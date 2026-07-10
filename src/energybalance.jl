@system EnergyBalance(Weather) begin
    gv ~ hold
    gh ~ hold
    PPFD ~ hold

    ϵ: leaf_thermal_emissivity => 0.97 ~ preserve(parameter)
    σ: stefan_boltzmann_constant => u"σ" ~ preserve(u"W/m^2/K^4")
    λ: latent_heat_of_vaporization_at_25 => 44 ~ preserve(u"kJ/mol", parameter)
    Cp: specific_heat_of_air => 29.3 ~ preserve(u"J/mol/K", parameter)

    k: radiation_conversion_factor => (1 / 4.55) ~ preserve(u"J/μmol")
    PAR(PPFD, k): photosynthetically_active_radiation => (PPFD * k) ~ track(u"W/m^2")

    # NIR(PAR): near_infrared_radiation => begin
    #     #FIXME: maybe δ or similar ratio supposed to be applied here?
    #     # If total solar radiation unavailable, assume NIR the same energy as PAR waveband
    #     PAR
    # end ~ track(u"W/m^2")

    # solar radiation absorptivity of leaves: =~ 0.5
    #FIXME: is α different from (1 - δ) in Irradiance?
    α_s: absorption_coefficient => 0.5 ~ preserve(parameter)

    #R_sw(PAR, NIR, α_s, δ): shortwave_radiation_absorbed => begin
    R_sw(PAR, α_s): shortwave_radiation_absorbed => begin
        #FIXME: why δ needed here? α should already take care of scattering
        # shortwave radiation (PAR (=0.85) + NIR (=0.15))
        #α_s*((1-δ)*PAR + δ*NIR)
        α_s*PAR
    end ~ track(u"W/m^2")

    R_wall(ϵ, σ, Tk_air): thermal_radiation_absorbed_from_wall => 2ϵ*σ*Tk_air^4 ~ track(u"W/m^2")
    R_leaf(ϵ, σ, Tk): thermal_radiation_emitted_by_leaf => 2ϵ*σ*Tk^4 ~ track(u"W/m^2")
    R_thermal(R_wall, R_leaf): thermal_radiation_absorbed => R_wall - R_leaf ~ track(u"W/m^2")
    R_net(R_sw, R_thermal): net_radiation_absorbed => R_sw + R_thermal ~ track(u"W/m^2")

    Δw(T, T_air, RH, #= P_air, =# ea=vp.ambient, es=vp.saturation): leaf_vapor_pressure_gradient => begin
        Es = es(T)
        Ea = ea(T_air, RH)
        Es - Ea # MAIZSIM: / (1 - (Es + Ea) / P_air)
    end ~ track(u"kPa")
    E(gv, Δw): transpiration => gv*Δw ~ track(u"mmol/m^2/s" #= H2O =#)

    H(Cp, gh, ΔT): sensible_heat_flux => Cp*gh*ΔT ~ track(u"W/m^2")
    λE(λ, E): latent_heat_flux => λ*E ~ track(u"W/m^2")

    ΔT(R_net, H, λE): temperature_adjustment => begin
        R_net ⩵ H + λE
    end ~ bisect(lower=-5, upper=5, u"K", evalunit=u"W/m^2")

    T(T_air, ΔT): leaf_temperature => (T_air + ΔT) ~ track(u"°C")
    Tk(T): absolute_leaf_temperature ~ track(u"K")

    # # Test taylor expansion
    # a_ξ(tc=nounit(T_air), tk=nounit(Tk_air), T_air, Tk_air, ϵ, σ, λ, gv, b=vp.b, c=vp.c, es=vp.es) => begin
    #     Esa = es(T_air)
    #     λ*gv*Esa * (b*c*tk * (b*c*tk - 2*tk*(c + tc)) / (c + tc)^4) + 24ϵ*σ*Tk_air^4
    # end ~ track(u"W/m^2")

    # b_ξ(tc=nounit(T_air), tk=nounit(Tk_air), T_air, Tk_air, ϵ, σ, λ, gv, Cp, gh, b=vp.b, c=vp.c, es=vp.es) => begin
    #     Esa = es(T_air)
    #     b = Cp*gh*Tk_air + λ*gv*Esa*(b*c*tk) / (c + tc)^2 + 8ϵ*σ*Tk_air^4
    # end ~ track(u"W/m^2")
    
    # c_ξ(T_air, λ, gv, RH, R_sw, es=vp.es) => begin
    #     Esa = es(T_air)
    #     c = λ*gv*Esa*(1 - RH) - R_sw
    # end ~ track(u"W/m^2")

    # ξ(a_ξ, b_ξ, c_ξ): temperature_adjustment_TE => begin
    #     0.5*a_ξ*ξ^2 + b_ξ*ξ + c_ξ
    # end ~ solve

    # # Tk(Tk_air, ξ) => Tk_air * (1 + ξ) ~ track(u"K")
    # # T(Tk) ~ track(u"°C")

    # Tk_test(Tk_air, ξ, Tk, i=context.clock.tick) => begin
    #     Tk_test = Tk_air * (1 + ξ)
    #     @info "[$i] updating temperature: bisect - $Tk, solve - $Tk_test"
    #     Tk_test
    # end ~ track(u"K")
end

@system EnergyBalanceDyn(Weather) begin
    gv ~ hold
    gh ~ hold
    PPFD ~ hold

    ϵ: leaf_thermal_emissivity => 0.97 ~ preserve(parameter)
    σ: stefan_boltzmann_constant => u"σ" ~ preserve(u"W/m^2/K^4")
    λ: latent_heat_of_vaporization_at_25 => 44 ~ preserve(u"kJ/mol", parameter)
    Cp: specific_heat_of_air => 29.3 ~ preserve(u"J/mol/K", parameter)

    k: radiation_conversion_factor => (1 / 4.55) ~ preserve(u"J/μmol")
    PAR(PPFD, k): photosynthetically_active_radiation => (PPFD * k) ~ track(u"W/m^2")

    # solar radiation absorptivity of leaves: =~ 0.5
    #FIXME: is α different from (1 - δ) in Irradiance?
    α_s: absorption_coefficient => 0.5 ~ preserve(parameter)

    #R_sw(PAR, NIR, α_s, δ): shortwave_radiation_absorbed => begin
    R_sw(PAR, α_s): shortwave_radiation_absorbed => begin
        #FIXME: why δ needed here? α should already take care of scattering
        # shortwave radiation (PAR (=0.85) + NIR (=0.15))
        #α_s*((1-δ)*PAR + δ*NIR)
        α_s*PAR
    end ~ track(u"W/m^2")

    # R_wall(ϵ, σ, Tk_air): thermal_radiation_absorbed_from_wall => 2ϵ*σ*Tk_air^4 ~ track(u"W/m^2")
    # R_leaf(ϵ, σ, Tk): thermal_radiation_emitted_by_leaf => 2ϵ*σ*Tk^4 ~ track(u"W/m^2")
    # R_thermal(R_wall, R_leaf): thermal_radiation_absorbed => R_wall - R_leaf ~ track(u"W/m^2")
    # R_net(R_sw, R_thermal): net_radiation_absorbed => R_sw + R_thermal ~ track(u"W/m^2")

    # Taylor expansion with Tk = Tk_air(1 + ξ) substituted into R_net ⩵ H + λE
    # f(ξ) = H + λE - R_net

    # f''(0)
    a_ξ(tc=nounit(T_air), tk=nounit(Tk_air), T_air, Tk_air, ϵ, σ, λ, gv, b=vp.b, c=vp.c, es=vp.es) => begin
        Esa = es(T_air)
        λ*gv*Esa * (b*c*tk * (b*c*tk - 2*tk*(c + tc)) / (c + tc)^4) + 24ϵ*σ*Tk_air^4
    end ~ track(u"W/m^2")

    # f'(0)
    b_ξ(tc=nounit(T_air), tk=nounit(Tk_air), T_air, Tk_air, ϵ, σ, λ, gv, Cp, gh, b=vp.b, c=vp.c, es=vp.es) => begin
        Esa = es(T_air)
        b = Cp*gh*Tk_air + λ*gv*Esa*(b*c*tk) / (c + tc)^2 + 8ϵ*σ*Tk_air^4
    end ~ track(u"W/m^2")
    
    # f(0)
    c_ξ(T_air, λ, gv, RH, R_sw, es=vp.es) => begin
        Esa = es(T_air)
        c = λ*gv*Esa*(1 - RH) - R_sw
    end ~ track(u"W/m^2")

    # second order Taylor expansion around ξ = 0, f(ξ) ≈ f(0) + f'(0)*ξ + 0.5*f''(0)*ξ^2
    ξ(a_ξ, b_ξ, c_ξ): temperature_adjustment_TE => begin
        0.5*a_ξ*ξ^2 + b_ξ*ξ + c_ξ
    end ~ solve
    
    Tk_target(Tk_air, ξ) => Tk_air * (1 + ξ) ~ track(u"K")

    τ_T: temperature_time_constant => 900 ~ preserve(u"s", parameter)

    # calculate actual change in T
    ΔTk(Tk_target, Tk, τ_T, h=context.clock.step) => begin
        (Tk_target - Tk) / τ_T * h
    end ~ track(u"K")

    # update T
    Tk_next(Tk, ΔTk, i=context.clock.tick): absolute_leaf_temperature_next => begin
        Tk + ΔTk
    end ~ track(u"K")

    # default initialization Tk = Tk_air, otherwise initialize in config
    Tk_init(Tk_air) ~ preserve(parameter, u"K")

    Tk_capture(Tk_next) ~ capture(u"K*hr", time=ztime)
    Tk(Tk_init, Tk_capture, h=context.clock.step, i=context.clock.tick): absolute_leaf_temperature => begin
        i == 0 ? Tk_init : Tk_capture / h 
    end ~ track(u"K")
    T(Tk): leaf_temperature ~ track(u"°C")
end
