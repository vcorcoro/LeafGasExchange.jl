@system StomataBase(Weather, Diffusion) begin
    gs: stomatal_conductance ~ hold
    gb: boundary_layer_conductance ~ hold
    A_net: net_photosynthesis ~ hold
    T: leaf_temperature ~ hold

    drb(Dw, Dc): diffusivity_ratio_boundary_layer => (Dw / Dc)^(2/3) ~ preserve(#= u"H2O/CO2", =# parameter)
    dra(Dw, Dc): diffusivity_ratio_air => (Dw / Dc) ~ preserve(#= u"H2O/CO2", =# parameter)

    Ca(CO2, P_air): co2_air => (CO2 * P_air) ~ track(u"μbar")
    Cs(Ca, A_net, gbc): co2_at_leaf_surface => begin
        Ca - A_net / gbc
    end ~ track(u"μbar")

    gv(gs, gb): total_conductance_h2o => (gs * gb / (gs + gb)) ~ track(u"mol/m^2/s/bar" #= H2O =#)

    rbc(gb, drb): boundary_layer_resistance_co2 => (drb / gb) ~ track(u"m^2*s/mol*bar")
    rsc(gs, dra): stomatal_resistance_co2 => (dra / gs) ~ track(u"m^2*s/mol*bar")
    rvc(rbc, rsc): total_resistance_co2 => (rbc + rsc) ~ track(u"m^2*s/mol*bar")

    gbc(rbc): boundary_layer_conductance_co2 => (1 / rbc) ~ track(u"mol/m^2/s/bar")
    gsc(rsc): stomatal_conductance_co2 => (1 / rsc) ~ track(u"mol/m^2/s/bar")
    gvc(rvc): total_conductance_co2 => (1 / rvc) ~ track(u"mol/m^2/s/bar")
end

@system StomataTuzet begin
    WP_leaf: leaf_water_potential => 0 ~ preserve(u"MPa", parameter)
    Ψv(WP_leaf): bulk_leaf_water_potential ~ track(u"MPa")
    Ψf: reference_leaf_water_potential => -2.0 ~ preserve(u"MPa", parameter)
    sf: stomata_sensitivity_param => 2.3 ~ preserve(u"MPa^-1", parameter)
    fΨv(Ψv, Ψf, sf): stomata_sensitivty => begin
        (1 + exp(sf*Ψf)) / (1 + exp(sf*(Ψf-Ψv)))
    end ~ track
end

@system StomataBallBerry(StomataBase, StomataTuzet) begin
    # Set default Ball-Berry model parameter values for C3 plants. Assume g0 is not different from 0. 
    # See Franks et al (2017) Plant Physiology and Miner et al (2017) Plant Cell Environ
    # For temperate C4 species, use g1 = 5.2.  
    g0 => 0.0 ~ preserve(u"mol/m^2/s/bar" #= H2O =#, parameter)
    g1 => 13.1 ~ preserve(parameter)

    # include g0 as minimum gs in hs equation, added lower=RH to preserve concentration gradient assumption
    hs1(g0, g1, gb, A_net, Cs, RH): relative_humidity_at_leaf_surface1 => begin
        gs = g0 + g1*(A_net*hs1/Cs) 
        (hs1 - RH)*gb ⩵ (1 - hs1)*gs
    end ~ solve(lower=RH, upper=1)
    # solution when gs = g0 (i.e. A_net ≈ 0)
    hs2(g0, gb, RH): relative_humidity_at_leaf_surface2 => begin
        # (hs2 - RH)*gb ⩵ (1 - hs2)*g0
        (g0 + RH*gb) / (gb + g0)
    end ~ track(min=RH, max=1)
    hs(hs1, hs2, nounit(A_net)): relative_humidity_at_leaf_surface => begin
        if (A_net + 1) ≈ 1    # check A_net ≈ 0
            hs2
        else 
            hs1
        end
    end ~ track

    #HACK: avoid scaling issue with dimensionless unit
    # hs(g0, g1, gb, A_net, Cs, fΨv, RH): relative_humidity_at_leaf_surface => begin
    #     gs = g0 + g1*(A_net*hs/Cs) * fΨv 
    #         (hs - RH)*gb ⩵ (1 - hs)*gs
    # end ~ solve(lower=0, upper=1) #, u"percent")
    Ds(D=vp.D, T, hs): vapor_pressure_deficit_at_leaf_surface => begin
        D(T, hs)
    end ~ track(u"kPa")

    gs(g0, g1, A_net, hs, Cs, fΨv): stomatal_conductance => begin
        g0 + g1*(A_net*hs/Cs) * fΨv
    end ~ track(u"mol/m^2/s/bar" #= H2O =#, min=g0)
end

@system StomataMedlyn(StomataBase, StomataTuzet) begin
    # Set default Medlyn model parameters for C3 plants. 
    # See Franks et al (2017) Plant Physiology (http://www.plantphysiol.org/cgi/doi/10.1104/pp.17.00287)
    # See also Lin et al. (2015) Nature Climate Change 
    g0 => 0.0 ~ preserve(u"mol/m^2/s/bar" #= H2O =#, parameter)
    g1 => 4.45 ~ preserve(u"√kPa", parameter)

    wa(ea=vp.ea, T_air, RH): vapor_pressure_at_air => ea(T_air, RH) ~ track(u"kPa")
    wi(es=vp.es, T): vapor_pressure_at_intercellular_space => es(T) ~ track(u"kPa")
    ws(Ds, wi): vapor_pressure_at_leaf_surface => (wi - Ds) ~ track(u"kPa")
    Ds¹ᐟ²(g0, g1, gb, A_net, Cs, fΨv, wi, wa) => begin
        #HACK: SymPy couldn't extract polynomial coeffs for ps inside √
        gs = g0 + (1 + g1 / Ds¹ᐟ²) * (A_net / Cs) * fΨv
        ws = wi - Ds¹ᐟ²^2
        (ws - wa)*gb ⩵ (wi - ws)*gs
    end ~ solve(lower=0, upper=√wi', u"√kPa")
    Ds(Ds¹ᐟ²): vapor_pressure_deficit_at_leaf_surface => Ds¹ᐟ²^2 ~ track(u"kPa", min=1u"Pa")
    hs(RH=vp.RH, T, Ds): relative_humidity_at_leaf_surface => RH(T, Ds) ~ track

    gs(g0, g1, A_net, Ds, Cs, fΨv): stomatal_conductance => begin
        g0 + (1 + g1/√Ds)*(A_net/Cs) * fΨv
    end ~ track(u"mol/m^2/s/bar" #= H2O =#, min=g0)
end

# three step mechanism for dynamic stomata response (Kirschbaum et al., 1988)
@system StomataKirschbaum(StomataBase, StomataTuzet) begin
    h(context.clock.step) ~ track(u"hr")
    gs_high => 0.03 ~ preserve(parameter, u"mol/m^2/s/bar")

    S_eq ~ hold # depends on photosynthesis type
    S_eq0: initial_stomata_equilibrium => 1 ~ preserve(parameter)

    τ_i => 0.37 ~ preserve(parameter, u"minute")
    τ_d => 7 ~ preserve(parameter, u"minute")
    τ_π => 13 ~ preserve(parameter, u"minute")
    τ_w => 15 ~ preserve(parameter, u"minute")
    
    dS(S_eq, S, τ_i, τ_d, h) => begin
        function f(S)
            if S_eq > S
                (S_eq - S) / τ_i
            else
                (S_eq - S) / τ_d
            end
        end
        k1 = f(S)
        k2 = f(S + k1*h/2)
        k3 = f(S + k2*h/2)
        k4 = f(S + h*k3)
        dS = (k1 + 2k2 + 2k3 + k4) / 6
    end ~ track(u"hr^-1")

    S(dS) ~ accumulate(init=S_eq0, max=1)

    dπ(π, S, τ_π, h) => begin
        function f(π)
            (S - π) / τ_π
        end
        k1 = f(π)
        k2 = f(π + k1*h/2)
        k3 = f(π + k2*h/2)
        k4 = f(π + h*k3)
        dπ = (k1 + 2k2 + 2k3 + k4) / 6
    end ~ track(u"hr^-1")

    π(dπ) ~ accumulate(init=S_eq0, max=1)

    dwc(wc, π, τ_w, h) => begin
        function f(wc)
            (π - wc) / τ_w
        end
        k1 = f(wc)
        k2 = f(wc + k1*h/2)
        k3 = f(wc + k2*h/2)
        k4 = f(wc + h*k3)
        dw = (k1 + 2k2 + 2k3 + k4) / 6
    end ~ track(u"hr^-1")

    wc(dwc) ~ accumulate(init=S_eq0, max=1)

    gs(wc, gs_high) => wc * gs_high ~ track(u"mol/m^2/s/bar")

    # gvc value for use in IntercellularSpaceDyn
    g(gvc) ~ track(u"mol/m^2/s/bar")
end

# CAM biochemical signal S as a function of circadian oscillator and leaf water potential
@system StomataKirschbaumCAM(StomataKirschbaum) begin
    z ~ hold
    S_eq(z, fΨv) => (1 - z) * fΨv ~ track # TODO: need a version of g0/S_min?
end

# C3 biochemical signal S as a function of light and leaf water potential
@system StomataKirschbaumC3(StomataKirschbaum) begin
    I2 ~ hold
    θs => 0.866 ~ preserve(parameter)
    α => 0.0088 ~ preserve(parameter, u"m^2*s/μmol")
    SI_min => 0 ~ preserve(parameter)

    SI_var(I2,θs,α,SI_min) => begin
        a = θs
        b = -(1 + SI_min + α*I2)
        c = (SI_min + α*I2)
        x = SI_var + SI_min
        a*x^2 + b*x + c
    end ~ solve(pick=:minimum)

    SI_eq(SI_var, SI_min) => SI_var + SI_min ~ track  

    SΨ_eq(fΨv) ~ track

    S_eq(SI_eq, SΨ_eq) => SI_eq * SΨ_eq ~ track
end

# dynamic stomata model - prognostic updates with time constant (Liu et al., 2024)
@system StomataDyn(StomataBase, StomataTuzet, StomataBallBerry) begin
    # initialize stomatal conductance
    gs_init => 0.001 ~ preserve(parameter, u"mol/m^2/s/bar") 
    gvc_init(gs_init, gb, drb, dra) => begin
        1 / (drb/gb + dra/gs_init)
    end ~ preserve(u"mol/m^2/s/bar")

    # capture stomatal conductance at last time step
    gvc_capture(gvc) ~ capture(u"mol/m^2/s/bar*hr", time=ztime)
    gvc_prev(gvc_init, gvc_capture, h=context.clock.step, i=context.clock.tick) => begin 
        i == 0 ? gvc_init : gvc_capture / h     
    end ~ track(u"mol/m^2/s/bar")
    # gvc value for use in IntercellularSpaceDyn
    g(gvc_prev) ~ track(u"mol/m^2/s/bar", min=g0)

    gs_capture(gs) ~ capture(u"mol/m^2/s/bar*hr", time=ztime)
    gs_prev(gs_init, gs_capture, h=context.clock.step, i=context.clock.tick) => begin
        i == 0 ? gs_init : gs_capture / h 
    end ~ track(u"mol/m^2/s/bar", min=g0)
    
    # calculate target steady state conductance from Ball-Berry model
    gs_target(g0, g1, A_net, hs, Cs, fΨv): stomatal_conductance_steady_state_target => begin
        if any(i -> isnan(i), [A_net, hs, Cs])
            @warn "NaN in equation for gs_target: A_net => $A_net, hs => $hs, Cs => $Cs"
        end
        g0 + g1*(A_net*hs/Cs) * fΨv
    end ~ track(u"mol/m^2/s/bar" #= H2O =#, min=g0)
    
    # calculate actual change in stomatal conductance 
    τ: stomata_time_constant => 900 ~ preserve(parameter, u"s")
    Δgs(gs_prev, gs_target, τ, h=context.clock.step) => begin
        (gs_target - gs_prev) / τ * h
    end ~ track(u"mol/m^2/s/bar")

    # update gs
    gs(gs_prev, Δgs): stomatal_conductance => begin
        gs_prev + Δgs
    end ~ track(u"mol/m^2/s/bar", min=g0)
end

###############
# From Photo3 #
###############
# @system StomataCAM(StomataBase, StomataTuzet, StomataMedlyn) begin
#     a1: stomatal_conductance_parameter => 0.8*15 ~ preserve(parameter)

#     fD(VPD): stomata_response_to_vpd => begin
#         if VPD < 0.01u"Pa"
#             1
#         else 
#             3/13/sqrt(Cropbox.deunitfy(VPD)) # no /1000 since in kPa instead of Pa
#         end
#     end ~ track

#     # in Photo3 Cs = Ca
#     Cs(Ca): co2_at_leaf_surface ~ track(u"μbar")

#     gsc(a1, A_net, Cs, fD): stomatal_conductance_co2 => begin
#         a1*A_net/Cs*fD
#     end ~ track(u"mol/m^2/s/bar" #= C02 =#)

#     gs(gsc): stomatal_conductance_to_water => begin
#         gsc*1.6
#     end ~ track(u"mol/m^2/s/bar" #= H2O =#)
# end
