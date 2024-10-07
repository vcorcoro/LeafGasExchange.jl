# CAM photosynthesis model derived from Photo3 model - Hartzell, Bartlett, & Porporato (2018)

@system CAMBase(CBase) begin
    # parameters from Photo3
    TH: high_temperature => 302.65 ~ preserve(parameter, u"K") 
    TW: low_temperature => 283.15 ~ preserve(parameter, u"K") 
    TOPT => 288.65 ~ preserve(u"K")
    Am_max: max_rate_of_malic_acid_storage_flux => 14 ~ preserve(parameter, u"μmol/m^2/s") # species dependent 
    TR: relaxation_time => 90 ~ preserve(u"minute") # Relaxation time for circadian oscillator (min)
    C0: parameter_for_decarboxylation_of_malic_acid => 3000 ~ preserve(u"μbar") # (umol/mol)

    Q: photosynthetic_active_radiation_conversion_factor => begin
        # 4.55 is a conversion factor from W to photons for solar radiation, Goudriaan and van Laar (1994)
        # some use 4.6 i.e., Amthor 1994, McCree 1981, Challa 1995.
        4.6
    end ~ preserve(u"μmol/J", parameter)

    solrad(PFD, Q): solar_radiation => begin
        PFD / Q
    end ~ track(u"W/m^2")

    # (Photo3)    
    Ci(Cs, fD, a1): intercellular_co2 => begin
        Cs * (1 - 1/(fD*a1))
    end ~ track(u"μbar")
    Cm(Ci): mesophyll_co2 ~ track(u"μbar")
    
    # malic acid decarboxylation increases co2 concentration in mesophyll (Photo3)
    Cc(Cm, f_c, C0): corrected_mesophyll_co2 => begin
        Cm + f_c*C0
    end ~ track(u"μbar")

end

# coppied from C3, parameters adjusted for CAM based on Photo3
# base temperature Tb = 20 instead of 25
@system CAMc(CAMBase) begin
    # 302 μmol/mol (at 293.2 K ~ 20.05 C) - Photo3
    Kc25: rubisco_constant_for_co2_at_25 => 302 ~ preserve(u"μbar", parameter)
    # 59430 J/mol (at 293.2 K ~ 20.05 C) - Photo3
    Eac: activation_energy_for_co2 => 59.43 ~ preserve(u"kJ/mol", parameter)
    Kc(Kc25, kT, Eac): rubisco_constant_for_co2 => begin
        Kc25 * kT(Eac)
    end ~ track(u"μbar")

    # 256. # Michaelis constant for 02 at TO (mmol/mol) - Photo3
    Ko25: rubisco_constant_for_o2_at_25 => 256 ~ preserve(u"mbar", parameter)
    # 36000. # Activation Energy for Ko (J/mol) - Photo3
    Eao: activation_energy_for_o2 => 36 ~ preserve(u"kJ/mol", parameter)
    Ko(Ko25, kT, Eao): rubisco_constant_for_o2 => begin
        Ko25 * kT(Eao)
    end ~ track(u"mbar")

    # mesophyll O2 partial pressure
    # OI = .209  # Oxygen Concentration (mol/mol) - Photo3
    Om: mesophyll_o2_partial_pressure => 209 ~ preserve(u"mbar", parameter)
    # effective M-M constant for Kc in the presence of O2
    Km(Kc, Om, Ko): rubisco_constant_for_co2_with_o2 => begin
        Kc * (1 + Om / Ko)
    end ~ track(u"μbar")

    # VCMAX0 = 18. # maximum carboxylation capacity - Photo3
    Vcm25: maximum_carboxylation_rate_at_25 => 18 ~ preserve(u"μmol/m^2/s" #= CO2 =#, parameter)
    # HKC =  59430. # Activation Energy for Kc (J/mol) - Photo3
    EaVc: activation_energy_for_carboxylation => 52.43 ~ preserve(u"kJ/mol", parameter)
    Vcmax(Vcm25, kT, EaVc, kN): maximum_carboxylation_rate => begin
        Vcm25 * kT(EaVc) * kN
    end ~ track(u"μmol/m^2/s" #= CO2 =#)
end

# coppied from C3, parameters adjusted for CAM based on Photo3
# base temperature Tb = 20 instead of 25
@system CAMj(CAMBase) begin
    # JMAX0 = 36. # maximum electron transport capacity - Photo3
    Jm25: maximum_electron_transport_rate_at_25 => 36 ~ preserve(u"μmol/m^2/s" #= Electron =#, parameter)
    # HAJ = 50000. # Activation Energy for Jmax (J/mol) - Photo3
    Eaj: activation_enthalpy_for_electron_transport => 50 ~ preserve(u"kJ/mol", parameter)
    # HDJ = 200000. # Deactivation Energy for Jmax (J/mol) - Photo3
    Hj: deactivation_enthalpy_for_electron_transport => 200 ~ preserve(u"kJ/mol", parameter)
    Sj: sensitivity_for_electron_transport => 616.4 ~ preserve(u"J/mol/K", parameter)
    # temperature dependent limitation
    Jmax(Jm25, kTpeak, Eaj, Hj, Sj, kN): maximum_electron_transport_rate => begin
        Jm25 * kTpeak(Eaj, Hj, Sj) * kN
    end ~ track(u"μmol/m^2/s" #= Electron =#)
    
    # θ: sharpness of transition from light limitation to light saturation
    # I2: effective irradiation absorbed by photosystem II, light dependent limitation
    # minh(I2, Jmax, θ)
    θ: light_transition_sharpness => 0.7 ~ preserve(parameter)
    J(I2, Jmax, θ): electron_transport_rate => begin
        a = θ
        b = -(I2+Jmax)
        c = I2*Jmax
        a*J^2 + b*J + c ⩵ 0
    end ~ solve(lower=0, upper=Jmax, u"μmol/m^2/s")
end

# coppied from C3, parameters adjusted for CAM based on Photo3
# base temperature Tb = 20 instead of 25
# addition of respiration fluxes - Photo3
@system CAMr(CAMBase) begin
    # 0.32 from Photo3 (at 293.2 K ~ 20.05 C)
    Rd25: dark_respiration_at_25 => 0.32 ~ preserve(u"μmol/m^2/s" #= CO2 =#, parameter)
    # 53000 J/mol from Photo3 (at 293.2 K ~ 20.05 C)
    Ear: activation_energy_for_respiration => 53 ~ preserve(u"kJ/mol", parameter)
    Rd(Rd25, kT, Ear): dark_respiration => begin
        Rd25 * kT(Ear)
    end ~ track(u"μmol/m^2/s")
    #Rm(Rd) => 0.5Rd ~ track(u"μmol/m^2/s")
    
    # equations from PHoto 3
    # Photo3 uses Φ which is incoming radiaion (W/m^2), same as solrad
    Rdv(Rd, solrad): dark_respiration_flux_to_calvin_cycle => begin
        Rd * exp(-Cropbox.deunitfy(solrad))
    end ~ track(u"μmol/m^2/s")
    Rdc(Rd, solrad): dark_respiration_flux_to_vacuole => begin
        Rd * (1 - exp(-Cropbox.deunitfy(solrad)))
    end ~ track(u"μmol/m^2/s")

    # CO2 compensation point in the absence of day respiration,
    # GAMMA_0 = 34.6 - Photo3
    Γ25: co2_compensation_point_at_25 => 34.6 ~ preserve(u"μbar", parameter)
    # different equation used, does not have this paramter - Photo3
    Eag: activation_energy_for_co2_compensation_point => 37.83 ~ preserve(u"kJ/mol", parameter)
    Γ(Γ25, kT, Eag): co2_compensation_point => begin
        Γ25 * kT(Eag)
    end ~ track(u"μbar")
end

###############
# From Photo3 #
###############
@system CircadianCycle(CAMBase) begin
   
    # Circadian oscillator constants
    MU => .5 ~ preserve
    BETA => 2.764 ~ preserve 
    CIRC_1 => .365 ~ preserve
    CIRC_2 => .55 ~ preserve
    CIRC_3 => 10 ~ preserve
    
    α_1 => 0.01 ~ preserve
    α_2 => 1/7 ~ preserve
    K => .003 ~ preserve
    
    m0: initial_malic_acid_concentration => 0 ~ preserve(u"mol/m^3")
    m_max: max_malic_acid_concentration => 230000000 ~ preserve(parameter, u"μmol/m^3") # species dependent
    vcm => 0.0027 ~ preserve(parameter, u"m") # Value controlling relative storage of malate (m)
    z0: initial_circadian_rhythm_order => 0.55 ~ preserve
    
    m_s(m_max,TH,Tk,TW,α_2): max_storage_concentration_of_malic_acid => begin
        m_max*((TH - Tk)/(TH - TW)*(1 - α_2) + α_2)
    end ~ track(u"mol/m^3")
    
    m_e(m_max,CIRC_1,TH,Tk,TW,BETA,z,MU,CIRC_2,f_o,m,α_1,solrad): equilibrium_concentration_of_malic_acid => begin
        if solrad > 0u"W/m^2"
            m_max*(CIRC_1*((TH - Tk)/(TH - TW) + 1)*(BETA*(z - MU))^3 - BETA*(TH - Tk)/(TH - TW)*(z - MU) + CIRC_2*(TH - Tk)/(TH - TW) -(1- f_o)*(1-m/(m+α_1*m_max)))
        else
            m_max*(CIRC_1*((TH - Tk)/(TH - TW) + 1)*(BETA*(z - MU))^3 - BETA*(TH - Tk)/(TH - TW)*(z - MU) + CIRC_2*(TH - Tk)/(TH - TW) + (1-f_o))
        end
    end ~ track(u"mol/m^3")
    
    f_o(z,MU,CIRC_3): circadian_order_function => begin
        exp(-(z/MU)^CIRC_3)
    end ~ track
    
    f_m(f_o,m_s,m,α_2): malic_acid_storage_function => begin
        f_o*(m_s - m)/(α_2*m_s + m_s - m)
    end ~ track
    
    f_c(f_o, m, α_1, m_max): carbon_circadian_control_function => begin
        (1 - f_o)*m/(α_1*m_max + m)
    end ~ track
    
    m(Asv, Avc, Rdv, vcm): malic_acid_concentration => begin
        (Asv - Avc + Rdv)/vcm
    end ~ accumulate(init=m0, min=0, u"μmol/m^3")

    z(m, m_e, m_max, TR): circadian_rhythm_order => begin
        (m - m_e)/(m_max*TR)
    end ~ accumulate(init=z0, min=0)
end

@system CAMRate(CAMc, CAMj, CAMr, CircadianCycle) begin
    
    # Same as C3, but not -Rd since respired carbon is recycled
    Ac_cc(Vcmax, Cc, Γ, Km, Rd): enzyme_limited_photosynthesis_rate_cc => begin
        Vcmax * (Cc - Γ) / (Cc + Km) # - Rd
    end ~ track(u"μmol/m^2/s" #= CO2 =#)
    
    Aj_cc(J, Cc, Γ, Rd): transport_limited_photosynthesis_rate_cc => begin
        J * (Cc - Γ) / 4(Cc + 2Γ) # - Rd
    end ~ track(u"μmol/m^2/s" #= CO2 =#)

    Ac_ci(Vcmax, Ci, Γ, Km, Rd): enzyme_limited_photosynthesis_rate_ci => begin
        Vcmax * (Ci - Γ) / (Ci + Km) # - Rd
    end ~ track(u"μmol/m^2/s" #= CO2 =#)
    
    Aj_ci(J, Ci, Γ, Rd): transport_limited_photosynthesis_rate_ci => begin
        J * (Ci - Γ) / 4(Ci + 2Γ) # - Rd
    end ~ track(u"μmol/m^2/s" #= CO2 =#)
    
    Ad_cc(Ac_cc, Aj_cc): co2_demand_function_cc => begin
        min(Ac_cc, Aj_cc)
    end ~ track(u"μmol/m^2/s" #= CO2 =#, min=0)

    Ad_ci(Ac_ci, Aj_ci): co2_demand_function_ci => begin
        min(Ac_ci, Aj_ci)
    end ~ track(u"μmol/m^2/s" #= CO2 =#, min=0)

    ###############
    # From Photo3 #
    ###############
    Asc(Ad_ci, Rdc, f_c, fΨv): co2_flux_stomata_to_calvin_cycle => begin
        (Ad_ci - Rdc)*(1-f_c) * fΨv
    end ~ track(min=0, u"μmol/m^2/s")
    
    Asv(m_s,m,K,Tk,T,TOPT,Am_max,Rdv,f_m,fΨv): co2_flux_stomata_to_vacuole => begin
        if m_s > m && 1 - K*(Cropbox.deunitfy(Tk - TOPT)^2) > 0
            (Am_max*(1 - K*(Cropbox.deunitfy(Tk - TOPT))^2) - Rdv)*f_m * fΨv
        else
            0
        end
    end ~ track(u"μmol/m^2/s")
    
    Avc(Ad_cc, Rdc, f_c): co2_flux_vacuole_to_calvin_cycle => begin
        (Ad_cc - Rdc)*f_c
    end ~ track(u"μmol/m^2/s")
    
    A_net(Asc, Asv): net_photosynthesis => begin
        Asc + Asv
    end ~ track(u"μmol/m^2/s" #= CO2 =#)
    
    # Respired carbon is recycled to calvin cycle or stored in vacuole, so A_net = A_gross
    A_gross(A_net): gross_photosynthesis ~ track(u"μmol/m^2/s" #= CO2 =#)
    # A_gross(A_net, Rd): gross_photosynthesis => begin
    #     A_net + Rd
    # end ~ track(u"μmol/m^2/s" #= CO2 =#)

end

@system CAM(CAMRate)