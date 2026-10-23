@system C4Base(CBase) begin
    Cm(Ci): mesophyll_co2 ~ track(u"μbar")

    gbs: bundle_sheath_conductance => 0.003 ~ preserve(u"mol/m^2/s/bar" #= CO2 =#, parameter) # bundle sheath conductance to CO2, mol m-2 s-1
    # gi => 1.0 ~ preserve(parameter) # conductance to CO2 from intercelluar to mesophyle, mol m-2 s-1, assumed
end

@system C4c(C4Base) begin
    # Kp25: Michaelis constant for PEP caboxylase for CO2
    Kp25: pep_carboxylase_constant_for_co2_at_25 => 80 ~ preserve(u"μbar", parameter)
    Kp(Kp25, kTQ10): pep_carboxylase_constant_for_co2 => begin
        Kp25 * kTQ10
    end ~ track(u"μbar")

    Vpm25: maximum_pep_carboxylation_rate_for_co2_at_25 => 70 ~ preserve(u"μmol/m^2/s" #= CO2 =#, parameter)
    EaVp: activation_energy_for_pep_carboxylation => 75.1 ~ preserve(u"kJ/mol", parameter)
    Vpmax(Vpm25, kT, EaVp, kN): maximum_pep_carboxylation_rate => begin
        Vpm25 * kT(EaVp) * kN
    end ~ track(u"μmol/m^2/s" #= CO2 =#)

    # PEP regeneration limited Vp, value adopted from vC book
    Vpr25: pep_regeneration_rate_for_co2_at_25 => 80 ~ preserve(u"μmol/m^2/s" #= CO2 =#, parameter)
    # Vpr(Vpr25, kTQ10): pep_regeneration_rate => begin
    #     Vpr25 * kTQ10
    # end ~ track(u"μmol/m^2/s" #= CO2 =#)
    # Vp_activation_type => 0 ~ preserve(parameter) # 0 for constant full activation, 1 for light dependent regeneration of PEP 
    Vp_activation_flag => false ~ flag(parameter) # true for dynamic light dependent regeneration of PEP - CC 10/8/26
    pa_init => 1 ~ preserve(parameter, min=p0, max=1)
    α_p: ATPase_activation_slope => 0.0017 ~ preserve(parameter, u"m^2*s/μmol") # Wang et al. (2021)
    p0: ATPase_activation_intercept => 0.16 ~ preserve(parameter) # Wang et al. (2021)
    τ_p: activation_time_constant_p => 30 ~ preserve(parameter, u"s") # Wang et al. (2021)
    pa_ss(α_p, p0, I2): ATPase_ss_activation_state => begin
        α_p * I2 + p0
    end ~ track(min=p0, max=1) # Wang et al. (2021)
    
    Spa(Spa, pa_ss, τ_p) => begin
        (pa_ss - Spa) / τ_p
    end ~ accumulate(init=pa_init, when=Vp_activation_flag)

    Vpr(Vpr25, kTQ10, Spa): pep_regeneration_rate => begin
        Vpr25 * kTQ10 * Spa
    end ~ track(u"μmol/m^2/s" #= CO2 =#)

    Vp(Vpmax, Vpr, Cm, Kp): pep_carboxylation_rate => begin
        # PEP carboxylation rate, that is the rate of C4 acid generation
        (Cm * Vpmax) / (Cm + Kp)
    end ~ track(u"μmol/m^2/s" #= CO2 =#, max=Vpr)

    Vcm25: maximum_carboxylation_rate_at_25 => 50 ~ preserve(u"μmol/m^2/s" #= CO2 =#, parameter)
    # EaVc: Sage (2002) JXB
    EaVc: activation_energy_for_carboxylation => 55.9 ~ preserve(u"kJ/mol", parameter)
    # Vcmax(Vcm25, kT, EaVc, kN): maximum_carboxylation_rate => begin
    #     Vcm25 * kT(EaVc) * kN
    # end ~ track(u"μmol/m^2/s" #= CO2 =#)

    ### Dyanmic light dependent activation of Rubisco -- CC 9/10/2026, updated 10/8/2026
    rubisco_activation_flag => false ~ flag(parameter) # true for dynamic rubisco activation
    # ra_init => 1 ~ preserve(parameter, min=r0, max=1)
    # α_r: initial_slope => 0.04 ~ preserve(parameter, u"m^2*s/μmol") # Gross et al. (1991) cited Seemann et al. (1988)
    # r0: minimum_state => 0.16 ~ preserve(parameter) # Gross et al. (1991) cited Seemann et al. (1988)
    # τ_r: activation_time_constant_r => 360 ~ preserve(parameter, u"s")
    # ra_ss(α_r, r0, I2): rubisco_ss_activation_state => begin
    #     α_r * I2 * (1 - r0) / (α_r * I2 + (1 - r0)) + r0      # Gross et al. (1991)
    # end ~ track(min=r0, max=1)

    # include light activation of Rubisco activase and subsequent activation or Rubisco (Wang et al., 2021) - CC 10/20/26
    α_R_init: initial_rubisco_activation_state => 1 ~ preserve(parameter, min=c0_r, max=1)
    α_RCA_init: initial_rubisco_activase_activation_state => 1 ~ preserve(parameter, min=c0_rca, max=1)
    k_rca: rubisco_activase_activation_slope => 0.0017 ~ preserve(parameter, u"m^2*s/μmol") # Wang et al. (2021)
    c0_rca: minimum_activation_state_rubisco_activase => 0.06 ~ preserve(parameter) # Wang et al. (2021)
    c0_r: minimum_activation_state_rubisco => 0 ~ preserve(parameter) # parameter value, source needed
    
    τ_r: activation_time_constant_rubisco => 360 ~ preserve(parameter, u"s")
    τ_RCA: activation_time_constant_rubisco_activase => 45.564 ~ preserve(parameter, u"s") # Wang et al. (2021)
    
    c_rac: rubisco_activase_concentration_constant => 216.9 ~ preserve(parameter, u"minute*mg/m^2") # Wang et al. (2021)
    K_activase: rubisco_activase_half_saturaton_constant => 12.3 ~ preserve(parameter, u"mg/m^2") # Wang et al. (2021)
 
    RCA(c_rac, τ_r): rubisco_activase_concentration => begin
        c_rac / τ_r
    end ~ track(u"mg/m^2")

    α_RCA_ss(k_rca, c0_rca, I2): rubisco_activase_ss_activation_state => begin
        k_rca * I2 + c0_rca      # Wang et al. (2021)
    end ~ track(min=c0_rca, max=1)

    α_RCA(α_RCA, α_RCA_ss, τ_RCA): rubisco_activase_activation_state => begin
        (α_RCA_ss - α_RCA) / τ_RCA
    end ~ accumulate(init=α_RCA_init, min=0, max=1)

    ActRCA(α_RCA, RCA): active_rubisco_activase_concentration => begin
        α_RCA * RCA
    end ~ track(u"mg/m^2")

    α_R_ss(ActRCA, K_activase): rubisco_activation_steady_state => begin
        ActRCA / (K_activase + ActRCA)
    end ~ track(min=0, max=1)

    α_R(α_R, α_R_ss, τ_r): rubisco_activation_state => begin
        (α_R_ss - α_R) / τ_r 
    end ~ accumulate(init=α_R_init, when=rubisco_activation_flag)

    Vcmax(Vcm25, kT, EaVc, kN, α_R): maximum_carboxylation_rate => begin
        Vcm25 * kT(EaVc) * kN * α_R
    end ~ track(u"μmol/m^2/s" #= CO2 =#)

end

@system C4j(C4Base) begin
    Jm25: maximum_electron_transport_rate_at_25 => 300 ~ preserve(u"μmol/m^2/s" #= Electron =#, parameter)
    Eaj: activation_enthalpy_for_electron_transport => 32.8 ~ preserve(u"kJ/mol", parameter)
    Hj: deactivation_enthalpy_for_electron_transport => 220 ~ preserve(u"kJ/mol", parameter)
    Sj: sensitivity_for_electron_transport => 702.6 ~ preserve(u"J/mol/K", parameter)
    Jmax(Jm25, kTpeak, Eaj, Hj, Sj, kN): maximum_electron_transport_rate => begin
        Jm25 * kTpeak(Eaj, Hj, Sj) * kN
    end ~ track(u"μmol/m^2/s" #= Electron =#)

    # θ: sharpness of transition from light limitation to light saturation
    θ: light_transition_sharpness => 0.5 ~ preserve(parameter)
    J(I2, Jmax, θ): electron_transport_rate => begin
        a = θ
        b = -(I2+Jmax)
        c = I2*Jmax
        a*J^2 + b*J + c ⩵ 0
    end ~ solve(lower=0, upper=Jmax, u"μmol/m^2/s")
end

@system C4r(C4Base) begin
    # Kc25: Michaelis constant of rubisco for CO2 of C4 plants (2.5 times that of tobacco), ubar, Von Caemmerer 2000
    Kc25: rubisco_constant_for_co2_at_25 => 650 ~ preserve(u"μbar", parameter)
    Eac: activation_energy_for_co2 => 59.4 ~ preserve(u"kJ/mol", parameter)
    Kc(kT, Kc25, Eac): rubisco_constant_for_co2 => begin
        Kc25 * kT(Eac)
    end ~ track(u"μbar")

    # Ko25: Michaelis constant of rubisco for O2 (2.5 times C3), mbar
    Ko25: rubisco_constant_for_o2_at_25 => 450 ~ preserve(u"mbar", parameter)
    # Activation energy for Ko, Bernacchi (2001)
    Eao: activation_energy_for_o2 => 36 ~ preserve(u"kJ/mol", parameter)
    Ko(Ko25, kT, Eao): rubisco_constant_for_o2 => begin
        Ko25 * kT(Eao)
    end ~ track(u"mbar")

    # mesophyll O2 partial pressure
    Om: mesophyll_o2_partial_pressure => 210 ~ preserve(u"mbar", parameter)
    Km(Kc, Om, Ko): rubisco_constant_for_co2_with_o2 => begin
        # effective M-M constant for Kc in the presence of O2
        Kc * (1 + Om / Ko)
    end ~ track(u"μbar")

    # Rd25: Values in Kim (2006) are for 31C, and the values here are normalized for 25C. SK
    Rd25: dark_respiration_at_25 => 2 ~ preserve(u"μmol/m^2/s" #= CO2 =#, parameter)
    Ear: activation_energy_for_respiration => 39.8 ~ preserve(u"kJ/mol", parameter)
    Rd(Rd25, kT, Ear): dark_respiration => begin
        Rd25 * kT(Ear)
    end ~ track(u"μmol/m^2/s")
    Rm(Rd) => 0.5Rd ~ track(u"μmol/m^2/s")
end

@system C4Rate(C4c, C4j, C4r) begin
    # Enzyme limited A (Rubisco or PEP carboxylation)
    Ac1(Vp, gbs, Cm, Rm) => (Vp + gbs*Cm - Rm) ~ track(u"μmol/m^2/s" #= CO2 =#)
    Ac2(Vcmax, Rd) => (Vcmax - Rd) ~ track(u"μmol/m^2/s" #= CO2 =#)
    Ac(Ac1, Ac2): enzyme_limited_photosynthesis_rate => begin
        #Ac1 = max(0, Ac1) # prevent Ac1 from being negative Yang 9/26/06
        min(Ac1, Ac2)
    end ~ track(u"μmol/m^2/s" #= CO2 =#)

    # x: Partitioning factor of J, yield maximal J at this value
    x: electron_transport_partitioning_factor => 0.4 ~ preserve(parameter)
    # Light and electron transport limited A mediated by J
    Aj1(x, J, Rm, gbs, Cm) => (x * J/2 - Rm + gbs*Cm) ~ track(u"μmol/m^2/s" #= CO2 =#)
    Aj2(x, J, Rd) => (1-x) * J/3 - Rd ~ track(u"μmol/m^2/s" #= CO2 =#)
    Aj(Aj1, Aj2): transport_limited_photosynthesis_rate => begin
        min(Aj1, Aj2)
    end ~ track(u"μmol/m^2/s" #= CO2 =#)

    # smoothing the transition between Ac and Aj
    β: photosynthesis_transition_factor => 0.99 ~ preserve(parameter)
    A_net(Ac, Aj, β): net_photosynthesis => begin
        x = A_net
        a = β
        b = -(Ac+Aj)
        c = Ac*Aj
        a*x^2 + b*x + c ⩵ 0
    end ~ solve(pick=:minimum, u"μmol/m^2/s")

    A_gross(A_net, Rd): gross_photosynthesis => begin
        A_gross = A_net + Rd
        # gets negative when PFD = 0, Rd needs to be examined, 10/25/04, SK
        #max(A_gross, zero(A_gross))
    end ~ track(u"μmol/m^2/s" #= CO2 =#)
end

@system C4(C4Rate) begin
    #FIXME: currently not used variables

    # FIXME are they even used?
    # beta_ABA => 1.48e2 # Tardieu-Davies beta, Dewar (2002) Need the references !?
    # delta => -1.0
    # alpha_ABA => 1.0e-4
    # lambda_r => 4.0e-12 # Dewar's email
    # lambda_l => 1.0e-12
    # K_max => 6.67e-3 # max. xylem conductance (mol m-2 s-1 MPa-1) from root to leaf, Dewar (2002)

    # alpha: fraction of PSII activity in the bundle sheath cell, very low for NADP-ME types
    α: bundle_sheath_PSII_activity_fraction => 0.0001 ~ preserve(parameter)
    a_o: ratio_of_solubility_and_diffusivity_O2_to_CO2 => 0.047 ~ preserve(parameter)
    # Bundle sheath O2 partial pressure, mbar
    Os(A_net, gbs, Om, α, a_o): bundle_sheath_o2 => begin
        α * A_net / (a_o * gbs) + Om
    end ~ track(u"mbar")

    Cbs(A_net, Vp, Cm, Rm, gbs): bundle_sheath_co2 => begin
        Cm + (Vp - A_net - Rm) / gbs # Bundle sheath CO2 partial pressure, ubar
    end ~ track(u"μbar")

    # half the reciprocal of rubisco specificity, to account for O2 dependence of CO2 comp point,
    # note that this become the same as that in C3 model when multiplied by [O2]
    Γ1 => 0.193 ~ preserve(parameter)
    Γ★(Γ1, Os) => Γ1 * Os ~ track(u"μbar")
    Γ(Rd, Km, Vcmax, Γ★): co2_compensation_point => begin
        (Rd*Km + Vcmax*Γ★) / (Vcmax - Rd)
    end ~ track(u"μbar")
end

@system C4Dyn(C4) begin
    A_net(Ac, Aj): net_photosynthesis => begin
        min(Ac, Aj)
    end ~ track(u"μmol/m^2/s")
end

@system C4Exact(C4) begin

    a_c(α, a_o, Kc, Ko) => begin
        1 - (α / a_o * Kc / Ko)
    end ~ track(u"μbar/μbar")

    b_c(Vp, Vcmax, Rm, Rd, gbs, Cm, Om, Kc, Ko, Γ1, α, a_o) => begin
        -((Vp - Rm + gbs*Cm) + (Vcmax - Rd) + gbs*(Kc*(1 + Om/Ko)) + α/a_o*(Γ1*Vcmax + Rd*Kc/Ko))
    end ~ track(u"μmol/m^2/s")

    c_c(Vp, Vcmax, Rm, Rd, gbs, Cm, Om, Kc, Ko, Γ1) => begin
        (Vcmax - Rd) * (Vp - Rm + gbs*Cm) - (Vcmax*gbs*Γ1*Om + Rd*gbs*Kc*(1+Om/Ko))
    end ~ track(u"μmol^2/m^4/s^2")
    
    Ac_min(Rd) => -Rd ~ track(u"μmol/m^2/s")
    Ac_max(Vp, gbs, Cm, Rm, Vcmax, Rd) => begin
        Ac1 = (Vp + gbs*Cm - Rm)
        Ac2 = (Vcmax - Rd)
        min(Ac1, Ac2)
    end ~ track(u"μmol/m^2/s" #= CO2 =#)

    Ac(a = a_c, b = b_c, c = c_c): enzyme_limited_photosynthesis_rate => begin
        a*Ac^2 + b*Ac + c
    end ~ solve(u"μmol/m^2/s", lower=Ac_min, upper=Ac_max)
    # other option: input quadratic equation directly taking the (-b-√(b^2 - 4ac))/(2a) root (von Caemmerer, 2000)
    
    # a_j(Γ1, α, a_o) => begin
    #     1 - (7*Γ1*α) / (3 * a_o)
    # end ~ track

    # b_j(z, x, J, Rm, Rd, Γ1, Cm, Om, α, a_o) => begin
    #     -((z/2*x*J - Rm + gbs*Cm) + (z/3*(1-x)*J - Rd) + gbs*(7*Γ1*Om/3) + α*Γ1/a_o*(z/3*(1-x)*J + 7/3*Rd))
    # end ~ track(u"μmol/m^2/s") 

    # c_j(z, x, J, Rm, Rd, Γ1, Cm, Om) => begin
    #     (z/2*x*J - Rm + gbs*Cm) * (z/3*(1-x)*J - Rd) - gbs*Γ1*Om*(z/3*(1-x)*J + 7/3*Rd)
    # end ~ track(u"μmol^2/m^4/s^2")

    # Aj(a = a_j, b = b_j, c = c_j): transport_limited_photosynthesis_rate => begin
    #     (-b + (b^2 - 4*a*c)^(-1/2))  / (2*a)
    # end ~ track(u"μmol/m^2/s")

    A_net(Ac, Aj): net_photosynthesis => begin
        min(Ac, Aj)
    end ~ track(u"μmol/m^2/s")
end