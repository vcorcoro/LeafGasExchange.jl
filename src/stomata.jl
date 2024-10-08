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

    #HACK: avoid scaling issue with dimensionless unit
    hs(g0, g1, gb, A_net, Cs, fΨv, RH): relative_humidity_at_leaf_surface => begin
        gs = g0 + g1*(A_net*hs/Cs) * fΨv
        (hs - RH)*gb ⩵ (1 - hs)*gs
    end ~ solve(lower=0, upper=1) #, u"percent")
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

###############
# From Photo3 #
###############
@system StomataCAM(StomataBase, StomataTuzet) begin
    a1: stomatal_conductance_parameter => 0.8*15 ~ preserve(parameter)

    fD(VPD): stomata_response_to_vpd => begin
        if VPD < 0.01u"Pa"
            1
        else 
            3/13/sqrt(Cropbox.deunitfy(VPD)) # no /1000 since in kPa instead of Pa
        end
    end ~ track

    # in Photo3 Cs = Ca
    Cs(Ca): co2_at_leaf_surface ~ track(u"μbar")

    Cc ~ hold

    gsc(a1, A_net, Cc, fD): stomatal_conductance_co2 => begin
        a1*A_net/Cc*fD
    end ~ track(u"mol/m^2/s/bar" #= C02 =#)

    gs(gsc): stomatal_conductance_to_water => begin
        gsc*1.6
    end ~ track(u"mol/m^2/s/bar" #= H2O =#)
end
