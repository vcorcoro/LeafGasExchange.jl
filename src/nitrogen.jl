@system Nitrogen begin
    SPAD: SPAD_greenness => 60 ~ preserve(parameter)
    SNa: SPAD_N_coeff_a => 0.0004 ~ preserve(u"g/m^2", parameter)
    SNb: SPAD_N_coeff_b => 0.012 ~ preserve(u"g/m^2", parameter)
    SNc: SPAD_N_coeff_c => 0 ~ preserve(u"g/m^2", parameter)
    N(SPAD, a=SNa, b=SNb, c=SNc): leaf_nitrogen_content => begin
        a*SPAD^2 + b*SPAD + c
    end ~ preserve(u"g/m^2", parameter)

    Np(N, SLA) => N * SLA ~ track(u"percent")
    SLA: specific_leaf_area => 200 ~ preserve(u"cm^2/g")
end
