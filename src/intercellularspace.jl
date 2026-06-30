@system IntercellularSpace(Weather) begin
    A_net ~ hold
    #TODO: interface between boundary/stomata/intercellular space (i.e. soil layers?)
    gvc ~ hold

    #FIXME: duplicate in Stomata
    Ca(CO2, P_air): co2_air => (CO2 * P_air) ~ track(u"μbar")

    #HACK: high temperature simulation requires higher upper bound
    Cimax(Ca): intercellular_co2_upper_limit => 2Ca ~ track(u"μbar")
    Cimin: intercellular_co2_lower_limit => 0 ~ preserve(u"μbar")
    Ci(Ca, Ci, A_net, gvc): intercellular_co2 => begin
        Ca - Ci ⩵ A_net / gvc
    end ~ bisect(min=Cimin, upper=Cimax, u"μbar")
end

@system IntercellularSpaceDynC3(Weather) begin
    g ~ hold
    Ca ~ hold

    # Use solve for analytical solution of Ci
    Ci_ac(Vcmax, Km, Rd, Γ, g, Ca): intercellular_co2_ac => begin
        a = 1
        b = g^-1 * (Vcmax - Rd) + Km - Ca
        c = -g^-1 * (Vcmax * Γ + Rd * Km) - Ca * Km
        x = Ci_ac
        a * x^2 + b * x + c
    end ~ solve(u"μbar")

    Ci_aj(J, Γ, Rd, g, Ca): intercellular_co2_aj => begin
        a = 1
        b = g^-1 * (J/4 - Rd) + 2Γ - Ca
        c = -J*Γ/4g - 2Rd*Γ/g -2Γ*Ca
        x = Ci_aj
        a * x^2 + b * x + c
    end ~ solve(u"μbar")

    Ci_ap(Tp, Rd, g, Ca): intercellular_co2_ap => begin
        Ca + (Rd - 3Tp) / g
    end ~ track(u"μbar")

    Ci(Ci_ac, Ci_aj, Ci_ap): intercellular_co2 => begin
        max(Ci_ac, Ci_aj, Ci_ap) 
    end ~ track(u"μbar", min=Cimin)

end

@system IntercellularSpaceDynC4(Weather) begin
    g ~ hold
    Ca ~ hold

    # Ci = Cm
    Ci_ac1(Ca, g, Vpmax, Kp, gbs, Rm) => begin
        Ci = Ci_ac1
        # (((Ci * Vpmax) / (Ci + Kp)) + gbs*Ci - Rm) == g * (Ca - Ci)
        a = gbs + g
        b = Kp*gbs - Rm + (Kp - Ca)*g + Vpmax
        c = -(Rm*Kp + g*Kp*Ca)  
        a*Ci^2 + b*Ci + c
    end ~ solve(u"μbar")

    Ci_ac2(Ca, g, Vcmax, Rd) => begin
        Ca - (Vcmax - Rd) / g
    end ~ track(u"μbar")

    Ci_aj1(Ca, g, x, J, Rm, gbs) => begin
        # (x * J/2 - Rm + gbs*Ci) == g * (Ca - Ci)
        (g*Ca - x * J/2 + Rm) / (gbs + g)
    end ~ track(u"μbar")

    Ci_aj2(Ca, g, x, J, Rd) => begin
        Ca - ((1-x) * J/3 - Rd)/g
    end ~ track(u"μbar")

    Ci(Ci_ac1, Ci_ac2, Ci_aj1, Ci_aj2): intercellular_co2 => begin
        max(Ci_ac1, Ci_ac2, Ci_aj1, Ci_aj2)
    end ~ track(u"μbar", min=Cimin)
end

@system IntercellularSpaceDynCAM(Weather) begin
    Ca ~ hold
    g ~ hold

    # Use solve for analytical solution of Ci
    # Ci(Vpmax, Kp, Rdv, fΨv, f_m, gvc, Ca): intercellular_co2 => begin
    #     a = 1
    #     b = fΨv * f_m / gvc * (Vpmax - Rdv) - (Ca - Kp)
    #     c = fΨv * f_m / gvc * Rdv * Kp - Ca * Kp
    #     a * Ci^2 + b * Ci + c
    # end ~ solve(u"μbar")

    ### Asv
    a1(f_c, C0) => f_c * C0 ~ track(u"μbar")
    b1(f_c, C0, Kp) => f_c * C0 + Kp ~ track(u"μbar")
    d1 => 1 ~ preserve
    f1(f_m) => f_m ~ track 
    v1(Vpmax) ~ track(u"μmol/m^2/s")
    R1(Rdv) ~ track(u"μmol/m^2/s")

    ### Asc
    a2(f_c, C0, Γ) => f_c*C0 - Γ ~ track(u"μbar")
    f2(f_c) => (1 - f_c) ~ track
    R2(Rdc) ~ track(u"μmol/m^2/s")
    #Ac
    b2(f_c, C0, Km) => (f_c*C0 + Km) ~ track(u"μbar")
    d2 => 1 ~ preserve
    v2(Vcmax) ~ track(u"μmol/m^2/s")
    #Aj
    b3(f_c, C0, Γ) => (4f_c*C0 + 8Γ) ~ track(u"μbar")
    d3 => 4 ~ preserve
    v3(J) ~ track(u"μmol/m^2/s")

    # xTODO: implement simpler solve for when Asv or Asc = 0
    # causes cyclic dependency because Ci -> Asc and Asc -> flag -> Ci
    # Ci_1() => begin
    #     x = Ci_1
    #     a = g*d
    #     b = f*(v - R*d) - g*(d*Ca - b)
    #     c = f*(v*a - R*b - g*b*Ca)
    #     a*x^2 + b*x + c
    # end ~ solve

    # Ci_ac(Ca, g, a1, b1, d1, f1, v1, R1, a2, b2, d2, f2, v2, R2): intercellular_co2_ac => begin
    #     (v1*(Ci_ac + a1)*f1*(d2 * Ci_ac + b2) 
    #         - R1*f1*(d1*Ci_ac + b1)*(d2*Ci_ac + b2) 
    #         + v2*(Ci_ac+a2)*f2*(d1*Ci_ac + b1) 
    #         - R2*f2*(d1*Ci_ac + b1)*(d2*Ci_ac + b2) 
    #         - g*(Ca - Ci_ac)*(d1*Ci_ac + b1)*(d2*Ci_ac + b2))
    # end ~ solve(u"μbar")
    
    Ci_ac(Ca, g, a1, b1, d1, f1, v1, R1, a2, b2, d2, f2, v2, R2): intercellular_co2_ac => begin
        x = Ci_ac
        a = d1*d2*g
        b = v1*f1*d2 - R1*f1*d1*d2 + v2*f2*d1 - R2*f2*d1*d2 + (d1*b2 + b1*d2 - d1*d2*Ca)*g
        c = v1*f1*(b2 + a1*d2) - R1*f1*(d1*b2 + b1*d2) + v2*f2*(b1 + a2*d1) - R2*f2*(d1*b2 + b1*d2) + (b1*b2 - (d1*b2 + b1*d2)*Ca)*g
        d = v1*a1*f1*b2 - R1*f1*b1*b2 + v2*a2*f2*b1 - R2*f2*b1*b2 - b1*b2*Ca*g
        a*x^3 + b*x^2 + c*x + d
    end ~ solve(u"μbar")

    # Ci_aj(Ca, g, a1, b1, d1, f1, v1, R1, a2, b3, d3, f2, v3, R2): intercellular_co2_aj => begin
    #     (v1*(Ci_aj + a1)*f1*(d3 * Ci_aj + b3) 
    #         - R1*f1*(d1*Ci_aj + b1)*(d3*Ci_aj + b3) 
    #         + v3*(Ci_aj+a2)*f2*(d1*Ci_aj + b1) 
    #         - R2*f2*(d1*Ci_aj + b1)*(d3*Ci_aj + b3) 
    #         - g*(Ca - Ci_aj)*(d1*Ci_aj + b1)*(d3*Ci_aj + b3))
    # end ~ solve(u"μbar")

    Ci_aj(Ca, g, a1, b1, d1, f1, v1, R1, a2, b3, d3, f2, v3, R2): intercellular_co2_aj => begin
        x = Ci_aj
        a = d1*d3*g
        b = v1*f1*d3 - R1*f1*d1*d3 + v3*f2*d1 - R2*f2*d1*d3 + (d1*b3 + b1*d3 - d1*d3*Ca)*g
        c = v1*f1*(b3 + a1*d3) - R1*f1*(d1*b3 + b1*d3) + v3*f2*(b1 + a2*d1) - R2*f2*(d1*b3 + b1*d3) + (b1*b3 - (d1*b3 + b1*d3)*Ca)*g
        d = v1*a1*f1*b3 - R1*f1*b1*b3 + v3*a2*f2*b1 - R2*f2*b1*b3 - b1*b3*Ca*g
        a*x^3 + b*x^2 + c*x + d
    end ~ solve(u"μbar")

    Cimin: intercellular_co2_lower_limit => 0 ~ preserve(u"μbar", parameter)
    Ci(Ci_ac, Ci_aj): intercellular_co2 => begin
        max(Ci_ac, Ci_aj) 
    end ~ track(u"μbar", min=Cimin)

    # Ci(Ca, gvc, a1, b1, d1, f1, v1, R1, a2, b2, d2, f2, v2, R2): intercellular_co2 => begin
    #     v1*(Ci + a1)*f1*(d2 * Ci + b2) 
    #         - R1*f1*(d1*Ci + b1)*(d2*Ci + b2) 
    #         + v2*(Ci+a2)*f2*(d1*Ci + b1) 
    #         - R2*f2*(d1*Ci + b1)*(d2*Ci + b2) 
    #         - gvc*(Ca - Ci)*(d1*Ci + b1)*(d2*Ci + b2)
    # end ~ solve(u"μbar")

    # Ci(AsvFlag, AscFlag, Ca, gvc, a1, b1, d1, f1, v1, R1, a2, b2, d2, f2, v2, R2): intercellular_co2 => begin
    #     if !AsvFlag || !AscFlag
    #         if !AscFlag
    #             v1 * (Ci + a1) - R1 * (d1*Ci + b1) - gvc/f1 *(Ca - Ci)*(d1 * Ci + b1)
    #         else
    #             v2 * (Ci + a2) - R2 * (d2*Ci + b2) - gvc/f2 *(Ca - Ci)*(d2 * Ci + b2)
    #         end
    #     else
    #         v1*(Ci + a1)*f1*(d2 * Ci + b2) 
    #             - R1*f1*(d1*Ci + b1)*(d2*Ci + b2) 
    #             + v2*(Ci+a2)*f2*(d1*Ci + b1) 
    #             - R2*f2*(d1*Ci + b1)*(d2*Ci + b2) 
    #             - gvc*(Ca - Ci)*(d1*Ci + b1)*(d2*Ci + b2)
    #     end
    # end ~ solve(u"μbar")

end