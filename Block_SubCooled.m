function [L_sc, P_out_sc, T_out_sc, Q_sc, dP_sc, dT_sc] = Block_SubCooled(L_sh, L_tp, P_in_sc, T_sat, m_ref, T_air_in, m_air_in, h_air, D_h, alpha, W_cond, A_tube_total, A_surface_total, eta_o, G)
    % coder.extrinsic('py.CoolProp.CoolProp.PropsSI');
    Refrig = 'R134a';

    L_total = W_cond;
    L_sc = L_total - L_sh - L_tp;

    % Check if no length left for subcooling, exit immediately
    if L_sc <= 1e-5
        Q_sc = 0;
        T_out_sc = T_sat;
        P_out_sc = P_in_sc;
        dP_sc = 0;
        dT_sc = 0;
        L_sc = 0;
        return;
    end

    w_sc = L_sc / L_total;

    mu_l = py.CoolProp.CoolProp.PropsSI('V', 'P', P_in_sc, 'Q', 0, Refrig);
    k_l = py.CoolProp.CoolProp.PropsSI('L', 'P', P_in_sc, 'Q', 0, Refrig);
    Pr_l = py.CoolProp.CoolProp.PropsSI('PRANDTL', 'P', P_in_sc, 'Q', 0, Refrig);
    rho_l = py.CoolProp.CoolProp.PropsSI('D', 'P', P_in_sc, 'Q', 0, Refrig);
    
    % Retrieve cp from single-phase properties slightly below T_sat to avoid singularity
    cp_l = py.CoolProp.CoolProp.PropsSI('C', 'P', P_in_sc, 'T', T_sat - 0.5, Refrig);

    Re_l = G * D_h / mu_l;

    % Q_eNTU_n from e-NTU method
    cp_air = 1005;
    C_air = w_sc * m_air_in * cp_air;
    C_r134a = m_ref * cp_l;
    C_min = min(C_air, C_r134a);
    C_max = max(C_air, C_r134a);
    C_r = C_min / C_max;

    if Re_l < 2300
        Nu_l = 7.541 * (1 - 2.610*alpha + 4.970*alpha^2 - 5.119*alpha^3 + 2.702*alpha^4 - 0.548*alpha^5);
    else
        Nu_l = 0.023 * (Re_l^0.8) * (Pr_l^0.4);
    end
    h_sc = Nu_l * k_l / D_h;

    R_ref = 1 / (h_sc * A_tube_total);
    R_air = 1 / (h_air * A_surface_total * eta_o);
    UA_sc = w_sc / (R_ref + R_air);

    if C_min > 0
        NTU_sc = UA_sc / C_min;
        epsilon = 1 - exp((NTU_sc^0.22 / C_r) * (exp(-C_r * NTU_sc^0.78) - 1));
        Q_sc = epsilon * C_min * (T_sat - T_air_in);
    else
        Q_sc = 0;
    end

    % Calculate final thermodynamic states
    T_out_sc = T_sat - Q_sc / (m_ref * cp_l);
    dT_sc = T_sat - T_out_sc;

    % 3. Pressure drop (Using Fanning friction factor)
    if Re_l < 2300
        f_F = (24 / Re_l) * (1 - 1.3553*alpha + 1.9467*alpha^2 - 1.7012*alpha^3 + 0.9564*alpha^4 - 0.2537*alpha^5);
    else
        f_F = 0.079 * Re_l^(-0.25);
    end
    dP_sc = 2 * L_sc * f_F * G^2 / (rho_l * D_h);
    P_out_sc = P_in_sc - dP_sc;
end

