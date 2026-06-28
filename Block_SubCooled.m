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

    % P_out_guess = P_in_sc - dP_sc;
    % dP_exit = PressureDrop_exit(P_out_guess, T_out_sc, Refrig, D_h, G, sigma);

    % dP_sc_total = dP_sc + dP_exit;
    % P_out_sc = P_in_sc - dP_sc_total;
    P_out_sc = P_in_sc - dP_sc;
end

% function dP_exit = PressureDrop_exit(P_out_sc, T_out_sc, Refrig, D_h, G, sigma)
    
%     rho_out = py.CoolProp.CoolProp.PropsSI('D', 'P', P_out_sc, 'T', T_out_sc, Refrig);
%     mu_out = py.CoolProp.CoolProp.PropsSI('V', 'P', P_out_sc, 'T', T_out_sc, Refrig);
%     Re_Dh = G * D_h / mu_out;

%     if Re_Dh < 2300
%         % Entrance pressure loss coefficients Kc for  multiple-tube flat-tube core From Kays and London, 1998
%         sigma_lami_array = [0.0034418, 0.013091, 0.026403, 0.03603, 0.050523, 0.063813, 0.077087, 0.089114, 0.10847, 0.12539, 0.14349, 0.16283, 0.187, 0.20272, 0.21963, 0.23897, 0.25589, 0.28009, 0.30066, 0.32123, 0.33814, 0.35993, 0.38051, 0.39623, 0.41922, 0.43617, 0.46401, 0.49309, 0.51611, 0.53793, 0.56219, 0.57915, 0.61308, 0.63856, 0.66888, 0.70043, 0.72227, 0.74777, 0.77692, 0.8097, 0.83276, 0.86801, 0.89595, 0.9239, 0.95914, 0.98953];
%         Ke_lami_array = [0.99813, 0.97407, 0.95183, 0.92223, 0.89076, 0.86298, 0.83152, 0.79268, 0.75932, 0.72598, 0.68156, 0.64452, 0.59452, 0.56488, 0.52969, 0.4908, 0.45746, 0.41669, 0.37963, 0.34626, 0.30923, 0.27401, 0.24249, 0.21285, 0.17209, 0.14613, 0.10163, 0.062667, 0.031126, 0.0014401, -0.02273, -0.048689, -0.093228, -0.11925, -0.15454, -0.18614, -0.20845, -0.22894, -0.25315, -0.28107, -0.30339, -0.3221, -0.34261, -0.36128, -0.38368, -0.40237];
        
%         % Interpolation for Kc corresponding sigma
%         Ke_interp = interp1(sigma_lami_array, Ke_lami_array, sigma, 'linear');
%     else
%         sigma_turb_array = [0.0034269, 0.01795, 0.030021, 0.048176, 0.051749, 0.073537, 0.090474, 0.11101, 0.12309, 0.14487, 0.16182, 0.18723, 0.21265, 0.23931, 0.27684, 0.30347, 0.31557, 0.34344, 0.36888, 0.40283, 0.43189, 0.47437, 0.51074, 0.54714, 0.58476, 0.61271, 0.65158, 0.69167, 0.73664, 0.77552, 0.82417, 0.86311, 0.90325, 0.94465, 0.97754, 0.99825];
%         Ke_turb_array = [0.99444, 0.97034, 0.94257, 0.91292, 0.89259, 0.85737, 0.82772, 0.78513, 0.75921, 0.7203, 0.69434, 0.65171, 0.61093, 0.57567, 0.51634, 0.47186, 0.45147, 0.41436, 0.37727, 0.33826, 0.29376, 0.25285, 0.20829, 0.16927, 0.13393, 0.11526, 0.087291, 0.057466, 0.029448, 0.0033222, -0.017343, -0.030555, -0.04562, -0.05147, -0.055411, -0.055568];

%         % Interpolation for Kc corresponding sigma
%         Ke_interp = interp1(sigma_turb_array, Ke_turb_array, sigma, 'linear');
%     end

%     % 4. Tính tổn thất áp suất Cửa ra theo Kays & London
%     % LƯU Ý: dP_exit là áp suất BỊ MẤT. Vì sự mở rộng làm TĂNG áp suất, 
%     % nên giá trị mất mát này thường mang dấu ÂM (Pressure Recovery).
%     dP_exit = - (G^2 / (2*rho_out)) * (1 - sigma^2 - Ke_interp);
% end
