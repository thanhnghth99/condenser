function [L_sh, P_out_sh, dP_sh, T_sat_v, h_sat_v, Q_sh] = Block_Superheated(P_ref_in, T_ref_in, m_ref, T_air_in, m_air_in, h_air, D_h, alpha, W_cond, A_tube_total, A_surface_total, eta_o, G)
    % coder.extrinsic('py.CoolProp.CoolProp.PropsSI');
    Refrig = 'R134a';
    
    % CoolProp library
    h_sat_v = py.CoolProp.CoolProp.PropsSI('H', 'P', P_ref_in, 'Q', 1, Refrig);
    h_ref_in = py.CoolProp.CoolProp.PropsSI('H', 'P', P_ref_in, 'T', T_ref_in, Refrig);
    cp_v = py.CoolProp.CoolProp.PropsSI('C', 'P', P_ref_in, 'T', T_ref_in, Refrig);
    mu_v = py.CoolProp.CoolProp.PropsSI('V', 'P', P_ref_in, 'T', T_ref_in, Refrig);
    k_v = py.CoolProp.CoolProp.PropsSI('L', 'P', P_ref_in, 'T', T_ref_in, Refrig);
    rho_v = py.CoolProp.CoolProp.PropsSI('D', 'P', P_ref_in, 'T', T_ref_in, Refrig);
    Pr_v = py.CoolProp.CoolProp.PropsSI('PRANDTL', 'P', P_ref_in, 'T', T_ref_in, Refrig);
    T_sat_v = py.CoolProp.CoolProp.PropsSI('T', 'P', P_ref_in, 'Q', 1, Refrig);

    % 1. Required rejected heat transfer rate to saturated vapor (Q_req)
    Q_req = m_ref * (h_ref_in - h_sat_v);

    % 2. Initial guess for Newton-Raphson method
    % Initial guess for Newton-Raphson method
    L_total = W_cond;
    % Initial guess for superheated region length
    L_n = 0.1 * L_total;
    % Absolute tolerance for heat transfer error [w]
    epsilon_tol = 0.01;
    % Maximun iterations to prevent infinite loops (Failed to convergence)
    iter_max = 25;

    dL = 1e-5;
    Q_eNTU_n = 0.0;

    % 3. Newton-Raphson loop
    for iter = 1:iter_max
        % Q_eNTU_n from e-NTU method
        Q_eNTU_n = HeatTransfer_eNTU_SH(L_n, m_ref, T_ref_in, T_air_in, m_air_in, cp_v, mu_v, k_v, Pr_v, D_h, alpha, W_cond, A_tube_total, h_air, A_surface_total, eta_o, G);
        
        % The error between required heat transfer and e-NTU calculated heat transfer
        f_Ln = Q_eNTU_n - Q_req;

        % Check for convergence
        if abs(f_Ln) < epsilon_tol
            break;
        end

        % Numerical derivative df/dL using central difference
        Ln_plus = L_n + dL;
        % Ensure Ln_plus does not exceed total condenser length to avoid unphysical results
        if Ln_plus >= L_total
            % If Ln_plus exceeds total length, use backward difference instead
            Ln_minus = L_n - dL;
            f_Ln_minus = HeatTransfer_eNTU_SH(Ln_minus, m_ref, T_ref_in, T_air_in, m_air_in, cp_v, mu_v, k_v, Pr_v, D_h, alpha, W_cond, A_tube_total, h_air, A_surface_total, eta_o, G) - Q_req;
            df_dLn = (f_Ln - f_Ln_minus) / dL;
        else
            f_Ln_plus = HeatTransfer_eNTU_SH(Ln_plus, m_ref, T_ref_in, T_air_in, m_air_in, cp_v, mu_v, k_v, Pr_v, D_h, alpha, W_cond, A_tube_total, h_air, A_surface_total, eta_o, G) - Q_req;
            df_dLn = (f_Ln_plus - f_Ln) / dL;
        end

        % Avoid division by zero (Singularity)
        if abs(df_dLn) < 1e-10
            df_dLn = sign(df_dLn + 1e-20) * 1e-10;
        end

        % Update L_n using Newton-Raphson formula
        L_next = L_n - 0.8 * (f_Ln / df_dLn); % Relaxation factor = 0.8

        % Ensure L_next is within physical bounds [0, L_total]
        if L_next <= 0
            L_next = 1e-5 * L_total;
        elseif L_next >= L_total
            L_next = 0.99 * L_total;
            if iter == iter_max
                L_next = L_total;
            end
        end
        L_n = L_next;
    end

    L_sh = L_n;
    Q_sh = Q_eNTU_n;

    % 4. Calculate outlet pressure and pressure drop
    Re_Dh = G * D_h / mu_v;
    if Re_Dh < 2300
        f_D = 4 * (24 / Re_Dh) * (1 - 1.3553*alpha + 1.9467*alpha^2 - 1.7012*alpha^3 + 0.9564*alpha^4 - 0.2537*alpha^5);
    else
        f_D = 1 / (0.79 * log(Re_Dh) - 1.64)^2;
    end
    dP_sh = f_D * L_sh * (G^2) / (2 * D_h * rho_v);
    P_out_sh = P_ref_in - dP_sh;
end

% Heat transfer function using e-NTU method
function Q_eNTU = HeatTransfer_eNTU_SH(L_sh, m_ref, T_ref_in, T_air_in, m_air_in, cp_v, mu_v, k_v, Pr_v, D_h, alpha, W_cond, A_tube_total, h_air, A_surface_total, eta_o, G)
    if L_sh <= 1e-6
        Q_eNTU = 0;
        return;
    end

    % Area fraction of superheated region
    w_sh = L_sh / W_cond;

    % Air side
    cp_air = 1005;
    C_air = w_sh * m_air_in * cp_air;
    C_r134a = m_ref * cp_v;
    C_min = min(C_air, C_r134a);
    C_max = max(C_air, C_r134a);
    C_r = C_min / C_max;

    Re_Dh = G * D_h / mu_v;

    if Re_Dh < 2300
        % Laminar flow in rectangular flat tube (Shah & London)
        Nu_sh = 7.541 * (1 - 2.610*alpha + 4.970*alpha^2 - 5.119*alpha^3 + 2.702*alpha^4 - 0.548*alpha^5);
        h_sh = Nu_sh * (k_v / D_h);
    else
        % Turbulent flow (Darcy friction factor)
        f_D = 1 / (0.79 * log(Re_Dh) - 1.64)^2;
        numerator = (f_D / 8) * Re_Dh * Pr_v;
        denominator = 1.07 + 12.7 * (f_D / 8)^0.5 * (Pr_v^(2/3) - 1);
        h_sh = (numerator / denominator) * (k_v / D_h);
    end

    R_ref = 1 / (h_sh * A_tube_total);
    R_air = 1 / (h_air * A_surface_total * eta_o);
    UA_sh = w_sh / (R_ref + R_air);

    if C_min > 0
        NTU_sh = UA_sh / C_min;
        epsilon = 1 - exp((NTU_sh^0.22 / C_r) * (exp(-C_r * NTU_sh^0.78) - 1));
        Q_eNTU = epsilon * C_min * (T_ref_in - T_air_in);
    else
        Q_eNTU = 0;
    end
end

