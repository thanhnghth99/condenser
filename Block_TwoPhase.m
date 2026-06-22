function [L_tp, P_out_tp, T_sat, h_sat_l, Q_tp] = Block_TwoPhase(L_sh, P_in_tp, m_ref, T_air_in, m_air_in, h_air, D_h, alpha, W_cond, A_tube_total, A_surface_total, eta_o, G)
    % coder.extrinsic('py.CoolProp.CoolProp.PropsSI');
    Refrig = 'R134a';

    T_sat = py.CoolProp.CoolProp.PropsSI('T', 'P', P_in_tp, 'Q', 1, Refrig);
    h_sat_v = py.CoolProp.CoolProp.PropsSI('H', 'P', P_in_tp, 'Q', 1, Refrig);
    h_sat_l = py.CoolProp.CoolProp.PropsSI('H', 'P', P_in_tp, 'Q', 0, Refrig);
    mu_l = py.CoolProp.CoolProp.PropsSI('V', 'P', P_in_tp, 'Q', 0, Refrig);
    k_l = py.CoolProp.CoolProp.PropsSI('L', 'P', P_in_tp, 'Q', 0, Refrig);
    Pr_l = py.CoolProp.CoolProp.PropsSI('PRANDTL', 'P', P_in_tp, 'Q', 0, Refrig);

    % 1. Required rejected heat transfer rate completely liquid condensation (Q_req)
    Q_req = m_ref * (h_sat_v - h_sat_l);

    % 2. Initial guess for Newton-Raphson method
    L_total = W_cond;
    L_avail = L_total - L_sh; % Available length for the two-phase region

    if L_avail <= 1e-5
        L_tp = 0;
        P_out_tp = P_in_tp;
        Q_tp = 0;
        return;
    end

    % Initial guess for two-phase region length
    L_n = 0.1 * L_avail;
    % Absolute tolerance for heat transfer error [w]
    epsilon_tol = 0.01;
    % Maximun iterations to prevent infinite loops (Failed to convergence)
    iter_max = 25;

    dL = 1e-5;
    % Q_eNTU_n = 0.0;

    %  3. Newton-Raphson loop
    for iter = 1:iter_max
        % Q_eNTU_n from e-NTU method
        Q_eNTU_n = HeatTransfer_eNTU_TP(L_n, P_in_tp, T_sat, m_air_in, T_air_in, mu_l, k_l, Pr_l, D_h, alpha, W_cond, A_tube_total, h_air, A_surface_total, eta_o, G);
        % The error between required heat transfer and e-NTU calculated heat transfer
        f_Ln = Q_eNTU_n - Q_req;

        % Check for convergence
        if abs(f_Ln) < epsilon_tol
            break;
        end

        % Numerical derivative df/dL using central difference
        Ln_plus = L_n + dL;

        % Ensure Ln_plus and Ln_minus are within physical bounds to avoid unphysical results
        if Ln_plus >= L_avail
            Ln_minus = L_n - dL;
            f_Ln_minus = HeatTransfer_eNTU_TP(Ln_minus, P_in_tp, T_sat, m_air_in, T_air_in, mu_l, k_l, Pr_l, D_h, alpha, W_cond, A_tube_total, h_air, A_surface_total, eta_o, G) - Q_req;
            df_dLn = (f_Ln - f_Ln_minus) / dL;
        else
            f_Ln_plus = HeatTransfer_eNTU_TP(Ln_plus, P_in_tp, T_sat, m_air_in, T_air_in, mu_l, k_l, Pr_l, D_h, alpha, W_cond, A_tube_total, h_air, A_surface_total, eta_o, G) - Q_req;
            df_dLn = (f_Ln_plus - f_Ln) / dL;
        end

         % Avoid division by zero (Singularity)
        if abs(df_dLn) < 1e-10
            df_dLn = sign(df_dLn + 1e-20) * 1e-10; % Add a small number to preserve the sign of the derivative
        end

        % Update L_n using Newton-Raphson formula
        L_next = L_n - 0.8 * (f_Ln / df_dLn);

        % Ensure L_next is within physical bounds [0, L_avail]
        if L_next <= 0
            L_next = 1e-5 * L_avail;
        elseif L_next >= L_avail
            L_next = 0.99 * L_avail;
            if iter == iter_max
                L_next = L_avail;
            end
        end
        
        L_n = L_next;
    end

    L_tp = L_n;
    Q_tp = Q_eNTU_n;
    P_out_tp = P_in_tp;
end

% Heat transfer function using e-NTU method
function Q_eNTU = HeatTransfer_eNTU_TP(L_tp, P, T_sat, m_air_in, T_air_in, mu_l, k_l, Pr_l, D_h, alpha, W_cond, A_tube_total, h_air, A_surface_total, eta_o, G)
    if L_tp <= 1e-6
        Q_eNTU = 0;
        return;
    end
    
    % Area fraction of the two-phase region
    w_tp = L_tp / W_cond;

    % Air side
    cp_air = 1005;
    C_air = w_tp * m_air_in * cp_air;

    % Reynolds number assuming all the mass flowing as liquid
    Re_l = G * D_h / mu_l;

    if Re_l < 2300
        % Nusselt number (Shah and London, 1978) for laminar flow in rectangular channels
        Nu_l = 7.541 * (1 - 2.610*alpha + 4.970*alpha^2 - 5.119*alpha^3 + 2.702*alpha^4 - 0.548*alpha^5);
    else
        % Turbulent flow (Dittus-Boelter)
        Nu_l = 0.023 * (Re_l^0.8) * (Pr_l^0.4);
    end

    h_l = Nu_l * (k_l / D_h);

    % Critical pressure of R134a [Pa]
    P_crit = 4.0593e6;
    P_r = P / P_crit;

    % Refrigerant heat transfer coefficient h_tpm [W/m^2.K]
    h_tpm = h_l * (0.55 + 2.09 / P_r^0.38);

    R_ref = 1 / (h_tpm * A_tube_total);
    R_air = 1 / (h_air * A_surface_total * eta_o);
    UA_tp = w_tp / (R_ref + R_air);

    if C_air > 0
        NTU_tp = UA_tp / C_air;
        epsilon = 1 - exp(-NTU_tp);
        Q_eNTU = epsilon * C_air * (T_sat - T_air_in);
    else
        Q_eNTU = 0;
    end
end