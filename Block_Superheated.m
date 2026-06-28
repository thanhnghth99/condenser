function [L_sh, P_out_sh, dP_sh_total, dP_entrance, T_sat_v, h_sat_v, Q_sh] = Block_Superheated(P_ref_in, T_ref_in, m_ref, T_air_in, m_air_in, h_air, D_h, alpha, W_cond, A_tube_total, A_surface_total, eta_o, G, sigma)
    % coder.extrinsic('py.CoolProp.CoolProp.PropsSI');
    Refrig = 'R134a';
    
    % CoolProp library
    h_sat_v = py.CoolProp.CoolProp.PropsSI('H', 'P', P_ref_in, 'Q', 1, Refrig);
    h_ref_in = py.CoolProp.CoolProp.PropsSI('H', 'P', P_ref_in, 'T', T_ref_in, Refrig);
    T_sat_v = py.CoolProp.CoolProp.PropsSI('T', 'P', P_ref_in, 'Q', 1, Refrig);

    T_rm_sh = (T_ref_in + T_sat_v) / 2;
    cp_v = py.CoolProp.CoolProp.PropsSI('C', 'P', P_ref_in, 'T', T_rm_sh, Refrig);
    mu_v = py.CoolProp.CoolProp.PropsSI('V', 'P', P_ref_in, 'T', T_rm_sh, Refrig);
    k_v = py.CoolProp.CoolProp.PropsSI('L', 'P', P_ref_in, 'T', T_rm_sh, Refrig);
    rho_v = py.CoolProp.CoolProp.PropsSI('D', 'P', P_ref_in, 'T', T_rm_sh, Refrig);
    Pr_v = py.CoolProp.CoolProp.PropsSI('PRANDTL', 'P', P_ref_in, 'T', T_rm_sh, Refrig);

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
    dP_sh = PressureDrop_SH(L_sh, rho_v, mu_v, D_h, alpha, G);
    
    % 5. Calculate minor losses at entrance
    dP_entrance = PressureDrop_entrance(P_ref_in, T_ref_in, Refrig, D_h, G, sigma);

    % 6. Total pressure drop of superheated region
    dP_sh_total = dP_entrance + dP_sh;
    P_out_sh = P_ref_in - dP_sh - dP_entrance;
end

function dP_entrance = PressureDrop_entrance(P_ref_in, T_ref_in, Refrig, D_h, G, sigma)
    rho_in = py.CoolProp.CoolProp.PropsSI('D', 'P', P_ref_in, 'T', T_ref_in, Refrig);
    mu_in = py.CoolProp.CoolProp.PropsSI('V', 'P', P_ref_in, 'T', T_ref_in, Refrig);
    
    Re_Dh = G * D_h / mu_in;

    if Re_Dh < 2300
        % Entrance pressure loss coefficients Kc for  multiple-tube flat-tube core From Kays and London, 1998
        sigma_lami_array = [0.017172, 0.063467, 0.12194, 0.17791, 0.2363, 0.30077, 0.35295, 0.40883, 0.46953, 0.51315, 0.54591, 0.581, 0.6149, 0.66586, 0.70447, 0.73717, 0.76859, 0.79638, 0.82541, 0.86284, 0.89911, 0.93416, 0.95944, 0.99205, 1.0077];
        Kc_lami_array = [0.79657, 0.79621, 0.79576, 0.78974, 0.78183, 0.77387, 0.75668, 0.7432, 0.72595, 0.70696, 0.69552, 0.6766, 0.65955, 0.64237, 0.61223, 0.59706, 0.5763, 0.55744, 0.54043, 0.51403, 0.49137, 0.46872, 0.44428, 0.42165, 0.40661];
        
        % Interpolation for Kc corresponding sigma
        Kc_interp = interp1(sigma_lami_array, Kc_lami_array, sigma, 'linear');
    else
        sigma_turb_array = [0.0092614, 0.043363, 0.085993, 0.11889, 0.16392, 0.21384, 0.2455, 0.28929, 0.33433, 0.37686, 0.4231, 0.46929, 0.49722, 0.53856, 0.57622, 0.60415, 0.64302, 0.67945, 0.71468, 0.74742, 0.78137, 0.81418, 0.85904, 0.90151, 0.93426, 0.95366, 0.98156, 1.0034];
        Kc_turb_array = [0.46548, 0.46337, 0.46119, 0.46094, 0.4532, 0.44727, 0.44518, 0.43375, 0.42785, 0.40903, 0.39943, 0.38243, 0.36742, 0.35416, 0.33722, 0.32221, 0.30342, 0.28464, 0.26773, 0.24343, 0.21728, 0.20408, 0.16859, 0.14237, 0.11808, 0.10313, 0.084422, 0.065758];

        % Interpolation for Kc corresponding sigma
        Kc_interp = interp1(sigma_turb_array, Kc_turb_array, sigma, 'linear');
    end

    % Pressure drop due to minor losses
    dP_entrance = (G^2 / (2*rho_in)) * (1 - sigma^2 + Kc_interp);
end

% Pressure drop function
function dP_sh = PressureDrop_SH(L_sh, rho_v, mu_v, D_h, alpha, G)
    if L_sh <= 1e-6
        dP_sh = 0;
        return;
    end

    Re_Dh = G * D_h / mu_v;

    if Re_Dh < 2300
        % Hydrodynamic entrance region for noncircular duct
        x_plus = L_sh / (D_h * Re_Dh);
        x_plus = max(x_plus, 1e-10); % Avoid division-by-zero singularity

        % Fully developed Fanning friction factor
        fRe_fd = 24 * (1 - 1.3553*alpha + 1.9467*alpha^2 - 1.7012*alpha^3 + 0.9564*alpha^4 - 0.2537*alpha^5);

        % Incremental pressure drop number K(infinity)
        K_inf = 0.6796 + 1.2197*alpha + 3.3089*alpha^2 - 9.5921*alpha^3 + 8.9089*alpha^4 - 2.9959*alpha^5;

        % Apparent Fanning friction factor times Re correlation
        fappRe = (3.44 / sqrt(x_plus)) + (fRe_fd + (K_inf / (4 * x_plus)) - (3.44 / sqrt(x_plus))) / (1 + 0.0002367 / (x_plus^2));

        % Convert Apparent Fanning to Darcy friction factor (f_darcy = 4 * f_fanning)
        f_D = 4 * (fappRe / Re_Dh);
    else
        % Turbulent flow (Darcy friction factor via Petukhov)
        f_D = 1 / (0.79 * log(Re_Dh) - 1.64)^2; 
    end

    dP_sh = f_D * L_sh * (G^2) / (2 * D_h * rho_v);
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
        % Fully developed Nusselt number
        Nu_fd = 7.541 * (1 - 2.610*alpha + 4.970*alpha^2 - 5.119*alpha^3 + 2.702*alpha^4 - 0.548*alpha^5);

        % Laminar flow in rectangular flat tube (Shah & London)
        x_star = L_sh / (D_h * Re_Dh * Pr_v);
        x_star = max(x_star, 1e-10); % Avoid singularity

        % Stephan correlation for developing Nusselt number
        Nu_sh = Nu_fd + 0.0668 / (x_star^(1/3) * (0.04 + x_star^(2/3)));
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

