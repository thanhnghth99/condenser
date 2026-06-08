classdef SuperheatedRegion
    properties
        % Reference to the condenser model to access inlet conditions and geometry
        Model
    end

    methods
        function obj = SuperheatedRegion(condenser_model)
            obj.Model = condenser_model;
        end
        
        function [L_sh, P_out_sh, dP_sh, T_sat_v, h_sat_v, Q_eNTU_sh] = defineRegion(obj)

            % Refrigerant inlet conditions
            P_ref_in = obj.Model.Inlet.P_ref_in;
            T_ref_in = obj.Model.Inlet.T_ref_in;
            h_ref_in = obj.Model.Inlet.h_ref_in;
            m_ref = obj.Model.Inlet.m_ref;
            Refrig = obj.Model.Inlet.Refrigerant;

            h_sat_v = ThermoProp.get_SatVaporProps(P_ref_in, Refrig).h_v;

            props = ThermoProp.get_SinglePhaseProps(P_ref_in, T_ref_in, Refrig);

            % 1. Required rejected heat transfer rate to saturated vapor (Q_req)
            Q_req = m_ref * (h_ref_in - h_sat_v);

            % 2. Initial guess for Newton-Raphson method
            L_total = obj.Model.W_cond;
            % Initial guess for superheated region length
            L_n = 0.1 * L_total;
            % Absolute tolerance for heat transfer error [w]
            epsilon_tol = 0.01;
            % Maximun iterations to prevent infinite loops (Failed to convergence)
            iter_max = 10;
            
            dL = 1e-5;
            converged = false;

            % 3. Newton-Raphson loop
            for iter = 1:iter_max
                fprintf('Iteration %d: L_sh = %.6f m\n', iter, L_n);

                % Q_eNTU_n from e-NTU method
                Q_eNTU_n = obj.HeatTransfer_eNTU(L_n, props);

                % The error between required heat transfer and e-NTU calculated heat transfer
                f_Ln = Q_eNTU_n - Q_req;

                % Check for convergence
                if abs(f_Ln) < epsilon_tol
                    converged = true;
                    break;
                end

                % Numerical derivative df/dL using central difference
                Ln_plus = L_n + dL;
                % Ensure Ln_plus does not exceed total condenser length to avoid unphysical results
                if Ln_plus >= L_total
                    % If Ln_plus exceeds total length, use backward difference instead
                    Ln_minus = L_n - dL;
                    f_Ln_minus = obj.HeatTransfer_eNTU(Ln_minus, props) - Q_req;
                    df_dLn = (f_Ln - f_Ln_minus) / dL;
                else
                    f_Ln_plus = obj.HeatTransfer_eNTU(Ln_plus, props) - Q_req;
                    df_dLn = (f_Ln_plus - f_Ln) / dL;
                end

                % Avoid division by zero (Singularity)
                if abs(df_dLn) < 1e-10
                    df_dLn = sign(df_dLn + 1e-20) * 1e-10;
                end

                % Update L_n using Newton-Raphson formula
                L_next = L_n - f_Ln / df_dLn;

                % Ensure L_next is within physical bounds [0, L_total]
                if L_next <= 0
                    warning('Superheated region length became negative. Setting to a small positive value.');
                    L_next = 1e-5 * L_total; % Set to a small positive value to avoid zero length
                elseif L_next >= L_total
                    warning('Superheated region length exceeded total condenser length. Setting to maximum possible length.');
                    L_next = 0.99 * L_total;
                    if iter == iter_max
                        L_next = L_total;
                    end
                end
                L_n = L_next;
            end

            if ~converged
                % Warning if the loop reachs the iteration limit without converging
                warning('Superheated region FAILED to converge after %d iterations!', iter_max);                
            end

            L_sh = L_n;
            Q_eNTU_sh = Q_eNTU_n;

            % 4. Calculate outlet pressure and pressure drop
            rho_sh = props.rho;
            mu_sh = props.mu;
            dP_sh = obj.PressureDropSH(L_sh, rho_sh, mu_sh);
            P_out_sh = P_ref_in - dP_sh;

            % 5. Saturated vapor temperature and enthalpy at outlet pressure
            T_sat_v = ThermoProp.get_T_sat(P_out_sh, 1, Refrig);
            h_sat_v = ThermoProp.get_SatVaporProps(P_out_sh, Refrig).h_v;
        end
    end

    methods (Access = private)

        % Pressure drop
        function dP_sh = PressureDropSH(obj, L_sh, rho, mu)
            % Reynolds number
            Re_Dh = obj.Model.G * obj.Model.D_h / mu;

            % % Friction factor f
            % f = 1 / (1.58 * log(Re_Dh) - 3.28)^2; % Turbulent flow
            % fprintf('f = %.4f\n', f)

            if Re_Dh < 2300
                % Laminar flow in rectangular flat tube
                alpha = obj.Model.alpha;

                f = 4 * (24 / Re_Dh) * (1 - 1.3553*alpha + 1.9467*alpha^2 - 1.7012*alpha^3 + 0.9564*alpha^4 - 0.2537*alpha^5);
            else
                f = 1 / (1.58 * log(Re_Dh) - 3.28)^2; % Turbulent flow
            end

            % Pressure drop dP_sh
            dP_sh = f * L_sh * (obj.Model.G^2) / (2 * obj.Model.D_h * rho);
        end

        % Heat transfer using e-NTU method
        function Q_eNTU = HeatTransfer_eNTU(obj, L_sh, props)
            if L_sh <= 1e-6
                Q_eNTU = 0;
                return;
            end

            w_sh = L_sh / obj.Model.W_cond; % Area fraction of superheated region

            % Air side
            % Specific heat capacity of air [J/kg.K]
            cp_air = 1005;
            T_air_in = obj.Model.Inlet.T_air_in;
            m_air_in = obj.Model.Inlet.m_air_in;
            C_air = w_sh * m_air_in * cp_air;

            % Refrigerant side
            % Specific heat capacity of refrigerant
            T_ref_in = obj.Model.Inlet.T_ref_in;
            m_ref = obj.Model.Inlet.m_ref;
            C_r134a = m_ref * props.cp;

            % Determine C_min, C_max and specific heat capacity ratio C_r
            C_min = min(C_air, C_r134a);
            C_max = max(C_air, C_r134a);
            C_r = C_min / C_max;

            % Reynolds number
            D_h = obj.Model.D_h;
            Re_Dh = obj.Model.G * D_h / props.mu;

            % Prandtl number
            Pr_sh = props.Pr;

            % Thermal conductivity
            k_sh = props.k;

            % Friction factor f
            f = 1 / (1.58 * log(Re_Dh) - 3.28)^2;
            % Refrigerant heat transfer coefficient h_sh [W/m^2.K]
            numerator = (f / 2) * Re_Dh * Pr_sh;
            denominator = 1.07 + 12.7 * (f / 2)^0.5 * (Pr_sh^(2/3) - 1);
            h_sh = (numerator / denominator) * (k_sh / D_h);


            UA_sh = obj.Model.UA(w_sh, h_sh);
            NTU_sh = UA_sh / C_min;
            epsilon = 1 - exp((NTU_sh^0.22 / C_r) * (exp(-C_r * NTU_sh^0.78) - 1));

            Q_eNTU = epsilon * C_min * (T_ref_in - T_air_in);
        end
    end
end