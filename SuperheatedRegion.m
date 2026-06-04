classdef SuperheatedRegion
    properties
        % Reference to the condenser model to access inlet conditions and geometry
        Model
    end

    methods
        function obj = SuperheatedRegion(condenser_model)
            obj.Model = condenser_model;
        end
        
        function [L_sh, P_out_sh, dP_sh, T_sat_v, h_v, Q_eNTU_sh] = defineRegion(obj)

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
            iter_max = 50;
            
            dL = 1e-6;
            converged = false;

            % 3. Newton-Raphson loop
            for iter = 1:iter_max
                fprintf('Iteration %d: L_sh = %.6f m\n', iter, L_n);

                % Q_eNTU_sh from e-NTU method
                Q_eNTU_sh = obj.HeatTransfer_eNTU(L_n, props);

                % The error between required heat transfer and e-NTU calculated heat transfer
                f_Ln = Q_eNTU_sh - Q_req;

                % Check for convergence
                if abs(f_Ln) < epsilon_tol
                    converged = true;
                    break;
                end

                % Numerical derivative df/dL using central difference
                f_Ln_plus = obj.HeatTransfer_eNTU(L_n + dL, props) - Q_req;
                df_dLn = (f_Ln_plus - f_Ln) / dL;

                % Avoid division by zero (Singularity)
                if abs(df_dLn) < 1e-10
                    df_dLn = sign(df_dLn) * 1e-10;
                end

                % Update L_n using Newton-Raphson formula
                L_next = L_n - f_Ln / df_dLn;

                % Ensure L_next is within physical bounds [0, L_total]
                if L_next < 0
                    L_next = 0;
                elseif L_next > L_total
                    warning('Superheated region length exceeded total condenser length. Stopping iteration.');
                    break;
                end
                L_n = L_next;
            end

            if iter == iter_max
                % Warning if the loop reachs the iteration limit without converging
                warning('Superheated region FAILED to converge!');                
            end

            L_sh = L_n;
            % Calculate outlet pressure and pressure drop
            rho_sh = props.rho;
            mu_sh = props.mu;
            dP_sh = obj.PressureDropSH(L_sh, rho_sh, mu_sh);
            P_out_sh = P_ref_in - dP_sh;
            % Saturated vapor temperature and enthalpy at outlet pressure
            T_sat_v = ThermoProp.get_T_sat(P_out_sh, 1, Refrig);
            h_v = ThermoProp.get_SatVaporProps(P_out_sh, Refrig).h_v;
        end
    end

    methods (Access = private)

        % Pressure drop
        function dP_sh = PressureDropSH(obj, L_sh, rho, mu)
            % Reynolds number
            Re_Dh = obj.Model.G * obj.Model.D_h / mu;

            % Friction factor f
            f = 1 / (1.58 * log(Re_Dh) - 3.28)^2; % Turbulent flow

            % Pressure drop dP_sh
            dP_sh = f * L_sh * (obj.Model.G^2) / (2 * obj.Model.D_h * rho);
        end

        % Heat transfer using e-NTU method
        function Q_eNTU = HeatTransfer_eNTU(obj, L_sh, props)
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