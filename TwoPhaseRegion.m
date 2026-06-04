classdef TwoPhaseRegion
    properties
        % Reference to the condenser model to access inlet conditions and geometry
        Model
    end

    methods
        % Constructor to initialize the TwoPhaseRegion
        function obj = TwoPhaseRegion(condenser_model)
            obj.Model = condenser_model;
        end

        % Function to define the two-phase region and calculate the area fraction
        function [L_tp, P_out_tp, T_sat, h_sat_v, h_sat_l, Q_eNTU_tp] = defineRegion(obj, L_sh, P_out_sh)
            P_in_tp = P_out_sh; % Initial pressure at the inlet of the two-phase region (outlet pressure of SH region)

            m_ref = obj.Model.Inlet.m_ref;
            Refrig = obj.Model.Inlet.Refrigerant;

            T_sat = ThermoProp.get_T_sat(P_in_tp, 1, Refrig);

            h_sat_v = ThermoProp.get_SatVaporProps(P_in_tp, Refrig).h_v;
            fprintf('h_sat_v: %.2f J/kg\n', h_sat_v);
            h_sat_l = ThermoProp.get_SatLiquidProps(P_in_tp, Refrig).h_l;
            fprintf('h_sat_l: %.2f J/kg\n\n', h_sat_l);

            props = ThermoProp.get_SatLiquidProps(P_in_tp, Refrig);

            % 1. Required rejected heat transfer rate completely liquid condensation (Q_req)
            Q_req = m_ref * (h_sat_v- h_sat_l);
            fprintf('Q_req: %.2f W\n\n', Q_req);

            % 2. Initial guess for Newton-Raphson method
            L_total = obj.Model.W_cond;
            L_avail = L_total - L_sh; % Available length for the two-phase region

            if L_avail <= 0
                error('No available length for two-phase region. Superheated region occupies the entire condenser length.');
            end

            % Initial guess for two-phase region length
            L_n = 0.1 * L_avail;
            % Absolute tolerance for heat transfer error [w]
            epsilon_tol = 0.01;
            % Maximun iterations to prevent infinite loops (Failed to convergence)
            iter_max = 50;
            
            dL = 1e-5;
            converged = false;

            %  3. Newton-Raphson loop
            for iter = 1:iter_max
                fprintf('Iteration %d: L_n = %.6f m\n', iter, L_n);

                % Q_eNTU_n from e-NTU method
                Q_eNTU_n = obj.HeatTransfer_eNTU(L_sh, L_n, P_in_tp, T_sat, props);
                fprintf('Q_eNTU_n: %.2f W\n', Q_eNTU_n);

                % The error between required heat transfer and e-NTU calculated heat transfer
                f_Ln = Q_eNTU_n - Q_req;
                fprintf('f_Ln: %.2f W\n', f_Ln);

                % Check for convergence
                if abs(f_Ln) < epsilon_tol
                    converged = true;
                    break;
                end

                % Numerical derivative df/dL using central difference
                Ln_plus = L_n + dL;
                fprintf('Ln_plus: %.6f m\n', Ln_plus);

                % Ensure Ln_plus and Ln_minus are within physical bounds to avoid unphysical results
                if Ln_plus >= L_avail
                    Ln_minus = L_n - dL;
                    f_Ln_minus = obj.HeatTransfer_eNTU(L_sh, Ln_minus, P_in_tp, T_sat, props) - Q_req;
                    df_dLn = (f_Ln - f_Ln_minus) / dL;
                    fprintf('Using BACKWARD difference: Ln_minus = %.6f m, f_Ln_minus = %.2f W\n', Ln_minus, f_Ln_minus);
                    fprintf('df_dLn: %.2f W/m\n', df_dLn);
                else
                    f_Ln_plus = obj.HeatTransfer_eNTU(L_sh, Ln_plus, P_in_tp, T_sat, props) - Q_req;
                    df_dLn = (f_Ln_plus - f_Ln) / dL;
                    fprintf('Using FORWARD difference: Ln_plus = %.6f m, f_Ln_plus = %.2f W\n', Ln_plus, f_Ln_plus);
                    fprintf('df_dLn: %.2f W/m\n', df_dLn);
                end

                % Avoid division by zero (Singularity)
                if abs(df_dLn) < 1e-10
                    df_dLn = sign(df_dLn + 1e-20) * 1e-10; % Add a small number to preserve the sign of the derivative
                end

                % Update L_n using Newton-Raphson formula
                L_next = L_n - f_Ln / df_dLn;
                fprintf('L_next (before bounds check): %.6f m\n', L_next);

                % step = f_Ln / df_dLn;
                % inner_iter = 0;
                % while (L_next <= 0 || L_next >= L_avail) && (inner_iter < 10)
                %     step = step / 2;
                %     L_next = L_n - step;
                %     inner_iter = inner_iter + 1;
                % end

                % Ensure L_next is within physical bounds [0, L_avail]
                if L_next <= 0
                    warning('Two-phase region length became negative. Setting to a small positive value.');
                    L_next = 1e-5 * L_avail; % Set to a small positive value to avoid zero length
                elseif L_next >= L_avail
                    warning('Two-phase region length exceeded total condenser length. Setting to maximum possible length.');
                    L_next = 0.99 * L_avail;
                end

                L_n = L_next;
            end

            if ~converged
                % Warning if the loop reachs the iteration limit without converging
                warning('Two-phase region FAILED to converge after %d iterations!', iter_max);                
            end

            L_tp = L_n;
            Q_eNTU_tp = Q_eNTU_n;

            P_out_tp = P_in_tp;
        end
    end

    methods (Access = private)
        % Function to calculate the actual heat transfer capacity (Q_eNTU) based on the NTU method for 1 element
        function Q_eNTU = HeatTransfer_eNTU(obj, L_sh, L_tp, P, T_sat, props)
            if L_tp <= 1e-6
                Q_eNTU = 0;
                return;
            end

            w_tp = L_tp / (obj.Model.W_cond - L_sh); % Area fraction of the two-phase region

            % Air side
            % Specific heat capacity of air [J/kg.K] for 1 element
            cp_air = 1005;
            T_air_in = obj.Model.Inlet.T_air_in;
            m_air = obj.Model.Inlet.m_air_in;
            C_air = w_tp * m_air * cp_air;

            % Reynolds number assuming all the mass flowing as liquid
            D_h = obj.Model.D_h;
            Re_l = obj.Model.G * D_h / props.mu_l;
            fprintf('Re_l: %.2f\n', Re_l);
            fprintf('mu_l: %.5e Pa.s\n', props.mu_l);

            % Heat transfer coefficient for liquid phase (h_l)
            if Re_l < 2300
                W_w = obj.Model.W_w;
                W_c = obj.Model.W_c;

                a_channel = min(W_w, W_c);
                b_channel = max(W_w, W_c);
                alpha = a_channel / b_channel;

                % Nusselt number (Shah and London, 1978) for laminar flow in rectangular channels
                Nu_l = 7.541 * (1 - 2.610*alpha + 4.970*alpha^2 - 5.119*alpha^3 + 2.702*alpha^4 - 0.548*alpha^5);
            else
                Nu_l = 0.023 * (Re_l^0.8) * (props.Pr_l^0.4); % Turbulent flow (Dittus-Boelter)
            end

            fprintf('props.Pr_l: %.4f\n', props.Pr_l);
            fprintf('Nu_l: %.2f\n', Nu_l);
            fprintf('k_l: %.6f W/m.K\n', props.k_l);
            h_l = Nu_l * (props.k_l / D_h);

            P_crit = 4.0593e6; % Critical pressure of R134a [Pa]
            P_r = P / P_crit;

            % Refrigerant heat transfer coefficient h_tpm [W/m^2.K]
            h_tpm = h_l * (0.55 + 2.09 / P_r^0.38); % Empirical correlation for two-phase flow

            UA_tp = obj.Model.UA(w_tp, h_tpm);
            fprintf('L_tp: %.4f\n', L_tp);
            fprintf('w_tp: %.4f\n', w_tp);
            fprintf('UA_tp: %.4f\n', UA_tp);
            NTU_tp = UA_tp / C_air;
            fprintf('NTU_tp: %.4f\n', NTU_tp);
            epsilon = 1 - exp(-NTU_tp);
            fprintf('epsilon: %.4f\n', epsilon);
            % Calculate the actual heat transfer
            Q_eNTU = epsilon * C_air * (T_sat - T_air_in);
        end

    end
end