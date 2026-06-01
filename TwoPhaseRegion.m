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
        function [w_tp, P_out_tp, T_sat_f, h_f] = defineRegion(obj, w_sh, P_in_tp, h_in_tp)

            % w is the ratio of the area of ​​each region to the total area of ​​the pipe
            % Lower bound of area fraction
            w_min = 0;
            % Upper bound of area fraction
            w_max = 1 - w_sh;

            % Absolute tolerance for heat transfer error [w]
            tol = 0.0001;

            % Maximun iterations to prevent infinite loops (Failed to convergence)
            iter_max = 100;

            % Loop counter
            iter = 0;

            Refrig = obj.Model.Inlet.Refrigerant;

            % START BISECTION LOOP
            while iter < iter_max
                iter = iter + 1;
                fprintf('Loop %d\n', iter);

                % Step 1: Guess two-phase region area fraction
                w_tp = (w_min + w_max) / 2;

                % Step 2: Calculate frictional pressure drop
                % tp_props = ThermoProp.get_SinglePhaseProps(P_out_sh, T_ref_in, Refrig);

                % Calculate pressure drop and update outlet pressure
                dP_tp = obj.PressureDropTP(w_tp, tp_props.rho, tp_props.mu);
                P_out_tp = P_out_sh - dP_tp;

                if P_out_tp < 101325
                    P_out_tp = 101325;
                    warning("The calculated pressure is negative. It's been forced down to 1 atm to prevent errors.");
                end
                
                % Step 3: Update required energy (Q_req)
                T_sat_f = ThermoProp.get_T_sat(P_out_tp, 0, Refrig);
                h_f = ThermoProp.get_h_sat(P_out_tp, 0, Refrig);

                % Required heat transfer to cool the refrigerant down to the dew point
                Q_req = obj.Model.Inlet.m_ref * (obj.Model.Inlet.h_ref_in - h_f);

                % Step 4: Calculate actual capacity (Q_eNTU)
            end
        end
    end

    methods (Access = private)
        % Function to calculate pressure drop in the two-phase region
        function dP_tp = PressureDropTP(obj, w_tp, P, rho_tp_l, rho_tp_v, mu_tp_l, x)
            L_tp = w_tp * obj.Model.H_cond; % Length of the two-phase region
            if L_tp == 0
                dP_tp = 0; % No two-phase region, so no pressure drop
                return;
            end

            % Reynolds number for two-phase flow (using mass flux and hydraulic diameter)
            G = obj.Model.G; % Mass flux [kg/m^2.s]
            D_h = obj.Model.D_h; % Hydraulic diameter [m]
            Re_tp = G * D_h / mu_tp_l;

            % Friction factor for two-phase flow
            f_lo = 0.079 * Re_tp^(-0.25); % Laminar flow friction factor

            % Pressure drop due to friction
            P_crit = 4.0593e6; % Critical pressure of R134a [Pa]
            P_r = P / P_crit; % Reduced pressure

            term_1 = (1 - x)^2;
            term_2 = 2.2 * x^2 * P_r^(-0.94);
            term_3 = 2.6 * x^0.8 * (1 - x)^0.25 * P_r^(-1.44);

            dP_fr = L_tp * 2 * (f_lo * G^2) / (rho_tp_l * D_h) * (term_1 + term_2 + term_3);

            % Void fraction (alpha) using homogeneous flow model
            alpha = 1 / (1 + ((1 - x) / x) * (rho_tp_v / rho_tp_l)^(2/3));

            % Refrigerant flows in vertically downward direction inside the tubes
            omega_deg = -90;
            g = 9.81; % Gravitational acceleration [m/s^2]

            % Pressure drop due to gravity
            dP_gravity = ((1 - alpha) * rho_tp_l + alpha * rho_tp_v) * g * L_tp * sind(omega_deg);

            % Pressure recovery due to deceleration of the flow
            dP_recovery = 0.5 * G^2 * (alpha / rho_tp_v + (1 - alpha) / rho_tp_l);


        end

        % Function to calculate the actual heat transfer capacity (Q_eNTU) based on the NTU method
        function Q_eNTU = CalculateQeNTU(obj, w_tp, T_ref_in, T_sat_f)
            % Placeholder for calculating the actual heat transfer capacity (Q_eNTU)
            % This should be based on the NTU method and the properties of the refrigerant
            Q_eNTU = 0; % Replace with actual calculation
        end

    end
end