classdef SuperheatedRegion
    properties
        % Reference to the condenser model to access inlet conditions and geometry
        Model
    end

    methods
        function obj = SuperheatedRegion(condenser_model)
            obj.Model = condenser_model;
        end
        
        function [w_sh, P_out_sh, dP_sh, T_sat_v, h_v] = defineRegion(obj)

            % w is the ratio of the area of ​​each region to the total area of ​​the pipe
            % Lower bound of area fraction
            w_min = 0;
            % Upper bound of area fraction
            w_max = 1;

            % Absolute tolerance for heat transfer error [w]
            tol = 0.0001;

            % Maximun iterations to prevent infinite loops (Failed to convergence)
            iter_max = 100;

            % Loop counter
            iter = 0;

            P_ref_in = obj.Model.Inlet.P_ref_in;
            T_ref_in = obj.Model.Inlet.T_ref_in;
            Refrig = obj.Model.Inlet.Refrigerant;

            % START BISECTION LOOP
            while iter < iter_max
                iter = iter + 1;
                fprintf('Loop %d\n', iter);

                % Step 1: Guess superheated region area fraction
                w_sh = (w_min + w_max) / 2;

                % Step 2: Calculate frictional pressure drop
                sh_props = ThermoProp.get_SinglePhaseProps(P_ref_in, T_ref_in, Refrig);

                % Calculate pressure drop and update outlet pressure
                dP_sh = obj.PressureDropSH(w_sh, sh_props.rho, sh_props.mu);
                P_out_sh = obj.Model.Inlet.P_ref_in - dP_sh;

                if P_out_sh < 101325
                    P_out_sh = 101325;
                    dP_sh = obj.Model.Inlet.P_ref_in - P_out_sh;
                    warning("The calculated pressure is negative. It's been forced down to 1 atm to prevent errors.");
                end
                
                % Step 3: Update required energy (Q_req)
                T_sat_v = ThermoProp.get_T_sat(P_out_sh, 1, Refrig);
                h_v = ThermoProp.get_SatVaporProps(P_out_sh, Refrig).h_v;

                % Required heat transfer to cool the refrigerant down to the dew point
                Q_req = obj.Model.Inlet.m_ref * (obj.Model.Inlet.h_ref_in - h_v);

                % Step 4: Calculate actual capacity (Q_eNTU)
                T_rm_sh = (obj.Model.Inlet.T_ref_in + T_sat_v) / 2;
                rm_props = ThermoProp.get_SinglePhaseProps(P_out_sh, T_rm_sh, Refrig);
                
                % Actual heat transfer capacity using e-NTU method
                Q_eNTU = obj.HeatTransfer_eNTU(w_sh, rm_props);

                % Step 5: Convergence check & Boundary update
                Error = Q_eNTU - Q_req;

                if abs(Error / Q_req) < tol
                    fprintf('Convergence after %d loop, w_sh = %.2f%% with error = %.3f W\n', iter, w_sh * 100, Error);
                    break;
                end
                
                % Bisection boundary update logic
                if Error < 0
                    % Actual < Target -> Under-sized -> Increase lower bound
                    w_min = w_sh;
                else
                    % Actual > Target -> Over-sized -> Decrease uppper bound
                    w_max = w_sh;
                end
            end % End of while loop
            
            if iter == iter_max
                % Warning if the loop reachs the iteration limit without converging
                warning('Superheated region FAILED to converge!');                
            end
        end
    end

    methods (Access = private)

        % Pressure drop
        function dP_sh = PressureDropSH(obj, w_sh, rho, mu)
            L_sh = w_sh * obj.Model.H_cond;
            if L_sh <= 1e-6
                dP_sh = 0; % No superheated region, so no pressure drop
                return;
            end

            % Reynolds number
            Re_Dh = obj.Model.G * obj.Model.D_h / mu;

            % Friction factor f
            f = 1 / (1.58 * log(Re_Dh) - 3.28)^2; % Turbulent flow


            % Tube length of superheated region
            L_sh = w_sh * obj.Model.H_cond;

            % Pressure drop dP_sh
            dP_sh = f * L_sh * (obj.Model.G^2) / (2 * obj.Model.D_h * rho);
        end

        % Heat transfer using e-NTU method
        function Q_eNTU = HeatTransfer_eNTU(obj, w_sh, props)
            L_sh = w_sh * obj.Model.H_cond;
            if L_sh <= 1e-6
                Q_eNTU = 0;
                return;
            end    

            % Air side
            % Specific heat capacity of air [J/kg.K]
            cp_air = 1005;
            C_air = w_sh * obj.Model.Inlet.m_air_in * cp_air;

            % Refrigerant side
            % Specific heat capacity of refrigerant
            C_r134a = obj.Model.Inlet.m_ref * props.cp;

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

            Q_eNTU = epsilon * C_min * (obj.Model.Inlet.T_ref_in - obj.Model.Inlet.T_air_in);
        end
    end
end