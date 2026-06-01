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
        function [w_tp, P_out_tp, T_sat_l, h_l, Q_eNTU_tp] = defineRegion(obj, w_sh, P_in_tp)

            Refrig = obj.Model.Inlet.Refrigerant;

            % 1. Geometry and Initialization
            dz_base = 1.5e-3; % 1.5 mm element length for numerical integration
            H_cond = obj.Model.H_cond; % Total length of the condenser [m]

            % dA_ref = obj.Model.Inlet.A_tube_total / H_cond * dz_base; % Refrigerant side area for 1 element [m^2]
            % dA_air = obj.Model.Inlet.A_surface_total / H_cond * dz_base; % Air side area for 1 element [m^2]

            m_ref = obj.Model.Inlet.m_ref; % Mass flow rate of refrigerant [kg/s]

            x_i = 1; % Initial vapor quality at the inlet of the two-phase region (saturated vapor)
            P_i = P_in_tp; % Initial pressure at the inlet of the two-phase region
            Q_eNTU_tp = 0; % Initialize actual heat transfer for the two-phase region
            L_acc = 0; % Accumulated length of the two-phase region
            step = 0; % Step counter for the iteration

            % START SPATIAL MARCHING LOOP
            while x_i > 0
                step = step + 1;
                fprintf('Loop %d\n', step);

                % --- BƯỚC 1: KIỂM TRA GIỚI HẠN CHIỀU DÀI VẬT LÝ ---
                L_avail = H_cond * (1 - w_sh) - L_acc;
                if L_avail <= 1e-6
                    warning('Condenser length is physically not enough to fully condense the refrigerant. Region truncated.');
                    break;
                end

                % Ép bước nhảy không được vượt quá chiều dài còn lại của ống
                dz_step = min(dz_base, L_avail);

                % Tính toán diện tích cho riêng bước nhảy này
                dA_ref_step = obj.Model.Inlet.A_tube_total / H_cond * dz_step; 
                dA_air_step = obj.Model.Inlet.A_surface_total / H_cond * dz_step; 

                if P_i < 101325
                    error("The calculated pressure is negative. Flow crashed!");
                end
                
                T_sat_i = ThermoProp.get_T_sat(P_i, x_i, Refrig);
                rho_l_i = ThermoProp.get_SatLiquidProps(P_i, Refrig).rho_l;
                rho_v_i = ThermoProp.get_SatVaporProps(P_i, Refrig).rho_v;
                mu_l_i = ThermoProp.get_SatLiquidProps(P_i, Refrig).mu_l;
                Pr_l_i = ThermoProp.get_SatLiquidProps(P_i, Refrig).Pr_l;
                k_l_i = ThermoProp.get_SatLiquidProps(P_i, Refrig).k_l;
                h_l_i = ThermoProp.get_SatLiquidProps(P_i, Refrig).h_l;
                h_v_i = ThermoProp.get_SatVaporProps(P_i, Refrig).h_v;

                % Truyền biến dz_step động vào hàm tính nhiệt
                d_QeNTU_i = obj.HeatTransfer_eNTU(x_i, P_i, T_sat_i, mu_l_i, Pr_l_i, k_l_i, dz_step, dA_ref_step, dA_air_step);

                dx_i = d_QeNTU_i / (m_ref * abs(h_v_i - h_l_i));
                
                % --- BƯỚC 2: OVERSHOOT PROTECTION (NỘI SUY KHI x < 0) ---
                if (x_i - dx_i) < 0
                    % Tính tỷ lệ chiều dài chính xác để x vừa chạm đúng 0
                    fraction = x_i / dx_i;
                    dz_exact = fraction * dz_step;
                    
                    % Cập nhật lại diện tích cho đoạn cắt nhỏ này
                    dA_ref_exact = obj.Model.Inlet.A_tube_total / H_cond * dz_exact;
                    dA_air_exact = obj.Model.Inlet.A_surface_total / H_cond * dz_exact;
                    
                    % Tính lại nhiệt lượng và sụt áp cho tỷ lệ fraction chính xác
                    d_QeNTU_i = obj.HeatTransfer_eNTU(x_i, P_i, T_sat_i, mu_l_i, Pr_l_i, k_l_i, dz_exact, dA_ref_exact, dA_air_exact);
                    x_next = 0; % Chốt hạ x = 0 một cách mượt mà
                    
                    [dP_total_i, dP_fr_i, dP_grav_i, dP_dec_i] = obj.PressureDropTP(P_i, rho_l_i, rho_v_i, mu_l_i, x_i, x_next, dz_exact, Refrig);
                    P_next = P_i - dP_total_i;
                    L_acc = L_acc + dz_exact;
                else
                    % Nếu chưa chạm biên, chạy bình thường
                    x_next = x_i - dx_i;
                    [dP_total_i, dP_fr_i, dP_grav_i, dP_dec_i] = obj.PressureDropTP(P_i, rho_l_i, rho_v_i, mu_l_i, x_i, x_next, dz_step, Refrig);
                    P_next = P_i - dP_total_i;
                    L_acc = L_acc + dz_step;
                end
                
                % Update variables
                x_i = x_next;
                P_i = P_next;
                Q_eNTU_tp = Q_eNTU_tp + d_QeNTU_i; 
            end
            
            w_tp = L_acc / H_cond; 
            P_out_tp = P_i; 
            T_sat_l = ThermoProp.get_T_sat(P_out_tp, 0, Refrig); 
            h_l = ThermoProp.get_SatLiquidProps(P_out_tp, Refrig).h_l; 
        end
        %     % START SPATIAL MARCHING LOOP
        %     while x_i > 0
        %         step = step + 1;
        %         fprintf('Loop %d\n', step);

        %         if L_acc >= H_cond * (1 - w_sh)
        %             warning('Condenser length is not enough to fully condense the refrigerant. The two-phase region is truncated at the end of the condenser.');
        %             break;
        %         end

        %         if P_i < 101325
        %             error("The calculated pressure is negative. Flow crashed! Please check the input conditions and model parameters.");
        %         end
                
        %         % Thermodynamic properties of the refrigerant at the current state
        %         T_sat_i = ThermoProp.get_T_sat(P_i, x_i, Refrig);

        %         rho_l_i = ThermoProp.get_SatLiquidProps(P_i, Refrig).rho_l;
        %         rho_v_i = ThermoProp.get_SatVaporProps(P_i, Refrig).rho_v;

        %         mu_l_i = ThermoProp.get_SatLiquidProps(P_i, Refrig).mu_l;

        %         Pr_l_i = ThermoProp.get_SatLiquidProps(P_i, Refrig).Pr_l;

        %         k_l_i = ThermoProp.get_SatLiquidProps(P_i, Refrig).k_l;

        %         h_l_i = ThermoProp.get_SatLiquidProps(P_i, Refrig).h_l;
        %         h_v_i = ThermoProp.get_SatVaporProps(P_i, Refrig).h_v;

        %         % Calculate the actual heat transfer for 1 element using e-NTU method
        %         d_QeNTU_i = obj.HeatTransfer_eNTU(x_i, P_i, T_sat_i, mu_l_i, Pr_l_i, k_l_i, dz, dA_ref, dA_air);

        %         % Update x vapor quality based on the actual heat transfer
        %         dx_i = d_QeNTU_i / (m_ref * abs(h_v_i - h_l_i));
        %         x_next = x_i - dx_i;

        %         % Calculate pressure drop and update outlet pressure
        %         [dP_total_i, dP_fr_i, dP_grav_i, dP_dec_i] = obj.PressureDropTP(P_i, rho_l_i, rho_v_i, mu_l_i, x_i, x_next, dz, Refrig);
        %         P_next = P_i - dP_total_i;
                
        %         % Update variables for the next iteration
        %         x_i = x_next;
        %         P_i = P_next;
        %         L_acc = L_acc + dz;
        %         Q_eNTU_tp = Q_eNTU_tp + d_QeNTU_i; % Store the actual heat transfer for the current element (can be used for post-processing or debugging)
        %     end
        %     % Output variables after convergence or reaching the end of the condenser
        %     w_tp = L_acc / H_cond; % Area fraction of the two-phase region
        %     P_out_tp = P_i; % Outlet pressure at the end of the two-phase region
        %     T_sat_l = ThermoProp.get_T_sat(P_out_tp, 0, Refrig); % Saturated temperature at the outlet of the two-phase region (saturated liquid)
        %     h_l = ThermoProp.get_SatLiquidProps(P_out_tp, Refrig).h_l; % Saturated liquid enthalpy at the outlet of the two-phase region
        % end
    end

    methods (Access = private)
        % Function to calculate pressure drop in the two-phase region for 1 element
        function [dP_total_ele, dP_fr_ele, dP_grav_ele, dP_dec_ele] = PressureDropTP(obj, P_1, rho_l_1, rho_v_1, mu_l_1, x_1, x_2, dz, Refrig)

            % Reynolds number for two-phase flow (using mass flux and hydraulic diameter)
            G = obj.Model.G; % Mass flux [kg/m^2.s]
            D_h = obj.Model.D_h; % Hydraulic diameter [m]
            Re_lo = G * D_h / mu_l_1;

            if Re_lo < 2300
                f_lo = 16 / Re_lo; % Laminar flow
            else
                f_lo = 0.079 * Re_lo^(-0.25); % Turbulent flow (Blasius)
            end

            % 1. Pressure drop due to friction
            P_crit = 4.0593e6; % Critical pressure of R134a [Pa]
            P_r = P_1 / P_crit; % Reduced pressure

            if x_1 >= 1
                term_1 = 0;
                term_2 = 2.2 * P_r^(-0.94);
                term_3 = 0;
            elseif x_1 <= 0
                term_1 = 1;
                term_2 = 0;
                term_3 = 0;
            else
                term_1 = (1 - x_1)^2;
                term_2 = 2.2 * x_1^2 * P_r^(-0.94);
                term_3 = 2.6 * x_1^0.8 * (1 - x_1)^0.25 * P_r^(-1.44);
            end

            dp_dz_fr = 2 * (f_lo * G^2) / (rho_l_1 * D_h) * (term_1 + term_2 + term_3);
            dP_fr_ele = dp_dz_fr * dz;

            % Void fraction and Momentum at 1st point
            if x_1 >= 1
                alpha_1 = 1;
                M_1 = G^2 / rho_v_1;
            elseif x_1 <= 0
                alpha_1 = 0;
                M_1 = G^2 / rho_l_1;
            else
                alpha_1 = 1 / (1 + ((1 - x_1) / x_1) * (rho_v_1 / rho_l_1)^(2/3));
                M_1 = G^2 * ((x_1^2)/(rho_v_1 * alpha_1) + ((1 - x_1)^2)/(rho_l_1 * (1 - alpha_1)));
            end

            % Refrigerant flows in vertically downward direction inside the tubes
            omega_deg = -90;
            g = 9.81; % Gravitational acceleration [m/s^2]

            % 2. Pressure drop due to gravity
            dp_dz_grav = ((1 - alpha_1) * rho_l_1 + alpha_1 * rho_v_1) * g * sind(omega_deg);
            dP_grav_ele = dp_dz_grav * dz;

            % 3. Pressure drop due to deceleration of the flow (Predictor-Corrector Method)
            % Predictor step: Assume densities at x_2 are identical to  x_1 to estimate tentative pressure
            if x_2 >= 1
                alpha_2_tent = 1;
                M_2_tent = G^2 / rho_v_1;
            elseif x_2 <= 0
                alpha_2_tent = 0;
                M_2_tent = G^2 / rho_l_1;
            else
                alpha_2_tent = 1 / (1 + ((1 - x_2) / x_2) * (rho_v_1 / rho_l_1)^(2/3));
                M_2_tent = G^2 * ((x_2^2)/(rho_v_1 * alpha_2_tent) + ((1 - x_2)^2)/(rho_l_1 * (1 - alpha_2_tent)));
            end

            % Estimate tentative pressure drop and resulting pressure
            dP_dec_tent = M_2_tent - M_1;
            dP_total_tent = dP_fr_ele + dP_grav_ele + dP_dec_tent;
            P_2_tent = P_1 - dP_total_tent;
            P_2_tent = max(P_2_tent, 101325); % Prevent internal CoolProp errors if tentative pressure drops below 1 atm

            % Corrector step: Retrieve accurate densities based on the tentative pressure (P_2_tent)
            rho_l_2 = ThermoProp.get_SatLiquidProps(P_2_tent, Refrig).rho_l;
            rho_v_2 = ThermoProp.get_SatVaporProps(P_2_tent, Refrig).rho_v;
            
            % Calculate the final, corrected Momentum at 2nd point
            if x_2 >= 1
                alpha_2 = 1;
                M_2 = G^2 / rho_v_2;
            elseif x_2 <= 0
                alpha_2 = 0;
                M_2 = G^2 / rho_l_2;
            else
                alpha_2 = 1 / (1 + ((1 - x_2) / x_2) * (rho_v_2 / rho_l_2)^(2/3)); % Update void fraction at 2nd point
                M_2 = G^2 * ((x_2^2)/(rho_v_2 * alpha_2) + ((1 - x_2)^2)/(rho_l_2 * (1 - alpha_2)));
            end

            % Final decelation pressure
            dP_dec_ele = M_2 - M_1;

            % Total pressure drop in the two-phase region
            dP_total_ele = dP_fr_ele + dP_grav_ele + dP_dec_ele;
        end

        % Function to calculate the actual heat transfer capacity (Q_eNTU) based on the NTU method for 1 element
        function d_QeNTU_ele = HeatTransfer_eNTU(obj, x, P, T_sat, mu_l, Pr_l, k_l, dz, dA_ref, dA_air)
            % Air side
            % Specific heat capacity of air [J/kg.K] for 1 element
            cp_air = 1005;
            m_air_ele = obj.Model.Inlet.m_air_in * (dz / obj.Model.H_cond); % Mass flow rate of air for 1 element
            C_air_ele = m_air_ele * cp_air;

            % Reynolds number assuming all the mass flowing as liquid
            D_h = obj.Model.D_h;
            Re_l = obj.Model.G * D_h / mu_l;

            if Re_l < 2300
                Nu_l = 3.66 * 1.5; % Laminar flow
            else
                Nu_l = 0.023 * (Re_l^0.8) * (Pr_l^0.4); % Turbulent flow (Dittus-Boelter)
            end
            % Heat transfer coefficient for liquid phase (h_l)
            h_l = Nu_l * (k_l / D_h);

            P_crit = 4.0593e6; % Critical pressure of R134a [Pa]
            P_r = P / P_crit;

            % At x = 1; Condensation just starts, h_ref cannot be zero. Not real in physics.
            if x >= 1
                x_calc = 1 - 1e-6; % Avoid division by zero or extremely small numbers
            elseif x <= 0
                x_calc = 1e-6; % Avoid division by zero or extremely small numbers
            else
                x_calc = x;
            end

            % Refrigerant heat transfer coefficient h_ref_ele [W/m^2.K]
            h_ref_ele = h_l * ((1 - x_calc)^0.8 + (3.8 * x_calc^0.76 * (1 - x_calc)^0.04) / (P_r^0.38));
            % Thermal resistance for refrigerant side for 1 element
            if h_ref_ele <= 1e-6
                R_ref_ele = Inf; % Avoid division by zero, set thermal resistance to infinity if heat transfer coefficient is very small
            else
                R_ref_ele = 1 / (h_ref_ele * dA_ref);
            end

            % Air side thermal resistance for 1 element
            h_air = obj.Model.Inlet.h_air;
            eta_o = obj.Model.Inlet.eta_o;
            R_air_ele = 1 / (h_air * dA_air * eta_o);

            % e-NTU method to calculate actual heat transfer for 1 element
            UA_ele = 1 / (R_ref_ele + R_air_ele);   

            NTU_ele = UA_ele / C_air_ele;
            epsilon = 1 - exp(-NTU_ele);

            % Calculate the actual heat transfer for 1 element
            d_QeNTU_ele = epsilon * C_air_ele * (T_sat - obj.Model.Inlet.T_air_in);
        end

    end
end