classdef SubCooledRegion
    properties
        Model
    end

    methods
        function obj = SubCooledRegion(condenser_model)
            obj.Model = condenser_model;
        end

        function [L_sc, P_out_sc, T_out_sc, Q_eNTU_sc, dP_sc, dT_sc] = defineRegion(obj, L_sh, L_tp, P_out_tp, T_sat)
            
            L_total = obj.Model.W_cond;
            L_sc = L_total - L_sh - L_tp;
            
            % Check if no length left for subcooling, exit immediately
            if L_sc <= 1e-5
                Q_eNTU_sc = 0;
                T_out_sc = T_sat;
                P_out_sc = P_out_tp;
                dP_sc = 0;
                dT_sc = 0;
                L_sc = 0;
                return;
            end

            w_sc = L_sc / L_total;

            P_in_sc = P_out_tp; % Initial pressure at the inlet of the two-phase region (outlet pressure of SH region)

            m_ref = obj.Model.Inlet.m_ref;
            Refrig = obj.Model.Inlet.Refrigerant;

            props_sat = ThermoProp.get_SatLiquidProps(P_in_sc, Refrig);
            
            % Retrieve cp from single-phase properties slightly below T_sat
            props_single = ThermoProp.get_SinglePhaseProps(P_in_sc, T_sat - 0.5, Refrig);
            cp_l = props_single.cp;

            % Calculate Reynolds number once to pass into sub-functions
            D_h = obj.Model.D_h;
            G = obj.Model.G;
            mu_l = props_sat.mu_l;
            Re_l = G * D_h / mu_l;

            % Calculate Heat Transfer Capacity using the private e-NTU method
            Q_eNTU_sc = obj.HeatTransfer_eNTU(w_sc, T_sat, props_sat, cp_l, Re_l);

            % Calculate final thermodynamic states
            T_out_sc = T_sat - Q_eNTU_sc / (m_ref * cp_l);
            dT_sc = T_sat - T_out_sc;

            % 4. Calculate Hydraulic Pressure Drop using the private method
            dP_sc = obj.PressureDropSC(L_sc, Re_l, props_sat);
            P_out_sc = P_in_sc - dP_sc;
        end
    end

    methods (Access = private)
        function dP_sc = PressureDropSC(obj, L_sc, Re_l, props_sat)
            D_h = obj.Model.D_h;
            G = obj.Model.G;
            rho_l = props_sat.rho_l;

            if Re_l < 2300
                % Laminar flow in rectangular flat tube (Fanning friction factor)
                alpha = obj.Model.alpha;
                f_lo = (24 / Re_l) * (1 - 1.3553*alpha + 1.9467*alpha^2 - 1.7012*alpha^3 + 0.9564*alpha^4 - 0.2537*alpha^5);
            else
                f_lo = 0.079 * Re_l^(-0.25); % Turbulent friction factor
            end
            
            % Pressure drop (Using Fanning friction factor)
            dP_sc = 2 * L_sc * f_lo * G^2 / (rho_l * D_h);
        end

        function Q_eNTU = HeatTransfer_eNTU(obj, w_sc, T_sat, props_sat, cp_l, Re_l)
            % Air side
            T_air_in = obj.Model.Inlet.T_air_in;
            m_air_in = obj.Model.Inlet.m_air_in;
            cp_air = 1005;
            C_air = w_sc * m_air_in * cp_air;

            % Refrigerant side
            m_ref = obj.Model.Inlet.m_ref;
            C_r134a = m_ref * cp_l;

            C_min = min(C_air, C_r134a);
            C_max = max(C_air, C_r134a);
            C_r = C_min / C_max;

            % Reynolds number
            D_h = obj.Model.D_h;

            % Prandtl number
            Pr_l = props_sat.Pr_l;

            % Thermal conductivity
            k_l = props_sat.k_l;

            % Heat transfer coefficient for liquid phase (h_l)
            if Re_l < 2300
                alpha = obj.Model.alpha;
                % Nusselt number (Shah and London, 1978) for laminar flow in rectangular channels
                Nu_l = 7.541 * (1 - 2.610*alpha + 4.970*alpha^2 - 5.119*alpha^3 + 2.702*alpha^4 - 0.548*alpha^5);
            else
                Nu_l = 0.023 * (Re_l^0.8) * (Pr_l^0.4); % Turbulent flow (Dittus-Boelter)
            end

            h_sc = Nu_l * k_l / D_h;
            UA_sc = obj.Model.UA(w_sc, h_sc);

            % Effectiveness-NTU calculation
            if C_min > 0
                NTU_sc = UA_sc / C_min;
                epsilon = 1 - exp((NTU_sc^0.22 / C_r) * (exp(-C_r * NTU_sc^0.78) - 1));
                Q_eNTU = epsilon * C_min * (T_sat - T_air_in);
            else
                Q_eNTU = 0;
            end
        end
    end
end