clear;
clc;
% R134a refrigerant side
Refrigerant = 'R134a';
P_ref_in = 1569064; % Pa
T_ref_in = 355.1500; % K
m_ref = 0.0227; % kg/s --> 900 - 1000 rpm
A_tube_total = 0.6977; % m^2

% Air side
T_air_in = 303.15; % K
m_air_in = 0.2014; % kg/s
h_air = 115.79; % W/m^2.K
eta_o = 0.94448; % Overall surface efficiency 
A_surface_total  = 2.642; % m^2

% Condenser model specifications
D_h = 0.00072; % m
N_c = 19; % Number of channels
N_t = 36; % Number of tubes
W_w = 0.0006; % Channel height [m]
W_c = 0.0009; % Channel width [m]
A_channel = 5.4e-7; % m^2
W_cond = 0.34; % m
sigma = 0.7315;

% 2. TÍNH TOÁN CÁC THÔNG SỐ TRUNG GIAN (Thay thế cho Class CondenserModel)
% =========================================================================
% Mass flux of refrigerant [kg/m^2.s]
G = m_ref / (A_channel * N_c * N_t);

% Aspect ratio
a_channel = min(2 * W_w, 2 * W_c);
b_channel = max(2 * W_w, 2 * W_c);
alpha = a_channel / b_channel;

% Initial enthalpy for cross-check (Gọi trực tiếp thư viện)
h_ref_in = py.CoolProp.CoolProp.PropsSI('H', 'P', P_ref_in, 'T', T_ref_in, Refrigerant);

% =========================================================================
% 3. CHẠY CÁC KHỐI MÔ PHỎNG (Sẵn sàng cho Simulink)
% =========================================================================

% --- SUPERHEATED (SH) REGION SIMULATION ---
disp('Running Superheated Block...');
[L_sh, P_out_sh, dP_sh_total, dP_entrance, T_sat_v, h_sat_v, Q_eNTU_sh] = Block_Superheated(...
    P_ref_in, T_ref_in, m_ref, T_air_in, m_air_in, h_air, D_h, alpha, W_cond, A_tube_total, A_surface_total, eta_o, G, sigma);

fprintf('\n=== SUPERHEATED (SH) REGION SIMULATION RESULTS ===\n');
fprintf('1. Superheated region length (L_sh)     : %.4f m\n', L_sh);
fprintf('2. Outlet pressure (P_out_sh)           : %.2f Pa\n', P_out_sh);
fprintf('3. Pressure drop (dP_sh)                : %.2f Pa\n', dP_sh_total);
fprintf('3. Minor losses at entrance (dP_entrance)   : %.2f Pa\n', dP_entrance);
fprintf('4. Saturated temperature (T_sat_v)      : %.2f K\n', T_sat_v);
fprintf('5. Saturated vapor enthalpy (h_sat_v)   : %.2f J/kg\n', h_sat_v);
fprintf('6. Heat Transfer Capacity (Q_sh)        : %.2f W\n', Q_eNTU_sh);
fprintf('==================================================\n\n');


% --- TWO-PHASE (TP) REGION SIMULATION ---
disp('Running Two-Phase Block...');
[L_tp, P_out_tp, T_sat, h_sat_l, Q_eNTU_tp] = Block_TwoPhase(...
    L_sh, P_out_sh, m_ref, T_air_in, m_air_in, h_air, D_h, alpha, W_cond, A_tube_total, A_surface_total, eta_o, G);

fprintf('\n=== TWO-PHASE (TP) REGION SIMULATION RESULTS ===\n');
fprintf('1. Two-phase region length (L_tp)       : %.4f m\n', L_tp);
fprintf('2. Outlet pressure (P_out_tp)           : %.2f Pa\n', P_out_tp);
fprintf('3. Saturated temperature (T_sat)        : %.2f K\n', T_sat);
fprintf('4. Saturated vapor enthalpy (h_sat_v)   : %.2f J/kg\n', h_sat_v);
fprintf('5. Saturated liquid enthalpy (h_sat_l)  : %.2f J/kg\n', h_sat_l);
fprintf('6. Heat Transfer Capacity (Q_tp)        : %.2f W\n', Q_eNTU_tp);
fprintf('==================================================\n\n');


% --- SUBCOOLED (SC) REGION SIMULATION ---
disp('Running Subcooled Block...');
[L_sc, P_out_sc, T_out_sc, Q_eNTU_sc, dP_sc, dT_sc] = Block_SubCooled(...
    L_sh, L_tp, P_out_tp, T_sat, m_ref, T_air_in, m_air_in, h_air, D_h, alpha, W_cond, A_tube_total, A_surface_total, eta_o, G);

fprintf('\n=== SUBCOOLED (SC) REGION SIMULATION RESULTS ===\n');
fprintf('1. Subcooled region length (L_sc)       : %.4f m\n', L_sc);
fprintf('2. Outlet pressure (P_out_sc)           : %.2f Pa\n', P_out_sc);
fprintf('3. Subcooled temperature (T_sc_out)     : %.2f K\n', T_out_sc);
fprintf('4. Pressure drop (dP_sc)                : %.2f Pa\n', dP_sc);
fprintf('5. Temperature drop (dT_sc)             : %.2f K\n', dT_sc);
fprintf('6. Heat Transfer Capacity (Q_sc)        : %.2f W\n', Q_eNTU_sc);
fprintf('==================================================\n\n');


% =========================================================================
% 4. OVERALL CONDENSER PERFORMANCE
% =========================================================================

% 1. Total Heat transfer capacity
Q_total_condenser = Q_eNTU_sh + Q_eNTU_tp + Q_eNTU_sc;

% 2. Total Pressure Drop
dP_total = dP_sh_total + (P_out_sh - P_out_tp) + dP_sc; % dP_tp ~ 0 or P_in - P_out

% Print the results to the console professionally
fprintf('\n==================================================\n');
fprintf('    OVERALL CONDENSER PERFORMANCE SUMMARY\n');
fprintf('==================================================\n');
fprintf('1. Total Heat Transfer Capacity (Q_tot): %.2f W\n', Q_total_condenser);
fprintf('2. Total Pressure Drop (dP_tot)        : %.2f Pa\n', dP_total);
fprintf('3. Final Outlet Temperature            : %.2f K (%.2f oC)\n', T_out_sc, T_out_sc - 273.15);
fprintf('4. Final Outlet Pressure               : %.2f Pa\n', P_out_sc);
fprintf('==================================================\n\n');


% =========================================================================
% 5. CROSS-CHECK USING ENTHALPY BALANCE
% =========================================================================
h_in_total = h_ref_in;

% SAFE CALCULATION OF OUTLET ENTHALPY (h_out_final)
if L_sc <= 1e-5
    % Case 1: Incomplete condensation (Outlet state is a two-phase mixture)
    h_out_final = h_sat_v - (Q_eNTU_tp / m_ref);
else
    % Case 2: Successful liquid subcooling (Subcooled liquid state)
    % Use CoolProp directly since we don't have ThermoProp class here anymore
    cp_l = py.CoolProp.CoolProp.PropsSI('C', 'P', P_out_sc, 'T', T_sat - 0.5, Refrigerant);
    h_out_final = h_sat_l - cp_l * (T_sat - T_out_sc);
end

% Calculate total heat transfer based on macroscopic energy balance
Q_enthalpy_balance = m_ref * (h_in_total - h_out_final);

fprintf('\n>>> CROSS-CHECK VALIDATION:\n');
fprintf('- Q calculated by Zone-by-Zone (e-NTU) : %.2f W\n', Q_total_condenser);
fprintf('- Q calculated by Global Enthalpy Drop : %.2f W\n', Q_enthalpy_balance);
fprintf('- Error difference                     : %.4f W\n', abs(Q_total_condenser - Q_enthalpy_balance));
