clear;
clc;
% R134a refrigerant side
Refrigerant = 'R134a';
P_ref_in = 1569064; % Pa
T_ref_in = 355.1500; % K
m_ref = 0.015; % kg/s --> 900 - 1000 rpm
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

inlet_data = CondenserInlet(P_ref_in, T_ref_in, m_ref, A_tube_total, T_air_in, m_air_in, h_air, eta_o, A_surface_total);
disp(inlet_data);

cond_specs = CondenserModel(inlet_data, D_h, A_channel, N_c, N_t, W_w, W_c, W_cond);
disp(cond_specs);

% SUPERHEATED (SH) REGION SIMULATION
% Initialize the Superheated Region solver
sh_solver = SuperheatedRegion(cond_specs);
disp(sh_solver);

[L_sh, P_out_sh, dP_sh, T_sat_v, h_sat_v, Q_eNTU_sh] = sh_solver.defineRegion();
% Print the results to the console professionally
fprintf('\n=== SUPERHEATED (SH) REGION SIMULATION RESULTS ===\n');
fprintf('1. Superheated region length (L_sh)     : %.4f m\n', L_sh);
fprintf('2. Outlet pressure (P_out_sh)           : %.2f Pa\n', P_out_sh);
fprintf('3. Pressure drop (dP_sh)                : %.2f Pa\n', dP_sh);
fprintf('4. Saturated temperature (T_sat_v)      : %.2f K\n', T_sat_v);
fprintf('5. Saturated vapor enthalpy (h_sat_v)   : %.2f J/kg\n', h_sat_v);
fprintf('6. Heat Transfer Capacity (Q_sh)        : %.2f W\n', Q_eNTU_sh);
fprintf('==================================================\n\n');


% TWO-PHASE (TP) REGION SIMULATION
% Initialize the Two-Phase Region solver
tp_solver = TwoPhaseRegion(cond_specs);
disp(tp_solver);

[L_tp, P_out_tp, T_sat, h_sat_v, h_sat_l, Q_eNTU_tp] = tp_solver.defineRegion(L_sh, P_out_sh);

% Print the results to the console professionally
fprintf('\n=== TWO-PHASE (TP) REGION SIMULATION RESULTS ===\n');
fprintf('1. Two-phase region length (L_tp)       : %.4f m\n', L_tp);
fprintf('2. Outlet pressure (P_out_tp)           : %.2f Pa\n', P_out_tp);
fprintf('3. Saturated temperature (T_sat)        : %.2f K\n', T_sat);
fprintf('4. Saturated vapor enthalpy (h_sat_v)   : %.2f J/kg\n', h_sat_v);
fprintf('5. Saturated liquid enthalpy (h_sat_l)  : %.2f J/kg\n', h_sat_l);
fprintf('6. Heat Transfer Capacity (Q_tp)        : %.2f W\n', Q_eNTU_tp);
fprintf('==================================================\n\n');


% SUBCOOLED (SC) REGION SIMULATION
% Initialize the Subcooled Region solver
sc_solver = SubCooledRegion(cond_specs);
disp(sc_solver);

[L_sc, P_out_sc, T_out_sc, Q_eNTU_sc, dP_sc, dT_sc] = sc_solver.defineRegion(L_sh, L_tp, P_out_tp, T_sat);

% Print the results to the console professionally
fprintf('\n=== SUBCOOLED (SC) REGION SIMULATION RESULTS ===\n');
fprintf('1. Subcooled region length (L_sc)       : %.4f m\n', L_sc);
fprintf('2. Outlet pressure (P_out_sc)           : %.2f Pa\n', P_out_sc);
fprintf('3. Subcooled temperature (T_sc_out)     : %.2f K\n', T_out_sc);
fprintf('4. Pressure drop (dP_sc)                : %.2f J/kg\n', dP_sc);
fprintf('5. Temperature drop (dT_sc)             : %.2f J/kg\n', dT_sc);
fprintf('6. Heat Transfer Capacity (Q_sc)        : %.2f W\n', Q_eNTU_sc);
fprintf('==================================================\n\n');


% =========================================================================
% TỔNG KẾT TOÀN BỘ DÀN NGƯNG (OVERALL CONDENSER PERFORMANCE)
% =========================================================================

% 1. Tính tổng công suất truyền nhiệt
Q_total_condenser = Q_eNTU_sh + Q_eNTU_tp + Q_eNTU_sc;

% 2. Tính tổng tổn thất áp suất (Total Pressure Drop)
dP_total = dP_sh + (P_out_sh - P_out_tp) + dP_sc; % Vùng 2 pha coi như dP ~ 0 hoặc lấy P_in - P_out

% In kết quả ra màn hình
fprintf('\n==================================================\n');
fprintf('    OVERALL CONDENSER PERFORMANCE SUMMARY\n');
fprintf('==================================================\n');
fprintf('1. Total Heat Transfer Capacity (Q_tot): %.2f W\n', Q_total_condenser);
fprintf('2. Total Pressure Drop (dP_tot)        : %.2f Pa\n', dP_total);
fprintf('3. Final Outlet Temperature            : %.2f K (%.2f oC)\n', T_out_sc, T_out_sc - 273.15);
fprintf('4. Final Outlet Pressure               : %.2f Pa\n', P_out_sc);
fprintf('==================================================\n\n');


% KIỂM TRA CHÉO BẰNG CÂN BẰNG ENTHALPY
% Lấy enthalpy tại cổng vào dàn ngưng
h_in_total = inlet_data.h_ref_in;

% Tính enthalpy tại cổng ra cuối cùng (trạng thái lỏng quá lạnh)
props_out_final = ThermoProp.get_SinglePhaseProps(P_out_sc, T_out_sc, Refrigerant);
h_out_final = props_out_final.h;

% Tính tổng nhiệt lượng theo cân bằng năng lượng vĩ mô
Q_enthalpy_balance = m_ref * (h_in_total - h_out_final);

fprintf('>> CROSS-CHECK VALIDATION:\n');
fprintf('- Q calculated by Zone-by-Zone (e-NTU) : %.2f W\n', Q_total_condenser);
fprintf('- Q calculated by Global Enthalpy Drop : %.2f W\n', Q_enthalpy_balance);
fprintf('- Error difference                     : %.4f W\n', abs(Q_total_condenser - Q_enthalpy_balance));
