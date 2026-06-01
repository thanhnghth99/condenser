clear;
clc;

inlet_data = CondenserInlet(1323900, 350.1500, 0.019551, 0.6977, 299.15, 0.204, 116.08, 0.94448, 2.642);
disp(inlet_data);

cond_specs = CondenserModel(inlet_data, 0.00072, 5.4e-7, 19, 36, 0.34);
disp(cond_specs);

sh_solver = SuperheatedRegion(cond_specs);
disp(sh_solver);

[w_sh, P_out_sh, dP_sh, T_sat_v, h_v] = sh_solver.defineRegion();
% Print the results to the console professionally
fprintf('\n=== SUPERHEATED (SH) REGION SIMULATION RESULTS ===\n');
fprintf('1. Area fraction (w_sh)            : %.2f %%\n', w_sh * 100);
fprintf('2. Outlet pressure (P_out_sh)      : %.2f Pa\n', P_out_sh);
fprintf('3. Pressure drop (dP_sh)           : %.2f Pa\n', dP_sh);
fprintf('4. Saturated temperature (T_sat_v) : %.2f K\n', T_sat_v);
fprintf('4. Saturated vapor enthalpy (h_v)  : %.2f J/kg\n', h_v);
fprintf('==================================================\n\n');


% TWO-PHASE (TP) REGION SIMULATION
% Initialize the Two-Phase solver
tp_solver = TwoPhaseRegion(cond_specs);
disp(tp_solver);

% Run the spatial marching solver using inputs from the SH region
% Note: The input pressure for TP is the outlet pressure of SH (P_out_sh)
[w_tp, P_out_tp, T_sat_l, h_l, Q_eNTU_tp] = tp_solver.defineRegion(w_sh, P_out_sh);

% Calculate pressure drop in Two-Phase region
dP_tp = P_out_sh - P_out_tp;

% Print the results to the console professionally
fprintf('\n=== TWO-PHASE (TP) REGION SIMULATION RESULTS ===\n');
fprintf('1. Area fraction (w_tp)            : %.2f %%\n', w_tp * 100);
fprintf('2. Outlet pressure (P_out_tp)      : %.2f Pa\n', P_out_tp);
fprintf('3. Pressure drop (dP_tp)           : %.2f Pa\n', dP_tp);
fprintf('4. Heat Transfer Capacity (Q_tp)   : %.2f W\n', Q_eNTU_tp);
fprintf('5. Saturated liquid temp (T_sat_l) : %.2f K\n', T_sat_l);
fprintf('6. Saturated liquid enthalpy (h_l) : %.2f J/kg\n', h_l);
fprintf('==================================================\n\n');

