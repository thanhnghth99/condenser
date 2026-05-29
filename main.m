clear;
clc;

inlet_data = CondenserInlet(980.6e3, 337.15, 0.8503e-2, 0.6977, 299.15, 0.204, 116.08, 0.94448, 2.642);
disp(inlet_data);

cond_specs = CondenserModel(inlet_data, 0.00072, 5.4e-7, 19, 36, 0.34);
disp(cond_specs);

sh_solver = SuperheatedRegion(cond_specs);
disp(sh_solver);

[w_sh, P_out_sh, T_sat_v, h_v] = sh_solver.defineRegion();
% Print the results to the console professionally
fprintf('\n=== SUPERHEATED (SH) REGION SIMULATION RESULTS ===\n');
fprintf('1. Area fraction (w_sh)            : %.2f %%\n', w_sh * 100);
fprintf('2. Outlet pressure (P_out_sh)      : %.2f Pa\n', P_out_sh);
fprintf('3. Saturated temperature (T_sat_v) : %.2f K\n', T_sat_v);
fprintf('4. Saturated vapor enthalpy (h_v)  : %.2f J/kg\n', h_v);
fprintf('==================================================\n\n');


