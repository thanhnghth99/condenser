clear;
clc;
% R134a refrigerant side
Refrigerant = 'R134a';
P_ref_in = 1569064; % Pa
T_ref_in = 355.1500; % K
m_ref = 0.02746; % kg/s
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
A_channel = 5.4e-7; % m^2
W_cond = 0.34; % m

inlet_data = CondenserInlet(P_ref_in, T_ref_in, m_ref, A_tube_total, T_air_in, m_air_in, h_air, eta_o, A_surface_total);
disp(inlet_data);

cond_specs = CondenserModel(inlet_data, D_h, A_channel, N_c, N_t, W_cond);
disp(cond_specs);

% SUPERHEATED (SH) REGION SIMULATION
% Initialize the Superheated Region solver
sh_solver = SuperheatedRegion(cond_specs);
disp(sh_solver);

[L_sh, P_out_sh, dP_sh, T_sat_v, h_sat_v, Q_eNTU_sh] = sh_solver.defineRegion();
% Print the results to the console professionally
fprintf('\n=== SUPERHEATED (SH) REGION SIMULATION RESULTS ===\n');
fprintf('1. Superheated region length (L_sh)     : %.2f m\n', L_sh);
fprintf('2. Outlet pressure (P_out_sh)           : %.2f Pa\n', P_out_sh);
fprintf('3. Pressure drop (dP_sh)                : %.2f Pa\n', dP_sh);
fprintf('4. Saturated temperature (T_sat_v)      : %.2f K\n', T_sat_v);
fprintf('5. Saturated vapor enthalpy (h_sat_v)   : %.2f J/kg\n', h_sat_v);
fprintf('6. Heat Transfer Capacity (Q_sh)        : %.2f W\n', Q_eNTU_sh);
fprintf('==================================================\n\n');


% % TWO-PHASE (TP) REGION SIMULATION
% % Initialize the Two-Phase solver
% tp_solver = TwoPhaseRegion(cond_specs);
% disp(tp_solver);

% % Run the spatial marching solver using inputs from the SH region
% % Note: The input pressure for TP is the outlet pressure of SH (P_out_sh)
% [w_tp, P_out_tp, T_sat_l, h_l, Q_eNTU_tp] = tp_solver.defineRegion(w_sh, P_out_sh);

% % Calculate pressure drop in Two-Phase region
% dP_tp = P_out_sh - P_out_tp;

% % Print the results to the console professionally
% fprintf('\n=== TWO-PHASE (TP) REGION SIMULATION RESULTS ===\n');
% fprintf('1. Area fraction (w_tp)            : %.2f %%\n', w_tp * 100);
% fprintf('2. Outlet pressure (P_out_tp)      : %.2f Pa\n', P_out_tp);
% fprintf('3. Pressure drop (dP_tp)           : %.2f Pa\n', dP_tp);
% fprintf('4. Heat Transfer Capacity (Q_tp)   : %.2f W\n', Q_eNTU_tp);
% fprintf('5. Saturated liquid temp (T_sat_l) : %.2f K\n', T_sat_l);
% fprintf('6. Saturated liquid enthalpy (h_l) : %.2f J/kg\n', h_l);
% fprintf('==================================================\n\n');

