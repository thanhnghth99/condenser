% Đoạn mã xuất mảng thông số Độ nhớt động lực học của hơi quá nhiệt R134a
% Sử dụng CoolProp qua môi trường Python trong MATLAB

clear; clc;

% 1. Khai báo các thông số đầu vào
% Mảng nhiệt độ (Breakpoints) từ 319.15 K đến 333.15 K, bước nhảy 1 K
T_array = 319.15 : 1 : 333.15; 

% Áp suất ngưng tụ không đổi (Đơn vị: Pa). Giả định 1 MPa.
P_condenser = 1e6; 

% Tên môi chất
Refrigerant = 'R134a';

% Khởi tạo mảng rỗng để chứa kết quả độ nhớt
mu_array = zeros(1, length(T_array));

% 2. Vòng lặp tính toán dùng CoolProp
fprintf('Đang lấy dữ liệu từ CoolProp...\n');
for i = 1:length(T_array)
    % Gọi hàm PropsSI từ Python
    % 'V' là mã của Dynamic Viscosity (Độ nhớt động lực học, đơn vị: Pa.s)
    % 'T' là Nhiệt độ, 'P' là Áp suất
    try
        mu_array(i) = py.CoolProp.CoolProp.PropsSI('V', 'T', T_array(i), 'P', P_condenser, Refrigerant);
    catch ME
        warning(['Lỗi tại nhiệt độ ', num2str(T_array(i)), ' K: ', ME.message]);
    end
end

% 3. In kết quả ra màn hình chuẩn định dạng copy vào Simulink 1-D Lookup Table
fprintf('\n--- KẾT QUẢ ĐỂ PASTE VÀO SIMULINK ---\n');

% In Breakpoints (Nhiệt độ)
breakpoints_str = sprintf('%.2f ', T_array);
fprintf('\nBreakpoints 1:\n');
fprintf('[%s]\n', strtrim(breakpoints_str));

% In Table Data (Độ nhớt) với 9 chữ số thập phân
table_data_str = sprintf('%.9f ', mu_array);
fprintf('\nTable data:\n');
fprintf('[%s]\n', strtrim(table_data_str));

a = "[0.000012682 0.000012725 0.000012769 0.000012812 0.000012855 0.000012898 0.000012941 0.000012984 0.000013027 0.00001307 0.000013113 0.000013155 0.000013198 0.00001324 0.000013283]";