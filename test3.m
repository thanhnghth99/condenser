% Đoạn mã xuất mảng thông số Nhiệt dung riêng (c_p) của hơi quá nhiệt R134a
% Sử dụng CoolProp qua môi trường Python trong MATLAB

clear; clc;

% 1. Khai báo các thông số đầu vào
% Mảng nhiệt độ (Breakpoints) từ 319.15 K đến 333.15 K, bước nhảy 1 K
T_array = 319.15 : 1 : 333.15; 

% Áp suất ngưng tụ không đổi (Đơn vị: Pa). Giả định 1 MPa (1,000,000 Pa).
P_condenser = 1e6; 

% Tên môi chất
Refrigerant = 'R134a';

% Khởi tạo mảng rỗng để chứa kết quả c_p
cp_array = zeros(1, length(T_array));

% 2. Vòng lặp tính toán dùng CoolProp
fprintf('Đang lấy dữ liệu c_p từ CoolProp...\n');
for i = 1:length(T_array)
    % Gọi hàm PropsSI từ Python
    % 'C' là mã của Specific Heat Capacity (Nhiệt dung riêng, đơn vị: J/kg.K)
    try
        cp_array(i) = py.CoolProp.CoolProp.PropsSI('C', 'T', T_array(i), 'P', P_condenser, Refrigerant);
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

% In Table Data (c_p)
% c_p thường có giá trị lớn (quanh mức 1000 - 1200 J/kg.K), nên lấy 4 chữ số thập phân là quá đủ độ chính xác
table_data_str = sprintf('%.4f ', cp_array);
fprintf('\nTable data (Đơn vị: J/kg.K):\n');
fprintf('[%s]\n', strtrim(table_data_str));

