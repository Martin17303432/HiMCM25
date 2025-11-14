function task3_T_vs_x_driver()
% 回字型教学楼：画出 T_total(x) 曲线，并观察边际收益

clc; close all;
fprintf('=== Task3 Hui-building: T_total(x) from 2 to 10 ===\n');

x_vals = 2:8;
T_vals = zeros(size(x_vals));
S_vals = zeros(size(x_vals));

for k = 1:numel(x_vals)
    x = x_vals(k);
    [T_vals(k), S_vals(k)] = task3_hui_sim_once(x, false);  % 不打印细节
    fprintf('x = %d, T_total = %.1f s (%.2f min), Stotal = %.1f\n', ...
        x, T_vals(k), T_vals(k)/60, S_vals(k));
end

%% 画时间曲线
figure('Name','Task3 Hui-building: T_{total}(x)');
plot(x_vals, T_vals/60, '-o', 'LineWidth', 1.5);
grid on;
xlabel('Firefighters x');
ylabel('Total time T_{total} (min)');
title('Hui-building scenario: T_{total} vs x');

%% 计算边际收益（多 1 人减少多少分钟）
dT = diff(T_vals)/60;  % 分钟
fprintf('\n--- Marginal time reduction (ΔT when x -> x+1) ---\n');
for k = 1:numel(dT)
    fprintf('from x = %d to %d: ΔT = %.2f min\n', ...
        x_vals(k), x_vals(k+1), dT(k));
end

end
