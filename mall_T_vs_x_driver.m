function mall_T_vs_x_driver()
% 画出商场场景中 T_total(x) 随消防员人数 x 变化的曲线
% 要求已经存在函数:  mall_sim_once(x, p_room, verbose)

clc; close all;
fprintf('=== Plot T_total(x) for mall scenario ===\n');

%% 固定每个房间的并行人数策略（用你现在认为合理/最优的那组）
% 例如：前两个房间 2 人并行，商店 1 人，儿童区 3 人：
p_room = [2; 4; 2; 2; 2];    % 这里可以改成你实际用的那组

%% 扫描的人数范围
x_vals = 5:10;               % 从 2 人到 10 人
T_vals = zeros(size(x_vals));

for k = 1:numel(x_vals)
    x = x_vals(k);
    % 第三个参数 verbose=false，这里只要 T_total 不打印细节
    T_vals(k) = mall_sim_once(x, p_room, false);
    fprintf('x = %d, T_total = %.1f s (%.2f min)\n', ...
            x, T_vals(k), T_vals(k)/60);
end

%% 画图
figure('Name','Mall T_{total}(x)');
plot(x_vals, T_vals/60, '-o', 'LineWidth', 1.5);
grid on;
xlabel('Number of firefighters x');
ylabel('Total time T_{total} (min)');
title('Mall scenario: T_{total} vs number of firefighters');

%% 计算边际收益（多一个人节省多少时间）
dT = diff(T_vals)/60;   % 单位：分钟
fprintf('\n--- Marginal time reduction (ΔT when x -> x+1) ---\n');
for k = 1:numel(dT)
    fprintf('from x = %d to %d: ΔT = %.2f min\n', ...
            x_vals(k), x_vals(k+1), dT(k));
end

end
