function optimize_mall_strategy_8()
% 在 x = 8 时, 穷举搜索每个房间并行人数 p_room,
% 找到总时间 T_total 最小的策略

x = 6;

% 每个房间并行人数的候选集合, 你可以改范围
cand = [2 4 2 2 2 ];    % 表示 1 人并行, 2 人并行, 3 人并行

best_T = inf;
best_p = [];

% 穷举 5 个房间的 p_room 组合
for p1 = cand
  for p2 = cand
    for p3 = cand
      for p4 = cand
        for p5 = cand
            p_room = [p1; p2; p3; p4; p5];

            % 只想看 T_total, 先关掉 verbose = false
            T_total = mall_sim_once(x, p_room, false);

            if T_total < best_T
                best_T = T_total;
                best_p = p_room;
            end
        end
      end
    end
  end
end

fprintf('\n================ Optimization result (x = %d) ================\n', x);
fprintf('Best p_room = [%s]\n', num2str(best_p(:).'));
fprintf('Best T_total = %.1f s (%.2f min)\n', best_T, best_T/60);

% 最后再用最佳方案跑一遍, 打印详细表
fprintf('\n--- Run simulation with best strategy ---\n');
mall_sim_once(x, best_p, true);

end
