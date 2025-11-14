function task2_full_sim()
% Task 2: 单层办公楼火灾应急排查模型（基于任务一公式）
% 输出：1) 详细计算表 2) 安全度表 3) 人数-总时间曲线 4) 边际参数自动标定

clc; close all;
fprintf('\n=== Task 2 Simulation: Single-floor Office Fire Inspection ===\n');

%% ================= 场景参数（与任务二设定一致） =================
m = 6;                                 % 房间数 R1..R6（普通办公室）
A  = 20*ones(m,1);                     % 面积 (m^2)
alphaTi = ones(m,1);                   % α(office) = 1.0
gammaV  = 0.1*ones(m,1);               % 遮挡率
k0 = floor(0.1*A.*alphaTi + 0.2*gammaV); % 预设扫描点（应为 2）

% 两名专业消防员
x = 2;
Sskill = [4;4];                        % 技能参数（专业=2）
v = [2;2];                         % 走动速度 (m/s)（按任务二示例）

% 面积/技能函数（任务一定义）
c1 = 5;                                % s/m^2 (A<=30)
farea  = @(A) c1*A;                    % 本场景 A=20<=30
b1 = 4.8; b0 = 0.8;
fskill = @(S) b1./S + b0;              % 专业S=2 => 3.2

% 烟雾线性模型（任务二）
Vis   = @(t) 15 - 0.02*t;              % 可视度(m), t单位:s
eta   = @(t) min(1, Vis(t)/8);         % 归一化能见度
kappa = @(e) min(1, e + 0.1);          % 技能发挥系数
room_time = @(alpha,A,S,tstart) alpha .* farea(A) .* fskill(S) .* ...
                               (kappa(eta(tstart)) ./ max(eta(tstart),1e-9));

% 边际效应（任务一）
c0=20; gamma=1.2;                      % 协调 c(x)=c0*x^gamma
theta0=0.3;                            % 冲突 theta(x)=1+theta0*ln(x)
delta0=1.5;                            % 碎片 delta(x)=max(1,delta0*x/m)
c_fun     = @(x) c0 * x.^gamma;
theta_fun = @(x) 1 + theta0*log(x);
delta_fun = @(x,m) max(1, delta0 * x/m);


% 安全度（任务一）
w = [0.4, 0.3, 0.3];                   % [完整性C, 可靠性V, 风险R]
V_low = 50;                            % 单人排查 V=50
Cscore = @(Ks,k0i) 100*min(1, Ks/max(k0i,1));
Rscore = @(alpha,ttot,t0,beta) 100*alpha .* max(0, 1 - max(0,(ttot-t0)./t0).*beta);

t0i   = 210*ones(m,1);                 % 安全时间上限（s）
betaT = ones(m,1);                     % 低风险 β=1
Smin_i= 80*ones(m,1);

% 布局：走廊30m；房门位置；两人从 E1/E2 对向清查
roomPos  = [5;10;15;15;20;25];         % R1..R6 门距 E1 的距离
startPos = [0;30];                     % F1 从 0m，F2 从 30m
assignF{1} = [1 2 3];                  % F1 清查 R1-3
assignF{2} = [4 5 6];                  % F2 清查 R4-6

theta = theta_fun(x);
delta = delta_fun(x,m);

%% ================= 核心计算（详细表格 + 单人耗时） =================
T_total = 0;                           % 预置
Tj = zeros(x,1);                       % 每人总时长
Ks = zeros(m,1);                       % 实际扫描点
V  = V_low*ones(m,1);                  % 可靠性（单人）
R  = zeros(m,1);                       % 风险规避分
tableRows = cell(0,5);                 % 房间详细表：Room, t_start, η, κ, t_ij

for j = 1:x
    pos = startPos(j);
    tcur = 0; moveDist = 0; Tfirst = 0;

    for r = assignF{j}
        % 移动到房门
        d = abs(roomPos(r) - pos);
        moveDist = moveDist + d;
        tcur = tcur + d / v(j);

        % 房间时间（含烟雾修正，按开始时刻 tcur）
        e_now = eta(tcur);
        k_now = kappa(e_now);
        tij   = room_time(alphaTi(r), A(r), Sskill(j), tcur);

        % 统计
        Tfirst = Tfirst + tij;
        Ks(r)  = max(Ks(r), k0(r));
        R(r)   =  Rscore(alphaTi(r), tij, t0i(r), betaT(r));

        % 记录详细表
        tableRows(end+1,:) = {sprintf('R%d', r), tcur, e_now, k_now, tij};

        % 更新位置与时间
        pos  = roomPos(r);
        tcur = tcur + tij;
    end

    % 单人总时长（含冲突/碎片修正）
    Tj(j) = Tfirst*delta + (moveDist / v(j)) * theta;
end

% 并行瓶颈 + 协调
T_total = max(Tj) + c_fun(x);

%% ================= 安全度计算与输出 =================
C = arrayfun(@(i) Cscore(Ks(i),k0(i)), 1:m)';     % 完整性
S = w(1)*C + w(2)*V + w(3)*R;                     % 单房间安全度
Stotal = mean(S);

fprintf('\n--- Detailed Room Table ---\n');
T = cell2table(tableRows, 'VariableNames', {'Room','t_start(s)','eta','kappa','t_ij(s)'});
disp(T);

fprintf('\n--- Room Safety Table ---\n');
RoomIdx = (1:m)';
T_safe = table(RoomIdx, k0, C, V, R, S, Smin_i, ...
    'VariableNames', {'Room','k0','C','V','R','S','Smin'});
disp(T_safe);

fprintf('\n=== Summary ===\n');
fprintf('Number of firefighters x = %d\n', x);
fprintf('===> Total Time = %.2f s (%.2f min)\n', T_total, T_total/60);
fprintf('Per-person Time = [%0.2f , %0.2f] s\n', Tj(1), Tj(2));
fprintf('Stotal = %.\n', Stotal);

%% ================= 人数—总时间曲线（x=1..6） =================
%Xmax = 6; Tcurve = nan(Xmax,1);
%for xx = 1:Xmax
%    Tcurve(xx) = simulate_Ttotal(xx, roomPos, A, alphaTi, Sskill, v, ...
%        room_time, c_fun, theta_fun, delta_fun);
%end
%figure('Name','Ttotal_vs_x'); 
%plot(1:Xmax, Tcurve/60, '-o', 'LineWidth', 1.5); grid on;
%xlabel('Responders (x)'); ylabel('Total time (min)');
%title('T_{total}(x) under Fire Scenario');



end

%% ================= 辅助函数：给定人数 x 计算 T_total =================
function Ttot = simulate_Ttotal(x, roomPos, A, alphaTi, Sskill_base, v_base, ...
                                room_time, c_fun, theta_fun, delta_fun)
m = numel(A);

% 起点（两端入口）；超过2人时，多余人员默认从左端进入（简化）
startPos = zeros(x,1);
startPos(1) = 0; 
if x >= 2, startPos(2) = 30; end
if x > 2, startPos(3:x) = 0; end

% 简单就近/均分分配：将房间序列切成 x 份（左右端交替分配方向，减少折返）
rooms = 1:m;
chunks = cell(x,1);
for j = 1:x
    idx = j:x:m;             % 轮转切分
    chunks{j} = rooms(idx);
end
% 对右端进入的人让其从远到近扫，减少折返
for j = 1:x
    if startPos(j) > 0
        chunks{j} = fliplr(chunks{j});
    end
end

% 人员参数（多于模板则复用第1/2人的能力）
Sskill = zeros(x,1); v = zeros(x,1);
Sskill(:) = Sskill_base(1);
v(:)      = v_base(1);
if numel(Sskill_base) >= 2
    Sskill(1:2) = Sskill_base(1:2);
    v(1:2)      = v_base(1:2);
end

theta = theta_fun(x);
delta = delta_fun(x,m);

Tj = zeros(x,1);
for j = 1:x
    pos = startPos(j); tcur = 0; moveDist = 0; Tfirst = 0;
    seq = chunks{j};
    for r = seq
        d = abs(roomPos(r) - pos);
        moveDist = moveDist + d;
        tcur = tcur + d / v(j);
        tij = room_time(alphaTi(r), A(r), Sskill(j), tcur);
        Tfirst = Tfirst + tij;
        pos = roomPos(r);
        tcur = tcur + tij;
    end
    Tj(j) = Tfirst*delta + (moveDist / v(j)) * theta;
end

Ttot = max(Tj) + c_fun(x);
end
