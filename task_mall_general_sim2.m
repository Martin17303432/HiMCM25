
function task_mall_general_sim(x)
% 通用大型商场火灾排查模拟
% 支持：
%   - 可变消防员人数 x（默认 4）
%   - 每个房间不同的并行人数 p_room(i)
%   - 每个房间必须由两名不同消防员各自检查一次（主查+复检）
%   - 房间 4 和 5 在二楼，含爬楼梯时间

clc; close all;
fprintf('\n=== General Mall Fire Inspection Simulation ===\n');

if nargin < 1
    x = 6;   % 默认 4 名消防员
end

%% ---------- 房间定义 ----------
% 1: Dining_1 (一楼)
% 2: Dining_2 (一楼)
% 3: Shops    (一楼)
% 4: Kids_1   (二楼)
% 5: Kids_2   (二楼)

m = 5;

A = [100; 100; 200; 150; 150];      % 面积
alphaTi = [1.5; 1.5; 1.0; 1.8; 1.8];
gammaV  = [0.2; 0.2; 0.1; 0.2; 0.2];

% 一楼 / 二楼 标记
floorRoom = [1; 1; 1; 2; 2];

% 可靠性、最小安全度阈值
V      = [75; 75; 75; 100; 100];
Smin_i = [80; 80; 80; 85; 85];

% t0(i) 安全时间窗（秒）
t0i    = [900; 900; 900; 1000; 1000];

% β(Ti) 风险系数
betaT  = [1.5; 1.5; 1.0; 2.0; 2.0];

% 预设扫描点 k0
k0 = floor(0.1*A.*alphaTi + 0.2*gammaV);

%% ---------- 每个房间允许的并行人数 p_room ----------
% 例子：前两个房间并行 2 人，商店 1 人，儿童区并行 3 人
p_room = [2; 4; 2; 2; 2];

%% ---------- 人员参数 ----------
Sskill = 4*ones(x,1);              % 技能值
b1 = 4.8; b0 = 0.8;
fskill = @(S) b1./S + b0;

v_walk  = 2*ones(x,1);             % 水平走动速度 (m/s)
v_stair = v_walk;                  % 上下楼梯速度（先取一样）

%% ---------- 烟雾模型 ----------
V0      = 10;
k_smoke = 0.001;

Vis   = @(t) V0*exp(-k_smoke*t);
eta   = @(t) min(1, Vis(t)/3);
kappa = @(e) min(1, e+0.2);

% f_area(A) = 5*A
farea = @(Ai) 5*Ai;

% 单人房间排查时间
room_time = @(alpha,Ai,S,tstart) ...
    alpha .* farea(Ai).*fskill(S).*(kappa(eta(tstart))./max(eta(tstart),1e-9));

%% ---------- 边际效应 ----------
c0 = 20; gamma = 1.2; theta0 = 0.3; delta0 = 1.5;

c_fun     = @(xx) c0 * xx.^gamma;
theta_fun = @(xx) 1 + theta0*log(xx);
delta_fun = @(xx,mm) max(1, delta0 * xx/mm);

% 安全度
w = [0.4, 0.3, 0.3];   % [C, V, R]
Cscore = @(Ks,k0i) 100*min(1, Ks/max(k0i,1));
Rscore = @(alpha,ttot,t0,beta) ...
    100*alpha .* max(0, 1 - max(0,(ttot-t0)./t0).*beta);

%% ---------- 布局：走廊 + 楼梯 ----------
doorPos   = [  5;  15;  25;  10;  20];  % 每个房间在本楼层走廊的位置
entrPos   = 0;                           % 一楼入口
stairsPos = 12;                          % 楼梯口在走廊上的位置
stairsLen = 3;                           % 楼梯等效长度（垂直）

% 所有人从一楼入口出发
startPosX  = entrPos * ones(x,1);
startFloor = ones(x,1);                  % 都在 1 楼

%% ---------- 任务分配：每个房间两名不同的人检查 ----------
%% ---------- 任务分配：结合 p_room 和 “至少两人复检” ----------
if x < 2
    error('至少需要 2 名消防员才能做到“不同的人复检一次”。');
end

% 每个房间需要的“不同消防员人数” = max(2, p_room(i))
visits_per_room = max(2, p_room(:));   % 列向量
total_visits = sum(visits_per_room);

% 构造一个房间序列，例如：
% p_room = [2;2;1;3;3] → visits_per_room = [2;2;2;3;3]
% roomList = [1 1 2 2 3 3 4 4 4 5 5 5]
roomList = [];
for r = 1:m
    roomList = [roomList,  r*ones(1, visits_per_room(r))];
end

% 把这些“房间任务”轮流分配给 x 个消防员，使负载尽量均衡
assignF = cell(x,1);
for k = 1:total_visits
    j = mod(k-1, x) + 1;     % 第 k 个任务给第 j 个消防员
    r = roomList(k);
    assignF{j}(end+1) = r;
end

% 打印每个消防员要去的房间
fprintf('\n--- Task assignment (room list per firefighter) ---\n');
for j = 1:x
    fprintf('Firefighter %d: ', j);
    if isempty(assignF{j})
        fprintf('(no rooms)');
    else
        fprintf('R%d ', assignF{j});
    end
    fprintf('\n');
end

% 额外打印：每个房间被分配到了哪些消防员（体现 p_room）
roomMembers = cell(m,1);
for j = 1:x
    for r = assignF{j}
        roomMembers{r}(end+1) = j;
    end
end

fprintf('\n--- Room -> firefighters (distinct heads per room) ---\n');
for r = 1:m
    u = unique(roomMembers{r});
    fprintf('Room %d: firefighters ', r);
    fprintf('%d ', u);
    fprintf('(count = %d, target p_room = %d)\n', numel(u), p_room(r));
end


fprintf('\n--- Task assignment (room list per firefighter) ---\n');
for j = 1:x
    fprintf('Firefighter %d: ', j);
    if isempty(assignF{j})
        fprintf('(no rooms)');
    else
        fprintf('R%d ', assignF{j});
    end
    fprintf('\n');
end

theta = theta_fun(x);
delta = delta_fun(x,m);

%% ---------- 核心模拟 ----------
Tj = zeros(x,1);
Ks = zeros(m,1);
R  = zeros(m,1);

tableRows = cell(0,5);   % {Room, t_start, eta, kappa, t_ij}

for j = 1:x
    posX     = startPosX(j);
    curFloor = startFloor(j);
    tcur     = 0;
    moveDist = 0;
    Tfirst   = 0;

    seq = assignF{j};
    for idx = 1:numel(seq)
        r = seq(idx);

        % --------- 计算移动时间（含楼梯） ---------
        targetFloor = floorRoom(r);
        targetPosX  = doorPos(r);

        if curFloor == targetFloor
            % 同楼层：水平走廊
            horiz = abs(targetPosX - posX);
            t_move = horiz / v_walk(j);
            moveDist = moveDist + horiz;
        else
            % 不同楼层：走到楼梯 → 上/下楼 → 再走到房门
            horiz1 = abs(stairsPos - posX);
            horiz2 = abs(targetPosX - stairsPos);
            t_stair = stairsLen / v_stair(j);

            t_move   = (horiz1 + horiz2)/v_walk(j) + t_stair;
            moveDist = moveDist + (horiz1 + horiz2);   % 楼梯长度不计入水平
            curFloor = targetFloor;
        end

        tcur = tcur + t_move;

        % --------- 房间开始时间 ---------
        t_start = tcur;

        % 烟雾 & 技能发挥
        e_now = eta(t_start);
        k_now = kappa(e_now);

        % 单人时间
        t_single = room_time(alphaTi(r), A(r), Sskill(j), t_start);

        % 并行加速
        p   = max(p_room(r), 1);     % 防止 0
        tij = t_single / p;

        % 时间 & 风险统计
        Tfirst = Tfirst + tij;
        Ks(r)  = max(Ks(r), k0(r));
        R(r)   = max(R(r), Rscore(alphaTi(r), tij, t0i(r), betaT(r)));

        % 记录详细表
        tableRows(end+1,:) = { ...
            sprintf('R%d', r), ...
            t_start, ...
            e_now, ...
            k_now, ...
            tij };

        % 结束后位置
        posX = targetPosX;
        tcur = tcur + tij;
    end

    Tj(j) = Tfirst*delta + (moveDist / v_walk(j)) * theta;
end

T_total = max(Tj) + c_fun(x);

%% ---------- 安全度 ----------
C = arrayfun(@(i) Cscore(Ks(i), k0(i)),1:m)';
S = w(1)*C + w(2)*V + w(3)*R;
Stotal = mean(S);

%% ---------- 输出 ----------
fprintf('\n=== Detailed Room Table ===\n');
T_detail = cell2table(tableRows, ...
    'VariableNames', {'Room','t_start_s','eta','kappa','t_ij_s'});
disp(T_detail);

fprintf('\n=== Safety Table ===\n');
RoomIdx = (1:m)';
T_safe = table(RoomIdx, k0, C, V, R, S, Smin_i, ...
    'VariableNames', {'Room','k0','C','V','R','S','Smin'});
disp(T_safe);

fprintf('\n=== SUMMARY ===\n');
fprintf('Firefighters x = %d\n', x);
fprintf('Per-person time Tj (s): '); disp(Tj');
fprintf('T_total = %.1f s (%.2f min)\n', T_total, T_total/60);
fprintf('Stotal  = %.1f\n', Stotal);

end
