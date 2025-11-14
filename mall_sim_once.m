function T_total = mall_sim_once(x, p_room, verbose)
% 单次商场火灾排查仿真
% 输入:
%   x       - 消防员人数
%   p_room  - 5x1 向量, 第 i 个元素是房间 i 的并行人数设定
%   verbose - (可选) 是否打印详细表, 默认 true
% 输出:
%   T_total - 协调后的总排查时间 (s)

if nargin < 3
    verbose = true;
end

fprintf('\n=== Mall Simulation: x = %d, p_room = [%s] ===\n', ...
        x, num2str(p_room(:).'));

%% ---------- 房间定义 ----------
m = 5;
A = [100; 100; 200; 150; 150];
alphaTi = [1.5; 1.5; 1.0; 1.8; 1.8];
gammaV  = [0.2; 0.2; 0.1; 0.2; 0.2];
floorRoom = [1; 1; 1; 2; 2];

V      = [75; 75; 75; 100; 100];
Smin_i = [70; 70; 70; 75; 75];

t0i    = [900; 900; 900; 600; 600];
betaT  = [1.5; 1.5; 1.0; 2.0; 2.0];

k0 = floor(0.1*A.*alphaTi + 0.2*gammaV);

%% ---------- 人员参数 ----------
Sskill = 4*ones(x,1);
b1 = 4.8; b0 = 0.8;
fskill = @(S) b1./S + b0;

v_walk  = 2*ones(x,1);
v_stair = v_walk;

%% ---------- 烟雾模型 ----------
V0      = 10;
k_smoke = 0.001;
Vis   = @(t) V0*exp(-k_smoke*t);
eta   = @(t) min(1, Vis(t)/3);
kappa = @(e) min(1, e+0.2);

farea      = @(Ai) 5*Ai;
room_time  = @(alpha,Ai,S,tstart) ...
    alpha .* farea(Ai).*fskill(S).*(kappa(eta(tstart))./max(eta(tstart),1e-9));

%% ---------- 边际效应 ----------
c0 = 20; gamma = 1.2; theta0 = 0.3; delta0 = 1.5;
c_fun     = @(xx) c0 * xx.^gamma;
theta_fun = @(xx) 1 + theta0*log(xx);
delta_fun = @(xx,mm) max(1, delta0 * xx/mm);

w      = [0.4, 0.3, 0.3];
Cscore = @(Ks,k0i) 100*min(1, Ks/max(k0i,1));
Rscore = @(alpha,ttot,t0,beta) ...
    100*alpha .* max(0, 1 - max(0,(ttot-t0)./t0).*beta);

%% ---------- 布局 ----------
doorPos   = [  5;  15;  25;  10;  20];
entrPos   = 0;
stairsPos = 12;
stairsLen = 3;

startPosX  = entrPos * ones(x,1);
startFloor = ones(x,1);

%% ---------- 任务分配, 用 p_room ----------
if x < 2
    error('至少需要 2 名消防员才能做到“不同的人复检一次”。');
end

visits_per_room = max(2, p_room(:));      % 每个房间至少 2 人
total_visits    = sum(visits_per_room);

roomList = [];
for r = 1:m
    roomList = [roomList, r*ones(1, visits_per_room(r))]; %#ok<AGROW>
end

assignF = cell(x,1);
for k = 1:total_visits
    j = mod(k-1, x) + 1;
    r = roomList(k);
    assignF{j}(end+1) = r;
end

theta = theta_fun(x);
delta = delta_fun(x,m);

%% ---------- 核心仿真 ----------
Tj = zeros(x,1);
Ks = zeros(m,1);
R  = zeros(m,1);

tableRows = cell(0,5);

for j = 1:x
    posX     = startPosX(j);
    curFloor = startFloor(j);
    tcur     = 0;
    moveDist = 0;
    Tfirst   = 0;

    seq = assignF{j};
    for idx = 1:numel(seq)
        r = seq(idx);

        targetFloor = floorRoom(r);
        targetPosX  = doorPos(r);

        if curFloor == targetFloor
            horiz   = abs(targetPosX - posX);
            t_move  = horiz / v_walk(j);
            moveDist = moveDist + horiz;
        else
            horiz1  = abs(stairsPos - posX);
            horiz2  = abs(targetPosX - stairsPos);
            t_stair = stairsLen / v_stair(j);

            t_move  = (horiz1 + horiz2)/v_walk(j) + t_stair;
            moveDist = moveDist + (horiz1 + horiz2);
            curFloor = targetFloor;
        end

        tcur = tcur + t_move;

        t_start = tcur;
        e_now   = eta(t_start);
        k_now   = kappa(e_now);

        t_single = room_time(alphaTi(r), A(r), Sskill(j), t_start);
        p        = max(p_room(r),1);
        tij      = t_single / p;

        Tfirst = Tfirst + tij;
        Ks(r)  = max(Ks(r), k0(r));
        R(r)   = max(R(r), Rscore(alphaTi(r), tij, t0i(r), betaT(r)));

        tableRows(end+1,:) = {sprintf('R%d', r), t_start, e_now, k_now, tij};

        posX = targetPosX;
        tcur = tcur + tij;
    end

    Tj(j) = Tfirst*delta + (moveDist / v_walk(j)) * theta;
end

T_total = max(Tj) + c_fun(x);

%% ---------- 输出 ----------
if verbose
    % 房间 -> 消防员
    roomMembers = cell(m,1);
    for j = 1:x
        for r = assignF{j}
            roomMembers{r}(end+1) = j; %#ok<AGROW>
        end
    end
    fprintf('\n--- Room -> firefighters ---\n');
    for r = 1:m
        u = unique(roomMembers{r});
        fprintf('Room %d: firefighters ', r);
        fprintf('%d ', u);
        fprintf('(count = %d, target p_room = %d)\n', numel(u), p_room(r));
    end

    fprintf('\n=== Detailed Room Table ===\n');
    T_detail = cell2table(tableRows, ...
        'VariableNames', {'Room','t_start_s','eta','kappa','t_ij_s'});
    disp(T_detail);

    C = arrayfun(@(i) Cscore(Ks(i), k0(i)),1:m)';
    S = w(1)*C + w(2)*V + w(3)*R;
    RoomIdx = (1:m)';
    T_safe = table(RoomIdx, k0, C, V, R, S, Smin_i, ...
        'VariableNames', {'Room','k0','C','V','R','S','Smin'});
    fprintf('\n=== Safety Table ===\n');
    disp(T_safe);

    Stotal = mean(S);
    fprintf('\n=== SUMMARY (single run) ===\n');
    fprintf('T_total = %.1f s (%.2f min)\n', T_total, T_total/60);
    fprintf('Stotal  = %.1f\n', Stotal);
end

end
