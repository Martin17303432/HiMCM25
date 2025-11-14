function task3_full_sim5()

% Task 3: Hui-type high school building fire inspection (5 responders)

clc; close all;
fprintf('\n=== Task 3 Simulation: Hui-type High School Building (5 firefighters) ===\n');

%% ============ Scene parameters =============
m = 10;                              % 8 classrooms + 2 labs

% 1..8: classrooms, 9: chemistry lab, 10: physics lab
isClass = true(m,1);
isClass(9:10) = false;

% Area A_i
A = zeros(m,1);
A(isClass)  = 20;        % classroom 20 m^2
A(~isClass) = 30;        % lab 40 m^2

% alpha(T_i)
alphaTi = zeros(m,1);
alphaTi(isClass)  = 1.0;
alphaTi(~isClass) = 1.8;

% gamma(V_i)
gammaV = zeros(m,1);
gammaV(isClass)  = 0.1;
gammaV(~isClass) = 0.2;

% preset scan points k0
k0 = floor(0.1*A.*alphaTi + 0.2*gammaV);   % classroom 2, lab 7

% safe time window t0(i) (s)
t0i = zeros(m,1);
t0i(isClass)  = 30;
t0i(~isClass) = 45;

% beta(T_i)
betaT = zeros(m,1);
betaT(isClass)  = 1.5;
betaT(~isClass) = 2.0;

% min safety S_min(i) (for info)
Smin_i = zeros(m,1);
Smin_i(isClass)  = 80;
Smin_i(~isClass) = 85;

% reliability V(i)
V = zeros(m,1);
V(isClass)  = 75;        % classroom: main+recheck
V(~isClass) = 100;       % lab: double-check

%% ============ Number of responders (3) =============
x = 5;                                % 3 firefighters

Sskill = 4*ones(x,1);                 % all skill=2
b1 = 4.8; b0 = 0.8;
fskill = @(S) b1./S + b0;             % skill factor
v_walk = 2*ones(x,1);               % walking speed (m/s)

%% ============ Smoke model ============
V0      = 10;                         % initial visibility
k_smoke = 0.001;                      % decay

Vis   = @(t) V0 .* exp(-k_smoke .* t);
eta   = @(t) min(1, Vis(t)/3);
kappa = @(e) min(1, e + 0.2);

farea = @(Ai) 5*Ai;                   % 5 * area
room_time = @(alpha,Ai,S,tstart) ...
    alpha .* farea(Ai) .* fskill(S) .* ...
    (kappa(eta(tstart)) ./ max(eta(tstart),1e-9));

%% ============ Marginal effects ============
c0 = 20; gamma = 1.2; theta0 = 0.3; delta0 = 1.5;

c_fun     = @(xx) c0 * xx.^gamma;
theta_fun = @(xx) 1 + theta0*log(xx);
delta_fun = @(xx,mm) max(1, delta0 * xx/mm);

%% ============ Safety score ============
w = [0.4, 0.3, 0.3];                  % [C,V,R] weights

Cscore = @(Ks,k0i) 100*min(1, Ks./max(k0i,1));
Rscore = @(alpha,ttot,t0,beta) ...
    100*max(0, 1-(ttot - t0)./t0.*beta);

t0i(isClass)   = 180; 
t0i(~isClass)   = 450; % 安全时间上限（s）



%% ============ Layout (corridor & doors) =============
Lcorridor = 40;                       % total corridor length
roomPos   = linspace(5, Lcorridor-5, m)';   % door positions

% start positions of 3 firefighters
startPos = zeros(x,1);
startPos(1) = 0;          % F1 from E1
startPos(2) = Lcorridor;  % F2 from E2
startPos(3) = 0;          % F3 from E1
startPos(4) = Lcorridor;          % F3 from E1
startPos(5) = 0;
%% ============ Assignment (3-person realistic) ============
% F1: classrooms 1~4
% F2: classrooms 8~5
% F3: labs 9 & 10
assignF = cell(x,1);
assignF{1} = 1:4;
assignF{2} = [5 6 7 10];
assignF{3} = [1 2 3 9];
assignF{4} = [4 5 8 9];
assignF{5} = [6 7 8 10];
theta = theta_fun(x);
delta = delta_fun(x,m);

%% ============ Core calculation (detailed table) ============
Tj = zeros(x,1);                      % total time per responder
Ks = zeros(m,1);                      % actual scan points
R  = zeros(m,1);                      % risk score
tableRows = cell(0,5);                % Room, t_start, eta, kappa, t_ij

V(isClass)  = 75;        % classroom: main+recheck
V(~isClass) = 100;       % lab: double-check                % 可靠性（单人）

for j = 1:x
    pos = startPos(j);
    tcur = 0;
    moveDist = 0;
    Tfirst = 0;
    seq = assignF{j};

    for r = seq
        % move to room door
        d = abs(roomPos(r) - pos);
        moveDist = moveDist + d;
        tcur = tcur + d / v_walk(j);

        % room inspection
        e_now = eta(tcur);
        k_now = kappa(e_now);
        tij   = room_time(alphaTi(r), A(r), Sskill(j), tcur);

        Tfirst = Tfirst + tij;
        Ks(r)  = max(Ks(r), k0(r));
        R(r)   = max(R(r), Rscore(alphaTi(r), tij, t0i(r), betaT(r)));

        tableRows(end+1,:) = {sprintf('R%d', r), tcur, e_now, k_now, tij};

        pos  = roomPos(r);
        tcur = tcur + tij;
    end

    Tj(j) = Tfirst*delta + (moveDist / v_walk(j)) * theta;
end

T_total = max(Tj) + c_fun(x);

%% ============ Safety result ============
C = arrayfun(@(i) Cscore(Ks(i), k0(i)), 1:m)';   % integrity
S = w(1)*C + w(2)*V + w(3)*R;                    % room safety
Stotal = mean(S);

%% ============ Output =====================
fprintf('\n--- Room Detailed Table ---\n');
T_detail = cell2table(tableRows, ...
    'VariableNames', {'Room','t_start_s','eta','kappa','t_ij_s'});
disp(T_detail);

fprintf('\n--- Room Safety Table ---\n');
RoomIdx = (1:m)';
T_safe = table(RoomIdx, k0, C, V, R, S, Smin_i, ...
    'VariableNames', {'Room','k0','C','V','R','S','Smin'});
disp(T_safe);

fprintf('\n=== Summary ===\n');
fprintf('Number of firefighters x = %d\n', x);
fprintf('Per-responder time Tj = [');
fprintf(' %.1f', Tj);
fprintf(' ] s\n');
fprintf('Total time T_total = %.1f s (%.2f min)\n', T_total, T_total/60);




end

%% ===== Helper: compute T_total for given x ======================
function Ttot = simulate_Ttotal_task3(x, roomPos, A, alphaTi, ...
                                      Sskill_base, v_base, ...
                                      room_time, c_fun, theta_fun, delta_fun)
m = numel(A);

Lcorridor = 40;
startPos = zeros(x,1);
startPos(1) = 0;
if x >= 2, startPos(2) = Lcorridor; end
if x > 2,  startPos(3:x) = 0;       end

rooms = 1:m;
chunks = cell(x,1);
for j = 1:x
    idx = j:x:m;
    chunks{j} = rooms(idx);
end
for j = 1:x
    if startPos(j) > 0
        chunks{j} = fliplr(chunks{j});
    end
end

Sskill = Sskill_base * ones(x,1);
v_walk = v_base * ones(x,1);

theta = theta_fun(x);
delta = delta_fun(x,m);

Tj = zeros(x,1);
for j = 1:x
    pos = startPos(j);
    tcur = 0; moveDist = 0; Tfirst = 0;
    seq = chunks{j};
    for r = seq
        d = abs(roomPos(r) - pos);
        moveDist = moveDist + d;
        tcur = tcur + d / v_walk(j);

        tij = room_time(alphaTi(r), A(r), Sskill(j), tcur);
        Tfirst = Tfirst + tij;

        pos  = roomPos(r);
        tcur = tcur + tij;
    end
    Tj(j) = Tfirst*delta + (moveDist / v_walk(j)) * theta;
end

Ttot = max(Tj) + c_fun(x);
end
