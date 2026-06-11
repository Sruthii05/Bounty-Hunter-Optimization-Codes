function [A_mat, P_mat, info] = GWO_three_layer_sepPower( ...
    H_sat, H_abs, H_tbs, ...
    Pmax_sat, Pmax_abs, Pmax_tbs, ...
    sigma2, B, maxIter)

%% ================= QoS =================
R_min = 0.5e6;

%% ================= DIMENSIONS =================
[K, Nsat] = size(H_sat);
[~, Nabs] = size(H_abs);
[~, Ntbs] = size(H_tbs);

if Nsat~=12 || Nabs~=12 || Ntbs~=12
    error('Each layer must have 12 channels');
end

N = Nsat + Nabs + Ntbs;

if K ~= N
    error('Requires K = N');
end

H = [H_sat, H_abs, H_tbs];
layer_of_channel = [ones(1,Nsat), 2*ones(1,Nabs), 3*ones(1,Ntbs)];
Pmax_vec = [Pmax_sat, Pmax_abs, Pmax_tbs];

%% ================= GWO PARAMETERS =================
numWolves = 50;

%% ================= INITIALIZATION =================
wolves = zeros(numWolves, N);

for i = 1:numWolves
    wolves(i,:) = randperm(N);
end

fitness = zeros(numWolves,1);

for i = 1:numWolves
    fitness(i) = compute_fitness_sepPower( ...
        wolves(i,:), H, K, sigma2, B, ...
        Pmax_vec, R_min, layer_of_channel);
end

%% Identify alpha, beta, delta
[fitness, idx_sort] = sort(fitness,'descend');
wolves = wolves(idx_sort,:);

alpha = wolves(1,:);
beta  = wolves(2,:);
delta = wolves(3,:);

alpha_fit = fitness(1);

%% ================= MAIN LOOP =================
for iter = 1:maxIter

    a = 2 - 2*(iter/maxIter);   % exploration → exploitation

    for i = 1:numWolves

        X = wolves(i,:);

        % --- GWO position updates ---
        X1 = gwo_update(X, alpha, a);
        X2 = gwo_update(X, beta,  a);
        X3 = gwo_update(X, delta, a);

        % --- Weighted combination (IMPORTANT IMPROVEMENT) ---
        temp = 0.5*X1 + 0.3*X2 + 0.2*X3;

        % --- Convert to valid permutation (CRITICAL FIX) ---
        [~, new_sol] = sort(temp);

        % --- Fitness evaluation ---
        new_fit = compute_fitness_sepPower( ...
            new_sol, H, K, sigma2, B, ...
            Pmax_vec, R_min, layer_of_channel);

        % --- Greedy selection ---
        if new_fit > fitness(i)
            wolves(i,:) = new_sol;
            fitness(i) = new_fit;
        end
    end

    % --- Update alpha, beta, delta ---
    [fitness, idx_sort] = sort(fitness,'descend');
    wolves = wolves(idx_sort,:);

    alpha = wolves(1,:);
    beta  = wolves(2,:);
    delta = wolves(3,:);

    alpha_fit = fitness(1);
end

best_sol = alpha;

%% ================= POWER ALLOCATION =================
[p_alloc, rate_alloc] = water_filling_sepPower( ...
    best_sol, H, K, sigma2, B, ...
    Pmax_vec, R_min, layer_of_channel);

%% ================= OUTPUT =================
A_mat = zeros(K,N);
P_mat = zeros(K,N);

for k = 1:K
    ch = best_sol(k);
    A_mat(k,ch) = 1;
    P_mat(k,ch) = p_alloc(k);
end

%% ================= INFO =================
info.best_fit = alpha_fit;
info.best_sol = best_sol;
info.sum_rate = sum(rate_alloc);
info.rate_per_user = rate_alloc;
info.assigned_layer = layer_of_channel(best_sol);

info.total_power_sat = sum(p_alloc(info.assigned_layer==1));
info.total_power_abs = sum(p_alloc(info.assigned_layer==2));
info.total_power_tbs = sum(p_alloc(info.assigned_layer==3));

end

function Xnew = gwo_update(X, leader, a)

N = length(X);

r1 = rand(1,N);
r2 = rand(1,N);

A = 2*a*r1 - a;
C = 2*r2;

D = abs(C .* leader - X);

Xnew = leader - A .* D;

end