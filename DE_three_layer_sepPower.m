function [A_mat, P_mat, info] = DE_three_layer_sepPower( ...
    H_sat, H_abs, H_tbs, ...
    Pmax_sat, Pmax_abs, Pmax_tbs, ...
    sigma2, B, maxIter)

%% ================= PARAMETERS =================
R_min = 0.1e6;        % QoS (bps)

NP = 60;              % population size
F  = 0.5;             % mutation factor
CR = 0.8;             % crossover rate

%% ================= DIMENSIONS =================
[K, Nsat] = size(H_sat);
[~, Nabs] = size(H_abs);
[~, Ntbs] = size(H_tbs);

N = Nsat + Nabs + Ntbs;

if K ~= N
    error('DE requires K = N (permutation-based)');
end

%% ================= MERGE CHANNELS =================
H = [H_sat, H_abs, H_tbs];

layer_of_channel = [ ...
    ones(1,Nsat), ...
    2*ones(1,Nabs), ...
    3*ones(1,Ntbs)];

Pmax_vec = [Pmax_sat, Pmax_abs, Pmax_tbs];

%% ================= INITIALIZATION =================
pop = zeros(NP, K);

for i = 1:NP
    pop(i,:) = randperm(N);
end

fit = zeros(NP,1);

for i = 1:NP
    fit(i) = compute_fitness_sepPower( ...
        pop(i,:), H, K, sigma2, B, ...
        Pmax_vec, R_min, layer_of_channel);
end

[best_fit, idx] = max(fit);
best_sol = pop(idx,:);

%% ================= MAIN LOOP =================
for iter = 1:maxIter

    for i = 1:NP

        % -------- SELECT r1, r2, r3 --------
        idxs = randperm(NP,3);
        while any(idxs == i)
            idxs = randperm(NP,3);
        end

        r1 = pop(idxs(1),:);
        r2 = pop(idxs(2),:);
        r3 = pop(idxs(3),:);

        % -------- MUTATION (Permutation-based) --------
        mutant = r1;

        for k = 1:K
            if rand < F
                mutant(k) = r2(k);
            end
        end

        % -------- REPAIR (important) --------
        mutant = repair_permutation(mutant, N);

        % -------- CROSSOVER --------
        trial = pop(i,:);

        for k = 1:K
            if rand < CR
                trial(k) = mutant(k);
            end
        end

        trial = repair_permutation(trial, N);

        % -------- FITNESS --------
        trial_fit = compute_fitness_sepPower( ...
            trial, H, K, sigma2, B, ...
            Pmax_vec, R_min, layer_of_channel);

        % -------- SELECTION --------
        if trial_fit > fit(i)
            pop(i,:) = trial;
            fit(i) = trial_fit;
        end

    end

    % -------- GLOBAL BEST --------
    [curr_best, idx] = max(fit);
    if curr_best > best_fit
        best_fit = curr_best;
        best_sol = pop(idx,:);
    end

end

%% ================= FINAL POWER =================
[p_alloc, rate_alloc] = water_filling_sepPower( ...
    best_sol, H, K, sigma2, B, ...
    Pmax_vec, R_min, layer_of_channel);

%% ================= OUTPUT =================
A_mat = zeros(K,N);
P_mat = zeros(K,N);

for k = 1:K
    n = best_sol(k);
    A_mat(k,n) = 1;
    P_mat(k,n) = p_alloc(k);
end

%% ================= INFO =================
info.best_fit = best_fit;
info.best_sol = best_sol;
info.rate_per_user = rate_alloc;
info.sum_rate = sum(rate_alloc);

info.layer_of_channel = layer_of_channel;
info.assigned_layer = layer_of_channel(best_sol);

info.total_power_sat = sum(p_alloc(info.assigned_layer == 1));
info.total_power_abs = sum(p_alloc(info.assigned_layer == 2));
info.total_power_tbs = sum(p_alloc(info.assigned_layer == 3));

info.A_sat = A_mat(:,1:12);
info.A_abs = A_mat(:,13:24);
info.A_tbs = A_mat(:,25:36);

info.P_sat = P_mat(:,1:12);
info.P_abs = P_mat(:,13:24);
info.P_tbs = P_mat(:,25:36);

info.qos_satisfied = sum(rate_alloc >= R_min);
info.qos_ratio = mean(rate_alloc >= R_min);

end