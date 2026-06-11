function [A_mat, P_mat, info] = BHO_three_layer_sepPower(H_sat, H_abs, H_tbs, Pmax_sat, Pmax_abs, Pmax_tbs, sigma2, B, maxIter)
% Three-layer BHO baseline with separate power budgets per layer.
% Inputs:
%   H_sat, H_abs, H_tbs : K x 12 channel gain matrices (nonnegative gains)
%   Pmax_sat, Pmax_abs, Pmax_tbs : scalar power budgets for each layer
%   sigma2 : noise power
%   B      : bandwidth per channel
%   maxIter: BHO iterations
%
% Outputs:
%   A_mat : K x 36 binary assignment matrix
%   P_mat : K x 36 power allocation matrix
%   info  : struct with assignment/power summaries

R_min = 0.1e6;

%R_min = 0.5e6;
numHunters = 60;
epsilon = 0.25;

[K, Nsat] = size(H_sat);
[K2, Nabs] = size(H_abs);
[K3, Ntbs] = size(H_tbs);

if K2 ~= K || K3 ~= K
    error('H_sat, H_abs, and H_tbs must have the same number of users (rows).');
end
if Nsat ~= 12 || Nabs ~= 12 || Ntbs ~= 12
    error('Each layer must have exactly 12 channels.');
end

N = Nsat + Nabs + Ntbs;
if K ~= N
    error('Current permutation-based BHO requires K = 36 users.');
end

H = [H_sat, H_abs, H_tbs];
layer_of_channel = [ones(1,Nsat), 2*ones(1,Nabs), 3*ones(1,Ntbs)];
Pmax_vec = [Pmax_sat, Pmax_abs, Pmax_tbs];

hunters = zeros(numHunters, K);
for i = 1:numHunters
    hunters(i,:) = randperm(N);
end

fitness = zeros(numHunters,1);
for i = 1:numHunters
    fitness(i) = compute_fitness_sepPower(hunters(i,:), H, K, sigma2, B, Pmax_vec, R_min, layer_of_channel);
end

[best_fit, idx] = max(fitness);
best_sol = hunters(idx,:);

for iter = 1:maxIter
    for i = 1:numHunters
        current = hunters(i,:);

        if rand < 0.5
            j = randi(numHunters);
            new_sol = current;
            pos = randperm(K,2);
            new_sol(pos) = hunters(j,pos);
        else
            new_sol = current;
            idx_replace = randperm(K, max(1, round(K/2)));
            new_sol(idx_replace) = best_sol(idx_replace);
        end

        if rand < epsilon
            pos = randperm(K,2);
            new_sol(pos) = new_sol(fliplr(pos));
        end

        new_sol = repair_permutation(new_sol, N);
        new_fit = compute_fitness_sepPower(new_sol, H, K, sigma2, B, Pmax_vec, R_min, layer_of_channel);

        if new_fit > fitness(i)
            hunters(i,:) = new_sol;
            fitness(i) = new_fit;
        end
    end

    [curr_best, idx] = max(fitness);
    if curr_best > best_fit
        best_fit = curr_best;
        best_sol = hunters(idx,:);
    end
end

[p_alloc, rate_alloc] = water_filling_sepPower(best_sol, H, K, sigma2, B, Pmax_vec, R_min, layer_of_channel);

A_mat = zeros(K,N);
P_mat = zeros(K,N);
for k = 1:K
    n = best_sol(k);
    A_mat(k,n) = 1;
    P_mat(k,n) = p_alloc(k);
end

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
end

function fit = compute_fitness_sepPower(sol, H, K, sigma2, B, Pmax_vec, R_min, layer_of_channel)
% [p_alloc, rate_alloc] = water_filling_sepPower(sol, H, K, sigma2, B, Pmax_vec, R_min, layer_of_channel);
% penalty = 0;
% for k = 1:K
%     if rate_alloc(k) < R_min
%         penalty = penalty + 100*(R_min - rate_alloc(k));
%     end
% end
% fit = sum(rate_alloc) - penalty;
[p_alloc, rate_alloc] = water_filling_sepPower(sol, H, K, sigma2, B, Pmax_vec, R_min, layer_of_channel);

% QoS penalty
penalty = 0;
for k = 1:K
    if rate_alloc(k) < R_min
        penalty = penalty + 100*(R_min - rate_alloc(k));
    end
end

% -------- WPF OBJECTIVE --------
epsilon_pf = 1e-9;        % for numerical stability (log)
w = ones(K,1);            % user weights (can be customized)

pf_term = sum(w .* log(rate_alloc + epsilon_pf));

% Final fitness
fit = pf_term - penalty;
end

function [p_alloc, rate_alloc] = water_filling_sepPower(sol, H, K, sigma2, B, Pmax_vec, R_min, layer_of_channel)
p_alloc = zeros(K,1);
rate_alloc = zeros(K,1);

for layer = 1:3
    idx_users = find(layer_of_channel(sol) == layer);
    if isempty(idx_users)
        continue;
    end

    gains = zeros(length(idx_users),1);
    for ii = 1:length(idx_users)
        k = idx_users(ii);
        n = sol(k);
        gains(ii) = max(H(k,n), eps);
    end

    P_layer = Pmax_vec(layer);
    p_layer = waterfill_single_layer(gains, sigma2, P_layer);

    for ii = 1:length(idx_users)
        k = idx_users(ii);
        p_alloc(k) = p_layer(ii);
        rate_alloc(k) = B * log2(1 + p_alloc(k) * gains(ii) / sigma2);
    end

    weak_users = idx_users(rate_alloc(idx_users) < R_min);
    for uu = 1:length(weak_users)
        k = weak_users(uu);
        n = sol(k);
        g = max(H(k,n), eps);
        p_req = (2^(R_min/B) - 1) * sigma2 / g;
        if p_req > p_alloc(k)
            delta = p_req - p_alloc(k);
            candidate_users = idx_users(idx_users ~= k);
            reducible = candidate_users(p_alloc(candidate_users) > 0);

            while delta > 1e-12 && ~isempty(reducible)
                [~, idx_max] = max(p_alloc(reducible));
                k_red = reducible(idx_max);
                take = min(delta, p_alloc(k_red));
                p_alloc(k_red) = p_alloc(k_red) - take;
                p_alloc(k) = p_alloc(k) + take;
                delta = delta - take;
                candidate_users = idx_users(idx_users ~= k);
                reducible = candidate_users(p_alloc(candidate_users) > 1e-12);
            end
        end
    end

    for ii = 1:length(idx_users)
        k = idx_users(ii);
        n = sol(k);
        g = max(H(k,n), eps);
        rate_alloc(k) = B * log2(1 + p_alloc(k) * g / sigma2);
    end

    p_sum = sum(p_alloc(idx_users));
    if p_sum > P_layer + 1e-9
        p_alloc(idx_users) = p_alloc(idx_users) * (P_layer / p_sum);
        for ii = 1:length(idx_users)
            k = idx_users(ii);
            n = sol(k);
            g = max(H(k,n), eps);
            rate_alloc(k) = B * log2(1 + p_alloc(k) * g / sigma2);
        end
    end
end
end

function p = waterfill_single_layer(gains, sigma2, Ptot)
M = length(gains);
if M == 0
    p = [];
    return;
end

inv_terms = sigma2 ./ gains;
[inv_sorted, order] = sort(inv_terms, 'ascend');

mu = 0;
for m = 1:M
    mu_candidate = (Ptot + sum(inv_sorted(1:m))) / m;
    if m == M || mu_candidate <= inv_sorted(m+1)
        mu = mu_candidate;
        break;
    end
end

p_sorted = max(mu - inv_sorted, 0);
p = zeros(M,1);
p(order) = p_sorted;
end

function sol = repair_permutation(sol, N)
used = false(1,N);
duplicates = [];
for i = 1:length(sol)
    if sol(i) >= 1 && sol(i) <= N && ~used(sol(i))
        used(sol(i)) = true;
    else
        duplicates(end+1) = i; %#ok<AGROW>
    end
end
missing = find(~used);
for t = 1:length(duplicates)
    sol(duplicates(t)) = missing(t);
end
end
