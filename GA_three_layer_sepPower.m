% function [A_mat, P_mat, info] = GA_three_layer_sepPower(H_sat, H_abs, H_tbs, Pmax_sat, Pmax_abs, Pmax_tbs, sigma2, B, maxIter)
% 
% %% ================= DIMENSIONS =================
% [K, N_sat] = size(H_sat);
% [~, N_abs] = size(H_abs);
% [~, N_tbs] = size(H_tbs);
% 
% N_total = N_sat + N_abs + N_tbs;
% 
% H = [H_sat, H_abs, H_tbs];   % 36 × 36
% 
% %% ================= GA PARAMETERS =================
% POP = 60;
% G   = maxIter;
% 
% pcross = 0.9;
% pmA = 0.02;
% pmP = 0.02;
% 
% ELITE = max(1, round(0.05 * POP));
% 
% %% ================= INITIALIZATION =================
% pop_assign = zeros(POP, N_total);
% for p = 1:POP
%     pop_assign(p,:) = randperm(K, N_total); % permutation
% end
% 
% pop_power = rand(POP, N_total);
% pop_power = project_power_layers(pop_power, Pmax_sat, Pmax_abs, Pmax_tbs);
% 
% best_fit = -inf;
% 
% %% ================= GA LOOP =================
% for g = 1:G
% 
%     %% ---- FITNESS ----
%     fit = zeros(POP,1);
%     for i = 1:POP
%         fit(i) = compute_fitness_sepPower( ...
%             pop_assign(i,:), pop_power(i,:), ...
%             H, sigma2, B);
%     end
% 
%     %% ---- BEST ----
%     [fmax, idx] = max(fit);
%     if fmax > best_fit
%         best_fit = fmax;
%         best_assign = pop_assign(idx,:);
%         best_power  = pop_power(idx,:);
%     end
% 
%     %% ---- ELITISM ----
%     [~, idx_sort] = sort(fit, 'descend');
%     eliteA = pop_assign(idx_sort(1:ELITE),:);
%     eliteP = pop_power(idx_sort(1:ELITE),:);
% 
%     %% ---- OFFSPRING ----
%     offA = zeros(POP-ELITE, N_total);
%     offP = zeros(POP-ELITE, N_total);
% 
%     for i = 1:2:(POP-ELITE)
%         j = min(i+1, POP-ELITE);
% 
%         p1 = idx_sort(i);
%         p2 = idx_sort(j);
% 
%         A1 = pop_assign(p1,:); A2 = pop_assign(p2,:);
%         P1 = pop_power(p1,:);  P2 = pop_power(p2,:);
% 
%         if rand < pcross
%             [C1A, C2A] = order_crossover_assign(A1, A2);
%             [C1P, C2P] = arithmetic_crossover_power(P1, P2);
%         else
%             C1A = A1; C2A = A2;
%             C1P = P1; C2P = P2;
%         end
% 
%         offA(i,:) = C1A; offP(i,:) = C1P;
%         offA(j,:) = C2A; offP(j,:) = C2P;
%     end
% 
%     %% ---- POWER REPAIR (CRITICAL FOR 3-LAYER) ----
%     offP = project_power_layers(offP, Pmax_sat, Pmax_abs, Pmax_tbs);
% 
%     %% ---- NEXT GEN ----
%     pop_assign = [eliteA; offA];
%     pop_power  = [eliteP; offP];
% end
% 
% %% ================= OUTPUT =================
% A_mat = zeros(K, N_total);
% P_mat = zeros(K, N_total);
% 
% for n = 1:N_total
%     u = best_assign(n);
%     A_mat(u,n) = 1;
%     P_mat(u,n) = best_power(n);
% end
% 
% %% Layer info (for plotting)
% assigned_layer = zeros(K,1);
% for u = 1:K
%     ch = find(A_mat(u,:) > 0.5,1);
%     if ch <= N_sat
%         assigned_layer(u) = 1;
%     elseif ch <= N_sat + N_abs
%         assigned_layer(u) = 2;
%     else
%         assigned_layer(u) = 3;
%     end
% end
% 
% info.best_fit = best_fit;
% info.assigned_layer = assigned_layer;
% 
% end
% 
% function fit = compute_fitness_sepPower(assign, power, H, sigma2, B)
% 
% N = length(assign);
% g = H(sub2ind(size(H), assign, 1:N));
% 
% snr = (power .* g) ./ sigma2;
% fit = sum(B * log2(1 + snr));
% 
% end
% 
% function P = project_power_layers(P, Pmax_sat, Pmax_abs, Pmax_tbs)
% 
% % Split layers
% P_sat = P(:,1:12);
% P_abs = P(:,13:24);
% P_tbs = P(:,25:36);
% 
% % Normalize each layer separately
% P_sat = normalize_layer(P_sat, Pmax_sat);
% P_abs = normalize_layer(P_abs, Pmax_abs);
% P_tbs = normalize_layer(P_tbs, Pmax_tbs);
% 
% P = [P_sat, P_abs, P_tbs];
% 
% end
% 
% function P = normalize_layer(P, Pmax)
% P = max(P,0);
% s = sum(P,2);
% idx = s > 0;
% P(idx,:) = P(idx,:) .* (Pmax ./ s(idx));
% end
% 
% function [C1, C2] = arithmetic_crossover_power(P1, P2)
% alpha = rand(size(P1));
% C1 = alpha .* P1 + (1-alpha).*P2;
% C2 = (1-alpha).*P1 + alpha.*P2;
% end
% 
% function [child1, child2] = order_crossover_assign(p1, p2)
% 
% n = length(p1);
% pts = sort(randperm(n,2));
% c1 = pts(1); c2 = pts(2);
% 
% child1 = zeros(1,n);
% child2 = zeros(1,n);
% 
% child1(c1:c2) = p1(c1:c2);
% child2(c1:c2) = p2(c1:c2);
% 
% rem1 = p2(~ismember(p2, child1));
% rem2 = p1(~ismember(p1, child2));
% 
% idx = [1:c1-1, c2+1:n];
% child1(idx) = rem1;
% child2(idx) = rem2;
% 
% end

function [A_mat, P_mat, info] = GA_three_layer_sepPower( ...
    H_sat, H_abs, H_tbs, ...
    Pmax_sat, Pmax_abs, Pmax_tbs, ...
    sigma2, B, maxIter)

%% ================= QoS =================
R_min = 0.1e6;   % 0.5 Mbps

%% ================= DIMENSIONS =================
[K, N_sat] = size(H_sat);
[~, N_abs] = size(H_abs);
[~, N_tbs] = size(H_tbs);

N_total = N_sat + N_abs + N_tbs;

if K ~= N_total
    error('GA requires K = N_total');
end

H = [H_sat, H_abs, H_tbs];

layer_of_channel = [ ...
    ones(1,N_sat), ...
    2*ones(1,N_abs), ...
    3*ones(1,N_tbs)];

Pmax_vec = [Pmax_sat, Pmax_abs, Pmax_tbs];

%% ================= GA PARAMETERS =================
POP = 60;
G   = maxIter;

pcross = 0.9;
pmA    = 0.02;

ELITE = max(1, round(0.05 * POP));

%% ================= INITIALIZATION =================
pop_assign = zeros(POP, N_total);

for p = 1:POP
    pop_assign(p,:) = randperm(K);
end

best_fit = -inf;

%% ================= GA LOOP =================
for g = 1:G

    %% ---- FITNESS ----
    fit = zeros(POP,1);

    for i = 1:POP
        fit(i) = compute_fitness_sepPower( ...
            pop_assign(i,:), ...
            H, K, sigma2, B, ...
            Pmax_vec, R_min, layer_of_channel);
    end

    %% ---- BEST ----
    [fmax, idx] = max(fit);
    if fmax > best_fit
        best_fit = fmax;
        best_assign = pop_assign(idx,:);
    end

    %% ---- ELITISM ----
    [~, idx_sort] = sort(fit, 'descend');
    eliteA = pop_assign(idx_sort(1:ELITE),:);

    %% ---- OFFSPRING ----
    offA = zeros(POP-ELITE, N_total);

    for i = 1:2:(POP-ELITE)

        j = min(i+1, POP-ELITE);

        p1 = idx_sort(i);
        p2 = idx_sort(j);

        A1 = pop_assign(p1,:);
        A2 = pop_assign(p2,:);

        % ---- CROSSOVER ----
        if rand < pcross
            [C1A, C2A] = order_crossover_assign(A1, A2);
        else
            C1A = A1;
            C2A = A2;
        end

        % ---- MUTATION ----
        C1A = mutate_assign(C1A, pmA);
        C2A = mutate_assign(C2A, pmA);

        offA(i,:) = C1A;
        offA(j,:) = C2A;
    end

    %% ---- REPAIR ----
    for i = 1:size(offA,1)
        offA(i,:) = repair_permutation(offA(i,:), N_total);
    end

    %% ---- NEXT GENERATION ----
    pop_assign = [eliteA; offA];

end

%% ================= FINAL POWER =================
[p_alloc, rate_alloc] = water_filling_sepPower( ...
    best_assign, H, K, sigma2, B, ...
    Pmax_vec, R_min, layer_of_channel);

%% ================= OUTPUT =================
A_mat = zeros(K, N_total);
P_mat = zeros(K, N_total);

for k = 1:K
    ch = best_assign(k);
    A_mat(k, ch) = 1;
    P_mat(k, ch) = p_alloc(k);
end

%% ================= INFO =================
assigned_layer = layer_of_channel(best_assign);

info.best_fit = best_fit;
info.best_assign = best_assign;

info.rate_per_user = rate_alloc;
info.sum_rate = sum(rate_alloc);

info.qos_satisfied = sum(rate_alloc >= R_min);
info.qos_ratio = mean(rate_alloc >= R_min);

info.assigned_layer = assigned_layer;

info.total_power_sat = sum(p_alloc(assigned_layer == 1));
info.total_power_abs = sum(p_alloc(assigned_layer == 2));
info.total_power_tbs = sum(p_alloc(assigned_layer == 3));

info.A_sat = A_mat(:,1:12);
info.A_abs = A_mat(:,13:24);
info.A_tbs = A_mat(:,25:36);

info.P_sat = P_mat(:,1:12);
info.P_abs = P_mat(:,13:24);
info.P_tbs = P_mat(:,25:36);

end

%% ================= HELPER FUNCTIONS =================

function fit = compute_fitness_sepPower(assign, H, K, sigma2, B, ...
                                        Pmax_vec, R_min, layer_of_channel)

[p_alloc, rates] = water_filling_sepPower( ...
    assign, H, K, sigma2, B, ...
    Pmax_vec, R_min, layer_of_channel);

qos_violations = sum(rates < R_min);
penalty = 1e6 * qos_violations;

fit = sum(rates) - penalty;

end

% ----------------------------------------------------

function [p_alloc, rates] = water_filling_sepPower( ...
    assign, H, K, sigma2, B, ...
    Pmax_vec, R_min, layer_of_channel)

p_alloc = zeros(K,1);
rates   = zeros(K,1);

for layer = 1:3

    idx = find(layer_of_channel(assign) == layer);
    if isempty(idx), continue; end

    g = zeros(length(idx),1);
    for i = 1:length(idx)
        u = idx(i);
        ch = assign(u);
        g(i) = H(u,ch);
    end

    Pmax = Pmax_vec(layer);

    % Initial equal allocation
    p = (Pmax / length(idx)) * ones(length(idx),1);

    % QoS iterative adjustment
    for iter = 1:10
        r = B * log2(1 + p .* g / sigma2);
        viol = r < R_min;

        if ~any(viol), break; end

        p(viol) = p(viol) * 1.2;

        % Normalize
        p = max(p,0);
        p = p * (Pmax / sum(p));
    end

    for i = 1:length(idx)
        u = idx(i);
        p_alloc(u) = p(i);
        rates(u) = B * log2(1 + p(i)*g(i)/sigma2);
    end

end

end

% ----------------------------------------------------

function A = mutate_assign(A, pm)

for i = 1:length(A)
    if rand < pm
        j = randi(length(A));
        tmp = A(i);
        A(i) = A(j);
        A(j) = tmp;
    end
end

end

% ----------------------------------------------------

function sol = repair_permutation(sol, N)

missing = setdiff(1:N, sol);
counts = histcounts(sol, 1:N+1);

m = 1;

for i = 1:length(sol)
    if counts(sol(i)) > 1
        sol(i) = missing(m);
        counts(sol(i)) = counts(sol(i)) + 1;
        m = m + 1;

        if m > length(missing)
            break;
        end
    end
end

end

% ----------------------------------------------------

function [child1, child2] = order_crossover_assign(p1, p2)

n = length(p1);

pts = sort(randperm(n,2));
c1 = pts(1);
c2 = pts(2);

child1 = zeros(1,n);
child2 = zeros(1,n);

child1(c1:c2) = p1(c1:c2);
child2(c1:c2) = p2(c1:c2);

rem1 = p2(~ismember(p2, child1));
rem2 = p1(~ismember(p1, child2));

idx = [1:c1-1, c2+1:n];

child1(idx) = rem1;
child2(idx) = rem2;

end

