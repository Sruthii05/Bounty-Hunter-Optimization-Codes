function [A_mat, P_mat, info] = PSO_three_layer_sepPower(H_sat, H_abs, H_tbs,Pmax_sat, Pmax_abs, Pmax_tbs, sigma2, B, maxIter)

%R_min = 0.1;
R_min = 0.1e6;
numParticles = 60;
w  = 0.7;
c1 = 1.0;
c2 = 1.0;

%% ================= DIMENSIONS =================
[K, Nsat] = size(H_sat);
[K2, Nabs] = size(H_abs);
[K3, Ntbs] = size(H_tbs);

if K2 ~= K || K3 ~= K
    error('All channel matrices must have same number of users.');
end

if Nsat ~= 12 || Nabs ~= 12 || Ntbs ~= 12
    error('Each layer must have 12 channels.');
end

N = Nsat + Nabs + Ntbs;

if K ~= N
    error('Permutation-based PSO requires K = N.');
end

%% ================= MERGE CHANNELS =================
H = [H_sat, H_abs, H_tbs];

layer_of_channel = [ ...
    ones(1,Nsat), ...
    2*ones(1,Nabs), ...
    3*ones(1,Ntbs)];

Pmax_vec = [Pmax_sat, Pmax_abs, Pmax_tbs];

%% ================= INITIALIZATION =================
particles  = zeros(numParticles, K);

for i = 1:numParticles
    particles(i,:) = randperm(N);
end

fitness = zeros(numParticles,1);

for i = 1:numParticles
    fitness(i) = compute_fitness_sepPower( ...
        particles(i,:), H, K, sigma2, B, ...
        Pmax_vec, R_min, layer_of_channel);
end

% Personal best
pbest = particles;
pbest_fit = fitness;

% Global best
[gbest_fit, idx] = max(fitness);
gbest = particles(idx,:);

%% ================= MAIN LOOP =================
for iter = 1:maxIter

     for i = 1:numParticles

        current = particles(i,:);
        new_sol = current;

        % -------- DISCRETE PSO UPDATE --------
          for k = 1:K 
        r1 = rand; 
        r2 = rand; 
        % Move toward pbest 
         if r1 < 0.5
         if new_sol(k) ~= pbest(i,k) 
         swap_idx = find(new_sol == pbest(i,k)); 
         new_sol([k swap_idx]) = new_sol([swap_idx k]); 
         end 
        % Move toward gbest 
         elseif r2 < 0.5 
         if new_sol(k) ~= gbest(k) 
         swap_idx = find(new_sol == gbest(k)); 
         new_sol([k swap_idx]) = new_sol([swap_idx k]); 
         end 
         end 
         end 
        % -------- INERTIA (exploration) -------- 
         if rand < w 
             pos = randperm(K,2); 
         new_sol(pos) = new_sol(fliplr(pos)); 
         end
        

        % -------- REPAIR --------
        new_sol = repair_permutation(new_sol, N);

        % -------- FITNESS --------
        new_fit = compute_fitness_sepPower( ...
            new_sol, H, K, sigma2, B, ...
            Pmax_vec, R_min, layer_of_channel);

        % -------- UPDATE --------
        if new_fit > pbest_fit(i)
            pbest(i,:) = new_sol;
            pbest_fit(i) = new_fit;
        end

        particles(i,:) = new_sol;
        fitness(i) = new_fit;

    end

    % -------- GLOBAL BEST --------
    [curr_best, idx] = max(fitness);
    if curr_best > gbest_fit
        gbest_fit = curr_best;
        gbest = particles(idx,:);
    end

end

%% ================= FINAL POWER =================
[p_alloc, rate_alloc] = water_filling_sepPower( ...
    gbest, H, K, sigma2, B, ...
    Pmax_vec, R_min, layer_of_channel);

%% ================= OUTPUT =================
A_mat = zeros(K,N);
P_mat = zeros(K,N);

for k = 1:K
    n = gbest(k);
    A_mat(k,n) = 1;
    P_mat(k,n) = p_alloc(k);
end

%% ================= INFO =================
info.best_fit = gbest_fit;
info.best_sol = gbest;
info.rate_per_user = rate_alloc;
info.sum_rate = sum(rate_alloc);

info.layer_of_channel = layer_of_channel;
info.assigned_layer = layer_of_channel(gbest);

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