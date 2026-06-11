% function fit = compute_fitness_sepPower(sol, H, K, sigma2, B, Pmax_vec, R_min, layer_of_channel)
% 
% [p_alloc, rate_alloc] = water_filling_sepPower(sol, H, K, sigma2, B, Pmax_vec, R_min, layer_of_channel);
% 
% penalty = 0;
% for k = 1:K
%     if rate_alloc(k) < R_min
%         penalty = penalty + 100*(R_min - rate_alloc(k));
%     end
% end
% 
% fit = sum(rate_alloc) - penalty;
% 
% end

function fit = compute_fitness_sepPower(sol, H, K, sigma2, B, Pmax_vec, R_min, layer_of_channel)

    % --- power allocation ---
    [p_alloc, rate] = water_filling_sepPower(sol, H, K, sigma2, B, Pmax_vec, R_min, layer_of_channel);

    % --- QoS penalty ---
    violation = max(0, R_min - rate);   % per-user deficit

    penalty = sum(violation);

    lambda = 200;   % tuning parameter

    % --- final fitness ---
    fit = sum(rate) - lambda * penalty;

end

