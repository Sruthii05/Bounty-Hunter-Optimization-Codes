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

