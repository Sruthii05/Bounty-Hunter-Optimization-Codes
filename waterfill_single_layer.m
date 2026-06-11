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

