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