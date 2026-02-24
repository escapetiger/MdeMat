function res = binom(k,j)
    m = numel(k); n = numel(j);
    k = k(:); j = j(:);
    diff = k-j'; diff(diff<0) = 0;
    fk = factorial(k);
    fj = factorial(j);
    fkj = factorial(diff);
    res = fk./(fj'.*fkj);
    if m>1 && n>1
        res = tril(res);
    elseif m == 1 || n == 1
        res = res(:).';
    end
end
