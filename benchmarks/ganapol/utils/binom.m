function res = binom(k, j)
%BINOM Pairwise binomial coefficients; j>k gives zero.
validateattributes(k, {'double'}, {'vector', 'integer', 'nonnegative', 'finite'});
validateattributes(j, {'double'}, {'vector', 'integer', 'nonnegative', 'finite'});
k = k(:); j = j(:).';
valid = j <= k;
res = zeros(numel(k), numel(j));
values = exp(gammaln(k+1)-gammaln(j+1)-gammaln(max(k-j, 0)+1));
res(valid) = round(values(valid));
if isvector(res), res = res(:).'; end
end
