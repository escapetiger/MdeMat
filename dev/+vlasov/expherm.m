function I = expherm(m, mu, sigma)
% EXPHERM Compute expectation of normalized Hermite polynomials.
%
% m     : vector of polynomial orders
% mu    : mean of Gaussian
% sigma : standard deviation


I = zeros(1, m+1);

I(1) = 1;

if m < 1, return; end

I(2) = mu;

if m < 2, return; end

for k = 2:m
    I(k+1) = (mu * I(k) + (sigma^2 - 1) * sqrt(k-1) * I(k-1)) / sqrt(k);
end

end
