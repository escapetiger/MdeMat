%% Function: 1D Guass-Lobatto quadrature weights and nodes
% Inputs:
%   n - number of quadrature nodes
% Outputs:
%   x - quadrature nodes
%   w - weights
function [x, w] = gaussloba(n, a, b, normalized)
    if nargin < 4
        normalized = true;
        if nargin < 3
            a = -1; b = 1;
        end
    end

    x = cos(pi * (0:(n-1)) / (n-1))';
    
    P = zeros(n, n);
    
    xold = ones(n, 1) * 2;
    while max(abs(x - xold)) > eps(1)
        xold = x;
        P(:, 1) = 1;
        P(:, 2) = x;
        for k = 3:n
            P(:, k) = ((2 * k - 3) * x .* P(:, k - 1) - (k - 2) * P(:, k - 2)) / (k-1);
        end
        x = xold - (x .* P(:, n) - P(:, n-1)) ./ (n * P(:, n));
    end
    w = 2 ./ ((n-1) * n * P(:, n) .^ 2);

    % Linear map from[-1,1] to [a,b]
    x=(a*(1-x)+b*(1+x))/2;      
    w=(b-a)/2*w;
    if normalized
        w = w / (b-a);
    end
end




