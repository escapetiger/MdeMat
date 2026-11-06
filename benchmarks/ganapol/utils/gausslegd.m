%% Function: 1D Guass-Legendre quadrature weights and nodes
% Inputs:
%   n - number of quadrature nodes
% Outputs:
%   x - quadrature nodes
%   w - weights
%   T - Jacobi matrix
function [x, w, T] = gausslegd(n, a, b, normalized)
    if nargin < 4
        normalized = true;
        if nargin < 3
            a = -1; b = 1;
        end
    end
    
    validateattributes(n, {'double'}, {'scalar', 'positive', 'integer'});
    validateattributes(a, {'double'}, {'scalar', 'finite', 'real'});
    validateattributes(b, {'double'}, {'scalar', 'finite', 'real'});
    if b < a || (normalized && b == a)
        error('mdemat:ganapol:Interval', 'Require an ordered interval of positive length when normalized.');
    end
    if n == 1
        x = (a+b)/2;
        w = b-a;
        if normalized, w = 1; end
        T = 0;
        return;
    end
    beta = 0.5*(1-(2*(1:n-1)).^(-2)).^(-0.5);
    T = diag(beta, 1) + diag(beta, -1);
    [V, D] = eig(T);
    x = diag(D); 
    w = 2*(V(1, :).^2)';

    % Linear map from[-1,1] to [a,b]
    x=(a*(1-x)+b*(1+x))/2;      
    w=(b-a)/2*w;

    if normalized
        w = w / (b-a);
    end
end
