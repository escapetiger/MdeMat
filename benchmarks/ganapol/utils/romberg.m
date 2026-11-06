function [I, report] = romberg(func, a, b, tol, kmax)
%ROMBERG Richardson-extrapolated trapezoids with a finite refinement budget.
if nargin < 4, tol = 1e-10; end
if nargin < 5, kmax = 15; end
validateattributes(tol, {'double'}, {'scalar', 'positive', '<', 1, 'finite'});
validateattributes(kmax, {'double'}, {'scalar', 'integer', 'positive', '<=', 25});
validateattributes(a, {'double'}, {'scalar', 'real', 'finite'});
validateattributes(b, {'double'}, {'scalar', 'real', 'finite'});
row = (b-a)*(func(a)+func(b))/2;
change = Inf; converged = a == b; level = 0;
for level = 1:kmax
    h = (b-a)/2^level;
    added = 0;
    for index = 1:2^(level-1)
        added = added+func(a+(2*index-1)*h);
    end
    next = zeros(1, level+1);
    next(1) = row(1)/2+h*added;
    for column = 2:level+1
        next(column) = next(column-1)+(next(column-1)-row(column-1))/(4^(column-1)-1);
    end
    change = abs(next(end)-row(end));
    row = next;
    if change <= tol*max(1, abs(row(end)))
        converged = true;
        break;
    end
end
I = row(end);
report = struct('Converged', converged, 'Levels', level, 'LastChange', change);
if ~converged
    warning('mdemat:ganapol:RombergBudget', 'Maximum refinement budget reached without convergence.');
end
end
