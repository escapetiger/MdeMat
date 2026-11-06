function [L, x, boundaryLift, h] = secondDerivative(n, bounds, options)
%SECONDDERIVATIVE Centered physical second derivative on a uniform 1D grid.
%   Dirichlet: n interior unknowns; u_xx = L*u + boundaryLift*[gLeft;gRight].
%   Periodic: n distinct nodes without a repeated endpoint; boundaryLift is zero.
%   L approximates +d^2/dx^2. A Poisson solve for -u_xx therefore uses -L.
arguments
    n (1, 1) double {mustBeInteger, mustBePositive}
    bounds (1, 2) double {mustBeFinite}
    options.Boundary (1, 1) string {mustBeMember(options.Boundary, ["dirichlet", "periodic"])} = "dirichlet"
end
if bounds(2) <= bounds(1)
    error('mdemat:fdm:InvalidDomain', 'Domain endpoints must be strictly increasing.');
end
if options.Boundary == "periodic" && n < 3
    error('mdemat:fdm:PeriodicSize', 'A periodic centered stencil requires at least three nodes.');
end
if options.Boundary == "periodic"
    h = diff(bounds) / n;
    x = bounds(1) + (0:n-1)' * h;
else
    h = diff(bounds) / (n + 1);
    x = bounds(1) + (1:n)' * h;
end
e = ones(n, 1);
L = spdiags([e, -2*e, e], -1:1, n, n) / h^2;
boundaryLift = sparse(n, 2);
if options.Boundary == "periodic"
    L(1, n) = 1 / h^2;
    L(n, 1) = 1 / h^2;
else
    boundaryLift(1, 1) = 1 / h^2;
    boundaryLift(n, 2) = 1 / h^2;
end
end
