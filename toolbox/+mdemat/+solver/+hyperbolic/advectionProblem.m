function problem = advectionProblem(velocity, initial, options)
%ADVECTIONPROBLEM Describe u_t+div(velocity*u)=0 independently of execution.
%   Velocity is a constant vector, with one entry per dimension.
%   Initial(x) uses components-by-points samples; boundaryValue(x,t) supplies
%   prescribed inflow data when boundaryType="dirichlet".
arguments
    velocity {mustBeNumeric, mustBeReal, mustBeFinite, mustBeNonempty, mustBeVector}
    initial {mustBeA(initial, {'numeric', 'function_handle'})}
    options.boundaryType (1, 1) string ...
        {mustBeMember(options.boundaryType, ["periodic", "dirichlet"])} = "periodic"
    options.boundaryValue = []
end
if options.boundaryType == "dirichlet" && ~isa(options.boundaryValue, 'function_handle')
    error('mdemat:hyperbolic:InvalidBoundary', 'Prescribe inflow through boundaryValue(x,t).');
elseif options.boundaryType == "periodic" && ~isempty(options.boundaryValue)
    error('mdemat:hyperbolic:InvalidBoundary', 'Periodic transport has no prescribed boundary values.');
end
problem = struct('Type', "mdemat.advectionProblem", 'Velocity', velocity, ...
    'Initial', initial, 'BoundaryType', options.boundaryType, ...
    'BoundaryValue', options.boundaryValue);
end
