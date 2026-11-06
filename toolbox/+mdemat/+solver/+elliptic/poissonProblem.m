function problem = poissonProblem(force, options)
%POISSONPROBLEM Describe -Delta(u)=force without execution or assessment data.
%   P = poissonProblem(F, boundaryType="dirichlet", boundaryValue=G)
%   prescribes u=G(x). The default boundaryType="periodic" prescribes paired
%   periodic traces and zero physical mean. F is scalar numeric or a scalar
%   field callback; G is a scalar field callback, both using d-by-points x.
%   Discretization, exact references, plotting and files are caller concerns.
arguments
    force {mustBeA(force, {'numeric', 'function_handle'})} = 0
    options.boundaryType (1, 1) string ...
        {mustBeMember(options.boundaryType, ["periodic", "dirichlet"])} = "periodic"
    options.boundaryValue = []
end
if isnumeric(force) && (~isscalar(force) || ~isfinite(force))
    error('mdemat:elliptic:InvalidForce', 'A numeric force must be a finite scalar.');
end
if options.boundaryType == "dirichlet"
    if ~isa(options.boundaryValue, 'function_handle')
        error('mdemat:elliptic:InvalidBoundary', ...
            'Dirichlet data must be supplied as a scalar field callback.');
    end
elseif ~isempty(options.boundaryValue)
    error('mdemat:elliptic:InvalidBoundary', ...
        'Periodic problems do not accept prescribed Dirichlet data.');
end
problem = struct('Type', "mdemat.poissonProblem", 'Force', force, ...
    'BoundaryType', options.boundaryType, 'BoundaryValue', options.boundaryValue);
end
