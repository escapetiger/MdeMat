function problem = diffusionProblem(diffusivity, initial, options)
%DIFFUSIONPROBLEM Describe u_t=div(D*grad(u)) with explicit physical data.
arguments
    diffusivity {mustBeNumeric, mustBeReal, mustBeFinite, mustBeNonempty, mustBeMatrix}
    initial {mustBeA(initial, {'numeric', 'function_handle'})}
    options.boundaryType (1, 1) string ...
        {mustBeMember(options.boundaryType, ["periodic", "dirichlet"])} = "periodic"
    options.boundaryValue = []
end
if options.boundaryType == "dirichlet" && ~isa(options.boundaryValue, 'function_handle')
    error('mdemat:parabolic:InvalidBoundary', 'Prescribe Dirichlet data through boundaryValue(x,t).');
elseif options.boundaryType == "periodic" && ~isempty(options.boundaryValue)
    error('mdemat:parabolic:InvalidBoundary', 'Periodic diffusion has no prescribed boundary values.');
end
problem = struct('Type', "mdemat.diffusionProblem", 'Diffusivity', diffusivity, ...
    'Initial', initial, 'BoundaryType', options.boundaryType, ...
    'BoundaryValue', options.boundaryValue);
end
