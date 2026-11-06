function [solution, report, history] = solveDiffusion(problem, space, integrator, options)
%SOLVEDIFFUSION Compose the original LDG elimination with generic ODE stages.
arguments
    problem (1, 1) struct
    space mdemat.common.fem.FiniteElementSpace
    integrator mdemat.common.ode.OdeIntegrator
    options.dt (1, 1) double {mustBeFinite, mustBePositive}
    options.penaltyType (2, 2) cell = {'boundary', 'right'; 'left', ''}
    options.storeHistory (1, 1) logical = false
    options.maxSteps (1, 1) double {mustBeNonnegative} = Inf
    options.checkpoint (1, 1) struct = struct()
end
expected = sort({'Type'; 'Diffusivity'; 'Initial'; 'BoundaryType'; 'BoundaryValue'});
if ~isequal(sort(fieldnames(problem)), expected) || ~isequal(problem.Type, "mdemat.diffusionProblem")
    error('mdemat:parabolic:InvalidProblem', 'Use diffusionProblem to construct model data.');
end
problem = mdemat.solver.parabolic.diffusionProblem(problem.Diffusivity, problem.Initial, ...
    boundaryType=problem.BoundaryType, boundaryValue=problem.BoundaryValue);
nd = space.NDims;
tensor = problem.Diffusivity;
if isscalar(tensor), tensor = tensor*eye(nd); end
if size(tensor, 1) ~= nd || size(tensor, 2) ~= nd
    error('mdemat:parabolic:DiffusivityShape', 'D must be a constant scalar or d-by-d tensor.');
end
for k = 1:size(tensor, 3)
    symmetric = (tensor(:, :, k)+tensor(:, :, k).')/2;
    if min(eig(symmetric)) <= 0
        error('mdemat:parabolic:InvalidDiffusivity', 'The symmetric part of D must be positive definite.');
    end
end
initial = space.fit(problem.Initial);
if size(initial, 2) ~= 1
    error('mdemat:parabolic:InitialShape', 'Supply one scalar initial field.');
end
assemblyClock = tic;
assembly = mdemat.common.fem.EllipticAssembly(space, char(problem.BoundaryType), options.penaltyType);
linear = assembly.assembleMatrix(tensor);
source = [];
if problem.BoundaryType == "dirichlet"
    source = @(~, t) assembly.assembleVector(tensor, ...
        @(x) boundary(problem.BoundaryValue, x, t));
end
assemblySeconds = toc(assemblyClock);
[coefficients, timeReport, history] = mdemat.common.ode.integrateLinear( ...
    initial, integrator, dt=options.dt, linear=linear, source=source, ...
    storeHistory=options.storeHistory, maxSteps=options.maxSteps, checkpoint=options.checkpoint);
weights = space.integrationWeights();
solution = struct('Type', "mdemat.diffusionSolution", 'U', coefficients, 'Time', timeReport.FinalTime);
report = struct('Time', timeReport, 'AssemblySeconds', assemblySeconds, ...
    'InitialMass', weights.'*initial, 'FinalMass', weights.'*coefficients, ...
    'MassChange', weights.'*(coefficients-initial), 'BoundaryType', problem.BoundaryType);
end

function values = boundary(callback, x, t)
if nargin(callback) == 1, values = callback(x); else, values = callback(x, t); end
end
