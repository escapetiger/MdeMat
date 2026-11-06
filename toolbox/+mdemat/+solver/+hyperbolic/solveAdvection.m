function [solution, report, history] = solveAdvection(problem, space, integrator, options)
%SOLVEADVECTION Compose DG transport kernels with an explicit ODE workspace.
%   method="uwdg" uses conservative upwind assembly. method="sldg" uses the
%   existing semi-Lagrangian overlap/transfer kernel and requires its SL space.
%   No execution controller, exact reference, observer or file is constructed.
arguments
    problem (1, 1) struct
    space mdemat.common.fem.FiniteElementSpace
    integrator mdemat.common.ode.OdeIntegrator
    options.dt (1, 1) double {mustBeFinite, mustBePositive}
    options.method (1, 1) string {mustBeMember(options.method, ["uwdg", "sldg"])} = "uwdg"
    options.storeHistory (1, 1) logical = false
    options.maxSteps (1, 1) double {mustBeNonnegative} = Inf
    options.checkpoint (1, 1) struct = struct()
end
expected = sort({'Type'; 'Velocity'; 'Initial'; 'BoundaryType'; 'BoundaryValue'});
if ~isequal(sort(fieldnames(problem)), expected) || ~isequal(problem.Type, "mdemat.advectionProblem")
    error('mdemat:hyperbolic:InvalidProblem', 'Use advectionProblem to construct model data.');
end
problem = mdemat.solver.hyperbolic.advectionProblem(problem.Velocity, problem.Initial, ...
    boundaryType=problem.BoundaryType, boundaryValue=problem.BoundaryValue);
velocity = problem.Velocity;
if isvector(velocity), velocity = velocity(:); end
if size(velocity, 1) ~= space.NDims || size(velocity, 2) ~= 1
    error('mdemat:hyperbolic:VelocityShape', 'Velocity must have one entry per spatial dimension.');
end
initial = space.fit(problem.Initial);
if size(initial, 2) ~= 1
    error('mdemat:hyperbolic:InitialShape', 'Supply one scalar initial field.');
end
assemblyClock = tic;
linear = [];
source = [];
mass = [];
if options.method == "uwdg"
    assembly = mdemat.common.fem.UpwindAssembly(space, char(problem.BoundaryType));
    linear = -assembly.assembleMatrix(velocity);
    if problem.BoundaryType == "dirichlet"
        source = @(~, t) -assembly.assembleVector(velocity, velocity, ...
            @(x) boundary(problem.BoundaryValue, x, t));
    end
else
    assembly = mdemat.common.fem.SemiLagrangianAssembly(space, char(problem.BoundaryType));
    mass = @transfer;
end
assemblySeconds = toc(assemblyClock);
[coefficients, timeReport, history] = mdemat.common.ode.integrateLinear( ...
    initial, integrator, dt=options.dt, mass=mass, linear=linear, source=source, ...
    storeHistory=options.storeHistory, maxSteps=options.maxSteps, checkpoint=options.checkpoint);
weights = space.integrationWeights();
solution = struct('Type', "mdemat.advectionSolution", 'U', coefficients, ...
    'Time', timeReport.FinalTime, 'Method', options.method);
report = struct('Time', timeReport, 'AssemblySeconds', assemblySeconds, ...
    'InitialMass', weights.'*initial, 'FinalMass', weights.'*coefficients, ...
    'MassChange', weights.'*(coefficients-initial), 'BoundaryType', problem.BoundaryType);

    function [lift, matrix] = transfer(tBegin, tEnd)
        assembly.reset();
        if nargout > 1
            matrix = assembly.assembleMatrix(velocity, tBegin, tEnd);
        end
        lift = [];
        if problem.BoundaryType == "dirichlet" && tBegin ~= tEnd
            lift = assembly.assembleVector(velocity, tBegin, tEnd, ...
                @(x, varargin) boundaryWithTime(problem.BoundaryValue, x, varargin{:}));
        end
    end
end

function values = boundary(callback, x, t)
if nargin(callback) == 1, values = callback(x); else, values = callback(x, t); end
end

function values = boundaryWithTime(callback, x, varargin)
% Legacy SL boundary evaluation may provide time positionally or as 't',t.
if isempty(varargin), t = 0; else, t = varargin{end}; end
values = boundary(callback, x, t);
end
