function [solution, report] = solvePoisson(problem, space, options)
%SOLVEPOISSON Solve a scalar Poisson problem in an explicit DG space.
%   [S,R] = solvePoisson(P,V,method="ldg") assembles and solves the mixed
%   system from docs/mathematics/elliptic.md (P2-P6). method="pdg" eliminates
%   the gradient before solving and then reconstructs every gradient column.
%   S.U is Nx-by-1 and S.Q is Nx-by-d, with Q=+grad(U). S/R are plain records
%   and retain no space, state, visualization, exact-reference or file handles.
%   Sample coefficients explicitly with V.eval(xRef,S.U) and V.eval(xRef,S.Q).
arguments
    problem (1, 1) struct
    space mdemat.common.fem.FiniteElementSpace
    options.method (1, 1) string {mustBeMember(options.method, ["ldg", "pdg"])} = "ldg"
    options.penaltyType (2, 2) cell = {'boundary', 'right'; 'left', ''}
    options.residualTolerance (1, 1) double {mustBeFinite, mustBePositive} = 1e-10
    options.compatibilityTolerance (1, 1) double {mustBeFinite, mustBeNonnegative} = 1e-10
end
validateProblem(problem);
assemblyClock = tic;
nd = space.NDims;
ng = space.NGlobalDofs;
assembly = mdemat.common.fem.AdjointAssembly(space, char(problem.BoundaryType), ...
    options.penaltyType);
source = mdemat.common.fem.SourceAssembly(space);
dual = assembly.assembleMatrix(eye(nd));
force = source.assembleVector(problem.Force);
if ~isequal(size(force), [ng, 1]) || any(~isfinite(force(:)))
    error('mdemat:elliptic:InvalidForce', 'The force must evaluate to one finite scalar field.');
end
boundaryPrimary = zeros(ng, 1);
boundaryGradient = zeros(ng, nd);
if problem.BoundaryType == "dirichlet"
    trace = assembly.assembleVector(eye(nd), problem.BoundaryValue);
    for d = 1:nd
        boundaryPrimary = boundaryPrimary + trace{1, 1}{d};
        boundaryGradient(:, d) = trace{2, 1}{d};
    end
end
physicalWeights = space.integrationWeights();
volume = sum(physicalWeights.'*source.assembleVector(1));
compatibilityDefect = physicalWeights.'*force;
periodic = problem.BoundaryType == "periodic";
compatibilityScale = max(1, norm(physicalWeights, 1)*norm(force, inf));
if periodic && abs(compatibilityDefect) > options.compatibilityTolerance*compatibilityScale
    error('mdemat:elliptic:IncompatiblePeriodicForce', ...
        'The discrete physical forcing integral is %g; periodic Poisson requires zero.', ...
        compatibilityDefect);
end
if options.method == "ldg"
    blocks = repmat({sparse(ng, ng)}, nd+1, nd+1);
    for d = 1:nd
        blocks{1, 1} = blocks{1, 1} - dual{1, 1}{d};
        blocks{1, d+1} = -dual{1, 2}{d};
        blocks{d+1, 1} = -dual{2, 1}{d};
        blocks{d+1, d+1} = speye(ng);
    end
    matrix = mdemat.common.linalg.block(blocks);
    rhs = [force+boundaryPrimary; boundaryGradient(:)];
    constraint = [physicalWeights; zeros(nd*ng, 1)];
    constant = [source.assembleVector(1); zeros(nd*ng, 1)];
else
    matrix = sparse(ng, ng);
    rhs = force+boundaryPrimary;
    for d = 1:nd
        matrix = matrix - dual{1, 1}{d} - dual{1, 2}{d}*dual{2, 1}{d};
        rhs = rhs + dual{1, 2}{d}*boundaryGradient(:, d);
    end
    constraint = physicalWeights;
    constant = source.assembleVector(1);
end
if periodic
    matrix = [matrix, constant; constraint.', 0];
    rhs = [rhs; 0];
end
assemblyTime = toc(assemblyClock);
linear = mdemat.common.linalg.LinearSolver(nnzTh=Inf, condTh=Inf, tol=options.residualTolerance);
solveClock = tic;
[packed, statistics] = linear.solve(A=matrix, b=rhs);
solveTime = toc(solveClock);
multiplier = 0;
if periodic
    multiplier = packed(end);
    packed = packed(1:end-1);
end
if options.method == "ldg"
    packed = reshape(packed, ng, nd+1);
    potential = packed(:, 1);
    gradient = packed(:, 2:end);
else
    potential = packed;
    gradient = zeros(ng, nd);
    for d = 1:nd
        gradient(:, d) = dual{2, 1}{d}*potential+boundaryGradient(:, d);
    end
end
meanValue = (physicalWeights.'*potential)/volume;
solution = struct('Type', "mdemat.poissonSolution", 'U', potential, ...
    'Q', gradient, 'Method', options.method, 'BoundaryType', problem.BoundaryType);
report = struct('Converged', statistics.converged, ...
    'RelativeResidual', statistics.relativeResidual, 'Method', options.method, ...
    'SystemSize', size(matrix, 1), 'MatrixNonzeros', nnz(matrix), ...
    'SpatialDofs', ng, 'Dimensions', nd, 'PhysicalMean', meanValue, ...
    'PeriodicMultiplier', multiplier, 'ForcingIntegral', compatibilityDefect, ...
    'AssemblySeconds', assemblyTime, 'SolveSeconds', solveTime);
end

function validateProblem(problem)
expected = sort({'Type'; 'Force'; 'BoundaryType'; 'BoundaryValue'});
if ~isequal(sort(fieldnames(problem)), expected) || ...
        ~isequal(problem.Type, "mdemat.poissonProblem")
    error('mdemat:elliptic:InvalidProblem', ...
        'Use poissonProblem to create a problem with the supported schema.');
end
% Revalidate caller-modified records, including boundary/force pairing.
mdemat.solver.elliptic.poissonProblem(problem.Force, ...
    boundaryType=problem.BoundaryType, boundaryValue=problem.BoundaryValue);
end
