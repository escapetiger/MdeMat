classdef TestEvolutionInterfaces < matlab.unittest.TestCase
    %TESTEVOLUTIONINTERFACES Verify independent evolution and restart contracts.
    methods (Test)
        function completedStepRestartMatchesUninterruptedTrajectories(testCase)
            constructors = {@mdemat.common.ode.BeIntegrator, ...
                @mdemat.common.ode.Bdf2Integrator, @mdemat.common.ode.Bdf3Integrator, ...
                @mdemat.common.ode.Sdirk2Integrator, @mdemat.common.ode.Exrk4Integrator};
            mass = diag([1, 2]);
            linear = [-1, 0.2; 0, -4];
            for i = 1:numel(constructors)
                integrator = constructors{i}(0.5);
                [whole, ~, fullHistory] = mdemat.common.ode.integrateLinear([1; 2], ...
                    integrator, dt=0.07, mass=mass, linear=linear, storeHistory=true);
                [~, partial, prefix] = mdemat.common.ode.integrateLinear([1; 2], ...
                    constructors{i}(0.5), dt=0.07, mass=mass, linear=linear, ...
                    maxSteps=3, storeHistory=true);
                [resumed, report, suffix] = mdemat.common.ode.integrateLinear([1; 2], ...
                    constructors{i}(0.5), dt=0.07, mass=mass, linear=linear, ...
                    checkpoint=partial.Checkpoint, storeHistory=true);
                testCase.verifyFalse(partial.Completed);
                testCase.verifyTrue(report.Completed);
                testCase.verifyTrue(report.Resumed);
                testCase.verifyEqual(resumed, whole, 'AbsTol', 2e-13);
                testCase.verifyEqual([prefix.Time; suffix.Time(2:end)], fullHistory.Time);
                testCase.verifyEqual([prefix.Coefficients, suffix.Coefficients(:, 2:end)], ...
                    fullHistory.Coefficients, 'AbsTol', 2e-13);
            end
        end

        function rejectsIntermediateStageRestart(testCase)
            [~, report] = mdemat.common.ode.integrateLinear(1, ...
                mdemat.common.ode.BeIntegrator(1), dt=0.1, linear=-1, maxSteps=1);
            checkpoint = report.Checkpoint;
            checkpoint.Phase = "inside-stage";
            testCase.verifyError(@() mdemat.common.ode.integrateLinear(1, ...
                mdemat.common.ode.BeIntegrator(1), dt=0.1, linear=-1, ...
                checkpoint=checkpoint), 'mdemat:ode:InvalidCheckpoint');
        end

        function upwindPeriodicTransportConvergesAndConservesMass(testCase)
            errors = zeros(3, 1);
            for level = 1:3
                space = TestEvolutionInterfaces.space(4*2^(level-1), false);
                initial = @(x) 1+0.1*sin(2*pi*x(1, :));
                exact = @(x) 1+0.1*sin(2*pi*(x(1, :)-0.1));
                problem = mdemat.solver.hyperbolic.advectionProblem(1, initial);
                [solution, report] = mdemat.solver.hyperbolic.solveAdvection( ...
                    problem, space, mdemat.common.ode.Exrk4Integrator(0.1), dt=0.001);
                errors(level) = TestEvolutionInterfaces.error(space, solution.U, exact);
                testCase.verifyEqual(report.MassChange, 0, 'AbsTol', 2e-12);
            end
            rate = log(errors(2)/errors(3))/log(2);
            testCase.verifyGreaterThan(rate, 2.5);
        end

        function semiLagrangianUsesTheExistingTransferKernel(testCase)
            space = TestEvolutionInterfaces.space(4, true);
            problem = mdemat.solver.hyperbolic.advectionProblem(1, 1.7);
            [solution, report] = mdemat.solver.hyperbolic.solveAdvection( ...
                problem, space, mdemat.common.ode.BeIntegrator(0.3), ...
                dt=0.15, method="sldg");
            testCase.verifyEqual(space.eval(space.Element.Volume.Nodes, solution.U), ...
                1.7*ones(space.NGlobalDofs, 1), 'AbsTol', 2e-12);
            testCase.verifyEqual(report.MassChange, 0, 'AbsTol', 2e-12);
        end

        function periodicDiffusionMatchesItsDecayAndMass(testCase)
            space = TestEvolutionInterfaces.space(12, false);
            initial = @(x) 1+0.1*sin(2*pi*x(1, :));
            exact = @(x) 1+0.1*exp(-4*pi^2*0.005)*sin(2*pi*x(1, :));
            problem = mdemat.solver.parabolic.diffusionProblem(1, initial);
            [solution, report] = mdemat.solver.parabolic.solveDiffusion( ...
                problem, space, mdemat.common.ode.Sdirk2Integrator(0.005), dt=0.0001);
            error = TestEvolutionInterfaces.error(space, solution.U, exact);
            testCase.verifyLessThan(error, 2e-4);
            testCase.verifyEqual(report.MassChange, 0, 'AbsTol', 2e-12);
        end

        function semiLagrangianLongStepsUseBothInflowDirections(testCase)
            % The whole initial footprint leaves the single cell. No
            % interior transfer block may be assembled from an empty set.
            for velocity = [-1, 1]
                space = TestEvolutionInterfaces.space(1, true);
                problem = mdemat.solver.hyperbolic.advectionProblem(velocity, 2, ...
                    boundaryType="dirichlet", boundaryValue=@(x, t) 1+0*x(1, :)+0*t);
                solution = mdemat.solver.hyperbolic.solveAdvection(problem, space, ...
                    mdemat.common.ode.BeIntegrator(2), dt=2, method="sldg");
                testCase.verifyEqual(space.eval(space.Element.Volume.Nodes, solution.U), ...
                    ones(space.NGlobalDofs, 1), 'AbsTol', 2e-12);
            end
        end

        function semiLagrangianUsesTheLastPhysicalEntryAtACorner(testCase)
            mesh = mdemat.data.mesh.UniformGrid([0, 2, -1, 3], [1, 1]);
            element = mdemat.common.fem.BH1SLOrthotopeElement.nodal(2, 3);
            space = mdemat.common.fem.FiniteElementSpace(mesh, element);
            problem = mdemat.solver.hyperbolic.advectionProblem([2; -4], 0, ...
                boundaryType="dirichlet", boundaryValue=@(x, t) t+0*x(1, :));
            solution = mdemat.solver.hyperbolic.solveAdvection(problem, space, ...
                mdemat.common.ode.BeIntegrator(1), dt=1, method="sldg");
            % Backward characteristics enter through the first face reached
            % from the destination: t_entry=1-min(x/2,(3-y)/4). The unequal
            % bounds and mixed velocity signs expose axis-broadcast errors.
            % Nodal projection
            % uses these same volume nodes; no claim of polynomial exactness
            % across the characteristic fan is needed here.
            nodes = space.Element.Volume.Nodes;
            exact = space.feval(@(x) 1-min(x(1, :)/2, (3-x(2, :))/4), nodes);
            testCase.verifyEqual(space.eval(nodes, solution.U), exact, 'AbsTol', 3e-12);
        end

        function advectionInterfaceResumesTheSameCoefficientTrajectory(testCase)
            space = TestEvolutionInterfaces.space(4, false);
            problem = mdemat.solver.hyperbolic.advectionProblem(1, ...
                @(x) 1+0.1*sin(2*pi*x(1, :)));
            whole = mdemat.solver.hyperbolic.solveAdvection(problem, space, ...
                mdemat.common.ode.Exrk4Integrator(0.1), dt=0.01);
            [~, partial] = mdemat.solver.hyperbolic.solveAdvection(problem, space, ...
                mdemat.common.ode.Exrk4Integrator(0.1), dt=0.01, maxSteps=3);
            resumed = mdemat.solver.hyperbolic.solveAdvection(problem, space, ...
                mdemat.common.ode.Exrk4Integrator(0.1), dt=0.01, ...
                checkpoint=partial.Time.Checkpoint);
            testCase.verifyEqual(resumed.U, whole.U, 'AbsTol', 2e-13);
            testCase.verifyEqual(resumed.Time, 0.1, 'AbsTol', eps);
        end
    end

    methods (Static, Access = private)
        function space = space(cells, semiLagrangian)
            mesh = mdemat.data.mesh.UniformGrid([0, 1], cells);
            if semiLagrangian
                element = mdemat.common.fem.BH1SLOrthotopeElement.nodal(1, 3);
            else
                element = mdemat.common.fem.BH1OrthotopeElement.nodal(1, 3);
            end
            space = mdemat.common.fem.FiniteElementSpace(mesh, element);
        end

        function error = error(space, coefficients, exact)
            quadrature = mdemat.common.quadrature.GaussLegendreRule(1);
            [nodes, weights] = quadrature.generate(space.Element.Geometry, 8);
            values = space.eval(nodes, coefficients)-space.feval(exact, nodes);
            determinants = space.Mesh.computeElementJacobianDeterminants();
            if isscalar(determinants)
                determinants = repmat(determinants, space.NMeshElements, 1);
            end
            physical = kron(determinants(:), weights(:));
            error = sqrt(physical.'*(values.^2));
        end
    end
end
