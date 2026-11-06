classdef TestPoissonInterface < matlab.unittest.TestCase
    %TESTPOISSONINTERFACE Test explicit mathematical inputs and sampled fields.
    methods (Test)
        function quadraticDirichletFieldsAreExactInOneAndTwoDimensions(testCase)
            for nd = 1:2
                space = TestPoissonInterface.space(nd, 2, 3);
                if nd == 1
                    exact = @(x) x(1, :).^2+0.4*x(1, :)+0.2;
                    derivative = @(x) 2*x(1, :)+0.4;
                    force = -2;
                else
                    exact = @(x) x(1, :).^2+2*x(2, :).^2 ...
                        +0.5*x(1, :)-0.75*x(2, :)+1;
                    derivative = @(x) [2*x(1, :)+0.5; 4*x(2, :)-0.75];
                    force = -6;
                end
                problem = mdemat.solver.elliptic.poissonProblem(force, ...
                    boundaryType="dirichlet", boundaryValue=exact);
                for method = ["ldg", "pdg"]
                    [solution, report] = mdemat.solver.elliptic.solvePoisson( ...
                        problem, space, method=method);
                    nodes = space.Element.Volume.Nodes;
                    testCase.verifyEqual(space.eval(nodes, solution.U), ...
                        space.feval(exact, nodes), 'AbsTol', 2e-10);
                    testCase.verifyEqual(space.eval(nodes, solution.Q), ...
                        space.feval(derivative, nodes), 'AbsTol', 2e-10);
                    testCase.verifySize(solution.Q, [space.NGlobalDofs, nd]);
                    testCase.verifyTrue(report.Converged);
                    testCase.verifyLessThan(report.RelativeResidual, 1e-10);
                end
            end
        end

        function nonuniformPeriodicSolutionHasThePhysicalGauge(testCase)
            mesh = mdemat.data.mesh.Grid({[0, 0.08, 0.31, 0.72, 1]});
            element = mdemat.common.fem.BH1OrthotopeElement.nodal(1, 5);
            space = mdemat.common.fem.FiniteElementSpace(mesh, element);
            exact = @(x) x(1, :).^2.*(1-x(1, :)).^2-1/30;
            gradient = @(x) 2*x(1, :)-6*x(1, :).^2+4*x(1, :).^3;
            problem = mdemat.solver.elliptic.poissonProblem( ...
                @(x) -2+12*x(1, :)-12*x(1, :).^2);
            for method = ["ldg", "pdg"]
                [solution, report] = mdemat.solver.elliptic.solvePoisson( ...
                    problem, space, method=method);
                nodes = element.Volume.Nodes;
                testCase.verifyEqual(space.eval(nodes, solution.U), ...
                    space.feval(exact, nodes), 'AbsTol', 2e-10);
                testCase.verifyEqual(space.eval(nodes, solution.Q), ...
                    space.feval(gradient, nodes), 'AbsTol', 2e-10);
                testCase.verifyEqual(report.PhysicalMean, 0, 'AbsTol', 1e-12);
                testCase.verifyEqual(report.PeriodicMultiplier, 0, 'AbsTol', 1e-10);
            end
        end

        function rejectsIncompatiblePeriodicDataAndUnknownProblemFields(testCase)
            space = TestPoissonInterface.space(1, 3, 3);
            problem = mdemat.solver.elliptic.poissonProblem(1);
            testCase.verifyError(@() mdemat.solver.elliptic.solvePoisson(problem, space), ...
                'mdemat:elliptic:IncompatiblePeriodicForce');
            problem = mdemat.solver.elliptic.poissonProblem(0);
            problem.Reference = @(x) 0;
            testCase.verifyError(@() mdemat.solver.elliptic.solvePoisson(problem, space), ...
                'mdemat:elliptic:InvalidProblem');
        end

        function smoothPeriodicFieldsConvergeInOneAndTwoDimensions(testCase)
            for nd = 1:2
                errors = zeros(3, 2);
                refinements = [3, 6, 12];
                if nd == 1
                    exact = @(x) sin(2*pi*x(1, :));
                    derivative = @(x) 2*pi*cos(2*pi*x(1, :));
                    force = @(x) 4*pi^2*sin(2*pi*x(1, :));
                    method = "ldg";
                else
                    exact = @(x) sin(2*pi*x(1, :)).*sin(2*pi*x(2, :));
                    derivative = @(x) 2*pi*[cos(2*pi*x(1, :)).*sin(2*pi*x(2, :)); ...
                        sin(2*pi*x(1, :)).*cos(2*pi*x(2, :))];
                    force = @(x) 8*pi^2*sin(2*pi*x(1, :)).*sin(2*pi*x(2, :));
                    method = "pdg";
                end
                problem = mdemat.solver.elliptic.poissonProblem(force);
                for level = 1:3
                    space = TestPoissonInterface.space(nd, refinements(level), 3);
                    solution = mdemat.solver.elliptic.solvePoisson(problem, space, method=method);
                    rule = mdemat.common.quadrature.GaussLegendreRule(nd);
                    [nodes, weights] = rule.generate(space.Element.Geometry, 8*ones(1, nd));
                    determinants = space.Mesh.computeElementJacobianDeterminants();
                    if isscalar(determinants)
                        determinants = repmat(determinants, space.NMeshElements, 1);
                    end
                    physicalWeights = kron(determinants(:), weights(:));
                    potentialError = space.eval(nodes, solution.U)-space.feval(exact, nodes);
                    gradientError = space.eval(nodes, solution.Q)-space.feval(derivative, nodes);
                    errors(level, 1) = sqrt(physicalWeights.'*(potentialError.^2));
                    errors(level, 2) = sqrt(sum(physicalWeights.'*(gradientError.^2)));
                end
                rates = log(errors(2, :)./errors(3, :))/log(2);
                testCase.verifyLessThan(errors(3, :), errors(2, :));
                testCase.verifyGreaterThan(rates(1), 2.5);
                testCase.verifyGreaterThan(rates(2), 1.8);
            end
        end
    end

    methods (Static, Access = private)
        function result = space(nd, cells, order)
            mesh = mdemat.data.mesh.UniformGrid(repmat([0, 1], 1, nd), cells*ones(1, nd));
            element = mdemat.common.fem.BH1OrthotopeElement.nodal(nd, order);
            result = mdemat.common.fem.FiniteElementSpace(mesh, element);
        end
    end
end
