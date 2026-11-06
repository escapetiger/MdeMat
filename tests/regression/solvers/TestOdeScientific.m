classdef TestOdeScientific < matlab.unittest.TestCase
    %TESTODESCIENTIFIC Independent order, mass, additive and nonlinear checks.
    methods (Test)
        function everyNamedMethodConvergesOnNonautonomousLinearData(testCase)
            [names, orders] = TestOdeScientific.catalog();
            exact = exp(-1)+1;
            for method = 1:numel(names)
                errors = zeros(1, 3);
                for level = 1:3
                    integrator = TestOdeScientific.create(names(method), 1);
                    value = TestOdeScientific.advance(integrator, 1, 8*2^(level-1), ...
                        L=-1, S=@(~, t) 2*t+t.^2, M=1);
                    errors(level) = abs(value-exact);
                end
                rate = log(errors(2)/errors(3))/log(2);
                testCase.verifyGreaterThan(rate, orders(method)-0.35, char(names(method)));
            end
        end

        function timeDependentMassPreservesTheConservedProduct(testCase)
            names = TestOdeScientific.catalog();
            for method = 1:numel(names)
                integrator = TestOdeScientific.create(names(method), 0.7);
                value = TestOdeScientific.advance(integrator, 1, 14, M=@(t) 1+t);
                testCase.verifyEqual(value, 1/1.7, 'AbsTol', 3e-11);
            end
        end

        function implicitStagesSolveTheActualNonlinearEquation(testCase)
            integrator = mdemat.common.ode.BeIntegrator(0.2);
            result = TestOdeScientific.advance(integrator, 1, 1, F=@(u) -u.^2, M=1);
            expected = (sqrt(1.8)-1)/0.4;
            testCase.verifyEqual(result, expected, 'AbsTol', 2e-11);
            testCase.verifyTrue(integrator.SolveReports{end}.Converged);
            testCase.verifyLessThanOrEqual(integrator.SolveReports{end}.ResidualNorm, ...
                integrator.SolveReports{end}.Tolerance);
            report = integrator.SolveReports{end};
            testCase.verifyNumElements(report.NonlinearStatistics.krylovStats.exitStatusId, report.Iterations);
            testCase.verifyNumElements(report.NonlinearStatistics.krylovStats.rawIterations, report.Iterations);
            testCase.verifyEqual(report.NonlinearStatistics.fNorm(end), report.ResidualNorm, 'AbsTol', 1e-13);

            % An explicit predictor is negative for this stiff positive
            % decay and would select the other algebraic quadratic root.
            stiff = mdemat.common.ode.BeIntegrator(10);
            result = TestOdeScientific.advance(stiff, 1, 1, F=@(u) -u.^2, M=1);
            testCase.verifyEqual(result, (sqrt(41)-1)/20, 'AbsTol', 2e-11);

            names = ["ImpIntegrator", "Sdirk2Integrator", "Sdirk3Integrator", ...
                "Sdirk4Integrator", "Esdirk3Integrator"];
            orders = [2, 2, 3, 4, 3];
            for method = 1:numel(names)
                errors = zeros(1, 3);
                for level = 1:3
                    integrator = TestOdeScientific.create(names(method), 0.5);
                    result = TestOdeScientific.advance(integrator, 1, 4*2^(level-1), ...
                        F=@(u) -u.^2, M=1);
                    errors(level) = abs(result-2/3);
                end
                testCase.verifyGreaterThan(log(errors(2)/errors(3))/log(2), ...
                    orders(method)-0.4, char(names(method)));
            end
        end

        function nonlinearFailureIsExplicit(testCase)
            integrator = mdemat.common.ode.BeIntegrator(1);
            integrator.NonlinearMaxIterations = 1;
            testCase.verifyError(@() TestOdeScientific.advance(integrator, 1, 1, ...
                F=@(u) -10*u.^3, M=1), 'mdemat:ode:NonlinearStageFailed');
        end

        function additiveImexTreatsItsNonlinearTermExplicitlyAtEveryStage(testCase)
            names = ["Ars111Integrator", "Ars222Integrator", "Ars443Integrator"];
            orders = [1, 2, 3];
            for method = 1:numel(names)
                errors = zeros(1, 3);
                for level = 1:3
                    integrator = TestOdeScientific.create(names(method), 1);
                    value = TestOdeScientific.advance(integrator, 1, 16*2^(level-1), ...
                        L=-5, F=@(u, t) 2*exp(t)*u.^2, ...
                        S=@(~, t) 2*exp(-t), M=1);
                    errors(level) = abs(value-exp(-1));
                end
                testCase.verifyGreaterThan(log(errors(2)/errors(3))/log(2), ...
                    orders(method)-0.4, char(names(method)));
            end
        end

        function noncommutingSplitsHaveTheirDeclaredOrders(testCase)
            first = [0, 1; -1, 0];
            second = [-1, 0; 0, -2];
            initial = [1; 0.2];
            exact = expm(0.5*(first+second))*initial;
            constructors = {@mdemat.common.ode.LieTrotterIntegrator, @mdemat.common.ode.StrangIntegrator};
            for method = 1:2
                errors = zeros(1, 3);
                for level = 1:3
                    integrator = constructors{method}(factors={ ...
                        mdemat.common.ode.Sdirk4Integrator(0.5), mdemat.common.ode.Sdirk4Integrator(0.5)});
                    value = TestOdeScientific.advance(integrator, initial, 4*2^(level-1), ...
                        L={first, second}, M={eye(2), eye(2)});
                    errors(level) = norm(value-exact);
                end
                testCase.verifyGreaterThan(log(errors(2)/errors(3))/log(2), method-0.25);
            end
        end

        function augmentedStagesSatisfyTheIndependentAlgebraicConstraint(testCase)
            names = ["BeIntegrator", "Sdirk2Integrator", "Ars222Integrator", "Bdf2Integrator"];
            for method = 1:numel(names)
                integrator = TestOdeScientific.create(names(method), 0.2);
                value = TestOdeScientific.advance(integrator, [1; 1], 20, ...
                    M=1, L=-1, P=1, A=1, B=-1, C=0, S=@(~, t) exp(t));
                testCase.verifyEqual(value(1)-value(2), 0, 'AbsTol', 1e-11);
                testCase.verifyEqual(value(1), exp(0.2), 'AbsTol', 2e-3);
            end
        end

        function omittedMassIsIdentityOnDifferentialAxesOnly(testCase)
            names = ["BeIntegrator", "Sdirk2Integrator", "Exrk4Integrator", ...
                "Ars222Integrator", "Bdf2Integrator"];
            for method = 1:numel(names)
                value = TestOdeScientific.advance(TestOdeScientific.create(names(method), 0.2), ...
                    [1; 1], 20, L=-1, P=1, A=1, B=-1, C=0, S=@(~, t) exp(t));
                testCase.verifyEqual(value(1)-value(2), 0, 'AbsTol', 1e-11);
                testCase.verifyEqual(value(1), exp(0.2), 'AbsTol', 2e-3);
            end
        end

        function splitImplicitFactorsRetainTheirNonlinearTerms(testCase)
            integrator = mdemat.common.ode.StrangIntegrator(factors={ ...
                mdemat.common.ode.BeIntegrator(0.2), mdemat.common.ode.Sdirk2Integrator(0.2)});
            value = TestOdeScientific.advance(integrator, 1, 1, ...
                F={@(u) ones(size(u)), @(u) 2*ones(size(u))});
            testCase.verifyEqual(value, 1.6, 'AbsTol', 2e-11);
        end

        function semiLagrangianInflowUsesBothTimeEndpoints(testCase)
            mesh = mdemat.data.mesh.UniformGrid([0, 1], 4);
            element = mdemat.common.fem.BH1SLOrthotopeElement.nodal(1, 3);
            space = mdemat.common.fem.FiniteElementSpace(mesh, element);
            problem = mdemat.solver.hyperbolic.advectionProblem(1, @(x) exp(-x(1, :)), ...
                boundaryType="dirichlet", boundaryValue=@(x, t) exp(t-x(1, :)));
            solution = mdemat.solver.hyperbolic.solveAdvection(problem, space, ...
                mdemat.common.ode.Bdf2Integrator(0.2), dt=0.1, method="sldg");
            rule = mdemat.common.quadrature.GaussLegendreRule(1);
            [nodes, weights] = rule.generate(element.Geometry, 8);
            error = space.eval(nodes, solution.U)-space.feval(@(x) exp(0.2-x(1, :)), nodes);
            testCase.verifyLessThan(sqrt(kron(ones(4, 1)/4, weights(:)).'*(error.^2)), 4e-4);
        end

        function allUnsplitMethodsRespectSemiLagrangianStageAndHistoryTransfers(testCase)
            [names, ~] = TestOdeScientific.catalog();
            for method = 1:numel(names)
                for velocity = [-0.9, 0.7]
                    mesh = mdemat.data.mesh.UniformGrid([-0.3, 1.2], 4);
                    element = mdemat.common.fem.BH1SLOrthotopeElement.nodal(1, 3);
                    space = mdemat.common.fem.FiniteElementSpace(mesh, element);
                    problem = mdemat.solver.hyperbolic.advectionProblem(velocity, ...
                        @(x) 1+2*x(1, :), boundaryType="dirichlet", ...
                        boundaryValue=@(x, t) 1+2*(x(1, :)-velocity*t));
                    solution = mdemat.solver.hyperbolic.solveAdvection(problem, space, ...
                        TestOdeScientific.create(names(method), 0.27), dt=0.09, method="sldg");
                    nodes = space.Element.Volume.Nodes;
                    expected = space.feval(@(x) 1+2*(x(1, :)-0.27*velocity), nodes);
                    testCase.verifyEqual(space.eval(nodes, solution.U), expected, ...
                        char(names(method)), 'AbsTol', 2e-11);
                end
            end
        end

        function varyingMassAndTimeDependentLinearOperatorRetainOrder(testCase)
            [names, orders] = TestOdeScientific.catalog();
            exact = (exp(-1)+1)/2;
            for method = 1:numel(names)
                errors = zeros(1, 3);
                for level = 1:3
                    integrator = TestOdeScientific.create(names(method), 1);
                    result = TestOdeScientific.advance(integrator, 1, 8*2^(level-1), ...
                        M=@(t) 1+t, L=@(~, t) -(1+t), S=@(~, t) 2*t+t.^2);
                    errors(level) = abs(result-exact);
                end
                testCase.verifyGreaterThan(log(errors(2)/errors(3))/log(2), ...
                    orders(method)-0.35, char(names(method)));
            end
        end

        function variableStepBdfUsesItsActualHistoryAndOrderMatchedStart(testCase)
            names = ["Bdf2Integrator", "Bdf3Integrator"];
            for method = 1:2
                errors = zeros(1, 3);
                for level = 1:3
                    integrator = TestOdeScientific.create(names(method), 1);
                    integrator.reset();
                    value = 1;
                    base = 1/(8*2^(level-1));
                    while ~integrator.Timeline.IsExhausted
                        integrator.setTimeStep(dt=base*(1+0.2*sin(2*pi*integrator.Timeline.Now)));
                        integrator.update(value);
                        value = integrator.step(L=-1, S=@(~, t) 2*t+t.^2, M=1);
                        integrator.advance();
                    end
                    errors(level) = abs(value-(exp(-1)+1));
                end
                testCase.verifyGreaterThan(log(errors(2)/errors(3))/log(2), method+1-0.4);
            end
        end

        function legacyBdf3StartupRetainsItsDocumentedGlobalOrder(testCase)
            errors = zeros(1, 3);
            for level = 1:3
                integrator = mdemat.common.ode.Bdf3Integrator(1);
                integrator.Startup = "legacy";
                value = TestOdeScientific.advance(integrator, 1, 16*2^(level-1), L=-1, M=1);
                errors(level) = abs(value-exp(-1));
            end
            rate = log(errors(2)/errors(3))/log(2);
            testCase.verifyGreaterThan(rate, 1.7);
            testCase.verifyLessThan(rate, 2.3);
        end

        function explicitStepOverrideInvalidatesImplicitFactorization(testCase)
            integrator = mdemat.common.ode.BeIntegrator(1);
            % Fill the complete timeline queue so HasStepSizeChanged is
            % false before the explicit override; a startup-only test misses
            % the stale factorization after accepted constant-size steps.
            for step = 1:3
                integrator.setTimeStep(dt=0.1);
                integrator.update(1);
                integrator.step(L=-1, M=1);
                integrator.advance();
            end
            integrator.setTimeStep(dt=0.1);
            integrator.update(1);
            value = integrator.step(L=-1, M=1, dt=0.25);
            testCase.verifyEqual(value, 1/1.25, 'AbsTol', 1e-13);
        end
    end

    methods (Static, Access = private)
        function [names, orders] = catalog()
            names = ["FeIntegrator", "EmpIntegrator", "HeunIntegrator", "Ssprk2Integrator", ...
                "Ssprk3Integrator", "Exrk4Integrator", "BeIntegrator", "ImpIntegrator", ...
                "Sdirk2Integrator", "Sdirk3Integrator", "Sdirk4Integrator", "Esdirk3Integrator", ...
                "Bdf2Integrator", "Bdf3Integrator", "Ars111Integrator", "Ars222Integrator", "Ars443Integrator"];
            orders = [1, 2, 2, 2, 3, 4, 1, 2, 2, 3, 4, 3, 2, 3, 1, 2, 3];
        end

        function integrator = create(name, final)
            factory = str2func("mdemat.common.ode."+name);
            integrator = factory(final);
        end

        function value = advance(integrator, initial, steps, varargin)
            value = initial;
            integrator.reset();
            dt = integrator.Timeline.Final/steps;
            while ~integrator.Timeline.IsExhausted
                integrator.setTimeStep(dt=dt);
                integrator.update(value);
                value = integrator.step(varargin{:});
                integrator.advance();
            end
        end
    end
end
