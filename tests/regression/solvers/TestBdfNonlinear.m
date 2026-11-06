classdef TestBdfNonlinear < matlab.unittest.TestCase
    %TESTBDFNONLINEAR Verify the public F contract, startup and saved history.
    methods (Test)
        function nonlinearDecayRetainsTheDeclaredGlobalOrder(testCase)
            for order = 2:3
                errors = zeros(1, 3);
                for level = 1:3
                    integrator = TestBdfNonlinear.create(order, 0.5);
                    value = TestBdfNonlinear.advance(integrator, 1, 8*2^(level-1), F=@(u) -u.^2);
                    errors(level) = abs(value-1/1.5);
                    testCase.verifyTrue(integrator.SolveReports{end}.Converged);
                end
                testCase.verifyGreaterThan(log(errors(2)/errors(3))/log(2), order-0.35);
            end
        end
        function nonlinearProductMassEquationUsesActualStageTimes(testCase)
            % Exact u=1/(1+t^2), M=1+t. Differentiate M*u independently.
            for order = 2:3
                errors = zeros(1, 3);
                for level = 1:3
                    integrator = TestBdfNonlinear.create(order, 0.5);
                    value = TestBdfNonlinear.advance(integrator, 1, 8*2^(level-1), ...
                        M=@(t) 1+t, F=@(u) -u.^2, S=@(~, t) (2-2*t-t.^2)./(1+t.^2).^2);
                    errors(level) = abs(value-0.8);
                end
                testCase.verifyGreaterThan(log(errors(2)/errors(3))/log(2), order-0.4);
            end
        end
        function nonlinearStagesSatisfyTheAugmentedConstraint(testCase)
            for order = 2:3
                integrator = TestBdfNonlinear.create(order, 0.2);
                value = TestBdfNonlinear.advance(integrator, [1; 1], 32, ...
                    M=2, F=@(u) -2*u(1).^2, P=0, A=-1, B=1, C=0);
                testCase.verifyEqual(value(1), 1/1.2, 'AbsTol', 5e-5);
                testCase.verifyEqual(value(2), value(1), 'AbsTol', 2e-11);
            end
        end
        function nonlinearContinuationRestoresTheSameBdfHistory(testCase)
            mass = diag([2, 3]);
            nonlinear = @(u) -mass*(u.^2);
            for order = 2:3
                whole = TestBdfNonlinear.advance(TestBdfNonlinear.create(order, 0.4), ...
                    [1; 2], 20, M=mass, F=nonlinear);
                integrator = TestBdfNonlinear.create(order, 0.4);
                integrator.reset();
                value = [1; 2];
                steps = 0;
                while ~integrator.Timeline.IsExhausted
                    integrator.setTimeStep(dt=0.02);
                    integrator.update(value);
                    value = integrator.step(M=mass, F=nonlinear);
                    integrator.advance();
                    steps = steps+1;
                    if steps == 7
                        record = mdemat.common.ode.snapshotIntegrator(integrator);
                        integrator = mdemat.common.ode.restoreIntegrator(record);
                    end
                end
                testCase.verifyEqual(value, whole, 'AbsTol', 2e-12);
            end
        end
    end
    methods (Static, Access = private)
        function integrator = create(order, final)
            if order == 2
                integrator = mdemat.common.ode.Bdf2Integrator(final);
            else
                integrator = mdemat.common.ode.Bdf3Integrator(final);
            end
        end
        function value = advance(integrator, value, steps, varargin)
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
