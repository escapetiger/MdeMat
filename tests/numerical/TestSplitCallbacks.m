classdef TestSplitCallbacks < matlab.unittest.TestCase
    methods (Test)
        function hooksReceiveTheActualFactorsAndEffectiveSubsteps(testCase)
            left = [-1, 1; 0, -2];
            right = [0, 0; 1, -3];
            initial = [1; 2];
            stepSize = 0.2;
            for symmetric = [false, true]
                factors = {mdemat.common.ode.BeIntegrator(1), mdemat.common.ode.FeIntegrator(1)};
                if symmetric
                    composition = mdemat.common.ode.StrangIntegrator(factors=factors);
                    half = eye(2)-stepSize*left/2;
                    expected = half\((eye(2)+stepSize*right)*(half\initial));
                    expectedIndices = [1, 2, 1];
                    expectedSteps = [stepSize/2, stepSize, stepSize/2];
                else
                    composition = mdemat.common.ode.LieTrotterIntegrator(factors=factors);
                    expected = (eye(2)+stepSize*right)*((eye(2)-stepSize*left)\initial);
                    expectedIndices = [1, 2];
                    expectedSteps = [stepSize, stepSize];
                end
                composition.add(initial);
                composition.setTimeStep(dt=stepSize);
                seenIndices = [];
                seenSteps = [];
                value = composition.step(L={left, []}, F={[], @(u) right*u}, ...
                    fn={@(factor, values) apply(factor, values, 1), ...
                        @(factor, values) apply(factor, values, 2)});
                testCase.verifyEqual(value, expected, 'AbsTol', 2e-14);
                testCase.verifyEqual(seenIndices, expectedIndices);
                testCase.verifyEqual(seenSteps, expectedSteps);
                testCase.verifyEqual(composition.U0{1}, initial);
                testCase.verifyEqual(factors{2}.U0{1}, initial);
            end

            function value = apply(factor, values, index)
                testCase.verifySameHandle(factor, factors{index});
                seenIndices(end+1) = index;
                seenSteps(end+1) = values.dt;
                % This composition has no algebraic constraints; FE accepts
                % only its numerical operators and the effective substep.
                value = factor.step(L=values.L, F=values.F, S=values.S, ...
                    M=values.M, dt=values.dt);
            end
        end
    end
end
