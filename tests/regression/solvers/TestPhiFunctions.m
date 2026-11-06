classdef TestPhiFunctions < matlab.unittest.TestCase
    %TESTPHIFUNCTIONS Independent entire-function and source-response checks.
    methods (Test)
        function zeroAndNilpotentMatricesHaveTheirExactFiniteSeries(testCase)
            sequence = mdemat.common.ode.expseq(zeros(3), 6);
            for order = 0:6
                testCase.verifyEqual(sequence{order+1}, eye(3)/factorial(order), 'AbsTol', 2e-16);
            end
            % Large norm with zero eigenvalues: division by this matrix
            % cannot recover the analytic functions on its kernel.
            matrix = [0, 100, 0; 0, 0, 0.3; 0, 0, 0];
            sequence = mdemat.common.ode.expseq(matrix, 5);
            for order = 0:5
                exact = eye(3)/factorial(order)+matrix/factorial(order+1) ...
                    +matrix^2/factorial(order+2);
                testCase.verifyEqual(sequence{order+1}, exact, 'AbsTol', 3e-12);
            end
        end

        function scalarComplexValuesAgreeWithIndependentIntegralDefinitions(testCase)
            values = [1e-12, -0.49, -0.51, 0.2+0.1i, -2+3i, -30];
            for value = values
                sequence = mdemat.common.ode.expseq(value, 5);
                testCase.verifyEqual(sequence{1}, exp(value), 'AbsTol', 2e-13);
                for order = 1:5
                    exact = integral(@(s) exp((1-s)*value).*s.^(order-1)/factorial(order-1), ...
                        0, 1, 'AbsTol', 2e-14, 'RelTol', 2e-13);
                    testCase.verifyEqual(sequence{order+1}, exact, 'AbsTol', 3e-13);
                end
            end
        end

        function singularAndNonnormalMatricesMatchTheirPolynomialSourceResponse(testCase)
            matrices = {[0, 4; 0, -2], [-8, 20; -1, -3], [-20, 50; 0, -20]};
            for index = 1:numel(matrices)
                matrix = matrices{index};
                sequence = mdemat.common.ode.expseq(matrix, 4);
                for order = 1:4
                    % Variation of constants for u'=Z*u+t^(k-1)*e_j/(k-1)!.
                    % Quadrature of the original small exponential is an
                    % independent reference, with no augmented exponential.
                    exact = integral(@(s) expm((1-s)*matrix)*s^(order-1)/factorial(order-1), ...
                        0, 1, 'ArrayValued', true, 'AbsTol', 1e-13, 'RelTol', 1e-12);
                    testCase.verifyEqual(sequence{order+1}, exact, 'AbsTol', 2e-12);
                    testCase.verifyEqual(matrix*sequence{order+1}, ...
                        sequence{order}-eye(2)/factorial(order-1), 'AbsTol', 3e-12);
                end
            end
        end

        function sparseSingleAndStiffScalarInputsKeepDefinedValues(testCase)
            sequence = mdemat.common.ode.expseq(sparse(diag([-1, 0])), 2);
            testCase.verifyEqual(sequence{2}, diag([1-exp(-1), 1]), 'AbsTol', 2e-14);
            sequence = mdemat.common.ode.expseq(single(0.1), 2);
            testCase.verifyClass(sequence{1}, 'single');
            testCase.verifyEqual(sequence{2}, single(expm1(0.1)/0.1), 'AbsTol', single(2e-7));
            sequence = mdemat.common.ode.expseq(-Inf, 4);
            testCase.verifyEqual(cell2mat(sequence), zeros(1, 5));
            sequence = mdemat.common.ode.expseq(1, 0);
            testCase.verifyEqual(sequence{1}, exp(1), 'AbsTol', 1e-14);
        end

        function invalidOrdersAndNonfiniteMatricesAreRejected(testCase)
            for order = {-1, 0.5, Inf, NaN, [1, 2], 1i, '2'}
                testCase.verifyError(@() mdemat.common.ode.expseq(1, order{1}), ...
                    'mdemat:ode:InvalidIndex');
            end
            testCase.verifyError(@() mdemat.common.ode.expseq([], 1), 'mdemat:ode:InvalidMatrix');
            testCase.verifyError(@() mdemat.common.ode.expseq('x', 1), 'mdemat:ode:InvalidMatrix');
            testCase.verifyError(@() mdemat.common.ode.expseq(ones(2, 3), 1), 'mdemat:ode:NonSquareMatrix');
            for matrix = {NaN, Inf, diag([-Inf, -Inf]), [0, NaN; 0, 0]}
                testCase.verifyError(@() mdemat.common.ode.expseq(matrix{1}, 1), ...
                    'mdemat:ode:NonFiniteMatrix');
            end
        end
    end
end
