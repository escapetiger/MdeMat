classdef TestFiniteDifferences < matlab.unittest.TestCase
    methods (Test)
        function quadraticAndBoundaryLift(test)
            [L, x, B] = mdemat.common.fdm.secondDerivative(8, [0, 1]);
            u = (-L) \ (2*ones(8, 1) + B*[1; 2]);
            test.verifyEqual(u, x.*(1-x) + 1+x, AbsTol=2e-14);
        end
        function manufacturedSecondOrder(test)
            sizes = [15, 31, 63]; errors = zeros(size(sizes));
            for index = 1:numel(sizes)
                [L, x] = mdemat.common.fdm.secondDerivative(sizes(index), [0, 1]);
                u = (-L) \ (pi^2 * sin(pi*x));
                errors(index) = norm(u-sin(pi*x), Inf);
            end
            rates = log2(errors(1:end-1)./errors(2:end));
            test.verifyGreaterThan(rates, [1.98, 1.98]);
            test.verifyLessThan(rates, [2.03, 2.03]);
        end
        function periodicNullspaceAndFourierMode(test)
            n = 24;
            [L, x, B, h] = mdemat.common.fdm.secondDerivative(n, [0, 2*pi], Boundary="periodic");
            test.verifyLessThan(norm(L*ones(n, 1)), 1e-13);
            test.verifyEqual(L*sin(x), -4*sin(h/2)^2/h^2*sin(x), AbsTol=2e-14);
            test.verifyEqual(nnz(B), 0);
            test.verifyError(@() mdemat.common.fdm.secondDerivative(2, [0, 1], ...
                Boundary="periodic"), 'mdemat:fdm:PeriodicSize');
        end
    end
end
