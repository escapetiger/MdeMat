classdef TestSymbolicFamilies < matlab.unittest.TestCase
    methods (TestClassSetup)
        function requireSymbolic(test)
            test.assertTrue(logical(license('test','Symbolic_Toolbox')), ...
                'This dedicated feature suite requires Symbolic Math Toolbox.');
        end
    end
    methods (Test)
        function classicalPolynomialValuesAndDerivativesAgree(test)
            variable = sym('x','real');
            points = [-0.7,-0.1,0.3,0.8];
            families = {'LegendreBasis','ChebyshevBasis','HermiteBasis'};
            for family = families
                for normalized = [false,true]
                    symbolic = mdemat.common.symbolic.(family{1})(variable,5,isNormalized=normalized);
                    direct = mdemat.common.spectral.(family{1})(5,isNormalized=normalized);
                    symbolic.compile(2);
                    values = evaluate(symbolic,points,1);
                    test.verifyEqual(values,direct.eval(points),AbsTol=3e-13);
                    h = 1e-5;
                    derivative = (direct.eval(points+h)-direct.eval(points-h))/(2*h);
                    test.verifyEqual(evaluate(symbolic,points,2),derivative,AbsTol=8e-8);
                end
            end
        end

        function shiftedIntervalsKeepEndpointsAndPhysicalNormalization(test)
            variable = sym('x','real');
            lower = 2; upper = 5; points = linspace(lower,upper,13);
            mapped = 2*(points-(lower+upper)/2)/(upper-lower);
            for family = {'LegendreBasis','ChebyshevBasis'}
                symbolic = mdemat.common.symbolic.(family{1})(variable,4, ...
                    isNormalized=true,lower=lower,upper=upper);
                symbolic.compile();
                direct = mdemat.common.spectral.(family{1})(4,isNormalized=true);
                test.verifyEqual(evaluate(symbolic,points,1), ...
                    sqrt(2/(upper-lower))*direct.eval(mapped),AbsTol=3e-13);
                test.verifyError(@() mdemat.common.symbolic.(family{1})(variable,2, ...
                    lower=3,upper=2),'MdeMat:Symbolic:Interval');
            end
            test.verifyError(@() mdemat.common.symbolic.LegendreBasis(variable,2, ...
                isNormalized=true,isMonic=true),'MdeMat:Symbolic:Normalization');
        end

        function fourierModesMatchAnIndependentPhysicalInterval(test)
            variable = sym('x','real'); lower = -0.7; upper = 2.1;
            points = linspace(lower,upper,29);
            for normalized = [false,true]
                symbolic = mdemat.common.symbolic.FourierBasis(variable,3, ...
                    lower=lower,upper=upper,isNormalized=normalized);
                direct = mdemat.common.spectral.FourierBasis(3, ...
                    lower=lower,upper=upper,isNormalized=normalized);
                symbolic.compile(1);
                test.verifyEqual(evaluate(symbolic,points,1),direct.eval(points),AbsTol=3e-13);
            end
        end
    end
end

function values = evaluate(basis,points,column)
values = zeros(basis.NCodims,numel(points));
for row = 1:basis.NCodims, values(row,:)=basis.Handles{row,column}(points); end
end
