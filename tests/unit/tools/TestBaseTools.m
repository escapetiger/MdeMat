classdef TestBaseTools < matlab.unittest.TestCase
    methods (Test)
        function fieldChecksVisitEveryAxisAndHaveConsistentOutputs(test)
            value = struct('a', 1, 'b', 2, 'c', 3);
            names = {'a', 'c'; 'b', 'missing'};
            single = mdemat.tools.base.hasFields(value, names);
            [multiple, missing] = mdemat.tools.base.hasFields(value, names);
            test.verifyClass(single, 'logical');
            test.verifySize(single, [1, 1]);
            test.verifyFalse(single);
            test.verifyEqual(multiple, single);
            test.verifyEqual(missing, 'missing');
            test.verifyTrue(mdemat.tools.base.hasFields(value, ["a", "b"]));
            test.verifyTrue(mdemat.tools.base.hasFields(value, {}));
            test.verifyError(@() mdemat.tools.base.hasFields(value, {''}), ...
                'mdemat:tools:base:hasFields:InvalidInput');
            test.verifyError(@() mdemat.tools.base.hasFields(value, {42}), ...
                'mdemat:tools:base:hasFields:InvalidInput');
        end

        function classChecksDoNotSkipTrailingElements(test)
            values = {1, 2; 3, 'bad'};
            test.verifyFalse(mdemat.tools.base.isAllClass(values, 'double'));
            test.verifyTrue(mdemat.tools.base.isAllClass(values, ["double", "char"]));
            test.verifyTrue(mdemat.tools.base.isAllClass([], {}));
            test.verifyFalse(mdemat.tools.base.isAllClass(1, {}));
            timer = mdemat.tools.analyze.Timer();
            test.verifyTrue(mdemat.tools.base.isAllClass({timer}, 'handle'));
            test.verifyTrue(mdemat.tools.base.isAllSameClass());
            test.verifyTrue(mdemat.tools.base.isAllSameClass(uint8([]), uint8(2)));
            test.verifyFalse(mdemat.tools.base.isAllSameClass(uint8([]), []));
        end

        function heterogeneousArraysUseAllowedClassesPerElement(test)
            location = fullfile(fileparts(mfilename('fullpath')), 'fixtures');
            test.applyFixture(matlab.unittest.fixtures.PathFixture(location));
            values = [foundationprobe.Alpha(), foundationprobe.Beta()];
            allowed = {'foundationprobe.Alpha', 'foundationprobe.Beta'};
            test.verifyTrue(mdemat.tools.base.isAllClass(values, allowed));
            test.verifyFalse(mdemat.tools.base.isAllClass(values, allowed(1)));
        end

        function validatorsCheckTheContainerRatherThanItsContents(test)
            mdemat.tools.base.mustBeCellOrNumeric({@sin, 'anything'});
            mdemat.tools.base.mustBeCellOrObject({42});
            test.verifyError(@() mdemat.tools.base.mustBeCellOrNumeric('text'), ...
                'mdemat:tools:base:mustBeCellOrNumeric:InvalidInput');
            test.verifyError(@() mdemat.tools.base.mustBeCellOrObject(42), ...
                'mdemat:tools:base:mustBeCellOrObject:InvalidInput');
        end

        function factorOwnersKeepEveryCheckedElementAndRejectBadFactories(test)
            factors = mdemat.data.variable.Separable(factors={1,3;2,4});
            test.verifyEqual(factors.NFactors, 4);
            test.verifyEqual(factors.Factors, 1:4);
            test.verifyError(@() mdemat.data.variable.Separable.createFactors(42, 1:2), ...
                'mdemat:data:variable:Separable:InvalidFunctionHandle');
            test.verifyError(@() mdemat.data.variable.Separable.createFactors( ...
                @(first,second) first+second, reshape(1:4,2,2), 1:2), ...
                'mdemat:data:variable:Separable:InvalidInput');
            solver = mdemat.common.linalg.LinearSolver();
            test.verifyError(@() solver.precompute(), ...
                'mdemat:common:linalg:LinearSolver:InvalidState');
        end

        function digestsMatchIndependentLittleEndianByteFixtures(test)
            arrays = {uint8([127, 128, 255]), int16([1, -2; 32767, -32768]), ...
                complex(single([1, -3]), single([2, 4])), sparse([0, 2; 3, 0]), zeros(0, 3)};
            % Python hashlib over documented headers, NUL, and struct.pack bytes.
            expected = { ...
                'b3f6d3d3e6cb2e4c77143a71f38818894c4950ddf37b406f14739ef6dc3ed94d', ...
                'd284ca143cb3d9f0572d8f91675a3021ff55986c37ab5a4b9d853c5d8ccc583d', ...
                'f9c090c4be46d63354efbfb5b995dae0547fc30e2b725fab7a561aa104d3944b', ...
                '386f342f5607d57b8335f86d624b0ff5df94ce52f9946c7ee0a55d9f6a811f3c', ...
                'b81a48834ba0956692a24641ae3a72212cbb6d1b6b3178fd66cd065aac47e473'};
            for index = 1:numel(arrays)
                hex = mdemat.tools.base.hash(arrays{index});
                raw = mdemat.tools.base.hash(arrays{index}, format='uint8');
                test.verifyEqual(hex, expected{index});
                test.verifyEqual(reshape(lower(dec2hex(raw, 2)).', 1, []), hex);
                test.verifyEqual(mdemat.tools.base.hash(arrays{index}, format='double'), double(raw));
            end
            test.verifyEqual(mdemat.tools.base.hash(arrays{1}, method='SHA-1'), ...
                '26d3ea4098e18fcecda16ae17a3c46466651799e');
            test.verifyEqual(mdemat.tools.base.hash(arrays{1}, method='MD5'), ...
                'cc7d93cd96685ebde38d0043c4052ffa');
        end

        function matrixReuseDoesNotRequireDiagnosticHashing(test)
            solver = mdemat.common.linalg.LinearSolver();
            solver.setLhs([3, 1; 1, 2]);
            solver.precompute();
            solver.setLhs([3, 1; 1, 2]);
            test.verifyTrue(solver.IsPrecomputed);
            test.verifyEmpty(solver.LhsHash);
            diagnostic = mdemat.common.linalg.LinearSolver(recordHash=true);
            diagnostic.setLhs(solver.Lhs);
            test.verifyEqual(diagnostic.LhsHash, mdemat.tools.base.hash(solver.Lhs));
            test.verifyError(@() mdemat.tools.base.hash(1, method='unknown'), ...
                'mdemat:tools:base:hash:UnsupportedMethod');
            test.verifyError(@() mdemat.tools.base.hash(1, format='unknown'), ...
                'mdemat:tools:base:hash:UnsupportedFormat');
        end
    end
end
