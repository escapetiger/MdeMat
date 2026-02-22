classdef TestInterval < matlab.unittest.TestCase
    % TESTINTERVAL Unit tests for the Interval class.
    %
    % This test class validates the functionality of the Interval class, which
    % represents one-dimensional intervals on the real line using the bbox
    % format [lower, upper]. Tests include construction with various
    % parameters, validation of boundary inclusion properties, testing points
    % for inside/boundary conditions, and verifying string representation.
    %
    % See also:
    %   core.geometry.Interval

    properties (Constant)
        Tolerance = sqrt(eps) % Numerical tolerance for floating-point comparisons
    end

    methods (Test)
        function testConstructor(testCase)
            % TESTCONSTRUCTOR Test interval constructor with bbox format.
            %
            % Verifies that the constructor properly handles the bbox input
            % format [lower, upper] and creates correct intervals with proper
            % boundary inclusion settings.

            % Default closed interval [0, 1]
            interval1 = core.geometry.Interval([0, 1]);
            testCase.verifyEqual(interval1.lower, 0);
            testCase.verifyEqual(interval1.upper, 1);
            testCase.verifyEqual(interval1.includeLower, true);
            testCase.verifyEqual(interval1.includeUpper, true);

            % Half-open interval [0, 1)
            interval2 = core.geometry.Interval([0, 1], 'IncludeUpper', false);
            testCase.verifyEqual(interval2.lower, 0);
            testCase.verifyEqual(interval2.upper, 1);
            testCase.verifyEqual(interval2.includeLower, true);
            testCase.verifyEqual(interval2.includeUpper, false);

            % Open interval (0, 1)
            interval3 = core.geometry.Interval([0, 1], 'IncludeLower', false, ...
                'IncludeUpper', false);
            testCase.verifyEqual(interval3.lower, 0);
            testCase.verifyEqual(interval3.upper, 1);
            testCase.verifyEqual(interval3.includeLower, false);
            testCase.verifyEqual(interval3.includeUpper, false);

            % Half-open interval (0, 1]
            interval4 = core.geometry.Interval([0, 1], 'IncludeLower', false);
            testCase.verifyEqual(interval4.lower, 0);
            testCase.verifyEqual(interval4.upper, 1);
            testCase.verifyEqual(interval4.includeLower, false);
            testCase.verifyEqual(interval4.includeUpper, true);

            % Test with non-zero bounds
            interval5 = core.geometry.Interval([-5, 10]);
            testCase.verifyEqual(interval5.lower, -5);
            testCase.verifyEqual(interval5.upper, 10);

            % Test with column vector input
            interval6 = core.geometry.Interval([2; 5]);
            testCase.verifyEqual(interval6.lower, 2);
            testCase.verifyEqual(interval6.upper, 5);
        end

        function testConstructorValidation(testCase)
            % TESTCONSTRUCTORVALIDATION Test constructor input validation.
            %
            % Verifies that the constructor properly validates input and
            % throws appropriate errors for invalid inputs.

            % Test empty input
            testCase.verifyError(@() core.geometry.Interval(), ?MException);

            % Test wrong number of elements
            testCase.verifyError(@() core.geometry.Interval([1]), ?MException);
            testCase.verifyError(@() core.geometry.Interval([1, 2, 3]), ?MException);

            % Test non-vector input
            testCase.verifyError(@() core.geometry.Interval([1, 2; 3, 4]), ?MException);

            % Test lower > upper bounds
            testCase.verifyError(@() core.geometry.Interval([1, 0]), ?MException);

            % Test invalid parameter values
            testCase.verifyError(@() core.geometry.Interval([0, 1], ...
                'IncludeLower', 'invalid'), ?MException);
            testCase.verifyError(@() core.geometry.Interval([0, 1], ...
                'IncludeUpper', [true, false]), ?MException);

            % Test valid edge cases
            % Equal bounds (degenerate interval)
            degenerateInterval = core.geometry.Interval([1, 1]);
            testCase.verifyEqual(degenerateInterval.lower, 1);
            testCase.verifyEqual(degenerateInterval.upper, 1);

            % Negative bounds
            negativeInterval = core.geometry.Interval([-2, -1]);
            testCase.verifyEqual(negativeInterval.lower, -2);
            testCase.verifyEqual(negativeInterval.upper, -1);

            % Infinite bounds
            infiniteInterval = core.geometry.Interval([-Inf, Inf]);
            testCase.verifyEqual(infiniteInterval.lower, -Inf);
            testCase.verifyEqual(infiniteInterval.upper, Inf);
        end

        function testInheritance(testCase)
            % TESTINHERITANCE Test inheritance hierarchy.
            %
            % Verifies that Interval correctly inherits from Orthotope
            % and Geometry classes.

            interval = core.geometry.Interval([0, 1]);
            testCase.verifyTrue(isa(interval, 'core.geometry.Orthotope'));
            testCase.verifyTrue(isa(interval, 'core.geometry.Geometry'));
            testCase.verifyEqual(interval.nDims, 1);
        end

        function testIsInside(testCase)
            % TESTISINSIDE Test interior point detection.
            %
            % Verifies that the isInside method correctly identifies points
            % strictly inside various interval types, accounting for
            % boundary inclusion settings.

            % Closed interval [0, 1]
            closedInterval = core.geometry.Interval([0, 1]);
            testCase.verifyEqual(closedInterval.isInside([0, 0.5, 1]), ...
                [false, true, false]);
            testCase.verifyEqual(closedInterval.isInside([-1, 2]), ...
                [false, false]);

            % Half-open interval [0, 1)
            leftClosedInterval = core.geometry.Interval([0, 1], 'IncludeUpper', false);
            testCase.verifyEqual(leftClosedInterval.isInside([0, 0.5, 1]), ...
                [false, true, false]);

            % Half-open interval (0, 1]
            rightClosedInterval = core.geometry.Interval([0, 1], 'IncludeLower', false);
            testCase.verifyEqual(rightClosedInterval.isInside([0, 0.5, 1]), ...
                [false, true, false]);

            % Open interval (0, 1)
            openInterval = core.geometry.Interval([0, 1], 'IncludeLower', false, ...
                'IncludeUpper', false);
            testCase.verifyEqual(openInterval.isInside([0, 0.5, 1]), ...
                [false, true, false]);

            % Test infinite upper bound
            infiniteInterval = core.geometry.Interval([0, Inf]);
            testCase.verifyEqual(infiniteInterval.isInside([0, 1, 1000]), ...
                [false, true, true]);
            testCase.verifyEqual(infiniteInterval.isInside(-1), false);

            % Test infinite lower bound
            negInfiniteInterval = core.geometry.Interval([-Inf, 0]);
            testCase.verifyEqual(negInfiniteInterval.isInside([-1000, -1, 0]), ...
                [true, true, false]);

            % Test dimension mismatch
            testCase.verifyError(@() closedInterval.isInside([0, 1; 2, 3]), ?MException);
        end

        function testIsOnBoundary(testCase)
            % TESTISONBOUNDARY Test boundary point detection.
            %
            % Verifies that the isOnBoundary method correctly identifies
            % points on the interval boundary, accounting for boundary
            % inclusion settings.

            % Closed interval [0, 1]
            closedInterval = core.geometry.Interval([0, 1]);
            testCase.verifyEqual(closedInterval.isOnBoundary([0, 0.5, 1]), ...
                [true, false, true]);
            testCase.verifyEqual(closedInterval.isOnBoundary([-1, 2]), ...
                [false, false]);

            % Half-open interval [0, 1)
            leftClosedInterval = core.geometry.Interval([0, 1], 'IncludeUpper', false);
            testCase.verifyEqual(leftClosedInterval.isOnBoundary([0, 0.5, 1]), ...
                [true, false, false]);

            % Half-open interval (0, 1]
            rightClosedInterval = core.geometry.Interval([0, 1], 'IncludeLower', false);
            testCase.verifyEqual(rightClosedInterval.isOnBoundary([0, 0.5, 1]), ...
                [false, false, true]);

            % Open interval (0, 1)
            openInterval = core.geometry.Interval([0, 1], 'IncludeLower', false, ...
                'IncludeUpper', false);
            testCase.verifyEqual(openInterval.isOnBoundary([0, 0.5, 1]), ...
                [false, false, false]);

            % Test with tolerance (points very close to boundary)
            tol = testCase.Tolerance / 2;
            testCase.verifyEqual(closedInterval.isOnBoundary([0 + tol, 1 - tol]), ...
                [true, true]);

            % Test infinite upper bound
            infiniteInterval = core.geometry.Interval([0, Inf]);
            testCase.verifyEqual(infiniteInterval.isOnBoundary([0, 1, 1000]), ...
                [true, false, false]);

            % Test infinite lower bound
            negInfiniteInterval = core.geometry.Interval([-Inf, 0]);
            testCase.verifyEqual(negInfiniteInterval.isOnBoundary([-1000, -1, 0]), ...
                [false, false, true]);

            % Test dimension mismatch
            testCase.verifyError(@() closedInterval.isOnBoundary([0, 1; 2, 3]), ?MException);
        end

        function testToString(testCase)
            % TESTTOSTRING Test string representation method.
            %
            % Verifies that the toString method produces correct mathematical
            % notation for various interval types and boundary conditions.

            % Closed interval [0, 1]
            closedInterval = core.geometry.Interval([0, 1]);
            testCase.verifyEqual(closedInterval.toString(), '[0,1]');

            % Half-open interval [0, 1)
            leftClosedInterval = core.geometry.Interval([0, 1], 'IncludeUpper', false);
            testCase.verifyEqual(leftClosedInterval.toString(), '[0,1)');

            % Half-open interval (0, 1]
            rightClosedInterval = core.geometry.Interval([0, 1], 'IncludeLower', false);
            testCase.verifyEqual(rightClosedInterval.toString(), '(0,1]');

            % Open interval (0, 1)
            openInterval = core.geometry.Interval([0, 1], 'IncludeLower', false, ...
                'IncludeUpper', false);
            testCase.verifyEqual(openInterval.toString(), '(0,1)');

            % Test with negative and non-integer bounds
            customInterval = core.geometry.Interval([-5.5, 10.75]);
            testCase.verifyEqual(customInterval.toString(), '[-5.5,10.75]');

            % Test infinite upper bound
            infiniteInterval = core.geometry.Interval([0, Inf]);
            testCase.verifyEqual(infiniteInterval.toString(), '[0,Inf]');

            % Test infinite lower bound
            negInfiniteInterval = core.geometry.Interval([-Inf, 0]);
            testCase.verifyEqual(negInfiniteInterval.toString(), '[-Inf,0]');

            % Test both bounds infinite
            bothInfiniteInterval = core.geometry.Interval([-Inf, Inf], ...
                'IncludeLower', false, 'IncludeUpper', false);
            testCase.verifyEqual(bothInfiniteInterval.toString(), '(-Inf,Inf)');
        end

        function testBoundaryAndInteriorConsistency(testCase)
            % TESTBOUNDARYANDINTERIORCONSISTENCY Test consistency between
            % inside and boundary detection.
            %
            % Verifies that a point cannot be both strictly inside and on
            % the boundary simultaneously for any interval type.

            % Test various interval types
            testIntervals = { ...
                core.geometry.Interval([0, 1]), ...
                core.geometry.Interval([0, 1], 'IncludeUpper', false), ...
                core.geometry.Interval([0, 1], 'IncludeLower', false), ...
                core.geometry.Interval([0, 1], 'IncludeLower', false, ...
                    'IncludeUpper', false) ...
            };

            testPoints = [0, 0.3, 0.5, 0.7, 1];

            for i = 1:length(testIntervals)
                interval = testIntervals{i};
                insideResults = interval.isInside(testPoints);
                boundaryResults = interval.isOnBoundary(testPoints);

                % Verify mutual exclusivity of inside and boundary
                for j = 1:length(testPoints)
                    testCase.verifyFalse(insideResults(j) && boundaryResults(j), ...
                        'Point cannot be both inside and on boundary for interval.');
                end
            end
        end

        function testOrthotopeMethods(testCase)
            % TESTORTHOTOPMETHODS Test inherited Orthotope functionality.
            %
            % Verifies that Interval properly inherits and overrides
            % Orthotope methods while maintaining interval-specific behavior.

            interval = core.geometry.Interval([0, 1]);

            % Test magnitude (should be length of interval)
            testCase.verifyEqual(interval.magnitude(), 1, 'AbsTol', testCase.Tolerance);

            % Test with different interval length
            interval2 = core.geometry.Interval([-2, 3]);
            testCase.verifyEqual(interval2.magnitude(), 5, 'AbsTol', testCase.Tolerance);

            % Test degenerate interval
            degenerateInterval = core.geometry.Interval([1, 1]);
            testCase.verifyEqual(degenerateInterval.magnitude(), 0, 'AbsTol', testCase.Tolerance);

            % Test that Interval retains Orthotope behavior for interior points
            testCase.verifyTrue(interval.isInside(0.5));
            testCase.verifyFalse(interval.isInside(2));

            % Test that boundary points are handled correctly by both classes
            testCase.verifyTrue(interval.isOnBoundary(0));
            testCase.verifyTrue(interval.isOnBoundary(1));
        end

        function testParameterParsing(testCase)
            % TESTPARAMETERPARSING Test parameter parsing functionality.
            %
            % Verifies that optional parameters are correctly parsed and
            % applied to interval boundary inclusion settings.

            % Test single parameter
            interval1 = core.geometry.Interval([0, 1], 'IncludeLower', false);
            testCase.verifyFalse(interval1.includeLower);
            testCase.verifyTrue(interval1.includeUpper);

            % Test both parameters in different orders
            interval2 = core.geometry.Interval([0, 1], 'IncludeUpper', false, ...
                'IncludeLower', true);
            testCase.verifyTrue(interval2.includeLower);
            testCase.verifyFalse(interval2.includeUpper);

            interval3 = core.geometry.Interval([0, 1], 'IncludeLower', false, ...
                'IncludeUpper', true);
            testCase.verifyFalse(interval3.includeLower);
            testCase.verifyTrue(interval3.includeUpper);

            % Test default values when no parameters specified
            interval4 = core.geometry.Interval([0, 1]);
            testCase.verifyTrue(interval4.includeLower);
            testCase.verifyTrue(interval4.includeUpper);
        end

        function testVectorInputHandling(testCase)
            % TESTVECTORINPUTHANDLING Test bbox vector input variations.
            %
            % Verifies that the constructor correctly handles both row and
            % column vectors for the bbox input.

            % Test row vector
            intervalRow = core.geometry.Interval([1, 2]);
            testCase.verifyEqual(intervalRow.lower, 1);
            testCase.verifyEqual(intervalRow.upper, 2);

            % Test column vector
            intervalCol = core.geometry.Interval([1; 2]);
            testCase.verifyEqual(intervalCol.lower, 1);
            testCase.verifyEqual(intervalCol.upper, 2);

            % Verify both produce identical results
            testCase.verifyEqual(intervalRow.lower, intervalCol.lower);
            testCase.verifyEqual(intervalRow.upper, intervalCol.upper);
            testCase.verifyEqual(intervalRow.includeLower, intervalCol.includeLower);
            testCase.verifyEqual(intervalRow.includeUpper, intervalCol.includeUpper);
        end
    end
end