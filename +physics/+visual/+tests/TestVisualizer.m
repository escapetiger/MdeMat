classdef TestVisualizer < matlab.unittest.TestCase
    % TESTVISUALIZER Unit tests for the Visualizer class.
    %
    %   TestVisualizer provides comprehensive test coverage for the
    %   Visualizer class functionality including constructor validation,
    %   setter methods, state-dependent plotting, and configuration
    %   management.
    %
    % See also:
    %   physics.visual.Visualizer

    properties (Constant)
        Tolerance = 1e-12
    end

    methods (Test)
        function testEmptyConstructor(testCase)
            % Test empty constructor functionality
            vis = physics.visual.Visualizer();
            testCase.verifyEmpty(vis.NDims);
            testCase.verifyEmpty(vis.Timeline);
            testCase.verifyEmpty(vis.Density);
            testCase.verifyClass(vis.Components, 'struct');
            testCase.verifyClass(vis.Exacts, 'struct');
            testCase.verifyClass(vis.Plotters, 'struct');
            testCase.verifyEqual(vis.TitleHead, '');
            testCase.verifyEqual(vis.Status, 0b000);
            testCase.verifyClass(vis.FigureManager, 'physics.visual.FigureManager');
            testCase.verifyClass(vis.Database, 'physics.visual.Database');
        end

        function testConstructorWithDimensions(testCase)
            % Test constructor with spatial dimensions
            vis = physics.visual.Visualizer(2);
            testCase.verifyEqual(vis.NDims, 2);
            testCase.verifyEmpty(vis.Timeline);
            testCase.verifyEmpty(vis.Density);
        end

        function testInvalidDimensionsValidation(testCase)
            % Test input validation for dimensions
            testCase.verifyError(@() physics.visual.Visualizer(-1), ...
                'MATLAB:validators:mustBePositive');
            testCase.verifyError(@() physics.visual.Visualizer(2.5), ...
                'MATLAB:validators:mustBeInteger');
        end

        function testSetDensity(testCase)
            % Test setDensity method
            vis = physics.visual.Visualizer();
            density = [10, 15, 20];
            vis = vis.setDensity(density);
            testCase.verifyEqual(vis.Density, density);
        end

        function testInvalidDensityValidation(testCase)
            % Test density validation
            vis = physics.visual.Visualizer();
            testCase.verifyError(@() vis.setDensity([-1, 5]), ...
                'MATLAB:validators:mustBePositive');
            testCase.verifyError(@() vis.setDensity([1; 2; 3]), ...
                'MATLAB:validators:mustBeVector');
        end

        function testSetTimeline(testCase)
            % Test setTimeline method
            vis = physics.visual.Visualizer();
            finalTime = 10.0;
            nTimeNodes = 50;
            vis = vis.setTimeline(finalTime, nTimeNodes);
            testCase.verifyClass(vis.Timeline, 'approx.mesh.StaticTimeline');
        end

        function testInvalidTimelineValidation(testCase)
            % Test timeline parameter validation
            vis = physics.visual.Visualizer();
            testCase.verifyError(@() vis.setTimeline(-1, 10), ...
                'MATLAB:validators:mustBeNonnegative');
            testCase.verifyError(@() vis.setTimeline(10, 0), ...
                'MATLAB:validators:mustBePositive');
            testCase.verifyError(@() vis.setTimeline(10, 2.5), ...
                'MATLAB:validators:mustBeInteger');
        end

        function testSetComponents(testCase)
            % Test setComponents method
            vis = physics.visual.Visualizer();
            components = struct('u', [1, 2], 'v', [1]);
            vis = vis.setComponents(components);
            testCase.verifyEqual(vis.Components, components);
        end

        function testSetTitleHead(testCase)
            % Test setTitleHead method
            vis = physics.visual.Visualizer();
            titleHead = "Test_Title_With_Underscores";
            vis = vis.setTitleHead(titleHead);
            testCase.verifyEqual(vis.TitleHead, "Test-Title-With-Underscores");
        end

        function testAddDataset(testCase)
            % Test addDataset method
            vis = physics.visual.Visualizer();
            vis = vis.addDataset("testData", 1);
            % Note: This tests the interface. Actual verification would
            % require checking the Database state.
            testCase.verifyClass(vis, 'physics.visual.Visualizer');
        end

        function testInvalidDatasetValidation(testCase)
            % Test dataset parameter validation
            vis = physics.visual.Visualizer();
            testCase.verifyError(@() vis.addDataset("test", -1), ...
                'MATLAB:validators:mustBeNonnegative');
            testCase.verifyError(@() vis.addDataset("test", 1.5), ...
                'MATLAB:validators:mustBeInteger');
        end

        function testAddExact(testCase)
            % Test addExact method
            vis = physics.visual.Visualizer();
            exactFunc = @(x, y) sin(x) .* cos(y);
            vis = vis.addExact("testExact", exactFunc);
            testCase.verifyEqual(vis.Exacts.testExact, exactFunc);
        end

        function testInvalidExactValidation(testCase)
            % Test exact solution validation
            vis = physics.visual.Visualizer();
            testCase.verifyError(@() vis.addExact("test", "not_a_function"), ...
                'MATLAB:validators:UnableToConvert');
        end

        function testIsEnabledProperty(testCase)
            % Test IsEnabled dependent property
            vis = physics.visual.Visualizer();
            testCase.verifyFalse(vis.IsEnabled);
            
            % Set required properties
            vis = vis.setDensity([10, 10]);
            vis = vis.setComponents(struct('u', [1]));
            
            % Should still be false until properly configured
            testCase.verifyFalse(vis.IsEnabled);
        end

        function testHasExactProperty(testCase)
            % Test HasExact dependent property
            vis = physics.visual.Visualizer();
            testCase.verifyFalse(vis.HasExact);
            
            exactFunc = @(x, y) x.^2 + y.^2;
            vis = vis.addExact("test", exactFunc);
            testCase.verifyTrue(vis.HasExact);
        end

        function testNComponentsProperty(testCase)
            % Test NComponents dependent property
            vis = physics.visual.Visualizer();
            testCase.verifyEmpty(vis.NComponents);
            
            components = struct('u', [1, 2], 'v', [1], 'p', [1, 2, 3]);
            vis = vis.setComponents(components);
            expectedComponents = [2, 1, 3];
            testCase.verifyEqual(vis.NComponents, expectedComponents);
        end

        function testPlotMethodValidation(testCase)
            % Test plot method input validation
            vis = physics.visual.Visualizer();
            
            % Should throw error for invalid state type
            invalidState = struct('invalid', true);
            testCase.verifyError(@() vis.plot(invalidState), ...
                'core:except:InvalidInput');
        end
    end
end