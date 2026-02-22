classdef TestTiledPlotter < matlab.unittest.TestCase
    % TESTTILEDPLOTTER Unit tests for the TiledPlotter class.
    %
    %   TestTiledPlotter provides test coverage for the TiledPlotter
    %   class functionality including constructor validation, grid
    %   configuration, and strategy management.
    %
    % See also:
    %   physics.visual.TiledPlotter

    methods (Test)
        function testConstructor(testCase)
            % Test TiledPlotter constructor
            plotter = physics.visual.TiledPlotter(2, 3);
            testCase.verifyClass(plotter, 'physics.visual.TiledPlotter');
            testCase.verifyEqual(plotter.Size, [2, 3]);
            testCase.verifyEmpty(plotter.Specs);
            testCase.verifyEmpty(plotter.FigObj);
            testCase.verifyEmpty(plotter.Axes);
            testCase.verifyEmpty(plotter.Strategy);
        end

        function testConstructorValidation(testCase)
            % Test constructor input validation
            testCase.verifyError(@() physics.visual.TiledPlotter(0, 2), ...
                'MATLAB:validators:mustBePositive');
            testCase.verifyError(@() physics.visual.TiledPlotter(2.5, 2), ...
                'MATLAB:validators:mustBeInteger');
            testCase.verifyError(@() physics.visual.TiledPlotter(2, 0), ...
                'MATLAB:validators:mustBePositive');
        end

        function testSetStrategy(testCase)
            % Test setStrategy method
            plotter = physics.visual.TiledPlotter(1, 1);
            strategy = physics.visual.Strategy();
            plotter = plotter.setStrategy(strategy);
            testCase.verifyEqual(plotter.Strategy, strategy);
        end

        function testAddAxisSpec(testCase)
            % Test addAxisSpec method
            plotter = physics.visual.TiledPlotter(3, 3);
            plotter = plotter.addAxisSpec(1, 2); % 1x2 panel
            plotter = plotter.addAxisSpec(2, 1); % 2x1 panel
            
            testCase.verifyEqual(length(plotter.Specs), 2);
            testCase.verifyEqual(plotter.Specs(1).rows, 1);
            testCase.verifyEqual(plotter.Specs(1).cols, 2);
            testCase.verifyEqual(plotter.Specs(2).rows, 2);
            testCase.verifyEqual(plotter.Specs(2).cols, 1);
        end

        function testBindVariable(testCase)
            % Test bind method for variable-to-axis mapping
            plotter = physics.visual.TiledPlotter(2, 2);
            plotter = plotter.bind('temperature', 1);
            plotter = plotter.bind('pressure', 2);
            
            testCase.verifyEqual(plotter.VarToAxis.temperature, 1);
            testCase.verifyEqual(plotter.VarToAxis.pressure, 2);
        end

        function testSetFigure(testCase)
            % Test setFigure method
            plotter = physics.visual.TiledPlotter(1, 1);
            fig = figure('Visible', 'off');
            plotter = plotter.setFigure(fig);
            
            testCase.verifyEqual(plotter.FigObj, fig);
            
            % Clean up
            close(fig);
        end
    end
end