classdef TestFigureManager < matlab.unittest.TestCase
    % TESTFIGUREMANAGER Unit tests for the FigureManager class.
    %
    %   TestFigureManager provides test coverage for the FigureManager
    %   class functionality including constructor, figure creation,
    %   and cleanup operations.
    %
    % See also:
    %   physics.visual.FigureManager

    methods (Test)
        function testConstructor(testCase)
            % Test FigureManager constructor
            manager = physics.visual.FigureManager();
            testCase.verifyClass(manager, 'physics.visual.FigureManager');
            testCase.verifyEmpty(manager.FigureHandles);
            testCase.verifyClass(manager.ScreenSize, 'double');
            testCase.verifyClass(manager.FigureSize, 'double');
            testCase.verifyEqual(length(manager.FigureSize), 2);
        end

        function testConstants(testCase)
            % Test constant properties
            manager = physics.visual.FigureManager();
            testCase.verifyEqual(manager.DefaultFigureWidth, 560);
            testCase.verifyEqual(manager.DefaultFigureHeight, 420);
            testCase.verifyGreaterThan(manager.FigureSpacing, 0);
            testCase.verifyGreaterThan(manager.ScreenMargin, 0);
        end

        function testCreateFiguresValidation(testCase)
            % Test input validation for createFigures
            manager = physics.visual.FigureManager();
            testCase.verifyError(@() manager.createFigures(0), ...
                'MATLAB:validators:mustBePositive');
            testCase.verifyError(@() manager.createFigures(1.5), ...
                'MATLAB:validators:mustBeInteger');
        end

        function testCreateSingleFigure(testCase)
            % Test creating a single figure (with cleanup)
            manager = physics.visual.FigureManager();
            figures = manager.createFigures(1);
            
            testCase.addTeardown(@() manager.closeAll());
            
            testCase.verifyClass(figures, 'cell');
            testCase.verifyEqual(length(figures), 1);
            testCase.verifyClass(figures{1}, 'matlab.ui.Figure');
            testCase.verifyEqual(length(manager.FigureHandles), 1);
        end

        function testCreateMultipleFigures(testCase)
            % Test creating multiple figures
            manager = physics.visual.FigureManager();
            figures = manager.createFigures(3);
            
            testCase.addTeardown(@() manager.closeAll());
            
            testCase.verifyEqual(length(figures), 3);
            testCase.verifyEqual(length(manager.FigureHandles), 3);
            
            % Verify all figures are valid
            for i = 1:3
                testCase.verifyClass(figures{i}, 'matlab.ui.Figure');
                testCase.verifyTrue(isvalid(figures{i}));
            end
        end

        function testIsValid(testCase)
            % Test isValid method
            manager = physics.visual.FigureManager();
            
            % Should be valid with no figures
            testCase.verifyTrue(manager.isValid());
            
            figures = manager.createFigures(2);
            testCase.addTeardown(@() manager.closeAll());
            
            % Should be valid with open figures
            testCase.verifyTrue(manager.isValid());
            
            % Close one figure manually
            close(figures{1});
            
            % Should now be invalid
            testCase.verifyFalse(manager.isValid());
        end

        function testBringToFront(testCase)
            % Test bringToFront method
            manager = physics.visual.FigureManager();
            figures = manager.createFigures(2);
            
            testCase.addTeardown(@() manager.closeAll());
            
            % Should not error when bringing valid figure to front
            manager.bringToFront(1);
            manager.bringToFront(2);
            
            % Should error for invalid index
            testCase.verifyError(@() manager.bringToFront(0), ...
                'MATLAB:validators:mustBePositive');
            testCase.verifyError(@() manager.bringToFront(3), ...
                'physics:visual:FigureManager:InvalidIndex');
        end

        function testBringAllToFront(testCase)
            % Test bringAllToFront method
            manager = physics.visual.FigureManager();
            figures = manager.createFigures(3);
            
            testCase.addTeardown(@() manager.closeAll());
            
            % Should not error
            manager.bringAllToFront();
        end

        function testCloseAll(testCase)
            % Test closeAll method
            manager = physics.visual.FigureManager();
            figures = manager.createFigures(2);
            
            % Verify figures exist
            testCase.verifyEqual(length(manager.FigureHandles), 2);
            testCase.verifyTrue(isvalid(figures{1}));
            testCase.verifyTrue(isvalid(figures{2}));
            
            % Close all
            manager.closeAll();
            
            % Verify cleanup
            testCase.verifyEmpty(manager.FigureHandles);
            testCase.verifyFalse(isvalid(figures{1}));
            testCase.verifyFalse(isvalid(figures{2}));
        end

        function testLayoutCalculation(testCase)
            % Test layout calculation for different figure counts
            manager = physics.visual.FigureManager();
            
            % Test various layouts (using protected method access)
            testCase.verifyEqual(manager.calculateLayout(1), [1, 1]);
            testCase.verifyEqual(manager.calculateLayout(2), [1, 2]);
            testCase.verifyEqual(manager.calculateLayout(4), [2, 2]);
            testCase.verifyEqual(manager.calculateLayout(6), [2, 3]);
            testCase.verifyEqual(manager.calculateLayout(9), [3, 3]);
            
            % Test larger numbers produce reasonable layouts
            layout = manager.calculateLayout(10);
            testCase.verifyGreaterThanOrEqual(layout(1) * layout(2), 10);
        end

        function testScreenSizeAccess(testCase)
            % Test screen size retrieval
            manager = physics.visual.FigureManager();
            
            screenSize = manager.getScreenSize();
            testCase.verifyClass(screenSize, 'double');
            testCase.verifyEqual(length(screenSize), 2);
            testCase.verifyGreaterThan(screenSize(1), 0);
            testCase.verifyGreaterThan(screenSize(2), 0);
        end
    end
end