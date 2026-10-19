classdef FigureManager < handle
    % FIGUREMANAGER Automatic figure positioning and lifecycle management.
    %
    %   FigureManager handles creation, positioning, and cleanup of figure
    %   windows with automatic screen centering. Figures are positioned
    %   in the center of the screen with default MATLAB styling.
    %
    %   The manager creates figures with consistent positioning and handles
    %   their lifecycle, ensuring proper cleanup when visualization is
    %   complete.
    %
    % Examples:
    %   % Create manager and figures
    %   manager = FigureManager();
    %   figures = manager.createFigures(2);
    %   
    %   % Later cleanup
    %   manager.closeAll();

    properties (SetAccess = protected, GetAccess = public)
        figureHandles % Cell array of created figure handles
        screenSize % Screen dimensions [width, height] in pixels
        figureSize % Standard figure size [width, height] in pixels
    end

    properties (Constant)
        DEFAULT_FIGURE_WIDTH = 560; % Default MATLAB figure width in pixels
        DEFAULT_FIGURE_HEIGHT = 420; % Default MATLAB figure height in pixels
        FIGURE_SPACING = 20; % Spacing between figures in pixels
        SCREEN_MARGIN = 50; % Margin from screen edges in pixels
    end

    methods
        function obj = FigureManager()
            % FIGUREMANAGER Constructor for FigureManager.
            %
            %   obj = FigureManager() creates a figure manager and
            %   initializes screen size information for automatic
            %   positioning calculations.
            %
            % Outputs:
            %   obj - Constructed FigureManager object

            obj.figureHandles = {};
            obj.screenSize = obj.getScreenSize();
            obj.figureSize = [obj.DEFAULT_FIGURE_WIDTH, obj.DEFAULT_FIGURE_HEIGHT];
        end

        function figureHandles = createFigures(obj, nFigures)
            % CREATEFIGURES Create non-overlapping figures with default styling.
            %
            %   figureHandles = createFigures(obj, nFigures) creates the
            %   specified number of figures positioned to avoid overlap.
            %   For 1-2 figures, they are placed side by side. For more
            %   figures, they are arranged in a simple grid pattern.
            %   Uses default MATLAB figure styling.
            %
            % Inputs:
            %   obj - The FigureManager object
            %   nFigures - Number of figures to create (positive integer)
            %
            % Outputs:
            %   figureHandles - Cell array of created figure handles

            core.except.assert(isnumeric(nFigures) && isscalar(nFigures) && nFigures > 0, ...
                'InvalidInput', 'Number of figures must be a positive integer');

            % Calculate layout for non-overlapping positioning
            layoutDims = obj.calculateLayout(nFigures);
            
            % Create positioned figures
            obj.figureHandles = cell(1, nFigures);
            for iFigure = 1:nFigures
                position = obj.calculateNonOverlappingPosition(iFigure, layoutDims);
                
                figHandle = figure('Position', position, ...
                    'NumberTitle', 'off', ...
                    'Name', sprintf('Figure %d', iFigure));
                
                obj.figureHandles{iFigure} = figHandle;
            end
            
            figureHandles = obj.figureHandles;
        end

        function obj = closeAll(obj)
            % CLOSEALL Close all managed figures.
            %
            %   obj = closeAll(obj) closes all figures created by this
            %   manager and clears the figure handle cache.
            %
            % Inputs:
            %   obj - The FigureManager object
            %
            % Outputs:
            %   obj - The FigureManager object

            for i = 1:length(obj.figureHandles)
                figHandle = obj.figureHandles{i};
                if isvalid(figHandle)
                    close(figHandle);
                end
            end
            obj.figureHandles = {};
        end

        function TF = isValid(obj)
            % ISVALID Check if all managed figures are still valid.
            %
            %   TF = isValid(obj) returns true if all managed figures
            %   are still open and valid.
            %
            % Inputs:
            %   obj - The FigureManager object
            %
            % Outputs:
            %   TF - True if all figures are valid (logical)

            TF = true;
            for i = 1:length(obj.figureHandles)
                if ~isvalid(obj.figureHandles{i})
                    TF = false;
                    return;
                end
            end
        end

        function obj = bringToFront(obj, figureIndex)
            % BRINGTOFRONT Bring specified figure to front.
            %
            %   obj = bringToFront(obj, figureIndex) brings the specified
            %   figure to the front of all windows and gives it focus.
            %
            % Inputs:
            %   obj - The FigureManager object
            %   figureIndex - Index of figure to bring to front
            %
            % Outputs:
            %   obj - The FigureManager object

            core.except.assert(figureIndex > 0 && figureIndex <= length(obj.figureHandles), ...
                'InvalidIndex', 'Figure index must be between 1 and %d', length(obj.figureHandles));

            figHandle = obj.figureHandles{figureIndex};
            if isvalid(figHandle)
                figure(figHandle);
            end
        end

        function obj = bringAllToFront(obj)
            % BRINGALLTOFRONT Bring all managed figures to front.
            %
            %   obj = bringAllToFront(obj) brings all managed figures
            %   to the front in creation order.
            %
            % Inputs:
            %   obj - The FigureManager object
            %
            % Outputs:
            %   obj - The FigureManager object

            for i = 1:length(obj.figureHandles)
                obj.bringToFront(i);
            end
        end
    end

    methods (Access = protected)
        function screenSize = getScreenSize(obj)
            % GETSCREENSIZE Get usable screen dimensions.
            %
            %   screenSize = getScreenSize(obj) returns the screen
            %   dimensions for positioning calculations.
            %
            % Outputs:
            %   screenSize - [width, height] in pixels

            screenSizeVec = get(0, 'ScreenSize'); % [left, bottom, width, height]
            screenSize = [screenSizeVec(3), screenSizeVec(4)];
        end

        function position = calculateCenteredPosition(obj, figureIndex)
            % CALCULATECENTEREDPOSITION Calculate centered position for figure.
            %
            %   position = calculateCenteredPosition(obj, figureIndex)
            %   calculates the screen position to center the figure with
            %   a small offset for multiple figures to prevent complete overlap.
            %
            % Inputs:
            %   obj - The FigureManager object
            %   figureIndex - Index of figure in creation order
            %
            % Outputs:
            %   position - [left, bottom, width, height] position vector

            figWidth = obj.figureSize(1);
            figHeight = obj.figureSize(2);
            
            % Calculate centered position
            screenWidth = obj.screenSize(1);
            screenHeight = obj.screenSize(2);
            
            % Center the first figure
            centerX = (screenWidth - figWidth) / 2;
            centerY = (screenHeight - figHeight) / 2;
            
            % Apply offset for subsequent figures to prevent complete overlap
            offset = (figureIndex - 1) * obj.FIGURE_OFFSET;
            left = centerX + offset;
            bottom = centerY - offset;
            
            % Ensure figure stays within screen bounds
            maxLeft = screenWidth - figWidth;
            maxBottom = screenHeight - figHeight;
            left = min(max(0, left), maxLeft);
            bottom = min(max(0, bottom), maxBottom);
            
            position = [left, bottom, figWidth, figHeight];
        end

        function layoutDims = calculateLayout(obj, nFigures)
            % CALCULATELAYOUT Determine optimal layout dimensions.
            %
            %   layoutDims = calculateLayout(obj, nFigures) calculates
            %   the optimal number of rows and columns for arranging
            %   figures without overlap.
            %
            % Inputs:
            %   obj - The FigureManager object
            %   nFigures - Number of figures to arrange
            %
            % Outputs:
            %   layoutDims - [nRows, nCols] layout dimensions

            if nFigures == 1
                layoutDims = [1, 1];
            elseif nFigures == 2
                layoutDims = [1, 2]; % Side by side
            elseif nFigures <= 4
                layoutDims = [2, 2]; % 2x2 grid
            elseif nFigures <= 6
                layoutDims = [2, 3]; % 2x3 grid
            elseif nFigures <= 9
                layoutDims = [3, 3]; % 3x3 grid
            else
                % For larger numbers, approximate square layout
                nCols = ceil(sqrt(nFigures));
                nRows = ceil(nFigures / nCols);
                layoutDims = [nRows, nCols];
            end
        end

        function position = calculateNonOverlappingPosition(obj, figureIndex, layoutDims)
            % CALCULATENONOVERLAPPINGPOSITION Calculate position without overlap.
            %
            %   position = calculateNonOverlappingPosition(obj, figureIndex, layoutDims)
            %   calculates the screen position for the figure to ensure no
            %   overlap with other figures in the layout.
            %
            % Inputs:
            %   obj - The FigureManager object
            %   figureIndex - Index of figure in creation order
            %   layoutDims - [nRows, nCols] layout dimensions
            %
            % Outputs:
            %   position - [left, bottom, width, height] position vector

            nRows = layoutDims(1);
            nCols = layoutDims(2);
            
            % Convert linear index to grid coordinates
            row = floor((figureIndex - 1) / nCols) + 1;
            col = mod(figureIndex - 1, nCols) + 1;
            
            % Use default figure size
            figWidth = obj.figureSize(1);
            figHeight = obj.figureSize(2);
            
            % Calculate total grid dimensions
            totalGridWidth = nCols * figWidth + (nCols - 1) * obj.FIGURE_SPACING;
            totalGridHeight = nRows * figHeight + (nRows - 1) * obj.FIGURE_SPACING;
            
            % Check if grid fits on screen, if not, scale down
            availableWidth = obj.screenSize(1) - 2 * obj.SCREEN_MARGIN;
            availableHeight = obj.screenSize(2) - 2 * obj.SCREEN_MARGIN;
            
            if totalGridWidth > availableWidth || totalGridHeight > availableHeight
                % Scale down figure size to fit
                scaleWidth = availableWidth / totalGridWidth;
                scaleHeight = availableHeight / totalGridHeight;
                scale = min(scaleWidth, scaleHeight);
                
                figWidth = floor(figWidth * scale);
                figHeight = floor(figHeight * scale);
                
                % Recalculate grid dimensions with scaled figures
                totalGridWidth = nCols * figWidth + (nCols - 1) * obj.FIGURE_SPACING;
                totalGridHeight = nRows * figHeight + (nRows - 1) * obj.FIGURE_SPACING;
            end
            
            % Center the entire grid on screen
            gridStartX = (obj.screenSize(1) - totalGridWidth) / 2;
            gridStartY = (obj.screenSize(2) - totalGridHeight) / 2;
            
            % Calculate figure position within grid
            left = gridStartX + (col - 1) * (figWidth + obj.FIGURE_SPACING);
            bottom = gridStartY + (nRows - row) * (figHeight + obj.FIGURE_SPACING);
            
            position = [left, bottom, figWidth, figHeight];
        end
    end
end