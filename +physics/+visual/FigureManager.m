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
        FigureHandles % Cell array of created figure handles
        ScreenSize % Screen dimensions [width, height] in pixels
        FigureSize % Standard figure size [width, height] in pixels
    end
    
    properties (Constant)
        DefaultFigureWidth = 560 % Default MATLAB figure width in pixels
        DefaultFigureHeight = 420 % Default MATLAB figure height in pixels
        FigureSpacing = 20 % Spacing between figures in pixels
        ScreenMargin = 50 % Margin from screen edges in pixels
    end
    
    properties (Dependent)
        NFigures % Number of managed figures
    end
    
    methods
        function obj = FigureManager()
            % FIGUREMANAGER Construct an instance of FigureManager.
            %
            %   obj = FigureManager() creates a figure manager and initializes
            %   screen size information for automatic positioning calculations.
            
            obj.FigureHandles = {};
            obj.ScreenSize = obj.getScreenSize();
            obj.FigureSize = [obj.DefaultFigureWidth, obj.DefaultFigureHeight];
        end
        
        function n = get.NFigures(obj)
            % GET.NFIGURES Get the number of managed figures.
            
            n = length(obj.FigureHandles);
        end
        
        function figureHandles = createFigures(obj, nFigures)
            % CREATEFIGURES Create non-overlapping figures with default styling.
            %
            %   figureHandles = createFigures(obj, nFigures) creates the
            %   specified number of figures positioned to avoid overlap. For
            %   1-2 figures, they are placed side by side. For more figures,
            %   they are arranged in a simple grid pattern.
            
            arguments
                obj physics.visual.FigureManager
                nFigures double {mustBePositive, mustBeInteger}
            end
            
            layoutDims = obj.calculateLayout(nFigures);
            
            obj.FigureHandles = cell(1, nFigures);
            for iFigure = 1:nFigures
                position = obj.calculateNonOverlappingPosition(iFigure, layoutDims);
                
                figHandle = figure(Position=position);
                
                obj.FigureHandles{iFigure} = figHandle;
            end
            
            figureHandles = obj.FigureHandles;
        end
        
        function obj = closeAll(obj)
            % CLOSEALL Close all managed figures.
            %
            %   obj = closeAll(obj) closes all figures created by this
            %   manager and clears the figure handle cache.
            
            for i = 1:length(obj.FigureHandles)
                figHandle = obj.FigureHandles{i};
                if isvalid(figHandle)
                    close(figHandle);
                end
            end
            obj.FigureHandles = {};
        end
        
        function TF = isValid(obj)
            % ISVALID Check if all managed figures are still valid.
            %
            %   TF = isValid(obj) returns true if all managed figures
            %   are still open and valid.
            
            TF = true;
            for i = 1:length(obj.FigureHandles)
                if ~isvalid(obj.FigureHandles{i})
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
            
            arguments
                obj physics.visual.FigureManager
                figureIndex double {mustBePositive, mustBeInteger}
            end
            
            core.except.assert(figureIndex <= obj.NFigures, ...
                'InvalidIndex', 'Figure index must be between 1 and %d', obj.NFigures);
            
            figHandle = obj.FigureHandles{figureIndex};
            if isvalid(figHandle)
                figure(figHandle);
            end
        end
        
        function obj = bringAllToFront(obj)
            % BRINGALLTOFRONT Bring all managed figures to front.
            %
            %   obj = bringAllToFront(obj) brings all managed figures
            %   to the front in creation order.
            
            for i = 1:length(obj.FigureHandles)
                obj.bringToFront(i);
            end
        end
    end
    
    methods (Access = protected)
        function position = calculateNonOverlappingPosition(obj, figureIndex, layoutDims)
            % CALCULATENONOVERLAPPINGPOSITION Calculate position without overlap.
            %
            %   position = calculateNonOverlappingPosition(obj, figureIndex, layoutDims)
            %   calculates the screen position for the figure to ensure no
            %   overlap with other figures in the layout.
            
            nRows = layoutDims(1);
            nCols = layoutDims(2);
            
            % Convert linear index to grid coordinates
            row = floor((figureIndex - 1) / nCols) + 1;
            col = mod(figureIndex - 1, nCols) + 1;
            
            % Use default figure size
            figWidth = obj.FigureSize(1);
            figHeight = obj.FigureSize(2);
            
            % Calculate total grid dimensions
            totalGridWidth = nCols * figWidth + (nCols - 1) * obj.FigureSpacing;
            totalGridHeight = nRows * figHeight + (nRows - 1) * obj.FigureSpacing;
            
            % Check if grid fits on screen, if not, scale down
            availableWidth = obj.ScreenSize(1) - 2 * obj.ScreenMargin;
            availableHeight = obj.ScreenSize(2) - 2 * obj.ScreenMargin;
            
            if totalGridWidth > availableWidth || totalGridHeight > availableHeight
                % Scale down figure size to fit
                scaleWidth = availableWidth / totalGridWidth;
                scaleHeight = availableHeight / totalGridHeight;
                scale = min(scaleWidth, scaleHeight);
                
                figWidth = floor(figWidth * scale);
                figHeight = floor(figHeight * scale);
                
                % Recalculate grid dimensions with scaled figures
                totalGridWidth = nCols * figWidth + (nCols - 1) * obj.FigureSpacing;
                totalGridHeight = nRows * figHeight + (nRows - 1) * obj.FigureSpacing;
            end
            
            % Center the entire grid on screen
            gridStartX = (obj.ScreenSize(1) - totalGridWidth) / 2;
            gridStartY = (obj.ScreenSize(2) - totalGridHeight) / 2;
            
            % Calculate figure position within grid
            left = gridStartX + (col - 1) * (figWidth + obj.FigureSpacing);
            bottom = gridStartY + (nRows - row) * (figHeight + obj.FigureSpacing);
            
            position = [left, bottom, figWidth, figHeight];
        end
    end
    
    methods (Static)
        function screenSize = getScreenSize()
            % GETSCREENSIZE Get usable screen dimensions.
            %
            %   screenSize = getScreenSize() returns the screen
            %   dimensions [width, height] in pixels for positioning.
            
            screenSizeVec = get(0, 'ScreenSize'); % [left, bottom, width, height]
            screenSize = [screenSizeVec(3), screenSizeVec(4)];
        end
        
        function layoutDims = calculateLayout(nFigures)
            % CALCULATELAYOUT Determine optimal layout dimensions.
            %
            %   layoutDims = calculateLayout(nFigures) calculates
            %   the optimal number of rows and columns for arranging
            %   figures without overlap.
            
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
    end
end