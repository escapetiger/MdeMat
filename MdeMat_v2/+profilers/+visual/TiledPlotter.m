classdef TiledPlotter < handle
    % TILEDPLOTTER Multi-panel figure manager for data visualization.
    %
    %   TiledPlotter provides infrastructure for creating and managing
    %   multi-panel figures using MATLAB's tiledlayout functionality.
    %   This class simplifies the creation of complex visualization layouts
    %   with automatic tile placement and flexible axis specifications.
    %
    %   The plotter supports both manual axis specification and automatic
    %   layout computation, making it suitable for both simple multi-plot
    %   figures and complex dashboard-style visualizations with varying
    %   panel sizes.
    %
    % Examples:
    %   % Create 2×3 tiled layout
    %   plotter = TiledPlotter(2, 3);
    %   
    %   % Add axis specifications for different panel sizes
    %   plotter.addAxisSpec(1, 1); % Standard 1×1 panel
    %   plotter.addAxisSpec(1, 2); % Wide 1×2 panel
    %   plotter.addAxisSpec(2, 1); % Tall 2×1 panel
    %   
    %   % Set up layout and strategy
    %   plotter.draft('TileSpacing', 'compact');
    %   plotter.setStrategy(strategy);
    %   
    %   % Render data
    %   plotter.render('My Results', dataStructure, styleStructure);
    %
    % See also:
    %   profilers.visual.Strategy, 
    %   profilers.visual.GeometryPlotter

    properties (SetAccess = protected, GetAccess = public)
        size        % [nRows, nCols] overall grid dimensions
        axes        % Array of axes handles for rendering
        specs       % Array of axis specifications with rows/cols fields
        tiledObj    % Handle to the tiledlayout container
        tileIndices % Computed tile indices for axis placement
        strategy    % Rendering strategy object
    end

    methods
        function obj = TiledPlotter(nRows, nCols)
            % TILEDPLOTTER Constructor for TiledPlotter.
            %
            %   obj = TiledPlotter(nRows, nCols) creates a tiled plotter
            %   with the specified grid dimensions. The actual tiledlayout
            %   is created when draft() is called.
            %
            % Inputs:
            %   nRows - Number of rows in the overall grid (positive integer)
            %   nCols - Number of columns in the overall grid (positive integer)
            %
            % Outputs:
            %   obj - Constructed TiledPlotter object

            obj.size = [max(1, nRows), max(1, nCols)];
            obj.specs = struct('rows', {}, 'cols', {});
            obj.axes = [];
            obj.tileIndices = [];
            obj.strategy = [];
        end

        function obj = setStrategy(obj, strategy)
            % SETSTRATEGY Configure the rendering strategy.
            %
            %   obj = setStrategy(obj, strategy) sets the strategy object
            %   that will be used to render data onto the axes. The
            %   strategy determines how data is mapped to visual elements
            %   and coordinates the rendering process.
            %
            % Inputs:
            %   obj - The TiledPlotter object
            %   strategy - Strategy object implementing the rendering logic
            %
            % Outputs:
            %   obj - The TiledPlotter object

            obj.strategy = strategy;
        end

        function obj = addAxisSpec(obj, nRows, nCols)
            % ADDAXISSPEC Add axis specification for a panel.
            %
            %   obj = addAxisSpec(obj, nRows, nCols) adds a specification
            %   for an axis that will span the specified number of rows and
            %   columns in the tiled layout. Axes are placed automatically
            %   in the order they are added.
            %
            % Inputs:
            %   obj - The TiledPlotter object
            %   nRows - Number of rows this axis should span (positive integer)
            %   nCols - Number of columns this axis should span (positive integer)
            %
            % Outputs:
            %   obj - The TiledPlotter object

            nTotalRows = obj.size(1);
            nTotalCols = obj.size(2);
            core.except.assert(nRows <= nTotalRows && nCols <= nTotalCols, ...
                'InvalidInput', 'Axis span exceeds grid size.');

            obj.specs(end+1) = struct('rows', nRows, 'cols', nCols);
        end

        function obj = draft(obj, varargin)
            % DRAFT Create the tiledlayout and axes.
            %
            %   obj = draft(obj) creates the MATLAB tiledlayout object and
            %   generates all axes according to the added specifications.
            %   This method must be called before rendering can occur.
            %
            %   obj = draft(obj, 'ParameterName', ParameterValue) passes
            %   additional arguments to the tiledlayout constructor.
            %
            % Inputs:
            %   obj - The TiledPlotter object
            %   varargin - Additional arguments for tiledlayout
            %
            % Outputs:
            %   obj - The TiledPlotter object

            if ~isempty(obj.axes), return; end

            nRows = obj.size(1);
            nCols = obj.size(2);
            obj.tiledObj = tiledlayout(nRows, nCols, varargin{:});
            obj.computeTileIndices();

            nAxes = length(obj.specs);
            obj.axes = gobjects(1, nAxes);
            for iAxis = 1:nAxes
                spec = obj.specs(iAxis);
                tileIdx = obj.tileIndices(iAxis);
                ax = nexttile(obj.tiledObj, tileIdx, [spec.rows, spec.cols]);
                obj.axes(iAxis) = ax;
            end
        end

        function obj = render(obj, titleText, dataset, styleset, renderFn)
            % RENDER Render data using the configured strategy.
            %
            %   obj = render(obj, titleText, dataset) renders the provided
            %   dataset using the configured strategy with default styling.
            %
            %   obj = render(obj, titleText, dataset, styleset) renders
            %   with custom styling specifications.
            %
            %   obj = render(obj, titleText, dataset, styleset, renderFn)
            %   uses a custom rendering function instead of the strategy.
            %
            % Inputs:
            %   obj - The TiledPlotter object
            %   titleText - Overall title for the figure (string)
            %   dataset - Structure containing data to render
            %   styleset - Rendering style specifications (optional)
            %   renderFn - Custom rendering function handle (optional)
            %
            % Outputs:
            %   obj - The TiledPlotter object

            core.except.assert(~isempty(obj.axes), ...
                'EmptyAxes', 'Axes must be drafted first.');

            %< Prepare default styleset if not provided
            if nargin < 4 || isempty(styleset)
                datasetNames = fieldnames(dataset);
                styleset = struct();
                for i = 1:length(datasetNames)
                    styleset.(datasetNames{i}) = {};
                end
            end

            %< Choose rendering method
            if nargin < 5 || isempty(renderFn)
                core.except.assert(~isempty(obj.strategy), ...
                    'EmptyStrategy', 'Strategy must be set first.');

                obj.strategy.render(obj.axes, dataset, styleset);
            else
                renderFn(obj.axes, dataset, styleset);
            end

            %< Set overall title
            if ~isempty(titleText)
                sgtitle(titleText);
            end

            drawnow;
        end
    end

    methods (Access = protected)
        function obj = computeTileIndices(obj)
            % COMPUTETILEINDICES Compute tile placement for axis
            % specifications.
            %
            %   computeTileIndices(obj) automatically determines where to
            %   place each axis specification in the tiled layout, ensuring
            %   no overlaps occur. Uses a greedy placement algorithm.
            % 
            % Inputs:
            %   obj - The TiledPlotter object
            %
            % Outputs:
            %   obj - The TiledPlotter object

            rows = obj.size(1);
            cols = obj.size(2);
            occupied = false(rows, cols);
            nSpecs = length(obj.specs);
            obj.tileIndices = zeros(1, nSpecs);

            for iSpec = 1:nSpecs
                spec = obj.specs(iSpec);
                placed = false;

                %< Try to place the current specification
                for idx = 1:(rows * cols)
                    row = floor((idx - 1) / cols) + 1;
                    col = mod(idx - 1, cols) + 1;

                    %< Check if specification fits at this position
                    if row + spec.rows - 1 <= rows && col + spec.cols - 1 <= cols
                        block = occupied(row:row+spec.rows-1, col:col+spec.cols-1);
                        if all(block(:) == false)
                            occupied(row:row+spec.rows-1, col:col+spec.cols-1) = true;
                            obj.tileIndices(iSpec) = idx;
                            placed = true;
                            break;
                        end
                    end
                end

                core.except.assert(placed, 'BadAxisSpec', ...
                    'Cannot place axis spec %d in the grid.', iSpec);
            end
        end
    end
end