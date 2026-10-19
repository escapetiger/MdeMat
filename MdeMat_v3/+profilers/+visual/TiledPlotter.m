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
    %   Each plotter instance receives its figure from a parent Visualizer,
    %   providing centralized figure management and better control over
    %   multiple visualization windows.
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
    %   % Assign figure and set up layout (managed by Visualizer)
    %   fig = figure('Name', 'Results', 'NumberTitle', 'off');
    %   plotter.setFigure(fig);
    %   plotter.draft('TileSpacing', 'compact');
    %   plotter.setStrategy(strategy);
    %
    %   % Render data
    %   plotter.render('My Results', dataStructure, styleStructure);
    %
    % See also:
    %   profilers.visual.Strategy, profilers.visual.GeometryPlotter

    properties (Constant)
        INCOMPATIABLE_STRATEGIES = {'Strategy2d', 'Strategy3d'}
    end

    properties (SetAccess = protected, GetAccess = public)
        size % [nRows, nCols] overall grid dimensions
        figObj % Figure handle for this plotter
        axes % Array of axes handles for rendering
        specs % Array of axis specifications with rows/cols fields
        tiledObj % Handle to the tiledlayout container
        tileIndices % Computed tile indices for axis placement
        varToAxis % Structure mapping varName names to axis indices
        strategy % Rendering strategy object
    end

    methods
        function obj = TiledPlotter(nRows, nCols)
            % TILEDPLOTTER Constructor for TiledPlotter.
            %
            %   obj = TiledPlotter(nRows, nCols) creates a tiled plotter
            %   with the specified grid dimensions. A figure must be
            %   assigned using setFigure() before the tiledlayout can be
            %   created with draft().
            %
            % Inputs:
            %   nRows - Number of rows in the overall grid (positive integer)
            %   nCols - Number of columns in the overall grid (positive integer)
            %
            % Outputs:
            %   obj - Constructed TiledPlotter object

            obj.size = [max(1, nRows), max(1, nCols)];
            obj.specs = struct('rows', {}, 'cols', {});
            obj.figObj = [];
            obj.axes = [];
            obj.tileIndices = [];
            obj.strategy = [];
            obj.varToAxis = struct();
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

        function obj = setFigure(obj, figureHandle)
            % SETFIGURE Assign a figure handle to this plotter.
            %
            %   obj = setFigure(obj, figureHandle) assigns the specified
            %   figure handle to this plotter. The figure will be used
            %   when draft() is called to create the tiledlayout.
            %
            % Inputs:
            %   obj - The TiledPlotter object
            %   figureHandle - Valid MATLAB figure handle
            %
            % Outputs:
            %   obj - The TiledPlotter object

            core.except.assert(isa(figureHandle, 'matlab.ui.Figure') && isvalid(figureHandle), ...
                'InvalidInput', 'Input must be a valid figure handle.');
            obj.figObj = figureHandle;
        end

        function obj = draft(obj, varargin)
            % DRAFT Create the tiledlayout and axes.
            %
            %   obj = draft(obj) creates the MATLAB tiledlayout object and
            %   generates all axes according to the added specifications.
            %   A figure must be assigned using setFigure() before calling
            %   this method.
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
            %
            % Examples:
            %   fig = figure('Name', 'Results');
            %   plotter.setFigure(fig);
            %   plotter.draft('TileSpacing', 'compact');

            if ~isempty(obj.axes), return; end

            core.except.assert(~isempty(obj.figObj) && isvalid(obj.figObj), ...
                'InvalidFigure', 'A valid figure must be assigned before drafting.');

            %< Create tiledlayout in the assigned figure
            nRows = obj.size(1);
            nCols = obj.size(2);
            figure(obj.figObj);
            obj.tiledObj = tiledlayout(nRows, nCols, varargin{:});
            obj.computeTileIndices();

            %< Create axes according to specifications
            nAxes = length(obj.specs);
            obj.axes = gobjects(1, nAxes);
            for iAxis = 1:nAxes
                spec = obj.specs(iAxis);
                tileIdx = obj.tileIndices(iAxis);
                ax = nexttile(obj.tiledObj, tileIdx, [spec.rows, spec.cols]);
                obj.axes(iAxis) = ax;
            end
        end

        function obj = bind(obj, varName, axisIdx)
            % BIND Map a varName to a specific axis index.
            %
            %   obj = bind(obj, varName, axisIdx) creates a mapping from
            %   the specified varName name to an axis index for use.
            %
            % Inputs:
            %   obj - The TiledPlotter object
            %   varName - Name of the varName to map (string)
            %   axisIdx - Target axis index (positive integer)
            %
            % Outputs:
            %   obj - The TiledPlotter object

            obj.varToAxis.(varName) = axisIdx;
        end

        function obj = render(obj, titleText, database, renderFn)
            % RENDER Render data using the configured strategy.
            %
            %   obj = render(obj, titleText, database) renders the provided
            %   database using the configured strategy with default styling.
            %
            %   obj = render(obj, titleText, database, renderFn) uses a
            %   custom rendering function instead of the strategy.
            %
            % Inputs:
            %   obj - The TiledPlotter object
            %   titleText - Overall title for the figure (string)
            %   database - Database object containing data to render
            %   renderFn - Custom rendering function handle (optional)
            %
            % Outputs:
            %   obj - The TiledPlotter object

            core.except.assert(~isempty(obj.axes), ...
                'EmptyAxes', 'Axes must be drafted first.');
            core.except.assert(~isempty(obj.figObj) && isvalid(obj.figObj), ...
                'InvalidFigure', 'Figure must be valid and drafted first.');

            %< Make sure the plotter's figure is current
            figure(obj.figObj);

            if nargin < 4 || isempty(renderFn)
                core.except.assert(~isempty(obj.strategy), ...
                    'EmptyStrategy', 'Strategy must be set first.');

                varNames = fieldnames(obj.varToAxis);
                for j = 1:length(varNames)
                    varName = varNames{j};
                    offset = obj.varToAxis.(varName);

                    if offset > length(obj.axes)
                        core.except.verify(false, 'InvalidAxisIdx', ...
                            'Axis index %d exceeds available axes.', offset);
                        continue;
                    end

                    varDataset = database.parse(varName);

                    obj.strategy.render(obj.axes, varDataset, varName, offset);
                end
            else
                renderFn(obj.axes, database);
            end

            if ~isempty(titleText)
                prefix = 'profilers.visual.';
                if any(cellfun(@(s) isa(obj.strategy, [prefix, s]), obj.INCOMPATIABLE_STRATEGIES))
                    obj.addAnnotationTitle(titleText);
                else
                    obj.addSgTitle(titleText);
                end
            end

            drawnow;
        end
        
        function obj = close(obj)
            % CLOSE Clear the plotter's graphics objects.
            %
            %   obj = close(obj) clears the tiledlayout and axes objects
            %   but does not close the figure (figure management is handled
            %   by the parent Visualizer). After calling this method,
            %   draft() must be called again before rendering.
            %
            % Inputs:
            %   obj - The TiledPlotter object
            %
            % Outputs:
            %   obj - The TiledPlotter object

            obj.axes = [];
            obj.tiledObj = [];
            obj.tileIndices = [];
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
                    row = floor((idx - 1)/cols) + 1;
                    col = mod(idx-1, cols) + 1;

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

        function addSgTitle(obj, titleText)
            % ADDSGTITLE Add built-in sgtitle.
            %
            % Inputs:
            %   obj - The TiledPlotter object
            %   titleText - Title text to display

            sgtitle(obj.tiledObj, titleText, ...
                'FontSize', 14, ...
                'FontWeight', 'bold');
        end

        function addAnnotationTitle(obj, titleText)
            % ADDANNOTATIONTITLE Add title using annotation for reliable positioning.
            %
            % Inputs:
            %   obj - The TiledPlotter object
            %   titleText - Title text to display

            % Remove any existing annotation title
            delete(findall(obj.figObj, 'Tag', 'AnnotationTitle'));

            % Create title using annotation (figure-normalized coordinates)
            annotation(obj.figObj, 'textbox', [0, 0.92, 1, 0.08], ...
                'String', titleText, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'EdgeColor', 'none', ...
                'BackgroundColor', 'none', ...
                'FontSize', 14, ...
                'FontWeight', 'bold', ...
                'Tag', 'AnnotationTitle', ...
                'FitBoxToText', 'off');
        end
    end
end