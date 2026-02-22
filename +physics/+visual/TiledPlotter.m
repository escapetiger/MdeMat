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
    % See also:
    %   physics.visual.Strategy, physics.visual.Visualizer

    properties (SetAccess = protected, GetAccess = public)
        Size % [nRows, nCols] overall grid dimensions
        FigObj % Figure handle for this plotter
        Axes % Array of axes handles for rendering
        Specs % Array of axis specifications with rows/cols fields
        TiledObj % Handle to the tiledlayout container
        TileIndices % Computed tile indices for axis placement
        VarToAxis % Structure mapping variable names to axis indices
        Strategy % Rendering strategy object
    end

    methods
        function obj = TiledPlotter(nRows, nCols)
            % TILEDPLOTTER Construct an instance of TiledPlotter.
            %
            %   obj = TiledPlotter(nRows, nCols) creates a tiled plotter with
            %   the specified grid dimensions. A figure must be assigned using
            %   setFigure() before the tiledlayout can be created with draft().

            arguments
                nRows{mustBeNonnegative, mustBeInteger}
                nCols{mustBeNonnegative, mustBeInteger}
            end

            obj.Size = [max(1, nRows), max(1, nCols)];
            obj.Specs = struct('rows', {}, 'cols', {});
            obj.FigObj = [];
            obj.Axes = [];
            obj.TileIndices = [];
            obj.Strategy = [];
            obj.VarToAxis = struct();
        end

        function obj = setStrategy(obj, strategy)
            % SETSTRATEGY Configure the rendering strategy.
            %
            %   obj = setStrategy(obj, strategy) sets the strategy object that
            %   will be used to render data onto the axes. The strategy
            %   determines how data is mapped to visual elements and coordinates
            %   the rendering process.

            arguments
                obj physics.visual.TiledPlotter
                strategy physics.visual.Strategy
            end

            obj.Strategy = strategy;
        end

        function obj = addAxisSpec(obj, nRows, nCols)
            % ADDAXISSPEC Add axis specification for a panel.
            %
            %   obj = addAxisSpec(obj, nRows, nCols) adds a specification
            %   for an axis that will span the specified number of rows and
            %   columns in the tiled layout. Axes are placed automatically
            %   in the order they are added.

            arguments
                obj physics.visual.TiledPlotter
                nRows{mustBePositive, mustBeInteger}
                nCols{mustBePositive, mustBeInteger}
            end

            nTotalRows = obj.Size(1);
            nTotalCols = obj.Size(2);
            core.except.assert(nRows <= nTotalRows && nCols <= nTotalCols, ...
                'InvalidInput', 'Axis span exceeds grid size.');

            obj.Specs(end+1) = struct('rows', nRows, 'cols', nCols);
        end

        function obj = setFigure(obj, figureHandle)
            % SETFIGURE Assign a figure handle to this plotter.
            %
            %   obj = setFigure(obj, figureHandle) assigns the specified
            %   figure handle to this plotter. The figure will be used
            %   when draft() is called to create the tiledlayout.

            arguments
                obj physics.visual.TiledPlotter
                figureHandle matlab.ui.Figure
            end

            obj.FigObj = figureHandle;
        end

        function obj = draft(obj, options)
            % DRAFT Create the tiledlayout and axes.
            %
            %   obj = draft(obj) creates the MATLAB tiledlayout object and
            %   generates all axes according to the added specifications.
            %   A figure must be assigned using setFigure() before calling
            %   this method.
            %
            %   obj = draft(obj, 'ParameterName', ParameterValue) passes
            %   additional arguments to the tiledlayout constructor.

            arguments
                obj physics.visual.TiledPlotter
                options struct = struct()
            end

            if ~isempty(obj.Axes), return; end

            core.except.assert(~isempty(obj.FigObj) && isvalid(obj.FigObj), ...
                'InvalidFigure', 'Figure must be valid and assigned before drafting.');

            %< Create tiledlayout in the assigned figure
            nRows = obj.Size(1);
            nCols = obj.Size(2);
            figure(obj.FigObj);
            args = namedargs2cell(options);
            obj.TiledObj = tiledlayout(nRows, nCols, args{:});
            obj.computeTileIndices();

            %< Create axes according to specifications
            nAxes = length(obj.Specs);
            obj.Axes = gobjects(1, nAxes);
            for iAxis = 1:nAxes
                spec = obj.Specs(iAxis);
                tileIdx = obj.TileIndices(iAxis);
                ax = nexttile(obj.TiledObj, tileIdx, [spec.rows, spec.cols]);
                obj.Axes(iAxis) = ax;
            end
        end

        function obj = bind(obj, variable, axisIdx)
            % BIND Map a variable to a specific axis index.
            %
            %   obj = bind(obj, variable, axisIdx) creates a mapping from
            %   the specified variable name to an axis index for use.

            arguments
                obj physics.visual.TiledPlotter
                variable{mustBeTextScalar}
                axisIdx{mustBePositive, mustBeInteger}
            end

            obj.VarToAxis.(variable) = axisIdx;
        end

        function obj = render(obj, database, titleText, options)
            % RENDER Render data using the configured strategy.
            %
            %   obj = render(obj, database, titleText) renders the provided
            %   database using the configured strategy with default styling.
            %
            %   obj = render(obj, database, titleText, renderFn=renderFn)
            %   uses a custom rendering function instead of the strategy.

            arguments
                obj physics.visual.TiledPlotter
                database physics.visual.Database
                titleText{mustBeTextScalar} = ''
                options.renderFn{mustBeFunctionOrEmpty} = []
            end

            core.except.assert(~isempty(obj.Axes), ...
                'EmptyAxes', 'Axes must be drafted first.');
            core.except.assert(~isempty(obj.FigObj) && isvalid(obj.FigObj), ...
                'InvalidFigure', 'Figure must be valid and drafted first.');

            figure(obj.FigObj);

            if ~isempty(options.renderFn)
                options.renderFn(obj.Axes, database);
            else
                core.except.assert(~isempty(obj.Strategy), ...
                    'EmptyStrategy', 'Strategy must be set first.');

                variables = fieldnames(obj.VarToAxis);
                for j = 1:length(variables)
                    variable = variables{j};
                    offset = obj.VarToAxis.(variable);

                    if offset > length(obj.Axes)
                        core.except.verify(false, 'InvalidAxisIdx', ...
                            'Axis index %d exceeds available axes.', offset);
                        continue;
                    end

                    group = database.groupBy(variable);

                    switch obj.Strategy.NDims
                        case 1
                            ax = obj.Axes(offset);
                            obj.Strategy.render(ax, group, variable);
                        case {2, 3}
                            axes = obj.Axes(offset+(0:group.NDatasets - 1));
                            obj.Strategy.render(axes, group, variable);
                    end
                end
            end

            if ~isempty(titleText)
                switch obj.Strategy.NDims
                    case 1
                        obj.addSgTitle(titleText);
                    case {2, 3}
                        obj.addAnnotationTitle(titleText);
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

            obj.Axes = [];
            obj.TiledObj = [];
            obj.TileIndices = [];
        end
    end

    methods (Access = protected)
        function obj = computeTileIndices(obj)
            % COMPUTETILEINDICES Compute tile placement for axis
            % specifications.

            rows = obj.Size(1);
            cols = obj.Size(2);
            occupied = false(rows, cols);
            nSpecs = length(obj.Specs);
            obj.TileIndices = zeros(1, nSpecs);

            for iSpec = 1:nSpecs
                spec = obj.Specs(iSpec);
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
                            obj.TileIndices(iSpec) = idx;
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

            sgtitle(obj.TiledObj, titleText, 'FontSize', 14, 'FontWeight', 'bold');
        end

        function addAnnotationTitle(obj, titleText)
            % ADDANNOTATIONTITLE Add title using annotation for reliable positioning.

            delete(findall(obj.FigObj, 'Tag', 'AnnotationTitle'));

            annotation(obj.FigObj, 'textbox', [0, 0.92, 1, 0.08], ...
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

function mustBeFunctionOrEmpty(x)
if ~isempty(x)
    mustBeA(x, 'function_handle');
end
end