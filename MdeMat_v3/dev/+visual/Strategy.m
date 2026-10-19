classdef Strategy < handle
    % STRATEGY Enhanced base class for data rendering strategies.
    %
    %   Strategy defines the interface and common functionality for all
    %   rendering strategies used in data visualization. This class
    %   provides mode-based rendering control, variable-to-axis mapping
    %   capabilities, default styling, and utility methods for coordinated 
    %   multi-panel visualizations.
    %
    %   The class includes merged helper methods for data preparation,
    %   rendering coordination, and style management, following a unified
    %   approach to visualization workflows.
    %
    % See also:
    %   profilers.visual.Strategy1d, profilers.visual.Strategy2d,
    %   profilers.visual.Strategy3d, profilers.visual.TiledPlotter

    properties (Constant)
        % Color palettes for different scenarios
        DEFAULT_MARKERS = {'o', 's', '^', 'v', 'd', 'x', '+'};
        DEFAULT_LINES = {'-', '--', ':', '-.', '-', '--', ':'};
        DEFAULT_COLOR_PALETTE = {'#000000', '#004488', '#BB5566', '#DDAA33', '#FFFFFF'};
%         DEFAULT_COLOR_PALETTE = {'#000000', '#AA3377', '#EE6677', '#CCBB44', '#228833', '#66CCEE', '#4477AA'};
        
        % Line and plot styling
        DEFAULT_LINE_WIDTH = 1;
        DEFAULT_COLORMAP_2D = 'turbo';
        DEFAULT_COLORMAP_3D = 'turbo';
        
        % Color limit and padding settings
        DEFAULT_COLOR_PADDING = 0.05;  % 5% padding for color limits
        USE_SYMMETRIC_COLOR_LIMITS = true;
        
        % Text and legend formatting
        DEFAULT_LEGEND_POSITION = 'best';
        DEFAULT_AXIS_LABEL_FONT_SIZE = 16;
        DEFAULT_LEGEND_FONT_SIZE = 14;
        REPLACE_UNDERSCORE_WITH_HYPHEN = true;
        
        % 1D plot options
        ENABLE_LOG_PLOT = true;
    end

    properties
        nDims   % Number of spatial dimensions handled by this strategy
        varMap  % Structure mapping variable names to axis indices
        mode    % Rendering mode ('sequential' | 'map')
        logScale % Structure controlling log scale for 1D plots (x, y)
    end

    methods
        function obj = Strategy(nDims, mode)
            % STRATEGY Constructor for Strategy base class.
            %
            %   obj = Strategy(nDims, mode) creates a rendering strategy
            %   for the specified number of dimensions with the given
            %   rendering mode.
            %
            % Inputs:
            %   nDims - Number of spatial dimensions (1, 2, or 3)
            %   mode - Rendering mode ('sequential' | 'map', default: 'sequential')
            %
            % Notes:
            %   mode = 'sequential': Variables rendered in order of appearance
            %   mode = 'map': Variables rendered according to bind() mapping
            %
            % Outputs:
            %   obj - Constructed Strategy object

            if nargin < 2, mode = 'sequential'; end

            core.except.assert(ismember(mode, {'sequential', 'map'}), ...
                'InvalidInput', 'Mode must be ''sequential'' or ''map''');

            obj.nDims = nDims;
            obj.varMap = struct();
            obj.mode = mode;
            obj.logScale = struct('x', false, 'y', false);
        end

        function obj = bind(obj, varName, axisIdx)
            % BIND Map a variable to a specific axis index.
            %
            %   obj = bind(obj, varName, axisIdx) creates a mapping from
            %   the specified variable name to an axis index for use in
            %   'map' rendering mode.
            %
            % Inputs:
            %   obj - The Strategy object
            %   varName - Name of the variable to map (string)
            %   axisIdx - Target axis index (positive integer)
            % 
            % Outputs:
            %   obj - The Strategy object

            obj.varMap.(varName) = axisIdx;
        end

        function obj = setLogScale(obj, axis, isEnabled)
            % SETLOGSCALE Configure logarithmic scaling for 1D plots.
            %
            %   obj = setLogScale(obj, axis, isEnabled) enables or disables
            %   logarithmic scaling for the specified axis in 1D plots.
            %
            % Inputs:
            %   obj - The Strategy object
            %   axis - Axis to configure ('x' or 'y')
            %   isEnabled - Boolean flag for log scale (logical)
            %
            % Outputs:
            %   obj - The Strategy object

            core.except.assert(ismember(axis, {'x', 'y'}), ...
                'InvalidInput', 'Axis must be ''x'' or ''y''');
            obj.logScale.(axis) = isEnabled;
        end

        function obj = render(obj, axes, dataset)
            % RENDER Render data with variable-to-axis coordination.
            %
            %   obj = render(obj, axes, dataset) renders the provided
            %   dataset onto the axes array using either sequential or
            %   mapped variable placement with automatic default styling.
            %
            % Inputs:
            %   obj - The Strategy object
            %   axes - Array of axes handles for rendering
            %   dataset - Structure containing data fields to render
            % 
            % Outputs:
            %   obj - The Strategy object

            if strcmp(obj.mode, 'sequential')
                varNames = obj.parseVarNames(dataset);
                obj.renderSequential(axes, dataset, varNames);
            else
                varNames = fieldnames(obj.varMap);
                obj.renderMapped(axes, dataset, varNames);
            end
        end

        function varNames = parseVarNames(obj, dataset)
            % PARSEVARNAMES Collect all variable names from datasets.
            %
            %   varNames = parseVarNames(obj, dataset) extracts all unique
            %   variable names from the provided dataset structure.
            %
            % Inputs:
            %   obj - The Strategy object
            %   dataset - Map of dataset names to dataset objects
            %
            % Outputs:
            %   varNames - Cell array of all unique variable names

            dataNames = fieldnames(dataset);
            varNames = {};

            for i = 1:length(dataNames)
                dataName = dataNames{i};
                data = dataset.(dataName);
                subVarNames = fieldnames(data.data);

                for j = 1:length(subVarNames)
                    if ~ismember(subVarNames{j}, varNames)
                        varNames{end+1} = subVarNames{j};
                    end
                end
            end
        end

        function subdata = prepareData(obj, dataset, varName)
            % PREPAREDATA Prepare data for rendering from datasets.
            %
            %   subdata = prepareData(obj, dataset, varName) extracts and
            %   formats data for the specified variable from all datasets,
            %   organizing coordinates, values, and data types for visualization.
            %
            % Inputs:
            %   obj - The Strategy object
            %   dataset - Map of dataset names to dataset objects
            %   varName - Name of the variable to prepare
            %
            % Outputs:
            %   subdata - Map of dataset names to prepared data structures
            %             Each structure contains: u, x1, x2, x3, dataType

            dataNames = fieldnames(dataset);
            [fieldName, varIdx] = obj.splitVarName(varName);
            subdata = struct();

            for i = 1:length(dataNames)
                dataName = dataNames{i};
                g = dataset.(dataName);

                if ~isfield(g.data, fieldName)
                    continue;
                end

                if isempty(g.coords), continue; end

                x = g.coords;
                u = g.data.(fieldName)(:, varIdx);

                %< Parse coordinates based on format and dimensions
                if iscell(x)
                    coords = x(1:obj.nDims);
                else
                    coords = cell(1, obj.nDims);
                    for iDim = 1:obj.nDims
                        coords{iDim} = unique(x(iDim, :));
                    end
                end

                %< Reshape data if needed
                if obj.nDims == 1
                    n = [1, length(coords{1})];
                else
                    n = cellfun(@length, coords);
                end
                if length(u) == prod(n)
                    u = reshape(u, n(:).');
                end

                %< Create data structure with type information
                dataStruct = struct('u', u, 'type', g.type);
                for iDim = 1:obj.nDims
                    coordName = sprintf('x%d', iDim);
                    dataStruct.(coordName) = coords{iDim};
                end

                subdata.(dataName) = dataStruct;
            end
        end

        function [fieldName, varIdx] = splitVarName(obj, str)
            % SPLITVARNAME Split variable name into field and index.
            %
            %   [fieldName, varIdx] = splitVarName(obj, str) parses a
            %   variable name string to extract the base field name and
            %   numerical index component.
            %
            % Inputs:
            %   obj - The Strategy object
            %   str - Variable name string (e.g., 'u1', 'velocity2')
            %
            % Outputs:
            %   fieldName - Base field name string
            %   varIdx - Numerical index (double)

            tokens = regexp(str, '^([a-zA-Z_]+)(\d+)$', 'tokens');
            core.except.assert(~isempty(tokens), 'InvalidInput', ...
                'Input does not match the expected format.');
            fieldName = tokens{1}{1};
            varIdx = str2double(tokens{1}{2});
        end

        function displayName = formatDisplayName(obj, name)
            % FORMATDISPLAYNAME Format name for display in legends/titles.
            %
            %   displayName = formatDisplayName(obj, name) converts a
            %   dataset or variable name to a display-friendly format by
            %   replacing underscores with hyphens.
            %
            % Inputs:
            %   obj - The Strategy object
            %   name - Original name string
            %
            % Outputs:
            %   displayName - Formatted name for display

            if obj.REPLACE_UNDERSCORE_WITH_HYPHEN
                displayName = strrep(name, '_', '-');
            else
                displayName = name;
            end
        end

        function style = getDefaultLineStyle(obj, colorIdx, plotIdx)
            % GETDEFAULTLINESTYLE Get default line style for dataset.
            %
            %   style = getDefaultLineStyle(obj, dataType, offset)
            %   returns appropriate default styling based on data type and
            %   the current offset for line plotting type.
            %
            % Inputs:
            %   obj - The Strategy object
            %   colorIdx - Current color index
            %   plotIdx - Current plot index
            %
            % Outputs:
            %   style - Cell array of style arguments for plot functions

            nColors = length(obj.DEFAULT_COLOR_PALETTE);
            nLines = length(obj.DEFAULT_LINES);
            colorIdx = mod(colorIdx - 1, nColors) + 1;
            plotIdx = mod(plotIdx - 1, nLines) + 1;
            
            style = {
                'Color', obj.DEFAULT_COLOR_PALETTE{colorIdx}, ...
                'LineStyle', obj.DEFAULT_LINES{plotIdx}, ...
                'Marker', 'none', ...
                'LineWidth', obj.DEFAULT_LINE_WIDTH
            };

        end

        function style = getDefaultScatterStyle(obj, colorIdx, plotIdx)
            % GETDEFAULTSCATTERSTYLE Get default scatter style for dataset.
            %
            %   style = getDefaultScatterStyle(obj, colorIdx, plotIdx)
            %   returns appropriate default styling based on data type and
            %   the current offset for scatter plotting type.
            %
            % Inputs:
            %   obj - The Strategy object
            %   colorIdx - Current color index
            %   plotIdx - Current plot index
            %
            % Outputs:
            %   style - Cell array of style arguments for plot functions

            nColors = length(obj.DEFAULT_COLOR_PALETTE);
            nMarkers = length(obj.DEFAULT_MARKERS);
            colorIdx = mod(colorIdx - 1, nColors) + 1;
            plotIdx = mod(plotIdx - 1, nMarkers) + 1;
            
            style = {
                'Color', obj.DEFAULT_COLOR_PALETTE{colorIdx}, ...
                'LineStyle', 'none', ...
                'Marker', obj.DEFAULT_MARKERS{plotIdx}, ...
                'LineWidth', obj.DEFAULT_LINE_WIDTH
            };
        end

        function setColorLimits(obj, ax, u)
            % SETCOLORLIMITS Set color axis limits based on data values.
            %
            %   setColorLimits(obj, ax, u) configures the color axis limits
            %   for the specified axes based on the data range, using
            %   symmetric limits for bipolar data.
            %
            % Inputs:
            %   obj - The Strategy object
            %   ax - Axes handle to configure
            %   u - Data array for limit calculation

            uMin = min(u(:));
            uMax = max(u(:));
            absMax = max(abs(uMin), abs(uMax));

            if obj.USE_SYMMETRIC_COLOR_LIMITS && uMin < 0 && uMax > 0
                cLimits = [-absMax, absMax];
            else
                padding = obj.DEFAULT_COLOR_PADDING * (uMax - uMin);
                if padding == 0
                    padding = obj.DEFAULT_COLOR_PADDING * abs(uMin);
                    if padding == 0
                        padding = 0.1;
                    end
                end
                cLimits = [uMin - padding, uMax + padding];
            end

            caxis(ax, cLimits);
        end

        function syncColorLimits(obj, axes, subdata)
            % SYNCCOLORLIMITS Set unified color limits for multiple plots.
            %
            %   syncColorLimits(obj, axes, subdata) synchronizes color
            %   limits across multiple axes based on the combined data range.
            %
            % Inputs:
            %   obj - The Strategy object
            %   axes - Array of axes handles to synchronize
            %   subdata - Map of data names to structures containing data

            if isempty(fieldnames(subdata))
                return;
            end

            uMin = inf;
            uMax = -inf;

            dataNames = fieldnames(subdata);
            for i = 1:length(dataNames)
                dataName = dataNames{i};
                u = subdata.(dataName).u;
                uMin = min(uMin, min(u(:)));
                uMax = max(uMax, max(u(:)));
            end

            absMax = max(abs(uMin), abs(uMax));

            if obj.USE_SYMMETRIC_COLOR_LIMITS && uMin < 0 && uMax > 0
                cLimits = [-absMax, absMax];
            else
                padding = obj.DEFAULT_COLOR_PADDING * (uMax - uMin);
                if padding == 0
                    padding = obj.DEFAULT_COLOR_PADDING * abs(uMin);
                    if padding == 0
                        padding = 0.1;
                    end
                end
                cLimits = [uMin - padding, uMax + padding];
            end

            for i = 1:length(axes)
                caxis(axes(i), cLimits);
            end
        end

        function finalizeAxis(obj, ax, xmin, xmax, ymin, ymax, xLabel, yLabel, legends)
            % FINALIZEAXIS Apply common formatting to an axis.
            %
            %   finalizeAxis(obj, ax, xmin, xmax, ymin, ymax, xLabel, yLabel, legends)
            %   applies standard formatting including limits, labels, and legends.
            %
            % Inputs:
            %   obj - The Strategy object
            %   ax - Axes handle to format
            %   xmin, xmax - X-axis limits
            %   ymin, ymax - Y-axis limits
            %   xLabel, yLabel - Axis labels
            %   legends - Cell array of legend entries

            epsilon = 1e-8;
            xmin = min(-epsilon, xmin);
            xmax = max(epsilon, xmax);
            ymin = min(-epsilon, ymin);
            ymax = max(epsilon, ymax);

            %< Set axis limits with log scale support
            if obj.ENABLE_LOG_PLOT && obj.logScale.x
                set(ax, 'XScale', 'log');
                xlim(ax, [max(epsilon, xmin), xmax]);
            else
                xlim(ax, [xmin, xmax]);
            end

            if obj.ENABLE_LOG_PLOT && obj.logScale.y
                set(ax, 'YScale', 'log');
                ylim(ax, [max(epsilon, ymin), ymax]);
            else
                ylim(ax, [ymin, ymax]);
            end

            %< Set labels with formatting
            xlabel(ax, xLabel, 'FontSize', obj.DEFAULT_AXIS_LABEL_FONT_SIZE);
            ylabel(ax, yLabel, 'FontSize', obj.DEFAULT_AXIS_LABEL_FONT_SIZE);

            %< Configure legend
            if ~isempty(legends)
                legend(ax, legends, 'Location', obj.DEFAULT_LEGEND_POSITION, ...
                    'FontSize', obj.DEFAULT_LEGEND_FONT_SIZE);
            end

            hold(ax, 'off');
        end
    end

    methods (Access = private)
        function obj = renderSequential(obj, axes, dataset, varNames)
            % RENDERSEQUENTIAL Render variables in sequential order.

            for j = 1:length(varNames)
                varName = varNames{j};

                nAxes = length(fieldnames(dataset)) + 1;
                offset = 1 + (j - 1) * nAxes;

                if offset + nAxes - 1 > length(axes)
                    core.except.verify(false, 'NotEnoughAxes', ...
                        'Not enough axes for all variables.');
                    break;
                end

                obj.renderVariable(axes, dataset, varName, offset);
            end
        end

        function renderMapped(obj, axes, dataset, varNames)
            % RENDERMAPPED Render variables using explicit axis mapping.

            for j = 1:length(varNames)
                varName = varNames{j};
                offset = obj.varMap.(varName);

                if offset > length(axes)
                    core.except.verify(false, 'InvalidAxisIdx', ...
                        'Axis index %d exceeds available axes.', offset);
                    continue;
                end

                obj.renderVariable(axes, dataset, varName, offset);
            end
        end
    
        function renderVariable(obj, axes, dataset, varName, offset)
            % RENDERVARIABLE Render a single variable onto axes.

            switch obj.nDims
                case 1
                    obj.renderVariable1d(axes, dataset, varName, offset);
                case 2
                    obj.renderVariable2d(axes, dataset, varName, offset);
                case 3
                    obj.renderVariable3d(axes, dataset, varName, offset);
            end
        end

        function renderVariable1d(obj, axes, dataset, varName, offset)
            % RENDERVARIABLE1D Render 1D variable data.

            subdata = obj.prepareData(dataset, varName);

            if isempty(fieldnames(subdata))
                return;
            end

            obj.render1d(axes, subdata, varName, offset);
        end

        function renderVariable2d(obj, axes, dataset, varName, offset)
            % RENDERVARIABLE2D Render 2D variable data.

            subdata = obj.prepareData(dataset, varName);

            if isempty(fieldnames(subdata))
                return;
            end

            varAxes = obj.render2d(axes, subdata, varName, offset);
            obj.renderSlices(axes, subdata, varName, offset, 2);

            if length(varAxes) > 1
                obj.syncColorLimits(varAxes, subdata);
            end
        end

        function renderVariable3d(obj, axes, dataset, varName, offset)
            % RENDERVARIABLE3D Render 3D variable data.

            subdata = obj.prepareData(dataset, varName);

            if isempty(fieldnames(subdata))
                return;
            end

            varAxes = obj.render3d(axes, subdata, varName, offset);
            obj.renderSlices(axes, subdata, varName, offset, 3);

            if length(varAxes) > 1
                obj.syncColorLimits(varAxes, subdata);
            end
        end

        function ax = render1d(obj, axes, subdata, varName, offset)
            % RENDER1D Render 1D line plot with default styling.

            dataNames = fieldnames(subdata);
            n = length(dataNames);
            ax = axes(offset);
            cla(ax);
            hold(ax, 'on');

            [xmin, xmax, umin, umax] = deal(inf, -inf, inf, -inf);
            legends = cell(1, n);

            lineIdx = 1;
            scatterIdx = 1;
            for i = 1:n
                dataName = dataNames{i};
                g = subdata.(dataName);

                xmin = min(xmin, min(g.x1));
                xmax = max(xmax, max(g.x1));
                umin = min(umin, min(g.u));
                umax = max(umax, max(g.u));

                if g.type < 2
                    style = obj.getDefaultLineStyle(i, lineIdx);
                    lineIdx = lineIdx + 1;
                else
                    style = obj.getDefaultScatterStyle(i, scatterIdx);
                    scatterIdx = scatterIdx + 1;
                end
                
                plot(ax, g.x1(:).', g.u(:).', style{:});
                legends{i} = obj.formatDisplayName(dataName);
            end

            obj.finalizeAxis(ax, xmin, xmax, umin, umax, 'x', varName, legends);
        end

        function varAxes = render2d(obj, axes, subdata, varName, offset)
            % RENDER2D Render 2D surface plots with default styling.

            dataNames = fieldnames(subdata);
            varAxes = [];
            
            for i = 1:length(dataNames)
                dataName = dataNames{i};
                g = subdata.(dataName);
                ax = axes(offset + i - 1);
                varAxes = [varAxes, ax];

                cla(ax);
                hold(ax, 'on');

                imagesc(ax, g.x1, g.x2, g.u');
                axis(ax, 'xy', 'equal', 'tight');
                xlabel(ax, 'x', 'FontSize', obj.DEFAULT_AXIS_LABEL_FONT_SIZE);
                ylabel(ax, 'y', 'FontSize', obj.DEFAULT_AXIS_LABEL_FONT_SIZE);
                colormap(ax, obj.DEFAULT_COLORMAP_2D);
                colorbar(ax);

                obj.setColorLimits(ax, g.u);

                titleText = sprintf('%s: %s', obj.formatDisplayName(dataName), varName);
                title(ax, titleText);

                hold(ax, 'off');
            end
        end

        function varAxes = render3d(obj, axes, subdata, varName, offset)
            % RENDER3D Render 3D slice plots with default styling.

            varAxes = [];
            dataNames = fieldnames(subdata);
            
            for i = 1:length(dataNames)
                dataName = dataNames{i};
                g = subdata.(dataName);
                ax = axes(offset + i - 1);
                varAxes = [varAxes, ax];

                cla(ax);
                hold(ax, 'on');

                z1 = (min(g.x1) + max(g.x1)) / 2;
                z2 = (min(g.x2) + max(g.x2)) / 2;
                z3 = (min(g.x3) + max(g.x3)) / 2;
                u = permute(g.u, [2, 1, 3]);
                h = slice(ax, g.x1, g.x2, g.x3, u, z1, z2, z3);
                set(h, 'EdgeColor', 'none');

                xlabel(ax, 'x', 'FontSize', obj.DEFAULT_AXIS_LABEL_FONT_SIZE);
                ylabel(ax, 'y', 'FontSize', obj.DEFAULT_AXIS_LABEL_FONT_SIZE);
                zlabel(ax, 'z', 'FontSize', obj.DEFAULT_AXIS_LABEL_FONT_SIZE);
                colormap(ax, obj.DEFAULT_COLORMAP_3D);
                colorbar(ax);
                axis(ax, 'vis3d', 'equal', 'tight');
                obj.setColorLimits(ax, u);
                
                titleText = sprintf('%s: %s', obj.formatDisplayName(dataName), varName);
                title(ax, titleText);
                view(ax, 3);

                hold(ax, 'off');
            end
        end

        function renderSlices(obj, axes, subdata, varName, offset, nDims)
            % RENDERSLICES Render diagonal slice plots for multi-D data.

            dataNames = fieldnames(subdata);
            n = length(dataNames);
            if n < 1
                return;
            end

            ax = axes(offset + n);
            cla(ax);
            hold(ax, 'on');

            [xmin, xmax, ymin, ymax] = deal(inf, -inf, inf, -inf);
            legends = cell(1, n);

            lineIdx = 1;
            scatterIdx = 1;
            for i = 1:n
                dataName = dataNames{i};
                g = subdata.(dataName);

                r = linspace(min(g.x1), max(g.x1), 64);
                if nDims == 2
                    R = sqrt(2);
                    s = interp2(g.x1, g.x2, g.u.', r/R, r/R);
                elseif nDims == 3
                    R = sqrt(3);
                    u = permute(g.u, [2, 1, 3]);
                    s = interp3(g.x1, g.x2, g.x3, u, r/R, r/R, r/R);
                end

                %< Use default style based on data type
                if g.type < 2
                    style = obj.getDefaultLineStyle(i, lineIdx);
                else
                    style = obj.getDefaultScatterStyle(i, scatterIdx);
                end
                plot(ax, r, s, style{:}, 'DisplayName', dataName);

                xmin = min(xmin, min(r(:)));
                xmax = max(xmax, max(r(:)));
                if ~isempty(s)
                    ymin = min(ymin, min(s(:)));
                    ymax = max(ymax, max(s(:)));
                end

                legends{i} = obj.formatDisplayName(dataName);
            end

            if ymin == ymax
                ymin = ymin - 1e-8;
                ymax = ymax + 1e-8;
            end

            obj.finalizeAxis(ax, xmin, xmax, ymin, ymax, 'x', varName, legends);
            title(ax, 'diagonal slice');
        end
    end
end