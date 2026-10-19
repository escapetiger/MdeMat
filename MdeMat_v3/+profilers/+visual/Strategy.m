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
    % Examples:
    %   % Create strategy and configure log scaling
    %   strategy = Strategy();
    %   strategy.setLogScale('y', true);
    %   
    %   % Get styling for different plot types
    %   lineStyle = strategy.getDefaultLineStyle(1, 1);
    %   scatterStyle = strategy.getDefaultScatterStyle(2, 1);
    %   
    %   % Apply color limits and formatting
    %   strategy.setColorLimits(ax, dataArray);
    %   strategy.finalizeAxis(ax, 0, 1, 0, 1, 'X', 'Y', {'Series 1'});
    %
    % See Also:
    %   profilers.visual.Strategy1d, profilers.visual.Strategy2d,
    %   profilers.visual.Strategy3d, profilers.visual.StrategySlice1d

    properties (Constant)
        %< Paul Tol discrete color schemes (embedded for offline use)
        TOL_BRIGHT_COLORS = {'#000000', '#EE6677', '#228833', '#4477AA', '#CCBB44', ...
            '#66CCEE', '#AA3377', '#BBBBBB'};
        TOL_MUTED_COLORS = {'#000000', '#CC6677', '#332288', '#DDCC77', '#117733', ...
            '#88CCEE', '#882255', '#44AA99', '#999933', ...
            '#AA4499', '#DDDDDD'};
        TOL_VIBRANT_COLORS = {'#000000', '#EE7733', '#0077BB', '#33BBEE', '#EE3377', ...
            '#CC3311', '#009988', '#BBBBBB'};
        TOL_LIGHT_COLORS = {'#000000', '#EEDD88', '#99DDFF', '#44BB99', '#FFAABB', ...
            '#EECC66', '#77AADD', '#EE8866', '#DDDDDD'};
        TOL_HIGH_CONTRAST_COLORS = {'#000000', '#004488', '#BB5566', '#DDAA33', '#FFFFFF'}

        %< Paul Tol colormap data (subset for key scientific visualization maps)
        %< Sunset colormap: dark blues (negative) to oranges/reds (positive)
        TOL_SUNSET_MAP = [; ...
            54, 75, 154; ... % dark blue (strong negative)
            64, 90, 168; ... % slightly lighter blue
            74, 105, 183; ... % medium-dark blue
            84, 120, 198; ... % medium blue
            94, 135, 205; ... % light blue
            110, 150, 215; ... % lighter blue
            126, 166, 220; ... % very light blue
            142, 182, 225; ... % pale blue
            158, 198, 230; ... % very pale blue
            174, 214, 235; ... % blue-white
            190, 225, 238; ... % very pale blue
            206, 232, 240; ... % almost white blue
            222, 238, 242; ... % white-blue
            238, 244, 244; ... % near white
            245, 240, 220; ... % warm white (transition)
            252, 236, 196; ... % pale warm
            254, 228, 172; ... % very light orange
            254, 218, 148; ... % light orange
            254, 208, 124; ... % medium-light orange
            254, 198, 100; ... % medium orange
            253, 188, 88; ... % orange
            252, 178, 76; ... % orange-red
            251, 168, 64; ... % red-orange
            248, 150, 52; ... % red
            245, 132, 40; ... % darker red
            240, 114, 28; ... % dark red
            235, 96, 16; ... % darker red
            225, 78, 8; ... % very dark red
            215, 60, 0; ... % darkest red
            180, 40, 0; ... % maroon
            145, 20, 0; ... % dark maroon (strong positive)
            ] / 255;

        %< Blue-red diverging colormap for bipolar data
        TOL_BURD_MAP = [; ...
            33, 102, 172; ... % dark blue (negative extreme)
            67, 147, 195; ... % blue
            146, 197, 222; ... % light blue
            209, 229, 240; ... % very light blue
            247, 247, 247; ... % white (neutral)
            253, 219, 199; ... % very light red
            244, 165, 130; ... % light red
            214, 96, 77; ... % red
            178, 24, 43; ... % dark red (positive extreme)
            ] / 255;

        %< Iridescent colormap for smooth gradations
        TOL_IRIDESCENT_MAP = [; ...
            254, 251, 233; ... % very light yellow (low values)
            252, 247, 213; ... % light yellow
            245, 243, 193; ... % yellow
            234, 240, 181; ... % yellow-green
            221, 236, 191; ... % light green
            208, 231, 202; ... % green
            194, 227, 210; ... % green-cyan
            181, 221, 216; ... % cyan
            168, 216, 220; ... % light cyan
            155, 210, 225; ... % cyan-blue
            141, 203, 228; ... % light blue
            129, 196, 231; ... % blue
            123, 188, 231; ... % medium blue
            126, 178, 228; ... % blue-purple
            136, 165, 221; ... % purple-blue
            147, 152, 210; ... % purple
            155, 138, 196; ... % light purple
            157, 125, 178; ... % purple-pink
            154, 112, 158; ... % pink-purple
            144, 99, 136; ... % dark purple
            128, 87, 112; ... % very dark purple
            104, 73, 87; ... % dark grey
            70, 53, 58; ... % very dark (high values)
            ] / 255;

        %< Color scheme configuration
        DISCRETE_COLOR_SCHEME = 'bright'; % Current discrete color scheme
        CONTINUOUS_COLOR_SCHEME = 'iridescent'; % Current continuous colormap

        %< Marker and line style patterns
        SCATTER_SCHEME = {'o', 's', '^', 'v', 'd', 'x', '+'}; % Marker styles for scatter plots
        LINE_SCHEME = {'-', '--', ':', '-.', '-', '--', ':'}; % Line styles for line plots

        %< Default styling parameters
        DEFAULT_LINE_WIDTH = 1; % Standard line width for plots

        %< Color mapping configuration
        DEFAULT_COLOR_PADDING = 0.05; % 5% padding for color axis limits
        USE_SYMMETRIC_COLOR_LIMITS = true; % Use symmetric limits for bipolar data

        %< Typography and layout settings
        DEFAULT_LEGEND_POSITION = 'best'; % Automatic legend positioning
        DEFAULT_AXIS_LABEL_FONT_SIZE = 16; % Font size for axis labels
        DEFAULT_LEGEND_FONT_SIZE = 14; % Font size for legend text

        %< 1D plot specific options
        ENABLE_LOG_PLOT = true; % Enable logarithmic scaling support
    end

    properties
        logScale % Structure controlling log scale for 1D plots (x, y)
    end

    properties (Dependent)
        cpal % Discrete color palette based on current scheme
        cmap % Continuous color map based on current scheme
    end

    methods
        function obj = Strategy()
            % STRATEGY Constructor for Strategy base class.
            %
            %   obj = Strategy() creates a strategy object with default
            %   configuration for log scaling and color schemes.
            %
            % Outputs:
            %   obj - Constructed Strategy object

            obj.logScale = struct('x', false, 'y', false);
        end

        function cmap = get.cmap(obj)
            % GET.CMAP Get continuous colormap based on current scheme.
            %
            %   cmap = get.cmap(obj) returns the colormap array or name
            %   corresponding to the current continuous color scheme.
            %
            % Outputs:
            %   cmap - Colormap array (Nx3) or colormap name string

            switch obj.CONTINUOUS_COLOR_SCHEME
                case 'sunset'
                    cmap = obj.TOL_SUNSET_MAP;
                case 'burd'
                    cmap = obj.TOL_BURD_MAP;
                case 'iridescent'
                    cmap = obj.TOL_IRIDESCENT_MAP;
                case 'turbo'
                    cmap = 'turbo';
                otherwise
                    cmap = 'viridis'; %< Fallback to MATLAB default
            end
        end

        function cpal = get.cpal(obj)
            switch obj.DISCRETE_COLOR_SCHEME
                case 'bright'
                    cpal = obj.TOL_BRIGHT_COLORS;
                case 'muted'
                    cpal = obj.TOL_MUTED_COLORS;
                case 'vibrant'
                    cpal = obj.TOL_VIBRANT_COLORS;
                case 'light'
                    cpal = obj.TOL_LIGHT_COLORS;
                case 'high_contrast'
                    cpal = obj.TOL_HIGH_CONTRAST_COLORS;
                otherwise
                    cpal = obj.TOL_BRIGHT_COLORS; %< Default fallback
            end
        end

        function obj = setLogScale(obj, axis, isEnabled)
            % SETLOGSCALE Configure logarithmic scaling for 1D plots.
            %
            %   obj = setLogScale(obj, axis, isEnabled) enables or disables
            %   logarithmic scaling for the specified axis in 1D plots.
            %   This setting is used by the finalizeAxis method when
            %   ENABLE_LOG_PLOT is true.
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

        function style = getDefaultLineStyle(obj, colorIdx, plotIdx)
            % GETDEFAULTLINESTYLE Get default line style for dataset.
            %
            %   style = getDefaultLineStyle(obj, colorIdx, plotIdx)
            %   returns appropriate default styling parameters for line
            %   plots, combining color and line style selections.
            %
            % Inputs:
            %   obj - The Strategy object
            %   colorIdx - Color index for discrete color selection
            %   plotIdx - Plot index for line style selection
            %
            % Outputs:
            %   style - Cell array of property-value pairs for plot functions

            color = obj.cpal{mod(colorIdx-1, length(obj.cpal))+1};
            lineStyle = obj.LINE_SCHEME{mod(plotIdx-1, length(obj.LINE_SCHEME))+1};

            style = {; ...
                'Color', color, ...
                'LineStyle', lineStyle, ...
                'Marker', 'none', ...
                'LineWidth', obj.DEFAULT_LINE_WIDTH; ...
                };
        end

        function style = getDefaultScatterStyle(obj, colorIdx, plotIdx)
            % GETDEFAULTSCATTERSTYLE Get default scatter style for dataset.
            %
            %   style = getDefaultScatterStyle(obj, colorIdx, plotIdx)
            %   returns appropriate default styling parameters for scatter
            %   plots, combining color and marker selections.
            %
            % Inputs:
            %   obj - The Strategy object
            %   colorIdx - Color index for discrete color selection
            %   plotIdx - Plot index for marker style selection
            %
            % Outputs:
            %   style - Cell array of property-value pairs for plot functions

            color = obj.cpal{mod(colorIdx-1, length(obj.cpal))+1};
            marker = obj.SCATTER_SCHEME{mod(plotIdx-1, length(obj.SCATTER_SCHEME))+1};

            style = {; ...
                'Color', color, ...
                'LineStyle', 'none', ...
                'Marker', marker, ...
                'LineWidth', obj.DEFAULT_LINE_WIDTH; ...
                };
        end

        function obj = setColorLimits(obj, ax, u)
            % SETCOLORLIMITS Set color axis limits based on data values.
            %
            %   setColorLimits(obj, ax, u) configures the color axis limits
            %   for the specified axes based on the data range. Uses
            %   symmetric limits for bipolar data when enabled and applies
            %   padding to prevent visual clipping.
            %
            % Inputs:
            %   obj - The Strategy object
            %   ax - Axes handle to configure
            %   u - Data array for limit calculation
            %
            % Outputs:
            %   obj - The Strategy object

            uMin = min(u(:));
            uMax = max(u(:));
            absMax = max(abs(uMin), abs(uMax));

            %< Use symmetric limits for bipolar data
            if obj.USE_SYMMETRIC_COLOR_LIMITS && uMin < 0 && uMax > 0
                cLimits = [-absMax, absMax];
            else
                %< Apply padding to data range
                padding = obj.DEFAULT_COLOR_PADDING * (uMax - uMin);
                if padding == 0
                    %< Handle constant data case
                    padding = obj.DEFAULT_COLOR_PADDING * abs(uMin);
                    if padding == 0
                        padding = 0.1; %< Minimum padding for zero data
                    end
                end
                cLimits = [uMin - padding, uMax + padding];
            end

            caxis(ax, cLimits);
        end

        function obj = syncColorLimits(obj, axes, dataset)
            % SYNCCOLORLIMITS Set unified color limits for multiple plots.
            %
            %   obj = syncColorLimits(obj, axes, dataset) synchronizes
            %   color limits across multiple axes based on the combined
            %   data range from all datasets. This ensures consistent color
            %   mapping across related visualizations.
            %
            % Inputs:
            %   obj - The Strategy object
            %   axes - Array of axes handles to synchronize
            %   dataset - Map of data names to structures
            %
            % Outputs:
            %   obj - The Strategy object

            if isempty(fieldnames(dataset))
                return;
            end

            %< Find global data range across all datasets
            uMin = inf;
            uMax = -inf;

            dataNames = fieldnames(dataset);
            for i = 1:length(dataNames)
                dataName = dataNames{i};
                u = dataset.(dataName).u;
                uMin = min(uMin, min(u(:)));
                uMax = max(uMax, max(u(:)));
            end

            absMax = max(abs(uMin), abs(uMax));

            %< Apply consistent limit calculation logic
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

            %< Apply limits to all axes
            for i = 1:length(axes)
                caxis(axes(i), cLimits);
            end
        end

        function obj = finalizeAxis(obj, ax, xmin, xmax, ymin, ymax, xLabel, yLabel, legends)
            % FINALIZEAXIS Apply common formatting to an axis.
            %
            %   obj = finalizeAxis(obj, ax, xmin, xmax, ymin, ymax, xLabel,
            %   yLabel, legends) applies standard formatting including
            %   limits, labels, legends, and log scaling support. This
            %   method provides consistent appearance across all
            %   visualization strategies.
            %
            % Inputs:
            %   obj - The Strategy object
            %   ax - Axes handle to format
            %   xmin, xmax - X-axis limits (scalars)
            %   ymin, ymax - Y-axis limits (scalars)
            %   xLabel, yLabel - Axis label strings
            %   legends - Cell array of legend entries (may be empty)
            %
            % Outputs:
            %   obj - The Strategy object

            epsilon = 1e-8;
            %< Ensure non-zero ranges for proper display
            xmin = min(-epsilon, xmin);
            xmax = max(epsilon, xmax);
            ymin = min(-epsilon, ymin);
            ymax = max(epsilon, ymax);

            %< Configure X-axis with optional log scaling
            if obj.ENABLE_LOG_PLOT && obj.logScale.x
                set(ax, 'XScale', 'log');
                xlim(ax, [max(epsilon, xmin), xmax]); %< Ensure positive limits for log scale
            else
                xlim(ax, [xmin, xmax]);
            end

            %< Configure Y-axis with optional log scaling
            if obj.ENABLE_LOG_PLOT && obj.logScale.y
                set(ax, 'YScale', 'log');
                ylim(ax, [max(epsilon, ymin), ymax]); %< Ensure positive limits for log scale
            else
                ylim(ax, [ymin, ymax]);
            end

            %< Apply consistent label formatting
            xlabel(ax, xLabel, 'FontSize', obj.DEFAULT_AXIS_LABEL_FONT_SIZE);
            ylabel(ax, yLabel, 'FontSize', obj.DEFAULT_AXIS_LABEL_FONT_SIZE);

            %< Add legend if entries are provided
            if ~isempty(legends)
                legend(ax, legends, 'Location', obj.DEFAULT_LEGEND_POSITION, ...
                    'FontSize', obj.DEFAULT_LEGEND_FONT_SIZE);
            end

            hold(ax, 'off');
        end
    end
end