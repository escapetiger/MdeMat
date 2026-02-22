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
    %   physics.visual.Strategy1d, physics.visual.Strategy2d,
    %   physics.visual.Strategy3d, physics.visual.StrategySlice1d

    properties (Constant)
        %< Paul Tol discrete color schemes (embedded for offline use)
        TolBrightColors = {'#000000', '#EE6677', '#228833', '#4477AA', '#CCBB44', ...
            '#66CCEE', '#AA3377', '#BBBBBB'};
        TolMutedColors = {'#000000', '#CC6677', '#332288', '#DDCC77', '#117733', ...
            '#88CCEE', '#882255', '#44AA99', '#999933', ...
            '#AA4499', '#DDDDDD'};
        TolVibrantColors = {'#000000', '#EE7733', '#0077BB', '#33BBEE', '#EE3377', ...
            '#CC3311', '#009988', '#BBBBBB'};
        TolLightColors = {'#000000', '#EEDD88', '#99DDFF', '#44BB99', '#FFAABB', ...
            '#EECC66', '#77AADD', '#EE8866', '#DDDDDD'};
        TolHighContrastColors = {'#000000', '#004488', '#BB5566', '#DDAA33', '#FFFFFF'}

        %< Francis Filbet's discrete color scheme
        FilbetColors = {
            '#000000', ... % black
            '#CBCBFF', ... % light blue
            '#D04231', ... % red
            '#EF9F93', ... % rose
            '#AE7AC5', ... % lavander
            '#5F94EE', ... % blue
            '#9C640D', ... % brown
            '#222F3D', ... % blue-gray
            '#8A8A8A', ... % gray
        };

        %< Paul Tol colormap data (subset for key scientific visualization maps)
        %< Sunset colormap: dark blues (negative) to oranges/reds (positive)
        TolSunsetMap = [; ...
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
        TolBurdMap = [; ...
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
        TolIridescentMap = [; ...
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
        DiscreteColorScheme = 'bright'; % Current discrete color scheme
        ContinuousColorScheme = 'jet'; % Current continuous color scheme

        %< Marker and line style patterns
        ScatterScheme = {'o', 's', '^', 'v', 'd', 'x', '+'}; % Marker styles for scatter plots
        % LineScheme = {'-', '--', ':', '-.', '-', '--', ':'}; % Line styles for line plots
        LineScheme = {'-'}

        %< Default styling parameters
        DefaultLineWidth = 2; % Standard line width for plots

        %< Color mapping configuration
        DefaultColorPadding = 0.05; % 5% padding for color axis limits
        UseSymmetricColorLimits = false; % Use symmetric limits for bipolar data

        %< Typography and layout settings
        DefaultFontSize = 18; % Font size for any elements
        DefaultLegendPosition = 'best'; % Automatic legend positioning
        DefaultAxisLabelFontSize = 28; % Font size for axis labels
        DefaultLegendFontSize = 18; % Font size for legend text
        DefaultTitleFontSize = 28; % Font size for title
    end

    properties (Dependent)
        ColorPalette % Discrete color palette based on current scheme
        ColorMap % Continuous color map based on current scheme
    end

    methods
        function cmap = get.ColorMap(obj)
            % GET.CMAP Get continuous colormap based on current scheme.

            switch obj.ContinuousColorScheme
                case 'sunset'
                    cmap = obj.TolSunsetMap;
                case 'burd'
                    cmap = obj.TolBurdMap;
                case 'iridescent'
                    cmap = obj.TolIridescentMap;
                case 'turbo'
                    cmap = 'turbo';
                case 'jet'
                    cmap = 'jet';
                otherwise
                    cmap = 'viridis'; %< Fallback to MATLAB default
            end
        end

        function cpal = get.ColorPalette(obj)
            % GET.CPAL Get discrete color palette based on current scheme.
            switch obj.DiscreteColorScheme
                case 'bright'
                    cpal = obj.TolBrightColors;
                case 'muted'
                    cpal = obj.TolMutedColors;
                case 'vibrant'
                    cpal = obj.TolVibrantColors;
                case 'light'
                    cpal = obj.TolLightColors;
                case 'high_contrast'
                    cpal = obj.TolHighContrastColors;
                case 'filbet'
                    cpal = obj.FilbetColors;
                otherwise
                    cpal = obj.TolBrightColors; %< Default fallback
            end
        end

        function style = getDefaultLineStyle(obj, colorIdx, plotIdx)
            % GETDEFAULTLINESTYLE Get default line style for dataset.
            %
            %   style = getDefaultLineStyle(obj, colorIdx, plotIdx)
            %   returns appropriate default styling parameters for line
            %   plots, combining color and line style selections.

            arguments
                obj physics.visual.Strategy
                colorIdx {mustBePositive, mustBeInteger}
                plotIdx {mustBePositive, mustBeInteger}
            end

            colorIdx = mod(colorIdx-1, length(obj.ColorPalette))+1;
            plotIdx = mod(plotIdx-1, length(obj.LineScheme))+1;
            if iscell(obj.ColorPalette)
                color = obj.ColorPalette{colorIdx};
            else
                color = obj.ColorPalette(colorIdx);
            end
            if iscell(obj.ColorPalette)
                lineStyle = obj.LineScheme{plotIdx};
            else
                lineStyle = obj.LineScheme(plotIdx);
            end

            style = {; ...
                'Color', color, ...
                'LineStyle', lineStyle, ...
                'Marker', 'none', ...
                'LineWidth', obj.DefaultLineWidth; ...
                };
        end

        function style = getDefaultScatterStyle(obj, colorIdx, plotIdx)
            % GETDEFAULTSCATTERSTYLE Get default scatter style for dataset.
            %
            %   style = getDefaultScatterStyle(obj, colorIdx, plotIdx)
            %   returns appropriate default styling parameters for scatter
            %   plots, combining color and marker selections.

            arguments
                obj physics.visual.Strategy
                colorIdx {mustBePositive, mustBeInteger}
                plotIdx {mustBePositive, mustBeInteger}
            end

            colorIdx = mod(colorIdx-1, length(obj.ColorPalette))+1;
            plotIdx = mod(plotIdx-1, length(obj.ScatterScheme))+1;
            if iscell(obj.ColorPalette)
                color = obj.ColorPalette{colorIdx};
            else
                color = obj.ColorPalette(colorIdx);
            end
            if iscell(obj.ColorPalette)
                marker = obj.ScatterScheme{plotIdx};
            else
                marker = obj.ScatterScheme(plotIdx);
            end

            style = {; ...
                'Color', color, ...
                'LineStyle', 'none', ...
                'Marker', marker, ...
                'LineWidth', obj.DefaultLineWidth; ...
                };
        end

        function obj = setColorLimits(obj, ax, u)
            % SETCOLORLIMITS Set color axis limits based on data values.
            %
            %   setColorLimits(obj, ax, u) configures the color axis limits
            %   for the specified axes based on the data range. Uses
            %   symmetric limits for bipolar data when enabled and applies
            %   padding to prevent visual clipping.

            arguments
                obj physics.visual.Strategy
                ax matlab.graphics.axis.Axes
                u double
            end

            uMin = min(u(:));
            uMax = max(u(:));
            absMax = max(abs(uMin), abs(uMax));

            %< Use symmetric limits for bipolar data
            if obj.UseSymmetricColorLimits && uMin < 0 && uMax > 0
                cLimits = [-absMax, absMax];
            else
                %< Apply padding to data range
                padding = obj.DefaultColorPadding * (uMax - uMin);
                if padding == 0
                    %< Handle constant data case
                    padding = obj.DefaultColorPadding * abs(uMin);
                    if padding == 0
                        padding = 0.1; %< Minimum padding for zero data
                    end
                end
                cLimits = [uMin - padding, uMax + padding];
            end

            caxis(ax, cLimits);
        end

        function obj = syncColorLimits(obj, axes, database)
            % SYNCCOLORLIMITS Set unified color limits for multiple plots.
            %
            %   obj = syncColorLimits(obj, axes, database) synchronizes
            %   color limits across multiple @a axes based on the combined
            %   data range from @a database. This ensures consistent color
            %   mapping across related visualizations.

            arguments
                obj physics.visual.Strategy
                axes matlab.graphics.axis.Axes
                database physics.visual.Database
            end

            if database.IsEmpty || isempty(axes)
                return;
            end

            uMin = inf;
            uMax = -inf;
            datasetNames = database.DatasetNames;
            for i = 1:length(datasetNames)
                dataset = database.getDataset(datasetNames{i});
                u = dataset.getData('u');
                uMin = min(uMin, min(u(:)));
                uMax = max(uMax, max(u(:)));
            end

            absMax = max(abs(uMin), abs(uMax));

            %< Apply consistent limit calculation logic
            if obj.UseSymmetricColorLimits && uMin < 0 && uMax > 0
                cLimits = [-absMax, absMax];
            else
                padding = obj.DefaultColorPadding * (uMax - uMin);
                if padding == 0
                    padding = obj.DefaultColorPadding * abs(uMin);
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
    end
end