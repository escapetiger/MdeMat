classdef StrategyTimeSeries1d < profilers.visual.Strategy
    % STRATEGYTIMESERIES1D Time series line plot rendering strategy.
    %
    %   StrategyTimeSeries1d extends the base Strategy class to provide
    %   specialized rendering for time series data using line plots.
    %   The strategy creates separate subplots for each score/property
    %   with consistent styling and automatic axis formatting.
    %
    %   Each score is rendered as a line plot with markers to show
    %   individual data points. The strategy supports multiple scores
    %   with automatic grid layout and synchronized time axes.
    %
    % Examples:
    %   % Create time series strategy and render data
    %   strategy = StrategyTimeSeries1d();
    %   strategy.render(axesArray, timeSeriesData, scoreNames, offset);
    %
    % See Also:
    %   profilers.visual.Strategy, profilers.visual.Strategy1d

    properties (Constant)
        %< Time series specific styling
        TIME_SERIES_LINE_WIDTH = 1.5; % Line width for time series plots
        TIME_SERIES_MARKER_SIZE = 4; % Marker size for data points
        TIME_SERIES_GRID_ALPHA = 0.3; % Grid transparency level
        TIME_SERIES_PADDING = 0.05; % Axis padding as fraction of range
    end

    methods
        function obj = render(obj, axes, timeSeries, scoreNames, offset)
            % RENDER Render time series plots for specified scores.
            %
            %   obj = render(obj, axes, timeSeries, scoreNames, offset)
            %   creates time series line plots for each score in the
            %   scoreNames list. Each score gets its own subplot with
            %   consistent styling and proper axis formatting.
            %
            % Inputs:
            %   obj - The StrategyTimeSeries1d object
            %   axes - Array of axes handles for rendering
            %   timeSeries - Structure containing 't' and score data
            %   scoreNames - Cell array of score names to plot
            %   offset - Starting axis offset for multi-panel layouts (scalar)
            %
            % Outputs:
            %   obj - The StrategyTimeSeries1d object

            core.except.assert(isstruct(timeSeries) && isfield(timeSeries, 't'), ...
                'InvalidInput', 'timeSeries must be a struct with ''t'' field');
            core.except.assert(iscell(scoreNames), ...
                'InvalidInput', 'scoreNames must be a cell array');

            if isempty(timeSeries.t)
                return; %< No data to plot
            end

            timeVector = timeSeries.t;
            nScores = length(scoreNames);

            %< Render each score on its designated axis
            for iScore = 1:nScores
                scoreName = scoreNames{iScore};

                %< Skip if score data not available
                if ~isfield(timeSeries, scoreName)
                    continue;
                end

                scoreData = timeSeries.(scoreName);

                %< Skip if score data is empty or wrong size
                if isempty(scoreData) || length(scoreData) ~= length(timeVector)
                    continue;
                end

                ax = axes(offset+iScore-1);

                %< Clear and prepare axes
                cla(ax);
                hold(ax, 'on');

                %< Get consistent styling from base strategy
                style = obj.getDefaultLineStyle(colorIndex, 1);

                %< Create line plot with markers for data points
                plot(ax, timeVector, scoreData, style{:}, ...
                    'LineWidth', obj.TIME_SERIES_LINE_WIDTH, ...
                    'Marker', 'o', ...
                    'MarkerSize', obj.TIME_SERIES_MARKER_SIZE, ...
                    'MarkerFaceColor', style{2}); %< Use line color for marker fill

                %< Format axes with consistent styling
                obj.formatScoreAxis(ax, scoreName, timeVector, scoreData);

                hold(ax, 'off');
            end
        end
    end

    methods (Access = protected)
        function obj = formatScoreAxis(obj, ax, scoreName, timeVector, scoreData)
            % FORMATSCOREAXIS Apply consistent formatting to score axis.
            %
            %   formatScoreAxis(obj, ax, scoreName, timeVector, scoreData)
            %   applies standard formatting including labels, grid, and
            %   limits to the specified score axis.
            %
            % Inputs:
            %   obj - The StrategyTimeSeries1d object
            %   ax - Axes handle to format
            %   scoreName - Name of the score for labeling
            %   timeVector - Time values for axis limits
            %   scoreData - Score values for axis limits
            %
            % Outputs:
            %   obj - The StrategyTimeSeries1d object

            %< Set axis labels with consistent font sizing
            xlabel(ax, 'Time', 'FontSize', obj.DEFAULT_AXIS_LABEL_FONT_SIZE);
            ylabel(ax, strrep(scoreName, '_', '-'), ...
                'FontSize', obj.DEFAULT_AXIS_LABEL_FONT_SIZE);

            %< Apply grid with transparency
            grid(ax, 'on');
            set(ax, 'GridAlpha', obj.TIME_SERIES_GRID_ALPHA);

            %< Set reasonable axis limits with padding
            if length(timeVector) > 1
                %< Time axis limits
                timeRange = [min(timeVector), max(timeVector)];
                timePadding = (timeRange(2) - timeRange(1)) * obj.TIME_SERIES_PADDING;
                if timePadding == 0
                    timePadding = 0.1; %< Minimum padding for constant time
                end
                xlim(ax, [timeRange(1) - timePadding, timeRange(2) + timePadding]);

                %< Score axis limits
                scoreRange = [min(scoreData), max(scoreData)];
                scorePadding = (scoreRange(2) - scoreRange(1)) * obj.TIME_SERIES_PADDING;
                if scorePadding == 0
                    scorePadding = 0.1 * abs(scoreRange(1)); %< Handle constant data
                    if scorePadding == 0
                        scorePadding = 0.1; %< Minimum padding for zero data
                    end
                end
                ylim(ax, [scoreRange(1) - scorePadding, scoreRange(2) + scorePadding]);
            else
                %< Single data point case
                xlim(ax, [timeVector(1) - 0.1, timeVector(1) + 0.1]);
                ylim(ax, [scoreData(1) - 0.1, scoreData(1) + 0.1]);
            end

            %< Ensure proper aspect ratio for time series
            axis(ax, 'normal');
        end
    end
end