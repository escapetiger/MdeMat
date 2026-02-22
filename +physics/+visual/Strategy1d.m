classdef Strategy1d < physics.visual.Strategy
    % STRATEGY1D One-dimensional line plot rendering strategy.
    %
    %   Strategy1d extends the base Strategy class to provide specialized
    %   rendering for 1D data using line and scatter plots. The strategy
    %   automatically selects appropriate styling based on dataset types,
    %   using lines for analytical/exact data and markers for numerical
    %   data points.
    %
    %   The strategy supports multiple datasets on the same axes with
    %   automatic legend generation and consistent styling. Log scaling
    %   can be configured for both x and y axes independently.
    %
    % See also:
    %   physics.visual.Strategy, physics.visual.Strategy2d,
    %   physics.visual.StrategySlice1d
    
    properties (Constant)
        NDims = 1 % Number of spatial dimensions
        EnableLogPlot = true % Enable logarithmic scaling support
    end
    
    properties
        LogScale % Structure controlling log scale for 1D plots (x, y)
    end
    
    methods
        function obj = Strategy1d()
            % STRATEGY1D Construct an instance of Strategy1d.
            %
            %   obj = Strategy1d() creates a Strategy1d object with default
            %   configuration for log scaling disabled on both axes.
            
            obj@physics.visual.Strategy();
            obj.LogScale = struct('x', false, 'y', false);
        end
        
        function obj = setLogScale(obj, axis, isEnabled)
            % SETLOGSCALE Configure logarithmic scaling for 1D plots.
            %
            %   obj = setLogScale(obj, axis, isEnabled) enables or disables
            %   logarithmic scaling for the specified axis in 1D plots.
            %   This setting is used by the finalizeAxis method.
            
            arguments
                obj physics.visual.Strategy1d
                axis {mustBeMember(axis, {'x', 'y'})}
                isEnabled logical {mustBeScalarOrEmpty}
            end
            
            obj.LogScale.(axis) = isEnabled;
        end
        
        function obj = render(obj, ax, database, variable)
            % RENDER Render 1D line plot with automatic dataset styling.
            %
            %   obj = render(obj, ax, database, variable) creates a 1D line
            %   plot on the specified @a ax, automatically selecting line
            %   styles for analytical data and scatter styles for numerical
            %   data. Multiple datasets are overlaid with distinct styling
            %   and legend entries.
            
            arguments
                obj physics.visual.Strategy1d
                ax matlab.graphics.axis.Axes
                database physics.visual.Database
                variable {mustBeTextScalar}
            end
            
            dataNames = fieldnames(database.Datasets);
            
            cla(ax);
            hold(ax, 'on');
            
            options = struct();
            options.xmin = inf;
            options.xmax = -inf;
            options.ymin = inf;
            options.ymax = -inf;
            options.xLabel = 'x';
            options.yLabel = variable;
            options.legends = cell(1, length(dataNames));

            lineIdx = 1;
            scatterIdx = 1;
            for i = 1:length(dataNames)
                dataName = dataNames{i};
                dataset = database.getDataset(dataName);
                g = dataset.Data;
                
                options.xmin = min(options.xmin, min(g.x{1}));
                options.xmax = max(options.xmax, max(g.x{1}));
                options.ymin = min(options.ymin, min(g.u));
                options.ymax = max(options.ymax, max(g.u));
                
                if strcmp(dataset.Type, 'fixed') || strcmp(dataset.Type, 'exact')
                    style = obj.getDefaultLineStyle(i, lineIdx);
                    lineIdx = lineIdx + 1;
                else
                    style = obj.getDefaultScatterStyle(i, scatterIdx);
                    scatterIdx = scatterIdx + 1;
                end
                
                plot(ax, g.x{1}(:).', g.u(:).', style{:});
                options.legends{i} = strrep(dataName, '_', '-');
            end
            
            args = namedargs2cell(options);
            obj.formatAxis(ax, args{:});
        end
        
        function obj = formatAxis(obj, ax, options)
            % FORMATAXIS Apply common formatting to 1D axis with log
            % support.
            %
            %   obj = formatAxis(obj, axes, options) applies standard
            %   formatting including limits, labels, legends, and log
            %   scaling support for 1D plots. This is a specialization of
            %   the base finalize method for Strategy1d.
            
            arguments
                obj physics.visual.Strategy1d
                ax matlab.graphics.axis.Axes
                options.xmin double
                options.xmax double
                options.ymin double
                options.ymax double
                options.xLabel {mustBeTextScalar}
                options.yLabel {mustBeTextScalar}
                options.legends {mustBeA(options.legends, 'cell')} = {}
            end
            
            %< Extract parameters from options structure
            xmin = options.xmin;
            xmax = options.xmax;
            ymin = options.ymin;
            ymax = options.ymax;
            xLabel = options.xLabel;
            yLabel = options.yLabel;
            legends = options.legends;
            
            %< Ensure non-zero ranges for proper display
            tol = 1e-8;
            xmin = min(-tol, xmin);
            xmax = max(tol, xmax);
            ymin = min(-tol, ymin);
            ymax = max(tol, ymax);
            
            %< Configure X-axis with optional log scaling
            if obj.EnableLogPlot && obj.LogScale.x
                set(ax, 'XScale', 'log');
                xlim(ax, [max(tol, xmin), xmax]);
            else
                xlim(ax, [xmin, xmax]);
            end
            
            %< Configure Y-axis with optional log scaling
            if obj.EnableLogPlot && obj.LogScale.y
                set(ax, 'YScale', 'log');
                ylim(ax, [max(tol, ymin), ymax]);
            else
                ylim(ax, [ymin, ymax]);
            end
            
            %< Apply consistent label formatting
            xl = xlabel(ax, xLabel);
            yl = ylabel(ax, yLabel);
            set(xl, 'FontSize', obj.DefaultAxisLabelFontSize);
            set(xl, 'Interpreter', 'tex');
            set(yl, 'FontSize', obj.DefaultAxisLabelFontSize);
            set(yl, 'Interpreter', 'tex');
            
            %< Add legend if entries are provided
            if ~isempty(legends)
                leg = legend(ax, legends{:});
                set(leg, 'Interpreter', 'tex');
                set(leg, 'Location', obj.DefaultLegendPosition);
                set(leg, 'FontSize', obj.DefaultLegendFontSize);
            end
            
            hold(ax, 'off');
        end
    end
end