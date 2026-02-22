classdef Strategy2d < physics.visual.Strategy
    % STRATEGY2D Two-dimensional surface plot rendering strategy.
    %
    %   Strategy2d extends the base Strategy class to provide specialized
    %   rendering for 2D data using image-based surface plots. The strategy
    %   creates separate subplots for each dataset with synchronized color
    %   limits and consistent formatting.
    %
    %   Each dataset is rendered as a colored image with proper axis
    %   scaling, colorbar, and title. When multiple datasets are present,
    %   color limits are synchronized across all plots to enable meaningful
    %   comparison.
    %
    % See also:
    %   physics.visual.Strategy, physics.visual.Strategy1d,
    %   physics.visual.Strategy3d

    properties (Constant)
        NDims = 2 % Number of spatial dimensions
    end

    methods
        function obj = render(obj, axes, database, variable)
            % RENDER Render 2D surface plots with synchronized styling.
            %
            %   obj = render(obj, axes, database, variable) creates 2D
            %   surface visualizations for each dataset using imagesc with
            %   proper coordinate mapping and color synchronization. Each
            %   dataset receives its own subplot with consistent formatting
            %   and colorbar display.

            arguments
                obj physics.visual.Strategy2d
                axes matlab.graphics.axis.Axes
                database physics.visual.Database
                variable{mustBeTextScalar}
            end

            dataNames = fieldnames(database.Datasets);

            state = struct();
            state.n = length(dataNames);

            for i = 1:state.n
                dataName = dataNames{i};
                dataset = database.getDataset(dataName);
                g = dataset.Data;
                ax = axes(i);

                cla(ax);
                hold(ax, 'on');

                imagesc(ax, g.x{1}, g.x{2}, g.u');

                axis(ax, 'xy', 'equal', 'tight');
                xlabel(ax, 'x', 'FontSize', obj.DefaultAxisLabelFontSize);
                ylabel(ax, 'y', 'FontSize', obj.DefaultAxisLabelFontSize);

                colormap(ax, obj.ColorMap);
                colorbar(ax);
                obj.setColorLimits(ax, g.u);

                titleText = sprintf('%s: %s', ...
                    strrep(dataName, '_', '-'), variable);
                title(ax, titleText);

                hold(ax, 'off');
            end

            if length(axes) > 1
                obj.syncColorLimits(axes, database);
            end
        end
    end
end