classdef Strategy3d < physics.visual.Strategy
    % STRATEGY3D Three-dimensional slice plot rendering strategy.
    %
    %   Strategy3d extends the base Strategy class to provide specialized
    %   rendering for 3D volumetric data using slice plots. The strategy
    %   creates cross-sectional views through the volume at the midpoint
    %   of each spatial dimension, providing insight into the 3D data
    %   structure.
    %
    %   Each dataset is rendered with orthogonal slices and synchronized
    %   color limits when multiple datasets are present. The visualization
    %   uses 3D perspective with proper axis labeling and colorbar display.
    %
    % See also:
    %   physics.visual.Strategy, physics.visual.Strategy2d,
    %   physics.visual.StrategySlice1d

    properties (Constant)
        NDims = 3 % Number of spatial dimensions
    end

    methods
        function obj = render(obj, axes, database, variable)
            % RENDER Render 3D slice plots with synchronized styling.
            %
            %   obj = render(obj, axes, database, variable) creates 3D slice
            %   visualizations for each dataset using orthogonal cross-sections
            %   through the volume center. Each dataset receives its own subplot
            %   with consistent 3D formatting and synchronized color limits.

            arguments
                obj physics.visual.Strategy3d
                axes matlab.graphics.axis.Axes
                database physics.visual.Database
                variable {mustBeTextScalar}
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

                z1 = (min(g.x{1}) + max(g.x{1})) / 2;
                z2 = (min(g.x{2}) + max(g.x{2})) / 2;
                z3 = (min(g.x{3}) + max(g.x{3})) / 2;
                
                u = permute(g.u, [2, 1, 3]);
                
                h = slice(ax, g.x{1}, g.x{2}, g.x{3}, u, z1, z2, z3);
                set(h, 'EdgeColor', 'none');

                xlabel(ax, 'x', 'FontSize', obj.DefaultAxisLabelFontSize);
                ylabel(ax, 'y', 'FontSize', obj.DefaultAxisLabelFontSize);
                zlabel(ax, 'z', 'FontSize', obj.DefaultAxisLabelFontSize);
                
                colormap(ax, obj.ColorMap);
                colorbar(ax);
                axis(ax, 'equal', 'tight');
                obj.setColorLimits(ax, u);

                titleText = sprintf('%s: %s', strrep(dataName, '_', '-'), variable);
                title(ax, titleText);
                view(ax, 3);

                hold(ax, 'off');
            end

            if length(axes) > 1
                obj.syncColorLimits(axes, database);
            end
        end
    end
end