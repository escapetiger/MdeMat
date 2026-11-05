classdef Strategy3d < profilers.visual.Strategy
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
    % Examples:
    %   % Create 3D strategy and render volumetric data
    %   strategy = Strategy3d();
    %   strategy.render(axesArray, datasetStruct, 'temperature3', 1);
    %
    % See Also:
    %   profilers.visual.Strategy, profilers.visual.Strategy2d,
    %   profilers.visual.StrategySlice1d

    methods
        function obj = render(obj, axes, varDataset, varName, offset)
            % RENDER Render 3D slice plots with synchronized styling.
            %
            %   render(obj, axes, varDataset, varName, offset) creates
            %   3D slice visualizations for each dataset using orthogonal
            %   cross-sections through the volume center. Each dataset
            %   receives its own subplot with consistent 3D formatting
            %   and synchronized color limits.
            %
            % Inputs:
            %   obj - The Strategy3d object
            %   axes - Array of axes handles for rendering
            %   varDataset - Structure mapping dataset names to data structures
            %   varName - Variable name for plot titles
            %   offset - Starting axis offset for multi-panel layouts (scalar)
            %
            % Outputs:
            %   obj - The Strategy3d object (for method chaining)

            varAxes = [];
            dataNames = fieldnames(varDataset);

            %< Render each dataset on separate axes
            for i = 1:length(dataNames)
                dataName = dataNames{i};
                g = varDataset.(dataName);
                ax = axes(offset+i-1);
                varAxes = [varAxes, ax];

                %< Prepare axes for 3D rendering
                cla(ax);
                hold(ax, 'on');

                %< Calculate slice positions at volume center
                z1 = (min(g.x1) + max(g.x1)) / 2; %< X-direction slice position
                z2 = (min(g.x2) + max(g.x2)) / 2; %< Y-direction slice position
                z3 = (min(g.x3) + max(g.x3)) / 2; %< Z-direction slice position
                
                %< Permute data for correct slice orientation
                u = permute(g.u, [2, 1, 3]);
                
                %< Create orthogonal slice plots through volume center
                h = slice(ax, g.x1, g.x2, g.x3, u, z1, z2, z3);
                set(h, 'EdgeColor', 'none'); %< Remove edge lines for cleaner appearance

                %< Apply standard 3D axis labels with consistent font sizing
                xlabel(ax, 'x', 'FontSize', obj.DEFAULT_AXIS_LABEL_FONT_SIZE);
                ylabel(ax, 'y', 'FontSize', obj.DEFAULT_AXIS_LABEL_FONT_SIZE);
                zlabel(ax, 'z', 'FontSize', obj.DEFAULT_AXIS_LABEL_FONT_SIZE);
                
                %< Apply colormap and add colorbar
                colormap(ax, obj.cmap);
                colorbar(ax);
                
                %< Set equal aspect ratio for proper 3D visualization
                axis(ax, 'equal');
                
                %< Set color limits based on data range
                obj.setColorLimits(ax, u);

                %< Create descriptive title with dataset and variable info
                titleText = sprintf('%s: %s', strrep(dataName, '_', '-'), varName);
                title(ax, titleText);
                
                %< Set 3D viewing angle
                view(ax, 3);

                hold(ax, 'off');
            end

            %< Synchronize color limits across multiple datasets for comparison
            if length(varAxes) > 1
                obj.syncColorLimits(varAxes, varDataset);
            end
        end
    end
end