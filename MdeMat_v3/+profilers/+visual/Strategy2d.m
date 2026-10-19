classdef Strategy2d < profilers.visual.Strategy
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
    % Examples:
    %   % Create 2D strategy and render data
    %   strategy = Strategy2d();
    %   strategy.render(axesArray, datasetStruct, 'pressure2', 1);
    %
    % See Also:
    %   profilers.visual.Strategy, profilers.visual.Strategy1d,
    %   profilers.visual.Strategy3d

    methods
        function obj = render(obj, axes, varDataset, varName, offset)
            % RENDER Render 2D surface plots with synchronized styling.
            %
            %   render(obj, axes, varDataset, varName, offset) creates
            %   2D surface visualizations for each dataset using imagesc
            %   with proper coordinate mapping and color synchronization.
            %   Each dataset receives its own subplot with consistent
            %   formatting and colorbar display.
            %
            % Inputs:
            %   obj - The Strategy2d object
            %   axes - Array of axes handles for rendering
            %   varDataset - Structure mapping dataset names to data structures
            %   varName - Variable name for plot titles
            %   offset - Starting axis offset for multi-panel layouts (scalar)
            %
            % Outputs:
            %   obj - The Strategy2d object (for method chaining)

            dataNames = fieldnames(varDataset);
            varAxes = [];
            
            %< Render each dataset on separate axes
            for i = 1:length(dataNames)
                dataName = dataNames{i};
                g = varDataset.(dataName);
                ax = axes(offset + i - 1);
                varAxes = [varAxes, ax];

                %< Prepare axes for 2D rendering
                cla(ax);
                hold(ax, 'on');

                %< Create image plot with proper coordinate mapping
                imagesc(ax, g.x1, g.x2, g.u');
                
                %< Configure axis properties for correct display
                axis(ax, 'xy', 'equal', 'tight');
                
                %< Apply standard axis labels with consistent font sizing
                xlabel(ax, 'x', 'FontSize', obj.DEFAULT_AXIS_LABEL_FONT_SIZE);
                ylabel(ax, 'y', 'FontSize', obj.DEFAULT_AXIS_LABEL_FONT_SIZE);
                
                %< Apply colormap and add colorbar
                colormap(ax, obj.cmap);
                colorbar(ax);

                %< Set color limits based on data range
                obj.setColorLimits(ax, g.u);

                %< Create descriptive title with dataset and variable info
                titleText = sprintf('%s: %s', strrep(dataName, '_', '-'), varName);
                title(ax, titleText);

                hold(ax, 'off');
            end

            %< Synchronize color limits across multiple datasets for comparison
            if length(varAxes) > 1
                obj.syncColorLimits(varAxes, varDataset);
            end
        end
    end
end