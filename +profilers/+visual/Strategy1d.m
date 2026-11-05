classdef Strategy1d < profilers.visual.Strategy
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
    %   can be configured through the base Strategy class methods.
    %
    % Examples:
    %   % Create and configure 1D strategy
    %   strategy = Strategy1d();
    %   strategy.setLogScale('y', true);
    %   
    %   % Render data to axes
    %   strategy.render(axesArray, datasetStruct, 'u1', 1);
    %
    % See Also:
    %   profilers.visual.Strategy, profilers.visual.Strategy2d,
    %   profilers.visual.StrategySlice1d

    methods        
        function obj = render(obj, axes, varDataset, varName, offset)
            % RENDER Render 1D line plot with automatic dataset styling.
            %
            %   obj = render(obj, axes, varDataset, varName, offset)
            %   creates a 1D line plot on the specified axes, automatically
            %   selecting line styles for analytical data (type < 2) and
            %   scatter styles for numerical data (type >= 2). Multiple
            %   datasets are overlaid with distinct styling and legend
            %   entries.
            %
            % Inputs:
            %   obj - The Strategy1d object
            %   axes - Array of axes handles for rendering
            %   varDataset - Structure mapping dataset names to data structures
            %   varName - Variable name for axis labeling
            %   offset - Axis offset for multi-panel layouts (scalar)
            %
            % Outputs:
            %   obj - The Strategy1d object

            dataNames = fieldnames(varDataset);
            n = length(dataNames);
            ax = axes(offset);
            
            %< Prepare axes for rendering
            cla(ax);
            hold(ax, 'on');

            %< Initialize axis limits and legend storage
            [xmin, xmax, umin, umax] = deal(inf, -inf, inf, -inf);
            legends = cell(1, n);

            %< Track style indices separately for lines and scatter plots
            lineIdx = 1;
            scatterIdx = 1;
            
            %< Process each dataset
            for i = 1:n
                dataName = dataNames{i};
                g = varDataset.(dataName);

                %< Update axis limits based on current dataset
                xmin = min(xmin, min(g.x1));
                xmax = max(xmax, max(g.x1));
                umin = min(umin, min(g.u));
                umax = max(umax, max(g.u));

                %< Select styling based on data type
                if g.type < 2
                    %< Analytical/exact data: use line styling
                    style = obj.getDefaultLineStyle(i, lineIdx);
                    lineIdx = lineIdx + 1;
                else
                    %< Numerical data: use scatter styling
                    style = obj.getDefaultScatterStyle(i, scatterIdx);
                    scatterIdx = scatterIdx + 1;
                end
                
                %< Render data with selected styling
                plot(ax, g.x1(:).', g.u(:).', style{:});
                
                %< Prepare legend entry with underscore replacement
                legends{i} = strrep(dataName, '_', '-');
            end

            %< Apply final axis formatting and limits
            obj.finalizeAxis(ax, xmin, xmax, umin, umax, 'x', varName, legends);
        end
    end
end