classdef StrategySlice1d < profilers.visual.Strategy
    % STRATEGYSLICE1D One-dimensional diagonal slice rendering strategy.
    %
    %   StrategySlice1d extends the base Strategy class to provide
    %   specialized rendering for 1D slices through multi-dimensional
    %   data. The strategy extracts diagonal slices from 2D or 3D datasets
    %   and renders them as line plots for comparative analysis.
    %
    %   The diagonal slice follows the path from one corner to the opposite
    %   corner of the domain, providing insight into the data variation
    %   along the main diagonal. This is particularly useful for analyzing
    %   symmetry properties and overall trends in multi-dimensional data.
    %
    % Examples:
    %   % Create slice strategy and render diagonal cuts
    %   strategy = StrategySlice1d();
    %   strategy.render(axesArray, datasetStruct, 'velocity2', 1);
    %
    % See Also:
    %   profilers.visual.Strategy, profilers.visual.Strategy1d,
    %   profilers.visual.Strategy2d, profilers.visual.Strategy3d

    methods
        function obj = render(obj, axes, varDataset, varName, offset)
            % RENDER Render diagonal slice plots for multi-dimensional
            % data.
            %
            %   obj = render(obj, axes, varDataset, varName, offset)
            %   extracts diagonal slices from 2D or 3D datasets and renders
            %   them as 1D line plots. The slice follows the main diagonal
            %   of the domain, scaled appropriately for each dimension.
            %
            % Inputs:
            %   obj - The StrategySlice1d object
            %   axes - Array of axes handles for rendering
            %   varDataset - Structure mapping dataset names to data structures
            %   varName - Variable name for axis labeling and plot title
            %   offset - Axis offset for multi-panel layouts (scalar)
            %
            % Outputs:
            %   obj - The StrategySlice1d object

            dataNames = fieldnames(varDataset);
            n = length(dataNames);
            if n < 1
                return; %< No data to render
            end

            ax = axes(offset);
            
            %< Prepare axes for slice rendering
            cla(ax);
            hold(ax, 'on');

            %< Initialize axis limits and legend storage
            [xmin, xmax, ymin, ymax] = deal(inf, -inf, inf, -inf);
            legends = cell(1, n);

            %< Track style indices separately for lines and scatter plots
            lineIdx = 1;
            scatterIdx = 1;
            
            %< Process each dataset
            for i = 1:n
                dataName = dataNames{i};
                g = varDataset.(dataName);

                %< Create diagonal slice coordinate vector
                r = linspace(min(g.x1), max(g.x1), 64); %< 64 points along diagonal
                
                %< Extract diagonal slice based on dimensionality
                if g.nDims == 2
                    %< 2D case: diagonal from (min,min) to (max,max)
                    R = sqrt(2); %< Diagonal scaling factor
                    s = interp2(g.x1, g.x2, g.u.', r/R, r/R);
                elseif g.nDims == 3
                    %< 3D case: space diagonal from corner to corner
                    R = sqrt(3); %< 3D diagonal scaling factor
                    u = permute(g.u, [2, 1, 3]); %< Adjust array orientation
                    s = interp3(g.x1, g.x2, g.x3, u, r/R, r/R, r/R);
                else
                    %< Unsupported dimensionality - skip this dataset
                    continue;
                end

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
                
                %< Render slice data with selected styling
                plot(ax, r, s, style{:}, 'DisplayName', dataName);

                %< Update axis limits based on current slice
                xmin = min(xmin, min(r(:)));
                xmax = max(xmax, max(r(:)));
                if ~isempty(s) && any(~isnan(s))
                    ymin = min(ymin, min(s(:)));
                    ymax = max(ymax, max(s(:)));
                end

                %< Prepare legend entry with underscore replacement
                legends{i} = strrep(dataName, '_', '-');
            end

            %< Handle degenerate case where all data is constant
            if ymin == ymax
                ymin = ymin - 1e-8;
                ymax = ymax + 1e-8;
            end

            %< Apply final axis formatting with descriptive title
            obj.finalizeAxis(ax, xmin, xmax, ymin, ymax, 'x', varName, legends);
            title(ax, 'Diagonal slice');
        end
    end
end