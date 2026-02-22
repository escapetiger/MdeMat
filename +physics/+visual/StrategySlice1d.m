classdef StrategySlice1d < physics.visual.Strategy1d
    % STRATEGYSLICE1D One-dimensional diagonal slice rendering strategy.
    %
    %   StrategySlice1d extends the Strategy1d class to provide specialized
    %   rendering for 1D slices through multi-dimensional data. The strategy
    %   extracts diagonal slices from 2D or 3D datasets and renders them as
    %   line plots for comparative analysis.
    %
    %   The diagonal slice follows the path from one corner to the opposite
    %   corner of the domain, providing insight into the data variation along
    %   the main diagonal. This is particularly useful for analyzing symmetry
    %   properties and overall trends in multi-dimensional data.
    %
    % See also:
    %   physics.visual.Strategy1d, physics.visual.Strategy2d,
    %   physics.visual.Strategy3d

    methods
        function obj = render(obj, ax, database, variable)
            % RENDER Render diagonal slice plots for multi-dimensional data.
            %
            %   obj = render(obj, ax, database, variable) extracts diagonal slices
            %   from 2D or 3D datasets and renders them as 1D line plots on the
            %   specified @a ax. The slice follows the main diagonal of the domain,
            %   scaled appropriately for each dimension.

            arguments
                obj physics.visual.StrategySlice1d
                ax matlab.graphics.axis.Axes
                database physics.visual.Database
                variable {mustBeTextScalar}
            end

            dataNames = fieldnames(database.Datasets);
            if length(dataNames) < 1
                return;
            end
            
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

                r = linspace(min(g.x{1}), max(g.x{1}), 64);
                nDims = length(g.x);
                
                if nDims == 2
                    R = sqrt(2);
                    s = interp2(g.x{1}, g.x{2}, full(g.u.'), r/R, r/R);
                elseif nDims == 3
                    R = sqrt(3);
                    u = permute(full(g.u), [2, 1, 3]);
                    s = interp3(g.x{1}, g.x{2}, g.x{3}, u, r/R, r/R, r/R);
                else
                    continue;
                end

                if strcmpi(dataset.Type, 'fixed') || strcmpi(dataset.Type, 'exact')
                    style = obj.getDefaultLineStyle(i, lineIdx);
                    lineIdx = lineIdx + 1;
                else
                    style = obj.getDefaultScatterStyle(i, scatterIdx);
                    scatterIdx = scatterIdx + 1;
                end

                plot(ax, r, s, style{:}, 'DisplayName', dataName);

                options.xmin = min(options.xmin, min(r(:)));
                options.xmax = max(options.xmax, max(r(:)));
                if ~isempty(s) && any(~isnan(s))
                    options.ymin = min(options.ymin, min(s(:)));
                    options.ymax = max(options.ymax, max(s(:)));
                end

                options.legends{i} = strrep(dataName, '_', '-');
            end

            tol = 1e-8;
            if options.ymin == options.ymax
                options.ymin = options.ymin - tol;
                options.ymax = options.ymax + tol;
            end

            title(ax, 'Diagonal slice');

            args = namedargs2cell(options);
            obj.formatAxis(ax, args{:});
        end
    end
end