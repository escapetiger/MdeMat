classdef ConvergenceProfiler < handle
    % CONVERGENCEPROFILER Convergence analysis and visualization for
    % numerical methods.
    %
    %   ConvergenceProfiler provides comprehensive functionality to analyze
    %   the convergence behavior of numerical methods by computing
    %   convergence orders and generating formatted error tables and plots.
    %   This class is essential for validating numerical schemes and
    %   understanding their theoretical convergence properties.
    %
    %   The profiler stores analysis configuration using MATLAB table
    %   structures and provides methods for building formatted error tables
    %   with convergence orders and creating log-log convergence
    %   visualizations that clearly show the relationship between mesh
    %   resolution and numerical errors.
    %
    % Examples:
    %   % Create profiler for 4 refinement levels
    %   profiler = ConvergenceProfiler(4);
    %   
    %   % Set mesh resolutions and descriptions
    %   profiler.resolutions = [0.1; 0.05; 0.025; 0.0125];
    %   profiler.descriptions = {'h=0.1'; 'h=0.05'; 'h=0.025'; 'h=0.0125'};
    %   
    %   % Build error table with convergence orders
    %   errors = [1e-2, 8e-3; 5e-3, 4e-3; 2.5e-3, 2e-3; 1.25e-3, 1e-3];
    %   names = {'L1', 'L2'};
    %   errorTable = profiler.buildTable(errors, names);
    %   
    %   % Create convergence plot with theoretical rates
    %   profiler.draw(errors, names, [1.0, 2.0]);
    %
    % See also:
    %   table, loglog, profilers.analysis.Analyzer

    properties
        nLevels     % Number of refinement levels (positive integer)
        resolutions % (nLevels×1) mesh resolution vector (positive values)
        descriptions % (nLevels×1) cell array of mesh level descriptions
    end

    methods
        function obj = ConvergenceProfiler(nLevels)
            % CONVERGENCEPROFILER Constructor for ConvergenceProfiler.
            %
            %   obj = ConvergenceProfiler(nLevels) creates a convergence
            %   profiler configured for the specified number of refinement
            %   levels. The profiler initializes storage for resolutions
            %   and descriptions that must be set before analysis.
            %
            % Inputs:
            %   nLevels - Number of refinement levels (positive integer)
            %
            % Outputs:
            %   obj - Constructed ConvergenceProfiler object
            %
            % Examples:
            %   % Create profiler for standard convergence study
            %   profiler = ConvergenceProfiler(5);
            %
            %   % Create profiler for quick verification
            %   profiler = ConvergenceProfiler(3);
            
            obj.nLevels = nLevels;
            obj.resolutions = zeros(nLevels, 1);
            obj.descriptions = cell(nLevels, 1);
        end

        function result = buildTable(obj, errors, names)
            % BUILDTABLE Build formatted convergence table with error
            % orders.
            %
            %   result = buildTable(obj, errors, names) creates a formatted
            %   table containing errors and computed convergence orders for
            %   each refinement level. Convergence orders are calculated
            %   using the formula: 
            %
            %   \f[
            %       order = \log(\frac{E_h}{E_{h/2}}) / \log(2) 
            %   \f]
            %
            %   where \f$E_h\f$ is the error at resolution \f$h\f$.
            %
            % Inputs:
            %   obj - The ConvergenceProfiler object
            %   errors - (nLevels×nTypes) numerical error matrix
            %   names - (1×nTypes) cell array of error type names (optional)
            %
            % Outputs:
            %   result - MATLAB table

            [~, nTypes] = size(errors);
            
            %< Default names if not provided
            if nargin < 3 || isempty(names)
                names = cell(1, nTypes);
                for i = 1:nTypes
                    names{i} = sprintf('Error %d', i);
                end
            end
                    
            %< Validate inputs
            core.except.assert(length(names) == nTypes, 'SizeMismatch', ...
                'Number of names (%d) must match number of error types (%d).', ...
                length(names), nTypes);
            
            %< Replace empty descriptions with default names
            for i = 1:obj.nLevels
                if isempty(obj.descriptions{i})
                    obj.descriptions{i} = sprintf('Level %d', i);
                end
            end

            %< Compute convergence orders
            orders = zeros(size(errors));
            for j = 1:nTypes
                for k = 2:obj.nLevels
                    hRatio = obj.resolutions(k-1) / obj.resolutions(k);
                    if errors(k, j) > 0 && errors(k-1, j) > 0
                        eRatio = errors(k-1, j) / errors(k, j);
                        orders(k, j) = log(eRatio) / log(hRatio);
                    else
                        orders(k, j) = 0;
                    end
                end
            end

            %< Prepare table data
            tableData = cell(obj.nLevels, 2*nTypes+1);
            tableData(:, 1) = obj.descriptions;
            for j = 1:nTypes
                tableData(:, 2*j) = num2cell(errors(:, j));
                tableData(:, 2*j+1) = num2cell(orders(:, j));
            end

            %< Generate column names
            columnNames = cell(1, 2*nTypes+1);
            columnNames{1} = 'Resolution';
            for j = 1:nTypes
                columnNames{2*j} = sprintf('%s_Error', names{j});
                columnNames{2*j+1} = sprintf('%s_Order', names{j});
            end

            %< Create table
            result = cell2table(tableData, 'VariableNames', columnNames);
        end

        function obj = draw(obj, errors, names, rates)
            % DRAW Create log-log convergence visualization.
            %
            %   obj = draw(obj, errors, names, rates) creates a log-log
            %   plot showing the convergence behavior of numerical errors
            %   as a function of mesh resolution. Optional theoretical
            %   convergence rates can be plotted as reference lines for
            %   comparison.
            %
            % Inputs:
            %   obj - The ConvergenceProfiler object
            %   errors - (nLevels×nTypes) numerical error matrix
            %   names - (1×nTypes) cell array of error type names (optional)
            %   rates - (1×nTypes) theoretical convergence rates (optional)
            %
            % Outputs:
            %   obj - The ConvergenceProfiler object

            [~, nTypes] = size(errors);
            
            %< Default names if not provided
            if nargin < 3 || isempty(names)
                names = cell(1, nTypes);
                for i = 1:nTypes
                    names{i} = sprintf('Error %d', i);
                end
            end
            
            %< Default rates if not provided
            if nargin < 4, rates = []; end

            %< Validate inputs
            core.except.assert(length(names) == nTypes, 'SizeMismatch', ...
                'Number of names (%d) must match number of error types (%d).', ...
                length(names), nTypes);
            core.except.assert(all(obj.resolutions > 0), 'InvalidResolutions', ...
                'All resolutions must be positive for log-log plot.');

            figure;
            hold on;
            plotColors = lines(nTypes);

            %< Plot actual errors
            for j = 1:nTypes
                loglog(obj.resolutions, errors(:, j), '-o', ...
                    'Color', plotColors(j, :), ...
                    'DisplayName', sprintf('%s Error', names{j}), ...
                    'LineWidth', 2, 'MarkerSize', 6);

                %< Plot theoretical convergence rates if provided
                if ~isempty(rates) && length(rates) >= j && rates(j) > 0
                    refErrors = errors(1, j) * (obj.resolutions / obj.resolutions(1)).^(rates(j));
                    loglog(obj.resolutions, refErrors, '--', ...
                        'Color', plotColors(j, :), ...
                        'DisplayName', sprintf('%s Rate %g', names{j}, rates(j)), ...
                        'LineWidth', 1.5);
                end
            end

            %< Format plot
            set(gca, 'XScale', 'log', 'YScale', 'log');
            xlabel('Resolution (mesh size)');
            ylabel('Error');
            legend('show', 'Location', 'best');
            grid on;
            title('Error Convergence Analysis');
            hold off;
        end
    end
end