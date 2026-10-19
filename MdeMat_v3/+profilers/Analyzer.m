classdef Analyzer < handle
    % ANALYZER Comprehensive error analysis and convergence study.
    %
    %   Analyzer provides sophisticated error analysis capabilities for
    %   numerical methods, including absolute error computation against
    %   exact solutions, relative error analysis using Richardson
    %   extrapolation, and convergence analysis with visualization.
    %
    %   This unified class handles the complete workflow from error
    %   computation through convergence table generation and plotting,
    %   eliminating the need for separate convergence profiler management.
    %
    % Examples:
    %   % Create analyzer for convergence study
    %   analyzer = Analyzer(4);
    %   analyzer.setDensity([10, 10]);
    %   analyzer.setReductions({'L1', 'L2'});
    %   analyzer.setComponents(struct('u', 1, 'v', 2));
    %   analyzer.addExact('u', @(x) sin(x(1,:)));
    %   analyzer.initialize();
    %
    %   % Perform multi-level error analysis
    %   for iLevel = 1:4
    %       analyzer.addLevel(mesh);
    %       analyzer.addAbsolute(space, dofs);
    %   end
    %
    %   % Generate results and visualization
    %   results = analyzer.analyze();
    %   analyzer.draw(errors, names, rates);
    %
    % See also:
    %   profilers.metrics.MeshSpaceAbsoluteMetric,
    %   profilers.metrics.GridSpaceRichardsonMetric

    properties
        % Convergence study configuration
        nLevels % Number of refinement levels (positive integer)
        resolutions % (nLevels×1) mesh resolution vector
        descriptions % (nLevels×1) cell array of mesh level descriptions

        % Error analysis configuration
        density % (1×nDims) error evaluation points per element per dimension
        reductions % Cell array of error reduction types ('L1', 'L2', 'Lx')
        components % Structure mapping field names to component indices
        exacts % Structure of exact solution function handles

        % Data storage
        errors % Structure storing computed errors for each level
        currentLevel % Current level count for error storage
        tiledLayout % Handle to tiled layout for plotting
    end

    properties (Dependent)
        hasExact % True if any exact solution is available
        isEnabled % True if error analysis is enabled
        nComponents % Number of components for each field
        nReductions % Number of reduction types
    end

    methods
        function obj = Analyzer(nLevels)
            % ANALYZER Constructor for unified error and convergence analyzer.
            %
            %   obj = Analyzer(nLevels) creates an analyzer configured for
            %   error analysis and convergence study with the specified
            %   number of refinement levels.
            %
            % Inputs:
            %   nLevels - Number of refinement levels (positive integer)
            %
            % Outputs:
            %   obj - Constructed Analyzer object

            obj.nLevels = nLevels;
            obj.resolutions = zeros(nLevels, 1);
            obj.descriptions = cell(nLevels, 1);
            obj.density = 0;
            obj.reductions = {};
            obj.components = struct();
            obj.exacts = struct();
            obj.errors = struct();
            obj.currentLevel = 1;
        end

        function obj = initialize(obj)
            % INITIALIZE Initialize error storage structure.
            %
            %   obj = initialize(obj) prepares the error storage structure
            %   based on the configured components. This method must be
            %   called after setting components and before error computation.
            %
            % Inputs:
            %   obj - The Analyzer object
            %
            % Outputs:
            %   obj - The Analyzer object

            fieldNames = fieldnames(obj.components);
            nFields = length(fieldNames);
            args = cell(2, nFields);
            args(1, :) = fieldNames;
            args(2, :) = {{cell(obj.nLevels, 1)}};
            obj.errors = struct(args{:});
            obj.currentLevel = 1;
        end

        function obj = setDensity(obj, density)
            % SETDENSITY Set the number of error evaluation points per element.
            %
            %   obj = setDensity(obj, density) sets the number of points
            %   per element per dimension for error evaluation sampling.
            %
            % Inputs:
            %   obj - The Analyzer object
            %   density - Number of points per element per dimension (vector)
            %
            % Outputs:
            %   obj - The Analyzer object

            obj.density = density;
        end

        function obj = setComponents(obj, components)
            % SETCOMPONENTS Set the solution components to analyze.
            %
            %   obj = setComponents(obj, components) sets the indices of
            %   solution components to analyze for each field.
            %
            % Inputs:
            %   obj - The Analyzer object
            %   components - Structure mapping field names to component indices
            %
            % Outputs:
            %   obj - The Analyzer object

            obj.components = components;
        end

        function obj = setReductions(obj, reductions)
            % SETREDUCTIONS Set the error reduction types.
            %
            %   obj = setReductions(obj, reductions) sets the types of
            %   error norms to compute (L1, L2, Lx, etc.).
            %
            % Inputs:
            %   obj - The Analyzer object
            %   reductions - Cell array of reduction types
            %
            % Outputs:
            %   obj - The Analyzer object

            obj.reductions = reductions;
        end

        function obj = addExact(obj, name, func)
            % ADDEXACT Add exact solution function.
            %
            %   obj = addExact(obj, name, func) registers an exact
            %   solution function for the specified field name.
            %
            % Inputs:
            %   obj - The Analyzer object
            %   name - Solution field name (string or char array)
            %   func - Exact solution function handle
            %
            % Outputs:
            %   obj - The Analyzer object

            obj.exacts.(name) = func;
        end

        function obj = addLevel(obj, mesh)
            % ADDLEVEL Configure analyzer for current refinement level.
            %
            %   obj = addLevel(obj, mesh) configures the analyzer for the
            %   current refinement level using mesh properties. This method
            %   extracts the mesh resolution and creates descriptive
            %   strings for convergence analysis.
            %
            % Inputs:
            %   obj - The Analyzer object
            %   mesh - Mesh object with resolution and measure properties
            %
            % Outputs:
            %   obj - The Analyzer object

            h = mesh.computeMeasure();
            hStr = strrep(num2str(mesh.resolution), '  ', 'x');
            obj.resolutions(obj.currentLevel) = h;
            obj.descriptions{obj.currentLevel} = hStr;
        end

        function obj = addAbsolute(obj, space, dofs, varargin)
            % ADDABSOLUTE Compute absolute error against exact solutions.
            %
            %   obj = addAbsolute(obj, space, dofs) computes absolute
            %   errors by comparing numerical solutions with exact
            %   analytical solutions. The method evaluates both numerical
            %   and exact solutions at the same sample points and computes
            %   specified error metrics.
            %
            % Inputs:
            %   obj - The Analyzer object
            %   space - Function space for solution evaluation
            %   dofs - Structure containing solution coefficient dofs
            %   varargin - Additional arguments for exact solution evaluation
            %
            % Outputs:
            %   obj - The Analyzer object

            metric = profilers.metrics.MeshSpaceAbsoluteMetric(obj.reductions, space);
            metric.setPoints(obj.density);

            X = metric.nodes;
            fieldNames = fieldnames(dofs);
            for iField = 1:length(fieldNames)
                fieldName = fieldNames{iField};
                component = obj.components.(fieldName);
                C = dofs.(fieldName);
                f = obj.exacts.(fieldName);
                U = space.evaluate([], X, C);
                U = U(:, component);
                U0 = space.functionEvaluate([], X, f, varargin{:});
                U0 = U0(:, component);
                R = metric.evaluate(U, U0);
                E = zeros(size(R, 1), 1);
                for j = 1:length(obj.reductions)
                    switch lower(obj.reductions{j})
                        case 'l1'
                            E(j) = norm(R(j, :), 1);
                        case 'l2'
                            E(j) = norm(R(j, :), 2);
                        case 'lx'
                            E(j) = norm(R(j, :), inf);
                    end
                end
                obj.errors.(fieldName){obj.currentLevel} = E;
            end
            obj.currentLevel = obj.currentLevel + 1;
        end

        function obj = addRelative(obj, fSpace, fDofs, cSpace, cDofs)
            % ADDRELATIVE Compute relative error using Richardson
            % extrapolation.
            %
            %   obj = addRelative(obj, fSpace, fDofs, cSpace, cDofs)
            %   computes relative errors using Richardson extrapolation
            %   between fine and coarse grid solutions. This method is
            %   useful when exact solutions are not available.
            %
            % Inputs:
            %   obj - The Analyzer object
            %   fSpace - Fine grid function space
            %   fDofs - Structure containing fine grid solution dofs
            %   cSpace - Coarse grid function space
            %   cDofs - Structure containing coarse grid solution dofs
            %
            % Outputs:
            %   obj - The Analyzer object

            metric = profilers.metrics.GridSpaceRichardsonMetric( ...
                obj.reductions, fSpace, cSpace);
            metric.setPoints(obj.density);

            X = metric.nodes{1};
            X0 = metric.nodes{2};

            fieldNames = fieldnames(fDofs);
            for iField = 1:length(fieldNames)
                fieldName = fieldNames{iField};
                component = obj.components.(fieldName);
                iNestDofs = metric.iNestDofs; % Fixed Hungarian notation
                C = fDofs.(fieldName)(iNestDofs, :);
                C0 = cDofs.(fieldName);
                U = fSpace.evaluate([], X, C);
                U = U(:, component);
                U0 = cSpace.evaluate([], X0, C0);
                U0 = U0(:, component);
                R = metric.evaluate(U, U0);
                E = zeros(size(R, 1), 1);
                for j = 1:length(obj.reductions)
                    switch lower(obj.reductions{j})
                        case 'l1'
                            E(j) = norm(R(j, :), 1);
                        case 'l2'
                            E(j) = norm(R(j, :), 2);
                        case 'lx'
                            E(j) = norm(R(j, :), inf);
                    end
                end
                obj.errors.(fieldName){obj.currentLevel} = E;
            end
            obj.currentLevel = obj.currentLevel + 1;
        end

        function tables = analyze(obj, verbose)
            % ANALYZE Generate comprehensive convergence analysis results.
            %
            %   tables = analyze(obj) processes all collected error data
            %   and generates structured results including convergence
            %   tables for each field component and error reduction type.
            %   The results can be used for convergence verification and
            %   method validation.
            %
            % Inputs:
            %   obj - The Analyzer object
            %   verbose - Verbose flag (scalar)
            %
            % Outputs:
            %   tables - Structure containing convergence tables

            tables = struct();
            fieldNames = fieldnames(obj.errors);
            for iField = 1:length(fieldNames)
                fieldName = fieldNames{iField};
                component = obj.components.(fieldName);
                if isempty(component), continue; end
                E = cell2mat(obj.errors.(fieldName));
                for iComponent = 1:size(E, 2)
                    Ej = reshape(E(:, iComponent), obj.nReductions, obj.nLevels).';
                    name = sprintf('%s%d', fieldName, component(iComponent));
                    tables.(name) = obj.analyzeImpl(Ej, obj.reductions);
                    if verbose > 0
                        disp(tables.(name));
                    end
                end
            end
        end

        function obj = draw(obj, varargin)
            % DRAW Create convergence visualization plots.
            %
            %   obj = draw(obj) creates log-log convergence plots using
            %   stored error data from previous error computations.
            %
            %   obj = draw(obj, rates) creates plots with optional
            %   theoretical convergence rates for comparison.
            %
            % Inputs:
            %   obj - The Analyzer object
            %   rates - Theoretical convergence rates (optional)
            %
            % Outputs:
            %   obj - The Analyzer object

            rates = [];
            if nargin >= 2
                rates = varargin{1};
            end

            fieldNames = fieldnames(obj.errors);
            if isempty(fieldNames)
                core.except.verify(0, 'NoData', 'No error data available for plotting.');
                return;
            end

            %< Validate that we have resolution data
            core.except.assert(all(obj.resolutions > 0), 'InvalidResolutions', ...
                'All resolutions must be positive for log-log plot.');

            %< Calculate total number of components for layout
            totalComponents = sum(obj.nComponents);

            %< Create tiled layout
            obj.setupTiledLayout(totalComponents);

            %< Process each field and component
            componentIndex = 1;
            rateIndex = 1;

            for iField = 1:length(fieldNames)
                fieldName = fieldNames{iField};
                component = obj.components.(fieldName);
                E = reshape(cell2mat(obj.errors.(fieldName)), [], obj.nLevels);

                for iComponent = 1:length(component)
                    %< Extract errors for this component across all reduction types
                    componentErrors = zeros(obj.nLevels, obj.nReductions);
                    for iReduction = 1:obj.nReductions
                        errorIndex = (iReduction - 1) * length(component) + iComponent;
                        componentErrors(:, iReduction) = E(errorIndex, :);
                    end

                    %< Extract rates for this component if provided
                    componentRates = [];
                    if ~isempty(rates)
                        endIndex = min(rateIndex+obj.nReductions-1, length(rates));
                        componentRates = rates(rateIndex:endIndex);
                        rateIndex = endIndex + 1;
                    end

                    %< Create subplot for this component
                    componentTitle = sprintf('%s%d', fieldName, component(iComponent));
                    obj.drawImpl(componentErrors, obj.reductions, componentRates, componentTitle);

                    componentIndex = componentIndex + 1;
                end
            end
        end
    end

    methods
        function TF = get.isEnabled(obj)
            % GET.ISENABLED Check if error analysis is enabled.

            TF = any(obj.density > 0) && sum(obj.nComponents) > 0;
        end

        function TF = get.hasExact(obj)
            % GET.HASEXACT Check if exact solutions are available.

            TF = ~isempty(fieldnames(obj.exacts));
        end

        function n = get.nComponents(obj)
            % GET.NCOMPONENTS Get number of components for each field.

            n = cellfun(@(f) length(obj.components.(f)), fieldnames(obj.components));
        end

        function n = get.nReductions(obj)
            % GET.NREDUCTIONS Get number of reduction types.

            n = length(obj.reductions);
        end
    end

    methods (Access = protected)
        function obj = setupTiledLayout(obj, nComponents)
            % SETUPTILEDLAYOUT Create optimally sized tiled layout.
            %
            %   obj = setupTiledLayout(obj, nComponents) creates a figure
            %   with tiled layout optimized for the given number of components.
            %
            % Inputs:
            %   obj - The Analyzer object
            %   nComponents - Number of components to display
            %
            % Outputs:
            %   obj - The Analyzer object

            %< Determine optimal layout dimensions
            if nComponents == 1
                nRows = 1;
                nCols = 1;
            elseif nComponents == 2
                nRows = 1;
                nCols = 2;
            elseif nComponents <= 4
                nRows = 2;
                nCols = 2;
            elseif nComponents <= 6
                nRows = 2;
                nCols = 3;
            else
                nRows = ceil(sqrt(nComponents));
                nCols = ceil(nComponents/nRows);
            end

            %< Create figure with tiled layout
            figure;
            obj.tiledLayout = tiledlayout(nRows, nCols, 'TileSpacing', 'compact', ...
                'Padding', 'compact');
            title(obj.tiledLayout, 'Error Convergence Analysis', ...
                'FontSize', 14, 'FontWeight', 'bold');
        end

        function T = analyzeImpl(obj, errors, names)
            % ANALYZEIMPL Build formatted convergence table with error orders.
            %
            %   T = analyzeImpl(obj, errors, names) creates a formatted
            %   table containing errors and computed convergence orders for
            %   each refinement level. Convergence orders are calculated
            %   using the formula:
            %
            %   order = log(E_h / E_{h/2}) / log(hRatio)
            %
            %   where E_h is the error at resolution h.
            %
            % Inputs:
            %   obj - The Analyzer object
            %   errors - Numerical error matrix (nLevels × nTypes)
            %   names - Cell array of error type names (1 × nTypes)
            %
            % Outputs:
            %   T - MATLAB table

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
            T = cell2table(tableData, 'VariableNames', columnNames);
        end

        function obj = drawImpl(obj, errors, names, rates, plotTitle)
            % DRAWIMPL Create single subplot in tiled layout.
            %
            %   obj = drawImpl(obj, errors, names, rates, plotTitle) creates
            %   a single subplot showing convergence behavior for one component.
            %
            % Inputs:
            %   obj - The Analyzer object
            %   errors - (nLevels×nTypes) numerical error matrix for one component
            %   names - (1×nTypes) cell array of error type names
            %   rates - (1×nTypes) theoretical convergence rates (optional)
            %   plotTitle - Title for this subplot
            %
            % Outputs:
            %   obj - The Analyzer object

            nexttile;

            [~, nTypes] = size(errors);
            plotColors = lines(nTypes);

            hold on;

            %< Plot actual errors for this component
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

            %< Format subplot
            set(gca, 'XScale', 'log', 'YScale', 'log');
            xlabel('Resolution');
            ylabel('Error');
            title(plotTitle);
            legend('show', 'Location', 'best');
            grid on;
            hold off;
        end
    end
end