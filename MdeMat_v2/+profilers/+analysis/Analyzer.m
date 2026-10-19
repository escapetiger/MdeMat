classdef Analyzer < handle
    % ANALYZER Comprehensive error analysis and metrics computation.
    %
    %   Analyzer provides sophisticated error analysis capabilities for
    %   numerical methods, including absolute error computation against
    %   exact solutions and relative error analysis using Richardson
    %   extrapolation. This class is essential for validating numerical
    %   schemes and quantifying their accuracy.
    %
    %   The analyzer supports multiple error reduction types (L1, L2,
    %   Lx), handles multi-component solutions, and can work with both
    %   exact solutions and grid-based Richardson extrapolation for
    %   problems where analytical solutions are unavailable.
    %
    % Examples:
    %   % Create analyzer for convergence study
    %   analyzer = Analyzer('nLevels', 4, ...
    %                      'density', [10, 10], ...
    %                      'reductions', {'L1', 'L2'}, ...
    %                      'components', struct('u', 1, 'v', 2), ...
    %                      'exacts', struct('u', @(x) sin(x(1,:))));
    %   
    %   % Compute absolute error against exact solution
    %   analyzer.setLevel(mesh);
    %   analyzer.absolute(space, dofs);
    %   
    %   % Compute relative error using Richardson extrapolation
    %   analyzer.relative(fineSpace, fineDofs, coarseSpace, coarseDofs);
    %   
    %   % Generate analysis results
    %   results = analyzer.analyze();
    %
    % See also:
    %   profilers.analysis.ConvergenceProfiler, profilers.metrics.MeshSpaceAbsoluteMetric

    properties
        profiler    % ConvergenceProfiler for managing convergence analysis
        density     % (1×nDims) number of error evaluation points per element per dimension
        reductions  % Cell array of error reduction types ('L1', 'L2', 'Lx')
        components  % Structure mapping field names to component indices
        exacts      % Structure of exact solution function handles
        errors      % Structure storing computed errors for each level
        count       % Current level count for error storage
    end

    properties (Dependent)
        hasExact     % True if any exact solution is available
        isEnabled    % True if error analysis is enabled
        nComponents  % Number of components for each field
    end
    
    methods
        function obj = Analyzer(varargin)
            % ANALYZER Constructor for Analyzer.
            %
            %   obj = Analyzer('Parameter', Value, ...) creates an analyzer
            %   with specified configuration for error analysis. The
            %   analyzer supports various error metrics and can handle
            %   multi-component solutions with optional exact solutions.
            %
            % Inputs:
            %   varargin - Input arguments
            %<   nLevels - Number of refinement levels (positive integer)
            %<   density - (1×nDims) error sampling density per element (optional)
            %<   reductions - Cell array of error reduction types (optional)
            %<   components - Structure mapping field names to indices (optional)
            %<   exacts - Structure of exact solution functions (optional)
            %
            % Outputs:
            %   obj - Constructed Analyzer object

            p = inputParser;
            addParameter(p, 'nLevels', [], @(x) isnumeric(x) && isscalar(x));
            addParameter(p, 'reductions', [], @(x) iscell(x));
            addParameter(p, 'density', [], @(x) isnumeric(x));
            addParameter(p, 'components', [], @(x) isstruct(x));
            parse(p, varargin{:});
            
            obj.profiler = profilers.analysis.ConvergenceProfiler(p.Results.nLevels);
            obj.density = p.Results.density;
            obj.reductions = p.Results.reductions;
            obj.components = p.Results.components;
            obj.exacts = struct();

            %< Initialize error storage structure
            fieldNames = fieldnames(p.Results.components);
            nFields = length(fieldNames);
            args = cell(2, nFields);
            args(1, :) = fieldNames;
            args(2, :) = {{cell(p.Results.nLevels, 1)}};
            obj.errors = struct(args{:});
            obj.count = 1;
        end

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
        
        function obj = addExact(obj, name, func)
            % ADDEXACT Add exact solution.
            %
            %   obj = addExact(obj, name, style) registers a new exact
            %   solution with the specified function handle.
            %
            % Inputs:
            %   obj - The Visualizer object
            %   name - Solution name (string or char array)
            %   func - Exact funtion (function handle)
            %
            % Outputs:
            %   obj - The Visualizer object          
            
            obj.exacts.(name) = func;
        end

        function obj = setLevel(obj, mesh)
            % SETLEVEL Set current refinement level information.
            %
            %   obj = setLevel(obj, mesh) configures the analyzer for the
            %   current refinement level using mesh properties. This method
            %   extracts the mesh resolution and creates descriptive
            %   strings for the convergence profiler.
            %
            % Inputs:
            %   obj - The Analyzer object
            %   mesh - Mesh object with resolution and measure properties
            %
            % Outputs:
            %   obj - The Analyzer object

            h = mesh.measure;
            hStr = strrep(num2str(mesh.resolution), '  ', 'x');
            obj.profiler.resolutions(obj.count) = h;
            obj.profiler.descriptions{obj.count} = hStr;
        end

        function obj = absolute(obj, space, dofs, varargin)
            % ABSOLUTE Compute absolute error against exact solutions.
            %
            %   obj = absolute(obj, space, dofs) computes absolute errors
            %   by comparing numerical solutions with exact analytical
            %   solutions. The method evaluates both numerical and exact
            %   solutions at the same sample points and computes specified
            %   error metrics.
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
                U0 = space.evaluate([], f, X, varargin{:});
                U0 = U0(:, component);
                obj.errors.(fieldName){obj.count} = metric.evaluate(U, U0);
            end
            obj.count = obj.count + 1;
        end

        function obj = relative(obj, fSpace, fDofs, cSpace, cDofs)
            % RELATIVE Compute relative error using Richardson
            % extrapolation.
            %
            %   obj = relative(obj, fSpace, fDofs, cSpace, cDofs) computes
            %   relative errors using Richardson extrapolation between fine
            %   and coarse grid solutions. This method is useful when exact
            %   solutions are not available.
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
                C = fDofs.(fieldName)(metric.iNestDofs, :);
                C0 = cDofs.(fieldName);
                U = fSpace.evaluate([], X, C);
                U = U(:, component);
                U0 = cSpace.evaluate([], X0, C0);
                U0 = U0(:, component);
                obj.errors.(fieldName){obj.count} = metric.evaluate(U, U0);
            end
            obj.count = obj.count + 1;
        end
    
        function results = analyze(obj)
            % ANALYZE Generate comprehensive error analysis results.
            %
            %   results = analyze(obj) processes all collected error data
            %   and generates structured results including convergence
            %   tables for each field component and error reduction type.
            %   The results can be used for convergence verification and
            %   method validation.
            %
            % Inputs:
            %   obj - The Analyzer object
            %
            % Outputs:
            %   results - Convergence tables

            results = struct();
            nLevels = obj.profiler.nLevels;
            nReductions = length(obj.reductions);
            fieldNames = fieldnames(obj.errors);
            for iField = 1:length(fieldNames)
                fieldName = fieldNames{iField};
                component = obj.components.(fieldName);
                E = cell2mat(obj.errors.(fieldName));
                for iComponent = 1:size(E, 2)
                    Ej = reshape(E(:, iComponent), nReductions, nLevels).';
                    name = sprintf('%s%d', fieldName, component(iComponent));
                    results.(name) = obj.profiler.buildTable(Ej, obj.reductions);
                end
            end
        end
    end
end