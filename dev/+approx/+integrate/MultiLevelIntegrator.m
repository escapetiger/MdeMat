classdef MultiLevelIntegrator < handle
    % MULTILEVELINTEGRATOR Multi-level integrator with Richardson
    % extrapolation.
    %
    %   MultiLevelIntegrator implements multi-level numerical integration
    %   using an array of Integrator objects at different refinement
    %   levels. It supports Richardson extrapolation to achieve higher
    %   accuracy by combining results from multiple levels with known
    %   convergence rates.
    %
    %   The class manages a collection of Integrator objects, each
    %   operating at different resolutions, and provides methods for
    %   setting up multi-level integration points, evaluating functions at
    %   all levels, and computing extrapolated results.
    %
    % See also:
    %   approx.integrate.Integrator, core.geometry.Orthotope
    
    properties
        Integrators % Array of Integrator objects, one per level
    end
    
    properties (Dependent)
        NLevels % Number of refinement levels
        Nodes % Integration nodes
        Weights % Integration weights
    end
    
    methods
        function n = get.NLevels(obj)
            % GET.NLEVELS Returns the number of refinement levels.
            
            n = length(obj.Integrators);
        end

        function X = get.Nodes(obj)
            n = obj.NLevels;
            X = cell(1, n);
            for i = 1:n
                X{i} = obj.Integrators(i).Nodes;
            end
        end

        function X = get.Weights(obj)
            n = obj.NLevels;
            X = cell(1, n);
            for i = 1:n
                X{i} = obj.Integrators(i).Weights;
            end
        end
    end
    
    methods
        function obj = MultiLevelIntegrator(rule, options)
            % MULTILEVELINTEGRATOR Construct a MultiLevelIntegrator.
            %
            %   obj = MultiLevelIntegrator(rule) creates an multi-level
            %   integrator.
            %
            %   obj = MultiLevelIntegrator(rule, nLevels=nLevels) creates a
            %   multi-level integrator with the specified number of levels.
            
            arguments
                rule approx.integrate.IntegrationRule
                options.nLevels(1, 1) {mustBePositive, mustBeInteger} = 2
            end
            
            obj.Integrators = arrayfun(@(i) approx.integrate.Integrator(rule), 1:options.nLevels);
        end
        
        function newObj = copy(obj)
            % COPY Create a deep copy of the multi-level integrator.
            %
            %   newObj = copy(obj) creates a new MultiLevelIntegrator
            %   object with the same geometry, number of levels, and
            %   integrator configurations as the original object. All
            %   integrator data including nodes, weights, and values are
            %   copied.
            
            arguments
                obj approx.integrate.MultiLevelIntegrator
            end
            
            newObj = approx.integrate.MultiLevelIntegrator(nLevels = obj.NLevels);
            
            for i = 1:obj.NLevels
                newObj.Integrators(i) = obj.Integrators(i).copy();
            end
        end
        
        function obj = setPoints(obj, geometry, np)
            % SETPOINTS Set up multi-level integration points.
            %
            %   obj = setPoints(obj, geometry, np) generates
            %   multi-level integration points using the specified number
            %   of points @a np and quadrature type @a tp. Currently
            %   supports Orthotope geometries.
            
            arguments
                obj approx.integrate.MultiLevelIntegrator
                geometry core.geometry.Geometry
                np {mustBeVector, mustBePositive, mustBeInteger}
            end
            
            obj.setPointsOnOrthotope(geometry, np);
        end
        
        function obj = setValues(obj, func)
            % SETVALUES Evaluate function at all integration levels.
            %
            %   obj = setValues(obj, func) evaluates the specified function
            %   at the integration nodes for all levels and stores the
            %   results in each integrator's Values property.
            
            arguments
                obj approx.integrate.MultiLevelIntegrator
                func {mustBeFunction}
            end
            
            for i = 1:obj.NLevels
                obj.Integrators(i).setValues(func);
            end
        end
        
        function obj = addWeightFunction(obj, weightFunc)
            % ADDWEIGHTFUNCTION Apply weight function to all levels.
            %
            %   obj = addWeightFunction(obj, weightFunc) multiplies the
            %   integration weights by the specified weight function
            %   evaluated at the integration nodes for all levels.
            
            arguments
                obj approx.integrate.MultiLevelIntegrator
                weightFunc {mustBeFunction}
            end
            
            for i = 1:obj.NLevels
                obj.Integrators(i).addWeightFunction(weightFunc);
            end
        end
        
        function I = integrate(obj, func, options)
            % INTEGRATE Compute integral using multi-level extrapolation.
            %
            %   I = integrate(obj) integrates using stored function values
            %   from all levels and applies Richardson extrapolation.
            %
            %   I = integrate(obj, func) evaluates the specified function
            %   at all integration levels and computes the extrapolated
            %   integral.
            %
            %   I = integrate(obj, func, levelIndex=levelIndex) computes
            %   the integral using only the specified level.
            
            arguments
                obj approx.integrate.MultiLevelIntegrator
                func {mustBeFunction}
                options.levelIndex(1, 1) {mustBePositive, mustBeInteger} = []
                options.useExtrapolation(1, 1) {mustBeNumericOrLogical} = false
            end
            
            if ~isempty(options.levelIndex)
                I = obj.Integrators(options.levelIndex).integrate(func);
                return;
            end
            
            results = zeros(1, obj.NLevels);
            for i = 1:obj.NLevels
                results(i) = obj.Integrators(i).integrate(func);
            end
            
            if obj.NLevels == 1 || ~options.useExtrapolation
                I = results(end);
            else
                n = length(results);
                R = zeros(n, n);
                R(:, 1) = results(:);
                
                for j = 2:n
                    for i = j:n
                        q = 4^(j-1);
                        R(i, j) = (q * R(i, j-1) - R(i-1, j-1)) / (q - 1);
                    end
                end
                
                I = R(n, n);
            end
        end
    end
    
    methods (Access = private)
        function obj = setPointsOnOrthotope(obj, orthotope, np)
            % SETPOINTSONORTHOTOPE Set points for Orthotope geometry.

            for level = 1:obj.NLevels
                [X, w] = obj.setPointsOnOrthotopeImpl(level, orthotope, np);
                obj.Integrators(level).Nodes = X;
                obj.Integrators(level).Weights = w;
            end
        end
        
        function [X, w] = setPointsOnOrthotopeImpl(obj, level, orthotope, np)
            % SETPOINTSONORTHOTOPEIMPL Implementation of refined
            % integration points for Orthotope geometry.
            
            nd = orthotope.NDims;
            if level == 1
                rule = obj.Integrators(level).Rule;
                [X, w] = rule.generate(orthotope, np);
            else
                a = orthotope.Lower;
                b = orthotope.Upper;
                m = 2^nd;
                c = (a + b) / 2;
                X0 = cell(1, m);
                w0 = cell(1, m);
                for j = 1:m
                    t = dec2bin(j - 1, nd) - '0';
                    p = a .* (1 - t) + c .* t;
                    q = c .* (1 - t) + b .* t;
                    bbox = [p(:).'; q(:).'];
                    g = core.geometry.Orthotope(bbox(:));
                    [X0{j}, w0{j}] = obj.setPointsOnOrthotopeImpl(level - 1, g, np);
                end
                X = [X0{:}];
                w = [w0{:}];
            end
        end
    end
end