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
    %   mdemat.common.quadrature.Integrator, mdemat.data.geometry.Orthotope
    
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
                rule mdemat.common.quadrature.IntegrationRule
                options.nLevels(1, 1) {mustBePositive, mustBeInteger} = 2
            end
            
            obj.Integrators = arrayfun(@(i) mdemat.common.quadrature.Integrator(rule), 1:options.nLevels);
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
                obj mdemat.common.quadrature.MultiLevelIntegrator
            end
            
            newObj = mdemat.common.quadrature.MultiLevelIntegrator( ...
                obj.Integrators(1).Rule, nLevels=obj.NLevels);
            
            for i = 1:obj.NLevels
                newObj.Integrators(i) = obj.Integrators(i).copy();
            end
        end
        
        function obj = setPoints(obj, geometry, np)
            % SETPOINTS Set up multi-level integration points.
            %
            %   obj = setPoints(obj, geometry, np) generates
            %   multi-level integration points using the specified number
            %   of points @a np. Currently
            %   supports Orthotope geometries.
            
            arguments
                obj mdemat.common.quadrature.MultiLevelIntegrator
                geometry mdemat.data.geometry.Geometry
                np {mustBeVector, mustBePositive, mustBeInteger}
            end
            
            assert(all(isa(geometry, 'mdemat.data.geometry.Orthotope'), 'all'), ...
                'mdemat:common:quadrature:MultiLevelIntegrator:InvalidGeometry', '%s', 'Multilevel refinement requires an orthotope.');
            obj.setPointsOnOrthotope(geometry, np);
        end
        
        function obj = setValues(obj, func)
            % SETVALUES Evaluate function at all integration levels.
            %
            %   obj = setValues(obj, func) evaluates the specified function
            %   at the integration nodes for all levels and stores the
            %   results in each integrator's Values property.
            
            arguments
                obj mdemat.common.quadrature.MultiLevelIntegrator
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
                obj mdemat.common.quadrature.MultiLevelIntegrator
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
            %   useExtrapolation=true assumes an even-power error expansion
            %   I(h)=I+c_2*h^2+c_4*h^4+... and successive mesh halving.
            
            arguments
                obj mdemat.common.quadrature.MultiLevelIntegrator
                func = []
                options.levelIndex {mustBePositive, mustBeInteger} = []
                options.useExtrapolation(1, 1) {mustBeNumericOrLogical} = false
            end
            
            if ~isempty(func)
                mustBeFunction(func);
            end
            if ~isempty(options.levelIndex)
                assert(all(isscalar(options.levelIndex) && ...
                    options.levelIndex <= obj.NLevels, 'all'), ...
                    'mdemat:common:quadrature:MultiLevelIntegrator:InvalidLevel', '%s', 'Requested integration level is unavailable.');
                I = obj.Integrators(options.levelIndex).integrate(func);
                return;
            end

            results = cell(obj.NLevels, obj.NLevels);
            for i = 1:obj.NLevels
                results{i, 1} = obj.Integrators(i).integrate(func);
            end
            if obj.NLevels == 1 || ~options.useExtrapolation
                I = results{obj.NLevels, 1};
                return;
            end

            % Richardson extrapolation acts componentwise. The refinement
            % factor is two and the leading error is assumed second order.
            for j = 2:obj.NLevels
                factor = 4^(j - 1);
                for i = j:obj.NLevels
                    results{i, j} = (factor * results{i, j - 1} - ...
                        results{i - 1, j - 1}) / (factor - 1);
                end
            end
            I = results{obj.NLevels, obj.NLevels};
        end
    end
    
    methods (Access = private)
        function obj = setPointsOnOrthotope(obj, orthotope, np)
            % SETPOINTSONORTHOTOPE Set points for Orthotope geometry.

            for level = 1:obj.NLevels
                [X, w] = obj.setPointsOnOrthotopeImpl(level, orthotope, np);
                obj.Integrators(level).Nodes = X;
                obj.Integrators(level).Weights = w;
                obj.Integrators(level).Values = [];
                obj.Integrators(level).Partials = [];
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
                    g = mdemat.data.geometry.Orthotope(bbox(:));
                    [X0{j}, w0{j}] = obj.setPointsOnOrthotopeImpl(level - 1, g, np);
                end
                X = [X0{:}];
                w = [w0{:}];
            end
        end
    end
end
