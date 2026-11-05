classdef PointSet < handle
    % POINTSET Point distributions within geometric domains.
    %
    %   PointSet generates and manages point distributions within various
    %   geometric domains using appropriate integration rules. It supports
    %   orthotope (rectangular) and hypersphere geometries with multiple
    %   point generation strategies including deterministic quadrature
    %   and Monte Carlo methods.
    %
    % Examples:
    %   % Create points in 2D rectangle using Gauss-Legendre rule
    %   orthotope = msclab.geometry.Orthotope([-1, 0], [1, 2]);
    %   pointSet = PointSet(orthotope, [5, 5], 'RuleType', 'gauss_legendre');
    %   
    %   % Create points on 3D sphere using Lebedev rule
    %   sphere = msclab.geometry.Hypersphere([0, 0, 0], 1);
    %   pointSet = PointSet(sphere, 110, 'RuleType', 'lebedev');
    %   
    %   % Monte Carlo points with custom seed
    %   pointSet = PointSet(sphere, 1000, 'RuleType', 'monte_carlo', 'Seed', 42);
    %
    % Notes:
    %   The class automatically selects appropriate integrators based on
    %   geometry type and provides a unified interface for accessing
    %   nodes, weights, and point set properties.
    %
    % See also:
    %   approx.integrate.OrthotopeIntegrator, approx.integrate.HypersphereIntegrator,
    %   msclab.geometry.Orthotope, msclab.geometry.Hypersphere
    
    properties (Access = public)
        integrator % Integration object for point generation
    end

    properties (Dependent)
        nodes   % Matrix of point coordinates (nDims × nPoints)
        weights % Vector of integration weights (1 × nPoints)
        nPoints % Number of points in the set
        nDims   % Number of spatial dimensions
    end
    
    methods
        function obj = PointSet(G, n, varargin)
            % POINTSET Constructor for PointSet.
            %
            %   obj = PointSet(G, n) creates a point set in geometry G
            %   with default resolution n and automatic rule selection.
            %
            %   obj = PointSet(G, n, Name, Value) creates a point set
            %   with additional options specified as name-value pairs.
            %
            % Inputs:
            %   G - Geometry object (Orthotope or Hypersphere)
            %   n - Point resolution (scalar or vector)
            %   varargin - Name-value pairs:
            %<     'RuleType' - Integration rule type string
            %<                  For orthotopes: 'gauss_legendre', 'gauss_lobatto'
            %<                  For spheres: 'lebedev', 'gauss_trapezoidal', 'monte_carlo'
            %<     'Seed' - Random seed for Monte Carlo (positive integer)
            %<     'CoordinateSystem' - Coordinate system for spheres (1=Cartesian, 2=Spherical)
            %<     'Normalization' - Weight normalization flag (logical)
            %
            % Outputs:
            %   obj - Constructed PointSet object
            
            parser = inputParser;
            addParameter(parser, 'RuleType', '', @(x) ischar(x) || isstring(x));
            addParameter(parser, 'Seed', [], @(x) isnumeric(x) && isscalar(x));
            addParameter(parser, 'CoordinateSystem', 1, @(x) isnumeric(x) && isscalar(x));
            addParameter(parser, 'Normalization', false, @(x) islogical(x) || isnumeric(x));
            parse(parser, varargin{:});
            
            ruleType = parser.Results.RuleType;
            seed = parser.Results.Seed;
            coordSys = parser.Results.CoordinateSystem;
            normalization = parser.Results.Normalization;
            
            if isa(G, 'msclab.geometry.Orthotope')
                obj.createOrthotopeIntegrator(G, n, ruleType, normalization);
            elseif isa(G, 'msclab.geometry.Hypersphere')
                obj.createHypersphereIntegrator(G, n, ruleType, seed, coordSys, normalization);
            else
                core.except.assert(0, 'InvalidInput', ...
                    'Geometry type %s is not supported.', class(G));
            end
        end
    
        function X = get.nodes(obj)
            % GET.NODES Get the point coordinates.
            
            X = obj.integrator.nodes;
        end
        
        function w = get.weights(obj)
            % GET.WEIGHTS Get the integration weights.
            
            w = obj.integrator.weights;
        end
        
        function n = get.nPoints(obj)
            % GET.NPOINTS Get the number of points.
            
            n = obj.integrator.nPoints;
        end
        
        function n = get.nDims(obj)
            % GET.NDIMS Get the number of dimensions.

            n = obj.integrator.nDims;
        end
    end
    
    methods (Access = private)
        function createOrthotopeIntegrator(obj, G, n, ruleType, normalization)
            % CREATEORTHOTOPEINTEGRATOR Create orthotope integrator.
            %
            %   createOrthotopeIntegrator(obj, G, n, ruleType,
            %   normalization) creates an integrator for generating points
            %   within an orthotope (rectangular domain).
            %
            % Inputs:
            %   obj - The PointSet object
            %   G - Orthotope geometry object
            %   n - Point resolution (scalar or vector)
            %   ruleType - Integration rule type string
            %   normalization - Weight normalization flag
            % 
            % Outputs:
            %   NULL
            
            if isempty(ruleType)
                ruleType = 'gauss_legendre';
            end
            
            I = msclab.integrate.OrthotopeIntegrator(G.nDims, ruleType);
            I.setPoints(n, normalization, G.lower, G.upper);
            obj.integrator = I;
        end
        
        function createHypersphereIntegrator(obj, G, n, ruleType, seed, coordSys, normalization)
            % CREATEHYPERSPHEREINTEGRATOR Create hypersphere integrator.
            %
            %   createHypersphereIntegrator(obj, G, n, ruleType, seed,
            %   coordSys, normalization) creates an integrator for
            %   generating points on or within a hypersphere.
            %
            % Inputs:
            %   obj - The PointSet object
            %   G - Hypersphere geometry object
            %   n - Point resolution (scalar or vector)
            %   ruleType - Integration rule type string
            %   seed - Random seed for Monte Carlo methods
            %   coordSys - Coordinate system specification
            %   normalization - Weight normalization flag
            %
            % Outputs:
            %   NULL
            
            if isempty(ruleType)
                if G.nDims == 3
                    ruleType = 'lebedev';
                else
                    ruleType = 'gauss_trapezoidal';
                end
            end
            
            if strcmpi(ruleType, 'monte_carlo') && ~isempty(seed)
                I = msclab.integrate.HypersphereIntegrator(G.nDims, ruleType, seed);
            else
                I = msclab.integrate.HypersphereIntegrator(G.nDims, ruleType);
            end
            
            I.setPoints(n, normalization, G.centroid, G.radius, coordSys);
            obj.integrator = I;
        end
    end
end