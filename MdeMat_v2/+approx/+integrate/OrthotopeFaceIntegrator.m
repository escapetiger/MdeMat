classdef OrthotopeFaceIntegrator < approx.integrate.Integrator
    % ORTHOTOPEFACEINTEGRATOR Numerical integration on orthotope faces.
    %
    %   OrthotopeFaceIntegrator provides specialized integration
    %   capabilities for functions defined on the faces (boundaries) of
    %   multidimensional orthotopes (hyperrectangles). Each face is
    %   identified by an index and corresponds to fixing one coordinate
    %   at either its minimum or maximum value.
    %
    % Examples:
    %   % Create integrator for face 1 of a 3D orthotope
    %   integrator = OrthotopeFaceIntegrator(3, 'gauss_legendre', 1);
    %   
    %   % Using existing rule object
    %   rule = approx.integrate.SeparableRule.gaussLobatto(2);
    %   integrator = OrthotopeFaceIntegrator(rule, 3);
    %
    % Notes:
    %   Face indexing: for d dimensions, faces 1,3,5,... correspond to
    %   minimum values and faces 2,4,6,... to maximum values of
    %   coordinates 1,2,3,... respectively.
    %
    % See also:
    %   approx.integrate.Integrator, approx.integrate.OrthotopeIntegrator,
    %   approx.integrate.SeparableRule

    properties (Access = public)
        faceIdx % Face index (positive integer)
    end

    methods
        function obj = OrthotopeFaceIntegrator(varargin)
            % ORTHOTOPEFACEINTEGRATOR Constructor for
            % OrthotopeFaceIntegrator.
            %
            %   obj = OrthotopeFaceIntegrator(rule, faceIdx) creates an
            %   integrator using the specified rule and face index.
            %
            %   obj = OrthotopeFaceIntegrator(nDims, ruleId, faceIdx)
            %   creates an integrator with the specified parameters.
            %
            % Inputs:
            %   varargin - Input arguments
            %<   rule - SeparableRule object for the face integration
            %<   faceIdx - Face index (positive integer)
            %<   nDims - Number of dimensions of the original orthotope
            %<   ruleId - Rule identifier string: {'gauss_legendre', 'gauss_lobatto'}
            %
            % Outputs:
            %   obj - Constructed OrthotopeFaceIntegrator object

            if nargin == 2
                [rule, faceIdx] = varargin{:};
            elseif nargin == 3
                [nDims, ruleId, faceIdx] = varargin{:};
                if strcmp(ruleId, 'gauss_legendre')
                    rule = approx.integrate.SeparableRule.gaussLegendre(nDims);
                elseif strcmp(ruleId, 'gauss_lobatto')
                    rule = approx.integrate.SeparableRule.gaussLobatto(nDims);
                else
                    error('Invalid rule identifier.');
                end
            else
                error('Invalid number of arguments.');
            end
            obj@approx.integrate.Integrator(rule);
            obj.faceIdx = faceIdx;
        end

        function newObj = copy(obj)
            % COPY Create a deep copy of the OrthotopeFaceintegrator.
            %
            %   newObj = copy(obj) creates a new OrthotopeFaceIntegrator
            %   object with the same rule, face index, nodes, and weights.
            %
            % Inputs:
            %   obj - The OrthotopeFaceIntegrator object to copy
            %
            % Outputs:
            %   newObj - Deep copy of the original integrator

            newObj = approx.integrate.OrthotopeFaceIntegrator(obj.rule, obj.faceIdx);
            newObj.nodes = obj.nodes;
            newObj.weights = obj.weights;
        end

        function obj = setPoints(obj, n, a, b)
            % SETPOINTS Generate integration points on the specified face.
            %
            %   obj = setPoints(obj, n) generates integration points on the
            %   face with @a n points per dimension for the unit orthotope.
            %
            %   obj = setPoints(obj, n, a, b) generates points for a
            %   custom orthotope with bounds [@a a, @a b].
            %
            % Inputs:
            %   obj - The OrthotopeFaceIntegrator object
            %   n - Number of points per dimension (vector or scalar)
            %   a - Lower bounds vector (optional, default: zeros)
            %   b - Upper bounds vector (optional, default: ones)
            %
            % Outputs:
            %   obj - The OrthotopeFaceIntegrator object

            d = obj.nDims;
            if nargin < 3
                a = zeros(1, d);
                b = ones(1, d);
            end
            i = obj.faceIdx;
            l = a(:).';
            u = b(:).';
            j = ceil(i/2);
            if mod(i, 2) == 0
                l(j) = u(j);
            else
                u(j) = l(j);
            end
            args = arrayfun(@(k) {n(k), l(k), u(k)}, 1:d, 'Un', 0);
            [obj.nodes, obj.weights] = obj.rule.generate(args{:});
        end
    end
end