classdef InterpolationPolynomialBasisFunction < core.symbolic.PolynomialBasisFunction
    % INTERPOLATIONPOLYNOMIALBASISFUNCTION Base class for all interpolation polynomial basis functions.
    %
    %   InterpolationPolynomialBasisFunction provides common functionality
    %   for interpolation-based polynomial basis functions such as Lagrange
    %   and Hermite polynomials. It handles node management and
    %   automatically calculates the polynomial degree based on the number
    %   of nodes.
    %
    % Examples:
    %   % Cannot instantiate directly - use concrete subclasses
    %   % nodes = [0, 1, 2];
    %   % basis = LagrangeBasisFunction(nodes, 1);
    %
    % See also:
    %   core.symbolic.PolynomialBasisFunction
    
    properties (Access = public)
        nodes % Array of interpolation nodes
    end
    
    methods
        function obj = InterpolationPolynomialBasisFunction(nodes, index, variable, monic)
            % INTERPOLATIONPOLYNOMIALBASISFUNCTION Constructor for
            % InterpolationPolynomialBasisFunction.
            %
            %   obj = InterpolationPolynomialBasisFunction(nodes, index)
            %   creates a new interpolation polynomial basis function with
            %   the specified nodes and index using default settings.
            %   
            %   obj = InterpolationPolynomialBasisFunction(nodes, index,
            %   variable, monic) creates a new interpolation polynomial
            %   basis function with full customization options.
            %
            % Inputs:
            %   nodes - Array of interpolation nodes (vector)
            %   index - Index of the basis function (positive integer)
            %   variable - Symbol(s) (optional, default: {sym('x')})
            %   monic - Monic form flag (optional, default: false)
            %
            % Outputs:
            %   obj - Constructed InterpolationPolynomialBasisFunction object
            
            core.except.assert(isvector(nodes) && ~isempty(nodes), ...
                'InvalidInput', 'Nodes must be a non-empty vector.');
            
            if nargin < 3
                variable = sym('x');
            end
            if nargin < 4
                monic = false;
            end
            
            degree = length(nodes) - 1;
            bbox = [min(nodes), max(nodes)];
            obj@core.symbolic.PolynomialBasisFunction(degree, index, bbox, variable, monic);
            obj.nodes = nodes(:).';
            obj.expression = obj.generate();
        end
    end
end