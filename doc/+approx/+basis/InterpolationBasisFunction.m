classdef InterpolationBasisFunction < core.function.CompiledFunction
    % INTERPOLATIONBASISFUNCTION Compiled interpolation basis functions.
    %
    %   InterpolationBasisFunction represents a set of interpolation basis
    %   functions (such as Lagrange or Hermite polynomials) with automatic
    %   metadata loading from the core.symbolic compilation system.
    %
    %   The class automatically discovers appropriate metadata files based
    %   on the specified basis type, interval type, node type, and number
    %   of basis functions, providing efficient evaluation of both function
    %   values and derivatives through precompiled function handles.
    %
    %   Supported configurations:
    %   - Basis types: Lagrange, Hermite polynomials
    %   - Intervals: unit [0,1], canonical [-1,1]
    %   - Node types: Gauss-Legendre, Gauss-Lobatto quadrature points
    %
    % Examples:
    %   % Lagrange basis on unit interval with Gauss-Legendre nodes
    %   basis = InterpolationBasisFunction(5, 'lagrange', 'unit', 'gauss_legendre');
    %   basis = basis.autoLoad();
    %   y = basis.evaluate([0.25]);
    %   dy = basis.derivative([0.25], [1]);
    %   
    %   % Hermite basis on canonical interval with Gauss-Lobatto nodes
    %   basis = InterpolationBasisFunction(3, 'hermite', 'canonical', 'gauss_lobatto');
    %   basis = basis.autoLoad();
    %   
    %   % Manual metadata loading
    %   basis = InterpolationBasisFunction(4, 'lagrange', 'unit', 'gauss_legendre');
    %   basis = basis.load('lagrange_unit_gauss_legendre_4');
    %
    % See Also:
    %   core.function.CompiledFunction,
    %   core.function.OrthogonalBasisFunction

    properties
        basisType     % Type of interpolation basis ('lagrange', 'hermite')
        intervalType  % Type of interval ('unit', 'canonical')
        nodeType      % Type of node distribution ('gauss_legendre', 'gauss_lobatto')
    end

    methods
        function obj = InterpolationBasisFunction(nCodims, basisType, intervalType, nodeType)
            % INTERPOLATIONBASISFUNCTION Constructor for InterpolationBasisFunction.
            %
            %   obj = InterpolationBasisFunction(nCodims, basisType,
            %   intervalType, nodeType) creates a new
            %   InterpolationBasisFunction object with the specified
            %   configuration. Validates all input parameters and
            %   initializes the parent CompiledFunction.
            %
            % Inputs:
            %   nCodims - Number of basis functions (positive integer)
            %   basisType - String specifying basis type ('lagrange', 'hermite')
            %   intervalType - String specifying interval type ('unit', 'canonical')
            %   nodeType - String specifying node type ('gauss_legendre', 'gauss_lobatto')
            %
            % Outputs:
            %   obj - The constructed InterpolationBasisFunction object

            core.except.assert(nCodims > 0 && isscalar(nCodims) && ...
                nCodims == floor(nCodims), 'InvalidInput', ...
                'Number of basis functions must be a positive integer.');

            supportedBasisTypes = {'lagrange', 'hermite'};
            core.except.assert( ...
                ismember(lower(basisType), supportedBasisTypes), ...
                'InvalidInput', ...
                'Supported basis types: %s', strjoin(supportedBasisTypes, ', '));

            supportedIntervalTypes = {'unit', 'canonical'};
            core.except.assert( ...
                ismember(lower(intervalType), supportedIntervalTypes), ...
                'InvalidInput', ...
                'Supported interval types: %s', strjoin(supportedIntervalTypes, ', '));

            supportedNodeTypes = {'gauss_legendre', 'gauss_lobatto'};
            core.except.assert( ...
                ismember(lower(nodeType), supportedNodeTypes), ...
                'InvalidInput', ...
                'Supported node types: %s', strjoin(supportedNodeTypes, ', '));

            obj@core.function.CompiledFunction(1, nCodims);
            obj.basisType = lower(basisType);
            obj.intervalType = lower(intervalType);
            obj.nodeType = lower(nodeType);
        end
    end

    methods (Access = protected)
        function filename = autoFilename(obj)
            % AUTOFILENAME Generate appropriate metadata filename for this basis.
            %
            %   filename = autoFilename(obj) generates the metadata filename
            %   based on the basis configuration. The filename follows the
            %   pattern: {basisType}_{intervalType}_{nodeType}_{nCodims}
            %
            % Inputs:
            %   obj - The InterpolationBasisFunction object
            %
            % Outputs:
            %   filename - String of metadata filename

            filename = sprintf('%s_%s_%s_%d', obj.basisType, obj.intervalType, ...
                obj.nodeType, obj.nCodims);

            core.except.assert( ...
                ismember(filename, obj.metasource), ...
                'FileNotFound', ...
                'Metadata file not found: %s. Available files: %s', ...
                filename, strjoin(obj.metasource, ', '));
        end
    end
end