classdef OrthogonalBasisFunction < core.function.CompiledFunction
    % ORTHOGONALBASISFUNCTION Compiled orthogonal polynomial basis
    % functions.
    %
    %   OrthogonalBasisFunction represents a set of orthogonal polynomial
    %   basis functions (such as Legendre or Chebyshev polynomials) with
    %   automatic metadata loading from the core.symbolic compilation
    %   system.
    %
    %   The class automatically discovers appropriate metadata files based
    %   on the specified basis type, interval type, and number of basis
    %   functions, providing efficient evaluation of both function values
    %   and derivatives through precompiled function handles.
    %
    %   Supported orthogonal polynomial families:
    %   - Legendre polynomials: orthogonal on \f$[-1,1]\f$ with weight
    %   \f$w(x) = 1\f$
    %   - Chebyshev polynomials: orthogonal on \f$[-1,1]\f$ with weight
    %   \f$w(x) = \frac{1}{\sqrt{1-x^2}}\f$
    %   - Monic Legendre: Legendre polynomials with leading coefficient = 1
    %
    % Examples:
    %   % Standard Legendre basis on unit interval
    %   basis = OrthogonalBasisFunction(5, 'legendre', 'unit');
    %   basis = basis.autoLoad();
    %   y = basis.evaluate([0.25]);
    %   dy = basis.derivative([0.25], [1]);
    %   
    %   % Monic Legendre polynomials on canonical interval
    %   basis = OrthogonalBasisFunction(4, 'monic_legendre', 'canonical');
    %   basis = basis.autoLoad();
    %   
    %   % Chebyshev polynomials with manual loading
    %   basis = OrthogonalBasisFunction(6, 'chebyshev', 'unit');
    %   basis = basis.load('chebyshev_unit_6');
    %
    % See Also:
    %   core.function.CompiledFunction,
    %   core.function.InterpolationBasisFunction

    properties
        basisType    % Type of orthogonal basis ('monic_legendre', 'legendre', 'chebyshev')
        intervalType % Type of interval ('unit', 'canonical')
    end

    methods
        function obj = OrthogonalBasisFunction(nCodims, basisType, intervalType)
            % ORTHOGONALBASISFUNCTION Constructor for
            % OrthogonalBasisFunction.
            %
            %   obj = OrthogonalBasisFunction(nCodims, basisType,
            %   intervalType) creates a new OrthogonalBasisFunction object
            %   with the specified configuration. Validates all input
            %   parameters and initializes the parent
            %   CompiledFunction.
            %
            % Inputs:
            %   nCodims - Number of basis functions (positive integer)
            %   basisType - Basis type ('monic_legendre', 'legendre', 'chebyshev', 'fourier', 'spherical_harmonic')
            %   intervalType - Interval type ('unit', 'canonical')
            %
            % Outputs:
            %   obj - The constructed OrthogonalBasisFunction object

            basisType = lower(basisType);
            switch basisType
                case {'monic_legendre', 'legendre', 'chebyshev', 'fourier'}
                    nDims = 1;
                case 'spherical_harmonic'
                    nDims = 2;
            end
            intervalType = lower(intervalType);

            obj@core.function.CompiledFunction(nDims, nCodims);
            obj.basisType = basisType;
            obj.intervalType = intervalType;
        end
    end

    methods (Access = protected)
        function filename = autoFilename(obj)
            % AUTOFILENAME Generate appropriate metadata filename for this
            % basis.
            %
            %   filename = autoFilename(obj) generates the metadata
            %   filename based on the basis configuration. The filename
            %   follows the pattern: {basisType}_{intervalType}_{nCodims}
            %
            % Inputs:
            %   obj - The OrthogonalBasisFunction object
            %
            % Outputs:
            %   filename - String of metadata filename

            switch obj.basisType
                case {'monic_legendre', 'legendre', 'chebyshev', 'fourier'}
                    filename = sprintf('%s_%s_%d', obj.basisType, obj.intervalType, obj.nCodims);
                case 'spherical_harmonic'
                    filename = sprintf('%s_%d', obj.basisType, obj.nCodims);
            end

            core.except.assert( ...
                ismember(filename, obj.metasource), ...
                'FileNotFound', ...
                'Metadata file not found: %s. Available files: %s', ...
                filename, strjoin(obj.metasource, ', '));
        end
    end
end