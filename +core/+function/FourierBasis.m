classdef FourierBasis < core.function.Function
    % FOURIERBASIS Fourier orthogonal basis functions.
    %
    %   FourierBasis generates trigonometric Fourier basis functions up to
    %   a specified maximum degree. The basis functions include constant,
    %   sine, and cosine terms, orthogonal on the interval \f$[0, 2\pi]\f$.
    %
    % See also:
    %   core.function.Function,
    %   core.function.SphericalHarmonicBasis
    
    properties
        MaxDegree % Maximum degree
        IsNormalized % Flag for orthonormal basis functions
        Lower % Lower bound
        Upper % Upper bound
    end
    
    methods
        function obj = FourierBasis(maxDegree, options)
            % FOURIERBASIS Construct Fourier basis functions.
            %
            %   obj = FourierBasis(maxDegree) creates Fourier basis
            %   functions  up to @a maxDegree.
            %
            %   obj = FourierBasis(maxDegree, isNormalized=true) creates
            %   orthonormal Fourier basis functions.
            
            arguments
                maxDegree(1, 1) {mustBeNonnegative, mustBeInteger}
                options.isNormalized(1, 1) {mustBeNumericOrLogical} = false
                options.lower(1, 1) {mustBeNumeric, mustBeFinite} = 0
                options.upper(1, 1) {mustBeNumeric, mustBeFinite} = 2 * pi
            end
            
            obj = obj@core.function.Function(nDims = 1, nCodims = 2*maxDegree+1);
            obj.MaxDegree = maxDegree;
            obj.IsNormalized = options.isNormalized;
            obj.Lower = options.lower;
            obj.Upper = options.upper;
        end
    end
    
    methods (Access = protected)
        function Y = evalImpl(obj, X)
            % EVALIMPL Implementation of function evaluation.
            
            L = obj.Upper - obj.Lower;
            Y = zeros(obj.NCodims, size(X, 2));
            if obj.IsNormalized
                Y(1, :) = 1 / sqrt(L);
            else
                Y(1, :) = 1;
            end
            for k = 1:obj.MaxDegree
                if obj.IsNormalized
                    normalizer = sqrt(2 / L);
                else
                    normalizer = 1;
                end
                Y(2*k, :) = normalizer * cos(k*X);
                Y(2*k+1, :) = normalizer * sin(k*X);
            end
        end
    end
end