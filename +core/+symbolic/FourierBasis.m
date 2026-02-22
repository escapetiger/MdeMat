classdef FourierBasis < core.symbolic.SymbolicFunction
    % FOURIERBASIS Fourier orthogonal basis functions.
    %
    %   FourierBasis generates trigonometric Fourier basis functions up to
    %   a specified maximum degree. The basis functions include constant,
    %   sine, and cosine terms, orthogonal on the interval \f$[0, 2\pi]\f$.
    %
    % See also:
    %   core.symbolic.SymbolicFunction,
    %   core.symbolic.SphericalHarmonicBasis

    properties
        MaxDegree % Maximum degree
        IsNormalized % Flag for orthonormal basis functions
        Lower % Lower bound
        Upper % Upper bound
    end

    methods
        function obj = FourierBasis(variables, maxDegree, options)
            % FOURIERBASIS Construct Fourier basis functions.
            %
            %   obj = FourierBasis(variables, maxDegree) creates Fourier
            %   basis functions using the specified @a variables up to @a
            %   maxDegree.
            %
            %   obj = FourierBasis(variables, maxDegree, isNormalized=true)
            %   creates orthonormal Fourier basis functions.

            arguments
                variables {mustBeA(variables, 'sym')}
                maxDegree {mustBeNonnegative, mustBeInteger}
                options.isNormalized {mustBeNumericOrLogical} = false
                options.lower {mustBeNumeric, mustBeFinite} = 0
                options.upper {mustBeNumeric, mustBeFinite} = 2 * pi
            end

            obj = obj@core.symbolic.SymbolicFunction(variables);
            obj.MaxDegree = maxDegree;
            obj.IsNormalized = options.isNormalized;
            obj.Lower = options.lower;
            obj.Upper = options.upper;
            obj.computeExpressions();
        end
    end

    methods (Access = private)
        function computeExpressions(obj)
            % COMPUTEEXPRESSIONS Compute Fourier basis function
            % expressions.

            obj.Expressions = repmat(sym([]), 2*obj.MaxDegree+1, 1);
            L = obj.Upper - obj.Lower;
            x = obj.Variables(1);
            x = 2 * sym(pi) * (x - obj.Lower) / L;

            if obj.IsNormalized
                obj.Expressions(1) = sqrt(2*sym(pi)) / sqrt(2*sym(pi)*L);
            else
                obj.Expressions(1) = sym(1);
            end

            for k = 1:obj.MaxDegree

                if obj.IsNormalized
                    normalizer = sqrt(2*sym(pi)) / sqrt(sym(pi)*L);
                else
                    normalizer = sym(1);
                end

                obj.Expressions(2*k) = normalizer * cos(k*x);
                obj.Expressions(2*k+1) = normalizer * sin(k*x);
            end
        end
    end
end