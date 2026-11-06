classdef SphericalHarmonicBasis < mdemat.common.symbolic.SymbolicFunction
    % SPHERICALHARMONICBASIS Spherical harmonic orthogonal basis functions.
    %
    %   SphericalHarmonicBasis generates spherical harmonic basis functions
    %   up to a specified maximum degree. The basis functions are
    %   orthogonal on the unit sphere with respect to the standard
    %   spherical measure.
    %   Modes are ordered by degree, then [m=0, cos(phi), sin(phi), ...].
    %   The real convention cancels the Condon--Shortley phase. For m>0,
    %   the orthonormal factor is sqrt(2)*N_lm, where
    %   N_lm = sqrt((2*l+1)/(4*pi)*(l-m)!/(l+m)!).
    %
    % See also:
    %   mdemat.common.symbolic.SymbolicFunction, mdemat.common.symbolic.FourierBasis

    properties
        MaxDegree % Maximum degree l
        IsNormalized % Flag for orthonormal basis functions
    end

    methods
        function obj = SphericalHarmonicBasis(variables, maxDegree, options)
            % SPHERICALHARMONICBASIS Construct spherical harmonic basis
            % functions.
            %
            %   obj = SphericalHarmonicBasis(variables, maxDegree) creates
            %   spherical harmonic basis functions using the specified
            %   @a variables (mu, phi) up to @a maxDegree.
            %
            %   obj = SphericalHarmonicBasis(variables, maxDegree,
            %   isNormalized=true) creates orthonormal spherical harmonic
            %   basis functions.

            arguments
                variables {mustBeA(variables, 'sym'), mustBeNonempty}
                maxDegree (1, 1) {mustBeNonnegative, mustBeInteger}
                options.isNormalized (1, 1) {mustBeNumericOrLogical} = false
            end

            assert(all(numel(variables) == 2, 'all'), 'mdemat:common:symbolic:SphericalHarmonicBasis:InvalidInput', ...
                '%s', 'Spherical harmonics require the variables [mu; phi].');
            obj = obj@mdemat.common.symbolic.SymbolicFunction(variables);
            obj.MaxDegree = maxDegree;
            obj.IsNormalized = options.isNormalized;
            obj.computeExpressions();
        end
    end

    methods (Access = private)
        function computeExpressions(obj)
            % COMPUTEEXPRESSIONS Compute spherical harmonic basis function
            % expressions.

            obj.Expressions = repmat(sym([]), (obj.MaxDegree + 1)^2, 1);

            mu = obj.Variables(1);
            phi = obj.Variables(2);

            idx = 1;

            for l = 0:obj.MaxDegree
                for m = 0:l
                    if m == 0
                        expr = obj.computeSphericalHarmonic(mu, phi, l, 0);
                        if obj.IsNormalized
                            normalizer = sqrt((2 * l + 1)/(4 * sym(pi))* ...
                                factorial(l-abs(m))/factorial(l+abs(m)));
                            expr = simplify(normalizer*expr);
                        end
                        obj.Expressions(idx) = expr;
                        idx = idx + 1;
                    else
                        expr1 = obj.computeSphericalHarmonic(mu, phi, l, m);
                        expr2 = obj.computeSphericalHarmonic(mu, phi, l, -m);
                        if obj.IsNormalized
                            normalizer = sqrt((2 * l + 1)/(2 * sym(pi))* ...
                                factorial(l-abs(m))/factorial(l+abs(m)));
                            expr1 = simplify(normalizer*expr1);
                            expr2 = simplify(normalizer*expr2);
                        end
                        obj.Expressions(idx) = (-1)^m * expr1;
                        obj.Expressions(idx+1) = (-1)^m * expr2;
                        idx = idx + 2;
                    end
                end
            end
        end

        function harmonic = computeSphericalHarmonic(obj, mu, phi, l, m)
            % COMPUTESPHERICALHARMONIC Compute spherical harmonic
            % Y_l^m(mu, phi).

            p = obj.computeAssociatedLegendre(mu, l, abs(m));

            if m > 0
                f = cos(m * phi);
            elseif m < 0
                f = sin(abs(m)*phi);
            else
                f = sym(1);
            end

            harmonic = simplify(p*f);
        end

        function p = computeAssociatedLegendre(obj, mu, l, m)
            % COMPUTEASSOCIATEDLEGENDRE Compute associated Legendre
            % polynomial P_l^m(mu).

            if m == 0
                p = obj.computeLegendre(mu, l);
            else
                p = obj.computeLegendre(mu, l);
                p = (-1)^m * (1 - mu^2)^(m / 2) * diff(p, mu, m);
                p = simplify(p);
                p = horner(p);
            end
        end

        function p = computeLegendre(obj, mu, n)
            % COMPUTELEGENDRE Compute Legendre polynomial P_n(mu).

            if n == 0
                p = sym(1);
            elseif n == 1
                p = mu;
            else
                p1 = obj.computeLegendre(mu, n-1);
                p2 = obj.computeLegendre(mu, n-2);
                p = ((2 * n - 1) * mu * p1 - (n - 1) * p2) / n;
                p = simplify(p);
                p = horner(p);
            end
        end
    end
end
