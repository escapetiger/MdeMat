classdef SphericalHarmonicBasis < core.function.Function
    % SPHERICALHARMONICBASIS Spherical harmonic orthogonal basis functions.
    %
    %   SphericalHarmonicBasis generates spherical harmonic basis functions
    %   up to a specified maximum degree. The basis functions are
    %   orthogonal on the unit sphere with respect to the standard
    %   spherical measure.
    %
    % See also:
    %   core.function.Function, core.function.FourierBasis

    properties
        MaxDegree % Maximum degree l
        IsNormalized % Flag for orthonormal basis functions
    end

    methods
        function obj = SphericalHarmonicBasis(maxDegree, options)
            % SPHERICALHARMONICBASIS Construct spherical harmonic basis
            % functions.
            %
            %   obj = SphericalHarmonicBasis(maxDegree) creates
            %   spherical harmonic basis functions up to @a maxDegree.
            %
            %   obj = SphericalHarmonicBasis(maxDegree, isNormalized=true)
            %   creates orthonormal spherical harmonic basis functions.

            arguments
                maxDegree(1, 1) {mustBeNonnegative, mustBeInteger}
                options.isNormalized(1, 1) {mustBeNumericOrLogical} = false
            end

            obj = obj@core.function.Function(nDims = 2, nCodims = (maxDegree + 1)^2);
            obj.MaxDegree = maxDegree;
            obj.IsNormalized = options.isNormalized;
        end
    end

    methods (Access = protected)
        function Y = evalImpl(obj, X)
            % EVALIMPL Implementation of function evaluation.

            np = size(X, 2);
            Y  = zeros(obj.NCodims, np);
        
            mu  = X(1, :);
            phi = X(2, :);
        
            idx = 1;
            for l = 0:obj.MaxDegree
                %< MATLAB legendre returns P_l^m(mu), m=0..l, with Condon–Shortley phase included.
                P = legendre(l, mu);                  %< size: (l+1) x np, rows m=0..l
                P = P([1, repelem(2:l+1, 2)], :);     %< size: (2*l+1) x np, order [0,1,1,2,2,...]
        
                %< Trig part for real basis: [1, cos(phi), sin(phi), cos(2phi), sin(2phi), ...]
                T = zeros(2*l+1, np);
                T(1, :) = 1;
                if l >= 1
                    mvec = (1:l).';
                    T(2:2:end, :) = cos(mvec * phi);
                    T(3:2:end, :) = sin(mvec * phi);
                end
        
                Z = P .* T;
        
                if obj.IsNormalized
                    %< Real orthonormal spherical harmonics:
                    %<   m=0:        N_{l0} P_l^0
                    %<   m=+k:  sqrt(2) N_{lk} P_l^k cos(k phi)
                    %<   m=-k:  sqrt(2) N_{lk} P_l^k sin(k phi), k>=1
                    %<
                    %< N_{lm} = sqrt((2l+1)/(4pi) * (l-m)!/(l+m)!)
                    %<
                    %< Use gammaln for stability (avoids factorial overflow).
                    m_real = [0, reshape([1:l; -(1:l)], 1, [])].';   % (2*l+1)x1
                    m_abs  = abs(m_real);
        
                    log_ratio = gammaln(l - m_abs + 1) - gammaln(l + m_abs + 1);
                    Nlm = sqrt((2*l + 1) / (4*pi) .* exp(log_ratio)); % (2*l+1)x1
        
                    S = ones(2*l+1, 1);  %< sqrt(2) for nonzero m in real basis
                    if l >= 1
                        S(2:end) = sqrt(2);
                    end
        
                    %< IMPORTANT:
                    %< Do NOT multiply by (-1).^m here: MATLAB's legendre already includes
                    %< the Condon–Shortley phase. Adding it again would flip signs for odd m.
                    Z = (-1).^m_abs .* Z .* (Nlm .* S);
                end
        
                Y(idx:idx+2*l, :) = Z;
                idx = idx + 2*l + 1;
            end
        end
    end
end