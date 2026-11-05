classdef SphericalHarmonicBasisFunction < core.symbolic.BasisFunction
    % SPHERICALHARMONICBASISFUNCTION Real-valued spherical harmonic basis functions.
    %
    %   SphericalHarmonicBasisFunction represents real-valued spherical
    %   harmonic basis functions on the canonical domain \f$\mu\in[-1,1]\f$
    %   and \f$\theta\in[0,2\pi]\f$. The functions are constructed using
    %   associated Legendre polynomials and trigonometric functions.
    %
    %   The spherical harmonic function with degree \f$l\f$ and order
    %   \f$m\f$ are: 
    %   - \f$(l^2+1)\f$-th function (zonal harmonics, \f$m=0\f$): 
    %     \f$P_l^0(\mu)\f$
    %   - \f$(l^2+2m-1)\f$-th function (tesseral harmonics, \f$m>0\f$): 
    %     \f$\cos(m\theta)P_l^m(\mu)\f$ 
    %   - \f$(l^2 + 2m)\f$-th function (tesseral harmonics, \f$m>0\f$): 
    %     \f$\sin(m\theta)P_l^m(\mu)\f$
    %
    % Examples:
    %   % Zonal harmonic (constant)
    %   basis1 = core.symbolic.SphericalHarmonicBasisFunction(1);  % P_0^0(mu) = 1
    %
    %   % First degree harmonics
    %   basis2 = core.symbolic.SphericalHarmonicBasisFunction(2);  % P_1^0(mu) = mu
    %   basis3 = core.symbolic.SphericalHarmonicBasisFunction(3);  % cos(theta)*P_1^1(mu)
    %   basis4 = core.symbolic.SphericalHarmonicBasisFunction(4);  % sin(theta)*P_1^1(mu)
    %
    % See also:
    %   core.symbolic.BasisFunction
    
    methods
        function obj = SphericalHarmonicBasisFunction(index)
            % SPHERICALHARMONICBASISFUNCTION Constructor for
            % SphericalHarmonicBasisFunction.
            %
            %   obj = SphericalHarmonicBasisFunction(index) creates a new
            %   spherical harmonic basis function with fixed domain
            %   variables \f$\mu\in[-1,1]\f$ and \f$\theta\in[0,2\pi]\f$.
            %
            % Inputs:
            %   index - Index of the spherical harmonic (positive integer)
            %
            % Outputs:
            %   obj - Constructed SphericalHarmonicBasisFunction object
            
            bbox = [-1, 1, 0, 2*pi];
            variables = {sym('mu'), sym('theta')};
            obj@core.symbolic.BasisFunction(index, bbox, variables);
            obj.expression = obj.generate();
        end
    end
    
    methods (Access = protected)
        function expr = generate(obj)
            % GENERATE Generates the spherical harmonic basis function.
            %
            %   expr = generate(obj) creates the symbolic expression for
            %   the real-valued spherical harmonic based on the index.
            %   Determines the appropriate degree l and order m, then
            %   constructs the expression using associated Legendre
            %   polynomials and trigonometric functions.
            %
            % Inputs:
            %   obj - The SphericalHarmonicBasisFunction object
            %
            % Outputs:
            %   expr - Symbolic expression for the spherical harmonic
            
            mu = obj.variables{1};
            theta = obj.variables{2};
            
            %< Determine degree l and order m from index
            [l, m] = obj.getDegreesAndOrder();
            
            %< Generate associated Legendre polynomial P_l^m(mu)
            P_lm = obj.associatedLegendre(l, m, mu);
            
            %< Generate the appropriate spherical harmonic
            if m == 0
                %< Zonal harmonic: P_l^0(mu)
                expr = P_lm;
            else
                %< Determine cosine vs sine for tesseral harmonics
                remainder = obj.index - l^2 - 1;
                isCos = mod(remainder, 2) == 1;
                
                if isCos
                    %< Tesseral harmonic (cosine): cos(m*theta)*P_l^m(mu)
                    expr = cos(m * theta) * P_lm;
                else
                    %< Tesseral harmonic (sine): sin(m*theta)*P_l^m(mu)
                    expr = sin(m * theta) * P_lm;
                end
            end
            
            % Apply Condon-Shortley phase factor
            expr = (-sym(1))^m * expr;

            % Simplify the expression
            expr = simplify(expr);
        end
    end
    
    methods (Access = private)
        function [l, m] = getDegreesAndOrder(obj)
            % GETDEGREESANDORDER Determines degree l and order m from index.
            %
            %   Calculates the spherical harmonic degree and order based on
            %   the standard indexing convention for real-valued spherical
            %   harmonics.
            %
            % Inputs:
            %   obj - The SphericalHarmonicBasisFunction object
            %
            % Outputs:
            %   l - Degree of the spherical harmonic
            %   m - Order of the spherical harmonic
            
            idx = obj.index;
            
            %< Find degree l: largest integer such that l^2 < idx
            l = floor(sqrt(idx - 1));
            
            %< Calculate order m
            if idx == l^2 + 1
                %< Zonal harmonic (m = 0)
                m = 0;
            else
                %< Tesseral harmonic (m > 0)
                remainder = idx - l^2 - 1;
                m = ceil(remainder / 2);
            end
        end
        
        function P_lm = associatedLegendre(~, l, m, mu)
            % ASSOCIATEDLEGENDRE Computes associated Legendre polynomial
            % P_l^m(mu).
            %
            %   Calculates the associated Legendre polynomial using the
            %   standard formula involving derivatives of Legendre
            %   polynomials. For m = 0, returns the standard Legendre
            %   polynomial.
            %
            % Inputs:
            %   l - Degree (non-negative integer)
            %   m - Order (non-negative integer, m <= l)
            %   mu - Symbolic variable
            %
            % Outputs:
            %   P_lm - Associated Legendre polynomial
            
            if m == 0
                %< Standard Legendre polynomial
                P_lm = legendreP(l, mu);
            else
                %< Associated Legendre polynomial
                %< Use the relation: P_l^m(x) = (-1)^m * (1-x^2)^(m/2) * d^m/dx^m[P_l(x)]
                P_l = legendreP(l, mu);
                
                P_l_deriv = P_l;
                for k = 1:m
                    P_l_deriv = diff(P_l_deriv, mu);
                end
                P_lm = (-1)^m * (1 - mu^2)^(m/2) * P_l_deriv;
                P_lm = simplify(P_lm);
            end
        end
    end
end