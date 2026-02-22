classdef FourierBasisFunction < core.symbolic.BasisFunction
    % FOURIERBASISFUNCTION Real-valued Fourier basis functions.
    %
    %   FourierBasisFunction represents real-valued Fourier basis functions
    %   on a specified interval. The basis functions include a constant
    %   term and sine/cosine harmonics. The functions are scaled from the
    %   canonical interval \f$[0, 2\pi]\f$ to the given bounding box.
    %
    %   The first five Fourier basis functions with wavelength \f$\omega =
    %   \frac{2\pi}{b-a}\f$ are: 
    %   \f[
    %        \phi_0(x) = 1,
    %        \phi_1(x) = \cos(\omega x), 
    %        \phi_2(x) = \sin(\omega x),
    %        \phi_3(x) = \cos(2\omega x),
    %        \phi_4(x) = \sin(2\omega x).
    %   \f]
    %
    % Examples:
    %   % Constant function
    %   basis1 = core.symbolic.FourierBasisFunction(1, [0, 2*pi], sym('x'));
    %
    %   % First harmonic cosine and sine
    %   basis2 = core.symbolic.FourierBasisFunction(2, [0, 2*pi], sym('x'));
    %   basis3 = core.symbolic.FourierBasisFunction(3, [0, 2*pi], sym('x'));
    %
    % See also:
    %   core.symbolic.BasisFunction
    
    methods
        function obj = FourierBasisFunction(index, bbox, variables)
            % FOURIERBASISFUNCTION Constructor for FourierBasisFunction.
            %
            %   obj = FourierBasisFunction(index, bbox) creates a new 
            %   FourierBasisFunction object with the specified index and
            %   interval using the default variable sym('x').
            %   
            %   obj = FourierBasisFunction(index, bbox, variables) creates
            %   a new FourierBasisFunction object with specified variables.
            %
            % Inputs:
            %   index - Index of the basis function (positive integer)
            %   bbox - Interval [a, b] for the function domain
            %   variables - Symbol(s) (optional, default: {sym('x')})
            %
            % Outputs:
            %   obj - Constructed FourierBasisFunction object
            
            obj@core.symbolic.BasisFunction(index, bbox, variables);
            obj.expression = obj.generate();
        end
    end

    methods (Access = protected)
        function expr = generate(obj)
            % GENERATE Generates the Fourier basis function.
            %
            %   expr = generate(obj) creates the symbolic expression for
            %   the Fourier basis function based on the index. Handles
            %   constant, cosine, and sine terms with appropriate scaling
            %   to the specified interval.
            %
            % Inputs:
            %   obj - The FourierBasisFunction object
            %
            % Outputs:
            %   expr - Symbolic expression for the Fourier basis function
            
            x = obj.variables{1};
            a = obj.bbox(1);
            b = obj.bbox(2);
            t = 2 * pi * (x - a) / (b - a);
            
            if obj.index == 1
                expr = sym(1);
            else
                harmonicNum = ceil((obj.index - 1) / 2);
                
                if mod(obj.index, 2) == 0
                    expr = cos(harmonicNum * t);
                else
                    expr = sin(harmonicNum * t);
                end
            end
            
            expr = simplify(expr);
        end
    end
end