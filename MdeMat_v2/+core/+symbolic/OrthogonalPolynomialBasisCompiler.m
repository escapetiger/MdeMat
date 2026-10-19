classdef OrthogonalPolynomialBasisCompiler < core.symbolic.Compiler
    % ORTHOGONALPOLYNOMIALBASISCOMPILER Compiler for orthogonal polynomial basis functions.
    %
    %   OrthogonalPolynomialBasisCompiler specializes in compiling basis
    %   matrices for orthogonal polynomial families such as Legendre and
    %   Chebyshev polynomials. It creates and manages arrays of orthogonal
    %   polynomial basis functions with optional monic normalization.
    %
    %   Supported orthogonal polynomial families:
    %   - Legendre: orthogonal on [-1,1] with weight w(x) = 1
    %   - Chebyshev: orthogonal on [-1,1] with weight w(x) = 1/√(1-x²)
    %   Both families support standard and monic forms.
    %
    % Examples:
    %   % Create compiler and set Legendre functions
    %   compiler = core.symbolic.OrthogonalPolynomialBasisCompiler();
    %   compiler.setFunctions(3, [0, 1], 'legendre');
    %   
    %   % Compile and load basis matrix
    %   compiler.compile('legendre_basis', 'Legendre', 2);
    %   metadata = compiler.load('legendre_basis');
    %   
    %   % Chebyshev polynomials in monic form
    %   compiler.setFunctions(4, [-1, 1], 'monic_chebyshev');
    %   compiler.compile('chebyshev_monic', 'Chebyshev_Monic', 1);
    %
    % See Also:
    %   core.symbolic.Compiler, core.symbolic.LegendreBasisFunction, 
    %   core.symbolic.ChebyshevBasisFunction
    
    methods
        function obj = setFunctions(obj, degree, bbox, basisType)
            % SETFUNCTIONS Sets orthogonal polynomial basis functions.
            %
            %   obj = setFunctions(obj, degree, bbox, basisType) creates an
            %   array of orthogonal polynomial basis functions from degree
            %   0 to the specified maximum degree. Dispatches to the
            %   appropriate specialized method based on basis type.
            %
            % Inputs:
            %   obj - The OrthogonalPolynomialBasisCompiler object
            %   degree - Maximum degree of polynomials (non-negative integer)
            %   bbox - Bounding box [lower, upper] defining the interval
            %   basisType - Type of basis ('legendre', 'monic_legendre', 'chebyshev', 'monic_chebyshev')
            %
            % Outputs:
            %   obj - The OrthogonalPolynomialBasisCompiler object
            
            core.except.assert(degree >= 0 && isscalar(degree) && ...
                degree == floor(degree), 'InvalidInput', ...
                'Degree must be a non-negative integer.');
            
            core.except.assert(isvector(bbox) && length(bbox) == 2, ...
                'InvalidInput', 'Bounding box must be [lower, upper].');
            
            core.except.assert(bbox(1) < bbox(2), ...
                'InvalidInput', 'Lower bound must be less than upper bound.');
            
            supportedTypes = {'legendre', 'monic_legendre', 'chebyshev', 'monic_chebyshev'};
            core.except.assert(ismember(lower(basisType), supportedTypes), ...
                'InvalidInput', 'Supported basis types: %s', strjoin(supportedTypes, ', '));
            
            switch lower(basisType)
                case {'legendre', 'monic_legendre'}
                    obj.setLegendre(degree, bbox, basisType);
                case {'chebyshev', 'monic_chebyshev'}
                    obj.setChebyshev(degree, bbox, basisType);
                otherwise
                    core.except.assert(0, 'UnsupportedType', ...
                        'Basis type "%s" not implemented', basisType);
            end
        end
        
        function setLegendre(obj, degree, bbox, basisType)
            % SETLEGENDRE Sets Legendre polynomial basis functions.
            %
            %   setLegendre(obj, degree, bbox, basisType) creates an array
            %   of Legendre polynomial basis functions from degree 0 to the
            %   specified maximum degree. Supports both standard and monic
            %   forms with interval scaling to the specified bounding box.
            %
            % Inputs:
            %   obj - The OrthogonalPolynomialBasisCompiler object
            %   degree - Maximum degree of Legendre polynomials
            %   bbox - Bounding box [lower, upper]
            %   basisType - Type ('legendre' or 'monic_legendre')
            % 
            % Outputs:
            %   NULL
            
            x = obj.variables;
            monic = strcmpi(basisType, 'monic_legendre');
            
            obj.funcs = arrayfun( ...
                @(d) core.symbolic.LegendreBasisFunction(d, bbox, x, monic), 0:degree);
            
            fprintf('Set %d %s basis functions, degree 0 to %d on [%.6g, %.6g]\n', ...
                degree + 1, basisType, degree, bbox(1), bbox(2));
        end
        
        function setChebyshev(obj, degree, bbox, basisType)
            % SETCHEBYSHEV Sets Chebyshev polynomial basis functions.
            %
            %   setChebyshev(obj, degree, bbox, basisType) creates an array
            %   of Chebyshev polynomial basis functions from degree 0 to the
            %   specified maximum degree. Supports both standard and monic
            %   forms with interval scaling to the specified bounding box.
            %
            % Inputs:
            %   obj - The OrthogonalPolynomialBasisCompiler object
            %   degree - Maximum degree of Chebyshev polynomials
            %   bbox - Bounding box [lower, upper]
            %   basisType - Type ('chebyshev' or 'monic_chebyshev')
            %
            % Outputs:
            %   NULL
            
            x = obj.variables;
            monic = strcmpi(basisType, 'monic_chebyshev');
            
            obj.funcs = arrayfun( ...
                @(d) core.symbolic.ChebyshevBasisFunction(d, bbox, x, monic), 0:degree);
            
            fprintf('Set %d %s basis functions, degree 0 to %d on [%.6g, %.6g]\n', ...
                degree + 1, basisType, degree, bbox(1), bbox(2));
        end
    end
end