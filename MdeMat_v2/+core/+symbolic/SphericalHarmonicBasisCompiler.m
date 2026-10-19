classdef SphericalHarmonicBasisCompiler < core.symbolic.Compiler
    % SPHERICALHARMONICBASISCOMPILER Compiler for spherical harmonic basis functions.
    %
    %   SphericalHarmonicBasisCompiler specializes in compiling basis
    %   matrices for real-valued spherical harmonic basis functions. It
    %   creates and manages arrays of spherical harmonic basis functions
    %   including zonal and tesseral harmonics.
    %
    % Examples:
    %   % Create compiler and set functions
    %   compiler = core.symbolic.SphericalHarmonicBasisCompiler();
    %   compiler.setFunctions(8, 'spherical_harmonic');
    %
    %   % Compile and load basis matrix
    %   compiler.compile('spherical_harmonic_basis', 'SphericalHarmonic', 1);
    %   metadata = compiler.load('spherical_harmonic_basis');
    %
    % See also:
    %   core.symbolic.Compiler, core.symbolic.SphericalHarmonicBasisFunction
    
    methods
        function obj = setFunctions(obj, maxIndex, basisType)
            % SETFUNCTIONS Sets spherical harmonic basis functions.
            %
            %   obj = setFunctions(obj, maxIndex) creates an array of
            %   spherical harmonic basis functions from index 1 to maxIndex
            %   using the default basis type.
            %   
            %   obj = setFunctions(obj, maxIndex, basisType) creates
            %   spherical harmonic basis functions with explicit type
            %   specification.
            %
            % Inputs:
            %   obj - The SphericalHarmonicBasisCompiler object
            %   maxIndex - Maximum index for basis functions (positive integer)
            %   basisType - Type of basis (optional, default: 'spherical_harmonic')
            %
            % Outputs:
            %   obj - The SphericalHarmonicBasisCompiler object
            
            core.except.assert(maxIndex >= 1, ...
                'InvalidInput', 'Maximum index must be positive.');
            
            if nargin < 3
                basisType = 'spherical_harmonic';
            end
            
            supportedTypes = {'spherical_harmonic'};
            core.except.assert(ismember(lower(basisType), supportedTypes), ...
                'InvalidInput', 'Supported basis types: %s', strjoin(supportedTypes, ', '));
            
            obj.funcs = arrayfun(@(i) core.symbolic.SphericalHarmonicBasisFunction(i), ...
                1:maxIndex);
            
            fprintf('Set %d %s basis functions, index 1 to %d\n', ...
                maxIndex, basisType, maxIndex);
        end
    end
end