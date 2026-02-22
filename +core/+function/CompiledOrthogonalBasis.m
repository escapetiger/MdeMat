classdef CompiledOrthogonalBasis < core.function.CompiledFunction
    % COMPILEDORTHOGONALBASIS Compiled orthogonal polynomial basis
    % functions.
    %
    %   CompiledOrthogonalBasis represents a set of orthogonal polynomial
    %   basis functions (such as Legendre or Chebyshev polynomials) with
    %   automatic metadata loading from the core.symbolic compilation
    %   system.
    %
    %   The class automatically discovers appropriate metadata files based
    %   on the specified basis type, providing efficient evaluation of both
    %   function values and derivatives through precompiled function
    %   handles.
    %
    % See Also:
    %   core.function.CompiledFunction
    
    properties
        basisName % Basis name
    end
    
    methods
        function obj = CompiledOrthogonalBasis(nCodims, basisName)
            % COMPILEDORTHOGONALBASIS Constructor for
            % CompiledOrthogonalBasis.
            %
            %   obj = CompiledOrthogonalBasis(nCodims, basisName) creates a
            %   new CompiledOrthogonalBasis object with the specified
            %   configuration. Validates all input parameters and
            %   initializes the parent CompiledFunction.
            
            arguments
                nCodims (1,1) {mustBeNonnegative, mustBeInteger}
                basisName {mustBeMember(basisName, { ...
                    'monic_unit_legendre', ...
                    'normal_legendre', ...
                    'normal_hermite', ...
                    'normal_fourier', ...
                    'normal_spherical_harmonic'})}
            end
            
            switch basisName
                case {'monic_unit_legendre', 'normal_legendre', 'normal_fourier', 'normal_hermite'}
                    nDims = 1;
                case 'spherical_harmonic'
                    nDims = 2;
            end
            
            obj@core.function.CompiledFunction(nDims, nCodims, basisName);
            obj.basisName = basisName;
        end
    end
end