classdef CompiledOrthogonalBasis < mdemat.common.approximation.CompiledFunction
    % COMPILEDORTHOGONALBASIS Compiled orthogonal polynomial basis
    % functions.
    %
    %   CompiledOrthogonalBasis represents a set of orthogonal polynomial
    %   basis functions (such as Legendre or Chebyshev polynomials) with
    %   automatic metadata loading from the mdemat.common.symbolic compilation
    %   system.
    %
    %   The class automatically discovers appropriate metadata files based
    %   on the specified basis type, providing efficient evaluation of both
    %   function values and derivatives through precompiled function
    %   handles.
    %
    % See Also:
    %   mdemat.common.approximation.CompiledFunction
    
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
                case 'normal_spherical_harmonic'
                    nDims = 2;
            end
            
            obj@mdemat.common.approximation.CompiledFunction(nDims, nCodims, basisName);
            obj.basisName = basisName;
        end
    end
end