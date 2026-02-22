classdef CompiledInterpolationBasis < core.function.CompiledFunction
    % COMPILEDINTERPOLATIONBASIS Compiled interpolation basis functions.
    %
    %   CompiledInterpolationBasis represents a set of interpolation basis
    %   functions (such as Lagrange or Hermite polynomials) with automatic
    %   metadata loading from the core.symbolic compilation system.
    %
    %   The class automatically discovers appropriate metadata files based
    %   on the specified basis type, interval type, node type, and number
    %   of basis functions, providing efficient evaluation of both function
    %   values and derivatives through precompiled function handles.
    %
    % See Also:
    %   core.function.CompiledFunction

    properties
        basisName % Basis name
    end

    methods
        function obj = CompiledInterpolationBasis(nCodims, basisName)
            % COMPILEDINTERPOLATIONBASIS Constructor for
            % CompiledInterpolationBasis.
            %
            %   obj = CompiledInterpolationBasis(nCodims, basisName)
            %   creates a new CompiledInterpolationBasis object with the
            %   specified configuration. Validates all input parameters and
            %   initializes the parent CompiledFunction.

            arguments
                nCodims (1,1) {mustBeNonnegative, mustBeInteger}
                basisName {mustBeMember(basisName, {
                    'unit_gauss_legendre_lagrange', ...
                    'unit_gauss_lobatto_lagrange'})}
            end

            basisName = lower(basisName);
            nDims = 1;
            basisName = sprintf('%s_%d', basisName, nCodims);

            obj@core.function.CompiledFunction(nDims, nCodims, basisName);
            obj.basisName = basisName;
        end
    end
end