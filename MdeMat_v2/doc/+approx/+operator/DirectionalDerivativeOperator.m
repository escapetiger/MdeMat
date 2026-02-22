classdef DirectionalDerivativeOperator < approx.operator.MeshOperator
    % DIRECTIONALDERIVATIVEOPERATOR

    properties (Constant)
        VALID_BC_TYPES = {'periodic', 'dirichlet'}
        VALID_FLUX_TYPES = {'upwind', 'central', 'outward', 'inward'}
        VALID_PERIODIC_FLUX_MODES = {'periodic'}
        VALID_DIRICHLET_FLUX_MODES = {'extrapolate', 'truncate', 'inflow'}
    end

    properties
        volume % 
        flux % 
        bcType % Boundary condition type
        fluxType % Flux type
        fluxMode % Flux mode
    end

    methods
        function obj = DirectionalDerivativeOperator(space, bcType, fluxType, fluxMode)
            obj@approx.operator.MeshOperator(space);
            obj.setBcType(bcType);
            obj.setFluxType(fluxType);
            obj.setFluxMode(fluxMode);
        end

        function obj = setBcType(obj, bcType)
            core.except.assert(obj.isValidBcType(bcType), ...
                'InvalidInput', 'BC type is not supported.');

            core.except.assert(obj.isCompitable(bcType, obj.fluxMode), ...
                'InvalidInput', 'Flux mode is incompatiable with BC.');

            obj.bcType = bcType;
        end

        function obj = setFluxType(obj, fluxType)
            core.except.assert(isempty(fluxType) || ...
                ismember(fluxType, obj.VALID_FLUX_TYPES), ...
                'InvalidInput', 'Flux type is not supported.');

            obj.fluxType = fluxType;
        end

        function obj = setFluxMode(obj, fluxMode)
            core.except.assert(obj.isValidFluxMode(fluxMode), ...
                'InvalidInput', 'Flux mode is not supported.');

            core.except.assert(obj.isCompitable(obj.bcType, fluxMode), ...
                'InvalidInput', 'Flux mode is incompatiable with BC.');

            obj.fluxMode = fluxMode;
        end

        function D = addBilinearVolume(obj, coe)
            d = obj.space.nDims;
            if ~isscalar(coe)
                C = reshape(coe, 1, 1, d, []);
            else
                C = coe;
            end

            invJ = obj.space.mesh.computeElementInverseJacobians();
            W = diag(obj.space.element.volume.weights);
            U = obj.space.element.volume.derivatives(:, :, 1:d);
            V = obj.space.element.volume.values;
            A = pagemtimes(pagemtimes(U, W), V.');
            A = reshape(A, [], d) * reshape(invJ, d, []);
            A = reshape(A, size(U, 1), size(V, 1), d, size(invJ, 3));
            A = A .* C;
            D = obj.addBilinear([], [], [], [], A);
        end

        function F = addBilinearInteriorFlux(obj, coe, vel)
            if nargin < 3
                core.except.assert(~strcmpi(obj.fluxType, 'upwind'), ...
                    'InvalidInput', 'Upwind flux requires velocity.');
                vel = [];
            end
        end

        function F = addBilinearBoundaryFlux(obj, coe, vel)
            if nargin < 3
                core.except.assert(~strcmpi(obj.fluxType, 'upwind'), ...
                    'InvalidInput', 'Upwind flux requires velocity.');
                vel = [];
            end
        end

        function F = addLinearBoundaryFlux(obj, coe, vel)
            if nargin < 3
                core.except.assert(~strcmpi(obj.fluxType, 'upwind'), ...
                    'InvalidInput', 'Upwind flux requires velocity.');
                vel = [];
            end
        end
    
        function alpha = computeStablizer(obj, coe, vel)
            switch obj.fluxType
                case 'upwind'
                case 'central'
                case 'outward'
                case 'inflow'
            end
        end
    end

    methods (Access = protected)
        function F = addBilinearInteriorAverage(obj, coe)
        end

        function F = addBiliearInteriorJump(obj, coe)
        end

        function F = addBilinearBoundaryAverage(obj, coe)
        end

        function F = addBilinearBoundaryJump(obj, coe)
        end

        function F = addLinearBoundaryAverage(obj, coe)
        end

        function F = addLinearBoundaryJump(obj, coe)
        end
    end

    methods (Access = private)
        function tf = isValidBcType(obj, bcType)
            tf = ismember(bcType, obj.VALID_BC_TYPES);
        end

        function tf = isValidFluxType(obj, fluxType)
            tf = isempty(fluxType) || ismember(fluxType, obj.VALID_FLUX_TYPES);
        end

        function tf = isValidFluxMode(obj, fluxMode)
            VALID_FLUX_MODES = [obj.VALID_PERIODIC_FLUX_MODES, ...
                obj.VALID_DIRICHLET_FLUX_MODES];
            tf = isempty(fluxMode) || ismember(fluxMode, VALID_FLUX_MODES);
        end

        function tf = isCompatible(obj, bcType, fluxMode)
            switch bcType
                case 'periodic'
                    tf = ismember(fluxMode, obj.VALID_PERIODIC_FLUX_MODES);
                case 'dirichlet'
                    tf = ismember(fluxMode, obj.VALID_DIRICHLET_FLUX_MODES);
                otherwise
                    tf = false;
            end
        end
    end
end