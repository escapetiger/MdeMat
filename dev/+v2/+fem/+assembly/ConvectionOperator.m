classdef ConvectionOperator < fem.assembly.Assembly
    % CONVECTION

    properties
        volume
        flux
    end

    methods
        function obj = ConvectionOperator(context, bcType)
            obj@fem.assembly.Assembly(context);
            obj.volume = fem.assembly.VolumeAssembly(context);
            obj.flux = fem.assembly.FluxAssembly(context, bcType, []);
            if strcmpi(bcType, 'dirichlet')
                obj.flux.setFluxMode('upwind');
            end
        end

        function A = linear(obj, coe, vel)
            if nargin < 3 || isempty(vel), vel = coe; end

            n = obj.context.nGlobalDofs;
            A = sparse(n, n);
            for dim = 1:obj.context.nDims
                if coe(dim) == 0, continue; end
                if vel(dim) > 0
                    obj.flux.setFluxType('left');
                else
                    obj.flux.setFluxType('right');
                end
                C = obj.volume.scaleConstant(dim, coe(dim));
                T1 = obj.volume.assembleVolumePartial(dim, C);
                T2 = obj.flux.assembleFluxPartial(dim, C);
                T = [T1; T2];
                B = core.linalg.sparseFromTriplets(T, n, n);
                A = A + B;
            end
        end

        function A = linearBc(obj, coe, vel, f, varargin)
            n = obj.context.nGlobalDofs;
            A = cell(1, obj.context.nDims);
            for dim = 1:obj.context.nDims
                if coe(dim) == 0, continue; end
                C = obj.volume.scaleConstant(dim, coe(dim));
                if vel(dim) > 0
                    T = obj.flux.boundary.assembleLeftTrace(dim, C, f, varargin{:});
                else
                    T = obj.flux.boundary.assembleRightTrace(dim, C, f, varargin{:});
                end
                A{dim} = core.linalg.sparseFromTriplets(T, n, max(T(:, 2)));
            end

            for dim = obj.context.nDims:-1:2
                A{1} = A{1} + A{dim};
                A{dim} = [];
            end
            A = A{1};
        end
    end
end
