classdef PrimalDivergenceOperator < approx.assembly.Assembly
    % PRIMALDIVERGENCEOPERATOR Assembly for primal divergence operator.

    properties
        volume   % Volume assembly for interior terms (VolumeAssembly)
        flux  % Primal flux assembly (FluxAssembly)
        prmMat   % Primal operator matrices (cell array)
        trMat    % Trace matrices (cell array)
    end

    methods
        function obj = PrimalDivergenceOperator(space, operator, bcType, fluxType)
            % PRIMALDIVERGENCEOPERATOR Constructor for
            % PrimalDivergenceOperator.
            %
            %   obj = PrimalDivergenceOperator(space, operator, bcType,
            %   fluxType) creates a primal divergence operator with the
            %   specified components for Local Discontinuous Galerkin
            %   discretization.
            %
            % Inputs:
            %   space - MeshSpace object
            %   operator - ElementOperator object
            %   bcType - Boundary condition type: {'periodic', 'dirichlet'}
            %   fluxType - Flux type: {'left', 'right', 'central'}
            %
            % Outputs:
            %   obj - Constructed PrimalDivergenceOperator object

            obj@approx.assembly.Assembly(space, operator);
            obj.volume = approx.assembly.VolumeAssembly(space, operator);
            obj.flux = approx.assembly.FluxAssembly(space, operator, bcType, fluxType);
            if strcmpi(bcType, 'dirichlet')
                obj.flux.setFluxMode('primal');
            end
            obj.prmMat = cell(1, obj.space.nDims);
            obj.trMat = cell(1, obj.space.element.nFluxes);
        end

        function A = linear(obj, coe)
            % LINEAR Assemble auxiliary linear divergence matrix.
            %
            %   A = linear(obj, coe) assembles the auxiliary divergence of
            %   linear field.
            %
            % Inputs:
            %   obj - The AuxiliaryDivergenceOperator object
            %   coe - Constant coefficients (nDims x 1 vector)
            %
            % Outputs:
            %   A - Assembled matrix (sparse)

            n = obj.space.nGlobalDofs;

            for dim = 1:obj.space.nDims
                fp = isempty(obj.prmMat{dim});

                if fp
                    C = obj.volume.scaleConstant(dim, 1);
                    T1 = obj.volume.assembleVolumePartial(dim, C);
                    T2 = obj.flux.assembleFluxPartial(dim, C);
                    T = [T1; T2];
                    obj.prmMat{dim} = core.linalg.sparseFromTriplets(T, n, n);
                end
            end

            A = sparse(n, n);
            for dim = 1:obj.space.nDims
                if coe(dim) == 0, continue; end
                A = A + coe(dim) * obj.prmMat{dim};
            end
        end
    
        function A = linearBc(obj, coe, f, varargin)
            % LINEARBC Assemble boundary condition for primal linear
            % divergence.
            %
            %   A = linearBc(obj, coe) assembles the boundary condition for
            %   primal linear divergence of linear field.
            %
            % Inputs:
            %   obj - The PrimalDivergenceOperator object
            %   coe - Constant coefficients (nDims x 1 vector)
            %   f - Boundary condition handle
            %   varargin - Input arguments passed to f
            %
            % Outputs:
            %   A - Assembled matrix (sparse)

            n = obj.space.nGlobalDofs;
            
            for dim = 1:obj.space.nDims
                C = obj.volume.scaleConstant(dim, 1);
                T = obj.flux.boundary.assembleTrace(2*dim-1, C, f, varargin{:});
                obj.trMat{2*dim-1} = core.linalg.sparseFromTriplets(T, n, max(T(:, 2)));
                T = obj.flux.boundary.assembleTrace(2*dim, C, f, varargin{:});
                obj.trMat{2*dim} = core.linalg.sparseFromTriplets(T, n, max(T(:, 2)));
            end

            A = coe(1) * (obj.trMat{1} + obj.trMat{2});
            for dim = obj.space.nDims:-1:2
                A = A + coe(dim) * (obj.trMat{2*dim} - obj.trMat{2*dim-1});
            end
        end
    
    end
end