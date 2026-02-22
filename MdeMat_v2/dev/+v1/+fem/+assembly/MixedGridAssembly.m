classdef MixedGridAssembly < fem.assembly.EulerianGridAssembly
    % MIXEDGRIDASSEMBLY Grid-based assembly for mixed finite element
    % systems.
    %
    %   MixedGridAssembly extends the Eulerian grid assembly with support
    %   for different flux types (left-biased, right-biased, central) used
    %   in mixed finite element methods and discontinuous Galerkin schemes.
    %
    %   This class provides the foundation for assembling operators in
    %   mixed formulations where different variables may require different
    %   flux treatments. It supports various numerical flux choices
    %   commonly used in discontinuous methods.
    %
    % Examples:
    %   % Create mixed assembly with central flux and Dirichlet BC
    %   assembly = MixedGridAssembly(fe, mesh, op, 1, 3);
    %   
    %   % Assemble divergence operator with dimension-specific coefficients
    %   coeffs = [1.0, 0.5]; % coefficients for each dimension
    %   divMatrix = assembly.divergence(coeffs);
    %
    %   % Assemble partial derivative in specific dimension
    %   partialMatrix = assembly.partial(1, 2.0);
    %
    % See also:
    %   fem.assembly.EulerianGridAssembly, fem.assembly.GridAssembly
    
    properties
        fluxType % Flux type (1=left-biased, 2=right-biased, 3=central)
    end

    methods
        function obj = MixedGridAssembly(fe, mesh, op, bcType, fluxType)
            % MIXEDGRIDASSEMBLY Constructor for MixedGridAssembly.
            %
            %   obj = MixedGridAssembly(fe, mesh, op, bcType, fluxType)
            %   creates an assembly for mixed finite element systems with
            %   specified flux type for discontinuous Galerkin methods.
            %
            % Inputs:
            %   fe - Finite element object
            %   mesh - Grid mesh object (must be approx.mesh.Grid)
            %   op - Finite element operator object
            %   bcType - Boundary condition type (0=periodic, 1=Dirichlet)
            %   fluxType - Flux type (1=left-biased, 2=right-biased, 3=central)
            %
            % Outputs:
            %   obj - Constructed MixedGridAssembly object

            core.except.assert(isa(mesh, 'approx.mesh.Grid'), ...
                'InvalidInput', ...
                'Finite element space must be based on grid.')

            obj@fem.assembly.EulerianGridAssembly(fe, mesh, op, bcType);
            obj.fluxType = fluxType;
        end

        function A = divergence(obj, c)
            % DIVERGENCE Assemble divergence operator with coefficients.
            %
            %   A = divergence(obj, c) creates the divergence operator
            %   matrix by summing partial derivative operators across all
            %   spatial dimensions, each weighted by the corresponding
            %   coefficient.
            %
            % Inputs:
            %   obj - The MixedGridAssembly object
            %   c - Coefficient vector (one per spatial dimension)
            %
            % Outputs:
            %   A - Assembled divergence operator sparse matrix

            n = obj.nGlobalDofs;
            A = sparse(n, n);
            for d = 1:obj.nDims
                if c(d) == 0, continue; end
                T = obj.partial(d, c(d));
                A = A + T;
            end
        end

        function A = partial(obj, d, c)
            % PARTIAL Assemble partial derivative operator with
            % coefficients.
            %
            %   A = partial(obj, d, c) creates the partial derivative
            %   operator matrix for the specified dimension, including
            %   volume and flux contributions based on the finite element
            %   type and flux choice.
            %
            % Inputs:
            %   obj - The MixedGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient value for this dimension
            %
            % Outputs:
            %   A - Assembled partial derivative sparse matrix

            ne = obj.nTotalElements;

            %< ------------------------------------------------------------
            %< COEFFICIENTS
            %< ------------------------------------------------------------
            if isa(obj.mesh, 'approx.mesh.UniformGrid')
                h = obj.mesh.spacings{d};
                c = ones(1, ne) .* c ./ h;
            end

            if isa(obj.mesh, 'approx.mesh.NonuniformGrid')
                h = obj.mesh.spacings;
                [h{:}] = ndgrid(h{:});
                c = c ./ h{d}(:);
            end

            c = c(:);

            %< ------------------------------------------------------------
            %< VOLUME CONTRIBUTION
            %< ------------------------------------------------------------
            T1 = obj.assembleVolumeTerms(d, c);

            %< ------------------------------------------------------------
            %< FLUX CONTRIBUTION
            %< ------------------------------------------------------------
            if isa(obj.fe, 'fem.element.L2FiniteElement')
                switch obj.fluxType
                    case 1
                        T2 = obj.assembleLeftFluxTerms(d, c);
                    case 2
                        T2 = obj.assembleRightFluxTerms(d, c);
                    case 3
                        T2 = obj.assembleCentralFluxTerms(d, c);
                end
            end

            %< ------------------------------------------------------------
            %< TRIPLETS TO SPARSE MATRIX
            %< ------------------------------------------------------------
            n = obj.nGlobalDofs;
            T = [T1; T2];
            A = obj.sparseFromTriplets(T, n, n);
        end
    end

    methods (Abstract, Access = protected)
        % ASSEMBLEINNERLEFTFLUXDIRICHLETBOUNDARYTERMS Inner contribution to
        % a left flux along a specific dimension on the boundary with
        % Dirichlet boundary conditions.
        T = assembleInnerLeftFluxDirichletBoundaryTerms(obj, d, c)

        % ASSEMBLEOUTERLEFTFLUXDIRICHLETBOUNDARYTERMS Outer contribution to
        % a left flux along a specific dimension on the boundary with
        % Dirichlet boundary conditions.
        T = assembleOuterLeftFluxDirichletBoundaryTerms(obj, d, c)
        
        % ASSEMBLEINNERRIGHTFLUXDIRICHLETBOUNDARYTERMS Inner contribution
        % to a right flux along a specific dimension on the boundary with
        % Dirichlet boundary conditions.
        T = assembleInnerRightFluxDirichletBoundaryTerms(obj, d, c)
        
        % ASSEMBLEOUTERRIGHTFLUXDIRICHLETBOUNDARYTERMS Outer contribution
        % to a right flux along a specific dimension on the boundary with
        % Dirichlet boundary conditions.
        T = assembleOuterRightFluxDirichletBoundaryTerms(obj, d, c)
    end
end