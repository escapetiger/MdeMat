classdef EulerianGridAssembly < fem.assembly.GridAssembly
    % EULERIANGRIDASSEMBLY Eulerian grid-based assembly.
    %
    %   EulerianGridAssembly provides assembly methods for finite element 
    %   operators on structured grid meshes using an Eulerian approach.
    %   This class handles volume and flux contributions with support for
    %   different boundary condition types including periodic and Dirichlet
    %   conditions.
    %
    %   The assembly supports various flux formulations (left-biased,
    %   right-biased, and central) and handles both interior and boundary
    %   contributions for each type.
    %
    % Examples:
    %   % Create Eulerian grid assembly
    %   assembly = EulerianGridAssembly(fe, mesh, op, 0); % Periodic BC
    %   
    %   % Assemble different flux types
    %   T_left = assembly.assembleLeftFluxTerms(1, coeffs);
    %   T_central = assembly.assembleCentralFluxTerms(1, coeffs);
    %
    % See also:
    %   fem.assembly.GridAssembly, fem.assembly.Assembly

    methods (Access = protected)
        function T = assembleVolumeTerms(obj, d, c)
            % ASSEMBLEVOLUMETERMS Assemble volume integral terms.
            %
            %   T = assembleVolumeTerms(obj, d, c) creates triplets for the
            %   volume contribution of differential operators along the 
            %   specified spatial dimension.
            %
            % Inputs:
            %   obj - The EulerianGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector, one per element
            %
            % Outputs:
            %   T - Triplet matrix [row, col, value] for volume terms

            ne = obj.nTotalElements;
            np = obj.nLocalDofs;
            P = obj.fe.projector;
            G = obj.op.gradient;
            F = G.volumeData(d).matrix;

            T = zeros(np^2*ne, 3);
            [i, j, v] = find(-P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (0:ne - 1)*np, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (0:ne - 1)*np, j(:)), [], 1);
            T(k, 3) = kron(c(:), v(:));
        end

        function T = assembleLeftFluxTerms(obj, d, c)
            % ASSEMBLELEFTFLUXTERMS Assemble left-biased flux terms.
            %
            %   T = assembleLeftFluxTerms(obj, d, c) combines interior and
            %   boundary contributions for left-biased flux terms along the
            %   specified dimension.
            %
            % Inputs:
            %   obj - The EulerianGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector, one per element
            %
            % Outputs:
            %   T - Combined triplet matrix for left flux terms

            T1 = obj.assembleInnerLeftFluxInteriorTerms(d, c);
            T2 = obj.assembleInnerLeftFluxBoundaryTerms(d, c);
            T3 = obj.assembleOuterLeftFluxInteriorTerms(d, c);
            T4 = obj.assembleOuterLeftFluxBoundaryTerms(d, c);
            T = [T1; T2; T3; T4];
        end

        function T = assembleRightFluxTerms(obj, d, c)
            % ASSEMBLYIGHTFLUXTERMS Assemble right-biased flux terms.
            %
            %   T = assembleRightFluxTerms(obj, d, c) combines interior and
            %   boundary contributions for right-biased flux terms along
            %   the specified dimension.
            %
            % Inputs:
            %   obj - The EulerianGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector, one per element
            %
            % Outputs:
            %   T - Combined triplet matrix for right flux terms

            T1 = obj.assembleInnerRightFluxInteriorTerms(d, c);
            T2 = obj.assembleInnerRightFluxBoundaryTerms(d, c);
            T3 = obj.assembleOuterRightFluxInteriorTerms(d, c);
            T4 = obj.assembleOuterRightFluxBoundaryTerms(d, c);
            T = [T1; T2; T3; T4];
        end

        function T = assembleCentralFluxTerms(obj, d, c)
            % ASSEMBLECENTRALFLUXTERMS Assemble central flux terms.
            %
            %   T = assembleCentralFluxTerms(obj, d, c) combines left and
            %   right flux contributions with equal weighting to create
            %   central difference-type operators.
            %
            % Inputs:
            %   obj - The EulerianGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector, one per element
            %
            % Outputs:
            %   T - Combined triplet matrix for central flux terms

            T1 = obj.assembleInnerLeftFluxTerms(d, c/2);
            T2 = obj.assembleOuterLeftFluxInteriorTerms(d, c/2);
            T3 = obj.assembleOuterLeftFluxBoundaryTerms(d, c/2);
            T4 = obj.assembleInnerRightFluxTerms(d, c/2);
            T5 = obj.assembleOuterRightFluxInteriorTerms(d, c/2);
            T6 = obj.assembleOuterRightFluxBoundaryTerms(d, c/2);
            T = [T1; T2; T3; T4; T5; T6];
        end

        function T = assembleInnerLeftFluxTerms(obj, d, c)
            % ASSEMBLEINNERLEFTFLUXTERMS Assemble inner left flux
            % contribution.
            %
            %   T = assembleInnerLeftFluxTerms(obj, d, c) computes the
            %   inner contribution to a left flux along the specified
            %   dimension.
            %
            % Inputs:
            %   obj - The EulerianGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector, one per element
            %
            % Outputs:
            %   T - Triplet matrix for inner left flux terms

            ne = obj.nTotalElements;
            np = obj.nLocalDofs;
            P = obj.fe.projector;
            G = obj.op.gradient;

            T = zeros(np^2*ne, 3);

            F = G.fluxData(1, 2*d).matrix;
            [i, j, v] = find(P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (0:ne - 1)*np, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (0:ne - 1)*np, j(:)), [], 1);
            T(k, 3) = kron(c(:), v(:));
        end

        function T = assembleInnerRightFluxTerms(obj, d, c)
            % ASSEMBLEINNERRIGHTFLUXTERMS Assemble inner right flux
            % contribution.
            %
            %   T = assembleInnerRightFluxTerms(obj, d, c) computes the
            %   inner contribution to a right flux along the specified
            %   dimension.
            %
            % Inputs:
            %   obj - The EulerianGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector, one per element
            %
            % Outputs:
            %   T - Triplet matrix for inner right flux terms

            ne = obj.nTotalElements;
            np = obj.nLocalDofs;
            P = obj.fe.projector;
            G = obj.op.gradient;

            T = zeros(np^2*ne, 3);

            F = G.fluxData(1, 2*d-1).matrix;
            [i, j, v] = find(-P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (0:ne - 1)*np, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (0:ne - 1)*np, j(:)), [], 1);
            T(k, 3) = kron(c(:), v(:));
        end

        function T = assembleInnerLeftFluxInteriorTerms(obj, d, c)
            % ASSEMBLEINNERLEFTFLUXINTERIORTERMS Assemble inner left flux
            % interior terms.
            %
            %   T = assembleInnerLeftFluxInteriorTerms(obj, d, c) computes
            %   the inner contribution to a left flux along the specified
            %   dimension for elements in the domain interior.
            %
            % Inputs:
            %   obj - The EulerianGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector, one per element
            %
            % Outputs:
            %   T - Triplet matrix for interior contribution

            np = obj.nLocalDofs;
            P = obj.fe.projector;
            G = obj.op.gradient;

            m = obj.mesh.allElementMultiIndices;
            e = m(:, d) > 1 & m(:, d) < obj.mesh.resolution(d);
            m = m(e, :);
            l = obj.mesh.indexer.multiToLinear(m, 0);
            ne = numel(l);

            T = zeros(np^2*ne, 3);

            F = G.fluxData(1, 2*d).matrix;
            [i, j, v] = find(P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l(:).' - 1)*np, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l(:).' - 1)*np, j(:)), [], 1);
            T(k, 3) = kron(c(l), v(:));
        end

        function T = assembleInnerRightFluxInteriorTerms(obj, d, c)
            % ASSEMBLEINNERRIGHTFLUXINTERIORTERMS Assemble inner right flux
            % interior terms.
            %
            %   T = assembleInnerRightFluxInteriorTerms(obj, d, c) computes
            %   the inner contribution to a right flux along the specified
            %   dimension for elements in the domain interior.
            %
            % Inputs:
            %   obj - The EulerianGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector, one per element
            %
            % Outputs:
            %   T - Triplet matrix for interior contribution

            np = obj.nLocalDofs;
            P = obj.fe.projector;
            G = obj.op.gradient;

            m = obj.mesh.allElementMultiIndices;
            e = m(:, d) > 1 & m(:, d) < obj.mesh.resolution(d);
            ne = numel(find(e));
            m = m(e, :);
            l = obj.mesh.indexer.multiToLinear(m, 0);

            T = zeros(np^2*ne, 3);

            F = G.fluxData(1, 2*d-1).matrix;
            [i, j, v] = find(-P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l(:).' - 1)*np, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l(:).' - 1)*np, j(:)), [], 1);
            T(k, 3) = kron(c(l), v(:));
        end
        
        function T = assembleInnerLeftFluxBoundaryTerms(obj, d, c)
            % ASSEMBLEINNERLEFTFLUXBOUNDARYTERMS Assemble inner left flux
            % boundary terms.
            %
            %   T = assembleInnerLeftFluxBoundaryTerms(obj, d, c) computes
            %   the inner contribution to a left flux along the specified
            %   dimension on the boundary, dispatching to appropriate
            %   boundary condition methods.
            %
            % Inputs:
            %   obj - The EulerianGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector, one per element
            %
            % Outputs:
            %   T - Triplet matrix for boundary contribution

            switch obj.bcType
                case 0
                    T = obj.assembleInnerLeftFluxPeriodicBoundaryTerms(d, c);
                case 1
                    T = obj.assembleInnerLeftFluxDirichletBoundaryTerms(d, c);
            end
        end

        function T = assembleInnerRightFluxBoundaryTerms(obj, d, c)
            % ASSEMBLEINNERRIGHTFLUXBOUNDARYTERMS Assemble inner right flux
            % boundary terms.
            %
            %   T = assembleInnerRightFluxBoundaryTerms(obj, d, c) computes
            %   the inner contribution to a right flux along the specified
            %   dimension on the boundary, dispatching to appropriate
            %   boundary condition methods.
            %
            % Inputs:
            %   obj - The EulerianGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector, one per element
            %
            % Outputs:
            %   T - Triplet matrix for boundary contribution

            switch obj.bcType
                case 0
                    T = obj.assembleInnerRightFluxPeriodicBoundaryTerms(d, c);
                case 1
                    T = obj.assembleInnerRightFluxDirichletBoundaryTerms(d, c);
            end
        end
        
        function T = assembleOuterLeftFluxInteriorTerms(obj, d, c)
            % ASSEMBLEOUTERLEFTFLUXINTERIORTERMS Assemble outer left flux
            % interior terms.
            %
            %   T = assembleOuterLeftFluxInteriorTerms(obj, d, c) computes
            %   the outer contribution to a left flux along the specified
            %   dimension for elements in the domain interior. This handles
            %   coupling between adjacent elements.
            %
            % Inputs:
            %   obj - The EulerianGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector, one per element
            %
            % Outputs:
            %   T - Triplet matrix for outer contribution

            np = obj.nLocalDofs;
            P = obj.fe.projector;
            G = obj.op.gradient;

            m = obj.mesh.allElementMultiIndices;
            e = m(:, d) > 1 & m(:, d) < obj.mesh.resolution(d);
            ne = numel(find(e));
            m = m(e, :);
            l1 = obj.mesh.indexer.multiToLinear(m, 0);
            m(:, d) = m(:, d) - 1;
            l2 = obj.mesh.indexer.multiToLinear(m, 0);

            T = zeros(np^2*ne, 3);

            F = G.fluxData(2, 2*d-1).matrix;
            [i, j, v] = find(-P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l1(:).' - 1)*np, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l2(:).' - 1)*np, j(:)), [], 1);
            T(k, 3) = kron(c(l1), v(:));
        end

        function T = assembleOuterRightFluxInteriorTerms(obj, d, c)
            % ASSEMBLEOUTERRIGHTFLUXINTERIORTERMS Assemble outer right flux
            % interior terms.
            %
            %   T = assembleOuterRightFluxInteriorTerms(obj, d, c) computes
            %   the outer contribution to a right flux along the specified
            %   dimension for elements in the domain interior. This handles
            %   coupling between adjacent elements.
            %
            % Inputs:
            %   obj - The EulerianGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector, one per element
            %
            % Outputs:
            %   T - Triplet matrix for outer contribution

            np = obj.nLocalDofs;
            P = obj.fe.projector;
            G = obj.op.gradient;

            m = obj.mesh.allElementMultiIndices;
            e = m(:, d) > 1 & m(:, d) < obj.mesh.resolution(d);
            ne = numel(find(e));
            m = m(e, :);
            l1 = obj.mesh.indexer.multiToLinear(m, 0);
            m(:, d) = m(:, d) + 1;
            l2 = obj.mesh.indexer.multiToLinear(m, 0);

            T = zeros(np^2*ne, 3);

            F = G.fluxData(2, 2*d).matrix;
            [i, j, v] = find(P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l1(:).' - 1)*np, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l2(:).' - 1)*np, j(:)), [], 1);
            T(k, 3) = kron(c(l1), v(:));
        end

        function T = assembleOuterLeftFluxBoundaryTerms(obj, d, c)
            % ASSEMBLEOUTERLEFTFLUXBOUNDARYTERMS Assemble outer left flux
            % boundary terms.
            %
            %   T = assembleOuterLeftFluxBoundaryTerms(obj, d, c) computes
            %   the outer contribution to a left flux along the specified
            %   dimension on the boundary, dispatching to appropriate
            %   boundary condition methods.
            %
            % Inputs:
            %   obj - The EulerianGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector, one per element
            %
            % Outputs:
            %   T - Triplet matrix for boundary contribution

            switch obj.bcType
                case 0
                    T = obj.assembleOuterLeftFluxPeriodicBoundaryTerms(d, c);
                case 1
                    T = obj.assembleOuterLeftFluxDirichletBoundaryTerms(d, c);
            end
        end

        function T = assembleOuterRightFluxBoundaryTerms(obj, d, c)
            % ASSEMBLEOUTERRIGHTFLUXBOUNDARYTERMS Assemble outer right flux
            % boundary terms.
            %
            %   T = assembleOuterRightFluxBoundaryTerms(obj, d, c) computes
            %   the outer contribution to a right flux along the specified
            %   dimension on the boundary, dispatching to appropriate
            %   boundary condition methods.
            %
            % Inputs:
            %   obj - The EulerianGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector, one per element
            %
            % Outputs:
            %   T - Triplet matrix for boundary contribution

            switch obj.bcType
                case 0
                    T = obj.assembleOuterRightFluxPeriodicBoundaryTerms(d, c);
                case 1
                    T = obj.assembleOuterRightFluxDirichletBoundaryTerms(d, c);
            end
        end
        
        function T = assembleInnerLeftFluxPeriodicBoundaryTerms(obj, d, c)
            % ASSEMBLEINNERLEFTFLUXPERIODICBOUNDARYTERMS Assemble inner
            % left flux periodic boundary terms.
            %
            %   T = assembleInnerLeftFluxPeriodicBoundaryTerms(obj, d, c)
            %   computes the inner contribution to a left flux along the
            %   specified dimension on boundaries with periodic boundary
            %   conditions.
            %
            % Inputs:
            %   obj - The EulerianGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector, one per element
            %
            % Outputs:
            %   T - Triplet matrix for periodic boundary contribution

            np = obj.nLocalDofs;
            P = obj.fe.projector;
            G = obj.op.gradient;

            m = obj.mesh.allElementMultiIndices;
            e = m(:, d) == 1 | m(:, d) == obj.mesh.resolution(d);
            ne = numel(find(e));
            m = m(e, :);
            l = obj.mesh.indexer.multiToLinear(m, 0);

            T = zeros(np^2*ne, 3);
            
            F = G.fluxData(1, 2*d).matrix;
            [i, j, v] = find(P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l(:).' - 1)*np, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l(:).' - 1)*np, j(:)), [], 1);
            T(k, 3) = kron(c(l), v(:));
        end

        function T = assembleInnerRightFluxPeriodicBoundaryTerms(obj, d, c)
            % ASSEMBLEINNERRIGHTFLUXPERIODICBOUNDARYTERMS Assemble inner
            % right flux periodic boundary terms.
            %
            %   T = assembleInnerRightFluxPeriodicBoundaryTerms(obj, d, c)
            %   computes the inner contribution to a right flux along the
            %   specified dimension on boundaries with periodic boundary
            %   conditions.
            %
            % Inputs:
            %   obj - The EulerianGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector, one per element
            %
            % Outputs:
            %   T - Triplet matrix for periodic boundary contribution

            np = obj.nLocalDofs;
            P = obj.fe.projector;
            G = obj.op.gradient;

            m = obj.mesh.allElementMultiIndices;
            e = m(:, d) == 1 | m(:, d) == obj.mesh.resolution(d);
            ne = numel(find(e));
            m = m(e, :);
            l = obj.mesh.indexer.multiToLinear(m, 0);

            T = zeros(np^2*ne, 3);
            
            F = G.fluxData(1, 2*d-1).matrix;
            [i, j, v] = find(-P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l(:).' - 1)*np, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l(:).' - 1)*np, j(:)), [], 1);
            T(k, 3) = kron(c(l), v(:));
        end
    
        function T = assembleOuterLeftFluxPeriodicBoundaryTerms(obj, d, c)
            % ASSEMBLEOUTERLEFTFLUXPERIODICBOUNDARYTERMS Assemble outer
            % left flux periodic boundary terms.
            %
            %   T = assembleOuterLeftFluxPeriodicBoundaryTerms(obj, d, c)
            %   computes the outer contribution to a left flux along the
            %   specified dimension on boundaries with periodic boundary
            %   conditions.
            %
            % Inputs:
            %   obj - The EulerianGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector, one per element
            %
            % Outputs:
            %   T - Triplet matrix for periodic boundary contribution

            np = obj.nLocalDofs;
            P = obj.fe.projector;
            G = obj.op.gradient;

            m = obj.mesh.allElementMultiIndices;
            e = m(:, d) == 1 | m(:, d) == obj.mesh.resolution(d);
            ne = numel(find(e));
            m = m(e, :);
            l1 = obj.mesh.indexer.multiToLinear(m, 0);
            m(:, d) = m(:, d) - 1;
            l2 = obj.mesh.indexer.multiToLinear(m, 0);

            T = zeros(np^2*ne, 3);
            
            F = G.fluxData(2, 2*d-1).matrix;
            [i, j, v] = find(-P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l1(:).' - 1)*np, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l2(:).' - 1)*np, j(:)), [], 1);
            T(k, 3) = kron(c(l1), v(:));
        end
        
        function T = assembleOuterRightFluxPeriodicBoundaryTerms(obj, d, c)
            % ASSEMBLEOUTERRIGHTFLUXPERIODICBOUNDARYTERMS Assemble outer
            % right flux periodic boundary terms.
            %
            %   T = assembleOuterRightFluxPeriodicBoundaryTerms(obj, d, c)
            %   computes the outer contribution to a right flux along the
            %   specified dimension on boundaries with periodic boundary
            %   conditions.
            %
            % Inputs:
            %   obj - The EulerianGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector, one per element
            %
            % Outputs:
            %   T - Triplet matrix for periodic boundary contribution

            np = obj.nLocalDofs;
            P = obj.fe.projector;
            G = obj.op.gradient;

            m = obj.mesh.allElementMultiIndices;
            e = m(:, d) == 1 | m(:, d) == obj.mesh.resolution(d);
            ne = numel(find(e));
            m = m(e, :);
            l1 = obj.mesh.indexer.multiToLinear(m, 0);
            m(:, d) = m(:, d) + 1;
            l2 = obj.mesh.indexer.multiToLinear(m, 0);

            T = zeros(np^2*ne, 3);
            
            F = G.fluxData(2, 2*d).matrix;
            [i, j, v] = find(P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l1(:).' - 1)*np, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l2(:).' - 1)*np, j(:)), [], 1);
            T(k, 3) = kron(c(l1), v(:));
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