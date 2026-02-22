classdef SemiLagrangianDirichletAssembly < approx.assembly.Assembly
    % SEMILAGRANGIANDIRICHLETASSEMBLY Assembly for semi-Lagrangian dynamic
    % with Dirichlet boundary conditions.
    %
    %   SemiLagrangianDirichletAssembly implements the assembly of
    %   transport operators for semi-Lagrangian methods with Dirichlet
    %   boundary conditions. It handles characteristic tracing, boundary
    %   intersections, and the incorporation of boundary data into the
    %   discrete system.
    %
    %   The class manages two key operations: volume assembly for interior
    %   characteristic tracing and boundary assembly for enforcing
    %   Dirichlet conditions. Characteristics that exit the computational
    %   domain are traced to boundary intersections where boundary data is
    %   evaluated and incorporated into the system.
    %
    % Notes:
    %   Boundary condition assembly involves geometric intersection
    %   computations and requires careful handling of characteristic
    %   trajectories that exit the domain.
    %
    % See also:
    %   approx.assembly.Assembly,
    %   approx.assembly.SemiLagrangianPeriodicAssembly,
    %   core.geometry.LineSegmentByPlaneIntersector

    properties (Constant)
        BC_TYPE = 'dirichlet' % Boundary condition type identifier
    end

    methods
        function T = assembleVolume(obj)
            % ASSEMBLEVOLUME Assemble volume terms for semi-Lagrangian
            % transport.
            %
            %   T = assembleVolume(obj) assembles the triplet
            %   representation of the transport operator matrix for
            %   interior elements. The method handles characteristic
            %   tracing and interpolation for semi-Lagrangian methods,
            %   excluding elements that intersect domain boundaries.
            %
            % Inputs:
            %   obj - The SemiLagrangianDirichletAssembly object
            %
            % Outputs:
            %   T - Triplet matrix representation [row, col, val]

            nl = obj.space.nLocalDofs;
            np = obj.operator.clipper.nVolumePieces;
            ne = obj.space.mesh.nTotalElements;

            %< Project mass matrices for each volume piece
            fM = @(i) obj.space.element.projector.project(obj.operator.mass(i).matrix);
            M = arrayfun(fM, 1:np, 'Un', 0);
            s = obj.operator.clipper.volumeShifts;
            d = obj.space.nDims;
            
            T = cell(np, 1);
            for p = 1:np
                T{p} = zeros(nl^2*ne, 3);
                [i, j, v] = find(M{p});
                nv = numel(v);
                m = obj.space.mesh.allElementMultiIndices;
                b = m + s(:, p).';
                
                %< Identify boundary interactions
                b = (b < 1) .* (2 * (1:d) - 1) + ...
                    (b > obj.space.mesh.resolution) .* (2 * (1:d));
                    
                %< Keep only interior elements
                m = m(sum(b ~= 0, 2) == 0, :);
                nb = size(m, 1);
                k = 1:(nv * nb);
                
                %< Compute global indices for assembly
                l1 = obj.space.mesh.indexer.multiToLinear(m, 1);
                l2 = obj.space.mesh.indexer.multiToLinear(m + s(:, p).', 1);
                T{p}(k, 1) = reshape(bsxfun(@plus, (l1(:).' - 1)*nl, i(:)), [], 1);
                T{p}(k, 2) = reshape(bsxfun(@plus, (l2(:).' - 1)*nl, j(:)), [], 1);
                T{p}(k, 3) = kron(ones(nb, 1), v(:));
            end
            T = vertcat(T{:});
        end
    
        function T = assembleBoundaryCondition(obj, tBegin, tEnd, f)
            % ASSEMBLEBOUNDARYCONDITION Assemble Dirichlet boundary
            % condition terms.
            %
            %   T = assembleBoundaryCondition(obj, tBegin, tEnd, f)
            %   assembles the boundary condition contribution for
            %   semi-Lagrangian transport. The method traces
            %   characteristics that exit the domain to their boundary
            %   intersections and evaluates the boundary condition function
            %   at those points.
            %
            % Inputs:
            %   obj - The SemiLagrangianDirichletAssembly object
            %   tBegin - Beginning time for characteristic integration (scalar)
            %   tEnd - End time for characteristic integration (scalar)
            %   f - Boundary condition function handle f(x, t)
            %
            % Outputs:
            %   T - Triplet matrix representation [row, col, val]
            %
            % Notes:
            %   Uses geometric intersection algorithms to find where
            %   characteristics meet domain boundaries. Handles
            %   multi-dimensional boundary intersections correctly.

            nl = obj.space.nLocalDofs;
            np = obj.operator.clipper.nVolumePieces;
            nd = obj.space.nDims;
            ht = tBegin - tEnd;
            v = obj.operator.clipper.velocity;

            T = cell(1, np);
            for p = 1:np
                %< Find elements whose pth volume piece exits domain
                m = obj.space.mesh.allElementMultiIndices;
                s = obj.operator.clipper.volumeShifts(:, p).';
                b = (m + s < 1) .* (2 * (1:nd) - 1) + ...
                    (m + s > obj.space.mesh.resolution) .* (2 * (1:nd));
                b = m(sum(b ~= 0, 2) > 0, :);
                nb = size(b, 1);
                
                %< Initialize triplet structure
                k = obj.space.mesh.indexer.multiToLinear(b, 1);
                T{p} = zeros(nb*nl, 3);
                T{p}(:, 1) = reshape(bsxfun(@plus, (k(:)-1).'*nl, (1:nl).'), [], 1);
                T{p}(:, 2) = 1;

                if nb == 0, continue; end

                %< Set up characteristic end points
                nq = obj.operator.clipper.volumeEulerianData(1, p).nPoints;
                X1 = zeros(nd + 1, nb*nq);
                X2 = zeros(nd + 1, nb*nq);
                Xq = obj.operator.clipper.volumeEulerianData(1, p).nodes;
                X2(1:nd, :) = obj.space.mesh.collocate(Xq, b);
                X2(nd+1, :) = tBegin;
                X1(1:nd, :) = X2(1:nd, :) - ht * v(:);
                X1(nd+1, :) = tEnd;

                %< Find boundary intersections of characteristics
                I = core.geometry.LineSegmentByPlaneIntersector(nd+1);
                bbox = reshape(obj.space.mesh.bbox, 2, []);
                bbox = -bbox(2*(1:nd) - (v > 0));
                Xb = zeros(nd + 1, nb*nq);
                
                for d = 1:nd
                    a = zeros(nd+1, 1);
                    a(d) = 1;
                    [X, TF] = I.intersect(X1, X2, a, bbox(d));
                    Xb(:, TF) = X;
                end

                %< Evaluate boundary function and project to element space
                U = f(Xb(1:nd, :), Xb(nd+1, :));
                U = reshape(U, nq, nb);
                P = obj.operator.volumePieceData(1, p).vector;
                T{p}(:, 3) = reshape(obj.space.element.projector.project(P * U), [], 1);
            end
            T = vertcat(T{:});
        end
    end
end