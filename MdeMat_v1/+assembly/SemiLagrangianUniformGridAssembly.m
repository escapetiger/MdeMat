classdef SemiLagrangianUniformGridAssembly < fem.assembly.GridAssembly
    % SEMILAGRANGIANUNIFORMGRIDASSEMBLY Grid-based assembly for
    % semi-Lagrangian methods.
    %
    %   SemiLagrangianUniformGridAssembly provides assembly methods for
    %   semi-Lagrangian finite element operators on uniform structured
    %   grids. This specialized assembly handles transport operators and
    %   boundary conditions for semi-Lagrangian discretizations.
    %
    %   The assembly supports both free transport operators and boundary
    %   condition assembly, with particular emphasis on handling
    %   characteristic tracing and interpolation required by
    %   semi-Lagrangian methods.
    %
    % Examples:
    %   % Create semi-Lagrangian assembly with periodic boundaries
    %   assembly = SemiLagrangianUniformGridAssembly(fe, mesh, op, 0);
    %   
    %   % Assemble transport operator only
    %   A = assembly.transport();
    %   
    %   % Assemble both transport and boundary conditions
    %   [b, A] = assembly.transport(t0, t1, boundaryFunction);
    %
    % See also:
    %   fem.assembly.GridAssembly, fem.assembly.Assembly

    methods
        function obj = SemiLagrangianUniformGridAssembly(fe, mesh, op, bcType)
            % SEMILAGRANGIANUNIFORMGRIDASSEMBLY Constructor for
            % semi-Lagrangian assembly.
            %
            %   obj = SemiLagrangianUniformGridAssembly(fe, mesh, op,
            %   bcType) creates an assembly for semi-Lagrangian finite
            %   element operators on uniform structured grids.
            %
            % Inputs:
            %   fe - Finite element object
            %   mesh - Uniform grid mesh object
            %   op - Semi-Lagrangian finite element operator object
            %   bcType - Boundary condition type (0=periodic, 1=Dirichlet)
            %
            % Outputs:
            %   obj - Constructed SemiLagrangianUniformGridAssembly object

            core.except.assert( ...
                isa(op, 'fem.element.SemiLagrangianFiniteElementOperator'), ...
                'InvalidInput', ...
                'Require a Semi-Lagrangian finite element operator.')

            obj@fem.assembly.GridAssembly(fe, mesh, op, bcType);
        end

        function [b, A] = transport(obj, varargin)
            % TRANSPORT Assemble semi-Lagrangian transport operator.
            %
            %   b = transport(obj) assembles boundary condition vector for
            %   the transport operator with current boundary conditions.
            %
            %   [b, A] = transport(obj) assembles both the boundary condition
            %   vector and the free transport operator matrix.
            %
            %   [b, A] = transport(obj, tBegin, tEnd, f) assembles transport
            %   components with specified time interval and boundary function.
            %
            % Inputs:
            %   obj - The SemiLagrangianUniformGridAssembly object
            %   varargin - Input arguments
            %<  tBegin - Beginning time for boundary conditions (optional)
            %<  tEnd - End time for boundary conditions (optional)
            %<  f - Boundary condition function handle (optional)
            %
            % Outputs:
            %   b - Boundary condition vector (sparse column vector)
            %   A - Free transport operator matrix (sparse matrix)
            
            %< ------------------------------------------------------------
            %< COMMON PARAMETERS
            %< ------------------------------------------------------------
            n = obj.nGlobalDofs;

            %< ------------------------------------------------------------
            %< FREE TRANSPORT
            %< ------------------------------------------------------------
            if nargout >= 2
                T = obj.assembleFreeTransport();
                A = obj.sparseFromTriplets(T, n, n);
            end

            %< ------------------------------------------------------------
            %< BOUNDARY CONDITION
            %< ------------------------------------------------------------
            if nargout >= 1
                b = [];
                if obj.bcType == 1
                    [tBegin, tEnd, f] = varargin{:};
                    if abs(tBegin - tEnd) < 1e-8, return; end
                    T = obj.assembleBoundaryCondition(tBegin, tEnd, f);
                    b = obj.sparseFromTriplets(T, n, 1);
                end
            end
        end
    end

    methods (Access = protected)
        function T = assembleFreeTransport(obj)
            % ASSEMBLEFREETRANSPORT Assemble free transport operator.
            %
            %   T = assembleFreeTransport(obj) assembles the triplet
            %   representation of the free transport operator matrix,
            %   handling characteristic tracing and interpolation for
            %   semi-Lagrangian methods.
            %
            % Inputs:
            %   obj - The SemiLagrangianUniformGridAssembly object
            %
            % Outputs:
            %   T - Triplet matrix [row, col, value] for transport operator

            nl = obj.nLocalDofs;
            np = obj.op.clipper.nVolumePieces;
            ne = obj.mesh.nTotalElements;

            M = arrayfun(@(i)obj.fe.projector.project(obj.op.mass(i).matrix), 1:np, 'Un', 0);
            s = obj.op.clipper.volumeShifts;
            d = obj.nDims;

            T = cell(np, 1);
            for p = 1:np
                T{p} = zeros(nl^2*ne, 3);
                [i, j, v] = find(M{p});
                nv = numel(v);
                m = obj.mesh.allElementMultiIndices;
                if obj.bcType == 0
                    m = m + s(:, p).';
                    l = obj.mesh.indexer.multiToLinear(m, 0);
                    k = 1:(nv * ne);
                    T{p}(k, 1) = reshape(bsxfun(@plus, (0:ne - 1)*nl, i(:)), [], 1);
                    T{p}(k, 2) = reshape(bsxfun(@plus, (l(:).' - 1)*nl, j(:)), [], 1);
                    T{p}(k, 3) = kron(ones(ne, 1), v(:));
                else
                    b = m + s(:, p).';
                    b = (b < 1) .* (2 * (1:d) - 1) + ...
                        (b > obj.mesh.resolution) .* (2 * (1:d));
                    m = m(sum(b ~= 0, 2) == 0, :);
                    nb = size(m, 1);
                    k = 1:(nv * nb);
                    l1 = obj.mesh.indexer.multiToLinear(m, 1);
                    l2 = obj.mesh.indexer.multiToLinear(m + s(:, p).', 1);
                    T{p}(k, 1) = reshape(bsxfun(@plus, (l1(:).' - 1)*nl, i(:)), [], 1);
                    T{p}(k, 2) = reshape(bsxfun(@plus, (l2(:).' - 1)*nl, j(:)), [], 1);
                    T{p}(k, 3) = kron(ones(nb, 1), v(:));
                end
            end
            T = vertcat(T{:});
        end
    
        function T = assembleBoundaryCondition(obj, tBegin, tEnd, f)
            % ASSEMBLEBOUNDARYCONDITION Assemble boundary condition terms.
            %
            %   T = assembleBoundaryCondition(obj, tBegin, tEnd, f)
            %   assembles the boundary condition contribution for
            %   semi-Lagrangian transport, handling characteristic tracing
            %   to domain boundaries and evaluation of boundary data.
            %
            % Inputs:
            %   obj - The SemiLagrangianUniformGridAssembly object
            %   tBegin - Beginning time for integration
            %   tEnd - End time for integration  
            %   f - Boundary condition function handle f(x, t)
            %
            % Outputs:
            %   T - Triplet matrix [row, col, value] for boundary conditions

            nl = obj.nLocalDofs;
            np = obj.op.clipper.nVolumePieces;
            nd = obj.nDims;
            ht = tBegin - tEnd;
            v = obj.op.clipper.velocity;

            T = cell(1, np);
            for p = 1:np
                %< Find all elements whose pth volume piece goes outside
                %< the computational region
                m = obj.mesh.allElementMultiIndices;
                s = obj.op.clipper.volumeShifts(:, p).';
                b = (m + s < 1) .* (2 * (1:nd) - 1) + (m + s > obj.mesh.resolution) .* (2 * (1:nd));
                b = m(sum(b ~= 0, 2) > 0, :);
                nb = size(b, 1);
                k = obj.mesh.indexer.multiToLinear(b, 1);
                T{p} = zeros(nb*nl, 3);
                T{p}(:, 1) = reshape(bsxfun(@plus, (k(:)-1).'*nl, (1:nl).'), [], 1);
                T{p}(:, 2) = 1;

                if nb == 0, continue; end

                %< Find all quadrature points to evaluate
                nq = obj.op.clipper.volumeEulerianData(1, p).nPoints;
                X1 = zeros(nd + 1, nb*nq);
                X2 = zeros(nd + 1, nb*nq);
                Xq = obj.op.clipper.volumeEulerianData(1, p).nodes;
                X2(1:nd, :) = obj.mesh.collocate(Xq, b);
                X2(nd+1, :) = tBegin;
                X1(1:nd, :) = X2(1:nd, :) - ht * v(:);
                X1(nd+1, :) = tEnd;

                %< Find intersection of characteristics and boundaries
                I = core.geometry.LineSegmentByPlaneIntersector(nd+1);
                b = reshape(obj.mesh.bbox, 2, []);
                b = -b(2*(1:nd) - (v > 0));
                Xb = zeros(nd + 1, nb*nq);
                for d = 1:nd
                    a = zeros(nd+1, 1);
                    a(d) = 1;
                    [X, TF] = I.intersect(X1, X2, a, b(d));
                    Xb(:, TF) = X;
                end

                %< Evaluate boundary function and project
                U = f(Xb(1:nd, :), Xb(nd+1, :));
                U = reshape(U, nq, nb);
                P = obj.op.volumePieceData(1, p).vector;
                T{p}(:, 3) = reshape(obj.fe.projector.project(P * U), [], 1);
            end
            T = vertcat(T{:});
        end
    end
end