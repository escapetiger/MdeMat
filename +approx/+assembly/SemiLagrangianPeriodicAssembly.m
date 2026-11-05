classdef SemiLagrangianPeriodicAssembly < approx.assembly.Assembly
    % SEMILAGRANGIANPERIODICASSEMBLY Assembly for semi-Lagrangian dynamic
    % with periodic boundary conditions.
    %
    %   SemiLagrangianPeriodicAssembly implements the assembly of transport
    %   operators for semi-Lagrangian methods with periodic boundary
    %   conditions. It handles characteristic tracing with automatic
    %   wrapping at domain boundaries, enabling stable advection
    %   computations without explicit boundary condition enforcement.
    %
    %   The periodic treatment allows characteristics that exit one side of
    %   the domain to re-enter from the opposite side, maintaining
    %   conservation and stability properties. This is particularly useful
    %   for problems on periodic domains or with circular advection
    %   patterns.
    %
    % Notes:
    %   Unlike Dirichlet assembly, periodic assembly only requires volume
    %   assembly as boundary conditions are handled implicitly through
    %   the periodic indexing scheme.
    %
    % See also:
    %   approx.assembly.Assembly,
    %   approx.assembly.SemiLagrangianDirichletAssembly,
    %   approx.mesh.UniformGrid

    properties (Constant)
        BC_TYPE = 'periodic' % Boundary condition type identifier
    end

    methods
        function T = assembleVolume(obj)
            % ASSEMBLEVOLUME Assemble volume terms for periodic
            % semi-Lagrangian transport.
            %
            %   T = assembleVolume(obj) assembles the triplet
            %   representation of the transport operator matrix for all
            %   elements with periodic boundary conditions. The method
            %   handles characteristic tracing with automatic wrapping at
            %   domain boundaries using periodic indexing.
            %
            % Inputs:
            %   obj - The SemiLagrangianPeriodicAssembly object
            %
            % Outputs:
            %   T - Triplet matrix representation [row, col, val]
            %
            % Examples:
            %   % Standard volume assembly for periodic transport
            %   triplets = assembly.assembleVolume();
            %   n = assembly.space.nGlobalDofs;
            %   A = core.linalg.sparseFromTriplets(triplets, n, n);
            %   
            %   % Matrix A now includes periodic wrapping automatically
            %   U_new = A * U_old;  % Time step with periodic BCs
            %
            % Notes:
            %   Periodic indexing (flag=0 in multiToLinear) automatically
            %   handles wrap-around at domain boundaries. All elements
            %   participate in assembly unlike Dirichlet case.

            nl = obj.space.nLocalDofs;
            np = obj.operator.clipper.nVolumePieces;
            ne = obj.space.mesh.nTotalElements;

            %< Project mass matrices for each volume piece
            fM = @(i) obj.space.element.projector.project(obj.operator.mass(i).matrix);
            M = arrayfun(fM, 1:np, 'Un', 0);
            s = obj.operator.clipper.volumeShifts;
            
            T = cell(np, 1);
            for p = 1:np
                T{p} = zeros(nl^2*ne, 3);
                [i, j, v] = find(M{p});
                nv = numel(v);
                
                %< Apply volume shifts with periodic wrapping
                m = obj.space.mesh.allElementMultiIndices;
                m = m + s(:, p).';
                
                %< Convert to linear indices with periodic boundary conditions
                l = obj.space.mesh.indexer.multiToLinear(m, 0);  % flag=0 for periodic
                k = 1:(nv * ne);
                
                %< Assemble triplet entries
                T{p}(k, 1) = reshape(bsxfun(@plus, (0:ne - 1)*nl, i(:)), [], 1);
                T{p}(k, 2) = reshape(bsxfun(@plus, (l(:).' - 1)*nl, j(:)), [], 1);
                T{p}(k, 3) = kron(ones(ne, 1), v(:));
            end
            T = vertcat(T{:});
        end
    end
end