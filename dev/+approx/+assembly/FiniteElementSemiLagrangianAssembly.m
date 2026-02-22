classdef FiniteElementSemiLagrangianAssembly < approx.assembly.FiniteElementAssembly
    % FINITEELEMENTSEMILAGRANGIANASSEMBLY Assembly for semi-Lagrangian finite element methods.
    %
    %   FiniteElementSemiLagrangianAssembly implements assembly operations for
    %   semi-Lagrangian finite element discretizations of transport equations.
    %   The semi-Lagrangian method follows characteristic curves backward in time
    %   to solve advection-dominated problems with excellent stability properties.
    %
    %   The assembly process involves clipping volume elements based on the
    %   advection field and time step to handle characteristics that may exit
    %   the computational domain. Status flags track whether clipping has been
    %   performed for matrix and vector assembly operations.
    %
    %   This implementation is specialized for L2 semi-Lagrangian orthotope
    %   elements on uniform grids, providing optimal efficiency for structured
    %   transport problems.
    %
    %   The assembly requires BH1SLOrthotopeElement and UniformGrid for proper
    %   functioning. Status flags use bitwise operations to track clipping state.
    %
    % See also:
    %   FiniteElementAssembly, BH1SLOrthotopeElement, UniformGrid
    
    properties (Constant)
        ClipMatrixMask = 0b10 % Clipped when assembling matrix
        ClipVectorMask = 0b01 % Clipped when assembling vector
    end
    
    properties
        BcType % Boundary condition type
        Status % Status mask
    end
    
    methods
        function obj = FiniteElementSemiLagrangianAssembly(space, bcType)
            % FINITEELEMENTSEMILAGRANGIANASSEMBLY Constructor for semi-Lagrangian assembly.
            %
            %   obj = FiniteElementSemiLagrangianAssembly(space, bcType) creates
            %   a semi-Lagrangian assembly object with specified finite element
            %   @a space and boundary condition @a bcType. Validates that the
            %   space uses compatible element and mesh types, then initializes
            %   the assembly status.
            
            arguments
                space {mustBeL2SLUniformGridSpace}
                bcType{mustBeMember(bcType, {'periodic', 'dirichlet'})}
            end
            
            obj@approx.assembly.FiniteElementAssembly(space);
            obj.BcType = bcType;
            obj.reset();
        end
        
        function obj = setBcType(obj, bcType)
            % SETBCTYPE Set boundary condition type for semi-Lagrangian assembly.
            %
            %   obj = setBcType(obj, bcType) sets the boundary condition type
            %   to @a bcType and validates it. Semi-Lagrangian methods support
            %   both periodic and Dirichlet boundary conditions with different
            %   handling of characteristics that exit the domain.
            
            arguments
                obj approx.assembly.FiniteElementSemiLagrangianAssembly
                bcType {mustBeMember(bcType, {'periodic', 'dirichlet'})}
            end
            
            obj.BcType = bcType;
        end
        
        function obj = reset(obj)
            % RESET Reset assembly status flags.
            %
            %   obj = reset(obj) clears all status flags to indicate that no
            %   clipping operations have been performed. This method should be
            %   called when starting a new time step or when the advection field
            %   changes to ensure proper clipping state management.
            
            arguments
                obj approx.assembly.FiniteElementSemiLagrangianAssembly
            end
            
            obj.Status = 0;
        end
        
        function G = assembleMatrix(obj, field, tBegin, tEnd)
            % ASSEMBLEMATRIX Assemble semi-Lagrangian transport matrix.
            %
            %   G = assembleMatrix(obj, field, tBegin, tEnd) assembles the mass
            %   matrix for semi-Lagrangian transport with advection @a field
            %   from time @a tBegin to @a tEnd. The method clips volume elements
            %   based on characteristic trajectories and assembles contributions
            %   from Eulerian and upstream integration points.
            %
            %   The time step \f$ \Delta t = tBegin - tEnd \f$ and spatial
            %   discretization determine the clipping parameters. Status flags
            %   are updated to track clipping operations.
            
            arguments
                obj approx.assembly.FiniteElementSemiLagrangianAssembly
                field {mustBeNumeric}
                tBegin {mustBeNumeric}
                tEnd {mustBeNumeric}
            end
            
            nd = obj.Space.NDims;
            
            % Parse advection field
            field = reshape(field, nd, []);
            
            % Clipping
            ht = tBegin - tEnd;
            hx = cellfun(@(h) h(1), obj.Space.Mesh.Spacings);
            obj.Space.Element.clipVolume(field, ht, hx);
            
            % Update status
            obj.Status = bitor(obj.Status, obj.ClipMatrixMask);
            
            % Assemble global mass matrix
            switch obj.BcType
                case 'periodic'
                    G = obj.assembleMatrixPeriodic();
                case 'dirichlet'
                    G = obj.assembleMatrixDirichlet();
            end
        end
        
        function G = assembleVector(obj, field, tBegin, tEnd, func)
            % ASSEMBLEVECTOR Assemble semi-Lagrangian boundary vector.
            %
            %   G = assembleVector(obj, field, tBegin, tEnd, func) assembles
            %   the boundary contribution vector for semi-Lagrangian transport
            %   with Dirichlet boundary conditions. Handles characteristics that
            %   exit the computational domain by evaluating boundary function
            %   @a func at intersection points.
            %
            %   For periodic boundary conditions, no vector assembly is needed.
            %   The method ensures clipping consistency with matrix assembly
            %   by performing clipping if not already done.
            
            arguments
                obj approx.assembly.FiniteElementSemiLagrangianAssembly
                field {mustBeNumeric}
                tBegin {mustBeNumeric}
                tEnd {mustBeNumeric}
                func {mustBeA(func, 'function_handle')}
            end
            
            if bitand(obj.Status, obj.ClipMatrixMask)
                nd = obj.Space.NDims;
                
                % Parse advection field
                field = reshape(field, nd, []);
                
                % Clipping
                ht = tBegin - tEnd;
                hx = cellfun(@(h) h(1), obj.Space.Mesh.Spacings);
                obj.Space.Element.clipVolume(field, ht, hx);
                
                % Update status
                obj.Status = bitor(obj.Status, obj.ClipMatrixMask);
            end
            
            % Skip for periodic boundary condition
            if strcmpi(obj.BcType, 'periodic')
                G = obj.createZeros(nCols=1);
                return;
            end
            
            G = obj.assembleVectorDirichlet(field, tBegin, tEnd, func);
        end
    end
    
    methods (Access = private)
        function G = assembleMatrixPeriodic(obj)
            % ASSEMBLEMATRIXPERIODIC Assemble matrix for periodic boundary conditions.
            
            G = obj.createZeros();
            for i = 1:obj.Space.Element.NVolumePieces
                %< Compute local mass matrix
                W = diag(obj.Space.Element.VolumePieces.eulerian(i).Weights); % (nq, nq)
                V = obj.Space.Element.VolumePieces.eulerian(i).Values; % (nl, nq)
                U = obj.Space.Element.VolumePieces.upstream(i).Values; % (nl, nq)
                L = V * W * U.'; % nl x nl
                
                %< Add contribution to global mass matrix
                IE = 1:obj.Space.NMeshElements;
                ME = obj.Space.Mesh.Indexer.linearToMulti(IE);
                ME = ME + obj.Space.Element.VolumePieces.indexShifts(:, i).';
                JE = obj.Space.Mesh.Indexer.multiToLinear(ME, pad='wrap');
                GP = obj.createMatrix(L, TEI=IE, SEI=JE);
                G = G + GP;
            end
        end
        
        function G = assembleMatrixDirichlet(obj)
            % ASSEMBLEMATRIXDIRICHLET Assemble matrix for Dirichlet boundary conditions.
            
            d = obj.Space.NDims;
            G = obj.createZeros();
            for i = 1:obj.Space.Element.NVolumePieces
                %< Compute local mass matrix
                W = diag(obj.Space.Element.VolumePieces.eulerian(i).Weights); % (nq, nq)
                V = obj.Space.Element.VolumePieces.eulerian(i).Values; % (nl, nq)
                U = obj.Space.Element.VolumePieces.upstream(i).Values; % (nl, nq)
                L = V * W * U.'; % (nl, nl)
                
                %< Find elements that have at least one interior upstream root
                ME = obj.Space.Mesh.Indexer.generate();
                SE = obj.Space.Element.VolumePieces.indexShifts(:, i).';
                MU = ME + SE;
                mask1 = (MU < 1) .* (2 * (1:d) - 1);
                mask2 = (MU > obj.Space.Mesh.Resolution) .* (2 * (1:d));
                mask = mask1 + mask2;
                mask = sum(mask ~= 0, 2) == 0;
                ME = ME(mask, :);
                
                %< Add contribution to global mass matrix
                IE = obj.Space.Mesh.Indexer.multiToLinear(ME, pad='empty');
                JE = obj.Space.Mesh.Indexer.multiToLinear(ME+SE, pad='empty');
                GP = obj.createMatrix(L, TEI=IE, SEI=JE);
                G = G + GP;
            end
        end
        
        function G = assembleVectorDirichlet(obj, field, tBegin, tEnd, func)
            % ASSEMBLEVECTORDIRICHLET Assemble Dirichlet boundary vector.
            
            np = obj.Space.Element.NVolumePieces;
            nd = obj.Space.NDims;
            ht = tBegin - tEnd;
            
            G = obj.createZeros(nCols=1);
            for i = 1:np
                %< Find elements whose pth volume piece exits domain
                ME = obj.Space.Mesh.Indexer.generate();
                SE = obj.Space.Element.VolumePieces.indexShifts(:, i).';
                MU = ME + SE;
                mask1 = (MU < 1) .* (2 * (1:nd) - 1);
                mask2 = (MU > obj.Space.Mesh.Resolution) .* (2 * (1:nd));
                mask = mask1 + mask2;
                mask = sum(mask ~= 0, 2) > 0;
                ME = ME(mask, :);
                IE = obj.Space.Mesh.Indexer.multiToLinear(ME, pad='empty');
                ne = size(ME, 1);
                
                if ne == 0, return; end
                
                %< Set up characteristic end points
                nq = obj.Space.Element.VolumePieces.eulerian(i).NPoints;
                X1 = zeros(nd+1, ne*nq);
                X2 = zeros(nd+1, ne*nq);
                Xq = obj.Space.Element.VolumePieces.eulerian(i).Nodes;
                X2(1:nd, :) = obj.Space.Mesh.collocate(Xq, IE);
                X2(nd+1, :) = tBegin;
                X1(1:nd, :) = X2(1:nd, :) - ht * field(:);
                X1(nd+1, :) = tEnd;
                
                %< Find boundary intersections of characteristics
                intersector = core.geometry.SegmentByPlaneIntersector(nd+1);
                bbox = reshape(obj.Space.Mesh.BBox, 2, []);
                bbox = -bbox(2*(1:nd)-(field > 0));
                Xb = zeros(nd+1, nq*ne);
                for d = 1:nd
                    a = zeros(nd+1, 1);
                    a(d) = 1;
                    [X, TF] = intersector.intersect(X1, X2, a, bbox(d));
                    Xb(:, TF) = X;
                end
                
                %< Evaluate boundary function and project to element space
                V = obj.Space.Element.VolumePieces.eulerian(i).Values; % (nl, nq)
                W = diag(obj.Space.Element.VolumePieces.eulerian(i).Weights); % (nq, nq)
                U = func(Xb(1:nd, :), t=Xb(nd+1, :)); % (1, nq*ne)
                U = reshape(U, nq, ne);
                L = V * W * U; % (nl, ne)
                GP = obj.createVector(L, TEI=IE);
                G = G + GP;
            end
        end
    end
end

function mustBeL2SLUniformGridSpace(x)
if ~isa(x, 'approx.space.FiniteElementSpace')
    error('Input must be a FiniteElementSpace object.');
end
tf1 = isa(x.Element, 'approx.element.BH1SLOrthotopeElement');
tf2 = isa(x.Mesh, 'approx.mesh.UniformGrid');
if ~(tf1 && tf2)
    error('SL assembly requires L2 SL orthotope element and uniform grid.');
end
end
