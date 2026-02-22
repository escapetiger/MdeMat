classdef UpwindAssembly < approx.assembly.FiniteElementAssembly
    % FINITEELEMENTUPWINDASSEMBLY Assembly for finite element upwind
    % operators.
    %
    %   UpwindAssembly handles the assembly of upwind finite
    %   element operators. These operators are mainly used for describing
    %   advection procedure.
    
    properties
        BcType % Boundary condition type
    end
    
    methods
        function obj = UpwindAssembly(space, bcType)
            % FINITEELEMENTUPWINDASSEMBLY Constructor for upwind assembly.
            %
            %   obj = UpwindAssembly(space, bcType) creates an
            %   upwind assembly object with specified finite element @a
            %   space and boundary condition @a bcType. The upwind method
            %   provides stable discretization of advection operators by
            %   biasing the flux calculation based on flow direction.
            
            arguments
                space approx.space.FiniteElementSpace
                bcType{mustBeMember(bcType, {'periodic', 'dirichlet'})}
            end
            
            obj@approx.assembly.FiniteElementAssembly(space);
            obj.setBcType(bcType);
        end
        
        function obj = setBcType(obj, bcType)
            % SETBCTYPE Set boundary condition type for upwind assembly.
            %
            %   obj = setBcType(obj, bcType) validates and sets the
            %   boundary condition type to @a bcType. Upwind methods handle
            %   boundary conditions by considering the flow direction at
            %   boundaries to determine whether inflow or outflow
            %   conditions apply.
            
            arguments
                obj approx.assembly.UpwindAssembly
                bcType{mustBeMember(bcType, {'periodic', 'dirichlet'})}
            end
            
            obj.BcType = bcType;
        end
        
        function G = assembleMatrix(obj, field, bias)
            % ASSEMBLEMATRIX Assemble upwind advection matrix.
            %
            %   G = assembleMatrix(obj, field) assembles the global matrix
            %   for the upwind discretization of the advection operator \f$
            %   \nabla \cdot (a u) \f$ where @a advection represents the
            %   advection velocity \f$ a \f$.
            %
            %   G = assembleMatrix(obj, field, bias) uses @a bias for
            %   upwind bias calculation, otherwise uses @a field.
            
            arguments
                obj approx.assembly.UpwindAssembly
                field{mustBeNumeric}
                bias{mustBeNumeric} = []
            end
            
            nd = obj.Space.NDims;
            
            field = reshape(field, nd, []);
            
            if isempty(bias)
                bias = field;
            else
                bias = reshape(bias, nd, []);
            end
            
            GV = obj.assembleVolumeMatrix(field);
            GF = obj.assembleFluxMatrix(field, bias);
            G = GF - GV;
        end
        
        function G = assembleVector(obj, field, bias, func, options)
            % ASSEMBLEVECTOR Assemble the vector contribution for
            % the divergence of piecewise constant advection field.
            %
            %   G = assembleVector(obj, field, bias, func) assembles the
            %   vector contribution for the upwind operator with advection
            %   @a field, upwind @a bias, and boundary function @a func.
            %
            %   G = assembleVector(obj, field, bias, func, options) allows
            %   additional arguments to be passed to function evaluation
            %   through @a options.args.
            
            arguments
                obj approx.assembly.UpwindAssembly
                field{mustBeNumeric}
                bias{mustBeNumeric}
                func{mustBeA(func, 'function_handle')}
                options.mode = ''
                options.args cell = {}
            end
            
            if strcmpi(obj.BcType, 'periodic')
                G = obj.createZeros(nCols=1);
                return;
            end
            
            nd = obj.Space.NDims;
            
            field = reshape(field, nd, []);
            
            if isempty(bias)
                bias = field;
            else
                bias = reshape(bias, nd, []);
            end
            
            G = obj.assembleFluxVector(field, bias, func, options.mode, options.args{:});
        end
    end
    
    methods (Access = private)
        function G = assembleVolumeMatrix(obj, field)
            % ASSEMBLEVOLUMEMATRIX Assemble the matrix contribution from
            % volume integrals for the divergence of piecewise constant
            % advection advection.
            
            %< Sizes
            nd = obj.Space.NDims;
            
            %< Inverse of Jacobian
            invJ = obj.Space.Mesh.computeElementInverseJacobians(); % (nd, nd, ne)
            
            %< Integration data
            W = diag(obj.Space.Element.Volume.Weights); % (nq, nq)
            V = obj.Space.Element.Volume.Partials{1}; % (nl, nq, nd)
            U = obj.Space.Element.Volume.Values; % (nl, nq)
            
            %< Reference matrices
            L = pagemtimes(pagemtimes(V, W), U.'); % (nl, nl, nd)
            L = obj.Space.Element.Approximator.project(L); % (nl, nl, nd)
            
            %< Physical marices
            L = reshape(L, [], nd) * reshape(invJ, nd, []); % (nl^2, nd, ne)
            L = reshape(L, size(V, 1), size(U, 1), nd, []);
            L = sum(L.*reshape(field, 1, 1, nd, []), 3);
            L = reshape(L, size(V, 1), size(U, 1), []);
            G = obj.createMatrix(L);
        end
        
        function G = assembleFluxMatrix(obj, field, bias)
            % ASSEMBLEFLUXMATRIX Assemble the matrix contribution from flux
            % integrals for the divergence of piecewise constant advection
            % advection.
            
            %< Initialization
            G = obj.createZeros();
            
            %< Add contribution from internal elements
            IEI = obj.Space.Mesh.getInternalElements(); % (:, 1)
            for LFI = 1:obj.Space.Element.NFluxes
                GF = obj.assembleInternalFluxMatrix(field, bias, IEI, LFI);
                G = G + GF;
            end
            
            %< Add contribution from boundary elements
            BEI = obj.Space.Mesh.getBoundaryElements(); % (:, 1)
            for LFI = 1:obj.Space.Element.NFluxes
                %< Find elements of which LFI-th face is not on boundary
                EI = obj.Space.Mesh.findInternalElements(BEI, LFI);
                if ~isempty(EI)
                    GF = obj.assembleInternalFluxMatrix(field, bias, EI, LFI);
                    G = G + GF;
                end
                
                %< Find elements of which LFI-th face is on boundary
                EI = obj.Space.Mesh.findBoundaryElements(BEI, LFI);
                if ~isempty(EI)
                    GF = obj.assembleInternalFluxMatrix(field, bias, EI, LFI);
                    G = G + GF;
                end
            end
        end
        
        function G = assembleInternalFluxMatrix(obj, field, bias, EI, LFI)
            % ASSEMBLEINTERNALFLUXMATRIX Assemble internal flux
            % contributions for upwind.
            
            % Sizes
            nl = obj.Space.NLocalDofs;
            
            %< Outward normal vectors
            normal = obj.Space.Mesh.computeOutwardNormals(EI, LFI); % (nd, :)
            
            %< Stabilizer
            eta = sum(bias.*normal, 1);
            mask = eta < 0;
            if size(bias, 2) == 1 && size(normal, 2) == 1
                mask = repmat(mask, 1, length(EI));
            end
            
            %< Test basis values and weights
            V = obj.Space.Element.Flux(LFI).Values; % (nl, nq)
            W = diag(obj.Space.Element.Flux(LFI).Weights); % (nq, nq)
            
            %< Face scaling factors
            ratio = obj.Space.Mesh.computeFaceScalingFactors(EI, LFI); % (1, :)
            ratio = reshape(ratio, 1, []);
            
            %< Face coefficients
            c = ratio .* sum(field.*normal, 1);
            
            %< Initialization
            G = obj.createZeros();
            
            %< Interior contribution
            if ~any(~mask)
                A = obj.createZeros();
            else
                U = obj.Space.Element.Flux(LFI).Values; % (nl, nq)
                L = V * W * U.'; % (nl, nl)
                L = obj.Space.Element.Approximator.project(L); % (nl, nl)
                L = reshape(L, nl, nl, 1);
                C = reshape(c(~mask), 1, 1, []);
                A = obj.createMatrix(L.*C, TEI=EI(~mask), SEI=EI(~mask));
            end
            G = G + A;
            
            %< Exterior contribution
            if ~any(mask)
                A = obj.createZeros();
            else
                [NEI, NLFI] = obj.Space.Mesh.findNeighborElements(EI(mask), LFI, obj.BcType); % (ne, 1), (ne, 1)
                U = arrayfun(@(j) obj.Space.Element.Flux(j).Values, NLFI, 'Un', 0);
                U = cat(3, U{:}); % (nl, nq, ne)
                L = pagemtimes(V*W, U.'); % (nl, nl, ne)
                L = obj.Space.Element.Approximator.project(L); % (nl, nl, ne)
                C = reshape(c(mask), 1, 1, []);
                A = obj.createMatrix(L.*C, TEI=EI(mask), SEI=NEI);
            end
            G = G + A;
        end
        
        function G = assembleFluxVector(obj, field, bias, func, mode, varargin)
            % ASSEMBLEFLUXVECTOR Assemble the vector contribution
            % from flux integrals for the divergence of piecewise constant
            % advection advection.
            
            % Initialization
            G = obj.createZeros(nCols=1);
            
            % Contribution from boundary elements
            BEI = obj.Space.Mesh.getBoundaryElements(); % (:, 1)
            for LFI = 1:obj.Space.Element.NFluxes
                % Find elements of which LFI-th face is on boundary
                EI = obj.Space.Mesh.findBoundaryElements(BEI, LFI);
                if ~isempty(EI)
                    GF = obj.assembleBoundaryFluxVector(field, bias, func, EI, LFI, mode, varargin{:});
                    G = G + GF;
                end
            end
        end
        
        function G = assembleBoundaryFluxVector(obj, field, bias, func, EI, LFI, mode, varargin)
            % ASSEMBLEBOUNDARYFLUXVECTOR Assemble boundary flux vector.
            
            % Outward normal vectors
            normal = obj.Space.Mesh.computeOutwardNormals(EI, LFI); % (nd, :)
            
            % Stabilizer
            eta = sum(bias.*normal, 1);
            mask = eta < 0;
            if size(bias, 2) == 1 && size(normal, 2) == 1
                mask = repmat(mask, 1, length(EI));
            end
            
            if ~any(mask)
                G = obj.createZeros(nCols=1);
                return;
            end
            
            % Face scaling factors
            ratio = obj.Space.Mesh.computeFaceScalingFactors(EI, LFI);
            ratio = reshape(ratio, 1, []); % (1, :)
            
            % Face coefficients
            C = ratio .* sum(field.*normal, 1); % (1, :)
            C = reshape(C(mask), 1, []); % (1, ne)
            
            % Test basis values and weights
            V = obj.Space.Element.Flux(LFI).Values; % (nl, nq)
            W = diag(obj.Space.Element.Flux(LFI).Weights); % (nq, nq)
            
            % Face integration points
            nq = obj.Space.Element.Flux(LFI).NPoints;
            xRef = obj.Space.Element.Flux(LFI).Nodes; % (nd, nq)
            switch mode
                case ''
                    U = obj.Space.feval(func, xRef, EI=EI(mask), args = varargin); % (nd, nq, ne)
                case 'face'
                    U = obj.Space.feval(func, xRef, EI=EI(mask), args = [LFI, varargin]); % (nd, nq, ne)
            end
            U = reshape(U, nq, []); % (nq, ne)
            
            nl = obj.Space.NLocalDofs;
            ne = length(EI);
            L = V * W * U; % (nl, ne)
            L = reshape(L, nl, []); % (nl, ne)
            L = obj.Space.Element.Approximator.project(L); % (nl, ne)
            L = reshape(L, nl, ne, []);
            L = L .* C; % (nl, ne)
            G = obj.createVector(L, TEI=EI(mask));
        end
    end
end
