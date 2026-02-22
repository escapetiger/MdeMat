classdef FiniteElementTraceAssembly < approx.assembly.FiniteElementAssembly
    % FINITEELEMENTTRACEASSEMBLY Assembly for finite element trace operations.
    %
    %   FiniteElementTraceAssembly provides specialized assembly methods for
    %   trace operations in finite element computations. These operations are
    %   particularly important for handling boundary conditions and interface
    %   problems where contributions must be accumulated from traces of elements
    %   on boundaries.
    %
    %   The class handles vectorized assembly operations (batch processing) for
    %   efficiency when multiple components need to be assembled simultaneously.
    %   All assembly operations respect periodic boundary conditions by skipping
    %   assembly when such conditions are specified.
    %
    % See also:
    %   FiniteElementAssembly, FiniteElementEllipticAssembly
    
    properties
        BcType % Boundary condition type
    end
    
    methods
        function G = assembleVectorBatch(obj, nc, func, options)
            % ASSEMBLEVECTORBATCH Assemble batch vector contributions from trace elements.
            %
            %   G = assembleVectorBatch(obj, nc, func) assembles vector contributions
            %   from trace elements for @a nc components using function @a func.
            %   Returns a cell array of sparse vectors, one for each spatial dimension.
            %
            %   G = assembleVectorBatch(obj, nc, func, options) allows additional
            %   arguments to be passed to the function evaluation through @a options.args.
            
            arguments
                obj approx.assembly.FiniteElementTraceAssembly
                nc {mustBePositive, mustBeInteger}
                func {mustBeA(func, 'function_handle')}
                options.args cell = {}
            end
            
            if strcmpi(obj.bcType, 'periodic')
                G = obj.createZeros(nCols=nc);
                return;
            end
            
            G = obj.assembleFluxVectorBatch(nc, func, options.args{:});
        end
    end
    
    methods (Access = private)
        function G = assembleFluxVectorBatch(obj, nc, func, varargin)
            % ASSEMBLEFLUXVECTORBATCH Assemble flux vector contributions in batch mode.
            
            % Initialization
            nd = obj.Space.NDims;
            G = repmat({obj.createZeros(nCols=nc)}, 1, nd);
            
            % Contribution from boundary elements
            EIB = obj.Space.Mesh.getBoundaryElements(); % (:, 1)
            for LFI = 1:obj.Space.Element.NFluxes
                % Find elements of which LFI-th face is on boundary
                EI = obj.Space.Mesh.findBoundaryElements(EIB, LFI);
                if ~isempty(EI)
                    GF = obj.assembleTraceVectorBatch(func, nc, EI, LFI, varargin{:});
                    for d = 1:nd
                        G{d} = G{d} + GF{d};
                    end
                end
            end
        end
        
        function G = assembleTraceVectorBatch(obj, func, nc, EI, LFI, varargin)
            % ASSEMBLETRACEVECTORBATCH Assemble trace contributions for specific elements.
            
            % Sizes
            nd = obj.Space.NDims;
            nl = obj.Space.NLocalDofs;
            
            % Outward normal vectors
            normal = obj.Space.Mesh.computeOutwardNormals(EI, LFI); % (nd, :)
            
            % Face scaling factors
            ratio = obj.Space.Mesh.computeFaceScalingFactors(EI, LFI); % (:, 1)
            ratio = reshape(ratio, 1, []); % (1, :)
            
            % Face coefficients
            FC = ratio .* normal;
            
            % Test basis values and weights
            V = obj.Space.Element.Flux(LFI).Values; % (nl, nq)
            W = diag(obj.Space.Element.Flux(LFI).Weights); % (nq, nq)
            
            % Face integration points
            nq = obj.Space.Element.Flux(LFI).NPoints;
            xRef = obj.Space.Element.Flux(LFI).Nodes; % (nd, nq)
            U = obj.Space.feval(func, xRef, EI=EI, args=varargin);
            U = reshape(U, nq, []); % (nq, ne*nc)
            
            ne = obj.Space.NMeshElements;
            G = repmat({obj.createZeros(nCols=nc)}, 1, nd);
            for d = 1:nd
                L = V * W * U; % (nl, ne*nc)
                L = obj.Space.Element.Approximator.project(L); % (nl, ne, nc)
                L = reshape(L, nl, [], nc); % (nl, ne)
                if size(L, 2) == 1, L = repmat(L, 1, ne, 1); end
                C = reshape(FC(d, :), 1, []); % (1, ne)
                G{d} = obj.createVectorBatch(L.*C, TEI=EI);
            end
        end
    end
end

