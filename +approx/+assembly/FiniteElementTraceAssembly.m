classdef FiniteElementTraceAssembly < approx.assembly.FiniteElementAssembly
    %FINITEELEMENTTRACEASSEMBLY 
    
    methods
        function G = assembleVectorBatch(obj, nc, func, varargin)
            %< Skip for periodic boundary condition
            if strcmpi(obj.bcType, 'periodic')
                G = obj.createZeros(nc);
                return;
            end

            G = obj.assembleFluxVectorBatch(func, varargin{:});
        end
    end

    methods (Access = private)
        function G = assembleFluxVectorBatch(obj, nc, func, varargin)
            %< Initialization
            nd = obj.space.nDims;
            G = repmat({obj.createZeros(nc)}, 1, nd);

            %< Contribution from boundary elements
            IEB = obj.space.mesh.getBoundaryElements(); % [] x 1
            for IFE = 1:obj.space.element.nFluxes
                %< Find elements of which IFE-th face is on boundary
                IE = obj.space.mesh.findBoundaryElements(IEB, IFE);
                if ~isempty(IE)
                    GF = obj.assembleTraceVectorBatch(nc, IE, IFE, func, varargin{:});
                    for d = 1:nd
                        G{d} = G{d} + GF{d};
                    end
                end
            end
        end

        function G = assembleTraceVectorBatch(obj, nc, IE, IFE, func, varargin)
            %< Sizes
            nd = obj.space.nDims;
            nl = obj.space.nLocalDofs;

            %< Outward normal vectors
            normal = obj.space.mesh.computeOutwardNormals(IE, IFE); % d x []

            %< Face scaling factors
            ratio = obj.space.mesh.computeFaceScalingFactors(IE, IFE); % [] x 1
            ratio = reshape(ratio, 1, []);

            %< Face coefficients
            FC = ratio .* normal;

            %< Test basis values and weights
            V = obj.space.element.flux(IFE).values; % nl x nq
            W = diag(obj.space.element.flux(IFE).weights); % nq x nq

            %< Face integration points
            nq = obj.space.element.flux(IFE).nPoints;
            xRef = obj.space.element.flux(IFE).nodes; % d x nq
            U = obj.space.functionEvaluate(IE, xRef, func, varargin{:});
            U = reshape(U, nq, []); % nq x ne*nc

            ne = obj.space.nMeshElements;
            G = repmat({obj.createZeros(nc)}, 1, nd);
            for d = 1:nd
                L = V * W * U; % nl x ne*nc
                L = obj.space.element.approximator.project(L); % nl x ne x nc
                L = reshape(L, nl, [], nc); % nl x ne
                if size(L, 2) == 1, L = repmat(L, 1, ne, 1); end
                C = reshape(FC(d, :), 1, []); % 1 x ne
                G{d} = obj.createVectorBatch(IE, [], L.*C);
            end
        end
    end
end

