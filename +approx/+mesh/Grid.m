classdef Grid < approx.mesh.Mesh
    % GRID Base class for multidimensional structured grids.
    %
    %   Grid provides a foundational framework for representing structured
    %   computational grids using memory-efficient tensor-based storage.
    %   Unlike unstructured meshes that store explicit connectivity,
    %   structured grids leverage their regular topology to achieve O(d×n)
    %   memory usage instead of O(n^d) for d-dimensional grids with n
    %   elements per dimension.
    %
    % See also:
    %   approx.mesh.Mesh, approx.mesh.SeparableGraph
    
    properties
        Resolution % Element count per dimension (1 × NDims vector)
        Centroids % Element center coordinates per dimension (1 × NDims cell)
        Spacings % Element spacing per dimension (1 × NDims cell)
        Boundary % Boundary element specifications (1 × 2*NDims cell)
        Indexer % Multi-index converter (core.linalg.MultiIndexer)
    end
    
    properties (Dependent)
        BBox % Bounding box
        AllElementMultiIndices % Multi-indices of all elements
    end
    
    methods
        function obj = Grid(vertices)
            % GRID Constructor for Grid base class.
            %
            %   obj = Grid(vertices) creates a non-uniform grid from vertex
            %   coordinate vectors that define element boundaries in each
            %   dimension, automatically computing element properties.

            arguments
                vertices {mustBeA(vertices, 'cell')}
            end
            
            elements = cellfun(@(x) 1:(length(x) - 1), vertices, 'Un', 0);
            obj@approx.mesh.Mesh(vertices, elements);
            obj.Resolution = cellfun(@(x) length(x)-1, vertices);
            obj.Centroids = cellfun(@(x) (x(1:end-1) + x(2:end))/2, vertices, 'Un', 0);
            obj.Spacings = cellfun(@(x) diff(x), vertices, 'Un', 0);
            obj.Indexer = core.linalg.MultiIndexer(shape=obj.Resolution);
            obj.Boundary = cell(1, 2*obj.NDims);
            for i = 1:obj.NDims
                f = arrayfun(@(n) 1:n, obj.Resolution, 'Un', 0);
                f{i} = 1;
                obj.Boundary{2*i-1} = f;
                f{i} = obj.Resolution(i);
                obj.Boundary{2*i} = f;
            end
        end
        
        function IEI = getInternalElements(obj)
            % GETINTERIORELEMENTS Get linear indices of internal elements.
            %
            %   IEI = getInternalElements(obj) computes linear indices of
            %   internal elements that do not touch any boundary.

            arguments
                obj approx.mesh.Grid
            end
            
            ME = cellfun(@(I) I(2:end-1), obj.Elements, 'Un', 0);
            ME = obj.Indexer.factorToMulti(ME);
            IEI = obj.Indexer.multiToLinear(ME);
        end
        
        function BEI = getBoundaryElements(obj, LFI)
            % GETBOUNDARYELEMENTS Get linear indices of boundary elements.
            %
            %   BEI = getBoundaryElements(obj) computes the linear indices
            %   of all elements that have at least one face on the
            %   boundary.
            %
            %   BEI = getBoundaryElements(obj, LFI) finds elements whose
            %   LFI-th face is on the boundary.

            arguments
                obj approx.mesh.Grid
                LFI {mustBeInteger, mustBePositive} = []
            end
            
            if isempty(LFI)
                nb = length(obj.Boundary);
                BEI = cell(1, nb);
                for i = 1:nb
                    M = obj.Indexer.factorToMulti(obj.Boundary{i});
                    BEI{i} = obj.Indexer.multiToLinear(M);
                end
                BEI = vertcat(BEI{:});
                BEI = unique(BEI, 'stable');
                return;
            end
            
            ME = obj.Indexer.factorToMulti(obj.Boundary{LFI});
            BEI = obj.Indexer.multiToLinear(ME);
        end
        
        function IEI = findInternalElements(obj, EI, LFI)
            % FINDINTERIORELEMENTS Find internal elements across faces.
            %
            %   IEI = findInternalElements(obj, EI, LFI) finds elements
            %   whose LFI-th face is not on the boundary within the
            %   specified elements @a EI.

            arguments
                obj approx.mesh.Grid
                EI {mustBeInteger}
                LFI {mustBeInteger, mustBePositive}
            end
            
            ME = obj.Indexer.linearToMulti(EI);
            dim = ceil(LFI/2);
            if mod(LFI, 2) == 1
                mask = ME(:, dim) > 1;
            else
                mask = ME(:, dim) < obj.Resolution(dim);
            end
            IEI = obj.Indexer.multiToLinear(ME(mask, :));
        end
        
        function BEI = findBoundaryElements(obj, EI, LFI)
            % FINDBOUNDARYELEMENTS Find boundary elements across faces.
            %
            %   BEI = findBoundaryElements(obj, EI, LFI) finds elements
            %   whose LFI-th face is on the boundary within the specified
            %   elements @a EI.

            arguments
                obj approx.mesh.Grid
                EI {mustBeInteger}
                LFI {mustBeInteger, mustBePositive}
            end
            
            ME = obj.Indexer.linearToMulti(EI);
            dim = ceil(LFI/2);
            if mod(LFI, 2) == 1
                mask = ME(:, dim) == 1;
            else
                mask = ME(:, dim) == obj.Resolution(dim);
            end
            BEI = obj.Indexer.multiToLinear(ME(mask, :));
        end
        
        function [NEI, NLFI] = findNeighborElements(obj, EI, LFI, BC)
            % FINDNEIGHBORELEMENTS Find neighboring elements across faces.
            %
            %   [NEI, NLFI] = findNeighborElements(obj, EI, LFI, BC) computes
            %   neighbor elements sharing the specified local face of given
            %   elements, with appropriate boundary condition treatment.

            arguments
                obj approx.mesh.Grid
                EI {mustBeInteger}
                LFI {mustBeInteger, mustBePositive}
                BC {mustBeMember(BC, {'', 'periodic', 'dirichlet'})} = ''
            end
            
            ME = obj.Indexer.linearToMulti(EI);
            
            dim = ceil(LFI/2);
            
            if mod(LFI, 2) == 1
                ME(:, dim) = ME(:, dim) - 1;
                NLFI = 2 * dim;
            else
                ME(:, dim) = ME(:, dim) + 1;
                NLFI = 2 * dim - 1;
            end
            
            switch lower(BC)
                case 'periodic'
                    NEI = obj.Indexer.multiToLinear(ME, pad='wrap');
                case 'dirichlet'
                    NEI = obj.Indexer.multiToLinear(ME, pad='empty');
                otherwise
                    NEI = obj.Indexer.multiToLinear(ME, pad='empty');
            end
        end
        
        function graphObj = graphify(obj)
            % GRAPHIFY Convert grid to separable graph representation.
            %
            %   graphObj = graphify(obj) creates a SeparableGraph that
            %   represents the grid's connectivity structure using 1D
            %   factors with variable vertex spacing.

            arguments
                obj approx.mesh.Grid
            end
            
            d = obj.NDims;
            n = obj.Resolution;
            h = obj.Spacings;
            c = obj.Centroids;
            V = arrayfun(@(i) [c{i} - h{i} / 2, c{i}(end) + h{i}(end) / 2].', 1:d, 'Un', 0);
            E = arrayfun(@(k) [(1:k).', (2:(k + 1)).'], n, 'Un', 0);
            graphObj = approx.mesh.SeparableGraph(V, E);
        end
        
        function Y = collocate(obj, X, EI)
            % COLLOCATE Map reference coordinates to physical coordinates.
            %
            %   Y = collocate(obj, X) maps reference coordinates from the
            %   standard hypercube $[-1/2,1/2]^nd$ to physical coordinates
            %   across all grid elements using element-specific
            %   transformations.
            %
            %   Y = collocate(obj, X, EI) maps coordinates to specific
            %   elements identified by linear indices EI.

            arguments
                obj approx.mesh.Grid
                X {mustBeNumericOrCell}
                EI {mustBeInteger} = []
            end
            
            nd = obj.NDims;
            if ~iscell(X)
                C = cell(1, nd);
                [C{:}] = ndgrid(obj.Centroids{:});
                C = reshape(cat(nd+1, C{:}), [], nd).';
                if ~isempty(EI), C = C(:, EI); end
                C = reshape(C, nd, 1, []);
                
                h = cell(1, nd);
                [h{:}] = ndgrid(obj.Spacings{:});
                h = reshape(cat(nd+1, h{:}), [], nd).';
                if ~isempty(EI), h = h(:, EI); end
                h = reshape(h, nd, 1, []);
                
                X = reshape(X, nd, [], 1);
                Y = C + h .* X;
                Y = reshape(Y, nd, []);
            else
                Y = cell(nd, 1);
                if isempty(EI)
                    for i = 1:nd
                        C = obj.Centroids{i};
                        h = obj.Spacings{i};
                        Y{i} = C(:).' + h(:).' .* X{i}(:);
                        Y{i} = Y{i}(:);
                    end
                else
                    M = obj.Indexer.linearToMulti(EI);
                    for i = 1:nd
                        C = obj.Centroids{i}(M(:, i));
                        h = obj.Spacings{i}(M(:, i));
                        Y{i} = C(:).' + h(:).' .* X{i}(:);
                        Y{i} = Y{i}(:);
                    end
                end
            end
        end
        
        function newObj = refine(obj, nLevels)
            % REFINE Create refined grid preserving spacing patterns.
            %
            %   newObj = refine(obj, nLevels) creates a refined grid by
            %   uniformly subdividing each element while preserving the
            %   relative spacing characteristics of the original grid.

            arguments
                obj approx.mesh.Grid
                nLevels {mustBeInteger, mustBeNonnegative} = 1
            end
            
            if nLevels == 0
                newObj = obj;
                return;
            end
            
            nd = obj.NDims;
            p = 2^nLevels;
            x = cell(1, nd);
            for i = 1:nd
                y = obj.Centroids{i};
                n = obj.Resolution(i);
                s = obj.Spacings{i};
                a = y - s / 2;
                b = y + s / 2;
                h = (b - a) ./ p;
                x{i} = zeros(1, 1+n*p);
                x{i}(1) = a(1);
                for j = 1:n
                    x{i}(1 + (j - 1) * p + (1:p)) = (a(j) + h(j)):h(j):b(j);
                end
            end
            newObj = approx.mesh.Grid(x);
        end
        
        function h = computeMeasure(obj)
            % COMPUTEMEASURE Compute minimum element spacing across the grid.
            %
            %   h = computeMeasure(obj) returns the minimum element spacing
            %   among all elements in all dimensions, providing a
            %   characteristic length scale for the grid.

            arguments
                obj approx.mesh.Grid
            end
            
            h = min(cellfun(@min, obj.Spacings));
        end
        
        function detJ = computeAllElementJacobianDeterminants(obj)
            % COMPUTEALLELEMENTJACOBIANDETERMINANTS Compute
            % element-specific Jacobian determinants.
            %
            %   detJ = computeAllElementJacobianDeterminants(obj) computes
            %   the Jacobian determinant for coordinate transformation for
            %   each element. For structured grids, determinants vary by
            %   element based on local spacing.

            arguments
                obj approx.mesh.Grid
            end
            
            h = cell(1, obj.NDims);
            [h{:}] = ndgrid(obj.Spacings{:});
            h = cellfun(@(x) x(:).', h, 'Un', 0);
            detJ = prod(cat(1, h{:}), 1);
            detJ = detJ(:);
        end
        
        function J = computeAllElementJacobians(obj)
            % COMPUTEALLELEMENTJACOBIANS Compute element-specific Jacobian
            % matrices.
            %
            %   J = computeAllElementJacobians(obj) computes the Jacobian
            %   matrix for coordinate transformation for each element. For
            %   structured grids, Jacobians vary by element but remain
            %   diagonal.

            arguments
                obj approx.mesh.Grid
            end
            
            nd = obj.NDims;
            h = cell(1, nd);
            [h{:}] = ndgrid(obj.Spacings{:});
            h = cellfun(@(x) x(:).', h, 'Un', 0);
            h = cat(1, h{:});
            J = reshape(h, nd, 1, []) .* eye(nd, nd);
        end
        
        function invJ = computeAllElementInverseJacobians(obj)
            % COMPUTEALLELEMENTINVERSEJACOBIANS Compute element-specific
            % inverse Jacobian matrices.
            %
            %   invJ = computeAllElementInverseJacobians(obj) computes the
            %   inverse Jacobian matrix for coordinate transformation for
            %   each element. For structured grids, inverse Jacobians vary
            %   by element.

            arguments
                obj approx.mesh.Grid
            end
            
            nd = obj.NDims;
            h = cell(1, nd);
            [h{:}] = ndgrid(obj.Spacings{:});
            h = cellfun(@(x) x(:).', h, 'Un', 0);
            h = 1 ./ cat(1, h{:});
            invJ = reshape(h, nd, 1, []) .* eye(nd, nd);
        end
        
        function detFJac = computeAllFaceJacobianDeterminants(obj)
            % COMPUTEALLFACEJACOBIANDETERMINANTS Compute element-specific
            % face Jacobian determinants.
            %
            %   detFJac = computeAllFaceJacobianDeterminants(obj) computes
            %   Jacobian determinants for coordinate transformation on
            %   element faces. For structured grids, face determinants vary
            %   by element.

            arguments
                obj approx.mesh.Grid
            end
            
            nd = obj.NDims;
            if nd == 1
                detFJac = repmat({ones(obj.NElements, 1)}, 1, 2);
                return;
            end
            
            detFJac = cell(1, 2*nd);
            for i = 1:2 * nd
                dim = ceil(i/2);
                h = cell(1, nd);
                [h{:}] = ndgrid(obj.Spacings{:});
                h = cellfun(@(x) x(:).', h(setdiff(1:nd, dim)), 'Un', 0);
                detFJac{i} = prod(cat(1, h{:}), 1).';
            end
        end
        
        function normal = computeAllOutwardNormals(obj)
            % COMPUTEALLOUTWARDNORMALS Compute outward normal vectors for
            % all faces.
            %
            %   normal = computeAllOutwardNormals(obj) computes unit outward
            %   normal vectors for each face type. For structured grids,
            %   normals remain axis-aligned and constant for faces of the
            %   same type.

            arguments
                obj approx.mesh.Grid
            end
            
            nd = obj.NDims;
            normal = cell(1, 2*nd);
            for i = 1:(2 * nd)
                dim = ceil(i/2);
                normal{i} = zeros(nd, 1);
                if mod(i, 2) == 0
                    normal{i}(dim) = +1;
                else
                    normal{i}(dim) = -1;
                end
            end
        end
        
        function v = computeUniqueDirection(obj)
            % COMPUTEUNIQUEDIRECTION Compute a unique direction vector.
            %
            %   v = computeUniqueDirection(obj) computes a unique unit
            %   direction vector to define the "plus" and "minus" side of a
            %   face.

            arguments
                obj approx.mesh.Grid
            end
            
            nd = obj.NDims;
            if nd == 1
                v = 1;
            else
                v = ones(nd, 1) / sqrt(nd);
            end
        end
    end
    
    methods
        function BBox = get.BBox(obj)
            % GET.BBOX Get the bounding box.

            xmin = cellfun(@(x) min(x), obj.Vertices);
            xmax = cellfun(@(x) max(x), obj.Vertices);
            BBox = reshape([xmin; xmax], 1, []);
        end

        function M = get.AllElementMultiIndices(obj)
            % GET.ALLELEMENTMULTIINDICES Get multi-indices of all elements.

            M = obj.Indexer.factorToMulti(obj.Elements);
        end
    end
end

function mustBeNumericOrCell(x)
    if ~(isnumeric(x) || iscell(x))
        error('Input must be numeric or cell.');
    end
end