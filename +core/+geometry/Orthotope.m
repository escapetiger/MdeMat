classdef Orthotope < core.geometry.Polytope
    % ORTHOTOPE Geometry for hyperrectangular domains.
    %
    %   Orthotope (hyperrectangle) represents a generalized geometric
    %   domain that extends the concept of rectangles and cuboids to
    %   arbitrary dimensions. It is defined by lower and upper bounds in
    %   each dimension, creating a fully rectangular geometric region.
    %
    % Examples:
    %   % Create a 2D rectangle from [1,3] to [2,4]
    %   rect = core.geometry.Orthotope([1, 2, 3, 4]);
    %   vol = rect.magnitude();  % Returns 2 (area)
    %
    %   % Create a unit square in 2D
    %   rect = core.geometry.Orthotope.unit(2);
    %
    % See also:
    %   core.geometry.Polytope, core.geometry.Simplex

    properties
        lower % Lower bounds in each dimension
        upper % Upper bounds in each dimension
    end

    properties (Dependent)
        nBoundaries % Number of boundaries (2 * nDims)
        outerNormals % Outer normals of boundaries
        bbox % Bounding box in alternating format
    end

    methods
        function obj = Orthotope(bbox)
            % ORTHOTOPE Constructor for Orthotope.
            %
            %   obj = Orthotope(bbox) creates an orthotope (hyperrectangle)
            %   from a bounding box specification. The bounding box is
            %   defined as alternating lower and upper bounds for each
            %   dimension [a1, b1, a2, b2, ...] where ai is the lower bound
            %   and bi is the upper bound for dimension i.
            %
            % Inputs:
            %   bbox - Bounding box vector [a1, b1, a2, b2, ...]
            %
            % Outputs:
            %   obj - Constructed Orthotope object

            core.except.assert(nargin >= 1, 'InvalidInput', ...
                'Bounding box must be specified.');

            core.except.assert(isvector(bbox), 'InvalidInput', ...
                'Bounding box must be a vector.');

            core.except.assert(mod(length(bbox), 2) == 0, 'InvalidInput', ...
                'Bounding box must have even number of elements.');

            core.except.assert(length(bbox) >= 2, 'InvalidInput', ...
                'Bounding box must have at least 2 elements.');

            bbox = bbox(:)';
            nDims = length(bbox) / 2;

           
            lower = bbox(1:2:end);
            upper = bbox(2:2:end);

            core.except.assert(all(lower <= upper), 'InvalidInput', ...
                'Lower bounds must not exceed upper bounds.');

            nVertices = 2^nDims;
            vertices = zeros(nDims, nVertices);
            for i = 0:(nVertices - 1)
                binary = dec2bin(i, nDims);
                for j = 1:nDims
                    if binary(j) == '0'
                        vertices(j, i+1) = lower(j);
                    else
                        vertices(j, i+1) = upper(j);
                    end
                end
            end

            nFaces = 2 * nDims;
            faces = cell(1, nFaces);
            binaries = false(nVertices, nDims);
            for i = 0:(nVertices - 1)
                binStr = dec2bin(i, nDims);
                for j = 1:nDims
                    binaries(i+1, j) = binStr(j) == '1';
                end
            end
            nVerticesPerFace = 2^(nDims - 1);
            for j = 1:nDims
                f1 = zeros(1, nVerticesPerFace);
                f2 = zeros(1, nVerticesPerFace);
                k1 = 1;
                k2 = 1;
                for i = 1:2^nDims
                    if binaries(i, j) == false
                        f1(k1) = i;
                        k1 = k1 + 1;
                    else
                        f2(k2) = i;
                        k2 = k2 + 1;
                    end
                end
                faces{2*j-1} = f1;
                faces{2*j} = f2;
            end

            obj@core.geometry.Polytope(vertices, faces);
            obj.lower = lower;
            obj.upper = upper;
        end

        function n = get.nBoundaries(obj)
            % GET.NBOUNDARIES Returns the number of boundaries.
            %
            %   n = get.nBoundaries(obj) returns the number of boundary
            %   faces of the orthotope. Each dimension contributes two
            %   boundaries (lower and upper), so the total number is
            %   2 * nDims.
            %
            % Outputs:
            %   n - Number of boundaries (2 * nDims)

            n = 2 * obj.nDims;
        end

        function N = get.outerNormals(obj)
            % GET.OUTERNORMALS Returns the outer normals of boundaries.
            %
            %   N = get.outerNormals(obj) returns a matrix where each
            %   column is the outward normal vector for a boundary face
            %   of the orthotope. The normals are unit vectors pointing
            %   outward from the orthotope.
            %
            % Outputs:
            %   N - Matrix of outer normal vectors (nDims × nBoundaries)

            N = kron(eye(obj.nDims), [-1, 1]);
        end

        function bbox = get.bbox(obj)
            % GET.BBOX Returns the bounding box in alternating format.
            %
            %   bbox = get.bbox(obj) returns the bounding box as
            %   [a1, b1, a2, b2, ...] where ai and bi are the lower and
            %   upper bounds for dimension i.
            %
            % Outputs:
            %   bbox - Bounding box vector in alternating format
            
            bbox = reshape([obj.lower; obj.upper], [], 1);
        end

        function m = magnitude(obj)
            % MAGNITUDE Returns the volume of the orthotope.
            %
            %   m = magnitude(obj) computes the n-dimensional volume of the
            %   orthotope as the product of the extent in each dimension.
            %   For 1D this is length, for 2D this is area, for 3D this is
            %   volume, and for higher dimensions this is hypervolume.
            %
            % Inputs:
            %   obj - The Orthotope object
            %
            % Outputs:
            %   m - Volume of the orthotope

            m = prod(obj.upper-obj.lower);
        end

        function TF = isInside(obj, X)
            % ISINSIDE Tests if points are inside the orthotope.
            %
            %   TF = isInside(obj, X) determines whether points lie strictly
            %   within the interior of the orthotope, excluding boundary
            %   faces. Points are inside if they satisfy ai < xi < bi for
            %   all dimensions i.
            %
            % Inputs:
            %   obj - The Orthotope object
            %   X - A matrix of size d×m where each column represents a point
            %
            % Outputs:
            %   TF - A logical vector indicating which points are inside

            d = obj.nDims;
            n = size(X);
            core.except.assert(n(1) == d, 'DimensionMismatch', ...
                'Point dimension must match geometry dimension.');

            a = obj.lower;
            b = obj.upper;
            X = reshape(X, n(1), []);
            TF = true(1, size(X, 2));
            for i = 1:d
                TF = TF & (X(i, :) > a(i)) & (X(i, :) < b(i));
            end

            if length(n) > 2, TF = reshape(TF, n(2:end)); end
        end

        function TF = isOnBoundary(obj, X)
            % ISONBOUNDARY Tests if points are on the orthotope boundary.
            %
            %   TF = isOnBoundary(obj, X) determines whether points lie
            %   exactly on the boundary faces of the orthotope. Points are
            %   on the boundary if they are within the orthotope bounds
            %   but not strictly inside.
            %
            % Inputs:
            %   obj - The Orthotope object
            %   X - A matrix of size d×m where each column represents a point
            %
            % Outputs:
            %   TF - A logical vector indicating which points are on boundary

            d = obj.nDims;
            n = size(X);
            core.except.assert(n(1) == d, 'DimensionMismatch', ...
                'Point dimension must match geometry dimension.');

            a = obj.lower;
            b = obj.upper;
            X = reshape(X, n(1), []);
            TF = true(1, size(X, 2));
            for i = 1:d
                TF = TF & (X(i, :) >= a(i)) & (X(i, :) <= b(i));
            end

            if length(n) > 2, TF = reshape(TF, n(2:end)); end

            TF = TF & ~obj.isInside(X);
        end

    end

    methods (Static)
        function obj = unit(nDims)
            % UNIT Creates a unit orthotope.
            %
            %   obj = unit(nDims) creates a unit orthotope centered at the
            %   origin with extent 1 in each dimension (from -1/2 to 1/2).
            %   This is useful for creating standard reference geometries.
            %
            % Inputs:
            %   nDims - Number of dimensions for the unit orthotope
            %
            % Outputs:
            %   obj - Unit orthotope with bounds [-1/2, 1/2] in each dimension
            
            obj = core.geometry.Orthotope(repmat([-1/2, 1/2], 1, nDims));
        end
    end
end