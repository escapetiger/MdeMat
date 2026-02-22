classdef Orthotope < core.geometry.Polytope
    % ORTHOTOPE Geometry for hyperrectangular domains.
    %
    %   Orthotope (hyperrectangle) represents a generalized geometric
    %   domain that extends the concept of rectangles and cuboids to
    %   arbitrary dimensions. It is defined by lower and upper bounds in
    %   each dimension, creating a fully rectangular geometric region.
    %
    % See also:
    %   core.geometry.Polytope, core.geometry.Simplex

    properties
        Lower % Lower bounds in each dimension
        Upper % Upper bounds in each dimension
    end

    properties (Dependent)
        NFaces % Number of boundaries (2 * nd)
        OutNormals % Outer normals of boundaries
        BBox % Bounding box in alternating format
    end

    methods
        function obj = Orthotope(bbox)
            % ORTHOTOPE Construct an instance of Orthotope.
            %
            %   obj = Orthotope(bbox) creates an orthotope (hyperrectangle)
            %   from a bounding box specification. The bounding box is
            %   defined as alternating lower and upper bounds for each
            %   dimension [a1, b1, a2, b2, ...] where ai is the lower bound
            %   and bi is the upper bound for dimension i.
            %
            %   The @a bbox parameter must have an even number of elements
            %   and all lower bounds must not exceed their corresponding
            %   upper bounds. Vertices and faces are automatically
            %   computed.

            arguments
                bbox{mustBeVector, mustBeNumeric}
            end

            core.except.assert(mod(length(bbox), 2) == 0, 'InvalidInput', ...
                'Bounding box must be a vector of even length.');

            bbox = bbox(:).';
            nd = length(bbox) / 2;

            lower = bbox(1:2:end);
            upper = bbox(2:2:end);

            core.except.assert(all(lower <= upper), 'InvalidInput', ...
                'Lower bounds must not exceed upper bounds.');

            nVertices = 2^nd;
            vertices = zeros(nd, nVertices);
            for i = 0:(nVertices - 1)
                binary = dec2bin(i, nd);
                for j = 1:nd
                    if binary(j) == '0'
                        vertices(j, i+1) = lower(j);
                    else
                        vertices(j, i+1) = upper(j);
                    end
                end
            end

            nFaces = 2 * nd;
            faces = cell(1, nFaces);
            binaries = false(nVertices, nd);
            for i = 0:(nVertices - 1)
                binStr = dec2bin(i, nd);
                for j = 1:nd
                    binaries(i+1, j) = binStr(j) == '1';
                end
            end
            nVerticesPerFace = 2^(nd - 1);
            for j = 1:nd
                f1 = zeros(1, nVerticesPerFace);
                f2 = zeros(1, nVerticesPerFace);
                k1 = 1;
                k2 = 1;
                for i = 1:2^nd
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
            obj.Lower = lower;
            obj.Upper = upper;
        end

        function n = get.NFaces(obj)
            % GET.NFACES Returns the number of boundaries.

            n = 2 * obj.NDims;
        end

        function N = get.OutNormals(obj)
            % GET.OUTNORMALS Returns the outward normals of boundaries.

            N = kron(eye(obj.NDims), [-1, 1]);
        end

        function bbox = get.BBox(obj)
            % GET.BBOX Returns the bounding box in alternating format.
            %
            %   bbox = get.BBox(obj) returns the bounding box as
            %   [a1, b1, a2, b2, ...] where ai and bi are the lower and
            %   upper bounds for dimension i.

            bbox = reshape([obj.Lower; obj.Upper], [], 1);
        end

        function m = magnitude(obj)
            % MAGNITUDE Returns the volume of the orthotope.
            %
            %   m = magnitude(obj) computes the n-dimensional volume of the
            %   orthotope as the product of the extent in each dimension.
            %   For 1D this is length, for 2D this is area, for 3D this is
            %   volume, and for higher dimensions this is hypervolume.

            arguments
                obj core.geometry.Orthotope
            end

            m = prod(obj.Upper-obj.Lower);
        end

        function TF = isInside(obj, X)
            % ISINSIDE Tests if points are inside the orthotope.
            %
            %   TF = isInside(obj, X) determines whether points lie
            %   strictly within the interior of the orthotope, excluding
            %   boundary faces. Points are inside if they satisfy ai < xi <
            %   bi for all dimensions i.

            arguments
                obj core.geometry.Orthotope
                X{mustBeNumeric}
            end

            nd = obj.NDims;
            nx = size(X);
            core.except.assert(nx(1) == nd, 'DimensionMismatch', ...
                'Point dimension must match geometry dimension.');

            a = obj.Lower;
            b = obj.Upper;
            X = reshape(X, nx(1), []);
            TF = true(1, size(X, 2));
            for i = 1:nd
                TF = TF & (X(i, :) > a(i)) & (X(i, :) < b(i));
            end

            if length(nx) > 2, TF = reshape(TF, nx(2:end)); end
        end

        function TF = isOnBoundary(obj, X)
            % ISONBOUNDARY Tests if points are on the orthotope boundary.
            %
            %   TF = isOnBoundary(obj, X) determines whether points lie
            %   exactly on the boundary faces of the orthotope. Points are
            %   on the boundary if they are within the orthotope bounds
            %   but not strictly inside.

            arguments
                obj core.geometry.Orthotope
                X{mustBeNumeric}
            end

            nd = obj.NDims;
            nx = size(X);
            core.except.assert(nx(1) == nd, 'DimensionMismatch', ...
                'Point dimension must match geometry dimension.');

            a = obj.Lower;
            b = obj.Upper;
            X = reshape(X, nx(1), []);
            TF = true(1, size(X, 2));
            for i = 1:nd
                TF = TF & (X(i, :) >= a(i)) & (X(i, :) <= b(i));
            end

            if length(nx) > 2, TF = reshape(TF, nx(2:end)); end

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

            arguments
                nDims(1, 1) {mustBePositive, mustBeInteger}
            end

            obj = core.geometry.Orthotope(repmat([-1 / 2, 1 / 2], 1, nDims));
        end
    end
end