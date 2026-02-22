classdef Geometry < handle
    % GEOMETRY Base class for all geometries.
    %
    %   Geometry defines the common properties and interfaces for all
    %   geometric objects in the MdeMat library. This abstract base class
    %   establishes the fundamental interface that all geometric entities
    %   must implement, including magnitude calculation and point
    %   containment testing. Concrete subclasses implement specific
    %   geometry types such as polytopes and hyperspherical geometries.
    %
    %   This is an abstract base class that cannot be instantiated
    %   directly. All geometric objects inherit from this class and must
    %   implement the abstract methods magnitude(), isInside(), and
    %   isOnBoundary().
    %
    % See also:
    %   core.geometry.Polytope, core.geometry.HypersphericalGeometry

    properties
        NDims % Number of dimensions
    end

    methods
        function obj = Geometry(nDims)
            % GEOMETRY Construct an instance of Geometry.
            %
            %   obj = Geometry(nDims) creates a geometric domain object
            %   representing an n-dimensional space with dimension @a nDims.

            arguments
                nDims (1,1) {mustBeNonnegative, mustBeInteger}
            end
            
            obj.NDims = nDims;
        end

        function m = magnitude(obj) 
            % MAGNITUDE Returns the measure (size) of the geometry.
            %
            %   m = magnitude(obj) returns the measure of the geometric
            %   domain, which depends on the dimension: length in 1D, area
            %   in 2D, volume in 3D, and hypervolume in higher dimensions.

            arguments
                obj core.geometry.Geometry
            end

            m = [];
            core.except.assert(0, 'NotImplemented', ...
                'Method isOnBoundary is not implemented by class %s.', class(obj));
        end

        function TF = isInside(obj, X)
            % ISINSIDE Determines if points are inside the geometry.
            %
            %   TF = isInside(obj, X) tests whether specified points @a X
            %   are contained within the interior of the geometric domain,
            %   excluding its boundary. @a X should be an nd-by-np matrix
            %   where nd is the number of dimensions and np is the number
            %   of points.

            arguments
                obj core.geometry.Geometry
                X {mustBeNumeric}
            end

            TF = logical(X);
            core.except.assert(0, 'NotImplemented', ...
                'Method isInside is not implemented by class %s.', class(obj));
        end

        function TF = isOnBoundary(obj, X)
            % ISONBOUNDARY Determines if points are on the boundary of the
            % geometry.
            %
            %   TF = isOnBoundary(obj, X) tests whether specified points @a
            %   X lie exactly on the boundary of the geometric domain. @a X
            %   should be an nd-by-np matrix where nd is the number of
            %   dimensions and np is the number of points.

            arguments
                obj core.geometry.Geometry
                X {mustBeNumeric}
            end

            TF = logical(X);
            core.except.assert(0, 'NotImplemented', ...
                'Method isOnBoundary is not implemented by class %s.', class(obj));
        end
    end
end