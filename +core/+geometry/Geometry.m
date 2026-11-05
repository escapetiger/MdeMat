classdef Geometry < handle
    % GEOMETRY Base class for all geometries.
    %
    %   Geometry defines common properties and interfaces of a geometry.
    %
    % See also:
    %   core.geometry.Polytope, core.geometry.HypersphericalGeometry

    properties
        nDims % Geometry dimension
    end

    methods
        function obj = Geometry(nDims)
            % GEOMETRY Constructor for Geometry.
            %
            %   obj = Geometry(nDims) creates a geometric domain object
            %   representing an n-dimensional space.
            %
            % Inputs:
            %   nDims - Geometry dimension (positive integer)
            %
            % Outputs:
            %   obj - Constructed Geometry object

            obj.nDims = nDims;
        end

        function m = magnitude(obj)
            % MAGNITUDE Returns the measure (size) of the geometry.
            %
            %   m = magnitude(obj) returns the measure of the geometric
            %   domain, which depends on the dimension: length in 1D, area
            %   in 2D, volume in 3D, and hypervolume in higher dimensions.
            %   Concrete subclasses must implement this method to define
            %   how the geometric measure is calculated.
            %
            % Inputs:
            %   obj - The Geometry object
            %
            % Outputs:
            %   m - The measure (size) of the geometric domain
        end

        function TF = isInside(obj, X)
            % ISINSIDE Determines if points are inside the geometry.
            %
            %   TF = isInside(obj, X) tests whether specified points are
            %   contained within the interior of the geometric domain,
            %   excluding its boundary. Concrete subclasses must implement
            %   this method to define the containment test for their
            %   specific geometry type.
            %
            % Inputs:
            %   obj - The Geometry object
            %   X - Coordinates (matrix)
            %
            % Outputs:
            %   TF - True if the point is inside the geometry
        end

        function TF = isOnBoundary(obj, X)
            % ISONBOUNDARY Determines if points are on the boundary of the
            % geometry.
            %
            %   TF = isOnBoundary(obj, X) tests whether specified points
            %   lie exactly on the boundary of the geometric domain, rather
            %   than in its interior or exterior. Concrete subclasses must
            %   implement this method to define the boundary test for their
            %   specific geometry type.
            %
            % Inputs:
            %   obj - The Geometry object
            %   X - Coordinates (matrix)
            %
            % Outputs:
            %   TF - True if the point is on the boundary of the geometry
        end
    end
end