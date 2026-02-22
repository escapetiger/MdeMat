classdef Intersector < handle
    % INTERSECTOR Base class for all geometric intersectors.
    %
    %   Intersector defines the interface and common functionality for
    %   computing intersections between geometric objects in n-dimensional
    %   space. This abstract class serves as the foundation for specific
    %   intersection algorithms.
    %
    % See also:
    %   core.geometry.LineSegmentByLineSegmentIntersector,
    %   core.geometry.LineSegmentByPlaneIntersector

    properties
        nDims % Dimension of the space (positive integer)
    end

    methods
        function obj = Intersector(nDims)
            % INTERSECTOR Constructor for Intersector object.
            %
            %   obj = Intersector(nDims) creates an Intersector object
            %   for operations in nDims-dimensional space.
            %
            % Inputs:
            %   nDims - Dimension of the space (positive integer)
            %
            % Outputs:
            %   obj - The constructed Intersector object
            %
            % Examples:
            %   % Create intersector for 3D space
            %   obj = ConcreteIntersector(3);
            
            core.except.assert(nDims >= 1, ...
                'InvalidInput', 'Dimension must be a positive integer.');

            obj.nDims = nDims;
        end
    end
    
    methods (Abstract)
        % INTERSECT Computes the intersection between geometric objects.
        %
        %   [X, TF] = intersect(obj, varargin) computes intersection points
        %   between geometric objects. The specific implementation depends
        %   on the concrete subclass.
        %
        % Inputs:
        %   obj - The Intersector object
        %   varargin - Variable input arguments specific to intersection type
        %
        % Outputs:
        %   X - Matrix with each column representing an intersection point
        %   TF - Logical array indicating valid intersections
        [X, TF] = intersect(obj, varargin)
    end
end