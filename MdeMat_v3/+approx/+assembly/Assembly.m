classdef Assembly < handle
    % ASSEMBLY Base class for assembly operations.
    %
    %   Assembly provides the foundation for assembling discrete operators
    %   from finite-dimensional function spaces. The assembly process
    %   involves iterating over mesh elements, computing local data, and
    %   accumulating contributions into global sparse data using
    %   appropriate connectivity information from the mesh space.
    %
    % See also:
    %   approx.space.MeshSpace

    properties
        space % MeshSpace object representing the computational space
    end

    methods
        function obj = Assembly(space)
            % ASSEMBLY Constructor for Assembly base class.
            %
            %   obj = Assembly(space) creates an assembly object with the
            %   specified mesh space. The constructor validates that the
            %   provided space is a valid MeshSpace object.
            %
            % Inputs:
            %   space - MeshSpace object
            %
            % Outputs:
            %   obj - Constructed Assembly object
            
            obj.space = space;
            obj.checkValidSpace();
        end
    end

    methods (Access = protected)
        function obj = checkValidSpace(obj)
            % CHECKVALIDSPACE Validate that the space is a MeshSpace
            % object.
            %
            %   obj = checkValidSpace(obj) verifies that the space property
            %   is a valid MeshSpace object. Throws an error if the space
            %   is not of the expected type.
            %
            % Outputs:
            %   obj - The Assembly object (for method chaining)
            
            tf = isa(obj.space, 'approx.space.MeshSpace');
            core.except.assert(tf, 'InvalidInput', ...
                'Assembly requires a MeshSpace object.');
        end
    end
end