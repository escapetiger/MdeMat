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
        Space % Mesh-based space
    end
    
    methods
        function obj = Assembly(space)
            % ASSEMBLY Constructor for Assembly base class.
            %
            %   obj = Assembly(space) creates an assembly object with the
            %   specified mesh space. The constructor validates that the
            %   provided space is a valid MeshSpace object.
            
            arguments
                space approx.space.MeshSpace
            end
            
            obj.Space = space;
        end
    end
end