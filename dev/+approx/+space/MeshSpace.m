classdef MeshSpace < handle
    % MESHSPACE Base class for all mesh-based function spaces.
    %
    %   MeshSpace represents a function space defined over a geometric
    %   mesh. It serves as the foundation of classical numerical methods
    %   leveraging the local approximation techniques, such as finite
    %   element method, finite difference method and finite volume method.
    
    properties
        Mesh % Mesh object
    end
    
    properties (Dependent)
        NDims % Number of dimensions
        NMeshElements % Total number of mesh elements
    end
    
    methods
        function obj = MeshSpace(mesh)
            % MESHSPACE Constructor for MeshSpace.
            %
            %   obj = MeshSpace(mesh) creates a function space defined on
            %   the specified geometric mesh.
            
            arguments
                mesh approx.mesh.Mesh
            end
            
            obj.Mesh = mesh;
        end
        
        function n = get.NDims(obj)
            % GET.NDIMS Returns the number of dimensions.
            
            n = obj.Mesh.NDims;
        end
        
        function n = get.NMeshElements(obj)
            % GET.NTOTALELEMENTS Get the total number of elements.
            
            n = obj.Mesh.NElements;
        end
        
        function newObj = refine(obj, nLevels)
            % REFINE Create refined mesh space.
            %
            %   newObj = refine(obj, nLevels) creates a refined mesh space.
            
            arguments
                obj approx.space.MeshSpace
                nLevels (1,1) {mustBeInteger, mustBePositive}
            end
            
            newMesh = obj.Mesh.refine(nLevels);
            newObj = approx.space.MeshSpace(newMesh);
        end
    end
end