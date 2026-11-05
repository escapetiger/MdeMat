classdef MeshSpace < handle
    % MESHSPACE Base class for all mesh-based function spaces.
    %
    %   MeshSpace represents a function space defined over a geometric
    %   mesh. It serves as the foundation of classical numerical methods
    %   leveraging the local approximation techniques, such as finite
    %   element method, finite difference method and finite volume method.

    properties
        mesh % Mesh object
    end

    properties (Dependent)
        nDims % Number of dimensions
        nMeshElements % Total number of mesh elements
    end

    methods
        function obj = MeshSpace(mesh)
            % MESHSPACE Constructor for MeshSpace.
            %
            %   obj = MeshSpace(mesh) creates a function space defined on
            %   the specified geometric mesh.
            %
            % Inputs:
            %   mesh - Mesh object
            %
            % Outputs:
            %   obj - Constructed MeshSpace object

            obj.mesh = mesh;
        end

        function n = get.nDims(obj)
            % GET.NDIMS Returns the number of dimensions.
            
            n = obj.mesh.nDims;
        end

        function n = get.nMeshElements(obj)
            % GET.NTOTALELEMENTS Get the total number of elements.

            n = obj.mesh.nElements;
        end
    
        function newObj = refine(obj, nLevels)
            % REFINE Create refined mesh space.
            %
            %   newObj = refine(obj, nLevels) creates a refined mesh space.
            %
            % Inputs:
            %   obj - The MeshSpace object
            %   nLevels - Number of refinement levels (non-negative integer)
            %
            % Outputs:
            %   newObj - Refined MeshSpace object

            newMesh = obj.mesh.refine(nLevels);
            newObj = approx.space.MeshSpace(newMesh);
        end
    end
end