classdef Assembly < handle
    % ASSEMBLY Base class for finite element assembly operations.
    %
    %   Assembly provides the foundation for assembling discrete operators
    %   from element-level contributions in finite element methods. It
    %   encapsulates the common interface and data structures needed to
    %   transform local element operators into global system matrices.
    %
    %   The assembly process involves iterating over mesh elements,
    %   extracting local operators from element data, and accumulating
    %   contributions into global sparse matrices using appropriate
    %   connectivity information from the mesh space.
    %
    % See also:
    %   approx.assembly.VolumeAssembly, approx.assembly.FluxAssembly,
    %   approx.space.MeshSpace, approx.element.ElementOperator
    
    properties
        space % Mesh space providing geometry and connectivity (MeshSpace)
        operator % Element operator providing local computations (ElementOperator)
    end
    
    methods
        function obj = Assembly(space, operator)
            % ASSEMBLY Constructor for Assembly.
            %
            %   obj = Assembly(space, operator) creates an assembly object
            %   with the specified mesh space and element operator. The
            %   mesh space provides geometric information and global
            %   connectivity, while the element operator provides local
            %   computational kernels.
            %
            % Inputs:
            %   space - MeshSpace object providing geometry and connectivity
            %   operator - ElementOperator object providing local operators
            %
            % Outputs:
            %   obj - Constructed Assembly object

            core.except.assert(isa(space, 'approx.space.MeshSpace'), ...
                'InvalidInput', 'Assembly requires a mesh space.');
            obj.space = space;
            obj.operator = operator;
        end
    end
end