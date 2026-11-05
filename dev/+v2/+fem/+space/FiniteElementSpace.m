classdef FiniteElementSpace < handle
    % FINITEELEMENTSPACE Base class for finite element function spaces.
    %
    %   FiniteElementSpace provides an abstract representation of a
    %   function space that is constructed from a finite element. This
    %   class serves as the foundation for specialized function spaces such
    %   as affine spaces, mesh spaces, and scaled variants.
    %
    %   The class encapsulates the finite element that generates the space
    %   and provides common properties and interfaces that all concrete
    %   finite element spaces must implement. Each space defines how
    %   functions are represented, evaluated, and manipulated within that
    %   space.
    %
    % See also:
    %   fem.space.AffineSpace, fem.space.MeshSpace, 
    %   fem.space.ScaledAffineSpace
    
    properties
        fe % Finite element
    end

    properties (Dependent)
        nDims % Number of spatial dimensions
    end
    
    methods (Abstract)
        % EVALUATE Evaluate functions in the finite element space.
        %
        %   Abstract method that must be implemented by concrete subclasses
        %   to define how functions in the space are evaluated at given
        %   points. The specific signature and behavior depends on the
        %   concrete space implementation.
        %
        % Inputs:
        %   obj - The FiniteElementSpace object
        %   varargin - Variable arguments depending on space type
        %
        % Outputs:
        %   Y - Evaluation results (format depends on space type)
        Y = evaluate(obj, varargin)
    end

    methods
        function obj = FiniteElementSpace(fe)
            % FINITEELEMENTSPACE Constructor for FiniteElementSpace.
            %
            %   obj = FiniteElementSpace(fe) creates an abstract finite
            %   element space using the specified finite element as the
            %   basis for function representation.
            %
            % Inputs:
            %   fe - Finite element object that defines the space basis
            %
            % Outputs:
            %   obj - Constructed FiniteElementSpace object
            
            obj.fe = fe;
        end

        function n = get.nDims(obj)
            % GET.NDIMS Get the number of spatial dimensions.
            
            n = obj.fe.nDims;
        end
    end
end