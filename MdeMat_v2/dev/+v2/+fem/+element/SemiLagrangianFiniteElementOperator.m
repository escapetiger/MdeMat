classdef SemiLagrangianFiniteElementOperator < fem.element.FiniteElementOperator
    % SEMILAGRANGIANFINITEELEMENTOPERATOR Base class for semi-Lagrangian
    % finite element operators.
    %
    %   SemiLagrangianFiniteElementOperator provides the foundation for
    %   operators used in semi-Lagrangian finite element methods. These
    %   methods track fluid parcels backward in time along characteristic
    %   curves, requiring specialized operators that handle the coupling
    %   between Eulerian and upstream (Lagrangian) data.
    %
    %   The class maintains a reference to a clipper object that manages
    %   the geometric intersection and integration over the overlap
    %   regions between Eulerian elements and upstream footprints.
    %   This is essential for accurate transport calculations in
    %   convection-dominated problems.
    %
    %   Semi-Lagrangian methods are particularly effective for problems
    %   with high Péclet numbers where traditional Eulerian methods
    %   suffer from numerical diffusion or instability.
    %
    % See also:
    %   fem.element.FiniteElementOperator,
    %   fem.element.C0SemiLagrangianFiniteElementOperator,
    %   fem.element.L2SemiLagrangianFiniteElementOperator

    properties
        clipper % Finite element clipper for geometric intersections
    end

    methods
        function obj = SemiLagrangianFiniteElementOperator(fe, clipper)
            % SEMILAGRANGIANFINITEELEMENTOPERATOR Constructor for
            % SemiLagrangianFiniteElementOperator.
            %
            %   obj = SemiLagrangianFiniteElementOperator(fe, clipper)
            %   creates a semi-Lagrangian finite element operator with
            %   the specified finite element and clipper. This constructor
            %   is called by concrete subclasses to initialize the base
            %   class functionality.
            %
            % Inputs:
            %   fe - Finite element object (FiniteElement)
            %   clipper - Geometric clipper for element intersections
            %
            % Outputs:
            %   obj - Constructed SemiLagrangianFiniteElementOperator object
            
            obj@fem.element.FiniteElementOperator(fe);
            obj.clipper = clipper;
        end
    end
end