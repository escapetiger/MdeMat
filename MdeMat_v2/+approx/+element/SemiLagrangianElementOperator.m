classdef SemiLagrangianElementOperator < approx.element.ElementOperator
    % SEMILAGRANGIANELEMENTOPERATOR Base class for semi-Lagrangian element
    % operators.
    %
    %   SemiLagrangianElementOperator provides the foundation for
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
    %   suffer from numerical diffusion or instability. By following
    %   characteristics backward, these methods naturally handle
    %   advection without the CFL time step restriction.
    %
    % See also:
    %   approx.element.ElementOperator,
    %   approx.element.C0SemiLagrangianElementOperator,
    %   approx.element.L2SemiLagrangianElementOperator

    properties
        clipper % Element clipper for geometric intersections
    end

    methods
        function obj = SemiLagrangianElementOperator(element, clipper)
            % SEMILAGRANGIANELEMENTOPERATOR Constructor for
            % SemiLagrangianElementOperator.
            %
            %   obj = SemiLagrangianElementOperator(element, clipper)
            %   creates a semi-Lagrangian element operator with the
            %   specified element and clipper. The clipper handles the
            %   geometric computations required for intersection and
            %   integration of Eulerian elements with upstream footprints.
            %
            % Inputs:
            %   element - Element object on which the operator acts
            %   clipper - Geometric clipper for element intersections
            %
            % Outputs:
            %   obj - Constructed SemiLagrangianElementOperator object
            
            obj@approx.element.ElementOperator(element);
            obj.clipper = clipper;
        end
    end
end