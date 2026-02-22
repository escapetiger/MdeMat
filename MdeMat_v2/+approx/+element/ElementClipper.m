classdef ElementClipper < handle
    % ELEMENTCLIPPER Base class for all element clippers.
    %
    %   ElementClipper provides the foundation for geometric clipping
    %   operations in semi-Lagrangian finite element methods. Clippers
    %   handle the complex geometric intersections that arise when tracking
    %   fluid parcels backward along characteristic curves, computing the
    %   overlap between Eulerian elements and upstream footprints.
    %
    %   The clipping process is essential for accurate semi-Lagrangian
    %   transport calculations, as it determines the integration domains
    %   for coupling Eulerian and Lagrangian data. This involves computing
    %   intersections between the current element and the region occupied
    %   by fluid parcels at the previous time step.
    %
    %   Clippers organize the resulting geometric pieces and set up
    %   appropriate quadrature rules for accurate integration over
    %   these potentially complex domains. The geometric decomposition
    %   enables efficient and accurate computation of transport operators
    %   while handling the complex geometries that arise from characteristic
    %   tracking in convection-dominated flows.
    %
    % See also:
    %   approx.element.C0ElementClipper, approx.element.L2ElementClipper,
    %   approx.element.SemiLagrangianElementOperator

    properties
        element % Element object for which clipping is performed
        velocity % Characteristic velocity (column vector)
    end

    methods
        function obj = ElementClipper(element, velocity)
            % ELEMENTCLIPPER Constructor for ElementClipper.
            %
            %   obj = ElementClipper(element, velocity) creates an element
            %   clipper with the specified element and velocity field.
            %   The velocity defines the characteristic directions for
            %   backward tracing in semi-Lagrangian methods.
            %
            % Inputs:
            %   element - Element object for clipping operations (Element)
            %   velocity - Characteristic velocity (column vector, nDims x 1)
            %
            % Outputs:
            %   obj - Constructed ElementClipper object

            obj.element = element;
            obj.velocity = velocity(:);
        end
    end
end