classdef FiniteElementClipper < handle
    % FINITEELEMENTCLIPPER Abstract base class for finite element clippers.
    %
    %   FiniteElementClipper provides the foundation for geometric clipping
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
    %   these potentially complex domains.
    %
    % See also:
    %   fem.element.C0FiniteElementClipper, fem.element.L2FiniteElementClipper,
    %   fem.element.SemiLagrangianFiniteElementOperator

    properties
        fe % Finite element object for which clipping is performed
        velocity % Velocity vector for characteristic tracking
    end

    methods
        function obj = FiniteElementClipper(fe, velocity)
            % FINITEELEMENTCLIPPER Constructor for FiniteElementClipper.
            %
            %   obj = FiniteElementClipper(fe, velocity) creates a finite
            %   element clipper with the specified finite element and
            %   velocity field. This constructor is called by concrete
            %   subclasses to initialize the base class functionality.
            %
            % Inputs:
            %   fe - Finite element object (FiniteElement)
            %   velocity - Velocity vector for characteristic tracking (column vector)
            %
            % Outputs:
            %   obj - Constructed FiniteElementClipper object

            obj.fe = fe;
            obj.velocity = velocity(:);
        end
    end
end