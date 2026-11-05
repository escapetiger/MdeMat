classdef C0Element < approx.element.Element
    % C0ELEMENT Continuous Galerkin element.
    %
    %   C0Element implements a continuous element function. This
    %   class extends the base Element class by adding volume integration
    %   capabilities through precomputed basis values at integration
    %   points.
    %
    % Examples:
    %   % Create C0 element with quadrature
    %   geometry = core.geometry.Orthotope.unit(2);
    %   basis = approx.basis.LagrangeBasis(2);
    %   projector = approx.project.NodalProjector(basis);
    %   integrator = approx.integrate.GaussLegendre(3);
    %   element = C0Element(geometry, projector, integrator);
    %
    %   % Access pre-computed data
    %   values = element.volume.values;
    %   weights = element.volume.weights;
    %
    % See also:
    %   approx.integrate.IntegrationData, approx.element.Element

    properties (Constant)
        TYPE = 'c0' % Element type identifier
    end

    properties
        volume % Volume integration data (IntegrationData)
    end

    methods
        function obj = C0Element(geometry, projector, integrator)
            % C0ELEMENT Constructor for C0Element.
            %
            %   obj = C0Element(geometry, projector, integrator)
            %   creates a continuous element with the specified geometry,
            %   projector, and integrator for volume computations.
            %
            % Inputs:
            %   geometry - Geometry object
            %   projector - GalerkinProjector object
            %   integrator - Integrator object
            %
            % Outputs:
            %   obj - Constructed C0Element object

            obj@approx.element.Element(geometry, projector);
            obj.volume = approx.integrate.IntegrationData(integrator, 0, projector.basis);
        end
    end
end