classdef BH1Element < approx.element.L2Element
    % BH1ELEMENT Broken \f$H^1\f$ element with volume and flux integration.
    %
    %   BH1Element implements a discontinuous finite element based on the
    %   broken \f$H^1\f$ function space. This class extends the L2Element class
    %   by adding flux integration capabilities through
    %   precomputed basis values and derivatives at integration points.
    %
    %   BH1Element assumes discontinuity occurs only at element boundaries,
    %   making them suitable for discontinuous Galerkin (DG) methods and
    %   finite volume schemes. The element supports both modal (orthogonal)
    %   and nodal (interpolation) basis functions.
    %
    % See also:
    %   approx.integrate.Integrator, approx.element.L2Element

    properties
        Flux % Flux integrators (Object array of Integrators)
    end

    properties (Dependent)
        NFluxes % Number of fluxes (integer)
    end

    methods
        function obj = BH1Element(geometry, approximator, volumeIntegrator, fluxIntegrators)
            % BH1ELEMENT Constructor for BH1Element.
            %
            %   obj = BH1Element(geometry, approximator, volumeIntegrator,
            %   fluxIntegrators) creates a BH1 element with both volume and
            %   flux integration capabilities using the specified components.

            arguments
                geometry core.geometry.Geometry
                approximator approx.linear.LinearApproximator
                volumeIntegrator approx.integrate.Integrator
                fluxIntegrators (1, :) {mustBeA(fluxIntegrators, 'approx.integrate.Integrator')}
            end

            obj@approx.element.L2Element(geometry, approximator, volumeIntegrator);
            obj.Volume = obj.Volume.setPartials(approximator.Basis, 1);
            obj.Flux = fluxIntegrators;
            for i = 1:length(fluxIntegrators)
                obj.Flux(i) = obj.Flux(i).setValues(approximator.Basis);
            end
        end

        function n = get.NFluxes(obj)
            % GET.NFLUXES Get the number of fluxes.

            n = length(obj.Flux);
        end
    end
end