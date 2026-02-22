classdef OFL2Element < approx.element.Element
    % OFL2Element Oscillation-free L2 Galerkin element.
    %
    %   OFL2Element extends L2Element by adding the flux derivatives
    %   evaluation capabilities for oscillation-free treatment.
    %
    % See also:
    %   approx.integrate.IntegrationData, approx.element.Element

    properties (Constant)
        TYPE = 'ofl2' % Element type identifier
    end

    properties
        volume % Volume integration data (IntegrationData)
        flux % Flux integration data (IntegrationData array)
    end

    properties (Dependent)
        nFluxes % Number of fluxes (integer)
    end

    methods
        function obj = OFL2Element(G, P, VI, FI, m)
            % OFL2Element Constructor for OFL2Element.
            %
            %   obj = OFL2Element(G, P, VI, FI) creates a OFL2 element with
            %   both volume and flux integration capabilities and flux
            %   derivative order @a m.
            %
            % Inputs:
            %   G - Geometry object
            %   P - GalerkinProjector object 
            %   VI - Integrator object
            %   FI - Array of Integrator objects
            %   m - Flux derivative order
            %
            % Outputs:
            %   obj - Constructed OFL2Element object

            obj@approx.element.Element(G, P);
            obj.volume = approx.integrate.IntegrationData(VI, 1, P.basis);
            obj.flux = arrayfun(@(I) approx.integrate.IntegrationData(I, m, P.basis), FI);
        end

        function n = get.nFluxes(obj)
            % GET.NFLUXES Get the number of fluxes.

            n = length(obj.flux);
        end
    end
end