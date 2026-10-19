classdef OFL2Element < approx.element.Element
    % OFL2ELEMENT Oscillation-free \f$L^2\f$ element.
    %
    %   OFL2Element extends the standard L2Element by adding flux
    %   derivative evaluation capabilities for oscillation-free treatment
    %   in Galerkin methods. This element type is particularly useful for
    %   problems where standard Galerkin methods may produce spurious
    %   oscillations, such as advection-dominated transport problems.
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
        function obj = OFL2Element(G, A, VI, FI, m)
            % OFL2Element Constructor for OFL2Element.
            %
            %   obj = OFL2Element(G, A, VI, FI) creates an OFL2Element with
            %   both volume and flux integration capabilities and flux
            %   derivative order @a m.
            %
            % Inputs:
            %   G - Geometry object
            %   A - LinearApproximator object 
            %   VI - Integrator object
            %   FI - Array of Integrator objects
            %   m - Flux derivative order
            %
            % Outputs:
            %   obj - Constructed OFL2Element object

            obj@approx.element.Element(G, A);
            obj.volume = approx.integrate.IntegrationData(VI, 1, A.basis);
            obj.flux = arrayfun(@(I) approx.integrate.IntegrationData(I, m, A.basis), FI);
        end

        function n = get.nFluxes(obj)
            % GET.NFLUXES Get the number of fluxes.

            n = length(obj.flux);
        end
    end
end