classdef L2Element < approx.element.Element
    % L2Element L2 Galerkin element.
    %
    %   L2Element implements a L2 function. This class extends the base
    %   Element class by adding volume and flux integration capabilities
    %   through precomputed basis values and derivatives at integration
    %   points. 
    % 
    %   L2Element assumes discontinuity occurs only at element boundaries.
    %
    % See also:
    %   approx.integrate.IntegrationData, approx.element.Element

    properties (Constant)
        TYPE = 'l2' % Element type identifier
    end

    properties
        volume % Volume integration data (IntegrationData)
        flux % Flux integration data (IntegrationData array)
    end

    properties (Dependent)
        nFluxes % Number of fluxes (integer)
    end

    methods
        function obj = L2Element(G, P, VI, FI)
            % L2Element Constructor for L2Element.
            %
            %   obj = L2Element(G, P, VI, FI) creates a L2 element with
            %   both volume and flux integration capabilities using the
            %   specified components.
            %
            % Inputs:
            %   G - Geometry object
            %   P - GalerkinProjector object 
            %   VI - Integrator object
            %   FI - Array of Integrator objects
            %
            % Outputs:
            %   obj - Constructed L2Element object

            obj@approx.element.Element(G, P);
            obj.volume = approx.integrate.IntegrationData(VI, 1, P.basis);
            obj.flux = arrayfun(@(I) approx.integrate.IntegrationData(I, 0, P.basis), FI);
        end

        function n = get.nFluxes(obj)
            % GET.NFLUXES Get the number of fluxes.

            n = length(obj.flux);
        end
    end
end