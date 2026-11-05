classdef DgElement < approx.element.C0Element
    % DGELEMENT Element-level DG function.
    %
    %   DgElement implements an element-level function that is generally
    %   discontinuous across element boundaries and supports flux
    %   integration.
    %
    %   The class manages both volume and flux integration data, making it
    %   suitable for discontinuous finite element methods where
    %   inter-element communication occurs through numerical fluxes
    %   computed at element boundaries. Each element boundary has its own
    %   integration data for efficient flux computations.
    %
    % See also:
    %   approx.element.C0Element, approx.element.Element,
    %   approx.element.ElementFunction

    properties
        fluxData % Array of function data for flux integrals (ElementFunction array)
    end

    properties (Dependent)
        nFluxes % Number of fluxes (integer)
    end

    methods
        function obj = DgElement(geometry, projector, volumeIntegrator, fluxIntegrators)
            % DGELEMENT Constructor for DgElement.
            %
            %   obj = DgElement(geometry, projector, volumeIntegrator,
            %   fluxIntegrators) creates a DG element with both volume and
            %   flux integration capabilities using the specified
            %   components.
            %
            % Inputs:
            %   geometry - Geometry object
            %   projector - GalerkinProjector object 
            %   volumeIntegrator - Integrator object
            %   fluxIntegrators - Array of Integrator objects
            %
            % Outputs:
            %   obj - Constructed DgElement object

            obj@approx.element.C0Element(geometry, projector, volumeIntegrator);
            obj.fluxData = arrayfun( ...
                @(integrator) approx.element.ElementFunction(integrator), ...
                fluxIntegrators);
        end

        function n = get.nFluxes(obj)
            % GET.NFLUXES Get the number of fluxes.

            n = obj.geometry.nBoundaries;
        end

        function obj = setFluxData(obj, k)
            % SETFLUXDATA Set function data for flux integrals.
            %
            %   obj = setFluxData(obj, k) computes and stores function
            %   values and derivatives up to order @a k at boundary
            %   integration points for all element faces. This method
            %   prepares the data needed for flux computations in
            %   discontinuous Galerkin methods.
            %
            % Inputs:
            %   obj - The DgElement object
            %   k - Maximum derivative order (non-negative integer)
            %
            % Outputs:
            %   obj - The DgElement object

            D = obj.fluxData;

            for i = 1:obj.nFluxes
                X = D(i).integrator.nodes(1:obj.nDims, :);
                D(i).setValues(obj.projector.basis, X);
                if k > 0, D(i).setDerivatives(obj.projector.basis, X, k); end
            end
        end
    end
end