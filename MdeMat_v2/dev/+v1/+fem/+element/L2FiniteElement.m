classdef L2FiniteElement < fem.element.C0FiniteElement
    % L2FINITEELEMENT Discontinuous finite element.
    %
    %   L2FiniteElement implements a discontinuous (L2) finite element that
    %   extends C0 elements with flux contributions at element boundaries.
    %   Functions in this space may be discontinuous across element
    %   boundaries, making them suitable for discontinuous Galerkin (DG)
    %   methods and problems with discontinuous solutions.
    %
    %   The L2 space allows for complete flexibility in function values
    %   across element interfaces, requiring numerical fluxes to couple
    %   neighboring elements. This approach is particularly effective for
    %   hyperbolic problems, convection-dominated flows, and problems with
    %   solution discontinuities.
    %
    %   Flux integration is performed separately on each element boundary
    %   face, with dedicated integrators for accurate computation of
    %   interface terms in the weak formulation.
    %
    % Examples:
    %   % Create L2 finite element with volume and flux integration
    %   geometry = core.geometry.Orthotope.unit(2);
    %   basis = approx.basis.LagrangeBasis(2);
    %   projector = approx.project.NodalProjector(basis);
    %   volumeIntegrator = approx.integrate.GaussLegendre(3);
    %   fluxIntegrators = arrayfun(@(i) approx.integrate.FaceIntegrator(i), 1:4);
    %   element = L2FiniteElement(geometry, projector, volumeIntegrator, fluxIntegrators);
    %
    %   % Set up data for DG method assembly
    %   element.setVolumeData(1);  % Volume terms with derivatives
    %   element.setFluxData(0);    % Flux terms with function values
    %
    %   % Access flux data for boundary face i
    %   faceValues = element.fluxData(i).values;
    %   faceWeights = element.fluxData(i).weights;
    %
    %   % Use static factory method for standard orthotope elements
    %   element = L2FiniteElement.Orthotope(2, 3, 2, 1);
    %
    % See also:
    %   fem.element.C0FiniteElement, fem.element.FiniteElement,
    %   fem.data.FunctionData

    properties
        fluxData % Array of function data for flux integrals (FunctionData array)
    end

    properties (Dependent)
        nFluxes % Number of flux boundaries (dependent property)
    end

    methods
        function obj = L2FiniteElement(K, P, I1, I2)
            % L2FINITEELEMENT Constructor for L2FiniteElement.
            %
            %   obj = L2FiniteElement(K, P, I1, I2) creates a discontinuous
            %   finite element with both volume and flux integration
            %   capabilities using the specified components.
            %
            % Inputs:
            %   K - Element geometry object (Geometry)
            %   P - Function space projector (Projector)
            %   I1 - Volume integrator object (Integrator)
            %   I2 - Array of flux integrator objects (Integrator array)
            %
            % Outputs:
            %   obj - Constructed L2FiniteElement object

            obj@fem.element.C0FiniteElement(K, P, I1);
            obj.fluxData = arrayfun(@(I) fem.data.FunctionData(I), I2);
        end

        function n = get.nFluxes(obj)
            % GET.NFLUXES Get the number of flux boundaries.
            %
            %   n = get.nFluxes(obj) returns the number of boundary faces
            %   of the element geometry, which equals the number of flux
            %   integrators required for discontinuous Galerkin methods.
            %
            % Inputs:
            %   obj - The L2FiniteElement object
            %
            % Outputs:
            %   n - Number of element boundary faces (positive integer)

            n = obj.geometry.nBoundaries;
        end

        function setFluxData(obj, r, varargin)
            % SETFLUXDATA Set function data for flux integrals.
            %
            %   setFluxData(obj, r) computes and stores function values
            %   and derivatives up to order r at boundary integration points
            %   for all element faces. This enables efficient computation
            %   of flux terms in discontinuous Galerkin formulations.
            %
            %   setFluxData(obj, r, integrationArgs) allows additional
            %   arguments to be passed to the integration point setup.
            %
            % Inputs:
            %   obj - The L2FiniteElement object
            %   r - Maximum derivative order to compute (non-negative integer)
            %   varargin - Optional arguments for integration point setup
            %
            % Outputs:
            %   NULL

            m = obj.nFluxes;
            
            %< Process each boundary face
            for i = 1:m
                D = obj.fluxData(i);
                
                %< Set integration points if additional arguments provided
                if ~isempty(varargin)
                    D.setPoints(varargin{:});
                end

                %< Extract integration and projection information
                I = D.integrator;
                n = obj.nDims;
                P = obj.projector;
                f = P.basis;
                X = I.nodes(1:n, :);

                %< Compute function values at boundary integration points
                D.setValues(f, X);
                
                %< Compute derivatives if requested
                if r > 0
                    D.setDerivatives(f, X, r);
                end
            end
        end
    end
end