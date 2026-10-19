classdef C0FiniteElement < fem.element.FiniteElement
    % C0FINITEELEMENT Continuous finite element.
    %
    %   C0FiniteElement implements a continuous (C0) finite element with
    %   volume integration capabilities. Functions in this space are
    %   continuous across element boundaries but may have discontinuous
    %   derivatives. This corresponds to the standard Galerkin finite
    %   element method for elliptic and parabolic problems.
    %
    %   The C0 continuity requirement ensures that the global finite element
    %   function is continuous across inter-element boundaries, which is
    %   essential for conforming finite element approximations of problems
    %   requiring H^1 regularity.
    %
    %   Volume integration is performed using numerical quadrature over the
    %   element domain, with function data pre-computed at integration points
    %   for efficient assembly of system matrices and vectors.
    %
    % Examples:
    %   % Create C0 finite element with quadrature
    %   geometry = core.geometry.Orthotope.unit(2);
    %   basis = approx.basis.LagrangeBasis(2);
    %   projector = approx.project.NodalProjector(basis);
    %   integrator = approx.integrate.GaussLegendre(3);
    %   element = C0FiniteElement(geometry, projector, integrator);
    %   
    %   % Set up function data for mass matrix assembly
    %   element.setVolumeData(0);  % Function values only
    %   
    %   % Set up function data for stiffness matrix assembly
    %   element.setVolumeData(1);  % Include first derivatives
    %
    %   % Access pre-computed data
    %   values = element.volumeData.values;
    %   weights = element.volumeData.weights;
    %
    % See also:
    %   fem.element.FiniteElement, fem.element.L2FiniteElement,
    %   fem.element.FunctionData

    properties
        volumeData % Function data for volume integrals (FunctionData object)
    end

    methods
        function obj = C0FiniteElement(K, P, I)
            % C0FINITEELEMENT Constructor for C0FiniteElement.
            %
            %   obj = C0FiniteElement(K, P, I) creates a continuous finite
            %   element with volume integration capabilities using the
            %   specified geometry, projector, and integrator.
            %
            % Inputs:
            %   K - Element geometry object (Geometry)
            %   P - Function space projector (Projector)
            %   I - Volume integrator object (Integrator)
            %
            % Outputs:
            %   obj - Constructed C0FiniteElement object

            obj@fem.element.FiniteElement(K, P);
            obj.volumeData = fem.element.FunctionData(I);
        end
                
        function obj = setVolumeData(obj, r, varargin)
            % SETVOLUMEDATA Set function data for volume integrals.
            %
            %   obj = setVolumeData(obj, r) computes and stores function
            %   values and derivatives up to order r at integration points
            %   for use in volume integral computations. This
            %   pre-computation enables efficient assembly of bilinear
            %   forms.
            %
            %   obj = setVolumeData(obj, r, args) allows additional
            %   arguments to be passed to the integration point setup.
            %
            % Inputs:
            %   obj - The C0FiniteElement object
            %   r - Maximum derivative order to compute (non-negative integer)
            %   varargin - Optional arguments for integration point setup
            %
            % Outputs:
            %   obj - The C0FiniteElement object

            D = obj.volumeData;
            
            %< Set integration points if additional arguments provided
            if ~isempty(varargin)
                D.setPoints(varargin{:});
            end

            %< Extract integration and projection information
            I = D.integrator;
            d = obj.nDims;
            P = obj.projector;
            f = P.basis;
            X = I.nodes(1:d, :);

            %< Compute function values at integration points
            D.setValues(f, X);
            
            %< Compute derivatives if requested
            if r > 0
                D.setDerivatives(f, X, r);
            end
        end
    end
end