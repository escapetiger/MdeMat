classdef DgElementOperator < approx.element.C0ElementOperator
    % DGELEMENTOPERATOR Operator for discontinuous finite elements.
    %
    %   DgElementOperator extends C0 operators with flux contributions
    %   for DG elements. It manages both volume and flux bilinear
    %   form data required for discontinuous Galerkin (DG) methods, where
    %   inter-element coupling is handled through numerical fluxes.
    %
    %   Flux data is organized in a tensor structure to accommodate
    %   different flow directions (inflow/outflow) and derivative orders
    %   across multiple element boundaries. This organization enables
    %   efficient computation of DG operators that account for the
    %   discontinuous nature of the solution space.
    %
    %   The operator constructs discrete representations of differential
    %   operators that include both volume integrals (as in continuous
    %   methods) and boundary flux integrals (specific to DG methods).
    %
    % Examples:
    %   % Create DG element and its operator
    %   element = approx.element.DgElement(geometry, projector, volInt, fluxInt);
    %   element.setVolumeData(1);  % Include first derivatives
    %   element.setFluxData(0);    % Function values on boundaries
    %   operator = DgElementOperator(element);
    %
    %   % Set up operator data and access gradient
    %   operator.setVolumeData();
    %   operator.setFluxData();
    %   gradOp = operator.gradient;
    %
    %   % Access volume and flux contributions
    %   volGrad = gradOp.volumeData;   % Volume gradient terms
    %   fluxGrad = gradOp.fluxData;    % Boundary flux terms
    %
    %   % Use in DG assembly
    %   for i = 1:element.nDims
    %       K_vol = volGrad(i).matrix;  % Volume stiffness contribution
    %   end
    %   for i = 1:element.nFluxes
    %       K_flux = fluxGrad(:, i);    % Flux contributions for face i
    %   end
    %
    % See also:
    %   approx.element.C0ElementOperator, approx.element.DgElement,
    %   approx.element.ElementBilinearForm

    properties
        fluxData % 3D tensor of bilinear form data for flux integrals
        %< Dimensions: (direction, derivative, boundary)
        %< direction: 1=inflow, 2=outflow
        %< derivative: 0=values, 1=first_deriv, etc.
        %< boundary: element face number
        %< (ElementBilinearForm array)
    end

    methods
        function obj = DgElementOperator(element)
            % DGELEMENTOPERATOR Constructor for DgElementOperator.
            %
            %   obj = DgElementOperator(element) creates an operator for
            %   the specified DG element. The element must be an DgElement
            %   with both volume and flux data available for operator
            %   construction.
            %
            % Inputs:
            %   element - DG element object (DgElement)
            %
            % Outputs:
            %   obj - Constructed DgElementOperator object

            cls = 'approx.element.DgElement';
            core.except.assert(isa(element, cls), ...
                'InvalidInput', 'Element must be a DG Element.');

            obj@approx.element.C0ElementOperator(element);
        end

        function obj = setFluxData(obj)
            % SETFLUXDATA Set up bilinear form data for flux integrals.
            %
            %   obj = setFluxData(obj) creates a tensor of bilinear form
            %   matrices for flux computations, organizing data by
            %   derivative order, flux boundary, and flow direction. This
            %   structure enables efficient assembly of DG operators with
            %   proper treatment of inflow and outflow contributions.
            %
            %   The flux data tensor accommodates:
            %   - Inflow contributions: coupling within the same element
            %   - Outflow contributions: coupling with neighboring elements
            %   - Multiple derivative orders for higher-order DG methods
            %   - All element boundaries for complete flux integration
            %
            % Inputs:
            %   obj - The DgElementOperator object
            %
            % Outputs:
            %   obj - The DgElementOperator object

            D = obj.element.fluxData;
            nf = obj.element.nFluxes;
            nd = 2; %< Number of directions (inflow/outflow)
            nr = max([D.nDerivatives]) + 1; %< Maximum derivative order + 1

            %< Initialize flux data tensor
            f = @(i) approx.element.ElementBilinearForm();
            obj.fluxData = arrayfun(f, 1:nd*nr*nf);
            obj.fluxData = reshape(obj.fluxData, nd, nr, nf);

            %< Populate flux data for each boundary
            for i = 1:nf
                nr = D(i).nDerivatives;
                k = 2 * ceil(i/2) - 1 + mod(i, 2); %< Neighboring element index

                for j = 0:nr
                    %< Inflow contribution (same element)
                    obj.fluxData(1, j+1, i).setData(D(i), D(i), j, 0);
                    %< Outflow contribution (neighboring element)
                    obj.fluxData(2, j+1, i).setData(D(i), D(k), j, 0);
                end
            end
        end
    end

    methods (Access = protected)
        function G = extractGradient(obj)
            % EXTRACTGRADIENT Extract gradient data with flux
            % contributions.
            %
            %   G = extractGradient(obj) extracts the gradient data
            %   including both volume and flux contributions. The volume
            %   data corresponds to interior gradient terms, while the
            %   flux data represents boundary gradient contributions in
            %   discontinuous Galerkin methods.
            %
            % Inputs:
            %   obj - The DgElementOperator object
            %
            % Outputs:
            %   G - Gradient data structure containing:
            %<       .volumeData - Volume gradient operators for each dimension
            %<       .fluxData - Flux gradient operators for each boundary

            G = struct( ...
                'volumeData', obj.volumeData(1:obj.element.nDims), ...
                'fluxData', squeeze(obj.fluxData(:, 1, :)));
        end
    end
end