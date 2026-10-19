classdef L2FiniteElementOperator < fem.element.C0FiniteElementOperator
    % L2FINITEELEMENTOPERATOR Operator for discontinuous finite elements.
    %
    %   L2FiniteElementOperator extends C0 operators with flux contributions
    %   for L2 finite elements. It manages both volume and flux bilinear
    %   form data required for discontinuous Galerkin (DG) methods, where
    %   inter-element coupling is handled through numerical fluxes.
    %
    %   The class handles the assembly of discrete differential operators
    %   that include both volume terms (from element interiors) and flux
    %   terms (from element boundaries). This dual structure is essential
    %   for DG methods where function discontinuities across element
    %   interfaces require special treatment.
    %
    %   Flux data is organized in a tensor structure to accommodate
    %   different flow directions (inflow/outflow) and derivative orders
    %   across multiple element boundaries.
    %
    % Examples:
    %   % Create L2 element and its operator
    %   element = fem.element.L2FiniteElement(geometry, projector, volInt, fluxInt);
    %   element.setVolumeData(1);  % Include first derivatives
    %   element.setFluxData(0);    % Function values on boundaries
    %   operator = L2FiniteElementOperator(element);
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
    % See also:
    %   fem.element.C0FiniteElementOperator, fem.element.L2FiniteElement,
    %   fem.element.BilinearFormData
    
    properties
        fluxData % 3D tensor of bilinear form data for flux integrals
                 %< Dimensions: (direction, derivative_order, boundary_index)
                 %< direction: 1=inflow, 2=outflow
                 %< derivative_order: 0=values, 1=first_deriv, etc.
                 %< boundary_index: element face number
    end

    methods
        function obj = L2FiniteElementOperator(fe)
            % L2FINITEELEMENTOPERATOR Constructor for L2FiniteElementOperator.
            %
            %   obj = L2FiniteElementOperator(fe) creates a discontinuous
            %   finite element operator for the specified L2 finite element.
            %   The operator provides access to discrete differential
            %   operators with both volume and flux contributions.
            %
            % Inputs:
            %   fe - L2 finite element object (L2FiniteElement)
            %
            % Outputs:
            %   obj - Constructed L2FiniteElementOperator object
            
            cls = 'fem.element.L2FiniteElement';
            core.except.assert(isa(fe, cls), ...
                'InvalidInput', 'fe must be a L2 FiniteElement.');

            obj@fem.element.C0FiniteElementOperator(fe);
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
            % Inputs:
            %   obj - The L2FiniteElementOperator object
            %
            % Outputs:
            %   obj - The L2FiniteElementOperator object
            
            D = obj.fe.fluxData;
            nf = obj.fe.nFluxes;
            nd = 2;  %< Number of directions (inflow/outflow)
            nr = max([D.nDerivatives]) + 1;  %< Maximum derivative order + 1

            %< Initialize flux data tensor
            obj.fluxData = arrayfun(@(i) fem.element.BilinearFormData(), 1:nd*nr*nf);
            obj.fluxData = reshape(obj.fluxData, nd, nr, nf);
            
            %< Populate flux data for each boundary
            for i = 1:nf
                nr = D(i).nDerivatives;
                k = 2 * ceil(i/2) - 1 + mod(i, 2);  %< Neighboring element index
                
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
        function G = buildGradient(obj)
            % BUILDGRADIENT Build gradient operator with flux contributions.
            %
            %   G = buildGradient(obj) constructs gradient operators that
            %   include both volume and flux contributions for
            %   discontinuous Galerkin methods. The flux terms handle
            %   discontinuities across element boundaries using numerical
            %   fluxes.
            %
            % Inputs:
            %   obj - The L2FiniteElementOperator object
            %
            % Outputs:
            %   G - Structure with volume and flux gradient data
            %       .volumeData - Volume gradient operators (array)
            %       .fluxData - Flux gradient operators (matrix)
            
            fe = obj.fe;
            d = fe.nDims;
            
            G = struct('volumeData', obj.volumeData(1:d), ...
                       'fluxData', squeeze(obj.fluxData(:, 1, :)));
        end
    end
end