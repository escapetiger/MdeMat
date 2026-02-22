classdef L2OrthotopeFiniteElement < fem.element.L2FiniteElement
    % L2ORTHOTOPEFINITEELEMENT Discontinuous finite element defined on an
    % orthotope.
    %
    % Examples:
    %   % Create L2 finite element on an orthotope
    %   element = L2OrthotopeFiniteElement(nDims, order, id, pattern);
    %
    % See also:
    %   fem.element.L2FiniteElement

    methods
        function obj = L2OrthotopeFiniteElement(varargin)
            % L2ORTHOTOPEFINITEELEMENT Create standard L2 finite element on
            % unit orthotope.
            %
            %   obj = L2OrthotopeFiniteElement(...) creates a discontinuous
            %   finite element on the unit orthotope with specified
            %   polynomial order and basis configuration.
            %
            % Inputs:
            %   varargin - Input arguments
            %<   nDims - Number of spatial dimensions (positive integer)
            %<   order - Polynomial order (non-negative integer)
            %<   id - Basis type identifier
            %<   pattern - Tensor product pattern
            %
            % Notes:
            %   id = 1: Orthogonal (Legendre) basis
            %   id = 2: Lagrange-Gauss-Legendre basis
            %   id = 3: Lagrange-Gauss-Lobatto basis
            %   pattern = 1: \f$L_\infty\f$ norm ordering
            %   pattern = 2: \f$L_1\f$ norm ordering
            %
            % Outputs:
            %   obj - Constructed L2OrthotopeFiniteElement object

            p = inputParser;
            addParameter(p, 'nDims', [], @(x) isnumeric(x) && isscalar(x));
            addParameter(p, 'order', 1, @(x) isnumeric(x) && isscalar(x));
            addParameter(p, 'id', 1, @(x) isnumeric(x) && isscalar(x));
            addParameter(p, 'pattern', 1, @(x) isnumeric(x));
            parse(p, varargin{:});

            %< Create unit orthotope geometry
            K = core.geometry.Orthotope.unit(p.Results.nDims);

            %< Create basis functions based on type
            switch p.Results.id
                case 1
                    U = approx.basis.OrthogonalBasisFunction( ...
                        p.Results.order, 'monic_legendre', 'unit');
                case 2
                    U = approx.basis.InterpolationBasisFunction( ...
                        p.Results.order, 'lagrange', 'unit', 'gauss_legendre');
                case 3
                    U = approx.basis.InterpolationBasisFunction( ...
                        p.Results.order, 'lagrange', 'unit', 'gauss_lobatto');
            end
            U = repmat(U, 1, p.Results.nDims);

            %< Create separable basis with specified pattern
            switch p.Results.pattern
                case 1
                    B = approx.basis.SeparableBasisFunction(U, 'lx');
                case 2
                    B = approx.basis.SeparableBasisFunction(U, 'l1');
            end
            B.autoLoad();

            %< Create volume integrator
            switch p.Results.id
                case {1, 2}
                    Q = 'gauss_legendre';
                case 3
                    Q = 'gauss_lobatto';
            end
            I1 = approx.integrate.OrthotopeIntegrator(K.nDims, Q);
            I1.setPoints(B.nFactorCodims, false, K.lower, K.upper);

            %< Create face integrators for all boundaries
            I2 = arrayfun(@(i) approx.integrate.OrthotopeFaceIntegrator( ...
                K.nDims, Q, i), 1:(2 * K.nDims));
            for i = 1:(2 * K.nDims)
                I2(i).setPoints(B.nFactorCodims, false, K.lower, K.upper);
            end

            %< Create projector based on basis type
            switch p.Results.id
                case 1
                    P = approx.project.ModalProjector(B);
                    P.setMass(I1.nodes, I1.weights);
                case {2, 3}
                    P = approx.project.NodalProjector(B);
                    P.setMass(I1.nodes);
            end

            %< Construct L2 orthotope finite element
            obj@fem.element.L2FiniteElement(K, P, I1, I2);
        end
    end
end