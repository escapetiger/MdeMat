classdef L2OrthotopeElement < approx.element.L2Element
    % L2ORTHOTOPEELEMENT \f$L^2\f$ element over the orthotope.
    %
    %   L2OrthotopeElement implements an \f$L^2\f$ finite element over the
    %   orthotope domain. The element supports both modal and nodal basis
    %   functions with various polynomial orders and ordering patterns.
    %
    % See also:
    %   approx.element.L2Element

    methods
        function obj = L2OrthotopeElement(nDims, basisOrder, basisType, basisPattern)
            % L2ORTHOTOPEELEMENT Create a \f$L^2\f$ element on the unit
            % orthotope.
            %
            %   obj = L2OrthotopeElement(nDims, basisOrder, basisType,
            %   basisPattern) creates a \f$L^2\f$ element on the unit
            %   orthotope with specified basis configuration.
            %
            % Inputs:
            %   nDims - Number of dimensions (positive integer)
            %   basisOrder - Basis order (positive integer)
            %   basisType - Basis function type (char)
            %   basisPattern - Basis ordering pattern (char, 'P' or 'Q')
            %
            % Outputs:
            %   obj - The L2OrthotopeElement object
            %
            % Notes:
            %   Basis types:
            %     'legendre': Monic Legendre basis
            %     'lagrange_gauss_legendre': Lagrange-Gauss-Legendre basis
            %     'lagrange_gauss_lobatto': Lagrange-Gauss-Lobatto basis
            %
            %   Basis patterns:
            %     'Q': L_infinity norm ordering (tensor product)
            %     'P': L_1 norm ordering (total degree)

            %< Validate compatibility between basis type and pattern
            validateComptiableBasis(basisType, basisPattern);

            %< Create unit orthotope geometry
            G = core.geometry.Orthotope.unit(nDims);

            %< Create basis functions based on type
            switch basisType
                case 'legendre'
                    U = approx.basis.CompiledOrthogonalBasis( ...
                        basisOrder, 'monic_legendre', 'unit');
                case 'lagrange_gauss_legendre'
                    U = approx.basis.CompiledInterpolationBasis( ...
                        basisOrder, 'lagrange', 'unit', 'gauss_legendre');
                case 'lagrange_gauss_lobatto'
                    U = approx.basis.CompiledInterpolationBasis( ...
                        basisOrder, 'lagrange', 'unit', 'gauss_lobatto');
            end
            U = repmat(U, 1, nDims);

            %< Create separable basis with specified pattern
            switch basisPattern
                case 'Q'
                    B = approx.basis.CompiledSeparableBasis(U, 'lx');
                    B.autoLoad();
                case 'P'
                    B = approx.basis.CompiledSeparableBasis(U, 'l1');
                    B.autoLoad();
            end

            %< Create volume integrator
            switch basisType
                case {'legendre', 'lagrange_gauss_legendre'}
                    Q = 'gauss_legendre';
                case 'lagrange_gauss_lobatto'
                    Q = 'gauss_lobatto';
            end
            VI = approx.integrate.OrthotopeIntegrator(G.nDims, Q);
            VI.setPoints(B.nFactorCodims, G.lower, G.upper);

            %< Create face integrators for all boundaries
            FI = arrayfun(@(i) approx.integrate.OrthotopeFaceIntegrator( ...
                G.nDims, Q, i), 1:(2 * G.nDims));
            for i = 1:(2 * G.nDims)
                FI(i).setPoints(B.nFactorCodims, G.lower, G.upper);
            end

            %< Create approximator based on basis type
            switch basisType
                case 'legendre'
                    A = approx.linear.ModalApproximator(B);
                    A.setMass(VI.nodes, VI.weights);
                    A.setDesign(VI.nodes, VI.weights);
                case {'lagrange_gauss_legendre', 'lagrange_gauss_lobatto'}
                    A = approx.linear.NodalApproximator(B);
                    A.setMass(VI.nodes, VI.weights);
                    A.setDesign(VI.nodes);
            end

            %< Call parent constructor
            obj@approx.element.L2Element(G, A, VI, FI);
        end
    end
end

function tf = validateComptiableBasis(basisType, basisPattern)
% CHECKCOMPATIBLEBASIS Check if basis configuration is compatible.
%
%   tf = validateComptiableBasis(basisType, basisPattern) returns true if
%   the basis type and the pattern are compatiable.
%
% Inputs:
%   basisType - Type of basis functions (char)
%   basisPattern - Basis ordering pattern (char)
%
% Outputs:
%   tf - True if basis configuration is valid (logical)

tf1 = isCompatibleModalBasis(basisType, basisPattern);
tf2 = isCompatibleNodalBasis(basisType, basisPattern);
tf = tf1 | tf2;
core.except.assert(tf, 'InvalidInput', 'Invalid basis configuration.');
end

function tf = isCompatibleModalBasis(basisType, basisPattern)
% ISCOMPATIBLEMODALBASIS Check if modal basis configuration is compatible.
%
%   tf = isCompatibleModalBasis(basisType, basisPattern) returns true
%   if the basis type is modal and the pattern is supported.
%
% Inputs:
%   basisType - Type of basis functions (char)
%   basisPattern - Basis ordering pattern (char)
%
% Outputs:
%   tf - True if modal configuration is valid (logical)

tf1 = ismember(basisType, {'legendre'});
tf2 = ismember(basisPattern, {'P', 'Q'});
tf = tf1 & tf2;
end

function tf = isCompatibleNodalBasis(basisType, basisPattern)
% ISCOMPATIBLENODALBASIS Check if nodal basis configuration is compatible.
%
%   tf = isCompatibleNodalBasis(basisType, basisPattern) returns true
%   if the basis type is nodal and the pattern is supported.
%
% Inputs:
%   basisType - Type of basis functions (char)
%   basisPattern - Basis ordering pattern (char)
%
% Outputs:
%   tf - True if nodal configuration is valid (logical)

tf1 = ismember(basisType, {'lagrange_gauss_legendre', 'lagrange_gauss_lobatto'});
tf2 = ismember(basisPattern, {'Q'});
tf = tf1 & tf2;
end