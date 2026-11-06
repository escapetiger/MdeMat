classdef FiniteElementSpace < mdemat.common.fem.MeshSpace
    % FINITEELEMENTSPACE Base class for all finite element spaces.
    %
    %   FiniteElementSpace implements the finite-dimensional function space
    %   used in the finite element method, It serves as the foundation for
    %   finite element computations on complex geometries.
    %
    % See also:
    %   mdemat.common.fem.MeshSpace
    
    properties
        Element % Element object
    end
    
    properties (Dependent)
        NLocalDofs % Number of local degrees of freedom
        NGlobalDofs % Number of global degrees of freedom
    end
    
    methods
        function obj = FiniteElementSpace(mesh, element)
            % FINITEELEMENTSPACE Constructor for FiniteElementSpace.
            %
            %   obj = FiniteElementSpace(mesh, element) creates a finite
            %   element space defined by the specified mesh and finite
            %   element.
            
            arguments
                mesh mdemat.data.mesh.Mesh
                element mdemat.common.fem.Element
            end
            
            obj@mdemat.common.fem.MeshSpace(mesh);
            obj.Element = element;
        end
        
        function n = get.NLocalDofs(obj)
            % GET.NLOCALDOFS Get the number of degrees of freedom per
            % element.
            
            n = obj.Element.NDofs;
        end
        
        function n = get.NGlobalDofs(obj)
            % GET.NGLOBALDOFS Get the total number of degrees of freedom.
            
            n = obj.NLocalDofs * obj.NMeshElements;
        end

        function weights = integrationWeights(obj)
            % INTEGRATIONWEIGHTS Physical integral of each global basis mode.
            %
            %   The ng-by-1 column follows localToGlobal ordering. For an
            %   affine element map F_K and basis phi_j on the reference
            %   element, ell_(K,j) = |det J_K| * sum_q w_q*phi_j(xi_q).
            %   Thus weights.'*C integrates each coefficient column of C;
            %   the transpose is not conjugated, even for a complex basis.

            arguments
                obj mdemat.common.fem.FiniteElementSpace
            end

            assert(all(isa(obj.Element, 'mdemat.common.fem.L2Element'), 'all'), ...
                'mdemat:common:fem:FiniteElementSpace:UnsupportedElement', '%s', 'Physical integration requires a volume quadrature.');
            volume = obj.Element.Volume;
            values = obj.Element.Approximator.Basis.eval(volume.Nodes);
            referenceIntegrals = values * volume.Weights(:);
            determinants = abs(obj.Mesh.computeElementJacobianDeterminants());
            if isscalar(determinants)
                determinants = repmat(determinants, obj.NMeshElements, 1);
            end
            weights = kron(determinants(:), referenceIntegrals);
        end
        
        function Y = eval(obj, xRef, C, options)
            % EVAL Evaluate finite element functions in the mesh space.
            %
            %   Y = eval(obj, xRef, C) evaluates the piecewise
            %   function defined by coefficients @a C at reference points
            %   @a xRef on all elements.
            %
            %   Y = eval(obj, xRef, C, EI=EI) evaluates on specified elements @a EI.
            %   C has global or selected-element rows. coefficientLayout='auto'
            %   recognizes the row count and favors selected rows when both
            %   interpretations coincide; 'global' explicitly reorders all rows.
            
            arguments
                obj mdemat.common.fem.FiniteElementSpace
                xRef {mustBeNumeric}
                C {mustBeNumeric}
                options.EI {mustBeInteger, mustBePositive} = []
                options.coefficientLayout {mustBeMember(options.coefficientLayout, ...
                    {'auto', 'global', 'selected'})} = 'auto'
            end
            
            if isempty(options.EI)
                EI = (1:obj.NMeshElements).';
            else
                EI = options.EI(:);
            end
            
            if isempty(xRef)
                xRef = obj.Element.Volume.Nodes;
            end
            
            ne = size(EI, 1);
            np = size(xRef, 2);
            nl = obj.NLocalDofs;
            C = selectCoefficients(obj, C, EI, options.coefficientLayout);
            C = reshape(C, nl, []);
            Y = obj.Element.eval(xRef, C); %< ne*nc x np
            Y = reshape(Y.', np*ne, []);
        end
        
        function Y = diff(obj, xRef, C, order, options)
            % DIFF Differentiate finite element functions in the
            % mesh space.
            %
            %   Y = diff(obj, xRef, C, order) differentiates
            %   piecewise function defined by coefficients @a C at
            %   reference points @a xRef on all elements.
            %
            %   Y = diff(obj, xRef, C, order, EI=EI) differentiates on
            %   specified elements @a EI.
            %   Derivatives are with respect to reference coordinates xRef.
            %   Coefficient row selection follows eval's coefficientLayout.
            
            arguments
                obj mdemat.common.fem.FiniteElementSpace
                xRef {mustBeNumeric}
                C {mustBeNumeric}
                order {mustBeInteger, mustBeNonnegative}
                options.EI {mustBeInteger, mustBePositive} = []
                options.coefficientLayout {mustBeMember(options.coefficientLayout, ...
                    {'auto', 'global', 'selected'})} = 'auto'
            end
            
            if isempty(options.EI)
                EI = (1:obj.NMeshElements).';
            else
                EI = options.EI(:);
            end
            
            if isempty(xRef)
                xRef = obj.Element.Volume.Nodes;
            end
            
            ne = size(EI, 1);
            np = size(xRef, 2);
            nl = obj.NLocalDofs;
            C = selectCoefficients(obj, C, EI, options.coefficientLayout);
            C = reshape(C, nl, []);
            Y = obj.Element.diff(xRef, C, order); %< ne*nc x np
            Y = reshape(Y.', np*ne, []);
        end
        
        function Y = feval(obj, f, xRef, options)
            % FEVAL Evaluate analytical function on mesh elements.
            %
            %   Y = feval(obj, f, xRef) evaluates the analytical function
            %   @a f at reference points @a xRef mapped to physical
            %   coordinates on all elements.
            %
            %   Y = feval(obj, f, xRef, EI=EI) evaluates on specified
            %   elements @a EI.
            %
            %   Y = feval(obj, f, xRef, args=args) passes additional
            %   arguments to the function @a f.
            
            arguments
                obj mdemat.common.fem.FiniteElementSpace
                f {mustBeA(f, {'numeric', 'function_handle'})}
                xRef {mustBeNumeric}
                options.EI {mustBeInteger, mustBePositive} = []
                options.args = {}
            end
            
            if isempty(options.EI)
                EI = (1:obj.NMeshElements).';
            else
                EI = options.EI(:);
            end
            
            X = obj.Mesh.collocate(xRef, EI);
            Y = evaluateField(f, X, options.args);
        end
        
        function A = mean(obj, f, options)
            % MEAN Compute weighted mean of a function over specified
            % elements.
            %
            %   A = mean(obj, f) computes the weighted mean of function
            %   @a f using integration weights over all elements.
            %
            %   A = mean(obj, f, EI=EI) computes over specified elements @a EI.
            %
            %   A = mean(obj, f, args=args) passes additional arguments to
            %   function @a f.
            
            arguments
                obj mdemat.common.fem.FiniteElementSpace
                f {mustBeA(f, 'function_handle')}
                options.EI {mustBeInteger, mustBePositive} = []
                options.args = {}
            end
            
            if isempty(options.EI)
                EI = (1:obj.NMeshElements).';
            else
                EI = options.EI(:);
            end
            
            ne = size(EI, 1);
            np = obj.Element.Volume.NPoints;
            wRef = obj.Element.Volume.Weights;
            xRef = obj.Element.Volume.Nodes;
            X = obj.Mesh.collocate(xRef, EI);
            F = f(X, options.args{:});
            F = reshape(F, np, []);
            A = wRef * F;
            A = reshape(A, ne, []);
        end
        
        function U = fit(obj, f, options)
            % FIT Fit a function within the finite element space.
            %
            %   U = fit(obj, f) fits the function @a f onto the finite
            %   element space over all elements using the appropriate
            %   fitting method.
            %
            %   U = fit(obj, f, EI=EI) fits over specified elements @a EI.
            %
            %   U = fit(obj, f, args=args) passes additional arguments to
            %   function @a f. X=X supplies physical sample coordinates in
            %   the configured quadrature/design order for selected elements.
            
            arguments
                obj mdemat.common.fem.FiniteElementSpace
                f {mustBeA(f, {'numeric', 'function_handle'})}
                options.EI {mustBeInteger, mustBePositive} = []
                options.X = []
                options.args = {}
            end
            
            if isempty(options.EI)
                EI = (1:obj.NMeshElements).';
            else
                EI = options.EI(:);
            end
            
            ne = size(EI, 1);
            nl = obj.NLocalDofs;
            if isempty(options.X)
                xRef = obj.Element.Volume.Nodes;
                X = obj.Mesh.collocate(xRef, EI);
            else
                X = options.X;
                assert(all(isnumeric(X) && size(X, 1) == obj.NDims, 'all'), ...
                    'mdemat:common:fem:FiniteElementSpace:InvalidInput', '%s', 'Explicit points must have NDims rows.');
                assert(all(size(X, 2) == obj.Element.Volume.NPoints * ne, 'all'), ...
                    'mdemat:common:fem:FiniteElementSpace:InvalidInput', '%s', 'Explicit points must match the configured element samples.');
            end

            U = evaluateField(f, X, options.args);
            
            A = obj.Element.Approximator;
            switch A.Type
                case 'modal'
                    D = obj.Element.Volume;
                    F = A.embed(U, D.Values, D.Weights);
                case 'nodal'
                    F = A.embed(U);
                otherwise
                    assert(all(false, 'all'), 'mdemat:common:fem:FiniteElementSpace:UnsupportedProjector', ...
                        'Projector type %s not supported.', A.Type);
            end

            U = A.fit(F);
            U = reshape(U, nl*ne, []);
        end
        
        function GI = localToGlobal(obj, EI, LI)
            % LOCALTOGLOBAL Map local DoF indices to global DoF indices.
            %
            %   GI = localToGlobal(obj, EI, LI) maps local DoF indices @a LI to
            %   global DoF indices for the specified elements @a EI.
            
            arguments
                obj mdemat.common.fem.FiniteElementSpace
                EI {mustBeInteger, mustBePositive}
                LI {mustBeInteger, mustBePositive}
            end
            
            nl = obj.NLocalDofs;
            
            if isempty(EI)
                EI = (1:obj.NMeshElements);
            end
            
            if isempty(LI)
                LI = (1:nl);
            end
            assert(all(EI <= obj.NMeshElements, 'all') && all(LI <= nl, 'all'), ...
                'MdeMat:Fem:IndexRange', 'Element or local degree-of-freedom index is out of range.');
            GI = LI(:) + (EI(:).' - 1) * nl;
            GI = GI(:);
        end
        
        function [EI, LI] = globalToLocal(obj, GI)
            % GLOBALTOLOCAL Map global DoF indices to local DoF indices.
            %
            %   [EI, LI] = globalToLocal(obj, GI) maps global DoF indices @a GI
            %   to local DoF indices.
            
            arguments
                obj mdemat.common.fem.FiniteElementSpace
                GI {mustBeInteger, mustBePositive}
            end
            
            nl = obj.NLocalDofs;
            assert(all(GI <= obj.NGlobalDofs, 'all'), 'MdeMat:Fem:IndexRange', ...
                'Global degree-of-freedom index is out of range.');
            EI = ceil(GI/nl);
            LI = GI - (EI - 1) * nl;
        end
        
        function newObj = refine(obj, nLevels)
            % REFINE Create refined finite element space.
            %
            %   newObj = refine(obj, nLevels) creates a refined finite
            %   element space.
            
            arguments
                obj mdemat.common.fem.FiniteElementSpace
                nLevels (1,1) {mustBeInteger, mustBePositive}
            end
            
            newMesh = obj.Mesh.refine(nLevels);
            newObj = mdemat.common.fem.FiniteElementSpace(newMesh, obj.Element);
        end
    end
end

function coefficients = selectCoefficients(space, coefficients, elements, layout)
assert(all(elements <= space.NMeshElements), 'MdeMat:Fem:IndexRange', ...
    'Selected element index is out of range.');
selectedRows = space.NLocalDofs*numel(elements);
globalRows = space.NGlobalDofs;
if isrow(coefficients) && (selectedRows > 1 || ...
        (strcmp(layout, 'global') && globalRows > 1))
    coefficients = coefficients(:);
end
useGlobal = strcmp(layout, 'global') || (strcmp(layout, 'auto') && ...
    size(coefficients, 1) == globalRows && globalRows ~= selectedRows);
if useGlobal
    assert(size(coefficients, 1) == globalRows, 'MdeMat:Fem:CoefficientShape', ...
        'Global coefficient data requires NGlobalDofs rows.');
    coefficients = coefficients(space.localToGlobal(elements, []), :);
else
    assert(size(coefficients, 1) == selectedRows, 'MdeMat:Fem:CoefficientShape', ...
        'Selected coefficient data requires NLocalDofs times selected-element rows.');
end
end

function values = evaluateField(field, coordinates, args)
% Map the public components-by-points convention to FE points-by-fields.
% Keep recognizable older points-by-fields/packed-vector callbacks readable.
% A square matrix follows the public convention, rather than guessing axes.
nPoints = size(coordinates, 2);
if isnumeric(field)
    if isempty(field) || ~isvector(field)
        error('mdemat:fem:InvalidField', 'Numeric fields must be scalar or vector constants.');
    end
    values = repmat(field(:).', nPoints, 1);
    return;
end
values = field(coordinates, args{:});
if ~isnumeric(values) || isempty(values) || ~ismatrix(values)
    error('mdemat:fem:InvalidField', 'A field must return a numeric sample matrix.');
end
if isscalar(values)
    values = repmat(values, nPoints, 1);
elseif isvector(values) && mod(numel(values), nPoints) == 0
    values = reshape(values, nPoints, []);
elseif size(values, 2) == nPoints
    values = values.';
elseif size(values, 1) ~= nPoints
    error('mdemat:fem:InvalidFieldShape', ...
        'Field samples require one column per point (or a recognizable legacy layout).');
end
end
