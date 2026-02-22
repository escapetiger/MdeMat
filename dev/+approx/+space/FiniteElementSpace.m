classdef FiniteElementSpace < approx.space.MeshSpace
    % FINITEELEMENTSPACE Base class for all finite element spaces.
    %
    %   FiniteElementSpace implements the finite-dimensional function space
    %   used in the finite element method, It serves as the foundation for
    %   finite element computations on complex geometries.
    %
    % See also:
    %   approx.space.MeshSpace
    
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
                mesh approx.mesh.Mesh
                element approx.element.Element
            end
            
            obj@approx.space.MeshSpace(mesh);
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
        
        function Y = eval(obj, xRef, C, options)
            % EVAL Evaluate finite element functions in the mesh space.
            %
            %   Y = eval(obj, xRef, C) evaluates the piecewise
            %   function defined by coefficients @a C at reference points
            %   @a xRef on all elements.
            %
            %   Y = eval(obj, xRef, C, EI=EI) evaluates on specified elements @a EI.
            
            arguments
                obj approx.space.FiniteElementSpace
                xRef {mustBeNumeric}
                C {mustBeNumeric}
                options.EI {mustBeInteger, mustBePositive} = []
            end
            
            if isempty(options.EI)
                EI = (1:obj.NMeshElements).';
            else
                EI = options.EI;
            end
            
            if isempty(xRef)
                xRef = obj.Element.Volume.Nodes;
            end
            
            ne = size(EI, 1);
            np = size(xRef, 2);
            nl = obj.NLocalDofs;
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
            
            arguments
                obj approx.space.FiniteElementSpace
                xRef {mustBeNumeric}
                C {mustBeNumeric}
                order {mustBeInteger, mustBeNonnegative}
                options.EI {mustBeInteger, mustBePositive} = []
            end
            
            if isempty(options.EI)
                EI = (1:obj.NMeshElements).';
            else
                EI = options.EI;
            end
            
            if isempty(xRef)
                xRef = obj.Element.Volume.Nodes;
            end
            
            ne = size(EI, 1);
            np = size(xRef, 2);
            nl = obj.NLocalDofs;
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
                obj approx.space.FiniteElementSpace
                f {mustBeA(f, 'function_handle')}
                xRef {mustBeNumeric}
                options.EI {mustBeInteger, mustBePositive} = []
                options.args = {}
            end
            
            if isempty(options.EI)
                EI = (1:obj.NMeshElements).';
            else
                EI = options.EI;
            end
            
            ne = size(EI, 1);
            np = size(xRef, 2);
            X = obj.Mesh.collocate(xRef, EI);
            Y = f(X, options.args{:});
            Y = reshape(Y, ne*np, []);
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
                obj approx.space.FiniteElementSpace
                f {mustBeA(f, 'function_handle')}
                options.EI {mustBeInteger, mustBePositive} = []
                options.args = {}
            end
            
            if isempty(options.EI)
                EI = (1:obj.NMeshElements).';
            else
                EI = options.EI;
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
            %   function @a f.
            
            arguments
                obj approx.space.FiniteElementSpace
                f {mustBeA(f, 'function_handle')}
                options.EI {mustBeInteger, mustBePositive} = []
                options.args = {}
            end
            
            if isempty(options.EI)
                EI = (1:obj.NMeshElements).';
            else
                EI = options.EI;
            end
            
            A = obj.Element.Approximator;
            
            ne = size(EI, 1);
            nl = obj.NLocalDofs;
            if isempty(options.args)
                xRef = obj.Element.Volume.Nodes;
                X = obj.Mesh.collocate(xRef, EI);
                U = f(X);
            else
                U = f(options.args{:});
            end
            
            switch A.Type
                case 'modal'
                    D = obj.Element.Volume;
                    F = A.embed(U, D.Values, D.Weights);
                case 'nodal'
                    F = A.embed(U);
                otherwise
                    core.except.error('UnsupportedProjector', ...
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
                obj approx.space.FiniteElementSpace
                EI {mustBeInteger, mustBeNonnegative}
                LI {mustBeInteger, mustBePositive}
            end
            
            nl = obj.NLocalDofs;
            
            if isempty(EI)
                EI = (1:obj.NMeshElements);
            end
            
            if isempty(LI)
                LI = (1:nl);
            end
            GI = LI(:) + (EI(:).' - 1) * nl;
            GI = GI(:);
        end
        
        function [EI, LI] = globalToLocal(obj, GI)
            % GLOBALTOLOCAL Map global DoF indices to local DoF indices.
            %
            %   [EI, LI] = globalToLocal(obj, GI) maps global DoF indices @a GI
            %   to local DoF indices.
            
            arguments
                obj approx.space.FiniteElementSpace
                GI {mustBeInteger, mustBePositive}
            end
            
            nl = obj.NLocalDofs;
            EI = ceil(GI/nl);
            LI = GI - (EI - 1) * nl;
        end
        
        function newObj = refine(obj, nLevels)
            % REFINE Create refined finite element space.
            %
            %   newObj = refine(obj, nLevels) creates a refined finite
            %   element space.
            
            arguments
                obj approx.space.FiniteElementSpace
                nLevels (1,1) {mustBeInteger, mustBePositive}
            end
            
            newMesh = obj.Mesh.refine(nLevels);
            newObj = approx.space.FiniteElementSpace(newMesh, obj.Element);
        end
    end
end