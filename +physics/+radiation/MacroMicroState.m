classdef MacroMicroState < physics.state.KineticState
    % MACROMICROSTATE State class for macro-micro radiation transport.
    %
    %   MacroMicroState represents the solution state for linear kinetic
    %   equations using orthogonal decomposition into macroscopic and
    %   microscopic components based on the paper "An Orthogonal Implicit
    %   Integrator for Multiscale Linear Kinetic Equations".
    %
    %   The macro-micro decomposition splits the distribution function:
    %
    %   \f[
    %       f(x,v,t) = u(x,v,t) + g(x,v,t)
    %   \f]
    %
    %   where \f$u\f$ is the macroscopic component and \f$g\f$ is the
    %   microscopic component.
    %
    % See also:
    %   physics.state.KineticState, approx.space.SumSpace

    properties (Constant)
        IsMacroComputed = 0b01 % Macro computation mask
        IsMicroComputed = 0b10 % Micro computation mask
        IsAllComputed = 0b11 % All computation mask
    end

    properties
        Levels % Basis levels
        Decomposition % Decomposition levels
        Cache % Precomputed values
        Status % Cache status
    end

    properties (Dependent)
        NXDims % Number of spatial dimensions
        NVDims % Number of velocity dimensions
        NXDofs % Number of spatial degrees of freedom
        NVDofs % Number of velocity degrees of freedom (total modes)
        NVMacroDofs % Number of velocity degrees of freedom (macroscopic modes)
        NVMicroDofs % Number of velocity degrees of freedom (microscopic modes)
        VLhsNodes % Velocity nodes in Cartesian coordinates
        VLhsWeights % Velocity weights in Cartesian coordinates
        VLhsBasisValues % Macroscopic basis function values at velocity nodes
        VLhsMass % Macroscopic basis function mass matrix
        VRhsNodes % Velocity nodes in Cartesian coordinates
        VRhsWeights % Velocity weights in Cartesian coordinates
        VRhsBasisValues % Microscopic basis function values at velocity nodes
        VNodes % All velocity nodes
        VInflowIndices % Velocity inflow indices for each boundary
        VOutflowIndices % Velocity outflow indices for each boundary
    end

    methods
        function obj = MacroMicroState(xDisc, vDisc)
            % MACROMICROSTATE Construct an instance of MacroMicroState.
            %
            %   obj = MacroMicroState(xDisc, vDisc) creates a macro-micro
            %   state with the specified spatial discretization @a xDisc
            %   and velocity SumSpace discretization @a vDisc.

            arguments
                xDisc approx.space.FiniteElementSpace
                vDisc approx.space.SumSpace
            end

            obj@physics.state.KineticState(xDisc, vDisc);

            %< Set degrees of freedom following equation: f = u + g
            obj.setDof('U', []); % Macroscopic component
            obj.setDof('G', []); % Microscopic component
            obj.setDof('F', []); % Kinetic component

            %< Set coefficients for radiation transport equation
            obj.setCoefficient('CS', []); % Scattering coefficient
            obj.setCoefficient('CA', []); % Absorption coefficient
            obj.setCoefficient('SU', []); % Macroscopic source
            obj.setCoefficient('SG', []); % Microscopic source

            %< Set history for free energy
            obj.setHistory('time', []);
            obj.setHistory('mass', []);
            obj.setHistory('freeEnergy', []); % Free energy

            %< Reset cache and status
            obj.setLevels();
            obj.reset();
            obj.Decomposition = [];
        end

        function newObj = refine(obj, nLevels)
            % REFINE Create refined macro-micro state.
            %
            %   newObj = refine(obj, nLevels) creates a refined macro-micro
            %   state by refining only the spatial discretization while
            %   keeping the velocity discretization unchanged.

            arguments
                obj physics.radiation.MacroMicroState
                nLevels{mustBeNonnegative, mustBeInteger}
            end

            newXDisc = obj.XDisc.refine(nLevels);
            newObj = physics.radiation.MacroMicroState(newXDisc, obj.VDisc);
        end
        
        function newObj = copy(obj)
            % COPY Create a shallow copy of the state.
            %
            %   newObj = copy(obj) creates a copy of the MacroMicroState
            %   with the same discretizations and copies of DOFs and
            %   coefficients.

            arguments
                obj physics.vlasov.MacroMicroState
            end

            newObj = physics.radiation.MacroMicroState(obj.XDisc, obj.VDisc);

            % Copy DOFs
            dofNames = fieldnames(obj.Dofs);
            for i = 1:length(dofNames)
                name = dofNames{i};
                newObj.setDof(name, obj.Dofs.(name));
            end

            % Copy coefficients
            coefNames = fieldnames(obj.Coefs);
            for i = 1:length(coefNames)
                name = coefNames{i};
                newObj.setCoefficient(name, obj.Coefs.(name));
            end

            % Copy history
            histNames = fieldnames(obj.History);
            for i = 1:length(histNames)
                name = histNames{i};
                newObj.setHistory(name, obj.History.(name));
            end
        end

        function obj = updateHistory(obj, t)
            % UPDATEHISTORY Update history data.
            %
            %   obj = updateHistory(obj, t) computes history data at
            %   time @a t.

            arguments
                obj physics.radiation.MacroMicroState
                t{mustBeNumeric}
            end

            obj.History.time(end+1) = t;

            m = obj.NVMacroDofs;
            n = obj.NVMicroDofs;
            ne = obj.XDisc.NMeshElements;
            nxq = obj.XDisc.Element.Volume.NPoints;
            wxq = obj.XDisc.Element.Volume.Weights;
            xRef = obj.XDisc.Element.Volume.Nodes;
            detJac = obj.XDisc.Mesh.computeElementJacobianDeterminants();

            %< Mass: rho
            if m > 0
                CU = obj.XDisc.eval(xRef, obj.Dofs.U(:, 1));
                CU = wxq * reshape(CU, nxq, []);
                CU = detJac(:).' * reshape(CU, ne, []);
                mass = sum(CU);
                obj.History.mass(end+1) = mass;
            end

            %< Free energy: (1/2) * int |f|^2 dx dv
            %  Assuming orthonormal velocity basis, this simplifies to
            %  (1/2) * sum_i ∫ |c_i(x)|^2 dx where c_i are modal coefficients
            freeEnergy = 0;

            if m > 0
                CU = obj.XDisc.eval(xRef, obj.Dofs.U);
                CU = wxq * reshape(CU.^2, nxq, []);
                CU = detJac(:).' * reshape(CU, ne, []);
                freeEnergy = freeEnergy + sum(CU) / 2;
            end

            if n > 0
                CG = obj.XDisc.eval(xRef, obj.Dofs.G);
                CG = wxq * reshape(CG.^2, nxq, []);
                CG = detJac(:).' * reshape(CG, ne, []);
                freeEnergy = freeEnergy + sum(CG) / 2;
            end

            obj.History.freeEnergy(end+1) = freeEnergy;
        end

        function obj = setScaling(obj, epsilon)
            % SETSCALING Set the scaling weights based on the scaling parameter.
            %
            %   obj = setScaling(obj, epsilon) sets the scaling weights
            %   based on the scaling parameter @a epsilon.

            arguments
                obj physics.radiation.MacroMicroState
                epsilon{mustBePositive}
            end

            m = obj.NVMacroDofs;
            n = obj.NVMicroDofs;

            w = epsilon .^ obj.Decomposition(1:m);
            obj.VDisc.setLhsWeight(w);
            w = epsilon .^ obj.Decomposition(m+(1:n));
            obj.VDisc.setRhsWeight(w);
        end

        function obj = setLevels(obj)
            % SETLEVELS Set the basis levels based on the velocity space.

            arguments
                obj physics.radiation.MacroMicroState
            end

            m = obj.NVMacroDofs;
            n = obj.NVMicroDofs;

            if m == 0
                obj.Levels = ones(1, n);
                return;
            end

            vBasis = obj.VDisc.Lhs.Element.Approximator.Basis.Lhs;

            switch class(vBasis)
                case 'core.function.LegendreBasis'
                    obj.Levels = [0:m-1, repmat(m, 1, n)];
                case 'core.function.FourierBasis'
                    k = ceil((m - 1) / 2);
                    l = [0, repelem(1:k, 2)];
                    obj.Levels(1:m) = l(1:m);
                    if n > 0
                        obj.Levels(m+(1:n)) = k + 1;
                    end
                case 'core.function.SphericalHarmonicBasis'
                    k = 0;
                    j = 0;
                    while j < m
                        s = repmat(k, 1, 2*k + 1);
                        i = j + 1:min(j + 2*k + 1, m);
                        obj.Levels(i) = s(1:length(i));
                        j = j + 2*k + 1;
                        k = k + 1;
                    end
                    if n > 0
                        obj.Levels(m+(1:n)) = k;
                    end
                otherwise
                    obj.Levels = [0:m-1, repmat(m, 1, n)];
            end
        end

        function obj = setDecomposition(obj, mode)
            % SETDECOMPOSITION Set the decomposition levels.
            %
            %   obj = setDecomposition(obj, mode) sets the decomposition
            %   levels based on the decomposition @a mode. Valid modes are
            %   'none', 'normal', 'diffusive', and 'hilbert'.

            arguments
                obj physics.radiation.MacroMicroState
                mode char{mustBeMember(mode, {'', 'partial', 'full'})}
            end

            m = obj.NVMacroDofs;
            n = obj.NVMicroDofs;

            switch mode
                case 'full'
                    obj.Decomposition = obj.Levels;
                case 'partial'
                    obj.Decomposition = [0, ones(1, m+n-1)];
                otherwise
                    obj.Decomposition = zeros(1, m + n);
            end
        end

        function F = distribution(obj, xRef, v)
            % DISTRIBUTION Compute distribution function.
            %
            %   F = distribution(obj, xRef, v) computes the distribution
            %   function at reference points @a xRef with velocity @a v.

            arguments
                obj physics.radiation.MacroMicroState
                xRef{mustBeNumeric}
                v{mustBeNumeric}
            end

            m = obj.NVMacroDofs;
            n = obj.NVMicroDofs;

            if m > 0 && n == 0
                U = obj.Dofs.U;
                U = obj.XDisc.eval(xRef, U); % (nx*nq, m)
                F = obj.VDisc.Lhs.eval(v, U.');
            elseif m == 0 && n > 0
                G = obj.Dofs.G;
                G = obj.XDisc.eval(xRef, G); % (nx*nq, n)
                F = obj.VDisc.Rhs.eval(v, G.');
            else
                U = obj.Dofs.U;
                U = obj.XDisc.eval(xRef, U); % (nx*nq, m)
                G = obj.Dofs.G;
                G = obj.XDisc.eval(xRef, G); % (nx*nq, n)
                F = obj.VDisc.eval(v, U.', G.'); % (nx*nq, nv)
            end
        end

        function obj = kineticReconstruct(obj, options)
            arguments
                obj physics.radiation.MacroMicroState
                options.v {mustBeNumeric} = []
            end

            if isempty(options.v)
                if ~isempty(obj.VDisc.Rhs)
                    v = obj.VDisc.Rhs.Element.Volume.Nodes;
                else
                    v = obj.VDisc.Lhs.Element.Volume.Nodes;
                end
            else
                v = options.v;
            end

            m = obj.NVMacroDofs;
            n = obj.NVMicroDofs;
            nl = obj.XDisc.NLocalDofs;
            ne = obj.XDisc.NMeshElements;
            xRef = obj.XDisc.Element.Volume.Nodes;
            
            if m > 0 && n == 0
                U = obj.Dofs.U;
                U = obj.XDisc.eval(xRef, U); % (nx*nq, m)
                F = obj.VDisc.Lhs.eval(v, U.');
            elseif m == 0 && n > 0
                G = obj.Dofs.G;
                G = obj.XDisc.eval(xRef, G); % (nx*nq, n)
                F = obj.VDisc.Rhs.eval(v, G.');
            else
                U = obj.Dofs.U;
                U = obj.XDisc.eval(xRef, U); % (nx*nq, m)
                G = obj.Dofs.G;
                G = obj.XDisc.eval(xRef, G); % (nx*nq, n)
                F = obj.VDisc.eval(v, U.', G.'); % (nx*nq, nv)
            end
            
            A = obj.XDisc.Element.Approximator;
            switch A.Type
                case 'modal'
                    D = obj.XDisc.Element.Volume;
                    F = A.embed(F, D.Values, D.Weights);
                case 'nodal'
                    F = A.embed(F);
                otherwise
                    core.except.error('UnsupportedProjector', ...
                        'Projector type %s not supported.', A.Type);
            end
            
            F = A.fit(F);
            F = reshape(F, nl*ne, []);
            obj.Dofs.F = F;
        end

        function obj = reset(obj)
            % RESET Resets the cache and status.

            arguments
                obj physics.radiation.MacroMicroState
            end

            obj.Cache = [];
            obj.Status = 0b00;
        end

        function obj = lazyFit(obj, f, options)
            % LAZYFIT Lazily fits the function handle to the state.
            %
            %   When both macro (m > 0) and micro (n > 0) modes are present,
            %   this function evaluates f at both macro and micro velocity
            %   nodes. Macro coefficients are computed using macro nodes,
            %   while micro coefficients need f evaluated at micro nodes.

            arguments
                obj physics.radiation.MacroMicroState
                f function_handle
                options.x double = []
                options.v double = []
                options.t double = 0
            end

            ng = obj.XDisc.NGlobalDofs;
            m = obj.NVMacroDofs;
            n = obj.NVMicroDofs;

            if nargin(f) == 1
                obj.Cache.C = obj.XDisc.fit(f);
                obj.Cache.C = reshape(obj.Cache.C, ng, []);
                return;
            end

            t = options.t;
            nArgs = nargin(f);

            %< For macro coefficients, use macro velocity nodes
            if m > 0
                VU = obj.VLhsNodes;
                BU = obj.VLhsBasisValues;
                wU = obj.VLhsWeights;
                nv = size(VU, 2);

                obj.Cache.F = zeros(ng, nv);
                for k = 1:nv
                    vk = VU(:, k);
                    if nArgs == 2
                        fk = @(x) f(x, repmat(vk, 1, size(x, 2)));
                    else
                        fk = @(x) f(x, repmat(vk, 1, size(x, 2)), t);
                    end
                    obj.Cache.F(:, k) = obj.XDisc.fit(fk);
                end

                S = obj.VDisc.Lhs.Element.Approximator.embed(obj.Cache.F.', BU, wU);
                obj.Cache.C = obj.VDisc.Lhs.Element.Approximator.project(S);
            end

            %< For micro coefficients, evaluate f at micro velocity nodes
            if n > 0
                VG = obj.VRhsNodes;
                nv = size(VG, 2);

                obj.Cache.F = zeros(ng, nv);
                for k = 1:nv
                    vk = VG(:, k);
                    if nArgs == 2
                        fk = @(x) f(x, repmat(vk, 1, size(x, 2)));
                    else
                        fk = @(x) f(x, repmat(vk, 1, size(x, 2)), t);
                    end
                    obj.Cache.F(:, k) = obj.XDisc.fit(fk);
                end
            end
        end

        function obj = lazyFEval(obj, f, options)
            % LAZYFEVAL Lazily evaluates the function handle to the state.
            %
            %   When both macro (m > 0) and micro (n > 0) modes are present,
            %   this function evaluates f at both macro and micro velocity
            %   nodes. Macro coefficients are computed using macro nodes,
            %   while micro coefficients need f evaluated at micro nodes.

            arguments
                obj physics.radiation.MacroMicroState
                f function_handle
                options.x double = []
                options.v double = []
                options.t double = 0
            end

            m = obj.NVMacroDofs;
            n = obj.NVMicroDofs;

            if isempty(options.x)
                I = obj.XDisc.Mesh.AllElementMultiIndices;
                xRef = obj.XDisc.Element.Volume.Nodes;
                XX = obj.XDisc.Mesh.collocate(xRef, I);
            else
                XX = options.x;
            end

            if nargin(f) == 1
                obj.Cache.C = f(XX);
                obj.Cache.C = reshape(obj.Cache.C, size(XX, 2), []);
                return;
            end

            t = options.t;

            %< When specific velocity is provided, just evaluate directly
            if ~isempty(options.v)
                VV = options.v;
                X = kron(ones(1, size(VV, 2)), XX);
                V = kron(VV, ones(1, size(XX, 2)));

                if nargin(f) == 2
                    F = f(X, V);
                else
                    F = f(X, V, t);
                end
                F = reshape(F, size(XX, 2), []);
                obj.Cache.F = F;
                obj.Cache.C = [];
                return;
            end

            %< For macro coefficients, use macro velocity nodes
            if m > 0
                VU = obj.VLhsNodes;
                BU = obj.VLhsBasisValues;
                wU = obj.VLhsWeights;

                X = kron(ones(1, size(VU, 2)), XX);
                V = kron(VU, ones(1, size(XX, 2)));

                if nargin(f) == 2
                    FU = f(X, V);
                else
                    FU = f(X, V, t);
                end
                FU = reshape(FU, size(XX, 2), []);
                obj.Cache.F = FU;

                S = obj.VDisc.Lhs.Element.Approximator.embed(FU.', BU, wU);
                obj.Cache.C = obj.VDisc.Lhs.Element.Approximator.project(S);
            end

            %< For micro coefficients, evaluate f at micro velocity nodes
            if n > 0
                VG = obj.VRhsNodes;

                X = kron(ones(1, size(VG, 2)), XX);
                V = kron(VG, ones(1, size(XX, 2)));

                if nargin(f) == 2
                    FG = f(X, V);
                else
                    FG = f(X, V, t);
                end
                FG = reshape(FG, size(XX, 2), []);
                obj.Cache.F = FG;
            end
        end

        function obj = lazyTraceEvaluate(obj, LFI, f, options)
            % LAZYTRACEEVALUATE Lazily evaluates the function handle on
            % boundary trace points.
            %
            %   obj = lazyTraceEvaluate(obj, LFI, f) evaluates function @a
            %   f on the boundary faces @LFI. For outflow velocities, uses
            %   the internal solution instead of the BC function.

            arguments
                obj physics.radiation.MacroMicroState
                LFI {mustBePositive, mustBeInteger}
                f function_handle
                options.x double = []
                options.v double = []
                options.t double = 0
            end

            m = obj.NVMacroDofs;
            n = obj.NVMicroDofs;
            nl = obj.XDisc.NLocalDofs;
            t = options.t;
            nArgs = nargin(f);

            %< Get boundary nodes
            EI = obj.XDisc.Mesh.getBoundaryElements(LFI);
            xRef = obj.XDisc.Element.Flux(LFI).Nodes;
            if ~isempty(options.x)
                XX = options.x;
            else
                XX = obj.XDisc.Mesh.collocate(xRef, EI);
            end
            nxx = size(XX, 2);

            core.except.assert(isa(obj.XDisc.Mesh, 'approx.mesh.Grid'), ...
                'InvalidMesh', 'Only support grid!');

            if m > 0 && n > 0 && isempty(options.v)
                VU = obj.VLhsNodes;
                VG = obj.VRhsNodes;

                %< Evaluate f at macro velocity nodes
                X = kron(ones(1, size(VU, 2)), XX);
                V = kron(VU, ones(1, nxx));
                if nArgs == 2
                    FU = f(X, V);
                else
                    FU = f(X, V, t);
                end
                FU = reshape(FU, nxx, []);

                %< Evaluate f at micro velocity nodes
                X = kron(ones(1, size(VG, 2)), XX);
                V = kron(VG, ones(1, nxx));
                if nArgs == 2
                    FG = f(X, V);
                else
                    FG = f(X, V, t);
                end
                FG = reshape(FG, nxx, []);

                %< Handle outflow
                normal = obj.XDisc.Mesh.computeOutwardNormals(EI, LFI);
                K = bsxfun(@plus, (EI(:).' - 1)*nl, (1:nl).');

                Co = obj.XDisc.eval(xRef, obj.Dofs.U(K(:), :), EI=EI);
                Uo = obj.VDisc.Lhs.Element.Approximator.project(Co.');

                JUo = find(sum(normal .* VU, 1) > 0);
                if ~isempty(JUo)
                    VUo = obj.VDisc.Lhs.Element.Volume.Nodes(:, JUo);
                    Go = obj.XDisc.eval(xRef, obj.Dofs.G(K(:), :), EI=EI);
                    Go = obj.VDisc.Rhs.eval(VUo, Go.');
                    FU(:, JUo) = obj.VDisc.Lhs.eval(VUo, Uo) + Go;
                end

                JGo = find(sum(normal .* VG, 1) > 0);
                if ~isempty(JGo)
                    VGo = obj.VDisc.Rhs.Element.Volume.Nodes(:, JGo);
                    Go = obj.XDisc.eval(xRef, obj.Dofs.G(K(:), JGo), EI=EI);
                    FG(:, JGo) = obj.VDisc.Lhs.eval(VGo, Uo) + Go;
                end

                %< Store micro values in Cache.F
                obj.Cache.F = FG;

                %< Project macro from FU
                BU = obj.VLhsBasisValues;
                wU = obj.VLhsWeights;
                S = obj.VDisc.Lhs.Element.Approximator.embed(FU.', BU, wU);
                obj.Cache.C = obj.VDisc.Lhs.Element.Approximator.project(S);
            else
                %< Get velocity nodes
                if ~isempty(options.v)
                    VV = options.v;
                elseif m > 0
                    VV = obj.VLhsNodes;
                else
                    VV = obj.VRhsNodes;
                end

                X = kron(ones(1, size(VV, 2)), XX);
                V = kron(VV, ones(1, nxx));

                if nArgs == 2
                    obj.Cache.F = f(X, V);
                else
                    obj.Cache.F = f(X, V, t);
                end
                obj.Cache.F = reshape(obj.Cache.F, nxx, []);

                %< Handle outflow
                normal = obj.XDisc.Mesh.computeOutwardNormals(EI, LFI);
                K = bsxfun(@plus, (EI(:).' - 1)*nl, (1:nl).');
                Jo = find(sum(normal .* VV, 1) > 0);
                if ~isempty(Jo) && m > 0 && n == 0 
                    Vo = obj.VDisc.Lhs.Element.Volume.Nodes(:, Jo);
                    Co = obj.XDisc.eval(xRef, obj.Dofs.U(K(:), :), EI=EI);
                    Uo = obj.VDisc.Lhs.Element.Approximator.project(Co.');
                    obj.Cache.F(:, Jo) = obj.VDisc.Lhs.eval(Vo, Uo);
                end

                %< Project to macro coefficients
                if m > 0
                    B = obj.VLhsBasisValues;
                    w = obj.VLhsWeights;
                    S = obj.VDisc.Lhs.Element.Approximator.embed(obj.Cache.F.', B, w);
                    obj.Cache.C = obj.VDisc.Lhs.Element.Approximator.project(S);
                end
            end
        end

        function U = macroFit(obj, f, options)
            % MACROFIT Fits the macroscopic component to the function
            % handle.
            %
            %   U = macroFit(obj, f) fits the macroscopic component to
            %   the function handle @a f.
            %
            %   U = macroFit(obj, f, x=X, v=V, t=t) uses specified
            %   evaluation points and time.

            arguments
                obj physics.radiation.MacroMicroState
                f function_handle
                options.x double = []
                options.v double = []
                options.t double = 0
            end

            m = obj.NVMacroDofs;

            if m == 0
                U = [];
                return;
            end

            if isempty(obj.Cache) || bitand(obj.Status, obj.IsMacroComputed)
                obj.lazyFit(f, x=options.x, v=options.v, t=options.t);
            end

            U = obj.macroImpl();

            obj.Status = bitor(obj.Status, obj.IsMacroComputed);
            if bitand(obj.Status, obj.IsAllComputed)
                obj.reset();
            end
        end

        function G = microFit(obj, f, options)
            % MICROFIT Fits the microscopic component to the function
            % handle.
            %
            %   G = microFit(obj, f) fits the microscopic component to
            %   the function handle @a f.
            %
            %   G = microFit(obj, f, x=X, v=V, t=t) uses specified
            %   evaluation points and time.

            arguments
                obj physics.radiation.MacroMicroState
                f function_handle
                options.x double = []
                options.v double = []
                options.t double = 0
            end

            n = obj.NVMicroDofs;

            if n == 0
                G = [];
                return;
            end

            if isempty(obj.Cache) || bitand(obj.Status, obj.IsMicroComputed)
                obj.lazyFit(f, x=options.x, v=options.v, t=options.t);
            end

            G = obj.microImpl();

            obj.Status = bitor(obj.Status, obj.IsMicroComputed);
            if bitand(obj.Status, obj.IsAllComputed)
                obj.reset();
            end
        end

        function U = macroFEval(obj, f, options)
            % MACROFEVAL Evaluates the macroscopic component from the
            % function handle.
            %
            %   U = macroFEval(obj, f) evaluates the macroscopic component
            %   from the function handle @a f.
            %
            %   U = macroFEval(obj, f, x=X, v=V, t=t) uses specified
            %   evaluation points and time.

            arguments
                obj physics.radiation.MacroMicroState
                f function_handle
                options.x double = []
                options.v double = []
                options.t double = 0
            end

            m = obj.NVMacroDofs;

            if m == 0
                U = [];
                return;
            end

            if isempty(obj.Cache) || bitand(obj.Status, obj.IsMacroComputed)
                obj.lazyFEval(f, x=options.x, v=options.v, t=options.t);
            end

            U = obj.macroImpl();

            obj.Status = bitor(obj.Status, obj.IsMacroComputed);
            if bitand(obj.Status, obj.IsAllComputed)
                obj.reset();
            end
        end

        function G = microFEval(obj, f, options)
            % MICROFEVAL Evaluates the microscopic component from the
            % function handle.
            %
            %   G = microFEval(obj, f) evaluates the microscopic component
            %   from the function handle @a f.
            %
            %   G = microFEval(obj, f, x=X, v=V, t=t) uses specified
            %   evaluation points and time.

            arguments
                obj physics.radiation.MacroMicroState
                f function_handle
                options.x double = []
                options.v double = []
                options.t double = 0
            end

            n = obj.NVMicroDofs;

            if n == 0 && isempty(options.v)
                G = [];
                return;
            end

            if isempty(obj.Cache) || bitand(obj.Status, obj.IsMicroComputed)
                obj.lazyFEval(f, x=options.x, v=options.v, t=options.t);
            end

            %< When specific velocity is provided, return raw function values
            if ~isempty(options.v)
                G = obj.Cache.F;
                obj.reset();
                return;
            end

            G = obj.microImpl();

            obj.Status = bitor(obj.Status, obj.IsMicroComputed);
            if bitand(obj.Status, obj.IsAllComputed)
                obj.reset();
            end
        end

        function F = kineticFEval(obj, f, options)
            % KINETICFEVAL Evaluate the distribution function from the
            % function handle.
            %
            %   F = kineticFEval(obj, f) evaluates the distribution function
            %   from the function handle @a f at the scheme's Rhs velocity
            %   nodes.
            %
            %   F = kineticFEval(obj, f, x=X, v=V, t=t) uses specified
            %   evaluation points, velocity nodes, and time.

            arguments
                obj physics.radiation.MacroMicroState
                f function_handle
                options.x double = []
                options.v double = []
                options.t double = 0
            end

            if isempty(options.v)
                options.v = obj.VRhsNodes;
            end

            obj.lazyFEval(f, x=options.x, v=options.v, t=options.t);
            F = obj.Cache.F;
            obj.reset();
        end

        function U = macroTraceEvaluate(obj, i, f, options)
            % MACROTRACEEVALUATE Evaluates the macroscopic component from
            % the function handle on boundary trace points.
            %
            %   U = macroTraceEvaluate(obj, i, f) evaluates the macroscopic
            %   component from function @a f on boundary face @a i.

            arguments
                obj physics.radiation.MacroMicroState
                i{mustBePositive, mustBeInteger}
                f function_handle
                options.x double = []
                options.v double = []
                options.t double = 0
                options.w double = []
            end

            m = obj.NVMacroDofs;

            if m == 0
                U = [];
                return;
            end

            if isempty(obj.Cache) || bitand(obj.Status, obj.IsMacroComputed)
                obj.lazyTraceEvaluate(i, f, x=options.x, v=options.v, t=options.t);
            end

            U = obj.macroImpl();

            obj.Status = bitor(obj.Status, obj.IsMacroComputed);
            if bitand(obj.Status, obj.IsAllComputed)
                obj.reset();
            end
        end

        function G = microTraceEvaluate(obj, i, f, options)
            % MICROTRACEEVALUATE Evaluates the microscopic component from
            % the function handle on boundary trace points.
            %
            %   G = microTraceEvaluate(obj, i, f) evaluates the microscopic
            %   component from function @a f on boundary face @a i.

            arguments
                obj physics.radiation.MacroMicroState
                i{mustBePositive, mustBeInteger}
                f function_handle
                options.x double = []
                options.v double = []
                options.t double = 0
                options.w double = []
            end

            n = obj.NVMicroDofs;

            if n == 0
                G = [];
                return;
            end

            if isempty(obj.Cache) || bitand(obj.Status, obj.IsMicroComputed)
                obj.lazyTraceEvaluate(i, f, x=options.x, v=options.v, t=options.t);
            end

            if ~isempty(options.v)
                G = obj.Cache.F;
                obj.reset();
                return;
            end

            G = obj.microImpl();
            
            obj.Status = bitor(obj.Status, obj.IsMicroComputed);
            if bitand(obj.Status, obj.IsAllComputed)
                obj.reset();
            end
        end
    
        function U = macroImpl(obj)
            % MACROIMPL Implements the macroscopic component from the
            % cached coefficients.

            arguments
                obj physics.radiation.MacroMicroState
            end

            M = obj.VDisc.Lhs.Element.Approximator.Mass;
            U = (M * obj.Cache.C).';
            if ~isempty(obj.VDisc.LhsWeight)
                U = U ./ obj.VDisc.LhsWeight(:).';
            end
        end

        function G = microImpl(obj)
            % MICROIMPL Implements the microscopic component from the
            % cached coefficients.
            %
            %   The micro component is computed as the residual after
            %   removing the macro contribution:
            %
            %   \f[
            %       g(x,v_j) = f(x,v_j) - \sum_k c_k(x) \phi_k(v_j)
            %   \f]
            %
            %   where \f$\phi_k(v)\f$ are macro basis functions evaluated at
            %   the micro velocity nodes \f$v_j\f$.

            arguments
                obj physics.radiation.MacroMicroState
            end

            m = obj.NVMacroDofs;
            n = obj.NVMicroDofs;

            if n == 0
                G = [];
                return;
            end

            V = obj.VDisc.Rhs.Element.Volume.Nodes;
            F = obj.Cache.F;
            if m > 0
                C = obj.Cache.C;
                G = obj.VDisc.invRhs(F, V, C);
            else
                G = F;
            end
        end

    end

    methods
        function n = get.NXDims(obj)
            % GET.NXDIMS Returns the number of spatial dimensions.

            n = obj.XDisc.NDims;
        end

        function n = get.NVDims(obj)
            % GET.NVDIMS Returns the number of velocity dimensions.

            n = obj.VDisc.NDims;
        end

        function n = get.NXDofs(obj)
            % GET.NXDOFS Returns the number of spatial degrees of freedom.

            n = obj.XDisc.NGlobalDofs;
        end

        function n = get.NVDofs(obj)
            % GET.NVDOFS Returns the total number of velocity degrees of
            % freedom.

            n = obj.VDisc.NDofs;
        end

        function n = get.NVMacroDofs(obj)
            % GET.NVMACRODOFS Returns the number of macroscopic velocity
            % modes.

            n = obj.VDisc.NLhsDofs;
        end

        function n = get.NVMicroDofs(obj)
            % GET.NVMICRODOFS Returns the number of microscopic velocity
            % modes.

            n = obj.VDisc.NRhsDofs;
        end

        function V = get.VLhsNodes(obj)
            % GET.VLHSNODES Returns the velocity nodes in Cartesian
            % coordinates.

            if isempty(obj.VDisc.Lhs)
                V = [];
                return;
            end

            nxd = obj.NXDims;
            V = obj.VDisc.Lhs.Element.Volume.Nodes;
            G = obj.VDisc.Lhs.Element.Geometry;
            reduction = obj.VDisc.Lhs.Element.DimReduction;
            if isa(G, 'core.geometry.Sphere')
                if strcmpi(reduction, 'topology')
                    V = G.sphericalToCartesian(V);
                else
                    V(1, :) = acos(-V(1, :));
                    V = G.sphericalToCartesian(V);
                end
            end
            V = V(1:nxd, :);
        end

        function w = get.VLhsWeights(obj)
            % GET.VLHSWEIGHTS Returns the velocity weights.

            if isempty(obj.VDisc.Lhs)
                w = [];
                return;
            end

            w = obj.VDisc.Lhs.Element.Volume.Weights;
        end

        function B = get.VLhsBasisValues(obj)
            % GET.VLHSBASISVALUES Returns the macroscopic basis function
            % values.

            if isempty(obj.VDisc.Lhs)
                B = [];
                return;
            end

            B = obj.VDisc.Lhs.Element.Volume.Values;
        end

        function M = get.VLhsMass(obj)
            % GET.VLHSMASS Returns the macroscopic basis function mass
            % matrix.

            if isempty(obj.VDisc.Lhs)
                M = [];
                return;
            end

            M = obj.VDisc.Lhs.Element.Approximator.Mass;
        end

        function V = get.VRhsNodes(obj)
            % GET.VRHSNODES Returns the velocity nodes in Cartesian
            % coordinates.

            if isempty(obj.VDisc.Rhs)
                V = [];
                return;
            end

            nxd = obj.NXDims;
            V = obj.VDisc.Rhs.Element.Volume.Nodes;
            G = obj.VDisc.Rhs.Element.Geometry;
            reduction = obj.VDisc.Rhs.Element.DimReduction;
            if isa(G, 'core.geometry.Sphere')
                if strcmpi(reduction, 'topology')
                    V = G.sphericalToCartesian(V);
                else
                    V(1, :) = acos(-V(1, :));
                    V = G.sphericalToCartesian(V);
                end
            end
            V = V(1:nxd, :);
        end

        function w = get.VRhsWeights(obj)
            % GET.VRHSWEIGHTS Returns the velocity weights.

            if isempty(obj.VDisc.Rhs)
                w = [];
                return;
            end

            w = obj.VDisc.Rhs.Element.Volume.Weights;
        end

        function B = get.VRhsBasisValues(obj)
            % GET.VRHSBASISVALUES Returns the microscopic basis function
            % values.

            if isempty(obj.VDisc.Rhs)
                B = [];
                return;
            end

            B = obj.VDisc.Rhs.Element.Volume.Values;
        end

        function V = get.VNodes(obj)
            % GET.VNODES Returns all velocity nodes (from Lhs if
            % available).

            if ~isempty(obj.VDisc.Lhs)
                V = obj.VLhsNodes;
            else
                V = obj.VRhsNodes;
            end
        end

        function I = get.VInflowIndices(obj)
            % GET.VINFLOWINDICES Returns the velocity inflow indices for
            % each boundary face.

            nxd = obj.NXDims;
            V = obj.VNodes;
            I = cell(1, 2*nxd);
            for d = 1:nxd
                I{2*d-1} = find(V(d, :) > 0);
                I{2*d} = find(V(d, :) < 0);
            end
        end

        function I = get.VOutflowIndices(obj)
            % GET.VOUTFLOWINDICES Returns the velocity outflow indices for
            % each boundary face.

            nxd = obj.NXDims;
            V = obj.VNodes;
            I = cell(1, 2*nxd);
            for d = 1:nxd
                I{2*d-1} = find(V(d, :) < 0);
                I{2*d} = find(V(d, :) > 0);
            end
        end
    end

    methods (Static)
        function obj = load(filename)
            % LOAD Load kinetic state from file.
            %
            %   obj = MacroMicroState.load(filename) loads a kinetic state
            %   object from the specified MAT file.

            arguments
                filename{mustBeTextScalar}
            end

            core.except.assert(exist(filename, 'file') == 2, ...
                'InvalidInput', 'File %s does not exist', filename);

            state = load(filename).state;

            requiredFields = {'xDisc', 'vDisc'};
            [tf, mf] = core.except.hasFields(state, requiredFields);
            core.except.assert(tf, 'InvalidInput', ...
                'Loaded file missing required field: %s', mf);

            obj = physics.radiation.MacroMicroState(state.xDisc, state.vDisc);

            obj.Dofs = state.dofs;
            obj.Coefs = state.coefs;
            obj.History = state.history;

            fprintf('[M] MacroMicroState loaded from: %s\n', filename);
        end
    end
end
