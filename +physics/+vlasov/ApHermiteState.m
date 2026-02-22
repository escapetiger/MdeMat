classdef ApHermiteState < physics.state.KineticState
    % APHERMITESTATE

    properties
        VBBox % Velocity bounding box
    end

    properties (Dependent)
        NXDofs % Number of spatial degrees of freedom
        NVDofs % Number of velocity degrees of freedom
    end

    methods
        function obj = ApHermiteState(xDisc, vDisc)
            % APHERMITESTATE Construct an instance of ApHermiteState.
            %
            %   obj = ApHermiteState(xDisc, vDisc) creates an AP Hermite
            %   state with the specified spatial discretization @a xDisc
            %   and velocity discretization @a vDisc.

            arguments
                xDisc approx.space.FiniteElementSpace
                vDisc approx.space.SpectralSpace
            end

            obj@physics.state.KineticState(xDisc, vDisc);
            obj.setDof('D', []); % Hermite coefficients
            obj.setDof('Psi', []); % Potential
            obj.setDof('Q', []); % Electric field
            obj.setCoefficient('T0', []); % Temperature parameter
            obj.setCoefficient('rhoi', []); % Ion density coefficient
            obj.setCoefficient('rhoInf', []); % Density equilibrium
            obj.setCoefficient('phiInf', []); % Potential equilibrium
            obj.setCoefficient('sqrtRhoInf', []); % Square root of density equilibrium
            obj.setCoefficient('EInf', []); % eletric equilibrium
            obj.setHistory('time', []); % Time history
            obj.setHistory('mass', []); % Mass
            obj.setHistory('momentum', []); % Momentum
            obj.setHistory('kineticEnergy', []); % Kinetic Energy
            obj.setHistory('potentialEnergy', []); % Potential Energy
            obj.setHistory('totalEnergy', []); % Total Energy = Kinetic Energy + Potential Energy
            obj.setHistory('L2Entropy', []); % |f-\rho M|
            obj.setHistory('L2FDistance', []); % |f-f_\infty|
            obj.setHistory('L2RhoDistance', []); % |\rho-\rho_\infty|
        end

        function bbox = get.VBBox(obj)
            % GET.VBBox Returns the velocity bounding box.

            bbox = obj.VDisc.Element.BBox;
        end

        function n = get.NXDofs(obj)
            % GET.NXDofs Returns the number of spatial degrees of freedom.

            n = obj.XDisc.NGlobalDofs;
        end

        function n = get.NVDofs(obj)
            % GET.NVDofs Returns the number of velocity degrees of freedom.

            n = obj.VDisc.Element.NDofs;
        end

        function newObj = copy(obj)
            % COPY Create a shallow copy of the state.
            %
            %   newObj = copy(obj) creates a copy of the ApHermiteState
            %   with the same discretizations and copies of DOFs and
            %   coefficients.

            arguments
                obj physics.vlasov.ApHermiteState
            end

            newObj = physics.vlasov.ApHermiteState(obj.XDisc, obj.VDisc);

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

        function newObj = refine(obj, nLevels)
            % REFINE Create refined Hermite state.
            %
            %   newObj = refine(obj, nLevels) creates a refined Hermite
            %   state by refining only the spatial discretization while
            %   keeping the velocity discretization unchanged.

            arguments
                obj physics.vlasov.ApHermiteState
                nLevels{mustBeNonnegative, mustBeInteger}
            end

            newXDisc = obj.XDisc.refine(nLevels);
            oldVElement = obj.VDisc.Element;
            newNVDofs = oldVElement.NDofs * 2;
            newT0 = oldVElement.Approximator.Basis.Rhs.Temperature;
            newVElement = approx.element.L2OrthotopeElement.hermite(1, newNVDofs, T = newT0);
            newVDisc = approx.space.SpectralSpace(newVElement);
            newObj = physics.vlasov.ApHermiteState(newXDisc, newVDisc);
        end

        function F = distribution(obj, xRef, v)
            % DISTRIBUTION Compute distribution function.
            %
            %   F = distribution(obj, xRef, v) computes the distribution
            %   function at reference points @a xRef with velocity @a v.

            arguments
                obj physics.vlasov.ApHermiteState
                xRef{mustBeNumeric}
                v{mustBeNumeric}
            end

            sqrtRhoInf = obj.Coefs.sqrtRhoInf;
            D = obj.Dofs.D;
            sqrtRhoInf = obj.XDisc.eval(xRef, sqrtRhoInf); % (nx*nq, 1)
            D = obj.XDisc.eval(xRef, D); % (nx*nq, nh)
            C = sqrtRhoInf .* D;
            F = obj.VDisc.eval(v, C.'); % (nx*nq, nv)
        end

        function F = equilibrium(obj, xRef, v)
            % EQUILIBRIUM Compute equilibrium.
            %
            %   F = equilibrium(obj, xRef, v) computes the equilibrium at
            %   reference points @a xRef with velocity @a v.

            arguments
                obj physics.vlasov.ApHermiteState
                xRef {mustBeNumeric}
                v {mustBeNumeric}    
            end

            % DInf = zeros(size(obj.Dofs.D));
            % DInf(:, 1) = obj.Coefs.sqrtRhoInf;
            % DInf = obj.XDisc.eval(xRef, DInf);
            % F = obj.VDisc.eval(v, DInf.');
            rhoInf = obj.XDisc.eval(xRef, obj.Coefs.rhoInf);
            M = obj.VDisc.Element.Approximator.Basis.Rhs.eval(v); % (1,nv)
            F = rhoInf .* M;
        end

        function obj = updateHistory(obj, t)
            % UPDATEHISTORY Update history data.
            %
            %   obj = updateHistory(obj, t) computes history data at
            %   time @a t.

            arguments
                obj physics.vlasov.ApHermiteState
                t{mustBeNumeric}
            end

            obj.History.time(end+1) = t;

            T0 = obj.Coefs.T0;
            sqrtRhoInf = obj.Coefs.sqrtRhoInf;
            D0 = obj.Dofs.D(:, 1);
            D1 = obj.Dofs.D(:, 2);
            D2 = obj.Dofs.D(:, 3);
            C0 = sqrtRhoInf .* D0;
            C1 = sqrtRhoInf .* D1;
            C2 = sqrtRhoInf .* D2;
            E = sqrt(T0) .* sqrtRhoInf .* obj.Dofs.Q;

            ne = obj.XDisc.NMeshElements;
            nq = obj.XDisc.Element.Volume.NPoints;
            wq = obj.XDisc.Element.Volume.Weights;
            xRef = obj.XDisc.Element.Volume.Nodes;
            detJac = obj.XDisc.Mesh.computeElementJacobianDeterminants();

            %< Mass: C_0
            mass = obj.XDisc.eval(xRef, C0);
            mass = wq * reshape(mass, nq, []);
            mass = detJac(:).' * reshape(mass, ne, []);
            obj.History.mass(end+1) = mass;

            %< Momentum: sqrt(T0) * C_1
            momentum = obj.XDisc.eval(xRef, C1);
            momentum = wq * reshape(momentum, nq, []);
            momentum = detJac(:).' * reshape(momentum, ne, []);
            obj.History.momentum(end+1) = sqrt(T0) * momentum;

            %< Kinetic energy: T0 * (\sqrt{2}C_2 + C0)/2
            kineticEnergy = obj.XDisc.eval(xRef, (sqrt(2) * C2 + C0));
            kineticEnergy = wq * reshape(kineticEnergy, nq, []);
            kineticEnergy = detJac(:).' * reshape(kineticEnergy, ne, []);
            obj.History.kineticEnergy(end+1) = T0 * kineticEnergy / 2;

            %< Potential energy: \norm{E}_2^2 / (2*T0)
            Ep = obj.XDisc.eval(xRef, E);
            potentialEnergy = wq * reshape(Ep.^2, nq, []);
            potentialEnergy = detJac(:).' * reshape(potentialEnergy, ne, []);
            obj.History.potentialEnergy(end+1) = potentialEnergy / (2 * T0);

            %< Total energy: Kinetic energy + Potential energy
            obj.History.totalEnergy(end+1) = kineticEnergy + potentialEnergy;

            %< Distance
            DInfDofs = zeros(size(obj.Dofs.D));
            DInfDofs(:, 1) = sqrtRhoInf;
            DFn = @(x) obj.XDisc.eval(xRef, obj.Dofs.D);
            DInfFn = @(x) obj.XDisc.eval(xRef, DInfDofs);
            L2D = obj.XDisc.feval(@(x) DFn(x) - DInfFn(x), xRef);
            L2D = wq * reshape(L2D.^2, nq, []);
            L2D = detJac(:).' * reshape(L2D, ne, []);
                
            %< L^2 entropy: \norm{f-\rho M}_{L^2(f_\infty^{-1})}^2 = \sum_{k>0} \norm{D_k}^2
            obj.History.L2Entropy(end+1) = sqrt(sum(L2D(2:end)));

            %< \norm{f-f_\infty}_{L^2(f_\infty^{-1})}^2 = \norm{D_0-\sqrt{\rho_\infty}}^2 + \sum_{k>0} \norm{D_k}^2
            obj.History.L2RhoDistance(end+1) = sqrt(L2D(1));
            obj.History.L2FDistance(end+1) = sqrt(sum(L2D));
        end
    end

    methods (Static)
        function obj = load(filename)
            % LOAD Load kinetic state from file.
            %
            %   obj = KineticState.load(filename) loads a kinetic state
            %   object from the specified MAT file and returns a new
            %   KineticState instance.

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

            obj = physics.vlasov.ApHermiteState(state.xDisc, state.vDisc);

            obj.Dofs = state.dofs;
            obj.Coefs = state.coefs;
            obj.History = state.history;

            fprintf('[M] ApHermiteState loaded from: %s\n', filename);
        end
    end
end
