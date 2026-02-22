classdef HermiteState < physics.state.KineticState
    % HERMITESTATE

    properties
        VBBox % Velocity bounding box
    end

    properties (Dependent)
        NXDofs % Number of spatial degrees of freedom
        NVDofs % Number of velocity degrees of freedom
    end

    methods
        function obj = HermiteState(xDisc, vDisc)
            % HERMITESTATE Construct an instance of HermiteState.
            %
            %   obj = HermiteState(xDisc, vDisc) creates a VpHermite state
            %   with the specified spatial discretization @a xDisc and
            %   velocity discretization @a vDisc.

            arguments
                xDisc approx.space.FiniteElementSpace
                vDisc approx.space.SpectralSpace
            end

            obj@physics.state.KineticState(xDisc, vDisc);
            obj.setDof('D', []); % Hermite coefficients
            obj.setDof('P', []); % Potential
            obj.setDof('E', []); % Electric field
            obj.setCoefficient('T0', []); % Temperature parameter
            obj.setCoefficient('rhoi', []); % Ion density coefficient
            obj.setHistory('time', []); % Time history
            obj.setHistory('mass', []); % Mass = D_0
            obj.setHistory('momentum', []); % Momentum = D_1
            obj.setHistory('kineticEnergy', []); % Kinetic Energy = (\sqrt(2)D_2 + D0)/2
            obj.setHistory('potentialEnergy', []); % Potential Energy = \norm{E}_2^2 / 2
            obj.setHistory('totalEnergy', []); % Total Energy = Kinetic Energy + Potential Energy
            obj.setHistory('L2Entropy', []); % \norm{f-\rho M}_{L^2(f_\infty^{-1})}
            obj.setHistory('L2FDistance', []); % \norm{f-f_\infty}_{L^2(f_\infty^{-1})}
            obj.setHistory('L2RhoDistance', []); % \norm{\rho-\rho_\infty}_{L^2(\rho_\infty^{-1})}
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

        function newObj = refine(obj, nLevels)
            % REFINE Create refined Hermite state.
            %
            %   newObj = refine(obj, nLevels) creates a refined Hermite
            %   state by refining only the spatial discretization while
            %   keeping the velocity discretization unchanged.

            arguments
                obj physics.vlasov.HermiteState
                nLevels{mustBeNonnegative, mustBeInteger}
            end

            newXDisc = obj.XDisc.refine(nLevels);
            oldVElement = obj.VDisc.Element;
            newNVDofs = oldVElement.NDofs * 2;
            newT0 = oldVElement.Approximator.Basis.Rhs.Temperature;
            newVElement = approx.element.L2OrthotopeElement.hermite(1, newNVDofs, T = newT0);
            newVDisc = approx.space.SpectralSpace(newVElement);
            newObj = physics.vlasov.HermiteState(newXDisc, newVDisc);
        end

        function F = distribution(obj, xRef, v)
            % DISTRIBUTION Compute distribution function.
            %
            %   F = distribution(obj, xRef, v) computes the distribution
            %   function at reference points @a xRef with velocity @a v.

            arguments
                obj physics.vlasov.HermiteState
                xRef{mustBeNumeric}
                v{mustBeNumeric}
            end

            D = obj.Dofs.D;
            D = obj.XDisc.eval(xRef, D); % (nx*nq, nh)
            F = obj.VDisc.eval(v, D.'); % (nx*nq, nv)
        end

        function obj = updateHistory(obj, t)
            % UPDATEHISTORY Update history data.
            %
            %   obj = updateHistory(obj, t) computes history data at
            %   time @a t.

            arguments
                obj physics.vlasov.HermiteState
                t{mustBeNumeric}
            end

            obj.History.time(end+1) = t;

            T0 = obj.Coefs.T0;
            D0 = obj.Dofs.D(:, 1);
            D1 = obj.Dofs.D(:, 2);
            D2 = obj.Dofs.D(:, 3);
            E = obj.Dofs.E;

            n = obj.XDisc.Element.Volume.NPoints;
            w = obj.XDisc.Element.Volume.Weights;
            xRef = obj.XDisc.Element.Volume.Nodes;
            detJac = obj.XDisc.Mesh.computeElementJacobianDeterminants();

            %< Mass: D_0
            mass = obj.XDisc.eval(xRef, D0);
            mass = w * reshape(mass, n, []);
            mass = sum(mass.*detJac(:).');
            obj.History.mass(end+1) = mass;

            %< Momentum: D_1
            momentum = obj.XDisc.eval(xRef, D1);
            momentum = w * reshape(momentum, n, []);
            momentum = sum(momentum.*detJac(:).');
            obj.History.momentum(end+1) = sqrt(T0) * momentum;

            %< Kinetic energy: (\sqrt{2}D_2 + D0)/2
            kineticEnergy = obj.XDisc.eval(xRef, (sqrt(2) * D2 + D0));
            kineticEnergy = w * reshape(kineticEnergy, n, []);
            kineticEnergy = sum(kineticEnergy.*detJac(:).');
            obj.History.kineticEnergy(end+1) = T0 * kineticEnergy / 2;

            %< Potential energy: \norm{E}_2^2 / 2
            Ep = obj.XDisc.eval(xRef, E);
            potentialEnergy = w * reshape(Ep.^2, n, []);
            potentialEnergy = sum(potentialEnergy.*detJac(:).');
            obj.History.potentialEnergy(end+1) = potentialEnergy / 2;

            %< Total energy: Kinetic energy + Potential energy
            obj.History.totalEnergy(end+1) = kineticEnergy + potentialEnergy;

            %< L^2 entropy: \norm{f-\rho M}_{L^2(f_\infty^{-1})}^2 = \sum_{k>0} \norm{D_k}^2
            L2Entropy = 0;
            for k = 2:obj.NVDofs
                DkFn = @(x) obj.XDisc.eval(xRef, obj.Dofs.D(:, k));
                L2Dk = obj.XDisc.feval(@(x) (DkFn(x)), xRef);
                L2Dk = w * reshape(L2Dk.^2, n, []);
                L2Dk = sum(L2Dk.*detJac(:).');
                L2Entropy = L2Entropy + L2Dk;
            end
            obj.History.L2Entropy(end+1) = sqrt(L2Entropy);

            %< \norm{f-f_\infty}_{L^2(f_\infty^{-1})}^2 = \norm{D_0-\sqrt{\rho_\infty}}^2 + \sum_{k>0} \norm{D_k}^2
            D0Fn = @(x) obj.XDisc.eval(xRef, D0);
            L2D0 = obj.XDisc.feval(@(x) (D0Fn(x) - 1), xRef);
            L2D0 = w * reshape(L2D0.^2, n, []);
            L2D0 = sum(L2D0.*detJac(:).');
            obj.History.L2RhoDistance(end+1) = sqrt(L2D0);
            obj.History.L2FDistance(end+1) = sqrt(L2D0 + L2Entropy);
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

            obj = physics.vlasov.HermiteState(state.xDisc, state.vDisc);

            obj.Dofs = state.dofs;
            obj.Coefs = state.coefs;
            obj.History = state.history;

            fprintf('[M] HermiteState loaded from: %s\n', filename);
        end
    end
end
