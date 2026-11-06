classdef ApHermiteState < mdemat.common.fem.KineticState
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
                xDisc mdemat.common.fem.FiniteElementSpace
                vDisc mdemat.common.fem.SpectralSpace
            end

            obj@mdemat.common.fem.KineticState(xDisc, vDisc);
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

            bbox = obj.VDisc.Element.Geometry.BBox;
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
                obj mdemat.application.plasma.ApHermiteState
            end

            newObj = mdemat.application.plasma.ApHermiteState(obj.XDisc, obj.VDisc);

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
            %   state by refining the spatial mesh by 2^nLevels in each
            %   direction and doubling the Hermite mode count. Temperature
            %   and basis normalization are preserved; coefficients are not
            %   transferred, since a convergence run initializes each level.

            arguments
                obj mdemat.application.plasma.ApHermiteState
                nLevels{mustBeNonnegative, mustBeInteger}
            end

            newXDisc = obj.XDisc.refine(nLevels);
            oldVElement = obj.VDisc.Element;
            newNVDofs = oldVElement.NDofs * 2;
            newT0 = oldVElement.Approximator.Basis.Rhs.Temperature;
            newVElement = mdemat.common.fem.L2OrthotopeElement.hermite(1, newNVDofs, T = newT0);
            newVDisc = mdemat.common.fem.SpectralSpace(newVElement);
            newObj = mdemat.application.plasma.ApHermiteState(newXDisc, newVDisc);
        end

        function F = distribution(obj, xRef, v)
            % DISTRIBUTION Compute distribution function.
            %
            %   F = distribution(obj, xRef, v) computes the distribution
            %   function at reference points @a xRef with velocity @a v.

            arguments
                obj mdemat.application.plasma.ApHermiteState
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
                obj mdemat.application.plasma.ApHermiteState
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
            % UPDATEHISTORY Integrate physical samples, never products of DOFs.
            % D are coefficients of f=sqrt(rhoInf)*M_T*sum D_k He_k/sqrt(k!).
            % The field diagnostic retains the scaled convention int E^2/(2*T).
            arguments
                obj mdemat.application.plasma.ApHermiteState
                t (1, 1) {mustBeNumeric, mustBeReal, mustBeFinite}
            end
            space = obj.XDisc;
            reference = space.Element.Volume.Nodes;
            determinants = abs(space.Mesh.computeElementJacobianDeterminants());
            if isscalar(determinants)
                determinants = repmat(determinants, space.NMeshElements, 1);
            end
            weights = kron(determinants(:), space.Element.Volume.Weights(:));
            values = space.eval(reference, obj.Dofs.D);
            equilibriumScale = space.eval(reference, obj.Coefs.sqrtRhoInf);
            coefficients = values.*equilibriumScale;
            density = coefficients(:, 1);
            firstMoment = zeros(size(density));
            secondMode = zeros(size(density));
            if size(values, 2) >= 2
                firstMoment = coefficients(:, 2);
            end
            if size(values, 2) >= 3
                secondMode = coefficients(:, 3);
            end
            temperature = obj.Coefs.T0;
            field = sqrt(temperature)*equilibriumScale.*space.eval(reference, obj.Dofs.Q);
            kinetic = temperature/2*(weights.'*(density+sqrt(2)*secondMode));
            potential = weights.'*abs(field).^2/(2*temperature);
            distance = values;
            distance(:, 1) = distance(:, 1)-equilibriumScale;
            squaredDistance = weights.'*abs(distance).^2;

            obj.History.time(end+1) = t;
            obj.History.mass(end+1) = weights.'*density;
            obj.History.momentum(end+1) = sqrt(temperature)*(weights.'*firstMoment);
            obj.History.kineticEnergy(end+1) = kinetic;
            obj.History.potentialEnergy(end+1) = potential;
            obj.History.totalEnergy(end+1) = kinetic+potential;
            obj.History.L2Entropy(end+1) = sqrt(sum(squaredDistance(2:end)));
            obj.History.L2RhoDistance(end+1) = sqrt(squaredDistance(1));
            obj.History.L2FDistance(end+1) = sqrt(sum(squaredDistance));
        end
    end

    methods
        function record = snapshot(obj, options)
            %SNAPSHOT Preserve physical state parameters; never persist derived caches.
            arguments
                obj mdemat.application.plasma.ApHermiteState
                options.omitCallbacks (1, 1) logical = false
            end
            record = snapshot@mdemat.common.fem.KineticState(obj, ...
                omitCallbacks=options.omitCallbacks);
            record.kind = 'ap-hermite';
        end
    end

    methods (Static)
        function obj = load(filename, options)
            %LOAD Load a versioned data-only state; callbacks are caller inputs.
            arguments
                filename {mustBeTextScalar}
                options.callbacks (1, 1) struct = struct()
            end
            record = mdemat.tools.io.readRecord(filename);
            obj = mdemat.application.plasma.ApHermiteState.fromRecord(record, callbacks=options.callbacks);
        end

        function obj = fromRecord(record, options)
            arguments
                record (1, 1) struct
                options.callbacks (1, 1) struct = struct()
            end
            record = mdemat.data.variable.validateStateRecord(record, 'ap-hermite');
            record = mdemat.data.variable.bindStateCallbacks(record, options.callbacks);
            if ~isempty(fieldnames(record.parameters))
                error('mdemat:persistence:InvalidState', 'Unexpected state parameters.');
            end
            [xDisc, vDisc] = mdemat.common.fem.SpatialState.restoreSpaces(record, true);
            obj = mdemat.application.plasma.ApHermiteState(xDisc, vDisc);
            obj.restoreRecordData(record);
        end
    end
end
