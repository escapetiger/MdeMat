classdef HermiteState < mdemat.common.fem.KineticState
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
                xDisc mdemat.common.fem.FiniteElementSpace
                vDisc mdemat.common.fem.SpectralSpace
            end

            obj@mdemat.common.fem.KineticState(xDisc, vDisc);
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

        function newObj = refine(obj, nLevels)
            % REFINE Create refined Hermite state.
            %
            %   newObj = refine(obj, nLevels) creates a refined Hermite
            %   state by refining the spatial mesh by 2^nLevels in each
            %   direction and doubling the Hermite mode count. Temperature
            %   and basis normalization are preserved; coefficients are not
            %   transferred, since a convergence run initializes each level.

            arguments
                obj mdemat.application.plasma.HermiteState
                nLevels{mustBeNonnegative, mustBeInteger}
            end

            newXDisc = obj.XDisc.refine(nLevels);
            oldVElement = obj.VDisc.Element;
            newNVDofs = oldVElement.NDofs * 2;
            newT0 = oldVElement.Approximator.Basis.Rhs.Temperature;
            newVElement = mdemat.common.fem.L2OrthotopeElement.hermite(1, newNVDofs, T = newT0);
            newVDisc = mdemat.common.fem.SpectralSpace(newVElement);
            newObj = mdemat.application.plasma.HermiteState(newXDisc, newVDisc);
        end

        function F = distribution(obj, xRef, v)
            % DISTRIBUTION Compute distribution function.
            %
            %   F = distribution(obj, xRef, v) computes the distribution
            %   function at reference points @a xRef with velocity @a v.

            arguments
                obj mdemat.application.plasma.HermiteState
                xRef{mustBeNumeric}
                v{mustBeNumeric}
            end

            D = obj.Dofs.D;
            D = obj.XDisc.eval(xRef, D); % (nx*nq, nh)
            F = obj.VDisc.eval(v, D.'); % (nx*nq, nv)
        end

        function obj = updateHistory(obj, t)
            % UPDATEHISTORY Integrate velocity moments and squared physical samples.
            arguments
                obj mdemat.application.plasma.HermiteState
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
            firstMoment = zeros(size(values, 1), 1);
            secondMode = firstMoment;
            if size(values, 2) >= 2
                firstMoment = values(:, 2);
            end
            if size(values, 2) >= 3
                secondMode = values(:, 3);
            end
            temperature = obj.Coefs.T0;
            field = space.eval(reference, obj.Dofs.E);
            kinetic = temperature/2*(weights.'*(values(:, 1)+sqrt(2)*secondMode));
            potential = weights.'*abs(field).^2/2;
            distance = values;
            distance(:, 1) = distance(:, 1)-1;
            squaredDistance = weights.'*abs(distance).^2;

            obj.History.time(end+1) = t;
            obj.History.mass(end+1) = weights.'*values(:, 1);
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
                obj mdemat.application.plasma.HermiteState
                options.omitCallbacks (1, 1) logical = false
            end
            record = snapshot@mdemat.common.fem.KineticState(obj, ...
                omitCallbacks=options.omitCallbacks);
            record.kind = 'hermite';
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
            obj = mdemat.application.plasma.HermiteState.fromRecord(record, callbacks=options.callbacks);
        end

        function obj = fromRecord(record, options)
            arguments
                record (1, 1) struct
                options.callbacks (1, 1) struct = struct()
            end
            record = mdemat.data.variable.validateStateRecord(record, 'hermite');
            record = mdemat.data.variable.bindStateCallbacks(record, options.callbacks);
            if ~isempty(fieldnames(record.parameters))
                error('mdemat:persistence:InvalidState', 'Unexpected state parameters.');
            end
            [xDisc, vDisc] = mdemat.common.fem.SpatialState.restoreSpaces(record, true);
            obj = mdemat.application.plasma.HermiteState(xDisc, vDisc);
            obj.restoreRecordData(record);
        end
    end
end
