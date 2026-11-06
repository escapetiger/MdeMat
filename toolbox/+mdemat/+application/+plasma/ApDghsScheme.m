classdef ApDghsScheme < mdemat.application.plasma.VlasovScheme
    % APDGHScheme Asymptotic preserving Discontinuous Galerkin Hermite
    % spectral method for Vlasov-Poisson-Fokker-Planck equations.
    %
    %   ApApDghsScheme implements an asymptotic preserving discontinuous
    %   Galerkin Hermite spectral method for solving
    %   Vlasov-Poisson-Fokker-Planck equations using discontinuous Galerkin
    %   spatial discretization and operator splitting temporal integration.
    %
    %   The scheme supports dealiasing to reduce aliasing errors in the
    %   nonlinear coupling terms using three filter options:
    %   1. Hou-Li high-order filter: preserves 12-15% more effective modes
    %   2. Traditional 2/3 dealiasing rule: sharp spectral cutoff
    %   3. Quasi time-consistent filter: preserves conservation laws and
    %      Galilean invariance while suppressing numerical recurrence
    %      (based on arxiv:1712.06433)
    %
    % See also:
    %   mdemat.application.plasma.VlasovScheme
    
    properties
        MultAmb % Assembly for multiplier
        xOpAAmb % Assembly for adjoint operators for A
        xOpBAmb % Assembly for adjoint operators for B
        SrcAmb % Assembly for source terms

        xOpA % adjoint operators for A
        xOpB % adjoint operators for B

        Poisson % Linear solver for the Poisson equation
        P
        A
        B
        C
        N
    end

    properties (Constant)
        Name = 'APDGHS' % Scheme name
    end
    
    methods
        function state = initialize(obj, state)
            % INITIALIZE Initialize the DG Hermite spectral scheme.
            %
            %   state = initialize(obj, state) sets up the DGHS scheme
            %   by assembling operators, fitting initial conditions,
            %   and preparing visualization components.
            
            arguments
                obj mdemat.application.plasma.ApDghsScheme
                state mdemat.application.plasma.ApHermiteState
            end
            
            %< Set up linear solver options
            linOpts = struct('nnzTh', inf, 'condTh', inf); % always use backslash
            obj.TDisc.reset(linOpts = linOpts);
            obj.Poisson = mdemat.common.linalg.LinearSolver(nnzTh=Inf, condTh=Inf);

            %< Assembly for weighted adjoint operators
            % obj.xOpAAmb = mdemat.common.fem.WeightedAdjointAssembly( ...
            %     state.XDisc, obj.Config.bcType, obj.Config.xPenaltyType);
            obj.xOpAAmb = mdemat.common.fem.AdjointAssembly( ...
                state.XDisc, obj.Config.bcType, obj.Config.xOpAType);
            obj.xOpBAmb = mdemat.common.fem.AdjointAssembly( ...
                state.XDisc, obj.Config.bcType, obj.Config.xOpBType);

            %< Assembly for multiplier
            obj.MultAmb = mdemat.common.fem.MultiplierAssembly( ...
                state.XDisc);
            
            %< Assembly for source
            obj.SrcAmb = mdemat.common.fem.SourceAssembly( ...
                state.XDisc);
            
            %< Set initial condition for Hermite coefficients
            state.setDof('D', state.XDisc.fit(obj.Config.ic, args={state.NVDofs}));

            %< Set reference temperature
            state.setCoefficient('T0', obj.Config.T0);

            %< Set ion density coefficient
            state.setCoefficient('rhoi', state.XDisc.fit(obj.Config.rhoi));

            %< Set electric potential equilibrium coefficient
            state.setCoefficient('phiInf', state.XDisc.fit(obj.Config.phiInf));

            %< Set density equilibrium coefficient
            state.setCoefficient('rhoInf', state.XDisc.fit(obj.Config.rhoInf));

            %< Set square root of density equilibrium coefficient
            state.setCoefficient('sqrtRhoInf', state.XDisc.fit(obj.Config.sqrtRhoInf));

            %< Set mass term
            obj.setMassTerm(state);
            
            %< Set transport term
            obj.setTransportTerm(state);

            %< Solve Poisson equation for initial potential
            obj.stepPoisson(state);
            
            %< Set initial history
            state.updateHistory(obj.TDisc.Timeline.Now);
        end
        
        function state = preStep(obj, state)
            % PRESTEP Preprocess operations before time step.
            %
            %   state = preStep(obj, state) performs preprocessing
            %   including time step size computation based on CFL condition
            %   or fixed time step size.
            
            arguments
                obj mdemat.application.plasma.ApDghsScheme
                state mdemat.application.plasma.ApHermiteState
            end
            
            if obj.hasConfig('cfl') && ~isempty(obj.Config.cfl)
                h = state.XDisc.Mesh.computeMeasure();
                obj.TDisc.setTimeStep(h = h, C = obj.Config.cfl, p = 1);
            end
            
            if obj.hasConfig('dt') && ~isempty(obj.Config.dt)
                obj.TDisc.setTimeStep(dt = obj.Config.dt);
            end
        end
        
        function state = step(obj, state)
            % STEP Perform one DG Hermite spectral time step.
            %
            %   state = step(obj, state) advances the solution by one time
            %   step using either Strang splitting or coupled IMEX scheme
            %   based on useSplitting configuration.

            arguments
                obj mdemat.application.plasma.ApDghsScheme
                state mdemat.application.plasma.ApHermiteState
            end

            dt = obj.TDisc.Timeline.StepSize;

            if obj.IsSplitting
                % Strang splitting scheme
                % state = obj.stepLinear(state, dt);
                % state = obj.stepLinearized(state, dt);
                % state = obj.stepLinearized(state, dt);
                % state = obj.stepNonlinear(state, dt);
                state = obj.stepLinearized(state, dt / 2);
                state = obj.stepNonlinear(state, dt);
                state = obj.stepLinearized(state, dt / 2);
            else
                % Coupled implicit scheme
                state = obj.stepCoupled(state, dt);
            end

            state = obj.applyDealiasingFilter(state);

            U = [state.Dofs.D, state.Dofs.Psi, state.Dofs.Q];
            obj.TDisc.add(U(:));
        end
        
        function state = postStep(obj, state)
            % POSTSTEP Postprocess operations after time step.
            %
            %   state = postStep(obj, state) performs postprocessing
            %   including visualization updates and timeline advancement.
            
            arguments
                obj mdemat.application.plasma.ApDghsScheme
                state mdemat.application.plasma.ApHermiteState
            end
            
            state.updateHistory(obj.TDisc.Timeline.Next);
            
            obj.TDisc.advance();
        end
        
        function state = finalize(obj, state)
            % FINALIZE Finalize simulation with scheme and state reporting.
            %
            %   state = finalize(obj, state) performs final reporting
            %   including scheme configuration details and state
            %   information summary.
            
            arguments
                obj mdemat.application.plasma.ApDghsScheme
                state mdemat.application.plasma.ApHermiteState
            end
            
            state.report();
            
            fprintf('[R] ODE Integrator: %s\n', class(obj.TDisc));
            fprintf('[R] Final Time: %.6g\n', obj.Config.tFinal);
            if obj.hasConfig('cfl') && ~isempty(obj.Config.cfl)
                fprintf('[R] CFL number: %.3f\n', obj.Config.cfl);
            end
            if obj.hasConfig('dt') && ~isempty(obj.Config.dt)
                fprintf('[R] Time step size: %.3f\n', obj.Config.dt);
            end
            fprintf('[R] Boundary Conditions: %s\n', obj.Config.bcType);
        end
       
    end
    
    methods
        function samples = samplePhysicalAssessment(obj,state,coarse,field,time)
            %SAMPLEPHYSICALASSESSMENT Physical F, including the AP equilibrium scale.
            samples = mdemat.application.plasma.assessmentSamples(obj,state,coarse,field,time);
        end
        function fields = samplePlotFields(obj, state, plan)
            fields = samplePlotFields@mdemat.internal.execution.MolScheme(obj, state, plan);
            physical = mdemat.application.plasma.plotDistribution(state, plan, ...
                obj.Observation.Plot.Components, obj.Config.vBBox);
            names = fieldnames(physical);
            for index = 1:numel(names), fields.(names{index}) = physical.(names{index}); end
        end
        function fields = samplePlotReferences(obj, state, plan, time)
            exacts = obj.Observation.Plot.Exacts;
            spatial = exacts;
            for name = {'F','Distance'}
                if isfield(spatial,name{1}),spatial=rmfield(spatial,name{1});end
            end
            fields = mdemat.internal.execution.plotReferenceFields(state.XDisc, spatial, ...
                obj.Observation.Plot.Components, plan, time);
            if isfield(exacts,'F') && isfield(obj.Observation.Plot.Components,'F')
                physical = mdemat.application.plasma.plotDistribution(state,plan, ...
                    struct('F',1),obj.Config.vBBox,reference=exacts.F,time=time);
                fields.F = physical.F;
            end
        end
    end

    methods (Access = protected)
        function obj = setup(obj)
            % SETUP Setup default configuration options for DGHS schemes.
            
            %< Call parent setup
            obj = setup@mdemat.application.plasma.VlasovScheme(obj);

            % AP potential/field retain their own physical field identities.
            obj.addConfig('components', default=struct('D', [], 'Psi', [], 'Q', [], 'F', 1), ...
                validator=@(x) isstruct(x) && all(structfun(@(y) isempty(y) || ...
                    isnumeric(y) && isscalar(y) && y > 0, x)));

            obj.addConfig('epsilon', default=1, ...
                validator=@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
            obj.addConfig('tau', default=Inf, ...
                validator=@(x) isnumeric(x) && isscalar(x) && x > 0);
            obj.addConfig('phiInf', default=@(x) zeros(1, size(x, 2)), ...
                validator=@(x) isa(x, 'function_handle'));
            obj.addConfig('rhoInf', default=@(x) ones(1, size(x, 2)), ...
                validator=@(x) isa(x, 'function_handle'));
            obj.addConfig('sqrtRhoInf', default=@(x) ones(1, size(x, 2)), ...
                validator=@(x) isa(x, 'function_handle'));
            obj.addConfig('EInf', default=@(x) zeros(1, size(x, 2)), ...
                validator=@(x) isa(x, 'function_handle'));
            obj.addConfig('ic', default=@(x, nh) [ones(1, size(x, 2)); zeros(nh-1, size(x, 2))], ...
                validator=@(x) isa(x, 'function_handle'));
            obj.addConfig('tSplitName', default='strang', ...
                validator=@mdemat.internal.execution.OdeScheme.validateSplitName);
            obj.addConfig('tOdeIntName', default={'sdirk2', 'sdirk2'}, ...
                validator=@mdemat.internal.execution.OdeScheme.validateOdeIntName);
            
            %< Add DGHS-specific required configuration options
            obj.addConfig('useSemiImplicit', default = true, ...
                validator = @(x) islogical(x) && isscalar(x));
            obj.addConfig('xBasisOrder', default = 1, ...
                validator = @(x) isnumeric(x) && isscalar(x) && (x > 0));
            obj.addConfig('xBasisType', default='nodal', ...
                validator=@(x) (ischar(x) || isstring(x)) && ismember(lower(x), {'nodal'}));
            obj.addConfig('xBasisPattern', default = 'Q', ...
                validator = @(x) (ischar(x) || isstring(x)) && ismember(upper(x), {'P', 'Q'}));
            obj.addConfig('xOpAType', default = {'', 'central'; 'central', ''}, ...
                validator = @(x) iscell(x) && isequal(size(x), [2, 2]));
            obj.addConfig('xOpBType', default = {'', 'right'; 'left', ''}, ...
                validator = @(x) iscell(x) && isequal(size(x), [2, 2]));
            obj.addConfig('vBBox', default = [-6, 6], ...
                validator = @(x) isnumeric(x) && isvector(x) && mod(length(x), 2) == 0);
            obj.addConfig('nh', default = 3, ...
                validator = @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 2 && fix(x) == x);
            obj.addConfig('useFilter', default = false, ...
                validator = @(x) islogical(x) && isscalar(x));
            obj.addConfig('dealiasingType', default = 'hou-li', ...
                validator = @(x) (ischar(x) || isstring(x)) && ismember(lower(x), {'hou-li', 'tc', 'quasi-tc'}));
            obj.addConfig('dealiasingBeta', default = 36, ...
                validator = @(x) isnumeric(x) && isscalar(x) && (x > 0));
            obj.addConfig('dealiasingGamma', default = 36, ...
                validator = @(x) isnumeric(x) && isscalar(x) && (x > 0));
            obj.addConfig('dealiasingCutoff', default = 2/3, ...
                validator = @(x) isnumeric(x) && isscalar(x) && (x > 0) && (x <= 1));
        end
        
        function obj = configure(obj)
            % CONFIGURE Generate DGHS-specific configuration dependencies.
            
            %< Call parent configure
            obj = configure@mdemat.application.plasma.VlasovScheme(obj);
            
            %< Determine spatial discretization identifier
            xId = sprintf('%s%d', upper(obj.Name), obj.Config.xBasisOrder);
            vId = sprintf('H%d', obj.Config.nh);
            obj.setDerived('xId', xId);
            obj.setDerived('vId', vId);
            obj.setDerived('sId', sprintf('%s%s', xId, vId));
            
            %< Set plot title
            eId = obj.Config.eId;
            obj.setDerived('titlePrefix', sprintf('%s with %s', eId, obj.Config.sId));
            
        end

    end
    
    methods (Access = private)
        function obj = setMassTerm(obj, state)
            % SETMASSTERM Set up mass term matrices.
            
            epsilon = obj.Config.epsilon;
            nh = state.NVDofs;
            ng = state.NXDofs;
                
            if obj.IsSplitting
                obj.M{1} = speye(ng * nh) * epsilon;
            else
                obj.M = speye(ng * nh) * epsilon;
            end
        end
        
        function obj = setTransportTerm(obj, state)
            % SETTRANSPORTTERM Set up transport term operators.

            T0 = obj.Config.T0;
            matD = obj.xOpAAmb.assembleMatrix(1);
            % sqrtRhoInf = state.Coefs.sqrtRhoInf;
            % rhs = 2 * T0 * matD{2, 1}{1} * sqrtRhoInf;            
            % f = @(x) reshape(state.XDisc.eval([], sqrtRhoInf), 1, []);
            % lhs = obj.MultAmb.assembleMatrix(f);
            % EInf = lhs \ rhs;
            EInf = obj.Config.EInf;
            matI = obj.MultAmb.assembleMatrix(EInf);
            obj.xOpA = cell(2, 2);
            obj.xOpA{1, 2} = -sqrt(T0) * matD{1, 2}{1} - matI ./ (2*sqrt(T0));
            obj.xOpA{2, 1} = +sqrt(T0) * matD{2, 1}{1} - matI ./ (2*sqrt(T0));
            obj.xOpA{1, 1} = +sqrt(T0) * matD{1, 1}{1};

            matD = obj.xOpBAmb.assembleMatrix(1);
            % sqrtRhoInf = state.Coefs.sqrtRhoInf;
            % rhs = 2 * T0 * matD{2, 1}{1} * sqrtRhoInf;
            % f = @(x) reshape(state.XDisc.eval([], sqrtRhoInf), 1, []);
            % lhs = obj.MultAmb.assembleMatrix(f);
            % EInf = lhs \ rhs;
            EInf = obj.Config.EInf;
            matI = obj.MultAmb.assembleMatrix(EInf);
            obj.xOpB = cell(2, 2);
            obj.xOpB{1, 2} = -sqrt(T0) * matD{1, 2}{1} - matI ./ (2*sqrt(T0));
            obj.xOpB{2, 1} = +sqrt(T0) * matD{2, 1}{1} - matI ./ (2*sqrt(T0));
            obj.xOpB{1, 1} = -sqrt(T0) * matD{1, 1}{1};
            
            %< Prestored weighted operators
            % T0 = obj.Config.T0;
            % weight = obj.Config.sqrtRhoInf;
            % dWeight = obj.Config.dSqrtRhoInf;
            % D = obj.xOpAAmb.assembleMatrix(sqrt(T0), weight, dWeight=dWeight);
            % obj.xOpA = cell(2, 2);
            % obj.xOpA{1, 2} = D{1, 2}{1};
            % obj.xOpA{2, 1} = D{2, 1}{1};
            % obj.xOpA{1, 1} = D{1, 1}{1};

            %< Linear stage: free transport + collision
            obj.addLinearizedTransport(state);
        end

        function obj = addLinearizedTransport(obj, state)
            % ADDLINEARIZEDTRANSPORT Add linearized transport terms.
            
            nh = state.NVDofs;
            ng = state.NXDofs;
            
            if obj.IsSplitting
                obj.L{1} = repmat({sparse(ng, ng)}, nh, nh);
                obj.P{1} = [
                    repmat({sparse(ng, ng)}, nh, 2), repmat({sparse(ng, 1)}, nh, 1)];
                obj.A{1} = [
                    repmat({sparse(ng, ng)}, 2, nh); 
                    repmat({sparse(1, ng)}, 1, nh)];
                obj.B{1} = [
                    repmat({sparse(ng, ng)}, 2, 2), repmat({sparse(ng, 1)}, 2, 1);
                    repmat({sparse(1, ng)}, 1, 2), {sparse(1, 1)}];
                obj.C{1} = [
                    repmat({sparse(ng, 1)}, 2, 1); 
                    {sparse(1, 1)}];
                obj.addFreeTransport(state);
                obj.addPoissonCoulping(state);
                obj.addCollision(state);
                obj.L{1} = mdemat.common.linalg.block(obj.L{1});
                obj.A{1} = mdemat.common.linalg.block(obj.A{1});
                obj.B{1} = mdemat.common.linalg.block(obj.B{1});
                obj.C{1} = mdemat.common.linalg.block(obj.C{1});
                obj.P{1} = mdemat.common.linalg.block(obj.P{1});
            else
                obj.L = repmat({sparse(ng, ng)}, nh, nh);
                obj.P = [
                    repmat({sparse(ng, ng)}, nh, 2), repmat({sparse(ng, 1)}, nh, 1)];
                obj.A = [
                    repmat({sparse(ng, ng)}, 2, nh);
                    repmat({sparse(1, ng)}, 1, nh)];
                obj.B = [
                    repmat({sparse(ng, ng)}, 2, 2), repmat({sparse(ng, 1)}, 2, 1);
                    repmat({sparse(1, ng)}, 1, 2), {sparse(1, 1)}];
                obj.C = [
                    repmat({sparse(ng, 1)}, 2, 1);
                    {sparse(1, 1)}];
                obj.addFreeTransport(state);
                obj.addPoissonCoulping(state);
                obj.addCollision(state);
                obj.L = mdemat.common.linalg.block(obj.L);
                obj.A = mdemat.common.linalg.block(obj.A);
                obj.B = mdemat.common.linalg.block(obj.B);
                obj.C = mdemat.common.linalg.block(obj.C);
                obj.P = mdemat.common.linalg.block(obj.P);
            end
        end

        function obj = addFreeTransport(obj, state)
            % ADDFREETRANSPORT Add free transport terms.
            
            nh = state.NVDofs;
            
            if obj.IsSplitting
                obj.L{1}{1, 2} = obj.xOpA{1, 2};
                for k = 2:nh - 1
                    obj.L{1}{k, k - 1} = -sqrt(k-1) * obj.xOpA{2, 1};
                    obj.L{1}{k, k + 1} = sqrt(k) * obj.xOpA{1, 2};
                end
                obj.L{1}{nh, nh - 1} = -sqrt(nh - 1) * obj.xOpA{2, 1};
            else
                obj.L{1, 2} = obj.xOpA{1, 2};
                for k = 2:nh - 1
                    obj.L{k, k - 1} = -sqrt(k-1) * obj.xOpA{2, 1};
                    obj.L{k, k + 1} = sqrt(k) * obj.xOpA{1, 2};
                end
                obj.L{nh, nh - 1} = -sqrt(nh - 1) * obj.xOpA{2, 1};
            end
        end

        function obj = addPoissonCoulping(obj, state)
            % ADDELECTRIC Add linearized Poisson coupling term.
            ng = state.NXDofs;
            rhoInf = state.Coefs.rhoInf;
            f = @(x) reshape(state.XDisc.eval([], rhoInf), 1, []);
            T = obj.MultAmb.assembleMatrix(f);
            [cv, ch] = obj.poissonGauge(state);

            if obj.IsSplitting
                obj.P{1}{2, 2} = T;
                obj.A{1}{1, 1} = -speye(ng);
                obj.B{1}{1, 1} = obj.xOpB{1, 1};
                obj.B{1}{1, 2} = -obj.xOpB{1, 2};
                obj.B{1}{2, 1} = obj.xOpB{2, 1};
                obj.B{1}{2, 2} = T;
                obj.B{1}{1, 3} = cv;
                obj.B{1}{3, 1} = ch.';
                obj.C{1}{1} = -state.Coefs.sqrtRhoInf;
            else
                obj.P{2, 2} = T;
                obj.A{1, 1} = -speye(ng);
                obj.B{1, 1} = obj.xOpB{1, 1};
                obj.B{1, 2} = -obj.xOpB{1, 2};
                obj.B{2, 1} = obj.xOpB{2, 1};
                obj.B{2, 2} = T;
                obj.B{1, 3} = cv;
                obj.B{3, 1} = ch.';
                obj.C{1} = -state.Coefs.sqrtRhoInf;
            end
        end
        
        function obj = addCollision(obj, state)
            % ADDCOLLISION Add collision terms.

            epsilon = obj.Config.epsilon;
            tau = obj.Config.tau;
            nh = state.NVDofs;
            T = obj.MultAmb.assembleMatrix(1);
            if obj.IsSplitting
                for k = 2:nh
                    obj.L{1}{k, k} = -(k - 1) * T * epsilon / tau;
                end
            else
                for k = 2:nh
                    obj.L{k, k} = -(k - 1) * T * epsilon / tau;
                end
            end
        end
       
        function obj = setPoissonLhs(obj, state)
            % SETPOISSONLHS Set left-hand side matrix for Poisson equation.

            if obj.Poisson.IsPrecomputed
                return;
            end

            ng = state.NXDofs;
            rhoInf = state.Coefs.rhoInf;
            f = @(x) reshape(state.XDisc.eval([], rhoInf), 1, []);
            T = obj.MultAmb.assembleMatrix(f);
            
            [cv, ch] = obj.poissonGauge(state);
            lhs = [
                repmat({sparse(ng, ng)}, 2, 2), repmat({sparse(ng, 1)}, 2, 1);
                repmat({sparse(1, ng)}, 1, 2), {sparse(1, 1)}];
            lhs{1, 1} = obj.xOpB{1, 1};
            lhs{1, 2} = -obj.xOpB{1, 2};
            lhs{2, 1} = obj.xOpB{2, 1};
            lhs{2, 2} = T;
            lhs{1, 3} = cv;
            lhs{3, 1} = ch.';
            lhs = mdemat.common.linalg.block(lhs);
            obj.Poisson.setLhs(lhs);
            obj.Poisson.precompute();
        end
        
        function obj = setPoissonRhs(obj, state)
            % SETPOISSONRHS Set right-hand side vector for Poisson equation.

            ng = state.XDisc.NGlobalDofs;
            rhs = repmat({sparse(ng, 1)}, 2, 1);
            rhs{1} = state.Dofs.D(:, 1) - state.Coefs.sqrtRhoInf;
            rhs = vertcat(rhs{:});
            rhs = [rhs; 0];
            obj.Poisson.setRhs(rhs);
        end
        
        function state = stepPoisson(obj, state)
            % STEPPOISSON Solve Poisson equation for potential and its
            % derivative.

            obj.setPoissonLhs(state);
            obj.setPoissonRhs(state);
            W = obj.Poisson.solve();
            W = reshape(W(1:end-1), [], 2);
            state.Dofs.Psi = W(:, 1);
            state.Dofs.Q = W(:, 2);
        end
        
        function state = stepLinear(obj, state, dt)
            % STEPLINEAR Perform linear transport and collision step.

            nh = state.NVDofs;
            W0 = state.Dofs.D(:);
            
            if obj.IsSplitting
                tDisc = obj.TDisc.getFactor(1);
                tDisc.update(W0);
                W = tDisc.step(dt=dt, L=obj.L{1}, M=obj.M{1});
            else
                tDisc = obj.TDisc;
                tDisc.update(W0);
                W = tDisc.step(dt=dt, ...
                    L = mdemat.common.linalg.block(obj.L), ...
                    M = mdemat.common.linalg.block(obj.M));
            end
            state.Dofs.D = reshape(W, [], nh);
        end

        function state = stepLinearized(obj, state, dt)
            % STEPLINEAR Perform linearized step.

            nh = state.NVDofs;
            W0 = [state.Dofs.D(:); state.Dofs.Psi(:); state.Dofs.Q(:); 0];
            
            if obj.IsSplitting
                tDisc = obj.TDisc.getFactor(1);
                tDisc.update(W0);
                W = tDisc.step(dt=dt, ...
                    L=obj.L{1}, M=obj.M{1}, ...
                    P = obj.P{1}, A = obj.A{1}, B = obj.B{1}, C = obj.C{1});
            else
                tDisc = obj.TDisc;
                tDisc.update(W0);
                W = tDisc.step(dt=dt, ...
                    L = mdemat.common.linalg.block(obj.L), ...
                    M = obj.M, P = obj.P, A = obj.A, B = obj.B, C = obj.C);
            end
            W = reshape(W(1:end-1), [], nh+2);
            state.Dofs.D = W(:, 1:nh);
            state.Dofs.Psi = W(:, nh+1);
            state.Dofs.Q = W(:, nh+2);
        end
        
        function state = stepNonlinear(obj, state, dt)
            % STEPNONLINEAR Perform nonlinear electric field coupling step.

            tDisc = obj.TDisc.getFactor(2);
            switch class(tDisc)
                case 'mdemat.common.ode.BeIntegrator'
                    state = obj.stepNonlinearBe(state, dt);
                case 'mdemat.common.ode.Sdirk2Integrator'
                    state = obj.stepNonlinearSdirk2(state, dt);
                otherwise
                    assert(all(0, 'all'), 'mdemat:application:plasma:ApDghsScheme:InvalidOdeInt', ...
                        '%s', 'Unsupported ODE integrator for nonlinear step.');
            end
        end
        
    end

    methods (Access = protected)
        function [constant, weightedMean] = poissonGauge(obj, state)
            % POISSONGAUGE Physical mean of Psi/sqrt(rhoInf), not a
            % Euclidean inner product of mass-normalized coefficients.
            constant = obj.SrcAmb.assembleVector(@(x) ones(1, size(x, 2)));
            multiplier = obj.MultAmb.assembleMatrix( ...
                @(x) 1 ./ obj.Config.sqrtRhoInf(x));
            weightedMean = multiplier.' * state.XDisc.integrationWeights();
        end

        function state = stepNonlinearBe(obj, state, dt)
            % STEPNONLINEARBE Perform nonlinear electric field coupling
            % step using backward Euler method.

            nh = state.NVDofs;
            epsilon = obj.Config.epsilon;

            D = state.Dofs.D;
            Q = state.Dofs.Q;
            sqrtRhoInf = state.Coefs.sqrtRhoInf;
            f = @(x) reshape(state.XDisc.eval([], sqrtRhoInf.*Q), 1, []);
            T = obj.MultAmb.assembleMatrix(f);
            DNew = zeros(size(D));
            DNew(:, 1) = D(:, 1);
            for k = 2:nh
                previous = DNew(:, k-1);
                if k == 2
                    previous = previous - sqrtRhoInf;
                end
                DNew(:, k) = D(:, k) + dt * sqrt(k-1) ...
                    * (T * previous) / epsilon;
            end
            state.Dofs.D = DNew;
        end

        function state = stepNonlinearSdirk2(obj, state, dt)
            % STEPNONLINEARSDIRK2 Perform nonlinear step using second-order
            % SDIRK method.

            nh = state.NVDofs;
            epsilon = obj.Config.epsilon;

            D = state.Dofs.D;
            Q = state.Dofs.Q;
            sqrtRhoInf = state.Coefs.sqrtRhoInf;
            DInf = zeros(size(D));
            DInf(:, 1) = sqrtRhoInf;
            f = @(x) reshape(state.XDisc.eval([], sqrtRhoInf.*Q), 1, []);
            T = obj.MultAmb.assembleMatrix(f);

            %< First stage
            gamma = 1 - sqrt(2)/2;
            a11 = gamma;
            F1 = zeros(size(D));
            D1 = zeros(size(D));
            D1(:, 1) = D(:, 1);
            for k = 2:nh
                F1(:, k) = sqrt(k-1) .* (T * (D1(:, k-1) - DInf(:, k-1))) / epsilon;
                D1(:, k) = D(:, k) + a11 * dt * F1(:, k);
            end

            %< Second stage
            a21 = 1-gamma; a22 = gamma;
            F2 = zeros(size(D));
            D2 = zeros(size(D));
            D2(:, 1) = D(:, 1);
            for k = 2:nh
                F2(:, k) = sqrt(k-1) .* (T * (D2(:, k-1) - DInf(:, k-1))) / epsilon;
                D2(:, k) = D(:, k) + a21 * dt * F1(:, k) + a22 * dt * F2(:, k);
            end

            %< Final stage
            b1 = 1-gamma; b2 = gamma;
            DNew = D + dt * (b1 * F1 + b2 * F2);
            state.Dofs.D = DNew;
        end

    end

    methods (Access = private)

        function state = stepCoupled(obj, state, dt)
            % STEPCOUPLED Perform coupled step with fixed-point iteration.
            %
            %   state = stepCoupled(obj, state, dt) performs a coupled
            %   time step using fixed-point iteration to solve the nonlinear
            %   system arising from implicit coupling. Only computes the
            %   time step integration once per iteration.

            if obj.IsSplitting, return; end

            nh = state.NVDofs;
            W0 = [state.Dofs.D(:); state.Dofs.Psi(:); state.Dofs.Q(:); 0];
            tDisc = obj.TDisc;
            tDisc.update(W0);

            if obj.IsFullyImplicit && obj.Config.useSemiImplicit
                obj.updatePoissonCoupling(state);
                
                W = tDisc.step(dt=dt, ...
                    L = obj.L + obj.N, S = obj.S, ...
                    M = obj.M, P = obj.P, A = obj.A, B = obj.B, C = obj.C);    
            elseif obj.IsImex
                obj.updatePoissonCoupling(state);
                
                FFn = @(W) obj.assembleF(state, W);

                W = tDisc.step(dt=dt, ...
                    L = obj.L, F = FFn, ...
                    M = obj.M, P = obj.P, A = obj.A, B = obj.B, C = obj.C);    
            else
                LFn = @(W) obj.assembleL(state, W);
                SFn = @(W) obj.assembleS(state, W);

                W = tDisc.step(dt=dt, ...
                    L = LFn, S = SFn, ...
                    M = obj.M, P = obj.P, A = obj.A, B = obj.B, C = obj.C);
            end

            W = reshape(W(1:end-1), [], nh+2);
            state.Dofs.D = W(:, 1:nh);
            state.Dofs.Psi = W(:, nh+1);
            state.Dofs.Q = W(:, nh+2);
        end

        function obj = updatePoissonCoupling(obj, state)
            nh = state.NVDofs;
            ng = state.NXDofs;

            Q = state.Dofs.Q;
            
            obj.N = repmat({sparse(ng, ng)}, nh, nh);
            sqrtRhoInf = state.Coefs.sqrtRhoInf;
            Tn = obj.MultAmb.assembleMatrix(sqrtRhoInf .* Q);
            for k = 2:nh
                obj.N{k, k-1} = sqrt(k-1) * Tn;
            end
            obj.N = mdemat.common.linalg.block(obj.N);

            obj.S = repmat({sparse(ng, 1)}, nh, 1);
            rhoInf = state.Coefs.rhoInf;
            obj.S{2} = -rhoInf .* Q;
            obj.S = mdemat.common.linalg.block(obj.S);
        end

        function L = assembleL(obj, state, W)
            nh = state.NVDofs;
            ng = state.NXDofs;

            Q = W((ng*(nh+1)+1):(ng*(nh+2)));

            obj.N = repmat({sparse(ng, ng)}, nh, nh);
            sqrtRhoInf = state.Coefs.sqrtRhoInf;
            Tn = obj.MultAmb.assembleMatrix(sqrtRhoInf .* Q);
            for k = 2:nh
                obj.N{k, k-1} = sqrt(k-1) * Tn;
            end
            L = obj.L + mdemat.common.linalg.block(obj.N);
        end

        function S = assembleS(obj, state, W)
            nh = state.NVDofs;
            ng = state.NXDofs;

            obj.S = repmat({sparse(ng, 1)}, nh, 1);
            Q = W((ng*(nh+1)+1):(ng*(nh+2)));
            rhoInf = state.Coefs.rhoInf;
            obj.S{2} = -rhoInf .* Q;
            S = mdemat.common.linalg.block(obj.S);
        end

        function F = assembleF(obj, state, W)
            nh = state.NVDofs;
            ng = state.NXDofs;

            D = reshape(W(1:ng*nh), [], nh);
            Q = state.Dofs.Q;

            obj.F = repmat({sparse(ng, 1)}, nh, 1);
            sqrtRhoInf = state.Coefs.sqrtRhoInf;
            rhoInf = state.Coefs.rhoInf;
            Tn = obj.MultAmb.assembleMatrix(sqrtRhoInf .* Q);
            for k = 2:nh
                obj.F{k} = sqrt(k-1) * Tn * D(:, k-1);
            end
            obj.F{2} = obj.F{2} - rhoInf .* Q;
            F = mdemat.common.linalg.block(obj.F);
        end

        function state = applyDealiasingFilter(obj, state)
            % APPLYDEALIASINGFILTER Apply dealiasing filter to Hermite
            % coefficients.
            %
            %   state = applyDealiasingFilter(obj, state) applies the specified
            %   dealiasing filter to the Hermite coefficient matrix in state to reduce
            %   aliasing errors in spectral computations.

            if ~obj.Config.useFilter
                return;
            end

            state.Dofs.D = obj.filter(state.Dofs.D);
        end
    
        function D = filter(obj, D, k)
            % FILTER Filter function.

            nh = size(D, 2);
            
            if nargin < 3, k = 1:nh; end

            dt = obj.TDisc.Timeline.StepSize;
            epsilon = obj.Config.epsilon;
            beta = obj.Config.dealiasingBeta;
            gamma = obj.Config.dealiasingGamma;
            c = obj.Config.dealiasingCutoff;
            k0 = k - 1;
            eta = k0 / nh;

            switch lower(obj.Config.dealiasingType)
                case 'hou-li'
                    g = @(eta, dt) 1;
                case 'tc'
                    g = @(eta, dt) epsilon * dt;
                case 'quasi-tc'
                    g = @(eta, dt) (epsilon * dt).^(1-eta.^gamma);
                otherwise
                    assert(all(0, 'all'), 'mdemat:application:plasma:ApDghsScheme:InvalidDealiasingType', ...
                        '%s', ['Invalid dealiasing type: ', obj.Config.dealiasingType]);
            end

            sigma = (eta <= c) .* 1 + (eta > c) .* exp(-beta * (eta).^gamma) .* g(eta, dt);
            D(:, k) = D(:, k) .* sigma;
        end
    end
end
