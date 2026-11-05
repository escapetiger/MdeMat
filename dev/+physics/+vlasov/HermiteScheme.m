classdef HermiteScheme < physics.Scheme
    properties
        constant
        auxiliary
        primal
        diffusion
        M
        L
        F
        S
        Ah
        Bh
    end

    methods
        function state = initialize(obj, state)
            % INITIALIZE Initialize the macro-micro DG scheme.

            %< Reset timer
            obj.timer.reset();

            if isempty(obj.config.bc)
                bcType = 'periodic';
            else
                bcType = 'dirichlet';
            end

            %< Set auxiliary operator
            obj.auxiliary = approx.assembly.AuxiliaryDivergenceOperator( ...
                state.xDisc.space, state.xDisc.operator, bcType, ...
                obj.config.xAuxFluxType, obj.config.xAuxBoundaryJumpType);

            %< Set primal operator
            obj.primal = approx.assembly.PrimalDivergenceOperator( ...
                state.xDisc.space, state.xDisc.operator, bcType, ...
                obj.config.xPrmFluxType);

            obj.diffusion = approx.assembly.DiffusionOperator( ...
                state.xDisc.space, state.xDisc.operator, bcType, ...
                obj.config.xAuxFluxType, ...
                obj.config.xPrmFluxType, ...
                obj.config.xAuxBoundaryJumpType);

            %< Set constant operator
            obj.constant = approx.assembly.ConstantOperator( ...
                state.xDisc.space, state.xDisc.operator);

            %< Set initial condition
            state.dofs.D = state.xDisc.space.project(obj.config.ic.D);
            state.dofs.omega = state.xDisc.space.project(obj.config.ic.omega);

            %< Set coefficients
            state.coefs.rhoInf = state.xDisc.space.average(obj.config.rhoInf);
            state.coefs.EInf = state.xDisc.space.average(obj.config.EInf);

            %< Set mass term
            obj.setMassTerm(state);

            %< Set linear term
            obj.setLinearTerm(state);

            %< Set source term
            obj.setSourceTerm(state);

            %< Set nonlinear term
            obj.F = cell(1, 2);
           
            %< Reset time discretization
            obj.tDisc.reset();

            %< Reset visualizer
            if ~isempty(obj.visualizer) && obj.visualizer.isEnabled
                obj.visualizer.reset();

                obj.timer.start(2, '[M] Visualization.');
                obj.visualizer.plot(state.xDisc.space, state.dofs, obj.tDisc);
                obj.timer.stop(2);
            end
        end

        function state = preStep(obj, state)
            % PRESTEP Preprocess operations before time step.

            %< Set time step
            h = state.xDisc.mesh.measure;
            obj.tDisc.timeline.setTimeStep(h, obj.config.cfl);
            if isa(obj.tDisc, 'approx.odeint.BdfIntegrator')
                obj.tDisc.setCoefficients();
            end
        end

        function state = step(obj, state)
            % STEP Perform one upwind DG time step.

            nh = state.nModes;

            %< Update history
            W = [state.dofs.D(:); state.dofs.omega(:)];
            obj.tDisc.update(W);

            %< Integration step
            W = obj.tDisc.step(obj.L, obj.F, obj.S, obj.M);
            
            W = reshape(W, [], nh + 2);
            state.dofs.D = W(:, 1:nh+1);
            state.dofs.omega = W(:, nh+2);
        end

        function state = postStep(obj, state)
            % POSTSTEP Postprocess operations after time step.

            if ~isempty(obj.visualizer) && obj.visualizer.isEnabled
                obj.timer.start(2, '[M] Visualization.');
                obj.visualizer.plot(state.xDisc.space, state.dofs, obj.tDisc);
                obj.timer.stop(2);
            end

            obj.tDisc.timeline.advance();
        end

        function state = finalize(obj, varargin)
            % FINALIZE Finalize simulation with error analysis.
            if length(varargin) == 1
                state = varargin{1};
                coarse = [];
            else
                [state, coarse] = varargin{:};
            end

            if ~isempty(obj.analyzer) && obj.analyzer.isEnabled
                obj.timer.start(3, '[M] Error Evaluation.');
                state.clean();
                if isempty(coarse)
                    core.except.assert(0, 'InvalidInput', ...
                        'Do not use absolute error.');
                else
                    obj.analyzer.setLevel(state.xDisc.mesh);
                    obj.analyzer.relative(state.xDisc.space, state.dofs, ...
                        coarse.xDisc.space, coarse.dofs);
                end
                obj.timer.stop(3);
            end

            obj.timer.report();
        end
    end

    methods (Access = private)
        function obj = setMassTerm(obj, state)
            nh = state.nModes;
            ng = state.xDisc.space.nGlobalDofs;
            obj.M = cell(nh+2, nh+2);
            for i = 1:nh+1
                obj.M{i, i} = speye(ng, ng) * obj.config.epsilon;
            end
            obj.M{nh+2, nh+2} = sparse(ng, ng);
            obj.M = core.linalg.block(obj.M);
        end

        function obj = setLinearTerm(obj, state)
            nh = state.nModes;
            ng = state.xDisc.space.nGlobalDofs;

            obj.L = cell(1, 2);
            obj.L{1} = repmat({sparse(ng, ng)}, nh+2, nh+2);
            obj.L{2} = @(U) computeL2(U);

            obj.assembleDualOperators(state);

            obj.addLinearizedFreeTransport(state);

            obj.addCollision(state);

            obj.L{1} = core.linalg.block(obj.L{1});

            function L = computeL2(U)
                ng = state.xDisc.space.nGlobalDofs;

                omega = reshape(U, ng, []);
                omega = omega(:, nh+2);
                eta = obj.Ah * omega;
                f = @(x) state.xDisc.space.evaluate([], x, eta);
                eta = state.xDisc.space.average(f);
                
                L = repmat({sparse(ng, ng)}, nh+2, nh+2);
                C = obj.constant.assemble(eta ./ sqrt(state.coefs.rhoInf));
                
                for k = 2:nh + 1
                    L{k, k-1} = sqrt(k) * C;
                end

                L{nh+2, nh + 2} = -speye(ng, ng);
                L = core.linalg.block(L);
            end
        end

        function obj = assembleDualOperators(obj, state)
            T0 = obj.config.T0;
            D1 = obj.auxiliary.linear(1);
            D2 = obj.primal.linear(1);
            C1 = obj.constant.assemble(state.coefs.EInf);
            obj.Ah = (sqrt(T0) * D1 - C1 / (2 * sqrt(T0)));
            obj.Bh = (-sqrt(T0) * D2 - C1 / (2 * sqrt(T0)));
        end

        function obj = addLinearizedFreeTransport(obj, state)
            nh = state.nModes;
            for k = 1:nh + 1
                if k > 1
                    obj.L{1}{k, k-1} = -obj.Ah;
                end
                if k < nh + 1
                    obj.L{1}{k, k+1} = obj.Bh;
                end
            end
            obj.L{1}{2, nh+2} = -obj.Ah;
            I = obj.constant.assemble(1);
            obj.L{1}{nh+2, 1} = I;

            ng = state.xDisc.space.nGlobalDofs;
            h = state.xDisc.space.mesh.spacings{1};
            flux1 = approx.assembly.FluxAssembly( ...
                state.xDisc.space, state.xDisc.operator, 'periodic', 'right');
            flux2 = approx.assembly.FluxAssembly( ...
                state.xDisc.space, state.xDisc.operator, 'periodic', 'left');
            coe = obj.auxiliary.volume.scaleConstant(1, 1);
            T1 = flux1.assembleFluxPartial(1, coe);
            T2 = flux2.assembleFluxPartial(1, coe);
            D3 = core.linalg.sparseFromTriplets([-T1; T2], ng, ng);
            C2 = obj.constant.assemble(1./state.coefs.rhoInf);
            obj.L{1}{nh+2, nh+2} = -(obj.Bh * C2 * (obj.Ah - D3));
        end

        function obj = addCollision(obj, state)
            nh = state.nModes;
            tau0 = obj.config.tau0;
            epsilon = obj.config.epsilon;
            for k = 2:nh + 1
                C = obj.constant.assemble(k-1);
                obj.L{1}{k, k} = -C / (tau0 * epsilon);
            end
        end
        
        function obj = setSourceTerm(obj, state)
            nh = state.nModes;
            ng = state.xDisc.space.nGlobalDofs;
            obj.S{1} = repmat({sparse(ng, 1)}, nh+2, 1);
            obj.S{1}{nh + 2} = -diag(obj.constant.assemble(sqrt(state.coefs.rhoInf)));
            obj.S{1} = vertcat(obj.S{1}{:});

            obj.S{2} = [];
            obj.S{2} = @(U) computeS2(U);

            function S = computeS2(U)
                omega = reshape(U, ng, []);
                omega = omega(:, nh+2);
                eta = obj.Ah * omega;
                S = repmat({sparse(ng, 1)}, nh+2, 1);
                S{2} = eta;
                S{nh+2} = omega;
                S = vertcat(S{:});
            end
        end
    end
end