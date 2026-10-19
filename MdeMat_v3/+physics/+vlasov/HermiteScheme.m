classdef HermiteScheme < physics.vlasov.VlasovScheme
    properties
        multiplier
        adjoint
        source

        dual
        dualw
        laplace

        P
        A
        B
        C
    end

    methods
        function obj = plotDistribution(obj, state)
            visualizer = obj.config.visualizer;
            xRef = visualizer.database.refNodes;
            X = state.xDisc.mesh.collocate(xRef);
            [F, V] = state.distribution(xRef, obj.config.T0, obj.config.sqrtRhoInf);

            figure(2);
            ax = gca;
            cla(ax);
            hold on;
            imagesc(ax, X, V, F.');
            xlabel(ax, 'x');
            ylabel(ax, 'v');

            axis(ax, 'xy', 'equal', 'tight');

            colormap('turbo');
            colorbar;

            uMin = min(F(:));
            uMax = max(F(:));

            %< Apply padding to data range
            padding = 0.05 * (uMax - uMin);
            if padding == 0
                %< Handle constant data case
                padding = 0.05 * abs(uMin);
                if padding == 0
                    padding = 0.1; %< Minimum padding for zero data
                end
            end
            cLimits = [uMin - padding, uMax + padding];

            caxis(ax, cLimits);
            hold off;
            title(sprintf('t = %.2f', obj.tDisc.timeline.next));
        end

        function state = initialize(obj, state)
            % INITIALIZE Initialize the macro-micro DG scheme.

            %< Assembly for adjoint operators
            obj.adjoint = approx.assembly.FiniteElementAdjointAssembly( ...
                state.xDisc, obj.config.bcType, obj.config.xPenaltyType);

            %< Assembly for multipliers
            obj.multiplier = approx.assembly.FiniteElementMultiplierAssembly( ...
                state.xDisc);

            %< Assembly for source
            obj.source = approx.assembly.FiniteElementSourceAssembly(state.xDisc);

            %< Set initial condition
            state.addDof('D', state.xDisc.fit([], obj.config.ic.D));
            state.addDof('omega', state.xDisc.fit([], obj.config.ic.omega));

            %< Set coefficients
            state.addCoefficient('rhoInf', state.xDisc.mean([], obj.config.rhoInf));
            state.addCoefficient('EInf', state.xDisc.mean([], obj.config.EInf));

            %< Set mass term
            obj.setMassTerm(state);

            %< Set linear term
            obj.setLinearTerm(state);

            %< Set source term
            obj.setSourceTerm(state);

            state.addHistory(obj.tDisc.timeline.now, obj.config.T0, obj.config.rhoInf);

            obj.plotDistribution(state);
        end

        function state = preStep(obj, state)
            % PRESTEP Preprocess operations before time step.

            if isfield(obj.config, 'cfl')
                h = state.xDisc.mesh.computeMeasure();
                obj.tDisc.setTimeStep(h, obj.config.cfl, []);
            end

            if isfield(obj.config, 'dt')
                obj.tDisc.setTimeStep(obj.config.dt, 1, []);
            end
        end

        function state = step(obj, state)
            % STEP Perform one upwind DG time step.

            nh = state.nModes - 1;

            %< Update history
            W0 = [state.dofs.D(:); state.dofs.omega(:); 0];

            %< VERSION 1
            obj.tDisc.update(W0);

            W = obj.tDisc.stepA(obj.L, obj.F, obj.S, obj.M, obj.P, ...
                obj.A, obj.B, obj.C);

            obj.tDisc.add(W0(1:end-1));

            %< VRESION 2

%             tDisc1 = obj.tDisc.getFactor(1);
%             tDisc1.update(W0);
% 
%             W = tDisc1.stepA(obj.L{1}, obj.S{1}, obj.M{1}, obj.P{1}, ...
%                 obj.A{1}, obj.B{1}, obj.C{1});
% 
%             tDisc1.add(W(1:end-1));

            W = reshape(W(1:end-1), [], nh+2);
            state.dofs.D = W(:, 1:nh+1);
            state.dofs.omega = W(:, nh+2);
        end

        function state = postStep(obj, state)
            % POSTSTEP Postprocess operations after time step.

            state.addHistory(obj.tDisc.timeline.next, obj.config.T0, obj.config.rhoInf);

            obj.plotDistribution(state);

            obj.tDisc.advance();
        end

        function state = finalize(obj, state)
            % FINALIZE Finalize simulation with error analysis.

            %< Report state information
            state.report();

            fprintf('[R] ODE Integrator: %s\n', class(obj.tDisc));
            fprintf('[R] Final Time: %.6g\n', obj.config.tFinal);
            if isfield(obj.config, 'cfl')
                fprintf('[R] CFL number: %.3f\n', obj.config.cfl);
            end
            if isfield(obj.config, 'dt')
                fprintf('[R] Time step size: %.3f\n', obj.config.dt);
            end
            %< Report boundary condition type
            if isempty(obj.config.bc)
                fprintf('[R] Boundary Conditions: Periodic\n');
            else
                fprintf('[R] Boundary Conditions: Dirichlet\n');
            end
        end

    end

    methods (Access = private)
        function obj = setMassTerm(obj, state)
            epsilon = obj.config.epsilon;
            tau = obj.config.tau;
            nh = state.nModes - 1;
            ng = state.xDisc.nGlobalDofs;
            obj.M{1} = cell(nh+1, nh+1);
            obj.M{1}{1, 1} = spdiags(repmat(epsilon, ng, 1), 0, ng, ng);
            for i = 2:nh + 1
                obj.M{1}{i, i} = spdiags(repmat(epsilon*tau, ng, 1), 0, ng, ng);
            end
            obj.M{1} = core.linalg.block(obj.M{1});

            obj.M{2} = cell(nh+1, nh+1);
            obj.M{2}{1, 1} = speye(ng, ng);
            for i = 2:nh + 1
                obj.M{2}{i, i} = spdiags(repmat(epsilon, ng, 1), 0, ng, ng);
            end
            obj.M{2} = core.linalg.block(obj.M{2});
        end

        function obj = setLinearTerm(obj, state)
            nh = state.nModes - 1;
            ng = state.xDisc.nGlobalDofs;

            obj.L{1} = cell(nh+1, nh+1);
            obj.L{1} = repmat({sparse(ng, ng)}, nh+1, nh+1);
            obj.P{1} = repmat({sparse(ng, ng+1)}, nh+1, 1);
            obj.A{1} = repmat({sparse(ng+1, ng)}, 1, nh+1);
            obj.B{1} = sparse(ng+1, ng+1);
            obj.assembleDualOperators(state);
            obj.addLinearizedFreeTransport(state);
            obj.addCollision(state);
            obj.L{1} = core.linalg.block(obj.L{1});
            obj.P{1} = core.linalg.block(obj.P{1});
            obj.A{1} = core.linalg.block(obj.A{1});

            obj.L{2} = @(U) obj.computeL2(state, U);
            obj.P{2} = [];
            %             obj.A{2} = repmat({sparse(ng+1, ng)}, 1, nh+1);
            %             obj.A{2} = core.linalg.block(obj.A{2});
            %             obj.B{2} = speye(ng+1, ng+1);
            obj.A{2} = obj.A{1};
            obj.B{2} = obj.B{1};
        end

        function obj = assembleDualOperators(obj, ~)
            T0 = obj.config.T0;
            D = obj.adjoint.assembleMatrix();
            %             C1 = obj.multiplier.assembleConstant(state.coefs.EInf);
            C1 = obj.multiplier.assemblePointwise(obj.config.EInf);
            obj.dual = cell(2, 2);
            obj.dual{1, 2} = -sqrt(T0) * D{1, 2}{1} - C1 / (2 * sqrt(T0));
            obj.dual{2, 1} = sqrt(T0) * D{2, 1}{1} - C1 / (2 * sqrt(T0));
            obj.dual{1, 1} = sqrt(T0) * D{1, 1}{1};
            %             C2 = obj.multiplier.assembleConstant(1./state.coefs.rhoInf);
            %             obj.laplace = obj.dual{1,2} * C2 * obj.dual{2,1} - C2 * obj.dual{1,1};
            DW = obj.adjoint.assemblePointwiseMatrix(obj.config.invRhoInf);
            obj.dualw{1, 2} = -sqrt(T0) * DW{1, 2}{1} - C1 / (2 * sqrt(T0));
            obj.dualw{1, 1} = sqrt(T0) * DW{1, 1}{1};
            obj.laplace = obj.dualw{1, 2} * obj.dual{2, 1} - obj.dualw{1, 1};
        end

        function obj = addLinearizedFreeTransport(obj, state)
            nh = state.nModes - 1;
            ng = state.xDisc.nGlobalDofs;
            tau = obj.config.tau;

            % kinetic
            obj.L{1}{1, 2} = obj.dual{1, 2};
            for k = 2:nh
                obj.L{1}{k, k - 1} = -sqrt(k-1) * obj.dual{2, 1} * tau;
                obj.L{1}{k, k + 1} = sqrt(k) * obj.dual{1, 2} * tau;
            end
            obj.L{1}{nh + 1, nh} = -sqrt(nh) * obj.dual{2, 1} * tau;
            obj.P{1}{2} = [-obj.dual{2, 1} * tau, sparse(ng, 1)];

            % poisson
            obj.A{1}{1} = [speye(ng, ng); sparse(1, ng)];
            c = obj.source.assemblePointwiseVector(obj.config.invSqrtRhoInf);
            obj.B{1} = -[obj.laplace, c; c.', 0];
        end

        function obj = addCollision(obj, state)
            nh = state.nModes - 1;
            for k = 2:nh + 1
                T = obj.multiplier.assembleConstant(1);
                obj.L{1}{k, k} = -(k - 1) * T;
            end
        end

        function L = computeL2(obj, state, U)
            nh = state.nModes - 1;
            ng = state.xDisc.nGlobalDofs;
            omega = reshape(U(1:end-1), ng, []);
            omega = omega(:, nh+2);

            EInf = obj.config.EInf;
            sqrtRhoInf = obj.config.sqrtRhoInf;
            T0 = obj.config.T0;
            dOmega = @(x) reshape(state.xDisc.differentiate([], [], omega, 1), 1, []);
            f = @(x) (sqrt(T0) .* dOmega(x) - EInf(x) ./ (2 * sqrt(T0))) ./ sqrtRhoInf(x);
            T = obj.multiplier.assemblePointwise(f);

            L = repmat({sparse(ng, ng)}, nh+1, nh+1);
            for k = 2:nh + 1
                L{k, k-1} = -sqrt(k-1) * T;
            end
            L = core.linalg.block(L);
        end

        function obj = setSourceTerm(obj, state)
            nh = state.nModes - 1;
            ng = state.xDisc.nGlobalDofs;
            sqrtRhoInf = obj.source.assemblePointwiseVector(obj.config.sqrtRhoInf);
            obj.C = cell(2, 1);
            obj.C{1} = [-sqrtRhoInf(:); 0];
            obj.C{2} = obj.C{1};
            %             obj.C{2} = @(U) -U(end-ng:end);
            obj.S{1} = sparse(ng*(nh + 1), 1);
            obj.S{2} = @(U) obj.computeS2(state, U);
        end

        function S = computeS2(obj, state, U)
            nh = state.nModes - 1;
            ng = state.xDisc.nGlobalDofs;
            omega = reshape(U(1:end-1), ng, []);
            omega = omega(:, nh+2);
            T0 = obj.config.T0;
            EInf = obj.config.EInf;

            dOmega = @(x) reshape(state.xDisc.differentiate([], [], omega, 1), 1, []);
            f = @(x) sqrt(T0) .* dOmega(x) - EInf(x) ./ (2 * sqrt(T0));
            b = obj.source.assemblePointwiseVector(f);

            S = repmat({sparse(ng, 1)}, nh+1, 1);
            S{2} = b;
            S = vertcat(S{:});
        end
    end
end