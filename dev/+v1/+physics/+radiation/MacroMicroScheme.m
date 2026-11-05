classdef MacroMicroScheme < physics.Scheme
    properties
        pattern
        assembler
        operator
    end

    methods
        function obj = MacroMicroScheme(config, tDisc, visualizer, analyzer, pattern)
            obj@physics.Scheme(config, tDisc, visualizer, analyzer);
            obj.pattern = pattern;
            obj.assembler = struct();
            obj.operator = struct();
        end

        function state = initialize(obj, state)
            % INITIALIZE Initialize the macro-micro DG scheme.

            %< Reset timer
            obj.timer.reset();

            %< Set local data
            state.disc.v.fe.setVolumeData(0);
            state.disc.x.fe.setVolumeData(1);
            state.disc.x.fe.setFluxData(0);
            state.disc.x.op.setVolumeData();
            state.disc.x.op.setFluxData();

            %< Set initial condition
            state.dofs.U = state.macroProject(obj.config.ic.f);
            state.dofs.G = state.microProject(obj.config.ic.f);

            %< Set scattering coefficients
            if ~isempty(obj.config.scattering)
                state.coefs.CS = state.mean(obj.config.scattering);
            end

            %< Set absorption coefficients
            if ~isempty(obj.config.absorption)
                state.coefs.CA = state.mean(obj.config.absorption);
            end

            %< Set source coefficients
            if ~isempty(obj.config.source) && nargin(obj.config.source) < 3
                state.coefs.QU = state.macroProject(obj.config.source);
                state.coefs.QG = state.microProject(obj.config.source);
            end

            %< Set free streaming pattern
            obj.pattern.setFreeStreamingPattern(state);

            %< Set scaling pattern
            obj.pattern.setScalingPattern(state);

            %< Set upwind advection assembly
            obj.assembler.advection = fem.assembly.UpwindGridAssembly( ...
                state.disc.x.fe, state.disc.x.mesh, state.disc.x.op, ...
                ~isempty(obj.config.bc.f));

            %< Set diffusion assembly
            obj.assembler.diffusion = fem.assembly.DiffusionGridAssembly( ...
                state.disc.x.fe, state.disc.x.mesh, state.disc.x.op, ...
                ~isempty(obj.config.bc.f), obj.config.xDiffusionFlux(1), ...
                obj.config.xDiffusionFlux(2), obj.config.xDiffusionBoundaryJump);

            %< Set constant assembly
            obj.assembler.constant = fem.assembly.ConstantAssembly( ...
                state.disc.x.fe, state.disc.x.mesh, state.disc.x.op);

            %< Set mass operator
            obj.operator.M = obj.computeMassOperator(state);

            %< Set linear operator
            obj.operator.L = obj.computeLinearOperator(state);

            %< Reset time discretization
            obj.tDisc.reset();

            %< Reset visualizer
            obj.visualizer.reset();
            obj.visualizer.addDataset('DG', ...
                {'Color', 'b', ...
                'Marker', 'o', ...
                'LineStyle', 'none'});
            if obj.visualizer.hasExact
                obj.visualizer.addExact('U', ...
                    @(x, t) state.macroEvaluate(obj.config.exacts.f, x, [], t));
                obj.visualizer.addExact('G', ...
                    @(x, t) state.microEvaluate(obj.config.exacts.f, x, [], t));
                obj.visualizer.addDataset('REF', ...
                    {'Color', 'r', ...
                    'Marker', 'none', ...
                    'LineStyle', '-'});
            end

            if obj.visualizer.isEnabled
                obj.timer.start(2, '[M] Visualization.');
                obj.visualizer.plot(state.disc.x.space, state.dofs, obj.tDisc);
                obj.timer.stop(2);
            end
        end

        function state = preStep(obj, state)
            % PRESTEP Preprocess operations before time step.

            %< Set time step
            h = state.disc.x.mesh.measure;
            obj.tDisc.timeline.setTimeStep(h, obj.config.cfl);
            if isa(obj.tDisc, 'approx.odeint.BdfIntegrator')
                obj.tDisc.setCoefficients();
            end

            %< Update history
            W = [state.dofs.U(:); state.dofs.G(:)];
            obj.tDisc.update(W);
        end

        function state = step(obj, state)
            % STEP Perform one upwind DG time step.

            m = state.nMacroModes;
            n = state.nMicroModes;

            %< Get mass operator
            M = obj.operator.M;

            %< Get linear operator
            L = -obj.operator.L;

            %< Set time-dependent source term
            if ~isempty(obj.config.source)
                S = @(t) obj.computeSourceOperator(state, t);
            else
                S = [];
            end
            
            %< Integration step
            W = obj.tDisc.step(L, S, M);
            W = reshape(W, [], m+n);
            state.dofs.U = W(:, 1:m);
            state.dofs.G = W(:, m+1:m+n);
        end

        function state = postStep(obj, state)
            % POSTSTEP Postprocess operations after time step.

            if obj.visualizer.isEnabled
                obj.timer.start(2, '[M] Visualization.');
                obj.visualizer.plot(state.disc.x.space, state.dofs, obj.tDisc);
                obj.timer.stop(2);
            end

            obj.tDisc.timeline.advance();
        end

        function state = finalize(obj, varargin)
            % FINALIZE Finalize simulation with error analysis.

            if obj.analyzer.isEnabled
                obj.timer.start(3, '[M] Error Evaluation.');
                if length(varargin) == 1
                    state = varargin{1};
                    obj.analyzer.addExact('U', ...
                        @(x, t) state.macroEvaluate(obj.config.exacts.f, x, [], t));
                    obj.analyzer.addExact('G', ...
                        @(x, t) state.microEvaluate(obj.config.exacts.f, x, [], t));
                    obj.analyzer.setLevel(state.disc.x.mesh);
                    obj.analyzer.absolute(state.disc.x.space, state.dofs, ...
                        obj.tDisc.timeline.now);
                else
                    [state, coarse] = varargin{:};
                    obj.analyzer.setLevel(state.disc.x.mesh);
                    obj.analyzer.relative(state.disc.x.space, state.dofs, ...
                        coarse.disc.x.space, coarse.dofs);
                end
                obj.timer.stop(3);
            end

            obj.timer.report();
        end
    end

    methods (Access = protected)
        function M = computeMassOperator(obj, state)
            m = state.nMacroModes;
            n = state.nMicroModes;
            p = state.disc.x.space.nDofs;
            S = obj.pattern.scaling;
            M = cell(m + n, m + n);
            for i = 1:(m + n)
                M{i, i} = speye(p, p) * S(i, 1);
            end
            M = core.linalg.block(M);
        end

        function L = computeLinearOperator(obj, state)
            % Assemble linear operator sparse matrix

            m = state.nMacroModes;
            n = state.nMicroModes;

            L = cell(m+n, m+n);

            % Set streaming operator
            L = obj.computeFreeStreamingOperator(state, L);

            % Set collision operator
            L = obj.computeCollisionOperator(state, L);

            L = core.linalg.block(L);
        end

        function L = computeFreeStreamingOperator(obj, state, L)
            % Assemble streaming part of linear operator

            m = state.nMacroModes;
            n = state.nMicroModes;
            p = state.disc.x.space.nDofs;

            asmAux = obj.assembler.diffusion.aux.asm;
            asmPrm = obj.assembler.diffusion.prm.asm;
            asmUw = obj.assembler.advection;

            V = state.microNodes;

            % Macro-macro streaming
            for i = 1:m
                for j = 1:m
                    a = squeeze(obj.pattern.streaming(i, j, :));
                    if all(a == 0)
                        L{i, j} = sparse(p, p);
                        continue; 
                    end
                    if i < j
                        D = asmAux.divergence(a);
                    else
                        D = asmPrm.divergence(a);
                    end
                    L{i, j} = D * obj.pattern.scaling(i, 1+j);
                end
            end

            % Macro-micro streaming
            for i = 1:m
                for j = 1:n
                    a = squeeze(obj.pattern.streaming(i, m+j, :));
                    if all(a == 0)
                        L{i, m + j} = sparse(p, p);
                        continue; 
                    end
                    D = asmAux.divergence(a);
                    L{i, m+j} = D * obj.pattern.scaling(i, 1+m+j);
                end
            end

            % Micro-macro streaming
            for i = 1:n
                for j = 1:m
                    a = squeeze(obj.pattern.streaming(m + i, j, :));
                    if all(a == 0)
                        L{m + i, j} = sparse(p, p);
                        continue; 
                    end
                    D = asmPrm.divergence(a);
                    L{m + i, j} = D * obj.pattern.scaling(m+i, 1+j);
                end
            end

            % Micro-micro streaming
            for i = 1:n
                for j = 1:n
                    a = squeeze(obj.pattern.streaming(m + i, m + j, :));
                    if all(a == 0)
                        L{m + i, m + j} = sparse(p, p);
                        continue; 
                    end
                    D = asmUw.divergence(a, V(:, j));
                    L{m + i, m + j} = D * obj.pattern.scaling(m+i, 1+m+j);
                end
            end
        end

        function L = computeCollisionOperator(obj, state, L)
            m = state.nMacroModes;
            n = state.nMicroModes;
            asm = obj.assembler.constant;

            % Scattering
            if ~isempty(state.coefs.CS)
                C = asm.assemble(state.coefs.CS);
                for i = 1:(m + n)
                    L{i, i} = L{i, i} + C * obj.pattern.scaling(i, m+n+2);
                end
            end

            % Absorption
            if ~isempty(state.coefs.CA)
                C = asm.assemble(state.coefs.CA);
                for i = 1:(m + n)
                    L{i, i} = L{i, i} + C * obj.pattern.scaling(i, m+n+3);
                end
            end
        end
    
        function S = computeSourceOperator(obj, state, t)
            state.coefs.QU = state.macroProject(obj.config.source, [], [], t);
            state.coefs.QG = state.microProject(obj.config.source, [], [], t);
            m = state.nMacroModes;
            n = state.nMicroModes;
            p = state.disc.x.fe.nDofs;
            C = reshape(repmat(state.coefs.CA(:).', p, 1), [], 1);
            QU = C .* state.coefs.QU .* obj.pattern.scaling(1:m, m+n+3).';
            QG = C .* state.coefs.QG .* obj.pattern.scaling(m+1:m+n, m+n+3).';
            S = [QU(:); QG(:)];
        end
    
    end
end