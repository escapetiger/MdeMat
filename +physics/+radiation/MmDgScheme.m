classdef MmDgScheme < physics.radiation.RadiationScheme
    % MMDGSCHEME MMDG scheme for linear kinetic equations.
    %
    %   MmDgScheme implements the macro-micro discontinuous Galerkin scheme
    %   for multiscale linear kinetic equations. This scheme uses
    %   macro-micro decomposition f = u + g where u is the macroscopic
    %   component and g is the microscopic component.
    %
    % See also:
    %   physics.radiation.RadiationScheme, physics.radiation.MacroMicroState

    properties
        Pattern % Kinetic pattern

        AdjAmb % Adjoint assembly
        MultAmb % Multiplier assembly
        SrcAmb % Source assembly
        UwAmb % Upwind assembly

        Dual % Dual operators
    end

    properties (Constant)
        Name = 'MMDG' % Scheme name
    end

    methods
        function state = initialize(obj, state)
            % INITIALIZE Initialize the macro-micro scheme.
            %
            %   state = initialize(obj, state) sets up the macro-micro
            %   scheme by setting up decomposition patterns, assembling
            %   operators, and preparing initial conditions.

            arguments
                obj physics.radiation.MmDgScheme
                state physics.radiation.MacroMicroState
            end

            %< Set up linear solver options
            obj.TDisc.reset();

            %< Set up decomposition pattern
            state.setDecomposition(obj.Config.decomposition);
            state.setScaling(obj.Config.epsilon);

            %< Set up free transport pattern
            obj.setFreeTransport(state);

            %< Set up scaling pattern
            obj.setScaling(state);

            %< Set initial condition
            state.setDof('U', state.macroFit(obj.Config.ic));
            state.setDof('G', state.microFit(obj.Config.ic));
            state.kineticReconstruct();

            %< Set scattering coefficients
            if ~isempty(obj.Config.scattering)
                state.setCoefficient('CS', state.XDisc.fit(obj.Config.scattering));
            end

            %< Set absorption coefficients
            if ~isempty(obj.Config.absorption)
                state.setCoefficient('CA', state.XDisc.fit(obj.Config.absorption));
            end

            %< Set source coefficients
            if ~isempty(obj.Config.source) && nargin(obj.Config.source) < 3
                state.setCoefficient('SU', ...
                    @(x, varargin) state.macroFit(obj.Config.source, 'x', x, varargin{:}));
                state.setCoefficient('SG', ...
                    @(x, varargin) state.microFit(obj.Config.source, 'x', x, varargin{:}));
            end

            %< Assembly for adjoint operators
            obj.AdjAmb = approx.assembly.AdjointAssembly( ...
                state.XDisc, obj.Config.bcType, obj.Config.xPenaltyType);

            %< Assembly for multipliers
            obj.MultAmb = approx.assembly.MultiplierAssembly( ...
                state.XDisc);

            %< Assembly for source
            obj.SrcAmb = approx.assembly.SourceAssembly( ...
                state.XDisc);

            %< Assembly for upwind fluxes
            obj.UwAmb = approx.assembly.UpwindAssembly( ...
                state.XDisc, obj.Config.bcType);

            %< Assemble precomputed operators

            %< Set mass term
            obj.setMassTerm(state);

            %< Set linear term
            obj.setLinearTerm(state);

            %< Set source term
            obj.setSourceTerm(state);

            %< Initialize visualizer with exact solution if available
            visualizer = obj.getConfig('visualizer');
            if ~isempty(visualizer) && visualizer.IsEnabled
                fnU = @(x, varargin) state.macroFEval(obj.Config.exact, 'x', x, varargin{:});
                fnG = @(x, varargin) state.microFEval(obj.Config.exact, 'x', x, varargin{:});
                visualizer.addExact('U', fnU);
                visualizer.addExact('G', fnG);
            end

            %< Initialize analyzer with exact solution if available
            analyzer = obj.getConfig('analyzer');
            if ~isempty(obj.Config.exact) && analyzer.IsEnabled
                fnU = @(x, varargin) state.macroFEval(obj.Config.exact, 'x', x, varargin{:});
                fnG = @(x, varargin) state.microFEval(obj.Config.exact, 'x', x, varargin{:});
                fnF = @(x, varargin) state.kineticFEval(obj.Config.exact, 'x', x, varargin{:});
                analyzer.addExact('U', fnU);
                analyzer.addExact('G', fnG);
                analyzer.addExact('F', fnF);
            end

            %< Set initial history
            state.updateHistory(obj.TDisc.Timeline.Now);
        end

        function state = preStep(obj, state)
            % PRESTEP Preprocess operations before time step.
            %
            %   state = preStep(obj, state) performs preprocessing
            %   including time step size computation based on CFL
            %   condition.

            arguments
                obj physics.radiation.MmDgScheme
                state physics.radiation.MacroMicroState
            end

            if obj.hasConfig('cfl') && ~isempty(obj.Config.cfl)
                h = state.XDisc.Mesh.computeMeasure();
                obj.TDisc.setTimeStep(h=h, C=obj.Config.cfl, p=1);
            end

            if obj.hasConfig('dt') && ~isempty(obj.Config.dt)
                obj.TDisc.setTimeStep(dt=obj.Config.dt);
            end
        end

        function state = step(obj, state)
            % STEP Perform one macro-micro time step.
            %
            %   state = step(obj, state) advances the solution by one time
            %   step using the orthogonal implicit integrator method.

            arguments
                obj physics.radiation.MmDgScheme
                state physics.radiation.MacroMicroState
            end

            m = state.NVMacroDofs;
            n = state.NVMicroDofs;

            if m > 0 && n == 0
                W0 = state.Dofs.U(:);
            elseif m == 0 && n > 0
                W0 = state.Dofs.G(:);
            else
                W0 = [state.Dofs.U(:); state.Dofs.G(:)];
            end
            obj.TDisc.update(W0);

            W = obj.TDisc.step(L=obj.L, S=obj.S, M=obj.M);

            if m > 0 && n == 0
                state.Dofs.U = reshape(W, [], m);
            elseif m == 0 && n > 0
                state.Dofs.G = reshape(W, [], n);
            else
                W = reshape(W, [], m+n);
                state.Dofs.U = W(:, 1:m);
                state.Dofs.G = W(:, m+1:m+n);
            end

            %< Apply filter if enabled
            if obj.Config.useFilter && m > 0
                state = obj.applyFilter(state);
            end

            %< Apply positivity limiter if enabled
            if obj.Config.usePositivityLimiter && m > 0
                state = obj.applyPositivityLimiter(state);
            end

            %< Reconstruct distribution
            state.kineticReconstruct();
        end

        function state = postStep(obj, state)
            % POSTSTEP Postprocess operations after time step.
            %
            %   state = postStep(obj, state) performs postprocessing
            %   including timeline advancement.

            arguments
                obj physics.radiation.MmDgScheme
                state physics.radiation.MacroMicroState
            end

            obj.TDisc.advance();
        end

        function state = finalize(obj, state)
            % FINALIZE Finalize simulation with scheme and state reporting.
            %
            %   state = finalize(obj, state) performs final reporting
            %   including scheme configuration and state information.

            arguments
                obj physics.radiation.MmDgScheme
                state physics.radiation.MacroMicroState
            end

            state.report();

            fprintf('[R] ODE Integrator: %s\n', class(obj.TDisc));
            fprintf('[R] Final Time: %.6g\n', obj.Config.tFinal);
            fprintf('[R] Epsilon: %.6g\n', obj.Config.epsilon);
            fprintf('[R] Macro Modes: %d\n', state.NVMacroDofs);
            fprintf('[R] Micro Modes: %d\n', state.NVMicroDofs);
            if obj.hasConfig('cfl') && ~isempty(obj.Config.cfl)
                fprintf('[R] CFL number: %.3f\n', obj.Config.cfl);
            end
            if obj.hasConfig('dt') && ~isempty(obj.Config.dt)
                fprintf('[R] Time step size: %.3f\n', obj.Config.dt);
            end
            fprintf('[R] Boundary Conditions: %s\n', obj.Config.bcType);
        end
    end

    methods (Access = protected)
        function obj = setup(obj)
            % SETUP Setup default configuration options for macro-micro
            % schemes.

            %< Call parent setup
            obj = setup@physics.radiation.RadiationScheme(obj);

            %< Add MMDG-specific required configuration options
            obj.addConfig('xBasisOrder', default=1, ...
                validator=@(x) isnumeric(x) && isscalar(x) && (x > 0));
            obj.addConfig('xBasisType', default='nodal', ...
                validator=@(x) (ischar(x) || isstring(x)) && ...
                    ismember(lower(x), {'modal', 'nodal'}));
            obj.addConfig('xBasisPattern', default='Q', ...
                validator=@(x) (ischar(x) || isstring(x)) && ...
                    ismember(upper(x), {'P', 'Q'}));
            obj.addConfig('xPenaltyType', default={'', 'right'; 'left', ''}, ...
                validator=@(x) iscell(x) && isequal(size(x), [2, 2]));
            obj.addConfig('decomposition', default='', ...
                validator=@(x) (ischar(x) || isstring(x)) && ...
                    ismember(lower(x), {'', 'partial', 'full'}));
            obj.addConfig('timeScale', default=0, ...
                validator=@(x) isnumeric(x) && isscalar(x));
            obj.addConfig('scatteringScale', default=0, ...
                validator=@(x) isnumeric(x) && isscalar(x));
            obj.addConfig('absorptionScale', default=0, ...
                validator=@(x) isnumeric(x) && isscalar(x));
            obj.addConfig('nu', default=0, ...
                validator=@(x) isnumeric(x) && isscalar(x) && (x >= 0));
            obj.addConfig('nv', default=0, ...
                validator=@(x) isnumeric(x));
            obj.addConfig('M', default=0, ...
                validator=@(x) isnumeric(x) && isscalar(x) && (x >= 0));
            obj.addConfig('N', default=0, ...
                validator=@(x) isnumeric(x) && isscalar(x) && (x >= 0));

            %< Add positivity limiter configuration options
            obj.addConfig('usePositivityLimiter', default=false, ...
                validator=@(x) islogical(x) || isnumeric(x));
            obj.addConfig('positivityLimiterType', default='zhang_shu', ...
                validator=@(x) (ischar(x) || isstring(x)) && ...
                    ismember(lower(x), {'cutoff', 'zhang_shu'}));

            %< Add filter configuration options
            obj.addConfig('useFilter', default=false, ...
                validator=@(x) islogical(x) || isnumeric(x));
            obj.addConfig('filterType', default='rot', ...
                validator=@(x) (ischar(x) || isstring(x)) && ...
                    ismember(lower(x), {'rot', 'exp', ''}));
            obj.addConfig('filterStrength', default=0.1, ...
                validator=@(x) isnumeric(x) && isscalar(x) && (x > 0) && (x <= 1));
            obj.addConfig('filterOrder', default=4, ...
                validator=@(x) isnumeric(x) && isscalar(x) && (x > 0));
            obj.addConfig('filterCutoff', default=0, ...
                validator=@(x) isnumeric(x) && isscalar(x) && (x >= 0) && (x <= 1));
            obj.addConfig('filterCollisionRate', default=0, ...
                validator=@(x) isnumeric(x) && isscalar(x) && (x >= 0));
        end

        function obj = configure(obj)
            % CONFIGURE Generate macro-micro specific configuration dependencies.

            %< Call parent configure
            obj = configure@physics.radiation.RadiationScheme(obj);

            %< Determine velocity discretization parameters
            nDims = obj.Config.nDims;
            vDimReduction = obj.Config.vDimReduction;
            if obj.Config.nu > 0
                if nDims == 1 && strcmpi(vDimReduction, 'topology')
                    obj.setConfig('M', obj.Config.nu);
                elseif nDims == 1 && strcmpi(vDimReduction, 'symmetry')
                    obj.setConfig('M', obj.Config.nu);
                elseif nDims == 2 && strcmpi(vDimReduction, 'topology')
                    obj.setConfig('M', 2 * obj.Config.nu - 1);
                else
                    obj.setConfig('M', obj.Config.nu^2);
                end
            else
                obj.setConfig('M', 0);
            end
            obj.setConfig('N', prod(obj.Config.nv));

            %< Determine spatial discretization identifier
            xId = sprintf('%s%d', upper(obj.Name), obj.Config.xBasisOrder);
            vId = sprintf('P%dS%d', obj.Config.M, obj.Config.N);
            obj.setConfig('xId', xId);
            obj.setConfig('vId', vId);
            obj.setConfig('sId', sprintf('%s%s', xId, vId));

            %< Set plot title
            eId = obj.Config.eId;
            obj.setConfig('titlePrefix', sprintf('%s with %s', eId, obj.Config.sId));

            %< Set visualization options
            visualizer = obj.Config.visualizer;
            visualizer.setTitlePrefix(obj.Config.titlePrefix);
        end
    
    end

    methods (Access = private)
        function obj = setFreeTransport(obj, state)
            % SETFREETRANSPORT Set up free transport operators.

            tol = 1e-13;
            m = state.NVMacroDofs;
            n = state.NVMicroDofs;
            nxd = state.NXDims;

            if m > 0 && n == 0
                EU = state.VLhsBasisValues.';
                VU = state.VLhsNodes;
                WU = diag(state.VLhsWeights);
                T = zeros(m, m, nxd);
                for d = 1:nxd
                    VUd = diag(VU(d, :));
                    T(:, :, d) = EU.' * WU * VUd * EU;
                end
            elseif m == 0 && n > 0
                VG = state.VRhsNodes;
                T = zeros(n, n, nxd);
                for d = 1:nxd
                    VGd = diag(VG(d, :));
                    T(:, :, d) = VGd;
                end
            else
                vG = state.VDisc.Rhs.Element.Volume.Nodes;
                BL = state.VDisc.Lhs.Element.Approximator.Basis;
                EU = state.VLhsBasisValues.';
                EG = BL.eval(vG).';
                VU = state.VLhsNodes;
                VG = state.VRhsNodes;
                WU = diag(state.VLhsWeights);
                M = BL.Rhs.eval(vG).^2;
                WG = diag(state.VRhsWeights(:) ./ M(:));
                L = state.Levels(1:m);
                S = ones(1, m);
                S(L+1 < max(L)+1) = 0;
                S = diag(S);
                T = zeros(m+n, m+n, nxd);
                for d = 1:nxd
                    VUd = diag(VU(d, :));
                    VGd = diag(VG(d, :));
                    TUUd = EU.' * WU * VUd * EU;
                    TUGd = EG.' * WG * VGd;
                    T(1:m, 1:m, d) = TUUd;
                    T(1:m, m+1:m+n, d) = S * TUGd;
                    T(m+1:m+n, 1:m, d) = (VGd * EG - EG * TUUd);
                    T(m+1:m+n, m+1:m+n, d) = (VGd - EG * S * TUGd);
                end
            end
            T(abs(T) < tol) = 0;
            obj.Pattern.freeTransport = T;
        end

        function obj = setScaling(obj, state)
            % SETSCALING Set up scaling weights based on epsilon.

            m = state.NVMacroDofs;
            n = state.NVMicroDofs;
            epsilon = obj.Config.epsilon;

            sct = obj.Config.timeScale;
            scs = obj.Config.scatteringScale;
            sca = obj.Config.absorptionScale;

            if m > 0
                scm = state.Decomposition(1:m);
            else
                scm = [];
            end

            if n > 0
                scn = state.Decomposition(m+1);
            else
                scn = [];
            end

            if m > 0 && n == 0
                sc = scm;
            elseif m == 0 && n > 0
                sc = repmat(scn, 1, n);
            else
                sc = [scm, repmat(scn, 1, n)];
            end
            S1 = sct + sc(:);
            S2 = repmat(sc, length(sc), 1);
            S2(all(obj.Pattern.freeTransport == 0, 3)) = inf;
            S3 = inf(length(sc), 1);
            if ~isempty(scs)
                if m > 0
                    S3(2:end) = scs + sc(2:end);
                else
                    S3 = scs + sc(:);
                end
            end
            S4 = inf(length(sc), 1);
            if ~isempty(sca), S4 = sca + sc(:); end
            S = [S1, S2, S3, S4];
            S = S - min(S, [], 2);
            obj.Pattern.scaling = (epsilon.^S) .* (S ~= inf);
        end

        function obj = setMassTerm(obj, state)
            % SETMASSTERM Assemble mass matrix for macro-micro scheme.

            m = state.NVMacroDofs;
            n = state.NVMicroDofs;
            ng = state.XDisc.NGlobalDofs;
            PS = obj.Pattern.scaling;
            obj.M = cell(m + n, m + n);
            for i = 1:(m + n)
                obj.M{i, i} = speye(ng, ng) * PS(i, 1);
            end
            obj.M = core.linalg.block(obj.M);
        end

        function obj = setLinearTerm(obj, state)
            % SETLINEARTERM Assemble linear operator sparse matrix.

            nxd = state.NXDims;
            obj.Dual = obj.AdjAmb.assembleMatrix(eye(nxd));

            m = state.NVMacroDofs;
            n = state.NVMicroDofs;

            obj.L = cell(m+n, m+n);

            obj.addFreeTransport(state);

            obj.addCollision(state);

            obj.L = core.linalg.block(obj.L);
        end

        function obj = addFreeTransport(obj, state)
            % ADDFREETRANSPORT Assemble free transport part of linear
            % operator.

            nxd = state.NXDims;
            m = state.NVMacroDofs;
            n = state.NVMicroDofs;
            ng = state.XDisc.NGlobalDofs;
            PFT = obj.Pattern.freeTransport;
            PS = obj.Pattern.scaling;

            %< Macro-macro transport
            for i = 1:m
                for j = 1:m
                    a = squeeze(PFT(i, j, :));
                    if all(a == 0)
                        obj.L{i, j} = sparse(ng, ng);
                        continue;
                    end
                    if i < j
                        D = sparse(ng, ng);
                        for d = 1:nxd
                            if a(d) == 0, continue; end
                            D = D + a(d) * obj.Dual{1, 2}{d};
                        end
                        if ~isempty(obj.Config.bc) && ~isempty(obj.Config.xPenaltyType{1, 1})
                            %< Jump term
                            Dj = sparse(ng, ng);
                            for d = 1:nxd
                                if a(d) == 0, continue; end
                                Dj = Dj - a(d) * obj.Dual{1, 1}{d};
                            end
                            obj.L{i, i} = obj.L{i, i} - Dj * PS(i, 1+j);
                        end
                    else
                        D = sparse(ng, ng);
                        for d = 1:nxd
                            if a(d) == 0, continue; end
                            D = D + a(d) * obj.Dual{2, 1}{d};
                        end
                    end
                    obj.L{i, j} = - D * PS(i, 1+j);
                end
            end

            %< Macro-micro transport
            for i = 1:m
                for j = 1:n
                    a = squeeze(PFT(i, m+j, :));
                    if all(a == 0)
                        obj.L{i, m+j} = sparse(ng, ng);
                        continue;
                    end
                    D = sparse(ng, ng);
                    for d = 1:nxd
                        if a(d) == 0, continue; end
                        D = D + a(d) * obj.Dual{1, 2}{d};
                    end
                    obj.L{i, m+j} = -D * PS(i, 1+m+j);

                    if ~isempty(obj.Config.bc) && ~isempty(obj.Config.xPenaltyType{1, 1})
                        %< Jump term
                        Dj = sparse(ng, ng);
                        if mod(state.Levels(i), 2) == 0
                            aa = squeeze(PFT(m+j, i, :));
                            for d = 1:nxd
                                if a(d) == 0 || aa(d) == 0, continue; end
                                Dj = Dj - aa(d) * a(d) * obj.Dual{1, 1}{d};
                            end
                        else
                            for d = 1:nxd
                                if a(d) == 0, continue; end
                                Dj = Dj - a(d) * obj.Dual{1, 1}{d};
                            end
                        end
                        obj.L{i, i} = obj.L{i, i} - Dj * PS(i, 1+m+j);
                    end
                end
            end

            %< Micro-macro transport
            for i = 1:n
                for j = 1:m
                    a = squeeze(PFT(m + i, j, :));
                    if all(a == 0)
                        obj.L{m + i, j} = sparse(ng, ng);
                        continue;
                    end
                    D = sparse(ng, ng);
                    for d = 1:nxd
                        if a(d) == 0, continue; end
                        D = D + a(d) * obj.Dual{2, 1}{d};
                    end
                    obj.L{m + i, j} = -D * PS(m+i, 1+j);
                end
            end

            %< Micro-micro transport
            VG = state.VRhsNodes;
            for i = 1:n
                for j = 1:n
                    a = squeeze(PFT(m + i, m + j, :));
                    if all(a == 0)
                        obj.L{m + i, m + j} = sparse(ng, ng);
                        continue;
                    end
                    D = obj.UwAmb.assembleMatrix(a, VG(:, j));
                    obj.L{m + i, m + j} = - D * PS(m+i, 1+m+j);
                end
            end
        end

        function obj = addCollision(obj, state)
            % ADDCOLLISION Add collision terms to linear operator.

            m = state.NVMacroDofs;
            n = state.NVMicroDofs;
            PS = obj.Pattern.scaling;

            %< Scattering
            if ~isempty(obj.Config.scattering)
                C = obj.MultAmb.assembleMatrix(obj.Config.scattering);
                for i = 1:(m + n)
                    obj.L{i, i} = obj.L{i, i} - C * PS(i, m+n+2);
                end
                
                if m == 0 && n > 0
                    M = obj.Config.E;
                    w = state.VRhsWeights * M;
                    for i = 1:(m+n)
                        for j = 1:(m+n)
                            obj.L{i, j} = obj.L{i, j} + C * w(j) * PS(i, m+n+2);
                        end
                    end
                end
            end

            %< Absorption
            if ~isempty(obj.Config.absorption)
                C = obj.MultAmb.assembleMatrix(obj.Config.absorption);
                for i = 1:(m + n)
                    obj.L{i, i} = obj.L{i, i} - C * PS(i, m+n+3);
                end
            end
        end

        function obj = setSourceTerm(obj, state)
            % SETSOURCETERM Set up source term function.

            m = state.NVMacroDofs;
            n = state.NVMicroDofs;
            ng = state.XDisc.NGlobalDofs;

            if ~isempty(obj.Config.source) || ~isempty(obj.Config.bc)
                obj.S = @(W, t) computeS(W, t);
            else
                obj.S = [];
            end

            function S = computeS(~, t)
                S = repmat({sparse(ng, 1)}, m+n, 1);

                if ~isempty(obj.Config.source)
                    S = obj.addSource(state, t, S);
                end

                if ~isempty(obj.Config.bc)
                    S = obj.addBc(state, t, S);
                end

                S = vertcat(S{:});
            end
        end

        function S = addSource(obj, state, t, S)
            % ADDSOURCE Add source contributions to right-hand side.

            m = state.NVMacroDofs;
            n = state.NVMicroDofs;
            PS = obj.Pattern.scaling;
            % nl = state.XDisc.NLocalDofs;

            state.Coefs.SU = state.macroFit(obj.Config.source, t=t);
            state.Coefs.SG = state.microFit(obj.Config.source, t=t);

            SU = state.Coefs.SU .* PS(1:m, m+n+3).';
            SG = state.Coefs.SG .* PS(m+1:m+n, m+n+3).';

            for i = 1:m
                S{i} = SU(:, i);
            end
            for i = 1:n
                S{m+i} = SG(:, i);
            end
        end

        function S = addBc(obj, state, t, S)
            % ADDBC Add boundary condition contributions to right-hand side.

            PFT = obj.Pattern.freeTransport;
            PS = obj.Pattern.scaling;
            m = state.NVMacroDofs;
            n = state.NVMicroDofs;
            ng = state.XDisc.NGlobalDofs;
            nd = state.NXDims;
            
            %< Macro Bcs
            fU = @(x, LFI) state.macroTraceEvaluate(LFI, obj.Config.bc, x=x, t=t);
            if m > 0
                bU = obj.AdjAmb.assembleVector(eye(nd), fU, mode='face');
            end

            %< Macro-macro transport BCs
            for i = 1:m
                for j = 1:m
                    a = squeeze(PFT(i, j, :));
                    if all(a == 0), continue; end
                    if ~isempty(obj.Config.xPenaltyType{1, 1})
                        if i < j
                            D = sparse(ng, 1);
                            for d = 1:nd
                                if a(d) == 0, continue; end
                                D = D - a(d) * bU{1, 1}{d}(:, i);
                            end
                            S{i} = S{i} - D * PS(i, 1+j);
                        end
                    end

                    if i > j 
                        D = sparse(ng, 1);
                        for d = 1:nd
                            if a(d) == 0, continue; end
                            D = D + a(d) * bU{2, 1}{d}(:, j);
                        end
                        S{i} = S{i} - D * PS(i, 1+j);
                    end
                end
            end

            %< Macro-micro transport BCs
            for i = 1:m
                for j = 1:n
                    a = squeeze(PFT(i, m + j, :));
                    if all(a == 0), continue; end
                    if ~isempty(obj.Config.xPenaltyType{1, 1})
                        D = sparse(ng, 1);
                        if mod(state.Levels(i), 2) == 0
                            aa = squeeze(PFT(m+j, i, :));
                            for d = 1:nd
                                if a(d) == 0 || aa(d) == 0, continue; end
                                D = D - aa(d) * a(d) * bU{1, 1}{d}(:, i);
                            end
                        else
                            for d = 1:nd
                                if a(d) == 0, continue; end
                                D = D - a(d) * bU{1, 1}{d}(:, i);
                            end
                        end
                        S{i} = S{i} - D * PS(i, 1+m+j);
                    end
                end
            end

            %< Micro-macro transport BCs
            for i = 1:n
                for j = 1:m
                    a = squeeze(PFT(m + i, j, :));
                    if all(a == 0), continue; end
                    D = sparse(ng, 1);
                    for d = 1:nd
                        if a(d) == 0, continue; end
                        D = D + a(d) * bU{2, 1}{d}(:, j);
                    end
                    S{m+i} = S{m+i} - D * PS(m+i, 1+j);
                end
            end

            %< Micro-micro transport BCs
            % state.reset();
            if n > 0
                V = state.VRhsNodes;
                fG = @(x, LFI) state.microTraceEvaluate(LFI, obj.Config.bc, x=x, t=t);
                bG = obj.AdjAmb.assembleVector(eye(nd), fG, mode='face_full');
            end

            for i = 1:n
                for j = 1:n
                    a = squeeze(PFT(m + i, m + j, :));
                    if all(a == 0), continue; end
                    if ~isa(state.XDisc.Mesh, 'approx.mesh.Grid')
                        D = sparse(ng, 1);
                        for d = 1:nd
                            if a(d) == 0, continue; end
                            if V(d, j) > 0
                                D = D + a(d) * bG{2, 1}{2*d-1}{d}(:, j);
                            else
                                D = D + a(d) * bG{2, 1}{2*d}{d}(:, j);
                            end
                        end
                        S{m+i} = S{m+i} - D * PS(m+i, 1+m+j);
                    else
                        nl = state.XDisc.NLocalDofs;
                        nf = state.XDisc.Element.Geometry.NFaces;
                        for LFI = 1:nf
                            BEI = state.XDisc.Mesh.getBoundaryElements(LFI);
                            normal = state.XDisc.Mesh.computeOutwardNormals(BEI, LFI);
                            BEI = BEI(sum(normal .* V(:, j), 1) < zeros(1, length(BEI)));
                            if ~isempty(BEI)
                                K = bsxfun(@plus, (BEI(:).' - 1)*nl, (1:nl).');
                                K = K(:);
                                D = sparse(length(K), 1);
                                for d = 1:nd
                                    if a(d) == 0, continue; end
                                    G = bG{2, 1}{LFI}{d}(:, j);
                                    D = D + a(d) * G(K);
                                end
                                D = sparse(K, 1, D, ng, 1);
                                S{m+i} = S{m+i} - D * PS(m+i, 1+m+j);
                            end
                        end
                    end
                end
            end
        end

        function state = applyFilter(obj, state)
            % APPLYFILTER Apply spectral filter to macro angular coefficients.
            %
            %   state = applyFilter(obj, state) damps the higher angular
            %   modes of the macro DOFs @a state.Dofs.U to suppress Gibbs
            %   oscillations from the truncated angular basis.
            %
            %   Two filter types are available:
            %
            %   'rot' — Rotational hyperviscosity (Padé form). Adds
            %   biharmonic angular diffusion \f$ -\nu_r (-\partial_\theta^2)^2 u \f$.
            %   Intermediate modes are time-consistent (\f$ \sigma_k \to 1 \f$
            %   as \f$ \Delta t \to 0 \f$); the top mode has fixed attenuation
            %   @a filterStrength per step.
            %
            %   'exp' — Vandeven exponential filter. Applies
            %   \f$ \sigma_k = \exp(-\alpha \eta_k^p) \f$ where
            %   \f$ \eta_k = (k - k_c)/(N - k_c) \f$. Preserves accuracy of
            %   order @a filterOrder below the cutoff.
            %
            %   Both types interpret @a filterStrength as the target
            %   attenuation \f$ \sigma_N^* \in (0,1] \f$ at the top mode
            %   \f$ k = N \f$. Set @a filterCollisionRate to \f$ \sigma_s +
            %   \sigma_a \f$ to compensate for physical collision damping
            %   that already acts on higher modes.

            m = state.NVMacroDofs;
            if m == 0
                return;
            end

            k = state.Levels(1:m);
            kmax = max(k);

            if kmax == 0
                return;
            end

            tp = obj.Config.filterType;
            sN = obj.Config.filterStrength;
            p = obj.Config.filterOrder;
            dt = obj.TDisc.Timeline.StepSize;
            kc = round(obj.Config.filterCutoff * kmax);

            if kmax <= kc
                return;
            end

            %< Adjust target for physical collision damping already present
            cr = obj.Config.filterCollisionRate;
            sN_eff = sN * exp(cr * dt);
            if sN_eff >= 1
                return;
            end

            %< Normalized level relative to cutoff
            eta = max(k - kc, 0) / (kmax - kc);

            switch lower(tp)
                case 'rot'
                    %< Biharmonic rotational diffusion (Pade form)
                    %< sigma_k = 1/(1 + c * (k-kc)^4)
                    %< where c is chosen so that sigma(kmax) = sN_eff
                    alpha = (1/sN_eff - 1) / (kmax - kc)^4;
                    % alpha = alpha .* dt.^(1 - eta.^p);
                    sigma = 1 ./ (1 + alpha .* max(k - kc, 0).^4);
                case 'exp'
                    %< Vandeven exponential filter (time-consistent)
                    %< sigma_k = exp(-alpha * dt * eta^p)
                    %< filterStrength = target cumulative attenuation per unit time
                    %< at k=kmax: cumulative sigma = sN_eff^T (dt-independent)
                    alpha = -log(sN_eff);
                    sigma = exp(-alpha .* dt .* eta.^p);
                otherwise
                    sigma = ones(1, m);
            end

            sigma(k <= kc) = 1;
            state.Dofs.U = state.Dofs.U .* sigma;
        end

        function state = applyPositivityLimiter(obj, state)
            % APPLYPOSITIVITYLIMITER Apply positivity limiter to density.
            %
            %   state = applyPositivityLimiter(obj, state) enforces
            %   non-negativity of the zeroth velocity moment (density) in
            %   the macro DOFs @a state.Dofs.U(:,1).
            %
            %   Two limiter types are available via config.positivityLimiterType:
            %
            %   - 'cutoff': Clamps cell-average DOFs to zero. Simple but
            %     NOT mass-conservative; adds spurious mass when cell averages
            %     are negative.
            %
            %   - 'zhang_shu' (default): Zhang-Shu (2010) limiter. Scales
            %     within-element spatial fluctuations to enforce
            %     \f$ u_0(x) \geq 0 \f$ at all quadrature nodes while
            %     preserving the cell average exactly. Mass-conservative
            %     whenever the cell average is non-negative. If the cell
            %     average itself is negative, it is clamped to zero (mass
            %     loss is unavoidable in that case).

            if strcmpi(obj.Config.xBasisType, 'modal')
                nl = state.XDisc.NLocalDofs;
                ne = state.XDisc.NMeshElements;
                C = reshape(state.Dofs.U(:, 1), nl, ne);  %< (nl x ne)

                switch lower(obj.Config.positivityLimiterType)
                    case 'zhang_shu'
                        %< Cell averages: first DOF per element (monic Legendre)
                        u_avg = C(1, :);  %< (1 x ne)

                        %< Evaluate density at element quadrature nodes
                        xRef = state.XDisc.Element.Volume.Nodes;
                        u_vals = state.XDisc.Element.eval(xRef, C);  %< (ne x np)

                        %< Minimum density per element over quadrature nodes
                        u_min = min(u_vals, [], 2).';  %< (1 x ne)

                        %< Zhang-Shu theta: largest scaling that ensures u >= 0
                        theta = ones(1, ne);
                        mask = (u_min < 0) & (u_avg > 0);
                        theta(mask) = u_avg(mask) ./ (u_avg(mask) - u_min(mask));
                        theta = min(1, theta);

                        %< Scale within-element fluctuations; preserve averages
                        C(2:end, :) = C(2:end, :) .* theta;

                        %< Clamp negative cell averages (accept unavoidable mass loss)
                        C(1, :) = max(0, u_avg);

                    otherwise  %< 'cutoff'
                        C(1, :) = max(0, C(1, :));
                end

                state.Dofs.U(:, 1) = C(:);
            else
                %< Nodal basis: DOFs are already values at quadrature nodes
                nl = state.XDisc.NLocalDofs;
                ne = state.XDisc.NMeshElements;
                U = reshape(state.Dofs.U(:, 1), nl, ne);  %< (nl x ne)

                switch lower(obj.Config.positivityLimiterType)
                    case 'zhang_shu'
                        %< Cell averages via quadrature weights (unit ref element)
                        w = state.XDisc.Element.Volume.Weights;  %< (1 x nl)
                        u_avg = w * U;  %< (1 x ne)

                        %< Minimum nodal value per element (DOFs = node values)
                        u_min = min(U, [], 1);  %< (1 x ne)

                        %< Zhang-Shu theta
                        theta = ones(1, ne);
                        mask = (u_min < 0) & (u_avg > 0);
                        theta(mask) = u_avg(mask) ./ (u_avg(mask) - u_min(mask));
                        theta = min(1, theta);

                        %< Scale fluctuations; preserve cell average
                        U = u_avg + theta .* (U - u_avg);

                        %< Clamp negative cell averages
                        mask_neg = u_avg < 0;
                        U(:, mask_neg) = 0;

                    otherwise  %< 'cutoff'
                        U = max(0, U);
                end

                state.Dofs.U(:, 1) = U(:);
            end
        end
    end

end
