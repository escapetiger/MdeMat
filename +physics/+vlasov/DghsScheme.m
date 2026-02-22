classdef DghsScheme < physics.vlasov.VlasovScheme
    % DGHScheme Discontinuous Galerkin Hermite spectral method for
    % Vlasov-Poisson equations.
    %
    %   DghsScheme implements a discontinuous Galerkin Hermite spectral
    %   method for solving Vlasov-Poisson-Fokker-Planck equations using
    %   discontinuous Galerkin spatial discretization and operator
    %   splitting temporal integration.
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
    %   physics.vlasov.VlasovScheme
    
    properties
        Multiplier % Assembly for multipliers
        Adjoint % Assembly for adjoint operators
        Source % Assembly for source terms

        Dual % Dual operators
        
        Poisson % Linear solver for the Poisson equation
    end
    
    properties (Constant)
        Name = 'DGHS' % Scheme name
    end
    
    methods
        function state = initialize(obj, state)
            % INITIALIZE Initialize the DG Hermite spectral scheme.
            %
            %   state = initialize(obj, state) sets up the DGHS scheme
            %   by assembling operators, fitting initial conditions,
            %   and preparing visualization components.
            
            arguments
                obj physics.vlasov.DghsScheme
                state physics.vlasov.HermiteState
            end
            
            %< Set up linear solver options
            linOpts = struct('nnzTh', inf, 'condTh', inf); % always use backslash
            obj.TDisc.reset(linOpts = linOpts);
            obj.Poisson = core.linalg.LinearSolver(nnzTh=Inf, condTh=Inf);
            
            %< Assembly for adjoint operators
            obj.Adjoint = approx.assembly.AdjointAssembly( ...
                state.XDisc, obj.Config.bcType, obj.Config.xPenaltyType);
            
            %< Assembly for multipliers
            obj.Multiplier = approx.assembly.MultiplierAssembly( ...
                state.XDisc);
            
            %< Assembly for source
            obj.Source = approx.assembly.SourceAssembly( ...
                state.XDisc);
            
            %< Set mass term
            obj.setMassTerm(state);
            
            %< Set transport term
            obj.setTransportTerm(state);
            
            %< Set initial condition for Hermite coefficients
            state.setDof('D', state.XDisc.fit(obj.Config.ic, args={state.NVDofs}));

            %< Set reference temperature
            state.setCoefficient('T0', obj.Config.T0);

            %< Set ion density coefficient
            state.setCoefficient('rhoi', state.XDisc.fit(obj.Config.rhoi));

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
                obj physics.vlasov.DghsScheme
                state physics.vlasov.HermiteState
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
            %   step using the DG Hermite spectral method with Strang
            %   splitting.
            
            arguments
                obj physics.vlasov.DghsScheme
                state physics.vlasov.HermiteState
            end
            
            dt = obj.TDisc.Timeline.StepSize;

            state = obj.stepLinear(state, dt / 2);

            state = obj.stepNonlinear(state, dt);

            state = obj.stepLinear(state, dt / 2);

            state = obj.applyDealiasingFilter(state);

            state = obj.stepPoisson(state);

            U = [state.Dofs.D, state.Dofs.P, state.Dofs.E];
            obj.TDisc.add(U(:));
        end
        
        function state = postStep(obj, state)
            % POSTSTEP Postprocess operations after time step.
            %
            %   state = postStep(obj, state) performs postprocessing
            %   including visualization updates and timeline advancement.
            
            arguments
                obj physics.vlasov.DghsScheme
                state physics.vlasov.HermiteState
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
                obj physics.vlasov.DghsScheme
                state physics.vlasov.HermiteState
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
    
    methods (Access = protected)
        function obj = setup(obj)
            % SETUP Setup default configuration options for DGHS schemes.
            
            %< Call parent setup
            obj = setup@physics.vlasov.VlasovScheme(obj);
            
            %< Add DGHS-specific required configuration options
            obj.addConfig('xBasisOrder', default = 1, ...
                validator = @(x) isnumeric(x) && isscalar(x) && (x > 0));
            obj.addConfig('xBasisType', default='nodal', ...
                validator=@(x) (ischar(x) || isstring(x)) && ismember(lower(x), {'modal', 'nodal'}));
            obj.addConfig('xBasisPattern', default = 'Q', ...
                validator = @(x) (ischar(x) || isstring(x)) && ismember(upper(x), {'P', 'Q'}));
            obj.addConfig('xPenaltyType', default = {'', 'right'; 'left', ''}, ...
                validator = @(x) iscell(x) && isequal(size(x), [2, 2]));
            obj.addConfig('vBBox', default = [-6, 6], ...
                validator = @(x) isnumeric(x) && isvector(x) && mod(length(x), 2) == 0);
            obj.addConfig('nh', default = 2, ...
                validator = @(x) isnumeric(x) && isscalar(x) && (x >= 0));
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
            obj = configure@physics.vlasov.VlasovScheme(obj);
            
            %< Determine spatial discretization identifier
            xId = sprintf('%s%d', upper(obj.Name), obj.Config.xBasisOrder);
            vId = sprintf('H%d', obj.Config.nh);
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
        
        function obj = plot(obj, state)
            % PLOT Generate visualization for current solution state.
            %
            %   obj = plot(obj, state) creates plots for the current
            %   solution state using default styling based on dataset
            %   types.
            
            visualizer = obj.Config.visualizer;

            xSpace = state.XDisc;
            dofs = state.Dofs;

            vBBox = obj.Config.vBBox;
            V = linspace(vBBox(1), vBBox(2), 512);

            visualizer.setRefXNodes(xSpace.Element.Geometry);

            if visualizer.IsEnabled
                visualizer.setPhyXNodes(xSpace.Mesh);
                
                [V0, n0] = obj.precompute(xSpace, dofs);
                
                while visualizer.Timeline.isBefore(obj.TDisc.Timeline.Next)
                    visualizer.update(xSpace, dofs, ...
                        timeline=obj.TDisc.Timeline, V0 = V0, n0 = n0);
                    visualizer.render();
                    visualizer.Timeline.advance();
                end
            end
            
            xRef = visualizer.Cache.RefXNodes;
            X = state.XDisc.Mesh.collocate(xRef);
            F = state.distribution(xRef, V);
            
            if obj.Config.verbose > 1
                strategy = physics.visual.Strategy2d();
                
                figure(2);
                ax = gca;
                cla(ax);
                hold(ax, 'on');
                imagesc(ax, X, V, F.', 'Interpolation', 'bilinear');
                
                axis(ax, 'xy', 'equal', 'tight');
                xl = xlabel(ax, 'x');
                yl = ylabel(ax, 'v');
                set(xl, 'FontSize', strategy.DefaultAxisLabelFontSize);
                set(yl, 'FontSize', strategy.DefaultAxisLabelFontSize);
                
                colormap(ax, strategy.ColorMap);
                colorbar(ax);
                strategy.setColorLimits(ax, F);
                
                hold(ax, 'off');
                title(sprintf('t = %.2f', obj.TDisc.Timeline.Next));
                drawnow;
                
                if obj.ShouldCheckpointing
                    obj.saveCkptFigure(gcf);
                end
            end
        end
    end
    
    methods (Access = private)
        function obj = setMassTerm(obj, state)
            % SETMASSTERM Set up mass term matrices.
            
            nh = state.NVDofs;
            ng = state.NXDofs;
            obj.M{1} = speye(ng * nh) * obj.Config.epsilon;
        end
        
        function obj = setTransportTerm(obj, state)
            % SETTRANSPORTTERM Set up transport term operators.
            
            %< Prestored operators
            D = obj.Adjoint.assembleMatrix(1);
            obj.Dual = cell(2, 2);
            obj.Dual{1, 2} = -D{1, 2}{1};
            obj.Dual{2, 1} = D{2, 1}{1};
            obj.Dual{1, 1} = -D{1, 1}{1};
            
            %< Linear stage: free transport + collision
            obj.addLinearizedTransport(state);
        end

        function obj = addLinearizedTransport(obj, state)
            % ADDLINEARIZEDTRANSPORT Add linearized transport terms.
            
            nh = state.NVDofs;
            ng = state.NXDofs;
            
            obj.L{1} = repmat({sparse(ng, ng)}, nh, nh);
            obj.addFreeTransport(state);
            obj.addCollision(state);
            obj.L{1} = core.linalg.block(obj.L{1});
        end
        
        function obj = addFreeTransport(obj, state)
            % ADDFREETRANSPORT Add free transport terms.
            
            nh = state.NVDofs;
            
            T0 = obj.Config.T0;
            obj.L{1}{1, 2} = sqrt(T0) * obj.Dual{1, 2};
            for k = 2:nh - 1
                obj.L{1}{k, k - 1} = -sqrt(k-1) * sqrt(T0) *  obj.Dual{2, 1};
                obj.L{1}{k, k + 1} = sqrt(k) * sqrt(T0) *  obj.Dual{1, 2};
            end
            obj.L{1}{nh, nh - 1} = -sqrt(nh - 1) * sqrt(T0) *  obj.Dual{2, 1};
        end
        
        function obj = addCollision(obj, state)
            % ADDCOLLISION Add collision terms.

            nh = state.NVDofs;
            T = obj.Multiplier.assembleMatrix(1);
            for k = 2:nh
                obj.L{1}{k, k} = -(k - 1) * T / obj.Config.tau0;
            end
        end

        function obj = setPoissonLhs(obj, state)
            % SETPOISSONLHS Set left-hand side matrix for Poisson equation.

            if obj.Poisson.IsPrecomputed
                return;
            end

            ng = state.NXDofs;
            lhs = repmat({sparse(ng, ng)}, 2, 2);
            lhs{1, 1} = obj.Dual{1, 1};
            lhs{1, 2} = -obj.Dual{1, 2};
            lhs{2, 1} = obj.Dual{2, 1};
            lhs{2, 2} = speye(ng, ng);
            lhs = core.linalg.block(lhs);
            c = obj.Source.assembleVector(@(x) ones(1, size(x, 2)));
            c = [c; sparse(ng, 1)];
            lhs = [lhs, c; c.', 0];
            obj.Poisson.setLhs(lhs);
            obj.Poisson.precompute();
        end
        
        function obj = setPoissonRhs(obj, state)
            % SETPOISSONRHS Set right-hand side vector for Poisson equation.

            ng = state.XDisc.NGlobalDofs;
            rhs = repmat({sparse(ng, 1)}, 2, 1);
            rhs{1} = state.Dofs.D(:, 1) - state.Coefs.rhoi;
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
            state.Dofs.P = W(:, 1);
            state.Dofs.E = W(:, 2);
        end
        
        function state = stepLinear(obj, state, dt)
            % STEPLINEAR Perform linear transport and collision step.

            nh = state.NVDofs;
            tDisc = obj.TDisc.getFactor(1);
            W0 = state.Dofs.D(:);
            tDisc.update(W0);
            W = tDisc.step(L=obj.L{1}, M=obj.M{1}, dt=dt);
            state.Dofs.D = reshape(W, [], nh);
        end
        
        function state = stepNonlinear(obj, state, dt)
            % STEPNONLINEAR Perform nonlinear electric field coupling step.

            state = obj.stepPoisson(state);

            tDisc = obj.TDisc.getFactor(2);
            switch class(tDisc)
                case 'approx.odeint.BeIntegrator'
                    state = obj.stepNonlinearBe(state, dt);
                case 'approx.odeint.ImpIntegrator'
                    state = obj.stepNonlinearImp(state, dt);
                case 'approx.odeint.Sdirk2Integrator'
                    state = obj.stepNonlinearSdirk2(state, dt);
                otherwise
                    core.except.assert(0, 'InvalidOdeInt', ...
                        'Unsupported ODE integrator for nonlinear step.');
            end
        end
        
        function state = stepNonlinearBe(obj, state, dt)
            % STEPNONLINEARBE Perform nonlinear electric field coupling
            % step using backward Euler method.

            nh = state.NVDofs;
            T0 = obj.Config.T0;

            D = state.Dofs.D;
            E = state.Dofs.E;
            f = @(x) reshape(state.XDisc.eval([], E), 1, []);
            T = obj.Multiplier.assembleMatrix(f);
            DNew = zeros(size(D));
            DNew(:, 1) = D(:, 1);
            k = 2:nh;
            DNew(:, k) = D(:, k) + dt * sqrt(k-1) / sqrt(T0) .* (T * DNew(:, k-1));
            state.Dofs.D = DNew;
        end

        function state = stepNonlinearImp(obj, state, dt)
            % STEPNONLINEARIMP Perform nonlinear electric field coupling
            % step using implicit midpoint method.

            nh = state.NVDofs;
            T0 = obj.Config.T0;
            epsilon = obj.Config.epsilon;

            D = state.Dofs.D;
            E = state.Dofs.E;
            f = @(x) reshape(state.XDisc.eval([], E), 1, []);
            T = obj.Multiplier.assembleMatrix(f);
            DNew = zeros(size(D));
            DNew(:, 1) = D(:, 1);
            k = 2:nh;
            DNew(:, k) = D(:, k);
            DNew(:, k) = DNew(:, k) + dt * sqrt(k-1) / sqrt(T0) .* (T * (D(:, k-1) + DNew(:, k-1)) / 2) / epsilon;
            state.Dofs.D = DNew;
        end

        function state = stepNonlinearSdirk2(obj, state, dt)
            % STEPNONLINEARSDIRK2 Perform nonlinear electric field coupling
            % step using second-order SDIRK method.

            nh = state.NVDofs;
            T0 = obj.Config.T0;
            epsilon = obj.Config.epsilon;

            D = state.Dofs.D;
            E = state.Dofs.E;
            f = @(x) reshape(state.XDisc.eval([], E), 1, []);
            T = obj.Multiplier.assembleMatrix(f);
            %< First stage
            gamma = 1 - sqrt(2)/2;
            a11 = gamma;
            D1 = zeros(size(D));
            D1(:, 1) = D(:, 1);
            k = 2:nh;
            D1(:, k) = D(:, k); 
            D1(:, k) = D1(:, k) + a11 * dt * sqrt(k-1) / sqrt(T0) .* (T * D1(:, k-1)) / epsilon;
            %< Second stage
            a21 = 1-gamma; a22 = gamma;
            D2 = zeros(size(D));
            D2(:, 1) = D(:, 1);
            k = 2:nh;
            D2(:, k) = D(:, k); 
            D2(:, k) = D2(:, k) + a21 * dt * sqrt(k-1) / sqrt(T0) .* (T * D1(:, k-1)) / epsilon;
            D2(:, k) = D2(:, k) + a22 * dt * sqrt(k-1) / sqrt(T0) .* (T * D2(:, k-1)) / epsilon;
            %< Final stage
            b1 = 1-gamma; b2 = gamma;
            DNew = zeros(size(D));
            DNew(:, 1) = D(:, 1);
            k = 2:nh;
            DNew(:, k) = D(:, k);
            DNew(:, k) = DNew(:, k) + b1 * dt * sqrt(k-1) / sqrt(T0) .* (T * D1(:, k-1)) / epsilon;
            DNew(:, k) = DNew(:, k) + b2 * dt * sqrt(k-1) / sqrt(T0) .* (T * D2(:, k-1)) / epsilon;
            state.Dofs.D = DNew;
        end

        function state = applyDealiasingFilter(obj, state)
            % APPLYDEALIASINGFILTER Apply dealiasing filter to Hermite coefficients.
            %
            %   state = applyDealiasingFilter(obj, state) applies the specified
            %   dealiasing filter to the Hermite coefficient matrix in state to reduce
            %   aliasing errors in spectral computations.

            if ~obj.Config.useFilter
                return;
            end

            nh = state.NVDofs;
            dt = obj.TDisc.Timeline.StepSize;
            beta = obj.Config.dealiasingBeta;
            gamma = obj.Config.dealiasingGamma;
            c = obj.Config.dealiasingCutoff;
            k = 1:nh;
            k0 = k - 1;
            eta = k0 / nh;

            switch lower(obj.Config.dealiasingType)
                case 'hou-li'
                    g = @(eta, dt) 1;
                case 'tc'
                    g = @(eta, dt) dt;
                case 'quasi-tc'
                    g = @(eta, dt) dt.^(1-eta.^gamma);
                otherwise
                    core.except.assert(0, 'InvalidDealiasingType', ...
                        ['Invalid dealiasing type: ', obj.Config.dealiasingType]);
            end

            sigma = (eta <= c) .* 1 + (eta > c) .* exp(-beta * (eta).^gamma) .* g(eta, dt);
            state.Dofs.D = state.Dofs.D .* sigma;
        end
    end
end