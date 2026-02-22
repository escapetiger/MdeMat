classdef LieTrotterIntegrator < approx.odeint.SeparableOdeIntegrator
    % LIETROTTERINTEGRATOR First-order Lie-Trotter operator splitting.
    %
    %   LieTrotterIntegrator implements the first-order Lie-Trotter
    %   operator splitting method for solving ordinary differential
    %   equations of the form:
    %
    %   \f[
    %     \frac{\mathrm{d}(Mu)}{\mathrm{d}t} = \sum_{i=1}^{n} (L_i u + F_i(u) + S_i)
    %   \f]
    %
    %   by decomposing it into n separate sub-problems:
    %   \f[
    %     \frac{\mathrm{d}(Mu)}{\mathrm{d}t} = L_i u + F_i(u) + S_i, \quad i = 1, 2, \ldots, n
    %   \f]
    %
    %   The Lie-Trotter method uses the splitting formula:
    %   \f[
    %     e^{(A_1+A_2+\cdots+A_n)\Delta t} \approx e^{A_1\Delta t} e^{A_2\Delta t} \cdots e^{A_n\Delta t}
    %   \f]
    %
    %   This achieves first-order accuracy with global error O(Δt).
    %   Each factor integrator solves its respective sub-problem using the
    %   most appropriate method (explicit, implicit, or IMEX).
    %
    % Examples:
    %   % Reaction-diffusion: implicit diffusion + explicit reaction
    %   integrators = {approx.odeint.BeIntegrator(1.0), ...
    %                  approx.odeint.FeIntegrator(1.0)};
    %   splitter = approx.odeint.LieTrotterIntegrator(integrators);
    %   
    %   % Step with operators [diffusion, reaction]
    %   L = {diffusionMatrix, []};
    %   F = {[], @(u,t) reactionTerm(u)};
    %   S = {[], []};
    %   U_new = splitter.step(L, F, S, M);
    %   
    %   % Three-way split: diffusion + advection + reaction
    %   integrators = {approx.odeint.BeIntegrator(1.0), ...
    %                  approx.odeint.FeIntegrator(1.0), ...
    %                  approx.odeint.FeIntegrator(1.0)};
    %   L = {diffusionMatrix, advectionMatrix, []};
    %   F = {[], [], @(u,t) reactionTerm(u)};
    %   S = {sourceVector, [], []};
    %
    % See also:
    %   approx.odeint.SeparableOdeIntegrator, 
    %   approx.odeint.StrangIntegrator,
    %   approx.odeint.BeIntegrator, 
    %   approx.odeint.FeIntegrator

    properties (Constant)
        ORDER = 1 % Accuracy order of the Lie-Trotter method
    end

    methods
        function U = step(obj, L, F, S, M)
            % STEP Advance one time step using Lie-Trotter splitting.
            %
            %   U = step(obj, L, F, S, M) performs one time step by
            %   sequentially applying each factor integrator to its
            %   corresponding sub-problem.
            %
            % Inputs:
            %   obj - The LieTrotterIntegrator object
            %   L - Cell array of linear operators {L1, L2, ..., Ln}
            %       Each Li can be matrix, function, or empty
            %   F - Cell array of nonlinear operators {F1, F2, ..., Fn}
            %       Each Fi can be function or empty
            %   S - Cell array of source terms {S1, S2, ..., Sn}
            %       Each Si can be vector, function, or empty
            %   M - Mass operator (matrix, function, or empty)
            %       Same mass operator used for all sub-problems
            %
            % Outputs:
            %   U - Solution at the next time step

            U0 = obj.getFactor(1).U0{1};
            dt = obj.timeline.dt;

            for i = 1:obj.nFactors
                U = obj.applyFactorStep(i, dt, L{i}, F{i}, S{i}, M);
                if i < obj.nFactors
                    obj.getFactor(i+1).U0{1} = U;
                end
            end
            
            for i = 1:obj.nFactors
                obj.getFactor(i).U0{1} = U0;
            end
        end
    end
end