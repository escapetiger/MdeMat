classdef StrangIntegrator < approx.odeint.SeparableOdeIntegrator
    % STRANGINTEGRATOR Second-order Strang operator splitting.
    %
    %   StrangIntegrator implements the second-order Strang operator
    %   splitting method for solving ordinary differential equations of
    %   the form:
    %
    %   \f[
    %     \frac{\mathrm{d}(Mu)}{\mathrm{d}t} = \sum_{i=1}^{n} (L_i u +
    %     F_i(u) + S_i)
    %   \f]
    %
    %   by decomposing it into n separate sub-problems:
    %
    %   \f[
    %     \frac{\mathrm{d}(Mu)}{\mathrm{d}t} = L_i u + F_i(u) + S_i, \quad
    %     i = 1, 2, \ldots, n
    %   \f]
    %
    %   The Strang method uses the symmetric splitting formula for n
    %   operators: 
    %
    %   \f[
    %     e^{(A_1+A_2+\cdots+A_n)\Delta t} \approx e^{A_1\Delta t/2}
    %     e^{A_2\Delta t/2} \cdots e^{A_n\Delta t/2} e^{A_n\Delta t/2}
    %     \cdots e^{A_2\Delta t/2} e^{A_1\Delta t/2}
    %   \f]
    %
    %   This achieves second-order accuracy with global error \f$O(\Delta
    %   t^2)\f$. The symmetric structure causes second-order commutator
    %   terms to cancel exactly, providing superior accuracy compared to
    %   first-order Lie-Trotter splitting.
    %
    %   The splitting sequence applies operators in forward order with
    %   half time steps, then in reverse order with half time steps.
    %
    % Examples:
    %   % Reaction-diffusion: implicit diffusion + explicit reaction
    %   integrators = {approx.odeint.BeIntegrator(1.0), ...
    %                  approx.odeint.FeIntegrator(1.0)};
    %   splitter = approx.odeint.StrangIntegrator(integrators);
    %   
    %   % Step with operators [diffusion, reaction]
    %   L = {diffusionMatrix, []};
    %   F = {[], @(u,t) reactionTerm(u)};
    %   S = {[], []};
    %   U_new = splitter.step(L, F, S, M);
    %   
    %   % Three-way Strang split: A_dt/2, B_dt/2, C_dt, B_dt/2, A_dt/2
    %   integrators = {approx.odeint.BeIntegrator(1.0), ...
    %                  approx.odeint.FeIntegrator(1.0), ...
    %                  approx.odeint.FeIntegrator(1.0)};
    %   L = {diffusionMatrix, advectionMatrix, []};
    %   F = {[], [], @(u,t) reactionTerm(u)};
    %   S = {sourceVector, [], []};
    %
    % Notes:
    %   For n operators, the method applies 2n-1 sub-steps per time step.
    %   The middle operator (when n is odd) gets applied once with full
    %   time step, while others get applied twice with half time steps.
    %
    % See also:
    %   approx.odeint.SeparableOdeIntegrator, 
    %   approx.odeint.LieTrotterIntegrator,
    %   approx.odeint.BeIntegrator, 
    %   approx.odeint.FeIntegrator

    properties (Constant)
        ORDER = 2 % Accuracy order of the Strang method
    end

    methods
        function U = step(obj, L, F, S, M)
            % STEP Advance one time step using Strang splitting.
            %
            %   U = step(obj, L, F, S, M) performs one time step using
            %   the symmetric Strang splitting: forward half-steps, then
            %   backward half-steps.
            %
            % Inputs:
            %   obj - The StrangIntegrator object
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

            n = obj.nFactors;
            U0 = obj.getFactor(1).U0{1};
            dt = obj.timeline.dt;

            %< Apply forward half-steps: A1_dt/2, A2_dt/2, ..., A_{n-1}_dt/2
            for i = 1:(n-1)
                U = obj.applyFactorStep(i, dt/2, L{i}, F{i}, S{i}, M);
                obj.getFactor(i+1).U0{1} = U;
            end
            
            %< Apply middle full step: An_dt
            U = obj.applyFactorStep(n, dt, L{n}, F{n}, S{n}, M);
            obj.getFactor(n-1).U0{1} = U;
            
            %< Apply backward half-steps: A_{n-1}_dt/2, ..., A2_dt/2, A1_dt/2
            for i = (n-1):-1:1
                U = obj.applyFactorStep(i, dt/2, L{i}, F{i}, S{i}, M);
                if i > 1
                    obj.getFactor(i-1).U0{1} = U;
                end
            end

            %< Recovery
            for i = 1:n
                obj.getFactor(i).U0{1} = U0;
            end
        end
    end
end