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
    %     e^{A_2\Delta t/2} \cdots e^{A_{n-1}\Delta t/2} e^{A_n\Delta t}
    %     e^{A_{n-1}\Delta t/2} \cdots e^{A_2\Delta t/2} e^{A_1\Delta t/2}
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
        Order = 2 % Accuracy order of the Strang method
    end

    methods
        function U = step(obj, options)
            % STEP Advance one time step using Strang splitting.
            %
            %   U = step(obj, L=L, F=F, S=S, M=M) performs one time step
            %   using the symmetric Strang splitting: forward half-steps,
            %   then backward half-steps.
            %
            %   U = step(obj, L=L, F=F, S=S, M=M, P=P, A=A, B=B, C=C)
            %   performs one time step by sequentially applying each factor
            %   integrator to its corresponding augmented sub-problem.

            arguments
                obj approx.odeint.StrangIntegrator
                options.L = {} % Linear operators (cell array)
                options.F = {} % Nonlinear operators (cell array)
                options.S = {} % Source terms (cell array)
                options.M = {} % Mass operators (cell array)
                options.P = {} % Coupling operators (cell array)
                options.A = {} % Constraint matrices (cell array)
                options.B = {} % Constraint coupling matrices (cell array)
                options.C = {} % Constraint source terms (cell array)
                options.fn = {} % Custom step functions for each factor (cell array)
            end

            n = obj.NFactors;
            U0 = obj.getFactor(1).U0{1};
            dt = obj.Timeline.StepSize;

            %< Apply forward half-steps: A1_dt/2, A2_dt/2, ..., A_{n-1}_dt/2
            for i = 1:(n-1)
                L = obj.getEntry(i, options.L);
                F = obj.getEntry(i, options.F);
                S = obj.getEntry(i, options.S);
                M = obj.getEntry(i, options.M);
                P = obj.getEntry(i, options.P);
                A = obj.getEntry(i, options.A);
                B = obj.getEntry(i, options.B);
                C = obj.getEntry(i, options.C);
                fn = obj.getEntry(i, options.fn);
                U = obj.factorStep(i, dt = dt / 2, L = L, F = F, ...
                    S = S, M = M, P = P, ...
                    A = A, B = B, C = C, fn = fn);
                obj.getFactor(i+1).add(U);
            end

            %< Apply middle full step: An_dt
            L = obj.getEntry(n, options.L);
            F = obj.getEntry(n, options.F);
            S = obj.getEntry(n, options.S);
            M = obj.getEntry(n, options.M);
            P = obj.getEntry(n, options.P);
            A = obj.getEntry(n, options.A);
            B = obj.getEntry(n, options.B);
            C = obj.getEntry(n, options.C);
            fn = obj.getEntry(n, options.fn);
            U = obj.factorStep(n, dt = dt, L = L, F = F, ...
                S = S, M = M, P = P, ...
                A = A, B = B, C = C, fn = fn);
            obj.getFactor(n-1).add(U);

            %< Apply backward half-steps: A_{n-1}_dt/2, ..., A2_dt/2, A1_dt/2
            for i = (n-1):-1:1
                L = obj.getEntry(i, options.L);
                F = obj.getEntry(i, options.F);
                S = obj.getEntry(i, options.S);
                M = obj.getEntry(i, options.M);
                P = obj.getEntry(i, options.P);
                A = obj.getEntry(i, options.A);
                B = obj.getEntry(i, options.B);
                C = obj.getEntry(i, options.C);
                fn = obj.getEntry(i, options.fn);
                U = obj.factorStep(i, dt = dt / 2, L = L, F = F, ...
                    S = S, M = M, P = P, ...
                    A = A, B = B, C = C, fn = fn);
                if i > 1
                    obj.getFactor(i-1).add(U);
                end
            end

            %< Recovery
            for i = 1:n
                obj.getFactor(i).add(U0);
            end
        end
    end
end
