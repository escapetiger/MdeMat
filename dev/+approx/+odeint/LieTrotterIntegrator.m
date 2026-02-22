classdef LieTrotterIntegrator < approx.odeint.SeparableOdeIntegrator
    % LIETROTTERINTEGRATOR First-order Lie-Trotter operator splitting.
    %
    %   LieTrotterIntegrator implements the first-order Lie-Trotter
    %   operator splitting method for solving ordinary differential
    %   equations of the form:
    %
    %   \f[
    %     \frac{\mathrm{d}(Mu)}{\mathrm{d}t} = \sum_{i=1}^{n} (L_i u +
    %     F_i(u) + S_i)
    %   \f]
    %
    %   by decomposing it into n separate sub-problems:
    %   \f[
    %     \frac{\mathrm{d}(Mu)}{\mathrm{d}t} = L_i u + F_i(u) + S_i, \quad
    %     i = 1, 2, \ldots, n
    %   \f]
    %
    %   The Lie-Trotter method uses the splitting formula:
    %   \f[
    %     e^{(A_1+A_2+\cdots+A_n)\Delta t} \approx e^{A_1\Delta t}
    %     e^{A_2\Delta t} \cdots e^{A_n\Delta t}
    %   \f]
    %
    %   This achieves first-order accuracy with global error \f$O(\Delta
    %   t)\f$. Each factor integrator solves its respective sub-problem
    %   using the most appropriate method (explicit, implicit, or IMEX).
    %
    % See also:
    %   approx.odeint.SeparableOdeIntegrator,
    %   approx.odeint.StrangIntegrator,
    %   approx.odeint.BeIntegrator,
    %   approx.odeint.FeIntegrator

    properties (Constant)
        Order = 1 % Accuracy order of the Lie-Trotter method
    end

    methods
        function U = step(obj, options)
            % STEP Advance one time step using Lie-Trotter splitting.
            %
            %   U = step(obj, L=L, F=F, S=S, M=M) performs one time step by
            %   sequentially applying each factor integrator to its
            %   corresponding regular sub-problem.
            %
            %   U = step(obj, L=L, F=F, S=S, M=M, P=P, A=A, B=B, C=C)
            %   performs one time step by sequentially applying each factor
            %   integrator to its corresponding augmented sub-problem.

            arguments
                obj approx.odeint.LieTrotterIntegrator
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

            U0 = obj.getFactor(1).U0{1};
            dt = obj.Timeline.StepSize;

            for i = 1:obj.NFactors
                L = obj.getEntry(i, options.L);
                F = obj.getEntry(i, options.F);
                S = obj.getEntry(i, options.S);
                M = obj.getEntry(i, options.M);
                P = obj.getEntry(i, options.P);
                A = obj.getEntry(i, options.A);
                B = obj.getEntry(i, options.B);
                C = obj.getEntry(i, options.C);
                fn = obj.getEntry(i, options.fn);
                U = obj.factorStep(i, dt = dt, L = L, F = F, ...
                    S = S, M = M, P = P, A = A, B = B, C = C, fn = fn);

                if i < obj.NFactors
                    obj.getFactor(i+1).add(U);
                end
            end

            for i = 1:obj.NFactors
                obj.getFactor(i).add(U0);
            end
        end
    end
end
