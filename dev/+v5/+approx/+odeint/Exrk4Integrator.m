classdef Exrk4Integrator < approx.odeint.ExrkIntegrator
    % EXRK4INTEGRATOR Classic fourth-order Runge-Kutta integrator.
    %
    %   Exrk4Integrator implements the classic RK4 method, a fourth-order
    %   explicit scheme for solving ordinary differential equations of the
    %   form:
    %
    %   \f[
    %     \frac{\mathrm{d}(Mu)}{\mathrm{d}t} = Lu + F(u) + S(t),
    %   \f]
    %
    %   where \f$M\f$ is a mass operator, \f$L\f$ is a linear operator,
    %   \f$F(u)\f$ is a nonlinear operator, and \f$S(t)\f$ is a source
    %   term. This is one of the most widely used explicit methods,
    %   providing excellent accuracy with four stages per time step. It
    %   offers a good balance of accuracy and computational efficiency.
    %
    %   The classic RK4 method evaluates derivatives at the beginning, two
    %   midpoints, and end of the time interval, then combines them with
    %   weights [1/6, 1/3, 1/3, 1/6] to achieve fourth-order accuracy.
    %   It's the gold standard for non-stiff problems.
    %
    % Examples:
    %   % Basic usage with all terms
    %   integrator = approx.odeint.Exrk4Integrator(1.0);
    %   U_new = integrator.step(L, F, S, M);
    %   integrator.update(U_new);
    %   integrator.advance();
    %
    %   % Pure diffusion (only linear term)
    %   U_new = integrator.step(L, [], [], []);
    %
    %   % Pure nonlinear dynamics
    %   U_new = integrator.step([], F, [], []);
    %
    %   % Linear system with forcing
    %   U_new = integrator.step(L, [], S, []);
    %
    % See also:
    %   approx.odeint.ExrkIntegrator, approx.odeint.HeunIntegrator,
    %   approx.odeint.Ssprk3Integrator
    
    properties (Constant)
        ORDER = 4 % Accuracy order of the classic RK4 method
    end

    methods
        function obj = Exrk4Integrator(final)
            % EXRK4INTEGRATOR Constructor for Exrk4Integrator.
            %
            %   obj = Exrk4Integrator(final) creates a classic fourth-order
            %   Runge-Kutta integrator with the specified final time.
            %
            % Inputs:
            %   final - Final time for integration (positive scalar)
            %
            % Outputs:
            %   obj - Constructed Exrk4Integrator object
            
            obj@approx.odeint.ExrkIntegrator(4, final);
            
            %< Classic RK4 coefficients
            obj.c = [0, 1/2, 1/2, 1];
            
            %< Butcher tableau coefficient matrix
            obj.A = [0, 0, 0, 0; ...
                     1/2, 0, 0, 0; ...
                     0, 1/2, 0, 0; ...
                     0, 0, 1, 0];
            
            %< Stage weights for final update
            obj.b = [1/6, 1/3, 1/3, 1/6];
        end
    end
end