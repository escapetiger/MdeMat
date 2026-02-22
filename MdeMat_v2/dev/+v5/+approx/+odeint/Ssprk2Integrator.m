classdef Ssprk2Integrator < approx.odeint.ExrkIntegrator
    % SSPRK2INTEGRATOR SSPRK2 integrator.
    %
    %   Ssprk2Integrator implements SSPRK2 method, a second-order explicit
    %   Runge-Kutta scheme for solving ordinary differential equations of
    %   the form:
    %
    %   \f[
    %     \frac{\mathrm{d}(Mu)}{\mathrm{d}t} = Lu + F(u) + S(t),
    %   \f]
    %
    %   where \f$M\f$ is a mass operator, \f$L\f$ is a linear operator,
    %   \f$F(u)\f$ is a nonlinear operator, and \f$S(t)\f$ is a source
    %   term. This method provides second-order accuracy with good
    %   stability properties. It's also known as the improved Euler method
    %   or RK2.
    %
    %   SSPRK2 method evaluates the derivative at \f$t\f$ and at \f$t +
    %   (2/3)\Delta t\f$, then takes a weighted average. The weights [1/4,
    %   3/4] give optimal second-order accuracy with this intermediate
    %   evaluation point.
    %
    % Examples:
    %   % Basic usage with all terms
    %   integrator = approx.odeint.Ssprk2Integrator(1.0);
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
    %   approx.odeint.ExrkIntegrator, approx.odeint.Exrk4Integrator,
    %   approx.odeint.FeIntegrator
    
    properties (Constant)
        ORDER = 2 % Accuracy order of SSPRK2 method
    end

    methods
        function obj = Ssprk2Integrator(final)
            % SSPRK2INTEGRATOR Constructor for Ssprk2Integrator.
            %
            %   obj = Ssprk2Integrator(final) creates a SSPRK2 method
            %   integrator with the specified final time.
            %
            % Inputs:
            %   final - Final time for integration (positive scalar)
            %
            % Outputs:
            %   obj - Constructed Ssprk2Integrator object
            
            obj@approx.odeint.ExrkIntegrator(2, final);
            
            %< SSPRK2 method coefficients
            obj.c = [0, 1];
            
            %< Butcher tableau coefficient matrix
            obj.A = [0, 0; 1, 0];
            
            %< Stage weights for final update
            obj.b = [1/2, 1/2];
        end
    end
end