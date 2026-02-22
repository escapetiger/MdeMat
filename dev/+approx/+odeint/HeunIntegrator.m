classdef HeunIntegrator < approx.odeint.ExrkIntegrator
    % HEUNINTEGRATOR Heun's integrator.
    %
    %   HeunIntegrator implements Heun's method, a second-order explicit
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
    %   Heun's method evaluates the derivative at \f$t\f$ and at \f$t +
    %   (2/3)\Delta t\f$, then takes a weighted average. The weights [1/4,
    %   3/4] give optimal second-order accuracy with this intermediate
    %   evaluation point.
    %
    % See also:
    %   approx.odeint.ExrkIntegrator, approx.odeint.Exrk4Integrator,
    %   approx.odeint.FeIntegrator
    
    properties (Constant)
        Order = 2 % Accuracy order of Heun's method
    end

    methods
        function obj = HeunIntegrator(final)
            % HEUNINTEGRATOR Constructor for HeunIntegrator.
            %
            %   obj = HeunIntegrator(final) creates a Heun's method
            %   integrator with the specified final time.

            arguments
                final {mustBeNumeric, mustBeNonnegative}
            end
            
            obj@approx.odeint.ExrkIntegrator(2, final);
            
            %< Heun's method coefficients
            obj.c = [0, 2/3];
            
            %< Butcher tableau coefficient matrix
            obj.A = [0, 0; 2/3, 0];
            
            %< Stage weights for final update
            obj.b = [1/4, 3/4];
        end
    end
end