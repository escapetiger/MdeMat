classdef Sdirk3Integrator < approx.odeint.DirkIntegrator
    % SDIRK3INTEGRATOR Third-order SDIRK integrator.
    %
    %   Sdirk3Integrator implements a third-order SDIRK method with three
    %   stages for solving ordinary differential equations of the form:
    %
    %   \f[
    %     \frac{\mathrm{d}(Mu)}{\mathrm{d}t} = Lu + S,
    %   \f]
    %
    %   where \f$M\f$ is a mass operator, \f$L\f$ is a linear operator, and
    %   \f$S\f$ is a source term. This method is A-stable and provides good
    %   stability properties for stiff problems. The diagonal elements of
    %   the Butcher tableau are identical.
    %
    %   The method uses \f$\gamma = 0.435866521508459\f$, which is chosen
    %   to optimize stability properties while maintaining third-order
    %   accuracy. The method requires three implicit solves per time step.
    %
    % See also:
    %   approx.odeint.DirkIntegrator, approx.odeint.Sdirk2Integrator,
    %   approx.odeint.Sdirk4Integrator

    properties (Constant)
        Order = 3 % Accuracy order of the SDIRK3 method
    end
        
    methods
        function obj = Sdirk3Integrator(final)
            % SDIRK3INTEGRATOR Constructor for Sdirk3Integrator.
            %
            %   obj = Sdirk3Integrator(final) creates a third-order SDIRK
            %   integrator with the specified final time.

            arguments
                final {mustBeNumeric, mustBeNonnegative}
            end
            
            obj@approx.odeint.DirkIntegrator(3, final);
            
            %< SDIRK3 coefficients for third-order accuracy
            gamma = 0.435866521508459;
            beta1 = -3 / 2 * gamma^2 + 4 * gamma - 1 / 4;
            beta2 = 3 / 2 * gamma^2 - 5 * gamma + 5 / 4;
            
            %< Stage time coefficients
            obj.c = [gamma, (1 + gamma) / 2, 1];
            
            %< Butcher tableau coefficient matrix
            obj.A = [gamma, 0, 0; 
                     (1 - gamma) / 2, gamma, 0; 
                     beta1, beta2, gamma];
            
            %< Stage weights for final update
            obj.b = [beta1, beta2, gamma];
        end
    end
end