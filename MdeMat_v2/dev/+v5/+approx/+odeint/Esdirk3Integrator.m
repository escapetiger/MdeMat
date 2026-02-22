classdef Esdirk3Integrator < approx.odeint.DirkIntegrator
    % ESDIRK3INTEGRATOR Third-order ESDIRK integrator.
    %
    %   Esdirk3Integrator implements a third-order ESDIRK method with three
    %   stages for solving ordinary differential equations of the form:
    %
    %   \f[
    %     \frac{\mathrm{d}(Mu)}{\mathrm{d}t} = Lu + S,
    %   \f]
    %
    %   where \f$M\f$ is a mass operator, \f$L\f$ is a linear operator, and
    %   \f$S\f$ is a source term. This method features an explicit first
    %   stage (\f$c_1 = 0\f$) followed by two implicit stages, providing
    %   computational efficiency while maintaining good stability properties.
    %
    %   The parameter \f$\gamma = (3 + \sqrt{3})/6 \approx 0.7886751346\f$
    %   ensures third-order accuracy and A-stability. The explicit first
    %   stage reduces computational cost compared to fully implicit methods.
    %
    % Examples:
    %   % Basic usage
    %   integrator = approx.odeint.Esdirk3Integrator(1.0);
    %   U_new = integrator.step(L, S, M);
    %   integrator.update(U_new);
    %   integrator.advance();
    %
    %   % Pure diffusion (only linear term)
    %   U_new = integrator.step(L, [], []);
    %
    %   % Pure forcing (only source term)
    %   U_new = integrator.step([], S, []);
    %
    %   % With mass matrix
    %   U_new = integrator.step(L, S, M);
    %
    % See also:
    %   approx.odeint.DirkIntegrator, approx.odeint.Sdirk2Integrator,
    %   approx.odeint.Sdirk3Integrator
    
    properties (Constant)
        ORDER = 3 % Accuracy order of the ESDIRK3 method
    end

    methods
        function obj = Esdirk3Integrator(final)
            % ESDIRK3INTEGRATOR Constructor for Esdirk3Integrator.
            %
            %   obj = Esdirk3Integrator(final) creates a third-order ESDIRK
            %   integrator with the specified final time.
            %
            % Inputs:
            %   final - Final time for integration (positive scalar)
            %
            % Outputs:
            %   obj - Constructed Esdirk3Integrator object
            
            obj@approx.odeint.DirkIntegrator(3, final);
            
            %< ESDIRK3 coefficients for third-order accuracy
            c3 = 1;
            gamma = (3 + sqrt(3)) / 6;
            a32 = c3 * (c3 - 2 * gamma) / (4 * gamma);
            b2 = (-2 + 3 * c3) / (12 * (c3 - 2 * gamma) * gamma);
            b3 = (1 - 3 * gamma) / (3 * c3 * (c3 - 2 * gamma));
            
            %< Stage time coefficients
            obj.c = [0, 2 * gamma, 1];
            
            %< Butcher tableau coefficient matrix
            obj.A = [0, 0, 0; 
                     gamma, gamma, 0; 
                     1 - gamma - a32, a32, gamma];
            
            %< Stage weights for final update
            obj.b = [1 - b2 - b3, b2, b3];
        end
    end
end