classdef Sdirk2Integrator < approx.odeint.DirkIntegrator
    % SDIRK2INTEGRATOR Second-order SDIRK integrator.
    %
    %   Sdirk2Integrator implements a second-order SDIRK method with two
    %   stages for solving ordinary differential equations of the form:
    %
    %   \f[
    %     \frac{\mathrm{d}(Mu)}{\mathrm{d}t} = Lu + S,
    %   \f]
    %
    %   where \f$M\f$ is a mass operator, \f$L\f$ is a linear operator, and
    %   \f$S\f$ is a source term. This method is L-stable and A-stable,
    %   making it suitable for stiff problems. The diagonal elements of the
    %   Butcher tableau are identical.
    %
    %   The method uses the parameter \f$\gamma = 1 - \sqrt{2}/2 \approx
    %   0.2928932188\f$, which ensures L-stability and optimal damping
    %   properties for stiff systems.
    %
    % Examples:
    %   % Basic usage
    %   integrator = approx.odeint.Sdirk2Integrator(1.0);
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
    %   approx.odeint.DirkIntegrator, approx.odeint.Sdirk3Integrator,
    %   approx.odeint.Esdirk3Integrator
        
    properties (Constant)
        ORDER = 2 % Accuracy order of the SDIRK2 method
    end
   
    methods
        function obj = Sdirk2Integrator(final)
            % SDIRK2INTEGRATOR Constructor for Sdirk2Integrator.
            %
            %   obj = Sdirk2Integrator(final) creates a second-order SDIRK
            %   integrator with the specified final time.
            %
            % Inputs:
            %   final - Final time for integration (positive scalar)
            %
            % Outputs:
            %   obj - Constructed Sdirk2Integrator object
            
            obj@approx.odeint.DirkIntegrator(2, final);
            
            %< SDIRK2 coefficients with L-stable parameter
            gamma = 1 - sqrt(2) / 2;
            
            %< Stage time coefficients
            obj.c = [gamma, 1];
            
            %< Butcher tableau coefficient matrix
            obj.A = [gamma, 0; 
                     1 - gamma, gamma];
            
            %< Stage weights for final update
            obj.b = [1 - gamma, gamma];
        end
    end
end