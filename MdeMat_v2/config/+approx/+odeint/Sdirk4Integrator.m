classdef Sdirk4Integrator < approx.odeint.DirkIntegrator
    % SDIRK4INTEGRATOR Fourth-order SDIRK integrator.
    %
    %   Sdirk4Integrator implements a fourth-order SDIRK method with five
    %   stages for solving ordinary differential equations of the form:
    %
    %   \f[
    %     \frac{\mathrm{d}(Mu)}{\mathrm{d}t} = Lu + S,
    %   \f]
    %
    %   where \f$M\f$ is a mass operator, \f$L\f$ is a linear operator, and
    %   \f$S\f$ is a source term. This method provides high-order accuracy
    %   with good stability properties for stiff problems. The diagonal
    %   elements of the Butcher tableau are identical.
    %
    %   The method uses \f$\gamma = 1/4\f$, which provides A-stability and
    %   good damping properties. The method requires five implicit solves
    %   per time step but achieves fourth-order accuracy, making it suitable
    %   for problems requiring high precision.
    %
    % Examples:
    %   % Basic usage
    %   integrator = approx.odeint.Sdirk4Integrator(1.0);
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
        ORDER = 4 % Accuracy order of the SDIRK4 method
    end

    methods
        function obj = Sdirk4Integrator(final)
            % SDIRK4INTEGRATOR Constructor for Sdirk4Integrator.
            %
            %   obj = Sdirk4Integrator(final) creates a fourth-order SDIRK
            %   integrator with the specified final time.
            %
            % Inputs:
            %   final - Final time for integration (positive scalar)
            %
            % Outputs:
            %   obj - Constructed Sdirk4Integrator object
            
            obj@approx.odeint.DirkIntegrator(5, final);
            
            %< Stage time coefficients
            obj.c = [1 / 4, 3 / 4, 11 / 20, 1 / 2, 1];
            
            %< Butcher tableau coefficient matrix (5x5)
            obj.A = [1 / 4, 0, 0, 0, 0; ...
                     1 / 2, 1 / 4, 0, 0, 0; ...
                     17 / 50, -1 / 25, 1 / 4, 0, 0; ...
                     371 / 1360, -137 / 2720, 15 / 544, 1 / 4, 0; ...
                     25 / 24, -49 / 48, 125 / 16, -85 / 12, 1 / 4];
            
            %< Stage weights for final update
            obj.b = [25 / 24, -49 / 48, 125 / 16, -85 / 12, 1 / 4];
        end
    end
end