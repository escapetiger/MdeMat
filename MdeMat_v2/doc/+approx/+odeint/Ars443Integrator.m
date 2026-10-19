classdef Ars443Integrator < approx.odeint.ImexrkIntegrator
    % ARS443INTEGRATOR Third-order ARS(4,4,3) IMEX Runge-Kutta integrator.
    %
    %   Ars443Integrator implements the ARS(4,4,3) method, a third-order
    %   IMEX scheme for solving ordinary differential equations of the
    %   form:
    %
    %   \f[
    %     \frac{\mathrm{d}(Mu)}{\mathrm{d}t} = Lu + F(u) + S(t),
    %   \f]
    %
    %   where \f$M\f$ is a mass operator, \f$L\f$ is a linear operator
    %   (treated implicitly), \f$F(u)\f$ is a nonlinear operator (treated
    %   explicitly), and \f$S(t)\f$ is a source term (treated explicitly).
    %   This method provides high accuracy with good stability properties,
    %   making it suitable for problems requiring precise solutions with
    %   moderate computational cost.
    %
    %   This method uses 5 stages to achieve third-order accuracy. The
    %   method maintains A-stability for the implicit part and provides
    %   excellent accuracy for smooth solutions. The explicit part is
    %   designed to minimize storage requirements and computational
    %   overhead.
    %
    % Examples:
    %   % Basic usage
    %   integrator = approx.odeint.Ars443Integrator(1.0);
    %   U_new = integrator.step(L, F, S, M);
    %   integrator.update(U_new);
    %   integrator.advance();
    %
    %   % Pure implicit diffusion
    %   U_new = integrator.step(L, [], [], []);
    %
    %   % Mixed implicit-explicit problem
    %   U_new = integrator.step(L, F, [], []);
    %
    %   % Explicit nonlinear dynamics with forcing
    %   U_new = integrator.step([], F, S, []);
    %
    % See also:
    %   approx.odeint.ImexrkIntegrator, approx.odeint.Ars111Integrator,
    %   approx.odeint.Ars222Integrator
    
    properties (Constant)
        ORDER = 3 % Accuracy order of the ARS(4,4,3) method
    end

    methods
        function obj = Ars443Integrator(final)
            % ARS443INTEGRATOR Constructor for Ars443Integrator.
            %
            %   obj = Ars443Integrator(final) creates an ARS(4,4,3) IMEX
            %   integrator with the specified final time.
            %
            % Inputs:
            %   final - Final time for integration (positive scalar)
            %
            % Outputs:
            %   obj - Constructed Ars443Integrator object
            
            obj@approx.odeint.ImexrkIntegrator(5, final);
            
            %< ARS(4,4,3) coefficients for third-order accuracy
            obj.c = [0, 1/2, 2/3, 1/2, 1];
            
            %< Explicit Butcher tableau (5x5 matrix)
            obj.AE = [0, 0, 0, 0, 0; ...
                      1/2, 0, 0, 0, 0; ...
                      11/18, 1/18, 0, 0, 0; ...
                      5/6, -5/6, 1/2, 0, 0; ...
                      1/4, 7/4, 3/4, -7/4, 0];
            obj.bE = [1/4, 7/4, 3/4, -7/4, 0];
            
            %< Implicit Butcher tableau (4x4 matrix)
            obj.AI = [1/2, 0, 0, 0; ...
                      1/6, 1/2, 0, 0; ...
                      -1/2, 1/2, 1/2, 0; ...
                      3/2, -3/2, 1/2, 1/2];
            obj.bI = [3/2, -3/2, 1/2, 1/2];
        end
    end
end