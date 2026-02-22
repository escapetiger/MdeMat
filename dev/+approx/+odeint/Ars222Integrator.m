classdef Ars222Integrator < approx.odeint.ImexrkIntegrator
    % ARS222INTEGRATOR Second-order ARS(2,2,2) IMEX Runge-Kutta integrator.
    %
    %   Ars222Integrator implements the ARS(2,2,2) method, a second-order
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
    %   This method provides a good balance of accuracy and computational
    %   efficiency for problems with both stiff linear and nonlinear
    %   components.
    %
    %   The method uses \f$\gamma = 1 - \sqrt{2}/2 \approx 0.2928932188\f$
    %   for optimal stability properties. This parameter ensures
    %   A-stability for the implicit part while maintaining second-order
    %   accuracy.
    %
    % See also:
    %   approx.odeint.ImexrkIntegrator, approx.odeint.Ars111Integrator,
    %   approx.odeint.Ars443Integrator
    
    properties (Constant)
        Order = 2 % Accuracy order of the ARS(2,2,2) method
    end

    methods
        function obj = Ars222Integrator(final)
            % ARS222INTEGRATOR Constructor for Ars222Integrator.
            %
            %   obj = Ars222Integrator(final) creates an ARS(2,2,2) IMEX
            %   integrator with the specified final time.

            arguments
                final {mustBeNumeric, mustBeNonnegative}
            end
            
            obj@approx.odeint.ImexrkIntegrator(3, final);
            
            %< ARS(2,2,2) coefficients with optimal stability parameter
            gamma = 1 - sqrt(2) / 2;
            delta = 1 - 1 / (2 * gamma);
            
            obj.c = [0, gamma, 1];
            
            %< Explicit Butcher tableau
            obj.AE = [0, 0, 0; ...
                      gamma, 0, 0; ...
                      delta, 1 - delta, 0];
            obj.bE = [delta, 1 - delta, 0];
            
            %< Implicit Butcher tableau
            obj.AI = [gamma, 0; ...
                      1 - gamma, gamma];
            obj.bI = [1 - gamma, gamma];
        end
    end
end