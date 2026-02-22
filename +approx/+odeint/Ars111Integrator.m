classdef Ars111Integrator < approx.odeint.ImexrkIntegrator
    % ARS111INTEGRATOR First-order ARS(1,1,1) IMEX Runge-Kutta integrator.
    %
    %   Ars111Integrator implements the ARS(1,1,1) method, a first-order
    %   IMEX scheme for solving ordinary differential equations of the form:
    %
    %   \f[
    %     \frac{\mathrm{d}(Mu)}{\mathrm{d}t} = Lu + F(u) + S(t),
    %   \f]
    %
    %   where \f$M\f$ is a mass operator, \f$L\f$ is a linear operator
    %   (treated implicitly), \f$F(u)\f$ is a nonlinear operator (treated
    %   explicitly), and \f$S(t)\f$ is a source term (treated explicitly).
    %   This is the simplest IMEX method, combining forward Euler for the
    %   explicit part and backward Euler for the implicit part.
    %
    %   The method is A-stable for the implicit part and provides
    %   first-order accuracy. It's suitable for initial testing and
    %   problems where simplicity is preferred over high accuracy.
    %
    % See also:
    %   approx.odeint.ImexrkIntegrator, approx.odeint.Ars222Integrator,
    %   approx.odeint.Ars443Integrator

    properties (Constant)
        Order = 1 % Accuracy order of the ARS(1,1,1) method
    end
        
    methods
        function obj = Ars111Integrator(final)
            % ARS111INTEGRATOR Constructor for Ars111Integrator.
            %
            %   obj = Ars111Integrator(final) creates an ARS(1,1,1) IMEX
            %   integrator with the specified final time.

            arguments
                final {mustBeNumeric, mustBeNonnegative}
            end
            
            obj@approx.odeint.ImexrkIntegrator(2, final);
            
            %< ARS(1,1,1) coefficients
            obj.c = [0, 1];
            
            %< Explicit Butcher tableau (forward Euler structure)
            obj.AE = [0, 0; 1, 0];
            obj.bE = [1, 0];
            
            %< Implicit Butcher tableau (backward Euler structure)
            obj.AI = 1;
            obj.bI = 1;
        end
    end
end