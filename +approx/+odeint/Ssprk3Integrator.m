classdef Ssprk3Integrator < approx.odeint.ExrkIntegrator
    % SSPRK3INTEGRATOR Third-order SSPRK integrator.
    %
    %   Ssprk3Integrator implements the SSPRK3 method, a third-order
    %   explicit scheme for solving ordinary differential equations of the
    %   form:
    %
    %   \f[
    %     \frac{\mathrm{d}(Mu)}{\mathrm{d}t} = Lu + F(u) + S(t),
    %   \f]
    %
    %   where \f$M\f$ is a mass operator, \f$L\f$ is a linear operator,
    %   \f$F(u)\f$ is a nonlinear operator, and \f$S(t)\f$ is a source
    %   term. This method is designed to preserve strong stability
    %   properties, making it particularly suitable for hyperbolic
    %   conservation laws and problems where monotonicity preservation is
    %   important.
    %
    %   SSPRK3 maintains the total variation diminishing (TVD) property
    %   under appropriate CFL conditions. The method uses three stages with
    %   specific coefficients designed to preserve stability properties of
    %   the spatial discretization.
    %
    %   SSPRK3 is particularly effective for problems involving
    %   discontinuous solutions or steep gradients where maintaining
    %   monotonicity and avoiding spurious oscillations is critical.
    %
    % See also:
    %   approx.odeint.ExrkIntegrator, approx.odeint.Exrk4Integrator,
    %   approx.odeint.HeunIntegrator
    
    properties (Constant)
        Order = 3 % Accuracy order of the SSPRK3 method
    end

    methods
        function obj = Ssprk3Integrator(final)
            % SSPRK3INTEGRATOR Constructor for Ssprk3Integrator.
            %
            %   obj = Ssprk3Integrator(final) creates a strong stability
            %   preserving third-order Runge-Kutta integrator with the
            %   specified final time.

            arguments
                final {mustBeNumeric, mustBeNonnegative}
            end
            
            obj@approx.odeint.ExrkIntegrator(3, final);
            
            %< SSPRK3 coefficients for strong stability preservation
            obj.c = [0, 1, 1/2];
            
            %< Butcher tableau coefficient matrix
            obj.A = [0, 0, 0; 
                     1, 0, 0; 
                     1/4, 1/4, 0];
            
            %< Stage weights for final update
            obj.b = [1/6, 1/6, 2/3];
        end
    end
end