classdef EmpIntegrator < approx.odeint.ExrkIntegrator
    % EMPINTEGRATOR Explicit midpoint integrator.
    %
    %   EmpIntegrator implements the explicit midpoint rule, a
    %   second-order explicit time integration scheme for solving ordinary
    %   differential equations of the form:
    %
    %   \f[
    %     \frac{\mathrm{d}(Mu)}{\mathrm{d}t} = Lu + F(u) + S,
    %   \f]
    %
    %   where \f$M\f$ is a mass operator, \f$L\f$ is a linear operator,
    %   \f$F(u)\f$ is a nonlinear operator, and \f$S\f$ is a source term.
    %   Any of these terms can be empty. For mass operator, it means it
    %   reduces to an identity. For linear, nonlinear and source terms,
    %   it means they vanish from the equation.
    %
    %   The explicit midpoint rule is given by:
    %
    %   \f[
    %     u_{n+1} = u_n + h \cdot f\left(t_n + \frac{h}{2}, u_n + \frac{h}{2} f(t_n, u_n)\right)
    %   \f]
    %
    %   This method is second-order accurate and has good stability
    %   properties for non-stiff problems. It requires only explicit
    %   evaluations and no nonlinear system solving, making it much
    %   faster than implicit methods.
    %
    % See also:
    %   approx.odeint.ExrkIntegrator, approx.odeint.HeunIntegrator,
    %   approx.odeint.ImpIntegrator

    properties (Constant)
        Order = 2 % Accuracy order of the explicit midpoint method
    end

    methods
        function obj = EmpIntegrator(final)
            % EMPINTEGRATOR Constructor for EmpIntegrator.
            %
            %   obj = EmpIntegrator(final) creates an explicit
            %   midpoint integrator with the specified final time.

            arguments
                final {mustBeNumeric, mustBeNonnegative}
            end

            obj@approx.odeint.ExrkIntegrator(2, final);

            %< Explicit midpoint method coefficients
            obj.c = [0, 1/2];

            %< Butcher tableau coefficient matrix
            obj.A = [0, 0; 1/2, 0];

            %< Stage weights for final update
            obj.b = [0, 1];
        end
    end
end