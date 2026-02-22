classdef ImpIntegrator < approx.odeint.DirkIntegrator
    % IMPINTEGRATOR Implicit midpoint integrator.
    %
    %   ImpIntegrator implements the implicit midpoint rule, a
    %   second-order implicit time integration scheme for solving ordinary
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
    %   The implicit midpoint rule is given by:
    %
    %   \f[
    %     u_{n+1} = u_n + h \cdot f\left(t_n + \frac{h}{2}, \frac{u_n + u_{n+1}}{2}\right)
    %   \f]
    %
    %   This method is A-stable, symplectic, and has second-order accuracy.
    %   It is particularly suitable for Hamiltonian systems and problems
    %   requiring energy conservation. The method requires solving a
    %   nonlinear system at each time step.
    %
    % See also:
    %   approx.odeint.OdeIntegrator, approx.odeint.BeIntegrator,
    %   approx.odeint.BdfIntegrator

    properties (Constant)
        Order = 2 % Accuracy order of the implicit midpoint method
    end

    methods
        function obj = ImpIntegrator(final)
            % IMPINTEGRATOR Constructor for ImpIntegrator.
            %
            %   obj = ImpIntegrator(final) creates an implicit
            %   midpoint integrator with the specified final time.

            arguments
                final {mustBeNumeric, mustBeNonnegative}
            end

            obj@approx.odeint.DirkIntegrator(1, final);
            obj.c = 0.5;
            obj.A = 0.5;
            obj.b = 1;
        end
    end
end