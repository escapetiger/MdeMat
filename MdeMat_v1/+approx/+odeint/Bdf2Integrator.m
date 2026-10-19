classdef Bdf2Integrator < approx.odeint.BdfIntegrator
    % BDF2INTEGRATOR Second-order BDF integrator.
    %
    %   Bdf2Integrator implements two-step Backward Differentiation
    %   Formula (BDF) methods for solving ordinary differential equations
    %   of the form:
    %
    %   \f[
    %     \frac{\mathrm{d}(Mu)}{\mathrm{d}t} = Lu + S,
    %   \f]
    %
    %   where \f$M\f$ is a mass operator, \f$L\f$ is a linear operator, and
    %   \f$S\f$ is a source term.
    %
    %   Bdf2Integrator uses the Backward Euler method for the first step.
    %   It acheives second-order accuracy and A-stability.
    %
    % Examples:
    %   integrator = approx.odeint.Bdf2Integrator(1.0);
    %   U_new = integrator.step(L, N);
    %   integrator.update(U_new);
    %
    % See also:
    %   approx.odeint.BdfIntegrator, approx.odeint.Bdf3Integrator
    
    properties (Constant)
        ORDER = 2 % Accuracy order
    end

    methods
        function obj = Bdf2Integrator(final)
            % BDF2INTEGRATOR Create a BDF2 integrator.
            %
            % Syntax:
            %   obj = Bdf2Integrator(final)
            %
            % Inputs:
            %   final - Final time for integration
            %
            % Outputs:
            %   obj - The constructed Bdf2Integrator object
            
            obj@approx.odeint.BdfIntegrator(2, final);
            obj.alpha = [-4/3, 1/3];
            obj.beta = 2/3;
        end

        function setCoefficients(obj)
            % SETCOEFFICIENTS Set BDF2 coefficients for variable step
            % sizes.
            %
            % Computes the alpha and beta coefficients based on the current
            % and previous step sizes. For uniform step sizes, uses the
            % standard BDF2 coefficients. For the first step, reduces to
            % Backward Euler coefficients.
            %
            % Syntax:
            %   setCoefficients(obj)
            %
            % Inputs:
            %   obj - The Bdf2Integrator object

            if obj.timeline.count == 1
                obj.alpha = -1;
                obj.beta = 1;
                return;
            end
            
            h1 = obj.timeline.h(1);
            h2 = obj.timeline.h(2);
            obj.alpha = [-(h1 + h2)^2 / (h2 * (2 * h1 + h2)), h1^2 / (h2 * (2 * h1 + h2))];
            obj.beta = (h1 + h2) / (2 * h1 + h2);
        end
    end
end