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
    %   It achieves second-order accuracy and A-stability.
    %
    % See also:
    %   approx.odeint.BdfIntegrator, approx.odeint.Bdf3Integrator
    
    properties (Constant)
        Order = 2 % Accuracy order
    end

    methods
        function obj = Bdf2Integrator(final)
            % BDF2INTEGRATOR Create a BDF2 integrator.
            %
            %   obj = Bdf2Integrator(final) creates a BDF2 integrator with
            %   the specified final time.

            arguments
                final {mustBeNumeric, mustBeNonnegative}
            end
            
            obj@approx.odeint.BdfIntegrator(2, final);
            obj.alpha = [-4/3, 1/3];
            obj.beta = 2/3;
        end

        function obj = setCoefficients(obj)
            % SETCOEFFICIENTS Set BDF2 coefficients for variable step
            % sizes.
            %
            %   obj = setCoefficients(obj) computes the alpha and beta
            %   coefficients based on the current and previous step sizes.

            arguments
                obj approx.odeint.Bdf2Integrator
            end

            if obj.Timeline.Count == 1
                obj.alpha = -1;
                obj.beta = 1;
                return;
            end
            
            h1 = obj.Timeline.StepSizeQueue(1);
            h2 = obj.Timeline.StepSizeQueue(2);
            obj.alpha = [-(h1 + h2)^2 / (h2 * (2 * h1 + h2)), h1^2 / (h2 * (2 * h1 + h2))];
            obj.beta = (h1 + h2) / (2 * h1 + h2);            
        end
    end
end