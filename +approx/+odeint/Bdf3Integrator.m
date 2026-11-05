classdef Bdf3Integrator < approx.odeint.BdfIntegrator
    % BDF3INTEGRATOR Third-order BDF integrator.
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
    %   Bdf3Integrator uses the Backward Euler method for the first step
    %   and BDF2 method for the second step. It acheives third-order
    %   accuracy and A-stability.
    %
    % Examples:
    %   integrator = approx.odeint.Bdf3Integrator(1.0);
    %   U_new = integrator.step(L, N);
    %   integrator.update(U_new);
    %
    % Notes:
    % For the first time step, reduces to Backward Euler (first-order).
    % For the second time step, reduces to BDF2 (second-order).
    % Full third-order accuracy is achieved from the third step onward.
    %
    % See also:
    %   approx.odeint.BdfIntegrator, approx.odeint.Bdf2Integrator
        
    properties (Constant)
        ORDER = 3 % Accuracy order
    end
   
    methods
        function obj = Bdf3Integrator(final)
            % BDF3INTEGRATOR Create a BDF3 integrator.
            %
            % Syntax:
            %   obj = Bdf3Integrator(final)
            %
            % Inputs:
            %   final - Final time for integration
            %
            % Outputs:
            %   obj - The constructed Bdf3Integrator object
            
            obj@approx.odeint.BdfIntegrator(3, final);
            obj.alpha = [-18/11, 9/11, -2/11];
            obj.beta = 6/11;
        end

        function setCoefficients(obj)
            % SETCOEFFICIENTS Set BDF3 coefficients for variable step
            % sizes.
            %
            % Computes the alpha and beta coefficients based on the current
            % and previous step sizes. Automatically handles cases with
            % insufficient history by reducing to lower-order methods.
            %
            % Syntax:
            %   setCoefficients(obj)
            %
            % Inputs:
            %   obj - The Bdf3Integrator object
            %
            % Notes:
            % The coefficients are derived from the BDF3 formula:
            % (11*U_n+1 - 18*U_n + 9*U_n-1 - 2*U_n-2)/(6*dt) = f(t_n+1, U_n+1)
            % For variable step sizes, the coefficients are modified
            % according to the step size ratios using divided differences.

            if obj.timeline.count == 1
                obj.alpha = -1;
                obj.beta = 1;
                return;
            end

            h = obj.timeline.h;

            if obj.timeline.count == 2
                h1 = h(1);
                h2 = h(2);
                d0 = (h2 * (2 * h1 + h2));
                d1 = -(h1 + h2)^2;
                d2 = h1^2;
                obj.alpha = [d1 / d0, d2 / d0];
                obj.beta = (h1 + h2) / (2 * h1 + h2);
                return;
            end

            d0 = sum(1./cumsum(h(1:3)));
            d1 = -sum(h(1:2)) * sum(h(1:3)) / (h(1) * h(2) * sum(h(2:3)));
            d2 = h(1) * sum(h(1:3)) / (sum(h(1:2)) * h(2) * h(3));
            d3 = -h(1) * sum(h(1:2)) / (sum(h(1:3)) * sum(h(2:3)) * h(3));
            obj.alpha = [d1 / d0, d2 / d0, d3 / d0];
            obj.beta = 1 / (d0 * h(1));
        end
    end
end