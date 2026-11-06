classdef Bdf3Integrator < mdemat.common.ode.BdfIntegrator
    % BDF3INTEGRATOR Third-order BDF integrator.
    %
    %   Bdf3Integrator implements three-step Backward Differentiation
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
    %   Third-order SDIRK startup supplies the first two accepted states.
    %   Startup="legacy" preserves the original BE/BDF2 startup; at fixed
    %   step size its initial errors generally limit global order to two.
    %
    % See also:
    %   mdemat.common.ode.BdfIntegrator, mdemat.common.ode.Bdf2Integrator
    
    properties (Constant)
        Order = 3 % Accuracy order
    end

    properties
        Startup (1, 1) string {mustBeMember(Startup, ["sdirk3", "legacy"])} = "sdirk3"
    end
    
    methods
        function obj = Bdf3Integrator(final)
            % BDF3INTEGRATOR Create a BDF3 integrator.
            %
            %   obj = Bdf3Integrator(final) creates a BDF3 integrator with
            %   the specified final time.

            arguments
                final {mustBeNumeric, mustBeNonnegative}
            end
            
            obj@mdemat.common.ode.BdfIntegrator(3, final);
            obj.alpha = [-18/11, 9/11, -2/11];
            obj.beta = 6/11;
        end
        
        function obj = setCoefficients(obj)
            % SETCOEFFICIENTS Set BDF3 coefficients for variable step
            % sizes.
            %
            %   obj = setCoefficients(obj) computes the alpha and beta
            %   coefficients based on the current and previous step sizes.

            arguments
                obj mdemat.common.ode.Bdf3Integrator
            end
            
            if obj.Timeline.Count == 1
                obj.alpha = -1;
                obj.beta = 1;
                return;
            end
            
            h = obj.Timeline.StepSizeQueue;
            
            if obj.Timeline.Count == 2
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

        function value = step(obj, options)
            % STEP Use an order-matched start before the multistep formula.
            arguments
                obj mdemat.common.ode.Bdf3Integrator
                options.L = []
                options.F = []
                options.S = []
                options.M = []
                options.P = []
                options.A = []
                options.B = []
                options.C = []
                options.dt = []
                options.fn = []
            end
            if ~isempty(options.fn)
                value = options.fn(obj, options);
                return;
            end
            args = namedargs2cell(options);
            if obj.Startup == "legacy" || obj.Timeline.Count >= 3
                value = step@mdemat.common.ode.BdfIntegrator(obj, args{:});
                return;
            end
            if isempty(options.dt), dt = obj.Timeline.StepSize; else, dt = options.dt; end
            starter = mdemat.common.ode.Sdirk3Integrator(obj.Timeline.Final);
            starter.NonlinearTolerance = obj.NonlinearTolerance;
            starter.NonlinearMaxIterations = obj.NonlinearMaxIterations;
            starter.Timeline.Now = obj.Timeline.Now;
            starter.setTimeStep(dt=dt);
            starter.update(obj.U0{1});
            value = starter.step(args{:});
            obj.SolveReports = starter.SolveReports;
        end
    end
end
