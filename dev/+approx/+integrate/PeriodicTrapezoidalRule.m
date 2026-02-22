classdef PeriodicTrapezoidalRule < approx.integrate.OrthotopeRule
    % PERIODICTRAPEZOIDALRULE Periodic trapezoidal numerical integration
    % rule.
    %
    %   PeriodicTrapezoidalRule implements numerical integration using the
    %   trapezoidal method on periodic domains. This method is particularly
    %   efficient for integrating periodic functions, providing exponential
    %   convergence for smooth periodic integrands.
    %
    %   The periodic trapezoidal rule assumes the integrand is periodic
    %   over the integration interval. It excludes the right endpoint
    %   to avoid double-counting due to periodicity, using equally spaced
    %   nodes over the half-open interval \f$[a, b)\f$.
    %
    % See also:
    %   approx.integrate.ClosedNewtonCotesRule, 
    %   approx.integrate.GaussLegendreRule,
    %   approx.integrate.GaussLobattoRule

    methods (Access = protected)
        function [x, w] = generate1d(~, np, options)
            % GENERATE1D Generate 1D periodic trapezoidal nodes and weights.
            %
            %   [x, w] = generate1d(obj, np) generates @a np equally spaced
            %   nodes @a x and uniform weights @a w for the interval \f$[0,
            %   1)\f$, excluding the right endpoint.
            %
            %   [x, w] = generate1d(obj, np, lower=a, upper=b) generates
            %   nodes and weights for the interval \f$[a, b)\f$.
            
            arguments
                ~
                np(1, 1) {mustBePositive, mustBeInteger}
                options.lower(1, 1) {mustBeReal} = 0
                options.upper(1, 1) {mustBeReal} = 1
            end
            
            a = options.lower;
            b = options.upper;
            
            core.except.assert(a <= b, 'InvalidInput', ...
                'a <= b is required.');
            
            if a == b
                [x, w] = deal(a, 1);
                return;
            end
            
            x = linspace(0, 1-1/np, np);
            w = ones(1, np) / np;
            x = a + (b - a) * x;
            w = w * (b - a);
        end
    end
end