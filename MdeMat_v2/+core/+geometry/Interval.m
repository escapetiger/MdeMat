classdef Interval < core.geometry.Orthotope
    % INTERVAL One-dimensional interval.
    %
    %   Interval represents a one-dimensional segment of the real line.
    %   It supports open, closed, and half-open intervals with configurable
    %   boundary inclusion properties.
    %
    % Examples:
    %   % Create a closed interval [0,1]
    %   interval1 = core.geometry.Interval([0, 1]);
    %
    %   % Create a half-open interval [0,1)
    %   interval2 = core.geometry.Interval([0, 1], 'IncludeUpper', false);
    %
    %   % Create an open interval (0,1)
    %   interval3 = core.geometry.Interval([0, 1], 'IncludeLower', false, 'IncludeUpper', false);
    %
    % See also:
    %   core.geometry.Orthotope

    properties
        includeLower % Whether the interval includes the lower bound
        includeUpper % Whether the interval includes the upper bound
    end

    methods
        function obj = Interval(bbox, varargin)
            % INTERVAL Constructor for Interval.
            %
            %   obj = Interval(bbox) creates a closed one-dimensional
            %   interval on the real line using the bounding box format
            %   [lower, upper]. The interval includes both endpoints by
            %   default.
            %
            %   obj = Interval(bbox, 'IncludeLower', includeLowerFlag)
            %   creates an interval with configurable lower boundary
            %   inclusion.
            %
            %   obj = Interval(bbox, 'IncludeUpper', includeUpperFlag)
            %   creates an interval with configurable upper boundary
            %   inclusion.
            %
            %   obj = Interval(bbox, 'IncludeLower', includeLowerFlag,
            %   'IncludeUpper', includeUpperFlag) creates an interval with
            %   both boundary inclusion properties specified.
            %
            % Inputs:
            %   bbox - Bounding box vector [lower, upper] defining the interval
            %   varargin - Input arguments
            %<   'IncludeLower' - Whether includes lower bound (logical, default: true)
            %<   'IncludeUpper' - Whether includes upper bound (logical, default: true)
            %
            % Outputs:
            %   obj - Constructed Interval object

            core.except.assert(nargin >= 1, 'InvalidInput', ...
                'Bounding box must be specified.');

            core.except.assert(isvector(bbox), 'InvalidInput', ...
                'Bounding box must be a vector.');

            core.except.assert(length(bbox) == 2, 'InvalidInput', ...
                'Interval bounding box must have exactly 2 elements [lower, upper].');

            bbox = bbox(:)';

            core.except.assert(bbox(1) <= bbox(2), 'InvalidInput', ...
                'Lower bound must be less than or equal to upper bound.');

            obj@core.geometry.Orthotope(bbox);

            obj.includeLower = true;
            obj.includeUpper = true;

            if nargin > 1
                p = inputParser;
                p.addParameter('IncludeLower', true, @(x) islogical(x) && isscalar(x));
                p.addParameter('IncludeUpper', true, @(x) islogical(x) && isscalar(x));
                p.parse(varargin{:});

                obj.includeLower = p.Results.IncludeLower;
                obj.includeUpper = p.Results.IncludeUpper;
            end
        end

        function TF = isInside(obj, X)
            % ISINSIDE Tests if points are strictly inside the interval.
            %
            %   TF = isInside(obj, X) determines if points are strictly
            %   inside the interval, accounting for whether the interval
            %   includes its endpoints. Points on included boundaries are
            %   not considered inside. The test considers infinite bounds
            %   appropriately.
            %
            % Inputs:
            %   obj - The Interval object
            %   X - A vector of points to test
            %
            % Outputs:
            %   TF - A logical vector of the same size as X where true
            %        indicates the corresponding point is inside

            n = size(X);
            core.except.assert(n(1) == 1, 'DimensionMismatch', ...
                'Point dimension must match geometry dimension.');

            if isinf(obj.lower)
                if obj.lower < 0
                    lowerCondition = true(size(X));
                else
                    lowerCondition = false(size(X));
                end
            else
                lowerCondition = X > obj.lower;
            end

            if isinf(obj.upper)
                if obj.upper > 0
                    upperCondition = true(size(X));
                else
                    upperCondition = false(size(X));
                end
            else
                upperCondition = X < obj.upper;
            end

            % Combined condition
            TF = lowerCondition & upperCondition;
        end

        function TF = isOnBoundary(obj, X)
            % ISONBOUNDARY Tests if points are on the interval boundary.
            %
            %   TF = isOnBoundary(obj, X) determines if points are exactly
            %   on the boundary of the interval. Points are on the boundary
            %   if they equal an endpoint that is included in the interval.
            %   Infinite boundaries are handled appropriately.
            %
            % Inputs:
            %   obj - The Interval object
            %   X - A vector of points to test
            %
            % Outputs:
            %   TF - A logical vector of the same size as X where true
            %        indicates the corresponding point is on the boundary

            n = size(X);
            core.except.assert(n(1) == 1, 'DimensionMismatch', ...
                'Point dimension must match geometry dimension.');

            tol = sqrt(eps);

            % Check if points are at included boundaries
            atLower = abs(X-obj.lower) < tol & obj.includeLower;

            if isinf(obj.upper)
                atUpper = false(size(X));
            else
                atUpper = abs(X-obj.upper) < tol & obj.includeUpper;
            end

            TF = atLower | atUpper;
        end

        function str = toString(obj)
            % TOSTRING Returns a string representation of the interval.
            %
            %   str = toString(obj) creates a string in mathematical
            %   notation format using appropriate bracket types to indicate
            %   boundary inclusion. Square brackets indicate included
            %   boundaries, parentheses indicate excluded boundaries.
            %
            % Inputs:
            %   obj - The Interval object
            %
            % Outputs:
            %   str - String representation of the interval in mathematical notation

            if obj.includeLower
                openBracket = '[';
            else
                openBracket = '(';
            end

            if obj.includeUpper
                closeBracket = ']';
            else
                closeBracket = ')';
            end

            str = [openBracket, num2str(obj.lower), ',', num2str(obj.upper), closeBracket];
        end
    end
end