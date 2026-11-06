function [flux, report] = ganapolReference(points, options)
%GANAPOLREFERENCE Callable independent infinite-medium transport reference.
%   Optical length/total cross section are nondimensionalized to one.
%   Phase(l+1) multiplies (2*l+1)*P_l(mu)*P_l(mu')/2. Phase(1)=1.
%   The finite Legendre system and Fourier inversion are independent of DG.
%   Origin singularities of point/line and isotropic plane delta sources are
%   returned as Inf. The plane uncollided term E1(abs(x))/2 is logarithmic;
%   finite angular quadrature must not regularize its value at x=0.
arguments
    points double {mustBeReal, mustBeFinite}
    options.Geometry (1, 1) string {mustBeMember(options.Geometry, ["plane", "sphere", "cylinder"])} = "plane"
    options.Source (1, 1) string {mustBeMember(options.Source, ["isotropic", "beam", "step"])} = "isotropic"
    options.ScatteringRatio (1, 1) double {mustBeNonnegative, mustBeLessThan(options.ScatteringRatio, 1)} = 0
    options.Phase (1, :) double {mustBeReal, mustBeFinite} = [1, 0.5, 0.3, 0.2]
    options.BeamDirection (1, 1) double {mustBeReal, mustBeFinite} = 1
    options.AngularPoints (1, 1) double {mustBeInteger, mustBePositive} = 30
    options.FourierOrder (1, 1) double {mustBeInteger, mustBeNonnegative} = 10
    options.EulerParameter (1, 1) double {mustBeFinite, mustBeNonnegative} = 0.45
    options.ProjectionPoints (1, 1) double {mustBeInteger, mustBePositive} = 32
end
if isempty(options.Phase) || options.Phase(1) ~= 1 || ...
        (options.Source == "beam" && (abs(options.BeamDirection)>1 || options.BeamDirection==0)) || ...
        (options.Geometry ~= "plane" && options.Source ~= "isotropic") || ...
        (options.Geometry ~= "plane" && any(points(:)<0))
    error('mdemat:ganapol:Problem', 'Invalid phase, beam direction or source/geometry combination.');
end
oldPath = path;
pathGuard = onCleanup(@() path(oldPath));
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));
[mu, w] = gausslegd(options.AngularPoints, -1, 1, false);
degrees = 0:numel(options.Phase)-1;
P = legdval(degrees, mu).';
collision = diag(options.ScatteringRatio*(2*degrees+1).*options.Phase);
started = tic;
flux = arrayfun(@evaluate, points);
report = struct('Status', 'COMPLETED', 'Assessment', 'UNVERIFIED', ...
    'Geometry', options.Geometry, 'Source', options.Source, ...
    'ScatteringRatio', options.ScatteringRatio, 'AngularPoints', options.AngularPoints, ...
    'FourierOrder', options.FourierOrder, 'ProjectionPoints', options.ProjectionPoints, ...
    'Duration', toc(started), 'CollisionlessAnalyticBranch', options.ScatteringRatio==0);

    function value = fourier(k, source)
        denominator = 1+1i*k*mu;
        T = (P.*(w./denominator).')*P.'/2;
        if source == "beam"
            Q = legdval(degrees, options.BeamDirection).'/(1+1i*k*options.BeamDirection);
        else
            Q = P*(w./denominator)/2;
        end
        moments = (eye(numel(degrees))-T*collision)\Q;
        value = moments(1);
    end
    function value = sphere(radius)
        if radius == 0
            value = Inf;
        elseif options.ScatteringRatio == 0
            value = exp(-radius)/(4*pi*radius^2);
        else
            transformed = @(k) fourier(k, "isotropic")*k/(2*pi*1i*radius);
            value = cifft(transformed, radius, options.FourierOrder, options.EulerParameter);
        end
    end
    function value = evaluate(x)
        if options.Geometry == "sphere"
            value = sphere(x);
        elseif options.Geometry == "cylinder"
            if x == 0
                value = Inf;
            else
                % Project over the entire missing coordinate, including both
                % halves. The spherical inverse uses its own radius hypot(x,z).
                value = 2*gaussint(@(t) sphere(hypot(x, t/(1-t)))/(1-t)^2, ...
                    0, 1, options.ProjectionPoints);
            end
        elseif options.Source == "isotropic" && x == 0
            value = Inf;
        elseif options.ScatteringRatio == 0
            if options.Source == "isotropic"
                value = expint(abs(x))/2;
            elseif options.Source == "beam"
                s = x/options.BeamDirection;
                value = double(s>0)*exp(-max(s, 0))/abs(options.BeamDirection);
                if x==0, value=1/(2*abs(options.BeamDirection)); end
            else
                a = abs(x);
                if a == 0
                    value = 0.5;
                else
                    tail = (exp(-a)-a*expint(a))/2;
                    value = tail;
                    if x>0, value=1-tail; end
                end
            end
        elseif options.Source == "step"
            transformed = @(k) fourier(k, "isotropic")/(1i*k);
            value = cifft(transformed, x, options.FourierOrder, options.EulerParameter) ...
                +1/(2*(1-options.ScatteringRatio));
        else
            value = cifft(@(k) fourier(k, options.Source), x, ...
                options.FourierOrder, options.EulerParameter);
        end
    end
end
