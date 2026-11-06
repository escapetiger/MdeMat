function [field, report] = evolveFixedRank(initial, rhs, interval, step)
%EVOLVEFIXEDRANK First-order K-S-L projector splitting for matrix dynamics.
%   RHS(t,F) returns a MatrixFactors coefficient field. Weighted projections
%   use the fixed SPD left/right metrics. Subproblems use RK4 without core
%   inversion. No complete phase-space field is formed by this implementation.
%   Global splitting order is one; inner RK4 does not imply fourth-order dynamics.
arguments
    initial mdemat.data.tensor.MatrixFactors
    rhs function_handle
    interval (1, 2) double {mustBeFinite}
    step (1, 1) double {mustBePositive, mustBeFinite}
end
if interval(2) < interval(1) || initial.Rank == 0
    error('mdemat:ode:FactorInitialState', 'Require an ordered interval and positive representation rank.');
end
r = initial.Rank;
[field, ~] = mdemat.common.lowrank.compress(initial.U, initial.S, initial.V, ...
    LeftMetric=initial.LeftMetric, RightMetric=initial.RightMetric, ...
    Axes=initial.Axes, RelativeTolerance=0, MinRank=r, MaxRank=r);
t = interval(1); count = 0;
while t < interval(2)
    h = min(step, interval(2) - t);
    V = field.V;
    K = rk4(@(time, value) kAction(time, value, V), t, field.U*field.S, h);
    [U, S] = mdemat.common.lowrank.weightedQr(K, field.LeftMetric);
    S = rk4(@(time, value) sAction(time, value, U, V), t, S, h);
    L = rk4(@(time, value) lAction(time, value, U), t, V*S', h);
    [V, R] = mdemat.common.lowrank.weightedQr(L, field.RightMetric);
    field = field.withFactors(U, R', V);
    t = t + h; count = count + 1;
end
report = struct('Method', 'KSL-RK4-substeps', 'SplittingOrder', 1, ...
    'Steps', count, 'Time', t, 'Rank', field.Rank, 'StoredEntries', field.StoredEntries);

    function value = kAction(time, K, frozenV)
        G = checkedRhs(time, field.withFactors(K, eye(r), frozenV));
        value = G.U * G.S * (G.V' * action(field.RightMetric, frozenV));
    end
    function value = sAction(time, S, frozenU, frozenV)
        G = checkedRhs(time, field.withFactors(frozenU, S, frozenV));
        value = -(frozenU' * action(field.LeftMetric, G.U)) * G.S ...
            * (G.V' * action(field.RightMetric, frozenV));
    end
    function value = lAction(time, L, frozenU)
        G = checkedRhs(time, field.withFactors(frozenU, eye(r), L));
        value = G.V * G.S' * (G.U' * action(field.LeftMetric, frozenU));
    end
    function G = checkedRhs(time, F)
        G = rhs(time, F);
        if ~isa(G, 'mdemat.data.tensor.MatrixFactors') || ~isequal(G.Shape, F.Shape) || ...
                ~isequal(G.Axes, F.Axes) || ~isequal(G.LeftMetric, F.LeftMetric) || ...
                ~isequal(G.RightMetric, F.RightMetric)
            error('mdemat:ode:FactorRhs', 'RHS must return factors on the original axes.');
        end
    end
end

function next = rk4(fun, t, value, h)
k1 = fun(t, value);
k2 = fun(t+h/2, value+h*k1/2);
k3 = fun(t+h/2, value+h*k2/2);
k4 = fun(t+h, value+h*k3);
next = value + h*(k1+2*k2+2*k3+k4)/6;
end

function value = action(M, X)
if isempty(M), value = X; else, value = M*X; end
end
