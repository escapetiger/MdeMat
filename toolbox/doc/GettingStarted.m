function result = GettingStarted()
%GETTINGSTARTED Use MdeMat without plotting or implicit file output.
%   The installed toolbox adds one path root for the mdemat namespace.
%   Explicit guide entry after installation:
%     addpath(fullfile(mdemat.tools.environment.toolboxRoot(), 'doc'))
%     result = GettingStarted();
%   Remove that optional documentation path afterwards if it is no longer used.
%   Core numerical methods use base MATLAB. Symbolic inputs and source
%   regeneration explicitly require the optional Symbolic Math Toolbox.
%
%   This example solves -u''=2 on (0,1), u(0)=u(1)=0. The exact continuous
%   solution is u=x(1-x). The centered derivative is exact for this quadratic.
%   Rows identify interior physical nodes; L acts on their column values.
%   To persist a run, construct createRunRecord(result, executionReport) and
%   pass it to writeRecord with an explicit destination. Missing assessment is
%   UNVERIFIED; UTC/seed/code version are explicit caller provenance inputs.
%
%   See also mdemat.common.fdm.secondDerivative,
%   mdemat.tools.analyze.weightedError, mdemat.tools.io.createRunRecord,
%   mdemat.tools.environment.describe, mdemat.tools.io.writeRecord.
[L, x, boundaryLift, h] = mdemat.common.fdm.secondDerivative(15, [0, 1]);
u = (-L) \ (2*ones(size(x)) + boundaryLift*[0; 0]);
exact = x .* (1-x);
assessment = mdemat.tools.analyze.weightedError(u.', exact.', h*ones(size(x)));
result = struct('Coordinates', x, 'Values', u, 'Reference', exact, ...
    'Assessment', assessment, 'Equation', '-u''''=2, u(0)=u(1)=0');
end
