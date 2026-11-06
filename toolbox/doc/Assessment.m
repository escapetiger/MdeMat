function Assessment
%ASSESSMENT Explicit sampling, physical norms and convergence observations.
%   plan = mdemat.common.fem.samplingPlan(space,[4,4]); % points on each axis
%   values = mdemat.common.fem.sampleField(space,coefficients,plan);
%   reference = exact(plan.Coordinates);
%   error = mdemat.tools.analyze.weightedError(values,reference,plan.Weights);
%
%   Values/reference are components-by-points; weights include abs(det J).
%   ComponentL1/ComponentL2/ComponentLInf retain each spatial component norm.
%   L1 sums component L1 norms; L2 combines component L2 norms in Euclidean
%   norm; LInf takes the maximum component supremum. Legacy sequence instead
%   applies Euclidean norm after every spatial norm, including L1 and LInf.
%   Thus use norm(error.ComponentL1,2) to request that particular mixed norm.
%
%   nestedSamplingPlan(fine,coarse,np) supplies common physical points, weights
%   and FineDofIndices for integer subdivisions. A difference between two levels
%   is not automatically a Richardson error estimate; division by r^p-1 needs
%   explicit order/asymptotic assumptions. Coarse modal zero padding is valid
%   only for proven compatible bases/scales/normalizations and component IDs.
%   Kinetic physical reconstruction belongs to the plasma/radiation model.
%   Plasma physical F assessment uses positive Gauss-Legendre velocity weights
%   on Config.vBBox and includes each AP State's sqrtRhoInf. This is a finite
%   velocity-domain dx dv norm, not a full-line or coefficient/Maxwellian norm.
%
%   summarizeConvergence accepts ordinary observation structs with Scale,
%   Field, ComponentIds, SpatialNorm, ComponentNorm, Comparison, Reference and
%   Error. Each identity is a separate series; no solver, callback or plot is run.
%   ComponentIds='all' is only a sequence identity and requires actual
%   ResolvedComponentIds=1:N per level; summaries retain every resolved cutoff.
%   Fixed numeric selections keep their order as part of the curve identity.
%   convergenceRates uses stable log ratios, retaining zero/negative-rate cases.
%
%   See also mdemat.tools.analyze.weightedError,
%   mdemat.tools.analyze.convergenceRates, mdemat.tools.analyze.summarizeConvergence,
%   mdemat.common.fem.samplingPlan, mdemat.common.fem.nestedSamplingPlan.
%   mdemat.application.plasma.assessmentSamples.
end
