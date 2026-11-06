# Sampled assessment

For finite samples $e_{cq}=u_{cq}-u^\star_{cq}$ and positive physical quadrature
weights $w_q=|\det J_K|\widehat w_q$, `weightedError` computes

\[
a_{1,c}=\sum_q w_q|e_{cq}|,\qquad
a_{2,c}=\left(\sum_q w_q|e_{cq}|^2\right)^{1/2},\qquad
a_{\infty,c}=\max_{q:w_q>0}|e_{cq}|.
\]

These are `ComponentL1`, `ComponentL2`, and `ComponentLInf`. Global `L1` is
$\sum_c a_{1,c}$, global `L2` is $\|a_2\|_2$, and global `LInf` is
$\max_c a_{\infty,c}$. The retained study `sequence` convention is instead
$\|a_p\|_2$ for each spatial norm $p$. Disjoint component supports distinguish
these mixed norms. No volume normalization is implicit.

A fixed ordered component selection, such as `[2,4]`, belongs to the curve
identity; `[4,2]` remains a distinct selection even when its sequence norm agrees.
For mode refinement, `ComponentIds='all'` with `ComponentNorm='sequence'` denotes
the complete current expansion, rather than a fixed cutoff. At level $j$ its
error is

\[
 E_{p,j}=\left(\sum_{c=1}^{N_j} a_{p,c,j}^{\,2}\right)^{1/2}.
\]

Each observation records `ResolvedComponentIds=1:N_j`; the summary retains these
IDs in a cell for every level while grouping the explicit `all` identity into one
curve. The cutoff is thus visible even when it changes between levels. This is a
sequence of truncated coefficient diagnostics, without an implicit claim that
it equals the physical velocity-space norm. Interlevel padding still requires
the common nested basis and normalization described below.

Real/imaginary differences and weighted sums use scaled logarithms when ordinary
arithmetic overflows; representable ordinary arithmetic retains its direct result.
Each output norm can still be mathematically outside double range. Relative L2
is declared defined only when both norms and their ratio are finite and the
reference exceeds the explicit floor. Zero-weight points do not enter arithmetic.

`samplingPlan` and `nestedSamplingPlan` contain numbers only. The latter aligns
fine local DOFs and coarse subcell reference points, preserving DG traces and
actual physical Jacobians. Interlevel difference $\|u_h-u_H\|$ is distinct from
the asymptotic fine-error estimate $\|u_h-u_H\|/(r^p-1)$, which requires an explicit
order and refinement ratio. Velocity-mode padding requires a common nested,
identically normalized basis and scale; general kinetic comparisons reconstruct
on common physical velocity samples or apply the model's specified Gram measure.

Observed rates are
$p_i=\log(E_i/E_{i+1})/\log(h_i/h_{i+1})$. Nearby positive inputs use `log1p`
of the relative difference; widely separated inputs use separate logarithms.
Zero errors give an undefined logarithmic rate; negative rates remain visible.
The tool does not infer whether the measured error is spatial, temporal or modal.

## Physical plasma distributions

`application.plasma.assessmentSamples` measures scalar `F` using physical
Lebesgue measure $dx\,dv$ on the explicit finite interval
$v\in[\mathrm{vBBox}_1,\mathrm{vBBox}_2]$. It is a truncated velocity-domain
assessment, distinct from both a full-line norm and a Maxwellian-weighted norm
of Hermite coefficients. The returned reference identity records these bounds.

Each State reconstructs using its own basis and temperature:

\[
f_h^{H}(x,v)=\sum_k D_k(x)\phi_k^T(v),\qquad
f_h^{AP}(x,v)=\sqrt{\rho_\infty(x)}\sum_k D_k(x)\phi_k^T(v).
\]

Positive Gauss-Legendre velocity weights combine with physical spatial weights
as `kron(velocityWeights,spatialWeights)`. Values are flattened with the spatial
point index fastest. Thus the physical L2 error is

\[
 E_2=\left(\sum_{j,q}w_j^v w_q^x
       |f_h(x_q,v_j)-f_\star(x_q,v_j)|^2\right)^{1/2}.
\]

The default velocity count is $\max(8,2\max(N_v^{fine},N_v^{coarse}))$; it is a
quadrature choice, not an exactness claim for Gaussian products. Independent
Gaussian-product integrals test the numerical measure. Interlevel samples use
the same physical points while each State evaluates its own basis and AP scale;
no coefficient padding or equality of temperatures is assumed. `D` remains an
explicitly selected coefficient diagnostic with its separately declared meaning.
