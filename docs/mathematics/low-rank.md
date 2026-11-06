# Matrix-factor dynamics and the slab transport example

The optional matrix representation is $F=USV^*$, with $U\in\mathbb C^{n_x\times r}$,
$S\in\mathbb C^{r\times r}$, $V\in\mathbb C^{n_v\times r}$. Arrays outside this
consumer continue to use ordinary dense/sparse MATLAB storage. `Rank` denotes
factor width. Dense expansion is explicit (`materialize`); `StoredEntries`
counts factor payload, not metric/operator/workspace storage.

## Physical metric, algebra and operator action

For fixed Hermitian positive definite $M_x,M_v$,

$$
\|F\|_M^2=\operatorname{tr}(F^*M_xFM_v),\qquad
U^*M_xU=I,\quad V^*M_vV=I.
$$

`common.lowrank.weightedQr` uses $M=C^*C$ and an ordinary QR of $CU$;
`compress` takes the SVD of the small core $R_xSR_v^*$. It selects $k$ with
$\sqrt{\sum_{j>k}\sigma_j^2}\leq\max(\varepsilon_a,\varepsilon_r\|F\|_M)$.
A binding rank cap may prevent that inequality: the report explicitly records
`ToleranceMet=false` and the discarded norm. No inverse of $S$ is formed.

For $\mathcal A(F)=\sum_\ell A_\ell F B_\ell^*$,
`common.operator.applySeparated` forms factors $[A_\ell U]$,
$\operatorname{blockdiag}(S,\ldots,S)$ and $[B_\ell V]$, then recompresses.
The algorithm does not build a global Kronecker matrix or full phase-space field.
Independent tests may build those small dense references explicitly.

## Fixed-rank time evolution

`common.ode.evolveFixedRank` integrates the projected equation

$$
\dot F=P_FG(t,F),\quad
P_FG=UU^*M_xG+GM_vVV^*-UU^*M_xGM_vVV^*.
$$

One K-S-L step freezes the indicated orthonormal factors and integrates

$$
\dot K=G(t,KV_0^*)M_vV_0,\quad K_1=U_1\widehat S,
$$
$$
\dot S=-U_1^*M_xG(t,U_1SV_0^*)M_vV_0,\quad S(t_n)=\widehat S,
$$
$$
\dot L=G(t,U_1L^*)^*M_xU_1,\quad L(t_n)=V_0\widetilde S^*,
\quad L_1=V_1S_1^*.
$$

Weighted QR supplies each reorthogonalization. Each subproblem uses RK4 over
the original step interval. The composition has **first-order splitting**;
fourth-order substeps do not make the projected evolution fourth order.
The method avoids division by small singular values. It does not claim general
mass, positivity, energy or asymptotic-preserving invariants.

## Rank adaptation and restart

`common.ode.evolveAdaptive` is an explicit Euler step $F+hG$, followed by
weighted compression. This is a distinct first-order retraction method.
`RepresentationChangeBounds` records each compression norm plus any subsequent
moment correction; it is not a bound on the full PDE solution error.

For $a=U^*M_x\mathbf1$, $b=V^*M_v\mathbf1$, the weighted integral is $a^*Sb$.
When `PreserveMass=true`, the small-core correction

$$
\Delta S=\frac{(m_{\rm target}-a^*Sb)ab^*}{\|a\|^2\|b\|^2}
$$

restores that one scalar moment. A subspace that cannot represent a nonzero
required correction causes an explicit error. This correction does not restore
pointwise positivity. `MatrixFactors.snapshot/fromSnapshot` stores a versioned
plain numeric record. Restart occurs at completed Euler steps and reconstructs
fresh factors/metrics rather than loading a handle graph.

## Selected kinetic problem and equation-to-code map

`application.radiation.solveLowRankTransport` solves periodic slab transport
with nondimensional total speed/scattering inputs:

$$
\partial_t f+v\partial_x f=\nu(\rho-f),\qquad
\rho_i=\frac{\sum_jw_jF_{ij}}{\sum_jw_j},\quad w_j>0.
$$

For cell width $h_x$, it uses $M_x=h_xI$, $M_v=\operatorname{diag}(w)$ and

$$
\dot F=-D_-F\operatorname{diag}(v^+)
       -D_+F\operatorname{diag}(v^-)
       +\nu F\left(\frac{w\mathbf1^T}{\mathbf1^Tw}-I\right),
$$

where $(D_-u)_i=(u_i-u_{i-1})/h_x$ and
$(D_+u)_i=(u_{i+1}-u_i)/h_x$, with periodic indexing. The `left`/`right` cell
arrays in that function correspond to those three terms in order. The method
reports spatial order one. The Euler monotonicity bound
$\Delta t(\max_j|v_j|/h_x+\nu)\leq1$ applies before compression; a compressed
solution has no general positivity guarantee and this is not an AP scheme.

`tests/numerical/TestLowRankDynamics.m` separates small-singular-value evolution,
comparison with independently projected `ode45`, rank growth/drop, numeric
restart and transport comparison with a dense matrix exponential. Space, time,
quadrature, rank truncation and operator approximation remain distinct errors.
