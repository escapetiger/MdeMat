# Independent infinite-medium transport reference

The developer reference `benchmarks/ganapol/ganapolReference.m` is independent
of DG assembly and is excluded from the installed toolbox. Optical length and
total cross section are scaled to one. The plane equation is

$$
\mu\partial_x\psi+\psi
=\frac c2\sum_{\ell=0}^{L}(2\ell+1)f_\ell P_\ell(\mu)
 \int_{-1}^1P_\ell(\mu')\psi(x,\mu')\,d\mu'+q(x,\mu).
$$

With $\widehat g(k)=\int g(x)e^{-ikx}dx$ and scalar moments
$\widehat\phi_m=\int P_m\widehat\psi\,d\mu$, define

$$
T_{m\ell}(k)=\frac12\int_{-1}^1
 \frac{P_m(\mu)P_\ell(\mu)}{1+ik\mu}\,d\mu,\quad
[I-T\operatorname{diag}(c(2\ell+1)f_\ell)]\widehat\phi=Q.
$$

The nested `fourier` function implements these terms directly. Isotropic plane
source uses $Q_m=\tfrac12\int P_m/(1+ik\mu)$; a beam at $\mu_0\ne0$ uses
$Q_m=P_m(\mu_0)/(1+ik\mu_0)$. The original drivers retain their phase/scattering
sweeps, including the separate coefficient sweep that does not impose the
normalized physical phase constraint of this new reference interface.

`utils/cifft` computes $\pi^{-1}\Re\int_0^\infty\widehat g(k)e^{ikx}dk$ by
oscillatory segments/Euler-Knopp acceleration. At $x=0$, it integrates the whole
non-oscillatory integral once. A plane step source has Fourier factor $1/(ik)$
plus its analytically treated delta contribution $1/[2(1-c)]$; no Symbolic
`dirac` evaluation is needed in the numerical integrand.

The three-dimensional point reference is
$\phi_3(r)=\operatorname{cifft}[k\widehat\phi_0(k)/(2\pi i r)](r)$.
The cylindrical line reference is its projection

$$
\phi_2(r)=2\int_0^\infty\phi_3(\sqrt{r^2+z^2})\,dz.
$$

The factor two and the radius inside the spherical inverse are essential.
The source origin is singular and returns `Inf`; it is not assigned a finite
value by a truncation artifact. For $c=0$, independent analytic limits are
$\phi_1(x)=E_1(|x|)/2$ and $\phi_3(r)=e^{-r}/(4\pi r^2)$.
`TestGanapolReferences` checks these limits, the cylindrical projection by a
different angular substitution, Fourier identities and refinement. A completed
reference call reports `UNVERIFIED` assessment until an experiment actually
compares/assesses it; finite truncation does not certify its own accuracy.
