# First finite-difference consumer

`common.fdm.secondDerivative(n,domain,Boundary=...)` builds the sparse central
approximation to **positive** $\partial_{xx}$ on a uniform one-dimensional mesh:

$$
(Lu)_i=\frac{u_{i-1}-2u_i+u_{i+1}}{h^2}.
$$

For Dirichlet boundaries, $n$ denotes interior points, $h=(b-a)/(n+1)$ and
`boundaryLift` is $B$ in

$$
u_{xx}\approx Lu+B\begin{bmatrix}g_a\\g_b\end{bmatrix}.
$$

Thus $-u_{xx}=f$ is solved by `(-L)\(f+B*[g_a;g_b])`. This sign convention is
covered by an independent quadratic solution with nonzero boundary values.
For periodic boundaries $n\geq3$ denotes distinct points (the endpoint is not
duplicated), $h=(b-a)/n$, and the two corner entries implement wraparound.
`L*ones(n,1)=0`; solving periodic Poisson additionally requires compatible
forcing and an explicit gauge, which this stencil generator does not select.

Taylor expansion gives $Lu=u_{xx}+h^2u_{xxxx}/12+O(h^4)$ for a smooth function.
`tests/unit/TestFiniteDifferences.m` checks the exact quadratic, independent
Dirichlet sine refinement and the periodic Fourier eigenvalue
$-4\sin^2(\pi k/n)/h^2$. Its refinement acceptance targets second-order error.
