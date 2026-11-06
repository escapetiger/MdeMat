# Lebedev quadrature source

The retained `source/` files are the exact vendor copy in v4 commit `ee0d143`.
Their headers identify Dmitri Laikov and, for several table helpers, John
Burkardt. Author/reference headers are preserved. The project MIT license does
not replace this component's upstream license.

The [original upstream MATLAB page](https://people.math.sc.edu/Burkardt/m_src/sphere_lebedev_rule/sphere_lebedev_rule.html)
explicitly distributes its code and data under GNU LGPL. `upstream.html` captures
that page; per-file source URLs and SHA-256 digests are in `provenance.json`.
The upstream page does not state a license version. The GNU LGPL version 3 text
linked by that page and its corresponding GPL version 3 text are reproduced as
`COPYING.LESSER` and `COPYING`. No MIT relicensing is inferred from newer pages.

Compared with that upstream copy, v4 changes headings/parameter comments.
`xyz_to_tp.m` also changes two scalar multiplication operators to elementwise
operators; this existing change is recorded, not presented as upstream code.
Runtime files under `mdemat.internal.lebedev` additionally qualify internal
function calls with that package name. Full editable MATLAB runtime source,
original attribution, these notices and both license texts are distributed.

The retained original copy is provenance material and is never added to the
MATLAB path. The runtime package is the implementation used by LebedevRule.
