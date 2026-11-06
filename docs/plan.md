# Development plan

The source is the pinned v4 snapshot; work continues on main. Scope decisions
are in [decisions.md](decisions.md). Current progress and next work are maintained
only in [README.md](README.md) and [tasks.csv](tasks.csv). Quantities and measured
results live in [status.csv](status.csv), with raw evidence linked there.

## Sequence and ownership

| Milestone | Tasks | Outcome |
|---|---|---|
| M0 Baseline repair | T000, T010, T011, T012 | Traceable audit/cleanup, repaired defects and intended unit contracts |
| M1 Foundation | T020-T022 | Actual toolbox/project entry, small tools and supporting data contracts |
| M2 Common methods | T030, T032-T035 | All active algebra/JFNK, bases, fitting, quadrature, geometry, elements/spaces and assembly |
| M3 Elliptic family | T040-T041 | Complete Poisson LDG/PDG and independent assessment |
| M4 Evolution | T050-T053 | All ODE methods/history, upwind transport, diffusion and SLDG transfer |
| M5 Applications | T060-T062 | Kinetic representations, radiation and plasma composition/diagnostics |
| M6 Tools/delivery | T070-T076 | Persistence/restart, visuals, help/examples, benchmarks, measurements, package and CI |
| M7 Acceptance | T081-T082 | Verified namespace cutover and complete in-scope migration |
| M8 Follow-on | T090, T100-T103 | First FDM method and bounded low rank progression |
| M9 Tools refinement | T200-T290 | Tools contracts, repository API cutover, correctness and measured optimization; implementation and acceptance in [tools-plan.md](tools-plan.md) |

Poisson was the first complete solver-family slice; every other in-scope family
is included in final acceptance. Implement tests, help and examples with each task. Do not
create empty class trees. Only required prerequisites delay a task; exact
dependencies are in the task CSV rather than a second diagram.

M0-M8 retain their accepted source-snapshot evidence. The later tools review
identified additional boundary defects; it does not retroactively extend those
checks or certify repairs. M9 implements the revised tools architecture; its
current gate status is recorded in the task register. Its per-file dispositions are in [tools-file-review.csv](tools-file-review.csv),
separate from the v4 provenance ledger. The target replaces redundant classes,
moves model-specific plot methods to data, and removes tools proxies after
repository consumers migrate. Source API changes and continued reading of
existing v1 files are separate contracts; see [tools-plan.md](tools-plan.md).

Milestones describe sequence, not priority. P0 repairs blockers; P1 protects
required numerical functionality/acceptance; P2 covers usability and delivery;
P3 adds new methods. Priority orders ready tasks after dependencies are met.
T075 references, T076 measurements and T082 final acceptance remain P1 even
though delivered late. Universal checks are in [standards.md](standards.md);
task rows contain only their specific work and acceptance.

## File ledger

[file-plan.csv](file-plan.csv) contains v4 implementation/resource/build/support
identities plus implementation/support additions. It excludes developer docs and PDF
references. Git sourceBlob, symbol and declaration metadata replace the separate
inventory; there are no working-document SHA-256 or repeated analyzer totals.

sourcePath is an origin key, not a promise that the file exists. currentPath and
sourceState distinguish live files from removed originals. Removed records stay
for provenance and are recoverable from sourceCommit. The Gaussian duplicate
points to its working replacement; its cleanup is complete, while that example's
toolbox relocation is tracked separately from numerical acceptance. New code is tracked as a
working-tree addition without inventing a v4 blob identity.

KEEP preserves behavior with required relocation/verification. OPTIMIZE changes
contracts, implementation or ownership deliberately. REMOVE is a completed
cleanup disposition only where status/evidence establish removal. Ordinary
pending rows retain their target, task and specific rationale. Semicolons join
multiple paths/tasks. Source state and task completion are separate.

`RELOCATED_UNVERIFIED` means a retained origin has an actual destination;
`IMPLEMENTED_UNVERIFIED` marks additions awaiting their task's full gate. Neither
status certifies a mathematical result. The reproducible developer command
`python3 buildUtilities/auditMigration.py --refresh` reconciles current source
paths and additions; without `--refresh` it verifies the ledger. It never hashes
developer documents. Native execution evidence is maintained separately.

`VERIFIED` records completion of the assigned task gates for the accepted source
snapshot in [final-verification.json](evidence/final-verification.json).
`VERIFIED_WITH_OPTIONAL_LIMIT` retains an explicitly unavailable external
execution environment documented in [issues.csv](issues.csv). Neither status
claims every branch or arbitrary parameter regime was tested. Relevant source
changes require renewed evidence; `--refresh` reconciles paths and never reruns
MATLAB. The source fingerprint in the audit identifies the accepted tree.

Native MATLAB may replace generated `resources/project/**/*.xml` when project
membership changes. Refresh retires only missing NEW XML rows with no pinned v4
identity and with an unchanged original/current path. Their exact paths remain
in the audit's generated-metadata retirement record. Other missing files remain
errors; this rule never retires numerical code or a fixed v4 origin.

## Completion

T082 requires every live in-scope source/capability to be migrated or consolidated
into a verified migrated target, including all tests, examples, configurations,
resources, reference workflows and documented contracts. A removed duplicate
does not exempt its retained scenario or replacement from migration.

For each origin, account for its constructors/factories, methods/local helpers,
options, public/dependent properties, persisted fields and resource dependencies.
Retain, relocate or consolidate each responsibility, or document and verify its
intended replacement. File-level assignment or a renamed class alone does not
complete that capability. Required defect corrections preserve intended behavior,
not the original bug. Keep this mapping with the implementation change rather
than maintaining a second permanent method inventory.

All required family studies and contract/failure checks pass with recorded
evidence. No unresolved source defect, weakened assertion, blanket skip or
unimplemented command can be hidden by a successful pilot. FDM/low rank follow-on
tasks do not block this migration or replace unfinished existing functionality.

Rebuild the final package and verify help/resources/examples and required tests
in a fresh foreign-directory process with no old source or removed-directory
dependency. T073 excludes repository documentation and reference PDFs explicitly.
Update ledger, task status and evidence together when work completes.
