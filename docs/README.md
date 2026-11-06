# MdeMat development

| Current state | Owner / next action |
|---|---|
| Milestone | Tools-layer implementation, repository API migration and native acceptance complete |
| Scope | Every tools-layer file plus affected data, execution, model, example, test and delivery entry points |
| Findings | Configuration, timing, sampling, component identity, graphics and persistence defects are repaired with independent regressions |
| Implementation | Data-owned plot methods; ordinary sample/style/result records; Config value plus pure parser; Timer/Visualizer resource owners; no Metric/Analyzer or old tools proxies |
| Next action | Use the current guides and measured profiles when extending numerical layers |
| Prior acceptance | The original v4 snapshot remains preserved separately; current tools gates accept the revised source |
| Delivery | Local MdeMat-2.0.0.mltbx built, exact payload/source bytes checked, fresh install and native project verified |
| Optional follow-up | Doxygen, other platforms and an authorized hosted CI run; DistMesh is bundled and its circle-mesh integration is tested |

Task status/dependencies are authoritative in [tasks.csv](tasks.csv). This table
is the current summary. Completion concerns the implemented scope and stated
numerical regimes; broader stability/AP/positivity claims require new evidence.
Repository quantities and recorded test/analyzer totals are maintained only in
[status.csv](status.csv), with scopes and evidence. Do not copy them into prose.

| Document | Purpose |
|---|---|
| [Architecture](architecture.md) | Layer responsibilities, numerical/data contracts and target layout |
| [Standard](standards.md) | Coding, tests, help, state/resources and release requirements |
| [Plan](plan.md) | Sequence, priority meaning, ledger semantics and final acceptance |
| [Tools refinement plan](tools-plan.md) | Findings, target contracts/boundaries, compatibility, staged work and verification |
| [Tools file review](tools-file-review.csv) | Every tools file's disposition, task, consumers and specific verification |
| [Tools review evidence](evidence/tools-layer-audit.json) | Source snapshot, reference candidates, native observations and limitations |
| [Tools design evidence](evidence/tools-design-review.json) | Second review, additional native probes, design choices and unchanged-source verification |
| [Tools implementation evidence](evidence/tools-implementation-verification.json) | Actual repairs, owning modules, independent regressions and delivery gates |
| [Tools capability map](evidence/tools-capability-map.json) | Explicit old methods/properties/options and data contracts mapped to implemented replacements |
| [Tools measurements](evidence/tools-study.json) | Validated native timings, IO bytes, graph artists and memory snapshots |
| [Task register](tasks.csv) | Work, dependencies, specific acceptance, status and evidence |
| [Source ledger](file-plan.csv) | Source identities, current/removed paths, targets and task ownership |
| [Baseline and defects](baseline.md) | Evidence interpretation and concrete source repairs |
| [Decisions](decisions.md) | Short accepted architecture/scope/cleanup decisions |
| [Low rank design](low-rank.md) | Implemented factor/operator contracts and bounded numerical progression |
| [Original evidence](evidence/baseline.json) | Historical analyzer, resources and unit outcomes |
| [Cleanup evidence](evidence/cleanup.json) | Actual removals and Gaussian preservation/native checks |
| [Plan coverage audit](evidence/plan-audit.json) | Independent pinned-v4 file, declaration, resource and task reconciliation |
| [Implementation audit](evidence/tools-migration-audit.json) | Current source/option/data-contract reconciliation; structural evidence only |
| [Prior v4 acceptance](evidence/v4-final-verification.json) | Preserved acceptance of the source snapshot before the tools redesign |
| [Mathematical formulations](mathematics/README.md) | Equations, numerical assumptions and implementation entry points |
| [Issue register](issues.csv) | Actual unresolved/optional environment limits and their follow-up |
| [Final acceptance](evidence/final-verification.json) | Native gates, artifact identity and scoped migration reconciliation |
| [Optimization review](optimization.md) | Measured costs, storage/error tradeoffs and bounded next improvements |

## Developer references

These files retain their authors'/publishers' rights. They are developer-only
references, explicitly excluded from the toolbox artifact by T073. Their presence
does not grant redistribution rights.

* [Johnson: The Elements of MATLAB Style](references/johnson_2010_elements_matlab_style.pdf).
* [Seibold and Frank: StaRMAP](references/seibold_frank_2014_starmap.pdf).
* [Peng et al.: oscillation-eliminating DG](references/peng_etal_2025_oedg.pdf).
* [Blaustein and Filbet: structure/AP VPFP scheme](references/blaustein_filbet_2024_vpfp_structure_ap.pdf).
