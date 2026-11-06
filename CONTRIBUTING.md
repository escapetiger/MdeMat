# Contributing

Start at [development status](docs/README.md), then follow the
[plan](docs/plan.md), [architecture](docs/architecture.md) and
[standard](docs/standards.md).

Find the affected source in [file-plan.csv](docs/file-plan.csv) and its task in
[tasks.csv](docs/tasks.csv). Add new implementation files to the source ledger;
developer documents do not need ledger entries or content hashes.

Make bounded changes with independent contract/numerical checks and native help.
Preserve valid test requirements; document intended API replacements and actual
PASS, FAIL or NOT RUN outcomes. Record source/task status and relevant evidence
together. Update the status entry when the next task changes.

The scope and cleanup decisions are in [decisions.md](docs/decisions.md).
Removed origins remain traceable to v4. After implementation changes, rerun the
affected gates and update the accepted source/evidence snapshot. Reference PDFs and all repository docs are excluded
from the distributed toolbox.
