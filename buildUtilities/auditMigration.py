#!/usr/bin/env python3
"""Reconcile pinned v4 origins and current implementation without hashing docs.

python3 buildUtilities/auditMigration.py --refresh updates live-path/new-file rows.
The normal command writes docs/evidence/tools-migration-audit.json. The earlier
migration-audit.json is retained as historical acceptance/source evidence.
This is a developer audit; MATLAB users do not need Python or Git at runtime.
"""
from __future__ import annotations

import argparse
import collections
import csv
import datetime as dt
import hashlib
import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
COMMIT = "ee0d143074c4e426e96086082b142dd31c5309d2"
LEDGER = ROOT / "docs/file-plan.csv"
ROOT_FILES = {".gitattributes", ".gitignore", "LICENSE", "README.md", "CONTRIBUTING.md",
              "MdeMat.prj", "setupMdeMat.m", "runTests.m", "buildfile.m"}
TREES = {"toolbox", "tests", "benchmarks", "buildUtilities", "third_party", "resources", ".github"}
FUNCTION = re.compile(r"^\s*function\s+(?:(?:\[[^\]]*\]|\w+)\s*=\s*)?([\w.]+)", re.M)
CLASS = re.compile(r"^\s*classdef\s+(?:\([^)]*\)\s*)?\w+\s*(?:<\s*([^\n%]+))?", re.M)
NATIVE = {"assert", "error", "warning", "MException", "cell", "struct", "fieldnames",
          "isfield", "numel", "isempty", "setfield", "rmfield", "onCleanup", "table",
          "matlab.unittest.fixtures.TemporaryFolderFixture", "plot", "loglog"}


def references(value: str) -> list[str]:
    return list(filter(None, value.split(';')))


def native_reference(path: str) -> bool:
    return path.startswith('builtin:')


def source_fingerprint(identities: list[dict[str, str]]) -> str:
    return hashlib.sha256(json.dumps(sorted(identities,key=lambda item:item['path']),
        sort_keys=True,separators=(',',':')).encode()).hexdigest()


def resolve_class(symbol: str) -> Path | None:
    parts = symbol.split('.')
    folder = ROOT/'toolbox'/Path(*('+'+part for part in parts[:-1]))
    for candidate in (folder/(parts[-1]+'.m'), folder/('@'+parts[-1])/(parts[-1]+'.m')):
        if candidate.is_file():
            return candidate
    return None


def implementation_sources(paths: list[str]) -> dict[str, str]:
    result = {}
    for path in paths:
        if native_reference(path) or not (ROOT/path).is_file():
            continue
        file = ROOT/path
        result[path] = file.read_text(errors='replace') if file.suffix == '.m' else ''
        if file.parent.name.startswith('@') and file.stem == file.parent.name[1:]:
            for external in file.parent.glob('*.m'):
                result[str(external.relative_to(ROOT))] = external.read_text(errors='replace')
    return result


def declarations_with_parents(text: str, visited: set[str] | None = None) -> set[str]:
    visited = set() if visited is None else visited
    names = set(FUNCTION.findall(text))
    for match in CLASS.finditer(text.replace('...\n', ' ')):
        for symbol in re.findall(r'\bmdemat\.[\w.]+', match[1] or ''):
            if symbol in visited:
                continue
            parent = resolve_class(symbol)
            if parent:
                inherited = implementation_sources([str(parent.relative_to(ROOT))])
                names |= declarations_with_parents('\n'.join(inherited.values()), visited|{symbol})
    return names


def tools_tasks(path: str) -> set[str]:
    lower = path.lower()
    result = set()
    if '/+tools/+base/' in path: result.add('T220' if path.endswith('/hash.m') else 'T210')
    if '/+tools/+environment/' in path:
        result.add('T230' if Path(path).stem in {'Config', 'readConfig'} else 'T220')
    if '/+tools/+visual/' in path or '/+data/+mesh/@' in path or '/+mesh/+internal/' in path:
        result.add('T260')
    if '/+tools/+analyze/' in path:
        result.add('T220' if Path(path).stem == 'Timer' else 'T250')
        result.add('T300')
    if '/+tools/+io/' in path or any(word in lower for word in (
            'checkpoint', 'bindstatecallbacks', 'validatestaterecord')): result.add('T240')
    if '/+fem/' in path and any(word in lower for word in (
            'samplingplan', 'samplefield', 'referencevalues', 'metric')): result.add('T250')
    if any(word in lower for word in ('plotsampling', 'plotfieldsamples')): result.add('T260')
    if '/+execution/' in path:
        stem = Path(path).stem
        if stem in {'Scheme','OdeScheme','defineStudyInputs','configureObservation'}: result.add('T230')
        if stem in {'assessState','runConvergence','convergenceTables'}: result.add('T250')
        if stem in {'configurePlotLayouts','sampleObservation','interpolateObservation',
                    'plotReferenceFields','sampledDataset'}: result.add('T260')
    if path.startswith('tests/'):
        if any(word in lower for word in ('config','schemehierarchy')): result.add('T230')
        if any(word in lower for word in ('basetools','foundationprobe')): result.add('T210')
        if 'timer' in lower: result.add('T220')
        if any(word in lower for word in ('analysis','convergence','sampling')): result.add('T250')
        if any(word in lower for word in ('render','visual','database','dataset','tiled','figure')): result.add('T260')
        if any(word in lower for word in ('persistence','recordio')): result.add('T240')
        result.add('T290')
    if path.endswith('migrateToolExceptions.py'): result.add('T210')
    if path.endswith('Configuration.m'): result.add('T230')
    if result or path.startswith('buildUtilities/'): result.add('T270')
    return result


def declared_tests(text: str) -> set[str]:
    """Only methods in a Test block are regression references, not helpers."""
    result = set()
    active_indent = None
    for line in text.splitlines():
        header = re.match(r'^(\s*)methods\s*(?:\(([^)]*)\))?', line)
        if header:
            active_indent = len(header[1]) if re.search(r'\bTest\b', header[2] or '') else None
            continue
        if active_indent is not None:
            if line.strip() == 'end' and len(line)-len(line.lstrip()) == active_indent:
                active_indent = None
            else:
                result.update(FUNCTION.findall(line))
    return result


def git(*args: str) -> bytes:
    return subprocess.check_output(["git", *args], cwd=ROOT)


def task_ids(path: str) -> str:
    lower = path.lower()
    if path == 'tests/numerical/TestModeConvergence.m':
        return 'T250;T270;T290'
    if "retainedstudies" in lower:
        return "T050;T051;T052;T061;T062;T070;T072"
    if "retainedutilitystudies" in lower:
        return "T034;T050;T072"
    if "bdfnonlinear" in lower or "splitcallbacks" in lower:
        return "T050;T070"
    if "retainedcontrollers" in lower:
        return "T040;T041;T050;T051;T052;T053"
    if "angularstudies" in lower:
        return "T032;T061;T072"
    if "phifunctions" in lower:
        return "T050"
    if "relativeanalysis" in lower:
        return "T041"
    if "studyconfig" in lower:
        return "T072"
    if "lowranktransport" in lower or "evolveadaptive" in lower:
        return "T103"
    if "evolvefixedrank" in lower:
        return "T102"
    if "/+lowrank/" in lower or "testfactor" in lower:
        return "T101"
    if "/+tensor/" in lower or "applyseparated" in lower:
        return "T100;T101"
    if "/+fdm/" in lower or "finitedifference" in lower:
        return "T090"
    if "ganapol" in lower:
        return "T075"
    if path.startswith(".github/"):
        return "T074"
    if "distmesh" in lower or "nonuniform_trimesh" in lower:
        return "T034;T072;T073"
    if "lebedev" in lower or "third_party" in lower:
        return "T033;T073"
    if "symbolic" in lower or "basisresource" in lower:
        return "T032"
    if "documentation" in lower or path.startswith("toolbox/doc/"):
        return "T072"
    if path.startswith("resources/") or path == "MdeMat.prj":
        return "T020"
    if path.startswith("buildUtilities/") or path in {"buildfile.m", "runTests.m", "setupMdeMat.m"}:
        return "T020;T073;T074"
    if "lowrankdynamics" in lower:
        return "T102;T103"
    if "evolution" in lower or "/+ode/" in lower:
        return "T050;T070"
    if "quadrature" in lower:
        return "T033"
    if "geometry" in lower:
        return "T034"
    if "assembly" in lower or "/+fem/" in lower:
        return "T035;T070"
    if "poisson" in lower or "/+elliptic/" in lower:
        return "T040"
    if "/+hyperbolic/" in lower:
        return "T051;T053"
    if "/+parabolic/" in lower:
        return "T052"
    if "application" in lower or "/+radiation/" in lower or "/+plasma/" in lower:
        return "T060;T061;T062"
    if "/+io/" in lower or "recordio" in lower or "snapshot" in lower:
        return "T070"
    if "/+analyze/" in lower or "/+analysis/" in lower or "testanalysis" in lower:
        return "T041;T076"
    if "/+visual/" in lower:
        return "T071"
    if "core" in lower or "spherical" in lower:
        return "T010;T032"
    if "contract" in lower:
        return "T011"
    if "linearalgebra" in lower:
        return "T030"
    return "T082"


def current_files() -> set[str]:
    result = {name for name in ROOT_FILES if (ROOT/name).is_file()}
    for tree in TREES:
        result.update(str(p.relative_to(ROOT)) for p in (ROOT/tree).rglob("*")
                      if p.is_file() and "__pycache__" not in p.parts)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true")
    options = parser.parse_args()
    report_path = ROOT/'docs/evidence/tools-migration-audit.json'
    prior_report = json.loads(report_path.read_text()) if report_path.is_file() else {}
    retired_metadata = set(prior_report.get('retiredGeneratedProjectMetadata',[]))
    retired_now = []
    with LEDGER.open(newline="") as stream:
        reader = csv.DictReader(stream)
        columns = reader.fieldnames
        rows = list(reader)
    if options.refresh:
        retained = []
        for row in rows:
            path = row['sourcePath']
            if (row['sourceState'] == 'NEW' and row['sourceCommit'] != COMMIT and
                    path.startswith('resources/project/') and Path(path).suffix.lower() == '.xml' and
                    row['currentPath'] == path and not (ROOT/path).exists()):
                # MATLAB regenerates its project metadata. This narrow exception
                # retires only absent, generated records without any v4 identity.
                retired_metadata.add(path)
                retired_now.append(path)
            else:
                retained.append(row)
        rows = retained
    mapping = json.loads((ROOT/"docs/evidence/namespace-migration.json").read_text())["files"]
    mapping.update({"Doxyfile": "buildUtilities/documentation/Doxyfile",
                    "m2cpp.pl": "buildUtilities/documentation/m2cpp.pl"})
    cutover = json.loads((ROOT/'docs/evidence/tools-cutover-map.json').read_text())['files']
    layout = json.loads((ROOT/'docs/evidence/tools-layout-map.json').read_text())['files']
    capability = json.loads((ROOT/'docs/evidence/tools-capability-map.json').read_text())
    previous = json.loads((ROOT/'docs/evidence/migration-audit.json').read_text())
    baseline = {item['path']: item['sha256'] for item in previous['identities']}
    accepted_file = ROOT/'docs/evidence/tools-accepted-identities.json'
    accepted = json.loads(accepted_file.read_text()) if accepted_file.is_file() else {}
    task_states = {row['id']:row['status'] for row in
                   csv.DictReader((ROOT/'docs/tasks.csv').open())}
    accepted_paths = {item['path']:item['sha256'] for item in accepted.get('identities',[])}
    acceptance_enabled = (accepted.get('schemaVersion') == 1 and
                          accepted.get('status') == 'ACCEPTED' and
                          accepted.get('task') == 'T290' and task_states.get('T290') == 'DONE' and
                          accepted.get('evidence') == 'docs/evidence/final-verification.json' and
                          (ROOT/accepted['evidence']).is_file() and
                          len(accepted_paths) == len(accepted.get('identities',[])) and
                          accepted.get('sourceTreeFingerprint') == source_fingerprint(
                              accepted.get('identities',[])))
    exception_changes = json.loads((ROOT/'docs/evidence/tools-exception-migration.json').read_text())
    exception_paths = {item['path'] for item in exception_changes['files']}
    target_tasks = collections.defaultdict(set)
    for item in cutover.values():
        for path in item['targets']:
            target_tasks[path].update(item['tasks'])

    def targets_for(row):
        origin = row['sourcePath']
        initial = mapping.get(origin, origin)
        if initial in cutover:
            return cutover[initial]['targets']
        if origin in cutover:
            return cutover[origin]['targets']
        # A current explicit cutover must take precedence over an older namespace map.
        current = references(row['currentPath'])
        if current and any(path in cutover for path in current):
            return list(dict.fromkeys(target for path in current
                for target in cutover.get(path, {'targets':[path]})['targets']))
        return [initial] if origin in mapping else current
    known = {row["sourcePath"] for row in rows}
    physical = current_files()
    accepted_whole_tree = acceptance_enabled and set(accepted_paths) == physical and all(
        hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == digest
        for path,digest in accepted_paths.items())
    if options.refresh:
        for row in rows:
            origin = row["sourcePath"]
            if row["sourceState"] == "REMOVED":
                if "gaussian-DESKTOP" in origin:
                    row["currentPath"] = "toolbox/examples/radiation/ex03_gaussian.m"
                    row["targetPath"] = row["currentPath"]
                continue
            targets = [layout.get(path, path) for path in targets_for(row)]
            row['currentPath'] = ';'.join(targets)
            row['targetPath'] = row['currentPath']
            changed = any(native_reference(path) or not (ROOT/path).is_file() or
                baseline.get(path) != hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
                for path in targets)
            accepted_targets = acceptance_enabled and all(native_reference(path) or
                ((ROOT/path).is_file() and accepted_paths.get(path) ==
                 hashlib.sha256((ROOT/path).read_bytes()).hexdigest()) for path in targets)
            # Native replacements have no individual source file. Their acceptance
            # is tied to the complete independent source snapshot below.
            if accepted_targets and any(native_reference(path) for path in targets):
                accepted_targets = accepted_whole_tree
            tasks_for_row = set(references(row['taskIds']))
            for path in targets:
                tasks_for_row.update(target_tasks[path])
                if changed: tasks_for_row.update(tools_tasks(path))
                if path in exception_paths: tasks_for_row.update({'T210','T270'})
            row['taskIds'] = ';'.join(sorted(tasks_for_row))
            if changed:
                row['migrationStatus'] = 'VERIFIED' if accepted_targets else 'IMPLEMENTED_UNVERIFIED'
                evidence = references(row['evidence'])
                evidence += ['evidence/tools-cutover-map.json', 'evidence/tools-capability-map.json',
                             'evidence/tools-migration-audit.json']
                row['evidence'] = ';'.join(dict.fromkeys(evidence))
                note = ('Tools cutover accepted against independent T290 source identities.' if
                        accepted_targets else
                        'Tools cutover implemented; current source API/acceptance awaits full T290 gates.')
                stale = ('Tools cutover implemented; current source API/acceptance awaits full T290 gates.'
                         if accepted_targets else
                         'Tools cutover accepted against independent T290 source identities.')
                row['notes'] = row['notes'].replace(stale,'').strip()
                if note not in row['notes']: row['notes'] += ' '+note
                if accepted_targets:
                    row['evidence'] += ';evidence/tools-accepted-identities.json;evidence/final-verification.json'
                    row['evidence'] = ';'.join(dict.fromkeys(references(row['evidence'])))
        owned = {p for row in rows for p in references(row['currentPath']) if not native_reference(p)}
        for path in sorted(physical-owned):
            if path in known:
                raise ValueError(f"Origin exists but is assigned elsewhere: {path}")
            source = ROOT/path
            text = source.read_text(errors="replace") if source.suffix == ".m" else ""
            kind = "class" if re.search(r"^\s*classdef\b", text, re.M) else (
                "function" if source.suffix == ".m" else "asset")
            row = dict.fromkeys(columns, "")
            row.update(sourcePath=path, kind=kind, role="new-implementation-or-support",
                       symbol=source.stem if source.suffix == ".m" else "",
                       functionDeclarations=str(len(FUNCTION.findall(text))), decision="KEEP",
                       sourceState="NEW", currentPath=path, targetPath=path,
                       taskIds=';'.join(sorted(set(references(task_ids(path))) | tools_tasks(path) |
                                             target_tasks[path])), migrationStatus="IMPLEMENTED_UNVERIFIED",
                       notes="Added during implementation; no invented v4 Git identity.")
            rows.append(row)
        with LEDGER.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    origins = set(git("ls-tree", "-r", "--name-only", COMMIT).decode().splitlines())
    origins = {p for p in origins if not p.startswith("docs/")}
    recorded_origins = {row["sourcePath"] for row in rows if row["sourceCommit"] == COMMIT}
    owned = {p for row in rows for p in references(row['currentPath']) if not native_reference(p)}
    tasks = {r["id"] for r in csv.DictReader((ROOT/"docs/tasks.csv").open())}
    errors = []
    for path in sorted(origins-recorded_origins): errors.append("Unaccounted v4 origin: "+path)
    for path in sorted(recorded_origins-origins): errors.append("Invalid v4 origin: "+path)
    for path in sorted(physical-owned): errors.append("Unregistered current source: "+path)
    for row in rows:
        if row["sourceState"] != "REMOVED":
            for path in row["currentPath"].split(";"):
                if native_reference(path):
                    if path.partition(':')[2] not in NATIVE: errors.append('Unknown native target: '+path)
                elif not path or not (ROOT/path).is_file(): errors.append("Missing current path: "+path)
        for task in filter(None, row["taskIds"].split(";")):
            if task not in tasks: errors.append("Unknown task: "+task)
        if row["sourceCommit"] == COMMIT:
            actual = git("rev-parse", COMMIT+":"+row["sourcePath"]).decode().strip()
            if actual != row["sourceBlob"]: errors.append("Source identity changed: "+row["sourcePath"])
    live = sorted(physical)
    # Hash code/resources/support, never developer documents or this evidence.
    identities = [{"path": p, "sha256": hashlib.sha256((ROOT/p).read_bytes()).hexdigest()} for p in live]
    declaration_gaps = []
    test_methods = collections.defaultdict(list)
    for path in (ROOT/'tests').rglob('*.m'):
        for name in declared_tests(path.read_text(errors='replace')):
            test_methods[f'{path.stem}/{name}'].append(str(path.relative_to(ROOT)))
    snapshot = json.loads((ROOT/'docs/evidence/tools-layer-audit.json').read_text())
    original_declarations = {}
    original_text = {}
    for origin in recorded_origins:
        if origin.endswith('.m'):
            original_text[origin] = git('show',COMMIT+':'+origin).decode(errors='replace')
            names = set(FUNCTION.findall(original_text[origin]))
            original_declarations[origin] = names
            original_declarations.setdefault(mapping.get(origin,origin),set()).update(names)
    for item in snapshot['files']:
        original_declarations.setdefault(item['path'],set()).update(
            name for declaration in item['declarations']
            for name in FUNCTION.findall(declaration['text']))
    explained = {}
    capability_errors = []

    def validate_evidence(entry, label, require_documentation=False):
        if not entry.get('replacements') or not entry.get('tests'):
            capability_errors.append('Incomplete capability evidence: '+label)
        for replacement in entry.get('replacements',[]):
            if 'native' in replacement:
                if replacement['native'] not in NATIVE or not replacement.get('operation'):
                    capability_errors.append('Unsupported native capability: '+label)
                continue
            path = replacement.get('path','')
            if not (ROOT/path).is_file():
                capability_errors.append('Missing capability replacement: '+label+' -> '+path)
                continue
            text = '\n'.join(implementation_sources([path]).values())
            if replacement.get('declaration'):
                if replacement['declaration'] not in declarations_with_parents(text):
                    capability_errors.append('Missing replacement declaration: '+label+' -> '+str(replacement))
            elif not replacement.get('operation') or not replacement.get('anchor') or replacement['anchor'] not in text:
                capability_errors.append('Value capability lacks a real source anchor: '+label)
        for test in entry.get('tests',[]):
            if len(test_methods[test]) != 1:
                capability_errors.append('Missing/ambiguous capability regression: '+label+' -> '+test)
        if require_documentation and not entry.get('documentation'):
            capability_errors.append('Missing capability documentation: '+label)
        for path in entry.get('documentation',[]):
            if not (ROOT/path).is_file():
                capability_errors.append('Missing capability documentation: '+label+' -> '+path)

    for entry in capability['methods']:
        label = entry['source']+'#'+entry['name']
        if not entry.get('capability') or not entry.get('replacements') or not entry.get('tests'):
            capability_errors.append('Incomplete exact capability mapping: '+label)
            continue
        if entry['name'] not in original_declarations.get(entry['source'],set()):
            capability_errors.append('Mapping is not an original declaration: '+label)
        validate_evidence(entry,label)
        for source in entry.get('aliases', [entry['source']]):
            if source != entry['source'] and mapping.get(source) != entry['source']:
                capability_errors.append('Unrelated capability source alias: '+source+' -> '+entry['source'])
            key = (source,entry['name'])
            if key in explained: capability_errors.append('Duplicate capability mapping: '+str(key))
            explained[key] = entry

    # Properties and name-value options are capabilities too. A missing method
    # check cannot establish retention of a value field, option or wire contract.
    baseline_interfaces = {item['source']:item for item in capability['interfacesBaseline']}
    expected_interfaces = set()
    if set(baseline_interfaces) != {item['path'] for item in snapshot['files']}:
        capability_errors.append('Interface baseline does not cover the exact original 44 files.')
    for item in snapshot['files']:
        before = baseline_interfaces.get(item['path'],{})
        if before.get('sha256Before') != item['sha256']:
            capability_errors.append('Interface source identity mismatch: '+item['path'])
        for kind,field in [('property','properties'),('option','options')]:
            for member in before.get(field,[]):
                if not member.get('name') or not member.get('line') or not (member.get('declaration') or member.get('text')):
                    capability_errors.append('Interface baseline lacks source declaration: '+item['path'])
                expected_interfaces.add((item['path'],kind,member.get('name')))
        for block in item.get('propertyBlocks',[]):
            for member in block['members']:
                name = re.match(r'^([A-Za-z]\w*)\b',member['text'])
                if name and (item['path'],'property',name[1]) not in expected_interfaces:
                    capability_errors.append('Original property missing from baseline: '+item['path']+'#'+name[1])
        for member in item.get('optionLines',[]):
            name = re.match(r'^options\.([A-Za-z]\w*)\b',member['text'])
            if name and (item['path'],'option',name[1]) not in expected_interfaces:
                capability_errors.append('Original name-value option missing from baseline: '+item['path']+'#'+name[1])
    covered_interfaces = set()
    for entry in capability['interfaceGroups']:
        label = entry['source']+'#'+entry['kind']+':'+','.join(entry['names'])
        if not entry.get('targetFieldsOrResponsibilities') or not entry['names']:
            capability_errors.append('Missing interface capability responsibility: '+label)
        validate_evidence(entry,label,True)
        for name in entry['names']:
            key = (entry['source'],entry['kind'],name)
            if key not in expected_interfaces:
                capability_errors.append('Interface mapping has no baseline declaration: '+str(key))
            if key in covered_interfaces:
                capability_errors.append('Duplicate interface mapping: '+str(key))
            covered_interfaces.add(key)
    for key in sorted(expected_interfaces-covered_interfaces):
        capability_errors.append('Unmapped property/name-value capability: '+str(key))
    for entry in capability['wireContracts']:
        label = 'wire:'+entry['name']
        for field in ['keys','owner','targets','compatibility','newChecks','tests']:
            if not entry.get(field): capability_errors.append('Incomplete wire contract: '+label+'#'+field)
        for path in entry.get('targets',[]):
            if not (ROOT/path).is_file(): capability_errors.append('Missing wire contract owner: '+label+' -> '+path)
        for test in entry.get('tests',[]):
            if len(test_methods[test]) != 1:
                capability_errors.append('Missing/ambiguous wire regression: '+label+' -> '+test)
    errors.extend(capability_errors)
    with (ROOT/'docs/tools-file-review.csv').open(newline='') as stream:
        reviews = list(csv.DictReader(stream))
    if len(reviews) != len(snapshot['files']) or {row['path'] for row in reviews} != set(baseline_interfaces):
        errors.append('Implementation review must retain exactly the original 44 tool files.')
    for row in reviews:
        path = row['path']
        if references(row.get('currentPaths','')) != cutover.get(path,{}).get('targets',[]):
            errors.append('Implementation review differs from current cutover: '+path)
        if row.get('implementationStatus') not in {'IMPLEMENTED','RETAINED'}:
            errors.append('Implementation review lacks actual disposition: '+path)
        for evidence in references(row.get('implementationEvidence','')):
            if not (ROOT/'docs'/evidence).is_file():
                errors.append('Implementation review evidence missing: '+path+' -> '+evidence)
    checked = 0
    inherited = []
    for row in rows:
        origin = row["sourcePath"]
        if row["sourceCommit"] != COMMIT or not origin.endswith('.m') or row["sourceState"] == "REMOVED":
            continue
        before = original_text[origin]
        after = '\n'.join(implementation_sources(references(row['currentPath'])).values())
        local_missing = set(FUNCTION.findall(before))-set(FUNCTION.findall(after))
        explained_names = {name for name in local_missing if (origin,name) in explained or
                           (mapping.get(origin,origin),name) in explained}
        inherited_names = (local_missing-explained_names) & declarations_with_parents(after,set())
        if inherited_names:
            inherited.append({'origin':origin, 'names':sorted(inherited_names),
                              'disposition':'Preserved through explicit current superclass implementation.'})
        missing = local_missing-inherited_names
        missing = {name for name in missing if (origin,name) not in explained and
                   (mapping.get(origin,origin),name) not in explained}
        if missing:
            declaration_gaps.append({"origin": origin, "missingNames": sorted(missing)})
        checked += 1
    tool_gaps = []
    for file in snapshot['files']:
        source = file['path']
        targets = cutover.get(source, {'targets':[source]})['targets']
        after = '\n'.join(implementation_sources(targets).values())
        before_names = set()
        for item in file['declarations']:
            before_names.update(FUNCTION.findall(item['text']))
        missing = before_names-declarations_with_parents(after)
        missing = sorted(name for name in missing if (source,name) not in explained)
        if missing: tool_gaps.append({'source':source,'missingNames':missing})
    if declaration_gaps:
        errors.append("Unexplained original method/helper declaration gaps remain.")
    if tool_gaps:
        errors.append('Unexplained pre-tools-upgrade declaration gaps remain.')
    report = {"schemaVersion": 2, "capturedUTC": dt.datetime.now(dt.UTC).isoformat(),
              "sourceCommit": COMMIT, "scope": "v4 non-doc origins, excluding removed legacy; current implementation/support closure",
              "status": "PASS" if not errors else "FAIL", "errors": errors,
              "originRows": len(recorded_origins), "currentFiles": len(live),
              "currentMatlabFiles": sum(p.endswith('.m') for p in live),
              "ledgerStates": dict(collections.Counter(r["migrationStatus"] for r in rows)),
              "declarationFilesCompared": checked, "unexplainedDeclarationGaps": declaration_gaps,
              "toolsBaselineFilesCompared":len(snapshot['files']),
              "unexplainedToolsDeclarationGaps":tool_gaps,
              "exactCapabilityMappings":len(capability['methods']),
              "interfaceCapabilityGroups":len(capability['interfaceGroups']),
              "originalPropertyAndOptionNames":len(expected_interfaces),
              "wireContractGroups":len(capability['wireContracts']),
              "implementationReviewRows":len(reviews),
              "retiredGeneratedProjectMetadata":sorted(retired_metadata),
              "retiredGeneratedProjectMetadataCount":len(retired_metadata),
              "newlyRetiredGeneratedProjectMetadataCount":len(retired_now),
              "projectMetadataRetirementPolicy":"Only absent NEW resources/project XML rows without the pinned v4 source identity and with unchanged currentPath are retired; all source gaps remain errors.",
              "capabilityValidationErrors":capability_errors,
              "implementationAcceptance":"SEPARATE_NATIVE_EVIDENCE_REQUIRED",
              "nativeAcceptanceEvidence":"docs/evidence/final-verification.json",
              "acceptedIdentitySnapshot":("MATCHED" if accepted_whole_tree else
                  "STALE" if acceptance_enabled else "NOT_FINALIZED" if accepted else "NOT_AVAILABLE"),
              "inheritedMethodConsolidations": inherited,
              "capabilityEvidence":"docs/evidence/tools-capability-map.json",
              "identities": identities,
              "sourceTreeFingerprint": source_fingerprint(identities),
              "limits": "Structural reconciliation only. Relocation does not establish mathematical or runtime acceptance."}
    report_path.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({k:v for k,v in report.items() if k not in {"identities"}}, indent=2))
    return int(bool(errors))


if __name__ == "__main__":
    raise SystemExit(main())
