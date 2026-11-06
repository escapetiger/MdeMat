#!/usr/bin/env python3
"""One-time, reviewed MATLAB exception migration (not a runtime dependency).

Usage: python3 buildUtilities/migrateToolExceptions.py --self-test
       python3 buildUtilities/migrateToolExceptions.py --apply --evidence FILE

The scanner preserves source spans, skips comments/quoted literals, recognizes
adjacent conjugate/nonconjugate transpose, and balances all MATLAB containers.
It fails before writing if a call has an unsupported arity or nonliteral ID.
Full source/ID mappings and before/after identities are retained in evidence.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


PREFIX = "mdemat.tools.except."
TARGETS = {PREFIX + "assert", PREFIX + "verify"}
NAMES = re.compile(r"[A-Za-z_]\w*(?:\.[A-Za-z_]\w*)*")
NUMBERS = re.compile(r"(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?[ij]?")


@dataclass(frozen=True)
class Token:
    kind: str
    start: int
    end: int
    text: str


def tokens(source: str) -> list[Token]:
    """Tokenize the source subset needed to locate and balance function calls."""
    result: list[Token] = []
    index = 0
    block_depth = 0
    while index < len(source):
        char = source[index]
        if char.isspace():
            index += 1
            continue
        if source.startswith("%{", index) and not source[source.rfind("\n", 0, index)+1:index].strip():
            block_depth += 1
            index += 2
            while index < len(source) and block_depth:
                line_end = source.find("\n", index)
                if line_end < 0:
                    line_end = len(source)
                line = source[index:line_end].strip()
                if line == "%{":
                    block_depth += 1
                elif line == "%}":
                    block_depth -= 1
                index = min(line_end + 1, len(source))
            if block_depth:
                raise ValueError("Unterminated MATLAB block comment")
            continue
        if char == "%" or source.startswith("...", index):
            line_end = source.find("\n", index)
            index = len(source) if line_end < 0 else line_end + 1
            continue
        previous = result[-1] if result else None
        adjacent_operand = previous is not None and previous.end == index and (
            previous.kind in {"name", "number", "string"}
            or previous.text in {")", "]", "}", "'", ".'"}
        )
        if source.startswith(".'", index) and previous is not None:
            result.append(Token("operator", index, index+2, ".'"))
            index += 2
            continue
        if char == "'" and adjacent_operand:
            result.append(Token("operator", index, index+1, "'"))
            index += 1
            continue
        if char in {"'", '"'}:
            start = index
            quote = char
            index += 1
            while index < len(source):
                if source[index] == "\n":
                    raise ValueError(f"Unterminated MATLAB literal at offset {start}")
                if source[index] == quote:
                    if index+1 < len(source) and source[index+1] == quote:
                        index += 2
                        continue
                    index += 1
                    break
                index += 1
            else:
                raise ValueError(f"Unterminated MATLAB literal at offset {start}")
            result.append(Token("string", start, index, source[start:index]))
            continue
        number = NUMBERS.match(source, index)
        if number:
            result.append(Token("number", index, number.end(), number.group()))
            index = number.end()
            continue
        name = NAMES.match(source, index)
        if name:
            result.append(Token("name", index, name.end(), name.group()))
            index = name.end()
            continue
        result.append(Token("operator", index, index+1, char))
        index += 1
    return result


def call_arguments(items: list[Token], opening: int) -> tuple[list[list[Token]], int]:
    stack: list[str] = []
    arguments: list[list[Token]] = []
    start = opening+1
    partners = {")": "(", "]": "[", "}": "{"}
    for index in range(opening+1, len(items)):
        text = items[index].text
        if items[index].kind == "string":
            continue
        if text in {"(", "[", "{"}:
            stack.append(text)
        elif text in partners:
            if not stack:
                if text != ")":
                    raise ValueError("Mismatched call delimiter")
                arguments.append(items[start:index])
                if any(not argument for argument in arguments):
                    raise ValueError("Empty exception argument")
                return arguments, index
            if stack.pop() != partners[text]:
                raise ValueError("Mismatched nested MATLAB delimiter")
        elif text == "," and not stack:
            arguments.append(items[start:index])
            start = index+1
    raise ValueError("Unclosed exception call")


def literal_id(argument: list[Token]) -> str:
    if len(argument) != 1 or argument[0].kind != "string":
        raise ValueError("Exception identifiers must be literal text scalars")
    text = argument[0].text
    return text[1:-1].replace(text[0]*2, text[0])


def namespace(path: Path) -> str:
    packages = [part[1:] for part in path.parts if part.startswith("+")]
    if packages:
        return ":".join(packages + [path.stem])
    category = "source"
    if "benchmarks" in path.parts:
        category = "benchmark"
    elif "buildUtilities" in path.parts:
        category = "build"
    elif "tests" in path.parts:
        category = "test"
    elif "examples" in path.parts:
        category = "example"
    return f"mdemat:{category}:{path.stem}"


def irrecoverable(path: Path, kind: str) -> bool:
    return kind == "verify" and path.name in {"Separable.m", "LinearSolver.m"}


def migrate(source: str, path: Path) -> tuple[str, list[dict]]:
    items = tokens(source)
    edits: list[tuple[int, int, str]] = []
    records: list[dict] = []
    for index, token in enumerate(items):
        if token.kind != "name" or token.text not in TARGETS:
            continue
        if index+1 >= len(items) or items[index+1].text != "(":
            raise ValueError("Unexpected qualified exception reference without a call")
        arguments, closing = call_arguments(items, index+1)
        kind = token.text.rsplit(".", 1)[1]
        if len(arguments) not in {1} and len(arguments) < 3:
            raise ValueError(f"Unsupported {kind} arity {len(arguments)}")
        reason = literal_id(arguments[1]) if len(arguments) >= 3 else (
            "AssertionFailed" if kind == "assert" else "VerificationFailed"
        )
        if not re.fullmatch(r"[A-Za-z]\w*(?::[A-Za-z]\w*)*", reason):
            raise ValueError(f"Invalid literal exception identifier {reason!r}")
        full_id = reason if ":" in reason else f"{namespace(path)}:{reason}"
        args = [source[argument[0].start:argument[-1].end] for argument in arguments]
        mode = "assert" if kind == "assert" or irrecoverable(path, kind) else "warning"
        if mode == "assert":
            if len(arguments) == 1:
                replacement = f"assert(all({args[0]}, 'all'), '{full_id}', '%s', 'Assertion failed.')"
                edits.append((token.start, items[closing].end, replacement))
            else:
                edits.append((token.start, token.end, "assert"))
                edits.append((arguments[0][0].start, arguments[0][-1].end, f"all({args[0]}, 'all')"))
                edits.append((arguments[1][0].start, arguments[1][-1].end, f"'{full_id}'"))
                if len(arguments) == 3:
                    edits.append((arguments[2][0].start, arguments[2][0].start, "'%s', "))
        else:
            line_start = source.rfind("\n", 0, token.start)+1
            indentation = source[line_start:token.start]
            if indentation.strip():
                raise ValueError("A warning call must be a standalone statement")
            if len(arguments) == 1:
                message_args = "'%s', 'Verification failed.'"
            elif len(arguments) == 3:
                message_args = "'%s', " + args[2]
            else:
                message_args = ", ".join(args[2:])
            replacement = (
                f"if ~all({args[0]}, 'all')\n"
                f"{indentation}    warning('{full_id}', {message_args});\n"
                f"{indentation}end"
            )
            edits.append((token.start, items[closing].end, replacement))
        records.append({
            "line": source.count("\n", 0, token.start)+1,
            "column": token.start-(source.rfind("\n", 0, token.start)+1)+1,
            "originalTool": kind,
            "nativeTool": mode,
            "reason": reason,
            "oldIdentifier": f"{namespace(path)}:{reason}",
            "newIdentifier": full_id,
            "flag": args[0],
            "originalCall": source[token.start:items[closing].end],
            "unrecoverableWarningPromotedToFailure": irrecoverable(path, kind),
        })
    # No call in this project takes another void exception call as an argument.
    ordered = sorted(edits)
    for previous, current in zip(ordered, ordered[1:]):
        if previous[1] > current[0]:
            raise ValueError("Overlapping/nested exception replacements require manual review")
    for start, end, replacement in reversed(ordered):
        source = source[:start] + replacement + source[end:]
    # Re-tokenize output before any filesystem mutation.
    if any(token.kind == "name" and token.text in TARGETS for token in tokens(source)):
        raise ValueError("Unconverted exception call remains")
    return source, records


def scalar_conditions(source: str, identifiers: set[str]) -> tuple[str, int]:
    """Repair only explicitly mapped native assertions, preserving empty truth."""
    items = tokens(source)
    edits = []
    for index, token in enumerate(items):
        if token.kind != 'name' or token.text != 'assert' or index+1 >= len(items) or items[index+1].text != '(':
            continue
        arguments, _ = call_arguments(items, index+1)
        if len(arguments) < 3 or literal_id(arguments[1]) not in identifiers:
            continue
        condition = arguments[0]
        if len(condition) > 2 and condition[0].text == 'all' and condition[1].text == '(':
            inner, closing = call_arguments(condition, 1)
            if (closing == len(condition)-1 and len(inner) == 2 and len(inner[1]) == 1
                    and inner[1][0].kind == 'string' and literal_id(inner[1]) == 'all'):
                continue
        text = source[condition[0].start:condition[-1].end]
        edits.append((condition[0].start, condition[-1].end, f"all({text}, 'all')"))
    for start, end, text in reversed(edits):
        source = source[:start]+text+source[end:]
    return source, len(edits)


def excluded(path: Path) -> bool:
    parts = path.parts
    return "third_party" in parts or "+except" in parts or (
        "+tools" in parts and "+visual" in parts
    ) or ("+data" in parts and "+mesh" in parts and (
        "+internal" in parts or (path.name == "plot.m" and any(p.startswith("@") for p in parts))
    ))


def self_test() -> None:
    source = """% mdemat.tools.except.assert(false, 'Ignore', 'comment');
%{
mdemat.tools.except.assert(false, 'Ignore', 'block comment');
%}
literal = "mdemat.tools.except.assert(false, 'Ignore', 'text')";
mdemat.tools.except.assert(norm(A' * B.' - f(x,y)) < tol, ...
    'Mismatch', ['can''t ', char("a,b")], x(1,2));
mdemat.tools.except.assert(all(isfield(s, {'a','b'})), 'Fields', '100% literal');
mdemat.tools.except.verify(cellfun(@(x) isa(x,'double'), values), ...
    'Recoverable', 'text');
"""
    result, records = migrate(source, Path("toolbox/+mdemat/+base/Probe.m"))
    assert len(records) == 3
    assert "norm(A' * B.' - f(x,y))" in result
    assert "'mdemat:base:Probe:Mismatch'" in result
    assert "'mdemat:base:Probe:Fields', '%s', '100% literal'" in result
    assert "if ~all(cellfun(@(x) isa(x,'double'), values), 'all')" in result
    assert records[-1]["nativeTool"] == "warning"
    promoted, records = migrate(
        "mdemat.tools.except.verify(iscell(x), 'Bad', 'Value must be a cell.');\n",
        Path("toolbox/+mdemat/+data/+variable/Separable.m"),
    )
    assert promoted.startswith("assert(")
    assert records[0]["unrecoverableWarningPromotedToFailure"]


def digest(source: str) -> str:
    return hashlib.sha256(source.encode("utf-8")).hexdigest()


def atomic_write(path: Path, source: str) -> None:
    mode = path.stat().st_mode
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", newline="", dir=path.parent,
                                     suffix=".exception-migration.tmp", delete=False) as output:
        temporary = Path(output.name)
        output.write(source)
    os.chmod(temporary, mode)
    os.replace(temporary, path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--apply", action="store_true")
    parser.add_argument("--remove-implementations", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--evidence", type=Path)
    options = parser.parse_args()
    if options.remove_implementations and not options.apply:
        parser.error("--remove-implementations requires --apply")
    self_test()
    if options.self_test:
        print("MATLAB lexical migration self-tests passed.")
        return
    root = options.root.resolve()
    files: set[Path] = set(root.glob("*.m"))
    for folder in ("toolbox", "tests", "benchmarks", "buildUtilities"):
        files.update((root/folder).rglob("*.m"))
    changes: list[tuple[Path, str, str]] = []
    reports: list[dict] = []
    for path in sorted(files):
        relative = path.relative_to(root)
        if excluded(relative):
            continue
        source = path.read_bytes().decode("utf-8")
        if PREFIX not in source:
            continue
        try:
            migrated, calls = migrate(source, relative)
        except Exception as exception:
            raise RuntimeError(f"{relative}: {exception}") from exception
        if not calls:
            continue
        changes.append((path, source, migrated))
        reports.append({"path": relative.as_posix(), "beforeSha256": digest(source),
                        "afterSha256": digest(migrated), "calls": calls})
    retired = []
    if options.remove_implementations:
        for relative in ("toolbox/+mdemat/+tools/+except/assert.m", 
                         "toolbox/+mdemat/+tools/+except/verify.m", 
                         "toolbox/+mdemat/+tools/+except/private/backtrack.m"):
            path = root/relative
            if path.is_file():
                retired.append({"path": relative,
                                "beforeSha256": hashlib.sha256(path.read_bytes()).hexdigest()})
        transformed = {path: source for path, _, source in changes}
        for path in files:
            if "third_party" in path.parts or "+except" in path.parts:
                continue
            source = transformed.get(path, path.read_bytes().decode("utf-8"))
            if PREFIX in source and any(t.kind == "name" and t.text in TARGETS for t in tokens(source)):
                raise RuntimeError(f"Cannot retire wrappers while excluded calls remain: {path}")
    report = {
        "schemaVersion": 1,
        "capturedUTC": datetime.now(timezone.utc).isoformat(),
        "mode": "apply" if options.apply else "preview",
        "applied": False,
        "scope": ["root MATLAB files", "toolbox", "tests", "benchmarks", "buildUtilities"],
        "exclusions": ["third_party", "tools.except implementations", "tools.visual", 
                       "data.mesh @Class/plot.m", "data.mesh.internal", "historical docs/evidence"],
        "identifierPolicy": "Explicit package:owner:reason, preserving valid preexisting IDs.",
        "legacyIdentifierEvidence": "Source package/file plus literal reason; legacy stack/path-dependent IDs may differ for local functions or installations. Native regression remains required.",
        "files": reports,
        "removedImplementations": retired,
        "fileCount": len(reports),
        "callCount": sum(len(item["calls"]) for item in reports),
        "verification": {"lexicalSelfTestsPassed": True, "nativeTestsPending": True},
    }
    if options.evidence:
        options.evidence.parent.mkdir(parents=True, exist_ok=True)
        options.evidence.write_text(json.dumps(report, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    if options.apply:
        for path, source, _ in changes:
            if path.read_bytes().decode("utf-8") != source:
                raise RuntimeError(f"Concurrent source change detected before migration: {path}")
        for path, source, migrated in changes:
            if path.read_bytes().decode("utf-8") != source:
                raise RuntimeError(f"Concurrent source change detected during migration: {path}")
            atomic_write(path, migrated)
        for item in retired:
            (root/item["path"]).unlink()
        for relative in ("toolbox/+mdemat/+tools/+except/private", "toolbox/+mdemat/+tools/+except"):
            path = root/relative
            if path.is_dir() and not any(path.iterdir()):
                path.rmdir()
        report["applied"] = True
        if options.evidence:
            options.evidence.write_text(json.dumps(report, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    print(json.dumps({"mode": report["mode"], "files": report["fileCount"],
                      "calls": report["callCount"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
