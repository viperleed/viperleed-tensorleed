"""Writing frozen test cases and editing the titles and namelists of their inputs.

Shared by the case generators; imports only the standard library and
baseline_tree, so it runs without viperleed.
"""
import re
import shutil

from .baseline_tree import GeneratorError


def write_case(case_dir, header, program, description, covers, runs, inputs):
    """Write case_dir/case.toml and the input files, replacing case_dir/inputs.

    runs: dicts with name, stdin, param, outputs and optionally files
    ({key: path}); inputs: {path: text}. Paths are relative to case_dir.
    """
    if (case_dir / "inputs").exists():
        shutil.rmtree(case_dir / "inputs")
    (case_dir / "inputs").mkdir(parents=True)
    for path, text in inputs.items():
        (case_dir / path).write_text(text)
    toml = [header, f'program = "{program}"', f'description = "{description}"',
            "covers = [" + ", ".join(f'"{c}"' for c in covers) + "]", ""]
    for run in runs:
        toml += ["[[runs]]", f'name = "{run["name"]}"', f'stdin = "{run["stdin"]}"',
                 f'param = "{run["param"]}"']
        if run.get("files"):
            toml.append("files = { " + ", ".join(f'{key} = "{path}"'
                                                for key, path in run["files"].items()) + " }")
        toml += ["outputs = [" + ", ".join(f'"{o}"' for o in run["outputs"]) + "]", ""]
    (case_dir / "case.toml").write_text("\n".join(toml))


def set_first_line_title(text, title):
    """Replace the first line (the title of a theory spectrum or a ref-calc FIN)."""
    _, rest = text.split("\n", maxsplit=1)
    return title + "\n" + rest


def set_wexpel_title(text, title):
    """Replace the WEXPEL title line (first line after the &NL3 ... &END block)."""
    lines = text.split("\n")
    i_nl3 = next(i for i, l in enumerate(lines) if l.strip() == "&NL3")
    i_end = next(i for i in range(i_nl3, len(lines)) if lines[i].strip() == "&END")
    lines[i_end + 1] = f" {title}"
    return "\n".join(lines)


def set_namelist_value(text, key, value):
    """Set KEY= value, in a namelist (viperleed layout: one key per line)."""
    if key == "WR":
        new, n = re.subn(r"(?m)^(\s*WR=\s*).*$", rf"\g<1>{value}", text)
    else:
        new, n = re.subn(rf"(?m)^(\s*{key}=\s*)[^,\n]*,", rf"\g<1>{value},", text)
    if n != 1:
        raise GeneratorError(f"expected one '{key}=' line in WEXPEL, found {n}")
    return new
