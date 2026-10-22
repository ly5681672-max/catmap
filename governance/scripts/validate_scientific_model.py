"""Fail-closed check of the actual CatMAP model selected for an authorized run."""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
import os
import shutil
import tempfile
from pathlib import Path


HERE = Path(__file__).resolve()
DEFAULT_ROOT = HERE.parents[2]
MANIFEST = HERE.parents[1] / "SCIENTIFIC_MODEL_BASELINE.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def literal_assignments(model_path: Path) -> dict:
    tree = ast.parse(model_path.read_text(encoding="utf-8-sig"), filename=str(model_path))
    result = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            try:
                result[node.targets[0].id] = ast.literal_eval(node.value)
            except (ValueError, TypeError):
                continue
    return result


def parse_with_catmap(model_path: Path, energy_path: Path, energy_filename: str) -> dict:
    """Parse in a disposable directory so CatMAP cannot write beside source inputs."""
    import catmap
    from catmap import ReactionModel

    with tempfile.TemporaryDirectory(prefix="catmap_model_preflight_") as tmp:
        tmpdir = Path(tmp)
        local_model = tmpdir / model_path.name
        shutil.copy2(model_path, local_model)
        shutil.copy2(energy_path, tmpdir / energy_filename)
        old_cwd = Path.cwd()
        try:
            os.chdir(tmpdir)
            model = ReactionModel(setup_file=local_model.name)
            model.parse()
            model.verify()
            adsorbates = list(model.adsorbate_names)
            ts_names = list(model.transition_state_names)
            gases = list(model.gas_names)
            site_counts = {
                species: model.species_definitions[species].get("n_sites", 1)
                for species in adsorbates
            }
            gas_inputs = {
                gas: {
                    key: model.species_definitions[gas][key]
                    for key in ("pressure", "concentration")
                    if key in model.species_definitions[gas]
                }
                for gas in gases
            }
            parsed_rxns = []
            for rxn in model.elementary_rxns:
                parsed_rxns.append({
                    "initial": list(rxn[0]),
                    "transition_state": list(rxn[1]) if len(rxn) > 2 else None,
                    "final": list(rxn[-1]),
                })
            return {
                "catmap_version": catmap.__version__,
                "parsed_reaction_count": len(parsed_rxns),
                "parsed_reactions": parsed_rxns,
                "adsorbates": sorted(adsorbates),
                "transition_states": sorted(ts_names),
                "gases": sorted(gases),
                "n_sites": site_counts,
                "temperature_K": model.temperature,
                "pressure_mode": model.pressure_mode,
                "gas_inputs": gas_inputs,
            }
        finally:
            os.chdir(old_cwd)


def validate(root: Path, manifest: dict | None = None, model_arg: Path | None = None,
             parsed: dict | None = None) -> list[str]:
    manifest = manifest or json.loads(MANIFEST.read_text(encoding="utf-8"))
    root = root.resolve()
    approved_model = (root / manifest["source_model"]).resolve()
    model_path = (model_arg if model_arg is not None else approved_model).resolve()
    energy_path = (root / manifest["energy_input"]).resolve()
    errors: list[str] = []

    if model_path != approved_model:
        errors.append(f"selected model path is not the approved formal input: {model_path}")
    if not model_path.is_file():
        errors.append(f"missing selected model: {model_path}")
        return errors
    if sha256(model_path).lower() != manifest["model_sha256"].lower():
        errors.append(f"selected model SHA256 differs from frozen baseline: {model_path}")
    if not energy_path.is_file():
        errors.append(f"missing energy input: {energy_path}")
    elif sha256(energy_path).lower() != manifest["energy_input_sha256"].lower():
        errors.append(f"energy input SHA256 differs from frozen baseline: {energy_path}")

    try:
        assignments = literal_assignments(model_path)
        expressions = assignments.get("rxn_expressions")
        if expressions != manifest["reaction_expressions"]:
            errors.append("source reaction expressions differ from frozen formal eight-stage mechanism")
        if assignments.get("temperature") != manifest["temperature_K"]:
            errors.append("source temperature differs from the frozen baseline")
        if assignments.get("pressure_mode") != manifest["pressure_mode"]:
            errors.append("source pressure_mode differs from the frozen baseline")
        definitions = assignments.get("species_definitions", {})
        gas_inputs = {}
        for gas in manifest["gas_inputs"]:
            definition = definitions.get(gas, {})
            if not definition:
                # The formal file may construct this mapping through indexed assignments.
                import re
                source = model_path.read_text(encoding="utf-8-sig")
                match = re.search(r"species_definitions\[['\"]" + re.escape(gas) + r"['\"]\]\s*=\s*\{([^}]*)\}", source)
                if match:
                    for key, value in re.findall(r"['\"](pressure|concentration)['\"]\s*:\s*([^,\n]+)", match.group(1)):
                        try:
                            definition[key] = ast.literal_eval(value.strip())
                        except (ValueError, SyntaxError):
                            definition[key] = value.strip()
            gas_inputs[gas] = {k: definition[k] for k in ("pressure", "concentration") if k in definition}
        if gas_inputs != manifest["gas_inputs"]:
            errors.append(f"source gas reservoir fields/values differ from baseline: {gas_inputs}")
        if "OOH*" not in " ".join(expressions or []):
            errors.append("explicit OOH* intermediate is absent")
    except (SyntaxError, ValueError) as exc:
        errors.append(f"cannot inspect selected model source: {exc}")

    if parsed is not None:
        if parsed["catmap_version"] != manifest["catmap_version"]:
            errors.append(f"CatMAP version differs from baseline: {parsed['catmap_version']}")
        if parsed["parsed_reaction_count"] != len(manifest["reaction_expressions"]):
            errors.append(f"CatMAP parsed {parsed['parsed_reaction_count']} reactions; expected 8")
        if set(parsed["adsorbates"]) != {s.replace("*", "_s") for s in manifest["stable_adsorbates"]}:
            errors.append(f"CatMAP parsed stable adsorbates differ from baseline: {parsed['adsorbates']}")
        if set(parsed["transition_states"]) != {"OOH-C6H12_s", "Hb-O_s"}:
            errors.append(f"CatMAP parsed transition states differ from baseline: {parsed['transition_states']}")
        if any(float(n) != 1 for n in parsed["n_sites"].values()):
            errors.append(f"CatMAP parsed non-baseline n_sites: {parsed['n_sites']}")
        if float(parsed["temperature_K"]) != float(manifest["temperature_K"]):
            errors.append(f"CatMAP parsed temperature differs from baseline: {parsed['temperature_K']}")
        if parsed["pressure_mode"] != manifest["pressure_mode"]:
            errors.append(f"CatMAP parsed pressure_mode differs from baseline: {parsed['pressure_mode']}")
        if parsed["gas_inputs"] != manifest["gas_inputs"]:
            errors.append(f"CatMAP parsed gas inputs differ from baseline: {parsed['gas_inputs']}")
        if not any("OOH_s" in (r["initial"] + (r["transition_state"] or []) + r["final"])
                   for r in parsed["parsed_reactions"]):
            errors.append("CatMAP parsed network does not contain OOH_s")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--model", type=Path, help="actual .mkm selected for the upcoming run")
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    args = parser.parse_args()
    root = args.project_root.resolve()
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    if args.model is None:
        model_path = (root / manifest["source_model"]).resolve()
    elif args.model.is_absolute():
        model_path = args.model.resolve()
    else:
        model_path = (root / args.model).resolve()
    energy_path = root / manifest["energy_input"]
    errors = validate(root, manifest, model_arg=model_path)
    parsed = None
    if model_path.is_file() and energy_path.is_file():
        try:
            parsed = parse_with_catmap(model_path, energy_path, Path(manifest["energy_input"]).name)
            errors.extend(validate(root, manifest, model_arg=model_path, parsed=parsed))
        except Exception as exc:  # fail closed: parser errors prevent the run
            errors.append(f"CatMAP parser/verify failed: {type(exc).__name__}: {exc}")
    if errors:
        print("MODEL_CONTRACT=BLOCKED")
        for error in dict.fromkeys(errors):
            print(f"- {error}")
        if parsed is not None:
            print("ACTUAL_CATMAP_PARSE=" + json.dumps(parsed, ensure_ascii=False, sort_keys=True))
        return 2
    print("MODEL_CONTRACT=PASS")
    print(f"selected_model={model_path}")
    print(f"model_sha256={manifest['model_sha256']}")
    print(f"energy_sha256={manifest['energy_input_sha256']}")
    print("ACTUAL_CATMAP_PARSE=" + json.dumps(parsed, ensure_ascii=False, sort_keys=True))
    print("No steady-state/TOF/DRC calculation was executed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
