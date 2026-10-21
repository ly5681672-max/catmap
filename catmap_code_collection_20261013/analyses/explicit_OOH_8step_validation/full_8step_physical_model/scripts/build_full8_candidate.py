"""Build an isolated, explicitly conditional full-state CatMAP candidate.

No source scientific data are edited. Composite-state site counts are formal
fragment-bookkeeping assumptions only and are not claimed as validated geometry.
"""
from __future__ import annotations

import csv
import hashlib
import json
import shutil
from collections import Counter
from fractions import Fraction
from functools import reduce
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parents[2]
SOURCE = PROJECT / "analyses" / "explicit_OOH_8step_validation" / "candidate_tico_kinetic_control" / "candidate_energies.txt"
ORIGINAL_MKM = PROJECT / "analyses" / "explicit_OOH_8step_validation" / "multisite_failure_and_full_model_diagnostics" / "inputs" / "original_8step_snapshot" / "alkene_epoxidation.mkm"
OLD_CANDIDATE_MKM = PROJECT / "analyses" / "explicit_OOH_8step_validation" / "candidate_tico_kinetic_control" / "candidate_tico.mkm"
MODEL_OUT = ROOT / "inputs" / "FULL_8STEP_MODEL.mkm"
ENERGY_OUT = ROOT / "inputs" / "FULL_8STEP_TICO_energies.txt"
DATA = ROOT / "data"
REPORTS = ROOT / "reports"

SPECIES = {
    "H2O2_s": ({"H": 2, "O": 2}, 1, {}),
    "Ha_s": ({"H": 1}, 1, {"Ha": 1}),
    "Hb_s": ({"H": 1}, 1, {"Hb": 1}),
    "OOH_s": ({"H": 1, "O": 2}, 1, {"OOH": 1}),
    "C6H12_s": ({"C": 6, "H": 12}, 1, {"C6H12": 1}),
    "C6H12O_s": ({"C": 6, "H": 12, "O": 1}, 1, {"C6H12O": 1}),
    "O_s": ({"O": 1}, 1, {"O": 1}),
    "OH_s": ({"H": 1, "O": 1}, 1, {"OH": 1}),
    "H2O_s": ({"H": 2, "O": 1}, 1, {"H2O": 1}),
    "NEB1IS_s": ({"C": 6, "H": 14, "O": 2}, 3, {"Ha": 1, "OOH": 1, "C6H12": 1}),
    "NEB1TS_s": ({"C": 6, "H": 14, "O": 2}, 3, {"Ha": 1, "OOH": 1, "C6H12": 1}),
    "NEB1FS_s": ({"C": 6, "H": 14, "O": 2}, 4, {"Ha": 1, "Hb": 1, "O": 1, "C6H12O": 1}),
    "NEB2IS_s": ({"H": 2, "O": 1}, 3, {"Ha": 1, "Hb": 1, "O": 1}),
    "NEB2TS_s": ({"H": 2, "O": 1}, 3, {"Ha": 1, "Hb": 1, "O": 1}),
    "NEB2FS_s": ({"H": 2, "O": 1}, 2, {"Ha": 1, "OH": 1}),
}
GAS = {
    "H2O2_g": {"H": 2, "O": 2},
    "C6H12_g": {"C": 6, "H": 12},
    "C6H12O_g": {"C": 6, "H": 12, "O": 1},
    "H2O_g": {"H": 2, "O": 1},
}
RXNS = [
    "* + H2O2_g <-> H2O2*",
    "H2O2* + * -> Ha* + OOH*",
    "* + C6H12_g <-> C6H12*",
    "Ha* + OOH* + C6H12* <-> NEB1IS*",
    "NEB1IS* + * <-> NEB1TS* + * -> NEB1FS*",
    "NEB1FS* <-> Ha* + C6H12O* + O* + Hb*",
    "C6H12O* <-> C6H12O_g + *",
    "Ha* + Hb* + O* <-> NEB2IS*",
    "NEB2IS* <-> NEB2TS* -> NEB2FS* + *",
    "NEB2FS* <-> Ha* + OH*",
    "Ha* + OH* -> H2O* + *",
    "H2O* <-> H2O_g + *",
]


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def nullspace(a: list[list[int]]) -> list[list[int]]:
    m = [[Fraction(x) for x in row] for row in a]
    nr, nc = len(m), len(m[0]) if m else 0
    pivots: list[int] = []
    r = 0
    for c in range(nc):
        p = next((i for i in range(r, nr) if m[i][c]), None)
        if p is None:
            continue
        m[r], m[p] = m[p], m[r]
        pivot = m[r][c]
        m[r] = [x / pivot for x in m[r]]
        for i in range(nr):
            if i != r and m[i][c]:
                q = m[i][c]
                m[i] = [x - q * y for x, y in zip(m[i], m[r])]
        pivots.append(c)
        r += 1
    free = [c for c in range(nc) if c not in pivots]
    basis = []
    for f in free:
        v = [Fraction(0) for _ in range(nc)]
        v[f] = 1
        for i, p in enumerate(pivots):
            v[p] = -m[i][f]
        d = reduce(lcm, (x.denominator for x in v), 1)
        ints = [int(x * d) for x in v]
        g = reduce(gcd, (abs(x) for x in ints if x), 0) or 1
        ints = [x // g for x in ints]
        if next(x for x in ints if x) < 0:
            ints = [-x for x in ints]
        basis.append(ints)
    return basis


def parse_energy_table() -> tuple[list[str], dict[str, float], dict[str, str]]:
    rows = list(csv.DictReader(SOURCE.open(encoding="utf-8-sig", newline=""), delimiter="\t"))
    selected = [
        r for r in rows
        if (r["surface_name"] == "tico" and r["site_name"] == "111")
        or (r["surface_name"] == "None" and r["site_name"] == "gas")
    ]
    energies = {}
    references = {}
    for r in selected:
        name = r["species_name"] + ("_g" if r["site_name"] == "gas" else "")
        if name in energies:
            raise ValueError(f"duplicate energy row for {name}")
        energies[name] = float(r["formation_energy"])
        references[name] = r["reference"]
    required = {n.removesuffix("_s") for n in SPECIES} | set(GAS)
    missing = sorted(required - energies.keys())
    if missing:
        raise ValueError(f"missing required energies: {missing}")
    ENERGY_OUT.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(SOURCE, ENERGY_OUT)
    return selected, energies, references


def main() -> None:
    ROOT.mkdir(parents=True, exist_ok=True)
    for d in (DATA, REPORTS, ROOT / "inputs", ROOT / "runs", ROOT / "figures"):
        d.mkdir(parents=True, exist_ok=True)
    selected, energies, references = parse_energy_table()
    mkm = """# CONDITIONAL FULL-STATE MODEL; not formally physically validated.
# Every NEB state keeps its Ha spectator in the mapped total formula.
# Composite n_sites values are bookkeeping counts from fragment projection,
# not DFT-validated geometric site footprints. Missing-link TSs use CatMAP's
# native endpoint treatment only as a declared conditional approximation.
rxn_expressions = %r
surface_names = ['tico']
temperature = 333.15
pressure_mode = 'concentration'
scaler = 'ThermodynamicScaler'
species_definitions = {
    's': {'site_names': ['111'], 'total': 1},
    'Ha_s': {'composition': {'H': 1}, 'n_sites': 1},
    'Hb_s': {'composition': {'H': 1}, 'n_sites': 1},
    'OOH_s': {'composition': {'H': 1, 'O': 2}, 'n_sites': 1},
    'NEB1IS_s': {'composition': {'C': 6, 'H': 14, 'O': 2}, 'n_sites': 3},
    'NEB1TS_s': {'composition': {'C': 6, 'H': 14, 'O': 2}, 'n_sites': 3},
    'NEB1FS_s': {'composition': {'C': 6, 'H': 14, 'O': 2}, 'n_sites': 4},
    'NEB2IS_s': {'composition': {'H': 2, 'O': 1}, 'n_sites': 3},
    'NEB2TS_s': {'composition': {'H': 2, 'O': 1}, 'n_sites': 3},
    'NEB2FS_s': {'composition': {'H': 2, 'O': 1}, 'n_sites': 2},
    'C6H12_g': {'concentration': 1.0},
    'H2O2_g': {'concentration': 0.185},
    'H2O_g': {'concentration': 0.815},
    'C6H12O_g': {'concentration': 1e-20},
}
descriptor_names = ['temperature', 'pressure']
descriptor_ranges = None
resolution = 1
input_file = 'FULL_8STEP_TICO_energies.txt'
data_file = 'FULL_8STEP_MODEL.pkl'
gas_thermo_mode = 'frozen_gas'
adsorbate_thermo_mode = 'frozen_adsorbate'
avoid_scaling = True
output_variables = ['turnover_frequency', 'forward_rate', 'reverse_rate', 'coverage']
use_numbers_solver = False
decimal_precision = 120
tolerance = 1e-50
max_rootfinding_iterations = 250
""" % RXNS
    MODEL_OUT.write_text(mkm, encoding="utf-8")

    # CatMAP 0.3.1 parser defines the actual elementary reaction sequence.
    import catmap
    from catmap import ReactionModel
    old_cwd = Path.cwd()
    import os
    os.chdir(MODEL_OUT.parent)
    try:
        model = ReactionModel(setup_file=MODEL_OUT.name)
        model.parse()
        reactions = list(model.elementary_rxns)
        names = list(model.adsorbate_names)
        species_lookup = {name: name.removesuffix("_s") for name in names}
        matrix = [[0 for _ in reactions] for _ in names]
        for j, rxn in enumerate(reactions):
            initial, final = Counter(rxn[0]), Counter(rxn[-1])
            for i, name in enumerate(names):
                matrix[i][j] = final.get(name, 0) - initial.get(name, 0)
        matrix_rows = []
        for i, name in enumerate(names):
            key = species_lookup[name]
            formula, nsites, _ = SPECIES[key + "_s"]
            matrix_rows.append({"species": name, "composition": json.dumps(formula, sort_keys=True),
                                "formal_n_sites_assumption": nsites,
                                "reaction_projection": json.dumps(SPECIES[key + "_s"][2], sort_keys=True),
                                **{f"r{j+1}": matrix[i][j] for j in range(len(reactions))}})
        write_csv(DATA / "FULL_8STEP_STOICHIOMETRY.csv", matrix_rows)

        # Parse the archived original eight-step input independently and emit
        # its actual CatMAP matrix for a direct, source-versioned comparison.
        if not ORIGINAL_MKM.is_file():
            raise FileNotFoundError(f"original eight-step snapshot missing: {ORIGINAL_MKM}")
        os.chdir(ORIGINAL_MKM.parent)
        original = ReactionModel(setup_file=ORIGINAL_MKM.name)
        original.parse()
        original_reactions = list(original.elementary_rxns)
        original_names = list(original.adsorbate_names)
        original_matrix = [[0 for _ in original_reactions] for _ in original_names]
        for j, rxn in enumerate(original_reactions):
            left, right = Counter(rxn[0]), Counter(rxn[-1])
            for i, name in enumerate(original_names):
                original_matrix[i][j] = right.get(name, 0) - left.get(name, 0)
        original_rows = []
        for i, name in enumerate(original_names):
            definition = original.species_definitions[name]
            original_rows.append({"species": name,
                                  "composition": json.dumps(definition.get("composition", {}), sort_keys=True),
                                  "n_sites": definition.get("n_sites", 1),
                                  **{f"r{j+1}": original_matrix[i][j] for j in range(len(original_reactions))}})
        write_csv(DATA / "ORIGINAL_8STEP_STOICHIOMETRY.csv", original_rows)
        original_inv = nullspace([list(row) for row in zip(*original_matrix)])
        write_csv(DATA / "ORIGINAL_8STEP_INVENTORY_BASIS.csv", [
            {"basis_id": f"I{i}", "expression": " + ".join(
                f"({c})theta({n})" for c, n in zip(vec, original_names) if c) + " = constant",
             **{n: c for n, c in zip(original_names, vec)}}
            for i, vec in enumerate(original_inv, 1)])

        # Parse the historical nine-reaction candidate so its lost inventory
        # dimension is quantified rather than only described.
        os.chdir(OLD_CANDIDATE_MKM.parent)
        old_candidate = ReactionModel(setup_file=OLD_CANDIDATE_MKM.name)
        old_candidate.parse()
        old_reactions = list(old_candidate.elementary_rxns)
        old_names = list(old_candidate.adsorbate_names)
        old_matrix = [[0 for _ in old_reactions] for _ in old_names]
        for j, rxn in enumerate(old_reactions):
            left, right = Counter(rxn[0]), Counter(rxn[-1])
            for i, name in enumerate(old_names):
                old_matrix[i][j] = right.get(name, 0) - left.get(name, 0)
        write_csv(DATA / "OLD_NINE_REACTION_STOICHIOMETRY.csv", [
            {"species": name,
             "composition": json.dumps(old_candidate.species_definitions[name].get("composition", {}), sort_keys=True),
             "n_sites": old_candidate.species_definitions[name].get("n_sites", 1),
             **{f"r{j+1}": old_matrix[i][j] for j in range(len(old_reactions))}}
            for i, name in enumerate(old_names)])
        old_inv = nullspace([list(row) for row in zip(*old_matrix)])
        write_csv(DATA / "OLD_NINE_REACTION_INVENTORY_BASIS.csv", [
            {"basis_id": f"J{i}", "expression": " + ".join(
                f"({c})theta({n})" for n, c in zip(old_names, vec) if c) + " = constant",
             **{n: c for n, c in zip(old_names, vec)}}
            for i, vec in enumerate(old_inv, 1)])
        os.chdir(MODEL_OUT.parent)

        # Left-nullspace is null(A^T), where A is species x elementary reaction.
        inv = nullspace([list(row) for row in zip(*matrix)])
        inv_rows = []
        for i, vec in enumerate(inv, 1):
            inv_rows.append({"basis_id": f"J{i}", "expression": " + ".join(
                f"({c})theta({n})" for c, n in zip(vec, names) if c) + " = constant",
                **{n: c for n, c in zip(names, vec)}})
        write_csv(DATA / "FULL_8STEP_INVENTORY_BASIS.csv", inv_rows)

        # Verify atom and formal site balance on each elementary CatMAP reaction.
        formula_by_ident = {n.removesuffix("_s"): f for n, (f, _, _) in SPECIES.items()}
        formula_by_ident.update(GAS)
        nsites_by_ident = {n.removesuffix("_s"): ns for n, (_, ns, _) in SPECIES.items()}
        balance = []
        for j, rxn in enumerate(reactions, 1):
            sides = []
            for side in (rxn[0], rxn[-1]):
                atom = Counter()
                sites = 0
                for ident in side:
                    if ident == "s":
                        sites += 1
                        continue
                    base = ident.removesuffix("_s").removesuffix("_g")
                    for element, count in formula_by_ident[base].items():
                        atom[element] += count
                    if not ident.endswith("_g"):
                        sites += nsites_by_ident[base]
                sides.append((atom, sites))
            atom_ok = sides[0][0] == sides[1][0]
            site_ok = sides[0][1] == sides[1][1]
            balance.append({"elementary_reaction": j, "equation": model.print_rxn(rxn, mode="text"),
                            "atom_balance": "PASS" if atom_ok else "FAIL",
                            "formal_site_balance": "PASS" if site_ok else "FAIL",
                            "left_formal_sites": sides[0][1], "right_formal_sites": sides[1][1],
                            "left_atoms": json.dumps(dict(sorted(sides[0][0].items())), sort_keys=True),
                            "right_atoms": json.dumps(dict(sorted(sides[1][0].items())), sort_keys=True)})
        write_csv(DATA / "FULL_8STEP_BALANCE_CHECKS.csv", balance)

        def catmap_state_energy(state: list[str]) -> float:
            total = 0.0
            for ident in state:
                if ident in model.site_names:
                    continue
                value = model.species_definitions[ident]["formation_energy"]
                total += float(value[0] if isinstance(value, (list, tuple)) else value)
            return total

        reaction_energy_rows = []
        for j, rxn in enumerate(reactions, 1):
            g_is = catmap_state_energy(list(rxn[0]))
            g_fs = catmap_state_energy(list(rxn[-1]))
            delta_g = g_fs - g_is
            if len(rxn) > 2 and rxn[1]:
                g_ts = catmap_state_energy(list(rxn[1]))
                barrier_f, barrier_r = g_ts - g_is, g_ts - g_fs
                barrier_basis = "EXPLICIT_DFT_NEB_TS"
            else:
                g_ts = max(g_is, g_fs)
                barrier_f, barrier_r = max(0.0, delta_g), max(0.0, -delta_g)
                barrier_basis = "CATMAP_NATIVE_ENDPOINT_MAX_APPROXIMATION"
            reaction_energy_rows.append({
                "reaction_index": j, "reaction": model.print_rxn(rxn, mode="text"),
                "G_initial_eV": g_is, "G_TS_eV": g_ts, "G_final_eV": g_fs,
                "deltaG_forward_eV": delta_g, "forward_barrier_eV": barrier_f,
                "reverse_barrier_eV": barrier_r, "barrier_basis": barrier_basis,
            })
        write_csv(DATA / "REACTION_ENERGY_AND_BARRIERS.csv", reaction_energy_rows)
    finally:
        os.chdir(old_cwd)

    # Ti-Co endpoint energy map; values are formation-energy sums under the
    # archived common elemental reference. Their structural matching remains
    # conditional and is stated in the status column.
    def d(names_l: list[str], names_r: list[str]) -> float:
        return sum(energies[n] for n in names_r) - sum(energies[n] for n in names_l)

    connection_rows = [
        {"connection": "independent Ha+OOH+C6H12 -> full NEB1-IS", "deltaG_eV": d(["Ha", "OOH", "C6H12"], ["NEB1IS"]), "data_class": "DERIVED_ENDPOINT_FREE_ENERGY", "kinetics": "missing connector TS; endpoint approximation only conditionally", "mapping": "Ha spectator included in NEB state formula; geometry/site footprint unverified"},
        {"connection": "NEB1-IS -> NEB1-FS", "deltaG_eV": energies["NEB1FS"]-energies["NEB1IS"], "data_class": "DIRECT_NEB_ENDPOINTS", "kinetics": "NEB1TS explicit; signed barriers from IS/TS/FS", "mapping": "full NEB state energies"},
        {"connection": "full NEB1-FS -> Ha+C6H12O+O+Hb", "deltaG_eV": d(["NEB1FS"], ["Ha", "C6H12O", "O", "Hb"]), "data_class": "DERIVED_ENDPOINT_FREE_ENERGY", "kinetics": "missing dissociation/reorganization TS; endpoint approximation only conditionally", "mapping": "product/O/Hb fragments plus retained Ha spectator"},
        {"connection": "NEB1-FS -> NEB2-IS + C6H12O_g", "deltaG_eV": d(["NEB1FS"], ["NEB2IS", "C6H12O_g"]), "data_class": "DERIVED_ENDPOINT_FREE_ENERGY", "kinetics": "no independent connector TS", "mapping": "published candidate endpoint delta"},
        {"connection": "Ha+Hb+O -> full NEB2-IS", "deltaG_eV": d(["Ha", "Hb", "O"], ["NEB2IS"]), "data_class": "DERIVED_ENDPOINT_FREE_ENERGY", "kinetics": "missing association/reorganization TS; endpoint approximation only conditionally", "mapping": "Ha spectator included; geometry/site footprint unverified"},
        {"connection": "NEB2-IS -> NEB2-FS", "deltaG_eV": energies["NEB2FS"]-energies["NEB2IS"], "data_class": "DIRECT_NEB_ENDPOINTS", "kinetics": "NEB2TS explicit; signed barriers from IS/TS/FS", "mapping": "full NEB state energies"},
        {"connection": "full NEB2-FS -> Ha+OH", "deltaG_eV": d(["NEB2FS"], ["Ha", "OH"]), "data_class": "DERIVED_ENDPOINT_FREE_ENERGY", "kinetics": "missing deaggregation TS", "mapping": "Ha spectator retained"},
        {"connection": "full NEB2-FS -> H2O*", "deltaG_eV": energies["H2O"]-energies["NEB2FS"], "data_class": "DERIVED_ENDPOINT_FREE_ENERGY", "kinetics": "missing water-formation/reorganization TS", "mapping": "existing endpoint delta"},
        {"connection": "Ha+OH -> H2O*", "deltaG_eV": d(["Ha", "OH"], ["H2O"]), "data_class": "DERIVED_ENDPOINT_FREE_ENERGY", "kinetics": "no explicit TS in original network; endpoint assumption", "mapping": "independent states"},
    ]
    write_csv(DATA / "NEB_CONNECTION_ENERGY_MAP.csv", connection_rows)

    # Formula/invariant/state map and reports are emitted after actual CatMAP parse.
    write_csv(DATA / "FULL_8STEP_STATE_PROJECTION.csv", [
        {"state": name, "formula": json.dumps(formula, sort_keys=True), "formal_n_sites": nsites,
         "mapped_constituent_moieties": json.dumps(mapping, sort_keys=True),
         "inventory_projection_status": "EXPLICIT_BOOKKEEPING_HYPOTHESIS_NOT_STRUCTURAL_PROOF" if name.startswith("NEB") else "independent adsorbate label"}
        for name, (formula, nsites, mapping) in SPECIES.items()
    ])

    hashes = []
    for path in (SOURCE, ORIGINAL_MKM, ORIGINAL_MKM.with_name("energies_dft_neb.txt"), OLD_CANDIDATE_MKM,
                 MODEL_OUT, ENERGY_OUT):
        hashes.append({"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "size_bytes": path.stat().st_size})
    write_csv(DATA / "SOURCE_HASHES.csv", hashes)

    balance_ok = all(r["atom_balance"] == "PASS" and r["formal_site_balance"] == "PASS" for r in balance)
    if len(inv) != 2:
        raise RuntimeError(f"expected two lifted inventory invariants, got {len(inv)}")
    if len(original_inv) != 2:
        raise RuntimeError(f"expected two original inventory invariants, got {len(original_inv)}")
    if len(old_inv) != 1:
        raise RuntimeError(f"expected one old-candidate inventory invariant, got {len(old_inv)}")
    (DATA / "PARSER_VALIDATION.json").write_text(json.dumps({
        "catmap_version": getattr(catmap, "__version__", "unknown"),
        "original_network_elementary_reactions": len(original_reactions),
        "old_candidate_elementary_reactions": len(old_reactions),
        "candidate_network_elementary_reactions": len(reactions),
        "original_stable_adsorbates": len(original_names),
        "old_candidate_stable_adsorbates": len(old_names),
        "candidate_stable_adsorbates": len(names),
        "original_inventory_nullity": len(original_inv),
        "old_candidate_inventory_nullity": len(old_inv),
        "candidate_inventory_nullity": len(inv),
        "candidate_atom_balance": "PASS" if all(r["atom_balance"] == "PASS" for r in balance) else "FAIL",
        "candidate_formal_site_balance": "PASS" if all(r["formal_site_balance"] == "PASS" for r in balance) else "FAIL",
        "candidate_state_footprints": "CONDITIONAL_FRAGMENT_COUNT_NOT_STRUCTURALLY_VALIDATED",
    }, indent=2), encoding="utf-8")
    report = f"""# 完整显式 OOH 状态映射：候选模型与门控

## 状态映射

原有八步的稳定吸附态保留；完整 NEB 状态不再被压成一个一有效位点物种，而以显式复合态标签表示。当前物种式和组分投影如下，均为根据已有 NEB 状态组成建立的 bookkeeping 映射，**不等同于结构文件已经证明的局部构型或几何占位**：

- NEB1 初态/TS：Ha + OOH + C6H12；NEB1 终态：Ha + C6H12O + O + Hb。
- NEB2 初态/TS：Ha + Hb + O；NEB2 终态：Ha + OH。
- 所以 Ha 在两条路径中均作为显式保留的旁观组分；NEB1 终态生成的 Hb/O 与 NEB2 初态之间有显式状态连接，环氧产物另走吸附态脱附；NEB2 终态的 Ha/OH 再进入原水生成步骤。

为做 CatMAP parser/守恒验证，模型输入对复合态暂用 fragment-count `n_sites`：NEB1 IS/TS=3、FS=4；NEB2 IS/TS=3、FS=2。此数字是可复算的形式占位假设，不是由结构几何验证的物理位点 footprint。因此 `FULL_8STEP_MODEL.mkm` 是 **CONDITIONAL** 候选，不是正式物理模型。

## 计量和库存

CatMAP {getattr(catmap, '__version__', 'unknown')} 对原始八步快照解析 {len(original_reactions)} 个基元步骤、旧九反应候选解析 {len(old_reactions)} 步、当前完整状态候选解析 {len(reactions)} 步。候选逐反应原子守恒与形式位点守恒：{"通过" if balance_ok else "存在失败，检查 FULL_8STEP_BALANCE_CHECKS.csv"}。左零空间维数分别为原始八步 {len(original_inv)}、旧九反应 {len(old_inv)}、完整状态候选 {len(inv)}。九反应候选确实少一个库存方向；完整状态方案在代数上恢复了第二个独立不变量。矩阵见 `data/ORIGINAL_8STEP_STOICHIOMETRY.csv`、`data/OLD_NINE_REACTION_STOICHIOMETRY.csv`、`data/FULL_8STEP_STOICHIOMETRY.csv`；相应基见对应 `*_INVENTORY_BASIS.csv`。

映射约定下的库存是：

```text
CatMAP候选矩阵的一组基为 `J1 = θ(Ha*) − θ(Hb*) − θ(OH*) − θ(OOH*)`、`J2 = θ(Hb*) − θ(O*)`。故原研究写法 `I1 = θ(Ha*) − θ(OOH*) − θ(OH*) − θ(O*)` 等于 `J1 + J2`，而 `I2=J2`。对洁净表面零库存类，两者均为0。四个复合态按下列片段投影计算的 `(I1,I2)` 均为 `(0,0)`：NEB1-IS/TS = Ha+OOH+C6H12；NEB1-FS = Ha+Hb+O+C6H12O；NEB2-IS = Ha+Hb+O；NEB2-FS = Ha+OH。该映射保持库存代数，但仍是化学片段记账，不证明真实局域结构或可达性。
```

复合态的 I1/I2 系数按其代表 moieties 投影求和。显式 NEB 反应的这些量保持不变；association/deaggregation 连接也按同一投影保持不变。洁净表面初值对应两常数均为0，但这并不证明模型的稳态根能从该初态动力学可达。

## 能量与动力学边界

Ti–Co 独立态和六个 NEB IS/TS/FS 能量均直接取自隔离的历史 CatMAP 输入快照；没有修改原始 Gibbs。连接能量差由相同能量表的端点求和得到，详见 `data/NEB_CONNECTION_ENERGY_MAP.csv`。现有数值包括 NEB1 `FS→NEB2IS+C6H12O_g = +0.811237 eV`、NEB2 `FS→H2O* = −0.601109 eV`；它们是端点自由能差，不是连接活化能。

NEB1、NEB2 的 TS 保持显式，CatMAP 采用的正/逆能垒分别为 NEB1 `0.556344/0.796030 eV`、NEB2 `0.225994/0.706574 eV`。对未显式 TS 的步骤，本机源码 `catmap/solvers/mean_field_solver.py` 的 `reaction_energy_equations()` 明确令 `G_TS=max(G_IS,G_FS)`；对应正逆端点能垒列在 `data/REACTION_ENERGY_AND_BARRIERS.csv`，故其仍是官方端点近似而非 DFT TS。解析候选总反应自由能为 `−2.088914 eV`，与气相净反应 `H2O2(g)+C6H12(g)→C6H12O(g)+H2O(g)` 的输入能量差闭合。未把端点 ΔG 加到 NEB TS，也未修改 NEB 有符号端点差。

## 结论

现有能量足以建立一个原子守恒、形式位点守恒、并保留两个库存不变量的**条件性拓扑候选**；因此九反应候选不是唯一可表达形式。它还不足以证明正式物理自洽，因为复合态结构投影、Ha 的局域位置、真实 `n_sites`、复合体解离/结合的过渡态及其速率尚未由数据给定。

Ti–Co 已在该候选下尝试零库存稳态，使用 CatMAP 官方速率常数与稳态函数、受库存投影的正对数覆盖度坐标、SciPy 仅作初猜以及 CatMAP `NewtonRoot` 高精度抛光；当前全部有限初猜均未通过残差门槛，无可报告的接受 TOF/覆盖度。尝试记录位于 `runs/tico_provisional_zero_inventory/`。因此，本候选已通过解析和反应化学计量检查，但尚未通过数值稳态门控；更重要的是，复合态几何 footprint 与连接 TS 缺失，不能升级为正式模型。最少缺失数据见 `reports/MINIMUM_MISSING_DFT_DATA.md`。
"""
    (REPORTS / "FULL_8STEP_STATE_MAPPING.md").write_text(report, encoding="utf-8")

    print(json.dumps({"elementary_reactions": len(reactions), "inventory_nullity": len(inv),
                      "atom_and_formal_site_balance": balance_ok,
                      "model": str(MODEL_OUT), "matrix": str(DATA / "FULL_8STEP_STOICHIOMETRY.csv")}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
