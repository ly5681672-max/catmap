"""Static CatMAP and thermodynamic-evaluability audit; no kinetics are run."""

from __future__ import annotations

import csv
import hashlib
import json
import runpy
import sys
from pathlib import Path


V2_ROOT = Path(__file__).resolve().parents[1]
ROOT = V2_ROOT.parent
V1_DATA = ROOT / "data"
V1_MODEL_DIR = ROOT / "models" / "v1_serial"
DATA = V2_ROOT / "data"
REPORTS = V2_ROOT / "reports"
MODEL_FILE = V1_MODEL_DIR / "ti_re_serial_neb.mkm"
EXPECTED_RXNS = [
    "* + H2O2_g + C6H12_g <-> Rpair*",
    "Rpair* <-> TS1* -> EpoxideHHO*",
    "EpoxideHHO* <-> HHO* + C6H12O_g",
    "HHO* <-> TS2* -> OHHa*",
    "OHHa* <-> H2O*",
    "H2O* <-> H2O_g + *",
]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest().upper()


def static_catmap_validation() -> dict:
    sys.path.insert(0, str(V1_MODEL_DIR))
    import validate_serial_network

    config = runpy.run_path(str(MODEL_FILE))
    if config["rxn_expressions"] != EXPECTED_RXNS:
        raise AssertionError("Frozen V1 six-step reaction network differs from the approved expression list")
    v1_manifest = json.loads((V1_DATA / "serial_v1_sha256.json").read_text(encoding="utf-8"))
    expected_hash = v1_manifest["output_sha256"].get(str(MODEL_FILE))
    observed_hash = sha256(MODEL_FILE)
    if not expected_hash or expected_hash != observed_hash:
        raise AssertionError("V1 model hash no longer matches the frozen V1 audit")
    model = validate_serial_network.build_model()
    model.verify()
    version = __import__("catmap").__version__
    if version != "0.3.1":
        raise AssertionError(f"Expected CatMAP 0.3.1, observed {version}")
    return {
        "catmap_version": version,
        "reaction_count": len(model.elementary_rxns),
        "adsorbate_count": len(model.adsorbate_names),
        "ts_count": len(model.transition_state_names),
        "verify_status": "PASS_STATIC_ELEMENT_AND_SITE_BALANCE_UNDER_V1_ASSUMPTION",
        "model_sha256_matches_v1": True,
        "formation_energies_present": False,
        "rate_or_steady_solver_called": False,
        "tof_calculated": False,
    }


def build_consistency_rows() -> list[dict[str, str]]:
    tire1 = next(row for row in read_csv(V1_DATA / "neb_to_catmap_energy_mapping.csv")
                 if row["catalyst"] == "tire" and row["path"] == "NEB1")
    tire2 = next(row for row in read_csv(V1_DATA / "neb_to_catmap_energy_mapping.csv")
                 if row["catalyst"] == "tire" and row["path"] == "NEB2")
    surface = next(row for row in read_csv(V1_DATA / "SERIAL_14_SURFACE_MAPPING.csv") if row["catalyst"] == "tire")
    steps = [
        ("step_1_feed_capture", "* + H2O2_g + C6H12_g <-> Rpair*", "", "MISSING_FEED_CAPTURE_FREE_ENERGY_AND_REFERENCE"),
        ("step_2_NEB1", "Rpair* <-> TS1* -> EpoxideHHO*", tire1["raw_FS_minus_IS"], "NEB1_LOCAL_RELATION_REPLAYED;GLOBAL_REFERENCE_NOT_ESTABLISHED"),
        ("step_3_epoxide_desorption", "EpoxideHHO* <-> HHO* + C6H12O_g", surface["epoxide_desorption_delta_G_eV"], "CONDITIONAL_ON_COMMON_SLAB_AND_SPECTATOR_MAPPING"),
        ("step_4_NEB2", "HHO* <-> TS2* -> OHHa*", tire2["raw_FS_minus_IS"], "NEB2_LOCAL_RELATION_REPLAYED;IS_TS_TEMPERATURE_MISSING"),
        ("step_5_water_formation", "OHHa* <-> H2O*", surface["water_formation_delta_G_eV"], "CONDITIONAL_ENDPOINT_ONLY;NO_TS"),
        ("step_6_water_desorption", "H2O* <-> H2O_g + *", surface["water_desorption_delta_G_eV"], "CONDITIONAL_ENDPOINT_ONLY;NO_TS"),
    ]
    return [{
        "catalyst": "tire", "step": name, "reaction": expression,
        "local_or_conditional_delta_G_eV": value,
        "global_reaction_delta_G_eV": "",
        "source_energy_evidence": evidence,
        "global_thermodynamic_status": "NOT_EVALUABLE",
        "reason": "At least one required global formation energy is blank; no common clean-slab reference is established.",
    } for name, expression, value, evidence in steps]


def build_surface_rows() -> list[dict[str, str]]:
    base = read_csv(V1_DATA / "SERIAL_14_SURFACE_MAPPING.csv")
    inventory = read_csv(V1_DATA / "original_neb_state_inventory.csv")
    temperatures: dict[str, dict[str, str]] = {}
    for row in inventory:
        if row["NEB_role"] in {"NEB1_IS", "NEB1_TS", "NEB1_FS", "NEB2_IS", "NEB2_TS", "NEB2_FS"}:
            temperatures.setdefault(row["catalyst"], {})[row["NEB_role"]] = row["temperature_K"]
    rows = []
    for item in base:
        cat = item["catalyst"]
        temp = temperatures[cat]
        neb1_status = "ALL_333.15_K_FIELDS_RECORDED" if all(temp.get(f"NEB1_{state}") == "333.15" for state in ("IS", "TS", "FS")) else "NEB1_TEMPERATURE_INCOMPLETE"
        neb2_status = "ALL_333.15_K_FIELDS_RECORDED" if all(temp.get(f"NEB2_{state}") == "333.15" for state in ("IS", "TS", "FS")) else "NEB2_TEMPERATURE_INCOMPLETE"
        rows.append({
            "catalyst": cat,
            "NEB1_IS_G_eV": item["NEB1_IS_G_eV"], "NEB1_TS_G_eV": item["NEB1_TS_G_eV"], "NEB1_FS_G_eV": item["NEB1_FS_G_eV"],
            "NEB1_raw_TS_minus_IS_eV": item["NEB1_raw_TS_minus_IS_eV"], "NEB1_raw_TS_minus_FS_eV": item["NEB1_raw_TS_minus_FS_eV"],
            "NEB1_CatMAP_Ea_forward_eV": item["NEB1_CatMAP_Ea_forward_eV"], "NEB1_CatMAP_Ea_reverse_eV": item["NEB1_CatMAP_Ea_reverse_eV"],
            "NEB2_IS_G_eV": item["NEB2_IS_G_eV"], "NEB2_TS_G_eV": item["NEB2_TS_G_eV"], "NEB2_FS_G_eV": item["NEB2_FS_G_eV"],
            "NEB2_raw_TS_minus_IS_eV": item["NEB2_raw_TS_minus_IS_eV"], "NEB2_raw_TS_minus_FS_eV": item["NEB2_raw_TS_minus_FS_eV"],
            "NEB2_CatMAP_Ea_forward_eV": item["NEB2_CatMAP_Ea_forward_eV"], "NEB2_CatMAP_Ea_reverse_eV": item["NEB2_CatMAP_Ea_reverse_eV"],
            "NEB1_temperature_status": neb1_status, "NEB2_temperature_status": neb2_status,
            "epoxide_desorption_conditional_delta_G_eV": item["epoxide_desorption_delta_G_eV"],
            "epoxide_desorption_status": "CONDITIONAL_NOT_ACTIVATION_BARRIER",
            "water_formation_conditional_delta_G_eV": item["water_formation_delta_G_eV"],
            "water_desorption_conditional_delta_G_eV": item["water_desorption_delta_G_eV"],
            "clean_slab_reference": "", "global_formation_energy_eV": "",
            "global_reference_status": "UNRESOLVED", "whole_cycle_status": "NOT_EVALUABLE",
            "missing_evidence": "clean slab and frequencies; atom-resolved compositions; cross-NEB reference; missing NEB2 IS/TS temperature metadata",
        })
    return rows


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def build() -> dict:
    DATA.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)
    static = static_catmap_validation()
    consistency = build_consistency_rows()
    surfaces = build_surface_rows()
    write_csv(DATA / "V2_THERMODYNAMIC_CONSISTENCY.csv", consistency)
    write_csv(DATA / "V2_14_SURFACE_GLOBAL_REFERENCE.csv", surfaces)
    consistency_report = [
        "# V2六步热力学一致性与CatMAP静态验证", "",
        "## 静态CatMAP审计", "",
        f"CatMAP {static['catmap_version']}解析得到{static['reaction_count']}步、{static['adsorbate_count']}个稳定吸附物种及{static['ts_count']}个显式TS；`verify()`结果：{static['verify_status']}。V1模型SHA256与冻结清单一致。该检查验证的是V1单有效中心粗粒化假设下的语法/元素/位点账面，不验证真实几何占位。", "",
        "未调用稳态、覆盖度、速率、TOF或动力学接口。无全局形成能输入，因此CatMAP反应能和全循环热力学闭合均不能评估。",
        "", "## 全循环与详细平衡", "",
        "六步中仅NEB1/NEB2内部有符号局部能差有源表支持；脱附及水形成/脱附端点是条件值，前置进料捕获反应仍缺能量。跨NEB统一能量零点未建立，因此六步反应能总和与气相净反应化学势的比较、全局详细平衡检查均为`NOT_EVALUABLE`，不是PASS。",
        "", "逐步状态见`../data/V2_THERMODYNAMIC_CONSISTENCY.csv`。NEB1与NEB2内部相对能差分别保留；NEB2 IS/TS温度字段缺失，故不得将全部NEB2热化学量称为已证实的333.15 K值。",
    ]
    (REPORTS / "V2_CATMAP_STATIC_VALIDATION.md").write_text("\n".join(consistency_report) + "\n", encoding="utf-8")
    surface_report = [
        "# V2 14位点全局参考状态", "",
        "所有位点均保留NEB1/NEB2源表能量和有符号差值。由于clean-slab、结构映射和共同参考证据缺失，14个体系的全局formation_energy均为空，参考状态为UNRESOLVED；不得使用本表跨位点拟合形成能或推导TOF。", "",
        "| 位点 | NEB1温度字段 | NEB2温度字段 | NEB1 TS−IS | NEB2 TS−IS | 脱附条件ΔG/eV | 全局形成能 | 状态 |", "|---|---|---|---:|---:|---:|---|---|",
    ]
    for row in surfaces:
        surface_report.append(f"| {row['catalyst']} | {row['NEB1_temperature_status']} | {row['NEB2_temperature_status']} | {row['NEB1_raw_TS_minus_IS_eV']} | {row['NEB2_raw_TS_minus_IS_eV']} | {row['epoxide_desorption_conditional_delta_G_eV']} |  | UNRESOLVED |")
    surface_report.extend([
        "", "NEB2原始负TS−IS值保留在CSV中；CatMAP有效Ea单独列出，不裁剪源数据。脱附/水端点值均为条件热力学比较，不是活化能。当前不具备将Ti–Re规则推广为VERIFIED全局形成能映射的依据；其它位点仅完成只读conditional/unresolved状态整理。",
    ])
    (REPORTS / "V2_14_SURFACE_REFERENCE_STATUS.md").write_text("\n".join(surface_report) + "\n", encoding="utf-8")
    return {"catmap": static, "thermodynamic_steps": len(consistency), "surface_rows": len(surfaces),
            "whole_cycle_status": "NOT_EVALUABLE", "all_global_formation_energies_blank": all(not row["global_formation_energy_eV"] for row in surfaces)}


if __name__ == "__main__":
    print(json.dumps(build(), ensure_ascii=False, indent=2))
