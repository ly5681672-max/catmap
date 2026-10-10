"""Independent 333.15 K root replay and branch-audit primitives.

This module deliberately starts with root self-consistency.  It replays the
cross-precision candidates in fresh CatMAP directories and captures the full
physical state before any branch clustering or stability classification.
"""
from __future__ import annotations

import argparse
import io
import math
from contextlib import redirect_stderr, redirect_stdout
from fractions import Fraction
from pathlib import Path

import pandas as pd

import run_ooh_analysis as core


AUDIT_ROOT = core.OUTPUT / "steady_state_branch_audit"
REPLAY_ROOT = AUDIT_ROOT / "replay"
AUDIT_LOG_ROOT = AUDIT_ROOT / "audit_logs"
TRIALS_CSV = core.OUTPUT / "recovery_single_temperature" / "seed_trials.csv"
TARGET_TEMPERATURE = 333.15
TARGET_PRESSURE = 1.0
RECOMPUTE_RELATIVE_TOL = 1e-8
COVERAGE_ABSOLUTE_TOL = "1e-8"
COVERAGE_LOG10_TOL = "0.05"
LOW_COVERAGE_FLOOR = "1e-30"
RATE_RELATIVE_TOL = "1e-6"
TOF_LOG10_TOL = "0.02"
JACOBIAN_RELATIVE_ERROR_TOL = "1e-4"
NEAR_NEUTRAL_RELATIVE_TOL = "1e-8"


def _as_bool(value) -> bool:
    return str(value).strip().lower() == "true"


def read_cross_precision_candidates(path: Path = TRIALS_CSV) -> list[dict]:
    """Read the accepted cross-precision stage rows, not single-stage rows."""
    # Do not let pandas infer 1e-340 as a binary float zero.  The tolerance
    # string is part of the numerical provenance and must reach CatMAP intact.
    frame = pd.read_csv(
        path,
        dtype={
            "run_tolerance": "string",
            "initial_seed_source": "string",
        },
    )
    required = {
        "surface", "seed_id", "refinement_stage", "quality_pass",
        "refinement_stable", "run_precision", "run_tolerance",
        "initial_seed_source",
    }
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"seed_trials.csv missing required columns: {sorted(missing)}")
    selected = frame.loc[
        frame["quality_pass"].map(_as_bool)
        & frame["refinement_stable"].map(_as_bool)
    ].copy()
    selected["refinement_stage_num"] = pd.to_numeric(
        selected["refinement_stage"], errors="raise"
    )
    selected = selected.sort_values(
        ["surface", "seed_id", "refinement_stage_num"], kind="stable"
    )
    # One accepted row per seed is required.  Duplicate accepted stages would
    # otherwise be counted as independent initial states.
    selected = selected.drop_duplicates(["surface", "seed_id"], keep="last")
    return selected.drop(columns=["refinement_stage_num"]).to_dict("records")


def _root_id(candidate: dict) -> str:
    stage = int(float(candidate["refinement_stage"]))
    return f"{candidate['surface']}__{candidate['seed_id']}__stage_{stage:02d}"


def replay_candidate(candidate: dict, output_root: Path = REPLAY_ROOT) -> dict:
    """Replay one candidate in a fresh directory and capture full state."""
    root_id = _root_id(candidate)
    root = output_root / root_id
    root.mkdir(parents=True, exist_ok=True)
    log_path = AUDIT_LOG_ROOT / f"{root_id}.txt"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    seed_file = Path(str(candidate["initial_seed_source"])).resolve()
    if not seed_file.is_file():
        raise FileNotFoundError(f"Candidate seed file missing: {seed_file}")

    with log_path.open("w", encoding="utf-8") as stream:
        with redirect_stdout(stream), redirect_stderr(stream):
            result = core.run_single_surface(
                str(candidate["surface"]),
                precision=int(float(candidate["run_precision"])),
                tolerance=str(candidate["run_tolerance"]),
                output_root=root,
                seed_file=seed_file,
                max_iterations=400,
                max_bisections=5,
                capture_state=True,
            )
    result.update({
        "root_id": root_id,
        "source_surface": str(candidate["surface"]),
        "source_seed_id": str(candidate["seed_id"]),
        "source_refinement_stage": int(float(candidate["refinement_stage"])),
        "source_precision": int(float(candidate["run_precision"])),
        "source_tolerance": str(candidate["run_tolerance"]),
        "source_seed_file": str(seed_file),
        "replay_root": str(root.resolve()),
        "replay_log": str(log_path.resolve()),
        "replay_rate_self_consistent": (
            float(result.get("rate_recompute_relative_error", math.inf))
            <= RECOMPUTE_RELATIVE_TOL
        ),
        "replay_tof_self_consistent": (
            float(result.get("tof_recompute_relative_error", math.inf))
            <= RECOMPUTE_RELATIVE_TOL
        ),
        "replay_quality_pass": bool(
            result.get("quality_pass")
            and result.get("recomputed_quality_pass")
        ),
    })
    return result


def run_existing_root_replay() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Replay all existing cross-precision candidates into the audit tree."""
    candidates = read_cross_precision_candidates()
    results = []
    for candidate in candidates:
        root_id = _root_id(candidate)
        print(f"[branch-audit] replay {root_id}")
        try:
            results.append(replay_candidate(candidate))
        except Exception as exc:
            results.append({
                "root_id": root_id,
                "source_surface": str(candidate["surface"]),
                "source_seed_id": str(candidate["seed_id"]),
                "source_refinement_stage": int(float(candidate["refinement_stage"])),
                "source_precision": int(float(candidate["run_precision"])),
                "source_tolerance": str(candidate["run_tolerance"]),
                "replay_status": f"FAILED: {type(exc).__name__}: {exc}",
                "replay_quality_pass": False,
            })
    all_roots = pd.DataFrame(results)
    validated = all_roots.loc[
        all_roots.get("replay_quality_pass", pd.Series(dtype=bool)).fillna(False)
        & all_roots.get("replay_rate_self_consistent", pd.Series(dtype=bool)).fillna(False)
        & all_roots.get("replay_tof_self_consistent", pd.Series(dtype=bool)).fillna(False)
    ].copy()
    AUDIT_ROOT.mkdir(parents=True, exist_ok=True)
    all_roots.to_csv(AUDIT_ROOT / "all_root_candidates.csv", index=False)
    validated.to_csv(AUDIT_ROOT / "validated_roots.csv", index=False)
    return all_roots, validated


def retry_failed_replays() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Retry only audit rows that failed, preserving successful replays."""
    result_path = AUDIT_ROOT / "all_root_candidates.csv"
    if not result_path.is_file():
        return run_existing_root_replay()
    previous = pd.read_csv(result_path, dtype={"source_tolerance": "string"})
    failed_ids = set(
        previous.loc[
            ~previous["replay_quality_pass"].map(_as_bool), "root_id"
        ].astype(str)
    )
    candidates = read_cross_precision_candidates()
    replacements = {}
    for candidate in candidates:
        if _root_id(candidate) in failed_ids:
            print(f"[branch-audit] retry {_root_id(candidate)}")
            try:
                replacements[_root_id(candidate)] = replay_candidate(candidate)
            except Exception as exc:
                replacements[_root_id(candidate)] = {
                    "root_id": _root_id(candidate),
                    "source_surface": str(candidate["surface"]),
                    "source_seed_id": str(candidate["seed_id"]),
                    "source_refinement_stage": int(float(candidate["refinement_stage"])),
                    "source_precision": int(float(candidate["run_precision"])),
                    "source_tolerance": str(candidate["run_tolerance"]),
                    "replay_status": f"FAILED: {type(exc).__name__}: {exc}",
                    "replay_quality_pass": False,
                }
    rows = [
        replacements.get(str(row["root_id"]), row.to_dict())
        for _, row in previous.iterrows()
    ]
    all_roots = pd.DataFrame(rows)
    validated = all_roots.loc[
        all_roots["replay_quality_pass"].map(_as_bool)
        & all_roots["replay_rate_self_consistent"].map(_as_bool)
        & all_roots["replay_tof_self_consistent"].map(_as_bool)
    ].copy()
    all_roots.to_csv(result_path, index=False)
    validated.to_csv(AUDIT_ROOT / "validated_roots.csv", index=False)
    return all_roots, validated


def _mp_vector(text: str) -> list:
    from mpmath import mp

    if not isinstance(text, str) or not text:
        raise ValueError("Missing high-precision vector")
    return [mp.mpf(value) for value in text.split(";")]


def coverage_distance_mp(left: dict, right: dict) -> dict:
    """Compare full physical coverages with an explicit low-value floor."""
    from mpmath import mp

    a = _mp_vector(left["full_coverage_vector_high_precision"])
    b = _mp_vector(right["full_coverage_vector_high_precision"])
    if left["coverage_order"] != right["coverage_order"] or len(a) != len(b):
        raise ValueError("Coverage order/length mismatch")
    abs_max = max(abs(x - y) for x, y in zip(a, b))
    floor = mp.mpf(LOW_COVERAGE_FLOOR)
    log_diffs = [
        abs(mp.log10(x) - mp.log10(y))
        for x, y in zip(a, b)
        if x > floor and y > floor
    ]
    log_max = max(log_diffs) if log_diffs else mp.mpf("0")
    return {
        "coverage_abs_max": abs_max,
        "coverage_log10_max_above_floor": log_max,
        "coverage_same": (
            abs_max <= mp.mpf(COVERAGE_ABSOLUTE_TOL)
            and log_max <= mp.mpf(COVERAGE_LOG10_TOL)
        ),
    }


def rate_distance_mp(left: dict, right: dict) -> dict:
    """Compare all elementary rates and epoxide TOF in high precision."""
    from mpmath import mp

    a = _mp_vector(left["net_rate_vector_high_precision"])
    b = _mp_vector(right["net_rate_vector_high_precision"])
    if len(a) != len(b):
        raise ValueError("Net-rate vector length mismatch")
    scale = max(max(abs(x) for x in a), max(abs(y) for y in b))
    abs_max = max(abs(x - y) for x, y in zip(a, b))
    relative = abs_max / scale if scale else mp.inf
    tof_delta = abs(
        mp.mpf(left["log10_net_C6H12O"])
        - mp.mpf(right["log10_net_C6H12O"])
    )
    return {
        "rate_relative_max": relative,
        "tof_log10_delta": tof_delta,
        "rate_same": (
            relative <= mp.mpf(RATE_RELATIVE_TOL)
            and tof_delta <= mp.mpf(TOF_LOG10_TOL)
        ),
    }


def _same_branch(left: dict, right: dict) -> tuple[bool, dict]:
    coverage = coverage_distance_mp(left, right)
    rates = rate_distance_mp(left, right)
    same = bool(coverage["coverage_same"] and rates["rate_same"])
    return same, {**coverage, **rates}


def cluster_roots(records: list[dict]) -> list[dict]:
    """Cluster validated roots deterministically within each surface."""
    grouped = {}
    for record in records:
        grouped.setdefault(str(record["source_surface"]), []).append(record)
    output = []
    for surface in sorted(grouped):
        rows = sorted(grouped[surface], key=lambda row: str(row["root_id"]))
        parent = list(range(len(rows)))

        def find(index):
            while parent[index] != index:
                parent[index] = parent[parent[index]]
                index = parent[index]
            return index

        def union(left, right):
            left, right = find(left), find(right)
            if left != right:
                parent[right] = left

        for i in range(len(rows)):
            for j in range(i + 1, len(rows)):
                same, _ = _same_branch(rows[i], rows[j])
                if same:
                    union(i, j)
        branch_roots = {}
        for index, row in enumerate(rows):
            branch_roots.setdefault(find(index), []).append(row)
        ordered = sorted(
            branch_roots.values(),
            key=lambda members: min(str(row["root_id"]) for row in members),
        )
        branch_ids = {
            id(members): f"{surface}_branch_{index:02d}"
            for index, members in enumerate(ordered, 1)
        }
        for members in ordered:
            branch_id = branch_ids[id(members)]
            for row in members:
                output.append({
                    "surface": surface,
                    "branch_id": branch_id,
                    "root_id": row["root_id"],
                    "seed_id": row["source_seed_id"],
                    "full_coverage_vector": row[
                        "full_coverage_vector_high_precision"
                    ],
                    "full_numbers_vector": row[
                        "full_numbers_vector_high_precision"
                    ],
                    "net_rate_vector": row["net_rate_vector_high_precision"],
                    "gas_tof_vector": row["gas_tof_vector_high_precision"],
                    "dominant_species": row["dominant_species"],
                    "log10_TOF": row["log10_net_C6H12O"],
                    "precision": row["run_precision"],
                    "tolerance": row["run_tolerance"],
                    "quality_status": row["quality_status"],
                    "replay_quality_pass": row["replay_quality_pass"],
                })
    return output


def cluster_existing_validated_roots() -> pd.DataFrame:
    validated_path = AUDIT_ROOT / "validated_roots.csv"
    if not validated_path.is_file():
        raise FileNotFoundError("先完成根重放并生成validated_roots.csv")
    frame = pd.read_csv(validated_path, dtype=str)
    records = frame.to_dict("records")
    clustered = pd.DataFrame(cluster_roots(records))
    clustered.to_csv(AUDIT_ROOT / "root_clusters.csv", index=False)
    return clustered


def _clear_solver_memos(solver) -> None:
    for memo_name in ("_rate_constant_memoize", "_steady_state_memoize"):
        memo = getattr(solver, memo_name, None)
        if hasattr(memo, "clear"):
            memo.clear()


def _matrix_max_abs(matrix) -> object:
    return max(abs(matrix[i, j]) for i in range(matrix.rows) for j in range(matrix.cols))


def _row_sum_norm(matrix) -> object:
    return max(
        sum(abs(matrix[i, j]) for j in range(matrix.cols))
        for i in range(matrix.rows)
    )


def _exact_rank(rows: list[list[int]]) -> int:
    """Return the exact rank of a small integer matrix."""
    if not rows or not rows[0]:
        return 0
    work = [[Fraction(value) for value in row] for row in rows]
    n_rows = len(work)
    n_cols = len(work[0])
    rank = 0
    for col in range(n_cols):
        pivot = next(
            (row for row in range(rank, n_rows) if work[row][col]), None
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        pivot_value = work[rank][col]
        work[rank] = [value / pivot_value for value in work[rank]]
        for row in range(n_rows):
            if row == rank or not work[row][col]:
                continue
            factor = work[row][col]
            work[row] = [
                value - factor * pivot_value2
                for value, pivot_value2 in zip(work[row], work[rank])
            ]
        rank += 1
        if rank == n_rows:
            break
    return rank


def _independent_column_indices(rows: list[list[int]]) -> list[int]:
    """Select an exact independent column basis from an integer matrix."""
    if not rows or not rows[0]:
        return []
    selected: list[int] = []
    rank = 0
    for column in range(len(rows[0])):
        candidate = [
            [row[index] for index in selected + [column]]
            for row in rows
        ]
        candidate_rank = _exact_rank(candidate)
        if candidate_rank > rank:
            selected.append(column)
            rank = candidate_rank
    return selected


def _surface_stoichiometric_rows(elementary_rxns, adsorbate_names) -> list[list[int]]:
    """Build d(theta_ads)/dt stoichiometry in solver adsorbate order."""
    from collections import Counter

    rows = []
    for species in adsorbate_names:
        row = []
        for reaction in elementary_rxns:
            initial = Counter(reaction[0])
            final = Counter(reaction[-1])
            row.append(final.get(species, 0) - initial.get(species, 0))
        rows.append(row)
    return rows


def _compatibility_class_jacobian(reduced, stoichiometric_rows):
    """Represent J on the stoichiometric compatibility tangent space."""
    from mpmath import mp

    basis_columns = _independent_column_indices(stoichiometric_rows)
    if not basis_columns:
        raise ValueError("Surface stoichiometric matrix has zero rank")
    basis = mp.matrix(len(stoichiometric_rows), len(basis_columns))
    for i, row in enumerate(stoichiometric_rows):
        for j, column in enumerate(basis_columns):
            basis[i, j] = mp.mpf(row[column])
    # L is a left inverse of B.  Therefore L * J * B is the Jacobian
    # represented in independent stoichiometric coordinates.
    left_inverse = (basis.T * basis) ** -1 * basis.T
    return left_inverse * reduced * basis, basis_columns


def _classify_linear_stability(eigenvalues, precision: int, derivative_error):
    """Classify sign after treating precision-scale zeros as neutral."""
    from mpmath import mp

    if derivative_error > mp.mpf(JACOBIAN_RELATIVE_ERROR_TOL):
        return "near_neutral_or_undetermined"
    if not eigenvalues:
        return "near_neutral_or_undetermined"
    # Using the largest eigenvalue as an absolute threshold would hide slow
    # but genuine modes behind fast adsorption modes.  The threshold is tied
    # to the working precision instead; exact conservation zeros are many
    # orders smaller than physical slow modes in this model.
    near_threshold = mp.power(10, -max(20, int(precision) // 2))
    real_parts = [mp.re(value) for value in eigenvalues]
    if any(value > near_threshold for value in real_parts):
        return "linearly_unstable"
    if all(value < -near_threshold for value in real_parts):
        return "linearly_stable"
    return "near_neutral_or_undetermined"


def _classify_full_simplex_stability(
    full_status: str, compatibility_status: str, extra_conserved_modes: int
) -> str:
    """Keep physical invariant directions in the full-simplex conclusion."""
    if extra_conserved_modes < 0:
        raise ValueError("Negative conserved-mode count")
    return compatibility_status if extra_conserved_modes == 0 else full_status


def _reduced_physical_jacobian(solver, theta_full, rxn_parameters):
    """Return CatMAP's coverage Jacobian with the vacancy coordinate removed."""
    from mpmath import mp

    adsorbates = list(solver.adsorbate_names)
    site_names = [name for name in solver.site_names if name != "g"]
    if len(site_names) != 1 or len(theta_full) != len(adsorbates) + 1:
        raise ValueError("Expected one site and adsorbates plus one vacancy coverage")
    _clear_solver_memos(solver)
    solver.get_rate(rxn_parameters, coverages=theta_full, verify_coverages=False)
    function_values = list(solver.ideal_steady_state_function(theta_full))
    full_jacobian = solver.ideal_steady_state_jacobian(theta_full)
    n_ads = len(adsorbates)
    if full_jacobian.rows != n_ads + 1 or full_jacobian.cols != n_ads + 1:
        raise ValueError(
            f"Unexpected CatMAP coverage Jacobian shape: "
            f"{full_jacobian.rows}x{full_jacobian.cols}"
        )
    vacancy_index = n_ads
    reduced = mp.matrix(n_ads, n_ads)
    for i in range(n_ads):
        for j in range(n_ads):
            # theta_vacant = 1 - sum(theta_adsorbates)
            reduced[i, j] = full_jacobian[i, j] - full_jacobian[i, vacancy_index]
    return function_values, full_jacobian, reduced


def _finite_difference_reduced_jacobian(solver, theta_full, step_hint, min_step):
    """Finite-difference df_ads/dtheta_ads on the constrained simplex."""
    from mpmath import mp

    n_ads = len(solver.adsorbate_names)
    base = [mp.mpf(value) for value in theta_full]
    vacancy = base[-1]
    if len(base) != n_ads + 1:
        raise ValueError("Coverage vector does not contain the expected vacancy")
    reduced = mp.matrix(n_ads, n_ads)
    used_steps = []
    for j in range(n_ads):
        h = min(
            mp.mpf(step_hint),
            base[j] / 4 if base[j] > 0 else mp.mpf("0"),
            vacancy / 4 if vacancy > 0 else mp.mpf("0"),
        )
        if h <= mp.mpf(min_step):
            raise ArithmeticError(
                f"feasible finite-difference step too small for coordinate {j}: {h}"
            )
        plus = list(base)
        minus = list(base)
        plus[j] += h
        plus[-1] -= h
        minus[j] -= h
        minus[-1] += h
        f_plus = list(solver.ideal_steady_state_function(plus))[:n_ads]
        f_minus = list(solver.ideal_steady_state_function(minus))[:n_ads]
        for i in range(n_ads):
            reduced[i, j] = (f_plus[i] - f_minus[i]) / (2 * h)
        used_steps.append(h)
    return reduced, used_steps


def analyze_branch_stability(root: dict) -> dict:
    """Analyze one validated root in physical coverage coordinates."""
    from mpmath import mp
    from catmap import ReactionModel

    surface = str(root["source_surface"])
    replay_root = Path(str(root["replay_root"]))
    config = replay_root / surface / f"baseline_{surface}.mkm"
    if not config.is_file():
        raise FileNotFoundError(f"Replay setup missing: {config}")
    model = ReactionModel(setup_file=str(config))
    # ReactionModel.__init__ only parses the setup file.  CatMAP installs
    # _mpfloat/_matrix and compiles the generated rate/Jacobian functions in
    # ReactionModel.run().  Re-run the isolated one-point replay from its
    # preserved candidate seed so the stability audit uses the official
    # initialization path without touching the source model or DFT inputs.
    captured = io.StringIO()
    with redirect_stdout(captured), redirect_stderr(captured):
        model.run(recalculate=True)
    model._descriptors = [TARGET_TEMPERATURE, TARGET_PRESSURE]
    model.solver._descriptors = [TARGET_TEMPERATURE, TARGET_PRESSURE]
    model.solver.compile()
    theta = _mp_vector(root["full_coverage_vector_high_precision"])
    params = model.scaler.get_rxn_parameters([TARGET_TEMPERATURE, TARGET_PRESSURE])
    function_values, full_jacobian, reduced = _reduced_physical_jacobian(
        model.solver, theta, params
    )
    precision = int(model.decimal_precision)
    finite, steps = _finite_difference_reduced_jacobian(
        model.solver,
        theta,
        step_hint="1e-8",
        min_step=mp.power(10, -(precision - 20)),
    )
    difference = max(
        abs(reduced[i, j] - finite[i, j])
        for i in range(reduced.rows) for j in range(reduced.cols)
    )
    scale = max(_matrix_max_abs(reduced), _matrix_max_abs(finite))
    derivative_error = difference / scale if scale else mp.inf
    try:
        eigenvalues = list(mp.eig(reduced, left=False, right=False))
    except Exception as exc:
        return {
            "root_id": root["root_id"],
            "surface": surface,
            "stability_status": "jacobian_evaluation_failed",
            "jacobian_error": f"eigenvalue_failure: {type(exc).__name__}: {exc}",
        }
    precision = int(model.decimal_precision)
    full_reduced_stability_status = _classify_linear_stability(
        eigenvalues, precision, derivative_error
    )

    stoichiometric_rows = _surface_stoichiometric_rows(
        model.elementary_rxns, list(model.solver.adsorbate_names)
    )
    compatibility, basis_columns = _compatibility_class_jacobian(
        reduced, stoichiometric_rows
    )
    compatibility_finite, _ = _compatibility_class_jacobian(
        finite, stoichiometric_rows
    )
    compatibility_difference = max(
        abs(compatibility[i, j] - compatibility_finite[i, j])
        for i in range(compatibility.rows)
        for j in range(compatibility.cols)
    )
    compatibility_scale = max(
        _matrix_max_abs(compatibility), _matrix_max_abs(compatibility_finite)
    )
    compatibility_derivative_error = (
        compatibility_difference / compatibility_scale
        if compatibility_scale else mp.inf
    )
    try:
        compatibility_eigenvalues = list(
            mp.eig(compatibility, left=False, right=False)
        )
    except Exception as exc:
        return {
            "root_id": root["root_id"],
            "surface": surface,
            "stability_status": "jacobian_evaluation_failed",
            "jacobian_error": (
                f"compatibility_eigenvalue_failure: {type(exc).__name__}: {exc}"
            ),
        }
    compatibility_stability_status = _classify_linear_stability(
        compatibility_eigenvalues, precision, compatibility_derivative_error
    )
    # Stoichiometric invariants constrain real trajectories but do not remove
    # physical coverage coordinates.  A root can be attractive inside its
    # invariant class while remaining neutral in the full coverage simplex.
    stability_status = _classify_full_simplex_stability(
        full_reduced_stability_status,
        compatibility_stability_status,
        reduced.rows - len(basis_columns),
    )
    compatibility_real_parts = [mp.re(value) for value in compatibility_eigenvalues]
    try:
        condition_number = _row_sum_norm(reduced) * _row_sum_norm(reduced ** -1)
    except Exception:
        condition_number = mp.inf
    try:
        compatibility_condition_number = (
            _row_sum_norm(compatibility)
            * _row_sum_norm(compatibility ** -1)
        )
    except Exception:
        compatibility_condition_number = mp.inf
    return {
        "root_id": root["root_id"],
        "surface": surface,
        "branch_id": root.get("branch_id", ""),
        "dominant_species": root["dominant_species"],
        "log10_TOF": root["log10_net_C6H12O"],
        "jacobian_dimension": reduced.rows,
        "jacobian_precision": precision,
        "jacobian_tolerance": str(root.get("run_tolerance", "")),
        "full_function_max_abs": mp.nstr(max(abs(value) for value in function_values), 80),
        "full_jacobian_has_redundant_vacancy": True,
        "eigenvalues_high_precision": ";".join(mp.nstr(value, 80) for value in eigenvalues),
        "max_real_part_high_precision": mp.nstr(max(mp.re(value) for value in eigenvalues), 80),
        "jacobian_condition_number_high_precision": mp.nstr(condition_number, 80),
        "finite_difference_steps_high_precision": ";".join(mp.nstr(value, 40) for value in steps),
        "derivative_check_relative_error_high_precision": mp.nstr(derivative_error, 80),
        "derivative_check_pass": derivative_error <= mp.mpf(JACOBIAN_RELATIVE_ERROR_TOL),
        "full_reduced_stability_status": full_reduced_stability_status,
        "stoichiometric_rank": len(basis_columns),
        "structural_zero_mode_count": reduced.rows - len(basis_columns),
        "stoichiometric_basis_column_indices": ";".join(str(i) for i in basis_columns),
        "compatibility_jacobian_dimension": compatibility.rows,
        "compatibility_eigenvalues_high_precision": ";".join(
            mp.nstr(value, 80) for value in compatibility_eigenvalues
        ),
        "compatibility_max_real_part_high_precision": mp.nstr(
            max(compatibility_real_parts), 80
        ),
        "compatibility_condition_number_high_precision": mp.nstr(
            compatibility_condition_number, 80
        ),
        "compatibility_derivative_check_relative_error_high_precision": mp.nstr(
            compatibility_derivative_error, 80
        ),
        "compatibility_derivative_check_pass": (
            compatibility_derivative_error
            <= mp.mpf(JACOBIAN_RELATIVE_ERROR_TOL)
        ),
        "compatibility_stability_status": compatibility_stability_status,
        "stability_status": stability_status,
        "jacobian_error": "",
    }


def analyze_existing_branch_stability() -> pd.DataFrame:
    validated = pd.read_csv(AUDIT_ROOT / "validated_roots.csv", dtype=str)
    clusters = pd.read_csv(AUDIT_ROOT / "root_clusters.csv", dtype=str)
    validated = validated.merge(
        clusters[["root_id", "branch_id"]], on="root_id", how="left"
    )
    rows = []
    for root in validated.to_dict("records"):
        print(f"[branch-audit] jacobian {root['root_id']}")
        try:
            rows.append(analyze_branch_stability(root))
        except Exception as exc:
            rows.append({
                "root_id": root["root_id"],
                "surface": root["source_surface"],
                "branch_id": root.get("branch_id", ""),
                "stability_status": "jacobian_evaluation_failed",
                "jacobian_error": f"{type(exc).__name__}: {exc}",
            })
    result = pd.DataFrame(rows)
    result.to_csv(AUDIT_ROOT / "jacobian_eigenvalues.csv", index=False)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--replay-existing",
        action="store_true",
        help="replay the 333.15 K cross-precision candidates",
    )
    parser.add_argument(
        "--retry-failed",
        action="store_true",
        help="retry only failed rows in an existing audit CSV",
    )
    parser.add_argument(
        "--cluster-existing",
        action="store_true",
        help="cluster validated roots using full coverage and rate vectors",
    )
    parser.add_argument(
        "--analyze-jacobian",
        action="store_true",
        help="analyze validated roots in reduced physical coverage space",
    )
    args = parser.parse_args()
    if args.analyze_jacobian:
        result = analyze_existing_branch_stability()
        print(f"jacobian_roots={len(result)}")
        print(result["stability_status"].value_counts().to_string())
        return
    if args.cluster_existing:
        clustered = cluster_existing_validated_roots()
        print(f"clustered_roots={len(clustered)}")
        print(f"branches={clustered[['surface', 'branch_id']].drop_duplicates().shape[0]}")
        return
    if args.retry_failed:
        all_roots, validated = retry_failed_replays()
    elif args.replay_existing:
        all_roots, validated = run_existing_root_replay()
    else:
        parser.error("先明确指定 --replay-existing、--retry-failed、--cluster-existing 或 --analyze-jacobian")
    print(f"all_root_candidates={len(all_roots)}")
    print(f"validated_roots={len(validated)}")


if __name__ == "__main__":
    main()
