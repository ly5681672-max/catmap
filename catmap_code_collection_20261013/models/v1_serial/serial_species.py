"""Evidence-aware species and energy-reference helpers for serial NEB V1."""

from __future__ import annotations

from typing import Any


STABLE_SURFACE_SPECIES = (
    "Rpair_s",
    "EpoxideHHO_s",
    "HHO_s",
    "OHHa_s",
    "H2O_s",
)

TRANSITION_STATE_SPECIES = ("TS1_s", "TS2_s")

_DEFINITIONS: dict[str, dict[str, Any]] = {
    "Rpair_s": {
        "composition": {"C": 6, "H": 14, "O": 2},
        "role": "stable_adsorbate",
        "physical_interpretation": "C6H12 plus activated H2O2-derived OOH and Ha in one complete coadsorbed state",
    },
    "TS1_s": {
        "composition": {"C": 6, "H": 14, "O": 2},
        "role": "transition_state",
        "physical_interpretation": "NEB1 oxygen-transfer transition state",
    },
    "EpoxideHHO_s": {
        "composition": {"C": 6, "H": 14, "O": 2},
        "role": "stable_adsorbate",
        "physical_interpretation": "adsorbed C6H12O plus surface Ha, Hb, O; no molecular H2O",
    },
    "HHO_s": {
        "composition": {"H": 2, "O": 1},
        "role": "stable_adsorbate",
        "physical_interpretation": "surface Ha, Hb, O fragments after epoxide desorption; not molecular H2O",
    },
    "TS2_s": {
        "composition": {"H": 2, "O": 1},
        "role": "transition_state",
        "physical_interpretation": "NEB2 surface Hb-O hydrogen-transfer transition state",
    },
    "OHHa_s": {
        "composition": {"H": 2, "O": 1},
        "role": "stable_adsorbate",
        "physical_interpretation": "surface OH plus retained Ha",
    },
    "H2O_s": {
        "composition": {"H": 2, "O": 1},
        "role": "stable_adsorbate",
        "physical_interpretation": "molecular water adsorbed on the surface; workbook post-NEB2 state",
    },
}


def species_definitions() -> dict[str, dict[str, Any]]:
    """Return explicit atom counts and a clearly provisional one-center occupancy.

    The user's confirmed element inventory is treated as structure evidence. The
    single effective-center ``n_sites`` is only a coarse-grained parsing assumption;
    it is not inferred from unavailable coordinates.
    """
    result: dict[str, dict[str, Any]] = {}
    for name, definition in _DEFINITIONS.items():
        result[name] = {
            **definition,
            "n_sites": 1,
            "n_sites_evidence": "ASSUMED_SINGLE_EFFECTIVE_CENTER",
        }
    result.update(
        {
            "H2O2_g": {"composition": {"H": 2, "O": 2}, "role": "gas_reservoir"},
            "C6H12_g": {"composition": {"C": 6, "H": 12}, "role": "gas_reservoir"},
            "C6H12O_g": {"composition": {"C": 6, "H": 12, "O": 1}, "role": "gas_reservoir"},
            "H2O_g": {"composition": {"H": 2, "O": 1}, "role": "gas_reservoir"},
        }
    )
    return result


def convert_state_energy(
    *,
    raw_g_eV: float,
    source_cell: str,
    source_formula: str,
    common_reference_eV: float | None,
    spectator_adjustments_eV: list[float],
    spectator_notes: list[str] | None = None,
) -> dict[str, Any]:
    """Apply explicit spectator terms and a stated common reference, if known.

    Spectator adjustments are signed additions to the source state G. A formation
    energy is intentionally unavailable when the common reference is not supplied.
    """
    adjustment = sum(spectator_adjustments_eV)
    adjusted_g = float(raw_g_eV) + adjustment
    resolved = common_reference_eV is not None
    return {
        "source_cell": source_cell,
        "source_formula": source_formula,
        "raw_G_eV": float(raw_g_eV),
        "spectator_adjustments_eV": list(spectator_adjustments_eV),
        "spectator_adjustment_sum_eV": adjustment,
        "spectator_notes": list(spectator_notes or []),
        "adjusted_G_eV": adjusted_g,
        "common_reference_eV": common_reference_eV,
        "formation_energy_eV": adjusted_g - common_reference_eV if resolved else None,
        "conversion_status": "RESOLVED" if resolved else "UNRESOLVED_REFERENCE",
        "reason": "" if resolved else "common slab/reference energy is not established from the available source records",
    }
