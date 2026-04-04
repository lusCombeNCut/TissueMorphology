#!/usr/bin/env python3
"""
Convert CryptBudding JSON config files between physical (SI) and
non-dimensional (simulation) units.

Conversion factors derived from CryptBuddingParams.hpp ComputeConversionFactors():
    L0 = cellDiameterMicrometers * 1e-6        [m]
    T0 = 3600                                   [s]
    eta_phys = (referenceForcenN * 1e-9) / (referenceVelocityUmPerH * 1e-6 / T0)  [N·s/m]
    kF = eta_phys / T0                          stiffness factor [N/m per sim unit]
    pF = eta_phys / (T0 * L0)                   pressure factor  [Pa per sim unit]
    bF = L0 ** 2                                bending factor   [N·m per sim unit]

Usage:
    python convert_units.py params-Node2d.json --to-physical -o params-Node2d-physical.json
    python convert_units.py params-Node2d-physical.json --to-sim -o params-Node2d.json
"""

import argparse
import json
import sys
from pathlib import Path

# ── Parameters that require unit conversion ──────────────────────────────
# Exactly matches ComputeConversionFactors() in CryptBuddingParams.hpp
#
# NOTE: gammaStemScale / gammaTransitScale / gammaDiffScale are dimensionless
# scaling factors and are NOT converted. Only gammaApical/Basal/Lateral are
# surface tensions in N/m and ARE converted.

STIFFNESS_PARAMS = {
    "cellGhostSpringStiffness",
    "springStiffness",
    "ghostE0",
    "ghostE1",
    "gammaApical",
    "gammaBasal",
    "gammaLateral",
    "nhMembraneSurface",
    "nhCellCellAdhesion",
    "nhBoundaryAdhesion",
    "nhStemStemAdhesion",
    "nhStemTransitAdhesion",
    "nhStemDiffAdhesion",
    "nhTransitTransitAdhesion",
    "nhTransitDiffAdhesion",
    "nhDiffDiffAdhesion",
    "nhStemBoundaryAdhesion",
    "nhTransitBoundaryAdhesion",
    "nhDiffBoundaryAdhesion",
}

PRESSURE_PARAMS = {
    "lumenPressure",
    # apicalConstrictionStrength is NOT converted: F = k * area_difference where
    # area is in [CD²], making the coefficient a geometric stiffness (η/CD·h), not pressure.
    "polarityBendingStrength",
}

BENDING_PARAMS = {
    "bendingStiffness",
}

ALL_CONVERTED = STIFFNESS_PARAMS | PRESSURE_PARAMS | BENDING_PARAMS

# ── C++ SetDefaults() values (simulation units) ───────────────────────────
# Must stay in sync with CryptBuddingParams.hpp SetDefaults().
# These are used when a param is absent from the source JSON: the physical
# output will include an explicit physical-unit value derived from the default
# sim value, preventing C++ from applying the conversion to its own default.
CPP_SIM_DEFAULTS = {
    # Stiffness [sim units, converted by /kF in C++]
    "cellGhostSpringStiffness": 5.0,
    "springStiffness":          30.0,
    "ghostE0":               5.0,
    "ghostE1":               1.0,
    "gammaApical":               0.85,
    "gammaBasal":                0.85,
    "gammaLateral":              0.7,
    "nhMembraneSurface":        10.0,
    "nhCellCellAdhesion":        1.0,
    "nhBoundaryAdhesion":        2.0,
    "nhStemStemAdhesion":        1.0,
    "nhStemTransitAdhesion":     1.0,
    "nhStemDiffAdhesion":        1.0,
    "nhTransitTransitAdhesion":  1.0,
    "nhTransitDiffAdhesion":     1.0,
    "nhDiffDiffAdhesion":        1.0,
    "nhStemBoundaryAdhesion":    2.0,
    "nhTransitBoundaryAdhesion": 2.0,
    "nhDiffBoundaryAdhesion":    2.0,
    # Pressure [sim units, converted by /pF in C++]
    "lumenPressure":             2.0,
    "polarityBendingStrength":   0.3,
    # Bending [sim units, converted by /bF in C++]
    "bendingStiffness":          5.0,
}

UNIT_LABELS = {}
for p in STIFFNESS_PARAMS:
    UNIT_LABELS[p] = "N/m"
for p in PRESSURE_PARAMS:
    UNIT_LABELS[p] = "Pa"
for p in BENDING_PARAMS:
    UNIT_LABELS[p] = "N·m"


def compute_factors(cell_diam_um: float, force_nN: float, vel_um_h: float):
    """Return (kF, pF, bF) conversion factors."""
    L0 = cell_diam_um * 1e-6
    T0 = 3600.0
    f_ref = force_nN * 1e-9
    v_ref = (vel_um_h * 1e-6) / T0
    eta = f_ref / v_ref

    kF = eta / T0
    pF = eta / (T0 * L0)
    bF = L0 * L0
    return kF, pF, bF


def get_factor(param_name: str, kF: float, pF: float, bF: float) -> float:
    if param_name in STIFFNESS_PARAMS:
        return kF
    if param_name in PRESSURE_PARAMS:
        return pF
    if param_name in BENDING_PARAMS:
        return bF
    raise ValueError(f"Unknown converted parameter: {param_name}")


def get_unit_refs(data: dict):
    uc = data.get("unitConversion", {})
    cell_diam = uc.get("cellDiameterMicrometers", 10.0)
    force = uc.get("referenceForcenN", 10.0)
    vel = uc.get("referenceVelocityUmPerH", 10.0)
    return cell_diam, force, vel


def is_physical(data: dict) -> bool:
    uc = data.get("unitConversion", {})
    return uc.get("usePhysicalUnits", False)


def collect_leaf_keys(obj, found=None):
    """Collect all leaf (non-object) keys from a nested dict."""
    if found is None:
        found = set()
    if isinstance(obj, dict):
        for k, v in obj.items():
            if isinstance(v, (dict, list)):
                collect_leaf_keys(v, found)
            else:
                found.add(k)
    elif isinstance(obj, list):
        for item in obj:
            collect_leaf_keys(item, found)
    return found


def walk_and_convert(obj, factor_fn, direction: str):
    """
    Recursively walk JSON, converting leaf values whose key is in ALL_CONVERTED.
    direction: 'to_physical' multiplies by factor, 'to_sim' divides by factor.
    Drops annotation keys starting with '_'.
    """
    if isinstance(obj, dict):
        out = {}
        for k, v in obj.items():
            if k.startswith("_"):
                continue
            if isinstance(v, (dict, list)):
                out[k] = walk_and_convert(v, factor_fn, direction)
            elif k in ALL_CONVERTED and isinstance(v, (int, float)):
                factor = factor_fn(k)
                if direction == "to_physical":
                    out[k] = v * factor
                else:
                    out[k] = v / factor
            else:
                out[k] = v
        return out
    elif isinstance(obj, list):
        return [walk_and_convert(item, factor_fn, direction) for item in obj]
    else:
        return obj


def format_value(v):
    if isinstance(v, float):
        return round(v, 12)
    return v


def clean_floats(obj):
    if isinstance(obj, dict):
        return {k: clean_floats(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [clean_floats(v) for v in obj]
    elif isinstance(obj, float):
        return format_value(obj)
    return obj


def main():
    parser = argparse.ArgumentParser(
        description="Convert CryptBudding JSON configs between physical and simulation units."
    )
    parser.add_argument("input", help="Input JSON config file")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--to-physical", action="store_true",
                       help="Convert simulation units → physical (SI) units")
    group.add_argument("--to-sim", action="store_true",
                       help="Convert physical (SI) units → simulation units")
    parser.add_argument("-o", "--output", help="Output file (default: stdout)")
    parser.add_argument("--cell-diameter", type=float, default=None,
                        help="Override cellDiameterMicrometers (µm)")
    parser.add_argument("--force", type=float, default=None,
                        help="Override referenceForcenN (nN)")
    parser.add_argument("--velocity", type=float, default=None,
                        help="Override referenceVelocityUmPerH (µm/h)")

    args = parser.parse_args()

    with open(args.input) as f:
        data = json.load(f)

    src_physical = is_physical(data)
    if args.to_physical and src_physical:
        print("Warning: input already has usePhysicalUnits=true. Converting anyway.", file=sys.stderr)
    if args.to_sim and not src_physical:
        print("Warning: input has usePhysicalUnits=false. Converting anyway.", file=sys.stderr)

    cell_diam, force_nN, vel = get_unit_refs(data)
    if args.cell_diameter is not None:
        cell_diam = args.cell_diameter
    if args.force is not None:
        force_nN = args.force
    if args.velocity is not None:
        vel = args.velocity

    kF, pF, bF = compute_factors(cell_diam, force_nN, vel)
    factor_fn = lambda name: get_factor(name, kF, pF, bF)
    direction = "to_physical" if args.to_physical else "to_sim"

    result = walk_and_convert(data, factor_fn, direction)

    # ── Ensure all converted params are explicit in physical output ───────
    # When usePhysicalUnits=true, C++ applies ComputeConversionFactors() to ALL
    # converted params — including those not in the JSON, using its sim-unit
    # defaults. To prevent those defaults from being wrongly reconverted, the
    # physical config must explicitly set every converted parameter.
    if args.to_physical:
        present_keys = collect_leaf_keys(data)
        missing = {p for p in ALL_CONVERTED if p not in present_keys}

        if missing:
            defaults_section = {}
            for param in sorted(missing):
                sim_default = CPP_SIM_DEFAULTS[param]
                phys_value = sim_default * factor_fn(param)
                defaults_section[param] = phys_value
                print(f"  Added default: {param} = {sim_default} sim → {phys_value:.6g} {UNIT_LABELS[param]}",
                      file=sys.stderr)

            # Inject into a dedicated section so the origin is clear.
            # The C++ parser flattens all nested keys, so placement doesn't matter.
            result["_defaults_from_cpp_SetDefaults"] = defaults_section

    # Update the usePhysicalUnits flag
    if "unitConversion" not in result:
        result["unitConversion"] = {}
    result["unitConversion"]["usePhysicalUnits"] = args.to_physical

    if args.cell_diameter is not None:
        result["unitConversion"]["cellDiameterMicrometers"] = cell_diam
    if args.force is not None:
        result["unitConversion"]["referenceForcenN"] = force_nN
    if args.velocity is not None:
        result["unitConversion"]["referenceVelocityUmPerH"] = vel

    result = clean_floats(result)

    # Print conversion summary
    L0 = cell_diam * 1e-6
    T0 = 3600.0
    eta = (force_nN * 1e-9) / ((vel * 1e-6) / T0)
    print(f"Reference: cellDiam={cell_diam} µm, F_ref={force_nN} nN, v_ref={vel} µm/h", file=sys.stderr)
    print(f"  L0={L0:.2e} m  η={eta:.4f} N·s/m", file=sys.stderr)
    print(f"  kF={kF:.6e} N/m per sim  (stiffness)", file=sys.stderr)
    print(f"  pF={pF:.6e} Pa  per sim  (pressure)", file=sys.stderr)
    print(f"  bF={bF:.6e} N·m per sim  (bending)", file=sys.stderr)
    print(f"Direction: {direction}", file=sys.stderr)

    output_str = json.dumps(result, indent=2)

    if args.output:
        Path(args.output).write_text(output_str + "\n")
        print(f"Written to {args.output}", file=sys.stderr)
    else:
        print(output_str)


if __name__ == "__main__":
    main()
