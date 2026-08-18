#!/usr/bin/env python3
"""completeness-aware-caller — three-state gene calls and ANI gating for
incomplete genomes.

Turns (gene search results, a completeness estimate, optional ANI
comparisons) into calls that refuse to overreach: a gene not found in a
genome that is only 70% complete is NOT absent — it is `cannot_conclude`.
Likewise, ranking two references by ANI is refused when the margin between
them is within the drift that incompleteness alone can induce.

Domain decisions (see SKILL.md for citations):
- Absence may be asserted only when completeness >= 0.95, so that under a
  random-loss model the probability of having missed a truly present gene
  stays at or below 5%.
- ANI drift is modelled as DRIFT_COEFF * (1 - completeness), calibrated on
  the STM815 degradation benchmark (0.41 ANI points of drift at 50%
  retention). A ranking must exceed SAFETY_FACTOR times that drift.
- The drift coefficient is specific to the tool that measured completeness.
  BUSCO and CheckM2 do not report the same number for the same assembly, so
  each has its own coefficient fitted on the same benchmark. Mixing one
  tool's completeness with the other's coefficient mis-scales the gate.
- ANI values inside the 94-96% species-boundary zone with completeness
  below 0.9 are flagged as uncertain species assignments.
- Detection ("present") is only trustworthy if the contig belongs to the
  genome, so `present` is gated on contamination using the MIMAG tiers.
"""

import argparse
import csv
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

MIN_COMPLETENESS_FOR_ABSENT = 0.95
SAFETY_FACTOR = 2.0
SPECIES_BOUNDARY = (94.0, 96.0)
BOUNDARY_MIN_COMPLETENESS = 0.90

# ANI points of drift per unit incompleteness, fitted by least squares through
# the origin on the STM815 degradation benchmark (4 retention levels x 2
# references). The value depends on which tool produced the completeness
# number, because the two tools disagree on the same assembly:
#
#   completeness scale   coefficient   R^2
#   true base retention      0.88      0.964   (physical reference)
#   BUSCO C%                 0.82      0.997
#   CheckM2 Completeness     0.98      0.984
#
# Derivation and raw tool output:
# https://github.com/Qihao-Duan/completeness-aware-caller (benchmark/)
DRIFT_COEFF = 0.82           # BUSCO C% (the original published calibration)
DRIFT_COEFF_CHECKM2 = 0.98   # CheckM2 Completeness
DRIFT_COEFFS = {"busco": DRIFT_COEFF, "checkm2": DRIFT_COEFF_CHECKM2}

# MIMAG quality tiers (Bowers et al. 2017): high-quality drafts allow <5%
# contamination, medium-quality <10%. Above 10% a contig cannot be assumed to
# belong to the genome under assessment, so detection stops being evidence
# about *this* organism.
CONTAM_HIGH_QUALITY = 5.0
CONTAM_MEDIUM_QUALITY = 10.0

DISCLAIMER = (
    "*ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare "
    "professional before making any medical decisions.*"
)


def parse_completeness(value):
    """Accept a fraction (0-1] or a percentage (1-100]; return a fraction."""
    c = float(value)
    if 1.0 < c <= 100.0:
        c /= 100.0
    if not 0.0 < c <= 1.0:
        raise ValueError(
            f"completeness {value} outside (0, 1] as fraction or (1, 100] as percent"
        )
    return c


def parse_checkm2_report(path, genome=None):
    """Read completeness and contamination from a CheckM2 quality_report.tsv.

    CheckM2 scores every genome it was given, so a report may hold many rows;
    `genome` selects one by its Name column. Returns percentages as CheckM2
    writes them (0-100), leaving normalisation to parse_completeness.
    """
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise ValueError(f"CheckM2 report {path} contains no rows")
    for required in ("Name", "Completeness", "Contamination"):
        if required not in rows[0]:
            raise ValueError(
                f"{path} is not a CheckM2 quality_report.tsv "
                f"(missing column {required!r})"
            )

    if genome is not None:
        matches = [r for r in rows if r["Name"] == genome]
        if not matches:
            available = ", ".join(sorted(r["Name"] for r in rows))
            raise ValueError(
                f"genome {genome!r} not found in {path}. Available: {available}"
            )
        row = matches[0]
    elif len(rows) > 1:
        available = ", ".join(sorted(r["Name"] for r in rows))
        raise ValueError(
            f"{path} holds {len(rows)} genomes; pass --genome to choose one. "
            f"Available: {available}"
        )
    else:
        row = rows[0]

    return {
        "genome": row["Name"],
        "completeness": float(row["Completeness"]),
        "contamination": float(row["Contamination"]),
        "model_used": row.get("Completeness_Model_Used"),
    }


def contamination_verdict(contamination):
    """Map a contamination percentage onto a MIMAG-anchored trust tier."""
    if contamination is None:
        return "unknown"
    if contamination < CONTAM_HIGH_QUALITY:
        return "clean"
    if contamination < CONTAM_MEDIUM_QUALITY:
        return "moderate"
    return "high"


def call_gene(found, completeness, contamination=None,
              min_completeness=MIN_COMPLETENESS_FOR_ABSENT):
    """Three-state call for one gene given a completeness estimate.

    `contamination` (percent, as CheckM2 reports it) gates positive calls:
    detection only supports a claim about *this* genome if the contig it was
    found on plausibly belongs to it. Absence calls are unaffected —
    contamination adds sequence, it does not remove it.
    """
    if found:
        # Only a supplied contamination figure can gate a positive call; with
        # none (the BUSCO path) the original ungated behaviour stands.
        if contamination is not None:
            if contamination >= CONTAM_MEDIUM_QUALITY:
                return {
                    "status": "cannot_conclude",
                    "confidence": 0.0,
                    "message": (
                        f"CANNOT CONCLUDE presence: the bin is "
                        f"{contamination:.1f}% contaminated, above the 10% "
                        f"MIMAG medium-quality bound, so a detected sequence "
                        f"cannot be attributed to this genome rather than to "
                        f"the contaminating one. Decontaminate or re-bin "
                        f"before asserting the gene is carried here."
                    ),
                }
            if contamination >= CONTAM_HIGH_QUALITY:
                return {
                    "status": "present",
                    "confidence": 1.0 - contamination / 100.0,
                    "message": (
                        f"Detected, but the bin carries {contamination:.1f}% "
                        f"contamination (MIMAG medium-quality): there is a "
                        f"residual chance the contig belongs to another "
                        f"organism."
                    ),
                }
        return {
            "status": "present",
            "confidence": 1.0,
            "message": "Detected in the assembly; detection is positive evidence.",
        }
    miss_probability = 1.0 - completeness
    if completeness >= min_completeness:
        return {
            "status": "absent",
            "confidence": completeness,
            "message": (
                f"Not found at {completeness:.1%} completeness; residual chance "
                f"of a missed gene is ~{miss_probability:.1%}."
            ),
        }
    return {
        "status": "cannot_conclude",
        "confidence": completeness,
        "message": (
            f"CANNOT CONCLUDE absence: the assembly is only "
            f"{completeness * 100:.0f}% complete, so a truly present gene "
            f"would be missed with ~{miss_probability:.0%} probability. "
            f"Re-sequence or close the assembly before asserting loss of "
            f"function."
        ),
    }


def ani_margin_gate(ani_a, ani_b, completeness, source="busco",
                    drift_coeff=None, safety=SAFETY_FACTOR):
    """Decide whether ANI to reference A vs B may be ranked at this
    completeness. Margin must exceed safety * modelled drift.

    `source` names the tool that measured completeness and selects the drift
    coefficient fitted for that tool's scale. An explicit `drift_coeff`
    overrides it. A completeness supplied by hand has no calibration of its
    own, so it falls back to the BUSCO coefficient and the result is marked
    "assumed".
    """
    if drift_coeff is None:
        drift_coeff = DRIFT_COEFFS.get(source, DRIFT_COEFF)
    calibration = "fitted" if source in DRIFT_COEFFS else "assumed"
    margin = abs(ani_a - ani_b)
    drift = drift_coeff * (1.0 - completeness)
    lo, hi = SPECIES_BOUNDARY
    boundary_uncertain = (
        completeness < BOUNDARY_MIN_COMPLETENESS
        and any(lo <= ani <= hi for ani in (ani_a, ani_b))
    )
    if margin >= safety * drift:
        decision = "rank"
        message = (
            f"Margin {margin:.2f} ANI points exceeds {safety:.0f}x modelled "
            f"drift ({drift:.2f}); ranking is supported."
        )
    else:
        decision = "cannot_conclude"
        message = (
            f"CANNOT CONCLUDE ranking: margin {margin:.2f} ANI points is "
            f"within {safety:.0f}x the drift ({drift:.2f}) that "
            f"{(1 - completeness):.0%} incompleteness alone can induce. "
            f"Do not report which reference is closer."
        )
    return {
        "decision": decision,
        "margin": margin,
        "drift": drift,
        "drift_coeff": drift_coeff,
        "completeness_source": source,
        "calibration": calibration,
        "species_boundary_uncertain": boundary_uncertain,
        "message": message,
    }


# ------------------------------------------------------------------ reporting
def write_outputs(outdir, completeness, gene_calls, ani_gate, argv,
                  contamination=None, source="busco", genome=None):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    result = {
        "completeness": completeness,
        "contamination": contamination,
        "completeness_source": source,
        "genome": genome,
        "gene_calls": gene_calls,
        "ani_gate": ani_gate,
        "thresholds": {
            "min_completeness_for_absent": MIN_COMPLETENESS_FOR_ABSENT,
            "drift_coeff": DRIFT_COEFFS.get(source, DRIFT_COEFF),
            "safety_factor": SAFETY_FACTOR,
            "contamination_high_quality": CONTAM_HIGH_QUALITY,
            "contamination_medium_quality": CONTAM_MEDIUM_QUALITY,
        },
        "generated": datetime.now(timezone.utc).isoformat(),
    }
    (outdir / "result.json").write_text(json.dumps(result, indent=2))

    lines = [
        "# Completeness-Aware Caller Report",
        "",
    ]
    if genome:
        lines.append(f"**Genome**: {genome}")
        lines.append("")
    lines.append(
        f"**Assembly completeness**: {completeness:.1%} "
        f"(source: {source})"
    )
    lines.append("")
    if contamination is not None:
        tier = contamination_verdict(contamination)
        note = {
            "clean": "below the 5% MIMAG high-quality bound",
            "moderate": "MIMAG medium-quality (5-10%)",
            "high": "above the 10% MIMAG medium-quality bound",
        }[tier]
        lines.append(f"**Contamination**: {contamination:.2f}% — {note}")
        lines.append("")
    lines += [
        "## Gene calls",
        "",
        "| Gene | Status | Confidence | Rationale |",
        "|------|--------|-----------|-----------|",
    ]
    for g in gene_calls:
        lines.append(
            f"| {g['gene']} | {g['status'].upper().replace('_', ' ')} | "
            f"{g['confidence']:.2f} | {g['message']} |"
        )
    if ani_gate is not None:
        lines += [
            "",
            "## ANI ranking gate",
            "",
            f"**Decision**: {ani_gate['decision'].upper().replace('_', ' ')}",
            "",
            ani_gate["message"],
        ]
        if ani_gate.get("calibration") == "assumed":
            lines.append(
                "\n⚠️ Completeness was supplied directly, so no fitted drift "
                f"coefficient applies; the BUSCO calibration "
                f"({ani_gate['drift_coeff']}) was assumed. Supply "
                "--busco-json or --checkm2-tsv for a calibrated gate."
            )
        if ani_gate["species_boundary_uncertain"]:
            lines.append(
                "\n⚠️ At least one ANI value falls in the 94-96% species-"
                "boundary zone while completeness is below 90%: species "
                "assignment is uncertain."
            )
    lines += ["", "---", "", DISCLAIMER, ""]
    (outdir / "report.md").write_text("\n".join(lines))

    (outdir / "commands.sh").write_text(
        "#!/bin/sh\n# Replay command\n" + " ".join(argv) + "\n"
    )


# ------------------------------------------------------------------ demo
# Mirrors the real frag70 level of the STM815 benchmark. nodB/nodC/nodS/nodU
# lose every copy and vanish from the assembly; nifH survives only because
# STM815 carries two copies of it on pBPHY02 — copy number, not importance,
# decides which genes a fragmented assembly can still show you.
DEMO_GENES = [
    {"gene": "nodC", "found": False},
    {"gene": "nodB", "found": False},
    {"gene": "nifH", "found": True},
    {"gene": "nifD", "found": True},
]
DEMO_COMPLETENESS = 0.70
DEMO_ANI = (81.45, 80.75)  # STM815@70% vs LB400 and vs J2315 (benchmark)


def run(args, argv):
    contamination = None
    genome = None
    source = "busco"

    if args.demo:
        completeness = DEMO_COMPLETENESS
        gene_specs = DEMO_GENES
        ani_pair = DEMO_ANI
    else:
        if args.busco_json and args.checkm2_tsv:
            raise SystemExit(
                "Use either --busco-json or --checkm2-tsv, not both: each has "
                "its own drift calibration and mixing them is undefined."
            )
        if args.checkm2_tsv:
            report = parse_checkm2_report(args.checkm2_tsv, args.genome)
            completeness = parse_completeness(report["completeness"])
            contamination = report["contamination"]
            genome = report["genome"]
            source = "checkm2"
        elif args.busco_json:
            busco = json.loads(Path(args.busco_json).read_text())
            completeness = parse_completeness(busco["C"])
        elif args.completeness:
            completeness = parse_completeness(args.completeness)
            source = "manual"
        else:
            raise SystemExit(
                "Provide --completeness, --busco-json, --checkm2-tsv, or --demo"
            )
        # An explicit --contamination wins, so a BUSCO run can still be gated
        # on a contamination figure obtained elsewhere.
        if args.contamination is not None:
            contamination = args.contamination
        gene_specs = (
            json.loads(Path(args.genes).read_text()) if args.genes else []
        )
        ani_pair = (
            (args.ani_a, args.ani_b)
            if args.ani_a is not None and args.ani_b is not None
            else None
        )

    gene_calls = [
        {"gene": spec["gene"],
         **call_gene(spec["found"], completeness, contamination)}
        for spec in gene_specs
    ]
    ani_gate = (
        ani_margin_gate(ani_pair[0], ani_pair[1], completeness, source=source)
        if ani_pair else None
    )
    write_outputs(args.output, completeness, gene_calls, ani_gate, argv,
                  contamination=contamination, source=source, genome=genome)
    print(f"Report written to: {Path(args.output) / 'report.md'}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--completeness", help="Fraction (0-1] or percent (1-100]")
    ap.add_argument("--busco-json",
                    help="Path to a busco-assessor result.json (reads C%%)")
    ap.add_argument("--checkm2-tsv",
                    help="Path to a CheckM2 quality_report.tsv (reads "
                         "Completeness and Contamination)")
    ap.add_argument("--genome",
                    help="Row to use when --checkm2-tsv holds several genomes")
    ap.add_argument("--contamination", type=float,
                    help="Contamination percent; overrides the CheckM2 value")
    ap.add_argument("--genes",
                    help='JSON file: [{"gene": "nifH", "found": false}, ...]')
    ap.add_argument("--ani-a", type=float,
                    help="ANI to reference A (e.g. from FastANI)")
    ap.add_argument("--ani-b", type=float, help="ANI to reference B")
    ap.add_argument("--output", required=True, help="Output directory")
    ap.add_argument("--demo", action="store_true",
                    help="Run the bundled STM815@70%% benchmark scenario")
    args = ap.parse_args()
    # Bad inputs are a user problem, not a bug: report them as a message
    # rather than a traceback.
    try:
        run(args, sys.argv)
    except (ValueError, KeyError) as exc:
        raise SystemExit(f"error: {exc}") from exc
    except FileNotFoundError as exc:
        raise SystemExit(f"error: no such file: {exc.filename}") from exc


if __name__ == "__main__":
    main()
