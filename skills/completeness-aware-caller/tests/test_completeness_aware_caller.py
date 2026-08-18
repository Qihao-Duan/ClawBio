"""Red/green TDD suite for completeness-aware-caller.

The skill turns (gene search result, completeness estimate, ANI comparisons)
into three-state calls: present / absent / cannot_conclude. Absence may only
be asserted when completeness makes a miss unlikely; ANI-based ranking of two
references may only be asserted when the margin exceeds completeness-induced
drift.
"""

import json
import subprocess
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "completeness_aware_caller.py"

sys.path.insert(0, str(SKILL_DIR))

from completeness_aware_caller import (
    DRIFT_COEFF,
    DRIFT_COEFF_CHECKM2,
    ani_margin_gate,
    call_gene,
    parse_checkm2_report,
    parse_completeness,
)


# ---------------------------------------------------------------- gene calls
def test_present_at_any_completeness():
    for c in (0.5, 0.7, 0.99):
        result = call_gene(found=True, completeness=c)
        assert result["status"] == "present"


def test_absent_requires_high_completeness():
    result = call_gene(found=False, completeness=0.98)
    assert result["status"] == "absent"
    assert result["confidence"] == pytest.approx(0.98)


def test_low_completeness_abstains():
    result = call_gene(found=False, completeness=0.70)
    assert result["status"] == "cannot_conclude"
    assert "70" in result["message"]
    assert "cannot" in result["message"].lower()


def test_boundary_exactly_at_threshold_allows_absent():
    result = call_gene(found=False, completeness=0.95)
    assert result["status"] == "absent"


def test_completeness_percent_normalised():
    assert parse_completeness("72.5") == pytest.approx(0.725)
    assert parse_completeness("0.725") == pytest.approx(0.725)


def test_invalid_completeness_raises():
    with pytest.raises(ValueError):
        parse_completeness("150")
    with pytest.raises(ValueError):
        parse_completeness("-0.1")


# ---------------------------------------------------------------- ANI gate
def test_ani_gate_refuses_when_margin_within_noise():
    # Real numbers from the STM815 degradation benchmark: at 50% retention the
    # LB400/J2315 margin (0.70 ANI pts) sits inside 2x the observed drift.
    result = ani_margin_gate(ani_a=81.28, ani_b=80.58, completeness=0.50)
    assert result["decision"] == "cannot_conclude"
    assert result["margin"] == pytest.approx(0.70, abs=0.01)


def test_ani_gate_ranks_when_margin_is_clear():
    result = ani_margin_gate(ani_a=97.0, ani_b=82.0, completeness=0.90)
    assert result["decision"] == "rank"


def test_ani_gate_flags_species_boundary_zone():
    result = ani_margin_gate(ani_a=95.2, ani_b=82.0, completeness=0.70)
    assert result["species_boundary_uncertain"] is True


# ---------------------------------------------------------------- CLI / demo
def test_demo_mode_produces_report_with_refusal(tmp_path):
    out = tmp_path / "demo_out"
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--demo", "--output", str(out)],
        capture_output=True, text=True, timeout=120,
    )
    assert proc.returncode == 0, proc.stderr
    report = (out / "report.md").read_text()
    assert "CANNOT CONCLUDE" in report
    assert "not a medical device" in report
    assert (out / "result.json").exists()
    assert (out / "commands.sh").exists()


def test_demo_result_json_schema(tmp_path):
    out = tmp_path / "demo_json"
    subprocess.run(
        [sys.executable, str(SCRIPT), "--demo", "--output", str(out)],
        capture_output=True, text=True, timeout=120,
    )
    data = json.loads((out / "result.json").read_text())
    assert data["completeness"] == pytest.approx(0.70)
    statuses = {g["gene"]: g["status"] for g in data["gene_calls"]}
    # nodC loses every copy at frag70 -> unprovable absence, must abstain.
    assert statuses["nodC"] == "cannot_conclude"
    assert statuses["nodB"] == "cannot_conclude"
    # nifH survives frag70 only because STM815 carries two copies of it.
    assert statuses["nifH"] == "present"
    assert statuses["nifD"] == "present"
    assert data["ani_gate"]["decision"] in ("rank", "cannot_conclude")


def test_cli_with_real_inputs(tmp_path):
    genes = tmp_path / "genes.json"
    genes.write_text(json.dumps([
        {"gene": "nifH", "found": False},
        {"gene": "nodC", "found": True},
    ]))
    out = tmp_path / "real_out"
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--completeness", "0.98",
         "--genes", str(genes), "--output", str(out)],
        capture_output=True, text=True, timeout=120,
    )
    assert proc.returncode == 0, proc.stderr
    data = json.loads((out / "result.json").read_text())
    statuses = {g["gene"]: g["status"] for g in data["gene_calls"]}
    assert statuses["nifH"] == "absent"
    assert statuses["nodC"] == "present"


# ------------------------------------------------------- CheckM2 integration
# CheckM2 supplies completeness AND contamination from one run, which is what
# lets the skill gate `present` calls as well as `absent` ones.

CHECKM2_HEADER = (
    "Name\tCompleteness\tContamination\tCompleteness_Model_Used\t"
    "Translation_Table_Used\n"
)


def _checkm2_tsv(tmp_path, rows, name="quality_report.tsv"):
    """rows: [(genome, completeness, contamination), ...]"""
    body = "".join(
        f"{g}\t{c}\t{x}\tNeural Network (Specific Model)\t11\n"
        for g, c, x in rows
    )
    path = tmp_path / name
    path.write_text(CHECKM2_HEADER + body)
    return path


def test_checkm2_report_single_row_parsed(tmp_path):
    tsv = _checkm2_tsv(tmp_path, [("bin1", 74.3, 0.18)])
    got = parse_checkm2_report(tsv)
    assert got["genome"] == "bin1"
    assert got["completeness"] == pytest.approx(74.3)
    assert got["contamination"] == pytest.approx(0.18)


def test_checkm2_report_multi_row_requires_genome(tmp_path):
    tsv = _checkm2_tsv(tmp_path, [("bin1", 74.3, 0.18), ("bin2", 90.0, 1.0)])
    with pytest.raises(ValueError, match="pass --genome"):
        parse_checkm2_report(tsv)


def test_checkm2_report_selects_named_genome(tmp_path):
    tsv = _checkm2_tsv(tmp_path, [("bin1", 74.3, 0.18), ("bin2", 90.0, 1.0)])
    got = parse_checkm2_report(tsv, genome="bin2")
    assert got["completeness"] == pytest.approx(90.0)


def test_checkm2_report_unknown_genome_lists_available(tmp_path):
    tsv = _checkm2_tsv(tmp_path, [("bin1", 74.3, 0.18)])
    with pytest.raises(ValueError, match="bin1"):
        parse_checkm2_report(tsv, genome="nope")


def test_checkm2_report_rejects_non_checkm2_tsv(tmp_path):
    path = tmp_path / "other.tsv"
    path.write_text("a\tb\n1\t2\n")
    with pytest.raises(ValueError, match="not a CheckM2"):
        parse_checkm2_report(path)


# --------------------------------------------------------- contamination gate
def test_clean_bin_leaves_present_untouched():
    """Below the MIMAG 5% bound, behaviour must match the pre-CheckM2 skill."""
    assert call_gene(True, 0.9, contamination=1.0) == call_gene(True, 0.9)


def test_absent_contamination_is_backwards_compatible():
    """No contamination supplied => exactly the original behaviour."""
    assert call_gene(True, 0.9, contamination=None)["confidence"] == 1.0


def test_moderate_contamination_downweights_present():
    result = call_gene(True, 0.99, contamination=7.0)
    assert result["status"] == "present"
    assert result["confidence"] == pytest.approx(0.93)


def test_high_contamination_refuses_presence():
    result = call_gene(True, 0.99, contamination=12.0)
    assert result["status"] == "cannot_conclude"
    assert "CANNOT CONCLUDE presence" in result["message"]


def test_contamination_does_not_affect_absence():
    """Contamination adds sequence; it cannot explain a missing gene."""
    clean = call_gene(False, 0.99, contamination=0.1)
    dirty = call_gene(False, 0.99, contamination=30.0)
    assert clean["status"] == dirty["status"] == "absent"


def test_contamination_boundaries_are_inclusive_below():
    assert call_gene(True, 0.99, contamination=4.99)["confidence"] == 1.0
    assert call_gene(True, 0.99, contamination=5.0)["status"] == "present"
    assert call_gene(True, 0.99, contamination=9.99)["status"] == "present"
    assert call_gene(True, 0.99, contamination=10.0)["status"] == "cannot_conclude"


# ------------------------------------------------- estimator-aware drift model
def test_drift_coefficient_differs_by_completeness_source():
    """BUSCO and CheckM2 disagree on the same assembly, so each has its own
    coefficient. Reusing one for the other silently mis-scales the gate."""
    busco = ani_margin_gate(81.45, 80.75, 0.681, source="busco")
    checkm2 = ani_margin_gate(81.45, 80.75, 0.743, source="checkm2")
    assert busco["drift_coeff"] == DRIFT_COEFF
    assert checkm2["drift_coeff"] == DRIFT_COEFF_CHECKM2
    assert checkm2["drift_coeff"] > busco["drift_coeff"]


def test_manual_completeness_is_marked_uncalibrated():
    gate = ani_margin_gate(81.45, 80.75, 0.70, source="manual")
    assert gate["calibration"] == "assumed"
    assert gate["drift_coeff"] == DRIFT_COEFF


def test_known_sources_are_marked_fitted():
    for src in ("busco", "checkm2"):
        assert ani_margin_gate(81.45, 80.75, 0.70, source=src)["calibration"] == "fitted"


def test_explicit_drift_coeff_overrides_source():
    gate = ani_margin_gate(81.45, 80.75, 0.70, source="checkm2", drift_coeff=0.5)
    assert gate["drift_coeff"] == 0.5


# ------------------------------------------------------------------- CLI paths
def test_cli_rejects_both_completeness_sources(tmp_path):
    tsv = _checkm2_tsv(tmp_path, [("bin1", 74.3, 0.18)])
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--checkm2-tsv", str(tsv),
         "--busco-json", "whatever.json", "--output", str(tmp_path / "o")],
        capture_output=True, text=True, timeout=120,
    )
    assert proc.returncode != 0
    assert "not both" in (proc.stdout + proc.stderr)


def test_cli_checkm2_end_to_end(tmp_path):
    tsv = _checkm2_tsv(tmp_path, [("bin1", 74.3, 0.18)])
    genes = tmp_path / "genes.json"
    genes.write_text(json.dumps([
        {"gene": "nodC", "found": False},
        {"gene": "nifH", "found": True},
    ]))
    out = tmp_path / "cm2_out"
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--checkm2-tsv", str(tsv),
         "--genes", str(genes), "--output", str(out)],
        capture_output=True, text=True, timeout=120,
    )
    assert proc.returncode == 0, proc.stderr
    data = json.loads((out / "result.json").read_text())
    assert data["completeness_source"] == "checkm2"
    assert data["contamination"] == pytest.approx(0.18)
    assert data["genome"] == "bin1"
    statuses = {g["gene"]: g["status"] for g in data["gene_calls"]}
    assert statuses["nodC"] == "cannot_conclude"
    assert statuses["nifH"] == "present"
    assert "Contamination" in (out / "report.md").read_text()


def test_cli_high_contamination_flips_present_to_refusal(tmp_path):
    tsv = _checkm2_tsv(tmp_path, [("dirty", 99.0, 12.0)])
    genes = tmp_path / "genes.json"
    genes.write_text(json.dumps([{"gene": "nifH", "found": True}]))
    out = tmp_path / "dirty_out"
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--checkm2-tsv", str(tsv),
         "--genes", str(genes), "--output", str(out)],
        capture_output=True, text=True, timeout=120,
    )
    assert proc.returncode == 0, proc.stderr
    data = json.loads((out / "result.json").read_text())
    assert data["gene_calls"][0]["status"] == "cannot_conclude"


def test_cli_missing_file_is_a_message_not_a_traceback(tmp_path):
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--checkm2-tsv", str(tmp_path / "nope.tsv"),
         "--output", str(tmp_path / "o")],
        capture_output=True, text=True, timeout=120,
    )
    assert proc.returncode != 0
    combined = proc.stdout + proc.stderr
    assert "Traceback" not in combined
    assert "error:" in combined


def test_shipped_checkm2_example_is_readable():
    """The CheckM2 evidence shipped in examples/ must stay parseable."""
    tsv = SKILL_DIR / "examples" / "checkm2_quality_report.tsv"
    got = parse_checkm2_report(tsv, genome="stm815_frag70")
    assert got["completeness"] == pytest.approx(74.3, abs=0.1)
