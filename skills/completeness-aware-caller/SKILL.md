---
name: completeness-aware-caller
description: >-
  Three-state gene-presence calls (present / absent / cannot-conclude) and
  ANI-ranking gates for incomplete genomes and MAGs. Refuses to call a gene
  absent, or to rank references by ANI, when assembly completeness cannot
  support the claim.
license: MIT
metadata:
  version: "0.1.0"
  author: Qihao Duann
  domain: genomics
  tags:
    - comparative-genomics
    - completeness
    - metagenomics
    - mag
    - ani
    - abstention
  inputs:
    - name: completeness
      type: value
      format:
        - float
      description: >-
        Assembly completeness as fraction (0-1] or percent (1-100], or via
        --busco-json pointing at a busco-assessor result.json (required
        unless --demo)
      required: true
    - name: checkm2_tsv
      type: file
      format:
        - tsv
      description: >-
        CheckM2 quality_report.tsv; supplies completeness AND contamination
        from one run. Mutually exclusive with --busco-json. Use --genome to
        pick a row when the report covers several genomes.
      required: false
    - name: genome
      type: value
      format:
        - string
      description: >-
        Row (Name column) to select when --checkm2-tsv covers several genomes
      required: false
    - name: contamination
      type: value
      format:
        - float
      description: >-
        Contamination percent (0-100), overriding the CheckM2 value; gates
        `present` calls against MIMAG tiers
      required: false
    - name: genes
      type: file
      format:
        - json
      description: 'Gene search results: [{"gene": "nifH", "found": false}, ...]'
      required: false
  outputs:
    - name: report
      type: file
      format: md
      description: Three-state call table with rationale and refusal messages
    - name: result
      type: file
      format: json
      description: Machine-readable calls, gate decision, and thresholds used
  dependencies:
    python: ">=3.10"
    packages: []
    external:
      - checkm2>=1.1.0 (optional; only for --checkm2-tsv, via scripts/install_checkm2.sh)
      - CheckM2 UniRef100/KO database, 3.1 GB (Zenodo record 14897628; not a conda package)
  demo_data:
    - path: examples/demo_genes.json
      description: STM815 degradation-benchmark gene set (nifH/nodC lost at 70%)
    - path: examples/checkm2_quality_report.tsv
      description: >-
        Real CheckM2 1.1.0 output for the four STM815 degradation levels
        (completeness 100.00/87.41/74.30/54.66, contamination <1%)
  endpoints:
    cli: python skills/completeness-aware-caller/completeness_aware_caller.py --busco-json {busco_result} --genes {genes_json} --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
    install:
      - kind: conda
        package: checkm2=1.1.0
        channels:
          - bioconda
          - conda-forge
    always: false
    emoji: "🚦"
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - darwin
      - linux
    trigger_keywords:
      - "is this gene really absent"
      - "completeness-aware"
      - "can I trust this MAG"
      - "gene absence incomplete genome"
      - "cannot conclude"
      - "abstention genomics"
      - "MAG interpretation"
      - "ANI confidence"
---

# 🚦 Completeness-Aware Caller

You are the **completeness-aware-caller**, a ClawBio skill that separates
taxonomic confidence from functional confidence in incomplete genomes. A
fragmented genome makes missing genes look absent and close relatives look
functionally identical; this skill's deliverable is the refusal — an explicit
`CANNOT CONCLUDE` printed as a first-class result, not a footnote.

## Trigger

**Fire when the user says any of:**
- "is this gene really absent", "gene absent or just missing"
- "can I trust this MAG", "MAG completeness interpretation"
- "completeness-aware call", "abstention", "cannot conclude"
- "is the ANI difference meaningful", "ANI confidence"
- "does incompleteness affect my conclusion"

**Do NOT fire when:**
- User wants to *measure* completeness → route to `busco-assessor`
  (eukaryotes, or marker-based) or CheckM2 via `scripts/run_checkm2.sh`
  (prokaryotic genomes and MAGs, and the only source of contamination)
- User wants to *compute* ANI or relatedness → run FastANI (or
  `galaxy-bridge`), then come back here with the values
- User wants variant-level interpretation → `clinical-variant-reporter`
- User wants microbiome community profiling → `claw-metagenomics`

## Why This Exists

- **Without it**: pipelines print "gene absent" from 70%-complete MAGs and
  rank references by ANI margins smaller than the noise incompleteness
  induces. Absence of evidence is silently converted into evidence of
  absence.
- **With it**: every absence claim and every ANI ranking carries either a
  confidence grounded in measured completeness, or an explicit refusal with
  an actionable message.
- **Why ClawBio**: chains directly after `busco-assessor` (reads its
  `result.json` verbatim) and complements `galaxy-bridge` ANI workflows.

## Scope

One skill, one task: **decide whether completeness supports a claim**. This
skill does not measure completeness, compute ANI, search for genes, or
assemble genomes. It consumes those results and gates the conclusions.

## Domain Decisions

1. **Absence needs ≥95% completeness** (`MIN_COMPLETENESS_FOR_ABSENT =
   0.95`). Under a random-loss model, P(gene missed | truly present) ≈
   1 − completeness; 0.95 caps that at 5%. Anchored to MIMAG quality tiers
   (Bowers et al. 2017) where "high-quality draft" starts at >90%.
2. **ANI drift model**: drift = k × (1 − completeness) ANI points,
   calibrated on the *P. phymatum* STM815 degradation benchmark bundled with
   this skill's provenance (0.41 points of observed drift at 50% retention
   against both *P. xenovorans* LB400 and *B. cenocepacia* J2315).
   **k depends on which tool measured completeness**, because BUSCO and
   CheckM2 do not report the same number for the same assembly:
   `DRIFT_COEFF = 0.82` for BUSCO C%, `DRIFT_COEFF_CHECKM2 = 0.98` for
   CheckM2 Completeness (0.88 against ground-truth base retention, for
   reference). Fitted by least squares through the origin on four retention
   levels × two references; derivation and raw tool output at
   https://github.com/Qihao-Duan/completeness-aware-caller
   (`benchmark/benchmark_summary.md`, `benchmark/checkm2/`).
5. **Contamination gates presence, not absence** (MIMAG tiers, Bowers et al.
   2017): <5% → `present` unchanged; 5–10% → `present` with confidence
   reduced to 1 − contamination; ≥10% → `cannot_conclude`, because a detected
   contig can no longer be attributed to this genome rather than the
   contaminating one. Contamination adds sequence, so it never affects an
   absence call.
3. **Safety factor 2**: a ranking margin must exceed 2× modelled drift —
   signal must beat twice the noise.
4. **Species-boundary zone 94–96%** (Jain et al. 2018): ANI values in this
   window with completeness <90% flag species assignment as uncertain.

## Workflow

1. Parse completeness from `--completeness` (fraction or percent), from a
   `busco-assessor` `result.json` via `--busco-json` (field `C`), or from a
   CheckM2 `quality_report.tsv` via `--checkm2-tsv` (columns `Completeness`
   and `Contamination`; `--genome` selects a row in a multi-genome report).
   Reject values outside (0, 1] / (1, 100]. Record which tool supplied the
   number — it selects the drift coefficient in step 3.
2. For each gene in `--genes`: found → `present`, subject to the
   contamination gate (≥10% contamination downgrades it to
   `cannot_conclude`); not found and completeness ≥ 0.95 → `absent` with
   confidence = completeness; otherwise → `cannot_conclude` with the miss
   probability spelled out.
3. If `--ani-a`/`--ani-b` supplied: compute margin and modelled drift using
   the coefficient fitted for the completeness source from step 1, then
   decide `rank` vs `cannot_conclude`; flag the species-boundary zone. A
   hand-supplied `--completeness` has no calibration of its own, so the gate
   falls back to the BUSCO coefficient and is marked `"calibration":
   "assumed"` in `result.json`.
4. Write `report.md` (call table + gate + disclaimer), `result.json`
   (calls, thresholds, decision), `commands.sh` (replay).

## CLI Reference

```bash
# One-time: install CheckM2 and its 3.1 GB reference database
sh skills/completeness-aware-caller/scripts/install_checkm2.sh

# Score genomes, then gate the calls on completeness AND contamination
sh skills/completeness-aware-caller/scripts/run_checkm2.sh \
  --input ./bins --output /tmp/checkm2_out
python skills/completeness-aware-caller/completeness_aware_caller.py \
  --checkm2-tsv /tmp/checkm2_out/quality_report.tsv --genome bin_042 \
  --genes genes.json --ani-a 81.45 --ani-b 80.75 --output /tmp/cac_out

# From a busco-assessor run + gene search results + FastANI values
python skills/completeness-aware-caller/completeness_aware_caller.py \
  --busco-json /path/to/busco_out/result.json \
  --genes genes.json --ani-a 81.45 --ani-b 80.75 --output /tmp/cac_out

# Direct completeness value
python skills/completeness-aware-caller/completeness_aware_caller.py \
  --completeness 0.98 --genes genes.json --output /tmp/cac_out

# Demo: STM815 degradation benchmark at 70% completeness
python skills/completeness-aware-caller/completeness_aware_caller.py \
  --demo --output /tmp/cac_demo
```

## Example Output

```markdown
# Completeness-Aware Caller Report

**Assembly completeness**: 68.1%

| Gene | Status | Confidence | Rationale |
|------|--------|-----------|-----------|
| nifH | CANNOT CONCLUDE | 0.68 | CANNOT CONCLUDE absence: the assembly is only 68% complete, so a truly present gene would be missed with ~32% probability. |
| nifD | PRESENT | 1.00 | Detected in the assembly; detection is positive evidence. |

## ANI ranking gate
**Decision**: RANK
Margin 0.70 ANI points exceeds 2x modelled drift (0.26); ranking is supported.
```

(Real run: BUSCO C:68.1% on STM815 degraded to 70% retention — where the
naive call for nifH would have been a false "absent".)

## Gotchas

1. **The model will want to treat "not found" as "absent".** Do not. The
   whole point of this skill is that at 70% completeness a missing gene is
   missing with ~30% probability even when truly present. Always route
   "not found" through `call_gene`, never report absence directly.
2. **The model will want to reuse the drift coefficient across clades.**
   The default 0.82 was calibrated on one Burkholderiaceae benchmark. For
   distant clades or eukaryotes, recalibrate (degrade a complete genome of
   the target clade and measure) or say the default is a Burkholderiaceae
   calibration.
3. **BUSCO C% is itself noisy at high retention.** In the calibration
   benchmark, dropping 9.7% of bases produced M=19.8% (marker loss is
   clustered, not uniform). Treat completeness as an estimate with variance,
   which is exactly why the thresholds are conservative.
4. **`present` is only contamination-safe when contamination is supplied.**
   Detection confidence 1.0 assumes the contig truly belongs to the genome.
   `--checkm2-tsv` supplies contamination and the gate applies automatically;
   the BUSCO path does not, because BUSCO does not measure contamination. If
   you gate on BUSCO completeness, pass `--contamination` explicitly or treat
   every `present` as ungated.
5. **Never mix one tool's completeness with the other's drift coefficient.**
   CheckM2 read 87.4% where BUSCO read 78.4% on the same assembly. Reusing
   `DRIFT_COEFF` (0.82, BUSCO) with CheckM2 completeness understates the
   threshold by ~16%, so the gate ranks when it should abstain — at frag50 it
   would wrongly rank any margin in [0.741, 0.885). The skill selects the
   coefficient from the input flag; do not override `drift_coeff` by hand
   unless you have fitted your own.
6. **CheckM2 runs optimistic on heavily fragmented assemblies.** Against
   ground-truth base retention it read +4.8 and +4.2 points high at 70% and
   50% retention. It is still the more accurate estimator overall (MAE 2.95
   vs BUSCO's 4.26), but a completeness just over 95% from a fragmented bin
   deserves scepticism before absence is asserted.

## Safety

- All processing is local; no sequence data leaves the machine.
- Every `report.md` ends with the ClawBio disclaimer: *"ClawBio is a
  research and educational tool. It is not a medical device and does not
  provide clinical diagnoses. Consult a healthcare professional before
  making any medical decisions."*
- The skill never invents completeness or ANI values; it only gates numbers
  produced by upstream tools.

## Agent Boundary

- **Agent dispatches**: completeness (or busco-assessor result.json path),
  gene search results, optional ANI pair.
- **Skill executes**: threshold logic, three-state calls, gate decision,
  report writing.
- **Agent explains**: what the refusal means for the user's question and
  what evidence would unlock a conclusion.
- **Agent must NOT**: override a `cannot_conclude` into a definitive call,
  tweak thresholds ad hoc, or report absence for genes the skill refused.

## Chaining Partners

| Upstream | Handoff | Downstream |
|----------|---------|-----------|
| CheckM2 (`scripts/run_checkm2.sh`) | `quality_report.tsv` (completeness + contamination) | `completeness-aware-caller` |
| `busco-assessor` | `result.json` (C%) | `completeness-aware-caller` |
| FastANI / `galaxy-bridge` | ANI values | `completeness-aware-caller` |
| `completeness-aware-caller` | `result.json` three-state calls | `profile-report`, `lit-synthesizer` |

## Maintenance

- **Review cadence**: quarterly, or when BUSCO/MIMAG standards revise.
- **Staleness signals**: busco-assessor changes its `result.json` schema;
  new consensus thresholds for MAG quality; drift recalibrations published.
- **Deprecation path**: fold into a broader confidence-gating skill if one
  emerges; keep the three-state vocabulary stable.

## Citations

- Bowers R.M. et al. (2017). Minimum information about a single amplified
  genome (MISAG) and a metagenome-assembled genome (MIMAG). *Nature
  Biotechnology*. https://doi.org/10.1038/nbt.3893
- Jain C. et al. (2018). High throughput ANI analysis of 90K prokaryotic
  genomes reveals clear species boundaries. *Nature Communications*.
  https://doi.org/10.1038/s41467-018-07641-9
- Manni M. et al. (2021). BUSCO Update. *Molecular Biology and Evolution*.
  https://doi.org/10.1093/molbev/msab199
- Chklovski A. et al. (2023). CheckM2: a rapid, scalable and accurate tool
  for assessing microbial genome quality using machine learning. *Nature
  Methods*. https://doi.org/10.1038/s41592-023-01940-w
