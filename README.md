# MicroScape

_A computational framework for integrative spatial modelling of microbial landscapes._

[![Status](https://img.shields.io/badge/status-alpha-informational)](#roadmap) [![OS](https://img.shields.io/badge/os-Linux%20%7C%20macOS-lightgrey)](#requirements) [![Python](https://img.shields.io/badge/python-%E2%89%A53.10-blue)](#requirements) [![License](https://img.shields.io/badge/license-See%20LICENSE-lightgrey)](./LICENSE)

> MicroScape is in **ALPHA DEVELOPMENTAL PHASE**

---

## ✨ What is MicroScape?
MicroScape helps you **map**, **simulate**, and **analyse** micro‑scale spatial microbiome data. It is designed for workflows such as 3D host–microbe mapping and spatial ecology where you need to:

- **Map** different spatial omic datasets, including spatial metagenomics, metatranscriptomics and metabolomics.
- **Simulate** microbial landscapes under realistic constraints.
- **Analyse** microbiome functioning, microbe–microbe and host–microbiota interactions.

---

## 🚀 Quick start

### Requirements
- **Python** ≥ 3.10 (Conda recommended).
- Unix‑like OS (**Linux** or **macOS**).
- ~2–4 GB RAM for the synthetic demo.

### Install (Conda, recommended)
```bash
conda env create -f https://raw.githubusercontent.com/anttonalberdi/microscape/refs/heads/main/envs/microscape.yaml
conda activate microscape
microscape --help
```

### Alternative install (pip)
```bash
pip install git+https://github.com/anttonalberdi/microscape.git
microscape --help
```

---

## 🗂️ Repository layout
```
microscape/
├─ envs/            # Conda environment definitions
├─ examples/        # Ready‑to‑run demos (synthetic first)
├─ microscape/      # Source code package
├─ pyproject.toml   # Build metadata
└─ README.md        # This file
```

## 🗂️ Example folder layout (after a full run)

```
test/
├─ ecology_summary.csv
├─ constraints__environmental.csv
├─ constraints__environmental__reactions.json
├─ constraints__environmental__debug.json
├─ constraints__transcriptional.csv
├─ constraints__transcriptional__reactions.json
├─ constraints__transcriptional__debug.json
├─ constraints__combined.csv
├─ constraints__combined__reactions.json
├─ constraints__combined__debug.json
├─ metabolism_unconstrained.csv
├─ metabolism_unconstrained.json
├─ metabolism_constrained_combined.csv
├─ metabolism_constrained_combined.json
├─ community_metrics_spot.csv
├─ community_pairwise.csv
└─ community_microbe.csv
```

---

### 🔁 Copy‑paste demo (synthetic example)

```bash
# Validate the configuration
microscape validate demo/system.yml

# Optional: ecological summaries
microscape ecology demo/system.yml --outdir test

# Build constraints (choose one or run all)
microscape constrain demo/system.yml --outdir test --mode environmental
microscape constrain demo/system.yml --outdir test --mode transcriptional
microscape constrain demo/system.yml --outdir test --mode combined

# Run metabolism WITHOUT constraints (baseline)
microscape metabolism demo/system.yml --outdir test

# Run metabolism WITH constraints (example: combined mode)
microscape metabolism demo/system.yml \
  --outdir test \
  --constraints test/constraints__combined__reactions.json

# Community‑level metrics from a named metabolism JSON
microscape community test/metabolism_unconstrained.json --outdir test
```

---

### 1) `validate` — sanity‑check your configuration

**What it does.**  
Parses your `system.yml` and referenced paths to ensure everything needed by downstream steps is findable and well‑formed.

**Use (standard example):**
```bash
microscape validate Github/microscape/examples/demo/system.yml
```

**Output highlights.**
- No files are written. You should see `OK` on success.
- On failure, the error message points to the exact missing/malformed field (fix here before proceeding).

---

### 2) *(optional)* `ecology` — quick per‑spot ecological snapshot

**What it does.**  
Builds simple per‑spot/per‑microbe summaries (presence/abundance and other ready‑to‑plot fields) to help you inspect the dataset before modelling.

**Use (standard example):**
```bash
microscape ecology Github/microscape/examples/demo/system.yml --outdir test
```

**Output highlights.**
- `test/ecology_summary.csv` — one row per **spot × microbe** with basic ecology fields.
- Use this file to check spot coverage and microbe representation before adding constraints.

---

### 3) `constrain` — turn measurements into model bounds

**What it does.**  
Builds **reaction bounds** for each **spot × microbe**. You can choose the source of constraints:
- **environmental** — extracellular concentrations → exchange bounds.
- **transcriptional** — TPM + GPR rules → tightens **internal** reaction bounds.
- **combined** — applies both environmental and transcriptional bounds.

**Use (standard examples):**
```bash
# environmental
microscape constrain Github/microscape/examples/demo/system.yml --outdir test --mode environmental

# transcriptional
microscape constrain Github/microscape/examples/demo/system.yml --outdir test --mode transcriptional

# combined
microscape constrain Github/microscape/examples/demo/system.yml --outdir test --mode combined
```

**Output highlights (per mode):**
- `test/constraints__{mode}.csv` — per **spot × microbe** summary with:
  - `changed_ex`: how many **exchange** reactions were altered;
  - `changed_internal`: how many **internal** reactions were altered;
  - `warnings`: parsing or data‑shape notes, if any.
- `test/constraints__{mode}__reactions.json` — **compact** map of **only changed** reactions and their **final** bounds.  
  *This is the file you pass to the metabolism step via `--constraints`.*
- `test/constraints__{mode}__debug.json` — **rich** per‑reaction audit trail (base bounds, env vs tx contributions, final).  
  Use this to understand *why* a reaction was constrained.

> Notes: Transcript inputs are auto‑detected whether provided per‑microbe, flat shared, or under common keys (e.g., `all`).

---

### 4) `metabolism` — steady‑state FBA per spot × microbe

**What it does.**  
Solves a COBRA model for each **spot × microbe**, applying:
1) per‑spot exchange bounds derived from your metabolism rules, and  
2) (optionally) the **final reaction bounds** from `constrain`.

You can run a **baseline** (unconstrained) and one or more **constrained** runs.

**Use (standard examples):**
```bash
# baseline (no constraints)
microscape metabolism Github/microscape/examples/demo/system.yml --outdir test

# with constraints (example: combined mode)
microscape metabolism Github/microscape/examples/demo/system.yml \
  --outdir test \
  --constraints test/constraints__combined__reactions.json
```

**Output highlights.**
- **Naming is automatic** and reflects constraint state:
  - `test/metabolism_unconstrained.{csv,json}` — no constraints applied.
  - `test/metabolism_constrained_{environmental|transcriptional|combined|custom}.{csv,json}` — constraints applied; `{mode}` read from the JSON’s `"mode"`.
- **CSV content (per row = spot × microbe):**
  - `spot_id`, `microbe`, `status`, `objective` (biomass or chosen objective).
  - A set of `flux:EX_*` columns for recorded exchanges. *Sign convention:* uptake < 0, secretion > 0.
- **JSON** mirrors the CSV per spot, and includes audit arrays:
  - `applied_env`: list of exchange bounds set from concentrations.
  - `applied_constraints`: list of reaction bounds overridden by the constraints file.
- **Strict by default:** if `--constraints` is malformed or irrelevant to the current run (no matching spot × microbe, or unknown reactions), the command errors out. Use `--no-constraints-strict` to relax.

---

### 5) `community` — competition, complementarity & cross‑feeding

**What it does.**  
Aggregates the **exchange fluxes** from a named metabolism JSON to quantify **niche overlap** (competition), **niche separation** (complementarity) and **cross‑feeding** within each spot.

**Use (standard example):**
```bash
# Point to any named metabolism JSON produced above:
microscape community test/metabolism_unconstrained.json --outdir test
# add --graph to also write a simple network JSON
```

**Output highlights.**
- `test/community_metrics_spot.csv` — one row per spot with:
  - mean/median **competition** (Schoener’s D over resource‑use profiles),
  - mean/median **complementarity** (= 1 − D),
  - **resource redundancy** (how shared resources are across the community),
  - cross‑feeding **connectance** (density of i→j feeding edges),
  - totals: uptake, secretion, and sum of objectives.
- `test/community_pairwise.csv` — directed pairs (i→j) per spot with:
  - `competition_schoenerD`, `complementarity`,
  - `crossfeeding_ij`, `crossfeeding_ji`, and `net_crossfeeding_i_minus_j`,
  - `uptake_jaccard` (overlap of resources consumed).
- `test/community_microbe.csv` — per **spot × microbe**:
  - `unique_resource_fraction`, `total_uptake`, `total_secretion`, `objective`.
- *(optional)* `test/community_graph.json` (with `--graph`) — bipartite microbe↔exchange resource edges and microbe→microbe cross‑feeding edges (weighted).

> Exchange flux interpretation: **uptake = max(0, −flux)**, **secretion = max(0, +flux)**.

---

## 🧪 Practical tips
- **Compare runs.** Run `metabolism` once without constraints and once with; diff objectives and key `flux:EX_*` columns to see how constraints reshape phenotype.
- **Audit quickly.** If something looks off for a pair, open the matching `constraints__*__debug.json` to see exactly how its bounds were computed.
- **Keep provenance.** Store the `constraints__{mode}__reactions.json` used alongside its corresponding metabolism outputs.

---

## 🙏 Acknowledgements
This project is developed in the context of international collaboration between the University of Copenhagen and Kiel University on spatial host–microbe research.

