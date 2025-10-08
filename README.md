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


### 6) `build` — assemble the modeling table

**What it does.** Joins per‑spot features with per‑microbe targets from a chosen metabolism JSON. Each **row = (spot, microbe)**.

**Standard use:**
```bash
# Predict growth (objective)
microscape build demo/system.yml   --metabolism-json test/metabolism_constrained_combined.json   --outdir out/model --target objective

# Or: predict a flux (e.g., uptake/secretion of C0011)
microscape build demo/system.yml   --metabolism-json test/metabolism_constrained_combined.json   --outdir out/model_flux --target flux:EX_C0011_e
```

**Also supported.**
- Multiple targets with combination (e.g., total uptake across metabolites):
  ```bash
  microscape build demo/system.yml     --metabolism-json test/metabolism_constrained_combined.json     --outdir out/model_multi     --target flux:EX_C0011_e --target flux:EX_C0012_e     --target-op uptake_sum
  ```
- Progress bar: `--no-progress` to disable.

**Outputs (in `--outdir`):**
- `table.csv` — tidy dataset for training (columns include `spot_id`, `microbe`, `treatment`, `cage`, `env_id`, `abundance`, `met:*`, and a single `target`).
- `schema.json` — provenance + column types.

---

### 7) `train` — fit a hierarchical Bayesian model

**What it does.** Learns a regression with:
- fixed effects for **features** (standardized; `abundance` is `log1p` before z‑scoring),
- optional fixed effect for **treatment**,
- random intercepts for **cage** and **environment**.

**Standard use:**
```bash
microscape train out/model/table.csv --outdir out/model_run
```

**Scaling & robustness options (recommended for large data):**
```bash
# fewer processes (avoid OOM), gentle init, light shrinkage on betas
microscape train out/model/table.csv --outdir out/model_run   --cores 1 --init adapt_diag --target-accept 0.9 --shrinkage
# or use NumPyro NUTS (single process) if JAX is available
microscape train out/model/table.csv --outdir out/model_run --jax
```

**Outputs (in `--outdir`):**
- `posterior.nc` — posterior draws (ArviZ NetCDF).
- `model_card.json` — features used, group encodings, training stats (means/SDs), and a summary table.

---

### 8) `evaluate` — check fit and interpret coefficients

**What it does.** Uses the saved posterior + the training table to compute:
- per‑row **posterior means/HDIs** and **residuals**,
- **RMSE/MAE/R²** and HDI coverage; WAIC/LOO when available,
- summaries of **β coefficients**, **treatment effects**, and **random intercepts**.

**Standard use:**
```bash
microscape evaluate out/model/table.csv out/model_run/posterior.nc out/model_run/model_card.json   --outdir out/eval
```

**Outputs (in `--outdir`):**
- `residuals.csv` — observed vs. predicted (+ HDIs) per row.
- `coef_summary.csv` — α and βs with mean/SD/HDIs and Pr(>0)/Pr(<0).
- `treatment_effects.csv` — per‑level fixed effects (if modeled).
- `random_intercepts_cage.csv`, `random_intercepts_env.csv` — RE summaries (if modeled).
- `metrics.json` — RMSE/MAE/R²/coverage (+ WAIC/LOO if available).

---

### 9) `predict` — posterior predictions on any table

**What it does.** Produces posterior **conditional means** and **HDIs** for new rows (or the training rows), reusing training standardization and encodings from the model card.

**Standard use (in‑sample):**
```bash
microscape predict out/model_run/posterior.nc out/model_run/model_card.json   --outdir out/pred
```

**Predict on new data:**
```bash
microscape predict out/model_run/posterior.nc out/model_run/model_card.json   --table out/model_new/table.csv   --outdir out/pred_new
```

**Notes.**
- Unknown group levels default to **zero random effect**; switch to strict with `--new-level-policy error`.
- Add predictive intervals (noise via `sigma`) with `--ppc`.

**Outputs (in `--outdir`):**
- `predictions.csv` — `y_hat_mean` + HDIs (and, with `--ppc`, predictive intervals).
- `predict_meta.json` — run configuration and feature list.

---


### 10) `microscape whatif` — counterfactual scenarios

**What it does.** Creates *what-if* versions of selected rows (spot–microbe pairs), applies edits to features and/or group labels, and returns posterior **means** and **HDIs** for **baseline**, **what‑if**, and their **delta** — without re‑running metabolism.

**Standard uses:**
```bash
# Change one metabolite for a specific spot & microbe
microscape whatif out/model_run/posterior.nc out/model_run/model_card.json   --select "spot_id=S0001" --select "microbe=M0002"   --set "met:C0010=8.0"   --outdir out/whatif_S1_M2

# "Remove" a microbe everywhere (set its abundance to 0 for all matching rows)
microscape whatif out/model_run/posterior.nc out/model_run/model_card.json   --select "microbe=M0004"   --set "abundance=0"   --limit 1000000   --outdir out/whatif_drop_M0004

# Zero a metabolite across all rows (global intervention)
microscape whatif out/model_run/posterior.nc out/model_run/model_card.json   --set "met:C0011=0"   --limit 1000000   --outdir out/whatif_zero_metC0011_all

# Apply changes only within a treatment group
microscape whatif out/model_run/posterior.nc out/model_run/model_card.json   --select "treatment=diet_B"   --set "met:C0011=0"   --outdir out/whatif_zero_metC0011_dietB
```

**Arguments you’ll care about.**
- `--set col=value` overwrites a column (e.g., `met:C0010=8.0`, `treatment=diet_B`).
- `--delta col=+/-x` adds a numeric delta (e.g., `abundance=+10`, `met:C0004=-1.5`).
- `--select col=value` filters rows; repeat to AND conditions.
- `--limit` caps how many rows are processed (default **50**). Use a large value (e.g., `--limit 1000000`) for global edits.
- `--new-level-policy zero|error` controls unseen group levels (default **zero** → no random effect).

**Outputs (in `--outdir`):**
- `whatif.csv` — baseline vs. what‑if (`y_base_mean/HDI`, `y_new_mean/HDI`) and `delta_mean/HDI` per row; includes before/after values for changed columns.
- `whatif_meta.json` — run configuration and feature list.

> Note: *what-if* operates on the **trained model** (counterfactual prediction). To recompute mechanistic steady states under removals/edits, change inputs and re‑run `microscape metabolism` (and downstream modules).

---

## 🧪 Practical tips
- **Compare runs.** Run `metabolism` once without constraints and once with; diff objectives and key `flux:EX_*` columns to see how constraints reshape phenotype.
- **Audit quickly.** If something looks off for a pair, open the matching `constraints__*__debug.json` to see exactly how its bounds were computed.
- **Keep provenance.** Store the `constraints__{mode}__reactions.json` used alongside its corresponding metabolism outputs.

---

## 🙏 Acknowledgements
This project is developed in the context of international collaboration between the University of Copenhagen and Kiel University on spatial host–microbe research.

