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
├─ docs/            # Reference & how‑to docs (in progress)
├─ envs/            # Conda environment definitions
├─ examples/        # Ready‑to‑run demos (synthetic first)
├─ microscape/      # Source code package
├─ tests/           # Unit tests (pytest)
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

## 📦 Command reference & expected outputs

> Paths are relative to your `--outdir`. Filenames are deterministic so they’re easy to script against.

| Command | Purpose | Key outputs | Notes |
|---|---|---|---|
| `microscape update` | Refresh internal registries/caches | — | Prints progress; no user‑facing files. |
| `microscape validate <system.yml>` | Sanity‑check structure & paths | — | Prints `OK` on success; otherwise a clear error. |
| `microscape ecology <system.yml> --outdir DIR` | Per‑spot ecological summary | `ecology_summary.csv` | Optional QC/inspection step. |
| `microscape constrain <system.yml> --outdir DIR --mode {environmental\|transcriptional\|combined}` | Convert measurements → **reaction bounds** | `constraints__{mode}.csv` • `constraints__{mode}__reactions.json` • `constraints__{mode}__debug.json` | Use `__reactions.json` later with `metabolism --constraints`. `__debug.json` is a full per‑reaction audit trail. |
| `microscape metabolism <system.yml> --outdir DIR` | FBA per **spot × microbe** (baseline) | `metabolism_unconstrained.{csv,json}` | Records objective & selected exchange fluxes. |
| `microscape metabolism <system.yml> --outdir DIR --constraints constraints__{mode}__reactions.json` | FBA with **final bounds** applied | `metabolism_constrained_{mode}.{csv,json}` | `mode` auto‑detected from JSON (`environmental`, `transcriptional`, `combined`, or `custom`). Strict by default (see below). |
| `microscape community <named_metabolism.json> --outdir DIR` | **Community metrics** per spot | `community_metrics_spot.csv` • `community_pairwise.csv` • `community_microbe.csv` • *(optional)* `community_graph.json` | Uses exchange fluxes: uptake = `max(0, −flux)`, secretion = `max(0, +flux)`. Add `--graph` to emit a bipartite + cross‑feeding graph. |

---

## 🔍 Outputs—at a glance

### Constraints (`constraints__{mode}*`)
- **CSV**: one row per `spot_id × microbe` with `changed_ex`, `changed_internal`, `warnings`.
- **`__reactions.json` (compact)**: only **changed** reactions with **final** bounds, consumable by `metabolism --constraints`.

---

## 🙏 Acknowledgements
This project is developed in the context of international collaboration between the University of Copenhagen and Kiel University on spatial host–microbe research.

