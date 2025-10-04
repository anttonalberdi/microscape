# MicroScape

_A computational framework for integrative spatial modelling of microbial landscapes._

[![Status](https://img.shields.io/badge/status-alpha-informational)](#roadmap) [![OS](https://img.shields.io/badge/os-Linux%20%7C%20macOS-lightgrey)](#requirements) [![Python](https://img.shields.io/badge/python-%E2%89%A53.10-blue)](#requirements) [![License](https://img.shields.io/badge/license-See%20LICENSE-lightgrey)](./LICENSE)

> MicroScape is in **ALPHA DEVELOPMENTAL PHASE**

---

## ✨ What is MicroScape?
MicroScape helps you **map**, **simulate**, and **analyse** micro‑scale spatial microbiome data. It is designed for workflows such as 3D host–microbe mapping and spatial ecology where you need to:

- **Map** different spatial omic datasets, including spatial metagenomics, metatranscriptomics and metabolomics.
- **Simulate** microbial landscapes under realistic constraints.
- **Analyse** microbiome functioning, microbe-microbe and host-microbiota interactions

---

## 🚀 Quick start

### Requirements
- **Python** ≥ 3.10 (via Conda recommended).
- Unix‑like OS (Linux or macOS).
- ~2–4 GB RAM for the synthetic demo.

### Install (Conda)
```bash
conda env create -f https://raw.githubusercontent.com/anttonalberdi/microscape/refs/heads/main/envs/microscape.yaml
conda activate microscape
microscape --help
```

### Run the synthetic demo
```bash
# 1) Validate an example Spatial Design Parameter (SDP) folder
microscape validate -i demo/system.yml

# 2) Simulate a spatially explicit microbial landscape
microscape profile -i demo/system.yml -o results/

# 2) Simulate a spatially explicit microbial landscape
microscape calibrate -i demo/system.yml -o results/

---

## 🧰 CLI overview
All functionality is exposed via the `microscape` command:

```
microscape --help
```

Common subcommands (run with `-h` to see full flags):

- `validate-sdp <PATH>` – checks your spatial design specification (folder or YAML) for completeness and consistency.
- `simulate <CONFIG.yml> --out <FILE.npz>` – runs a stochastic simulation using the provided configuration; accepts options like `--seed`, `--n-replicates`, and `--progress`.
- `export <FILE.npz> --format csv` – converts simulation artifacts to analysis‑ready tables.
- `plot <FILE.npz> --what abundance` – produces quick‑look figures for sanity checks.

> Tip: Keep configs under version control (`examples/` shows suggested structure). Use `--seed` for reproducibility.

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

---

## 📘 Configuration guide (SDP & simulation)

### 1) Spatial Design Parameters (SDP)
An SDP captures your planned sampling layout and constraints. A minimal SDP includes:
- **grid / coordinates**: spatial grid or point set
- **sample plan**: number of samples per region/stratum
- **biological priors**: optional community/gradient priors

Validate with:
```bash
microscape validate-sdp path/to/sdp
```

### 2) Simulation config (`sim_config.yml`)
Key fields you’ll likely set:
- `random_seed`: integer for reproducibility
- `n_replicates`: number of independent simulations
- `landscape`: model name and parameters
- `measurement`: noise model and depth
- `export`: which tables to write out

Run with:
```bash
microscape simulate sim_config.yml --out outputs/run.npz
```

---

## 📈 Example: quick visualization

```bash
microscape plot outputs/run_001.npz --what abundance --save fig_abundance.png
```
---

## 🧾 License
Distributed under the terms described in [LICENSE](./LICENSE).

## 🆘 Support & contact
- Issues & bugs: please open a GitHub issue
- Questions & ideas: open a discussion or reach out to the maintainers

---

## 🙏 Acknowledgements
This project is developed in the context of international collaboration between the University of Copenhagen and Kiel University on spatial host–microbe research.

