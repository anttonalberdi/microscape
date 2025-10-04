# MicroScape

_A computational framework for integrative spatial modelling of microbial landscapes._

[![Status](https://img.shields.io/badge/status-alpha-informational)](#roadmap) [![OS](https://img.shields.io/badge/os-Linux%20%7C%20macOS-lightgrey)](#requirements) [![Python](https://img.shields.io/badge/python-%E2%89%A53.10-blue)](#requirements) [![License](https://img.shields.io/badge/license-See%20LICENSE-lightgrey)](./LICENSE)

---

## ✨ What is MicroScape?
MicroScape helps you **simulate**, **validate** and **analyze** micro‑scale spatial microbiome data. It is designed for workflows such as 3D host–microbe mapping and spatial ecology where you need to:

- Define **spatial design parameters (SDP)** and validate them before running experiments
- **Simulate** microbial landscapes under realistic constraints
- Export results in analysis‑ready formats (e.g., `.npz`, `.csv`) and visualize basic diagnostics

> If you are looking for a practical way to prototype spatial sampling designs and benchmark downstream analyses on synthetic ground truth, MicroScape is for you.

---

## 🚀 Quick start

### Requirements
- **Python** ≥ 3.10 (via Conda recommended)
- Unix‑like OS (Linux or macOS)
- ~2–4 GB RAM for the synthetic demo

### Install (Conda)
```bash
conda env create -f https://raw.githubusercontent.com/anttonalberdi/microscape/refs/heads/main/envs/microscape.yaml
conda activate microscape
microscape --help
```

### Run the synthetic demo
```bash
# 1) Validate an example Spatial Design Parameter (SDP) folder
microscape validate-sdp examples/00_synthetic/sdp_demo

# 2) Simulate a spatially explicit microbial landscape
microscape simulate examples/00_synthetic/sim_config.yml --out outputs/run_001.npz
```

> The `outputs/run_001.npz` file can be loaded in Python or exported with the CLI subcommands below for plotting and tabular summaries.

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

## 📦 Using MicroScape in Python
Although the CLI is the recommended interface, you can load outputs programmatically:

```python
import numpy as np
z = np.load("outputs/run_001.npz")
print(list(z.keys()))  # e.g., coordinates, abundance, taxonomy, metadata
```

---

## 🧪 Development

### Set up a dev environment
```bash
# Clone your fork
git clone https://github.com/<you>/microscape
cd microscape

# Create the environment used for development
conda env create -f envs/microscape.yaml
conda activate microscape

# Install in editable mode (if needed)
pip install -e .
```

### Run tests & linters
```bash
pytest -q
```

> Consider enabling pre‑commit hooks to enforce style/quality. We follow standard Python typing and docstrings.

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

## 🤝 Contributing
Contributions are welcome! Please open an issue first to discuss the change you’d like to make. Typical contributions include:
- New landscape or measurement models
- Additional exporters and validators
- Documentation improvements & examples

### Pull request checklist
- [ ] Tests added/updated
- [ ] Docs updated (`README`, `docs/`, or `examples/`)
- [ ] CI passes locally (`pytest`)

---

## 🗺️ Roadmap
- [ ] Public API stabilization of core models
- [ ] Expanded demo datasets and tutorials
- [ ] Benchmarks against empirical datasets
- [ ] Optional GPU‑accelerated kernels where relevant

Have suggestions? Please open a discussion or an issue.

---

## 📣 Citing MicroScape
If you use MicroScape in a publication, please cite this repository and the associated preprint/article when available. For now:

```
Alberdi A., et al. MicroScape: integrative spatial modelling of microbial landscapes. GitHub repository: https://github.com/anttonalberdi/microscape
```

> Also consider citing relevant methods describing micro‑scale spatial metagenomics and 3D host–microbiota mapping.

---

## 🧾 License
Distributed under the terms described in [LICENSE](./LICENSE).

## 🆘 Support & contact
- Issues & bugs: please open a GitHub issue
- Questions & ideas: open a discussion or reach out to the maintainers

---

## 🙏 Acknowledgements
This project is developed in the context of spatial host–microbe research. We thank contributors and testers who provide feedback via issues and PRs.

