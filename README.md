# MicroScape

**MicroScape** is a computational framework for integrative spatial modelling of microbial landscapes.

## Quick start

### Installation

The best way to install **MicroScape** is to create a conda environment containing **MicroScape** and all its dependencies.

```bash
conda env create -f https://raw.githubusercontent.com/anttonalberdi/microscape/refs/heads/main/envs/microscape.yaml
conda activate microscape
microscape --help
```

### Run the synthetic demo

```bash
# Download data
wget https://raw.githubusercontent.com/anttonalberdi/microscape/refs/heads/main/examples/demo1/demo1.yml
wget https://raw.githubusercontent.com/anttonalberdi/microscape/refs/heads/main/examples/demo1/BP.xml
wget https://raw.githubusercontent.com/anttonalberdi/microscape/refs/heads/main/examples/demo1/FD.xml
wget https://raw.githubusercontent.com/anttonalberdi/microscape/refs/heads/main/examples/demo1/LU.xml

microscape simulate demo1.yml
```

https://github.com/anttonalberdi/microscape/tree/branch/foldername

`microscape demo` runs a tiny, biologically motivated simulation on a small tissue slice:

- Creates a synthetic **fibre patch** in the gut lumen.
- Uses a simple **metapathway proxy** to produce **butyrate** near fibre (a stand-in for “fibres → SCFAs”).
- Diffuses butyrate through tissue with **stable no-flux boundaries** (epithelium doesn’t teleport metabolites).
- Computes an **SCFA-at-mucosa** score and a **radial profile** (butyrate vs distance from the mucosa band).
- Saves **human-readable plots** and a **CSV summary**, plus the raw arrays.

