# MicroScape

**MicroScape** is a computational framework for  integrative spatial modelling of microbial landscapes.

## Quick start

### Conda environments with dependencies
```bash
# minimal stack
conda env create -f envs/microscape.yaml
conda activate microscape
microscape --help

```

### Installation
```bash
conda activate microscape-min
pip install git+https://github.com/anttonalberdi/microscape.git
```

## Run the synthetic demo
```bash
microscape validate-sdp examples/00_synthetic/sdp_demo
microscape simulate examples/00_synthetic/sim_config.yml --out outputs/run_001.npz
```
