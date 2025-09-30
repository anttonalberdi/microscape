# MicroScape

**MicroScape** is a computational framework for integrative spatial modelling of microbial landscapes.

## Quick start

### Installation
```bash
conda env create -f https://raw.githubusercontent.com/anttonalberdi/microscape/refs/heads/main/envs/microscape.yaml
conda activate microscape
microscape --help
```

### Run the synthetic demo
```bash
microscape validate-sdp examples/00_synthetic/sdp_demo
microscape simulate examples/00_synthetic/sim_config.yml --out outputs/run_001.npz
```
