# microscape/io/spot_loader.py
from pathlib import Path
import yaml

def load_spot(p: Path) -> dict:
    with p.open("r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh) or {}
    if "spot" in data:
        return data["spot"]
    return data
