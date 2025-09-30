from importlib.resources import files
from pathlib import Path

def get_packaged_path(rel_path: str) -> str:
    """
    Return an absolute path to a resource shipped inside the 'microscape' package.
    Falls back to repo-relative path when running from source.
    """
    # 1) Try package data (installed wheel)
    try:
        res = files("microscape").joinpath(rel_path)
        # Traversable returned by importlib.resources
        if hasattr(res, "is_file") and res.is_file():
            return str(res)
    except Exception:
        pass

    # 2) Fallback: resolve relative to the repo root when running from source
    here = Path(__file__).resolve().parents[1]  # microscape/
    candidate = (here / rel_path).resolve()
    return str(candidate)

# Backwards-compatible alias used elsewhere in the codebase
def get_packaged_example(rel_path: str) -> str:
    return get_packaged_path(rel_path)