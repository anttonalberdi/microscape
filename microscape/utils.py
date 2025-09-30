# microscape/utils.py
from __future__ import annotations
from contextlib import ExitStack
from importlib.resources import files, as_file
from pathlib import Path

def get_demo_dir(subpath: str = "00_synthetic/sdp_demo") -> Path:
    """
    Return a real filesystem path to the bundled demo directory.
    Works from both source and installed wheels.
    """
    data_root = files("microscape").joinpath("examples").joinpath(subpath)
    # as_file handles zip/packed wheels gracefully
    with ExitStack() as stack:
        p = stack.enter_context(as_file(data_root))
        # We return the path itself; callers should use it immediately (don’t cache across processes).
        return Path(str(p))
