import importlib.resources as ir

def get_packaged_example(rel_path: str) -> str:
    """
    Return the full path to an example file bundled with microscape.

    Example:
        get_packaged_example("examples/00_synthetic/community_sbml.yml")
    """
    return str(ir.files("microscape").joinpath(rel_path))
