from microscape.coupling.loop import run_minimal

def test_run_minimal():
    out = run_minimal(sim_steps=5)
    assert "butyrate" in out and out["butyrate"].shape[1] > 0
