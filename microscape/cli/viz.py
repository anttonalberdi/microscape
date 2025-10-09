# microscape/cli/viz.py
from __future__ import annotations
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
import json, datetime as dt

import typer
import pandas as pd
import yaml

# Re-use your existing helpers
from ..io.system_loader import load_system, iter_spot_files_for_env
from ..io.spot_loader import load_spot

app = typer.Typer(add_completion=False, no_args_is_help=True)

def _safe(val, cast=float):
    try:
        return cast(val)
    except Exception:
        return None

def _read_json_maybe(path: Optional[Path]) -> Optional[Dict[str, Any]]:
    if not path:
        return None
    try:
        return json.loads(Path(path).read_text())
    except Exception:
        return None

def _read_table_maybe(path: Optional[Path]) -> Optional[pd.DataFrame]:
    if not path:
        return None
    p = Path(path)
    try:
        if p.suffix.lower() == ".parquet":
            return pd.read_parquet(p)
        return pd.read_csv(p)
    except Exception:
        return None

def _index_metabolism(meta: Dict[str, Any]) -> Dict[str, Any]:
    """Normalize metabolism JSON into spots -> microbes -> {objective, fluxes{}}."""
    if not meta:
        return {}
    if "spots" in meta and isinstance(meta["spots"], dict):
        return meta["spots"]
    # fallback to flat dict keyed by spot_id
    return meta

def _collect_from_system(system_yml: Path) -> Tuple[List[Dict[str, Any]], Dict[str, Dict[str, float]], Dict[str, float]]:
    """
    Returns:
      spots: [{id, env_id, x, y, z}]
      abund_by_microbe: microbe -> {spot_id -> value}
      mets_by_id: metabolite_id -> {spot_id -> value}
    """
    sys_info = load_system(system_yml)
    spots_meta: List[Dict[str, Any]] = []
    abund_by_microbe: Dict[str, Dict[str, float]] = {}
    mets_by_id: Dict[str, Dict[str, float]] = {}

    for env_file in sys_info["environment_files"]:
        env_doc = yaml.safe_load(Path(env_file).read_text())
        env = env_doc.get("environment", {})
        env_id = env.get("id")
        for sid, spot_path in iter_spot_files_for_env(env_file, sys_info["paths"]):
            sobj = load_spot(spot_path) or {}
            s = sobj.get("spot") or sobj
            spot_id = str(s.get("id") or s.get("name") or sid)
            pos = (s.get("position") or {})
            sx = _safe(pos.get("x_um"))
            sy = _safe(pos.get("y_um"))
            sz = _safe(pos.get("z_um"))
            spots_meta.append({"id": spot_id, "env_id": env_id, "x": sx, "y": sy, "z": sz})

            meas = (s.get("measurements") or {})
            # microbes abundances (dict of microbe_id -> value)
            mics = (meas.get("microbes") or {}).get("values") or {}
            if isinstance(mics, dict):
                for mid, val in mics.items():
                    abund_by_microbe.setdefault(str(mid), {})[spot_id] = _safe(val)
            # metabolites (dict of metabolite_id -> value)
            mets = (meas.get("metabolites") or {}).get("values") or {}
            if isinstance(mets, dict):
                for cid, val in mets.items():
                    mets_by_id.setdefault(str(cid), {})[spot_id] = _safe(val)

    return spots_meta, abund_by_microbe, mets_by_id

def _collect_from_metabolism(metabolism_json: Optional[Path]) -> Tuple[Dict[str, Dict[str, float]], Dict[str, Dict[str, Dict[str, float]]], List[str]]:
    """
    Returns:
      objective_by_microbe: microbe -> {spot_id -> objective}
      flux_by_exchange: exchange_id -> microbe -> {spot_id -> flux}
      exchange_ids: list of exchange ids present
    """
    meta = _read_json_maybe(metabolism_json) or {}
    per_spot = _index_metabolism(meta)
    obj_by_microbe: Dict[str, Dict[str, float]] = {}
    flux_by_ex: Dict[str, Dict[str, Dict[str, float]]] = {}
    ex_ids = set()

    for sid, node in (per_spot or {}).items():
        microbes = (node or {}).get("microbes") or {}
        for mid, rec in microbes.items():
            if "objective" in rec:
                obj_by_microbe.setdefault(mid, {})[sid] = _safe(rec.get("objective"))
            # fluxes: prefer curated; also accept broader set if present
            fx = (rec.get("fluxes") or {})
            if not fx and "fluxes_all_exchanges" in rec:
                fx = rec.get("fluxes_all_exchanges") or {}
            for rid, val in (fx or {}).items():
                ex_ids.add(rid)
                flux_by_ex.setdefault(rid, {}).setdefault(mid, {})[sid] = _safe(val)

    return obj_by_microbe, flux_by_ex, sorted(ex_ids)

def _collect_from_eval(residuals_csv: Optional[Path]) -> Dict[str, Dict[str, float]]:
    """Residuals per (spot, microbe) from evaluate/residuals.csv."""
    df = _read_table_maybe(residuals_csv)
    out: Dict[str, Dict[str, float]] = {}
    if df is None or df.empty:
        return out
    need = {"spot_id", "microbe", "residual"}
    if not need.issubset(set(df.columns)):
        return out
    for _, r in df.iterrows():
        out.setdefault(str(r["microbe"]), {})[str(r["spot_id"])] = _safe(r["residual"])
    return out

def _collect_from_predict(pred_csv: Optional[Path]) -> Dict[str, Dict[str, float]]:
    """Predicted y_hat_mean per (spot, microbe) from predict/predictions.csv."""
    df = _read_table_maybe(pred_csv)
    out: Dict[str, Dict[str, float]] = {}
    if df is None or df.empty:
        return out
    if not {"spot_id", "microbe"}.issubset(set(df.columns)):
        return out
    val_col = "y_hat_mean" if "y_hat_mean" in df.columns else None
    if not val_col:
        return out
    for _, r in df.iterrows():
        out.setdefault(str(r["microbe"]), {})[str(r["spot_id"])] = _safe(r[val_col])
    return out

def _collect_from_whatif(whatif_path: Optional[Path], baseline_preds_df: Optional[pd.DataFrame]) -> Dict[str, Dict[str, float]]:
    """
    What-if delta per (spot, microbe).
    Accepts either:
      - legacy CSV with columns ['spot_id','microbe','delta_mean']
      - predictions CSV with columns ['spot_id','microbe','y_hat_mean']; in this case,
        if baseline_preds_df is provided and has y_hat_mean for the same keys,
        we compute delta = whatif - baseline.
    """
    df = _read_table_maybe(whatif_path)
    out: Dict[str, Dict[str, float]] = {}
    if df is None or df.empty:
        return out

    cols = set(df.columns)

    # Case A: explicit delta file
    if {"spot_id", "microbe", "delta_mean"}.issubset(cols):
        for _, r in df.iterrows():
            out.setdefault(str(r["microbe"]), {})[str(r["spot_id"])] = _safe(r["delta_mean"])
        return out

    # Case B: it's a predictions file; try to compute delta against baseline
    if {"spot_id", "microbe", "y_hat_mean"}.issubset(cols) and baseline_preds_df is not None:
        if not {"spot_id", "microbe", "y_hat_mean"}.issubset(set(baseline_preds_df.columns)):
            return out
        base = baseline_preds_df[["spot_id","microbe","y_hat_mean"]].copy()
        base["spot_id"] = base["spot_id"].astype(str)
        base["microbe"] = base["microbe"].astype(str)
        base = base.set_index(["spot_id","microbe"])["y_hat_mean"].to_dict()

        for _, r in df.iterrows():
            key = (str(r["spot_id"]), str(r["microbe"]))
            if key in base:
                delta = _safe(r["y_hat_mean"]) - _safe(base[key])
                out.setdefault(key[1], {})[key[0]] = delta
        return out

    # Otherwise, can't interpret as delta
    return out

def _build_payload(
    system_yml: Path,
    metabolism_json: Optional[Path],
    residuals_csv: Optional[Path],
    predictions_csv: Optional[Path],
    whatif_csv: Optional[Path],
) -> Dict[str, Any]:
    spots, abund_by_microbe, mets_by_id = _collect_from_system(system_yml)
    obj_by_microbe, flux_by_ex, exchange_ids = _collect_from_metabolism(metabolism_json)
    resid_by_microbe = _collect_from_eval(residuals_csv)

    # Baseline predictions table if provided (used to compute what-if deltas when needed)
    pred_df = _read_table_maybe(predictions_csv)
    yhat_by_microbe = _collect_from_predict(predictions_csv)
    delta_by_microbe = _collect_from_whatif(whatif_csv, pred_df)

    microbes = sorted(set(list(abund_by_microbe.keys()) +
                          list(obj_by_microbe.keys()) +
                          list(yhat_by_microbe.keys()) +
                          list(resid_by_microbe.keys()) +
                          list(delta_by_microbe.keys())))
    metabolites = sorted(mets_by_id.keys())
    envs = sorted({s["env_id"] for s in spots if s.get("env_id")})

    payload = {
        "created_utc": dt.datetime.utcnow().isoformat(),
        "system": str(system_yml),
        "spots": spots,               # [{id, env_id, x,y,z}]
        "envs": envs,
        "microbes": microbes,
        "metabolites": metabolites,
        "exchanges": exchange_ids,
        "layers": {
            # measured
            "abundance": abund_by_microbe,   # microbe -> {spot -> value}
            "metabolite": mets_by_id,        # metabolite -> {spot -> value}
            # metabolism
            "objective": obj_by_microbe,     # microbe -> {spot -> value}
            "flux": flux_by_ex,              # exch -> microbe -> {spot -> value}
            # model outputs
            "residual": resid_by_microbe,    # microbe -> {spot -> value}
            "yhat": yhat_by_microbe,         # microbe -> {spot -> value}
            # what-if
            "delta": delta_by_microbe,       # microbe -> {spot -> value}
        },
        "has": {
            "metabolism": bool(obj_by_microbe or flux_by_ex),
            "evaluate": bool(resid_by_microbe),
            "predict": bool(yhat_by_microbe),
            "whatif": bool(delta_by_microbe),
        }
    }
    return payload

_HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en"><head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width,initial-scale=1"/>
<title>MicroScape Viz Report</title>
<style>
  body { font-family: system-ui, -apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; margin: 0; display: grid; grid-template-columns: 320px 1fr; height: 100vh; }
  aside { padding: 16px; border-right: 1px solid #eee; overflow: auto; }
  main { padding: 8px; overflow: auto; }
  h1 { font-size: 18px; margin: 0 0 8px; }
  .section { margin-bottom: 16px; }
  label { display: block; font-size: 12px; color: #444; margin: 8px 0 4px; }
  select, button { width: 100%; padding: 6px 8px; }
  .row { display: grid; grid-template-columns: 1fr 1fr; gap: 8px; }
  .muted { color:#666; font-size: 12px; }
  .badge { display:inline-block; padding:2px 6px; border-radius: 10px; background:#eef; color:#334; font-size: 12px; margin-right:4px; }
  .footer { font-size: 12px; color: #777; margin-top: 16px; }
</style>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
</head>
<body>
<aside>
  <h1>MicroScape Viz</h1>
  <div class="section">
    <div>System: <span id="sys"></span></div>
    <div class="muted">Created: <span id="created"></span></div>
  </div>

  <div class="section">
    <label for="envSel">Environment</label>
    <select id="envSel"></select>

    <label for="layerSel">Layer</label>
    <select id="layerSel">
      <optgroup label="Measured">
        <option value="abundance">Abundance (per microbe)</option>
        <option value="metabolite">Metabolite (per compound)</option>
      </optgroup>
      <optgroup label="Metabolism">
        <option value="objective">Objective (per microbe)</option>
        <option value="flux">Exchange flux (per microbe + exchange)</option>
      </optgroup>
      <optgroup label="Model">
        <option value="yhat">Prediction ŷ (per microbe)</option>
        <option value="residual">Residual (per microbe)</option>
      </optgroup>
      <optgroup label="What-if">
        <option value="delta">Delta (per microbe)</option>
      </optgroup>
    </select>

    <div class="row">
      <div>
        <label for="microbeSel">Microbe</label>
        <select id="microbeSel"></select>
      </div>
      <div>
        <label for="exSel">Exchange</label>
        <select id="exSel"></select>
      </div>
    </div>

    <label for="metSel">Metabolite</label>
    <select id="metSel"></select>

    <div class="section">
      <button id="refreshBtn">Update</button>
      <div class="muted" style="margin-top:8px">
        <span class="badge" id="badgeMetabolism">metabolism</span>
        <span class="badge" id="badgeEvaluate">evaluate</span>
        <span class="badge" id="badgePredict">predict</span>
        <span class="badge" id="badgeWhatif">what-if</span>
      </div>
    </div>

    <div class="footer">
      Tip: choose a layer, then a microbe or metabolite. Missing values are hidden.<br/>
      Colors: Viridis. Negative (flux uptake) vs positive (secretion) shown on the same scale.<br/>
      Δ (what-if) is computed as (scenario − baseline) when both predictions are provided.
    </div>
  </div>
</aside>

<main>
  <div id="plot" style="width:100%;height:60vh;"></div>
  <div id="hist" style="width:100%;height:32vh;"></div>
</main>

<script id="DATA" type="application/json">{DATA_JSON}</script>
<script>
(function(){
  const DATA = JSON.parse(document.getElementById('DATA').textContent);

  // DOM
  const sysEl = document.getElementById('sys');
  const createdEl = document.getElementById('created');
  const envSel = document.getElementById('envSel');
  const layerSel = document.getElementById('layerSel');
  const microbeSel = document.getElementById('microbeSel');
  const exSel = document.getElementById('exSel');
  const metSel = document.getElementById('metSel');
  const refreshBtn = document.getElementById('refreshBtn');

  const bMet = document.getElementById('badgeMetabolism');
  const bEval = document.getElementById('badgeEvaluate');
  const bPred = document.getElementById('badgePredict');
  const bWI  = document.getElementById('badgeWhatif');

  sysEl.textContent = DATA.system;
  createdEl.textContent = DATA.created_utc.replace('T',' ').split('.')[0];

  // populate envs (add "all")
  envSel.innerHTML = '<option value="__ALL__">All environments</option>' + DATA.envs.map(e => `<option>${e}</option>`).join('');
  // microbes
  microbeSel.innerHTML = DATA.microbes.map(m => `<option>${m}</option>`).join('');
  // exchanges
  exSel.innerHTML = DATA.exchanges.map(r => `<option>${r}</option>`).join('');
  // metabolites
  metSel.innerHTML = DATA.metabolites.map(m => `<option>${m}</option>`).join('');

  // badges
  function setBadge(el, on){ el.style.opacity = on ? 1 : 0.3; }
  setBadge(bMet, DATA.has.metabolism);
  setBadge(bEval, DATA.has.evaluate);
  setBadge(bPred, DATA.has.predict);
  setBadge(bWI,  DATA.has.whatif);

  function getLayerValues(layer, env, microbe, exch, metabolite) {
    const spots = DATA.spots.filter(s => env==='__ALL__' || s.env_id===env);

    const result = [];
    for (const s of spots) {
      let v = null;
      if (layer === 'metabolite') {
        if (!metabolite) continue;
        const bySpot = (DATA.layers.metabolite[metabolite] || {});
        v = bySpot[s.id];
      } else if (layer === 'flux') {
        if (!exch || !microbe) continue;
        const exBlock = DATA.layers.flux[exch] || {};
        const bySpot = exBlock[microbe] || {};
        v = bySpot[s.id];
      } else {
        // layers that are microbe->spot maps
        if (!microbe) continue;
        const block = (DATA.layers[layer] || {})[microbe] || {};
        v = block[s.id];
      }
      if (v !== null && v !== undefined) {
        result.push({x: s.x, y: s.y, z: s.z, spot: s.id, env: s.env_id, val: v});
      }
    }
    return result;
  }

  function render(){
    const env = envSel.value || '__ALL__';
    const layer = layerSel.value;
    const microbe = microbeSel.value || null;
    const exch = exSel.value || null;
    const metabolite = metSel.value || null;

    // enable/disable controls based on layer
    microbeSel.disabled = (layer === 'metabolite');
    exSel.disabled = (layer !== 'flux');
    metSel.disabled = (layer !== 'metabolite');

    const rows = getLayerValues(layer, env, microbe, exch, metabolite);
    const xs = rows.map(r=>r.x), ys = rows.map(r=>r.y), cs = rows.map(r=>r.val);
    const text = rows.map(r=>`${r.spot}<br>env=${r.env}<br>value=${r.val}`);

    const titleMap = {
      'abundance': `Abundance (microbe=${microbe})`,
      'metabolite': `Metabolite ${metabolite}`,
      'objective': `Objective (microbe=${microbe})`,
      'flux': `Flux ${exch} (microbe=${microbe})`,
      'yhat': `Prediction ŷ (microbe=${microbe})`,
      'residual': `Residual (microbe=${microbe})`,
      'delta': `What-if Δ (microbe=${microbe})`,
    };
    const title = titleMap[layer] || layer;

    const plotDiv = document.getElementById('plot');
    const histDiv = document.getElementById('hist');

    // Scatter (2D XY)
    const scatter = [{
      x: xs, y: ys, text: text, mode: 'markers',
      type: 'scattergl',
      marker: { size: 8, color: cs, colorscale: 'Viridis', showscale: true, colorbar: {title: 'value'} },
      hovertemplate: '%{text}<extra></extra>'
    }];
    Plotly.newPlot(plotDiv, scatter, {
      title, xaxis:{title:'x (µm)', zeroline:false}, yaxis:{title:'y (µm)', zeroline:false, scaleanchor:'x', scaleratio:1},
      margin:{l:48,r:24,t:36,b:48}, paper_bgcolor:'#fff', plot_bgcolor:'#fff'
    }, {responsive:true, displaylogo:false});

    // Histogram
    const hist = [{
      x: cs, type:'histogram', nbinsx: 40
    }];
    Plotly.newPlot(histDiv, hist, {
      title: 'Value distribution', xaxis:{title:'value'}, yaxis:{title:'count'},
      margin:{l:48,r:24,t:36,b:48}, paper_bgcolor:'#fff', plot_bgcolor:'#fff'
    }, {responsive:true, displaylogo:false});
  }

  refreshBtn.addEventListener('click', render);
  layerSel.addEventListener('change', render);
  envSel.addEventListener('change', render);
  microbeSel.addEventListener('change', render);
  exSel.addEventListener('change', render);
  metSel.addEventListener('change', render);

  // initial defaults: pick first values if present
  if (DATA.microbes.length) microbeSel.value = DATA.microbes[0];
  if (DATA.metabolites.length) metSel.value = DATA.metabolites[0];
  if (DATA.exchanges.length) exSel.value = DATA.exchanges[0];
  render();
})();
</script>
</body></html>
"""

@app.command("viz")
def viz_cmd(
    system_yml: Path = typer.Argument(..., help="Path to system.yml"),
    outdir: Path = typer.Option("outputs/viz", help="Output folder for the report"),
    # optional inputs (pass what you have)
    metabolism_json: Optional[Path] = typer.Option(None, "--metabolism-json", help="metabolism_* .json to show objective/fluxes"),
    eval_residuals: Optional[Path] = typer.Option(None, "--residuals", help="Path to evaluate/residuals.csv"),
    predictions: Optional[Path] = typer.Option(None, "--predictions", help="Path to predict/predictions.csv"),
    whatif_csv: Optional[Path] = typer.Option(None, "--whatif", help="Path to what-if output: either delta CSV or what-if predictions CSV"),
    separate_data: bool = typer.Option(False, "--separate-data", help="Also write viz_data.json next to HTML (HTML still embeds data)."),
    filename: str = typer.Option("viz_report.html", help="Output HTML filename"),
):
    """
    Build a self-contained interactive HTML report to explore spatial layers:
      - Measured: microbe abundance (per microbe), metabolites (per compound)
      - Metabolism: objective (per microbe), exchange fluxes (per microbe + exchange)
      - Model: residuals (per microbe)
      - Predict: y_hat (per microbe)
      - What-if: delta (per microbe). If --whatif is a what-if predictions CSV,
                 and --predictions is also given as baseline predictions, Δ is computed as (scenario − baseline).
    """
    outdir = Path(outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    payload = _build_payload(system_yml, metabolism_json, eval_residuals, predictions, whatif_csv)

    # Write HTML (embed JSON)
    html = _HTML_TEMPLATE.replace("{DATA_JSON}", json.dumps(payload, separators=(",", ":")))
    html_path = outdir / filename
    html_path.write_text(html)

    if separate_data:
        (outdir / "viz_data.json").write_text(json.dumps(payload, indent=2))

    typer.echo(f"Viz report: {html_path}")
    if separate_data:
        typer.echo(f"Data JSON : {outdir / 'viz_data.json'}")
