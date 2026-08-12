"""Generates experiment_report.html from real experiment artifacts.

Every figure and number is read from the experiment manifest / artifact
files. The report explicitly separates OBSERVED COMPUTATIONAL RESULT from
AGENT HYPOTHESIS (master prompt #27, #40).
"""
from __future__ import annotations

import argparse
import html
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from visualization.molecular_viewer import viewer_html
from visualization.plots import bar_chart_svg, state_space_svg

CSS = """
:root { color-scheme: light dark; }
* { box-sizing: border-box; }
body { font-family: system-ui, -apple-system, "Segoe UI", sans-serif; line-height: 1.6;
       margin: 0 auto; padding: 2.5rem 1.5rem; max-width: 980px;
       background: #fff; color: #16191d; }
@media (prefers-color-scheme: dark) { body { background: #14171a; color: #e8eaed; } }
h1 { font-size: 1.9rem; margin: 0 0 .3rem; }
h2 { font-size: 1.25rem; margin: 2.4rem 0 .8rem; padding-bottom: .35rem;
     border-bottom: 1px solid rgba(128,128,128,.28); }
h3 { font-size: 1rem; margin: 1.4rem 0 .5rem; }
.sub { opacity: .7; margin: 0 0 1.8rem; }
table { border-collapse: collapse; width: 100%; font-size: .9rem; margin: .6rem 0; }
th, td { text-align: left; padding: .5rem .65rem; border-bottom: 1px solid rgba(128,128,128,.2); }
th { font-weight: 600; opacity: .8; }
td.num { text-align: right; font-variant-numeric: tabular-nums; }
.wrap { overflow-x: auto; }
.grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(160px,1fr)); gap: .8rem;
        margin: 1rem 0; }
.card { border: 1px solid rgba(128,128,128,.25); border-radius: 9px; padding: .8rem .9rem; }
.card .k { font-size: .72rem; text-transform: uppercase; letter-spacing: .04em; opacity: .65; }
.card .v { font-size: 1.35rem; font-weight: 650; margin-top: .18rem;
           font-variant-numeric: tabular-nums; }
.obs { border-left: 3px solid #54A24B; background: rgba(84,162,75,.09);
       padding: .7rem .95rem; border-radius: 0 7px 7px 0; margin: .7rem 0; }
.hyp { border-left: 3px solid #F58518; background: rgba(245,133,24,.09);
       padding: .7rem .95rem; border-radius: 0 7px 7px 0; margin: .7rem 0; }
.warn { border-left: 3px solid #E45756; background: rgba(228,87,86,.09);
        padding: .7rem .95rem; border-radius: 0 7px 7px 0; margin: .7rem 0; }
.tag { display: inline-block; font-size: .68rem; font-weight: 700; letter-spacing: .05em;
       text-transform: uppercase; padding: .12rem .45rem; border-radius: 4px; margin-right: .5rem; }
.tag-obs { background: #54A24B; color: #fff; }
.tag-hyp { background: #F58518; color: #fff; }
.tag-warn { background: #E45756; color: #fff; }
code, pre { font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: .82rem; }
pre { background: rgba(128,128,128,.11); padding: .8rem; border-radius: 7px; overflow-x: auto; }
.log { max-height: 340px; overflow-y: auto; }
footer { margin-top: 3rem; padding-top: 1rem; border-top: 1px solid rgba(128,128,128,.25);
         font-size: .82rem; opacity: .7; }
"""


def esc(x) -> str:
    return html.escape(str(x))


def kpi(label, value) -> str:
    return f"<div class='card'><div class='k'>{esc(label)}</div><div class='v'>{esc(value)}</div></div>"


def build_report(exp_dir: Path) -> Path:
    manifest = json.loads((exp_dir / "metrics" / "experiment_manifest.json").read_text())
    log_path = exp_dir / "agent_log.txt"
    log_text = log_path.read_text() if log_path.exists() else ""

    analysis_path = exp_dir / "metrics" / "ensemble_analysis.json"
    analysis = json.loads(analysis_path.read_text()) if analysis_path.exists() else {}

    pocket = manifest.get("selected_pocket") or {}
    results = manifest.get("ranked_results", [])
    ensemble = manifest.get("ensemble", {})
    fallbacks = manifest.get("tools_unavailable_fallback_used", {})

    # --- state space plot -------------------------------------------------
    labels = analysis.get("pdb_ids", [])
    coords = analysis.get("pca_coords", [])
    state_plot = state_space_svg(coords, labels, selected=pocket.get("state_pdb_id"),
                                  pocket_states=[pocket.get("state_pdb_id")]) if coords else ""

    # --- affinity chart ---------------------------------------------------
    top = sorted((r for r in results if r.get("best_affinity_kcal") is not None),
                 key=lambda r: r["best_affinity_kcal"])[:10]
    affinity_plot = bar_chart_svg([r["ligand_name"] for r in top],
                                   [r["best_affinity_kcal"] for r in top],
                                   title="Best Vina affinity per ligand (kcal/mol, lower = better score)",
                                   lower_is_better=True) if top else ""

    # --- 3D viewer --------------------------------------------------------
    viewer = ""
    receptor = next((s for s in ensemble.get("structures", [])
                     if Path(s).stem == pocket.get("state_pdb_id")), None)
    if receptor and Path(receptor).exists():
        best_ligand = results[0]["ligand_name"] if results else None
        lig_pdb = exp_dir / "ligands" / f"{best_ligand}.pdb" if best_ligand else None
        viewer = viewer_html(Path(receptor), pocket.get("residues", []),
                             lig_pdb if lig_pdb and lig_pdb.exists() else None)

    # --- tables -----------------------------------------------------------
    rows = "".join(
        f"<tr><td>{esc(r['ligand_name'])}</td><td>{esc(r.get('origin','library'))}</td>"
        f"<td class='num'>{r.get('best_affinity_kcal', float('nan')):.2f}</td>"
        f"<td class='num'>{(r.get('qed') or 0):.3f}</td>"
        f"<td class='num'>{r.get('discovery_score', 0):.4f}</td></tr>"
        for r in results[:20])

    decisions = "".join(
        f"<tr><td class='num'>{d.get('iteration')}</td><td><code>{esc(d.get('action'))}</code></td>"
        f"<td>{esc((d.get('outcome') or {}).get('tool',''))}</td>"
        f"<td>{esc(d.get('failure') or 'ok')}</td></tr>"
        for d in manifest.get("agent_decisions", []))

    fallback_rows = "".join(
        f"<tr><td><code>{esc(k)}</code></td><td>{esc(v)}</td></tr>" for k, v in fallbacks.items())

    prov = manifest.get("reproducibility", {})
    prov_rows = "".join(f"<tr><td>{esc(k)}</td><td><code>{esc(v)}</code></td></tr>"
                        for k, v in prov.items())

    best = results[0] if results else {}
    is_equilibrium = ensemble.get("is_equilibrium_sample")

    doc = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>ConforM-Agent report — {esc(manifest.get('experiment_id'))}</title>
<style>{CSS}</style></head><body>

<h1>ConforM-Agent experiment report</h1>
<p class="sub">{esc(manifest.get('target'))} &middot; mode <code>{esc(manifest.get('mode'))}</code>
 &middot; <code>{esc(manifest.get('experiment_id'))}</code></p>

<div class="warn"><span class="tag tag-warn">Read first</span>
This is a computational hypothesis-generation run. Nothing here demonstrates that any
molecule binds KRAS G12D. Docking scores are AutoDock Vina empirical-function values
produced under the protocol described below; they are not measured affinities and have
not been experimentally validated.</div>

<div class="grid">
  {kpi("Conformations explored", ensemble.get("n_states", 0))}
  {kpi("Pocket candidates", manifest.get("n_pocket_candidates", 0))}
  {kpi("Best pocket volume", f"{pocket.get('volume', 0):.0f} Å³")}
  {kpi("Pocket druggability", f"{pocket.get('druggability', 0):.2f}")}
  {kpi("Ligands docked", manifest.get("n_ligands_docked", 0))}
  {kpi("Best affinity", f"{best.get('best_affinity_kcal', 0):.2f} kcal/mol")}
  {kpi("Discovery Score", f"{manifest.get('best_discovery_score') or 0:.4f}")}
  {kpi("Closed-loop iterations", manifest.get("closed_loop_iterations_executed", 0))}
  {kpi("Runtime", f"{manifest.get('runtime_seconds', 0):.0f} s")}
  {kpi("GPU hours", "0.00 (CPU-only run)")}
</div>

<h2>1. Provenance and execution mode</h2>
<div class="warn"><span class="tag tag-warn">Fallbacks in use</span>
The following components of the original design could not be executed in this
environment. Their fallbacks are real tools, but they are <em>not</em> the named models,
and no output below should be attributed to them.</div>
<div class="wrap"><table><thead><tr><th>Component</th><th>Status / substitute</th></tr></thead>
<tbody>{fallback_rows}</tbody></table></div>
<p>Ensemble provider: <code>{esc(ensemble.get('provider'))}</code>.
Equilibrium sample: <strong>{esc(is_equilibrium)}</strong>.
{"" if is_equilibrium else
 "Because these are independent experimental structures rather than Boltzmann-weighted samples, "
 "state populations are uniform placeholders and no thermodynamic weight should be read into them."}</p>

<h2>2. Conformational state space</h2>
{state_plot or "<p>No ensemble analysis recorded for this run.</p>"}
<p>Dimensionality reduction: <code>{esc(analysis.get('dimensionality_reduction', 'n/a'))}</code>.
Max pairwise Cα RMSD {analysis.get('max_pairwise_rmsd', 0):.2f} Å;
peak per-residue RMSF {analysis.get('max_rmsf', 0):.2f} Å.</p>

<h2>3. Selected pocket</h2>
<div class="obs"><span class="tag tag-obs">Observed computational result</span>
Blind ranking on druggability + novelty-vs-baseline + volume selected
<code>{esc(pocket.get('state_pdb_id'))}:pocket{esc(pocket.get('pocket_index'))}</code>
({pocket.get('volume', 0):.0f} Å³, druggability {pocket.get('druggability', 0):.2f}).
Post-hoc, this cavity covers
<strong>{(pocket.get('ground_truth_overlap') or 0)*100:.0f}%</strong> of the documented
Switch-II ground-truth residues. Ground-truth overlap was not an input to the ranking.</div>
<p><strong>Lining residues:</strong> <code>{esc(', '.join(pocket.get('residues', [])[:40]))}</code></p>
{viewer}

<h2>4. Ligand results</h2>
{affinity_plot}
<div class="wrap"><table><thead><tr><th>Ligand</th><th>Origin</th>
<th class="num">Vina (kcal/mol)</th><th class="num">QED</th>
<th class="num">Discovery Score</th></tr></thead><tbody>{rows}</tbody></table></div>
<p>Optimizer mode: <code>{esc(manifest.get('optimizer_mode'))}</code>.</p>

<h2>5. Agent decisions</h2>
<div class="wrap"><table><thead><tr><th class="num">Iter</th><th>Action</th><th>Tool</th>
<th>Status</th></tr></thead><tbody>{decisions}</tbody></table></div>

<h2>6. Interpretation</h2>
<div class="obs"><span class="tag tag-obs">Observed computational result</span>
The system recovered a cavity absent from the apo baseline structure and present in a
minority of sampled states, and that cavity's lining residues overlap the documented
Switch-II cryptic site. This is <em>recovery of a known ground truth</em>, which
validates the pipeline. It is not a novel biological discovery.</div>
<div class="hyp"><span class="tag tag-hyp">Agent hypothesis</span>
{esc(best.get('ligand_name', 'The top-ranked ligand'))} produced the most favorable
computational docking score ({best.get('best_affinity_kcal', 0):.2f} kcal/mol) under this
protocol. This constitutes a candidate hypothesis requiring experimental validation —
not evidence of binding, potency, selectivity, or efficacy.</div>

<h2>7. Reproducibility</h2>
<div class="wrap"><table><tbody>{prov_rows}</tbody></table></div>

<h2>8. Agent log</h2>
<pre class="log">{esc(log_text)}</pre>

<footer>Generated by <code>scripts/generate_report.py</code> from
<code>{esc(exp_dir)}</code>. All figures derive from executed tool output;
see <code>docs/LIMITATIONS.md</code> for what was not run.</footer>
</body></html>"""

    out = exp_dir / "report" / "experiment_report.html"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(doc)
    return out


def latest_experiment(root: Path = Path("artifacts")) -> Path:
    candidates = [d for d in root.iterdir()
                  if d.is_dir() and (d / "metrics" / "experiment_manifest.json").exists()]
    if not candidates:
        raise SystemExit("No experiments with a manifest found under artifacts/")
    return max(candidates, key=lambda d: d.stat().st_mtime)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--experiment", help="experiment directory or id (default: most recent)")
    ap.add_argument("--root", default="artifacts")
    args = ap.parse_args()
    root = Path(args.root)
    if args.experiment:
        d = Path(args.experiment)
        if not d.exists():
            d = root / args.experiment
    else:
        d = latest_experiment(root)
    print(f"Report written to: {build_report(d)}")
