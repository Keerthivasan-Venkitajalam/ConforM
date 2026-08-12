"""Streamlit dashboard.

Consumes REAL experiment artifacts from artifacts/<experiment_id>/. There are
no hard-coded demo events: if no experiment has been run, the dashboard says
so rather than showing invented output.

Run:  streamlit run visualization/dashboard.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import streamlit as st
import streamlit.components.v1 as components

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from visualization.molecular_viewer import viewer_html
from visualization.plots import bar_chart_svg, state_space_svg

ARTIFACTS = Path("artifacts")


def list_experiments() -> list[Path]:
    if not ARTIFACTS.exists():
        return []
    dirs = [d for d in ARTIFACTS.iterdir()
            if d.is_dir() and (d / "metrics" / "experiment_manifest.json").exists()]
    return sorted(dirs, key=lambda d: d.stat().st_mtime, reverse=True)


st.set_page_config(page_title="ConforM-Agent", layout="wide", page_icon="🧬")
st.title("ConforM-Agent")
st.caption("Closed-Loop Molecular State Space Explorer")

experiments = list_experiments()
if not experiments:
    st.warning(
        "No experiments found. Run one first:\n\n"
        "```bash\npython scripts/run_experiment.py --target kras-g12d\n```\n\n"
        "This dashboard only displays real executed-experiment artifacts; "
        "it does not ship canned demo output."
    )
    st.stop()

choice = st.sidebar.selectbox("Experiment", experiments, format_func=lambda d: d.name)
manifest = json.loads((choice / "metrics" / "experiment_manifest.json").read_text())
analysis_path = choice / "metrics" / "ensemble_analysis.json"
analysis = json.loads(analysis_path.read_text()) if analysis_path.exists() else {}

pocket = manifest.get("selected_pocket") or {}
ensemble = manifest.get("ensemble", {})
results = manifest.get("ranked_results", [])
best = results[0] if results else {}

st.error(
    "**Scientific honesty notice** — docking scores below are AutoDock Vina empirical "
    "values under this protocol. They are not measured binding affinities and do not "
    "show that any molecule binds KRAS G12D.", icon="⚠️")

# --- metrics row ---------------------------------------------------------
cols = st.columns(5)
cols[0].metric("Conformations explored", ensemble.get("n_states", 0))
cols[1].metric("Pocket candidates", manifest.get("n_pocket_candidates", 0))
cols[2].metric("Best pocket volume", f"{pocket.get('volume', 0):.0f} Å³")
cols[3].metric("Best affinity", f"{best.get('best_affinity_kcal', 0):.2f} kcal/mol")
cols[4].metric("Discovery Score", f"{manifest.get('best_discovery_score') or 0:.3f}")

cols = st.columns(5)
cols[0].metric("Pocket druggability", f"{pocket.get('druggability', 0):.2f}")
cols[1].metric("Ligands docked", manifest.get("n_ligands_docked", 0))
cols[2].metric("Closed-loop iterations", manifest.get("closed_loop_iterations_executed", 0))
cols[3].metric("GPU time", "0.00 h (CPU-only)")
cols[4].metric("Runtime", f"{manifest.get('runtime_seconds', 0):.0f} s")

# --- fallback banner -----------------------------------------------------
fallbacks = manifest.get("tools_unavailable_fallback_used", {})
with st.expander("Execution mode / fallbacks actually used", expanded=True):
    st.write(f"Ensemble provider: `{ensemble.get('provider')}` — "
             f"equilibrium sample: **{ensemble.get('is_equilibrium_sample')}**")
    st.write(f"Optimizer mode: `{manifest.get('optimizer_mode')}`")
    for k, v in fallbacks.items():
        st.write(f"- **{k}**: {v}")

left, right = st.columns([1, 1])

with left:
    st.subheader("Agent event log")
    log_path = choice / "agent_log.txt"
    st.code(log_path.read_text() if log_path.exists() else "(no log)", language="text")

    st.subheader("Agent decisions")
    st.dataframe([
        {"iter": d.get("iteration"), "action": d.get("action"),
         "tool": (d.get("outcome") or {}).get("tool", ""),
         "status": d.get("failure") or "ok"}
        for d in manifest.get("agent_decisions", [])
    ], use_container_width=True)

with right:
    st.subheader("Structure / pocket / ligand")
    receptor = next((s for s in ensemble.get("structures", [])
                     if Path(s).stem == pocket.get("state_pdb_id")), None)
    if receptor and Path(receptor).exists():
        lig = choice / "ligands" / f"{best.get('ligand_name')}.pdb"
        components.html(
            viewer_html(Path(receptor), pocket.get("residues", []),
                        lig if lig.exists() else None, height=430),
            height=450)
        st.caption("Grey cartoon: receptor. Gold: Switch I. Magenta: Switch II. "
                   "Red surface: selected cavity. Cyan sticks: top-scoring ligand pose.")
    else:
        st.info("Receptor structure artifact not found for this experiment.")

    st.subheader("Conformational state space")
    if analysis.get("pca_coords"):
        st.markdown(state_space_svg(analysis["pca_coords"], analysis.get("pdb_ids", []),
                                     selected=pocket.get("state_pdb_id"),
                                     pocket_states=[pocket.get("state_pdb_id")]),
                    unsafe_allow_html=True)
        st.caption(analysis.get("dimensionality_reduction", ""))
    else:
        st.info("No ensemble analysis recorded for this run.")

st.subheader("Ligand results")
top = sorted((r for r in results if r.get("best_affinity_kcal") is not None),
             key=lambda r: r["best_affinity_kcal"])[:12]
if top:
    st.markdown(bar_chart_svg([r["ligand_name"] for r in top],
                               [r["best_affinity_kcal"] for r in top],
                               title="Best Vina affinity (kcal/mol, lower = better score)",
                               lower_is_better=True), unsafe_allow_html=True)
st.dataframe([
    {"ligand": r["ligand_name"], "origin": r.get("origin", "library"),
     "vina_kcal": r.get("best_affinity_kcal"), "QED": r.get("qed"),
     "discovery_score": r.get("discovery_score")}
    for r in results
], use_container_width=True)
