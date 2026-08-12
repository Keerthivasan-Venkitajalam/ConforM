"""py3Dmol viewers for the dashboard and HTML report.

Renders the real receptor PDB and the real docked-pose coordinates from the
experiment artifacts -- never a stand-in structure.
"""
from __future__ import annotations

from pathlib import Path

SWITCH_I = "30-40"
SWITCH_II = "60-76"


def viewer_html(receptor_pdb: Path, pocket_residues: list[str] | None = None,
                ligand_pdb: Path | None = None, height: int = 460,
                div_id: str = "conform_viewer") -> str:
    """Self-contained py3Dmol/3Dmol.js HTML fragment.

    Requires network access to load 3Dmol.js from a CDN when rendered in a
    browser; the Streamlit dashboard embeds this directly.
    """
    receptor_text = Path(receptor_pdb).read_text() if Path(receptor_pdb).exists() else ""
    ligand_text = (Path(ligand_pdb).read_text()
                   if ligand_pdb and Path(ligand_pdb).exists() else "")
    resnums = sorted({int("".join(filter(str.isdigit, r))) for r in (pocket_residues or [])
                      if any(ch.isdigit() for ch in r)})

    return f"""
<div id="{div_id}" style="width:100%;height:{height}px;position:relative;"></div>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<script>
(function() {{
  var el = document.getElementById("{div_id}");
  if (!el || typeof $3Dmol === "undefined") return;
  var viewer = $3Dmol.createViewer(el, {{backgroundColor: "white"}});
  viewer.addModel({receptor_text!r}, "pdb");
  viewer.setStyle({{}}, {{cartoon: {{color: "lightgrey", opacity: 0.85}}}});
  viewer.setStyle({{resi: "{SWITCH_I}"}}, {{cartoon: {{color: "gold"}}}});
  viewer.setStyle({{resi: "{SWITCH_II}"}}, {{cartoon: {{color: "magenta"}}}});
  var pocketResi = {resnums};
  if (pocketResi.length) {{
    viewer.addStyle({{resi: pocketResi}},
                    {{stick: {{colorscheme: "redCarbon", radius: 0.18}}}});
    viewer.addSurface($3Dmol.SurfaceType.VDW,
                      {{opacity: 0.72, color: "red"}}, {{resi: pocketResi}});
  }}
  var ligText = {ligand_text!r};
  if (ligText.length) {{
    viewer.addModel(ligText, "pdb");
    viewer.setStyle({{model: -1}}, {{stick: {{colorscheme: "cyanCarbon", radius: 0.24}}}});
  }}
  viewer.zoomTo();
  viewer.render();
}})();
</script>
"""
