"""Ligand optimization.

Fallback hierarchy (master prompt #39):
    REINVENT4 -> RDKit scaffold/R-group mutation -> existing library only

REINVENT4 (MolecularAI/REINVENT4, Apache 2.0) is an RL framework whose prior
models are distributed as separate checkpoints and whose RL loop is designed
around GPU execution. It is NOT installed or executed in this build -- see
docs/LIMITATIONS.md. `ReinventOptimizer` therefore raises, and the pipeline
falls back to `RDKitMutationOptimizer`, a deterministic R-group enumeration.

The mode actually used is reported everywhere as `optimizer_mode` so the UI
and reports can never present RDKit enumeration as REINVENT RL output.
"""
from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors, QED, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

# Deterministic R-group vocabulary applied to aromatic C-H positions.
R_GROUPS = [
    ("F", "[cH1:1]>>[c:1]F"),
    ("Cl", "[cH1:1]>>[c:1]Cl"),
    ("CH3", "[cH1:1]>>[c:1]C"),
    ("OH", "[cH1:1]>>[c:1]O"),
    ("OMe", "[cH1:1]>>[c:1]OC"),
    ("CN", "[cH1:1]>>[c:1]C#N"),
    ("NH2", "[cH1:1]>>[c:1]N"),
]


@dataclass
class Analog:
    name: str
    smiles: str
    parent: str
    modification: str
    qed: float
    mol_weight: float
    lipinski_violations: int


class ReinventOptimizer:
    mode = "reinvent4"

    def generate(self, *_args, **_kwargs):
        raise NotImplementedError(
            "REINVENT4 is not installed in this environment (no prior checkpoints, "
            "RL loop is GPU-oriented). See docs/LIMITATIONS.md."
        )


class RDKitMutationOptimizer:
    """Deterministic R-group enumeration on the seed scaffold.

    This is combinatorial chemistry enumeration, NOT reinforcement learning.
    It proposes analogs; the docking engine (not this module) decides whether
    any of them score better.
    """

    mode = "rdkit_rgroup_enumeration_fallback"

    def __init__(self, max_analogs: int = 12, max_mol_weight: float = 600.0,
                 max_lipinski_violations: int = 1):
        self.max_analogs = max_analogs
        self.max_mol_weight = max_mol_weight
        self.max_lipinski_violations = max_lipinski_violations

    def generate(self, seed_name: str, seed_smiles: str) -> list[Analog]:
        seed = Chem.MolFromSmiles(seed_smiles)
        if seed is None:
            raise ValueError(f"Seed SMILES for {seed_name} could not be parsed by RDKit")

        seen: set[str] = {Chem.MolToSmiles(seed)}
        analogs: list[Analog] = []

        for label, smarts in R_GROUPS:
            rxn = AllChem.ReactionFromSmarts(smarts)
            products = rxn.RunReactants((seed,))
            for i, (product,) in enumerate(products):
                try:
                    Chem.SanitizeMol(product)
                except Exception:  # noqa: BLE001
                    continue
                smiles = Chem.MolToSmiles(product)
                if smiles in seen:
                    continue
                mw = Descriptors.MolWt(product)
                violations = sum([
                    mw > 500,
                    Descriptors.MolLogP(product) > 5,
                    rdMolDescriptors.CalcNumHBD(product) > 5,
                    rdMolDescriptors.CalcNumHBA(product) > 10,
                ])
                if mw > self.max_mol_weight or violations > self.max_lipinski_violations:
                    continue
                seen.add(smiles)
                analogs.append(Analog(
                    name=f"{seed_name}_{label}_{i}", smiles=smiles, parent=seed_name,
                    modification=f"aromatic C-H -> {label} (site {i})",
                    qed=QED.qed(product), mol_weight=mw, lipinski_violations=violations,
                ))

        # Deterministic selection: highest QED first, tie-broken by name.
        analogs.sort(key=lambda a: (-a.qed, a.name))
        return analogs[: self.max_analogs]


def get_optimizer(prefer_reinvent: bool = False):
    """Returns (optimizer, mode). Falls back to RDKit enumeration if REINVENT is unavailable."""
    if prefer_reinvent:
        try:
            opt = ReinventOptimizer()
            opt.generate("probe", "CCO")
            return opt, opt.mode
        except NotImplementedError:
            pass
    opt = RDKitMutationOptimizer()
    return opt, opt.mode


if __name__ == "__main__":
    optimizer, mode = get_optimizer(prefer_reinvent=True)
    print(f"optimizer mode: {mode}")
    for a in optimizer.generate("seed", "COc1cc(N2CCN(C)CC2)c(F)cc1Nc1ncc(Cl)cn1"):
        print(f"  {a.name:28s} QED={a.qed:.3f} MW={a.mol_weight:.1f} {a.smiles}")
