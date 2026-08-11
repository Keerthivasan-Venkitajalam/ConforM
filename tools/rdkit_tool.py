"""Ligand validation and preparation using RDKit (real, installed via conda-forge)."""
from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, QED, rdMolDescriptors


@dataclass
class ValidatedLigand:
    name: str
    smiles_in: str
    canonical_smiles: str | None
    valid: bool
    error: str | None
    mol_weight: float | None = None
    logp: float | None = None
    qed: float | None = None
    lipinski_violations: int | None = None
    num_rotatable_bonds: int | None = None
    embedded_3d: bool = False
    pdb_path: str | None = None


def lipinski_violations(mol) -> int:
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    return violations


def validate_and_prepare(name: str, smiles: str, out_dir=None) -> ValidatedLigand:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ValidatedLigand(name=name, smiles_in=smiles, canonical_smiles=None,
                                valid=False, error="RDKit could not parse/sanitize SMILES")
    try:
        Chem.SanitizeMol(mol)
    except Exception as exc:  # noqa: BLE001
        return ValidatedLigand(name=name, smiles_in=smiles, canonical_smiles=None,
                                valid=False, error=f"Sanitization failed: {exc}")

    canonical = Chem.MolToSmiles(mol)
    result = ValidatedLigand(
        name=name, smiles_in=smiles, canonical_smiles=canonical, valid=True, error=None,
        mol_weight=Descriptors.MolWt(mol),
        logp=Descriptors.MolLogP(mol),
        qed=QED.qed(mol),
        lipinski_violations=lipinski_violations(mol),
        num_rotatable_bonds=rdMolDescriptors.CalcNumRotatableBonds(mol),
    )

    if out_dir is not None:
        mol_h = Chem.AddHs(mol)
        embed_status = AllChem.EmbedMolecule(mol_h, randomSeed=42, useRandomCoords=True)
        if embed_status == 0:
            AllChem.MMFFOptimizeMolecule(mol_h, maxIters=500)
            out_dir.mkdir(parents=True, exist_ok=True)
            pdb_path = out_dir / f"{name}.pdb"
            Chem.MolToPDBFile(mol_h, str(pdb_path))
            result.embedded_3d = True
            result.pdb_path = str(pdb_path)
        else:
            result.error = "3D embedding failed (ETKDG could not find valid conformer)"

    return result


if __name__ == "__main__":
    import csv
    import sys
    from pathlib import Path

    with open(sys.argv[1]) as f:
        reader = csv.DictReader(f)
        for row in reader:
            r = validate_and_prepare(row["name"], row["smiles"], Path("artifacts/ligands"))
            status = "OK" if r.valid and r.embedded_3d else ("INVALID" if not r.valid else "NO-3D")
            print(f"{status:8s} {row['name']:45s} err={r.error}")
