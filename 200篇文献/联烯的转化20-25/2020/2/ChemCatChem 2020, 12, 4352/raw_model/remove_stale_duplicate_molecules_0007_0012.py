#!/usr/bin/env python3
"""Remove three superseded machine-only registry duplicates before re-merge.

The standard merge deliberately preserves absent old molecules.  These three
IDs were created by the first page-7/10 pass before cross-page deduplication;
their canonical records are mol_19, mol_pmosnco, and mol_pnosnco.  No
source-reviewed fact is removed.
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
path = ROOT / "runs/current/dataset.json"
data = json.loads(path.read_text())
stale = {"mol_benzoyl_acetonitrile", "mol_s9", "mol_s10"}
standard = data.get("standard", {})
removed = []
kept = []
for molecule in standard.get("molecules", []):
    if molecule["molecule_id"] not in stale:
        kept.append(molecule)
        continue
    source_checked = any(
        node.get("verification") == "source_checked"
        for node in [molecule, *molecule.values()]
        if isinstance(node, dict)
    )
    if source_checked:
        raise RuntimeError(f"Refusing to remove source-reviewed molecule {molecule['molecule_id']}")
    removed.append(molecule["molecule_id"])
standard["molecules"] = kept
path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"removed_machine_only_duplicate_ids": removed}, ensure_ascii=False))
