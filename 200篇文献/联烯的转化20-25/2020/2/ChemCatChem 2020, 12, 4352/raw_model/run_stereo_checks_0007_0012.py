#!/usr/bin/env python3
"""Persist standard stereo graph diagnostics for page 7--12 source molecules."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
from standard_records import stereo_inspect  # noqa: E402

bundle = json.loads((ROOT / "runs/current/submitted.json").read_text())
page_evidence = {
    f"doc_2401da5b5fbe4542:page{page}:{kind}"
    for page in range(7, 13) for kind in ("text", "image")
}
results = []
for molecule in bundle["molecules"]:
    if not any(page_evidence.intersection(occurrence["evidence_ids"])
               for occurrence in molecule["occurrences"]):
        continue
    smiles = molecule["smiles"]["normalized_value"]
    item = {
        "molecule_id": molecule["molecule_id"],
        "source_evidence_ids": sorted(set().union(*(set(o["evidence_ids"]) for o in molecule["occurrences"])) & page_evidence),
        "registry_status": {
            "chirality_status": molecule["stereochemistry"]["chirality_status"]["normalized_value"],
            "geometry_status": molecule["stereochemistry"]["geometry_status"]["normalized_value"],
            "representation_status": molecule["stereochemistry"]["representation_status"]["normalized_value"],
        },
        "smiles": smiles,
    }
    if smiles:
        item["standard_stereo_diagnostic"] = stereo_inspect(smiles)
    else:
        item["standard_stereo_diagnostic"] = None
        item["limitation"] = "No safely encoded structure; source stereochemical statement/drawing is retained in the molecule record."
    results.append(item)

out = {
    "document_id": "doc_2401da5b5fbe4542",
    "page_range": [7, 12],
    "source_verified": False,
    "molecule_count": len(results),
    "checks": results,
    "limitations": "Toolkit graph diagnostics do not establish source identity, absolute configuration, ee, dr, E/Z, ring cis/trans, or axial assignment.",
}
(ROOT / "runs/current/raw_model/stereo_checks_pages_0007_0012.json").write_text(json.dumps(out, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"molecule_count": len(results), "encoded": sum(x["smiles"] is not None for x in results), "unencoded": sum(x["smiles"] is None for x in results)}))
