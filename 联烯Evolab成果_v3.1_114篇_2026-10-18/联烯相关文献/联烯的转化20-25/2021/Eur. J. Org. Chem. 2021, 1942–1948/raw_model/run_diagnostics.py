#!/usr/bin/env python3
"""Run required literal-mention and standard-stereo diagnostics as raw artifacts."""
import json
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
PYTHON = "/EvoLab/student-workspaces/xllin/OpenEvo-chem-integration/.venv/bin/python"
TOOL = ROOT / "skills/chemical-reaction-extraction/scripts/reaction_tools.py"
DATASET = ROOT / "runs/current/dataset.json"
SUBMITTED = ROOT / "runs/current/submitted.json"
OUT = ROOT / "runs/current/raw_model"


def call(args):
    result = subprocess.run([PYTHON, str(TOOL), *args], cwd=ROOT, check=True,
                            capture_output=True, text=True)
    return json.loads(result.stdout)


mentions = call([
    "standard", "mentions", "--dataset", str(DATASET),
    "--term", "1a", "--term", "1b", "--term", "1c", "--term", "2a",
    "--term", "3", "--term", "4a", "--term", "4b", "--term", "4e",
    "--term", "8", "--term", "9", "--term", "10", "--term", "β-ICD",
    "--context", "180",
])
(OUT / "source_mentions.json").write_text(json.dumps(mentions, ensure_ascii=False, indent=2) + "\n")

bundle = json.loads(SUBMITTED.read_text())
diagnostics = []
seen = set()
for molecule in bundle["molecules"]:
    smiles = molecule["smiles"]["normalized_value"]
    if not smiles or smiles in seen:
        continue
    seen.add(smiles)
    result = call(["standard", "stereo", "--smiles", smiles])
    diagnostics.append({
        "molecule_id": molecule["molecule_id"],
        "smiles": smiles,
        "diagnostic": result,
    })
(OUT / "stereo_diagnostics.json").write_text(json.dumps({
    "source_verified": False,
    "diagnostics": diagnostics,
    "limitations": "Graph diagnostics only; parsing/canonicalization does not verify source identity or stereochemistry.",
}, ensure_ascii=False, indent=2) + "\n")

print(json.dumps({"mention_hits": len(mentions["matches"]), "stereo_graphs": len(diagnostics)}))
