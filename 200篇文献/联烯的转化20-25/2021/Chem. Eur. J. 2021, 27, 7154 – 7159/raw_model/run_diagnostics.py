#!/usr/bin/env python3
"""Run the required CLI stereo and mention diagnostics and retain raw output."""
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
PY = "/EvoLab/student-workspaces/xllin/OpenEvo-chem-integration/.venv/bin/python"
TOOL = ROOT / "skills/chemical-reaction-extraction/scripts/reaction_tools.py"
RUN = ROOT / "runs/current"

bundle = json.loads((RUN / "submitted.json").read_text())
checks = {}
for mol in bundle["molecules"]:
    smi = mol["smiles"]["normalized_value"]
    if not smi:
        continue
    proc = subprocess.run(
        [PY, str(TOOL), "standard", "stereo", "--smiles", smi],
        cwd=ROOT, check=True, text=True, capture_output=True,
    )
    checks[mol["molecule_id"]] = json.loads(proc.stdout)
(RUN / "raw_model/stereo_checks.json").write_text(json.dumps(checks, ensure_ascii=False, indent=2) + "\n")

mention_cmd = [PY, str(TOOL), "standard", "mentions", "--dataset", str(RUN / "dataset.json")]
for term in ["1a", "2a", "3a", "4a", "6p", "7a"]:
    mention_cmd += ["--term", term]
proc = subprocess.run(mention_cmd, cwd=ROOT, check=True, text=True, capture_output=True)
(RUN / "raw_model/source_mentions.json").write_text(proc.stdout)
print(json.dumps({"stereo_checks": len(checks), "mention_matches": len(json.loads(proc.stdout)["matches"])}, ensure_ascii=False))
