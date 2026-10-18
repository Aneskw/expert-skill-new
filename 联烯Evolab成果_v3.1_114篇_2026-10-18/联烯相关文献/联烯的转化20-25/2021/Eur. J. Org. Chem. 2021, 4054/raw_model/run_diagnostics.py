"""Run literal source mention and CLI-equivalent stereo diagnostics for retained candidates."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
import standard_records as st

run = Path(__file__).resolve().parents[1]
dataset = json.loads((run / "dataset.json").read_text())
bundle = json.loads((run / "submitted.json").read_text())

terms = ["1a", "4a", "5a", "4i", "5i", "11", "12", "N.R."]
mentions = st.source_mentions(dataset, terms, context=180)
(run / "raw_model/source_mentions.json").write_text(json.dumps(mentions, ensure_ascii=False, indent=2) + "\n")

diagnostics = []
for molecule in bundle["molecules"]:
    smiles = molecule["smiles"]["normalized_value"]
    if smiles:
        diagnostics.append({"molecule_id": molecule["molecule_id"], "smiles": smiles,
                            "diagnostic": st.stereo_inspect(smiles)})
(run / "raw_model/stereo_diagnostics.json").write_text(json.dumps(diagnostics, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"mention_hits": len(mentions["matches"]), "stereo_graphs": len(diagnostics)}, ensure_ascii=False))
