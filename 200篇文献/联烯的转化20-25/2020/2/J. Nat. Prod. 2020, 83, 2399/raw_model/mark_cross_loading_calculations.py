#!/usr/bin/env python3
"""Mark cross-cyclomagnesiation ratio calculations as inferred in submitted v3 facts."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
import standard_records as st

path = ROOT / "runs/current/submitted.json"
bundle = json.loads(path.read_text())
derivations = {
    "0.1 relative to 10 of 15a/15b = 1 mol% (page 3 ratio)": "Calculated as 0.1/10 × 100 mol% using the printed ratio and 15a/15b as denominator.",
    "0.4 mmol relative to 4.0 mmol 15a/15b = 10 mol% (Experimental Section)": "Calculated as 0.4/4.0 × 100 mol% from the printed charges.",
    "36:10 = 3.6 equiv (page 3 ratio)": "Calculated as 36/10 using 15a/15b as the one-equivalent denominator.",
    "40.0 mmol relative to 4.0 mmol 15a/15b = 10 equiv (Experimental Section)": "Calculated as 40.0/4.0 using 15a/15b as the one-equivalent denominator.",
    "24:10 = 2.4 equiv (page 3 ratio)": "Calculated as 24/10 using 15a/15b as the one-equivalent denominator.",
    "30.0 mmol relative to 4.0 mmol 15a/15b = 7.5 equiv (Experimental Section)": "Calculated as 30.0/4.0 using 15a/15b as the one-equivalent denominator.",
}
changed = 0
for reaction in bundle["reactions"]:
    if reaction["reaction_id"] not in {
        "doc_2992316b1e64fc22:scheme4_10_15a_to_17",
        "doc_2992316b1e64fc22:scheme4_10_15b_to_18",
    }:
        continue
    for participant in reaction["participants"]:
        for charge in participant["charges"].values():
            for alternative in charge.get("alternatives", []):
                raw = alternative.get("raw_value")
                if raw in derivations:
                    alternative["acquisition_method"] = "inference"
                    alternative["verification"] = "inferred"
                    alternative["derivation"] = derivations[raw]
                    changed += 1
assert changed == 12, changed
st.write(path, bundle)
print({"updated_calculated_alternatives": changed})
