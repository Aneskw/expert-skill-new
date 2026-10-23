#!/usr/bin/env python3
"""Deterministic acceptance checks for the PDF-page 7-11 extraction."""
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
import reaction_tools as rt
import standard_records as st

RUN = ROOT / "runs/current"
data = rt.read(RUN / "dataset.json")
standard = data["standard"]
coverage = rt.read(RUN / "coverage/pages_0007_0011.json")
inventory = rt.read(RUN / "raw_model/pages_0007_0011_extraction.json")
updated = set(inventory["updated_reaction_ids"])
strict_ids = {r["record_id"] for r in data["experiments"]}
errors = rt.validate(data, RUN, strict_record_ids=strict_ids) + st.validate(standard, RUN)
by_id = {r["reaction_id"]: r for r in standard["reactions"]}

def value(f):
    return f["normalized_value"]

conflict35 = by_id["doc_c724145241bb10a8:scheme10_14_to_35"]["outcomes"][0]["value"]
hwe_temp = next(s for s in by_id["doc_c724145241bb10a8:scheme10_36_to_37"]["steps"] if s["step_id"] == "hwe_hold")["temperature"]
side = next(o for o in by_id["doc_c724145241bb10a8:scheme8_28_to_29"]["outcomes"] if o["outcome_id"] == "phosphate_adduct_yield")
md_22 = (RUN / "records/doc_c724145241bb10a8__scheme6_3a_to_22.md").read_text()
md_24 = (RUN / "records/doc_c724145241bb10a8__scheme6_23_to_24.md").read_text()
md_conflict = (RUN / "records/doc_c724145241bb10a8__scheme10_36_to_37.md").read_text()

checks = {
    "validation_passed": not errors,
    "coverage_pages_exact": [p["page"] for p in coverage["pages"]] == [7, 8, 9, 10, 11],
    "coverage_regions_nonempty": all(p["regions"] for p in coverage["pages"]),
    "coverage_record_ids_exist": all(rid in by_id for p in coverage["pages"] for rid in p["record_ids"]),
    "updated_record_count": len(updated) == 15,
    "all_updated_records_need_review": all(by_id[rid]["review"]["status"] == "needs_review" for rid in updated),
    "all_updated_records_partial": all(by_id[rid]["completeness"] == "partial" for rid in updated),
    "compound35_conflict_preserved": conflict35["data_status"] == "conflict" and value(conflict35) is None and len(conflict35["alternatives"]) == 2,
    "hwe_temperature_conflict_preserved": hwe_temp["data_status"] == "conflict" and value(hwe_temp) is None and len(hwe_temp["alternatives"]) == 2,
    "footnote39_targeted_yield": value(side["value"])["number"] == 64 and side["target_molecule_ids"] == ["phosphate_adduct_pair"],
    "footnote39_cc_label_present": any(x["code"] == "C-C" for x in value(by_id["doc_c724145241bb10a8:scheme8_28_to_29"]["classification"]["labels"])),
    "no_si_association_claimed": data["documents"][0]["supplementary_material"]["status"] == "unknown",
    "title_present": data["documents"][0]["title"] == "A stereotetrad-centered approach toward pironetin: Dead ends, Detour, and evolution of the synthetic strategy",
    "markdown_reagent_amounts_present": all(x in md_22 for x in ("35 mg", "0.186 mmol", "1.3 mL", "7.45 mmol")) and all(x in md_24 for x in ("562 mg", "1.78 mmol", "315 mg, 65% yield")),
    "markdown_conflict_present": "45 °C" in md_conflict and "room temperature" in md_conflict and "conflict" in md_conflict,
    "total_reaction_count": len(standard["reactions"]),
    "total_molecule_count": len(standard["molecules"]),
    "source_reviewed_records": sum(r["review"]["status"] == "source_reviewed" for r in standard["reactions"]),
}
boolean_failures = [k for k, v in checks.items() if k not in {"total_reaction_count", "total_molecule_count", "source_reviewed_records"} and v is not True]
result = {"checks": checks, "errors": sorted(set(errors)), "boolean_failures": boolean_failures,
          "chemical_accuracy_verified": False}
st.write(RUN / "raw_model/postflight_pages_0007_0011.json", result)
if errors or boolean_failures:
    raise SystemExit(json.dumps(result, ensure_ascii=False, indent=2))
print(json.dumps(result, ensure_ascii=False, indent=2))
