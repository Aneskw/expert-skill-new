#!/usr/bin/env python3
"""Final source-faithful review patches for PDF pages 7--12.

This script patches the already merged bundle without re-running the batch
constructor, so page-7--12 records are not duplicated.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
DOC = "doc_036d1c8f4ae80aff"
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
import standard_records as sr  # noqa: E402

bundle = json.loads((RUN / "standard_records.json").read_text())
molecules = {m["molecule_id"]: m for m in bundle["molecules"]}
reactions = {r["reaction_id"]: r for r in bundle["reactions"]}


def ev(page: int) -> list[str]:
    return [f"{DOC}:page{page}:image", f"{DOC}:page{page}:text"]


def fact(raw=None, normalized=None, evidence=(), *, method="original_extraction",
         status=None, verification=None, derivation=None, alternatives=()):
    if normalized is None and raw is not None:
        normalized = raw
    return sr.fact(raw=raw, normalized=normalized, evidence=evidence, method=method,
                   status=status, verification=verification, derivation=derivation,
                   alternatives=alternatives)


def inferred(normalized, evidence, derivation, raw=None):
    return fact(raw, normalized, evidence, method="inference", status="pending_review",
                verification="inferred", derivation=derivation)


def quantity(number, unit, raw):
    return {"number": number, "unit": unit, "raw_text": raw}


def find_label(label: str):
    for m in bundle["molecules"]:
        for occ in m.get("occurrences", []):
            if sr.value(occ.get("compound_label")) == label:
                return m["molecule_id"]
    return None


def ensure_agent(label: str, smiles: str | None = None):
    existing = find_label(label)
    if existing:
        return existing
    evidence = ev(11)
    mid = f"{DOC}:mol_{label.replace(' ', '_').replace('/', '_')}"
    m = sr.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = fact(label, label, evidence)
    m["names"]["original"] = fact(label, label, evidence)
    if smiles:
        diag = sr.stereo_inspect(smiles)
        canon = diag["canonical_isomeric_smiles"]
        m["representation_type"] = "smiles"
        m["smiles"] = fact(label, canon, evidence, method="structure_derivation",
                           verification="inferred",
                           derivation="Common workup-agent graph; not a source identity assertion.")
        m["stereochemistry"]["isomeric_smiles"] = fact(
            smiles, canon, evidence, method="structure_derivation", verification="inferred",
            derivation="Canonicalized with standard stereo diagnostics.")
        m["stereochemistry"]["identity_key"] = fact(
            canon, diag["stereo_graph_key"], evidence, method="structure_derivation",
            verification="inferred", derivation="Standard stereo graph key.")
        m["stereochemistry"]["chirality_status"] = inferred(
            "achiral", evidence, "The workup agent has no stereogenic element in the encoded graph.")
        m["stereochemistry"]["geometry_status"] = inferred(
            "not_applicable", evidence, "No E/Z or ring cis/trans descriptor applies.")
        m["stereochemistry"]["representation_status"] = inferred(
            "encoded", evidence, "Standard stereo diagnostics found no omitted stereo element.")
    else:
        m["identity_note"] = fact(
            "Source workup label retained; no molecular graph asserted.",
            "Source workup label retained; no molecular graph asserted.", evidence)
        m["stereochemistry"]["chirality_status"] = inferred(
            "unresolved", evidence, "No molecular graph was encoded.")
        m["stereochemistry"]["geometry_status"] = inferred(
            "unresolved", evidence, "No molecular graph was encoded.")
        m["stereochemistry"]["representation_status"] = inferred(
            "unresolved", evidence, "No molecular graph was encoded.")
    bundle["molecules"].append(m)
    molecules[mid] = m
    return mid


def add_participant(rid: str, pid: str, label: str, role: str, *, smiles=None,
                    volume=None, stage="workup"):
    r = reactions[rid]
    if any(p["participant_id"] == pid for p in r["participants"]):
        return
    p = sr.new_participant(pid, ensure_agent(label, smiles), role, ev(11), stage)
    if volume:
        p["charges"]["volume"] = fact(volume["raw_text"], volume, ev(11))
    r["participants"].insert(-1, p)


# Scheme 4 reports 16 h, while the detailed procedure reports 15 h.  Keep both
# source values rather than selecting one.
rid5 = f"{DOC}:scheme4_2e_to_5"
r5 = reactions[rid5]
hold = next(s for s in r5["steps"] if s["step_id"] == "reaction_hold")
hold["time"] = fact(
    None, None, ev(4) + ev(11), status="conflict", verification="missing",
    alternatives=[
        fact("16 h (Scheme 4)", quantity(16, "h", "16 h"), ev(4)),
        fact("15 h (detailed procedure)", quantity(15, "h", "15 h"), ev(11)),
    ],
)
time_note = "Reaction time is 16 h in Scheme 4 but 15 h in the detailed page-11 procedure; both are retained as a source conflict."
if not any(sr.value(n) == time_note for n in r5["notes"]):
    r5["notes"].append(fact(time_note, time_note, ev(4) + ev(11)))
if "Reaction time conflict: 16 h in Scheme 4 versus 15 h in the detailed procedure." not in r5["review"]["unresolved"]:
    r5["review"]["unresolved"].append(
        "Reaction time conflict: 16 h in Scheme 4 versus 15 h in the detailed procedure.")

# Workup materials named in the detailed procedures are participants even when
# the source does not give amounts.  The one explicit quench volume is retained.
rid4 = f"{DOC}:scheme4_2e_to_4"
add_participant(rid4, "Et2O_workup", "Et2O", "workup_agent", smiles="CCOCC")
add_participant(rid4, "Na2CO3_workup", "Na2CO3", "workup_agent",
                smiles="[Na+].[Na+].[O-]C([O-])=O")
add_participant(rid4, "MgSO4_drying", "MgSO4", "workup_agent",
                smiles="[Mg+2].[O-]S([O-])(=O)=O")

add_participant(rid5, "MeOH_quench", "MeOH", "workup_agent", smiles="CO",
                volume=quantity(10, "mL", "10 mL"))
add_participant(rid5, "NaHCO3_workup", "NaHCO3", "workup_agent")
add_participant(rid5, "Et2O_workup", "Et2O", "workup_agent", smiles="CCOCC")
add_participant(rid5, "water_workup", "water", "workup_agent", smiles="O")
add_participant(rid5, "Na2SO4_drying", "Na2SO4", "workup_agent")

rid6 = f"{DOC}:scheme4_4_to_6"
add_participant(rid6, "water_workup", "water", "workup_agent", smiles="O")
add_participant(rid6, "HCl_5pct_wash", "HCl", "workup_agent")
add_participant(rid6, "brine_wash", "brine", "workup_agent")
add_participant(rid6, "Na2SO4_drying", "Na2SO4", "workup_agent")

# 2e and its page-11 derivatives are enantioenriched samples, not pure single
# stereoisomers.  Headings assign the major product as (R); the structures were
# intentionally not encoded because the paper drawing has not been converted to
# source-verified stereochemical SMILES.
stereo_mids = [
    f"{DOC}:mol_2e_deriv_start",
    f"{DOC}:mol_4",
    f"{DOC}:mol_5",
    f"{DOC}:mol_6",
]
for mid in stereo_mids:
    m = molecules[mid]
    m["stereochemistry"]["chirality_status"] = inferred(
        "mixture", ev(3) + ev(11),
        "The starting 2e is enantioenriched and the paper states that derivatization proceeds without erosion; this is not a pure enantiomer record.")
    m["stereochemistry"]["geometry_status"] = inferred(
        "not_applicable", ev(11),
        "The printed name/drawing has no E/Z-eligible alkene and no ring cis/trans descriptor.")
    m["stereochemistry"]["representation_status"] = inferred(
        "not_encoded", ev(11),
        "The textual (R) assignment is preserved but a source-verified stereochemical SMILES was not asserted.")
    m["stereochemistry"]["assignment"] = fact(
        "major enantiomer assigned (R) in the source heading; configuration retained through derivatization",
        "major enantiomer assigned (R) in the source heading; configuration retained through derivatization",
        ev(3) + ev(11))
    m["stereochemistry"]["mixture"] = fact(
        "enantioenriched sample; exact derivative-specific er/ee not reported",
        "enantioenriched sample; exact derivative-specific er/ee not reported",
        ev(3) + ev(11))

# Correct the structure-derived propargyl candidates: prop-2-yn-1-yl is
# HC#C-CH2-, not an ynamine HC#C-N linkage.
propargyl_smiles = {
    "S10h": "C#CCN(c1c(C(=O)OC)cc(C(F)(F)F)cc1)S(=O)(=O)c1ccc(C)cc1",
    "S10i": "C#CCN(c1c(C(=O)OC)ccc(OC)c1)S(=O)(=O)c1ccc(C)cc1",
    "S10j": "C#CCN(c1c(C(=O)OC)ccc(Br)c1)S(=O)(=O)c1ccc(C)cc1",
    "S10n": "C#CCN(c1c(C(=O)OC)cc(OC)c(OC)c1)S(=O)(=O)c1ccc(C)cc1",
    "S14g": "C#CCN(c1c(CO)ccc(Cl)c1)S(=O)(=O)c1ccc(C)cc1",
    "S14k": "C#CCN(c1c(CO)cc(Cl)cc1)S(=O)(=O)c1ccc(C)cc1",
    "S14l": "C#CCN(c1c(CO)c(F)ccc1)S(=O)(=O)c1ccc(C)cc1",
    "S14m": "C#CCN(c1c(CO)c(C)ccc1)S(=O)(=O)c1ccc(C)cc1",
}
for label, smiles in propargyl_smiles.items():
    mid = find_label(label)
    m = molecules[mid]
    evidence = ev(7 if label.startswith("S10") else (8 if label == "S14g" else 9))
    diag = sr.stereo_inspect(smiles)
    canon = diag["canonical_isomeric_smiles"]
    m["smiles"] = fact(label, canon, evidence, method="structure_derivation",
                       verification="inferred",
                       derivation="Candidate corrected from the printed prop-2-yn-1-yl systematic name; parsing is not source verification.")
    m["stereochemistry"]["isomeric_smiles"] = fact(
        smiles, canon, evidence, method="structure_derivation", verification="inferred",
        derivation="Canonical isomeric SMILES generated by standard stereo diagnostics.")
    m["stereochemistry"]["identity_key"] = fact(
        canon, diag["stereo_graph_key"], evidence, method="structure_derivation",
        verification="inferred", derivation="Standard stereo graph key; not a sample identity.")
    m["stereochemistry"]["chirality_status"] = inferred(
        "achiral", evidence, "No stereogenic center or axis is present in the source-derived propargyl graph.")
    m["stereochemistry"]["geometry_status"] = inferred(
        "not_applicable", evidence, "No E/Z-eligible alkene or ring cis/trans assignment is present.")
    m["stereochemistry"]["representation_status"] = inferred(
        "encoded", evidence, "The structure-derived candidate contains no stereochemical element requiring encoding.")

# Standard stereo reports the cumulated double bonds as a potential stereo
# feature, but every source allene in this batch is terminal (CH2=C=), so axial
# chirality and E/Z geometry are not applicable.
terminal_allene_labels = [
    "S11d", "S11e", "S11f", "S11h", "S11i", "S11j", "S11n",
    "3a", "S12d", "S12f", "S12h", "S12i", "S12j", "S12n",
    "3d", "3e", "3f", "3g", "3h", "3i", "3j", "3k", "3l", "3m", "3n",
    "S15g", "S15k", "S15l", "S15m",
]
for label in terminal_allene_labels:
    mid = find_label(label)
    if not mid:
        continue
    m = molecules[mid]
    evidence = m["smiles"].get("evidence_ids", ev(9))
    m["stereochemistry"]["chirality_status"] = inferred(
        "achiral", evidence,
        "The source name/graph is a terminal CH2=C= allene and contains no other stereogenic element.")
    m["stereochemistry"]["geometry_status"] = inferred(
        "not_applicable", evidence,
        "E/Z is not defined for the terminal allene and no ring cis/trans descriptor is present.")
    m["stereochemistry"]["representation_status"] = inferred(
        "partial", evidence,
        "Chemically no axial descriptor is required for the terminal allene, but standard stereo flags the cumulated double bond as an unspecified potential stereo feature.")
    m["stereochemistry"]["axial"] = inferred(
        {"status": "not_applicable", "reason": "terminal allene end is CH2"}, evidence,
        "The terminal allene carbon bears two identical hydrogens.")

# Record the newly recognized time conflict in the raw inventory and coverage.
raw_path = RUN / "raw_model/pages_0007_0012_extraction.json"
raw = json.loads(raw_path.read_text())
conflict = {
    "field": "scheme4_2e_to_5.reaction_hold.time",
    "values": ["16 h (Scheme 4, page 4)", "15 h (detailed procedure, page 11)"],
    "decision": "retained as conflict; no source value selected",
}
if conflict not in raw.setdefault("source_conflicts", []):
    raw["source_conflicts"].append(conflict)
raw_path.write_text(json.dumps(raw, ensure_ascii=False, indent=2) + "\n")

coverage_path = RUN / "coverage/pages_0007_0012.json"
coverage = json.loads(coverage_path.read_text())
p11 = next(p for p in coverage["pages"] if p["page"] == 11)
message = "Scheme 4 gives 16 h for 2e→5, whereas the detailed procedure gives 15 h; retained as a conflict."
if message not in p11["unresolved"]:
    p11["unresolved"].append(message)
coverage_path.write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")

# Persist raw standard-stereo diagnostics after the propargyl corrections.
diagnostics = {}
for m in bundle["molecules"]:
    smiles = sr.value(m.get("smiles"))
    if smiles:
        diagnostics[m["molecule_id"]] = sr.stereo_inspect(smiles)
(RUN / "raw_model/pages_0007_0012_stereo_diagnostics.json").write_text(
    json.dumps(diagnostics, ensure_ascii=False, indent=2) + "\n")

(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"reactions": len(bundle["reactions"]), "molecules": len(bundle["molecules"]),
                  "patched": [rid4, rid5, rid6], "stereo_molecules": stereo_mids}))
