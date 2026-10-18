#!/usr/bin/env python3
"""Enrich existing v3.1 records with visually reviewed PDF pages 7-9.

Pages 7-8 complete experiments already inventoried on pages 2-3/6, so this
script merges evidence and procedure details into those records.  Page 9 has
references/X-ray context only.  Machine records remain partial/needs_review.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
import standard_records as st

DOC = "doc_41522308b99275e3"


def E(*pages: int) -> list[str]:
    out = []
    for p in pages:
        out.extend([f"{DOC}:page{p}:text", f"{DOC}:page{p}:image"])
    return out


def F(raw=None, normalized=None, pages=(), method="original_extraction",
      status=None, derivation=None, alternatives=()):
    return st.fact(raw, normalized, E(*pages), method=method, status=status,
                   derivation=derivation, alternatives=alternatives)


def missing(pages=()):
    return st.fact(evidence=E(*pages))


def q(number, unit, raw_text=None):
    value = {"number": number, "unit": unit}
    if raw_text is not None:
        value["raw_text"] = raw_text
    return value


def add_evidence(target, *pages):
    target["evidence_ids"] = list(dict.fromkeys(target.get("evidence_ids", []) + E(*pages)))


bundle = json.loads((RUN / "submitted.json").read_text())
# Refresh the document registry from the mechanical dataset so page-image
# evidence 7-9 registered by prepare/init is present in the submission.
bundle["documents"] = json.loads((RUN / "dataset.json").read_text())["documents"]
reactions = {r["reaction_id"].split(":")[-1]: r for r in bundle["reactions"]}
molecules = {m["molecule_id"]: m for m in bundle["molecules"]}


def add_molecule(mid, label, name, smiles, pages, *, concrete=True):
    if mid in molecules:
        return molecules[mid]
    m = st.new_molecule(mid, DOC, E(*pages))
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, pages)
    m["names"]["original"] = F(name, name, pages)
    m["names"]["english"] = F(name, name, pages)
    m["cas_numbers"] = missing(pages)
    if concrete:
        diag = st.stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = F(
            {"source_name": name, "candidate_smiles": smiles}, smiles, pages,
            method="structure_derivation",
            derivation="Candidate connectivity derived from the unambiguous source name; parsing is diagnostic, not source verification.")
        m["stereochemistry"]["atom_centers"] = F(
            "No stereogenic atom center in the named structure.", diag["atom_centers"], pages,
            method="structure_derivation", derivation="Standard stereo graph diagnostic.")
        m["stereochemistry"]["double_bonds"] = F(
            "No stereogenic E/Z alkene in the named structure.", diag["double_bonds"], pages,
            method="structure_derivation", derivation="Standard stereo graph diagnostic.")
        m["stereochemistry"]["axial"] = F(
            "No stereogenic axis in the named structure.", [], pages,
            method="structure_derivation", derivation="Structure contains no allene or stereogenic biaryl axis.")
        m["stereochemistry"]["mixture"] = missing(pages)
        m["stereochemistry"]["assignment"] = missing(pages)
        m["stereochemistry"]["chirality_status"] = F(
            None, "achiral", pages, method="inference",
            derivation="Assigned from the complete candidate graph, not from absence of a source descriptor.")
        m["stereochemistry"]["geometry_status"] = F(
            None, "not_applicable", pages, method="inference",
            derivation="No E/Z-applicable alkene is present in the candidate graph.")
        m["stereochemistry"]["representation_status"] = F(
            None, "not_applicable", pages, method="inference",
            derivation="No stereochemical encoding is applicable to the candidate graph.")
        m["stereochemistry"]["isomeric_smiles"] = F(
            {"candidate_smiles": smiles}, diag["canonical_isomeric_smiles"], pages,
            method="structure_derivation", derivation="Standard stereo canonicalization; no stereochemical elements are present.")
        m["stereochemistry"]["identity_key"] = F(
            {"candidate_smiles": smiles}, diag["stereo_graph_key"], pages,
            method="structure_derivation", derivation="Standard stereo graph key for comparison only.")
    else:
        m["representation_type"] = "unresolved"
        m["stereochemistry"]["chirality_status"] = F(
            label, "not_applicable", pages, method="inference",
            derivation="Bulk inorganic stationary phase is not treated as a discrete stereochemical molecule.")
        m["stereochemistry"]["geometry_status"] = F(
            label, "not_applicable", pages, method="inference",
            derivation="Bulk inorganic stationary phase has no molecular E/Z assignment.")
        m["stereochemistry"]["representation_status"] = F(
            label, "unresolved", pages, method="inference",
            derivation="The source gives a material name, not a unique molecular structure.")
    m["identity_note"] = F(
        f"{name} is source-named; CAS not reported.",
        "Source identity retained without external database enrichment.", pages)
    molecules[mid] = m
    bundle["molecules"].append(m)
    return m


# Explicit page-7 workup/purification materials that were absent from the prior registry.
add_molecule("hydrochloric_acid", "0.1 N HCl", "hydrochloric acid", "Cl", (7,))
add_molecule("sodium_bicarbonate", "NaHCO3", "sodium bicarbonate", "O=C([O-])O.[Na+]", (7,))
add_molecule("hexane", "hexane", "hexane", "CCCCCC", (7, 8))
add_molecule("neutral_silica_gel", "neutral silica gel", "neutral silica gel", None, (7, 8), concrete=False)


def participant(r, pid, mid, role, pages, stage="reaction"):
    found = next((p for p in r["participants"] if p["participant_id"] == pid), None)
    if found is None:
        found = st.new_participant(pid, mid, role, E(*pages), stage)
        r["participants"].append(found)
    return found


def charge(p, key, raw, normalized, pages, method="original_extraction", alternatives=()):
    p["charges"][key] = F(raw, normalized, pages, method=method, alternatives=alternatives,
                           status="conflict" if alternatives else None)


def step(sid, pages, stage, operation, *, temperature=None, time=None,
         atmosphere=None, pids=()):
    s = st.new_step(sid, E(*pages), stage)
    s["participant_ids"] = list(pids)
    s["operation"] = F(operation, operation, pages)
    if temperature is not None:
        s["temperature"] = F(temperature, temperature, pages)
    if time is not None:
        s["time"] = F(time.get("raw_text", time) if isinstance(time, dict) else time, time, pages)
    if atmosphere is not None:
        s["atmosphere"] = F(atmosphere, atmosphere, pages)
    return s


def outcome_basis(r, pages, product_mass, product_mid, yield_value):
    p = next(p for p in r["participants"] if p["molecule_id"] == product_mid and p["role"] == "product")
    charge(p, "mass", f"{product_mass} mg", q(product_mass, "mg", f"{product_mass} mg"), pages)
    o = next(o for o in r["outcomes"] if product_mid in o["target_molecule_ids"] and o["kind"] == "yield")
    o["yield_type"] = F(f"{yield_value}% after general procedure/column", "isolated yield", pages)
    o["measurement_basis"] = F(
        f"to give compound {product_mid} ({product_mass} mg, {yield_value}%)",
        "isolated after neutral silica gel column chromatography; product mass reported", pages)


def merge_occurrence(mid, *pages):
    m = molecules[mid]
    occ = next(o for o in m["occurrences"] if o["document_id"] == DOC)
    occ["evidence_ids"] = list(dict.fromkeys(occ["evidence_ids"] + E(*pages)))


# Page 7 completes the Procedure-B preparations 2m, 2n and 2o.
for suffix, substrate, product, input_mass, product_mass, yld in (
    ("table2_2m", "1m", "2m", 55, 40, 70),
    ("table2_2n", "1n", "2n", 63, 30, 45),
    ("table2_2o", "1o", "2o", 58, 46, 76),
):
    r = reactions[suffix]
    add_evidence(r, 7)
    r["location"] = F(
        f"Table 2 {product}; general procedure B; specific preparation on PDF page 7",
        {"pages": sorted({int(x.split(":page")[1].split(":")[0]) for x in r["evidence_ids"]}),
         "location": f"Table 2 {product}; general procedure B; specific preparation"},
        (2, 3, 6, 7))
    sub = next(p for p in r["participants"] if p["participant_id"] == "substrate")
    alternatives = [
        F("50 mg (general procedure B)", q(50, "mg", "50 mg"), (6,), method="inherited_general_condition"),
        F(f"{input_mass} mg (specific preparation)", q(input_mass, "mg", f"{input_mass} mg"), (7,)),
    ]
    sub["charges"]["mass"] = F(
        [a["raw_value"] for a in alternatives], None, (6, 7), status="conflict", alternatives=alternatives)
    charge(sub, "amount", "0.16 mmol", q(0.16, "mmol", "0.16 mmol"), (7,))
    r["reaction_scale"] = F(
        f"{substrate} (0.16 mmol)",
        {"reference_molecule_id": substrate, "amount": q(0.16, "mmol")}, (7,))
    outcome_basis(r, (7,), product_mass, product, yld)
    r["procedure"]["general_procedure_reference"] = F(
        "Following general procedure B", f"{DOC}:general_procedure_B", (7,))
    r["procedure"]["original_text"] = F(
        f"Following general procedure B, compound {substrate} ({input_mass} mg, 0.16 mmol) was converted to {product} ({product_mass} mg, {yld}%).",
        "Specific preparation follows general procedure B; characterization follows the synthesis statement.", (7,))
    r["notes"] = [
        n for n in r["notes"]
        if "General/specific substrate-mass discrepancy" not in str(n.get("normalized_value"))
        and "specific preparations" not in str(n.get("raw_value"))
    ]
    r["notes"].append(F(
        f"General procedure B states 50 mg, whereas the specific {product} preparation states {input_mass} mg.",
        "General/specific substrate-mass discrepancy retained as a conflict; neither value selected.", (6, 7)))
    r["review"]["unresolved"] = [
        f"General-procedure 50 mg versus specific {input_mass} mg input is retained as a conflict.",
        "Complete stereochemical assignment of the potentially stereogenic cyclobutene-ring carbon is not reported."
    ]
    merge_occurrence(substrate, 7)
    merge_occurrence(product, 7)


def add_common_workup(r, pages, *, hcl=False, bicarbonate=False, chromatography="20:1"):
    pids = []
    if hcl:
        pids.append(participant(r, "hcl_quench", "hydrochloric_acid", "workup_agent", pages, "workup")["participant_id"])
    else:
        pids.append(participant(r, "water_quench", "water", "workup_agent", pages, "workup")["participant_id"])
    pids.extend([
        participant(r, "acoet_extract", "ethyl_acetate", "workup_agent", pages, "workup")["participant_id"],
        participant(r, "brine_wash", "brine", "workup_agent", pages, "workup")["participant_id"],
        participant(r, "na2so4_dry", "sodium_sulfate", "workup_agent", pages, "workup")["participant_id"],
    ])
    if bicarbonate:
        pids.insert(2, participant(r, "nahco3_wash", "sodium_bicarbonate", "workup_agent", pages, "workup")["participant_id"])
    pur = [
        participant(r, "hexane_eluent", "hexane", "purification_agent", pages, "purification")["participant_id"],
        participant(r, "acoet_eluent", "ethyl_acetate", "purification_agent", pages, "purification")["participant_id"],
        participant(r, "neutral_silica", "neutral_silica_gel", "purification_agent", pages, "purification")["participant_id"],
    ]
    return pids, pur, chromatography


# Full desilylation procedure for 2f.
r = reactions["scheme3_desilylation_2f"]
add_evidence(r, 7)
sub = next(p for p in r["participants"] if p["participant_id"] == "substrate")
charge(sub, "mass", "150 mg", q(150, "mg", "150 mg"), (7,))
charge(sub, "amount", "0.45 mmol", q(0.45, "mmol", "0.45 mmol"), (7,))
r["reaction_scale"] = F("2a (0.45 mmol)", {"reference_molecule_id": "2a", "amount": q(0.45, "mmol")}, (7,))
r["solvent_conditions"]["composition"] = F("DMF", "DMF", (7,))
r["solvent_conditions"]["concentration"] = F("0.06 M", q(0.06, "M", "0.06 M"), (7,))
r["solvent_conditions"]["concentration_basis"] = F("compound 2a", "2a", (7,))
work, pur, ratio = add_common_workup(r, (7,), hcl=True, bicarbonate=True, chromatography="3:1")
next(p for p in r["participants"] if p["participant_id"] == "hcl_quench")["charges"]["volume"] = missing((7,))
r["steps"] = [
    step("add_tbaf", (7,), "reaction", "add tetrabutylammonium fluoride trihydrate to 2a in DMF", pids=("substrate", "solvent", "tba")),
    step("heat_hold", (7,), "reaction", "stir", temperature=q(80, "°C", "80 °C"), time=q(4, "h", "4 h"), atmosphere="nitrogen", pids=("substrate", "solvent", "tba")),
    step("quench_extract_wash_dry", (7,), "workup", "quench with 0.1 N HCl; extract with AcOEt; wash with NaHCO3 and brine; dry over Na2SO4; concentrate in vacuo", pids=work),
    step("column", (7,), "purification", "neutral silica gel column chromatography; hexane/AcOEt = 3:1", pids=pur),
]
outcome_basis(r, (7,), 74, "2f", 94)
r["procedure"]["original_text"] = F(
    "2a (150 mg, 0.45 mmol) in DMF (0.06 M); n-Bu4NF·3H2O (3.0 equiv.); N2, 80 °C, 4 h; 0.1 N HCl quench; AcOEt extraction; NaHCO3/brine wash; Na2SO4; neutral silica, hexane/AcOEt 3:1; 2f (74 mg, 94%).",
    "Full page-7 desilylation preparation and isolation sequence.", (7,))
r["procedure"]["addition_order"] = F("TBAF trihydrate was added to 2a in DMF.", ["2a in DMF", "add TBAF trihydrate"], (7,))
r["review"]["unresolved"] = ["The volume of 0.1 N HCl and all extraction/wash volumes are unreported."]
merge_occurrence("2f", 7)


# Full Wittig procedure for compound 4; preserve qualitative overnight duration.
r = reactions["scheme3_wittig_4"]
add_evidence(r, 7)
sub = next(p for p in r["participants"] if p["participant_id"] == "substrate")
charge(sub, "mass", "100 mg", q(100, "mg", "100 mg"), (7,))
charge(sub, "amount", "0.3 mmol", q(0.3, "mmol", "0.3 mmol"), (7,))
r["reaction_scale"] = F("2a (0.3 mmol)", {"reference_molecule_id": "2a", "amount": q(0.3, "mmol")}, (7,))
r["solvent_conditions"]["composition"] = F("THF", "THF", (7,))
r["solvent_conditions"]["concentration"] = F("0.2 M", q(0.2, "M", "0.2 M"), (7,))
r["solvent_conditions"]["concentration_basis"] = F("methyltriphenylphosphonium bromide", "methyltriphenylphosphonium bromide", (7,))
work, pur, ratio = add_common_workup(r, (7,), chromatography="20:1")
r["steps"] = [
    step("add_base_at_zero", (7,), "reaction", "add t-BuOK to methyltriphenylphosphonium bromide in THF", temperature=q(0, "°C", "0 °C"), pids=("ylide", "base", "solvent")),
    step("ylide_hold", (7,), "reaction", "stir", temperature="room temperature", time=q(30, "min", "30 minutes"), pids=("ylide", "base", "solvent")),
    step("cool_again", (7,), "reaction", "cool to 0 °C again", temperature=q(0, "°C", "0 °C"), pids=("ylide", "base", "solvent")),
    step("add_2a", (7,), "reaction", "add compound 2a", temperature=q(0, "°C", "0 °C"), pids=("substrate", "ylide", "base", "solvent")),
    step("overnight_hold", (7,), "reaction", "stir", time="overnight", atmosphere="nitrogen", pids=("substrate", "ylide", "base", "solvent")),
    step("quench_extract_wash_dry", (7,), "workup", "quench with water; extract with AcOEt; wash with brine; dry over Na2SO4; concentrate in vacuo", pids=work),
    step("column", (7,), "purification", "neutral silica gel column chromatography; hexane/AcOEt = 20:1", pids=pur),
]
outcome_basis(r, (7,), 75, "4", 76)
r["procedure"]["original_text"] = F(
    "Ph3PMeBr (2.2 equiv.) in THF (0.2 M) at 0 °C; add t-BuOK (2.2 equiv.); room temperature 30 min; cool to 0 °C; add 2a (100 mg, 0.3 mmol); N2 overnight; aqueous workup; neutral silica, hexane/AcOEt 20:1; 4 (75 mg, 76%).",
    "Full page-7 Wittig preparation and isolation sequence.", (7,))
r["procedure"]["addition_order"] = F(
    "t-BuOK was added to phosphonium salt at 0 °C; after 30 min at room temperature the mixture was recooled to 0 °C before adding 2a.",
    ["phosphonium salt in THF", "add t-BuOK at 0 °C", "30 min at room temperature", "cool to 0 °C", "add 2a", "stir overnight under N2"], (7,))
r["review"]["unresolved"] = ["The temperature during the overnight hold and the exact overnight duration are unreported.", "Workup volumes are unreported."]
merge_occurrence("4", 7)


# Full Pinnick oxidation/cyclization procedure for compound 5.
r = reactions["scheme3_pinnick_5"]
add_evidence(r, 7)
sub = next(p for p in r["participants"] if p["participant_id"] == "substrate")
charge(sub, "mass", "100 mg", q(100, "mg", "100 mg"), (7,))
charge(sub, "amount", "0.3 mmol", q(0.3, "mmol", "0.3 mmol"), (7,))
water = next(p for p in r["participants"] if p["participant_id"] == "water")
charge(water, "volume", "1.0 ml of H2O", q(1.0, "mL", "1.0 ml of H2O"), (7,))
r["reaction_scale"] = F("2a (0.3 mmol)", {"reference_molecule_id": "2a", "amount": q(0.3, "mmol")}, (7,))
r["solvent_conditions"]["composition"] = F("tert-butanol with 1.0 mL H2O", ["tert-butanol", "water"], (7,))
r["solvent_conditions"]["concentration"] = F("0.05 M", q(0.05, "M", "0.05 M"), (7,))
r["solvent_conditions"]["concentration_basis"] = F("compound 2a in tert-butanol", "2a", (7,))
work, pur, ratio = add_common_workup(r, (7,), chromatography="20:1")
r["steps"] = [
    step("add_reagents", (7,), "reaction", "add 2-methyl-2-butene and an aqueous NaClO2/NaH2PO4 solution to 2a in tert-butanol", temperature="ambient temperature", pids=("substrate", "solvent", "scavenger", "chlorite", "phosphate", "water")),
    step("ambient_hold", (7,), "reaction", "stir", temperature="ambient temperature", time=q(1, "h", "1 h"), atmosphere="nitrogen", pids=("substrate", "solvent", "scavenger", "chlorite", "phosphate", "water")),
    step("quench_extract_wash_dry", (7,), "workup", "quench with water; extract with AcOEt; wash with brine; dry over Na2SO4; concentrate in vacuo", pids=work),
    step("column", (7,), "purification", "neutral silica gel column chromatography; hexane/AcOEt = 20:1", pids=pur),
]
outcome_basis(r, (7,), 73, "5", 74)
r["procedure"]["original_text"] = F(
    "2a (100 mg, 0.3 mmol) in tert-butanol (0.05 M); add 2-methyl-2-butene (30 equiv.) and NaClO2/NaH2PO4 (5 equiv. each) in 1.0 ml H2O; N2, ambient temperature, 1 h; aqueous workup; neutral silica, hexane/AcOEt 20:1; 5 (73 mg, 74%).",
    "Full page-7 Pinnick oxidation/cyclization preparation and isolation sequence.", (7,))
r["procedure"]["addition_order"] = F(
    "2-methyl-2-butene and the aqueous NaClO2/NaH2PO4 solution were added to 2a in tert-butanol.",
    ["2a in tert-butanol", "add 2-methyl-2-butene", "add aqueous NaClO2/NaH2PO4 solution"], (7,))
r["review"]["unresolved"] = ["The volume of tert-butanol and all workup volumes are unreported."]
merge_occurrence("5", 7)


def add_scaleup_workup(r, pages, cond):
    if cond == "A":
        pur = [
            participant(r, "hexane_eluent", "hexane", "purification_agent", pages, "purification")["participant_id"],
            participant(r, "acoet_eluent", "ethyl_acetate", "purification_agent", pages, "purification")["participant_id"],
            participant(r, "neutral_silica", "neutral_silica_gel", "purification_agent", pages, "purification")["participant_id"],
        ]
        return [], pur
    work = [
        participant(r, "acoet_extract", "ethyl_acetate", "workup_agent", pages, "workup")["participant_id"],
        participant(r, "brine_wash", "brine", "workup_agent", pages, "workup")["participant_id"],
        participant(r, "na2so4_dry", "sodium_sulfate", "workup_agent", pages, "workup")["participant_id"],
    ]
    pur = [
        participant(r, "hexane_eluent", "hexane", "purification_agent", pages, "purification")["participant_id"],
        participant(r, "acoet_eluent", "ethyl_acetate", "purification_agent", pages, "purification")["participant_id"],
        participant(r, "neutral_silica", "neutral_silica_gel", "purification_agent", pages, "purification")["participant_id"],
    ]
    return work, pur


# Scheme 4 condition A, fully described on page 7.
r = reactions["scheme4_scaleup_a"]
add_evidence(r, 7)
sub = next(p for p in r["participants"] if p["participant_id"] == "substrate")
charge(sub, "mass", "625 mg", q(625, "mg", "625 mg"), (7,))
charge(sub, "amount", "2 mmol", q(2, "mmol", "2 mmol"), (7,))
r["reaction_scale"] = F("1a (2 mmol)", {"reference_molecule_id": "1a", "amount": q(2, "mmol")}, (7,))
r["solvent_conditions"]["composition"] = F("dry toluene with H2O (2 equiv.)", ["dry toluene", "H2O"], (7,))
r["solvent_conditions"]["concentration"] = F("0.03 M", q(0.03, "M", "0.03 M"), (7,))
r["solvent_conditions"]["concentration_basis"] = F("compound 1a", "1a", (7,))
work, pur = add_scaleup_workup(r, (7,), "A")
r["steps"] = [
    step("add_water", (7,), "reaction", "add H2O to 1a in dry toluene", pids=("substrate", "solvent", "water")),
    step("sonicate", (7,), "reaction", "sonicate", time=q(10, "min", "10 minutes"), pids=("substrate", "solvent", "water")),
    step("reflux", (7,), "reaction", "reflux", temperature="reflux", time=q(2, "h", "2 hours"), atmosphere="nitrogen", pids=("substrate", "solvent", "water")),
    step("cool_concentrate", (7,), "workup", "cool to room temperature and concentrate in vacuo", temperature="room temperature"),
    step("column", (7,), "purification", "neutral silica gel column chromatography; hexane/AcOEt = 5:1", pids=pur),
]
outcome_basis(r, (7,), 549, "2a", 83)
r["procedure"]["original_text"] = F(
    "Condition A: 1a (625 mg, 2 mmol) in dry toluene (0.03 M), H2O (2 equiv.); sonicate 10 min; reflux under N2 2 h; cool, concentrate; neutral silica, hexane/AcOEt 5:1; 2a (549 mg, 83%).",
    "Full page-7 scale-up condition A.", (7,))
r["procedure"]["pretreatment"] = F("the mixture was sonicated for 10 minutes", q(10, "min", "10 minutes"), (7,))
r["review"]["unresolved"] = ["The source does not report solvent or workup volumes beyond the 0.03 M concentration."]


# Scheme 4 condition B spans pages 7-8.
r = reactions["scheme4_scaleup_b"]
add_evidence(r, 7, 8)
sub = next(p for p in r["participants"] if p["participant_id"] == "substrate")
charge(sub, "mass", "625 mg", q(625, "mg", "625 mg"), (7,))
charge(sub, "amount", "2 mmol", q(2, "mmol", "2 mmol"), (7,))
r["reaction_scale"] = F("1a (2 mmol)", {"reference_molecule_id": "1a", "amount": q(2, "mmol")}, (7,))
r["solvent_conditions"]["composition"] = F("H2O", "water", (7,))
r["solvent_conditions"]["concentration"] = F("0.03 M", q(0.03, "M", "0.03 M"), (7,))
r["solvent_conditions"]["concentration_basis"] = F("compound 1a", "1a", (7,))
work, pur = add_scaleup_workup(r, (7, 8), "B")
r["steps"] = [
    step("sonicate", (7,), "reaction", "sonicate the solution of 1a in water", time=q(10, "min", "10 minutes"), pids=("substrate", "solvent")),
    step("reflux", (7, 8), "reaction", "reflux", temperature="reflux", time=q(2, "h", "2 hours"), atmosphere="air", pids=("substrate", "solvent")),
    step("cool_extract_wash_dry", (8,), "workup", "cool to room temperature; extract with AcOEt; wash with brine; dry over Na2SO4; concentrate in vacuo", temperature="room temperature", pids=work),
    step("column", (8,), "purification", "neutral silica gel column chromatography; hexane/AcOEt = 5:1", pids=pur),
]
outcome_basis(r, (8,), 450, "2a", 68)
r["procedure"]["original_text"] = F(
    "Condition B: 1a (625 mg, 2 mmol) in H2O (0.03 M); sonicate 10 min; reflux under air 2 h; cool; AcOEt extraction, brine wash, Na2SO4, concentrate; neutral silica, hexane/AcOEt 5:1; 2a (450 mg, 68%).",
    "Scale-up condition B begins on page 7 and continues on page 8.", (7, 8))
r["procedure"]["pretreatment"] = F("the solution was sonicated for 10 minutes", q(10, "min", "10 minutes"), (7,))
r["review"]["unresolved"] = ["Extraction and wash volumes are unreported."]


# Complete page-8 D2O-labeling procedure; analytical basis of the 1:1 ratio remains missing.
r = reactions["scheme5_d2o_labeling"]
add_evidence(r, 8)
sub = next(p for p in r["participants"] if p["participant_id"] == "substrate")
charge(sub, "mass", "50 mg", q(50, "mg", "50 mg"), (8,))
charge(sub, "amount", "0.16 mmol", q(0.16, "mmol", "0.16 mmol"), (8,))
r["reaction_scale"] = F("1a (0.16 mmol)", {"reference_molecule_id": "1a", "amount": q(0.16, "mmol")}, (8,))
r["solvent_conditions"]["composition"] = F("dioxane", "dioxane", (8,))
r["solvent_conditions"]["concentration"] = F("0.03 M", q(0.03, "M", "0.03 M"), (8,))
r["solvent_conditions"]["concentration_basis"] = F("compound 1a", "1a", (8,))
pur = [
    participant(r, "hexane_eluent", "hexane", "purification_agent", (8,), "purification")["participant_id"],
    participant(r, "acoet_eluent", "ethyl_acetate", "purification_agent", (8,), "purification")["participant_id"],
    participant(r, "neutral_silica", "neutral_silica_gel", "purification_agent", (8,), "purification")["participant_id"],
]
r["steps"] = [
    step("flame_dry_purge_seal", (8,), "reaction", "place 1a/dioxane solution in a flame-dried test tube, purge with nitrogen and seal", atmosphere="nitrogen", pids=("substrate", "solvent")),
    step("add_d2o", (8,), "reaction", "add D2O using a microsyringe", pids=("substrate", "solvent", "d2o")),
    step("sonicate", (8,), "reaction", "sonicate", time=q(10, "min", "10 minutes"), pids=("substrate", "solvent", "d2o")),
    step("reflux", (8,), "reaction", "reflux", temperature="reflux", time=q(2, "h", "2 hours"), atmosphere="nitrogen (sealed after purge)", pids=("substrate", "solvent", "d2o")),
    step("cool_concentrate", (8,), "workup", "cool to room temperature and concentrate in vacuo", temperature="room temperature"),
    step("column", (8,), "purification", "neutral silica gel column chromatography; hexane/AcOEt = 5:1", pids=pur),
]
r["procedure"]["original_text"] = F(
    "1a (50 mg, 0.16 mmol) in dioxane (0.03 M), flame-dried tube, N2 purge/seal; add D2O (2 equiv.) by microsyringe; sonicate 10 min; reflux 2 h; cool, concentrate; neutral silica, hexane/AcOEt 5:1; 1:1 deuterium-labeled 2a:unlabeled 2a.",
    "Full page-8 deuterium-labeling procedure and reported mixture ratio.", (8,))
r["procedure"]["addition_order"] = F("D2O was added by microsyringe to the sealed 1a/dioxane solution.", ["1a in dioxane", "purge and seal", "add D2O", "sonicate", "reflux"], (8,))
r["procedure"]["pretreatment"] = F("flame-dried test tube; sonicated for 10 minutes", ["flame-dried test tube", q(10, "min", "10 minutes")], (8,))
o = next(o for o in r["outcomes"] if o["outcome_id"] == "labeled_unlabeled_ratio")
o["value"] = F("1:1 mixture of deuterium-labeled 2a and unlabeled 2a", "1:1", (8,))
o["scope"] = F("deuterium-labeled 2a : unlabeled 2a", {"metric": "isotopologue ratio", "order": ["2a_d_labeled", "2a"]}, (8,))
o["measurement_basis"] = missing((8,))
r["review"]["unresolved"] = [
    "The analytical method and uncertainty for the 1:1 isotopologue ratio are not reported.",
    "No total yield or isolated mass of the mixture is reported.",
    "The drawn D/OD positions remain textually preserved because a verified isotope/stereo SMILES was not established."
]
merge_occurrence("1a", 8)
merge_occurrence("2a", 8)
merge_occurrence("2a_d_labeled", 8)


# Page-7 characterization and page-9 X-ray note are source context, not reactions.
ctx = {x["context_id"]: x for x in bundle["source_context"]}
# Keep the earlier Figure-1 context byte-for-byte stable because standard merge
# preserves old source context.  The page-9 note is a distinct occurrence.
ctx[f"{DOC}:figure1_xray"]["evidence_ids"] = [f"{DOC}:page2:text"]
ctx[f"{DOC}:figure1_xray"]["description"] = (
    "Figure 1 X-ray structure of 2a; text states the structure was determined by X-ray crystallography (CCDC 1958998)."
)
page7_context = {
    "context_id": f"{DOC}:characterization_page7",
    "document_id": DOC,
    "kind": "characterization",
    "description": "Page 7 contains 1H/13C NMR, HRMS and melting-point data for 2m, 2n, 2o, 2f, 4, 5 and 3a; peak lists are not separate synthesis experiments.",
    "evidence_ids": E(7),
}
if page7_context["context_id"] not in ctx:
    bundle["source_context"].append(page7_context)
reference11_context = {
    "context_id": f"{DOC}:reference11_xray",
    "document_id": DOC,
    "kind": "structure_characterization",
    "description": "Reference note 11 repeats that 2a was determined by X-ray structure analysis and gives CCDC 1958998; this is not a separate synthesis.",
    "evidence_ids": E(9),
}
if reference11_context["context_id"] not in ctx:
    bundle["source_context"].append(reference11_context)


# Raw page inventory and exact requested coverage artifact.
inventory = {
    "document_id": DOC,
    "requested_pages": [7, 8, 9],
    "visual_review": "All three registered page images were opened at original detail.",
    "new_reaction_records": 0,
    "enriched_existing_reaction_records": 9,
    "new_molecule_records": 4,
    "decisions": [
        "Pages 7-8 complete existing cross-page experiments; no duplicate reactions were created.",
        "Wittig duration is preserved as qualitative 'overnight', not converted to 12 h.",
        "Page-7 spectra and page-9 X-ray note are source_context, not synthesis records.",
        "Page 9 contains no new reaction arrow, table row or execution.",
    ],
}
st.write(RUN / "raw_model/pages_0007_0009_extraction.json", inventory)

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 7,
            "state": "attempted",
            "regions": [
                {"region": "upper left, continuation of 2m followed by preparations 2n and 2o", "finding": "The 2m cross-page sentence is completed (40 mg, 70%); full specific charges/masses for 2n and 2o are visible and merged into their existing Procedure-B records. NMR/HRMS/melting points are characterization context.", "evidence_ids": E(7)},
                {"region": "lower left, Preparation of 2f from 2a", "finding": "Complete desilylation execution, workup, column conditions, isolated mass and 94% yield; merged into the existing Scheme 3 record.", "evidence_ids": E(7)},
                {"region": "bottom left to upper right, preparation of compound 4", "finding": "Complete Wittig execution. The timed stages are 0 °C addition, room-temperature 30 min hold, recooling to 0 °C, 2a addition, then an overnight nitrogen hold of unreported temperature/duration; merged into the existing record.", "evidence_ids": E(7)},
                {"region": "middle right, preparation of compound 5", "finding": "Complete Pinnick oxidation/cyclization execution, including the 1.0 mL aqueous reagent solution, workup, isolation and 74% yield; merged into the existing record.", "evidence_ids": E(7)},
                {"region": "lower right, 3a analytical block", "finding": "Only 1H/13C NMR characterization of 3a; retained as source_context and not made into a synthesis experiment.", "evidence_ids": E(7)},
                {"region": "lower right, Scheme 4 scale-up procedure", "finding": "Condition A is complete on this page; condition B begins here and continues on page 8. They remain two distinct pre-existing scale-up records.", "evidence_ids": E(7)},
            ],
            "record_ids": [
                f"{DOC}:table2_2m", f"{DOC}:table2_2n", f"{DOC}:table2_2o",
                f"{DOC}:scheme3_desilylation_2f", f"{DOC}:scheme3_wittig_4", f"{DOC}:scheme3_pinnick_5",
                f"{DOC}:scheme4_scaleup_a", f"{DOC}:scheme4_scaleup_b",
            ],
            "unresolved": [
                "General procedure B states 50 mg, while the specific 2m/2n/2o preparations state 55/63/58 mg; both values are retained as conflicts.",
                "The Wittig overnight-hold temperature and exact duration are not reported.",
                "Workup/extraction volumes are mostly unreported.",
                "No absolute configuration is assigned for potentially stereogenic cyclobutene-ring centers.",
            ],
        },
        {
            "page": 8,
            "state": "attempted",
            "regions": [
                {"region": "upper left, continuation of Scheme 4 condition B", "finding": "Completes the water-solvent scale-up workup, column isolation, 450 mg product mass and 68% isolated yield; merged with the condition-B record begun on page 7.", "evidence_ids": E(8)},
                {"region": "upper left, Procedure for deuterium labeling using D2O", "finding": "Complete mechanistic isotope-labeling execution with scale, concentration, D2O loading, vessel preparation, sonication, reflux and column; reports a 1:1 labeled:unlabeled mixture without analytical method or total yield.", "evidence_ids": E(8)},
                {"region": "remaining left column and full right column", "finding": "Conflicts statement, acknowledgements and references 1-10 only; no current-paper reaction execution or reaction arrow.", "evidence_ids": E(8)},
            ],
            "record_ids": [f"{DOC}:scheme4_scaleup_b", f"{DOC}:scheme5_d2o_labeling"],
            "unresolved": [
                "The D2O experiment does not report the analytical method/uncertainty for the 1:1 ratio or a total isolated yield.",
                "Extraction and wash volumes for condition B are unreported.",
            ],
        },
        {
            "page": 9,
            "state": "attempted",
            "regions": [
                {"region": "top left and top right", "finding": "Continuation of reference 10 and reference note 11 only. Note 11 states the X-ray assignment of 2a (CCDC 1958998); this is characterization context, not an experimental synthesis record.", "evidence_ids": E(9)},
                {"region": "remainder of page", "finding": "Blank journal page area; no scheme, table, arrow, experimental paragraph or reaction execution.", "evidence_ids": [f"{DOC}:page9:image"]},
            ],
            "record_ids": [],
            "unresolved": ["No reaction record is associated with this page because it contains references/characterization note only."],
        },
    ],
}
st.write(RUN / "coverage/pages_0007_0009.json", coverage)

# New-molecule standard-stereo diagnostics are raw model artifacts.
stereo = {mid: st.stereo_inspect(molecules[mid]["smiles"]["normalized_value"])
          for mid in ("hydrochloric_acid", "sodium_bicarbonate", "hexane")}
st.write(RUN / "raw_model/stereo_checks_pages_0007_0009.json", stereo)

# Guard the batch invariants before writing the submission.
assert len(bundle["reactions"]) == 57
assert all(r["review"]["status"] == "needs_review" and r["completeness"] == "partial" for r in bundle["reactions"])
assert all(not r["classification"]["model_votes"] for r in bundle["reactions"])
st.write(RUN / "submitted.json", bundle)
print(json.dumps({"reactions": len(bundle["reactions"]), "molecules": len(bundle["molecules"]), "source_context": len(bundle["source_context"]), "enriched": 9}, ensure_ascii=False))
