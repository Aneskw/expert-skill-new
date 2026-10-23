#!/usr/bin/env python3
"""Build the page-7--12 v3.1 evidence update after direct page-image inspection."""
from __future__ import annotations

import copy
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DOC = "doc_2878b0a078d500b2"
EV = {p: f"{DOC}:page{p}:image" for p in range(7, 13)}
TXT = {p: f"{DOC}:page{p}:text" for p in range(7, 13)}
GPB_ID = f"{DOC}:GP_B_generic"


def dump(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def fact(raw=None, normalized=None, evidence=(), method="original_extraction", derivation=None,
         status=None, verification=None, alternatives=()):
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation,
                   status=status, verification=verification, alternatives=alternatives)


base = json.loads((RUN / "standard_records.json").read_text())
template = json.loads((RUN / "raw_model/submitted_pages_0007_0012_template.json").read_text())
old_molecules = {m["molecule_id"]: m for m in base["molecules"]}
old_reactions = {r["reaction_id"]: r for r in base["reactions"]}
incoming_molecules: dict[str, dict] = {}
incoming_reactions: list[dict] = []
stereo_diagnostics = []


def add_unique(items, values):
    for value in values:
        if value not in items:
            items.append(value)


def touch_molecule(mid: str, evidence=(), name: str | None = None, note: str | None = None):
    if mid in incoming_molecules:
        m = incoming_molecules[mid]
    elif mid in old_molecules:
        m = copy.deepcopy(old_molecules[mid])
        incoming_molecules[mid] = m
    else:
        m = st.new_molecule(mid, DOC, list(evidence))
        m["identity_status"] = "source_identified"
        label = mid.rsplit(":", 1)[-1]
        m["occurrences"][0]["compound_label"] = fact(label, label, evidence)
        m["names"]["original"] = fact(name or label, name or label, evidence)
        m["names"]["english"] = fact(name or label, name or label, evidence)
        m["names"]["abbreviations"] = fact(label, [label], evidence)
        m["stereochemistry"]["chirality_status"] = fact(
            "unresolved", "unresolved", evidence,
            method="inference", derivation="The inspected mechanism drawing does not establish a sample stereochemical identity.")
        m["stereochemistry"]["geometry_status"] = fact(
            "unresolved", "unresolved", evidence,
            method="inference", derivation="The inspected mechanism drawing was not promoted to a complete geometry-resolved structure.")
        m["stereochemistry"]["representation_status"] = fact(
            "unresolved", "unresolved", evidence,
            method="inference", derivation="No complete source-verified structure representation was assigned.")
        incoming_molecules[mid] = m
    if evidence:
        occ = next((o for o in m["occurrences"] if o["document_id"] == DOC), None)
        if occ is None:
            label = mid.rsplit(":", 1)[-1]
            occ = {"document_id": DOC, "compound_label": fact(label, label, evidence),
                   "evidence_ids": list(evidence)}
            m["occurrences"].append(occ)
        else:
            add_unique(occ["evidence_ids"], evidence)
            add_unique(occ["compound_label"]["evidence_ids"], evidence)
    if name:
        m["names"]["original"] = fact(name, name, evidence)
        m["names"]["english"] = fact(name, name, evidence)
    if note:
        m["identity_note"] = fact(note, note, evidence)
    return m


def make_named_unresolved(label: str, name: str, evidence, note: str):
    return touch_molecule(f"{DOC}:mol:{label}", evidence, name, note)


def get_reaction(short: str):
    rid = short if short.startswith(DOC + ":") else f"{DOC}:{short}"
    r = copy.deepcopy(old_reactions[rid])
    return r


def add_reaction_evidence(r, evidence):
    add_unique(r["evidence_ids"], evidence)
    loc = r["location"]
    add_unique(loc["evidence_ids"], evidence)
    pages = loc.get("normalized_value", {}).get("pages")
    if isinstance(pages, list):
        for eid in evidence:
            match = re.search(r":page(\d+):", eid)
            if match and int(match.group(1)) not in pages:
                pages.append(int(match.group(1)))
        pages.sort()


def participant(r, pid):
    return next(p for p in r["participants"] if p["participant_id"] == pid)


def set_charge(p, key, raw, normalized, evidence, method="original_extraction", derivation=None):
    p["charges"][key] = fact(raw, normalized, evidence, method=method, derivation=derivation)


def set_equiv_reference(p, molecule_id, evidence):
    p["equivalence_reference"] = fact(molecule_id, molecule_id, evidence)


def new_solvent(pid, molecule_id, volume, evidence):
    p = st.new_participant(pid, molecule_id, "solvent", evidence)
    set_charge(p, "volume", f"{volume:g} mL", {"number": volume, "unit": "mL"}, evidence)
    return p


def new_step(step_id, stage, pids, operation_raw, evidence, temperature=None, time=None,
             atmosphere=None):
    s = st.new_step(step_id, evidence, stage)
    s["participant_ids"] = list(pids)
    s["operation"] = fact(operation_raw, operation_raw, evidence)
    if temperature is not None:
        s["temperature"] = fact(str(temperature), temperature, evidence)
    if time is not None:
        s["time"] = fact(str(time), time, evidence)
    if atmosphere is not None:
        s["atmosphere"] = fact(str(atmosphere), atmosphere, evidence)
    return s


def find_outcome(r, kind):
    return next((o for o in r["outcomes"] if o["kind"] == kind), None)


def set_isolated_outcome(r, product_mid, yield_value, mass_mg, evidence, purification,
                         measurement_raw=None):
    out = find_outcome(r, "yield")
    if out is None:
        out = st.new_outcome("yield", "yield", [product_mid], evidence)
        r["outcomes"].append(out)
    out["target_molecule_ids"] = [product_mid]
    out["value"] = fact(f"{yield_value}%", {"number": yield_value, "unit": "%"}, evidence)
    out["yield_type"] = fact("obtained after purification", "isolated", evidence,
                             method="inference",
                             derivation="The product paragraph says the material was purified and obtained with a reported mass and yield; isolated is a normalization of that wording.")
    raw = measurement_raw or f"purified by {purification} and obtained in {yield_value}% yield"
    out["measurement_basis"] = fact(raw, f"yield of product obtained after {purification}", evidence)
    out["product_role"] = fact("title compound", "main product", evidence)
    amount = find_outcome(r, "product_amount")
    if amount is None:
        amount = st.new_outcome("product_mass", "product_amount", [product_mid], evidence)
        r["outcomes"].append(amount)
    amount["target_molecule_ids"] = [product_mid]
    amount["value"] = fact(f"{mass_mg:g} mg", {"number": mass_mg, "unit": "mg"}, evidence)
    amount["measurement_basis"] = fact("obtained as the title compound", "reported mass of purified product", evidence)
    amount["product_role"] = fact("title compound", "main product", evidence)


def update_product_stereo(mid: str, evidence, geometry=None, mixture=None, assignment=None,
                          chirality=None, representation=None):
    m = touch_molecule(mid, evidence)
    if geometry:
        m["stereochemistry"]["double_bonds"] = fact(
            geometry, [{"position": 1, "configuration": "Z"}, {"position": 7, "configuration": "E"}], evidence)
        m["stereochemistry"]["geometry_status"] = fact("(1Z,7E)", "specified", evidence)
        m["stereochemistry"]["representation_status"] = fact(
            "geometry appears in the source name but no complete isomeric SMILES is registered",
            representation or "not_encoded", evidence)
    if mixture:
        m["stereochemistry"]["mixture"] = fact(mixture, mixture, evidence)
        m["stereochemistry"]["chirality_status"] = fact("mixture", "mixture", evidence)
        m["stereochemistry"]["representation_status"] = fact(
            "mixture is not encoded as one structure", representation or "unresolved", evidence)
    if assignment:
        m["stereochemistry"]["assignment"] = fact(assignment, assignment, evidence)
    if chirality:
        m["stereochemistry"]["chirality_status"] = fact(chirality, chirality, evidence)
    return m


def use_gpb_small(r, evidence, purification_text, product_mid, yield_value, mass_mg):
    add_reaction_evidence(r, evidence)
    if GPB_ID not in r["general_scheme_ids"]:
        r["general_scheme_ids"].append(GPB_ID)
    # Preserve the two explicitly reported 1.0 mL solvent portions instead of calling their sum
    # an author-reported total volume.
    r["participants"] = [p for p in r["participants"] if p["participant_id"] != "solvent"]
    substrate = participant(r, "substrate")
    diazo = participant(r, "diazo")
    catalyst = participant(r, "catalyst")
    ref_mid = substrate["molecule_id"]
    set_charge(substrate, "amount", "0.2 mmol", {"number": 0.2, "unit": "mmol"}, evidence)
    set_charge(substrate, "equivalents", "1.0 equiv", {"number": 1.0, "unit": "equiv"}, evidence)
    set_equiv_reference(substrate, ref_mid, evidence)
    set_charge(diazo, "amount", "0.4 mmol", {"number": 0.4, "unit": "mmol"}, evidence)
    set_charge(diazo, "equivalents", "2.0 equiv", {"number": 2.0, "unit": "equiv"}, evidence)
    set_equiv_reference(diazo, ref_mid, evidence)
    set_charge(catalyst, "mass", "18.6 mg", {"number": 18.6, "unit": "mg"}, evidence)
    set_charge(catalyst, "amount", "0.02 mmol", {"number": 0.02, "unit": "mmol"}, evidence)
    set_charge(catalyst, "mol_percent", "10 mol %", {"number": 10, "unit": "mol%"}, evidence)
    set_equiv_reference(catalyst, ref_mid, evidence)
    ch2cl2 = f"{DOC}:mol:CH2Cl2"
    r["participants"].append(new_solvent("diazo_solution_solvent", ch2cl2, 1.0, evidence))
    r["participants"].append(new_solvent("reaction_solution_solvent", ch2cl2, 1.0, evidence))
    reaction_pids = [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"]
    r["steps"] = [
        new_step("syringe_pump_addition", "reaction", reaction_pids,
                 "a solution of aryl diazoacetate in CH2Cl2 was injected into the allenyne/Cu(PPh3)3Br solution via a syringe pump over 2 h",
                 evidence, "r.t.", {"number": 2, "unit": "h"}),
        new_step("post_addition_hold", "reaction", reaction_pids,
                 "the resulting solution was continued at room temperature for 12 h",
                 evidence, "r.t.", {"number": 12, "unit": "h"}),
        new_step("concentration", "workup", [], "the mixture was concentrated under reduced pressure", evidence),
        new_step("flash_chromatography", "purification", ["product"], purification_text, evidence),
    ]
    r["reaction_scale"] = fact("1,n-allenyne (0.2 mmol, 1.0 equiv)",
                               {"reference_molecule_id": ref_mid, "amount": {"number": 0.2, "unit": "mmol"}}, evidence)
    r["solvent_conditions"]["composition"] = fact(
        "CH2Cl2 (1.0 mL) for diazo solution; CH2Cl2 (1.0 mL) for allenyne/catalyst solution",
        [{"use": "diazo solution", "volume": {"number": 1.0, "unit": "mL"}},
         {"use": "allenyne/catalyst solution", "volume": {"number": 1.0, "unit": "mL"}}], evidence)
    r["solvent_conditions"]["total_volume"] = fact()
    r["procedure"]["general_procedure_reference"] = fact("following GP-B", GPB_ID, evidence)
    r["procedure"]["original_text"] = fact(
        "A solution of aryl diazoacetate (0.4 mmol, 2.0 equiv) in CH2Cl2 (1.0 mL) was injected into a solution of 1,n-allenyne (0.2 mmol, 1.0 equiv) and Cu(PPh3)3Br (18.6 mg, 0.02 mmol, 10 mol %) in CH2Cl2 (1.0 mL) via a syringe pump over 2 h at rt; the resulting solution was continued at rt for 12 h, concentrated, and purified by flash chromatography.",
        "GP-B small-scale procedure", evidence)
    r["procedure"]["addition_order"] = fact(
        "diazo solution injected into allenyne/Cu(PPh3)3Br solution",
        "diazo solution -> allenyne/catalyst solution", evidence)
    set_isolated_outcome(r, product_mid, yield_value, mass_mg, evidence, purification_text)
    r["review"]["status"] = "needs_review"
    r["review"]["method"] = "model"
    # The individual preparation paragraphs on pages 7--10 resolve warnings that
    # were appropriate when only the scheme captions on pages 4--6 had been read.
    r["review"]["unresolved"] = [
        u for u in r["review"]["unresolved"]
        if "yield measurement basis" not in u.lower()
        and "yield basis" not in u.lower()
        and "alkene e/z assignments are not textually stated" not in u.lower()
    ]
    add_unique(r["review"]["unresolved"], [
        "Machine record remains needs_review; structure drawings and stereochemical identity require expert confirmation.",
        "The source gives purification wording and product mass, but no independent analytical yield method or uncertainty."
    ])


# Exact product names and isolated product details from the individual paragraphs.
products = {
    "3aa": (7, "Methyl 8-Phenyl-3-tosyl-3-azabicyclo[5.2.0]nona-1,6-diene-8-carboxylate", 83, 68.0, "flash column chromatography on silica gel (hexanes/EtOAc = 7:1)"),
    "3ab": (7, "Ethyl 8-Phenyl-3-tosyl-3-azabicyclo[5.2.0]nona-1,6-diene-8-carboxylate", 81, 68.6, "flash column chromatography on silica gel (hexanes/EtOAc = 7:1)"),
    "3ac": (7, "Benzyl 8-Phenyl-3-tosyl-3-azabicyclo[5.2.0]nona-1,6-diene-8-carboxylate", 70, 68.0, "flash column chromatography on silica gel (hexanes/EtOAc = 7:1)"),
    "3ad": (7, "Allyl 8-Phenyl-3-tosyl-3-azabicyclo[5.2.0]nona-1,6-diene-8-carboxylate", 65, 56.6, "flash column chromatography on silica gel (hexanes/EtOAc = 7:1)"),
    "3ae": (7, "Methyl 8-(p-Tolyl)-3-tosyl-3-azabicyclo[5.2.0]nona-1,6-diene-8-carboxylate", 87, 73.7, "flash column chromatography on silica gel (hexanes/EtOAc = 7:1)"),
    "3af": (7, "Methyl 8-(4-(tert-Butyl)phenyl)-3-tosyl-3-azabicyclo[5.2.0]nona-1,6-diene-8-carboxylate", 85, 79.2, "flash column chromatography on silica gel (hexanes/EtOAc = 7:1)"),
    "3ag": (7, "Methyl 8-([1,1'-Biphenyl]-4-yl)-3-tosyl-3-azabicyclo[5.2.0]nona-1,6-diene-8-carboxylate", 78, 75.8, "flash column chromatography on silica gel (hexanes/EtOAc = 7:1)"),
    "3ah": (7, "Methyl 8-(4-Methoxyphenyl)-3-tosyl-3-azabicyclo[5.2.0]nona-1,6-diene-8-carboxylate", 91, 80.0, "flash column chromatography on silica gel (hexanes/EtOAc = 7:1)"),
    "3ai": ((7, 8), "Methyl 8-(4-Fluorophenyl)-3-tosyl-3-azabicyclo[5.2.0]nona-1,6-diene-8-carboxylate", 71, 60.7, "flash column chromatography on silica gel (hexanes/EtOAc = 7:1)"),
    "3aj": (8, "Methyl 8-(4-Chlorophenyl)-3-tosyl-3-azabicyclo[5.2.0]nona-1,6-diene-8-carboxylate", 68, 60.4, "flash column chromatography on silica gel (hexanes/EtOAc = 7:1)"),
    "3ak": (8, "Methyl 8-(4-Bromophenyl)-3-tosyl-3-azabicyclo[5.2.0]nona-1,6-diene-8-carboxylate", 66, 64.5, "flash column chromatography on silica gel (hexanes/EtOAc = 7:1)"),
    "3al": (8, "Methyl 8-(4-Iodophenyl)-3-tosyl-3-azabicyclo[5.2.0]nona-1,6-diene-8-carboxylate", 61, 65.3, "flash column chromatography on silica gel (hexanes/EtOAc = 7:1)"),
    "3am": (8, "Methyl 8-(3-Chlorophenyl)-3-tosyl-3-azabicyclo[5.2.0]nona-1,6-diene-8-carboxylate", 69, 61.3, "flash column chromatography on silica gel (hexanes/EtOAc = 7:1)"),
    "3an": (8, "Methyl 3-Tosyl-8-(3-(trifluoromethyl)phenyl)-3-azabicyclo[5.2.0]nona-1,6-diene-8-carboxylate", 53, 50.6, "flash column chromatography on silica gel (hexanes/EtOAc = 7:1)"),
    "3ao": (8, "Methyl 8-(o-Tolyl)-3-tosyl-3-azabicyclo[5.2.0]nona-1,6-diene-8-carboxylate", 23, 19.4, "flash column chromatography on silica gel (hexanes/EtOAc = 10:1 to 7:1)"),
    "3ap": (8, "Methyl 8-(3,4-Dichlorophenyl)-3-tosyl-3-azabicyclo[5.2.0]nona-1,6-diene-8-carboxylate", 60, 57.4, "flash column chromatography on silica gel (hexanes/EtOAc = 7:1)"),
    "3aq": (8, "Methyl 8-(3-Bromo-4-fluorophenyl)-3-tosyl-3-azabicyclo[5.2.0]nona-1,6-diene-8-carboxylate", 63, 63.8, "flash column chromatography on silica gel (hexanes/EtOAc = 7:1)"),
    "3ar": (8, "Methyl 8-(Naphthalen-2-yl)-3-tosyl-3-azabicyclo[5.2.0]nona-1,6-diene-8-carboxylate", 58, 53.3, "flash column chromatography on silica gel (hexanes/EtOAc = 7:1)"),
    "3as": (8, "Methyl 8-(Benzofuran-3-yl)-3-tosyl-3-azabicyclo[5.2.0]nona-1,6-diene-8-carboxylate", 69, 62.0, "flash column chromatography on silica gel (hexanes/EtOAc = 7:1)"),
    "3ba": (8, "Methyl 8-Phenyl-3-(phenylsulfonyl)-3-azabicyclo[5.2.0]nona-1,6-diene-8-carboxylate", 80, 63.2, "flash column chromatography on silica gel (hexanes/EtOAc = 7:1)"),
    "3ca": ((8, 9), "Methyl 3-((4-Bromophenyl)sulfonyl)-8-phenyl-3-azabicyclo[5.2.0]nona-1,6-diene-8-carboxylate", 76, 72.0, "flash column chromatography on silica gel (hexanes/EtOAc = 7:1)"),
    "3da": (9, "Methyl 3-((4-Nitrophenyl)sulfonyl)-8-phenyl-3-azabicyclo[5.2.0]nona-1,6-diene-8-carboxylate", 71, 62.5, "flash column chromatography on silica gel (hexanes/EtOAc = 7:1)"),
    "3ea": (9, "Methyl 3-(Methylsulfonyl)-8-phenyl-3-azabicyclo[5.2.0]nona-1,6-diene-8-carboxylate", 73, 48.7, "flash column chromatography on silica gel (hexanes/EtOAc = 7:1)"),
    "3fa": (9, "Methyl 5-Methyl-8-phenyl-3-tosyl-3-azabicyclo[5.2.0]nona-1,6-diene-8-carboxylate", 82, 69.4, "flash column chromatography on silica gel (hexanes/EtOAc = 7:1)"),
    "3ga": (9, "Methyl 4,8-Diphenyl-3-tosyl-3-azabicyclo[5.2.0]nona-1,6-diene-8-carboxylate", 73, 70.9, "flash column chromatography on silica gel (hexanes/EtOAc = 7:1)"),
    "3ha": (9, "Methyl 1-Phenyl-4-tosyl-2,4,4a,5,6,7,8,8a-octahydro-1H-benzo[b]cyclobuta[e]azepine-1-carboxylate", 35, 32.4, "flash column chromatography on silica gel (hexanes/EtOAc = 7:1), then recrystallization in hexanes/EtOAc"),
    "3ia": (9, "Methyl 9-Methyl-8-phenyl-3-tosyl-3-azabicyclo[5.2.0]nona-1,6-diene-8-carboxylate", 74, 62.7, "flash column chromatography on silica gel (hexanes/EtOAc = 7:1)"),
    "5a": (9, "Methyl (1Z,7E)-9-Phenyl-3-tosyl-3-azabicyclo[6.2.0]deca-1,7-diene-9-carboxylate", 61, 51.7, "flash column chromatography on silica gel (hexanes/EtOAc = 15:1 to 7:1)"),
    "5b": (9, "Methyl (1Z,7E)-9-(p-Tolyl)-3-tosyl-3-azabicyclo[6.2.0]deca-1,7-diene-9-carboxylate", 65, 57.0, "flash column chromatography on silica gel (hexanes/EtOAc = 15:1 to 7:1)"),
    "5c": (9, "Methyl (1Z,7E)-9-([1,1'-Biphenyl]-4-yl)-3-tosyl-3-azabicyclo[6.2.0]deca-1,7-diene-9-carboxylate", 68, 68.0, "flash column chromatography on silica gel (hexanes/EtOAc = 15:1 to 7:1)"),
    "5d": ((9, 10), "Methyl (1Z,7E)-9-(4-Fluorophenyl)-3-tosyl-3-azabicyclo[6.2.0]deca-1,7-diene-9-carboxylate", 60, 53.0, "flash column chromatography on silica gel (hexanes/EtOAc = 15:1 to 7:1)"),
    "5e": (10, "Methyl (1Z,7E)-9-(4-Chlorophenyl)-3-tosyl-3-azabicyclo[6.2.0]deca-1,7-diene-9-carboxylate", 58, 53.0, "flash column chromatography on silica gel (hexanes/EtOAc = 15:1 to 7:1)"),
    "5f": (10, "Methyl (1Z,7E)-9-(4-Bromophenyl)-3-tosyl-3-azabicyclo[6.2.0]deca-1,7-diene-9-carboxylate", 55, 55.3, "flash column chromatography on silica gel (hexanes/EtOAc = 15:1 to 7:1)"),
    "5g": (10, "Methyl (1Z,7E)-9-(4-Iodophenyl)-3-tosyl-3-azabicyclo[6.2.0]deca-1,7-diene-9-carboxylate", 52, 57.0, "flash column chromatography on silica gel (hexanes/EtOAc = 15:1 to 7:1)"),
    "5h": (10, "Methyl (1Z,7E)-9-(3-Chlorophenyl)-3-tosyl-3-azabicyclo[6.2.0]deca-1,7-diene-9-carboxylate", 38, 34.8, "flash column chromatography on silica gel (hexanes/EtOAc = 15:1 to 7:1)"),
    "5i": (10, "Methyl (1Z,7E)-9-(Naphthalen-2-yl)-3-tosyl-3-azabicyclo[6.2.0]deca-1,7-diene-9-carboxylate", 47, 44.5, "flash column chromatography on silica gel (hexanes/EtOAc = 15:1 to 7:1)"),
}

scheme_for = {}
for label in products:
    scheme_for[label] = "scheme2_" + label if label.startswith("3a") and label not in {"3ba"} else None
scheme_for.update({x: "scheme3_" + x for x in ("3ba", "3ca", "3da", "3ea", "3fa", "3ga", "3ha", "3ia")})
scheme_for.update({x: "scheme4_" + x for x in ("5a", "5b", "5c", "5d", "5e", "5f", "5g", "5h", "5i")})

for label, (pages, name, yld, mass, purification) in products.items():
    if isinstance(pages, int):
        pages = (pages,)
    evidence = [EV[p] for p in pages]
    mid = f"{DOC}:mol:{label}"
    touch_molecule(mid, evidence, name,
                   "Full source name and characterization paragraph inspected; no concrete SMILES was promoted without an independently checked drawing-to-graph transcription.")
    if label.startswith("5"):
        update_product_stereo(mid, evidence, geometry="(1Z,7E)")
    if label in {"3fa", "3ga"}:
        ratio = "approximately 1.5:1" if label == "3fa" else "approximately 2:1"
        update_product_stereo(mid, evidence,
                              mixture=f"mixture of two inseparable regioisomers; source also reports dr {ratio}")
    if label == "3ha":
        update_product_stereo(mid, evidence,
                              mixture="flash chromatography gave four inseparable regioisomers; recrystallization afforded one major regioisomer",
                              assignment="racemic product; relative stereochemistry reported as established by X-ray in the Results and Discussion",
                              chirality="racemic", representation="not_encoded")
    if label == "3ia":
        update_product_stereo(mid, evidence,
                              assignment="racemic product; relative stereochemistry reported as established by X-ray in the Results and Discussion",
                              chirality="racemic", representation="not_encoded")
    r = get_reaction(scheme_for[label])
    use_gpb_small(r, evidence, purification, mid, yld, mass)
    if label in {"3fa", "3ga"}:
        sel = find_outcome(r, "selectivity")
        if sel:
            add_unique(sel["value"]["evidence_ids"], evidence)
            add_unique(sel["measurement_basis"]["evidence_ids"], evidence)
            sel["target_molecule_ids"] = [mid]
    if label == "3ha":
        r["notes"].append(fact(
            "flash chromatography gave four inseparable regioisomers; recrystallization afforded one major regioisomer",
            "four-regioisomer mixture followed by isolation of one major regioisomer", evidence))
    if label == "3ia":
        r["notes"].append(fact("one regioisomer", "one regioisomer", evidence))
        r["notes"].append(fact(
            "Scheme 3 footnote d says 'One diastereomer was obtained', whereas the individual paragraph says 'one regioisomer'",
            None, [f"{DOC}:page5:image", EV[9]], status="conflict", alternatives=[
                fact("One diastereomer was obtained", "one diastereomer", [f"{DOC}:page5:image"]),
                fact("one regioisomer", "one regioisomer", [EV[9]])
            ]))
        add_unique(r["review"]["unresolved"], [
            "Source terminology conflict: Scheme 3 footnote d says one diastereomer, while the individual 3ia paragraph says one regioisomer; neither wording was silently repaired."
        ])
    incoming_reactions.append(r)


# Product 7 paragraph is the same Scheme 5b control experiment, now with the isolated mass and GP-B details.
ev10 = [EV[10]]
name7 = "Methyl 8-((4-Methyl-N-(propa-1,2-dien-1-yl)phenyl)sulfonamido)-2-phenylocta-2,3-dienoate"
mid7 = f"{DOC}:mol:7"
touch_molecule(mid7, ev10, name7,
               "The named product contains two allene units; axial configurations are unreported and no axial isomeric encoding was promoted.")
m7 = incoming_molecules[mid7]
m7["stereochemistry"]["chirality_status"] = fact("allene axial configurations are not reported", "unspecified", ev10)
m7["stereochemistry"]["geometry_status"] = fact("allene geometry/axial configuration not assigned", "unspecified", ev10)
m7["stereochemistry"]["representation_status"] = fact("no complete stereo-resolved SMILES", "unresolved", ev10)
r7 = get_reaction("scheme5_bisallene7")
use_gpb_small(r7, ev10, "flash column chromatography on silica gel (hexanes/EtOAc = 15:1 to 7:1)", mid7, 56, 49.0)
incoming_reactions.append(r7)


# GP-B general procedure record.
gpb_evidence = [EV[7]]
gpb = st.new_reaction(GPB_ID, DOC, gpb_evidence)
gpb["record_type"] = "generic_scheme"
gpb["experiment_categories"] = fact(["other"], ["other"], gpb_evidence)
gpb["location"] = fact("page 7, General Procedure for the Synthesis of Compounds 3, 5, and 7 (GP-B)",
                       {"pages": [7], "region": "General Procedure GP-B below Scheme 6"}, gpb_evidence)
gpb["reaction_name"] = fact("General Procedure for the Synthesis of Compounds 3, 5, and 7 (GP-B)",
                            "General Procedure for the Synthesis of Compounds 3, 5, and 7 (GP-B)", gpb_evidence)
gpb["transformation_type"] = fact(
    None, "copper-catalyzed allenyne/diazo cross-coupling followed by cycloaddition, or bisallene isolation",
    gpb_evidence, method="inference", derivation="Summary of the products explicitly covered by GP-B and the paper's reaction description.")
generic_allenyne = f"{DOC}:mol:generic_1n_allenyne"
generic_diazo = f"{DOC}:mol:generic_aryl_diazoacetate"
generic_product = f"{DOC}:mol:generic_azabicyclo_product"
cu = f"{DOC}:mol:CuPPh3Br"
ch2cl2 = f"{DOC}:mol:CH2Cl2"
for mid in (generic_allenyne, generic_diazo, generic_product, cu, ch2cl2):
    touch_molecule(mid, gpb_evidence)
gpb["participants"] = [
    st.new_participant("allenyne", generic_allenyne, "reactant", gpb_evidence),
    st.new_participant("diazo", generic_diazo, "reactant", gpb_evidence),
    st.new_participant("catalyst", cu, "catalyst", gpb_evidence),
    new_solvent("diazo_solution_solvent", ch2cl2, 1.0, gpb_evidence),
    new_solvent("reaction_solution_solvent", ch2cl2, 1.0, gpb_evidence),
    st.new_participant("product", generic_product, "product", gpb_evidence),
]
set_charge(participant(gpb, "allenyne"), "amount", "0.2 mmol", {"number": 0.2, "unit": "mmol"}, gpb_evidence)
set_charge(participant(gpb, "allenyne"), "equivalents", "1.0 equiv", {"number": 1.0, "unit": "equiv"}, gpb_evidence)
set_equiv_reference(participant(gpb, "allenyne"), generic_allenyne, gpb_evidence)
set_charge(participant(gpb, "diazo"), "amount", "0.4 mmol", {"number": 0.4, "unit": "mmol"}, gpb_evidence)
set_charge(participant(gpb, "diazo"), "equivalents", "2.0 equiv", {"number": 2.0, "unit": "equiv"}, gpb_evidence)
set_equiv_reference(participant(gpb, "diazo"), generic_allenyne, gpb_evidence)
set_charge(participant(gpb, "catalyst"), "mass", "18.6 mg", {"number": 18.6, "unit": "mg"}, gpb_evidence)
set_charge(participant(gpb, "catalyst"), "amount", "0.02 mmol", {"number": 0.02, "unit": "mmol"}, gpb_evidence)
set_charge(participant(gpb, "catalyst"), "mol_percent", "10 mol %", {"number": 10, "unit": "mol%"}, gpb_evidence)
set_equiv_reference(participant(gpb, "catalyst"), generic_allenyne, gpb_evidence)
gpb["reaction_scale"] = fact("1,n-allenyne (0.2 mmol, 1.0 equiv)",
                             {"reference_molecule_id": generic_allenyne, "amount": {"number": 0.2, "unit": "mmol"}}, gpb_evidence)
gpb["solvent_conditions"]["composition"] = fact(
    "CH2Cl2 (1.0 mL) in each of two explicitly reported solutions",
    [{"use": "diazo solution", "volume": {"number": 1.0, "unit": "mL"}},
     {"use": "allenyne/catalyst solution", "volume": {"number": 1.0, "unit": "mL"}}], gpb_evidence)
gpb["steps"] = [
    new_step("syringe_pump_addition", "reaction",
             ["allenyne", "diazo", "catalyst", "diazo_solution_solvent", "reaction_solution_solvent"],
             "inject the diazo solution into the allenyne/catalyst solution via syringe pump over 2 h",
             gpb_evidence, "r.t.", {"number": 2, "unit": "h"}),
    new_step("post_addition_hold", "reaction", ["allenyne", "diazo", "catalyst", "product"],
             "continue the resulting solution at room temperature for 12 h",
             gpb_evidence, "r.t.", {"number": 12, "unit": "h"}),
    new_step("concentration", "workup", [], "concentrate under reduced pressure", gpb_evidence),
    new_step("purification", "purification", ["product"],
             "flash column chromatography on silica gel (hexanes/EtOAc)", gpb_evidence),
]
gpb["procedure"]["original_text"] = fact(
    "A solution of aryl diazoacetate (0.4 mmol, 2.0 equiv) in CH2Cl2 (1.0 mL) was injected into a solution of 1,n-allenyne (0.2 mmol, 1.0 equiv) and Cu(PPh3)3Br (18.6 mg, 0.02 mmol, 10 mol %) in CH2Cl2 (1.0 mL) via a syringe pump over 2 h at rt. Then the resulting solution was continued at rt for 12 h. The mixture was concentrated under reduced pressure and purified by flash column chromatography (eluent: hexanes/EtOAc) on silica gel.",
    "GP-B", gpb_evidence)
gpb["procedure"]["addition_order"] = fact(
    "diazo solution injected into allenyne/catalyst solution",
    "diazo solution -> allenyne/catalyst solution", gpb_evidence)
gpb["classification"] = copy.deepcopy(old_reactions[f"{DOC}:scheme1e_generic"]["classification"])
gpb["classification"]["coverage"] = fact(
    "GP-B covers bicyclic products 3/5 and bisallene 7", "partial", gpb_evidence,
    method="inference", derivation="Products governed by one procedure have different net bond inventories; no single complete mapping applies.")
gpb["classification"]["limitations"] = [
    "The general procedure covers both cyclized products and the noncyclized bisallene control; formed-bond count varies by instance.",
    "No complete atom mapping was asserted."
]
for bf in gpb["classification"]["bond_changes"]:
    bf["evidence_ids"] = gpb_evidence
for key in ("labels",):
    gpb["classification"][key]["evidence_ids"] = gpb_evidence
gpb["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
                 "unresolved": ["Generic products 3, 5, and 7 do not share one fully mapped bond inventory; instance records carry the source-specific outcomes."]}
incoming_reactions.append(gpb)


# Detailed gram-scale repeat.
scale = get_reaction("scheme2_3aa_scaleup")
add_reaction_evidence(scale, ev10)
scale["general_scheme_ids"] = list(dict.fromkeys(scale["general_scheme_ids"] + [GPB_ID]))
scale["participants"] = [p for p in scale["participants"] if p["participant_id"] != "solvent"]
sub = participant(scale, "substrate"); dz = participant(scale, "diazo"); cat = participant(scale, "catalyst")
set_charge(sub, "mass", "1306.5 mg", {"number": 1306.5, "unit": "mg"}, ev10)
set_charge(sub, "amount", "5.0 mmol", {"number": 5.0, "unit": "mmol"}, ev10)
set_charge(sub, "equivalents", "1.0 equiv", {"number": 1.0, "unit": "equiv"}, ev10); set_equiv_reference(sub, sub["molecule_id"], ev10)
set_charge(dz, "mass", "1761.8 mg", {"number": 1761.8, "unit": "mg"}, ev10)
set_charge(dz, "amount", "10.0 mmol", {"number": 10.0, "unit": "mmol"}, ev10)
set_charge(dz, "equivalents", "2.0 equiv", {"number": 2.0, "unit": "equiv"}, ev10); set_equiv_reference(dz, sub["molecule_id"], ev10)
set_charge(cat, "mass", "465.2 mg", {"number": 465.2, "unit": "mg"}, ev10)
set_charge(cat, "amount", "0.5 mmol", {"number": 0.5, "unit": "mmol"}, ev10)
set_charge(cat, "mol_percent", "10 mol %", {"number": 10, "unit": "mol%"}, ev10); set_equiv_reference(cat, sub["molecule_id"], ev10)
scale["participants"].append(new_solvent("diazo_solution_solvent", f"{DOC}:mol:CH2Cl2", 25.0, ev10))
scale["participants"].append(new_solvent("reaction_solution_solvent", f"{DOC}:mol:CH2Cl2", 25.0, ev10))
pids = [p["participant_id"] for p in scale["participants"] if p["stage"] == "reaction"]
scale["steps"] = [
    new_step("dropwise_addition", "reaction", pids,
             "phenyl diazoacetate solution was added dropwise into the allenyne/catalyst solution via a constant pressure dropping funnel over 2 h",
             ev10, "r.t.", {"number": 2, "unit": "h"}, "N2"),
    new_step("post_addition_hold", "reaction", pids,
             "continue at room temperature for 12 h (monitored by TLC)", ev10, "r.t.", {"number": 12, "unit": "h"}, "N2"),
    new_step("concentration", "workup", [], "concentrate under reduced pressure", ev10),
    new_step("purification", "purification", ["product"],
             "flash column chromatography on silica gel (hexanes/EtOAc = 7:1)", ev10),
]
scale["reaction_scale"] = fact("allenyne 1a: 5.0 mmol", {"reference_molecule_id": sub["molecule_id"], "amount": {"number": 5.0, "unit": "mmol"}}, ev10)
scale["solvent_conditions"]["composition"] = fact(
    "CH2Cl2 (25.0 mL) for diazo solution; CH2Cl2 (25.0 mL) for allenyne/catalyst solution",
    [{"use": "diazo solution", "volume": {"number": 25.0, "unit": "mL"}},
     {"use": "allenyne/catalyst solution", "volume": {"number": 25.0, "unit": "mL"}}], ev10)
scale["solvent_conditions"]["total_volume"] = fact()
scale["procedure"]["original_text"] = fact(
    "A solution of phenyl diazoacetate 2a (1761.8 mg, 10.0 mmol, 2.0 equiv) in CH2Cl2 (25.0 mL) was added dropwise into a solution of allenyne 1a (1306.5 mg, 5.0 mmol, 1.0 equiv) and Cu(PPh3)3Br (465.2 mg, 0.5 mmol, 10 mol %) in CH2Cl2 (25.0 mL) via a constant pressure dropping funnel over 2 h at rt under N2; the mixture was continued at rt for 12 h, concentrated, and purified by flash chromatography (hexanes/EtOAc = 7:1).",
    "gram-scale procedure", ev10)
scale["procedure"]["modifications"] = fact("constant pressure dropping funnel under N2", "constant-pressure dropwise addition under N2", ev10)
scale["procedure"]["addition_order"] = fact("2a solution added to 1a/catalyst solution", "2a solution -> 1a/catalyst solution", ev10)
set_isolated_outcome(scale, f"{DOC}:mol:3aa", 72, 1474.2, ev10,
                     "flash column chromatography on silica gel (hexanes/EtOAc = 7:1)")
scale["review"]["unresolved"] = [u for u in scale["review"]["unresolved"] if "pages 1-6" not in u and "Yield basis" not in u]
add_unique(scale["review"]["unresolved"], ["Machine record remains needs_review; the source does not report an independent analytical yield method or uncertainty."])
incoming_reactions.append(scale)


# Detailed BHT/TEMPO radical-inhibition repeats.
for additive, mass, yld, product_mass in (("BHT", 44.1, 80, 65.5), ("TEMPO", 31.2, 84, 68.8)):
    r = get_reaction(f"scheme5_{additive}")
    use_gpb_small(r, ev10, "flash column chromatography on silica gel (hexanes/EtOAc = 7:1)",
                  f"{DOC}:mol:3aa", yld, product_mass)
    sub = participant(r, "substrate"); dz = participant(r, "diazo"); cat = participant(r, "catalyst")
    set_charge(sub, "mass", "52.3 mg", {"number": 52.3, "unit": "mg"}, ev10)
    set_charge(dz, "mass", "70.5 mg", {"number": 70.5, "unit": "mg"}, ev10)
    addp = participant(r, "additive1")
    set_charge(addp, "mass", f"{mass:g} mg", {"number": mass, "unit": "mg"}, ev10)
    set_charge(addp, "amount", "0.2 mmol", {"number": 0.2, "unit": "mmol"}, ev10)
    set_charge(addp, "equivalents", "1.0 equiv", {"number": 1.0, "unit": "equiv"}, ev10)
    set_equiv_reference(addp, sub["molecule_id"], ev10)
    r["procedure"]["modifications"] = fact(f"{additive} ({mass:g} mg, 0.2 mmol, 1.0 equiv) was present in the allenyne/catalyst solution",
                                             f"{additive}, 1.0 equiv", ev10)
    add_unique(r["review"]["unresolved"], ["The experiment supports persistence of product formation in the presence of the inhibitor; no radical-trapping detection limit is reported."])
    incoming_reactions.append(r)


# Scheme 6 proposed mechanism: update the existing prose record rather than creating experiments for arrows.
mech = get_reaction("proposed_mechanism_text")
mech_ev = [EV[7]]
add_reaction_evidence(mech, mech_ev)
intermediates = {
    "A": "author-proposed copper-carbene intermediate A",
    "B": "author-proposed migratory-insertion intermediate B",
    "C": "author-proposed bisallene intermediate C",
    "D": "author-proposed cupracyclopentane intermediate D",
    "CuI_generic": "Cu(I) species in Scheme 6",
    "Hplus": "H+ in Scheme 6",
    "HBr": "HBr shown in Scheme 6",
}
for label, name in intermediates.items():
    make_named_unresolved(label, name, mech_ev,
                          "Mechanistic species from the author-proposed cycle; not isolated and not promoted to a verified concrete structure.")
touch_molecule(f"{DOC}:mol:N2", mech_ev)
existing = {p["participant_id"] for p in mech["participants"]}
for pid, label, role in (("A", "A", "unknown"), ("B", "B", "unknown"), ("C", "C", "unknown"),
                         ("D", "D", "unknown"), ("CuI", "CuI_generic", "unknown"),
                         ("Hplus", "Hplus", "unknown"), ("N2", "N2", "product"), ("HBr", "HBr", "product")):
    if pid not in existing:
        mech["participants"].append(st.new_participant(pid, f"{DOC}:mol:{label}", role, mech_ev))
mech["steps"] = [
    new_step("carbene_formation", "reaction", ["1a", "2a", "CuI", "A", "N2", "HBr"],
             "Scheme 6 proposes formation of copper-carbene intermediate A from 1a and 2a, with N2 and HBr shown", mech_ev),
    new_step("migratory_insertion", "reaction", ["A", "B"],
             "migratory insertion of the carbenic carbon into the alkynyl group gives intermediate B", mech_ev),
    new_step("protonation", "reaction", ["B", "Hplus", "C"],
             "protonation of B gives bisallene intermediate C", mech_ev),
    new_step("cupracycle_formation", "reaction", ["C", "CuI", "D"],
             "regioselective coordination followed by C-C bond formation gives cupracyclopentane D", mech_ev),
    new_step("reductive_elimination", "reaction", ["D", "3aa", "CuI"],
             "reductive elimination gives 3aa and regenerates Cu(I)", mech_ev),
]
mech["procedure"]["original_text"] = fact("Scheme 6. Proposed Reaction Mechanism", "author-proposed mechanism", mech_ev)
mech["notes"].append(fact("Proposed Reaction Mechanism", "author-proposed; not an isolated experimental sequence", mech_ev))
mech["review"]["unresolved"] = [
    "Author-proposed mechanism only; intermediates A-D were not reported as isolated experiments.",
    "Mechanistic structures and copper oxidation-state drawings were retained as labels/descriptions rather than promoted to complete verified SMILES."
]
incoming_reactions.append(mech)


# Add every registry molecule referenced by an incoming reaction to the import bundle.
for reaction in incoming_reactions:
    for p in reaction["participants"]:
        touch_molecule(p["molecule_id"], reaction["evidence_ids"])


template["molecules"] = list(incoming_molecules.values())
template["reactions"] = incoming_reactions
template["source_context"] = [
    {"context_id": f"{DOC}:context:pages7_10_characterization", "document_id": DOC,
     "evidence_ids": [EV[7], EV[8], EV[9], EV[10]], "kind": "characterization",
     "note": "Individual product paragraphs contain 1H/13C NMR, HRMS, melting-point and physical-form characterization. These observations support the associated preparation records but are not separate synthesis experiments."},
    {"context_id": f"{DOC}:context:page10_unattached_si_notice", "document_id": DOC,
     "evidence_ids": [EV[10]], "kind": "supporting_information_notice",
     "note": "The article states that NMR spectra and crystal data are available as Supporting Information, but no SI file was supplied in this workspace; no supplementary document was associated or inferred."},
    {"context_id": f"{DOC}:context:pages11_12_back_matter", "document_id": DOC,
     "evidence_ids": [EV[11], EV[12]], "kind": "nonreaction_back_matter",
     "note": "Author information, acknowledgments and cited references/footnotes. Cited experiments were not imported as experiments of this paper."},
]

dump(RUN / "submitted.json", template)


coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 7, "state": "attempted", "regions": [
            {"region": "upper-left Scheme 6", "finding": "Author-proposed Cu(I) catalytic cycle A-D; merged with the existing proposed-mechanism record and not treated as five experimental reactions.", "evidence_ids": [EV[7]]},
            {"region": "left column below Scheme 6", "finding": "GP-B general procedure with two explicit CH2Cl2 portions, 2 h syringe-pump addition, 12 h room-temperature hold, concentration and flash chromatography; added as one generic-procedure record.", "evidence_ids": [EV[7]]},
            {"region": "remaining columns, product paragraphs 3aa-3ai", "finding": "Individual preparations give isolated material wording, masses, yields and chromatography conditions; merged into the corresponding Scheme 2 experiments. 3ai continues on page 8.", "evidence_ids": [EV[7]]},
            {"region": "NMR/HRMS blocks", "finding": "Characterization only; retained in source_context rather than creating synthesis records.", "evidence_ids": [EV[7]]}
        ], "record_ids": [f"{DOC}:proposed_mechanism_text", GPB_ID] + [f"{DOC}:scheme2_3a{x}" for x in "abcdefghi"],
         "unresolved": ["Mechanistic intermediates A-D are proposed and not isolated; their complete structures/stereochemistry were not promoted to verified SMILES.", "Product characterization was not transcribed peak-by-peak because it is not a separate synthesis experiment."]},
        {"page": 8, "state": "attempted", "regions": [
            {"region": "top-left continuation of 3ai", "finding": "Completes the 3ai characterization paragraph begun on page 7; same experiment, merged evidence.", "evidence_ids": [EV[8]]},
            {"region": "product paragraphs 3aj-3as", "finding": "Individual Scheme 2 preparations with reported product masses, yields and flash chromatography conditions; merged with existing scope records.", "evidence_ids": [EV[8]]},
            {"region": "lower-right 3ba and start of 3ca", "finding": "Scheme 3 product preparations; 3ca continues on page 9 and is one experiment, not duplicated.", "evidence_ids": [EV[8]]},
            {"region": "NMR/HRMS and melting points", "finding": "Characterization context only; no independent reaction records.", "evidence_ids": [EV[8]]}
        ], "record_ids": [f"{DOC}:scheme2_3ai"] + [f"{DOC}:scheme2_3a{x}" for x in "jklmnopqrs"] + [f"{DOC}:scheme3_3ba", f"{DOC}:scheme3_3ca"],
         "unresolved": ["Concrete product SMILES remain blank because the named/drawn structures were not independently graph-transcribed and verified."]},
        {"page": 9, "state": "attempted", "regions": [
            {"region": "upper-left continuation of 3ca", "finding": "Completes the 3ca paragraph begun on page 8; merged with the same Scheme 3 experiment.", "evidence_ids": [EV[9]]},
            {"region": "left column 3da-3ha", "finding": "Individual preparations with isolated masses/yields; 3fa/3ga retain inseparable-regioisomer wording and reported dr, while 3ha retains the four-regioisomer/recrystallization sequence.", "evidence_ids": [EV[9]]},
            {"region": "right column 3ia and 5a-5d", "finding": "3ia terminology conflict (one regioisomer in paragraph versus one diastereomer in Scheme 3 footnote) retained; 5a-5d names explicitly specify (1Z,7E), with 5d continuing on page 10.", "evidence_ids": [EV[9]]},
            {"region": "spectral blocks", "finding": "NMR/HRMS/melting-point characterization retained as source context, not new synthesis experiments.", "evidence_ids": [EV[9]]}
        ], "record_ids": [f"{DOC}:scheme3_{x}" for x in ("3ca", "3da", "3ea", "3fa", "3ga", "3ha", "3ia")] + [f"{DOC}:scheme4_5{x}" for x in "abcd"],
         "unresolved": ["The source's regioisomer/diastereomer terminology for 3ia is internally inconsistent and remains unresolved.", "Absolute configurations for racemic 3ha/3ia are not encoded in registry SMILES."]},
        {"page": 10, "state": "attempted", "regions": [
            {"region": "top-left continuation of 5d and product paragraphs 5e-5i", "finding": "Completes 5d and records individual 1,7-allenyne products 5e-5i; all names explicitly retain (1Z,7E), masses, yields and chromatography ranges.", "evidence_ids": [EV[10]]},
            {"region": "upper-right product 7 paragraph", "finding": "Bisallene 7 preparation under GP-B with 49.0 mg and 56% yield; merged with the existing Scheme 5b mechanistic control.", "evidence_ids": [EV[10]]},
            {"region": "right-column Gram-Scale Synthesis", "finding": "Detailed 5.0 mmol repeat of 3aa with exact masses, two 25.0 mL CH2Cl2 portions, N2, 2 h dropwise addition, 12 h hold, chromatography, 1474.2 mg and 72% yield; merged with existing scale-up record.", "evidence_ids": [EV[10]]},
            {"region": "right-column Radical Inhibition Experiments", "finding": "Distinct BHT and TEMPO experiments with additive masses/amounts and isolated product masses/yields; merged with the two existing mechanistic-study controls.", "evidence_ids": [EV[10]]},
            {"region": "Associated Content", "finding": "SI availability notice only. No SI file is present in the workspace, so no supplementary source was associated.", "evidence_ids": [EV[10]]}
        ], "record_ids": [f"{DOC}:scheme4_5{x}" for x in "defghi"] + [f"{DOC}:scheme5_bisallene7", f"{DOC}:scheme2_3aa_scaleup", f"{DOC}:scheme5_BHT", f"{DOC}:scheme5_TEMPO"],
         "unresolved": ["The source reports yields after purification but does not provide an independent analytical yield method or uncertainty.", "No supplied SI file was associated despite the article's SI notice."]},
        {"page": 11, "state": "attempted", "regions": [
            {"region": "left column author information, contributions, acknowledgments", "finding": "Nonreaction back matter; no experiment record created.", "evidence_ids": [EV[11]]},
            {"region": "references 1-13", "finding": "Cited literature and a note listing other screened metal catalysts. These are background citations/footnote context, not newly detailed experiments of this paper on this page.", "evidence_ids": [EV[11]]}
        ], "record_ids": [],
         "unresolved": ["The listed alternative catalysts refer back to Table 1/background screening but do not supply row-specific conditions or outcomes beyond the already extracted main-paper records."]},
        {"page": 12, "state": "attempted", "regions": [
            {"region": "references and footnotes 14-22", "finding": "Other-solvent list, CCDC identifiers, mechanistic literature citations and precursor literature references; no new synthesis experiment of this paper is reported.", "evidence_ids": [EV[12]]}
        ], "record_ids": [],
         "unresolved": ["CCDC and cited external papers were not imported as source experiments; no CIF or SI files were supplied in this workspace."]},
    ]
}
dump(RUN / "coverage/pages_0007_0012.json", coverage)


# Mark only this bounded page task range attempted.
tasks_path = RUN / "tasks.json"
if tasks_path.exists():
    tasks = json.loads(tasks_path.read_text())
    for task in tasks:
        match = re.search(r":page(\d+):", task.get("evidence_id", ""))
        if match and 7 <= int(match.group(1)) <= 12:
            task["state"] = "attempted"
            task["reason"] = "Original page image directly inspected for the bounded PDF pages 7-12 pass; see coverage/pages_0007_0012.json."
    dump(tasks_path, tasks)


# Literal mention hints are retained as raw workflow support, not used as automatic identity merges.
mentions = st.source_mentions(json.loads((RUN / "dataset.json").read_text()),
                              ["3aa", "3ia", "5a", "7", "GP-B", "Gram-Scale Synthesis"], context=180)
dump(RUN / "raw_model/source_mentions_pages_0007_0012.json", mentions)
dump(RUN / "raw_model/visual_transcription_pages_0007_0012.json", {
    "document_id": DOC,
    "pages": [7, 8, 9, 10, 11, 12],
    "source": "direct visual inspection of the supplied original-detail page images",
    "new_reaction_ids": [GPB_ID],
    "updated_reaction_ids": [r["reaction_id"] for r in incoming_reactions if r["reaction_id"] != GPB_ID],
    "updated_reaction_count": len(incoming_reactions) - 1,
    "incoming_bundle_reaction_count": len(incoming_reactions),
    "inventory": {
        "scheme6_proposed_mechanism": 1,
        "GP_B_generic": 1,
        "individual_product_preparations_3aa_to_3as": 19,
        "individual_product_preparations_scheme3": 8,
        "individual_product_preparations_5a_to_5i": 9,
        "bisallene7_control": 1,
        "gram_scale_repeat": 1,
        "radical_inhibition_controls": 2,
        "new_specific_experiments_after_cross_page_deduplication": 0,
        "nonreaction_pages": [11, 12]
    },
    "limitations": [
        "No SI/CIF file was supplied or associated.",
        "Full product SMILES were not promoted without independently checked drawing-to-graph transcription.",
        "All imported machine records remain needs_review."
    ]
})

print(json.dumps({
    "incoming_reactions": len(incoming_reactions),
    "new_reactions": 1,
    "updated_reactions": len(incoming_reactions) - 1,
    "incoming_molecules": len(incoming_molecules),
    "coverage_pages": [7, 8, 9, 10, 11, 12]
}, ensure_ascii=False))
