#!/usr/bin/env python3
"""Build source-derived v3.1 records for JACS 2021, 143, 10374-10381, pages 1-6.

This is a transcription/build helper, not a chemistry verifier. Dense structures remain
unencoded when the page drawing cannot be transcribed with stereochemical confidence.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DATASET = RUN / "dataset.json"
DOC = "doc_9508ed0ca7d13d60"
EV = {p: f"{DOC}:page{p}:image" for p in range(1, 7)}
TXT = {p: f"{DOC}:page{p}:text" for p in range(1, 7)}


def dump(path: Path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n")


def f(raw=None, normalized=None, evidence=(), method="original_extraction",
      status=None, derivation=None):
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   derivation=derivation)


def missing(status="unreported"):
    return f(status=status, method="not_applicable" if status == "not_applicable" else "original_extraction")


data = json.loads(DATASET.read_text())
doc = data["documents"][0]
doc["title"] = "Allenone-Mediated Racemization/Epimerization-Free Peptide Bond Formation and Its Application in Peptide Synthesis"
doc["supplementary_material"] = {
    "status": "mentioned_not_supplied",
    "note": "Pages 1, 3, and 7 mention Supporting Information, but no SI file was supplied or explicitly associated in this session."
}
if not any(x.get("event") == "source_metadata_transcribed_pages_1_6" for x in data["run_log"]):
    data["run_log"].append({"event": "source_metadata_transcribed_pages_1_6", "document_id": DOC})
dump(DATASET, data)

bundle = json.loads((RUN / "submitted.json").read_text())
bundle["documents"] = data["documents"]
bundle["molecules"] = []
bundle["reactions"] = []
bundle["source_context"] = []

molecules = {}


def unresolved_molecule(label, description, evidence, stereo="unresolved"):
    mid = f"{DOC}:mol:{label}"
    if mid in molecules:
        occ = molecules[mid]["occurrences"][0]
        occ["evidence_ids"] = sorted(set(occ["evidence_ids"] + list(evidence)))
        occ["compound_label"]["evidence_ids"] = sorted(set(occ["compound_label"]["evidence_ids"] + list(evidence)))
        return mid
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    m["names"]["original"] = f(description, description, evidence)
    m["names"]["english"] = f(description, description, evidence)
    m["representation_type"] = "unresolved"
    m["smiles"] = missing()
    m["generic_template"] = None
    m["cas_numbers"] = missing()
    for key in ("atom_centers", "double_bonds", "axial", "mixture", "assignment",
                "isomeric_smiles", "identity_key"):
        m["stereochemistry"][key] = missing()
    if stereo == "drawn":
        m["stereochemistry"]["chirality_status"] = f(
            "stereochemical bonds are shown in the source drawing; absolute descriptors are not printed",
            "partially_specified", evidence)
        m["stereochemistry"]["assignment"] = f(
            "wedge/dash stereochemistry shown; descriptor assignment not transcribed",
            {"source_representation": "wedge/dash", "absolute_descriptors": None}, evidence)
    else:
        m["stereochemistry"]["chirality_status"] = f(
            None, "unresolved", evidence, method="inference",
            derivation="The source label/drawing was retained, but a full stereochemical graph was not encoded.")
    m["stereochemistry"]["geometry_status"] = f(
        None, "unresolved", evidence, method="inference",
        derivation="E/Z and ring-relative stereochemistry were not assigned from the dense page drawing.")
    m["stereochemistry"]["representation_status"] = f(
        None, "not_encoded", evidence, method="inference",
        derivation="Complete stereochemically faithful SMILES was not asserted from the page image.")
    m["identity_note"] = f(
        description,
        "Source-identified label/drawing retained; complete structure requires expert review.",
        evidence)
    molecules[mid] = m
    return mid


def simple_molecule(label, description, smiles, evidence, chirality="achiral", geometry="not_applicable"):
    mid = f"{DOC}:mol:{label}"
    if mid in molecules:
        occ = molecules[mid]["occurrences"][0]
        occ["evidence_ids"] = sorted(set(occ["evidence_ids"] + list(evidence)))
        return mid
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    m["names"]["original"] = f(description, description, evidence)
    m["names"]["english"] = f(description, description, evidence)
    m["cas_numbers"] = missing()
    m["representation_type"] = "smiles"
    m["smiles"] = f(
        f"structure shown/named as {description}", smiles, evidence,
        method="structure_derivation",
        derivation="Candidate connectivity was transcribed from the source drawing/name; parsing is diagnostic and does not verify identity.")
    m["generic_template"] = None
    for key in ("atom_centers", "double_bonds", "axial", "mixture", "assignment",
                "isomeric_smiles", "identity_key"):
        m["stereochemistry"][key] = missing()
    m["stereochemistry"]["chirality_status"] = f(
        None, chirality, evidence, method="structure_derivation",
        derivation="Candidate graph and source drawing contain no assigned tetrahedral or axial stereochemistry." if chirality == "achiral" else "Candidate graph/source drawing leaves potential stereochemistry unspecified.")
    m["stereochemistry"]["geometry_status"] = f(
        None, geometry, evidence, method="structure_derivation",
        derivation="Candidate graph/source drawing was checked for specified alkene geometry.")
    m["stereochemistry"]["representation_status"] = f(
        None, "not_applicable" if chirality == "achiral" and geometry == "not_applicable" else "partial",
        evidence, method="structure_derivation",
        derivation="No stereochemical encoding is applicable to the candidate graph." if chirality == "achiral" and geometry == "not_applicable" else "Connectivity is encoded, but unassigned stereogenic elements remain.")
    m["identity_note"] = f(
        description,
        "Source identity with a machine candidate SMILES; structure remains needs_review.",
        evidence)
    molecules[mid] = m
    return mid


def generic_molecule(label, description, fixed, variables, evidence):
    mid = f"{DOC}:mol:{label}"
    if mid in molecules:
        return mid
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    m["names"]["original"] = f(description, description, evidence)
    m["names"]["english"] = f(description, description, evidence)
    m["cas_numbers"] = missing()
    m["representation_type"] = "generic_template"
    m["smiles"] = missing()
    var_objs = []
    for var in variables:
        var_objs.append({
            "label": var,
            "attachment_points": f(f"source attachment for {var}", {"label": var}, evidence),
            "definition": missing(), "allowed_values": missing(), "constraints": missing()
        })
    m["generic_template"] = {
        "fixed_structure": f(fixed, {"source_description": fixed}, evidence),
        "variables": var_objs, "instance_mapping": [],
        "reason": f(description, "Source is an R-group/general structure, not one concrete molecule.", evidence)
    }
    for key in ("atom_centers", "double_bonds", "axial", "mixture", "assignment",
                "isomeric_smiles", "identity_key"):
        m["stereochemistry"][key] = missing()
    for key in ("chirality_status", "geometry_status", "representation_status"):
        m["stereochemistry"][key] = f(
            None, "unresolved", evidence, method="inference",
            derivation="Generic R-group structure cannot establish a single stereochemical state.")
    m["identity_note"] = f(description, "Generic source template retained without wildcard SMILES.", evidence)
    molecules[mid] = m
    return mid


# Shared identities and templates.
M2A = simple_molecule("2a", "allenone 2a", "C=C=CC(=O)c1ccccc1", [EV[1], EV[2], EV[3], EV[4], EV[5], TXT[6]], geometry="unresolved")
M1A = simple_molecule("1a", "benzoic acid 1a", "O=C(O)c1ccccc1", [EV[2]])
MPHEN = simple_molecule("phenethylamine", "2-phenylethylamine", "NCCc1ccccc1", [EV[2], EV[3]])
M3A = simple_molecule("3a", "alpha-carbonyl vinyl ester 3a", "CC(OC(=O)c1ccccc1)=CC(=O)c1ccccc1", [EV[2], EV[3]], chirality="achiral", geometry="unspecified")
M4A = simple_molecule("4a", "N-(2-phenylethyl)benzamide 4a", "O=C(NCCc1ccccc1)c1ccccc1", [EV[2], EV[3]])
M5 = simple_molecule("5", "1-phenyl-1,3-butanedione byproduct 5", "CC(=O)CC(=O)c1ccccc1", [EV[2], EV[3]])
M5P = unresolved_molecule("5prime", "trace condensation product 5'", [EV[2]])
MDCE = simple_molecule("DCE", "1,2-dichloroethane", "ClCCCl", [EV[2], EV[3], EV[4], EV[5]])
MDMF = simple_molecule("DMF", "N,N-dimethylformamide", "CN(C)C=O", [EV[2], EV[3], EV[4], EV[5], EV[6]])
MHOBT = unresolved_molecule("HOBt", "1-hydroxybenzotriazole (HOBt)", [EV[3], EV[4], EV[6]])
MTFA = simple_molecule("TFA", "trifluoroacetic acid", "O=C(O)C(F)(F)F", [EV[5], EV[6]])
MDCM = simple_molecule("DCM", "dichloromethane", "ClCCl", [EV[5]])
MTHF = simple_molecule("THF", "tetrahydrofuran", "C1CCOC1", [EV[5]])
MKI = simple_molecule("KI", "potassium iodide", "[K+].[I-]", [EV[5]])
MMOR = simple_molecule("morpholine", "morpholine", "C1COCCN1", [EV[5]])
MTIS = unresolved_molecule("TIS", "triisopropylsilane (TIS)", [EV[6]])
MPIP = simple_molecule("piperidine", "piperidine", "N1CCCCC1", [EV[6]])

GEN_ACID = generic_molecule("generic_acid", "generic carboxylic acid 1", "R1-C(=O)-OH", ["R1"], [EV[2], EV[3]])
GEN_AMINE = generic_molecule("generic_amine", "generic amine 6", "H-N(R2)-R3", ["R2", "R3"], [EV[3]])
GEN_ESTER = generic_molecule("generic_active_ester", "generic alpha-carbonyl vinyl ester 3", "R1-C(=O)-O-C(Me)=CH-C(=O)-Ph", ["R1"], [EV[2], EV[3]])
GEN_AMIDE = generic_molecule("generic_amide", "generic amide 4", "R1-C(=O)-N(R2)-R3", ["R1", "R2", "R3"], [EV[2], EV[3]])


def participant(pid, mid, role, evidence, amount=None, volume=None, equiv=None, molpct=None, ref=None,
                stage="reaction", product_role=None):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    if product_role:
        p["product_role"] = f(product_role, product_role, evidence)
    if amount is not None:
        p["charges"]["amount"] = f(amount[0], {"number": amount[1], "unit": amount[2]}, evidence)
    if volume is not None:
        p["charges"]["volume"] = f(volume[0], {"number": volume[1], "unit": volume[2]}, evidence)
    if equiv is not None:
        p["charges"]["equivalents"] = f(equiv[0], {"number": equiv[1], "unit": "equiv"}, evidence)
    if molpct is not None:
        p["charges"]["mol_percent"] = f(molpct[0], {"number": molpct[1], "unit": "mol%"}, evidence)
    if ref is not None:
        p["equivalence_reference"] = f(ref, ref, evidence)
    return p


def step(step_id, pids, evidence, operation, temperature=None, time=None, stage="reaction"):
    s = st.new_step(step_id, evidence, stage=stage)
    s["participant_ids"] = list(pids)
    s["operation"] = f(operation, operation, evidence)
    if temperature is not None:
        norm = temperature if isinstance(temperature, dict) else temperature
        s["temperature"] = f(str(temperature), norm, evidence)
    if time is not None:
        norm = time if isinstance(time, dict) else time
        s["time"] = f(str(time), norm, evidence)
    return s


def classification(label, evidence, count=1, limitation=None):
    if label is None:
        return {
            "standard_id": st.STANDARD, "mapped_reaction_smiles": missing(),
            "bond_changes": [], "labels": f(None, [], evidence, method="inference", derivation="No observed bond formation is classified for this contextual/control record."),
            "coverage": f(None, "not_applicable", evidence, method="inference", derivation="No applicable observed transformation was assigned."),
            "limitations": [limitation or "No atom-mapped observed transformation."], "model_votes": []
        }
    elems = {"C-O": ("C", "O", 3), "C-N": ("C", "N", 5)}[label]
    changes = []
    ids = []
    for i in range(count):
        cid = f"b{i+1}"
        ids.append(cid)
        changes.append(f(None, {
            "change_id": cid, "atoms": [{"atom_map": None, "element": elems[0]}, {"atom_map": None, "element": elems[1]}],
            "before_order": 0, "after_order": 1, "change_type": "formed"
        }, evidence, method="structure_derivation",
        derivation="Net bond type inferred from the source reactant/product drawings; atoms are not fully mapped."))
    labels = [{"category_number": elems[2], "code": label, "supporting_change_ids": ids}]
    return {
        "standard_id": st.STANDARD, "mapped_reaction_smiles": missing(),
        "bond_changes": changes,
        "labels": f(None, labels, evidence, method="structure_derivation",
                    derivation="Taxonomy label derived from the source transformation drawing; full atom mapping is unresolved."),
        "coverage": f(None, "partial", evidence, method="inference",
                      derivation="Principal formed bond retained; bond-order changes, leaving groups, hydrogens, and complete atom mapping remain unresolved."),
        "limitations": [limitation or "Principal net bond only; no complete atom mapping or hydrogen-source assignment."],
        "model_votes": []
    }


def reaction(rid, evidence, location, name, transformation, categories, cls_label,
             participants=(), steps=(), outcomes=(), notes=(), unresolved=(), record_type="specific_experiment"):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, evidence)
    r["location"] = f(location, location, evidence)
    r["reaction_name"] = f(name, name, evidence)
    r["transformation_type"] = f(None, transformation, evidence, method="inference",
                                 derivation="Transformation type summarized from the source arrow and product identity.")
    r["participants"] = list(participants)
    r["steps"] = list(steps)
    r["outcomes"] = list(outcomes)
    r["classification"] = classification(cls_label, evidence)
    r["notes"] = list(notes)
    r["review"]["unresolved"] = list(unresolved)
    r["completeness"] = "partial"
    bundle["reactions"].append(r)
    return r


def outcome(oid, kind, target_mid, raw, normalized, evidence, basis=None, scope=None, yield_type=None):
    o = st.new_outcome(oid, kind, [target_mid], evidence)
    o["value"] = f(raw, normalized, evidence)
    if basis is not None:
        o["measurement_basis"] = f(basis, basis, evidence)
    if scope is not None:
        o["scope"] = f(scope, scope, evidence)
    if yield_type is not None:
        o["yield_type"] = f(yield_type, yield_type, evidence)
    return o


def set_common_solvent_and_scale(r, reference_mid, amount, solvent_mid, volume, evidence):
    r["reaction_scale"] = f(
        f"{reference_mid} {amount['number']} {amount['unit']}",
        {"reference_molecule_id": reference_mid, "amount": amount}, evidence)
    r["solvent_conditions"]["composition"] = f(solvent_mid, [solvent_mid], evidence)
    if volume is not None:
        r["solvent_conditions"]["total_volume"] = f(
            f"{volume['number']} {volume['unit']}", volume, evidence)


# Page 2: three cited design precedents and this-work generic concept.
for idx, desc in [(1, "Sheehan 1955 carbodiimide coupling concept"),
                  (2, "Stevens 1958 ketenimine coupling concept"),
                  (3, "Zhao 2016 ynamide coupling concept")]:
    reaction(f"scheme1_literature_eq{idx}", [EV[2]], f"Scheme 1, eq {idx}", desc,
             "cited peptide/amide coupling precedent", ["other"], None,
             notes=[f(desc, desc, [EV[2]])],
             unresolved=["Cited literature concept only; not a present-paper experiment."],
             record_type="literature_example")

r = reaction("scheme1_current_generic", [EV[1], EV[2]], "Graphical abstract and Scheme 1, eq 4",
             "Allenone-mediated amide bond formation", "generic allenone-mediated amide formation",
             ["other"], "C-N",
             participants=[participant("acid", GEN_ACID, "reactant", [EV[2]]),
                           participant("amine", GEN_AMINE, "reactant", [EV[2]]),
                           participant("allenone", M2A, "reagent", [EV[1], EV[2]]),
                           participant("amide", GEN_AMIDE, "product", [EV[1], EV[2]], product_role="target")],
             notes=[f("This work", "General design concept; no experiment-specific yield.", [EV[2]])],
             unresolved=["Generic R-group scheme; no experiment-specific quantities."],
             record_type="generic_scheme")


# Scheme 3 transcription. 3a/4a are also shown in Scheme 2 and use both evidence pages.
scheme3 = [
    ("3a", "8 h", 8, 95, "4a", "5 min", 5, "min", 99, 94, False, False),
    ("3b", "8 h", 8, 93, "4b", "10 min", 10, "min", 96, 94, False, False),
    ("3c", "12 h", 12, 92, "4c", "5 min", 5, "min", 92, 93, False, False),
    ("3d", "19 h", 19, 96, "4d", "5 min", 5, "min", 94, 95, False, False),
    ("3e", "4 h", 4, 93, "4e", "5 min", 5, "min", 97, 97, False, False),
    ("3f", "2 h", 2, 95, "4f", "5 min", 5, "min", 92, 97, False, False),
    ("3j", "12 h", 12, 96, "4g", "15 min", 15, "min", 95, 94, True, True),
    ("3h", "10 h", 10, 96, "4h", "6 h", 6, "h", 95, 94, True, True),
    ("3e", "4 h", 4, 93, "4i", "24 h", 24, "h", 99, 96, True, True),
]
seen_activation = set()
scheme3_raw = []
for ester, t_raw, t_num, y1, product, t2_raw, t2_num, t2_unit, y2, yop, hobt, onepot_hobt in scheme3:
    evi = [EV[3]] + ([EV[2]] if ester == "3a" or product == "4a" else [])
    acid = M1A if ester == "3a" else unresolved_molecule(f"acid_for_{ester}", f"carboxylic acid precursor drawn for {ester}", evi)
    ester_mid = M3A if ester == "3a" else unresolved_molecule(ester, f"alpha-carbonyl vinyl ester {ester}", evi)
    amine = MPHEN if product in ("4a", "4b", "4c", "4d", "4f") else unresolved_molecule(f"amine_for_{product}", f"amine precursor drawn for {product}", evi)
    product_mid = M4A if product == "4a" else unresolved_molecule(product, f"amide product {product}", evi)
    if ester not in seen_activation:
        seen_activation.add(ester)
        parts = [participant("acid", acid, "reactant", evi, amount=("0.22 mmol", .22, "mmol")),
                 participant("allenone", M2A, "reactant", evi, amount=("0.2 mmol", .2, "mmol")),
                 participant("DCE", MDCE, "solvent", evi, volume=("2 mL", 2, "mL")),
                 participant("ester", ester_mid, "product", evi, product_role="target")]
        stp = step("activation", [p["participant_id"] for p in parts], evi,
                   "1,4-addition/isomerization activation", "rt", {"number": t_num, "unit": "h"})
        rr = reaction(f"scheme3_activation_{ester}", evi, f"Scheme 3, {ester}",
                      "Allenone-mediated carboxylic acid activation", "alpha-carbonyl vinyl ester formation",
                      ["substrate_scope"], "C-O", parts, [stp],
                      [outcome("yield", "yield", ester_mid, f"{y1}%", {"number": y1, "unit": "%"}, evi,
                               basis="isolated", yield_type="isolated yield")],
                      unresolved=["Complete source structure and atom mapping remain unresolved."])
        set_common_solvent_and_scale(rr, M2A, {"number": .2, "unit": "mmol"}, MDCE, {"number": 2, "unit": "mL"}, evi)
    parts = [participant("ester", ester_mid, "reactant", evi, amount=("0.2 mmol", .2, "mmol")),
             participant("amine", amine, "reactant", evi, amount=("0.22 mmol", .22, "mmol")),
             participant("DMF", MDMF, "solvent", evi, volume=("1 mL", 1, "mL")),
             participant("product", product_mid, "product", evi, product_role="target")]
    if hobt:
        parts.append(participant("HOBt", MHOBT, "catalyst", evi, molpct=("10 mol%", 10), ref=ester_mid))
    rr = reaction(f"scheme3_aminolysis_{product}", evi, f"Scheme 3, {product}",
                  "Aminolysis of alpha-carbonyl vinyl ester", "amide bond formation",
                  ["substrate_scope"], "C-N", parts,
                  [step("aminolysis", [p["participant_id"] for p in parts], evi, "aminolysis", "rt", {"number": t2_num, "unit": t2_unit})],
                  [outcome("yield", "yield", product_mid, f"{y2}%", {"number": y2, "unit": "%"}, evi,
                           basis="isolated", yield_type="isolated yield")],
                  unresolved=["Complete source structure and atom mapping remain unresolved."])
    set_common_solvent_and_scale(rr, ester_mid, {"number": .2, "unit": "mmol"}, MDMF, {"number": 1, "unit": "mL"}, evi)
    parts = [participant("acid", acid, "reactant", evi, amount=("0.22 mmol", .22, "mmol")),
             participant("allenone", M2A, "reagent", evi, amount=("0.2 mmol", .2, "mmol")),
             participant("amine", amine, "reactant", evi, amount=("0.22 mmol", .22, "mmol")),
             participant("DCE", MDCE, "solvent", evi, volume=("2 mL", 2, "mL")),
             participant("DMF", MDMF, "solvent", evi, volume=("1 mL", 1, "mL")),
             participant("product", product_mid, "product", evi, product_role="target")]
    if onepot_hobt:
        parts.append(participant("HOBt", MHOBT, "catalyst", evi, molpct=("10 mol%", 10), ref=acid))
    rr = reaction(f"scheme3_onepot_{product}", evi, f"Scheme 3, {product}, footnote {'d' if onepot_hobt else 'c'}",
                  "Two-step one-pot allenone-mediated amide formation", "one-pot amide bond formation",
                  ["substrate_scope"], "C-N", parts,
                  [step("activation", ["acid", "allenone", "DCE"], evi, "activation in DCE", "rt"),
                   step("aminolysis", ["amine", "DMF", "product"] + (["HOBt"] if onepot_hobt else []), evi, "one-pot aminolysis in DMF", "rt")],
                  [outcome("yield", "yield", product_mid, f"{yop}%", {"number": yop, "unit": "%"}, evi,
                           basis="isolated", yield_type="isolated yield")],
                  unresolved=["One-pot stage times are not reported on pages 1-6."])
    set_common_solvent_and_scale(rr, M2A, {"number": .2, "unit": "mmol"}, MDCE, None, evi)
    scheme3_raw.append({"ester": ester, "activation_time": t_raw, "activation_yield_percent": y1,
                        "amide": product, "aminolysis_time": t2_raw, "aminolysis_yield_percent": y2,
                        "onepot_yield_percent": yop, "HOBt_aminolysis": hobt, "HOBt_onepot": onepot_hobt})

# Scheme 2 side-condensation control.
parts = [participant("diketone", M5, "reactant", [EV[2]]),
         participant("amine", MPHEN, "reactant", [EV[2]]),
         participant("DMF", MDMF, "solvent", [EV[2]]),
         participant("trace_product", M5P, "product", [EV[2]], product_role="trace product")]
reaction("scheme2_control_eq3", [EV[2]], "Scheme 2, eq 3", "Side-condensation control",
         "1,3-dicarbonyl/amine side condensation", ["control"], None, parts,
         [step("control", [p["participant_id"] for p in parts], [EV[2]], "mix diketone 5 and phenethylamine in DMF", "rt", {"number": 24, "unit": "h"})],
         [outcome("trace", "trace", M5P, "trace", "trace", [EV[2]])],
         unresolved=["Trace-product structure and measurement/detection basis are not fully reported on pages 1-6."])


# Scheme 4: 26 unique active-esters and 36 sequential + 36 one-pot peptide couplings.
active = {
"9a":(8,92,"ee"),"9b":(3,99,"ee"),"9c":(7,99,"de"),"9d":(8,95,"ee"),"9e":(6,90,"ee"),"9f":(8,99,"ee"),
"9g":(6,91,"ee"),"9h":(6,91,"ee"),"9i":(3,99,"ee"),"9j":(6,94,"de"),"9k":(12,90,"de"),"9l":(14,92,"de"),
"9m":(5,97,"ee"),"9n":(12,90,"ee"),"9o":(3,99,"ee"),"9p":(12,24,"ee"),"9q":(6,98,"ee"),"9r":(4,90,"ee"),
"9s":(10,91,"ee"),"9t":(6,95,"ee"),"9u":(10,95,"ee"),"9v":(10,92,"ee"),"9w":(12,94,"ee"),"9x":(12,92,"ee"),
"9y":(8,93,"ee"),"9z":(8,93,"ee")}
finals = [
("10a","9a",1.5,"h",97,95,False,False),("10b","9b",5,"min",94,93,False,False),("10c","9c",15,"min",98,99,True,True),
("10d","9d",10,"min",99,92,False,False),("10e","9e",20,"min",79,90,False,False),("10f","9f",5,"min",95,94,False,False),
("10g","9f",30,"min",78,77,True,True),("10h","9g",15,"min",93,90,True,True),("10i","9g",1.5,"h",92,90,True,True),
("10j","9a",2,"h",93,95,False,False),("10k","9h",30,"min",90,90,True,True),("10l","9h",1.5,"h",98,95,True,True),
("10m","9i",10,"min",96,92,False,False),("10n","9b",1.5,"h",94,95,True,True),("10o","9i",1.5,"h",94,95,True,True),
("10p","9j",25,"min",99,91,False,False),("10q","9k",1.5,"h",90,90,False,False),("10r","9l",25,"min",99,91,False,False),
("10s","9m",15,"min",99,97,False,False),("10t","9n",5,"min",99,92,False,False),("10u","9o",4.5,"h",91,90,False,False),
("10v","9b",2,"h",96,90,True,False),("10w","9p",20,"min",99,90,False,False),("10x","9a",2,"h",90,92,False,False),
("10y","9q",30,"min",95,96,False,False),("10z","9r",15,"min",92,93,False,False),("10aa","9s",25,"min",98,91,False,False),
("10ab","9t",15,"min",95,94,False,False),("10ac","9u",1.5,"h",99,92,False,False),("10ad","9v",10,"min",96,91,False,False),
("10ae","9w",1,"h",95,92,False,False),("10af","9x",1.5,"h",94,96,False,False),("10ag","9y",2,"h",95,94,False,False),
("10ah","9z",1,"h",95,94,False,False),("10ai","9b",15,"min",92,96,False,False),("10aj","9a",15,"min",93,96,False,False)]
scheme4_raw = {"active_esters": [], "couplings": []}
for ester, (hours, yld, metric) in active.items():
    acid = unresolved_molecule(f"aa_for_{ester}", f"amino-acid precursor drawn for {ester}", [EV[4]], stereo="drawn")
    emid = unresolved_molecule(ester, f"amino-acid alpha-carbonyl vinyl ester {ester}", [EV[4]], stereo="drawn")
    parts = [participant("amino_acid", acid, "reactant", [EV[4]], amount=("0.22 mmol", .22, "mmol")),
             participant("allenone", M2A, "reactant", [EV[4]], amount=("0.2 mmol", .2, "mmol")),
             participant("DCE", MDCE, "solvent", [EV[4]], volume=("2 mL", 2, "mL")),
             participant("active_ester", emid, "product", [EV[4]], product_role="target")]
    outs = [outcome("yield", "yield", emid, f"{yld}%", {"number": yld, "unit": "%"}, [EV[4]], basis="isolated", yield_type="isolated yield"),
            outcome(metric, "selectivity", emid, f"{metric} >99%", {"number": 99, "unit": "%", "qualifier": "gt", "raw_text": f"{metric} >99%"}, [EV[4]],
                    basis="HPLC analysis", scope={"metric": metric, "reported_relation": ">99%", "configuration_assignment": None})]
    rr = reaction(f"scheme4_activation_{ester}", [EV[4]], f"Scheme 4, {ester}",
                  "Allenone-mediated amino-acid activation", "alpha-carbonyl vinyl ester formation",
                  ["substrate_scope"], "C-O", parts,
                  [step("activation", [p["participant_id"] for p in parts], [EV[4]], "activation in DCE", "rt", {"number": hours, "unit": "h"})],
                  outs, unresolved=["Dense source structure is not encoded as a complete stereochemical SMILES."])
    set_common_solvent_and_scale(rr, M2A, {"number": .2, "unit": "mmol"}, MDCE, {"number": 2, "unit": "mL"}, [EV[4]])
    scheme4_raw["active_esters"].append({"label": ester, "time_h": hours, "yield_percent": yld, "selectivity": f"{metric} >99%"})

for product, ester, tnum, tunit, yseq, yone, hobt_seq, hobt_one in finals:
    emid = molecules[f"{DOC}:mol:{ester}"]["molecule_id"]
    acid = molecules[f"{DOC}:mol:aa_for_{ester}"]["molecule_id"]
    amine = unresolved_molecule(f"amine_for_{product}", f"amino coupling partner drawn for {product}", [EV[4]], stereo="drawn")
    pmid = unresolved_molecule(product, f"dipeptide product {product}", [EV[4]], stereo="drawn")
    parts = [participant("active_ester", emid, "reactant", [EV[4]], amount=("0.2 mmol", .2, "mmol")),
             participant("amine", amine, "reactant", [EV[4]], amount=("0.22 mmol", .22, "mmol")),
             participant("DMF", MDMF, "solvent", [EV[4]], volume=("1.5 mL", 1.5, "mL")),
             participant("product", pmid, "product", [EV[4]], product_role="target")]
    if hobt_seq:
        parts.append(participant("HOBt", MHOBT, "catalyst", [EV[4]], molpct=("10 mol%", 10), ref=emid))
    temperature = {"number": -55, "unit": "°C"} if product == "10u" else "rt"
    outs = [outcome("yield", "yield", pmid, f"{yseq}%", {"number": yseq, "unit": "%"}, [EV[4]], basis="isolated", yield_type="isolated yield"),
            outcome("de", "selectivity", pmid, "de >99%", {"number": 99, "unit": "%", "qualifier": "gt", "raw_text": "de >99%"}, [EV[4]],
                    basis="HPLC analysis", scope={"metric": "de", "diastereomer_order": None, "configuration_assignment": None})]
    rr = reaction(f"scheme4_aminolysis_{product}", [EV[4]], f"Scheme 4, {product}",
                  "Aminolysis of amino-acid alpha-carbonyl vinyl ester", "peptide bond formation",
                  ["substrate_scope"], "C-N", parts,
                  [step("aminolysis", [p["participant_id"] for p in parts], [EV[4]], "aminolysis in DMF", temperature, {"number": tnum, "unit": tunit})],
                  outs, unresolved=["Absolute configuration descriptors and complete source structures are not encoded."])
    set_common_solvent_and_scale(rr, emid, {"number": .2, "unit": "mmol"}, MDMF, {"number": 1.5, "unit": "mL"}, [EV[4]])
    parts = [participant("amino_acid", acid, "reactant", [EV[4]], amount=("0.22 mmol", .22, "mmol")),
             participant("allenone", M2A, "reagent", [EV[4]], amount=("0.2 mmol", .2, "mmol")),
             participant("amine", amine, "reactant", [EV[4]], amount=("0.22 mmol", .22, "mmol")),
             participant("DCE", MDCE, "solvent", [EV[4]], volume=("2 mL", 2, "mL")),
             participant("DMF", MDMF, "solvent", [EV[4]], volume=("1.5 mL", 1.5, "mL")),
             participant("product", pmid, "product", [EV[4]], product_role="target")]
    if hobt_one:
        parts.append(participant("HOBt", MHOBT, "catalyst", [EV[4]], molpct=("10 mol%", 10), ref=acid))
    rr = reaction(f"scheme4_onepot_{product}", [EV[4]], f"Scheme 4, {product}, footnote {'d' if hobt_one else 'c'}",
                  "Two-step one-pot allenone-mediated peptide bond formation", "one-pot peptide bond formation",
                  ["substrate_scope"], "C-N", parts,
                  [step("activation", ["amino_acid", "allenone", "DCE"], [EV[4]], "activation in DCE", "rt"),
                   step("aminolysis", ["amine", "DMF", "product"] + (["HOBt"] if hobt_one else []), [EV[4]], "one-pot aminolysis in DMF", "rt")],
                  [outcome("yield", "yield", pmid, f"{yone}%", {"number": yone, "unit": "%"}, [EV[4]], basis="isolated", yield_type="isolated yield"),
                   outcome("de", "selectivity", pmid, "de >99%", {"number": 99, "unit": "%", "qualifier": "gt", "raw_text": "de >99%"}, [EV[4]],
                           basis="HPLC analysis", scope={"metric": "de", "diastereomer_order": None, "configuration_assignment": None})],
                  unresolved=["One-pot stage times are not reported on page 4; dense source structures are not fully encoded."])
    set_common_solvent_and_scale(rr, M2A, {"number": .2, "unit": "mmol"}, MDCE, None, [EV[4]])
    scheme4_raw["couplings"].append({"product": product, "active_ester": ester, "time": f"{tnum} {tunit}",
                                      "sequential_yield_percent": yseq, "onepot_yield_percent": yone,
                                      "sequential_HOBt": hobt_seq, "onepot_HOBt": hobt_one,
                                      "temperature": "-55 °C" if product == "10u" else "rt", "de": ">99%"})


# Scheme 5 peptide fragment condensation.
scheme5 = [
("13a",10,99,"14a",20,"min",90,"[2+1]",False),("13b",18,95,"14b",30,"min",93,"[2+1]",False),
("13c",36,92,"14c",1,"h",95,"[2+1]",False),("13d",10,95,"14d",20,"min",93,"[2+1]",False),
("13e",12,99,"14e",20,"min",95,"[2+2]",False),("13f",96,91,"14f",20,"min",90,"[2+3]",False),
("13g",48,96,"14g",2,"h",90,"[3+1]",False),("13h",48,99,"14h",1.5,"h",94,"[3+1]",False),
("13i",36,99,"14i",30,"min",95,"[3+1]",True),("13j",48,92,"14j",40,"min",93,"[4+2]",False),
("13k",48,98,"14k",45,"min",92,"[4+2]",False)]
scheme5_raw = []
for ester, th, y1, prod, tn, tu, y2, span, zero in scheme5:
    acid = unresolved_molecule(f"peptide_acid_for_{ester}", f"peptide acid 11 precursor drawn for {ester}", [EV[5]], stereo="drawn")
    emid = unresolved_molecule(ester, f"peptide alpha-carbonyl vinyl ester {ester}", [EV[5]], stereo="drawn")
    amine = unresolved_molecule(f"peptide_amine_for_{prod}", f"amine/peptide fragment 12 drawn for {prod}", [EV[5]], stereo="drawn")
    pmid = unresolved_molecule(prod, f"peptide product {prod}", [EV[5]], stereo="drawn")
    parts = [participant("peptide_acid", acid, "reactant", [EV[5]], amount=("0.1 mmol", .1, "mmol")),
             participant("allenone", M2A, "reactant", [EV[5]], amount=("0.2 mmol", .2, "mmol")),
             participant("DCE", MDCE, "solvent", [EV[5]], volume=("2 mL", 2, "mL")),
             participant("active_ester", emid, "product", [EV[5]], product_role="target")]
    rr = reaction(f"scheme5_activation_{ester}", [EV[5]], f"Scheme 5, {ester}",
                  "Allenone-mediated peptide-acid activation", "peptide alpha-carbonyl vinyl ester formation",
                  ["substrate_scope"], "C-O", parts,
                  [step("activation", [p["participant_id"] for p in parts], [EV[5]], "activation in DCE", "rt", {"number": th, "unit": "h"})],
                  [outcome("yield", "yield", emid, f"{y1}%", {"number": y1, "unit": "%"}, [EV[5]], basis="isolated", yield_type="isolated yield"),
                   outcome("de", "selectivity", emid, "de >99%", {"number": 99, "unit": "%", "qualifier": "gt", "raw_text": "de >99%"}, [EV[5]],
                           basis="HPLC analysis", scope={"metric": "de", "configuration_assignment": None})],
                  unresolved=["Dense peptide structure is not encoded as a complete stereochemical SMILES."])
    set_common_solvent_and_scale(rr, acid, {"number": .1, "unit": "mmol"}, MDCE, {"number": 2, "unit": "mL"}, [EV[5]])
    parts = [participant("active_ester", emid, "reactant", [EV[5]], amount=("0.1 mmol", .1, "mmol")),
             participant("amine", amine, "reactant", [EV[5]], amount=("0.11 mmol", .11, "mmol")),
             participant("DMF", MDMF, "solvent", [EV[5]], volume=("1 mL", 1, "mL")),
             participant("product", pmid, "product", [EV[5]], product_role="target")]
    temp = {"number": 0, "unit": "°C"} if zero else "rt"
    rr = reaction(f"scheme5_condensation_{prod}", [EV[5]], f"Scheme 5, {prod}, fragment span {span}",
                  "Allenone-mediated peptide fragment condensation", "peptide bond formation",
                  ["substrate_scope"], "C-N", parts,
                  [step("aminolysis", [p["participant_id"] for p in parts], [EV[5]], "fragment aminolysis in DMF", temp, {"number": tn, "unit": tu})],
                  [outcome("yield", "yield", pmid, f"{y2}%", {"number": y2, "unit": "%"}, [EV[5]], basis="isolated", yield_type="isolated yield"),
                   outcome("de", "selectivity", pmid, "de >99%", {"number": 99, "unit": "%", "qualifier": "gt", "raw_text": "de >99%"}, [EV[5]],
                           basis="HPLC analysis", scope={"metric": "de", "configuration_assignment": None, "fragment_span": span})],
                  unresolved=["Dense peptide structure and complete atom mapping remain unresolved."])
    set_common_solvent_and_scale(rr, emid, {"number": .1, "unit": "mmol"}, MDMF, {"number": 1, "unit": "mL"}, [EV[5]])
    scheme5_raw.append({"active_ester": ester, "activation_time_h": th, "activation_yield_percent": y1,
                        "product": prod, "fragment_span": span, "condensation_time": f"{tn} {tu}",
                        "condensation_temperature": "0 °C" if zero else "rt", "yield_percent": y2, "de": ">99%"})


# Scheme 6 carfilzomib route: nine displayed arrows.
for label in ["chloroacetic_acid", "15", "Phe_OtBu", "16", "17", "Leu_OtBu", "18", "19", "20", "21", "warhead_amine", "22", "23"]:
    if label == "chloroacetic_acid":
        simple_molecule(label, "chloroacetic acid", "O=C(O)CCl", [EV[5]])
    else:
        unresolved_molecule(label, "carfilzomib" if label == "23" else f"Scheme 6 material {label}", [EV[5]], stereo="drawn" if label not in ("15",) else "unresolved")

route = [
    ("15", "chloroacetic_acid", M2A, "activation", "C-O", 10, "h", 94, "DCE", "rt"),
    ("16", "15", "Phe_OtBu", "coupling", "C-N", 15, "min", 95, "DMF", "rt"),
    ("17", "16", M2A, "deprotection_activation", "C-O", None, None, 93, "TFA/DCM then 2a/DCE", None),
    ("18", "17", "Leu_OtBu", "coupling", "C-N", 20, "min", 99, "DMF", "rt"),
    ("19", "18", M2A, "deprotection_activation", "C-O", None, None, 98, "TFA/DCM then 2a/DCE", None),
    ("20", "19", "Phe_OtBu", "coupling", "C-N", 20, "min", 92, "DMF", "rt"),
    ("21", "20", M2A, "deprotection_activation", "C-O", None, None, 99, "TFA/DCM then 2a/DCE", None),
    ("22", "21", "warhead_amine", "coupling", "C-N", 20, "min", 93, "DMF", "rt"),
    ("23", "22", MMOR, "substitution", "C-N", 12, "h", 99, "THF/KI", "rt"),
]
scheme6_raw = []
for prod, start, partner, kind, clabel, tn, tu, yld, medium, temp in route:
    smid = molecules[f"{DOC}:mol:{start}"]["molecule_id"] if not start.startswith(DOC) else start
    pmid = molecules[f"{DOC}:mol:{prod}"]["molecule_id"]
    partner_mid = partner if isinstance(partner, str) and partner.startswith(DOC) else molecules[f"{DOC}:mol:{partner}"]["molecule_id"]
    parts = [participant("start", smid, "reactant", [EV[5]]),
             participant("partner", partner_mid, "reagent" if kind == "deprotection_activation" else "reactant", [EV[5]]),
             participant("product", pmid, "product", [EV[5]], product_role="target")]
    if kind == "activation":
        parts += [participant("DCE", MDCE, "solvent", [EV[5]])]
    elif kind == "coupling":
        parts += [participant("DMF", MDMF, "solvent", [EV[5]])]
    elif kind == "deprotection_activation":
        parts += [participant("TFA", MTFA, "reagent", [EV[5]]), participant("DCM", MDCM, "solvent", [EV[5]]), participant("DCE", MDCE, "solvent", [EV[5]])]
    else:
        parts += [participant("KI", MKI, "reagent", [EV[5]]), participant("THF", MTHF, "solvent", [EV[5]])]
    steps_here = []
    if kind == "deprotection_activation":
        steps_here = [step("deprotection", ["start", "TFA", "DCM"], [EV[5]], "TFA deprotection in DCM"),
                      step("activation", ["partner", "DCE", "product"], [EV[5]], "allenone activation in DCE")]
    else:
        steps_here = [step(kind, [p["participant_id"] for p in parts], [EV[5]], kind.replace("_", " "), temp,
                           {"number": tn, "unit": tu} if tn is not None else None)]
    rr = reaction(f"scheme6_{start}_to_{prod}", [EV[5]], f"Scheme 6, {start} to {prod}",
                  "Carfilzomib synthesis with allenone 2a", kind.replace("_", " "),
                  ["other"], clabel, parts, steps_here,
                  [outcome("yield", "yield", pmid, f"{yld}%", {"number": yld, "unit": "%"}, [EV[5]])],
                  notes=[f("two steps" if kind == "deprotection_activation" else medium,
                           {"reported_medium": medium, "yield_scope": "two-step" if kind == "deprotection_activation" else "single arrow"}, [EV[5]])],
                  unresolved=["Yield measurement basis is not stated on page 5."] + (["Individual deprotection and activation times and per-step yields are not reported; the displayed yield spans both steps."] if kind == "deprotection_activation" else []))
    scheme6_raw.append({"start": start, "product": prod, "operation": kind, "medium": medium,
                        "time": None if tn is None else f"{tn} {tu}", "temperature": temp,
                        "yield_percent": yld, "yield_scope": "two-step" if kind == "deprotection_activation" else "single arrow"})


# Page 6: one multicycle SPPS experiment plus three comparator syntheses.
MRESIN = unresolved_molecule("Fmoc_Gly_CTC_resin", "Fmoc-Gly-CTC resin (0.77 mmol/g loading)", [EV[6]])
MBLOCKS = generic_molecule("aa_vinyl_ester_building_blocks", "proteinogenic amino-acid alpha-carbonyl vinyl ester building blocks", "Fmoc-AA-O-C(Me)=CH-C(=O)-Ph", ["AA"], [EV[6]])
MACP = unresolved_molecule("ACP_65_74", "H-Val-Gln-Ala-Ala-Ile-Asp-Tyr-Ile-Asn-Gly-OH (ACP 65-74)", [EV[6]], stereo="drawn")
for reagent in ["PyBOP", "HBTU", "pentafluorophenyl_ester_building_blocks"]:
    unresolved_molecule(reagent, reagent.replace("_", " "), [EV[6]])

parts = [participant("resin", MRESIN, "reactant", [EV[6]]),
         participant("building_blocks", MBLOCKS, "reactant", [EV[6]], equiv=("2 equiv per coupling", 2), ref=MRESIN),
         participant("HOBt", MHOBT, "catalyst", [EV[6]], equiv=("0.1 equiv", .1), ref=MRESIN),
         participant("DMF", MDMF, "solvent", [EV[6]]),
         participant("piperidine", MPIP, "reagent", [EV[6]]),
         participant("TFA", MTFA, "workup_agent", [EV[6]], stage="workup"),
         participant("TIS", MTIS, "workup_agent", [EV[6]], stage="workup"),
         participant("product", MACP, "product", [EV[6]], product_role="target")]
steps = []
residues = ["Asn", "Ile", "Tyr", "Asp", "Ile", "Ala", "Ala", "Gln", "Val"]
for i, aa in enumerate(residues, 1):
    trange = {"minimum": 1, "maximum": 2, "unit": "h", "raw_text": "1-2 h"} if aa in ("Val", "Ile") else {"minimum": 20, "maximum": 30, "unit": "min", "raw_text": "20-30 min"}
    steps.append(step(f"coupling_{i}_{aa}", ["resin", "building_blocks", "HOBt", "DMF"], [EV[6]], f"couple {aa} alpha-carbonyl vinyl ester", None, trange))
    steps.append(step(f"deprotection_{i}_{aa}", ["resin", "piperidine", "DMF"], [EV[6]], "Fmoc removal with piperidine/DMF (1/4)"))
steps.append(step("cleavage", ["resin", "TFA", "TIS", "product"], [EV[6]], "TFA/TIS/H2O cleavage (95/2.5/2.5)", stage="workup"))
rr = reaction("scheme7_ACP_vinyl_ester_SPPS", [EV[6], TXT[6]], "Scheme 7 and Figure 1, vinyl ester method",
              "Solid-phase synthesis of ACP (65-74)", "iterative peptide bond formation",
              ["other"], "C-N", parts, steps,
              [outcome("purity", "other", MACP, "purity 98%", {"number": 98, "unit": "%"}, [EV[6]], basis="HPLC profile of crude ACP (65-74)", scope={"metric": "crude purity"})],
              notes=[f("Fmoc-Gly-CTC resin (0.77 mmol/g resin loading)", {"resin_loading": {"number": .77, "unit": "mmol/g"}}, [EV[6]]),
                     f("coupling 2 equiv vinyl ester, 0.1 equiv HOBt, DMF; Fmoc removal piperidine/DMF (1/4); cleavage TFA/TIS/H2O=95/2.5/2.5",
                       "Common SPPS conditions", [EV[6]])],
              unresolved=["No isolated yield is reported on page 6; 98% is crude HPLC purity, not product yield.",
                          "Per-cycle conversions and exact times for each non-hindered residue are not individually tabulated."])
rr["classification"] = classification("C-N", [EV[6]], count=9,
                                      limitation="Nine iterative peptide-bond formations are represented at residue-sequence level; full atom mapping is unresolved.")

for method_label, reagent_label, purity in [
    ("PyBOP", "PyBOP", 86), ("HBTU", "HBTU", 77),
    ("aryl_ester", "pentafluorophenyl_ester_building_blocks", 60)]:
    reagent_mid = molecules[f"{DOC}:mol:{reagent_label}"]["molecule_id"]
    parts = [participant("resin", MRESIN, "reactant", [EV[6]]),
             participant("coupling_material", reagent_mid, "reagent", [EV[6]]),
             participant("product", MACP, "product", [EV[6]], product_role="target")]
    reaction(f"figure1_ACP_control_{method_label}", [EV[6], TXT[6]], f"Figure 1, {method_label} comparison",
             f"ACP (65-74) SPPS comparison using {method_label}", "comparative solid-phase peptide synthesis",
             ["control"], "C-N", parts,
             [step("SPPS", [p["participant_id"] for p in parts], [EV[6]], f"SPPS using {method_label}")],
             [outcome("purity", "other", MACP, f"purity {purity}%", {"number": purity, "unit": "%"}, [EV[6]],
                      basis="HPLC profile of crude ACP (65-74)", scope={"metric": "crude purity"})],
             unresolved=["Pages 1-6 do not report the complete comparator procedure, quantities, temperatures, or cycle times.",
                         "Reported percentage is crude HPLC purity, not yield."])


bundle["molecules"] = list(molecules.values())
bundle["source_context"] = [
    {"context_id": f"{DOC}:graphical_abstract", "document_id": DOC, "evidence_ids": [EV[1]],
     "kind": "graphical_abstract", "note": "Generic peptide-acid/amine coupling through an alpha-carbonyl vinyl ester; duplicate conceptual presentation, not a separate experiment."},
    {"context_id": f"{DOC}:figure1_hplc", "document_id": DOC, "evidence_ids": [EV[6]],
     "kind": "characterization_comparison", "note": "HPLC profiles of crude ACP (65-74): PyBOP 86%, HBTU 77%, aryl ester 60%, vinyl ester 98%. These are purity measurements, not synthetic yields."},
    {"context_id": f"{DOC}:storage_discussion", "document_id": DOC, "evidence_ids": [TXT[6]],
     "kind": "stability_statement", "note": "Authors state allenone/active esters can be refrigerated for 12 months and stock solutions remain intact in DMF at room temperature for 5 days; these are prose observations, not separate synthesis records."}
]

dump(RUN / "submitted.json", bundle)
diagnostics = []
for molecule in bundle["molecules"]:
    smiles = molecule["smiles"]["normalized_value"]
    if smiles:
        diagnostics.append({"molecule_id": molecule["molecule_id"], "smiles": smiles,
                            "diagnostic": st.stereo_inspect(smiles)})
dump(RUN / "raw_model/stereo_diagnostics.json", diagnostics)
dump(RUN / "raw_model/visual_transcription.json", {
    "document_id": DOC, "pages_visually_inspected": [1, 2, 3, 4, 5, 6],
    "scheme3": scheme3_raw, "scheme4": scheme4_raw, "scheme5": scheme5_raw,
    "scheme6": scheme6_raw,
    "limitations": ["Machine transcription remains needs_review.", "Dense peptide structures were not converted to complete SMILES.", "No SI was supplied or associated."]
})

page_records = {p: [] for p in range(1, 7)}
for r in bundle["reactions"]:
    for eid in r["evidence_ids"]:
        for p in range(1, 7):
            if eid in (EV[p], TXT[p]):
                page_records[p].append(r["reaction_id"])
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title/authors/abstract", "finding": "True article title, authors, DOI context, and method summary; no experiment-specific quantities.", "evidence_ids": [EV[1], TXT[1]]},
            {"region": "graphical abstract", "finding": "Generic allenone-mediated peptide coupling through an alpha-carbonyl vinyl ester; duplicate concept, not a separate experiment.", "evidence_ids": [EV[1]]}],
         "record_ids": sorted(set(page_records[1])), "unresolved": ["Graphical abstract peptide identities and quantities are generic/unreported."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Scheme 1", "finding": "Three cited literature coupling concepts and one present-work generic design; cited work is not promoted to a current-paper experiment.", "evidence_ids": [EV[2], TXT[2]]},
            {"region": "Scheme 2 eqs 1-2", "finding": "Specific 1a to 3a activation and 3a to 4a aminolysis, linked to repeated Scheme 3 evidence rather than duplicated.", "evidence_ids": [EV[2], TXT[2]]},
            {"region": "Scheme 2 eq 3", "finding": "24 h room-temperature side-condensation control gives trace 5'; retained as a control.", "evidence_ids": [EV[2], TXT[2]]}],
         "record_ids": sorted(set(page_records[2])), "unresolved": ["Scheme 2 yield measurement basis is supplied only by the repeated Scheme 3 footnote for 3a/4a; 5' trace detection basis is unreported."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Scheme 3 general arrow and footnotes", "finding": "Common activation/aminolysis conditions, isolated-yield basis, HOBt exceptions, and one-pot variants transcribed.", "evidence_ids": [EV[3], TXT[3]]},
            {"region": "Scheme 3 product grid", "finding": "Eight unique active-ester preparations, nine aminolyses, and nine one-pot experiments; repeated 3e activation is not duplicated.", "evidence_ids": [EV[3]]},
            {"region": "right-column discussion", "finding": "Scheme 4/5 stereochemical and scope discussion; no additional standalone experiment beyond linked schemes.", "evidence_ids": [EV[3], TXT[3]]}],
         "record_ids": sorted(set(page_records[3])), "unresolved": ["Dense structures other than safely transcribed 3a/4a are retained without complete SMILES."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "Scheme 4 general arrow and footnotes", "finding": "Activation and aminolysis quantities, DCE/DMF volumes, isolated yields, HPLC basis, HOBt flags, one-pot variants, and -55 °C exception transcribed.", "evidence_ids": [EV[4], TXT[4]]},
            {"region": "Scheme 4 six-by-six product grid", "finding": "26 unique active-esters plus 36 sequential aminolyses and 36 distinct one-pot experiments; repeated active-esters are linked rather than duplicated.", "evidence_ids": [EV[4]]},
            {"region": "bottom prose", "finding": "Carfilzomib overall nine-step/68% statement begins; the actual route arrows are on page 5 and recorded there.", "evidence_ids": [EV[4], TXT[4]]}],
         "record_ids": sorted(set(page_records[4])), "unresolved": ["Full stereochemical SMILES for dense peptide structures are not safely transcribed; absolute descriptors are not printed.", "One-pot stage times are not reported."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "Scheme 5", "finding": "Eleven peptide-acid activations and eleven fragment aminolyses, including fragment spans, times, yields, HPLC de values, and the 0 °C exception.", "evidence_ids": [EV[5], TXT[5]]},
            {"region": "Scheme 6", "finding": "Nine displayed arrows in the carfilzomib route; three deprotection/activation arrows retain two-step yield scope.", "evidence_ids": [EV[5], TXT[5]]},
            {"region": "SPPS discussion opening", "finding": "Historical active-ester discussion only; no additional present-paper experiment beyond Scheme 7 on page 6.", "evidence_ids": [EV[5], TXT[5]]}],
         "record_ids": sorted(set(page_records[5])), "unresolved": ["Scheme 6 yield measurement basis and quantities for most arrows are unreported on supplied pages.", "Complete peptide structures and atom maps remain unresolved."]},
        {"page": 6, "state": "attempted", "regions": [
            {"region": "Scheme 7", "finding": "Nine iterative vinyl-ester couplings from Fmoc-Gly-CTC resin, repeated Fmoc deprotection, and final cleavage recorded as one multistep SPPS experiment.", "evidence_ids": [EV[6], TXT[6]]},
            {"region": "Scheme 7 footnote/prose", "finding": "2 equiv vinyl ester, 0.1 equiv HOBt, DMF, piperidine/DMF 1/4, TFA/TIS/H2O 95/2.5/2.5, coupling-time ranges, and resin loading retained.", "evidence_ids": [EV[6], TXT[6]]},
            {"region": "Figure 1 HPLC profiles", "finding": "Crude ACP purities 86/77/60/98% for PyBOP/HBTU/aryl ester/vinyl ester; retained as purity outcomes and characterization context, not yields.", "evidence_ids": [EV[6], TXT[6]]}],
         "record_ids": sorted(set(page_records[6])), "unresolved": ["Comparator procedures are not reported on pages 1-6.", "No isolated ACP yield is reported."]}
    ]
}
dump(RUN / "coverage/pages_0001_0006.json", coverage)

print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
                  "source_context": len(bundle["source_context"])}, ensure_ascii=False))
