#!/usr/bin/env python3
"""Source-derived v3.1 extraction for JACS 2021, 143, 12622-12632, PDF pages 1-6.

The values below are a visual transcription of the supplied page images.  The
script intentionally leaves dense structures as unresolved identities rather
than manufacturing SMILES from the drawings.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_f2a3e97f1c7189a7"
DATASET = ROOT / "runs/current/dataset.json"
SUBMITTED = ROOT / "runs/current/submitted.json"


def ev(page: int) -> str:
    return f"{DOC}:page{page}:image"


def F(raw=None, normalized=None, evidence=(), method="original_extraction",
      status=None, derivation=None):
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   derivation=derivation)


def q(number, unit, qualifier=None, raw_text=None):
    value = {"number": number, "unit": unit}
    if qualifier:
        value["qualifier"] = qualifier
    if raw_text:
        value["raw_text"] = raw_text
    return value


def safe(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_")


dataset = st.read(DATASET)
doc = dataset["documents"][0]
doc["title"] = ("Enantio- and Diastereodivergent Construction of 1,3-Nonadjacent "
                "Stereocenters Bearing Axial and Central Chirality through Synergistic "
                "Pd/Cu Catalysis")
doc["supplementary_material"] = {
    "status": "mentioned_not_supplied",
    "note": ("The main article explicitly mentions Supporting Information on PDF page 8, "
             "but no SI file was supplied or explicitly associated in source_mapping.json.")
}
st.write(DATASET, dataset)
bundle = st.empty(dataset["documents"])
molecules = {}


def molecule(label, evidence, key=None, name=None, stereo=None, note=None):
    mid = f"{DOC}:mol:{key or safe(label)}"
    evidence = list(dict.fromkeys(evidence))
    if mid in molecules:
        m = molecules[mid]
        occ = m["occurrences"][0]
        occ["evidence_ids"] = list(dict.fromkeys(occ["evidence_ids"] + evidence))
        occ["compound_label"]["evidence_ids"] = list(dict.fromkeys(
            occ["compound_label"]["evidence_ids"] + evidence))
        return mid
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, evidence)
    m["names"]["original"] = F(label, label, evidence)
    m["names"]["english"] = F(name or label, name or label, evidence)
    m["representation_type"] = "unresolved"
    m["smiles"] = F(evidence=evidence)
    m["generic_template"] = None
    m["cas_numbers"] = F(evidence=evidence)
    for field in ("atom_centers", "double_bonds", "axial", "mixture", "assignment",
                  "isomeric_smiles", "identity_key"):
        m["stereochemistry"][field] = F(evidence=evidence)
    m["stereochemistry"]["geometry_status"] = F(
        "geometry not independently resolved from the dense drawing", "unresolved", evidence)
    if stereo and stereo.get("kind") == "racemic":
        m["stereochemistry"]["mixture"] = F(
            stereo["raw"], {"type": "racemic", "axis": "allene"}, evidence)
        m["stereochemistry"]["axial"] = F(
            stereo["raw"], {"axis": "allene", "configuration": "racemic/paired enantiomers"}, evidence)
        m["stereochemistry"]["chirality_status"] = F(stereo["raw"], "racemic", evidence)
        m["stereochemistry"]["representation_status"] = F(
            "source axial chirality not encoded as SMILES", "not_encoded", evidence)
    elif stereo and stereo.get("kind") == "assigned":
        assignment = stereo["assignment"]
        m["stereochemistry"]["assignment"] = F(assignment, assignment, evidence)
        m["stereochemistry"]["chirality_status"] = F(assignment, "specified", evidence)
        m["stereochemistry"]["representation_status"] = F(
            "source configuration retained textually; no complete isomeric SMILES", "not_encoded", evidence)
        if stereo.get("axial"):
            m["stereochemistry"]["axial"] = F(
                stereo["axial"], {"axis": "allene", "configuration": stereo["axial"]}, evidence)
        if stereo.get("centers"):
            m["stereochemistry"]["atom_centers"] = F(
                stereo["centers"], stereo["centers"], evidence)
    else:
        m["stereochemistry"]["chirality_status"] = F(
            "stereochemistry not safely resolved", "unresolved", evidence)
        m["stereochemistry"]["representation_status"] = F(
            "complete structure not encoded", "unresolved", evidence)
    m["identity_note"] = F(
        note or "Source label/drawing retained; complete connectivity and stereochemical SMILES were not safely resolved from the main-text image.",
        note or "Source label/drawing retained; complete connectivity and stereochemical SMILES were not safely resolved from the main-text image.",
        evidence)
    m["external_properties"] = []
    molecules[mid] = m
    return mid


def assigned(label, evidence, key=None):
    config = label.split(")", 1)[0].lstrip("(") if label.startswith("(") else None
    parts = config.split(",") if config else []
    axial = None
    centers = None
    if len(parts) == 2:
        axial, centers = parts[0], [{"source_order": 2, "configuration": parts[1]}]
    elif len(parts) == 3:
        # Scheme 4 color coding: pre-existing alcohol center, allene axis, newly formed center.
        centers = [{"source_order": 1, "configuration": parts[0], "role": "pre-existing alcohol center"},
                   {"source_order": 3, "configuration": parts[2], "role": "new central stereocenter"}]
        axial = parts[1]
    return molecule(label, evidence, key=key,
                    stereo={"kind": "assigned", "assignment": config,
                            "axial": axial, "centers": centers})


def rac(label, evidence, key=None):
    return molecule(label, evidence, key=key,
                    stereo={"kind": "racemic", "raw": f"rac-{label.removeprefix('rac-')}"})


# Shared source-labelled materials.
M = {}
M["1a"] = rac("rac-1a", [ev(1), ev(2), ev(3), ev(5), ev(6)], "rac_1a")
for letter in "bcdefghijklmnop":
    M[f"1{letter}"] = rac(f"rac-1{letter}", [ev(4)], f"rac_1{letter}")
M["1q"] = assigned("(R)-1q", [ev(5)], "R_1q")
M["1qprime"] = assigned("(S)-1q′", [ev(5)], "S_1q_prime")
for label in ["2a′", "2a", "2a″", "2a‴"]:
    M[label] = molecule(label, [ev(2), ev(3), ev(5), ev(6)] if label == "2a" else [ev(2)])
for letter in "bcdefghijklmnopqrs":
    M[f"2{letter}"] = molecule(f"2{letter}", [ev(5)])

common_materials = {
    "Pd_allyl": ("[Pd(η3-C3H5)Cl]2", "allylpalladium chloride dimer"),
    "Cu_salt": ("[Cu(CH3CN)4]PF6", "tetrakis(acetonitrile)copper(I) hexafluorophosphate"),
    "base": ("Cs2CO3", "cesium carbonate"),
    "THF": ("THF", "tetrahydrofuran"),
    "mesitylene": ("mesitylene", "mesitylene"),
}
for key, (label, name) in common_materials.items():
    M[key] = molecule(label, [ev(2), ev(5), ev(6)], key=key, name=name)


def set_concrete(mid, smiles, evidence, *, partial=False, preserve_chirality=False):
    """Attach only visually/source-identifiable candidate connectivity and graph diagnostics."""
    m = molecules[mid]
    diag = st.stereo_inspect(smiles)
    derivation = ("Candidate connectivity transcribed from the source label/drawing and checked by RDKit; "
                  "parsing/graph diagnostics do not establish source identity correctness.")
    m["representation_type"] = "smiles"
    m["smiles"] = F("candidate connectivity from source label/drawing", smiles, evidence,
                    method="structure_derivation", derivation=derivation)
    m["stereochemistry"]["isomeric_smiles"] = F(
        None, diag["canonical_isomeric_smiles"], evidence, method="structure_derivation", derivation=derivation)
    m["stereochemistry"]["identity_key"] = F(
        None, diag["stereo_graph_key"], evidence, method="structure_derivation", derivation=derivation)
    m["stereochemistry"]["representation_status"] = F(
        "axial configuration is not encoded" if partial else "graph stereochemistry fully represented for this candidate",
        "partial" if partial else "encoded", evidence,
        method="structure_derivation", derivation=derivation)
    if not preserve_chirality:
        m["stereochemistry"]["chirality_status"] = F(
            "no potential stereocenter in candidate graph", "achiral", evidence,
            method="structure_derivation", derivation=derivation)
        m["stereochemistry"]["geometry_status"] = F(
            "no applicable E/Z assignment", "not_applicable", evidence,
            method="structure_derivation", derivation=derivation)
    return diag


stereo_diagnostics = {
    "rac-1a_connectivity_candidate": set_concrete(M["1a"], "CC(C)(C)C=C=COC(C)=O", [ev(1), ev(2), ev(6)], partial=True, preserve_chirality=True),
    "THF": set_concrete(M["THF"], "C1CCOC1", [ev(2), ev(4), ev(5), ev(6)]),
    "mesitylene": set_concrete(M["mesitylene"], "Cc1cc(C)cc(C)c1", [ev(5), ev(6)]),
    "cesium_carbonate": set_concrete(M["base"], "O=C([O-])[O-].[Cs+].[Cs+]", [ev(2), ev(5), ev(6)]),
}
ligands = {
    "R_L1": "(R)-BINAP (L1)", "R_L2": "(R)-MeO-BIPHEP (L2)",
    "R_L3": "(R)-SegPhos (L3)", "S_L3": "(S)-SegPhos (L3)",
    "R_L4": "(R)-DTBM-SegPhos (L4)", "S_L4": "(S)-DTBM-SegPhos (L4)",
    "S_DM_L3": "(S)-DM-SegPhos", "RRp_L5": "(R,Rp)-L5",
    "RRp_L6": "(R,Rp)-L6", "SSp_L6": "(S,Sp)-L6",
}
for key, label in ligands.items():
    M[key] = molecule(label, [ev(2), ev(3), ev(4), ev(5), ev(6)], key=key,
                      stereo={"kind": "assigned", "assignment": label.split("-")[0].strip("()"),
                              "centers": [{"source_assignment": label}]})


def product(label, pages, key=None):
    return assigned(label, [ev(p) for p in pages], key=key or safe(label))


def addp(r, mid, role, evidence, pid=None, amount=None, equiv=None, molpct=None,
         ref=None, stage="reaction", product_role=None):
    pid = pid or f"p{len(r['participants'])+1}_{safe(mid.split(':')[-1])}"
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    if product_role:
        p["product_role"] = F(product_role, product_role, evidence)
    if amount is not None:
        p["charges"]["amount"] = F(f"{amount} mmol", q(amount, "mmol"), evidence)
    if equiv is not None:
        p["charges"]["equivalents"] = F(f"{equiv} equiv", q(equiv, "equiv"), evidence)
    if molpct is not None:
        p["charges"]["mol_percent"] = F(f"{molpct} mol %", q(molpct, "mol%"), evidence)
    if ref is not None and (equiv is not None or molpct is not None):
        p["equivalence_reference"] = F(ref, ref, evidence)
    r["participants"].append(p)
    return pid


def reaction(rid, pages, record_type="specific_experiment", categories=None,
             location=None, name="dynamic kinetic asymmetric allenylation"):
    evidence = [ev(p) for p in pages]
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = F(categories or [], categories or [], evidence)
    r["location"] = F(location or {"pdf_pages": pages}, location or {"pdf_pages": pages}, evidence)
    r["reaction_name"] = F(name, name, evidence)
    r["transformation_type"] = F(
        "allenylic substitution forming a C-C bond",
        "allenylic C-C bond formation", evidence, method="structure_derivation",
        derivation="Derived from the displayed reactant/product connectivity; not an author-supplied taxonomy label.")
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [],
                   "unresolved": ["Supporting Information was mentioned but not supplied; full procedure/workup and dense structures remain unresolved."]}
    r["completeness"] = "partial"
    return r


def set_cc(r, evidence, reason="Displayed substitution joins the carbon nucleophile to the allenylic fragment."):
    change = {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"},
                                               {"atom_map": None, "element": "C"}],
              "before_order": 0, "after_order": 1, "change_type": "formed"}
    label = {"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1"]}
    r["classification"]["bond_changes"] = [F(None, change, evidence,
        method="structure_derivation", derivation=reason)]
    r["classification"]["labels"] = F(None, [label], evidence,
        method="structure_derivation", derivation=reason)
    r["classification"]["coverage"] = F(None, "partial", evidence,
        method="inference", derivation="Only the visually evident net new C-C connection is classified; complete atom mapping and hydrogen balance are unresolved.")
    r["classification"]["limitations"] = [
        "Dense source structures were not converted to complete mapped SMILES.",
        "Implicit hydrogen and leaving-group balance are unresolved; classification is partial."
    ]


def set_not_applicable(r, evidence, reason):
    r["classification"]["bond_changes"] = []
    r["classification"]["labels"] = F(None, [], evidence, method="inference", derivation=reason)
    r["classification"]["coverage"] = F(None, "not_applicable", evidence,
        method="inference", derivation=reason)
    r["classification"]["limitations"] = [reason]


def standard_setup(r, substrate, nucleophile, product_mid, evidence, pd_lig="R_L3",
                   cu_lig="RRp_L6", temp=-10, time=24, omit=()):
    omit = set(omit)
    p_sub = addp(r, substrate, "reactant", evidence, amount=0.30, equiv=1.5, ref=nucleophile)
    p_nuc = addp(r, nucleophile, "reactant", evidence, amount=0.20, equiv=1.0, ref=nucleophile)
    if product_mid:
        addp(r, product_mid, "product", evidence, product_role="major/desired product")
    if "Pd" not in omit:
        addp(r, M["Pd_allyl"], "catalyst", evidence, molpct=2.5, ref=nucleophile)
        addp(r, M[pd_lig], "ligand", evidence, molpct=5.5, ref=nucleophile)
    if "Cu" not in omit:
        addp(r, M["Cu_salt"], "catalyst", evidence, molpct=5.0, ref=nucleophile)
        addp(r, M[cu_lig], "ligand", evidence, molpct=5.5, ref=nucleophile)
    if "base" not in omit:
        addp(r, M["base"], "reagent", evidence, amount=0.24, equiv=1.2, ref=nucleophile)
    addp(r, M["THF"], "solvent", evidence)
    r["reaction_scale"] = F("2: 0.20 mmol", {"reference_molecule_id": nucleophile,
                                               "amount": q(0.20, "mmol")}, evidence)
    r["solvent_conditions"]["composition"] = F("THF", [M["THF"]], evidence)
    step = st.new_step("reaction", evidence)
    step["participant_ids"] = [p["participant_id"] for p in r["participants"]]
    step["operation"] = F("combined under the tabulated reaction conditions",
                          "combined under tabulated conditions", evidence)
    step["temperature"] = F(f"{temp} °C", q(temp, "°C"), evidence)
    step["time"] = F(f"{time} h", q(time, "h"), evidence)
    r["steps"] = [step]
    r["procedure"]["general_procedure_reference"] = F(
        "Table 1 footnote a / entry 7 as cited by scope tables", "Table 1 entry 7", evidence)
    set_cc(r, evidence)


def add_outcome(r, oid, kind, target, raw, norm, evidence, basis=None, scope=None,
                yield_type=None, status=None):
    o = st.new_outcome(oid, kind, [target] if target else [], evidence)
    o["value"] = F(raw, norm, evidence, status=status)
    if target and kind in ("yield", "selectivity", "product_amount"):
        o["product_role"] = F("reported product", "reported product", evidence)
    if yield_type:
        o["yield_type"] = F(yield_type, yield_type, evidence)
    if basis:
        o["measurement_basis"] = F(basis, basis, evidence)
    if scope:
        o["scope"] = F(scope, scope, evidence)
    r["outcomes"].append(o)
    return o


def pct_value(raw):
    if raw.startswith(">"):
        return q(float(raw[1:]), "%", "gt", raw)
    return q(float(raw), "%", raw_text=raw + "%")


GENERIC = f"{DOC}:scheme1C_generic"
g = reaction("scheme1C_generic", [1, 2], record_type="generic_scheme", categories=[],
             location={"pdf_pages": [1, 2], "regions": ["graphical abstract", "Scheme 1C"]},
             name="stereodivergent synergistic Pd/Cu-catalyzed allenylic substitution")
g["general_scheme_ids"] = []
g["review"]["unresolved"] = ["Generic R-group structures and exact experimental conditions are intentionally unresolved."]
set_cc(g, [ev(1), ev(2)])
bundle["reactions"].append(g)

# Scheme 1B: explicitly cited prior work, retained but not counted as this paper's experiment.
lit_sub = molecule("racemic allenyl acetate (Scheme 1B)", [ev(2)], "scheme1B_allenyl_acetate",
                   stereo={"kind": "racemic", "raw": "racemic"})
lit_nuc = molecule("1,3-dioxaborole α-hydroxyketone equivalent", [ev(2)], "scheme1B_boronate")
lit_prod = molecule("chiral allene product (Scheme 1B)", [ev(2)], "scheme1B_product",
                    stereo={"kind": "assigned", "assignment": "enantioenriched; absolute configuration not printed in this panel"})
M["Pd2dba3"] = molecule("Pd2dba3·CHCl3", [ev(2)], "Pd2dba3")
M["Trost"] = molecule("(R,R)-Trost ligand", [ev(2)], "RR_Trost",
                      stereo={"kind": "assigned", "assignment": "R,R",
                              "centers": [{"configuration": "R"}, {"configuration": "R"}]})
lr = reaction("scheme1B_literature_example", [2], record_type="literature_example", categories=[],
              location={"pdf_page": 2, "region": "Scheme 1B"},
              name="Pd-catalyzed asymmetric allenylation of α-hydroxyketones, Trost (2019)")
addp(lr, lit_sub, "reactant", [ev(2)])
addp(lr, lit_nuc, "reactant", [ev(2)])
addp(lr, lit_prod, "product", [ev(2)], product_role="reported literature product")
addp(lr, M["Pd2dba3"], "catalyst", [ev(2)], molpct=2.5, ref=lit_nuc)
addp(lr, M["Trost"], "ligand", [ev(2)], molpct=7.5, ref=lit_nuc)
addp(lr, M["base"], "reagent", [ev(2)])
addp(lr, M["THF"], "solvent", [ev(2)])
lr["solvent_conditions"]["composition"] = F("THF", [M["THF"]], [ev(2)])
s = st.new_step("reaction", [ev(2)]); s["participant_ids"] = [p["participant_id"] for p in lr["participants"]]
s["temperature"] = F("rt", "rt", [ev(2)]); lr["steps"] = [s]
add_outcome(lr, "yield", "yield", lit_prod, "50-77% yield", {"minimum": 50, "maximum": 77, "unit": "%", "raw_text": "50-77%"}, [ev(2)])
add_outcome(lr, "dr", "selectivity", lit_prod, "6:1->20:1 dr", "6:1 to >20:1", [ev(2)], scope={"metric": "dr"})
add_outcome(lr, "ee", "selectivity", lit_prod, "98% ee", q(98, "%"), [ev(2)], scope={"metric": "ee"})
set_cc(lr, [ev(2)])
bundle["reactions"].append(lr)

# Table 1 optimization. Entry 7 is also the top-left Scheme 3 experiment.
opt = [
    (1, "2a′", "Me", "R_L1", "RRp_L5", 0, "88", "1:1", ">99", "R,R"),
    (2, "2a′", "Me", "R_L2", "RRp_L5", 0, "92", "1:1", ">99", "R,R"),
    (3, "2a′", "Me", "R_L3", "RRp_L5", 0, "39", "15:1", ">99", "R,R"),
    (4, "2a′", "Me", "R_L4", "RRp_L5", 0, "48", "1:6", ">99", "R,S"),
    (5, "2a′", "Me", "R_L3", "RRp_L6", 0, "55", "11:1", ">99", "R,R"),
    (6, "2a′", "Me", "R_L3", "RRp_L6", -10, "50", ">20:1", ">99", "R,R"),
    (7, "2a", "Et", "R_L3", "RRp_L6", -10, "80", ">20:1", ">99", "R,R"),
    (8, "2a″", "Bn", "R_L3", "RRp_L6", -10, "62", "12:1", ">99", "R,R"),
    (9, "2a‴", "tBu", "R_L3", "RRp_L6", -10, "86", "13:1", ">99", "R,R"),
    (10, "2a", "Et", "R_L3", None, -10, "14", "1:2", ">99", "R,S"),
    (11, "2a", "Et", None, "RRp_L6", -10, "nr", None, None, None),
    (12, "2a", "Et", "R_L3", "RRp_L6", -10, "nr", None, None, None),
]
for entry, nuc_label, rgroup, pdlig, culig, temp, yld, dr, ee, config in opt:
    pages = [2, 3] if entry == 7 else [2]
    evidence = [ev(p) for p in pages]
    prod_mid = None
    if yld != "nr":
        suffix = {"2a′": "3aa_prime", "2a": "3aa", "2a″": "3aa_doubleprime", "2a‴": "3aa_tripleprime"}[nuc_label]
        prod_mid = product(f"({config})-{suffix.replace('_prime','′').replace('_doubleprime','″').replace('_tripleprime','‴')}",
                           pages, key=f"{safe(config)}_{suffix}")
    r = reaction(f"table1_entry_{entry:02d}", pages, categories=["optimization"] + (["control"] if entry >= 10 else []),
                 location={"pdf_page": 2, "table": 1, "entry": entry,
                           "duplicate_evidence": "Scheme 3 top-left" if entry == 7 else None})
    r["general_scheme_ids"] = [GENERIC]
    omit = []
    if entry == 10: omit = ["Cu"]
    if entry == 11: omit = ["Pd"]
    if entry == 12: omit = ["base"]
    standard_setup(r, M["1a"], M[nuc_label], prod_mid, evidence,
                   pd_lig=pdlig or "R_L3", cu_lig=culig or "RRp_L6", temp=temp, omit=omit)
    if yld == "nr":
        r["outcomes"] = []
        add_outcome(r, "no_reaction", "no_reaction", None, "nr", True, [ev(2)],
                    basis="Table 1 notation; detection limit not reported")
        set_not_applicable(r, [ev(2)], "No reaction was reported; no expected-product bond label is assigned.")
    else:
        add_outcome(r, "isolated_yield", "yield", prod_mid, yld + "%", q(float(yld), "%"), [ev(2)],
                    basis="isolated yield of all diastereoisomers", yield_type="isolated")
        add_outcome(r, "dr", "selectivity", prod_mid, dr + " dr", dr, [ev(2)],
                    basis="HPLC analysis", scope={"metric": "dr", "reported_major_configuration": config})
        add_outcome(r, "ee", "selectivity", prod_mid, ee + "% ee", pct_value(ee), [ev(2)],
                    basis="HPLC analysis", scope={"metric": "ee", "reported_configuration": config})
    bundle["reactions"].append(r)

# Remaining three Scheme 3 stereoisomers (the fourth is Table 1 entry 7).
scheme3 = [
    ("S,R", 77, "S_L3", "RRp_L6"),
    ("R,S", 88, "R_L3", "SSp_L6"),
    ("S,S", 84, "S_L3", "SSp_L6"),
]
for config, yld, pdlig, culig in scheme3:
    pm = product(f"({config})-3aa", [3], key=f"{safe(config)}_3aa")
    r = reaction(f"scheme3_{safe(config)}_3aa", [3], categories=["substrate_scope"],
                 location={"pdf_page": 3, "scheme": 3, "product": f"({config})-3aa"})
    r["general_scheme_ids"] = [GENERIC]
    standard_setup(r, M["1a"], M["2a"], pm, [ev(3)], pd_lig=pdlig, cu_lig=culig)
    add_outcome(r, "yield", "yield", pm, f"{yld}% yield", q(yld, "%"), [ev(3)])
    add_outcome(r, "dr", "selectivity", pm, ">20:1 dr", ">20:1", [ev(3)], scope={"metric": "dr"})
    add_outcome(r, "ee", "selectivity", pm, ">99% ee", pct_value(">99"), [ev(3)], scope={"metric": "ee"})
    bundle["reactions"].append(r)

# Table 2, two catalyst-controlled stereochemical series.
t2_rs = [("b",63,">20:1",">99",None),("c",66,">20:1",">99",None),("d",73,"14:1",">99",None),
         ("e",73,"10:1",">99",None),("f",79,"20:1","98",None),("g",94,">20:1",">99",None),
         ("h",82,">20:1","96",None),("i",48,"6:1",">99","b"),("j",52,"5:1",">99","b"),
         ("k",62,"7:1",">99","b"),("l",43,"9:1",">99",None),("m",83,"5:1",">99","b"),
         ("n",81,"5:1",">99","b"),("o",51,"12:1",">99","b"),("p",82,">20:1",">99",None)]
t2_ss = [("b",55,">20:1",">99",None),("c",72,">20:1",">99",None),("d",80,">20:1",">99",None),
         ("e",71,"15:1",">99",None),("f",74,">20:1","99",None),("g",94,">20:1",">99",None),
         ("h",79,">20:1",">99",None),("i",88,"8:1",">99",None),("j",75,"8:1","87",None),
         ("k",89,"5:1",">99",None),("l",77,"5:1",">99",None),("m",89,"3:1","96",None),
         ("n",85,"4:1",">99",None),("o",87,"4:1","99",None),("p",50,"9:1",">99","c")]
for config, rows in [("R,S", t2_rs), ("S,S", t2_ss)]:
    for letter, yld, dr, ee, foot in rows:
        pm = product(f"({config})-3{letter}a", [4], key=f"{safe(config)}_3{letter}a")
        pdlig = "R_L3" if config == "R,S" else "S_L3"
        if foot == "b": pdlig = "R_L4" if config == "R,S" else "S_L4"
        if foot == "c": pdlig = "S_DM_L3"
        r = reaction(f"table2_{safe(config)}_3{letter}a", [4], categories=["substrate_scope"],
                     location={"pdf_page": 4, "table": 2, "product": f"({config})-3{letter}a", "footnote": foot})
        r["general_scheme_ids"] = [GENERIC]
        standard_setup(r, M[f"1{letter}"], M["2a"], pm, [ev(4)], pd_lig=pdlig, cu_lig="SSp_L6")
        if foot == "b": r["notes"].append(F("Using DTBM-SegPhos instead of SegPhos.", "DTBM-SegPhos substitution", [ev(4)]))
        if foot == "c": r["notes"].append(F("Using DM-SegPhos instead of SegPhos.", "DM-SegPhos substitution", [ev(4)]))
        add_outcome(r, "yield", "yield", pm, f"{yld}% yield", q(yld, "%"), [ev(4)])
        add_outcome(r, "dr", "selectivity", pm, f"{dr} dr", dr, [ev(4)], scope={"metric": "dr"})
        add_outcome(r, "ee", "selectivity", pm, f"{ee}% ee", pct_value(ee), [ev(4)], scope={"metric": "ee"})
        bundle["reactions"].append(r)

# Table 3 imino ester scope.
t3_rs = [("b",66,">20:1","98"),("c",80,">20:1",">99"),("d",74,"12:1",">99"),("e",82,"20:1",">99"),
         ("f",63,">20:1",">99"),("g",63,"15:1",">99"),("h",66,"15:1","98"),("i",82,"17:1","93"),
         ("j",76,"14:1","99"),("k",91,"16:1",">99"),("l",60,"19:1",">99"),("m",83,">20:1",">99")]
t3_ss = [("b",78,">20:1","97"),("c",77,">20:1",">99"),("d",80,"11:1",">99"),("e",74,">20:1","97"),
         ("f",76,">20:1","95"),("g",78,">20:1",">99"),("h",66,"15:1",">99"),("i",69,">20:1",">99"),
         ("j",90,">20:1",">99"),("k",95,">20:1",">99"),("l",65,">20:1",">99"),("m",89,">20:1",">99")]
for config, rows in [("R,S", t3_rs), ("S,S", t3_ss)]:
    for letter, yld, dr, ee in rows:
        pm = product(f"({config})-3a{letter}", [5], key=f"{safe(config)}_3a{letter}")
        r = reaction(f"table3_{safe(config)}_3a{letter}", [5], categories=["substrate_scope"],
                     location={"pdf_page": 5, "table": 3, "product": f"({config})-3a{letter}"})
        r["general_scheme_ids"] = [GENERIC]
        standard_setup(r, M["1a"], M[f"2{letter}"], pm, [ev(5)],
                       pd_lig="R_L3" if config == "R,S" else "S_L3", cu_lig="SSp_L6")
        add_outcome(r, "yield", "yield", pm, f"{yld}% yield", q(yld, "%"), [ev(5)])
        add_outcome(r, "dr", "selectivity", pm, f"{dr} dr", dr, [ev(5)], scope={"metric": "dr"})
        add_outcome(r, "ee", "selectivity", pm, f"{ee}% ee", pct_value(ee), [ev(5)], scope={"metric": "ee"})
        if letter == "k":
            r["notes"].append(F("Absolute configuration determined by X-ray crystallography.",
                                  "X-ray assignment reported", [ev(5)]))
        bundle["reactions"].append(r)

# Table 4 aldimine/ketimine ester scope.
t4_rs = [("n",77,"10:1",">99"),("o",46,">20:1","99"),("p",72,"10:1","92"),
         ("q",69,">20:1","97"),("r",50,"12:1","98"),("s",58,"11:1",">99")]
t4_ss = [("n",71,"6:1",">99"),("o",78,">20:1","99"),("p",72,"11:1",">99"),
         ("q",56,">20:1",">99"),("r",70,">20:1","99"),("s",61,">20:1","96")]
for config, rows in [("R,S", t4_rs), ("S,S", t4_ss)]:
    for letter, yld, dr, ee in rows:
        pm = product(f"({config})-3a{letter}", [5], key=f"{safe(config)}_3a{letter}")
        r = reaction(f"table4_{safe(config)}_3a{letter}", [5], categories=["substrate_scope"],
                     location={"pdf_page": 5, "table": 4, "product": f"({config})-3a{letter}"})
        r["general_scheme_ids"] = [GENERIC]
        standard_setup(r, M["1a"], M[f"2{letter}"], pm, [ev(5)],
                       pd_lig="R_L3" if config == "R,S" else "S_L3", cu_lig="SSp_L6")
        add_outcome(r, "yield_after_hydrolysis", "yield", pm, f"{yld}% yield", q(yld, "%"), [ev(5)],
                    basis="yield determined after hydrolysis",
                    scope={"reported_scope": "after hydrolysis; exact assayed derivative identity unresolved"})
        add_outcome(r, "dr_after_TsCl", "selectivity", pm, f"{dr} dr", dr, [ev(5)],
                    basis="measured after N-protecting with TsCl", scope={"metric": "dr"})
        add_outcome(r, "ee_after_TsCl", "selectivity", pm, f"{ee}% ee", pct_value(ee), [ev(5)],
                    basis="measured after N-protecting with TsCl", scope={"metric": "ee"})
        r["notes"].append(F("The yield of product was determined after hydrolysis. The dr and ee values were measured after N-protecting with TsCl.",
                              "Table 4 footnote b", [ev(5)]))
        bundle["reactions"].append(r)

# Scheme 4: all eight stereoisomers from 1q/1q' and 2a.
scheme4 = [
    ("R,R,R", "1q", 81, ">20:1", "R_L3", "RRp_L6"),
    ("S,R,R", "1qprime", 66, "12:1", "R_L3", "RRp_L6"),
    ("R,R,S", "1q", 88, ">20:1", "R_L3", "SSp_L6"),
    ("S,R,S", "1qprime", 66, "18:1", "R_L3", "SSp_L6"),
    ("R,S,R", "1q", 68, ">20:1", "S_L3", "RRp_L6"),
    ("S,S,R", "1qprime", 63, ">20:1", "S_L3", "RRp_L6"),
    ("R,S,S", "1q", 60, ">20:1", "S_L3", "SSp_L6"),
    ("S,S,S", "1qprime", 74, ">20:1", "S_L3", "SSp_L6"),
]
for config, substrate_key, yld, dr, pdlig, culig in scheme4:
    pm = product(f"({config})-3qa", [5], key=f"{safe(config)}_3qa")
    r = reaction(f"scheme4_{safe(config)}_3qa", [5], categories=["substrate_scope"],
                 location={"pdf_page": 5, "scheme": 4, "product": f"({config})-3qa"})
    r["general_scheme_ids"] = [GENERIC]
    standard_setup(r, M[substrate_key], M["2a"], pm, [ev(5)], pd_lig=pdlig, cu_lig=culig)
    add_outcome(r, "yield", "yield", pm, f"{yld}% yield", q(yld, "%"), [ev(5)])
    add_outcome(r, "dr", "selectivity", pm, f"{dr} dr", dr, [ev(5)], scope={"metric": "dr"})
    r["notes"].append(F("No ee value is printed for this Scheme 4 arrow.",
                          "ee unreported in Scheme 4", [ev(5)]))
    bundle["reactions"].append(r)

# One monitored reaction sampled at six time points (Figure 1), not six independent syntheses.
pm_rs_3aa = molecules[f"{DOC}:mol:R_S_3aa"]["molecule_id"]
tc = reaction("figure1_time_course", [5, 6], categories=["mechanistic_study"],
              location={"pdf_pages": [5, 6], "figure": 1},
              name="time-course monitoring of rac-1a with 2a under standard conditions")
tc["general_scheme_ids"] = [GENERIC]
standard_setup(tc, M["1a"], M["2a"], pm_rs_3aa, [ev(5), ev(6)], pd_lig="R_L3", cu_lig="SSp_L6")
tc["steps"][0]["time"] = F("monitored through 10.0 h", q(10.0, "h"), [ev(5), ev(6)])
tc["notes"].append(F("Samples were analyzed at 0.5 h, 1.0 h, 3.0 h, 6.0 h, 8.0 h, and 10.0 h.",
                       {"sampling_times": [0.5, 1.0, 3.0, 6.0, 8.0, 10.0], "sampling_time_unit": "h"},
                       [ev(5), ev(6)]))
addp(tc, M["mesitylene"], "reagent", [ev(5), ev(6)], pid="p_characterization_mesitylene", stage="characterization")
for oid, kind, target, raw, basis, metric in [
    ("recovered_1a_yield_series", "recovery", M["1a"], "Figure 1A plotted recovered yields at six times; exact values not safely legible", "Figure 1A / HPLC context", "recovered yield"),
    ("recovered_1a_ee_series", "selectivity", M["1a"], "Figure 1A plotted ee values of recovered 1a; exact values not safely legible", "chiral HPLC, CHIRALPAK IC-3", "ee"),
    ("2a_conversion_series", "conversion", M["2a"], "Figure 1B plotted conversion of 2a; exact values not safely legible", "GC-MS using mesitylene as internal standard", "conversion"),
    ("3aa_ee_series", "selectivity", pm_rs_3aa, "Figure 1B plotted product ee at six times; exact values not safely legible", "chiral HPLC, CHIRALPAK IC-3", "ee"),
    ("3aa_dr_series", "selectivity", pm_rs_3aa, "Figure 1B plotted product dr at six times; exact values not safely legible", "chiral HPLC, CHIRALPAK IC-3", "dr"),
]:
    add_outcome(tc, oid, kind, target, raw, None, [ev(5), ev(6)], basis=basis,
                scope={"metric": metric, "timepoints_h": [0.5,1,3,6,8,10]}, status="unreadable")
tc["review"]["unresolved"].append("Exact numeric values in Figure 1 curves/bars were not transcribed because the plotted values are not printed as data labels.")
bundle["reactions"].append(tc)

# Scheme 5 is an author-proposed mechanism, not an experimental record.
mech = reaction("scheme5_proposed_mechanism", [5, 6], record_type="proposed_mechanism", categories=[],
                location={"pdf_pages": [5, 6], "scheme": 5},
                name="proposed mechanism of stereodivergent Pd/Cu-catalyzed DyKAT")
mech["general_scheme_ids"] = [GENERIC]
addp(mech, M["1a"], "reactant", [ev(5), ev(6)])
addp(mech, M["2a"], "reactant", [ev(5), ev(6)])
addp(mech, pm_rs_3aa, "product", [ev(5), ev(6)], product_role="illustrated product")
mech["steps"] = []
mech["outcomes"] = []
mech["transformation_type"] = F("proposed catalytic cycles and intermediates 4A-4D, 5, and 6",
                                  "proposed mechanism", [ev(5), ev(6)])
mech["classification"]["coverage"] = F(None, "unresolved", [ev(5), ev(6)], method="inference",
    derivation="The scheme is a proposed mechanism rather than a single observed net experiment; no elementary-step taxonomy is asserted.")
mech["review"]["unresolved"] = ["Mechanistic intermediates are proposed, not experimentally isolated in this scheme."]
bundle["reactions"].append(mech)

# Scheme 6: two stereoconvergent reactions and two no-racemization controls.
for enant, start_ee, yld in [("R", 96, 82), ("S", 94, 92)]:
    start = assigned(f"({enant})-1a", [ev(6)], key=f"{enant}_1a")
    r = reaction(f"scheme6_{enant}_1a_to_RS_3aa", [6], categories=["mechanistic_study"],
                 location={"pdf_page": 6, "scheme": 6, "row": f"({enant})-1a product experiment"})
    r["general_scheme_ids"] = [GENERIC]
    standard_setup(r, start, M["2a"], pm_rs_3aa, [ev(6)], pd_lig="R_L3", cu_lig="SSp_L6")
    add_outcome(r, "starting_material_ee", "selectivity", start, f"{start_ee}% ee", q(start_ee, "%"), [ev(6)], scope={"metric": "starting-material ee"})
    add_outcome(r, "yield", "yield", pm_rs_3aa, f"{yld}% yield", q(yld, "%"), [ev(6)])
    add_outcome(r, "dr", "selectivity", pm_rs_3aa, ">20:1 dr", ">20:1", [ev(6)], scope={"metric": "dr"})
    add_outcome(r, "ee", "selectivity", pm_rs_3aa, ">99% ee", pct_value(">99"), [ev(6)], scope={"metric": "ee"})
    bundle["reactions"].append(r)

for enant, start_ee, recovered_ee in [("R", 96, 93), ("S", 94, 92)]:
    start = molecules[f"{DOC}:mol:{enant}_1a"]["molecule_id"]
    r = reaction(f"scheme6_{enant}_1a_no_racemization_control", [6], categories=["control", "mechanistic_study"],
                 location={"pdf_page": 6, "scheme": 6, "row": f"({enant})-1a 10 mol% 2a control"},
                 name="no-racemization control under standard conditions")
    r["general_scheme_ids"] = [GENERIC]
    pstart = addp(r, start, "reactant", [ev(6)])
    addp(r, M["2a"], "reagent", [ev(6)], molpct=10, ref=start)
    addp(r, M["Pd_allyl"], "catalyst", [ev(6)], molpct=2.5, ref=start)
    addp(r, M["R_L3"], "ligand", [ev(6)], molpct=5.5, ref=start)
    addp(r, M["Cu_salt"], "catalyst", [ev(6)], molpct=5.0, ref=start)
    addp(r, M["SSp_L6"], "ligand", [ev(6)], molpct=5.5, ref=start)
    addp(r, M["base"], "reagent", [ev(6)], equiv=1.2, ref=start)
    addp(r, M["THF"], "solvent", [ev(6)])
    r["solvent_conditions"]["composition"] = F("THF", [M["THF"]], [ev(6)])
    step = st.new_step("reaction", [ev(6)]); step["participant_ids"] = [p["participant_id"] for p in r["participants"]]
    step["temperature"] = F("standard conditions (-10 °C)", q(-10, "°C"), [ev(6)], method="inherited_general_condition")
    step["time"] = F("standard conditions (24 h)", q(24, "h"), [ev(6)], method="inherited_general_condition")
    r["steps"] = [step]
    add_outcome(r, "starting_ee", "selectivity", start, f"{start_ee}% ee", q(start_ee, "%"), [ev(6)], scope={"metric": "starting-material ee"})
    add_outcome(r, "recovered_ee", "selectivity", start, f"recovered ({enant})-1a, {recovered_ee}% ee", q(recovered_ee, "%"), [ev(6)], scope={"metric": "recovered-material ee"})
    r["notes"].append(F("No obvious racemization was found.", "no obvious racemization", [ev(6)]))
    set_not_applicable(r, [ev(6)], "This is a substrate-racemization control with no reported synthetic product or newly formed bond.")
    bundle["reactions"].append(r)

bundle["molecules"] = list(molecules.values())
bundle["source_context"] = [
    {"context_id": f"{DOC}:page1_graphical_abstract", "document_id": DOC,
     "evidence_ids": [ev(1)], "kind": "duplicate_generic_presentation",
     "note": "Graphical abstract repeats the generic Pd/Cu stereodivergent coupling concept; linked to scheme1C_generic, not counted as another experiment."},
    {"context_id": f"{DOC}:page2_scheme1A_concept", "document_id": DOC,
     "evidence_ids": [ev(2)], "kind": "conceptual_background",
     "note": "Scheme 1A is a conceptual stereogenic-center difficulty chart, not a reaction experiment."},
    {"context_id": f"{DOC}:page2_scheme2_bioactive", "document_id": DOC,
     "evidence_ids": [ev(2)], "kind": "background_structures",
     "note": "Scheme 2 shows representative bioactive allene motifs and general allene chemistry; no present-paper experimental arrow is extracted."},
    {"context_id": f"{DOC}:figure1_numeric_plot", "document_id": DOC,
     "evidence_ids": [ev(5), ev(6)], "kind": "mechanistic_measurement_plot",
     "note": "Figure 1 values are represented in the time-course experiment, but exact plotted numeric series remain unreadable because data labels are absent."},
]

st.write(SUBMITTED, bundle)

visual = {
    "document_id": DOC,
    "scope": "original PDF pages 1-6",
    "pages_visually_opened": [1, 2, 3, 4, 5, 6],
    "page_inventory": {
        "1": ["title/abstract", "graphical abstract generic coupling", "introduction"],
        "2": ["Scheme 1A conceptual", "Scheme 1B cited literature example", "Scheme 1C current generic scheme", "Scheme 2 background", "Table 1 entries 1-12"],
        "3": ["discussion of Table 1", "Scheme 3 four stereoisomer arrows", "scope discussion"],
        "4": ["Table 2: 30 product entries", "scope prose"],
        "5": ["Table 3: 24 entries", "Table 4: 12 entries", "Scheme 4: eight arrows", "mechanistic/time-course prose"],
        "6": ["Scheme 5 proposed mechanism", "Figure 1 time-course plot", "Scheme 6 four arrows", "DFT discussion"],
    },
    "notes": ["All values are machine records and remain needs_review.",
              "No Supporting Information file was supplied or associated.",
              "Dense structures were not converted to guessed SMILES."],
}
st.write(ROOT / "runs/current/raw_model/visual_transcription.json", visual)
st.write(ROOT / "runs/current/raw_model/stereo_diagnostics.json", stereo_diagnostics)
st.write(ROOT / "runs/current/raw_model/source_mentions.json",
         st.source_mentions(dataset, ["1a", "2a", "3aa", "3ak", "3qa"], context=220))

# Coverage is generated from the actual page evidence and reaction assignments.
by_page = {p: [] for p in range(1, 7)}
for r in bundle["reactions"]:
    for p in range(1, 7):
        if ev(p) in r["evidence_ids"]:
            by_page[p].append(r["reaction_id"])
regions = {
    1: [{"region": "title, abstract, and graphical abstract", "finding": "Current-work generic Pd/Cu stereodivergent coupling; duplicate presentation, not a separate experiment.", "evidence_ids": [ev(1)]},
        {"region": "introduction columns", "finding": "Background discussion only; no specific present-paper experiment.", "evidence_ids": [ev(1)]}],
    2: [{"region": "Scheme 1A-C", "finding": "Conceptual background, one cited Trost literature reaction, and one current-work generic scheme were distinguished.", "evidence_ids": [ev(2)]},
        {"region": "Scheme 2", "finding": "Bioactive examples/background; not present-paper synthesis records.", "evidence_ids": [ev(2)]},
        {"region": "Table 1 entries 1-12 and footnotes", "finding": "Twelve optimization/control experiments, including two no-reaction controls, were recorded row by row.", "evidence_ids": [ev(2)]}],
    3: [{"region": "Scheme 3", "finding": "Four stereoisomer arrows inspected; (R,R)-3aa is the same experiment as Table 1 entry 7 and evidence was merged, while three additional experiments were created.", "evidence_ids": [ev(3)]},
        {"region": "results prose", "finding": "Catalyst-control and scope discussion mapped to records; no additional independent reaction arrow.", "evidence_ids": [ev(3)]}],
    4: [{"region": "Table 2 upper (R,S) series", "finding": "Fifteen distinct substrate-scope entries transcribed with yield, dr, ee and ligand footnotes.", "evidence_ids": [ev(4)]},
        {"region": "Table 2 lower (S,S) series", "finding": "Fifteen distinct substrate-scope entries transcribed with yield, dr, ee and ligand footnotes.", "evidence_ids": [ev(4)]}],
    5: [{"region": "Table 3", "finding": "Twenty-four imino-ester scope entries transcribed (12 per stereochemical series).", "evidence_ids": [ev(5)]},
        {"region": "Table 4", "finding": "Twelve aldimine/ketimine scope entries transcribed; hydrolysis/TsCl measurement scope retained.", "evidence_ids": [ev(5)]},
        {"region": "Scheme 4", "finding": "Eight stereoisomer-forming arrows transcribed independently.", "evidence_ids": [ev(5)]},
        {"region": "mechanistic/time-course prose", "finding": "Proposed mechanism kept separate; one monitored reaction linked across pages 5-6.", "evidence_ids": [ev(5)]}],
    6: [{"region": "Scheme 5", "finding": "Author-proposed catalytic mechanism; recorded as proposed_mechanism, not experimental proof.", "evidence_ids": [ev(6)]},
        {"region": "Figure 1", "finding": "One time-course experiment at six sampling times; exact plotted numeric series lacked printed labels and remain unreadable.", "evidence_ids": [ev(6)]},
        {"region": "Scheme 6", "finding": "Two enantiopure-substrate product experiments and two no-racemization controls transcribed arrow by arrow.", "evidence_ids": [ev(6)]},
        {"region": "DFT discussion", "finding": "Computational transition-state energy discussion is source context, not a synthetic experiment.", "evidence_ids": [ev(6)]}],
}
unresolved = {
    1: ["Generic graphical structures contain R groups; no concrete SMILES assigned."],
    2: ["Complete substrate/product connectivity is drawing-only; SI procedures unavailable."],
    3: ["Scheme 3 yield measurement basis is not printed."],
    4: ["Dense scope structures and yield/selectivity measurement bases are not fully reported on this page."],
    5: ["Dense scope structures remain without complete SMILES; Table 4 assayed derivative identities after hydrolysis/TsCl are not fully specified in the main text."],
    6: ["Exact Figure 1 numeric series are not printed as labels; Scheme 6 yield/selectivity measurement bases are unreported; DFT pathways are not experimental reactions."],
}
coverage = {"document_id": DOC, "semantic_accuracy_verified": False,
            "pages": [{"page": p, "state": "attempted", "regions": regions[p],
                       "record_ids": by_page[p], "unresolved": unresolved[p]}
                      for p in range(1, 7)]}
st.write(ROOT / "runs/current/coverage/pages_0001_0006.json", coverage)

print(json.dumps({"molecules": len(bundle["molecules"]),
                  "reactions": len(bundle["reactions"]),
                  "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in bundle["reactions"]),
                  "needs_review": sum(r["review"]["status"] == "needs_review" for r in bundle["reactions"])},
                 ensure_ascii=False, indent=2))
