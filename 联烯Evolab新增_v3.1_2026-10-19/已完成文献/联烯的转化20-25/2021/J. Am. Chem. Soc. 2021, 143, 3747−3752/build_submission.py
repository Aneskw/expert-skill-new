#!/usr/bin/env python3
"""Build the source-derived v3.1 submission for JACS 2021, 143, 3747-3752.

The script is intentionally document-specific.  It transcribes the six supplied
main-text pages and does not use an SI or any previous paper's reaction records.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / "runs/current"
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_bf6cf3c425610cb9"
EV = {i: f"{DOC}:page{i}:image" for i in range(1, 7)}
EVT = {i: f"{DOC}:page{i}:text" for i in range(1, 7)}


def dump(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def F(raw=None, normalized=None, evidence=(), method="original_extraction",
      derivation=None, status=None):
    return st.fact(raw, normalized, evidence, method=method,
                   derivation=derivation, status=status)


dataset = json.loads((RUN / "dataset.json").read_text())
bundle = st.empty(dataset["documents"])
molecules = {}


def slug(value):
    return re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")


def molecule(key, label, evidence, name=None, smiles=None, chirality="unresolved",
             geometry="unresolved", representation=None, assignment=None,
             mixture=None, double_bonds=None, identity_note=None):
    """Create or extend one paper-local molecule registry item."""
    mid = f"{DOC}:mol_{slug(key)}"
    evidence = list(dict.fromkeys(evidence))
    if mid in molecules:
        occ = molecules[mid]["occurrences"][0]
        occ["evidence_ids"] = list(dict.fromkeys(occ["evidence_ids"] + evidence))
        return mid
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, evidence)
    m["names"]["original"] = F(name or label, name or label, evidence)
    m["names"]["english"] = F(name or label, name or label, evidence)
    m["names"]["abbreviations"] = F([label], [label], evidence)
    m["names"]["aliases"] = F()
    m["cas_numbers"] = F()
    if smiles:
        diag = st.stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = F(
            {"source": "structure drawing", "candidate": smiles},
            diag["canonical_isomeric_smiles"], evidence,
            method="structure_derivation",
            derivation="Candidate connectivity was transcribed from the supplied drawing; RDKit parsing/normalization does not verify source identity.",
        )
        rep = representation or "partial"
        m["stereochemistry"]["isomeric_smiles"] = F(
            smiles, diag["canonical_isomeric_smiles"], evidence,
            method="structure_derivation",
            derivation="Canonical isomeric SMILES generated from the source-derived candidate graph; allene axial support is limited.",
        )
        m["stereochemistry"]["identity_key"] = F(
            None, diag["stereo_graph_key"], evidence,
            method="structure_derivation",
            derivation="RDKit stereo graph key from the candidate SMILES; not a sample identity.",
        )
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = F()
        rep = representation or "not_encoded"
    m["stereochemistry"]["chirality_status"] = F(
        chirality, chirality, evidence, method="inference",
        derivation="Status preserves the source label/drawing; lack of a stereo tag is not treated as proof of achirality.",
    )
    m["stereochemistry"]["geometry_status"] = F(
        geometry, geometry, evidence, method="inference",
        derivation="Status records whether alkene/allene or ring geometry is actually assigned in the supplied source.",
    )
    m["stereochemistry"]["representation_status"] = F(
        rep, rep, evidence, method="inference",
        derivation="Complete stereochemical encoding was not claimed unless supported by the candidate graph and source assignment.",
    )
    if assignment is not None:
        m["stereochemistry"]["assignment"] = F(assignment, assignment, evidence)
    if mixture is not None:
        m["stereochemistry"]["mixture"] = F(mixture, mixture, evidence)
    if double_bonds is not None:
        m["stereochemistry"]["double_bonds"] = F(double_bonds, double_bonds, evidence)
    m["identity_note"] = F(
        identity_note or "Identity is paper-local; CAS was not reported and no external identity completion was performed.",
        identity_note or "Identity is paper-local; CAS was not reported and no external identity completion was performed.",
        evidence,
    )
    molecules[mid] = m
    return mid


def generic_molecule(key, label, evidence, fixed, variables):
    mid = molecule(key, label, evidence, name=label, representation="not_encoded")
    m = molecules[mid]
    m["representation_type"] = "generic_template"
    m["smiles"] = F()
    m["generic_template"] = {
        "fixed_structure": F(label, fixed, evidence),
        "variables": [
            {
                "label": v,
                "attachment_points": F(f"Scheme attachment point for {v}", f"source-defined attachment point for {v}", evidence),
                "definition": F(f"{v} as drawn", None, evidence, status="unreadable"),
                "allowed_values": F(),
                "constraints": F(),
            }
            for v in variables
        ],
        "instance_mapping": [],
        "reason": F("generic R-group drawing", "Source is a generic scheme rather than a single concrete molecule.", evidence),
    }
    return mid


# Generic context from Scheme 1.
g_allyl = generic_molecule("generic_allyltriflamide", "racemic alpha-allyl amine template", [EV[1]],
                           "allylamine/triflamide framework shown in Scheme 1B", ["R", "R'", "R''", "PG"])
g_allene = generic_molecule("generic_allene", "allene template", [EV[1]], "C=C=C with source-defined substituents", ["R1", "R2"])
g_product = generic_molecule("generic_tetrahydropyridine", "chiral tetrahydropyridine template", [EV[1]],
                             "six-membered N-heterocycle from formal (4+2) annulation", ["R", "R'", "R''", "PG"])
g_recovered = generic_molecule("generic_chiral_allyltriflamide", "chiral recovered alpha-allyl amine template", [EV[1]],
                               "enantioenriched allylamine/triflamide framework shown in Scheme 1B", ["R", "R'", "R''", "PG"])
for mid, description in [(g_product, "generic product is explicitly described as chiral"),
                         (g_recovered, "generic recovered starting material is explicitly described as chiral")]:
    molecules[mid]["stereochemistry"]["chirality_status"] = F(
        "partially_specified", "partially_specified", [EV[1]], method="inference",
        derivation="Scheme 1B specifies chirality but not a complete absolute configuration for the generic R-group template.")
    molecules[mid]["stereochemistry"]["assignment"] = F(description, {"source_assignment": description, "absolute_configuration": None}, [EV[1]])
g_prev_s = generic_molecule("prior_alpha_branched_amine", "prior-work alpha-branched amine template", [EV[1]],
                            "alpha-branched benzyl/alkyl amine template", ["R", "NHR"])
g_prev_p = generic_molecule("prior_functionalized_amine", "prior-work functionalized amine template", [EV[1]],
                            "iodinated or arylated kinetic-resolution product template", ["R", "NHR", "X"])

# Concrete allene candidates that are legible on page 2.  Product/substrate dense
# structures remain label-based below rather than forcing unsafe full SMILES.
allene_smiles = {
    "2a": "C=C=C1CCCCC1",
    "2b": "C=C=C(CCCC)CCCC",
    "2c": "C=C=CC1CCCCC1",
    "2d": "C=C=CCO[Si](C)(C)C(C)(C)C",
    "2e": "C=C=CCN1C(=O)c2ccccc2C1=O",
    "2f": "C=C=CCOC(C)=O",
    "2g": "C=C=Cc1ccccc1",
    "2h": "C=C=Cc1ccc(C(F)(F)F)cc1",
    "2i": "CC(=C=Cc1ccccc1)C",
    "free_allenol": "C=C=CCO",
}
allene_names = {
    "2a": "vinylidenecyclohexane (paper label 2a)",
    "2b": "5-vinylidenenonane",
    "2c": "propa-1,2-dien-1-ylcyclohexane",
    "2d": "silyl-protected buta-2,3-dien-1-ol 2d",
    "2e": "phthalimide-functionalized allene 2e",
    "2f": "acetate-functionalized allene 2f",
    "2g": "phenyl-substituted allene 2g",
    "2h": "trifluoromethylphenyl-substituted allene 2h",
    "2i": "trisubstituted allene 2i",
    "free_allenol": "buta-2,3-dien-1-ol",
}
for key, smi in allene_smiles.items():
    molecule(key, key if key != "free_allenol" else "buta-2,3-dien-1-ol", [EV[2], EVT[3]],
             allene_names[key], smiles=smi, chirality="achiral", geometry="not_applicable",
             assignment={"basis": "candidate graph has duplicate substituents at at least one allene terminus; no source axial assignment"},
             identity_note="Candidate SMILES is structure-derived from the main-text drawing/name and remains unverified without SI.")


def rac_substrate(label, evidence):
    return molecule(f"rac_{label}", f"rac-{label}", evidence, f"racemic allyltriflamide {label}",
                    chirality="racemic", geometry="unspecified", representation="not_encoded",
                    mixture={"type": "racemate", "ratio": "not stated"},
                    identity_note="The source explicitly marks the input racemic; dense connectivity and alkene geometry are not encoded without SI.")


def recovered_substrate(label, evidence, explicit_s=False):
    assignment = {"configuration": "S", "basis": "Scheme 3 generic recovered material is labeled (S)-1"}
    if label == "1b":
        assignment["basis"] = "X-ray structure on page 3 is labeled (S)-1b"
    if label == "1a":
        assignment["basis"] = "Page 4 prose and Scheme 5 explicitly label (S)-1a"
    return molecule(f"s_{label}", f"(S)-{label}" if explicit_s or label in ("1a", "1b") else label,
                    evidence, f"enantioenriched recovered allyltriflamide {label}",
                    chirality="specified", geometry="unspecified", representation="not_encoded",
                    assignment=assignment,
                    identity_note="Absolute configuration is retained textually; no unsafe @/@@ encoding is supplied for the dense source drawing.")


for label in [f"1{x}" for x in "abcdefghijklm"]:
    rac_substrate(label, [EV[2], EV[3]])
    recovered_substrate(label, [EV[2], EV[3], EVT[4]], explicit_s=(label in ("1a", "1b")))


def product(label, evidence, scheme3=False, assignment=None, geometry="unspecified"):
    chirality = "specified" if scheme3 or assignment else "partially_specified"
    assn = assignment or ({"configuration": "R", "basis": "Scheme 3 generic product is labeled (R)-3"} if scheme3 else
                          {"basis": "source drawing uses a stereobond but gives no explicit R/S label for this panel"})
    return molecule(label, label, evidence, f"tetrahydropyridine {label}", chirality=chirality,
                    geometry=geometry, representation="not_encoded", assignment=assn,
                    identity_note="Product drawing and stereochemical annotation are retained textually; complete stereochemical SMILES is not assigned without SI.")


for label in ["3aa", "3ab", "3ac", "3ad", "3ae", "3af", "3ag", "3ah", "3ai"]:
    geom = "specified" if label in {"3ac", "3ad", "3ae", "3af", "3ag", "3ah", "3ai"} else "unspecified"
    db = None
    if geom == "specified":
        db = {"source_metric": "Z:E ratio reported for the product panel", "exact double-bond atom indices": "not encoded"}
    molecule(label, label, [EV[2], EVT[3], EVT[4]], f"tetrahydropyridine {label}",
             chirality="partially_specified", geometry=geom, representation="not_encoded",
             assignment={"basis": "source product drawing/ee; absolute R/S is not explicitly assigned in Scheme 2"},
             double_bonds=db,
             identity_note="Scheme 2 structure is source-identified, but dense connectivity/stereo is not encoded without SI.")

scheme3_products = ["3ba", "3ca", "3da", "3ea", "3fa", "3ga", "3ha", "3ia", "3ja", "3ka", "3la", "3ma", "3hc"]
for label in scheme3_products:
    assn = {"configuration": "R", "basis": "Scheme 3 generic product is labeled (R)-3"}
    if label == "3ba": assn["basis"] = "X-ray structure on page 3 is explicitly labeled (R)-3ba"
    product(label, [EV[3], EVT[4]], scheme3=True, assignment=assn,
            geometry="specified" if label == "3hc" else "unspecified")

for label, name, chi, geom, assn in [
    ("4", "hydrogenation/isomerization product 4", "partially_specified", "unspecified", {"basis": "stereobonds shown; prose says fully diastereoselective"}),
    ("5", "oxidation product 5", "partially_specified", "unspecified", {"basis": "stereobond shown; no R/S label"}),
    ("6", "detriflylated product 6", "partially_specified", "unspecified", {"basis": "stereobond shown; no R/S label"}),
    ("7", "chiral keto amino product 7", "specified", "not_applicable", {"configuration": "S", "basis": "derived from explicitly labeled (S)-1a; source drawing retains stereobond"}),
    ("8", "N-propargyl enyne 8", "specified", "unspecified", {"configuration": "S", "basis": "source drawing and sequence from (S)-1a"}),
    ("9", "piperidine 9", "partially_specified", "unspecified", {"basis": "source drawing shows relative stereochemistry; prose reports a single diastereoisomer"}),
    ("s_3aa", "(S)-3aa", "specified", "unspecified", {"configuration": "S", "basis": "Scheme 5 explicitly labels (S)-3aa"}),
]:
    molecule(label, "(S)-3aa" if label == "s_3aa" else label, [EV[3], EVT[4]], name,
             chirality=chi, geometry=geom, representation="not_encoded", assignment=assn)


# Reagents, ligands, oxidants and solvents.  They remain name-resolved and do not
# receive guessed CAS values or structures.
chemical_names = {
    "pd_oac2": "palladium(II) acetate", "cu_oac2_h2o": "copper(II) acetate monohydrate",
    "cs2co3": "cesium carbonate", "dmso": "dimethyl sulfoxide", "toluene": "toluene",
    "boc_l_phe_nhome": "Boc-L-Phe-NHOMe", "boc_d_phe_nhome": "Boc-D-Phe-NHOMe",
    "boc_val_nhome": "Boc-Val-NHOMe", "boc_leu_nhome": "Boc-Leu-NHOMe",
    "boc_phe_oh": "Boc-Phe-OH", "boc_val_oh": "Boc-Val-OH", "boc_leu_oh": "Boc-Leu-OH",
    "boc_tleu_oh": "Boc-tert-Leu-OH", "boc_ile_oh": "Boc-Ile-OH",
    "tcboc_phe_oh": "TcBoc-Phe-OH", "ac_phe_oh": "Ac-Phe-OH",
    "ffbz_phe_oh": "2,6-F,F-Bz-Phe-OH", "agoac": "silver acetate", "ag2co3": "silver carbonate",
    "thf": "tetrahydrofuran", "dce": "1,2-dichloroethane", "iproh": "isopropanol",
    "tamyl_oh": "tert-amyl alcohol", "h2": "hydrogen", "pd_c": "palladium on carbon",
    "etoac": "ethyl acetate", "rucl3": "ruthenium trichloride", "naio4": "sodium periodate",
    "water": "water", "acn": "acetonitrile", "red_al": "Red-Al", "ccl4": "carbon tetrachloride",
    "k2co3": "potassium carbonate", "propargyl_bromide": "propargyl bromide",
    "ir_catalyst": "[IrCp*Cl2]2",
}
for key, name in chemical_names.items():
    molecule(key, name, [EV[2], EV[3]], name, chirality="unresolved", geometry="unresolved",
             representation="not_encoded")


def qty(number, unit, raw=None, qualifier=None):
    value = {"number": number, "unit": unit}
    if qualifier: value["qualifier"] = qualifier
    if raw: value["raw_text"] = raw
    return value


def participant(pid, mid, role, evidence, *, amount=None, mass=None, volume=None,
                equiv=None, molpct=None, reference=None, product_role=None,
                stage="reaction"):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    if product_role:
        p["product_role"] = F(product_role, product_role, evidence)
    fields = {"amount": amount, "mass": mass, "volume": volume,
              "equivalents": equiv, "mol_percent": molpct}
    for key, item in fields.items():
        if item is not None:
            raw = item.get("raw_text") if isinstance(item, dict) else str(item)
            p["charges"][key] = F(raw or item, item, evidence)
    if reference is not None:
        p["equivalence_reference"] = F(reference, reference, evidence)
    return p


def base_reaction(suffix, evidence, record_type, categories, location, name, transformation):
    rid = f"{DOC}:{suffix}"
    r = st.new_reaction(rid, DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = F(categories, categories, evidence)
    r["location"] = F(location, location, evidence)
    r["reaction_name"] = F(name, name, evidence)
    r["transformation_type"] = F(
        transformation, transformation, evidence, method="inference",
        derivation="Concise transformation label derived from the source scheme/title; not asserted as an author-defined ontology term.",
    )
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    return r


def reaction_step(r, evidence, temperature=None, time=None, atmosphere=None, operation="reaction under the displayed conditions"):
    s = st.new_step(f"{r['reaction_id'].split(':')[-1]}_reaction", evidence)
    s["participant_ids"] = [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"]
    s["operation"] = F(operation, operation, evidence)
    if temperature is not None:
        s["temperature"] = F(temperature[0], temperature[1], evidence)
    if time is not None:
        s["time"] = F(time[0], time[1], evidence)
    if atmosphere is not None:
        s["atmosphere"] = F(atmosphere, atmosphere, evidence)
    r["steps"].append(s)


def outcome(r, oid, kind, targets, evidence, raw, normalized, *, role=None,
            yield_type=None, basis=None, scope=None, method="original_extraction", derivation=None):
    o = st.new_outcome(oid, kind, targets, evidence)
    o["value"] = F(raw, normalized, evidence, method=method, derivation=derivation)
    if role is not None: o["product_role"] = F(role, role, evidence)
    if yield_type is not None: o["yield_type"] = F(yield_type, yield_type, evidence)
    if basis is not None: o["measurement_basis"] = F(basis, basis, evidence)
    if scope is not None: o["scope"] = F(scope, scope, evidence)
    r["outcomes"].append(o)


def classify_changes(r, changes, evidence, note):
    labels = {}
    facts = []
    for i, ch in enumerate(changes, 1):
        ch = dict(ch)
        ch["change_id"] = f"b{i}"
        facts.append(F(None, ch, evidence, method="structure_derivation", derivation=note))
        if ch["change_type"] in {"formed", "order_increased", "hydrogen_gain"}:
            lab = st.bond_label(ch["atoms"][0]["element"], ch["atoms"][1]["element"], ch["after_order"])
            labels.setdefault(lab["code"], dict(lab, supporting_change_ids=[]))["supporting_change_ids"].append(ch["change_id"])
    r["classification"]["bond_changes"] = facts
    r["classification"]["labels"] = F(None, list(labels.values()), evidence, method="structure_derivation", derivation=note)
    r["classification"]["coverage"] = F(None, "partial", evidence, method="inference",
                                              derivation="Only visually secure net bond classes are asserted; complete atom mapping is unavailable.")
    r["classification"]["limitations"] = [note, "No atom-mapped reaction SMILES or three-model vote was fabricated."]


def formal_annulation_classification(r, evidence):
    classify_changes(r, [
        {"atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
         "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "N"}],
         "before_order": 0, "after_order": 1, "change_type": "formed"},
    ], evidence, "Formal (4+2) annulation drawing supports one new C-C and one new C-N connection; other bond-order/H changes remain unresolved.")


def no_reaction_classification(r, evidence):
    r["classification"]["bond_changes"] = []
    r["classification"]["labels"] = F("no observed target transformation", [], evidence,
                                             method="inference", derivation="Negative/complex-mixture result; expected bonds are not classified as observed.")
    r["classification"]["coverage"] = F("no observed target transformation", "not_applicable", evidence,
                                               method="inference", derivation="No target annulation product was reported.")
    r["classification"]["limitations"] = ["Complex mixture components were not identified in the main text."]


def common_annulation_participants(r, substrate, allene, product_mid, recovered_mid, evidence,
                                   ligand="boc_l_phe_nhome", temperature=70, solvent="toluene",
                                   cu_equiv=2, include_base=True, include_dmso=True):
    ref = substrate
    r["participants"] = [
        participant("substrate", substrate, "reactant", evidence, amount=qty(0.1, "mmol", "0.1 mmol")),
        participant("allene", allene, "reactant", evidence, amount=qty(0.1, "mmol", "0.1 mmol")),
        participant("product", product_mid, "product", evidence, product_role="major cycloadduct"),
        participant("recovered_starting_material", recovered_mid, "reactant", evidence, product_role="recovered starting material"),
        participant("pd_catalyst", molecules[f"{DOC}:mol_pd_oac2"]["molecule_id"], "catalyst", evidence,
                    molpct=qty(10, "mol%", "10 mol%"), reference=ref),
        participant("ligand", molecules[f"{DOC}:mol_{slug(ligand)}"]["molecule_id"], "ligand", evidence,
                    molpct=qty(40, "mol%", "40 mol%"), reference=ref),
        participant("oxidant", molecules[f"{DOC}:mol_cu_oac2_h2o"]["molecule_id"], "reagent", evidence,
                    equiv=qty(cu_equiv, "equiv", f"{cu_equiv} equiv"), reference=ref),
    ]
    if include_base:
        r["participants"].append(participant("base", molecules[f"{DOC}:mol_cs2co3"]["molecule_id"], "reagent", evidence,
                                                 equiv=qty(1.5, "equiv", "1.5 equiv"), reference=ref))
    if include_dmso:
        r["participants"].append(participant("dmso", molecules[f"{DOC}:mol_dmso"]["molecule_id"], "solvent", evidence,
                                                 equiv=qty(15, "equiv", "15 equiv"), reference=ref))
    r["participants"].append(participant("main_solvent", molecules[f"{DOC}:mol_{slug(solvent)}"]["molecule_id"], "solvent", evidence,
                                             volume=qty(1.2, "mL", "1.2 mL")))
    r["reaction_scale"] = F("substrate 0.1 mmol", {"reference_molecule_id": substrate, "amount": qty(0.1, "mmol")}, evidence)
    r["solvent_conditions"]["composition"] = F(
        f"{solvent}; DMSO" if include_dmso else solvent,
        [molecules[f"{DOC}:mol_{slug(solvent)}"]["molecule_id"]] + ([molecules[f"{DOC}:mol_dmso"]["molecule_id"]] if include_dmso else []), evidence)
    if solvent == "toluene" and include_dmso:
        r["solvent_conditions"]["mixture_ratio"] = F("toluene/DMSO (10:1) in model-condition prose", "10:1", [EVT[2]])
    reaction_step(r, evidence, (f"{temperature} °C", qty(temperature, "°C")), ("24 h", qty(24, "h")), "air")
    formal_annulation_classification(r, evidence)


reactions = []

# Scheme 1A: prior work only.
r = base_reaction("scheme1a_prior_work", [EV[1]], "literature_example", ["other"],
                  "Scheme 1A", "Previous work: kinetic resolution of alpha-branched amines",
                  "prior-work C-H iodination/cross-coupling")
r["participants"] = [participant("prior_substrate", g_prev_s, "reactant", [EV[1]]),
                     participant("prior_product", g_prev_p, "product", [EV[1]])]
r["review"]["unresolved"] = ["Contextual prior-work scheme is not a current-paper experiment and has no experimental quantities on this page."]
reactions.append(r)

# Scheme 1B: this work's generic transformation.
r = base_reaction("scheme1b_generic_annulation", [EV[1]], "generic_scheme", ["other"],
                  "Scheme 1B", "Kinetic resolution of allyl amines via C-H functionalization with allenes",
                  "formal (4+2) annulation")
r["participants"] = [participant("generic_allyl", g_allyl, "reactant", [EV[1]]),
                     participant("generic_allene", g_allene, "reactant", [EV[1]]),
                     participant("generic_product", g_product, "product", [EV[1]]),
                     participant("generic_recovered", g_recovered, "reactant", [EV[1]], product_role="recovered chiral allyl amine")]
formal_annulation_classification(r, [EV[1]])
r["review"]["unresolved"] = ["Generic R groups and conditions are intentionally not expanded into a concrete experiment."]
reactions.append(r)
generic_id = r["reaction_id"]


# Table 1 optimization rows.
table_rows = [
    (1, "none", "boc_l_phe_nhome", 2, "cu_oac2_h2o", "toluene", True, True, 49, 94, 90, 100),
    (2, "Boc-Val-NHOMe", "boc_val_nhome", 2, "cu_oac2_h2o", "toluene", True, True, 45, 82, 66, 20),
    (3, "Boc-Leu-NHOMe", "boc_leu_nhome", 2, "cu_oac2_h2o", "toluene", True, True, 37, 88, 52, 26),
    (4, "Boc-Phe-OH", "boc_phe_oh", 2, "cu_oac2_h2o", "toluene", True, True, 63, 58, 99, 18),
    (5, "Boc-Val-OH", "boc_val_oh", 2, "cu_oac2_h2o", "toluene", True, True, 59, 58, 84, 10),
    (6, "Boc-Leu-OH", "boc_leu_oh", 2, "cu_oac2_h2o", "toluene", True, True, 64, 54, 98, 14),
    (7, "Boc-t-Leu-OH", "boc_tleu_oh", 2, "cu_oac2_h2o", "toluene", True, True, 58, 50, 70, 6),
    (8, "Boc-Ile-OH", "boc_ile_oh", 2, "cu_oac2_h2o", "toluene", True, True, 65, 50, 92, 9),
    (9, "TcBoc-Phe-OH", "tcboc_phe_oh", 2, "cu_oac2_h2o", "toluene", True, True, 64, 20, 36, 2),
    (10, "Ac-Phe-OH", "ac_phe_oh", 2, "cu_oac2_h2o", "toluene", True, True, 31, 48, 22, 4),
    (11, "2,6-F,F-Bz-Phe-OH", "ffbz_phe_oh", 2, "cu_oac2_h2o", "toluene", True, True, 47, 60, 54, 7),
    (12, "0.5 equiv Cu(OAc)2.H2O", "boc_l_phe_nhome", 0.5, "cu_oac2_h2o", "toluene", True, True, 36, 94, 52, 54),
    (13, "AgOAc as oxidant", "boc_l_phe_nhome", 2, "agoac", "toluene", True, True, 9, 86, 8, 14),
    (14, "Ag2CO3 as oxidant", "boc_l_phe_nhome", 2, "ag2co3", "toluene", True, True, 6, 92, 6, 25),
    (15, "THF as solvent", "boc_l_phe_nhome", 2, "cu_oac2_h2o", "thf", True, True, 25, 84, 28, 15),
    (16, "DCE as solvent", "boc_l_phe_nhome", 2, "cu_oac2_h2o", "dce", True, True, 39, 86, 56, 23),
    (17, "i-PrOH as solvent", "boc_l_phe_nhome", 2, "cu_oac2_h2o", "iproh", True, True, 51, 78, 80, 20),
    (18, "t-AmylOH as solvent", "boc_l_phe_nhome", 2, "cu_oac2_h2o", "tamyl_oh", True, True, 55, 80, 99, 40),
    (19, "no base", "boc_l_phe_nhome", 2, "cu_oac2_h2o", "toluene", False, True, 4, 44, 2, 3),
    (20, "no DMSO", "boc_l_phe_nhome", 2, "cu_oac2_h2o", "toluene", True, False, 12, 72, 10, 7),
]
rac1a = molecules[f"{DOC}:mol_rac_1a"]["molecule_id"]
s1a = molecules[f"{DOC}:mol_s_1a"]["molecule_id"]
m2a = molecules[f"{DOC}:mol_2a"]["molecule_id"]
m3aa = molecules[f"{DOC}:mol_3aa"]["molecule_id"]
for row in table_rows:
    entry, deviation, ligand, ox_eq, oxidant, solvent, base, dmso, conv, ee_p, ee_s, sf = row
    evidence = [EV[2], EVT[2]]
    r = base_reaction(f"table1_entry_{entry:02d}", evidence, "specific_experiment", ["optimization"],
                      f"Table 1 entry {entry}", "Optimization of conditions", "formal (4+2) annulation / kinetic resolution")
    common_annulation_participants(r, rac1a, m2a, m3aa, s1a, evidence, ligand=ligand,
                                   solvent=solvent, cu_equiv=ox_eq, include_base=base, include_dmso=dmso)
    # Replace Cu oxidant in entries 13/14.
    if oxidant != "cu_oac2_h2o":
        p = next(p for p in r["participants"] if p["participant_id"] == "oxidant")
        p["molecule_id"] = molecules[f"{DOC}:mol_{slug(oxidant)}"]["molecule_id"]
    r["general_scheme_ids"] = [generic_id]
    r["procedure"]["modifications"] = F(deviation, deviation, evidence)
    outcome(r, "conversion", "conversion", [rac1a], evidence, f"C = {conv}%", qty(conv, "%"),
            basis="C = eeSM/(eeSM + eePR)", scope={"metric": "calculated conversion"})
    outcome(r, "product_ee", "selectivity", [m3aa], evidence, f"{ee_p}% ee", qty(ee_p, "%"),
            basis="chiral HPLC analysis of the reaction crude", scope={"metric": "ee", "entity": "3aa"})
    outcome(r, "starting_material_ee", "selectivity", [s1a], evidence, f"{ee_s}% ee", qty(ee_s, "%"),
            basis="chiral HPLC analysis of the reaction crude", scope={"metric": "ee", "entity": "recovered 1a"})
    outcome(r, "selectivity_factor", "selectivity", [m3aa, s1a], evidence, f"s = {sf}", sf,
            basis="s = ln[(1-C)(1-eeSM)]/ln[(1-C)(1+eeSM)]", scope={"metric": "kinetic-resolution selectivity factor s"})
    r["review"]["unresolved"] = ["No product/recovery yield is reported in Table 1; no SI was supplied for workup or isolation details."]
    reactions.append(r)


# Scheme 2 allene scope.
scheme2_rows = [
    ("3aa", "2a", 46, 94, 35, 90, 49, 100, 70, None, None),
    ("3ab", "2b", 40, 90, 36, 82, 48, 48, 70, None, None),
    ("3ac", "2c", 45, 92, 36, 90, 49, 74, 70, "Z:E > 20:1", None),
    ("3ad", "2d", 41, 90, 35, 96, 52, 74, 75, "Z:E > 10:1", None),
    ("3ae", "2e", 35, 80, 25, 88, 52, 26, 80, "Z:E > 20:1", None),
    ("3af", "2f", 30, 94, 32, 74, 44, 72, 75, "Z:E > 10:1", None),
    ("3ag", "2g", 40, 88, 38, 78, 47, 37, 80, "Z:E > 20:1", None),
    ("3ah", "2h", 31, 94, 36, 62, 40, 61, 80, "Z:E > 20:1", None),
    ("3ai", "2i", 35, 86, 39, 88, 51, 38, 70, "Z:E > 1:15", "4.5:1 rr"),
]
for prod, allene, py, pee, ry, ree, conv, sf, temp, ze, rr in scheme2_rows:
    evidence = [EV[2], EVT[3]] + ([EVT[4]] if prod in {"3ai"} else [])
    pmid = molecules[f"{DOC}:mol_{prod}"]["molecule_id"]
    amid = molecules[f"{DOC}:mol_{allene}"]["molecule_id"]
    r = base_reaction(f"scheme2_{prod}", evidence, "specific_experiment", ["substrate_scope"],
                      f"Scheme 2 panel {prod}", "Scope of the formal (4+2) cycloaddition", "formal (4+2) annulation / kinetic resolution")
    common_annulation_participants(r, rac1a, amid, pmid, s1a, evidence, temperature=temp)
    r["general_scheme_ids"] = [generic_id]
    outcome(r, "product_yield", "yield", [pmid], evidence, f"{py}% yield", qty(py, "%"), role="major product")
    outcome(r, "product_ee", "selectivity", [pmid], evidence, f"{pee}% ee", qty(pee, "%"),
            scope={"metric": "ee", "entity": prod})
    outcome(r, "recovered_starting_material", "recovery", [s1a], evidence, f"1a, {ry}% yield", qty(ry, "%"),
            role="recovered starting material", yield_type="recovery")
    outcome(r, "starting_material_ee", "selectivity", [s1a], evidence, f"{ree}% ee", qty(ree, "%"),
            scope={"metric": "ee", "entity": "recovered 1a"})
    outcome(r, "conversion", "conversion", [rac1a], evidence, f"C = {conv}%", qty(conv, "%"),
            basis="C reported in panel; formula defined in Table 1", scope={"metric": "conversion"})
    outcome(r, "selectivity_factor", "selectivity", [pmid, s1a], evidence, f"s = {sf}", sf,
            basis="selectivity-factor formula defined in Table 1", scope={"metric": "kinetic-resolution selectivity factor s"})
    if ze:
        qual, ratio = ("gt", ze.split(">", 1)[1].strip()) if ">" in ze else (None, ze.split()[-1])
        outcome(r, "z_e_ratio", "selectivity", [pmid], evidence, ze,
                {"ratio": ratio, "qualifier": qual} if qual else {"ratio": ratio}, scope={"metric": "Z:E", "order": "Z:E"})
    if rr:
        outcome(r, "regioisomeric_ratio", "selectivity", [pmid], evidence, rr, {"ratio": "4.5:1"},
                scope={"metric": "rr", "order": "major regioisomer:minor regioisomer"})
    r["review"]["unresolved"] = ["Scheme 2 does not state yield or selectivity measurement methods; complete structures/workup require the unavailable SI."]
    reactions.append(r)


# Prose-only failed free allenol experiment.
free_allenol = molecules[f"{DOC}:mol_free_allenol"]["molecule_id"]
r = base_reaction("prose_free_allenol_failure", [EVT[3]], "specific_experiment", ["substrate_scope", "control"],
                  "Page 3 prose below Scheme 2", "Attempted annulation with buta-2,3-dien-1-ol", "failed annulation")
r["participants"] = [participant("substrate", rac1a, "reactant", [EVT[3]]),
                     participant("allene", free_allenol, "reactant", [EVT[3]])]
reaction_step(r, [EVT[3]], None, None, None,
              operation="attempted formal annulation; exact main-text conditions not separately stated")
outcome(r, "complex_mixture", "decomposition", [], [EVT[3]], "failed to participate; complex mixture of products", "complex mixture")
no_reaction_classification(r, [EVT[3]])
r["review"]["unresolved"] = ["Exact temperature, charges, product identities and analytical method are not reported in the main-text sentence."]
reactions.append(r)


# Scheme 3 allyltriflamide scope.
scheme3_rows = [
    ("3ba", "1b", 46, 90, 26, 92, 51, 62, 75, "boc_l_phe_nhome", "2a"),
    ("3ca", "1c", 49, 94, 42, 94, 50, 115, 70, "boc_l_phe_nhome", "2a"),
    ("3da", "1d", 36, 92, 22, 92, 50, 79, 70, "boc_l_phe_nhome", "2a"),
    ("3ea", "1e", 45, 90, 42, 99, 52, 99, 70, "boc_l_phe_nhome", "2a"),
    ("3fa", "1f", 47, 94, 50, 96, 51, 127, 70, "boc_l_phe_nhome", "2a"),
    ("3ga", "1g", 33, 92, 41, 99, 52, 126, 70, "boc_l_phe_nhome", "2a"),
    ("3ha", "1h", 49, 92, 40, 90, 49, 74, 80, "boc_l_phe_nhome", "2a"),
    ("3ia", "1i", 48, 90, 40, 88, 49, 55, 80, "boc_l_phe_nhome", "2a"),
    ("3ja", "1j", 46, 90, 30, 76, 46, 44, 70, "boc_l_phe_nhome", "2a"),
    ("3ka", "1k", 37, 84, 31, 82, 49, 29, 40, "boc_phe_oh", "2a"),
    ("3la", "1l", 28, 86, 22, 52, 38, 22, 40, "boc_phe_oh", "2a"),
    ("3ma", "1m", 49, 90, 33, 96, 52, 74, 40, "boc_phe_oh", "2a"),
    ("3hc", "1h", 40, 90, 38, 76, 46, 44, 80, "boc_l_phe_nhome", "2c"),
]
for prod, sub, py, pee, ry, ree, conv, sf, temp, ligand, allene in scheme3_rows:
    evidence = [EV[3], EVT[4]]
    pmid = molecules[f"{DOC}:mol_{prod}"]["molecule_id"]
    smid = molecules[f"{DOC}:mol_rac_{sub}"]["molecule_id"]
    recmid = molecules[f"{DOC}:mol_s_{sub}"]["molecule_id"]
    amid = molecules[f"{DOC}:mol_{allene}"]["molecule_id"]
    r = base_reaction(f"scheme3_{prod}", evidence, "specific_experiment", ["substrate_scope"],
                      f"Scheme 3 panel {prod}", "Scope of the formal (4+2) cycloaddition", "formal (4+2) annulation / kinetic resolution")
    common_annulation_participants(r, smid, amid, pmid, recmid, evidence, ligand=ligand, temperature=temp)
    r["general_scheme_ids"] = [generic_id]
    outcome(r, "product_yield", "yield", [pmid], evidence, f"{py}% yield", qty(py, "%"), role="major product")
    outcome(r, "product_ee", "selectivity", [pmid], evidence, f"{pee}% ee", qty(pee, "%"), scope={"metric": "ee", "entity": prod})
    outcome(r, "recovered_starting_material", "recovery", [recmid], evidence, f"{sub}, {ry}% yield", qty(ry, "%"),
            role="recovered starting material", yield_type="recovery")
    outcome(r, "starting_material_ee", "selectivity", [recmid], evidence, f"{ree}% ee", qty(ree, "%"),
            scope={"metric": "ee", "entity": f"recovered {sub}"})
    outcome(r, "conversion", "conversion", [smid], evidence, f"C = {conv}%", qty(conv, "%"),
            basis="C reported in panel; formula defined in Table 1", scope={"metric": "conversion"})
    outcome(r, "selectivity_factor", "selectivity", [pmid, recmid], evidence, f"s = {sf}", sf,
            basis="selectivity-factor formula defined in Table 1", scope={"metric": "kinetic-resolution selectivity factor s"})
    if prod == "3hc":
        outcome(r, "z_e_ratio", "selectivity", [pmid], evidence, "Z:E > 20:1", {"ratio": "20:1", "qualifier": "gt"},
                scope={"metric": "Z:E", "order": "Z:E"})
    r["review"]["unresolved"] = ["Scheme 3 does not state yield/ee measurement methods; complete structures and workup require the unavailable SI."]
    reactions.append(r)


def simple_deriv(suffix, evidence, start_mid, product_mid, name, participants_extra,
                 temp, time, yield_value, ee_value, classification_changes, unresolved):
    r = base_reaction(suffix, evidence, "specific_experiment", ["derivatization"],
                      f"Scheme {'4' if suffix.startswith('scheme4') else '5'} {name}", name, name)
    r["participants"] = [participant("starting_material", start_mid, "reactant", evidence),
                         participant("product", product_mid, "product", evidence, product_role="reported product")]
    r["participants"].extend(participants_extra)
    reaction_step(r, evidence, temp, time, None, operation=name)
    if yield_value is not None:
        outcome(r, "product_yield", "yield", [product_mid], evidence, f"{yield_value}% yield", qty(yield_value, "%"), role="product")
    if ee_value is not None:
        outcome(r, "product_ee", "selectivity", [product_mid], evidence, f"{ee_value}% ee", qty(ee_value, "%"), scope={"metric": "ee"})
    classify_changes(r, classification_changes, evidence, unresolved)
    r["review"]["unresolved"] = [unresolved, "Yield/ee measurement bases and unshown workup are not reported in the main text."]
    reactions.append(r)
    return r


m3ac = molecules[f"{DOC}:mol_3ac"]["molecule_id"]
for k in ["4", "5", "6", "7", "8", "9", "s_3aa"]:
    assert f"{DOC}:mol_{slug(k)}" in molecules
m4, m5, m6, m7, m8, m9, ms3aa = [molecules[f"{DOC}:mol_{slug(k)}"]["molecule_id"] for k in ["4", "5", "6", "7", "8", "9", "s_3aa"]]

r4 = simple_deriv("scheme4_3ac_to_4", [EV[3], EVT[4]], m3ac, m4, "hydrogenation/isomerization to 4", [
    participant("hydrogen", molecules[f"{DOC}:mol_h2"]["molecule_id"], "reagent", [EV[3]]),
    participant("pd_c", molecules[f"{DOC}:mol_pd_c"]["molecule_id"], "catalyst", [EV[3]], molpct=qty(10, "mol%", "10 mol%"), reference=m3ac),
    participant("solvent", molecules[f"{DOC}:mol_etoac"]["molecule_id"], "solvent", [EV[3]]),
], ("rt", "rt"), ("24 h", qty(24, "h")), 91, 92, [
    {"atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
    {"atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 2},
], "Hydrogenation gives net C-H formation; the accompanying alkene isomerization cannot be completely atom-mapped from the main-text drawing.")
r4["solvent_conditions"]["composition"] = F("AcOEt", [molecules[f"{DOC}:mol_etoac"]["molecule_id"]], [EV[3]])
r4["notes"].append(F("starting material 3ac is shown as 92% ee", {"starting_material": "3ac", "ee": qty(92, "%")}, [EV[3]]))

r5 = simple_deriv("scheme4_3ac_to_5", [EV[3], EVT[4]], m3ac, m5, "RuCl3/NaIO4 oxidation to 5", [
    participant("rucl3", molecules[f"{DOC}:mol_rucl3"]["molecule_id"], "catalyst", [EV[3]], molpct=qty(15, "mol%", "15 mol%"), reference=m3ac),
    participant("naio4", molecules[f"{DOC}:mol_naio4"]["molecule_id"], "reagent", [EV[3]], equiv=qty(3.5, "equiv", "3.5 equiv"), reference=m3ac),
    participant("water", molecules[f"{DOC}:mol_water"]["molecule_id"], "solvent", [EV[3]]),
    participant("acn", molecules[f"{DOC}:mol_acn"]["molecule_id"], "solvent", [EV[3]]),
    participant("etoac", molecules[f"{DOC}:mol_etoac"]["molecule_id"], "solvent", [EV[3]]),
], ("rt", "rt"), ("0.5 h", qty(0.5, "h")), 66, 96, [
    {"atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}], "before_order": 0, "after_order": 2, "change_type": "formed"},
], "Source drawing supports carbonyl formation (DB-C-O); complete oxidative-cleavage atom correspondence is unresolved.")
r5["solvent_conditions"]["composition"] = F("H2O/CH3CN/AcOEt", [molecules[f"{DOC}:mol_water"]["molecule_id"], molecules[f"{DOC}:mol_acn"]["molecule_id"], molecules[f"{DOC}:mol_etoac"]["molecule_id"]], [EV[3]])
r5["solvent_conditions"]["mixture_ratio"] = F("H2O/CH3CN/AcOEt 0.75:1:1", "0.75:1:1", [EV[3]])
r5["notes"].append(F("starting material 3ac is shown as 92% ee", {"starting_material": "3ac", "ee": qty(92, "%")}, [EV[3]]))

r6 = simple_deriv("scheme4_3ac_to_6", [EV[3], EVT[4]], m3ac, m6, "Red-Al detriflylation to 6", [
    participant("red_al", molecules[f"{DOC}:mol_red_al"]["molecule_id"], "reagent", [EV[3]], equiv=qty(10, "equiv", "10 equiv"), reference=m3ac),
    participant("solvent", molecules[f"{DOC}:mol_toluene"]["molecule_id"], "solvent", [EV[3]]),
], ("50 °C", qty(50, "°C")), ("20 h", qty(20, "h")), 99, 92, [
    {"atoms": [{"atom_map": None, "element": "N"}, {"atom_map": None, "element": "S"}], "before_order": 1, "after_order": 0, "change_type": "broken"},
    {"atoms": [{"atom_map": None, "element": "N"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain"},
], "Detriflylation drawing supports N-S cleavage and net N-H formation; hydride/proton sources are not atom-mapped.")
r6["notes"].append(F("Scheme 4: 99% yield; page 4 prose: quantitative yield", {"scheme": "99%", "prose": "quantitative"}, [EV[3], EVT[4]]))
r6["solvent_conditions"]["composition"] = F("toluene", [molecules[f"{DOC}:mol_toluene"]["molecule_id"]], [EV[3]])
r6["notes"].append(F("starting material 3ac is shown as 92% ee", {"starting_material": "3ac", "ee": qty(92, "%")}, [EV[3]]))

r7 = simple_deriv("scheme5_s1a_to_7", [EV[3], EVT[4]], s1a, m7, "oxidation to keto amino product 7", [
    participant("rucl3", molecules[f"{DOC}:mol_rucl3"]["molecule_id"], "catalyst", [EV[3]], molpct=qty(3, "mol%", "3 mol%"), reference=s1a),
    participant("naio4", molecules[f"{DOC}:mol_naio4"]["molecule_id"], "reagent", [EV[3]]),
    participant("ccl4", molecules[f"{DOC}:mol_ccl4"]["molecule_id"], "solvent", [EV[3]]),
    participant("acn", molecules[f"{DOC}:mol_acn"]["molecule_id"], "solvent", [EV[3]]),
    participant("water", molecules[f"{DOC}:mol_water"]["molecule_id"], "solvent", [EV[3]]),
], ("rt", "rt"), None, 71, 99, [
    {"atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}], "before_order": 0, "after_order": 2, "change_type": "formed"},
], "Oxidative transformation drawing supports carbonyl formation; cleavage products and full atom mapping are not shown.")
r7["solvent_conditions"]["composition"] = F("CCl4/CH3CN/H2O", [molecules[f"{DOC}:mol_ccl4"]["molecule_id"], molecules[f"{DOC}:mol_acn"]["molecule_id"], molecules[f"{DOC}:mol_water"]["molecule_id"]], [EV[3]])
r7["solvent_conditions"]["mixture_ratio"] = F("CCl4/CH3CN/H2O 1:1:1.6", "1:1:1.6", [EV[3]])
r7["notes"].append(F("starting material (S)-1a is shown as 99% ee", {"starting_material": "(S)-1a", "ee": qty(99, "%")}, [EV[3]]))

r8 = simple_deriv("scheme5_s1a_to_8", [EV[3], EVT[4]], s1a, m8, "N-propargylation to enyne 8", [
    participant("base", molecules[f"{DOC}:mol_k2co3"]["molecule_id"], "reagent", [EV[3]]),
    participant("propargyl_bromide", molecules[f"{DOC}:mol_propargyl_bromide"]["molecule_id"], "reagent", [EV[3]]),
], ("60 °C", qty(60, "°C")), None, None, None, [
    {"atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "N"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
], "N-propargylation drawing supports new C-N formation; bromide and N-H coproduct changes are not mapped.")
r8["notes"].append(F("starting material (S)-1a is shown as 99% ee", {"starting_material": "(S)-1a", "ee": qty(99, "%")}, [EV[3]]))

r9 = simple_deriv("scheme5_8_to_9", [EV[3], EVT[4]], m8, m9, "Ir-catalyzed cycloisomerization to piperidine 9", [
    participant("ir_catalyst", molecules[f"{DOC}:mol_ir_catalyst"]["molecule_id"], "catalyst", [EV[3]], molpct=qty(7.5, "mol%", "7.5 mol%"), reference=m8),
    participant("solvent", molecules[f"{DOC}:mol_dce"]["molecule_id"], "solvent", [EV[3]]),
], ("reflux", "reflux"), None, 50, 99, [
    {"atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
    {"atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 1, "after_order": 2, "change_type": "order_increased"},
], "Cycloisomerization drawing supports new C-C ring closure and a new C=C bond; complete enyne atom mapping is unresolved.")
r9["solvent_conditions"]["composition"] = F("dichloroethane", [molecules[f"{DOC}:mol_dce"]["molecule_id"]], [EV[3]])
reactions[-1]["outcomes"].append(st.new_outcome("diastereoselectivity", "selectivity", [m9], [EVT[4]]))
reactions[-1]["outcomes"][-1]["value"] = F("single diastereoisomer", "single diastereoisomer", [EVT[4]])
reactions[-1]["outcomes"][-1]["scope"] = F("diastereoselectivity; isomer order not reported", {"metric": "dr", "reported": "single diastereoisomer", "isomer_order": None}, [EVT[4]])

# Enantioenriched starting material undergoes the annulation with the D ligand.
r = base_reaction("scheme5_s1a_2a_to_s3aa", [EV[3], EVT[4]], "specific_experiment", ["derivatization"],
                  "Scheme 5 arrow d", "Formal (4+2) annulation of (S)-1a using Boc-D-Phe-NHOMe", "formal (4+2) annulation")
common_annulation_participants(r, s1a, m2a, ms3aa, s1a, [EV[3], EVT[4]], ligand="boc_d_phe_nhome", temperature=70)
# Remove the generic recovered-material output: this experiment reports product only.
r["participants"] = [p for p in r["participants"] if p["participant_id"] != "recovered_starting_material"]
r["steps"][0]["participant_ids"] = [p["participant_id"] for p in r["participants"]]
outcome(r, "product_yield", "yield", [ms3aa], [EV[3], EVT[4]], "88% yield", qty(88, "%"), role="product")
outcome(r, "product_ee", "selectivity", [ms3aa], [EV[3], EVT[4]], "99% ee", qty(99, "%"),
        scope={"metric": "ee"})
r["review"]["unresolved"] = ["Yield/ee measurement bases and workup are not reported in the main text; conditions are inherited from Table 1 except the explicitly stated D ligand."]
reactions.append(r)


bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": f"{DOC}:scheme1c_challenges", "document_id": DOC, "evidence_ids": [EV[1]],
     "kind": "conceptual_illustration", "summary": "Scheme 1C illustrates challenges of allylamine C-H functionalization; it is not a discrete experiment."},
    {"context_id": f"{DOC}:scheme1d_natural_products", "document_id": DOC, "evidence_ids": [EV[1]],
     "kind": "contextual_structures", "summary": "Scheme 1D shows palustrine and lepadin B as motivation, not syntheses reported in this paper."},
    {"context_id": f"{DOC}:z_selectivity_explanation", "document_id": DOC, "evidence_ids": [EVT[4]],
     "kind": "author_interpretation", "summary": "Authors attribute preferred Z formation with monosubstituted allenes to steric effects."},
    {"context_id": f"{DOC}:si_notice", "document_id": DOC, "evidence_ids": [EV[4], EVT[4]],
     "kind": "source_availability", "summary": "Main text advertises SI, but no SI file was explicitly supplied or associated in this workspace."},
    {"context_id": f"{DOC}:pages5_6_references", "document_id": DOC, "evidence_ids": [EV[5], EV[6]],
     "kind": "non_experimental_pages", "summary": "Pages 5-6 contain author/acknowledgment/reference material only; cited reactions were not imported as current-paper experiments."},
]

dump(RUN / "submitted.json", bundle)

# Save raw visual transcription separately from the formal submission.
raw = {
    "document_id": DOC,
    "pages_reviewed": [1, 2, 3, 4, 5, 6],
    "source": "assistant visual inspection of supplied full-page JPEGs plus extracted page text",
    "table1_rows": [
        {"entry": x[0], "deviation": x[1], "conversion_percent": x[8], "product_ee_percent": x[9],
         "starting_material_ee_percent": x[10], "selectivity_factor": x[11]} for x in table_rows
    ],
    "scheme2_panels": [dict(product=x[0], allene=x[1], product_yield=x[2], product_ee=x[3],
                                   recovered_1a_yield=x[4], recovered_1a_ee=x[5], conversion=x[6],
                                   selectivity_factor=x[7], temperature_C=x[8], z_e=x[9], rr=x[10]) for x in scheme2_rows],
    "scheme3_panels": [dict(product=x[0], substrate=x[1], product_yield=x[2], product_ee=x[3],
                                   recovered_yield=x[4], recovered_ee=x[5], conversion=x[6],
                                   selectivity_factor=x[7], temperature_C=x[8], ligand=x[9], allene=x[10]) for x in scheme3_rows],
    "derivatization_arrows": ["3ac->4", "3ac->5", "3ac->6", "(S)-1a->7", "(S)-1a->8", "8->9", "(S)-1a+2a->(S)-3aa"],
    "negative_scope": "buta-2,3-dien-1-ol failed; complex mixture",
    "non_experimental": {"page5": "author information, acknowledgments, references", "page6": "references"},
    "limitations": ["No SI supplied", "Dense structures were not assigned unsafe complete SMILES", "All formal records remain needs_review"],
}
dump(RUN / "raw_model/visual_transcription.json", raw)

diagnostics = {}
for mid, m in molecules.items():
    smi = m["smiles"]["normalized_value"]
    if smi:
        diagnostics[mid] = st.stereo_inspect(smi)
dump(RUN / "raw_model/stereo_diagnostics.json", diagnostics)
dump(RUN / "raw_model/source_mentions.json",
     st.source_mentions(dataset, ["1a", "3ac", "3ba", "3hc", "product 4", "product 7", "piperidine 9"], context=140))

# Exact per-page coverage requested by the batch controller.
ids = [r["reaction_id"] for r in reactions]
def matching(*needles):
    return [rid for rid in ids if any(n in rid for n in needles)]

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title/abstract", "finding": "Article identity and current formal (4+2) annulation/kinetic-resolution claim.", "evidence_ids": [EV[1], EVT[1]]},
            {"region": "Scheme 1A", "finding": "Prior literature context; classified as literature_example, not a current experiment.", "evidence_ids": [EV[1]]},
            {"region": "Scheme 1B", "finding": "Generic current-work transformation; classified as generic_scheme.", "evidence_ids": [EV[1]]},
            {"region": "Scheme 1C-D", "finding": "Conceptual challenge illustration and natural-product motivation; no current synthesis experiment.", "evidence_ids": [EV[1]]},
        ], "record_ids": matching("scheme1a_", "scheme1b_"), "unresolved": ["Generic R groups are not expanded; no SI supplied."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Table 1 entries 1-20 and footnotes", "finding": "Twenty independent optimization experiments transcribed row-by-row with conversion, both ee values and s.", "evidence_ids": [EV[2], EVT[2]]},
            {"region": "Scheme 2 panels 3aa-3ai and footnotes", "finding": "Nine allene-scope experiments transcribed with product/recovery values and panel selectivities.", "evidence_ids": [EV[2]]},
            {"region": "model-condition prose", "finding": "Duplicates Table 1 entry 1; merged as supporting evidence rather than a second experiment.", "evidence_ids": [EVT[2]]},
        ], "record_ids": matching("table1_entry_", "scheme2_3"), "unresolved": ["Scope yield and selectivity measurement methods are not printed; SI unavailable."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Scheme 3 panels 3ba-3hc and footnotes", "finding": "Thirteen allyltriflamide-scope experiments transcribed independently.", "evidence_ids": [EV[3]]},
            {"region": "Scheme 4 arrows a-c", "finding": "Three distinct derivatization experiments from 3ac.", "evidence_ids": [EV[3]]},
            {"region": "Scheme 5 arrows a-d", "finding": "Four distinct experiments: oxidation, N-propargylation, cycloisomerization, and enantioenriched annulation.", "evidence_ids": [EV[3]]},
            {"region": "lower-column scope prose", "finding": "Buta-2,3-dien-1-ol negative/complex-mixture experiment and narrative corroboration of Scheme 2.", "evidence_ids": [EVT[3]]},
        ], "record_ids": matching("scheme3_", "scheme4_", "scheme5_", "free_allenol"), "unresolved": ["Dense structures and methods not printed in main text remain unresolved without SI."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "left-column scope discussion", "finding": "Narrative corroboration of Scheme 2/3; X-ray assignments (R)-3ba and (S)-1b; no duplicate experiments created.", "evidence_ids": [EVT[4]]},
            {"region": "derivatization discussion", "finding": "Narrative corroboration of Schemes 4/5, including quantitative prose for the 99% Red-Al result and single-diastereoisomer statement for 9.", "evidence_ids": [EVT[4]]},
            {"region": "Associated Content", "finding": "SI is advertised but was not supplied/associated; no SI facts imported.", "evidence_ids": [EV[4], EVT[4]]},
        ], "record_ids": matching("scheme2_", "scheme3_", "scheme4_", "scheme5_"), "unresolved": ["Main-text qualitative 'quantitative yield' and scheme 99% are retained together; no attempt made to inspect absent SI."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "author information/acknowledgments/references", "finding": "No current-paper experiment, reaction arrow, table row, mass-spectrum table, molecular-weight table, or fragment-abundance table.", "evidence_ids": [EV[5], EVT[5]]},
        ], "record_ids": [], "unresolved": []},
        {"page": 6, "state": "attempted", "regions": [
            {"region": "reference list continuation", "finding": "Cited literature only; no current-paper experiment or independent characterization table.", "evidence_ids": [EV[6], EVT[6]]},
        ], "record_ids": [], "unresolved": []},
    ],
}
dump(RUN / "coverage/pages_0001_0006.json", coverage)

# When a merged dataset is already present, persist its deterministic advisory
# audit as a separate artifact.  It never promotes review status.
if dataset.get("experiments"):
    from reaction_tools import dataset_audit
    dump(RUN / "audit.json", dataset_audit(dataset))

print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
                  "pages": len(coverage["pages"])}, ensure_ascii=False))
