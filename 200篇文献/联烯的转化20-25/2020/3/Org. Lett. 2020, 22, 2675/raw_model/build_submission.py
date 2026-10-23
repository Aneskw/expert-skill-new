#!/usr/bin/env python3
"""Build the evidence-linked v3.1 submission for the six-page main article.

This script is a reproducible transcription of the page-image review.  It does
not modify the extraction tools, schemas, or scientific validation rules.
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
DATASET_PATH = RUN / "dataset.json"
DOC = "doc_f3ddb3d9d361e80d"
EV = {page: f"{DOC}:page{page}:image" for page in range(1, 7)}
TEXT = {page: f"{DOC}:page{page}:text" for page in range(1, 7)}
TITLE = "Phosphine-Catalyzed Cascade Michael Addition/[4+2] Cycloaddition Reaction of Allenoates and 2-Arylidene-1,3-indanediones"


def dump(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def F(raw, normalized=None, evidence=(), method="original_extraction", derivation=None):
    if normalized is None and raw is not None:
        normalized = raw
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation)


def D(normalized, evidence, derivation):
    return st.fact(None, normalized, evidence, method="structure_derivation", derivation=derivation)


def I(normalized, evidence, derivation):
    return st.fact(None, normalized, evidence, method="inference", derivation=derivation)


def missing(evidence=()):
    return st.fact(evidence=evidence)


def set_status(m, key, normalized, evidence, raw=None, method="inference", derivation=None):
    if method == "original_extraction":
        m["stereochemistry"][key] = F(raw or normalized, normalized, evidence)
    else:
        m["stereochemistry"][key] = st.fact(
            raw if method == "structure_derivation" else None,
            normalized,
            evidence,
            method=method,
            derivation=derivation or "Status assigned from the displayed structure and current representation limits.",
        )


dataset = json.loads(DATASET_PATH.read_text())
doc = dataset["documents"][0]
assert doc["document_id"] == DOC and doc["pages"] == 6
doc["title"] = TITLE
doc["supplementary_material"] = {
    "status": "unknown",
    "note": "The main article states that SI is available, but no SI file was explicitly supplied or associated with this run.",
}
DATASET_PATH.write_text(json.dumps(dataset, ensure_ascii=False, indent=2) + "\n")
bundle = st.empty(dataset["documents"])

molecules = {}
stereo_diagnostics = {}


def add_molecule(mid, label, name, evidence, smiles=None, stereo=None, note=None):
    m = st.new_molecule(mid, DOC, evidence)
    m["occurrences"][0]["compound_label"] = F(label, label, evidence)
    m["names"]["original"] = F(name, name, evidence)
    m["names"]["english"] = F(name, name, evidence)
    m["identity_note"] = F(note, note, evidence) if note else missing(evidence)
    if smiles:
        m["identity_status"] = "source_identified"
        m["representation_type"] = "smiles"
        m["smiles"] = st.fact(
            f"Candidate connectivity transcribed from the displayed structure for {label}",
            smiles,
            evidence,
            method="structure_derivation",
            derivation="Candidate SMILES transcribed from the page drawing/name; parsing does not establish source correctness.",
        )
        diag = st.stereo_inspect(smiles)
        stereo_diagnostics[mid] = diag
        atom_potential = [p for p in diag["potential_stereo"] if str(p["type"]).startswith("Atom_")]
        bond_potential = [p for p in diag["potential_stereo"] if str(p["type"]).startswith("Bond_")]
        set_status(
            m,
            "chirality_status",
            "unspecified" if atom_potential else "achiral",
            evidence,
            method="structure_derivation",
            raw=f"Displayed structure/name for {label}",
            derivation="RDKit graph diagnostics on the candidate transcription; absence of a toolkit flag is not used for complex unresolved structures.",
        )
        set_status(
            m,
            "geometry_status",
            "unspecified" if bond_potential else "not_applicable",
            evidence,
            method="structure_derivation",
            raw=f"Displayed structure/name for {label}",
            derivation="Potential double-bond stereochemistry assessed only on the candidate graph; no unprinted E/Z assignment was added.",
        )
        set_status(
            m,
            "representation_status",
            "partial" if atom_potential or bond_potential else "encoded",
            evidence,
            method="structure_derivation",
            raw=f"Candidate SMILES for {label}",
            derivation="Representation status reflects whether the candidate graph contains unspecified potential stereochemistry.",
        )
        m["stereochemistry"]["atom_centers"] = D(diag["atom_centers"], evidence, "RDKit diagnostics on the candidate SMILES; not source verification.")
        m["stereochemistry"]["double_bonds"] = D(diag["double_bonds"], evidence, "RDKit diagnostics on the candidate SMILES; carbonyl bonds are retained as graph output and are not E/Z claims.")
        m["stereochemistry"]["isomeric_smiles"] = D(diag["canonical_isomeric_smiles"], evidence, "Canonical isomeric SMILES generated from the candidate graph.")
        m["stereochemistry"]["identity_key"] = D(diag["stereo_graph_key"], evidence, "Graph key generated from the candidate isomeric SMILES; not a sample identity.")
        m["stereochemistry"]["axial"] = missing(evidence)
        m["stereochemistry"]["mixture"] = missing(evidence)
        m["stereochemistry"]["assignment"] = missing(evidence)
    else:
        m["identity_status"] = "source_identified"
        m["representation_type"] = "unresolved"
        set_status(m, "chirality_status", (stereo or {}).get("chirality", "unresolved"), evidence,
                   raw=(stereo or {}).get("chirality_raw"),
                   method="original_extraction" if (stereo or {}).get("chirality_raw") else "inference")
        set_status(m, "geometry_status", (stereo or {}).get("geometry", "unresolved"), evidence,
                   raw=(stereo or {}).get("geometry_raw"),
                   method="original_extraction" if (stereo or {}).get("geometry_raw") else "inference")
        set_status(m, "representation_status", (stereo or {}).get("representation", "not_encoded"), evidence,
                   method="inference", derivation="No safe complete SMILES was transcribed; source stereochemical markings are retained separately.")
        if (stereo or {}).get("assignment"):
            m["stereochemistry"]["assignment"] = F((stereo or {})["assignment"], (stereo or {})["assignment"], evidence)
        else:
            m["stereochemistry"]["assignment"] = missing(evidence)
        m["stereochemistry"]["atom_centers"] = missing(evidence)
        m["stereochemistry"]["double_bonds"] = missing(evidence)
        m["stereochemistry"]["axial"] = missing(evidence)
        m["stereochemistry"]["mixture"] = missing(evidence)
        m["stereochemistry"]["isomeric_smiles"] = missing(evidence)
        m["stereochemistry"]["identity_key"] = missing(evidence)
    molecules[mid] = m
    bundle["molecules"].append(m)
    return m


def add_generic(mid, label, name, evidence, fixed, variables, note):
    m = add_molecule(mid, label, name, evidence, note=note)
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": F(fixed, fixed, evidence),
        "variables": [],
        "instance_mapping": [],
        "reason": F(note, note, evidence),
    }
    for v in variables:
        m["generic_template"]["variables"].append({
            "label": v,
            "attachment_points": F(f"Attachment point shown for {v} in the page drawing", f"source-drawn attachment for {v}", evidence),
            "definition": F(f"{v} is a variable substituent in the source scheme", f"source-defined variable {v}", evidence),
            "allowed_values": missing(evidence),
            "constraints": missing(evidence),
        })
    return m


# Core source molecules and easy-to-check connectivity candidates.
add_molecule("m_1a", "1a", "tetrahydrobenzofuranone-derived allenoate 1a", [EV[1], EV[2], EV[3]],
             stereo={"chirality": "unspecified", "geometry": "unspecified", "representation": "not_encoded"},
             note="Complex allenoate drawing is legible, but a complete candidate SMILES was not forced from the compact main-text figure; no SI was provided.")
for suffix, alkyl in [("1b", "Et"), ("1c", "n-Pr"), ("1d", "i-Pr")]:
    add_molecule(f"m_{suffix}", suffix, f"alkyl-substituted allenoate {suffix} (R={alkyl})", [EV[3]],
                 stereo={"chirality": "unspecified", "geometry": "unspecified", "representation": "not_encoded"},
                 note=f"Table 3 defines R={alkyl}; the compact allenoate structure was not forced into a complete SMILES.")

aryl_smiles = {
    "2a": "O=C1C(=Cc2ccccc2)C(=O)c2ccccc21",
    "2b": "O=C1C(=Cc2ccccc2F)C(=O)c2ccccc21",
    "2c": "O=C1C(=Cc2cccc(F)c2)C(=O)c2ccccc21",
    "2d": "O=C1C(=Cc2ccc(F)cc2)C(=O)c2ccccc21",
    "2e": "O=C1C(=Cc2ccc(Cl)cc2)C(=O)c2ccccc21",
    "2f": "O=C1C(=Cc2ccc(Br)cc2)C(=O)c2ccccc21",
    "2g": "O=C1C(=Cc2ccc(I)cc2)C(=O)c2ccccc21",
    "2h": "O=C1C(=Cc2ccc(C(F)(F)F)cc2)C(=O)c2ccccc21",
    "2i": "O=C1C(=Cc2ccc([N+](=O)[O-])cc2)C(=O)c2ccccc21",
    "2j": "O=C1C(=Cc2c(Cl)cc(Cl)cc2)C(=O)c2ccccc21",
    "2k": "O=C1C(=Cc2ccccc2C)C(=O)c2ccccc21",
    "2l": "O=C1C(=Cc2cccc(C)c2)C(=O)c2ccccc21",
    "2m": "O=C1C(=Cc2ccc(C)cc2)C(=O)c2ccccc21",
    "2n": "O=C1C(=Cc2ccc(CC)cc2)C(=O)c2ccccc21",
    "2o": "O=C1C(=Cc2cccc(OC)c2)C(=O)c2ccccc21",
    "2p": "O=C1C(=Cc2ccc(OC)cc2)C(=O)c2ccccc21",
    "2q": "O=C1C(=Cc2c(C)ccc(C)c2)C(=O)c2ccccc21",
    "2r": "O=C1C(=Cc2cccc3ccccc23)C(=O)c2ccccc21",
    "2s": "O=C1C(=CC=Cc2ccccc2)C(=O)c2ccccc21",
}
names_2 = {
    "2a": "2-benzylidene-1H-indene-1,3(2H)-dione",
    "2b": "2-(2-fluorobenzylidene)-1,3-indanedione", "2c": "2-(3-fluorobenzylidene)-1,3-indanedione",
    "2d": "2-(4-fluorobenzylidene)-1,3-indanedione", "2e": "2-(4-chlorobenzylidene)-1,3-indanedione",
    "2f": "2-(4-bromobenzylidene)-1,3-indanedione", "2g": "2-(4-iodobenzylidene)-1,3-indanedione",
    "2h": "2-(4-trifluoromethylbenzylidene)-1,3-indanedione", "2i": "2-(4-nitrobenzylidene)-1,3-indanedione",
    "2j": "2-(2,4-dichlorobenzylidene)-1,3-indanedione", "2k": "2-(2-methylbenzylidene)-1,3-indanedione",
    "2l": "2-(3-methylbenzylidene)-1,3-indanedione", "2m": "2-(4-methylbenzylidene)-1,3-indanedione",
    "2n": "2-(4-ethylbenzylidene)-1,3-indanedione", "2o": "2-(3-methoxybenzylidene)-1,3-indanedione",
    "2p": "2-(4-methoxybenzylidene)-1,3-indanedione", "2q": "2-(2,5-dimethylbenzylidene)-1,3-indanedione",
    "2r": "2-(1-naphthylvinyl)-1,3-indanedione", "2s": "styryl-substituted 2-arylidene-1,3-indanedione",
}
for key, smi in aryl_smiles.items():
    ev = [EV[2], EV[3]] if key != "2a" else [EV[2], EV[3]]
    add_molecule(f"m_{key}", key, names_2[key], ev, smiles=smi,
                 note="Candidate connectivity follows the displayed Table 1/2 scaffold and printed R assignment; unprinted alkene configuration was not added.")

# Catalysts and solvents. Dense chiral phosphines P1-P3 are retained without unsafe SMILES.
simple = {
    "ph3p": ("Ph3P", "triphenylphosphine", "P(c1ccccc1)(c1ccccc1)c1ccccc1"),
    "ph2pme": ("Ph2PMe", "methyldiphenylphosphine", "CP(c1ccccc1)c1ccccc1"),
    "phpme2": ("PhPMe2", "dimethylphenylphosphine", "CP(C)c1ccccc1"),
    "nbu3p": ("n-Bu3P", "tributylphosphine", "CCCCP(CCCC)CCCC"),
    "ptol3p": ("(4-MeC6H4)3P", "tris(4-methylphenyl)phosphine", "P(c1ccc(C)cc1)(c1ccc(C)cc1)c1ccc(C)cc1"),
    "pf3p": ("(4-FC6H4)3P", "tris(4-fluorophenyl)phosphine", "P(c1ccc(F)cc1)(c1ccc(F)cc1)c1ccc(F)cc1"),
    "ph2pbn": ("Ph2PBn", "benzyldiphenylphosphine", "P(c1ccccc1)(c1ccccc1)Cc1ccccc1"),
    "dmap": ("DMAP", "4-dimethylaminopyridine", "CN(C)c1ccncc1"),
    "dabco": ("DABCO", "1,4-diazabicyclo[2.2.2]octane", "C1CN2CCN1CC2"),
    "dbu": ("DBU", "1,8-diazabicyclo[5.4.0]undec-7-ene", "C1CCC2=NCCCN2CC1"),
    "dcm": ("CH2Cl2", "dichloromethane", "ClCCl"),
    "dce": ("DCE", "1,2-dichloroethane", "ClCCCl"),
    "toluene": ("toluene", "toluene", "Cc1ccccc1"),
    "mecn": ("CH3CN", "acetonitrile", "CC#N"),
}
for key, (label, name, smi) in simple.items():
    add_molecule(f"m_{key}", label, name, [EV[2], EV[3]], smiles=smi)
for p in ("P1", "P2", "P3"):
    add_molecule(f"m_{p.lower()}", p, f"chiral phosphine {p}", [EV[3]],
                 stereo={"chirality": "specified", "chirality_raw": "stereochemical wedges/dashes shown in Scheme 2 catalyst drawing", "geometry": "unresolved", "representation": "not_encoded",
                         "assignment": "Scheme 2 depicts a single chiral catalyst structure with wedge/dash stereochemistry; no absolute descriptor is printed on the supplied page."},
                 note="Dense chiral catalyst drawing retained without a forced SMILES; absolute descriptor is not printed on the main-text page.")

add_molecule("m_allyl_reagent", "ethyl 2-(acetoxymethyl)acrylate", "ethyl 2-(acetoxymethyl)acrylate", [EV[3]],
             smiles="CCOC(=O)C(=C)COC(C)=O")


def major_product(label, evidence, xray=False):
    assignment = "Relative configuration shown by wedge/dash drawing"
    if xray:
        assignment += "; relative configuration reported as confirmed by X-ray crystallography"
    return add_molecule(
        f"m_{label}", label, f"major diastereomer of cascade product {label}", evidence,
        stereo={"chirality": "specified", "chirality_raw": "relative stereochemistry shown by wedge/dash bonds",
                "geometry": "unresolved", "representation": "not_encoded", "assignment": assignment},
        note="No absolute R/S descriptors or safely transcribed isomeric SMILES are present on the supplied main-text pages.",
    )


major_product("3aa", [EV[2], EV[3], EV[6]], xray=True)
add_molecule("m_4aa", "4aa", "[4+2] cycloadduct 4aa", [EV[2], EV[4], EV[6]],
             stereo={"chirality": "specified", "chirality_raw": "relative stereochemistry shown by wedge/dash bonds",
                     "geometry": "unresolved", "representation": "not_encoded",
                     "assignment": "Relative configuration shown; reported as confirmed by X-ray crystallography (CCDC 1921150)."},
             note="No absolute R/S descriptors or safely transcribed isomeric SMILES are present on the supplied main-text pages.")

table2 = [
    (1, "2b", "3ab", 1.5, 95, "8:1"), (2, "2c", "3ac", 1, 88, "15:1"),
    (3, "2d", "3ad", 1, 80, "15:1"), (4, "2e", "3ae", 1, 84, "14:1"),
    (5, "2f", "3af", 1, 87, "13:1"), (6, "2g", "3ag", 2, 84, "12:1"),
    (7, "2h", "3ah", 1, 83, ">19:1"), (8, "2i", "3ai", 3, 72, ">19:1"),
    (9, "2j", "3aj", 2, 59, "15:1"), (10, "2k", "3ak", 1, 78, "8:1"),
    (11, "2l", "3al", 1, 83, "11:1"), (12, "2m", "3am", 1, 78, "12:1"),
    (13, "2n", "3an", 2, 90, "15:1"), (14, "2o", "3ao", 2, 81, "9:1"),
    (15, "2p", "3ap", 1.5, 48, "10:1"), (16, "2q", "3aq", 2, 56, "12:1"),
    (17, "2r", "3ar", 11, 69, "10:1"), (18, "2s", "3as", 24, "trace", None),
]
for _, substrate, product, *_ in table2:
    if product == "3as":
        add_molecule("m_3as", "3as", "desired cascade product 3as", [EV[3]],
                     stereo={"chirality": "unresolved", "geometry": "unresolved", "representation": "not_encoded"},
                     note="Only a trace is reported and no dr or stereoisomer assignment is printed; a major diastereomer was not invented.")
    else:
        major_product(product, [EV[3]])
    add_molecule(f"m_4_from_{substrate}", f"4 (from {substrate})", f"corresponding trace [4+2] cycloadduct from {substrate}", [EV[3]],
                 stereo={"chirality": "unspecified", "geometry": "unresolved", "representation": "not_encoded"},
                 note="Table 2 calls this only product 4; an exact compound label and stereochemical assignment are not printed.")
for product in ("3ba", "3ca", "3da"):
    major_product(product, [EV[3]])
for substrate in ("1b", "1c", "1d"):
    add_molecule(f"m_4_from_{substrate}", f"4 (from {substrate})", f"corresponding trace [4+2] cycloadduct from {substrate}", [EV[3]],
                 stereo={"chirality": "unspecified", "geometry": "unresolved", "representation": "not_encoded"},
                 note="Table 3 calls this only product 4; an exact compound label and stereochemical assignment are not printed.")
add_molecule("m_5a", "5", "allylation adduct 5a", [EV[3]],
             stereo={"chirality": "specified", "chirality_raw": "relative stereochemistry shown by wedge/dash bonds",
                     "geometry": "unresolved", "representation": "not_encoded",
                     "assignment": "Relative configuration shown in Scheme 3b; no absolute descriptors printed."})

# Generic templates for page-1 schemes and the generic phosphine in Scheme 4.
add_generic("m_gen_allenoate", "generic allenoate 1", "tetrahydrobenzofuranone-derived allenoate 1", [EV[1]],
            "Source-drawn tetrahydrobenzofuranone-derived allenoate skeleton", ["R1"], "R1 is variable in the source scheme; no unique molecule exists.")
add_generic("m_gen_arylidene", "generic arylidene indanedione 2", "2-arylidene-1,3-indanedione 2", [EV[1], EV[3]],
            "Source-drawn 2-arylidene-1,3-indanedione skeleton", ["R2"], "R2 is variable in the source scheme; no unique molecule exists.")
add_generic("m_gen_product3", "generic cascade product 3", "spirocyclic 1,3-indanedione derivative 3", [EV[1], EV[3]],
            "Source-drawn cascade-product skeleton with relative wedge/dash bonds", ["R1", "R2"], "Generic product represents a family, not a concrete stereoisomeric molecule.")
add_generic("m_lit_lu_allenoate", "Lu allenoate", "generic allenoate used in Lu's prior work", [EV[1]],
            "Source-drawn allenoate with n-Bu and labeled alpha/beta/gamma/gamma-prime positions", [], "Prior-literature structure retained as a generic page drawing.")
add_generic("m_lit_lu_imine", "Lu imine", "generic alpha,beta-unsaturated imine", [EV[1]],
            "Source-drawn alpha,beta-unsaturated imine bearing TsN and CO2Me", ["Ph shown in example"], "Prior-literature example; compact structure not forced into a concrete registry identity.")
add_generic("m_lit_lu_product", "Lu product", "tetrahydropyridine product from prior work", [EV[1]],
            "Source-drawn tetrahydropyridine product with relative wedge/dash bonds", [], "Prior-literature example with no numerical result printed in Scheme 1a.")
add_generic("m_lit_guo_allenoate", "previous-work allenoate", "generic allenoate used in the authors' prior work", [EV[1]],
            "Source-drawn allenoate with labeled beta-prime and gamma-prime positions", [], "Prior-literature structure retained as a generic page drawing.")
add_generic("m_lit_cyclic_imine", "cyclic imine", "sulfamate-derived cyclic imine", [EV[1]],
            "Source-drawn sulfamate-derived cyclic imine", ["R"], "R is variable in the prior-work scheme.")
add_generic("m_lit_guo_product", "prior-work product", "annulation product from sulfamate-derived cyclic imine", [EV[1]],
            "Source-drawn fused annulation product with relative wedge/dash bonds", ["R"], "Prior-literature example retained as a generic structure.")
add_generic("m_pr3", "PR3", "generic phosphine catalyst", [EV[1], EV[4]], "Generic trivalent phosphine PR3", ["R"], "Mechanistic catalyst is generic PR3, not a unique molecule.")
for letter in "ABCDEFGHIJ":
    add_molecule(f"m_mech_{letter}", letter, f"proposed mechanistic intermediate {letter}", [EV[4]],
                 stereo={"chirality": "unresolved", "geometry": "unresolved", "representation": "not_encoded"},
                 note="Proposed intermediate in Scheme 4; not an isolated or experimentally established compound in the supplied pages.")


def charge(raw, number, unit, evidence, inherited=False):
    method = "inherited_general_condition" if inherited else "original_extraction"
    return F(raw, {"number": number, "unit": unit}, evidence, method=method)


def participant(r, pid, mid, role, evidence, charges=None, eqref=None, product_role=None):
    p = st.new_participant(pid, mid, role, evidence)
    if product_role:
        p["product_role"] = F(product_role, product_role, evidence)
    for key, val in (charges or {}).items():
        p["charges"][key] = val
    if eqref:
        p["equivalence_reference"] = F(eqref, eqref, evidence)
    r["participants"].append(p)
    return p


def cc_classification(r, evidence, count=3, limitation=None):
    changes = []
    ids = []
    for i in range(count):
        cid = f"b{i+1}"
        ids.append(cid)
        change = {"change_id": cid, "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
                  "before_order": 0, "after_order": 1, "change_type": "formed"}
        changes.append(D(change, evidence, "New C-C sigma connection inferred from the source-drawn reactant/product connectivity; atom mapping is unresolved."))
    r["classification"]["bond_changes"] = changes
    r["classification"]["labels"] = D([{"category_number": 2, "code": "C-C", "supporting_change_ids": ids}], evidence,
                                              "Taxonomy label derived from the displayed new carbon-carbon connections.")
    r["classification"]["coverage"] = I("partial", evidence, "No complete atom mapping or hydrogen-source accounting was asserted.")
    r["classification"]["limitations"] = [
        limitation or "Compact drawings support C-C bond formation, but exact atom mapping, bond-order shifts, and implicit-H changes remain unresolved.",
        "No three-model votes were produced.",
    ]


def no_transform_classification(r, evidence):
    r["classification"]["labels"] = I([], evidence, "No observed transformation; intended product was not used to infer a bond class.")
    r["classification"]["coverage"] = I("not_applicable", evidence, "The source reports NR/no observed product.")
    r["classification"]["limitations"] = ["Negative result; no intended bond formation is asserted."]


def unresolved_classification(r, evidence, why):
    r["classification"]["labels"] = I([], evidence, why)
    r["classification"]["coverage"] = I("unresolved", evidence, why)
    r["classification"]["limitations"] = [why]


def make_reaction(rid, evidence, record_type="specific_experiment", categories=None, location=None,
                  name="Phosphine-catalyzed cascade Michael addition/[4+2] cycloaddition"):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = F(categories or [], categories or [], evidence) if record_type == "specific_experiment" else F([], [], evidence)
    r["location"] = F(location or rid, location or rid, evidence)
    r["reaction_name"] = F(name, name, evidence)
    r["transformation_type"] = I("cascade Michael addition/[4+2] cycloaddition", evidence,
                                  "Normalized from the article's reported reaction name.")
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
                   "unresolved": ["No explicitly supplied SI; workup, purification, characterization, and some measurement details remain unavailable."]}
    r["completeness"] = "partial"
    bundle["reactions"].append(r)
    return r


def set_solvent(r, mid, label, evidence, volume_ml=None, inherited=False):
    charges = {}
    if volume_ml is not None:
        charges["volume"] = charge(f"{volume_ml:g} mL", volume_ml, "mL", evidence, inherited)
    participant(r, f"p_{label}", mid, "solvent", evidence, charges=charges)
    method = "inherited_general_condition" if inherited else "original_extraction"
    r["solvent_conditions"]["composition"] = F(label, [label], evidence, method=method)
    if volume_ml is not None:
        r["solvent_conditions"]["total_volume"] = F(f"{volume_ml:g} mL", {"number": volume_ml, "unit": "mL"}, evidence, method=method)


def add_step(r, evidence, temp=None, time=None, inherited=False, operation="reaction"):
    s = st.new_step(f"step_{len(r['steps'])+1}", evidence)
    s["participant_ids"] = [p["participant_id"] for p in r["participants"] if p["role"] != "product"]
    method = "inherited_general_condition" if inherited else "original_extraction"
    s["operation"] = F(operation, operation, evidence, method=method)
    if temp is not None:
        if isinstance(temp, (int, float)):
            s["temperature"] = F(f"{temp:g} °C", {"number": temp, "unit": "°C"}, evidence, method=method)
        else:
            s["temperature"] = F(temp, temp, evidence, method=method)
    if time is not None:
        s["time"] = F(f"{time:g} h", {"number": time, "unit": "h"}, evidence, method=method)
    r["steps"].append(s)


def outcome(r, oid, kind, targets, raw, normalized, evidence, basis=None, yield_type=None, scope=None, inherited=False):
    o = st.new_outcome(oid, kind, targets, evidence)
    method = "inherited_general_condition" if inherited else "original_extraction"
    o["value"] = F(raw, normalized, evidence, method=method)
    if yield_type:
        o["yield_type"] = F(yield_type, yield_type, evidence, method=method)
    if basis:
        o["measurement_basis"] = F(basis, basis, evidence, method=method)
    if scope:
        o["scope"] = F(scope if isinstance(scope, str) else json.dumps(scope, ensure_ascii=False), scope, evidence, method=method)
    r["outcomes"].append(o)
    return o


def add_dr(r, product_mid, dr, evidence, basis=None, inherited=False):
    outcome(r, f"o_dr_{len(r['outcomes'])+1}", "selectivity", [product_mid], dr, dr, evidence, basis=basis,
            scope={"metric": "dr", "order": "major:minor; stereoisomer identities unassigned"}, inherited=inherited)


# Scheme 1: two prior-work examples and one current-paper generic scheme.
r = make_reaction("scheme1a_lu_prior_work", [EV[1]], record_type="literature_example", location="Scheme 1a",
                  name="Lu group's phosphine-catalyzed [4+2] annulation")
participant(r, "p_lu_allenoate", "m_lit_lu_allenoate", "reactant", [EV[1]])
participant(r, "p_lu_imine", "m_lit_lu_imine", "reactant", [EV[1]])
participant(r, "p_pr3", "m_pr3", "catalyst", [EV[1]])
participant(r, "p_lu_product", "m_lit_lu_product", "product", [EV[1]], product_role="reported prior-work product")
r["notes"].append(F("The prose describes high yield and excellent enantioselectivity but Scheme 1a prints no numerical values.",
                    "qualitative prior-literature summary only", [EV[1]]))
unresolved_classification(r, [EV[1]], "Prior-literature generic drawing was not atom-mapped in this current-paper extraction.")

r = make_reaction("scheme1b_authors_prior_work", [EV[1]], record_type="literature_example", location="Scheme 1b",
                  name="Authors' prior phosphine-catalyzed [4+2] cycloaddition")
participant(r, "p_prev_allenoate", "m_lit_guo_allenoate", "reactant", [EV[1]])
participant(r, "p_cyclic_imine", "m_lit_cyclic_imine", "reactant", [EV[1]])
participant(r, "p_pr3", "m_pr3", "catalyst", [EV[1]])
participant(r, "p_prev_product", "m_lit_guo_product", "product", [EV[1]], product_role="reported prior-work product")
unresolved_classification(r, [EV[1]], "Prior-literature generic drawing was not atom-mapped in this current-paper extraction.")

r = make_reaction("scheme1c_generic", [EV[1], EV[3]], record_type="generic_scheme", location="Graphical abstract, Scheme 1c, and generic schemes above Tables 2-3")
participant(r, "p_gen_1", "m_gen_allenoate", "reactant", [EV[1], EV[3]])
participant(r, "p_gen_2", "m_gen_arylidene", "reactant", [EV[1], EV[3]])
participant(r, "p_gen_cat", "m_ph2pbn", "catalyst", [EV[1], EV[3]])
participant(r, "p_gen_dce", "m_dce", "solvent", [EV[1], EV[3]])
participant(r, "p_gen_3", "m_gen_product3", "product", [EV[1], EV[3]], product_role="generic cascade product family")
r["notes"].append(F("Ph2PBn (20 mol%), DCE, 0 °C; up to 95% yield and up to >19:1 dr", "generic scheme summary", [EV[1]]))
cc_classification(r, [EV[1], EV[3]], 3)

# Initial attempt before Table 1.
r = make_reaction("initial_attempt", [EV[2]], categories=["optimization"], location="Page 2 prose before Table 1")
participant(r, "p_1a", "m_1a", "reactant", [EV[2]], charges={"amount": charge("0.12 mmol", 0.12, "mmol", [EV[2]])})
participant(r, "p_2a", "m_2a", "reactant", [EV[2]], charges={"amount": charge("0.1 mmol", 0.1, "mmol", [EV[2]])})
participant(r, "p_cat", "m_ph3p", "catalyst", [EV[2]], charges={"mol_percent": charge("20 mol %", 20, "mol%", [EV[2]])}, eqref="m_2a")
set_solvent(r, "m_dcm", "CH2Cl2", [EV[2]])
participant(r, "p_3aa", "m_3aa", "product", [EV[2]], product_role="main tandem-reaction product")
participant(r, "p_4aa", "m_4aa", "product", [EV[2]], product_role="minor [4+2] cycloadduct")
r["reaction_scale"] = F("2a, 0.1 mmol", {"reference_molecule_id": "m_2a", "amount": {"number": 0.1, "unit": "mmol"}}, [EV[2]])
add_step(r, [EV[2]], temp=0, operation="model reaction investigated")
outcome(r, "o_yield_3aa", "yield", ["m_3aa"], "25%", {"number": 25, "unit": "%"}, [EV[2]], basis="isolated", yield_type="isolated")
outcome(r, "o_trace_4aa", "trace", ["m_4aa"], "trace amount", "trace", [EV[2]], scope="minor [4+2] cycloadduct")
r["procedure"]["original_text"] = F("allenoate 1a (0.12 mmol), alkene 2a (0.1 mmol), CH2Cl2, 0 °C, 20 mol % Ph3P", "initial-attempt conditions", [EV[2]])
r["review"]["unresolved"].append("Reaction time, solvent volume, dr, trace-product detection method, workup, and purification are not reported in the supplied main text.")
cc_classification(r, [EV[2]], 3)

table1 = [
    (1, "ph3p", "dcm", 62, 22, 1, "8:1", None),
    (2, "ph2pme", "dcm", 72, "trace", 24, "7:1", None),
    (3, "phpme2", "dcm", "NR", None, 24, None, None),
    (4, "nbu3p", "dcm", "NR", None, 24, None, None),
    (5, "ptol3p", "dcm", 60, 22, 1, "7:1", None),
    (6, "pf3p", "dcm", 63, 24, 1, "6:1", None),
    (7, "ph2pbn", "dcm", 84, "trace", 1, "13:1", None),
    (8, "dmap", "dcm", "NR", None, 24, None, None),
    (9, "dabco", "dcm", "NR", None, 24, None, None),
    (10, "dbu", "dcm", "NR", None, 24, None, None),
    (11, "ph2pbn", "dce", 88, "trace", 1, "15:1", None),
    (12, "ph2pbn", "toluene", 82, "trace", 1, "11:1", None),
    (13, "ph2pbn", "mecn", 45, "trace", 72, "10:1", None),
    (14, "ph2pbn", "dce", 82, "trace", 0.5, "10:1", {"temp": 25}),
    (15, "ph2pbn", "dce", 88, "trace", 24, "14:1", {"temp": -10}),
    (16, "ph2pbn", "dce", 85, "trace", 24, "13:1", {"loading": 10}),
    (17, "ph2pbn", "dce", 35, "trace", 1, "13:1", {"amount_2a": 0.1}),
]
for entry, cat, solvent, y3, y4, time_h, dr, exception in table1:
    ev = [EV[2]]
    r = make_reaction(f"table1_entry_{entry:02d}", ev, categories=["optimization"], location=f"Table 1 entry {entry}")
    participant(r, "p_1a", "m_1a", "reactant", ev, charges={"amount": charge("0.10 mmol", 0.10, "mmol", ev, True)})
    amount_2a = (exception or {}).get("amount_2a", 0.2)
    participant(r, "p_2a", "m_2a", "reactant", ev, charges={"amount": charge(f"{amount_2a:g} mmol", amount_2a, "mmol", ev, True)})
    loading = (exception or {}).get("loading", 20)
    cat_charges = {"mol_percent": charge(f"{loading:g} mol%", loading, "mol%", ev, True)}
    if entry != 16:
        cat_charges["amount"] = charge("0.02 mmol", 0.02, "mmol", ev, True)
    participant(r, "p_cat", f"m_{cat}", "catalyst", ev, charges=cat_charges, eqref="m_1a")
    solvent_label = simple[solvent][0]
    set_solvent(r, f"m_{solvent}", solvent_label, ev, volume_ml=1, inherited=True)
    r["reaction_scale"] = F("1a, 0.10 mmol", {"reference_molecule_id": "m_1a", "amount": {"number": 0.10, "unit": "mmol"}}, ev, method="inherited_general_condition")
    temp = (exception or {}).get("temp", 0)
    if y3 != "NR":
        participant(r, "p_3aa", "m_3aa", "product", ev, product_role="main tandem-reaction product")
        participant(r, "p_4aa", "m_4aa", "product", ev, product_role="minor [4+2] cycloadduct")
    add_step(r, ev, temp=temp, time=time_h, inherited=True)
    if y3 == "NR":
        outcome(r, "o_no_reaction", "no_reaction", [], "NR", True, ev)
        r["review"]["unresolved"].append("NR detection method, conversion, recovery, and any byproducts are not reported.")
        no_transform_classification(r, ev)
    else:
        outcome(r, "o_yield_3aa", "yield", ["m_3aa"], f"{y3}%", {"number": y3, "unit": "%"}, ev,
                basis="isolated", yield_type="isolated")
        if isinstance(y4, int):
            outcome(r, "o_yield_4aa", "yield", ["m_4aa"], f"{y4}%", {"number": y4, "unit": "%"}, ev,
                    basis="isolated", yield_type="isolated")
        elif y4 == "trace":
            outcome(r, "o_trace_4aa", "trace", ["m_4aa"], "trace", "trace", ev, scope="minor [4+2] cycloadduct")
        add_dr(r, "m_3aa", dr, ev, basis="1H NMR analysis")
        r["review"]["unresolved"].append("The dr stereoisomer identities/order beyond major:minor are not assigned; trace-product detection method is not reported.")
        cc_classification(r, ev, 3)
    r["procedure"]["general_procedure_reference"] = F("Table 1 footnote a with entry-specific footnotes d-g", "Table 1 general conditions plus row exception", ev)

# Table 2 scope entries.
for entry, substrate, product, time_h, yld, dr in table2:
    ev = [EV[2], EV[3]]
    r = make_reaction(f"table2_entry_{entry:02d}", ev, categories=["substrate_scope"], location=f"Table 2 entry {entry}")
    participant(r, "p_1a", "m_1a", "reactant", ev, charges={"amount": charge("0.1 mmol", 0.1, "mmol", [EV[3]], True)})
    participant(r, f"p_{substrate}", f"m_{substrate}", "reactant", ev, charges={"amount": charge("0.2 mmol", 0.2, "mmol", [EV[3]], True)})
    participant(r, "p_cat", "m_ph2pbn", "catalyst", ev, charges={"mol_percent": charge("20 mol %", 20, "mol%", [EV[3]], True)}, eqref="m_1a")
    set_solvent(r, "m_dce", "DCE", [EV[3]], volume_ml=1, inherited=True)
    participant(r, f"p_{product}", f"m_{product}", "product", ev, product_role="desired cascade product")
    byp = f"m_4_from_{substrate}"
    participant(r, "p_product4", byp, "product", [EV[3]], product_role="trace corresponding [4+2] cycloadduct")
    r["reaction_scale"] = F("1a, 0.1 mmol", {"reference_molecule_id": "m_1a", "amount": {"number": 0.1, "unit": "mmol"}}, [EV[3]], method="inherited_general_condition")
    add_step(r, [EV[3]], temp=0, time=time_h, inherited=True)
    if yld == "trace":
        outcome(r, f"o_trace_{product}", "trace", [f"m_{product}"], "trace", "trace", [EV[3]], scope="desired cascade product")
    else:
        outcome(r, f"o_yield_{product}", "yield", [f"m_{product}"], f"{yld}%", {"number": yld, "unit": "%"}, [EV[3]], basis="isolated", yield_type="isolated")
    outcome(r, "o_trace_product4", "trace", [byp], "Only a trace of product 4 was detected", "trace", [EV[3]], scope="corresponding [4+2] cycloadduct")
    if dr:
        add_dr(r, f"m_{product}", dr, [EV[3]], basis="1H NMR analysis")
    r["procedure"]["general_procedure_reference"] = F("Table 2 footnote a", "1a (0.1 mmol), 2 (0.2 mmol), Ph2PBn (20 mol%), 1 mL DCE, 0 °C", [EV[3]])
    r["review"]["unresolved"].append("Exact identity/stereochemistry and detection basis of trace product 4 are not printed; workup and purification are unavailable without SI.")
    cc_classification(r, [EV[3]], 3)

# Table 3 allenoate scope.
table3 = [(1, "1b", "3ba", 84, "12:1"), (2, "1c", "3ca", 59, "6:1"), (3, "1d", "3da", 89, "10:1")]
for entry, substrate, product, yld, dr in table3:
    ev = [EV[3]]
    r = make_reaction(f"table3_entry_{entry:02d}", ev, categories=["substrate_scope"], location=f"Table 3 entry {entry}")
    participant(r, f"p_{substrate}", f"m_{substrate}", "reactant", ev, charges={"amount": charge("0.1 mmol", 0.1, "mmol", ev, True)})
    participant(r, "p_2a", "m_2a", "reactant", ev, charges={"amount": charge("0.2 mmol", 0.2, "mmol", ev, True)})
    participant(r, "p_cat", "m_ph2pbn", "catalyst", ev, charges={"mol_percent": charge("20 mol %", 20, "mol%", ev, True)}, eqref=f"m_{substrate}")
    set_solvent(r, "m_dce", "DCE", ev, volume_ml=1, inherited=True)
    participant(r, f"p_{product}", f"m_{product}", "product", ev, product_role="desired cascade product")
    byp = f"m_4_from_{substrate}"
    participant(r, "p_product4", byp, "product", ev, product_role="trace corresponding [4+2] cycloadduct")
    r["reaction_scale"] = F(f"{substrate}, 0.1 mmol", {"reference_molecule_id": f"m_{substrate}", "amount": {"number": 0.1, "unit": "mmol"}}, ev, method="inherited_general_condition")
    add_step(r, ev, temp=0, inherited=True)
    outcome(r, f"o_yield_{product}", "yield", [f"m_{product}"], f"{yld}%", {"number": yld, "unit": "%"}, ev, basis="isolated", yield_type="isolated")
    outcome(r, "o_trace_product4", "trace", [byp], "A trace of product 4 was detected", "trace", ev, scope="corresponding [4+2] cycloadduct")
    add_dr(r, f"m_{product}", dr, ev, basis="1H NMR analysis")
    r["procedure"]["general_procedure_reference"] = F("Table 3 footnote a", "1 (0.1 mmol), 2a (0.2 mmol), Ph2PBn (20 mol%), 1 mL DCE, 0 °C", ev)
    r["review"]["unresolved"].append("Reaction time is not printed in Table 3; trace-product detection basis and exact identity are not reported.")
    cc_classification(r, ev, 3)

# Scheme 2 asymmetric attempts.
for catalyst, result in [("p1", "full"), ("p2", "trace"), ("p3", "trace")]:
    ev = [EV[3]]
    r = make_reaction(f"scheme2_{catalyst}_asymmetric_attempt", ev, categories=["control"], location=f"Scheme 2 / prose, catalyst {catalyst.upper()}")
    participant(r, "p_1a", "m_1a", "reactant", ev)
    participant(r, "p_2a", "m_2a", "reactant", ev)
    participant(r, "p_cat", f"m_{catalyst}", "catalyst", ev, eqref="m_1a")
    participant(r, "p_3aa", "m_3aa", "product", ev, product_role="cascade product")
    if result == "full":
        set_solvent(r, "m_dce", "DCE", ev)
        add_step(r, ev, temp=0, time=1.5)
        outcome(r, "o_yield_3aa", "yield", ["m_3aa"], "85% yield", {"number": 85, "unit": "%"}, ev)
        add_dr(r, "m_3aa", ">20:1", ev)
        outcome(r, "o_ee_3aa", "selectivity", ["m_3aa"], "3% ee", {"number": 3, "unit": "%"}, ev,
                scope={"metric": "ee", "major enantiomer": "unassigned"})
        r["review"]["unresolved"].append("Substrate/catalyst amounts, catalyst loading, yield basis, dr method, ee method, and absolute configuration are not printed for Scheme 2.")
    else:
        add_step(r, ev, operation="asymmetric-variant attempt; conditions not individually printed")
        outcome(r, "o_trace_3aa", "trace", ["m_3aa"], "trace amount of the product", "trace", ev)
        r["review"]["unresolved"].append("P2/P3-specific solvent, temperature, time, loading, amounts, and trace detection method are not reported; Scheme 2 prints conditions only over the P1 arrow.")
    cc_classification(r, ev, 3)

# Scheme 3a scale-up.
r = make_reaction("scheme3a_scale_up", [EV[3]], categories=["scale_up"], location="Scheme 3a")
participant(r, "p_1a", "m_1a", "reactant", [EV[3]], charges={"amount": charge("1 mmol", 1, "mmol", [EV[3]]), "mass": charge("0.19 g", 0.19, "g", [EV[3]])})
participant(r, "p_2a", "m_2a", "reactant", [EV[3]], charges={"amount": charge("2 mmol", 2, "mmol", [EV[3]]), "mass": charge("0.47 g", 0.47, "g", [EV[3]])})
participant(r, "p_cat", "m_ph2pbn", "catalyst", [EV[3]], charges={"mol_percent": charge("20 mol %", 20, "mol%", [EV[3]])}, eqref="m_1a")
set_solvent(r, "m_dce", "DCE", [EV[3]])
participant(r, "p_3aa", "m_3aa", "product", [EV[3]], product_role="scale-up cascade product")
r["reaction_scale"] = F("1a, 1 mmol (0.19 g)", {"reference_molecule_id": "m_1a", "amount": {"number": 1, "unit": "mmol"}}, [EV[3]])
add_step(r, [EV[3]], temp=0, time=1)
outcome(r, "o_amount_3aa", "product_amount", ["m_3aa"], "0.55 g", {"number": 0.55, "unit": "g"}, [EV[3]])
outcome(r, "o_yield_3aa", "yield", ["m_3aa"], "83% yield", {"number": 83, "unit": "%"}, [EV[3]])
add_dr(r, "m_3aa", "15:1", [EV[3]])
r["procedure"]["general_procedure_reference"] = F("under the optimal condition", "Ph2PBn (20 mol%), DCE, 0 °C, 1 h", [EV[3]])
r["review"]["unresolved"].append("DCE volume, catalyst absolute amount, yield basis, dr method, workup, and purification are not printed for the scale-up.")
cc_classification(r, [EV[3]], 3)

# Scheme 3b derivatization.
r = make_reaction("scheme3b_allylation", [EV[3]], categories=["derivatization"], location="Scheme 3b",
                  name="Allylation of product 3aa")
participant(r, "p_3aa", "m_3aa", "reactant", [EV[3]])
participant(r, "p_allyl", "m_allyl_reagent", "reactant", [EV[3]])
participant(r, "p_cat", "m_ph2pbn", "catalyst", [EV[3]], charges={"mol_percent": charge("20 mol %", 20, "mol%", [EV[3]])}, eqref="m_3aa")
set_solvent(r, "m_dce", "DCE", [EV[3]])
participant(r, "p_5a", "m_5a", "product", [EV[3]], product_role="allylation adduct")
add_step(r, [EV[3]], temp="rt", time=24)
outcome(r, "o_yield_5a", "yield", ["m_5a"], "66% yield", {"number": 66, "unit": "%"}, [EV[3]])
add_dr(r, "m_5a", "15:1", [EV[3]])
r["review"]["unresolved"].append("Reactant amounts, reaction scale, DCE volume, yield basis, dr method, addition order, workup, and purification are not printed.")
cc_classification(r, [EV[3]], 1, "Allylation visibly forms a new C-C connection; leaving-group departure, atom mapping, and any proton transfers are unresolved.")

# Scheme 4 proposed mechanism, explicitly not an experiment.
r = make_reaction("scheme4_proposed_mechanism", [EV[3], EV[4]], record_type="proposed_mechanism", location="Scheme 4",
                  name="Plausible reaction mechanism")
participant(r, "p_1a", "m_1a", "reactant", [EV[4]])
participant(r, "p_2a", "m_2a", "reactant", [EV[4]])
participant(r, "p_pr3", "m_pr3", "catalyst", [EV[4]])
participant(r, "p_3aa", "m_3aa", "product", [EV[4]], product_role="proposed cascade-path product")
participant(r, "p_4aa", "m_4aa", "product", [EV[4]], product_role="proposed alternate-path product")
for letter in "ABCDEFGHIJ":
    participant(r, f"p_{letter}", f"m_mech_{letter}", "unknown", [EV[4]])
r["notes"].append(F(
    "Phosphine attack gives A; A attacks 2a to B; proton migration gives C; C attacks a second 2a to D; intramolecular cycloaddition gives E; phosphine elimination gives F; F isomerizes to 3aa. The alternate A→G→H→I pathway gives J and then 4aa.",
    "author-proposed sequence A-J with two product branches", [EV[3], EV[4]]))
r["review"]["unresolved"] = ["Author-proposed arrow-pushing mechanism; intermediates A-J are not treated as experimentally verified."]
cc_classification(r, [EV[4]], 3, "Overall product drawings support C-C formation, but individual elementary-step mapping and proton transfers are proposed rather than experimentally established.")

# Source context: no standalone mass-spectral tables occur on the six pages.
bundle["source_context"] = [
    {"context_id": f"{DOC}:article_metadata", "document_id": DOC, "kind": "bibliographic_metadata",
     "title": TITLE, "doi": "10.1021/acs.orglett.0c00637", "evidence_ids": [TEXT[1], EV[1]]},
    {"context_id": f"{DOC}:si_notice", "document_id": DOC, "kind": "supporting_information_notice",
     "value": "Experimental procedures, characterization data, HPLC analysis data, and NMR spectra are stated to be available as SI, but no SI file was supplied in this run.",
     "evidence_ids": [TEXT[4], EV[4]]},
    {"context_id": f"{DOC}:crystallography", "document_id": DOC, "kind": "crystallographic_assignment",
     "value": {"3aa": "CCDC 1898871", "4aa": "CCDC 1921150", "claim": "relative configurations confirmed by X-ray crystallography"},
     "evidence_ids": [TEXT[2], EV[2], TEXT[6], EV[6]]},
    {"context_id": f"{DOC}:pages5_6_references", "document_id": DOC, "kind": "non_experimental_content",
     "value": "Pages 5-6 contain bibliography and the crystallographic deposition footnote; no new current-paper synthesis experiment is shown.",
     "evidence_ids": [EV[5], EV[6]]},
]

# Page-level coverage contract.
page1_ids = [f"{DOC}:scheme1a_lu_prior_work", f"{DOC}:scheme1b_authors_prior_work", f"{DOC}:scheme1c_generic"]
page2_ids = [f"{DOC}:initial_attempt"] + [f"{DOC}:table1_entry_{i:02d}" for i in range(1, 18)]
page3_ids = [f"{DOC}:scheme1c_generic"] + [f"{DOC}:table2_entry_{i:02d}" for i in range(1, 19)] + \
            [f"{DOC}:table3_entry_{i:02d}" for i in range(1, 4)] + \
            [f"{DOC}:scheme2_{p}_asymmetric_attempt" for p in ("p1", "p2", "p3")] + \
            [f"{DOC}:scheme3a_scale_up", f"{DOC}:scheme3b_allylation", f"{DOC}:scheme4_proposed_mechanism"]
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "graphical abstract and Scheme 1c", "finding": "Current-paper generic cascade scheme; repeated generic depiction merged into one generic_scheme record.", "evidence_ids": [EV[1]]},
            {"region": "Scheme 1a", "finding": "Lu group's cited prior work; classified as literature_example, not a current-paper experiment.", "evidence_ids": [EV[1]]},
            {"region": "Scheme 1b", "finding": "Authors' cited prior work; classified as literature_example, not a current-paper experiment.", "evidence_ids": [EV[1]]},
        ], "record_ids": page1_ids, "unresolved": ["Compact generic/prior-work structures were not forced into unsafe concrete SMILES."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "prose before Table 1", "finding": "Separate initial 1.2:1 attempt with 25% 3aa and trace 4aa.", "evidence_ids": [EV[2], TEXT[2]]},
            {"region": "Table 1 entries 1-17 and footnotes a-g", "finding": "All 17 optimization rows represented separately; NR rows have no fabricated products.", "evidence_ids": [EV[2], TEXT[2]]},
            {"region": "right-column Table 2 discussion", "finding": "Cross-region mentions for Table 2 scope; linked as supporting prose to the page-3 row records.", "evidence_ids": [EV[2], TEXT[2]]},
        ], "record_ids": page2_ids + [f"{DOC}:table2_entry_{i:02d}" for i in range(1, 19)],
         "unresolved": ["No workup/purification details without SI; trace/NR detection methods are not printed."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Table 2 entries 1-18 and footnotes", "finding": "All 18 arylidene scope rows extracted, including desired-product trace in entry 18 and trace product 4 in every row.", "evidence_ids": [EV[3], TEXT[3]]},
            {"region": "Table 3 entries 1-3 and footnotes", "finding": "All three allenoate scope rows extracted; reaction times are not printed.", "evidence_ids": [EV[3], TEXT[3]]},
            {"region": "Scheme 2", "finding": "P1, P2, and P3 attempts represented separately; only P1 arrow prints DCE/0 °C/1.5 h.", "evidence_ids": [EV[3], TEXT[3]]},
            {"region": "Scheme 3a-b", "finding": "Scale-up and allylation are distinct current-paper experiments.", "evidence_ids": [EV[3], TEXT[3]]},
            {"region": "mechanism discussion", "finding": "Author-proposed sequence begins here and is linked to the Scheme 4 record, not treated as experimental proof.", "evidence_ids": [EV[3], TEXT[3]]},
        ], "record_ids": page3_ids, "unresolved": ["Dense product/catalyst stereostructures lack safe isomeric SMILES and absolute descriptors."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "Scheme 4", "finding": "Plausible mechanism A-J with branches to 3aa and 4aa; classified proposed_mechanism.", "evidence_ids": [EV[4], TEXT[4]]},
            {"region": "Supporting Information notice and accession codes", "finding": "SI is mentioned but not supplied; CCDC identifiers retained as source_context.", "evidence_ids": [EV[4], TEXT[4]]},
            {"region": "conclusion and references", "finding": "No additional current-paper experiment.", "evidence_ids": [EV[4], TEXT[4]]},
        ], "record_ids": [f"{DOC}:scheme4_proposed_mechanism"], "unresolved": ["Mechanistic intermediates are proposed, not experimentally verified."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "bibliography, references 4-9", "finding": "Reference list only; no current-paper synthesis experiment or standalone characterization table.", "evidence_ids": [EV[5], TEXT[5]]},
        ], "record_ids": [], "unresolved": []},
        {"page": 6, "state": "attempted", "regions": [
            {"region": "references 10-11", "finding": "Bibliography and CCDC footnote only; no new synthesis experiment. Footnote supports crystallographic context for 3aa/4aa.", "evidence_ids": [EV[6], TEXT[6]]},
        ], "record_ids": [], "unresolved": ["The supplied main text gives no absolute R/S descriptors for 3aa or 4aa."]},
    ],
}

# Raw page decision file and literal source mentions.
raw_decisions = {
    "document_id": DOC,
    "title": TITLE,
    "pages_reviewed": [1, 2, 3, 4, 5, 6],
    "image_evidence_ids": [EV[i] for i in range(1, 7)],
    "record_ids": [r["reaction_id"] for r in bundle["reactions"]],
    "decisions": [
        "Optimization rows are separate experiments; initial 1.2:1 attempt is not merged with Table 1.",
        "Table 2 and Table 3 trace product 4 statements are retained with row-specific unresolved product identities.",
        "NR rows omit products and receive not_applicable bond classification.",
        "P2/P3 conditions are not inherited from the P1 arrow.",
        "Scheme 4 is a proposed mechanism, not mechanistic proof.",
        "Pages 5-6 are references/context, not experiments.",
        "No SI was associated from the article's availability notice.",
    ],
    "chemical_accuracy_verified": False,
}

mentions = st.source_mentions(dataset, ["1a", "2a", "3aa", "4aa", "5a", "trace", "NR"], context=180)
tasks_path = RUN / "tasks.json"
tasks = json.loads(tasks_path.read_text())
for task in tasks:
    task["state"] = "attempted"
    task["reason"] = "Page image and extracted text reviewed for current-paper experiments, tables, schemes, footnotes, references, and missing details; see coverage/pages_0001_0006.json."
dump(RUN / "submitted.json", bundle)
dump(RUN / "coverage/pages_0001_0006.json", coverage)
dump(RUN / "raw_model/pages_0001_0006_extraction.json", raw_decisions)
dump(RUN / "raw_model/stereo_checks.json", stereo_diagnostics)
dump(RUN / "raw_model/source_mentions.json", mentions)
dump(tasks_path, tasks)

print(json.dumps({
    "molecules": len(bundle["molecules"]),
    "reactions": len(bundle["reactions"]),
    "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in bundle["reactions"]),
    "records_by_type": {kind: sum(r["record_type"] == kind for r in bundle["reactions"])
                        for kind in sorted({r["record_type"] for r in bundle["reactions"]})},
}, ensure_ascii=False, indent=2))
