#!/usr/bin/env python3
"""Extend doc_1ba437add780cf5e with source-derived PDF pages 7-8."""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st


DOC = "doc_1ba437add780cf5e"
RUN = ROOT / "runs/current"
CURRENT = RUN / "standard_records.json"
INIT = RUN / "raw_model/standard_init_pages_0007_0008.json"
OUT = RUN / "submitted.json"


def T(page: int) -> str:
    return f"{DOC}:page{page}:text"


def I(page: int) -> str:
    return f"{DOC}:page{page}:image"


def ev(*pages: int) -> list[str]:
    return [eid for page in pages for eid in (I(page), T(page))]


def F(raw=None, normalized=None, evidence=(), method="original_extraction", status=None,
      derivation=None, alternatives=()):
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   derivation=derivation, alternatives=alternatives)


def q(number, unit, qualifier=None, raw_text=None):
    value = {"number": number, "unit": unit}
    if qualifier:
        value["qualifier"] = qualifier
    if raw_text:
        value["raw_text"] = raw_text
    return value


base = json.loads(CURRENT.read_text())
init = json.loads(INIT.read_text())
assert len(init["documents"]) == 1 and init["documents"][0]["document_id"] == DOC
for key, value in base["documents"][0].items():
    if key != "evidence":
        assert value == init["documents"][0][key]
assert len(init["documents"][0]["evidence"]) == 16
bundle = st.empty(init["documents"])
bundle["molecules"] = copy.deepcopy(base["molecules"])
bundle["reactions"] = copy.deepcopy(base["reactions"])
bundle["source_context"] = copy.deepcopy(base["source_context"])

molecules = {m["molecule_id"]: m for m in bundle["molecules"]}
reactions = {r["reaction_id"]: r for r in bundle["reactions"]}


def add_occurrence(mid: str, label: str, pages=(7,)):
    evidence = ev(*pages)
    m = molecules[mid]
    if not any(o["evidence_ids"] == evidence and st.value(o["compound_label"]) == label for o in m["occurrences"]):
        m["occurrences"].append({
            "document_id": DOC,
            "compound_label": F(label, label, evidence),
            "evidence_ids": evidence,
        })


def set_stereo(m, evidence, chirality, geometry, representation, assignment=None,
               atom_centers=None, double_bonds=None, axial=None, mixture=None):
    derivation = "Status assigned from the explicit Scheme 6 drawing/label; absent descriptors are not interpreted as racemic or achiral."
    m["stereochemistry"]["chirality_status"] = F("Scheme 6 drawing/label", chirality, evidence, method="inference", derivation=derivation)
    m["stereochemistry"]["geometry_status"] = F("Scheme 6 drawing/label", geometry, evidence, method="inference", derivation=derivation)
    m["stereochemistry"]["representation_status"] = F("unencoded Scheme 6 drawing", representation, evidence, method="inference", derivation=derivation)
    if assignment is not None:
        m["stereochemistry"]["assignment"] = F(assignment, assignment, evidence)
    if atom_centers is not None:
        m["stereochemistry"]["atom_centers"] = F(atom_centers, atom_centers, evidence)
    if double_bonds is not None:
        m["stereochemistry"]["double_bonds"] = F(double_bonds, double_bonds, evidence)
    if axial is not None:
        m["stereochemistry"]["axial"] = F(axial, axial, evidence)
    if mixture is not None:
        m["stereochemistry"]["mixture"] = F(mixture, mixture, evidence)


def new_molecule(mid: str, label: str, name: str, *, smiles=None, chirality="unresolved",
                 geometry="unresolved", representation=None, assignment=None, note=None):
    evidence = ev(7)
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, evidence)
    m["names"]["original"] = F(name, name, evidence)
    m["names"]["english"] = F(name, name, evidence)
    m["cas_numbers"] = F()
    if smiles:
        diag = st.stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = F(
            f"identity named in Scheme 6 caption: {name}",
            smiles,
            evidence,
            method="structure_derivation",
            derivation="Conservative candidate graph derived from the explicit source chemical name; parsing does not establish source identity.",
        )
        representation = representation or "encoded"
        m["stereochemistry"]["isomeric_smiles"] = F(
            "candidate graph canonicalization",
            diag["canonical_isomeric_smiles"],
            evidence,
            method="structure_derivation",
            derivation="Canonical isomeric SMILES computed from the candidate graph.",
        )
        m["stereochemistry"]["identity_key"] = F(
            "candidate graph stereo key",
            diag["stereo_graph_key"],
            evidence,
            method="structure_derivation",
            derivation="Stereo graph key computed from the candidate graph; it is not a sample identity.",
        )
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = F()
        representation = representation or "not_encoded"
    set_stereo(m, evidence, chirality, geometry, representation, assignment=assignment)
    if note:
        m["identity_note"] = F(note, note, evidence)
    bundle["molecules"].append(m)
    molecules[mid] = m
    return mid


# Scheme 6 reagents/solvents not present in the page-1-6 registry.
for args in [
    ("m_DDQ", "DDQ", "2,3-dicyano-5,6-dichlorobenzoquinone", "N#CC1=C(Cl)C(=O)C(C#N)=C(Cl)C1=O"),
    ("m_H2O", "H2O", "water", "O"),
    ("m_benzene", "benzene", "benzene", "c1ccccc1"),
    ("m_EtOAc", "EtOAc", "ethyl acetate", "CCOC(C)=O"),
    ("m_EtOH", "EtOH", "ethanol", "CCO"),
    ("m_NaOH", "NaOH", "sodium hydroxide", "[Na+].[OH-]"),
]:
    if args[0] not in molecules:
        new_molecule(args[0], args[1], args[2], smiles=args[3], chirality="achiral", geometry="not_applicable", representation="encoded",
                     assignment="No stereogenic element in the named candidate graph.")

if "m_HCl" not in molecules:
    new_molecule("m_HCl", "HCl", "hydrogen chloride", chirality="achiral",
                 geometry="not_applicable", representation="not_encoded",
                 assignment="Diatomic hydrogen chloride is not stereogenic.",
                 note="The source specifies concentrated HCl in EtOH (6 M) and a later HCl (1 M) workup; volumes are not reported.")

# Add page-7 occurrences and improve only the Scheme-6-specific stereo facts.
for mid, label in [
    ("m_2h", "2h"),
    ("m_2S_3S_3ha", "(2S,3S)-3ha"),
    ("m_6a", "6a"),
    ("m_7", "7"),
    ("m_8a", "8a"),
    ("m_mycF", "mycestericin F·HCl"),
    ("m_3hb", "(2S,3R)-4hb"),
    ("m_mycG", "mycestericin G·HCl"),
]:
    add_occurrence(mid, label)

set_stereo(
    molecules["m_2h"], ev(7), "unspecified", "not_applicable", "not_encoded",
    assignment="Scheme 6 draws a wavy bond at the oxazoline ester-bearing center; no single configuration is assigned.",
    atom_centers=[{"center": "oxazoline ester-bearing carbon", "assignment": "unspecified (wavy bond)"}],
)
molecules["m_2h"]["identity_note"] = F(
    "Oxazoline 2h is explicitly drawn in Scheme 6 with a wavy bond; no complete candidate SMILES or unique configuration is asserted.",
    "Source-identified oxazoline 2h with an unspecified center; complete structure remains unencoded.",
    ev(7),
)
set_stereo(
    molecules["m_6a"], ev(7), "specified", "not_applicable", "not_encoded",
    assignment="Source-drawn stereochemistry inherited through PMB deprotection from labelled (2S,3S)-3ha; exact center indices/R,S were not independently encoded.",
    atom_centers=[{"assignment": "two source-drawn stereocenters; wedge/dash shown"}],
)
molecules["m_6a"]["identity_note"] = F(
    "Scheme 6 supplies a complete drawing for alcohol 6a, but a source-verified isomeric SMILES was not asserted.",
    "Complete Scheme 6 drawing retained; structure encoding remains unresolved.",
    ev(7),
)
set_stereo(
    molecules["m_7"], ev(7), "achiral", "not_applicable", "not_encoded",
    assignment="The source drawing is an acyclic ketone with a terminal alkene and no stereogenic center; no E/Z element applies to the terminal alkene.",
)
molecules["m_7"]["identity_note"] = F(
    "Scheme 6 draws ketoalkene 7 as C6H13-C(O)-(CH2)8-terminal alkene; exact candidate-SMILES transcription was not asserted.",
    "Source-identified ketoalkene 7; complete structure remains unencoded.",
    ev(7),
)
set_stereo(
    molecules["m_8a"], ev(7), "specified", "unspecified", "not_encoded",
    assignment="Two source-drawn stereocenters are retained from 6a; the newly formed internal alkene is drawn but not assigned E/Z in Scheme 6.",
    atom_centers=[{"assignment": "two source-drawn stereocenters; wedge/dash shown"}],
    double_bonds=[{"bond": "new internal alkene from metathesis", "assignment": "E/Z not explicitly named"}],
)
molecules["m_8a"]["identity_note"] = F(
    "Scheme 6 supplies the 8a drawing and stereocenters, but the internal-alkene E/Z descriptor and a complete isomeric SMILES are not reported/encoded.",
    "Source-identified 8a with specified centers and unresolved alkene geometry.",
    ev(7),
)
set_stereo(
    molecules["m_mycF"], ev(7), "specified", "not_applicable", "not_encoded",
    assignment="Scheme 6 draws the stereochemistry of mycestericin F·HCl; a full isomeric SMILES was not asserted.",
    atom_centers=[{"assignment": "source-drawn stereocenters; exact atom indexing unresolved"}],
)
molecules["m_mycF"]["names"]["original"] = F(
    ["mycestericin F (page-6 prose)", "mycestericin F·HCl (Scheme 6)"],
    "mycestericin F hydrochloride (Scheme 6 product; prose short name mycestericin F)",
    ev(6, 7),
)
molecules["m_mycF"]["names"]["english"] = copy.deepcopy(molecules["m_mycF"]["names"]["original"])
molecules["m_mycF"]["identity_note"] = F(
    "Scheme 6 explicitly labels the isolated endpoint as mycestericin F·HCl; page-6 prose uses the unsalted short name mycestericin F.",
    "Scheme endpoint is the hydrochloride salt; prose uses the short natural-product name.",
    ev(6, 7),
)

set_stereo(
    molecules["m_3hb"], ev(6, 7), "specified", "not_applicable", "not_encoded",
    assignment="(2S,3R) in both page-6 prose and Scheme 6; source compound number conflicts (3hb versus 4hb).",
    atom_centers=[{"assignment": "(2S,3R)", "source_label_conflict": ["3hb", "4hb"]}],
)
molecules["m_3hb"]["names"]["original"] = F(
    ["(2S,3R)-3hb (page-6 prose)", "(2S,3R)-4hb (Scheme 6)"],
    "source-labelled syn-coupling product (page-6 3hb / Scheme-6 4hb)",
    ev(6, 7),
    method="inference",
    derivation="Neutral registry name retains the explicit source label disagreement without selecting one label.",
)
molecules["m_3hb"]["names"]["english"] = copy.deepcopy(molecules["m_3hb"]["names"]["original"])
molecules["m_3hb"]["identity_note"] = F(
    ["page-6 prose: (2S,3R)-3hb", "page-7 Scheme 6: (2S,3R)-4hb"],
    None,
    ev(6, 7),
    status="conflict",
    alternatives=[
        F("(2S,3R)-3hb", "(2S,3R)-3hb", ev(6)),
        F("(2S,3R)-4hb", "(2S,3R)-4hb", ev(7)),
    ],
)

set_stereo(
    molecules["m_mycG"], ev(7), "specified", "not_applicable", "not_encoded",
    assignment="Scheme 6 draws the stereochemistry of mycestericin G·HCl; a full isomeric SMILES was not asserted.",
    atom_centers=[{"assignment": "source-drawn stereocenters; exact atom indexing unresolved"}],
)
molecules["m_mycG"]["names"]["original"] = F(
    ["mycestericin G (page-6 prose)", "mycestericin G·HCl (Scheme 6)"],
    "mycestericin G hydrochloride (Scheme 6 product; prose short name mycestericin G)",
    ev(6, 7),
)
molecules["m_mycG"]["names"]["english"] = copy.deepcopy(molecules["m_mycG"]["names"]["original"])
molecules["m_mycG"]["identity_note"] = F(
    "Scheme 6 explicitly labels the isolated endpoint as mycestericin G·HCl; page-6 prose uses the unsalted short name mycestericin G.",
    "Scheme endpoint is the hydrochloride salt; prose uses the short natural-product name.",
    ev(6, 7),
)

# New intermediates revealed only by the page-7 scheme.
if "m_6b" not in molecules:
    new_molecule(
        "m_6b", "6b", "alcohol 6b", chirality="specified", geometry="not_applicable",
        representation="not_encoded",
        assignment="Source-drawn stereochemistry retained from labelled (2S,3R)-4hb through PMB deprotection; exact atom indices/R,S are not encoded.",
        note="The complete constitution and stereochemistry were not converted to a source-verified isomeric SMILES; the Scheme 6 drawing is retained as evidence.",
    )
    molecules["m_6b"]["stereochemistry"]["atom_centers"] = F(
        "Scheme 6 wedge/dash drawing", [{"assignment": "two source-drawn stereocenters"}], ev(7)
    )
if "m_8b" not in molecules:
    new_molecule(
        "m_8b", "8b", "metathesis product 8b", chirality="specified", geometry="unspecified",
        representation="not_encoded",
        assignment="Two source-drawn stereocenters are retained from 6b; the new internal alkene is not explicitly assigned E/Z.",
        note="The complete constitution and stereochemistry were not converted to a source-verified isomeric SMILES; the Scheme 6 drawing is retained as evidence.",
    )
    molecules["m_8b"]["stereochemistry"]["atom_centers"] = F(
        "Scheme 6 wedge/dash drawing", [{"assignment": "two source-drawn stereocenters"}], ev(7)
    )
    molecules["m_8b"]["stereochemistry"]["double_bonds"] = F(
        "Scheme 6 internal alkene drawing", [{"bond": "metathesis alkene", "assignment": "E/Z not explicitly named"}], ev(7)
    )


TAX = {"C-H": 1, "C-C": 2, "O-H": 24, "N-H": 29, "DB-C-C": 13}


def classification(r, codes, evidence, limitation, extra_changes=()):
    changes = []
    labels = []
    for code in codes:
        a, b = code.removeprefix("DB-").split("-")
        order = 2 if code.startswith("DB-") else 1
        cid = f"b{len(changes) + 1}"
        change = {
            "change_id": cid,
            "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
            "before_order": 0,
            "after_order": order,
            "change_type": "hydrogen_gain" if b == "H" else "formed",
        }
        changes.append(change)
        labels.append({"category_number": TAX[code], "code": code, "supporting_change_ids": [cid]})
    for item in extra_changes:
        change = dict(item)
        change["change_id"] = f"b{len(changes) + 1}"
        changes.append(change)
    r["classification"]["bond_changes"] = [
        F("Scheme 6 transformation drawing", change, evidence, method="structure_derivation", derivation=limitation)
        for change in changes
    ]
    r["classification"]["labels"] = F(
        "Scheme 6 transformation drawing", labels, evidence, method="structure_derivation", derivation=limitation
    )
    r["classification"]["coverage"] = F(
        "partial class-level changes", "partial", evidence, method="inference", derivation=limitation
    )
    r["classification"]["limitations"] = [
        limitation,
        "No complete atom mapping is asserted; exact atom identities, coproducts, and additional bond changes may be absent.",
    ]


def reaction(rid, name, transform, location, categories, evidence, unresolved):
    r = st.new_reaction(rid, DOC, evidence)
    r["record_type"] = "specific_experiment"
    r["experiment_categories"] = F(categories, categories, evidence)
    r["location"] = F(location, {"pdf_pages": sorted({int(x.split(":page")[1].split(":")[0]) for x in evidence}), "region": location}, evidence)
    r["reaction_name"] = F(name, name, evidence)
    r["transformation_type"] = F(transform, transform, evidence)
    r["general_scheme_ids"] = []
    r["review"] = {
        "status": "needs_review",
        "method": "assistant_visual",
        "reviewer": None,
        "checked_fields": [],
        "unresolved": unresolved,
    }
    r["completeness"] = "partial"
    return r


def addp(r, pid, mid, role, evidence, *, stage="reaction", product_role=None, molpct=None, ref=None):
    p = st.new_participant(pid, mid, role, evidence, stage)
    if product_role:
        p["product_role"] = F(product_role, product_role, evidence)
    if molpct is not None:
        p["charges"]["mol_percent"] = F(f"{molpct} mol%", q(molpct, "mol%"), evidence)
        p["equivalence_reference"] = F(ref, ref, evidence)
    r["participants"].append(p)
    return pid


def add_step(r, sid, evidence, pids, operation, *, temperature=None, time=None,
             atmosphere=None, pressure=None, stage="reaction"):
    s = st.new_step(sid, evidence, stage)
    s["participant_ids"] = pids
    s["operation"] = F(operation, operation, evidence)
    if temperature is not None:
        temperature_raw = temperature if isinstance(temperature, str) else temperature.get(
            "raw_text", f"{temperature['number']} {temperature['unit']}"
        )
        s["temperature"] = F(temperature_raw, temperature, evidence)
    if time is not None:
        s["time"] = F(time.get("raw_text", f"{time['number']} {time['unit']}") if isinstance(time, dict) else time, time, evidence)
    if atmosphere is not None:
        s["atmosphere"] = F(atmosphere, atmosphere, evidence)
    if pressure is not None:
        pressure_raw = pressure if isinstance(pressure, str) else pressure.get(
            "raw_text", f"{pressure['number']} {pressure['unit']}"
        )
        s["pressure"] = F(pressure_raw, pressure, evidence)
    r["steps"].append(s)


def add_outcome(r, oid, kind, targets, evidence, raw, normalized, *, product_role=None,
                yield_type=None, basis=None, scope=None):
    o = st.new_outcome(oid, kind, targets, evidence)
    o["value"] = F(raw, normalized, evidence)
    if product_role:
        o["product_role"] = F(product_role, product_role, evidence)
    if yield_type:
        o["yield_type"] = F(yield_type, yield_type, evidence)
    if basis:
        o["measurement_basis"] = F(basis, basis, evidence)
    if scope is not None:
        o["scope"] = F(scope, scope, evidence)
    r["outcomes"].append(o)


def put(r):
    old = reactions.get(r["reaction_id"])
    if old is None:
        bundle["reactions"].append(r)
    else:
        bundle["reactions"][bundle["reactions"].index(old)] = r
    reactions[r["reaction_id"]] = r


# Two 2 mmol-scale coupling arrows.
for rid, catalyst, ligand, metal, catalyst_system, product, yld, dr, ee_value, label_note in [
    ("doc_1ba437add780cf5e:scheme6_scaled_3ha", "m_Pd_L13", "m_L1", "m_CuSalt", "L13-Pd/L1-Cu", "m_2S_3S_3ha", 60, "6:1", 99, None),
    ("doc_1ba437add780cf5e:scheme6_syn_3hb", "m_Pd_L3", "m_L3", "m_AgOAc", "L3-Pd/L3-Ag", "m_3hb", 63, "14:1", 96,
     "Page-6 prose labels the product (2S,3R)-3hb, whereas Scheme 6 labels the same route/drawing (2S,3R)-4hb."),
]:
    evidence = ev(6, 7)
    r = reaction(
        rid,
        "scaled coupling for mycestericin synthesis",
        "diastereodivergent alkoxyallene/oxazoline coupling",
        "page-6 prose and page-7 Scheme 6 coupling arrow",
        ["scale_up"],
        evidence,
        [
            "Scheme 6 states only '2 mmol scale'; it does not explicitly identify the scale denominator, individual charges, solvent, temperature, time, workup, purification, or yield/selectivity measurement methods.",
            *([label_note] if label_note else []),
        ],
    )
    addp(r, "p_1a", "m_1a", "reactant", evidence)
    addp(r, "p_2h", "m_2h", "reactant", evidence)
    addp(r, "p_product", product, "product", evidence, product_role="major coupling product")
    addp(r, "p_pd", catalyst, "catalyst", evidence)
    addp(r, "p_ligand", ligand, "ligand", evidence)
    addp(r, "p_metal", metal, "catalyst", evidence)
    # The scheme does not explicitly identify which reactant the scale is based on;
    # retain the source quantity without inventing the required reference molecule.
    r["reaction_scale"] = F("2 mmol scale", None, ev(7))
    add_step(r, "coupling", evidence, ["p_1a", "p_2h", "p_pd", "p_ligand", "p_metal"],
             f"couple 1a and 2h using the displayed {catalyst_system} catalyst system on the stated 2 mmol scale")
    add_outcome(r, "yield", "yield", [product], ev(7), f"{yld}%", q(yld, "%"), product_role="coupling product")
    add_outcome(r, "dr", "selectivity", [product], ev(7), f"dr = {dr}", dr,
                scope={"metric": "dr", "major": st.value(molecules[product]["stereochemistry"]["assignment"]), "minor_configuration": "unreported"})
    add_outcome(r, "ee", "selectivity", [product], ev(7), f"{ee_value}% ee", q(ee_value, "%"),
                scope={"metric": "ee", "configuration": st.value(molecules[product]["stereochemistry"]["assignment"])})
    if label_note:
        r["notes"].append(F(
            ["page-6 prose: (2S,3R)-3hb", "page-7 Scheme 6: (2S,3R)-4hb"],
            None,
            evidence,
            status="conflict",
            alternatives=[F("(2S,3R)-3hb", "(2S,3R)-3hb", ev(6)), F("(2S,3R)-4hb", "(2S,3R)-4hb", ev(7))],
        ))
    classification(r, ["C-C", "N-H"], evidence,
                   "Partial net changes inferred from the coupling drawing and primary-amine-forming workup; exact oxazoline/product mapping is unresolved.")
    put(r)


# PMB deprotection arrows a.
for rid, start, product, yld in [
    ("doc_1ba437add780cf5e:scheme6_deprotection_6a", "m_2S_3S_3ha", "m_6a", 91),
    ("doc_1ba437add780cf5e:scheme6_deprotection_6b", "m_3hb", "m_6b", 90),
]:
    evidence = ev(6, 7) if product == "m_6a" else ev(7)
    r = reaction(
        rid,
        "PMB deprotection",
        "oxidative para-methoxybenzyl ether deprotection",
        f"page-7 Scheme 6 arrow a to {product.removeprefix('m_')}",
        ["derivatization"],
        evidence,
        ["Scheme 6 does not report substrate/reagent amounts, solvent ratio/volume, workup, purification, or yield measurement basis."],
    )
    addp(r, "p_start", start, "reactant", evidence)
    addp(r, "p_product", product, "product", evidence, product_role="deprotected alcohol")
    addp(r, "p_DDQ", "m_DDQ", "reagent", ev(7))
    addp(r, "p_DCM", "m_DCM", "solvent", ev(7))
    addp(r, "p_water", "m_H2O", "solvent", ev(7))
    r["solvent_conditions"]["composition"] = F("CH2Cl2/H2O", ["m_DCM", "m_H2O"], ev(7))
    r["procedure"]["original_text"] = F(
        "(a) 2,3-dicyano-5,6-dichlorobenzoquinone, CH2Cl2/H2O, rt, 0.5 h",
        "DDQ in CH2Cl2/H2O at room temperature for 0.5 h",
        ev(7),
    )
    add_step(r, "deprotection", ev(7), ["p_start", "p_DDQ", "p_DCM", "p_water"],
             "treat with DDQ in CH2Cl2/H2O", temperature="rt", time=q(0.5, "h"))
    add_outcome(r, "yield", "yield", [product], ev(7), f"{yld}%", q(yld, "%"), product_role="deprotected alcohol")
    classification(
        r,
        ["O-H"],
        ev(7),
        "Scheme 6 shows conversion of the PMB ether to an alcohol; O-H formation and C-O cleavage are retained only at class level.",
        extra_changes=[{
            "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}],
            "before_order": 1,
            "after_order": 0,
            "change_type": "broken",
        }],
    )
    put(r)


# Cross-metathesis arrows b.
for rid, start, product, yld in [
    ("doc_1ba437add780cf5e:scheme6_metathesis_8a", "m_6a", "m_8a", 65),
    ("doc_1ba437add780cf5e:scheme6_metathesis_8b", "m_6b", "m_8b", 45),
]:
    evidence = ev(6, 7) if product == "m_8a" else ev(7)
    r = reaction(
        rid,
        "olefin metathesis",
        "cross metathesis",
        f"page-7 Scheme 6 arrow b to {product.removeprefix('m_')}",
        ["derivatization"],
        evidence,
        ["Scheme 6 does not report reactant/catalyst amounts, benzene volume, workup, purification, alkene E/Z assignment, or yield measurement basis."],
    )
    addp(r, "p_start", start, "reactant", evidence)
    addp(r, "p_7", "m_7", "reactant", evidence)
    addp(r, "p_product", product, "product", evidence, product_role="cross-metathesis product")
    addp(r, "p_HGII", "m_HGII", "catalyst", evidence)
    addp(r, "p_benzene", "m_benzene", "solvent", ev(7))
    r["solvent_conditions"]["composition"] = F("benzene", ["m_benzene"], ev(7))
    r["procedure"]["original_text"] = F(
        "(b) Hoveyda-Grubbs II catalyst, benzene, 80 °C, 24 h",
        "Hoveyda-Grubbs II catalyst in benzene at 80 °C for 24 h",
        ev(7),
    )
    add_step(r, "cross_metathesis", ev(7), ["p_start", "p_7", "p_HGII", "p_benzene"],
             "cross-metathesis with ketoalkene 7 and Hoveyda-Grubbs II catalyst",
             temperature=q(80, "°C"), time=q(24, "h"))
    add_outcome(r, "yield", "yield", [product], ev(7), f"{yld}%", q(yld, "%"), product_role="cross-metathesis product")
    classification(r, ["DB-C-C"], ev(7),
                   "A new internal C=C connection is inferred from the two terminal-alkene reactants and Scheme 6 metathesis product; coproduct and exact atom mapping are unresolved.")
    put(r)


# Three-stage c-e conversions. The printed yields are overall from 8a/8b to the HCl salts.
for rid, start, product, yld in [
    ("doc_1ba437add780cf5e:scheme6_to_mycestericin_F", "m_8a", "m_mycF", 56),
    ("doc_1ba437add780cf5e:scheme6_to_mycestericin_G", "m_8b", "m_mycG", 60),
]:
    evidence = ev(6, 7)
    letter = "F" if product == "m_mycF" else "G"
    r = reaction(
        rid,
        f"completion of mycestericin {letter} hydrochloride synthesis",
        "hydrogenation, oxazoline/ester hydrolysis, saponification, and acidification",
        f"page-7 Scheme 6 arrows c-e to mycestericin {letter}·HCl",
        ["derivatization"],
        evidence,
        [
            "The 56%/60% value spans Scheme 6 steps c-e from 8a/8b to the HCl salt; it is not assigned to the last workup step alone.",
            "Scheme 6 does not report substrate/reagent volumes or amounts except Pd/C loading, solution concentrations, H2 pressure, stage-c/d times, and temperatures; stage-e time/temperature, intermediates, purification, and yield measurement basis remain unreported.",
        ],
    )
    addp(r, "p_start", start, "reactant", evidence)
    addp(r, "p_product", product, "product", evidence, product_role=f"mycestericin {letter} hydrochloride")
    addp(r, "p_PdC", "m_PdC", "catalyst", ev(7), molpct=10, ref=start)
    addp(r, "p_H2", "m_H2", "reagent", ev(7))
    addp(r, "p_EtOAc", "m_EtOAc", "solvent", ev(7))
    addp(r, "p_HCl6M", "m_HCl", "reagent", ev(7))
    addp(r, "p_EtOH", "m_EtOH", "solvent", ev(7))
    addp(r, "p_NaOH", "m_NaOH", "reagent", ev(7))
    addp(r, "p_MeOH", "m_MeOH", "solvent", ev(7))
    addp(r, "p_HCl1M", "m_HCl", "workup_agent", ev(7), stage="workup")
    r["solvent_conditions"]["composition"] = F(
        "EtOAc (c); EtOH (d); MeOH (e)",
        {"c": ["m_EtOAc"], "d": ["m_EtOH"], "e": ["m_MeOH"]},
        ev(7),
    )
    r["solvent_conditions"]["concentration"] = F(
        "HCl in EtOH (6 M); NaOH (1 M); HCl workup (1 M)",
        None,
        ev(7),
    )
    r["solvent_conditions"]["concentration_basis"] = F(
        "caption solution concentrations",
        {"d_HCl": "m_HCl in m_EtOH", "e_NaOH": "m_NaOH solution", "e_HCl_workup": "m_HCl workup solution"},
        ev(7),
    )
    r["procedure"]["original_text"] = F(
        "(c) Pd/C (10 mol%), H2 (1 atm balloon), EtOAc, rt, 12 h. (d) Concentrated HCl in EtOH (6 M), reflux, 12 h. (e) NaOH (1 M), MeOH then HCl (1 M) workup.",
        "Scheme 6 steps c-e",
        ev(7),
    )
    add_step(r, "c_hydrogenation", ev(7), ["p_start", "p_PdC", "p_H2", "p_EtOAc"],
             "hydrogenate with Pd/C in EtOAc", temperature="rt", time=q(12, "h"),
             atmosphere="H2", pressure="1 atm balloon")
    add_step(r, "d_acid_hydrolysis", ev(7), ["p_start", "p_HCl6M", "p_EtOH"],
             "heat with concentrated HCl in EtOH (6 M)", temperature="reflux", time=q(12, "h"))
    add_step(r, "e_base_hydrolysis", ev(7), ["p_start", "p_NaOH", "p_MeOH"],
             "treat with NaOH (1 M) in MeOH")
    add_step(r, "e_hcl_workup", ev(7), ["p_HCl1M"], "work up with HCl (1 M)", stage="workup")
    add_outcome(
        r,
        "overall_yield_c_to_e",
        "yield",
        [product],
        ev(7),
        f"c-e, {yld}%",
        q(yld, "%"),
        product_role=f"mycestericin {letter} hydrochloride",
        yield_type="overall yield across Scheme 6 steps c-e",
        scope={"metric": "overall_yield", "from_molecule_id": start, "to_molecule_id": product, "stages": ["c", "d", "e"]},
    )
    r["notes"].append(F(
        f"The printed c-e arrow gives {yld}% from {start.removeprefix('m_')} to mycestericin {letter}·HCl.",
        {"yield_span": {"from": start, "to": product, "stages": ["c", "d", "e"]}, "not_last_step_yield": True},
        ev(7),
        method="inference",
        derivation="The yield is printed beside a single arrow explicitly labelled c-e, so it spans all three displayed operations.",
    ))
    r["notes"].append(F(
        "HCl in EtOH (6 M); NaOH (1 M); HCl workup (1 M)",
        [
            {"stage": "d", "molecule_id": "m_HCl", "concentration": q(6, "M")},
            {"stage": "e", "molecule_id": "m_NaOH", "concentration": q(1, "M")},
            {"stage": "e_workup", "molecule_id": "m_HCl", "concentration": q(1, "M")},
        ],
        ev(7),
    ))
    classification(
        r,
        ["C-H", "O-H", "N-H"],
        ev(7),
        "Partial class-level changes inferred from hydrogenation plus oxazoline/ester hydrolysis in Scheme 6; exact atoms and full hydrolysis accounting are unresolved.",
        extra_changes=[{
            "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
            "before_order": 2,
            "after_order": 1,
            "change_type": "order_decreased",
        }],
    )
    put(r)


context_id = f"{DOC}:pages7_8_reference_and_negative_inventory"
bundle["source_context"] = [c for c in bundle["source_context"] if c.get("context_id") != context_id]
bundle["source_context"].append({
    "context_id": context_id,
    "document_id": DOC,
    "kind": "coverage_context",
    "evidence_ids": ev(7, 8),
    "content": {
        "page7": "Scheme 6 contains eight current-paper reaction arrows; conflict-of-interest, keywords, and references [1]-[14] are non-experimental context.",
        "page8": "References [14]-[30], CCDC availability notice, and manuscript dates only; no current-paper reaction arrow, table, or procedure.",
        "characterization_inventory": "No standalone mass-spectrum peak table, molecular-weight table, or fragment-relative-abundance block appears on pages 7-8.",
        "ccdc_note": "Reference [25] mentions CCDC 2036502 for (2S,3S)-3ha; the external crystallographic file was not supplied and is not treated as associated SI.",
    },
})

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 7,
            "state": "attempted",
            "regions": [
                {
                    "region": "left column, Scheme 6 upper route",
                    "finding": "Four current-paper arrows: 2 mmol L13-Pd/L1-Cu coupling to (2S,3S)-3ha; DDQ deprotection to 6a; Hoveyda-Grubbs-II cross-metathesis with 7 to 8a; three-stage c-e conversion to mycestericin F·HCl.",
                    "evidence_ids": ev(7),
                },
                {
                    "region": "left column, Scheme 6 lower route",
                    "finding": "Four current-paper arrows: 2 mmol L3-Pd/L3-Ag coupling to scheme-labelled (2S,3R)-4hb; DDQ deprotection to 6b; cross-metathesis with 7 to 8b; three-stage c-e conversion to mycestericin G·HCl. The 4hb scheme label conflicts with page-6 prose label 3hb.",
                    "evidence_ids": ev(7),
                },
                {
                    "region": "Scheme 6 caption and footnote conditions a-e",
                    "finding": "Explicit reagents, solvents, temperatures/times, Pd/C loading, H2 pressure, solution concentrations, and PMB definition were assigned to the corresponding arrows; c-e yields were retained as overall multistep yields.",
                    "evidence_ids": ev(7),
                },
                {
                    "region": "conflict-of-interest, keywords, and references [1]-[14]",
                    "finding": "Non-experimental context and cited literature; no additional current-paper experiment was created.",
                    "evidence_ids": ev(7),
                },
            ],
            "record_ids": [
                f"{DOC}:scheme6_scaled_3ha",
                f"{DOC}:scheme6_deprotection_6a",
                f"{DOC}:scheme6_metathesis_8a",
                f"{DOC}:scheme6_to_mycestericin_F",
                f"{DOC}:scheme6_syn_3hb",
                f"{DOC}:scheme6_deprotection_6b",
                f"{DOC}:scheme6_metathesis_8b",
                f"{DOC}:scheme6_to_mycestericin_G",
            ],
            "unresolved": [
                "Page-6 prose labels the syn product (2S,3R)-3hb, while page-7 Scheme 6 labels it (2S,3R)-4hb; both are retained without choosing.",
                "Scheme 6 does not state most charges/volumes, coupling temperature/time, purification, yield measurement methods, dr/ee methods, or E/Z assignments for 8a/8b.",
                "Complex Scheme 6 structures were not forced into complete isomeric SMILES; source drawings and explicit stereo assignments remain primary evidence.",
            ],
        },
        {
            "page": 8,
            "state": "attempted",
            "regions": [
                {
                    "region": "two-column references [14]-[30]",
                    "finding": "Bibliography only; cited papers are not current-paper experiments.",
                    "evidence_ids": ev(8),
                },
                {
                    "region": "references [24]-[28] mechanistic/crystallographic notes",
                    "finding": "Author discussion and external CCDC 2036502 availability notice; no supplied CIF/SI and no new executed experiment on this page.",
                    "evidence_ids": ev(8),
                },
                {
                    "region": "manuscript received/revised/accepted/version dates",
                    "finding": "Publication metadata; no reaction arrow, table row, experimental procedure, or standalone mass-spectrum/fragment table.",
                    "evidence_ids": ev(8),
                },
            ],
            "record_ids": [],
            "unresolved": [
                "External references and CCDC 2036502 were not supplied as associated documents and were not imported.",
            ],
        },
    ],
}
(RUN / "coverage/pages_0007_0008.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")

raw = {
    "document_id": DOC,
    "page_range": [7, 8],
    "visual_inspection": "Original rendered page 7 and page 8 images opened at original detail; Scheme 6 was additionally inspected in a source-lineage crop.",
    "page7_arrow_inventory": [
        "1a + 2h -> (2S,3S)-3ha (60%, dr 6:1, 99% ee; L13-Pd/L1-Cu; 2 mmol scale)",
        "(2S,3S)-3ha -> 6a (a; 91%)",
        "6a + 7 -> 8a (b; 65%)",
        "8a -> mycestericin F·HCl (c-e; 56% overall)",
        "1a + 2h -> (2S,3R)-4hb (63%, dr 14:1, 96% ee; L3-Pd/L3-Ag; 2 mmol scale)",
        "(2S,3R)-4hb -> 6b (a; 90%)",
        "6b + 7 -> 8b (b; 45%)",
        "8b -> mycestericin G·HCl (c-e; 60% overall)",
    ],
    "caption_conditions": {
        "a": "2,3-dicyano-5,6-dichlorobenzoquinone, CH2Cl2/H2O, rt, 0.5 h",
        "b": "Hoveyda-Grubbs II catalyst, benzene, 80 °C, 24 h",
        "c": "Pd/C (10 mol%), H2 (1 atm balloon), EtOAc, rt, 12 h",
        "d": "concentrated HCl in EtOH (6 M), reflux, 12 h",
        "e": "NaOH (1 M), MeOH then HCl (1 M) workup",
    },
    "deduplication": [
        "The upper coupling/deprotection/metathesis/final-route arrows enrich the same records introduced by page-6 prose.",
        "The lower page-6 prose route is split into the three distinct page-7 arrows 4hb->6b, 6b+7->8b, and 8b->mycestericin G·HCl; the broad prose route is not retained as a duplicate experiment.",
    ],
    "conflicts": [
        "Page-6 prose: (2S,3R)-3hb; page-7 Scheme 6: (2S,3R)-4hb. Both identify the same syn route context but disagree on compound number.",
    ],
    "page8": "References and publication metadata only; zero reaction records.",
}
(RUN / "raw_model/page_extraction_pages_0007_0008.json").write_text(json.dumps(raw, ensure_ascii=False, indent=2) + "\n")

stereo_checks = {}
for mid in ["m_DDQ", "m_H2O", "m_benzene", "m_EtOAc", "m_EtOH", "m_NaOH"]:
    candidate = st.value(molecules[mid]["smiles"])
    stereo_checks[mid] = st.stereo_inspect(candidate) if candidate else {"candidate_smiles": None}
(RUN / "raw_model/stereo_checks_pages_0007_0008.json").write_text(
    json.dumps(stereo_checks, ensure_ascii=False, indent=2) + "\n"
)

OUT.write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({
    "reactions": len(bundle["reactions"]),
    "molecules": len(bundle["molecules"]),
    "source_context": len(bundle["source_context"]),
    "page7_record_ids": len(coverage["pages"][0]["record_ids"]),
}, ensure_ascii=False))
