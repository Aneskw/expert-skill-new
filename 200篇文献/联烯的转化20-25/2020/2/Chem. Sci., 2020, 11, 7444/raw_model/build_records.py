#!/usr/bin/env python3
"""Build evidence-linked v3.1 records for DOI 10.1039/d0sc02878a, PDF pages 1-6.

This is the raw machine extraction recipe.  It deliberately leaves complex
source drawings unresolved rather than manufacturing exact structures.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st  # noqa: E402


RUN = ROOT / "runs/current"
DOC = "doc_48a246e692c52136"
E = {p: f"{DOC}:page{p}:image" for p in range(1, 7)}
T = {p: f"{DOC}:page{p}:text" for p in range(1, 7)}


def f(raw=None, normalized=None, evidence=(), method="original_extraction",
      status=None, derivation=None, alternatives=()):
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   derivation=derivation, alternatives=alternatives)


def missing(evidence=()):
    return f(evidence=evidence)


def na(evidence=(), reason=None):
    return f(evidence=evidence, method="not_applicable", status="not_applicable",
             derivation=reason)


bundle = json.loads((RUN / "submitted.json").read_text())
molecules = {}


def base_molecule(mid, label, evidence, name=None, chirality="unresolved",
                  geometry="unresolved", representation="unresolved",
                  assignment=None, mixture=None, identity_note=None):
    evidence = list(dict.fromkeys(evidence))
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    m["names"]["original"] = f(name or label, name or label, evidence)
    m["names"]["english"] = f(name or label, name or label, evidence)
    if name and name != label:
        m["names"]["abbreviations"] = f(label, [label], evidence)
    m["cas_numbers"] = missing(evidence)
    m["stereochemistry"]["chirality_status"] = f(
        chirality, chirality, evidence, method="inference",
        derivation="Status records what the supplied drawing/text specifies; absence of an assignment is not treated as achiral or racemic.")
    m["stereochemistry"]["geometry_status"] = f(
        geometry, geometry, evidence, method="inference",
        derivation="Status records the supplied drawing/text without inventing E/Z assignments.")
    m["stereochemistry"]["representation_status"] = f(
        representation, representation, evidence, method="inference",
        derivation="No complete stereochemical SMILES was promoted unless separately documented.")
    if assignment is not None:
        m["stereochemistry"]["assignment"] = f(assignment, assignment, evidence)
    if mixture is not None:
        m["stereochemistry"]["mixture"] = f(mixture, mixture, evidence)
    m["identity_note"] = f(
        identity_note or "Identity is retained by the paper label and source drawing; exact connectivity/stereochemical SMILES remains unencoded.",
        identity_note or "Identity is retained by the paper label and source drawing; exact connectivity/stereochemical SMILES remains unencoded.",
        evidence)
    molecules[mid] = m
    return m


def concrete_smiles(mid, label, name, smiles, evidence):
    m = base_molecule(mid, label, evidence, name=name, chirality="achiral",
                      geometry="not_applicable", representation="encoded",
                      identity_note="Connectivity candidate was derived from the explicit source name/drawing; parsing is not source verification.")
    m["representation_type"] = "smiles"
    m["smiles"] = f(name, smiles, evidence, method="structure_derivation",
                     derivation="Connectivity SMILES derived from the explicit source name/drawing and checked with the bundled RDKit diagnostic.")
    diag = st.stereo_inspect(smiles)
    m["stereochemistry"]["isomeric_smiles"] = f(
        name, diag["canonical_isomeric_smiles"], evidence, method="structure_derivation",
        derivation="Bundled standard stereo graph diagnostic; it does not independently verify source identity.")
    m["stereochemistry"]["identity_key"] = f(
        name, diag["stereo_graph_key"], evidence, method="structure_derivation",
        derivation="SHA-256 graph key from the canonical isomeric SMILES; not a sample identity.")
    return m


def generic_molecule(mid, label, fixed, variables, evidence, name=None):
    m = base_molecule(mid, label, evidence, name=name, chirality="unresolved",
                      geometry="unresolved", representation="not_encoded",
                      identity_note="The source uses a variable scaffold; it is retained as a generic template and not expanded to a concrete SMILES.")
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": f(fixed, {"description": fixed}, evidence),
        "variables": [],
        "instance_mapping": [],
        "reason": f("Variable/R-group source drawing", "Variable/R-group source drawing", evidence),
    }
    for variable in variables:
        m["generic_template"]["variables"].append({
            "label": variable,
            "attachment_points": f(f"Source attachment marked {variable}", f"Source attachment marked {variable}", evidence),
            "definition": missing(evidence),
            "allowed_values": missing(evidence),
            "constraints": missing(evidence),
        })
    molecules[mid] = m
    return m


# Explicitly named small reagents/solvent for which a conservative graph candidate is possible.
concrete_smiles("m_cpme", "CPME", "cyclopentyl methyl ether", "COC1CCCC1", [E[2], E[3], E[4]])
concrete_smiles("m_tfa", "TFA", "trifluoroacetic acid", "O=C(O)C(F)(F)F", [E[3], E[4]])
concrete_smiles("m_4brba", "4-BrPhCO2H", "4-bromobenzoic acid", "O=C(O)c1ccc(Br)cc1", [E[2]])

# Catalysts and salts. Their paper identities are preserved, but coordination/absolute
# stereochemistry is not forced into unsupported SMILES.
for label, name, ev, stereo in [
    ("P1", "prolinamide P1", [E[2]], "specified"),
    ("P2", "prolinamide P2", [E[2]], "specified"),
    ("P3", "prolinamide P3", [E[2]], "specified"),
    ("P4", "prolinamide P4", [E[2]], "specified"),
    ("P5", "prolinamide P5", [E[2]], "specified"),
    ("P6", "2-naphthylprolinamide P6", [E[2], E[3], E[4]], "specified"),
    ("Me-P6", "N-methylated derivative of P6", [E[4]], "specified"),
]:
    base_molecule(f"m_{label.lower().replace('-', '_')}", label, ev, name=name,
                  chirality=stereo, geometry="not_applicable", representation="not_encoded",
                  assignment={"source": "wedge-bonded prolinamide stereocenter shown; absolute descriptor not printed"})

for label, name, ev in [
    ("Cu(OTf)2", "copper(II) triflate", [E[1], E[2]]),
    ("CuI", "copper(I) iodide", [E[2]]),
    ("Cu(OAc)2", "copper(II) acetate", [E[2]]),
    ("Cu(MeCN)4PF6", "tetrakis(acetonitrile)copper(I) hexafluorophosphate", [E[2], E[3], E[4]]),
    ("Cu(MeCN)4BF4", "tetrakis(acetonitrile)copper(I) tetrafluoroborate", [E[2]]),
    ("Cu(acac)2", "copper(II) acetylacetonate", [E[2]]),
]:
    base_molecule("m_" + label.lower().replace("(", "_").replace(")", "").replace("/", "_").replace(" ", "").replace("+", ""),
                  label, ev, name=name, chirality="unresolved", geometry="not_applicable",
                  representation="unresolved")

CU_IDS = {
    "Cu(OTf)2": "m_cu_otf2", "CuI": "m_cui", "Cu(OAc)2": "m_cu_oac2",
    "Cu(MeCN)4PF6": "m_cu_mecn4pf6", "Cu(MeCN)4BF4": "m_cu_mecn4bf4",
    "Cu(acac)2": "m_cu_acac2",
}

# Paper-labelled current-work substrates/products.
base_molecule("m_1a", "1a", [E[1], E[2], E[3], E[4]],
              name="N-tethered allene-linked cyclohexanone 1a",
              chirality="unspecified", geometry="unspecified", representation="not_encoded")
base_molecule("m_2a", "2a", [E[1], E[2], E[3], E[4]],
              name="enantioenriched 4-vinyl morphan 2a", chirality="specified",
              geometry="unresolved", representation="not_encoded",
              assignment={"source": "bicyclic product drawing uses wedge bonds; absolute R/S not printed for 2a"})
base_molecule("m_2a_rac", "2a-rac", [E[4]], name="racemic morphan product 2a-rac",
              chirality="racemic", geometry="unresolved", representation="not_encoded",
              assignment={"source": "product explicitly labelled 2a-rac"}, mixture={"type": "racemate"})
base_molecule("m_5a", "5a", [E[4]], name="alkyne control substrate 5a",
              chirality="unspecified", geometry="not_applicable", representation="not_encoded")

for letter in "bcdefghijklmnop":
    ev = [E[3]]
    if letter in "gj":
        ev = [E[2], E[3]]
    base_molecule(f"m_1{letter}", f"1{letter}", ev,
                  name=f"N-linked allene substrate 1{letter}", chirality="unspecified",
                  geometry="unspecified", representation="not_encoded")
    assignment = {"source": "Scheme 2 wedge/dash product drawing"}
    if letter == "g":
        assignment = {"absolute_configuration": "1S,4R,5S", "method": "single-crystal X-ray diffraction", "CCDC": "2000237"}
    base_molecule(f"m_2{letter}", f"2{letter}", ev,
                  name=f"N-linked 4-vinyl morphan product 2{letter}", chirality="specified",
                  geometry="unresolved", representation="not_encoded", assignment=assignment)

for letter in "abcdefghijkl":
    base_molecule(f"m_3{letter}", f"3{letter}", [E[3]],
                  name=f"O-linked allene substrate 3{letter}", chirality="unspecified",
                  geometry="unspecified", representation="not_encoded")
    base_molecule(f"m_4{letter}", f"4{letter}", [E[3]],
                  name=f"4-vinyl oxamorphan product 4{letter}", chirality="specified",
                  geometry="unresolved", representation="not_encoded",
                  assignment={"source": "Scheme 2 wedge/dash product drawing"})

# Generic/literature-scheme identities.
generic_molecule("m_lit_dixon_s", "alkyne-tethered cyclohexanone", "generic alkyne-tethered cyclohexanone", ["R", "R′"], [E[2]])
generic_molecule("m_lit_dixon_p", "morphan core (Dixon 2017)", "generic bridged morphan product with exocyclic alkene", ["R", "R′"], [E[2]])
generic_molecule("m_lit_dong_s", "Dong alkyne-linked substrate", "generic alkyne-linked cyclohexanone", ["R", "R′", "Ar"], [E[2]])
generic_molecule("m_lit_dong_endo", "endo [3.2.1] morphan analogue", "generic endo morphan product", ["R", "R′", "Ar"], [E[2]])
generic_molecule("m_lit_dong_exo", "exo [3.2.1] morphan analogue", "generic exo morphan product", ["R", "R′", "Ar"], [E[2]])
generic_molecule("m_generic_s", "allene-tethered cyclohexanone", "allene-tethered cyclohexanone generic scaffold", ["R", "X"], [E[2], E[3]])
generic_molecule("m_generic_p", "4-vinyl (oxa)morphan core", "4-vinyl bridged (oxa)morphan generic scaffold", ["R", "X"], [E[2], E[3]])
for mid, label, name in [
    ("m_en", "En", "generic enamine organocatalyst"), ("m_ag", "Ag", "generic silver catalyst"),
    ("m_cu_generic", "Cu", "generic copper catalyst"), ("m_pd", "Pd", "generic palladium catalyst"),
    ("m_dtbmpp", "DTBMPP", "DTBMPP ligand"), ("m_diop", "(−)-DIOP", "(−)-DIOP ligand"),
]:
    base_molecule(mid, label, [E[2]], name=name, chirality="unresolved",
                  geometry="not_applicable", representation="unresolved")


def set_fact(obj, key, raw, norm, ev, method="original_extraction", derivation=None, status=None):
    obj[key] = f(raw, norm, ev, method=method, derivation=derivation, status=status)


def reaction(rid, ev, record_type, categories, location, name, transformation,
             general=(), unresolved=()):
    ev = list(dict.fromkeys(ev))
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, ev)
    r["location"] = f(location, location, ev)
    r["general_scheme_ids"] = list(general)
    r["reaction_name"] = f(name, name, ev)
    r["transformation_type"] = f(
        transformation, transformation, ev, method="inference",
        derivation="Transformation label derived from the depicted reactant/product relationship and article wording.")
    r["review"] = {
        "status": "needs_review", "method": "model", "reviewer": None,
        "checked_fields": [],
        "unresolved": list(unresolved) or ["Exact procedures and workup are unavailable because no SI was supplied."],
    }
    r["completeness"] = "partial"
    return r


def participant(r, pid, mid, role, ev, charges=None, eqref=None, product_role=None,
                stage="reaction"):
    p = st.new_participant(pid, mid, role, ev, stage=stage)
    if product_role is not None:
        p["product_role"] = f(product_role, product_role, ev)
    for key, value in (charges or {}).items():
        raw = value.get("raw", value.get("normalized"))
        p["charges"][key] = f(raw, value.get("normalized"), ev,
                               method=value.get("method", "original_extraction"),
                               derivation=value.get("derivation"))
    if eqref is not None:
        p["equivalence_reference"] = f(eqref, eqref, ev)
    r["participants"].append(p)
    return p


def one_step(r, ev, temperature=None, time=None, atmosphere=None, operation="reaction under the depicted conditions"):
    ids = [p["participant_id"] for p in r["participants"] if p["role"] != "product"]
    s = st.new_step("reaction", ev)
    s["participant_ids"] = ids
    s["operation"] = f(operation, operation, ev)
    if temperature is not None:
        s["temperature"] = f(temperature[0], temperature[1], ev)
    if time is not None:
        s["time"] = f(time[0], time[1], ev)
    if atmosphere is not None:
        s["atmosphere"] = f(atmosphere, atmosphere, ev)
    r["steps"].append(s)


def solvent(r, ev, composition=None, concentration=None, volume=None):
    if composition is not None:
        r["solvent_conditions"]["composition"] = f(composition, [composition], ev)
    if concentration is not None:
        r["solvent_conditions"]["concentration"] = f(concentration[0], concentration[1], ev)
        r["solvent_conditions"]["concentration_basis"] = f("substrate 1a or depicted scope substrate", "substrate", ev)
    if volume is not None:
        r["solvent_conditions"]["total_volume"] = f(volume[0], volume[1], ev)


def procedure(r, ev, text, general_ref=None):
    r["procedure"]["original_text"] = f(text, text, ev)
    if general_ref:
        r["procedure"]["general_procedure_reference"] = f(general_ref, general_ref, ev)


def outcome(r, oid, kind, targets, ev, raw, normalized, product_role=None,
            yield_type=None, basis=None, scope=None, value_fact=None):
    o = st.new_outcome(oid, kind, targets, ev)
    o["value"] = value_fact if value_fact is not None else f(raw, normalized, ev)
    if product_role is not None:
        o["product_role"] = f(product_role, product_role, ev)
    if yield_type is not None:
        o["yield_type"] = f(yield_type, yield_type, ev)
    if basis is not None:
        o["measurement_basis"] = f(basis, basis, ev)
    if scope is not None:
        o["scope"] = f(scope, scope, ev)
    r["outcomes"].append(o)
    return o


def formation_classification(r, ev, include_h=True, limitation=None):
    changes = [{
        "change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
        "before_order": 0, "after_order": 1, "change_type": "formed",
    }]
    labels = [{"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1"]}]
    if include_h:
        changes.append({
            "change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}],
            "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1,
        })
        labels.append({"category_number": 1, "code": "C-H", "supporting_change_ids": ["b2"]})
    c = r["classification"]
    c["mapped_reaction_smiles"] = missing(ev)
    c["bond_changes"] = [f("source-drawing net change", x, ev, method="structure_derivation",
                               derivation="Partial net-bond interpretation from the allene cyclization drawing; exact atom mapping is unavailable.") for x in changes]
    c["labels"] = f("drawing-derived partial labels", labels, ev, method="structure_derivation",
                    derivation="Taxonomy labels correspond only to the explicitly retained partial bond changes.")
    c["coverage"] = f("partial", "partial", ev, method="inference",
                      derivation="Complex source structures were not atom-mapped; allene bond-order changes and atom correspondence remain unresolved.")
    c["limitations"] = [
        limitation or "Partial drawing-derived classification; no complete atom mapping was asserted.",
        "The C-H gain is a net transformation and does not identify the elementary proton source.",
    ]
    c["model_votes"] = []


def no_reaction_classification(r, ev, reason):
    c = r["classification"]
    c["mapped_reaction_smiles"] = missing(ev)
    c["bond_changes"] = []
    c["labels"] = f(reason, [], ev, method="inference", derivation="No observed transformation is not assigned the intended bond-formation class.")
    c["coverage"] = f(reason, "not_applicable", ev, method="inference", derivation="No observed transformation.")
    c["limitations"] = [reason]
    c["model_votes"] = []


def unresolved_classification(r, ev, reason):
    c = r["classification"]
    c["mapped_reaction_smiles"] = missing(ev)
    c["bond_changes"] = []
    c["labels"] = f(reason, [], ev, method="inference", derivation=reason)
    c["coverage"] = f(reason, "unresolved", ev, method="inference", derivation=reason)
    c["limitations"] = [reason]
    c["model_votes"] = []


reactions = []

# Scheme 1: cited literature and a current-paper generic design.  The two Dong
# arrows are kept separate because they depict different ligand/product outcomes.
literature_specs = [
    ("scheme1_dixon_2017", "m_lit_dixon_s", "m_lit_dixon_p", [("en", "m_en", "catalyst"), ("ag", "m_ag", "catalyst")],
     "Dixon 2017 enantioselective alkyne-linked cyclohexanone desymmetrization"),
    ("scheme1_dong_endo", "m_lit_dong_s", "m_lit_dong_endo", [("pd", "m_pd", "catalyst"), ("lig", "m_dtbmpp", "ligand")],
     "Dong 2019 endo-selective, diastereoselective racemic annulation"),
    ("scheme1_dong_exo", "m_lit_dong_s", "m_lit_dong_exo", [("pd", "m_pd", "catalyst"), ("lig", "m_diop", "ligand")],
     "Dong 2019 exo-selective, diastereoselective racemic annulation"),
]
for rid, smid, pmid, agents, title in literature_specs:
    r = reaction(rid, [E[2]], "literature_example", ["other"], "Scheme 1 (cited prior work)", title,
                 "cited intramolecular desymmetrization/annulation",
                 unresolved=["Scheme 1 is a citation-level summary; full experimental conditions belong to the cited paper and were not imported."])
    participant(r, "substrate", smid, "reactant", [E[2]])
    participant(r, "product", pmid, "product", [E[2]], product_role="depicted cited product")
    for pid, mid, role in agents:
        participant(r, pid, mid, role, [E[2]])
    unresolved_classification(r, [E[2]], "Citation-level generic structures do not support a complete net-bond derivation here.")
    reactions.append(r)

GENERIC_ID = f"{DOC}:scheme1_current_generic"
r = reaction("scheme1_current_generic", [E[2], E[3]], "generic_scheme", ["other"],
             "Scheme 1 this work and Scheme 2 general transformation",
             "dual catalytic enantioselective desymmetrization of allene-tethered cyclohexanones",
             "intramolecular cycloisomerization/desymmetrization to a 4-vinyl (oxa)morphan",
             unresolved=["Generic R/X structures are templates; no single concrete atom mapping applies."])
participant(r, "generic_substrate", "m_generic_s", "reactant", [E[2], E[3]])
participant(r, "generic_product", "m_generic_p", "product", [E[2], E[3]], product_role="generic product")
participant(r, "en", "m_en", "catalyst", [E[2]])
participant(r, "cu", "m_cu_generic", "catalyst", [E[2]])
formation_classification(r, [E[2], E[3]], limitation="Generic scheme supports only the net C-C closure and protonation-derived C-H gain.")
reactions.append(r)


def q(number, unit, raw=None, qualifier=None):
    d = {"number": number, "unit": unit}
    if raw is not None:
        d["raw_text"] = raw
    if qualifier is not None:
        d["qualifier"] = qualifier
    return d


table_rows = [
    (1, "Cu(OTf)2", "P1", 76, "8:1", "82:18"),
    (2, "CuI", "P1", 44, ">20:1", "72:28"),
    (3, "Cu(OAc)2", "P1", 40, "7:1", "70:30"),
    (4, "Cu(MeCN)4PF6", "P1", 98, "12:1", "87:13"),
    (5, "Cu(MeCN)4BF4", "P1", 76, "11:1", "78:22"),
    (6, "Cu(acac)2", "P1", 93, "10:1", "80:20"),
    (7, "Cu(MeCN)4PF6", "P2", 99, ">20:1", "82:18"),
    (8, "Cu(MeCN)4PF6", "P3", 90, "18:1", "80:20"),
    (9, "Cu(MeCN)4PF6", "P4", 99, ">20:1", "89.5:10.5"),
    (10, "Cu(MeCN)4PF6", "P5", 90, "18:1", "89:11"),
    (11, "Cu(MeCN)4PF6", "P4", None, None, None),
    (12, "Cu(MeCN)4PF6", "P4", 91, ">20:1", "90:10"),
    (13, "Cu(MeCN)4PF6", "P4", 93, ">20:1", "91:9"),
    (14, "Cu(MeCN)4PF6", "P6", 82, ">20:1", "92.5:7.5"),
    (15, "Cu(MeCN)4PF6", "P6", 81, ">20:1", "96:4"),
]


def ratio_value(text, metric):
    qualifier = "gt" if text.startswith(">") else None
    body = text.lstrip(">")
    a, b = body.split(":")
    return {"metric": metric, "ratio": [float(a), float(b)], "order": "major:minor", **({"qualifier": qualifier} if qualifier else {})}


for entry, cu_label, pro, yld, dr, er in table_rows:
    ev = [E[2]]
    if entry == 1:
        ev = [E[1], E[2]]
    if entry == 15:
        ev = [E[2], E[3]]  # same optimized experiment is repeated as Scheme 2 product 2a
    categories = ["optimization"] + (["substrate_scope"] if entry == 15 else [])
    location = f"Table 1, entry {entry}" + ("; repeated as Scheme 2, 2a" if entry == 15 else "")
    r = reaction(f"table1_entry_{entry:02d}", ev, "specific_experiment", categories, location,
                 "dual catalytic enantioselective desymmetrization",
                 "allene cycloisomerization/desymmetrization to 4-vinyl morphan 2a",
                 general=[GENERIC_ID],
                 unresolved=["No SI was supplied; workup/purification and full operation order remain unavailable."])
    participant(r, "substrate_1a", "m_1a", "reactant", ev,
                charges={"amount": {"normalized": q(0.2, "mmol"), "raw": "1a (0.2 mmol)"}})
    participant(r, "product_2a", "m_2a", "product", ev, product_role="major bicyclic product")
    participant(r, "cu", CU_IDS[cu_label], "catalyst", ev,
                charges={"amount": {"normalized": q(0.02, "mmol"), "raw": "[Cu] catalyst (0.02 mmol)"},
                         "mol_percent": {"normalized": q(10, "mol%"), "raw": "10 mol%"}}, eqref="m_1a")
    p_mid = "m_" + pro.lower()
    participant(r, "prolinamide", p_mid, "catalyst", ev,
                charges={"amount": {"normalized": q(0.06, "mmol"), "raw": "prolinamide catalyst (0.06 mmol)"},
                         "mol_percent": {"normalized": q(30, "mol%"), "raw": "30 mol%"}}, eqref="m_1a")
    acid_mid = "m_tfa" if entry == 15 else "m_4brba"
    acid_raw = "TFA (50 mol%) used instead of 4-BrPhCO2H" if entry == 15 else "4-bromobenzoic acid (0.1 mmol, 50 mol%)"
    acid_charges = {"mol_percent": {"normalized": q(50, "mol%"), "raw": "50 mol%"}}
    if entry != 15:
        acid_charges["amount"] = {"normalized": q(0.1, "mmol"), "raw": "0.1 mmol"}
    participant(r, "acid", acid_mid, "additive", ev, charges=acid_charges, eqref="m_1a")
    solvent_charge = {} if entry >= 12 else {"volume": {"normalized": q(2, "mL"), "raw": "CPME (2 mL)"}}
    participant(r, "cpme", "m_cpme", "solvent", ev, charges=solvent_charge)
    concentration = 0.04 if entry == 12 else 0.02 if entry >= 13 else 0.1
    volume = None if entry >= 12 else ("2 mL", q(2, "mL"))
    solvent(r, ev, "CPME", (f"{concentration:g} M", q(concentration, "M")), volume)
    r["reaction_scale"] = f("1a (0.2 mmol)", {"reference_molecule_id": "m_1a", "amount": q(0.2, "mmol")}, ev)
    temperature = 100 if entry == 11 else 120
    one_step(r, ev, (f"{temperature} °C", q(temperature, "°C")), ("48 h", q(48, "h")), "Ar")
    procedure(r, ev,
              "General conditions: 1a (0.2 mmol), [Cu] catalyst (0.02 mmol, 10 mol%), 4-bromobenzoic acid (0.1 mmol, 50 mol%), prolinamide catalyst (0.06 mmol, 30 mol%), in CPME (0.1 M, 2 mL) at 120 °C under Ar for 48 h; listed entries apply stated changes.",
              "Table 1 general conditions and entry-specific change")
    if entry == 11:
        outcome(r, "no_reaction", "no_reaction", ["m_2a"], ev, "—; nd; nd", True,
                product_role="intended product not observed", scope={"entry": 11, "temperature": "100 °C"})
        no_reaction_classification(r, ev, "Entry 11 reports extinguished reactivity at 100 °C.")
    else:
        outcome(r, "yield_2a", "yield", ["m_2a"], ev, f"{yld}%", q(yld, "%"),
                product_role="major product", yield_type="isolated yield", basis="isolated")
        outcome(r, "dr_2a", "selectivity", ["m_2a"], ev, dr, ratio_value(dr, "dr"),
                product_role="major product", basis="1H NMR analysis of the crude reaction mixture",
                scope={"metric": "dr", "order": "major:minor"})
        outcome(r, "er_2a", "selectivity", ["m_2a"], ev, er, ratio_value(er, "er"),
                product_role="major product", basis="chiral HPLC analysis of the pure product",
                scope={"metric": "er", "order": "major:minor", "absolute configuration of major enantiomer": "unassigned for 2a"})
        formation_classification(r, ev)
    reactions.append(r)


scope_n = {
    "b": (75, ">20:1", "96:4"), "c": (88, ">20:1", "94:6"),
    "d": (83, ">20:1", "94:6"), "e": (84, ">20:1", "94:6"),
    "f": (90, ">20:1", "95:5"), "g": (84, ">20:1", "94:6"),
    "h": (86, "18:1", "95:5"), "i": (85, "9:1", "96:4"),
    "j": (86, ">20:1", "96:4"), "k": (85, "10:1", "91:9"),
    "l": (80, ">20:1", "95:5"), "m": (83, "18:1", "94:6"),
    "n": ("trace", None, None), "o": ("trace", None, None),
    "p": (89, ">20:1", "91:9"),
}
scope_o = {
    "a": ("trace", None, None), "b": (97, ">20:1", "92:8"),
    "c": (98, ">20:1", "90:10"), "d": (84, ">20:1", "93:7"),
    "e": (94, ">20:1", "94:6"), "f": (91, ">20:1", "94:6"),
    "g": (67, ">20:1", "81:19"), "h": (89, ">20:1", "93:7"),
    "i": (91, ">20:1", "94:6"), "j": (88, ">20:1", "93:7"),
    "k": (81, ">20:1", "92:8"), "l": (85, ">20:1", "94:6"),
}


def scope_record(series, letter, values):
    substrate_prefix, product_prefix = ("1", "2") if series == "N" else ("3", "4")
    ev = [E[3]]
    if series == "N" and letter in "gj":
        ev = [E[2], E[3]]
    smid, pmid = f"m_{substrate_prefix}{letter}", f"m_{product_prefix}{letter}"
    product_label = f"{product_prefix}{letter}"
    unresolved = [
        "No SI was supplied; substrate charge, workup, purification, and yield/selectivity measurement basis are unavailable.",
        "Complex source structures retain labels/drawings without promoted exact SMILES or complete atom mapping.",
    ]
    if series == "N" and letter == "g":
        unresolved.append("Yield conflict retained: Scheme 2 prints 84%, while page-2 prose states 78%.")
    r = reaction(f"scheme2_{product_label}", ev, "specific_experiment", ["substrate_scope"],
                 f"Scheme 2, product {product_label}",
                 "dual catalytic enantioselective desymmetrization",
                 "allene cycloisomerization/desymmetrization to a 4-vinyl (oxa)morphan",
                 general=[GENERIC_ID], unresolved=unresolved)
    participant(r, "substrate", smid, "reactant", ev)
    participant(r, "product", pmid, "product", ev, product_role="depicted bicyclic product")
    participant(r, "p6", "m_p6", "catalyst", ev,
                charges={"mol_percent": {"normalized": q(30, "mol%"), "raw": "P6 (30 mol%)"}}, eqref=smid)
    participant(r, "cu", "m_cu_mecn4pf6", "catalyst", ev,
                charges={"mol_percent": {"normalized": q(10, "mol%"), "raw": "Cu(MeCN)4PF6 (10 mol%)"}}, eqref=smid)
    participant(r, "tfa", "m_tfa", "additive", ev,
                charges={"mol_percent": {"normalized": q(50, "mol%"), "raw": "TFA (50 mol%)"}}, eqref=smid)
    participant(r, "cpme", "m_cpme", "solvent", ev)
    solvent(r, ev, "CPME", ("0.02 M", q(0.02, "M")))
    one_step(r, ev, ("120 °C", q(120, "°C")), ("48 h", q(48, "h")), "Ar")
    procedure(r, ev, "Scheme 2 arrow: P6 (30 mol%), Cu(MeCN)4PF6 (10 mol%), TFA (50 mol%), CPME (0.02 M), 120 °C, 48 h, Ar.",
              "Scheme 2 general arrow")
    yld, dr, er = values
    if yld == "trace":
        outcome(r, "trace_product", "trace", [pmid], ev, "trace", "trace",
                product_role="trace depicted product", scope={"metric": "product formation"})
    else:
        value_fact = None
        if series == "N" and letter == "g":
            alt_scheme = f("Scheme 2: 84%", q(84, "%"), [E[3]])
            alt_prose = f("page-2 prose: 78%", q(78, "%"), [E[2]])
            value_fact = f(["Scheme 2: 84%", "page-2 prose: 78%"], None, [E[2], E[3]],
                           status="conflict", alternatives=[alt_scheme, alt_prose])
        outcome(r, "yield_product", "yield", [pmid], ev, f"{yld}%", q(yld, "%"),
                product_role="depicted product", yield_type="reported yield", value_fact=value_fact)
        outcome(r, "dr_product", "selectivity", [pmid], ev, dr, ratio_value(dr, "dr"),
                product_role="depicted product", scope={"metric": "dr", "order": "major:minor"})
        er_scope = {"metric": "er", "order": "major:minor", "major absolute configuration": "unassigned"}
        if series == "N" and letter == "g":
            er_scope["major absolute configuration"] = "1S,4R,5S"
        outcome(r, "er_product", "selectivity", [pmid], ev, er, ratio_value(er, "er"),
                product_role="depicted product", scope=er_scope)
    formation_classification(r, ev)
    return r


for letter, vals in scope_n.items():
    reactions.append(scope_record("N", letter, vals))
for letter, vals in scope_o.items():
    reactions.append(scope_record("O", letter, vals))


# Scheme 3 controls.
def add_control_common(r, substrate_mid, product_mid=None, p6=True, copper=True, methylated=False):
    ev = [E[4]]
    participant(r, "substrate", substrate_mid, "reactant", ev)
    if product_mid:
        participant(r, "product", product_mid, "product", ev, product_role="depicted/target product")
    if p6:
        participant(r, "prolinamide", "m_me_p6" if methylated else "m_p6", "catalyst", ev,
                    charges={"mol_percent": {"normalized": q(30, "mol%"), "raw": "30 mol%"}}, eqref=substrate_mid)
    if copper:
        participant(r, "cu", "m_cu_mecn4pf6", "catalyst", ev,
                    charges={"mol_percent": {"normalized": q(10, "mol%"), "raw": "[Cu] cat. (10 mol%)"}}, eqref=substrate_mid)
    participant(r, "tfa", "m_tfa", "additive", ev,
                charges={"mol_percent": {"normalized": q(50, "mol%"), "raw": "TFA (50 mol%)"}}, eqref=substrate_mid)
    participant(r, "cpme", "m_cpme", "solvent", ev)
    solvent(r, ev, "CPME", ("0.02 M", q(0.02, "M")))
    one_step(r, ev, ("120 °C", q(120, "°C")), ("48 h", q(48, "h")), "Ar")


r = reaction("scheme3_A_alkyne_control", [E[4]], "specific_experiment", ["control", "mechanistic_study"],
             "Scheme 3A", "alkyne-intermediate control experiment", "test of an alkyne-mediated pathway",
             unresolved=["No SI was supplied; substrate amount, workup, purification and detection limit are unavailable."])
add_control_common(r, "m_5a", "m_2a")
procedure(r, [E[4]], "Alkyne derivative 5a was subjected to the standard conditions; 2a was not observed.", "standard conditions")
outcome(r, "no_cyclized_product", "no_reaction", ["m_2a"], [E[4]], "2a (not observed)", True,
        product_role="intended cyclized product", scope={"test": "alkyne intermediate"})
no_reaction_classification(r, [E[4]], "No cyclized product was observed from alkyne 5a.")
reactions.append(r)

r = reaction("scheme3_B_methylated_catalyst", [E[4]], "specific_experiment", ["control", "mechanistic_study"],
             "Scheme 3B", "methylated-prolinamide deletion control", "test of the prolinamide N-H requirement",
             general=[GENERIC_ID], unresolved=["No SI was supplied; workup and measurement methods are unavailable."])
add_control_common(r, "m_1a", "m_2a", methylated=True)
procedure(r, [E[4]], "Methylated derivative of the optimized catalyst (30 mol%) under standard reaction conditions.", "standard conditions")
outcome(r, "conversion_1a", "conversion", ["m_1a"], [E[4]], "<5% conversion",
        {"number": 5, "unit": "%", "qualifier": "lt", "raw_text": "<5% conversion"},
        scope={"metric": "conversion of 1a"})
outcome(r, "recovery_1a", "recovery", ["m_1a"], [E[4]], "95% SM recovered", q(95, "%"),
        scope={"metric": "starting material recovery"})
unresolved_classification(r, [E[4]], "Less than 5% conversion is reported, but the identity/amount of any formed product is not established.")
reactions.append(r)

r = reaction("scheme3_C_no_copper", [E[4]], "specific_experiment", ["control", "mechanistic_study"],
             "Scheme 3C", "copper deletion control", "test of the copper-catalyst requirement",
             general=[GENERIC_ID], unresolved=[
                 "Figure prints 2a: <5%, while prose says no formation was observed; both statements are retained as a source conflict.",
                 "No SI was supplied; measurement basis and detection limit are unavailable."])
add_control_common(r, "m_1a", "m_2a", copper=False)
procedure(r, [E[4]], "Standard reaction conditions without copper catalyst.", "standard conditions with copper omitted")
alt_fig = f("Scheme 3C: 2a <5%", {"number": 5, "unit": "%", "qualifier": "lt", "raw_text": "<5%"}, [E[4]])
alt_text = f("prose: no formation of the cyclized product was observed", "not observed", [T[4]])
conflict = f(["2a: <5%", "no formation observed"], None, [E[4], T[4]], status="conflict", alternatives=[alt_fig, alt_text])
outcome(r, "product_2a_source_conflict", "yield", ["m_2a"], [E[4], T[4]], None, None,
        product_role="target product", yield_type="source expression", value_fact=conflict)
no_reaction_classification(r, [E[4], T[4]], "Prose reports no observed cyclized product; the figure's <5% statement is retained as an unresolved conflict.")
reactions.append(r)

r = reaction("scheme3_D_no_prolinamide", [E[4]], "specific_experiment", ["control", "mechanistic_study"],
             "Scheme 3D", "prolinamide deletion control", "background cyclization without prolinamide",
             general=[GENERIC_ID], unresolved=["Yield/selectivity measurement bases, exact substrate charge and workup are unreported in the supplied main text."])
add_control_common(r, "m_1a", "m_2a_rac", p6=False)
procedure(r, [E[4]], "Standard reaction conditions without prolinamide catalyst.", "standard conditions with prolinamide omitted")
outcome(r, "yield_2a_rac", "yield", ["m_2a_rac"], [E[4]], "16%", q(16, "%"),
        product_role="racemic bicyclic product", yield_type="reported yield")
outcome(r, "dr_2a_rac", "selectivity", ["m_2a_rac"], [E[4]], "1:1 dr", ratio_value("1:1", "dr"),
        product_role="racemic bicyclic product", scope={"metric": "dr", "order": "diastereomer 1:diastereomer 2"})
formation_classification(r, [E[4]])
reactions.append(r)


# Computed mechanism is not an experiment and remains explicitly proposed.
r = reaction("scheme4_proposed_mechanism", [E[4], E[5]], "proposed_mechanism", ["other"],
             "Scheme 4 and accompanying DFT discussion", "computed reaction pathway for the dual catalytic cyclization",
             "proposed enamine attack on a copper-coordinated allene followed by protonation and hydrolysis",
             general=[GENERIC_ID], unresolved=["Computed intermediates/transition states are not experimentally observed products."])
participant(r, "substrate", "m_1a", "reactant", [E[4], E[5]])
participant(r, "product", "m_2a", "product", [E[4], E[5]], product_role="computed pathway product E1/2a scaffold")
participant(r, "p6", "m_p6", "catalyst", [E[4], E[5]])
participant(r, "cu", "m_cu_mecn4pf6", "catalyst", [E[4], E[5]])
participant(r, "tfa", "m_tfa", "additive", [E[4], E[5]])
one_step(r, [E[4], E[5]], operation="computed pathway, not an experimental operation")
procedure(r, [E[4], E[5]], "DFT pathway: enamine formation; copper-allene complexation; C-C bond formation; protonation; enamine hydrolysis.")
formation_classification(r, [E[4], E[5]], limitation="Classification describes the proposed net transformation, not experimental proof of the elementary pathway.")
reactions.append(r)


bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
bundle["source_context"] = [
    {
        "context_id": f"{DOC}:metadata",
        "kind": "bibliographic_metadata",
        "evidence_ids": [E[1]],
        "title": "Dual catalytic enantioselective desymmetrization of allene-tethered cyclohexanones",
        "doi": "10.1039/d0sc02878a",
        "citation": "Chem. Sci., 2020, 11, 7444-7450",
    },
    {
        "context_id": f"{DOC}:xray_2g",
        "kind": "structure_assignment",
        "evidence_ids": [E[2], E[3], E[6]],
        "observation": "2g assigned 1S,4R,5S by single-crystal X-ray diffraction; CCDC 2000237.",
    },
    {
        "context_id": f"{DOC}:dft",
        "kind": "computational_mechanism",
        "evidence_ids": [E[4], E[5]],
        "observation": "Scheme 4 is a computed energy profile, not a synthetic experiment; TSB1-C1 is favored by 2.2 kcal mol-1 over TSB2-C2 and the reported activation strains are 27.1 and 29.0 kcal mol-1.",
    },
    {
        "context_id": f"{DOC}:si_footnote_16",
        "kind": "unresolved_si_reference",
        "evidence_ids": [E[6]],
        "observation": "Footnote 16 says attempted preparation of allene S6b under standard conditions gave competitive cycloisomerized S6c, isolated as a racemate; Scheme S2 and the actual procedure are in SI that was not supplied, so no complete experiment record was fabricated.",
    },
]


def save(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n")


save(RUN / "submitted.json", bundle)

record_ids = [r["reaction_id"] for r in reactions]
by_prefix = lambda prefix: [x for x in record_ids if prefix in x]
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 1, "state": "attempted",
            "regions": [
                {"region": "title, abstract, bibliographic block", "finding": "Paper title/DOI and an explicit statement that ESI exists; no SI file was supplied or associated.", "evidence_ids": [E[1]]},
                {"region": "Results & discussion, lower right", "finding": "Model reaction/Table 1 entry 1 is described: 1a, P1, Cu(OTf)2, 4-bromobenzoic acid, CPME, 120 °C; 2a 76%, 82:18 er.", "evidence_ids": [E[1], T[1]]},
            ],
            "record_ids": [f"{DOC}:table1_entry_01"],
            "unresolved": ["Page 1 points to ESI and CCDC data, but no SI/CIF was provided."],
        },
        {
            "page": 2, "state": "attempted",
            "regions": [
                {"region": "Scheme 1, left column", "finding": "Three citation-level prior-work arrows (Dixon 2017; Dong 2019 endo and exo) and one current-paper generic transformation were separated by record type.", "evidence_ids": [E[2]]},
                {"region": "Table 1, right column", "finding": "All 15 optimization rows and footnotes were read row by row; entry 15 is the same experiment repeated as Scheme 2 product 2a.", "evidence_ids": [E[2], T[2]]},
                {"region": "scope prose, bottom", "finding": "Cross-page scope text identifies 2g absolute configuration and creates a yield conflict: prose 78% versus Scheme 2 84%; other repeated scope values were linked, not duplicated.", "evidence_ids": [E[2], T[2]]},
            ],
            "record_ids": by_prefix("scheme1_") + by_prefix("table1_") + [f"{DOC}:scheme2_2g", f"{DOC}:scheme2_2j"],
            "unresolved": ["Full procedures/workups and scope measurement bases are not present in the supplied main text."],
        },
        {
            "page": 3, "state": "attempted",
            "regions": [
                {"region": "Scheme 2 general arrow", "finding": "General scope conditions read as P6 30 mol%, Cu(MeCN)4PF6 10 mol%, TFA 50 mol%, CPME 0.02 M, 120 °C, 48 h, Ar.", "evidence_ids": [E[3]]},
                {"region": "Scheme 2A, nitrogen-linked allenes", "finding": "Products 2a-2p were inventoried individually; 2a was merged with Table 1 entry 15, 2n/2o are trace, and no unshown 2q experiment was invented despite the range label.", "evidence_ids": [E[3]]},
                {"region": "Scheme 2B, oxygen-linked allenes", "finding": "Products 4a-4l were inventoried individually; 4a is trace and all printed yield/dr/er values were retained.", "evidence_ids": [E[3]]},
                {"region": "scope prose", "finding": "2p and oxygen-series summary statements were cross-checked against the figure; no separate duplicate experiments were created.", "evidence_ids": [E[3], T[3]]},
            ],
            "record_ids": [f"{DOC}:table1_entry_15"] + by_prefix("scheme2_"),
            "unresolved": ["Exact complex structures and most absolute configurations were not promoted to SMILES; yield/dr/er measurement bases are unreported on the supplied page."],
        },
        {
            "page": 4, "state": "attempted",
            "regions": [
                {"region": "Scheme 3A-D, lower left", "finding": "Four independent mechanistic control experiments were recorded: alkyne substrate, methylated catalyst, no copper, and no prolinamide.", "evidence_ids": [E[4]]},
                {"region": "control prose, upper left", "finding": "Figure/prose result differences were retained, including Scheme 3C 2a <5% versus prose no formation observed.", "evidence_ids": [E[4], T[4]]},
                {"region": "DFT methods/discussion, right column", "finding": "Computational settings and mechanistic interpretation are source context/proposed mechanism, not synthetic experiments.", "evidence_ids": [E[4], T[4]]},
            ],
            "record_ids": by_prefix("scheme3_") + [f"{DOC}:scheme4_proposed_mechanism"],
            "unresolved": ["Control-experiment yield/conversion measurement bases and exact detection limits are not reported."],
        },
        {
            "page": 5, "state": "attempted",
            "regions": [
                {"region": "Scheme 4 energy profile", "finding": "Computed intermediates/transition states and energies were treated as a proposed mechanism and source context, not observed reaction records.", "evidence_ids": [E[5]]},
                {"region": "mechanistic discussion and conclusion", "finding": "Authors propose enamine attack on copper-bound allene, protonation, and hydrolysis; no additional synthetic experiment is present.", "evidence_ids": [E[5], T[5]]},
            ],
            "record_ids": [f"{DOC}:scheme4_proposed_mechanism"],
            "unresolved": ["Computational structures are not experimental molecule identities and were not promoted as observed products."],
        },
        {
            "page": 6, "state": "attempted",
            "regions": [
                {"region": "Acknowledgements and references", "finding": "No independent synthetic experiment, table, mass-spectrum peak list, or reaction arrow appears on this page.", "evidence_ids": [E[6]]},
                {"region": "reference note 14", "finding": "CCDC 2000237 supports the source's 2g crystallographic assignment and points to unavailable ESI for numbering.", "evidence_ids": [E[6]]},
                {"region": "reference note 16", "finding": "S6b/S6c racemate is mentioned only through an SI Scheme S2 reference; because SI and the starting precursor/conditions are unavailable, it remains unresolved source context rather than a fabricated full experiment.", "evidence_ids": [E[6], T[6]]},
            ],
            "record_ids": [],
            "unresolved": ["The S6b/S6c procedure and yield require explicitly supplied SI; no SI association was inferred."],
        },
    ],
}
save(RUN / "coverage/pages_0001_0006.json", coverage)

inventory = {
    "document_id": DOC,
    "covered_pages": [1, 2, 3, 4, 5, 6],
    "molecule_count": len(bundle["molecules"]),
    "reaction_count": len(bundle["reactions"]),
    "reaction_types": {kind: sum(r["record_type"] == kind for r in reactions) for kind in sorted({r["record_type"] for r in reactions})},
    "notes": [
        "Raw recipe was derived anew from this paper; no prior workflow record was imported.",
        "Entry 15/2a was deduplicated as one experiment with two evidence pages.",
        "2g yield conflict and Scheme 3C figure/prose conflict are preserved.",
        "No SI was associated.",
    ],
}
save(RUN / "raw_model/page_inventory.json", inventory)

# Mark only the requested six pages attempted. Page 7 remains pending/context-only.
tasks_path = RUN / "tasks.json"
tasks = json.loads(tasks_path.read_text())
for task in tasks:
    page = int(task["evidence_id"].split(":page", 1)[1].split(":", 1)[0])
    task["state"] = "attempted" if 1 <= page <= 6 else "pending"
save(tasks_path, tasks)

# Literal text hits are review leads only; they are saved separately and were
# reconciled against the page images above rather than merged automatically.
dataset_for_mentions = json.loads((RUN / "dataset.json").read_text())
save(RUN / "raw_model/source_mentions.json",
     st.source_mentions(dataset_for_mentions, ["1a", "2a", "2g", "2p", "4b", "P6", "TFA"], context=180))

print(json.dumps(inventory, ensure_ascii=False, indent=2))
