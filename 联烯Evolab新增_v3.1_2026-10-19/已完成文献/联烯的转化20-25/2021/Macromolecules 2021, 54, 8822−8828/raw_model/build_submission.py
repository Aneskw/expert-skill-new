#!/usr/bin/env python3
"""Build the source-derived v3.1 submission for PDF pages 1--6.

This file records model extraction decisions. It does not alter the schema or
promote any record beyond needs_review.
"""
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_56ca04973542741c"
RUN = ROOT / "runs/current"


def ev(page):
    return [f"{DOC}:page{page}:text", f"{DOC}:page{page}:image"]


def f(raw=None, normalized=None, evidence=(), method="original_extraction", derivation=None,
      status=None, alternatives=()):
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation,
                   status=status, alternatives=alternatives)


def set_name(m, name, evidence, english=None, aliases=None, abbreviations=None):
    m["identity_status"] = "source_identified"
    m["names"]["original"] = f(name, name, evidence)
    m["names"]["english"] = f(english or name, english or name, evidence)
    if aliases:
        m["names"]["aliases"] = f(aliases, aliases, evidence)
    if abbreviations:
        m["names"]["abbreviations"] = f(abbreviations, abbreviations, evidence)
    m["cas_numbers"] = f()


def set_concrete(m, smiles, evidence, chirality="achiral", geometry="not_applicable",
                 stereo_note=None):
    diag = st.stereo_inspect(smiles)
    deriv = "Candidate connectivity/isomeric SMILES was derived from the cited source drawing or unambiguous source name; toolkit parsing does not establish source correctness."
    m["representation_type"] = "smiles"
    m["smiles"] = f("source structure/name", smiles, evidence, "structure_derivation", deriv)
    m["stereochemistry"]["chirality_status"] = f(None, chirality, evidence, "structure_derivation", deriv)
    m["stereochemistry"]["geometry_status"] = f(None, geometry, evidence, "structure_derivation", deriv)
    m["stereochemistry"]["representation_status"] = f(None, "encoded", evidence, "structure_derivation", deriv)
    m["stereochemistry"]["isomeric_smiles"] = f(None, diag["canonical_isomeric_smiles"], evidence, "structure_derivation", deriv)
    m["stereochemistry"]["identity_key"] = f(None, diag["stereo_graph_key"], evidence, "structure_derivation", deriv)
    if stereo_note:
        m["stereochemistry"]["assignment"] = f(stereo_note, stereo_note, evidence)
    m["identity_note"] = f(None, "Candidate structure; remains machine-derived and needs expert review.", evidence,
                           "inference", "Structure candidate was transcribed from the cited source region.")


def set_unresolved(m, evidence, note, chirality="unresolved", geometry="unresolved"):
    m["representation_type"] = "unresolved"
    m["stereochemistry"]["chirality_status"] = f(note, chirality, evidence)
    m["stereochemistry"]["geometry_status"] = f(note, geometry, evidence)
    m["stereochemistry"]["representation_status"] = f(note, "not_encoded", evidence)
    m["identity_note"] = f(note, note, evidence)


def set_template(m, evidence, fixed, reason, chirality="unspecified", geometry="unspecified",
                 assignment=None):
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": f(fixed, {"description": fixed}, evidence),
        "variables": [],
        "instance_mapping": [],
        "reason": f(reason, reason, evidence),
    }
    m["stereochemistry"]["chirality_status"] = f(reason, chirality, evidence)
    m["stereochemistry"]["geometry_status"] = f(reason, geometry, evidence)
    m["stereochemistry"]["representation_status"] = f(reason, "not_encoded", evidence)
    if assignment:
        m["stereochemistry"]["assignment"] = f(assignment, assignment, evidence)
    m["identity_note"] = f(reason, reason, evidence)


def molecule(mid, label, name, page=2):
    evidence = ev(page)
    m = st.new_molecule(mid, DOC, evidence)
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    set_name(m, name, evidence)
    return m


M = {}


def add_molecule(mid, label, name, page=2, smiles=None, template=None, unresolved=None,
                 chirality="achiral", geometry="not_applicable", assignment=None):
    m = molecule(mid, label, name, page)
    evidence = ev(page)
    if smiles:
        set_concrete(m, smiles, evidence, chirality, geometry, assignment)
    elif template:
        set_template(m, evidence, template[0], template[1], chirality, geometry, assignment)
    else:
        set_unresolved(m, evidence, unresolved or "The source identity is named but a complete structure was not encoded.", chirality, geometry)
    M[mid] = m
    return m


# Core monomer/product/catalysts.
m1 = add_molecule(f"{DOC}:mol:1", "1", "allenylboronic acid pinacol ester", smiles="C=C=CB1OC(C)(C)C(C)(C)O1")
m1["stereochemistry"]["axial"] = f(
    "terminal CH2-substituted allene drawing", {"axis": "allene", "stereogenic": False}, ev(2),
    "structure_derivation", "One terminal allene carbon is CH2, so the drawn allene is not axially stereogenic.")
m1["stereochemistry"]["representation_status"] = f(
    "RDKit reports an unspecified potential allene stereo element", "partial", ev(2), "structure_derivation",
    "Connectivity is encoded, but toolkit axial/allene stereo handling is limited; no false complete-encoding claim is made.")
m1["allene_features"] = [{
    "fragment_id": "allene-1",
    "carbon_atom_maps": f("source allene C=C=C", None, ev(2), status="pending_review"),
    "substitution": f("CH2=C=CH-BPin", "terminal CH2 and BPin-substituted terminus", ev(2)),
    "axial_stereochemistry": f("no axial descriptor shown", "not stereogenic for the drawn substitution", ev(2),
                                "structure_derivation", "One terminal allene carbon bears two hydrogens."),
    "local_environment": f("BPin directly bonded to terminal allene carbon", "allenyl-BPin", ev(2)),
}]

m2 = add_molecule(
    f"{DOC}:mol:2", "2", "poly(allenylboronic acid pinacol ester)",
    template=("2,3-vinyl-addition polyallene repeat unit bearing a pendant vinyl-BPin group; chain length/end groups unspecified",
              "Macromolecular chain length and end groups are unspecified, and the source tentatively assigns minor NMR resonances to Z insertion."),
    chirality="unspecified", geometry="mixture",
    assignment="Major/minor repeat-unit resonances retain a 1 vinyl-H:2 methylene-H ratio; minor signals are tentatively attributed by the authors to Z insertion, without a quantified E/Z ratio.")

add_molecule(f"{DOC}:mol:3", "3", "[(η3-C3H5)Ni(OC(O)CF3)]2", unresolved="Coordination complex 3 is drawn, but a complete machine SMILES was not encoded.")
add_molecule(f"{DOC}:mol:4", "4", "cationic Ni-allyl catalyst generated from 3 and AgSbF6", unresolved="Catalyst 4 is specified operationally as [3] + AgSbF6; complete coordination structure/speciation is unresolved.")
add_molecule(f"{DOC}:mol:5", "5", "trans-(C6F5)2Ni(SbPh3)2", unresolved="Coordination complex 5 is drawn, but a complete machine SMILES was not encoded.")

# Solvents/workup agents and comonomers/reagents.
add_molecule(f"{DOC}:mol:toluene", "toluene", "toluene", smiles="Cc1ccccc1")
add_molecule(f"{DOC}:mol:thf", "THF", "tetrahydrofuran", smiles="C1CCOC1")
add_molecule(f"{DOC}:mol:ether", "Et2O", "diethyl ether", smiles="CCOCC")
add_molecule(f"{DOC}:mol:cdcl3", "CDCl3", "deuterated chloroform", smiles="[2H]C(Cl)(Cl)Cl")
add_molecule(f"{DOC}:mol:dcm", "DCM", "dichloromethane", smiles="ClCCl")
add_molecule(f"{DOC}:mol:meoh", "MeOH", "methanol", smiles="CO")
add_molecule(f"{DOC}:mol:alumina", "basic alumina", "basic alumina", unresolved="Heterogeneous basic alumina is a material rather than a single molecular species.")

nonadiene = add_molecule(f"{DOC}:mol:nonadiene", "1,2-nonadiene", "1,2-nonadiene", page=4, smiles="C=C=CCCCCCC")
nonadiene["stereochemistry"]["axial"] = f(
    "terminal allene drawing", {"axis": "allene", "stereogenic": False}, ev(4), "structure_derivation",
    "One terminal allene carbon is CH2, so the allene is not axially stereogenic.")
nonadiene["stereochemistry"]["representation_status"] = f(
    "RDKit reports an unspecified potential allene stereo element", "partial", ev(4), "structure_derivation",
    "Connectivity is encoded, but toolkit axial/allene stereo handling is limited; no false complete-encoding claim is made.")

add_molecule(
    f"{DOC}:mol:beta_pinadiene", "β-pinadiene", "β-pinadiene", page=4,
    template=("Bicyclic terpene-derived terminal allene drawn in Figure 4A; bridge connectivity retained descriptively.",
              "The bicyclic structure is legible, but absolute/relative ring stereochemistry is not explicitly marked in the supplied figure and was not encoded."),
    chirality="unresolved", geometry="unresolved",
    assignment="No wedge/dash or named absolute configuration is visible in the supplied Figure 4 drawing; ring stereochemistry is therefore unresolved.")

add_molecule(
    f"{DOC}:mol:6", "6", "block copolymer of 1 and 1,2-nonadiene", page=4,
    template=("Block copolymer with a vinyl-BPin-bearing polyallene block and a poly(1,2-nonadiene) block; n and m unspecified.",
              "Polymer chain lengths, block lengths, end groups, and repeat-unit stereogeometry are not fully specified."),
    chirality="unspecified", geometry="unspecified")
add_molecule(
    f"{DOC}:mol:7", "7", "block copolymer of 1 and β-pinadiene", page=4,
    template=("Block copolymer with a vinyl-BPin-bearing polyallene block and a β-pinadiene-derived block; n and m unspecified.",
              "Polymer chain lengths, end groups, repeat-unit stereogeometry, and β-pinadiene ring stereochemistry are not fully specified."),
    chirality="unresolved", geometry="unspecified")

add_molecule(f"{DOC}:mol:phbr", "Ph-Br", "bromobenzene", page=5, smiles="Brc1ccccc1")
add_molecule(f"{DOC}:mol:pdpph3", "Pd(PPh3)4", "tetrakis(triphenylphosphine)palladium(0)", page=5,
             unresolved="The named Pd(PPh3)4 catalyst is explicit, but a complete coordination SMILES was not encoded.")
add_molecule(f"{DOC}:mol:k2co3", "K2CO3", "potassium carbonate", page=5, smiles="[K+].[K+].[O-]C([O-])=O")
add_molecule(f"{DOC}:mol:water", "H2O", "water", page=5, smiles="O")
add_molecule(
    f"{DOC}:mol:8", "8", "phenyl-substituted polyallene product 8", page=5,
    template=("Polyallene derived from 2 by replacing pendant vinyl-BPin groups with pendant vinyl-phenyl groups; m/end groups unspecified.",
              "Macromolecular chain length, end groups, and pendant alkene geometry are not fully specified."),
    chirality="unspecified", geometry="unspecified")


def charge(p, key, raw, normalized, evidence, method="original_extraction", derivation=None):
    p["charges"][key] = f(raw, normalized, evidence, method, derivation)


def participant(pid, mid, role, evidence, stage="reaction", product_role=None):
    p = st.new_participant(pid, mid, role, evidence, stage)
    if product_role:
        p["product_role"] = f(product_role, product_role, evidence)
    return p


def polymer_changes(evidence):
    changes = [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
         "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
         "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
    ]
    return changes, [{"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1"]}]


def set_classification(r, evidence, kind="polymerization"):
    if kind == "cross_coupling":
        changes = [
            {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
             "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "B"}],
             "before_order": 1, "after_order": 0, "change_type": "broken"},
        ]
        labels = [{"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1"]}]
        limitation = "Pendant C-B replacement is represented at motif level; no complete polymer atom mapping or chain-end inventory is asserted."
    else:
        changes, labels = polymer_changes(evidence)
        limitation = "Polymerization is represented at repeat-unit level; exact atom mapping, chain ends, and all repeat-unit stereogeometry are unresolved."
    r["classification"]["bond_changes"] = [f(None, c, evidence, "structure_derivation",
        "Bond-change motif inferred from the source reaction drawing and prose.") for c in changes]
    r["classification"]["labels"] = f(None, labels, evidence, "structure_derivation",
        "Taxonomy label derived from the motif-level bond changes; order decrease is not labeled as a new sigma bond.")
    r["classification"]["coverage"] = f(None, "partial", evidence, "inference", limitation)
    r["classification"]["limitations"] = [limitation]
    r["classification"]["model_votes"] = []


def reaction(suffix, evidence, record_type="specific_experiment", categories=None, name="Vinyl-addition polymerization", transformation="2,3-vinyl-addition polymerization"):
    rid = f"{DOC}:{suffix}"
    r = st.new_reaction(rid, DOC, evidence)
    r["record_type"] = record_type
    cats = categories or (["other"] if record_type == "specific_experiment" else None)
    if cats is not None:
        r["experiment_categories"] = f(cats, cats, evidence)
    r["location"] = f("PDF pages " + ", ".join(sorted({x.split(":page")[1].split(":")[0] for x in evidence})),
                      {"evidence_ids": evidence}, evidence)
    r["reaction_name"] = f(name, name, evidence)
    r["transformation_type"] = f(None, transformation, evidence, "inference",
                                  "Transformation class inferred from the source scheme and text.")
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
                   "unresolved": ["No explicitly associated SI was supplied; missing experimental details remain unreported."]}
    r["completeness"] = "partial"
    return r


def add_polymer_participants(r, catalyst_mid, solvent_mid=None, evidence=None, monomer_mid=None, product_mid=None):
    evidence = evidence or r["evidence_ids"]
    monomer_mid = monomer_mid or m1["molecule_id"]
    product_mid = product_mid or m2["molecule_id"]
    ps = participant("monomer", monomer_mid, "reactant", evidence)
    pc = participant("catalyst", catalyst_mid, "catalyst", evidence)
    pp = participant("polymer", product_mid, "product", evidence, product_role="reported polymer product")
    r["participants"] = [ps, pc, pp]
    if solvent_mid:
        r["participants"].append(participant("solvent", solvent_mid, "solvent", evidence))
    return ps, pc, pp


def set_solvent(r, raw, mid, evidence, concentration=None, inherited=False):
    method = "inherited_general_condition" if inherited else "original_extraction"
    r["solvent_conditions"]["composition"] = f(raw, [mid], evidence, method)
    if concentration is not None:
        r["solvent_conditions"]["concentration"] = f(f"[1]0 = {concentration} M", {"number": concentration, "unit": "M"}, evidence, method)
        r["solvent_conditions"]["concentration_basis"] = f("[1]0", m1["molecule_id"], evidence, method)


def add_poly_step(r, evidence, temperature=23, time_h=1.5, operation="polymerize allene monomer"):
    s = st.new_step("reaction", evidence)
    s["participant_ids"] = [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"]
    s["operation"] = f(operation, operation, evidence)
    if temperature is not None:
        s["temperature"] = f(f"{temperature} °C", {"number": temperature, "unit": "°C"}, evidence)
    if time_h is not None:
        s["time"] = f(f"{time_h} h", {"number": time_h, "unit": "h"}, evidence)
    r["steps"].append(s)


def outcome(r, oid, kind, targets, evidence, raw=None, normalized=None, basis_raw=None, basis_norm=None,
            yield_type=None, scope=None, status=None, alternatives=()):
    o = st.new_outcome(oid, kind, targets, evidence)
    o["value"] = f(raw, normalized, evidence, status=status, alternatives=alternatives)
    if basis_raw is not None or basis_norm is not None:
        o["measurement_basis"] = f(basis_raw, basis_norm, evidence)
    if yield_type:
        o["yield_type"] = f(yield_type, yield_type, evidence)
    if scope:
        o["scope"] = f(scope, scope, evidence)
    r["outcomes"].append(o)
    return o


R = []
generic_id = f"{DOC}:scheme:vinyl_addition"
g = reaction("scheme:vinyl_addition", ev(1) + ev(2), record_type="generic_scheme", categories=None)
add_polymer_participants(g, f"{DOC}:mol:3", evidence=g["evidence_ids"])
add_poly_step(g, g["evidence_ids"], 23, 1.5)
g["steps"][0]["atmosphere"] = f("N2", "N2", ev(2))
set_classification(g, g["evidence_ids"])
g["review"]["unresolved"] = ["Generic scheme only; it must not be treated as a single quantified experiment."]
R.append(g)


table = [
    (1, f"{DOC}:mol:3", f"{DOC}:mol:toluene", 92, "27,260 Da", "1.53"),
    (2, f"{DOC}:mol:4", f"{DOC}:mol:toluene", None, "70,600 Da", "1.20"),
    (3, f"{DOC}:mol:5", f"{DOC}:mol:toluene", 8, "N.D.", "N.D."),
    (4, f"{DOC}:mol:3", f"{DOC}:mol:thf", 82, "28,000 Da", "2.08"),
    (5, f"{DOC}:mol:3", f"{DOC}:mol:ether", 90, "35,180 Da", "1.49"),
    (6, f"{DOC}:mol:3", f"{DOC}:mol:cdcl3", 80, "34,420 Da", "1.52"),
    (7, f"{DOC}:mol:3", f"{DOC}:mol:dcm", 90, "29,450 Da", "1.38"),
    (8, f"{DOC}:mol:3", f"{DOC}:mol:dcm", 96, "39,970 Da", "1.51"),
]

for entry, cat, solv, yld, mn, disp in table:
    evidence = ev(2) + (ev(3) if 2 <= entry <= 7 else []) + ([f"{DOC}:page7:text"] if entry == 8 else [])
    r = reaction(f"table1:entry{entry}", evidence, categories=["optimization"] + (["scale_up"] if entry == 8 else []))
    r["general_scheme_ids"] = [generic_id]
    ps, pc, pp = add_polymer_participants(r, cat, solv, evidence)
    if entry == 8:
        charge(ps, "mass", "2.0 g", {"number": 2.0, "unit": "g"}, ev(2))
        charge(ps, "amount", "12.05 mmol", {"number": 12.05, "unit": "mmol"}, ev(2))
        r["reaction_scale"] = f("1 (2.0 g; 12.05 mmol)", {"reference_molecule_id": m1["molecule_id"], "mass": {"number": 2.0, "unit": "g"}, "amount": {"number": 12.05, "unit": "mmol"}}, ev(2))
        r["review"]["unresolved"].append("Table footnote e changes the monomer scale but does not explicitly report the scaled catalyst amount.")
    else:
        charge(ps, "amount", "0.6 mmol", {"number": 0.6, "unit": "mmol"}, ev(2))
        charge(pc, "amount", "0.0033 mmol Ni", {"number": 0.0033, "unit": "mmol"}, ev(2))
        charge(pc, "mol_percent", None, {"number": 0.55, "unit": "mol%"}, ev(2), "inference",
               "0.0033 mmol Ni / 0.6 mmol monomer × 100 = 0.55 mol% Ni.")
        pc["equivalence_reference"] = f(None, m1["molecule_id"], ev(2), "inference", "Loading denominator is monomer 1.")
        r["reaction_scale"] = f("1 (0.6 mmol)", {"reference_molecule_id": m1["molecule_id"], "amount": {"number": 0.6, "unit": "mmol"}}, ev(2))
    set_solvent(r, M[solv]["names"]["original"]["normalized_value"], solv, ev(2), 0.5)
    add_poly_step(r, ev(2), 23, 1.5)
    r["steps"][0]["atmosphere"] = f("N2", "N2", ev(2))
    if entry == 1:
        pa = participant("basic_alumina", f"{DOC}:mol:alumina", "workup_agent", ev(2), "workup")
        r["participants"].append(pa)
        s = st.new_step("workup_filter", ev(2), "workup")
        s["participant_ids"] = ["basic_alumina"]
        s["operation"] = f("filtered through basic alumina and concentrated to dryness", "filter through basic alumina; concentrate to dryness", ev(2))
        r["steps"].append(s)
    if entry == 2:
        a1 = f("Table 1 entry 2: 56%", {"number": 56, "unit": "%"}, ev(2))
        a2 = f("page 3 prose: 52%", {"number": 52, "unit": "%"}, ev(3))
        outcome(r, "yield", "yield", [m2["molecule_id"]], evidence, raw="56% (table) versus 52% (prose)",
                normalized=None, status="conflict", alternatives=[a1, a2], yield_type="reported polymer yield")
        r["review"]["unresolved"].append("Yield conflict: Table 1 gives 56%, whereas page 3 prose gives 52%; neither value was selected.")
    else:
        basis = "isolated polymer/material" if entry in (1, 3, 5, 8) else None
        outcome(r, "yield", "yield", [m2["molecule_id"]], evidence, raw=f"{yld}%", normalized={"number": yld, "unit": "%"},
                basis_raw=basis, basis_norm=basis, yield_type="reported polymer yield")
        if basis is None:
            r["review"]["unresolved"].append("The exact yield measurement/isolation basis is not explicitly stated for this table row in the supplied main text.")
    if entry == 1:
        outcome(r, "monomer_conversion", "conversion", [m1["molecule_id"]], ev(2), raw="99% consumption of 1",
                normalized={"number": 99, "unit": "%"}, basis_raw="1H NMR spectroscopy (toluene-d8)",
                basis_norm="1H NMR", scope="consumption of monomer 1 after 1.5 h at 23 °C")
    outcome(r, "target_mn", "other", [m2["molecule_id"]], ev(2), raw="30,000 Da", normalized="30,000 Da", scope="target number-average molecular weight")
    outcome(r, "measured_mn", "other", [m2["molecule_id"]], ev(2), raw=mn, normalized=mn,
            basis_raw="GPC with triple detection", basis_norm="GPC with triple detection", scope="measured number-average molecular weight")
    outcome(r, "dispersity", "other", [m2["molecule_id"]], ev(2), raw=disp, normalized=disp,
            basis_raw="Mw/Mn", basis_norm="calculated from Mw/Mn", scope="dispersity")
    if entry == 3:
        r["notes"].append(f("Isolated material was insoluble in THF solvent.", "Isolated material was insoluble in THF; GPC Mn and dispersity were not determined.", ev(2)))
    set_classification(r, evidence)
    R.append(r)


for suffix, conc, measured in [("high_mn_0p5M", 0.5, "130,800 Da"), ("high_mn_0p1M", 0.1, "90,350 Da")]:
    evidence = ev(3)
    r = reaction(suffix, evidence, categories=["optimization"])
    r["general_scheme_ids"] = [generic_id]
    add_polymer_participants(r, f"{DOC}:mol:3", f"{DOC}:mol:dcm", evidence)
    set_solvent(r, "DCM", f"{DOC}:mol:dcm", evidence, conc, inherited=True)
    add_poly_step(r, evidence, None, None, "polymerize 1 while targeting Mn = 100,000 Da")
    outcome(r, "target_mn", "other", [m2["molecule_id"]], evidence, "100,000 Da", "100,000 Da", scope="target number-average molecular weight")
    outcome(r, "measured_mn", "other", [m2["molecule_id"]], evidence, measured, measured,
            basis_raw="measured Mn; method referenced to Figure S5/S6", basis_norm="GPC (details in unavailable SI)", scope="measured number-average molecular weight")
    r["review"]["unresolved"].append("Exact charges, temperature, time, yield, and dispersity for this SI-referenced experiment are not reported in supplied pages 1–6.")
    set_classification(r, evidence)
    R.append(r)


# Living-character experiments.
r = reaction("mechanistic:mn_conversion", ev(3) + ev(4), categories=["mechanistic_study"])
r["general_scheme_ids"] = [generic_id]
add_polymer_participants(r, f"{DOC}:mol:3", evidence=r["evidence_ids"])
add_poly_step(r, r["evidence_ids"], None, None, "monitor Mn and dispersity as functions of monomer conversion")
outcome(r, "mn_conversion_relation", "other", [m2["molecule_id"]], r["evidence_ids"],
        "linear correlation between Mn and monomer conversion (Figure 3A, R2 = 0.9972)",
        {"relationship": "linear", "R2": 0.9972}, basis_raw="GPC Mn and monomer conversion plot", basis_norm="GPC/monomer-conversion monitoring")
outcome(r, "dispersity_trend", "other", [m2["molecule_id"]], ev(4),
        "dispersity narrowed with increasing monomer conversion", "decreasing dispersity trend", basis_raw="Figure 3A", basis_norm="plotted dispersity")
r["review"]["unresolved"].append("Individual plotted points and conversion measurement method are not tabulated in the supplied main text.")
set_classification(r, r["evidence_ids"])
R.append(r)

r = reaction("mechanistic:chain_extension", ev(4), categories=["mechanistic_study"])
r["general_scheme_ids"] = [generic_id]
add_polymer_participants(r, f"{DOC}:mol:3", evidence=ev(4))
for i, label in enumerate(["polymerize first portion of 1 until complete consumption", "add second portion of 1 and polymerize until complete consumption", "add third portion of 1 and polymerize until complete consumption"], 1):
    s = st.new_step(f"extension_{i}", ev(4))
    s["participant_ids"] = ["monomer", "catalyst", "polymer"]
    s["operation"] = f(label, label, ev(4))
    r["steps"].append(s)
for i, (mn, d) in enumerate([("40,160 Da", "1.40"), ("67,550 Da", "1.39"), ("93,910 Da", "1.42")], 1):
    outcome(r, f"mn_after_extension_{i}", "other", [m2["molecule_id"]], ev(4), mn, mn, basis_raw="GPC trace in Figure 3B", basis_norm="GPC", scope=f"Mn after stage {i}")
    outcome(r, f"dispersity_after_extension_{i}", "other", [m2["molecule_id"]], ev(4), d, d, basis_raw="GPC trace in Figure 3B", basis_norm="GPC", scope=f"dispersity after stage {i}")
r["review"]["unresolved"].append("Portion sizes, stage times, temperatures, solvent, and absolute conversion method are in unavailable SI or unreported here.")
set_classification(r, ev(4))
R.append(r)

r = reaction("sample:low_mn_end_group", ev(4), categories=["other"])
r["general_scheme_ids"] = [generic_id]
add_polymer_participants(r, f"{DOC}:mol:3", evidence=ev(4))
add_poly_step(r, ev(4), None, None, "prepare low-Mn sample of 2 for end-group analysis")
outcome(r, "gpc_mn", "other", [m2["molecule_id"]], ev(4), "4865 Da", "4865 Da", basis_raw="GPC", basis_norm="GPC", scope="number-average molecular weight")
outcome(r, "dispersity", "other", [m2["molecule_id"]], ev(4), "1.45", "1.45", basis_raw="GPC", basis_norm="GPC", scope="dispersity")
outcome(r, "nmr_mn", "other", [m2["molecule_id"]], ev(4), "3363 Da", "3363 Da", basis_raw="1H NMR allyl end-group resonances", basis_norm="1H NMR end-group analysis", scope="number-average molecular weight")
r["review"]["unresolved"].append("Preparation conditions and isolated yield for this low-Mn sample are not reported in supplied pages 1–6.")
set_classification(r, ev(4))
R.append(r)


def block_record(suffix, comonomer_mid, product_mid, yield_pct, amount_mg, initial_mn, initial_d, final_mn, final_d, theo, actual):
    evidence = ev(4) + ev(5)
    r = reaction(suffix, evidence, categories=["other"], name="Block copolymerization", transformation="sequential vinyl-addition block copolymerization")
    r["general_scheme_ids"] = [generic_id]
    p1 = participant("monomer_1", m1["molecule_id"], "reactant", evidence)
    p2 = participant("comonomer", comonomer_mid, "reactant", evidence)
    pc = participant("catalyst", f"{DOC}:mol:3", "catalyst", evidence)
    pp = participant("block_copolymer", product_mid, "product", evidence, product_role="isolated block copolymer")
    r["participants"] = [p1, p2, pc, pp]
    charge(p1, "mass", "300 mg", {"number": 300, "unit": "mg"}, ev(5))
    charge(p2, "mass", "300 mg", {"number": 300, "unit": "mg"}, ev(5))
    r["reaction_scale"] = f("1 (300 mg) and comonomer (300 mg)", {"reference_molecule_id": m1["molecule_id"], "mass_each": {"number": 300, "unit": "mg"}}, ev(5))
    s1 = st.new_step("first_block", evidence)
    s1["participant_ids"] = ["monomer_1", "catalyst", "block_copolymer"]
    s1["operation"] = f("1) [3], 23 °C, 2 h", "polymerize monomer 1 with catalyst 3", ev(4))
    s1["temperature"] = f("23 °C", {"number": 23, "unit": "°C"}, ev(4))
    s1["time"] = f("2 h", {"number": 2, "unit": "h"}, ev(4))
    s2 = st.new_step("second_block", evidence)
    s2["participant_ids"] = ["comonomer", "catalyst", "block_copolymer"]
    s2["operation"] = f("Upon complete consumption of 1, comonomer was added; final polymer precipitated from MeOH after 2 h.", "add comonomer after complete consumption of 1; polymerize second block", evidence)
    s2["temperature"] = f("23 °C (same Figure 4 sequence)", {"number": 23, "unit": "°C"}, ev(4), "inherited_general_condition")
    s2["time"] = f("2 h", {"number": 2, "unit": "h"}, ev(4))
    r["steps"] = [s1, s2]
    pm = participant("methanol", f"{DOC}:mol:meoh", "workup_agent", evidence, "workup")
    r["participants"].append(pm)
    sw = st.new_step("precipitation", evidence, "workup")
    sw["participant_ids"] = ["block_copolymer", "methanol"]
    sw["operation"] = f("precipitated from MeOH", "precipitate polymer from methanol", evidence)
    r["steps"].append(sw)
    outcome(r, "isolated_yield", "yield", [product_mid], evidence, f"{yield_pct}%", {"number": yield_pct, "unit": "%"}, basis_raw="isolated block copolymer", basis_norm="isolated", yield_type="isolated yield")
    outcome(r, "isolated_amount", "product_amount", [product_mid], ev(5), f"{amount_mg} mg", {"number": amount_mg, "unit": "mg"}, basis_raw="Table 2 isolated yield (mg)", basis_norm="isolated mass")
    for oid, raw, scope in [("initial_mn", initial_mn, "first-block Mn"), ("initial_dispersity", initial_d, "first-block dispersity"), ("final_mn", final_mn, "final block-copolymer Mn"), ("final_dispersity", final_d, "final block-copolymer dispersity")]:
        outcome(r, oid, "other", [product_mid], ev(5), raw, raw, basis_raw="GPC with triple detection" if "mn" in oid else "calculated from Mw/Mn", basis_norm="GPC with triple detection" if "mn" in oid else "Mw/Mn", scope=scope)
    outcome(r, "theoretical_incorporation_1", "other", [product_mid], ev(5), f"{theo}%", {"number": theo, "unit": "%"}, basis_raw="calculated from feed ratio", basis_norm="feed-ratio calculation", scope="theoretical monomer 1 incorporation")
    outcome(r, "actual_incorporation_1", "other", [product_mid], ev(5), f"{actual}%", {"number": actual, "unit": "%"}, basis_raw="1H NMR spectroscopy", basis_norm="1H NMR", scope="actual monomer 1 incorporation")
    set_classification(r, evidence)
    return r


R.append(block_record("block:6", f"{DOC}:mol:nonadiene", f"{DOC}:mol:6", 46, 275.2, "37,360 Da", "1.30", "80,510 Da", "1.43", 43, 31))
R.append(block_record("block:7", f"{DOC}:mol:beta_pinadiene", f"{DOC}:mol:7", 71, 424.5, "35,330 Da", "1.34", "66,460 Da", "1.77", 49, 37))


def recovery_probe(suffix, comonomer_mid, product_mid):
    evidence = ev(5)
    r = reaction(suffix, evidence, categories=["mechanistic_study"], name="Block-copolymer precipitation/supernatant probe", transformation="sequential vinyl-addition polymerization followed by fractionation")
    r["general_scheme_ids"] = [generic_id]
    r["participants"] = [
        participant("monomer_1", m1["molecule_id"], "reactant", evidence),
        participant("comonomer", comonomer_mid, "reactant", evidence),
        participant("catalyst", f"{DOC}:mol:3", "catalyst", evidence),
        participant("block_copolymer", product_mid, "product", evidence, product_role="precipitated block copolymer fraction"),
        participant("homopolymer_2", m2["molecule_id"], "product", evidence, product_role="polymer recovered from methanol supernatant"),
        participant("methanol", f"{DOC}:mol:meoh", "workup_agent", evidence, "workup"),
    ]
    s = st.new_step("repeat_block_synthesis", evidence)
    s["participant_ids"] = ["monomer_1", "comonomer", "catalyst", "block_copolymer", "homopolymer_2"]
    s["operation"] = f("conducted a synthesis as described above", "repeat block-copolymer synthesis", evidence)
    w = st.new_step("fractionation", evidence, "workup")
    w["participant_ids"] = ["block_copolymer", "homopolymer_2", "methanol"]
    w["operation"] = f("precipitated solid polymer from MeOH, collected and concentrated supernatant to dryness", "methanol precipitation; collect supernatant; concentrate to dryness", evidence)
    r["steps"] = [s, w]
    outcome(r, "homopolymer_identification", "other", [m2["molecule_id"]], evidence, "supernatant solid was homopolymer 2", "homopolymer 2 detected", basis_raw="1H NMR spectral analysis", basis_norm="1H NMR")
    outcome(r, "total_polymer_mass_recovery", "other", [product_mid, m2["molecule_id"]], evidence, ">74% (statement applies to both probe cases)", {"number": 74, "unit": "%", "qualifier": "gt", "raw_text": ">74%"}, basis_raw="recovered mass of all polymeric products", basis_norm="total recovered polymer mass", scope="combined precipitate plus supernatant polymers")
    r["review"]["unresolved"].append("The >74% statement is shared across both probe cases; case-specific values and full conditions are not given in the supplied main text.")
    set_classification(r, evidence)
    return r


R.append(recovery_probe("mechanistic:recovery_probe_6", f"{DOC}:mol:nonadiene", f"{DOC}:mol:6"))
R.append(recovery_probe("mechanistic:recovery_probe_7", f"{DOC}:mol:beta_pinadiene", f"{DOC}:mol:7"))


def reverse_order(suffix, first_mid, block_mid):
    evidence = ev(5)
    r = reaction(suffix, evidence, categories=["mechanistic_study"], name="Reverse-order block copolymerization", transformation="sequential vinyl-addition polymerization")
    r["general_scheme_ids"] = [generic_id]
    r["participants"] = [
        participant("first_monomer", first_mid, "reactant", evidence),
        participant("monomer_1", m1["molecule_id"], "reactant", evidence),
        participant("catalyst", f"{DOC}:mol:3", "catalyst", evidence),
        participant("nominal_block_product", block_mid, "product", evidence, product_role="nominal block-copolymer product/fraction"),
        participant("homopolymer_2", m2["molecule_id"], "product", evidence, product_role="observed homopolymer side fraction"),
    ]
    s1 = st.new_step("first_block", evidence)
    s1["participant_ids"] = ["first_monomer", "catalyst", "nominal_block_product"]
    s1["operation"] = f("comonomer was polymerized before adding 1", "polymerize comonomer first", evidence)
    s2 = st.new_step("add_1", evidence)
    s2["participant_ids"] = ["monomer_1", "catalyst", "nominal_block_product", "homopolymer_2"]
    s2["operation"] = f("then adding 1 to the reaction", "add monomer 1 after first block", evidence)
    r["steps"] = [s1, s2]
    outcome(r, "homopolymer_2_formation", "other", [m2["molecule_id"]], evidence, "resulted in formation of 2", "homopolymer 2 formed", basis_raw="authors' main-text observation; Table S4/Figures S30-S31 cited", basis_norm="source observation; quantitative details unavailable")
    r["review"]["unresolved"].append("Exact conditions, product distributions, and recovered masses are only referenced to unavailable SI Table S4/Figures S30-S31.")
    set_classification(r, evidence)
    return r


R.append(reverse_order("mechanistic:reverse_nonadiene_first", f"{DOC}:mol:nonadiene", f"{DOC}:mol:6"))
R.append(reverse_order("mechanistic:reverse_beta_pinadiene_first", f"{DOC}:mol:beta_pinadiene", f"{DOC}:mol:7"))


# Suzuki--Miyaura post-polymerization modification.
r = reaction("derivatization:8", ev(5), categories=["derivatization"], name="Pd-catalyzed Suzuki–Miyaura cross-coupling", transformation="post-polymerization Suzuki–Miyaura cross-coupling")
r["participants"] = [
    participant("polymer_2", m2["molecule_id"], "reactant", ev(5)),
    participant("phenyl_bromide", f"{DOC}:mol:phbr", "reactant", ev(5)),
    participant("pd_catalyst", f"{DOC}:mol:pdpph3", "catalyst", ev(5)),
    participant("base", f"{DOC}:mol:k2co3", "reagent", ev(5)),
    participant("thf", f"{DOC}:mol:thf", "solvent", ev(5)),
    participant("water", f"{DOC}:mol:water", "solvent", ev(5)),
    participant("product_8", f"{DOC}:mol:8", "product", ev(5), product_role="cross-coupled polymer product"),
]
charge(r["participants"][1], "equivalents", "3 equiv.", {"number": 3, "unit": "equiv"}, ev(5))
r["participants"][1]["equivalence_reference"] = f("relative to polymer 2 repeat units", m2["molecule_id"], ev(5))
charge(r["participants"][2], "mol_percent", "3 mol %", {"number": 3, "unit": "mol%"}, ev(5))
r["participants"][2]["equivalence_reference"] = f("relative to polymer 2 reactive repeat units", m2["molecule_id"], ev(5))
charge(r["participants"][3], "equivalents", "6 equiv.", {"number": 6, "unit": "equiv"}, ev(5))
r["participants"][3]["equivalence_reference"] = f("relative to polymer 2 reactive repeat units", m2["molecule_id"], ev(5))
r["solvent_conditions"]["composition"] = f("THF/H2O", [f"{DOC}:mol:thf", f"{DOC}:mol:water"], ev(5))
s = st.new_step("cross_coupling", ev(5))
s["participant_ids"] = [p["participant_id"] for p in r["participants"]]
s["operation"] = f("Pd(PPh3)4, K2CO3, THF/H2O, N2, 80 °C, 20 h", "Suzuki–Miyaura cross-coupling", ev(5))
s["temperature"] = f("80 °C", {"number": 80, "unit": "°C"}, ev(5))
s["time"] = f("20 h", {"number": 20, "unit": "h"}, ev(5))
s["atmosphere"] = f("N2", "N2", ev(5))
r["steps"] = [s]
outcome(r, "conversion", "conversion", [m2["molecule_id"]], ev(5), "96%", {"number": 96, "unit": "%"}, scope="conversion of pendant vinyl-BPin motifs to coupled product")
outcome(r, "yield", "yield", [f"{DOC}:mol:8"], ev(5), "93%", {"number": 93, "unit": "%"}, yield_type="reported product yield")
outcome(r, "product_mn", "other", [f"{DOC}:mol:8"], ev(5), "52,000 Da", "52,000 Da", basis_raw="GPC Figure 5", basis_norm="GPC", scope="product 8 Mn")
outcome(r, "product_dispersity", "other", [f"{DOC}:mol:8"], ev(5), "1.86", "1.86", basis_raw="GPC Figure 5", basis_norm="GPC", scope="product 8 dispersity")
outcome(r, "starting_polymer_mn", "other", [m2["molecule_id"]], ev(5), "61,410 Da", "61,410 Da", basis_raw="GPC Figure 5", basis_norm="GPC", scope="starting polymer 2 Mn")
outcome(r, "starting_polymer_dispersity", "other", [m2["molecule_id"]], ev(5), "1.87", "1.87", basis_raw="GPC Figure 5", basis_norm="GPC", scope="starting polymer 2 dispersity")
r["review"]["unresolved"].append("The main text does not explicitly state measurement methods for the 96% conversion or 93% yield; solvent ratio/volume and polymer charge are also unreported.")
set_classification(r, ev(5), "cross_coupling")
R.append(r)


bundle = json.loads((RUN / "submitted.json").read_text())
bundle["molecules"] = list(M.values())
bundle["reactions"] = R
bundle["source_context"] = [
    {"context_id": f"{DOC}:page1", "page": 1, "kind": "graphical_abstract_and_background", "evidence_ids": ev(1), "finding": "Graphical abstract summarizes the current polymerization/copolymerization/cross-coupling; introduction and Figure 1A-C prior-work panels are not new experiments."},
    {"context_id": f"{DOC}:page2", "page": 2, "kind": "scheme_table_prose", "evidence_ids": ev(2), "finding": "Current-work generic scheme and Table 1 entries 1-8; catalyst drawings and footnotes inspected."},
    {"context_id": f"{DOC}:page3", "page": 3, "kind": "characterization_and_discussion", "evidence_ids": ev(3), "finding": "Figure 2 is NMR characterization, not a separate synthesis. It supports 2,3-insertion and a tentative minor Z-insertion assignment. Prose/table entry-2 yield conflict retained."},
    {"context_id": f"{DOC}:page4", "page": 4, "kind": "kinetics_chain_extension_block_schemes", "evidence_ids": ev(4), "finding": "Figure 3 contains kinetic/chain-extension experiments; Figure 4 contains two block-copolymer syntheses. GPC plots are characterization outcomes, not standalone syntheses."},
    {"context_id": f"{DOC}:page5", "page": 5, "kind": "block_table_mechanistic_probes_and_derivatization", "evidence_ids": ev(5), "finding": "Table 2 duplicates Figure 4 block experiments and adds quantities/properties. Main-text recovery and reverse-order probes are separate mechanistic experiments. Figure 5 is the Suzuki–Miyaura derivatization."},
    {"context_id": f"{DOC}:page6", "page": 6, "kind": "associated_content_and_references", "evidence_ids": ev(6), "finding": "No new current-paper experiment. SI is stated to exist online but was not supplied/associated; references are prior work."},
]
st.write(RUN / "submitted.json", bundle)

raw = {
    "document_id": DOC,
    "claimed_page_scope": [1, 2, 3, 4, 5, 6],
    "semantic_accuracy_verified": False,
    "molecule_ids": list(M),
    "reaction_ids": [r["reaction_id"] for r in R],
    "notes": [
        "All six supplied page images were opened at original detail; local crops were used for Figures/Table details.",
        "Page 7 was not claimed as visually covered; its extracted footnote 45 text was used only as necessary context confirming entry 8 isolated yield/scale-up interpretation.",
        "No SI file was associated. SI-referenced details remain unreported.",
        "Machine-derived structures and classifications remain needs_review.",
        "Table 1 entry 2 yield conflict (56% table vs 52% prose) is retained without selecting a value.",
    ],
}
st.write(RUN / "raw_model/extracted_pages_0001_0006.json", raw)
st.write(RUN / "raw_model/source_mentions.json", st.source_mentions(
    json.loads((RUN / "dataset.json").read_text()),
    ["polymer 2", "block copolymer 6", "block copolymer 7", "phenyl bromide", "1,2-nonadiene", "β-pinadiene"],
    context=180,
))
diagnostics = []
for mol in M.values():
    smiles = mol["smiles"]["normalized_value"]
    if smiles:
        diagnostics.append({"molecule_id": mol["molecule_id"], "input_smiles": smiles,
                            "diagnostic": st.stereo_inspect(smiles)})
(RUN / "raw_model/stereo_diagnostics.jsonl").write_text(
    "".join(json.dumps(row, ensure_ascii=False) + "\n" for row in diagnostics)
)
print(json.dumps({"molecules": len(M), "reactions": len(R)}, ensure_ascii=False))
