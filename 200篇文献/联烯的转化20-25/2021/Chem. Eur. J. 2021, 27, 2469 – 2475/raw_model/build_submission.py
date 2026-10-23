#!/usr/bin/env python3
"""Build the source-derived v3.1 machine bundle for PDF pages 1--6.

This is a transparent extraction artifact, not a source-review certificate.  It
uses only the prepared dataset and the six visually inspected page images.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_f15761c9ba9d1af6"
RUN = ROOT / "runs/current"


def ev(page: int):
    return [f"{DOC}:page{page}:text", f"{DOC}:page{page}:image"]


def f(raw=None, normalized=None, evidence=(), method="original_extraction", status=None, derivation=None, alternatives=()):
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   derivation=derivation, alternatives=alternatives)


def missing(status="unreported"):
    method = "not_applicable" if status == "not_applicable" else "original_extraction"
    return f(method=method, status=status)


def q(number, unit, raw_text=None, qualifier=None):
    out = {"number": number, "unit": unit}
    if raw_text is not None:
        out["raw_text"] = raw_text
    if qualifier is not None:
        out["qualifier"] = qualifier
    return out


def set_fact(container, key, raw, normalized, evidence, method="original_extraction", derivation=None, status=None):
    container[key] = f(raw, normalized, evidence, method, status=status, derivation=derivation)


stereo_checks = json.loads((RUN / "raw_stereo_checks.json").read_text())


def molecule(mid, label, name, evidence, smiles=None, stereo_kind="unresolved", identity_note=None):
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    m["names"]["original"] = f(name, name, evidence)
    m["names"]["english"] = f(name, name, evidence)
    m["names"]["abbreviations"] = missing()
    m["names"]["aliases"] = missing()
    m["cas_numbers"] = missing()
    m["external_properties"] = []
    m["identity_note"] = (f(None, identity_note, evidence, "inference", derivation=identity_note)
                          if identity_note else missing())
    if smiles is None:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing()
        m["generic_template"] = None
        reason = "The page provides a label or drawing, but a complete stereochemically reliable SMILES was not established from pages 1-6."
        for key in ("atom_centers", "double_bonds", "axial", "mixture", "assignment"):
            m["stereochemistry"][key] = missing()
        for key in ("chirality_status", "geometry_status", "representation_status"):
            m["stereochemistry"][key] = f(None, "unresolved", evidence, "inference", derivation=reason)
        m["stereochemistry"]["isomeric_smiles"] = missing()
        m["stereochemistry"]["identity_key"] = missing()
        return m

    m["representation_type"] = "smiles"
    derivation = "Candidate connectivity transcribed from the page drawing/name; RDKit parsing and stereo diagnostics do not establish source correctness."
    m["smiles"] = f(None, smiles, evidence, "structure_derivation", derivation=derivation)
    m["generic_template"] = None
    diag = None
    for row in stereo_checks.values():
        if row["smiles"] == smiles:
            diag = row["result"]
            break
    if diag is None:
        diag = st.stereo_inspect(smiles)
    stereo = m["stereochemistry"]
    if stereo_kind == "achiral":
        stereo["atom_centers"] = f(None, [], evidence, "structure_derivation", derivation="No tetrahedral stereocenter is present in the candidate graph.")
        stereo["double_bonds"] = f(None, [], evidence, "structure_derivation", derivation="No E/Z-configurable double bond is present in the candidate graph.")
        stereo["axial"] = f(None, [], evidence, "structure_derivation", derivation="No stereogenic axis is present in the candidate graph.")
        stereo["assignment"] = f(None, "achiral candidate graph", evidence, "structure_derivation", derivation="Symmetry/connectivity of the candidate graph was inspected; this does not validate the source identity.")
        stereo["chirality_status"] = f(None, "achiral", evidence, "structure_derivation", derivation="Candidate graph has no stereogenic center or axis.")
        stereo["geometry_status"] = f(None, "not_applicable", evidence, "structure_derivation", derivation="Candidate graph has no E/Z-configurable double bond.")
        stereo["representation_status"] = f(None, "encoded", evidence, "structure_derivation", derivation="No applicable stereochemical feature remains unencoded in the candidate graph.")
    elif stereo_kind == "allene_achiral_partial":
        stereo["atom_centers"] = f(None, [], evidence, "structure_derivation", derivation="Candidate graph contains no tetrahedral stereocenter.")
        stereo["double_bonds"] = f("terminal allene C=C=C", [{"feature": "terminal allene", "configuration": "not assigned"}], evidence,
                                    "structure_derivation", derivation="Allene drawn with a terminal CH2 end; ordinary E/Z assignment is not applicable.")
        stereo["axial"] = f("terminal allene; no axial label", [], evidence, "structure_derivation",
                              derivation="One allene terminus is CH2, so the drawn allene is not stereogenic about its axis.")
        stereo["assignment"] = f(None, "achiral terminal allene candidate", evidence, "structure_derivation",
                                   derivation="Terminal CH2 and duplicate phenyl substituents preclude an assigned stereogenic axis/center in the candidate.")
        stereo["chirality_status"] = f(None, "achiral", evidence, "structure_derivation", derivation="Candidate terminal allene is not axially chiral.")
        stereo["geometry_status"] = f(None, "not_applicable", evidence, "structure_derivation", derivation="Ordinary E/Z notation is not applicable to the terminal allene.")
        stereo["representation_status"] = f(None, "partial", evidence, "structure_derivation", derivation="Toolkit axial/allene support is limited; the source allene description is retained separately.")
    else:
        centers = diag.get("atom_centers", [])
        stereo["atom_centers"] = f("drawing has no wedge/dash or R/S label", centers, evidence, "structure_derivation",
                                     derivation="Potential tetrahedral centers were detected from the candidate graph; configuration is not specified in the page drawing.")
        stereo["double_bonds"] = f("no E/Z label shown", [], evidence, "structure_derivation",
                                    derivation="No source E/Z assignment was found on the inspected drawing.")
        stereo["axial"] = missing()
        stereo["assignment"] = f("configuration not specified", "unspecified", evidence)
        stereo["chirality_status"] = f("configuration not specified", "unspecified", evidence)
        stereo["geometry_status"] = f("no E/Z assignment shown", "not_applicable", evidence,
                                      "structure_derivation", derivation="Only carbonyl/terminal alkene bonds without an E/Z assignment are present in the candidate graph.")
        stereo["representation_status"] = f(None, "partial", evidence, "structure_derivation",
                                              derivation="Potential tetrahedral configuration is not encoded because it is not specified by the source.")
    stereo["mixture"] = missing()
    stereo["isomeric_smiles"] = f(None, diag["canonical_isomeric_smiles"], evidence, "structure_derivation",
                                    derivation="RDKit canonical isomeric SMILES of the candidate graph; unassigned features remain unassigned.")
    stereo["identity_key"] = f(None, diag["stereo_graph_key"], evidence, "structure_derivation",
                                derivation="SHA-256 stereo graph key from the candidate SMILES; not a sample identity.")
    return m


def participant(pid, mid, role, evidence, stage="reaction", mass=None, amount=None, volume=None,
                equivalents=None, mol_percent=None, reference=None, product_role=None):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    if product_role is not None:
        p["product_role"] = f(product_role, product_role, evidence)
    for key, value in (("mass", mass), ("amount", amount), ("volume", volume),
                       ("equivalents", equivalents), ("mol_percent", mol_percent)):
        if value is not None:
            raw = value.get("raw_text") or f"{value.get('number')} {value.get('unit')}"
            p["charges"][key] = f(raw, value, evidence)
    if reference is not None:
        p["equivalence_reference"] = f(reference, reference, evidence)
    return p


def step(step_id, stage, evidence, participants=(), operation=None, temperature=None, time=None, atmosphere=None, pressure=None,
         method="original_extraction"):
    s = st.new_step(step_id, evidence, stage=stage)
    s["participant_ids"] = list(participants)
    if operation is not None:
        s["operation"] = f(operation, operation, evidence, method)
    if temperature is not None:
        raw = temperature if isinstance(temperature, str) else temperature.get("raw_text") or f"{temperature['number']} {temperature['unit']}"
        s["temperature"] = f(raw, temperature, evidence, method)
    if time is not None:
        raw = time if isinstance(time, str) else time.get("raw_text") or f"{time['number']} {time['unit']}"
        s["time"] = f(raw, time, evidence, method)
    if atmosphere is not None:
        s["atmosphere"] = f(atmosphere, atmosphere, evidence, method)
    if pressure is not None:
        s["pressure"] = f(pressure, pressure, evidence, method)
    return s


def outcome(oid, kind, targets, evidence, raw, normalized, basis=None, role=None, scope=None, yield_type=None):
    o = st.new_outcome(oid, kind, targets, evidence)
    o["value"] = f(raw, normalized, evidence)
    if basis is not None:
        o["measurement_basis"] = f(basis, basis, evidence)
    if role is not None:
        o["product_role"] = f(role, role, evidence)
    if scope is not None:
        o["scope"] = f(scope, scope, evidence)
    if yield_type is not None:
        o["yield_type"] = f(yield_type, yield_type, evidence)
    return o


def classification(kind, evidence):
    specs = {
        "complexation": [
            {"change_id": "b1", "atoms": [{"element": "P"}, {"element": "Au"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        ],
        "hydration": [
            {"change_id": "b1", "atoms": [{"element": "C"}, {"element": "O"}], "before_order": 0, "after_order": 2, "change_type": "formed"},
            {"change_id": "b2", "atoms": [{"element": "C"}, {"element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 2},
            {"change_id": "b3", "atoms": [{"element": "C"}, {"element": "C"}], "before_order": 3, "after_order": 1, "change_type": "order_decreased"},
        ],
        "cycloisomerization": [
            {"change_id": "b1", "atoms": [{"element": "C"}, {"element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        ],
        "alkoxycyclization": [
            {"change_id": "b1", "atoms": [{"element": "C"}, {"element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b2", "atoms": [{"element": "C"}, {"element": "O"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b3", "atoms": [{"element": "C"}, {"element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain"},
        ],
        "hydroalkoxylation": [
            {"change_id": "b1", "atoms": [{"element": "C"}, {"element": "O"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b2", "atoms": [{"element": "C"}, {"element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain"},
        ],
    }
    if kind == "none":
        return {"standard_id": st.STANDARD, "mapped_reaction_smiles": missing(), "bond_changes": [],
                "labels": f(None, [], evidence, "inference", derivation="No observable reaction; no bond-formation label assigned."),
                "coverage": f(None, "not_applicable", evidence, "inference", derivation="Negative control/no observable transformation."),
                "limitations": ["No expected-product bonds were imputed for the negative control."], "model_votes": []}
    if kind not in specs:
        return {"standard_id": st.STANDARD, "mapped_reaction_smiles": missing(), "bond_changes": [],
                "labels": missing(), "coverage": f(None, "unresolved", evidence, "inference", derivation="Insufficient source structure/mapping for defensible classification."),
                "limitations": ["No atom mapping was asserted."], "model_votes": []}
    changes = specs[kind]
    labels = {}
    for ch in changes:
        if ch["change_type"] in {"formed", "hydrogen_gain", "order_increased"}:
            lab = st.bond_label(ch["atoms"][0]["element"], ch["atoms"][1]["element"], ch["after_order"])
            labels.setdefault(lab["code"], dict(lab, supporting_change_ids=[]))["supporting_change_ids"].append(ch["change_id"])
    deriv = "Partial bond-change interpretation from the page reaction drawing; atom mapping and complete atom balance were not established."
    return {"standard_id": st.STANDARD, "mapped_reaction_smiles": missing(),
            "bond_changes": [f(None, ch, evidence, "structure_derivation", derivation=deriv) for ch in changes],
            "labels": f(None, list(labels.values()), evidence, "structure_derivation", derivation=deriv),
            "coverage": f(None, "partial", evidence, "inference", derivation=deriv),
            "limitations": ["No atom-mapped reaction SMILES was asserted.", "Hydrogen changes are net changes; donor/acceptor identity and full atom balance remain unresolved.", "Only confidently visible net formed bonds are classified."],
            "model_votes": []}


def reaction(rid, evidence, record_type="specific_experiment", categories=None, location=None, name=None,
             transformation=None, general=(), class_kind="unresolved"):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories or ["other"], categories or ["other"], evidence)
    r["location"] = f(location, location, evidence)
    r["reaction_name"] = f(name, name, evidence) if name else missing()
    r["transformation_type"] = (f(None, transformation, evidence, "inference", derivation="Normalized from the source reaction title/drawing.")
                                    if transformation else missing())
    r["general_scheme_ids"] = [f"{DOC}:{x}" for x in general]
    r["classification"] = classification(class_kind, evidence)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
                   "unresolved": ["Machine extraction requires expert review; schema validation is not chemical verification."]}
    r["completeness"] = "partial"
    return r


def solvent_fields(r, evidence, composition=None, ratio=None, total=None, basis=None):
    if composition is not None:
        r["solvent_conditions"]["composition"] = f(composition, composition, evidence)
    if ratio is not None:
        r["solvent_conditions"]["mixture_ratio"] = f(ratio, ratio, evidence)
    if total is not None:
        r["solvent_conditions"]["total_volume"] = f(total.get("raw_text"), total, evidence)
    if basis is not None:
        r["solvent_conditions"]["concentration_basis"] = f(basis, basis, evidence)


def procedure_fields(r, evidence, original=None, addition=None, pretreatment=None, general_ref=None):
    if original:
        r["procedure"]["original_text"] = f(original, original, evidence)
    if addition:
        r["procedure"]["addition_order"] = f(addition, addition, evidence)
    if pretreatment:
        r["procedure"]["pretreatment"] = f(pretreatment, pretreatment, evidence)
    if general_ref:
        r["procedure"]["general_procedure_reference"] = f(general_ref, general_ref, evidence)


dataset = json.loads((RUN / "dataset.json").read_text())
bundle = st.empty(dataset["documents"])

# Molecule registry.  Organogold and visually ambiguous products deliberately
# remain unresolved structures; simple and clearly drawn organic species have
# candidate SMILES with explicit stereo status.
molecules = [
    molecule("m_1a", "1a", "phosphanoxy-substituted phosphaalkene 1a (R = nBu)", ev(2) + ev(5)),
    molecule("m_1b", "1b", "phosphanoxy-substituted phosphaalkene 1b (R = Ph)", ev(2) + ev(5)),
    molecule("m_2a", "2a", "mono(chlorogold) complex 2a (R = nBu)", ev(2) + ev(5)),
    molecule("m_2b", "2b", "mono(chlorogold) complex 2b (R = Ph)", ev(2) + ev(5)),
    molecule("m_3a", "3a", "bis(chlorogold) complex 3a (R = nBu)", ev(2) + ev(3) + ev(6)),
    molecule("m_3b", "3b", "bis(chlorogold) complex 3b (R = Ph)", ev(2) + ev(3) + ev(6)),
    molecule("m_4", "4", "ethoxydiphenylphosphine-chlorogold(I) complex 4", ev(4)),
    molecule("m_tht_aucl", "(tht)AuCl", "tetrahydrothiophene gold(I) chloride", ev(2) + ev(5) + ev(6)),
    molecule("m_dcm", "CH2Cl2", "dichloromethane", ev(2) + ev(6), "ClCCl", "achiral"),
    molecule("m_nhexane", "n-hexane", "n-hexane", ev(5) + ev(6), "CCCCCC", "achiral"),
    molecule("m_acetone", "acetone", "acetone", ev(3) + ev(6), "CC(C)=O", "achiral"),
    molecule("m_water", "H2O", "water", ev(3) + ev(6), "O", "achiral"),
    molecule("m_methanol", "MeOH", "methanol", ev(4) + ev(6), "CO", "achiral"),
    molecule("m_tfe", "CF3CH2OH", "2,2,2-trifluoroethan-1-ol", ev(5), "OCC(F)(F)F", "achiral"),
    molecule("m_toluene", "toluene", "toluene", ev(4) + ev(6), "Cc1ccccc1", "achiral"),
    molecule("m_hexane", "hexane", "hexane (column eluent)", ev(6), "CCCCCC", "achiral"),
    molecule("m_ethyl_acetate", "ethyl acetate", "ethyl acetate", ev(6), "CCOC(C)=O", "achiral"),
    molecule("m_silica", "silica gel", "silica gel", ev(6)),
    molecule("m_agsbf6", "AgSbF6", "silver hexafluoroantimonate", ev(4)),
    molecule("m_hyd_sub", "drawn 2-propynyl/propargyl acetate substrate", "1-phenylprop-2-yn-1-yl acetate (candidate from drawing)", ev(3) + ev(6),
             "CC(=O)OC(C#C)c1ccccc1", "unspecified",
             "The table/procedure calls the substrate 2-propynyl or propargyl acetate, while the drawing and 17.4 mg per 0.10 mmol show a phenyl-substituted substrate; identity remains pending review."),
    molecule("m_hyd_prod", "hydration product", "2-oxo-1-phenylpropyl acetate", ev(3) + ev(6),
             "CC(=O)OC(C(C)=O)c1ccccc1", "unspecified"),
    molecule("m_enyne", "1,6-enyne", "drawn dimethyl malonate-tethered 1,6-enyne (candidate)", ev(3) + ev(4) + ev(5) + ev(6),
             "C#CCC(C(=O)OC)(C(=O)OC)CC=C(C)C", "achiral",
             "Candidate connectivity was transcribed from the repeated drawing. The reported 25.3 mg/0.10 mmol charge and the drawing-derived candidate molecular mass require expert reconciliation."),
    molecule("m_cyclo_prod", "cycloisomerization product", "vinylcyclopentene cycloisomerization product (drawn; exact connectivity pending)", ev(3) + ev(5)),
    molecule("m_alkoxy_prod", "alkoxycyclization product", "methoxy-substituted methylenecyclopentene product (drawn; exact connectivity pending)", ev(4)),
    molecule("m_allenyl_alcohol", "allenyl alcohol", "2,2-diphenylhexa-4,5-dien-1-ol", ev(4) + ev(6),
             "OCC(c1ccccc1)(c1ccccc1)CC=C=C", "allene_achiral_partial"),
    molecule("m_vinyl_thf", "2-vinyltetrahydrofuran", "2-vinyl-4,4-diphenyltetrahydrofuran (candidate from drawing)", ev(4) + ev(6),
             "C=CC1CC(c2ccccc2)(c2ccccc2)CO1", "unspecified"),
]

# Preserve the 3b precursor conflict as an unresolved identity rather than
# silently choosing Scheme 2 or the experimental paragraph.
conflict = molecule("m_3b_precursor_conflict", "1a/1b conflict", "precursor reported inconsistently for preparation of 3b", ev(2) + ev(6))
conflict["occurrences"][0]["compound_label"] = f(
    ["Scheme 2: 1b", "Preparation of 3b: 1a"], None, ev(2) + ev(6), status="conflict",
    alternatives=[f("Scheme 2: 1b", "1b", ev(2)), f("Preparation of 3b: 1a", "1a", ev(6))])
conflict["names"]["original"] = conflict["occurrences"][0]["compound_label"]
conflict["names"]["english"] = missing()
molecules.append(conflict)
bundle["molecules"] = molecules


# Generic schemes shown over grouped rows/compound variants.
generic_defs = [
    ("scheme2_mono_generic", ev(2), "Scheme 2 Eq. (1)", "Synthesis of mono(chlorogold) complexes 2", "complexation"),
    ("scheme2_bis_generic", ev(2), "Scheme 2 Eq. (2)", "Synthesis of bis(chlorogold) complexes 3", "complexation"),
    ("table1_generic", ev(3), "Table 1 header scheme", "Gold-catalyzed hydration of 2-propynyl acetate", "hydration"),
    ("table2_generic", ev(3), "Table 2 header scheme", "Gold-catalyzed cycloisomerization of 1,6-enyne", "cycloisomerization"),
    ("table3_generic", ev(4), "Table 3 header scheme", "Gold-catalyzed alkoxycyclization of 1,6-enyne", "alkoxycyclization"),
    ("table4_generic", ev(4), "Table 4 header scheme", "Gold-catalyzed intramolecular hydroalkoxylation of allenyl alcohol", "hydroalkoxylation"),
    ("table5_generic", ev(4), "Table 5 header scheme", "Alkoxycyclization of 1,6-enyne catalyzed by mononuclear complexes 2", "alkoxycyclization"),
]
reactions = []
for rid, evidence, loc, name, kind in generic_defs:
    g = reaction(rid, evidence, record_type="generic_scheme", categories=["other"], location=loc,
                 name=name, transformation=kind, class_kind=kind)
    g["review"]["unresolved"].append("Generic scheme is retained separately from its individual compound/table-row experiments.")
    reactions.append(g)


def prep_record(label, precursor_mid, reagent_mass, product_mid, product_mass, product_yield, page, generic, bis=False, conflict_note=None):
    evidence = ev(2) + ev(page)
    r = reaction(f"prep_{label}", evidence, location=f"Scheme 2 and Preparation of {label}",
                 name=f"Preparation of {label}", transformation="chlorogold complexation",
                 general=[generic], class_kind="complexation")
    precursor_mass = {"2a": 71.5, "2b": 73.5, "3a": 71.5, "3b": 73.5}[label]
    r["participants"] = [
        participant("precursor", precursor_mid, "reactant", evidence, mass=q(precursor_mass, "mg"), amount=q(0.10, "mmol")),
        participant("tht_aucl", "m_tht_aucl", "reagent", evidence, mass=q(reagent_mass, "mg"), amount=q(0.20 if bis else 0.10, "mmol"),
                    equivalents=q(2 if bis else 1, "equiv"), reference=precursor_mid),
        participant("dcm_initial", "m_dcm", "solvent", evidence, volume=q(1, "mL")),
        participant("nhexane_wash", "m_nhexane", "purification_agent", ev(page), stage="purification", volume=q(5, "mL")),
        participant("product", product_mid, "product", evidence, product_role="main product"),
    ]
    if bis:
        r["participants"].append(participant("dcm_redissolution", "m_dcm", "solvent", ev(page), volume=q(1, "mL")))
        r["steps"] = [
            step("reaction_1", "reaction", ev(page), ["precursor", "tht_aucl", "dcm_initial"], "stir", "room temperature", q(1, "h"), "inert atmosphere (nitrogen or argon)"),
            step("remove_volatiles_1", "workup", ev(page), ["precursor", "tht_aucl", "dcm_initial"], "remove volatile materials under reduced pressure"),
            step("redissolve", "reaction", ev(page), ["dcm_redissolution"], "dissolve residue in dichloromethane", atmosphere="nitrogen"),
            step("reaction_2", "reaction", ev(page), ["dcm_redissolution"], "stir", "room temperature", q(1, "h"), "nitrogen"),
            step("remove_volatiles_2", "workup", ev(page), ["dcm_redissolution"], "remove volatile materials under reduced pressure"),
            step("wash", "purification", ev(page), ["nhexane_wash", "product"], "wash residue with n-hexane"),
        ]
    else:
        r["steps"] = [
            step("reaction", "reaction", ev(page), ["precursor", "tht_aucl", "dcm_initial"], "stir", "room temperature", q(5, "h"), "inert atmosphere (nitrogen or argon)"),
            step("remove_volatiles", "workup", ev(page), ["dcm_initial"], "remove volatile materials in vacuo"),
            step("wash", "purification", ev(page), ["nhexane_wash", "product"], "wash residue with n-hexane"),
        ]
    r["reaction_scale"] = f("precursor 0.10 mmol", {"reference_molecule_id": precursor_mid, "amount": q(0.10, "mmol")}, ev(page))
    solvent_fields(r, evidence, composition="dichloromethane", total=q(1, "mL", "1 mL"))
    procedure_fields(r, ev(page), original=f"Preparation of {label} paragraph; see page {page}.")
    appearance = {"2a": "pink", "2b": "white", "3a": "pink", "3b": "pink"}[label]
    isolation_phrase = f"to afford {label} as a {appearance} solid ({product_mass} mg, {product_yield}%)"
    r["outcomes"] = [
        outcome("product_mass", "product_amount", [product_mid], ev(page), f"{product_mass} mg", q(product_mass, "mg"), role="main product"),
        outcome("isolated_yield", "yield", [product_mid], evidence, f"{product_yield}%", q(product_yield, "%"), basis=isolation_phrase, role="main product"),
    ]
    if conflict_note:
        r["notes"].append(f(conflict_note, conflict_note, evidence))
        r["review"]["unresolved"].append(conflict_note)
    return r


reactions += [
    prep_record("2a", "m_1a", 32.1, "m_2a", 61.1, 65, 5, "scheme2_mono_generic"),
    prep_record("2b", "m_1b", 32.1, "m_2b", 62.2, 64, 5, "scheme2_mono_generic"),
    prep_record("3a", "m_1a", 64.1, "m_3a", 98.2, 80, 6, "scheme2_bis_generic", bis=True),
    prep_record("3b", "m_3b_precursor_conflict", 64.1, "m_3b", 102.1, 85, 6, "scheme2_bis_generic", bis=True,
                conflict_note="Source conflict: Scheme 2 assigns R = Ph precursor 1b to 3b, but the Preparation of 3b paragraph says 1a (73.5 mg, 0.10 mmol); no value was silently selected."),
]


def catalytic_base(rid, evidence, location, name, transformation, generic, class_kind, catalyst_mid, catalyst_loading,
                   substrate_mid, product_mid=None, substrate_mass=None, substrate_amount=None):
    r = reaction(rid, evidence, location=location, name=name, transformation=transformation,
                 general=[generic], class_kind=class_kind)
    r["participants"] = []
    if substrate_mid:
        r["participants"].append(participant("substrate", substrate_mid, "reactant", evidence,
                                             mass=q(substrate_mass, "mg") if substrate_mass is not None else None,
                                             amount=q(substrate_amount, "mmol") if substrate_amount is not None else None))
    r["participants"].append(participant("catalyst", catalyst_mid, "catalyst", evidence,
                                         mol_percent=q(catalyst_loading, "mol%"), reference=substrate_mid))
    if product_mid:
        r["participants"].append(participant("product", product_mid, "product", evidence, product_role="main product"))
    return r


def add_column_workup(r, evidence):
    r["participants"] += [
        participant("silica", "m_silica", "purification_agent", evidence, stage="purification"),
        participant("hexane_eluent", "m_hexane", "purification_agent", evidence, stage="purification"),
        participant("ethyl_acetate_eluent", "m_ethyl_acetate", "purification_agent", evidence, stage="purification"),
    ]
    r["steps"].append(step("column", "purification", evidence, ["silica", "hexane_eluent", "ethyl_acetate_eluent"],
                           "load reaction mixture directly on silica gel and elute with hexane/ethyl acetate 6:1"))


# Table 1: two hydration rows.
for catalyst, yld, rec in (("m_3a", 69, 31), ("m_3b", 60, 40)):
    label = catalyst[-2:]
    evidence = ev(3) + ev(6)
    r = catalytic_base(f"table1_hydration_{label}", evidence, f"Table 1 row {label}; hydration procedure",
                       "Gold-catalyzed hydration of 2-propynyl acetate", "terminal alkyne hydration",
                       "table1_generic", "hydration", catalyst, 3, "m_hyd_sub", "m_hyd_prod", 17.4, 0.10)
    r["participants"] += [
        participant("water", "m_water", "reagent", evidence, amount=q(1.0, "mmol"), equivalents=q(10, "equiv"), reference="m_hyd_sub"),
        participant("acetone", "m_acetone", "solvent", evidence, volume=q(1, "mL")),
    ]
    r["reaction_scale"] = f("substrate 0.10 mmol", {"reference_molecule_id": "m_hyd_sub", "amount": q(0.10, "mmol")}, ev(6))
    solvent_fields(r, evidence, composition="acetone", total=q(1, "mL", "1 mL"))
    r["steps"] = [step("reaction", "reaction", evidence, ["substrate", "catalyst", "water", "acetone"],
                       "add substrate solution and water to chlorogold complex; stir", "room temperature", q(48, "h"),
                       "inert atmosphere (nitrogen or argon)")]
    add_column_workup(r, ev(6))
    r["outcomes"] = [
        outcome("product_yield", "yield", ["m_hyd_prod"], ev(3), f"{yld}%", q(yld, "%"), "1H NMR spectroscopy", "main product", yield_type="NMR yield"),
        outcome("substrate_recovery", "recovery", ["m_hyd_sub"], ev(3), f"{rec}%", q(rec, "%"), "1H NMR spectroscopy", scope="recovered starting material"),
    ]
    procedure_fields(r, ev(6), original="Hydration of propargyl acetate procedure on page 6.")
    r["review"]["unresolved"].append("Substrate naming (2-propynyl/propargyl acetate) conflicts with the phenyl-substituted drawing; candidate identity needs expert review.")
    reactions.append(r)


# Table 2: two cycloisomerization rows.
for catalyst, yld, rec in (("m_3a", 95, 0), ("m_3b", 93, 7)):
    label = catalyst[-2:]
    evidence = ev(3) + ev(6)
    r = catalytic_base(f"table2_cycloisomerization_{label}", evidence, f"Table 2 row {label}; cycloisomerization procedure",
                       "Gold-catalyzed cycloisomerization of 1,6-enyne", "1,6-enyne cycloisomerization",
                       "table2_generic", "cycloisomerization", catalyst, 3, "m_enyne", "m_cyclo_prod", 25.3, 0.10)
    r["participants"].append(participant("dcm", "m_dcm", "solvent", evidence, volume=q(1, "mL")))
    r["reaction_scale"] = f("1,6-enyne 0.10 mmol", {"reference_molecule_id": "m_enyne", "amount": q(0.10, "mmol")}, ev(6))
    solvent_fields(r, evidence, composition="dichloromethane", total=q(1, "mL", "1 mL"))
    r["steps"] = [step("reaction", "reaction", evidence, ["substrate", "catalyst", "dcm"],
                       "add 1,6-enyne solution to chlorogold complex; stir", "room temperature", q(24, "h"),
                       "inert atmosphere (nitrogen or argon)")]
    add_column_workup(r, ev(6))
    r["outcomes"] = [
        outcome("product_yield", "yield", ["m_cyclo_prod"], ev(3), f"{yld}%", q(yld, "%"), "1H NMR spectroscopy", "main product", yield_type="NMR yield"),
        outcome("substrate_recovery", "recovery", ["m_enyne"], ev(3), f"{rec}%", q(rec, "%"), "1H NMR spectroscopy", scope="recovered 1,6-enyne"),
    ]
    procedure_fields(r, ev(6), original="Cycloisomerization of 1,6-enyne procedure on page 6.")
    r["review"]["unresolved"].append("Exact product connectivity and the 1,6-enyne candidate identity require expert structure review.")
    reactions.append(r)


# Table 3: two alkoxycyclization rows with catalyst recovery/decoordination outcomes.
for catalyst, yld, subrec, catrec, decoord in (("m_3a", 73, 23, 84, 16), ("m_3b", 97, 3, 67, 33)):
    label = catalyst[-2:]
    mono = "m_2a" if catalyst == "m_3a" else "m_2b"
    evidence = ev(4) + ev(6)
    r = catalytic_base(f"table3_alkoxycyclization_{label}", evidence, f"Table 3 row {label}; alkoxycyclization procedure",
                       "Gold-catalyzed alkoxycyclization of 1,6-enyne", "1,6-enyne alkoxycyclization",
                       "table3_generic", "alkoxycyclization", catalyst, 3, "m_enyne", "m_alkoxy_prod", 25.3, 0.10)
    r["participants"] += [
        participant("dcm", "m_dcm", "solvent", evidence, volume=q(0.5, "mL")),
        participant("methanol", "m_methanol", "reagent", evidence, volume=q(0.5, "mL")),
        participant("decoordination_product", mono, "product", evidence, product_role="catalyst de-coordination product"),
    ]
    r["reaction_scale"] = f("1,6-enyne 0.10 mmol", {"reference_molecule_id": "m_enyne", "amount": q(0.10, "mmol")}, ev(6))
    solvent_fields(r, evidence, composition=["dichloromethane", "methanol"], ratio="1:1")
    r["steps"] = [step("reaction", "reaction", evidence, ["substrate", "catalyst", "dcm", "methanol"],
                       "add 1,6-enyne solution in dichloromethane/methanol to chlorogold complex; stir", "room temperature", q(24, "h"),
                       "inert atmosphere (nitrogen or argon)")]
    add_column_workup(r, ev(6))
    r["outcomes"] = [
        outcome("product_yield", "yield", ["m_alkoxy_prod"], ev(4), f"{yld}%", q(yld, "%"), "1H NMR spectroscopy", "main product", yield_type="NMR yield"),
        outcome("substrate_recovery", "recovery", ["m_enyne"], ev(4), f"{subrec}%", q(subrec, "%"), "1H NMR spectroscopy", scope="recovered 1,6-enyne"),
        outcome("catalyst_recovery", "recovery", [catalyst], ev(4), f"{catrec}%", q(catrec, "%"), "31P NMR spectroscopy", scope=f"recovery of {label}"),
        outcome("decoordination_product_yield", "yield", [mono], ev(4), f"{decoord}%", q(decoord, "%"), "31P NMR spectroscopy", "catalyst de-coordination product", yield_type="NMR yield"),
    ]
    procedure_fields(r, ev(6), original="Alkoxycyclization of 1,6-enyne procedure on page 6.")
    r["review"]["unresolved"].append("Exact connectivity/stereochemistry of the drawn alkoxycyclization product remains unresolved.")
    reactions.append(r)


# Table 4: additive-free and AgSbF6 entries, with two separately introduced toluene portions.
for entry, additive, yld, rec in ((1, False, 11, 89), (2, True, 91, 9)):
    evidence = ev(4) + ev(6)
    r = catalytic_base(f"table4_hydroalkoxylation_entry{entry}", evidence, f"Table 4 entry {entry}; intramolecular hydroalkoxylation procedure",
                       "Gold-catalyzed intramolecular hydroalkoxylation of allenyl alcohol", "intramolecular allene hydroalkoxylation",
                       "table4_generic", "hydroalkoxylation", "m_3b", 5, "m_allenyl_alcohol", "m_vinyl_thf", 25.3, 0.10)
    r["participants"] += [
        participant("toluene_premix", "m_toluene", "solvent", ev(6), volume=q(0.4, "mL")),
        participant("toluene_substrate_solution", "m_toluene", "solvent", ev(6), volume=q(0.6, "mL")),
        participant("silica", "m_silica", "purification_agent", ev(6), stage="purification"),
        participant("hexane_eluent", "m_hexane", "purification_agent", ev(6), stage="purification"),
        participant("ethyl_acetate_eluent", "m_ethyl_acetate", "purification_agent", ev(6), stage="purification"),
    ]
    if additive:
        r["participants"].append(participant("additive", "m_agsbf6", "additive", evidence,
                                             mol_percent=q(10, "mol%"), reference="m_allenyl_alcohol"))
    r["reaction_scale"] = f("allenyl alcohol 0.10 mmol", {"reference_molecule_id": "m_allenyl_alcohol", "amount": q(0.10, "mmol")}, ev(6))
    solvent_fields(r, evidence, composition="toluene", basis="0.4 mL catalyst/additive premix plus 0.6 mL substrate solution; no author-reported measured total volume")
    premix_ids = ["catalyst", "toluene_premix"] + (["additive"] if additive else [])
    r["steps"] = [
        step("premix", "reaction", ev(6), premix_ids, "stir catalyst and additive (if present) in toluene", "room temperature", q(10, "min"),
             "inert atmosphere (nitrogen or argon)"),
        step("reaction", "reaction", ev(6), ["substrate", "product", "toluene_substrate_solution"] + premix_ids,
             "add allenyl alcohol solution in toluene to the premix; stir", "room temperature", q(24, "h"),
             "inert atmosphere (nitrogen or argon)"),
        step("column", "purification", ev(6), ["silica", "hexane_eluent", "ethyl_acetate_eluent"], "load directly on a short silica gel column and elute with hexane/ethyl acetate (ratio unreported)"),
    ]
    r["outcomes"] = [
        outcome("product_yield", "yield", ["m_vinyl_thf"], ev(4), f"{yld}%", q(yld, "%"), "1H NMR spectroscopy", "main product", yield_type="NMR yield"),
        outcome("substrate_recovery", "recovery", ["m_allenyl_alcohol"], ev(4), f"{rec}%", q(rec, "%"), "1H NMR spectroscopy", scope="recovered allenyl alcohol"),
    ]
    procedure_fields(r, ev(6), original="Intramolecular hydroalkoxylation of allenyl alcohol procedure on page 6.",
                     addition="stir 3b and additive in 0.4 mL toluene for 10 min, then add substrate in 0.6 mL toluene")
    r["review"]["unresolved"].append("Product tetrahedral configuration is not specified; candidate SMILES intentionally leaves it unassigned.")
    reactions.append(r)


# Table 5: seven distinct entries.  Entry 7 is a real negative control and does
# not invent an obtained product participant.
table5 = [
    (1, "m_2a", 1, "room temperature", 19, 81),
    (2, "m_2b", 1, "room temperature", 15, 85),
    (3, "m_2a", 6, "room temperature", 62, 38),
    (4, "m_2b", 6, "room temperature", 95, "<5"),
    (5, "m_2b", 1, q(40, "°C"), 25, 75),
    (6, "m_2b", 3, q(40, "°C"), 80, 20),
    (7, "m_4", 1, "room temperature", 0, ">99"),
]
for entry, catalyst, days, temp, yld, rec in table5:
    evidence = ev(4) + ev(6)
    product_mid = None if entry == 7 else "m_alkoxy_prod"
    r = catalytic_base(f"table5_entry{entry}", evidence, f"Table 5 entry {entry}; alkoxycyclization procedure",
                       "Alkoxycyclization of 1,6-enyne catalyzed by mononuclear complex", "1,6-enyne alkoxycyclization",
                       "table5_generic", "none" if entry == 7 else "alkoxycyclization", catalyst, 3, "m_enyne", product_mid, 25.3, 0.10)
    r["experiment_categories"] = f(["control"] if entry == 7 else ["optimization"], ["control"] if entry == 7 else ["optimization"], evidence)
    r["participants"] += [
        participant("dcm", "m_dcm", "solvent", evidence, volume=q(0.5, "mL")),
        participant("methanol", "m_methanol", "reagent", evidence, volume=q(0.5, "mL")),
    ]
    r["reaction_scale"] = f("1,6-enyne 0.10 mmol", {"reference_molecule_id": "m_enyne", "amount": q(0.10, "mmol")}, ev(6))
    solvent_fields(r, evidence, composition=["dichloromethane", "methanol"], ratio="1:1")
    r["steps"] = [step("reaction", "reaction", evidence, [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"],
                       "stir under the Table 5 conditions", temp, q(days, "d"), "inert atmosphere (nitrogen or argon)")]
    add_column_workup(r, ev(6))
    if entry == 7:
        r["outcomes"] = [
            outcome("no_reaction", "no_reaction", ["m_enyne"], ev(4), "0% product yield", True, "1H NMR spectroscopy", scope="no alkoxycyclization product observed in Table 5"),
            outcome("substrate_recovery", "recovery", ["m_enyne"], ev(4), ">99%", q(99, "%", ">99%", "gt"), "1H NMR spectroscopy", scope="recovered 1,6-enyne"),
        ]
        r["review"]["unresolved"].append("Negative control has no obtained product participant; 0% yield and >99% starting-material recovery are retained as no-reaction evidence.")
    else:
        yraw = f"{yld}%"
        yval = q(yld, "%")
        if rec == "<5":
            rraw, rval = "<5%", q(5, "%", "<5%", "lt")
        else:
            rraw, rval = f"{rec}%", q(rec, "%")
        r["outcomes"] = [
            outcome("product_yield", "yield", ["m_alkoxy_prod"], ev(4), yraw, yval, "1H NMR spectroscopy", "main product", yield_type="NMR yield"),
            outcome("substrate_recovery", "recovery", ["m_enyne"], ev(4), rraw, rval, "1H NMR spectroscopy", scope="recovered 1,6-enyne"),
        ]
        r["review"]["unresolved"].append("Exact connectivity/stereochemistry of the drawn alkoxycyclization product remains unresolved.")
    procedure_fields(r, ev(6), original="Alkoxycyclization general procedure supplies substrate charge, solvent portions, column workup; Table 5 supplies the entry-specific time/temperature.",
                     general_ref="Alkoxycyclization of 1,6-enyne procedure on page 6 with Table 5 modifications")
    reactions.append(r)


# Scheme 3 is a separate explicitly yielded experiment; its scale, solvent
# volumes, workup and measurement basis are not reported on pages 1-6.
evidence = ev(5)
r = catalytic_base("scheme3_tfe_cycloisomerization", evidence, "Scheme 3", "Cycloisomerization of 1,6-enyne promoted by 2b in 2,2,2-trifluoroethan-1-ol",
                   "1,6-enyne cycloisomerization", "table2_generic", "cycloisomerization", "m_2b", 3, "m_enyne", "m_cyclo_prod")
r["participants"] += [
    participant("dcm", "m_dcm", "solvent", evidence),
    participant("tfe", "m_tfe", "additive", evidence),
]
solvent_fields(r, evidence, composition=["dichloromethane", "2,2,2-trifluoroethan-1-ol"], ratio="1:1")
r["steps"] = [step("reaction", "reaction", evidence, ["substrate", "catalyst", "dcm", "tfe", "product"],
                   "stir under Scheme 3 conditions", "room temperature", q(24, "h"))]
r["outcomes"] = [outcome("product_yield", "yield", ["m_cyclo_prod"], evidence, "96%", q(96, "%"), role="main product")]
r["review"]["unresolved"] += ["Scheme 3 does not report reaction scale, solvent volumes, isolation/workup, or yield measurement basis on pages 1-6.",
                               "Exact product connectivity remains pending expert structure review."]
reactions.append(r)

bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": f"{DOC}:scheme1_prior_work", "document_id": DOC, "kind": "cited_prior_work",
     "evidence_ids": ev(1), "summary": "Scheme 1 is explicitly described as the authors' previous report [14]: P=P bond cleavage of Mes*P=PMes* to phosphanoxy-substituted phosphaalkene 1. It is not counted as a current-paper specific experiment."},
    {"context_id": f"{DOC}:figures1_3_characterization", "document_id": DOC, "kind": "characterization_and_computation",
     "evidence_ids": ev(2) + ev(3), "summary": "Figures 1-3 contain ORTEP/AIM, qualitative crystal structures, 31P NMR and DFT/conformational observations. These are retained as source context, not synthetic experiments."},
    {"context_id": f"{DOC}:figure4_proposed_interaction", "document_id": DOC, "kind": "proposed_interaction",
     "evidence_ids": ev(4), "summary": "Figure 4 is a plausible methanol/chloride interaction proposed to increase Lewis acidity; it is not represented as an experimentally proven mechanism."},
    {"context_id": f"{DOC}:unresolved_negative_statements", "document_id": DOC, "kind": "unresolved_experimental_statements",
     "evidence_ids": ev(4), "summary": "The prose reports no hydration/cycloisomerization with complexes 2 and no allenyl-alcohol hydroalkoxylation with 2b/AgSbF6, but variant-specific conditions/results are absent or deferred to unavailable SI; no fabricated specific records were created."},
    {"context_id": f"{DOC}:characterization_sections", "document_id": DOC, "kind": "compound_characterization",
     "evidence_ids": ev(5) + ev(6), "summary": "NMR, HRMS and crystallographic data for 2a, 2b, 3a and 3b are characterization context attached to their preparations; standalone spectral peak lists are not separate reactions."},
    {"context_id": f"{DOC}:page7_references_and_publication_history", "document_id": DOC, "kind": "bibliographic_context",
     "evidence_ids": ev(7), "summary": "PDF page 7 contains only the continuation of references [5]-[28] and the manuscript received, revised, accepted-online, and version-of-record dates. It contains no current-paper reaction arrow, structure scheme, experimental procedure, table row, yield, selectivity, or standalone characterization peak table, so no molecule or reaction record is created from this page."},
]

st.write(RUN / "submitted.json", bundle)

inventory = {
    "document_id": DOC,
    "title": "Gold(I) Complexation of Phosphanoxy-Substituted Phosphaalkenes for Activation-Free LAuCl Catalysis",
    "pages_visually_inspected": [1, 2, 3, 4, 5, 6],
    "specific_experiments": [r["reaction_id"] for r in reactions if r["record_type"] == "specific_experiment"],
    "generic_schemes": [r["reaction_id"] for r in reactions if r["record_type"] == "generic_scheme"],
    "source_context_ids": [x["context_id"] for x in bundle["source_context"]],
    "model_status": "needs_review",
    "limitations": ["No SI file was supplied or associated.", "Pages 1-6 were inspected; page 7 was used only as bibliographic context and is not claimed covered.",
                    "Organogold structures and two drawn cyclization products remain unresolved; candidate organic SMILES are not structure verification."]
}
st.write(RUN / "raw_model/visual_inventory.json", inventory)
print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
                  "specific": len(inventory["specific_experiments"]), "generic": len(inventory["generic_schemes"])}, ensure_ascii=False))
