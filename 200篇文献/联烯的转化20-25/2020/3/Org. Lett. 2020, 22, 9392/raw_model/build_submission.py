#!/usr/bin/env python3
"""Build the source-derived v3.1 submission for Org. Lett. 2020, 22, 9392.

This generator contains only document-specific extraction decisions. It does not
modify the extraction tools, schemas, tests, or scientific quality rules.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st
import reaction_tools


RUN = ROOT / "runs/current"
DATASET_PATH = RUN / "dataset.json"
DOC = "doc_c17630b33b36c6c9"
data = st.read(DATASET_PATH)


def ev(page: int, image: bool = True, text: bool = True) -> list[str]:
    result = []
    if image:
        result.append(f"{DOC}:page{page}:image")
    if text:
        result.append(f"{DOC}:page{page}:text")
    return result


def f(raw=None, normalized=None, evidence=(), method="original_extraction",
      status=None, derivation=None, alternatives=()):
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   derivation=derivation, alternatives=alternatives)


def missing(status="unreported"):
    method = "not_applicable" if status == "not_applicable" else "original_extraction"
    return f(method=method, status=status)


def q(number, unit, **extra):
    result = {"number": number, "unit": unit}
    result.update(extra)
    return result


def mid(label: str) -> str:
    return f"{DOC}:mol:{label}"


bundle = st.empty(data["documents"])
molecules: dict[str, dict] = {}
stereo_checks = {}


def add_molecule(label: str, evidence, name=None, *, chirality="unresolved",
                 geometry="unresolved", representation="not_encoded",
                 atom_centers=None, double_bonds=None, axial=None, mixture=None,
                 assignment=None, note=None, smiles=None) -> str:
    identity = mid(label)
    if identity in molecules:
        return identity
    m = st.new_molecule(identity, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    m["names"]["original"] = f(name or label, name or label, evidence)
    m["names"]["english"] = f(name or label, name or label, evidence)
    m["names"]["abbreviations"] = f(label, [label], evidence)
    m["names"]["aliases"] = missing()
    m["cas_numbers"] = missing()
    m["identity_note"] = f(note, note, evidence) if note else missing()
    m["external_properties"] = []
    if smiles:
        diag = st.stereo_inspect(smiles)
        stereo_checks[label] = diag
        deriv = "RDKit graph diagnostics on a source-named, machine-entered candidate; parsing does not certify source identity."
        m["representation_type"] = "smiles"
        m["smiles"] = f(smiles, diag["canonical_isomeric_smiles"], evidence,
                          method="structure_derivation", derivation=deriv)
        m["generic_template"] = None
        m["stereochemistry"]["atom_centers"] = f(None, diag["atom_centers"], evidence,
            method="structure_derivation", derivation=deriv)
        m["stereochemistry"]["double_bonds"] = f(None, diag["double_bonds"], evidence,
            method="structure_derivation", derivation=deriv)
        m["stereochemistry"]["axial"] = f(None, [], evidence,
            method="structure_derivation", derivation=deriv)
        m["stereochemistry"]["mixture"] = f(None, [], evidence,
            method="structure_derivation", derivation=deriv)
        m["stereochemistry"]["assignment"] = f(None, [], evidence,
            method="structure_derivation", derivation=deriv)
        m["stereochemistry"]["chirality_status"] = f(None, "achiral", evidence,
            method="structure_derivation", derivation=deriv)
        m["stereochemistry"]["geometry_status"] = f(None, "not_applicable", evidence,
            method="structure_derivation", derivation=deriv)
        m["stereochemistry"]["representation_status"] = f(None, "encoded", evidence,
            method="structure_derivation", derivation=deriv)
        m["stereochemistry"]["isomeric_smiles"] = f(None, diag["canonical_isomeric_smiles"], evidence,
            method="structure_derivation", derivation=deriv)
        m["stereochemistry"]["identity_key"] = f(None, diag["stereo_graph_key"], evidence,
            method="structure_derivation", derivation=deriv)
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing()
        m["generic_template"] = None
        values = {
            "atom_centers": atom_centers,
            "double_bonds": double_bonds,
            "axial": axial,
            "mixture": mixture,
            "assignment": assignment,
        }
        for key, value in values.items():
            m["stereochemistry"][key] = f(value, value, evidence) if value is not None else missing()
        m["stereochemistry"]["chirality_status"] = f(chirality, chirality, evidence)
        m["stereochemistry"]["geometry_status"] = f(geometry, geometry, evidence)
        m["stereochemistry"]["representation_status"] = f(representation, representation, evidence)
        m["stereochemistry"]["isomeric_smiles"] = missing()
        m["stereochemistry"]["identity_key"] = missing()
    molecules[identity] = m
    return identity


def add_template(label: str, evidence, name: str, fixed: str, variables: list[dict]) -> str:
    identity = add_molecule(label, evidence, name=name, chirality="unresolved",
                            geometry="unresolved", representation="not_encoded",
                            note="Generic source drawing; R/Ar variables are not expanded as a concrete molecule.")
    m = molecules[identity]
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": f(fixed, {"source_description": fixed}, evidence),
        "variables": [
            {
                "label": row["label"],
                "attachment_points": f(row["attachment"], row["attachment"], evidence),
                "definition": f(row["definition"], row["definition"], evidence),
                "allowed_values": f(row.get("allowed"), row.get("allowed"), evidence) if row.get("allowed") else missing(),
                "constraints": missing(),
            }
            for row in variables
        ],
        "instance_mapping": [],
        "reason": f("Generic R/Ar drawing in the source", "Generic template required; no wildcard SMILES asserted.", evidence),
    }
    return identity


def add_complex_product(label: str, evidence, *, only_major=False, special=None, name=None):
    if special == "multi_isomer":
        chirality = "mixture"
        mixture = "At least three pairs of isomers; source attributes extra axial chirality to the 2,4-dichlorophenyl group."
        assignment = "Axial and central configurations are not assigned."
    elif only_major:
        chirality = "partially_specified"
        mixture = "Only the major diastereomer was obtained; absolute configuration is unassigned."
        assignment = "Major relative isomer shown; absolute configuration not reported."
    else:
        chirality = "mixture"
        mixture = "Combined anti/syn diastereomeric product mixture unless otherwise noted; ratio is experiment-specific."
        assignment = "Major isomer is anti and minor is syn; absolute configurations are unassigned."
    return add_molecule(label, evidence, name=name or f"tetrahydroisoquinoline {label}",
        chirality=chirality, geometry="not_applicable", representation="not_encoded",
        atom_centers="Source drawing contains stereogenic ring centers; atom-level assignments are not encoded.",
        mixture=mixture, assignment=assignment,
        note="Dense source structure retained by compound label; complete stereochemical SMILES was not safely established.")


# Generic current-work structures and page-1 prior-work placeholders.
g1 = add_template("generic_1", ev(1), "cross-conjugated azatriene 1",
    "N-tosyl cross-conjugated azatriene with two aryl/alkyl substituent positions",
    [{"label": "Ar/R", "attachment": "two terminal carbon substituent positions", "definition": "source-variable aryl or alkyl groups"}])
g2 = add_template("generic_2", ev(1), "delta-sulfonamido-allenoate 2",
    "delta-sulfonamido-allenoate with ester and delta-substituent variables",
    [{"label": "R1", "attachment": "delta carbon", "definition": "aryl or alkyl group"},
     {"label": "R2", "attachment": "ester alkyl group", "definition": "ethyl, isopropyl, or tert-butyl in displayed examples"}])
g3 = add_template("generic_3", ev(1), "tetrahydroisoquinoline product 3",
    "fused tetrahydroisoquinoline bearing ester, aryl/alkyl, phenyl, and N-tosyl substituents",
    [{"label": "Ar/R", "attachment": "displayed ring substituent positions", "definition": "source-variable aryl or alkyl groups"},
     {"label": "R2", "attachment": "ester alkyl group", "definition": "source-variable ester substituent"}])
for tag in ("prior_a_reactant", "prior_a_product", "prior_b_reactant", "prior_b_product",
            "prior_c_reactant", "prior_c_product"):
    add_molecule(tag, ev(1), name=tag.replace("_", " "), chirality="unresolved",
                 geometry="unresolved", representation="not_encoded",
                 note="Generic structure shown only as a prior-work example in Scheme 1.")


# Source-labelled substrates and products.
for letter in list("abcdefghijklmnopqrst") + ["y", "z"]:
    add_molecule(f"1{letter}", ev(2), name=f"cross-conjugated azatriene 1{letter}",
        chirality="unspecified", geometry="unspecified", representation="not_encoded",
        double_bonds="Cross-conjugated alkene/imine geometry is drawn but no E/Z assignment is stated.",
        note="Source-labelled azatriene; complete structure is retained by the drawing and label, not a guessed SMILES.")
for letter in list("abcdefghijklmn"):
    evidence = ev(2) if letter == "a" else ev(3)
    add_molecule(f"2{letter}", evidence, name=f"delta-sulfonamido-allenoate 2{letter}",
        chirality="unspecified", geometry="not_applicable", representation="not_encoded",
        axial="Potential allene axial stereochemistry is not assigned in the source.",
        note="Source-labelled allenoate; no axial configuration or complete SMILES is asserted.")

scheme2_products = [f"3{x}a" for x in "abcdefghijklmnopqrst"] + ["3ya"]
for label in scheme2_products:
    add_complex_product(label, ev(2), only_major=(label == "3ta"))
scheme3_success = [f"3a{x}" for x in "bcdefghijkl"]
for label in scheme3_success:
    add_complex_product(label, ev(3), special="multi_isomer" if label == "3ag" else None)

# Mechanistic and derivatization structures.
add_molecule("3u", ev(3), name="cross-conjugated azatriene 3u (4-bromobenzenesulfonyl analog)",
    chirality="unspecified", geometry="unspecified", representation="not_encoded",
    note="Source drawing retained; the leaving sulfonamide-derived group differs from the N-tosyl group in product 3aa.")
add_molecule("D-3ab", ev(3), name="deuterium-labelled product D-3ab", chirality="mixture",
    geometry="not_applicable", representation="not_encoded",
    atom_centers="Tetrahydroisoquinoline stereogenic centers are not atom-level assigned.",
    mixture="Diastereomeric/isotopologue composition is not fully specified.",
    assignment="Deuterium incorporation is shown at four positions; exact atom mapping remains unresolved.")
add_molecule("11", ev(3), name="imine 11", chirality="unspecified", geometry="unspecified",
    representation="not_encoded", note="Source-drawn imine; complete structure candidate not asserted.")
add_molecule("12", ev(3), name="benzannulation product 12", chirality="unspecified",
    geometry="not_applicable", representation="not_encoded")
for label, name, evidence, chirality in [
    ("5", "alcohol derivative 5", ev(4), "partially_specified"),
    ("6na", "dihydroisoquinoline 6na", ev(4), "unspecified"),
    ("7na", "isoquinoline 7na", ev(4), "unspecified"),
    ("8", "carboxylic-acid intermediate 8", ev(4), "unspecified"),
    ("9", "fused pentacyclic product 9", ev(4), "unspecified"),
    ("6ma", "dihydroisoquinoline 6ma", ev(4), "unspecified"),
    ("10", "fused tetracyclic product 10", ev(4), "unspecified"),
]:
    add_molecule(label, evidence, name=name, chirality=chirality,
        geometry="not_applicable", representation="not_encoded",
        assignment="Relative/absolute configuration is not fully assigned in the main-text drawing.",
        note="Dense source structure retained by label; complete stereochemical SMILES not safely established.")


# Reagents and solvents. Candidate SMILES are checked but remain machine-entered.
simple = {
    "P(p-Tol)3": ("tri(p-tolyl)phosphine", "Cc1ccc(P(c2ccc(C)cc2)c2ccc(C)cc2)cc1", ev(2)),
    "CH2Cl2": ("dichloromethane", "ClCCl", ev(2)),
    "PCy3": ("tricyclohexylphosphine", "P(C1CCCCC1)(C1CCCCC1)C1CCCCC1", ev(2, text=True)),
    "CHCl3": ("chloroform", "ClC(Cl)Cl", ev(2, image=False)),
    "D2O": ("deuterium oxide", "[2H]O[2H]", ev(3)),
    "THF": ("tetrahydrofuran", "C1CCOC1", ev(3)),
    "DDQ": ("2,3-dichloro-5,6-dicyano-1,4-benzoquinone", "N#CC1=C(Cl)C(=O)C(Cl)=C(C#N)C1=O", ev(4)),
    "toluene": ("toluene", "Cc1ccccc1", ev(4)),
    "KOH": ("potassium hydroxide", "[K+].[OH-]", ev(4)),
    "EtOH": ("ethanol", "CCO", ev(4)),
    "TfOH": ("trifluoromethanesulfonic acid", "OS(=O)(=O)C(F)(F)F", ev(4)),
    "BBr3": ("boron tribromide", "B(Br)(Br)Br", ev(4)),
}
for label, (name, smiles, evidence) in simple.items():
    add_molecule(label, evidence, name=name, smiles=smiles,
                 note="Machine-entered candidate for a source-named reagent/solvent; RDKit parsing is not source certification.")
for label, name, evidence in [
    ("K3PO4.3H2O", "tripotassium phosphate trihydrate", ev(2)),
    ("K2CO3", "potassium carbonate", ev(2)),
    ("LiAlH4", "lithium aluminium hydride", ev(4)),
    ("SmI2", "samarium(II) iodide", ev(4)),
]:
    add_molecule(label, evidence, name=name, chirality="not_applicable",
                 geometry="not_applicable", representation="not_applicable",
                 note="Source-named reagent; salt/solvation representation left unresolved rather than guessed.")


def participant(pid, molecule, role, evidence, *, stage="reaction", product_role=None,
                amount=None, volume=None, equivalents=None, mol_percent=None, reference=None,
                method="original_extraction"):
    p = st.new_participant(pid, molecule, role, evidence, stage=stage)
    p["product_role"] = (f(product_role, product_role, evidence, method=method)
                         if product_role else missing("not_applicable"))
    for key, value in (("amount", amount), ("volume", volume),
                       ("equivalents", equivalents), ("mol_percent", mol_percent)):
        if value is not None:
            p["charges"][key] = f(value[0], value[1], evidence, method=method)
    if reference:
        p["equivalence_reference"] = f(reference[0], reference[1], evidence, method=method)
    else:
        p["equivalence_reference"] = missing("not_applicable")
    return p


def reaction_base(suffix, evidence, record_type="specific_experiment", categories=None,
                  location=None, name=None, transformation=None):
    r = st.new_reaction(f"{DOC}:{suffix}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = (f(categories or [], categories or [], evidence)
                                    if record_type == "specific_experiment" else missing("not_applicable"))
    r["location"] = f(location or suffix, location or suffix, evidence)
    r["reaction_name"] = f(name, name, evidence) if name else missing()
    r["transformation_type"] = f(transformation, transformation, evidence) if transformation else missing()
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    return r


def add_step(r, evidence, *, operation=None, temperature=None, time=None, participants=None,
             method="original_extraction", step_id="reaction"):
    s = st.new_step(step_id, evidence)
    s["participant_ids"] = participants or [p["participant_id"] for p in r["participants"]]
    if operation:
        s["operation"] = f(operation, operation, evidence, method=method)
    if temperature:
        s["temperature"] = f(temperature[0], temperature[1], evidence, method=method)
    if time:
        s["time"] = f(time[0], time[1], evidence, method=method)
    r["steps"].append(s)


def add_outcome(r, oid, kind, targets, evidence, raw, normalized, *, product_role=None,
                yield_type=None, basis=None, scope=None, method="original_extraction",
                alternatives=()):
    o = st.new_outcome(oid, kind, targets, evidence)
    o["value"] = f(raw, normalized, evidence, method=method, alternatives=alternatives,
                   status="conflict" if alternatives else None)
    o["product_role"] = (f(product_role, product_role, evidence, method=method)
                         if product_role else missing("not_applicable"))
    o["yield_type"] = (f(yield_type, yield_type, evidence, method=method)
                       if yield_type else missing("not_applicable" if kind not in ("yield", "recovery") else "unreported"))
    o["measurement_basis"] = (f(basis, basis, evidence, method=method) if basis else missing())
    o["scope"] = f(scope, scope, evidence, method=method) if scope is not None else missing()
    r["outcomes"].append(o)


def set_solvent(r, evidence, molecule, raw_name, *, volume=None, method="original_extraction", dry=False):
    p = participant("p_solvent", molecule, "solvent", evidence,
                    volume=(f"{volume} mL", q(volume, "mL")) if volume is not None else None,
                    method=method)
    r["participants"].append(p)
    text = ("dry " if dry else "") + raw_name
    r["solvent_conditions"]["composition"] = f(text, [text], evidence, method=method)
    if volume is not None:
        r["solvent_conditions"]["total_volume"] = f(f"{volume} mL", q(volume, "mL"), evidence, method=method)


def set_classification(r, evidence, changes=None, coverage="partial", limitation=None):
    c = r["classification"]
    c["mapped_reaction_smiles"] = missing()
    if coverage in ("unresolved", "not_applicable"):
        c["bond_changes"] = []
        c["labels"] = (f(None, [], evidence, method="inference",
                         derivation="No observed transformation is classified." )
                       if coverage == "not_applicable" else missing())
        c["coverage"] = f(None, coverage, evidence, method="inference",
                          derivation="No net observed transformation" if coverage == "not_applicable" else "Source structures are insufficient for a defensible net classification.")
        c["limitations"] = [limitation or "No classification asserted."]
        c["model_votes"] = []
        return
    norm_changes = []
    labels = {}
    for idx, row in enumerate(changes or [], 1):
        change = {"change_id": f"b{idx}", "atoms": [{"atom_map": None, "element": row[0]},
                  {"atom_map": None, "element": row[1]}], "before_order": row[2],
                  "after_order": row[3], "change_type": row[4]}
        if len(row) > 5:
            change["hydrogen_count"] = row[5]
        norm_changes.append(f(None, change, evidence, method="structure_derivation",
                              derivation="Partial net change inferred from the source drawings; atom correspondence is unresolved."))
        if row[4] in ("formed", "order_increased", "hydrogen_gain"):
            lab = st.bond_label(row[0], row[1], row[3])
            labels.setdefault(lab["code"], dict(lab, supporting_change_ids=[]))["supporting_change_ids"].append(f"b{idx}")
    c["bond_changes"] = norm_changes
    c["labels"] = f(None, list(labels.values()), evidence, method="structure_derivation",
                    derivation="Taxonomy labels cover only visually defensible new connections/order increases in the source drawings.")
    c["coverage"] = f(None, coverage, evidence, method="inference",
                      derivation="Complete atom mapping, hydrogen accounting, leaving groups, and coproducts are unresolved.")
    c["limitations"] = [limitation or "Partial heavy-atom classification only; no complete atom mapping or hydrogen/coproduct balance."]
    c["model_votes"] = []


ANNULATION_CHANGES = [("C", "C", 0, 1, "formed"), ("C", "C", 0, 1, "formed"),
                      ("C", "N", 0, 1, "formed")]


reactions = []

# Page-1 generic current reaction (also covers the graphical abstract) and prior work.
r = reaction_base("scheme1_current_generic", ev(1), "generic_scheme", location="Page 1 graphical abstract and Scheme 1d",
                  name="Phosphine-catalyzed sequential [3+3]/aza-6pi-electrocyclization",
                  transformation="sequential annulation to tetrahydroisoquinolines")
r["participants"] = [participant("p_1", g1, "reactant", ev(1)), participant("p_2", g2, "reactant", ev(1)),
                     participant("p_3", g3, "product", ev(1), product_role="target product")]
set_classification(r, ev(1), ANNULATION_CHANGES)
r["review"]["unresolved"] = ["Generic R/Ar structures and exact atom correspondence are unresolved."]
reactions.append(r)
for key, mode in (("a", "[5+1] annulation"), ("b", "[3+2] annulation"), ("c", "[4+2] prior-work reactions")):
    rr = reaction_base(f"scheme1_prior_{key}", ev(1), "literature_example",
                       location=f"Page 1 Scheme 1{key} (prior work)", name=mode,
                       transformation="prior literature example")
    rr["participants"] = [participant("p_reactant", mid(f"prior_{key}_reactant"), "reactant", ev(1)),
                          participant("p_product", mid(f"prior_{key}_product"), "product", ev(1), product_role="literature-example product")]
    set_classification(rr, ev(1), coverage="unresolved",
                       limitation="Prior-work generic drawing is not classified as a current-paper experiment.")
    rr["review"]["unresolved"] = ["Bibliographic prior-work example; experimental details are outside the supplied current paper."]
    reactions.append(rr)

# Generic scope schemes.
for suffix, page, loc in (("scheme2_general", 2, "Page 2 Scheme 2"), ("scheme3_general", 3, "Page 3 Scheme 3")):
    rr = reaction_base(suffix, ev(page), "generic_scheme", location=loc,
        name="Phosphine-catalyzed sequential [3+3]/aza-6pi-electrocyclization",
        transformation="substrate-scope annulation")
    rr["participants"] = [participant("p_1", g1, "reactant", ev(page)), participant("p_2", g2, "reactant", ev(page)),
                          participant("p_3", g3, "product", ev(page), product_role="scope product")]
    set_classification(rr, ev(page), ANNULATION_CHANGES)
    rr["review"]["unresolved"] = ["Generic structures and full atom mapping unresolved."]
    reactions.append(rr)


def scope_annulation(product, substrate, allenoate, page, yield_value=None, dr=None,
                     dr_nmr=False, scheme="scheme2_general", only_major=False,
                     decomposition=False, multi_isomer=False):
    evidence = ev(page)
    rr = reaction_base(f"{scheme}_{product}", evidence, categories=["substrate_scope"],
        location=f"Page {page} {'Scheme 2' if page == 2 else 'Scheme 3'} entry {product}",
        name="sequential [3+3]/aza-6pi-electrocyclization",
        transformation="annulation to tetrahydroisoquinoline")
    rr["general_scheme_ids"] = [f"{DOC}:{scheme}"]
    inherited = page == 3
    method = "inherited_general_condition" if inherited else "original_extraction"
    condition_evidence = evidence + (ev(2) if inherited else [])
    s_mid, a_mid, p_mid = mid(substrate), mid(allenoate), mid(product)
    rr["participants"] = [
        participant("p_substrate", s_mid, "reactant", condition_evidence,
                    amount=("0.1 mmol", q(0.1, "mmol")), method=method),
        participant("p_allenoate", a_mid, "reactant", condition_evidence,
                    amount=("0.12 mmol", q(0.12, "mmol")), method=method),
        participant("p_catalyst", mid("P(p-Tol)3"), "catalyst", condition_evidence,
                    mol_percent=("20 mol%", q(20, "mol%")), reference=(substrate, s_mid), method=method),
        participant("p_base", mid("K3PO4.3H2O"), "reagent", condition_evidence,
                    mol_percent=("50 mol%", q(50, "mol%")), reference=(substrate, s_mid), method=method),
        participant("p_product", p_mid, "product", evidence, product_role="target product"),
    ]
    set_solvent(rr, condition_evidence, mid("CH2Cl2"), "CH2Cl2", volume=1.0, method=method)
    rr["reaction_scale"] = f("substrate 1: 0.1 mmol",
        {"reference_molecule_id": s_mid, "amount": q(0.1, "mmol")}, condition_evidence, method=method)
    rr["solvent_conditions"]["concentration"] = f(None, q(0.1, "M"), condition_evidence,
        method="inference", derivation="0.1 mmol limiting azatriene divided by 1.0 mL reported solvent volume.")
    rr["solvent_conditions"]["concentration_basis"] = f(None, substrate, condition_evidence,
        method="inference", derivation="Concentration calculated on the limiting azatriene.")
    add_step(rr, condition_evidence, operation="sealed tube immersed in an oil bath",
             temperature=("60 degrees C", q(60, "°C")), time=("24 h", q(24, "h")), method=method)
    rr["procedure"]["general_procedure_reference"] = f(
        "Reaction carried out under optimal conditions" if inherited else "Optimal conditions in Scheme 2 footnote a",
        "optimal conditions", condition_evidence, method=method)
    if yield_value is not None:
        add_outcome(rr, "yield", "yield", [p_mid], evidence, f"{yield_value}%", q(yield_value, "%"),
                    product_role="combined anti/syn product", yield_type="combined yield of diastereomeric isomers",
                    scope={"starting_material": substrate, "product": product})
    if dr is not None:
        add_outcome(rr, "dr", "selectivity", [p_mid], evidence, dr, dr,
                    product_role="major anti / minor syn",
                    basis="1H NMR" if dr_nmr else None,
                    scope={"metric": "dr", "order": ["anti", "syn"], "major": "anti", "minor": "syn"})
    if only_major:
        add_outcome(rr, "only_major", "selectivity", [p_mid], evidence,
                    "Only the major isomer was obtained", "only major isomer",
                    product_role="major isomer; absolute configuration unassigned",
                    scope={"metric": "diastereoselectivity", "minor": "not reported/detected"})
    if decomposition:
        add_outcome(rr, "decomposition", "decomposition", [p_mid], evidence,
                    "3ra decomposed during purification on silica gel", "decomposed during silica-gel purification",
                    product_role="unstable target product", scope={"stage": "purification"})
    if multi_isomer:
        add_outcome(rr, "isomer_mixture", "selectivity", [p_mid], evidence,
                    "mixture with at least three pairs of isomers", "at least three pairs of isomers",
                    product_role="isomer mixture", scope={"metric": "isomer composition", "axial_chirality": "unassigned"})
    set_classification(rr, evidence, ANNULATION_CHANGES)
    rr["review"]["unresolved"] = ["Complete substrate/product SMILES, atom mapping, absolute configuration, yield measurement basis, workup, and purification are not reported in the supplied main text."]
    if dr is not None and not dr_nmr:
        rr["review"]["unresolved"].append("dr measurement method is not stated for this entry.")
    if decomposition:
        rr["review"]["unresolved"].append("No stable isolated yield is reported because the product decomposed on silica gel.")
    if multi_isomer:
        rr["review"]["unresolved"].append("The identities and ratios of the at-least-three isomer pairs are not assigned.")
    reactions.append(rr)


# Scheme 2: 20 displayed scope entries.
scheme2_rows = [
    ("3aa", "1a", 84, "8.3:1", False), ("3ba", "1b", 80, "5.2:1", False),
    ("3ca", "1c", 76, "4.8:1", True), ("3da", "1d", 81, "6.8:1", False),
    ("3ea", "1e", 71, "5.4:1", False), ("3fa", "1f", 79, "7.4:1", False),
    ("3ga", "1g", 66, "5.2:1", False), ("3ha", "1h", 85, "8.3:1", True),
    ("3ia", "1i", 84, "5.2:1", False), ("3ja", "1j", 80, "7.2:1", False),
    ("3ka", "1k", 74, "9.1:1", True), ("3la", "1l", 92, "8.8:1", False),
    ("3ma", "1m", 77, "7.1:1", False), ("3na", "1n", 88, "6:1", False),
    ("3oa", "1o", 81, "5.9:1", True), ("3pa", "1p", 85, "4.8:1", False),
    ("3qa", "1q", 77, "7.2:1", False),
    ("3sa", "1s", 91, "1.2:1", True),
]
for product, substrate, yld, dr, nmr in scheme2_rows:
    scope_annulation(product, substrate, "2a", 2, yld, dr, nmr)
scope_annulation("3ra", "1r", "2a", 2, decomposition=True)
scope_annulation("3ta", "1t", "2a", 2, 47, only_major=True)

# Page-2 equations 1 and 2.
rr = reaction_base("eq1_3ya", ev(2), categories=["substrate_scope"], location="Page 2 eq 1",
    name="sequential [3+3]/aza-6pi-electrocyclization", transformation="aliphatic-substituted annulation")
rr["general_scheme_ids"] = [f"{DOC}:scheme2_general"]
rr["participants"] = [participant("p_1y", mid("1y"), "reactant", ev(2)),
    participant("p_2a", mid("2a"), "reactant", ev(2)),
    participant("p_cat", mid("P(p-Tol)3"), "catalyst", ev(2), mol_percent=("20 mol%", q(20, "mol%")), reference=("1y", mid("1y"))),
    participant("p_base", mid("K3PO4.3H2O"), "reagent", ev(2), mol_percent=("50 mol%", q(50, "mol%")), reference=("1y", mid("1y"))),
    participant("p_product", mid("3ya"), "product", ev(2), product_role="target product")]
set_solvent(rr, ev(2), mid("CH2Cl2"), "CH2Cl2")
add_step(rr, ev(2), temperature=("60 degrees C", q(60, "°C")), time=("24 h", q(24, "h")))
add_outcome(rr, "yield", "yield", [mid("3ya")], ev(2), "15%", q(15, "%"), product_role="combined product", yield_type="reported yield")
add_outcome(rr, "dr", "selectivity", [mid("3ya")], ev(2), "4.7:1", "4.7:1", product_role="major/minor unassigned", scope={"metric":"dr","order":"major:minor"})
set_classification(rr, ev(2), ANNULATION_CHANGES)
rr["review"]["unresolved"] = ["Substrate charges, solvent volume, dr/yield measurement methods, workup, and absolute configuration are unreported."]
reactions.append(rr)

rr = reaction_base("eq2_3za_ND", ev(2), categories=["substrate_scope", "control"], location="Page 2 eq 2",
    name="attempted sequential annulation", transformation="no detected annulation product")
rr["general_scheme_ids"] = [f"{DOC}:scheme2_general"]
rr["participants"] = [participant("p_1z", mid("1z"), "reactant", ev(2)),
    participant("p_2a", mid("2a"), "reactant", ev(2)),
    participant("p_cat", mid("P(p-Tol)3"), "catalyst", ev(2), mol_percent=("20 mol%", q(20, "mol%")), reference=("1z", mid("1z"))),
    participant("p_base", mid("K3PO4.3H2O"), "reagent", ev(2), mol_percent=("50 mol%", q(50, "mol%")), reference=("1z", mid("1z")))]
set_solvent(rr, ev(2), mid("CH2Cl2"), "CH2Cl2")
add_step(rr, ev(2), temperature=("60 degrees C", q(60, "°C")), time=("24 h", q(24, "h")))
add_outcome(rr, "no_product", "no_reaction", [], ev(2), "3za ND", "product not detected", scope={"expected_label":"3za"})
add_outcome(rr, "recovery", "recovery", [mid("1z")], ev(2), "67% recovered", q(67, "%"), product_role="recovered starting material", yield_type="recovery")
set_classification(rr, ev(2), coverage="not_applicable", limitation="No product detected; expected bonds are not classified.")
rr["review"]["unresolved"] = ["Detection and recovery measurement methods, substrate charges, and solvent volume are unreported."]
reactions.append(rr)

# Initial non-optimal experiment from prose.
rr = reaction_base("initial_1a_2b_3ab", ev(2), categories=["optimization"], location="Page 2 optimization prose",
    name="initial phosphine-catalyzed annulation", transformation="annulation to 3ab")
rr["general_scheme_ids"] = [f"{DOC}:scheme1_current_generic"]
rr["participants"] = [participant("p_1a", mid("1a"), "reactant", ev(2)), participant("p_2b", mid("2b"), "reactant", ev(2)),
    participant("p_cat", mid("PCy3"), "catalyst", ev(2)), participant("p_base", mid("K2CO3"), "reagent", ev(2)),
    participant("p_product", mid("3ab"), "product", ev(2), product_role="target product")]
set_solvent(rr, ev(2), mid("CHCl3"), "chloroform")
add_step(rr, ev(2), operation="reaction in chloroform; other conditions referred to unavailable Table S1")
add_outcome(rr, "yield", "yield", [mid("3ab")], ev(2), "24% yield", q(24, "%"), product_role="combined product", yield_type="reported yield")
add_outcome(rr, "dr", "selectivity", [mid("3ab")], ev(2), "10:1", "10:1", product_role="major/minor unassigned", scope={"metric":"dr","order":"major:minor"})
set_classification(rr, ev(2), ANNULATION_CHANGES)
rr["review"]["unresolved"] = ["Table S1 was not supplied; catalyst/base loadings, charges, temperature, time, yield/dr methods, workup, and purification are unavailable."]
reactions.append(rr)

# Gram-scale 3na experiment is distinct from the 0.1 mmol scope entry.
rr = reaction_base("gram_scale_3na", ev(2), categories=["scale_up", "substrate_scope"], location="Page 2 gram-scale prose and Scheme 2 parenthesis",
    name="gram-scale sequential annulation", transformation="annulation to 3na")
rr["general_scheme_ids"] = [f"{DOC}:scheme2_general"]
rr["participants"] = [participant("p_1n", mid("1n"), "reactant", ev(2)), participant("p_2a", mid("2a"), "reactant", ev(2)),
    participant("p_cat", mid("P(p-Tol)3"), "catalyst", ev(2), mol_percent=("20 mol%", q(20, "mol%")), reference=("1n", mid("1n")), method="inherited_general_condition"),
    participant("p_base", mid("K3PO4.3H2O"), "reagent", ev(2), mol_percent=("50 mol%", q(50, "mol%")), reference=("1n", mid("1n")), method="inherited_general_condition"),
    participant("p_product", mid("3na"), "product", ev(2), product_role="target product")]
set_solvent(rr, ev(2), mid("CH2Cl2"), "CH2Cl2", method="inherited_general_condition")
add_step(rr, ev(2), operation="under the optimal condition", temperature=("60 degrees C", q(60, "°C")),
         time=("24 h", q(24, "h")), method="inherited_general_condition")
add_outcome(rr, "mass", "product_amount", [mid("3na")], ev(2), "1.04 g", q(1.04, "g"), product_role="isolated target product", scope={"scale":"gram scale"})
add_outcome(rr, "yield", "yield", [mid("3na")], ev(2), "76%", q(76, "%"), product_role="target product", yield_type="reported gram-scale yield")
add_outcome(rr, "catalyst_recovery", "recovery", [mid("P(p-Tol)3")], ev(2), "nearly 50% catalyst could be recovered", q(50, "%", qualifier="approx", raw_text="nearly 50%"), product_role="recovered catalyst", yield_type="recovery")
set_classification(rr, ev(2), ANNULATION_CHANGES)
rr["review"]["unresolved"] = ["Gram-scale substrate charges, solvent volume, product-yield basis, catalyst-recovery method, dr, workup, and purification are unreported."]
reactions.append(rr)

# Scheme 3 successful products and two no-product experiments.
scheme3_rows = [
    ("3ab", "2b", 58, "10:1", True), ("3ac", "2c", 80, "5.6:1", False),
    ("3ad", "2d", 88, "6.2:1", False), ("3ae", "2e", 88, "9.5:1", False),
    ("3af", "2f", 78, "6.3:1", True), ("3ag", "2g", 69, None, False),
    ("3ah", "2h", 77, "5:1", False), ("3ai", "2i", 63, "7.7:1", True),
    ("3aj", "2j", 25, "1.7:1", False), ("3ak", "2k", 78, "7.7:1", True),
    ("3al", "2l", 86, "7.5:1", False),
]
for product, allenoate, yld, dr, nmr in scheme3_rows:
    scope_annulation(product, "1a", allenoate, 3, yld, dr, nmr,
                     scheme="scheme3_general", multi_isomer=(product == "3ag"))
for product, allenoate, substituent in (("3am", "2m", "Cy"), ("3an", "2n", "H")):
    evidence = ev(3)
    rr = reaction_base(f"scheme3_{product}_ND", evidence, categories=["substrate_scope", "control"],
        location=f"Page 3 Scheme 3 entry {product}", name="attempted sequential annulation",
        transformation="no detected annulation product")
    rr["general_scheme_ids"] = [f"{DOC}:scheme3_general"]
    rr["participants"] = [participant("p_1a", mid("1a"), "reactant", evidence,
        amount=("0.1 mmol", q(0.1,"mmol")), method="inherited_general_condition"),
        participant("p_allenoate", mid(allenoate), "reactant", evidence,
        amount=("0.12 mmol", q(0.12,"mmol")), method="inherited_general_condition"),
        participant("p_cat", mid("P(p-Tol)3"), "catalyst", evidence,
        mol_percent=("20 mol%", q(20,"mol%")), reference=("1a",mid("1a")), method="inherited_general_condition"),
        participant("p_base", mid("K3PO4.3H2O"), "reagent", evidence,
        mol_percent=("50 mol%", q(50,"mol%")), reference=("1a",mid("1a")), method="inherited_general_condition")]
    set_solvent(rr, evidence, mid("CH2Cl2"), "CH2Cl2", volume=1.0, method="inherited_general_condition")
    add_step(rr, evidence, operation="reaction carried out under optimal conditions",
             temperature=("60 degrees C",q(60,"°C")), time=("24 h",q(24,"h")), method="inherited_general_condition")
    add_outcome(rr, "no_product", "no_reaction", [], evidence, f"{product}, R = {substituent}, ND", "product not detected", scope={"expected_label":product})
    set_classification(rr, evidence, coverage="not_applicable", limitation="No corresponding product detected; no expected bonds classified.")
    rr["review"]["unresolved"] = ["ND detection method and limit are unreported; no product participant is invented."]
    reactions.append(rr)

# Scheme 4a: source conflict 71% (scheme) vs 73% (prose).
e4 = ev(3)
rr = reaction_base("scheme4a_3u_2a_3aa", e4, categories=["mechanistic_study"], location="Page 3 Scheme 4a and prose",
    name="mechanistic sulfonamide-origin experiment", transformation="annulation to 3aa")
rr["general_scheme_ids"] = [f"{DOC}:scheme1_current_generic"]
rr["participants"] = [participant("p_3u", mid("3u"), "reactant", e4, amount=("0.1 mmol",q(0.1,"mmol")), method="inherited_general_condition"),
    participant("p_2a", mid("2a"), "reactant", e4, amount=("0.12 mmol",q(0.12,"mmol")), method="inherited_general_condition"),
    participant("p_cat",mid("P(p-Tol)3"),"catalyst",e4,mol_percent=("20 mol%",q(20,"mol%")),reference=("3u",mid("3u")),method="inherited_general_condition"),
    participant("p_base",mid("K3PO4.3H2O"),"reagent",e4,mol_percent=("50 mol%",q(50,"mol%")),reference=("3u",mid("3u")),method="inherited_general_condition"),
    participant("p_product",mid("3aa"),"product",e4,product_role="target product")]
set_solvent(rr,e4,mid("CH2Cl2"),"CH2Cl2",volume=1.0,method="inherited_general_condition")
add_step(rr,e4,operation="optimal condition",temperature=("60 degrees C",q(60,"°C")),time=("24 h",q(24,"h")),method="inherited_general_condition")
alts = [f("Scheme 4a: 71%",q(71,"%"),[f"{DOC}:page3:image"]), f("prose: 73% yield",q(73,"%"),[f"{DOC}:page3:text"])]
add_outcome(rr,"yield","yield",[mid("3aa")],e4,["71% (scheme)","73% (prose)"],None,product_role="target product",yield_type="reported yield",alternatives=alts)
add_outcome(rr,"dr","selectivity",[mid("3aa")],e4,"8.3:1","8.3:1",product_role="major anti / minor syn",scope={"metric":"dr","order":["anti","syn"]})
set_classification(rr,e4,ANNULATION_CHANGES)
rr["notes"].append(f("Scheme 4a reports 71%, whereas the adjacent prose reports 73%; neither value is selected.",
                      "Unresolved source conflict: 71% vs 73%",e4))
rr["review"]["unresolved"] = ["Yield conflict 71% versus 73% is unresolved.","Yield/dr measurement methods and full structures are unreported."]
reactions.append(rr)

# Scheme 4b deuterium-labelling experiment.
rr = reaction_base("scheme4b_deuterium", e4, categories=["mechanistic_study"], location="Page 3 Scheme 4b",
    name="deuterium-labeling experiment", transformation="annulation with deuterium incorporation")
rr["general_scheme_ids"] = [f"{DOC}:scheme1_current_generic"]
rr["participants"] = [participant("p_1a",mid("1a"),"reactant",e4,amount=("0.1 mmol",q(0.1,"mmol")),method="inherited_general_condition"),
    participant("p_2b",mid("2b"),"reactant",e4,amount=("0.12 mmol",q(0.12,"mmol")),method="inherited_general_condition"),
    participant("p_D2O",mid("D2O"),"reagent",e4,equivalents=("1.5 equiv D2O",q(1.5,"equiv")),reference=("1a",mid("1a"))),
    participant("p_cat",mid("P(p-Tol)3"),"catalyst",e4,mol_percent=("20 mol%",q(20,"mol%")),reference=("1a",mid("1a")),method="inherited_general_condition"),
    participant("p_base",mid("K3PO4.3H2O"),"reagent",e4,mol_percent=("50 mol%",q(50,"mol%")),reference=("1a",mid("1a")),method="inherited_general_condition"),
    participant("p_product",mid("D-3ab"),"product",e4,product_role="deuterium-labelled target product")]
set_solvent(rr,e4,mid("CH2Cl2"),"CH2Cl2",volume=1.0,method="inherited_general_condition")
add_step(rr,e4,operation="optimal condition with 1.5 equiv D2O",temperature=("60 degrees C",q(60,"°C")),time=("24 h",q(24,"h")),method="inherited_general_condition")
add_outcome(rr,"yield","yield",[mid("D-3ab")],e4,"D-3ab, 47%",q(47,"%"),product_role="deuterium-labelled product",yield_type="reported yield")
inc = {"displayed_position_1":"35% D","displayed_position_2":"3% D","displayed_position_3":"25% D","displayed_position_4":"2% D"}
add_outcome(rr,"D_incorporation","other",[mid("D-3ab")],e4,
            "Scheme labels 35% D, 3% D, 25% D, and 2% D at four displayed positions",inc,
            product_role="deuterium-labelled product",basis="1H NMR analysis",
            scope={"metric":"deuterium incorporation","atom_mapping":"unresolved","author_summary":"incorporation at C2 and nearly none at C1"})
set_classification(rr,e4,ANNULATION_CHANGES)
rr["review"]["unresolved"] = ["Exact atom mapping of the four displayed D percentages is not safely established from the main-text diagram.","Yield basis and full isotopologue distribution are unreported."]
reactions.append(rr)

# Scheme 4c imine experiment.
rr = reaction_base("scheme4c_imine11", e4, categories=["mechanistic_study"], location="Page 3 Scheme 4c",
    name="imine 11/allenoate mechanistic experiment", transformation="benzannulation to 12")
rr["participants"] = [participant("p_11",mid("11"),"reactant",e4),
    participant("p_2a",mid("2a"),"reactant",e4,equivalents=("1.5 equiv",q(1.5,"equiv")),reference=("11",mid("11"))),
    participant("p_cat",mid("PCy3"),"catalyst",e4,mol_percent=("40 mol%",q(40,"mol%")),reference=("11",mid("11"))),
    participant("p_base",mid("K2CO3"),"reagent",e4,mol_percent=("40 mol%",q(40,"mol%")),reference=("11",mid("11"))),
    participant("p_product",mid("12"),"product",e4,product_role="benzannulation product")]
set_solvent(rr,e4,mid("THF"),"THF")
add_step(rr,e4,temperature=("60 degrees C",q(60,"°C")),time=("1 d",q(1,"d")))
add_outcome(rr,"yield","yield",[mid("12")],e4,"41%",q(41,"%"),product_role="benzannulation product",yield_type="reported yield")
set_classification(rr,e4,[("C","C",0,1,"formed"),("C","C",0,1,"formed")],
                   limitation="Only the two [3+3] C-C ring-forming connections are retained; aromatization order changes and hydrogen balance are unresolved.")
rr["review"]["unresolved"] = ["Reactant charge, solvent volume, yield basis, workup, purification, complete structures, and aromatization atom mapping are unreported."]
reactions.append(rr)

# Figure 2 proposed mechanism.
rr = reaction_base("figure2_proposed_mechanism", ev(3), "proposed_mechanism", location="Page 3 Figure 2",
    name="possible mechanism", transformation="proposed [3+3] annulation followed by aza-6pi-electrocyclization")
rr["participants"] = [participant("p_2a",mid("2a"),"reactant",ev(3)),
    participant("p_1",g1,"reactant",ev(3)),participant("p_3aa",mid("3aa"),"product",ev(3),product_role="proposed pathway product"),
    participant("p_12",mid("12"),"product",ev(3),product_role="alternative proposed aromatization product")]
set_classification(rr,ev(3),coverage="unresolved",limitation="Elementary arrows are author-proposed and are not treated as independently demonstrated net experiments.")
rr["notes"].append(f("Intermediate A through K and proton shifts are author-proposed; Scheme 4 experiments are separately recorded.",
                     "author-proposed mechanism",ev(3)))
rr["review"]["unresolved"] = ["Intermediates A-K and elementary atom mappings are proposed rather than independently observed."]
reactions.append(rr)


def derivatization(suffix, start, product, reagent, solvent, equiv, yield_value, label,
                   changes, *, dry=False):
    evidence = ev(4)
    rr = reaction_base(suffix,evidence,categories=["derivatization"],location=f"Page 4 Scheme 5 {label}",
        name="transformation of tetrahydroisoquinoline",transformation=label)
    rr["participants"] = [participant("p_start",mid(start),"reactant",evidence),
        participant("p_reagent",mid(reagent),"reagent",evidence,equivalents=(f"{equiv} equiv",q(equiv,"equiv")),reference=(start,mid(start))),
        participant("p_product",mid(product),"product",evidence,product_role="derivatization product")]
    set_solvent(rr,evidence,mid(solvent),solvent,dry=dry)
    add_step(rr,evidence,operation=f"Scheme 5 condition {label}; temperature and time not reported")
    add_outcome(rr,"yield","yield",[mid(product)],evidence,f"{yield_value}%",q(yield_value,"%"),product_role="derivatization product",yield_type="reported yield")
    set_classification(rr,evidence,changes)
    rr["review"]["unresolved"] = ["Temperature, time, reagent/substrate charges, solvent volume, yield basis, workup, purification, and complete atom mapping are unreported."]
    reactions.append(rr)


derivatization("scheme5a_3na_to_5","3na","5","LiAlH4","THF",2,80,"a: ester reduction",
    [("C","O",2,1,"order_decreased"),("C","O",1,0,"broken"),("C","H",0,1,"hydrogen_gain",2),("O","H",0,1,"hydrogen_gain",1)],dry=True)
derivatization("scheme5b_3na_to_6na","3na","6na","SmI2","THF",6,85,"b: detosylation/spontaneous aromatization",
    [("C","N",1,2,"order_increased"),("C","C",1,2,"order_increased")],dry=True)
derivatization("scheme5c_6na_to_7na","6na","7na","DDQ","toluene",1.25,89,"c: DDQ oxidation",
    [("C","C",1,2,"order_increased")])
derivatization("scheme5d_3na_to_7na","3na","7na","DDQ","toluene",3,65,"d: direct DDQ oxidation",
    [("C","N",1,2,"order_increased"),("C","C",1,2,"order_increased"),("C","C",1,2,"order_increased")])
derivatization("scheme5e_i_6na_to_8","6na","8","KOH","EtOH",10,84,"e(i): ester hydrolysis",
    [("C","O",1,0,"broken"),("O","H",0,1,"hydrogen_gain",1)])
derivatization("scheme5e_ii_8_to_9","8","9","TfOH","CH2Cl2",5,93,"e(ii): Friedel-Crafts acylation",
    [("C","C",0,1,"formed")],dry=True)
derivatization("scheme5f_6ma_to_10","6ma","10","BBr3","CH2Cl2",2.2,73,"f: lactonization/demethylation",
    [("C","O",0,1,"formed"),("O","H",0,1,"hydrogen_gain",1)],dry=True)


# Source context that must not become synthetic experiments.
bundle["source_context"] = [
    {"context_id":f"{DOC}:figure1_bioactive","kind":"non_synthetic_figure","description":"Representative bioactive (hydro)isoquinoline structures; no current-paper synthesis arrow.","evidence_ids":ev(2)},
    {"context_id":f"{DOC}:xray_3fa","kind":"characterization","description":"Relative configuration of major anti-3fa confirmed by X-ray; CCDC 2012531.","evidence_ids":ev(2)+ev(4)},
    {"context_id":f"{DOC}:asymmetric_attempts","kind":"qualitative_mechanistic_context","description":"Main text says chiral-phosphine trials resulted in racemization; SI Scheme S3 is not supplied, so individual experiments are not enumerated.","evidence_ids":ev(3)},
    {"context_id":f"{DOC}:asymmetric_azatrienes","kind":"qualitative_scope_context","description":"Asymmetrically substituted azatrienes gave mixtures of isomers; SI Scheme S1 is not supplied, so individual experiments are not enumerated.","evidence_ids":ev(3)},
    {"context_id":f"{DOC}:si_notice","kind":"availability_notice","description":"The article advertises SI, but no SI file was supplied or explicitly associated with this run.","evidence_ids":ev(4)},
    {"context_id":f"{DOC}:references_pages5_6","kind":"bibliography","description":"Pages 5-6 contain bibliography/CCDC note only; cited reactions are not imported as current-paper experiments.","evidence_ids":ev(5)+ev(6)},
]

bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
st.write(RUN / "submitted.json", bundle)
st.write(RUN / "raw_model/stereo_checks.json", stereo_checks)
st.write(RUN / "raw_model/source_mentions.json", st.source_mentions(data,
    ["3aa","3ab","3na","3ra","3ya","3za","3ag","3am","3an","6na","7na","12"], context=220))

page_records = {page: [] for page in range(1,7)}
for reaction in reactions:
    pages = set()
    for eid in reaction["evidence_ids"]:
        if ":page" in eid:
            pages.add(int(eid.split(":page",1)[1].split(":",1)[0]))
    for page in pages:
        page_records[page].append(reaction["reaction_id"])

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page":1,"state":"attempted","regions":[
            {"region":"graphical abstract and Scheme 1d","finding":"One generic current-work sequential annulation; duplicate depictions merged.","evidence_ids":ev(1)},
            {"region":"Scheme 1a-c","finding":"Three prior-work generic examples classified as literature_example, not current-paper experiments.","evidence_ids":ev(1)},
            {"region":"remaining title/abstract/body","finding":"Article metadata and motivation; no additional specific current-paper experiment.","evidence_ids":ev(1)}],
         "record_ids":page_records[1],"unresolved":["Generic R/Ar structures and exact atom mapping remain unresolved."]},
        {"page":2,"state":"attempted","regions":[
            {"region":"optimization prose","finding":"Initial 1a+2b/PCy3/K2CO3/chloroform experiment gives 3ab in 24%, dr 10:1; unavailable Table S1 holds missing details.","evidence_ids":ev(2)},
            {"region":"Scheme 2 and footnotes","finding":"Twenty displayed scope outcomes (3aa-3ta, including decomposed 3ra) plus one generic scheme; all rows enumerated.","evidence_ids":ev(2)},
            {"region":"eq 1 and eq 2","finding":"3ya is obtained in 15%, dr 4.7:1; 3za is ND and 1z recovery is 67%.","evidence_ids":ev(2)},
            {"region":"gram-scale prose and Scheme 2 parenthesis","finding":"Distinct gram-scale 3na outcome: 1.04 g, 76%, with nearly 50% catalyst recovery.","evidence_ids":ev(2)},
            {"region":"Figure 1/X-ray discussion","finding":"Non-synthetic bioactive figure and 3fa characterization context retained outside reaction count.","evidence_ids":ev(2)}],
         "record_ids":page_records[2],"unresolved":["No SI/Table S1 supplied; measurement bases, full structures, workups, and absolute configurations are incomplete."]},
        {"page":3,"state":"attempted","regions":[
            {"region":"Scheme 3 and footnotes","finding":"Eleven successful scope products 3ab-3al plus two ND experiments 3am/3an; 3ag is an at-least-three-pair isomer mixture.","evidence_ids":ev(3)},
            {"region":"Scheme 4a","finding":"Mechanistic sulfonamide-origin experiment; scheme yield 71% conflicts with prose 73%, dr 8.3:1.","evidence_ids":ev(3)},
            {"region":"Scheme 4b","finding":"D2O labelling gives D-3ab in 47% with four displayed D incorporations; atom mapping remains unresolved.","evidence_ids":ev(3)},
            {"region":"Scheme 4c","finding":"Imine 11/allenoate experiment gives benzannulation product 12 in 41%.","evidence_ids":ev(3)},
            {"region":"Figure 2","finding":"Author-proposed mechanism recorded separately from mechanistic experiments.","evidence_ids":ev(3)}],
         "record_ids":page_records[3],"unresolved":["SI-only asymmetric/scope trials cannot be enumerated; exact D-position mapping and complete atom mapping remain unresolved."]},
        {"page":4,"state":"attempted","regions":[
            {"region":"Scheme 5 and footnote","finding":"Seven specific derivatization stages recorded: a, b, c, d, e(i), e(ii), and f, preserving intermediate 8 and both stage yields.","evidence_ids":ev(4)},
            {"region":"conclusion/Associated Content","finding":"No additional current-paper reaction; SI advertised but not supplied.","evidence_ids":ev(4)},
            {"region":"CCDC/author information/references start","finding":"Characterization/accession and bibliography context only.","evidence_ids":ev(4)}],
         "record_ids":page_records[4],"unresolved":["Scheme 5 omits temperature, time, scale, solvent volume, measurement basis, workup, and purification."]},
        {"page":5,"state":"attempted","regions":[
            {"region":"full page","finding":"Bibliography only; cited literature not imported as experiments newly reported in this paper.","evidence_ids":ev(5)}],
         "record_ids":[],"unresolved":[]},
        {"page":6,"state":"attempted","regions":[
            {"region":"full page","finding":"Bibliography and CCDC note continuation only; no current-paper synthesis experiment.","evidence_ids":ev(6)}],
         "record_ids":[],"unresolved":[]},
    ]
}
st.write(RUN / "coverage/pages_0001_0006.json", coverage)

inventory = {
    "document_id":DOC,
    "pages_attempted":[1,2,3,4,5,6],
    "reaction_records":len(reactions),
    "specific_experiments":sum(r["record_type"]=="specific_experiment" for r in reactions),
    "generic_schemes":sum(r["record_type"]=="generic_scheme" for r in reactions),
    "literature_examples":sum(r["record_type"]=="literature_example" for r in reactions),
    "proposed_mechanisms":sum(r["record_type"]=="proposed_mechanism" for r in reactions),
    "molecules":len(molecules),
    "page_record_ids":page_records,
    "limitations":["No SI was supplied or explicitly associated.","All machine records remain needs_review.","Page coverage is attempted, not semantic accuracy verification."],
}
st.write(RUN / "raw_model/pages_0001_0006_extraction.json", inventory)

tasks = st.read(RUN / "tasks.json")
for task in tasks:
    task["state"] = "attempted"
    task["reason"] = "Original page image and extracted text reviewed; see coverage/pages_0001_0006.json. Attempted does not certify completeness."
st.write(RUN / "tasks.json", tasks)
st.write(RUN / "raw_model/audit.json", reaction_tools.dataset_audit(st.read(DATASET_PATH)))

print(json.dumps(inventory, ensure_ascii=False, indent=2))
