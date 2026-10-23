#!/usr/bin/env python3
"""Build the evidence-linked v3.1 submission for PDF pages 1--6.

This is a source-specific extraction script, not a chemistry inference engine.
Uncertain structures and atom mappings remain explicitly unresolved/partial.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st  # noqa: E402


RUN = ROOT / "runs/current"
DATASET = RUN / "dataset.json"
SUBMITTED = RUN / "submitted.json"
DOC = "doc_b7f7182105e32edc"
TITLE = (
    "Regio- and Stereoselective Synthesis of Rotationally Hindered "
    "C12-Naphthylated Tribenzo[a,c,j]xanthenes through Catalytic "
    "Tricyclization of Yne-Allenones"
)


def ev(page: int, kind: str = "text") -> str:
    return f"{DOC}:page{page}:{kind}"


def fact(raw=None, normalized=None, evidence=(), method="original_extraction",
         status=None, derivation=None):
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   derivation=derivation)


def q(number, unit, **extra):
    return {"number": number, "unit": unit, **extra}


def set_source(f, raw, normalized, evidence, method="original_extraction",
               derivation=None, status=None):
    f.clear()
    f.update(fact(raw, normalized, evidence, method=method,
                  derivation=derivation, status=status))


def set_na(f, evidence=()):
    f.clear()
    f.update(fact(evidence=evidence, method="not_applicable",
                  status="not_applicable"))


def stereo_candidate(m, evidence, smiles, allene=False):
    diag = st.stereo_inspect(smiles)
    deriv = "RDKit graph diagnostics on the source-derived candidate SMILES; this does not verify source identity."
    for key, value in (
        ("atom_centers", diag["atom_centers"]),
        ("double_bonds", diag["double_bonds"]),
        ("axial", []),
        ("mixture", None),
        ("assignment", None),
    ):
        if value is None:
            continue
        m["stereochemistry"][key] = fact(
            "candidate graph diagnostics", value, evidence,
            method="structure_derivation", derivation=deriv)
    m["stereochemistry"]["chirality_status"] = fact(
        "source drawing/name contains no stereogenic center assignment",
        "achiral", evidence, method="structure_derivation",
        derivation=(
            "The yne-allenone has a terminal CH2 allene end and no other stereogenic element in the candidate graph."
            if allene else "No stereogenic element is present in the source-derived candidate graph."
        ))
    m["stereochemistry"]["geometry_status"] = fact(
        "no applicable E/Z assignment", "not_applicable", evidence,
        method="structure_derivation",
        derivation=(
            "The terminal allene is not assigned as an alkene E/Z unit."
            if allene else "The candidate graph contains no applicable alkene E/Z element."
        ))
    m["stereochemistry"]["representation_status"] = fact(
        "candidate SMILES stereo representation",
        "partial" if allene else "not_applicable", evidence,
        method="structure_derivation",
        derivation=(
            "Toolkit diagnostics report an unspecified cumulene/double-bond feature; the source shows a terminal allene and no axial assignment."
            if allene else "No stereochemical encoding is applicable to this candidate graph."
        ))
    m["stereochemistry"]["isomeric_smiles"] = fact(
        "candidate graph", diag["canonical_isomeric_smiles"], evidence,
        method="structure_derivation", derivation=deriv)
    m["stereochemistry"]["identity_key"] = fact(
        "candidate graph", diag["stereo_graph_key"], evidence,
        method="structure_derivation", derivation=deriv)
    return diag


def molecule(mid, label, evidence, name=None, smiles=None, allene=False,
             identity_note=None):
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = fact(label, label, evidence)
    m["names"]["original"] = fact(name or label, name or label, evidence)
    m["names"]["english"] = fact(name or label, name or label, evidence)
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = fact(
            name or label, smiles, evidence, method="structure_derivation",
            derivation="Candidate SMILES derived from the printed structure and/or explicit chemical name; parsing is not source verification.")
        stereo_candidate(m, evidence, smiles, allene=allene)
    else:
        m["representation_type"] = "unresolved"
        m["stereochemistry"]["chirality_status"] = fact(
            "stereochemistry not fully identifiable from retained label",
            "unresolved", evidence, method="inference",
            derivation="No complete candidate graph was asserted.")
        m["stereochemistry"]["geometry_status"] = fact(
            "geometry not fully identifiable from retained label",
            "unresolved", evidence, method="inference",
            derivation="No complete candidate graph was asserted.")
        m["stereochemistry"]["representation_status"] = fact(
            "structure not encoded", "unresolved", evidence,
            method="inference", derivation="No complete candidate graph was asserted.")
    if identity_note:
        m["identity_note"] = fact(identity_note, identity_note, evidence)
    return m


def product_mixture(mid, label, evidence):
    m = molecule(mid, label, evidence,
                 name=f"reported racemic Z/E product mixture {label}",
                 identity_note=(
                     "Complex tribenzo[a,c,j]xanthene structure retained from the drawing without a concrete SMILES; "
                     "the source reports a rotationally hindered chiral axis and Z/E mixtures."
                 ))
    m["stereochemistry"]["chirality_status"] = fact(
        "rac-3 in Schemes 1--2", "racemic", evidence)
    m["stereochemistry"]["geometry_status"] = fact(
        "Z/E ratios reported for scope products", "mixture", evidence)
    m["stereochemistry"]["representation_status"] = fact(
        "axial and Z/E stereochemistry not encoded as SMILES", "not_encoded", evidence)
    m["stereochemistry"]["axial"] = fact(
        "rotationally hindered chiral axis; no R_a/S_a assignment reported",
        [{"axis": "hindered aryl axis", "configuration": "unassigned", "sample": "racemic"}], evidence)
    m["stereochemistry"]["double_bonds"] = fact(
        "(Z)-products are major; E is minor", [{"configuration": "Z/E mixture", "major": "Z"}], evidence)
    m["stereochemistry"]["mixture"] = fact(
        "rac-3; experiment-specific Z/E ratios", {"axial": "racemic", "geometry": "Z/E mixture"}, evidence)
    m["stereochemistry"]["assignment"] = fact(
        "3b confirmed by single-crystal X-ray; derivatives assigned by analogy",
        {"3b": "single-crystal X-ray, CCDC 2054923", "others": "assigned by analogy"}, evidence)
    return m


def generic_molecule(mid, label, evidence, fixed, variables, reason,
                     product=False):
    m = molecule(mid, label, evidence, name=label, identity_note=reason)
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": fact(fixed, fixed, evidence),
        "variables": [
            {
                "label": v,
                "attachment_points": fact(f"printed attachment point for {v}", f"printed attachment point for {v}", evidence),
                "definition": fact(f"Scheme variable {v}", f"Scheme variable {v}", evidence),
                "allowed_values": fact(evidence=evidence),
                "constraints": fact(evidence=evidence),
            }
            for v in variables
        ],
        "instance_mapping": [],
        "reason": fact(reason, reason, evidence),
    }
    if product:
        m["stereochemistry"]["chirality_status"] = fact("rac-3", "racemic", evidence)
        m["stereochemistry"]["geometry_status"] = fact("generic product contains Z/E element", "unspecified", evidence)
        m["stereochemistry"]["representation_status"] = fact("R-group and axial stereo not encoded", "not_encoded", evidence)
        m["stereochemistry"]["axial"] = fact("rotationally hindered axis in rac-3", [{"configuration": "unassigned", "sample": "racemic"}], evidence)
    return m


def participant(pid, mid, role, evidence, stage="reaction", product_role=None,
                mass=None, amount=None, volume=None, equivalents=None,
                mol_percent=None, reference=None, inherited=False):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    method = "inherited_general_condition" if inherited else "original_extraction"
    deriv = "Inherited from Table 1 entry 15 / the stated standard operation." if inherited else None
    if product_role is not None:
        p["product_role"] = fact(product_role, product_role, evidence)
    if mass is not None:
        p["charges"]["mass"] = fact(f"{mass['number']} {mass['unit']}", mass, evidence, method=method, derivation=deriv)
    if amount is not None:
        p["charges"]["amount"] = fact(f"{amount['number']} {amount['unit']}", amount, evidence, method=method, derivation=deriv)
    if volume is not None:
        p["charges"]["volume"] = fact(f"{volume['number']} {volume['unit']}", volume, evidence, method=method, derivation=deriv)
    if equivalents is not None:
        p["charges"]["equivalents"] = fact(f"{equivalents['number']} equiv", equivalents, evidence, method=method, derivation=deriv)
    if mol_percent is not None:
        p["charges"]["mol_percent"] = fact(f"{mol_percent['number']} mol%", mol_percent, evidence, method=method, derivation=deriv)
    if reference is not None:
        p["equivalence_reference"] = fact(reference, reference, evidence, method=method, derivation=deriv)
    return p


def step(step_id, evidence, pids, operation=None, temperature=None, time=None,
         atmosphere=None, stage="reaction", inherited=False):
    s = st.new_step(step_id, evidence, stage=stage)
    s["participant_ids"] = list(pids)
    method = "inherited_general_condition" if inherited else "original_extraction"
    deriv = "Inherited from Table 1 entry 15 / the stated standard operation." if inherited else None
    if operation is not None:
        s["operation"] = fact(operation, operation, evidence, method=method, derivation=deriv)
    if temperature is not None:
        raw = temperature if isinstance(temperature, str) else f"{temperature['number']} {temperature['unit']}"
        s["temperature"] = fact(raw, temperature, evidence, method=method, derivation=deriv)
    if time is not None:
        raw = time if isinstance(time, str) else f"{time['number']} {time['unit']}"
        s["time"] = fact(raw, time, evidence, method=method, derivation=deriv)
    if atmosphere is not None:
        s["atmosphere"] = fact(atmosphere, atmosphere, evidence, method=method, derivation=deriv)
    return s


def outcome(outcome_id, kind, targets, evidence, value=None, raw=None,
            product_role=None, yield_type=None, measurement=None, scope=None,
            inherited=False):
    o = st.new_outcome(outcome_id, kind, targets, evidence)
    method = "inherited_general_condition" if inherited else "original_extraction"
    deriv = "Inherited from the cited general condition." if inherited else None
    if value is not None:
        o["value"] = fact(raw if raw is not None else value, value, evidence, method=method, derivation=deriv)
    if product_role is not None:
        o["product_role"] = fact(product_role, product_role, evidence)
    if yield_type is not None:
        o["yield_type"] = fact(yield_type, yield_type, evidence)
    if measurement is not None:
        o["measurement_basis"] = fact(measurement, measurement, evidence)
    if scope is not None:
        o["scope"] = fact(scope, scope, evidence)
    return o


def reactive_classification(evidence, c_o=True, c_c=True, limitation=None):
    changes = []
    labels = []
    if c_c:
        change = {"change_id": "class_C-C", "atoms": [{"element": "C"}, {"element": "C"}],
                  "before_order": 0, "after_order": 1, "change_type": "formed"}
        changes.append(fact(None, change, evidence, method="inference",
                            derivation="Class-level bond formation inferred from the source cycloaddition/annulation description; atom mapping is unresolved."))
        labels.append({"category_number": 2, "code": "C-C", "supporting_change_ids": ["class_C-C"]})
    if c_o:
        change = {"change_id": "class_C-O", "atoms": [{"element": "C"}, {"element": "O"}],
                  "before_order": 0, "after_order": 1, "change_type": "formed"}
        changes.append(fact(None, change, evidence, method="inference",
                            derivation="Class-level C-O formation inferred from the source oxa-cyclization description; atom mapping is unresolved."))
        labels.append({"category_number": 3, "code": "C-O", "supporting_change_ids": ["class_C-O"]})
    return {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": fact(evidence=evidence),
        "bond_changes": changes,
        "labels": fact(None, labels, evidence, method="inference",
                       derivation="Taxonomy labels are inferred only at bond-class level from the reported transformation."),
        "coverage": fact(None, "partial", evidence, method="inference",
                         derivation="Exact product atom mapping and the full set/count of new bonds are unresolved."),
        "limitations": [
            limitation or "Class-level C-C/C-O labels are retained; exact atom identities, bond counts, hydrogen changes, and complete atom mapping are unresolved."
        ],
        "model_votes": [],
    }


def no_reaction_classification(evidence):
    return {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": fact(evidence=evidence),
        "bond_changes": [],
        "labels": fact(None, [], evidence, method="inference",
                       derivation="No observed target transformation; no bond-formation label asserted."),
        "coverage": fact(None, "not_applicable", evidence, method="inference",
                         derivation="No observed target transformation."),
        "limitations": ["No bond changes are inferred from an intended but unobserved reaction."],
        "model_votes": [],
    }


def reaction(rid, evidence, record_type, categories, location, name,
             transformation, classification, unresolved):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = fact(categories, categories, evidence)
    r["location"] = fact(location, location, evidence)
    r["reaction_name"] = fact(name, name, evidence)
    r["transformation_type"] = fact(
        transformation, transformation, evidence, method="inference",
        derivation="Normalized from the authors' reaction description.")
    r["classification"] = classification
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": list(unresolved)}
    r["completeness"] = "partial"
    return r


def set_solvent(r, name, volume, evidence, inherited=False):
    method = "inherited_general_condition" if inherited else "original_extraction"
    deriv = "Inherited from Table 1 entry 15 / the stated standard operation." if inherited else None
    r["solvent_conditions"]["composition"] = fact(name, [name], evidence, method=method, derivation=deriv)
    if volume is not None:
        r["solvent_conditions"]["total_volume"] = fact(
            f"{volume['number']} {volume['unit']}", volume, evidence,
            method=method, derivation=deriv)


def add_selectivity(r, product_mid, ratio, evidence, scope_text="Z/E ratio"):
    r["outcomes"].append(outcome(
        f"{r['reaction_id']}:selectivity", "selectivity", [product_mid], evidence,
        value={"metric": "Z/E", "reported_ratio": ratio, "order": ["Z", "E"]},
        raw=f"Z/E = {ratio}", product_role="reported product mixture",
        scope=scope_text))


def main():
    data = st.read(DATASET)
    if data["experiments"] or data["candidates"] or "standard" in data:
        raise RuntimeError("Expected the mechanically prepared empty dataset; refusing to reuse prior records.")
    data["documents"][0]["title"] = TITLE
    st.write(DATASET, data)

    bundle = st.read(SUBMITTED)
    bundle["documents"] = data["documents"]
    mols = []
    stereo_checks = {}

    # Concrete structures directly supported by names/drawings.
    m1a = f"{DOC}:mol:1a"
    m2a = f"{DOC}:mol:2a"
    one_a = molecule(m1a, "1a", [ev(2, "image"), ev(5)],
                     "1-(2-(phenylethynyl)phenyl)buta-2,3-dien-1-one",
                     "C=C=CC(=O)c1ccccc1C#Cc1ccccc1", allene=True)
    two_a = molecule(m2a, "2a", [ev(2, "image"), ev(5)],
                     "1-(phenylethynyl)naphthalen-2-ol",
                     "Oc1ccc2ccccc2c1C#Cc1ccccc1")
    mols += [one_a, two_a]
    stereo_checks["1a"] = st.stereo_inspect(st.value(one_a["smiles"]))
    stereo_checks["2a"] = st.stereo_inspect(st.value(two_a["smiles"]))

    # Generic drawings (R groups are templates, not wildcard SMILES).
    tpl1 = f"{DOC}:mol:template:yne_allenone_1"
    tpl2 = f"{DOC}:mol:template:alkynylnaphthol_2"
    tpl3 = f"{DOC}:mol:template:rac_product_3"
    mols += [
        generic_molecule(tpl1, "generic yne-allenone 1", [ev(2, "image"), ev(3, "image")],
                         "printed yne-allenone scaffold with terminal allene and arylalkyne",
                         ["R1", "R2"], "Scheme uses variable substituents R1/R2; no wildcard concrete SMILES asserted."),
        generic_molecule(tpl2, "generic 1-alkynylnaphthalen-2-ol 2", [ev(2, "image"), ev(3, "image")],
                         "printed 1-alkynylnaphthalen-2-ol scaffold",
                         ["R3", "Ar"], "Scheme uses variable substituents R3/Ar; no wildcard concrete SMILES asserted."),
        generic_molecule(tpl3, "generic rac-3", [ev(2, "image"), ev(3, "image")],
                         "printed tribenzo[a,c,j]xanthene product scaffold",
                         ["R1", "R2", "R3", "Ar"], "Complex R-group product and axial stereochemistry are retained as a template.", product=True),
    ]

    # Yne-allenone scope labels 1b--1t.
    one_descriptions = {
        "1b": "yne-allenone with R1 = 4-fluorophenyl", "1c": "yne-allenone with R1 = 4-chlorophenyl",
        "1d": "yne-allenone with R1 = 3-chlorophenyl", "1e": "yne-allenone with R1 = 2-chlorophenyl",
        "1f": "yne-allenone with R1 = 4-bromophenyl", "1g": "yne-allenone with R1 = p-tolyl",
        "1h": "yne-allenone with R1 = m-tolyl", "1i": "yne-allenone with R1 = 4-tert-butylphenyl",
        "1j": "yne-allenone with R1 = 1-naphthalenyl", "1k": "yne-allenone with R1 = 2-thienyl",
        "1l": "yne-allenone with R1 = cyclohex-1-en-1-yl",
        "1m": "fluoro-substituted internal-arene yne-allenone 1m",
        "1n": "chloro-substituted internal-arene yne-allenone 1n",
        "1o": "chloro-substituted internal-arene yne-allenone 1o",
        "1p": "methyl-substituted internal-arene yne-allenone 1p",
        "1q": "methyl-substituted internal-arene yne-allenone 1q",
        "1r": "methyl-substituted internal-arene yne-allenone 1r",
        "1s": "chloro-substituted internal-arene yne-allenone 1s",
        "1t": "methoxy-substituted internal-arene yne-allenone 1t",
    }
    mids_1 = {"1a": m1a}
    for label, desc in one_descriptions.items():
        mid = f"{DOC}:mol:{label}"
        mids_1[label] = mid
        candidate = "C=C=CC(=O)c1ccccc1C#Cc1ccc(Cl)cc1" if label == "1c" else None
        mols.append(molecule(mid, label, [ev(3, "image"), ev(3)], desc,
                             candidate, allene=bool(candidate),
                             identity_note=(
                                 "Candidate SMILES encoded because both the drawing and prose specify the para-chlorophenyl substituent."
                                 if candidate else
                                 "Full candidate SMILES withheld because the small scope drawing/positional assignment was not encoded with sufficient confidence."
                             )))

    two_descriptions = {
        "2b": "1-alkynylnaphthalen-2-ol with R3 = 4-fluorophenyl",
        "2c": "1-alkynylnaphthalen-2-ol with R3 = 4-chlorophenyl",
        "2d": "1-alkynylnaphthalen-2-ol with R3 = 3-chlorophenyl",
        "2e": "1-alkynylnaphthalen-2-ol with R3 = 4-bromophenyl",
        "2f": "1-alkynylnaphthalen-2-ol with R3 = 4-tert-butylphenyl",
        "2g": "1-alkynylnaphthalen-2-ol with R3 = 4-methoxyphenyl (PMP)",
        "2h": "1-alkynylnaphthalen-2-ol with R3 = 4-biphenyl",
        "2-cyclopropyl": "cyclopropyl-substituted 1-alkynylnaphthalen-2-ol",
        "2-tert-butyl": "tert-butyl-substituted 1-alkynylnaphthalen-2-ol",
    }
    mids_2 = {"2a": m2a}
    for label, desc in two_descriptions.items():
        safe = label.replace("-", "_")
        mid = f"{DOC}:mol:{safe}"
        mids_2[label] = mid
        candidate = "Oc1ccc2ccccc2c1C#Cc1ccc(C(C)(C)C)cc1" if label == "2f" else None
        mols.append(molecule(mid, label, [ev(3, "image"), ev(3)], desc,
                             candidate,
                             identity_note=(
                                 "Candidate SMILES encoded because both the drawing and prose specify the para-tert-butylphenyl substituent."
                                 if candidate else
                                 "Substituent identity is retained; a full candidate SMILES was not asserted from the small scope drawing."
                             )))

    # Products actually reported in nonzero yield.
    product_labels = [
        "3a", "3b", "3c", "3d", "3e", "3f", "3g", "3h", "3i", "3j", "3k",
        "3m", "3n", "3o", "3p", "3q", "3r", "3s", "3t",
        "3u", "3v", "3w", "3x", "3y", "3z", "3aa", "3dd",
    ]
    mids_3 = {}
    for label in product_labels:
        evidence = [ev(3, "image"), ev(3)]
        if label == "3a": evidence += [ev(2, "image"), ev(4, "image"), ev(5)]
        if label == "3b": evidence += [ev(6)]
        mid = f"{DOC}:mol:{label}"
        mids_3[label] = mid
        mols.append(product_mixture(mid, label, evidence))

    # Reagents, solvents, controls, purification agents, and prior-work entities.
    reagent_specs = {
        "BF3Et2O": ("BF3·Et2O", None), "AgOTf": ("silver triflate", None),
        "ScOTf3": ("scandium triflate", None), "SnOTf2": ("tin(II) triflate", None),
        "FeCl3": ("iron(III) chloride", None), "YOTf3": ("yttrium triflate", None),
        "DCE": ("1,2-dichloroethane", "ClCCCl"), "DCM": ("dichloromethane", "ClCCl"),
        "MeCN": ("acetonitrile", "CC#N"), "dioxane": ("1,4-dioxane", "O1CCOCC1"),
        "EtOH": ("ethanol", "CCO"), "petroleum_ether": ("petroleum ether", None),
        "ethyl_acetate": ("ethyl acetate", "CCOC(C)=O"), "silica_gel": ("silica gel", None),
        "diphenylacetylene": ("1,2-diphenylethyne", "c1ccccc1C#Cc1ccccc1"),
        "control4": ("2-(phenylethynyl)phenol (4)", "Oc1ccccc1C#Cc1ccccc1"),
        "control5": ("1-(phenylethynyl)-2-methoxynaphthalene (5)", "COc1ccc2ccccc2c1C#Cc1ccccc1"),
    }
    rmids = {}
    for key, (name, smiles) in reagent_specs.items():
        mid = f"{DOC}:mol:{key}"
        rmids[key] = mid
        evidence = [ev(2, "image"), ev(2)]
        if key in {"petroleum_ether", "ethyl_acetate", "silica_gel"}: evidence = [ev(5)]
        if key in {"control4", "control5"}: evidence = [ev(4, "image"), ev(4)]
        mol = molecule(mid, key, evidence, name, smiles,
                       identity_note=("Named mixture; composition not expanded." if key == "petroleum_ether" else None))
        mols.append(mol)
        if smiles:
            stereo_checks[key] = st.stereo_inspect(smiles)

    # Prior-work generic partner and product.
    lit_alkyne = f"{DOC}:mol:literature:alkyne"
    lit_product = f"{DOC}:mol:literature:phenanthrenol"
    mols += [
        generic_molecule(lit_alkyne, "generic prior-work alkyne/alkene partner", [ev(1), ev(2, "image")],
                         "printed generic unsaturated partner", ["R2"],
                         "Prior-work Scheme 1a is generic and cited to reference 15."),
        generic_molecule(lit_product, "generic functionalized phenanthren-9-ol", [ev(1), ev(2, "image")],
                         "printed generic phenanthren-9-ol product", ["R1", "R2"],
                         "Prior-work product is generic and was not re-encoded as a concrete molecule."),
    ]

    bundle["molecules"] = mols
    reactions = []

    # Scheme 1a: successful cited prior work and the explicitly unsuccessful diphenylacetylene case.
    lit_success = reaction(
        "scheme1a_prior_success", [ev(1), ev(2, "image")], "literature_example", [],
        "Scheme 1a, upper arrow (previous work)", "double [2+2] cycloaddition relay",
        "prior-work cycloaddition/annulation", reactive_classification([ev(1), ev(2, "image")], c_o=False),
        ["Cited prior work only; substrate scope, quantities, time, temperature, yield, and exact atom mapping are not supplied here."])
    lit_success["participants"] = [
        participant("lit_yne_allenone", tpl1, "reactant", lit_success["evidence_ids"]),
        participant("lit_partner", lit_alkyne, "reactant", lit_success["evidence_ids"]),
        participant("lit_bf3", rmids["BF3Et2O"], "catalyst", lit_success["evidence_ids"]),
        participant("lit_product", lit_product, "product", lit_success["evidence_ids"], product_role="reported prior-work product"),
    ]
    reactions.append(lit_success)

    lit_fail = reaction(
        "scheme1a_prior_diphenylacetylene_failure", [ev(1), ev(2, "image")], "literature_example", [],
        "Scheme 1a, crossed lower arrow (previous work)", "attempted prior-work cycloaddition",
        "prior-work unreactive internal alkyne example", no_reaction_classification([ev(1), ev(2, "image")]),
        ["Cited prior-work negative result; quantities and recovery are not reported on these pages."])
    lit_fail["participants"] = [
        participant("lit_fail_yne_allenone", tpl1, "reactant", lit_fail["evidence_ids"]),
        participant("lit_fail_diphenylacetylene", rmids["diphenylacetylene"], "reactant", lit_fail["evidence_ids"]),
        participant("lit_fail_bf3", rmids["BF3Et2O"], "catalyst", lit_fail["evidence_ids"]),
    ]
    lit_fail["outcomes"] = [outcome("lit_fail_no_reaction", "no_reaction", [], lit_fail["evidence_ids"], True, "crossed reaction arrow / remained unreactive")]
    reactions.append(lit_fail)

    # One generic current-work scheme; Scheme 2 repeats the same general arrow as scope header.
    general_id = f"{DOC}:scheme1b_current_generic"
    general = reaction(
        "scheme1b_current_generic", [ev(2, "image"), ev(2), ev(3, "image")], "generic_scheme", [],
        "Scheme 1b and repeated scope header in Scheme 2", "tricyclization of yne-allenones",
        "BF3·Et2O-catalyzed tricyclization cascade", reactive_classification([ev(2, "image"), ev(2), ev(3, "image")]),
        ["R-group scheme; exact structures, quantities beyond shown equivalents, yield, time, and complete atom mapping are unresolved."])
    general["participants"] = [
        participant("gen_1", tpl1, "reactant", general["evidence_ids"], equivalents=q(1, "equiv"), reference=tpl1),
        participant("gen_2", tpl2, "reactant", general["evidence_ids"], equivalents=q(2, "equiv"), reference=tpl1),
        participant("gen_bf3", rmids["BF3Et2O"], "catalyst", general["evidence_ids"]),
        participant("gen_dce", rmids["DCE"], "solvent", general["evidence_ids"]),
        participant("gen_3", tpl3, "product", general["evidence_ids"], product_role="racemic product family"),
    ]
    set_solvent(general, "DCE", None, general["evidence_ids"])
    general["steps"] = [step("generic_reaction", general["evidence_ids"], [p["participant_id"] for p in general["participants"]], temperature=q(60, "°C"), atmosphere="Ar")]
    reactions.append(general)

    # Table 1 optimization rows.
    opt_rows = [
        (1, "BF3Et2O", 10, "DCE", 50, 75), (2, "BF3Et2O", 10, "DCM", 50, 65),
        (3, "BF3Et2O", 10, "MeCN", 50, 43), (4, "BF3Et2O", 10, "dioxane", 50, 25),
        (5, "BF3Et2O", 10, "EtOH", 50, None), (6, "AgOTf", 10, "DCE", 50, 67),
        (7, "ScOTf3", 10, "DCE", 50, 39), (8, "SnOTf2", 10, "DCE", 50, 46),
        (9, "FeCl3", 10, "DCE", 50, 35), (10, "YOTf3", 10, "DCE", 50, 55),
        (11, "BF3Et2O", 10, "DCE", 60, 78), (12, "BF3Et2O", 10, "DCE", 40, 64),
        (13, "BF3Et2O", 20, "DCE", 60, 72), (14, "BF3Et2O", 10, "DCE", 60, 83),
        (15, "BF3Et2O", 10, "DCE", 60, 87), (16, None, None, "DCE", 60, None),
    ]
    for entry, catalyst, loading, solvent_key, temp, yld in opt_rows:
        evidence = [ev(2, "image"), ev(2)]
        dry = entry in (15, 16)
        ratio = 2.2 if entry >= 14 else 2.0
        r = reaction(
            f"table1_entry_{entry:02d}", evidence, "specific_experiment", ["optimization"],
            f"Table 1 entry {entry}", "optimization of reaction conditions for 3a",
            "tricyclization of 1a and 2a" if yld is not None else "attempted tricyclization of 1a and 2a",
            reactive_classification(evidence) if yld is not None else no_reaction_classification(evidence),
            ["No SI was supplied; complete workup, purification, product structure encoding, and complete atom mapping remain unavailable."])
        r["general_scheme_ids"] = [general_id]
        p1 = participant("opt_1a", m1a, "reactant", evidence, amount=q(0.2, "mmol"), equivalents=q(1, "equiv"), reference=m1a)
        p2 = participant("opt_2a", m2a, "reactant", evidence,
                         amount=q(0.44 if ratio == 2.2 else 0.4, "mmol"), equivalents=q(ratio, "equiv"), reference=m1a)
        ps = participant("opt_solvent", rmids[solvent_key], "solvent", evidence, volume=q(2.5, "mL"))
        r["participants"] = [p1, p2, ps]
        if catalyst:
            r["participants"].append(participant("opt_catalyst", rmids[catalyst], "catalyst", evidence,
                                                   mol_percent=q(loading, "mol%"), reference=m1a))
        if yld is not None:
            r["participants"].append(participant("opt_3a", mids_3["3a"], "product", evidence, product_role="main product mixture"))
        set_solvent(r, ("dry " if dry else "") + reagent_specs[solvent_key][0], q(2.5, "mL"), evidence)
        r["steps"] = [step("reaction", evidence, [p["participant_id"] for p in r["participants"]],
                           operation="stir under the tabulated conditions", temperature=q(temp, "°C"),
                           time=q(12, "h"), atmosphere="argon")]
        if yld is not None:
            r["outcomes"] = [outcome(f"entry{entry}_yield", "yield", [mids_3["3a"]], evidence,
                                     q(yld, "%"), f"{yld}%", product_role="main product mixture",
                                     yield_type="isolated yield based on substrate 1a", measurement="isolated")]
            if entry == 1:
                add_selectivity(r, mids_3["3a"], "4:1", evidence)
        else:
            r["outcomes"] = [outcome(f"entry{entry}_no_reaction", "no_reaction", [], evidence, True, "N.R.")]
            if entry == 16:
                r["outcomes"].append(outcome("entry16_recovery", "recovery", [m1a, m2a], evidence,
                                                     q(100, "%"), "1a and 2a in 100% recovery",
                                                     scope="combined statement for both starting materials"))
        if entry == 12:
            r["notes"].append(fact("starting materials were not completely consumed, even after prolonging the reaction time",
                                    "incomplete consumption; prolonged duration not quantified", evidence))
            r["review"]["unresolved"].append("Table footnote states 12 h, while prose says entry 12 was prolonged; the prolonged duration is unreported.")
        reactions.append(r)

    # Scheme 2 scope: exact rows transcribed from the original image.
    scope_rows = [
        ("3a", "1a", "2a", 81, "4:1"), ("3b", "1b", "2a", 80, "10:1"),
        ("3c", "1c", "2a", 74, "10:1"), ("3d", "1d", "2a", 90, "10:1"),
        ("3e", "1e", "2a", 83, "4:1"), ("3f", "1f", "2a", 85, "4:1"),
        ("3g", "1g", "2a", 82, "5:1"), ("3h", "1h", "2a", 83, "10:1"),
        ("3i", "1i", "2a", 76, "4:1"), ("3j", "1j", "2a", 72, "10:1"),
        ("3k", "1k", "2a", 61, ">19:1"), ("3l", "1l", "2a", 0, None),
        ("3m", "1m", "2a", 55, "4:1"), ("3n", "1n", "2a", 80, ">19:1"),
        ("3o", "1o", "2a", 75, "6:1"), ("3p", "1p", "2a", 90, "5:1"),
        ("3q", "1q", "2a", 84, "5:1"), ("3r", "1r", "2a", 71, "10:1"),
        ("3s", "1s", "2a", 68, ">19:1"), ("3t", "1t", "2a", 70, "5:1"),
        ("3u", "1a", "2b", 70, ">19:1"), ("3v", "1a", "2c", 50, ">19:1"),
        ("3w", "1a", "2d", 70, "10:1"), ("3x", "1a", "2e", 65, "5:1"),
        ("3y", "1a", "2f", 76, "10:1"), ("3z", "1a", "2g", 50, "2:1"),
        ("3aa", "1a", "2h", 62, "3:1"), ("3bb", "1a", "2-cyclopropyl", 0, None),
        ("3cc", "1a", "2-tert-butyl", 0, None), ("3dd", "1c", "2f", 65, "10:1"),
    ]
    for product_label, one_label, two_label, yld, ratio in scope_rows:
        evidence = [ev(3, "image"), ev(3), ev(2)]
        if product_label == "3a": evidence.append(ev(5))
        success = yld > 0
        r = reaction(
            f"scheme2_scope_{product_label}", evidence, "specific_experiment", ["substrate_scope"],
            f"Scheme 2, product {product_label}", "substrate scope for the synthesis of 3",
            "BF3·Et2O-catalyzed tricyclization" if success else "attempted BF3·Et2O-catalyzed tricyclization",
            reactive_classification(evidence) if success else no_reaction_classification(evidence),
            ["No SI was supplied; scope yield measurement basis and full product SMILES/atom mapping are unreported."])
        r["general_scheme_ids"] = [general_id]
        one_mid, two_mid = mids_1[one_label], mids_2[two_label]
        r["participants"] = [
            participant("scope_1", one_mid, "reactant", evidence, amount=q(0.2, "mmol"), equivalents=q(1, "equiv"), reference=one_mid, inherited=True),
            participant("scope_2", two_mid, "reactant", evidence, amount=q(0.44, "mmol"), equivalents=q(2.2, "equiv"), reference=one_mid, inherited=True),
            participant("scope_bf3", rmids["BF3Et2O"], "catalyst", evidence, mol_percent=q(10, "mol%"), reference=one_mid, inherited=True),
            participant("scope_dce", rmids["DCE"], "solvent", evidence, volume=q(2.5, "mL"), inherited=True),
        ]
        if success:
            r["participants"].append(participant("scope_product", mids_3[product_label], "product", evidence, product_role="reported product mixture"))
        set_solvent(r, "dry 1,2-dichloroethane", q(2.5, "mL"), evidence, inherited=True)
        r["steps"] = [step("reaction", evidence, [p["participant_id"] for p in r["participants"]],
                           operation="stir under the standard conditions", temperature=q(60, "°C"),
                           time=q(12, "h"), atmosphere="argon", inherited=True)]
        if success:
            r["outcomes"] = [outcome(f"{product_label}_yield", "yield", [mids_3[product_label]], evidence,
                                     q(yld, "%"), f"{yld}%", product_role="reported product mixture",
                                     scope=f"Scheme 2 product {product_label}")]
            add_selectivity(r, mids_3[product_label], ratio, evidence)
        else:
            r["outcomes"] = [outcome(f"{product_label}_zero_target_yield", "other", [], evidence,
                                     q(0, "%"), f"{product_label} ({yld}%)",
                                     scope="reported desired-product yield; no produced molecule instantiated")]
            r["notes"].append(fact(
                "failed to give the desired product in tractable amount / did not give access to the desired product",
                "target product not obtained; this does not establish complete absence of all reaction", evidence))
            r["review"]["unresolved"].append("The source reports 0% desired product but does not fully characterize conversion, recovery, or byproducts.")
        if product_label == "3a":
            # Page-5 concrete procedure belongs to 3a and is not counted again.
            direct = {p["participant_id"]: p for p in r["participants"]}
            direct["scope_1"]["charges"]["mass"] = fact("48.8 mg", q(48.8, "mg"), [ev(5)])
            direct["scope_1"]["charges"]["amount"] = fact("0.2 mmol", q(0.2, "mmol"), [ev(5)])
            direct["scope_1"]["charges"]["equivalents"] = fact("1 equiv.", q(1, "equiv"), [ev(5)])
            direct["scope_1"]["equivalence_reference"] = fact("1a", one_mid, [ev(5)])
            direct["scope_2"]["charges"]["mass"] = fact("107.36 mg", q(107.36, "mg"), [ev(5)])
            direct["scope_2"]["charges"]["amount"] = fact("0.44 mmol", q(0.44, "mmol"), [ev(5)])
            direct["scope_2"]["charges"]["equivalents"] = fact("2.2 equiv.", q(2.2, "equiv"), [ev(5)])
            direct["scope_2"]["equivalence_reference"] = fact("1a", one_mid, [ev(5)])
            direct["scope_bf3"]["charges"]["mass"] = fact("2.82 mg", q(2.82, "mg"), [ev(5)])
            direct["scope_bf3"]["charges"]["mol_percent"] = fact("10 mol%", q(10, "mol%"), [ev(5)])
            direct["scope_bf3"]["equivalence_reference"] = fact("1a", one_mid, [ev(5)])
            direct["scope_dce"]["charges"]["volume"] = fact("2.5 mL", q(2.5, "mL"), [ev(5)])
            r["solvent_conditions"]["composition"] = fact("dry 1,2-dichloroethane", ["dry 1,2-dichloroethane"], [ev(5)])
            r["solvent_conditions"]["total_volume"] = fact("2.5 mL", q(2.5, "mL"), [ev(5)])
            r["steps"][0] = step(
                "reaction", [ev(5)], [p["participant_id"] for p in r["participants"]],
                operation="successively add 1a, 2a, BF3·Et2O, and dry DCE; stir",
                temperature=q(60, "°C"), time=q(12, "h"), atmosphere="argon")
            r["procedure"]["original_text"] = fact(
                "To a 10 mL Schlenk tube under Ar conditions, 1a, 2a, BF3·Et2O, and dry 1,2-dichloroethane were successively added; stirred at 60 °C for 12 h; concentrated in vacuum; column chromatography (petroleum ether/ethyl acetate = 5:1).",
                "page-5 general procedure for 3a", [ev(5)])
            r["procedure"]["addition_order"] = fact(
                "1a, 2a, BF3·Et2O, and dry DCE were successively added",
                ["1a", "2a", "BF3·Et2O", "dry DCE"], [ev(5)])
            r["special_conditions"]["reactor"] = fact("10 mL Schlenk tube", "10 mL Schlenk tube", [ev(5)])
            r["participants"] += [
                participant("scope_petroleum_ether", rmids["petroleum_ether"], "purification_agent", [ev(5)], stage="purification"),
                participant("scope_ethyl_acetate", rmids["ethyl_acetate"], "purification_agent", [ev(5)], stage="purification"),
                participant("scope_silica", rmids["silica_gel"], "purification_agent", [ev(5)], stage="purification"),
            ]
            r["steps"] += [
                step("concentrate", [ev(5)], ["scope_product"], operation="concentrate the reaction mixture in vacuum", stage="workup"),
                step("column", [ev(5)], ["scope_product", "scope_petroleum_ether", "scope_ethyl_acetate", "scope_silica"],
                     operation="silica-gel column chromatography; petroleum ether/ethyl acetate = 5:1", stage="purification"),
            ]
            r["notes"].append(fact("product 3a as light-yellow solid", "light-yellow solid", [ev(5)]))
        reactions.append(r)

    # Scale-up experiment.
    evidence = [ev(4, "image"), ev(4), ev(2)]
    scale = reaction(
        "scheme3_scale_up_3a", evidence, "specific_experiment", ["scale_up"],
        "Scheme 3", "scale-up transformation of 3a", "BF3·Et2O-catalyzed tricyclization",
        reactive_classification(evidence),
        ["Dry-DCE volume, yield measurement basis, Z/E ratio, workup, purification, and complete atom mapping are unreported for the scale-up."])
    scale["general_scheme_ids"] = [general_id]
    scale["participants"] = [
        participant("scale_1a", m1a, "reactant", evidence, amount=q(4, "mmol"), equivalents=q(1, "equiv"), reference=m1a),
        participant("scale_2a", m2a, "reactant", evidence, amount=q(8.8, "mmol"), equivalents=q(2.2, "equiv"), reference=m1a),
        participant("scale_bf3", rmids["BF3Et2O"], "catalyst", evidence, mol_percent=q(10, "mol%"), reference=m1a, inherited=True),
        participant("scale_dce", rmids["DCE"], "solvent", evidence),
        participant("scale_3a", mids_3["3a"], "product", evidence, product_role="reported product mixture"),
    ]
    scale["reaction_scale"] = fact("1a (4 mmol)", {"reference_molecule_id": m1a, "amount": q(4, "mmol")}, evidence)
    set_solvent(scale, "dry 1,2-dichloroethane", None, evidence)
    scale["steps"] = [step("scale_reaction", evidence, [p["participant_id"] for p in scale["participants"]],
                           operation="perform under the standard conditions", temperature=q(60, "°C"), time=q(12, "h"),
                           atmosphere="argon", inherited=True)]
    scale["outcomes"] = [
        outcome("scale_mass", "product_amount", [mids_3["3a"]], evidence, q(2.29, "g"), "2.29 g", product_role="reported product mixture"),
        outcome("scale_yield", "yield", [mids_3["3a"]], evidence, q(80, "%"), "80%", product_role="reported product mixture"),
    ]
    reactions.append(scale)

    # Scheme 4 concrete mechanistic controls.
    for suffix, partner_key, panel in (("phenol4", "control4", "a"), ("protected5", "control5", "b")):
        evidence = [ev(4, "image"), ev(4), ev(2)]
        r = reaction(
            f"scheme4_control_{suffix}", evidence, "specific_experiment", ["control", "mechanistic_study"],
            f"Scheme 4{panel}", "preliminary mechanistic control", "attempted tricyclization control",
            no_reaction_classification(evidence),
            ["The pages do not report exact recovery, conversion, analytical detection method, workup, or byproducts."])
        partner_mid = rmids[partner_key]
        r["general_scheme_ids"] = [general_id]
        r["participants"] = [
            participant("control_1a", m1a, "reactant", evidence, amount=q(0.2, "mmol"), equivalents=q(1, "equiv"), reference=m1a, inherited=True),
            participant("control_partner", partner_mid, "reactant", evidence, amount=q(0.44, "mmol"), equivalents=q(2.2, "equiv"), reference=m1a, inherited=True),
            participant("control_bf3", rmids["BF3Et2O"], "catalyst", evidence, mol_percent=q(10, "mol%"), reference=m1a, inherited=True),
            participant("control_dce", rmids["DCE"], "solvent", evidence, volume=q(2.5, "mL"), inherited=True),
        ]
        set_solvent(r, "dry 1,2-dichloroethane", q(2.5, "mL"), evidence, inherited=True)
        r["steps"] = [step("control_reaction", evidence, [p["participant_id"] for p in r["participants"]],
                           operation="stir under the standard conditions", temperature=q(60, "°C"), time=q(12, "h"),
                           atmosphere="argon", inherited=True)]
        r["outcomes"] = [outcome("control_no_reaction", "no_reaction", [], evidence, True, "no reaction")]
        if partner_key == "control4":
            r["outcomes"].append(outcome("control_recovery", "recovery", [m1a, partner_mid], evidence,
                                                 "starting materials recovered", "starting materials were recovered",
                                                 scope="qualitative combined recovery statement"))
        reactions.append(r)

    # Scheme 5 is a proposed pathway, not a series of executed experiments.
    evidence = [ev(4, "image"), ev(4)]
    mech = reaction(
        "scheme5_proposed_pathway", evidence, "proposed_mechanism", [],
        "Scheme 5", "plausible reaction pathway", "proposed double [2+2] relay / ring opening / oxa-cyclization / dehydration",
        reactive_classification(evidence),
        ["Intermediates A--F and E' are proposed, not isolated experiments; exact structures and atom mapping are not encoded."])
    mech["participants"] = [
        participant("mech_1a", m1a, "reactant", evidence),
        participant("mech_2a", m2a, "reactant", evidence),
        participant("mech_bf3", rmids["BF3Et2O"], "catalyst", evidence),
        participant("mech_3a", mids_3["3a"], "product", evidence, product_role="proposed pathway endpoint"),
    ]
    mech["notes"] = [
        fact("1a first undergoes intramolecular [2+2] cycloaddition to A", "proposed A formation", evidence),
        fact("A and 2a give B; ring opening gives C", "proposed B/C sequence", evidence),
        fact("VQM D from 2a is intercepted by C; oxa-cyclization and dehydration give major Z-3a", "proposed path a to major Z isomer", evidence),
        fact("E' gives the minor E stereoisomer via path b", "proposed path b to minor E isomer", evidence),
    ]
    reactions.append(mech)

    bundle["reactions"] = reactions
    bundle["source_context"] = [
        {"context_id": f"{DOC}:context:paper_metadata", "document_id": DOC,
         "evidence_ids": [ev(1)], "kind": "paper_metadata",
         "content": {"title": TITLE, "doi": "10.1002/adsc.202100061", "pages_reviewed": [1, 2, 3, 4, 5, 6]}},
        {"context_id": f"{DOC}:context:si_notice", "document_id": DOC,
         "evidence_ids": [ev(1)], "kind": "availability_notice",
         "content": "The main paper advertises supporting information online, but no SI file was supplied or associated in this workspace."},
        {"context_id": f"{DOC}:context:xray_3b", "document_id": DOC,
         "evidence_ids": [ev(3, "image"), ev(3), ev(6)], "kind": "structure_assignment",
         "content": "The stereo-structure of 3b was confirmed by single-crystal X-ray diffraction (CCDC 2054923); derivatives were assigned by analogy."},
        {"context_id": f"{DOC}:context:no_spectra", "document_id": DOC,
         "evidence_ids": [ev(1), ev(2), ev(3), ev(4), ev(5), ev(6)], "kind": "negative_inventory",
         "content": "No standalone mass-spectrum peak table, molecular-weight table, or fragment-relative-abundance block occurs on pages 1--6."},
    ]

    # Raw model inventory and literal mention hints are kept separate from submitted records.
    raw_inventory = {
        "document_id": DOC,
        "pages": [
            {"page": 1, "inventory": ["title/metadata", "abstract", "introduction", "prose describing cited prior work"],
             "decisions": [f"{DOC}:scheme1a_prior_success", f"{DOC}:scheme1a_prior_diphenylacetylene_failure"]},
            {"page": 2, "inventory": ["Scheme 1a: two prior-work arrows", "Scheme 1b: one current generic arrow", "Table 1: 16 optimization rows"],
             "decisions": [f"{DOC}:scheme1a_prior_success", f"{DOC}:scheme1a_prior_diphenylacetylene_failure", general_id] + [f"{DOC}:table1_entry_{i:02d}" for i in range(1, 17)]},
            {"page": 3, "inventory": ["Scheme 2 generic header", "30 scope outcomes 3a--3dd", "Figure 1 X-ray of 3b"],
             "decisions": [general_id] + [f"{DOC}:scheme2_scope_{row[0]}" for row in scope_rows]},
            {"page": 4, "inventory": ["Scheme 3 scale-up", "Scheme 4a/4b controls", "Scheme 5 proposed pathway"],
             "decisions": [f"{DOC}:scheme3_scale_up_3a", f"{DOC}:scheme4_control_phenol4", f"{DOC}:scheme4_control_protected5", f"{DOC}:scheme5_proposed_pathway"]},
            {"page": 5, "inventory": ["general procedure for 3a", "acknowledgements", "references 1--15"],
             "decisions": [f"{DOC}:scheme2_scope_3a"]},
            {"page": 6, "inventory": ["references 16--18", "CCDC 2054923 notice"], "decisions": []},
        ],
        "limitations": ["No SI supplied.", "Complex scope products and small positional variants were not assigned concrete SMILES.",
                        "Machine records remain needs_review; page inspection does not establish chemistry accuracy."],
    }
    st.write(RUN / "raw_model/page_extraction.json", raw_inventory)
    st.write(RUN / "raw_model/stereo_checks.json", stereo_checks)
    st.write(RUN / "raw_model/source_mentions.json", st.source_mentions(data, ["1 a", "2 a", "3 a", "3 b"], 220))

    coverage = {
        "document_id": DOC, "semantic_accuracy_verified": False,
        "pages": [
            {"page": 1, "state": "attempted", "regions": [
                {"region": "title, abstract, and both introduction columns", "finding": "Paper metadata and reaction overview; prose identifies successful cited double-[2+2] work and the cited diphenylacetylene limitation; no new concrete experiment arrow on this page.", "evidence_ids": [ev(1), ev(1, "image")]}
            ], "record_ids": [f"{DOC}:scheme1a_prior_success", f"{DOC}:scheme1a_prior_diphenylacetylene_failure"],
             "unresolved": ["The advertised supporting information was not supplied and was not associated."]},
            {"page": 2, "state": "attempted", "regions": [
                {"region": "Scheme 1a, upper-left", "finding": "Two cited prior-work outcomes: a successful generic cycloaddition relay and an explicitly crossed/unreactive diphenylacetylene case.", "evidence_ids": [ev(2, "image"), ev(2)]},
                {"region": "Scheme 1b, left-middle", "finding": "Current-work R-group generic tricyclization; retained as one generic scheme and linked to the repeated Scheme 2 header.", "evidence_ids": [ev(2, "image"), ev(2)]},
                {"region": "Table 1, upper-right and footnotes", "finding": "All 16 optimization rows recorded separately, including entries 5 and 16 as N.R. and entry 16 as 100% starting-material recovery.", "evidence_ids": [ev(2, "image"), ev(2)]},
                {"region": "optimization prose below table", "finding": "Entry 1 Z/E ratio and entry 12 incomplete consumption/prolonged-time statement associated with their rows; standard entry 15 conditions retained.", "evidence_ids": [ev(2, "image"), ev(2)]}
            ], "record_ids": [f"{DOC}:scheme1a_prior_success", f"{DOC}:scheme1a_prior_diphenylacetylene_failure", general_id] + [f"{DOC}:table1_entry_{i:02d}" for i in range(1, 17)],
             "unresolved": ["Entry 12 prose says the time was prolonged but does not state the prolonged duration."]},
            {"page": 3, "state": "attempted", "regions": [
                {"region": "Scheme 2 header and product grid", "finding": "All 30 printed scope entries 3a--3dd transcribed individually; 3l, 3bb, and 3cc retain 0% desired-product outcomes without invented product molecules.", "evidence_ids": [ev(3, "image"), ev(3)]},
                {"region": "scope discussion in both text columns", "finding": "Substrate-label relationships and the special 1c + 2f combination for 3dd associated with the corresponding records.", "evidence_ids": [ev(3, "image"), ev(3)]},
                {"region": "Figure 1 and adjacent prose", "finding": "3b X-ray/CCDC assignment retained as structure-assignment context, not a synthetic experiment.", "evidence_ids": [ev(3, "image"), ev(3), ev(6)]}
            ], "record_ids": [general_id] + [f"{DOC}:scheme2_scope_{row[0]}" for row in scope_rows],
             "unresolved": ["Scope yield measurement bases and Z/E determination methods are not stated on the supplied pages.", "Small complex product drawings and several positional isomers were not encoded as concrete SMILES."]},
            {"page": 4, "state": "attempted", "regions": [
                {"region": "Scheme 3 and opening prose", "finding": "One 4 mmol scale-up experiment recorded with 2.29 g product and 80% yield.", "evidence_ids": [ev(4, "image"), ev(4)]},
                {"region": "Scheme 4a and 4b", "finding": "Two real mechanistic control experiments recorded separately; both report no reaction, and panel a additionally reports recovered starting materials without a percentage.", "evidence_ids": [ev(4, "image"), ev(4)]},
                {"region": "Scheme 5 and mechanism prose", "finding": "One proposed pathway record covers intermediates A--F/E' and paths a/b; intermediates were not misrepresented as executed experiments.", "evidence_ids": [ev(4, "image"), ev(4)]}
            ], "record_ids": [f"{DOC}:scheme3_scale_up_3a", f"{DOC}:scheme4_control_phenol4", f"{DOC}:scheme4_control_protected5", f"{DOC}:scheme5_proposed_pathway"],
             "unresolved": ["Scale-up solvent volume, yield basis, Z/E ratio, and workup are unreported.", "Control recovery percentages and analytical detection methods are unreported."]},
            {"page": 5, "state": "attempted", "regions": [
                {"region": "Experimental Section, General Procedure", "finding": "Concrete 1a/2a/BF3·Et2O charges, 10 mL Schlenk reactor, dry DCE, 60 °C/12 h, TLC monitoring, vacuum concentration, and 5:1 petroleum ether/ethyl acetate silica chromatography linked to scope product 3a rather than duplicated.", "evidence_ids": [ev(5), ev(5, "image")]},
                {"region": "acknowledgements and references 1--15", "finding": "No additional experiment; references were not treated as this paper's experiments.", "evidence_ids": [ev(5), ev(5, "image")]}
            ], "record_ids": [f"{DOC}:scheme2_scope_3a"],
             "unresolved": ["The procedure does not state the isolated mass or yield of the described 3a preparation on this page."]},
            {"page": 6, "state": "attempted", "regions": [
                {"region": "references 16--18 and blank lower page", "finding": "CCDC 2054923 notice and literature references only; no synthesis experiment, table, or reaction arrow.", "evidence_ids": [ev(6), ev(6, "image")]}
            ], "record_ids": [], "unresolved": ["Crystallographic data are referenced externally but were not supplied as an associated source file."]},
        ],
    }
    st.write(RUN / "coverage/pages_0001_0006.json", coverage)

    # Mark only the assigned page tasks attempted; this does not claim semantic completeness.
    tasks = st.read(RUN / "tasks.json")
    for task in tasks:
        if task["document_id"] == DOC and task["task_id"] in {ev(i) for i in range(1, 7)}:
            task["state"] = "attempted"
            task["reason"] = "Original page image and extracted text inspected; findings mapped in coverage/pages_0001_0006.json."
    st.write(RUN / "tasks.json", tasks)

    st.write(SUBMITTED, bundle)
    print(json.dumps({"molecules": len(mols), "reactions": len(reactions),
                      "source_context": len(bundle["source_context"])}, ensure_ascii=False))


if __name__ == "__main__":
    main()
