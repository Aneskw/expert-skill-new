#!/usr/bin/env python3
"""Source-specific v3.1 extraction builder for doc_739c1c0b69cb742a.

This file records model transcription/interpretation separately from the merged
dataset.  It imports the project helpers but does not modify tool or schema code.
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

RUN = ROOT / "runs/current"
DOC = "doc_739c1c0b69cb742a"
E = {p: f"{DOC}:page{p}:image" for p in range(1, 5)}
T = {p: f"{DOC}:page{p}:text" for p in range(1, 5)}


def f(raw=None, normalized=None, evidence=(), method="original_extraction", derivation=None,
      status=None):
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   derivation=derivation)


def missing(evidence=()):
    return st.fact(evidence=evidence)


def derived(normalized, evidence, derivation, raw=None):
    return f(raw, normalized, evidence, "structure_derivation", derivation)


def inferred(normalized, evidence, derivation, raw=None):
    return f(raw, normalized, evidence, "inference", derivation)


def not_applicable(evidence, reason):
    return f(evidence=evidence, method="not_applicable", status="not_applicable",
             derivation=reason)


def mid(label):
    safe = re.sub(r"[^A-Za-z0-9]+", "_", label).strip("_").lower()
    return f"{DOC}:mol:{safe}"


bundle = json.loads((RUN / "submitted.json").read_text())
bundle["molecules"] = []
bundle["reactions"] = []
bundle["source_context"] = []
molecules = {}


def add_occurrence(m, label, evidence):
    evidence = list(dict.fromkeys(evidence))
    existing = next((o for o in m["occurrences"] if o["document_id"] == DOC and
                     o["compound_label"]["normalized_value"] == label), None)
    if existing:
        existing["evidence_ids"] = list(dict.fromkeys(existing["evidence_ids"] + evidence))
        existing["compound_label"]["evidence_ids"] = list(dict.fromkeys(
            existing["compound_label"]["evidence_ids"] + evidence))
    else:
        m["occurrences"].append({
            "document_id": DOC,
            "compound_label": f(label, label, evidence),
            "evidence_ids": evidence,
        })


def set_stereo(m, evidence, chirality="unresolved", geometry="unresolved",
               representation="unresolved", atom_centers=None, double_bonds=None,
               axial=None, mixture=None, assignment=None):
    ev = list(evidence)
    reason = "Status assigned conservatively from the supplied source drawing; no absent descriptor is treated as achiral or racemic."
    m["stereochemistry"]["chirality_status"] = inferred(chirality, ev, reason)
    m["stereochemistry"]["geometry_status"] = inferred(geometry, ev, reason)
    m["stereochemistry"]["representation_status"] = inferred(representation, ev, reason)
    if atom_centers is not None:
        m["stereochemistry"]["atom_centers"] = f(atom_centers, atom_centers, ev)
    if double_bonds is not None:
        m["stereochemistry"]["double_bonds"] = f(double_bonds, double_bonds, ev)
    if axial is not None:
        m["stereochemistry"]["axial"] = f(axial, axial, ev)
    if mixture is not None:
        m["stereochemistry"]["mixture"] = f(mixture, mixture, ev)
    if assignment is not None:
        m["stereochemistry"]["assignment"] = f(assignment, assignment, ev)


def add_molecule(label, evidence, name=None, smiles=None, stereo=None, identity_note=None):
    identity = mid(label)
    if identity in molecules:
        add_occurrence(molecules[identity], label, evidence)
        return identity
    m = st.new_molecule(identity, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"] = []
    add_occurrence(m, label, evidence)
    m["names"]["original"] = f(label, name or label, evidence)
    if name:
        m["names"]["english"] = f(name, name, evidence)
    m["cas_numbers"] = missing(evidence)
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = derived(smiles, evidence,
            "Candidate connectivity transcribed from the supplied drawing/name and checked with the standard stereo command; parsing is not source-identity certification.",
            raw={"source_label": label, "source_representation": "drawing or explicit chemical name"})
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing(evidence)
    set_stereo(m, evidence, **(stereo or {}))
    m["identity_note"] = f(identity_note or f"Source-labelled molecule {label}; complete structure was not asserted where the small drawing or stereochemistry could not be safely encoded.",
                                identity_note or f"Source-labelled molecule {label}; complete structure was not asserted where the small drawing or stereochemistry could not be safely encoded.", evidence)
    molecules[identity] = m
    bundle["molecules"].append(m)
    return identity


def add_template_molecule(label, evidence, fixed, reason, stereo=None):
    identity = mid(label)
    if identity in molecules:
        add_occurrence(molecules[identity], label, evidence)
        return identity
    m = st.new_molecule(identity, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"] = []
    add_occurrence(m, label, evidence)
    m["names"]["original"] = f(label, label, evidence)
    m["cas_numbers"] = missing(evidence)
    m["representation_type"] = "generic_template"
    m["smiles"] = missing(evidence)
    m["generic_template"] = {
        "fixed_structure": f(fixed, fixed, evidence),
        "variables": [],
        "instance_mapping": [],
        "reason": f(reason, reason, evidence),
    }
    set_stereo(m, evidence, **(stereo or {}))
    m["identity_note"] = f(reason, reason, evidence)
    molecules[identity] = m
    bundle["molecules"].append(m)
    return identity


# Recurring concrete/candidate identities.
add_molecule("1a", [E[1], E[2], E[3], E[4]], "N-phenoxyacetamide",
             "CC(=O)NOc1ccccc1",
             {"chirality": "achiral", "geometry": "not_applicable",
              "representation": "not_applicable"},
             "The complete 1a drawing is legible. The candidate SMILES is a graph transcription; no CAS was supplied.")
add_molecule("2a", [E[1], E[2], E[3], E[4]], "α-allenol 2a",
             "CC(O)C=C=CC1CCCCC1",
             {"chirality": "unspecified", "geometry": "not_applicable",
              "representation": "not_encoded",
              "atom_centers": [{"feature": "secondary alcohol carbon", "assignment": "unspecified"}],
              "axial": [{"feature": "allene axis", "assignment": "unspecified"}]},
             "Candidate connectivity retains the allene and secondary alcohol. The source does not assign the alcohol center or allene axis.")
add_molecule("3aa", [E[2], E[3], E[4]], "4-alkylidene chroman-2-ol 3aa",
             "CC1(O)CC(=CC2CCCCC2)c2ccccc2O1",
             {"chirality": "unspecified", "geometry": "specified",
              "representation": "not_encoded",
              "atom_centers": [{"feature": "chroman-2-ol hemiketal carbon", "assignment": "unspecified"}],
              "double_bonds": [{"feature": "exocyclic alkene", "configuration": "E"}],
              "assignment": "Only E-selective products were reported; absolute configuration at the hemiketal carbon was not assigned."},
             "Connectivity candidate is retained without an E tag because the exact directional encoding was not independently established; the source E assignment remains explicit.")
add_molecule("3ea", [E[2], E[3]], "chlorinated 4-alkylidene chroman-2-ol 3ea",
             "CC1(O)CC(=CC2CCCCC2)c2cc(Cl)ccc2O1",
             {"chirality": "unspecified", "geometry": "specified",
              "representation": "not_encoded",
              "atom_centers": [{"feature": "chroman-2-ol hemiketal carbon", "assignment": "unspecified"}],
              "double_bonds": [{"feature": "exocyclic alkene", "configuration": "E"}]},
             "Candidate connectivity is used to retain the deuterium-study substrate; exact positional atom mapping and E encoding remain unresolved.")

simple = [
    ("DCE", "1,2-dichloroethane", "ClCCCl"),
    ("KOPiv", "potassium pivalate", "CC(C)(C)C(=O)[O-].[K+]"),
    ("CsOAc", "cesium acetate", "CC(=O)[O-].[Cs+]"),
    ("D2O", "deuterium oxide", "[2H]O[2H]"),
    ("CD3OD", "methanol-d4", "[2H]OC([2H])([2H])[2H]"),
    ("MeOH", "methanol", "CO"),
    ("LiAlH4", "lithium aluminium hydride", "[Li+].[AlH4-]"),
]
for label, name, smiles in simple:
    pages = [E[2], E[3], E[4]] if label in ("DCE", "KOPiv") else ([E[2]] if label == "CsOAc" else [E[3]] if label in ("D2O", "CD3OD") else [E[4]])
    chirality = "unresolved" if label == "LiAlH4" else "achiral"
    representation = "not_encoded" if label == "LiAlH4" else "not_applicable"
    add_molecule(label, pages, name, smiles,
                 {"chirality": chirality, "geometry": "not_applicable",
                  "representation": representation},
                 f"Identity is stated by the source abbreviation/name; graph candidate checked with the standard stereo command.")

add_molecule("[Cp*RhCl2]2", [E[1], E[2], E[3], E[4]],
             "pentamethylcyclopentadienyl rhodium dichloride dimer", None,
             {"chirality": "unresolved", "geometry": "unresolved", "representation": "unresolved"})
add_molecule("N-phenoxytosylamide", [E[2]], "N-phenoxytosylamide",
             "Cc1ccc(S(=O)(=O)NOc2ccccc2)cc1",
             {"chirality": "achiral", "geometry": "not_applicable", "representation": "not_applicable"},
             "Candidate graph expands the source name; the ineffective optimization details are unavailable without SI.")

# Scope-specific identities.  Structures remain source-labelled if exact graph/stereo encoding was not safe.
for letter in "bcdefghijklmnopqr":
    stereo = {"chirality": "partially_specified" if letter == "r" else "unresolved",
              "geometry": "unresolved", "representation": "not_encoded" if letter == "r" else "unresolved"}
    if letter == "r":
        stereo["assignment"] = "L-tyrosine-derived side-chain stereochemistry is indicated by the source; the absolute descriptor is not encoded here."
    add_molecule(f"1{letter}", [E[2]], f"N-phenoxyacetamide 1{letter}", stereo=stereo)
for letter in "bcdefghij":
    add_molecule(f"2{letter}", [E[2]], f"α-allenol 2{letter}", stereo={
        "chirality": "unspecified", "geometry": "not_applicable", "representation": "not_encoded",
        "atom_centers": [{"feature": "secondary alcohol carbon", "assignment": "unspecified"}],
        "axial": [{"feature": "allene axis", "assignment": "unspecified or unresolved from the drawing"}],
    })

scope_rows = [
    ("3aa", "1a", "2a", 73), ("3ba", "1b", "2a", 63), ("3ca", "1c", "2a", 79),
    ("3da", "1d", "2a", 68), ("3ea", "1e", "2a", 54), ("3fa", "1f", "2a", 67),
    ("3ga", "1g", "2a", 55), ("3ha", "1h", "2a", 62), ("3ia", "1i", "2a", 47),
    ("3ja", "1j", "2a", 33), ("3ka", "1k", "2a", 38), ("3la", "1l", "2a", 68),
    ("3ma", "1m", "2a", 46), ("3na", "1n", "2a", 69), ("3oa", "1o", "2a", 57),
    ("3pa", "1p", "2a", 54), ("3qa", "1q", "2a", 53), ("3ra", "1r", "2a", 47),
    ("3ab", "1a", "2b", 51), ("3ac", "1a", "2c", 43), ("3ad", "1a", "2d", 39),
    ("3ae", "1a", "2e", 37), ("3af", "1a", "2f", 65), ("3ag", "1a", "2g", 83),
    ("3ah", "1a", "2h", 43), ("3ai", "1a", "2i", 42), ("3aj", "1a", "2j", 39),
]
for product, _, _, _ in scope_rows:
    if product in molecules:
        continue
    p_stereo = {
        "chirality": "partially_specified" if product == "3ra" else "unspecified",
        "geometry": "specified", "representation": "not_encoded",
        "atom_centers": [{"feature": "chroman-2-ol hemiketal carbon", "assignment": "unspecified"}],
        "double_bonds": [{"feature": "exocyclic alkene", "configuration": "E"}],
        "assignment": ("L-tyrosine-derived side-chain configuration is depicted; the new hemiketal center remains unspecified."
                       if product == "3ra" else "Only E-selective product reported; no absolute configuration assigned at the hemiketal carbon."),
    }
    add_molecule(product, [E[2], E[3]], f"4-alkylidene chroman-2-ol {product}", stereo=p_stereo,
                 identity_note=("X-ray structure reported with CCDC 2084815; CIF was not supplied in this session."
                                if product == "3ai" else None))


def rx(suffix, evidence, record_type="specific_experiment", categories=None,
       location=None, reaction_name=None, transformation=None):
    rid = f"{DOC}:{suffix}"
    r = st.new_reaction(rid, DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, evidence) if categories is not None else missing(evidence)
    r["location"] = f(location, location, evidence) if location else missing(evidence)
    r["reaction_name"] = f(reaction_name, reaction_name, evidence) if reaction_name else missing(evidence)
    r["transformation_type"] = inferred(transformation, evidence,
        "Transformation class inferred from the source reaction arrow and product class.") if transformation else missing(evidence)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    return r


def participant(r, label, role, evidence, pid=None, stage="reaction"):
    identity = add_molecule(label, evidence)
    p = st.new_participant(pid or f"p_{re.sub(r'[^A-Za-z0-9]+', '_', label).strip('_').lower()}",
                           identity, role, evidence, stage=stage)
    if role == "product":
        p["product_role"] = f("product", "main product", evidence)
    r["participants"].append(p)
    return p


def charge(p, key, raw, normalized, evidence, reference=None, method="original_extraction"):
    p["charges"][key] = f(raw, normalized, evidence, method=method)
    if reference is not None:
        p["equivalence_reference"] = f(reference, mid(reference), evidence, method=method)


def reaction_step(r, step_id, evidence, participant_ids, operation, temperature=None,
                  time=None, atmosphere=None, stage="reaction"):
    s = st.new_step(step_id, evidence, stage=stage)
    s["participant_ids"] = participant_ids
    s["operation"] = f(operation, operation, evidence)
    if temperature is not None:
        s["temperature"] = f(temperature[0], temperature[1], evidence)
    if time is not None:
        s["time"] = f(time[0], time[1], evidence)
    if atmosphere is not None:
        s["atmosphere"] = f(atmosphere, atmosphere, evidence)
    r["steps"].append(s)
    return s


def outcome(r, oid, kind, targets, evidence, raw, normalized, product_role=None,
            yield_type=None, basis=None, scope=None):
    o = st.new_outcome(oid, kind, [mid(x) for x in targets], evidence)
    if product_role:
        o["product_role"] = f(product_role, product_role, evidence)
    o["value"] = f(raw, normalized, evidence)
    if yield_type:
        o["yield_type"] = f(yield_type, yield_type, evidence)
    if basis:
        o["measurement_basis"] = f(basis, basis, evidence)
    if scope:
        o["scope"] = f(scope, scope, evidence)
    r["outcomes"].append(o)
    return o


def set_classification(r, bonds, evidence, limitation):
    changes = []
    labels = {}
    for idx, (a, b, order, kind) in enumerate(bonds, 1):
        change = {
            "change_id": f"b{idx}",
            "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
            "before_order": 0 if kind in ("formed", "hydrogen_gain") else 1,
            "after_order": order,
            "change_type": kind,
        }
        if kind == "hydrogen_gain":
            change["hydrogen_count"] = 1
        changes.append(change)
        label = st.bond_label(a, b, order)
        labels.setdefault(label["code"], {**label, "supporting_change_ids": []})["supporting_change_ids"].append(change["change_id"])
    derivation = "Class-level net bond formation inferred from the supplied reactant/product drawings; atom identities and full mapping were not asserted."
    r["classification"]["bond_changes"] = [derived(c, evidence, derivation) for c in changes]
    r["classification"]["labels"] = derived(list(labels.values()), evidence, derivation)
    r["classification"]["coverage"] = inferred("partial", evidence, derivation)
    r["classification"]["limitations"] = [limitation,
        "No atom-mapped reaction SMILES was asserted; implicit-H source and omitted coproducts remain unresolved."]


def no_classification(r, evidence, coverage="not_applicable", reason="No observed source transformation supports a formed-bond label."):
    r["classification"]["bond_changes"] = []
    r["classification"]["labels"] = inferred([], evidence, reason)
    r["classification"]["coverage"] = inferred(coverage, evidence, reason)
    r["classification"]["limitations"] = [reason]


CAT = "[Cp*RhCl2]2"

# Scheme 1a: cited prior reaction families, not experiments of this paper.
r = rx("scheme1a_prior_work", [E[1]], "literature_example", None,
       "PDF page 1, Scheme 1a", "classic reaction modes of TM-catalyzed C-H couplings with α-allenols")
r["notes"] = [f("Scheme 1a, paths I and II and oxidative coupling are cited prior-work families.",
                  "Cited prior-work families; no current-paper experimental row or yield is claimed.", [E[1]])]
no_classification(r, [E[1]], "not_applicable", "Grouped literature examples are not executed experiments in this paper.")
r["review"]["unresolved"] = ["Individual cited examples and conditions were not expanded from the reference list."]
bundle["reactions"].append(r)

# Scheme 1c generic current transformation.
generic1 = add_template_molecule("generic N-phenoxyacetamide 1", [E[1]],
    {"scaffold": "N-phenoxyacetamide", "variables": ["Ar/R1"]},
    "Scheme 1c uses an Ar/R-group generic substrate; it is not a concrete molecule.")
generic2 = add_template_molecule("generic α-allenol 2", [E[1]],
    {"scaffold": "α-allenol", "variables": ["R1", "R2"]},
    "Scheme 1c uses R-group variables and cannot be a concrete SMILES.",
    {"chirality": "unresolved", "geometry": "not_applicable", "representation": "not_encoded",
     "axial": [{"feature": "allene axis", "assignment": "variable-dependent and unspecified"}]})
generic3 = add_template_molecule("generic chroman-2-ol 3", [E[1]],
    {"scaffold": "4-alkylidene chroman-2-ol", "variables": ["Ar", "R1", "R2"]},
    "Scheme 1c is a generic product family; E selectivity is discussed for concrete scope products on page 3.")
r = rx("scheme1c_general", [E[1]], "generic_scheme", None, "PDF page 1, Scheme 1c",
       "C-H [3+3] annulation", "Rh(III)-catalyzed redox-neutral [3+3] annulation")
for label, role in [("generic N-phenoxyacetamide 1", "reactant"), ("generic α-allenol 2", "reactant"),
                    (CAT, "catalyst"), ("KOPiv", "additive"), ("DCE", "solvent"),
                    ("generic chroman-2-ol 3", "product")]:
    participant(r, label, role, [E[1]])
catp = next(p for p in r["participants"] if p["molecule_id"] == mid(CAT))
charge(catp, "mol_percent", "5 mol%", {"number": 5, "unit": "mol%"}, [E[1]], "generic N-phenoxyacetamide 1")
kop = next(p for p in r["participants"] if p["molecule_id"] == mid("KOPiv"))
charge(kop, "equivalents", "1 equiv", {"number": 1, "unit": "equiv"}, [E[1]], "generic N-phenoxyacetamide 1")
r["solvent_conditions"]["composition"] = f("DCE", [mid("DCE")], [E[1]])
reaction_step(r, "scheme1c_arrow", [E[1]], [p["participant_id"] for p in r["participants"]],
              "generic reaction arrow", ("100 °C", {"number": 100, "unit": "°C"}),
              ("4 h", {"number": 4, "unit": "h"}))
set_classification(r, [("C", "C", 1, "formed"), ("C", "O", 1, "formed")], [E[1]],
                   "The generic product family shows one new aryl-allene C-C connection and one ring-closing C-O connection; exact variable atoms are unresolved.")
r["review"]["unresolved"] = ["Generic R groups and atom mapping are unresolved; this is not counted as a concrete experiment."]
bundle["reactions"].append(r)
GENERAL = r["reaction_id"]

# Page-2 prose-only initial CsOAc experiment.
r = rx("tableS1_entry1_CsOAc", [T[2], E[2]], categories=["optimization"],
       location="PDF page 2 prose; Table S1 entry 1 is cited but SI was not supplied",
       reaction_name="regioselective [3+3] annulation", transformation="redox-neutral [3+3] annulation")
r["general_scheme_ids"] = [GENERAL]
for label, role in [("1a", "reactant"), ("2a", "reactant"), (CAT, "catalyst"),
                    ("CsOAc", "additive"), ("DCE", "solvent"), ("3aa", "product")]:
    participant(r, label, role, [T[2], E[2]])
r["solvent_conditions"]["composition"] = f("DCE", [mid("DCE")], [T[2]])
reaction_step(r, "reaction", [T[2]], [p["participant_id"] for p in r["participants"]],
              "N-phenoxyacetamide 1a was reacted with α-allenol 2a in the presence of [Cp*RhCl2]2 and CsOAc in DCE")
outcome(r, "product_observed", "other", ["3aa"], [T[2]],
        "delivering an intriguing 4-alkylidene chroman-2-ol derivative 3aa", "3aa obtained",
        product_role="desired product", scope="qualitative outcome; numerical yield not reported on supplied page")
set_classification(r, [("C", "C", 1, "formed"), ("C", "O", 1, "formed")], [T[2], E[2]],
                   "Class-level [3+3] annulation assignment; exact atoms and bond-order changes await SI/source mapping.")
r["review"]["unresolved"] = ["Table S1/SI was not supplied; charges, temperature, time, yield, workup and outcome basis remain unreported here."]
bundle["reactions"].append(r)

# Prose-reported ineffective N-phenoxytosylamide optimization.
r = rx("N_phenoxytosylamide_ineffective", [T[2]], categories=["optimization", "control"],
       location="PDF page 2 optimization prose", reaction_name="N-phenoxytosylamide test")
r["general_scheme_ids"] = [GENERAL]
for label, role in [("N-phenoxytosylamide", "reactant"), ("2a", "reactant")]:
    participant(r, label, role, [T[2]])
reaction_step(r, "reported_test", [T[2]], [p["participant_id"] for p in r["participants"]],
              "N-phenoxytosylamide was tested and proved ineffective")
outcome(r, "ineffective", "no_reaction", [], [T[2]], "proved to be ineffective", "ineffective",
        scope="qualitative optimization result; expected product and detection limit not reported")
no_classification(r, [T[2]], "not_applicable", "Ineffective test; no observed product transformation is specified.")
r["review"]["unresolved"] = ["Exact conditions, conversion/recovery, expected product identity and detection method are only in unavailable SI/optimization data."]
bundle["reactions"].append(r)


def add_scope_record(product, aryl_sub, allene, yld):
    ev = [E[2], T[2], E[3], T[3]]
    r = rx(f"scheme2_{product}", ev, categories=["substrate_scope"],
           location=f"PDF page 2, Scheme 2, product {product}",
           reaction_name="C-H [3+3] annulation", transformation="redox-neutral [3+3] annulation")
    r["general_scheme_ids"] = [GENERAL]
    for label, role in [(aryl_sub, "reactant"), (allene, "reactant"), (CAT, "catalyst"),
                        ("KOPiv", "additive"), ("DCE", "solvent"), (product, "product")]:
        participant(r, label, role, ev)
    p1 = next(p for p in r["participants"] if p["molecule_id"] == mid(aryl_sub))
    p2 = next(p for p in r["participants"] if p["molecule_id"] == mid(allene))
    pc = next(p for p in r["participants"] if p["molecule_id"] == mid(CAT))
    pk = next(p for p in r["participants"] if p["molecule_id"] == mid("KOPiv"))
    charge(p1, "amount", "1 (0.2 mmol)", {"number": 0.2, "unit": "mmol"}, [E[2]], method="inherited_general_condition")
    charge(p2, "amount", "2 (0.24 mmol)", {"number": 0.24, "unit": "mmol"}, [E[2]], method="inherited_general_condition")
    charge(pc, "mol_percent", "[Cp*RhCl2]2 (5 mol%)", {"number": 5, "unit": "mol%"}, [E[2]], aryl_sub, "inherited_general_condition")
    charge(pk, "equivalents", "KOPiv (1 equiv.)", {"number": 1, "unit": "equiv"}, [E[2]], aryl_sub, "inherited_general_condition")
    r["reaction_scale"] = f("1 (0.2 mmol)", {"reference_molecule_id": mid(aryl_sub), "amount": {"number": 0.2, "unit": "mmol"}}, [E[2]], "inherited_general_condition")
    r["solvent_conditions"]["composition"] = f("in DCE (0.2 M)", [mid("DCE")], [E[2]], "inherited_general_condition")
    r["solvent_conditions"]["concentration"] = f("DCE (0.2 M)", {"number": 0.2, "unit": "M"}, [E[2]], "inherited_general_condition")
    r["solvent_conditions"]["concentration_basis"] = f("1 (0.2 mmol) in DCE (0.2 M)", mid(aryl_sub), [E[2]], "inherited_general_condition")
    reaction_step(r, "reaction", [E[2]], [p["participant_id"] for p in r["participants"]],
                  "Scheme 2 standard conditions", ("100 °C", {"number": 100, "unit": "°C"}),
                  ("4 h", {"number": 4, "unit": "h"}), "without exclusion of air or moisture")
    r["procedure"]["general_procedure_reference"] = f("Scheme 2 standard conditions", "Scheme 2 standard conditions", [E[2]], "inherited_general_condition")
    outcome(r, "isolated_yield", "yield", [product], [E[2]], f"{product}, {yld}%",
            {"number": yld, "unit": "%"}, "main product", "isolated", "isolated yield",
            f"yield of {product}")
    outcome(r, "E_selectivity", "selectivity", [product], [T[3], E[3]],
            "only E-selective products were obtained regardless of the nature of both substrates",
            "E only", "main product", basis=None,
            scope="exocyclic-alkene geometry; E:Z ratio and analytical method not reported")
    set_classification(r, [("C", "C", 1, "formed"), ("C", "O", 1, "formed")], [E[2]],
                       "Partial net classification of the annulation; retained allene/product bond orders and implicit-H transfers were not atom-mapped.")
    unresolved = ["No SI was supplied; addition order, workup, purification and product characterization are unavailable.",
                  "E-only statement has no numerical E:Z ratio or analytical method on the supplied pages.",
                  "Product SMILES/absolute configuration remain unencoded except for retained source E assignment."]
    if product in {"3ia", "3ja", "3ka"}:
        unresolved.append("The prose reports partial starting-material recovery for electron-withdrawing substrates but gives no product-specific recovery percentage.")
    if product in {"3ae", "3ah", "3ai", "3aj"}:
        unresolved.append("The prose reports possible dehydration side products for aryl- or linear-alkyl-tethered α-allenols, without product-specific amounts or measurement basis.")
    r["review"]["unresolved"] = unresolved
    bundle["reactions"].append(r)


for row in scope_rows:
    add_scope_record(*row)


def add_standard_agents(r, evidence, reference, solvent="DCE", include_d2o=False):
    for label, role in [(CAT, "catalyst"), ("KOPiv", "additive"), (solvent, "solvent")]:
        if not any(p["molecule_id"] == mid(label) for p in r["participants"]):
            participant(r, label, role, evidence)
    if include_d2o and not any(p["molecule_id"] == mid("D2O") for p in r["participants"]):
        participant(r, "D2O", "reagent", evidence)
    pc = next(p for p in r["participants"] if p["molecule_id"] == mid(CAT))
    pk = next(p for p in r["participants"] if p["molecule_id"] == mid("KOPiv"))
    charge(pc, "mol_percent", "[Cp*RhCl2]2 (5 mol%)", {"number": 5, "unit": "mol%"}, evidence, reference)
    charge(pk, "equivalents", "KOPiv (1 equiv)", {"number": 1, "unit": "equiv"}, evidence, reference)


# Scheme 3a, deuterium exchange on recovered 1a.
r = rx("scheme3a_1a_D2O_exchange", [E[3], T[3]], categories=["mechanistic_study"],
       location="PDF page 3, Scheme 3a, left arrow", reaction_name="deuterium-labeling experiment",
       transformation="H/D exchange")
p1 = participant(r, "1a", "reactant", [E[3]])
add_standard_agents(r, [E[3]], "1a", include_d2o=True)
participant(r, "1a", "product", [E[3]], pid="p_1a_recovered")
r["solvent_conditions"]["composition"] = f("DCE/D2O", [mid("DCE"), mid("D2O")], [E[3]])
reaction_step(r, "exchange", [E[3]], [p["participant_id"] for p in r["participants"]],
              "deuterium-labeling experiment", ("100 °C", {"number": 100, "unit": "°C"}),
              ("0.5 h", {"number": 0.5, "unit": "h"}))
outcome(r, "recovery", "recovery", ["1a"], [E[3]], "76% recovery", {"number": 76, "unit": "%"},
        product_role="recovered starting material", basis=None, scope="recovery of 1a")
outcome(r, "D_incorporation", "other", ["1a"], [E[3]], "30% D at the H/D positions drawn",
        {"number": 30, "unit": "%", "positions": "two aryl H/D positions as drawn"},
        scope="deuterium incorporation; analytical method not reported on the page")
set_classification(r, [("C", "H", 1, "hydrogen_gain")], [E[3]],
                   "C-D formation is indexed under the C-H taxonomy at isotope level; exact atom mapping and exchanged-H loss are unresolved.")
r["review"]["unresolved"] = ["The page does not report D2O ratio/volume, atmosphere, workup, recovery basis or deuterium-analysis method."]
bundle["reactions"].append(r)

# Scheme 3a, annulation in D2O.
add_molecule("3aa-d-labeled", [E[3]], "deuterium-labelled 3aa",
             stereo={"chirality": "unspecified", "geometry": "specified", "representation": "not_encoded",
                     "double_bonds": [{"feature": "exocyclic alkene", "configuration": "E"}],
                     "assignment": "Positional D incorporation is shown, but a complete isotopomer graph is not encoded."})
r = rx("scheme3a_3aa_D2O_annulation", [E[3], T[3]], categories=["mechanistic_study"],
       location="PDF page 3, Scheme 3a, right arrow", reaction_name="deuterium-labeling annulation",
       transformation="deuterium-labelled [3+3] annulation")
for label, role in [("1a", "reactant"), ("2a", "reactant"), ("3aa-d-labeled", "product")]:
    participant(r, label, role, [E[3]])
add_standard_agents(r, [E[3]], "1a", include_d2o=True)
r["solvent_conditions"]["composition"] = f("DCE/D2O", [mid("DCE"), mid("D2O")], [E[3]])
reaction_step(r, "reaction", [E[3]], [p["participant_id"] for p in r["participants"]],
              "deuterium-labeling annulation", ("100 °C", {"number": 100, "unit": "°C"}),
              ("0.5 h", {"number": 0.5, "unit": "h"}))
outcome(r, "yield", "yield", ["3aa-d-labeled"], [E[3]], "3aa, 58%", {"number": 58, "unit": "%"},
        product_role="deuterium-labelled product", scope="reported yield; method not reported")
outcome(r, "D_incorporation", "other", ["3aa-d-labeled"], [E[3]],
        "~30% D at the aryl H/D site and 81% D at the allylic H/D site as drawn",
        [{"position": "aryl H/D site as drawn", "value": {"number": 30, "unit": "%", "qualifier": "approx"}},
         {"position": "allylic H/D site adjacent to Cy as drawn", "value": {"number": 81, "unit": "%"}}],
        scope="positional deuterium incorporation; analytical method not reported")
set_classification(r, [("C", "C", 1, "formed"), ("C", "O", 1, "formed"),
                       ("C", "H", 1, "hydrogen_gain")], [E[3]],
                   "Annulation plus isotopic C-D formation; isotope positions and full atom mapping remain partial.")
r["review"]["unresolved"] = ["D2O ratio/volume, atmosphere, workup, yield basis and deuterium-analysis method are unreported on supplied pages."]
bundle["reactions"].append(r)

# Scheme 3a, 3ea in CD3OD.
add_molecule("3ea-d2", [E[3]], "deuterated 3ea-d2",
             stereo={"chirality": "unspecified", "geometry": "specified", "representation": "not_encoded",
                     "double_bonds": [{"feature": "exocyclic alkene", "configuration": "E"}],
                     "assignment": "Two H/D positions are drawn; exact isotopomer graph is not encoded."})
r = rx("scheme3a_3ea_CD3OD_exchange", [E[3], T[3]], categories=["mechanistic_study"],
       location="PDF page 3, Scheme 3a, lower arrow", reaction_name="CD3OD deuterium-incorporation experiment",
       transformation="H/D exchange")
for label, role in [("3ea", "reactant"), ("CD3OD", "solvent"), ("3ea-d2", "product")]:
    participant(r, label, role, [E[3]])
r["solvent_conditions"]["composition"] = f("CD3OD", [mid("CD3OD")], [E[3]])
reaction_step(r, "exchange", [E[3]], [p["participant_id"] for p in r["participants"]],
              "treatment of 3ea in CD3OD", ("100 °C", {"number": 100, "unit": "°C"}),
              ("3 h", {"number": 3, "unit": "h"}))
outcome(r, "D_incorporation", "other", ["3ea-d2"], [E[3]], "~30% D at the two H/D positions drawn",
        {"number": 30, "unit": "%", "qualifier": "approx", "positions": "two allylic H/D sites as drawn"},
        scope="deuterium incorporation; product amount and analytical method not reported")
set_classification(r, [("C", "H", 1, "hydrogen_gain")], [E[3]],
                   "C-D formation indexed as C-H; exact exchanged atoms and isotope balance are unresolved.")
r["review"]["unresolved"] = ["CD3OD volume, atmosphere, recovery/product amount, workup and deuterium-analysis method are unreported."]
bundle["reactions"].append(r)

# KIE substrates/products.
add_molecule("1a-d5", [E[3]], "deuterated N-phenoxyacetamide 1a-d5",
             stereo={"chirality": "achiral", "geometry": "not_applicable", "representation": "not_encoded"})
add_molecule("3aa-d4", [E[3]], "deuterated product 3aa-d4",
             stereo={"chirality": "unspecified", "geometry": "specified", "representation": "not_encoded",
                     "double_bonds": [{"feature": "exocyclic alkene", "configuration": "E"}]})
r = rx("scheme3b_competitive_KIE", [E[3], T[3]], categories=["mechanistic_study"],
       location="PDF page 3, Scheme 3b, competitive experiment", reaction_name="competitive kinetic isotope experiment",
       transformation="competitive [3+3] annulation")
for label, role in [("1a", "reactant"), ("1a-d5", "reactant"), ("2a", "reactant"),
                    ("3aa", "product"), ("3aa-d4", "product")]:
    participant(r, label, role, [E[3]])
add_standard_agents(r, [E[3]], "1a")
r["solvent_conditions"]["composition"] = f("DCE", [mid("DCE")], [E[3]])
reaction_step(r, "competition", [E[3]], [p["participant_id"] for p in r["participants"]],
              "competitive KIE experiment", ("100 °C", {"number": 100, "unit": "°C"}),
              ("0.5 h", {"number": 0.5, "unit": "h"}))
outcome(r, "KIE", "selectivity", ["3aa", "3aa-d4"], [E[3]], "KH/KD = 1.6",
        {"ratio": 1.6, "numerator": "kH", "denominator": "kD"}, basis="competitive experiment",
        scope="primary kinetic isotope effect")
set_classification(r, [("C", "C", 1, "formed"), ("C", "O", 1, "formed")], [E[3]],
                   "Products are the protiated/deuterated annulation products; isotopic balance and exact mapping are unresolved.")
r["review"]["unresolved"] = ["Individual substrate amounts, product ratios/yields, kinetic measurement details, workup and atmosphere are unreported."]
bundle["reactions"].append(r)

r = rx("scheme3b_parallel_KIE", [E[3], T[3]], categories=["mechanistic_study"],
       location="PDF page 3, Scheme 3b, two parallel side-by-side reactions", reaction_name="parallel kinetic isotope experiment",
       transformation="parallel [3+3] annulations")
for label, role in [("1a", "reactant"), ("1a-d5", "reactant"), ("2a", "reactant"),
                    ("3aa", "product"), ("3aa-d4", "product")]:
    participant(r, label, role, [E[3]])
add_standard_agents(r, [E[3]], "1a")
r["solvent_conditions"]["composition"] = f("standard conditions (DCE shown for the paired KIE study)", [mid("DCE")], [E[3]], "inherited_general_condition")
reaction_step(r, "parallel_H_run", [E[3]], ["p_1a", "p_2a", "p_3aa"],
              "1a + 2a under standard conditions", time=("0.5 h", {"number": 0.5, "unit": "h"}))
reaction_step(r, "parallel_D_run", [E[3]], ["p_1a_d5", "p_2a", "p_3aa_d4"],
              "1a-d5 + 2a under standard conditions", time=("0.5 h", {"number": 0.5, "unit": "h"}))
outcome(r, "KIE", "selectivity", ["3aa", "3aa-d4"], [E[3]], "KH/KD = 1.3",
        {"ratio": 1.3, "numerator": "kH", "denominator": "kD"}, basis="two parallel, side-by-side reactions",
        scope="primary kinetic isotope effect")
set_classification(r, [("C", "C", 1, "formed"), ("C", "O", 1, "formed")], [E[3]],
                   "Paired product-forming reactions share class-level annulation changes; exact isotope mapping is unresolved.")
r["review"]["unresolved"] = ["The arrow says standard conditions but does not repeat temperature on the two lower arrows; temperature is not silently inherited into the step facts.",
                              "Individual charges, product yields, kinetic calculation details, workup and atmosphere are unreported."]
bundle["reactions"].append(r)

# Scheme 3c/d substrate and product labels.
for label, name in [("4", "modified α-allenol 4"), ("5", "2H-chromen-2-one 5"),
                    ("6", "modified α-allenol 6"), ("7", "diene derivative 7"),
                    ("8", "allene substrate 8 without hydroxyl group"),
                    ("9", "NHTs-bearing allene 9"), ("10", "4-alkylidene chroman-2-amide 10")]:
    add_molecule(label, [E[3]], name, stereo={"chirality": "unresolved", "geometry": "unresolved", "representation": "unresolved"})
add_template_molecule("disordered complex product mixture", [E[3]],
    {"composition": "unknown", "source_result": "disordered / inseparable mixture of complex products"},
    "The source reports an inseparable complex mixture and does not define unique components.",
    {"chirality": "mixture", "geometry": "mixture", "representation": "unresolved",
     "mixture": "inseparable complex mixture; composition unreported"})

for suffix, substrate, product, yld, labels in [
    ("scheme3c_allenol4_to_5", "4", "5", 64, [("C", "C", 1, "formed"), ("C", "O", 1, "formed")]),
    ("scheme3c_allenol6_to_7", "6", "7", 42, [("C", "C", 1, "formed")]),
    ("scheme3d_NHTs_allene9_to_10", "9", "10", 31, [("C", "C", 1, "formed"), ("C", "O", 1, "formed")]),
]:
    r = rx(suffix, [E[3], T[3]], categories=["control", "mechanistic_study"],
           location=f"PDF page 3, Scheme 3 {'c' if '3c' in suffix else 'd'}",
           reaction_name="control experiment", transformation="control transformation")
    for label, role in [("1a", "reactant"), (substrate, "reactant"), (product, "product")]:
        participant(r, label, role, [E[3]])
    add_standard_agents(r, [E[3]], "1a")
    r["solvent_conditions"]["composition"] = f("DCE", [mid("DCE")], [E[3]])
    reaction_step(r, "reaction", [E[3]], [p["participant_id"] for p in r["participants"]],
                  "Scheme 3 control reaction", ("100 °C", {"number": 100, "unit": "°C"}),
                  ("4 h", {"number": 4, "unit": "h"}))
    outcome(r, "yield", "yield", [product], [E[3]], f"{yld}%", {"number": yld, "unit": "%"},
            product_role="reported product", scope="reported yield; method not reported")
    set_classification(r, labels, [E[3]], "Class-level bond formation inferred from the control-experiment drawings; exact mapping is unresolved.")
    r["review"]["unresolved"] = ["Reactant charges, atmosphere, workup, purification and yield basis are not reported on the supplied page.",
                                  "Exact product SMILES/stereochemistry and atom mapping remain unresolved."]
    bundle["reactions"].append(r)

r = rx("scheme3d_noOH_allene8_mixture", [E[3], T[3]], categories=["control", "mechanistic_study"],
       location="PDF page 3, Scheme 3d, upper arrow", reaction_name="hydroxyl-removal control")
for label, role in [("1a", "reactant"), ("8", "reactant"), ("disordered complex product mixture", "product")]:
    participant(r, label, role, [E[3]])
add_standard_agents(r, [E[3]], "1a")
r["solvent_conditions"]["composition"] = f("DCE", [mid("DCE")], [E[3]])
reaction_step(r, "reaction", [E[3]], [p["participant_id"] for p in r["participants"]],
              "hydroxyl-removal coordination control", ("100 °C", {"number": 100, "unit": "°C"}),
              ("4 h", {"number": 4, "unit": "h"}))
outcome(r, "complex_mixture", "decomposition", ["disordered complex product mixture"], [E[3]],
        "disordered; an inseparable mixture of complex products", "inseparable complex product mixture",
        product_role="unresolved product mixture", scope="qualitative outcome; components and measurement method unreported")
no_classification(r, [E[3]], "unresolved", "Mixture components are undefined, so observed net bond changes cannot be assigned.")
r["review"]["unresolved"] = ["Mixture identities, quantities, conversion/recovery, measurement method, atmosphere and workup are unreported."]
bundle["reactions"].append(r)

# Scheme 4 proposed cycle, not executed experiments.
r = rx("scheme4_proposed_mechanism", [E[3], T[3]], "proposed_mechanism", None,
       "PDF page 3, Scheme 4", "proposed catalytic cycle")
r["notes"] = [f("Proposed intermediates A-F and C′; reversible C-H cleavage, allene insertion, O-N oxidative addition, beta-H/reductive elimination, protonolysis and cyclization.",
                  "Author-proposed catalytic cycle; intermediates are not treated as isolated products.", [E[3], T[3]])]
no_classification(r, [E[3]], "not_applicable", "Proposed elementary steps are not separate observed experiments.")
r["review"]["unresolved"] = ["Intermediate structures and elementary-step mappings are mechanistic proposals, not source-verified isolated species."]
bundle["reactions"].append(r)

# Page 4 scale-up and derivatizations.
for label, name in [("11", "reductive ortho-alkenyl phenol 11"),
                    ("12", "ring-opening product 12 (R = Me)"),
                    ("13", "ring-opening product 13 (R = C≡CH)"),
                    ("MeMgBr", "methylmagnesium bromide"),
                    ("HC≡CMgBr", "ethynylmagnesium bromide")]:
    add_molecule(label, [E[4], T[3], T[4]], name, stereo={
        "chirality": "unspecified" if label in ("12", "13") else "unresolved",
        "geometry": "unspecified" if label in ("11", "12", "13") else "not_applicable",
        "representation": "not_encoded" if label in ("11", "12", "13") else "unresolved",
    })

r = rx("eq1_scaleup_3aa", [E[4], T[3]], categories=["scale_up"],
       location="PDF page 4, equation (1); scale-up discussion begins on page 3",
       reaction_name="scaled-up [3+3] annulation", transformation="redox-neutral [3+3] annulation")
r["general_scheme_ids"] = [GENERAL]
for label, role in [("1a", "reactant"), ("2a", "reactant"), (CAT, "catalyst"),
                    ("KOPiv", "additive"), ("DCE", "solvent"), ("3aa", "product")]:
    participant(r, label, role, [E[4], T[3]])
pc = next(p for p in r["participants"] if p["molecule_id"] == mid(CAT))
pk = next(p for p in r["participants"] if p["molecule_id"] == mid("KOPiv"))
charge(pc, "mol_percent", "5 mol%", {"number": 5, "unit": "mol%"}, [E[4]], "1a")
charge(pk, "equivalents", "1 equiv", {"number": 1, "unit": "equiv"}, [E[4]], "1a")
r["reaction_scale"] = inferred({"reference_molecule_id": mid("1a"), "amount": {"number": 10, "unit": "mmol"}},
    [E[4], T[3]], "The source states 'on 10 mmol scale'; 1a is selected as the scale reference by continuity with Scheme 2, but the arrow does not explicitly name the denominator.",
    raw="68% on 10 mmol scale")
r["solvent_conditions"]["composition"] = f("DCE", [mid("DCE")], [E[4]])
reaction_step(r, "reaction", [E[4]], [p["participant_id"] for p in r["participants"]],
              "scale-up reaction", ("100 °C", {"number": 100, "unit": "°C"}),
              ("4 h", {"number": 4, "unit": "h"}))
outcome(r, "yield", "yield", ["3aa"], [E[4]], "68% on 10 mmol scale", {"number": 68, "unit": "%"},
        product_role="main product", scope="scale-up yield; method not explicitly stated")
outcome(r, "product_mass", "product_amount", ["3aa"], [E[4]], "3aa, 1.755 g", {"number": 1.755, "unit": "g"},
        product_role="main product", scope="isolated/collected mass basis not explicitly stated")
set_classification(r, [("C", "C", 1, "formed"), ("C", "O", 1, "formed")], [E[4]],
                   "Class-level scale-up annulation assignment; exact atom mapping remains unresolved.")
r["review"]["unresolved"] = ["The 10 mmol scale denominator is inferred as 1a; exact substrate charges, DCE volume, atmosphere, workup, purification and yield basis are not reported."]
bundle["reactions"].append(r)

r = rx("eq2_3aa_to_11", [E[4], T[3]], categories=["derivatization"],
       location="PDF page 4, equation (2), left arrow", reaction_name="LiAlH4 reduction/ring opening",
       transformation="reductive ring opening")
for label, role in [("3aa", "reactant"), ("LiAlH4", "reagent"), ("MeOH", "reagent"),
                    ("THF", "solvent"), ("11", "product")]:
    if label == "THF":
        add_molecule("THF", [E[4]], "tetrahydrofuran", "C1CCOC1",
                     {"chirality": "achiral", "geometry": "not_applicable", "representation": "not_applicable"})
    participant(r, label, role, [E[4]])
pla = next(p for p in r["participants"] if p["molecule_id"] == mid("LiAlH4"))
pmeoh = next(p for p in r["participants"] if p["molecule_id"] == mid("MeOH"))
charge(pla, "equivalents", "LiAlH4 (3 equiv)", {"number": 3, "unit": "equiv"}, [E[4]], "3aa")
charge(pmeoh, "equivalents", "MeOH (4 equiv)", {"number": 4, "unit": "equiv"}, [E[4]], "3aa")
r["solvent_conditions"]["composition"] = f("THF", [mid("THF")], [E[4]])
reaction_step(r, "reaction", [E[4]], [p["participant_id"] for p in r["participants"]],
              "LiAlH4/MeOH reduction of 3aa", ("rt", "rt"), ("3 h", {"number": 3, "unit": "h"}))
outcome(r, "yield", "yield", ["11"], [E[4]], "11, 61%", {"number": 61, "unit": "%"},
        product_role="reductive ortho-alkenyl phenol", scope="reported yield; method not stated")
set_classification(r, [("C", "H", 1, "hydrogen_gain"), ("O", "H", 1, "hydrogen_gain")], [E[4]],
                   "Partial net reduction/protonation classification; ring C-O cleavage and exact hydrogen count/source are unresolved.")
r["review"]["unresolved"] = ["3aa charge, THF volume, atmosphere, workup, purification, yield basis and full hydrogen mapping are unreported."]
bundle["reactions"].append(r)

for product, reagent, yld, rval in [("12", "MeMgBr", 51, "Me"), ("13", "HC≡CMgBr", 54, "C≡CH")]:
    r = rx(f"eq2_3aa_to_{product}", [E[4], T[3], T[4]], categories=["derivatization"],
           location=f"PDF page 4, equation (2), R = {rval}", reaction_name="Grignard ring opening",
           transformation="nucleophilic ring opening")
    for label, role in [("3aa", "reactant"), (reagent, "reagent"), ("THF", "solvent"), (product, "product")]:
        participant(r, label, role, [E[4]])
    pr = next(p for p in r["participants"] if p["molecule_id"] == mid(reagent))
    charge(pr, "equivalents", "RMgBr (2 equiv)", {"number": 2, "unit": "equiv"}, [E[4]], "3aa")
    r["solvent_conditions"]["composition"] = f("THF", [mid("THF")], [E[4]])
    reaction_step(r, "reaction", [E[4]], [p["participant_id"] for p in r["participants"]],
                  f"3aa treated with RMgBr, R = {rval}", ("rt", "rt"),
                  ("12 h", {"number": 12, "unit": "h"}), "under N2")
    outcome(r, "yield", "yield", [product], [E[4]], f"R = {rval}: {product}, {yld}%",
            {"number": yld, "unit": "%"}, product_role="ring-opening product",
            scope="reported yield; method not stated")
    set_classification(r, [("C", "C", 1, "formed"), ("O", "H", 1, "hydrogen_gain")], [E[4]],
                       "Partial Grignard-addition/ring-opening classification; Mg-C cleavage, proton source and full atom mapping remain unresolved.")
    r["review"]["unresolved"] = ["3aa charge, THF volume, workup, purification, yield basis, proton source and product stereochemistry are unreported."]
    bundle["reactions"].append(r)


bundle["source_context"] = [
    {"context_id": f"{DOC}:context:paper_metadata", "document_id": DOC,
     "evidence_ids": [E[1], T[1]], "kind": "paper_metadata",
     "value": {"title": "Rh(III)-Catalyzed and synergistic dual directing group-enabled redox-neutral [3+3] annulation of N-phenoxyacetamides with α-allenols",
               "doi": "10.1039/d1cc03206b", "journal": "Chem. Commun.", "year": 2021, "volume": 57, "pages": "9284-9287"}},
    {"context_id": f"{DOC}:context:si_unavailable", "document_id": DOC,
     "evidence_ids": [E[1], T[1]], "kind": "source_limitation",
     "value": "The article advertises ESI, experimental procedures, characterization, spectra and CCDC 2084815, but no SI/CIF was supplied or associated in this session."},
    {"context_id": f"{DOC}:context:optimization_summary", "document_id": DOC,
     "evidence_ids": [E[2], T[2]], "kind": "discussion_not_row_resolved",
     "value": "Solvent/catalyst/N-substituent screening is summarized in prose, but individual Table S1 rows and numerical outcomes are unavailable; only the explicit CsOAc entry-1 statement and ineffective N-phenoxytosylamide statement were made concrete records."},
    {"context_id": f"{DOC}:context:scope_limitations", "document_id": DOC,
     "evidence_ids": [E[2], T[2], E[3], T[3]], "kind": "discussion",
     "value": "Electron-withdrawing substrates partly recover starting material; aryl/linear-alkyl allenols may dehydrate; ester-, alkenyl-, heterocycle-substituted and trisubstituted α-allenols were incompatible. Product-specific quantities/identities are not supplied."},
    {"context_id": f"{DOC}:context:xray_3ai", "document_id": DOC,
     "evidence_ids": [E[2], T[2]], "kind": "structure_assignment",
     "value": "The structure of 3ai was confirmed by X-ray crystallography, CCDC 2084815; no CIF was supplied."},
    {"context_id": f"{DOC}:context:no_standalone_spectral_table", "document_id": DOC,
     "evidence_ids": [E[1], E[2], E[3], E[4]], "kind": "negative_inventory",
     "value": "No standalone mass-spectrum peak table, molecular-weight table or fragment-relative-abundance block occurs on inspected pages 1-4."},
]

# Save submission and raw extraction inventory.
st.write(RUN / "submitted.json", bundle)

page1_ids = [f"{DOC}:scheme1a_prior_work", f"{DOC}:scheme1c_general"]
scope_ids = [f"{DOC}:scheme2_{p}" for p, _, _, _ in scope_rows]
page2_ids = [f"{DOC}:tableS1_entry1_CsOAc", f"{DOC}:N_phenoxytosylamide_ineffective"] + scope_ids
page3_mech = [
    f"{DOC}:scheme3a_1a_D2O_exchange", f"{DOC}:scheme3a_3aa_D2O_annulation",
    f"{DOC}:scheme3a_3ea_CD3OD_exchange", f"{DOC}:scheme3b_competitive_KIE",
    f"{DOC}:scheme3b_parallel_KIE", f"{DOC}:scheme3c_allenol4_to_5",
    f"{DOC}:scheme3c_allenol6_to_7", f"{DOC}:scheme3d_noOH_allene8_mixture",
    f"{DOC}:scheme3d_NHTs_allene9_to_10", f"{DOC}:scheme4_proposed_mechanism",
]
page4_ids = [f"{DOC}:eq1_scaleup_3aa", f"{DOC}:eq2_3aa_to_11",
             f"{DOC}:eq2_3aa_to_12", f"{DOC}:eq2_3aa_to_13"]

inventory = {
    "document_id": DOC,
    "model_output_status": "needs_review",
    "source_pages_inspected": [1, 2, 3, 4],
    "deduplication": [
        "Scheme 2 3aa and the general optimized example share the same standard-condition experiment; the 73% Scheme 2 row is retained once.",
        "The prose-only CsOAc/Table S1 entry 1 differs by additive and is retained separately.",
        "Page-3 prose introducing scale-up/derivatization is linked to the page-4 arrows, not duplicated.",
    ],
    "pages": {
        "1": {"regions": ["title/abstract", "Scheme 1a prior work", "Scheme 1b natural products", "Scheme 1c current generic transformation", "ESI notice"], "record_ids": page1_ids},
        "2": {"regions": ["optimization prose", "Scheme 2 generic header/caption", "Scheme 2a 18 N-phenoxyacetamide scope products", "Scheme 2b 9 α-allenol scope products"], "record_ids": page2_ids},
        "3": {"regions": ["scope selectivity/limitations prose", "Scheme 3a deuterium experiments", "Scheme 3b KIE", "Scheme 3c/d controls", "Scheme 4 proposed mechanism", "scale-up/derivatization lead-in"], "record_ids": scope_ids + page3_mech + [f"{DOC}:eq1_scaleup_3aa", f"{DOC}:eq2_3aa_to_11", f"{DOC}:eq2_3aa_to_12", f"{DOC}:eq2_3aa_to_13"]},
        "4": {"regions": ["equation 1 scale-up", "equation 2 LiAlH4 derivatization", "equation 2 two RMgBr variants", "conclusion/references"], "record_ids": page4_ids},
    },
    "specific_experiment_count": sum(r["record_type"] == "specific_experiment" for r in bundle["reactions"]),
    "reaction_count_total": len(bundle["reactions"]),
    "molecule_count": len(bundle["molecules"]),
    "limitations": [
        "No SI/CIF was provided; all records remain partial and needs_review.",
        "Many compact scope structures are source-labelled but not assigned concrete SMILES.",
        "Scope E selectivity has no numerical E:Z ratio or analytical method.",
        "Mechanistic and equation-2 yield methods, workups and purifications are unreported on the supplied pages.",
    ],
}
st.write(RUN / "raw_model/page_extraction.json", inventory)

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "lower-right Scheme 1a", "finding": "grouped cited prior-work reaction modes; retained as literature_example, not current experiments", "evidence_ids": [E[1]]},
            {"region": "lower-right Scheme 1b", "finding": "representative natural products/bioactive molecules; non-experimental context", "evidence_ids": [E[1]]},
            {"region": "lower-right Scheme 1c", "finding": "current generic [3+3] transformation with catalyst/additive/DCE/100 °C/4 h; one generic_scheme record", "evidence_ids": [E[1]]},
            {"region": "footnote/metadata", "finding": "ESI and CCDC are advertised but were not supplied or associated", "evidence_ids": [E[1], T[1]]},
        ], "record_ids": page1_ids, "unresolved": ["Prior-work families are not expanded to individual cited experiments.", "No SI/CIF available."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "upper-left optimization prose", "finding": "explicit CsOAc/Table S1 entry-1 experiment and ineffective N-phenoxytosylamide test; other missing Table S1 rows retained as context only", "evidence_ids": [E[2], T[2]]},
            {"region": "bottom Scheme 2a", "finding": "18 independent N-phenoxyacetamide scope experiments 3aa-3ra with printed isolated yields", "evidence_ids": [E[2], T[2]]},
            {"region": "bottom Scheme 2b", "finding": "9 independent α-allenol scope experiments 3ab-3aj with printed isolated yields; 3ai X-ray annotation retained as context", "evidence_ids": [E[2], T[2]]},
            {"region": "Scheme 2 caption/footnote", "finding": "common charges, catalyst loading, KOPiv, DCE concentration, temperature, time, open-air/moisture statement and isolated-yield basis inherited by each scope row", "evidence_ids": [E[2], T[2]]},
        ], "record_ids": page2_ids, "unresolved": ["Individual SI optimization rows and procedures unavailable.", "Compact scope structures are not all encoded as SMILES."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "upper-left scope discussion", "finding": "E-only selectivity applies to the Scheme 2 products; dehydration/recovery/incompatible-substrate statements lack product-specific quantities", "evidence_ids": [E[3], T[3]]},
            {"region": "lower-left Scheme 3a", "finding": "three deuterium-labeling/exchange experiments with 76% recovery, 58% yield and positional D incorporation", "evidence_ids": [E[3], T[3]]},
            {"region": "lower-left Scheme 3b", "finding": "competitive KIE 1.6 and paired parallel KIE 1.3 retained as two mechanistic-study records", "evidence_ids": [E[3], T[3]]},
            {"region": "lower-left Scheme 3c/d", "finding": "two product-forming controls (5, 7), one unresolved complex-mixture control, and one NHTs-allene product 10", "evidence_ids": [E[3], T[3]]},
            {"region": "upper-right Scheme 4", "finding": "author-proposed catalytic cycle; one proposed_mechanism record, not executed syntheses", "evidence_ids": [E[3], T[3]]},
            {"region": "lower-right lead-in", "finding": "scale-up and derivatization prose linked to page-4 equations", "evidence_ids": [E[3], T[3]]},
        ], "record_ids": scope_ids + page3_mech + page4_ids, "unresolved": ["Mechanistic-study analytical methods, charges, atmospheres and workups are mostly unreported.", "Mixture composition and exact isotope atom mapping remain unresolved."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "upper-left equation (1)", "finding": "10 mmol scale 3aa synthesis, 68%, 1.755 g", "evidence_ids": [E[4], T[3]]},
            {"region": "upper-left equation (2), left arrow", "finding": "3aa to 11 with LiAlH4 3 equiv/MeOH 4 equiv in THF at rt for 3 h, 61%", "evidence_ids": [E[4], T[3]]},
            {"region": "upper-left equation (2), right arrow/table", "finding": "two separate RMgBr ring-opening variants: 12 (R=Me, 51%) and 13 (R=C≡CH, 54%), 2 equiv, THF, rt, 12 h, N2", "evidence_ids": [E[4], T[4]]},
            {"region": "remaining conclusion/references", "finding": "non-experimental summary and bibliography; no additional current reaction arrow", "evidence_ids": [E[4], T[4]]},
        ], "record_ids": page4_ids, "unresolved": ["Yield methods, substrate charges, solvent volumes, workups, purifications and product stereochemistry are unreported.", "The 10 mmol scale denominator is inferred as 1a and explicitly flagged."]},
    ],
}
st.write(RUN / "coverage/pages_0001_0004.json", coverage)

stereo_checks = {}
for molecule in bundle["molecules"]:
    smiles = molecule["smiles"]["normalized_value"]
    if smiles:
        stereo_checks[molecule["molecule_id"]] = st.stereo_inspect(smiles)
st.write(RUN / "raw_model/stereo_checks.json", {
    "source_verified": False,
    "note": "Standard stereo diagnostics for every molecule carrying a candidate SMILES; graph diagnostics do not certify transcription or source identity.",
    "results": stereo_checks,
})

source_mentions = st.source_mentions(json.loads((RUN / "dataset.json").read_text()),
    ["1a", "2a", "3aa", "3ea", "N-phenoxytosylamide", "LiAlH4", "Grignard reagents"],
    context=120)
st.write(RUN / "raw_model/source_mentions.json", source_mentions)

print(json.dumps({"reactions": len(bundle["reactions"]),
                  "specific_experiments": inventory["specific_experiment_count"],
                  "molecules": len(bundle["molecules"]),
                  "source_context": len(bundle["source_context"])}, ensure_ascii=False))
