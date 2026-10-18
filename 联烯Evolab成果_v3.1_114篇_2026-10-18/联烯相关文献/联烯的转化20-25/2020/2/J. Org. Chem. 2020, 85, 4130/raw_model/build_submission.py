#!/usr/bin/env python3
"""Build the page-1--6 v3.1 submission from the visually transcribed source regions.

This is a source-extraction artifact, not a chemistry-certification program.  It
uses only the prepared evidence registry and leaves unverified structures empty.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))

from standard_records import (  # noqa: E402
    bond_label,
    fact,
    new_molecule,
    new_outcome,
    new_participant,
    new_reaction,
    new_step,
    stereo_inspect,
)

DOC = "doc_5909d9a1e9253598"
OUT = ROOT / "runs/current"
TITLE = (
    "Lewis Base-Switched [3 + 3] and [4 + 2] Annulation Reactions of "
    "δ-Acetoxy Allenoates with Cyclic N-Sulfonyl Imines: Divergent Synthesis "
    "of Functionalized α-Pyridyl Acetates and Teraryl Scaffolds"
)


def ev(page: int) -> list[str]:
    return [f"{DOC}:page{page}:image"]


def F(raw, normalized=None, evidence=(), method="original_extraction", status=None,
      derivation=None, alternatives=()):
    if normalized is None and raw is not None:
        normalized = raw
    return fact(raw=raw, normalized=normalized, evidence=evidence, method=method,
                status=status, derivation=derivation, alternatives=alternatives)


def missing(evidence=(), status="unreported"):
    return fact(evidence=evidence, status=status)


def q(number, unit, raw=None, qualifier=None):
    value = {"number": number, "unit": unit}
    if raw is not None:
        value["raw_text"] = raw
    if qualifier:
        value["qualifier"] = qualifier
    return value


def set_fact(obj, key, raw, normalized=None, evidence=(), method="original_extraction",
             status=None, derivation=None, alternatives=()):
    obj[key] = F(raw, normalized, evidence, method, status, derivation, alternatives)


bundle = json.loads((OUT / "submitted.json").read_text())
bundle["molecules"] = []
bundle["reactions"] = []
bundle["source_context"] = []
bundle["documents"][0]["title"] = TITLE

molecules: dict[str, dict] = {}


def add_occurrence(m, label, evidence):
    existing = {(o["document_id"], o["compound_label"]["normalized_value"], tuple(o["evidence_ids"]))
                for o in m["occurrences"]}
    key = (DOC, label, tuple(evidence))
    if key not in existing:
        m["occurrences"].append({
            "document_id": DOC,
            "compound_label": F(label, label, evidence),
            "evidence_ids": list(evidence),
        })


def mol(label: str, name: str | None = None, pages=(2,), kind="ordinary", smiles=None):
    mid = f"{DOC}:mol:{label.replace('/', '_').replace(' ', '_')}"
    evidence = sorted({e for p in pages for e in ev(p)})
    if mid in molecules:
        for p in pages:
            add_occurrence(molecules[mid], label, ev(p))
        if name and molecules[mid]["names"]["original"]["normalized_value"] is None:
            molecules[mid]["names"]["original"] = F(name, name, evidence)
        return mid
    m = new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, evidence)
    m["names"]["original"] = F(name or label, name or label, evidence)
    m["names"]["abbreviations"] = F(label, [label], evidence)
    m["identity_note"] = F(
        {"source_label": label, "visual_review_scope": "PDF pages 1-6"},
        "Source label registered; full drawn structure was not independently structure-verified.",
        evidence,
    )
    stereo = m["stereochemistry"]
    if kind == "allene":
        stereo["axial"] = F(
            "allene drawn without R_a/S_a or P/M assignment",
            [{"type": "allene_axis", "assignment": "unspecified"}], evidence)
        stereo["chirality_status"] = F("no axial/absolute descriptor shown", "unspecified", evidence)
        stereo["geometry_status"] = F(None, "not_applicable", evidence, method="inference",
                                      derivation="Allene axial stereochemistry is tracked under axial, not forced into E/Z.")
        stereo["representation_status"] = F("source drawing not encoded as a verified SMILES",
                                             "not_encoded", evidence)
    elif kind == "diastereomeric_allene":
        ratio = "~1:0.9" if label == "2d" else "1:0.9"
        stereo["mixture"] = F(f"diastereomeric mixture ({ratio})",
                               {"type": "diastereomeric_mixture", "ratio": ratio,
                                "isomer_assignments": "unreported"}, evidence)
        stereo["axial"] = F("allene; individual axial configurations not assigned",
                             [{"type": "allene_axis", "assignment": "unassigned"}], evidence)
        stereo["chirality_status"] = F("diastereomeric mixture", "mixture", evidence)
        stereo["geometry_status"] = F(None, "not_applicable", evidence, method="inference",
                                      derivation="Reported isomerism is retained as a diastereomeric mixture, not converted to E/Z.")
        stereo["representation_status"] = F("mixture not encoded as a single stereoisomeric SMILES",
                                             "not_encoded", evidence)
    else:
        for field_name, value in (("chirality_status", "unresolved"),
                                  ("geometry_status", "unresolved"),
                                  ("representation_status", "unresolved")):
            stereo[field_name] = F(
                None, value, evidence, method="inference",
                derivation=("The source label/drawing was inspected, but no complete independently verified "
                            "structure representation is asserted; absence of stereo labels is not treated as achiral."))
    if smiles:
        diag = stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = F(label, diag["canonical_isomeric_smiles"], evidence,
                         method="structure_derivation",
                         derivation="Expanded from the unambiguous source abbreviation/name and checked with standard stereo.")
        stereo["isomeric_smiles"] = F(label, diag["canonical_isomeric_smiles"], evidence,
                                       method="structure_derivation",
                                       derivation="Canonical isomeric SMILES from standard stereo graph diagnostics.")
        stereo["identity_key"] = F(label, diag["stereo_graph_key"], evidence,
                                    method="structure_derivation",
                                    derivation="SHA-256 stereo graph key returned by standard stereo.")
        stereo["atom_centers"] = F(label, diag["atom_centers"], evidence,
                                    method="structure_derivation",
                                    derivation="Computed from the candidate graph; source identity remains needs-review.")
        stereo["double_bonds"] = F(label, diag["double_bonds"], evidence,
                                    method="structure_derivation",
                                    derivation="Computed from the candidate graph; source identity remains needs-review.")
        stereo["chirality_status"] = F(None, "achiral", evidence, method="structure_derivation",
                                        derivation="No potential stereocenters in the checked concrete graph.")
        stereo["geometry_status"] = F(None, "not_applicable", evidence, method="structure_derivation",
                                       derivation="No stereogenic E/Z double bond in the checked concrete graph.")
        stereo["representation_status"] = F(None, "encoded", evidence, method="structure_derivation",
                                             derivation="The checked concrete graph contains no unencoded stereo elements.")
    molecules[mid] = m
    bundle["molecules"].append(m)
    return mid


# Source-labelled reaction families and products.
for label in ["1a", "1b", "1c", "1d", "1e"]:
    mol(label, pages=(1, 2, 3, 4, 5, 6))
for label in ["2a", "2b", "2c", "2f"]:
    mol(label, pages=(1, 2, 3, 4, 5), kind="allene")
mol("2d", "Ethyl 5-acetoxy-5-(4-nitrophenyl)penta-2,3-dienoate", pages=(3, 4, 6),
    kind="diastereomeric_allene")
mol("2e", "tert-Butyl 3-(1-acetoxy-5-ethoxy-5-oxopenta-2,3-dien-1-yl)-1H-indole-1-carboxylate",
    pages=(3, 6), kind="diastereomeric_allene")
for label in ["6a", "6b"]:
    mol(label, pages=(3, 4))

scope3 = {
    "3aa": ("1a", "2a", 79), "3ab": ("1a", "2b", 73), "3ac": ("1a", "2c", 80),
    "3ad": ("1a", "2d", 71), "3ae": ("1a", "2e", 67), "3af": ("1a", "2f", 69),
    "3ba": ("1b", "2a", 81), "3bb": ("1b", "2b", 75), "3bc": ("1b", "2c", 83),
    "3bd": ("1b", "2d", 72), "3ca": ("1c", "2a", 74), "3cb": ("1c", "2b", 67),
    "3cc": ("1c", "2c", 70), "3da": ("1d", "2a", 75), "3db": ("1d", "2b", 72),
    "3dc": ("1d", "2c", 76), "3ea": ("1e", "2a", 69),
}
product_names = {
    "3aa": "Ethyl 2-(6-(2-hydroxyphenyl)-4-phenylpyridin-2-yl)acetate",
    "3ab": "Ethyl 2-(6-(2-hydroxyphenyl)-4-(4-methoxyphenyl)pyridin-2-yl)acetate",
    "3ac": "Ethyl 2-(4-(4-chlorophenyl)-6-(2-hydroxyphenyl)pyridin-2-yl)acetate",
    "3ad": "Ethyl 2-(6-(2-hydroxyphenyl)-4-(4-nitrophenyl)pyridin-2-yl)acetate",
    "3ae": "tert-Butyl 3-(2-(2-ethoxy-2-oxoethyl)-6-(2-hydroxyphenyl)pyridin-4-yl)-1H-indole-1-carboxylate",
}
for label in scope3:
    mol(label, product_names.get(label), pages=(3, 6) if label in product_names else (3,))

scope4 = {"4aa": ("1a", "2a", 73), "4ba": ("1b", "2a", 76),
          "4ca": ("1c", "2a", 71), "4af": ("1a", "2f", 65)}
scope5 = {
    "5aa": ("1a", "2a", 81), "5ab": ("1a", "2b", 75), "5ac": ("1a", "2c", 83),
    "5ad": ("1a", "2d", 85), "5ae": ("1a", "2e", 65), "5af": ("1a", "2f", 68),
    "5ba": ("1b", "2a", 74), "5bb": ("1b", "2b", 70), "5bc": ("1b", "2c", 75),
    "5bd": ("1b", "2d", 77), "5ca": ("1c", "2a", 73), "5cb": ("1c", "2b", 69),
    "5cc": ("1c", "2c", 73), "5cd": ("1c", "2d", 76), "5da": ("1d", "2a", 72),
}
for label in list(scope4) + list(scope5):
    mol(label, pages=(3, 4))

scope7 = {
    "7aa": ("6a", "2a", 71), "7ab": ("6a", "2b", 69), "7ac": ("6a", "2c", 73),
    "7ba": ("6b", "2a", 69), "7bb": ("6b", "2b", 64), "7bc": ("6b", "2c", 70),
    "7ae": ("6a", "2e", 57), "7af": ("6a", "2f", 62),
}
for label in scope7:
    mol(label, pages=(3, 4))
mol("8bb", pages=(3, 4))
scope9 = {"9aa": ("6a", "2a", 76), "9ab": ("6a", "2b", 71),
          "9ac": ("6a", "2c", 79), "9ad": ("6a", "2d", 82),
          "9af": ("6a", "2f", 70), "9ba": ("6b", "2a", 69)}
for label in scope9:
    mol(label, pages=(4,))

# Other source-mentioned chemicals.  Only two unambiguous abbreviations are
# expanded to concrete graphs; all others remain structure-unresolved.
chem = {}
def reagent(key, name=None, pages=(2, 3, 4, 6), smiles=None):
    chem[key] = mol(key, name or key, pages=pages, smiles=smiles)
    return chem[key]

reagent("toluene", "toluene (PhCH3)", pages=(2, 3, 4, 5, 6), smiles="Cc1ccccc1")
reagent("PPh3", "triphenylphosphine", pages=(1, 2, 3, 4, 5),
        smiles="P(c1ccccc1)(c1ccccc1)c1ccccc1")
for key, name, pages in [
    ("DBU", "1,8-diazabicyclo[5.4.0]undec-7-ene", (1, 2, 3, 4, 5, 6)),
    ("DABCO", "1,4-diazabicyclo[2.2.2]octane", (2,)),
    ("Na2CO3", "sodium carbonate", (1, 2, 3, 4, 5, 6)),
    ("K2CO3", "potassium carbonate", (2,)), ("t-BuOK", "potassium tert-butoxide", (2,)),
    ("Et3N", "triethylamine", (6,)), ("acetyl_chloride", "acetyl chloride", (6,)),
    ("DCM", "dry dichloromethane", (6,)), ("water", "water", (6,)),
    ("EtOAc", "ethyl acetate", (6,)), ("brine", "brine", (6,)),
    ("Na2SO4", "anhydrous sodium sulfate", (6,)), ("silica_gel", "silica gel", (6,)),
    ("hexane", "hexane", (6,)),
    ("2d_alcohol", "ethyl 5-hydroxy-5-(4-nitrophenyl)penta-2,3-dienoate", (6,)),
    ("2e_alcohol", "tert-butyl 3-(5-ethoxy-1-hydroxy-5-oxopenta-2,3-dien-1-yl)-1H-indole-1-carboxylate", (6,)),
]:
    reagent(key, name, pages=pages)


def add_charge(p, key, raw, normalized, evidence, method="original_extraction"):
    p["charges"][key] = F(raw, normalized, evidence, method=method)


def participant(r, label, role, evidence, amount=None, mass=None, volume=None,
                equivalents=None, mol_percent=None, reference=None, stage="reaction"):
    mid = chem.get(label) or mol(label, pages=tuple(int(x.split("page")[1].split(":")[0]) for x in evidence))
    pid = f"p{len(r['participants'])+1}_{label.replace('/', '_')}"
    p = new_participant(pid, mid, role, evidence, stage=stage)
    if amount is not None:
        add_charge(p, "amount", f"{amount} mmol", q(amount, "mmol"), evidence)
    if mass is not None:
        add_charge(p, "mass", f"{mass} mg", q(mass, "mg"), evidence)
    if volume is not None:
        add_charge(p, "volume", f"{volume} mL", q(volume, "mL"), evidence)
    if equivalents is not None:
        add_charge(p, "equivalents", f"{equivalents} equiv", q(equivalents, "equiv"), evidence)
    if mol_percent is not None:
        add_charge(p, "mol_percent", f"{mol_percent} mol %", q(mol_percent, "mol%"), evidence)
    if equivalents is not None or mol_percent is not None:
        if reference is None:
            raise ValueError("equivalents/mol% requires reference")
        p["equivalence_reference"] = F(reference, reference, evidence)
    r["participants"].append(p)
    return p


def add_step(r, step_id, evidence, operation=None, temperature=None, time=None,
             stage="reaction", participant_ids=None, method="original_extraction"):
    s = new_step(step_id, evidence, stage=stage)
    s["participant_ids"] = participant_ids or [p["participant_id"] for p in r["participants"] if p["stage"] == stage]
    if operation is not None:
        s["operation"] = F(operation, operation, evidence, method=method)
    if temperature is not None:
        raw = temperature if isinstance(temperature, str) else f"{temperature} °C"
        norm = temperature if isinstance(temperature, str) else q(temperature, "°C")
        s["temperature"] = F(raw, norm, evidence, method=method)
    if time is not None:
        raw = time if isinstance(time, str) else f"{time} h"
        norm = time if isinstance(time, str) else q(time, "h")
        s["time"] = F(raw, norm, evidence, method=method)
    r["steps"].append(s)
    return s


def outcome(r, kind, target_label, value, evidence, basis=None, yield_type=None,
            product_role="product", scope=None, raw=None):
    mid = chem.get(target_label) or mol(target_label, pages=tuple(int(x.split("page")[1].split(":")[0]) for x in evidence))
    o = new_outcome(f"o{len(r['outcomes'])+1}_{kind}_{target_label}", kind, [mid], evidence)
    o["value"] = F(raw if raw is not None else value, value, evidence)
    o["product_role"] = F(product_role, product_role, evidence)
    if basis is not None:
        o["measurement_basis"] = F(basis, basis, evidence)
    if yield_type is not None:
        o["yield_type"] = F(yield_type, yield_type, evidence)
    if scope is not None:
        o["scope"] = F(scope, scope, evidence)
    r["outcomes"].append(o)
    return o


def classify(r, pattern, evidence):
    changes = []
    def change(a, b, before, after, kind, hydrogen_count=None):
        value = {"change_id": f"b{len(changes)+1}",
                 "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
                 "before_order": before, "after_order": after, "change_type": kind}
        if hydrogen_count is not None:
            value["hydrogen_count"] = hydrogen_count
        changes.append(value)
    if pattern == "annulation33":
        change("C", "C", 0, 1, "formed")
        change("C", "N", 0, 1, "formed")
        change("C", "C", 1, 2, "order_increased")
    elif pattern == "annulation42":
        change("C", "C", 0, 1, "formed")
        change("C", "C", 0, 1, "formed")
        change("C", "C", 1, 2, "order_increased")
        change("C", "N", 1, 0, "broken")
    elif pattern == "acetylation":
        change("C", "O", 0, 1, "formed")
        change("O", "H", 1, 0, "hydrogen_loss", 1)
    elif pattern == "desulfamoylation":
        change("O", "S", 1, 0, "broken")
        change("O", "H", 0, 1, "hydrogen_gain", 1)
    elif pattern == "aromatization":
        change("C", "C", 1, 2, "order_increased")
        change("C", "S", 1, 0, "broken")
    elif pattern == "no_reaction":
        r["classification"]["labels"] = F("no observed transformation", [], evidence)
        r["classification"]["coverage"] = F(None, "not_applicable", evidence, method="inference",
                                                derivation="No reaction was observed; expected-product bonds were not assigned.")
        r["classification"]["limitations"] = ["No observed transformation to classify."]
        return
    else:
        r["classification"]["limitations"] = ["Source structure or atom correspondence is insufficient for bond classification."]
        return
    r["classification"]["bond_changes"] = [
        F(None, c, evidence, method="structure_derivation",
          derivation="Net bond-change class inferred from the source reaction drawing; atom mapping is unavailable.")
        for c in changes
    ]
    labels = {}
    for c in changes:
        if c["change_type"] in ("formed", "order_increased", "hydrogen_gain"):
            lab = bond_label(c["atoms"][0]["element"], c["atoms"][1]["element"], c["after_order"])
            labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(c["change_id"])
    r["classification"]["labels"] = F(None, list(labels.values()), evidence, method="structure_derivation",
                                             derivation="Taxonomy labels are derived from the partial net bond-change list.")
    r["classification"]["coverage"] = F(None, "partial", evidence, method="inference",
                                            derivation="No reviewed atom-mapped reaction SMILES; only source-visible net bond classes are recorded.")
    r["classification"]["limitations"] = [
        "Atom identities are not mapped; multiplicity and hydrogen bookkeeping may be incomplete.",
        "Labels are inferred from source drawings and remain needs-review.",
    ]


def reaction(suffix, page, record_type="specific_experiment", categories=(),
             name=None, transformation=None, location=None, evidence=None):
    evidence = evidence or ev(page)
    r = new_reaction(f"{DOC}:{suffix}", DOC, evidence)
    r["record_type"] = record_type
    if categories:
        r["experiment_categories"] = F(list(categories), list(categories), evidence)
    r["location"] = F(location or f"PDF page {page}", location or f"PDF page {page}", evidence)
    if name:
        r["reaction_name"] = F(name, name, evidence)
    if transformation:
        r["transformation_type"] = F(transformation, transformation, evidence)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    bundle["reactions"].append(r)
    return r


def solvent(r, evidence, volume=2.0, method="original_extraction"):
    participant(r, "toluene", "solvent", evidence, volume=volume)
    r["solvent_conditions"]["composition"] = F("toluene", [chem["toluene"]], evidence, method=method)
    r["solvent_conditions"]["total_volume"] = F(f"{volume} mL", q(volume, "mL"), evidence, method=method)


def mark_yield(r, product, percent, evidence, raw=None, basis="isolated yield"):
    value = q(percent, "%") if not isinstance(percent, dict) else percent
    outcome(r, "yield", product, value, evidence, basis=basis,
            yield_type="isolated" if basis == "isolated yield" else None,
            raw=raw or (f"{percent}%" if not isinstance(percent, dict) else percent.get("raw_text")))


def add_common_scope_workup(r, evidence):
    # Explicitly inherited from the page-6 General Procedure for compounds 3aa-3ea.
    for label, role, volume in [("water", "workup_agent", 10), ("EtOAc", "workup_agent", 5),
                                ("brine", "workup_agent", None), ("Na2SO4", "workup_agent", None),
                                ("silica_gel", "purification_agent", None),
                                ("hexane", "purification_agent", None)]:
        p = participant(r, label, role, evidence, volume=volume,
                        stage="purification" if role == "purification_agent" else "workup")
        if label == "EtOAc":
            p["charges"]["volume"] = F("3 × 5 mL", {"number": 15, "unit": "mL", "raw_text": "3 × 5 mL"}, evidence,
                                               method="inherited_general_condition")
    participant(r, "EtOAc", "purification_agent", evidence, stage="purification")
    add_step(r, "workup_quench_extract", evidence,
             operation="quench with water (10 mL); extract aqueous layer with ethyl acetate (3 × 5 mL)",
             stage="workup", method="inherited_general_condition")
    add_step(r, "workup_wash_dry_concentrate", evidence,
             operation="wash combined organic layer with brine; dry over anhydrous Na2SO4; concentrate under reduced pressure",
             stage="workup", method="inherited_general_condition")
    add_step(r, "purification_column", evidence,
             operation="silica gel column chromatography; ethyl acetate/hexane 10:90",
             stage="purification", method="inherited_general_condition")


# Four current-work generic branches (graphical abstract + Scheme 2c).
generic_ids = {}
for suffix, temp, base, transformation in [
    ("generic_dbu_50", 50, "DBU/Na2CO3", "[3 + 3] annulation to hydroxy 2-pyridyl acetate"),
    ("generic_dbu_100", 100, "DBU/Na2CO3", "[3 + 3] annulation to 2-pyridyl acetate"),
    ("generic_pph3_rt", "rt", "PPh3", "[4 + 2] benzannulation to sulfamoyloxy/sulfonamide teraryl"),
    ("generic_pph3_80", 80, "PPh3", "[4 + 2] benzannulation to phenolic/sulfonamide teraryl"),
]:
    r = reaction(suffix, 1, "generic_scheme", name="Lewis base-switched annulation",
                 transformation=transformation, location="Graphical abstract (page 1) and Scheme 2c (page 2)",
                 evidence=ev(1)+ev(2))
    r["procedure"]["original_text"] = F(f"{base}, toluene, {temp} °C" if temp != "rt" else f"{base}, toluene, rt",
                                           {"base": base, "solvent": "toluene", "temperature": temp}, ev(1)+ev(2))
    add_step(r, "generic_reaction", ev(1)+ev(2), temperature=temp,
             operation="generic reaction arrow; substrate identities vary")
    classify(r, "annulation33" if "dbu" in suffix else "annulation42", ev(1)+ev(2))
    r["review"]["unresolved"] = ["Generic scheme only; it is not counted as a specific experiment or assigned a yield."]
    generic_ids[suffix] = r["reaction_id"]

# Prior-work arrows are retained as literature examples, not attributed to this paper.
for page, panel, count in [(1, "Scheme 1", 3), (2, "Scheme 2a", 3), (2, "Scheme 2b", 2)]:
    for i in range(1, count+1):
        r = reaction(f"literature_{panel.lower().replace(' ', '_')}_{i}", page,
                     "literature_example", name=f"Prior-work arrow {i} in {panel}",
                     location=f"{panel}, arrow {i}, PDF page {page}")
        r["review"]["unresolved"] = [
            "Prior literature example; detailed substrate identity and conditions are not imported as this paper's experiment."
        ]
        r["classification"]["limitations"] = ["Prior-work drawing retained without atom mapping or a new-paper experiment claim."]

# Table 1: 15 optimization rows.
table1 = [
    (1, "DABCO", None, 25, "<5", 0, 0), (2, "DABCO", None, 25, 31, 0, 0),
    (3, "DABCO", "Na2CO3", 25, 43, 0, 0), (4, "DABCO", "Na2CO3", 50, 56, 0, 0),
    (5, "DABCO", "Na2CO3", 50, 62, 0, 0), (6, "DABCO", "Na2CO3", 50, 62, 0, 0),
    (7, "DBU", "Na2CO3", 50, 79, 0, 0), (8, "DBU", "K2CO3", 50, 72, 0, 0),
    (9, "DBU", "t-BuOK", 50, 75, 0, 0), (10, "DBU", "Na2CO3", 50, 79, 0, 0),
    (11, "PPh3", "Na2CO3", 50, 0, 53, 0), (12, "PPh3", "Na2CO3", 80, 0, 0, 73),
    (13, "PPh3", None, 80, 0, 0, 81), (14, "PPh3", None, 25, 0, 73, 0),
    (15, "PPh3", None, 80, 0, 0, 81),
]
for entry, base, additive, temp, y3, y4, y5 in table1:
    r = reaction(f"table1_entry_{entry:02d}", 2, categories=("optimization",),
                 name="Optimization of reactions giving 3aa, 4aa, and 5aa",
                 transformation="Lewis base-switched annulation", location=f"Table 1, entry {entry}")
    p1 = participant(r, "1a", "reactant", ev(2), amount=0.20)
    participant(r, "2a", "reactant", ev(2), amount=0.24)
    for prod in ("3aa", "4aa", "5aa"):
        participant(r, prod, "product", ev(2))
    loading = 20 if entry in (1, 13, 14, 15) else 50
    participant(r, base, "catalyst", ev(2), mol_percent=loading, reference=mol("1a"))
    if additive:
        equiv = 1.5 if 5 <= entry <= 10 else None
        participant(r, additive, "additive", ev(2), equivalents=equiv, reference=mol("1a") if equiv else None)
    solvent(r, ev(2), 2.0)
    hours = 24 if entry in (10, 15) else 12
    add_step(r, "reaction", ev(2), operation="stir in oil bath", temperature=temp, time=hours)
    values = [("3aa", y3), ("4aa", y4), ("5aa", y5)]
    for prod, value in values:
        if value == "<5":
            mark_yield(r, prod, {"number": 5, "unit": "%", "qualifier": "lt", "raw_text": "<5%"}, ev(2), raw="<5%")
        else:
            mark_yield(r, prod, value, ev(2))
    pattern = "annulation33" if (isinstance(y3, int) and y3 > 0) or y3 == "<5" else "annulation42"
    classify(r, pattern, ev(2))
    r["general_scheme_ids"] = [generic_ids["generic_dbu_50"] if pattern == "annulation33" else
                               generic_ids["generic_pph3_rt" if y4 else "generic_pph3_80"]]
    if entry in (1, 2):
        note = ("Table cells show no additive, while footnote a states additive (1.0 equiv) for entries 1 and 2; "
                "the additive identity is not supplied and is not invented.")
        r["notes"].append(F(note, note, ev(2)))
        r["review"]["unresolved"].append("Table/footnote additive conflict for entry %d." % entry)
    if entry in (3, 4):
        note = ("Na2CO3 is listed in the additive column, but footnote a does not assign an amount to entries 3 and 4.")
        r["notes"].append(F(note, note, ev(2)))
        r["review"]["unresolved"].append("Na2CO3 amount is unreported for entry %d." % entry)

# A genuine negative control stated across the page-2/3 continuation.
r = reaction("control_no_tpp_no_na2co3", 3, categories=("control",),
             name="Control without TPP and Na2CO3", transformation="no reaction",
             location="Results text continued from page 2 to page 3", evidence=ev(2)+ev(3))
participant(r, "1a", "reactant", ev(2)+ev(3), amount=0.20)
participant(r, "2a", "reactant", ev(2)+ev(3), amount=0.24)
solvent(r, ev(2)+ev(3), 2.0)
add_step(r, "reaction", ev(2)+ev(3), operation="reaction attempted without TPP and Na2CO3",
         time=24)
outcome(r, "no_reaction", "1a", True, ev(2)+ev(3), basis=None, product_role="input material",
        raw="did not observe any reaction even after 24 h")
classify(r, "no_reaction", ev(2)+ev(3))
r["review"]["unresolved"] = ["Temperature and analytical detection basis for the no-reaction statement are not reported."]

# Table 2 substrate scope, with the page-6 General Procedure inherited explicitly.
product_masses = {"3aa": 53.0, "3ab": 53.0, "3ac": 59.0, "3ad": 54.0, "3ae": 63.2}
for prod, (imine, allene, yld) in scope3.items():
    evidence = ev(3) + ev(6)
    r = reaction(f"scope_table2_{prod}", 3, categories=("substrate_scope",),
                 name="[3 + 3] annulation of cyclic sulfamidate imine with δ-acetoxy allenoate",
                 transformation="[3 + 3] annulation to 2-pyridinyl acetate",
                 location=f"Table 2, product {prod}; General Procedure ii on page 6", evidence=evidence)
    participant(r, imine, "reactant", evidence, amount=0.20,
                mass={"1a":39.4}.get(imine) if prod in product_masses else None)
    participant(r, allene, "reactant", evidence, amount=0.24,
                mass={"2a":62.4,"2b":69.6,"2c":70.6,"2d":73.2,"2e":95.8}.get(allene)
                if prod in product_masses else None)
    participant(r, prod, "product", evidence)
    dbu_amount = 0.20 if imine == "1d" else 0.10
    participant(r, "DBU", "catalyst", evidence, amount=dbu_amount,
                mass=(30.4 if imine == "1d" else 15.2) if prod in product_masses else None)
    participant(r, "Na2CO3", "additive", evidence, amount=0.30,
                mass=31.8 if prod in product_masses else None)
    solvent(r, evidence, 2.0, method="inherited_general_condition")
    add_step(r, "reaction", evidence, operation="stir at 50 °C in an oil bath; monitor by TLC",
             temperature=50, time=12, method="inherited_general_condition")
    add_common_scope_workup(r, ev(6))
    r["procedure"]["general_procedure_reference"] = F(
        "General Procedure ii for compounds 3aa-3ea", "General Procedure ii", ev(6))
    r["procedure"]["addition_order"] = F(
        "charge imine, allenoate, Na2CO3 and 1.5 mL toluene; then add DBU in 0.5 mL toluene",
        ["imine", "allenoate", "Na2CO3", "1.5 mL toluene", "DBU in 0.5 mL toluene"], ev(6))
    mark_yield(r, prod, yld, evidence)
    if prod in product_masses:
        outcome(r, "product_amount", prod, q(product_masses[prod], "mg"), ev(6),
                basis="isolated product mass", raw=f"{product_masses[prod]} mg")
    classify(r, "annulation33", evidence)
    r["general_scheme_ids"] = [generic_ids["generic_dbu_50"]]
    r["review"]["unresolved"] = [
        "Primary reactant/product structures and all stereo elements remain unverified; source labels and drawings are retained."
    ]

# Scheme 3 scale-up experiments.
for prod, route, yld, mass in [("3ab", "dbu", 72, 261.0), ("5ab", "pph3", 73, 254.0)]:
    r = reaction(f"scaleup_{prod}", 3, categories=("scale_up",),
                 name=f"Scale-up synthesis of {prod}", transformation="scale-up annulation",
                 location=f"Scheme 3, scale-up to {prod}")
    participant(r, "1a", "reactant", ev(3), amount=1.0, mass=197.0)
    participant(r, "2b", "reactant", ev(3), amount=1.2, mass=348.0)
    participant(r, prod, "product", ev(3))
    if route == "dbu":
        participant(r, "DBU", "catalyst", ev(3), amount=0.5, mol_percent=50, reference=mol("1a"))
        participant(r, "Na2CO3", "additive", ev(3), amount=1.5, equivalents=1.5, reference=mol("1a"))
        temp = 50; pattern = "annulation33"; generic = "generic_dbu_50"
    else:
        participant(r, "PPh3", "catalyst", ev(3), amount=0.20, mol_percent=20, reference=mol("1a"))
        temp = 80; pattern = "annulation42"; generic = "generic_pph3_80"
    solvent(r, ev(3), 10.0)
    add_step(r, "reaction", ev(3), operation="scale-up reaction", temperature=temp, time=12)
    mark_yield(r, prod, yld, ev(3))
    outcome(r, "product_amount", prod, q(mass, "mg"), ev(3), basis="isolated product mass", raw=f"{mass} mg")
    classify(r, pattern, ev(3)); r["general_scheme_ids"] = [generic_ids[generic]]
    r["review"]["unresolved"] = ["Scale-up workup is discussed later in the paper, outside the page-1--6 visual coverage claim."]

# Narrative optimization and mechanistic/derivatization experiments.
for suffix, temp, yld, dbu_amount in [("sulfonyl_6a_2a_trace", 50, "<5", 0.10),
                                      ("sulfonyl_6a_2a_100c_25pct", 100, 25, 0.10)]:
    r = reaction(suffix, 3, categories=("optimization",), name="Preliminary [3 + 3] annulation of 6a with 2a",
                 transformation="[3 + 3] annulation to 7aa", location="Results text, page 3")
    participant(r, "6a", "reactant", ev(3), amount=0.20)
    participant(r, "2a", "reactant", ev(3), amount=0.24)
    participant(r, "7aa", "product", ev(3))
    participant(r, "DBU", "catalyst", ev(3), amount=dbu_amount)
    participant(r, "Na2CO3", "additive", ev(3), amount=0.30)
    solvent(r, ev(3), 2.0)
    add_step(r, "reaction", ev(3), operation="preliminary annulation", temperature=temp,
             time=12 if yld == "<5" else None)
    if yld == "<5":
        mark_yield(r, "7aa", {"number":5,"unit":"%","qualifier":"lt","raw_text":"<5%"}, ev(3), raw="trace amount (<5%)")
    else:
        mark_yield(r, "7aa", yld, ev(3), basis=None)
        r["review"]["unresolved"].append("Reaction time and yield measurement basis are not stated for the 25% observation.")
    classify(r, "annulation33", ev(3)); r["general_scheme_ids"] = [generic_ids["generic_dbu_100" if temp == 100 else "generic_dbu_50"]]

r = reaction("derivatization_4aa_to_5aa", 4, categories=("derivatization",),
             name="Thermal conversion of 4aa to 5aa", transformation="S-O cleavage/desulfamoylation",
             location="Results text page 3 and Table 3 arrow page 4", evidence=ev(3)+ev(4))
participant(r, "4aa", "reactant", ev(3)+ev(4), amount=0.20)
participant(r, "5aa", "product", ev(3)+ev(4))
solvent(r, ev(3)+ev(4), 2.0)
add_step(r, "reaction", ev(3)+ev(4), operation="heat 4aa in toluene", temperature=80, time=6)
outcome(r, "yield", "5aa", "quantitative", ev(3)+ev(4), basis=None, yield_type=None,
        raw="afforded 5aa quantitatively")
classify(r, "desulfamoylation", ev(3)+ev(4))
r["review"]["unresolved"] = ["A numeric yield and measurement basis are not reported for 'quantitatively'."]

# Table 3 scope.
for prod, (imine, allene, yld) in {**scope4, **scope5}.items():
    family4 = prod.startswith("4")
    temp = "rt" if family4 else 80
    r = reaction(f"scope_table3_{prod}", 4, categories=("substrate_scope",),
                 name="[4 + 2] annulation of cyclic sulfamidate imine with δ-acetoxy allenoate",
                 transformation="[4 + 2] benzannulation to teraryl",
                 location=f"Table 3, product {prod}")
    participant(r, imine, "reactant", ev(4), amount=0.20)
    participant(r, allene, "reactant", ev(4), amount=0.24)
    participant(r, prod, "product", ev(4))
    participant(r, "PPh3", "catalyst", ev(4), amount=0.04, mol_percent=20, reference=mol(imine))
    solvent(r, ev(4), 2.0)
    add_step(r, "reaction", ev(4), operation="[4 + 2] annulation", temperature=temp, time=12)
    mark_yield(r, prod, yld, ev(4)); classify(r, "annulation42", ev(4))
    r["general_scheme_ids"] = [generic_ids["generic_pph3_rt" if family4 else "generic_pph3_80"]]
    r["review"]["unresolved"] = ["Full structure and stereochemical representation remain source-drawing candidates, not verified SMILES."]

# Table 4 scope.
for prod, (imine, allene, yld) in scope7.items():
    r = reaction(f"scope_table4_{prod}", 4, categories=("substrate_scope",),
                 name="[3 + 3] annulation of cyclic sulfonyl imine with δ-acetoxy allenoate",
                 transformation="[3 + 3] annulation to 2-pyridinyl acetate", location=f"Table 4, product {prod}")
    participant(r, imine, "reactant", ev(4), amount=0.20)
    participant(r, allene, "reactant", ev(4), amount=0.24)
    participant(r, prod, "product", ev(4))
    participant(r, "DBU", "catalyst", ev(4), amount=0.20)
    participant(r, "Na2CO3", "additive", ev(4), amount=0.30)
    solvent(r, ev(4), 2.0)
    add_step(r, "reaction", ev(4), operation="one-pot [3 + 3] annulation in oil bath", temperature=100, time=12)
    mark_yield(r, prod, yld, ev(4)); classify(r, "annulation33", ev(4))
    r["general_scheme_ids"] = [generic_ids["generic_dbu_100"]]
    r["review"]["unresolved"] = ["Full structure and stereochemical representation remain source-drawing candidates, not verified SMILES."]

# Isolated intermediate 8bb and its conversion to 7bb.
r = reaction("mechanistic_8bb_isolation", 4, categories=("mechanistic_study",),
             name="Isolation of fused dihydropyridine 8bb", transformation="DBU-mediated annulation",
             location="Table 4 mechanistic arrow and results text")
participant(r, "6b", "reactant", ev(4), amount=0.20)
participant(r, "2b", "reactant", ev(4), amount=0.24)
participant(r, "8bb", "product", ev(4))
participant(r, "DBU", "reagent", ev(4), amount=0.30, equivalents=1.5, reference=mol("6b"))
solvent(r, ev(4), 2.0)
add_step(r, "reaction", ev(4), operation="stir in oil bath", temperature=60, time=12)
mark_yield(r, "8bb", 74, ev(4)); classify(r, "annulation33", ev(4))
r["notes"].append(F("structure confirmed by X-ray crystallography", "X-ray structure reported", ev(3)+ev(4)))
r["review"]["unresolved"] = ["The partial bond labels do not encode the complete fused-ring atom mapping."]

r = reaction("mechanistic_8bb_to_7bb", 4, categories=("mechanistic_study", "derivatization"),
             name="Conversion of 8bb to 7bb", transformation="SO2 elimination/aromatization",
             location="Table 4 mechanistic arrow and results text")
participant(r, "8bb", "reactant", ev(4), amount=0.10)
participant(r, "7bb", "product", ev(4))
participant(r, "DBU", "reagent", ev(4), amount=0.10, equivalents=1.0, reference=mol("8bb"))
participant(r, "Na2CO3", "additive", ev(4), amount=0.15, equivalents=1.5, reference=mol("8bb"))
solvent(r, ev(4), 1.0)
add_step(r, "reaction", ev(4), operation="heat under basic conditions", temperature=100, time=4)
outcome(r, "other", "7bb", "product identity supported by 1H NMR/HRMS; yield unreported", ev(3)+ev(4),
        basis="1H NMR and HRMS", raw="7bb (1H NMR/HRMS evidence)")
classify(r, "aromatization", ev(4))
r["review"]["unresolved"] = ["No isolated or spectroscopic yield is reported for conversion of 8bb to 7bb."]

# Table 5 scope.
for prod, (imine, allene, yld) in scope9.items():
    r = reaction(f"scope_table5_{prod}", 4, categories=("substrate_scope",),
                 name="[4 + 2] annulation of cyclic sulfonyl imine with δ-acetoxy allenoate",
                 transformation="[4 + 2] benzannulation to sulfonamide teraryl", location=f"Table 5, product {prod}")
    participant(r, imine, "reactant", ev(4), amount=0.20)
    participant(r, allene, "reactant", ev(4), amount=0.24)
    participant(r, prod, "product", ev(4))
    participant(r, "PPh3", "catalyst", ev(4), amount=0.04, mol_percent=20, reference=mol(imine))
    solvent(r, ev(4), 2.0)
    add_step(r, "reaction", ev(4), operation="[4 + 2] annulation", temperature=80, time=12)
    mark_yield(r, prod, yld, ev(4)); classify(r, "annulation42", ev(4))
    r["general_scheme_ids"] = [generic_ids["generic_pph3_80"]]
    r["review"]["unresolved"] = ["Full structure and stereochemical representation remain source-drawing candidates, not verified SMILES."]

# Scheme 4 is a proposed mechanism, not a synthetic experiment.
r = reaction("proposed_mechanism_scheme4", 5, "proposed_mechanism",
             name="Probable reaction pathways for the [3 + 3] and [4 + 2] annulations",
             transformation="author-proposed SN2'/annulation/aromatization pathways",
             location="Scheme 4, PDF page 5")
r["procedure"]["original_text"] = F(
    "Scheme 4: A1-A3 route to 3aa and B1-B7 routes to 4aa/5aa",
    {"[3+3]": ["A1", "A2", "A3", "3aa"],
     "[4+2]": ["B1", "B2", "B3", "B4", "B5", "4aa", "5aa", "B3'", "B4'", "B5'", "B6", "B7"]}, ev(5))
r["notes"].append(F("Probable Reaction Pathways", "Author labels the scheme as probable", ev(5)))
r["classification"]["limitations"] = ["Proposed elementary arrows are not classified as observed synthetic experiments."]
r["review"]["unresolved"] = ["Mechanistic intermediates are proposed and were not registered as isolated molecules on pages 1-6."]

# Page-6 preparation of allenoates 2d and 2e, with separately staged operations.
for prod, precursor, precursor_mass, yield_mass, yield_pct in [
    ("2d", "2d_alcohol", 2630.0, 1700.0, 55),
    ("2e", "2e_alcohol", 3570.0, 1910.0, 48),
]:
    r = reaction(f"preparation_{prod}", 6, categories=("other",),
                 name=f"Synthesis of δ-acetoxy allenoate {prod}",
                 transformation="acetylation of allenylic alcohol", location=f"Experimental section i, {prod}")
    participant(r, precursor, "reactant", ev(6), amount=10.0, mass=precursor_mass)
    participant(r, "Et3N", "reagent", ev(6), amount=20.0, volume=2.80, equivalents=2.0,
                reference=chem[precursor])
    participant(r, "acetyl_chloride", "reagent", ev(6), amount=12.0, volume=0.85,
                equivalents=1.2, reference=chem[precursor])
    participant(r, "DCM", "solvent", ev(6), volume=20.0)
    participant(r, prod, "product", ev(6))
    r["solvent_conditions"]["composition"] = F("dry dichloromethane", [chem["DCM"]], ev(6))
    r["solvent_conditions"]["total_volume"] = F("20 mL", q(20, "mL"), ev(6))
    add_step(r, "cool_and_pre_stir", ev(6),
             operation="add triethylamine to allenylic alcohol solution at 0 °C and stir",
             temperature=0, time={"number":30,"unit":"min"})
    # Override generic time formatting for minutes.
    r["steps"][-1]["time"] = F("30 min", q(30, "min"), ev(6))
    add_step(r, "acetyl_chloride_addition", ev(6),
             operation="slowly add acetyl chloride", temperature=0)
    r["steps"][-1]["time"] = F("over 5 min", q(5, "min"), ev(6))
    add_step(r, "hold_at_zero", ev(6), operation="stir after addition", temperature=0)
    r["steps"][-1]["time"] = F("40 min", q(40, "min"), ev(6))
    for label, role, volume in [("water", "workup_agent", None), ("DCM", "workup_agent", 25),
                                ("brine", "workup_agent", 20), ("Na2SO4", "workup_agent", None),
                                ("silica_gel", "purification_agent", None),
                                ("EtOAc", "purification_agent", None), ("hexane", "purification_agent", None)]:
        participant(r, label, role, ev(6), volume=volume,
                    stage="purification" if role == "purification_agent" else "workup")
    add_step(r, "workup_extract", ev(6), operation="extract aqueous layer with dichloromethane (3 × 25 mL)", stage="workup")
    add_step(r, "workup_wash", ev(6), operation="wash combined organic layer with brine (2 × 20 mL)", stage="workup")
    add_step(r, "workup_dry_concentrate", ev(6), operation="dry over anhydrous Na2SO4 and concentrate under reduced pressure", stage="workup")
    add_step(r, "purification", ev(6), operation="silica gel column chromatography with ethyl acetate/hexane 10:90", stage="purification")
    mark_yield(r, prod, yield_pct, ev(6)); outcome(r, "product_amount", prod, q(yield_mass, "mg"), ev(6),
                                                   basis="isolated product mass", raw=f"{yield_mass/1000:g} g")
    classify(r, "acetylation", ev(6))
    r["review"]["unresolved"] = [
        "Product is a diastereomeric mixture; individual configurations are not assigned and no single isomeric SMILES is asserted."
    ]

# Source context separates non-synthetic characterization/general-method content.
bundle["source_context"] = [
    {"context_id": f"{DOC}:context:title_abstract", "document_id": DOC, "pages": [1],
     "kind": "bibliographic_and_abstract", "title": TITLE,
     "evidence_ids": ev(1), "note": "Title and abstract; not counted as an experiment."},
    {"context_id": f"{DOC}:context:general_methods", "document_id": DOC, "pages": [5, 6],
     "kind": "general_methods", "evidence_ids": ev(5)+ev(6),
     "note": "Instrument, chromatography and crystallization statements; not standalone synthesis records."},
    {"context_id": f"{DOC}:context:characterization_page6", "document_id": DOC, "pages": [6],
     "kind": "characterization", "evidence_ids": ev(6),
     "note": "NMR/IR/HRMS/mp data are linked to compounds 2d, 2e and 3aa-3ae; no separate synthesis records were created for spectral peak lists."},
    {"context_id": f"{DOC}:context:general_procedure_ii", "document_id": DOC, "pages": [6],
     "kind": "general_procedure", "evidence_ids": ev(6),
     "note": "Inherited only by matching Table 2 scope experiments; not counted as an additional experiment."},
]

# Add an explicit note to every record about SI and machine review status.
for r in bundle["reactions"]:
    r["notes"].append(F(
        "Supporting information was not explicitly provided in this batch.",
        "No SI facts imported; source mapping status is unestablished.", r["evidence_ids"]))
    if not r["review"]["unresolved"]:
        r["review"]["unresolved"] = ["Machine-extracted record requires expert source/structure review."]

# Raw visual transcription retained separately from formal records.
raw = {
    "batch_id": "allene-stereo-v31-200-20261017",
    "document_id": DOC,
    "visual_pages": [1, 2, 3, 4, 5, 6],
    "table1": table1,
    "table2": scope3,
    "table3_products4": scope4,
    "table3_products5": scope5,
    "table4": scope7,
    "table5": scope9,
    "stereo_diagnostics": {
        "toluene": stereo_inspect("Cc1ccccc1"),
        "PPh3": stereo_inspect("P(c1ccccc1)(c1ccccc1)c1ccccc1"),
    },
    "limitations": [
        "Raw transcription is not a source-reviewed chemistry record.",
        "Structures not independently verified are intentionally left without SMILES.",
        "No supporting-information association was inferred.",
    ],
}
(OUT / "raw_model/visual_transcription_pages_0001_0006.json").write_text(
    json.dumps(raw, ensure_ascii=False, indent=2) + "\n")
(OUT / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")

coverage_regions = {
    1: [
        {"region": "title, abstract, and graphical abstract", "finding": "current-work Lewis-base-switched generic [3+3]/[4+2] branches; generic schemes only, not specific-yield experiments", "evidence_ids": ev(1)},
        {"region": "Scheme 1", "finding": "three prior-literature reaction arrows retained as literature_example records", "evidence_ids": ev(1)},
        {"region": "introduction text", "finding": "background and cited work; no additional present-paper experiment", "evidence_ids": ev(1)},
    ],
    2: [
        {"region": "Scheme 2a-b", "finding": "five prior-work arrows retained as literature examples, separate from the present work", "evidence_ids": ev(2)},
        {"region": "Scheme 2c", "finding": "present-work divergent generic DBU and PPh3 approaches; merged with the graphical-abstract generic records", "evidence_ids": ev(2)},
        {"region": "Table 1 rows 1-15", "finding": "15 independent optimization experiments transcribed row-by-row, with 3aa/4aa/5aa yields and isolated-yield basis", "evidence_ids": ev(2)},
        {"region": "Table 1 footnote", "finding": "base loadings, additive statements, 2.0 mL toluene, times, and isolated-yield footnote; additive conflicts preserved", "evidence_ids": ev(2)},
    ],
    3: [
        {"region": "Table 2", "finding": "17 independent [3+3] substrate-scope experiments mapped to products 3aa-3ea", "evidence_ids": ev(3)},
        {"region": "Scheme 3", "finding": "two independent 1.0 mmol scale-up experiments to 3ab and 5ab", "evidence_ids": ev(3)},
        {"region": "results text above and below Table 2", "finding": "no-catalyst control, 4aa-to-5aa conversion, and two preliminary 6a+2a optimization observations", "evidence_ids": ev(3)},
        {"region": "right-column mechanistic-study text", "finding": "8bb isolation and subsequent conversion to 7bb identified; quantitative conditions cross-linked to page 4", "evidence_ids": ev(3)},
    ],
    4: [
        {"region": "Table 3", "finding": "19 independent [4+2] substrate-scope rows (4 products of type 4 and 15 products of type 5), plus distinct 4aa-to-5aa arrow", "evidence_ids": ev(4)},
        {"region": "Table 4", "finding": "8 independent [3+3] scope rows and two mechanistic-study arrows involving 8bb", "evidence_ids": ev(4)},
        {"region": "Table 5", "finding": "6 independent [4+2] scope rows", "evidence_ids": ev(4)},
        {"region": "mechanistic discussion text", "finding": "author interpretation of C-S/S-O cleavage; not duplicated as an additional synthesis experiment", "evidence_ids": ev(4)},
    ],
    5: [
        {"region": "Scheme 4", "finding": "author-labelled probable A1-A3 and B1-B7 pathways retained as one proposed_mechanism record, not experimental reactions", "evidence_ids": ev(5)},
        {"region": "conclusions and General Methods", "finding": "summary and instrument/general handling context; no additional specific experiment", "evidence_ids": ev(5)},
    ],
    6: [
        {"region": "Experimental section i, preparation of 2d", "finding": "specific acetylation with separate 0 °C pre-stir, 5 min addition, and 40 min hold steps; 55% isolated yield", "evidence_ids": ev(6)},
        {"region": "Experimental section i, preparation of 2e", "finding": "specific acetylation following 2d procedure; 48% isolated yield and dr 1:0.9", "evidence_ids": ev(6)},
        {"region": "General Procedure ii", "finding": "inherited only into matching Table 2 experiments, including addition order, workup, and column eluent", "evidence_ids": ev(6)},
        {"region": "individual entries 3aa-3ae", "finding": "product-specific charges, isolated masses/yields, and characterization linked to existing Table 2 experiment records rather than duplicated", "evidence_ids": ev(6)},
        {"region": "NMR/IR/HRMS/mp blocks", "finding": "characterization context linked to molecules; no standalone synthesis records created from spectral peak lists", "evidence_ids": ev(6)},
    ],
}
page_unresolved = {
    1: ["Generic drawings do not provide complete instance structures or yields.", "Prior-work conditions are not imported beyond the visible arrow classification."],
    2: ["Table 1 footnote says additive 1.0 equiv for entries 1-2 although the additive cells are blank.", "Na2CO3 amount is not assigned by the footnote for entries 3-4."],
    3: ["The 25% preliminary 7aa observation lacks an explicit time and yield measurement basis.", "The word 'quantitatively' for 4aa to 5aa has no numeric value or measurement basis."],
    4: ["The 8bb-to-7bb conversion has NMR/HRMS identity evidence but no reported yield.", "Complete atom-mapped structures are unavailable, so bond classifications remain partial."],
    5: ["Scheme 4 is proposed; intermediates were not treated as isolated compounds without explicit evidence."],
    6: ["2d and 2e are unassigned diastereomeric mixtures; individual configurations are not reported.", "No SI was explicitly provided or associated."],
}
coverage_pages = []
for page in range(1, 7):
    eid = ev(page)[0]
    record_ids = [r["reaction_id"] for r in bundle["reactions"] if eid in r["evidence_ids"]]
    coverage_pages.append({"page": page, "state": "attempted", "regions": coverage_regions[page],
                           "record_ids": record_ids, "unresolved": page_unresolved[page]})
coverage = {"document_id": DOC, "semantic_accuracy_verified": False, "pages": coverage_pages}
(OUT / "coverage").mkdir(exist_ok=True)
(OUT / "coverage/pages_0001_0006.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
                  "source_context": len(bundle["source_context"])}, ensure_ascii=False))
