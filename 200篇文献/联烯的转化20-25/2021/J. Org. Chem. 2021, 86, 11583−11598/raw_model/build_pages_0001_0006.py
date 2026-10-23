#!/usr/bin/env python3
"""Build the v3.1 submission for source PDF pages 1--6.

This is a session-local extraction script.  It does not change the shared
schema, validators, or scientific rules.
"""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))

import standard_records as sr  # noqa: E402

DOC = "doc_ee191f89470cf89e"
RUN = ROOT / "runs/current"


def ev(page: int, kind: str = "image") -> str:
    return f"{DOC}:page{page}:{kind}"


def q(number, unit, qualifier=None, raw_text=None):
    value = {"number": number, "unit": unit}
    if qualifier:
        value["qualifier"] = qualifier
    if raw_text:
        value["raw_text"] = raw_text
    return value


def rng(low, high, unit, raw_text):
    return {"minimum": low, "maximum": high, "unit": unit, "raw_text": raw_text}


def sf(raw=None, normalized=None, evidence=(), method="original_extraction",
       status=None, derivation=None, alternatives=()):
    return sr.fact(raw=raw, normalized=normalized, evidence=evidence, method=method,
                   status=status, derivation=derivation, alternatives=alternatives)


template = json.loads((RUN / "raw_model/pages_0001_0006_init.json").read_text())
bundle = copy.deepcopy(template)
bundle["molecules"] = []
bundle["reactions"] = []
bundle["source_context"] = []
molecules = {}


def merge_evidence(fact, evidence):
    for eid in evidence:
        if eid not in fact["evidence_ids"]:
            fact["evidence_ids"].append(eid)


def add_molecule(label, evidence, profile="unresolved", smiles=None, name=None):
    mid = f"{DOC}:mol:{label}"
    evidence = list(dict.fromkeys(evidence))
    if mid in molecules:
        m = molecules[mid]
        merge_evidence(m["occurrences"][0]["compound_label"], evidence)
        for eid in evidence:
            if eid not in m["occurrences"][0]["evidence_ids"]:
                m["occurrences"][0]["evidence_ids"].append(eid)
        return mid
    m = sr.new_molecule(mid, DOC, evidence)
    m["occurrences"][0]["compound_label"] = sf(label, label, evidence)
    m["names"]["original"] = sf(name or label, name or label, evidence)
    m["identity_note"] = sf(
        "Identity is registered from the paper label/drawing; machine records require expert review.",
        "Identity is registered from the paper label/drawing; machine records require expert review.",
        evidence,
    )
    st = m["stereochemistry"]
    if smiles:
        diag = sr.stereo_inspect(smiles)
        m["identity_status"] = "source_identified"
        m["representation_type"] = "smiles"
        m["smiles"] = sf(smiles, diag["canonical_isomeric_smiles"], evidence,
                         method="structure_derivation",
                         derivation="Candidate SMILES derived from the unambiguous common reagent/solvent identity; RDKit syntax/stereo diagnostics only, not source structure verification.")
        st["atom_centers"] = sf(None, diag["atom_centers"], evidence,
                                method="structure_derivation",
                                derivation="RDKit standard stereo inspection of the candidate SMILES.")
        st["double_bonds"] = sf(None, diag["double_bonds"], evidence,
                                method="structure_derivation",
                                derivation="RDKit standard stereo inspection of the candidate SMILES.")
        st["chirality_status"] = sf(None, "achiral", evidence, method="structure_derivation",
                                    derivation="Candidate graph has no tetrahedral or axial stereogenic element.")
        st["geometry_status"] = sf(None, "not_applicable", evidence, method="structure_derivation",
                                   derivation="Candidate graph has no source-relevant configurable alkene geometry.")
        st["representation_status"] = sf(None, "not_applicable", evidence, method="structure_derivation",
                                         derivation="No stereochemical encoding is applicable to this candidate graph.")
        st["isomeric_smiles"] = sf(None, diag["canonical_isomeric_smiles"], evidence,
                                   method="structure_derivation",
                                   derivation="Canonical isomeric SMILES returned by standard stereo inspection.")
        st["identity_key"] = sf(None, diag["stereo_graph_key"], evidence,
                                method="structure_derivation",
                                derivation="Stereo graph key from standard stereo inspection; not a sample identity.")
    else:
        m["identity_status"] = "source_identified"
        m["representation_type"] = "unresolved"
        m["smiles"] = sf()
        st["representation_status"] = sf(None, "not_encoded", evidence, method="inference",
                                         derivation="No source-verified complete SMILES was assigned from the drawing.")
        if profile == "allene_unspecified":
            st["axial"] = sf(
                "Allene shown; no axial descriptor or enantiomeric composition is reported.",
                {"axis_type": "allene", "assignment": "unspecified"}, evidence)
            st["chirality_status"] = sf(None, "unspecified", evidence, method="inference",
                                        derivation="The allene drawing has no axial stereochemical assignment; absence of a label is not treated as achiral or racemic.")
            st["geometry_status"] = sf(None, "not_applicable", evidence, method="inference",
                                       derivation="Allene axial stereochemistry is tracked under axial, not E/Z.")
        elif profile == "center_unspecified":
            st["atom_centers"] = sf(
                "Drawing contains a potentially stereogenic tetrahedral center without R/S, wedge, racemate, or composition assignment.",
                [{"assignment": "unspecified", "location": "source drawing"}], evidence)
            st["chirality_status"] = sf(None, "unspecified", evidence, method="inference",
                                        derivation="Potential stereocenter is present but configuration is not assigned in the inspected region.")
            st["geometry_status"] = sf(None, "unresolved", evidence, method="inference",
                                       derivation="Complete geometry assessment requires a source-verified structure encoding.")
        elif profile == "relative_spiro":
            st["atom_centers"] = sf(
                "Scheme 5 uses wedge/dashed bonds for the relative arrangement at the spirocyclopentene product.",
                [{"assignment": "relative configuration shown", "absolute_configuration": "unassigned"}], evidence)
            st["assignment"] = sf(
                "Relative configuration shown by wedge/dashed bonds; absolute configuration is not assigned.",
                {"type": "relative", "absolute": "unassigned"}, evidence)
            st["chirality_status"] = sf(None, "partially_specified", evidence, method="inference",
                                        derivation="Relative stereochemistry is drawn, while absolute configuration remains unassigned.")
            st["geometry_status"] = sf(None, "unresolved", evidence, method="inference",
                                       derivation="No explicit E/Z label is visible in Scheme 5; the ring alkene is not converted to an E/Z claim.")
        elif profile == "dienolate_ze":
            st["double_bonds"] = sf(
                "Experimental heading identifies (2Z,4E)-7ha.",
                [{"position": 2, "assignment": "Z"}, {"position": 4, "assignment": "E"}], [ev(13, "text")])
            st["chirality_status"] = sf(None, "unresolved", evidence, method="inference",
                                        derivation="No absolute tetrahedral assignment is made in the inspected structure.")
            st["geometry_status"] = sf(None, "specified", [ev(13, "text")], method="inference",
                                       derivation="Normalized from the explicit (2Z,4E) experimental heading.")
        else:
            st["chirality_status"] = sf(None, "unresolved", evidence, method="inference",
                                        derivation="Source drawing has not been converted into a source-verified complete structure.")
            st["geometry_status"] = sf(None, "unresolved", evidence, method="inference",
                                       derivation="Source drawing has not been converted into a source-verified complete structure.")
    molecules[mid] = m
    bundle["molecules"].append(m)
    return mid


AGENTS = {
    "PPh3": "P(c1ccccc1)(c1ccccc1)c1ccccc1",
    "pClP": "P(c1ccc(Cl)cc1)(c1ccc(Cl)cc1)c1ccc(Cl)cc1",
    "Bu3P": "P(CCCC)(CCCC)CCCC",
    "pFP": "P(c1ccc(F)cc1)(c1ccc(F)cc1)c1ccc(F)cc1",
    "Ph2PMe": "CP(c1ccccc1)c1ccccc1",
    "pTol3P": "P(c1ccc(C)cc1)(c1ccc(C)cc1)c1ccc(C)cc1",
    "Cs2CO3": "[Cs+].[Cs+].[O-]C([O-])=O",
    "K2CO3": "[K+].[K+].[O-]C([O-])=O",
    "tBuOK": "[K+].CC(C)(C)[O-]",
    "PhMe": "Cc1ccccc1",
    "DCE": "ClCCCl",
    "DMF": "CN(C)C=O",
    "MeOH": "CO",
    "Dioxane": "O1CCOCC1",
    "H2O": "O",
    "MeCN": "CC#N",
    "TsCl": "Cc1ccc(S(=O)(=O)Cl)cc1",
    "DCM": "ClCCl",
    "pyridine": "n1ccccc1",
    "acetyl_chloride": "CC(=O)Cl",
    "ethanol": "CCO",
    "EtOAc": "CCOC(C)=O",
    "Na2SO4": "[Na+].[Na+].[O-]S(=O)(=O)[O-]",
    "AgOAc": "[Ag+].CC(=O)[O-]",
}


def agent(label, evidence):
    return add_molecule(label, evidence, smiles=AGENTS.get(label), name=label)


def compound(label, evidence, profile="unresolved"):
    return add_molecule(label, evidence, profile=profile, name=label)


def participant(pid, mid, role, evidence, amount=None, mass=None, volume=None,
                equiv=None, molpct=None, ref=None, stage="reaction", product_role=None):
    p = sr.new_participant(pid, mid, role, evidence, stage=stage)
    if amount is not None:
        p["charges"]["amount"] = sf(str(amount), amount, evidence)
    if mass is not None:
        p["charges"]["mass"] = sf(str(mass), mass, evidence)
    if volume is not None:
        p["charges"]["volume"] = sf(str(volume), volume, evidence)
    if equiv is not None:
        p["charges"]["equivalents"] = sf(str(equiv), equiv, evidence)
    if molpct is not None:
        p["charges"]["mol_percent"] = sf(str(molpct), molpct, evidence)
    if ref is not None:
        p["equivalence_reference"] = sf(ref, ref, evidence)
    if product_role:
        p["product_role"] = sf(product_role, product_role, evidence)
    return p


def classification(reaction, evidence, bonds=(), coverage="partial"):
    c = reaction["classification"]
    c["mapped_reaction_smiles"] = sf()
    if coverage == "not_applicable":
        c["bond_changes"] = []
        c["labels"] = sf([], [], evidence, method="inference",
                         derivation="No product-forming transformation was observed in this control.")
        c["coverage"] = sf(None, "not_applicable", evidence, method="inference",
                           derivation="No observed product-forming transformation to classify.")
        c["limitations"] = ["Negative control; no formed-bond taxonomy asserted."]
        c["model_votes"] = []
        return
    changes = []
    labels = {}
    for idx, (a, b, order) in enumerate(bonds, 1):
        change_id = f"b{idx}"
        change = {"change_id": change_id,
                  "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
                  "before_order": 0, "after_order": order, "change_type": "formed"}
        changes.append(sf(None, change, evidence, method="inference",
                          derivation="Bond type inferred from the inspected reaction drawing; atom mapping is not claimed."))
        label = sr.bond_label(a, b, order)
        labels.setdefault(label["code"], dict(label, supporting_change_ids=[]))["supporting_change_ids"].append(change_id)
    c["bond_changes"] = changes
    c["labels"] = sf(None, list(labels.values()), evidence, method="inference",
                     derivation="Multi-label taxonomy derived from the source reaction drawing without atom mapping.")
    c["coverage"] = sf(None, coverage, evidence, method="inference",
                       derivation="Only visually supportable newly formed heavy-atom bonds are listed; implicit-H and full atom mapping remain unresolved.")
    c["limitations"] = [
        "No source-verified atom mapping was available.",
        "Implicit-hydrogen changes and complete atom balance are unresolved; unchanged or order-decreased bonds are not labeled as newly formed.",
    ]
    c["model_votes"] = []


def base_reaction(rid, evidence, record_type, categories, location, name, transformation,
                  unresolved=()):
    r = sr.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = sf(categories, categories, evidence)
    r["location"] = sf(location, location, evidence)
    r["reaction_name"] = sf(name, name, evidence)
    r["transformation_type"] = sf(transformation, transformation, evidence)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": list(unresolved)}
    r["completeness"] = "partial"
    bundle["reactions"].append(r)
    return r


def set_solvent(r, mid, evidence, volume=None, concentration=None):
    r["solvent_conditions"]["composition"] = sf("PhMe" if mid.endswith(":PhMe") else mid.split(":")[-1], [mid], evidence)
    if volume is not None:
        r["solvent_conditions"]["total_volume"] = sf(str(volume), volume, evidence)
    if concentration is not None:
        r["solvent_conditions"]["concentration"] = sf(str(concentration), concentration, evidence)
        r["solvent_conditions"]["concentration_basis"] = sf(
            "reported reaction concentration; denominator not explicitly defined",
            "reported reaction concentration; denominator not explicitly defined", evidence)


def add_step(r, sid, evidence, pids, operation, temperature=None, time=None, stage="reaction", atmosphere=None):
    s = sr.new_step(sid, evidence, stage=stage)
    s["participant_ids"] = list(pids)
    s["operation"] = sf(operation, operation, evidence)
    if temperature is not None:
        s["temperature"] = sf(str(temperature), temperature, evidence)
    if time is not None:
        s["time"] = sf(str(time), time, evidence)
    if atmosphere is not None:
        s["atmosphere"] = sf(atmosphere, atmosphere, evidence)
    r["steps"].append(s)


def add_outcome(r, oid, kind, target, evidence, value=None, raw=None, basis=None,
                yield_type=None, scope=None, product_role="desired product"):
    o = sr.new_outcome(oid, kind, [target] if target else [], evidence)
    if value is not None or raw is not None:
        o["value"] = sf(raw if raw is not None else str(value), value, evidence)
    if basis is not None:
        o["measurement_basis"] = sf(basis, basis, evidence)
    if yield_type is not None:
        o["yield_type"] = sf(yield_type, yield_type, evidence)
    if scope is not None:
        o["scope"] = sf(scope, scope, evidence)
    if target:
        o["product_role"] = sf(product_role, product_role, evidence)
    r["outcomes"].append(o)


def add_note(r, text, evidence, method="original_extraction", derivation=None):
    r["notes"].append(sf(text if method == "original_extraction" else None, text, evidence,
                         method=method, derivation=derivation))


def conflict_fact(values):
    alts = [sf(raw, norm, evidence) for raw, norm, evidence in values]
    return sf(None, None, (), status="conflict", alternatives=alts)


# Common source compounds and agents.
E1, E2, E3, E4, E5, E6 = ([ev(i, "image"), ev(i, "text")] for i in range(1, 7))
for lbl, evidence in [("PPh3", E3), ("pClP", E3), ("Bu3P", E3), ("pFP", E3),
                      ("Ph2PMe", E3), ("pTol3P", E4 + E5), ("Cs2CO3", E3 + E4 + E6),
                      ("K2CO3", E3), ("tBuOK", E3), ("PhMe", E1 + E3 + E4 + E5 + E6),
                      ("DCE", E3 + E6), ("DMF", E3 + E4), ("MeOH", E3 + E4),
                      ("Dioxane", E3), ("H2O", E3 + E4), ("MeCN", E4), ("TsCl", E6)]:
    agent(lbl, evidence)


# Three current-work overview schemes (pages 1-2). They are not counted as specific experiments.
for suffix, title, bonds, page in [
    ("generic_3_series", "Phosphine-catalyzed [3+3] hetero-annulation to dihydro-alpha-carbolines", (("C", "C", 1), ("C", "N", 1)), 1),
    ("generic_4_series", "[3+3] annulation with tosyl migration and aromatization", (("C", "C", 1), ("C", "N", 1), ("C", "S", 1)), 1),
    ("generic_6_series", "[4+1] spiro-carbo-annulation", (("C", "C", 1), ("C", "C", 1)), 1),
]:
    evidence = [ev(page, "image"), ev(page, "text"), ev(2, "image")]
    r = base_reaction(suffix, evidence, "generic_scheme", [], "Graphical abstract and Scheme 2c",
                      title, title, ["Generic family drawing; individual substrate identities and full atom mappings are represented in specific records."])
    classification(r, evidence, bonds)

GEN3 = f"{DOC}:generic_3_series"
GEN4 = f"{DOC}:generic_4_series"
GEN6 = f"{DOC}:generic_6_series"


# Table 1: 23 independent optimization entries.
table1 = [
    (1, "PPh3", None, "PhMe", 25, 24, "45:00"),
    (2, "PPh3", None, "PhMe", 25, 24, "67:00"),
    (3, "PPh3", None, "PhMe", 25, 36, "83:00"),
    (4, "pClP", None, "PhMe", 25, 36, "61:00"),
    (5, "Bu3P", None, "PhMe", 25, 36, "57:00"),
    (6, "pFP", None, "PhMe", 25, 36, "58:00"),
    (7, "Ph2PMe", None, "PhMe", 25, 36, "65:00"),
    (8, "PPh3", None, "DCE", 25, 36, "26:00"),
    (9, "PPh3", None, "DMF", 25, 36, "29:00"),
    (10, "PPh3", None, "MeOH", 25, 36, None),
    (11, "PPh3", None, "Dioxane", 25, 36, "31:00"),
    (12, "PPh3", None, "H2O", 25, 36, None),
    (13, "PPh3", None, "PhMe", 80, 36, "61:<10"),
    (14, "PPh3", "Cs2CO3", "PhMe", 80, 36, "≤5:75"),
    (15, "PPh3", "K2CO3", "PhMe", 80, 36, "≤5:61"),
    (16, "PPh3", "tBuOK", "PhMe", 80, 36, "≤5:64"),
    (17, "PPh3", "Cs2CO3", "PhMe", 80, 36, "Complex"),
    (18, "PPh3", "Cs2CO3", "PhMe", 80, 36, "≤5:53"),
    (19, "PPh3", None, "PhMe", 25, 36, "61:00"),
    (20, "PPh3", None, "PhMe", 25, 36, "83:00"),
    (21, "PPh3", "Cs2CO3", "PhMe", 80, 36, "≤5:75"),
    (22, "PPh3", None, "PhMe", 25, 36, "83:00"),
    (23, "PPh3", "Cs2CO3", "PhMe", 80, 36, "≤5:75"),
]

for entry, cat, additive, solv, temp, hours, result in table1:
    evidence = E3 + [ev(2, "text")]
    r = base_reaction(f"table1_entry_{entry:02d}", evidence, "specific_experiment", ["optimization"],
                      f"Table 1, entry {entry}", "Optimization of 3aa/4aa formation",
                      "phosphine-catalyzed annulation",
                      ["Primary substrate/product structures remain drawing-based and need expert structure encoding."])
    r["general_scheme_ids"] = [GEN3, GEN4]
    m1 = compound("1a", evidence)
    m2 = compound("2a", evidence, "allene_unspecified")
    m3 = compound("3aa", evidence, "center_unspecified")
    m4 = compound("4aa", evidence)
    mc = agent(cat, evidence)
    ms = agent(solv, evidence)
    parts = [
        participant("p_1a", m1, "reactant", evidence, amount=q(0.20, "mmol")),
        participant("p_2a", m2, "reactant", evidence, amount=q(0.24, "mmol")),
        participant("p_cat", mc, "catalyst", evidence, ref=m1),
        participant("p_solvent", ms, "solvent", evidence, volume=q(2.0, "mL")),
    ]
    amount_map = {**{i: q(0.04, "mmol") for i in range(1, 18)}, 19: q(0.02, "mmol"),
                  20: q(0.10, "mmol"), 21: q(0.10, "mmol"), 22: q(0.20, "mmol"), 23: q(0.20, "mmol")}
    if entry == 18:
        parts[2]["charges"]["amount"] = conflict_fact([
            ("0.04 mmol for entries 1-18", q(0.04, "mmol"), E3),
            ("0.02 mmol for entries 18 and 19", q(0.02, "mmol"), E3),
        ])
        r["review"]["unresolved"].append("Table 1 footnote assigns entry 18 to both the 0.04 mmol and 0.02 mmol phosphine ranges.")
    else:
        parts[2]["charges"]["amount"] = sf(str(amount_map[entry]), amount_map[entry], evidence)
        molpct = amount_map[entry]["number"] / 0.20 * 100
        parts[2]["charges"]["mol_percent"] = sf(None, q(molpct, "mol%"), evidence,
                                                  method="inference",
                                                  derivation="Calculated from the stated phosphine amount and 0.20 mmol 1a denominator.")
    if additive:
        ma = agent(additive, evidence)
        parts.append(participant("p_additive", ma, "additive", evidence, amount=q(0.40, "mmol"),
                                 equiv=q(2.0, "equiv"), ref=m1))
    include3 = result not in (None, "Complex") and result.split(":")[0] not in ("00",)
    include4 = result not in (None, "Complex") and ":" in result and result.split(":")[1] not in ("00",)
    if include3:
        parts.append(participant("p_3aa", m3, "product", evidence, product_role="product 3aa"))
    if include4:
        parts.append(participant("p_4aa", m4, "product", evidence, product_role="product 4aa"))
    r["participants"] = parts
    r["reaction_scale"] = sf("1a (0.20 mmol)", {"reference_molecule_id": m1, "amount": q(0.20, "mmol")}, evidence)
    set_solvent(r, ms, evidence, volume=q(2.0, "mL"))
    add_step(r, "reaction", evidence, [p["participant_id"] for p in parts if p["stage"] == "reaction"],
             "combine and stir", q(temp, "°C"), q(hours, "h"))
    if entry == 2:
        r["procedure"]["pretreatment"] = sf("degassed toluene", "degassed toluene", [ev(2, "text")])
    if result and result not in ("Complex",):
        left, right = result.split(":")
        if include3:
            val = q(float(left[1:]), "%", "le", left + "%") if left.startswith("≤") else q(float(left), "%")
            add_outcome(r, "yield_3aa", "yield", m3, evidence, val, raw=left + "%", basis="isolated yield", yield_type="isolated")
        if include4:
            val = q(float(right[1:]), "%", "lt", right + "%") if right.startswith("<") else q(float(right), "%")
            add_outcome(r, "yield_4aa", "yield", m4, evidence, val, raw=right + "%", basis="isolated yield", yield_type="isolated")
    elif result == "Complex":
        add_outcome(r, "complex_mixture", "other", None, evidence, "Complex", raw="Complex", scope="unresolved product mixture")
    else:
        r["review"]["unresolved"].append("Table 1 yield cell is blank; no zero yield is inferred.")
    if result in (None, "Complex"):
        classification(r, evidence, (), coverage="partial")
        r["classification"]["labels"] = sf(None, None, evidence, method="inference",
                                             derivation="No identified product outcome is available for formed-bond classification.")
    else:
        bonds = [("C", "C", 1), ("C", "N", 1)] + ([("C", "S", 1)] if include4 else [])
        classification(r, evidence, bonds)


# Scheme 3 substrate scope.
scheme3 = [
    ("3aa", "1a", "2a", 83), ("3ab", "1a", "2b", 85), ("3ac", "1a", "2c", 78),
    ("3ad", "1a", "2d", 77), ("3ae", "1a", "2e", 76), ("3af", "1a", "2f", 75),
    ("3ag", "1a", "2g", 72), ("3ah", "1a", "2h", 73), ("3ai", "1a", "2i", 75),
    ("3aj", "1a", "2j", 71), ("3ba", "1b", "2a", 84), ("3ca", "1c", "2a", 78),
    ("3da", "1d", "2a", 73), ("3ea", "1e", "2a", 74), ("3fa", "1f", "2a", 71),
]
for prod, sub1, sub2, yld in scheme3:
    evidence = E4
    r = base_reaction(f"scheme3_{prod}", evidence, "specific_experiment", ["substrate_scope"],
                      f"Scheme 3, {prod}", "[3+3] hetero-annulation", "[3+3] hetero-annulation",
                      ["Complete source-verified structures and absolute stereochemistry are unresolved."])
    r["general_scheme_ids"] = [GEN3]
    m1, m2 = compound(sub1, evidence), compound(sub2, evidence, "allene_unspecified")
    mp, mc, ms = compound(prod, evidence, "center_unspecified"), agent("PPh3", evidence), agent("PhMe", evidence)
    r["participants"] = [participant("p_sub1", m1, "reactant", evidence, amount=q(0.20, "mmol")),
                         participant("p_sub2", m2, "reactant", evidence, amount=q(0.24, "mmol")),
                         participant("p_cat", mc, "catalyst", evidence, amount=q(0.04, "mmol"), molpct=q(20, "mol%"), ref=m1),
                         participant("p_solvent", ms, "solvent", evidence),
                         participant("p_product", mp, "product", evidence, product_role="desired product")]
    r["reaction_scale"] = sf("1a-f (0.20 mmol)", {"reference_molecule_id": m1, "amount": q(0.20, "mmol")}, evidence)
    set_solvent(r, ms, evidence, concentration=q(0.1, "M"))
    add_step(r, "reaction", evidence, [p["participant_id"] for p in r["participants"]], "stir", q(25, "°C"), q(36, "h"))
    add_outcome(r, "yield", "yield", mp, evidence, q(yld, "%"), raw=f"{yld}%", basis="isolated yield", yield_type="isolated")
    if prod == "3fa":
        add_note(r, "X-ray structure reported for 3fa.", evidence)
    classification(r, evidence, (("C", "C", 1), ("C", "N", 1)))


# Scheme 4 substrate scope.
scheme4 = [("4aa", "1a", "2a", 75), ("4ab", "1a", "2b", 78), ("4ac", "1a", "2c", 73),
           ("4ad", "1a", "2d", 70), ("4af", "1a", "2f", 68), ("4ah", "1a", "2h", 64),
           ("4ak", "1a", "2k", 72), ("4aj", "1a", "2j", 68), ("4ba", "1b", "2a", 76),
           ("4ea", "1e", "2a", 65)]
for prod, sub1, sub2, yld in scheme4:
    evidence = E4
    r = base_reaction(f"scheme4_{prod}", evidence, "specific_experiment", ["substrate_scope"],
                      f"Scheme 4, {prod}", "[3+3] annulation accompanied by tosyl migration",
                      "annulation/tosyl migration/aromatization",
                      ["Full atom mapping, implicit-H changes, and complete structure encoding remain unresolved."])
    r["general_scheme_ids"] = [GEN4]
    m1, m2 = compound(sub1, evidence), compound(sub2, evidence, "allene_unspecified")
    mp, mc, ma, ms = compound(prod, evidence), agent("PPh3", evidence), agent("Cs2CO3", evidence), agent("PhMe", evidence)
    r["participants"] = [participant("p_sub1", m1, "reactant", evidence, amount=q(0.20, "mmol")),
                         participant("p_sub2", m2, "reactant", evidence, amount=q(0.24, "mmol")),
                         participant("p_cat", mc, "catalyst", evidence, amount=q(0.04, "mmol"), molpct=q(20, "mol%"), ref=m1),
                         participant("p_base", ma, "additive", evidence, amount=q(0.40, "mmol"), equiv=q(2.0, "equiv"), ref=m1),
                         participant("p_solvent", ms, "solvent", evidence),
                         participant("p_product", mp, "product", evidence, product_role="desired product")]
    r["reaction_scale"] = sf("1a,b or 1e (0.20 mmol)", {"reference_molecule_id": m1, "amount": q(0.20, "mmol")}, evidence)
    set_solvent(r, ms, evidence, concentration=q(0.1, "M"))
    add_step(r, "reaction", evidence, [p["participant_id"] for p in r["participants"]], "stir", q(80, "°C"), q(36, "h"))
    add_outcome(r, "yield", "yield", mp, evidence, q(yld, "%"), raw=f"{yld}%", basis="isolated yield", yield_type="isolated")
    if prod in ("4aa", "4ad"):
        add_note(r, f"X-ray structure reported for {prod}.", evidence)
    classification(r, evidence, (("C", "C", 1), ("C", "N", 1), ("C", "S", 1)))


# Table 2 optimization. The source footnote says 2a although the drawing and prose identify 5a.
table2 = [
    (1, "PPh3", "PhMe", 25, 47, "20:1"), (2, "PPh3", "PhMe", 80, 69, "20:1"),
    (3, "Ph2PMe", "PhMe", 80, 54, "20:1"), (4, "pClP", "PhMe", 80, 49, "20:1"),
    (5, "Bu3P", "PhMe", 80, 52, "20:1"), (6, "pFP", "PhMe", 80, 45, "20:1"),
    (7, "pTol3P", "PhMe", 80, 79, "20:1"), (8, "pTol3P", "MeCN", 80, 26, "20:1"),
    (9, "pTol3P", "DMF", 80, 29, "20:1"), (10, "pTol3P", "MeOH", 80, "trace", None),
    (11, "pTol3P", "H2O", 80, "trace", None), (12, "pTol3P", "PhMe", 80, 79, "20:1"),
]
for entry, cat, solv, temp, yld, dr in table2:
    evidence = E4
    r = base_reaction(f"table2_entry_{entry:02d}", evidence, "specific_experiment", ["optimization"],
                      f"Table 2, entry {entry}", "Optimization of [4+1] spiro-annulation to 6aa",
                      "[4+1] spiro-carbo-annulation",
                      ["Table 2 footnote names 2a, while the table drawing and narrative identify beta-prime-acetoxyallenoate 5a; 5a is used as the participant but the source conflict is retained.",
                       "Absolute configuration is unassigned; dr order follows the drawn major product versus other diastereomer(s)."])
    r["general_scheme_ids"] = [GEN6]
    m1, m5, mp = compound("1a", evidence), compound("5a", evidence, "center_unspecified"), compound("6aa", evidence, "relative_spiro")
    mc, ms = agent(cat, evidence), agent(solv, evidence)
    r["participants"] = [participant("p_1a", m1, "reactant", evidence, amount=q(0.20, "mmol")),
                         participant("p_5a", m5, "reactant", evidence, amount=q(0.24, "mmol")),
                         participant("p_cat", mc, "catalyst", evidence, amount=q(0.04, "mmol"), molpct=q(20, "mol%"), ref=m1),
                         participant("p_solvent", ms, "solvent", evidence, volume=q(2.0, "mL")),
                         participant("p_product", mp, "product", evidence, product_role="6aa")]
    if entry == 12:
        ma = agent("Cs2CO3", evidence)
        r["participants"].append(participant("p_additive", ma, "additive", evidence, amount=q(0.40, "mmol"), equiv=q(2.0, "equiv"), ref=m1))
    r["reaction_scale"] = sf("1a (0.20 mmol)", {"reference_molecule_id": m1, "amount": q(0.20, "mmol")}, evidence)
    set_solvent(r, ms, evidence, volume=q(2.0, "mL"))
    r["procedure"]["modifications"] = conflict_fact([
        ("table drawing and narrative: 5a (0.24 mmol)", {"substrate": m5, "amount": q(0.24, "mmol")}, evidence),
        ("Table 2 footnote: 2a (0.24 mmol)", {"substrate": compound("2a", evidence, "allene_unspecified"), "amount": q(0.24, "mmol")}, evidence),
    ])
    add_step(r, "reaction", evidence, [p["participant_id"] for p in r["participants"]], "stir", q(temp, "°C"), q(36, "h"))
    if yld == "trace":
        add_outcome(r, "trace_6aa", "trace", mp, evidence, "trace", raw="trace", basis="isolated yield footnote; numeric amount unreported")
    else:
        add_outcome(r, "yield_6aa", "yield", mp, evidence, q(yld, "%"), raw=f"{yld}%", basis="isolated yield", yield_type="isolated")
    if dr:
        add_outcome(r, "dr_6aa", "selectivity", mp, evidence, dr, raw=dr,
                    basis=None, scope={"metric": "dr", "isomer_order": "drawn major 6aa : other diastereomer(s)", "configuration": "absolute unassigned"})
        r["review"]["unresolved"].append("Diastereomeric-ratio measurement method is not reported in Table 2.")
    classification(r, evidence, (("C", "C", 1), ("C", "C", 1)))


# Scheme 5 scope with product-specific dr outcomes.
scheme5 = [
    ("6aa", "1a", "5a", 79, "20:1", True), ("6ab", "1a", "5b", 83, ">20:1", False),
    ("6ac", "1a", "5c", 76, "20:1", False), ("6ad", "1a", "5d", 71, "20:1", False),
    ("6ae", "1a", "5e", 65, ">20:1", True), ("6af", "1a", "5f", 70, ">20:1", False),
    ("6ag", "1a", "5g", 69, "20:1", False), ("6ba", "1b", "5a", 82, ">20:1", False),
    ("6ca", "1c", "5a", 79, ">20:1", False), ("6da", "1d", "5a", 75, "20:1", False),
    ("6bf", "1b", "5f", 73, "20:1", False), ("6ga", "1g", "5a", 71, "20:1", False),
    ("6hb", "1h", "5b", 84, ">20:1", False),
]
for prod, sub1, sub2, yld, dr, xray in scheme5:
    evidence = E5 + [ev(1, "image")]
    r = base_reaction(f"scheme5_{prod}", evidence, "specific_experiment", ["substrate_scope"],
                      f"Scheme 5, {prod}", "[4+1] spiro-carbo-annulation", "[4+1] spiro-carbo-annulation",
                      ["Absolute configuration is not assigned; wedge/dashed bonds and dr are retained as relative stereochemical evidence.",
                       "Diastereomeric-ratio measurement method is not reported in Scheme 5."])
    r["general_scheme_ids"] = [GEN6]
    m1, m5, mp = compound(sub1, evidence), compound(sub2, evidence, "center_unspecified"), compound(prod, evidence, "relative_spiro")
    mc, ms = agent("pTol3P", evidence), agent("PhMe", evidence)
    r["participants"] = [participant("p_sub1", m1, "reactant", evidence, amount=q(0.20, "mmol")),
                         participant("p_sub2", m5, "reactant", evidence, amount=q(0.24, "mmol")),
                         participant("p_cat", mc, "catalyst", evidence, amount=q(0.04, "mmol"), molpct=q(20, "mol%"), ref=m1),
                         participant("p_solvent", ms, "solvent", evidence),
                         participant("p_product", mp, "product", evidence, product_role="desired spiro product")]
    r["reaction_scale"] = sf("1a-d, 1g-h (0.20 mmol)", {"reference_molecule_id": m1, "amount": q(0.20, "mmol")}, evidence)
    set_solvent(r, ms, evidence, concentration=q(0.1, "M"))
    add_step(r, "reaction", evidence, [p["participant_id"] for p in r["participants"]], "stir", q(80, "°C"), q(36, "h"))
    add_outcome(r, "yield", "yield", mp, evidence, q(yld, "%"), raw=f"{yld}%", basis="isolated yield", yield_type="isolated")
    add_outcome(r, "dr", "selectivity", mp, evidence, dr, raw=dr, scope={"metric": "dr", "isomer_order": f"drawn major {prod} : other diastereomer(s)", "configuration": "absolute unassigned"})
    if xray:
        add_note(r, f"Relative configuration was ascertained by X-ray crystal structure of {prod}.", evidence)
    classification(r, evidence, (("C", "C", 1), ("C", "C", 1)))


# Scheme 6: noncyclized alpha-umpolung addition product 7ha.
evidence = E5
r = base_reaction("scheme6_7ha", evidence, "specific_experiment", ["substrate_scope"],
                  "Scheme 6", "Isolation of dienolate 7ha", "alpha-umpolung addition",
                  ["Complete source-verified structures remain unresolved; 2Z,4E assignment is retained from the experimental heading as text evidence."])
m1h, m2a, mp = compound("1h", evidence), compound("2a", evidence, "allene_unspecified"), compound("7ha", evidence, "dienolate_ze")
mc, ms = agent("PPh3", evidence), agent("PhMe", evidence)
r["participants"] = [participant("p_1h", m1h, "reactant", evidence, amount=q(0.20, "mmol")),
                     participant("p_2a", m2a, "reactant", evidence, amount=q(0.24, "mmol")),
                     participant("p_cat", mc, "catalyst", evidence, molpct=q(20, "mol%"), ref=m1h),
                     participant("p_solvent", ms, "solvent", evidence),
                     participant("p_product", mp, "product", evidence, product_role="sole noncyclized product")]
r["reaction_scale"] = sf("1h (0.20 mmol)", {"reference_molecule_id": m1h, "amount": q(0.20, "mmol")}, evidence)
set_solvent(r, ms, evidence, concentration=q(0.1, "M"))
add_step(r, "reaction", evidence, [p["participant_id"] for p in r["participants"]], "stir", q(25, "°C"), q(36, "h"))
add_outcome(r, "yield", "yield", mp, evidence, q(78, "%"), raw="78%", basis="isolated product", yield_type="isolated")
classification(r, evidence, (("C", "C", 1),))


# Scheme 7 control experiments: six combinations in (a), two temperatures in (b), one TsCl experiment in (c).
for agent_name in ("Cs2CO3", "PPh3", "PPh3+Cs2CO3"):
    for temp in (80, 110):
        evidence = E6
        rid = f"scheme7a_{agent_name.replace('+', '_')}_{temp}C"
        r = base_reaction(rid, evidence, "specific_experiment", ["control", "mechanistic_study"],
                          "Scheme 7a", "Attempted conversion of 3aa to 4aa", "control experiment",
                          ["Starting amount of 3aa is not reported in Scheme 7a."])
        m3, m4, ms = compound("3aa", evidence, "center_unspecified"), compound("4aa", evidence), agent("PhMe", evidence)
        parts = [participant("p_3aa", m3, "reactant", evidence), participant("p_target", m4, "product", evidence, product_role="not observed target"), participant("p_solvent", ms, "solvent", evidence)]
        if "PPh3" in agent_name:
            parts.append(participant("p_cat", agent("PPh3", evidence), "catalyst", evidence, molpct=q(20, "mol%"), ref=m3))
        if "Cs2CO3" in agent_name:
            parts.append(participant("p_base", agent("Cs2CO3", evidence), "additive", evidence, equiv=q(2.0, "equiv"), ref=m3))
        r["participants"] = parts
        set_solvent(r, ms, evidence, concentration=q(0.1, "M"))
        add_step(r, "reaction", evidence, [p["participant_id"] for p in parts], "heat", q(temp, "°C"), q(36, "h"))
        add_outcome(r, "no_4aa", "no_reaction", m4, evidence, "4aa not obtained", raw="did not yield compound 4aa", scope="attempted 3aa to 4aa conversion")
        classification(r, evidence, coverage="not_applicable")

for temp in (25, 80):
    evidence = E6
    r = base_reaction(f"scheme7b_no_PPh3_{temp}C", evidence, "specific_experiment", ["control", "mechanistic_study"],
                      "Scheme 7b", "Annulation without PPh3", "negative control",
                      ["No product was observed; no product identity is invented."])
    m1, m2, ma, ms = compound("1a", evidence), compound("2a", evidence, "allene_unspecified"), agent("Cs2CO3", evidence), agent("PhMe", evidence)
    r["participants"] = [participant("p_1a", m1, "reactant", evidence, amount=q(0.20, "mmol")),
                         participant("p_2a", m2, "reactant", evidence, amount=q(0.24, "mmol")),
                         participant("p_base", ma, "additive", evidence, equiv=q(2.0, "equiv"), ref=m1),
                         participant("p_solvent", ms, "solvent", evidence)]
    r["reaction_scale"] = sf("1a (0.20 mmol)", {"reference_molecule_id": m1, "amount": q(0.20, "mmol")}, evidence)
    set_solvent(r, ms, evidence, concentration=q(0.1, "M"))
    add_step(r, "reaction", evidence, [p["participant_id"] for p in r["participants"]], "stir without PPh3", q(temp, "°C"), q(36, "h"))
    add_outcome(r, "no_product", "no_reaction", None, evidence, "no product observed", raw="did not observe any product", scope="absence of PPh3 control")
    classification(r, evidence, coverage="not_applicable")

evidence = E6
r = base_reaction("scheme7c_TsCl", evidence, "specific_experiment", ["control", "mechanistic_study"],
                  "Scheme 7c", "TsCl control for tosyl migration", "control experiment",
                  ["The <10% yield measurement basis is not reported in Scheme 7c."])
m1, m2, mp = compound("1a", evidence), compound("2a", evidence, "allene_unspecified"), compound("4aa", evidence)
parts = [participant("p_1a", m1, "reactant", evidence, amount=q(0.20, "mmol")),
         participant("p_2a", m2, "reactant", evidence, amount=q(0.24, "mmol")),
         participant("p_TsCl", agent("TsCl", evidence), "additive", evidence, amount=q(0.24, "mmol"), ref=m1),
         participant("p_cat", agent("PPh3", evidence), "catalyst", evidence, molpct=q(20, "mol%"), ref=m1),
         participant("p_base", agent("Cs2CO3", evidence), "additive", evidence, equiv=q(2.0, "equiv"), ref=m1),
         participant("p_solvent", agent("PhMe", evidence), "solvent", evidence),
         participant("p_product", mp, "product", evidence, product_role="4aa")]
r["participants"] = parts
r["reaction_scale"] = sf("1a (0.20 mmol)", {"reference_molecule_id": m1, "amount": q(0.20, "mmol")}, evidence)
set_solvent(r, agent("PhMe", evidence), evidence, concentration=q(0.1, "M"))
add_step(r, "reaction", evidence, [p["participant_id"] for p in parts], "stir", q(80, "°C"), q(36, "h"))
add_outcome(r, "yield_4aa", "yield", mp, evidence, q(10, "%", "lt", "<10%"), raw="<10%", basis=None, yield_type=None)
classification(r, evidence, (("C", "C", 1), ("C", "N", 1), ("C", "S", 1)))


# Scale-up experiments: page 6 prose anchors them; page 13 extracted text supplies the explicit charges and yields.
scaleups = [
    ("3aa", "2a", "PPh3", None, 25, 79, q(406.2, "mg")),
    ("4aa", "2a", "PPh3", "Cs2CO3", 80, 72, q(369.0, "mg")),
    ("6aa", "5a", "pTol3P", None, 80, 75, q(432.0, "mg")),
]
for prod, sub2, cat, additive, temp, yld, prod_mass in scaleups:
    evidence = E6 + [ev(13, "text")]
    r = base_reaction(f"scaleup_{prod}", evidence, "specific_experiment", ["scale_up"],
                      f"Scheme 9 referenced on page 6; scale-up procedure on PDF page 13 text, {prod}",
                      f"1.0 mmol scale synthesis of {prod}", "scale-up annulation",
                      ["PDF page 13 was used as extracted-text context only and is not included in the visually reviewed page-1-to-6 coverage claim.",
                       "Primary compound structures remain source-drawing identities without source-verified SMILES."])
    m1, m2 = compound("1a", evidence), compound(sub2, evidence, "allene_unspecified" if sub2 == "2a" else "center_unspecified")
    profile = "relative_spiro" if prod == "6aa" else "center_unspecified" if prod == "3aa" else "unresolved"
    mp, mc, ms = compound(prod, evidence, profile), agent(cat, evidence), agent("PhMe", evidence)
    m2_mass = q(312.0, "mg") if sub2 == "2a" else q(386.5, "mg")
    cat_mass = q(52.5, "mg") if cat == "PPh3" else q(61.0, "mg")
    parts = [participant("p_1a", m1, "reactant", evidence, amount=q(1.0, "mmol"), mass=q(314.0, "mg")),
             participant("p_sub2", m2, "reactant", evidence, amount=q(1.20, "mmol"), mass=m2_mass),
             participant("p_cat", mc, "catalyst", evidence, amount=q(0.20, "mmol"), mass=cat_mass, molpct=q(20, "mol%"), ref=m1),
             participant("p_solvent", ms, "solvent", evidence, volume=q(10.0, "mL")),
             participant("p_product", mp, "product", evidence, mass=prod_mass, product_role="desired product")]
    if additive:
        parts.append(participant("p_base", agent(additive, evidence), "additive", evidence,
                                 amount=q(2.0, "mmol"), mass=q(651.6, "mg"), equiv=q(2.0, "equiv"), ref=m1))
    r["participants"] = parts
    r["reaction_scale"] = sf("1a (1.0 mmol)", {"reference_molecule_id": m1, "amount": q(1.0, "mmol")}, evidence)
    set_solvent(r, ms, evidence, volume=q(10.0, "mL"))
    if prod == "3aa":
        add_step(r, "reaction", evidence, ["p_1a", "p_sub2", "p_cat", "p_solvent"], "combine and stir", "room temperature", q(36, "h"))
    else:
        add_step(r, "addition", evidence, ["p_sub2"], "add substrate solution gradually", q(80, "°C"), q(45, "min"))
        add_step(r, "reaction", evidence, [p["participant_id"] for p in parts if p["stage"] == "reaction"], "hold after addition", q(80, "°C"), q(36, "h"))
    add_step(r, "workup", evidence, [], "quench with water; extract with ethyl acetate; wash with brine; dry over Na2SO4; concentrate", stage="workup")
    add_step(r, "purification", evidence, ["p_product"], "silica gel chromatography", stage="purification")
    add_outcome(r, "yield", "yield", mp, evidence, q(yld, "%"), raw=f"{yld}%", basis="isolated after chromatography", yield_type="isolated")
    bonds = (("C", "C", 1), ("C", "C", 1)) if prod == "6aa" else (("C", "C", 1), ("C", "N", 1), *(((("C", "S", 1),) if prod == "4aa" else ())))
    classification(r, evidence, bonds)


# Page 6 precursor syntheses 1f and 1h.
for key in ("ethanol", "EtOAc", "Na2SO4", "Dioxane", "DCE", "AgOAc"):
    agent(key, E6)

evidence = E6
r = base_reaction("page6_synthesis_1f", evidence, "specific_experiment", ["other"],
                  "Experimental section i, page 6", "Synthesis of 1f", "sulfonamidoindole precursor synthesis",
                  ["Ag-free azide transfer mechanism and full atom mapping are not inferred from the procedure."])
start = compound("6-bromo-1-ethyl-1H-indole", evidence)
azide = compound("p-toluenesulfonyl azide", evidence)
prod = compound("1f", evidence)
solv = agent("Dioxane", evidence)
r["participants"] = [participant("p_indole", start, "reactant", evidence, amount=q(5.0, "mmol"), mass=q(1.11, "g")),
                     participant("p_azide", azide, "reactant", evidence, amount=q(7.0, "mmol"), mass=q(1.38, "g")),
                     participant("p_solvent", solv, "solvent", evidence, volume=q(10.0, "mL")),
                     participant("p_ethanol", agent("ethanol", evidence), "workup_agent", evidence, volume=q(40, "mL"), stage="workup"),
                     participant("p_product", prod, "product", evidence, mass=q(1.17, "g"), product_role="1f")]
r["reaction_scale"] = sf("6-bromo-1-ethyl-1H-indole (5.0 mmol)", {"reference_molecule_id": start, "amount": q(5.0, "mmol")}, evidence)
set_solvent(r, solv, evidence, volume=q(10.0, "mL"))
add_step(r, "reaction", evidence, ["p_indole", "p_azide", "p_solvent"], "stir with heating", rng(75, 80, "°C", "75-80 °C"), rng(18, 24, "h", "18-24 h"), atmosphere="nitrogen")
add_step(r, "cooling", evidence, ["p_product"], "cool", "room temperature", None)
add_step(r, "precipitation", evidence, ["p_ethanol", "p_product"], "add ethanol to crystallize/precipitate product", stage="workup")
add_step(r, "purification", evidence, ["p_product"], "silica gel column chromatography, ethyl acetate/hexane 1:4", stage="purification")
add_outcome(r, "yield", "yield", prod, evidence, q(60, "%"), raw="60%", basis="isolated product", yield_type="isolated")
add_outcome(r, "product_amount", "product_amount", prod, evidence, q(1.17, "g"), raw="1.17 g", basis="reported isolated mass")
classification(r, evidence, (("C", "N", 1),))

evidence = E6
r = base_reaction("page6_synthesis_1h", evidence, "specific_experiment", ["other"],
                  "Experimental section ii, page 6", "Synthesis of 1h", "sulfonamidoindole precursor synthesis",
                  ["Catalyst identities are source labels; AgNTf2 structure remains unresolved."])
start = compound("tert-butyl 1H-indole-1-carboxylate", evidence)
azide = compound("p-toluenesulfonyl azide", evidence)
agntf2 = compound("AgNTf2", evidence)
prod = compound("1h", evidence)
solv = agent("DCE", evidence)
r["participants"] = [participant("p_indole", start, "reactant", evidence, amount=q(2.30, "mmol"), mass=q(500, "mg")),
                     participant("p_azide", azide, "reactant", evidence, amount=q(2.54, "mmol"), mass=q(500, "mg")),
                     participant("p_AgNTf2", agntf2, "catalyst", evidence, amount=q(0.23, "mmol"), mass=q(89.2, "mg"), molpct=q(10, "mol%"), ref=start),
                     participant("p_AgOAc", agent("AgOAc", evidence), "catalyst", evidence, amount=q(0.26, "mmol"), mass=q(76.7, "mg"), ref=start),
                     participant("p_solvent", solv, "solvent", evidence, volume=q(5.0, "mL")),
                     participant("p_product", prod, "product", evidence, mass=q(100.0, "mg"), product_role="1h")]
r["reaction_scale"] = sf("tert-butyl 1H-indole-1-carboxylate (2.30 mmol)", {"reference_molecule_id": start, "amount": q(2.30, "mmol")}, evidence)
set_solvent(r, solv, evidence, volume=q(5.0, "mL"))
add_step(r, "reaction", evidence, [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"], "stir with heating", q(40, "°C"), q(12, "h"))
add_step(r, "workup", evidence, [], "cool to room temperature; filter through Celite; wash with EtOAc; remove solvent", stage="workup")
add_step(r, "purification", evidence, ["p_product"], "silica gel column chromatography, ethyl acetate/hexane 1:4", stage="purification")
add_outcome(r, "yield", "yield", prod, evidence, q(13, "%"), raw="13%", basis="isolated product", yield_type="isolated")
add_outcome(r, "product_amount", "product_amount", prod, evidence, q(100.0, "mg"), raw="100.0 mg", basis="reported isolated mass")
classification(r, evidence, (("C", "N", 1),))


# Cross-page general acetylation procedure begins on page 6; retain it as a partial generic scheme.
evidence = E6
r = base_reaction("page6_general_synthesis_5a_5g", evidence, "generic_scheme", [],
                  "Experimental section iii, page 6 (procedure continues on page 7)",
                  "General procedure for the synthesis of compounds 5a-g", "acetylation of allenic alcohols",
                  ["Procedure continues beyond the visually reviewed range; products, workup, purification, and individual yields are not claimed from page 6 alone.",
                   "The generic hydroxy-allenoate substrate is not expanded into a concrete SMILES."])
substrate = compound("benzyl 2-(hydroxy-4-aryl/heteroaryl)buta-2,3-dienoate generic", evidence, "allene_unspecified")
solv, pyr, ac = agent("DCM", evidence), agent("pyridine", evidence), agent("acetyl_chloride", evidence)
r["participants"] = [participant("p_allenol", substrate, "reactant", evidence, amount=q(10.0, "mmol"), equiv=q(1.0, "equiv"), ref=substrate),
                     participant("p_pyridine", pyr, "reagent", evidence, amount=q(15.0, "mmol"), equiv=q(1.5, "equiv"), ref=substrate),
                     participant("p_AcCl", ac, "reagent", evidence, amount=q(15.0, "mmol"), equiv=q(1.5, "equiv"), ref=substrate),
                     participant("p_solvent", solv, "solvent", evidence, volume=q(20, "mL"))]
r["reaction_scale"] = sf("generic allenol (10.0 mmol)", {"reference_molecule_id": substrate, "amount": q(10.0, "mmol")}, evidence)
set_solvent(r, solv, evidence, volume=q(20, "mL"))
add_step(r, "charge", evidence, ["p_allenol", "p_solvent"], "prepare solution", q(-5, "°C"), None)
add_step(r, "pyridine_hold", evidence, ["p_pyridine"], "add pyridine and stir", q(-5, "°C"), q(5.0, "min"))
add_step(r, "acetyl_chloride_addition", evidence, ["p_AcCl"], "slowly add acetyl chloride", q(-5, "°C"), q(5.0, "min"))
add_step(r, "post_addition_hold", evidence, ["p_allenol", "p_pyridine", "p_AcCl"], "stir after addition", None, q(45, "min"))
classification(r, evidence, (("C", "O", 1),))


# Source context: bibliography, prior-art schemes, proposed mechanism discussion, and analytical-only text.
bundle["source_context"].extend([
    {
        "context_id": f"{DOC}:bibliography",
        "kind": "bibliographic_metadata",
        "title": "Divergent Reactivity of δ- and β′-Acetoxy Allenoates with 2-Sulfonamidoindoles via Phosphine Catalysis: Entry to Dihydro-α-carboline, α-Carboline, and Spiro-cyclopentene Motifs",
        "doi": "10.1021/acs.joc.1c01137",
        "evidence_ids": E1,
    },
    {
        "context_id": f"{DOC}:prior_art_schemes_1_2ab",
        "kind": "literature_examples",
        "description": "Scheme 1 and Scheme 2a-b summarize cited prior work and are not counted as experiments in this paper.",
        "evidence_ids": E2,
    },
    {
        "context_id": f"{DOC}:proposed_mechanism_discussion",
        "kind": "proposed_mechanism",
        "description": "Pages 5-6 discuss proposed intermediates A1-A6, B3-B7 and C1-C3 and refer to Scheme 8/Scheme S1; these are hypotheses, not independent experiments. Scheme 8 artwork is not visible on the provided page images 1-6.",
        "evidence_ids": E5 + E6,
    },
    {
        "context_id": f"{DOC}:analytical_general_methods",
        "kind": "characterization_context",
        "description": "General NMR/IR/HRMS/X-ray methods and standalone characterization text are retained as source context, not synthetic experiments.",
        "evidence_ids": E6,
    },
])


# Coverage contains exactly the visually inspected PDF pages 1-6.
page_records = {i: [] for i in range(1, 7)}
for rxn in bundle["reactions"]:
    for page in range(1, 7):
        if ev(page, "image") in rxn["evidence_ids"] or ev(page, "text") in rxn["evidence_ids"]:
            page_records[page].append(rxn["reaction_id"])

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title, abstract, and graphical abstract", "finding": "Current-work overview of three divergent annulation families; recorded as generic schemes, not duplicated as specific experiments.", "evidence_ids": [ev(1, "image"), ev(1, "text")]},
            {"region": "introduction and Scheme 1a lead-in", "finding": "Background and cited prior work; no new specific experiment assigned on this page.", "evidence_ids": [ev(1, "image"), ev(1, "text")]},
        ], "record_ids": page_records[1], "unresolved": ["Graphical-abstract structures are family-level drawings; full source-verified SMILES and atom mappings were not assigned."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Scheme 1 and Scheme 2a-b", "finding": "Cited previous work; retained as literature context and not counted as current-paper experiments.", "evidence_ids": [ev(2, "image"), ev(2, "text")]},
            {"region": "Scheme 2c and optimization narrative", "finding": "Current-work generic overview plus prose cross-check for Table 1 entries, including degassed toluene for entry 2.", "evidence_ids": [ev(2, "image"), ev(2, "text")]},
        ], "record_ids": page_records[2], "unresolved": ["Narrative reference to earlier use of one equivalent PPh3 is not cleanly assigned to a unique Table 1 row; table footnote values are retained and the prose is not used to overwrite them."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Table 1, rows 1-23 and footnotes", "finding": "Twenty-three independent optimization experiments were transcribed row by row, including blank cells, qualitative Complex, bounds, and isolated-yield basis.", "evidence_ids": [ev(3, "image"), ev(3, "text")]},
            {"region": "scope discussion below Table 1", "finding": "Narrative used only to cross-check the Scheme 3/4 substrate-scope records on page 4.", "evidence_ids": [ev(3, "image"), ev(3, "text")]},
        ], "record_ids": page_records[3], "unresolved": ["Entry 18 phosphine amount is internally conflicting because the footnote assigns it to both 0.04 and 0.02 mmol ranges.", "Entries 10 and 12 have blank yield cells; blank is not converted to zero."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "Scheme 3", "finding": "Fifteen [3+3] substrate-scope experiments (3aa-3aj, 3ba-3da, 3ea-3fa) were recorded with isolated yields.", "evidence_ids": [ev(4, "image"), ev(4, "text")]},
            {"region": "Scheme 4", "finding": "Ten annulation/tosyl-migration scope experiments were recorded with isolated yields and X-ray notes where shown.", "evidence_ids": [ev(4, "image"), ev(4, "text")]},
            {"region": "Table 2, rows 1-12 and footnotes", "finding": "Twelve [4+1] optimization experiments were recorded, including trace outcomes, product-specific dr, and the entry-12 Cs2CO3 addition.", "evidence_ids": [ev(4, "image"), ev(4, "text")]},
        ], "record_ids": page_records[4], "unresolved": ["Table 2 footnote names 2a while the structure and narrative identify 5a; both source statements are preserved as a conflict.", "Table 2 does not report a method for dr determination."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "Scheme 5", "finding": "Thirteen spiro-annulation scope experiments were recorded with isolated yields and product-targeted dr; wedge/dashed relative stereochemistry was retained without inventing absolute configuration.", "evidence_ids": [ev(5, "image"), ev(5, "text")]},
            {"region": "Scheme 6", "finding": "Specific alpha-umpolung addition giving noncyclized 7ha (78%) was recorded separately.", "evidence_ids": [ev(5, "image"), ev(5, "text")]},
            {"region": "mechanistic discussion", "finding": "Proposed pathways are retained as source context, not counted as observed experiments.", "evidence_ids": [ev(5, "image"), ev(5, "text")]},
        ], "record_ids": page_records[5], "unresolved": ["Scheme 5 does not state the analytical method for dr.", "Absolute configurations of the spiro products are not assigned; source-verified isomeric SMILES are therefore left empty."]},
        {"page": 6, "state": "attempted", "regions": [
            {"region": "Scheme 7a-c", "finding": "Six 3aa-to-4aa condition combinations, two no-PPh3 temperature controls, and one TsCl control were recorded as mechanistic-study experiments.", "evidence_ids": [ev(6, "image"), ev(6, "text")]},
            {"region": "scale-up statement", "finding": "Three 1.0 mmol scale-up experiments referenced here were linked to explicit extracted-text procedures on PDF page 13; page 13 is context, not visually claimed coverage.", "evidence_ids": [ev(6, "image"), ev(6, "text"), ev(13, "text")]},
            {"region": "experimental syntheses i-ii", "finding": "Specific syntheses of precursors 1f and 1h were recorded with reaction, workup, purification, masses, amounts, and isolated yields.", "evidence_ids": [ev(6, "image"), ev(6, "text")]},
            {"region": "general procedure iii opening", "finding": "The page-6 portion of the general acetylation procedure for 5a-g was recorded as a partial generic scheme with separate -5 °C hold and addition steps.", "evidence_ids": [ev(6, "image"), ev(6, "text")]},
            {"region": "general analytical methods", "finding": "NMR/IR/HRMS/X-ray method text is characterization context, not a synthetic experiment.", "evidence_ids": [ev(6, "image"), ev(6, "text")]},
        ], "record_ids": page_records[6], "unresolved": ["Scheme 8 and Scheme 9 artwork are not visible on the supplied page-1-to-6 images; proposed mechanism is not converted into experimental records.", "The 5a-g general procedure continues on page 7, so the page-6 generic record remains partial."]},
    ],
}


inventory = {
    "document_id": DOC,
    "visually_inspected_pages": [1, 2, 3, 4, 5, 6],
    "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in bundle["reactions"]),
    "generic_schemes": sum(r["record_type"] == "generic_scheme" for r in bundle["reactions"]),
    "reaction_ids": [r["reaction_id"] for r in bundle["reactions"]],
    "molecule_ids": [m["molecule_id"] for m in bundle["molecules"]],
    "notes": [
        "Previous-work schemes and analytical-only text are source_context, not experiments.",
        "No SI was associated.",
        "All records are machine submissions and must remain needs_review after merge.",
    ],
}

(RUN / "coverage").mkdir(parents=True, exist_ok=True)
(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
(RUN / "coverage/pages_0001_0006.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/pages_0001_0006_extraction.json").write_text(json.dumps(inventory, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
                  "specific": inventory["specific_experiments"], "generic": inventory["generic_schemes"]}, ensure_ascii=False))
