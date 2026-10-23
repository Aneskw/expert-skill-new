#!/usr/bin/env python3
"""Build the source-bounded v3.1 submission, including visual review of pages 1--10."""
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
DATASET = RUN / "dataset.json"
DOC = "doc_40554c21bc96f78b"
TITLE = "Copper(I)-Catalyzed Boryl Substitution of 1-Trifluoromethyl Allenes for the Synthesis of 3-Boryl-Substituted 1,1-gem-Difluorodienes"


def ev(page: int, kind: str = "text") -> str:
    return f"{DOC}:page{page}:{kind}"


def q(number, unit, raw=None, qualifier=None):
    value = {"number": number, "unit": unit}
    if raw is not None:
        value["raw_text"] = raw
    if qualifier:
        value["qualifier"] = qualifier
    return value


def sf(raw=None, norm=None, evidence=(), method="original_extraction", derivation=None,
       status=None):
    return st.fact(raw, norm, evidence, method=method, derivation=derivation, status=status)


data = json.loads(DATASET.read_text())
data["documents"][0]["title"] = TITLE
tmp = DATASET.with_suffix(".json.tmp")
tmp.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n")
tmp.replace(DATASET)
bundle = st.empty(data["documents"])
molecules = {}


def molecule(key, label, name=None, evidence=(ev(1, "image"),), stereo="unresolved"):
    if key in molecules:
        occ = molecules[key]["occurrences"][0]
        for eid in evidence:
            if eid not in occ["evidence_ids"]:
                occ["evidence_ids"].append(eid)
            if eid not in occ["compound_label"]["evidence_ids"]:
                occ["compound_label"]["evidence_ids"].append(eid)
        return key
    m = st.new_molecule(key, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = sf(label, label, evidence)
    m["names"]["original"] = sf(name or label, name or label, evidence)
    m["names"]["english"] = sf(name or label, name or label, evidence)
    m["identity_note"] = sf(
        "Source label/name registered; complete structure was not independently verified.",
        "Source label/name registered; complete structure was not independently verified.", evidence,
        method="inference", derivation="Identity is limited to the printed label/name; no CAS or external lookup was used.")
    if stereo == "allene":
        m["stereochemistry"]["axial"] = sf(
            "Allene C=C=C axis shown or named; no R_a/S_a or P/M assignment is reported.",
            [{"type": "allene_axis", "configuration": "unspecified"}], evidence)
        m["stereochemistry"]["chirality_status"] = sf("unassigned allene axis", "unspecified", evidence)
        m["stereochemistry"]["geometry_status"] = sf("E/Z is not used to assign the allene axis", "not_applicable", evidence,
                                                        method="inference", derivation="Allene axial configuration is separate from alkene E/Z.")
        m["stereochemistry"]["representation_status"] = sf("axis not encoded", "not_encoded", evidence)
    elif stereo in {"E", "Z"}:
        m["stereochemistry"]["double_bonds"] = sf(
            f"Compound label/drawing specifies {stereo}",
            [{"configuration": stereo, "atom_locations": "unmapped source drawing"}], evidence)
        m["stereochemistry"]["chirality_status"] = sf("no atom/axial assignment reported", "unresolved", evidence)
        m["stereochemistry"]["geometry_status"] = sf(stereo, "specified", evidence)
        m["stereochemistry"]["representation_status"] = sf("geometry retained textually; no verified SMILES", "not_encoded", evidence)
    else:
        # new_molecule already supplies explicit unresolved status facts.
        pass
    molecules[key] = m
    return key


def ensure(key, label=None, name=None, evidence=(ev(1, "text"),), stereo="unresolved"):
    return molecule(key, label or key, name, evidence, stereo)


def part(r, key, role, evidence, stage="reaction", mass=None, amount=None, volume=None,
         equiv=None, molpct=None, ref=None):
    ensure(key, evidence=evidence)
    pid = "p_" + re.sub(r"[^A-Za-z0-9_]+", "_", key)
    if any(p["participant_id"] == pid for p in r["participants"]):
        pid += "_" + str(len(r["participants"]) + 1)
    p = st.new_participant(pid, key, role, evidence, stage)
    for field, value in (("mass", mass), ("amount", amount), ("volume", volume),
                         ("equivalents", equiv), ("mol_percent", molpct)):
        if value is not None:
            raw, norm = value
            p["charges"][field] = sf(raw, norm, evidence)
    if ref:
        p["equivalence_reference"] = sf(ref, ref, evidence, method="inference",
                                           derivation="The table/procedure defines equivalents or mol% relative to this substrate.")
    r["participants"].append(p)
    return p["participant_id"]


def step(r, sid, evidence, operation, temperature=None, time=None, participants=(), stage="reaction", atmosphere=None):
    s = st.new_step(sid, evidence, stage)
    s["participant_ids"] = list(participants)
    s["operation"] = sf(operation, operation, evidence)
    if temperature is not None:
        raw, norm = temperature
        s["temperature"] = sf(raw, norm, evidence)
    if time is not None:
        raw, norm = time
        s["time"] = sf(raw, norm, evidence)
    if atmosphere is not None:
        s["atmosphere"] = sf(atmosphere, atmosphere, evidence)
    r["steps"].append(s)
    return s


def outcome(r, oid, kind, targets, evidence, raw, norm, yield_type=None, basis=None, scope=None,
            status=None):
    o = st.new_outcome(oid, kind, targets, evidence)
    o["value"] = sf(raw, norm, evidence, status=status)
    if yield_type is not None:
        o["yield_type"] = sf(yield_type, yield_type, evidence)
    if basis is not None:
        if isinstance(basis, tuple):
            braw, bnorm, bevidence = basis
        else:
            braw = bnorm = basis
            bevidence = evidence
        o["measurement_basis"] = sf(braw, bnorm, bevidence)
    if scope is not None:
        o["scope"] = sf(scope, scope, evidence)
    r["outcomes"].append(o)
    return o


def classification(r, evidence, formed=(), increased=(), broken=(), decreased=(), coverage="partial"):
    changes = []
    labels = {}
    def add(a, b, before, after, kind):
        cid = f"b{len(changes)+1}"
        ch = {"change_id": cid, "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
              "before_order": before, "after_order": after, "change_type": kind}
        changes.append(ch)
        if kind in {"formed", "order_increased", "hydrogen_gain"}:
            lab = st.bond_label(a, b, after)
            labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(cid)
    for a, b, order in formed:
        add(a, b, 0, order, "formed")
    for a, b, before, after in increased:
        add(a, b, before, after, "order_increased")
    for a, b, order in broken:
        add(a, b, order, 0, "broken")
    for a, b, before, after in decreased:
        add(a, b, before, after, "order_decreased")
    deriv = "Partial net-bond inventory inferred from the printed reactant/product drawings; atom mapping and omitted coproducts are unresolved."
    r["classification"]["bond_changes"] = [sf(None, c, evidence, method="structure_derivation", derivation=deriv) for c in changes]
    r["classification"]["labels"] = sf(None, list(labels.values()), evidence, method="structure_derivation", derivation=deriv)
    r["classification"]["coverage"] = sf(None, coverage, evidence, method="inference", derivation=deriv)
    r["classification"]["limitations"] = [deriv, "No complete atom-mapped reaction SMILES was asserted."]


def no_classification(r, evidence, not_applicable=False):
    cov = "not_applicable" if not_applicable else "unresolved"
    deriv = "No observed transformation to classify." if not_applicable else "Source structures are insufficient for a defensible bond inventory."
    r["classification"]["labels"] = sf(None, [] if not_applicable else None, evidence,
                                            method="inference", derivation=deriv)
    r["classification"]["coverage"] = sf(None, cov, evidence, method="inference", derivation=deriv)
    r["classification"]["limitations"] = [deriv]


def reaction(suffix, evidence, record_type="specific_experiment", categories=("other",), location=None,
             name=None, transform=None):
    rid = f"{DOC}:{suffix}"
    r = st.new_reaction(rid, DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = sf(list(categories), list(categories), evidence)
    r["location"] = sf(location or suffix, location or suffix, evidence)
    if name:
        r["reaction_name"] = sf(name, name, evidence)
    if transform:
        r["transformation_type"] = sf(None, transform, evidence, method="inference",
                                           derivation="Normalized from the authors' title/caption and the printed reaction arrow.")
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
                   "unresolved": ["Machine record; complete molecular structures and atom mapping remain unverified."]}
    r["completeness"] = "partial"
    bundle["reactions"].append(r)
    return r


# Molecule registry: principal substrates/products and stereoisomers.
for x in "abcdefghijk":
    ensure(f"m_1{x}", f"1{x}", evidence=(ev(2, "image") if x != "k" else ev(5, "image"),), stereo="allene")
for x in "abcdef":
    ensure(f"m_3{x}", f"3{x}", evidence=(ev(2, "image"), ev(6, "image")))
for x in "ghij":
    for geom in ("Z", "E"):
        ensure(f"m_3{x}_{geom}", f"({geom})-3{x}", evidence=(ev(2, "image"), ev(6, "text")), stereo=geom)
for x in "abcdeghi":
    ensure(f"m_4{x}", f"4{x}", evidence=(ev(3, "image"),))
ensure("m_4f_E", "(E)-4f", evidence=(ev(3, "image"),), stereo="E")
ensure("m_4f_Z", "(Z)-4f", evidence=(ev(3, "image"),), stereo="Z")
ensure("m_2", "2", "bis(pinacolato)diboron", (ev(2, "image"),))
ensure("m_5", "5", "4-phenyl-1,2,4-triazoline-3,5-dione", (ev(3, "image"),))
ensure("m_6", "6", "difluoro-containing heterocyclic vinyl boronate 6", (ev(3, "image"),))


# Page 1: current-work generic transformations. Scheme 1a bioactive structures are context, not reactions.
g1 = reaction("scheme1_boryl_substitution", (ev(1, "image"),), "generic_scheme", (), "Scheme 1b / graphical abstract",
              "Copper(I)-catalyzed synthesis of 3-boryl-1,1-gem-difluorodienes", "defluoroborylation of CF3-substituted allenes")
classification(g1, (ev(1, "image"),), formed=(("C", "B", 1),), increased=(("C", "C", 1, 2),),
               broken=(("C", "F", 1),), decreased=(("C", "C", 2, 1),))
g2 = reaction("scheme1_cross_coupling", (ev(1, "image"),), "generic_scheme", (), "Scheme 1c cross-coupling arrow",
              "Cross-coupling utility of 3-boryl-1,1-gem-difluorodienes", "Suzuki-Miyaura cross-coupling")
classification(g2, (ev(1, "image"),), formed=(("C", "C", 1),), broken=(("C", "B", 1), ("C", "Br", 1)))
g3 = reaction("scheme1_diels_alder", (ev(1, "image"),), "generic_scheme", (), "Scheme 1c Diels-Alder arrow",
              "Diels-Alder utility of 3-boryl-1,1-gem-difluorodienes", "[4+2] cycloaddition")
classification(g3, (ev(1, "image"),), formed=(("C", "N", 1), ("C", "N", 1)),
               increased=(("C", "C", 1, 2),), decreased=(("C", "C", 2, 1), ("C", "C", 2, 1)))


# Page 2 Table 1: every optimization row is an independent experiment.
table1 = [
    (1, "m_PPh3", "triphenylphosphine (PPh3)", "m_NaOMe", "sodium methoxide", "m_THF", "THF", 92),
    (2, "m_Xantphos", "Xantphos", "m_NaOMe", "sodium methoxide", "m_THF", "THF", 87),
    (3, "m_dppe", "dppe", "m_NaOMe", "sodium methoxide", "m_THF", "THF", 89),
    (4, "m_dppf", "dppf", "m_NaOMe", "sodium methoxide", "m_THF", "THF", 89),
    (5, "m_IMesHCl", "IMes·HCl", "m_NaOMe", "sodium methoxide", "m_THF", "THF", "trace"),
    (6, "m_PPh3", "triphenylphosphine (PPh3)", "m_NaOtBu", "sodium tert-butoxide", "m_THF", "THF", 92),
    (7, "m_PPh3", "triphenylphosphine (PPh3)", "m_KOtBu", "potassium tert-butoxide", "m_THF", "THF", 87),
    (8, "m_PPh3", "triphenylphosphine (PPh3)", "m_LiOMe", "lithium methoxide", "m_THF", "THF", "trace"),
    (9, "m_PPh3", "triphenylphosphine (PPh3)", "m_NaOMe", "sodium methoxide", "m_DMI", "DMI", 95),
    (10, "m_PPh3", "triphenylphosphine (PPh3)", "m_NaOMe", "sodium methoxide", "m_toluene", "toluene", 78),
    (11, "m_PPh3", "triphenylphosphine (PPh3)", "m_NaOMe", "sodium methoxide", "m_THF", "THF", 27),
    (12, None, None, "m_NaOMe", "sodium methoxide", "m_THF", "THF", 69),
]
for entry, lig, ligname, base, basename, solv, solvname, yld in table1:
    evidence = (ev(2, "image"), ev(2, "text"))
    r = reaction(f"table1_entry_{entry:02d}", evidence, categories=("optimization",),
                 location=f"Table 1 entry {entry}", name="Optimization of the boryl substitution",
                 transform="defluoroborylation of CF3-substituted allene 1a")
    ensure("m_CuCl", "CuCl", "copper(I) chloride", evidence)
    ensure(base, name=basename, evidence=evidence)
    ensure(solv, name=solvname, evidence=evidence)
    ids = [part(r, "m_1a", "reactant", evidence, amount=("0.25 mmol", q(0.25, "mmol"))),
           part(r, "m_2", "reagent", evidence, amount=("0.375 mmol", q(0.375, "mmol")), equiv=("1.5 equiv", q(1.5, "equiv")), ref="m_1a"),
           part(r, "m_CuCl", "catalyst", evidence, amount=("0.0125 mmol", q(0.0125, "mmol")), molpct=("5 mol%", q(5, "mol%")), ref="m_1a")]
    if lig:
        ensure(lig, name=ligname, evidence=evidence)
        ids.append(part(r, lig, "ligand", evidence, amount=("0.0125 mmol", q(0.0125, "mmol")), molpct=("5 mol%", q(5, "mol%")), ref="m_1a"))
    if entry == 11:
        ids.append(part(r, base, "reagent", evidence, molpct=("10 mol%", q(10, "mol%")), ref="m_1a"))
    else:
        ids.append(part(r, base, "reagent", evidence, amount=("0.375 mmol", q(0.375, "mmol")), equiv=("1.5 equiv", q(1.5, "equiv")), ref="m_1a"))
    ids += [part(r, solv, "solvent", evidence, volume=("500 μL", q(500, "uL"))),
            part(r, "m_3a", "product", evidence)]
    r["reaction_scale"] = sf("1a (0.25 mmol)", {"reference_molecule_id": "m_1a", "amount": q(0.25, "mmol")}, evidence)
    r["solvent_conditions"]["composition"] = sf(solvname, [solvname], evidence)
    r["solvent_conditions"]["total_volume"] = sf("500 μL", q(500, "uL"), evidence)
    step(r, "reaction", evidence, "stir", ("30 °C", q(30, "°C")), ("2 h", q(2, "h")), ids)
    if yld == "trace":
        outcome(r, "product_trace", "trace", ["m_3a"], evidence, "trace", "trace",
                basis="19F NMR analysis of the crude mixture with fluorobenzene as internal standard")
    else:
        outcome(r, "yield_3a", "yield", ["m_3a"], evidence, f"{yld}%", q(yld, "%"), "NMR yield",
                "19F NMR analysis of the crude mixture with fluorobenzene as internal standard")
    classification(r, evidence, formed=(("C", "B", 1),), increased=(("C", "C", 1, 2),),
                   broken=(("C", "F", 1),), decreased=(("C", "C", 2, 1),))


# Page 2 Table 2, enriched by the individual paragraphs on pages 5--7.
scope = {
    "a": (3, 83, None, 0.50), "b": (2.5, 89, None, 0.49), "c": (8, 79, None, 0.25),
    "d": (4, 91, None, 0.50), "e": (5, 84, None, 0.50), "f": (25, 85, None, 0.50),
    "g": (5, 50, ("Z/E", 94, 6), 0.50), "h": (5.5, 49, ("Z/E", 69, 31), 0.50),
    "i": (6, 87, ("Z/E", 58, 42), 0.50), "j": (7, 45, ("Z/E", 79, 21), 0.50),
}
prod_amounts = {
    "a": (150.0, 0.415, 6), "b": (136.9, 0.44, 6), "c": (64.9, None, 6), "d": (144.9, 0.45, 6),
    "e": (146.3, 0.42, 6), "f": (158.1, 0.43, 6), "g": (80.1, 0.25, 6), "h": (73.8, 0.25, 6),
    "i": (145.7, 0.44, 6), "j": (90.4, 0.23, 7),
}
substrate_masses = {"j": 144.8}
for x, (hours, yld, sel, scale_amt) in scope.items():
    detail_page = 6 if x in "abcdefghi" else 7
    evidence = (ev(2, "image"), ev(2, "text"), ev(detail_page, "image"), ev(detail_page, "text"))
    r = reaction(f"scope_table2_3{x}", evidence, categories=("substrate_scope",),
                 location=f"Table 2, 1{x} to 3{x}; individual paragraph page {detail_page}",
                 name="General borylation procedure", transform="defluoroborylation of a CF3-substituted allene")
    ids = [part(r, f"m_1{x}", "reactant", evidence,
                mass=(f"{substrate_masses[x]} mg", q(substrate_masses[x], "mg")) if x in substrate_masses else None,
                amount=(f"{scale_amt} mmol", q(scale_amt, "mmol"))),
           part(r, "m_CuCl", "catalyst", evidence, amount=("0.025 mmol", q(0.025, "mmol")), molpct=("5 mol%", q(5, "mol%")), ref=f"m_1{x}"),
           part(r, "m_PPh3", "ligand", evidence, amount=("0.025 mmol (table general condition)", q(0.025, "mmol")), molpct=("5 mol%", q(5, "mol%")), ref=f"m_1{x}"),
           part(r, "m_2", "reagent", evidence, amount=("0.75 mmol", q(0.75, "mmol")), equiv=("1.5 equiv", q(1.5, "equiv")), ref=f"m_1{x}"),
           part(r, "m_NaOMe", "reagent", evidence, amount=("0.75 mmol", q(0.75, "mmol")), equiv=("1.5 equiv", q(1.5, "equiv")), ref=f"m_1{x}"),
           part(r, "m_THF", "solvent", evidence, volume=("1.0 mL", q(1.0, "mL")))]
    products = [f"m_3{x}"] if x in "abcdef" else [f"m_3{x}_Z", f"m_3{x}_E"]
    ids += [part(r, p, "product", evidence) for p in products]
    r["reaction_scale"] = sf(f"1{x} ({scale_amt} mmol)", {"reference_molecule_id": f"m_1{x}", "amount": q(scale_amt, "mmol")}, evidence)
    r["solvent_conditions"]["composition"] = sf("THF", ["THF"], evidence)
    r["solvent_conditions"]["total_volume"] = sf("1.0 mL", q(1.0, "mL"), evidence)
    if x == "j":
        alt1 = sf("7 h (Table 2)", q(7, "h"), (ev(2, "image"),))
        alt2 = sf("6.5 h (individual paragraph)", q(6.5, "h"), (ev(7, "image"), ev(7, "text")))
        s = step(r, "reaction", evidence, "stir at 30 °C", ("30 °C", q(30, "°C")), None, ids)
        s["time"] = sf(["7 h", "6.5 h"], None, evidence, status="conflict")
        s["time"]["alternatives"] = [alt1, alt2]
        step(r, "purification", (ev(7, "image"), ev(7, "text")), "purify product 3j by GPC", stage="purification")
        r["review"]["unresolved"].append("3j reaction time conflicts: Table 2 gives 7 h; the later individual paragraph gives 6.5 h.")
    else:
        step(r, "reaction", evidence, "stir at 30 °C", ("30 °C", q(30, "°C")), (f"{hours} h", q(hours, "h")), ids)
    r["procedure"]["general_procedure_reference"] = sf("General Borylation Procedure", "General Borylation Procedure", evidence)
    outcome(r, f"yield_3{x}", "yield", products, evidence, f"{yld}%", q(yld, "%"), "isolated yield", "isolated yield")
    if x in prod_amounts:
        mg, mmol, amount_page = prod_amounts[x]
        amount_evidence = (ev(amount_page, "image"), ev(amount_page, "text"))
        outcome(r, f"mass_3{x}", "product_amount", products, amount_evidence, f"{mg} mg", q(mg, "mg"), basis="weighed isolated product")
        if mmol is not None:
            outcome(r, f"amount_3{x}", "product_amount", products, amount_evidence, f"{mmol} mmol", q(mmol, "mmol"), basis="reported product amount")
        elif x == "c":
            o = outcome(r, "amount_3c_unreadable", "product_amount", products, (ev(6, "image"), ev(6, "text")),
                        "printed/extracted as '395 mmol' next to 64.9 mg and 79% from 0.25 mmol scale", None,
                        basis="reported product amount", status="unreadable")
            r["review"]["unresolved"].append("The 3c product amount is printed/extracted as '395 mmol'; it is chemically inconsistent and was not silently corrected.")
    if x == "a":
        pph3 = next(p for p in r["participants"] if p["molecule_id"] == "m_PPh3")
        pph3["charges"]["amount"] = sf(["0.025 mmol", "0.026 mmol"], None, evidence, status="conflict")
        pph3["charges"]["amount"]["alternatives"] = [sf("0.025 mmol (Table 2)", q(0.025, "mmol"), (ev(2, "image"),)), sf("0.026 mmol (3a paragraph)", q(0.026, "mmol"), (ev(5, "text"),))]
        r["review"]["unresolved"].append("PPh3 amount for 3a is 0.025 mmol in Table 2 and 0.026 mmol in the individual procedure.")
    if x == "b":
        substrate = next(p for p in r["participants"] if p["molecule_id"] == "m_1b")
        scale_alternatives = [
            sf("0.50 mmol (Table 2)", q(0.50, "mmol"), (ev(2, "image"),)),
            sf("100.8 mg (0.49 mmol; individual paragraph)", q(0.49, "mmol"), (ev(6, "text"),)),
        ]
        substrate["charges"]["amount"] = sf(["0.50 mmol", "0.49 mmol"], None, evidence, status="conflict")
        substrate["charges"]["amount"]["alternatives"] = scale_alternatives
        r["reaction_scale"] = sf(["0.50 mmol", "0.49 mmol"], None, evidence, status="conflict")
        r["reaction_scale"]["alternatives"] = scale_alternatives
        r["review"]["unresolved"].append("1b charge is 0.50 mmol in Table 2 but 0.49 mmol in the individual procedure; neither value was selected as authoritative.")
        r["notes"].append(sf("Table 2: 1b 0.50 mmol; individual paragraph: 100.8 mg (0.49 mmol).", "Source conflict retained", evidence))
    if x in {"a", "g", "j"}:
        r["notes"].append(sf(f"3{x} contained a small amount of byproducts.", f"3{x} contained a small amount of byproducts.", evidence))
    if sel:
        metric, major, minor = sel
        basis = "19F NMR analysis of the crude mixture; configuration assignment also references GC and 2D NOESY in Table 2 footnote c"
        outcome(r, f"selectivity_3{x}", "selectivity", products, evidence, f"{metric} = {major}/{minor}",
                {"metric": metric, "numerator": major, "denominator": minor}, basis=basis,
                scope={"metric": metric, "order": metric.split("/")})
    classification(r, evidence, formed=(("C", "B", 1),), increased=(("C", "C", 1, 2),),
                   broken=(("C", "F", 1),), decreased=(("C", "C", 2, 1),))

# Table 2 failed terminal allene case 1k.
evidence = (ev(2, "image"), ev(2, "text"))
r = reaction("scope_table2_1k_complex_mixture", evidence, categories=("substrate_scope", "control"),
             location="Table 2, 1k", name="Attempted borylation of terminal allene 1k",
             transform="attempted defluoroborylation")
ids = [part(r, "m_1k", "reactant", evidence, amount=("0.50 mmol", q(0.50, "mmol"))),
       part(r, "m_CuCl", "catalyst", evidence, molpct=("5 mol%", q(5, "mol%")), ref="m_1k"),
       part(r, "m_PPh3", "ligand", evidence, molpct=("5 mol%", q(5, "mol%")), ref="m_1k"),
       part(r, "m_2", "reagent", evidence, equiv=("1.5 equiv", q(1.5, "equiv")), ref="m_1k"),
       part(r, "m_NaOMe", "reagent", evidence, equiv=("1.5 equiv", q(1.5, "equiv")), ref="m_1k"),
       part(r, "m_THF", "solvent", evidence, volume=("1.0 mL", q(1.0, "mL")))]
step(r, "reaction", evidence, "stir", ("30 °C", q(30, "°C")), None, ids)
outcome(r, "complex_mixture", "decomposition", ["m_1k"], evidence, "complex mixture", "complex mixture", scope="No product identity or yield reported")
no_classification(r, evidence, True)


# Page 3 Table 3: nine cross-coupling experiments.
cross = {
    "a": ("m_3a", "m_4bromoanisole", "4-bromoanisole", 25, 75, False, "isolated yield"),
    "b": ("m_3b", "m_1bromonaphthalene", "1-bromonaphthalene", 18, 61, False, "isolated yield"),
    "c": ("m_3b", "m_1bromo2methoxybenzene", "1-bromo-2-methoxybenzene", 27, 50, False, "isolated yield"),
    "d": ("m_3b", "m_1bromo2CF3benzene", "1-bromo-2-(trifluoromethyl)benzene", 21, 62, True, "isolated yield"),
    "e": ("m_3b", "m_5bromo123trimethoxybenzene", "5-bromo-1,2,3-trimethoxybenzene", 20, 83, True, "isolated yield"),
    "f": ("m_3g_Z", "m_4bromoanisole", "4-bromoanisole", 20, 59, True, "isolated yield"),
    "g": ("m_3d", "m_4bromoanisole", "4-bromoanisole", 20, 54, True, "isolated yield"),
    "h": ("m_3f", "m_4bromoanisole", "4-bromoanisole", 6.5, 76, True, "isolated yield"),
    "i": ("m_3f", "m_1bromopyrene", "1-bromopyrene", 23, 55, True, "GC analysis"),
}
cross_details = {
    "a": (7, 36.2, 0.10, 25.8, 0.076),
    "b": (7, 31.4, 0.10, 19.1, 0.061),
    "c": (7, 31.3, 0.100, 14.8, 0.0506),
    "d": (7, 31.9, 0.102, 21.0, 0.064),
    "e": (7, 31.2, 0.10, 29.1, 0.083),
    "f": (7, 32.0, 0.10, 17.6, 0.059),
    "g": (8, 32.4, 0.10, 16.3, 0.054),
    "h": (8, 36.9, 0.0997, 26.4, 0.075),
    "i": (8, 37.6, 0.102, 24.9, 0.056),
}
for x, (substrate, arbr, arname, hours, yld, low_arbr, basis) in cross.items():
    detail_page, substrate_mg, substrate_mmol, product_mg, product_mmol = cross_details[x]
    detail_pages = (7, 8) if x == "f" else (detail_page,)
    detail_evidence = tuple(eid for page in detail_pages for eid in (ev(page, "image"), ev(page, "text")))
    evidence = (ev(3, "image"), ev(3, "text"), *detail_evidence)
    ensure(arbr, name=arname, evidence=evidence)
    products = [f"m_4{x}"] if x != "f" else ["m_4f_E", "m_4f_Z"]
    r = reaction(f"scope_table3_4{x}", evidence, categories=("substrate_scope", "derivatization"),
                 location=f"Table 3, product 4{x}", name="Cross-coupling reactions of 3 with aryl bromides",
                 transform="Suzuki-Miyaura cross-coupling")
    ar_equiv = 1.05 if low_arbr else 1.5
    input_ids = [part(r, substrate, "reactant", evidence, mass=(f"{substrate_mg} mg", q(substrate_mg, "mg")),
                      amount=(f"{substrate_mmol} mmol", q(substrate_mmol, "mmol"))),
           part(r, arbr, "reactant", evidence, amount=(("0.105 mmol" if low_arbr else "0.15 mmol"), q(0.105 if low_arbr else 0.15, "mmol")),
                equiv=(f"{ar_equiv} equiv", q(ar_equiv, "equiv")), ref=substrate),
           part(r, "m_PdOAc2", "catalyst", evidence, amount=("0.005 mmol", q(0.005, "mmol")), molpct=("5 mol%", q(5, "mol%")), ref=substrate),
           part(r, "m_SPhos", "ligand", evidence, amount=("0.010 mmol", q(0.010, "mmol")), molpct=("10 mol%", q(10, "mol%")), ref=substrate),
           part(r, "m_NaOH", "reagent", evidence, amount=("0.30 mmol", q(0.30, "mmol")), equiv=("3.0 equiv", q(3.0, "equiv")), ref=substrate),
           part(r, "m_THF", "solvent", evidence, volume=("210 μL", q(210, "uL")))]
    product_ids = [part(r, p, "product", evidence) for p in products]
    ids = input_ids + product_ids
    r["reaction_scale"] = sf(f"{substrate_mmol} mmol of {substrate}",
                              {"reference_molecule_id": substrate, "amount": q(substrate_mmol, "mmol")}, evidence)
    r["solvent_conditions"]["composition"] = sf("THF and aqueous NaOH", ["THF", "2.5 M aqueous NaOH"], evidence)
    if x == "a":
        next(p for p in r["participants"] if p["molecule_id"] == arbr)["charges"]["mass"] = sf("28.5 mg", q(28.5, "mg"), detail_evidence)
        next(p for p in r["participants"] if p["molecule_id"] == "m_PdOAc2")["charges"]["mass"] = sf("1.1 mg", q(1.1, "mg"), detail_evidence)
        next(p for p in r["participants"] if p["molecule_id"] == "m_SPhos")["charges"]["mass"] = sf("4.0 mg", q(4.0, "mg"), detail_evidence)
        next(p for p in r["participants"] if p["molecule_id"] == "m_NaOH")["charges"]["volume"] = sf("0.12 mL of 2.5 M aqueous solution", q(0.12, "mL"), detail_evidence)
        r["solvent_conditions"]["total_volume"] = sf(
            "THF 0.21 mL plus 0.12 mL aqueous NaOH; author does not state a combined total volume",
            {"THF": q(0.21, "mL"), "aqueous_NaOH": q(0.12, "mL")}, detail_evidence)
        step(r, "setup", detail_evidence, "place Pd(OAc)2, SPhos, and 3a in an oven-dried vial; evacuate and backfill with nitrogen three times",
             participants=[input_ids[0], input_ids[2], input_ids[3]], stage="reaction", atmosphere="nitrogen")
        step(r, "addition", detail_evidence, "add THF, 2.5 M aqueous NaOH, and 4-bromoanisole",
             participants=[input_ids[1], input_ids[4], input_ids[5]], stage="reaction", atmosphere="nitrogen")
    else:
        r["solvent_conditions"]["total_volume"] = sf(
            "THF 210 μL; row-specific aqueous volume is not stated in Table 3",
            None, evidence, status="unreadable")
    step(r, "reaction", evidence, "stir under the Table 3/cross-coupling conditions", ("60 °C", q(60, "°C")),
         (f"{hours} h", q(hours, "h")), input_ids, atmosphere="nitrogen" if x == "a" else None)
    if x == "a":
        step(r, "workup", detail_evidence, "extract with Et2O, dry the organic layer over MgSO4, filter, and remove volatiles by rotary evaporation", stage="workup")
        step(r, "purification", detail_evidence, "flash silica gel chromatography, Et2O/hexane 0:100–2:98", stage="purification")
    r["procedure"]["general_procedure_reference"] = sf(
        "Cross-Coupling Reaction Procedure for 4a and Table 3 conditions",
        "Cross-Coupling Reaction Procedure for 4a and Table 3 conditions", evidence)
    y = outcome(r, f"yield_4{x}", "yield", products, evidence, f"{yld}%", q(yld, "%"), basis, basis)
    outcome(r, f"mass_4{x}", "product_amount", products, detail_evidence, f"{product_mg} mg", q(product_mg, "mg"),
            basis="reported product amount in the individual product paragraph")
    outcome(r, f"amount_4{x}", "product_amount", products, detail_evidence, f"{product_mmol} mmol", q(product_mmol, "mmol"),
            basis="reported product amount in the individual product paragraph")
    if substrate_mmol != 0.10:
        table_alt = sf("0.10 mmol (Table 3)", q(0.10, "mmol"), (ev(3, "image"), ev(3, "text")))
        paragraph_alt = sf(f"{substrate_mg} mg ({substrate_mmol} mmol; individual paragraph)", q(substrate_mmol, "mmol"), detail_evidence)
        substrate_part = next(p for p in r["participants"] if p["molecule_id"] == substrate)
        substrate_part["charges"]["amount"] = sf(["0.10 mmol", f"{substrate_mmol} mmol"], None, evidence, status="conflict")
        substrate_part["charges"]["amount"]["alternatives"] = [table_alt, paragraph_alt]
        r["reaction_scale"] = sf(["0.10 mmol", f"{substrate_mmol} mmol"], None, evidence, status="conflict")
        r["reaction_scale"]["alternatives"] = [table_alt, paragraph_alt]
        r["notes"].append(sf(
            f"Table 3 gives 0.10 mmol of substrate; the 4{x} paragraph gives {substrate_mg} mg ({substrate_mmol} mmol).",
            "Source scale conflict retained", evidence))
        r["review"]["unresolved"].append(
            f"4{x} substrate scale conflicts between Table 3 (0.10 mmol) and the individual paragraph ({substrate_mmol} mmol).")
    if x == "i":
        gc_alt = sf("Determined by GC analysis (Table 3 footnote c)", "GC analysis", (ev(3, "image"), ev(3, "text")))
        amount_alt = sf("Product 4i was obtained as 24.9 mg (0.056 mmol) in the individual paragraph", "reported product amount after the general isolation workflow", detail_evidence)
        y["yield_type"] = sf(["GC analysis", "reported product amount"], None, evidence, status="conflict")
        y["yield_type"]["alternatives"] = [gc_alt, amount_alt]
        y["measurement_basis"] = sf(["Table 3 footnote c: GC analysis", "individual paragraph: 24.9 mg, 0.056 mmol"], None, evidence, status="conflict")
        y["measurement_basis"]["alternatives"] = [gc_alt, amount_alt]
        r["notes"].append(sf(
            "Table 3 marks the 55% value as determined by GC analysis, while the 4i paragraph also reports 24.9 mg (0.056 mmol), numerically 55% of the stated charge.",
            "Yield-basis conflict retained without selecting one source", evidence))
        r["review"]["unresolved"].append("4i yield basis is conflicting/ambiguous: Table 3 says GC analysis, whereas the individual paragraph reports a matching product mass and amount.")
    if x == "f":
        outcome(r, "selectivity_4f", "selectivity", products, evidence, "E/Z = 94/6",
                {"metric": "E/Z", "numerator": 94, "denominator": 6}, basis="19F NMR analysis of the crude mixture",
                scope={"metric": "E/Z", "order": ["E", "Z"]})
        r["notes"].append(sf("4f contained a small amount of byproducts.", "4f contained a small amount of byproducts.", detail_evidence))
    classification(r, evidence, formed=(("C", "C", 1),), broken=(("C", "B", 1), ("C", "Br", 1)))


# Page 3 Scheme 2, enriched by the full individual procedure on page 8, and the textual negative observation.
evidence = (ev(3, "image"), ev(3, "text"), ev(8, "image"), ev(8, "text"))
r = reaction("scheme2_3b_5_to_6", evidence, categories=("derivatization",), location="Scheme 2",
             name="Diels-Alder reaction of 3b with a dienophile", transform="[4+2] cycloaddition")
input_ids = [part(r, "m_3b", "reactant", evidence, mass=("31.0 mg", q(31.0, "mg")), amount=("0.10 mmol", q(0.10, "mmol"))),
             part(r, "m_5", "reactant", evidence, mass=("69.3 mg", q(69.3, "mg")), amount=("0.40 mmol", q(0.40, "mmol")),
                  equiv=("4.0 equiv", q(4.0, "equiv")), ref="m_3b"),
             part(r, "m_toluene", "solvent", evidence, volume=("0.48 mL", q(0.48, "mL")))]
part(r, "m_6", "product", evidence)
r["reaction_scale"] = sf("3b (31.0 mg, 0.10 mmol)", {"reference_molecule_id": "m_3b", "amount": q(0.10, "mmol")}, evidence)
r["solvent_conditions"]["composition"] = sf("toluene", ["toluene"], evidence)
r["solvent_conditions"]["total_volume"] = sf("0.48 mL", q(0.48, "mL"), evidence)
r["procedure"]["general_procedure_reference"] = sf(
    "Diels-Alder reaction was performed according to reported procedure 9c",
    "reported procedure 9c (external cited literature, not supplied SI)", (ev(8, "image"), ev(8, "text")))
step(r, "setup", (ev(8, "image"), ev(8, "text")),
     "place dienophile 5 in an oven-dried vial; evacuate and backfill with nitrogen three times",
     participants=[input_ids[1]], stage="reaction", atmosphere="nitrogen")
step(r, "addition", (ev(8, "image"), ev(8, "text")), "add toluene and 3b",
     participants=[input_ids[2], input_ids[0]], stage="reaction", atmosphere="nitrogen")
step(r, "reaction", evidence, "stir", ("room temperature", "room temperature"), ("21 h", q(21, "h")), input_ids, atmosphere="nitrogen")
step(r, "workup", (ev(8, "image"), ev(8, "text")), "pass the reaction mixture through a short silica gel column eluting with Et2O", stage="workup")
step(r, "purification", (ev(8, "image"), ev(8, "text")), "flash silica gel chromatography, EtOAc/hexane 0:100–15:85", stage="purification")
outcome(r, "yield_6", "yield", ["m_6"], evidence, "30%", q(30, "%"), "isolated yield", "isolated after flash chromatography")
outcome(r, "mass_6", "product_amount", ["m_6"], (ev(8, "image"), ev(8, "text")), "14.7 mg", q(14.7, "mg"), basis="reported isolated product amount")
outcome(r, "amount_6", "product_amount", ["m_6"], (ev(8, "image"), ev(8, "text")), "0.030 mmol", q(0.030, "mmol"), basis="reported isolated product amount")
classification(r, evidence, formed=(("C", "N", 1), ("C", "N", 1)), increased=(("C", "C", 1, 2),),
               decreased=(("C", "C", 2, 1), ("C", "C", 2, 1)))

ensure("m_maleic_anhydride", "maleic anhydride", "maleic anhydride", evidence)
r = reaction("maleic_anhydride_unreactive", (ev(3, "text"),), categories=("control",),
             location="Page 3 discussion below Table 3", name="Attempted Diels-Alder reaction with maleic anhydride",
             transform="attempted [4+2] cycloaddition")
ids = [part(r, "m_3b", "reactant", (ev(3, "text"),)), part(r, "m_maleic_anhydride", "reactant", (ev(3, "text"),))]
step(r, "reaction_conditions_unreported", (ev(3, "text"),), "reaction attempted; conditions not reported", participants=ids)
outcome(r, "no_reaction", "no_reaction", ["m_3b", "m_maleic_anhydride"], (ev(3, "text"),), "unreactive", "unreactive", scope="No product or detection limit reported")
no_classification(r, (ev(3, "text"),), True)

r = reaction("figure1_proposed_mechanism", (ev(3, "image"), ev(3, "text")), "proposed_mechanism", (),
             "Figure 1a–b", "Proposed reaction mechanism for copper(I)-catalyzed borylation", "proposed catalytic cycle")
no_classification(r, (ev(3, "image"),), True)
r["notes"].append(sf("Authors propose Cu-alkoxide A, borylcopper B, sigma-allyl intermediates C/D, beta-fluoro elimination, and eta1–eta3 interconversion.",
                         "Author-proposed pathway, not a directly observed synthesis experiment.", (ev(3, "image"), ev(3, "text"))))


# Pages 4--5: Procedure A for 1a as three isolated transformations.
common4 = (ev(4, "text"), ev(4, "image"))
for key, label, name in [
    ("m_diisopropylamine", "diisopropylamine", "diisopropylamine"), ("m_nBuLi", "n-BuLi", "n-butyllithium"),
    ("m_2bromo333trifluoropropene", "2-bromo-3,3,3-trifluoropropene", "2-bromo-3,3,3-trifluoropropene"),
    ("m_4phenylcyclohexanone", "4-phenylcyclohexanone", "4-phenylcyclohexanone"),
    ("m_1a_alcohol", "alcohol intermediate toward 1a", "unnamed alcohol intermediate toward 1a"),
    ("m_Et3N", "Et3N", "triethylamine"), ("m_MsCl", "MsCl", "methanesulfonyl chloride"),
    ("m_DMAP", "DMAP", "4-dimethylaminopyridine"), ("m_CH2Cl2", "CH2Cl2", "dichloromethane"),
    ("m_1a_mesylate", "methanesulfonyl compound toward 1a", "unnamed methanesulfonyl intermediate toward 1a"),
    ("m_Et2Zn", "Et2Zn", "diethylzinc"), ("m_PdPPh3_4", "Pd(PPh3)4", "tetrakis(triphenylphosphine)palladium(0)"),
    ("m_EtOAc", "EtOAc", "ethyl acetate"), ("m_Et2O", "Et2O", "diethyl ether"),
    ("m_MgSO4", "MgSO4", "magnesium sulfate"), ("m_NaCl", "NaCl", "sodium chloride"),
]: ensure(key, label, name, common4)

r = reaction("prep1a_alcohol", common4, categories=("other",), location="Page 4 Procedure A, first operation",
             name="Preparation of the alcohol intermediate toward 1a", transform="carbonyl addition")
ids = [part(r, "m_diisopropylamine", "reagent", common4, volume=("5.6 mL", q(5.6, "mL")), amount=("40.0 mmol", q(40.0, "mmol"))),
       part(r, "m_nBuLi", "reagent", common4, volume=("25 mL", q(25, "mL")), amount=("40.0 mmol", q(40.0, "mmol"))),
       part(r, "m_2bromo333trifluoropropene", "reactant", common4, volume=("2.1 mL", q(2.1, "mL")), amount=("20.0 mmol", q(20.0, "mmol"))),
       part(r, "m_4phenylcyclohexanone", "reactant", common4, mass=("3.49 g", q(3.49, "g")), amount=("20.0 mmol", q(20.0, "mmol"))),
       part(r, "m_THF", "solvent", common4, volume=("30 mL + 10 mL + 10 mL", q(50, "mL"))),
       part(r, "m_1a_alcohol", "product", common4)]
step(r, "generate_LDA", common4, "add n-BuLi dropwise to diisopropylamine in THF and stir", ("0 °C", q(0, "°C")), ("1 h", q(1, "h")), ids, atmosphere="nitrogen")
step(r, "add_bromotrifluoropropene", common4, "cool and slowly add 2-bromo-3,3,3-trifluoropropene solution", ("−78 °C", q(-78, "°C")), ("1 h", q(1, "h")), ids, atmosphere="nitrogen")
step(r, "add_ketone", common4, "add 4-phenylcyclohexanone solution and stir", ("−78 °C", q(-78, "°C")), ("2.5 h", q(2.5, "h")), ids, atmosphere="nitrogen")
step(r, "workup", common4, "quench with 1 M HCl at 0 °C; extract with EtOAc three times; brine wash; MgSO4 dry; filter; evaporate", stage="workup")
step(r, "purification", common4, "flash silica gel chromatography, EtOAc/hexane 0:100–11:89", stage="purification")
outcome(r, "yield_alcohol", "yield", ["m_1a_alcohol"], common4, "67%", q(67, "%"), "isolated yield",
        ("product obtained after chromatography", "isolated after chromatography", common4))
outcome(r, "amount_alcohol", "product_amount", ["m_1a_alcohol"], common4, "3.57 g, 13.3 mmol", {"mass": q(3.57, "g"), "amount": q(13.3, "mmol")}, basis="reported isolated product amount")
classification(r, common4, formed=(("C", "C", 1),))

r = reaction("prep1a_mesylate", common4, categories=("other",), location="Page 4 Procedure A, second operation",
             name="Mesylation of the alcohol intermediate toward 1a", transform="mesylation")
ids = [part(r, "m_1a_alcohol", "reactant", common4, mass=("3.57 g", q(3.57, "g")), amount=("13.3 mmol", q(13.3, "mmol"))),
       part(r, "m_Et3N", "reagent", common4, volume=("2.2 mL", q(2.2, "mL")), amount=("16.0 mmol", q(16.0, "mmol"))),
       part(r, "m_MsCl", "reagent", common4, volume=("1.1 mL", q(1.1, "mL")), amount=("14.7 mmol", q(14.7, "mmol"))),
       part(r, "m_DMAP", "catalyst", common4, mass=("81.8 mg", q(81.8, "mg")), amount=("0.67 mmol", q(0.67, "mmol"))),
       part(r, "m_CH2Cl2", "solvent", common4, volume=("53 mL", q(53, "mL"))),
       part(r, "m_1a_mesylate", "product", common4)]
step(r, "addition", common4, "add Et3N dropwise to alcohol/MsCl/DMAP solution", ("0 °C", q(0, "°C")), None, ids, atmosphere="nitrogen")
step(r, "hold", common4, "stir", ("room temperature", "room temperature"), ("16 h", q(16, "h")), ids, atmosphere="nitrogen")
step(r, "workup", common4, "quench with saturated NH4Cl at 0 °C; extract with CH2Cl2 three times; brine wash; MgSO4 dry; filter; evaporate", stage="workup")
step(r, "purification", common4, "flash silica gel chromatography, EtOAc/hexane 0:100–15:85", stage="purification")
outcome(r, "yield_mesylate", "yield", ["m_1a_mesylate"], common4, "61%", q(61, "%"), "isolated yield", ("product obtained after chromatography", "isolated after chromatography", common4))
outcome(r, "amount_mesylate", "product_amount", ["m_1a_mesylate"], common4, "2.79 g, 8.1 mmol", {"mass": q(2.79, "g"), "amount": q(8.1, "mmol")}, basis="reported isolated product amount")
classification(r, common4, formed=(("O", "S", 1),))

r = reaction("prep1a_allene", common4, categories=("other",), location="Page 4 Procedure A, third operation",
             name="Preparation of CF3-substituted allene 1a", transform="palladium-catalyzed allene formation")
ids = [part(r, "m_1a_mesylate", "reactant", common4, mass=("2.79 g", q(2.79, "g")), amount=("8.1 mmol", q(8.1, "mmol"))),
       part(r, "m_Et2Zn", "reagent", common4, volume=("24.2 mL", q(24.2, "mL")), amount=("24.2 mmol", q(24.2, "mmol"))),
       part(r, "m_PdPPh3_4", "catalyst", common4, mass=("468.5 mg", q(468.5, "mg")), amount=("0.41 mmol", q(0.41, "mmol"))),
       part(r, "m_THF", "solvent", common4, volume=("82 mL", q(82, "mL"))), part(r, "m_1a", "product", common4)]
step(r, "addition", common4, "add Et2Zn dropwise to mesylate/Pd(PPh3)4 solution", ("0 °C", q(0, "°C")), None, ids, atmosphere="nitrogen")
step(r, "hold", common4, "stir", ("room temperature", "room temperature"), ("2.5 h", q(2.5, "h")), ids, atmosphere="nitrogen")
step(r, "workup", common4, "quench with saturated NH4Cl at 0 °C; filter through Celite; extract with EtOAc three times; brine wash; MgSO4 dry; filter; evaporate", stage="workup")
step(r, "purification", common4, "flash silica gel chromatography, hexane", stage="purification")
outcome(r, "yield_1a", "yield", ["m_1a"], common4, "45%", q(45, "%"), "isolated yield", ("1a obtained after chromatography", "isolated after chromatography", common4))
outcome(r, "amount_1a", "product_amount", ["m_1a"], common4, "924.5 mg, 3.7 mmol", {"mass": q(924.5, "mg"), "amount": q(3.7, "mmol")}, basis="reported isolated product amount")
classification(r, common4, formed=(("C", "C", 1),), increased=(("C", "C", 1, 2),), broken=(("O", "S", 1),))


# 1b--1i: source explicitly invokes Procedure A; product-specific yields are not printed on pages 4--5.
prep_overall = {
    "b": (4, "corresponding ketone for 1b", "Et2Zn"), "c": (4, "corresponding ketone for 1c", "Me3SiCH2ZnCl (0.49 M in THF)"),
    "d": (4, "corresponding ketone for 1d", "PhZnCl (0.29 M in THF)"), "e": (4, "corresponding ketone for 1e", "p-aryl zinc chloride (0.50 M in THF)"),
    "f": (4, "corresponding ketone for 1f", "(1-naphthyl)ZnCl (0.50 M in THF)"),
    "g": (5, "corresponding aldehyde for 1g", "Procedure A organozinc reagent"), "h": (5, "corresponding aldehyde for 1h", "Procedure A organozinc reagent"),
    "i": (5, "corresponding ketone for 1i", "Procedure A organozinc reagent"),
}
for x, (page, precursor_name, organozinc_name) in prep_overall.items():
    evidence = tuple(dict.fromkeys((ev(page, "text"), ev(page, "image"), ev(4, "text"))))
    precursor = f"m_precursor_1{x}"
    organozinc = f"m_organozinc_1{x}"
    ensure(precursor, precursor_name, precursor_name, evidence)
    ensure(organozinc, organozinc_name, organozinc_name, evidence)
    r = reaction(f"prep1{x}_procedure_A", evidence, categories=("other",), location=f"Page {page}, preparation of 1{x}",
                 name=f"Preparation of CF3-substituted allene 1{x} according to Procedure A", transform="multistep allene preparation")
    ids = [part(r, precursor, "reactant", evidence), part(r, "m_nBuLi", "reagent", evidence),
           part(r, "m_diisopropylamine", "reagent", evidence), part(r, "m_2bromo333trifluoropropene", "reactant", evidence),
           part(r, "m_MsCl", "reagent", evidence), part(r, "m_Et3N", "reagent", evidence), part(r, "m_DMAP", "catalyst", evidence),
           part(r, organozinc, "reagent", evidence), part(r, "m_PdPPh3_4", "catalyst", evidence),
           part(r, "m_THF", "solvent", evidence), part(r, "m_CH2Cl2", "solvent", evidence), part(r, f"m_1{x}", "product", evidence)]
    step(r, "procedure_A_stage1", evidence, "Procedure A: LDA generation, trifluoropropenyl reagent formation, and carbonyl addition", participants=ids)
    step(r, "procedure_A_stage2", evidence, "Procedure A: mesylation of the alcohol intermediate", participants=ids)
    step(r, "procedure_A_stage3", evidence, f"Procedure A: palladium-catalyzed substitution using {organozinc_name}", participants=ids)
    r["procedure"]["general_procedure_reference"] = sf("Procedure A described for 1a", "Procedure A", evidence)
    r["review"]["unresolved"].append("Product-specific scale and isolated yield are not reported in the inspected paragraph; exact substitutions to Procedure A beyond those printed remain unreported.")
    classification(r, evidence, formed=(("C", "C", 1),), increased=(("C", "C", 1, 2),), coverage="partial")


# Page 5: three-step preparation of 1j.
e5 = (ev(5, "text"), ev(5, "image"))
for key, label, name in [
    ("m_alkyne_1j", "alkyne", "unnamed alkyne used for 1j"), ("m_3phenylpropionaldehyde", "3-phenylpropionaldehyde", "3-phenylpropionaldehyde"),
    ("m_propargyl_alcohol_1j", "propargyl alcohol", "unnamed propargyl alcohol toward 1j"), ("m_SOCl2", "SOCl2", "thionyl chloride"),
    ("m_pyridine", "pyridine", "pyridine"), ("m_propargyl_chloride_1j", "propargyl chloride", "unnamed propargyl chloride toward 1j"),
    ("m_CF3SiMe3", "CF3SiMe3", "trimethyl(trifluoromethyl)silane"), ("m_CuTC", "copper(I) 2-thiophenecarboxylate", "copper(I) 2-thiophenecarboxylate"),
    ("m_KF", "KF", "potassium fluoride"),
]: ensure(key, label, name, e5)

r = reaction("prep1j_propargyl_alcohol", e5, categories=("other",), location="Page 5, preparation of 1j step 1", transform="carbonyl addition")
ids = [part(r, "m_nBuLi", "reagent", e5, volume=("13.5 mL", q(13.5, "mL")), amount=("21.0 mmol", q(21.0, "mmol"))),
       part(r, "m_alkyne_1j", "reactant", e5, volume=("2.41 mL", q(2.41, "mL")), amount=("22.0 mmol", q(22.0, "mmol"))),
       part(r, "m_3phenylpropionaldehyde", "reactant", e5, volume=("2.63 mL", q(2.63, "mL")), amount=("20.0 mmol", q(20.0, "mmol"))),
       part(r, "m_THF", "solvent", e5, volume=("80 mL", q(80, "mL"))), part(r, "m_propargyl_alcohol_1j", "product", e5)]
step(r, "lithiation", e5, "add n-BuLi dropwise to alkyne solution and stir", ("−78 °C", q(-78, "°C")), ("15 min", q(15, "min")), ids, atmosphere="nitrogen")
step(r, "aldehyde_addition", e5, "add 3-phenylpropionaldehyde; warm gradually to 0 °C and stir", ("−78 to 0 °C", {"minimum": -78, "maximum": 0, "unit": "°C", "raw_text": "−78 °C to 0 °C"}), ("1 h after reaching 0 °C", q(1, "h")), ids)
step(r, "workup", e5, "quench with 1 M HCl at 0 °C; extract with Et2O three times; brine wash; MgSO4 dry; filter; evaporate", stage="workup")
step(r, "purification", e5, "flash chromatography, EtOAc/hexane 0:100–15:85", stage="purification")
outcome(r, "yield_propargyl_alcohol", "yield", ["m_propargyl_alcohol_1j"], e5, "97%", q(97, "%"), "isolated yield", "isolated after chromatography")
outcome(r, "amount_propargyl_alcohol", "product_amount", ["m_propargyl_alcohol_1j"], e5, "4.60 g, 19.5 mmol", {"mass": q(4.60, "g"), "amount": q(19.5, "mmol")}, basis="reported isolated amount")
classification(r, e5, formed=(("C", "C", 1),))

r = reaction("prep1j_propargyl_chloride", e5, categories=("other",), location="Page 5, preparation of 1j step 2", transform="alcohol chlorination")
ids = [part(r, "m_propargyl_alcohol_1j", "reactant", e5, mass=("4.60 g", q(4.60, "g")), amount=("19.5 mmol", q(19.5, "mmol"))),
       part(r, "m_SOCl2", "reagent", e5, volume=("1.6 mL", q(1.6, "mL")), amount=("21.5 mmol", q(21.5, "mmol"))),
       part(r, "m_pyridine", "reagent", e5, volume=("1.9 mL", q(1.9, "mL")), amount=("23.4 mmol", q(23.4, "mmol"))),
       part(r, "m_CH2Cl2", "solvent", e5, volume=("30 mL", q(30, "mL"))), part(r, "m_propargyl_chloride_1j", "product", e5)]
step(r, "addition", e5, "add SOCl2 dropwise to alcohol/pyridine solution", ("0 °C", q(0, "°C")), ("10 min", q(10, "min")), ids, atmosphere="nitrogen")
step(r, "hold", e5, "warm to room temperature and stir overnight", ("room temperature", "room temperature"), ("overnight", "overnight"), ids)
step(r, "workup", e5, "dilute with Et2O; wash with 1 M HCl; extract water layer with Et2O; NaHCO3 wash; MgSO4 dry; filter; evaporate", stage="workup")
step(r, "purification", e5, "flash chromatography, EtOAc/hexane 0:100–3:97", stage="purification")
outcome(r, "yield_propargyl_chloride", "yield", ["m_propargyl_chloride_1j"], e5, "38%", q(38, "%"), "isolated yield", "isolated after chromatography")
outcome(r, "amount_propargyl_chloride", "product_amount", ["m_propargyl_chloride_1j"], e5, "1.90 g, 7.46 mmol", {"mass": q(1.90, "g"), "amount": q(7.46, "mmol")}, basis="reported isolated amount")
classification(r, e5, formed=(("C", "Cl", 1),), broken=(("C", "O", 1),))

r = reaction("prep1j_allene", e5, categories=("other",), location="Page 5, preparation of 1j step 3", transform="trifluoromethylative allene formation")
ids = [part(r, "m_propargyl_chloride_1j", "reactant", e5, mass=("1.90 g", q(1.90, "g")), amount=("7.46 mmol", q(7.46, "mmol"))),
       part(r, "m_CF3SiMe3", "reagent", e5, volume=("1.67 mL", q(1.67, "mL")), amount=("11.3 mmol", q(11.3, "mmol"))),
       part(r, "m_CuTC", "catalyst", e5, mass=("71.5 mg", q(71.5, "mg")), amount=("0.375 mmol", q(0.375, "mmol"))),
       part(r, "m_KF", "reagent", e5, mass=("658.6 mg", q(658.6, "mg")), amount=("11.3 mmol", q(11.3, "mmol"))),
       part(r, "m_THF", "solvent", e5, volume=("45 mL", q(45, "mL"))), part(r, "m_1j", "product", e5)]
step(r, "reaction", e5, "add propargyl chloride and CF3SiMe3 dropwise to CuTC/KF in THF and stir", ("60 °C", q(60, "°C")), ("23 h", q(23, "h")), ids, atmosphere="nitrogen")
step(r, "workup", e5, "pour into water; extract with Et2O three times; brine wash; filter; evaporate", stage="workup")
step(r, "purification", e5, "flash chromatography with hexane, then GPC", stage="purification")
outcome(r, "yield_1j", "yield", ["m_1j"], e5, "9%", q(9, "%"), "isolated yield", "isolated after chromatography and GPC")
outcome(r, "amount_1j", "product_amount", ["m_1j"], e5, "201.5 mg, 0.70 mmol", {"mass": q(201.5, "mg"), "amount": q(0.70, "mmol")}, basis="reported isolated amount")
classification(r, e5, formed=(("C", "C", 1),), increased=(("C", "C", 1, 2),), broken=(("C", "Cl", 1),))


# Page 5: three-step preparation of 1k.
for key, label, name in [
    ("m_CBr4", "CBr4", "carbon tetrabromide"), ("m_trifluoroacetophenone", "2,2,2-trifluoroacetophenone", "2,2,2-trifluoroacetophenone"),
    ("m_gemdibromo_1k", "gem-dibromo compound", "unnamed gem-dibromo intermediate toward 1k"),
    ("m_paraformaldehyde", "paraformaldehyde", "paraformaldehyde"), ("m_mesylate_1k", "methanesulfonyl compound", "unnamed methanesulfonyl intermediate toward 1k"),
    ("m_LiBr", "LiBr", "lithium bromide"), ("m_DMF", "dimethylformamide", "dimethylformamide"), ("m_Zn", "Zn", "zinc powder"),
]: ensure(key, label, name, e5)

r = reaction("prep1k_gem_dibromo", e5, categories=("other",), location="Page 5, preparation of 1k step 1", transform="Corey-Fuchs-type gem-dibromoolefination")
ids = [part(r, "m_CBr4", "reagent", e5, mass=("3.78 g", q(3.78, "g")), amount=("11.4 mmol", q(11.4, "mmol"))),
       part(r, "m_PPh3", "reagent", e5, mass=("5.80 g", q(5.80, "g")), amount=("22.1 mmol", q(22.1, "mmol"))),
       part(r, "m_trifluoroacetophenone", "reactant", e5, volume=("1.36 mL", q(1.36, "mL")), amount=("10.0 mmol", q(10.0, "mmol"))),
       part(r, "m_toluene", "solvent", e5, volume=("28 mL", q(28, "mL"))), part(r, "m_gemdibromo_1k", "product", e5)]
step(r, "reagent_preformation", e5, "add CBr4 quickly to PPh3 solution and stir", ("room temperature", "room temperature"), ("30 min", q(30, "min")), ids, atmosphere="nitrogen")
step(r, "substrate_addition", e5, "add trifluoroacetophenone dropwise over 15 min and stir", ("room temperature", "room temperature"), ("30 min after addition", q(30, "min")), ids)
step(r, "reflux", e5, "reflux", ("reflux", "reflux"), ("18 h", q(18, "h")), ids)
step(r, "workup", e5, "cool; add hexane; filter salts through Celite; quench filtrate with water; extract with hexane three times; brine wash; MgSO4 dry; filter; evaporate", stage="workup")
step(r, "purification", e5, "flash chromatography with hexane", stage="purification")
outcome(r, "yield_gemdibromo", "yield", ["m_gemdibromo_1k"], e5, "76%", q(76, "%"), "isolated yield", "isolated after chromatography")
outcome(r, "amount_gemdibromo", "product_amount", ["m_gemdibromo_1k"], e5, "2.51 g, 7.61 mmol", {"mass": q(2.51, "g"), "amount": q(7.61, "mmol")}, basis="reported isolated amount")
classification(r, e5, formed=(("C", "Br", 1), ("C", "Br", 1)), increased=(("C", "C", 1, 2),), decreased=(("C", "O", 2, 1),))

r = reaction("prep1k_mesylate", e5, categories=("other",), location="Page 5, preparation of 1k step 2", transform="homologation and mesylation")
ids = [part(r, "m_gemdibromo_1k", "reactant", e5, mass=("2.51 g", q(2.51, "g")), amount=("7.61 mmol", q(7.61, "mmol"))),
       part(r, "m_nBuLi", "reagent", e5, volume=("4.76 mL", q(4.76, "mL")), amount=("7.61 mmol", q(7.61, "mmol"))),
       part(r, "m_paraformaldehyde", "reactant", e5, mass=("933.4 mg", q(933.4, "mg")), amount=("31.1 mmol", q(31.1, "mmol"))),
       part(r, "m_Et3N", "reagent", e5, volume=("1.6 mL", q(1.6, "mL")), amount=("11.4 mmol", q(11.4, "mmol"))),
       part(r, "m_MsCl", "reagent", e5, volume=("1.2 mL", q(1.2, "mL")), amount=("15.2 mmol", q(15.2, "mmol"))),
       part(r, "m_THF", "solvent", e5, volume=("76 mL", q(76, "mL"))), part(r, "m_mesylate_1k", "product", e5)]
step(r, "lithiation", e5, "add n-BuLi dropwise to gem-dibromo compound solution and stir", ("−78 °C", q(-78, "°C")), ("40 min", q(40, "min")), ids, atmosphere="nitrogen")
step(r, "formaldehyde_addition", e5, "add paraformaldehyde; warm to room temperature with stirring", ("−78 °C to room temperature", {"minimum": -78, "maximum": 25, "unit": "°C", "raw_text": "−78 °C to room temperature"}), ("17 h", q(17, "h")), ids)
step(r, "base_addition", e5, "cool and add Et3N", ("0 °C", q(0, "°C")), ("30 min", q(30, "min")), ids)
step(r, "mesylation", e5, "add MsCl dropwise and stir", ("0 °C", q(0, "°C")), ("2.5 h", q(2.5, "h")), ids)
step(r, "workup", e5, "quench with 1 M HCl at 0 °C; extract with Et2O three times; brine wash; MgSO4 dry; filter; evaporate", stage="workup")
step(r, "purification", e5, "flash chromatography, EtOAc/hexane 3:97–20:80", stage="purification")
outcome(r, "yield_mesylate_1k", "yield", ["m_mesylate_1k"], e5, "60%", q(60, "%"), "isolated yield", "isolated after chromatography")
outcome(r, "amount_mesylate_1k", "product_amount", ["m_mesylate_1k"], e5, "1.65 g, 4.59 mmol", {"mass": q(1.65, "g"), "amount": q(4.59, "mmol")}, basis="reported isolated amount")
classification(r, e5, formed=(("C", "C", 1), ("O", "S", 1)), broken=(("C", "Br", 1),))

r = reaction("prep1k_allene", e5, categories=("other",), location="Page 5, preparation of 1k step 3", transform="allene formation")
ids = [part(r, "m_mesylate_1k", "reactant", e5, mass=("1.65 g", q(1.65, "g")), amount=("4.59 mmol", q(4.59, "mmol"))),
       part(r, "m_LiBr", "reagent", e5, mass=("397.2 mg", q(397.2, "mg")), amount=("4.57 mmol", q(4.57, "mmol"))),
       part(r, "m_DMF", "solvent", e5, volume=("9.2 mL", q(9.2, "mL"))),
       part(r, "m_Zn", "reagent", e5, mass=("334.8 mg", q(334.8, "mg")), amount=("5.12 mmol", q(5.12, "mmol"))),
       part(r, "m_1k", "product", e5)]
step(r, "bromide_stage", e5, "stir mesylate and LiBr in DMF", ("50 °C", q(50, "°C")), ("5 h", q(5, "h")), ids, atmosphere="nitrogen")
step(r, "zinc_stage", e5, "cool to room temperature, add zinc powder, and stir", ("room temperature", "room temperature"), ("22 h", q(22, "h")), ids)
step(r, "workup", e5, "quench with 1 M HCl at 0 °C; extract with Et2O three times; brine wash; MgSO4 dry; filter; evaporate", stage="workup")
step(r, "purification", e5, "flash chromatography with pentane, then Kugelrohr distillation", stage="purification")
outcome(r, "yield_1k", "yield", ["m_1k"], e5, "47%", q(47, "%"), "isolated yield", "isolated after chromatography and Kugelrohr distillation")
outcome(r, "amount_1k", "product_amount", ["m_1k"], e5, "396.5 mg, 2.15 mmol", {"mass": q(396.5, "mg"), "amount": q(2.15, "mmol")}, basis="reported isolated amount")
classification(r, e5, increased=(("C", "C", 1, 2),), broken=(("O", "S", 1), ("C", "Br", 1)))


# Source context: non-synthetic observations and characterization are not reactions.
bundle["source_context"] = [
    {"context_id": f"{DOC}:context:article_metadata", "kind": "article_metadata", "evidence_ids": [ev(1, "image"), ev(1, "text")],
     "title": TITLE, "authors": ["Sota Akiyama", "Syogo Nomura", "Koji Kubota", "Hajime Ito"],
     "citation": "J. Org. Chem. 2020, 85, 4172–4181", "doi": "10.1021/acs.joc.9b03353"},
    {"context_id": f"{DOC}:context:scheme1a_bioactive_examples", "kind": "nonreaction_illustration", "evidence_ids": [ev(1, "image")],
     "description": "Scheme 1a shows difluorovinyl-containing bioactive compounds and a carbonyl-bioisostere comparison; no synthetic arrow is reported."},
    {"context_id": f"{DOC}:context:general_methods", "kind": "general_methods", "evidence_ids": [ev(4, "text"), ev(4, "image")],
     "description": "General solvent handling, analytical instrumentation, chromatographic methods, and mass-spectrometry methods."},
    {"context_id": f"{DOC}:context:characterization_pages4_6", "kind": "characterization", "evidence_ids": [ev(4, "text"), ev(5, "text"), ev(6, "text")],
     "description": "NMR, HRMS, physical-state, and related characterization blocks for 1a–1k and 3a–3i. These blocks were not converted into synthesis experiments."},
    {"context_id": f"{DOC}:context:characterization_pages7_8", "kind": "characterization",
     "evidence_ids": [ev(7, "image"), ev(7, "text"), ev(8, "image"), ev(8, "text")],
     "description": "NMR, HRMS, physical-state, and related characterization blocks for 3i–3j, 4a–4i, and 6. These analytical peak lists were not converted into synthesis experiments."},
    {"context_id": f"{DOC}:context:associated_content_notice", "kind": "associated_content_notice",
     "evidence_ids": [ev(8, "image"), ev(8, "text")],
     "description": "The article advertises Supporting Information, but no SI file was supplied or explicitly associated in source_mapping.json; no SI facts were imported."},
    {"context_id": f"{DOC}:context:references_pages8_10", "kind": "references",
     "evidence_ids": [ev(8, "image"), ev(8, "text"), ev(9, "image"), ev(9, "text"), ev(10, "image"), ev(10, "text")],
     "description": "Reference-list entries and cited prior work only; not converted into current-paper synthesis experiments."},
]

bundle["molecules"] = list(molecules.values())
(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")

# Raw visual inventory kept separate from the formal submission.
raw_inventory = {
    "document_id": DOC,
    "page_range": [1, 6],
    "semantic_accuracy_verified": False,
    "pages": {
        "1": ["title/abstract", "graphical abstract", "Scheme 1a bioactive examples", "Scheme 1b generic current-work borylation", "Scheme 1c two generic utility arrows"],
        "2": ["Table 1: 12 optimization rows", "Table 2: 10 successful products 3a–3j and failed/complex-mixture 1k case", "Table footnotes and discussion"],
        "3": ["Table 3: nine cross-coupling rows", "Scheme 2 Diels-Alder experiment", "maleic-anhydride unreactive observation", "Figure 1 proposed mechanism"],
        "4": ["General methods", "Procedure A: three isolated operations toward 1a", "preparation references and characterization for 1b–1f"],
        "5": ["characterization/preparation references for 1g–1i", "three isolated operations toward 1j", "three isolated operations toward 1k", "start of general borylation procedure"],
        "6": ["completion of 3a general borylation procedure", "individual 3a–3i paragraphs", "NMR/HRMS characterization context"],
    },
    "limitations": ["Only pages 1–6 were visually reviewed in this session.", "Later page text was used only to preserve the explicit 3j time conflict; no later-page visual coverage is claimed.", "Supporting Information was not supplied or associated."],
}
(RUN / "raw_model/visual_transcription_pages_0001_0006.json").write_text(json.dumps(raw_inventory, ensure_ascii=False, indent=2) + "\n")

records_by_page = {i: [] for i in range(1, 7)}
for rxn in bundle["reactions"]:
    for eid in rxn["evidence_ids"]:
        m = re.search(r":page(\d+):", eid)
        if m and 1 <= int(m.group(1)) <= 6 and rxn["reaction_id"] not in records_by_page[int(m.group(1))]:
            records_by_page[int(m.group(1))].append(rxn["reaction_id"])
coverage = {
    "document_id": DOC, "semantic_accuracy_verified": False, "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title, abstract, graphical abstract", "finding": "Article metadata and a generic current-work transformation; not a concrete experiment.", "evidence_ids": [ev(1, "image"), ev(1, "text")]},
            {"region": "Scheme 1a", "finding": "Bioactive structures and carbonyl-bioisostere comparison; no synthetic experiment.", "evidence_ids": [ev(1, "image")]},
            {"region": "Scheme 1b–c", "finding": "Three generic arrows: boryl substitution, cross-coupling utility, and Diels–Alder utility.", "evidence_ids": [ev(1, "image")]},
        ], "record_ids": records_by_page[1], "unresolved": ["Generic R-group structures were not expanded into concrete SMILES."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Table 1 including scheme, rows 1–12, ligand drawings and footnotes", "finding": "Twelve independent optimization experiments, including two trace outcomes and a no-ligand row.", "evidence_ids": [ev(2, "image"), ev(2, "text")]},
            {"region": "Table 2 scheme, products 3a–3k and footnotes", "finding": "Ten successful substrate-scope experiments 3a–3j and one complex-mixture control 1k; Z/E ratios retained for 3g–3j.", "evidence_ids": [ev(2, "image"), ev(2, "text")]},
            {"region": "Results discussion", "finding": "Narrative repeats Table 1/2 outcomes; merged as additional evidence rather than duplicated.", "evidence_ids": [ev(2, "text")]},
        ], "record_ids": records_by_page[2], "unresolved": ["Many drawn structures lack independently verified SMILES.", "The Table 2 3j time conflicts with a later individual paragraph (7 h vs 6.5 h)."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Table 3 rows 4a–4i and footnotes", "finding": "Nine independent cross-coupling experiments; 4i yield is by GC, the others are reported under isolated-yield heading.", "evidence_ids": [ev(3, "image"), ev(3, "text")]},
            {"region": "Scheme 2", "finding": "Specific Diels–Alder reaction of 3b with dienophile 5 to product 6 (30%).", "evidence_ids": [ev(3, "image"), ev(3, "text")]},
            {"region": "discussion below Table 3", "finding": "Maleic anhydride is reported unreactive; retained as a no-reaction control without inventing a product.", "evidence_ids": [ev(3, "text")]},
            {"region": "Figure 1a–b", "finding": "Author-proposed catalytic cycle and Z/E isomerization pathway; recorded as proposed mechanism, not observed synthesis.", "evidence_ids": [ev(3, "image"), ev(3, "text")]},
        ], "record_ids": records_by_page[3], "unresolved": ["Scheme 2 does not state the 30% yield measurement basis on pages 1–6."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "General", "finding": "Analytical and handling methods; retained as source_context, not reactions.", "evidence_ids": [ev(4, "image"), ev(4, "text")]},
            {"region": "Procedure A for 1a", "finding": "Three isolated transformations: alcohol formation (67%), mesylation (61%), and Pd/organozinc allene formation (45%).", "evidence_ids": [ev(4, "image"), ev(4, "text")]},
            {"region": "1a–1f characterization/preparation paragraphs", "finding": "Characterization kept as context; 1b–1f preparation references to Procedure A recorded without inventing missing yields.", "evidence_ids": [ev(4, "image"), ev(4, "text")]},
        ], "record_ids": records_by_page[4], "unresolved": ["Unnamed intermediates and corresponding ketones remain identity-unresolved.", "Product-specific scales/yields for 1b–1f are unreported here."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "1g–1i paragraphs", "finding": "Procedure-A preparation references and characterization; no product-specific yields printed.", "evidence_ids": [ev(5, "image"), ev(5, "text")]},
            {"region": "Preparation of 1j", "finding": "Three concrete isolated operations: propargyl alcohol 97%, chloride 38%, allene 1j 9%.", "evidence_ids": [ev(5, "image"), ev(5, "text")]},
            {"region": "Preparation of 1k", "finding": "Three concrete isolated operations: gem-dibromo intermediate 76%, mesylate 60%, allene 1k 47%.", "evidence_ids": [ev(5, "image"), ev(5, "text")]},
            {"region": "start of General Borylation Procedure", "finding": "Charges and glovebox setup for 3a; merged into the Table 2 3a record.", "evidence_ids": [ev(5, "image"), ev(5, "text")]},
        ], "record_ids": records_by_page[5], "unresolved": ["Identity of the unnamed alkyne and several intermediates is not independently resolved."]},
        {"page": 6, "state": "attempted", "regions": [
            {"region": "completion of 3a procedure", "finding": "Reaction, workup, purification, 83% yield, product amount and trace-byproduct note merged into scope record 3a.", "evidence_ids": [ev(6, "image"), ev(6, "text")]},
            {"region": "individual 3b–3i paragraphs", "finding": "Product-specific substrate charges, times, yields, product masses/amounts, and Z/E results merged into existing Table 2 records.", "evidence_ids": [ev(6, "image"), ev(6, "text")]},
            {"region": "NMR/HRMS blocks", "finding": "Characterization retained in source_context and not counted as synthesis experiments.", "evidence_ids": [ev(6, "image"), ev(6, "text")]},
        ], "record_ids": records_by_page[6], "unresolved": ["3c product amount is printed/extracted as '395 mmol' next to 64.9 mg and was left unreadable rather than corrected.", "3i characterization continues to page 7; no page-7 visual coverage is claimed."]},
    ]
}
(RUN / "coverage/pages_0001_0006.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")

# Separate, exact visual inventory and semantic coverage for the requested pages 7--10.
raw_inventory_7_10 = {
    "document_id": DOC,
    "page_range": [7, 10],
    "semantic_accuracy_verified": False,
    "pages": {
        "7": [
            "continuation of 3i characterization",
            "individual 3j borylation paragraph and characterization",
            "full cross-coupling procedure for 4a",
            "individual 4b–4f paragraphs and characterization",
        ],
        "8": [
            "continuation of 4f characterization",
            "individual 4g–4i paragraphs and characterization",
            "full Diels–Alder procedure and product 6 characterization",
            "Associated Content, author information, acknowledgments, and start of references",
        ],
        "9": ["reference list entries 2–14 only"],
        "10": ["reference list entries 15–21 only"],
    },
    "limitations": [
        "Pages 7–10 were individually opened from the rendered original PDF images at original detail.",
        "Characterization peak lists and reference-list entries were retained as source context rather than synthesis experiments.",
        "Supporting Information was advertised but not supplied or explicitly associated.",
        "Machine extraction remains unreviewed; semantic_accuracy_verified is false.",
    ],
}
(RUN / "raw_model/visual_transcription_pages_0007_0010.json").write_text(
    json.dumps(raw_inventory_7_10, ensure_ascii=False, indent=2) + "\n")

prefix = DOC + ":"
coverage_7_10 = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 7,
            "state": "attempted",
            "regions": [
                {"region": "top-left continuation from page 6", "finding": "NMR/HRMS characterization for 3i; no new synthesis experiment.", "evidence_ids": [ev(7, "image"), ev(7, "text")]},
                {"region": "left column, 3j paragraph", "finding": "Individual borylation experiment 1j to 3j: 144.8 mg (0.50 mmol), 6.5 h, GPC, 45%, 90.4 mg (0.23 mmol), Z/E 79/21 by 19F NMR of the crude mixture; the 6.5 h value conflicts with Table 2's 7 h.", "evidence_ids": [ev(7, "image"), ev(7, "text")]},
                {"region": "left-to-right columns, Cross-Coupling Reaction Procedure for 4a", "finding": "Full setup, charges, nitrogen handling, 60 °C/25 h reaction, workup, chromatography, and 4a product amount/yield.", "evidence_ids": [ev(7, "image"), ev(7, "text")]},
                {"region": "remaining columns, individual 4b–4f paragraphs", "finding": "Product-specific substrate charges, times, product amounts, yields, and 4f E/Z ratio merged into the existing Table 3 records; spectral blocks are characterization context.", "evidence_ids": [ev(7, "image"), ev(7, "text")]},
            ],
            "record_ids": [
                prefix + "scope_table2_3i", prefix + "scope_table2_3j",
                prefix + "scope_table3_4a", prefix + "scope_table3_4b", prefix + "scope_table3_4c",
                prefix + "scope_table3_4d", prefix + "scope_table3_4e", prefix + "scope_table3_4f",
            ],
            "unresolved": [
                "3j reaction time remains a source conflict: Table 2 gives 7 h, while the page-7 paragraph gives 6.5 h.",
                "Complete source-verified structures/atom mappings for these labeled compounds remain unavailable; concrete SMILES were not invented.",
            ],
        },
        {
            "page": 8,
            "state": "attempted",
            "regions": [
                {"region": "top-left continuation", "finding": "Continuation of 4f characterization; no separate experiment.", "evidence_ids": [ev(8, "image"), ev(8, "text")]},
                {"region": "left and middle columns, 4g–4i paragraphs", "finding": "Individual cross-coupling times, substrate/product amounts, and yields merged into existing Table 3 records. The 4i paragraph's product mass/amount is retained alongside Table 3's GC footnote as a basis conflict.", "evidence_ids": [ev(8, "image"), ev(8, "text")]},
                {"region": "middle/right columns, Diels–Alder procedure", "finding": "Full experiment for 3b plus dienophile 5: exact charges, toluene volume, nitrogen handling, room temperature for 21 h, workup, chromatography, and isolated 6 (30%, 14.7 mg, 0.030 mmol).", "evidence_ids": [ev(8, "image"), ev(8, "text")]},
                {"region": "right column lower section", "finding": "Product 6 characterization followed by Associated Content, author information, acknowledgments, and reference-list material; retained as context, not new experiments.", "evidence_ids": [ev(8, "image"), ev(8, "text")]},
            ],
            "record_ids": [
                prefix + "scope_table3_4f", prefix + "scope_table3_4g", prefix + "scope_table3_4h",
                prefix + "scope_table3_4i", prefix + "scheme2_3b_5_to_6",
            ],
            "unresolved": [
                "4i yield basis remains conflicting/ambiguous between the Table 3 GC footnote and the matching product mass/amount in the individual paragraph.",
                "The advertised Supporting Information was not supplied or explicitly associated, so it was not used.",
            ],
        },
        {
            "page": 9,
            "state": "attempted",
            "regions": [
                {"region": "full page, references 2–14", "finding": "Cited prior literature only; no current-paper experimental arrow, table row, or procedure was converted into a record.", "evidence_ids": [ev(9, "image"), ev(9, "text")]},
            ],
            "record_ids": [],
            "unresolved": ["Reference-list descriptions were not treated as evidence of experiments performed in this paper."],
        },
        {
            "page": 10,
            "state": "attempted",
            "regions": [
                {"region": "full page, references 15–21", "finding": "Cited prior literature only; no current-paper experiment is present on this page.", "evidence_ids": [ev(10, "image"), ev(10, "text")]},
            ],
            "record_ids": [],
            "unresolved": ["Reference-list descriptions were not imported as this article's reaction records."],
        },
    ],
}
(RUN / "coverage/pages_0007_0010.json").write_text(
    json.dumps(coverage_7_10, ensure_ascii=False, indent=2) + "\n")

workflow = f"""# Workflow log — allene-stereo-v31-200-20261017

## Source and bounded scope

- Document: `{DOC}`; DOI `10.1021/acs.joc.9b03353`.
- Title: *{TITLE}*.
- Source mapping reports `presentation_conversion=false`; this is the journal PDF, not a PPTX-derived presentation.
- SI association remains unestablished. The article advertises Supporting Information, but no SI file was supplied or associated, so no SI facts were imported.
- Original rendered images for PDF pages 1–10 were individually opened at original detail across the two bounded passes. The present requested pass is restricted to exactly pages 7–10 and is documented separately.
- The pre-existing workflow described a different paper (`doc_5909d9a1e9253598`) and was treated as stale contamination. It was not used as source evidence and this file replaces it.

## Page and record decisions

- Page 1: Scheme 1b and the two Scheme 1c arrows are generic current-work transformations. Scheme 1a is a nonreaction bioactive/bioisostere illustration.
- Page 2: Table 1 contributes 12 independent optimization records. Table 2 contributes 10 successful scope experiments and one complex-mixture control. Table and individual-procedure evidence are merged by experiment, not duplicated.
- Page 3: Table 3 contributes nine cross-coupling experiments; Scheme 2 contributes one Diels–Alder experiment; maleic anhydride is a real no-reaction observation; Figure 1 is one proposed-mechanism record.
- Page 4: Procedure A for 1a is split into three isolated transformations. 1b–1f are separate source-specific preparations that explicitly reference Procedure A; missing yields/scales remain unreported. Spectral blocks are source context.
- Page 5: 1g–1i likewise reference Procedure A. The preparations of 1j and 1k each contain three isolated transformations and were split accordingly. The beginning of the 3a procedure is merged with the Table 2 record.
- Page 6: the 3a–3i individual paragraphs enrich their Table 2 records. NMR/HRMS peak lists are source context, not synthesis experiments.
- Page 7: the continuation of 3i is characterization only. The 3j paragraph enriches the existing Table 2 record and visually confirms the 7 h versus 6.5 h conflict. The full 4a procedure and individual 4b–4f paragraphs enrich the existing Table 3 records.
- Page 8: the 4f continuation is characterization only. Individual 4g–4i paragraphs enrich Table 3 records; 4i retains the table-GC versus matching product-amount ambiguity. The full Diels–Alder procedure resolves the prior missing yield basis for 6 as isolated after chromatography.
- Pages 9–10: reference lists only. Cited prior work was retained as source context and was not converted into this paper's experiments.
- Exact regions and record IDs are in `runs/current/coverage/pages_0001_0006.json` and `runs/current/coverage/pages_0007_0010.json`; the latter lists exactly pages 7–10 and every page has nonempty regions.

## Stereochemistry and structures

- Every molecule has explicit chirality, geometry, and representation status facts. CF3-substituted allenes 1a–1k retain an unspecified allene axis; absence of an assignment is not converted to achiral or racemic.
- Reported Z/E pairs for 3g–3j and E/Z for 4f are registered as separate geometric-isomer molecule IDs. Combined yields target both product isomers, and selectivity outcomes preserve the printed order and ratio.
- No ee/er/dr or absolute R/S/Ra/Sa assignment is reported on pages 1–10. No configuration was invented.
- Primary source drawings were inspected, but complete concrete SMILES were not independently verified. They remain empty and `representation_status=not_encoded/unresolved`; schema-validity is not treated as structure correctness.
- Graph-only `standard stereo` diagnostics are saved separately under `raw_model/`; `source_verified=false` remains explicit and the diagnostics were not promoted into source structures.
- In this pass, `standard stereo` was rerun for the printed 4-bromoanisole name candidate and dienophile 5 name candidate. Both outputs remain graph diagnostics only; the allene-axis limitation and absent source-verified product graphs prevent blanket SMILES promotion.

## Quantities, conflicts, and bond classification

- All quantities use explicit unit objects; qualitative `trace`, `complex mixture`, `r.t.`, `reflux`, and `overnight` remain qualitative.
- 3a PPh3 charge conflict is retained: Table 2 says 0.025 mmol, while the individual paragraph says 0.026 mmol.
- 1b scale differs: Table 2 gives 0.50 mmol; the paragraph gives 100.8 mg (0.49 mmol). Both are retained as explicit conflict alternatives; neither is silently selected.
- 3j time conflict is retained using Table 2 (7 h) and the visually inspected page-7 individual paragraph (6.5 h).
- 3c product amount appears as `395 mmol` beside 64.9 mg and a 0.25 mmol reaction scale. It is marked unreadable and was not silently repaired.
- Scheme 2's 30% yield basis was absent from pages 1–6 but is resolved by the page-8 procedure as isolated after silica chromatography; the reported 14.7 mg and 0.030 mmol are separately targeted to product 6.
- Table 3 gives 0.10 mmol of substrate for every cross-coupling row, whereas the individual paragraphs give 0.102 mmol for 4d, 0.0997 mmol for 4h, and 0.102 mmol for 4i. Both values and evidence are retained as conflicts.
- Table 3 footnote c states that 4i's 55% was determined by GC, while the individual paragraph reports 24.9 mg (0.056 mmol), numerically matching 55% of its stated charge. The yield basis is kept conflicting/ambiguous rather than normalized to either GC or isolated yield.
- Bond classifications are partial structure-derived net inventories. Defluoroborylation records include C–B formation, C–C bond-order increase, C–F cleavage, and a C–C bond-order decrease; cross-couplings include C–C formation with C–B/C–Br cleavage; cycloaddition records include two C–N formations and diene bond-order changes. No complete atom mapping or three-model vote was fabricated.

## Audit disposition

- Strict validation and standard validation pass. The audit reports 51 warned records, all with `structure_gap`; it reports no `reaction_identity`, `outcome_target`, or `measurement_basis` warnings after the page-7/8 enrichment.
- The `structure_gap` warnings remain unresolved because source-verified concrete SMILES/atom mappings are unavailable for the drawn/labeled reactant-product pairs. Raw graph candidates are kept in `raw_model/stereo_diagnostics.json` with `source_verified=false`; inventing structures merely to silence the audit was rejected.
- All 55 reaction records remain `needs_review` and `partial`. Validation establishes schema consistency only, not chemical correctness.

## Artifacts and validation boundary

- Fresh init templates: `runs/current/raw_model/standard_init_template.json` and `runs/current/raw_model/standard_init_template_pages_0007_0010.json`.
- Raw visual inventories: `runs/current/raw_model/visual_transcription_pages_0001_0006.json` and `runs/current/raw_model/visual_transcription_pages_0007_0010.json`; literal mention hints: `runs/current/raw_model/source_mentions.json`; reproducible builder: `runs/current/raw_model/build_submission.py`.
- Formal submission: `runs/current/submitted.json`; coverage: `runs/current/coverage/pages_0001_0006.json` and `runs/current/coverage/pages_0007_0010.json`.
- All machine reactions remain `needs_review` and `partial`. `semantic_accuracy_verified=false` and `chemical_accuracy_verified=false`; validation establishes schema consistency only.
"""
(ROOT / "runs/workflow.md").write_text(workflow)
print(json.dumps({
    "molecules": len(bundle["molecules"]),
    "reactions": len(bundle["reactions"]),
    "requested_pages": [7, 8, 9, 10],
    "enriched_existing_records": 11,
    "new_reaction_records": 0,
}, ensure_ascii=False))
