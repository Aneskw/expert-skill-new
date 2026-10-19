#!/usr/bin/env python3
"""Build the source-derived v3.1 machine submission for doc_ad5aa762c966de92.

This file is a reproducible extraction artifact, not a source-review certificate.
Dense drawn structures are deliberately left unresolved rather than guessed.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DOC = "doc_ad5aa762c966de92"
E = {p: f"{DOC}:page{p}:image" for p in range(1, 6)}
T = {p: f"{DOC}:page{p}:text" for p in range(1, 6)}


def f(raw=None, norm=None, ev=(), method="original_extraction", derivation=None,
      status=None):
    return st.fact(raw, norm, ev, method=method, derivation=derivation, status=status)


def missing(status="unreported"):
    x = st.fact()
    x["data_status"] = status
    x["acquisition_method"] = "not_applicable" if status == "not_applicable" else "original_extraction"
    return x


def q(number, unit, raw=None, qualifier=None):
    value = {"number": number, "unit": unit}
    if qualifier:
        value["qualifier"] = qualifier
    if raw:
        value["raw_text"] = raw
    return value


bundle = json.loads((RUN / "submitted.json").read_text())
bundle["molecules"] = []
bundle["reactions"] = []
bundle["source_context"] = []
molecules = {}


def add_molecule(mid, label, ev, *, name=None, smiles=None, chirality="unresolved",
                 geometry="unresolved", representation=None, stereo_raw=None,
                 mixture=None, axial=None, identity_note=None):
    m = st.new_molecule(mid, DOC, list(ev))
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    m["names"]["original"] = f(name or label, name or label, ev)
    m["names"]["english"] = f(name or label, name or label, ev)
    m["names"]["abbreviations"] = missing()
    m["names"]["aliases"] = missing()
    m["cas_numbers"] = missing()
    m["external_properties"] = []
    if smiles:
        diag = st.stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = f({"source": "candidate transcribed or named structure", "smiles": smiles},
                         smiles, ev, method="structure_derivation",
                         derivation="Candidate structure derived from the legible drawing or unambiguous common name; RDKit parsing does not certify source identity.")
        representation = representation or "encoded"
        m["stereochemistry"]["isomeric_smiles"] = f(None, diag["canonical_isomeric_smiles"], ev,
            method="structure_derivation", derivation="RDKit canonical isomeric SMILES for the candidate graph.")
        m["stereochemistry"]["identity_key"] = f(None, diag["stereo_graph_key"], ev,
            method="structure_derivation", derivation="RDKit stereo graph key for the candidate graph; not a sample identity.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing()
        representation = representation or "not_encoded"
        m["stereochemistry"]["isomeric_smiles"] = missing()
        m["stereochemistry"]["identity_key"] = missing()
    m["stereochemistry"]["chirality_status"] = f(stereo_raw or chirality, chirality, ev,
        method="inference", derivation="Status assigned from the source drawing/text while preserving any unassigned or racemic state.")
    m["stereochemistry"]["geometry_status"] = f(geometry, geometry, ev,
        method="inference", derivation="E/Z or other alkene geometry status inspected in the cited source region.")
    m["stereochemistry"]["representation_status"] = f(representation, representation, ev,
        method="inference", derivation="Whether the candidate machine representation retains the source-reported stereochemical state.")
    m["stereochemistry"]["atom_centers"] = f(stereo_raw, stereo_raw, ev) if stereo_raw else missing()
    m["stereochemistry"]["double_bonds"] = missing()
    m["stereochemistry"]["axial"] = f(axial, axial, ev) if axial else missing()
    m["stereochemistry"]["mixture"] = f(mixture, mixture, ev) if mixture else missing()
    m["stereochemistry"]["assignment"] = f(stereo_raw, stereo_raw, ev) if stereo_raw else missing()
    m["identity_note"] = f(identity_note, identity_note, ev) if identity_note else missing()
    molecules[mid] = m
    bundle["molecules"].append(m)
    return mid


# Solvents and simple, unambiguously named nucleophiles.
add_molecule("mol_chcl3", "CHCl3", (E[2], E[3]), name="chloroform", smiles="ClC(Cl)Cl",
             chirality="achiral", geometry="not_applicable")
add_molecule("mol_dcm", "DCM", (E[2],), name="dichloromethane", smiles="ClCCl",
             chirality="achiral", geometry="not_applicable")
add_molecule("mol_dce", "DCE", (E[2],), name="1,2-dichloroethane", smiles="ClCCCl",
             chirality="achiral", geometry="not_applicable")
for mid, label, name, smi in [
    ("mol_meoh", "MeOH", "methanol", "CO"),
    ("mol_etoh", "EtOH", "ethanol", "CCO"),
    ("mol_nproh", "nPrOH", "1-propanol", "CCCO"),
    ("mol_iproh", "iPrOH", "2-propanol", "CC(O)C"),
    ("mol_tbuoh", "tBuOH", "tert-butanol", "CC(C)(C)O"),
    ("mol_propargyl_alcohol", "propargyl alcohol", "prop-2-yn-1-ol", "C#CCO"),
    ("mol_4methoxyphenol", "4-methoxyphenol", "4-methoxyphenol", "COc1ccc(O)cc1"),
    ("mol_135trimethoxybenzene", "1,3,5-trimethoxybenzene", "1,3,5-trimethoxybenzene", "COc1cc(OC)cc(OC)c1"),
]:
    add_molecule(mid, label, (E[3], T[3]), name=name, smiles=smi,
                 chirality="achiral", geometry="not_applicable")

# Substrates and allenene intermediates. Exact dense graphs are not guessed.
for label in ["7a", "7b", "7c", "7d"]:
    add_molecule(f"mol_{label}", label, (E[2], E[3]),
                 chirality="achiral", geometry="not_applicable", representation="not_encoded",
                 identity_note="Fully drawn naphthol-derived propargyl ether; exact machine SMILES not asserted in this run.")
add_molecule("mol_7e", "7e", (E[3], T[3]), chirality="unspecified", geometry="not_applicable",
             stereo_raw="a stereogenic center is present in the propargylic unit; absolute configuration is not reported on supplied pages",
             identity_note="Structure appears only by textual reference to Scheme S2; no associated SI was supplied.")
for label in ["5a", "5b", "5c", "5d"]:
    add_molecule(f"mol_{label}", label, (E[2], E[3]), chirality="unspecified",
                 geometry="not_applicable", representation="not_encoded",
                 axial="potential allene axial stereogenicity; R_a/S_a or P/M not assigned in the supplied main-text pages",
                 identity_note="Drawn 1,5-allenene; axial configuration is not assigned and is not machine encoded.")

# Fused products. The schemes explicitly state racemic mixtures and show one enantiomer.
racemic_text = "products were isolated as racemic mixtures; only one enantiomer is shown"
relative_text = "source drawing specifies relative wedge/dash relationships for the fused product; absolute configuration is not assigned"
product_labels = ["6a", "6b", "6c", "6d", "6e", "6f", "6g", "6h", "6i", "6i_prime",
                  "6j", "6j_prime", "6k", "6l", "6m", "6n", "6o", "6p", "6q", "6r",
                  "6s", "6t", "6u", "6v", "6w"]
for label in product_labels:
    shown = label.replace("_prime", "′")
    add_molecule(f"mol_{label}", shown, (E[3],), chirality="racemic", geometry="specified",
                 representation="not_encoded", stereo_raw=relative_text, mixture=racemic_text,
                 identity_note="Dense 6-6-5 fused product is visually identified by label; exact connectivity/stereochemical SMILES was not guessed.")
for label in ["8a", "8b"]:
    add_molecule(f"mol_{label}", label, (E[3],), chirality="racemic", geometry="specified",
                 representation="not_encoded", stereo_raw=relative_text, mixture=racemic_text,
                 identity_note="One member of the source-reported pair of diastereomeric ether products; exact graph not encoded.")
for label in ["6x", "6x_prime"]:
    add_molecule(f"mol_{label}", label.replace("_prime", "′"), (E[3], T[3]),
                 chirality="unspecified", geometry="unresolved", representation="not_encoded",
                 stereo_raw="one of two source-reported diastereomers; configurations are not shown on supplied pages",
                 identity_note="Referenced to unavailable Scheme S2; exact structure and configuration remain unresolved.")

# Nucleophiles that are drawn but not safely nameable from the supplied main-text resolution.
for key, desc in [
    ("cresol", "cresol nucleophile used for 6i/6i′; positional identity not normalized beyond the source drawing"),
    ("chlorophenol", "chlorophenol nucleophile used for 6k; positional identity retained only in the drawing"),
    ("iodophenol", "iodophenol nucleophile used for 6l; positional identity retained only in the drawing"),
    ("phenol_6m", "alkenyl-substituted phenol drawn as the nucleophile leading to 6m"),
]:
    add_molecule(f"mol_nu_{key}", desc, (E[3], T[3]), chirality="unresolved",
                 geometry="unresolved", representation="not_encoded", identity_note=desc)
for label in ["6o", "6p", "6q", "6r", "6s", "6t", "6v", "6w"]:
    add_molecule(f"mol_nu_{label}", f"carbon nucleophile leading to {label} (Scheme 4 drawing)",
                 (E[3],), chirality="unresolved", geometry="unresolved", representation="not_encoded",
                 identity_note="Input identity is retained by product linkage and source drawing; exact standalone SMILES was not guessed.")

# Optimization catalyst systems, preserved exactly as table-cell identities.
cat_names = [
    "[JohnPhosAu(NCMe)]SbF6", "XPhosAu(NCMe)SbF6", "BrettPhosAuNTf2", "tBuXPhosAuNTf2",
    "CyJohnPhosAuCl/AgSbF6", "Ph3PAuCl/AgSbF6", "Me3PAuCl/AgSbF6", "IPrAuCl/AgNTf2",
    "IPrAuCl/NaBArF24", "IPrAuCl/AgSbF6", "Dichloro(2-picolinato)gold(III)",
    "JohnPhosAuCl/AgSbF6", "JohnPhosAuCl/AgNTf2", "JohnPhosAuCl/NaBArF24",
]
for i, name in enumerate(cat_names, 1):
    add_molecule(f"mol_cat_{i}", name, (E[2],), name=name, chirality="unresolved",
                 geometry="not_applicable", representation="not_encoded",
                 identity_note="Catalyst or catalyst/activator system transcribed literally from Table 1; multicomponent slash entries are not collapsed into a guessed single molecule.")
for key, name in [
    ("sbf6", "SbF6 counterion/component indicated by the Table 1 catalyst cell"),
    ("ntf2", "NTf2 counterion/component indicated by the Table 1 catalyst cell"),
    ("barf24", "BArF24 counterion/component indicated by the Table 1 catalyst cell"),
    ("unknown", "counterion reported by the Table 1 heading; identity unresolved for this catalyst cell"),
]:
    add_molecule(f"mol_counterion_{key}", name, (E[2],), name=name,
                 chirality="unresolved", geometry="not_applicable", representation="not_encoded",
                 identity_note="Registered separately to preserve the explicit 7 mol% table-heading quantity; whether it is a separately dosed activator or intrinsic counterion remains under review.")


def participant(pid, mid, role, ev, *, equiv=None, molpct=None, stage="reaction", reference="mol_7a",
                charge_method="original_extraction"):
    p = st.new_participant(pid, mid, role, ev, stage=stage)
    if role == "product":
        p["product_role"] = f("product", "product", ev)
    if equiv is not None:
        p["charges"]["equivalents"] = f(f"{equiv} equiv", q(equiv, "equiv"), ev, method=charge_method)
        p["equivalence_reference"] = f("substrate", reference, ev, method=charge_method)
    if molpct is not None:
        p["charges"]["mol_percent"] = f(f"{molpct} mol%", q(molpct, "mol%"), ev, method=charge_method)
        p["equivalence_reference"] = f("substrate", reference, ev, method=charge_method)
    return p


def classification(ev, codes=(), counts=None, limitation=None):
    counts = counts or {}
    changes = []
    labels = []
    n = 0
    atom_pairs = {"C-C": ("C", "C", 1), "C-O": ("C", "O", 1), "C-H": ("C", "H", 1),
                  "DB-C-C": ("C", "C", 2)}
    for code in codes:
        support = []
        for _ in range(counts.get(code, 1)):
            n += 1
            a, b, order = atom_pairs[code]
            cid = f"b{n}"
            change_type = "hydrogen_gain" if code == "C-H" else ("order_increased" if code == "DB-C-C" else "formed")
            before = 1 if code == "DB-C-C" else 0
            changes.append(f(None, {"change_id": cid, "atoms": [{"atom_map": None, "element": a},
                                      {"atom_map": None, "element": b}], "before_order": before,
                                      "after_order": order, "change_type": change_type}, ev,
                             method="structure_derivation",
                             derivation="Partial bond-change assertion from the source reaction drawing; atom mapping is unresolved."))
            support.append(cid)
        lab = st.bond_label(*atom_pairs[code])
        lab["supporting_change_ids"] = support
        labels.append(lab)
    return {"standard_id": st.STANDARD,
            "mapped_reaction_smiles": missing(),
            "bond_changes": changes,
            "labels": f(None, labels, ev, method="structure_derivation",
                        derivation="Taxonomy labels derived only from the explicitly asserted partial bond changes."),
            "coverage": f(None, "partial" if codes else "unresolved", ev, method="inference",
                          derivation="Dense structures lack reviewed atom mapping; only visually explicit new connections are represented."),
            "limitations": [limitation or "Partial visual classification only; no complete atom mapping, atom balance, or coproduct accounting."],
            "model_votes": []}


def make_reaction(rid, ev, *, categories, location, rname, transform, parts, step_temp=None,
                  step_time=None, operation=None, solvent=None, concentration=None, outcomes=(),
                  class_codes=(), class_counts=None, notes=(), unresolved=(), record_type="specific_experiment"):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, list(ev))
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, ev)
    r["location"] = f(location, location, ev)
    r["reaction_name"] = f(rname, rname, ev)
    r["transformation_type"] = f(None, transform, ev, method="inference",
                                  derivation="Transformation category inferred from the source drawing and prose.")
    r["participants"] = parts
    if solvent:
        r["solvent_conditions"]["composition"] = f(solvent[1], [solvent[0]], ev)
    if concentration is not None:
        r["solvent_conditions"]["concentration"] = f(f"{concentration} M", q(concentration, "M"), ev)
        r["solvent_conditions"]["concentration_basis"] = f("7a", "mol_7a", ev)
    if operation or step_temp is not None or step_time is not None:
        s = st.new_step("reaction", ev)
        s["participant_ids"] = [p["participant_id"] for p in parts]
        if operation:
            s["operation"] = f(operation, operation, ev)
        if step_temp is not None:
            s["temperature"] = f(step_temp, step_temp if isinstance(step_temp, str) else q(step_temp, "°C"), ev)
        if step_time is not None:
            s["time"] = f(f"{step_time} h", q(step_time, "h"), ev)
        r["steps"] = [s]
    r["outcomes"] = list(outcomes)
    r["classification"] = classification(ev, class_codes, class_counts)
    r["notes"] = [f(x, x, ev) for x in notes]
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": list(unresolved) or [
                       "No associated SI was supplied; experimental quantities, workup, purification, and analytical yield basis may be absent."]}
    r["completeness"] = "partial"
    bundle["reactions"].append(r)
    return r


def outcome(oid, kind, targets, value, ev, *, raw=None, basis=None, role="main product", scope=None):
    o = st.new_outcome(oid, kind, targets, ev)
    o["product_role"] = f(role, role, ev)
    o["value"] = f(raw if raw is not None else value, value, ev)
    if kind == "yield":
        o["yield_type"] = f("reported yield", "reported", ev)
    if basis:
        o["measurement_basis"] = f(basis, basis, ev)
    if scope:
        o["scope"] = f(scope, scope, ev)
    return o


# Table 1: one optimization record per row. Column assignment was checked visually.
table_rows = [
    (1, 1, None, 87, "CHCl3", 2, 5), (2, 2, None, 78, "CHCl3", 2, 5),
    (3, 3, None, 80, "CHCl3", 2, 5), (4, 4, None, 79, "CHCl3", 2, 5),
    (5, 5, None, 76, "CHCl3", 2, 5), (6, 6, None, 65, "CHCl3", 2, 5),
    (7, 7, 57, None, "CHCl3", 2, 5), (8, 8, None, 55, "CHCl3", 2, 5),
    (9, 9, 36, 31, "CHCl3", 2, 5), (10, 10, None, 61, "CHCl3", 2, 5),
    (11, 11, 88, None, "CHCl3", 2, 5), (12, 12, None, 75, "CHCl3", 2, 5),
    (13, 13, None, 64, "CHCl3", 2, 5), (14, 14, 46, 36, "CHCl3", 2, 5),
    (15, 1, None, 82, "DCM", 2, 5), (16, 1, None, 53, "DCE", 2, 5),
    (17, 1, None, 56, "CHCl3", 1, 5), (18, 1, None, 84, "CHCl3", 3, 5),
    (19, 1, None, 83, "CHCl3", 2, 2.5), (20, 1, 60, 22, "CHCl3", 2, 1),
]
solvent_mids = {"CHCl3": "mol_chcl3", "DCM": "mol_dcm", "DCE": "mol_dce"}
counterion_by_cat = {1: "sbf6", 2: "sbf6", 3: "ntf2", 4: "ntf2", 5: "sbf6", 6: "sbf6",
                     7: "sbf6", 8: "ntf2", 9: "barf24", 10: "sbf6", 11: "unknown",
                     12: "sbf6", 13: "ntf2", 14: "barf24"}
for row, cat, y5, y6, solv, meoh_eq, cat_pct in table_rows:
    ev = (E[2], T[2])
    parts = [participant("substrate", "mol_7a", "reactant", ev),
             participant("methanol", "mol_meoh", "reagent", ev, equiv=meoh_eq),
             participant("catalyst", f"mol_cat_{cat}", "catalyst", ev, molpct=cat_pct),
             participant("counterion_component", f"mol_counterion_{counterion_by_cat[cat]}", "unknown", ev, molpct=7),
             participant("solvent", solvent_mids[solv], "solvent", ev)]
    outs = []
    codes = ["C-C", "DB-C-C"]
    if y5 is not None:
        parts.append(participant("product_5a", "mol_5a", "product", ev))
        outs.append(outcome("yield_5a", "yield", ["mol_5a"], q(y5, "%"), ev, raw=f"5a {y5}%", role="allenene product/intermediate"))
    if y6 is not None:
        parts.append(participant("product_6a", "mol_6a", "product", ev))
        outs.append(outcome("yield_6a", "yield", ["mol_6a"], q(y6, "%"), ev, raw=f"6a {y6}%"))
        codes += ["C-O", "C-H"]
    make_reaction(f"table1_entry_{row:02d}", ev, categories=["optimization"],
        location={"page": 2, "table": "Table 1", "entry": row},
        rname="Optimization of the gold-catalyzed formation of 5a and 6a",
        transform="in situ propargyl-ether rearrangement and/or 1,5-allenene cyclization/trapping",
        parts=parts, step_temp="room temperature", step_time=2,
        operation="Treat 7a with the listed catalyst system and methanol under the table conditions.",
        solvent=(solvent_mids[solv], solv), concentration=0.1, outcomes=outs,
        class_codes=codes, notes=("Table heading states 5 mol% catalyst and 7 mol% counterion; the catalyst-cell salt/system is retained literally. Applicability of a separately dosed 7 mol% activator to every preformed salt is unresolved.",),
        unresolved=("Yield measurement basis is not reported in Table 1.",
                    "The table heading reports 7 mol% counterion and that charge is structured; whether the cell-specific component is separately dosed or intrinsic to a preformed salt cannot be resolved safely.",
                    "No associated SI was supplied; reaction scale, volumes, workup and purification are unreported."))

# Scheme 2: four independently reported Au(III)-catalyzed rearrangements.
for label, y in [("a", 98), ("b", 51), ("c", 98), ("d", 75)]:
    ev = (E[2], T[2])
    parts = [participant("substrate", f"mol_7{label}", "reactant", ev),
             participant("catalyst", "mol_cat_11", "catalyst", ev),
             participant("solvent", "mol_chcl3", "solvent", ev),
             participant("product", f"mol_5{label}", "product", ev)]
    basis = "isolated" if label == "a" else None
    make_reaction(f"scheme2_7{label}_to_5{label}", ev, categories=["substrate_scope"],
        location={"page": 2, "scheme": "Scheme 2", "arrow": f"7{label} to 5{label}"},
        rname="Au(III)-catalyzed rearrangement to a 1,5-allenene",
        transform="propargyl aryl ether rearrangement to 1,5-allenene",
        parts=parts, step_temp="room temperature", step_time=2,
        operation="Au(III)-catalyzed rearrangement in CHCl3.", solvent=("mol_chcl3", "CHCl3"),
        outcomes=[outcome("yield", "yield", [f"mol_5{label}"], q(y, "%"), ev,
                          raw=f"5{label} ({y}%)", basis=basis)],
        class_codes=["C-C", "DB-C-C"],
        unresolved=(("Catalyst loading, reaction scale, workup and purification are not reported on the supplied page.",) +
                    (() if basis else ("Yield measurement basis is not explicitly reported for this Scheme 2 entry.",))))

# Scheme 3 oxygen-nucleophile scope. 6i/6i′ and 6j/6j′ are paired outcomes of one experiment each.
oxygen_rows = [
    ("6a", "mol_7a", "mol_meoh", [("mol_6a", 86)]),
    ("6b", "mol_7b", "mol_meoh", [("mol_6b", 35)]),
    ("6c", "mol_7c", "mol_meoh", [("mol_6c", 76)]),
    ("6d", "mol_7d", "mol_meoh", [("mol_6d", 34)]),
    ("6e", "mol_7a", "mol_etoh", [("mol_6e", 67)]),
    ("6f", "mol_7a", "mol_nproh", [("mol_6f", 64)]),
    ("6g", "mol_7a", "mol_iproh", [("mol_6g", 84)]),
    ("6h", "mol_7a", "mol_tbuoh", [("mol_6h", 72)]),
    ("6i", "mol_7a", "mol_nu_cresol", [("mol_6i", 46), ("mol_6i_prime", 38)]),
    ("6j", "mol_7a", "mol_4methoxyphenol", [("mol_6j", 17), ("mol_6j_prime", 42)]),
    ("6k", "mol_7a", "mol_nu_chlorophenol", [("mol_6k", 81)]),
    ("6l", "mol_7a", "mol_nu_iodophenol", [("mol_6l", 70)]),
    ("6m", "mol_7a", "mol_nu_phenol_6m", [("mol_6m", 12)]),
    ("6n", "mol_7a", "mol_propargyl_alcohol", [("mol_6n", 68)]),
]
for tag, substrate, nu, products in oxygen_rows:
    ev = (E[3], T[3], E[2], T[2])
    parts = [participant("substrate", substrate, "reactant", ev),
             participant("nucleophile", nu, "reagent", ev, equiv=2, reference=substrate,
                         charge_method="inherited_general_condition"),
             participant("catalyst", "mol_cat_1", "catalyst", ev, molpct=5, reference=substrate,
                         charge_method="inherited_general_condition"),
             participant("solvent", "mol_chcl3", "solvent", ev)]
    outs = []
    for j, (mid, y) in enumerate(products, 1):
        parts.append(participant(f"product_{j}", mid, "product", ev))
        outs.append(outcome(f"yield_{j}", "yield", [mid], q(y, "%"), ev,
                            raw=f"{mid.removeprefix('mol_').replace('_prime', '′')} ({y}%)",
                            role="main product" if j == 1 else "competing product"))
    main_targets = [p[0] for p in products]
    outs.append(outcome("diastereoselectivity", "selectivity", main_targets,
                        "one diastereomer isolated as an enantiomeric pair; numerical dr not reported", ev,
                        scope={"metric": "diastereoselectivity/racemic state", "order": "not numerically reported"}))
    make_reaction(f"scheme3_{tag}", ev, categories=["substrate_scope"],
        location={"page": 3, "scheme": "Scheme 3", "product_label": tag},
        rname="Gold(I)-catalyzed 1,5-allenene cyclization using an oxygen nucleophile",
        transform="in situ allenene formation, 5-exo cyclization and oxygen-nucleophile trapping",
        parts=parts, step_temp="room temperature", step_time=3,
        operation="Generate allenene 5 in situ from alkyne 7 and cyclize/trap with the depicted oxygen nucleophile.",
        solvent=("mol_chcl3", "CHCl3"), outcomes=outs,
        class_codes=["C-C", "DB-C-C", "C-O", "C-H"],
        notes=(racemic_text,),
        unresolved=("All products are reported as racemic mixtures; absolute configurations and numerical dr are not reported.",
                    "Yield measurement basis is not explicitly stated in the Scheme 3 footnote.",
                    "Exact dense product SMILES and full atom mapping remain unresolved.",
                    "The 2 equiv nucleophile and 5 mol% catalyst loadings are inherited from the optimized conditions described on page 2; the main-paper Scheme 3 footnote does not restate them."))

# No externally added nucleophile: combined yield of two explicitly distinct diastereomers.
ev = (E[3], T[3])
parts = [participant("substrate", "mol_7a", "reactant", ev),
         participant("catalyst", "mol_cat_1", "catalyst", ev, molpct=5,
                     charge_method="inherited_general_condition"),
         participant("solvent", "mol_chcl3", "solvent", ev),
         participant("product_8a", "mol_8a", "product", ev),
         participant("product_8b", "mol_8b", "product", ev)]
make_reaction("scheme3_no_added_nucleophile_8a_8b", ev, categories=["control", "mechanistic_study"],
    location={"page": 3, "scheme": "Scheme 3", "products": "8a + 8b"},
    rname="Gold(I)-catalyzed cyclization without added nucleophile",
    transform="cyclization followed by source-proposed water/alcohol-mediated ether formation",
    parts=parts, step_temp="room temperature", step_time=3,
    operation="Run the Au(I)-catalyzed reaction without adding an external nucleophile.",
    solvent=("mol_chcl3", "CHCl3"),
    outcomes=[outcome("combined_yield", "yield", ["mol_8a", "mol_8b"], q(77, "%"), ev,
                      raw="8a + 8b (77%)", role="mixture of two diastereomeric ether products")],
    class_codes=["C-C", "DB-C-C", "C-O", "C-H"],
    notes=("The authors propose that traces of water form an alcohol which reacts with another benzylic cation; this mechanistic explanation is not a separately measured reagent dose.",),
    unresolved=("Yield measurement basis and the 8a:8b ratio are not reported.",
                "Water is described only as probable trace contamination and is not assigned a measured participant quantity.",
                "Exact structures and atom mapping remain unresolved."))

# Scheme 4 carbon-nucleophile scope.
carbon_rows = [("6o", 68), ("6p", 31), ("6q", 51), ("6r", 79), ("6s", 89),
               ("6t", 65), ("6u", 95), ("6v", 80), ("6w", 81)]
for tag, y in carbon_rows:
    ev = (E[3], T[3])
    nu = "mol_135trimethoxybenzene" if tag == "6u" else f"mol_nu_{tag}"
    parts = [participant("substrate", "mol_7a", "reactant", ev),
             participant("nucleophile", nu, "reagent", ev),
             participant("catalyst", "mol_cat_1", "catalyst", ev, molpct=5,
                         charge_method="inherited_general_condition"),
             participant("solvent", "mol_chcl3", "solvent", ev),
             participant("product", f"mol_{tag}", "product", ev)]
    outs = [outcome("yield", "yield", [f"mol_{tag}"], q(y, "%"), ev, raw=f"{tag} ({y}%)"),
            outcome("diastereoselectivity", "selectivity", [f"mol_{tag}"],
                    "one diastereomer isolated as an enantiomeric pair; numerical dr not reported", ev,
                    scope={"metric": "diastereoselectivity/racemic state", "order": "not numerically reported"})]
    make_reaction(f"scheme4_{tag}", ev, categories=["substrate_scope"],
        location={"page": 3, "scheme": "Scheme 4", "product_label": tag},
        rname="Gold(I)-catalyzed 1,5-allenene cyclization using a carbon nucleophile",
        transform="in situ allenene formation, 5-exo cyclization and carbon-nucleophile trapping",
        parts=parts, step_temp="room temperature", step_time=3,
        operation="Generate 5a in situ from 7a and trap the cyclized benzylic cation with the depicted carbon nucleophile.",
        solvent=("mol_chcl3", "CHCl3"), outcomes=outs,
        class_codes=["C-C", "DB-C-C", "C-H"], class_counts={"C-C": 2}, notes=(racemic_text,),
        unresolved=("Carbon-nucleophile amount, reaction scale, workup, purification and yield measurement basis are not reported in the supplied main-paper pages.",
                    "All products are racemic; absolute configurations and numerical dr are not reported.",
                    "Exact dense structures and complete atom mapping remain unresolved."))

# Text-only 7e result, retained as a partial specific experiment rather than importing unavailable SI.
ev = (E[3], T[3])
parts = [participant("substrate", "mol_7e", "reactant", ev),
         participant("product_6x", "mol_6x", "product", ev),
         participant("product_6x_prime", "mol_6x_prime", "product", ev)]
make_reaction("text_7e_to_6x_6x_prime", ev, categories=["substrate_scope"],
    location={"page": 3, "prose": "paragraph below Scheme 3; Scheme S2 referenced but not supplied"},
    rname="Cyclization of stereogenic propargylic substrate 7e",
    transform="gold-catalyzed cyclization giving two diastereomers",
    parts=parts, operation="React propargylic substrate 7e under conditions not fully stated on the supplied page.",
    outcomes=[outcome("diastereomer_ratio", "selectivity", ["mol_6x", "mol_6x_prime"],
                      {"first": 60, "second": 40, "metric": "dr"}, ev, raw="6x and 6x′; ratio 60:40",
                      scope={"metric": "dr", "order": ["6x", "6x′"]})],
    class_codes=(),
    unresolved=("Scheme S2/SI is not supplied; substrate/product structures, nucleophile, catalyst loading, solvent, time, temperature, yield and measurement method are unresolved.",
                "The 60:40 ratio is retained without assigning absolute configurations."))

# Scheme 1 generic transformations and Figure 1 calculated mechanism are not counted as specific experiments.
make_reaction("scheme1a_literature_allenene_cyclization", (E[2], T[2]), categories=["other"],
    location={"page": 2, "scheme": "Scheme 1a"}, rname="Literature Au(I)-catalyzed cyclization of 1,5-allenenes",
    transform="literature cyclization to bicycle/polycycle", parts=[], class_codes=(), record_type="literature_example",
    notes=("This is cited prior work and is not a current-paper experiment.",),
    unresolved=("Generic R/R′ structures and cited-work conditions are not expanded in this dataset.",))
make_reaction("scheme1b_generic_current_transformation", (E[1], E[2], T[1], T[2]), categories=["other"],
    location={"pages": [1, 2], "scheme": "graphical abstract and Scheme 1b"},
    rname="Presented Au(I)-catalyzed 5-exo cyclization and nucleophile trapping",
    transform="generic 1,5-allenene cyclization with oxygen or carbon nucleophile", parts=[], class_codes=(),
    record_type="generic_scheme", notes=("Generic NuH includes oxygen and carbon nucleophiles; it is not a single executed experiment.",),
    unresolved=("Generic R/R′/Nu variables are not expanded as concrete molecules.",))
fig_ev = (E[3], T[3], E[4], T[4])
fig_parts = [participant("allenene", "mol_5a", "reactant", fig_ev),
             participant("methanol", "mol_meoh", "reagent", fig_ev),
             participant("catalyst", "mol_cat_1", "catalyst", fig_ev),
             participant("product", "mol_6a", "product", fig_ev)]
make_reaction("figure1_calculated_mechanism", fig_ev, categories=["other"],
    location={"page": 4, "figure": "Figure 1"},
    rname="B3LYP-D3BJ(SMD)-calculated catalytic cycle",
    transform="proposed/calculated 5-exo cyclization and methanol trapping", parts=fig_parts,
    class_codes=["C-C", "C-O", "C-H"], record_type="proposed_mechanism",
    notes=("Calculated barriers: 5-exo pathway +11.4 kcal/mol; competing [2+3] pathway +17.4 kcal/mol; nucleophilic-attack transition states differ by >8 kcal/mol.",),
    unresolved=("This is a calculated/proposed mechanism, not a separate preparative experiment; intermediate SMILES and complete atom mapping are not asserted.",))

# Non-reaction source context, including the SI notice without claiming association.
bundle["source_context"] = [
    {"context_id": f"{DOC}:metadata", "document_id": DOC, "kind": "bibliographic_metadata",
     "evidence_ids": [E[1], T[1]], "title": "Controlling the Gold(I)-Catalyzed 1,5-Allenene Reaction: Construction of Fused Rings with Excellent Diastereoselectivity",
     "authors": ["Nina Semleit", "Gebhard Haberhauer"], "doi": "10.1021/acs.orglett.1c03886"},
    {"context_id": f"{DOC}:stereo_assignment", "document_id": DOC, "kind": "stereochemical_discussion",
     "evidence_ids": [E[3], T[3]], "text": "NOESY of 6e and 6g was used to determine product stereochemistry; Figures S7/S8 are not supplied."},
    {"context_id": f"{DOC}:calculation", "document_id": DOC, "kind": "computational_mechanism",
     "evidence_ids": [E[3], E[4], T[3], T[4]], "text": "B3LYP-D3BJ(SMD) mechanism and stereochemical rationale; not a synthetic experiment."},
    {"context_id": f"{DOC}:si_notice", "document_id": DOC, "kind": "supplement_notice",
     "evidence_ids": [E[4], T[4]], "text": "The article advertises a Supporting Information PDF, but no SI file was supplied or explicitly associated in this workspace."},
    {"context_id": f"{DOC}:references", "document_id": DOC, "kind": "bibliography",
     "evidence_ids": [E[4], E[5], T[4], T[5]], "text": "References are bibliographic context and were not imported as current-paper experiments."},
]

(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")

transcription = {
    "document_id": DOC,
    "status": "raw_machine_transcription_needs_review",
    "pages_visually_inspected": [1, 2, 3, 4, 5],
    "table1": [{"entry": r, "catalyst_cell": cat_names[c-1], "yield_5a_percent": y5,
                "yield_6a_percent": y6, "solvent": solv, "methanol_equiv": eq,
                "catalyst_mol_percent": pct} for r, c, y5, y6, solv, eq, pct in table_rows],
    "scheme2": [{"substrate": f"7{x}", "product": f"5{x}", "yield_percent": y}
                for x, y in [("a", 98), ("b", 51), ("c", 98), ("d", 75)]],
    "scheme3": [{"record": tag, "products": [{"molecule_id": m, "yield_percent": y} for m, y in ps]}
                for tag, _, _, ps in oxygen_rows] + [{"record": "8a+8b", "combined_yield_percent": 77}],
    "scheme4": [{"product": tag, "yield_percent": y} for tag, y in carbon_rows],
    "limitations": ["Exact dense structure graphs were not guessed from raster drawings.",
                    "Machine transcription is not source certification."]
}
(RUN / "raw_model/visual_transcription.json").write_text(json.dumps(transcription, ensure_ascii=False, indent=2) + "\n")

stereo_diagnostics = {}
for mid, m in molecules.items():
    smi = st.value(m["smiles"])
    if smi:
        stereo_diagnostics[mid] = st.stereo_inspect(smi)
(RUN / "raw_model/stereo_diagnostics.json").write_text(
    json.dumps(stereo_diagnostics, ensure_ascii=False, indent=2) + "\n")

all_ids = [r["reaction_id"] for r in bundle["reactions"]]
def ids_with(*needles):
    return [rid for rid in all_ids if any(n in rid for n in needles)]

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title, abstract and graphical abstract", "finding": "Current-paper generic Au(I)-catalyzed 1,5-allenene 5-exo cyclization/trapping summary; not a separate executed experiment.", "evidence_ids": [E[1], T[1]]},
            {"region": "two-column introduction", "finding": "Background and cited prior chemistry; no optimization row or current specific reaction arrow.", "evidence_ids": [E[1], T[1]]}],
         "record_ids": [f"{DOC}:scheme1b_generic_current_transformation"],
         "unresolved": ["Graphical-abstract dense structures are not assigned exact machine SMILES."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Scheme 1a-c", "finding": "Scheme 1a is literature precedent; Scheme 1b is the current generic transformation; Scheme 1c contains natural-product structures rather than experiments.", "evidence_ids": [E[2], T[2]]},
            {"region": "Table 1 entries 1-20 and footnotes a-f", "finding": "Twenty independent optimization experiments; 5a and 6a yield columns, solvent/MeOH/catalyst-loading exceptions checked row by row.", "evidence_ids": [E[2], T[2]]},
            {"region": "Scheme 2 arrows 7a-d to 5a-d", "finding": "Four current Au(III)-catalyzed rearrangement experiments at RT, 2 h, CHCl3, with displayed yields.", "evidence_ids": [E[2], T[2]]},
            {"region": "prose beneath Table 1 and Scheme 2", "finding": "Optimization interpretation, separate 98% isolation of 5a, and link from allenenes 5 to later cyclization; used as corroborating evidence without merging distinct conditions.", "evidence_ids": [E[2], T[2]]}],
         "record_ids": ids_with("table1_entry_", "scheme2_", "scheme1a_", "scheme1b_"),
         "unresolved": ["Table 1 yield measurement basis is not stated.", "The heading states 7 mol% counterion, but separate dosing/applicability to each catalyst-cell salt or slash system is not safely resolvable from this page."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Scheme 3 and footnote", "finding": "Fourteen oxygen-nucleophile scope experiments plus one no-added-nucleophile experiment; 6i/6i′ and 6j/6j′ are paired outcomes of their respective single experiments; common RT, 3 h, CHCl3 and racemic-product statement retained.", "evidence_ids": [E[3], T[3]]},
            {"region": "Scheme 4 and footnote", "finding": "Nine carbon-nucleophile scope experiments (6o-w, with separate 6r and 6s entries), common RT, 3 h, CHCl3, and racemic-product statement retained.", "evidence_ids": [E[3], T[3]]},
            {"region": "prose below/around schemes", "finding": "NOESY assignment for 6e/6g, text-only 7e to 6x/6x′ 60:40 result, water-mediated explanation for 8a/8b, competing O/C trapping for 6i/6j, and start of computational-mechanism discussion.", "evidence_ids": [E[3], T[3]]}],
         "record_ids": ids_with("scheme3_", "scheme4_", "text_7e_", "figure1_"),
         "unresolved": ["Scheme S2 and Figures S7-S9 are in an unavailable SI; exact 6x/6x′ structures and stereochemical assignments are unresolved.", "Scope yield measurement bases and most reagent quantities are not stated on the supplied page.", "Dense product and nucleophile structures are not assigned speculative exact SMILES."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "Figure 1a-b and continuation prose", "finding": "B3LYP-D3BJ(SMD) calculated catalytic cycle and stereochemical rationale; stored as proposed mechanism, not as a preparative experiment.", "evidence_ids": [E[4], T[4]]},
            {"region": "conclusion and Associated Content", "finding": "Summary only; SI availability notice observed, but no SI file is explicitly supplied or associated in the workspace.", "evidence_ids": [E[4], T[4]]},
            {"region": "author information, acknowledgments and references 1-6", "finding": "No additional current-paper reaction experiment or independent spectral peak table.", "evidence_ids": [E[4], T[4]]}],
         "record_ids": [f"{DOC}:figure1_calculated_mechanism"],
         "unresolved": ["Calculated intermediates and transition states are not encoded as source-verified molecule graphs."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "references 7-46 across both columns", "finding": "Bibliography only; cited reactions were not imported as current-paper experiments.", "evidence_ids": [E[5], T[5]]}],
         "record_ids": [],
         "unresolved": []}
    ]
}
(RUN / "coverage/pages_0001_0005.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")

print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
                  "source_context": len(bundle["source_context"])}, ensure_ascii=False))
