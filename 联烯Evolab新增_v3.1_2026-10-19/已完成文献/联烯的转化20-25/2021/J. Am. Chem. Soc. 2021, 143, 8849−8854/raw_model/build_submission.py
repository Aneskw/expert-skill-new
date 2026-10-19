#!/usr/bin/env python3
"""Build the page-1--6 source transcription into the v3.1 submission bundle.

This file is a reproducible raw-model artifact.  It deliberately leaves structures
unresolved when the main-text drawing is not safe to transcribe without the SI.
"""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_d6517b57d3476a7a"
RUN = ROOT / "runs/current"
EV = {p: f"{DOC}:page{p}:image" for p in range(1, 7)}
EVT = {p: f"{DOC}:page{p}:text" for p in range(1, 7)}


def f(raw=None, norm=None, pages=(), method="original_extraction", derivation=None,
      status=None):
    evidence = []
    for p in pages:
        evidence.extend([EV[p], EVT[p]])
    return st.fact(raw, norm, evidence, method=method, derivation=derivation,
                   status=status)


def missing():
    return st.fact()


bundle = json.loads((RUN / "submitted.json").read_text())
bundle["molecules"] = []
bundle["reactions"] = []
bundle["source_context"] = []
molecules = {}
stereo_diagnostics = {}


def molecule(mid, label, pages, name=None, *, product=False, chiral=False,
             smiles=None, structure_note=None):
    if mid in molecules:
        existing = molecules[mid]
        occurrence = existing["occurrences"][0]
        occurrence["evidence_ids"] = sorted(set(occurrence["evidence_ids"] +
                                                    [x for p in pages for x in (EV[p], EVT[p])]))
        return mid
    evs = [x for p in pages for x in (EV[p], EVT[p])]
    m = st.new_molecule(mid, DOC, evs)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, pages)
    m["names"]["original"] = f(name or label, name or label, pages)
    m["names"]["english"] = f(name or label, name or label, pages)
    m["cas_numbers"] = missing()
    m["external_properties"] = []
    if smiles:
        diag = st.stereo_inspect(smiles)
        stereo_diagnostics[mid] = {"input_smiles": smiles, **diag}
        m["representation_type"] = "smiles"
        m["smiles"] = f(label, smiles, pages, method="structure_derivation",
                         derivation="SMILES transcribed from the explicitly drawn or named main-text structure; graph parsing is not source verification.")
        m["stereochemistry"]["isomeric_smiles"] = f(label, diag["canonical_isomeric_smiles"], pages,
            method="structure_derivation", derivation="Canonical isomeric SMILES computed with the bundled RDKit stereo diagnostic.")
        m["stereochemistry"]["identity_key"] = f(label, diag["stereo_graph_key"], pages,
            method="structure_derivation", derivation="Stereo graph key computed by standard stereo; it is not a sample identity.")
        m["stereochemistry"]["chirality_status"] = f(label, "unresolved" if mid == "m_1a" else "achiral", pages,
            method="inference", derivation="Graph diagnostic and visible source structure; allene axial support is limited, so 1a remains unresolved.")
        m["stereochemistry"]["geometry_status"] = f(label, "not_applicable", pages,
            method="inference", derivation="No source-specified E/Z element is applicable in this explicitly transcribed molecule.")
        m["stereochemistry"]["representation_status"] = f(label, "partial" if mid == "m_1a" else "encoded", pages,
            method="inference", derivation="Toolkit graph encoding assessed with standard stereo; allene axial stereochemistry support is limited.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing()
        m["generic_template"] = None
        if product or chiral:
            m["stereochemistry"]["atom_centers"] = f(
                "wedge/dash stereochemical structure in the source drawing",
                [{"assignment": "source drawing specifies relative/absolute stereochemistry; CIP labels not transcribed"}],
                pages)
            m["stereochemistry"]["assignment"] = f(
                "enantioenriched syn product or chiral substrate drawn with stereobonds",
                "specified in source drawing; atom-level assignment unresolved", pages)
            m["stereochemistry"]["chirality_status"] = f("specified stereochemical drawing", "specified", pages)
            m["stereochemistry"]["geometry_status"] = f("alkene geometry not assigned as E/Z in the displayed label", "unresolved", pages)
            m["stereochemistry"]["representation_status"] = f("stereochemistry not encoded as SMILES", "not_encoded", pages)
        else:
            m["stereochemistry"]["chirality_status"] = f(label, "unresolved", pages,
                method="inference", derivation="Molecule is identified by the paper label, but its stereochemical graph was not safely transcribed.")
            m["stereochemistry"]["geometry_status"] = f(label, "unresolved", pages,
                method="inference", derivation="No E/Z or ring cis/trans assignment was safely transcribed.")
            m["stereochemistry"]["representation_status"] = f(label, "unresolved", pages,
                method="inference", derivation="Main-text drawing retained as evidence; concrete SMILES intentionally omitted.")
    note = structure_note or "Identity follows the paper label/drawing. Missing CAS and unresolved structural details were not filled from external sources."
    m["identity_note"] = f(note, note, pages)
    molecules[mid] = m
    bundle["molecules"].append(m)
    return mid


# Reusable explicitly named/drawn species.
molecule("m_1a", "1a", [1, 2, 3, 4], "O-benzhydryl 3-alkoxy-1,2-butadiene",
         smiles="C=C=C(C)OC(c1ccccc1)c1ccccc1")
molecule("m_1b", "1b", [4], "O-benzyl 3-alkoxy-1,2-butadiene")
molecule("m_2a", "2a", [2], "p-bromobenzyl alcohol", smiles="OCC1=CC=C(Br)C=C1")
molecule("m_3a", "3a", [2, 4], "p-bromobenzaldehyde", smiles="O=CC1=CC=C(Br)C=C1")
for mid, label, pages, name, smi in [
    ("m_h2ru", "H2Ru(CO)(PPh3)3", [2, 3, 4], "H2Ru(CO)(PPh3)3", None),
    ("m_hclru", "HClRu(CO)(PPh3)3", [2, 3], "HClRu(CO)(PPh3)3", None),
    ("m_rbinap", "(R)-BINAP", [1, 2, 3, 4], "(R)-BINAP", None),
    ("m_sbinap", "(S)-BINAP", [3, 4], "(S)-BINAP", None),
    ("m_lii", "LiI", [1, 2, 3, 4], "lithium iodide", "[Li+].[I-]"),
    ("m_licl", "LiCl", [2], "lithium chloride", "[Li+].[Cl-]"),
    ("m_libr", "LiBr", [2], "lithium bromide", "[Li+].[Br-]"),
    ("m_iproh", "2-PrOH", [1, 2, 3, 4], "2-propanol", "CC(O)C"),
    ("m_d8iproh", "d8-2-propanol", [4], "d8-2-propanol", None),
    ("m_thf", "THF", [1, 2, 3, 4], "tetrahydrofuran", "C1CCOC1"),
    ("m_cpme", "CPME", [2, 3], "cyclopentyl methyl ether", "COC1CCCC1"),
    ("m_dipe", "DIPE", [3, 4], "diisopropyl ether", "CC(C)OC(C)C"),
]:
    molecule(mid, label, pages, name, smiles=smi)


def participant(pid, mid, role, pages, *, molpct=None, equiv=None, ref=None, stage="reaction"):
    p = st.new_participant(pid, mid, role, [x for pg in pages for x in (EV[pg], EVT[pg])], stage=stage)
    if role == "product":
        p["product_role"] = f("main product", "main_product", pages)
    if molpct is not None:
        p["charges"]["mol_percent"] = f(f"{molpct} mol%", {"number": molpct, "unit": "mol%"}, pages)
        p["equivalence_reference"] = f(ref, ref, pages)
    if equiv is not None:
        p["charges"]["equivalents"] = f(f"{equiv} equiv", {"number": equiv, "unit": "equiv"}, pages)
        p["equivalence_reference"] = f(ref, ref, pages)
    return p


def add_changes(r, mode, pages):
    changes = []
    def change(cid, atoms, before, after, kind):
        obj = {"change_id": cid, "atoms": [{"atom_map": None, "element": x} for x in atoms],
               "before_order": before, "after_order": after, "change_type": kind}
        if kind.startswith("hydrogen_"):
            obj["hydrogen_count"] = 1
        changes.append(f(None, obj, pages, method="structure_derivation",
                         derivation="Net bond change inferred from the displayed reactant/product connectivity; atom identities remain unmapped."))
    if mode in {"alcohol", "aldehyde", "coupling"}:
        change("b1", ["C", "C"], 0, 1, "formed")
        change("b2", ["C", "H"], 0, 1, "hydrogen_gain")
        if mode == "aldehyde":
            change("b3", ["O", "H"], 0, 1, "hydrogen_gain")
            change("b4", ["C", "O"], 2, 1, "order_decreased")
        r["classification"]["bond_changes"] = changes
        labels = [
            {"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1"]},
            {"category_number": 1, "code": "C-H", "supporting_change_ids": ["b2"]},
        ]
        if mode == "aldehyde":
            labels.append({"category_number": 24, "code": "O-H", "supporting_change_ids": ["b3"]})
        r["classification"]["labels"] = f(None, labels, pages, method="structure_derivation",
            derivation="Taxonomy labels computed from the explicitly listed partial net changes; no complete atom mapping is claimed.")
        r["classification"]["coverage"] = f(None, "partial", pages, method="inference",
            derivation="Unresolved product structures, implicit-H sources, and omitted coproducts preclude complete atom mapping.")
        r["classification"]["limitations"] = [
            "Partial principal-product classification only; no complete atom mapping is asserted.",
            "Hydrogen donor/acceptor and coproduct atom identities are unresolved without the SI.",
        ]
    else:
        r["classification"]["coverage"] = f(None, "not_applicable" if mode == "no_reaction" else "unresolved", pages,
            method="inference", derivation="No observed bond-forming product" if mode == "no_reaction" else "Insufficient atom-level structure for safe classification.")
        r["classification"]["labels"] = f(None, [], pages, method="inference",
            derivation="No formation label asserted." if mode == "no_reaction" else "Classification deferred pending structure review.")
        r["classification"]["limitations"] = ["No bond-formation labels inferred from an unobserved or unresolved transformation."]


def reaction(rid, pages, record_type="specific_experiment", categories=None, location=None,
             name=None, transformation=None, class_mode="unresolved"):
    evs = [x for p in pages for x in (EV[p], EVT[p])]
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evs)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, pages) if categories is not None else missing()
    r["location"] = f(location or f"PDF page(s) {pages}", {"pages": pages, "region": location}, pages)
    r["reaction_name"] = f(name, name, pages) if name else missing()
    if transformation:
        r["transformation_type"] = f(transformation, transformation, pages, method="inference",
            derivation="Transformation category summarized from the displayed reactants and products.")
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    add_changes(r, class_mode, pages)
    bundle["reactions"].append(r)
    return r


def step(r, sid, pages, participant_ids, temp=None, time=None, operation=None, stage="reaction"):
    s = st.new_step(sid, [x for p in pages for x in (EV[p], EVT[p])], stage=stage)
    s["participant_ids"] = participant_ids
    if operation:
        s["operation"] = f(operation, operation, pages)
    if temp is not None:
        s["temperature"] = f(f"{temp} °C", {"number": temp, "unit": "°C"}, pages)
    if time is not None:
        s["time"] = f(f"{time} h", {"number": time, "unit": "h"}, pages)
    r["steps"].append(s)


def outcome(r, oid, kind, target, pages, value, *, basis=None, scope=None, yield_type=None):
    o = st.new_outcome(oid, kind, [target] if target else [], [x for p in pages for x in (EV[p], EVT[p])])
    raw = value
    if kind in {"yield", "conversion", "recovery"} and isinstance(value, (int, float)):
        norm = {"number": value, "unit": "%"}
        raw = f"{value}%"
    elif isinstance(value, tuple) and value[0] == "percent":
        norm = {"number": value[1], "unit": "%"}
        raw = f"{value[1]}%"
    else:
        norm = value
    o["value"] = f(raw, norm, pages)
    o["product_role"] = f("main product", "main_product", pages) if target else missing()
    o["yield_type"] = f(yield_type, yield_type, pages) if yield_type else missing()
    o["measurement_basis"] = f(basis, basis, pages) if basis else missing()
    o["scope"] = f(scope, scope, pages) if scope is not None else missing()
    r["outcomes"].append(o)


def set_solvent(r, pages, solvent, conc=None, basis=None, ratio=None):
    r["solvent_conditions"]["composition"] = f(solvent, solvent, pages)
    if conc is not None:
        r["solvent_conditions"]["concentration"] = f(f"{conc} M", {"number": conc, "unit": "M"}, pages)
        r["solvent_conditions"]["concentration_basis"] = f(basis, basis, pages)
    if ratio is not None:
        r["solvent_conditions"]["mixture_ratio"] = f(ratio, ratio, pages)


# Page 1: contextual/generic reaction arrows.  Prior literature is not a current-paper experiment.
r = reaction("figure1_prior_allylboration", [1], "literature_example", None,
             "Figure 1 upper auxiliary-based approach", "allene hydroboration–aldehyde allylboration sequence",
             "literature allene hydroboration/aldehyde allylboration", "unresolved")
r["notes"].append(f("Auxiliary-based prior work; cited JACS 2009, 14602. Reported up to >20:1 dr and 85–92% ee.",
                    "prior literature example, not an experiment of this paper", [1]))
r["review"]["unresolved"] = ["Detailed identities and individual prior-literature experiments were not imported from the citation."]

r = reaction("figure1_current_generic", [1], "generic_scheme", None,
             "Figure 1 lower catalytic approach", "Ru-BINAP-catalyzed carbonyl (α-alkoxy)allylation",
             "generic alcohol/aldehyde reductive coupling", "coupling")
r["notes"].append(f("1a + ROH or RCHO; LiI; 2-PrOH for RCHO; THF or DIPE, 75 °C; up to 19:1 dr and 86–99% ee.",
                    "generic overview; not one specific experiment", [1]))
r["review"]["unresolved"] = ["Generic scheme values summarize multiple experiments and are not assigned as a single yield."]

r = reaction("scheme1_transition_state_options", [1], "proposed_mechanism", None,
             "Scheme 1", "chair-like transition-structure alternatives", None, "unresolved")
r["notes"].append(f("Scheme 1 compares σ/π/σ allyl isomerization and anti/syn chair-like transition structures; these are proposed alternatives, not measured reactions.",
                    "proposed mechanistic alternatives", [1]))
r["review"]["unresolved"] = ["No individual experimental outcome belongs to the proposed transition structures."]


# Generic table arrows used by the individual experiment records.
for rid, page, loc, mode in [
    ("table2_alcohol_general", 2, "Table 2 general arrow, alcohol route", "alcohol"),
    ("table2_aldehyde_general", 2, "Table 2 general arrow, aldehyde route", "aldehyde"),
    ("table3_aldehyde_general", 3, "Table 3 general arrow", "aldehyde"),
]:
    rr = reaction(rid, [page], "generic_scheme", None, loc,
                  "ruthenium-catalyzed C–C coupling of alkoxyallene 1a", None, mode)
    rr["review"]["unresolved"] = ["Generic scope arrow; individual yields and selectivities are recorded in separate rows."]


def coupling_participants(r, pages, partner_mid, product_mid, route, *, cat=5, ligand=5,
                          additive="m_lii", additive_pct=10, solvent="m_thf"):
    r["participants"] = [
        participant("p_1a", "m_1a", "reactant", pages, molpct=200, ref=partner_mid),
        participant("p_partner", partner_mid, "reactant", pages),
        participant("p_product", product_mid, "product", pages),
        participant("p_ru", "m_h2ru", "catalyst", pages, molpct=cat, ref=partner_mid),
        participant("p_ligand", "m_rbinap", "ligand", pages, molpct=ligand, ref=partner_mid),
    ]
    if additive:
        r["participants"].append(participant("p_additive", additive, "additive", pages,
                                             molpct=additive_pct, ref=partner_mid))
    if route == "aldehyde":
        r["participants"].append(participant("p_iproh", "m_iproh", "reagent", pages,
                                             molpct=200 if solvent == "m_thf" else 300, ref=partner_mid))
    r["participants"].append(participant("p_solvent", solvent, "solvent", pages))


# Table 1: ten optimization experiments, one row per experiment.
opt = [
    (1, "CPME", 0.4, None, 70, 18, "4:1", 73, "m_h2ru"),
    (2, "CPME", 0.4, "m_licl", 70, 50, "6:1", 59, "m_h2ru"),
    (3, "CPME", 0.4, "m_libr", 70, 76, "7:1", 79, "m_h2ru"),
    (4, "CPME", 0.4, "m_lii", 70, 80, "8.5:1", 86, "m_h2ru"),
    (5, "CPME", 0.4, None, 70, 56, "6:1", 58, "m_hclru"),
    (6, "CPME", 0.4, "m_lii", 70, 65, "8:1", 87, "m_hclru"),
    (7, "THF", 0.4, "m_lii", 70, 73, "9.5:1", 87, "m_h2ru"),
    (8, "THF", 0.4, "m_lii", 75, 75, "9:1", 89, "m_h2ru"),
    (9, "THF", 0.3, "m_lii", 75, 78, "9:1 (10:1 after chromatography)", 90, "m_h2ru"),
    (10, "THF", 0.3, "m_lii", 75, 72, "7.5:1", 88, "m_hclru"),
]
molecule("m_4a", "4a", [2, 4], "syn-sec,tert-diol 4a", product=True)
for n, solv, conc, add, temp, yld, dr, ee, ru in opt:
    rid = f"table1_entry_{n:02d}"
    rr = reaction(rid, [2, 3], "specific_experiment", ["optimization"],
                  f"Table 1 entry {n}", "enantioselective ruthenium-catalyzed C–C coupling of 1a with alcohol 2a",
                  "alcohol-mediated carbonyl reductive coupling", "alcohol")
    rr["general_scheme_ids"] = [f"{DOC}:table2_alcohol_general"]
    rr["participants"] = [
        participant("p_1a", "m_1a", "reactant", [2], molpct=200, ref="m_2a"),
        participant("p_2a", "m_2a", "reactant", [2]),
        participant("p_4a", "m_4a", "product", [2]),
        participant("p_ru", ru, "catalyst", [2], molpct=5, ref="m_2a"),
        participant("p_ligand", "m_rbinap", "ligand", [2], molpct=5, ref="m_2a"),
        participant("p_solvent", "m_cpme" if solv == "CPME" else "m_thf", "solvent", [2]),
    ]
    if add:
        rr["participants"].append(participant("p_additive", add, "additive", [2], molpct=10, ref="m_2a"))
    set_solvent(rr, [2], solv, conc, "2a")
    step(rr, "reaction", [2], [p["participant_id"] for p in rr["participants"]], temp, 24,
         "combine under the Table 1 conditions and heat")
    outcome(rr, "yield", "yield", "m_4a", [2], yld, basis="silica gel chromatography",
            yield_type="isolated")
    outcome(rr, "dr", "selectivity", "m_4a", [2], dr,
            basis="1H NMR of crude reaction mixture" if n != 9 else "1H NMR; 9:1 crude and 10:1 after chromatographic purification",
            scope={"metric": "dr", "major": "syn", "ratio": dr})
    outcome(rr, "ee", "selectivity", "m_4a", [2], ("percent", ee),
            basis="chiral stationary phase HPLC", scope={"metric": "ee", "major_configuration": "unassigned in Table 1"})
    rr["notes"].append(f("Table column header reads '3a (Yield)', while the reaction drawing and discussion identify product 4a.",
                         "source-label conflict retained; outcome targets the drawn/discussed product 4a", [2]))
    rr["review"]["unresolved"] = [
        "Table 1 header/product-label conflict ('3a (Yield)' versus drawn and discussed 4a).",
        "Absolute configuration is drawn but not named in the table; product SMILES remains unresolved.",
        "Detailed operation/workup is in unprovided SI.",
    ]


# Table 2: every alcohol and aldehyde row is a distinct experiment.
table2 = {
    "a": ((78, "10:1", 90), (78, "11:1", 90)),
    "b": ((77, "6:1", 98), (78, "7:1", 94)),
    "c": ((78, "7:1", 94), (84, "9:1", 92)),
    "d": ((75, "8:1", 96), (88, "10:1", 97)),
    "e": ((85, "11:1", 91), (95, "15:1", 90)),
    "f": ((76, "17:1", 94), (82, "17:1", 95)),
    "g": ((82, "7:1", 87), (84, "9:1", 86)),
    "h": ((80, "17:1", 96), (92, "19:1", 98)),
    "i": ((85, "4:1", 90), (82, "8:1", 91)),
    "j": ((76, "10:1", 96), (83, "19:1", 98)),
    "k": ((76, "10:1", 94), (74, "17:1", 97)),
    "l": ((72, "10:1", 89), (79, "10:1", 90)),
    "m": ((51, "19:1", 97), (84, "19:1", 97)),
    "n": ((70, "11:1", 96), (85, "14:1", 95)),
    "o": ((52, "9:1", 95), (87, "10:1", 98)),
    "p": ((91, "13:1", 96), (93, "13:1", 96)),
}
slow = {"g": (10, 48), "m": (10, 48), "n": (5, 48), "o": (10, 48)}
for letter, (alcohol_vals, aldehyde_vals) in table2.items():
    product = "m_4" + letter
    molecule(product, "4" + letter, [2], f"syn-sec,tert-diol 4{letter}", product=True)
    for precursor in ("2", "3"):
        partner = f"m_{precursor}{letter}"
        if partner not in molecules:
            molecule(partner, precursor + letter, [2], f"scope substrate {precursor}{letter}")
        route = "alcohol" if precursor == "2" else "aldehyde"
        vals = alcohol_vals if precursor == "2" else aldehyde_vals
        cat, hours = slow.get(letter, (5, 24)) if precursor == "2" else (5, 24)
        rr = reaction(f"table2_{precursor}{letter}_to_4{letter}", [2], "specific_experiment",
                      ["substrate_scope"], f"Table 2, from {precursor}{letter}",
                      "diastereo- and enantioselective ruthenium-catalyzed C–C coupling",
                      "hydrogen autotransfer" if route == "alcohol" else "2-propanol-mediated reductive coupling",
                      route)
        rr["general_scheme_ids"] = [f"{DOC}:table2_{route}_general"]
        coupling_participants(rr, [2], partner, product, route, cat=cat, ligand=cat)
        rr["reaction_scale"] = f("standard conditions: 0.2 mmol scale",
                                 {"reference_molecule_id": partner, "amount": {"number": 0.2, "unit": "mmol"}}, [2])
        set_solvent(rr, [2], "THF")
        step(rr, "reaction", [2], [p["participant_id"] for p in rr["participants"]], 75, hours,
             "apply Table 2 standard conditions")
        outcome(rr, "yield", "yield", product, [2], vals[0], basis="silica gel chromatography", yield_type="isolated")
        outcome(rr, "dr", "selectivity", product, [2], vals[1], basis="1H NMR of purified material",
                scope={"metric": "dr", "major": "syn", "ratio": vals[1]})
        outcome(rr, "ee", "selectivity", product, [2], ("percent", vals[2]), basis="chiral stationary phase HPLC",
                scope={"metric": "ee", "major_configuration": "drawn product; CIP label not printed in table"})
        rr["review"]["unresolved"] = [
            "Detailed charges, atmosphere, workup and purification parameters are in unprovided SI.",
            "Product stereochemistry is drawn but atom-level CIP/SMILES transcription remains unresolved.",
        ]
        if precursor == "2" and letter in slow:
            rr["notes"].append(f(f"Footnotes: {'10% catalyst; ' if cat == 10 else ''}48 h.",
                                 {"catalyst_loading_mol_percent": cat, "time_h": hours}, [2]))
            if cat == 10:
                rr["review"]["unresolved"].append("Footnote says '10% catalyst'; equal 10 mol% Ru/ligand loading is inferred from the in-situ catalyst assembly and needs SI confirmation.")


# Eq. 1 repeats the 1a/3a experiment from Table 2; add its evidence rather than duplicate it.
dup = next(x for x in bundle["reactions"] if x["reaction_id"] == f"{DOC}:table2_3a_to_4a")
dup["evidence_ids"] = sorted(set(dup["evidence_ids"] + [EV[4], EVT[4]]))
dup["notes"].append(f("Equation (1) repeats the 1a + 3a result: 78% yield, 11:1 dr, 90% ee.",
                      "duplicate source mention for the same experiment", [4]))


# Table 3: 15 achiral-aldehyde scope rows plus 10 catalyst-diastereomer pairs.
table3 = [
    ("q", 72, "10:1", 97, False, False, None), ("r", 77, "15:1", 87, False, False, "1 mmol"),
    ("s", 71, "13:1", 98, False, False, None), ("t", 65, "8:1", 97, True, False, None),
    ("u", 62, "11:1", 96, True, True, None), ("v", 60, "10:1", 99, True, True, None),
    ("w", 64, "15:1", 94, False, False, None), ("x", 66, "14:1", 90, False, False, None),
    ("y", 62, "8:1", 96, False, False, None), ("z", 65, "9:1", 92, False, False, None),
    ("aa", 59, "11:1", 90, False, False, None), ("bb", 67, "10:1", 92, False, False, None),
    ("cc", 59, "12:1", 94, False, False, None), ("dd", 63, "10:1", 91, False, False, None),
    ("ee", 61, "8:1", 89, False, False, None),
    ("ff", 76, "6:1", None, False, False, "(R)-BINAP"), ("iso-ff", 70, "11:1", None, False, False, "(S)-BINAP"),
    ("gg", 77, "10:1", None, False, False, "(R)-BINAP"), ("iso-gg", 66, "13:1", None, False, False, "(S)-BINAP"),
    ("hh", 72, "13:1", None, False, False, "(R)-BINAP; X-ray"), ("iso-hh", 64, "8:1", None, False, False, "(S)-BINAP"),
    ("ii", 59, "12:1", None, False, False, "(R)-BINAP"), ("iso-ii", 64, "14:1", None, False, False, "(S)-BINAP"),
    ("jj", 62, "12:1", None, True, True, "(R)-BINAP"), ("iso-jj", 65, "11:1", None, True, True, "(S)-BINAP"),
]
for label, yld, dr, ee, highcat, slow48, note in table3:
    base = label.removeprefix("iso-")
    partner = "m_3" + base
    product = "m_4" + label.replace("-", "_")
    molecule(partner, "3" + base, [3], f"aliphatic aldehyde 3{base}", chiral=base in {"ff", "gg", "hh", "ii", "jj"})
    molecule(product, "4" + label, [3], f"scope product 4{label}", product=True)
    ligand_mid = "m_sbinap" if label.startswith("iso-") else "m_rbinap"
    cat = 10 if highcat else 7.5
    hours = 48 if slow48 else 24
    rr = reaction(f"table3_3{base}_to_4{label.replace('-', '_')}", [3], "specific_experiment",
                  ["substrate_scope"], f"Table 3, product 4{label}",
                  "diastereo- and enantioselective ruthenium-catalyzed reductive C–C coupling",
                  "2-propanol-mediated reductive coupling", "aldehyde")
    rr["general_scheme_ids"] = [f"{DOC}:table3_aldehyde_general"]
    rr["participants"] = [
        participant("p_1a", "m_1a", "reactant", [3], molpct=200, ref=partner),
        participant("p_partner", partner, "reactant", [3]),
        participant("p_product", product, "product", [3]),
        participant("p_ru", "m_h2ru", "catalyst", [3], molpct=cat, ref=partner),
        participant("p_ligand", ligand_mid, "ligand", [3], molpct=cat, ref=partner),
        participant("p_lii", "m_lii", "additive", [3], molpct=15, ref=partner),
        participant("p_iproh", "m_iproh", "reagent", [3], molpct=300, ref=partner),
        participant("p_dipe", "m_dipe", "solvent", [3]),
    ]
    scale = 1.0 if label == "r" else 0.2
    rr["reaction_scale"] = f(f"{scale:g} mmol scale", {"reference_molecule_id": partner, "amount": {"number": scale, "unit": "mmol"}}, [3])
    set_solvent(rr, [3], "DIPE")
    step(rr, "reaction", [3], [p["participant_id"] for p in rr["participants"]], 75, hours,
         "apply Table 3 conditions")
    outcome(rr, "yield", "yield", product, [3], yld, basis="silica gel chromatography", yield_type="isolated")
    outcome(rr, "dr", "selectivity", product, [3], dr, basis="1H NMR of purified material",
            scope={"metric": "dr", "major": "drawn diastereomer", "ratio": dr})
    if ee is not None:
        outcome(rr, "ee", "selectivity", product, [3], ("percent", ee), basis="chiral stationary phase HPLC",
                scope={"metric": "ee", "major_configuration": "drawn product; CIP label not printed in table"})
    if note:
        rr["notes"].append(f(note, note, [3]))
    rr["review"]["unresolved"] = [
        "Detailed charges, atmosphere, workup and purification parameters are in unprovided SI.",
        "Product stereochemistry is drawn but atom-level CIP/SMILES transcription remains unresolved.",
    ]
    if highcat:
        rr["review"]["unresolved"].append("Footnote says '10% catalyst'; equal 10 mol% Ru/ligand loading is inferred and needs SI confirmation.")
    if ee is None:
        rr["review"]["unresolved"].append("No ee is printed for the chiral-aldehyde matched/mismatched pair; no value was inferred.")


# Scheme 2, step 1: benzhydryl deprotection of 4r.
for mid, label, name in [
    ("m_linaph", "lithium naphthalide", "lithium naphthalide"),
    ("m_4r_deprotected", "deprotected 4r intermediate", "deprotected 4r intermediate"),
    ("m_hgii", "HG-II", "second-generation Hoveyda–Grubbs catalyst"),
    ("m_acrolein", "acrolein", "acrolein"),
    ("m_wittig", "Ph3P=C(Me)CO2Me", "Ph3P=C(Me)CO2Me"),
    ("m_dcm", "DCM", "dichloromethane"),
    ("m_benzene", "PhH", "benzene"),
    ("m_citreodiol", "(−)-citreodiol", "(−)-citreodiol"),
]:
    molecule(mid, label, [3], name, product=mid in {"m_4r_deprotected", "m_citreodiol"}, chiral=mid == "m_citreodiol")
rr = reaction("scheme2_step1_deprotection", [3], "specific_experiment", ["derivatization"],
              "Scheme 2, Step 1", "conversion of 4r toward (−)-citreodiol", "benzhydryl ether deprotection", "unresolved")
rr["participants"] = [participant("p_4r", "m_4r", "reactant", [3]), participant("p_linaph", "m_linaph", "reagent", [3], molpct=500, ref="m_4r"),
                      participant("p_thf", "m_thf", "solvent", [3]), participant("p_product", "m_4r_deprotected", "product", [3])]
set_solvent(rr, [3], "THF", 0.2, "4r")
step(rr, "reaction", [3], [p["participant_id"] for p in rr["participants"]], -20, None, "treat 4r with lithium naphthalide")
outcome(rr, "yield", "yield", "m_4r_deprotected", [3], 78, basis=None, yield_type=None,
        scope={"span": "Scheme 2 Step 1 only"})
rr["review"]["unresolved"] = ["Step-1 intermediate is not separately labeled/drawn; identity is retained as unresolved.",
                               "Time, atmosphere, workup, purification, and yield measurement basis are not reported in the main text."]

rr = reaction("scheme2_step2_citreodiol", [3], "specific_experiment", ["derivatization"],
              "Scheme 2, Step 2", "total synthesis of (−)-citreodiol", "cross-metathesis followed by Wittig olefination", "unresolved")
rr["participants"] = [
    participant("p_start", "m_4r_deprotected", "reactant", [3]), participant("p_hgii", "m_hgii", "catalyst", [3], molpct=7, ref="m_4r_deprotected"),
    participant("p_acrolein", "m_acrolein", "reagent", [3], molpct=1000, ref="m_4r_deprotected"), participant("p_dcm", "m_dcm", "solvent", [3]),
    participant("p_thf", "m_thf", "solvent", [3]), participant("p_wittig", "m_wittig", "reagent", [3], molpct=200, ref="m_4r_deprotected"),
    participant("p_benzene", "m_benzene", "solvent", [3]), participant("p_product", "m_citreodiol", "product", [3]),
]
set_solvent(rr, [3], "DCM–THF then PhH", ratio="4:1 DCM–THF; each stage 0.04 M")
step(rr, "metathesis", [3], ["p_start", "p_hgii", "p_acrolein", "p_dcm", "p_thf"], 60, None, "react with HG-II and acrolein in DCM–THF (4:1, 0.04 M)")
step(rr, "olefination", [3], ["p_wittig", "p_benzene", "p_product"], 90, None, "then add/use Ph3P=C(Me)CO2Me in PhH (0.04 M)")
outcome(rr, "yield", "yield", "m_citreodiol", [3], 86, basis=None, yield_type=None,
        scope={"span": "Scheme 2 Step 2 only; not an overall yield"})
rr["review"]["unresolved"] = ["Stage times, atmosphere, workup, purification, and yield measurement basis are not reported in the main text.",
                               "Atom-mapped multi-operation bond changes remain unresolved; no unsupported full classification was asserted."]


# Scheme 3 isotope-labeling experiments.
for mid, label, name, prod in [
    ("m_d2_2h", "deuterio-2h", "d2-3-furfuryl alcohol deuterio-2h", False),
    ("m_3h", "3h", "3-furfural", False),
    ("m_d4h_i", "deuterio-4h-I", "deuterio-4h-I", True),
    ("m_d4h_ii", "deuterio-4h-II", "deuterio-4h-II", True),
]:
    molecule(mid, label, [4], name, product=prod, chiral=prod)

for tag, partner, product, donor, yld, dr, ee, positions in [
    ("i", "m_d2_2h", "m_d4h_i", None, 90, "11:1", 97, {"Ha": "13% 2H", "Hb": "0% 2H", "Hc": "56% 2H", "Hd": "100% 2H"}),
    ("ii", "m_3h", "m_d4h_ii", "m_d8iproh", 74, "17:1", 98, {"Ha": "7% 2H", "Hb": "0% 2H", "Hc": "71% 2H", "Hd": "0% 2H"}),
]:
    rr = reaction(f"scheme3_labeling_{tag}", [4], "specific_experiment", ["mechanistic_study"],
                  f"Scheme 3 isotope-labeling experiment {tag.upper()}", "isotopic labeling study",
                  "deuterium-labeling carbonyl reductive coupling", "alcohol" if tag == "i" else "aldehyde")
    rr["participants"] = [
        participant("p_1a", "m_1a", "reactant", [4], molpct=200, ref=partner), participant("p_partner", partner, "reactant", [4]),
        participant("p_product", product, "product", [4]), participant("p_ru", "m_h2ru", "catalyst", [4], molpct=5, ref=partner),
        participant("p_ligand", "m_rbinap", "ligand", [4], molpct=5, ref=partner), participant("p_lii", "m_lii", "additive", [4], molpct=10, ref=partner),
        participant("p_thf", "m_thf", "solvent", [4]),
    ]
    if donor:
        rr["participants"].append(participant("p_donor", donor, "reagent", [4], molpct=200, ref=partner))
    set_solvent(rr, [4], "THF")
    step(rr, "reaction", [4], [p["participant_id"] for p in rr["participants"]], 75, 24, "run isotope-labeling coupling under the displayed standard conditions")
    outcome(rr, "yield", "yield", product, [4], yld, basis=None, yield_type=None)
    outcome(rr, "dr", "selectivity", product, [4], dr, basis=None, scope={"metric": "dr", "major": "drawn product", "ratio": dr})
    outcome(rr, "ee", "selectivity", product, [4], ("percent", ee), basis=None, scope={"metric": "ee", "major_configuration": "drawn product"})
    outcome(rr, "deuterium", "other", product, [4], positions, basis=None, scope={"metric": "site-specific deuterium incorporation"})
    rr["review"]["unresolved"] = ["Analytical methods for yield, dr, ee, and site-specific deuterium incorporation are not printed in Scheme 3.",
                                   "Detailed operation/workup and exact isotope-structure SMILES are unresolved without SI."]


# Proposed catalytic cycle is kept distinct from the two real isotope experiments.
rr = reaction("scheme3_proposed_cycle", [4], "proposed_mechanism", None, "Scheme 3 right-hand catalytic cycle",
              "general catalytic mechanism", None, "unresolved")
rr["notes"].append(f("Hydroruthenation gives proposed (Z)-σ-allylruthenium I; aldehyde addition via II gives III; alcohol exchange and β-hydride elimination proceed through IV/V.",
                     "author-proposed mechanism supported by labeling data, not a discrete isolated experiment", [4]))
rr["review"]["unresolved"] = ["Intermediates I–V are proposed and were not isolated; no experimental yield is assigned."]


# Equation (1), new 1b comparison.  The 1a result was merged with Table 2 above.
molecule("m_4a_1b", "4a (from 1b; RO = Bn)", [4], "benzyl-protected 4a analogue from 1b", product=True,
         structure_note="Kept separate from benzhydryl-protected 4a because eq. (1) changes RO from benzhydryl to benzyl despite reusing the 4a label.")
rr = reaction("eq1_1b_3a_to_benzyl_4a", [4], "specific_experiment", ["control"], "Equation (1), from 1b",
              "comparison of alkoxyallene ether substituent", "aldehyde reductive coupling", "aldehyde")
coupling_participants(rr, [4], "m_3a", "m_4a_1b", "aldehyde")
step(rr, "reaction", [4], [p["participant_id"] for p in rr["participants"]], 75, 24, "apply Table 2 standard conditions")
set_solvent(rr, [4], "THF")
outcome(rr, "yield", "yield", "m_4a_1b", [4], 73, basis=None, yield_type=None)
outcome(rr, "dr", "selectivity", "m_4a_1b", [4], "15:1", basis=None, scope={"metric": "dr", "major": "syn"})
outcome(rr, "ee", "selectivity", "m_4a_1b", [4], ("percent", 77), basis=None, scope={"metric": "ee", "major_configuration": "drawn product"})
rr["review"]["unresolved"] = ["Equation (1) does not repeat measurement bases; methods were not inherited without an explicit statement.",
                               "Paper reuses 4a label although RO differs; a separate molecule ID is retained pending SI review."]


# Explicit failed/no-reaction observations in the page-4 prose.
molecule("m_bisnaph_allene", "bis(1-naphthyl)-alkoxyallene", [4], "bis(1-naphthyl)-alkoxyallene")
rr = reaction("bisnaph_allene_no_coupling", [4], "specific_experiment", ["control"], "page 4 prose after eq. (1)",
              "standard-condition coupling test", "no coupling product observed", "no_reaction")
rr["participants"] = [participant("p_allene", "m_bisnaph_allene", "reactant", [4], molpct=200, ref="m_3a"), participant("p_3a", "m_3a", "reactant", [4]),
                      participant("p_ru", "m_h2ru", "catalyst", [4], molpct=5, ref="m_3a"), participant("p_ligand", "m_rbinap", "ligand", [4], molpct=5, ref="m_3a"),
                      participant("p_lii", "m_lii", "additive", [4], molpct=10, ref="m_3a"), participant("p_iproh", "m_iproh", "reagent", [4], molpct=200, ref="m_3a"),
                      participant("p_thf", "m_thf", "solvent", [4])]
set_solvent(rr, [4], "THF")
step(rr, "reaction", [4], [p["participant_id"] for p in rr["participants"]], 75, 24,
     "expose bis(1-naphthyl)-alkoxyallene to 3a under the Table 2 standard conditions")
outcome(rr, "no_reaction", "no_reaction", None, [4], "coupling product was not observed", basis=None,
        scope={"conditions": "standard conditions"})
rr["review"]["unresolved"] = ["Exact charges, detection method, detection limit, and recovered-material balance are not reported."]

molecule("m_mono_allene", "monosubstituted alkoxyallene", [4], "monosubstituted alkoxyallene")
molecule("m_et_allene", "ethyl-substituted allene target", [4], "ethyl-substituted alkoxyallene target")
molecule("m_eti", "ethyl iodide", [4], "ethyl iodide", smiles="CCI")
rr = reaction("ethyl_allene_incomplete_ethylation", [4], "specific_experiment", ["control"], "page 4 prose after eq. (1)",
              "attempted synthesis of ethyl-substituted allene", "lithiation followed by ethyl iodide", "unresolved")
rr["participants"] = [participant("p_start", "m_mono_allene", "reactant", [4]), participant("p_eti", "m_eti", "reagent", [4]), participant("p_target", "m_et_allene", "product", [4])]
step(rr, "attempted_reaction", [4], ["p_start", "p_eti", "p_target"], None, None,
     "lithiate the monosubstituted alkoxyallene, then react with ethyl iodide")
outcome(rr, "result", "other", "m_et_allene", [4], "incomplete ethylation", basis=None,
        scope={"author_comment": "possibly due to competing elimination"})
rr["review"]["unresolved"] = ["Lithium base, quantities, solvent, temperature, time, conversion, product yield, and analytical method are unreported in the main text."]


# Nonreaction source context and page inventory facts.
bundle["source_context"].extend([
    {"context_id": f"{DOC}:page1_natural_product_and_prior_work_context", "document_id": DOC,
     "kind": "background_and_literature_examples", "evidence_ids": [EV[1], EVT[1]],
     "note": "Figure 1 includes natural-product structures and prior-work reaction cartoons; these were not counted as current-paper experiments."},
    {"context_id": f"{DOC}:page4_tertiary_allene_preparation_failure", "document_id": DOC,
     "kind": "unresolved_failed_preparation_statement", "evidence_ids": [EV[4], EVT[4]],
     "note": "Preparation of tertiary allenic ethers could not be achieved because lithiation occurred predominantly at the γ-position; no discrete conditions/product row is provided."},
    {"context_id": f"{DOC}:page5_si_and_crystallography_notice", "document_id": DOC,
     "kind": "associated_content_notice", "evidence_ids": [EV[5], EVT[5]],
     "note": "The article states that SI exists online and lists CCDC data, but no SI file was explicitly supplied to this workspace; it was not associated or imported."},
    {"context_id": f"{DOC}:pages5_6_bibliography", "document_id": DOC,
     "kind": "bibliography", "evidence_ids": [EV[5], EVT[5], EV[6], EVT[6]],
     "note": "References are cited literature rather than experiments performed in this paper."},
])


# Save raw table transcription separately from the schema bundle.
raw = {
    "document_id": DOC,
    "pages_viewed": [1, 2, 3, 4, 5, 6],
    "table1": opt,
    "table2": table2,
    "table2_alcohol_footnotes": slow,
    "table3": table3,
    "decisions": [
        "Eq. (1) 1a/3a result was linked to Table 2 instead of duplicated.",
        "The 1b product was assigned a separate molecule ID because the RO protecting group differs despite label reuse.",
        "SI was not imported; all SI-dependent details remain unresolved.",
        "Main-text scope structures without safe atom-level transcription retain null SMILES and explicit stereo statuses.",
    ],
}
(RUN / "raw_model/visual_transcription.json").write_text(json.dumps(raw, ensure_ascii=False, indent=2) + "\n")

# Literal cross-region mention hints, saved as raw diagnostics rather than promoted facts.
dataset = json.loads((RUN / "dataset.json").read_text())
mentions = st.source_mentions(dataset, ["1a", "1b", "2a", "3a", "4a", "4r", "deuterio-2h", "deuterio-4h"], 220)
(RUN / "raw_model/mentions.json").write_text(json.dumps(mentions, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/stereo_checks.json").write_text(json.dumps(stereo_diagnostics, ensure_ascii=False, indent=2) + "\n")

regions = {
    1: [
        {"region": "title/abstract and graphical abstract", "finding": "Current work overview; no single experimental row.", "evidence_ids": [EV[1], EVT[1]]},
        {"region": "Figure 1 upper and lower reaction cartoons", "finding": "One prior-literature generic sequence and one current-work generic method; natural-product structures are context, not current experiments.", "evidence_ids": [EV[1], EVT[1]]},
        {"region": "Scheme 1", "finding": "Proposed chair-like transition-structure alternatives, not discrete experiments.", "evidence_ids": [EV[1], EVT[1]]},
    ],
    2: [
        {"region": "Table 1", "finding": "Ten optimization experiments transcribed row by row, including footnotes and the 3a/4a header conflict.", "evidence_ids": [EV[2], EVT[2]]},
        {"region": "Table 2", "finding": "Sixteen products, each reported from an alcohol and an aldehyde route (32 distinct experiments), with yield/dr/ee and footnotes.", "evidence_ids": [EV[2], EVT[2]]},
        {"region": "right-column discussion", "finding": "Method context and catalyst composition; detailed experimental procedures are deferred to unavailable SI.", "evidence_ids": [EV[2], EVT[2]]},
    ],
    3: [
        {"region": "Table 3", "finding": "Twenty-five aldehyde scope experiments, including matched/mismatched catalyst pairs and footnotes.", "evidence_ids": [EV[3], EVT[3]]},
        {"region": "Scheme 2", "finding": "Two-step conversion of 4r toward (−)-citreodiol retained as two step-specific records with 78% and 86% outcomes.", "evidence_ids": [EV[3], EVT[3]]},
        {"region": "lower-column discussion", "finding": "Table 1 discussion and Table 2 scope commentary; no duplicate records created.", "evidence_ids": [EV[3], EVT[3]]},
    ],
    4: [
        {"region": "Scheme 3 left", "finding": "Two real deuterium-labeling experiments with yield, dr, ee, and site-specific D incorporation.", "evidence_ids": [EV[4], EVT[4]]},
        {"region": "Scheme 3 right", "finding": "Author-proposed catalytic cycle I–V, kept distinct from experimental labeling records.", "evidence_ids": [EV[4], EVT[4]]},
        {"region": "Equation (1)", "finding": "1a/3a duplicates the Table 2 4a experiment and was merged as additional evidence; 1b comparison is a new record.", "evidence_ids": [EV[4], EVT[4]]},
        {"region": "prose after Equation (1)", "finding": "Bis(1-naphthyl) allene no-coupling control, tertiary-allene preparation failure context, and incomplete ethylation attempt.", "evidence_ids": [EV[4], EVT[4]]},
    ],
    5: [
        {"region": "conclusions/associated content/author information/references 1–6", "finding": "No new current-paper experiment; SI availability and crystallographic accessions are notices only, and cited references were not imported.", "evidence_ids": [EV[5], EVT[5]]},
    ],
    6: [
        {"region": "references 6–24", "finding": "Bibliography only; cited literature was not treated as experiments of this paper.", "evidence_ids": [EV[6], EVT[6]]},
    ],
}
unresolved_by_page = {
    1: ["Generic overview structures and prior-literature details were not expanded into specific experiments."],
    2: ["SI procedures are not supplied; many scope structures and atom-level stereochemical assignments remain unresolved."],
    3: ["Scope-product SMILES/CIP assignments and Scheme 2 times/measurement bases remain unresolved without SI."],
    4: ["Analytical methods for isotope outcomes and exact conditions for failed preparation attempts are not reported in the main text."],
    5: ["The article advertises SI, but no SI file was explicitly supplied or associated."],
    6: ["No current-paper experimental content appears on this page."],
}
coverage = {"document_id": DOC, "semantic_accuracy_verified": False, "pages": []}
for page in range(1, 7):
    evidence_page_ids = {EV[page], EVT[page]}
    record_ids = [r["reaction_id"] for r in bundle["reactions"] if evidence_page_ids.intersection(r["evidence_ids"])]
    coverage["pages"].append({"page": page, "state": "attempted", "regions": regions[page],
                              "record_ids": record_ids, "unresolved": unresolved_by_page[page]})
(RUN / "coverage").mkdir(exist_ok=True)
(RUN / "coverage/pages_0001_0006.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")

tasks_path = RUN / "tasks.json"
tasks = json.loads(tasks_path.read_text())
for task in tasks:
    if task.get("document_id") == DOC and task.get("evidence_id") in set(EVT.values()):
        task["state"] = "attempted"
        task["reason"] = "Full page image and extracted text inspected; detailed region decisions are in coverage/pages_0001_0006.json."
tasks_path.write_text(json.dumps(tasks, ensure_ascii=False, indent=2) + "\n")

st.write(RUN / "submitted.json", bundle)
print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"])}, ensure_ascii=False))
