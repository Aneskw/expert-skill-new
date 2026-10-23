#!/usr/bin/env python3
"""Build the source-derived v3.1 submission for DOI 10.1021/acs.orglett.1c02007.

This is a machine extraction.  It deliberately keeps dense Scheme 2 structures,
unreported procedures, and unreported measurement bases unresolved.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_33fc8f32de8ab299"
RUN = ROOT / "runs/current"
P = {n: [f"{DOC}:page{n}:text", f"{DOC}:page{n}:image"] for n in range(1, 6)}


def f(raw=None, normalized=None, ev=(), method="original_extraction", derivation=None,
      status=None):
    return st.fact(raw, normalized, ev, method=method, derivation=derivation, status=status)


def missing(ev=()):
    return f(ev=ev)


def inferred(normalized, ev, derivation):
    return f(None, normalized, ev, method="inference", derivation=derivation)


def derived(normalized, ev, derivation):
    return f(None, normalized, ev, method="structure_derivation", derivation=derivation)


def q(number, unit, raw_text=None, qualifier=None):
    value = {"number": number, "unit": unit}
    if raw_text is not None:
        value["raw_text"] = raw_text
    if qualifier is not None:
        value["qualifier"] = qualifier
    return value


def make_molecule(mid, label, name, ev, smiles=None, known_achiral=False, identity_note=None,
                  axial_note=None):
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"] = [{
        "document_id": DOC,
        "compound_label": f(label, label, ev),
        "evidence_ids": list(ev),
    }]
    m["names"]["original"] = f(name or label, name or label, ev)
    m["names"]["english"] = f(name or label, name or label, ev)
    m["cas_numbers"] = missing(ev)
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = derived(
            smiles, ev,
            "Candidate connectivity transcribed from the source drawing/name; RDKit parsing is diagnostic only.",
        )
        diag = st.stereo_inspect(smiles)
        m["stereochemistry"]["atom_centers"] = derived(
            diag["atom_centers"], ev, "RDKit graph diagnostic on the candidate SMILES."
        )
        m["stereochemistry"]["double_bonds"] = derived(
            diag["double_bonds"], ev, "RDKit graph diagnostic on the candidate SMILES."
        )
        m["stereochemistry"]["axial"] = (
            f(axial_note, axial_note, ev) if axial_note else
            derived([], ev, "No axial descriptor is encoded by the candidate graph; toolkit support is limited.")
        )
        m["stereochemistry"]["mixture"] = missing(ev)
        m["stereochemistry"]["assignment"] = derived(
            "no stereochemical descriptor encoded", ev,
            "No wedge/dash, E/Z, or axial descriptor is present in the source identity used for this candidate.",
        )
        m["stereochemistry"]["chirality_status"] = derived(
            "achiral" if known_achiral else "unresolved", ev,
            "Achirality is asserted only for simple fully transcribed graphs; otherwise source stereochemistry remains unresolved.",
        )
        m["stereochemistry"]["geometry_status"] = derived(
            "not_applicable" if known_achiral else "unresolved", ev,
            "No applicable specified E/Z geometry is represented in the candidate graph.",
        )
        m["stereochemistry"]["representation_status"] = derived(
            "encoded" if known_achiral else "partial", ev,
            "Connectivity is encoded; dense or toolkit-limited stereochemistry is not promoted beyond the source.",
        )
        if known_achiral:
            m["stereochemistry"]["isomeric_smiles"] = derived(
                diag["canonical_isomeric_smiles"], ev, "Canonical isomeric SMILES from RDKit diagnostic."
            )
            m["stereochemistry"]["identity_key"] = derived(
                diag["stereo_graph_key"], ev, "SHA-256 key of the canonical isomeric graph."
            )
    else:
        m["representation_type"] = "unresolved"
        m["identity_note"] = f(
            identity_note or "Source label/drawing retained; a complete structure was not safely encoded.",
            identity_note or "Source label/drawing retained; a complete structure was not safely encoded.", ev,
        )
        if axial_note:
            m["stereochemistry"]["axial"] = f(axial_note, axial_note, ev)
        for key in ("chirality_status", "geometry_status", "representation_status"):
            m["stereochemistry"][key] = inferred(
                "unresolved", ev,
                "The complete source graph/configuration was not safely encoded; absence of a descriptor is not achirality.",
            )
    return m


def make_template(mid, label, name, ev, fixed, variables):
    m = make_molecule(mid, label, name, ev, identity_note="Generic R-group structure in the source.")
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": f(fixed, fixed, ev),
        "variables": [{
            "label": var,
            "attachment_points": f("source drawing attachment", "source drawing attachment", ev),
            "definition": missing(ev),
            "allowed_values": missing(ev),
            "constraints": missing(ev),
        } for var in variables],
        "instance_mapping": [],
        "reason": f("generic R-group drawing", "Generic source drawing cannot be represented as a concrete SMILES.", ev),
    }
    return m


def participant(pid, mid, role, ev, amounts=None, reference=None, product_role=None):
    p = st.new_participant(pid, mid, role, ev)
    if product_role:
        p["product_role"] = f(product_role, product_role, ev)
    for key, raw, normalized in amounts or []:
        p["charges"][key] = f(raw, normalized, ev)
    if reference:
        p["equivalence_reference"] = f(reference, reference, ev)
    return p


def classification(ev, codes=None, not_applicable=False, limitations=None):
    codes = codes or []
    changes = []
    labels = []
    pairs = {
        "C-C": ("C", "C", 2), "C-O": ("C", "O", 3), "C-H": ("C", "H", 1),
        "C-Cl": ("C", "Cl", 8), "C-N": ("C", "N", 5),
    }
    for index, code in enumerate(codes, 1):
        a, b, number = pairs[code]
        change_id = f"b{index}"
        change = {"change_id": change_id, "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
                  "before_order": 0, "after_order": 1, "change_type": "formed"}
        changes.append(derived(change, ev, "Category-level net bond change inferred from the source reactant/product drawings."))
        labels.append({"category_number": number, "code": code, "supporting_change_ids": [change_id]})
    coverage = "not_applicable" if not_applicable else ("partial" if codes else "unresolved")
    return {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": missing(ev),
        "bond_changes": changes,
        "labels": derived(labels, ev, "Labels follow the supplied 41-category taxonomy from category-level source drawings.") if not not_applicable else derived([], ev, "No observed transformation to classify."),
        "coverage": inferred(coverage, ev, "No atom-complete mapping was asserted; category-level changes only."),
        "limitations": limitations or ["No atom-complete mapping; source drawings support category-level classification only."],
        "model_votes": [],
    }


def base_reaction(rid, ev, record_type, categories, location, name, transformation, codes=None,
                  not_applicable=False):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, ev)
    r["location"] = f(location, {"description": location}, ev)
    r["reaction_name"] = f(name, name, ev)
    r["transformation_type"] = inferred(transformation, ev, "Normalized from the source reaction title/scheme.")
    r["classification"] = classification(ev, codes, not_applicable)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
                   "unresolved": []}
    r["completeness"] = "partial"
    return r


def set_reaction_conditions(r, ev, scale=0.1, tmb_equiv=5, catalyst="m_aucl3",
                            catalyst_mol=5, oxidant="m_pida", oxidant_equiv=1,
                            concentration=0.1, time_h=None, include_1a=True,
                            substrate="m_1a", arene="m_tmb", product=None):
    parts = []
    if include_1a:
        parts.append(participant("p_substrate", substrate, "reactant", ev,
                                 [("amount", f"{scale} mmol", q(scale, "mmol")),
                                  ("equivalents", "1 equiv", q(1, "equiv"))], substrate))
    if arene:
        parts.append(participant("p_arene", arene, "reactant", ev,
                                 [("equivalents", f"{tmb_equiv} equiv", q(tmb_equiv, "equiv"))], substrate))
    if oxidant:
        parts.append(participant("p_oxidant", oxidant, "reagent", ev,
                                 [("equivalents", f"{oxidant_equiv} equiv", q(oxidant_equiv, "equiv"))], substrate))
    if catalyst:
        parts.append(participant("p_catalyst", catalyst, "catalyst", ev,
                                 [("mol_percent", f"{catalyst_mol} mol %", q(catalyst_mol, "mol%"))], substrate))
    parts.append(participant("p_solvent", "m_ch2cl2", "solvent", ev))
    if product:
        parts.append(participant("p_product", product, "product", ev, product_role="desired product"))
    r["participants"] = parts
    r["reaction_scale"] = f(
        f"{scale} mmol {substrate}",
        {"reference_molecule_id": substrate, "amount": q(scale, "mmol")}, ev,
    )
    r["solvent_conditions"]["composition"] = f("CH2Cl2", ["m_ch2cl2"], ev)
    r["solvent_conditions"]["concentration"] = f(f"{concentration} M", q(concentration, "M"), ev)
    r["solvent_conditions"]["concentration_basis"] = f(substrate, substrate, ev)
    step = st.new_step("reaction", ev)
    step["participant_ids"] = [p["participant_id"] for p in parts]
    step["operation"] = f("reaction under the stated conditions", "reaction", ev)
    step["temperature"] = f("room temperature", "room temperature", ev)
    step["time"] = f(f"{time_h} h", q(time_h, "h"), ev) if time_h is not None else missing(ev)
    r["steps"] = [step]
    r["procedure"]["general_procedure_reference"] = f("optimized conditions / Scheme footnote", "main-text scheme/table conditions", ev)
    return r


def add_outcome(r, oid, kind, targets, raw, normalized, ev, role, basis=None, yield_type=None, scope=None):
    o = st.new_outcome(oid, kind, targets, ev)
    o["value"] = f(raw, normalized, ev)
    o["product_role"] = f(role, role, ev)
    o["measurement_basis"] = f(basis, basis, ev) if basis else missing(ev)
    o["yield_type"] = f(yield_type, yield_type, ev) if yield_type else missing(ev)
    o["scope"] = f(scope, scope, ev) if scope else missing(ev)
    r["outcomes"].append(o)
    return o


def register_molecules():
    mols = []
    # Candidate structures used repeatedly.  These are intentionally limited to clear, simple identities.
    mols += [
        make_molecule("m_1a", "1a", "1-phenylpenta-3,4-dien-1-one", P[1] + P[2] + P[3],
                      "C=C=CCC(=O)c1ccccc1", False,
                      axial_note="Terminal allene is drawn without an axial descriptor; its CH2 terminus makes axial chirality inapplicable for this graph."),
        make_molecule("m_2a", "2a", "2-(2,4,6-trimethoxybenzyl)-5-phenylfuran", P[2] + P[3],
                      "COc1cc(OC)c(Cc2ccc(-c3ccccc3)o2)c(OC)c1", True),
        make_molecule("m_3a", "3a", "2-methyl-5-phenylfuran", P[1] + P[2] + P[3],
                      "Cc1ccc(o1)-c1ccccc1", True),
        make_molecule("m_tmb", "arene", "1,3,5-trimethoxybenzene", P[1] + P[2] + P[3],
                      "COc1cc(OC)cc(OC)c1", True),
        make_molecule("m_ch2cl2", "CH2Cl2", "dichloromethane", P[2] + P[3], "ClCCl", True),
    ]
    simple = [
        ("m_aucl3", "AuCl3", "gold(III) chloride"), ("m_ph3paucl", "Ph3PAuCl", "triphenylphosphinegold(I) chloride"),
        ("m_agsbf6", "AgSbF6", "silver hexafluoroantimonate"), ("m_pida", "PIDA", "phenyliodine diacetate"),
        ("m_pifa", "PIFA", "phenyliodine bis(trifluoroacetate)"), ("m_phio", "PhIO", "iodosobenzene"),
        ("m_selectfluor", "Selectfluor", "Selectfluor"), ("m_tbuooh", "t-BuOOH", "tert-butyl hydroperoxide"),
        ("m_d2o", "D2O", "deuterium oxide"),
    ]
    mols += [make_molecule(mid, lab, name, P[2] + P[3]) for mid, lab, name in simple]
    # Scope identities from Scheme 2. Full dense structures are left unresolved rather than guessed.
    for lab in ["1b", "1c", "1d", "1e", "1f", "1g", "1h", "1i"]:
        mols.append(make_molecule(f"m_{lab}", lab, f"beta-ketoallene {lab}", P[3],
                                  identity_note="Scheme 2 beta-ketoallene drawing retained; complete candidate SMILES not asserted.",
                                  axial_note="Terminal allene is drawn without an axial descriptor; other stereochemical features remain unresolved."))
    for lab in ["2b", "2c", "2d", "2e", "2f", "2g", "2h", "2i", "2j", "2k", "2l", "2m", "2n",
                "2o", "2p", "2q", "2r", "2s", "2t", "2u", "2v", "2w", "2x", "2y", "2z", "2aa", "2ab", "2ac"]:
        mols.append(make_molecule(f"m_{lab}", lab, f"2-furylmethylarene {lab}", P[3],
                                  identity_note="Scheme 2 product drawing/label retained; complete graph and stereochemistry not safely encoded."))
    arene_names = {
        "2j": "1,3-dimethoxybenzene", "2k": "2,6-dimethoxytoluene", "2l": "1-bromo-3,5-dimethoxybenzene",
        "2m": "3-methoxy-N,N-dimethylaniline", "2n": "5-methoxycarbazole", "2o": "indole",
        "2p": "N-methylindole", "2q": "N-phenylindole", "2r": "N-benzylindole",
        "2s": "substituted N-methylindole (R3=Me)", "2t": "substituted N-methylindole (R3=OMe)",
        "2u": "substituted N-methylindole (R3=F)", "2v": "substituted N-methylindole (R3=Br)",
        "2w": "substituted N-methylindole (R3=Bpin)", "2x": "phenyl-substituted N-methylindole",
        "2y": "methyl-substituted N-methylindole for 2y", "2z": "methyl-substituted N-methylindole for 2z",
        "2aa": "methyl-substituted N-methylindole for 2aa", "2ab": "3-methyl-N-methylindole",
        "2ac": "3-phenyl-N-methylindole",
    }
    for product, name in arene_names.items():
        mols.append(make_molecule(f"m_arene_{product}", f"Ar-H for {product}", name, P[3],
                                  identity_note="Arene drawing identified in Scheme 2; exact positional graph is retained only by source image where uncertain."))
    mols += [
        make_molecule("m_4", "4", "diaryliodonium salt 4", P[3], identity_note="Scheme 3C proposed/undetected salt drawing."),
        make_molecule("m_5", "5", "2,4,6-trimethoxychlorobenzene", P[3] + P[4]),
        make_molecule("m_3a_d1", "3a-d1", "deuterium-incorporated 2-methyl-5-phenylfuran", P[3] + P[4]),
        make_molecule("m_tmb_d3", "deuterated arene", "deuterium-labeled 1,3,5-trimethoxybenzene (92% D)", P[3] + P[5]),
        make_molecule("m_2a_d4", "2a-d4", "deuterium-labeled product 2a-d4 (91% D)", P[3]),
        make_template("m_generic_1", "1", "generic beta-ketoallene 1", P[2] + P[3],
                      {"source_scaffold": "beta-ketoallene", "defined_connections": "carbonyl/allenyl scaffold"}, ["R1", "R2"]),
        make_template("m_generic_arene", "(Het)Ar-H", "generic electron-rich arene or indole", P[2] + P[3],
                      {"source_scaffold": "electron-rich (hetero)arene bearing C-H"}, ["R3", "R4", "R5"]),
        make_template("m_generic_2", "2", "generic 2-furylmethylarene 2", P[2] + P[3],
                      {"source_scaffold": "2-furylmethylarene product"}, ["R1", "R2", "R3", "R4", "R5"]),
    ]
    return mols


def optimization_records():
    rows = [
        # entry, oxidant, equiv, time, recovered 1a, 2a, 3a, catalyst, mol%, AgSbF6, concentration, arene equiv, scale
        (1, None, None, 24, 0, 38, 41, "m_aucl3", 10, False, 0.4, 5, 0.1),
        (2, None, None, 24, 99, 0, 0, "m_ph3paucl", 10, False, 0.4, 5, 0.1),
        (3, None, None, 24, 0, 0, 85, "m_ph3paucl", 10, True, 0.4, 5, 0.1),
        (4, "m_pida", 1.2, 1.5, 0, 75, 3, "m_aucl3", 10, False, 0.4, 5, 0.1),
        (5, "m_pifa", 1.2, 3, 3, 57, 0, "m_aucl3", 10, False, 0.4, 5, 0.1),
        (6, "m_phio", 1.2, 12, 22, 55, 0, "m_aucl3", 10, False, 0.4, 5, 0.1),
        (7, "m_selectfluor", 1.2, 3, 0, 32, 25, "m_aucl3", 10, False, 0.4, 5, 0.1),
        (8, "m_tbuooh", 1.2, 12, 3, 39, 4, "m_aucl3", 10, False, 0.4, 5, 0.1),
        (9, "m_pida", 1.0, 1.5, 0, 73, 5, "m_aucl3", 10, False, 0.4, 5, 0.1),
        (10, "m_pida", 0.75, 3, 2, 62, 12, "m_aucl3", 10, False, 0.4, 5, 0.1),
        (11, "m_pida", 1.0, 1.5, 0, 86, 4, "m_aucl3", 10, False, 0.1, 5, 0.1),
        (12, "m_pida", 1.0, 3, 0, 87, 5, "m_aucl3", 5, False, 0.1, 5, 0.1),
        (13, "m_pida", 1.0, 1.5, 0, 75, 1, "m_aucl3", 5, False, 0.1, 2, 0.1),
        (14, "m_pida", 1.0, 12, 98, 0, 0, None, None, False, 0.1, 5, 0.1),
        (15, "m_pida", 1.0, 3, 4, 79, 5, "m_aucl3", 5, False, 0.1, 5, 1.0),
    ]
    out = []
    for entry, ox, oxeq, time, rec, y2, y3, cat, catmol, silver, conc, areq, scale in rows:
        ev = P[2] + (P[3] if entry == 12 else [])
        codes = []
        if y2 > 0:
            codes += ["C-O", "C-C"]
        if y3 > 0:
            if "C-O" not in codes:
                codes.append("C-O")
            codes.append("C-H")
        r = base_reaction(f"table1_entry_{entry:02d}", ev, "specific_experiment", ["optimization"],
                          f"PDF page 2, Table 1 entry {entry}",
                          "Optimization of the reaction conditions",
                          "gold-catalyzed tandem cyclization/oxidative arene coupling",
                          codes, not_applicable=not codes)
        set_reaction_conditions(r, ev, scale, areq, cat, catmol, ox, oxeq, conc, time, product="m_2a")
        # 3a is an explicitly tabulated byproduct; 1a is the recovered/remaining substrate column.
        r["participants"].append(participant("p_byproduct", "m_3a", "product", ev, product_role="protonated byproduct"))
        if silver:
            r["participants"].append(participant("p_agsbf6", "m_agsbf6", "additive", ev,
                                                  [("mol_percent", "10 mol %", q(10, "mol%"))], "m_1a"))
        add_outcome(r, "recovery_1a", "recovery", ["m_1a"], f"{rec}%", q(rec, "%"), ev,
                    "starting material column 1a", basis="isolated yield", scope="Table 1 column 1a")
        add_outcome(r, "yield_2a", "yield", ["m_2a"], f"{y2}%", q(y2, "%"), ev,
                    "desired product", basis="isolated yield", yield_type="isolated", scope="Table 1 column 2a")
        add_outcome(r, "yield_3a", "yield", ["m_3a"], f"{y3}%", q(y3, "%"), ev,
                    "protonated byproduct", basis="isolated yield", yield_type="isolated", scope="Table 1 column 3a")
        r["procedure"]["original_text"] = f(
            "Table 1 footnote a and entry row; footnotes c-i applied where marked.",
            "Table 1 entry with footnote-specific catalyst, concentration, arene loading, catalyst omission, or scale.", ev,
        )
        r["review"]["unresolved"] = [
            "Supporting Information was not supplied; addition order, workup, purification, and characterization remain unreported.",
            "The table footnote states isolated yield for the yield columns; no separate analytical details are provided in the main text.",
        ]
        if entry == 12:
            r["notes"].append(f(
                "Scheme 2 repeats 2a at 87% (3 h) under the same standard conditions; evidence merged rather than duplicated.",
                "Scheme 2 repeats 2a at 87% (3 h) under the same standard conditions; evidence merged rather than duplicated.", P[3]))
        out.append(r)
    return out


def scope_records():
    entries = [
        ("2b", 78, 3, "m_1b", "m_tmb"), ("2c", 86, 3, "m_1c", "m_tmb"),
        ("2d", 60, 6, "m_1d", "m_tmb"), ("2e", 81, 3, "m_1e", "m_tmb"),
        ("2f", 84, 3, "m_1f", "m_tmb"), ("2g", 74, 6, "m_1g", "m_tmb"),
        ("2h", 74, 9, "m_1h", "m_tmb"), ("2i", 74, 12, "m_1i", "m_tmb"),
        ("2j", 77, 3, "m_1a", "m_arene_2j"), ("2k", 80, 18, "m_1a", "m_arene_2k"),
        ("2l", 56, 9, "m_1a", "m_arene_2l"), ("2m", 59, 24, "m_1a", "m_arene_2m"),
        ("2n", 80, 9, "m_1a", "m_arene_2n"), ("2o", 63, 24, "m_1a", "m_arene_2o"),
        ("2p", 86, 12, "m_1a", "m_arene_2p"), ("2q", 82, 12, "m_1a", "m_arene_2q"),
        ("2r", 63, 18, "m_1a", "m_arene_2r"), ("2s", 81, 12, "m_1a", "m_arene_2s"),
        ("2t", 58, 12, "m_1a", "m_arene_2t"), ("2u", 75, 12, "m_1a", "m_arene_2u"),
        ("2v", 71, 18, "m_1a", "m_arene_2v"), ("2w", 59, 15, "m_1a", "m_arene_2w"),
        ("2x", 82, 9, "m_1a", "m_arene_2x"), ("2y", 83, 12, "m_1a", "m_arene_2y"),
        ("2z", 79, 15, "m_1a", "m_arene_2z"), ("2aa", 82, 9, "m_1a", "m_arene_2aa"),
        ("2ab", 72, 18, "m_1a", "m_arene_2ab"), ("2ac", 52, 24, "m_1a", "m_arene_2ac"),
    ]
    out = []
    for product, yld, time, substrate, arene in entries:
        r = base_reaction(f"scheme2_scope_{product}", P[3], "specific_experiment", ["substrate_scope"],
                          f"PDF page 3, Scheme 2 product {product}",
                          "Scope of substrates", "gold-catalyzed tandem cyclization/oxidative arene coupling",
                          ["C-O", "C-C"])
        catalyst_mol = 10 if product == "2m" else 5
        set_reaction_conditions(r, P[3], 0.1, 5, "m_aucl3", catalyst_mol, "m_pida", 1,
                                0.1, time, substrate=substrate, arene=arene, product=f"m_{product}")
        add_outcome(r, "yield", "yield", [f"m_{product}"], f"{yld}%", q(yld, "%"), P[3],
                    "desired product", yield_type=None, basis=None, scope=f"Scheme 2 product {product}")
        r["procedure"]["original_text"] = f(
            "Scheme 2 footnote a; footnote b for 2m.",
            "1 (0.1 mmol), arene (5 equiv), PIDA (1 equiv), AuCl3 (5 mol%; 10 mol% for 2m), CH2Cl2 (0.1 M), room temperature.", P[3])
        r["review"]["unresolved"] = [
            "Yield measurement basis is not stated in the main-text Scheme 2 footnote.",
            "Supporting Information was not supplied; addition order, workup, purification, and characterization remain unreported.",
            "Dense source structures remain label/drawing identities where complete candidate SMILES was not safely transcribed.",
        ]
        out.append(r)
    return out


def mechanism_experiments():
    out = []
    # Scheme 3A
    r = base_reaction("scheme3A_ph3paucl_pida", P[3], "specific_experiment", ["mechanistic_study"],
                      "PDF page 3, Scheme 3A", "Mechanistic study: Au(I) catalyst with PIDA",
                      "gold-catalyzed tandem cyclization/oxidative arene coupling", ["C-O", "C-C"])
    set_reaction_conditions(r, P[3], 0.1, 5, "m_ph3paucl", 5, "m_pida", 1, 0.1, 3, product="m_2a")
    add_outcome(r, "yield_2a", "yield", ["m_2a"], "73%", q(73, "%"), P[3], "desired product")
    r["review"]["unresolved"] = ["Yield basis, workup, purification, and full procedure are unreported without SI."]
    out.append(r)
    # Scheme 3B
    r = base_reaction("scheme3B_3a_control", P[3] + P[4], "specific_experiment", ["control", "mechanistic_study"],
                      "PDF page 3, Scheme 3B and page 4 discussion", "Control: attempted conversion of 3a to 2a",
                      "negative control", not_applicable=True)
    set_reaction_conditions(r, P[3], 0.1, 5, "m_aucl3", 5, "m_pida", 1, 0.1, 3,
                            substrate="m_3a", arene="m_tmb", product="m_2a")
    add_outcome(r, "recovery_3a", "recovery", ["m_3a"], "89%", q(89, "%"), P[3] + P[4], "recovered starting material")
    add_outcome(r, "no_2a", "no_reaction", ["m_2a"], "not transformed to 2a", "not transformed to 2a", P[3] + P[4], "expected product not formed")
    r["review"]["unresolved"] = ["Recovery measurement basis and detection limit for absent 2a are unreported."]
    out.append(r)
    # Scheme 3C
    r = base_reaction("scheme3C_arene_pida_aucl3", P[3] + P[4], "specific_experiment", ["control", "mechanistic_study"],
                      "PDF page 3, Scheme 3C and page 4 continuation", "Control: arene/PIDA/AuCl3 without beta-ketoallene",
                      "arene chlorination control", ["C-Cl"])
    r["participants"] = [
        participant("p_arene", "m_tmb", "reactant", P[3], [("amount", "0.2 mmol", q(0.2, "mmol"))]),
        participant("p_pida", "m_pida", "reagent", P[3], [("amount", "0.1 mmol", q(0.1, "mmol"))]),
        participant("p_aucl3", "m_aucl3", "reagent", P[3], [("amount", "0.1 mmol", q(0.1, "mmol"))]),
        participant("p_solvent", "m_ch2cl2", "solvent", P[3]),
        participant("p_4", "m_4", "product", P[3], product_role="undetected diaryliodonium salt"),
        participant("p_5", "m_5", "product", P[3] + P[4], product_role="observed chlorinated product"),
    ]
    r["reaction_scale"] = f("0.2 mmol arene", {"reference_molecule_id": "m_tmb", "amount": q(0.2, "mmol")}, P[3])
    r["solvent_conditions"]["composition"] = f("CH2Cl2", ["m_ch2cl2"], P[3])
    r["solvent_conditions"]["concentration"] = f("0.1 M", q(0.1, "M"), P[3])
    r["solvent_conditions"]["concentration_basis"] = f("1,3,5-trimethoxybenzene", "m_tmb", P[3])
    step = st.new_step("reaction", P[3]); step["participant_ids"] = [p["participant_id"] for p in r["participants"]]
    step["operation"] = f("reaction without 1a", "reaction without 1a", P[3]); step["temperature"] = f("room temperature", "room temperature", P[3])
    r["steps"] = [step]
    add_outcome(r, "product_5_amount", "product_amount", ["m_5"], "0.077 mmol", q(0.077, "mmol"), P[3], "observed product")
    add_outcome(r, "yield_5", "yield", ["m_5"], "77%", q(77, "%"), P[3] + P[4], "observed product")
    add_outcome(r, "salt_4_not_obtained", "no_reaction", ["m_4"], "4, not obtained", "not obtained", P[3] + P[4], "undetected proposed product")
    r["review"]["unresolved"] = ["Reaction time, yield basis, and detection method/limit for 4 are not reported in the main text."]
    out.append(r)
    # Scheme 3D
    r = base_reaction("scheme3D_d2o_labeling", P[3] + P[4], "specific_experiment", ["mechanistic_study"],
                      "PDF page 3, Scheme 3D and page 4 discussion", "D2O labeling/deuterodeauration experiment",
                      "gold-catalyzed cyclization and deuterodeauration", ["C-O", "C-H"])
    set_reaction_conditions(r, P[3], 0.1, 0, "m_aucl3", 5, None, None, 0.1, 3,
                            substrate="m_1a", arene=None, product="m_3a_d1")
    r["participants"].insert(1, participant("p_d2o", "m_d2o", "reagent", P[3],
                                             [("equivalents", "10 equiv", q(10, "equiv"))], "m_1a"))
    r["steps"][0]["participant_ids"] = [p["participant_id"] for p in r["participants"]]
    add_outcome(r, "yield_3a_d1", "yield", ["m_3a_d1"], "27%", q(27, "%"), P[3], "deuterated product")
    add_outcome(r, "d_incorporation", "selectivity", ["m_3a_d1"], "84% D", q(84, "%"), P[3], "deuterium incorporation", scope="D incorporation at the labeled product position")
    r["review"]["unresolved"] = ["Yield and D-incorporation measurement methods are referred to unavailable SI."]
    out.append(r)
    # Scheme 3E is a paired parallel KIE study; both arrows belong to the same comparative experiment.
    r = base_reaction("scheme3E_parallel_kie", P[3] + P[4], "specific_experiment", ["mechanistic_study"],
                      "PDF page 3, Scheme 3E and page 4 discussion", "Parallel H/D kinetic isotope effect experiment",
                      "parallel isotope-effect comparison", ["C-O", "C-C"])
    r["participants"] = [
        participant("p_1a_h", "m_1a", "reactant", P[3], [("amount", "0.1 mmol", q(0.1, "mmol")), ("equivalents", "1 equiv", q(1, "equiv"))], "m_1a"),
        participant("p_arene_h", "m_tmb", "reactant", P[3], [("equivalents", "5 equiv", q(5, "equiv"))], "m_1a"),
        participant("p_1a_d", "m_1a", "reactant", P[3], [("amount", "0.1 mmol", q(0.1, "mmol")), ("equivalents", "1 equiv", q(1, "equiv"))], "m_1a"),
        participant("p_arene_d", "m_tmb_d3", "reactant", P[3], [("equivalents", "5 equiv", q(5, "equiv"))], "m_1a"),
        participant("p_pida", "m_pida", "reagent", P[3], [("equivalents", "1 equiv", q(1, "equiv"))], "m_1a"),
        participant("p_aucl3", "m_aucl3", "catalyst", P[3], [("mol_percent", "5 mol %", q(5, "mol%"))], "m_1a"),
        participant("p_solvent", "m_ch2cl2", "solvent", P[3]),
        participant("p_2a_h", "m_2a", "product", P[3], product_role="H-arm product"),
        participant("p_2a_d", "m_2a_d4", "product", P[3], product_role="D-arm product"),
    ]
    r["reaction_scale"] = f("0.1 mmol 1a per parallel arm", {"reference_molecule_id": "m_1a", "amount": q(0.1, "mmol")}, P[3])
    r["solvent_conditions"]["composition"] = f("CH2Cl2", ["m_ch2cl2"], P[3])
    r["solvent_conditions"]["concentration"] = f("0.1 M", q(0.1, "M"), P[3])
    r["solvent_conditions"]["concentration_basis"] = f("1a", "m_1a", P[3])
    step1 = st.new_step("parallel_reactions_1h", P[3]); step1["participant_ids"] = [p["participant_id"] for p in r["participants"]]
    step1["operation"] = f("parallel H and D reactions; 1 h sampling point", "parallel reactions sampled at 1 h", P[3]); step1["temperature"] = f("room temperature", "room temperature", P[3]); step1["time"] = f("1 h", q(1, "h"), P[3])
    step3 = st.new_step("parallel_reactions_3h", P[3]); step3["participant_ids"] = [p["participant_id"] for p in r["participants"]]
    step3["operation"] = f("parallel H and D reactions; 3 h sampling point", "parallel reactions sampled at 3 h", P[3]); step3["temperature"] = f("room temperature", "room temperature", P[3]); step3["time"] = f("3 h", q(3, "h"), P[3])
    r["steps"] = [step1, step3]
    for oid, kind, target, raw, val, role, scope in [
        ("h_yield_1h", "yield", "m_2a", "time = 1 h: 30%", q(30, "%"), "H-arm product", "1 h"),
        ("h_recovery_1h", "recovery", "m_1a", "53% of 1a", q(53, "%"), "H-arm recovered substrate", "1 h"),
        ("h_yield_3h", "yield", "m_2a", "3 h: 87%", q(87, "%"), "H-arm product", "3 h"),
        ("d_yield_1h", "yield", "m_2a_d4", "time = 1 h: 29%", q(29, "%"), "D-arm product", "1 h"),
        ("d_recovery_1h", "recovery", "m_1a", "55% of 1a", q(55, "%"), "D-arm recovered substrate", "1 h"),
        ("d_yield_3h", "yield", "m_2a_d4", "3 h: 84%", q(84, "%"), "D-arm product", "3 h"),
    ]:
        add_outcome(r, oid, kind, [target], raw, val, P[3], role, scope=scope)
    add_outcome(r, "kie", "other", ["m_2a", "m_2a_d4"], "kH/kD = 1.0", {"metric": "kH/kD", "value": 1.0}, P[3] + P[4], "parallel KIE comparison", scope="H/D kinetic isotope effect")
    add_outcome(r, "arene_d_level", "selectivity", ["m_tmb_d3"], "92% D", q(92, "%"), P[3], "input isotopic enrichment", scope="deuterated arene")
    add_outcome(r, "product_d_level", "selectivity", ["m_2a_d4"], "91% D", q(91, "%"), P[3], "product isotopic enrichment", scope="2a-d4")
    r["review"]["unresolved"] = ["Individual rate data, fitting method, analytical method, and yield/recovery bases are in unavailable SI or unreported in the main text."]
    out.append(r)
    return out


def contextual_records():
    out = []
    r = base_reaction("scheme1D_generic_current_work", P[1] + P[2] + P[3], "generic_scheme", [],
                      "PDF pages 1-3, graphical abstract and Scheme 1D/Scheme 2 header",
                      "Oxidative coupling of beta-ketoallenes and arenes or indoles",
                      "generic gold-catalyzed tandem cyclization/oxidative coupling", ["C-O", "C-C"])
    r["participants"] = [participant("p_g1", "m_generic_1", "reactant", P[2]), participant("p_ga", "m_generic_arene", "reactant", P[2]),
                         participant("p_gp", "m_generic_2", "product", P[2], product_role="generic product"),
                         participant("p_cat", "m_aucl3", "catalyst", P[2]), participant("p_ox", "m_pida", "reagent", P[2])]
    r["review"]["unresolved"] = ["Generic R-group scheme; not a concrete experiment and no yield is assigned."]
    out.append(r)
    for suffix, title, codes in [
        ("scheme1B_a_ref7", "Prior-work Au-catalyzed tandem reaction from ref 7", ["C-O", "C-C"]),
        ("scheme1B_b_ref8", "Prior-work Au-catalyzed tandem reaction from ref 8", ["C-O", "C-C"]),
        ("scheme1C_ref11", "Prior-work oxidative oxi/amino-arylation from ref 11", ["C-C", "C-O", "C-N"]),
    ]:
        lr = base_reaction(suffix, P[2], "literature_example", [], f"PDF page 2, Scheme 1 {suffix}", title,
                           "cited-work schematic transformation", codes)
        lr["review"]["unresolved"] = ["Cited prior work shown as a generic schematic; no current-paper experiment or yield assigned."]
        out.append(lr)
    r = base_reaction("scheme3F_proposed_mechanism", P[3] + P[4], "proposed_mechanism", [],
                      "PDF page 3, Scheme 3F and page 4 discussion", "Possible catalytic cycles (routes A and B)",
                      "author-proposed Au(I)/Au(III) catalytic mechanism", not_applicable=True)
    r["participants"] = [participant("p_g1", "m_generic_1", "reactant", P[3]),
                         participant("p_ga", "m_generic_arene", "reactant", P[3]),
                         participant("p_gp", "m_generic_2", "product", P[3], product_role="generic product"),
                         participant("p_cat", "m_aucl3", "catalyst", P[3]), participant("p_ox", "m_pida", "reagent", P[3])]
    r["notes"] = [f(
        "Authors propose enolization, Au-allene coordination, 5-exocyclization, olefin isomerization, arene C-H auration, reductive elimination, and Au(I)/Au(III) turnover; route B may predominate but route A is not ruled out.",
        "Authors propose routes A and B and state that route B might be major while route A cannot be ruled out.", P[3] + P[4])]
    r["review"]["unresolved"] = ["Proposed intermediates I-V/A/III are mechanistic hypotheses, not observed products; no observed-transformation labels assigned."]
    out.append(r)
    return out


def main():
    data_path = RUN / "dataset.json"
    data = json.loads(data_path.read_text())
    assert not data.get("experiments"), "Refuse to mix pre-existing experiments into this source extraction."
    doc = next(d for d in data["documents"] if d["document_id"] == DOC)
    doc["title"] = "Gold-Catalyzed Tandem Oxidative Coupling Reaction between β-Ketoallenes and Electron-Rich Arenes to 2-Furylmethylarenes"
    doc["supplementary_material"] = {
        "status": "unknown",
        "note": "The main text states SI is available, but no SI file was supplied or explicitly associated in this workspace."
    }
    st.write(data_path, data)
    bundle = json.loads((RUN / "submitted.json").read_text())
    bundle["documents"] = data["documents"]
    bundle["molecules"] = register_molecules()
    bundle["reactions"] = optimization_records() + scope_records() + mechanism_experiments() + contextual_records()
    bundle["source_context"] = [
        {"context_id": f"{DOC}:page1:graphical_abstract", "document_id": DOC, "evidence_ids": P[1],
         "kind": "duplicate_generic_summary", "note": "Graphical abstract duplicates the current generic transformation and is not a separate experiment."},
        {"context_id": f"{DOC}:page2:scheme1A", "document_id": DOC, "evidence_ids": P[2],
         "kind": "application_examples", "note": "Bioactivity examples motivate the scaffold; they are not current synthesis experiments."},
        {"context_id": f"{DOC}:page4:si_notice", "document_id": DOC, "evidence_ids": P[4],
         "kind": "si_cross_reference", "note": "SI is mentioned but no SI file is supplied/associated."},
        {"context_id": f"{DOC}:page5:references", "document_id": DOC, "evidence_ids": P[5],
         "kind": "references_only", "note": "Reference list; cited reactions are not imported as current-paper experiments."},
    ]
    st.write(RUN / "submitted.json", bundle)
    diagnostics = {}
    for m in bundle["molecules"]:
        smi = st.value(m["smiles"])
        if smi:
            diagnostics[m["molecule_id"]] = st.stereo_inspect(smi)
    st.write(RUN / "raw_model/stereo_diagnostics.json", diagnostics)
    inventory = {
        "document_id": DOC,
        "source_pages_reviewed": [1, 2, 3, 4, 5],
        "table1_rows": 15,
        "scheme2_displayed_products": 29,
        "scheme2_new_scope_records": 28,
        "scheme2_duplicate": {"product": "2a", "merged_into": f"{DOC}:table1_entry_12"},
        "scheme3_real_studies": ["3A", "3B", "3C", "3D", "3E"],
        "scheme3_proposed_mechanism": "3F",
        "prior_work_schematics": ["1B-a", "1B-b", "1C"],
        "limitations": ["No SI supplied", "Dense Scheme 2 structures not guessed", "Most Scheme 2 yield bases unreported"],
    }
    st.write(RUN / "raw_model/extraction_inventory.json", inventory)
    print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"])}, ensure_ascii=False))


if __name__ == "__main__":
    main()
