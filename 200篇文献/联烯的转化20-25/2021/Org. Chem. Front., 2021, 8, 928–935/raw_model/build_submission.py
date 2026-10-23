#!/usr/bin/env python3
"""Session-local transcription of PDF pages 1--6 into the v3.1 bundle.

This file records model decisions; it does not modify the extraction tools or schemas.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st  # noqa: E402


DOC = "doc_755164270d1e5816"
RUN = ROOT / "runs/current"
P = {n: f"{DOC}:page{n}:image" for n in range(1, 7)}
PT = {n: f"{DOC}:page{n}:text" for n in range(1, 7)}


def q(number, unit, **extra):
    return {"number": number, "unit": unit, **extra}


def f(raw=None, normalized=None, ev=(), method="original_extraction", status=None,
      derivation=None, alternatives=()):
    return st.fact(raw, normalized, ev, method=method, status=status,
                   derivation=derivation, alternatives=alternatives)


def inferred(value, ev, reason, raw=None):
    return f(raw, value, ev, method="inference", derivation=reason)


def missing(ev=(), state="unreported"):
    return f(ev=ev, status=state)


bundle = json.loads((RUN / "submitted.json").read_text())
bundle["molecules"] = []
bundle["reactions"] = []
bundle["source_context"] = []
molecules = {}
stereo_diagnostics = {}


SMILES = {
    "1a": "O=C(Nc1cccc2ncccc12)c1ccccc1",
    "2a": "C=C=CC(=O)OCc1ccccc1",
    "2b": "C=C=CC(=O)OC",
    "2c": "C=C=CC(=O)OCC",
    "2d": "C=C=CC(=O)Oc1ccccc1",
    "2e": "C=C=CP(=O)(c1ccccc1)c1ccccc1",
    "2f": "C=C=CCC(=O)OCC",
    "2g": "C=C=C(C)C(=O)OCC",
    "2i": "C=C=Cc1ccccc1",
    "2j": "C=C=Cc1ccc(F)cc1",
    "2k": "C=C=Cc1ccc(Cl)cc1",
    "2l": "C=C=Cc1ccc(OC)cc1",
    "TFE": "OCC(F)(F)F",
    "O2": "O=O",
    "KOAc": "CC(=O)[O-].[K+]",
    "KOTf": "O=S(=O)([O-])C(F)(F)F.[K+]",
    "NaOPiv": "CC(C)(C)C(=O)[O-].[Na+]",
    "CD3OD": "[2H]OC([2H])([2H])[2H]",
    "TEMPO": "CC1(C)CCCC(C)(C)N1[O]",
    "NaH": "[Na+].[H-]",
    "benzaldehyde": "O=Cc1ccccc1",
    "DMF": "CN(C)C=O",
}


NAMES = {
    "1a": "N-(quinolin-8-yl)benzamide",
    "2a": "benzyl buta-2,3-dienoate",
    "2b": "methyl buta-2,3-dienoate",
    "2c": "ethyl buta-2,3-dienoate",
    "2d": "phenyl buta-2,3-dienoate",
    "2e": "diphenyl(propa-1,2-dien-1-yl)phosphine oxide",
    "2f": "ethyl penta-3,4-dienoate",
    "2g": "terminal disubstituted allene 2g",
    "2i": "propa-1,2-dien-1-ylbenzene (phenylallene)",
    "2j": "fluoro-substituted phenylallene 2j",
    "2k": "chloro-substituted phenylallene 2k",
    "2l": "methoxy-substituted phenylallene 2l",
    "Co(acac)2": "cobalt(II) acetylacetonate",
    "Eosin Y": "Eosin Y",
    "RuPS": "[Ru(bpy)3]Cl2",
    "IrPS": "[Ir(dF(CF3)ppy)2(dtbbpy)]PF6",
    "NaOPiv": "sodium pivalate",
    "KOAc": "potassium acetate",
    "KOTf": "potassium triflate",
    "TFE": "2,2,2-trifluoroethanol",
    "O2": "oxygen",
    "air": "air",
    "N2": "nitrogen",
    "CD3OD": "methanol-d4",
    "TEMPO": "2,2,6,6-tetramethyl-1-piperidinyloxy (TEMPO)",
    "NaH": "sodium hydride",
    "benzaldehyde": "benzaldehyde",
    "DMF": "N,N-dimethylformamide",
}


def mol(mid, label=None, ev=(P[2],), name=None, kind="unresolved", stereo=None,
        identity_note=None):
    if mid in molecules:
        m = molecules[mid]
        existing = {(o["document_id"], o["compound_label"]["normalized_value"], tuple(o["evidence_ids"]))
                    for o in m["occurrences"]}
        marker = (DOC, label or mid, tuple(ev))
        if marker not in existing:
            m["occurrences"].append({"document_id": DOC,
                                     "compound_label": f(label or mid, label or mid, ev),
                                     "evidence_ids": list(ev)})
        return mid
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label or mid, label or mid, ev)
    original_name = name or NAMES.get(label or mid) or NAMES.get(mid)
    if original_name:
        m["names"]["original"] = f(original_name, original_name, ev)
        m["names"]["english"] = f(original_name, original_name, ev)
    else:
        m["names"]["original"] = f(label or mid, label or mid, ev)
    m["cas_numbers"] = missing(ev)
    smi = SMILES.get(label or mid) or SMILES.get(mid)
    if smi:
        m["representation_type"] = "smiles"
        m["smiles"] = f({"source": "structure drawing/name", "compound_label": label or mid}, smi, ev,
                         method="structure_derivation",
                         derivation="Candidate connectivity transcribed from the supplied main-text drawing/name; parsing does not verify identity.")
        diag = st.stereo_inspect(smi)
        stereo_diagnostics[mid] = diag
        m["stereochemistry"]["atom_centers"] = inferred(diag["atom_centers"], ev, "RDKit graph diagnostic on the candidate SMILES.")
        m["stereochemistry"]["double_bonds"] = inferred(diag["double_bonds"], ev, "RDKit graph diagnostic on the candidate SMILES; cumulated bonds do not by themselves establish source E/Z.")
        m["stereochemistry"]["axial"] = inferred([], ev, "No axial descriptor is printed; candidate allenes here have a terminal CH2 and the candidate graph does not encode an allene stereogenic axis.")
        m["stereochemistry"]["mixture"] = missing(ev)
        m["stereochemistry"]["assignment"] = inferred("candidate graph has no assigned stereogenic element", ev, "Structure-derived diagnostic, not a source stereochemical assignment.")
        m["stereochemistry"]["chirality_status"] = inferred("achiral", ev, "Candidate graph has no potential atom center and, for candidate allenes, contains a terminal CH2; source identity still needs review.")
        m["stereochemistry"]["geometry_status"] = inferred("not_applicable", ev, "No independently configurable E/Z element is present in the candidate graph.")
        m["stereochemistry"]["representation_status"] = inferred("not_applicable", ev, "No stereochemical element was identified for encoding in this candidate graph.")
        m["stereochemistry"]["isomeric_smiles"] = f({"candidate_smiles": smi}, diag["canonical_isomeric_smiles"], ev,
                                                           method="structure_derivation",
                                                           derivation="Canonical isomeric SMILES calculated from the candidate graph; no source stereo assignment is implied.")
        m["stereochemistry"]["identity_key"] = f(None, diag["stereo_graph_key"], ev,
                                                       method="structure_derivation",
                                                       derivation="Hash of canonical candidate isomeric graph; not a sample identity.")
    else:
        m["representation_type"] = "unresolved"
        if stereo == "product_unspecified":
            m["stereochemistry"]["chirality_status"] = inferred("unspecified", ev, "No R/S or axial assignment is printed in the inspected source drawing.")
            m["stereochemistry"]["geometry_status"] = inferred("unspecified", ev, "An alkene is drawn but no E/Z descriptor is printed; geometry was not invented.")
            m["stereochemistry"]["representation_status"] = inferred("not_encoded", ev, "Definitive source structure/stereo was not encoded as SMILES.")
        else:
            m["stereochemistry"]["chirality_status"] = inferred("unresolved", ev, "Structure or stereochemical status is not safely established from the inspected main-text raster.")
            m["stereochemistry"]["geometry_status"] = inferred("unresolved", ev, "Structure or alkene geometry is not safely established from the inspected main-text raster.")
            m["stereochemistry"]["representation_status"] = inferred("unresolved", ev, "No definitive stereochemical structure representation was created.")
        m["identity_note"] = f(identity_note or "Source label/drawing retained; definitive SMILES, CAS, and stereochemical encoding remain unresolved.",
                               identity_note or "Source label/drawing retained; definitive SMILES, CAS, and stereochemical encoding remain unresolved.", ev)
    molecules[mid] = m
    bundle["molecules"].append(m)
    return mid


def participant(pid, mid, role, ev, stage="reaction"):
    return st.new_participant(pid, mid, role, ev, stage)


def charge(p, field, raw, normalized, ev, reference=None, method="original_extraction"):
    p["charges"][field] = f(raw, normalized, ev, method=method,
                             derivation=("Calculated/inherited only as stated." if method != "original_extraction" else None))
    if field in ("equivalents", "mol_percent"):
        p["equivalence_reference"] = f(reference, reference, ev)


def classify_annulation(ev):
    changes = [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
         "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "N"}],
         "before_order": 0, "after_order": 1, "change_type": "formed"},
    ]
    return changes, [
        {"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1"]},
        {"category_number": 5, "code": "C-N", "supporting_change_ids": ["b2"]},
    ]


def set_classification(r, ev, mode="annulation"):
    if mode == "annulation":
        changes, labels = classify_annulation(ev)
        r["classification"]["bond_changes"] = [f(None, c, ev, method="structure_derivation",
                                                       derivation="Partial net bond-change assignment from the displayed annulation connectivity; exact atom mapping and bond-order shifts remain unresolved.") for c in changes]
        r["classification"]["labels"] = f(None, labels, ev, method="structure_derivation",
                                                derivation="Taxonomy labels calculated from the two partial formed-bond assignments.")
        r["classification"]["coverage"] = inferred("partial", ev, "No complete atom mapping; implicit-H losses and all bond-order changes were not safely resolved.")
        r["classification"]["limitations"] = ["Partial source-drawing classification only; no complete atom map.",
                                                   "C-H activation implies hydrogen loss, but proton destination and complete atom balance are not reported."]
    elif mode == "wittig":
        ch = {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
              "before_order": 0, "after_order": 2, "change_type": "formed"}
        lab = [{"category_number": 13, "code": "DB-C-C", "supporting_change_ids": ["b1"]}]
        r["classification"]["bond_changes"] = [f(None, ch, ev, method="structure_derivation",
                                                       derivation="Displayed olefination forms a new C=C connection; other bond changes and atom origins remain unresolved.")]
        r["classification"]["labels"] = f(None, lab, ev, method="structure_derivation",
                                                derivation="Taxonomy label from the partial formed double bond.")
        r["classification"]["coverage"] = inferred("partial", ev, "Partial displayed-product classification without atom mapping.")
        r["classification"]["limitations"] = ["Byproducts and complete atom mapping are unresolved."]
    else:
        r["classification"]["coverage"] = inferred("unresolved", ev, "No observed net transformation was classified for this contextual mechanism record.")
        r["classification"]["limitations"] = ["Proposed intermediates are not treated as observed products."]


def base_reaction(rid, ev, category, location, reactants, product_ids, record_type="specific_experiment"):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = f(category, category, ev)
    r["location"] = f(location, {"page": min(int(e.split(":page")[1].split(":")[0]) for e in ev), "region": location}, ev)
    r["reaction_name"] = f("annulation of benzamides with allenes", "annulation of benzamides with allenes", ev)
    r["transformation_type"] = inferred("cobalt/photoredox C-H activation and allene annulation", ev,
                                             "Normalized transformation description from title and displayed transformation.")
    for idx, mid in enumerate(reactants):
        r["participants"].append(participant(f"r{idx+1}", mid, "reactant", ev))
    for idx, mid in enumerate(product_ids):
        p = participant(f"p{idx+1}", mid, "product", ev)
        p["product_role"] = f("desired product", "desired product", ev)
        r["participants"].append(p)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
                       "unresolved": ["Machine transcription requires expert source review.",
                                      "No associated SI was supplied; product characterization and detailed procedure are unavailable."]}
    set_classification(r, ev, "annulation")
    bundle["reactions"].append(r)
    return r


def add_agent(r, mid, role, ev, mol_percent=None, equiv=None, volume=None, ref=None, pid=None):
    p = participant(pid or f"a{len(r['participants'])+1}", mid, role, ev)
    if mol_percent is not None:
        charge(p, "mol_percent", f"{mol_percent} mol%", q(mol_percent, "mol%"), ev, ref)
    if equiv is not None:
        charge(p, "equivalents", f"{equiv} equiv", q(equiv, "equiv"), ev, ref)
    if volume is not None:
        charge(p, "volume", f"{volume} mL", q(volume, "mL"), ev)
    r["participants"].append(p)
    return p


def set_standard_conditions(r, ev, time_h=24, atmosphere="O2", light=True, solvent_volume=2.0,
                            co=20, ps="Eosin Y", base="KOTf", base_mol=20):
    ref = r["participants"][0]["molecule_id"]
    if co is not None:
        add_agent(r, mol("Co(acac)2", "Co(acac)2", ev), "catalyst", ev, mol_percent=co, ref=ref, pid="cat_co")
    if ps is not None:
        add_agent(r, mol(ps, ps, ev), "catalyst", ev, mol_percent=5, ref=ref, pid="photocat")
    if base is not None:
        add_agent(r, mol(base, base, ev), "additive", ev, mol_percent=base_mol, ref=ref, pid="base")
    add_agent(r, mol("TFE", "TFE", ev), "solvent", ev, volume=solvent_volume, pid="solvent_tfe")
    add_agent(r, mol(atmosphere, atmosphere, ev), "reagent", ev, pid="atmosphere_agent")
    r["solvent_conditions"]["composition"] = f("2,2,2-trifluoroethanol (TFE)", ["TFE"], ev)
    r["solvent_conditions"]["total_volume"] = f(f"{solvent_volume} mL", q(solvent_volume, "mL"), ev)
    s = st.new_step("reaction", ev)
    s["participant_ids"] = [p["participant_id"] for p in r["participants"]]
    s["operation"] = f("irradiated/stirred under the stated atmosphere", "irradiation" if light else "dark control", ev)
    s["temperature"] = f("room temperature", "room temperature", ev)
    s["time"] = f(f"{time_h} h", q(time_h, "h"), ev)
    s["atmosphere"] = f(atmosphere, atmosphere, ev)
    r["steps"] = [s]
    r["special_conditions"]["light_source"] = f("15 W CFL" if light else "without light",
                                                       {"type": "CFL", "power": q(15, "W")} if light else "none (dark control)", ev)
    r["special_conditions"]["power"] = f("15 W" if light else None, q(15, "W") if light else None, ev,
                                                status=None if light else "not_applicable")


def outcome(r, oid, kind, targets, raw, normalized, ev, basis=None, yield_type=None, scope=None,
            status=None, alternatives=()):
    o = st.new_outcome(oid, kind, targets, ev)
    o["product_role"] = f("desired product", "desired product", ev) if targets else missing(ev)
    o["value"] = f(raw, normalized, ev, status=status, alternatives=alternatives)
    o["yield_type"] = f(yield_type, yield_type, ev) if yield_type else missing(ev)
    o["measurement_basis"] = f(basis, basis, ev) if basis else missing(ev)
    o["scope"] = f(scope, scope, ev) if scope else missing(ev)
    r["outcomes"].append(o)
    return o


# Shared source-identified molecules and generic Scheme 1c entities.
for key in ("1a", "2a", "2b", "2c", "2d", "2e", "2f", "2g", "2i", "2j", "2k", "2l"):
    page = P[4] if key.startswith("2i") or key in ("2j", "2k", "2l") else P[3] if key != "1a" and key != "2a" else P[2]
    mol(key, key, (page,))
generic_benzamide = mol("generic_benzamide_q8", "generic benzamide, Q = 8-quinolyl", (P[2],))
generic_allene = mol("generic_allene_scheme1c", "generic substituted allene", (P[2],))
generic_product = mol("generic_isoquinolinone_scheme1c", "generic isoquinolinone/isoquinolidinone product", (P[2],), stereo="product_unspecified")

r = base_reaction("scheme1c_generic_annulation", (P[2],), ["other"], "Scheme 1c, this work",
                  [generic_benzamide, generic_allene], [generic_product], "generic_scheme")
set_standard_conditions(r, (P[2],), time_h=24, solvent_volume=2.0)
r["procedure"]["original_text"] = f("Co(acac)2, Eosin Y, O2; without metal oxidant; room temperature; 15 W CFL",
                                           "Co(acac)2/Eosin Y/O2, room temperature, 15 W CFL", (P[2],))


# Table 1: every optimization row is a distinct experiment.
table = [
    (1, 10, "RuPS", "O2", "NaOPiv", 20, 32, None, True),
    (2, 10, "IrPS", "O2", "NaOPiv", 20, 5, "lt", True),
    (3, 10, "Eosin Y", "O2", "NaOPiv", 20, 39, None, True),
    (4, 10, "Eosin Y", "O2", "KOAc", 20, 41, None, True),
    (5, 10, "Eosin Y", "O2", "KOTf", 20, 55, None, True),
    (6, 20, "Eosin Y", "O2", "KOTf", 20, 87, None, True),
    (7, 20, "Eosin Y", "air", "KOTf", 20, 53, None, True),
    (8, 20, "Eosin Y", "N2", "KOTf", 20, 5, "lt", True),
    (9, None, "Eosin Y", "O2", "KOTf", 20, 5, "lt", True),
    (10, 20, None, "O2", "KOTf", 20, 16, None, True),
    (11, 20, "Eosin Y", "O2", "KOTf", 20, 14, None, False),
]
mol("3a", "3a", (P[2], P[3]), name="isoquinolinone 3a", stereo="product_unspecified")
for entry, co, ps, atm, base, base_mol, y, qualifier, light in table:
    ev = (P[2],)
    r = base_reaction(f"table1_entry{entry:02d}", ev,
                      ["optimization"] + (["control"] if entry >= 7 else []),
                      f"Table 1 entry {entry}", ["1a", "2a"], ["3a"])
    charge(r["participants"][0], "amount", "0.2 mmol", q(0.2, "mmol"), ev)
    charge(r["participants"][1], "amount", "0.4 mmol", q(0.4, "mmol"), ev)
    r["reaction_scale"] = f("1a (0.2 mmol)", {"reference_molecule_id": "1a", "amount": q(0.2, "mmol")}, ev)
    set_standard_conditions(r, ev, co=co, ps=ps, atmosphere=atm, base=base,
                            base_mol=base_mol, light=light, solvent_volume=2.0)
    raw_y = f"<{y}%" if qualifier == "lt" else f"{y}%"
    val = q(y, "%", qualifier="lt", raw_text=raw_y) if qualifier else q(y, "%")
    basis = "1H NMR using 4-nitrobenzaldehyde as internal standard" if entry in (2, 8, 9) else "isolated"
    outcome(r, "yield_3a", "yield", ["3a"], raw_y, val, ev, basis=basis, yield_type=basis)
    r["procedure"]["original_text"] = f("Table 1 footnote: 1a (0.2 mmol), 2a (0.4 mmol), catalyst/additive as row, TFE (2.0 mL), room temperature, 15 W CFL, 24 h; oxygen unless otherwise stated.",
                                               "Table 1 row plus footnote", ev)
    if entry == 11:
        r["notes"].append(f("c Without light.", "dark control", ev))
    if entry == 6:
        r["evidence_ids"].append(P[3])
        r["notes"].append(f("The same 1a + 2a experiment, conditions, product 3a, and 87% result reappear as Scheme 2 product 3a.",
                               "Table 1 entry 6 and Scheme 2/3a are treated as duplicate source locations for one experiment.",
                               (P[2], P[3])))


# Scheme 2 substrate scope.
scheme2_yields = {"3a": 87, "3b": 85, "3c": 83, "3d": 89, "3e": 73, "3f": 88,
                  "3g": 91, "3h": 86, "3i": 90, "3j": 93, "3k": 83, "3l": 81,
                  "3m": 79, "3n": 77, "3o": 67, "3p": 71}
for i, (prod, y) in enumerate(scheme2_yields.items()):
    substrate = f"1{chr(ord('a') + i)}"
    if prod == "3a":
        continue
    mol(substrate, substrate, (P[3],), name=f"benzamide substrate {substrate}")
    mol(prod, prod, (P[3],), name=f"isoquinolinone {prod}", stereo="product_unspecified")
    r = base_reaction(f"scheme2_{prod}", (P[3],), ["substrate_scope"], f"Scheme 2 product {prod}",
                      [substrate, "2a"], [prod])
    set_standard_conditions(r, (P[3],))
    scope = "two regioisomers, 3:2" if prod == "3n" else "single displayed product; stereochemistry not specified"
    outcome(r, f"yield_{prod}", "yield", [prod], f"{y}%" + (" (3:2)" if prod == "3n" else ""), q(y, "%"),
            (P[3],), scope=scope)
    if prod == "3n":
        outcome(r, "regioisomer_ratio", "selectivity", [prod], "3:2", {"ratio": [3, 2], "metric": "regioisomer ratio"},
                (P[3],), scope="two unresolved regioisomers; printed order not safely encoded")
    r["procedure"]["original_text"] = f("Scheme 2 arrow: Co(acac)2 20 mol%, Eosin Y 5 mol%, TfOK 20 mol%, O2, TFE, r.t., 24 h, 15 W CFL.",
                                               "Scheme 2 displayed conditions", (P[3],))


# Scheme 3 allene scope; arylamide identities are retained as product-specific drawings when labels are not printed per row.
scheme3 = [
    ("4a", 75, "2b"), ("4b", 72, "2b"), ("4c", 65, "2c"), ("4d", 72, "2c"), ("4e", 76, "2c"),
    ("4f", 85, "2d"), ("4g", 82, "2d"), ("4h", 62, "2d"),
    ("4i", 83, "2e"), ("4j", 92, "2e"), ("4k", 85, "2e"), ("4l", 79, "2e"), ("4m", 83, "2e"), ("4n", 91, "2e"),
    ("4o", 61, "2f"), ("4p", 65, "2f"), ("4q", 73, "2g"), ("4r", 65, "2g"), ("4s", 78, "2g"),
]
for prod, y, allene in scheme3:
    amid = f"arylamide_for_{prod}"
    mol(amid, f"drawn arylamide precursor of {prod}", (P[3],), name=f"drawn arylamide precursor of {prod}")
    mol(prod, prod, (P[3], P[4]) if prod in ("4m", "4n") else (P[3],),
        name=f"isoquinolinone {prod}", stereo="product_unspecified",
        identity_note=("Regioselectivity/structure of 4n is reported as confirmed by X-ray (CCDC 2032156); definitive SMILES was not transcribed." if prod == "4n" else None))
    r = base_reaction(f"scheme3_{prod}", (P[3],), ["substrate_scope"], f"Scheme 3 product {prod}",
                      [amid, allene], [prod])
    set_standard_conditions(r, (P[3],))
    outcome(r, f"yield_{prod}", "yield", [prod], f"{y}%", q(y, "%"), (P[3],),
            scope="single displayed regioisomer; no E/Z or enantiomeric descriptor printed")
    r["procedure"]["original_text"] = f("Scheme 3 arrow: Co(acac)2 20 mol%, Eosin Y 5 mol%, TfOK 20 mol%, O2, TFE, r.t., 24 h, 15 W CFL.",
                                               "Scheme 3 displayed conditions", (P[3],))


# Scheme 4 phenylallene scope.
scheme4 = [("5a", 63, "2i"), ("5b", 61, "2i"), ("5c", 55, "2i"), ("5d", 57, "2i"), ("5e", 45, "2i"), ("5f", 45, "2i"),
           ("5g", 56, "2j"), ("5h", 52, "2j"), ("5i", 57, "2k"), ("5j", 51, "2k"), ("5k", 60, "2l")]
for prod, y, allene in scheme4:
    amid = f"arylamide_for_{prod}"
    mol(amid, f"drawn arylamide precursor of {prod}", (P[4],), name=f"drawn arylamide precursor of {prod}")
    mol(prod, prod, (P[4],), name=f"exo-alkene isoquinolinone {prod}", stereo="product_unspecified")
    r = base_reaction(f"scheme4_{prod}", (P[4],), ["substrate_scope"], f"Scheme 4 product {prod}",
                      [amid, allene], [prod])
    set_standard_conditions(r, (P[4],))
    outcome(r, f"yield_{prod}", "yield", [prod], f"{y}%", q(y, "%"), (P[4],),
            scope="single displayed exo-double-bond product; E/Z not printed")
    r["procedure"]["original_text"] = f("Scheme 4 arrow: Co(acac)2 20 mol%, Eosin Y 5 mol%, TfOK 20 mol%, O2, TFE, r.t., 24 h, 15 W CFL.",
                                               "Scheme 4 displayed conditions", (P[4],))


# Scheme 5 scale-up 4m and 4n. The 4n yield conflicts between scheme (81%) and prose (87%).
for prod, substrate, mmol, mass_g, volume, hours, y, amount in [
    ("4m", "arylamide_for_4m", 1.9, 0.6, 12, 20, 79, q(864.6, "mg")),
    ("4n", "arylamide_for_4n", 4.0, 1.17, 20, 40, None, q(1.19, "g")),
]:
    ev = (P[4], PT[4])
    r = base_reaction(f"scheme5_scale_{prod}", ev, ["scale_up"], f"Scheme 5 gram-scale synthesis of {prod}",
                      [substrate, "2e"], [prod])
    charge(r["participants"][0], "amount", f"{mmol} mmol", q(mmol, "mmol"), (P[4],))
    charge(r["participants"][0], "mass", f"{mass_g} g", q(mass_g, "g"), (P[4],))
    charge(r["participants"][1], "equivalents", "2.0 equiv", q(2.0, "equiv"), (PT[4],), reference=substrate)
    r["reaction_scale"] = f(f"{substrate}: {mmol} mmol, {mass_g} g", {"reference_molecule_id": substrate, "amount": q(mmol, "mmol"), "mass": q(mass_g, "g")}, (P[4],))
    set_standard_conditions(r, (P[4],), time_h=hours, solvent_volume=volume)
    if prod == "4m":
        outcome(r, "yield_4m", "yield", ["4m"], "Isolated Yield: 79%", q(79, "%"), (P[4],), basis="isolated", yield_type="isolated")
    else:
        a_scheme = f("Isolated Yield: 81%", q(81, "%"), (P[4],))
        a_prose = f("79% and 87% yield, respectively", q(87, "%"), (PT[4],))
        outcome(r, "yield_4n", "yield", ["4n"], ["Scheme 5: 81%", "prose: 87%"], None, ev,
                basis="isolated (explicit in Scheme 5)", yield_type="isolated", status="conflict",
                alternatives=(a_scheme, a_prose))
        r["notes"].append(f("Source conflict: Scheme 5 prints 81% (1.19 g), while the accompanying prose says the two products were obtained in 79% and 87%, respectively.",
                               "yield conflict retained without choosing", ev))
        r["review"]["unresolved"].append("Resolve the 4n scale-up yield conflict (Scheme 5 81% versus prose 87%) from associated SI or author correction; neither is selected.")
    outcome(r, f"amount_{prod}", "product_amount", [prod],
            "864.6 mg" if prod == "4m" else "1.19 g", amount, (P[4],), basis="isolated", scope="isolated product mass")
    r["procedure"]["original_text"] = f("Scheme 5 displayed scale-up conditions; text states 2.0 equiv preformed allene 2e.",
                                               "Scheme 5 scale-up conditions plus prose allene loading", ev)


# Scheme 5 post-annulation olefination.
mol("6", "6", (P[4],), name="olefin 6", stereo="product_unspecified")
for key in ("NaH", "benzaldehyde", "DMF"):
    mol(key, key, (P[4],))
r = base_reaction("scheme5_4m_to_6", (P[4],), ["derivatization"], "Scheme 5 post-annulation Wittig reaction",
                  ["4m", "benzaldehyde"], ["6"])
r["reaction_name"] = f("post-annulation Wittig reaction", "post-annulation Wittig reaction", (P[4],))
r["transformation_type"] = inferred("olefination", (P[4],), "Normalized from the author's Wittig-reaction description and displayed alkene product.")
add_agent(r, "NaH", "reagent", (P[4],), pid="nah")
add_agent(r, "DMF", "solvent", (P[4],), pid="dmf")
r["solvent_conditions"]["composition"] = f("DMF", ["DMF"], (P[4],))
s = st.new_step("olefination", (P[4],)); s["participant_ids"] = [p["participant_id"] for p in r["participants"]]
s["operation"] = f("NaH, DMF, benzaldehyde", "olefination", (P[4],)); r["steps"] = [s]
outcome(r, "yield_6", "yield", ["6"], "65%", q(65, "%"), (P[4],))
set_classification(r, (P[4],), "wittig")
r["review"]["unresolved"].append("Reagent quantities, temperature, time, and yield measurement basis are not printed in the main-text scheme.")


# Scheme 6 mechanistic experiments.
mol("D5-1a", "[D5]-1a", (P[4],), name="perdeuterated aryl-ring benzamide [D5]-1a")
mol("D5-4i", "[D5]-4i", (P[4],), name="deuterated product [D5]-4i", stereo="product_unspecified")
mol("1q", "1q", (P[4],), name="nitro-substituted benzamide 1q")

r = base_reaction("scheme6_kie", (P[4], PT[4]), ["mechanistic_study"], "Scheme 6a KIE experiment",
                  ["1a", "D5-1a", "2e"], ["4i", "D5-4i"])
set_standard_conditions(r, (P[4],), time_h=(40 / 60), solvent_volume=2.0)
r["steps"][0]["time"] = f("40 min", q(40, "min"), (P[4],))
outcome(r, "kie", "other", ["4i", "D5-4i"], "kH/kD = 1.1", {"metric": "kH/kD", "number": 1.1}, (P[4], PT[4]),
        scope="intermolecular competition of 1a/[D5]-1a (1:1)")
r["notes"].append(f("[D5]-1a/1a (1:1)", {"substrate_ratio": [1, 1]}, (PT[4],)))
r["review"]["unresolved"].append("The method used to derive kH/kD is not reported in the supplied main text.")

r = base_reaction("scheme6_hd_exchange", (P[4], PT[4]), ["mechanistic_study"], "Scheme 6a H/D-exchange control",
                  ["D5-1a", "2e"], ["D5-4i"])
set_standard_conditions(r, (P[4],), time_h=24, solvent_volume=2.0)
outcome(r, "yield_D5_4i", "yield", ["D5-4i"], "80%", q(80, "%"), (P[4],))
outcome(r, "deuterium_retention", "other", ["D5-4i"], "95% D; no H/D exchange observed", q(95, "%"), (P[4], PT[4]),
        scope="aryl deuterium retention/no H-D exchange; analytical method unreported")

mol("CD3OD", "CD3OD", (P[4],))
r = base_reaction("scheme6_cd3od", (P[4], PT[4]), ["mechanistic_study"], "Scheme 6a CD3OD incorporation control",
                  ["1q", "2e"], ["4n"])
set_standard_conditions(r, (P[4],), time_h=24, solvent_volume=2.0)
add_agent(r, "CD3OD", "solvent", (P[4],), pid="cosolvent_cd3od")
r["solvent_conditions"]["composition"] = f("TFE:CD3OD = 4:1", ["TFE", "CD3OD"], (P[4],))
r["solvent_conditions"]["mixture_ratio"] = f("TFE:CD3OD = 4:1", {"components": ["TFE", "CD3OD"], "ratio": [4, 1]}, (P[4],))
outcome(r, "yield_4n", "yield", ["4n"], "65%", q(65, "%"), (P[4],))
outcome(r, "deuterium_incorporation", "other", ["4n"], "no deuterium incorporation", "not observed", (PT[4],),
        scope="product 4n; analytical method and detection limit unreported")

mol("TEMPO", "TEMPO", (P[4],))
r = base_reaction("scheme6_tempo", (P[4], PT[5]), ["mechanistic_study", "control"], "Scheme 6b radical-quenching experiment",
                  ["1a", "2a"], ["3a"])
set_standard_conditions(r, (P[4],), time_h=24, solvent_volume=2.0)
add_agent(r, "TEMPO", "additive", (P[4],), equiv=1.0, ref="1a", pid="tempo")
outcome(r, "yield_3a", "yield", ["3a"], "less than 5%", q(5, "%", qualifier="lt", raw_text="less than 5%"), (P[4],))
r["notes"].append(f("TEMPO 1.0 equiv", q(1.0, "equiv"), (P[4],)))

for key in ("1d", "1e"):
    mol(key, key, (P[4],), name=f"benzamide {key}")
for key, name in (("competition_product_1e", "product derived from electron-rich amide 1e"),
                  ("competition_product_1d", "product derived from electron-poor amide 1d")):
    mol(key, key, (P[4],), name=name, stereo="product_unspecified")
r = base_reaction("scheme6_competition", (P[4], PT[5]), ["mechanistic_study"], "Scheme 6c intermolecular competition",
                  ["1d", "1e", "2a"], ["competition_product_1e", "competition_product_1d"])
set_standard_conditions(r, (P[4],), time_h=8, solvent_volume=2.0)
outcome(r, "product_ratio", "selectivity", ["competition_product_1e", "competition_product_1d"],
        "electron-rich-derived : electron-poor-derived = 1.2 : 1.0", {"ratio": [1.2, 1.0], "metric": "intermolecular product ratio"},
        (P[4], PT[5]), scope="ordered as 1e-derived product : 1d-derived product; printed compound labels are too small to transcribe safely")
r["review"]["unresolved"].append("The product-label glyphs and ratio measurement method are not safely readable/reported; identities are retained by precursor relationship.")


# Page 5 proposed mechanism is formalized separately; CV is source context, not a synthesis reaction.
r = st.new_reaction(f"{DOC}:scheme8_proposed_mechanism", DOC, (P[5], PT[5], P[6], PT[6]))
r["record_type"] = "proposed_mechanism"
r["experiment_categories"] = f("proposed mechanism", ["mechanistic_study"], (P[5], PT[5]))
r["location"] = f("Scheme 8", {"page": 5, "region": "Scheme 8 proposed reaction mechanism"}, (P[5],))
r["reaction_name"] = f("Proposed reaction mechanism", "proposed reaction mechanism", (P[5],))
r["transformation_type"] = inferred("Co(II/III/I)-photoredox annulation cycle", (P[5], PT[5]), "Normalized from author-proposed cycle; not an observed intermediate sequence.")
r["participants"] = [participant("substrate", "1a", "reactant", (P[5],)), participant("allene", "2a", "reactant", (P[5],)), participant("product", "3a", "product", (P[5],))]
r["notes"] = [f("Intermediates A-I and SET/reoxidation cycle are author-proposed; only intermediate A is said to have been observed under stoichiometric conditions, with details deferred to unavailable ESI.",
                    "author-proposed mechanism; intermediate-A evidence unavailable without SI", (P[5], PT[5]))]
r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
                   "unresolved": ["Proposed intermediates are not treated as isolated or fully observed species.",
                                  "The asserted observation of intermediate A cannot be checked because the referenced ESI was not supplied."]}
set_classification(r, (P[5], PT[5]), "unresolved")
bundle["reactions"].append(r)


bundle["source_context"] = [
    {"context_id": f"{DOC}:bibliographic_metadata", "document_id": DOC, "kind": "bibliographic_metadata",
     "evidence_ids": [P[1], PT[1]], "title": "Metallaphotoredox-catalyzed C-H activation: regio-selective annulation of allenes with benzamide",
     "authors": ["Tianlei Li", "Jishun Li", "Zihao Zhu", "Yuyuan Chen", "Xueyao Li", "Qingyun Yang", "Jie Xia", "Wenxuan Zhang", "Chao Zhang", "Weidong Pan", "Song Wu"],
     "journal": "Organic Chemistry Frontiers", "year": 2021, "volume": 8, "pages": "928-935", "doi": "10.1039/d0qo01127d"},
    {"context_id": f"{DOC}:esi_notice", "document_id": DOC, "kind": "supplement_notice", "evidence_ids": [P[1], PT[1]],
     "finding": "The paper says ESI is available and cites CCDC 2032156, but no SI/CIF was supplied or explicitly associated in this workspace."},
    {"context_id": f"{DOC}:scheme1_prior_work", "document_id": DOC, "kind": "literature_context", "evidence_ids": [P[2]],
     "finding": "Scheme 1a summarizes cited cobalt allene-annulation strategies (metal oxidant/electrolysis); these drawings are not experiments of the current paper."},
    {"context_id": f"{DOC}:scheme7_cv", "document_id": DOC, "kind": "mechanistic_measurement", "evidence_ids": [P[5], PT[5]],
     "finding": "Cyclic-voltammetry study: 1a 1.52 V vs SCE; 2a no relevant oxidation potential; in-situ cobalt catalyst 0.94 V vs SCE. General conditions: MeCN, 0.1 M n-Bu4NPF6, 5 mM TfOK, 5 mM Co(acac)2, 5 mM substrates as applicable, 100 mV s-1. This is retained as measurement context, not a synthesis reaction."},
    {"context_id": f"{DOC}:xray_4n", "document_id": DOC, "kind": "structure_statement", "evidence_ids": [P[4], PT[4]],
     "finding": "The structure and regioselectivity of nitro product 4n are stated to be confirmed by X-ray crystallography; CCDC 2032156 is printed. No CIF was associated."},
    {"context_id": f"{DOC}:page6_nonexperimental", "document_id": DOC, "kind": "nonreaction_context", "evidence_ids": [P[6], PT[6]],
     "finding": "Page 6 contains the end of the proposed-mechanism prose, conclusions, author contributions, conflicts, acknowledgements, and references 1-5; no new current-paper experiment is displayed."},
]


# Persist the bundle and raw model inventory.
st.write(RUN / "submitted.json", bundle)
st.write(RUN / "raw_model/stereo_diagnostics.json", stereo_diagnostics)
inventory = {
    "document_id": DOC,
    "pages_in_scope": [1, 2, 3, 4, 5, 6],
    "visual_review": "Every full page image was opened at original detail; source-lineage crops retained for dense regions.",
    "molecule_count": len(bundle["molecules"]),
    "reaction_count": len(bundle["reactions"]),
    "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in bundle["reactions"]),
    "generic_schemes": sum(r["record_type"] == "generic_scheme" for r in bundle["reactions"]),
    "proposed_mechanisms": sum(r["record_type"] == "proposed_mechanism" for r in bundle["reactions"]),
    "source_context_count": len(bundle["source_context"]),
    "record_ids": [r["reaction_id"] for r in bundle["reactions"]],
}
st.write(RUN / "raw_model/extraction_inventory.json", inventory)
prepared_dataset = json.loads((RUN / "dataset.json").read_text())
st.write(RUN / "raw_model/source_mentions.json",
         st.source_mentions(prepared_dataset, ["1a", "2a", "3a", "4n", "isolated yield"], context=180))

ids = inventory["record_ids"]
by_prefix = lambda prefix: [rid for rid in ids if rid.split(":", 1)[1].startswith(prefix)]
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title, abstract, affiliations and ESI/CCDC notice", "finding": "Bibliographic metadata and an explicit ESI/CCDC notice; no SI/CIF was supplied or associated.", "evidence_ids": [P[1], PT[1]]},
            {"region": "Introduction, both columns", "finding": "Background and cited prior cobalt/allene work only; no current-paper experimental arrow or result table.", "evidence_ids": [P[1], PT[1]]}],
         "record_ids": [], "unresolved": ["The paper advertises ESI and CCDC 2032156, but neither file is available in this session."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Scheme 1a", "finding": "Prior cobalt allene-annulation strategies are literature context, not current-paper experiments.", "evidence_ids": [P[2]]},
            {"region": "Scheme 1b-c", "finding": "Design rationale and one generic current-work annulation scheme; no scheme-specific yield.", "evidence_ids": [P[2], PT[2]]},
            {"region": "Table 1, entries 1-11 and footnotes a-c", "finding": "Eleven independent optimization experiments, including atmosphere/catalyst/light controls; row-specific isolated versus NMR bases retained.", "evidence_ids": [P[2], PT[2]]}],
         "record_ids": [f"{DOC}:scheme1c_generic_annulation"] + by_prefix("table1_"),
         "unresolved": ["Definitive structures of the Ru/Ir photocatalysts and product 3a were not encoded from the main-text raster."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "upper discussion text", "finding": "Optimization interpretation, directing-group failures, and substrate-scope prose; failures with weaker directing groups lack individually identified substrates/products and remain prose context.", "evidence_ids": [P[3], PT[3]]},
            {"region": "Scheme 2, products 3a-3p", "finding": "Sixteen independent benzamide-scope experiments; 3n retains the printed 3:2 regioisomer ratio.", "evidence_ids": [P[3]]},
            {"region": "Scheme 3, products 4a-4s", "finding": "Nineteen independent allene/arylamide scope experiments with printed yields and single displayed regioisomers.", "evidence_ids": [P[3], PT[3]]}],
         "record_ids": [f"{DOC}:table1_entry06"] + by_prefix("scheme2_") + by_prefix("scheme3_"),
         "unresolved": ["Scope yields do not state a measurement/isolation basis on this page.", "Most product and product-specific arylamide structures remain drawing-linked without definitive SMILES; no E/Z or enantiomeric descriptor is printed."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "Scheme 4, products 5a-5k", "finding": "Eleven phenylallene-scope experiments forming exo-double-bond products.", "evidence_ids": [P[4], PT[4]]},
            {"region": "Scheme 5", "finding": "Two gram-scale annulations and one post-annulation olefination. The 4n scale-up yield conflicts: scheme 81% versus prose 87%.", "evidence_ids": [P[4], PT[4]]},
            {"region": "Scheme 6a-c", "finding": "Five mechanistic experiments: KIE, H/D-exchange, CD3OD incorporation, TEMPO quench, and intermolecular competition.", "evidence_ids": [P[4], PT[4]]}],
         "record_ids": by_prefix("scheme4_") + by_prefix("scheme5_") + by_prefix("scheme6_"),
         "unresolved": ["The SI-referenced gram-scale procedures and analytical methods are unavailable.", "4n scale-up yield remains an unresolved source conflict (81% figure versus 87% prose).", "The competition-product label glyphs and ratio method are not safely readable/reported; identities are retained by precursor relationship."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "top-left continuation of Scheme 6 discussion", "finding": "Author interpretation of TEMPO and competition results; linked to the existing page-4 experiment records, not duplicated.", "evidence_ids": [P[5], PT[5]]},
            {"region": "Scheme 7 cyclic voltammograms and caption", "finding": "Mechanistic electrochemical measurements retained as source_context, not as synthesis reactions.", "evidence_ids": [P[5], PT[5]]},
            {"region": "Scheme 8 and mechanism prose", "finding": "One author-proposed catalytic mechanism record; proposed intermediates are not treated as isolated synthesis experiments.", "evidence_ids": [P[5], PT[5]]}],
         "record_ids": [f"{DOC}:scheme6_tempo", f"{DOC}:scheme6_competition", f"{DOC}:scheme8_proposed_mechanism"],
         "unresolved": ["Intermediate A is said to have been observed, but the supporting details are in unavailable ESI.", "Cyclic-voltammetry raw data/replicates and detailed analysis are not supplied."]},
        {"page": 6, "state": "attempted", "regions": [
            {"region": "mechanism continuation and Conclusions", "finding": "Continuation of Scheme 8 interpretation and summary; no new arrow or independent experiment.", "evidence_ids": [P[6], PT[6]]},
            {"region": "author contributions, conflicts, acknowledgements and references 1-5", "finding": "Nonexperimental metadata and cited literature only; no current-paper reaction record created.", "evidence_ids": [P[6], PT[6]]}],
         "record_ids": [f"{DOC}:scheme8_proposed_mechanism"],
         "unresolved": ["Pages 7-8 are outside this pass and are not claimed as visually covered."]},
    ]
}
st.write(RUN / "coverage/pages_0001_0006.json", coverage)
print(json.dumps(inventory, ensure_ascii=False, indent=2))
