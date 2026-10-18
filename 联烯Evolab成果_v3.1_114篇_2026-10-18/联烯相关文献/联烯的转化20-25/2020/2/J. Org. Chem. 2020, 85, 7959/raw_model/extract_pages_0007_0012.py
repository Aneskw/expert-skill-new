#!/usr/bin/env python3
"""Extend the v3.1 machine bundle with source PDF pages 7--12.

The six original page images were inspected in the current session.  This
program preserves the resulting machine transcription/interpretation while
leaving every reaction in needs_review state.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SKILL = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SKILL))

from standard_records import (  # noqa: E402
    fact,
    new_molecule,
    new_outcome,
    new_participant,
    new_reaction,
    new_step,
    source_mentions,
    stereo_inspect,
)

DOC = "doc_f3d393b861fb548d"
RUN = ROOT / "runs/current"
RAW = RUN / "raw_model"
OUT = RUN / "submitted.json"


def tids(*pages: int) -> list[str]:
    out = []
    for page in pages:
        out.extend([f"{DOC}:page{page}:text", f"{DOC}:page{page}:image"])
    return out


def F(raw, normalized=None, ev=(), method="original_extraction", derivation=None,
      status=None, verification=None, alternatives=()):
    if normalized is None and raw is not None:
        normalized = raw
    return fact(raw=raw, normalized=normalized, evidence=ev, method=method,
                derivation=derivation, status=status, verification=verification,
                alternatives=alternatives)


def Q(number, unit, **extra):
    q = {"number": number, "unit": unit}
    q.update(extra)
    return q


def inferred(raw, normalized, ev, derivation):
    return F(raw, normalized, ev, method="inference", derivation=derivation)


# Frozen v3.1 base produced by the prior page-1--6 pass.  Keeping this source
# bundle separate makes the page-7--12 extraction reproducible after merge.
BASE = RAW / "base_pages_0001_0006.json"
bundle = json.loads(BASE.read_text())
# The fresh standard-init template registers the page-7--12 image evidence
# blocks.  Retain the prior records while replacing the document registry with
# this superset; otherwise valid new image evidence would appear dangling.
fresh_template = json.loads((RAW / "init_template_pages_0007_0012.json").read_text())
bundle["documents"] = fresh_template["documents"]
molecules = {m["molecule_id"]: m for m in bundle["molecules"]}
reactions = {r["reaction_id"]: r for r in bundle["reactions"]}


def merge_evidence(obj, ev):
    if isinstance(obj, dict):
        if "evidence_ids" in obj and isinstance(obj["evidence_ids"], list):
            for e in ev:
                if e not in obj["evidence_ids"]:
                    obj["evidence_ids"].append(e)
        for v in obj.values():
            merge_evidence(v, ev)
    elif isinstance(obj, list):
        for v in obj:
            merge_evidence(v, ev)


def extend_molecule_evidence(mid, *pages):
    ev = tids(*pages)
    m = molecules[mid]
    occ = m["occurrences"][0]
    for e in ev:
        if e not in occ["evidence_ids"]:
            occ["evidence_ids"].append(e)
    for key in ("compound_label",):
        for e in ev:
            if e not in occ[key]["evidence_ids"]:
                occ[key]["evidence_ids"].append(e)
    return ev


def add_molecule(mid, label, name, pages, smiles=None, *, chirality_note=None,
                 mixture=None, source_assignment=None, axial_note=None,
                 representation="smiles"):
    ev = tids(*pages)
    if mid in molecules:
        extend_molecule_evidence(mid, *pages)
        return molecules[mid]
    m = new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, ev)
    m["names"]["original"] = F(name, name, ev)
    m["names"]["english"] = F(name, name, ev)
    m["cas_numbers"] = fact(raw=None, normalized=None, evidence=ev, status="unreported")
    m["identity_note"] = F(
        "Machine candidate derived from the inspected source name/structure; CAS was not reported.",
        "machine candidate; CAS unreported", ev)
    if representation == "unresolved" or not smiles:
        m["representation_type"] = "unresolved"
        m["smiles"] = fact(raw=None, normalized=None, evidence=ev, status="unreported")
        m["stereochemistry"]["chirality_status"] = inferred(
            "source entity is not fully encoded",
            "mixture" if mixture else "unresolved", ev,
            "The inspected source does not support one complete stereochemical molecular graph.")
        m["stereochemistry"]["geometry_status"] = inferred(
            "source entity is not fully encoded", "unresolved", ev,
            "No complete geometry assignment is asserted.")
        m["stereochemistry"]["representation_status"] = inferred(
            "no complete stereochemical encoding", "not_encoded", ev,
            "The source entity is retained without a complete candidate isomeric SMILES.")
        if mixture:
            m["stereochemistry"]["mixture"] = F(mixture, mixture, ev)
        if source_assignment:
            m["stereochemistry"]["assignment"] = F(source_assignment, source_assignment, ev)
    else:
        diag = stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = F(
            "candidate derived from the inspected source name/structure", smiles, ev,
            method="structure_derivation",
            derivation="Connectivity and any explicit isotope/stereo tags were derived from the source name and drawing; RDKit parsing is diagnostic only.")
        m["stereochemistry"]["identity_key"] = F(
            "standard stereo graph diagnostic", diag["stereo_graph_key"], ev,
            method="structure_derivation", derivation="Computed with standard_records.stereo_inspect.")
        m["stereochemistry"]["isomeric_smiles"] = F(
            "canonical isomeric candidate", diag["canonical_isomeric_smiles"], ev,
            method="structure_derivation", derivation="Computed from the candidate SMILES; source verification is not implied.")
        m["stereochemistry"]["atom_centers"] = F(
            "stereo diagnostic", diag["atom_centers"], ev, method="structure_derivation",
            derivation="Zero-based atom indices refer to the candidate SMILES.")
        m["stereochemistry"]["double_bonds"] = F(
            "double-bond diagnostic", diag["double_bonds"], ev, method="structure_derivation",
            derivation="Candidate-graph double bonds were inspected; carbonyls are not E/Z assignments.")
        potential_atom = any(x["type"].startswith("Atom_") and x["specified"] == "Unspecified"
                             for x in diag["potential_stereo"])
        potential_bond = any(x["type"].startswith("Bond_") and x["specified"] == "Unspecified"
                             for x in diag["potential_stereo"])
        specified_atom = any(x["type"].startswith("Atom_") and x["specified"] == "Specified"
                             for x in diag["potential_stereo"])
        if mixture:
            chir = "mixture"
            m["stereochemistry"]["mixture"] = F(mixture, mixture, ev)
        elif source_assignment or specified_atom:
            chir = "specified"
        elif chirality_note or potential_atom:
            chir = "unspecified"
        else:
            chir = "achiral"
        if source_assignment:
            m["stereochemistry"]["assignment"] = F(source_assignment, source_assignment, ev)
        elif chirality_note:
            m["stereochemistry"]["assignment"] = F(chirality_note, chirality_note, ev,
                method="structure_derivation",
                derivation="Potential stereogenicity was assessed from the source connectivity; no configuration was assigned.")
        else:
            m["stereochemistry"]["assignment"] = inferred(
                "no source stereochemical assignment", "no assigned tetrahedral configuration", ev,
                "The candidate graph and source naming do not specify a tetrahedral configuration.")
        if axial_note:
            m["stereochemistry"]["axial"] = F(axial_note, axial_note, ev,
                method="structure_derivation", derivation="Allene-axis substitution was inspected separately from toolkit atom-center output.")
        else:
            m["stereochemistry"]["axial"] = inferred(
                "no allene/atropisomeric axis in candidate", "not_applicable", ev,
                "The candidate contains no stereogenic allene or atropisomeric axis.")
        geom = "unspecified" if potential_bond else "not_applicable"
        if source_assignment and "E" in source_assignment:
            geom = "specified"
        rep = "encoded" if specified_atom and not mixture and not chirality_note else (
            "partial" if mixture or chirality_note or potential_atom or potential_bond else "not_applicable")
        m["stereochemistry"]["chirality_status"] = inferred(
            "source/stereo diagnostic", chir, ev,
            "Status combines the inspected source name/drawing with the standard stereo diagnostic; missing tags are never interpreted as racemic.")
        m["stereochemistry"]["geometry_status"] = inferred(
            "alkene geometry diagnostic", geom, ev,
            "Only potentially stereogenic non-carbonyl double bonds require E/Z assignment.")
        m["stereochemistry"]["representation_status"] = inferred(
            "candidate encoding status", rep, ev,
            "Potential unassigned or mixture stereochemistry makes the representation partial.")
    bundle["molecules"].append(m)
    molecules[mid] = m
    return m


def replace_existing_molecule(mid, name, label, pages, smiles, assignment=None,
                              chirality_note=None, axial_note=None):
    old = molecules[mid]
    ev = list(dict.fromkeys(old["occurrences"][0]["evidence_ids"] + tids(*pages)))
    idx = bundle["molecules"].index(old)
    del molecules[mid]
    new = add_molecule(mid, label, name, pages, smiles,
                       source_assignment=assignment,
                       chirality_note=chirality_note, axial_note=axial_note)
    new["occurrences"][0]["evidence_ids"] = ev
    for f in (new["occurrences"][0]["compound_label"], new["names"]["original"],
              new["names"]["english"], new["smiles"], new["identity_note"]):
        f["evidence_ids"] = ev
    bundle["molecules"].remove(new)
    bundle["molecules"][idx] = new
    molecules[mid] = new
    return new


def base_reaction(tag, pages, location, name, transformation, categories=("other",),
                  unresolved=()):
    ev = tids(*pages)
    rid = f"{DOC}:{tag}"
    if rid in reactions:
        raise ValueError(f"duplicate reaction {rid}")
    r = new_reaction(rid, DOC, ev)
    r["record_type"] = "specific_experiment"
    r["experiment_categories"] = F(list(categories), list(categories), ev)
    r["location"] = F(location, location, ev)
    r["reaction_name"] = F(name, name, ev)
    r["transformation_type"] = F(transformation, transformation, ev)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [],
                   "unresolved": list(unresolved) or ["Machine extraction requires expert source review."]}
    bundle["reactions"].append(r)
    reactions[rid] = r
    return r, ev


def participant(r, mid, role, ev, *, mass=None, amount=None, volume=None,
                eq=None, molpct=None, ref=None, stage="reaction", product_role=None,
                method="original_extraction"):
    p = new_participant(f"{r['reaction_id']}:p{len(r['participants'])+1}", mid, role, ev, stage=stage)
    if product_role:
        p["product_role"] = F(product_role, product_role, ev, method=method)
    for key, val in (("mass", mass), ("amount", amount), ("volume", volume)):
        if val:
            p["charges"][key] = F(val[0], Q(val[1], val[2]), ev, method=method)
    if eq is not None:
        p["charges"]["equivalents"] = F(f"{eq} equiv", Q(eq, "equiv"), ev, method=method)
        p["equivalence_reference"] = F(f"relative to {ref}", ref, ev, method=method)
    if molpct is not None:
        p["charges"]["mol_percent"] = F(f"{molpct} mol%", Q(molpct, "mol%"), ev, method=method)
        p["equivalence_reference"] = F(f"relative to {ref}", ref, ev, method=method)
    r["participants"].append(p)
    return p


def step(r, ev, operation, *, temp=None, time=None, atmosphere=None, pressure=None,
         stage="reaction", pids=None, method="original_extraction"):
    s = new_step(f"{r['reaction_id']}:s{len(r['steps'])+1}", ev, stage=stage)
    s["participant_ids"] = pids if pids is not None else [
        p["participant_id"] for p in r["participants"] if p["stage"] == stage]
    s["operation"] = F(operation, operation, ev, method=method)
    if temp:
        s["temperature"] = F(temp[0], temp[1], ev, method=method)
    if time:
        s["time"] = F(time[0], time[1], ev, method=method)
    if atmosphere:
        s["atmosphere"] = F(atmosphere, atmosphere, ev, method=method)
    if pressure:
        s["pressure"] = F(pressure, pressure, ev, method=method)
    r["steps"].append(s)
    return s


def outcome(r, ev, kind, targets, raw, normalized, *, role=None, yield_type=None,
            basis=None, scope=None, basis_missing=False):
    o = new_outcome(f"{r['reaction_id']}:o{len(r['outcomes'])+1}", kind, targets, ev)
    o["value"] = F(raw, normalized, ev)
    if role:
        o["product_role"] = F(role, role, ev)
    if yield_type:
        o["yield_type"] = F(yield_type, yield_type, ev)
    if basis:
        o["measurement_basis"] = F(basis, basis, ev)
    elif basis_missing:
        o["measurement_basis"] = fact(raw=None, normalized=None, evidence=ev, status="unreported")
    if scope:
        o["scope"] = F(scope, scope, ev)
    r["outcomes"].append(o)
    return o


def set_solvent(r, ev, mids, raw, volume=None, concentration=None, basis=None,
                method="original_extraction"):
    r["solvent_conditions"]["composition"] = F(raw, list(mids), ev, method=method)
    if volume:
        r["solvent_conditions"]["total_volume"] = F(volume[0], Q(volume[1], volume[2]), ev, method=method)
    if concentration:
        r["solvent_conditions"]["concentration"] = F(concentration[0], Q(concentration[1], concentration[2]), ev, method=method)
    if basis:
        r["solvent_conditions"]["concentration_basis"] = F(basis, basis, ev, method=method)


def classify(r, ev, formed, limitation):
    changes = []
    labels = []
    category = {"C-H": 1, "C-C": 2, "C-O": 3, "C-S": 4, "C-N": 5,
                "DB-C-C": 13, "DB-O-S": 18, "O-H": 24, "O-S": 26,
                "N-H": 29, "S-H": 34}
    for i, (a, b, order, code) in enumerate(formed, 1):
        cid = f"b{i}"
        chg_type = "formed" if order == 1 else "order_increased"
        changes.append(F(
            "partial source-derived bond-change candidate",
            {"change_id": cid, "atoms": [{"atom_map": None, "element": a},
                                           {"atom_map": None, "element": b}],
             "before_order": 0 if order == 1 else order - 1,
             "after_order": order, "change_type": chg_type},
            ev, method="inference", derivation=limitation))
        labels.append({"category_number": category[code], "code": code,
                       "supporting_change_ids": [cid]})
    r["classification"]["bond_changes"] = changes
    r["classification"]["labels"] = inferred(
        "formed-bond labels from partial source-derived model", labels, ev, limitation)
    r["classification"]["coverage"] = inferred(
        "partial source-derived classification", "partial", ev, limitation)
    r["classification"]["limitations"] = [limitation, "No complete atom-mapped reaction SMILES is asserted."]


def classify_cn(r, ev):
    classify(r, ev, [("C", "N", 1, "C-N")],
             "The named product and starting identities support C-N formation, but activation, leaving groups, proton transfers, and coproducts are unmapped.")


def classify_cc_ch(r, ev):
    classify(r, ev, [("C", "C", 1, "C-C"), ("C", "H", 1, "C-H")],
             "The substrate/product drawings support a new C-C connection and net C-H gain; allene atom correspondence and the proton source remain unresolved.")


def add_common_workup_molecules():
    specs = [
        ("m_diethyl_ether", "Et2O", "diethyl ether", [8, 9, 10], "CCOCC"),
        ("m_methanol", "MeOH", "methanol", [10], "CO"),
        ("m_triethylamine", "Et3N", "triethylamine", [8, 9, 10], "CCN(CC)CC"),
        ("m_nahso3", "NaHSO3", "sodium bisulfite", [10], None),
        ("m_mcpba", "mCPBA", "meta-chloroperoxybenzoic acid", [10], "OOC(=O)c1cccc(Cl)c1"),
    ]
    for mid, label, name, pages, smi in specs:
        add_molecule(mid, label, name, pages, smi, representation="smiles" if smi else "unresolved")


add_common_workup_molecules()


# Source molecules first appearing in pages 7--10.
new_specs = [
    ("m_1naphthoic_acid", "1-naphthoic acid", "1-naphthoic acid", [7], "O=C(O)c1cccc2ccccc12", None),
    ("m_picolinic_acid", "picolinic acid", "picolinic acid", [7], "O=C(O)c1ccccn1", None),
    ("m_thiophene2_carboxylic_acid", "thiophene-2-carboxylic acid", "thiophene-2-carboxylic acid", [7], "O=C(O)c1cccs1", None),
    ("m_acetic_acid", "acetic acid", "acetic acid", [7], "CC(=O)O", None),
    ("m_phenylacetic_acid", "2-phenylacetic acid", "phenylacetic acid", [7], "O=C(O)Cc1ccccc1", None),
    ("m_acrylic_acid", "acrylic acid", "acrylic acid", [7], "C=CC(=O)O", None),
    ("m_boc2o", "di-tert-butyl dicarbonate", "di-tert-butyl dicarbonate", [7], "CC(C)(C)OC(=O)OC(=O)OC(C)(C)C", None),
    ("m_cbz_glycine", "((benzyloxy)carbonyl)glycine", "N-Cbz-glycine", [7], "O=C(O)CNC(=O)OCc1ccccc1", None),
    ("m_pivalic_acid", "pivalic acid", "pivalic acid", [7], "CC(C)(C)C(=O)O", None),
    ("m_methylaminomethylpyridine", "2-[(methylamino)methyl]pyridine", "N-methyl-1-(pyridin-2-yl)methanamine", [8], "CNCc1ccccn1", None),
    ("m_r_benzyloxazolidinone", "(R)-4-benzyloxazolidin-2-one", "(R)-4-benzyloxazolidin-2-one", [8], "O=C1OC[C@@H](Cc2ccccc2)N1", "R configuration explicitly named by the source"),
    ("m_2bromomethylpyridine_hbr", "2-(bromomethyl)pyridine hydrobromide", "2-(bromomethyl)pyridine hydrobromide", [8], None, None),
    ("m_4picolylamine", "pyridin-4-ylmethanamine", "pyridin-4-ylmethanamine", [8], "NCc1ccncc1", None),
    ("m_benzylamine", "benzyl amine", "benzylamine", [8], "NCc1ccccc1", None),
    ("m_tosyl_chloride", "tosyl chloride", "4-methylbenzenesulfonyl chloride", [8, 9, 10], "Cc1ccc(S(=O)(=O)Cl)cc1", None),
    ("m_triphenylphosphine", "triphenylphosphine", "triphenylphosphine", [8, 9, 10], "P(c1ccccc1)(c1ccccc1)c1ccccc1", None),
    ("m_propargyl_alcohol", "propargyl alcohol", "prop-2-yn-1-ol", [8], "C#CCO", None),
    ("m_13", "13", "prop-2-yn-1-yl 4-methylbenzenesulfinate", [8, 9], "C#CCOS(=O)c1ccc(C)cc1", "Potential sulfur stereochemistry is not assigned by the source"),
    ("m_agsbf6", "AgSbF6", "silver hexafluoroantimonate(V)", [9, 10], None, None),
    ("m_liald4", "LiAlD4", "lithium aluminum deuteride", [9], "[Li+].[Al-]([2H])([2H])([2H])[2H]", None),
    ("m_ethyl_propiolate", "ethyl propiolate", "ethyl propiolate", [9], "CCOC(=O)C#C", None),
    ("m_14", "14", "prop-2-yn-1,1-d2-1-ol", [9], "[2H]C([2H])(O)C#C", None),
    ("m_15", "15", "prop-2-yn-1-yl-1,1-d2 4-methylbenzenesulfinate", [9], "C#CC([2H])([2H])OS(=O)c1ccc(C)cc1", "Potential sulfur stereochemistry is not assigned by the source"),
    ("m_16", "16", "ethyl 4-methylbenzenesulfinate", [9], "CCOS(=O)c1ccc(C)cc1", "Potential sulfur stereochemistry is not assigned by the source"),
    ("m_3butyn2ol", "3-butyn-2-ol", "but-3-yn-2-ol", [9], "C#CC(C)O", "The tetrahedral alcohol center is not assigned by the source"),
    ("m_17_mix", "17", "but-3-yn-2-yl 4-methylbenzenesulfinate, 1:1 diastereomer mixture", [9, 10], "C#CC(C)OS(=O)c1ccc(C)cc1", "Tetrahedral carbon and sulfur stereochemistry are unassigned"),
    ("m_but2yn1ol", "but-2-yn-1-ol", "but-2-yn-1-ol", [10], "CC#CCO", None),
    ("m_18", "18", "but-2-yn-1-yl 4-methylbenzenesulfinate", [10], "CC#CCOS(=O)c1ccc(C)cc1", "Potential sulfur stereochemistry is not assigned by the source"),
    ("m_19", "19", "p-tolyl 4-methylbenzenesulfinate", [10], "Cc1ccc(S(=O)Oc2ccc(C)cc2)cc1", "Potential sulfur stereochemistry is not assigned by the source"),
    ("m_thiophenol", "thiophenol", "benzenethiol", [10], "Sc1ccccc1", None),
]
for mid, label, name, pages, smi, stereo_note in new_specs:
    explicit_assignment = stereo_note if stereo_note and stereo_note.startswith("R configuration") else None
    add_molecule(mid, label, name, pages, smi,
                 chirality_note=None if explicit_assignment else stereo_note,
                 source_assignment=explicit_assignment,
                 representation="smiles" if smi else "unresolved")

# The isolated 17 material is a named mixture.  Its two source-unassigned
# diastereomer components remain distinct unresolved identities.
molecules["m_17_mix"]["stereochemistry"]["mixture"] = F(
    "1:1 mixture of diastereomers", {"ratio": [1, 1], "components": ["m_17_d1", "m_17_d2"]},
    tids(9, 10))
molecules["m_17_mix"]["stereochemistry"]["chirality_status"] = inferred(
    "1:1 mixture of diastereomers", "mixture", tids(9, 10),
    "The source explicitly characterizes 17 as a 1:1 diastereomer mixture without configurations.")
molecules["m_17_mix"]["stereochemistry"]["representation_status"] = inferred(
    "connectivity-only candidate for a diastereomer mixture", "partial", tids(9, 10),
    "The connectivity SMILES cannot encode the two source-unassigned diastereomers.")
add_molecule("m_17_d1", "17 diastereomer 1", "source-unassigned diastereomer 1 of 17",
             [9, 10], representation="unresolved",
             source_assignment="one component of the reported 1:1 diastereomer mixture; configuration unknown")
add_molecule("m_17_d2", "17 diastereomer 2", "source-unassigned diastereomer 2 of 17",
             [9, 10], representation="unresolved",
             source_assignment="second component of the reported 1:1 diastereomer mixture; configuration unknown")

# Page 8 gives the complete name for 1o and page 9 gives the exact isotope
# placement for 2a-d2, allowing their earlier unresolved registry entries to be
# upgraded without claiming source verification.
replace_existing_molecule(
    "m_1o", "(R)-4-benzyl-3-(pyridin-2-ylmethyl)oxazolidin-2-one", "1o", [8],
    "O=C1OC[C@@H](Cc2ccccc2)N1Cc1ccccn1",
    assignment="(R) configuration explicitly named; candidate atom index 4 is R by the standard diagnostic")
replace_existing_molecule(
    "m_2a_d2", "1-methyl-4-((propa-1,2-dien-1-yl-3,3-d2)sulfonyl)benzene", "2a-d2", [9],
    "[2H]C([2H])=C=CS(=O)(=O)c1ccc(C)cc1",
    axial_note="terminal CD2 end has two identical isotopic substituents; the allene axis is not stereogenic for this substitution pattern")

# Add source occurrences to existing products/reagents used in this page span.
for mid, pages in {
    "m_1d": (7, 11), "m_1e": (7, 11), "m_1f": (7, 11, 12), "m_1g": (7, 12),
    "m_1h": (7, 12), "m_1i": (7, 12), "m_1j": (7, 12), "m_1k": (7, 12),
    "m_1l": (7, 12), "m_1m": (7, 8), "m_1n": (8,), "m_1p": (8,), "m_11": (8,),
    "m_2a": (9, 10, 11, 12), "m_2b": (9, 10), "m_2c": (10,), "m_5a": (10,),
    "m_5b": (10,), "m_3a": (10, 11), "m_4a": (10, 11),
    "m_3b": (11,), "m_4b": (11,), "m_3c": (11,), "m_4c": (11,),
    "m_3d": (11,), "m_4d": (11,), "m_3e": (11,), "m_4e": (11,),
    "m_3f": (11, 12), "m_4f": (11, 12), "m_3g": (12,), "m_4g": (12,),
    "m_3h": (12,), "m_4h": (12,), "m_3i": (12,), "m_4i": (12,),
    "m_3j": (12,), "m_4j": (12,), "m_3k": (12,), "m_4k": (12,),
    "m_3l": (12,), "m_4l": (12,),
}.items():
    extend_molecule_evidence(mid, *pages)


def inherited_amidation_steps(r, specific_ev, procedure_letter, purification_text):
    gp_ev = tids(6)
    for e in gp_ev:
        if e not in r["evidence_ids"]:
            r["evidence_ids"].append(e)
    base_mid = "m_nahco3" if procedure_letter == "A" else "m_naoh"
    participant(r, base_mid, "workup_agent", specific_ev, stage="workup")
    participant(r, "m_water", "workup_agent", specific_ev, stage="workup")
    participant(r, "m_mgso4", "workup_agent", specific_ev, stage="workup")
    participant(r, "m_silica", "purification_agent", specific_ev, stage="purification")
    step(r, gp_ev, "charge acid, dichloromethane, and DMF; place in an ice/water bath",
         temp=("ice/water bath", "ice/water bath"), method="inherited_general_condition")
    step(r, gp_ev, "add oxalyl chloride dropwise by syringe",
         temp=("ice/water bath", "ice/water bath"), time=("over 5 min", Q(5, "min")),
         method="inherited_general_condition")
    step(r, gp_ev, "after 5 min remove the bath and warm to room temperature",
         temp=("room temperature", "r.t."), time=("5 min before warming", Q(5, "min")),
         method="inherited_general_condition")
    step(r, gp_ev, "maintain at room temperature and monitor by TLC and 1H NMR",
         temp=("room temperature", "r.t."),
         time=("2-4 h", {"minimum": 2, "maximum": 4, "unit": "h", "raw_text": "2-4 h"}),
         method="inherited_general_condition")
    step(r, gp_ev, "return to an ice/water bath and add the amine dropwise",
         temp=("ice/water bath", "ice/water bath"), method="inherited_general_condition")
    step(r, gp_ev, "remove bath, warm to room temperature, and stir vigorously",
         temp=("room temperature", "r.t."), time=("overnight", "overnight"),
         method="inherited_general_condition")
    work = ("add saturated aqueous sodium bicarbonate, separate, wash organic layer with water, dry over magnesium sulfate, filter, and concentrate"
            if procedure_letter == "A" else
            "add saturated aqueous sodium hydroxide, separate, extract aqueous layer twice with dichloromethane, combine organic layers, dry over magnesium sulfate, filter, and concentrate")
    step(r, gp_ev, work, stage="workup", method="inherited_general_condition")
    step(r, specific_ev, f"silica-gel flash chromatography ({purification_text})", stage="purification")


def add_amidation(tag, page, acid_mid, product_mid, letter, acid_mass, acid_amount,
                   ch2cl2, dmf_vol, dmf_amount, ox_vol, ox_amount, amine_vol,
                   amine_amount, product_mass, yld, purification, amine_mid="m_2picolylamine",
                   note=None):
    r, ev = base_reaction(
        f"synthesis_{tag}", [page], f"Compound synthesis paragraph for {tag}",
        f"synthesis of {tag}", "amide formation", unresolved=(
            "The experiment explicitly follows a general procedure; staged temperatures/times are inherited from page 6.",
        ) + ((note,) if note else ()))
    r["general_scheme_ids"] = [f"{DOC}:general_procedure_{letter.lower()}"]
    participant(r, acid_mid, "reactant", ev, mass=acid_mass, amount=acid_amount,
                eq=1.0, ref=acid_mid, method="original_extraction")
    participant(r, "m_ch2cl2", "solvent", ev, volume=ch2cl2)
    participant(r, "m_dmf", "reagent", ev, volume=dmf_vol, amount=dmf_amount,
                eq=0.05, ref=acid_mid, method="inherited_general_condition")
    participant(r, "m_oxalyl_chloride", "reagent", ev, volume=ox_vol, amount=ox_amount,
                eq=1.1, ref=acid_mid, method="inherited_general_condition")
    participant(r, amine_mid, "reactant", ev, volume=amine_vol, amount=amine_amount,
                eq=1.1, ref=acid_mid, method="inherited_general_condition")
    participant(r, product_mid, "product", ev, product_role="main")
    r["reaction_scale"] = F(acid_amount[0], {"reference_molecule_id": acid_mid,
                                              "amount": Q(acid_amount[1], acid_amount[2])}, ev)
    set_solvent(r, ev, ["m_ch2cl2"], "dichloromethane", ch2cl2,
                concentration=("0.5 M (general procedure)", 0.5, "M"),
                basis=f"relative to {acid_mid}", method="inherited_general_condition")
    inherited_amidation_steps(r, ev, letter, purification)
    outcome(r, ev, "product_amount", [product_mid], product_mass[0],
            Q(product_mass[1], product_mass[2]), role="isolated main product",
            basis="after silica-gel flash chromatography")
    outcome(r, ev, "yield", [product_mid], f"{yld}%", Q(yld, "%"),
            role="isolated main product", yield_type="isolated yield",
            basis="after silica-gel flash chromatography")
    classify_cn(r, ev)
    return r


# Complete the cross-page 1d record rather than duplicating the same experiment.
r = reactions[f"{DOC}:synthesis_1d"]
ev7 = tids(7)
merge_evidence(r, ev7)
step(r, ev7, "silica-gel flash chromatography (20-100% ethyl acetate/hexane)", stage="purification")
outcome(r, ev7, "product_amount", ["m_1d"], "445 mg", Q(445, "mg"),
        role="isolated main product", basis="after silica-gel flash chromatography")
outcome(r, ev7, "yield", ["m_1d"], "80%", Q(80, "%"), role="isolated main product",
        yield_type="isolated yield", basis="after silica-gel flash chromatography")
r["review"]["unresolved"] = [
    "The experiment follows General Procedure A; staged temperatures/times are inherited from page 6.",
    "Spectral data are characterization context and were not converted into outcomes.",
]

amidation_rows = [
    ("1e", 7, "m_1naphthoic_acid", "m_1e", "A", ("300 mg", 300, "mg"), ("1.7 mmol", 1.7, "mmol"), ("4 mL", 4, "mL"), ("0.01 mL", 0.01, "mL"), ("0.1 mmol", 0.1, "mmol"), ("0.17 mL", 0.17, "mL"), ("2.0 mmol", 2.0, "mmol"), ("0.20 mL", 0.20, "mL"), ("2.0 mmol", 2.0, "mmol"), ("353 mg", 353, "mg"), 77, "20-100% ethyl acetate/hexane"),
    ("1f", 7, "m_picolinic_acid", "m_1f", "A", ("300 mg", 300, "mg"), ("2.4 mmol", 2.4, "mmol"), ("5 mL", 5, "mL"), ("0.01 mL", 0.01, "mL"), ("0.12 mmol", 0.12, "mmol"), ("0.24 mL", 0.24, "mL"), ("2.8 mmol", 2.8, "mmol"), ("0.28 mL", 0.28, "mL"), ("2.8 mmol", 2.8, "mmol"), ("355 mg", 355, "mg"), 68, "20-100% ethyl acetate/hexane"),
    ("1g", 7, "m_thiophene2_carboxylic_acid", "m_1g", "A", ("300 mg", 300, "mg"), ("2.3 mmol", 2.3, "mmol"), ("5 mL", 5, "mL"), ("0.01 mL", 0.01, "mL"), ("0.1 mmol", 0.1, "mmol"), ("0.23 mL", 0.23, "mL"), ("2.6 mmol", 2.6, "mmol"), ("0.27 mL", 0.27, "mL"), ("2.6 mmol", 2.6, "mmol"), ("331 mg", 331, "mg"), 65, "0-100% ethyl acetate/hexane"),
    ("1h", 7, "m_acetic_acid", "m_1h", "B", ("250 mg", 250, "mg"), ("4.2 mmol", 4.2, "mmol"), ("8 mL", 8, "mL"), ("0.02 mL", 0.02, "mL"), ("0.2 mmol", 0.2, "mmol"), ("0.40 mL", 0.40, "mL"), ("4.7 mmol", 4.7, "mmol"), ("0.48 mL", 0.48, "mL"), ("4.7 mmol", 4.7, "mmol"), ("420 mg", 420, "mg"), 67, "0-100% acetone/ethyl acetate"),
    ("1i", 7, "m_phenylacetic_acid", "m_1i", "A", ("350 mg", 350, "mg"), ("2.6 mmol", 2.6, "mmol"), ("5 mL", 5, "mL"), ("0.01 mL", 0.01, "mL"), ("0.13 mmol", 0.13, "mmol"), ("0.25 mL", 0.25, "mL"), ("2.9 mmol", 2.9, "mmol"), ("0.30 mL", 0.30, "mL"), ("2.9 mmol", 2.9, "mmol"), ("484 mg", 484, "mg"), 83, "0-100% ethyl acetate/hexane"),
    ("1j", 7, "m_acrylic_acid", "m_1j", "B", ("315 mg", 315, "mg"), ("4.4 mmol", 4.4, "mmol"), ("9 mL", 9, "mL"), ("0.02 mL", 0.02, "mL"), ("0.2 mmol", 0.2, "mmol"), ("0.42 mL", 0.42, "mL"), ("4.9 mmol", 4.9, "mmol"), ("0.51 mL", 0.51, "mL"), ("4.9 mmol", 4.9, "mmol"), ("460 mg", 460, "mg"), 65, "0-100% acetone/ethyl acetate"),
    ("1l", 7, "m_cbz_glycine", "m_1l", "B", ("350 mg", 350, "mg"), ("1.7 mmol", 1.7, "mmol"), ("3 mL", 3, "mL"), ("0.01 mL", 0.01, "mL"), ("0.1 mmol", 0.1, "mmol"), ("0.16 mL", 0.16, "mL"), ("1.9 mmol", 1.9, "mmol"), ("0.19 mL", 0.19, "mL"), ("1.9 mmol", 1.9, "mmol"), ("127 mg", 127, "mg"), 25, "50-100% ethyl acetate/hexane, then 20-70% acetone/ethyl acetate"),
    ("1m", 7, "m_pivalic_acid", "m_1m", "B", ("300 mg", 300, "mg"), ("2.9 mmol", 2.9, "mmol"), ("6 mL", 6, "mL"), ("0.01 mL", 0.01, "mL"), ("0.1 mmol", 0.1, "mmol"), ("0.28 mL", 0.28, "mL"), ("3.3 mmol", 3.3, "mmol"), ("0.34 mL", 0.34, "mL"), ("3.3 mmol", 3.3, "mmol"), ("491 mg", 491, "mg"), 87, "0-100% ethyl acetate/hexane"),
]
for row in amidation_rows:
    add_amidation(*row)

# 1n follows General Procedure A with a different amine.  The acid amount is
# preserved exactly as printed even though 400 mg benzoic acid / 23.3 mmol is
# internally inconsistent.
add_amidation(
    "1n", 8, "m_benzoic_acid", "m_1n", "A",
    ("400 mg", 400, "mg"), ("23.3 mmol", 23.3, "mmol"), ("6 mL", 6, "mL"),
    ("0.013 mL", 0.013, "mL"), ("0.16 mmol", 0.16, "mmol"),
    ("0.31 mL", 0.31, "mL"), ("3.6 mmol", 3.6, "mmol"),
    ("0.44 mL", 0.44, "mL"), ("3.6 mmol", 3.6, "mmol"),
    ("573 mg", 573, "mg"), 77, "30-100% ethyl acetate/hexane",
    amine_mid="m_methylaminomethylpyridine",
    note="The printed benzoic-acid values are 400 mg and 23.3 mmol; the inconsistency is retained and not silently corrected.")


# Solvent-free Boc protection of 2-picolylamine (no solvent is stated).
r, ev = base_reaction("synthesis_1k", [7], "Compound synthesis paragraph for 1k",
                      "synthesis of tert-butyl (pyridin-2-ylmethyl)carbamate 1k",
                      "carbamate formation", unresolved=(
                          "No reaction solvent is reported; absence of a solvent charge is not converted into a solvent-free claim.",))
participant(r, "m_boc2o", "reactant", ev, mass=("250 mg", 250, "mg"), amount=("1.15 mmol", 1.15, "mmol"))
participant(r, "m_2picolylamine", "reactant", ev, volume=("0.12 mL", 0.12, "mL"), amount=("1.15 mmol", 1.15, "mmol"))
participant(r, "m_1k", "product", ev, product_role="main")
participant(r, "m_silica", "purification_agent", ev, stage="purification")
r["reaction_scale"] = F("1.15 mmol 2-picolylamine", {"reference_molecule_id": "m_2picolylamine", "amount": Q(1.15, "mmol")}, ev)
step(r, ev, "add di-tert-butyl dicarbonate and 2-picolylamine to a flame-dried flask under nitrogen",
     atmosphere="nitrogen")
step(r, ev, "stir", temp=("80 °C", Q(80, "°C")), time=("15 min", Q(15, "min")), atmosphere="nitrogen")
step(r, ev, "silica-gel flash chromatography (20-100% diethyl ether/hexane)", stage="purification")
outcome(r, ev, "product_amount", ["m_1k"], "230 mg", Q(230, "mg"), role="isolated main product", basis="after silica-gel flash chromatography")
outcome(r, ev, "yield", ["m_1k"], "96%", Q(96, "%"), role="isolated main product", yield_type="isolated yield", basis="after silica-gel flash chromatography")
classify_cn(r, ev)


add_molecule("m_nbutyllithium", "n-BuLi", "n-butyllithium", [8], "[Li]CCCC")

# Stereospecific N-alkylation furnishing named (R)-1o.
r, ev = base_reaction("synthesis_1o", [8], "Compound synthesis paragraph for 1o",
                      "synthesis of (R)-4-benzyl-3-(pyridin-2-ylmethyl)oxazolidin-2-one 1o",
                      "N-alkylation of an oxazolidinone", unresolved=(
                          "The source names the retained R configuration but reports no ee or optical purity.",))
participant(r, "m_r_benzyloxazolidinone", "reactant", ev, mass=("425 mg", 425, "mg"), amount=("2.4 mmol", 2.4, "mmol"))
participant(r, "m_thf", "solvent", ev, volume=("8 mL", 8, "mL"))
participant(r, "m_nbutyllithium", "reagent", ev, volume=("1.5 mL", 1.5, "mL"), amount=("2.4 mmol", 2.4, "mmol"))
participant(r, "m_2bromomethylpyridine_hbr", "reactant", ev, mass=("300 mg", 300, "mg"), amount=("1.2 mmol", 1.2, "mmol"))
participant(r, "m_1o", "product", ev, product_role="main; source names R configuration")
participant(r, "m_ch2cl2", "workup_agent", ev, stage="workup")
participant(r, "m_naoh", "workup_agent", ev, stage="workup")
participant(r, "m_mgso4", "workup_agent", ev, stage="workup")
participant(r, "m_silica", "purification_agent", ev, stage="purification")
r["reaction_scale"] = F("1.2 mmol 2-(bromomethyl)pyridine hydrobromide", {"reference_molecule_id": "m_2bromomethylpyridine_hbr", "amount": Q(1.2, "mmol")}, ev)
set_solvent(r, ev, ["m_thf"], "THF", ("8 mL", 8, "mL"))
step(r, ev, "charge (R)-4-benzyloxazolidin-2-one and THF; cool in a dry ice/acetone bath", temp=("-78 °C", Q(-78, "°C")), atmosphere="nitrogen")
step(r, ev, "add n-BuLi dropwise by syringe", temp=("-78 °C", Q(-78, "°C")), time=("over 5 min", Q(5, "min")), atmosphere="nitrogen")
step(r, ev, "maintain after n-BuLi addition", temp=("-78 °C", Q(-78, "°C")), time=("15 min", Q(15, "min")), atmosphere="nitrogen")
step(r, ev, "add 2-(bromomethyl)pyridine hydrobromide by temporarily removing the septum and stir", temp=("-78 °C", Q(-78, "°C")), time=("3 h", Q(3, "h")), atmosphere="nitrogen")
step(r, ev, "remove the dry ice/acetone bath, warm to room temperature, and stir", temp=("room temperature", "r.t."), time=("22 h", Q(22, "h")), atmosphere="nitrogen")
step(r, ev, "dilute with dichloromethane, wash with 10% sodium hydroxide, extract the aqueous layer with dichloromethane, combine, dry over magnesium sulfate, filter, and concentrate", stage="workup")
step(r, ev, "silica-gel flash chromatography (30-100% ethyl acetate/hexane)", stage="purification")
outcome(r, ev, "product_amount", ["m_1o"], "313 mg", Q(313, "mg"), role="isolated main product", basis="after silica-gel flash chromatography")
outcome(r, ev, "yield", ["m_1o"], "98%", Q(98, "%"), role="isolated main product", yield_type="isolated yield", basis="after silica-gel flash chromatography")
classify_cn(r, ev)


def add_modified_benzamide(tag, product_mid, amine_mid, amine_name, activation_h,
                           amine_vol, product_mass, yld, purification):
    r, ev = base_reaction(f"synthesis_{tag}", [8], f"Compound synthesis paragraph for {tag}",
                          f"synthesis of {tag}", "amide formation", unresolved=(
                              "The paragraph states a modified literature procedure and supplies its concrete staged operation.",))
    participant(r, "m_benzoic_acid", "reactant", ev, mass=("300 mg", 300, "mg"), amount=("2.5 mmol", 2.5, "mmol"))
    participant(r, "m_ch2cl2", "solvent", ev, volume=("5 mL", 5, "mL"))
    participant(r, "m_dmf", "reagent", ev, volume=("0.01 mL", 0.01, "mL"), amount=("0.12 mmol", 0.12, "mmol"))
    participant(r, "m_oxalyl_chloride", "reagent", ev, volume=("0.42 mL", 0.42, "mL"), amount=("4.9 mmol", 4.9, "mmol"))
    participant(r, amine_mid, "reactant", ev, volume=amine_vol, amount=("4.9 mmol", 4.9, "mmol"))
    participant(r, "m_triethylamine", "reagent", ev, volume=("0.68 mL", 0.68, "mL"), amount=("4.9 mmol", 4.9, "mmol"))
    participant(r, product_mid, "product", ev, product_role="main")
    participant(r, "m_nahco3", "workup_agent", ev, stage="workup")
    participant(r, "m_water", "workup_agent", ev, stage="workup")
    participant(r, "m_mgso4", "workup_agent", ev, stage="workup")
    participant(r, "m_silica", "purification_agent", ev, stage="purification")
    r["reaction_scale"] = F("2.5 mmol benzoic acid", {"reference_molecule_id": "m_benzoic_acid", "amount": Q(2.5, "mmol")}, ev)
    set_solvent(r, ev, ["m_ch2cl2"], "dichloromethane", ("5 mL initial charge; additional dichloromethane carries the amine/base but its volume is unreported", 5, "mL"))
    step(r, ev, "charge benzoic acid, dichloromethane, and DMF; cool in an ice/water bath", temp=("ice/water bath", "ice/water bath"), atmosphere="nitrogen")
    step(r, ev, "add oxalyl chloride dropwise", temp=("ice/water bath", "ice/water bath"), time=("30 min before warming", Q(30, "min")), atmosphere="nitrogen")
    step(r, ev, "remove bath, warm to room temperature, and maintain while monitoring by TLC", temp=("room temperature", "r.t."), time=(f"{activation_h} h", Q(activation_h, "h")), atmosphere="nitrogen")
    step(r, ev, "remove excess oxalyl chloride under high vacuum", pressure="high vacuum")
    step(r, ev, f"cool in an ice/water bath and add {amine_name} plus triethylamine in dichloromethane dropwise", temp=("ice/water bath", "ice/water bath"), atmosphere="nitrogen")
    step(r, ev, "remove bath, warm to room temperature, and stir vigorously", temp=("room temperature", "r.t."), time=("24 h", Q(24, "h")), atmosphere="nitrogen")
    step(r, ev, "add saturated aqueous sodium bicarbonate, dilute with dichloromethane, separate, wash organic layer with water, dry over magnesium sulfate, filter, and concentrate", stage="workup")
    step(r, ev, f"silica-gel flash chromatography ({purification})", stage="purification")
    outcome(r, ev, "product_amount", [product_mid], product_mass[0], Q(product_mass[1], product_mass[2]), role="isolated main product", basis="after silica-gel flash chromatography")
    outcome(r, ev, "yield", [product_mid], f"{yld}%", Q(yld, "%"), role="isolated main product", yield_type="isolated yield", basis="after silica-gel flash chromatography")
    classify_cn(r, ev)


add_modified_benzamide("1p", "m_1p", "m_4picolylamine", "pyridin-4-ylmethanamine",
                       5, ("0.50 mL", 0.50, "mL"), ("300 mg", 300, "mg"), 58,
                       "20-100% ethyl acetate/hexane, then 5% methanol/ethyl acetate")
add_modified_benzamide("11", "m_11", "m_benzylamine", "benzylamine",
                       8, ("0.54 mL", 0.54, "mL"), ("429 mg", 429, "mg"), 83,
                       "10-40% ethyl acetate/hexane")


def add_sulfinate_ester(tag, pages, alcohol_mid, alcohol_charge, product_mid,
                        product_mass, yld, addition_minutes, reaction_hours,
                        bath_range, dcm_second_volume, purification,
                        product_note=None):
    """Concrete PPh3-mediated preparation of a propargylic sulfinate."""
    r, ev = base_reaction(
        f"synthesis_{tag}", pages, f"Compound synthesis paragraph for {tag}",
        f"synthesis of sulfinate {tag}", "formation of a sulfinic ester",
        unresolved=((product_note,) if product_note else ()))
    participant(r, "m_tosyl_chloride", "reactant", ev,
                mass=("3.00 g", 3.00, "g") if tag == "13" else ("1.00 g", 1.00, "g") if tag in ("17", "18") else ("985 mg", 985, "mg"),
                amount=("15.7 mmol", 15.7, "mmol") if tag == "13" else ("5.3 mmol", 5.3, "mmol") if tag in ("17", "18") else ("5.2 mmol", 5.2, "mmol"))
    initial_dcm = 39 if tag == "13" else 13
    participant(r, "m_ch2cl2", "solvent", ev, volume=(f"{initial_dcm} mL initial charge", initial_dcm, "mL"))
    et3n = ("2.4 mL", 2.4, "mL", "17.3 mmol", 17.3) if tag == "13" else (
        ("0.82 mL", 0.82, "mL", "5.9 mmol", 5.9) if tag == "15" else
        ("0.8 mL", 0.8, "mL", "5.8 mmol", 5.8))
    participant(r, "m_triethylamine", "reagent", ev,
                volume=(et3n[0], et3n[1], et3n[2]), amount=(et3n[3], et3n[4], "mmol"))
    pph3_mass = ("4.12 g", 4.12, "g") if tag == "13" else (
        ("1355 mg", 1355, "mg") if tag == "15" else ("1.38 g", 1.38, "g"))
    pph3_amt = 15.7 if tag == "13" else 5.2 if tag == "15" else 5.3
    participant(r, "m_triphenylphosphine", "reagent", ev, mass=pph3_mass,
                amount=(f"{pph3_amt} mmol", pph3_amt, "mmol"))
    participant(r, "m_ch2cl2", "solvent", ev,
                volume=(f"{dcm_second_volume} mL in the alcohol/PPh3 solution", dcm_second_volume, "mL"))
    participant(r, alcohol_mid, "reactant", ev, **alcohol_charge)
    participant(r, product_mid, "product", ev, product_role="main")
    if tag == "15":
        participant(r, "m_16", "product", ev, product_role="co-isolated side product")
    participant(r, "m_diethyl_ether", "workup_agent", ev, stage="workup")
    participant(r, "m_hexane", "workup_agent", ev, stage="workup")
    participant(r, "m_silica", "purification_agent", ev, stage="purification")
    set_solvent(r, ev, ["m_ch2cl2"],
                f"dichloromethane, {initial_dcm} mL initial plus {dcm_second_volume} mL in the added solution",
                (f"{initial_dcm + dcm_second_volume} mL total stated DCM", initial_dcm + dcm_second_volume, "mL"))
    scale_amt = 15.7 if tag == "13" else 3.44 if tag == "15" else 5.3
    scale_mid = "m_tosyl_chloride" if tag == "13" else alcohol_mid
    r["reaction_scale"] = F(
        f"{scale_amt} mmol {'tosyl chloride' if tag == '13' else molecules[alcohol_mid]['names']['english']['normalized_value']}",
        {"reference_molecule_id": scale_mid, "amount": Q(scale_amt, "mmol")}, ev)
    step(r, ev, "charge tosyl chloride and dichloromethane, add triethylamine, and place in an ice/water bath",
         temp=(bath_range, {"minimum": int(bath_range.split('-')[0]), "maximum": int(bath_range.split('-')[1].split()[0]), "unit": "°C", "raw_text": bath_range}),
         atmosphere="nitrogen")
    step(r, ev, "prepare a separate solution of triphenylphosphine, dichloromethane, and the alcohol substrate")
    step(r, ev, "add the alcohol/triphenylphosphine solution dropwise while maintaining the bath range",
         temp=(bath_range, {"minimum": int(bath_range.split('-')[0]), "maximum": int(bath_range.split('-')[1].split()[0]), "unit": "°C", "raw_text": bath_range}),
         time=(f"over {addition_minutes} min", Q(addition_minutes, "min")) if addition_minutes is not None else None,
         atmosphere="nitrogen")
    step(r, ev, "continue the reaction while monitoring by TLC" + (" and 1H NMR" if tag == "15" else ""),
         temp=(bath_range, {"minimum": int(bath_range.split('-')[0]), "maximum": int(bath_range.split('-')[1].split()[0]), "unit": "°C", "raw_text": bath_range}),
         time=(f"{reaction_hours} h", Q(reaction_hours, "h")), atmosphere="nitrogen")
    step(r, ev, "transfer to a beaker, add 20% diethyl ether/hexane to precipitate solids, filter through a silica pad under reduced pressure, and rinse with ether", stage="workup")
    step(r, ev, "collect the filtrate and concentrate under reduced pressure", stage="workup", pressure="reduced pressure")
    step(r, ev, f"silica-gel flash chromatography ({purification})", stage="purification")
    outcome(r, ev, "product_amount", [product_mid], product_mass[0], Q(product_mass[1], product_mass[2]),
            role="isolated/co-isolated main material", basis="after silica-gel flash chromatography")
    if yld is not None:
        outcome(r, ev, "yield", [product_mid], f"{yld}%", Q(yld, "%"),
                role="isolated main product", yield_type="isolated yield",
                basis="after silica-gel flash chromatography")
    if tag == "15":
        outcome(r, ev, "selectivity", [product_mid, "m_16"],
                "title compound 15 : side product 16 = 85:15",
                {"ratio": [85, 15], "order": [product_mid, "m_16"], "metric": "composition ratio"},
                basis="1H NMR characterization of the isolated mixture",
                scope="co-isolated product composition; configurations not assigned")
    classify(r, ev, [("O", "S", 1, "O-S")],
             "The named sulfinic ester supports new O-S connectivity, but reagent-derived sulfur/oxygen correspondence and coproducts are not fully mapped.")


add_sulfinate_ester(
    "13", [8, 9], "m_propargyl_alcohol",
    {"volume": ("0.91 mL", 0.91, "mL"), "amount": ("15.7 mmol", 15.7, "mmol")},
    "m_13", ("2.35 g", 2.35, "g"), 77, 10, 3, "10-15 °C", 39,
    "0-15% ethyl acetate/hexane")
add_sulfinate_ester(
    "15", [9], "m_14",
    {"mass": ("approximately 200 mg of 14 in a 14/Et2O/EtOH 1:0.8:0.2 mixture", 200, "mg"),
     "amount": ("3.44 mmol 14", 3.44, "mmol")},
    "m_15", ("701 mg mixture containing 15 and 16", 701, "mg"), None, None, 2.5,
    "10-20 °C", 9, "0-20% ethyl acetate/hexane",
    "The isolated 701 mg material contains 15 and ethyl sulfinate 16 in an 85:15 ratio; no isolated yield for 15 is stated.")
add_sulfinate_ester(
    "17", [9, 10], "m_3butyn2ol",
    {"volume": ("0.4 mL", 0.4, "mL"), "amount": ("5.3 mmol", 5.3, "mmol")},
    "m_17_mix", ("0.94 g", 0.94, "g"), 86, 15, 3, "10-15 °C", 13,
    "0-20% ethyl acetate/hexane",
    "The product is explicitly a 1:1 mixture of source-unassigned diastereomers.")
add_sulfinate_ester(
    "18", [10], "m_but2yn1ol",
    {"volume": ("0.39 mL", 0.39, "mL"), "amount": ("5.3 mmol", 5.3, "mmol")},
    "m_18", ("0.80 g", 0.80, "g"), 73, 20, 5.5, "10-15 °C", 13,
    "0-20% ethyl acetate/hexane")


def add_silver_rearrangement(tag, page, precursor_mid, precursor_mass, precursor_amount,
                             ag_mass, ag_amount, dcm_volume, add_minutes, hours,
                             product_mid, product_mass, yld, purification):
    r, ev = base_reaction(
        f"synthesis_{tag}", [page], f"Compound synthesis paragraph for {tag}",
        f"silver-promoted rearrangement to allenyl sulfone {tag}",
        "rearrangement of a propargylic sulfinate to an allenyl sulfone",
        unresolved=("The source calls this a modified literature procedure; bond classification is a partial source-derived net-change model.",))
    participant(r, "m_agsbf6", "reagent", ev, mass=ag_mass, amount=ag_amount)
    participant(r, precursor_mid, "reactant", ev, mass=precursor_mass, amount=precursor_amount)
    participant(r, "m_ch2cl2", "solvent", ev, volume=(f"{dcm_volume} mL", dcm_volume, "mL"))
    participant(r, product_mid, "product", ev, product_role="main")
    participant(r, "m_silica", "purification_agent", ev, stage="purification")
    participant(r, "m_diethyl_ether", "purification_agent", ev, stage="purification")
    r["reaction_scale"] = F(precursor_amount[0],
        {"reference_molecule_id": precursor_mid, "amount": Q(precursor_amount[1], precursor_amount[2])}, ev)
    set_solvent(r, ev, ["m_ch2cl2"], "dichloromethane", (f"{dcm_volume} mL", dcm_volume, "mL"))
    step(r, ev, "charge silver hexafluoroantimonate in a flame-dried flask in a glovebox; equip the flask and place under nitrogen", atmosphere="nitrogen")
    step(r, ev, "add the propargylic sulfinate solution in dichloromethane dropwise by syringe/addition funnel",
         time=(f"over {add_minutes} min", Q(add_minutes, "min")), atmosphere="nitrogen")
    step(r, ev, "continue under nitrogen while monitoring by TLC",
         time=(f"{hours} h", Q(hours, "h")), atmosphere="nitrogen")
    step(r, ev, "vacuum-filter through a silica pad using diethyl ether, collect the filtrate, and concentrate under reduced pressure", stage="workup", pressure="reduced pressure")
    step(r, ev, f"silica-gel flash chromatography ({purification})", stage="purification")
    outcome(r, ev, "product_amount", [product_mid], product_mass[0], Q(product_mass[1], product_mass[2]),
            role="isolated main product", basis="after silica-gel flash chromatography")
    outcome(r, ev, "yield", [product_mid], f"{yld}%", Q(yld, "%"),
            role="isolated main product", yield_type="isolated yield",
            basis="after silica-gel flash chromatography")
    classify(r, ev, [("C", "S", 1, "C-S")],
             "The propargylic-sulfinate/allenyl-sulfone structures support a rearranged C-S bond, but a complete atom mapping and bond-order accounting are not asserted.")
    return r, ev


add_silver_rearrangement("2a", 9, "m_13", ("2.35 g", 2.35, "g"),
                         ("12.1 mmol", 12.1, "mmol"), ("83 mg", 83, "mg"),
                         ("0.24 mmol", 0.24, "mmol"), 24, 10, 2, "m_2a",
                         ("2.22 g", 2.22, "g"), 94, "0-30% ethyl acetate/hexane")
add_silver_rearrangement("2a_d2", 9, "m_15", ("100.9 mg", 100.9, "mg"),
                         ("0.51 mmol", 0.51, "mmol"), ("3.5 mg", 3.5, "mg"),
                         ("0.01 mmol", 0.01, "mmol"), 1, 5, 3, "m_2a_d2",
                         ("71 mg", 71, "mg"), 70, "5-30% ethyl acetate/hexane")
add_silver_rearrangement("2b", 10, "m_17_mix", ("0.93 g", 0.93, "g"),
                         ("4.5 mmol", 4.5, "mmol"), ("31 mg", 31, "mg"),
                         ("0.09 mmol", 0.09, "mmol"), 9, 5, 4, "m_2b",
                         ("0.80 g", 0.80, "g"), 86, "0-30% ethyl acetate/hexane")
r_2c, ev_2c = add_silver_rearrangement(
    "2c", 10, "m_18", ("0.79 g", 0.79, "g"), ("3.8 mmol", 3.8, "mmol"),
    ("26 mg", 26, "mg"), ("0.08 mmol", 0.08, "mmol"), 8, 5, 4, "m_2c",
    ("0.63 g", 0.63, "g"), 80, "0-30% ethyl acetate/hexane")
participant(r_2c, "m_19", "product", ev_2c, product_role="chromatographically separated side product")
outcome(r_2c, ev_2c, "product_amount", ["m_19"], "64 mg", Q(64, "mg"),
        role="isolated side product", basis="chromatographic separation")
outcome(r_2c, ev_2c, "yield", ["m_19"], "7%", Q(7, "%"),
        role="isolated side product", yield_type="isolated yield",
        basis="chromatographic separation")
r_2c["review"]["unresolved"].append(
    "Formation of side product 19 is preserved, but its distinct atom mapping is not represented by the partial main-product classification.")


# Selective deuteride reduction of ethyl propiolate to labeled propargyl alcohol 14.
r, ev = base_reaction(
    "synthesis_14", [9], "Compound synthesis paragraph for 14",
    "synthesis of prop-2-yn-1,1-d2-1-ol 14", "deuteride reduction of an ester",
    unresolved=(
        "Product 14 was split into two solvent-containing vial fractions; the source estimates the amount of 14 in each rather than reporting a single purified mass.",
        "The 25% yield is an overall yield for the estimated combined amount of 14."))
participant(r, "m_liald4", "reagent", ev, mass=("642 mg", 642, "mg"), amount=("15.3 mmol", 15.3, "mmol"))
participant(r, "m_diethyl_ether", "solvent", ev, volume=("26 mL initial", 26, "mL"))
participant(r, "m_ethyl_propiolate", "reactant", ev, volume=("2.07 mL", 2.07, "mL"), amount=("20.4 mmol", 20.4, "mmol"))
participant(r, "m_diethyl_ether", "solvent", ev, volume=("13 mL with ethyl propiolate", 13, "mL"))
participant(r, "m_14", "product", ev, product_role="main, isolated as solvent-containing fractions")
participant(r, "m_water", "workup_agent", ev, volume=("0.7 mL first addition", 0.7, "mL"), stage="workup")
participant(r, "m_naoh", "workup_agent", ev, mass=("0.7 g NaOH in 4.7 mL water", 0.7, "g"), volume=("4.7 mL water in 15% solution", 4.7, "mL"), stage="workup")
participant(r, "m_water", "workup_agent", ev, volume=("2 mL final addition", 2, "mL"), stage="workup")
participant(r, "m_mgso4", "workup_agent", ev, stage="workup")
r["reaction_scale"] = F("20.4 mmol ethyl propiolate", {"reference_molecule_id": "m_ethyl_propiolate", "amount": Q(20.4, "mmol")}, ev)
set_solvent(r, ev, ["m_diethyl_ether"], "diethyl ether, 26 mL initial plus 13 mL substrate solution", ("39 mL total stated Et2O", 39, "mL"))
step(r, ev, "charge lithium aluminum deuteride and diethyl ether; cool in an ethanol bath maintained by cryo-cool",
     temp=("-70 °C", Q(-70, "°C")))
step(r, ev, "add ethyl propiolate in diethyl ether by syringe pump at 25 mL/h while maintaining the bath between -70 and -60 °C",
     temp=("-70 to -60 °C", {"minimum": -70, "maximum": -60, "unit": "°C", "raw_text": "-70 to -60 °C"}),
     time=("over 30 min", Q(30, "min")))
step(r, ev, "raise the bath temperature after addition",
     temp=("-45 °C", Q(-45, "°C")))
step(r, ev, "maintain at -45 °C while monitoring by TLC until complete consumption",
     temp=("-45 °C", Q(-45, "°C")), time=("overnight", "overnight"))
step(r, ev, "add water, then 15% sodium hydroxide solution, then water to form an off-white precipitate; gravity-filter and rinse with diethyl ether", stage="workup")
step(r, ev, "dry the filtrate over magnesium sulfate, gravity-filter, and distill at atmospheric pressure to remove diethyl ether", stage="workup", pressure="atmospheric pressure")
step(r, ev, "transfer the bulk residue to vial 1; rinse the flask into vial 2 and further concentrate vial 2 by rotary evaporation",
     temp=("10 °C bath for vial 2 concentration", Q(10, "°C")), stage="workup", pressure="reduced pressure")
outcome(r, ev, "product_amount", ["m_14"],
        "vial 1: 445 mg mixture, 14/Et2O/EtOH = 1:0.8:0.2, approximately 200 mg 14; vial 2: 114 mg mixture, 1:0.11:0.08, approximately 95 mg 14",
        {"fractions": [
            {"mixture_mass": Q(445, "mg"), "composition_ratio": [1, 0.8, 0.2], "order": ["m_14", "m_diethyl_ether", "m_ethanol"], "estimated_product_mass": Q(200, "mg")},
            {"mixture_mass": Q(114, "mg"), "composition_ratio": [1, 0.11, 0.08], "order": ["m_14", "m_diethyl_ether", "m_ethanol"], "estimated_product_mass": Q(95, "mg")}]},
        role="two solvent-containing product fractions", basis="1H NMR composition and estimated product mass")
outcome(r, ev, "yield", ["m_14"], "overall yield 25%", Q(25, "%"),
        role="combined estimated product", yield_type="overall yield",
        basis="estimated combined amount of 14 in the two reported vial fractions")
classify(r, ev, [("C", "H", 1, "C-H"), ("O", "H", 1, "O-H")],
         "The reduction/deuteration creates C-D bonds (classified under C-H taxonomy) and an alcohol O-H bond; ester bond-order changes and atom/proton sources are not completely mapped.")


# Two-stage preparation and oxidation of crossover disulfone 5b.
r, ev = base_reaction(
    "synthesis_5b", [10], "Compound synthesis paragraph for 5b",
    "synthesis of 1-methyl-4-((2-(phenylsulfonyl)allyl)sulfonyl)benzene 5b",
    "thiol addition to allenyl sulfone followed by oxidation",
    unresolved=(
        "The source reports a single 71% isolated yield from 2a to 5b across the addition/oxidation sequence; it is not assigned as a yield of the oxidation stage alone.",
        "The intermediate sulfide is not isolated or structurally registered as a source product."))
participant(r, "m_2a", "reactant", ev, mass=("100 mg", 100, "mg"), amount=("0.52 mmol", 0.52, "mmol"))
participant(r, "m_thiophenol", "reactant", ev, volume=("0.06 mL", 0.06, "mL"), amount=("0.6 mmol", 0.6, "mmol"))
participant(r, "m_triethylamine", "reagent", ev, volume=("0.01 mL", 0.01, "mL"), amount=("0.05 mmol", 0.05, "mmol"))
participant(r, "m_methanol", "solvent", ev, volume=("2 mL", 2, "mL"))
participant(r, "m_ch2cl2", "solvent", ev, volume=("5 mL for oxidation", 5, "mL"))
participant(r, "m_mcpba", "reagent", ev, mass=("300 mg", 300, "mg"), amount=("1.6 mmol", 1.6, "mmol"))
participant(r, "m_5b", "product", ev, product_role="main two-stage product")
participant(r, "m_nahso3", "workup_agent", ev, stage="workup")
participant(r, "m_nahco3", "workup_agent", ev, stage="workup")
participant(r, "m_mgso4", "workup_agent", ev, stage="workup")
participant(r, "m_silica", "purification_agent", ev, stage="purification")
r["reaction_scale"] = F("0.52 mmol allenyl sulfone 2a", {"reference_molecule_id": "m_2a", "amount": Q(0.52, "mmol")}, ev)
set_solvent(r, ev, ["m_methanol", "m_ch2cl2"], "2 mL methanol for addition; residue dissolved in 5 mL dichloromethane for oxidation", ("7 mL total stated solvent across stages", 7, "mL"))
step(r, ev, "charge allenyl sulfone 2a, thiophenol, triethylamine, and methanol and react",
     temp=("room temperature", "r.t."), time=("30 min", Q(30, "min")))
step(r, ev, "concentrate the first-stage reaction under reduced pressure", stage="workup", pressure="reduced pressure")
step(r, ev, "dissolve the residue in dichloromethane and add meta-chloroperoxybenzoic acid in one portion",
     temp=("room temperature", "r.t."))
step(r, ev, "continue oxidation; TLC indicated the reaction was not complete at the stated endpoint",
     temp=("room temperature", "r.t."), time=("1 h", Q(1, "h")))
step(r, ev, "dilute with dichloromethane, wash with NaHSO3 solution and saturated NaHCO3 solution, dry over magnesium sulfate, filter, and concentrate", stage="workup")
step(r, ev, "silica-gel flash chromatography (0-50% ethyl acetate/hexane)", stage="purification")
outcome(r, ev, "product_amount", ["m_5b"], "124 mg", Q(124, "mg"), role="isolated main product", basis="after silica-gel flash chromatography")
outcome(r, ev, "yield", ["m_5b"], "71%", Q(71, "%"), role="isolated main product", yield_type="overall isolated yield", basis="overall sequence from charged 2a after silica-gel flash chromatography")
classify(r, ev, [("C", "S", 1, "C-S"), ("S", "O", 2, "DB-O-S")],
         "The final disulfone structure supports C-S formation and oxidation to S=O connectivity, but the unisolated intermediate and oxygen atom correspondence are not completely mapped.")


def set_participant_charge(p, key, raw, number, unit, ev, method="original_extraction"):
    p["charges"][key] = F(raw, Q(number, unit), ev, method=method)


def replace_missing_fact(container, key, raw, normalized, ev, method="original_extraction"):
    container[key] = F(raw, normalized, ev, method=method)


# Page 10/11 supplies the concrete preparative paragraph for Table 1 entry 7.
entry7 = reactions[f"{DOC}:table1_entry_07"]
ev7 = tids(10, 11)
for p in entry7["participants"]:
    if p["molecule_id"] == "m_1a":
        set_participant_charge(p, "mass", "40 mg", 40, "mg", ev7)
        set_participant_charge(p, "amount", "0.19 mmol", 0.19, "mmol", ev7)
    elif p["molecule_id"] == "m_2a":
        set_participant_charge(p, "mass", "74 mg", 74, "mg", ev7)
        set_participant_charge(p, "amount", "0.38 mmol", 0.38, "mmol", ev7)
    elif p["molecule_id"] == "m_thf":
        set_participant_charge(p, "volume", "1 mL", 1, "mL", ev7)
entry7["outcomes"][0]["yield_type"] = F("isolated after silica-gel flash chromatography", "isolated yield", ev7)
entry7["outcomes"][0]["measurement_basis"] = F(
    "51.7 mg isolated after silica-gel flash chromatography; 67% combined 3a/4a",
    "isolated after silica-gel flash chromatography", ev7)
entry7["outcomes"][1]["measurement_basis"] = F(
    "product selectivity measured after chromatography using 1H NMR integration ratios",
    "post-chromatography 1H NMR integration", ev7)
entry7["notes"].append(F(
    "The concrete paragraph reports a red-brown solid, 51.7 mg, 67%, 3a/4a = 98:2.",
    {"product_mass": Q(51.7, "mg"), "yield": Q(67, "%"), "ratio": [98, 2]}, ev7))
entry7["review"]["unresolved"] = [
    x for x in entry7["review"]["unresolved"] if "measurement basis" not in x]
entry7["review"]["unresolved"].append(
    "The isolated material is a 98:2 3a/4a mixture; product-specific isolated masses are not reported.")
merge_evidence(entry7, ev7)


# Page 11's 4a characterization paragraph explicitly says THF, whereas Table 1
# entry 6 on page 2 says EtOH.  Both source values are preserved as a conflict.
entry6 = reactions[f"{DOC}:table1_entry_06"]
ev6 = tids(11)
for p in entry6["participants"]:
    if p["molecule_id"] == "m_1a":
        set_participant_charge(p, "mass", "40 mg", 40, "mg", ev6)
        set_participant_charge(p, "amount", "0.19 mmol", 0.19, "mmol", ev6)
    elif p["molecule_id"] == "m_2a":
        set_participant_charge(p, "mass", "74 mg", 74, "mg", ev6)
        set_participant_charge(p, "amount", "0.38 mmol", 0.38, "mmol", ev6)
participant(entry6, "m_thf", "solvent", ev6, volume=("1 mL in the 4a characterization paragraph", 1, "mL"))
entry6["solvent_conditions"]["composition"] = F(
    "Table 1 entry 6: EtOH; 4a characterization paragraph: THF",
    None, tids(2, 11), status="conflict")
entry6["solvent_conditions"]["composition"]["normalized_value"] = None
entry6["solvent_conditions"]["composition"]["alternatives"] = [
    F("EtOH (Table 1 entry 6)", ["m_ethanol"], tids(2)),
    F("THF (page 11 characterization paragraph)", ["m_thf"], tids(11))]
entry6["solvent_conditions"]["total_volume"] = F(
    "1 mL in Table 1 general footnote and 1 mL in the characterization paragraph",
    None, tids(2, 11), status="conflict")
entry6["solvent_conditions"]["total_volume"]["normalized_value"] = None
entry6["solvent_conditions"]["total_volume"]["alternatives"] = [
    F("1 mL EtOH", Q(1, "mL"), tids(2)),
    F("1 mL THF", Q(1, "mL"), tids(11))]
entry6["outcomes"][0]["yield_type"] = F("isolated after silica-gel flash chromatography", "isolated yield", ev6)
entry6["outcomes"][0]["measurement_basis"] = F(
    "46.2 mg isolated after silica-gel flash chromatography; 53% combined 3a/4a",
    "isolated after silica-gel flash chromatography", ev6)
entry6["outcomes"][1]["measurement_basis"] = F(
    "product selectivity measured after chromatography using 1H NMR integration ratios",
    "post-chromatography 1H NMR integration", ev6)
entry6["notes"].append(F(
    "Source conflict retained: Table 1 entry 6 assigns EtOH, while the page-11 4a paragraph says THF for the same 53%, 63:37 experiment.",
    {"table_1_solvent": "m_ethanol", "characterization_paragraph_solvent": "m_thf"}, tids(2, 11)))
entry6["notes"].append(F(
    "The characterization paragraph reports a red-brown solid, 46.2 mg, 53%, 3a/4a = 63:37.",
    {"product_mass": Q(46.2, "mg"), "yield": Q(53, "%"), "ratio": [63, 37]}, ev6))
entry6["review"]["unresolved"] = [
    x for x in entry6["review"]["unresolved"] if "measurement basis" not in x]
entry6["review"]["unresolved"].extend([
    "Unresolved source conflict: Table 1 entry 6 reports EtOH, but the 4a characterization paragraph reports THF.",
    "The isolated material is a 63:37 3a/4a mixture; product-specific isolated masses are not reported."])
merge_evidence(entry6, ev6)


def inherited_scope_steps(r, current_ev, procedure, duration_h, purification,
                          extra_purification=None):
    # Page 3 supplies the Figure 1 default/footnoted scope duration; page 6
    # supplies the explicitly invoked General Procedure C/D operations.
    proc_ev = list(dict.fromkeys(current_ev + tids(3, 6)))
    for e in tids(3, 6):
        if e not in r["evidence_ids"]:
            r["evidence_ids"].append(e)
    opening = ("charge liquid picolyl amide, THF, and allenyl sulfone 2a open to air and seal the screw-top tube"
               if procedure == "D" else
               "charge solid picolyl amide and allenyl sulfone 2a in air and seal the screw-top tube")
    step(r, proc_ev, opening, method="inherited_general_condition")
    step(r, proc_ev, "evacuate for 1-2 s and fill with nitrogen; repeat three times",
         time=("1-2 s evacuation per cycle", {"minimum": 1, "maximum": 2, "unit": "s", "raw_text": "1-2 s"}),
         atmosphere="nitrogen", method="inherited_general_condition")
    if procedure == "C":
        step(r, proc_ev, "add THF by syringe, wrap the cap with parafilm, and lower the tube into a preheated oil bath",
             temp=("50 °C", Q(50, "°C")), atmosphere="nitrogen",
             method="inherited_general_condition")
    else:
        step(r, proc_ev, "wrap the cap with parafilm and lower the tube into a preheated oil bath",
             temp=("50 °C", Q(50, "°C")), atmosphere="nitrogen",
             method="inherited_general_condition")
    step(r, proc_ev, "stir and monitor reaction progress by 1H NMR",
         temp=("50 °C", Q(50, "°C")), time=(f"{duration_h} h", Q(duration_h, "h")),
         atmosphere="nitrogen", method="inherited_general_condition")
    step(r, proc_ev, "dilute with dichloromethane, transfer, and concentrate under reduced pressure",
         stage="workup", pressure="reduced pressure", method="inherited_general_condition")
    step(r, current_ev, f"silica-gel flash chromatography ({purification})", stage="purification")
    if extra_purification:
        step(r, current_ev, extra_purification, stage="purification")


scope_run_data = {
    "b": {"page": [11], "procedure": "C", "sub_mass": 45.3, "sub_mmol": 0.2, "allene_mass": 78, "allene_mmol": 0.4, "thf": 1.0,
          "purification": "20-70% ethyl acetate/hexane", "duration": 48,
          "runs": [(1, 55.0, 65, (96, 4)), (2, 52.6, 63, (96, 4))]},
    "c": {"page": [11], "procedure": "C", "sub_mass": 46.1, "sub_mmol": 0.2, "allene_mass": 78, "allene_mmol": 0.4, "thf": 1.0,
          "purification": "20-60% ethyl acetate/hexane", "duration": 48,
          "runs": [(1, 54.5, 64, (95, 5)), (2, 54.1, 64, (97, 3))]},
    "d": {"page": [11], "procedure": "D", "sub_mass": 49.0, "sub_mmol": 0.2, "allene_mass": 78, "allene_mmol": 0.4, "thf": 1.0,
          "purification": "20-70% ethyl acetate/hexane", "duration": 48,
          "runs": [(1, 49.6, 56, (96, 4), 7.3, 66), (2, 45.9, 52, (97, 3), 10.1, 65)]},
    "e": {"page": [11], "procedure": "C", "sub_mass": 52.5, "sub_mmol": 0.2, "allene_mass": 78, "allene_mmol": 0.4, "thf": 1.0,
          "purification": "20-70% ethyl acetate/hexane", "duration": 48,
          "runs": [(1, 52.8, 58, (94, 6)), (2, 56.7, 62, (95, 5))]},
    "f": {"page": [11, 12], "procedure": "C", "sub_mass": 42.6, "sub_mmol": 0.2, "allene_mass": 78, "allene_mmol": 0.4, "thf": 1.0,
          "purification": "20-70% ethyl acetate/hexane", "duration": 48,
          "runs": [(1, 54.0, 66, (99, 1)), (2, 49.7, 61, (98, 2))]},
    "g": {"page": [12], "procedure": "C", "sub_mass": 43.7, "sub_mmol": 0.2, "allene_mass": 78, "allene_mmol": 0.4, "thf": 1.0,
          "purification": "20-70% ethyl acetate/hexane", "duration": 48,
          "runs": [(1, 47.7, 58, (94, 6)), (2, 46.1, 56, (94, 6))]},
    "h": {"page": [12], "procedure": "C", "sub_mass": 30.0, "sub_mmol": 0.2, "allene_mass": 78, "allene_mmol": 0.4, "thf": 1.0,
          "purification": "50-100% ethyl acetate/hexane, then 10-40% acetone/ethyl acetate", "duration": 48,
          "runs": [(1, 41.4, 60, (94, 6)), (2, 39.1, 57, (94, 6))]},
    "i": {"page": [12], "procedure": "C", "sub_mass": 45.3, "sub_mmol": 0.2, "allene_mass": 78, "allene_mmol": 0.4, "thf": 1.0,
          "purification": "20-70% ethyl acetate/hexane", "duration": 48,
          "runs": [(1, 52.5, 62, (92, 8)), (2, 50.0, 59, (92, 8))]},
    "j": {"page": [12], "procedure": "C", "sub_mass": 32.5, "sub_mmol": 0.2, "allene_mass": 78, "allene_mmol": 0.4, "thf": 1.0,
          "purification": "20-100% ethyl acetate/hexane", "duration": 48,
          "runs": [(1, 49.9, 70, (94, 6)), (2, 49.3, 69, (92, 8))]},
    "k": {"page": [12], "procedure": "D", "sub_mass": 41.7, "sub_mmol": 0.2, "allene_mass": 78, "allene_mmol": 0.4, "thf": 1.0,
          "purification": "20-70% diethyl ether/hexane", "duration": None,
          "runs": [(1, 44.3, 55, (96, 4), 48), (2, 47.5, 59, (96, 4), 72)]},
    "l": {"page": [12], "procedure": "C", "purification": "20-100% ethyl acetate/hexane", "duration": 24,
          "runs": [
              (1, 44.9, 0.15, 58.3, 0.30, 0.8, 44.2, 60, (85, 15)),
              (2, 44.9, 0.15, 58.3, 0.30, 0.8, 41.4, 56, (94, 6)),
              (3, 10.0, 0.033, 13.0, 0.067, 0.18, 10.6, 64, (95, 5)),
              (4, 7.0, 0.024, 9.0, 0.047, 0.15, 7.0, 58, (96, 4))]},
}


def add_scope_run(suffix, run_tuple, data):
    pages = data["page"]
    ev = tids(*pages)
    if suffix == "l":
        run_no, sub_mass, sub_mmol, allene_mass, allene_mmol, thf_vol, prod_mass, yld, ratio = run_tuple
        duration = data["duration"]
    else:
        run_no, prod_mass, yld, ratio, *extras = run_tuple
        sub_mass, sub_mmol = data["sub_mass"], data["sub_mmol"]
        allene_mass, allene_mmol, thf_vol = data["allene_mass"], data["allene_mmol"], data["thf"]
        duration = extras[0] if suffix == "k" else data["duration"]
    unresolved = [
        "The isolated material contains major product 3 and minor product 4; product-specific masses are not reported.",
        "The newly generated tetrahedral center is not assigned; no racemic composition is inferred."]
    if suffix == "k":
        unresolved.append("The isolated material was contaminated with starting material 1k; the reported yield was calculated from 1H NMR molar ratios.")
    r, current_ev = base_reaction(
        f"scope_3{suffix}_run{run_no}", pages,
        f"Compound paragraph for 3{suffix}, run {run_no}",
        f"individual reported run {run_no}: formation of 3{suffix}/4{suffix}",
        "C(sp3)-H functionalization with allenyl sulfone 2a",
        categories=("substrate_scope",), unresolved=tuple(unresolved))
    r["general_scheme_ids"] = [f"{DOC}:figure1_generic"]
    participant(r, f"m_1{suffix}", "reactant", current_ev,
                mass=(f"{sub_mass:g} mg", sub_mass, "mg"),
                amount=(f"{sub_mmol:g} mmol", sub_mmol, "mmol"), eq=1.0, ref=f"m_1{suffix}")
    participant(r, "m_2a", "reactant", current_ev,
                mass=(f"{allene_mass:g} mg", allene_mass, "mg"),
                amount=(f"{allene_mmol:g} mmol", allene_mmol, "mmol"), eq=2.0, ref=f"m_1{suffix}")
    participant(r, "m_thf", "solvent", current_ev, volume=(f"{thf_vol:g} mL", thf_vol, "mL"))
    participant(r, f"m_3{suffix}", "product", current_ev, product_role="major product in isolated mixture")
    participant(r, f"m_4{suffix}", "product", current_ev, product_role="minor product in isolated mixture")
    participant(r, "m_ch2cl2", "workup_agent", current_ev, stage="workup")
    participant(r, "m_silica", "purification_agent", current_ev, stage="purification")
    if suffix == "h":
        participant(r, "m_acetone", "purification_agent", current_ev, stage="purification")
    r["reaction_scale"] = F(f"{sub_mmol:g} mmol of 1{suffix}",
        {"reference_molecule_id": f"m_1{suffix}", "amount": Q(sub_mmol, "mmol")}, current_ev)
    concentration = sub_mmol / thf_vol
    set_solvent(r, current_ev, ["m_thf"], "THF", (f"{thf_vol:g} mL", thf_vol, "mL"),
                (f"{concentration:.3g} M from stated charge and volume", concentration, "M"),
                basis=f"{sub_mmol:g} mmol 1{suffix} in {thf_vol:g} mL THF")
    r["solvent_conditions"]["concentration"] = inferred(
        f"{concentration:.3g} M from {sub_mmol:g} mmol in {thf_vol:g} mL THF",
        Q(concentration, "M"), current_ev,
        "Calculated from the source-reported substrate amount and THF volume.")
    r["solvent_conditions"]["concentration_basis"] = inferred(
        f"{sub_mmol:g} mmol 1{suffix} in {thf_vol:g} mL THF",
        f"{sub_mmol:g} mmol 1{suffix} in {thf_vol:g} mL THF", current_ev,
        "This basis records the source values used for the calculated concentration.")
    extra_purification = None
    if suffix == "k":
        extra_purification = "for characterization only, perform a second silica-gel flash chromatography (20-70% ethyl acetate/hexane)"
    inherited_scope_steps(r, current_ev, data["procedure"], duration,
                          data["purification"], extra_purification)
    outcome(r, current_ev, "product_amount", [f"m_3{suffix}", f"m_4{suffix}"],
            f"{prod_mass:g} mg isolated product mixture", Q(prod_mass, "mg"),
            role=f"combined 3{suffix}/4{suffix} mixture",
            basis="isolated after silica-gel flash chromatography")
    if suffix == "k":
        ybasis = "yield calculated using molar ratios from 1H NMR integration because isolated material contained 1k"
        ytype = "NMR-calculated yield of chromatographed mixture"
    else:
        ybasis = "isolated after silica-gel flash chromatography"
        ytype = "isolated yield"
    outcome(r, current_ev, "yield", [f"m_3{suffix}", f"m_4{suffix}"],
            f"{yld}%", Q(yld, "%"), role=f"combined 3{suffix}/4{suffix} mixture",
            yield_type=ytype, basis=ybasis)
    outcome(r, current_ev, "selectivity", [f"m_3{suffix}", f"m_4{suffix}"],
            f"3{suffix}/4{suffix} = {ratio[0]}:{ratio[1]}",
            {"ratio": list(ratio), "order": [f"m_3{suffix}", f"m_4{suffix}"], "metric": "product ratio"},
            basis="post-chromatography 1H NMR integration ratios",
            scope=f"3{suffix}/4{suffix} product selectivity; stereochemical configuration remains unassigned")
    if suffix == "d":
        recovered_mass, brsm = run_tuple[4], run_tuple[5]
        outcome(r, current_ev, "product_amount", ["m_1d"],
                f"{recovered_mass:g} mg recovered 1d", Q(recovered_mass, "mg"),
                role="recovered starting material", basis="recovered after chromatography")
        outcome(r, current_ev, "other", ["m_3d", "m_4d"],
                f"{brsm}% yield based on recovered starting material", Q(brsm, "%"),
                role="combined 3d/4d mixture", basis="yield based on recovered starting material",
                scope="BRSM")
    if suffix == "k":
        r["notes"].append(F(
            "The source states that the chromatographed product was contaminated with 1k; yield was calculated from 1H NMR molar ratios.",
            "product contaminated with 1k; NMR-calculated yield", current_ev))
    classify_cc_ch(r, current_ev)


for suffix, data in scope_run_data.items():
    for run_tuple in data["runs"]:
        add_scope_run(suffix, run_tuple, data)


# Characterization-only passages and named-compound spectra remain source
# context.  They are deliberately not promoted to synthesis experiments.
bundle["source_context"].extend([
    {
        "context_id": f"{DOC}:context_pages7_8_precursor_characterization",
        "document_id": DOC,
        "evidence_ids": tids(7, 8),
        "kind": "characterization_or_analysis",
        "description": "NMR, HRMS, TLC, physical-state, and literature-match statements for 1d-1p and 11 characterize compounds associated with preparative records; spectra are not separate reactions.",
    },
    {
        "context_id": f"{DOC}:context_pages8_10_sulfinate_allene_characterization",
        "document_id": DOC,
        "evidence_ids": tids(8, 9, 10),
        "kind": "characterization_or_analysis",
        "description": "NMR, HRMS, TLC, and literature-match data for 13-19, 2a, 2a-d2, 2b, 2c, and 5a/5b are characterization context. Compound 16, compound 19, and 5a have characterization paragraphs on these pages but no separate preparation beyond the linked mixture/byproduct records.",
    },
    {
        "context_id": f"{DOC}:context_pages10_12_product_characterization",
        "document_id": DOC,
        "evidence_ids": tids(10, 11, 12),
        "kind": "characterization_or_analysis",
        "description": "NMR, IR, HRMS, melting-point, and TLC data for 3a-3l and distinguishable minor-product resonances for 4a-4l support identity/selectivity but are not standalone syntheses or mass-spectrum reactions.",
    },
])


coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 7,
            "state": "attempted",
            "regions": [
                {"region": "top-left continuation of 1d paragraph", "finding": "The cross-page preparation of 1d was completed with its purification, 445 mg isolated mass, and 80% isolated yield; spectral data remain characterization context.", "evidence_ids": tids(7)},
                {"region": "two-column synthesis paragraphs for 1e-1j", "finding": "Six concrete amide preparations were recorded separately with page-specific charges and the ordered operations inherited from explicitly cited General Procedures A or B.", "evidence_ids": tids(7)},
                {"region": "lower-half paragraphs for 1k-1m", "finding": "The Boc preparation of 1k and amide preparations of 1l and 1m were recorded separately; the 1m spectra continue on page 8 and do not create another reaction.", "evidence_ids": tids(7)},
            ],
            "record_ids": [
                f"{DOC}:synthesis_1d", f"{DOC}:synthesis_1e", f"{DOC}:synthesis_1f",
                f"{DOC}:synthesis_1g", f"{DOC}:synthesis_1h", f"{DOC}:synthesis_1i",
                f"{DOC}:synthesis_1j", f"{DOC}:synthesis_1k", f"{DOC}:synthesis_1l",
                f"{DOC}:synthesis_1m"],
            "unresolved": [
                "The concrete amidation paragraphs explicitly invoke General Procedures A/B; their staged times and temperatures are inherited from page 6 and remain machine-interpreted.",
                "No stereochemical configuration is inferred for molecules whose source drawings/names do not assign one."],
        },
        {
            "page": 8,
            "state": "attempted",
            "regions": [
                {"region": "top continuation of 1m characterization", "finding": "The remaining 1m spectra are characterization context for the preparation on page 7, not a second synthesis.", "evidence_ids": tids(8)},
                {"region": "synthesis paragraphs for 1n, (R)-1o, 1p, and 11", "finding": "Four concrete preparations were recorded. The explicit R assignment for 1o was preserved; the printed 1n benzoic-acid values 400 mg/23.3 mmol were retained as a source conflict rather than corrected.", "evidence_ids": tids(8)},
                {"region": "bottom paragraph beginning synthesis of 13", "finding": "The preparation of sulfinate 13 begins here and continues on page 9; it is one merged cross-page experiment.", "evidence_ids": tids(8)},
            ],
            "record_ids": [
                f"{DOC}:synthesis_1m", f"{DOC}:synthesis_1n", f"{DOC}:synthesis_1o",
                f"{DOC}:synthesis_1p", f"{DOC}:synthesis_11", f"{DOC}:synthesis_13"],
            "unresolved": [
                "The source does not resolve the 400 mg/23.3 mmol inconsistency in the 1n paragraph.",
                "No ee or optical-purity result is reported for named (R)-1o."],
        },
        {
            "page": 9,
            "state": "attempted",
            "regions": [
                {"region": "top continuation of 13 and preparation of 2a", "finding": "The cross-page isolation of 13 and the AgSbF6-mediated preparation of 2a were recorded as distinct experiments with ordered addition and purification steps.", "evidence_ids": tids(9)},
                {"region": "middle-column preparations of 14, 15/16, and 2a-d2", "finding": "The low-temperature deuteride reduction to 14, co-isolation of 15/16 at 85:15, and rearrangement to 2a-d2 were recorded separately. Fraction compositions for 14 and the absence of an isolated yield for 15 were preserved.", "evidence_ids": tids(9)},
                {"region": "bottom paragraph beginning preparation of 17", "finding": "The preparation of 17 begins here and continues on page 10; the stated 1:1 diastereomer mixture is retained with unresolved configurations.", "evidence_ids": tids(9)},
            ],
            "record_ids": [
                f"{DOC}:synthesis_13", f"{DOC}:synthesis_2a", f"{DOC}:synthesis_14",
                f"{DOC}:synthesis_15", f"{DOC}:synthesis_2a_d2", f"{DOC}:synthesis_17"],
            "unresolved": [
                "The two diastereomers of 17 are not assigned configurations and cannot be represented as unique source-verified isomeric SMILES.",
                "The isolated 701 mg for 15 is a mixture with 16; no product-specific mass or yield is reported."],
        },
        {
            "page": 10,
            "state": "attempted",
            "regions": [
                {"region": "top continuation of 17 and preparations of 2b, 18, and 2c/19", "finding": "The isolation of 17, rearrangements to 2b and 2c, preparation of 18, and chromatographically separated 19 side product were recorded. Compound 19's following spectra are characterization only.", "evidence_ids": tids(10)},
                {"region": "middle characterization of 5a and preparation of 5b", "finding": "The 5a paragraph contains only prior byproduct characterization and was not duplicated; the two-stage thiol-addition/oxidation preparation of 5b was recorded with its overall 71% yield.", "evidence_ids": tids(10)},
                {"region": "bottom paragraph beginning 3a Table 1 entry 7", "finding": "The preparative paragraph continuing on page 11 was merged into the existing Table 1 entry 7 record rather than duplicated.", "evidence_ids": tids(10)},
            ],
            "record_ids": [
                f"{DOC}:synthesis_17", f"{DOC}:synthesis_2b", f"{DOC}:synthesis_18",
                f"{DOC}:synthesis_2c", f"{DOC}:synthesis_5b", f"{DOC}:table1_entry_07"],
            "unresolved": [
                "The 7% side-product 19 branch is preserved but not completely atom-mapped.",
                "The intermediate sulfide in the two-stage 5b preparation is not isolated or fully specified."],
        },
        {
            "page": 11,
            "state": "attempted",
            "regions": [
                {"region": "top continuation of 3a and 4a characterization paragraph", "finding": "Table 1 entries 7 and 6 were enriched with concrete masses, isolated-yield basis, and post-chromatography NMR selectivity. The entry-6 EtOH/THF source conflict was retained explicitly.", "evidence_ids": tids(11)},
                {"region": "individual run paragraphs for 3b-3e", "finding": "Both reported runs for each substrate were recorded separately, including 3d recovered-starting-material masses and per-run BRSM values.", "evidence_ids": tids(11)},
                {"region": "bottom paragraph beginning 3f", "finding": "The two 3f runs continue on page 12 and were merged across the page boundary; spectral peak tables remain characterization context.", "evidence_ids": tids(11)},
            ],
            "record_ids": [
                f"{DOC}:table1_entry_07", f"{DOC}:table1_entry_06",
                *[f"{DOC}:scope_3{x}_run{n}" for x in "bcde" for n in (1, 2)],
                f"{DOC}:scope_3f_run1", f"{DOC}:scope_3f_run2"],
            "unresolved": [
                "Unresolved source conflict: Table 1 entry 6 names EtOH, while the page-11 4a paragraph names THF for the same 53%, 63:37 experiment.",
                "All isolated scope materials are major/minor product mixtures; product-specific masses and absolute configurations are unreported."],
        },
        {
            "page": 12,
            "state": "attempted",
            "regions": [
                {"region": "top continuation of 3f and individual run paragraphs for 3g-3j", "finding": "The two runs for each substrate were recorded separately with isolated masses, yields, and post-chromatography 1H NMR product ratios.", "evidence_ids": tids(12)},
                {"region": "3k paragraph", "finding": "The 48 h and 72 h runs were separated; contamination by 1k and NMR-calculated yield basis were retained, together with the second purification used only for characterization.", "evidence_ids": tids(12)},
                {"region": "3l paragraph and four run statements", "finding": "All four differently scaled runs were recorded independently, including the 85:15 run that the Figure 1 average excluded; spectral data continuing to page 13 were not imported into this coverage pass.", "evidence_ids": tids(12)},
            ],
            "record_ids": [
                f"{DOC}:scope_3f_run1", f"{DOC}:scope_3f_run2",
                *[f"{DOC}:scope_3{x}_run{n}" for x in "ghijk" for n in (1, 2)],
                *[f"{DOC}:scope_3l_run{n}" for n in (1, 2, 3, 4)]],
            "unresolved": [
                "Product-specific masses for the 3/4 mixtures and configurations at the newly formed tetrahedral centers are unreported.",
                "The 3l characterization paragraph continues onto page 13, outside this pass; no page-13 facts are claimed as reviewed here."],
        },
    ],
}


page_inventory = {
    "document_id": DOC,
    "pages": list(range(7, 13)),
    "source_title": "Metal-Free C-C Coupling of an Allenyl Sulfone with Picolyl Amides to Access Vinyl Sulfones via Pyridine-Initiated In Situ Generation of Sulfinate Anion",
    "presentation_conversion": False,
    "notes": [
        "All six original PDF page images were inspected at original detail in this session.",
        "Characterization-only spectra and mass peaks were retained as source context rather than synthesis records.",
        "Supporting Information was not associated because none was explicitly supplied in the workspace mapping."],
    "new_reaction_ids": [r["reaction_id"] for r in bundle["reactions"] if r["reaction_id"] not in set(reactions.keys())],
}
# The dictionary above was updated as records were added, so derive new IDs by
# page evidence and an explicit tag set rather than relying on dictionary age.
new_prefixes = ("synthesis_1e", "synthesis_1f", "synthesis_1g", "synthesis_1h",
                "synthesis_1i", "synthesis_1j", "synthesis_1k", "synthesis_1l",
                "synthesis_1m", "synthesis_1n", "synthesis_1o", "synthesis_1p",
                "synthesis_11", "synthesis_13", "synthesis_14", "synthesis_15",
                "synthesis_17", "synthesis_18", "synthesis_2a", "synthesis_2a_d2",
                "synthesis_2b", "synthesis_2c", "synthesis_5b", "scope_")
page_inventory["new_reaction_ids"] = [
    r["reaction_id"] for r in bundle["reactions"]
    if r["reaction_id"].split(":", 1)[1].startswith(new_prefixes)]

stereo_checks = {}
for m in bundle["molecules"]:
    if any(e.startswith(f"{DOC}:page{p}:") for p in range(7, 13)
           for e in m["occurrences"][0]["evidence_ids"]):
        smi = m["smiles"]["normalized_value"]
        if smi:
            stereo_checks[m["molecule_id"]] = stereo_inspect(smi)

dataset = json.loads((RUN / "dataset.json").read_text())
mentions = source_mentions(
    dataset,
    ["1d", "1e", "1f", "1g", "1h", "1i", "1j", "1k", "1l", "1m",
     "1n", "1o", "1p", "11", "13", "14", "15", "16", "17", "18", "19",
     "2a", "2b", "2c", "3a", "4a", "3b", "3l", "5a", "5b"],
    context=220)

(RAW / "page_inventory_0007_0012.json").write_text(
    json.dumps(page_inventory, ensure_ascii=False, indent=2) + "\n")
(RAW / "source_mentions_pages_0007_0012.json").write_text(
    json.dumps(mentions, ensure_ascii=False, indent=2) + "\n")
(RAW / "stereo_checks_pages_0007_0012.json").write_text(
    json.dumps(stereo_checks, ensure_ascii=False, indent=2) + "\n")
(RUN / "coverage/pages_0007_0012.json").write_text(
    json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")
OUT.write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({
    "molecules": len(bundle["molecules"]),
    "reactions": len(bundle["reactions"]),
    "new_reactions": len(page_inventory["new_reaction_ids"]),
    "coverage_pages": [p["page"] for p in coverage["pages"]]}, ensure_ascii=False))
