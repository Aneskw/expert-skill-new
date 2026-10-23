#!/usr/bin/env python3
"""Build the source-derived v3.1 bundle for PDF pages 1–6.

This script is a reproducibility artifact. It uses only the prepared dataset and
the page evidence already present in this workspace. It does not certify the
machine structures or chemistry.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_1d33e0721b132f96"
RUN = ROOT / "runs/current"
DATASET = json.loads((RUN / "dataset.json").read_text())
EVT = {n: f"{DOC}:page{n}:text" for n in range(1, 13)}
EVI = {n: f"{DOC}:page{n}:image" for n in range(1, 7)}


def ev(*pages, image=True):
    out = []
    for page in pages:
        out.append(EVT[page])
        if image and page in EVI:
            out.append(EVI[page])
    return list(dict.fromkeys(out))


def src(raw, normalized=None, evidence=()):
    return st.fact(raw, raw if normalized is None else normalized, evidence)


def inf(normalized, evidence, derivation):
    return st.fact(None, normalized, evidence, method="inference", derivation=derivation)


def derived(normalized, evidence, derivation, raw="source structure drawing"):
    return st.fact(raw, normalized, evidence, method="structure_derivation", derivation=derivation)


def missing(evidence=()):
    return st.fact(evidence=evidence)


def conflict(alternatives, evidence):
    f = st.fact(None, None, evidence)
    f["data_status"] = "conflict"
    f["alternatives"] = alternatives
    return f


bundle = st.empty(DATASET["documents"])
molecules = {}
diagnostics = {}


def molecule(mid, label, name, evidence, smiles=None, stereo="achiral", note=None):
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = src(label, label, evidence)
    m["names"]["original"] = src(name, name, evidence)
    m["names"]["english"] = src(name, name, evidence)
    m["cas_numbers"] = missing(evidence)
    if smiles:
        diag = st.stereo_inspect(smiles)
        diagnostics[mid] = diag
        m["representation_type"] = "smiles"
        m["smiles"] = derived(
            diag["canonical_isomeric_smiles"], evidence,
            "Connectivity candidate transcribed from the source drawing/name and canonicalized by RDKit; it remains machine-inferred and not source-verified.",
            raw=f"structure drawing/name for {label}",
        )
        m["stereochemistry"]["isomeric_smiles"] = derived(
            diag["canonical_isomeric_smiles"], evidence,
            "Canonical isomeric SMILES from the same machine connectivity candidate; unassigned source stereochemistry remains unassigned.",
            raw=f"structure drawing/name for {label}",
        )
        m["stereochemistry"]["identity_key"] = inf(
            diag["stereo_graph_key"], evidence,
            "SHA-256 stereo graph key calculated from the canonical isomeric SMILES; not a sample identity or source verification.",
        )
        rep = "encoded" if stereo == "achiral" else "partial"
    else:
        m["representation_type"] = "unresolved"
        rep = "not_encoded"
    if stereo == "allene_unassigned":
        m["stereochemistry"]["axial"] = src(
            "2,3-allenoic-acid allene axis shown without R_a/S_a or P/M assignment",
            [{"type": "allene axis", "configuration": "unassigned"}], evidence,
        )
        chirality = "unspecified"
        geometry = "not_applicable"
    elif stereo == "center_unassigned":
        m["stereochemistry"]["atom_centers"] = src(
            "tetrahedral product center shown without wedge/dash or R/S assignment",
            [{"configuration": "unassigned", "source": "product drawing"}], evidence,
        )
        chirality = "unspecified"
        geometry = "not_applicable"
    elif stereo == "p_and_center_unassigned":
        m["stereochemistry"]["atom_centers"] = src(
            "tetrahedral carbon and potentially stereogenic phosphorus are not assigned",
            [{"configuration": "unassigned", "elements": ["C", "P"]}], evidence,
        )
        chirality = "unspecified"
        geometry = "not_applicable"
    elif stereo == "unresolved":
        chirality = "unresolved"
        geometry = "unresolved"
    else:
        chirality = "achiral"
        geometry = "not_applicable"
        m["stereochemistry"]["assignment"] = inf(
            "no stereogenic center/axis identified in the machine connectivity candidate",
            evidence,
            "Graph inspection of the machine connectivity candidate; absence of source stereochemical labels alone was not used.",
        )
    m["stereochemistry"]["chirality_status"] = inf(chirality, evidence, "Status assigned from the inspected source drawing plus the stated structure-candidate limitations.")
    m["stereochemistry"]["geometry_status"] = inf(geometry, evidence, "No independently assignable E/Z element is asserted; allene axial stereochemistry is tracked separately.")
    m["stereochemistry"]["representation_status"] = inf(rep, evidence, "Machine representation status relative to the stereochemical information visible in the source.")
    m["identity_note"] = src(note, note, evidence) if note else src(
        "Machine structure candidate; parsing does not establish source-level structural correctness." if smiles else
        "Source label/drawing retained; complete machine SMILES withheld because the structure was not encoded with sufficient confidence.",
        evidence=evidence,
    )
    molecules[mid] = m
    return mid


def generic_molecule(mid, label, fixed, variables, evidence):
    molecule(mid, label, label, evidence, stereo="unresolved")
    m = molecules[mid]
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": src(fixed, fixed, evidence),
        "variables": [
            {
                "label": v,
                "attachment_points": src(f"source-labelled attachment for {v}", [v], evidence),
                "definition": src(f"{v} is variable as drawn", f"source variable {v}", evidence),
                "allowed_values": missing(evidence),
                "constraints": missing(evidence),
            } for v in variables
        ],
        "instance_mapping": [],
        "reason": src("R-group structure in Scheme 1", "Generic source drawing; not a concrete molecule.", evidence),
    }
    return mid


# Generic Scheme 1 identities.
generic_molecule(f"{DOC}:mol:generic_allenoic_acid", "generic 2,3-allenoic acid 1", "2,3-allenoic acid scaffold", ["R1", "R2", "R3"], ev(1, 2))
generic_molecule(f"{DOC}:mol:generic_p_reagent", "generic H–P(O) reagent 2", "H–P(O)(R4)(R5)", ["R4", "R5"], ev(1, 2))
generic_molecule(f"{DOC}:mol:generic_phosphorylated_butenolide", "generic product 3/4", "phosphorylated furan-2(5H)-one scaffold", ["R1", "R2", "R3", "R4", "R5"], ev(1, 2))
generic_molecule(f"{DOC}:mol:generic_halobutenolide", "Scheme 1a halobutenolide", "beta-halobutenolide scaffold", ["R1", "R2", "R3", "X"], ev(2))
generic_molecule(f"{DOC}:mol:generic_cc_butenolide", "Scheme 1b C–C-functionalized butenolide", "C-substituted butenolide scaffold", ["R1", "R2", "R3", "FG"], ev(2))
generic_molecule(f"{DOC}:mol:generic_cs_butenolide", "Scheme 1c C–S-functionalized butenolide", "S-substituted butenolide scaffold", ["R1", "R2", "R3", "R"], ev(2))

# Scope substrates/products. These are evidence-linked candidates, not source-verified structures.
aryl = {
    "a": "c1ccccc1", "b": "c1ccc(C)cc1", "c": "c1ccc(C(C)(C)C)cc1",
    "d": "c1ccc(OC)cc1", "e": "c1ccc(F)cc1", "f": "c1ccc(Cl)cc1",
    "g": "c1ccc(Br)cc1", "h": "c1ccc(C(F)(F)F)cc1", "i": "c1cc(Cl)ccc1",
    "j": "c1ccccc1C", "k": "c1ccccc1OC", "l": "c1ccccc1Cl",
    "m": "c1c(C)ccc(C)c1", "n": "c1cc(C(F)(F)F)cc(C(F)(F)F)c1",
    "o": "c1ccc2ccccc2c1", "p": "c1ccccc1",
}
scope3_yields = {"a": 80, "b": 68, "c": 42, "d": 70, "e": 48, "f": 82, "g": 78, "h": 62, "i": 77, "j": 37, "k": 70, "l": 66, "m": 77, "n": 52, "o": 51, "p": 72, "q": 28, "r": 33}
for label in scope3_yields:
    evidence = ev(3, 4) if label in ("a", "b") else ev(4, 6) if label in ("c", "d") else ev(4)
    # Keep branch ring-closure digits distinct from the still-open lactone ring.
    # Reusing digit 1 here is syntactically valid but incorrectly creates a fused ring.
    product_aryl = aryl.get(label, "").translate(str.maketrans({"1": "3", "2": "4"}))
    if label == "p":
        s1 = f"CCCC(=C=C{aryl[label]})C(=O)O"
        prod = f"CCCC1=C(P(=O)(c2ccccc2)c2ccccc2)C({product_aryl})OC1=O"
    elif label == "q":
        s1 = "CC(=C=C(C)c1ccccc1)C(=O)O"
        prod = "CC1=C(P(=O)(c2ccccc2)c2ccccc2)C(C)(c2ccccc2)OC1=O"
    elif label == "r":
        s1 = "CC(=C=C(c1ccccc1)c1ccccc1)C(=O)O"
        prod = "CC1=C(P(=O)(c2ccccc2)c2ccccc2)C(c2ccccc2)(c2ccccc2)OC1=O"
    else:
        s1 = f"CC(=C=C{aryl[label]})C(=O)O"
        prod = f"CC1=C(P(=O)(c2ccccc2)c2ccccc2)C({product_aryl})OC1=O"
    molecule(f"{DOC}:mol:1{label}", f"1{label}", f"2,3-allenoic acid 1{label}", evidence, s1, "allene_unassigned")
    product_name = {
        "a": "4-(diphenylphosphoryl)-3-methyl-5-phenylfuran-2(5H)-one",
        "b": "4-(diphenylphosphoryl)-3-methyl-5-(p-tolyl)furan-2(5H)-one",
        "c": "5-(4-(tert-butyl)phenyl)-4-(diphenylphosphoryl)-3-methylfuran-2(5H)-one",
        "d": "4-(diphenylphosphoryl)-5-(4-methoxyphenyl)-3-methylfuran-2(5H)-one",
    }.get(label, f"phosphorylated butenolide 3{label}")
    stereo = "center_unassigned" if label != "r" else "achiral"
    molecule(f"{DOC}:mol:3{label}", f"3{label}", product_name, evidence, prod, stereo)

# Phosphorus reagents and Scheme 3 products.
p_groups = {
    "a": ("2a", "diphenylphosphine oxide", "O=[PH](c1ccccc1)c1ccccc1", "P(=O)(c2ccccc2)c2ccccc2"),
    "b": ("2b", "bis(4-methylphenyl)phosphine oxide", "O=[PH](c1ccc(C)cc1)c1ccc(C)cc1", "P(=O)(c2ccc(C)cc2)c2ccc(C)cc2"),
    "c": ("2c", "bis(4-methoxyphenyl)phosphine oxide", "O=[PH](c1ccc(OC)cc1)c1ccc(OC)cc1", "P(=O)(c2ccc(OC)cc2)c2ccc(OC)cc2"),
    "d": ("2d", "bis([1,1'-biphenyl]-4-yl)phosphine oxide", "O=[PH](c1ccc(-c2ccccc2)cc1)c1ccc(-c2ccccc2)cc1", "P(=O)(c2ccc(-c3ccccc3)cc2)c2ccc(-c3ccccc3)cc2"),
    "e": ("2e", "bis(4-fluorophenyl)phosphine oxide", "O=[PH](c1ccc(F)cc1)c1ccc(F)cc1", "P(=O)(c2ccc(F)cc2)c2ccc(F)cc2"),
    "f": ("2f", "bis(4-chlorophenyl)phosphine oxide", "O=[PH](c1ccc(Cl)cc1)c1ccc(Cl)cc1", "P(=O)(c2ccc(Cl)cc2)c2ccc(Cl)cc2"),
    "g": ("2g", "bis(4-bromophenyl)phosphine oxide", "O=[PH](c1ccc(Br)cc1)c1ccc(Br)cc1", "P(=O)(c2ccc(Br)cc2)c2ccc(Br)cc2"),
    "h": ("2h", "bis(3,5-dimethylphenyl)phosphine oxide", "O=[PH](c1cc(C)cc(C)c1)c1cc(C)cc(C)c1", "P(=O)(c2cc(C)cc(C)c2)c2cc(C)cc(C)c2"),
    "i": ("2i", "di-o-tolylphosphine oxide", "O=[PH](c1ccccc1C)c1ccccc1C", "P(=O)(c2ccccc2C)c2ccccc2C"),
    "j": ("2j", "diethyl phosphite", "O=[PH](OCC)OCC", "P(=O)(OCC)OCC"),
    "k": ("2k", "diisopropyl phosphite", "O=[PH](OC(C)C)OC(C)C", "P(=O)(OC(C)C)OC(C)C"),
    "l": ("2l", "ethyl phenylphosphinate", "O=[PH](OCC)c1ccccc1", "P(=O)(OCC)c2ccccc2"),
}
for key, (label, name, smi, _) in p_groups.items():
    stereo = "p_and_center_unassigned" if key == "l" else "achiral"
    molecule(f"{DOC}:mol:{label}", label, name, ev(1, 3, 5), smi, stereo)

scheme3 = [
    ("4a", "b", 84), ("4b", "c", 88), ("4c", "d", 58), ("4d", "e", 69),
    ("4e", "f", 61), ("4f", "g", 63), ("4g", "i", 80), ("4h", "h", 43),
    ("4i", "j", 31), ("4j", "k", 33), ("4k", "l", 41),
]
for product_label, reagent_key, _yield in scheme3:
    pg = p_groups[reagent_key][3]
    smi = f"CC1=C({pg})C(c2ccccc2)OC1=O"
    product_stereo = "p_and_center_unassigned" if product_label == "4k" else "center_unassigned"
    molecule(f"{DOC}:mol:{product_label}", product_label, f"phosphorylated butenolide {product_label}", ev(5), smi, product_stereo)

# Common reagents, solvents, controls, and derivatization identities.
common = [
    ("cu_oac2", "Cu(OAc)2", "copper(II) acetate"), ("fecl2", "FeCl2", "iron(II) chloride"),
    ("cui", "CuI", "copper(I) iodide"), ("fe_oac2", "Fe(OAc)2", "iron(II) acetate"),
    ("cucl", "CuCl", "copper(I) chloride"), ("cubr2", "CuBr2", "copper(II) bromide"),
    ("cu_me_cn_bf4", "Cu(CH3CN)4BF4", "tetrakis(acetonitrile)copper(I) tetrafluoroborate"),
    ("k2s2o8", "K2S2O8", "potassium persulfate"), ("nh4_s2o8", "(NH4)2S2O8", "ammonium persulfate"),
    ("na2s2o8", "Na2S2O8", "sodium persulfate"), ("tbhp", "TBHP", "tert-butyl hydroperoxide"),
    ("mecn", "CH3CN", "acetonitrile"), ("hoac", "HOAc", "acetic acid"), ("hcooh", "HCOOH", "formic acid"),
    ("tempo", "TEMPO", "2,2,6,6-tetramethylpiperidin-1-yl oxyl"), ("bht", "BHT", "butylated hydroxytoluene"),
    ("zncl2", "ZnCl2", "zinc chloride"), ("alcl3", "AlCl3", "aluminum chloride"),
    ("tbuoh", "t-BuOH", "tert-butanol"), ("tbuok", "t-BuOK", "potassium tert-butoxide"),
    ("pe", "PE", "petroleum ether"), ("ea", "EA", "ethyl acetate"),
]
simple_smiles = {"mecn": "CC#N", "hoac": "CC(=O)O", "hcooh": "O=CO", "tempo": "CC1(C)CCCC(C)(C)N1[O]", "bht": "CC(C)(C)c1cc(C)c(O)c(C(C)(C)C)c1", "tbuoh": "CC(C)(C)O", "tbuok": "CC(C)(C)[O-].[K+]", "ea": "CCOC(C)=O"}
for code, label, name in common:
    molecule(f"{DOC}:mol:{code}", label, name, ev(3, 5, 6), simple_smiles.get(code), "achiral" if simple_smiles.get(code) else "unresolved")

molecule(f"{DOC}:mol:5", "5", "tert-butyl 2-methyl-4-oxo-4-phenylbutanoate", ev(2, 5, 9, image=False), "CC(C(=O)OC(C)(C)C)CC(=O)c1ccccc1", "center_unassigned")
molecule(f"{DOC}:mol:6", "6", "3-methyl-5-phenyl-2(5H)-furanone", ev(3, 6), "CC1=CC(c2ccccc2)OC1=O", "center_unassigned")
molecule(f"{DOC}:mol:bht_2a_adduct", "BHT–2a adduct", "BHT–diphenylphosphine oxide adduct", ev(2, 3, 6), stereo="unresolved")


def set_classification(r, evidence, labels=None, not_applicable=False):
    c = r["classification"]
    c["mapped_reaction_smiles"] = missing(evidence)
    c["model_votes"] = []
    if not_applicable:
        c["bond_changes"] = []
        c["labels"] = inf([], evidence, "No observed transformation; no formation label assigned.")
        c["coverage"] = inf("not_applicable", evidence, "Negative/no-reaction experiment.")
        c["limitations"] = ["No bond formation inferred from an intended but unobserved product."]
        return
    if not labels:
        c["bond_changes"] = []
        c["labels"] = inf([], evidence, "No reliable atom-level net change was asserted.")
        c["coverage"] = inf("unresolved", evidence, "Complete atom mapping was not established.")
        c["limitations"] = ["Source-faithful complete atom mapping was not established."]
        return
    atoms = {"C-P": ("C", "P", 6), "C-O": ("C", "O", 3)}
    changes, out_labels = [], []
    for index, label in enumerate(labels, 1):
        a, b, number = atoms[label]
        cid = f"b{index}"
        changes.append(inf({"change_id": cid, "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}], "before_order": 0, "after_order": 1, "change_type": "formed"}, evidence, "Net bond type inferred from the inspected substrate/product drawings; atom mapping and atom indices are unresolved."))
        out_labels.append({"category_number": number, "code": label, "supporting_change_ids": [cid]})
    c["bond_changes"] = changes
    c["labels"] = inf(out_labels, evidence, "Taxonomy labels generated from the partial inferred net bond changes.")
    c["coverage"] = inf("partial", evidence, "Only visually supported new C–P/C–O connections are represented; bond-order, leaving-group, and hydrogen accounting is incomplete.")
    c["limitations"] = ["No atom-mapped reaction SMILES was asserted.", "Hydrogen provenance and complete bond-order changes remain unresolved."]


def participant(pid, mid, role, evidence, stage="reaction", product_role=None, mass=None, amount=None, volume=None, equiv=None, molpct=None, ref=None, method="original_extraction"):
    p = st.new_participant(pid, mid, role, evidence, stage)
    if product_role is not None:
        p["product_role"] = src(product_role, product_role, evidence)
    vals = {"mass": mass, "amount": amount, "volume": volume, "equivalents": equiv, "mol_percent": molpct}
    for key, value in vals.items():
        if value is not None:
            raw = value.get("raw_text") if isinstance(value, dict) and value.get("raw_text") else value
            p["charges"][key] = st.fact(raw, value, evidence, method=method)
    if ref is not None:
        p["equivalence_reference"] = src(ref, ref, evidence) if method == "original_extraction" else st.fact(ref, ref, evidence, method=method)
    return p


def step(step_id, stage, pids, evidence, operation, temp=None, time=None, atmosphere=None):
    s = st.new_step(step_id, evidence, stage)
    s["participant_ids"] = pids
    s["operation"] = src(operation, operation, evidence)
    if temp is not None:
        s["temperature"] = src(temp.get("raw_text", temp) if isinstance(temp, dict) else temp, temp, evidence)
    if time is not None:
        s["time"] = src(time.get("raw_text", time) if isinstance(time, dict) else time, time, evidence)
    if atmosphere is not None:
        s["atmosphere"] = src(atmosphere, atmosphere, evidence)
    return s


def outcome(oid, kind, targets, evidence, value, basis=None, role="main product", yield_type=None, scope=None):
    o = st.new_outcome(oid, kind, targets, evidence)
    o["value"] = src(value.get("raw_text", value) if isinstance(value, dict) else value, value, evidence)
    o["product_role"] = src(role, role, evidence)
    if basis is not None:
        o["measurement_basis"] = src(basis, basis, evidence)
    if yield_type is not None:
        o["yield_type"] = src(yield_type, yield_type, evidence)
    if scope is not None:
        o["scope"] = src(scope, scope, evidence)
    return o


def reaction(rid, evidence, record_type="specific_experiment", categories=None, location=None, name=None, transformation=None):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = src(categories or [], categories or [], evidence)
    r["location"] = src(location or rid, location or rid, evidence)
    if name:
        r["reaction_name"] = src(name, name, evidence)
    if transformation:
        r["transformation_type"] = inf(transformation, evidence, "Machine description based on the inspected reaction arrow/prose.")
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    return r


def set_solvent(r, evidence, composition, ratio=None, volume=None, method="original_extraction"):
    raw = composition if isinstance(composition, str) else "/".join(composition)
    r["solvent_conditions"]["composition"] = st.fact(raw, composition, evidence, method=method)
    if ratio is not None:
        r["solvent_conditions"]["mixture_ratio"] = st.fact(ratio, ratio, evidence, method=method)
    if volume is not None:
        r["solvent_conditions"]["total_volume"] = st.fact(volume.get("raw_text", volume), volume, evidence, method=method)


reactions = []

# Scheme 1: three literature summaries and one present-work generic transformation.
for suffix, product, title in [
    ("scheme1a_literature_cx", "generic_halobutenolide", "Prior cyclization to beta-halobutenolides (C–X formation)"),
    ("scheme1b_literature_cc", "generic_cc_butenolide", "Prior C–C-forming allenoic-acid cyclizations"),
    ("scheme1c_literature_cs", "generic_cs_butenolide", "Prior C–S-forming allenoic-acid cyclizations"),
]:
    r = reaction(suffix, ev(2), "literature_example", [], f"PDF page 2, Scheme 1{suffix[7]}", title, "literature-reported butenolide synthesis")
    r["participants"] = [participant("generic_input", f"{DOC}:mol:generic_allenoic_acid", "reactant", ev(2)), participant("generic_product", f"{DOC}:mol:{product}", "product", ev(2), product_role="generic product")]
    set_classification(r, ev(2))
    r["review"]["unresolved"] = ["Literature-summary arrow only; cited primary procedures and complete structures were not extracted as experiments of this paper."]
    reactions.append(r)

general_id = f"{DOC}:scheme1d_general_annulation"
r = reaction("scheme1d_general_annulation", ev(1, 2), "generic_scheme", [], "Graphical abstract and PDF page 2 Scheme 1d", "Copper-catalyzed phosphorylation of 2,3-allenoic acids", "annulation/phosphorylation to phosphorylated butenolides")
r["participants"] = [
    participant("generic_1", f"{DOC}:mol:generic_allenoic_acid", "reactant", ev(1, 2)),
    participant("generic_2", f"{DOC}:mol:generic_p_reagent", "reactant", ev(1, 2)),
    participant("generic_product", f"{DOC}:mol:generic_phosphorylated_butenolide", "product", ev(1, 2), product_role="generic main product"),
]
set_classification(r, ev(1, 2), ["C-P", "C-O"])
r["review"]["unresolved"] = ["Generic R-group transformation; example count/up-to yield is a scope summary, not an experiment-specific outcome."]
reactions.append(r)

# Table 1 optimization rows.
table = [
    (1, "fecl2", "k2s2o8", ["CH3CN", "HOAc"], "3:1", 100, 17, None),
    (2, "cui", "k2s2o8", ["CH3CN", "HOAc"], "3:1", 100, None, "n. r."),
    (3, "fe_oac2", "k2s2o8", ["CH3CN", "HOAc"], "3:1", 100, 27, None),
    (4, "cucl", "k2s2o8", ["CH3CN", "HOAc"], "3:1", 100, None, "trace"),
    (5, "cu_oac2", "k2s2o8", ["CH3CN", "HOAc"], "3:1", 100, 80, None),
    (6, "cubr2", "k2s2o8", ["CH3CN", "HOAc"], "3:1", 100, 37, None),
    (7, "cu_me_cn_bf4", "k2s2o8", ["CH3CN", "HOAc"], "3:1", 100, 66, None),
    (8, "cu_oac2", "nh4_s2o8", ["CH3CN", "HOAc"], "3:1", 100, 57, None),
    (9, "cu_oac2", "na2s2o8", ["CH3CN", "HOAc"], "3:1", 100, 45, None),
    (10, "cu_oac2", "tbhp", ["CH3CN", "HOAc"], "3:1", 100, 59, None),
    (11, "cu_oac2", "k2s2o8", ["CH3CN"], None, 100, 47, None),
    (12, "cu_oac2", "k2s2o8", ["CH3CN", "HCOOH"], "3:1", 100, 68, None),
    (13, "cu_oac2", "k2s2o8", ["CH3CN", "HOAc"], "1:1", 100, 73, None),
    (14, "cu_oac2", "k2s2o8", ["CH3CN", "HOAc"], "5:1", 100, 69, None),
    (15, "cu_oac2", "k2s2o8", ["CH3CN", "HOAc"], "3:1", 80, 69, None),
    (16, "cu_oac2", "k2s2o8", ["CH3CN", "HOAc"], "3:1", 100, 65, None),
]
for n, catalyst, oxidant, solvents, ratio, temp_c, yld, qualitative in table:
    evidence = ev(3, 4) if n == 5 else ev(3)
    r = reaction(f"table1_entry_{n:02d}", evidence, categories=["optimization"], location=f"PDF page 3, Table 1 entry {n}", name="Optimization of phosphorylation conditions", transformation="2,3-allenoic-acid phosphorylation/cyclization")
    r["general_scheme_ids"] = [general_id]
    mid1, mid2, mid3 = f"{DOC}:mol:1a", f"{DOC}:mol:2a", f"{DOC}:mol:3a"
    equiv2 = {"number": 1.1, "unit": "equiv", "raw_text": "1.1 equiv"} if n == 16 else {"number": 2, "unit": "equiv", "raw_text": "2 equiv"}
    amount2 = None if n == 16 else {"number": 0.4, "unit": "mmol", "raw_text": "0.4 mmol"}
    r["participants"] = [
        participant("substrate_1a", mid1, "reactant", evidence, amount={"number": 0.2, "unit": "mmol", "raw_text": "0.2 mmol"}, equiv={"number": 1, "unit": "equiv", "raw_text": "1 equiv"}, ref=mid1),
        participant("reagent_2a", mid2, "reactant", evidence, amount=amount2, equiv=equiv2, ref=mid1),
        participant("catalyst", f"{DOC}:mol:{catalyst}", "catalyst", evidence, molpct={"number": 20, "unit": "mol%", "raw_text": "20 mol %"}, ref=mid1),
        participant("oxidant", f"{DOC}:mol:{oxidant}", "reagent", evidence, equiv={"number": 2, "unit": "equiv", "raw_text": "2 equiv"}, ref=mid1),
        participant("product_3a", mid3, "product", evidence, product_role="desired product"),
    ]
    for ix, solv in enumerate(solvents, 1):
        code = {"CH3CN": "mecn", "HOAc": "hoac", "HCOOH": "hcooh"}[solv]
        vol = None
        if n == 5:
            vol = {"number": 1.5 if solv == "CH3CN" else 0.5, "unit": "mL", "raw_text": "1.5 mL" if solv == "CH3CN" else "0.5 mL"}
        r["participants"].append(participant(f"solvent_{ix}", f"{DOC}:mol:{code}", "solvent", evidence, volume=vol))
    set_solvent(r, evidence, solvents, ratio, {"number": 2, "unit": "mL", "raw_text": "2.0 mL"})
    all_pids = [p["participant_id"] for p in r["participants"]]
    r["steps"] = [step("reaction", "reaction", all_pids, evidence, "mix/heat and stir", {"number": temp_c, "unit": "°C", "raw_text": f"{temp_c} °C"}, {"number": 12, "unit": "h", "raw_text": "12 h"}, "Ar")]
    if n == 5:
        r["steps"] += [
            step("cool_and_concentrate", "workup", all_pids, evidence, "cool to room temperature and concentrate under reduced pressure", "room temperature"),
            step("flash_chromatography", "purification", ["product_3a"], evidence, "flash chromatography, petroleum ether/ethyl acetate = 2:1"),
        ]
        r["outcomes"].append(outcome("product_amount", "product_amount", [mid3], evidence, {"number": 59.9, "unit": "mg", "raw_text": "59.9 mg"}, "isolated after flash chromatography"))
        r["procedure"]["original_text"] = src("2,3-Allenoic Acid 1a (34.8 mg, 0.2 mmol), 2a (0.4 mmol, 2 equiv), Cu(OAc)2 (0.04 mmol, 20 mol %), K2S2O8 (0.4 mmol, 2 equiv), CH3CN (1.5 mL), and HOAc (0.5 mL) were mixed ... 100 °C ... 12 h under argon ... flash chromatography (petroleum ether/ethyl acetate = 2:1).", "representative general procedure for 3a", evidence)
        # Add detailed masses/amounts without changing the table row identity.
        r["participants"][0]["charges"]["mass"] = src("34.8 mg", {"number": 34.8, "unit": "mg"}, evidence)
        r["participants"][2]["charges"]["amount"] = src("0.04 mmol", {"number": 0.04, "unit": "mmol"}, evidence)
        r["participants"][3]["charges"]["amount"] = src("0.4 mmol", {"number": 0.4, "unit": "mmol"}, evidence)
    else:
        r["procedure"]["general_procedure_reference"] = src("Table 1 footnote a", "Table 1 common conditions with row-specific modifications", evidence)
    if yld is not None:
        r["outcomes"].insert(0, outcome("yield", "yield", [mid3], evidence, {"number": yld, "unit": "%", "raw_text": f"{yld}%"}, "isolated yield", yield_type="isolated"))
        set_classification(r, evidence, ["C-P", "C-O"])
    elif qualitative == "n. r.":
        r["outcomes"] = [outcome("no_reaction", "no_reaction", [mid3], evidence, "n. r.", None, role="intended product")]
        set_classification(r, evidence, not_applicable=True)
    else:
        r["outcomes"] = [outcome("trace", "trace", [mid3], evidence, "trace", None, role="desired product")]
        set_classification(r, evidence, ["C-P", "C-O"])
    r["review"]["unresolved"] = ["Structures and all quantitative facts remain machine-extracted and require expert review."]
    reactions.append(r)

# Scheme 2: 3a is the same experiment as Table 1 entry 5; add only 3b–3r.
for label, yld in scope3_yields.items():
    if label == "a":
        continue
    evidence = ev(4, 6) if label in ("c", "d") else ev(3, 4) if label == "b" else ev(4)
    r = reaction(f"scheme2_3{label}", evidence, categories=["substrate_scope"], location=f"PDF page 4, Scheme 2 product 3{label}", name="Substrate scope relative to 2,3-allenoic acids", transformation="phosphorylative annulation to a butenolide")
    r["general_scheme_ids"] = [general_id]
    m1, m2, mp = f"{DOC}:mol:1{label}", f"{DOC}:mol:2a", f"{DOC}:mol:3{label}"
    r["participants"] = [
        participant(f"substrate_1{label}", m1, "reactant", evidence, amount={"number": 0.2, "unit": "mmol", "raw_text": "0.2 mmol"}, equiv={"number": 1, "unit": "equiv", "raw_text": "1 equiv"}, ref=m1, method="inherited_general_condition"),
        participant("reagent_2a", m2, "reactant", evidence, amount={"number": 0.4, "unit": "mmol", "raw_text": "0.4 mmol"}, equiv={"number": 2, "unit": "equiv", "raw_text": "2 equiv"}, ref=m1, method="inherited_general_condition"),
        participant("cu_oac2", f"{DOC}:mol:cu_oac2", "catalyst", evidence, molpct={"number": 20, "unit": "mol%", "raw_text": "20 mol %"}, ref=m1, method="inherited_general_condition"),
        participant("k2s2o8", f"{DOC}:mol:k2s2o8", "reagent", evidence, equiv={"number": 2, "unit": "equiv", "raw_text": "2 equiv"}, ref=m1, method="inherited_general_condition"),
        participant("mecn", f"{DOC}:mol:mecn", "solvent", evidence), participant("hoac", f"{DOC}:mol:hoac", "solvent", evidence),
        participant(f"product_3{label}", mp, "product", evidence, product_role="main product"),
    ]
    set_solvent(r, evidence, ["CH3CN", "HOAc"], "3:1", {"number": 2, "unit": "mL", "raw_text": "2.0 mL"}, "inherited_general_condition")
    pids = [p["participant_id"] for p in r["participants"]]
    r["steps"] = [step("reaction", "reaction", pids, evidence, "apply Scheme 2 common conditions", {"number": 100, "unit": "°C", "raw_text": "100 °C"}, {"number": 12, "unit": "h", "raw_text": "12 h"}, "Ar")]
    r["procedure"]["general_procedure_reference"] = src("Scheme 2 footnote a", "Scheme 2 common conditions", evidence)
    r["outcomes"] = [outcome("yield", "yield", [mp], evidence, {"number": yld, "unit": "%", "raw_text": f"{yld}%"}, "isolated yield", yield_type="isolated")]
    amounts = {"b": 52.8, "c": 36.1, "d": 56.6}
    if label in amounts:
        r["outcomes"].append(outcome("product_amount", "product_amount", [mp], evidence, {"number": amounts[label], "unit": "mg", "raw_text": f"{amounts[label]} mg"}, "isolated product after chromatography"))
    set_classification(r, evidence, ["C-P", "C-O"])
    r["review"]["unresolved"] = ["Absolute/relative configuration at the product tetrahedral center is not reported; no single stereoisomer was invented."]
    reactions.append(r)

# Scheme 3 phosphorus-reagent scope.
for product_label, reagent_key, yld in scheme3:
    evidence = ev(5)
    r = reaction(f"scheme3_{product_label}", evidence, categories=["substrate_scope"], location=f"PDF page 5, Scheme 3 product {product_label}", name="Substrate scope relative to phosphorus reagent", transformation="phosphorylative annulation to a butenolide")
    r["general_scheme_ids"] = [general_id]
    m1, m2, mp = f"{DOC}:mol:1a", f"{DOC}:mol:{p_groups[reagent_key][0]}", f"{DOC}:mol:{product_label}"
    r["participants"] = [
        participant("substrate_1a", m1, "reactant", evidence, amount={"number": 0.2, "unit": "mmol", "raw_text": "0.2 mmol"}, equiv={"number": 1, "unit": "equiv", "raw_text": "1 equiv"}, ref=m1, method="inherited_general_condition"),
        participant("p_reagent", m2, "reactant", evidence, amount={"number": 0.4, "unit": "mmol", "raw_text": "0.4 mmol"}, equiv={"number": 2, "unit": "equiv", "raw_text": "2 equiv"}, ref=m1, method="inherited_general_condition"),
        participant("cu_oac2", f"{DOC}:mol:cu_oac2", "catalyst", evidence, molpct={"number": 20, "unit": "mol%", "raw_text": "20 mol %"}, ref=m1, method="inherited_general_condition"),
        participant("k2s2o8", f"{DOC}:mol:k2s2o8", "reagent", evidence, equiv={"number": 2, "unit": "equiv", "raw_text": "2 equiv"}, ref=m1, method="inherited_general_condition"),
        participant("mecn", f"{DOC}:mol:mecn", "solvent", evidence), participant("hoac", f"{DOC}:mol:hoac", "solvent", evidence),
        participant("product", mp, "product", evidence, product_role="main product"),
    ]
    set_solvent(r, evidence, ["CH3CN", "HOAc"], "3:1", {"number": 2, "unit": "mL", "raw_text": "2.0 mL"}, "inherited_general_condition")
    r["steps"] = [step("reaction", "reaction", [p["participant_id"] for p in r["participants"]], evidence, "apply Scheme 3 common conditions", {"number": 100, "unit": "°C", "raw_text": "100 °C"}, {"number": 12, "unit": "h", "raw_text": "12 h"}, "Ar")]
    r["procedure"]["general_procedure_reference"] = src("Scheme 3 footnote a", "Scheme 3 common conditions", evidence)
    r["outcomes"] = [outcome("yield", "yield", [mp], evidence, {"number": yld, "unit": "%", "raw_text": f"{yld}%"}, "isolated yield", yield_type="isolated")]
    set_classification(r, evidence, ["C-P", "C-O"])
    r["review"]["unresolved"] = ["Product stereocenter configuration is unreported; 4k may also have unassigned phosphorus stereochemistry."]
    reactions.append(r)

# Scheme 4 gram scale.
evidence = ev(3, 4, 5)
r = reaction("scheme4_gram_scale_3a", evidence, categories=["scale_up"], location="PDF page 5 Scheme 4; procedure spans PDF pages 3–4", name="Gram-scale synthesis of 3a", transformation="scale-up phosphorylative annulation")
r["general_scheme_ids"] = [general_id]
m1, m2, mp = f"{DOC}:mol:1a", f"{DOC}:mol:2a", f"{DOC}:mol:3a"
r["participants"] = [
    participant("substrate_1a", m1, "reactant", evidence, mass={"number": 1.134, "unit": "g", "raw_text": "1.134 g"}, amount={"number": 6, "unit": "mmol", "raw_text": "6 mmol"}, equiv={"number": 1, "unit": "equiv", "raw_text": "1 equiv"}, ref=m1),
    participant("reagent_2a", m2, "reactant", evidence, mass={"number": 2.426, "unit": "g", "raw_text": "2.426 g"}, amount={"number": 12, "unit": "mmol", "raw_text": "12 mmol"}, equiv={"number": 2, "unit": "equiv", "raw_text": "2 equiv"}, ref=m1),
    participant("cu_oac2", f"{DOC}:mol:cu_oac2", "catalyst", evidence, mass={"number": 0.218, "unit": "g", "raw_text": "0.218 g"}, molpct={"number": 20, "unit": "mol%", "raw_text": "20 mol %"}, ref=m1),
    participant("k2s2o8", f"{DOC}:mol:k2s2o8", "reagent", evidence, mass={"number": 3.243, "unit": "g", "raw_text": "3.243 g"}, equiv={"number": 2, "unit": "equiv", "raw_text": "2 equiv"}, ref=m1),
    participant("mecn", f"{DOC}:mol:mecn", "solvent", evidence, volume={"number": 45, "unit": "mL", "raw_text": "45.0 mL"}),
    participant("hoac", f"{DOC}:mol:hoac", "solvent", evidence, volume={"number": 15, "unit": "mL", "raw_text": "15.0 mL"}),
    participant("product_3a", mp, "product", evidence, product_role="main product"),
]
r["reaction_scale"] = src("1a: 6 mmol", {"reference_molecule_id": m1, "amount": {"number": 6, "unit": "mmol"}}, evidence)
set_solvent(r, evidence, ["CH3CN", "HOAc"], "3:1", {"number": 60, "unit": "mL", "raw_text": "60 mL"})
pids = [p["participant_id"] for p in r["participants"]]
r["steps"] = [step("reaction", "reaction", pids, evidence, "stir in oil bath", {"number": 100, "unit": "°C", "raw_text": "100 °C"}, {"number": 12, "unit": "h", "raw_text": "12 h"}, "Ar"), step("concentrate", "workup", pids, evidence, "concentrate under reduced pressure"), step("purify", "purification", ["product_3a"], evidence, "flash chromatography, PE/EA = 2:1")]
r["outcomes"] = [outcome("yield", "yield", [mp], evidence, {"number": 62, "unit": "%", "raw_text": "62%"}, "isolated after flash chromatography", yield_type="isolated"), outcome("product_amount", "product_amount", [mp], evidence, {"number": 1.392, "unit": "g", "raw_text": "1.392 g (Scheme: 1.39 g)"}, "isolated after flash chromatography")]
set_classification(r, evidence, ["C-P", "C-O"])
r["review"]["unresolved"] = ["Machine structure candidates and reagent identities require expert review; 1.392 g and scheme-rounded 1.39 g are treated as consistent expressions, not a conflict."]
reactions.append(r)

# Scheme 4 ring opening; later extracted text is contextual and conflicts with the page-5 drawing.
evidence = ev(2, 5) + ev(9, image=False)
r = reaction("scheme4_ring_opening_3a_to_5", evidence, categories=["derivatization"], location="PDF page 5 Scheme 4; later extracted procedure on PDF page 9 used as context only", name="Open ring of 3a", transformation="butenolide ring opening to a 1,4-dicarbonyl compound")
r["participants"] = [
    participant("reactant_3a", mp, "reactant", evidence, amount={"number": 0.4, "unit": "mmol", "raw_text": "0.4 mmol"}, equiv={"number": 1, "unit": "equiv", "raw_text": "1 equiv"}, ref=mp),
    participant("zncl2", f"{DOC}:mol:zncl2", "reagent", evidence, amount={"number": 0.4, "unit": "mmol", "raw_text": "0.4 mmol"}, equiv={"number": 1, "unit": "equiv", "raw_text": "1 equiv"}, ref=mp),
    participant("tbuoh", f"{DOC}:mol:tbuoh", "solvent", ev(5), volume={"number": 4, "unit": "mL", "raw_text": "4 mL"}),
    participant("tbuok_text", f"{DOC}:mol:tbuok", "unknown", ev(9, image=False), volume={"number": 0.4, "unit": "mL", "raw_text": "0.4 mL (extracted text)"}),
    participant("product_5", f"{DOC}:mol:5", "product", evidence, product_role="main product"),
    participant("pe", f"{DOC}:mol:pe", "purification_agent", ev(9, image=False), stage="purification"),
    participant("ea", f"{DOC}:mol:ea", "purification_agent", ev(9, image=False), stage="purification"),
]
r["solvent_conditions"]["composition"] = conflict([src("t-BuOH", "t-BuOH", ev(5)), src("t-BuOK", "t-BuOK", ev(9, image=False))], evidence)
r["solvent_conditions"]["total_volume"] = conflict([src("4 mL", {"number": 4, "unit": "mL"}, ev(5)), src("0.4 mL", {"number": 0.4, "unit": "mL"}, ev(9, image=False))], evidence)
r["steps"] = [step("reaction", "reaction", ["reactant_3a", "zncl2", "tbuoh", "tbuok_text"], evidence, "heat/stir", {"number": 135, "unit": "°C", "raw_text": "135 °C"}), step("concentrate", "workup", ["product_5"], ev(9, image=False), "concentrate under reduced pressure"), step("purify", "purification", ["product_5", "pe", "ea"], ev(9, image=False), "flash chromatography, PE/EA = 4:1")]
r["steps"][0]["time"] = conflict([src("12 h", {"number": 12, "unit": "h"}, ev(5)), src("overnight", "overnight", ev(9, image=False))], evidence)
r["outcomes"] = [outcome("yield", "yield", [f"{DOC}:mol:5"], evidence, {"number": 42, "unit": "%", "raw_text": "42%"}, "isolated after flash chromatography", yield_type="isolated"), outcome("product_amount", "product_amount", [f"{DOC}:mol:5"], ev(9, image=False), {"number": 42, "unit": "mg", "raw_text": "42 mg"}, "isolated after flash chromatography")]
set_classification(r, evidence)
r["notes"] = [src("Scheme 4 visually reports t-BuOH (4 mL), 135 °C, 12 h; later extracted procedure text reads t-BuOK (0.4 mL) and overnight.", "source conflict retained without correction", evidence)]
r["review"]["unresolved"] = ["The t-BuOH/t-BuOK and 4/0.4 mL discrepancy requires visual review of the later procedure page; that page is outside this batch's claimed visual coverage.", "Complete bond changes for ring opening were not atom mapped."]
reactions.append(r)

# Scheme 5 controls. The BHT record also preserves the detected adduct.
def control_record(suffix, extra_mid, extra_label, outcome_kind, outcome_value, labels, catalyst=True, oxidant=True, substrate="1a", product="3a", extra_equiv=None, panel=None):
    evidence = ev(2, 3, 6)
    panel = panel or suffix[-1]
    r = reaction(suffix, evidence, categories=["control", "mechanistic_study"], location=f"PDF page 6 Scheme 5 ({panel})", name="Control experiment", transformation="mechanistic control of phosphorylative annulation")
    ms = f"{DOC}:mol:{substrate}"; m2 = f"{DOC}:mol:2a"; mp2 = f"{DOC}:mol:{product}"
    r["participants"] = [participant("substrate", ms, "reactant", evidence), participant("reagent_2a", m2, "reactant", evidence), participant("intended_product", mp2, "product", evidence, product_role="intended product")]
    if catalyst:
        r["participants"].append(participant("cu_oac2", f"{DOC}:mol:cu_oac2", "catalyst", evidence, molpct={"number": 20, "unit": "mol%", "raw_text": "20 mol %"}, ref=ms))
    if oxidant:
        r["participants"].append(participant("k2s2o8", f"{DOC}:mol:k2s2o8", "reagent", evidence, equiv={"number": 2, "unit": "equiv", "raw_text": "2 equiv"}, ref=ms))
    if extra_mid:
        r["participants"].append(participant(extra_label, f"{DOC}:mol:{extra_mid}", "additive" if extra_mid in ("tempo", "bht") else "catalyst", evidence, equiv={"number": extra_equiv, "unit": "equiv", "raw_text": f"{extra_equiv} equiv"} if extra_equiv else None, molpct={"number": 20, "unit": "mol%", "raw_text": "20 mol %"} if extra_mid in ("zncl2", "alcl3") else None, ref=ms))
    r["participants"] += [participant("mecn", f"{DOC}:mol:mecn", "solvent", evidence), participant("hoac", f"{DOC}:mol:hoac", "solvent", evidence)]
    set_solvent(r, evidence, ["CH3CN", "HOAc"], "3:1", {"number": 2, "unit": "mL", "raw_text": "2 mL"}, "inherited_general_condition")
    r["steps"] = [step("reaction", "reaction", [p["participant_id"] for p in r["participants"]], evidence, "apply depicted control conditions", {"number": 100, "unit": "°C", "raw_text": "100 °C"}, {"number": 12, "unit": "h", "raw_text": "12 h"}, "Ar")]
    r["procedure"]["general_procedure_reference"] = src("standard conditions, with depicted omission/substitution", "standard conditions modified as shown in Scheme 5", evidence)
    r["outcomes"] = [outcome("control_outcome", outcome_kind, [mp2], evidence, outcome_value, None, role="intended product")]
    set_classification(r, evidence, labels, not_applicable=outcome_kind == "no_reaction")
    r["review"]["unresolved"] = ["Detection/measurement basis for the Scheme 5 product outcome is not stated."]
    return r

r = control_record("scheme5_a", "tempo", "tempo", "no_reaction", {"number": 0, "unit": "%", "raw_text": "0%"}, None, extra_equiv=2)
reactions.append(r)
r = control_record("scheme5_b", "bht", "bht", "no_reaction", "completely suppressed", None, extra_equiv=2)
r["participants"].append(participant("bht_2a_adduct", f"{DOC}:mol:bht_2a_adduct", "product", ev(2, 3, 6), product_role="detected radical-trapping adduct"))
r["outcomes"].append(outcome("adduct_detected", "other", [f"{DOC}:mol:bht_2a_adduct"], ev(2, 3, 6), "detected", "HRMS", role="radical-trapping adduct"))
r["review"]["unresolved"].append("No adduct yield or amount is reported; the HRMS observation is qualitative.")
reactions.append(r)
r = control_record("scheme5_c", None, "", "no_reaction", {"number": 0, "unit": "%", "raw_text": "0%"}, None, catalyst=True, oxidant=False)
reactions.append(r)
r = control_record("scheme5_d", None, "", "yield", {"number": 14, "unit": "%", "raw_text": "14%"}, ["C-P", "C-O"], catalyst=False, oxidant=True)
reactions.append(r)
for acid in ("zncl2", "alcl3"):
    r = control_record(f"scheme5_e_{acid}", acid, acid, "no_reaction", {"number": 0, "unit": "%", "raw_text": "0%"}, None, catalyst=False, oxidant=True, panel="e")
    reactions.append(r)
r = control_record("scheme5_f", None, "", "no_reaction", {"number": 0, "unit": "%", "raw_text": "0%"}, None, catalyst=True, oxidant=True, substrate="6")
reactions.append(r)

# Proposed radical pathway is described in page-3 prose; the Scheme-6 drawing is outside visual coverage.
r = reaction("scheme6_proposed_radical_pathway", ev(3), "proposed_mechanism", [], "PDF page 3 mechanism discussion referring to Scheme 6", "Proposed free-radical pathway", "proposed P-centered radical addition/oxidation/cyclization")
r["participants"] = [participant("substrate_1a", f"{DOC}:mol:1a", "reactant", ev(3)), participant("reagent_2a", f"{DOC}:mol:2a", "reactant", ev(3)), participant("product_3", f"{DOC}:mol:3a", "product", ev(3), product_role="overall product")]
set_classification(r, ev(3), ["C-P", "C-O"])
r["procedure"]["original_text"] = src("The authors propose P-centered radical generation, radical addition to 1a, oxidation to an allyl cation, and intramolecular nucleophilic attack to product 3.", "author-proposed pathway summary", ev(3))
r["review"]["unresolved"] = ["This is an author-proposed mechanism, not a set of isolated intermediates; the Scheme 6 drawing lies outside the claimed visual page range."]
reactions.append(r)

bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": f"{DOC}:page1_graphical_abstract", "document_id": DOC, "kind": "duplicate_generic_scheme", "evidence_ids": ev(1), "note": "Graphical abstract duplicates the present-work generic annulation and its 29-example/up-to-88% scope summary."},
    {"context_id": f"{DOC}:page2_figure1", "document_id": DOC, "kind": "nonreaction_figure", "evidence_ids": ev(2), "note": "Commercial drugs and bioactive furanones; not synthesis experiments of this paper."},
    {"context_id": f"{DOC}:pages4_6_characterization", "document_id": DOC, "kind": "characterization", "evidence_ids": ev(4, 6), "note": "NMR/HRMS/melting-point blocks support compound identity/amount context and are not separate synthesis reactions."},
    {"context_id": f"{DOC}:si_notice", "document_id": DOC, "kind": "source_availability", "evidence_ids": [EVT[9]], "note": "Article states that SI exists online, but no SI file was explicitly provided or associated in this workspace."},
    {"context_id": f"{DOC}:page9_ring_opening_context", "document_id": DOC, "kind": "outside_batch_text_context", "evidence_ids": [EVT[9]], "note": "Later extracted text was used only to supplement the page-5 ring-opening identity; visual coverage remains pages 1–6."},
]

(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/pages_0001_0006_structure_diagnostics.json").write_text(json.dumps(diagnostics, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/pages_0001_0006_inventory.json").write_text(json.dumps({
    "document_id": DOC,
    "claimed_visual_pages": [1, 2, 3, 4, 5, 6],
    "reaction_count": len(reactions),
    "molecule_count": len(molecules),
    "deduplication": {"scheme2_3a": "merged with table1_entry_05"},
    "counts": {
        "literature_examples": 3,
        "generic_schemes": 1,
        "optimization_rows": 16,
        "new_scheme2_scope_rows": 17,
        "scheme3_scope_rows": 11,
        "scale_up": 1,
        "derivatization": 1,
        "control_experiments": 7,
        "proposed_mechanisms": 1,
    },
}, ensure_ascii=False, indent=2) + "\n")

rid = lambda suffix: f"{DOC}:{suffix}"
table_ids = [rid(f"table1_entry_{n:02d}") for n in range(1, 17)]
scheme2_ids = [rid(f"scheme2_3{x}") for x in "bcdefghijklmnopqr"]
scheme3_ids = [rid(f"scheme3_4{x}") for x in "abcdefghijk"]
control_ids = [rid("scheme5_a"), rid("scheme5_b"), rid("scheme5_c"), rid("scheme5_d"), rid("scheme5_e_zncl2"), rid("scheme5_e_alcl3"), rid("scheme5_f")]
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 1,
            "state": "attempted",
            "regions": [
                {"region": "title, abstract, and graphical abstract", "finding": "Present-work generic Cu-catalyzed annulation; 29 examples/up to 88% is a scope summary, not an individual experiment.", "evidence_ids": [EVI[1], EVT[1]]},
                {"region": "Introduction columns", "finding": "Background/prior work and opening optimization discussion; no additional row-level experiment beyond the later table.", "evidence_ids": [EVI[1], EVT[1]]},
            ],
            "record_ids": [rid("scheme1d_general_annulation"), rid("table1_entry_01"), rid("table1_entry_05")],
            "unresolved": ["Graphical abstract R-group structures are generic and remain templates.", "No SI file is associated."],
        },
        {
            "page": 2,
            "state": "attempted",
            "regions": [
                {"region": "Figure 1", "finding": "Commercial drugs/bioactive furanones; nonreaction context, not experiments of this paper.", "evidence_ids": [EVI[2], EVT[2]]},
                {"region": "Scheme 1a–d", "finding": "Three prior-work generic arrows (C–X/C–C/C–S) and one present-work generic C–P annulation; prior work not counted as current experiments.", "evidence_ids": [EVI[2], EVT[2]]},
                {"region": "Results/discussion prose", "finding": "Scope, gram-scale, derivatization, and control summaries cross-reference later schemes; used only as duplicate/context evidence.", "evidence_ids": [EVI[2], EVT[2]]},
            ],
            "record_ids": [rid("scheme1a_literature_cx"), rid("scheme1b_literature_cc"), rid("scheme1c_literature_cs"), rid("scheme1d_general_annulation"), rid("scheme4_ring_opening_3a_to_5")] + control_ids,
            "unresolved": ["Primary procedures for cited prior work were not supplied and were not reconstructed.", "Control-summary prose does not state yield measurement methods."],
        },
        {
            "page": 3,
            "state": "attempted",
            "regions": [
                {"region": "Table 1 including reaction drawing, all 16 rows, and footnotes a–d", "finding": "Sixteen independent optimization experiments; entry 2 is n.r., entry 4 trace, entry 15 is 80 °C, and entry 16 uses 1.1 equiv 2a.", "evidence_ids": [EVI[3], EVT[3]]},
                {"region": "left-column control/mechanism discussion", "finding": "Author discussion of Scheme 5 controls and proposed radical pathway; proposed intermediates are not isolated experiments.", "evidence_ids": [EVI[3], EVT[3]]},
                {"region": "Experimental section and gram-scale procedure start", "finding": "Representative preparation of 3a and start of the distinct 6 mmol scale-up; characterization text is not a separate reaction.", "evidence_ids": [EVI[3], EVT[3]]},
            ],
            "record_ids": table_ids + [rid("scheme4_gram_scale_3a"), rid("scheme6_proposed_radical_pathway")] + control_ids,
            "unresolved": ["Scheme 5/6 drawings themselves occur later; page-3 prose is retained as discussion evidence.", "Control outcome measurement bases are not stated."],
        },
        {
            "page": 4,
            "state": "attempted",
            "regions": [
                {"region": "Scheme 2 reaction arrow, 18 products 3a–3r, and footnotes", "finding": "All 18 scope arrows inventoried; 3a is the same 0.2 mmol/80% experiment as Table 1 entry 5 and is merged rather than duplicated.", "evidence_ids": [EVI[4], EVT[4]]},
                {"region": "gram-scale continuation", "finding": "Continuation of the 6 mmol 3a scale-up: 100 °C/Ar/12 h, concentration, PE/EA 2:1, 1.392 g/62%.", "evidence_ids": [EVI[4], EVT[4]]},
                {"region": "3a/3b characterization", "finding": "Product identity/amount context only; NMR/HRMS peak lists are not separate synthesis experiments.", "evidence_ids": [EVI[4], EVT[4]]},
            ],
            "record_ids": [rid("table1_entry_05")] + scheme2_ids + [rid("scheme4_gram_scale_3a")],
            "unresolved": ["Absolute/relative configuration at product tetrahedral centers is not reported.", "Machine SMILES are candidates, not source-verified structures."],
        },
        {
            "page": 5,
            "state": "attempted",
            "regions": [
                {"region": "Scheme 3 arrow, 11 products 4a–4k, and footnotes", "finding": "Eleven independent phosphorus-reagent scope experiments, all with isolated yields.", "evidence_ids": [EVI[5], EVT[5]]},
                {"region": "Scheme 4 upper arrow", "finding": "Distinct 6 mmol scale-up to 3a, 62% and scheme-rounded 1.39 g.", "evidence_ids": [EVI[5], EVT[5]]},
                {"region": "Scheme 4 lower arrow", "finding": "3a ring opening to 5 at 135 °C/12 h, ZnCl2 1 equiv, t-BuOH 4 mL, 42%; later extracted procedure text conflicts and is retained separately.", "evidence_ids": [EVI[5], EVT[5]]},
            ],
            "record_ids": scheme3_ids + [rid("scheme4_gram_scale_3a"), rid("scheme4_ring_opening_3a_to_5")],
            "unresolved": ["Later page-9 extracted text says t-BuOK 0.4 mL/overnight versus page-5 t-BuOH 4 mL/12 h; page 9 was not visually covered in this batch.", "Complete atom mapping for ring opening remains unresolved."],
        },
        {
            "page": 6,
            "state": "attempted",
            "regions": [
                {"region": "Scheme 5 controls a–f", "finding": "TEMPO, BHT/adduct, no-oxidant, no-copper, two Lewis-acid substitutions, and precyclized-6 controls; ZnCl2 and AlCl3 are separate conditions, yielding seven records.", "evidence_ids": [EVI[6], EVT[6]]},
                {"region": "3c/3d characterization", "finding": "Isolated masses and yields support existing Scheme 2 identities; spectra are characterization context, not reactions.", "evidence_ids": [EVI[6], EVT[6]]},
            ],
            "record_ids": control_ids + [rid("scheme2_3c"), rid("scheme2_3d")],
            "unresolved": ["Scheme 5 does not state the product-yield measurement basis.", "The BHT–2a adduct is qualitatively detected by HRMS; no amount/yield is reported."],
        },
    ],
}
(RUN / "coverage/pages_0001_0006.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")

# Mark only the six claimed page tasks attempted; later pages remain pending.
tasks_path = RUN / "tasks.json"
tasks = json.loads(tasks_path.read_text())
for task in tasks:
    if task.get("document_id") == DOC and task.get("evidence_id") in {EVT[n] for n in range(1, 7)}:
        task["state"] = "attempted"
        task["reason"] = "Original page image and extracted text inspected; detailed region inventory is in coverage/pages_0001_0006.json."
tasks_path.write_text(json.dumps(tasks, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"reactions": len(reactions), "molecules": len(molecules)}, ensure_ascii=False))
