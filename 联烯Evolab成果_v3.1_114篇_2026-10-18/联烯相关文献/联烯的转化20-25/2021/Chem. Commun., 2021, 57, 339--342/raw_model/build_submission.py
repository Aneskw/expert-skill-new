#!/usr/bin/env python3
"""Source-specific v3.1 bundle builder for Chem. Commun. 2021, 57, 339–342."""
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_84632c1bd8ad7a3e"
RUN = ROOT / "runs/current"
E = {p: [f"{DOC}:page{p}:image", f"{DOC}:page{p}:text"] for p in range(1, 5)}
EI = {p: [f"{DOC}:page{p}:image"] for p in range(1, 5)}


def F(raw=None, norm=None, ev=(), method="original_extraction", status=None, derivation=None):
    return st.fact(raw, norm, ev, method=method, status=status, derivation=derivation)


def Q(number, unit, **extra):
    return {"number": number, "unit": unit, **extra}


def mid(label):
    return f"{DOC}:mol_{label}"


def rid(label):
    return f"{DOC}:{label}"


bundle = json.loads((RUN / "raw_model/standard_init.json").read_text())
molecules = {}
reactions = []


def add_molecule(label, name, page, *, smiles=None, chirality="unspecified",
                 geometry="unspecified", representation=None, assignment=None,
                 mixture=None, axial=None, note=None):
    ev = E[page]
    m = st.new_molecule(mid(label), DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, ev)
    m["names"]["original"] = F(name, name, ev)
    m["names"]["english"] = F(name, name, ev)
    m["names"]["abbreviations"] = F([], [], ev)
    m["names"]["aliases"] = F([], [], ev)
    m["cas_numbers"] = F(status="unreported")
    m["identity_note"] = F(note, note, ev) if note else F(status="unreported")
    m["external_properties"] = []
    if smiles:
        diag = st.stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = F(None, smiles, ev, method="structure_derivation",
                         derivation="Candidate SMILES transcribed from the source drawing/name; parsing does not verify source identity.")
        rep = representation or ("partial" if any(x["specified"] == "Unspecified" for x in diag["potential_stereo"]) else "encoded")
        m["stereochemistry"]["isomeric_smiles"] = F(
            None, diag["canonical_isomeric_smiles"], ev, method="structure_derivation",
            derivation="RDKit canonical isomeric SMILES calculated from the unverified candidate graph.")
        m["stereochemistry"]["identity_key"] = F(
            None, diag["stereo_graph_key"], ev, method="structure_derivation",
            derivation="RDKit stereo graph key calculated from the unverified candidate graph.")
        m["stereochemistry"]["representation_status"] = F(
            rep, rep, ev, method="inference", derivation="Assessed from the candidate graph and source stereo annotations.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = F(status="unreported")
        m["stereochemistry"]["isomeric_smiles"] = F(status="unreported")
        m["stereochemistry"]["identity_key"] = F(status="unreported")
        rep = representation or "not_encoded"
        m["stereochemistry"]["representation_status"] = F(rep, rep, ev)
    m["stereochemistry"]["chirality_status"] = F(chirality, chirality, ev)
    m["stereochemistry"]["geometry_status"] = F(geometry, geometry, ev)
    m["stereochemistry"]["atom_centers"] = F(
        "not assigned in the supplied pages" if chirality in ("unspecified", "partially_specified") else None,
        [] if chirality == "achiral" else None, ev,
        method="structure_derivation" if chirality == "achiral" else "original_extraction",
        derivation="No tetrahedral stereogenic center is present in the candidate graph." if chirality == "achiral" else None)
    m["stereochemistry"]["double_bonds"] = F(
        assignment if geometry in ("specified", "mixture", "partially_specified") else None,
        assignment if geometry in ("specified", "mixture", "partially_specified") else None, ev)
    m["stereochemistry"]["axial"] = F(
        axial, axial, ev) if axial is not None else F(status="unreported")
    m["stereochemistry"]["mixture"] = F(mixture, mixture, ev) if mixture is not None else F(status="unreported")
    m["stereochemistry"]["assignment"] = F(assignment, assignment, ev) if assignment is not None else F(status="unreported")
    molecules[label] = m
    return m


def add_template(label, name, page, fixed, variables):
    ev = E[page]
    m = add_molecule(label, name, page, chirality="unresolved", geometry="unresolved", representation="not_encoded")
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": F(fixed, {"description": fixed}, ev),
        "variables": [{
            "label": v,
            "attachment_points": F(f"source drawing attachment for {v}", f"source drawing attachment for {v}", ev),
            "definition": F("R-group as drawn", "R-group as drawn", ev),
            "allowed_values": F(status="unreported"),
            "constraints": F(status="unreported"),
        } for v in variables],
        "instance_mapping": [],
        "reason": F("Generic R-group drawing", "Source depicts a reaction family rather than one concrete molecule.", ev),
    }
    return m


def add_participant(r, label, role, ev, *, amount=None, mass=None, volume=None,
                    equiv=None, molpct=None, reference=None, stage="reaction", product_role=None):
    p = st.new_participant(f"p{len(r['participants'])+1}_{label}", mid(label), role, ev, stage=stage)
    if product_role:
        p["product_role"] = F(product_role, product_role, ev)
    for key, val in (("amount", amount), ("mass", mass), ("volume", volume),
                     ("equivalents", equiv), ("mol_percent", molpct)):
        if val is not None:
            p["charges"][key] = F(val.get("raw_text", f"{val['number']} {val['unit']}"), val, ev)
    if reference:
        p["equivalence_reference"] = F(reference, mid(reference), ev)
    r["participants"].append(p)
    return p["participant_id"]


def set_common(r, page, categories, location, name, transformation):
    ev = E[page]
    r["record_type"] = "specific_experiment"
    r["experiment_categories"] = F(categories, categories, ev)
    r["location"] = F(location, {"page": page, "location": location}, ev)
    r["reaction_name"] = F(name, name, ev)
    r["transformation_type"] = F(None, transformation, ev, method="inference",
                                  derivation="Normalized from the source reaction title/caption and depicted transformation.")
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"


def add_step(r, ev, *, temp=None, time=None, atmosphere=None, operation="reaction hold", stage="reaction", participants=None):
    s = st.new_step(f"step_{len(r['steps'])+1}", ev, stage=stage)
    s["participant_ids"] = participants or [p["participant_id"] for p in r["participants"] if p["stage"] == stage]
    s["operation"] = F(operation, operation, ev)
    if temp is not None:
        s["temperature"] = F(temp if isinstance(temp, str) else f"{temp['number']} {temp['unit']}", temp, ev)
    if time is not None:
        s["time"] = F(time if isinstance(time, str) else f"{time['number']} {time['unit']}", time, ev)
    if atmosphere is not None:
        s["atmosphere"] = F(atmosphere, atmosphere, ev)
    r["steps"].append(s)


def add_outcome(r, ev, kind, target, value, *, basis=None, yield_type=None, scope=None, product_role="product"):
    o = st.new_outcome(f"o{len(r['outcomes'])+1}_{kind}", kind, [mid(target)] if target else [], ev)
    o["product_role"] = F(product_role, product_role, ev)
    raw = value.get("raw_text", f"{value['number']} {value['unit']}") if isinstance(value, dict) else value
    o["value"] = F(raw, value, ev)
    if yield_type:
        o["yield_type"] = F(yield_type, yield_type, ev)
    if basis:
        o["measurement_basis"] = F(basis, basis, ev)
    if scope is not None:
        o["scope"] = F(scope, scope, ev)
    r["outcomes"].append(o)


def classify(r, ev, codes):
    changes = []
    labels = {}
    pairs = {"C-C": ("C", "C", 0, 1, "formed"), "C-P": ("C", "P", 0, 1, "formed"),
             "DB-C-C": ("C", "C", 1, 2, "order_increased"), "C-S": ("C", "S", 0, 1, "formed"),
             "C-H": ("C", "H", 0, 1, "hydrogen_gain")}
    for code in codes:
        a, b, before, after, kind = pairs[code]
        cid = f"b{len(changes)+1}"
        change = {"change_id": cid, "atoms": [{"element": a}, {"element": b}],
                  "before_order": before, "after_order": after, "change_type": kind}
        if kind == "hydrogen_gain":
            change["hydrogen_count"] = 1
        changes.append(change)
        lab = st.bond_label(a, b, after)
        labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(cid)
    deriv = "Class-level net changes inferred from the source drawing/text; exact atom identities and full mapping are unresolved."
    r["classification"]["bond_changes"] = [F(None, c, ev, method="inference", derivation=deriv) for c in changes]
    r["classification"]["labels"] = F(None, list(labels.values()), ev, method="inference", derivation=deriv)
    r["classification"]["coverage"] = F(None, "partial", ev, method="inference", derivation=deriv)
    r["classification"]["limitations"] = ["Exact atom mapping, coproducts, proton source, and complete hydrogen accounting are unresolved."]
    r["classification"]["model_votes"] = []


def no_classification(r, ev):
    r["classification"]["bond_changes"] = []
    r["classification"]["labels"] = F(None, [], ev, method="inference", derivation="No observed transformation to classify.")
    r["classification"]["coverage"] = F(None, "not_applicable", ev, method="inference", derivation="Source reports 0%/inert outcome.")
    r["classification"]["limitations"] = ["No intended bond formation is assigned to an unsuccessful experiment."]


# Generic/literature molecule families from Scheme 1.
for label, name, fixed, vars_ in [
    ("lit_allenol_family", "generic allenic alcohol family", "allenic alcohol scaffold", ["R1", "R2"]),
    ("lit_rearrangement_product", "generic 2-phosphinoyl-1,3-butadiene family", "phosphinoyl diene scaffold", ["R1", "R2"]),
    ("lit_propargyl_alcohol", "generic terminal propargyl alcohol family", "terminal propargyl alcohol scaffold", ["R"]),
    ("lit_addition_products", "Markovnikov/anti-Markovnikov product mixture family", "regioisomeric phosphinoyl diene products", ["R", "P"]),
    ("lit_allenyl_phosphine_oxide", "generic allenylphosphine oxide family", "allenyl-P(O)Ph2 scaffold", ["R1", "R2"]),
    ("lit_aryl_or_disulfide", "arylboronic-acid or diaryl-disulfide partner family", "alternative partners in Scheme 1c", ["Ar"]),
    ("generic_allenol", "generic unprotected multisubstituted allenic alcohol", "unprotected allenic alcohol scaffold", ["R1", "R2", "R3"]),
    ("generic_POH", "generic P(O)H compound", "H-phosphonate or H-phosphinate scaffold", ["R", "Rprime"]),
    ("generic_product", "generic multisubstituted 2-phosphinoyl-1,3-butadiene", "2-phosphinoyl-1,3-butadiene scaffold", ["R1", "R2", "R3", "P"]),
]:
    add_template(label, name, 1, fixed, vars_)

# Concrete recurrent reagents/solvents.
add_molecule("1a", "1-phenylbuta-2,3-dien-1-ol", 2, smiles="C=C=CC(O)c1ccccc1", chirality="unspecified", geometry="unspecified", representation="partial", axial="terminal-CH2 allene axis unassigned; local axis is not stereogenic on the candidate graph")
add_molecule("1u", "buta-2,3-dien-1-ol", 2, smiles="C=C=CCO", chirality="achiral", geometry="not_applicable", representation="partial", axial="terminal-CH2 allene; no configurationally stable axis on the candidate graph")
add_molecule("2a", "diisopropyl phosphonate", 2, smiles="[PH](=O)(OC(C)C)OC(C)C", chirality="achiral", geometry="not_applicable")
add_molecule("2b", "dimethyl phosphonate", 3, smiles="[PH](=O)(OC)OC", chirality="achiral", geometry="not_applicable")
add_molecule("2c", "diethyl phosphonate", 3, smiles="[PH](=O)(OCC)OCC", chirality="achiral", geometry="not_applicable")
add_molecule("2d", "dibutyl phosphonate", 3, smiles="[PH](=O)(OCCCC)OCCCC", chirality="achiral", geometry="not_applicable")
add_molecule("2e", "diphenyl phosphonate", 3, smiles="[PH](=O)(Oc1ccccc1)Oc1ccccc1", chirality="achiral", geometry="not_applicable")
add_molecule("2f", "ethyl phenylphosphinate", 3, smiles="[PH](=O)(OCC)c1ccccc1", chirality="unspecified", geometry="not_applicable", representation="partial")
add_molecule("2g", "diphenylphosphine oxide", 3, smiles="[PH](=O)(c1ccccc1)c1ccccc1", chirality="achiral", geometry="not_applicable")
add_molecule("2h", "dienophile 2h as drawn (diethyl alkene dicarboxylate)", 3, chirality="unspecified", geometry="unspecified", note="The compact drawing is retained without choosing E/Z or a unique name.")
add_molecule("2i", "4-methoxythiophenol", 3, smiles="COc1ccc(S)cc1", chirality="achiral", geometry="not_applicable")
add_molecule("pd_acac2", "Pd(acac)2", 2, chirality="not_applicable", geometry="not_applicable")
add_molecule("xantphos", "xantphos", 2, chirality="unspecified", geometry="not_applicable")
add_molecule("h3po3", "phosphorous acid", 2, smiles="[PH](=O)(O)O", chirality="achiral", geometry="not_applicable")
add_molecule("thf", "tetrahydrofuran", 2, smiles="C1CCOC1", chirality="achiral", geometry="not_applicable")
add_molecule("toluene", "toluene", 3, smiles="Cc1ccccc1", chirality="achiral", geometry="not_applicable")
add_molecule("water", "water", 3, smiles="O", chirality="achiral", geometry="not_applicable")
add_molecule("ethanol", "ethanol", 3, smiles="CCO", chirality="achiral", geometry="not_applicable")
add_molecule("tBuONa", "sodium tert-butoxide", 3, smiles="[Na+].[O-]C(C)(C)C", chirality="achiral", geometry="not_applicable")
add_molecule("ph2pcl", "chlorodiphenylphosphine", 1, chirality="unspecified", geometry="not_applicable")
add_molecule("et3n", "triethylamine", 1, smiles="CCN(CC)CC", chirality="achiral", geometry="not_applicable")
add_molecule("ni_cat", "nickel catalyst (unspecified)", 1, chirality="not_applicable", geometry="not_applicable")
add_molecule("pd0", "palladium(0) catalytic species", 4, chirality="not_applicable", geometry="not_applicable")

# Table 1 substrate/product identities.
t1 = [
    ("a", "H", 92, None), ("b", "4-F", 86, None), ("c", "4-Cl", 82, None),
    ("d", "3-Cl", 80, None), ("e", "2-Cl", 60, None), ("f", "4-CF3", 69, None),
    ("g", "4-CN", 56, None), ("h", "4-NO2", 67, None), ("i", "4-Me", 81, None),
    ("j", "4-SMe", 72, None), ("k", "4-CO2Me", 82, None), ("l", "4-OMe", 75, None),
    ("m", "4-tBu", 83, None), ("n", "4-Ph", 77, None), ("o", "2-naphthyl", 84, None),
    ("p", "2-furyl", 67, None), ("q", "benzothiophenyl", 76, None), ("r", "2-pyridyl", 75, None),
    ("s", "n-pentyl", 85, ("E/Z", [72, 28])), ("t", "phenylethynyl", 80, ("E/Z", [41, 59])),
    ("u", "H (unsubstituted buta-2,3-dien-1-ol)", 72, None),
]
for suffix, desc, _, ratio in t1:
    sl = f"1{suffix}"
    if sl not in molecules:
        add_molecule(sl, f"1-substituted allenic alcohol {sl} ({desc})", 2,
                     chirality="unspecified", geometry="unspecified", representation="not_encoded",
                     axial="allene axial configuration not assigned in the source")
    pl = f"3{suffix}"
    if suffix == "a":
        add_molecule(pl, "(E)-3a, diisopropyl 2-(1-phenylbuta-1,3-dien-2-yl)phosphonate", 2,
                     smiles="C=C(P(=O)(OC(C)C)OC(C)C)/C=C/c1ccccc1", chirality="achiral",
                     geometry="specified", assignment="E major/isolate; trace Z not isolated")
    elif suffix == "u":
        add_molecule(pl, "diisopropyl buta-1,3-dien-2-ylphosphonate (3u)", 2,
                     smiles="C=C(P(=O)(OC(C)C)OC(C)C)C=C", chirality="achiral", geometry="not_applicable")
    elif ratio:
        add_molecule(pl, f"2-phosphinoyl-1,3-butadiene {pl} ({desc}) E/Z mixture", 2,
                     chirality="unspecified", geometry="mixture", representation="not_encoded",
                     mixture={"metric": ratio[0], "ratio": ratio[1], "order": ["E", "Z"]},
                     assignment=f"{ratio[0]}={ratio[1][0]}:{ratio[1][1]}")
    else:
        add_molecule(pl, f"(E)-2-phosphinoyl-1,3-butadiene {pl} ({desc})", 2,
                     chirality="unspecified", geometry="specified", representation="not_encoded",
                     assignment="E stereoisomer reported as highly selective; trace Z not isolated")

# Table 2 identities and ratios.
t2 = [
    ("v", "C6H13-substituted tertiary allenol", 86, "E/Z", [50, 50], True),
    ("w", "cyclohexyl-substituted tertiary allenol", 79, None, None, True),
    ("x", "tertiary aryl allenol, R=H", 82, "E/Z", [83, 17], True),
    ("y", "tertiary aryl allenol, R=Cl", 71, "E/Z", [85, 15], True),
    ("z", "tertiary aryl allenol, R=C6H5", 89, "E/Z", [90, 10], True),
    ("aa", "4-chlorophenyl-substituted tertiary allenol", 78, "E/Z", [87, 13], True),
    ("ab", "1-aryl-4-alkyl allenol, R=H", 94, "E,E/E,Z", [90, 10], False),
    ("ac", "1-aryl-4-alkyl allenol, R=4-F", 90, "E,E/E,Z", [86, 14], False),
    ("ad", "1-aryl-4-alkyl allenol, R=4-CF3", 89, "E,E/E,Z", [85, 15], False),
    ("ae", "1-aryl-4-alkyl allenol, R=4-Me", 91, "E,E/E,Z", [90, 10], False),
    ("af", "alkyl-substituted allenic alcohol", 81, None, None, False),
    ("ag", "alkyl-substituted allenic alcohol", 89, "E,E/E,Z", [50, 50], False),
    ("ah", "C4-phenyl allenic alcohol", 41, None, None, False),
    ("ai", "C4-phenyl alkyl allenic alcohol", 51, "E,Z/Z,Z", [50, 50], False),
    ("aj", "1,1,4-trisubstituted allenic alcohol", 74, "E,Z/E,E/Z,Z", [23, 70, 7], True),
]
for suffix, desc, _, metric, ratio, _ in t2:
    add_molecule(f"1{suffix}", f"{desc} (1{suffix})", 3, chirality="unspecified", geometry="unspecified",
                 representation="not_encoded", axial="allene axis/configuration not assigned in the supplied pages")
    if ratio:
        add_molecule(f"3{suffix}", f"multisubstituted 2-phosphinoyl-1,3-butadiene 3{suffix} mixture", 3,
                     chirality="unspecified", geometry="mixture", representation="not_encoded",
                     mixture={"metric": metric, "ratio": ratio, "order": metric.split("/")},
                     assignment=f"{metric}=" + ":".join(map(str, ratio)))
    else:
        add_molecule(f"3{suffix}", f"multisubstituted 2-phosphinoyl-1,3-butadiene 3{suffix}", 3,
                     chirality="unspecified", geometry="unspecified", representation="not_encoded")
    if suffix == "aa":
        molecules[f"3{suffix}"]["stereochemistry"]["assignment"] = F(
            "X-ray of E-3aa, CCDC 2034498", "E assignment supported by the reported X-ray structure of 3aa", E[3])

# Table 3 and Scheme 2 products.
for suffix, name, geom in [
    ("ak", "dimethyl-phosphonate diene 3ak", "unspecified"),
    ("al", "diethyl-phosphonate diene 3al", "unspecified"),
    ("am", "dibutyl-phosphonate diene 3am", "unspecified"),
    ("an", "diphenyl-phosphonate diene 3an", "unspecified"),
    ("ao", "ethyl phenylphosphinate diene 3ao", "unspecified"),
    ("ap", "desired diphenylphosphine-oxide diene 3ap", "unspecified"),
]:
    add_molecule(f"3{suffix}", name, 3, chirality="unspecified", geometry=geom, representation="not_encoded")
for label, name, chir, geom, mix in [
    ("4a", "Diels–Alder adduct 4a", "unspecified", "unspecified", None),
    ("4b", "Diels–Alder adduct 4b", "unspecified", "unspecified", None),
    ("5a", "Michael addition product 5a", "unspecified", "mixture", {"metric": "E/Z", "ratio": [71, 29], "order": ["E", "Z"]}),
    ("5b", "Michael addition product 5b", "unspecified", "unspecified", None),
]:
    add_molecule(label, name, 3, chirality=chir, geometry=geom, representation="not_encoded", mixture=mix,
                 assignment="E/Z=71:29" if label == "5a" else None)


# Scheme 1 literature examples and current generic scheme.
def literature(label, name, participants, classification_codes=None):
    ev = E[1]
    r = st.new_reaction(rid(label), DOC, ev)
    r["record_type"] = "literature_example"
    r["experiment_categories"] = F(["other"], ["other"], ev)
    r["location"] = F(name, {"page": 1, "location": name}, ev)
    r["reaction_name"] = F(name, name, ev)
    r["transformation_type"] = F(status="unreported")
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
                   "unresolved": ["Cited prior-work family; exact substrates, yields, and procedures are not reported here."]}
    for lab, role in participants:
        add_participant(r, lab, role, ev)
    if classification_codes:
        classify(r, ev, classification_codes)
    else:
        r["classification"]["coverage"] = F(None, "unresolved", ev, method="inference", derivation="Alternative generic pathways are grouped in the source panel.")
        r["classification"]["labels"] = F(status="unreported")
    reactions.append(r)


literature("scheme1a_prior_sigmatropic", "Scheme 1a prior [2,3]-sigmatropic rearrangement", [("lit_allenol_family", "reactant"), ("ph2pcl", "reactant"), ("et3n", "reagent"), ("lit_rearrangement_product", "product")], ["C-P", "C-C", "DB-C-C"])
literature("scheme1b_prior_ni_addition", "Scheme 1b prior Ni-catalyzed P(O)H addition", [("lit_propargyl_alcohol", "reactant"), ("generic_POH", "reactant"), ("ni_cat", "catalyst"), ("lit_addition_products", "product")], ["C-P", "C-H"])
literature("scheme1c_prior_allenyl_phosphine_oxide", "Scheme 1c prior allenylphosphine-oxide coupling families", [("lit_allenyl_phosphine_oxide", "reactant"), ("lit_aryl_or_disulfide", "reactant"), ("lit_rearrangement_product", "product")])

ev = E[1]
g = st.new_reaction(rid("scheme1d_generic"), DOC, ev)
g["record_type"] = "generic_scheme"
g["experiment_categories"] = F(["other"], ["other"], ev)
g["location"] = F("Scheme 1d, this work", {"page": 1, "location": "Scheme 1d"}, ev)
g["reaction_name"] = F("Pd-catalyzed coupling of unprotected allenic alcohols with P(O)H compounds", "Pd-catalyzed coupling of unprotected allenic alcohols with P(O)H compounds", ev)
g["transformation_type"] = F(None, "C–P cross-coupling/rearrangement to 2-phosphinoyl-1,3-butadienes", ev, method="inference", derivation="Normalized from title and Scheme 1d.")
for lab, role in [("generic_allenol", "reactant"), ("generic_POH", "reactant"), ("generic_product", "product")]:
    add_participant(g, lab, role, ev)
classify(g, ev, ["C-C", "C-P", "DB-C-C"])
g["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [], "unresolved": ["R-group family only; no single concrete experiment."]}
reactions.append(g)


def add_coupling_participants(r, page, substrate, p_reagent, additive_pct, ev):
    substrate_pid = add_participant(r, substrate, "reactant", ev, amount=Q(0.33, "mmol"))
    p_pid = add_participant(r, p_reagent, "reactant", ev, amount=Q(0.30, "mmol"))
    add_participant(r, "pd_acac2", "catalyst", ev, molpct=Q(2, "mol%"), reference=p_reagent)
    add_participant(r, "xantphos", "ligand", ev, molpct=Q(3, "mol%"), reference=p_reagent)
    add_participant(r, "h3po3", "additive", ev, molpct=Q(additive_pct, "mol%"), reference=p_reagent)
    add_participant(r, "thf", "solvent", ev, volume=Q(4.5, "mL"))
    return substrate_pid, p_pid


# Table 1: 21 concrete experiments.
for suffix, desc, yld, ratio in t1:
    ev = E[2]
    r = st.new_reaction(rid(f"table1_3{suffix}"), DOC, ev)
    cats = ["substrate_scope"] + (["optimization"] if suffix == "a" else [])
    set_common(r, 2, cats, f"Table 1, product 3{suffix}", "Palladium-catalyzed coupling of 1-substituted allenic alcohol with 2a", "allene-alcohol/P(O)H C–P cross-coupling")
    r["general_scheme_ids"] = [rid("scheme1d_generic")]
    add_coupling_participants(r, 2, f"1{suffix}", "2a", 20, ev)
    add_participant(r, f"3{suffix}", "product", ev, product_role="isolated coupling product")
    r["reaction_scale"] = F("2a 0.30 mmol", {"reference_molecule_id": mid("2a"), "amount": Q(0.30, "mmol")}, ev)
    r["solvent_conditions"]["composition"] = F("THF", [mid("thf")], ev)
    r["solvent_conditions"]["total_volume"] = F("4.5 mL", Q(4.5, "mL"), ev)
    add_step(r, ev, temp=Q(80, "°C"), time=Q(0.5, "h"), atmosphere="argon")
    add_outcome(r, ev, "yield", f"3{suffix}", Q(yld, "%"), basis="isolated", yield_type="isolated yield")
    if suffix <= "r":
        add_outcome(r, ev, "selectivity", f"3{suffix}", "high E selectivity; trace Z not isolated",
                    basis="in situ 31P NMR and 1H NMR", scope={"metric": "E/Z", "major": "E", "ratio": "qualitative"})
    if ratio:
        metric, vals = ratio
        add_outcome(r, ev, "selectivity", f"3{suffix}", f"{metric} = {vals[0]}:{vals[1]}",
                    basis="31P NMR and 1H NMR", scope={"metric": metric, "order": metric.split("/"), "ratio": vals})
    classify(r, ev, ["C-C", "C-P", "DB-C-C"])
    r["procedure"]["general_procedure_reference"] = F("Table 1 footnote a; ESI Table S1 cited but not supplied", "Table 1 footnote a", ev)
    r["review"]["unresolved"] = ["Workup/purification details and full procedure are unavailable because the advertised ESI was not supplied.", "Exact atom mapping and complete product-structure encoding remain unresolved."]
    if suffix == "a":
        r["notes"].append(F("The prose optimized model experiment and Table 1 product 3a have the same reported conditions and are represented once.", "Cross-page duplicate evidence merged into this record.", ev, method="inference", derivation="Matching substrates, conditions, and 92% Table 1 outcome."))
    reactions.append(r)


# Table 2: 15 concrete experiments.
for suffix, desc, yld, metric, ratio, harsh in t2:
    ev = E[3]
    r = st.new_reaction(rid(f"table2_3{suffix}"), DOC, ev)
    set_common(r, 3, ["substrate_scope"], f"Table 2, product 3{suffix}", "Pd-catalyzed coupling of multisubstituted allenic alcohol with 2a", "multisubstituted allene-alcohol/P(O)H C–P cross-coupling")
    r["general_scheme_ids"] = [rid("scheme1d_generic")]
    add_coupling_participants(r, 3, f"1{suffix}", "2a", 60 if harsh else 20, ev)
    add_participant(r, f"3{suffix}", "product", ev, product_role="isolated coupling product")
    r["reaction_scale"] = F("2a 0.30 mmol", {"reference_molecule_id": mid("2a"), "amount": Q(0.30, "mmol")}, ev)
    r["solvent_conditions"]["composition"] = F("THF", [mid("thf")], ev)
    r["solvent_conditions"]["total_volume"] = F("4.5 mL", Q(4.5, "mL"), ev)
    add_step(r, ev, temp=Q(90 if harsh else 80, "°C"), time=Q(10, "h"), atmosphere="argon")
    add_outcome(r, ev, "yield", f"3{suffix}", Q(yld, "%"), basis="isolated", yield_type="isolated yield")
    if ratio:
        add_outcome(r, ev, "selectivity", f"3{suffix}", f"{metric} = " + ":".join(map(str, ratio)),
                    basis="31P NMR and 1H NMR", scope={"metric": metric, "order": metric.split("/"), "ratio": ratio})
    classify(r, ev, ["C-C", "C-P", "DB-C-C"])
    r["procedure"]["general_procedure_reference"] = F("Table 2 footnote a" + (" with footnote d modification" if harsh else ""), "Table 2 footnote conditions", ev)
    r["review"]["unresolved"] = ["Workup/purification and individual substrate preparation details are unavailable without ESI.", "Exact atom mapping and safe complete structures are unresolved."]
    reactions.append(r)


# Table 3: six P(O)H-scope experiments.
t3 = [
    ("ak", "1a", "2b", 70, False), ("al", "1a", "2c", 85, False),
    ("am", "1a", "2d", 90, False), ("an", "1w", "2e", 43, True),
    ("ao", "1w", "2f", 72, True), ("ap", "1a", "2g", 0, False),
]
for suffix, substrate, preg, yld, harsh in t3:
    ev = E[3]
    r = st.new_reaction(rid(f"table3_3{suffix}"), DOC, ev)
    set_common(r, 3, ["substrate_scope"], f"Table 3, product 3{suffix}", "Substrate scope of P(O)H compounds", "allene-alcohol/P(O)H C–P cross-coupling")
    r["general_scheme_ids"] = [rid("scheme1d_generic")]
    add_coupling_participants(r, 3, substrate, preg, 60 if harsh else 20, ev)
    product = f"3{suffix}"
    add_participant(r, product, "product", ev, product_role="desired coupling product" if yld == 0 else "isolated coupling product")
    r["reaction_scale"] = F(f"{preg} 0.30 mmol", {"reference_molecule_id": mid(preg), "amount": Q(0.30, "mmol")}, ev)
    r["solvent_conditions"]["composition"] = F("THF", [mid("thf")], ev)
    r["solvent_conditions"]["total_volume"] = F("4.5 mL", Q(4.5, "mL"), ev)
    add_step(r, ev, temp=Q(90 if harsh else 80, "°C"), time=Q(10 if harsh else 0.5, "h"), atmosphere="argon")
    add_outcome(r, ev, "yield", product, Q(yld, "%"), basis="isolated", yield_type="isolated yield")
    if yld:
        classify(r, ev, ["C-C", "C-P", "DB-C-C"])
    else:
        no_classification(r, ev)
        r["notes"].append(F("Secondary phosphine oxide 2g was inert under similar reaction conditions.", "2g inert; desired 3ap reported at 0%.", ev))
    r["review"]["unresolved"] = ["Workup/purification and complete structures are unavailable without ESI."]
    reactions.append(r)


# Scheme 2 application experiments.
def scheme2_base(label, title):
    ev = E[3]
    r = st.new_reaction(rid(label), DOC, ev)
    set_common(r, 3, [], f"Scheme 2, {title}", title, "synthetic application of phosphinoyl diene")
    return r, ev


r, ev = scheme2_base("scheme2a_scaleup_3u", "gram-scale preparation of 3u")
r["experiment_categories"] = F(["scale_up"], ["scale_up"], ev)
add_participant(r, "1u", "reactant", ev, mass=Q(1.54, "g"))
add_participant(r, "2a", "reactant", ev, mass=Q(3.32, "g"))
add_participant(r, "pd_acac2", "catalyst", ev, molpct=Q(2, "mol%"), reference="2a")
add_participant(r, "xantphos", "ligand", ev, molpct=Q(3, "mol%"), reference="2a")
add_participant(r, "h3po3", "additive", ev, molpct=Q(20, "mol%"), reference="2a")
add_participant(r, "thf", "solvent", ev)
add_participant(r, "3u", "product", ev, mass=Q(2.84, "g"), product_role="scale-up product")
r["reaction_scale"] = F("20 mmol of 2a (prose); 2a 3.32 g (scheme)", {"reference_molecule_id": mid("2a"), "amount": Q(20, "mmol")}, ev)
r["solvent_conditions"]["composition"] = F("THF", [mid("thf")], ev)
add_step(r, ev, temp=Q(80, "°C"), time=Q(2, "h"), atmosphere="argon")
add_outcome(r, ev, "product_amount", "3u", Q(2.84, "g"), scope={"stage": "isolated product amount"})
add_outcome(r, ev, "yield", "3u", Q(65, "%"))
classify(r, ev, ["C-C", "C-P", "DB-C-C"])
r["review"]["unresolved"] = ["Yield measurement basis and THF volume are not stated; other details are only described as similar conditions.", "Workup/purification is unavailable without ESI."]
reactions.append(r)

r, ev = scheme2_base("scheme2b_diels_alder_4a", "[4+2] cycloaddition of 3u with 2h")
r["experiment_categories"] = F(["derivatization"], ["derivatization"], ev)
add_participant(r, "3u", "reactant", ev, amount=Q(0.3, "mmol"))
add_participant(r, "2h", "reactant", ev, amount=Q(0.6, "mmol"))
add_participant(r, "toluene", "solvent", ev)
add_participant(r, "4a", "product", ev, product_role="cycloaddition product")
r["solvent_conditions"]["composition"] = F("toluene", [mid("toluene")], ev)
add_step(r, ev, temp=Q(120, "°C"), time=Q(4, "h"))
add_outcome(r, ev, "yield", "4a", Q(72, "%"))
classify(r, ev, ["C-C", "C-C"])
r["review"]["unresolved"] = ["Yield measurement basis, atmosphere, reagent stereochemistry, workup, and purification are not stated."]
reactions.append(r)

r, ev = scheme2_base("scheme2b_diels_alder_4b", "formal [4+2] cycloaddition/dimerization of 3a")
r["experiment_categories"] = F(["derivatization"], ["derivatization"], ev)
add_participant(r, "3a", "reactant", ev, amount=Q(0.4, "mmol"))
add_participant(r, "thf", "solvent", ev)
add_participant(r, "4b", "product", ev, product_role="cycloaddition product")
r["solvent_conditions"]["composition"] = F("THF", [mid("thf")], ev)
add_step(r, ev, temp=Q(80, "°C"), time=Q(12, "h"))
add_outcome(r, ev, "yield", "4b", Q(65, "%"))
classify(r, ev, ["C-C", "C-C"])
r["review"]["unresolved"] = ["Yield measurement basis, stoichiometric interpretation, atmosphere, workup, and purification are not stated."]
reactions.append(r)

r, ev = scheme2_base("scheme2c_michael_5a", "Michael addition of 4-methoxythiophenol to 3u")
r["experiment_categories"] = F(["derivatization"], ["derivatization"], ev)
add_participant(r, "3u", "reactant", ev, amount=Q(0.6, "mmol"))
add_participant(r, "2i", "reactant", ev, amount=Q(1.2, "mmol"))
add_participant(r, "water", "solvent", ev)
add_participant(r, "5a", "product", ev, product_role="Michael addition product")
r["solvent_conditions"]["composition"] = F("H2O", [mid("water")], ev)
add_step(r, ev, temp="rt", time=Q(6, "h"))
add_outcome(r, ev, "yield", "5a", Q(51, "%"))
add_outcome(r, ev, "selectivity", "5a", "E/Z = 71:29", scope={"metric": "E/Z", "order": ["E", "Z"], "ratio": [71, 29]})
classify(r, ev, ["C-S", "C-H"])
r["review"]["unresolved"] = ["Yield and E/Z measurement bases, solvent volume, atmosphere, workup, and purification are not stated."]
reactions.append(r)

r, ev = scheme2_base("scheme2c_michael_5b", "Michael addition of diphenylphosphine oxide to 3u")
r["experiment_categories"] = F(["derivatization"], ["derivatization"], ev)
add_participant(r, "3u", "reactant", ev, amount=Q(0.6, "mmol"))
add_participant(r, "2g", "reactant", ev, amount=Q(1.2, "mmol"))
add_participant(r, "tBuONa", "reagent", ev)
add_participant(r, "ethanol", "solvent", ev)
add_participant(r, "5b", "product", ev, product_role="Michael addition product")
r["solvent_conditions"]["composition"] = F("EtOH", [mid("ethanol")], ev)
add_step(r, ev, temp="rt", time=Q(1, "h"))
add_outcome(r, ev, "yield", "5b", Q(77, "%"))
classify(r, ev, ["C-P", "C-H"])
r["review"]["unresolved"] = ["Yield measurement basis, t-BuONa loading, solvent volume, atmosphere, workup, and purification are not stated."]
reactions.append(r)


# Scheme 3 proposed mechanism.
for label, name in [("mech_A", "hydrogen-bonded allenic alcohol species A"), ("mech_B", "allylpalladium hydroxide B"), ("mech_C", "palladium phosphonate intermediate C")]:
    add_molecule(label, name, 4, chirality="unresolved", geometry="unresolved", representation="not_encoded")
ev = E[4]
r = st.new_reaction(rid("scheme3_proposed_mechanism"), DOC, ev)
r["record_type"] = "proposed_mechanism"
r["experiment_categories"] = F(["mechanistic_study"], ["mechanistic_study"], ev)
r["location"] = F("Scheme 3", {"page": 4, "location": "Scheme 3"}, ev)
r["reaction_name"] = F("Plausible reaction mechanism", "Plausible reaction mechanism", ev)
r["transformation_type"] = F(None, "author-proposed Pd(0)/allylpalladium/transmetalation/reductive-elimination cycle", ev, method="inference", derivation="Normalized from Scheme 3 and adjacent prose.")
for lab, role in [("1a", "reactant"), ("h3po3", "additive"), ("pd0", "catalyst"), ("mech_A", "unknown"), ("mech_B", "unknown"), ("2a", "reactant"), ("mech_C", "unknown"), ("3a", "product"), ("water", "product")]:
    add_participant(r, lab, role, ev, product_role="proposed product" if role == "product" else None)
classify(r, ev, ["C-C", "C-P", "DB-C-C"])
r["notes"].append(F("Authors propose hydrogen-bond activation, Pd(0) attack/C–O cleavage, transmetalation, and reductive elimination.", "Proposed mechanism; intermediates A–C are not treated as isolated experiments.", ev))
r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
               "unresolved": ["Mechanistic intermediates are proposed, not independently isolated or quantified on the supplied pages."]}
reactions.append(r)


bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": rid("context_metadata_esi"), "document_id": DOC, "kind": "metadata_and_source_limit",
     "evidence_ids": E[1], "value": "Title/DOI/authors and ESI availability notice. ESI and CCDC 2034498 data are advertised but not supplied or associated."},
    {"context_id": rid("context_xray_3aa"), "document_id": DOC, "kind": "structure_assignment",
     "evidence_ids": E[3], "value": "X-ray structure of E-3aa and CCDC 2034498 are reported; no CIF is supplied."},
    {"context_id": rid("context_spectra_inventory"), "document_id": DOC, "kind": "negative_inventory",
     "evidence_ids": E[1] + E[2] + E[3] + E[4], "value": "No standalone mass-spectrum peak table, molecular-weight table, or fragment-relative-abundance block occurs on pages 1–4."},
    {"context_id": rid("context_references"), "document_id": DOC, "kind": "bibliography",
     "evidence_ids": E[4], "value": "Page 4 bibliography is not converted into current-paper experiments beyond the three explicitly depicted Scheme 1 precedents."},
]

(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")

raw_inventory = {
    "document_id": DOC,
    "pages": [
        {"page": 1, "image_evidence_id": EI[1][0], "inventory": ["title/metadata/abstract", "Scheme 1a–c: three cited prior-work arrows/families", "Scheme 1d: one current generic scheme", "ESI/CCDC availability notice"], "record_ids": [rid("scheme1a_prior_sigmatropic"), rid("scheme1b_prior_ni_addition"), rid("scheme1c_prior_allenyl_phosphine_oxide"), rid("scheme1d_generic")]},
        {"page": 2, "image_evidence_id": EI[2][0], "inventory": ["Table 1 header and footnotes", "21 product entries 3a–3u", "prose high-E assignment for 3a–3r", "optimized-model prose merged with 3a"], "record_ids": [rid(f"table1_3{x[0]}") for x in t1]},
        {"page": 3, "image_evidence_id": EI[3][0], "inventory": ["Table 2: 15 entries 3v–3aj", "Table 3: six entries 3ak–3ap", "Scheme 2: scale-up plus four derivatizations", "X-ray E-3aa context"], "record_ids": [rid(f"table2_3{x[0]}") for x in t2] + [rid(f"table3_3{x[0]}") for x in t3] + [rid(x) for x in ["scheme2a_scaleup_3u", "scheme2b_diels_alder_4a", "scheme2b_diels_alder_4b", "scheme2c_michael_5a", "scheme2c_michael_5b"]]},
        {"page": 4, "image_evidence_id": EI[4][0], "inventory": ["Scheme 3 proposed mechanism", "conclusion/acknowledgements/conflict statement", "references"], "record_ids": [rid("scheme3_proposed_mechanism")]},
    ],
    "deduplication": ["Table 1 product 3a is the optimized model reaction described in prose and is counted once.", "The scale-up 3u is a different experiment from Table 1 3u because scale/time/yield differ.", "Scheme 3 intermediates are one proposed mechanism, not executed experiments."],
}
(RUN / "raw_model/page_extraction.json").write_text(json.dumps(raw_inventory, ensure_ascii=False, indent=2) + "\n")

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title/abstract and ESI footnote", "finding": "Paper metadata and an ESI/CCDC availability notice; no supplied SI is associated.", "evidence_ids": EI[1]},
            {"region": "Scheme 1a–c", "finding": "Three explicitly cited prior-work transformation families, retained as literature examples.", "evidence_ids": EI[1]},
            {"region": "Scheme 1d", "finding": "Current-work R-group generic Pd-catalyzed coupling scheme, not one concrete experiment.", "evidence_ids": EI[1]}],
         "record_ids": raw_inventory["pages"][0]["record_ids"], "unresolved": ["Prior-work exact conditions/yields and advertised ESI are unavailable."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "left-column model-reaction and scope prose", "finding": "Optimized conditions and qualitative high-E assignment; matched to Table 1 3a and product records.", "evidence_ids": EI[2]},
            {"region": "Table 1 header/body", "finding": "All 21 entries 3a–3u inventoried; yields and two numerical E/Z ratios transcribed.", "evidence_ids": EI[2]},
            {"region": "Table 1 footnotes", "finding": "Common charges, catalyst/ligand/additive, THF volume, 80 °C/0.5 h/argon, isolated-yield basis, and NMR E/Z method extracted.", "evidence_ids": EI[2]}],
         "record_ids": raw_inventory["pages"][1]["record_ids"], "unresolved": ["No SI optimization rows, procedures, workup, or purification were supplied; compact product structures remain machine candidates/unencoded."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Table 2", "finding": "All 15 multisubstituted-allenol entries 3v–3aj and footnote-d condition changes inventoried.", "evidence_ids": EI[3]},
            {"region": "Table 3", "finding": "All six P(O)H-scope entries 3ak–3ap, including inert 2g/0% 3ap, inventoried.", "evidence_ids": EI[3]},
            {"region": "Scheme 2", "finding": "One gram-scale preparation and four separate derivatization arrows recorded; all printed charges/conditions/results retained.", "evidence_ids": EI[3]},
            {"region": "X-ray and adjacent prose", "finding": "E-3aa/CCDC assignment retained as structure context; external CIF not supplied.", "evidence_ids": EI[3]}],
         "record_ids": raw_inventory["pages"][2]["record_ids"], "unresolved": ["Scheme 2 yield/selectivity methods, several reagent loadings, solvent volumes, atmospheres, workups, and purification are unreported."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "Scheme 3", "finding": "One author-proposed catalytic mechanism; A–C are not expanded into executed experiments.", "evidence_ids": EI[4]},
            {"region": "conclusion/acknowledgements/conflicts", "finding": "Non-experimental summary and funding/conflict statements.", "evidence_ids": EI[4]},
            {"region": "references", "finding": "Bibliography only; no current-paper reaction arrow, table row, procedure, or standalone mass-spectrum table.", "evidence_ids": EI[4]}],
         "record_ids": raw_inventory["pages"][3]["record_ids"], "unresolved": ["Proposed intermediates are not isolated or quantified; external literature and CCDC files are not supplied."]},
    ],
}
(RUN / "coverage").mkdir(exist_ok=True)
(RUN / "coverage/pages_0001_0004.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")

stereo_checks = {}
for label, molecule in molecules.items():
    smiles = molecule["smiles"]["normalized_value"]
    if smiles:
        stereo_checks[mid(label)] = st.stereo_inspect(smiles)
(RUN / "raw_model/stereo_checks.json").write_text(json.dumps(stereo_checks, ensure_ascii=False, indent=2) + "\n")

print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]), "source_context": len(bundle["source_context"])}, ensure_ascii=False))
