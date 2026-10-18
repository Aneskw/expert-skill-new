#!/usr/bin/env python3
"""Build the source-derived v3.1 bundle for Org. Lett. 2020, 22, 6750–6755.

The script is intentionally document-specific.  It records only facts visible in
the supplied six-page main article; no Supporting Information is associated.
"""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st  # noqa: E402


RUN = ROOT / "runs/current"
DATASET_PATH = RUN / "dataset.json"
DOC = "doc_c35696f4709ebd07"
P = {n: f"{DOC}:page{n}:image" for n in range(1, 7)}
T = {n: f"{DOC}:page{n}:text" for n in range(1, 7)}


def f(raw=None, normalized=None, ev=(), method="original_extraction", derivation=None,
      status=None):
    return st.fact(raw, normalized, ev, method=method, derivation=derivation,
                   status=status)


def inferred(normalized, ev, why):
    return f(None, normalized, ev, method="inference", derivation=why)


def derived(normalized, ev, why, raw=None):
    return f(raw, normalized, ev, method="structure_derivation", derivation=why)


def na():
    return f(method="not_applicable", status="not_applicable")


def q(number, unit, **extra):
    return {"number": number, "unit": unit, **extra}


dataset = json.loads(DATASET_PATH.read_text())
doc = dataset["documents"][0]
doc["title"] = "Construction of CF3-Containing Tetrahydropyrano[3,2-b]indoles through DMAP-Catalyzed [4+1]/[3+3] Domino Sequential Annulation"
doc["supplementary_material"] = {
    "status": "unknown",
    "note": "The main article advertises Supporting Information, but no SI file was supplied or explicitly associated with this extraction."
}
st.write(DATASET_PATH, dataset)
bundle = st.empty(copy.deepcopy(dataset["documents"]))


MOLECULES = {}
STEREO_CHECKS = {}


def molecule(mid, label, name=None, ev=(P[2],), smiles=None, stereo_kind="unresolved",
             note=None, axial_note=None):
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    m["names"]["original"] = f(name or label, name or label, ev)
    m["names"]["english"] = f(name, name, ev) if name else f()
    m["cas_numbers"] = f()
    m["identity_note"] = f(note, note, ev) if note else f()
    if smiles:
        diag = st.stereo_inspect(smiles)
        canonical = diag["canonical_isomeric_smiles"]
        STEREO_CHECKS[mid] = {"input_smiles": smiles, **diag}
        m["representation_type"] = "smiles"
        m["smiles"] = derived(canonical, ev,
            "Candidate connectivity transcribed from the supplied drawing/name and canonicalized by the standard stereo graph diagnostic; parsing is not source identity certification.", raw=smiles)
        m["stereochemistry"]["atom_centers"] = derived(diag["atom_centers"], ev,
            "RDKit standard stereo diagnostic on the candidate SMILES.")
        m["stereochemistry"]["double_bonds"] = derived(diag["double_bonds"], ev,
            "RDKit standard stereo diagnostic on the candidate SMILES.")
        m["stereochemistry"]["axial"] = derived([], ev,
            "No axial descriptor is encoded in the candidate SMILES; allene axes still require explicit source review.")
        m["stereochemistry"]["mixture"] = f()
        m["stereochemistry"]["assignment"] = derived(
            {"candidate_graph_only": True, "source_absolute_assignment": "not reported"}, ev,
            "Graph diagnostic only; it does not establish sample stereochemical identity.")
        has_unspecified = any(x["specified"] == "Unspecified" for x in diag["potential_stereo"])
        m["stereochemistry"]["chirality_status"] = derived("unresolved" if has_unspecified else "achiral", ev,
            "The standard stereo diagnostic exposes an unspecified potential stereoelement; no source assignment is invented." if has_unspecified else
            "No potential atom stereocenter or stereogenic allene axis is present in this candidate graph.")
        m["stereochemistry"]["geometry_status"] = derived("unresolved" if has_unspecified else "not_applicable", ev,
            "Potential graph stereochemistry is unspecified and is not converted into a source E/Z assignment." if has_unspecified else
            "No independently assignable E/Z alkene is present in this candidate graph.")
        m["stereochemistry"]["representation_status"] = derived("partial" if has_unspecified else "encoded", ev,
            "The graph contains an unspecified potential stereoelement and therefore cannot claim fully encoded stereochemistry." if has_unspecified else
            "The candidate graph has no source-specified stereoelement omitted from its isomeric SMILES.")
        m["stereochemistry"]["isomeric_smiles"] = derived(canonical, ev,
            "Canonical isomeric SMILES returned by the standard stereo diagnostic.")
        m["stereochemistry"]["identity_key"] = derived(diag["stereo_graph_key"], ev,
            "Stereo graph key returned by the standard stereo diagnostic; not a sample identity key.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = f()
        if stereo_kind in ("product_relative", "l_menthyl", "product_l_menthyl"):
            raw = "Source drawing depicts relative stereochemical bonds" if stereo_kind != "l_menthyl" else "L-menthyl substituent"
            norm = {"relative_configuration": "depicted", "absolute_configuration": "unreported"}
            if stereo_kind in ("l_menthyl", "product_l_menthyl"):
                norm["named_chiral_substituent"] = "L-menthyl; center mapping unresolved"
            m["stereochemistry"]["atom_centers"] = f(raw, norm, ev)
            m["stereochemistry"]["assignment"] = f(raw, norm, ev)
            m["stereochemistry"]["chirality_status"] = f(raw, "partially_specified", ev)
            m["stereochemistry"]["representation_status"] = inferred(
                "not_encoded", ev, "The dense source structure/relative configuration was not converted to a complete isomeric SMILES.")
        else:
            m["stereochemistry"]["chirality_status"] = inferred(
                "unresolved", ev, "No complete candidate graph was encoded; absence of a source descriptor is not treated as achiral or racemic.")
            m["stereochemistry"]["representation_status"] = inferred(
                "unresolved", ev, "The complete source structure was not safely encoded.")
        m["stereochemistry"]["double_bonds"] = f()
        if axial_note:
            m["stereochemistry"]["axial"] = derived(axial_note, ev,
                "The drawing contains a terminal allene; one terminal carbon is CH2, so that allene axis is not stereogenic. This does not establish whole-molecule achirality.")
        else:
            m["stereochemistry"]["axial"] = f()
        m["stereochemistry"]["mixture"] = f()
        if stereo_kind not in ("product_relative", "l_menthyl", "product_l_menthyl"):
            m["stereochemistry"]["assignment"] = f()
        m["stereochemistry"]["geometry_status"] = inferred(
            "unresolved", ev, "Complete graph and all source geometric assignments were not established.")
        m["stereochemistry"]["isomeric_smiles"] = f()
        m["stereochemistry"]["identity_key"] = f()
    MOLECULES[mid] = m
    return mid


def template_molecule(mid, label, fixed, variables, ev):
    molecule(mid, label, label, ev=ev, note="Generic R-group structure; not a concrete compound identity.")
    m = MOLECULES[mid]
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": f(fixed, fixed, ev),
        "variables": [{
            "label": v,
            "attachment_points": f(f"source drawing attachment for {v}", {"label": v, "attachment": "as drawn"}, ev),
            "definition": f(f"{v} as shown in the source generic scheme", f"source-defined {v}", ev),
            "allowed_values": f(),
            "constraints": f()
        } for v in variables],
        "instance_mapping": [],
        "reason": f("R-group scheme", "The source uses unexpanded R/PG variables.", ev)
    }
    m["stereochemistry"]["chirality_status"] = inferred("unresolved", ev,
        "Variable substituents and product stereochemistry prevent a unique assignment.")
    m["stereochemistry"]["geometry_status"] = inferred("unresolved", ev,
        "Variable structure prevents a unique geometry assignment.")
    m["stereochemistry"]["representation_status"] = inferred("not_encoded", ev,
        "The source is a generic template, not a concrete isomeric SMILES.")
    MOLECULES[mid] = m
    return mid


# Core source compounds.  Only the two fully legible model substrates are given
# candidate structures; dense products and scope structures remain unresolved.
molecule("m_1a", "1a", "o-(tosylamino)trifluoroacetophenone derivative 1a", (P[1], P[2], P[3], P[4]),
         "Cc1ccc(S(=O)(=O)Nc2ccccc2C(=O)C(F)(F)F)cc1",
         note="Candidate structure from the model-substrate drawing; complete identity remains needs_review.")
molecule("m_2a", "2a", "benzyl 2-(acetoxymethyl)buta-2,3-dienoate", (P[1], P[2], P[3], P[4]),
         "C=C=C(COC(C)=O)C(=O)OCc1ccccc1",
         note="Terminal beta-prime-acetoxy allenoate candidate; the allene axis is non-stereogenic because one terminus is CH2.")

product_labels = [f"3a{x}" for x in "abcdefghijklm"] + [f"3{x}a" for x in "bcdefghijklmnop"]
for label in product_labels:
    ev = (P[3],) if label != "3aa" else (P[1], P[2], P[3], P[4])
    kind = "product_l_menthyl" if label == "3la" else "product_relative"
    note = None
    if label == "3ag":
        note = "Single-crystal X-ray analysis reported; CCDC 1951502. Absolute descriptors are not stated in the supplied main text."
    molecule(f"m_{label}", label, f"tetrahydropyrano[3,2-b]indole {label}", ev,
             stereo_kind=kind, note=note)

for label in [f"2{x}" for x in "bcdefghijklm"]:
    molecule(f"m_{label}", label, f"beta-prime-acetoxy allene scope substrate {label}", (P[3],),
             axial_note={"axis": "terminal allene", "stereogenic": False, "basis": "one terminal CH2"})
for label in [f"1{x}" for x in "bcdefghijklmnop"]:
    kind = "l_menthyl" if label == "1l" else "unresolved"
    molecule(f"m_{label}", label, f"o-aminotrifluoroacetophenone scope substrate {label}", (P[3],),
             stereo_kind=kind)
for label in ["4a", "4b", "4c", "4d", "4e", "d-3aa"]:
    note = "Single-crystal X-ray analysis reported; CCDC 1934884." if label == "4b" else None
    molecule(f"m_{label.replace('-', '_')}", label, label, (P[3], P[4]) if label == "4b" else (P[3],),
             stereo_kind="product_relative", note=note)

template_molecule("m_generic_1", "generic 1", "o-aminoaryl trifluoromethyl ketone with R1 and N-protecting group PG variables", ["R1", "PG"], (P[1], P[2]))
template_molecule("m_generic_2", "generic 2", "beta-prime-acetoxy allene bearing CO2R", ["R"], (P[1], P[2], P[3]))
template_molecule("m_generic_3", "generic 3", "CF3-containing tetrahydropyrano[3,2-b]indole with R/PG variables", ["R", "R1", "PG"], (P[1], P[2], P[3]))


simple = {
    "PPh3": ("triphenylphosphine", "P(c1ccccc1)(c1ccccc1)c1ccccc1"),
    "DABCO": ("1,4-diazabicyclo[2.2.2]octane", "C1CN2CCN1CC2"),
    "DMAP": ("4-dimethylaminopyridine", "CN(C)c1ccncc1"),
    "Cs2CO3": ("cesium carbonate", "[Cs+].[Cs+].[O-]C([O-])=O"),
    "Na2CO3": ("sodium carbonate", "[Na+].[Na+].[O-]C([O-])=O"),
    "K2CO3": ("potassium carbonate", "[K+].[K+].[O-]C([O-])=O"),
    "Et3N": ("triethylamine", "CCN(CC)CC"),
    "CH3ONa": ("sodium methoxide", "C[O-].[Na+]"),
    "CHCl3": ("chloroform", "ClC(Cl)Cl"),
    "acetone": ("acetone", "CC(C)=O"),
    "THF": ("tetrahydrofuran", "C1CCOC1"),
    "toluene": ("toluene", "Cc1ccccc1"),
    "MCPBA": ("m-chloroperbenzoic acid", "O=C(OO)c1cccc(Cl)c1"),
    "DCE": ("1,2-dichloroethane", "ClCCCl"),
    "EA": ("ethyl acetate", "CCOC(C)=O"),
    "MeCN": ("acetonitrile", "CC#N"),
    "H2O": ("water", "O"),
    "MeOH": ("methanol", "CO"),
    "PhNH2": ("aniline", "Nc1ccccc1"),
    "EDCI": ("1-ethyl-3-(3-dimethylaminopropyl)carbodiimide", "CCN=C=NCCCN(C)C"),
    "DCM": ("dichloromethane", "ClCCl"),
    "D2O": ("deuterium oxide", "[2H]O[2H]"),
}
for key, (name, smi) in simple.items():
    ev = (P[2], P[3]) if key in {"DMAP", "K2CO3", "CHCl3", "THF", "toluene"} else ((P[2],) if key in {"PPh3", "DABCO", "Cs2CO3", "Na2CO3", "Et3N", "CH3ONa", "acetone"} else (P[3],))
    molecule(f"m_{key}", key, name, ev, smi)
for key, name in {
    "RuCl3": "ruthenium(III) chloride", "NaIO4": "sodium periodate",
    "LiOH": "lithium hydroxide", "LiAlH4": "lithium aluminium hydride",
    "CH3MgBr": "methylmagnesium bromide", "PhMgBr": "phenylmagnesium bromide"
}.items():
    molecule(f"m_{key}", key, name, (P[3],), note="Named reagent; complete salt/solvation structure not encoded.")


bundle["molecules"] = list(MOLECULES.values())
REACTIONS = []


def participant(pid, mid, role, ev, *, amount=None, equiv=None, molpct=None,
                reference=None, product_role=None, method="original_extraction"):
    p = st.new_participant(pid, mid, role, ev)
    if product_role:
        p["product_role"] = f(product_role, product_role, ev)
    if amount:
        p["charges"]["amount"] = f(amount[0], amount[1], ev, method=method,
                                      derivation="Inherited from the displayed general conditions." if method == "inherited_general_condition" else None)
    if equiv:
        p["charges"]["equivalents"] = f(equiv[0], equiv[1], ev, method=method,
                                          derivation="Inherited from the displayed general conditions." if method == "inherited_general_condition" else None)
    if molpct:
        p["charges"]["mol_percent"] = f(molpct[0], molpct[1], ev, method=method,
                                         derivation="Inherited from the displayed general conditions." if method == "inherited_general_condition" else None)
    if reference:
        p["equivalence_reference"] = inferred(reference, ev,
            "The displayed equivalents/mol% are interpreted relative to the limiting substrate named by the scheme/table.")
    return p


def add_solvent(r, key, ev):
    pid = f"p_{key}_{len(r['participants'])+1}"
    r["participants"].append(participant(pid, f"m_{key}", "solvent", ev))
    return pid


def base_reaction(rid, ev, record_type, categories, location, name=None, transformation=None):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, ev) if categories is not None else f()
    r["location"] = f(location, location, ev)
    r["reaction_name"] = f(name, name, ev) if name else f()
    if transformation:
        r["transformation_type"] = inferred(transformation, ev, "Transformation type summarized from the source reaction arrow/prose.")
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    return r


def set_step(r, ev, pids, temp=None, time=None, operation="reaction"):
    s = st.new_step(f"step_{len(r['steps'])+1}", ev)
    s["participant_ids"] = list(pids)
    s["operation"] = f(operation, operation, ev)
    if temp is not None:
        s["temperature"] = f(temp[0], temp[1], ev)
    if time is not None:
        s["time"] = f(time[0], time[1], ev)
    r["steps"].append(s)


def solvent_conditions(r, solvent, ev, volume=None, concentration=None, inherited=False):
    method = "inherited_general_condition" if inherited else "original_extraction"
    derivation = "Inherited from the scheme's displayed general conditions." if inherited else None
    r["solvent_conditions"]["composition"] = f(solvent, [solvent], ev, method=method, derivation=derivation)
    if volume:
        r["solvent_conditions"]["total_volume"] = f(volume[0], volume[1], ev, method=method, derivation=derivation)
    if concentration:
        r["solvent_conditions"]["concentration"] = f(concentration[0], concentration[1], ev, method=method, derivation=derivation)
        r["solvent_conditions"]["concentration_basis"] = inferred("limiting substrate", ev,
            "The article discusses system concentration relative to the 0.1 mmol limiting substrate.")


def changes_and_labels(kinds):
    changes = []
    labels = {}
    for index, (a, b, change_type, before, after) in enumerate(kinds, 1):
        cid = f"b{index}"
        ch = {"change_id": cid, "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
              "before_order": before, "after_order": after, "change_type": change_type}
        changes.append(ch)
        if change_type in ("formed", "hydrogen_gain", "order_increased"):
            lab = st.bond_label(a, b, after)
            labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(cid)
    return changes, list(labels.values())


def classify_partial(r, ev, kinds, raw_basis, limitations):
    changes, labels = changes_and_labels(kinds)
    why = "Partial net bond-type assignment from the supplied source drawing/prose; exact atom mapping and complete hydrogen/coproduct accounting were not reconstructed."
    r["classification"]["bond_changes"] = [derived(ch, ev, why) for ch in changes]
    r["classification"]["labels"] = derived(labels, ev, why, raw=raw_basis)
    r["classification"]["coverage"] = inferred("partial", ev, why)
    r["classification"]["limitations"] = list(limitations)


def classify_main(r, ev):
    classify_partial(r, ev,
        [("C", "C", "formed", 0, 1), ("C", "N", "formed", 0, 1), ("C", "O", "formed", 0, 1)],
        "The article states that one C-C, one C-N, and one C-O bond are sequentially built.",
        ["Product atom correspondence is not encoded.", "Eliminated acetate/coproduct and all hydrogen transfers are unresolved."])


def classify_na(r, ev, why):
    r["classification"]["labels"] = inferred([], ev, why)
    r["classification"]["coverage"] = inferred("not_applicable", ev, why)
    r["classification"]["limitations"] = [why]


def outcome(r, oid, kind, targets, ev, raw, normalized, *, yield_type=None,
            basis=None, product_role=None, scope=None):
    o = st.new_outcome(oid, kind, targets, ev)
    o["value"] = f(raw, normalized, ev)
    if yield_type:
        o["yield_type"] = f(yield_type, yield_type, ev)
    if basis:
        o["measurement_basis"] = f(basis, basis, ev)
    if product_role:
        o["product_role"] = f(product_role, product_role, ev)
    if scope:
        o["scope"] = f(scope, scope, ev)
    r["outcomes"].append(o)


# Generic overview and prior-work context.
r = base_reaction("generic_domino_annulation", (P[1], P[2]), "generic_scheme", None,
                  "Graphical abstract and Scheme 1c", "[4+1]/[3+3] domino sequential annulation")
r["participants"] = [participant("p_g1", "m_generic_1", "reactant", (P[1], P[2])),
                     participant("p_g2", "m_generic_2", "reactant", (P[1], P[2])),
                     participant("p_g3", "m_generic_3", "product", (P[1], P[2]), product_role="desired product"),
                     participant("p_dmap", "m_DMAP", "catalyst", (P[1], P[2]), molpct=("30 mol%", q(30, "mol%")), reference="m_generic_1"),
                     participant("p_k2co3", "m_K2CO3", "reagent", (P[1], P[2]), equiv=("1.0 equiv", q(1, "equiv")), reference="m_generic_1")]
add_solvent(r, "CHCl3", (P[1], P[2]))
set_step(r, (P[1], P[2]), [p["participant_id"] for p in r["participants"]], temp=("rt", "room temperature"))
solvent_conditions(r, "CHCl3", (P[1], P[2]))
classify_main(r, (P[1], P[2]))
r["review"]["unresolved"] = ["Generic R/PG variables are not concrete molecule identities.", "Exact atom mapping and generic yield are not applicable."]
REACTIONS.append(r)

r = base_reaction("scheme1_prior_modes", (P[2],), "literature_example", None,
                  "Scheme 1a (previous work)", transformation="Lewis-base-bound 1,3-diene nucleophilic addition modes")
r["participants"] = [participant("p_2a", "m_2a", "reactant", (P[2],))]
r["notes"] = [f("Previous-work phosphonium/ammonium diene pathways; the beta-prime addition is described as well reported.",
                  "Previous-work context, not a current-paper experiment.", (P[2],))]
r["classification"]["coverage"] = inferred("unresolved", (P[2],), "Generic literature mechanism lacks a concrete mapped product.")
r["review"]["unresolved"] = ["Cited prior-work scheme is generic and has no concrete experiment conditions in this paper."]
REACTIONS.append(r)

for rid, loc in [("generic_scheme2", "Scheme 2 generic arrow"), ("generic_scheme3", "Scheme 3 generic arrow")]:
    ev = (P[3],)
    r = base_reaction(rid, ev, "generic_scheme", None, loc, "[4+1]/[3+3] domino sequential annulation")
    reactants = ("m_1a", "m_generic_2") if rid.endswith("2") else ("m_generic_1", "m_2a")
    r["participants"] = [participant("p_r1", reactants[0], "reactant", ev), participant("p_r2", reactants[1], "reactant", ev),
                         participant("p_prod", "m_generic_3", "product", ev, product_role="desired product"),
                         participant("p_dmap", "m_DMAP", "catalyst", ev, molpct=("30 mol%", q(30, "mol%")), reference=reactants[0]),
                         participant("p_base", "m_K2CO3", "reagent", ev, equiv=("1.0 equiv", q(1, "equiv")), reference=reactants[0])]
    add_solvent(r, "CHCl3", ev)
    set_step(r, ev, [p["participant_id"] for p in r["participants"]], temp=("rt", "room temperature"))
    solvent_conditions(r, "CHCl3", ev, concentration=("0.05 M", q(0.05, "M")))
    classify_main(r, ev)
    r["review"]["unresolved"] = ["Generic scope arrow; specific products are recorded separately."]
    REACTIONS.append(r)


# Table 1: each optimization/control row is a separate experiment.
table_rows = [
    (1, "PPh3", "Cs2CO3", 1.2, "CHCl3", "0", None, None),
    (2, "DABCO", "Cs2CO3", 1.2, "CHCl3", "0", None, None),
    (3, "DMAP", "Cs2CO3", 1.2, "CHCl3", "44", 44, ">20:1"),
    (4, None, "Cs2CO3", 1.2, "CHCl3", "0", None, None),
    (5, "DMAP", "Na2CO3", 1.2, "CHCl3", "57", 57, ">20:1"),
    (6, "DMAP", "K2CO3", 1.2, "CHCl3", "61", 61, ">20:1"),
    (7, "DMAP", "Et3N", 1.2, "CHCl3", "48", 48, ">20:1"),
    (8, "DMAP", "CH3ONa", 1.2, "CHCl3", "trace", "trace", ">20:1"),
    (9, "DMAP", None, None, "CHCl3", "45", 45, ">20:1"),
    (10, "DMAP", "K2CO3", 1.2, "acetone", "28", 28, ">20:1"),
    (11, "DMAP", "K2CO3", 1.2, "THF", "35", 35, ">20:1"),
    (12, "DMAP", "K2CO3", 1.2, "toluene", "11", 11, ">20:1"),
    (13, "DMAP", "K2CO3", 1.0, "CHCl3", "64", 64, ">20:1"),
    (14, "DMAP", "K2CO3", 1.0, "CHCl3", "73", 73, ">20:1"),
    (15, "DMAP", "K2CO3", 1.0, "CHCl3", "68", 68, ">20:1"),
    (16, "DMAP", "K2CO3", 1.0, "CHCl3", "58", 58, ">20:1"),
    (17, "DMAP", "K2CO3", 1.0, "CHCl3", "48", 48, ">20:1"),
    (18, "DMAP", "K2CO3", 1.0, "CHCl3", "71", 71, ">20:1"),
]
for entry, catalyst, base, base_eq, solvent, raw_yield, yval, dr in table_rows:
    ev = (P[2], T[2]) if entry != 14 else (P[2], T[2], P[3], T[3])
    r = base_reaction(f"table1_entry{entry:02d}", ev, "specific_experiment",
                      ["optimization"] + (["control"] if entry in (1, 2, 4) else []),
                      f"Table 1, entry {entry}", "[4+1]/[3+3] domino sequential annulation")
    r["general_scheme_ids"] = [f"{DOC}:generic_domino_annulation"]
    r["participants"] = [participant("p_1a", "m_1a", "reactant", ev, amount=("0.1 mmol", q(0.1, "mmol"))),
                         participant("p_2a", "m_2a", "reactant", ev,
                                     amount=None if entry == 16 else ("0.2 mmol", q(0.2, "mmol")),
                                     equiv=("1.5 equiv", q(1.5, "equiv")) if entry == 16 else ("2.0 equiv", q(2, "equiv")),
                                     reference="m_1a")]
    if yval is not None:
        r["participants"].append(participant("p_3aa", "m_3aa", "product", ev, product_role="desired product"))
    if catalyst:
        r["participants"].append(participant("p_cat", f"m_{catalyst}", "catalyst", ev,
                                             molpct=("30 mol%", q(30, "mol%")), reference="m_1a"))
    if base:
        r["participants"].append(participant("p_base", f"m_{base}", "reagent", ev,
                                             equiv=(f"{base_eq} equiv", q(base_eq, "equiv")), reference="m_1a"))
    add_solvent(r, solvent, ev)
    volume = 2.0 if entry in (14, 16, 17, 18) else 4.0 if entry == 15 else 1.0
    temp = q(0, "°C") if entry == 17 else q(40, "°C") if entry == 18 else "room temperature"
    temp_raw = "0 °C" if entry == 17 else "40 °C" if entry == 18 else "room temperature"
    set_step(r, ev, [p["participant_id"] for p in r["participants"]], temp=(temp_raw, temp), time=("18 h", q(18, "h")))
    conc = None
    if entry == 14:
        conc = ("0.05 M", q(0.05, "M"))
    solvent_conditions(r, solvent, ev, volume=(f"{volume:g} mL", q(volume, "mL")), concentration=conc)
    r["reaction_scale"] = f("1a (0.1 mmol)", {"reference_molecule_id": "m_1a", "amount": q(0.1, "mmol")}, ev)
    if yval is None:
        outcome(r, "o_no_reaction", "no_reaction", [], ev,
                {"yield_cell": "0", "prose": "no desired product observed"}, True,
                scope={"desired_product": "3aa"})
        classify_na(r, ev, "No desired product was observed; expected bond formation is not assigned to a negative control.")
        r["review"]["unresolved"] = ["Detection method and limit for the zero-yield/no-product observation are unreported."]
    elif yval == "trace":
        outcome(r, "o_trace", "trace", ["m_3aa"], ev, "trace", "trace", basis="yield column; detection method unreported", product_role="desired product")
        outcome(r, "o_dr", "selectivity", ["m_3aa"], ev, dr, dr, basis="1H NMR spectroscopy",
                scope={"metric": "dr", "major_minor_order": "unassigned"})
        classify_main(r, ev)
        r["review"]["unresolved"] = ["Trace-yield detection/measurement basis is unreported.", "Major/minor diastereomer identities are unassigned."]
    else:
        outcome(r, "o_yield", "yield", ["m_3aa"], ev, f"{raw_yield}%", q(yval, "%"), yield_type="isolated", basis="isolated yield", product_role="desired product")
        outcome(r, "o_dr", "selectivity", ["m_3aa"], ev, dr, dr, basis="1H NMR spectroscopy",
                scope={"metric": "dr", "major_minor_order": "unassigned"})
        classify_main(r, ev)
        r["review"]["unresolved"] = ["Exact product atom mapping and absolute stereochemistry are unresolved.", "Major/minor diastereomer identities are unassigned."]
    if entry == 14:
        r["general_scheme_ids"].append(f"{DOC}:generic_scheme2")
        r["notes"].append(f("Scheme 2 repeats 3aa under the same displayed 18 h/73%/>20:1 dr standard conditions.",
                            "Scheme 2 3aa is treated as duplicate evidence for Table 1 entry 14, not a second experiment.", (P[2], P[3])))
    REACTIONS.append(r)


def scope_record(rid, substrate_mid, product_mid, hours, yield_pct, generic_id, substrate_side):
    ev = (P[3],)
    r = base_reaction(rid, ev, "specific_experiment", ["substrate_scope"],
                      f"Scheme {'2' if substrate_side == 'allenoate' else '3'}, product {product_mid[2:]}",
                      "[4+1]/[3+3] domino sequential annulation")
    r["general_scheme_ids"] = [f"{DOC}:{generic_id}"]
    limiting = "m_1a" if substrate_side == "allenoate" else substrate_mid
    other = substrate_mid if substrate_side == "allenoate" else "m_2a"
    r["participants"] = [participant("p_lim", limiting, "reactant", ev, amount=("0.1 mmol", q(0.1, "mmol"))),
                         participant("p_other", other, "reactant", ev, amount=("0.2 mmol", q(0.2, "mmol")), equiv=("2.0 equiv", q(2, "equiv")), reference=limiting),
                         participant("p_prod", product_mid, "product", ev, product_role="desired product"),
                         participant("p_dmap", "m_DMAP", "catalyst", ev, molpct=("30 mol%", q(30, "mol%")), reference=limiting),
                         participant("p_base", "m_K2CO3", "reagent", ev, amount=("0.1 mmol", q(0.1, "mmol")), equiv=("1.0 equiv", q(1, "equiv")), reference=limiting)]
    add_solvent(r, "CHCl3", ev)
    set_step(r, ev, [p["participant_id"] for p in r["participants"]], temp=("room temperature", "room temperature"), time=(f"{hours} h", q(hours, "h")))
    solvent_conditions(r, "CHCl3", ev, volume=("2.0 mL", q(2, "mL")), concentration=("0.05 M", q(0.05, "M")))
    r["reaction_scale"] = f("limiting substrate (0.1 mmol)", {"reference_molecule_id": limiting, "amount": q(0.1, "mmol")}, ev)
    outcome(r, "o_yield", "yield", [product_mid], ev, f"{yield_pct}% yield", q(yield_pct, "%"), product_role="desired product")
    outcome(r, "o_dr", "selectivity", [product_mid], ev, ">20:1 dr", ">20:1", basis="1H NMR spectroscopy", scope={"metric": "dr", "major_minor_order": "unassigned"})
    classify_main(r, ev)
    r["review"]["unresolved"] = ["Yield measurement basis is not stated in the supplied main-text Scheme caption.", "Complete structures, atom mapping, and absolute stereochemistry remain unresolved."]
    return r


scheme2 = {
    "3ab": (16, 69), "3ac": (16, 73), "3ad": (36, 68), "3ae": (36, 81),
    "3af": (36, 89), "3ag": (36, 98), "3ah": (24, 71), "3ai": (16, 61),
    "3aj": (16, 64), "3ak": (16, 68), "3al": (16, 71), "3am": (24, 74),
}
for label, (hours, ypct) in scheme2.items():
    allene = f"m_2{label[-1]}"
    REACTIONS.append(scope_record(f"scheme2_{label}", allene, f"m_{label}", hours, ypct, "generic_scheme2", "allenoate"))

scheme3 = {
    "3ba": (18, 68), "3ca": (18, 71), "3da": (20, 90), "3ea": (16, 69),
    "3fa": (16, 77), "3ga": (18, 68), "3ha": (16, 78), "3ia": (16, 74),
    "3ja": (16, 73), "3ka": (48, 70), "3la": (18, 79), "3ma": (16, 72),
    "3na": (16, 70), "3oa": (12, 64), "3pa": (8, 51),
}
for label, (hours, ypct) in scheme3.items():
    substrate = f"m_1{label[1]}"
    REACTIONS.append(scope_record(f"scheme3_{label}", substrate, f"m_{label}", hours, ypct, "generic_scheme3", "aminoketone"))


# Scheme 4 scale-up.
ev = (P[3], T[3])
r = base_reaction("scheme4_gram_scale_3aa", ev, "specific_experiment", ["scale_up"], "Scheme 4, gram-scale synthesis", "[4+1]/[3+3] domino sequential annulation")
r["general_scheme_ids"] = [f"{DOC}:generic_domino_annulation"]
r["participants"] = [participant("p_1a", "m_1a", "reactant", ev, amount=("3.0 mmol", q(3, "mmol"))),
                     participant("p_2a", "m_2a", "reactant", ev, equiv=("standard conditions: 2.0 equiv", q(2, "equiv")), reference="m_1a", method="inherited_general_condition"),
                     participant("p_3aa", "m_3aa", "product", ev, product_role="desired product"),
                     participant("p_dmap", "m_DMAP", "catalyst", ev, molpct=("30 mol%", q(30, "mol%")), reference="m_1a"),
                     participant("p_base", "m_K2CO3", "reagent", ev, equiv=("1.0 equiv", q(1, "equiv")), reference="m_1a")]
add_solvent(r, "CHCl3", ev)
set_step(r, ev, [p["participant_id"] for p in r["participants"]], temp=("rt", "room temperature"))
solvent_conditions(r, "CHCl3", ev)
r["reaction_scale"] = f("1a (3.0 mmol)", {"reference_molecule_id": "m_1a", "amount": q(3, "mmol")}, ev)
outcome(r, "o_mass", "product_amount", ["m_3aa"], ev, "1.32 g", q(1.32, "g"), product_role="desired product")
outcome(r, "o_yield", "yield", ["m_3aa"], ev, "83% yield", q(83, "%"), product_role="desired product")
outcome(r, "o_dr", "selectivity", ["m_3aa"], ev, "excellent diastereoselectivity", "excellent diastereoselectivity", scope={"metric": "dr", "major_minor_order": "unassigned"})
classify_main(r, ev)
r["review"]["unresolved"] = ["Reaction time, isolated-yield basis, and numerical dr are not reported in the supplied main-text scale-up arrow.", "The 2a loading is inherited from 'standard reaction conditions'; its scaled amount is not calculated into a reported charge."]
REACTIONS.append(r)


def derivatization(rid, product, reagents, solvents, temp, yield_pct, dr, kinds, location, steps=None, conflict=None):
    ev = (P[3], T[3])
    r = base_reaction(rid, ev, "specific_experiment", ["derivatization"], location, transformation=rid.replace("scheme4_", ""))
    r["participants"] = [participant("p_3aa", "m_3aa", "reactant", ev), participant("p_product", product, "product", ev, product_role="desired product")]
    for i, key in enumerate(reagents, 1):
        r["participants"].append(participant(f"p_reagent_{i}", f"m_{key}", "reagent", ev))
    for key in solvents:
        add_solvent(r, key, ev)
    if steps:
        for op, step_reagents, step_solvents, step_temp in steps:
            ids = [p["participant_id"] for p in r["participants"] if p["molecule_id"] in {"m_3aa", product, *(f"m_{x}" for x in step_reagents), *(f"m_{x}" for x in step_solvents)}]
            set_step(r, ev, ids, temp=step_temp, operation=op)
    else:
        set_step(r, ev, [p["participant_id"] for p in r["participants"]], temp=temp)
    solvent_conditions(r, "/".join(solvents), ev)
    outcome(r, "o_yield", "yield", [product], ev, f"{yield_pct}% yield", q(yield_pct, "%"), product_role="desired product",
            scope={"span": "3aa to product"})
    if dr:
        outcome(r, "o_dr", "selectivity", [product], ev, dr, dr, scope={"metric": "dr", "major_minor_order": "unassigned"})
    classify_partial(r, ev, kinds, location, ["No complete atom mapping was reconstructed.", "Hydrogen/proton sources and coproducts are incomplete unless explicitly represented by the partial net assignment."])
    r["review"]["unresolved"] = ["Reagent amounts, reaction time, workup, purification, and yield measurement basis are not reported in the supplied main text.", "Complete atom mapping and absolute stereochemistry remain unresolved."]
    if conflict:
        c = f([conflict[0], conflict[1]], None, ev, status="conflict")
        c["alternatives"] = [f(conflict[0], conflict[0], (P[3],)), f(conflict[1], conflict[1], (T[3],))]
        r["procedure"]["modifications"] = c
        r["notes"].append(f("Scheme 4 labels CH3MgBr; the prose says PhMgBr.", "Unresolved source conflict in reagent identity.", ev))
        r["review"]["unresolved"].append("Source conflict: Scheme 4 shows CH3MgBr, whereas the prose says PhMgBr; neither was selected as the definitive reagent.")
    return r


REACTIONS.append(derivatization("scheme4_epoxidation_4a", "m_4a", ["MCPBA"], ["DCE"], ("50 °C", q(50, "°C")), 46,
    "excellent diastereoselectivity", [("C", "O", "formed", 0, 1), ("C", "O", "formed", 0, 1), ("C", "C", "order_decreased", 2, 1)], "Scheme 4, 3aa to 4a"))
REACTIONS.append(derivatization("scheme4_dihydroxylation_4b", "m_4b", ["RuCl3", "NaIO4"], ["EA", "MeCN", "H2O"], ("0 °C", q(0, "°C")), 80,
    ">20:1", [("C", "O", "formed", 0, 1), ("C", "O", "formed", 0, 1), ("O", "H", "hydrogen_gain", 0, 1), ("O", "H", "hydrogen_gain", 0, 1), ("C", "C", "order_decreased", 2, 1)], "Scheme 4, 3aa to 4b"))
REACTIONS.append(derivatization("scheme4_hydrolysis_amidation_4c", "m_4c", ["LiOH", "PhNH2", "EDCI", "DMAP"], ["H2O", "MeOH", "THF", "DCM"], None, 46,
    ">20:1", [("C", "N", "formed", 0, 1), ("C", "O", "broken", 1, 0)], "Scheme 4, two-step 3aa to 4c (46% total yield)",
    steps=[("hydrolysis", ["LiOH"], ["H2O", "MeOH", "THF"], ("65 °C", q(65, "°C"))), ("amidation", ["PhNH2", "EDCI", "DMAP"], ["DCM"], ("0 °C", q(0, "°C")))]))
REACTIONS.append(derivatization("scheme4_ester_reduction_4d", "m_4d", ["LiAlH4"], ["THF"], ("rt", "room temperature"), 74,
    ">20:1", [("C", "O", "order_decreased", 2, 1), ("C", "O", "broken", 1, 0), ("C", "H", "hydrogen_gain", 0, 1), ("C", "H", "hydrogen_gain", 0, 1), ("O", "H", "hydrogen_gain", 0, 1)], "Scheme 4, 3aa to 4d"))
REACTIONS.append(derivatization("scheme4_grignard_4e", "m_4e", ["CH3MgBr", "PhMgBr"], ["toluene"], ("0 °C", q(0, "°C")), 77,
    ">20:1", [("C", "C", "formed", 0, 1), ("C", "C", "formed", 0, 1), ("C", "O", "order_decreased", 2, 1), ("C", "O", "broken", 1, 0), ("O", "H", "hydrogen_gain", 0, 1)], "Scheme 4, 3aa to 4e", conflict=("CH3MgBr (Scheme 4)", "PhMgBr (prose)")))


# Scheme 5 deuterium-labeling experiment.
ev = (P[3], T[3])
r = base_reaction("scheme5_deuterium_labeling", ev, "specific_experiment", ["mechanistic_study"], "Scheme 5", "deuterium-labeling domino annulation")
r["participants"] = [participant("p_1a", "m_1a", "reactant", ev), participant("p_2a", "m_2a", "reactant", ev),
                     participant("p_d3aa", "m_d_3aa", "product", ev, product_role="deuterated product"),
                     participant("p_dmap", "m_DMAP", "catalyst", ev, molpct=("30 mol%", q(30, "mol%")), reference="m_1a"),
                     participant("p_base", "m_K2CO3", "reagent", ev, equiv=("1.0 equiv", q(1, "equiv")), reference="m_1a"),
                     participant("p_d2o", "m_D2O", "reagent", ev, equiv=("20 equiv", q(20, "equiv")), reference="m_1a")]
add_solvent(r, "CHCl3", ev)
set_step(r, ev, [p["participant_id"] for p in r["participants"]], temp=("rt", "room temperature"))
solvent_conditions(r, "CHCl3", ev)
outcome(r, "o_yield", "yield", ["m_d_3aa"], ev, "46%", q(46, "%"), product_role="deuterated product")
outcome(r, "o_beta_d", "other", ["m_d_3aa"], ev, "beta D 68%", q(68, "%"), scope={"metric": "deuterium incorporation", "position": "beta"})
outcome(r, "o_gamma_d", "other", ["m_d_3aa"], ev, "gamma D 20%", q(20, "%"), scope={"metric": "deuterium incorporation", "position": "gamma"})
classify_main(r, ev)
r["review"]["unresolved"] = ["Substrate amounts, reaction time, yield basis, and deuterium-incorporation measurement method are unreported in the supplied main text.", "Exact isotope atom mapping and the extent of H/D exchange in recovered/unreacted material are unresolved."]
REACTIONS.append(r)


# Scheme 6 is an author-proposed pathway, not a demonstrated experiment.
ev = (P[4], T[4])
r = base_reaction("scheme6_proposed_mechanism", ev, "proposed_mechanism", None, "Scheme 6", "plausible reaction mechanism")
r["participants"] = [participant("p_1a", "m_1a", "reactant", ev), participant("p_2a", "m_2a", "reactant", ev), participant("p_3aa", "m_3aa", "product", ev),
                     participant("p_dmap", "m_DMAP", "catalyst", ev), participant("p_base", "m_K2CO3", "reagent", ev)]
r["notes"] = [f("Intermediates I–VI and paths a/b are proposed; intermediate VI was not isolated, and the authors favor the concerted Diels–Alder pathway based on diastereoselectivity.",
                  {"intermediates": ["I", "II", "III", "IV", "V", "VI"], "paths": ["a", "b"], "VI_isolated": False, "favored": "concerted Diels-Alder"}, ev)]
classify_main(r, ev)
r["review"]["unresolved"] = ["Proposed elementary steps are not experimental records and are not claimed as verified.", "Intermediate VI was not isolated; path preference is the authors' interpretation."]
REACTIONS.append(r)


bundle["reactions"] = REACTIONS
bundle["source_context"] = [
    {"context_id": f"{DOC}:figure1_pharmaceuticals", "document_id": DOC, "kind": "nonreaction_figure",
     "evidence_ids": [P[1]], "summary": "Figure 1 shows example pharmaceuticals/scaffolds; it is not a synthesis performed in this paper."},
    {"context_id": f"{DOC}:scheme1b_unexplored", "document_id": DOC, "kind": "research_gap",
     "evidence_ids": [P[2]], "summary": "Scheme 1b marks gamma nucleophilic addition of an amine-bound diene as unexplored, not an experiment."},
    {"context_id": f"{DOC}:xray", "document_id": DOC, "kind": "characterization",
     "evidence_ids": [P[3], P[4]], "summary": "Single-crystal X-ray analyses are reported for 3ag (CCDC 1951502) and 4b (CCDC 1934884); these are characterization, not extra syntheses."},
    {"context_id": f"{DOC}:si_notice", "document_id": DOC, "kind": "supplement_notice",
     "evidence_ids": [P[4]], "summary": "The article advertises SI, but no SI file was supplied/associated in this workspace."},
    {"context_id": f"{DOC}:bibliography", "document_id": DOC, "kind": "bibliography",
     "evidence_ids": [P[4], P[5], P[6]], "summary": "Pages 5–6 and the lower part of page 4 contain cited literature only; no current-paper experiments were extracted from the bibliography."},
]


# Raw page-level extraction snapshot, intentionally separate from the formal bundle.
raw = {
    "document_id": DOC,
    "pages_reviewed": [1, 2, 3, 4, 5, 6],
    "table1_rows": table_rows,
    "scheme2_scope": scheme2,
    "scheme3_scope": scheme3,
    "scheme4_specific_records": ["gram_scale_3aa", "4a", "4b", "4c", "4d", "4e"],
    "scheme5": {"product": "d-3aa", "yield_percent": 46, "beta_D_percent": 68, "gamma_D_percent": 20},
    "source_conflicts": [{"scheme": "CH3MgBr", "prose": "PhMgBr", "record": f"{DOC}:scheme4_grignard_4e"}],
    "deduplication": {"Scheme 2 3aa": f"{DOC}:table1_entry14"},
    "record_count": len(REACTIONS),
    "specific_experiment_count": sum(r["record_type"] == "specific_experiment" for r in REACTIONS),
    "limitations": ["No SI supplied.", "Dense source structures remain unresolved rather than receiving unsafe SMILES.", "All records remain needs_review."]
}
st.write(RUN / "raw_model/pages_0001_0006_extraction.json", raw)
st.write(RUN / "raw_model/stereo_checks.json", STEREO_CHECKS)
st.write(RUN / "submitted.json", bundle)


coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title, graphical abstract, abstract", "finding": "Current-paper generic [4+1]/[3+3] annulation; overview conditions and one C-N/one C-C/one C-O claim.", "evidence_ids": [P[1], T[1]]},
            {"region": "Figure 1 and body introduction", "finding": "Pharmaceutical examples and background only; no current synthesis experiment.", "evidence_ids": [P[1], T[1]]}
        ], "record_ids": [f"{DOC}:generic_domino_annulation"], "unresolved": ["Generic R1/R/PG identities and exact atom mapping are unresolved."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Scheme 1a–c", "finding": "Prior generic Lewis-base modes, an explicitly unexplored mode, and the current generic annulation design.", "evidence_ids": [P[2], T[2]]},
            {"region": "Table 1, entries 1–18 and footnotes a–h", "finding": "Eighteen independent optimization/control experiments, including three zero-product controls, one trace result, solvent/volume/allene-loading/temperature exceptions, isolated-yield basis, and 1H NMR dr method.", "evidence_ids": [P[2], T[2]]}
        ], "record_ids": [f"{DOC}:scheme1_prior_modes", f"{DOC}:generic_domino_annulation"] + [f"{DOC}:table1_entry{i:02d}" for i in range(1, 19)], "unresolved": ["Zero-yield/trace detection limits are unreported.", "Exact product atom mapping and absolute stereochemistry are unresolved."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Scheme 2 and footnote", "finding": "Generic allenoate scope plus 13 products; 3aa duplicates Table 1 entry 14 and 12 other products are independent experiments.", "evidence_ids": [P[3], T[3]]},
            {"region": "Scheme 3 and footnote", "finding": "Generic aminoketone scope plus 15 independent product experiments with individual times/yields and >20:1 dr.", "evidence_ids": [P[3], T[3]]},
            {"region": "Scheme 4 and adjacent prose", "finding": "One gram-scale synthesis and five independent derivatizations; 4c is a two-step overall 46% transformation; Scheme/prose conflict CH3MgBr versus PhMgBr retained.", "evidence_ids": [P[3], T[3]]},
            {"region": "Scheme 5 and adjacent prose", "finding": "One deuterium-labeling mechanistic experiment: d-3aa 46%, beta D 68%, gamma D 20%.", "evidence_ids": [P[3], T[3]]}
        ], "record_ids": [f"{DOC}:generic_scheme2", f"{DOC}:generic_scheme3", f"{DOC}:table1_entry14"] + [f"{DOC}:scheme2_{x}" for x in scheme2] + [f"{DOC}:scheme3_{x}" for x in scheme3] + [f"{DOC}:scheme4_gram_scale_3aa", f"{DOC}:scheme4_epoxidation_4a", f"{DOC}:scheme4_dihydroxylation_4b", f"{DOC}:scheme4_hydrolysis_amidation_4c", f"{DOC}:scheme4_ester_reduction_4d", f"{DOC}:scheme4_grignard_4e", f"{DOC}:scheme5_deuterium_labeling"], "unresolved": ["Scheme 2/3/4/5 yield measurement bases are not stated in the supplied main text.", "Most complete structures and absolute stereochemical assignments remain unresolved.", "CH3MgBr/PhMgBr source conflict remains unresolved."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "Scheme 6 and mechanism discussion", "finding": "Author-proposed intermediates I–VI and paths a/b; VI was not isolated; mechanism stored separately from experiments.", "evidence_ids": [P[4], T[4]]},
            {"region": "Associated Content and Accession Codes", "finding": "SI is advertised but not supplied; CCDC entries characterize 3ag and 4b and do not create new syntheses.", "evidence_ids": [P[4], T[4]]},
            {"region": "lower-right references", "finding": "Bibliographic context only.", "evidence_ids": [P[4], T[4]]}
        ], "record_ids": [f"{DOC}:scheme6_proposed_mechanism"], "unresolved": ["Proposed pathway is not experimentally verified.", "Unavailable SI prevents extraction of procedures/workups/purifications."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "two-column references (continuation)", "finding": "Bibliography only; cited prior allene chemistry was not reclassified as current-paper experiments.", "evidence_ids": [P[5], T[5]]}
        ], "record_ids": [], "unresolved": []},
        {"page": 6, "state": "attempted", "regions": [
            {"region": "references 16d–g and CCDC notes 17–18", "finding": "Bibliography and crystallographic accession notes only; no additional current-paper reaction arrow or procedure.", "evidence_ids": [P[6], T[6]]}
        ], "record_ids": [], "unresolved": []},
    ]
}
st.write(RUN / "coverage/pages_0001_0006.json", coverage)

tasks = json.loads((RUN / "tasks.json").read_text())
for task in tasks:
    task["state"] = "attempted"
    task["reason"] = "Original page image and extracted text were reviewed; see coverage/pages_0001_0006.json. Attempted does not certify semantic accuracy."
st.write(RUN / "tasks.json", tasks)

print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
                  "specific_experiments": raw["specific_experiment_count"], "stereo_checks": len(STEREO_CHECKS)}, ensure_ascii=False))
