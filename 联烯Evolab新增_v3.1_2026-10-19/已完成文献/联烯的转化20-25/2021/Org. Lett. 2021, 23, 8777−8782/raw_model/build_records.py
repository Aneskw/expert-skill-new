#!/usr/bin/env python3
"""Build the v3.1 machine bundle from the visually transcribed main-text pages."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DATA = json.loads((RUN / "dataset.json").read_text())
DOC = "doc_0c91dd7389d7f800"
P = {n: f"{DOC}:page{n}:image" for n in range(1, 7)}
PT = {n: f"{DOC}:page{n}:text" for n in range(1, 7)}
GENERAL = f"{DOC}:general_thiocarbocyclization"
_UNSET = object()


def f(raw=None, norm=_UNSET, ev=(), method="original_extraction", derivation=None,
      status=None):
    if norm is _UNSET:
        norm = raw
    return st.fact(raw, norm, ev, method=method, derivation=derivation, status=status)


def inferred(norm, ev, why, raw=None):
    return f(raw, norm, ev, "inference", why)


def derived(norm, ev, why, raw=None):
    return f(raw, norm, ev, "structure_derivation", why)


def set_stereo(m, ev, chirality="unresolved", geometry="unresolved",
               representation="unresolved", reason="The source identity has not been fully encoded."):
    m["stereochemistry"]["chirality_status"] = inferred(chirality, ev, reason)
    m["stereochemistry"]["geometry_status"] = inferred(geometry, ev, reason)
    m["stereochemistry"]["representation_status"] = inferred(representation, ev, reason)


def base_molecule(suffix, label, ev, name=None):
    mid = f"{DOC}:mol_{suffix}"
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    m["occurrences"][0]["evidence_ids"] = list(ev)
    m["names"]["original"] = f(name or label, name or label, ev)
    m["names"]["english"] = f(name or label, name or label, ev)
    m["names"]["abbreviations"] = f(label, [label], ev)
    m["cas_numbers"] = f()
    m["external_properties"] = []
    return m


def unresolved_molecule(suffix, label, ev, name=None, kind="other", stereo=None):
    m = base_molecule(suffix, label, ev, name)
    m["representation_type"] = "unresolved"
    m["smiles"] = f()
    reason = "Dense or source-labelled structure was not converted to a complete stereochemically reliable SMILES from the supplied main-text figure."
    m["identity_note"] = inferred(reason, ev, reason)
    if stereo:
        set_stereo(m, ev, **stereo)
    elif kind == "allene":
        set_stereo(m, ev, "unspecified", "not_applicable", "not_encoded",
                   "The allene drawing has no axial descriptor; absence of a descriptor is not an achiral assignment.")
        m["stereochemistry"]["axial"] = inferred(
            {"axis": "C=C=C", "configuration": "unspecified"}, ev,
            "The source drawing contains an allene axis but gives no R_a/S_a or P/M assignment.")
    elif kind == "product":
        set_stereo(m, ev, "unspecified", "unresolved", "not_encoded",
                   "The product drawing carries no explicit stereochemical assignment; potential centers were not collapsed to achiral.")
    else:
        set_stereo(m, ev, "unresolved", "unresolved", "unresolved", reason)
    return m


def smiles_molecule(suffix, label, ev, smiles, name=None, allene_racemate=False):
    m = base_molecule(suffix, label, ev, name)
    diag = st.stereo_inspect(smiles)
    m["representation_type"] = "smiles"
    m["smiles"] = derived(smiles, ev,
        "Candidate SMILES transcribed from the visible structure/name and checked for parsing; parsing does not verify source identity.",
        raw={"source_label": label, "source_location": "main-text drawing or explicit name"})
    m["stereochemistry"]["atom_centers"] = derived(diag["atom_centers"], ev, "RDKit graph diagnostic on the candidate SMILES.")
    m["stereochemistry"]["double_bonds"] = derived(diag["double_bonds"], ev, "RDKit graph diagnostic on the candidate SMILES.")
    m["stereochemistry"]["identity_key"] = derived(diag["stereo_graph_key"], ev, "SHA-256 key of canonical isomeric graph; not a sample identity.")
    if allene_racemate:
        m["stereochemistry"]["axial"] = f("racemic allene", {"axis": "C=C=C", "configuration": "racemic, absolute configurations unassigned"}, ev)
        m["stereochemistry"]["mixture"] = f("racemic", "racemic", ev)
        m["stereochemistry"]["assignment"] = f("racemic allene; no R_a/S_a assignment", "racemic axial allene, configurations unassigned", ev)
        set_stereo(m, ev, "racemic", "not_applicable", "partial",
                   "The source explicitly calls the allene racemic, while the candidate SMILES/toolkit does not encode its allene axis.")
        m["stereochemistry"]["isomeric_smiles"] = f()
    else:
        no_potential = not diag["potential_stereo"]
        chir = "achiral" if no_potential else "unspecified"
        rep = "not_applicable" if no_potential else "partial"
        set_stereo(m, ev, chir, "not_applicable", rep,
                   "Graph diagnostics found no assignable stereogenic element." if no_potential else "Candidate graph contains unspecified potential stereochemistry.")
        if no_potential:
            m["stereochemistry"]["assignment"] = derived("no stereogenic element detected in candidate graph", ev, "RDKit potential-stereo diagnostic; source identity remains machine-extracted.")
        m["stereochemistry"]["isomeric_smiles"] = derived(diag["canonical_isomeric_smiles"], ev, "Canonical isomeric SMILES computed from the candidate graph.")
    m["identity_note"] = inferred("machine candidate; source-labelled identity retained and chemistry remains needs_review", ev,
                                  "SMILES parsing and graph diagnostics do not establish that the drawing was transcribed correctly.")
    return m


def generic_molecule(suffix, label, ev, fixed, variables):
    m = base_molecule(suffix, label, ev, label)
    m["representation_type"] = "generic_template"
    m["smiles"] = f()
    m["generic_template"] = {
        "fixed_structure": f(fixed, {"description": fixed}, ev),
        "variables": [{
            "label": v,
            "attachment_points": inferred("source drawing attachment point", ev, "Attachment is retained descriptively because the template is not atom mapped."),
            "definition": f(v, v, ev),
            "allowed_values": f(),
            "constraints": f(),
        } for v in variables],
        "instance_mapping": [],
        "reason": f("generic R-group drawing", "Source is a generic drawing rather than one concrete molecule.", ev),
    }
    set_stereo(m, ev, "unresolved", "unresolved", "not_encoded",
               "Generic R-group structure cannot support a unique stereochemical identity.")
    m["identity_note"] = inferred("generic source template; no concrete SMILES", ev,
                                  "Unexpanded R groups are retained as variables rather than wildcards in a concrete SMILES.")
    return m


bundle = st.empty(DATA["documents"])
molecules = []


def add(m):
    molecules.append(m)
    return m["molecule_id"]


# Core concrete candidates and common conditions.
M = {}
M["1a"] = add(smiles_molecule("1a", "1a", [P[2], PT[2], P[3], PT[3]],
    "CC(=C=Cc1ccccc1)c1ccccc1", "racemic diphenyl allene 1a", allene_racemate=True))
M["2a"] = add(smiles_molecule("2a", "2a", [P[2], PT[2], P[3]],
    "CCCSN1C(=O)CCC1=O", "1-(propylthio)pyrrolidine-2,5-dione"))
M["C1"] = add(smiles_molecule("C1", "C1", [P[2], PT[2]],
    "c1ccc([Se]c2ccccc2)cc1", "diphenyl selenide"))
M["TMSOTf"] = add(smiles_molecule("TMSOTf", "TMSOTf", [P[2], P[3], PT[2], PT[3]],
    "C[Si](C)(C)OS(=O)(=O)C(F)(F)F", "trimethylsilyl trifluoromethanesulfonate"))
M["CH2Cl2"] = add(smiles_molecule("CH2Cl2", "CH2Cl2", [P[2], P[3]], "ClCCl", "dichloromethane"))
M["std"] = add(smiles_molecule("nmr_internal_standard", "1,3,5-trimethoxybenzene", [P[2], PT[2]],
    "COc1cc(OC)cc(OC)c1", "1,3,5-trimethoxybenzene"))

for label, name in {
    "C2": "selenide catalyst C2", "C3": "selenide catalyst C3", "C4": "selenide catalyst C4",
    "TESOTf": "triethylsilyl triflate", "BF3OEt2": "boron trifluoride diethyl etherate",
    "TfOH": "triflic acid", "Tf2NH": "bis(trifluoromethanesulfonyl)imide",
    "PhSPh": "diphenyl sulfide", "Ph3PS": "triphenylphosphine sulfide",
    "Ph3PSe": "triphenylphosphine selenide",
}.items():
    M[label] = add(unresolved_molecule(label, label, [P[2], PT[2]], name))

M["3a"] = add(unresolved_molecule("3a", "3a", [P[1], P[2], P[3]], "indene-based propyl sulfide 3a", "product"))

# Sulfur-reagent scope candidates. 2f is retained unresolved because the footnote changes scaffold to saccharin.
scope2_reagents = {
    "3b": ("2b", "CCCCCCSN1C(=O)CCC1=O", "n-hexylthio electrophile shown for 3b"),
    "3c": ("2c", "CC(C)CCSN1C(=O)CCC1=O", "branched alkylthio electrophile shown for 3c"),
    "3d": ("2d", "ClCCCSN1C(=O)CCC1=O", "3-chloropropylthio electrophile shown for 3d"),
    "3e": ("2e", "c1ccccc1CSN2C(=O)CCC2=O", "benzylthio electrophile shown for 3e"),
    "3g": ("2g", "c1ccc(SN2C(=O)CCC2=O)cc1", "phenylthio electrophile shown for 3g"),
    "3h": ("2h", "Fc1ccccc1SN2C(=O)CCC2=O", "fluorophenylthio electrophile shown for 3h"),
    "3i": ("2i", "Cc1ccc(SN2C(=O)CCC2=O)cc1", "p-tolylthio electrophile 2i"),
    "3j": ("2j", "Cc1ccccc1SN2C(=O)CCC2=O", "o-tolylthio electrophile shown for 3j"),
}
for product, (label, smi, name) in scope2_reagents.items():
    M[label] = add(smiles_molecule(label, label, [P[2]], smi, name))
M["2f"] = add(unresolved_molecule("2f", "N-trifluoromethylthiosaccharin", [P[2]], "N-trifluoromethylthiosaccharin"))

for label in ["3b", "3c", "3d", "3e", "3f", "3g", "3h", "3i", "3j"]:
    M[label] = add(unresolved_molecule(label, label, [P[2]], f"indene-based sulfide {label}", "product"))

# Allene scope identities and their products. 1f also occurs in the kinetic-resolution discussion.
scope3_products = ["3k", "3l", "3m", "3n", "3o", "3p", "3q", "3r", "3s", "3t", "3u", "3v", "3w", "3x", "3y"]
scope3_substrates = [f"1{c}" for c in "bcdefghijklmnop"]
for substrate in scope3_substrates:
    ev = [P[3], P[4], PT[4]] if substrate == "1f" else [P[3]]
    m = unresolved_molecule(substrate, substrate, ev, f"allene substrate {substrate}", "allene")
    if substrate == "1f":
        m["stereochemistry"]["axial"] = f("racemic allene 1f", {"axis": "C=C=C", "configuration": "racemic, unassigned"}, [P[3], P[4], PT[4]])
        m["stereochemistry"]["mixture"] = f("racemic allene 1f", "racemic", [P[4], PT[4]])
        m["stereochemistry"]["chirality_status"] = f("racemic allene 1f", "racemic", [P[4], PT[4]])
    M[substrate] = add(m)
for product in scope3_products:
    M[product] = add(unresolved_molecule(product, product, [P[3]], f"indene-based sulfide {product}", "product"))

for label, name, kind in [
    ("1q", "meta-methoxy-substituted allene 1q", "allene"),
    ("3z", "regioisomeric indene-based sulfide 3z", "product"),
    ("1r", "allene substrate 1r", "allene"),
    ("4", "azidothiolation product 4", "product"),
    ("4prime", "thiocarbocyclization byproduct 4-prime", "product"),
    ("3aa", "kinetic-resolution cyclized product 3aa", "product"),
]:
    ev = [P[3], P[4], PT[4]] if label == "3aa" else [P[3], PT[3]]
    m = unresolved_molecule(label, "4′" if label == "4prime" else label, ev, name, kind)
    if label == "3aa":
        m["stereochemistry"]["mixture"] = f("0% ee / without ee", "racemic", [P[3], P[4], PT[4]])
        m["stereochemistry"]["chirality_status"] = f("0% ee / without ee", "racemic", [P[3], P[4], PT[4]])
        m["stereochemistry"]["assignment"] = f("stereocenter marked with an asterisk; absolute configuration unassigned", "racemic product with unassigned enantiomer identities", [P[3]])
    M[label] = add(m)

M["TMSN3"] = add(unresolved_molecule("TMSN3", "TMSN3", [P[3]], "trimethylsilyl azide"))
for cat in ["C5", "C6", "C7"]:
    m = unresolved_molecule(cat, cat, [P[3], P[4]], f"chiral catalyst {cat}", stereo={
        "chirality": "specified", "geometry": "not_applicable", "representation": "not_encoded",
        "reason": "The source drawing specifies relative/absolute-looking wedge stereochemistry, but no textual R/S assignment or reliable complete SMILES is available."
    })
    m["stereochemistry"]["atom_centers"] = f("wedge/dash stereochemistry shown in catalyst drawing", [{"configuration": "source-drawn; R/S unassigned"}], [P[3]])
    m["stereochemistry"]["assignment"] = f("source-drawn wedge/dash; no R/S descriptor", "specified in drawing, absolute label unresolved", [P[3]])
    M[cat] = add(m)

# Generic identities used only for current-family, prior-work, and mechanism records.
M["g_allene"] = add(generic_molecule("generic_allene", "generic allene", [P[1]], "Ar-C(R1)=C=C(R2)- framework", ["Ar", "R1", "R2"]))
M["g_alkene"] = add(generic_molecule("generic_alkene", "generic alkene", [P[1]], "R1-substituted C=C-R2 framework", ["R1", "R2"]))
M["g_alkyne"] = add(generic_molecule("generic_alkyne", "generic alkyne", [P[1]], "R1-C#C-R2 framework", ["R1", "R2"]))
M["g_SR"] = add(generic_molecule("generic_sulfur_electrophile", "RS+", [P[1]], "electrophilic sulfur fragment RS+", ["R"]))
M["g_sulfide"] = add(generic_molecule("generic_sulfide_product", "sulfides", [P[1]], "unspecified sulfide product family", ["R", "R1", "R2"]))
M["g_indene"] = add(generic_molecule("generic_indene_product", "generic indene-based sulfide 3", [P[1], P[3], P[4]], "indene scaffold bearing SR, R1 and R2", ["Ar", "R", "R1", "R2"]))

bundle["molecules"] = molecules


def participant(pid, mid, role, ev, stage="reaction", product_role=None):
    p = st.new_participant(pid, mid, role, ev, stage)
    if product_role:
        p["product_role"] = f(product_role, product_role, ev)
    return p


def charge(p, key, raw, norm, ev, ref=None, method="original_extraction"):
    p["charges"][key] = f(raw, norm, ev, method)
    if key in ("equivalents", "mol_percent"):
        p["equivalence_reference"] = f(ref, ref, ev, method)


def classification(ev, codes=None, kind="partial"):
    if not codes:
        return {
            "standard_id": st.STANDARD, "mapped_reaction_smiles": f(), "bond_changes": [],
            "labels": inferred([], ev, "No source-supported atom-complete bond changes were assigned."),
            "coverage": inferred(kind, ev, "The schematic is generic, contextual, or has no observed transformation to classify."),
            "limitations": ["No atom-complete mapping was created from the generic/source-context drawing."], "model_votes": []
        }
    changes, labels = [], []
    for i, code in enumerate(codes, 1):
        a, b = code.split("-")
        change = {"change_id": f"b{i}", "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
                  "before_order": 0, "after_order": 1, "change_type": "formed"}
        changes.append(change)
        lab = st.bond_label(a, b, 1)
        lab["supporting_change_ids"] = [f"b{i}"]
        labels.append(lab)
    why = "Net connectivity was inferred from the visible reactant/product drawings; atom correspondence, bond-order changes, hydrogens, leaving groups, and coproducts remain unresolved."
    return {
        "standard_id": st.STANDARD, "mapped_reaction_smiles": f(),
        "bond_changes": [derived(c, ev, why) for c in changes],
        "labels": derived(labels, ev, why), "coverage": inferred("partial", ev, why),
        "limitations": [why, "No three-model vote was performed."], "model_votes": []
    }


def reaction(rid, ev, record_type="specific_experiment", categories=None, location="", general=True,
             name="Catalytic electrophilic thiocarbocyclization of allenes", transform="thiocarbocyclization"):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, ev) if categories is not None else f()
    r["location"] = f(location, {"description": location}, ev)
    r["general_scheme_ids"] = [GENERAL] if general and record_type == "specific_experiment" else []
    r["reaction_name"] = f(name, name, ev)
    r["transformation_type"] = inferred(transform, ev, "Normalized transformation label inferred from the authors' reaction name and drawings.")
    r["review"]["unresolved"] = [
        "No explicitly associated Supporting Information was supplied; addition order, atmosphere, workup, purification, and characterization details remain unavailable.",
        "Machine extraction and visual transcription remain needs_review; schema validation is not chemistry verification."
    ]
    return r


def step(r, ev, temp=None, time=None, inherited=False):
    s = st.new_step("reaction", ev)
    s["participant_ids"] = [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"]
    method = "inherited_general_condition" if inherited else "original_extraction"
    s["operation"] = f("reaction under the displayed arrow conditions", "combine/react", ev, method)
    if temp is not None:
        s["temperature"] = f(f"{temp} °C", {"number": temp, "unit": "°C"}, ev, method)
    if time is not None:
        s["time"] = f(f"{time} h", {"number": time, "unit": "h"}, ev, method)
    r["steps"] = [s]


def outcome(r, oid, kind, targets, raw, norm, ev, basis=None, yield_type=None, scope=None):
    o = st.new_outcome(oid, kind, targets, ev)
    o["product_role"] = f("main product" if kind in ("yield", "product_amount") else kind,
                            "main_product" if kind in ("yield", "product_amount") else kind, ev)
    o["value"] = f(raw, norm, ev)
    if yield_type:
        o["yield_type"] = f(yield_type, yield_type, ev)
    if basis:
        o["measurement_basis"] = f(basis, basis, ev)
    if scope:
        o["scope"] = f(scope, scope, ev)
    r["outcomes"].append(o)


def set_solvent(r, ev, volume=None, inherited=False):
    method = "inherited_general_condition" if inherited else "original_extraction"
    r["solvent_conditions"]["composition"] = f("CH2Cl2", [M["CH2Cl2"]], ev, method)
    if volume is not None:
        r["solvent_conditions"]["total_volume"] = f(f"{volume} mL", {"number": volume, "unit": "mL"}, ev, method)


# Generic current-paper family.
r = reaction("general_thiocarbocyclization", [P[1], PT[1], P[2], P[3]], "generic_scheme", None,
             "graphical abstract and Schemes 2-3 general transformation", False)
r["participants"] = [
    participant("allene", M["g_allene"], "reactant", [P[1]]),
    participant("sulfur", M["g_SR"], "reactant", [P[1]]),
    participant("product", M["g_indene"], "product", [P[1]], product_role="product family"),
    participant("catalyst", M["C1"], "catalyst", [P[1]]),
    participant("additive", M["TMSOTf"], "additive", [P[1]]),
    participant("solvent", M["CH2Cl2"], "solvent", [P[1]])]
charge(r["participants"][3], "mol_percent", "catalytic selenide / 10 mol% in schemes", {"number": 10, "unit": "mol%"}, [P[1], P[2]], M["g_allene"])
step(r, [P[1]], -60, 22)
set_solvent(r, [P[1]])
r["classification"] = classification([P[1]], ["C", "S"] if False else ["C-S", "C-C"])
r["review"]["unresolved"].append("Generic drawing has variable R/Ar groups and is not a concrete experiment.")
bundle["reactions"].append(r)

# Scheme 1 prior-work/context arrows, retained as literature examples rather than current experiments.
for rid, substrate, loc in [
    ("literature_scheme1_alkene", "g_alkene", "Scheme 1 top-left alkene thiolation family"),
    ("literature_scheme1_alkyne", "g_alkyne", "Scheme 1 top-right alkyne thiolation family"),
    ("literature_scheme1_allene", "g_allene", "Scheme 1 bottom allene thiolation marked underdeveloped"),
]:
    q = reaction(rid, [P[1], PT[1]], "literature_example", None, loc, False,
                 "Prior electrophilic thiolation context", "literature-context thiolation")
    q["participants"] = [participant("substrate", M[substrate], "reactant", [P[1]]),
                         participant("electrophile", M["g_SR"], "reagent", [P[1]]),
                         participant("product_family", M["g_sulfide"], "product", [P[1]], product_role="generic literature product family")]
    q["classification"] = classification([P[1]], None, "unresolved")
    q["review"]["unresolved"].append("Scheme 1 is a generic literature-context schematic without a specific cited experimental row, conditions, or yield.")
    bundle["reactions"].append(q)

# Table 1 optimization rows. Entry 6 is split into NMR screening and the explicit 0.1 mmol isolated repeat.
table_rows = [
    (1, "C1", "TMSOTf", 25, 12, 0), (2, "C1", "TMSOTf", 0, 12, 18),
    (3, "C1", "TMSOTf", -20, 12, 41), (4, "C1", "TMSOTf", -60, 12, 66),
    (5, "C1", "TMSOTf", -78, 12, 61), (6, "C1", "TMSOTf", -60, 22, 89),
    (7, "C1", "TESOTf", -60, 22, 70), (8, "C1", "BF3OEt2", -60, 22, 41),
    (9, "C1", "TfOH", -60, 22, 56), (10, "C1", "Tf2NH", -60, 22, 55),
    (11, "C1", "TMSOTf", -60, 22, 63), (12, "C1", "TMSOTf", -60, 22, 69),
    (13, "C2", "TMSOTf", -60, 22, 70), (14, "C3", "TMSOTf", -60, 22, 67),
    (15, "C4", "TMSOTf", -60, 22, 68), (16, "PhSPh", "TMSOTf", -60, 22, 68),
    (17, "Ph3PS", "TMSOTf", -60, 22, 49), (18, "Ph3PSe", "TMSOTf", -60, 22, 18),
    (19, "C1", None, -60, 22, 0), (20, None, "TMSOTf", -60, 22, 61),
]
for entry, cat, additive, temp, hours, yld in table_rows:
    ev = [P[2], PT[2]]
    q = reaction(f"table1_entry{entry:02d}_nmr", ev, categories=["optimization"], location=f"Table 1 entry {entry}")
    p1 = participant("allene_1a", M["1a"], "reactant", ev); charge(p1, "amount", "0.05 mmol", {"number": 0.05, "unit": "mmol"}, ev)
    p2 = participant("sulfur_2a", M["2a"], "reactant", ev); charge(p2, "amount", "0.05 mmol", {"number": 0.05, "unit": "mmol"}, ev)
    q["participants"] = [p1, p2]
    if cat:
        pc = participant("catalyst", M[cat], "catalyst", ev); charge(pc, "mol_percent", "10 mol%", {"number": 10, "unit": "mol%"}, ev, M["1a"]); q["participants"].append(pc)
    if additive:
        pa = participant("additive", M[additive], "additive", ev)
        amount = 0.025 if entry == 11 else 0.1 if entry == 12 else 0.05
        charge(pa, "amount", f"{amount:g} mmol", {"number": amount, "unit": "mmol"}, ev); q["participants"].append(pa)
    ps = participant("solvent", M["CH2Cl2"], "solvent", ev); charge(ps, "volume", "1 mL", {"number": 1, "unit": "mL"}, ev); q["participants"].append(ps)
    q["participants"].append(participant("nmr_standard", M["std"], "reagent", ev, "characterization"))
    q["participants"].append(participant("product_3a", M["3a"], "product", ev, product_role="desired product"))
    q["reaction_scale"] = f("1a (0.05 mmol)", {"reference_molecule_id": M["1a"], "amount": {"number": 0.05, "unit": "mmol"}}, ev)
    set_solvent(q, ev, 1)
    step(q, ev, temp, hours)
    q["procedure"]["original_text"] = f("1a (0.05 mmol), 2a (0.05 mmol), catalyst (10 mol%), additive (0.05 mmol unless footnoted), CH2Cl2 (1 mL)",
                                               "Table 1 footnote conditions", ev)
    outcome(q, "yield_nmr", "yield", [M["3a"]], f"{yld}%", {"number": yld, "unit": "%"}, ev,
            "NMR yield using 1,3,5-trimethoxybenzene as internal standard", "NMR yield")
    if yld == 0:
        outcome(q, "no_desired_product", "no_reaction", [], "0% / no desired product observed", "no desired product observed", ev,
                scope="desired product 3a")
    q["classification"] = classification(ev, None, "not_applicable") if yld == 0 else classification(ev, ["C-S", "C-C"])
    if entry in (11, 12):
        q["notes"].append(f("entry-specific TMSOTf amount in footnote", f"TMSOTf amount {0.025 if entry == 11 else 0.1:g} mmol", ev))
    bundle["reactions"].append(q)

# Explicit isolated repeat associated with entry 6.
ev = [P[2], PT[2]]
q = reaction("table1_entry06_isolated_0p1mmol", ev, categories=["optimization", "scale_up"], location="Table 1 entry 6 parenthetical isolated repeat")
p1 = participant("allene_1a", M["1a"], "reactant", ev); charge(p1, "amount", "0.1 mmol scale", {"number": 0.1, "unit": "mmol"}, ev)
q["participants"] = [p1, participant("sulfur_2a", M["2a"], "reactant", ev), participant("catalyst", M["C1"], "catalyst", ev),
                      participant("additive", M["TMSOTf"], "additive", ev), participant("solvent", M["CH2Cl2"], "solvent", ev),
                      participant("product_3a", M["3a"], "product", ev, product_role="desired product")]
charge(q["participants"][2], "mol_percent", "10 mol%", {"number": 10, "unit": "mol%"}, ev, M["1a"])
q["reaction_scale"] = f("isolated yield in parentheses on a 0.1 mmol scale", {"reference_molecule_id": M["1a"], "amount": {"number": 0.1, "unit": "mmol"}}, ev)
set_solvent(q, ev); step(q, ev, -60, 22)
outcome(q, "yield_isolated", "yield", [M["3a"]], "85% (parentheses)", {"number": 85, "unit": "%"}, ev, "isolated yield", "isolated yield")
q["classification"] = classification(ev, ["C-S", "C-C"])
q["review"]["unresolved"].append("The footnote identifies a 0.1 mmol isolated repeat but does not reprint scaled absolute amounts for 2a, TMSOTf, or CH2Cl2; they remain unreported rather than multiplied by assumption.")
bundle["reactions"].append(q)

# Scheme 2 sulfur-electrophile scope.
scope2 = [("3b", "2b", 83), ("3c", "2c", 79), ("3d", "2d", 75), ("3e", "2e", 83),
          ("3f", "2f", 79), ("3g", "2g", 83), ("3h", "2h", 73), ("3i", "2i", 79), ("3j", "2j", 67)]
for product, reagent, yld in scope2:
    ev = [P[2], PT[2]]
    q = reaction(f"scheme2_{product}", ev, categories=["substrate_scope"], location=f"Scheme 2 product {product}")
    p1 = participant("allene_1a", M["1a"], "reactant", ev); charge(p1, "amount", "0.10 mmol", {"number": 0.10, "unit": "mmol"}, ev)
    p2 = participant("sulfur_reagent", M[reagent], "reactant", ev); charge(p2, "equivalents", "1.0 equiv", {"number": 1.0, "unit": "equiv"}, ev, M["1a"])
    pc = participant("catalyst_C1", M["C1"], "catalyst", ev); charge(pc, "mol_percent", "10 mol%", {"number": 10, "unit": "mol%"}, ev, M["1a"])
    pa = participant("TMSOTf", M["TMSOTf"], "additive", ev); charge(pa, "equivalents", "1.0 equiv", {"number": 1.0, "unit": "equiv"}, ev, M["1a"])
    ps = participant("CH2Cl2", M["CH2Cl2"], "solvent", ev); charge(ps, "volume", "2.0 mL", {"number": 2.0, "unit": "mL"}, ev)
    q["participants"] = [p1, p2, pc, pa, ps, participant("product", M[product], "product", ev, product_role="scope product")]
    q["reaction_scale"] = f("1a (0.10 mmol)", {"reference_molecule_id": M["1a"], "amount": {"number": 0.10, "unit": "mmol"}}, ev)
    set_solvent(q, ev, 2.0); step(q, ev, -60, 22)
    q["procedure"]["original_text"] = f("1a (0.10 mmol), 2 (1.0 equiv), C1 (10 mol%), TMSOTf (1.0 equiv), CH2Cl2 (2.0 mL), -60 °C, 22 h", "Scheme 2 footnote conditions", ev)
    outcome(q, "yield", "yield", [M[product]], f"{yld}%", {"number": yld, "unit": "%"}, ev, scope=f"product {product}")
    q["classification"] = classification(ev, ["C-S", "C-C"])
    q["review"]["unresolved"].append("Scheme 2 does not state the yield measurement basis in the supplied main article; it is not assumed to be isolated.")
    bundle["reactions"].append(q)

# Scheme 3 allene scope.
scope3_yields = [72, 83, 79, 82, 64, 62, 73, 68, 81, 82, 86, 61, 67, 99, 60]
for substrate, product, yld in zip(scope3_substrates, scope3_products, scope3_yields):
    ev = [P[3], PT[3]]
    q = reaction(f"scheme3_{product}", ev, categories=["substrate_scope"], location=f"Scheme 3 product {product}")
    pa1 = participant("allene", M[substrate], "reactant", ev); charge(pa1, "amount", "0.10 mmol", {"number": 0.10, "unit": "mmol"}, ev)
    reagent = "2i" if product in ("3x", "3y") else "2a"
    psr = participant("sulfur_reagent", M[reagent], "reactant", ev); charge(psr, "equivalents", "1.0 equiv", {"number": 1.0, "unit": "equiv"}, ev, M[substrate])
    pc = participant("catalyst_C1", M["C1"], "catalyst", ev); charge(pc, "mol_percent", "10 mol%", {"number": 10, "unit": "mol%"}, ev, M[substrate])
    pad = participant("TMSOTf", M["TMSOTf"], "additive", ev); charge(pad, "equivalents", "1.0 equiv", {"number": 1.0, "unit": "equiv"}, ev, M[substrate])
    psol = participant("CH2Cl2", M["CH2Cl2"], "solvent", ev); charge(psol, "volume", "2.0 mL", {"number": 2.0, "unit": "mL"}, ev)
    q["participants"] = [pa1, psr, pc, pad, psol, participant("product", M[product], "product", ev, product_role="scope product")]
    q["reaction_scale"] = f(f"{substrate} (0.10 mmol)", {"reference_molecule_id": M[substrate], "amount": {"number": 0.10, "unit": "mmol"}}, ev)
    set_solvent(q, ev, 2.0); step(q, ev, -60, 22)
    q["procedure"]["original_text"] = f("substrates 1b-1p (0.10 mmol), 2 (1.0 equiv), PhSePh (10 mol%), TMSOTf (1.0 equiv), CH2Cl2 (2.0 mL), -60 °C, 22 h", "Scheme 3 footnote conditions", ev)
    outcome(q, "yield", "yield", [M[product]], f"{yld}%", {"number": yld, "unit": "%"}, ev, scope=f"product {product}")
    q["classification"] = classification(ev, ["C-S", "C-C"])
    q["review"]["unresolved"].append("Scheme 3 does not state the yield measurement basis in the supplied main article; it is not assumed to be isolated.")
    bundle["reactions"].append(q)

# Eq 1 regioselectivity switch.
ev = [P[3], PT[3]]
q = reaction("eq1_3z", ev, categories=["substrate_scope"], location="page 3 equation 1 / product 3z")
q["participants"] = [participant("allene_1q", M["1q"], "reactant", ev), participant("sulfur_2a", M["2a"], "reactant", ev),
                      participant("catalyst_C1", M["C1"], "catalyst", ev), participant("TMSOTf", M["TMSOTf"], "additive", ev),
                      participant("CH2Cl2", M["CH2Cl2"], "solvent", ev), participant("product_3z", M["3z"], "product", ev, product_role="regioselective product")]
charge(q["participants"][2], "mol_percent", "10 mol%", {"number": 10, "unit": "mol%"}, ev, M["1q"])
set_solvent(q, ev); step(q, ev, -60, 22)
outcome(q, "yield", "yield", [M["3z"]], "61%", {"number": 61, "unit": "%"}, ev, scope="product 3z")
q["classification"] = classification(ev, ["C-S", "C-C"])
q["notes"].append(f("electron-rich aryl group acted as nucleophile", "authors attribute switched regioselectivity to the methoxy-substituted aryl group", ev))
q["review"]["unresolved"].append("Yield basis and absolute reagent amounts are not stated next to equation 1 in the supplied main article.")
bundle["reactions"].append(q)

# Eq 2 gram-scale reaction.
q = reaction("eq2_gram_scale_3a", ev, categories=["scale_up"], location="page 3 equation 2")
p1 = participant("allene_1a", M["1a"], "reactant", ev); charge(p1, "amount", "5.0 mmol", {"number": 5.0, "unit": "mmol"}, ev); charge(p1, "mass", "1.03 g", {"number": 1.03, "unit": "g"}, ev)
pc = participant("catalyst_C1", M["C1"], "catalyst", ev); charge(pc, "mol_percent", "10 mol%", {"number": 10, "unit": "mol%"}, ev, M["1a"])
q["participants"] = [p1, participant("sulfur_2a", M["2a"], "reactant", ev), pc, participant("TMSOTf", M["TMSOTf"], "additive", ev),
                      participant("CH2Cl2", M["CH2Cl2"], "solvent", ev), participant("product_3a", M["3a"], "product", ev, product_role="scale-up product")]
q["reaction_scale"] = f("1a: 5.0 mmol, 1.03 g", {"reference_molecule_id": M["1a"], "amount": {"number": 5.0, "unit": "mmol"}}, ev)
set_solvent(q, ev); step(q, ev, -60, 22)
outcome(q, "product_mass", "product_amount", [M["3a"]], "1.15 g", {"number": 1.15, "unit": "g"}, ev, scope="product 3a")
outcome(q, "yield", "yield", [M["3a"]], "82%", {"number": 82, "unit": "%"}, ev, scope="product 3a")
q["classification"] = classification(ev, ["C-S", "C-C"])
q["review"]["unresolved"].append("Equation 2 says similar conditions but does not reprint the absolute charges of 2a/TMSOTf or solvent volume; no scaling calculation was inserted. Yield basis is unreported.")
bundle["reactions"].append(q)

# Eq 3 intermolecular azidothiolation with an explicitly shown cyclization byproduct.
q = reaction("eq3_azidothiolation", ev, categories=["substrate_scope"], location="page 3 equation 3",
             name="intermolecular azidothiolation of allene 1r", transform="intermolecular azidothiolation")
pc = participant("catalyst_C1", M["C1"], "catalyst", ev); charge(pc, "mol_percent", "10 mol%", {"number": 10, "unit": "mol%"}, ev, M["1r"])
q["participants"] = [participant("allene_1r", M["1r"], "reactant", ev), participant("sulfur_2a", M["2a"], "reactant", ev),
                      participant("TMSN3", M["TMSN3"], "reagent", ev), pc, participant("TMSOTf", M["TMSOTf"], "additive", ev),
                      participant("CH2Cl2", M["CH2Cl2"], "solvent", ev), participant("product_4", M["4"], "product", ev, product_role="desired azidothiolation product"),
                      participant("product_4prime", M["4prime"], "product", ev, product_role="thiocarbocyclization byproduct")]
set_solvent(q, ev); step(q, ev, -60, 22)
outcome(q, "yield_4", "yield", [M["4"]], "52%", {"number": 52, "unit": "%"}, ev, scope="desired product 4")
outcome(q, "yield_4prime", "yield", [M["4prime"]], "20%", {"number": 20, "unit": "%"}, ev, scope="thiocarbocyclization byproduct 4-prime")
q["classification"] = classification(ev, ["C-S", "C-N"])
q["review"]["unresolved"].append("Yield bases, reagent amounts, and product ratio measurement method are unreported in the supplied main article; byproduct 4-prime is kept as a distinct product outcome.")
bundle["reactions"].append(q)


def kinetic_record(rid, catalyst, hours, yld=None, ee="no enantioselectivity", explicit=False):
    ev = [P[3], P[4], PT[4]]
    q = reaction(rid, ev, categories=["mechanistic_study"], location=f"equation 4 / kinetic-resolution test with {catalyst}",
                 name="kinetic resolution test of racemic allene 1f", transform="thiocarbocyclization / kinetic resolution test")
    pc = participant("catalyst", M[catalyst], "catalyst", ev); charge(pc, "mol_percent", "10 mol%", {"number": 10, "unit": "mol%"}, ev, M["1f"], "original_extraction" if explicit else "inherited_general_condition")
    q["participants"] = [participant("racemic_allene_1f", M["1f"], "reactant", ev), participant("sulfur_2i", M["2i"], "reactant", ev), pc,
                          participant("TMSOTf", M["TMSOTf"], "additive", ev), participant("CH2Cl2", M["CH2Cl2"], "solvent", ev),
                          participant("product_3aa", M["3aa"], "product", ev, product_role="cyclized product")]
    set_solvent(q, ev, inherited=not explicit); step(q, ev, -60, hours, inherited=not explicit)
    if yld is not None:
        outcome(q, "yield", "yield", [M["3aa"]], f"{yld}%", {"number": yld, "unit": "%"}, ev, scope="product 3aa")
    if ee == "0% ee":
        norm = {"number": 0, "unit": "%", "metric": "ee"}
    else:
        norm = {"metric": "ee", "result": "no enantioselectivity"}
    outcome(q, "enantioselectivity", "selectivity", [M["3aa"]], ee, norm, ev,
            scope="ee of product 3aa; enantiomer identities/order and analytical method unreported")
    q["classification"] = classification(ev, ["C-S", "C-C"])
    q["review"]["unresolved"].append("The ee analytical method and enantiomer elution/configuration order are not stated in the main article; yield basis is also unreported where a yield is shown.")
    return q


bundle["reactions"].append(kinetic_record("eq4_C5_22h", "C5", 22, 60, "0% ee", True))
bundle["reactions"].append(kinetic_record("kinetic_C5_1h", "C5", 1, 40, "without ee", True))
bundle["reactions"].append(kinetic_record("kinetic_C6_similar", "C6", 22, None, "no enantioselectivity", False))
bundle["reactions"].append(kinetic_record("kinetic_C7_similar", "C7", 22, None, "no enantioselectivity", False))

# Scheme 4 is the authors' proposal, not an observed experiment.
ev4 = [P[4], PT[4]]
q = reaction("scheme4_proposed_mechanism", ev4, "proposed_mechanism", None, "Scheme 4 proposed mechanism", False,
             "Proposed mechanism for catalytic electrophilic thiocarbocyclization", "proposed thiiranium-ion pathway")
q["participants"] = [participant("allene", M["g_allene"], "reactant", ev4), participant("sulfur", M["g_SR"], "reactant", ev4),
                      participant("catalyst", M["C1"], "catalyst", ev4), participant("additive", M["TMSOTf"], "additive", ev4),
                      participant("product", M["g_indene"], "product", ev4, product_role="generic product family")]
q["classification"] = classification(ev4, ["C-S", "C-C"])
q["notes"].append(f("possible mechanism is proposed", "authors propose activation, sulfonium intermediate II, thiiranium ion III, and aromatic-ring attack", ev4))
q["review"]["unresolved"].append("Intermediates I-IV and transition states are proposed; they are not isolated products or independently verified experiments.")
bundle["reactions"].append(q)

bundle["source_context"] = [
    {"context_id": f"{DOC}:graphical_abstract_duplicate", "document_id": DOC, "kind": "duplicate_schematic",
     "evidence_ids": [P[1]], "note": "Graphical abstract repeats the current generic transformation; it does not create an extra experiment."},
    {"context_id": f"{DOC}:si_advertisement", "document_id": DOC, "kind": "supplementary_material_notice",
     "evidence_ids": [P[4], PT[4]], "note": "The article advertises experimental details/characterization/NMR in SI, but no SI PDF was explicitly supplied or associated."},
    {"context_id": f"{DOC}:references_pages5_6", "document_id": DOC, "kind": "bibliography",
     "evidence_ids": [P[5], P[6]], "note": "Pages 5-6 are references; cited papers and titles were not converted into current-paper synthesis records."},
    {"context_id": f"{DOC}:spectra_absence", "document_id": DOC, "kind": "characterization_scope",
     "evidence_ids": [P[4]], "note": "No standalone mass-spectrum peak table or characterization table occurs in pages 1-6; characterization is only advertised as SI."},
]

# Raw transcription summary, kept separate from the formal bundle.
raw_rows = {
    "document_id": DOC,
    "source": "visual inspection of original rendered PDF pages 1-6",
    "table1": [{"entry": e, "catalyst": c or "dash", "additive": a or "dash", "temperature_C": t,
                "time_h": h, "yield_percent": y, "basis": "NMR; entry 6 also 85% isolated at 0.1 mmol" if e == 6 else "NMR"}
               for e, c, a, t, h, y in table_rows],
    "scheme2": [{"product": p, "yield_percent": y} for p, _, y in scope2],
    "scheme3": [{"substrate": s, "product": p, "yield_percent": y} for s, p, y in zip(scope3_substrates, scope3_products, scope3_yields)],
    "equations": [
        {"equation": 1, "product": "3z", "yield_percent": 61},
        {"equation": 2, "product": "3a", "scale": "5.0 mmol / 1.03 g 1a", "product_mass_g": 1.15, "yield_percent": 82},
        {"equation": 3, "products": [{"label": "4", "yield_percent": 52}, {"label": "4-prime", "yield_percent": 20}]},
        {"equation": 4, "catalyst": "C5", "time_h": 22, "yield_percent": 60, "ee_percent": 0},
        {"prose": "C5 1 h", "yield_percent": 40, "selectivity": "without ee"},
        {"prose": "C6 similar", "selectivity": "no enantioselectivity", "yield": None},
        {"prose": "C7 similar", "selectivity": "no enantioselectivity", "yield": None},
    ],
    "limitations": ["Machine visual transcription remains needs_review.", "No SI association was inferred."]
}

record_ids = [r["reaction_id"] for r in bundle["reactions"]]
by_prefix = lambda starts: [rid for rid in record_ids if any(rid.startswith(f"{DOC}:{s}") for s in starts)]
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "graphical abstract at top", "finding": "current generic thiocarbocyclization; duplicate family schematic, not a separate experiment", "evidence_ids": [P[1]]},
            {"region": "Scheme 1 lower-left", "finding": "three generic literature-context arrows for alkene, alkyne, and underdeveloped allene thiolation", "evidence_ids": [P[1]]},
            {"region": "two-column introduction", "finding": "background and motivation; no current specific experiment", "evidence_ids": [P[1], PT[1]]}],
         "record_ids": [GENERAL] + by_prefix(["literature_scheme1_"]),
         "unresolved": ["Generic R-group structures and possible intermediate are not concrete molecule identities."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Table 1 left column", "finding": "20 optimization rows; entry 6 additionally reports an isolated 85% repeat at 0.1 mmol", "evidence_ids": [P[2], PT[2]]},
            {"region": "Scheme 2 right column", "finding": "nine sulfur-reagent scope experiments 3b-3j with common conditions", "evidence_ids": [P[2], PT[2]]},
            {"region": "surrounding prose", "finding": "duplicates/interprets Table 1 and Scheme 2; no additional conditions beyond recorded entries", "evidence_ids": [P[2], PT[2]]}],
         "record_ids": by_prefix(["table1_", "scheme2_"]),
         "unresolved": ["Scheme 2 yield basis is not stated in the main article.", "Exact dense product structures remain source-labelled candidates."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Scheme 3 left column", "finding": "15 allene-scope experiments 3k-3y with common conditions", "evidence_ids": [P[3], PT[3]]},
            {"region": "equation 1 upper-right", "finding": "regioselectivity-switched product 3z, 61%", "evidence_ids": [P[3], PT[3]]},
            {"region": "equation 2 middle-right", "finding": "5.0 mmol/1.03 g scale-up to 1.15 g 3a, 82%", "evidence_ids": [P[3], PT[3]]},
            {"region": "equation 3 lower-right", "finding": "azidothiolation product 4 (52%) plus cyclization byproduct 4-prime (20%)", "evidence_ids": [P[3], PT[3]]},
            {"region": "equation 4 bottom-right", "finding": "C5 kinetic-resolution test: 3aa 60%, 0% ee; discussion continues on page 4", "evidence_ids": [P[3], P[4], PT[4]]}],
         "record_ids": by_prefix(["scheme3_", "eq1_", "eq2_", "eq3_", "eq4_"]),
         "unresolved": ["Yield bases for schemes/equations are not printed.", "Absolute configurations and ee analytical method are unreported."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "top-left continuation", "finding": "C6 and C7 no-enantioselectivity tests plus C5 1 h/40%/without ee", "evidence_ids": [P[4], PT[4]]},
            {"region": "Scheme 4 and adjacent prose", "finding": "authors' proposed sulfonium/thiiranium mechanism, not a verified experiment", "evidence_ids": [P[4], PT[4]]},
            {"region": "Associated Content", "finding": "SI is advertised but no SI file is supplied or associated", "evidence_ids": [P[4], PT[4]]}],
         "record_ids": by_prefix(["kinetic_", "eq4_", "scheme4_"]),
         "unresolved": ["C6/C7 yields are not reported.", "Proposed intermediates are not isolated.", "SI procedures and characterization are unavailable."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "full two-column page", "finding": "references 3-15 only; cited literature is context, not current-paper experiments", "evidence_ids": [P[5], PT[5]]}],
         "record_ids": [], "unresolved": ["No current-paper experimental arrow, table row, or characterization table appears on this page."]},
        {"page": 6, "state": "attempted", "regions": [
            {"region": "upper-left reference continuation", "finding": "references 15-19 only; remainder of page blank", "evidence_ids": [P[6], PT[6]]}],
         "record_ids": [], "unresolved": ["No current-paper experimental arrow, table row, or characterization table appears on this page."]},
    ]
}

stereo_checks = {}
for m in molecules:
    smi = st.value(m["smiles"])
    if smi:
        stereo_checks[m["molecule_id"]] = st.stereo_inspect(smi)

# Keep the standard-init output path, but replace its empty arrays with this source-derived bundle.
st.write(RUN / "submitted.json", bundle)
st.write(RUN / "raw_model/extracted_rows.json", raw_rows)
st.write(RUN / "raw_model/stereo_checks.json", stereo_checks)
st.write(RUN / "coverage/pages_0001_0006.json", coverage)

tasks = json.loads((RUN / "tasks.json").read_text())
for task in tasks:
    if task.get("document_id") == DOC and task.get("state") == "pending":
        task["state"] = "attempted"
        task["reason"] = "Original page image inspected; detailed regions and record mapping are in coverage/pages_0001_0006.json."
st.write(RUN / "tasks.json", tasks)

print(json.dumps({"molecules": len(molecules), "reactions": len(bundle["reactions"]),
                  "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in bundle["reactions"]),
                  "stereo_checks": len(stereo_checks)}, ensure_ascii=False))
