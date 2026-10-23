#!/usr/bin/env python3
"""Build the source-derived v3.1 bundle for ejoc.202000076, PDF pages 1--5.

The page images were inspected in this session.  Candidate structures below are
machine records and remain needs_review; unclear drawings are deliberately left
unresolved rather than expanded by convention.
"""
from __future__ import annotations

import copy
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_08b3dd4eacda5aa8"
E = {p: f"{DOC}:page{p}:image" for p in range(1, 6)}
T = {p: f"{DOC}:page{p}:text" for p in range(1, 6)}
TITLE = "Revisiting the Addition of in-situ Nucleophiles to Allenic Ketones: An Entry Towards Synthesis of Benzodioxins"


def f(raw=None, value=None, ev=(), method="original_extraction", derivation=None,
      status=None, alternatives=()):
    return st.fact(raw, value, ev, method=method, derivation=derivation,
                   status=status, alternatives=alternatives)


def inferred(value, ev, why):
    return f(None, value, ev, method="structure_derivation", derivation=why)


def q(number, unit, raw=None):
    return {"number": number, "unit": unit, "raw_text": raw or f"{number} {unit}"}


dataset = st.read(RUN / "dataset.json")
doc = dataset["documents"][0]
assert doc["document_id"] == DOC and doc["pages"] == 5
# This builder is the reproducible replacement for the mechanically prepared
# zero-record state.  Clear only machine records previously emitted by this
# same builder so removed candidates cannot linger across reruns.
dataset["experiments"] = []
dataset.pop("standard", None)
doc["title"] = TITLE
doc["supplementary_material"] = {
    "status": "mentioned_not_provided",
    "note": "Page 1 advertises Supporting Information online, but no SI file was supplied or explicitly associated in this session."
}
st.write(RUN / "dataset.json", dataset)
bundle = st.empty(dataset["documents"])
molecules = {}


def slug(text):
    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")


def add_molecule(key, label, name, smiles, ev, *, geometry="not_applicable",
                 chirality=None, assignment=None, axial=None, unresolved_reason=None):
    if key in molecules:
        m = molecules[key]
        for eid in ev:
            if eid not in m["occurrences"][0]["evidence_ids"]:
                m["occurrences"][0]["evidence_ids"].append(eid)
        return m["molecule_id"]
    mid = f"{DOC}:mol:{key}"
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    m["names"]["original"] = f(name, name, ev)
    m["names"]["english"] = f(name, name, ev)
    m["names"]["abbreviations"] = f()
    m["names"]["aliases"] = f(label, [label], ev)
    m["cas_numbers"] = f()
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = inferred(smiles, ev, "Candidate connectivity transcribed from the displayed structure; parsing does not verify the drawing.")
        try:
            diag = st.stereo_inspect(smiles)
        except Exception:
            diag = None
        if chirality is None:
            potential = bool(diag and any(x["type"].startswith("Atom_") for x in diag["potential_stereo"]))
            chirality = "unspecified" if potential else "achiral"
        m["stereochemistry"]["chirality_status"] = inferred(chirality, ev, "Assigned from the candidate graph; allene axial cases are handled explicitly and are not inferred from RDKit silence.")
        m["stereochemistry"]["geometry_status"] = inferred(geometry, ev, "Geometry state follows the source label/drawing; not all source E labels are encoded in the candidate SMILES.")
        rep = "partial" if geometry == "specified" or chirality in {"specified", "partially_specified", "unspecified"} else "encoded"
        if chirality == "achiral" and geometry == "not_applicable":
            rep = "not_applicable"
        m["stereochemistry"]["representation_status"] = inferred(rep, ev, "Candidate SMILES preserves connectivity; explicit source stereochemistry is only encoded when unambiguous.")
        m["stereochemistry"]["isomeric_smiles"] = inferred(diag["canonical_isomeric_smiles"], ev, "RDKit canonical isomeric form of the candidate graph; source verification remains pending.") if diag else f()
        m["stereochemistry"]["identity_key"] = inferred(diag["stereo_graph_key"], ev, "Graph key computed from the candidate isomeric SMILES; it is not a sample identity.") if diag else f()
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = f()
        chirality = chirality or "unresolved"
        m["stereochemistry"]["chirality_status"] = inferred(chirality, ev, unresolved_reason or "The displayed identity is not safely expressible as a concrete candidate structure.")
        m["stereochemistry"]["geometry_status"] = inferred(geometry if geometry else "unresolved", ev, unresolved_reason or "Geometry cannot be safely encoded from the available drawing.")
        m["stereochemistry"]["representation_status"] = inferred("unresolved", ev, unresolved_reason or "No concrete SMILES asserted.")
    if axial is not None:
        m["stereochemistry"]["axial"] = f(axial, axial, ev)
    if geometry == "specified":
        m["stereochemistry"]["double_bonds"] = f("E-isomer", [{"source_assignment": "E", "encoding": "not encoded in candidate SMILES"}], ev)
    if assignment is not None:
        m["stereochemistry"]["assignment"] = f(assignment, assignment, ev)
    m["identity_note"] = f(unresolved_reason, unresolved_reason, ev) if unresolved_reason else f()
    molecules[key] = m
    return mid


def mol(key):
    return molecules[key]["molecule_id"]


# Common reagents, solvents and work-up materials.
common = [
    ("li2co3", "Li2CO3", "lithium carbonate", "[Li+].[Li+].[O-]C([O-])=O"),
    ("k2co3", "K2CO3", "potassium carbonate", "[K+].[K+].[O-]C([O-])=O"),
    ("na2co3", "Na2CO3", "sodium carbonate", "[Na+].[Na+].[O-]C([O-])=O"),
    ("cs2co3", "Cs2CO3", "caesium carbonate", "[Cs+].[Cs+].[O-]C([O-])=O"),
    ("et3n", "Et3N", "triethylamine", "CCN(CC)CC"),
    ("pyridine", "pyridine", "pyridine", "n1ccccc1"),
    ("thf", "THF", "tetrahydrofuran", "C1CCOC1"),
    ("acetonitrile", "CH3CN", "acetonitrile", "CC#N"),
    ("dioxane", "1,4-dioxane", "1,4-dioxane", "O1CCOCC1"),
    ("benzene", "benzene", "benzene", "c1ccccc1"),
    ("chloroform", "CHCl3", "chloroform", "ClC(Cl)Cl"),
    ("dcm", "DCM", "dichloromethane", "ClCCl"),
    ("dce", "DCE", "1,2-dichloroethane", "ClCCCl"),
    ("cdcl3", "CDCl3", "deuterated chloroform", "[2H]C(Cl)(Cl)Cl"),
    ("dmf", "DMF", "N,N-dimethylformamide", "CN(C)C=O"),
    ("toluene", "toluene", "toluene", "Cc1ccccc1"),
    ("water", "water", "water", "O"),
    ("etoac", "ethyl acetate", "ethyl acetate", "CCOC(C)=O"),
    ("na2so4", "Na2SO4", "anhydrous sodium sulfate", "[Na+].[Na+].[O-]S(=O)(=O)[O-]"),
    ("silica", "silica gel", "silica gel", None),
]
for key, label, name, smi in common:
    add_molecule(key, label, name, smi, [E[2], E[3], E[4]], unresolved_reason="Named formulation/material; exact molecular composition is not represented as a single concrete SMILES." if smi is None else None)
add_molecule("brine", "brine", "brine solution", None, [E[4]], unresolved_reason="Named aqueous mixture; exact composition/concentration is unreported.")
add_molecule("oxone", "Oxone", "Oxone", None, [E[3], E[4]], unresolved_reason="Named triple-salt oxidant formulation; exact batch composition is not expanded from the article.")
add_molecule("tmsn3", "TMS-N3", "trimethylsilyl azide", "C[Si](C)(C)N=[N+]=[N-]", [E[1], E[3]])


ACYL = {
    "p_anisyl": "COc1ccc(C(=O)X)cc1",
    "p_tolyl": "Cc1ccc(C(=O)X)cc1",
    "o_anisyl": "COc1ccccc1C(=O)X",
    "benzyl": "c1ccccc1CC(=O)X",
    "phenethyl": "c1ccccc1CCC(=O)X",
    "phenyl": "O=C(X)c1ccccc1",
    "o_chlorophenyl": "O=C(X)c1ccccc1Cl",
    "p_fluorophenyl": "O=C(X)c1ccc(F)cc1",
    "p_bromophenyl": "O=C(X)c1ccc(Br)cc1",
    "dichlorophenyl": "O=C(X)c1c(Cl)cccc1Cl",
    "naphthyl": "O=C(X)c1cccc2ccccc12",
    "hexyl": "CCCCCCC(=O)X",
    "o_tolyl": "O=C(X)c1ccccc1C",
}


def core(acyl, r2, xgroup=None, allene=False):
    if allene:
        tail = "C=C=C" + ({"Me": "C", "Ph": "c1ccccc1", "Et": "CC", "H": ""}[r2])
    else:
        branch = {"Me": "(C)", "Ph": "(c1ccccc1)", "Et": "(CC)", "H": ""}[r2]
        tail = "C=C" + branch + (xgroup or "")
    return ACYL[acyl].replace("X", tail)


def allene_mid(acyl, r2, ev):
    key = f"allene_{acyl}_{r2}"
    if key not in molecules:
        ch = "unspecified" if r2 != "H" else "achiral"
        axial = [{"axis": "C=C=C", "configuration": "unspecified in source"}] if r2 != "H" else None
        add_molecule(key, "1", f"allenic ketone ({acyl}, R2={r2})", core(acyl, r2, allene=True), ev,
                     chirality=ch, geometry="not_applicable", axial=axial,
                     unresolved_reason="The substituted allene can possess axial chirality, but the paper does not assign it.")
    return mol(key)


def add_named(key, label, name, smiles, ev):
    return add_molecule(key, label, name, smiles, ev)


add_named("p_thiocresol", "2", "p-thiocresol", "Cc1ccc(S)cc1", [E[1], E[2]])


def base_reaction(rid, ev, record_type, categories, location, name, transformation):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, ev)
    r["location"] = f(location, location, ev)
    r["reaction_name"] = f(name, name, ev)
    r["transformation_type"] = inferred(transformation, ev, "Concise transformation description derived from the displayed reactant/product connectivity.")
    r["review"]["unresolved"] = []
    return r


def participant(pid, mid, role, ev, stage="reaction", equiv=None, volume=None, eqref=None):
    p = st.new_participant(pid, mid, role, ev, stage=stage)
    if role == "product":
        p["product_role"] = f("main product", "main_product", ev)
    if equiv is not None:
        p["charges"]["equivalents"] = f(f"{equiv} equiv.", q(equiv, "equiv"), ev)
        p["equivalence_reference"] = f(eqref, eqref, ev)
    if volume is not None:
        p["charges"]["volume"] = f(f"{volume} mL", q(volume, "mL"), ev)
    return p


def step(step_id, ev, stage, pids, operation, temp=None, time=None, atmosphere=None, inherited=False):
    s = st.new_step(step_id, ev, stage=stage)
    s["participant_ids"] = list(pids)
    method = "inherited_general_condition" if inherited else "original_extraction"
    s["operation"] = f(operation, operation, ev, method=method)
    if temp is not None:
        raw, val = temp
        s["temperature"] = f(raw, val, ev, method=method)
    if time is not None:
        raw, val = time
        s["time"] = f(raw, val, ev, method=method)
    if atmosphere is not None:
        s["atmosphere"] = f(atmosphere, atmosphere, ev, method=method)
    return s


def set_classification(r, kind, ev):
    patterns = {
        "thio": [("C", "S", 0, 1, "formed"), ("C", "H", 0, 1, "hydrogen_gain"), ("C", "C", 2, 1, "order_decreased")],
        "ether": [("C", "O", 0, 1, "formed"), ("C", "H", 0, 1, "hydrogen_gain"), ("C", "C", 2, 1, "order_decreased")],
        "azide": [("C", "N", 0, 1, "formed"), ("C", "H", 0, 1, "hydrogen_gain"), ("C", "C", 2, 1, "order_decreased")],
        "sulfoxide": [("O", "S", 0, 2, "formed")],
        "sulfone": [("O", "S", 0, 2, "formed"), ("O", "S", 0, 2, "formed")],
    }
    if kind not in patterns:
        r["classification"]["coverage"] = inferred("unresolved", ev, "The displayed cascade/rearrangement does not support a defensible complete atom mapping from these pages alone.")
        r["classification"]["limitations"] = ["No atom mapping asserted; bond-level classification remains unresolved."]
        return
    changes = []
    labels = {}
    for i, (a, b, before, after, ctype) in enumerate(patterns[kind], 1):
        cid = f"b{i}"
        ch = {"change_id": cid, "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
              "before_order": before, "after_order": after, "change_type": ctype}
        if ctype.startswith("hydrogen_"):
            ch["hydrogen_count"] = 1
        changes.append(inferred(ch, ev, "Net bond change derived from the generic displayed transformation; atom identities are not fully mapped."))
        if ctype in {"formed", "hydrogen_gain", "order_increased"}:
            lab = st.bond_label(a, b, after)
            labels.setdefault(lab["code"], dict(lab, supporting_change_ids=[]))["supporting_change_ids"].append(cid)
    r["classification"]["bond_changes"] = changes
    r["classification"]["labels"] = inferred(list(labels.values()), ev, "Taxonomy labels calculated from the explicitly listed partial net bond changes.")
    r["classification"]["coverage"] = inferred("partial", ev, "Only defensible net changes from the displayed generic transformation are recorded; full atom/source mapping is not claimed.")
    r["classification"]["limitations"] = ["Partial generic mapping only.", "Implicit proton source and coproducts are unresolved."]


def yield_outcome(oid, product_mids, y, ev, scope):
    o = st.new_outcome(oid, "yield", product_mids, ev)
    o["product_role"] = f("desired product", "main_product", ev)
    o["value"] = f(f"{y}%", q(y, "%"), ev)
    o["scope"] = f(scope, scope, ev)
    return o


def add_gp_a(r, ev, allene, nuc, base="li2co3", time_value="4-7 h"):
    inherited = [E[4], T[4]]
    pa = participant("allene", allene, "reactant", ev, equiv=1.0, eqref=allene)
    pn = participant("nucleophile", nuc, "reactant", ev, equiv=1.1, eqref=allene)
    parts = [pa, pn]
    if base:
        parts.append(participant("base", mol(base), "reagent", ev, equiv=1.1, eqref=allene))
    parts.append(participant("solvent", mol("thf"), "solvent", inherited, volume=5.0))
    parts += [participant("water", mol("water"), "workup_agent", inherited, stage="workup"),
              participant("etoac", mol("etoac"), "workup_agent", inherited, stage="workup", volume=20.0),
              participant("brine", mol("brine"), "workup_agent", inherited, stage="workup"),
              participant("na2so4", mol("na2so4"), "workup_agent", inherited, stage="workup"),
              participant("silica", mol("silica"), "purification_agent", inherited, stage="purification")]
    r["participants"] = parts
    r["reaction_scale"] = f("1 equiv. allenic ketone; absolute amount unreported", {"reference_molecule_id": allene, "equivalents": q(1.0, "equiv")}, ev)
    r["solvent_conditions"]["composition"] = f("THF", [mol("thf")], inherited, method="inherited_general_condition")
    r["solvent_conditions"]["total_volume"] = f("5 mL", q(5.0, "mL"), inherited, method="inherited_general_condition")
    tval = {"minimum": 4, "maximum": 7, "unit": "h", "raw_text": "4-7 h"}
    r["steps"] = [
        step("reaction", inherited, "reaction", [p["participant_id"] for p in parts if p["stage"] == "reaction"],
             "stir in an oven-dried 10 mL single-neck round-bottom flask; monitor by TLC",
             ("35 °C", q(35, "°C")), (time_value, tval), "open atmosphere", inherited=True),
        step("quench_extract", inherited, "workup", ["water", "etoac"], "quench with water and extract with ethyl acetate (3 × 20 mL)", inherited=True),
        step("wash_dry", inherited, "workup", ["brine", "na2so4"], "wash organic phase with brine and dry over anhydrous sodium sulfate", inherited=True),
        step("evaporate", inherited, "workup", [], "evaporate solvent under reduced pressure", inherited=True),
        step("chromatography", inherited, "purification", ["silica"], "purify by silica-gel column chromatography; solvent combination not specified", inherited=True),
    ]
    r["procedure"]["original_text"] = f("General Procedure (A) for Vinyl α,β-Unsaturated Ketone", "General Procedure A", inherited, method="inherited_general_condition")
    r["procedure"]["general_procedure_reference"] = f("General Procedure (A)", "General Procedure A", inherited, method="inherited_general_condition")
    r["procedure"]["addition_order"] = f("Allenic ketone in THF, then Li2CO3 (except azide), then nucleophile", ["allene in THF", "base except azide", "nucleophile"], inherited, method="inherited_general_condition")
    r["review"]["unresolved"] += ["Yield measurement basis is not stated in the provided main-text pages.", "Absolute reaction scale is unreported; SI is not supplied.", "Chromatography solvent combination is unspecified."]


def add_gp_b(r, ev, allene, aldehyde):
    inherited = [E[4], T[4]]
    parts = [participant("allene", allene, "reactant", ev, equiv=1.0, eqref=allene),
             participant("aldehyde", aldehyde, "reactant", ev, equiv=1.1, eqref=allene),
             participant("base", mol("li2co3"), "reagent", ev, equiv=1.2, eqref=allene),
             participant("oxidant", mol("oxone"), "reagent", ev, equiv=1.2, eqref=allene),
             participant("solvent", mol("dmf"), "solvent", inherited, volume=5.0),
             participant("water", mol("water"), "workup_agent", inherited, stage="workup"),
             participant("etoac", mol("etoac"), "workup_agent", inherited, stage="workup", volume=20.0),
             participant("brine", mol("brine"), "workup_agent", inherited, stage="workup"),
             participant("na2so4", mol("na2so4"), "workup_agent", inherited, stage="workup"),
             participant("silica", mol("silica"), "purification_agent", inherited, stage="purification")]
    r["participants"] = parts
    r["reaction_scale"] = f("1 equiv. allenic ketone; absolute amount unreported", {"reference_molecule_id": allene, "equivalents": q(1.0, "equiv")}, ev)
    r["solvent_conditions"]["composition"] = f("DMF", [mol("dmf")], ev)
    r["solvent_conditions"]["total_volume"] = f("5 mL", q(5.0, "mL"), inherited, method="inherited_general_condition")
    r["steps"] = [
        step("reaction", list(dict.fromkeys(ev + inherited)), "reaction", ["allene", "aldehyde", "base", "oxidant", "solvent"], "stir and monitor by TLC", ("room temperature (25 °C)", q(25, "°C")), ("12 h", q(12, "h")), "open atmosphere"),
        step("quench_extract", inherited, "workup", ["water", "etoac"], "quench with water and extract with ethyl acetate (3 × 20 mL)", inherited=True),
        step("wash_dry", inherited, "workup", ["brine", "na2so4"], "wash organic phase with brine and dry over anhydrous sodium sulfate", inherited=True),
        step("evaporate", inherited, "workup", [], "evaporate solvent under reduced pressure", inherited=True),
        step("chromatography", inherited, "purification", ["silica"], "purify by silica-gel column chromatography; solvent combination not specified", inherited=True),
    ]
    r["procedure"]["original_text"] = f("General Procedure (B) for Benzodioxin Synthesis", "General Procedure B", inherited, method="inherited_general_condition")
    r["procedure"]["general_procedure_reference"] = f("General Procedure (B)", "General Procedure B", inherited, method="inherited_general_condition")
    r["procedure"]["addition_order"] = f("Allenic ketone in DMF, followed by Li2CO3, substituted ortho-hydroxybenzaldehyde and oxone", ["allene in DMF", "Li2CO3", "aldehyde", "oxone"], inherited, method="inherited_general_condition")
    r["review"]["unresolved"] += ["Yield measurement basis is not stated in the provided main-text pages.", "Absolute reaction scale is unreported; SI is not supplied.", "Chromatography solvent combination is unspecified.", "Cascade atom mapping/classification is unresolved."]


reactions = []

# Four current-paper generic transformations.  Scheme 2 on page 1 and the table
# headers on pages 2--3 are duplicate presentations of these general arrows.
for rid, pages, name, trans, ckind in [
    ("scheme_thioether", [1, 2], "synthesis of (E)-β-vinyl thioether", "nucleophilic hydrothiolation of an allenic ketone", "thio"),
    ("scheme_ether", [1, 3], "synthesis of (E)-β-vinyl ether", "nucleophilic hydroalkoxylation of an allenic ketone", "ether"),
    ("scheme_azide", [1, 3], "synthesis of (E)-β-vinyl azide", "hydroazidation of an allenic ketone", "azide"),
    ("scheme_benzodioxin", [3], "benzodioxin skeleton synthesis", "oxidative annulation cascade", "unresolved"),
]:
    ev = [E[p] for p in pages]
    r = base_reaction(rid, ev, "generic_scheme", [], {"pages": pages, "region": "Scheme 2/table-header generic arrow"}, name, trans)
    r["review"]["unresolved"] = ["Generic R-group structures are not expanded as concrete molecules."]
    set_classification(r, ckind, ev)
    reactions.append(r)


# Table 1 optimization: 17 separate experiments.
opt = [
    (1, "acetonitrile", "k2co3", "r.t.", 60, "60:40"), (2, "thf", "k2co3", "r.t.", 70, "65:35"),
    (3, "thf", "na2co3", "r.t.", 75, "65:35"), (4, "thf", "na2co3", 40, 80, "85:15"),
    (5, "thf", "cs2co3", 40, 70, "85:15"), (6, "thf", "li2co3", 40, 85, "90:10"),
    (7, "thf", "li2co3", 35, 96, "100:0"), (8, "thf", "li2co3", "r.t.", 96, "55:45"),
    (9, "dioxane", "li2co3", 35, 85, "100:0"), (10, "benzene", "li2co3", 35, 70, "65:35"),
    (11, "chloroform", "li2co3", 35, 78, "100:0"), (12, "dcm", "li2co3", 35, 75, "100:0"),
    (13, "dce", "li2co3", 35, 80, "100:0"), (14, "thf", "et3n", 35, 75, "83:17"),
    (15, "thf", "pyridine", 35, 70, "65:35"), (16, "cdcl3", "et3n", 35, 72, "80:20"),
    (17, "cdcl3", "li2co3", 35, 80, "100:0"),
]
opt_ev = [E[2], T[2]]
allene_1f = allene_mid("benzyl", "Me", [E[1], E[2]])
prod_2f_e = add_molecule("product_2f_E", "2f", "(E)-2f β-vinyl thioether", core("benzyl", "Me", "Sc1ccc(C)cc1"), [E[1], E[2], E[3]], geometry="specified", assignment="Relative alkene geometry designated E by the authors; Figure 2 gives nOe support.")
prod_2f_z = add_molecule("product_2f_Z", "2f′", "(Z)-2f′ β-vinyl thioether", core("benzyl", "Me", "Sc1ccc(C)cc1"), [E[2]], geometry="specified", assignment="Alkene geometry designated Z in Table 1.")
for n, solvent, base, temp, y, ratio in opt:
    rid = f"table1_entry_{n:02d}"
    r = base_reaction(rid, opt_ev, "specific_experiment", ["optimization"], {"page": 2, "table": 1, "entry": n}, "optimization of 2f/2f′ formation", "hydrothiolation of allenic ketone")
    r["general_scheme_ids"] = [f"{DOC}:scheme_thioether"]
    r["participants"] = [participant("allene_1f", allene_1f, "reactant", opt_ev, equiv=1.0, eqref=allene_1f),
                         participant("thiol_2", mol("p_thiocresol"), "reactant", opt_ev, equiv=1.1, eqref=allene_1f),
                         participant("base", mol(base), "reagent", opt_ev, equiv=1.1, eqref=allene_1f),
                         participant("solvent", mol(solvent), "solvent", opt_ev, volume=5.0),
                         participant("product_E", prod_2f_e, "product", opt_ev),
                         participant("product_Z", prod_2f_z, "product", opt_ev)]
    r["reaction_scale"] = f("1f (1 equiv.); absolute amount unreported", {"reference_molecule_id": allene_1f, "equivalents": q(1.0, "equiv")}, opt_ev)
    r["solvent_conditions"]["composition"] = f(molecules[solvent]["names"]["original"]["normalized_value"], [mol(solvent)], opt_ev)
    r["solvent_conditions"]["total_volume"] = f("dry solvent (5 mL)", q(5.0, "mL"), opt_ev)
    tv = ("r.t.", "r.t.") if temp == "r.t." else (f"{temp} °C", q(temp, "°C"))
    r["steps"] = [step("reaction", opt_ev, "reaction", ["allene_1f", "thiol_2", "base", "solvent"], "react under the Table 1 conditions", tv, ("5 h", q(5, "h")), "air atmosphere")]
    r["outcomes"] = [yield_outcome("yield_total", [prod_2f_e, prod_2f_z], y, opt_ev, {"table": 1, "entry": n, "target": "2f + 2f′"})]
    sel = st.new_outcome("selectivity_E_Z", "selectivity", [prod_2f_e, prod_2f_z], opt_ev)
    sel["value"] = f(ratio, ratio, opt_ev)
    sel["measurement_basis"] = f("E/Z is determined by 1H-NMR", "1H NMR", opt_ev)
    sel["scope"] = f("2f:2f′", {"metric": "E/Z ratio", "order": [prod_2f_e, prod_2f_z]}, opt_ev)
    r["outcomes"].append(sel)
    r["procedure"]["original_text"] = f("Table 1 footnote: 1f (1 equiv.), 2 (1.1 equiv.), base (1.1 equiv.), dry solvent (5 mL), 5 h under air", "Table 1 footnote conditions", opt_ev)
    r["review"]["unresolved"] = ["Yield measurement basis is unreported; only the E/Z ratio method is specified as 1H NMR.", "Absolute reaction scale and work-up are not reported for the optimization row."]
    set_classification(r, "thio", opt_ev)
    reactions.append(r)


# Scope tables.  Candidate connectivity is included when the compact drawing is
# unambiguous; otherwise the product/nucleophile remains explicitly unresolved.
thio_rows = [
    ("2a", 96, "p_anisyl", "Me", None, None, "aliphatic thiol shown for 2a"),
    ("2b", 90, "p_tolyl", "Me", "SCCO", "OCCS", "2-mercaptoethanol"),
    ("2c", 94, "o_anisyl", "Me", "SCc1ccccc1", "c1ccccc1CS", "benzyl mercaptan"),
    ("2d", 88, "p_anisyl", "Me", "SCCO", "OCCS", "2-mercaptoethanol"),
    ("2e", 95, "phenyl", "Me", "SC1CCCCC1", "SC1CCCCC1", "cyclohexanethiol"),
    ("2f", 96, "benzyl", "Me", "Sc1ccc(C)cc1", "Cc1ccc(S)cc1", "p-thiocresol"),
    ("2g", 93, "benzyl", "Me", "Sc1ccccn1", "Sc1ccccn1", "2-mercaptopyridine"),
    ("2h", 93, "phenethyl", "H", "Sc1ccccc1", "Sc1ccccc1", "thiophenol"),
    ("2i", 90, "benzyl", "Me", "Sc1ccc(OC)cc1", "COc1ccc(S)cc1", "p-methoxythiophenol"),
    ("2j", 90, "benzyl", "Me", "Sc1ccc(Br)cc1", "Brc1ccc(S)cc1", "p-bromothiophenol"),
    ("2k", 70, "phenethyl", "Me", "Sc1ccc(N)cc1", "Nc1ccc(S)cc1", "p-aminothiophenol"),
    ("2l", 88, "benzyl", "Me", "Sc1ccc(F)cc1", "Fc1ccc(S)cc1", "p-fluorothiophenol"),
    ("2m", 94, "p_tolyl", "Me", "Sc1ccc(C)cc1", "Cc1ccc(S)cc1", "p-thiocresol"),
    ("2n", 92, "p_tolyl", "Me", "Sc1ccc(Br)cc1", "Brc1ccc(S)cc1", "p-bromothiophenol"),
    ("2o", 90, "dichlorophenyl", "Me", "Sc1ccc(F)cc1", "Fc1ccc(S)cc1", "p-fluorothiophenol"),
    ("2p", 88, "o_anisyl", "Me", "Sc1ccc(OC)cc1", "COc1ccc(S)cc1", "p-methoxythiophenol"),
    ("2q", 90, "dichlorophenyl", "Me", "Sc1ccc(Cl)cc1", "Clc1ccc(S)cc1", "p-chlorothiophenol"),
    ("2r", 92, "p_bromophenyl", "Me", "Sc1ccc(C)cc1", "Cc1ccc(S)cc1", "p-thiocresol"),
    ("2s", 92, "o_anisyl", "Me", "Sc1ccccn1", "Sc1ccccn1", "2-mercaptopyridine"),
    ("2t", 89, "dichlorophenyl", "Me", "Sc1ccc(Br)cc1", "Brc1ccc(S)cc1", "p-bromothiophenol"),
    ("2u", 88, "p_tolyl", "Ph", "Sc1ccc(C)cc1", "Cc1ccc(S)cc1", "p-thiocresol"),
    ("2v", 90, "hexyl", "Me", "Sc1ccc(C)cc1", "Cc1ccc(S)cc1", "p-thiocresol"),
    ("2w", 85, "hexyl", "Me", "SCCO", "OCCS", "2-mercaptoethanol"),
    ("2x", 86, "benzyl", "Me", "Sc1ccccc1N", "Nc1ccccc1S", "o-aminothiophenol"),
]
for label, y, acyl, r2, xgroup, nuc_smi, nuc_name in thio_rows:
    ev = [E[2], E[4]]
    allene = allene_mid(acyl, r2, ev)
    nkey = f"thiol_{slug(nuc_name)}"
    nuc = add_molecule(nkey, "R-SH", nuc_name, nuc_smi, [E[2]], unresolved_reason="Compact side-chain drawing is not safely resolved to a unique concrete thiol." if nuc_smi is None else None)
    if label == "2f":
        product = prod_2f_e
    else:
        psmiles = core(acyl, r2, xgroup) if xgroup else None
        product = add_molecule(f"product_{label}_E", label, f"(E)-β-vinyl thioether {label}", psmiles, [E[2]] + ([E[3]] if label == "2s" else []), geometry="specified", assignment="E-isomer stated in Table 2; Figure 2 provides nOe support." if label == "2s" else "E-isomer stated in Table 2.", unresolved_reason="The compact substituent drawing is not safely resolved; the source label and E assignment are retained.")
    r = base_reaction(f"table2_{label}", ev, "specific_experiment", ["substrate_scope"], {"page": 2, "table": 2, "entry": label}, "synthesis of (E)-β-vinyl thioether", "hydrothiolation of an allenic ketone")
    r["general_scheme_ids"] = [f"{DOC}:scheme_thioether"]
    add_gp_a(r, ev, allene, nuc)
    r["participants"].append(participant("product", product, "product", [E[2]]))
    r["outcomes"] = [yield_outcome("reported_yield", [product], y, [E[2]], {"table": 2, "entry": label})]
    set_classification(r, "thio", [E[2]])
    reactions.append(r)


ether_rows = [
    ("3a", 96, "p_tolyl", "Me", "OCC", "CCO", "ethanol"),
    ("3b", 90, "o_chlorophenyl", "Me", "OCC", "CCO", "ethanol"),
    ("3c", 88, "o_chlorophenyl", "Me", "OC(C)C", "CC(C)O", "2-propanol"),
    ("3d", 92, "o_chlorophenyl", "Me", "OC", "CO", "methanol"),
    ("3e", 94, "phenethyl", "Me", "OC", "CO", "methanol"),
    ("3f", 92, "p_tolyl", "Me", None, None, "formylphenol isomer shown for 3f"),
    ("3g", 92, "p_tolyl", "Me", "Oc1ccc([N+](=O)[O-])cc1", "Oc1ccc([N+](=O)[O-])cc1", "p-nitrophenol"),
    ("3h", 94, "p_tolyl", "Me", "Oc1ccccc1[N+](=O)[O-]", "Oc1ccccc1[N+](=O)[O-]", "o-nitrophenol"),
    ("3i", 90, "p_tolyl", "Me", "Oc1ccccn1", "Oc1ccccn1", "2-hydroxypyridine"),
    ("3j", 94, "p_tolyl", "Me", None, None, "methoxy-formyl phenol shown for 3j"),
    ("3k", 95, "p_fluorophenyl", "Me", "Oc1ccccc1C=O", "O=Cc1ccccc1O", "salicylaldehyde"),
    ("3l", 91, "p_fluorophenyl", "Me", None, None, "formylphenol isomer shown for 3l"),
    ("3m", 92, "p_tolyl", "Me", "Oc1ccccc1C=O", "O=Cc1ccccc1O", "salicylaldehyde"),
    ("3n", 90, "p_tolyl", "Me", "Oc1ccc(C=O)cc1", "O=Cc1ccc(O)cc1", "p-hydroxybenzaldehyde"),
    ("3o", 86, "hexyl", "Me", "OC", "CO", "methanol"),
]
for label, y, acyl, r2, xgroup, nuc_smi, nuc_name in ether_rows:
    ev = [E[3], E[4]]
    allene = allene_mid(acyl, r2, ev)
    nkey = f"alcohol_{slug(nuc_name)}"
    nuc = add_molecule(nkey, "R-OH", nuc_name, nuc_smi, [E[3]], unresolved_reason="The compact aromatic substitution pattern cannot be safely resolved from the provided drawing." if nuc_smi is None else None)
    psmiles = core(acyl, r2, xgroup) if xgroup else None
    product = add_molecule(f"product_{label}_E", label, f"(E)-β-vinyl ether {label}", psmiles, [E[3]], geometry="specified", assignment="E-isomer stated in Table 3; 3d is additionally shown in the nOe figure." if label == "3d" else "E-isomer stated in Table 3.", unresolved_reason="The compact aromatic substitution pattern is retained by label but not forced into a concrete SMILES.")
    r = base_reaction(f"table3_{label}", ev, "specific_experiment", ["substrate_scope"], {"page": 3, "table": 3, "entry": label}, "synthesis of (E)-β-vinyl ether", "hydroalkoxylation of an allenic ketone")
    r["general_scheme_ids"] = [f"{DOC}:scheme_ether"]
    add_gp_a(r, ev, allene, nuc)
    r["participants"].append(participant("product", product, "product", [E[3]]))
    r["outcomes"] = [yield_outcome("reported_yield", [product], y, [E[3]], {"table": 3, "entry": label})]
    set_classification(r, "ether", [E[3]])
    reactions.append(r)


azide_rows = [
    ("4a", 86, "o_tolyl", "Me"), ("4b", 88, "phenyl", "Me"), ("4c", 84, "p_fluorophenyl", "Me"),
    ("4d", 85, None, "Me"), ("4e", 83, None, "Me"), ("4f", 82, "dichlorophenyl", "Me"),
    ("4g", 90, "naphthyl", "Me"), ("4h", 85, None, "Me"), ("4i", 88, None, "Me"),
    ("4j", 83, None, "Et"), ("4k", 75, "benzyl", "Me"), ("4l", 78, "phenethyl", "Me"),
]
for label, y, acyl, r2 in azide_rows:
    ev = [E[3], E[4]]
    if acyl:
        allene = allene_mid(acyl, r2, ev)
        psmiles = core(acyl, r2, "N=[N+]=[N-]")
    else:
        allene = add_molecule(f"allene_for_{label}", "1", f"allenic ketone precursor of {label}", None, [E[3]], chirality="unresolved", geometry="not_applicable", unresolved_reason="Fused/heteroaryl compact drawing is not safely resolved to a unique concrete structure.")
        psmiles = None
    product = add_molecule(f"product_{label}_E", label, f"(E)-β-vinyl azide {label}", psmiles, [E[3]] + ([E[4]] if label in {"4j", "4l"} else []), geometry="specified", assignment="E-isomer stated in Table 4; Figure 2 gives nOe support." if label in {"4j", "4l"} else "E-isomer stated in Table 4.", unresolved_reason="Fused/heteroaryl compact drawing is retained by label without a forced SMILES.")
    r = base_reaction(f"table4_{label}", ev, "specific_experiment", ["substrate_scope"], {"page": 3, "table": 4, "entry": label}, "synthesis of (E)-β-vinyl azide", "hydroazidation of an allenic ketone")
    r["general_scheme_ids"] = [f"{DOC}:scheme_azide"]
    add_gp_a(r, ev, allene, mol("tmsn3"), base=None)
    r["participants"].append(participant("product", product, "product", [E[3]]))
    r["outcomes"] = [yield_outcome("reported_yield", [product], y, [E[3]], {"table": 4, "entry": label})]
    set_classification(r, "azide", [E[3]])
    reactions.append(r)


# Table 5: exact fused product SMILES are left unresolved; reactants are retained
# where their compact identities are legible.
add_named("salicylaldehyde", "ortho-hydroxy benzaldehyde", "salicylaldehyde", "O=Cc1ccccc1O", [E[3], E[4]])
add_named("methoxy_salicylaldehyde", "methoxy ortho-hydroxy benzaldehyde", "methoxy-substituted salicylaldehyde shown for 5c", None, [E[3]])
add_named("naphthol_aldehyde", "naphthol aldehyde", "naphthol-derived aldehyde shown for 5d", None, [E[3]])
add_named("methyl_salicylaldehyde", "methyl ortho-hydroxy benzaldehyde", "methyl-substituted salicylaldehyde shown for 5e", None, [E[3]])
benz_rows = [
    ("5a", 94, "p_fluorophenyl", "salicylaldehyde"), ("5b", 91, "p_tolyl", "salicylaldehyde"),
    ("5c", 92, "p_tolyl", "methoxy_salicylaldehyde"), ("5d", 90, "p_anisyl", "naphthol_aldehyde"),
    ("5e", 93, "naphthyl", "methyl_salicylaldehyde"),
]
for label, y, acyl, aldehyde_key in benz_rows:
    ev = [E[3], E[4]]
    allene = allene_mid(acyl, "Me", ev)
    product = add_molecule(f"product_{label}", label, f"benzodioxin {label}", None, [E[3]], geometry="not_applicable", unresolved_reason="The fused benzodioxin connectivity/substitution is preserved by the source label and image; no unverified concrete SMILES is asserted.")
    r = base_reaction(f"table5_{label}", ev, "specific_experiment", ["substrate_scope"], {"page": 3, "table": 5, "entry": label}, "benzodioxin skeleton synthesis", "oxidative annulation cascade")
    r["general_scheme_ids"] = [f"{DOC}:scheme_benzodioxin"]
    add_gp_b(r, ev, allene, mol(aldehyde_key))
    r["participants"].append(participant("product", product, "product", [E[3]]))
    r["outcomes"] = [yield_outcome("reported_yield", [product], y, [E[3]], {"table": 5, "entry": label})]
    set_classification(r, "unresolved", [E[3]])
    reactions.append(r)


# Scheme 4 post-functionalizations.  The 3m -> 5b route is a distinct experiment
# from the Table 5 one-pot formation (same product identity, different inputs).
scheme4_ev = [E[4], T[4]]
post = [
    ("scheme4_6a", "2m", "product_6a", "β-vinyl sulfoxide 6a", 90, 1.1, "sulfoxide"),
    ("scheme4_6b", "2n", "product_6b", "β-vinyl sulfoxide 6b", 85, 1.1, "sulfoxide"),
    ("scheme4_7", "2m", "product_7", "β-vinyl sulfone 7", 90, 2.2, "sulfone"),
]
for rid, source_label, pkey, pname, y, oxeq, ckind in post:
    source_mid = mol(f"product_{source_label}_E")
    product = add_molecule(pkey, pkey.replace("product_", ""), pname, None, [E[4]], geometry="specified", unresolved_reason="Oxidation state and E-alkene are visible, but a full concrete structure is not duplicated from the compact scheme without source review.")
    r = base_reaction(rid, scheme4_ev, "specific_experiment", ["derivatization"], {"page": 4, "scheme": 4, "arrow": rid}, "oxidation of β-vinyl thioether", "thioether oxidation")
    r["participants"] = [participant("substrate", source_mid, "reactant", scheme4_ev, equiv=1.0, eqref=source_mid),
                         participant("oxone", mol("oxone"), "reagent", scheme4_ev, equiv=oxeq, eqref=source_mid),
                         participant("dmf", mol("dmf"), "solvent", scheme4_ev),
                         participant("product", product, "product", scheme4_ev)]
    r["reaction_scale"] = f("substrate shown without absolute charge", {"reference_molecule_id": source_mid, "equivalents": q(1.0, "equiv")}, scheme4_ev)
    r["solvent_conditions"]["composition"] = f("DMF", [mol("dmf")], scheme4_ev)
    r["steps"] = [step("reaction", scheme4_ev, "reaction", ["substrate", "oxone", "dmf"], "treat with Oxone in DMF", ("r.t. (25 °C)", q(25, "°C")), ("12 h", q(12, "h")), "open-air")]
    r["outcomes"] = [yield_outcome("reported_yield", [product], y, scheme4_ev, {"scheme": 4, "product": pkey.replace("product_", "")})]
    r["review"]["unresolved"] = ["Yield measurement basis, absolute scale, solvent volume, work-up and purification are unreported in the provided pages."]
    set_classification(r, ckind, scheme4_ev)
    reactions.append(r)

source_3m = mol("product_3m_E")
product_5b = mol("product_5b")
r = base_reaction("scheme4_3m_to_5b", scheme4_ev, "specific_experiment", ["derivatization"], {"page": 4, "scheme": 4, "arrow": "3m to 5b"}, "oxidative cyclization of β-vinyl ether 3m", "oxidative cyclization")
r["participants"] = [participant("substrate", source_3m, "reactant", scheme4_ev, equiv=1.0, eqref=source_3m),
                     participant("oxone", mol("oxone"), "reagent", scheme4_ev, equiv=1.2, eqref=source_3m),
                     participant("dmf", mol("dmf"), "solvent", scheme4_ev),
                     participant("product", product_5b, "product", scheme4_ev)]
r["reaction_scale"] = f("3m shown without absolute charge", {"reference_molecule_id": source_3m, "equivalents": q(1.0, "equiv")}, scheme4_ev)
r["solvent_conditions"]["composition"] = f("DMF", [mol("dmf")], scheme4_ev)
r["steps"] = [step("reaction", scheme4_ev, "reaction", ["substrate", "oxone", "dmf"], "treat 3m with Oxone in DMF", ("r.t. (25 °C)", q(25, "°C")), ("12 h", q(12, "h")), "open-air")]
r["outcomes"] = [yield_outcome("reported_yield", [product_5b], 94, scheme4_ev, {"scheme": 4, "route": "3m to 5b"})]
r["review"]["unresolved"] = ["Yield measurement basis, absolute scale, solvent volume, work-up and purification are unreported.", "Cyclization atom mapping/classification is unresolved."]
set_classification(r, "unresolved", scheme4_ev)
reactions.append(r)

for label, y, substituent, precursor_label in [("8a", 75, "2Me", "4a"), ("8b", 74, "4F", "4c")]:
    # Scheme 4's R=2Me and R=4F drawings match the concrete 4a and 4c
    # connectivity shown directly above in Table 4; the cross-page association
    # remains a machine inference and does not promote review status.
    substrate = mol(f"product_{precursor_label}_E")
    product = add_molecule(f"product_{label}", label, f"2,5-disubstituted oxazole {label}", None, [E[4]], geometry="not_applicable", unresolved_reason="Compact R-group product drawing is preserved without an unverified concrete SMILES.")
    r = base_reaction(f"scheme4_{label}", scheme4_ev, "specific_experiment", ["derivatization"], {"page": 4, "scheme": 4, "product": label}, "thermolysis of β-vinyl azide", "2H-azirine-mediated oxazole formation")
    r["participants"] = [participant("azide", substrate, "reactant", scheme4_ev), participant("toluene", mol("toluene"), "solvent", scheme4_ev), participant("product", product, "product", scheme4_ev)]
    r["solvent_conditions"]["composition"] = f("toluene", [mol("toluene")], scheme4_ev)
    r["steps"] = [step("thermolysis", scheme4_ev, "reaction", ["azide", "toluene"], "heat β-vinyl azide via a proposed 2H-azirine", ("120 °C", q(120, "°C")), None, "nitrogen")]
    r["outcomes"] = [yield_outcome("reported_yield", [product], y, scheme4_ev, {"scheme": 4, "product": label, "R": substituent})]
    r["notes"].append(inferred({"scheme4_R": substituent, "matched_table4_precursor": precursor_label}, [E[3], E[4]], "Cross-page structure match between the Scheme 4 R substituent and Table 4 product drawing."))
    r["review"]["unresolved"] = ["Reaction time, absolute scale, solvent volume, work-up, purification and yield basis are unreported.", "Product structure and rearrangement atom mapping remain unresolved."]
    set_classification(r, "unresolved", scheme4_ev)
    reactions.append(r)


bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": f"{DOC}:context:scheme1_prior_work", "document_id": DOC, "kind": "literature_context", "evidence_ids": [E[1]], "description": "Scheme 1 summarizes Sugita's prior allenic-ketone/thiol work and four possible isomers; it is not a current-paper experiment."},
    {"context_id": f"{DOC}:context:scheme3_prior_routes", "document_id": DOC, "kind": "literature_context", "evidence_ids": [E[2]], "description": "Scheme 3 collects cited prior benzodioxin strategies (refs. 9-16); arrows were inventoried as prior work, not current experiments."},
    {"context_id": f"{DOC}:context:figure1_bioactive_structures", "document_id": DOC, "kind": "background_structures", "evidence_ids": [E[1]], "description": "Figure 1 contains biologically important benzodioxin examples, not reaction experiments."},
    {"context_id": f"{DOC}:context:figure2_noe", "document_id": DOC, "kind": "stereochemical_characterization", "evidence_ids": [E[3]], "description": "Figure 2 reports multiple nOe relationships for 2f, 2s, 3d, 4j and 4l; no separate synthesis record was fabricated."},
    {"context_id": f"{DOC}:context:table5_cascade_drawing", "document_id": DOC, "kind": "author_proposed_pathway", "evidence_ids": [E[3]], "description": "Table 5 includes a stepwise ether/oxidation/cyclization depiction supporting the one-pot cascade; it is not treated as three independently measured experiments."},
    {"context_id": f"{DOC}:context:references_page5", "document_id": DOC, "kind": "bibliography", "evidence_ids": [E[5]], "description": "Page 5 contains references and the received date only; no current-paper experiment, spectrum or mass table occurs."},
]

# Save a raw model/extraction artifact separately from the formal submission.
raw_path = RUN / "raw_model/extracted_bundle_pages_0001_0005.json"
st.write(raw_path, bundle)
st.write(RUN / "submitted.json", bundle)

by_page = {p: [] for p in range(1, 6)}
for rxn in reactions:
    cited_pages = {int(e.split(":page", 1)[1].split(":", 1)[0]) for e in rxn["evidence_ids"]}
    for p in cited_pages:
        by_page[p].append(rxn["reaction_id"])

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title, abstract and bibliographic header", "finding": "True title, DOI, authors and explicit online-SI notice; SI file not supplied/associated.", "evidence_ids": [E[1], T[1]]},
            {"region": "Scheme 1", "finding": "Sugita prior-work literature example with four isomers; retained as source_context, not a current experiment.", "evidence_ids": [E[1]]},
            {"region": "Scheme 2", "finding": "Current generic thioether/ether and azide arrows; linked to generic records and later scope tables.", "evidence_ids": [E[1]]},
            {"region": "Figure 1", "finding": "Biologically important structures; non-synthetic background context.", "evidence_ids": [E[1]]},
        ], "record_ids": by_page[1], "unresolved": ["No supplied SI despite explicit availability notice.", "Scheme 2 uses R groups; individual identities come from later tables."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Scheme 3", "finding": "Cited prior benzodioxin strategies; retained as literature context.", "evidence_ids": [E[2]]},
            {"region": "Table 1 including footnotes", "finding": "All 17 optimization rows extracted separately, including solvent/base/temperature/yield/E:Z ratio, 5 h, 5 mL dry solvent, air and 1H-NMR ratio method.", "evidence_ids": [E[2], T[2]]},
            {"region": "Table 2 including generic arrow", "finding": "All 24 thioether scope entries 2a-2x extracted separately with reported yields and shared conditions.", "evidence_ids": [E[2]]},
        ], "record_ids": by_page[2], "unresolved": ["Scope-yield measurement basis and absolute scale are not stated on this page.", "Compact 2a side chain is not forced into a unique structure."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Table 3", "finding": "All 15 ether scope entries 3a-3o extracted separately.", "evidence_ids": [E[3]]},
            {"region": "Table 4", "finding": "All 12 azide scope entries 4a-4l extracted separately; no external base in the arrow.", "evidence_ids": [E[3]]},
            {"region": "Figure 2", "finding": "nOe characterization for 2f, 2s, 3d, 4j and 4l; retained as stereochemical source_context, not synthesis.", "evidence_ids": [E[3]]},
            {"region": "Table 5", "finding": "All 5 benzodioxin scope entries 5a-5e and the one-pot/cascade depiction extracted; proposed intermediate path not split into fabricated experiments.", "evidence_ids": [E[3]]},
        ], "record_ids": by_page[3], "unresolved": ["Several compact aromatic/fused substitution patterns and all benzodioxin products remain structure-unresolved.", "Scope-yield basis is unreported.", "Cascade atom mapping is unresolved."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "Scheme 4", "finding": "Five distinct post-functionalization experiments extracted: 2→6a, 2→6b, 2m→7, 3m→5b, and two azide thermolyses to 8a/8b (six product rows but five arrow families; 8a and 8b retained as separate experiments).", "evidence_ids": [E[4], T[4]]},
            {"region": "General Procedure A", "finding": "Shared Table 2-4 reaction, work-up and purification sequence inherited into 51 scope records; azide rows exclude Li2CO3 as stated.", "evidence_ids": [E[4], T[4]]},
            {"region": "General Procedure B", "finding": "Shared Table 5 reaction, work-up and purification sequence inherited into five scope records.", "evidence_ids": [E[4], T[4]]},
        ], "record_ids": by_page[4], "unresolved": ["Absolute scales and chromatography solvent combinations are absent without SI.", "Scheme 4 generally omits work-up, scale and yield basis.", "Thermolysis reaction time is unreported."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "both columns", "finding": "References [7]-[18] and received date only; no current reaction arrow, experimental row, spectrum, mass table or molecular-weight table.", "evidence_ids": [E[5], T[5]]},
        ], "record_ids": by_page[5], "unresolved": []},
    ],
}
cov = RUN / "coverage/pages_0001_0005.json"
st.write(cov, coverage)

tasks = st.read(RUN / "tasks.json")
for task in tasks:
    task["state"] = "attempted"
    task["reason"] = "Original page image inspected; reaction inventory and unresolved items are recorded in coverage/pages_0001_0005.json."
st.write(RUN / "tasks.json", tasks)

# Stereo diagnostics are generated for every concrete candidate, not used as
# proof of source correctness.
diagnostics = []
for m in bundle["molecules"]:
    smi = m["smiles"]["normalized_value"]
    if smi:
        try:
            diagnostics.append({"molecule_id": m["molecule_id"], "smiles": smi, "diagnostics": st.stereo_inspect(smi)})
        except Exception as exc:
            diagnostics.append({"molecule_id": m["molecule_id"], "smiles": smi, "error": str(exc)})
st.write(RUN / "raw_model/stereo_diagnostics_pages_0001_0005.json", diagnostics)

workflow = f"""# Workflow log — allene-stereo-v31-200-20261017

## Scope and source controls

- Document: `{DOC}`; DOI `10.1002/ejoc.202000076`.
- Title: *{TITLE}*.
- Original PDF pages 1–5 were each opened at original detail and inspected over the full page. Dense Tables 1–5 and Scheme 4 were additionally cropped from the supplied page images for legibility; crop lineage receipts are under `runs/current/raw_model/`.
- Source mapping identifies a journal PDF (`presentation_conversion: false`), not a PPTX-derived presentation.
- The prepared dataset began with zero candidates and zero experiments. The pre-existing `runs/workflow.md` concerned another document and was not used for any chemistry; this log replaces it.
- Page 1 explicitly advertises online Supporting Information, but no SI file was supplied or explicitly associated. SI-only details remain missing.

## Page inventory and decisions

- Page 1: title/abstract; prior Sugita Scheme 1; current generic Scheme 2; biological Figure 1. Scheme 1 and Figure 1 are source context, not current experiments. Three Scheme 2 transformations are represented by generic records shared with Tables 2–4.
- Page 2: prior-method Scheme 3 retained as literature context; Table 1 yielded 17 separate optimization records; Table 2 yielded 24 separate thioether scope records.
- Page 3: Table 3 yielded 15 ether records; Table 4 yielded 12 azide records; Table 5 yielded 5 benzodioxin records. Figure 2 is nOe characterization, not a synthesis record. The Table 5 stepwise depiction is an author pathway for the one-pot experiment and was not split into fictitious measurements.
- Page 4: Scheme 4 yielded five transformation families and six specific product experiments (6a, 6b, 7, 5b, 8a, 8b). The 3m→5b route is distinct from Table 5's one-pot 5b route. General Procedures A and B were merged into the corresponding scope records instead of duplicated as reactions.
- Page 5: references and received date only; no current reaction, independent spectrum, mass peak table, molecular-weight table or fragment-abundance observation.

## Record accounting and evidence

- Formal bundle contains {len(reactions)} reactions: 4 generic schemes and {len(reactions)-4} specific experiments. The specific total is 17 optimization + 24 thioether + 15 ether + 12 azide + 5 benzodioxin + 6 post-functionalization records.
- Every printed yield and every Table 1 E/Z ratio targets its explicit product molecule(s). Table 1 ratio basis is 1H NMR. Yield bases are left unreported because these pages do not explicitly say isolated/NMR/GC yield.
- General Procedure A contributes allene 1 equiv, nucleophile 1.1 equiv, Li2CO3 1.1 equiv except azides, THF 5 mL, open atmosphere, 35 °C, 4–7 h, TLC, water quench, EtOAc extraction 3 × 20 mL, brine wash, Na2SO4 drying, reduced-pressure evaporation and silica chromatography. General Procedure B contributes allene 1 equiv, aldehyde 1.1 equiv, Li2CO3/Oxone 1.2 equiv each, DMF 5 mL, open atmosphere, 25 °C, 12 h and the same work-up/purification sequence.
- Table 1 rows retain their own dry solvent (5 mL), 5 h and air conditions rather than silently inheriting the later scope work-up.
- The generated Markdown is checked after merge for reagent doses, work-up quantities and explicit unresolved notes.

## Structures, stereochemistry and classification

- Molecules are independently registered and reactions reference `molecule_id`. CAS numbers are not guessed.
- Every molecule has explicit chirality, geometry and representation states. Candidate SMILES are generated only for drawings with a defensible connectivity; compact fused/R-group or uncertain substitution patterns remain unresolved.
- Substituted allenes are marked with unassigned axial chirality where applicable; missing axial assignment is not called achiral or racemic. Scope products retain the source E designation, while candidate SMILES intentionally do not force an uncertain CIP slash encoding; representation is partial.
- Figure 2 nOe evidence is associated with 2f, 2s, 3d, 4j and 4l. No ee/er/dr values occur in these five pages. Table 1 E/Z ratios remain experiment-specific selectivity outcomes.
- `standard stereo` diagnostics are saved for every concrete candidate. Toolkit parsing is not treated as source verification.
- Hydrothiolation/hydroalkoxylation/hydroazidation classifications are partial and multi-label (C–S/C–O/C–N plus implicit C–H formation; allene C=C order decrease retained but not mislabeled as a new C–C bond). Sulfoxide/sulfone oxidations retain S=O formation. Cascade/rearrangement classifications remain unresolved without defensible atom mapping. No model votes were fabricated.

## Outputs and unresolved matters

- Reproducible builder: `runs/current/raw_model/build_pages_0001_0005.py`.
- Raw extraction bundle: `runs/current/raw_model/extracted_bundle_pages_0001_0005.json`; formal submission: `runs/current/submitted.json`.
- Coverage: `runs/current/coverage/pages_0001_0005.json`, exactly pages 1–5, with nonempty inspected regions and `semantic_accuracy_verified: false`.
- All machine records remain `needs_review`; schema validation is not chemistry certification or recall proof.
- Unresolved: absent SI; absolute scope scales; chromatography solvent combinations; yield measurement bases; Scheme 4 work-up/scale details and thermolysis time; several compact substituent/fused structures; CIP encoding of the paper's E labels; complete atom mapping for cascades/rearrangements.
"""
(ROOT / "runs/workflow.md").write_text(workflow)
print(json.dumps({"reactions": len(reactions), "molecules": len(molecules), "coverage_pages": len(coverage["pages"])}, ensure_ascii=False))
