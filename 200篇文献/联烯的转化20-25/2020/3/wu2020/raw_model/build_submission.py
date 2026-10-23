#!/usr/bin/env python3
"""Author the page-1--6 extraction bundle from inspected source images.

This script contains source transcriptions and explicit extraction decisions. It
does not alter the schema/tooling and does not certify chemical correctness.
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

RUN = ROOT / "runs/current"
DOC = "doc_439fe3cfd49f9dfc"


def ev(page, kind="text"):
    return f"{DOC}:page{page}:{kind}"


def f(raw=None, normalized=None, pages=(), method="original_extraction", derivation=None,
      status=None):
    evidence = [eid for p in pages for eid in (ev(p), ev(p, "image"))]
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation,
                   status=status)


def inferred(value, pages, why, raw=None):
    return f(raw, value, pages, method="inference", derivation=why)


def derived(value, pages, why, raw=None):
    return f(raw, value, pages, method="structure_derivation", derivation=why)


def set_fact(obj, key, raw, normalized, pages, method="original_extraction", derivation=None,
             status=None):
    obj[key] = f(raw, normalized, pages, method, derivation, status)


dataset = json.loads((RUN / "dataset.json").read_text())
for document in dataset["documents"]:
    if document["document_id"] == DOC:
        document["title"] = "Phosphine-catalyzed [3 + 2] annulation of 2-aminoacrylates with allenoates and mechanistic studies"
(RUN / "dataset.json").write_text(json.dumps(dataset, ensure_ascii=False, indent=2) + "\n")
bundle = st.empty(dataset["documents"])


def molecule(mid, label, pages, *, name=None, smiles=None, note=None,
             chiral="unresolved", geometry="unresolved", representation="unresolved",
             assignment=None, mixture=None, axial=None):
    evidence = [eid for p in pages for eid in (ev(p), ev(p, "image"))]
    m = st.new_molecule(mid, DOC, evidence)
    m["occurrences"][0]["compound_label"] = f(label, label, pages)
    m["names"]["original"] = f(name or label, name or label, pages)
    m["names"]["english"] = f(name or label, name or label, pages)
    m["cas_numbers"] = st.fact(evidence=evidence)
    m["identity_note"] = f(note, note, pages) if note else st.fact(evidence=evidence)
    m["external_properties"] = []
    m["stereochemistry"]["chirality_status"] = inferred(
        chiral, pages, "Status assigned from the inspected source drawing/name; absence of a label is not treated as racemic.")
    m["stereochemistry"]["geometry_status"] = inferred(
        geometry, pages, "Status assigned from the inspected source drawing and candidate molecular graph.")
    m["stereochemistry"]["representation_status"] = inferred(
        representation, pages, "Reports whether the candidate representation preserves the inspected stereochemical information.")
    if assignment is not None:
        m["stereochemistry"]["assignment"] = f(assignment, assignment, pages)
    if mixture is not None:
        m["stereochemistry"]["mixture"] = f(mixture, mixture, pages)
    if axial is not None:
        m["stereochemistry"]["axial"] = f(axial, axial, pages)
    if smiles:
        diag = st.stereo_inspect(smiles)
        m["identity_status"] = "source_identified"
        m["representation_type"] = "smiles"
        m["smiles"] = derived(smiles, pages,
            "Machine candidate transcribed/derived from the inspected source structure or unambiguous formula; RDKit parsing is not source verification.",
            raw={"source_label": label, "candidate": smiles})
        m["stereochemistry"]["isomeric_smiles"] = derived(
            diag["canonical_isomeric_smiles"], pages,
            "RDKit canonical isomeric SMILES from the machine candidate; source_verified=false.", raw=smiles)
        m["stereochemistry"]["identity_key"] = derived(
            diag["stereo_graph_key"], pages,
            "RDKit stereo graph key from the machine candidate; not a sample identity.", raw=smiles)
    bundle["molecules"].append(m)
    return mid


# Concrete low-risk candidates used repeatedly. Each was passed through standard stereo above.
M = {}
M["1a"] = molecule("mol_1a", "1a", [3, 4, 5], name="methyl 2-(tosylamino)acrylate",
    smiles="C=C(NS(=O)(=O)c1ccc(C)cc1)C(=O)OC", chiral="not_applicable",
    geometry="not_applicable", representation="not_applicable",
    note="Candidate connectivity derived from the repeatedly drawn 1a structure; machine-only identity.")
M["2a"] = molecule("mol_2a", "2a", [3, 4, 5, 6], name="benzyl 2,3-butadienoate",
    smiles="C=C=CC(=O)OCc1ccccc1", chiral="not_applicable", geometry="not_applicable",
    representation="not_applicable", note="Terminal CH2 allene; no allene axial configuration applies to this candidate graph.")
M["2b"] = molecule("mol_2b", "2b", [4, 5, 6], name="ethyl 2,3-butadienoate",
    smiles="C=C=CC(=O)OCC", chiral="not_applicable", geometry="not_applicable", representation="not_applicable")
M["2c"] = molecule("mol_2c", "2c", [4, 5, 6], name="phenyl 2,3-butadienoate",
    smiles="C=C=CC(=O)Oc1ccccc1", chiral="not_applicable", geometry="not_applicable", representation="not_applicable")
M["2d"] = molecule("mol_2d", "2d", [4], name="methyl 2,3-butadienoate",
    smiles="C=C=CC(=O)OC", chiral="not_applicable", geometry="not_applicable", representation="not_applicable")
for key, name, smi, pages in [
    ("PPh3", "triphenylphosphine", "P(c1ccccc1)(c1ccccc1)c1ccccc1", [3, 4, 5, 6]),
    ("PBu3", "tributylphosphine", "CCCCP(CCCC)CCCC", [3]),
    ("PPhMe2", "dimethylphenylphosphine", "CP(C)c1ccccc1", [3]),
    ("PPh2Me", "methyldiphenylphosphine", "CP(c1ccccc1)c1ccccc1", [3]),
    ("PMe3", "trimethylphosphine", "CP(C)C", [3]),
    ("THF", "tetrahydrofuran", "C1CCOC1", [3, 4, 5]),
    ("dioxane", "1,4-dioxane", "C1COCCO1", [3]),
    ("PhMe", "toluene", "Cc1ccccc1", [3]),
    ("CH3CN", "acetonitrile", "CC#N", [3]),
    ("CDCl3", "chloroform-d", "[2H]C(Cl)(Cl)Cl", [3, 5, 6]),
    ("DCM", "dichloromethane", "ClCCl", [3]),
    ("Cs2CO3", "cesium carbonate", "[Cs+].[Cs+].[O-]C([O-])=O", [3]),
    ("H2", "hydrogen", "[H][H]", [5]),
    ("KOH", "potassium hydroxide", "[K+].[OH-]", [5]),
    ("H2O", "water", "O", [5]),
    ("HCl", "hydrochloric acid", "Cl", [5]),
]:
    M[key] = molecule("mol_" + key.replace("2", "_2").replace("3", "_3"), key, pages,
        name=name, smiles=smi, chiral="not_applicable", geometry="not_applicable", representation="not_applicable")

M["PdC"] = molecule("mol_palladium_on_carbon", "Pd/C", [5], name="palladium on carbon",
    note="Heterogeneous catalyst; no unique finite molecular structure.")
M["L"] = molecule("mol_chiral_catalyst_L", "L", [4], name="thiourea-derived chiral phosphine catalyst L",
    note="A specific chiral catalyst drawing is shown, but the absolute descriptor and a reliable complete machine structure are not established.",
    chiral="specified", geometry="unresolved", representation="not_encoded",
    assignment="A single wedge-defined catalyst stereoisomer is drawn; absolute configuration is not stated in the inspected page.")


def unresolved(label, pages, name=None, note=None, **stereo):
    key = label
    if key not in M:
        M[key] = molecule("mol_" + label.lower().replace("'", "prime").replace("/", "_").replace(" ", "_"),
                          label, pages, name=name, note=note, **stereo)
    return M[key]


# Product/sample identities. Product drawings contain a stereogenic quaternary carbon but
# standard entries do not assign absolute configuration; no racemate is invented.
for label in ["3aa", "3ba", "3ca", "3ea", "3fa", "3ga", "3ha", "3ia", "3ja", "3ka", "3la", "3ma", "3na", "3da", "3oa", "3pa",
              "3ab", "3bb", "3cb", "3db", "3ac", "3bc", "3cc", "3dc", "3ec", "3ad", "3an", "3ao", "3eb", "3fc"]:
    unresolved(label, [4] + ([3, 5] if label == "3aa" else []),
        name=f"3-pyrroline product {label}",
        note="Source-drawn product; quaternary stereogenic center has no absolute configuration assigned in the standard, achiral-catalyst entry.",
        chiral="unspecified", geometry="not_applicable", representation="not_encoded")
unresolved("3aa_prime", [4], name="enantioenriched 3aa'", chiral="mixture",
    geometry="not_applicable", representation="not_encoded",
    mixture="Enantiomeric mixture reported as 52% ee; major absolute configuration is unassigned.",
    note="Kept separate from standard 3aa because the sample is enantioenriched and the major enantiomer is not assigned.")
unresolved("4a", [5], name="hydrogenation product 4a", chiral="mixture",
    geometry="not_applicable", representation="partial",
    mixture="Diastereomeric mixture, dr >20:1; major relative configuration shown and assigned by single-crystal X-ray diffraction.",
    assignment="Major diastereomer relative configuration is shown with wedge/dash bonds and supported by Fig. 2 X-ray; absolute descriptors are not stated.")
unresolved("4b", [5], name="pyrrole carboxylic acid 4b", chiral="not_applicable",
    geometry="not_applicable", representation="unresolved")

# Scope reactants whose complete structures/substitution identity are retained by source description.
scope_substrates = {
    "sub_Ts_Me":"R1=Ts, R3=Me (1a)", "sub_4tBu_Me":"R1=SO2Ph(4-tBu), R3=Me",
    "sub_4OMe_Me":"R1=SO2Ph(4-OMe), R3=Me", "sub_Ph_Me":"R1=SO2Ph, R3=Me",
    "sub_246Me_Me":"R1=SO2Ph(2,4,6-Me), R3=Me", "sub_4Cl_Me":"R1=SO2Ph(4-Cl), R3=Me",
    "sub_Bs_Me":"R1=Bs, R3=Me", "sub_4F_Me":"R1=SO2Ph(4-F), R3=Me",
    "sub_3Br4NO2_Me":"R1=SO2Ph(3-Br,4-NO2), R3=Me", "sub_Ns_Me":"R1=Ns, R3=Me",
    "sub_4OCF3_Me":"R1=SO2Ph(4-OCF3), R3=Me", "sub_4CF3_Me":"R1=SO2Ph(4-CF3), R3=Me",
    "sub_4CN_Me":"R1=SO2Ph(4-CN), R3=Me", "sub_Ms_Me":"R1=Ms, R3=Me",
    "sub_thienylSO2_Me":"R1=thiophenesulfonyl, R3=Me", "sub_pyridylSO2_Me":"R1=methylpyridinesulfonyl, R3=Me",
    "sub_Boc_Me":"R1=Boc, R3=Me", "sub_Ac_Me":"R1=Ac, R3=Me",
    "sub_Ts_Et":"R1=Ts, R3=Et", "sub_Ts_Bn":"R1=Ts, R3=Bn",
    "1q":"1q (beta-substituted 2-aminoacrylate)", "1r":"1r (beta-substituted 2-aminoacrylate)",
}
for key, desc in scope_substrates.items():
    if key == "sub_Ts_Me":
        M[key] = M["1a"]
    else:
        M[key] = molecule("mol_" + key.lower(), desc, [4], name=desc,
            note="Complete candidate SMILES withheld because the page uses substituent notation or the exact substituent connectivity/identity is not fully secure.",
            chiral="unresolved", geometry="unresolved", representation="unresolved")
M["2e"] = molecule("mol_2e", "2e", [4], name="gamma-acetoxy-substituted allenoate 2e",
    note="No-reaction allenoate drawing; full connectivity retained as unresolved rather than guessed.",
    chiral="unresolved", geometry="unresolved", representation="unresolved")
M["2f"] = molecule("mol_2f", "2f", [4], name="gamma-phenyl-substituted allenoate 2f",
    note="The substituted allene may be axially stereogenic; no axial configuration is specified and no isomeric SMILES is invented.",
    chiral="unspecified", geometry="unresolved", representation="not_encoded",
    axial="Potential allene axis; configuration not specified in the source drawing.")


def participant(pid, mid, role, pages, *, stage="reaction", amount=None, equivalents=None,
                mol_percent=None, volume=None, mass=None, reference=None, product_role=None):
    p = st.new_participant(pid, mid, role, [ev(x) for x in pages], stage=stage)
    if product_role:
        p["product_role"] = f(product_role, product_role, pages)
    for key, val in [("amount", amount), ("equivalents", equivalents), ("mol_percent", mol_percent),
                     ("volume", volume), ("mass", mass)]:
        if val is not None:
            p["charges"][key] = f(val[0], val[1], pages)
    if reference:
        p["equivalence_reference"] = f(reference, reference, pages)
    return p


def change(change_id, a, b, before, after, kind):
    return {"change_id": change_id, "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
            "before_order": before, "after_order": after, "change_type": kind}


def classification(changes, pages, limitations, coverage="partial"):
    labels = {}
    for c in changes:
        if c["change_type"] in ("formed", "hydrogen_gain", "order_increased"):
            lab = st.bond_label(c["atoms"][0]["element"], c["atoms"][1]["element"], c["after_order"])
            labels.setdefault(lab["code"], dict(lab, supporting_change_ids=[]))["supporting_change_ids"].append(c["change_id"])
    why = "Partial source-drawing transformation analysis; atom identities are symbolic because a complete atom mapping was not asserted."
    return {"standard_id": st.STANDARD, "mapped_reaction_smiles": st.fact(evidence=[ev(p) for p in pages]),
            "bond_changes": [derived(c, pages, why) for c in changes],
            "labels": derived(list(labels.values()), pages, why),
            "coverage": inferred(coverage, pages, "Coverage limited by unresolved full structures, implicit-H/proton transfer, or omitted coproducts."),
            "limitations": limitations, "model_votes": []}


ANN_CHANGES = [change("b1", "C", "C", 0, 1, "formed"), change("b2", "C", "N", 0, 1, "formed")]


def reaction(rid, pages, record_type="specific_experiment", categories=None, location=None,
             name=None, transformation=None, changes=None, classification_coverage="partial",
             limitations=None, unresolved_items=None):
    evidence = [eid for p in pages for eid in (ev(p), ev(p, "image"))]
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, pages) if categories is not None else st.fact(evidence=evidence)
    r["location"] = f(location, location, pages) if location else f({"pages": pages}, {"pages": pages}, pages)
    r["reaction_name"] = f(name, name, pages) if name else st.fact(evidence=evidence)
    r["transformation_type"] = inferred(transformation, pages, "Normalized transformation description from the source title/scheme.") if transformation else st.fact(evidence=evidence)
    r["classification"] = classification(changes or [], pages, limitations or ["Full atom mapping not established."], classification_coverage)
    r["review"]["unresolved"] = unresolved_items or ["Machine-authored extraction requires expert source review."]
    r["review"]["status"] = "needs_review"
    r["review"]["method"] = "model"
    r["completeness"] = "partial"
    bundle["reactions"].append(r)
    return r


def step(r, sid, pages, pids, operation, temp=None, time=None, atmosphere=None, stage="reaction"):
    s = st.new_step(sid, [ev(p) for p in pages], stage=stage)
    s["participant_ids"] = pids
    s["operation"] = f(operation, operation, pages)
    if temp is not None: s["temperature"] = f(temp[0], temp[1], pages)
    if time is not None: s["time"] = f(time[0], time[1], pages)
    if atmosphere is not None: s["atmosphere"] = f(atmosphere, atmosphere, pages)
    r["steps"].append(s)
    return s


def solvent(r, raw, normalized, pages, volume=None):
    r["solvent_conditions"]["composition"] = f(raw, normalized, pages)
    if volume is not None:
        r["solvent_conditions"]["total_volume"] = f(volume[0], volume[1], pages)


def outcome(r, oid, kind, targets, pages, raw, normalized, *, basis=None, yield_type=None, scope=None,
            method="original_extraction", derivation=None):
    o = st.new_outcome(oid, kind, targets, [ev(p) for p in pages])
    o["value"] = f(raw, normalized, pages, method=method, derivation=derivation)
    if basis is not None: o["measurement_basis"] = f(basis, basis, pages)
    if yield_type is not None: o["yield_type"] = f(yield_type, yield_type, pages)
    if scope is not None: o["scope"] = f(scope, scope, pages)
    r["outcomes"].append(o)
    return o


def set_procedure(r, raw, pages, *, general=None, addition=None):
    r["procedure"]["original_text"] = f(raw, raw, pages)
    if general: r["procedure"]["general_procedure_reference"] = f(general, general, pages)
    if addition: r["procedure"]["addition_order"] = f(addition, addition, pages)


def annulation_participants(r, pages, sub_mid, allene_mid, product_mid=None, *, a_amount=None, b_amount=None,
                            catalyst_mid=None, catalyst_loading=None, solvent_mid=None, solvent_volume=None,
                            base_mid=None, base_amount=None):
    ps = []
    ps.append(participant("substrate", sub_mid, "reactant", pages, amount=a_amount))
    ps.append(participant("allenoate", allene_mid, "reactant", pages, amount=b_amount))
    if catalyst_mid:
        ps.append(participant("catalyst", catalyst_mid, "catalyst", pages, mol_percent=catalyst_loading, reference=sub_mid))
    if base_mid:
        ps.append(participant("base", base_mid, "additive", pages, amount=base_amount))
    if solvent_mid:
        ps.append(participant("solvent", solvent_mid, "solvent", pages, volume=solvent_volume))
    if product_mid:
        ps.append(participant("product", product_mid, "product", pages, product_role="main product"))
    r["participants"] = ps
    return [p["participant_id"] for p in ps]


# Scheme 1: two cited previous-work generic reactions and the authors' generic scheme.
for tag, desc in [("prev_alkene", "Previous work a: allenoate plus activated alkene"),
                  ("prev_ketimine", "Previous work b: allenoate plus ketimine")]:
    a = unresolved(tag + "_allenoate", [3], name="generic allenoate", note="Generic R-group structure in cited previous work.")
    b = unresolved(tag + "_partner", [3], name=desc, note="Generic R-group structure in cited previous work.")
    p = unresolved(tag + "_product", [3], name=desc + " cycloadduct", note="Generic cited-product drawing.")
    r = reaction("scheme1_" + tag, [3], record_type="literature_example", location="Scheme 1, Previous work",
                 name="phosphine-catalyzed [3+2] cycloaddition", transformation="[3+2] annulation",
                 changes=ANN_CHANGES, unresolved_items=["Cited prior-work generic scheme; not an experiment performed in this article."])
    r["participants"] = [participant("allenoate", a, "reactant", [3]), participant("partner", b, "reactant", [3]),
                         participant("product", p, "product", [3], product_role="generic cycloadduct")]
    r["notes"].append(f("Previous work", "Previous work", [3]))

generic_sub = unresolved("generic_2_aminoacrylate", [3], name="generic 2-aminoacrylate 1",
                         note="R1/R3 generic substrate in Scheme 1.")
generic_allene = unresolved("generic_allenoate_2", [3], name="generic allenoate 2",
                            note="R2 generic allenoate in Scheme 1.")
generic_product = unresolved("generic_3_pyrroline_3", [3], name="generic 3-pyrroline 3",
                             note="R1/R2/R3 generic product in Scheme 1.", chiral="unspecified",
                             geometry="not_applicable", representation="not_encoded")
r = reaction("scheme1_this_work_generic", [3], record_type="generic_scheme", location="Scheme 1, This work",
             name="Phosphine-Catalyzed [3 + 2] Cycloaddition Reactions", transformation="[3+2] annulation",
             changes=ANN_CHANGES, unresolved_items=["Generic R groups are not expanded to concrete molecules.",
                                                    "Scheme states N2/r.t.; model prose later states argon for the initial run."])
r["participants"] = [participant("substrate", generic_sub, "reactant", [3]), participant("allenoate", generic_allene, "reactant", [3]),
                     participant("catalyst", M["PPh3"], "catalyst", [3]), participant("product", generic_product, "product", [3], product_role="generic main product")]
step(r, "generic_reaction", [3], ["substrate", "allenoate", "catalyst", "product"],
     "phosphine-catalyzed annulation", ("r.t.", "r.t."), atmosphere="N2")
r["notes"].append(f("In situ generated ketoimine", "In situ generated ketoimine", [3]))


# Table 1 optimization. Entry 4 is the same standard 3aa experiment repeated in Table 2.
opt_rows = [
    (1,"PPh3",None,"dioxane","trace"),(2,"PPh3",None,"PhMe",47),(3,"PPh3",None,"CH3CN",27),
    (4,"PPh3",None,"THF",92),(5,"PPh3",None,"CDCl3",43),(6,"PPh3",None,"DCM",32),
    (7,"PBu3",None,"THF","trace"),(8,"PPhMe2",None,"THF",80),(9,"PPh2Me",None,"THF",88),
    (10,"PMe3",None,"THF","trace"),(11,"PPh3",None,"THF",82),(12,"PPh3",None,"THF",76),
    (13,"PPh3","Cs2CO3","THF",56),(14,None,None,"THF",0),
]
standard_r = None
for entry, cat, base, solv, y in opt_rows:
    rid = "rxn_model_3aa_opt_scope" if entry == 4 else f"table1_entry_{entry:02d}"
    pages = [3,4] if entry == 4 else [3]
    categories = ["optimization", "substrate_scope"] if entry == 4 else ["optimization"]
    changes = [] if entry == 14 else ANN_CHANGES
    coverage = "not_applicable" if entry == 14 else "partial"
    r = reaction(rid, pages, categories=categories, location=f"Table 1 entry {entry}" + ("; Table 2 3aa" if entry == 4 else ""),
                 name="phosphine-catalyzed [3+2] annulation of 1a with 2a", transformation="[3+2] annulation",
                 changes=changes, classification_coverage=coverage,
                 limitations=["No complete atom mapping; proton transfers and omitted species are unresolved."] if entry != 14 else ["No reaction observed; intended product bonds are not classified."],
                 unresolved_items=["Yield measurement basis is unreported." if entry != 4 else "Only the 3aa Table 2 footnote explicitly identifies isolated yield.",
                                   "Machine candidate structures and bond classification require expert review."])
    a_amt = ("0.20 mmol", {"number":0.20,"unit":"mmol"}) if entry == 12 else ("0.10 mmol", {"number":0.10,"unit":"mmol"})
    b_amt = ("0.20 mmol", {"number":0.20,"unit":"mmol"}) if entry == 11 else ("0.10 mmol", {"number":0.10,"unit":"mmol"}) if entry == 12 else ("0.12 mmol", {"number":0.12,"unit":"mmol"})
    cat_load = ("10 mol%", {"number":10,"unit":"mol%"}) if cat else None
    product_mid = None if entry == 14 else M["3aa"]
    pids = annulation_participants(r, pages, M["1a"], M["2a"], product_mid,
        a_amount=a_amt, b_amount=b_amt, catalyst_mid=M[cat] if cat else None, catalyst_loading=cat_load,
        solvent_mid=M[solv], solvent_volume=("2.0 mL", {"number":2.0,"unit":"mL"}),
        base_mid=M[base] if base else None,
        base_amount=("0.2 mmol", {"number":0.2,"unit":"mmol"}) if base else None)
    r["reaction_scale"] = f("1a (0.10 mmol)" if entry != 12 else "1a (0.20 mmol)",
                            {"reference_molecule_id":M["1a"],"amount":a_amt[1]}, pages)
    solvent(r, solv, [M[solv]], pages, ("2.0 mL", {"number":2.0,"unit":"mL"}))
    atmosphere = "argon" if entry == 1 else None
    step(r, "reaction", pages, pids, "carried out in a reaction tube", ("room temperature", "r.t."),
         ("4 h", {"number":4,"unit":"h"}), atmosphere)
    set_procedure(r, "Unless otherwise indicated: 1a, 2a and catalyst in solvent in a reaction tube for 4 h at room temperature.", pages,
                  addition="1a was treated with 2a and phosphine catalyst")
    if y == "trace":
        outcome(r, "trace_product", "trace", [M["3aa"]], pages, "trace", "trace")
    elif entry == 14:
        outcome(r, "no_reaction", "no_reaction", [], pages, "yield 0%; no PPh3", True,
                scope="absence of PPh3")
    else:
        basis = "isolated yield" if entry == 4 else None
        outcome(r, "product_yield", "yield", [M["3aa"]], pages, f"{y}%", {"number":y,"unit":"%"}, basis=basis,
                yield_type="isolated" if entry == 4 else None)
    if entry == 1:
        r["notes"].append(f("Prose states argon for the initial 1,4-dioxane run; Scheme 1 generic conditions state N2.",
                            "Atmosphere statements have distinct source scopes and were not silently reconciled.", [3]))
        r["review"]["unresolved"].append("Atmosphere differs between initial-run prose (argon) and generic Scheme 1 (N2); Table 1 footnote is silent.")
    if entry == 4: standard_r = r


# Table 2 positive scope entries. 3aa is already represented by the merged standard record.
scope_rows = [
    ("sub_4tBu_Me","2a","3ba",90),("sub_4OMe_Me","2a","3ca",92),("sub_Ph_Me","2a","3ea",70),
    ("sub_246Me_Me","2a","3fa",83),("sub_4Cl_Me","2a","3ga",95),("sub_Bs_Me","2a","3ha",70),
    ("sub_4F_Me","2a","3ia",60),("sub_3Br4NO2_Me","2a","3ja",64),("sub_Ns_Me","2a","3ka",60),
    ("sub_4OCF3_Me","2a","3la",42),("sub_4CF3_Me","2a","3ma",56),("sub_4CN_Me","2a","3na",55),
    ("sub_Ms_Me","2a","3da",75),("sub_thienylSO2_Me","2a","3oa",50),("sub_pyridylSO2_Me","2a","3pa",80),
    ("sub_Ts_Me","2b","3ab",92),("sub_4tBu_Me","2b","3bb",74),("sub_4OMe_Me","2b","3cb",85),("sub_4CF3_Me","2b","3db",83),
    ("sub_Ts_Me","2c","3ac",79),("sub_4tBu_Me","2c","3bc",74),("sub_4OMe_Me","2c","3cc",71),("sub_4Cl_Me","2c","3dc",84),("sub_4CF3_Me","2c","3ec",87),
    ("sub_Ts_Me","2d","3ad",88),("sub_Ts_Et","2a","3an",86),("sub_Ts_Bn","2a","3ao",94),
    ("sub_Ts_Bn","2b","3eb",80),("sub_Ts_Bn","2c","3fc",80),
]
for sub, allene, prod, y in scope_rows:
    r = reaction(f"table2_{prod}", [4], categories=["substrate_scope"], location=f"Table 2, {prod}",
        name="phosphine-catalyzed [3+2] annulation", transformation="[3+2] annulation", changes=ANN_CHANGES,
        limitations=["Partial bond classification from generic source drawings; no complete atom mapping."],
        unresolved_items=["Yield measurement basis is not explicitly assigned to this row.",
                          "Full source-drawn substrate/product structure requires expert confirmation."])
    pids = annulation_participants(r, [4], M[sub], M[allene], M[prod],
        a_amount=("0.1 mmol", {"number":0.1,"unit":"mmol"}), b_amount=("0.12 mmol", {"number":0.12,"unit":"mmol"}),
        catalyst_mid=M["PPh3"], catalyst_loading=("10 mol%", {"number":10,"unit":"mol%"}),
        solvent_mid=M["THF"], solvent_volume=("2.0 mL", {"number":2.0,"unit":"mL"}))
    r["reaction_scale"] = f("1 (0.1 mmol)", {"reference_molecule_id":M[sub],"amount":{"number":0.1,"unit":"mmol"}}, [4])
    solvent(r, "THF", [M["THF"]], [4], ("2.0 mL", {"number":2.0,"unit":"mL"}))
    step(r, "reaction", [4], pids, "standard Table 2 reaction in a reaction tube", ("room temperature", "r.t."),
         ("4 h", {"number":4,"unit":"h"}))
    set_procedure(r, "Unless otherwise indicated, all reactions were carried out with 1 (0.1 mmol), 2 (0.12 mmol), PPh3 (10 mol%) in THF (2.0 mL) at room temperature in a reaction tube for 4 h.", [4])
    outcome(r, "product_yield", "yield", [M[prod]], [4], f"{y}%", {"number":y,"unit":"%"})


# Table 2 no-reaction controls/scope failures.
no_rows = [
    ("Boc", "sub_Boc_Me", "2a"), ("Ac", "sub_Ac_Me", "2a"), ("1q", "1q", "2a"),
    ("1r", "1r", "2a"), ("2e", "sub_Ts_Me", "2e"), ("2f", "sub_Ts_Me", "2f"),
]
for label, sub, allene in no_rows:
    r = reaction(f"table2_no_reaction_{label.lower()}", [4], categories=["substrate_scope", "control"],
        location=f"Table 2, {label}, no reaction", name="attempted phosphine-catalyzed [3+2] annulation",
        transformation="no observed annulation", changes=[], classification_coverage="not_applicable",
        limitations=["No reaction observed; intended product bonds are not classified."],
        unresolved_items=["Observation/detection limit and measurement basis are not reported.",
                          "No product molecule was invented for a failed experiment."])
    pids = annulation_participants(r, [4], M[sub], M[allene], None,
        a_amount=("0.1 mmol", {"number":0.1,"unit":"mmol"}), b_amount=("0.12 mmol", {"number":0.12,"unit":"mmol"}),
        catalyst_mid=M["PPh3"], catalyst_loading=("10 mol%", {"number":10,"unit":"mol%"}),
        solvent_mid=M["THF"], solvent_volume=("2.0 mL", {"number":2.0,"unit":"mL"}))
    r["reaction_scale"] = f("1 (0.1 mmol)", {"reference_molecule_id":M[sub],"amount":{"number":0.1,"unit":"mmol"}}, [4])
    solvent(r, "THF", [M["THF"]], [4], ("2.0 mL", {"number":2.0,"unit":"mL"}))
    step(r, "reaction", [4], pids, "standard Table 2 attempted reaction", ("room temperature", "r.t."),
         ("4 h", {"number":4,"unit":"h"}))
    set_procedure(r, "Standard Table 2 conditions; source reports no reaction.", [4])
    outcome(r, "no_reaction", "no_reaction", [], [4], "no reaction", True, scope=label)


# Scheme 2 asymmetric example.
r = reaction("scheme2_asymmetric_3aa_prime", [4], categories=["substrate_scope"], location="Scheme 2",
    name="Asymmetric reaction of 2-aminoacrylate 1a with allenoate 2a", transformation="enantioselective [3+2] annulation",
    changes=ANN_CHANGES, limitations=["No complete atom mapping; major absolute configuration is unassigned."],
    unresolved_items=["Catalyst loading, allenoate charge and reaction time are not reported on the inspected page.",
                      "52% ee does not establish the major absolute configuration."])
r["participants"] = [participant("substrate", M["1a"], "reactant", [4], amount=("0.1 mmol scale", {"number":0.1,"unit":"mmol"})),
                     participant("allenoate", M["2a"], "reactant", [4]), participant("catalyst", M["L"], "catalyst", [4]),
                     participant("solvent", M["THF"], "solvent", [4]),
                     participant("product", M["3aa_prime"], "product", [4], product_role="enantioenriched main product")]
r["reaction_scale"] = f("Reaction was run on 0.1 mmol scale", {"reference_molecule_id":M["1a"],"amount":{"number":0.1,"unit":"mmol"}}, [4])
solvent(r, "THF", [M["THF"]], [4])
step(r, "reaction", [4], [p["participant_id"] for p in r["participants"]], "asymmetric annulation", ("r.t.", "r.t."))
set_procedure(r, "Reaction was run on 0.1 mmol scale in THF at r.t.", [4])
outcome(r, "isolated_yield", "yield", [M["3aa_prime"]], [4], "yield: 92%", {"number":92,"unit":"%"}, basis="isolated yield", yield_type="isolated")
outcome(r, "enantioselectivity", "selectivity", [M["3aa_prime"]], [4], "ee: 52%", {"number":52,"unit":"%"},
        basis="HPLC analysis on a chiral stationary phase", scope="enantiomeric excess; major enantiomer absolute configuration unassigned")


# Scheme 3 derivatizations and gram-scale repetition.
r = reaction("scheme3_hydrogenation_4a", [5], categories=["derivatization"], location="Scheme 3a, hydrogenation",
    name="Hydrogenation of 3-pyrroline 3aa", transformation="hydrogenation and benzyl ester hydrogenolysis",
    changes=[change("b1","C","C",2,1,"order_decreased"), change("b2","C","H",0,1,"hydrogen_gain"),
             change("b3","C","H",0,1,"hydrogen_gain"), change("b4","O","H",0,1,"hydrogen_gain")],
    limitations=["Hydrogen changes are net gains; proton/hydrogen source mapping and benzyl coproducts are incomplete."],
    unresolved_items=["Reaction time, temperature, pressure, solvent and catalyst loading are not reported on the inspected page.",
                      "Product reported as dr >20:1; absolute descriptors are not stated."])
r["participants"] = [participant("substrate", M["3aa"], "reactant", [5]), participant("hydrogen", M["H2"], "reagent", [5]),
                     participant("catalyst", M["PdC"], "catalyst", [5]), participant("product", M["4a"], "product", [5], product_role="main product")]
step(r, "reaction", [5], ["substrate","hydrogen","catalyst","product"], "H2, Pd/C hydrogenation")
set_procedure(r, "3aa --H2, Pd/C--> 4a", [5])
outcome(r, "product_yield", "yield", [M["4a"]], [5], "yield: 100%", {"number":100,"unit":"%"})
outcome(r, "diastereoselectivity", "selectivity", [M["4a"]], [5], "dr > 20:1", {"qualifier":"gt","ratio":[20,1],"raw_text":"dr > 20:1"},
        scope="diastereomeric ratio; major/minor absolute identities not named")

r = reaction("scheme3_koh_hydrolysis_4b", [5], categories=["derivatization"], location="Scheme 3a, KOH/HCl sequence",
    name="Base-promoted conversion of 3aa to pyrrole carboxylic acid 4b", transformation="ester cleavage/aromatization",
    changes=[change("b1","C","C",1,2,"order_increased"), change("b2","O","H",0,1,"hydrogen_gain")],
    limitations=["Partial net transformation only; benzyl fragment, proton transfers and full aromatic bond redistribution are unresolved."],
    unresolved_items=["KOH/HCl quantities, time, temperature and purification are not reported on the inspected page.",
                      "The prose calls the event debenzylation; the drawing also indicates pyrrole unsaturation, so classification is partial."])
r["participants"] = [participant("substrate", M["3aa"], "reactant", [5]), participant("koh", M["KOH"], "reagent", [5]),
                     participant("thf", M["THF"], "solvent", [5]), participant("water", M["H2O"], "solvent", [5]),
                     participant("hcl", M["HCl"], "workup_agent", [5], stage="workup"),
                     participant("product", M["4b"], "product", [5], product_role="main product")]
step(r, "base_treatment", [5], ["substrate","koh","thf","water"], "1. KOH, THF, H2O")
step(r, "acidification", [5], ["hcl","product"], "2. HCl aq.", stage="workup")
solvent(r, "THF, H2O", [M["THF"],M["H2O"]], [5])
set_procedure(r, "1. KOH, THF, H2O; 2. HCl aq.", [5])
outcome(r, "product_yield", "yield", [M["4b"]], [5], "yield: 56%", {"number":56,"unit":"%"})

r = reaction("scheme3_scaleup_3aa", [5], categories=["scale_up"], location="Scheme 3b",
    name="Gram-scale synthesis of 3aa", transformation="[3+2] annulation", changes=ANN_CHANGES,
    limitations=["Standard conditions are inherited qualitatively; scaled solvent volume is not calculated or invented."],
    unresolved_items=["The Scheme 3 arrow says Standard Reaction Conditions; exact scaled catalyst and solvent charges are not stated here.",
                      "Yield measurement basis is not explicitly stated for the scale-up."])
r["participants"] = [participant("substrate", M["1a"], "reactant", [5], amount=("12 mmol", {"number":12,"unit":"mmol"})),
                     participant("allenoate", M["2a"], "reactant", [5], amount=("14.5 mmol", {"number":14.5,"unit":"mmol"})),
                     participant("catalyst", M["PPh3"], "catalyst", [5], mol_percent=("standard conditions: 10 mol%", {"number":10,"unit":"mol%"}), reference=M["1a"]),
                     participant("solvent", M["THF"], "solvent", [5]),
                     participant("product", M["3aa"], "product", [5], mass=("3.53 g", {"number":3.53,"unit":"g"}), product_role="main product")]
r["participants"][2]["charges"]["mol_percent"] = f(
    {"scheme3":"Standard Reaction Conditions","table2":"PPh3 (10 mol%)"},
    {"number":10,"unit":"mol%"}, [4,5], method="inherited_general_condition")
r["reaction_scale"] = f("1a (12 mmol)", {"reference_molecule_id":M["1a"],"amount":{"number":12,"unit":"mmol"}}, [5])
r["solvent_conditions"]["composition"] = f(
    {"scheme3":"Standard Reaction Conditions","table2":"THF"}, [M["THF"]], [4,5], method="inherited_general_condition")
scale_step = step(r, "reaction", [5], [p["participant_id"] for p in r["participants"]], "Standard Reaction Conditions")
scale_step["temperature"] = f(
    {"scheme3":"Standard Reaction Conditions","table2":"room temperature"}, "r.t.", [4,5], method="inherited_general_condition")
scale_step["time"] = f(
    {"scheme3":"Standard Reaction Conditions","table2":"4 h"}, {"number":4,"unit":"h"}, [4,5], method="inherited_general_condition")
set_procedure(r, "Standard Reaction Conditions; 1a (12 mmol), 2a (14.5 mmol)", [5], general="Table 2 standard conditions")
outcome(r, "product_amount", "product_amount", [M["3aa"]], [5], "3.53 g", {"number":3.53,"unit":"g"})
outcome(r, "product_yield", "yield", [M["3aa"]], [5], "yield: 83%", {"number":83,"unit":"%"})


# Scheme 4 is an author-proposed catalytic cycle, not a performed experiment.
r = reaction("scheme4_proposed_mechanism", [5], record_type="proposed_mechanism", location="Scheme 4",
    name="Proposed Reaction Mechanism", transformation="proposed phosphine-catalyzed [3+2] cycle",
    changes=ANN_CHANGES, limitations=["Arrow-pushing proposal; intermediates I-IV are not treated as experimentally isolated products."],
    unresolved_items=["Proposed intermediates and proton shift are mechanistic claims, not independently verified synthesis records."])
r["participants"] = [participant("substrate", M["1a"], "reactant", [5]), participant("allenoate", M["2a"], "reactant", [5]),
                     participant("catalyst", M["PPh3"], "catalyst", [5]), participant("product", M["3aa"], "product", [5], product_role="proposed cycle product")]
step(r, "proposed_cycle", [5], ["substrate","allenoate","catalyst","product"], "Proposed zwitterion I, ylide II, cyclization III, [1,2]-proton shift IV and elimination",
     ("r.t.", "r.t."))
solvent(r, "THF", [M["THF"]], [5])
r["notes"].append(f("An active intermediate could be detected by 31P NMR during the reaction (details said to be in Fig. S7 in SI).",
                    "Author statement; SI was not supplied or associated in this session.", [5]))


# Figures 3/4: six concrete time-dependent 1H-NMR experiments. Values are approximate graph readings.
curve10 = {
    "2a":[[0,0],[2,1],[4,3],[6,4],[8,5],[10,8],[12,9],[14,11],[16,13]],
    "2b":[[0,0],[2,9],[4,12],[6,14],[8,16],[10,19],[12,20],[14,23],[16,25]],
    "2c":[[0,0],[2,11],[4,16],[6,20],[8,23],[10,26],[12,30],[14,31],[16,33]],
}
curve30 = {
    "2a":[[0,0],[2,1],[4,2],[6,3],[8,4],[10,5],[12,5.5],[14,6],[16,7]],
    "2b":[[0,0],[2,21],[4,25],[6,29],[8,30],[10,32],[12,33.5],[14,34.5],[16,35.5]],
    "2c":[[0,0],[2,30],[4,33],[6,37],[8,40],[10,42],[12,43],[14,44],[16,45]],
}
product_for = {"2a":"3aa","2b":"3ab","2c":"3ac"}
for loading, curves, figure, pages in [(10,curve10,"Fig. 3",[5,6]), (30,curve30,"Fig. 4",[5,6])]:
    for allene, points in curves.items():
        prod = product_for[allene]
        r = reaction(f"kinetics_{allene}_{loading}molpct_pph3", pages, categories=["mechanistic_study"],
            location=f"{figure}, {allene}, {loading} mol% PPh3", name="Time-dependent 1H-NMR study",
            transformation="monitored [3+2] annulation", changes=ANN_CHANGES,
            limitations=["Curve values are approximate visual digitizations, not author-tabulated values; no full atom mapping."],
            unresolved_items=["1a/allenoate charges and temperature are not reported in the inspected main-text pages.",
                              "Approximate curve readings require independent digitization if numerical precision is needed."])
        r["participants"] = [participant("substrate", M["1a"], "reactant", pages), participant("allenoate", M[allene], "reactant", pages),
                             participant("catalyst", M["PPh3"], "catalyst", pages, mol_percent=(f"{loading} mol%", {"number":loading,"unit":"mol%"}), reference=M["1a"]),
                             participant("solvent", M["CDCl3"], "solvent", pages),
                             participant("product", M[prod], "product", pages, product_role="NMR-monitored product")]
        solvent(r, "CDCl3", [M["CDCl3"]], pages)
        step(r, "nmr_monitoring", pages, [p["participant_id"] for p in r["participants"]],
             "time-dependent 1H-NMR monitoring", time=("0-16 min", {"minimum":0,"maximum":16,"unit":"min","raw_text":"0-16 min"}))
        set_procedure(r, f"Reaction of 1a with {allene} in the presence of {loading} mol% PPh3 in CDCl3; time-dependent 1H-NMR study.", pages)
        outcome(r, "yield_time-curve", "other", [M[prod]], pages, f"{figure} plotted yield/time curve",
                [{"time_min":t,"yield_percent_approx":y} for t,y in points], basis="time-dependent 1H NMR",
                scope=f"Approximate graph readings for {allene}; 0-16 min",
                method="structure_derivation", derivation="Approximate visual reading from the plotted markers; values are not source-tabulated.")


# Non-reaction observations and analytical context remain source_context.
bundle["source_context"] = [
    {"context_id":f"{DOC}:page1_cover", "document_id":DOC, "evidence_ids":[ev(1)],
     "kind":"publication_front_matter", "summary":"Accepted-manuscript cover; no reaction experiment."},
    {"context_id":f"{DOC}:page2_title", "document_id":DOC, "evidence_ids":[ev(2)],
     "kind":"bibliographic", "title":"Phosphine-catalyzed [3 + 2] annulation of 2-aminoacrylates with allenoates and mechanistic studies",
     "authors":["Xiao-Yun Wu","Hou-Ze Gui","Harish Jangra","Yin Wei","Hendrik Zipse","Min Shi"], "doi":"10.1039/D0CY00092B"},
    {"context_id":f"{DOC}:page3_figure1", "document_id":DOC, "evidence_ids":[ev(3)],
     "kind":"background_figure", "summary":"Fig. 1 shows representative natural products/therapeutic 3-pyrrolines; not experiments performed here."},
    {"context_id":f"{DOC}:page5_xray", "document_id":DOC, "evidence_ids":[ev(5)],
     "kind":"characterization", "summary":"Fig. 2 single-crystal X-ray structure of 4a; supports structural assignment, not a separate synthesis."},
    {"context_id":f"{DOC}:page5_31p_nmr", "document_id":DOC, "evidence_ids":[ev(5)],
     "kind":"mechanistic_observation", "summary":"Text states an active intermediate was detected by 31P NMR; details refer to unavailable/unassociated SI Fig. S7."},
    {"context_id":f"{DOC}:page6_maa", "document_id":DOC, "evidence_ids":[ev(6)],
     "kind":"calculation", "summary":"Fig. 5 calculated Methyl Anion Affinity values (kJ/mol): 2a 163.1, 2b 163.4, 2c 182.5; not synthesis yields."},
    {"context_id":f"{DOC}:page6_general_methods", "document_id":DOC, "evidence_ids":[ev(6)],
     "kind":"analytical_methods", "summary":"General NMR/TLC/column/IR/MS instrumentation and solvent-drying statements; no product-specific procedure on page 6."},
]


# Raw transcription is deliberately separate from the normalized submission.
raw = {
    "document_id": DOC,
    "pages_inspected": [1,2,3,4,5,6],
    "semantic_accuracy_verified": False,
    "table1_rows": [{"entry":e,"catalyst":c,"base":b,"solvent":s,"yield":y} for e,c,b,s,y in opt_rows],
    "table2_positive_rows": [{"substrate_key":s,"allenoate":a,"product":p,"yield_percent":y} for s,a,p,y in scope_rows],
    "table2_no_reaction_rows": [{"label":l,"substrate_key":s,"allenoate":a} for l,s,a in no_rows],
    "figure3_10molpct_approx": curve10,
    "figure4_30molpct_approx": curve30,
    "notes": [
        "Table 1 entry 4 and Table 2 3aa were treated as the same standard experiment and merged with both page citations.",
        "Only 3aa carries the Table 2 b footnote 'Isolated yield'; other scope yield bases remain unreported.",
        "Curve values are approximate image readings and are encoded as derived, pending-review data.",
        "No SI was associated; SI-only details were not imported.",
    ],
}


# Coverage maps every inspected page and locatable region to records/decisions.
record_ids_by_page = {p:[] for p in range(1,7)}
for r in bundle["reactions"]:
    for eid in r["evidence_ids"]:
        p = int(eid.split(":page",1)[1].split(":",1)[0])
        if p in record_ids_by_page and r["reaction_id"] not in record_ids_by_page[p]:
            record_ids_by_page[p].append(r["reaction_id"])
coverage = {"document_id":DOC,"semantic_accuracy_verified":False,"pages":[
    {"page":1,"state":"attempted","regions":[
        {"region":"full-page accepted-manuscript cover","finding":"Publication front matter only; no chemical experiment.","evidence_ids":[ev(1),ev(1,"image")]}
    ],"record_ids":record_ids_by_page[1],"unresolved":[]},
    {"page":2,"state":"attempted","regions":[
        {"region":"title, authors and abstract","finding":"Bibliographic metadata and article-level summary; no concrete reaction row.","evidence_ids":[ev(2),ev(2,"image")]},
        {"region":"introduction and ESI footnote","finding":"Background/cited work; footnote says ESI available, but no SI file was supplied or associated.","evidence_ids":[ev(2),ev(2,"image")]}
    ],"record_ids":record_ids_by_page[2],"unresolved":["Supporting information is mentioned but not provided/associated."]},
    {"page":3,"state":"attempted","regions":[
        {"region":"Fig. 1 upper left","finding":"Background example molecules and bioactivity; not experiments performed in this article.","evidence_ids":[ev(3),ev(3,"image")]},
        {"region":"Scheme 1 lower left","finding":"Two cited previous-work generic schemes and one this-work generic scheme; classified separately from concrete experiments.","evidence_ids":[ev(3),ev(3,"image")]},
        {"region":"optimization prose and Table 1 right","finding":"Fourteen optimization rows transcribed; entry 4 is merged with repeated Table 2 3aa standard result.","evidence_ids":[ev(3),ev(3,"image")]}
    ],"record_ids":record_ids_by_page[3],"unresolved":["Generic Scheme 1 says N2 while initial-run prose says argon; source scopes retained separately.","Most optimization yield measurement bases are unreported."]},
    {"page":4,"state":"attempted","regions":[
        {"region":"Table 2 upper right and continuation prose","finding":"Thirty positive scope products plus six no-reaction substrate/allenoate cases inventoried; 3aa is duplicate evidence for the standard experiment.","evidence_ids":[ev(4),ev(4,"image")]},
        {"region":"Scheme 2 lower right","finding":"One asymmetric 1a/2a experiment: 92% isolated yield, 52% ee by chiral-stationary-phase HPLC.","evidence_ids":[ev(4),ev(4,"image")]},
        {"region":"synthetic-utility prose lower right","finding":"Points forward to Scheme 3 on page 5; not duplicated as separate experiments here.","evidence_ids":[ev(4),ev(4,"image")]}
    ],"record_ids":record_ids_by_page[4],"unresolved":["Many scope structures use substituent notation and remain without asserted concrete SMILES.","Only 3aa is explicitly footnoted as isolated yield."]},
    {"page":5,"state":"attempted","regions":[
        {"region":"Scheme 3 and Fig. 2 left upper","finding":"Hydrogenation to 4a, KOH/HCl conversion to 4b, and 12 mmol scale-up recorded; X-ray is characterization evidence, not a separate synthesis.","evidence_ids":[ev(5),ev(5,"image")]},
        {"region":"Scheme 4 left lower","finding":"Author-proposed mechanism recorded as proposed_mechanism; intermediates are not treated as isolated experiments.","evidence_ids":[ev(5),ev(5,"image")]},
        {"region":"kinetics prose and lower-right plot","finding":"Three 10 mol% PPh3 time-dependent NMR experiments (2a/2b/2c); curve readings retained as approximate derived values.","evidence_ids":[ev(5),ev(5,"image")]}
    ],"record_ids":record_ids_by_page[5],"unresolved":["Derivatization time/temperature/charges are absent from the inspected page.","31P-NMR intermediate details are referred to unavailable SI."]},
    {"page":6,"state":"attempted","regions":[
        {"region":"Fig. 3/4 plots left upper","finding":"Captions complete the 10 mol% series and show three 30 mol% PPh3 NMR kinetic experiments; graph readings are approximate.","evidence_ids":[ev(6),ev(6,"image")]},
        {"region":"Fig. 5 left lower","finding":"Calculated MAA values retained as source_context, not reaction experiments.","evidence_ids":[ev(6),ev(6,"image")]},
        {"region":"conclusion/general experimental method right","finding":"Article-level conclusion and analytical instrumentation; no new product-specific synthesis.","evidence_ids":[ev(6),ev(6,"image")]}
    ],"record_ids":record_ids_by_page[6],"unresolved":["Exact numeric kinetic values are not tabulated and require dedicated digitization for precision."]},
]}

(RUN / "raw_model/pages_0001_0006_extraction.json").write_text(json.dumps(raw, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/stereo_diagnostics.json").write_text(json.dumps({
    m["molecule_id"]: st.stereo_inspect(st.value(m["smiles"]))
    for m in bundle["molecules"] if m["representation_type"] == "smiles"
}, ensure_ascii=False, indent=2) + "\n")
(RUN / "coverage/pages_0001_0006.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")
(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"molecules":len(bundle["molecules"]),"reactions":len(bundle["reactions"]),
                  "source_context":len(bundle["source_context"]),"records_by_page":{k:len(v) for k,v in record_ids_by_page.items()}},
                 ensure_ascii=False, indent=2))
