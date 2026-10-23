#!/usr/bin/env python3
"""Incremental v3.1 extraction for visually inspected PDF pages 7--12.

This script updates the page-1--6 bundle.  It does not treat characterization
paragraphs as new syntheses, does not associate the linked SI, and preserves
the printed solvent/name/number inconsistencies for review.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_5bda04d0a1e86544"
RUN = ROOT / "runs/current"
BUNDLE_PATH = RUN / "submitted.json"
BASE_BUNDLE_PATH = RUN / "standard_records.json"
DATASET_PATH = RUN / "dataset.json"


def ev(*pages: int) -> list[str]:
    out: list[str] = []
    for page in pages:
        out += [f"{DOC}:page{page}:image", f"{DOC}:page{page}:text"]
    return out


def uniq(values):
    return list(dict.fromkeys(values))


def src(raw, normalized=None, pages=(7,), status="pending_review", alternatives=()):
    if normalized is None:
        normalized = copy.deepcopy(raw)
    return st.fact(raw, normalized, ev(*pages), method="original_extraction", status=status, alternatives=alternatives)


def infer(normalized, pages=(7,), derivation="Inferred from the cited source drawing/text; expert review remains required."):
    return st.fact(None, normalized, ev(*pages), method="inference", derivation=derivation)


def derived(raw, normalized, pages=(7,), derivation="Candidate graph derived from the cited source; parsing is not source verification."):
    return st.fact(raw, normalized, ev(*pages), method="structure_derivation", derivation=derivation)


def inherited(raw, normalized=None, pages=(8,), derivation="Inherited from the explicitly stated general procedure."):
    if normalized is None:
        normalized = copy.deepcopy(raw)
    return st.fact(raw, normalized, ev(*pages), method="inherited_general_condition", derivation=derivation)


def missing(pages=(7,), status="unreported"):
    return st.fact(evidence=ev(*pages), status=status)


def q(number, unit, **extra):
    value = {"number": number, "unit": unit}
    value.update(extra)
    return value


def mid(key):
    return f"{DOC}:mol_{key}"


# The incremental source is the last successfully merged bundle.  The output is
# kept separate until validation succeeds.
bundle = json.loads(BASE_BUNDLE_PATH.read_text())
dataset = json.loads(DATASET_PATH.read_text())
bundle["documents"] = copy.deepcopy(dataset["documents"])
molecules = {m["molecule_id"]: m for m in bundle["molecules"]}
reactions = {r["reaction_id"]: r for r in bundle["reactions"]}

# Re-running after a successful merge is intentionally a no-op: the merged
# standard bundle is already the authoritative result, and records must not be
# duplicated merely because this provenance script is invoked again.
if f"{DOC}:scheme5_proposed_mechanism" in reactions:
    st.write(BUNDLE_PATH, bundle)
    print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]), "already_applied": True}, ensure_ascii=False))
    raise SystemExit(0)


def add_evidence(obj, *pages):
    obj["evidence_ids"] = uniq(obj.get("evidence_ids", []) + ev(*pages))


def merge_fact_evidence(fact, *pages):
    fact["evidence_ids"] = uniq(fact.get("evidence_ids", []) + ev(*pages))


def add_molecule(key, label, name, pages, *, smiles=None, representation_type="unresolved",
                 chirality="unresolved", geometry="not_applicable", representation=None,
                 mixture=None, identity_note=None):
    molecule_id = mid(key)
    if molecule_id in molecules:
        add_evidence(molecules[molecule_id], *pages)
        return molecule_id
    evidence = ev(*pages)
    m = st.new_molecule(molecule_id, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"] = [{"document_id": DOC, "compound_label": src(label, label, pages), "evidence_ids": evidence}]
    m["names"]["original"] = src(name, name, pages)
    m["names"]["english"] = src(name, name, pages)
    m["names"]["abbreviations"] = src([label], [label], pages)
    m["names"]["aliases"] = src([], [], pages)
    m["cas_numbers"] = missing(pages)
    m["generic_template"] = None
    if smiles:
        diag = st.stereo_inspect(smiles)
        canonical = diag["canonical_isomeric_smiles"]
        m["representation_type"] = "smiles"
        m["smiles"] = derived(f"candidate SMILES for {label}", canonical, pages)
        m["stereochemistry"]["isomeric_smiles"] = derived(f"canonical isomeric SMILES for {label}", canonical, pages)
        m["stereochemistry"]["identity_key"] = derived(f"standard stereo graph key for {label}", diag["stereo_graph_key"], pages)
        representation = representation or ("partial" if any(x["specified"] == "Unspecified" for x in diag["potential_stereo"]) else "encoded")
        m["stereochemistry"]["atom_centers"] = infer(diag["atom_centers"], pages, "Computed graph diagnostic; not a source assignment.")
        m["stereochemistry"]["double_bonds"] = infer(diag["double_bonds"], pages, "Computed graph diagnostic; not a source assignment.")
    else:
        m["representation_type"] = representation_type
        m["smiles"] = missing(pages)
        m["stereochemistry"]["isomeric_smiles"] = missing(pages)
        m["stereochemistry"]["identity_key"] = missing(pages)
        m["stereochemistry"]["atom_centers"] = missing(pages)
        m["stereochemistry"]["double_bonds"] = missing(pages)
        representation = representation or "not_encoded"
    m["stereochemistry"]["axial"] = missing(pages)
    m["stereochemistry"]["mixture"] = src(mixture, mixture, pages) if mixture is not None else missing(pages)
    m["stereochemistry"]["assignment"] = missing(pages)
    m["stereochemistry"]["chirality_status"] = infer(chirality, pages, "Center/axis status reviewed from the source and candidate graph; omission is not interpreted as racemic.")
    m["stereochemistry"]["geometry_status"] = infer(geometry, pages, "E/Z and ring cis/trans status reviewed from the source and candidate graph.")
    m["stereochemistry"]["representation_status"] = infer(representation, pages, "Reports whether all source stereochemical information is represented in the retained graph.")
    m["identity_note"] = src(identity_note or ("Candidate structure retained for review." if smiles else "Structure not encoded; source label/name retained."),
                             identity_note or ("Candidate structure retained for review." if smiles else "Structure not encoded; source label/name retained."), pages)
    m["external_properties"] = []
    bundle["molecules"].append(m)
    molecules[molecule_id] = m
    return molecule_id


# Explicit materials from the page-8 procedures.  CAS numbers remain missing.
M_PTOL = add_molecule("ptoluidine", "p-toluidine", "p-toluidine", (8,), smiles="Cc1ccc(N)cc1", chirality="achiral")
M_TSCL = add_molecule("tscl", "TsCl", "p-toluenesulfonyl chloride", (8,), smiles="Cc1ccc(S(=O)(=O)Cl)cc1", chirality="achiral")
M_TEA = add_molecule("tea", "triethylamine", "triethylamine", (8,), smiles="CCN(CC)CC", chirality="achiral")
M_DCM = add_molecule("dcm", "DCM", "dichloromethane", (8,), smiles="ClCCl", chirality="achiral")
M_WATER = add_molecule("water", "water", "water", (8,), smiles="O", chirality="achiral")
M_NA2SO4 = add_molecule("na2so4", "Na2SO4", "sodium sulfate", (8, 10), smiles="[Na+].[Na+].[O-]S(=O)(=O)[O-]", chirality="achiral")
M_6A = mid("6a")
M_PROPBROM = add_molecule("propargyl_bromide", "3-bromoprop-1-yne", "3-bromoprop-1-yne", (8,), smiles="C#CCBr", chirality="achiral")
M_KI = add_molecule("ki", "KI", "potassium iodide", (8,), smiles="[K+].[I-]", chirality="achiral")
M_ETOA = add_molecule("etoac", "EtOAc", "ethyl acetate", (8, 10), smiles="CCOC(C)=O", chirality="achiral")
M_PROPARGYL = add_molecule("propargyl_sulfonamide", "propargyl precursor", "4-methyl-N-(prop-2-yn-1-yl)-N-(p-tolyl)benzenesulfonamide", (8,),
                           smiles="C#CCN(c1ccc(C)cc1)S(=O)(=O)c1ccc(C)cc1", chirality="achiral")
M_TBUOK = add_molecule("tbuok", "tBuOK", "potassium tert-butoxide", (8,), smiles="CC(C)(C)[O-].[K+]", chirality="achiral")
M_THF = add_molecule("thf", "THF", "tetrahydrofuran", (8,), smiles="C1CCOC1", chirality="achiral")
M_SILICA = add_molecule("silica", "silica gel", "silica gel", (7, 8, 10), representation_type="unresolved", chirality="not_applicable", mixture="Merck silica gel, 40-63 μm / unspecified flash silica")
M_PE = add_molecule("petroleum_ether", "PE", "petroleum ether", (8, 10), representation_type="unresolved", chirality="unresolved", mixture="petroleum ether; composition not reported")
M_OHBOH = add_molecule("2hydroxyphenylboronic_acid", "(2-hydroxyphenyl)boronic acid", "(2-hydroxyphenyl)boronic acid", (8,), smiles="OB(O)c1ccccc1O", chirality="achiral")
M_PINACOL = add_molecule("pinacol", "pinacol", "pinacol", (8,), smiles="CC(O)(C)C(C)(C)O", chirality="achiral")
M_BPINPHENOL = add_molecule("bpin_phenol", "2-(Bpin)phenol", "2-(4,4,5,5-tetramethyl-1,3,2-dioxaborolan-2-yl)phenol", (8,),
                            smiles="CC1(C)OB(c2ccccc2O)OC1(C)C", chirality="achiral")
M_DIPEA = mid("dipea")
M_TF2O = add_molecule("triflic_anhydride", "trifluoromethanesulfonic anhydride", "trifluoromethanesulfonic anhydride", (8,),
                      smiles="O=S(=O)(OS(=O)(=O)C(F)(F)F)C(F)(F)F", chirality="achiral")
M_BRINE = add_molecule("brine", "brine", "brine", (8,), representation_type="unresolved", chirality="not_applicable", mixture="aqueous brine; concentration not reported")
M_HEXANE = add_molecule("nhexane", "n-hexane", "n-hexane", (8,), smiles="CCCCCC", chirality="achiral")
M_1G = add_molecule("1g_preparation", "1g", "2-(4,4,5,5-tetramethyl-1,3,2-dioxaborolan-2-yl)phenyl trifluoromethanesulfonate", (8,),
                    smiles="CC1(C)OB(c2ccccc2OS(=O)(=O)C(F)(F)F)OC1(C)C", chirality="achiral")


def participant(pid, molecule_id, role, pages, *, stage="reaction", mass=None, amount=None, volume=None, equiv=None, molpct=None, ref=None, raw_volume=None):
    p = st.new_participant(pid, molecule_id, role, ev(*pages), stage=stage)
    if mass is not None:
        p["charges"]["mass"] = src(f"{mass} g" if mass >= 0.1 else f"{mass} g", q(mass, "g"), pages)
    if amount is not None:
        p["charges"]["amount"] = src(f"{amount} mmol", q(amount, "mmol"), pages)
    if volume is not None:
        raw = raw_volume or f"{volume} mL"
        p["charges"]["volume"] = src(raw, q(volume, "mL"), pages)
    if equiv is not None:
        p["charges"]["equivalents"] = src(f"{equiv} equiv", q(equiv, "equiv"), pages)
    if molpct is not None:
        p["charges"]["mol_percent"] = src(f"{molpct} mol %", q(molpct, "mol%"), pages)
    if ref is not None:
        p["equivalence_reference"] = infer(ref, pages, "Denominator follows the source-stated substrate scale.")
    return p


def step(step_id, stage, pids, pages, operation, *, temp=None, time=None, atmosphere=None, pressure=None, inherited_values=False):
    s = st.new_step(step_id, ev(*pages), stage=stage)
    s["participant_ids"] = pids
    fact = inherited if inherited_values else src
    s["operation"] = src(operation, operation, pages)
    s["temperature"] = fact(temp[0], q(temp[1], temp[2]), pages) if temp else missing(pages)
    s["time"] = fact(time[0], q(time[1], time[2]), pages) if time else missing(pages)
    s["atmosphere"] = fact(atmosphere, atmosphere, pages) if atmosphere else missing(pages)
    # The controlled vocabulary has no pressure units.  Preserve numeric torr
    # expressions as unresolved source text instead of silently converting units.
    s["pressure"] = fact(pressure[0], pressure[0], pages) if pressure else missing(pages)
    return s


def classification(kind, pages):
    table = {
        "ns": ([{"change_id": "b1", "atoms": [{"element": "N", "source_role": "p-toluidine"}, {"element": "S", "source_role": "TsCl"}], "before_order": 0, "after_order": 1, "change_type": "formed"}],
               [{"category_number": 31, "code": "N-S", "supporting_change_ids": ["b1"]}]),
        "cn": ([{"change_id": "b1", "atoms": [{"element": "C", "source_role": "propargyl bromide"}, {"element": "N", "source_role": "sulfonamide"}], "before_order": 0, "after_order": 1, "change_type": "formed"}],
               [{"category_number": 5, "code": "C-N", "supporting_change_ids": ["b1"]}]),
        "isomerization": ([
            {"change_id": "b1", "atoms": [{"element": "C", "source_role": "terminal propargyl carbon"}, {"element": "H", "source_role": "intramolecularly transferred hydrogen"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain"},
            {"change_id": "b2", "atoms": [{"element": "C", "source_role": "N-adjacent propargyl carbon"}, {"element": "H", "source_role": "transferred hydrogen"}], "before_order": 1, "after_order": 0, "change_type": "hydrogen_loss"},
            {"change_id": "b3", "atoms": [{"element": "C"}, {"element": "C"}], "before_order": 3, "after_order": 2, "change_type": "order_decreased"},
            {"change_id": "b4", "atoms": [{"element": "C"}, {"element": "C"}], "before_order": 1, "after_order": 2, "change_type": "order_increased"}],
            [{"category_number": 1, "code": "C-H", "supporting_change_ids": ["b1"]}, {"category_number": 13, "code": "DB-C-C", "supporting_change_ids": ["b4"]}]),
        "bos": ([
            {"change_id": "b1", "atoms": [{"element": "B", "source_role": "boronic acid"}, {"element": "O", "source_role": "pinacol"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b2", "atoms": [{"element": "B", "source_role": "boronic acid"}, {"element": "O", "source_role": "pinacol"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b3", "atoms": [{"element": "O", "source_role": "phenol"}, {"element": "S", "source_role": "triflic anhydride"}], "before_order": 0, "after_order": 1, "change_type": "formed"}],
            [{"category_number": 41, "code": "other-B-O", "supporting_change_ids": ["b1", "b2"]}, {"category_number": 26, "code": "O-S", "supporting_change_ids": ["b3"]}]),
        "hydroarylation": ([
            {"change_id": "b1", "atoms": [{"element": "C", "source_role": "aryl boronic acid"}, {"element": "C", "source_role": "allene carbon"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b2", "atoms": [{"element": "C", "source_role": "allene carbon"}, {"element": "H", "source_role": "unresolved hydrogen donor"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain"}],
            [{"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1"]}, {"category_number": 1, "code": "C-H", "supporting_change_ids": ["b2"]}]),
    }
    changes, labels = table.get(kind, ([], []))
    coverage = "partial" if changes else "unresolved"
    return {"standard_id": st.STANDARD, "mapped_reaction_smiles": missing(pages),
            "bond_changes": [derived(f"inferred net bond change {c['change_id']}", c, pages, "Net change inferred from source reactant/product identity; no atom mapping asserted.") for c in changes],
            "labels": infer(labels, pages, "Labels computed only from listed newly formed bonds; unchanged/lowered bonds are not labels."),
            "coverage": infer(coverage, pages, "No source-verified atom-mapped reaction is available."),
            "limitations": ["Bond changes are inferred and require expert review.", "Coproducts and complete atom correspondence are not asserted."], "model_votes": []}


def new_reaction(key, pages, name, transformation, categories, class_kind, record_type="specific_experiment", unresolved=()):
    rid = f"{DOC}:{key}"
    r = st.new_reaction(rid, DOC, ev(*pages))
    r["record_type"] = record_type
    r["experiment_categories"] = src(list(categories), list(categories), pages)
    r["location"] = src({"pdf_pages": list(pages), "source_region": key}, {"pdf_pages": list(pages), "source_region": key}, pages)
    r["reaction_name"] = src(name, name, pages)
    r["transformation_type"] = infer(transformation, pages)
    r["classification"] = classification(class_kind, pages)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
                   "unresolved": list(unresolved) or ["Complete source-verified structures and atom mapping remain unavailable."]}
    r["completeness"] = "partial"
    bundle["reactions"].append(r)
    reactions[rid] = r
    return r


def add_yield(r, oid, target, raw, number, pages, *, basis="isolated", product_mass=None, scope="single_experiment"):
    o = st.new_outcome(oid, "yield", [target], ev(*pages))
    o["product_role"] = src("product", "product", pages)
    o["value"] = src(raw, q(number, "%"), pages)
    o["yield_type"] = src("isolated yield", "isolated", pages) if basis == "isolated" else missing(pages)
    o["measurement_basis"] = src("isolated yield", "isolated", pages) if basis == "isolated" else missing(pages)
    o["scope"] = src(scope, scope, pages)
    r["outcomes"].append(o)
    if product_mass:
        m = st.new_outcome(oid + "_product_mass", "product_amount", [target], ev(*pages))
        m["product_role"] = src("isolated product", "isolated_product", pages)
        m["value"] = src(f"{product_mass[0]} {product_mass[1]}", q(product_mass[0], product_mass[1]), pages)
        m["measurement_basis"] = src("isolated product mass", "isolated_mass", pages)
        m["scope"] = src(scope, scope, pages)
        r["outcomes"].append(m)


def null_conditions(r, pages):
    r["reaction_scale"] = missing(pages)
    r["solvent_conditions"]["composition"] = missing(pages)
    r["solvent_conditions"]["total_volume"] = missing(pages)


# Page 7: Scheme 5 is explicitly proposed, not an observed experiment.
mech = new_reaction("scheme5_proposed_mechanism", (7,), "Proposed ligand-controlled palladium hydroarylation mechanism",
                    "proposed catalytic cycles for regioselective hydroarylation", ("mechanistic_study",), "hydroarylation",
                    record_type="proposed_mechanism",
                    unresolved=("Intermediates B-E and I-III are author-proposed and are not reported as isolated or directly observed species.",
                                "The boronic-acid hydrogen-source assignment is tentative."))
mech["general_scheme_ids"] = [f"{DOC}:generic_anti_markovnikov", f"{DOC}:generic_markovnikov"]
mech["procedure"]["original_text"] = src("Scheme 5. Possible Mechanism", "author-proposed mechanism", (7,))
mech["notes"].append(src("Cycle I invokes MeCN/phenyl/allene coordination and gives 3a; cycle II invokes dppf and gives 4a.",
                              "Conceptual catalytic cycles; no separate experimental outcome.", (7,)))
null_conditions(mech, (7,))


# Page 8: three discrete preparations leading to allenamide 1a.
M_1A = mid("1a")
prep1 = new_reaction("prep_1a_stage1_sulfonamide", (8,), "Preparation of 4-methyl-N-(p-tolyl)benzenesulfonamide",
                     "sulfonamide formation", ("other",), "ns",
                     unresolved=("The phrase 'organic phase was combined' is retained although only one organic phase is described before washing.",))
prep1["participants"] = [
    participant("p1_ptol", M_PTOL, "reactant", (8,), mass=3.2, amount=30),
    participant("p1_tscl", M_TSCL, "reactant", (8,), mass=6.3, amount=33),
    participant("p1_tea", M_TEA, "reagent", (8,), mass=9.1, amount=90),
    participant("p1_dcm1", M_DCM, "solvent", (8,), volume=15),
    participant("p1_dcm2", M_DCM, "solvent", (8,), volume=15),
    participant("p1_water", M_WATER, "workup_agent", (8,), stage="workup"),
    participant("p1_drier", M_NA2SO4, "workup_agent", (8,), stage="workup"),
    participant("p1_prod", M_6A, "product", (8,)),
]
prep1["reaction_scale"] = src("p-toluidine 30 mmol", {"reference_molecule_id": M_PTOL, "amount": q(30, "mmol")}, (8,))
prep1["solvent_conditions"]["composition"] = src("DCM", [M_DCM], (8,))
prep1["solvent_conditions"]["total_volume"] = src("15 mL + 15 mL", q(30, "mL", components=[15, 15]), (8,))
prep1["steps"] = [
    step("initial_charge", "reaction", ["p1_ptol", "p1_tscl", "p1_dcm1"], (8,), "TsCl added to p-toluidine solution", temp=("room temperature", 25, "°C")),
    step("dropwise_base_addition", "reaction", ["p1_tea", "p1_dcm2"], (8,), "triethylamine solution added dropwise through a constant-pressure funnel under vigorous stirring", temp=("room temperature", 25, "°C")),
    step("overnight_hold", "reaction", [], (8,), "stirred overnight", temp=("room temperature", 25, "°C")),
    step("aqueous_wash", "workup", ["p1_water"], (8,), "organic phase washed three times with water"),
    step("drying", "workup", ["p1_drier"], (8,), "organic phase dried over Na2SO4"),
]
add_yield(prep1, "sulfonamide_yield", M_6A, "91%", 91, (8,), product_mass=(7.1, "g"))
prep1["notes"].append(src("white solid", "white solid", (8,)))
prep1["notes"].append(src("overnight", "qualitative duration: overnight", (8,)))

prep2 = new_reaction("prep_1a_stage2_propargylation", (8,), "Preparation of the propargyl sulfonamide precursor to 1a",
                     "N-propargylation", ("other",), "cn",
                     unresolved=("The source says 'stirred for 2 h and refluxed at 60 °C'; the ordering/relationship of these two holds is grammatically ambiguous.",))
prep2["participants"] = [
    participant("p2_6a", M_6A, "reactant", (8,), mass=1.83, amount=7),
    participant("p2_mecn", mid("mecn"), "solvent", (8,), volume=40),
    participant("p2_br", M_PROPBROM, "reactant", (8,), volume=0.7, amount=9.1),
    participant("p2_k2co3", mid("k2co3"), "reagent", (8,), mass=1.94, amount=14),
    participant("p2_ki", M_KI, "reagent", (8,), mass=0.58, amount=3.5),
    participant("p2_water", M_WATER, "workup_agent", (8,), stage="workup"),
    participant("p2_etoac", M_ETOA, "workup_agent", (8,), stage="workup"),
    participant("p2_drier", M_NA2SO4, "workup_agent", (8,), stage="workup"),
    participant("p2_silica", M_SILICA, "purification_agent", (8,), stage="purification"),
    participant("p2_pe", M_PE, "purification_agent", (8,), stage="purification"),
    participant("p2_ea", M_ETOA, "purification_agent", (8,), stage="purification"),
    participant("p2_prod", M_PROPARGYL, "product", (8,)),
]
prep2["reaction_scale"] = src("sulfonamide 7 mmol", {"reference_molecule_id": M_6A, "amount": q(7, "mmol")}, (8,))
prep2["solvent_conditions"]["composition"] = src("MeCN", [mid("mecn")], (8,))
prep2["solvent_conditions"]["total_volume"] = src("40 mL", q(40, "mL"), (8,))
prep2["steps"] = [
    step("successive_charge", "reaction", ["p2_6a", "p2_mecn", "p2_br", "p2_k2co3", "p2_ki"], (8,), "3-bromoprop-1-yne, K2CO3, and KI added successively"),
    step("two_hour_hold", "reaction", [], (8,), "stirred for 2 h", time=("2 h", 2, "h")),
    step("heated_hold", "reaction", [], (8,), "refluxed at 60 °C", temp=("60 °C", 60, "°C")),
    step("quench_extract", "workup", ["p2_water", "p2_etoac"], (8,), "quenched with water; aqueous layer extracted three times with EtOAc"),
    step("dry_concentrate", "workup", ["p2_drier"], (8,), "dried over Na2SO4 and concentrated in vacuo"),
    step("column", "purification", ["p2_silica", "p2_pe", "p2_ea"], (8,), "silica gel column chromatography, PE/EA = 5:1"),
]
add_yield(prep2, "propargyl_precursor_yield", M_PROPARGYL, "89%", 89, (8,), product_mass=(1.86, "g"))
prep2["notes"].append(src("white solid", "white solid", (8,)))

prep3 = new_reaction("prep_1a_stage3_isomerization", (8,), "Preparation of allenamide 1a from the propargyl precursor",
                     "base-promoted propargyl-to-allene isomerization", ("other",), "isomerization",
                     unresolved=("No post-addition reaction-hold duration, quench, extraction solvent, or workup volumes are reported in the inspected text.",
                                 "The text jumps from 30 min dropwise addition to combined organic phase/drying."))
prep3["participants"] = [
    participant("p3_tbuok", M_TBUOK, "reagent", (8,), mass=0.67, amount=6, equiv=1.5, ref=M_PROPARGYL),
    participant("p3_thf1", M_THF, "solvent", (8,), volume=10),
    participant("p3_sub", M_PROPARGYL, "reactant", (8,), mass=1.2, amount=4),
    participant("p3_thf2", M_THF, "solvent", (8,), volume=5),
    participant("p3_drier", M_NA2SO4, "workup_agent", (8,), stage="workup"),
    participant("p3_silica", M_SILICA, "purification_agent", (8,), stage="purification"),
    participant("p3_pe", M_PE, "purification_agent", (8,), stage="purification"),
    participant("p3_ea", M_ETOA, "purification_agent", (8,), stage="purification"),
    participant("p3_prod", M_1A, "product", (8,)),
]
prep3["reaction_scale"] = src("propargyl sulfonamide 4 mmol", {"reference_molecule_id": M_PROPARGYL, "amount": q(4, "mmol")}, (8,))
prep3["solvent_conditions"]["composition"] = src("dry THF", [M_THF], (8,))
prep3["solvent_conditions"]["total_volume"] = src("10 mL + 5 mL", q(15, "mL", components=[10, 5]), (8,))
prep3["steps"] = [
    step("base_prestir", "reaction", ["p3_tbuok", "p3_thf1"], (8,), "tBuOK solution stirred under nitrogen", temp=("-5 °C", -5, "°C"), atmosphere="N2"),
    step("dropwise_addition", "reaction", ["p3_sub", "p3_thf2"], (8,), "substrate solution added dropwise", temp=("-5 °C", -5, "°C"), time=("over 30 min", 30, "min"), atmosphere="N2"),
    step("dry_evaporate", "workup", ["p3_drier"], (8,), "organic phase combined, dried over Na2SO4, and evaporated"),
    step("flash", "purification", ["p3_silica", "p3_pe", "p3_ea"], (8,), "flash chromatography, PE/EA = 10:1"),
]
add_yield(prep3, "allenamide_1a_yield", M_1A, "78%", 78, (8,), product_mass=(0.93, "g"))


# Page 8: a two-stage procedure with one explicitly overall yield.
prep1g = new_reaction("prep_1g_two_stage", (8,), "Two-stage preparation of aryl triflate 1g",
                      "boronate ester formation followed by phenolic triflation", ("other",), "bos",
                      unresolved=("The 86% yield is explicitly overall from (2-hydroxyphenyl)boronic acid and is not assigned to the second stage alone.",
                                  "Brine concentration and crude intermediate amount are not reported."))
prep1g["participants"] = [
    participant("g_acid", M_OHBOH, "reactant", (8,), mass=1.10, amount=8.0),
    participant("g_pinacol", M_PINACOL, "reactant", (8,), mass=1.04, amount=8.8),
    participant("g_dcm1", M_DCM, "solvent", (8,), volume=15),
    participant("g_water1", M_WATER, "workup_agent", (8,), stage="workup", volume=10, raw_volume="ca. 10 mL"),
    participant("g_dcm_extract1", M_DCM, "workup_agent", (8,), stage="workup", volume=30, raw_volume="ca. 30 mL × 3"),
    participant("g_brine", M_BRINE, "workup_agent", (8,), stage="workup", volume=30, raw_volume="ca. 30 mL × 1"),
    participant("g_drier1", M_NA2SO4, "workup_agent", (8,), stage="workup"),
    participant("g_intermediate", M_BPINPHENOL, "unknown", (8,)),
    participant("g_dcm2", M_DCM, "solvent", (8,), volume=25),
    participant("g_dipea", M_DIPEA, "reagent", (8,), volume=1.57, amount=9),
    participant("g_tf2o", M_TF2O, "reactant", (8,), volume=1.21, amount=7.2),
    participant("g_water2", M_WATER, "workup_agent", (8,), stage="workup", volume=10, raw_volume="ca. 10 mL"),
    participant("g_dcm_extract2", M_DCM, "workup_agent", (8,), stage="workup", volume=10, raw_volume="ca. 10 mL × 3"),
    participant("g_drier2", M_NA2SO4, "workup_agent", (8,), stage="workup"),
    participant("g_silica", M_SILICA, "purification_agent", (8,), stage="purification"),
    participant("g_hexane", M_HEXANE, "purification_agent", (8,), stage="purification"),
    participant("g_ea", M_ETOA, "purification_agent", (8,), stage="purification"),
    participant("g_prod", M_1G, "product", (8,)),
]
prep1g["reaction_scale"] = src("(2-hydroxyphenyl)boronic acid 8.0 mmol", {"reference_molecule_id": M_OHBOH, "amount": q(8.0, "mmol")}, (8,))
prep1g["solvent_conditions"]["composition"] = src("DCM (stage 1: 15 mL; stage 2: 25 mL)", [M_DCM], (8,))
prep1g["solvent_conditions"]["total_volume"] = missing((8,), status="not_applicable")
prep1g["steps"] = [
    step("stage1_charge", "reaction", ["g_acid", "g_pinacol", "g_dcm1"], (8,), "pinacol added to boronic acid solution", temp=("room temperature", 25, "°C")),
    step("stage1_hold", "reaction", [], (8,), "stirred at the same temperature", temp=("room temperature", 25, "°C"), time=("2 h", 2, "h")),
    step("stage1_workup", "workup", ["g_water1", "g_dcm_extract1", "g_brine", "g_drier1"], (8,), "water addition, DCM extraction (ca. 30 mL × 3), brine wash, drying, filtration, concentration"),
    step("stage1_crude", "reaction", ["g_intermediate"], (8,), "crude 2-(Bpin)phenol carried forward without further purification"),
    step("stage2_charge", "reaction", ["g_intermediate", "g_dcm2", "g_dipea", "g_tf2o"], (8,), "DIPEA and triflic anhydride added", temp=("-78 °C", -78, "°C")),
    step("stage2_hold", "reaction", [], (8,), "stirred at the same temperature", temp=("-78 °C", -78, "°C"), time=("3 h", 3, "h")),
    step("stage2_workup", "workup", ["g_water2", "g_dcm_extract2", "g_drier2"], (8,), "water addition, DCM extraction (ca. 10 mL × 3), drying, filtration, concentration"),
    step("short_pad", "purification", ["g_silica", "g_hexane", "g_ea"], (8,), "short silica pad, n-hexane/ethyl acetate = 10:1"),
    step("kugelrohr", "purification", [], (8,), "Kugelrohr distillation", temp=("ca. 150 °C", 150, "°C"), pressure=("2.0 torr", 2.0, "torr")),
]
add_yield(prep1g, "overall_yield", M_1G, "86% from (2-hydroxyphenyl)boronic acid", 86, (8,), product_mass=(1.82, "g"), scope="overall_two_stage_yield")
prep1g["notes"].append(src("colorless oil", "colorless oil", (8,)))


# General procedures: attach exact charges/ordered operations to the existing parents.
def replace_general_procedure(rid, branch):
    r = reactions[rid]
    pages = (8,) if branch == "anti" else (10,)
    add_evidence(r, *pages)
    merge_fact_evidence(r["location"], *pages)
    substrate, aryl = mid("1a"), mid("2a")
    prefix = "gp3" if branch == "anti" else "gp4"
    catalyst = mid("pdmecn") if branch == "anti" else mid("pdoac2")
    cat_mass = 5.0
    ligand = None if branch == "anti" else mid("dppf")
    time_h = 1.0 if branch == "anti" else 0.5
    r["participants"] = [
        participant(prefix+"_sub", substrate, "reactant", pages, amount=0.20),
        participant(prefix+"_aryl", aryl, "reactant", pages, amount=0.24),
        participant(prefix+"_cat", catalyst, "catalyst", pages, amount=0.02, molpct=10, ref=substrate),
        participant(prefix+"_base", mid("k3po4"), "reagent", pages, mass=0.085, amount=0.40, equiv=2.0, ref=substrate),
        participant(prefix+"_solv", mid("dioxane"), "solvent", pages, volume=2.0),
        participant(prefix+"_etoac", M_ETOA, "workup_agent", pages, stage="workup", volume=10, raw_volume="3 × 10 mL"),
        participant(prefix+"_drier", M_NA2SO4, "workup_agent", pages, stage="workup"),
        participant(prefix+"_silica", M_SILICA, "purification_agent", pages, stage="purification"),
    ]
    # Catalyst mass is explicitly mg, not the default g used by the helper.
    r["participants"][2]["charges"]["mass"] = src("5.0 mg", q(cat_mass, "mg"), pages)
    r["participants"][3]["charges"]["mass"] = src("85.0 mg", q(85.0, "mg"), pages)
    if ligand:
        lp = participant(prefix+"_lig", ligand, "ligand", pages, amount=0.02, molpct=10, ref=substrate)
        lp["charges"]["mass"] = src("11.0 mg", q(11.0, "mg"), pages)
        r["participants"].insert(3, lp)
    r["reaction_scale"] = src("allenamide 0.20 mmol", {"reference_molecule_id": substrate, "amount": q(0.20, "mmol")}, pages)
    if branch == "anti":
        conflict = st.fact("toluene 1,4-dioxane (2 mL)", None, ev(3, 8), method="original_extraction", status="conflict")
        conflict["alternatives"] = [src("1,4-dioxane (2 mL)", [mid("dioxane")], (3,)), st.fact("toluene 1,4-dioxane (2 mL)", None, ev(8), method="original_extraction", status="unreadable")]
        r["solvent_conditions"]["composition"] = conflict
        r["review"]["unresolved"] = uniq(r["review"]["unresolved"] + ["Page 8 prints 'Then, toluene 1,4-dioxane (2 mL) was added', whereas Table 2 specifies 1,4-dioxane; the solvent wording is retained as a conflict."])
        r["notes"].append(src("Table 2: 1,4-dioxane; page-8 procedure: 'toluene 1,4-dioxane'.", "Solvent source conflict retained without editorial correction.", (3, 8)))
    else:
        r["solvent_conditions"]["composition"] = src("1,4-dioxane", [mid("dioxane")], pages)
    r["solvent_conditions"]["total_volume"] = src("2 mL", q(2, "mL"), pages)
    reaction_ids = [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"]
    r["steps"] = [
        step("charge", "reaction", reaction_ids, pages, "components added to a reaction tube; solvent added under nitrogen", atmosphere="N2"),
        step("heated_hold", "reaction", [], pages, "stirred in an oil bath", temp=("50 °C", 50, "°C"), time=(f"{time_h:g} h", time_h, "h"), atmosphere="N2"),
        step("cool_extract", "workup", [prefix+"_etoac"], pages, "cooled to room temperature; aqueous layer extracted with EtOAc (3 × 10 mL)", temp=("room temperature", 25, "°C")),
        step("dry_concentrate", "workup", [prefix+"_drier"], pages, "combined organic phase dried over Na2SO4 and concentrated in vacuo"),
        step("column", "purification", [prefix+"_silica"], pages, "silica gel column chromatography"),
    ]
    r["procedure"]["original_text"] = src("General Procedures for Products 3" if branch == "anti" else "General Procedures for Products 4",
                                             "general procedure with ordered reaction/workup/purification steps", pages)
    r["special_conditions"]["reactor"] = src("reaction tube; oil bath", "reaction tube heated in oil bath", pages)


replace_general_procedure(f"{DOC}:generic_anti_markovnikov", "anti")
replace_general_procedure(f"{DOC}:generic_markovnikov", "mark")


# Product-specific isolation/characterization paragraphs update existing scope
# experiments; spectra remain source context rather than new reaction records.
anti_data = {
    "3a": (8, 70.9, 86, "15:1", "yellow oil"), "3b": (8, 58.5, 71, "15:1", "yellow oil"),
    "3c": (8, 55.4, 70, "15:1", "yellow solid"), "3d": (8, 52.9, 64, "15:1", "yellow solid"),
    "3e": (8, 59.3, 75, "15:1", "yellow solid"), "3f": (9, 46.6, 51, "15:1", "yellow oil"),
    "3g": (9, 54.7, 52, "15:1", "yellow solid"), "3h": (9, 60.2, 55, "5:1", "yellow oil"),
    "3i": (9, 54.6, 67, "15:1", "yellow oil"), "3j": (9, 44.8, 53, "10:1", "yellow solid"),
    "3k": (9, 38.6, 48, "10:1", "yellow solid"), "3l": (9, 63.3, 71, "15:1", "yellow solid"),
    "3m": (9, 50.9, 61, "10:1", "yellow solid"), "3n": (9, 53.4, 48, "15:1", "yellow oil"),
    "3o": (9, 56.4, 63, "10:1", "yellow oil"), "3p": (9, 79.1, 81, "15:1", "yellow oil"),
    "3q": (10, 45.0, 52, "15:1", "yellow oil"), "3r": (10, 48.7, 55, "15:1", "yellow oil"),
    "3s": (10, 55.6, 65, "15:1", "yellow oil"), "3t": (10, 29.6, 46, "10:1", "yellow oil"),
}
mark_data = {
    "4a": (10, 66.7, 81, "15:1", "yellow oil"), "4b": (10, 56.0, 68, "15:1", "yellow oil"),
    "4c": (10, 66.9, 93, "15:1", "yellow solid"), "4d": (10, 56.8, 79, "15:1", "yellow solid"),
    "4e": (10, 64.5, 78, "15:1", "yellow solid"), "4f": (11, 69.4, 76, "10:1", "yellow oil"),
    "4g": (11, 46.5, 57, "10:1", "yellow oil"), "4h": (11, 68.5, 84, "10:1", "yellow oil"),
    "4i": (11, 77.2, 89, "15:1", "yellow oil"), "4j": (11, 61.3, 56, "10:1", "yellow oil"),
    "4k": (11, 47.0, 54, "5:1", "yellow oil"), "4l": (11, 68.0, 81, "5:1", "yellow oil"),
    "4m": (11, 45.6, 51, "10:1", "yellow oil"), "4n": (11, 41.0, 48, "15:1", "yellow solid"),
    "4o": (11, 81.9, 85, "15:1", "yellow solid"), "4p": (11, 76.9, 89, "10:1", "yellow oil"),
    "4q": (11, 79.6, 93, "10:1", "yellow oil"),
}


def update_scope_record(label, data, branch):
    page, mass_mg, yield_pct, ratio, form = data
    rid = f"{DOC}:table2_{label}" if branch == "anti" else f"{DOC}:table3_{label}"
    r = reactions[rid]
    general_page = 8 if branch == "anti" else 10
    add_evidence(r, general_page, page)
    location_pages = uniq(r["location"]["normalized_value"]["pdf_pages"] + [general_page, page])
    r["location"] = src({"pdf_pages": location_pages, "source_region": r["location"]["normalized_value"]["source_region"]},
                        {"pdf_pages": location_pages, "source_region": r["location"]["normalized_value"]["source_region"]}, tuple(location_pages))
    for s in r["steps"]:
        if s["stage"] == "reaction":
            s["time"] = inherited("1 h" if branch == "anti" else "0.5 h", q(1.0 if branch == "anti" else 0.5, "h"), (general_page,))
    prod = next(p["molecule_id"] for p in r["participants"] if p["role"] == "product")
    # Attach the product page to the existing isolated-yield assertion.
    for o in r["outcomes"]:
        if o["kind"] == "yield" and prod in o["target_molecule_ids"]:
            merge_fact_evidence(o["value"], page)
            o["yield_type"] = src("isolated yield", "isolated", (3 if branch == "anti" else 4, page))
            o["measurement_basis"] = src("isolated yield", "isolated", (3 if branch == "anti" else 4, page))
    amount = st.new_outcome("isolated_product_mass", "product_amount", [prod], ev(page))
    amount["product_role"] = src(label, "isolated_product", (page,))
    amount["value"] = src(f"{mass_mg} mg", q(mass_mg, "mg"), (page,))
    amount["measurement_basis"] = src("isolated product mass", "isolated_mass", (page,))
    amount["scope"] = src("product-specific isolation paragraph", "single_experiment", (page,))
    r["outcomes"].append(amount)
    r["notes"].append(src(f"column chromatography (PE/EA = {ratio}) on silica gel; {form}",
                          {"purification": "silica gel column chromatography", "eluent_ratio": f"PE/EA {ratio}", "appearance": form}, (page,)))
    r["procedure"]["general_procedure_reference"] = src("standard procedure", f"{DOC}:generic_anti_markovnikov" if branch == "anti" else f"{DOC}:generic_markovnikov", (page,))
    r["review"]["unresolved"] = [x for x in r["review"]["unresolved"] if "Reaction time is not reported" not in x]
    r["review"]["unresolved"] = uniq(r["review"]["unresolved"] + ["Product-specific NMR/HRMS data are characterization evidence and were not represented as separate synthetic experiments."])


for label, data in anti_data.items():
    update_scope_record(label, data, "anti")
for label, data in mark_data.items():
    update_scope_record(label, data, "mark")


# Product 3u in the page-10 characterization belongs to Scheme 4a, not the
# different Table-2 structure also printed as 3u.  Likewise 4v belongs to 4d.
def update_scheme_product(rid, page, product_mid, mass_mg, form, ratio):
    r = reactions[f"{DOC}:{rid}"]
    add_evidence(r, page)
    for o in r["outcomes"]:
        if o["kind"] == "yield" and product_mid in o["target_molecule_ids"]:
            o["yield_type"] = src("isolated yield", "isolated", (page,))
            o["measurement_basis"] = src("isolated yield", "isolated", (page,))
            merge_fact_evidence(o["value"], page)
    amount = st.new_outcome("isolated_product_mass", "product_amount", [product_mid], ev(page))
    amount["product_role"] = src("isolated product", "isolated_product", (page,))
    amount["value"] = src(f"{mass_mg} mg", q(mass_mg, "mg"), (page,))
    amount["measurement_basis"] = src("isolated product mass", "isolated_mass", (page,))
    amount["scope"] = src("product-specific isolation paragraph", "single_experiment", (page,))
    r["outcomes"].append(amount)
    r["notes"].append(src(f"column chromatography (PE/EA = {ratio}); {form}",
                          {"eluent_ratio": f"PE/EA {ratio}", "appearance": form}, (page,)))


update_scheme_product("scheme4_a", 10, mid("scheme4_3u"), 30.7, "yellow oil", "15:1")
update_scheme_product("scheme4_d", 11, mid("scheme4_4v"), 79.4, "yellow oil", "15:1")
add_evidence(reactions[f"{DOC}:scheme4_d"], 12)

# Page 7 completes Scheme 4c's negative observation.
s4c = reactions[f"{DOC}:scheme4_c"]
add_evidence(s4c, 7)
s4c["notes"].append(src("reaction of 1u did not perform well, and we could not isolate an obvious product", "trace was depicted on page 6; no obvious product was isolated per page 7 prose", (6, 7)))
s4c["review"]["unresolved"] = uniq(s4c["review"]["unresolved"] + ["Trace depiction and failure to isolate an obvious product are both retained; analytical basis for trace is not stated."])

# Preserve source naming/numbering problems rather than silently correcting.
reactions[f"{DOC}:table2_3q"]["notes"].append(src("The printed page-10 compound name for 3q describes a 2,6-dichloropyridin-3-yl product, duplicating the identity printed for 3o, while the scope drawing/analytical formula context differs.", "source name/structure inconsistency", (9, 10)))
reactions[f"{DOC}:table2_3q"]["review"]["unresolved"].append("Printed 3q name conflicts with the scope identity and appears to duplicate 3o; not editorially corrected.")
reactions[f"{DOC}:table3_4d"]["notes"].append(src("The printed page-10 name for 4d repeats the 2-fluorophenyl name used for 4c although the scope row distinguishes the products.", "source name/structure inconsistency", (10,)))
reactions[f"{DOC}:table3_4d"]["review"]["unresolved"].append("Printed 4d name duplicates 4c; retained without correction.")
reactions[f"{DOC}:table3_4o"]["notes"].append(src("The paragraph says 'provide 40' (zero) while the product is consistently labeled 4o elsewhere.", "source typographical product-label inconsistency", (11,)))
reactions[f"{DOC}:table3_4o"]["review"]["unresolved"].append("Printed '40' versus 4o label retained as a source typo/conflict.")


# Page 12: a distinct scale-up experiment.  Only explicitly stated scale and
# outcome are populated; other charges are not multiplied from the small scale.
scale = new_reaction("scaleup_3a_1mmol", (12,), "Larger-scale synthesis of 3a",
                     "anti-Markovnikov hydroarylation scale-up", ("scale_up",), "hydroarylation",
                     unresolved=("The paragraph reports a 1 mmol operation but does not restate reagent/catalyst/base/solvent quantities, temperature, or time.",
                                 "Conditions were not arithmetically scaled from the 0.2 mmol general procedure."))
scale["general_scheme_ids"] = [f"{DOC}:generic_anti_markovnikov"]
scale["participants"] = [
    participant("scale_1a", M_1A, "reactant", (12,), amount=1.0),
    participant("scale_2a", mid("2a"), "reactant", (12,)),
    participant("scale_cat", mid("pdmecn"), "catalyst", (12,)),
    participant("scale_base", mid("k3po4"), "reagent", (12,)),
    participant("scale_solv", mid("dioxane"), "solvent", (12,)),
    participant("scale_prod", mid("3a"), "product", (12,)),
]
scale["reaction_scale"] = src("1 mmol", {"reference_molecule_id": M_1A, "amount": q(1.0, "mmol")}, (12,))
scale["solvent_conditions"]["composition"] = missing((12,))
scale["solvent_conditions"]["total_volume"] = missing((12,))
scale["steps"] = [step("larger_scale_operation", "reaction", ["scale_1a", "scale_2a", "scale_cat", "scale_base", "scale_solv"], (12,), "larger-scale operation performed", inherited_values=False)]
scale["procedure"]["general_procedure_reference"] = src("larger-scale synthesis of 3a", f"{DOC}:generic_anti_markovnikov", (12,))
add_yield(scale, "scaleup_3a_yield", mid("3a"), "82% isolated yield", 82, (12,), product_mass=(337.8, "mg"))


# Page-level nonreaction material is retained as source context.
bundle["source_context"].extend([
    {"context_id": f"{DOC}:page7_general_information", "document_id": DOC, "kind": "experimental_general_information", "evidence_ids": ev(7),
     "note": "Flash silica specification and NMR/MS instrumentation context; not a reaction experiment."},
    {"context_id": f"{DOC}:pages8_12_characterization_inventory", "document_id": DOC, "kind": "characterization_inventory", "evidence_ids": ev(8, 9, 10, 11, 12),
     "note": "Product-specific 1H/13C/19F NMR and HRMS paragraphs support identity/characterization but were not converted into standalone syntheses, mass-spectrum experiments, or product yields beyond the explicitly printed isolation statements."},
    {"context_id": f"{DOC}:page12_si_notice_unassociated", "document_id": DOC, "kind": "supplementary_information_notice", "evidence_ids": ev(12),
     "note": "The paper states that SI is available online, but no SI file was explicitly supplied or associated in this workspace; no SI facts were imported."},
    {"context_id": f"{DOC}:pages12_references_start", "document_id": DOC, "kind": "bibliography", "evidence_ids": ev(12),
     "note": "Author information, acknowledgments, and the beginning of references are not current-paper experiments."},
])


# Raw inventory and standard-stereo diagnostics are deliberately separate from
# the formal bundle.
new_ids = [f"{DOC}:{x}" for x in ["scheme5_proposed_mechanism", "prep_1a_stage1_sulfonamide", "prep_1a_stage2_propargylation", "prep_1a_stage3_isomerization", "prep_1g_two_stage", "scaleup_3a_1mmol"]]
inventory = {
    "document_id": DOC,
    "batch_id": "allene-stereo-v31-200-20261017",
    "visually_inspected_pages": [7, 8, 9, 10, 11, 12],
    "new_reaction_ids": new_ids,
    "updated_scope_records": list(anti_data) + list(mark_data),
    "updated_existing_records": [f"{DOC}:generic_anti_markovnikov", f"{DOC}:generic_markovnikov", f"{DOC}:scheme4_a", f"{DOC}:scheme4_c", f"{DOC}:scheme4_d"],
    "notes": [
        "All machine records remain needs_review.",
        "Pages 8-12 characterization spectra were not counted as synthesis experiments.",
        "No SI document was associated.",
        "Printed solvent, name, and product-label conflicts were retained.",
    ],
}
st.write(RUN / "raw_model/pages_0007_0012_extraction.json", inventory)

diagnostics = {}
for key in ["ptoluidine", "tscl", "tea", "dcm", "water", "propargyl_bromide", "ki", "etoac", "propargyl_sulfonamide", "tbuok", "thf", "2hydroxyphenylboronic_acid", "pinacol", "bpin_phenol", "triflic_anhydride", "nhexane", "1g_preparation"]:
    m = molecules[mid(key)]
    if st.value(m["smiles"]):
        diagnostics[key] = st.stereo_inspect(st.value(m["smiles"]))
st.write(RUN / "raw_model/pages_0007_0012_stereo_diagnostics.json", diagnostics)
st.write(RUN / "raw_model/pages_0007_0012_mentions.json", st.source_mentions(dataset, ["General Procedures for Products 3", "General Procedures for Products 4", "Larger-Scale Synthesis", "3q", "4d", "4o"], context=300))


anti_page = {p: [f"{DOC}:table2_{label}" for label, data in anti_data.items() if data[0] == p] for p in range(8, 11)}
mark_page = {p: [f"{DOC}:table3_{label}" for label, data in mark_data.items() if data[0] == p] for p in range(10, 12)}
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 7, "state": "attempted", "regions": [
            {"region": "top-left continuation of Scheme 4 discussion", "finding": "Scheme 4c continuation states that the 1u/dppf reaction did not perform well and no obvious product could be isolated; merged with the page-6 experiment.", "evidence_ids": ev(7)},
            {"region": "Scheme 5 and mechanism discussion", "finding": "Author-proposed anti-Markovnikov and Markovnikov catalytic cycles; classified as proposed mechanism, not an observed experiment.", "evidence_ids": ev(7)},
            {"region": "Experimental Section, General Information", "finding": "Chromatography and analytical-instrument context; no standalone synthesis or spectral-peak experiment.", "evidence_ids": ev(7)},
        ], "record_ids": [f"{DOC}:scheme4_c", f"{DOC}:scheme5_proposed_mechanism"],
         "unresolved": ["Proposed intermediates are not directly observed or source-verified structures.", "Trace analytical basis for Scheme 4c is not stated."]},
        {"page": 8, "state": "attempted", "regions": [
            {"region": "Synthesis of 1a", "finding": "Three discrete operations recorded: sulfonamide formation, N-propargylation, and low-temperature propargyl-to-allene isomerization, with ordered additions/workups and isolated yields.", "evidence_ids": ev(8)},
            {"region": "Synthesis of 1g", "finding": "Two-stage boronate formation/triflation recorded as one overall experiment with the explicitly overall 86% yield.", "evidence_ids": ev(8)},
            {"region": "General Procedure for Products 3", "finding": "Exact reaction, 1 h hold, workup, and purification added to the generic anti-Markovnikov parent; solvent wording conflict preserved.", "evidence_ids": ev(8)},
            {"region": "Products 3a-3e", "finding": "Isolation masses, column ratios, appearances, and characterization context linked to existing scope records; spectra not made into new syntheses.", "evidence_ids": ev(8)},
        ], "record_ids": [f"{DOC}:prep_1a_stage1_sulfonamide", f"{DOC}:prep_1a_stage2_propargylation", f"{DOC}:prep_1a_stage3_isomerization", f"{DOC}:prep_1g_two_stage", f"{DOC}:generic_anti_markovnikov"] + anti_page[8],
         "unresolved": ["1a stage-3 post-addition hold and complete workup are not reported.", "Products-3 solvent wording conflicts with the Table-2 footnote."]},
        {"page": 9, "state": "attempted", "regions": [
            {"region": "Products 3e continuation and 3f-3p", "finding": "Isolation/characterization paragraphs linked to existing anti-Markovnikov scope experiments; no independent spectral experiments created.", "evidence_ids": ev(9)},
        ], "record_ids": [f"{DOC}:table2_3e"] + anti_page[9],
         "unresolved": ["Complete source-verified structures and atom mapping remain unavailable for most compact scope products."]},
        {"page": 10, "state": "attempted", "regions": [
            {"region": "Products 3p continuation and 3q-3t", "finding": "Isolation/characterization paragraphs linked to existing Table-2 records; printed 3q name conflict retained.", "evidence_ids": ev(10)},
            {"region": "Product 3u", "finding": "30.7 mg/35% isolation belongs to Scheme 4a cinnamyl-tethered 3u and remains distinct from the different Table-2 3u.", "evidence_ids": ev(10)},
            {"region": "General Procedure for Products 4", "finding": "Exact charges, 0.5 h hold, workup, and purification added to the generic Markovnikov parent.", "evidence_ids": ev(10)},
            {"region": "Products 4a-4e", "finding": "Isolation masses and characterization context linked to existing scope records; printed 4d naming conflict retained.", "evidence_ids": ev(10)},
        ], "record_ids": [f"{DOC}:table2_3p"] + anti_page[10] + [f"{DOC}:scheme4_a", f"{DOC}:generic_markovnikov"] + mark_page[10],
         "unresolved": ["Printed 3q identity duplicates 3o name while other context differs.", "Printed 4d name duplicates 4c."]},
        {"page": 11, "state": "attempted", "regions": [
            {"region": "Products 4f-4q", "finding": "Isolation masses, column ratios, appearances, and characterization context linked to existing Markovnikov scope records.", "evidence_ids": ev(11)},
            {"region": "Product 4v", "finding": "79.4 mg/81% isolation linked to Scheme 4d; characterization continues on page 12.", "evidence_ids": ev(11)},
        ], "record_ids": mark_page[11] + [f"{DOC}:scheme4_d"],
         "unresolved": ["Printed '40' in the 4o paragraph is retained as a product-label conflict.", "Most scope structures remain source-labeled rather than source-verified SMILES."]},
        {"page": 12, "state": "attempted", "regions": [
            {"region": "Product 4v continuation", "finding": "Remaining characterization data linked as source context to the existing Scheme 4d experiment.", "evidence_ids": ev(12)},
            {"region": "Larger-Scale Synthesis", "finding": "Distinct 1 mmol synthesis of 3a recorded with 337.8 mg and 82% isolated yield; unreported charges/conditions were not scaled by inference.", "evidence_ids": ev(12)},
            {"region": "Associated Content/author information/acknowledgments/references", "finding": "SI availability notice and bibliographic material only; no SI file associated and no new current-paper experiment created.", "evidence_ids": ev(12)},
        ], "record_ids": [f"{DOC}:scheme4_d", f"{DOC}:scaleup_3a_1mmol"],
         "unresolved": ["Scale-up paragraph omits full reagent quantities and operating conditions.", "SI is mentioned but not supplied/associated."]},
    ],
}
st.write(RUN / "coverage/pages_0007_0012.json", coverage)

st.write(BUNDLE_PATH, bundle)
print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]), "new_reactions": len(new_ids), "coverage_pages": [7, 8, 9, 10, 11, 12]}, ensure_ascii=False))
