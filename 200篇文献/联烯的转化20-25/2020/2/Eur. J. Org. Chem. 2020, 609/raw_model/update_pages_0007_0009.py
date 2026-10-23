#!/usr/bin/env python3
"""Incrementally enrich the v3.1 extraction from source PDF pages 7--9.

The pre-existing page-1--6 records are used only to preserve stable IDs and
avoid duplicate experiments.  Pages 7--9 are the newly inspected source
range.  Complex structures remain unresolved when an evidence-checked graph
cannot be encoded safely; every machine record remains needs_review.
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
DOC = "doc_ef264899de17791f"
E = {p: f"{DOC}:page{p}:image" for p in range(1, 10)}
T = {p: f"{DOC}:page{p}:text" for p in range(1, 10)}


def evidence(*pages: int) -> list[str]:
    out: list[str] = []
    for page in pages:
        out.extend([E[page], T[page]])
    return out


def f(raw=None, normalized=None, ev=(), method="original_extraction", *,
      derivation=None, status=None, verification=None, alternatives=()):
    return st.fact(raw, normalized, ev, method=method, derivation=derivation,
                   status=status, verification=verification,
                   alternatives=alternatives)


def missing(ev=(), status="unreported"):
    return f(ev=ev, status=status)


def q(number, unit, raw_text=None, **extra):
    out = {"number": number, "unit": unit}
    if raw_text is not None:
        out["raw_text"] = raw_text
    out.update(extra)
    return out


def inferred(raw, normalized, ev, why):
    return f(raw, normalized, ev, method="inference", derivation=why)


def structure_fact(raw, normalized, ev, why):
    return f(raw, normalized, ev, method="structure_derivation", derivation=why)


bundle = st.read(RUN / "standard_records.json")
# The prior formal bundle was prepared before page-image evidence for pages
# 7--9 was registered.  Keep the records, but take the immutable document and
# evidence registry from the current mechanical dataset so every cited image
# is a real, hash-bound source asset.
bundle["documents"] = copy.deepcopy(st.read(RUN / "dataset.json")["documents"])
molecules = {m["molecule_id"]: m for m in bundle["molecules"]}
reactions = {r["reaction_id"]: r for r in bundle["reactions"]}


def set_molecule_evidence(m, ev):
    occ = m["occurrences"][0]
    occ["evidence_ids"] = sorted(set(occ["evidence_ids"]) | set(ev))


def new_named_molecule(mid, name, ev, *, label=None, identity_note=None):
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label or name, label or name, ev)
    m["names"]["original"] = f(name, name, ev)
    m["names"]["english"] = missing(ev)
    m["names"]["abbreviations"] = missing(ev)
    m["names"]["aliases"] = missing(ev)
    m["cas_numbers"] = missing(ev)
    m["representation_type"] = "unresolved"
    m["smiles"] = missing(ev)
    m["generic_template"] = None
    for key in ("atom_centers", "double_bonds", "axial", "mixture", "assignment",
                "isomeric_smiles", "identity_key"):
        m["stereochemistry"][key] = missing(ev)
    why = "Source name or role was read on pages 7-8, but no complete evidence-checked stereochemical graph was encoded."
    m["stereochemistry"]["chirality_status"] = inferred(name, "unresolved", ev, why)
    m["stereochemistry"]["geometry_status"] = inferred(name, "unresolved", ev, why)
    m["stereochemistry"]["representation_status"] = inferred(name, "unresolved", ev, why)
    note = identity_note or "Source-named process component; structure not encoded."
    m["identity_note"] = f(note, note, ev)
    m["external_properties"] = []
    molecules[mid] = m
    return m


def ensure_component(mid, name, ev):
    if mid not in molecules:
        new_named_molecule(mid, name, ev, label=name,
                           identity_note="Named process component; no CAS or structure was inferred from the name.")
    else:
        set_molecule_evidence(molecules[mid], ev)
    return mid


# New named process components on pages 7-8.  The 35% peroxide is not merged
# with the pre-existing 30% peroxide solution because concentration is part of
# the reported material identity.
for mid, name, pages in [
    ("H2O2_35", "35% (w/w) aqueous hydrogen peroxide", (7,)),
    ("HCl_1N", "1 N hydrochloric acid", (7,)),
    ("NaBH4", "sodium borohydride", (7,)),
]:
    ensure_component(mid, name, evidence(*pages))


# Evidence-linked candidate structures for the two simple, fully named
# cyclohexanone/dione compounds.  RDKit diagnostics are graph checks only.
ev20 = evidence(3, 7)
m20 = molecules["m20"]
set_molecule_evidence(m20, ev20)
m20_name = "2-(but-3-yn-1-yl)-2-methylcyclohexane-1,3-dione (20)"
m20["names"]["original"] = f(m20_name, m20_name, ev20)
smi20 = "O=C1C(C)(CCC#C)C(=O)CCC1"
diag20 = st.stereo_inspect(smi20)
m20["representation_type"] = "smiles"
m20["smiles"] = structure_fact(None, smi20, ev20,
    "Candidate connectivity derived from the complete source name and Scheme 4 drawing; parsing does not establish source correctness.")
m20["stereochemistry"]["atom_centers"] = structure_fact(
    "standard stereo diagnostic: no potential tetrahedral stereocenters",
    [], ev20, "RDKit standard stereo diagnostic on the candidate graph found no potential tetrahedral stereocenters.")
m20["stereochemistry"]["double_bonds"] = structure_fact(
    "candidate graph contains carbonyl double bonds and a terminal alkyne, with no E/Z-applicable alkene",
    diag20["double_bonds"], ev20, "Computed from the candidate graph; carbonyl bonds are not E/Z stereoelements.")
m20["stereochemistry"]["axial"] = f("no allene/atropisomeric axis in candidate graph", [], ev20)
m20["stereochemistry"]["mixture"] = missing(ev20)
m20["stereochemistry"]["assignment"] = missing(ev20, "not_applicable")
m20["stereochemistry"]["chirality_status"] = structure_fact(
    "candidate graph has no potential stereocenters", "achiral", ev20,
    "Achiral status follows the source-derived candidate graph and standard stereo diagnostic; it remains a machine candidate.")
m20["stereochemistry"]["geometry_status"] = structure_fact(
    "no E/Z-applicable alkene in candidate graph", "not_applicable", ev20,
    "The candidate graph contains only carbonyl double bonds and a terminal alkyne.")
m20["stereochemistry"]["representation_status"] = structure_fact(
    "candidate graph contains no stereochemical elements requiring tags", "encoded", ev20,
    "No stereochemical tags are required for the candidate graph.")
m20["stereochemistry"]["isomeric_smiles"] = structure_fact(
    None, diag20["canonical_isomeric_smiles"], ev20,
    "Generated by standard stereo from the source-derived candidate graph.")
m20["stereochemistry"]["identity_key"] = structure_fact(
    None, diag20["stereo_graph_key"], ev20,
    "Generated by standard stereo; this is not a sample identity or source verification.")
m20["identity_note"] = f(
    "Complete source name and Scheme 4 connectivity used for a machine candidate; formula/HRMS characterization appears on page 7.",
    "Complete source name and Scheme 4 connectivity used for a machine candidate; formula/HRMS characterization appears on page 7.", ev20)


ev21 = evidence(3, 7)
m21 = molecules["m21"]
set_molecule_evidence(m21, ev21)
m21_name = "(rel-2S,3S)-2-(but-3-yn-1-yl)-3-[(tert-butyldimethylsilyl)oxy]-2-methylcyclohexan-1-one (21)"
m21["names"]["original"] = f(m21_name, m21_name, ev21)
smi21 = "O=C1C(C)(CCC#C)C(O[Si](C)(C)C(C)(C)C)CCC1"
diag21 = st.stereo_inspect(smi21)
m21["representation_type"] = "smiles"
m21["smiles"] = structure_fact(None, smi21, ev21,
    "Connectivity-only candidate derived from the complete source name and Scheme 4; relative stereochemistry and mixture composition are not encoded in this SMILES.")
m21["stereochemistry"]["atom_centers"] = f(
    "major product heading: rel-2S,3S; absolute configuration is not assigned",
    [{"source_positions": "2,3", "configuration": "relative rel-2S,3S", "absolute_assignment": "unassigned"}], ev21)
m21["stereochemistry"]["double_bonds"] = structure_fact(
    "no E/Z-applicable alkene in connectivity candidate", diag21["double_bonds"], ev21,
    "Computed from the connectivity-only candidate; the remaining double bond is a carbonyl.")
m21["stereochemistry"]["axial"] = f("no allene axis in compound 21", [], ev21)
m21["stereochemistry"]["mixture"] = f(
    "Scheme 4 reports dr 5:1",
    {"type": "diastereomeric_mixture", "ratio": "5:1", "major": "rel-2S,3S", "minor": "configuration unassigned"},
    evidence(3))
m21["stereochemistry"]["assignment"] = f(
    "relative stereochemistry of the siloxy group was determined by NOESY and comparison with prior literature; the major hydroxy/siloxy group is anti to the adjacent methyl and syn to the larger adjacent substituent",
    {"assignment_type": "relative", "major": "rel-2S,3S", "method": ["NOESY", "comparison with prior literature"],
     "relationship": "siloxy anti to adjacent methyl and syn to the larger adjacent substituent"}, ev21)
m21["stereochemistry"]["chirality_status"] = inferred(
    "rel-2S,3S major; dr 5:1; minor configuration unassigned", "mixture", ev21,
    "The source reports a diastereomeric ratio and a relative assignment for the major component, not an absolute single-isomer identity.")
m21["stereochemistry"]["geometry_status"] = structure_fact(
    "no E/Z-applicable alkene in candidate graph", "not_applicable", ev21,
    "The candidate connectivity contains only a carbonyl double bond and a terminal alkyne.")
m21["stereochemistry"]["representation_status"] = inferred(
    "connectivity encoded; relative rel-2S,3S and the 5:1 mixture are not encoded", "partial", ev21,
    "The SMILES intentionally omits absolute @/@@ tags because the source reports relative configuration and a mixture.")
m21["stereochemistry"]["isomeric_smiles"] = missing(ev21)
m21["stereochemistry"]["identity_key"] = missing(ev21)
m21["identity_note"] = f(
    "Compound 21 is represented as the reported 5:1 diastereomeric material; the major component is rel-2S,3S and the minor configuration is not assigned.",
    "Compound 21 is represented as the reported 5:1 diastereomeric material; the major component is rel-2S,3S and the minor configuration is not assigned.", ev21)


# Crude ketol is an explicitly used intermediate and the product class of the
# parallel NaBH4 experiment.  Its connectivity is derivable, but its relative
# mixture is not encoded in SMILES.
ev_ketol = evidence(7)
if "m_21_ketol" not in molecules:
    ketol = new_named_molecule("m_21_ketol", "crude ketol precursor to 21", ev_ketol,
                               label="crude ketol precursor")
else:
    ketol = molecules["m_21_ketol"]
smi_ketol = "O=C1C(C)(CCC#C)C(O)CCC1"
diag_ketol = st.stereo_inspect(smi_ketol)
ketol["representation_type"] = "smiles"
ketol["smiles"] = structure_fact(None, smi_ketol, ev_ketol,
    "Connectivity candidate derived from reduction of one carbonyl of 20 and the named structure of silylated product 21; stereochemistry is not encoded.")
ketol["stereochemistry"]["atom_centers"] = f(
    "LTBA reduction gives the relative configuration later named rel-2S,3S after silylation; NaBH4 gives the same stereochemical outcome at lower dr",
    [{"source_positions": "2,3", "configuration": "relative major as for 21", "absolute_assignment": "unassigned"}], ev_ketol)
ketol["stereochemistry"]["double_bonds"] = structure_fact(
    "candidate graph has no E/Z-applicable alkene", diag_ketol["double_bonds"], ev_ketol,
    "Computed from the connectivity candidate; the remaining double bond is a carbonyl.")
ketol["stereochemistry"]["axial"] = f("no allene axis in the ketol candidate", [], ev_ketol)
ketol["stereochemistry"]["mixture"] = f(
    "LTBA reduction: dr 5:1 in Scheme 4; NaBH4 parallel experiment: lower dr (number unreported)",
    {"type": "diastereomeric_mixture", "LTBA_ratio": "5:1", "NaBH4_ratio": "lower dr; numeric value unreported"},
    evidence(3, 7))
ketol["stereochemistry"]["assignment"] = f(
    "major hydroxy group anti to adjacent methyl and syn to the larger adjacent substituent",
    {"assignment_type": "relative", "relationship": "OH anti to methyl and syn to larger substituent"}, ev_ketol)
ketol["stereochemistry"]["chirality_status"] = inferred(
    "relative major configuration with reported diastereomeric mixtures", "mixture", ev_ketol,
    "Both reductions are described in terms of stereochemical outcome/dr; the minor configuration is not reported.")
ketol["stereochemistry"]["geometry_status"] = structure_fact(
    "no E/Z-applicable alkene", "not_applicable", ev_ketol,
    "The candidate contains a carbonyl and terminal alkyne only.")
ketol["stereochemistry"]["representation_status"] = inferred(
    "connectivity encoded; relative mixture not encoded", "partial", ev_ketol,
    "The source gives relative, mixture-level information that is not represented by the untagged connectivity SMILES.")
ketol["stereochemistry"]["isomeric_smiles"] = missing(ev_ketol)
ketol["stereochemistry"]["identity_key"] = missing(ev_ketol)
ketol["identity_note"] = f(
    "Crude, unisolated ketol intermediate; connectivity is a machine candidate and relative mixture composition is retained separately.",
    "Crude, unisolated ketol intermediate; connectivity is a machine candidate and relative mixture composition is retained separately.", ev_ketol)


# Improve source names and stereo states for compounds whose complete headings
# occur on pages 7-8.  Complex graphs remain intentionally unencoded.
name_updates = {
    "m22": ("tert-butyldimethyl({(rel-1S,2R)-2-methyl-2-[5-methyl-7-(2-methyl-1,3-dithian-2-yl)hepta-3,4-dien-1-yl]-3-methylenecyclohexyl}oxy)silane (22)", evidence(3, 7)),
    "m23": ("tert-butyl{[(rel-1S,8aR)-6-(2,5-dimethylcyclopent-1-en-1-yl)-8a-methyl-1,2,3,4,6,7,8,8a-octahydronaphthalen-1-yl]oxy}dimethylsilane (23)", evidence(3, 8)),
    "m24": ("tert-butyldimethyl{[(rel-4R,6aR,7S,10aS,10bS)-3,6a,10b-trimethyl-2,4,5,6,6a,7,8,9,10,10b-decahydro-1H-4,10a-methanobenzo[e]azulen-7-yl]oxy}silane (24)", evidence(3, 8)),
}
for mid, (name, ev) in name_updates.items():
    m = molecules[mid]
    set_molecule_evidence(m, ev)
    m["names"]["original"] = f(name, name, ev)
    m["representation_type"] = "unresolved"
    m["smiles"] = missing(ev)
    m["stereochemistry"]["chirality_status"] = inferred(
        name, "partially_specified", ev,
        "The heading gives relative configurations (rel), but no absolute sample assignment is stated; the complex graph is not encoded.")
    m["stereochemistry"]["representation_status"] = inferred(
        name, "not_encoded", ev,
        "Relative stereochemistry is preserved textually because an evidence-checked isomeric SMILES was not encoded.")
    m["stereochemistry"]["isomeric_smiles"] = missing(ev)
    m["stereochemistry"]["identity_key"] = missing(ev)

m22 = molecules["m22"]
m22["stereochemistry"]["atom_centers"] = f(
    "compound heading gives rel-1S,2R", [{"source_positions": "1,2", "configuration": "relative rel-1S,2R"}], evidence(7))
m22["stereochemistry"]["axial"] = f(
    "allene C=C=C is present; axial configuration is not assigned",
    [{"type": "allene", "configuration": "unspecified"}], evidence(3, 7))
m22["stereochemistry"]["geometry_status"] = inferred(
    "no E/Z assignment is reported; allene axial configuration is handled separately", "unresolved", evidence(3, 7),
    "The complex source drawing was not encoded, so E/Z applicability is not claimed.")
m22["stereochemistry"]["assignment"] = f(
    "rel-1S,2R in the compound heading; allene axial configuration not assigned",
    {"central_relative": "rel-1S,2R", "allene_axis": "unassigned"}, evidence(3, 7))

for mid, rel in [("m23", "rel-1S,8aR"), ("m24", "rel-4R,6aR,7S,10aS,10bS")]:
    m = molecules[mid]
    m["stereochemistry"]["atom_centers"] = f(
        f"compound heading gives {rel}", [{"configuration": rel, "assignment_type": "relative"}], evidence(8))
    m["stereochemistry"]["axial"] = missing(evidence(3, 8), "not_applicable")
    m["stereochemistry"]["assignment"] = f(
        f"relative configuration {rel}; no absolute configuration is stated",
        {"assignment_type": "relative", "configuration": rel}, evidence(3, 8))
    m["stereochemistry"]["geometry_status"] = inferred(
        "source contains unsaturated rings but no explicit E/Z assignment", "unresolved", evidence(3, 8),
        "Ring/double-bond geometry was not converted to E/Z without an encoded source-checked graph.")


def reaction(rid, pages, *, categories=("other",), name=None, transformation=None,
             location=None):
    ev = evidence(*pages)
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = "specific_experiment"
    r["experiment_categories"] = f(list(categories), list(categories), ev)
    r["location"] = f(location or {"pages": list(pages)}, location or {"pages": list(pages)}, ev)
    r["reaction_name"] = f(name, name, ev) if name else missing(ev)
    r["transformation_type"] = inferred(
        transformation, transformation, ev,
        "Transformation type inferred from the source procedure and Scheme 4 drawing.") if transformation else missing(ev)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    return r, ev


def participant(r, pid, mid, role, ev, *, stage="reaction", product_role=None,
                mass=None, amount=None, volume=None, equiv=None, mol_percent=None,
                reference=None):
    p = st.new_participant(pid, mid, role, ev, stage=stage)
    if product_role is not None:
        p["product_role"] = f(product_role, product_role, ev)
    for key, item in {"mass": mass, "amount": amount, "volume": volume,
                      "equivalents": equiv, "mol_percent": mol_percent}.items():
        if item is not None:
            raw, normalized = item
            p["charges"][key] = f(raw, normalized, ev)
    if reference is not None:
        p["equivalence_reference"] = f(reference, reference, ev)
    r["participants"].append(p)
    set_molecule_evidence(molecules[mid], ev)
    return pid


def step(r, sid, ev, operation, *, temp=None, time=None, pids=(), stage="reaction",
         atmosphere=None, atmosphere_ev=None, method="original_extraction"):
    s = st.new_step(sid, ev, stage=stage)
    s["participant_ids"] = list(pids)
    s["operation"] = f(operation, operation, ev, method=method,
                         derivation="Condition inherited from the immediately preceding procedure explicitly said to be repeated." if method == "inherited_general_condition" else None)
    if temp is not None:
        s["temperature"] = f(temp[0], temp[1], ev, method=method,
                              derivation="Condition inherited from the immediately preceding procedure explicitly said to be repeated." if method == "inherited_general_condition" else None)
    if time is not None:
        s["time"] = f(time[0], time[1], ev, method=method,
                       derivation="Condition inherited from the immediately preceding procedure explicitly said to be repeated." if method == "inherited_general_condition" else None)
    if atmosphere is not None:
        aev = atmosphere_ev or ev
        s["atmosphere"] = f(atmosphere, atmosphere, aev,
                              method="inherited_general_condition",
                              derivation="Inherited from the paper's General Information statement that all reactions were carried out under Ar.")
    r["steps"].append(s)
    return s


def outcome(r, oid, kind, targets, ev, raw, normalized, *, role=None,
            yield_type=None, basis=None, scope=None, basis_status=None):
    o = st.new_outcome(oid, kind, targets, ev)
    if role is not None:
        o["product_role"] = f(role, role, ev)
    o["value"] = f(raw, normalized, ev)
    if yield_type is not None:
        o["yield_type"] = f(yield_type, yield_type, ev)
    if basis is not None:
        o["measurement_basis"] = f(basis, basis, ev)
    elif basis_status:
        o["measurement_basis"] = missing(ev, basis_status)
    if scope is not None:
        o["scope"] = f(scope, scope, ev)
    r["outcomes"].append(o)
    return o


CAT = {row["code"]: row["number"] for row in st.read(
    ROOT / "skills/chemical-reaction-extraction/references/bond_taxonomy.json")["categories"]}


def classification(r, ev, codes, note):
    changes = []
    labels = {}
    for i, code in enumerate(codes, 1):
        if code.startswith("DB-"):
            elems, order, ctype = code[3:].split("-"), 2, "formed"
        elif code.startswith("TB-"):
            elems, order, ctype = code[3:].split("-"), 3, "formed"
        else:
            elems, order = code.split("-"), 1
            ctype = "hydrogen_gain" if "H" in elems else "formed"
        cid = f"b{i}"
        change = {"change_id": cid,
                  "atoms": [{"atom_map": None, "element": elems[0]},
                            {"atom_map": None, "element": elems[1]}],
                  "before_order": 0, "after_order": order, "change_type": ctype}
        if ctype == "hydrogen_gain":
            change["hydrogen_count"] = 1
        changes.append(structure_fact(None, change, ev, note))
        labels.setdefault(code, {"category_number": CAT[code], "code": code,
                                 "supporting_change_ids": []})["supporting_change_ids"].append(cid)
    r["classification"]["bond_changes"] = changes
    r["classification"]["labels"] = structure_fact(None, list(labels.values()), ev, note)
    r["classification"]["coverage"] = inferred(None, "partial", ev, note)
    r["classification"]["limitations"] = [
        "Partial element-level bond-change inference from the source procedure/drawing; atom mapping was not asserted.",
        "Implicit-H/proton source and omitted coproduct mappings are unresolved.",
    ]
    r["classification"]["model_votes"] = []


def set_solvent(r, ev, composition, *, ratio=None, volume=None, concentration=None,
                basis=None, method="original_extraction"):
    derivation = "Condition inherited from the immediately preceding procedure explicitly said to be repeated." if method == "inherited_general_condition" else None
    r["solvent_conditions"]["composition"] = f(composition, composition, ev, method=method, derivation=derivation)
    if ratio is not None:
        r["solvent_conditions"]["mixture_ratio"] = f(ratio, ratio, ev, method=method, derivation=derivation)
    if volume is not None:
        r["solvent_conditions"]["total_volume"] = f(volume[0], volume[1], ev, method=method, derivation=derivation)
    if concentration is not None:
        r["solvent_conditions"]["concentration"] = f(concentration[0], concentration[1], ev, method=method, derivation=derivation)
    if basis is not None:
        r["solvent_conditions"]["concentration_basis"] = f(basis, basis, ev, method=method, derivation=derivation)


def set_procedure(r, ev, raw_text, *, addition=None, pretreatment=None,
                  general_reference=None, modifications=None, method="original_extraction"):
    derivation = "Procedure inherited because the source explicitly states that the immediately preceding procedure was repeated." if method == "inherited_general_condition" else None
    r["procedure"]["original_text"] = f(raw_text, raw_text, ev, method=method, derivation=derivation)
    if addition:
        r["procedure"]["addition_order"] = f(addition, addition, ev, method=method, derivation=derivation)
    if pretreatment:
        r["procedure"]["pretreatment"] = f(pretreatment, pretreatment, ev)
    if general_reference:
        r["procedure"]["general_procedure_reference"] = f(general_reference, general_reference, ev)
    if modifications:
        r["procedure"]["modifications"] = f(modifications, modifications, ev)


AR_EV = evidence(5)


# r23: first independently reported epoxidation batch on page 7.
r, ev = reaction("r23_1_to_s4_epoxide_batch1", (3, 5, 7),
                 name="nucleophilic epoxidation (Scheme 4 batch 1)",
                 transformation="nucleophilic enone epoxidation",
                 location="Scheme 4 stage i and experimental section for compound 20, first epoxidation batch")
participant(r, "r23_p1", "m01", "reactant", ev,
            mass=("2.0 g", q(2.0, "g")), amount=("10.997 mmol", q(10.997, "mmol")),
            equiv=("1 equiv.", q(1, "equiv")), reference="m01")
participant(r, "r23_p2", "MeOH", "solvent", ev, volume=("55 mL", q(55, "mL")))
participant(r, "r23_p3", "H2O2_35", "reagent", ev,
            volume=("2.8 mL", q(2.8, "mL")), amount=("32.991 mmol", q(32.991, "mmol")),
            equiv=("3 equiv.", q(3, "equiv")), reference="m01")
participant(r, "r23_p4", "NaOH_10", "reagent", ev,
            volume=("2.0 mL", q(2.0, "mL")), amount=("5.498 mmol", q(5.498, "mmol")),
            equiv=("0.5 equiv.", q(0.5, "equiv")), reference="m01")
participant(r, "r23_p5", "m_s4_epoxide", "product", ev,
            product_role="nucleophilic epoxidation intermediate", mass=("1700 mg", q(1700, "mg")))
for pid, mid, role, vol, stage_name in [
    ("r23_w1", "Na2S2O3_aq", "workup_agent", ("200 mL", q(200, "mL")), "workup"),
    ("r23_w2", "EtOAc", "workup_agent", ("3 × 200 mL", q(200, "mL", raw_text="3 × 200 mL", multiplicity=3)), "workup"),
    ("r23_w3", "MgSO4", "workup_agent", None, "workup"),
    ("r23_w4", "silica", "purification_agent", None, "purification"),
    ("r23_w5", "hexane", "purification_agent", None, "purification"),
]:
    participant(r, pid, mid, role, ev, stage=stage_name, volume=vol)
r["reaction_scale"] = f(
    "Wieland–Miescher ketone: 10.997 mmol",
    {"reference_molecule_id": "m01", **q(10.997, "mmol")}, ev)
set_solvent(r, ev, "methanol", volume=("55 mL", q(55, "mL")),
            concentration=("0.2 M", q(0.2, "M")), basis="compound 1")
step(r, "r23_s1", ev, "Add 35% H2O2 solution and 10% NaOH solution dropwise to 1 in methanol.",
     temp=("0 °C", q(0, "°C")), pids=("r23_p1", "r23_p2", "r23_p3", "r23_p4"),
     atmosphere="Ar", atmosphere_ev=AR_EV)
step(r, "r23_s2", ev, "Stir the resulting mixture.", temp=("0 °C", q(0, "°C")),
     time=("12 h", q(12, "h")), pids=("r23_p1", "r23_p3", "r23_p4"),
     atmosphere="Ar", atmosphere_ev=AR_EV)
step(r, "r23_s3", ev, "Quench with saturated sodium thiosulfate (or remove residual H2O2 under reduced pressure), extract with ethyl acetate, dry, filter and concentrate.",
     pids=("r23_w1", "r23_w2", "r23_w3"), stage="workup")
step(r, "r23_s4", ev, "Flash chromatograph on silica gel (ethyl acetate/hexane 1:3 to 1:2).",
     pids=("r23_w4", "r23_w5"), stage="purification")
outcome(r, "r23_o1", "yield", ["m_s4_epoxide"], ev, "80 %", q(80, "%"),
        role="nucleophilic epoxidation intermediate", yield_type="isolated yield",
        basis="isolated after flash column chromatography", scope="compound 1 to the nucleophilic epoxidation intermediate; first batch")
set_procedure(r, ev,
    "Wieland–Miescher ketone in methanol was treated dropwise with 35% H2O2 and 10% NaOH at 0 °C, stirred 12 h, worked up, and chromatographed to give 1700 mg (80%) epoxide intermediate.",
    addition="H2O2 solution and NaOH solution were added dropwise at 0 °C; individual addition durations were not reported.")
classification(r, ev, ["C-O", "C-O"],
               "The named nucleophilic epoxidation forms the two C–O bonds of an epoxide; the unlabelled intermediate lacks a complete source structure/atom map.")
r["review"]["unresolved"].append("The epoxide intermediate is unlabelled and not structurally depicted as a standalone compound; its SMILES and atom mapping remain unresolved.")
reactions[r["reaction_id"]] = r


# r24: explicitly repeated epoxidation batch.  Only its scale, isolated mass,
# and yield are independently stated; reagent charges are not back-calculated.
r, ev = reaction("r24_1_to_s4_epoxide_batch2", (7,),
                 name="repeated nucleophilic epoxidation (Scheme 4 batch 2)",
                 transformation="nucleophilic enone epoxidation",
                 location="Experimental section for compound 20, repeated 2.9 g epoxidation batch")
participant(r, "r24_p1", "m01", "reactant", ev,
            mass=("2.9 gram scale", q(2.9, "g", raw_text="2.9 gram scale")), reference="m01")
participant(r, "r24_p2", "MeOH", "solvent", ev)
participant(r, "r24_p3", "H2O2_35", "reagent", ev)
participant(r, "r24_p4", "NaOH_10", "reagent", ev)
participant(r, "r24_p5", "m_s4_epoxide", "product", ev,
            product_role="nucleophilic epoxidation intermediate, repeated batch",
            mass=("2262 mg", q(2262, "mg")))
for pid, mid, role, stage_name in [
    ("r24_w1", "Na2S2O3_aq", "workup_agent", "workup"),
    ("r24_w2", "EtOAc", "workup_agent", "workup"),
    ("r24_w3", "MgSO4", "workup_agent", "workup"),
    ("r24_w4", "silica", "purification_agent", "purification"),
    ("r24_w5", "hexane", "purification_agent", "purification"),
]:
    participant(r, pid, mid, role, ev, stage=stage_name)
r["reaction_scale"] = f(
    "procedure repeated on 2.9 gram scale",
    {"reference_molecule_id": "m01", **q(2.9, "g", raw_text="2.9 gram scale")}, ev)
set_solvent(r, ev, "methanol (procedure repeated; volume/concentration for this batch unreported)",
            method="inherited_general_condition")
step(r, "r24_s1", ev, "Repeat the immediately preceding H2O2/NaOH epoxidation procedure on a 2.9 g scale.",
     temp=("0 °C", q(0, "°C")), pids=("r24_p1", "r24_p2", "r24_p3", "r24_p4"),
     atmosphere="Ar", atmosphere_ev=AR_EV, method="inherited_general_condition")
step(r, "r24_s2", ev, "Stir according to the repeated procedure.", temp=("0 °C", q(0, "°C")),
     time=("12 h", q(12, "h")), pids=("r24_p1", "r24_p3", "r24_p4"),
     atmosphere="Ar", atmosphere_ev=AR_EV, method="inherited_general_condition")
step(r, "r24_s3", ev, "Use the workup and flash chromatography of the immediately preceding procedure.",
     pids=("r24_w1", "r24_w2", "r24_w3", "r24_w4", "r24_w5"), stage="workup",
     method="inherited_general_condition")
outcome(r, "r24_o1", "yield", ["m_s4_epoxide"], ev, "72 %", q(72, "%"),
        role="nucleophilic epoxidation intermediate, repeated batch", yield_type="isolated yield",
        basis="2262 mg obtained after the explicitly repeated preceding procedure", scope="compound 1 to epoxide intermediate; repeated 2.9 g batch")
set_procedure(r, ev,
    "The source states: 'This procedure has been repeated in 2.9 gram scale and yielded 2262 mg (72%)'.",
    general_reference="immediately preceding detailed epoxidation procedure on page 7",
    modifications="starting-material scale 2.9 g; batch-specific reagent charges and solvent volume are not reported",
    method="inherited_general_condition")
classification(r, ev, ["C-O", "C-O"],
               "The explicitly repeated nucleophilic epoxidation forms the two C–O bonds of the same epoxide intermediate; no atom mapping is asserted.")
r["review"]["unresolved"].extend([
    "The source reports only the repeated 2.9 g scale, 2262 mg product and 72% yield; batch-specific mmol, reagent charges and solvent volume remain unreported and were not scaled by calculation.",
    "The epoxide intermediate is unlabelled and lacks a standalone complete source structure.",
])
reactions[r["reaction_id"]] = r


# Replace r10 with the actual pooled epoxide-fragmentation execution.  The two
# upstream epoxidation batches are now independent records r23/r24.
r, ev = reaction("r10_1_to_20_scheme4", (3, 5, 7),
                 name="Eschenmoser–Tanabe fragmentation of pooled epoxide batches",
                 transformation="epoxide fragmentation to gamma-ethynyl diketone",
                 location="Scheme 4 stage ii and experimental section for compound 20")
participant(r, "r10_p1", "m_s4_epoxide", "reactant", ev,
            mass=("3962 mg (1700 mg + 2262 mg)", q(3962, "mg", raw_text="3962 mg (1700 mg + 2262 mg)")),
            amount=("20.398 mmol", q(20.398, "mmol")), equiv=("1 equiv.", q(1, "equiv")),
            reference="m_s4_epoxide")
participant(r, "r10_p2", "EtOH", "solvent", ev, volume=("81.6 mL", q(81.6, "mL")))
participant(r, "r10_p3", "TsNHNH2", "reagent", ev,
            mass=("4001 mg", q(4001, "mg")), amount=("20.806 mmol", q(20.806, "mmol")),
            equiv=("1.02 equiv.", q(1.02, "equiv")), reference="m_s4_epoxide")
participant(r, "r10_p4", "AcOH", "reagent", ev,
            volume=("1.18 mL", q(1.18, "mL")), amount=("20.398 mmol", q(20.398, "mmol")),
            equiv=("1 equiv.", q(1, "equiv")), reference="m_s4_epoxide")
participant(r, "r10_p5", "m20", "product", ev, product_role="diketone 20",
            mass=("2189 mg", q(2189, "mg")))
for pid, mid, role, vol, stage_name in [
    ("r10_w1", "K2CO3_aq", "workup_agent", None, "workup"),
    ("r10_w2", "EtOAc", "workup_agent", ("3 × 200 mL", q(200, "mL", raw_text="3 × 200 mL", multiplicity=3)), "workup"),
    ("r10_w3", "MgSO4", "workup_agent", None, "workup"),
    ("r10_w4", "silica", "purification_agent", None, "purification"),
    ("r10_w5", "hexane", "purification_agent", None, "purification"),
]:
    participant(r, pid, mid, role, ev, stage=stage_name, volume=vol)
r["reaction_scale"] = f(
    "pooled epoxide precursor: 20.398 mmol",
    {"reference_molecule_id": "m_s4_epoxide", **q(20.398, "mmol")}, ev)
set_solvent(r, ev, "ethanol", volume=("81.6 mL", q(81.6, "mL")),
            concentration=("0.25 M", q(0.25, "M")), basis="pooled epoxide precursor")
step(r, "r10_s1", ev, "Add p-toluenesulfonyl hydrazide to the pooled epoxide precursor in ethanol and stir.",
     temp=("room temperature", "room temperature"), time=("3 h", q(3, "h")),
     pids=("r10_p1", "r10_p2", "r10_p3"), atmosphere="Ar", atmosphere_ev=AR_EV)
step(r, "r10_s2", ev, "Add acetic acid and continue stirring.",
     temp=("room temperature", "room temperature"), time=("12 h", q(12, "h")),
     pids=("r10_p4",), atmosphere="Ar", atmosphere_ev=AR_EV)
step(r, "r10_s3", ev, "Quench into saturated potassium carbonate, extract with ethyl acetate, dry, filter and concentrate.",
     pids=("r10_w1", "r10_w2", "r10_w3"), stage="workup")
step(r, "r10_s4", ev, "Flash chromatograph on silica gel (ethyl acetate/hexane 1:4).",
     pids=("r10_w4", "r10_w5"), stage="purification")
outcome(r, "r10_o1", "yield", ["m20"], ev, "60 %", q(60, "%"),
        role="diketone 20", yield_type="isolated yield",
        basis="2189 mg isolated after flash column chromatography",
        scope="pooled epoxide precursor (3962 mg, 20.398 mmol) to compound 20")
set_procedure(r, ev,
    "The pooled epoxide precursor from two batches was treated with TsNHNH2 for 3 h at room temperature, then AcOH for 12 h at room temperature, worked up and chromatographed to give 20 (2189 mg, 60%).")
classification(r, ev, ["TB-C-C"],
               "The Eschenmoser–Tanabe fragmentation produces the terminal C≡C bond in compound 20; C–C cleavage and coproduct mappings are unresolved.")
r["review"]["unresolved"].append(
    "The stable legacy reaction ID retains '1_to_20_scheme4', but the detailed execution represented here starts from pooled isolated epoxide material; upstream batches are separate records r23 and r24.")
reactions[r["reaction_id"]] = r


# Detailed 20 -> 21 reduction/protection sequence.
r, ev = reaction("r11_20_to_21_scheme4", (2, 3, 5, 7),
                 name="LTBA reduction and TBS protection",
                 transformation="diastereoselective ketone reduction followed by silylation",
                 location="Scheme 4 and experimental section for compound 21")
participant(r, "r11_p1", "m20", "reactant", ev,
            mass=("510 mg", q(510, "mg")), amount=("2.863 mmol", q(2.863, "mmol")),
            equiv=("1 equiv.", q(1, "equiv")), reference="m20")
participant(r, "r11_p2", "THF", "solvent", ev, volume=("19 mL", q(19, "mL")))
participant(r, "r11_p3", "LTBA", "reagent", ev,
            volume=("3.58 mL of 1 M solution in THF", q(3.58, "mL", raw_text="1 M solution in THF, 3.58 mL")),
            amount=("3.578 mmol", q(3.578, "mmol")), equiv=("1.25 equiv.", q(1.25, "equiv")), reference="m20")
participant(r, "r11_p4", "m_21_ketol", "product", ev, product_role="crude ketol intermediate")
participant(r, "r11_p5", "HCl_1N", "workup_agent", ev, stage="workup", volume=("19 mL", q(19, "mL")))
participant(r, "r11_p6", "EtOAc", "workup_agent", ev, stage="workup",
            volume=("3 × 150 mL", q(150, "mL", raw_text="3 × 150 mL", multiplicity=3)))
participant(r, "r11_p7", "MgSO4", "workup_agent", ev, stage="workup")
participant(r, "r11_p8", "m_21_ketol", "reactant", ev,
            product_role="crude ketol used without further purification", equiv=("1 equiv.", q(1, "equiv")), reference="m_21_ketol")
participant(r, "r11_p9", "DCM", "solvent", ev, volume=("14.3 mL", q(14.3, "mL")))
participant(r, "r11_p10", "lutidine", "reagent", ev,
            volume=("0.67 mL", q(0.67, "mL")), amount=("5.725 mmol", q(5.725, "mmol")),
            equiv=("2 equiv.", q(2, "equiv")), reference="m_21_ketol")
participant(r, "r11_p11", "TBSOTf", "reagent", ev,
            volume=("1.1 mL", q(1.1, "mL")), amount=("4.580 mmol", q(4.580, "mmol")),
            equiv=("1.6 equiv.", q(1.6, "equiv")), reference="m_21_ketol")
participant(r, "r11_p12", "m21", "product", ev, product_role="compound 21, 5:1 diastereomeric material",
            mass=("719 mg", q(719, "mg")))
participant(r, "r11_w1", "NH4Cl_aq", "workup_agent", ev, stage="workup", volume=("150 mL", q(150, "mL")))
participant(r, "r11_w2", "EtOAc", "workup_agent", ev, stage="workup",
            volume=("3 × 150 mL", q(150, "mL", raw_text="3 × 150 mL", multiplicity=3)))
participant(r, "r11_w3", "MgSO4", "workup_agent", ev, stage="workup")
participant(r, "r11_w4", "silica", "purification_agent", ev, stage="purification")
participant(r, "r11_w5", "hexane", "purification_agent", ev, stage="purification")
r["reaction_scale"] = f(
    "compound 20: 2.863 mmol",
    {"reference_molecule_id": "m20", **q(2.863, "mmol")}, ev)
r["solvent_conditions"]["composition"] = f(
    ["THF in reduction stage", "dichloromethane in silylation stage"],
    [{"stage": "reduction", "solvent": "THF"}, {"stage": "silylation", "solvent": "dichloromethane"}], ev)
r["solvent_conditions"]["total_volume"] = f(
    ["19 mL THF", "14.3 mL dichloromethane"],
    [{"stage": "reduction", **q(19, "mL")}, {"stage": "silylation", **q(14.3, "mL")}], ev)
r["solvent_conditions"]["concentration"] = f(
    ["0.15 M reduction", "0.2 M silylation"],
    [{"stage": "reduction", **q(0.15, "M")}, {"stage": "silylation", **q(0.2, "M")}], ev)
r["solvent_conditions"]["concentration_basis"] = f(
    ["compound 20", "crude ketol precursor"],
    [{"stage": "reduction", "basis": "compound 20"}, {"stage": "silylation", "basis": "crude ketol precursor"}], ev)
step(r, "r11_s1", ev, "Add LTBA solution in THF to 20 in THF.", temp=("−60 °C", q(-60, "°C")),
     pids=("r11_p1", "r11_p2", "r11_p3"), atmosphere="Ar", atmosphere_ev=AR_EV)
step(r, "r11_s2", ev, "Hold the reduction mixture.", temp=("−60 °C", q(-60, "°C")),
     time=("1 h", q(1, "h")), pids=("r11_p1", "r11_p3"), atmosphere="Ar", atmosphere_ev=AR_EV)
step(r, "r11_s3", ev, "Add 1 N HCl, extract with ethyl acetate, dry, filter and concentrate; use the crude ketol without further purification.",
     pids=("r11_p5", "r11_p6", "r11_p7"), stage="workup")
step(r, "r11_s4", ev, "Add 2,6-lutidine and TBSOTf to the crude ketol in dichloromethane.",
     temp=("0 °C", q(0, "°C")), pids=("r11_p8", "r11_p9", "r11_p10", "r11_p11"),
     atmosphere="Ar", atmosphere_ev=AR_EV)
step(r, "r11_s5", ev, "Stir after the additions.", temp=("0 °C", q(0, "°C")),
     time=("30 min", q(30, "min")), pids=("r11_p8", "r11_p10", "r11_p11"),
     atmosphere="Ar", atmosphere_ev=AR_EV)
step(r, "r11_s6", ev, "Quench with saturated NH4Cl, extract with ethyl acetate, dry, filter and concentrate.",
     pids=("r11_w1", "r11_w2", "r11_w3"), stage="workup")
step(r, "r11_s7", ev, "Flash chromatograph on silica gel (ethyl acetate/hexane 1:30).",
     pids=("r11_w4", "r11_w5"), stage="purification")
outcome(r, "r11_o1", "yield", ["m21"], ev, "719 mg, 85 % over two steps",
        q(85, "%", raw_text="719 mg, 85 % over two steps"), role="compound 21",
        yield_type="isolated yield over two steps", basis="719 mg isolated after flash column chromatography",
        scope="compound 20 to compound 21 over LTBA reduction and TBS protection")
o = outcome(r, "r11_o2", "selectivity", ["m21"], evidence(3), "dr 5:1", "5:1",
            role="diastereomeric material 21", scope="diastereomeric ratio, major rel-2S,3S : minor configuration unassigned",
            basis_status="unreported")
set_procedure(r, ev,
    "20 was reduced with LTBA for 1 h at −60 °C, the crude ketol was carried directly into TBS protection for 30 min at 0 °C, and chromatography afforded 21 in 85% over two steps (Scheme 4: dr 5:1).",
    addition="LTBA solution was added to 20 at −60 °C; 2,6-lutidine and TBSOTf were added to the crude ketol at 0 °C.")
classification(r, ev, ["C-H", "O-Si"],
               "Reduction gives a net C–H bond at the reduced carbonyl carbon and silylation forms O–Si; proton/hydride atom mapping is incomplete.")
r["review"]["unresolved"].append("The numerical dr 5:1 is printed in Scheme 4, but its measurement method is not reported on pages 7-9.")
reactions[r["reaction_id"]] = r


# Parallel NaBH4 stereochemical comparison, explicitly described as an
# experiment but without detailed conditions or a numerical dr.
r, ev = reaction("r25_parallel_nabh4_reduction", (7,),
                 categories=("mechanistic_study",),
                 name="parallel NaBH4 reduction stereochemical comparison",
                 transformation="ketone reduction",
                 location="Discussion following the experimental procedure for compound 21")
participant(r, "r25_p1", "m20", "reactant", ev)
participant(r, "r25_p2", "NaBH4", "reagent", ev)
participant(r, "r25_p3", "m_21_ketol", "product", ev,
            product_role="ketol with same major relative stereochemical outcome as LTBA reduction")
step(r, "r25_s1", ev,
     "Perform a parallel NaBH4 reduction experiment by reference to the cited literature; detailed conditions are not reported here.",
     pids=("r25_p1", "r25_p2"))
outcome(r, "r25_o1", "selectivity", ["m_21_ketol"], ev,
        "same stereochemical outcome with the LTBA case, but with a lower dr",
        "same major relative stereochemical outcome; dr lower than 5:1; numeric value unreported",
        role="crude ketol diastereomeric mixture", scope="qualitative diastereoselectivity comparison with the LTBA reduction",
        basis_status="unreported")
set_procedure(r, ev,
    "The authors state that a parallel NaBH4 reduction experiment, referring to prior literature, gave the same stereochemical outcome as LTBA but a lower dr.",
    general_reference="references 25a and 25b; no supporting-information file was supplied")
classification(r, ev, ["C-H"],
               "Net ketone reduction creates a C–H bond at the reduced carbonyl carbon; hydride/proton sources are not atom-mapped.")
r["review"]["unresolved"].extend([
    "Substrate/reagent charges, solvent, temperature, time, workup, isolation and yield are not reported in the page-7 statement.",
    "The NaBH4 experiment is described only qualitatively as giving a lower dr; neither the numerical ratio nor measurement method is reported.",
])
reactions[r["reaction_id"]] = r


# Detailed 21 -> 22 two-stage sequence.
r, ev = reaction("r12_21_to_22_scheme4", (2, 3, 5, 7, 8),
                 name="Wittig methylenation and Wang allenation",
                 transformation="carbonyl methylenation followed by copper-mediated allene formation",
                 location="Scheme 4 and experimental section for compound 22")
participant(r, "r12_p1", "PPh3MeBr", "reagent", ev,
            mass=("4347 mg", q(4347, "mg")), amount=("11.927 mmol", q(11.927, "mmol")),
            equiv=("3 equiv.", q(3, "equiv")), reference="m21")
participant(r, "r12_p2", "KOtBu", "reagent", ev,
            mass=("1340 mg", q(1340, "mg")), amount=("11.927 mmol", q(11.927, "mmol")),
            equiv=("3 equiv.", q(3, "equiv")), reference="m21")
participant(r, "r12_p3", "THF", "solvent", ev,
            volume=("6 mL + 4 mL substrate solution", q(10, "mL", raw_text="6 mL + 4 mL")))
participant(r, "r12_p4", "m21", "reactant", ev,
            mass=("1171 mg", q(1171, "mg")), amount=("3.976 mmol", q(3.976, "mmol")),
            equiv=("1 equiv.", q(1, "equiv")), reference="m21")
participant(r, "r12_p5", "m_21_methylene", "product", ev,
            product_role="isolated methylenated intermediate", mass=("1125 mg", q(1125, "mg")),
            amount=("3.845 mmol used in stage 2", q(3.845, "mmol")))
for pid, mid, role, vol, stage_name in [
    ("r12_w1", "NH4Cl_aq", "workup_agent", ("250 mL", q(250, "mL")), "workup"),
    ("r12_w2", "EtOAc", "workup_agent", ("3 × 250 mL", q(250, "mL", raw_text="3 × 250 mL", multiplicity=3)), "workup"),
    ("r12_w3", "MgSO4", "workup_agent", None, "workup"),
    ("r12_w4", "silica", "purification_agent", None, "purification"),
    ("r12_w5", "hexane", "purification_agent", None, "purification"),
]:
    participant(r, pid, mid, role, ev, stage=stage_name, volume=vol)
participant(r, "r12_p6", "m12", "reactant", ev,
            mass=("2043 mg", q(2043, "mg")), amount=("9.998 mmol", q(9.998, "mmol")),
            equiv=("2.6 equiv.", q(2.6, "equiv")), reference="m_21_methylene")
participant(r, "r12_p7", "TsNHNH2", "reagent", ev,
            mass=("1993 mg", q(1993, "mg")), amount=("10.383 mmol", q(10.383, "mmol")),
            equiv=("2.7 equiv.", q(2.7, "equiv")), reference="m_21_methylene")
participant(r, "r12_p8", "sieves_4A", "additive", ev)
participant(r, "r12_p9", "dioxane", "solvent", ev,
            volume=("10 mL then 38 mL carrier", q(48, "mL", raw_text="10 mL + 38 mL")))
participant(r, "r12_p10", "m_21_methylene", "reactant", ev,
            mass=("1125 mg", q(1125, "mg")), amount=("3.845 mmol", q(3.845, "mmol")),
            equiv=("1 equiv.", q(1, "equiv")), reference="m_21_methylene")
participant(r, "r12_p11", "CuI", "reagent", ev,
            mass=("754 mg", q(754, "mg")), amount=("3.960 mmol", q(3.960, "mmol")),
            equiv=("1 equiv.", q(1, "equiv")), reference="m_21_methylene")
participant(r, "r12_p12", "LiOtBu", "reagent", ev,
            mass=("1700 mg", q(1700, "mg")), amount=("21.148 mmol", q(21.148, "mmol")),
            equiv=("5.5 equiv.", q(5.5, "equiv")), reference="m_21_methylene")
participant(r, "r12_p13", "m22", "product", ev, product_role="allene 22",
            mass=("1380 mg", q(1380, "mg")))
participant(r, "r12_w6", "silica", "purification_agent", ev, stage="purification")
participant(r, "r12_w7", "EtOAc", "purification_agent", ev, stage="purification")
participant(r, "r12_w8", "hexane", "purification_agent", ev, stage="purification")
r["reaction_scale"] = f(
    ["compound 21: 3.976 mmol", "methylenated intermediate in stage 2: 3.845 mmol"],
    {"reference_molecule_id": "m21", **q(3.976, "mmol"),
     "stage_scales": [
         {"stage": "methylenation", "reference_molecule_id": "m21", **q(3.976, "mmol")},
         {"stage": "allenation", "reference_molecule_id": "m_21_methylene", **q(3.845, "mmol")},
     ]}, ev)
r["solvent_conditions"]["composition"] = f(
    ["THF in methylenation", "1,4-dioxane in hydrazone generation/allenation"],
    [{"stage": "methylenation", "solvent": "THF"}, {"stage": "allenation", "solvent": "1,4-dioxane"}], ev)
r["solvent_conditions"]["total_volume"] = f(
    ["6 mL THF plus 4 mL substrate solution", "10 mL 1,4-dioxane plus 38 mL substrate carrier"],
    [{"stage": "methylenation", **q(10, "mL", raw_text="6 mL + 4 mL")},
     {"stage": "allenation", **q(48, "mL", raw_text="10 mL + 38 mL")}], ev)
r["solvent_conditions"]["concentration"] = f(
    ["final concentration 0.4 M in methylenation", "12 in 1,4-dioxane: 1 M before coupling", "final concentration 0.08 M in coupling"],
    [{"stage": "methylenation", **q(0.4, "M")}, {"stage": "hydrazone generation", **q(1, "M")},
     {"stage": "coupling", **q(0.08, "M")}], ev)
r["solvent_conditions"]["concentration_basis"] = f(
    ["compound 21", "compound 12", "methylenated intermediate"],
    [{"stage": "methylenation", "basis": "compound 21"}, {"stage": "hydrazone generation", "basis": "compound 12"},
     {"stage": "coupling", "basis": "methylenated intermediate"}], ev)
step(r, "r12_s0", ev, "Heat methyltriphenylphosphonium bromide before reaction setup.",
     temp=("110 °C", q(110, "°C")), time=("overnight", "overnight"), pids=("r12_p1",),
     atmosphere="Ar", atmosphere_ev=AR_EV)
step(r, "r12_s1", ev, "Add KOtBu to the dry phosphonium salt in THF and stir.",
     temp=("0 °C", q(0, "°C")), time=("30 min", q(30, "min")),
     pids=("r12_p1", "r12_p2", "r12_p3"), atmosphere="Ar", atmosphere_ev=AR_EV)
step(r, "r12_s2", ev, "Add 21 in THF by cannula at 0 °C, then stir after warming to room temperature.",
     temp=("room temperature after addition at 0 °C", "room temperature"), time=("24 h", q(24, "h")),
     pids=("r12_p4",), atmosphere="Ar", atmosphere_ev=AR_EV)
step(r, "r12_s3", ev, "Quench with saturated NH4Cl, extract with ethyl acetate, dry, filter and concentrate.",
     pids=("r12_w1", "r12_w2", "r12_w3"), stage="workup")
step(r, "r12_s4", ev, "Flash chromatograph on silica gel with hexane to isolate the methylenated intermediate.",
     pids=("r12_w4", "r12_w5"), stage="purification")
step(r, "r12_s5", ev, "Combine 12, p-toluenesulfonyl hydrazide and 4 Å molecular sieves in 1,4-dioxane and reflux.",
     temp=("reflux", "reflux"), time=("12 h", q(12, "h")),
     pids=("r12_p6", "r12_p7", "r12_p8", "r12_p9"), atmosphere="Ar", atmosphere_ev=AR_EV)
step(r, "r12_s6", ev, "Cool to room temperature, then add the methylenated intermediate in 1,4-dioxane by cannula, CuI, additional 4 Å sieves and LiOtBu.",
     temp=("room temperature", "room temperature"), pids=("r12_p10", "r12_p11", "r12_p12", "r12_p8"),
     atmosphere="Ar", atmosphere_ev=AR_EV)
step(r, "r12_s7", ev, "Reflux the coupling mixture.", temp=("reflux", "reflux"),
     time=("24 h", q(24, "h")), pids=("r12_p10", "r12_p11", "r12_p12"),
     atmosphere="Ar", atmosphere_ev=AR_EV)
step(r, "r12_s8", ev, "Cool, filter through a silica pad, and flash chromatograph (ethyl acetate/hexane 1:50 to 1:30).",
     pids=("r12_w6", "r12_w7", "r12_w8"), stage="purification")
outcome(r, "r12_o1", "yield", ["m_21_methylene"], ev, "1125 mg, 97 %",
        q(97, "%", raw_text="1125 mg, 97 %"), role="isolated methylenated intermediate",
        yield_type="isolated yield", basis="1125 mg isolated after flash chromatography",
        scope="Wittig methylenation of 21")
outcome(r, "r12_o2", "yield", ["m22"], ev, "1380 mg, 75 %",
        q(75, "%", raw_text="1380 mg, 75 %"), role="allene 22",
        yield_type="isolated yield", basis="1380 mg isolated after flash chromatography",
        scope="coupling of isolated methylenated intermediate with in-situ reagent generated from 12")
set_procedure(r, ev,
    "Dry phosphonium salt and KOtBu were used for Wittig methylenation of 21 (97% isolated intermediate); a hydrazone generated from 12 was then coupled using CuI/LiOtBu in refluxing 1,4-dioxane to give allene 22 (75%).",
    pretreatment="Methyltriphenylphosphonium bromide was heated at 110 °C overnight.")
classification(r, ev, ["DB-C-C", "DB-C-C"],
               "The Wittig stage and Wang allenation stage each form a new C=C connection; full atom mapping and allene axial assignment remain unresolved.")
r["review"]["unresolved"].append("Allene 22 has an unassigned allene axis; no evidence-checked full isomeric SMILES was encoded.")
reactions[r["reaction_id"]] = r


# Detailed 22 -> 23 + 24 three-stage sequence.
r, ev = reaction("r13_22_to_23_24_scheme4", (2, 3, 5, 8),
                 name="modified TMM diyl-mediated tandem cycloaddition",
                 transformation="dithiane deprotection, hydrazone formation and TMM cyclization with competing 1,4-HAT",
                 location="Scheme 4 and experimental section for compounds 23 and 24")
participant(r, "r13_p1", "m22", "reactant", ev,
            mass=("1380 mg", q(1380, "mg")), amount=("2.870 mmol", q(2.870, "mmol")),
            equiv=("1 equiv.", q(1, "equiv")), reference="m22")
participant(r, "r13_p2", "MeCN", "solvent", ev)
participant(r, "r13_p3", "water", "solvent", ev)
participant(r, "r13_p4", "CaCO3", "reagent", ev,
            mass=("1160 mg", q(1160, "mg")), amount=("11.479 mmol", q(11.479, "mmol")),
            equiv=("4 equiv.", q(4, "equiv")), reference="m22")
participant(r, "r13_p5", "MeI", "reagent", ev,
            volume=("1.8 mL", q(1.8, "mL")), amount=("28.697 mmol", q(28.697, "mmol")),
            equiv=("10 equiv.", q(10, "equiv")), reference="m22")
participant(r, "r13_p6", "m_22_ketone", "product", ev,
            product_role="isolated allenyl ketone intermediate", mass=("1031 mg", q(1031, "mg")))
for pid, mid, role, vol, stage_name in [
    ("r13_w1", "NH4Cl_aq", "workup_agent", ("200 mL", q(200, "mL")), "workup"),
    ("r13_w2", "EtOAc", "workup_agent", ("3 × 200 mL", q(200, "mL", raw_text="3 × 200 mL", multiplicity=3)), "workup"),
    ("r13_w3", "MgSO4", "workup_agent", None, "workup"),
    ("r13_w4", "silica", "purification_agent", None, "purification"),
    ("r13_w5", "hexane", "purification_agent", None, "purification"),
]:
    participant(r, pid, mid, role, ev, stage=stage_name, volume=vol)
participant(r, "r13_p7", "m_22_ketone", "reactant", ev,
            mass=("437 mg", q(437, "mg")), amount=("1.117 mmol", q(1.117, "mmol")),
            equiv=("1 equiv.", q(1, "equiv")), reference="m_22_ketone")
participant(r, "r13_p8", "MeOH", "solvent", ev, volume=("3.72 mL", q(3.72, "mL")))
participant(r, "r13_p9", "TsNHNH2", "reagent", ev,
            mass=("257 mg", q(257, "mg")), amount=("1.341 mmol", q(1.341, "mmol")),
            equiv=("1.2 equiv.", q(1.2, "equiv")), reference="m_22_ketone")
participant(r, "r13_p10", "m_22_hydrazone", "product", ev,
            product_role="isolated hydrazone precursor", mass=("409 mg", q(409, "mg")))
participant(r, "r13_w6", "silica", "purification_agent", ev, stage="purification")
participant(r, "r13_w7", "EtOAc", "purification_agent", ev, stage="purification")
participant(r, "r13_w8", "hexane", "purification_agent", ev, stage="purification")
participant(r, "r13_p11", "NaH", "reagent", ev,
            mass=("36.64 mg of 60% (w/w) dispersion in mineral oil", q(36.64, "mg", raw_text="60% (w/w) dispersion in mineral oil, 36.64 mg")),
            amount=("0.916 mmol", q(0.916, "mmol")), equiv=("1.25 equiv.", q(1.25, "equiv")), reference="m_22_hydrazone")
participant(r, "r13_p12", "m_22_hydrazone", "reactant", ev,
            mass=("409.2 mg", q(409.2, "mg")), amount=("0.733 mmol", q(0.733, "mmol")),
            equiv=("1 equiv.", q(1, "equiv")), reference="m_22_hydrazone")
participant(r, "r13_p13", "toluene", "solvent", ev,
            volume=("110 mL + 12 mL substrate solution", q(122, "mL", raw_text="110 mL + 12 mL")))
participant(r, "r13_p14", "m23", "product", ev, product_role="1,4-hydrogen atom transfer product",
            mass=("43 mg", q(43, "mg")))
participant(r, "r13_p15", "m24", "product", ev, product_role="bridged cycloaddition product",
            mass=("165 mg", q(165, "mg")))
participant(r, "r13_w9", "NH4Cl_aq", "workup_agent", ev, stage="workup", volume=("200 mL", q(200, "mL")))
participant(r, "r13_w10", "EtOAc", "workup_agent", ev, stage="workup",
            volume=("3 × 200 mL", q(200, "mL", raw_text="3 × 200 mL", multiplicity=3)))
participant(r, "r13_w11", "MgSO4", "workup_agent", ev, stage="workup")
participant(r, "r13_w12", "silica", "purification_agent", ev, stage="purification")
r["reaction_scale"] = f(
    ["22: 2.870 mmol", "allenyl ketone portion: 1.117 mmol", "hydrazone precursor: 0.733 mmol"],
    {"reference_molecule_id": "m22", **q(2.870, "mmol"),
     "stage_scales": [
         {"stage": "dithiane deprotection", "reference_molecule_id": "m22", **q(2.870, "mmol")},
         {"stage": "hydrazone formation", "reference_molecule_id": "m_22_ketone", **q(1.117, "mmol")},
         {"stage": "cyclization", "reference_molecule_id": "m_22_hydrazone", **q(0.733, "mmol")},
     ]}, ev)
r["solvent_conditions"]["composition"] = f(
    ["MeCN/H2O in deprotection", "methanol in hydrazone formation", "toluene in cyclization"],
    [{"stage": "deprotection", "solvent": "MeCN/H2O"}, {"stage": "hydrazone formation", "solvent": "methanol"},
     {"stage": "cyclization", "solvent": "toluene"}], ev)
r["solvent_conditions"]["mixture_ratio"] = f("MeCN/H2O = 3:1 (v/v)", "3:1 (v/v)", ev)
r["solvent_conditions"]["total_volume"] = f(
    ["28.7 mL MeCN/H2O", "3.72 mL methanol", "110 mL toluene plus 12 mL substrate solution"],
    [{"stage": "deprotection", **q(28.7, "mL")}, {"stage": "hydrazone formation", **q(3.72, "mL")},
     {"stage": "cyclization", **q(122, "mL", raw_text="110 mL + 12 mL")}], ev)
r["solvent_conditions"]["concentration"] = f(
    ["0.1 M deprotection", "0.3 M hydrazone formation", "final concentration 0.006 M cyclization"],
    [{"stage": "deprotection", **q(0.1, "M")}, {"stage": "hydrazone formation", **q(0.3, "M")},
     {"stage": "cyclization", **q(0.006, "M")}], ev)
r["solvent_conditions"]["concentration_basis"] = f(
    ["compound 22", "allenyl ketone portion", "hydrazone precursor"],
    [{"stage": "deprotection", "basis": "compound 22"}, {"stage": "hydrazone formation", "basis": "allenyl ketone"},
     {"stage": "cyclization", "basis": "hydrazone precursor"}], ev)
step(r, "r13_s1", ev, "Add CaCO3 and iodomethane to 22 in MeCN/H2O 3:1 at room temperature.",
     temp=("room temperature", "room temperature"), pids=("r13_p1", "r13_p2", "r13_p3", "r13_p4", "r13_p5"),
     atmosphere="Ar", atmosphere_ev=AR_EV)
step(r, "r13_s2", ev, "Heat and stir the deprotection mixture.", temp=("80 °C", q(80, "°C")),
     time=("4 h", q(4, "h")), pids=("r13_p1", "r13_p4", "r13_p5"),
     atmosphere="Ar", atmosphere_ev=AR_EV)
step(r, "r13_s3", ev, "Cool, quench with saturated NH4Cl, extract with ethyl acetate, dry, filter and concentrate.",
     pids=("r13_w1", "r13_w2", "r13_w3"), stage="workup")
step(r, "r13_s4", ev, "Flash chromatograph (silica gel, hexane to ethyl acetate/hexane 1:20) to isolate the allenyl ketone.",
     pids=("r13_w4", "r13_w5"), stage="purification")
step(r, "r13_s5", ev, "Add p-toluenesulfonyl hydrazide to the allenyl ketone portion in methanol and stir.",
     temp=("room temperature", "room temperature"), time=("1.5 h", q(1.5, "h")),
     pids=("r13_p7", "r13_p8", "r13_p9"), atmosphere="Ar", atmosphere_ev=AR_EV)
step(r, "r13_s6", ev, "Concentrate and flash chromatograph (silica gel, ethyl acetate/hexane 1:10 to 1:5) to isolate the hydrazone precursor.",
     pids=("r13_w6", "r13_w7", "r13_w8"), stage="purification")
step(r, "r13_s7", ev, "Add freshly prepared hydrazone in toluene by cannula to NaH suspension at 0 °C.",
     temp=("0 °C", q(0, "°C")), pids=("r13_p11", "r13_p12", "r13_p13"),
     atmosphere="Ar", atmosphere_ev=AR_EV)
step(r, "r13_s8", ev, "Degas the resulting solution using an Ar balloon.",
     temp=("0 °C", q(0, "°C")), time=("25 min", q(25, "min")),
     pids=("r13_p11", "r13_p12", "r13_p13"), atmosphere="Ar", atmosphere_ev=AR_EV)
step(r, "r13_s9", ev, "Reflux the reaction mixture.", temp=("reflux", "reflux"),
     time=("24 h", q(24, "h")), pids=("r13_p11", "r13_p12", "r13_p13"),
     atmosphere="Ar", atmosphere_ev=AR_EV)
step(r, "r13_s10", ev, "Cool, quench with saturated NH4Cl, extract with ethyl acetate, dry, filter and concentrate.",
     pids=("r13_w9", "r13_w10", "r13_w11"), stage="workup")
step(r, "r13_s11", ev, "Flash chromatograph on silica gel to separate products 23 and 24.",
     pids=("r13_w12",), stage="purification")
outcome(r, "r13_o1", "yield", ["m_22_ketone"], ev, "1031 mg, 92 %",
        q(92, "%", raw_text="1031 mg, 92 %"), role="isolated allenyl ketone intermediate",
        yield_type="isolated yield", basis="1031 mg isolated after flash chromatography",
        scope="22 to allenyl ketone intermediate")
outcome(r, "r13_o2", "yield", ["m_22_hydrazone"], ev, "409 mg, 66 %",
        q(66, "%", raw_text="409 mg, 66 %"), role="isolated hydrazone precursor",
        yield_type="isolated yield", basis="409 mg isolated after flash chromatography",
        scope="437 mg allenyl ketone portion to hydrazone precursor")
outcome(r, "r13_o3", "yield", ["m23"], ev, "43 mg, 16 %",
        q(16, "%", raw_text="43 mg, 16 %"), role="1,4-hydrogen atom transfer product 23",
        yield_type="isolated yield", basis="43 mg isolated after flash chromatography",
        scope="hydrazone precursor to product 23")
outcome(r, "r13_o4", "yield", ["m24"], ev, "165 mg, 60 %",
        q(60, "%", raw_text="165 mg, 60 %"), role="bridged cycloaddition product 24",
        yield_type="isolated yield", basis="165 mg isolated after flash chromatography",
        scope="hydrazone precursor to product 24")
outcome(r, "r13_o5", "yield", ["m23", "m24"], evidence(3), "76 %", q(76, "%"),
        role="combined products 23 and 24", yield_type="combined isolated yield",
        basis="Scheme 4 reports the combined 76%; page 8 separately reports 16% and 60% after flash chromatography",
        scope="hydrazone cyclization to combined products 23 and 24")
outcome(r, "r13_o6", "selectivity", ["m23", "m24"], evidence(3),
        "23:24 = 4:15", "4:15", role="product ratio",
        scope="product ratio in the explicit order 23:24", basis_status="unreported")
set_procedure(r, ev,
    "22 was deprotected with MeI/CaCO3 (92% allenyl ketone), a portion was converted to isolated hydrazone (66%), and NaH/toluene cyclization furnished isolated 23 (16%) and 24 (60%); Scheme 4 reports 23:24=4:15 and 76% combined.")
classification(r, ev, ["C-C", "C-C", "C-H"],
               "The cycloaddition branch forms two C–C ring connections, while the competing 1,4-HAT branch has a net C–H formation; classification is the union of observed branches without product-specific atom maps.")
r["review"]["unresolved"].extend([
    "The 23:24 ratio is printed in Scheme 4, but its measurement method is not reported.",
    "Bond-change classification combines the cycloaddition and 1,4-HAT branches because complete product-specific atom maps were not asserted.",
])
reactions[r["reaction_id"]] = r


# Ensure all modified/new molecules have page evidence from their uses.
for r in reactions.values():
    for p in r["participants"]:
        if p["molecule_id"] in molecules:
            set_molecule_evidence(molecules[p["molecule_id"]], r["evidence_ids"])


# Non-synthetic source content seen on the requested pages.
new_context = [
    {"context_id": f"{DOC}:ctx_page7_characterization", "document_id": DOC,
     "kind": "characterization",
     "evidence_ids": evidence(7),
     "description": "Page 7 begins with continuation of 18/19 NMR/MS data and includes NMR/HRMS blocks for 20 and 21 plus the beginning of 22 characterization; these peak tables are not separate syntheses."},
    {"context_id": f"{DOC}:ctx_page8_characterization", "document_id": DOC,
     "kind": "characterization",
     "evidence_ids": evidence(8),
     "description": "Page 8 completes 22 HRMS and gives NMR/HRMS data for isolated products 23 and 24; characterization is linked context, not independent synthetic records."},
    {"context_id": f"{DOC}:ctx_pages8_9_backmatter", "document_id": DOC,
     "kind": "references_and_backmatter",
     "evidence_ids": evidence(8, 9),
     "description": "The right side of page 8 and all of page 9 contain acknowledgments, keywords, references and the received date. Cited papers are background and were not converted into current-paper experiments."},
]
existing_context_ids = {c.get("context_id") for c in bundle.get("source_context", [])}
for context in new_context:
    if context["context_id"] not in existing_context_ids:
        bundle["source_context"].append(context)


bundle["molecules"] = list(molecules.values())
bundle["reactions"] = list(reactions.values())


# Raw extraction inventory: source observations and decisions, not a second
# hand-authored formal record format.
raw_inventory = {
    "batch_id": "allene-stereo-v31-200-20261017",
    "document_id": DOC,
    "title": bundle["documents"][0]["title"],
    "pages_attempted": [7, 8, 9],
    "visual_inspection": {
        "page7": "Viewed full rendered page. Characterization continuation for 18/19; detailed preparations of 20 and 21; beginning and full operational details of 22 preparation.",
        "page8": "Viewed full rendered page. Completion of 22 characterization; detailed three-stage preparation of 23 and 24; acknowledgments and references.",
        "page9": "Viewed full rendered page. References 21-25 and received date only; no reaction arrow, table, or experimental procedure.",
    },
    "new_reaction_ids": [
        f"{DOC}:r23_1_to_s4_epoxide_batch1",
        f"{DOC}:r24_1_to_s4_epoxide_batch2",
        f"{DOC}:r25_parallel_nabh4_reduction",
    ],
    "enriched_reaction_ids": [
        f"{DOC}:r10_1_to_20_scheme4",
        f"{DOC}:r11_20_to_21_scheme4",
        f"{DOC}:r12_21_to_22_scheme4",
        f"{DOC}:r13_22_to_23_24_scheme4",
    ],
    "deduplicated_existing_context": [
        f"{DOC}:r09_14_to_18_19",
        "Page 7 characterization continuation for products 18/19 was not treated as a new synthesis.",
    ],
    "standard_stereo_checks": {"m20": diag20, "m21_connectivity": diag21, "m_21_ketol_connectivity": diag_ketol},
    "limitations": [
        "No supporting-information file was supplied or associated.",
        "Complex structures 22-24 retain source names and relative/allene-axis descriptions without speculative SMILES.",
        "The repeated 2.9 g epoxidation batch does not report batch-specific reagent charges; no scaled values were calculated.",
        "The parallel NaBH4 experiment reports only qualitative lower dr and no detailed conditions.",
        "All machine records remain needs_review; validation is not chemical verification.",
    ],
}
(RUN / "raw_model/pages_0007_0009_extraction.json").write_text(
    json.dumps(raw_inventory, ensure_ascii=False, indent=2) + "\n")
(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")


coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 7, "state": "attempted", "regions": [
            {"region": "upper-left continuation", "finding": "Continuation of 18/19 NMR and MS characterization; linked as source_context and not counted as a new synthesis.", "evidence_ids": evidence(7)},
            {"region": "left column, procedure for 20", "finding": "Two independent nucleophilic epoxidation batches (2.0 g, 80%; repeated 2.9 g, 72%) followed by a pooled Eschenmoser–Tanabe fragmentation to 20 (60%); represented by r23, r24 and enriched r10.", "evidence_ids": evidence(7)},
            {"region": "lower-left to upper-right, procedure for 21", "finding": "Detailed LTBA reduction/TBS protection sequence with −60 °C hold, separate workups, 85% over two steps, relative assignment discussion, and a separate qualitative NaBH4 comparison; represented by enriched r11 and new r25.", "evidence_ids": evidence(7)},
            {"region": "right column, procedure for 22", "finding": "Detailed Wittig methylenation (97%) and CuI/LiOtBu Wang allenation (75%) including pretreatment, addition order, times, concentrations, workup and purification; enriches r12.", "evidence_ids": evidence(7)},
            {"region": "NMR/HRMS blocks for 20-22", "finding": "Characterization retained in source_context; no synthetic records created from spectral peaks.", "evidence_ids": evidence(7)},
        ], "record_ids": [
            f"{DOC}:r09_14_to_18_19", f"{DOC}:r10_1_to_20_scheme4", f"{DOC}:r11_20_to_21_scheme4",
            f"{DOC}:r12_21_to_22_scheme4", f"{DOC}:r23_1_to_s4_epoxide_batch1",
            f"{DOC}:r24_1_to_s4_epoxide_batch2", f"{DOC}:r25_parallel_nabh4_reduction",
        ], "unresolved": [
            "The repeated 2.9 g epoxidation batch does not provide batch-specific reagent amounts or solvent volume.",
            "The NaBH4 comparison gives no numerical dr, measurement method, detailed conditions or yield.",
            "The epoxide intermediate is unlabelled and has no standalone complete structure drawing.",
            "Allene 22 has no axial assignment and no evidence-checked full SMILES was encoded.",
        ]},
        {"page": 8, "state": "attempted", "regions": [
            {"region": "upper-left continuation", "finding": "Completion of 22 characterization followed by the detailed three-stage 22→allenyl ketone→hydrazone→23+24 procedure; enriches r13 and provides product-specific isolated yields 16% and 60%.", "evidence_ids": evidence(8)},
            {"region": "middle-left NMR/HRMS blocks", "finding": "Characterization of 23 and 24 retained in source_context; spectral peaks are not separate reactions.", "evidence_ids": evidence(8)},
            {"region": "lower-left and right column", "finding": "Acknowledgments, keywords and references 1-21; cited works are background, not current-paper experiments.", "evidence_ids": evidence(8)},
        ], "record_ids": [f"{DOC}:r12_21_to_22_scheme4", f"{DOC}:r13_22_to_23_24_scheme4"],
         "unresolved": [
             "The Scheme 4 ratio 23:24=4:15 has no reported measurement method.",
             "Complex relative stereochemistry for 23/24 is preserved from headings but not encoded as isomeric SMILES.",
         ]},
        {"page": 9, "state": "attempted", "regions": [
            {"region": "entire page", "finding": "References 21-25 and received date only; no scheme, table, reaction arrow, experimental procedure, control, or new characterization record.", "evidence_ids": evidence(9)},
        ], "record_ids": [],
         "unresolved": ["No synthetic experiment is present on this page; cited references were not reclassified as executions in the current paper."]},
    ],
}
(RUN / "coverage/pages_0007_0009.json").write_text(
    json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")

print(json.dumps({
    "molecules": len(bundle["molecules"]),
    "reactions": len(bundle["reactions"]),
    "new_reactions": 3,
    "enriched_reactions": 4,
    "pages": [7, 8, 9],
}, ensure_ascii=False))
