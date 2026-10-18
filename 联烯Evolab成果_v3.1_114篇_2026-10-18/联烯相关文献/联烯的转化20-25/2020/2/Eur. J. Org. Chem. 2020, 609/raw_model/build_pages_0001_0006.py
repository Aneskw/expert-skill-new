#!/usr/bin/env python3
"""Build the source-derived v3.1 submission for original PDF pages 1--6.

This is a session extraction artifact, not a reusable chemistry rule.  Complex
drawn structures remain unresolved instead of being converted to speculative
SMILES.  All machine records remain needs_review after standard merge.
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DOC = "doc_ef264899de17791f"
E = {p: f"{DOC}:page{p}:image" for p in range(1, 7)}
T = {p: f"{DOC}:page{p}:text" for p in range(1, 7)}


def f(raw=None, normalized=None, ev=(), method="original_extraction", *, derivation=None,
      status=None):
    return st.fact(raw, normalized, ev, method=method, derivation=derivation, status=status)


def missing(ev=(), status="unreported"):
    return f(ev=ev, status=status)


def q(number, unit, raw_text=None, **extra):
    out = {"number": number, "unit": unit}
    if raw_text is not None:
        out["raw_text"] = raw_text
    out.update(extra)
    return out


def derived(raw, normalized, ev, why):
    return f(raw, normalized, ev, method="inference", derivation=why)


mol_evidence = defaultdict(set)
molecules = {}


def molecule(mid, label, pages, *, name=None, stereo="unresolved", geometry="unresolved",
             repr_status="unresolved", stereo_raw=None, identity_note=None):
    ev = [E[p] for p in pages]
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    m["names"]["original"] = f(name or label, name or label, ev)
    m["names"]["english"] = missing(ev)
    m["names"]["abbreviations"] = missing(ev)
    m["names"]["aliases"] = missing(ev)
    m["cas_numbers"] = missing(ev)
    m["representation_type"] = "unresolved"
    m["smiles"] = missing(ev)
    m["generic_template"] = None
    why = "Status assigned conservatively from the page image; no complete source-checked molecular graph was encoded."
    m["stereochemistry"]["chirality_status"] = derived(stereo_raw or label, stereo, ev, why)
    m["stereochemistry"]["geometry_status"] = derived(label, geometry, ev, why)
    m["stereochemistry"]["representation_status"] = derived(label, repr_status, ev, why)
    for key in ("atom_centers", "double_bonds", "axial", "mixture", "assignment",
                "isomeric_smiles", "identity_key"):
        m["stereochemistry"][key] = missing(ev)
    if stereo_raw:
        m["stereochemistry"]["assignment"] = f(stereo_raw, stereo_raw, ev)
    m["identity_note"] = f(identity_note or f"Source label/drawing {label}; complete structure not encoded.",
                            identity_note or f"Source label/drawing {label}; complete structure not encoded.", ev)
    m["external_properties"] = []
    molecules[mid] = m
    for eid in ev:
        mol_evidence[mid].add(eid)
    return mid


def ensure(mid, name, pages=(5,)):
    if mid not in molecules:
        molecule(mid, name, pages, name=name, stereo="unresolved", geometry="unresolved",
                 repr_status="unresolved", identity_note="Named process component; no source structure drawing was checked in pages 1-6.")
    return mid


# Source compounds and drawn intermediates.
molecule("m01", "1", (1, 2, 3, 5), name="Wieland–Miescher ketone 1",
         stereo="partially_specified", repr_status="not_encoded",
         stereo_raw="Solid-wedge bridgehead methyl is drawn; absolute assignment is not stated on pages 1-6.")
molecule("m02", "2", (2, 6), name="hexa-2,5-dione (compound 2)")
# Compound 2 is the one simple source-drawn substrate for which a complete
# candidate graph can be encoded without guessing a complex ring system.
_smi2 = "CC(=O)CCC(=O)C"
_diag2 = st.stereo_inspect(_smi2)
molecules["m02"]["representation_type"] = "smiles"
molecules["m02"]["smiles"] = f(None, _smi2, [E[2], E[6]], method="structure_derivation",
                                  derivation="Candidate SMILES derived from the source name and complete line drawing of compound 2.")
molecules["m02"]["stereochemistry"]["chirality_status"] = derived("standard stereo diagnostic: no potential stereocenters", "achiral", [E[2], E[6]], "RDKit standard stereo diagnostic on the candidate graph found no potential stereocenters; graph remains a machine candidate.")
molecules["m02"]["stereochemistry"]["geometry_status"] = derived("standard stereo diagnostic: carbonyl double bonds only", "not_applicable", [E[2], E[6]], "No E/Z-applicable alkene is present in the candidate graph.")
molecules["m02"]["stereochemistry"]["representation_status"] = derived("standard stereo candidate graph", "encoded", [E[2], E[6]], "The candidate graph contains no stereochemical elements requiring tags.")
molecules["m02"]["stereochemistry"]["isomeric_smiles"] = f(None, _diag2["canonical_isomeric_smiles"], [E[2], E[6]], method="structure_derivation", derivation="Generated by standard stereo from the candidate graph.")
molecules["m02"]["stereochemistry"]["identity_key"] = f(None, _diag2["stereo_graph_key"], [E[2], E[6]], method="structure_derivation", derivation="Generated by standard stereo; not a sample identity.")
molecule("m03", "3", (1, 2, 6), name="6-(but-3-yn-1-yl)-6-methyl-7-methylene-1,4-dioxaspiro[4.5]decane (3)")
molecule("m04", "4", (2,), name="in-situ hydrazone reagent 4")
molecule("m05", "5", (1, 2, 6), name="hydrazone precursor 5")
molecule("m06", "6", (1, 2), name="anticipated tetracyclic compound 6",
         stereo="partially_specified", repr_status="not_encoded", stereo_raw="Relative wedges/dashes are shown in Figure 2.")
molecule("m07", "7", (1, 2), name="known conidiogenone B intermediate 7",
         stereo="partially_specified", repr_status="not_encoded", stereo_raw="Relative wedges/dashes are shown in Figure 2.")
molecule("m08", "8", (1, 2), name="conidiogenone B (8)", stereo="partially_specified",
         repr_status="not_encoded", stereo_raw="Figure 2 gives a stereochemical drawing without a pages-1-6 absolute-name assignment.")
molecule("m09", "9", (1, 2, 5), name="8a-methyl-3,4,8,8a-tetrahydro-2H-spiro[naphthalene-1,2′-[1,3]dioxolan]-6(7H)-one (9)",
         stereo="partially_specified", repr_status="not_encoded", stereo_raw="Bridgehead methyl is drawn as a solid wedge.")
molecule("m10", "10", (1, 2, 5, 6), name="(rel-1aS,4aR,8aS)-epoxide 10",
         stereo="specified", repr_status="not_encoded", stereo_raw="The experimental heading states rel-1aS,4aR,8aS and the drawing uses stereobonds.")
molecule("m11", "11", (1, 2, 5, 6), name="6-(but-3-yn-1-yl)-6-methyl-1,4-dioxaspiro[4.5]decan-7-one (11)")
molecule("m12", "12", (1, 2, 6), name="4-(2-methyl-1,3-dithian-2-yl)butan-2-one (12)")
molecule("m13", "13", (1, 2, 6), name="allene 13", stereo="unspecified", repr_status="not_encoded",
         stereo_raw="An allene is drawn, but no axial assignment is given on pages 1-6.")
molecules["m13"]["stereochemistry"]["axial"] = f("allene C=C=C in source drawing; configuration not assigned",
                                                        [{"type": "allene", "configuration": "unspecified"}], [E[2], E[6]])
molecule("m14", "14", (1, 2, 6), name="allenyl ketone 14", stereo="unspecified", repr_status="not_encoded",
         stereo_raw="An allene is drawn, but no axial assignment is given on pages 1-6.")
molecules["m14"]["stereochemistry"]["axial"] = f("allene C=C=C in source drawing; configuration not assigned",
                                                        [{"type": "allene", "configuration": "unspecified"}], [E[2], E[6]])
for n in (15, 16, 17):
    molecule(f"m{n:02d}", str(n), (2,), name=f"proposed TMM-pathway intermediate {n}",
             stereo="unresolved", repr_status="not_encoded")
molecule("m18", "18", (2, 6), name="bridged cycloaddition product 18",
         stereo="specified", repr_status="not_encoded",
         stereo_raw="Scheme 2 shows relative stereobonds; single-crystal X-ray is stated to confirm the structure.")
molecule("m19", "19", (2, 6), name="1,4-HAT product 19", stereo="partially_specified",
         repr_status="not_encoded", stereo_raw="Reported as a single diastereomer, but relative stereochemistry was not determined.")
molecule("m20", "20", (2, 3), name="2-(but-3-yn-1-yl)-2-methylcyclohexane-1,3-dione (20)")
molecule("m21", "21", (2, 3), name="protected alcohol 21", stereo="specified", repr_status="not_encoded",
         stereo_raw="Scheme 4 shows relative stereobonds and reports dr 5:1.")
molecule("m22", "22", (3,), name="allenyl dithiane 22", stereo="partially_specified", repr_status="not_encoded",
         stereo_raw="Scheme 4 preserves the stereobond at the protected alcohol center; allene axis is unassigned.")
molecules["m22"]["stereochemistry"]["axial"] = f("allene C=C=C in source drawing; configuration not assigned",
                                                        [{"type": "allene", "configuration": "unspecified"}], [E[3]])
molecule("m23", "23", (2, 3), name="1,4-HAT product 23", stereo="partially_specified", repr_status="not_encoded",
         stereo_raw="Scheme 4 gives a stereochemical drawing, without a pages-1-6 absolute assignment.")
molecule("m24", "24", (2, 3), name="bridged cycloaddition product 24", stereo="specified", repr_status="not_encoded",
         stereo_raw="Scheme 4 gives relative stereobonds and key HMBC correlations.")
for n in (25, 27):
    molecule(f"m{n}", str(n), (3, 4), name=f"model TMM diyl {n}", stereo="unresolved", repr_status="not_encoded")
molecule("m26", "26", (3, 4), name="bridged model product 26", stereo="partially_specified", repr_status="not_encoded",
         stereo_raw="Table 1 footnote states that relative stereochemistry of the bridgehead methyl group is not determined.")
molecule("m28", "28", (3, 4), name="linear model product 28", stereo="partially_specified", repr_status="not_encoded",
         stereo_raw="Table 1 gives a stereochemical drawing; no absolute assignment is stated.")
molecule("m31", "31", (4,), name="vinyl-sulfone TMM precursor 31", stereo="unspecified", repr_status="not_encoded",
         stereo_raw="Allene precursor drawing has no axial assignment.")
molecule("m32", "32", (4,), name="TMM diyl 32", stereo="unresolved", repr_status="not_encoded")
molecule("m33", "33", (4,), name="linear cycloaddition product 33", stereo="partially_specified", repr_status="not_encoded",
         stereo_raw="Scheme 5 gives relative stereobonds; no absolute assignment is stated.")

# Unlabelled intermediates needed to attach stage-specific reported yields.
molecule("m_s4_epoxide", "unlabelled epoxidation intermediate in Scheme 4", (3,),
         name="unlabelled epoxidation intermediate in Scheme 4")
molecule("m_21_methylene", "unlabelled methylenated intermediate from 21", (3,),
         name="unlabelled methylenated intermediate from 21")
molecule("m_22_ketone", "unlabelled allenyl ketone from 22", (3,), name="unlabelled allenyl ketone from 22",
         stereo="unspecified", repr_status="not_encoded", stereo_raw="Allene axis is not assigned in Scheme 4.")
molecule("m_22_hydrazone", "unlabelled hydrazone from the allenyl ketone", (3,), name="unlabelled hydrazone from the allenyl ketone")
molecule("m_31_hydrazone", "unlabelled hydrazone from 31", (4,), name="unlabelled hydrazone from 31")

# Literature-example substrates/products in Table 2.  They are kept separate
# from this paper's experiments and are not merged by the generic drawings.
RVALUES = {"a": "Ac", "b": "CH2OH", "c": "CH(OMe)2", "d": "CMe(OMe)2", "e": "CMe(OCH2CH2O)"}
for suffix, rval in RVALUES.items():
    molecule(f"m29{suffix}", f"29{suffix}", (3,), name=f"literature precursor 29{suffix} (R = {rval})")
    molecule(f"m30{suffix}_tmm", f"30{suffix}", (3,), name=f"literature TMM diyl 30{suffix} (R = {rval})")
    molecule(f"m30{suffix}_linear", f"linear product from 30{suffix}", (3,), name=f"linear literature product, R = {rval}",
             stereo="partially_specified", repr_status="not_encoded", stereo_raw="Generic product drawing gives relative H stereobonds.")
    molecule(f"m30{suffix}_bridged", f"bridged product from 30{suffix}", (3,), name=f"bridged literature product, R = {rval}",
             stereo="partially_specified", repr_status="not_encoded", stereo_raw="Generic product drawing gives relative H stereobonds.")


# Named process components.  No CAS values are supplied or inferred.
for mid, name, pages in [
    ("ethylene_glycol", "ethylene glycol", (2, 5)), ("sieves_4A", "4 Å molecular sieves", (5, 6)),
    ("pTsOH_H2O", "p-toluenesulfonic acid monohydrate", (5,)), ("NaHCO3_aq", "saturated aqueous sodium bicarbonate", (5, 6)),
    ("EtOAc", "ethyl acetate", (5, 6)), ("MgSO4", "anhydrous magnesium sulfate", (5, 6)),
    ("silica", "silica gel", (5, 6)), ("Et2O", "diethyl ether", (5, 6)), ("hexane", "hexane", (5, 6)),
    ("MeOH", "methanol", (2, 3, 4, 5, 6)), ("H2O2_30", "30% (w/w) aqueous hydrogen peroxide", (5,)),
    ("NaOH_10", "10% (w/w) aqueous sodium hydroxide", (5,)), ("Na2S2O3_aq", "saturated aqueous sodium thiosulfate", (5,)),
    ("EtOH", "ethanol", (2, 3, 5)), ("TsNHNH2", "p-toluenesulfonyl hydrazide", (2, 3, 4, 5, 6)),
    ("AcOH", "acetic acid", (2, 3, 5, 6)), ("K2CO3_aq", "saturated aqueous potassium carbonate", (5,)),
    ("DCM", "dichloromethane", (2, 3, 5, 6)), ("Na2CO3", "sodium carbonate", (5,)), ("water", "water", (2, 3, 5, 6)),
    ("PPh3MeBr", "methyltriphenylphosphonium bromide", (2, 3, 6)), ("KOtBu", "potassium tert-butoxide", (2, 3, 6)),
    ("THF", "tetrahydrofuran", (2, 3, 6)), ("NH4Cl_aq", "saturated aqueous ammonium chloride", (6,)),
    ("propanedithiol", "1,3-propanedithiol", (2, 6)), ("BF3OEt2", "boron trifluoride diethyl etherate", (2, 6)),
    ("dioxane", "1,4-dioxane", (1, 2, 3, 6)), ("CuI", "copper(I) iodide", (1, 2, 3, 6)),
    ("LiOtBu", "lithium tert-butoxide", (1, 2, 3, 6)), ("MeCN", "acetonitrile", (2, 3, 6)),
    ("CaCO3", "calcium carbonate", (2, 3, 6)), ("MeI", "iodomethane", (2, 3, 6)),
    ("NaH", "sodium hydride (60% w/w dispersion in mineral oil)", (2, 3, 4, 6)), ("toluene", "toluene", (2, 3, 4, 6)),
    ("Ar", "argon", (5, 6)), ("LTBA", "lithium tri-tert-butoxyaluminum hydride", (2, 3)),
    ("TBSOTf", "tert-butyldimethylsilyl trifluoromethanesulfonate", (2, 3)), ("lutidine", "2,6-lutidine", (3,)),
]:
    ensure(mid, name, pages)


def charge(raw, normalized, ev):
    return f(raw, normalized, ev)


def participant(r, pid, mid, role, ev, *, stage="reaction", product_role=None,
                mass=None, amount=None, volume=None, equiv=None, mol_percent=None,
                reference=None):
    p = st.new_participant(pid, mid, role, ev, stage=stage)
    if product_role is not None:
        p["product_role"] = f(product_role, product_role, ev)
    values = {"mass": mass, "amount": amount, "volume": volume,
              "equivalents": equiv, "mol_percent": mol_percent}
    for key, item in values.items():
        if item is not None:
            raw, norm = item
            p["charges"][key] = charge(raw, norm, ev)
    if reference is not None:
        p["equivalence_reference"] = f(reference, reference, ev)
    r["participants"].append(p)
    mol_evidence[mid].update(ev)
    return pid


def step(r, sid, ev, operation, *, temp=None, time=None, pids=(), stage="reaction", atmosphere=None):
    s = st.new_step(sid, ev, stage=stage)
    s["participant_ids"] = list(pids)
    s["operation"] = f(operation, operation, ev)
    if temp is not None:
        s["temperature"] = f(temp[0], temp[1], ev)
    if time is not None:
        s["time"] = f(time[0], time[1], ev)
    if atmosphere is not None:
        s["atmosphere"] = f(atmosphere, atmosphere, [E[5]], method="inherited_general_condition")
    r["steps"].append(s)
    return s


def outcome(r, oid, kind, targets, ev, raw, normalized, *, role=None, yield_type=None,
            basis=None, scope=None):
    o = st.new_outcome(oid, kind, targets, ev)
    if role is not None:
        o["product_role"] = f(role, role, ev)
    o["value"] = f(raw, normalized, ev)
    if yield_type is not None:
        o["yield_type"] = f(yield_type, yield_type, ev)
    if basis is not None:
        o["measurement_basis"] = f(basis, basis, ev)
    if scope is not None:
        o["scope"] = f(scope, scope, ev)
    r["outcomes"].append(o)
    return o


CAT = {row["code"]: row["number"] for row in st.read(ROOT / "skills/chemical-reaction-extraction/references/bond_taxonomy.json")["categories"]}


def classification(r, ev, codes, note):
    changes = []
    labels = {}
    for i, code in enumerate(codes, 1):
        if code.startswith("DB-"):
            elems, after, ctype = code[3:].split("-"), 2, "formed"
        elif code.startswith("TB-"):
            elems, after, ctype = code[3:].split("-"), 3, "formed"
        else:
            elems, after, ctype = code.split("-"), 1, "hydrogen_gain" if "H" in code.split("-") else "formed"
        cid = f"b{i}"
        change = {"change_id": cid,
                  "atoms": [{"atom_map": None, "element": elems[0]}, {"atom_map": None, "element": elems[1]}],
                  "before_order": 0, "after_order": after, "change_type": ctype}
        if ctype == "hydrogen_gain":
            change["hydrogen_count"] = 1
        changes.append(f(None, change, ev, method="structure_derivation", derivation=note))
        labels.setdefault(code, {"category_number": CAT[code], "code": code, "supporting_change_ids": []})["supporting_change_ids"].append(cid)
    r["classification"]["bond_changes"] = changes
    r["classification"]["labels"] = f(None, list(labels.values()), ev, method="structure_derivation", derivation=note)
    r["classification"]["coverage"] = f(None, "partial", ev, method="inference", derivation=note)
    r["classification"]["limitations"] = [
        "Partial element-level bond-change inference from the source drawing; atom mapping was not asserted.",
        "Implicit-H/proton source and omitted coproduct mappings are unresolved.",
    ]
    r["classification"]["model_votes"] = []


def reaction(suffix, pages, *, record_type="specific_experiment", categories=("other",),
             name=None, transformation=None, location=None):
    ev = []
    for p in pages:
        ev.extend([E[p], T[p]])
    r = st.new_reaction(f"{DOC}:{suffix}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = f(list(categories), list(categories), ev)
    r["location"] = f(location or {"pages": list(pages)}, location or {"pages": list(pages)}, ev)
    if name:
        r["reaction_name"] = f(name, name, ev)
    if transformation:
        r["transformation_type"] = derived(transformation, transformation, ev, "Transformation type inferred from the source arrow and accompanying prose.")
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    return r, ev


def set_solvent(r, ev, composition, *, ratio=None, volume=None, concentration=None, basis=None):
    r["solvent_conditions"]["composition"] = f(composition, composition, ev)
    if ratio is not None:
        r["solvent_conditions"]["mixture_ratio"] = f(ratio, ratio, ev)
    if volume is not None:
        r["solvent_conditions"]["total_volume"] = f(volume[0], volume[1], ev)
    if concentration is not None:
        r["solvent_conditions"]["concentration"] = f(concentration[0], concentration[1], ev)
    if basis is not None:
        r["solvent_conditions"]["concentration_basis"] = f(basis, basis, ev)


def set_procedure(r, ev, text, addition=None, pretreatment=None):
    r["procedure"]["original_text"] = f(text, text, ev)
    if addition:
        r["procedure"]["addition_order"] = f(addition, addition, ev)
    if pretreatment:
        r["procedure"]["pretreatment"] = f(pretreatment, pretreatment, ev)


reactions = []

# 1 -> 9, detailed procedure on page 5 and arrow on page 2.
r, ev = reaction("r01_1_to_9", (1, 2, 5), name="chemoselective ketal protection", transformation="ketal formation", location="Scheme 1 and experimental section for 9")
p1 = participant(r, "r01_p1", "m01", "reactant", ev, mass=("1843 mg", q(1843, "mg")), amount=("10.134 mmol", q(10.134, "mmol")), equiv=("1 equiv.", q(1, "equiv")), reference="m01")
participant(r, "r01_p2", "ethylene_glycol", "reagent", ev, volume=("50.6 mL", q(50.6, "mL")))
participant(r, "r01_p3", "sieves_4A", "additive", ev)
participant(r, "r01_p4", "pTsOH_H2O", "reagent", ev, mass=("1928 mg", q(1928, "mg")), amount=("10.134 mmol", q(10.134, "mmol")), equiv=("1 equiv.", q(1, "equiv")), reference="m01")
participant(r, "r01_p5", "m09", "product", ev, product_role="main product", mass=("2071 mg", q(2071, "mg")))
for pid, mid, vol in [("r01_w1", "NaHCO3_aq", ("200 mL", q(200, "mL"))), ("r01_w2", "EtOAc", ("3 × 200 mL", q(200, "mL", multiplicity=3, raw_text="3 × 200 mL"))), ("r01_w3", "MgSO4", None), ("r01_w4", "silica", None), ("r01_w5", "Et2O", None), ("r01_w6", "hexane", None)]:
    participant(r, pid, mid, "workup_agent" if pid != "r01_w4" else "purification_agent", ev, stage="workup" if pid != "r01_w4" else "purification", volume=vol)
set_solvent(r, ev, "ethylene glycol", volume=("50.6 mL", q(50.6, "mL")), concentration=("0.2 M", q(0.2, "M")), basis="compound 1")
step(r, "r01_s1", ev, "Add molecular sieves and p-TsOH·H2O to 1 in ethylene glycol and stir.", temp=("room temperature", "room temperature"), time=("23 min", q(23, "min")), pids=(p1, "r01_p2", "r01_p3", "r01_p4"), atmosphere="Ar")
step(r, "r01_s2", ev, "Quench into saturated NaHCO3 with ice; extract, dry, filter and concentrate.", pids=("r01_w1", "r01_w2", "r01_w3", "r01_w4"), stage="workup")
step(r, "r01_s3", ev, "Dissolve crude syrup in minimal ether, add hexane until murky, and crystallize in the fridge overnight.", temp=("0 °C during hexane addition", q(0, "°C")), time=("overnight", "overnight"), pids=("r01_w5", "r01_w6"), stage="purification")
outcome(r, "r01_o1", "yield", ["m09"], ev, "92 %", q(92, "%"), role="main product", yield_type="isolated yield", basis="white crystals obtained by filtration", scope="1 to 9")
set_procedure(r, ev, "Ground 1 in ethylene glycol was treated with 4 Å sieves and p-TsOH·H2O for 23 min at room temperature, followed by aqueous workup and crystallization.")
classification(r, ev, ["C-O", "C-O"], "Ketal formation creates two C–O connections at the protected carbonyl carbon; complete atom mapping is not asserted.")
reactions.append(r)

# 9 -> 10.
r, ev = reaction("r02_9_to_10", (1, 2, 5), name="nucleophilic epoxidation", transformation="enone epoxidation", location="Scheme 1 and experimental section for 10")
participant(r, "r02_p1", "m09", "reactant", ev, mass=("459 mg", q(459, "mg")), amount=("2.067 mmol", q(2.067, "mmol")), equiv=("1 equiv.", q(1, "equiv")), reference="m09")
participant(r, "r02_p2", "MeOH", "solvent", ev, volume=("2.67 mL", q(2.67, "mL")))
participant(r, "r02_p3", "H2O2_30", "reagent", ev, volume=("0.61 mL", q(0.61, "mL")), amount=("5.994 mmol", q(5.994, "mmol")), equiv=("2.9 equiv.", q(2.9, "equiv")), reference="m09")
participant(r, "r02_p4", "NaOH_10", "reagent", ev, volume=("0.36 mL", q(0.36, "mL")), amount=("0.992 mmol", q(0.992, "mmol")), equiv=("0.48 equiv.", q(0.48, "equiv")), reference="m09")
participant(r, "r02_p5", "m10", "product", ev, product_role="main product", mass=("376 mg", q(376, "mg")))
for pid, mid, vol, role in [("r02_w1", "Na2S2O3_aq", ("200 mL", q(200, "mL")), "workup_agent"), ("r02_w2", "EtOAc", ("3 × 200 mL", q(200, "mL", multiplicity=3, raw_text="3 × 200 mL")), "workup_agent"), ("r02_w3", "MgSO4", None, "workup_agent"), ("r02_w4", "silica", None, "purification_agent"), ("r02_w5", "hexane", None, "purification_agent")]:
    participant(r, pid, mid, role, ev, stage="workup" if role == "workup_agent" else "purification", volume=vol)
set_solvent(r, ev, "methanol", volume=("2.67 mL", q(2.67, "mL")), concentration=("0.775 M", q(0.775, "M")), basis="compound 9")
step(r, "r02_s1", ev, "Add H2O2 solution and NaOH solution dropwise over 30 min.", temp=("0 °C", q(0, "°C")), time=("30 min addition", q(30, "min")), pids=("r02_p1", "r02_p2", "r02_p3", "r02_p4"), atmosphere="Ar")
step(r, "r02_s2", ev, "Stir after addition.", temp=("0 °C", q(0, "°C")), time=("6 h", q(6, "h")), pids=("r02_p1", "r02_p3", "r02_p4"), atmosphere="Ar")
step(r, "r02_s3", ev, "Hold in the cold room; add an additional 1 equiv H2O2 after each 24 h period.", temp=("1–5 °C", {"minimum": 1, "maximum": 5, "unit": "°C", "raw_text": "1–5 °C"}), time=("3 d", q(3, "d")), pids=("r02_p3",), atmosphere="Ar")
step(r, "r02_s4", ev, "Quench with saturated sodium thiosulfate, extract with ethyl acetate, dry, concentrate and flash chromatograph (ethyl acetate/hexane 1:3).", pids=("r02_w1", "r02_w2", "r02_w3", "r02_w4", "r02_w5"), stage="workup")
outcome(r, "r02_o1", "yield", ["m10"], ev, "76 %", q(76, "%"), role="main product", yield_type="isolated yield", basis="isolated after flash chromatography", scope="9 to 10")
outcome(r, "r02_o2", "yield", ["m10"], ev, "BRSM 80 %", q(80, "%"), role="main product", yield_type="yield based on recovered starting material", scope="9 to 10")
set_procedure(r, ev, "9 was treated with 30% H2O2 and 10% NaOH at 0 °C, held 6 h at 0 °C and 3 d at 1–5 °C with periodic additional H2O2, then worked up and chromatographed.")
classification(r, ev, ["C-O", "C-O"], "Epoxidation forms two C–O connections across the source alkene; the accompanying C=C order decrease is not a new sigma-bond label.")
reactions.append(r)

# 10 -> 11, primary procedure.
r, ev = reaction("r03_10_to_11_primary", (1, 2, 5), name="Eschenmoser–Tanabe fragmentation", transformation="epoxide fragmentation to gamma-ethynyl ketone", location="Scheme 1 and first experimental procedure for 11")
participant(r, "r03_p1", "m10", "reactant", ev, mass=("140 mg", q(140, "mg")), amount=("0.589 mmol", q(0.589, "mmol")), equiv=("1 equiv.", q(1, "equiv")), reference="m10")
participant(r, "r03_p2", "EtOH", "solvent", ev, volume=("5.9 mL", q(5.9, "mL")))
participant(r, "r03_p3", "TsNHNH2", "reagent", ev)
participant(r, "r03_p4", "AcOH", "reagent", ev, volume=("33.7 μL", q(33.7, "uL")), amount=("0.589 mmol", q(0.589, "mmol")), equiv=("1 equiv.", q(1, "equiv")), reference="m10")
participant(r, "r03_p5", "m11", "product", ev, product_role="main product", mass=("101 mg", q(101, "mg")))
for pid, mid, vol, role in [("r03_w1", "K2CO3_aq", ("100 mL", q(100, "mL")), "workup_agent"), ("r03_w2", "EtOAc", ("3 × 100 mL", q(100, "mL", multiplicity=3, raw_text="3 × 100 mL")), "workup_agent"), ("r03_w3", "MgSO4", None, "workup_agent"), ("r03_w4", "silica", None, "purification_agent"), ("r03_w5", "hexane", None, "purification_agent")]:
    participant(r, pid, mid, role, ev, stage="workup" if role == "workup_agent" else "purification", volume=vol)
set_solvent(r, ev, "ethanol", volume=("5.9 mL", q(5.9, "mL")), concentration=("0.1 M", q(0.1, "M")), basis="compound 10")
step(r, "r03_s1", ev, "Add p-toluenesulfonyl hydrazide to 10 in ethanol and stir.", temp=("room temperature", "room temperature"), time=("3 h", q(3, "h")), pids=("r03_p1", "r03_p2", "r03_p3"), atmosphere="Ar")
step(r, "r03_s2", ev, "Add acetic acid dropwise and continue stirring.", temp=("room temperature", "room temperature"), time=("12 h", q(12, "h")), pids=("r03_p4",), atmosphere="Ar")
step(r, "r03_s3", ev, "Pour into saturated K2CO3, extract, dry, concentrate and flash chromatograph (ethyl acetate/hexane 1:7).", pids=("r03_w1", "r03_w2", "r03_w3", "r03_w4", "r03_w5"), stage="workup")
outcome(r, "r03_o1", "yield", ["m11"], ev, "77 %", q(77, "%"), role="main product", yield_type="isolated yield", basis="isolated after flash chromatography", scope="10 to 11")
set_procedure(r, ev, "10 and p-toluenesulfonyl hydrazide were stirred in ethanol for 3 h, acetic acid was added, and the mixture was stirred 12 h before workup.")
classification(r, ev, ["TB-C-C"], "The Eschenmoser–Tanabe fragmentation produces the terminal alkyne shown in 11; C–C cleavage and other changes are not fully atom-mapped.")
r["review"]["unresolved"].append("The amount of p-toluenesulfonyl hydrazide is not reported for the primary procedure and remains unreported.")
reactions.append(r)

# Alternative 10 -> 11 procedure crossing pages 5-6.
r, ev = reaction("r04_10_to_11_alternative", (5, 6), name="alternative Eschenmoser–Tanabe fragmentation", transformation="epoxide fragmentation to gamma-ethynyl ketone", location="Alternative procedure for 11, spanning pages 5-6")
participant(r, "r04_p1", "m10", "reactant", ev, mass=("108 mg", q(108, "mg")), amount=("0.453 mmol", q(0.453, "mmol")), equiv=("1 equiv.", q(1, "equiv")), reference="m10")
participant(r, "r04_p2", "DCM", "solvent", ev, volume=("0.36 mL", q(0.36, "mL")))
participant(r, "r04_p3", "TsNHNH2", "reagent", ev, mass=("89 mg", q(89, "mg")), amount=("0.462 mmol", q(0.462, "mmol")), equiv=("1 equiv.", q(1, "equiv")), reference="m10")
participant(r, "r04_p4", "AcOH", "reagent", ev, volume=("34 μL", q(34, "uL")), amount=("0.589 mmol", q(0.589, "mmol")), equiv=("1.35 equiv.", q(1.35, "equiv")), reference="m10")
participant(r, "r04_p5", "Na2CO3", "reagent", ev, mass=("115 mg", q(115, "mg")), amount=("1.359 mmol", q(1.359, "mmol")), equiv=("3 equiv.", q(3, "equiv")), reference="m10")
participant(r, "r04_p6", "MeOH", "solvent", ev, volume=("0.54 mL", q(0.54, "mL")))
participant(r, "r04_p7", "m11", "product", ev, product_role="main product", mass=("71 mg", q(71, "mg")))
for pid, mid, vol, role in [("r04_w1", "NaHCO3_aq", ("100 mL", q(100, "mL")), "workup_agent"), ("r04_w2", "DCM", ("3 × 100 mL", q(100, "mL", multiplicity=3, raw_text="3 × 100 mL")), "workup_agent"), ("r04_w3", "water", ("100 mL", q(100, "mL")), "workup_agent"), ("r04_w4", "MgSO4", None, "workup_agent"), ("r04_w5", "silica", None, "purification_agent"), ("r04_w6", "hexane", None, "purification_agent")]:
    participant(r, pid, mid, role, ev, stage="workup" if role == "workup_agent" else "purification", volume=vol)
set_solvent(r, ev, "dichloromethane, then methanol added", basis="compound 10")
r["solvent_conditions"]["concentration"] = f(
    ["initial 1.26 M in dichloromethane", "final concentration 0.5 M after methanol addition"],
    [{"stage": "initial", **q(1.26, "M")}, {"stage": "after methanol addition", **q(0.5, "M")}], ev)
step(r, "r04_s1", ev, "Add p-toluenesulfonyl hydrazide and acetic acid dropwise to 10 in dichloromethane.", temp=("0 °C", q(0, "°C")), time=("30 min", q(30, "min")), pids=("r04_p1", "r04_p2", "r04_p3", "r04_p4"), atmosphere="Ar")
step(r, "r04_s2", ev, "Warm and stir.", temp=("room temperature", "room temperature"), time=("1.5 h", q(1.5, "h")), pids=("r04_p1", "r04_p3", "r04_p4"), atmosphere="Ar")
step(r, "r04_s3", ev, "Add sodium carbonate and methanol, then stir.", temp=("40 °C", q(40, "°C")), time=("2 h", q(2, "h")), pids=("r04_p5", "r04_p6"), atmosphere="Ar")
step(r, "r04_s4", ev, "Pour into saturated NaHCO3; extract with dichloromethane, wash with water, dry, concentrate and flash chromatograph (ethyl acetate/hexane 1:7).", pids=("r04_w1", "r04_w2", "r04_w3", "r04_w4", "r04_w5", "r04_w6"), stage="workup")
outcome(r, "r04_o1", "yield", ["m11"], ev, "71 %", q(71, "%"), role="main product", yield_type="isolated yield", basis="isolated after flash chromatography", scope="10 to 11; alternative procedure")
set_procedure(r, ev, "Alternative procedure: 10, p-toluenesulfonyl hydrazide and acetic acid were held at 0 °C then room temperature; Na2CO3 and MeOH were added and the mixture was held at 40 °C before workup.")
classification(r, ev, ["TB-C-C"], "The alternative fragmentation produces the terminal alkyne shown in 11; complete atom mapping is unresolved.")
r["review"]["unresolved"].append("The source reports two stage-specific concentrations (1.26 M initial and final 0.5 M); they are retained as separate stage values rather than collapsed into one concentration.")
reactions.append(r)

# 11 -> 3.
r, ev = reaction("r05_11_to_3", (1, 2, 6), name="Wittig methylenation", transformation="carbonyl methylenation", location="Scheme 1 and experimental section for 3")
participant(r, "r05_p1", "PPh3MeBr", "reagent", ev, mass=("1410 mg", q(1410, "mg")), amount=("3.867 mmol", q(3.867, "mmol")), equiv=("3 equiv.", q(3, "equiv")), reference="m11")
participant(r, "r05_p2", "KOtBu", "reagent", ev, mass=("434 mg", q(434, "mg")), amount=("3.867 mmol", q(3.867, "mmol")), equiv=("3 equiv.", q(3, "equiv")), reference="m11")
participant(r, "r05_p3", "m11", "reactant", ev, mass=("287 mg", q(287, "mg")), amount=("1.289 mmol", q(1.289, "mmol")), equiv=("1 equiv.", q(1, "equiv")), reference="m11")
participant(r, "r05_p4", "THF", "solvent", ev, volume=("5 mL + 1.4 mL substrate solution", q(6.4, "mL", raw_text="5 mL + 1.4 mL")))
participant(r, "r05_p5", "m03", "product", ev, product_role="main product", mass=("223 mg", q(223, "mg")))
for pid, mid, vol, role in [("r05_w1", "NH4Cl_aq", ("150 mL", q(150, "mL")), "workup_agent"), ("r05_w2", "EtOAc", ("3 × 150 mL", q(150, "mL", multiplicity=3, raw_text="3 × 150 mL")), "workup_agent"), ("r05_w3", "MgSO4", None, "workup_agent"), ("r05_w4", "silica", None, "purification_agent"), ("r05_w5", "hexane", None, "purification_agent")]:
    participant(r, pid, mid, role, ev, stage="workup" if role == "workup_agent" else "purification", volume=vol)
set_solvent(r, ev, "THF", concentration=("final concentration 0.2 M", q(0.2, "M")), basis="compound 11")
step(r, "r05_s0", ev, "Dry methyltriphenylphosphonium bromide by heating.", temp=("110 °C", q(110, "°C")), time=("overnight", "overnight"), pids=("r05_p1",), atmosphere="Ar")
step(r, "r05_s1", ev, "Add KOtBu to the dry phosphonium salt in THF and stir.", temp=("0 °C", q(0, "°C")), time=("30 min", q(30, "min")), pids=("r05_p1", "r05_p2", "r05_p4"), atmosphere="Ar")
step(r, "r05_s2", ev, "Add 11 in THF by cannula at 0 °C, then stir at room temperature.", temp=("room temperature after addition at 0 °C", "room temperature"), time=("24 h", q(24, "h")), pids=("r05_p3",), atmosphere="Ar")
step(r, "r05_s3", ev, "Quench with saturated NH4Cl, extract, dry, concentrate and flash chromatograph (ethyl acetate/hexane 1:20).", pids=("r05_w1", "r05_w2", "r05_w3", "r05_w4", "r05_w5"), stage="workup")
outcome(r, "r05_o1", "yield", ["m03"], ev, "78 %", q(78, "%"), role="main product", yield_type="isolated yield", basis="isolated after flash chromatography", scope="11 to 3")
set_procedure(r, ev, "A Wittig reagent was generated from dried methyltriphenylphosphonium bromide and KOtBu in THF; 11 was added at 0 °C and the mixture was stirred 24 h at room temperature.", pretreatment="Methyltriphenylphosphonium bromide was heated at 110 °C overnight.")
classification(r, ev, ["DB-C-C"], "Wittig methylenation forms the exocyclic C=C bond shown in product 3; atom mapping and coproducts are unresolved.")
reactions.append(r)

# 2 -> 12.
r, ev = reaction("r06_2_to_12", (2, 6), name="dithiane formation", transformation="carbonyl thioacetalization", location="Scheme 1 and experimental section for 12")
participant(r, "r06_p1", "m02", "reactant", ev, volume=("0.286 mL", q(0.286, "mL")), amount=("1.972 mmol", q(1.972, "mmol")), equiv=("2 equiv.", q(2, "equiv")), reference="propanedithiol")
participant(r, "r06_p2", "propanedithiol", "reactant", ev, volume=("0.1 mL", q(0.1, "mL")), amount=("0.986 mmol", q(0.986, "mmol")), equiv=("1 equiv.", q(1, "equiv")), reference="propanedithiol")
participant(r, "r06_p3", "DCM", "solvent", ev, volume=("9.9 mL", q(9.9, "mL")))
participant(r, "r06_p4", "BF3OEt2", "reagent", ev, volume=("0.262 mL", q(0.262, "mL")), amount=("0.986 mmol", q(0.986, "mmol")), equiv=("1 equiv.", q(1, "equiv")), reference="propanedithiol")
participant(r, "r06_p5", "m12", "product", ev, product_role="main product", mass=("214 mg", q(214, "mg")))
for pid, mid, vol, role in [("r06_w1", "NaHCO3_aq", ("200 mL", q(200, "mL")), "workup_agent"), ("r06_w2", "DCM", ("3 × 200 mL", q(200, "mL", multiplicity=3, raw_text="3 × 200 mL")), "workup_agent"), ("r06_w3", "MgSO4", None, "workup_agent"), ("r06_w4", "silica", None, "purification_agent"), ("r06_w5", "hexane", None, "purification_agent")]:
    participant(r, pid, mid, role, ev, stage="workup" if role == "workup_agent" else "purification", volume=vol)
set_solvent(r, ev, "dichloromethane", volume=("9.9 mL", q(9.9, "mL")), concentration=("0.1 M", q(0.1, "M")), basis="1,3-propanedithiol")
step(r, "r06_s1", ev, "Add BF3·OEt2 dropwise to hexa-2,5-dione and 1,3-propanedithiol in dichloromethane and stir.", temp=("0 °C", q(0, "°C")), time=("30 min", q(30, "min")), pids=("r06_p1", "r06_p2", "r06_p3", "r06_p4"), atmosphere="Ar")
step(r, "r06_s2", ev, "Quench with sodium bicarbonate solution, extract, dry, concentrate and flash chromatograph (ethyl acetate/hexane 1:7).", pids=("r06_w1", "r06_w2", "r06_w3", "r06_w4", "r06_w5"), stage="workup")
outcome(r, "r06_o1", "yield", ["m12"], ev, "quant.", "quantitative", role="main product", yield_type="isolated yield", basis="isolated after flash chromatography", scope="2 plus 1,3-propanedithiol to 12")
set_procedure(r, ev, "Hexa-2,5-dione and 1,3-propanedithiol were treated with BF3·OEt2 for 30 min at 0 °C, then worked up and chromatographed.")
classification(r, ev, ["C-S", "C-S"], "Dithiane formation creates two C–S connections at the protected carbonyl center; full mapping is unresolved.")
reactions.append(r)

# 3 + 12 -> 13.
r, ev = reaction("r07_3_plus_12_to_13", (1, 2, 6), name="Wang allenation", transformation="copper-mediated allene formation", location="Scheme 1 and experimental section for 13")
participant(r, "r07_p1", "m12", "reactant", ev, mass=("205 mg", q(205, "mg")), amount=("1.004 mmol", q(1.004, "mmol")), equiv=("2.7 equiv.", q(2.7, "equiv")), reference="m03")
participant(r, "r07_p2", "TsNHNH2", "reagent", ev, mass=("200 mg", q(200, "mg")), amount=("1.042 mmol", q(1.042, "mmol")), equiv=("2.8 equiv.", q(2.8, "equiv")), reference="m03")
participant(r, "r07_p3", "sieves_4A", "additive", ev)
participant(r, "r07_p4", "m03", "reactant", ev, mass=("81.9 mg", q(81.9, "mg")), amount=("0.372 mmol", q(0.372, "mmol")), equiv=("1 equiv.", q(1, "equiv")), reference="m03")
participant(r, "r07_p5", "dioxane", "solvent", ev, volume=("1 mL then 5 mL carrier", q(6, "mL", raw_text="1 mL + 5 mL")))
participant(r, "r07_p6", "CuI", "reagent", ev, mass=("73 mg", q(73, "mg")), amount=("0.383 mmol", q(0.383, "mmol")), equiv=("1 equiv.", q(1, "equiv")), reference="m03")
participant(r, "r07_p7", "LiOtBu", "reagent", ev, mass=("164 mg", q(164, "mg")), amount=("2.046 mmol", q(2.046, "mmol")), equiv=("5.5 equiv.", q(5.5, "equiv")), reference="m03")
participant(r, "r07_p8", "m13", "product", ev, product_role="main product", mass=("123 mg", q(123, "mg")))
participant(r, "r07_w1", "silica", "purification_agent", ev, stage="purification")
participant(r, "r07_w2", "EtOAc", "purification_agent", ev, stage="purification")
participant(r, "r07_w3", "hexane", "purification_agent", ev, stage="purification")
set_solvent(r, ev, "1,4-dioxane", concentration=("final concentration 0.08 M", q(0.08, "M")), basis="compound 3")
step(r, "r07_s1", ev, "Combine 12, p-toluenesulfonyl hydrazide and 4 Å sieves in 1,4-dioxane, reflux, then cool.", temp=("reflux", "reflux"), time=("12 h", q(12, "h")), pids=("r07_p1", "r07_p2", "r07_p3", "r07_p5"), atmosphere="Ar")
step(r, "r07_s2", ev, "Add 3 in dioxane by cannula, CuI, additional sieves and LiOtBu; reflux and cool.", temp=("reflux", "reflux"), time=("24 h", q(24, "h")), pids=("r07_p4", "r07_p6", "r07_p7"), atmosphere="Ar")
step(r, "r07_s3", ev, "Filter through a silica pad and flash chromatograph (ethyl acetate/hexane 1:30 to 1:20).", pids=("r07_w1", "r07_w2", "r07_w3"), stage="purification")
outcome(r, "r07_o1", "yield", ["m13"], ev, "81 %", q(81, "%"), role="main product", yield_type="isolated yield", basis="isolated after flash chromatography", scope="3 plus in-situ 4 (from 12) to allene 13")
set_procedure(r, ev, "One-pot hydrazone formation from 12 followed by CuI/LiOtBu coupling with 3 in refluxing 1,4-dioxane afforded allene 13.")
classification(r, ev, ["DB-C-C"], "The Wang allenation creates the new allene C=C connection shown in 13; no full atom mapping or axial assignment is asserted.")
reactions.append(r)

# 13 -> 14.
r, ev = reaction("r08_13_to_14", (1, 2, 6), name="dithiane deprotection", transformation="dithiane-to-ketone conversion", location="Scheme 1 and experimental section for 14")
participant(r, "r08_p1", "m13", "reactant", ev, mass=("1368 mg", q(1368, "mg")), amount=("3.348 mmol", q(3.348, "mmol")), equiv=("1 equiv.", q(1, "equiv")), reference="m13")
participant(r, "r08_p2", "MeCN", "solvent", ev)
participant(r, "r08_p3", "water", "solvent", ev)
participant(r, "r08_p4", "CaCO3", "reagent", ev, mass=("1460 mg", q(1460, "mg")), amount=("14.416 mmol", q(14.416, "mmol")), equiv=("4 equiv.", q(4, "equiv")), reference="m13")
participant(r, "r08_p5", "MeI", "reagent", ev, volume=("2.25 mL", q(2.25, "mL")), amount=("36.036 mmol", q(36.036, "mmol")), equiv=("10 equiv.", q(10, "equiv")), reference="m13")
participant(r, "r08_p6", "m14", "product", ev, product_role="main product", mass=("989 mg", q(989, "mg")))
for pid, mid, vol, role in [("r08_w1", "NH4Cl_aq", ("250 mL", q(250, "mL")), "workup_agent"), ("r08_w2", "EtOAc", ("3 × 250 mL", q(250, "mL", multiplicity=3, raw_text="3 × 250 mL")), "workup_agent"), ("r08_w3", "MgSO4", None, "workup_agent"), ("r08_w4", "silica", None, "purification_agent"), ("r08_w5", "hexane", None, "purification_agent")]:
    participant(r, pid, mid, role, ev, stage="workup" if role == "workup_agent" else "purification", volume=vol)
set_solvent(r, ev, "MeCN/H2O", ratio="3:1 (v/v)", volume=("36 mL", q(36, "mL")), concentration=("0.1 M", q(0.1, "M")), basis="compound 13")
step(r, "r08_s1", ev, "Add CaCO3 and iodomethane to 13 in MeCN/H2O at room temperature, then heat and stir.", temp=("80 °C", q(80, "°C")), time=("4 h", q(4, "h")), pids=("r08_p1", "r08_p2", "r08_p3", "r08_p4", "r08_p5"), atmosphere="Ar")
step(r, "r08_s2", ev, "Cool, quench with saturated NH4Cl, extract, dry, concentrate and flash chromatograph (ethyl acetate/hexane 1:5).", pids=("r08_w1", "r08_w2", "r08_w3", "r08_w4", "r08_w5"), stage="workup")
outcome(r, "r08_o1", "yield", ["m14"], ev, "93 %", q(93, "%"), role="main product", yield_type="isolated yield", basis="isolated after flash chromatography", scope="13 to 14")
set_procedure(r, ev, "13, CaCO3 and iodomethane were heated in MeCN/H2O (3:1) at 80 °C for 4 h, then worked up and chromatographed.")
classification(r, ev, ["DB-C-O"], "Dithiane deprotection forms the ketone C=O connection in 14; sulfur-containing coproducts and full mapping are unresolved.")
reactions.append(r)

# 14 -> 18 + 19, including hydrazone and TMM stages.
r, ev = reaction("r09_14_to_18_19", (1, 2, 6), name="TMM diyl-mediated tandem cycloaddition", transformation="hydrazone-to-TMM tandem cyclization with competing 1,4-HAT", location="Scheme 2 and experimental section beginning on page 6")
participant(r, "r09_p1", "m14", "reactant", ev, mass=("120 mg", q(120, "mg")), amount=("0.376 mmol", q(0.376, "mmol")), equiv=("1 equiv.", q(1, "equiv")), reference="m14")
participant(r, "r09_p2", "MeOH", "solvent", ev, volume=("0.75 mL", q(0.75, "mL")))
participant(r, "r09_p3", "TsNHNH2", "reagent", ev, mass=("87 mg", q(87, "mg")), amount=("0.451 mmol", q(0.451, "mmol")), equiv=("1.2 equiv.", q(1.2, "equiv")), reference="m14")
participant(r, "r09_p4a", "m05", "product", ev, product_role="isolated stage-1 hydrazone intermediate", mass=("183 mg", q(183, "mg")))
participant(r, "r09_p4", "m05", "reactant", ev, product_role="isolated intermediate", mass=("183 mg", q(183, "mg")), amount=("0.376 mmol in second stage", q(0.376, "mmol")), equiv=("1 equiv. in second stage", q(1, "equiv")), reference="m05")
participant(r, "r09_p5", "NaH", "reagent", ev, mass=("18.8 mg", q(18.8, "mg")), amount=("0.470 mmol", q(0.470, "mmol")), equiv=("1.25 equiv.", q(1.25, "equiv")), reference="m05")
participant(r, "r09_p6", "toluene", "solvent", ev, volume=("44 mL + 10 mL substrate solution", q(54, "mL", raw_text="44 mL + 10 mL")))
participant(r, "r09_p7", "m18", "product", ev, product_role="major bridged product", mass=("60 mg", q(60, "mg")))
participant(r, "r09_p8", "m19", "product", ev, product_role="1,4-HAT product", mass=("49 mg", q(49, "mg")))
for pid, mid, vol, role in [("r09_w1", "silica", None, "purification_agent"), ("r09_w2", "EtOAc", None, "purification_agent"), ("r09_w3", "hexane", None, "purification_agent"), ("r09_w4", "NH4Cl_aq", ("100 mL", q(100, "mL")), "workup_agent"), ("r09_w5", "EtOAc", ("3 × 100 mL", q(100, "mL", multiplicity=3, raw_text="3 × 100 mL")), "workup_agent"), ("r09_w6", "MgSO4", None, "workup_agent"), ("r09_w7", "Et2O", None, "purification_agent")]:
    participant(r, pid, mid, role, ev, stage="workup" if role == "workup_agent" else "purification", volume=vol)
set_solvent(r, ev, "methanol in hydrazone stage; toluene in cyclization stage", basis="14 then hydrazone 5")
r["solvent_conditions"]["concentration"] = f(
    ["0.5 M in hydrazone stage", "final concentration 0.007 M in cyclization stage"],
    [{"stage": "hydrazone formation", **q(0.5, "M")}, {"stage": "cyclization", **q(0.007, "M")}], ev)
step(r, "r09_s1", ev, "Add p-toluenesulfonyl hydrazide to 14 in methanol and stir.", temp=("room temperature", "room temperature"), time=("1.5 h", q(1.5, "h")), pids=("r09_p1", "r09_p2", "r09_p3"), atmosphere="Ar")
step(r, "r09_s2", ev, "Concentrate and flash chromatograph (ethyl acetate/hexane 1:4 to 1:3) to isolate hydrazone 5.", pids=("r09_w1", "r09_w2", "r09_w3"), stage="purification")
step(r, "r09_s3", ev, "Add freshly prepared hydrazone in toluene by cannula to NaH suspension at 0 °C; degas with an Ar balloon.", temp=("0 °C", q(0, "°C")), time=("15 min degassing", q(15, "min")), pids=("r09_p4", "r09_p5", "r09_p6"), atmosphere="Ar")
step(r, "r09_s4", ev, "Reflux the reaction mixture.", temp=("reflux", "reflux"), time=("24 h", q(24, "h")), pids=("r09_p4", "r09_p5", "r09_p6"), atmosphere="Ar")
step(r, "r09_s5", ev, "Cool, quench with saturated NH4Cl, extract, dry, concentrate and flash chromatograph (diethyl ether/hexane 1:150).", pids=("r09_w4", "r09_w5", "r09_w6", "r09_w7", "r09_w3"), stage="workup")
outcome(r, "r09_o0", "yield", ["m05"], ev, "quantitative", "quantitative", role="isolated hydrazone intermediate", yield_type="isolated intermediate yield", basis="isolated after flash chromatography", scope="14 to hydrazone 5")
outcome(r, "r09_o1", "yield", ["m18"], ev, "54 %", q(54, "%"), role="major bridged product", yield_type="isolated yield", basis="isolated after flash chromatography", scope="hydrazone 5 to 18")
outcome(r, "r09_o2", "yield", ["m19"], ev, "44 %", q(44, "%"), role="1,4-HAT product", yield_type="isolated yield", basis="isolated after flash chromatography", scope="hydrazone 5 to 19")
outcome(r, "r09_o3", "selectivity", ["m19"], ev, "single diastereomer; relative stereochemistry not determined", "single diastereomer; relative stereochemistry not determined", role="1,4-HAT product", scope="diastereomeric outcome for 19; configuration order unavailable")
set_procedure(r, ev, "14 was converted to isolated hydrazone 5 in methanol; freshly made 5 was added to NaH in toluene at 0 °C, degassed, refluxed 24 h, worked up and chromatographed to give 18 and 19.")
classification(r, ev, ["C-C", "C-C", "C-H"], "The tandem cycloaddition forms new C–C ring connections; the competing 1,4-HAT outcome includes a newly formed C–H connection. Product-specific atom mapping is unresolved.")
r["review"]["unresolved"].append("A single reaction execution gives products 18 and 19 through different branches; element-level classification is the union of both observed outcomes and lacks product-specific atom maps.")
reactions.append(r)

# Scheme 4: 1 -> 20 multistage sequence.
r, ev = reaction("r10_1_to_20_scheme4", (2, 3), name="direct nucleophilic epoxidation and Eschenmoser–Tanabe fragmentation", transformation="two-stage epoxidation/fragmentation", location="Scheme 4")
participant(r, "r10_p1", "m01", "reactant", ev, reference="m01")
participant(r, "r10_p2", "H2O2_30", "reagent", ev)
participant(r, "r10_p3", "NaOH_10", "reagent", ev)
participant(r, "r10_p4", "MeOH", "solvent", ev)
participant(r, "r10_p5", "m_s4_epoxide", "product", ev, product_role="stage-1 intermediate")
participant(r, "r10_p6", "TsNHNH2", "reagent", ev)
participant(r, "r10_p7", "AcOH", "reagent", ev)
participant(r, "r10_p8", "EtOH", "solvent", ev)
participant(r, "r10_p9", "m20", "product", ev, product_role="stage-2 product")
set_solvent(r, ev, "MeOH in stage i; EtOH in stage ii")
step(r, "r10_s1", ev, "Treat 1 with H2O2 and NaOH in MeOH.", temp=("0 °C", q(0, "°C")), pids=("r10_p1", "r10_p2", "r10_p3", "r10_p4"), atmosphere="Ar")
step(r, "r10_s2", ev, "Treat the epoxidation product with TsNHNH2 and AcOH in EtOH.", temp=("RT", "room temperature"), pids=("r10_p5", "r10_p6", "r10_p7", "r10_p8"), atmosphere="Ar")
outcome(r, "r10_o1", "yield", ["m_s4_epoxide"], ev, "80 %", q(80, "%"), role="stage-1 intermediate", scope="stage i: epoxidation of 1")
outcome(r, "r10_o2", "yield", ["m20"], ev, "60 %", q(60, "%"), role="stage-2 product", scope="stage ii: fragmentation to 20")
set_procedure(r, ev, "Scheme 4: i) H2O2, NaOH, MeOH, 0 °C, 80%; ii) TsNHNH2, AcOH, EtOH, RT, 60%.")
classification(r, ev, ["TB-C-C"], "The overall fragmentation sequence produces the terminal alkyne drawn in 20; transient epoxide C–O bonds are stage-specific and not asserted as final-product bonds.")
r["review"]["unresolved"].append("Scheme 4 does not report stage times, reagent charges, workup, purification, or yield measurement bases on pages 1-6.")
reactions.append(r)

# 20 -> 21.
r, ev = reaction("r11_20_to_21_scheme4", (2, 3), name="ketone reduction and TBS protection", transformation="diastereoselective reduction/silylation", location="Scheme 4")
participant(r, "r11_p1", "m20", "reactant", ev, reference="m20")
participant(r, "r11_p2", "LTBA", "reagent", ev)
participant(r, "r11_p3", "THF", "solvent", ev)
participant(r, "r11_p4", "TBSOTf", "reagent", ev)
participant(r, "r11_p5", "lutidine", "reagent", ev)
participant(r, "r11_p6", "DCM", "solvent", ev)
participant(r, "r11_p7", "m21", "product", ev, product_role="protected alcohol")
step(r, "r11_s1", ev, "Reduce 20 with LiAlH(OtBu)3 in THF.", temp=("−60 °C", q(-60, "°C")), pids=("r11_p1", "r11_p2", "r11_p3"), atmosphere="Ar")
step(r, "r11_s2", ev, "Protect the resulting alcohol with TBSOTf and 2,6-lutidine in DCM.", temp=("0 °C", q(0, "°C")), pids=("r11_p4", "r11_p5", "r11_p6"), atmosphere="Ar")
outcome(r, "r11_o1", "yield", ["m21"], ev, "85 % over two steps", q(85, "%"), role="protected alcohol", scope="20 to 21 over reduction and TBS protection")
outcome(r, "r11_o2", "selectivity", ["m21"], ev, "dr 5:1", "5:1", role="protected alcohol", scope="diastereomeric ratio; major/minor configurations are not fully assigned on pages 1-6")
set_solvent(r, ev, "THF in stage i; DCM in stage ii")
set_procedure(r, ev, "Scheme 4: i) LiAlH(OtBu)3, THF, −60 °C (dr 5:1); ii) TBSOTf, 2,6-lutidine, DCM, 0 °C; 85% over two steps.")
classification(r, ev, ["C-H", "O-Si"], "Net reduction introduces a C–H bond at the former carbonyl carbon and silylation forms O–Si; hydride/proton sources and full maps are unresolved.")
r["review"]["unresolved"].append("Scheme 4 reports no times, charges, workup, purification, selectivity method, or yield measurement basis on pages 1-6.")
reactions.append(r)

# 21 -> 22.
r, ev = reaction("r12_21_to_22_scheme4", (2, 3), name="methylenation and Wang allenation sequence", transformation="two-stage methylenation/allene formation", location="Scheme 4")
participant(r, "r12_p1", "m21", "reactant", ev, reference="m21")
participant(r, "r12_p2", "PPh3MeBr", "reagent", ev)
participant(r, "r12_p3", "KOtBu", "reagent", ev)
participant(r, "r12_p4", "THF", "solvent", ev)
participant(r, "r12_p5", "m_21_methylene", "product", ev, product_role="stage-1 intermediate")
participant(r, "r12_p6", "m12", "reactant", ev)
participant(r, "r12_p7", "TsNHNH2", "reagent", ev)
participant(r, "r12_p8", "dioxane", "solvent", ev)
participant(r, "r12_p9", "CuI", "reagent", ev)
participant(r, "r12_p10", "LiOtBu", "reagent", ev)
participant(r, "r12_p11", "m22", "product", ev, product_role="stage-2 allene product")
step(r, "r12_s1", ev, "Methylenate 21 with PPh3MeBr/KOtBu in THF.", temp=("RT", "room temperature"), pids=("r12_p1", "r12_p2", "r12_p3", "r12_p4"), atmosphere="Ar")
step(r, "r12_s2", ev, "Generate the hydrazone from 12 and TsNHNH2 in 1,4-dioxane, then couple with the methylenated intermediate using CuI/LiOtBu.", temp=("reflux", "reflux"), pids=("r12_p5", "r12_p6", "r12_p7", "r12_p8", "r12_p9", "r12_p10"), atmosphere="Ar")
outcome(r, "r12_o1", "yield", ["m_21_methylene"], ev, "97 %", q(97, "%"), role="stage-1 intermediate", scope="Wittig methylenation of 21")
outcome(r, "r12_o2", "yield", ["m22"], ev, "75 %", q(75, "%"), role="stage-2 allene product", scope="coupling of methylenated 21 with reagent generated from 12")
set_solvent(r, ev, "THF in stage i; 1,4-dioxane in stage ii")
set_procedure(r, ev, "Scheme 4: i) PPh3MeBr, KOtBu, THF, RT, 97%; ii) 12, TsNHNH2, 1,4-dioxane, reflux; CuI, LiOtBu, reflux, 75%.")
classification(r, ev, ["DB-C-C", "DB-C-C"], "The Wittig stage and Wang allenation stage each form a new C=C connection; exact atom maps are unresolved.")
r["review"]["unresolved"].append("Scheme 4 reports no times, reagent charges, workup, purification, or yield measurement bases on pages 1-6.")
reactions.append(r)

# 22 -> 23 + 24.
r, ev = reaction("r13_22_to_23_24_scheme4", (2, 3), name="modified TMM tandem cycloaddition", transformation="dithiane deprotection/hydrazone formation/TMM cyclization", location="Scheme 4")
for pid, mid, role, prole in [("r13_p1", "m22", "reactant", None), ("r13_p2", "MeI", "reagent", None), ("r13_p3", "CaCO3", "reagent", None), ("r13_p4", "MeCN", "solvent", None), ("r13_p5", "water", "solvent", None), ("r13_p6", "m_22_ketone", "product", "stage-1 intermediate"), ("r13_p7", "TsNHNH2", "reagent", None), ("r13_p8", "MeOH", "solvent", None), ("r13_p9", "m_22_hydrazone", "product", "stage-2 intermediate"), ("r13_p10", "NaH", "reagent", None), ("r13_p11", "toluene", "solvent", None), ("r13_p12", "m23", "product", "1,4-HAT product"), ("r13_p13", "m24", "product", "bridged cycloaddition product")]:
    participant(r, pid, mid, role, ev, product_role=prole, reference="m22" if pid == "r13_p1" else None)
step(r, "r13_s1", ev, "Convert dithiane 22 to the allenyl ketone using MeI/CaCO3 in MeCN/H2O 3:1.", temp=("80 °C", q(80, "°C")), pids=("r13_p1", "r13_p2", "r13_p3", "r13_p4", "r13_p5"), atmosphere="Ar")
step(r, "r13_s2", ev, "Form the hydrazone with TsNHNH2 in MeOH.", temp=("RT", "room temperature"), pids=("r13_p6", "r13_p7", "r13_p8"), atmosphere="Ar")
step(r, "r13_s3", ev, "Treat the hydrazone with NaH in toluene under reflux.", temp=("reflux", "reflux"), pids=("r13_p9", "r13_p10", "r13_p11"), atmosphere="Ar")
outcome(r, "r13_o1", "yield", ["m_22_ketone"], ev, "92 %", q(92, "%"), role="stage-1 intermediate", scope="22 to allenyl ketone")
outcome(r, "r13_o2", "yield", ["m_22_hydrazone"], ev, "66 %", q(66, "%"), role="stage-2 intermediate", scope="allenyl ketone to hydrazone")
outcome(r, "r13_o3", "yield", ["m23", "m24"], ev, "76 %", q(76, "%"), role="combined products", scope="hydrazone cyclization; combined 23 and 24")
outcome(r, "r13_o4", "selectivity", ["m23", "m24"], ev, "23:24 = 4:15", "4:15", role="product ratio", scope="product ratio in order 23:24")
set_solvent(r, ev, "MeCN/H2O 3:1 in stage i; MeOH in stage ii; toluene in stage iii", ratio="MeCN/H2O = 3:1")
set_procedure(r, ev, "Scheme 4: i) MeI, CaCO3, MeCN/H2O=3:1, 80 °C, 92%; ii) TsNHNH2, MeOH, RT, 66%; iii) NaH, toluene, reflux, 76%; products 23:24 = 4:15.")
classification(r, ev, ["C-C", "C-C", "C-H"], "The cycloaddition branch forms new C–C ring connections, while the competing 1,4-HAT product includes a new C–H connection; stage/product atom maps are unresolved.")
r["review"]["unresolved"].append("Scheme 4 does not report times, reagent charges, workup, purification, yield/selectivity measurement bases, or separate isolated yields for 23 and 24 on pages 1-6.")
reactions.append(r)

# Table 1 current model experiments.
for suffix, src, product, yld, result in [
    ("r14_25_to_26_table1", "m25", "m26", 57, "bridged product only"),
    ("r15_27_to_28_table1", "m27", "m28", 88, "linear product only"),
]:
    r, ev = reaction(suffix, (3, 4), categories=("mechanistic_study",), name="intramolecular [3+2] TMM diyl–olefin cycloaddition", transformation="model TMM cycloaddition", location="Table 1; duplicated as entries 6-7 in Table 2 and discussed with Table 3")
    participant(r, suffix+"_p1", src, "reactant", ev)
    participant(r, suffix+"_p2", product, "product", ev, product_role=result)
    step(r, suffix+"_s1", ev, "Intramolecular [3+2] cycloaddition of the drawn model TMM diyl; experimental conditions are not shown in Table 1.", pids=(suffix+"_p1",))
    outcome(r, suffix+"_o1", "yield", [product], ev, f"{yld} %", q(yld, "%"), role=result, scope=f"model TMM diyl to {result}")
    outcome(r, suffix+"_o2", "selectivity", [product], ev, result, result, role=result, scope="reported regioselectivity; no ratio beyond 'only'")
    set_procedure(r, ev, f"Table 1 reports {src[1:]} to {product[1:]} in {yld}% with {result}; conditions and measurement basis are not printed in the table.")
    classification(r, ev, ["C-C", "C-C"], "The intramolecular [3+2] cycloaddition forms two new C–C ring connections in the drawn product; atom mapping is unresolved.")
    r["review"]["unresolved"].append("Table 1 does not state precursor generation, solvent, temperature, time, reagent charges, workup, purification, yield method, or selectivity method.")
    reactions.append(r)

# Scheme 5 current experiment.
r, ev = reaction("r16_31_to_33_scheme5", (4,), categories=("mechanistic_study",), name="EWG-modified TMM tandem cycloaddition", transformation="hydrazone formation and TMM cycloaddition", location="Scheme 5")
participant(r, "r16_p1", "m31", "reactant", ev, reference="m31")
participant(r, "r16_p2", "TsNHNH2", "reagent", ev)
participant(r, "r16_p3", "MeOH", "solvent", ev)
participant(r, "r16_p4", "m_31_hydrazone", "product", ev, product_role="stage-1 hydrazone intermediate")
participant(r, "r16_p5", "NaH", "reagent", ev)
participant(r, "r16_p6", "toluene", "solvent", ev)
participant(r, "r16_p7", "m33", "product", ev, product_role="linear product")
step(r, "r16_s1", ev, "Treat 31 with TsNHNH2 in MeOH.", temp=("RT", "room temperature"), pids=("r16_p1", "r16_p2", "r16_p3"), atmosphere="Ar")
step(r, "r16_s2", ev, "Treat the hydrazone with NaH in toluene under reflux to generate TMM diyl 32 and cyclize.", temp=("reflux", "reflux"), pids=("r16_p4", "r16_p5", "r16_p6"), atmosphere="Ar")
outcome(r, "r16_o1", "yield", ["m_31_hydrazone"], ev, "80 %", q(80, "%"), role="stage-1 hydrazone intermediate", scope="31 to unlabelled hydrazone")
outcome(r, "r16_o2", "yield", ["m33"], ev, "16 %", q(16, "%"), role="linear product", scope="hydrazone via TMM diyl 32 to 33")
outcome(r, "r16_o3", "selectivity", ["m33"], ev, "linear product", "linear product", role="linear product", scope="regioselectivity; no linear/bridged ratio reported")
set_solvent(r, ev, "MeOH in stage i; toluene in stage ii")
set_procedure(r, ev, "Scheme 5: i) TsNHNH2, MeOH, RT, 80%; ii) NaH, toluene, reflux, 16%; TMM diyl 32 leads to linear product 33.")
classification(r, ev, ["C-C", "C-C"], "The intramolecular [3+2] cycloaddition forms two new C–C ring connections in linear product 33; atom mapping is unresolved.")
r["review"]["unresolved"].append("Scheme 5 does not state times, charges, workup, purification, yield measurement bases, or a numerical regioisomer ratio on page 4.")
reactions.append(r)

# Table 2 prior-work rows, explicitly literature examples rather than this paper's executions.
LIT = [
    ("a", "7.3:1", 98, None),
    ("b", "1:1.9", 91, None),
    ("c", "1:2.0", 86, "yield over two steps: acetalization and cycloaddition"),
    ("d", "1:15.7", 68, "yield over two steps: cycloaddition and ketal deprotection"),
    ("e", "1:10", 89, None),
]
for idx, (sfx, ratio, yld, foot) in enumerate(LIT, 17):
    rid = f"r{idx:02d}_literature_29{sfx}"
    r, ev = reaction(rid, (3,), record_type="literature_example", categories=("other",), name="reported intramolecular TMM diyl cycloaddition", transformation="literature TMM cycloaddition", location=f"Table 2 entry {idx-16}; cited prior work [14]")
    participant(r, rid+"_p1", f"m29{sfx}", "reactant", ev)
    participant(r, rid+"_p2", "THF", "solvent", ev)
    participant(r, rid+"_p2a", f"m30{sfx}_tmm", "unknown", ev, product_role="transient TMM diyl named in Table 2")
    participant(r, rid+"_p3", f"m30{sfx}_linear", "product", ev, product_role="linear regioisomer")
    participant(r, rid+"_p4", f"m30{sfx}_bridged", "product", ev, product_role="bridged regioisomer")
    step(r, rid+"_s1", ev, "Heat the depicted literature precursor in THF under reflux; Table 2 attributes the example to reference 14.", temp=("reflux", "reflux"), pids=(rid+"_p1", rid+"_p2"))
    outcome(r, rid+"_o1", "yield", [f"m30{sfx}_linear", f"m30{sfx}_bridged"], ev, f"{yld} %" + (f" ({foot})" if foot else ""), q(yld, "%"), role="combined regioisomers", scope=foot or "combined linear and bridged products")
    outcome(r, rid+"_o2", "selectivity", [f"m30{sfx}_linear", f"m30{sfx}_bridged"], ev, f"linear/bridged = {ratio}", ratio, role="regioisomer ratio", scope="ratio order is linear:bridged")
    set_solvent(r, ev, "THF")
    set_procedure(r, ev, f"Table 2 entry {idx-16}: precursor 29{sfx}, R={RVALUES[sfx]}, THF reflux; linear/bridged {ratio}; {yld}%" + (f"; {foot}." if foot else "."))
    classification(r, ev, ["C-C", "C-C"], "The literature [3+2] cycloaddition drawings form two new C–C ring connections; atom mapping is unresolved.")
    r["review"]["unresolved"].append("This is a summarized cited-literature example, not an execution in the current paper; time, charges, workup, purification and measurement methods are not given in Table 2.")
    reactions.append(r)

# Proposed mechanistic pathway is a separate non-experimental record.
r, ev = reaction("r22_proposed_tmm_mechanism", (2, 3), record_type="proposed_mechanism", categories=("mechanistic_study",), name="proposed singlet/triplet TMM pathway", transformation="proposed TMM formation, intersystem crossing, cyclization and 1,4-HAT", location="Schemes 2 and 3")
for pid, mid, role in [("r22_p1", "m14", "reactant"), ("r22_p2", "m05", "unknown"), ("r22_p3", "m15", "unknown"), ("r22_p4", "m16", "unknown"), ("r22_p5", "m17", "unknown"), ("r22_p6", "m18", "product"), ("r22_p7", "m19", "product")]:
    participant(r, pid, mid, role, ev, product_role="proposed pathway species" if role in ("unknown", "product") else None)
step(r, "r22_s1", ev, "Authors propose hydrazone-derived diazo/TMM formation, singlet/triplet pathways, 5-exo/6-endo closure, radical recombination or 1,4-HAT.", pids=tuple(f"r22_p{i}" for i in range(1, 8)))
set_procedure(r, ev, "Schemes 2-3 are a mechanistic layout, not an independently executed synthetic procedure.")
classification(r, ev, ["C-C", "C-C", "C-H"], "Proposed product-forming branches contain C–C cyclization and a C–H-forming 1,4-HAT branch; these changes are mechanistic proposals and not atom-mapped proof.")
r["review"]["unresolved"].append("Proposed intermediates and spin-state pathways are author interpretation, not independently verified experiments; no quantitative conditions/outcomes are assigned here.")
reactions.append(r)


# Expand molecule occurrence evidence to every page on which it is used.
for r in reactions:
    for p in r["participants"]:
        mol_evidence[p["molecule_id"]].update(r["evidence_ids"])
for mid, m in molecules.items():
    evs = sorted(mol_evidence[mid] or set(m["occurrences"][0]["evidence_ids"]))
    m["occurrences"][0]["evidence_ids"] = evs


dataset = st.read(RUN / "dataset.json")
bundle = st.empty(dataset["documents"])
bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": f"{DOC}:ctx_page1_prior_syntheses", "document_id": DOC, "kind": "cited_prior_work",
     "evidence_ids": [E[1], T[1]], "description": "Figure 1 shows six natural products from cited prior TMM methods; it does not supply concrete experimental conditions for new executions in this paper."},
    {"context_id": f"{DOC}:ctx_pages1_2_retrosynthesis", "document_id": DOC, "kind": "retrosynthetic_analysis",
     "evidence_ids": [E[1], T[1], E[2], T[2]], "description": "Figure 2 is a retrosynthetic plan linking conidiogenone B (8), known intermediate 7, anticipated 6, and precursors; retrosynthetic arrows are not counted as experiments."},
    {"context_id": f"{DOC}:ctx_page3_fmo", "document_id": DOC, "kind": "computational_discussion",
     "evidence_ids": [E[3], T[3]], "description": "Table 2 HOMO(TMM)-LUMO(olefin) gaps are calculated descriptors; entries 1-5 are cited literature summaries and entries 6-7 repeat Table 1 experiments."},
    {"context_id": f"{DOC}:ctx_page4_barriers", "document_id": DOC, "kind": "computational_results",
     "evidence_ids": [E[4], T[4]], "description": "Table 3 and Figure 3 report calculated singlet cycloaddition barriers and transition-state distances; they are retained as computational context, not synthetic executions."},
    {"context_id": f"{DOC}:ctx_page5_general", "document_id": DOC, "kind": "general_experimental_information",
     "evidence_ids": [E[5], T[5]], "description": "General information states oven-dried glassware, Ar atmosphere, default room temperature, solvent handling, chromatography and analytical instrumentation. Ar is inherited only for this paper's experimental records."},
    {"context_id": f"{DOC}:ctx_pages5_6_characterization", "document_id": DOC, "kind": "characterization",
     "evidence_ids": [E[5], T[5], E[6], T[6]], "description": "NMR and HRMS data following preparations of 9, 10, 11, 3, 12, 13 and 14 are characterization, not separate synthetic experiments."},
]

raw_inventory = {
    "batch_id": "allene-stereo-v31-200-20261017",
    "document_id": DOC,
    "title": dataset["documents"][0]["title"],
    "pages_attempted": [1, 2, 3, 4, 5, 6],
    "reaction_inventory": [{"reaction_id": r["reaction_id"], "record_type": r["record_type"],
                            "evidence_ids": r["evidence_ids"], "unresolved": r["review"]["unresolved"]}
                           for r in reactions],
    "source_context": bundle["source_context"],
    "stereo_checks": {"m02": _diag2},
    "notes": [
        "Complex source drawings were retained by labels and stereochemical observations; speculative SMILES and atom maps were not created.",
        "No supporting-information file was supplied or associated.",
        "Page 6 begins the 14-to-18/19 procedure; its reaction conditions and isolated yields are visible on page 6, while later characterization is outside the claimed page range.",
    ],
}

(RUN / "raw_model").mkdir(parents=True, exist_ok=True)
(RUN / "raw_model/pages_0001_0006_extraction.json").write_text(json.dumps(raw_inventory, ensure_ascii=False, indent=2) + "\n")
(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title and abstract", "finding": "Article identity and the reported unexpected bridged selectivity; no standalone experiment arrow.", "evidence_ids": [E[1], T[1]]},
            {"region": "Figure 1", "finding": "Six cited prior syntheses presented as background; no conditions/yields sufficient for concrete executions in this paper.", "evidence_ids": [E[1], T[1]]},
            {"region": "Results and Discussion, lower right", "finding": "Narrative inventory of the Scheme 1 sequence 1→9→10→11→3, coupling with in-situ 4 from 12 to 13, and deprotection to 14.", "evidence_ids": [E[1], T[1]]}
        ], "record_ids": [f"{DOC}:r01_1_to_9", f"{DOC}:r02_9_to_10", f"{DOC}:r03_10_to_11_primary", f"{DOC}:r05_11_to_3", f"{DOC}:r06_2_to_12", f"{DOC}:r07_3_plus_12_to_13", f"{DOC}:r08_13_to_14"],
         "unresolved": ["Figure 1 does not provide enough experimental detail to reconstruct cited executions; it remains source_context."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Figure 2", "finding": "Retrosynthetic arrows and anticipated/known intermediates; planning context, not experimental executions.", "evidence_ids": [E[2], T[2]]},
            {"region": "Scheme 1", "finding": "Seven current-work synthetic arrows with reported conditions/yields, linked to detailed records and later procedures.", "evidence_ids": [E[2], T[2]]},
            {"region": "Scheme 2", "finding": "Actual tandem cycloaddition of 14 gives 18 (54%) and 19 (44%); intermediate arrows 15–17 are proposed mechanism.", "evidence_ids": [E[2], T[2]]},
            {"region": "Scheme 3", "finding": "Author-proposed singlet/triplet and radical pathway, including 5-exo/6-endo closure and 1,4-HAT; separately classified as proposed_mechanism.", "evidence_ids": [E[2], T[2]]},
            {"region": "lower-right prose", "finding": "Narrative preview of modified Scheme 4 and reported 21 dr 5:1 / products 23 and 24.", "evidence_ids": [E[2], T[2]]}
        ], "record_ids": [f"{DOC}:r01_1_to_9", f"{DOC}:r02_9_to_10", f"{DOC}:r03_10_to_11_primary", f"{DOC}:r05_11_to_3", f"{DOC}:r06_2_to_12", f"{DOC}:r07_3_plus_12_to_13", f"{DOC}:r08_13_to_14", f"{DOC}:r09_14_to_18_19", f"{DOC}:r10_1_to_20_scheme4", f"{DOC}:r11_20_to_21_scheme4", f"{DOC}:r12_21_to_22_scheme4", f"{DOC}:r13_22_to_23_24_scheme4", f"{DOC}:r22_proposed_tmm_mechanism"],
         "unresolved": ["The desired linear product in Scheme 3 is hypothetical and was not created as an observed-product experiment."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Scheme 4", "finding": "Four multistage current-work arrows: 1→20, 20→21, 21→22, and 22→23+24, with stage yields, dr/product ratio and several unlabelled intermediates.", "evidence_ids": [E[3], T[3]]},
            {"region": "Table 1", "finding": "Two current model experiments retained separately: 25→26 (57%, bridged only) and 27→28 (88%, linear only).", "evidence_ids": [E[3], T[3]]},
            {"region": "Table 2", "finding": "Entries 1–5 are cited literature examples (separate literature_example records); entries 6–7 duplicate Table 1 and are not double-counted. HOMO–LUMO gaps are computational context.", "evidence_ids": [E[3], T[3]]}
        ], "record_ids": [f"{DOC}:r10_1_to_20_scheme4", f"{DOC}:r11_20_to_21_scheme4", f"{DOC}:r12_21_to_22_scheme4", f"{DOC}:r13_22_to_23_24_scheme4", f"{DOC}:r14_25_to_26_table1", f"{DOC}:r15_27_to_28_table1"] + [f"{DOC}:r{i:02d}_literature_29{s}" for i, s in zip(range(17, 22), "abcde")] + [f"{DOC}:r22_proposed_tmm_mechanism"],
         "unresolved": ["Scheme 4 omits charges, times, workups and measurement bases on this page.", "Table 1 omits experimental conditions and measurement methods."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "Table 3", "finding": "Calculated singlet cycloaddition barriers and observed/predicted regioselectivity; computational context, not five new executions.", "evidence_ids": [E[4], T[4]]},
            {"region": "Figure 3", "finding": "Calculated transition-state structures/distances for 25-TS and 27-TS; not a synthetic experiment.", "evidence_ids": [E[4], T[4]]},
            {"region": "Scheme 5", "finding": "Current two-stage experiment 31→hydrazone→33 via proposed TMM diyl 32, with 80% stage yield and 16% final yield.", "evidence_ids": [E[4], T[4]]}
        ], "record_ids": [f"{DOC}:r14_25_to_26_table1", f"{DOC}:r15_27_to_28_table1", f"{DOC}:r16_31_to_33_scheme5"],
         "unresolved": ["Scheme 5 omits charges, times, workup, purification, yield measurement basis and a numerical regioisomer ratio."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "Experimental Section general information", "finding": "General Ar atmosphere, glassware, default room temperature, solvent/purification and analytical context; not a reaction record.", "evidence_ids": [E[5], T[5]]},
            {"region": "procedure for 9", "finding": "Detailed protection of 1 to 9 including charges, 23 min reaction, workup, crystallization, mass and 92% yield.", "evidence_ids": [E[5], T[5]]},
            {"region": "procedure for 10", "finding": "Detailed epoxidation of 9 to 10 with 30 min addition, 6 h at 0 °C, 3 d at 1–5 °C, periodic H2O2, workup and 76% (BRSM 80%).", "evidence_ids": [E[5], T[5]]},
            {"region": "procedure for 11", "finding": "Primary 10→11 procedure and beginning of alternative procedure; alternative continues on page 6.", "evidence_ids": [E[5], T[5]]}
        ], "record_ids": [f"{DOC}:r01_1_to_9", f"{DOC}:r02_9_to_10", f"{DOC}:r03_10_to_11_primary", f"{DOC}:r04_10_to_11_alternative"],
         "unresolved": ["Primary 10→11 procedure does not report the p-toluenesulfonyl hydrazide charge."]},
        {"page": 6, "state": "attempted", "regions": [
            {"region": "top-left continuation", "finding": "Completion of the alternative 10→11 procedure and its 71% yield.", "evidence_ids": [E[6], T[6]]},
            {"region": "procedures for 3 and 12", "finding": "Detailed Wittig methylenation 11→3 and dithiane formation 2→12, with all printed charges, ordered conditions, workups and yields.", "evidence_ids": [E[6], T[6]]},
            {"region": "procedure for 13", "finding": "Detailed one-pot hydrazone generation/copper-mediated allene formation 3+12→13, with two reflux stages and 81% yield.", "evidence_ids": [E[6], T[6]]},
            {"region": "procedures for 14 and 18/19", "finding": "Detailed 13→14 deprotection (93%) and the beginning/complete reaction-operation portion of 14→18+19 with 54%/44% outcomes; later characterization is outside claimed coverage.", "evidence_ids": [E[6], T[6]]},
            {"region": "NMR/HRMS blocks", "finding": "Characterization retained as source_context; no independent synthetic records created.", "evidence_ids": [E[6], T[6]]}
        ], "record_ids": [f"{DOC}:r04_10_to_11_alternative", f"{DOC}:r05_11_to_3", f"{DOC}:r06_2_to_12", f"{DOC}:r07_3_plus_12_to_13", f"{DOC}:r08_13_to_14", f"{DOC}:r09_14_to_18_19"],
         "unresolved": ["The detailed 18/19 characterization continues beyond page 6 and is not claimed as covered."]}
    ]
}
(RUN / "coverage").mkdir(parents=True, exist_ok=True)
(RUN / "coverage/pages_0001_0006.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(reactions),
                  "source_context": len(bundle["source_context"])}, ensure_ascii=False))
