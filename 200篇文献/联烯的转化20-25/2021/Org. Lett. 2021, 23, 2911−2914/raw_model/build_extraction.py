#!/usr/bin/env python3
"""Build the source-derived v3.1 submission for DOI 10.1021/acs.orglett.1c00554.

This is a session-local transcription program. It does not alter the bundled
schema, tools, taxonomy, or quality rules. All records remain needs_review.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DOC = "doc_0e80a55b596be81e"
P = {n: [f"{DOC}:page{n}:text", f"{DOC}:page{n}:image"] for n in range(1, 5)}
MISSING = object()


def f(raw=None, norm=MISSING, ev=(), method="original_extraction", status=None,
      derivation=None):
    if norm is MISSING:
        norm = raw
    return st.fact(raw, norm, ev, method=method, status=status,
                   derivation=derivation)


def setf(obj, key, raw, norm=MISSING, ev=(), method="original_extraction",
         status=None, derivation=None):
    obj[key] = f(raw, norm, ev, method, status, derivation)


def source_molecule(mid, label, ev, *, name=None, stereo="unresolved",
                    geometry="unresolved", representation="unresolved",
                    smiles=None, identity_note=None, assignment=None, axial=None,
                    double_bonds=None):
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    m["names"]["original"] = f(name or label, name or label, ev)
    m["names"]["english"] = f()
    m["names"]["abbreviations"] = f(label, [label], ev) if label else f()
    m["names"]["aliases"] = f()
    m["cas_numbers"] = f()
    m["stereochemistry"]["chirality_status"] = f(
        stereo, stereo, ev, "inference",
        derivation="Status follows the inspected source drawing/text; no absent stereochemical mark was converted into achiral or racemic."
    )
    m["stereochemistry"]["geometry_status"] = f(
        geometry, geometry, ev, "inference",
        derivation="Geometry status follows the inspected source drawing/text."
    )
    m["stereochemistry"]["representation_status"] = f(
        representation, representation, ev, "inference",
        derivation="Encoding status reflects the submitted molecular representation."
    )
    if assignment is not None:
        m["stereochemistry"]["assignment"] = f(assignment, assignment, ev)
    if axial is not None:
        m["stereochemistry"]["axial"] = f(axial, axial, ev)
    if double_bonds is not None:
        m["stereochemistry"]["double_bonds"] = f(double_bonds, double_bonds, ev)
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = f(name or label, smiles, ev, "structure_derivation",
                         derivation="Candidate connectivity derived from the unambiguous source name/abbreviation; parsing is not source verification.")
        diag = st.stereo_inspect(smiles)
        m["stereochemistry"]["isomeric_smiles"] = f(
            None, diag["canonical_isomeric_smiles"], ev, "structure_derivation",
            derivation="RDKit canonical isomeric SMILES calculated from the submitted candidate connectivity."
        )
        m["stereochemistry"]["identity_key"] = f(
            None, diag["stereo_graph_key"], ev, "structure_derivation",
            derivation="RDKit stereo graph key calculated from the submitted candidate connectivity."
        )
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = f()
    if identity_note:
        m["identity_note"] = f(identity_note, identity_note, ev)
    return m


def template_molecule(mid, label, ev, fixed, variables, *, stereo="unresolved",
                      geometry="unresolved", note=None):
    m = source_molecule(mid, label, ev, stereo=stereo, geometry=geometry,
                        representation="not_encoded", identity_note=note)
    m["representation_type"] = "generic_template"
    if stereo in ("specified", "partially_specified"):
        m["stereochemistry"]["axial"] = f(
            "source-drawn allene axis; absolute descriptor not printed",
            [{"type": "allene axis", "configuration": "source-drawn; R_a/S_a not named"}], ev
        )
    if geometry in ("specified", "partially_specified"):
        m["stereochemistry"]["double_bonds"] = f(
            "alkene geometry shown by the generic drawing without E/Z labels",
            [{"configuration": "as drawn; E/Z descriptor not printed"}], ev
        )
    m["generic_template"] = {
        "fixed_structure": f(fixed, fixed, ev),
        "variables": [
            {
                "label": v,
                "attachment_points": f(f"source drawing attachment for {v}",
                                       f"source drawing attachment for {v}", ev),
                "definition": f(f"variable {v} as drawn", f"variable {v} as drawn", ev),
                "allowed_values": f(),
                "constraints": f(),
            }
            for v in variables
        ],
        "instance_mapping": [],
        "reason": f("R-group/general structure in Scheme 1", "R-group/general structure in Scheme 1", ev),
    }
    return m


def q(number, unit, raw_text=None, qualifier=None):
    x = {"number": number, "unit": unit}
    if raw_text:
        x["raw_text"] = raw_text
    if qualifier:
        x["qualifier"] = qualifier
    return x


def participant(pid, mid, role, ev, *, amount=None, volume=None, equiv=None,
                molpct=None, ref=None, product_role=None):
    p = st.new_participant(pid, mid, role, ev)
    if product_role:
        p["product_role"] = f(product_role, product_role, ev)
    for key, val in (("amount", amount), ("volume", volume),
                     ("equivalents", equiv), ("mol_percent", molpct)):
        if val is not None:
            p["charges"][key] = f(val.get("raw_text", str(val)), val, ev)
    if ref:
        p["equivalence_reference"] = f(f"relative to {ref}", ref, ev)
    return p


def step(step_id, ev, pids, *, operation=None, temperature=None, time=None,
         stage="reaction", atmosphere=None, pressure=None):
    s = st.new_step(step_id, ev, stage)
    s["participant_ids"] = list(pids)
    if operation is not None:
        s["operation"] = f(operation, operation, ev)
    if temperature is not None:
        s["temperature"] = f(str(temperature), temperature, ev)
    if time is not None:
        s["time"] = f(str(time), time, ev)
    if atmosphere is not None:
        s["atmosphere"] = f(atmosphere, atmosphere, ev)
    if pressure is not None:
        s["pressure"] = f(pressure, pressure, ev)
    return s


def outcome(oid, kind, targets, ev, value, *, role=None, yield_type=None,
            basis=None, scope=None):
    o = st.new_outcome(oid, kind, targets, ev)
    o["value"] = f(value if isinstance(value, str) else str(value), value, ev)
    if role:
        o["product_role"] = f(role, role, ev)
    if yield_type:
        o["yield_type"] = f(yield_type, yield_type, ev)
    if basis:
        o["measurement_basis"] = f(basis, basis, ev)
    if scope:
        o["scope"] = f(scope, scope, ev)
    return o


def ratio_value(text, metric="er", order="major:minor; configurations unassigned"):
    a, b = text.replace(">", "").split(":")
    value = {"ratio": [int(a), int(b)], "metric": metric, "order": order,
             "raw_text": text}
    if text.startswith(">"):
        value["qualifier"] = "gt"
    return value


def unresolved_classification(r, ev, note):
    r["classification"]["mapped_reaction_smiles"] = f()
    r["classification"]["bond_changes"] = []
    r["classification"]["labels"] = f(None, [], ev, "inference",
                                          derivation=note)
    r["classification"]["coverage"] = f(
        None, "unresolved", ev, "inference", derivation=note
    )
    r["classification"]["limitations"] = [note]
    r["classification"]["model_votes"] = []


def cyclo_classification(r, ev):
    changes = [
        {"change_id": "b1", "atoms": [{"element": "C", "source_role": "diene terminus"}, {"element": "C", "source_role": "allene terminus"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b2", "atoms": [{"element": "C", "source_role": "other diene terminus"}, {"element": "C", "source_role": "allene central carbon"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b3", "atoms": [{"element": "C"}, {"element": "C"}], "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
        {"change_id": "b4", "atoms": [{"element": "C"}, {"element": "C"}], "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
        {"change_id": "b5", "atoms": [{"element": "C"}, {"element": "C"}], "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
        {"change_id": "b6", "atoms": [{"element": "C", "source_role": "internal diene atoms"}, {"element": "C", "source_role": "internal diene atoms"}], "before_order": 1, "after_order": 2, "change_type": "order_increased"},
    ]
    labels = [
        {"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1", "b2"]},
        {"category_number": 13, "code": "DB-C-C", "supporting_change_ids": ["b6"]},
    ]
    deriv = "Partial net bond changes inferred from the inspected intramolecular allene–diene [4+2] drawings; exact atom maps, hydrogen changes, and full coproduct balance are unresolved."
    r["classification"]["mapped_reaction_smiles"] = f()
    r["classification"]["bond_changes"] = [f(None, c, ev, "structure_derivation", derivation=deriv) for c in changes]
    r["classification"]["labels"] = f(None, labels, ev, "structure_derivation", derivation=deriv)
    r["classification"]["coverage"] = f(None, "partial", ev, "inference", derivation=deriv)
    r["classification"]["limitations"] = [
        "No exact atom-mapped reaction SMILES was asserted because the complex stereochemical drawings were not safely transcribed.",
        "Hydrogen and catalyst/coproduct atom sources are not mapped.",
    ]
    r["classification"]["model_votes"] = []


def simple_classification(r, ev, changes, labels, note):
    r["classification"]["mapped_reaction_smiles"] = f()
    r["classification"]["bond_changes"] = [f(None, c, ev, "structure_derivation", derivation=note) for c in changes]
    r["classification"]["labels"] = f(None, labels, ev, "structure_derivation", derivation=note)
    r["classification"]["coverage"] = f(None, "partial", ev, "inference", derivation=note)
    r["classification"]["limitations"] = [note, "No complete atom mapping or coproduct balance was asserted."]
    r["classification"]["model_votes"] = []


def reaction(rid, ev, *, record_type="specific_experiment", categories=None,
             location=None, name=None, transformation=None):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories or ["other"], categories or ["other"], ev)
    r["location"] = f(location or rid, location or rid, ev)
    if name:
        r["reaction_name"] = f(name, name, ev)
    if transformation:
        r["transformation_type"] = f(None, transformation, ev, "inference",
                                       derivation="Transformation type assigned from the inspected source arrow and accompanying text.")
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    return r


def set_solvent(r, ev, composition, volume=None):
    r["solvent_conditions"]["composition"] = f(composition, composition, ev)
    if volume is not None:
        r["solvent_conditions"]["total_volume"] = f(volume.get("raw_text", str(volume)), volume, ev)


def add_note(r, ev, text):
    r["notes"].append(f(text, text, ev))


def write_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n")


dataset = json.loads((RUN / "dataset.json").read_text())
bundle = json.loads((RUN / "submitted.json").read_text())
bundle["documents"] = dataset["documents"]

# Molecule registry. Generic Scheme 1 structures use templates; complex scope
# drawings remain source-identified but unresolved rather than unsafe SMILES.
molecules = []
num_ev = {n: P[1] if n <= 12 else P[2] if n <= 41 else P[3] for n in range(1, 70)}
for n in [1, 3, 4, 5, 6, 7, 8, 9]:
    molecules.append(template_molecule(
        f"{DOC}:mol:{n}", str(n), P[1], f"Scheme 1 generic structure labeled {n}", ["R1", "R2", "X"],
        stereo="specified" if n in (3, 4, 5, 7, 9) else "unresolved",
        geometry="partially_specified" if n in (6, 7, 8, 9) else "unresolved",
        note="Generic R-group drawing; no concrete SMILES was asserted."
    ))

specific_numbers = sorted(set([2, 10, 11, 12, 13] + list(range(14, 50)) + list(range(51, 70, 2))))
for n in specific_numbers:
    ev = num_ev[n]
    if n == 11:
        molecules.append(source_molecule(f"{DOC}:mol:11", "11", ev,
            name="3-methyl-2-butanone (structure candidate from drawing)",
            stereo="achiral", geometry="not_applicable", representation="not_applicable",
            smiles="CC(=O)C(C)C",
            identity_note="Simple ketone candidate; source drawing is legible but machine record is not structure-checked."))
        continue
    is_substrate = n in [10, 12] or (14 <= n <= 48 and n % 2 == 0)
    is_product = n == 13 or (15 <= n <= 49 and n % 2 == 1) or n >= 51
    stereo = "specified" if (is_substrate or is_product or n == 2) else "unresolved"
    geometry = "partially_specified" if is_substrate else "specified" if is_product else "unresolved"
    rep = "not_encoded" if (is_substrate or is_product or n == 2) else "unresolved"
    assignment = None
    axial = None
    double = None
    note = "Complex source drawing retained by label; definitive source-safe SMILES and absolute descriptors were not asserted."
    if is_substrate:
        axial = [{"type": "allene axis", "configuration": "source-drawn enantiomer; R_a/S_a not named"}]
        assignment = "Optically active allene is drawn; absolute axial descriptor is not printed."
    if is_product:
        assignment = "Source drawing specifies relative stereochemistry at the fused-ring junction; absolute center descriptors are not printed."
        double = [{"type": "exocyclic alkene", "configuration": "as drawn; E/Z descriptor not printed"}]
    if n == 2:
        assignment = "Scheme 1 draws a stereodefined diphenylprolinol catalyst; center-by-center absolute descriptors are not printed."
    if n == 29:
        assignment = "Single-crystal X-ray: cis 5/6 ring junction; exocyclic double bond as drawn with the allene-derived group toward the five-membered ring."
        note += " Figure 2 and prose provide X-ray stereochemical evidence for product 29."
        ev = P[2] + P[3]
    molecules.append(source_molecule(f"{DOC}:mol:{n}", str(n), ev,
        stereo=stereo, geometry=geometry, representation=rep,
        identity_note=note, assignment=assignment, axial=axial, double_bonds=double))

# Table 3 inputs 1-10 are represented independently because the source table
# prints the product compound number but not a substrate compound number.
for i, xdesc in enumerate([
    "p-methylphenyl", "p-methoxyphenyl", "p-fluorophenyl", "p-chlorophenyl", "p-trifluoromethylphenyl",
    "m-chlorophenyl", "m-fluorophenyl", "3,5-difluorophenyl", "o-fluorophenyl", "2,6-difluorophenyl"
], 1):
    molecules.append(source_molecule(f"{DOC}:mol:table3_input_{i}", f"Table 3 entry {i} allene", P[3],
        name=f"aromatic allenyl diene substrate, X = {xdesc}", stereo="specified",
        geometry="partially_specified", representation="not_encoded",
        axial=[{"type": "allene axis", "configuration": "source-drawn enantiomer; R_a/S_a not named"}],
        assignment="Optically active aromatic allene is drawn; absolute descriptor is not printed.",
        identity_note="Table row defines the aryl substituent, but a definitive full SMILES was not asserted."))

common_specs = [
    ("aldehyde_R1", "R1CHO", None, "unresolved", "unresolved"),
    ("CuBr2", "CuBr2", "[Br-].[Br-].[Cu+2]", "achiral", "not_applicable"),
    ("dioxane", "dioxane", "O1CCOCC1", "achiral", "not_applicable"),
    ("MsCl", "methanesulfonyl chloride", "CS(=O)(=O)Cl", "achiral", "not_applicable"),
    ("Et3N", "triethylamine", "CCN(CC)CC", "achiral", "not_applicable"),
    ("CH2Cl2", "dichloromethane", "ClCCl", "achiral", "not_applicable"),
    ("PPh3", "triphenylphosphine", "P(c1ccccc1)(c1ccccc1)c1ccccc1", "achiral", "not_applicable"),
    ("DIAD", "diisopropyl azodicarboxylate", "CC(C)OC(=O)N=NC(=O)OC(C)C", "achiral", "unspecified"),
    ("THF", "tetrahydrofuran", "C1CCOC1", "achiral", "not_applicable"),
    ("NaH", "sodium hydride", "[H-].[Na+]", "achiral", "not_applicable"),
    ("toluene", "toluene", "Cc1ccccc1", "achiral", "not_applicable"),
]
for key, name, smiles, chiral, geom in common_specs:
    ev = P[1] + P[2] + P[3]
    molecules.append(source_molecule(f"{DOC}:mol:{key}", key, ev, name=name,
        stereo=chiral, geometry=geom,
        representation="partial" if key == "DIAD" else "not_applicable" if smiles else "unresolved",
        smiles=smiles, identity_note="Identity follows the printed reagent/solvent label; CAS was not reported."))

for key, name, role_stereo in [
    ("RhcodCl2", "[Rh(cod)Cl]2", "unresolved"),
    ("Rh_SS_Boz", "[Rh(nbd){(S,S)-Me-BozPHOS}]SbF6", "specified"),
    ("Rh_RR_Boz", "[Rh(nbd){(R,R)-Me-BozPHOS}]SbF6", "specified"),
    ("BINAPO", "BINAP(O)", "unspecified"), ("BINAP", "BINAP", "unspecified"),
    ("dppe", "dppe", "unresolved"), ("dppeO", "dppeO", "unresolved"),
    ("dppmO", "dppmO", "unresolved"), ("dpppO", "dpppO", "unresolved"),
    ("RhCODdppeOSbF6", "[Rh(COD)dppeO]SbF6", "unresolved"),
]:
    molecules.append(source_molecule(f"{DOC}:mol:{key}", key, P[2] + P[3], name=name,
        stereo=role_stereo, geometry="not_applicable", representation="not_encoded",
        assignment=("Source complex name explicitly specifies (S,S)." if key == "Rh_SS_Boz" else
                    "Source complex name explicitly specifies (R,R)." if key == "Rh_RR_Boz" else None),
        identity_note="Complex/ligand retained by the source abbreviation; no definitive SMILES or CAS was asserted."))

# Scheme 1 catalyst 2 is already a numbered molecule; ensure it can be used as a catalyst.
mol_by_id = {m["molecule_id"]: m for m in molecules}

reactions = []

# Scheme 1: four generic preparation arrows and one specific 10 + 11 -> 12 example.
r = reaction("scheme1_1_aldehyde_to_3", P[1], record_type="generic_scheme",
             location="page 1, Scheme 1, top-left arrow", name="prolinol-controlled copper-catalyzed asymmetric allene synthesis",
             transformation="generic asymmetric allene formation")
r["participants"] = [participant("s1_1", f"{DOC}:mol:1", "reactant", P[1]),
                     participant("s1_aldehyde", f"{DOC}:mol:aldehyde_R1", "reactant", P[1]),
                     participant("s1_cat2", f"{DOC}:mol:2", "catalyst", P[1]),
                     participant("s1_cubr2", f"{DOC}:mol:CuBr2", "catalyst", P[1]),
                     participant("s1_dioxane", f"{DOC}:mol:dioxane", "solvent", P[1]),
                     participant("s1_3", f"{DOC}:mol:3", "product", P[1], product_role="main product")]
set_solvent(r, P[1], "dioxane")
r["steps"] = [step("reaction", P[1], [p["participant_id"] for p in r["participants"]],
                   operation="catalytic reaction shown in Scheme 1",
                   temperature="alternative temperatures as printed")]
r["steps"][0]["temperature"] = f("70 or 130 °C", "alternative temperatures as printed", P[1])
r["steps"][0]["temperature"]["alternatives"] = [
    f("70 °C", q(70, "°C", "70 °C"), P[1]),
    f("130 °C", q(130, "°C", "130 °C"), P[1]),
]
r["outcomes"] = [outcome("yield_range", "yield", [f"{DOC}:mol:3"], P[1],
                         {"minimum": 50, "maximum": 70, "unit": "%", "raw_text": "50–70%"}, role="main product"),
                 outcome("ee_range", "selectivity", [f"{DOC}:mol:3"], P[1],
                         {"minimum": 93, "maximum": 99, "unit": "%", "raw_text": "93–99% ee"},
                         scope={"metric": "ee", "configuration": "unassigned"})]
r["review"]["unresolved"] = ["Generic R groups, catalyst loadings, time, yield basis, and ee method are not reported on this page."]
unresolved_classification(r, P[1], "Multicomponent allene-forming sequence is generic and not safely atom-mapped from the main-text drawing.")
reactions.append(r)

r = reaction("scheme1_3_to_4_mesylation", P[1], record_type="generic_scheme",
             location="page 1, Scheme 1, top-right arrow", name="mesylation of allenol 3", transformation="alcohol mesylation")
r["participants"] = [participant("s1_3", f"{DOC}:mol:3", "reactant", P[1]),
                     participant("s1_mscl", f"{DOC}:mol:MsCl", "reagent", P[1]),
                     participant("s1_et3n", f"{DOC}:mol:Et3N", "reagent", P[1]),
                     participant("s1_ch2cl2", f"{DOC}:mol:CH2Cl2", "solvent", P[1]),
                     participant("s1_4", f"{DOC}:mol:4", "product", P[1], product_role="main product")]
set_solvent(r, P[1], "CH2Cl2")
r["steps"] = [step("reaction", P[1], [p["participant_id"] for p in r["participants"]], operation="MsCl/Et3N treatment")]
r["review"]["unresolved"] = ["No quantities, temperature, time, yield, workup, or yield basis are printed in the main text."]
simple_classification(r, P[1],
    [{"change_id": "b1", "atoms": [{"element": "O"}, {"element": "S"}], "before_order": 0, "after_order": 1, "change_type": "formed"}],
    [{"category_number": 26, "code": "O-S", "supporting_change_ids": ["b1"]}],
    "O–S formation inferred from conversion of the drawn allenol to its mesylate; other reagent/coproduct changes unresolved.")
reactions.append(r)

r = reaction("scheme1_5_plus_6_to_7", P[1], record_type="generic_scheme",
             location="page 1, Scheme 1, middle arrow", name="Mitsunobu substitution", transformation="Mitsunobu N-alkylation")
r["participants"] = [participant("s1_5", f"{DOC}:mol:5", "reactant", P[1]), participant("s1_6", f"{DOC}:mol:6", "reactant", P[1]),
                     participant("s1_pph3", f"{DOC}:mol:PPh3", "reagent", P[1]), participant("s1_diad", f"{DOC}:mol:DIAD", "reagent", P[1]),
                     participant("s1_thf", f"{DOC}:mol:THF", "solvent", P[1]), participant("s1_7", f"{DOC}:mol:7", "product", P[1], product_role="main product")]
set_solvent(r, P[1], "THF")
r["steps"] = [step("reaction", P[1], [p["participant_id"] for p in r["participants"]], operation="Mitsunobu substitution", temperature="rt", time=q(20, "h", "20 h"))]
r["review"]["unresolved"] = ["No quantities, yield, workup, purification, or yield basis are printed in the main text."]
simple_classification(r, P[1],
    [{"change_id": "b1", "atoms": [{"element": "C"}, {"element": "N"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
     {"change_id": "b2", "atoms": [{"element": "C"}, {"element": "O"}], "before_order": 1, "after_order": 0, "change_type": "broken"}],
    [{"category_number": 5, "code": "C-N", "supporting_change_ids": ["b1"]}],
    "C–N formation/C–O cleavage inferred from the inspected Mitsunobu arrow; exact atom maps and byproducts unresolved.")
reactions.append(r)

r = reaction("scheme1_8_plus_4_to_9", P[1], record_type="generic_scheme",
             location="page 1, Scheme 1, lower-middle arrow", name="alkylation of 1,3-diesterdiene 8 with mesylate 4", transformation="C-alkylation")
r["participants"] = [participant("s1_8", f"{DOC}:mol:8", "reactant", P[1]), participant("s1_4", f"{DOC}:mol:4", "reactant", P[1]),
                     participant("s1_nah", f"{DOC}:mol:NaH", "reagent", P[1]), participant("s1_thf", f"{DOC}:mol:THF", "solvent", P[1]),
                     participant("s1_9", f"{DOC}:mol:9", "product", P[1], product_role="main product")]
set_solvent(r, P[1], "THF")
r["steps"] = [step("deprotonation", P[1], ["s1_8", "s1_nah", "s1_thf"], operation="treat 8 with NaH", temperature=q(0, "°C", "0 °C"), time=q(30, "min", "30 min")),
              step("mesylate_addition", P[1], ["s1_8", "s1_4", "s1_9"], operation="add mesylate allene 4")]
r["review"]["unresolved"] = ["The second-stage temperature/time, all quantities, yield, workup, purification, and yield basis are not printed."]
simple_classification(r, P[1],
    [{"change_id": "b1", "atoms": [{"element": "C"}, {"element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
     {"change_id": "b2", "atoms": [{"element": "C"}, {"element": "O"}], "before_order": 1, "after_order": 0, "change_type": "broken"}],
    [{"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1"]}],
    "C–C formation and leaving-group C–O cleavage inferred from the inspected alkylation arrow; exact maps and proton transfers unresolved.")
reactions.append(r)

r = reaction("scheme1_10_plus_11_to_12", P[1], categories=["other"],
             location="page 1, Scheme 1, bottom arrow", name="direct conversion of dienyne 10 to allenyl diene 12", transformation="asymmetric allene formation")
r["participants"] = [participant("s1_10", f"{DOC}:mol:10", "reactant", P[1]), participant("s1_11", f"{DOC}:mol:11", "reactant", P[1]),
                     participant("s1_cat2", f"{DOC}:mol:2", "catalyst", P[1]), participant("s1_cubr2", f"{DOC}:mol:CuBr2", "catalyst", P[1]),
                     participant("s1_dioxane", f"{DOC}:mol:dioxane", "solvent", P[1]), participant("s1_12", f"{DOC}:mol:12", "product", P[1], product_role="main product")]
set_solvent(r, P[1], "dioxane")
r["steps"] = [step("reaction", P[1], [p["participant_id"] for p in r["participants"]], operation="catalytic two-step sequence shown in Scheme 1", temperature=q(130, "°C", "130 °C"))]
r["outcomes"] = [outcome("overall_yield", "yield", [f"{DOC}:mol:12"], P[1], q(32, "%", "32% over 2 steps"), role="main product", scope={"span": "10 + 11 to 12", "stages": 2}),
                 outcome("product_er", "selectivity", [f"{DOC}:mol:12"], P[1], ratio_value("99:1"), scope={"metric": "er", "order": "major:minor; configurations unassigned"})]
r["review"]["unresolved"] = ["Catalyst loadings, scale, time, individual step boundary/yields, yield basis, er method, workup, and purification are absent from the main text."]
unresolved_classification(r, P[1], "The two-step allene-forming sequence lacks an explicit intermediate and cannot be safely atom-mapped from the main text alone.")
reactions.append(r)


def cyclo_record(rid, ev, substrate_mid, product_mid, *, location, categories,
                 catalyst_mid, catalyst_pct, ligand_mid=None, ligand_pct=None,
                 solvent="toluene", solvent_mid="toluene", solvent_volume=None,
                 temp=None, duration=None, yield_value=None, yield_basis=None,
                 product_er=None, product_er_basis=None, substrate_er=None,
                 trace=False, unstable=False, sealed=False, extra_evidence=None):
    ev = list(dict.fromkeys(ev + (extra_evidence or [])))
    r = reaction(rid, ev, categories=categories, location=location,
                 name="rhodium-catalyzed intramolecular [4 + 2] cycloisomerization",
                 transformation="intramolecular allene–diene [4+2] cycloisomerization")
    s_mid = f"{DOC}:mol:{substrate_mid}" if isinstance(substrate_mid, int) else substrate_mid
    p_mid = f"{DOC}:mol:{product_mid}" if isinstance(product_mid, int) else product_mid
    parts = [participant("substrate", s_mid, "reactant", ev),
             participant("product", p_mid, "product", ev, product_role="main cycloisomerization product"),
             participant("catalyst", f"{DOC}:mol:{catalyst_mid}", "catalyst", ev,
                         molpct=q(catalyst_pct, "mol%", f"{catalyst_pct} mol%"), ref=s_mid)]
    if ligand_mid:
        parts.append(participant("ligand", f"{DOC}:mol:{ligand_mid}", "ligand", ev,
                                 molpct=q(ligand_pct, "mol%", f"{ligand_pct} mol%"), ref=s_mid))
    parts.append(participant("solvent", f"{DOC}:mol:{solvent_mid}", "solvent", ev,
                             volume=solvent_volume))
    r["participants"] = parts
    if isinstance(substrate_mid, int) and substrate_mid == 12 and ("Table 1" in location or "Table 2" in location):
        r["participants"][0]["charges"]["amount"] = f("0.1 mmol", q(0.1, "mmol", "0.1 mmol"), ev)
        r["reaction_scale"] = f("12 (0.1 mmol)", {"reference_molecule_id": s_mid, "amount": q(0.1, "mmol", "0.1 mmol")}, ev)
    set_solvent(r, ev, solvent, solvent_volume)
    pids = [p["participant_id"] for p in parts]
    r["steps"] = [step("cycloisomerization", ev, pids, operation="cycloisomerization under displayed conditions", temperature=temp, time=duration)]
    if sealed:
        r["special_conditions"]["reactor"] = f("sealed tube", "sealed tube", ev)
    outs = []
    if yield_value is not None:
        outs.append(outcome("product_yield", "yield", [p_mid], ev, yield_value,
                            role="main product", yield_type="isolated" if yield_basis == "isolated yield" else None,
                            basis=yield_basis))
    if trace:
        outs.append(outcome("trace_product", "trace", [p_mid], ev, "trace", role="main product"))
    if unstable:
        outs.append(outcome("product_stability", "other", [p_mid], ev, "product unstable", role="main product",
                            scope={"observation": "product unstable; no isolated yield printed"}))
    if substrate_er:
        outs.append(outcome("starting_allene_er", "selectivity", [s_mid], ev, ratio_value(substrate_er),
                            basis=None,
                            scope={"metric": "starting_material_er", "order": "major:minor; configurations unassigned"}))
    if product_er:
        outs.append(outcome("product_er", "selectivity", [p_mid], ev, ratio_value(product_er),
                            basis=product_er_basis,
                            scope={"metric": "product_er", "order": "major:minor; configurations unassigned"}))
    if "Table 1" in location:
        outs.append(outcome("diastereomer_observation", "selectivity", [p_mid], ev,
                            "Within the limits of detection only one diastereomer was observed.",
                            scope={"metric": "diastereomer observation", "ratio": "not quantified"}))
    r["outcomes"] = outs
    unresolved = ["Substrate/product definitive SMILES, exact atom maps, workup, and purification are unavailable from the supplied main text."]
    if yield_value is not None and not yield_basis:
        unresolved.append("Yield measurement/isolation basis is not printed for this scope value.")
    if product_er and not product_er_basis:
        unresolved.append("Product er measurement method is not printed for this scope value.")
    if substrate_er:
        unresolved.append("Starting-material er measurement method is not separately printed for this value.")
    r["review"]["unresolved"] = unresolved
    cyclo_classification(r, ev)
    return r


# Table 1, eight catalyst-screening rows.
table1 = [
    (1, "RhcodCl2", 2.5, None, None, 24, 7, "50:50"),
    (2, "Rh_SS_Boz", 5, None, None, 12, 37, "54:46"),
    (3, "Rh_SS_Boz", 5, None, None, 24, 54, "54:46"),
    (4, "Rh_RR_Boz", 5, None, None, 24, 57, "85:15"),
    (5, "RhcodCl2", 2.5, "BINAPO", 5, 24, 43, "51:49"),
    (6, "RhcodCl2", 2.5, "BINAP", 5, 24, 18, "95:5"),
    (7, "RhcodCl2", 2.5, "dppe", 5, 24, 57, "92:8"),
    (8, "RhcodCl2", 2.5, "dppeO", 5, 24, 83, "97:3"),
]
for entry, cat, catpct, lig, ligpct, hours, yld, er in table1:
    r = cyclo_record(f"table1_entry_{entry}", P[1] + P[2], 12, 13,
        location=f"page 2, Table 1, entry {entry}", categories=["optimization"],
        catalyst_mid=cat, catalyst_pct=catpct, ligand_mid=lig, ligand_pct=ligpct,
        solvent_volume=q(1.0, "mL", "1.0 mL"), temp=q(80, "°C", "80 °C"), duration=q(hours, "h", f"{hours} h"),
        yield_value=q(yld, "%", f"{yld}%"), yield_basis="isolated yield",
        product_er=er, product_er_basis="chiral HPLC", substrate_er="97:3")
    add_note(r, P[1] + P[2], "Source-label conflict retained: page-1 prose says the catalyst probe used 'substrate 13', whereas the page-2 Table 1 reaction drawing and footnote explicitly show substrate 12 giving product 13.")
    r["review"]["unresolved"].append("Page-1 prose calls the screening material substrate 13; Table 1 drawing/footnote identify it as substrate 12. The structured participants follow the explicit table drawing and preserve the prose conflict in notes.")
    reactions.append(r)

# Table 2, six independent ligand/temperature/time rows. Entry 2 is also the
# Scheme 2 12 -> 13 scope point and therefore receives both evidence roles.
table2 = [
    (1, 12, "dppeO", 80, 82, "97:3", False),
    (2, 6, "dppeO", 80, 83, "97:3", False),
    (3, 6, "dppmO", 80, 88, "97:3", False),
    (4, 6, "dpppO", 80, 67, "96:4", False),
    (5, 6, "dppeO", 55, 30, "97:3", False),
    (6, 6, "dppeO", "rt", None, None, True),
]
for entry, hours, lig, temp, yld, er, is_trace in table2:
    extra = P[2] if entry == 2 else []
    r = cyclo_record(f"table2_entry_{entry}", P[2], 12, 13,
        location=f"page 2, Table 2, entry {entry}" + ("; same 6 h experiment also displayed as Scheme 2 12→13" if entry == 2 else ""),
        categories=["optimization"] + (["substrate_scope"] if entry == 2 else []),
        catalyst_mid="RhcodCl2", catalyst_pct=2.5, ligand_mid=lig, ligand_pct=5,
        solvent_volume=q(1.0, "mL", "1.0 mL"), temp=q(temp, "°C", f"{temp} °C") if isinstance(temp, (int, float)) else "rt",
        duration=q(hours, "h", f"{hours} h"), yield_value=q(yld, "%", f"{yld}%") if yld is not None else None,
        yield_basis="isolated yield" if yld is not None else None, product_er=er,
        product_er_basis="chiral HPLC" if er else None, substrate_er="97:3", trace=is_trace,
        extra_evidence=extra)
    if entry == 2:
        add_note(r, P[2], "The 6 h, 80 °C, dppeO result (83%, er 97:3) is displayed both in Table 2 entry 2 and Scheme 2; it is one experiment, not two records.")
        add_note(r, P[2], "Yield conflict retained: the surrounding prose states 6 h, 82%, er 97:3, while Table 2 entry 2 and Scheme 2 state 6 h, 83%, er 97:3.")
        conflict = f(["prose: 82%", "Table 2 entry 2 / Scheme 2: 83%"], None, P[2], status="conflict")
        conflict["alternatives"] = [
            f("82% in surrounding prose", q(82, "%", "82%"), P[2]),
            f("83% in Table 2 entry 2 and Scheme 2", q(83, "%", "83%"), P[2]),
        ]
        r["outcomes"][0]["value"] = conflict
        r["review"]["unresolved"].append("The 6 h isolated yield conflicts between prose (82%) and Table 2/Scheme 2 (83%); no value was selected over the other.")
    reactions.append(r)

# Scheme 2 scope points except 12 -> 13, already represented by Table 2 entry 2.
scheme2 = [
    (14, 15, "97:3", 96, "90:10", False), (16, 17, "99:1", 87, "96:4", False),
    (18, 19, "96:4", 85, "93:7", False), (20, 21, "98:2", 81, "87:13", False),
    (22, 23, "99:1", 85, "95:5", False), (24, 25, "96:4", 75, "93:7", False),
    (26, 27, "99:1", 68, "93:7", False), (28, 29, "99:1", 67, "82:18", False),
    (30, 31, "99:1", None, None, True), (32, 33, "99:1", 80, "98:2", False),
    (34, 35, "97:3", 82, "96:4", False), (36, 37, "98:2", 83, "91:9", False),
    (38, 39, "99:1", 75, "96:4", False), (40, 41, "99:1", 75, "99:1", False),
]
for s, p, ser, yld, per, unstable in scheme2:
    reactions.append(cyclo_record(f"scheme2_{s}_to_{p}", P[2] + (P[3] if s == 30 else []), s, p,
        location=f"page 2, Scheme 2, {s}→{p}", categories=["substrate_scope"],
        catalyst_mid="RhcodCl2", catalyst_pct=2.5, ligand_mid="dppeO", ligand_pct=5,
        temp=q(80, "°C", "80 °C"), duration=q(6, "h", "6 h"),
        yield_value=q(yld, "%", f"{yld}%") if yld is not None else None,
        product_er=per, substrate_er=ser, unstable=unstable))

# Scheme 3 scope.
for s, p, ser, yld, per in [
    (42, 43, "99:1", 71, ">99:1"), (44, 45, "96:4", 93, "97:3"),
    (46, 47, "97:3", 78, "98:2"), (48, 49, "96:4", 81, "95:5")
]:
    reactions.append(cyclo_record(f"scheme3_{s}_to_{p}", P[3], s, p,
        location=f"page 3, Scheme 3, {s}→{p}", categories=["substrate_scope"],
        catalyst_mid="RhcodCl2", catalyst_pct=2.5, ligand_mid="dppeO", ligand_pct=5,
        temp=q(80, "°C", "80 °C"), duration=q(6, "h", "6 h"),
        yield_value=q(yld, "%", f"{yld}%"), product_er=per, substrate_er=ser))

# Table 3 aromatic substrates; substrate labels are deliberately row-local and
# not inferred as even compound numbers.
table3 = [
    (1, 51, "99:1", "90:10", 61), (2, 53, "90:10", "85:15", 61),
    (3, 55, "98:2", "85:15", 53), (4, 57, "98:2", "84:16", 46),
    (5, 59, "98:2", "87:13", 67), (6, 61, "97:3", "90:10", 51),
    (7, 63, "99:1", "92:8", 53), (8, 65, "97:3", "92:8", 58),
    (9, 67, "98:2", "95:5", 53), (10, 69, "96:4", "90:10", 59),
    (11, 29, "99:1", "99:1", 51),
]
for entry, product, ser, per, yld in table3:
    substrate = f"{DOC}:mol:28" if entry == 11 else f"{DOC}:mol:table3_input_{entry}"
    extra = P[2] if entry == 11 else []
    reactions.append(cyclo_record(f"table3_entry_{entry}_product_{product}", P[3], substrate, product,
        location=f"page 3, Table 3, entry {entry}, product {product}", categories=["substrate_scope"],
        catalyst_mid="RhCODdppeOSbF6", catalyst_pct=5,
        solvent="CH2Cl2", solvent_mid="CH2Cl2", temp=q(40, "°C", "40 °C"), duration=q(18, "h", "18 h"),
        yield_value=q(yld, "%", f"{yld}%"), product_er=per, substrate_er=ser, sealed=True,
        extra_evidence=extra))

bundle["molecules"] = molecules
bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": f"{DOC}:metadata", "kind": "bibliographic_metadata", "evidence_ids": P[1],
     "title": "Diastereoselective Rhodium Catalyzed [4 + 2] Cycloisomerization of Allenes",
     "authors": ["Jun Li", "Scott R. Gilbertson"], "journal": "Organic Letters", "year": 2021,
     "volume": 23, "pages": "2911–2914", "doi": "10.1021/acs.orglett.1c00554"},
    {"context_id": f"{DOC}:graphical_abstract", "kind": "summary_not_extra_experiment", "evidence_ids": P[1],
     "note": "Graphical abstract summarizes the study and is not counted as an additional experiment."},
    {"context_id": f"{DOC}:prior_work", "kind": "literature_context", "evidence_ids": P[1] + P[4],
     "note": "Introduction and references describe prior Ni/Rh/Au/Pt/Ru allene chemistry; cited work is not promoted to current-paper experiments."},
    {"context_id": f"{DOC}:aromatic_isomerization", "kind": "author_observation", "evidence_ids": P[3],
     "note": "Authors state that aromatic allene substrates isomerize over time under reaction conditions; no row-level time-course values or complete standalone experiment are printed."},
    {"context_id": f"{DOC}:catalyst_xray", "kind": "characterization", "evidence_ids": P[3],
     "note": "Figure 1 is the X-ray structure of [Rh(COD)dppeO]SbF6; characterization, not a synthetic experiment."},
    {"context_id": f"{DOC}:product29_xray", "kind": "characterization", "evidence_ids": P[3],
     "note": "Figure 2 is the X-ray structure of product 29 and supports cis ring fusion/exocyclic alkene assignment; no separate synthesis record was created."},
    {"context_id": f"{DOC}:four_atom_tether_limitation", "kind": "reaction_scope_limitation", "evidence_ids": P[3],
     "note": "Authors report that 4-atom tethers intended to form [4.4.0] systems do not react; no specific substrate, conditions, or row is supplied, so no fabricated control record was created."},
    {"context_id": f"{DOC}:si_notice", "kind": "supporting_information_notice", "evidence_ids": P[4],
     "note": "Page 4 explicitly states that SI contains procedures, characterization, and NMR spectra, but no SI file was supplied or associated in this workspace."},
    {"context_id": f"{DOC}:ccdc_notice", "kind": "crystallographic_data_notice", "evidence_ids": P[4],
     "note": "CCDC 2062291 and 2062293 are cited; no CIF was supplied or reviewed."},
    {"context_id": f"{DOC}:page4_bibliography", "kind": "bibliography", "evidence_ids": P[4],
     "note": "Page 4 references and author/acknowledgment metadata contain no new current-paper synthetic experiment."},
]

write_json(RUN / "submitted.json", bundle)

record_ids = [r["reaction_id"] for r in reactions]
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title, abstract, graphical abstract, introduction", "finding": "bibliographic metadata, study summary, and cited prior work; graphical abstract is not an extra experiment", "evidence_ids": P[1]},
            {"region": "Scheme 1 top row", "finding": "generic 1 + R1CHO -> 3 and 3 -> 4 preparation arrows", "evidence_ids": P[1]},
            {"region": "Scheme 1 middle/lower rows", "finding": "generic Mitsunobu 5 + 6 -> 7, alkylation 8 + 4 -> 9, and specific two-step 10 + 11 -> 12 (32%, er 99:1)", "evidence_ids": P[1]},
            {"region": "lower-column catalyst-screen introduction", "finding": "prose context for Table 1; no additional row beyond page-2 table", "evidence_ids": P[1]},
        ], "record_ids": [rid for rid in record_ids if any(x in rid for x in ["scheme1_", "table1_"])],
         "unresolved": ["No associated SI; Scheme 1 omits most loadings, scales, times, workup, purification, and measurement bases.", "Generic R-group structures and allene absolute descriptors remain unresolved."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Table 1 catalyst screening", "finding": "eight independent optimization rows, with isolated yields, product er by chiral HPLC, and one-diastereomer footnote", "evidence_ids": P[2]},
            {"region": "Table 2 ligand/temperature screen", "finding": "six independent optimization rows including room-temperature trace result", "evidence_ids": P[2]},
            {"region": "Scheme 2 substrate scope", "finding": "fifteen displayed scope arrows; 12->13 at 6 h duplicates Table 2 entry 2 and is merged, leaving fourteen additional records", "evidence_ids": P[2]},
            {"region": "prose surrounding Tables 1-2 and Scheme 2", "finding": "discussion cross-checks ligand, time, temperature, substrate classes, and product instability; no extra numbered experiment", "evidence_ids": P[2]},
        ], "record_ids": [rid for rid in record_ids if any(x in rid for x in ["table1_", "table2_", "scheme2_"])],
         "unresolved": ["Complex scope structures and exact atom mapping are not safely encoded.", "Scheme 2 yield/er measurement bases are not footnoted on the supplied page; only Table 1/2 methods were assigned where explicit."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Scheme 3 reaction with other tethers", "finding": "four substrate-scope arrows with yields and input/product er", "evidence_ids": P[3]},
            {"region": "Table 3 aromatic allene substrates", "finding": "eleven independent scope rows under discrete [Rh(COD)dppeO]SbF6/CH2Cl2/40 °C/18 h conditions", "evidence_ids": P[3]},
            {"region": "Figure 1 and catalyst-complex prose", "finding": "X-ray characterization of catalyst complex and equivalence of discrete/in-situ catalyst systems; no separate synthesis record", "evidence_ids": P[3]},
            {"region": "Figure 2 and stereochemical prose", "finding": "X-ray characterization of product 29 supports cis ring fusion and drawn exocyclic alkene geometry", "evidence_ids": P[3]},
            {"region": "scope limitations/conclusion", "finding": "aromatic-substrate isomerization observation and generic 4-atom-tether failure; insufficient detail for standalone mechanistic/control records", "evidence_ids": P[3]},
        ], "record_ids": [rid for rid in record_ids if any(x in rid for x in ["scheme3_", "table3_"]) or rid.endswith("scheme2_30_to_31")],
         "unresolved": ["Table 3 does not explicitly state yield or er measurement methods.", "Aromatic-substrate time-course values and specific 4-atom-tether identity/conditions are absent.", "No CIF was supplied; X-ray assignments rely on the paper's figure/prose only."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "Associated Content / Supporting Information", "finding": "SI availability notice only; SI was not supplied or associated", "evidence_ids": P[4]},
            {"region": "Accession Codes", "finding": "CCDC 2062291 and 2062293 notice; no CIF reviewed", "evidence_ids": P[4]},
            {"region": "author information, acknowledgments, references", "finding": "metadata and cited prior work; no current-paper reaction arrow, table row, spectrum, or new experiment", "evidence_ids": P[4]},
        ], "record_ids": [],
         "unresolved": ["Supporting Information and CIF files were not provided; their experimental details were not inferred."]},
    ],
}
write_json(RUN / "coverage/pages_0001_0004.json", coverage)

inventory = {
    "document_id": DOC,
    "scope": "original PDF pages 1-4",
    "source_title": "Diastereoselective Rhodium Catalyzed [4 + 2] Cycloisomerization of Allenes",
    "source_sha256": dataset["documents"][0]["sha256"],
    "visual_review": {"full_pages_opened": [1, 2, 3, 4],
                      "source_lineage_crops": [p.name for p in sorted((RUN / "raw_model").glob("page*.png"))]},
    "units": [
        {"region": "page1 Scheme1", "records": [rid for rid in record_ids if "scheme1_" in rid]},
        {"region": "page2 Table1", "records": [rid for rid in record_ids if "table1_" in rid]},
        {"region": "page2 Table2", "records": [rid for rid in record_ids if "table2_" in rid]},
        {"region": "page2 Scheme2", "records": [rid for rid in record_ids if "scheme2_" in rid], "deduplication": "12->13/6 h is Table2 entry2"},
        {"region": "page3 Scheme3", "records": [rid for rid in record_ids if "scheme3_" in rid]},
        {"region": "page3 Table3", "records": [rid for rid in record_ids if "table3_" in rid]},
        {"region": "page3 Figures1-2", "disposition": "source_context / stereochemical evidence, not extra syntheses"},
        {"region": "page4", "disposition": "SI/CCDC/metadata/bibliography; no new experiment"},
    ],
    "reaction_count": len(reactions),
    "molecule_count": len(molecules),
    "machine_status": "needs_review",
}
write_json(RUN / "raw_model/extraction_inventory.json", inventory)

tasks_path = RUN / "tasks.json"
tasks = json.loads(tasks_path.read_text())
for task in tasks:
    if task["document_id"] == DOC:
        task["state"] = "attempted"
        task["reason"] = "Full original page image inspected and mapped in coverage/pages_0001_0004.json; attempted does not mean semantic accuracy verified."
write_json(tasks_path, tasks)

print(json.dumps({"molecules": len(molecules), "reactions": len(reactions),
                  "source_context": len(bundle["source_context"])}, ensure_ascii=False))
