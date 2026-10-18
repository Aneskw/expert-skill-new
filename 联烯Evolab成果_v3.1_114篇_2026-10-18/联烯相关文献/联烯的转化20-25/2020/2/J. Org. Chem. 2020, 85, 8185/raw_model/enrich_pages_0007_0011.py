#!/usr/bin/env python3
"""Enrich the v3.1 bundle from original PDF images/pages 7--11.

Pages 7--9 continue product characterization and give three full procedures.
Pages 10--11 are references.  Existing cross-page experiments are updated in
place; characterization and bibliography are retained as source context rather
than fabricated reaction records.  All imported records remain needs_review.
"""

from __future__ import annotations

import copy
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills" / "chemical-reaction-extraction" / "scripts"
sys.path.insert(0, str(SCRIPTS))

import standard_records as st


DOC = "doc_4e4afb29c8afe53f"
CURRENT = ROOT / "runs" / "current"
RAW = CURRENT / "raw_model"
# Immutable same-document pages-1--6 construction bundle is the baseline.
# Using it makes this enrichment reproducible after the live dataset has
# already been merged once, without re-appending outcomes/participants.
SOURCE = RAW / "extracted_bundle_pages_0001_0006.json"
SUBMITTED = CURRENT / "submitted.json"


def ev(page: int) -> list[str]:
    return [f"{DOC}:page{page}:image", f"{DOC}:page{page}:text"]


def q(number, unit, raw_text=None):
    value = {"number": number, "unit": unit}
    if raw_text is not None:
        value["raw_text"] = raw_text
    return value


def fact(raw=None, normalized=None, evidence=(), method="original_extraction",
         status=None, derivation=None, alternatives=()):
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   derivation=derivation, alternatives=alternatives)


def missing(evidence=(), status="unreported"):
    return fact(evidence=evidence, status=status)


def append_unique(seq: list, values: list):
    for value in values:
        if value not in seq:
            seq.append(value)


bundle = json.loads(SOURCE.read_text())
# The earlier bounded same-document bundle contains only the image evidence that was in scope
# for pages 1--6.  The immutable prepared dataset now contains the full 1--11
# evidence registry; use it verbatim, exactly as a fresh `standard init` does.
dataset = json.loads((CURRENT / "dataset.json").read_text())
bundle["documents"] = copy.deepcopy(dataset["documents"])
molecules = {m["molecule_id"]: m for m in bundle["molecules"]}
reactions = {r["reaction_id"]: r for r in bundle["reactions"]}


def mid(label: str) -> str:
    return f"{DOC}:mol_{label}"


def rid(suffix: str) -> str:
    return f"{DOC}:{suffix}"


def update_fact_evidence(obj: dict, evidence: list[str]):
    append_unique(obj["evidence_ids"], evidence)


def add_molecule_occurrence_evidence(label: str, page: int, *, name_visible=True):
    m = molecules[mid(label)]
    evidence = ev(page)
    occurrence = m["occurrences"][0]
    append_unique(occurrence["evidence_ids"], evidence)
    update_fact_evidence(occurrence["compound_label"], evidence)
    if name_visible:
        for key in ("original", "english", "aliases"):
            update_fact_evidence(m["names"][key], evidence)
        update_fact_evidence(m["identity_note"], evidence)


def add_reaction_page(reaction: dict, page: int):
    append_unique(reaction["evidence_ids"], ev(page))


def participant_by_id(reaction: dict, participant_id: str) -> dict:
    return next(p for p in reaction["participants"] if p["participant_id"] == participant_id)


def outcome_by_id(reaction: dict, outcome_id: str) -> dict:
    return next(o for o in reaction["outcomes"] if o["outcome_id"] == outcome_id)


def new_participant(pid: str, molecule_id: str, role: str, evidence: list[str],
                    stage="reaction", *, amount=None, mass=None, volume=None,
                    equivalents=None, mol_percent=None, eq_ref=None,
                    product_role=None) -> dict:
    p = st.new_participant(pid, molecule_id, role, evidence, stage=stage)
    p["product_role"] = (fact(product_role, product_role, evidence)
                         if product_role else missing(evidence))
    for key, value in (("amount", amount), ("mass", mass), ("volume", volume),
                       ("equivalents", equivalents), ("mol_percent", mol_percent)):
        p["charges"][key] = fact(value, value, evidence) if value is not None else missing(evidence)
    p["equivalence_reference"] = (
        fact(eq_ref, eq_ref, evidence, method="structure_derivation",
             derivation="Loading/equivalent denominator is the explicitly charged reference substrate.")
        if eq_ref else missing(evidence)
    )
    return p


def new_step(step_id: str, stage: str, pids: list[str], operation: str,
             evidence: list[str], temperature=None, time=None, atmosphere=None) -> dict:
    s = st.new_step(step_id, evidence, stage=stage)
    s["participant_ids"] = pids
    s["operation"] = fact(operation, operation, evidence)
    s["temperature"] = fact(temperature, temperature, evidence) if temperature is not None else missing(evidence)
    s["time"] = fact(time, time, evidence) if time is not None else missing(evidence)
    s["atmosphere"] = fact(atmosphere, atmosphere, evidence) if atmosphere is not None else missing(evidence)
    s["pressure"] = missing(evidence)
    return s


def product_amount_outcome(oid: str, product_label: str, amount: dict,
                           evidence: list[str], basis: str) -> dict:
    o = st.new_outcome(oid, "product_amount", [mid(product_label)], evidence)
    o["product_role"] = fact("target", "target", evidence)
    o["value"] = fact(amount, amount, evidence)
    o["yield_type"] = missing(evidence, status="not_applicable")
    o["measurement_basis"] = fact(basis, basis, evidence)
    o["scope"] = missing(evidence)
    return o


def set_isolated_yield(outcome: dict, evidence: list[str], wording: str):
    combined = list(dict.fromkeys(evidence + ev(4)))
    outcome["yield_type"] = fact(
        "product paragraph yield after the general flash-chromatography procedure",
        "isolated", combined, method="inherited_general_condition")
    outcome["measurement_basis"] = fact(
        wording, "isolated product after chromatographic purification", combined,
        method="inherited_general_condition")


# Product paragraphs on pages 7 and 8 add reported isolated masses.  They do
# not create new reactions: these are the same Table 2/Table 3 rows.
scope_page7 = {
    "x": 82, "y": 117, "z": 142, "aa": 121, "ab": 126,
    "ac": 112, "ad": 42, "ae": 65, "af": 141, "ag": 127,
}
scope_page8 = {"ah": 139, "ai": 104.5}
table3_page8 = {"aj": 149, "ak": 159, "al": 154, "am": 83,
                "an": 16, "ap": 106}

# 3w characterization begins on page 6 and merely continues on page 7.
r3w = reactions[rid("table2_scope_3w")]
add_reaction_page(r3w, 7)
add_molecule_occurrence_evidence("3w", 7, name_visible=False)
r3w["notes"].append(fact(
    "The 3w product heading, yield, and mass are on page 6; page 7 only continues its characterization.",
    "Page 7 supplies characterization continuation only, not a second synthesis experiment.", ev(7)))

for letter, mass in scope_page7.items():
    reaction = reactions[rid(f"table2_scope_3{letter}")]
    add_reaction_page(reaction, 7)
    add_molecule_occurrence_evidence(f"3{letter}", 7)
    update_fact_evidence(reaction["location"], ev(7))
    reaction["location"]["raw_value"] = f"Table 2 scope row 3{letter}; product paragraph on PDF page 7"
    reaction["location"]["normalized_value"] = f"PDF page 2/3 Table 2 and page 7 product 3{letter}"
    y = outcome_by_id(reaction, "yield")
    update_fact_evidence(y["value"], ev(7))
    set_isolated_yield(y, ev(7), f"3{letter} product paragraph reports yield and mass; general procedure states flash chromatography")
    reaction["outcomes"].append(product_amount_outcome(
        "product_amount", f"3{letter}", q(mass, "mg"), ev(7),
        "reported product mass in the experimental product paragraph"))
    reaction["review"]["unresolved"] = [
        x for x in reaction["review"]["unresolved"]
        if "Yield measurement basis" not in x
    ]

# The 3ag heading/yield/mass are on page 7 and its spectrum continues on page 8.
r3ag = reactions[rid("table2_scope_3ag")]
add_reaction_page(r3ag, 8)
add_molecule_occurrence_evidence("3ag", 8, name_visible=False)
r3ag["notes"].append(fact(
    "The 3ag product heading, yield, and mass are on page 7; page 8 only continues its characterization.",
    "Page 8 supplies characterization continuation only, not a second synthesis experiment.", ev(8)))

for letter, mass in scope_page8.items():
    reaction = reactions[rid(f"table2_scope_3{letter}")]
    add_reaction_page(reaction, 8)
    add_molecule_occurrence_evidence(f"3{letter}", 8)
    update_fact_evidence(reaction["location"], ev(8))
    reaction["location"]["raw_value"] = f"Table 2 scope row 3{letter}; product paragraph on PDF page 8"
    reaction["location"]["normalized_value"] = f"PDF pages 2-3 Table 2 and page 8 product 3{letter}"
    y = outcome_by_id(reaction, "yield")
    update_fact_evidence(y["value"], ev(8))
    set_isolated_yield(y, ev(8), f"3{letter} product paragraph reports yield and mass; general procedure states flash chromatography")
    reaction["outcomes"].append(product_amount_outcome(
        "product_amount", f"3{letter}", q(mass, "mg"), ev(8),
        "reported product mass in the experimental product paragraph"))
    reaction["review"]["unresolved"] = [
        x for x in reaction["review"]["unresolved"]
        if "Yield measurement basis" not in x
    ]

for suffix, mass in table3_page8.items():
    reaction = reactions[rid(f"table3_scope_3{suffix}")]
    add_reaction_page(reaction, 8)
    add_molecule_occurrence_evidence(f"3{suffix}", 8)
    update_fact_evidence(reaction["location"], ev(8))
    reaction["location"]["raw_value"] = f"Table 3 row 3{suffix}; product paragraph on PDF page 8"
    reaction["location"]["normalized_value"] = f"PDF page 3 Table 3 and page 8 product 3{suffix}"
    y = outcome_by_id(reaction, "yield")
    update_fact_evidence(y["value"], ev(8))
    set_isolated_yield(y, ev(8), f"3{suffix} product paragraph reports yield and mass; general procedure states flash chromatography")
    reaction["outcomes"].append(product_amount_outcome(
        "product_amount", f"3{suffix}", q(mass, "mg"), ev(8),
        "reported product mass in the experimental product paragraph"))
    reaction["review"]["unresolved"] = [
        x for x in reaction["review"]["unresolved"]
        if "Yield measurement basis" not in x
    ]


# Correct a ring-number transcription error in the prior 3aa candidate.  The
# corrected candidate follows the explicit page-7 product name; it remains an
# inferred, needs-review graph rather than a source-certified structure.
m3aa = molecules[mid("3aa")]
smiles_3aa = "CC(=O)CC1(P(=O)(c2ccccc2)c2ccccc2)CCCCC1"
diag_3aa = st.stereo_inspect(smiles_3aa)
e3aa = list(dict.fromkeys(ev(2) + ev(7)))
m3aa["smiles"] = fact(
    "1-(1-diphenylphosphinylcyclohexyl)propan-2-one candidate",
    smiles_3aa, e3aa, method="structure_derivation",
    derivation="Corrected candidate derived from the page-7 product name; distinct ring digits avoid the prior parser-induced fused graph.")
m3aa["stereochemistry"]["atom_centers"] = fact(
    "standard stereo diagnostic", diag_3aa["atom_centers"], e3aa,
    method="structure_derivation", derivation="Calculated with standard stereo; the candidate has no potential tetrahedral center.")
m3aa["stereochemistry"]["double_bonds"] = fact(
    "standard stereo diagnostic", diag_3aa["double_bonds"], e3aa,
    method="structure_derivation", derivation="Calculated with standard stereo; carbonyl/P=O bonds are not E/Z assignments.")
m3aa["stereochemistry"]["axial"] = missing(e3aa, status="not_applicable")
for key, value in (("chirality_status", "achiral"),
                   ("geometry_status", "not_applicable"),
                   ("representation_status", "not_applicable")):
    m3aa["stereochemistry"][key] = fact(
        value, value, e3aa, method="structure_derivation",
        derivation="Status assigned from the corrected candidate graph and standard stereo diagnostics; no sample-level source verification is claimed.")
m3aa["stereochemistry"]["isomeric_smiles"] = fact(
    smiles_3aa, diag_3aa["canonical_isomeric_smiles"], e3aa,
    method="structure_derivation", derivation="Canonical isomeric graph calculated by standard stereo.")
m3aa["stereochemistry"]["identity_key"] = fact(
    diag_3aa["canonical_isomeric_smiles"], diag_3aa["stereo_graph_key"], e3aa,
    method="structure_derivation", derivation="SHA-256 stereo graph key calculated by standard stereo; not a sample identity.")
m3aa["identity_note"] = fact(
    "Page 7 names 3aa as 1-(1-diphenylphosphinylcyclohexyl)propan-2-one.",
    "Corrected candidate has a cyclohexane carbon bearing phosphinyl and acetonyl substituents; no stereogenic center is present in that graph.",
    e3aa, method="structure_derivation",
    derivation="Candidate interpretation of the explicit product name; chemistry correctness remains needs_review.")


# Water/ice is explicitly used in the page-8 workup and therefore gets an
# independent registry entry rather than being hidden in free text.
if mid("water") not in molecules:
    mw = st.new_molecule(mid("water"), DOC, ev(8))
    mw["identity_status"] = "source_identified"
    mw["occurrences"] = [{
        "document_id": DOC,
        "compound_label": fact("ice; water", "water", ev(8)),
        "evidence_ids": ev(8),
    }]
    mw["names"]["original"] = fact("ice; water", "water", ev(8))
    mw["names"]["english"] = fact("water", "water", ev(8))
    mw["names"]["abbreviations"] = fact("H2O", "H2O", ev(8))
    mw["names"]["aliases"] = fact(["ice", "water"], ["ice", "water"], ev(8))
    mw["cas_numbers"] = missing(ev(8))
    mw["representation_type"] = "smiles"
    mw["smiles"] = fact("water", "O", ev(8), method="structure_derivation",
                         derivation="Standard molecular graph for explicitly named water/ice; no external lookup used.")
    mw["generic_template"] = None
    d = st.stereo_inspect("O")
    mw["stereochemistry"]["atom_centers"] = fact("standard stereo diagnostic", d["atom_centers"], ev(8), method="structure_derivation", derivation="Calculated by standard stereo.")
    mw["stereochemistry"]["double_bonds"] = fact("standard stereo diagnostic", d["double_bonds"], ev(8), method="structure_derivation", derivation="Calculated by standard stereo.")
    for key in ("axial", "mixture", "assignment"):
        mw["stereochemistry"][key] = missing(ev(8), status="not_applicable")
    for key, value in (("chirality_status", "achiral"), ("geometry_status", "not_applicable"), ("representation_status", "not_applicable")):
        mw["stereochemistry"][key] = fact(value, value, ev(8), method="structure_derivation", derivation="No molecular stereogenic element in the standard graph.")
    mw["stereochemistry"]["isomeric_smiles"] = fact("O", d["canonical_isomeric_smiles"], ev(8), method="structure_derivation", derivation="Calculated by standard stereo.")
    mw["stereochemistry"]["identity_key"] = fact(d["canonical_isomeric_smiles"], d["stereo_graph_key"], ev(8), method="structure_derivation", derivation="Calculated by standard stereo; not a sample identity.")
    mw["identity_note"] = fact("ice; water", "Same H2O registry molecule used for the ice quench and water wash.", ev(8))
    mw["external_properties"] = []
    molecules[mid("water")] = mw
    bundle["molecules"].append(mw)


# Page 8/9: full bromination procedure and isolated product data.
brom = reactions[rid("scheme2b_bromination_4")]
for page in (8, 9):
    add_reaction_page(brom, page)
add_molecule_occurrence_evidence("4", 8)
add_molecule_occurrence_evidence("4", 9)
participant_by_id(brom, "p_hbr")["charges"]["volume"] = fact(q(0.44, "mL"), q(0.44, "mL"), ev(8))
participant_by_id(brom, "p_hbr")["charges"]["amount"] = fact(q(3, "mmol"), q(3, "mmol"), ev(8))
participant_by_id(brom, "p_dmso")["charges"]["volume"] = fact(q(1.7, "mL"), q(1.7, "mL"), ev(8))
brom["participants"].extend([
    new_participant("p_water", mid("water"), "workup_agent", ev(8), stage="workup"),
    new_participant("p_silica", mid("silica_gel"), "purification_agent", ev(8), stage="purification"),
])
brom["steps"] = [
    new_step("dissolve", "reaction", ["p_3a", "p_dmso"], "stir 3a in DMSO", ev(8)),
    new_step("add_hbr", "reaction", ["p_3a", "p_dmso", "p_hbr"], "add 48% aqueous HBr slowly", ev(8)),
    new_step("hold_55", "reaction", ["p_3a", "p_dmso", "p_hbr"], "stir in an oil bath in an open flask", ev(8), q(55, "°C"), q(4, "h"), "open flask"),
    new_step("ice_quench", "workup", ["p_4", "p_water"], "pour the solution onto ice", ev(8)),
    new_step("isolate", "workup", ["p_4", "p_water"], "filter the solid, wash with water, and dry under vacuum", ev(8)),
    new_step("chromatography", "purification", ["p_4", "p_silica"], "purify the resulting mixture by silica gel chromatography", ev(8)),
]
brom["reaction_scale"] = fact(
    "3a (1 mmol)", {"reference_molecule_id": mid("3a"), "amount": q(1, "mmol")}, ev(8))
brom["solvent_conditions"]["composition"] = fact("DMSO", ["DMSO"], ev(8))
brom["solvent_conditions"]["total_volume"] = fact("1.7 mL", q(1.7, "mL"), ev(8))
brom["procedure"]["original_text"] = fact(
    "3a (1 mmol) in DMSO (1.7 mL); slowly add 0.44 mL (3 mmol) of 48% aqueous HBr; open flask, 55 °C, 4 h; ice, filtration, water wash, vacuum drying, silica chromatography.",
    "Page-8 bromination procedure transcription", ev(8))
brom["procedure"]["addition_order"] = fact(
    "3a in DMSO, then slowly 48% aqueous HBr", ["3a", "DMSO", "48% aqueous HBr"], ev(8))
brom["notes"].append(fact(
    "48% aqueous HBr", {"material": "aqueous HBr", "concentration": q(48, "%")}, ev(8)))
yb = outcome_by_id(brom, "yield")
update_fact_evidence(yb["value"], ev(9))
yb["yield_type"] = fact("56% yield after silica gel chromatography", "isolated", ev(8) + ev(9))
yb["measurement_basis"] = fact("56% yield, 239 mg; product isolated after chromatography", "isolated product", ev(8) + ev(9))
brom["outcomes"].append(product_amount_outcome("product_amount", "4", q(239, "mg"), ev(9), "reported isolated product mass"))
brom["location"] = fact("Scheme 2b and experimental procedure/characterization", "PDF pages 3, 8-9", ev(3) + ev(8) + ev(9))
brom["review"]["unresolved"] = [
    "Absolute/relative configurations of the two potential stereocenters in 4 are not assigned."
]


# Page 9: full allenol-to-enone control procedure and isolated product data.
ctrl = reactions[rid("scheme3a_1a_to_5_with_tfoh")]
add_reaction_page(ctrl, 9)
add_molecule_occurrence_evidence("5", 9)
participant_by_id(ctrl, "p_tfoh")["charges"]["amount"] = fact(q(0.3, "mmol"), q(0.3, "mmol"), ev(9))
ctrl["participants"].append(new_participant("p_silica", mid("silica_gel"), "purification_agent", ev(9), stage="purification"))
ctrl["steps"] = [
    new_step("charge", "reaction", ["p_1a", "p_tfoh", "p_thf"], "add 1a and TfOH in THF to a Schlenk tube equipped with a magnetic stir bar", ev(9), atmosphere="Ar"),
    new_step("hold_100", "reaction", ["p_1a", "p_tfoh", "p_thf"], "stir in an oil bath", ev(9), q(100, "°C"), q(24, "h"), "Ar"),
    new_step("chromatography", "purification", ["p_5", "p_silica"], "directly purify the resulting mixture by silica gel chromatography", ev(9)),
]
ctrl["reaction_scale"] = fact("1a (1 mmol)", {"reference_molecule_id": mid("1a"), "amount": q(1, "mmol")}, ev(9))
ctrl["solvent_conditions"]["composition"] = fact("THF", ["THF"], ev(9))
ctrl["solvent_conditions"]["total_volume"] = fact("5 mL", q(5, "mL"), ev(9))
ctrl["procedure"]["original_text"] = fact(
    "1-phenyl-2,3-butadien-1-ol (1 mmol), TfOH (0.3 mmol), THF (5 mL), Ar, 100 °C oil bath, 24 h; direct silica chromatography.",
    "Page-9 procedure transcription", ev(9))
ctrl["procedure"]["addition_order"] = fact("1a + TfOH in THF, then heat", ["1a", "TfOH", "THF", "heat"], ev(9))
yc = outcome_by_id(ctrl, "yield_5")
update_fact_evidence(yc["value"], ev(9))
yc["yield_type"] = fact("51% yield after direct silica chromatography", "isolated", ev(9))
yc["measurement_basis"] = fact("51% yield, 75 mg; product isolated after chromatography", "isolated product", ev(9))
ctrl["outcomes"].append(product_amount_outcome("product_amount", "5", q(75, "mg"), ev(9), "reported isolated product mass"))
ctrl["location"] = fact("Scheme 3a and Procedure for the Synthesis of 5", "PDF pages 3 and 9", ev(3) + ev(9))
ctrl["review"]["unresolved"] = [
    "Water/oxygen provenance in the net allenol-to-enone rearrangement is not resolved by the experiment."
]


# Page 9: phospha-Michael procedure.  The heading and Scheme 3b identify
# enone 5, but the procedure sentence literally says compound 4.  Preserve the
# conflict and do not silently replace one source reading with the other.
addn = reactions[rid("scheme3b_enone_5_to_3a")]
add_reaction_page(addn, 9)
participant_by_id(addn, "p_tfoh")["charges"]["amount"] = fact(q(0.12, "mmol"), q(0.12, "mmol"), ev(9))
addn["participants"].extend([
    new_participant("p_silica", mid("silica_gel"), "purification_agent", ev(9), stage="purification"),
    new_participant("p_pe", mid("petroleum_ether"), "purification_agent", ev(9), stage="purification"),
    new_participant("p_ea", mid("EtOAc"), "purification_agent", ev(9), stage="purification"),
])
addn["steps"] = [
    new_step("charge", "reaction", ["p_5", "p_2a", "p_tfoh", "p_thf"], "charge the substrate named by the heading/Scheme 3b as enone 5 with 2a and TfOH in THF", ev(9), atmosphere="Ar"),
    new_step("hold_100", "reaction", ["p_5", "p_2a", "p_tfoh", "p_thf"], "stir in an oil bath", ev(9), q(100, "°C"), q(24, "h"), "Ar"),
    new_step("chromatography", "purification", ["p_3a", "p_silica", "p_pe", "p_ea"], "directly purify by silica gel chromatography using petroleum ether/ethyl acetate", ev(9)),
]
heading = fact("Procedure for the Synthesis of 3a via the Addition of 2a to 5", "enone 5", ev(9))
sentence = fact("was added 4 (0.4 mmol)", "compound 4", ev(9))
addn["procedure"]["original_text"] = fact(
    evidence=ev(9), status="conflict", alternatives=[heading, sentence])
addn["procedure"]["addition_order"] = fact(
    evidence=ev(9), status="conflict", alternatives=[
        fact("heading/Scheme 3b: 5 (0.4 mmol), 2a, TfOH, THF", ["5", "2a", "TfOH", "THF"], ev(3) + ev(9)),
        fact("procedure sentence: 4 (0.4 mmol), 2a, TfOH, THF", ["4", "2a", "TfOH", "THF"], ev(9)),
    ])
addn["reaction_scale"] = fact(
    evidence=ev(3) + ev(9), status="conflict", alternatives=[
        fact("Scheme 3b/heading: 5 (0.4 mmol)", {"reference_molecule_id": mid("5"), "amount": q(0.4, "mmol")}, ev(3) + ev(9)),
        fact("procedure sentence: 4 (0.4 mmol)", {"reference_molecule_id": mid("4"), "amount": q(0.4, "mmol")}, ev(9)),
    ])
addn["solvent_conditions"]["composition"] = fact("THF", ["THF"], ev(9))
addn["solvent_conditions"]["total_volume"] = fact("2 mL", q(2, "mL"), ev(9))
ya = outcome_by_id(addn, "yield")
ya["yield_type"] = fact("96% after the page-9 chromatographic procedure", "isolated", ev(3) + ev(9), method="inherited_general_condition")
ya["measurement_basis"] = fact("page-3 96% yield paired with page-9 direct chromatography", "isolated product after chromatography", ev(3) + ev(9), method="inherited_general_condition")
addn["notes"].append(fact(
    "Heading and Scheme 3b identify enone 5, but the experimental sentence literally reads compound 4 (0.4 mmol).",
    "Reactant identity conflict retained; the reaction participant remains 5 because the heading, Scheme 3b, product, and mechanistic discussion all identify 5, but this needs expert review.",
    ev(3) + ev(9)))
addn["location"] = fact("Scheme 3b and conflicting experimental paragraph", "PDF pages 3 and 9", ev(3) + ev(9))
addn["review"]["unresolved"] = [
    "Reaction identity conflict: Scheme 3b and the page-9 heading name enone 5, while the procedure sentence says compound 4 (0.4 mmol).",
    "Complete proton-transfer/coproduct balance remains unresolved."
]


# Source context: spectral blocks and reference lists are not reactions.  The
# page-9 SI notice is recorded without associating an SI file.
new_context = [
    {
        "context_id": f"{DOC}:product_characterization_pages_7_9",
        "document_id": DOC,
        "kind": "characterization_only",
        "description": "NMR, IR, HRMS, melting-point and physical-form data on pages 7-9 were associated with existing products where relevant; individual peaks, exact masses and fragment-like spectral entries were not promoted to synthesis experiments.",
        "evidence_ids": ev(7) + ev(8) + ev(9),
    },
    {
        "context_id": f"{DOC}:si_notice_page9_unassociated",
        "document_id": DOC,
        "kind": "source_limitation",
        "description": "Page 9 advertises crystallographic/NMR Supporting Information, but source_mapping.json does not explicitly associate an SI file; no SI content was inferred or linked.",
        "evidence_ids": ev(9),
    },
    {
        "context_id": f"{DOC}:reference_list_pages_10_11",
        "document_id": DOC,
        "kind": "bibliography_only",
        "description": "Pages 10-11 contain cited-reference text only. They do not report additional experiments of the current paper and were not converted into current-paper reaction records.",
        "evidence_ids": ev(10) + ev(11),
    },
]
for context in new_context:
    if not any(x.get("context_id") == context["context_id"] for x in bundle["source_context"]):
        bundle["source_context"].append(context)


# Save explicit standard-stereo diagnostics for every molecule product/reagent
# whose identity is materially used on pages 7-9.  Unencoded mixtures retain
# their source assignments instead of receiving invented single-isomer SMILES.
review_labels = [
    "3w", "3x", "3y", "3z", "3aa", "3ab", "3ac", "3ad", "3ae",
    "3af", "3ag", "3ah", "3ai", "3aj", "3ak", "3al", "3am", "3an",
    "3ap", "4", "5", "3a", "1a", "2a", "TfOH", "THF", "DMSO",
    "HBr", "water", "silica_gel", "petroleum_ether", "EtOAc",
]
stereo_review = {}
for label in review_labels:
    m = molecules[mid(label)]
    smi = m["smiles"]["normalized_value"]
    stereo_review[mid(label)] = {
        "source_pages": [7, 8, 9],
        "chirality_status": m["stereochemistry"]["chirality_status"]["normalized_value"],
        "geometry_status": m["stereochemistry"]["geometry_status"]["normalized_value"],
        "representation_status": m["stereochemistry"]["representation_status"]["normalized_value"],
        "diagnostic": st.stereo_inspect(smi) if smi else None,
        "note": "No single-isomer graph asserted" if smi is None else "Graph diagnostic only; source_verified remains false",
    }


coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 7,
            "state": "attempted",
            "regions": [
                {
                    "region": "two-column experimental text, top-left continuation through bottom-right",
                    "finding": "Continuation of 3w characterization followed by product/characterization paragraphs for 3x-3ag; reported masses and yields were merged into the corresponding existing Table 2 experiments, while spectral peaks were retained as characterization context.",
                    "evidence_ids": ev(7),
                }
            ],
            "record_ids": [rid(f"table2_scope_3{x}") for x in ["w", "x", "y", "z", "aa", "ab", "ac", "ad", "ae", "af", "ag"]],
            "unresolved": [
                "The 3ag characterization paragraph continues onto page 8.",
                "Potential stereocenters remain configuration-unassigned unless explicitly stated; successful SMILES parsing is not chemistry verification."
            ],
        },
        {
            "page": 8,
            "state": "attempted",
            "regions": [
                {
                    "region": "left and upper-right product paragraphs",
                    "finding": "Continuation of 3ag and product/characterization paragraphs for 3ah-3ai and 3aj-3ap; masses and isolated-yield basis were merged into existing Table 2/Table 3 records, not made into duplicate experiments.",
                    "evidence_ids": ev(8),
                },
                {
                    "region": "bottom-right Procedure for the Synthesis of 4",
                    "finding": "Full bromination procedure: 3a/DMSO, slow addition of 48% aqueous HBr, open-flask 55 °C/4 h hold, ice quench, filtration/wash/vacuum dry, and silica chromatography; merged with Scheme 2b record.",
                    "evidence_ids": ev(8),
                },
            ],
            "record_ids": [
                rid("table2_scope_3ag"), rid("table2_scope_3ah"), rid("table2_scope_3ai"),
                *[rid(f"table3_scope_3{x}") for x in ["aj", "ak", "al", "am", "an", "ap"]],
                rid("scheme2b_bromination_4"),
            ],
            "unresolved": [
                "3an is reported as two isomers but their configurations/order and the dr measurement method remain unassigned/unreported.",
                "3ap is reported as a racemic form with dr 3.4:1; individual stereoisomer configurations remain unassigned."
            ],
        },
        {
            "page": 9,
            "state": "attempted",
            "regions": [
                {
                    "region": "top-left product 4 characterization",
                    "finding": "Isolated product 4 amount/yield and characterization; amount/yield were linked to the existing bromination experiment, while spectra were kept as context.",
                    "evidence_ids": ev(9),
                },
                {
                    "region": "left-column Procedure for the Synthesis of 5 and product 5 paragraph",
                    "finding": "Full mechanistic control procedure and isolated 5 amount/yield; merged with the existing Scheme 3a experiment.",
                    "evidence_ids": ev(9),
                },
                {
                    "region": "left-column Procedure for the Synthesis of 3a via Addition of 2a to 5",
                    "finding": "The heading and Scheme 3b name enone 5, but the procedure sentence literally says compound 4 (0.4 mmol); both readings are preserved as a source conflict in the existing Scheme 3b record.",
                    "evidence_ids": ev(9),
                },
                {
                    "region": "Associated Content and author/reference text",
                    "finding": "Supporting Information is advertised but no SI file is explicitly mapped; author metadata/references are source context, not experiments.",
                    "evidence_ids": ev(9),
                },
            ],
            "record_ids": [rid("scheme2b_bromination_4"), rid("scheme3a_1a_to_5_with_tfoh"), rid("scheme3b_enone_5_to_3a")],
            "unresolved": [
                "Reactant identity conflict in the 3a-from-enone procedure requires expert resolution (5 in heading/Scheme 3b versus literal 4 in procedure sentence).",
                "Advertised Supporting Information was not supplied/mapped and remains unassociated."
            ],
        },
        {
            "page": 10,
            "state": "attempted",
            "regions": [
                {
                    "region": "entire two-column page",
                    "finding": "Reference list entries only (continuation of citations 3-31); no current-paper reaction arrow, table row, experimental procedure, or new synthesis record.",
                    "evidence_ids": ev(10),
                }
            ],
            "record_ids": [],
            "unresolved": ["Cited literature was not supplied and was not re-extracted as experiments of the current paper."],
        },
        {
            "page": 11,
            "state": "attempted",
            "regions": [
                {
                    "region": "entire page",
                    "finding": "Final reference list entries 31-33 only; no current-paper reaction arrow, table row, experimental procedure, or new synthesis record.",
                    "evidence_ids": ev(11),
                }
            ],
            "record_ids": [],
            "unresolved": ["Cited literature was not supplied and was not re-extracted as experiments of the current paper."],
        },
    ],
}


# The standard merge accepts this complete, internally consistent bundle and
# automatically regenerates legacy compatibility views.  No second hand-built
# legacy format is written here.
for reaction in bundle["reactions"]:
    reaction["review"]["status"] = "needs_review"
    reaction["review"]["method"] = "model"
    reaction["review"]["reviewer"] = None
    reaction["review"]["checked_fields"] = []

st.write(RAW / "extracted_bundle_pages_0007_0011.json", bundle)
st.write(RAW / "stereo_checks_pages_0007_0011.json", stereo_review)
st.write(CURRENT / "coverage" / "pages_0007_0011.json", coverage)
st.write(SUBMITTED, bundle)

print(json.dumps({
    "submitted_reactions": len(bundle["reactions"]),
    "updated_reactions_in_page_range": 22,
    "new_reaction_ids": 0,
    "molecules": len(bundle["molecules"]),
    "new_molecule_ids": [mid("water")],
    "source_context": len(bundle["source_context"]),
}, ensure_ascii=False, indent=2))
