#!/usr/bin/env python3
"""Rebuild the current document and add source-derived records for PDF pages 13-18.

Pages 13-15 contain the experimental section.  Existing same-document records
are enriched with image evidence instead of duplicated.  Three preparations
first fully described on page 15 are added.  Pages 16-18 are analytical-method
and reference matter and therefore contribute source context/coverage only.
"""

from __future__ import annotations

import json
import runpy
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = ROOT / "runs" / "current"
DOC = "doc_837adcd3d5d8cb7b"

# Reconstruct all source-derived pages 1-12 records.  This chain starts from a
# standard-init template and does not import legacy compatibility experiments.
ns = runpy.run_path(str(HERE / "extend_pages_7_12.py"))
M = ns["M"]
R = ns["R"]
add_mol = ns["add_mol"]
participant = ns["participant"]
step = ns["step"]
outcome = ns["outcome"]
reaction = ns["reaction"]
f = ns["f"]
q = ns["q"]
write = ns["write"]
stereo_inspect = ns["stereo_inspect"]
source_mentions = ns["source_mentions"]


def et(page: int) -> str:
    return f"{DOC}:page{page}:text"


def ei(page: int) -> str:
    return f"{DOC}:page{page}:image"


def ev(*pages: int, text: bool = True, image: bool = True) -> list[str]:
    result: list[str] = []
    for page in pages:
        if text:
            result.append(et(page))
        if image:
            result.append(ei(page))
    return result


def pids(parts):
    return [part["participant_id"] for part in parts]


def simple_parts(rid, evidence, specs):
    return [
        participant(
            rid,
            index,
            slug,
            role,
            evidence,
            amounts=amounts,
            reference=reference,
            product_role=product_role,
            stage=stage_name,
        )
        for index, (slug, role, amounts, reference, product_role, stage_name) in enumerate(specs, 1)
    ]


def inherited_argon(evidence):
    return f(
        "All manipulations were carried out in an oxygen- and moisture-free argon atmosphere",
        "argon",
        evidence,
        method="inherited_general_condition",
    )


# The prior extraction used pages 13-15 text as matching procedure context but
# could not yet claim those images as inspected.  Add the now-inspected image
# evidence to every existing reaction/molecule occurrence already supported by
# the corresponding page text.  This merges repeated descriptions rather than
# creating duplicate experiments.
for page in (13, 14, 15):
    text_id = et(page)
    image_id = ei(page)
    for record in R:
        if text_id in record["evidence_ids"] and image_id not in record["evidence_ids"]:
            record["evidence_ids"].append(image_id)
    for molecule in M.values():
        occurrence = molecule["occurrences"][0]
        if text_id in occurrence["evidence_ids"] and image_id not in occurrence["evidence_ids"]:
            occurrence["evidence_ids"].append(image_id)

# Correct the page lineage of the two printed charge conflicts now that the
# experimental pages themselves are in scope.  The values are preserved; only
# their exact page/image evidence is repaired.
complex6_record = next(item for item in R if item["reaction_id"].endswith(":isolation_complex6"))
complex6_volume = complex6_record["participants"][1]["charges"]["volume"]
complex6_volume["evidence_ids"] = ev(14)
for alternative in complex6_volume["alternatives"]:
    alternative["evidence_ids"] = ev(14)

complex9a_record = next(item for item in R if item["reaction_id"].endswith(":synthesis_complex9_from7"))
for evidence_id in ev(14, 15):
    if evidence_id not in complex9a_record["evidence_ids"]:
        complex9a_record["evidence_ids"].append(evidence_id)
complex9a_record["participants"][0]["charges"]["mass"]["evidence_ids"] = ev(14)
complex9a_amount = complex9a_record["participants"][0]["charges"]["amount"]
complex9a_amount["evidence_ids"] = ev(14)
for alternative in complex9a_amount["alternatives"]:
    alternative["evidence_ids"] = ev(14)
for key in ("mass", "amount"):
    complex9a_record["participants"][1]["charges"][key]["evidence_ids"] = ev(14)
complex9a_record["steps"][0]["operation"]["evidence_ids"] = ev(14, 15)
complex9a_record["steps"][0]["time"]["evidence_ids"] = ev(14, 15)


e15 = ev(15)
e13_15 = ev(13, 15)

# Molecules first required by the page-15 preparations.  Simple ionic/organic
# candidates are graph-checked; the organozirconium precursor remains drawing-
# based because forcing a Cp/Zr coordination graph would be unsafe.
add_mol(
    "cp2zr_me_cl",
    "Cp2Zr(Me)Cl",
    "Cp2Zr(Me)Cl",
    e15,
    chirality="unspecified",
    geometry="unresolved",
    rep="not_encoded",
    identity_note="Source-named organozirconium precursor; complete Cp/Zr coordination is retained from the name rather than forced into SMILES.",
)
add_mol(
    "li_c2sime3",
    "lithium trimethylsilylacetylide",
    "LiC≡CSiMe3",
    e15,
    smiles="[Li+].C[Si](C)(C)C#[C-]",
    chirality="achiral",
    geometry="not_applicable",
    rep="encoded",
)
add_mol(
    "n_hexane",
    "n-hexane",
    "n-hexane",
    e15,
    smiles="CCCCCC",
    chirality="achiral",
    geometry="not_applicable",
    rep="encoded",
)
add_mol(
    "n_buli",
    "n-butyllithium",
    "n-BuLi",
    e15,
    smiles="[Li+].[CH2-]CCC",
    abbrev="n-BuLi",
    chirality="achiral",
    geometry="not_applicable",
    rep="encoded",
)
add_mol(
    "bh3_sme2",
    "borane-dimethyl sulfide complex",
    "H3B·SMe2",
    e15,
    smiles="[BH3-][S+](C)C",
    abbrev="BH3·SMe2",
    chirality="achiral",
    geometry="not_applicable",
    rep="encoded",
)

# Existing products/reagents gain the inspected page-15 image occurrence.
for slug in (
    "complex11",
    "licl",
    "toluene",
    "dimethylamine_borane",
    "li_amidoborane",
    "thf",
    "pentane",
    "mu_aminodiborane",
    "benzene",
    "toluene_d8",
):
    occurrence = M[slug]["occurrences"][0]
    for evidence_id in e15:
        if evidence_id not in occurrence["evidence_ids"]:
            occurrence["evidence_ids"].append(evidence_id)


# Page 15: preparation of complex 11.
rid = f"{DOC}:synthesis_complex11"
parts = simple_parts(
    rid,
    e15,
    [
        ("cp2zr_me_cl", "reactant", {"mass": q(2.88, "g"), "amount": q(10.59, "mmol"), "equivalents": q(1, "equiv")}, "cp2zr_me_cl", None, "reaction"),
        ("li_c2sime3", "reagent", {"mass": q(1.12, "g"), "amount": q(12.71, "mmol"), "equivalents": q(1.2, "equiv")}, "cp2zr_me_cl", None, "reaction"),
        ("toluene", "solvent", {"volume": q(100, "mL")}, None, None, "reaction"),
        ("n_hexane", "purification_agent", {}, None, None, "purification"),
        ("complex11", "product", {}, None, "isolated product", "reaction"),
        ("licl", "product", {}, None, "salt coproduct", "reaction"),
    ],
)
steps = [
    step(rid, 1, e15, pids(parts[:2]), "weigh the two solids into a 250 mL Schlenk flask in the glove box", None, None, "glovebox"),
    step(rid, 2, e15, pids(parts[:3]), "add 100 mL of pre-cooled toluene to the white mixture", q(-78, "°C")),
    step(rid, 3, e15, pids(parts[:3]), "allow the white suspension to slowly warm from -78 °C to room temperature", "warming to room temperature", "overnight"),
    step(rid, 4, e15, pids(parts), "filter off the supernatant and remove solvent from the filtrate under vacuum", stage="workup"),
    step(rid, 5, e15, pids(parts), "repeatedly recrystallize from saturated n-hexane solutions until the filter cake is grey", stage="purification"),
    step(rid, 6, e15, pids(parts), "dry the grey needles under vacuum", stage="purification"),
]
for item in steps[1:]:
    item["atmosphere"] = inherited_argon(ev(13))
outs = [
    outcome(rid, 1, "product_amount", ["complex11"], e15, q(1.08, "g"), basis="isolated grey needles", role="isolated product"),
    outcome(rid, 2, "yield", ["complex11"], e15, q(32, "%"), basis="isolated grey needles", role="isolated product", yield_type="isolated"),
]
r11 = reaction(
    "synthesis_complex11",
    "specific_experiment",
    e15,
    ["other"],
    "synthesis of complex 11",
    "alkynyl substitution at zirconium",
    parts,
    steps,
    outs,
    class_changes=[
        ("C", "Zr", 0, 1, "formed", "A Zr-C(alkynyl) bond is formed in complex 11."),
        ("Cl", "Zr", 1, 0, "broken", "The Zr-Cl bond is displaced during salt metathesis."),
    ],
    notes=["All synthesis steps were explicitly performed under exclusion of light."],
    unresolved=["The organozirconium structures are not safely encoded as complete SMILES; exact atom mapping remains unresolved."],
)
r11["reaction_scale"] = f(
    "Cp2Zr(Me)Cl 10.59 mmol",
    {"reference_molecule_id": M["cp2zr_me_cl"]["molecule_id"], "amount": q(10.59, "mmol")},
    e15,
)
r11["solvent_conditions"]["composition"] = f("toluene", [M["toluene"]["molecule_id"]], e15)
r11["solvent_conditions"]["total_volume"] = f("100 mL", q(100, "mL"), e15)
r11["special_conditions"]["light_source"] = f("under exclusion of light", "excluded", e15)
r11["special_conditions"]["reactor"] = f("250 mL Schlenk flask", "250 mL Schlenk flask", e15)
r11["procedure"]["original_text"] = f(
    "Cp2Zr(Me)Cl (2.88 g, 10.59 mmol, 1 eq) and LiC≡CSiMe3 (1.12 g, 12.71 mmol, 1.2 eq) were weighed into a 250 mL Schlenk flask; 100 mL pre-cooled (-78 °C) toluene was added and the suspension was allowed to slowly warm to room temperature overnight.",
    "Cp2Zr(Me)Cl and LiC≡CSiMe3 were combined under exclusion of light, treated with pre-cooled toluene, and warmed from -78 °C to room temperature overnight.",
    e15,
)
r11["procedure"]["addition_order"] = f(
    "solids first; pre-cooled toluene added at -78 °C; then slow warming",
    ["Cp2Zr(Me)Cl", "LiC≡CSiMe3", "pre-cooled toluene", "slow warming"],
    e15,
)
r11["procedure"]["pretreatment"] = f("all synthesis steps under exclusion of light; toluene pre-cooled to -78 °C", "exclude light and pre-cool toluene to -78 °C", e15)


# Page 15: synthesis of Li[NMe2BH3].  The evolved gas is not identified by the
# source, so no gas product is invented.
rid = f"{DOC}:synthesis_li_amidoborane"
parts = simple_parts(
    rid,
    e15,
    [
        ("dimethylamine_borane", "reactant", {"mass": q(400, "mg"), "amount": q(6.80, "mmol")}, None, None, "reaction"),
        ("n_buli", "reagent", {"volume": q(2.90, "mL"), "amount": q(7.50, "mmol")}, "dimethylamine_borane", None, "reaction"),
        ("thf", "solvent", {"volume": q(15, "mL")}, None, None, "reaction"),
        ("n_hexane", "solvent", {}, None, "carrier solvent for n-BuLi", "reaction"),
        ("pentane", "purification_agent", {"volume": q(10, "mL", raw_text="2 × 5 mL")}, None, None, "purification"),
        ("li_amidoborane", "product", {}, None, "isolated product", "reaction"),
    ],
)
steps = [
    step(rid, 1, e15, pids(parts[:3]), "dissolve H3B·NMe2H in THF", q(-78, "°C")),
    step(rid, 2, e15, pids(parts[:4]), "add n-BuLi (2.50 M in n-hexane) dropwise; gas evolution was noticed", q(-78, "°C")),
    step(rid, 3, e15, pids(parts[:4]), "warm the reaction from -78 °C to room temperature", "warming to room temperature"),
    step(rid, 4, e15, pids(parts[:4]), "stir at room temperature", "room temperature", "overnight"),
    step(rid, 5, e15, pids(parts), "remove all volatiles in vacuo", stage="workup"),
    step(rid, 6, e15, pids(parts), "wash the precipitate with pentane 2 × 5 mL", stage="purification"),
]
for item in steps:
    item["atmosphere"] = inherited_argon(ev(13))
outs = [
    outcome(rid, 1, "product_amount", ["li_amidoborane"], e15, q(350, "mg"), basis="isolated white solid", role="isolated product"),
    outcome(rid, 2, "yield", ["li_amidoborane"], e15, q(87, "%"), basis="isolated white solid", role="isolated product", yield_type="isolated"),
]
rli = reaction(
    "synthesis_li_amidoborane",
    "specific_experiment",
    e13_15,
    ["other"],
    "synthesis of Li[NMe2BH3]",
    "deprotonation of dimethylamine borane",
    parts,
    steps,
    outs,
    class_changes=None,
    notes=["The n-BuLi solution concentration is 2.50 M in n-hexane; the identity of the evolved gas is not stated."],
    unresolved=["The evolved gas is unidentified, and Li/N ion pairing is not treated as a confidently mapped covalent bond change."],
)
rli["reaction_scale"] = f(
    "H3B·NMe2H 6.80 mmol",
    {"reference_molecule_id": M["dimethylamine_borane"]["molecule_id"], "amount": q(6.80, "mmol")},
    e15,
)
rli["solvent_conditions"]["composition"] = f(
    "15 mL THF plus n-hexane carrier from the n-BuLi solution",
    [M["thf"]["molecule_id"], M["n_hexane"]["molecule_id"]],
    e15,
)
rli["solvent_conditions"]["concentration"] = f("n-BuLi 2.50 M in n-hexane", q(2.50, "M"), e15)
rli["solvent_conditions"]["concentration_basis"] = f("n-BuLi solution", M["n_buli"]["molecule_id"], e15)
rli["special_conditions"]["reactor"] = f("Schlenk flask", "Schlenk flask", e15)
rli["procedure"]["original_text"] = f(
    "H3B·NMe2H (400 mg, 6.80 mmol) in 15 mL THF at -78 °C was treated dropwise with n-BuLi (2.90 mL, 7.50 mmol, 2.50 M in n-hexane), warmed to room temperature, stirred overnight, evaporated, and washed with pentane (2 × 5 mL).",
    "deprotonation in THF followed by evaporation and pentane washing",
    e15,
)
rli["procedure"]["addition_order"] = f(
    "H3B·NMe2H/THF first; n-BuLi added dropwise at -78 °C; warm to room temperature",
    ["H3B·NMe2H in THF", "dropwise n-BuLi", "warm to room temperature"],
    e15,
)
rli["procedure"]["pretreatment"] = f("THF solution held at -78 °C before n-BuLi addition", "cool to -78 °C", e15)


# Page 15: NMR-scale preparation/observation of mu-aminodiborane.
rid = f"{DOC}:synthesis_mu_aminodiborane"
parts = simple_parts(
    rid,
    e15,
    [
        ("dimethylamine_borane", "reactant", {"mass": q(30, "mg"), "amount": q(0.51, "mmol")}, None, None, "reaction"),
        ("bh3_sme2", "reagent", {"volume": q(48, "uL"), "amount": q(0.51, "mmol")}, "dimethylamine_borane", None, "reaction"),
        ("benzene", "solvent", {"volume": q(0.6, "mL")}, None, None, "reaction"),
        ("mu_aminodiborane", "product", {}, None, "NMR-observed product", "reaction"),
        ("toluene_d8", "solvent", {}, None, "NMR solvent", "characterization"),
    ],
)
steps = [
    step(rid, 1, e15, pids(parts[:3]), "dissolve H3B·NMe2H in benzene in a Young NMR tube and add H3B·SMe2", None),
    step(rid, 2, e15, pids(parts[:4]), "heat the colourless solution", q(80, "°C"), "overnight"),
    step(rid, 3, e15, pids(parts), "observe the product by 11B NMR in toluene-d8", q(25, "°C"), None, stage="characterization"),
]
for item in steps:
    item["atmosphere"] = inherited_argon(ev(13))
outs = [
    outcome(
        rid,
        1,
        "other",
        ["mu_aminodiborane"],
        e15,
        "product observed; 11B NMR delta -19.6 ppm (m)",
        basis="11B NMR, 25 °C, toluene-d8, 96.32 MHz",
        role="NMR-observed product",
    )
]
rmu = reaction(
    "synthesis_mu_aminodiborane",
    "specific_experiment",
    e13_15,
    ["mechanistic_study", "other"],
    "synthesis of BH2(mu-Me2N)(mu-H)BH2",
    "aminodiborane formation",
    parts,
    steps,
    outs,
    class_changes=[
        ("B", "N", 0, 1, "formed", "A new B-N connection is inferred on forming the dimethylamino-bridged diborane."),
    ],
    notes=["The source reports NMR observation only; no isolated mass or yield is given."],
    unresolved=["Product isolation, yield, coproduct identity, and complete atom/hydrogen mapping are unreported."],
)
rmu["reaction_scale"] = f(
    "H3B·NMe2H 0.51 mmol",
    {"reference_molecule_id": M["dimethylamine_borane"]["molecule_id"], "amount": q(0.51, "mmol")},
    e15,
)
rmu["solvent_conditions"]["composition"] = f("benzene", [M["benzene"]["molecule_id"]], e15)
rmu["solvent_conditions"]["total_volume"] = f("0.6 mL benzene", q(0.6, "mL"), e15)
rmu["special_conditions"]["reactor"] = f("Young NMR tube", "Young NMR tube", e15)
rmu["procedure"]["original_text"] = f(
    "H3B·NMe2H (30 mg, 0.51 mmol) was dissolved in 0.6 mL benzene in a Young NMR tube; H3B·SMe2 (48 μL, 0.51 mmol) was added and the mixture was heated to 80 °C overnight.",
    "equimolar H3B·NMe2H/H3B·SMe2 heated in benzene",
    e15,
)
rmu["procedure"]["addition_order"] = f(
    "dissolve H3B·NMe2H in benzene, then add H3B·SMe2",
    ["H3B·NMe2H in benzene", "H3B·SMe2"],
    e15,
)


# Add page-bounded non-reaction context.  Characterization peak lists and
# analytical/computational methods are explicitly kept out of the reaction set.
base = json.loads((OUT / "submitted.json").read_text())
contexts = list(base.get("source_context", []))
contexts.extend(
    [
        {
            "context_id": f"{DOC}:context_pages13_15_general_and_characterization",
            "document_id": DOC,
            "evidence_ids": ev(13, 14, 15),
            "kind": "experimental_methods_and_characterization",
            "note": "General inert-atmosphere/solvent/analytical methods and the NMR, IR, Raman, MS, melting-point, CHN, and XRD characterization blocks support the associated preparations; standalone peaks and fragment abundances are not synthetic experiments.",
        },
        {
            "context_id": f"{DOC}:context_page13_dft_methods",
            "document_id": DOC,
            "evidence_ids": ev(13),
            "kind": "computational_methods",
            "note": "BP86/LANL2DZ/TZVP preliminary calculations and B3LYP/def2-TZVP/GD3BJ production calculations, frequency analyses, and IRC confirmation are computational methods, not executed chemical reactions.",
        },
        {
            "context_id": f"{DOC}:context_page15_general_dehydrocoupling_procedure",
            "document_id": DOC,
            "evidence_ids": ev(15),
            "kind": "general_procedure",
            "note": "The general dehydrocoupling procedure supplies 78.4 mg/1.33 mmol H3B·NMe2H, 4 mL toluene, glovebox transfer, an Ar-connected gas buret, MeLi solution (1.6 M in Et2O), immediate acquisition, GC-TCD gas analysis, and 11B NMR aliquot analysis to the existing catalytic records; it is not a separate reaction.",
        },
        {
            "context_id": f"{DOC}:context_page16_uvvis_svd_method",
            "document_id": DOC,
            "evidence_ids": ev(16),
            "kind": "spectroscopic_data_analysis",
            "note": "The extinction-diagram and singular-value-decomposition discussion explains UV-vis data analysis and refers back to the already recorded activation experiment; it is not an additional experiment. The printed phrase 'reaction of 5 to 2' is retained as source wording and not used to create a new chemical arrow.",
        },
        {
            "context_id": f"{DOC}:context_pages16_18_references",
            "document_id": DOC,
            "evidence_ids": ev(16, 17, 18),
            "kind": "references_and_notes",
            "note": "Conflicts statement, acknowledgements, notes, and references only. Cited literature and ESI pointers are not imported as current-paper experiments, and no SI is associated.",
        },
    ]
)

# Build from the fresh template generated specifically for this segment.
bundle = json.loads((HERE / "standard_init_pages13_18.json").read_text())
bundle["molecules"] = list(M.values())
bundle["reactions"] = R
bundle["source_context"] = contexts
write(OUT / "submitted.json", bundle)


# Required exact page coverage for pages 13-18.
page_records = {page: [] for page in range(13, 19)}
for record in R:
    for page in page_records:
        if et(page) in record["evidence_ids"] or ei(page) in record["evidence_ids"]:
            page_records[page].append(record["reaction_id"])

regions = {
    13: [
        {"region": "left column: Experimental details—General", "finding": "General inert-atmosphere, solvent-purification, NMR/Raman/IR/MS/CHN/melting-point methods; characterization context, not separate syntheses.", "evidence_ids": ev(13)},
        {"region": "right column top: X-ray and DFT analysis", "finding": "XRD/CCDC and computational methods only; no executed chemical experiment added.", "evidence_ids": ev(13)},
        {"region": "right column lower half: Synthesis of complex 4", "finding": "Concrete preparation of complex 4; merged with existing synthesis_complex4 and image evidence added.", "evidence_ids": ev(13)},
    ],
    14: [
        {"region": "left column top: Synthesis of complex 5", "finding": "Concrete preparation and characterization of complex 5; merged with existing synthesis_complex5.", "evidence_ids": ev(14)},
        {"region": "left column middle: Synthesis of complex 6", "finding": "Concrete preparation of complex 6 with the printed MeLi volume/amount conflict; merged with existing isolation_complex6 without correcting the source.", "evidence_ids": ev(14)},
        {"region": "left-bottom through right-middle: Synthesis of complex 7, Methods A/B", "finding": "Two distinct preparation methods for the same product 7; retained as two existing experiments, not merged by product identity.", "evidence_ids": ev(14)},
        {"region": "right column lower-middle: Synthesis of complex 8", "finding": "Concrete light-excluded low-temperature preparation of 8; merged with existing synthesis_complex8.", "evidence_ids": ev(14)},
        {"region": "right column bottom: Synthesis of complex 9, Method A begins", "finding": "Method A begins here and continues on page 15; merged with existing synthesis_complex9_from7. Printed 25 mg/0.83 mmol conflict remains unresolved.", "evidence_ids": ev(14)},
        {"region": "NMR/IR/Raman/MS/CHN blocks", "finding": "Characterization and fragment-abundance data retained as source context, not counted as additional syntheses.", "evidence_ids": ev(14)},
    ],
    15: [
        {"region": "left column top: complex 9 Methods A/B", "finding": "Completes Method A and gives isolated Method B (97 mg, 72%); merged into the two existing complex-9 preparation records.", "evidence_ids": ev(15)},
        {"region": "left column middle: Synthesis of complexes 10 and 11", "finding": "Complex 10 preparation is merged with its existing record; complex 11 is added as a new concrete experiment with separate low-temperature addition, overnight warming, workup, and recrystallization steps.", "evidence_ids": ev(15)},
        {"region": "right column middle: Synthesis of Li[NMe2BH3]", "finding": "New concrete deprotonation/preparation record with all printed charges, solvent carrier concentration, staged temperatures, pentane washing, 350 mg product, and 87% isolated yield.", "evidence_ids": ev(15)},
        {"region": "right column lower-middle: Synthesis of BH2(mu-Me2N)(mu-H)BH2", "finding": "New NMR-scale preparation/observation record; no isolated yield is invented.", "evidence_ids": ev(15)},
        {"region": "right column bottom: General dehydrocoupling procedure", "finding": "General apparatus/scale/analysis context associated with existing catalytic records; not a separate experiment.", "evidence_ids": ev(15)},
        {"region": "individual NMR/IR/Raman/MS/CHN blocks", "finding": "Characterization-only observations retained in source_context; no peak table or fragment abundance is converted to a synthesis.", "evidence_ids": ev(15)},
    ],
    16: [
        {"region": "left column top: numerical-rank/SVD continuation", "finding": "UV-vis data-analysis method referring to an already recorded activation experiment; no new reaction record.", "evidence_ids": ev(16)},
        {"region": "left column lower half and right column: conflicts, acknowledgements, references 1-17", "finding": "Administrative and cited-literature matter only; cited work is not imported as current-paper experiments.", "evidence_ids": ev(16)},
    ],
    17: [
        {"region": "both columns: notes and references 17-60", "finding": "Reference list and explanatory footnotes only. ESI pointers remain unassociated; no current synthesis, table row, or reaction arrow occurs.", "evidence_ids": ev(17)},
    ],
    18: [
        {"region": "upper half: references 61-70", "finding": "References and SVD literature notes only; no current-paper experiment, reaction arrow, or table.", "evidence_ids": ev(18)},
        {"region": "lower half blank page area", "finding": "No chemical content or experimental record.", "evidence_ids": ev(18)},
    ],
}
unresolved = {
    13: ["Complex-4 characterization includes unsatisfactory CHN values; these are retained as characterization context rather than interpreted as reaction yield/purity."],
    14: ["The printed complex-6 MeLi volume (0.27 μL) conflicts with 0.43 mmol at 1.6 M; the printed complex-9 Method-A charge (25 mg, 0.83 mmol) is also internally inconsistent. Both remain unresolved conflicts."],
    15: ["The Li[NMe2BH3] gas identity and the aminodiborane coproducts/yield are unreported; complete organozirconium structures and atom mappings remain unresolved."],
    16: ["The printed analytical statement 'reaction of 5 to 2' is potentially inconsistent with the earlier 1→6→2 discussion; it is preserved as wording and not converted into a new reaction."],
    17: ["References point to ESI and prior syntheses, but no SI was supplied or associated."],
    18: ["No experimental content is present on this page."],
}
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": page, "state": "attempted", "regions": regions[page], "record_ids": page_records[page], "unresolved": unresolved[page]}
        for page in range(13, 19)
    ],
}
write(OUT / "coverage" / "pages_0013_0018.json", coverage)

inventory = {
    "document_id": DOC,
    "pages_visually_inspected": list(range(13, 19)),
    "source_images": {str(page): f"assets/{DOC}/page{page}.jpg" for page in range(13, 19)},
    "page_inventory": [
        {"page": page, "regions": regions[page], "record_ids": page_records[page], "unresolved": unresolved[page]}
        for page in range(13, 19)
    ],
    "deduplication": [
        "Preparations of complexes 4-10 already represented from schemes/prose were enriched with pages 13-15 image evidence and detailed procedures, not duplicated.",
        "Complex 7 Methods A and B and complex 9 Methods A and B remain separate because they are distinct experimental routes.",
        "The page-15 general dehydrocoupling procedure is shared context for existing catalytic runs, not a standalone reaction.",
        "Characterization spectra, MS fragment abundances, XRD/CCDC data, DFT/SVD methods, and references remain source context.",
    ],
    "limitations": [
        "Exactly PDF pages 13-18 were visually inspected for this segment.",
        "No supporting-information file was supplied or associated.",
        "Unsafe complete organozirconium SMILES and unreported products/methods were not invented.",
    ],
}
write(HERE / "page_extraction_0013_0018.json", inventory)

# Raw standard-stereo diagnostics for every parseable candidate molecule.
stereo = {}
for molecule in M.values():
    smiles = molecule["smiles"]["normalized_value"]
    if smiles:
        stereo[molecule["molecule_id"]] = stereo_inspect(smiles)
write(HERE / "stereo_checks.json", stereo)

prepared_dataset = json.loads((OUT / "dataset.json").read_text())
mentions = source_mentions(
    prepared_dataset,
    ["Synthesis of complex 11", "Synthesis of Li[NMe", "Synthesis of BH", "General procedure for dehydrocoupling", "reaction of 5 to 2"],
    context=420,
)
write(HERE / "source_mentions_pages13_18.json", mentions)

summary = {
    "document_id": DOC,
    "molecules": len(M),
    "reactions": len(R),
    "new_reactions_pages13_18": 3,
    "specific_experiments": sum(item["record_type"] == "specific_experiment" for item in R),
    "all_needs_review": all(item["review"]["status"] == "needs_review" for item in R),
    "pages_covered_this_segment": list(range(13, 19)),
}
write(HERE / "verification_summary.json", summary)
print(json.dumps(summary, indent=2, ensure_ascii=False))
