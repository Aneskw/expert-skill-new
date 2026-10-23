#!/usr/bin/env python3
"""Rebuild the current document and add source-derived records for PDF pages 7-12.

The pages 1-6 builder is executed from a fresh standard-init template so this
extension does not import compatibility records.  Only current-document source
facts are added here; mechanisms and standalone characterization remain
separate from executed experiments.
"""

from __future__ import annotations

import json
import runpy
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = ROOT / "runs" / "current"
DOC = "doc_837adcd3d5d8cb7b"

# Reconstruct the pages 1-6 state from source-specific code and the freshly
# generated standard template, then extend that in-memory registry.
ns = runpy.run_path(str(HERE / "build_submission.py"))
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
fact = ns["fact"]
stereo_inspect = ns["stereo_inspect"]
source_mentions = ns["source_mentions"]
DEHYDRO = ns["DEHYDRO"]


def et(page: int) -> str:
    return f"{DOC}:page{page}:text"


def ei(page: int) -> str:
    return f"{DOC}:page{page}:image"


def ev(*pages: int, text: bool = True, image: bool = True) -> list[str]:
    result: list[str] = []
    for page in pages:
        if text:
            result.append(et(page))
        if image and page <= 12:
            result.append(ei(page))
    return result


def add_evidence(slug: str, evidence: list[str]) -> None:
    """Extend the one current-document occurrence without changing identity."""
    m = M[slug]
    occurrence = m["occurrences"][0]
    for eid in evidence:
        if eid not in occurrence["evidence_ids"]:
            occurrence["evidence_ids"].append(eid)


def pids(parts):
    return [p["participant_id"] for p in parts]


def simple_parts(rid, evidence, specs):
    return [
        participant(
            rid,
            i,
            slug,
            role,
            evidence,
            amounts=amounts,
            reference=reference,
            product_role=product_role,
            stage=stage_name,
        )
        for i, (slug, role, amounts, reference, product_role, stage_name) in enumerate(specs, 1)
    ]


def concentration_conflict(evidence_main: list[str]):
    """Preserve the paper's page-11 1.33 M value and page-15 0.33 M implication."""
    return fact(
        raw="[H3B·NMe2H] = 1.33 M (page 11); general procedure charges 1.33 mmol in 4 mL toluene (page 15)",
        normalized=None,
        evidence=evidence_main + [et(15)],
        status="conflict",
        verification="unreviewed",
        alternatives=[
            f("[H3B·NMe2H] = 1.33 M", q(1.33, "M"), evidence_main),
            f(
                "1.33 mmol in 4 mL implies approximately 0.33 M",
                q(0.3325, "M", qualifier="approx", raw_text="1.33 mmol / 4 mL"),
                [et(15)],
                method="inference",
                derivation="Arithmetic from the printed general-procedure amount and volume; not selected over the page-11 concentration.",
            ),
        ],
    )


# Extend evidence for molecules already introduced on pages 1-6.
for slug, pages in {
    "complex1": (7, 8, 9, 10, 11, 12),
    "complex2": (9, 12),
    "complex4": (8, 9),
    "complex5": (8, 9),
    "complex6": (7, 9, 12),
    "me_li": (7, 8, 9, 12),
    "dimethylamine_borane": (7, 8, 9, 10, 11, 12),
    "d3b_nme2h": (10,),
    "li_amidoborane": (7, 8, 9, 12),
    "cyclic_diborazane": (7, 8, 9, 10, 11, 12),
    "linear_diborazane": (7, 12),
    "diaminoborane": (7, 8, 12),
    "aminoborane": (7, 9, 10, 12),
    "mu_aminodiborane": (8,),
    "hydrogen": (7, 10, 11, 12),
    "methane": (8, 9, 10),
    "bh3_thf": (8,),
    "toluene": (7, 8, 9, 10, 11, 12),
    "toluene_d8": (7,),
    "benzene": (7, 8, 9, 10),
    "diethyl_ether": (7, 10),
    "licl": (7, 9),
    "nacl": (7,),
}.items():
    add_evidence(slug, ev(*pages))


# Molecules first introduced on pages 7-12.  Complex organometallic graphs are
# retained from source drawings and names rather than forced into unsafe SMILES.
add_mol("complex3b", "postulated dinuclear zirconocene bis(amidoborane) complex", "3b", ev(8, 9), abbrev="3b",
        chirality="unspecified", geometry="unresolved", rep="not_encoded",
        identity_note="Proposed intermediate only; two alternative Zr-N and multicentre Zr-H-B binding modes are drawn.")
add_mol("complex7", "Cp2Zr(C2SiMe3)(mu-CH2SiMe3)ZrCp2", "7", ev(7, 8, 9, 10, 11, 14), abbrev="7",
        chirality="unspecified", geometry="unresolved", rep="not_encoded",
        identity_note="Dinuclear alkynyl/trimethylsilylmethyl zirconocene; source drawing and SC-XRD support a Zr(IV/II) solid-state assignment.")
add_mol("complex8", "Cp2Zr(CH2SiMe3)(C2SiMe3)", "8", ev(7, 14), abbrev="8",
        chirality="unspecified", geometry="unresolved", rep="not_encoded")
add_mol("complex9", "Cp2Zr(C2SiMe3)(mu-H)ZrCp2", "9", ev(7, 8, 9, 10, 11, 12, 15), abbrev="9",
        chirality="unspecified", geometry="unresolved", rep="not_encoded",
        identity_note="Dinuclear hydride complex assigned by spectroscopy/literature and independently prepared; exact organometallic graph is not encoded.")
add_mol("complex9_d", "deuteride isotopologue of complex 9", "9-D", ev(10), abbrev="9-D",
        chirality="unspecified", geometry="unresolved", rep="not_encoded")
add_mol("complex10", "Cp2Zr(C2SiMe3)(mu-Me)ZrCp2", "10", ev(10, 11, 12, 15), abbrev="10",
        chirality="unspecified", geometry="unresolved", rep="not_encoded",
        identity_note="Methyl-bridged dinuclear zirconocene; source drawing and SC-XRD support a Zr(III/III) solid-state assignment.")
add_mol("complex11", "Cp2Zr(Me)(C2SiMe3)", "11", ev(10, 15), abbrev="11",
        chirality="unspecified", geometry="unresolved", rep="not_encoded")
add_mol("generic_zrr_zrh", "generic dinuclear zirconocene hydride ZrR-ZrH", "ZrR-ZrH", ev(9, 10),
        chirality="unresolved", geometry="unresolved", rep="not_encoded",
        identity_note="Generic proposed active-species family with R = H or H3BNMe2; not a single isolated molecular identity.")
add_mol("cp2zr_cl_alkynyl", "Cp2Zr(Cl)(C2SiMe3)", "Cp2Zr(Cl)(C2SiMe3)", ev(7, 14),
        chirality="unspecified", geometry="unresolved", rep="not_encoded")
add_mol("cp2zr_butadiene", "Cp2Zr(eta4-butadiene)", "[Cp2Zr]", ev(7, 14), abbrev="[Cp2Zr]",
        chirality="unspecified", geometry="unresolved", rep="not_encoded")
add_mol("cp2zr_py_bistmsacetylene", "Cp2Zr(py)(eta2-Me3SiC2SiMe3)", "Cp2Zr(py)(eta2-Me3SiC2SiMe3)", ev(7, 8, 10, 15),
        chirality="unspecified", geometry="unresolved", rep="not_encoded")
add_mol("li_ch2sime3", "lithium trimethylsilylmethyl", "LiCH2SiMe3", ev(7, 14),
        chirality="achiral", geometry="not_applicable", rep="not_encoded")
add_mol("lih", "lithium hydride", "LiH", ev(7), smiles="[Li+].[H-]",
        chirality="achiral", geometry="not_applicable", rep="encoded")
add_mol("lialh4", "lithium aluminium hydride", "LiAlH4", ev(7),
        chirality="unresolved", geometry="unresolved", rep="not_encoded",
        identity_note="The formula is source-explicit, but the ionic/coordination graph is not forced into SMILES because the toolkit treats the tetrahedral Al centre as potential stereochemistry.")
add_mol("nabeth3h", "sodium triethylborohydride", "NaBEt3H", ev(7, 14),
        chirality="achiral", geometry="not_applicable", rep="not_encoded")
add_mol("mehb_nme2", "methyl(dimethylamino)borane", "MeHB=NMe2", ev(8),
        chirality="unspecified", geometry="unspecified", rep="not_encoded")
add_mol("me2b_nme2", "dimethyl(dimethylamino)borane", "Me2B=NMe2", ev(8),
        chirality="unspecified", geometry="unspecified", rep="not_encoded")
add_mol("sime4", "tetramethylsilane", "SiMe4", ev(10), smiles="C[Si](C)(C)C",
        chirality="achiral", geometry="not_applicable", rep="encoded")
add_mol("bistmsacetylene", "bis(trimethylsilyl)acetylene", "Me3SiC2SiMe3", ev(10),
        smiles="C[Si](C)(C)C#C[Si](C)(C)C", chirality="achiral", geometry="not_applicable", rep="encoded")
add_mol("pyridine", "pyridine", "py", ev(7, 10), smiles="n1ccccc1", abbrev="py",
        chirality="achiral", geometry="not_applicable", rep="encoded")
add_mol("butadiene", "1,3-butadiene", "butadiene", ev(7), smiles="C=CC=C",
        chirality="achiral", geometry="unspecified", rep="encoded")
add_mol("pentane", "pentane", "pentane", ev(14, 15), smiles="CCCCC",
        chirality="achiral", geometry="not_applicable", rep="encoded")
add_mol("thf", "tetrahydrofuran", "THF", ev(15), smiles="C1CCOC1",
        chirality="achiral", geometry="not_applicable", rep="encoded")
add_mol("benzene_d6", "benzene-d6", "benzene-d6", ev(7, 8, 9, 10, 15),
        chirality="achiral", geometry="not_applicable", rep="not_encoded")


# Page 7, Fig. 6: one time-resolved in situ experiment, distinct from the
# post-reaction spectrum in Fig. 1 and from the cited H3B.NMeH2 spectrum.
e7 = ev(7)
rid = f"{DOC}:in_situ_nmr_dehydrocoupling"
parts = simple_parts(rid, e7, [
    ("dimethylamine_borane", "reactant", {}, None, None, "reaction"),
    ("complex1", "catalyst", {"mol_percent": q(2.5, "mol%")}, "dimethylamine_borane", None, "reaction"),
    ("me_li", "reagent", {"mol_percent": q(5.2, "mol%")}, "dimethylamine_borane", None, "reaction"),
    ("toluene_d8", "solvent", {}, None, None, "reaction"),
    ("linear_diborazane", "product", {}, None, "dominant species at low conversion", "reaction"),
    ("cyclic_diborazane", "product", {}, None, "main product at late stage", "reaction"),
    ("aminoborane", "product", {}, None, "very-low-concentration early species", "reaction"),
    ("diaminoborane", "product", {}, None, "very-low-concentration early species", "reaction"),
    ("hydrogen", "product", {}, None, "gas coproduct", "reaction"),
    ("complex7", "product", {}, None, "tentative assignment of a high-field signal", "characterization"),
    ("complex9", "product", {}, None, "tentative hydride-species assignment", "characterization"),
])
st = [step(rid, 1, e7, pids(parts), "initiate dehydrocoupling and acquire time-resolved 11B and 1H{11B} NMR spectra", q(25, "°C"), q(8, "h"))]
outs = [
    outcome(rid, 1, "other", ["linear_diborazane"], e7, "dominant at low conversion", basis="time-resolved 11B NMR (96.3 MHz)"),
    outcome(rid, 2, "other", ["cyclic_diborazane"], e7, "large amounts only at a very late stage", basis="time-resolved 11B NMR (96.3 MHz)"),
    outcome(rid, 3, "trace", ["aminoborane", "diaminoborane"], e7, "very low concentration at the onset", basis="time-resolved 11B NMR (96.3 MHz)"),
    outcome(rid, 4, "other", ["complex7", "complex9"], e7, "high-field singlets at -5.12 and -5.89 ppm; assignments are tentative", basis="time-resolved 1H{11B} NMR (400.1 MHz)"),
]
reaction("in_situ_nmr_dehydrocoupling", "specific_experiment", e7, ["mechanistic_study"],
         "in situ NMR monitoring of dehydrocoupling", "amine-borane dehydrocoupling and catalyst speciation",
         parts, st, outs, general=("scheme1_current_dehydrocoupling",), class_changes=DEHYDRO,
         notes=["The Fig. 6 bottom H3B·NMeH2 spectrum is explicitly previous work and is not counted as a current-paper experiment."],
         unresolved=["The -5.12/-5.89 ppm resonances are possible Zr-hydride/speciation signals, not unambiguous product identities."])


# Hydride-source attempts and the three Scheme 3 preparation arrows.
for slug, reagent_slug, reagent_name in [
    ("hydride_source_lih", "lih", "LiH"),
    ("hydride_source_lialh4", "lialh4", "LiAlH4"),
]:
    rid = f"{DOC}:{slug}"
    parts = simple_parts(rid, e7, [
        ("complex1", "reactant", {}, None, None, "reaction"),
        (reagent_slug, "reagent", {}, None, None, "reaction"),
    ])
    st = [step(rid, 1, e7, pids(parts), f"add {reagent_name} to complex 1")]
    outs = [outcome(rid, 1, "decomposition", ["complex1"], e7, "ill-defined mixture", basis="source discussion; analytical method not separately stated")]
    reaction(slug, "specific_experiment", e7, ["mechanistic_study", "control"],
             f"attempted hydride synthesis using {reagent_name}", "hydride-source reaction",
             parts, st, outs, class_coverage="not_applicable",
             unresolved=["Reagent amount, solvent, temperature, time, and mixture identities are not reported on the supplied pages."])

e7_14 = ev(7) + [et(14)]
rid = f"{DOC}:synthesis_complex7_method_a"
parts = simple_parts(rid, e7_14, [
    ("complex1", "reactant", {"mass": q(300, "mg"), "amount": q(0.43, "mmol")}, None, None, "reaction"),
    ("nabeth3h", "reagent", {"volume": q(1.08, "mL"), "amount": q(1.08, "mmol")}, "complex1", None, "reaction"),
    ("toluene", "solvent", {"volume": q(20, "mL")}, None, None, "reaction"),
    ("pentane", "purification_agent", {"volume": q(6, "mL", raw_text="3 x 2 mL")}, None, None, "purification"),
    ("complex7", "product", {}, None, "isolated product", "reaction"),
    ("nacl", "product", {}, None, "salt coproduct shown in Scheme 3a", "reaction"),
])
st = [
    step(rid, 1, e7_14, pids(parts[:3]), "cool the toluene solution of 1", q(-78, "°C")),
    step(rid, 2, e7_14, pids(parts[:3]), "add NaBEt3H dropwise", q(-78, "°C"), q(10, "min")),
    step(rid, 3, e7_14, pids(parts[:3]), "stir after addition", None, "overnight"),
    step(rid, 4, e7_14, pids(parts), "cannula-filter and concentrate until precipitation", stage="workup"),
    step(rid, 5, e7_14, pids(parts), "filter red crystals and wash with pentane 3 x 2 mL", stage="purification"),
]
outs = [
    outcome(rid, 1, "product_amount", ["complex7"], e7_14, q(78, "mg"), basis="isolated dark-red solid"),
    outcome(rid, 2, "yield", ["complex7"], e7_14, q(26, "%"), basis="isolated dark-red solid", yield_type="isolated"),
]
reaction("synthesis_complex7_method_a", "specific_experiment", e7_14, ["mechanistic_study", "other"],
         "synthesis of complex 7, method A", "hydride-induced allenediide activation",
         parts, st, outs, class_changes=None,
         notes=["Page 14 procedure text is used as matching context; page 14 is not claimed as visually covered."],
         unresolved=["The detailed carbon/zirconium/hydrogen atom mapping for the rearrangement is unresolved."])

rid = f"{DOC}:synthesis_complex8"
parts = simple_parts(rid, e7_14, [
    ("cp2zr_cl_alkynyl", "reactant", {"mass": q(64, "mg"), "amount": q(0.18, "mmol")}, None, None, "reaction"),
    ("li_ch2sime3", "reagent", {"mass": q(19, "mg"), "amount": q(0.18, "mmol")}, "cp2zr_cl_alkynyl", None, "reaction"),
    ("diethyl_ether", "solvent", {"volume": q(10, "mL")}, None, None, "reaction"),
    ("pentane", "purification_agent", {"volume": q(1.5, "mL", raw_text="3 x 0.5 mL")}, None, None, "purification"),
    ("complex8", "product", {}, None, "isolated product", "reaction"),
    ("licl", "product", {}, None, "salt coproduct", "reaction"),
])
st = [
    step(rid, 1, e7_14, pids(parts[:3]), "combine the two solids with pre-cooled Et2O under exclusion of light", q(-78, "°C")),
    step(rid, 2, e7_14, pids(parts[:3]), "stir at -78 °C", q(-78, "°C"), "overnight"),
    step(rid, 3, e7_14, pids(parts), "remove Et2O, wash the solid with pentane 3 x 0.5 mL at -78 °C, and filter", q(-78, "°C"), stage="workup"),
    step(rid, 4, e7_14, pids(parts), "store at -78 °C to obtain colourless crystals", q(-78, "°C"), stage="purification"),
]
outs = [
    outcome(rid, 1, "product_amount", ["complex8"], e7_14, q(53, "mg"), basis="isolated colourless crystals"),
    outcome(rid, 2, "yield", ["complex8"], e7_14, q(83, "%"), basis="isolated colourless crystals", yield_type="isolated"),
]
reaction("synthesis_complex8", "specific_experiment", e7_14, ["other"], "synthesis of complex 8",
         "alkyl substitution at zirconium", parts, st, outs,
         class_changes=[("C", "Zr", 0, 1, "formed", "A Zr-CH2SiMe3 bond is formed."),
                        ("Cl", "Zr", 1, 0, "broken", "The Zr-Cl bond is displaced with LiCl formation.")],
         notes=["Page 14 procedure text is used as matching context; page 14 is not claimed as visually covered."],
         unresolved=["Complete organometallic atom mapping is not asserted."])

rid = f"{DOC}:synthesis_complex7_method_b"
parts = simple_parts(rid, e7_14, [
    ("complex8", "reactant", {"mass": q(16, "mg"), "amount": q(0.04, "mmol")}, None, None, "reaction"),
    ("cp2zr_butadiene", "reactant", {"mass": q(10.8, "mg"), "amount": q(0.04, "mmol")}, "complex8", None, "reaction"),
    ("benzene", "solvent", {"volume": q(0.6, "mL")}, None, None, "reaction"),
    ("complex7", "product", {}, None, "isolated product", "reaction"),
    ("butadiene", "product", {}, None, "coproduct shown in Scheme 3b", "reaction"),
])
st = [
    step(rid, 1, e7_14, pids(parts[:3]), "dissolve the two precursors in benzene and stir", q(75, "°C"), q(3, "d")),
    step(rid, 2, e7_14, pids(parts), "remove benzene in vacuo and dry the red solid", None, q(2, "h"), stage="workup"),
]
outs = [
    outcome(rid, 1, "product_amount", ["complex7"], e7_14, q(12, "mg"), basis="isolated red solid"),
    outcome(rid, 2, "yield", ["complex7"], e7_14, q(75, "%"), basis="isolated red solid", yield_type="isolated"),
]
reaction("synthesis_complex7_method_b", "specific_experiment", e7_14, ["other"], "synthesis of complex 7, method B",
         "dinuclear zirconocene assembly", parts, st, outs, class_changes=None,
         notes=["The paper calls this route poorly reproducible; page 14 procedure text is context and is not claimed as visually covered."],
         unresolved=["Exact changed bonds and the fate/amount of butadiene are not atom-mapped."])

e78 = ev(7, 8)
rid = f"{DOC}:failed_alkyne_complex_route_to7"
parts = simple_parts(rid, e78, [
    ("complex8", "reactant", {}, None, None, "reaction"),
    ("cp2zr_py_bistmsacetylene", "reactant", {}, None, None, "reaction"),
    ("complex7", "product", {}, None, "candidate product explicitly reported as not obtained", "reaction"),
])
st = [step(rid, 1, e78, pids(parts), "attempt alternative preparation of complex 7 using Cp2Zr(py)(eta2-Me3SiC2SiMe3)")]
outs = [outcome(rid, 1, "no_reaction", ["complex7"], e78, "attempts were not successful", basis="source statement; detection method not reported")]
reaction("failed_alkyne_complex_route_to7", "specific_experiment", e78, ["control", "mechanistic_study"],
         "unsuccessful alternative route to complex 7", "alternative dinuclear-complex assembly",
         parts, st, outs, class_coverage="not_applicable",
         unresolved=["The sentence does not fully state the co-reactant/conditions; complex 8 is retained as the likely sequence substrate but requires expert review."])


# Page 8: catalytic UV-vis and stoichiometric speciation experiments.
e8 = ev(8)
rid = f"{DOC}:complex4_bh3thf_reaction"
parts = simple_parts(rid, e8, [
    ("complex4", "reactant", {}, None, None, "reaction"),
    ("bh3_thf", "reagent", {}, None, None, "reaction"),
    ("complex5", "product", {}, None, "only identified complex", "reaction"),
    ("complex7", "product", {}, None, "candidate product not observed", "characterization"),
    ("complex9", "product", {}, None, "candidate hydride product not observed", "characterization"),
])
st = [step(rid, 1, e8, pids(parts), "react complex 4 with BH3·THF and analyze high-field resonances")]
outs = [
    outcome(rid, 1, "other", ["complex5"], e8, "only complex 5 identified, along with several unidentified metallocene species", basis="NMR spectroscopy (Fig. S66)"),
    outcome(rid, 2, "no_reaction", ["complex7", "complex9"], e8, "high-field resonances for 7 or a Zr hydride were not observed", basis="NMR spectroscopy (Fig. S66)"),
]
reaction("complex4_bh3thf_reaction", "specific_experiment", e8, ["mechanistic_study", "control"],
         "reaction of complex 4 with BH3·THF", "borate formation control",
         parts, st, outs, class_changes=None,
         unresolved=["Amounts, solvent, temperature, time, and unidentified metallocene products are in unavailable SI Fig. S66."])

rid = f"{DOC}:catalytic_uvvis_monitoring"
parts = simple_parts(rid, e8, [
    ("dimethylamine_borane", "reactant", {}, None, None, "reaction"),
    ("complex1", "catalyst", {"mol_percent": q(2.5, "mol%")}, "dimethylamine_borane", None, "reaction"),
    ("me_li", "reagent", {"mol_percent": q(5.2, "mol%")}, "dimethylamine_borane", None, "reaction"),
    ("hydrogen", "product", {}, None, "gas coproduct", "reaction"),
])
st = [step(rid, 1, e8, pids(parts), "monitor catalytic dehydrocoupling by operando UV-vis", None, q(480, "min"))]
outs = [outcome(rid, 1, "other", ["complex1"], e8,
                "at least three UV-vis-active species; maxima at 385 and 625 nm evolve between 0 and 480 min; reaction is not kinetically uniform",
                basis="UV-vis spectroscopy, singular-value decomposition, and extinction diagram")]
reaction("catalytic_uvvis_monitoring", "specific_experiment", e8, ["mechanistic_study"],
         "UV-vis monitoring of catalytic dehydrocoupling", "catalyst speciation during dehydrocoupling",
         parts, st, outs, general=("scheme1_current_dehydrocoupling",), class_changes=DEHYDRO,
         unresolved=["Solvent, concentration, and exact species concentrations are in unavailable SI Fig. S83; UV-vis cannot provide reliable kinetics here."])

rid = f"{DOC}:stoichiometric_1_meli_dmab"
parts = simple_parts(rid, e8, [
    ("complex1", "reactant", {}, None, None, "reaction"),
    ("me_li", "reagent", {}, "complex1", None, "reaction"),
    ("dimethylamine_borane", "reactant", {"equivalents": q(1, "equiv")}, "complex1", None, "reaction"),
    ("benzene_d6", "solvent", {}, None, None, "reaction"),
    ("complex7", "product", {}, None, "dominant zirconium species", "reaction"),
    ("methane", "product", {}, None, "observed coproduct", "reaction"),
    ("cyclic_diborazane", "product", {}, None, "dehydrocoupling product", "reaction"),
    ("diaminoborane", "product", {}, None, "dehydrocoupling product", "reaction"),
    ("mehb_nme2", "product", {}, None, "methyl-scrambling product", "reaction"),
    ("me2b_nme2", "product", {}, None, "methyl-scrambling product", "reaction"),
])
st = [step(rid, 1, e8, pids(parts), "react 1 with MeLi and equimolar H3B·NMe2H in benzene-d6; analyze by 1H and 11B NMR")]
outs = [
    outcome(rid, 1, "other", ["complex7"], e8, "dominant zirconium species", basis="1H NMR (Fig. S68)"),
    outcome(rid, 2, "other", ["methane"], e8, "observed at delta 0.16 ppm", basis="1H NMR (Fig. S68)"),
    outcome(rid, 3, "other", ["cyclic_diborazane", "diaminoborane", "mehb_nme2", "me2b_nme2"], e8,
            "dehydrocoupling and methyl-scrambling products observed", basis="11B NMR (Fig. S69)"),
]
reaction("stoichiometric_1_meli_dmab", "specific_experiment", e8, ["mechanistic_study"],
         "stoichiometric reaction of 1/MeLi with H3B·NMe2H", "precatalyst activation and methyl scrambling",
         parts, st, outs, class_changes=None,
         unresolved=["The MeLi amount and full stoichiometric recipe are not stated on page 8; exact product quantities are unreported."])

for slug, equivalents, primary, result, basis in [
    ("stoichiometric_1_liamidoborane_2eq", 2, "complex7",
     "blue-green mixture; signal at -5.09 ppm is dominant; cyclic diborazane, typical dehydrocoupling products, mu-aminodiborane and Me2B=NMe2 observed",
     "1H/11B NMR (Figs. S62/S64)"),
    ("stoichiometric_1_liamidoborane_4eq", 4, "complex9",
     "the -5.09 ppm signal is absent and a new -7.21 ppm resonance assigned to complex 9 appears",
     "1H NMR (Fig. S62)"),
]:
    rid = f"{DOC}:{slug}"
    parts = simple_parts(rid, e8, [
        ("complex1", "reactant", {}, None, None, "reaction"),
        ("li_amidoborane", "reagent", {"equivalents": q(equivalents, "equiv")}, "complex1", None, "reaction"),
        (primary, "product", {}, None, "assigned dominant/diagnostic zirconium species", "reaction"),
        ("cyclic_diborazane", "product", {}, None, "observed dehydrocoupling product", "reaction"),
        ("mu_aminodiborane", "product", {}, None, "observed minor boron species", "reaction"),
        ("me2b_nme2", "product", {}, None, "observed methylated boron species", "reaction"),
    ])
    st = [step(rid, 1, e8, pids(parts), f"react complex 1 with {equivalents} equivalents Li[NMe2BH3] and analyze by NMR")]
    outs = [outcome(rid, 1, "other", [primary, "cyclic_diborazane", "mu_aminodiborane", "me2b_nme2"], e8, result, basis=basis)]
    reaction(slug, "specific_experiment", e8, ["mechanistic_study"],
             f"stoichiometric 1/Li[NMe2BH3] ({equivalents} equivalents)", "amidoborane activation and zirconium speciation",
             parts, st, outs, class_changes=None,
             unresolved=["Exact quantities, solvent, time, temperature, and product amounts are not stated on page 8."])

rid = f"{DOC}:variable_temperature_1_liamidoborane"
parts = simple_parts(rid, e8, [
    ("complex1", "reactant", {}, None, None, "reaction"),
    ("li_amidoborane", "reagent", {}, "complex1", None, "reaction"),
    ("complex7", "product", {}, None, "assigned product", "reaction"),
])
st = [
    step(rid, 1, e8, pids(parts), "analyze 1/Li[NMe2BH3] by variable-temperature 1H and 11B NMR", q(-60, "°C")),
    step(rid, 2, e8, pids(parts), "warm the NMR sample to room temperature and compare resonances", "room temperature"),
]
outs = [outcome(rid, 1, "other", ["complex7"], e8,
                "7 and an unassigned ca. -6.0 ppm signal form at low temperature; 11B broad singlet at 4 ppm at -60 °C resolves into a triplet at room temperature",
                basis="variable-temperature 1H/11B NMR")]
reaction("variable_temperature_1_liamidoborane", "specific_experiment", e8, ["mechanistic_study"],
         "variable-temperature NMR of 1/Li[NMe2BH3]", "low-temperature catalyst speciation",
         parts, st, outs, class_changes=None,
         unresolved=["Li[NMe2BH3] ratio, solvent, hold times, and identities of the -6.0/-5.0 ppm species are unresolved."])


# Page 9 controls plus two explicitly proposed speciation mechanisms.
e9 = ev(9)
rid = f"{DOC}:complex6_liamidoborane"
parts = simple_parts(rid, e9, [
    ("complex6", "reactant", {}, None, None, "reaction"),
    ("li_amidoborane", "reagent", {}, "complex6", None, "reaction"),
    ("complex7", "product", {}, None, "trace identified species", "reaction"),
])
st = [step(rid, 1, e9, pids(parts), "react complex 6 with Li[NMe2BH3] and analyze by NMR")]
outs = [outcome(rid, 1, "trace", ["complex7"], e9,
                "complex mixture containing traces of 7 and unidentified resonances at -5.82 and -6.33 ppm",
                basis="NMR spectroscopy (Fig. S74)")]
reaction("complex6_liamidoborane", "specific_experiment", e9, ["mechanistic_study"],
         "reaction of complex 6 with Li[NMe2BH3]", "monomethyl-complex activation",
         parts, st, outs, class_changes=None,
         unresolved=["Amounts, conditions, unidentified products, and exact conversion are in unavailable SI Fig. S74."])

rid = f"{DOC}:complex2_dmab_stoichiometric"
parts = simple_parts(rid, e9, [
    ("complex2", "reactant", {}, None, None, "reaction"),
    ("dimethylamine_borane", "reactant", {}, "complex2", None, "reaction"),
    ("complex7", "product", {}, None, "candidate product not observed", "characterization"),
    ("complex9", "product", {}, None, "candidate product not observed", "characterization"),
])
st = [step(rid, 1, e9, pids(parts), "perform stoichiometric reaction of isolated 2 with H3B·NMe2H and analyze by NMR")]
outs = [outcome(rid, 1, "no_reaction", ["complex7", "complex9"], e9,
                "no high-field resonances corresponding to complexes 7 or 9", basis="NMR spectroscopy (Fig. S72)")]
reaction("complex2_dmab_stoichiometric", "specific_experiment", e9, ["mechanistic_study", "control"],
         "stoichiometric reaction of isolated 2 with H3B·NMe2H", "dimethyl-complex activation control",
         parts, st, outs, class_coverage="not_applicable",
         unresolved=["Other products, amounts, solvent, temperature, and time are not reported on page 9."])

rid = f"{DOC}:scheme4_speciation_summary"
parts = simple_parts(rid, e9, [
    ("complex1", "reactant", {}, None, None, "reaction"),
    ("me_li", "reagent", {}, None, None, "reaction"),
    ("li_amidoborane", "reagent", {}, None, None, "reaction"),
    ("complex2", "product", {}, None, "experimentally observed branch species", "reaction"),
    ("complex6", "product", {}, None, "experimentally observed branch species", "reaction"),
    ("complex7", "product", {}, None, "experimentally confirmed species", "reaction"),
    ("complex9", "product", {}, None, "experimentally confirmed species", "reaction"),
    ("complex3b", "product", {}, None, "postulated intermediate", "reaction"),
])
reaction("scheme4_speciation_summary", "proposed_mechanism", e9, ["mechanistic_study"],
         "catalyst speciation summary", "interconversion of experimentally observed and postulated zirconocene species",
         parts, [], [], class_changes=None,
         notes=["Scheme 4 deliberately mixes experimentally confirmed and postulated species; its arrows are not counted as executed experiments."],
         unresolved=["Several arrows and the ZrR-ZrH family are postulated and are not assigned quantitative conditions or yields."])

rid = f"{DOC}:scheme5_amidoborane_hydride_mechanism"
parts = simple_parts(rid, e9, [
    ("complex1", "reactant", {}, None, None, "reaction"),
    ("li_amidoborane", "reagent", {"equivalents": q(2, "equiv")}, "complex1", None, "reaction"),
    ("toluene", "solvent", {}, None, None, "reaction"),
    ("complex3b", "product", {}, None, "postulated amidoborane intermediate", "reaction"),
    ("generic_zrr_zrh", "product", {}, None, "proposed active hydride family", "reaction"),
    ("aminoborane", "product", {}, None, "proposed eliminated product", "reaction"),
    ("licl", "product", {}, None, "salt coproduct", "reaction"),
])
st = [step(rid, 1, e9, pids(parts), "proposed reaction of 1 with two equivalents Li[NMe2BH3]", q(-78, "°C"))]
reaction("scheme5_amidoborane_hydride_mechanism", "proposed_mechanism", e9, ["mechanistic_study"],
         "proposed formation of Zr amidoborane and hydride complexes", "amidoborane-mediated precatalyst activation",
         parts, st, [], class_changes=None,
         unresolved=["Intermediate 3b and the elimination sequence are proposed; no isolated yield or complete mapping is reported."])


# Page 10, Fig. 7/Table 1 and Scheme 6 plus isotope/hydrogen-transfer controls.
e10 = ev(10)
rid = f"{DOC}:fig7_formal_hydride_cycles"
parts = simple_parts(rid, e10, [
    ("generic_zrr_zrh", "catalyst", {}, None, None, "reaction"),
    ("dimethylamine_borane", "reactant", {}, None, None, "reaction"),
    ("aminoborane", "product", {}, None, "proposed released intermediate", "reaction"),
    ("hydrogen", "product", {}, None, "gas coproduct", "reaction"),
    ("cyclic_diborazane", "product", {}, None, "off-metal dimerization product", "reaction"),
])
outs = [outcome(rid, 1, "other", ["hydrogen", "aminoborane"], e10,
                "formal catalytic cycles have overall DeltaRG = -13.0 kcal mol-1; Table 1 reports seven calculated steps for R = H and R = H3BNMe2",
                basis="DFT calculation (Fig. 7 and Table 1)")]
reaction("fig7_formal_hydride_cycles", "proposed_mechanism", e10, ["mechanistic_study"],
         "proposed dehydrocoupling mechanism with ZrR-ZrH", "DFT catalytic cycle",
         parts, [], outs, class_changes=DEHYDRO,
         unresolved=["The cycle is computational/proposed; intermediates Int1-Int4 and TS1/TS2 are not isolated experiments."])

e10_15 = ev(10) + [et(15)]
rid = f"{DOC}:synthesis_complex10"
parts = simple_parts(rid, e10_15, [
    ("complex11", "reactant", {"mass": q(240, "mg"), "amount": q(0.74, "mmol")}, None, None, "reaction"),
    ("cp2zr_py_bistmsacetylene", "reactant", {"mass": q(350, "mg"), "amount": q(0.74, "mmol")}, "complex11", None, "reaction"),
    ("benzene", "solvent", {"volume": q(50, "mL")}, None, None, "reaction"),
    ("pentane", "purification_agent", {"volume": q(20, "mL")}, None, None, "purification"),
    ("thf", "purification_agent", {"volume": q(0.1, "mL")}, None, None, "purification"),
    ("complex10", "product", {}, None, "isolated product", "reaction"),
    ("bistmsacetylene", "product", {}, None, "coproduct shown in Scheme 6", "reaction"),
    ("pyridine", "product", {}, None, "coproduct shown in Scheme 6", "reaction"),
])
st = [
    step(rid, 1, e10_15, pids(parts[:3]), "dissolve the two complexes in benzene and stir", q(60, "°C"), "overnight"),
    step(rid, 2, e10_15, pids(parts), "remove volatiles, dissolve the red solid in pentane and THF, and slowly cool", q(-78, "°C"), stage="workup"),
    step(rid, 3, e10_15, pids(parts), "filter off the orange supernatant at -78 °C and dry the red residue", q(-78, "°C"), stage="purification"),
]
outs = [
    outcome(rid, 1, "product_amount", ["complex10"], e10_15, q(130, "mg", raw_text="0.13 g"), basis="isolated red crystalline precipitate"),
    outcome(rid, 2, "yield", ["complex10"], e10_15, q(33, "%"), basis="isolated red crystalline precipitate", yield_type="isolated"),
]
reaction("synthesis_complex10", "specific_experiment", e10_15, ["other"], "synthesis of complex 10",
         "dinuclear zirconocene assembly", parts, st, outs,
         class_changes=[("C", "Zr", 0, 1, "formed", "A bridging methyl carbon forms a connection to the second Zr centre.")],
         notes=["Page 15 procedure text is used as matching context; page 15 is not claimed as visually covered."],
         unresolved=["Exact migration atom mapping and coproduct quantities are unresolved."])

rid = f"{DOC}:synthesis_complex9_from10"
parts = simple_parts(rid, e10_15, [
    ("complex10", "reactant", {"mass": q(133, "mg"), "amount": q(0.24, "mmol")}, None, None, "reaction"),
    ("dimethylamine_borane", "reactant", {"mass": q(600, "mg"), "amount": q(10.18, "mmol")}, "complex10", None, "reaction"),
    ("benzene", "solvent", {"volume": q(15, "mL")}, None, None, "reaction"),
    ("pentane", "purification_agent", {"volume": q(12, "mL", raw_text="4 x 3 mL")}, None, None, "purification"),
    ("complex9", "product", {}, None, "isolated product", "reaction"),
    ("methane", "product", {}, None, "coproduct shown in Scheme 6", "reaction"),
])
st = [
    step(rid, 1, e10_15, pids(parts[:3]), "dissolve complex 10 and excess H3B·NMe2H in benzene and stir", None, "overnight"),
    step(rid, 2, e10_15, pids(parts), "remove volatiles and dry the residue in a warm water bath in vacuo", "warm water bath", stage="workup"),
    step(rid, 3, e10_15, pids(parts), "wash with pentane 4 x 3 mL and keep the product cold", None, stage="purification"),
]
outs = [
    outcome(rid, 1, "product_amount", ["complex9"], e10_15, q(97, "mg"), basis="isolated blue-green solid"),
    outcome(rid, 2, "yield", ["complex9"], e10_15, q(72, "%"), basis="isolated blue-green solid", yield_type="isolated"),
    outcome(rid, 3, "decomposition", ["complex9"], e10_15, "decomposition starts immediately during drying or room-temperature storage", basis="source observation"),
]
reaction("synthesis_complex9_from10", "specific_experiment", e10_15, ["mechanistic_study", "other"],
         "synthesis of complex 9, method B", "conversion of methyl bridge to hydride bridge",
         parts, st, outs,
         class_changes=[("H", "Zr", 0, 1, "formed", "A bridging Zr-H interaction is formed in complex 9."),
                        ("C", "H", 0, 1, "formed", "Methane formation creates C-H bond(s) from the methyl bridge.")],
         notes=["Page 15 procedure text is used as matching context; page 15 is not claimed as visually covered."],
         unresolved=["Hydrogen atom origins and complete atom mapping are not established; complex 9 is highly unstable."])

rid = f"{DOC}:d3b_transfer_probe_1_meli"
parts = simple_parts(rid, e10, [
    ("d3b_nme2h", "reactant", {}, None, None, "reaction"),
    ("complex1", "catalyst", {}, None, None, "reaction"),
    ("me_li", "reagent", {}, None, None, "reaction"),
    ("complex7", "product", {}, None, "candidate isotopologue signal not observed", "characterization"),
])
st = [step(rid, 1, e10, pids(parts), "dehydrocouple D3B·NMe2H with 1/MeLi and inspect the high-field 1H NMR signal")]
outs = [outcome(rid, 1, "other", ["complex7"], e10,
                "no signal at -5.09 ppm", basis="1H NMR; used by authors to infer transfer from the borane moiety to the allenediide ligand")]
reaction("d3b_transfer_probe_1_meli", "specific_experiment", e10, ["mechanistic_study"],
         "B-deuterated hydrogen-transfer probe", "isotope-tracing dehydrocoupling",
         parts, st, outs, class_changes=None,
         notes=["This NMR tracing experiment is separate from the page-4 volumetric D3B·NMe2H comparison and is therefore not merged with that record."],
         unresolved=["Exact catalyst loading, solvent, temperature, time, and deuterated product identity are not reported on page 10."])

e10_15 = ev(10) + [et(15)]
printed_complex7_amount = fact(
    raw="25 mg, 0.83 mmol",
    normalized=None,
    evidence=[et(15)],
    status="conflict",
    verification="unreviewed",
    alternatives=[
        f("25 mg; 0.83 mmol as printed", {"mass": q(25, "mg"), "amount": q(0.83, "mmol")}, [et(15)]),
        f("25 mg at the reported M = 626 implies about 0.040 mmol", q(0.040, "mmol", qualifier="approx", raw_text="25 mg / 626 g mol-1"), [et(14), et(15)],
          method="inference", derivation="Arithmetic consistency check only; not selected over the printed amount."),
    ],
)
rid = f"{DOC}:synthesis_complex9_from7"
parts = simple_parts(rid, e10_15, [
    ("complex7", "reactant", {"amount": {"_fact": printed_complex7_amount}}, None, None, "reaction"),
    ("dimethylamine_borane", "reactant", {"mass": q(45, "mg"), "amount": q(0.75, "mmol")}, "complex7", None, "reaction"),
    ("benzene_d6", "solvent", {"volume": q(0.6, "mL")}, None, None, "reaction"),
    ("complex9", "product", {}, None, "NMR-confirmed product", "reaction"),
    ("sime4", "product", {}, None, "observed coproduct", "reaction"),
    ("cyclic_diborazane", "product", {}, None, "observed dehydrocoupling product", "reaction"),
])
# Attach the separately reported mass while leaving the amount conflict intact.
parts[0]["charges"]["mass"] = f("25 mg", q(25, "mg"), [et(15)])
st = [step(rid, 1, e10_15, pids(parts), "dissolve 7 and H3B·NMe2H in benzene-d6 in a Young NMR tube, close immediately, and hold", None, q(1, "d"))]
st[0]["pressure"] = f("closed tube; H2 pressure allowed to build", "closed tube; H2 pressure allowed to build", e10_15)
outs = [outcome(rid, 1, "other", ["complex9", "sime4", "cyclic_diborazane"], e10_15,
                "solution changes from dark red to blue; 9 confirmed by NMR, with SiMe4 and cyclic diborazane observed",
                basis="NMR spectroscopy (Fig. S30 and method A)")]
reaction("synthesis_complex9_from7", "specific_experiment", e10_15, ["mechanistic_study", "other"],
         "synthesis of complex 9, method A", "hydrogen transfer from amine borane to complex 7",
         parts, st, outs,
         class_changes=[("H", "Zr", 0, 1, "formed", "A bridging Zr-H interaction is formed in complex 9."),
                        ("C", "H", 0, 1, "formed", "Hydrogenation of the CH2SiMe3 fragment gives SiMe4.")],
         notes=["The page-15 printed charge '25 mg, 0.83 mmol' is internally inconsistent with the reported molecular ion/molar mass; both the printed amount and arithmetic alternative are retained with no correction selected.",
                "Full conversion requires large excess amine borane or a closed flask that permits H2 pressure buildup."],
         unresolved=["Printed complex-7 amount conflict remains unresolved.", "The complete hydrogen atom mapping is not established."])

rid = f"{DOC}:deuterated_transfer_7_to9d"
parts = simple_parts(rid, e10, [
    ("complex7", "reactant", {}, None, None, "reaction"),
    ("d3b_nme2h", "reactant", {}, "complex7", None, "reaction"),
    ("complex9_d", "product", {}, None, "deuteride product", "reaction"),
])
st = [step(rid, 1, e10, pids(parts), "react 7 with D3B·NMe2H and analyze the hydride region by 1H NMR")]
outs = [outcome(rid, 1, "other", ["complex9_d"], e10,
                "the product 9-D shows no resonance at -7.21 ppm", basis="1H NMR; authors assign B-D to Zr-D transfer")]
reaction("deuterated_transfer_7_to9d", "specific_experiment", e10, ["mechanistic_study"],
         "deuterium-transfer experiment from D3B·NMe2H to 7", "B-D to Zr-D transfer probe",
         parts, st, outs, class_changes=None,
         unresolved=["Amounts, solvent, temperature, time, deuterium incorporation level, and full isotopomer mapping are unreported."])

rid = f"{DOC}:hydrogenation_complex10_to9"
parts = simple_parts(rid, e10, [
    ("complex10", "reactant", {}, None, None, "reaction"),
    ("hydrogen", "reagent", {}, None, None, "reaction"),
    ("complex9", "product", {}, None, "slowly formed product", "reaction"),
])
st = [step(rid, 1, e10, pids(parts), "pressurize a J. Young NMR tube containing 10 with H2 and monitor conversion", None, None)]
st[0]["pressure"] = f(None, None, e10, status="unreported")
outs = [outcome(rid, 1, "other", ["complex9"], e10, "complex 9 slowly formed", basis="NMR monitoring")]
reaction("hydrogenation_complex10_to9", "specific_experiment", e10, ["mechanistic_study", "control"],
         "hydrogenation of complex 10 to complex 9", "methyl-bridge hydrogenolysis",
         parts, st, outs,
         class_changes=[("H", "Zr", 0, 1, "formed", "A bridging Zr-H interaction is formed in complex 9."),
                        ("C", "H", 0, 1, "formed", "Hydrogenolysis of the methyl bridge is expected to form C-H bond(s).")],
         unresolved=["H2 pressure, solvent, temperature, time, methane observation, and conversion are not reported."])


# Pages 11-12, Fig. 9: three separate catalyst runs.  The parent 1/MeLi curve
# already exists as catalyst_loading_2_50 and receives only repeated evidence.
for eid in ev(11):
    baseline = next(r for r in R if r["reaction_id"].endswith("catalyst_loading_2_50"))
    if eid not in baseline["evidence_ids"]:
        baseline["evidence_ids"].append(eid)

cat_defs = [
    ("catalysis_complex7", "complex7", q(0.2, "equiv", qualifier="approx", raw_text="approximately 0.2 equivalents after 140 hours"), q(140, "h"),
     None, None, "very slow H2 evolution; no significant contribution to parent-system activity"),
    ("catalysis_complex9", "complex9", q(0.8, "equiv", raw_text="0.8 equivalents within approximately two days"), q(2, "d", qualifier="approx", raw_text="approximately two days"),
     "cyclic_diborazane", "main product", "no induction period; crude catalyst approximately 90% pure; residual substrate remains"),
    ("catalysis_complex10", "complex10", "close to full conversion", None,
     "cyclic_diborazane", "main product", "pronounced induction period followed by steady H2 evolution; traces of diaminoborane"),
]
for slug, catalyst, gas_value, duration, borazane, borazane_role, note in cat_defs:
    evidence = ev(11, 12) + [et(15)] if slug != "catalysis_complex7" else ev(11) + [et(15)]
    rid = f"{DOC}:{slug}"
    specs = [
        ("dimethylamine_borane", "reactant", {"mass": q(78.4, "mg"), "amount": q(1.33, "mmol")}, None, None, "reaction"),
        (catalyst, "catalyst", {"mol_percent": q(2.5, "mol%")}, "dimethylamine_borane", None, "reaction"),
        ("toluene", "solvent", {"volume": q(4, "mL")}, None, None, "reaction"),
        ("hydrogen", "product", {}, None, "measured gas product", "reaction"),
    ]
    if borazane:
        specs.append((borazane, "product", {}, None, borazane_role, "reaction"))
    if slug == "catalysis_complex10":
        specs.append(("diaminoborane", "product", {}, None, "trace product", "reaction"))
    parts = simple_parts(rid, evidence, specs)
    st = [step(rid, 1, evidence, pids(parts), "perform isobaric dehydrocoupling and monitor H2 evolution", q(25, "°C"), duration, "Ar")]
    outs = [outcome(rid, 1, "product_amount" if isinstance(gas_value, dict) else "conversion", ["hydrogen"] if isinstance(gas_value, dict) else ["dimethylamine_borane"],
                    evidence, gas_value, basis="volumetric gas evolution" if isinstance(gas_value, dict) else "volumetric gas evolution; qualitative close-to-full statement")]
    if borazane:
        outs.append(outcome(rid, 2, "other", [borazane], evidence, borazane_role, basis="11B NMR at maximum conversion"))
    if slug == "catalysis_complex9":
        outs.append(outcome(rid, 3, "recovery", ["dimethylamine_borane"], evidence, "traces of residual substrate", basis="11B NMR at maximum conversion"))
    if slug == "catalysis_complex10":
        outs.append(outcome(rid, 3, "trace", ["diaminoborane"], evidence, "traces", basis="11B NMR (Fig. S57)"))
    r = reaction(slug, "specific_experiment", evidence, ["mechanistic_study"],
                 f"dehydrocoupling with complex {catalyst.replace('complex', '')}", "amine-borane dehydrocoupling",
                 parts, st, outs, general=("scheme1_current_dehydrocoupling",),
                 class_changes=DEHYDRO if borazane else [("H", "H", 0, 1, "formed", "Net H-H formation is reported by H2 evolution.")],
                 notes=[note, "Page 11 reports [H3B·NMe2H] = 1.33 M, whereas the page-15 general procedure charges 1.33 mmol in 4 mL (about 0.33 M); both are retained as a conflict."],
                 unresolved=["The 1.33 M versus 1.33 mmol/4 mL concentration conflict is unresolved.", "Individual borazane yields are not reported."])
    r["reaction_scale"] = f("H3B·NMe2H 1.33 mmol", {"reference_molecule_id": M["dimethylamine_borane"]["molecule_id"], "amount": q(1.33, "mmol")}, evidence)
    r["solvent_conditions"]["composition"] = f("toluene", [M["toluene"]["molecule_id"]], evidence)
    r["solvent_conditions"]["total_volume"] = f("4 mL", q(4, "mL"), [et(15)], method="inherited_general_condition")
    r["solvent_conditions"]["concentration"] = concentration_conflict(ev(11))
    r["solvent_conditions"]["concentration_basis"] = f("[H3B·NMe2H]", M["dimethylamine_borane"]["molecule_id"], evidence)


e12 = ev(12)
rid = f"{DOC}:fig10_complex9_cycle"
parts = simple_parts(rid, e12, [
    ("complex9", "catalyst", {}, None, None, "reaction"),
    ("dimethylamine_borane", "reactant", {}, None, None, "reaction"),
    ("aminoborane", "product", {}, None, "proposed released intermediate", "reaction"),
    ("hydrogen", "product", {}, None, "gas coproduct", "reaction"),
    ("cyclic_diborazane", "product", {}, None, "assumed off-metal dimerization product", "reaction"),
])
outs = [outcome(rid, 1, "other", ["hydrogen", "aminoborane"], e12,
                "calculated cycle DeltaRG = -12.2 kcal mol-1; reported barriers 10.8 and 9.2 kcal mol-1",
                basis="DFT calculation (Fig. 10)")]
reaction("fig10_complex9_cycle", "proposed_mechanism", e12, ["mechanistic_study"],
         "proposed dehydrocoupling mechanism with complex 9", "DFT catalytic cycle",
         parts, [], outs, class_changes=DEHYDRO,
         unresolved=["The cycle is proposed/computational; Int1-Int3 and TS1-TS3 are not isolated experiments and B-N coupling is assumed off-metal."])


# Build the final v3 bundle.  The clean standard-init document/evidence lineage
# is retained, while only source-derived molecules/reactions/context are added.
bundle = json.loads((HERE / "standard_init.json").read_text())
bundle["molecules"] = list(M.values())
bundle["reactions"] = R
bundle["source_context"] = list(json.loads((OUT / "submitted.json").read_text()).get("source_context", []))
bundle["source_context"].extend([
    {"context_id": f"{DOC}:context_pages7_12_spectroscopy", "document_id": DOC, "evidence_ids": ev(7, 8, 9, 10, 11, 12),
     "kind": "spectroscopy_and_structure_assignment",
     "note": "NMR peak assignments, MS fragments, EPR non-resonance, HMBC correlations, UV-vis maxima, SC-XRD structures, and bond distances support records but are not independent syntheses."},
    {"context_id": f"{DOC}:context_pages8_12_dft", "document_id": DOC, "evidence_ids": ev(8, 9, 10, 11, 12),
     "kind": "computational_mechanistic_context",
     "note": "DFT pathways, Table 1 thermodynamic values, migration barriers, and Figures 7/10 are computational/proposed; executed chemical experiments are recorded separately."},
    {"context_id": f"{DOC}:context_pages7_12_negative_spectral_inventory", "document_id": DOC, "evidence_ids": ev(7, 8, 9, 10, 11, 12),
     "kind": "negative_inventory",
     "note": "No standalone mass-spectrum peak table, molecular-weight table, or fragment-relative-abundance block on pages 7-12 was converted into a synthesis record; the prose MS fragments for 7 remain characterization context."},
])
write(OUT / "submitted.json", bundle)


# Required page-bounded coverage for exactly pages 7-12.
page_records = {page: [] for page in range(7, 13)}
for record in R:
    for page in page_records:
        if et(page) in record["evidence_ids"] or ei(page) in record["evidence_ids"]:
            page_records[page].append(record["reaction_id"])

regions = {
    7: [
        {"region": "top-left continuation and Fig. 6", "finding": "Li[NMe2BH3] activation conclusion plus a distinct time-resolved 11B/1H NMR dehydrocoupling/speciation experiment; cited H3B·NMeH2 spectrum is prior work.", "evidence_ids": ev(7)},
        {"region": "right-column hydride-source prose", "finding": "LiH and LiAlH4 give separate ill-defined-mixture attempts; NaBEt3H affords complex 7 in low yield.", "evidence_ids": ev(7)},
        {"region": "Scheme 3a and 3b", "finding": "Concrete arrows for 1→7, Cp2Zr(Cl)(C2SiMe3)→8, and 8+[Cp2Zr]→7; matching page-14 procedures supply quantities without extending visual coverage.", "evidence_ids": ev(7)},
    ],
    8: [
        {"region": "left-column computational/complex-4 discussion", "finding": "DFT routes to 7 are proposed context; the real 4+BH3·THF control gives 5 plus unidentified metallocenes and no 7/hydride signal.", "evidence_ids": ev(8)},
        {"region": "left-column UV-vis section", "finding": "One operando catalytic UV-vis mechanistic experiment (2.5 mol% 1, 5.2 mol% MeLi, 0-480 min).", "evidence_ids": ev(8)},
        {"region": "right-column stoichiometric reactions", "finding": "Separate 1/MeLi/DMAB experiment and 1/Li[NMe2BH3] experiments at 1:2 and 1:4 ratios; VT-NMR follow-up retained separately.", "evidence_ids": ev(8)},
    ],
    9: [
        {"region": "Scheme 4", "finding": "Mixed experimental/postulated catalyst-speciation summary retained as proposed mechanism, not split into fictitious executed arrows.", "evidence_ids": ev(9)},
        {"region": "lower-left prose and Scheme 5", "finding": "Real 6+Li[NMe2BH3] and isolated-2+DMAB controls separated from the proposed amidoborane/hydride sequence.", "evidence_ids": ev(9)},
        {"region": "right-column DFT discussion", "finding": "Computational amidoborane/hydride interpretation only; no additional executed experiment counted.", "evidence_ids": ev(9)},
    ],
    10: [
        {"region": "Fig. 7 and Table 1", "finding": "Proposed formal ZrR-ZrH cycles and calculated thermodynamic values; one proposed-mechanism record, not seven experiments.", "evidence_ids": ev(10)},
        {"region": "Scheme 6 and improved-synthesis prose", "finding": "Concrete synthesis of 10 and conversion of 10 to 9; page-15 procedures supply charges/yields and the source instability note.", "evidence_ids": ev(10)},
        {"region": "right-column control-experiment paragraph", "finding": "B-deuterated 1/MeLi NMR probe, 7+DMAB→9, 7+D3B·NMe2H→9-D, and H2 pressurization of 10 are distinct real experiments.", "evidence_ids": ev(10)},
    ],
    11: [
        {"region": "left-column EPR/HMBC/DFT and Fig. 8", "finding": "Characterization and bond-assignment context for 7 and 10; no synthesis record created from spectra, bond distances, or X-ray images.", "evidence_ids": ev(11)},
        {"region": "right-column catalytic section and Fig. 9", "finding": "Three separate 2.5 mol% catalyst runs using 7, crude 9, and 10; parent 1/MeLi curve is an existing duplicate and only gains evidence.", "evidence_ids": ev(11)},
    ],
    12: [
        {"region": "top-left continuation of Fig. 9 discussion", "finding": "Completes the complex-10 catalytic outcome: close to full conversion, cyclic diborazane, and trace diaminoborane.", "evidence_ids": ev(12)},
        {"region": "Fig. 10 and DFT text", "finding": "Proposed catalytic cycle for 9 retained as a proposed mechanism; computed intermediates/transition states are not experiments.", "evidence_ids": ev(12)},
        {"region": "Conclusions", "finding": "Summary of already recorded experiments and author interpretation; no additional experiment created.", "evidence_ids": ev(12)},
    ],
}
unresolved = {
    7: ["LiH/LiAlH4 conditions and mixture identities are not stated; Scheme 3 procedures are completed using extracted page-14 text only."],
    8: ["Most SI figure details and stoichiometric product quantities are unavailable; proposed intermediates are not treated as isolated species."],
    9: ["Scheme 4 explicitly mixes confirmed and postulated species; exact atom mapping and conditions for several arrows remain unresolved."],
    10: ["The printed Method-A charge for complex 7 (25 mg, 0.83 mmol) is internally inconsistent and is preserved as a conflict; isotope-control quantities are mostly unreported."],
    11: ["Catalytic concentration conflicts with the page-15 general-procedure amount/volume; exact individual borazane yields are absent."],
    12: ["Fig. 10 is computational/proposed and does not establish isolated intermediates; complex-10 endpoint is qualitative close-to-full conversion."],
}
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": page, "state": "attempted", "regions": regions[page], "record_ids": page_records[page], "unresolved": unresolved[page]}
        for page in range(7, 13)
    ],
}
write(OUT / "coverage" / "pages_0007_0012.json", coverage)

inventory = {
    "document_id": DOC,
    "pages_visually_inspected": list(range(7, 13)),
    "source_images": {str(p): f"assets/{DOC}/page{p}.jpg" for p in range(7, 13)},
    "page_inventory": [
        {"page": page, "regions": regions[page], "record_ids": page_records[page], "unresolved": unresolved[page]}
        for page in range(7, 13)
    ],
    "deduplication": [
        "Fig. 9 parent 1/MeLi curve was linked to catalyst_loading_2_50 rather than duplicated.",
        "The Fig. 6 current H3B·NMe2H time series is distinct from Fig. 1 post-reaction NMR; the bottom H3B·NMeH2 spectrum is cited prior work.",
        "Schemes 4/5 and Figures 7/10 are proposed/mechanistic summaries, not executed arrow-by-arrow experiments.",
        "EPR, HMBC, XRD, MS fragments, DFT energies, and bond distances remain source context rather than synthesis records.",
    ],
    "limitations": [
        "Only PDF pages 7-12 were visually inspected for this batch segment.",
        "Extracted text from pages 14-15 was used only as matching procedural context; those pages are not claimed as visually covered.",
        "ESI was advertised but not supplied or associated.",
    ],
}
write(HERE / "page_extraction_0007_0012.json", inventory)

# Save all parseable graph diagnostics as raw machine output and literal mention
# hints for cross-page review; neither establishes structure correctness.
stereo = {}
for molecule in M.values():
    smi = molecule["smiles"]["normalized_value"]
    if smi:
        stereo[molecule["molecule_id"]] = stereo_inspect(smi)
write(HERE / "stereo_checks.json", stereo)

prepared_dataset = json.loads((OUT / "dataset.json").read_text())
mentions = source_mentions(
    prepared_dataset,
    ["complex 7", "complex 9", "complex 10", "Li[NMe2BH3]", "NaBEt3H", "D3B·NMe2H", "H2 pressure"],
    context=360,
)
write(HERE / "source_mentions_pages7_12.json", mentions)

summary = {
    "document_id": DOC,
    "molecules": len(M),
    "reactions": len(R),
    "new_reactions_pages7_12": len(R) - 28,
    "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in R),
    "proposed_mechanisms": sum(r["record_type"] == "proposed_mechanism" for r in R),
    "all_needs_review": all(r["review"]["status"] == "needs_review" for r in R),
    "pages_covered_this_segment": list(range(7, 13)),
}
write(HERE / "verification_summary.json", summary)
print(json.dumps(summary, indent=2, ensure_ascii=False))
