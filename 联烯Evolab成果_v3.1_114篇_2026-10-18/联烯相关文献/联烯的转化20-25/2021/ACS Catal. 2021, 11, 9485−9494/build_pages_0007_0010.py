#!/usr/bin/env python3
"""Append source-derived v3.1 records for PDF pages 7-10 only."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / "runs/current"
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st  # noqa: E402
from reaction_tools import dataset_audit  # noqa: E402

PREV_PATH = RUN / "build_pages_0001_0006.py"
spec = importlib.util.spec_from_file_location("pages_0001_0006", PREV_PATH)
prev = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(prev)

DOC = "doc_7aafc6041d0daf10"
TITLE = "Scope and Mechanistic Investigations of Pd-Catalyzed Coupling/Cyclization and Cycloisomerization of Allenyl Malonates"
Builder = prev.Builder
sf = prev.sf


def ev(page: int, kind: str = "image") -> str:
    return f"{DOC}:page{page}:{kind}"


def conflict_fact(alternatives, evidence):
    fact = sf()
    fact.update({
        "raw_value": None,
        "normalized_value": None,
        "acquisition_method": "original_extraction",
        "evidence_ids": list(evidence),
        "data_status": "conflict",
        "verification": "unreviewed",
        "derivation": None,
        "alternatives": [sf(raw, norm, ids) for raw, norm, ids in alternatives],
    })
    return fact


def set_equiv(p, raw, number, evidence, reference="m_2a_racemic"):
    p["charges"]["equivalents"] = sf(raw, {"number": number, "unit": "equiv"}, evidence)
    p["equivalence_reference"] = sf("relative to allenyl malonate 2a", reference, evidence)


def set_molpct(p, raw, number, evidence, reference="m_2a_racemic", basis=None):
    norm = {"number": number, "unit": "mol%"}
    if basis:
        norm["basis"] = basis
    p["charges"]["mol_percent"] = sf(raw, norm, evidence)
    p["equivalence_reference"] = sf("relative to allenyl malonate 2a", reference, evidence)


def add_named_molecules(b: Builder):
    e7 = [ev(7)]
    for mid, label, name, note in [
        ("m_pdtfa2", "Pd(TFA)2", "palladium(II) trifluoroacetate", "Named Pd source in Table 3 entry 13; no complete source-verified structure encoded."),
        ("m_pdmecn2cl2", "Pd(MeCN)2Cl2", "bis(acetonitrile)palladium(II) dichloride", "Named Pd source in Table 3 entry 14; no complete source-verified structure encoded."),
        ("m_kh_paraffin_50", "KH (50% w/w in paraffin)", "potassium hydride dispersion in paraffin, 50% w/w", "Named reagent formulation from Table 3 footnote f; the active KH equivalents are row-specific and the mixture has no unique SMILES."),
        ("m_unspecified_base_page7", "base", "base used in the split-substrate mechanistic experiment", "The prose says only 'base'; its identity and charge are not restated for this experiment."),
        ("m_unspecified_catalyst_system_page7", "catalyst", "catalyst system used in the split-substrate mechanistic experiment", "The prose says only 'the catalyst'; catalyst, ligand and loadings are not restated for this experiment."),
    ]:
        b.molecule(mid, label, name, e7, identity="unresolved", chirality="not_applicable",
                   geometry="not_applicable", representation="not_applicable", note=note)

    # Isotopologues in Scheme 5 are kept separate from the unrelated compound 6
    # used for the page-4 dihydroxylation. No unchecked concrete SMILES are assigned.
    scheme5 = [
        ("m_scheme5_6", "6", "partially C1-deuterated allenyl malonate 6", "unresolved", "unresolved",
         {"type": "isotopologue mixture", "site": "C1 malonate carbon", "deuterium_incorporation": ">75%"},
         "Allene axis is present, but Scheme 5 gives no R_a/S_a assignment. This compound-number 6 is structurally distinct from the page-4 dihydroxylated product 6 and is not merged with it."),
        ("m_scheme5_7", "7", "cycloisomerized product 7", "unspecified", "not_applicable",
         {"type": "isotopologue mixture", "site": "product C4", "deuterium_incorporation": "26%"},
         "Cyclopentene product with 26% D incorporation at C4; absolute/relative configuration and isotope measurement method are not stated."),
        ("m_scheme5_8", "8", "C2-dideuterated allenyl malonate 8", "unresolved", "unresolved",
         {"type": "isotopologue", "site": "C2", "deuterium_count": 2},
         "The source drawing places two D labels at the tether carbon assigned C2 by the Scheme 5 numbering; allene axial configuration is unassigned."),
        ("m_scheme5_9", "9", "cycloisomerized product 9", "unspecified", "not_applicable",
         {"type": "isotopologue", "site": "product C2", "deuterium_count": 2},
         "The two D labels are retained at product C2; the figure explicitly notes no D at C4. Product stereochemistry is not assigned."),
        ("m_scheme5_10", "10", "C3-deuterated allenyl malonate 10", "unresolved", "unresolved",
         {"type": "isotopologue", "site": "allene carbon C3", "deuterium_count": 1},
         "The source drawing places D at allene carbon C3; allene axial configuration is unassigned."),
        ("m_scheme5_11", "11", "cycloisomerized product 11", "unspecified", "not_applicable",
         {"type": "isotopologue", "site": "product C3", "deuterium_count": 1},
         "The D label is retained at product C3; the figure explicitly notes no D at C4. Product stereochemistry is not assigned."),
        ("m_scheme5_12", "12", "C5-dideuterated allenyl malonate 12", "unresolved", "unresolved",
         {"type": "isotopologue", "site": "C5 side-chain methylene", "deuterium_count": 2},
         "The source drawing shows a C4H9-CD2 substituent at C5; allene axial configuration is unassigned."),
        ("m_scheme5_13", "13", "cycloisomerized product 13", "unspecified", "not_applicable",
         {"type": "isotopologue", "site": "C5 side-chain methylene", "deuterium_count": 2},
         "The C4H9-CD2 substituent is retained; the figure explicitly notes no D at C4. Product stereochemistry is not assigned."),
    ]
    for mid, label, name, chirality, geometry, mixture, note in scheme5:
        b.molecule(mid, label, name, [ev(8)], identity="source_identified", chirality=chirality,
                   geometry=geometry, representation="unresolved", mixture=mixture, note=note)


def add_table3(b: Builder):
    e = [ev(7)]
    rows = [
        (1, "m_pd2dba3", None, None, 10, "rt", 0, "quant", None),
        (2, "m_pd2dba3", "m_nah", 1.2, 0, "rt", 8, 40, None),
        (3, "m_pd2dba3", "m_nah", 1.2, 10, "rt", 41, 11, None),
        (4, "m_pd2dba3", "m_khmds", 1.2, 10, "rt", 18, 70, "one_portion"),
        (5, "m_pd2dba3", "m_khmds", 1.2, 10, "rt", 48, 32, "over_2h"),
        (6, "m_pd2dba3", "m_khmds", 1.2, 10, "rt", 47, 40, "over_7h"),
        (7, "m_pd2dba3", "m_kh_paraffin_50", 1.0, 10, "rt", 67, None, "kh_paraffin"),
        (8, "m_pd2dba3", "m_kh_paraffin_50", 1.2, 10, "rt", 67, None, "kh_paraffin"),
        (9, "m_pd2dba3", "m_kh_paraffin_50", 0.25, 10, "rt", 53, None, "kh_paraffin"),
        (10, "m_pd2dba3", "m_kh_paraffin_50", 1.0, 10, "rt", 0, "quant", "pre_stir_conflict"),
        (11, "m_pd2dba3", "m_kh_paraffin_50", 1.2, 5, "rt", 55, None, "kh_paraffin"),
        (12, "m_pd2dba3", "m_kh_paraffin_50", 1.2, 20, "rt", 27, 38, "kh_paraffin"),
        (13, "m_pdtfa2", "m_kh_paraffin_50", 1.0, 10, "rt", 36, None, "kh_paraffin"),
        (14, "m_pdmecn2cl2", "m_kh_paraffin_50", 1.0, 10, "rt", "trace", None, "kh_paraffin"),
        (15, "m_pd2dba3", "m_kh_paraffin_50", 1.0, 10, "0 °C", 50, 24, "kh_paraffin"),
        (16, "m_pd2dba3", "m_kh_paraffin_50", 1.0, 10, "60 °C", 32, None, "kh_paraffin"),
    ]
    ids = []
    for n, pd, base, base_eq, ligand_pct, temp, result, recovery, special in rows:
        cats = ["optimization", "control"] if n in {1, 2, 10} else ["optimization"]
        r = b.base_reaction(f"table3_entry_{n:02d}", e, "specific_experiment", cats,
                            f"Table 3, entry {n}",
                            "Optimization of cycloisomerization of 2a",
                            "Pd-catalyzed cycloisomerization of allenyl malonate to cyclopentene")
        r["general_scheme_ids"] = [f"{DOC}:scheme1_current_cycloisomerization_generic"]
        parts = [
            b.participant("substrate", "m_2a_racemic", "reactant", e),
            b.participant("pd", pd, "catalyst", e),
            b.participant("dmf", "m_dmf", "solvent", e),
            b.participant("product", "m_5", "product", e, product_role="cycloisomerization product"),
        ]
        set_molpct(parts[1], "5 mol% [Pd]", 5, e, basis="Pd")
        if ligand_pct:
            ligand = b.participant("ligand", "m_poph3", "ligand", e)
            set_molpct(ligand, f"{ligand_pct:g} mol% P(OPh)3", ligand_pct, e)
            parts.append(ligand)
        if base:
            bp = b.participant("base", base, "reagent", e)
            set_equiv(bp, f"{base_eq:g} equiv", base_eq, e)
            parts.append(bp)
        if special in {"one_portion", "over_2h", "over_7h"}:
            parts.append(b.participant("thf_carrier", "m_thf", "solvent", e))
        r["participants"] = parts
        solvents = [{"solvent": "DMF", "role": "reaction solvent"}]
        if special in {"one_portion", "over_2h", "over_7h"}:
            solvents.append({"solvent": "THF", "role": "carrier for 1 M KHMDS solution", "volume": "unreported"})
        r["solvent_conditions"]["composition"] = sf("DMF" if len(solvents) == 1 else "DMF plus KHMDS solution in THF", solvents, e)
        tnorm = {"number": float(temp.split()[0]), "unit": "°C"} if "°C" in temp else temp

        if special == "one_portion":
            b.add_step(r, e, "add 1 M KHMDS in THF in one portion", (temp, tnorm), participants=["substrate", "base", "thf_carrier"])
        elif special in {"over_2h", "over_7h"}:
            hours = 2 if special == "over_2h" else 7
            b.add_step(r, e, f"add 1 M KHMDS in THF over {hours} h", (temp, tnorm),
                       (f"{hours} h", {"number": hours, "unit": "h"}),
                       participants=["substrate", "base", "thf_carrier"])
        elif special == "pre_stir_conflict":
            step = b.add_step(r, e, "pre-stir substrate 2a with KH dispersion before catalyst exposure",
                              (temp, tnorm), participants=["substrate", "base"])
            step["time"] = conflict_fact([
                ("45 min", {"number": 45, "unit": "min"}, [ev(7, "text")]),
                ("60 min", {"number": 60, "unit": "min"}, [ev(7)]),
            ], [ev(7), ev(7, "text")])
            r["procedure"]["pretreatment"] = conflict_fact([
                ("Pre-stirring 2a with KH for 45 min", "pre-stir 2a with KH for 45 min", [ev(7, "text")]),
                ("Substrate and base pre-stirred for 60 min", "pre-stir substrate and base for 60 min", [ev(7)]),
            ], [ev(7), ev(7, "text")])
            r["notes"].append(sf("Table 3 footnote g gives 60 min, while the discussion text gives 45 min.",
                                 "source conflict: pretreatment time 45 min versus 60 min", [ev(7), ev(7, "text")]))
            r["review"]["unresolved"].append("Pretreatment duration conflicts within page 7: prose reports 45 min, Table 3 footnote g reports 60 min; neither value was selected.")

        b.add_step(r, e, "hold under Table 3 cycloisomerization conditions", (temp, tnorm),
                   ("overnight", "overnight"))
        r["procedure"]["original_text"] = sf(
            {"table_header": "5 mol% [Pd]; X mol% P(OPh)3; Y equiv base; DMF; temp; overnight", "entry": n},
            f"Table 3 entry {n}", e)
        if n == 1:
            r["notes"].append(sf("A control experiment with no PhI or base (entry 1) gave no 5",
                                 "entry 1 control omits PhI and base and gives no product 5", [ev(7, "text")]))
        if special == "one_portion":
            r["procedure"]["addition_order"] = sf("1 M in THF added in one portion", "add KHMDS solution in one portion", e)
        elif special in {"over_2h", "over_7h"}:
            hours = 2 if special == "over_2h" else 7
            r["procedure"]["addition_order"] = sf(f"1 M in THF added over {hours} h", f"add KHMDS solution over {hours} h", e)

        nmr_basis = "1H NMR yield with PhSiMe3 as an internal standard"
        if result == "trace":
            b.outcome(r, "trace_5", "trace", ["m_5"], e, "trace", "trace", nmr_basis,
                      product_role="cycloisomerization product")
        else:
            basis = "isolated yield" if n == 7 else nmr_basis
            ytype = "isolated" if n == 7 else "NMR"
            b.outcome(r, "yield_5", "yield", ["m_5"], e, f"{result}%", {"number": result, "unit": "%"}, basis,
                      product_role="cycloisomerization product", yield_type=ytype)
        if recovery is not None:
            raw = "quant recovered 2a" if recovery == "quant" else f"{recovery}% recovered 2a"
            norm = "quantitative" if recovery == "quant" else {"number": recovery, "unit": "%"}
            b.outcome(r, "recovery_2a", "recovery", ["m_2a_racemic"], e, raw, norm, nmr_basis,
                      product_role="recovered starting material")

        if result == 0:
            b.classify_none(r, e, "The row reports zero product and quantitative recovered 2a; no intended-product bond formation is assigned to this negative control.")
        else:
            b.classify_changes(r, e, prev.cycloisomerization_changes(f"t3e{n}", "5"), [
                "The drawing supports ring-forming C-C connection and net C-H redistribution, but exact atom and proton-source mapping is unresolved.",
                "Trace product in entry 14 supports only qualitative occurrence of the same net product class.",
            ])
        r["review"]["unresolved"] += [
            "Reaction scale, concentration, workup, purification, and product characterization are not reported on page 7.",
            "The Table 3 header gives overnight but no numeric duration.",
        ]
        if special in {"one_portion", "over_2h", "over_7h"}:
            r["review"]["unresolved"].append("The THF volume of the 1 M KHMDS carrier solution is not reported, so total solvent volume and concentration remain missing.")
        ids.append(r["reaction_id"])
    return ids


def add_split_substrate_experiment(b: Builder):
    e = [ev(7), ev(7, "text")]
    r = b.base_reaction("page7_split_substrate_test", e, "specific_experiment", ["mechanistic_study", "control"],
                        "Page 7 discussion below Table 3",
                        "Split-substrate test of neutral and deprotonated 2a roles",
                        "Pd-catalyzed cycloisomerization under incompletely restated standard conditions")
    r["general_scheme_ids"] = [f"{DOC}:scheme1_current_cycloisomerization_generic"]
    pre = b.participant("substrate_pre", "m_2a_racemic", "reactant", e)
    set_equiv(pre, "0.2 equiv of malonate 2a", 0.2, e)
    remaining = b.participant("substrate_remaining", "m_2a_racemic", "reactant", e)
    set_equiv(remaining, "remaining 0.8 equiv of 2a", 0.8, e)
    r["participants"] = [
        pre,
        remaining,
        b.participant("base", "m_unspecified_base_page7", "unknown", e),
        b.participant("catalyst", "m_unspecified_catalyst_system_page7", "unknown", e),
        b.participant("product", "m_5", "product", e, product_role="cycloisomerization product"),
    ]
    b.add_step(r, e, "pre-stir 0.2 equiv of 2a with base; duration not stated for this experiment",
               participants=["substrate_pre", "base"])
    b.add_step(r, e, "transfer the pre-stirred portion to a flask containing catalyst and the remaining 0.8 equiv of 2a",
               participants=["substrate_pre", "base", "catalyst", "substrate_remaining"])
    r["procedure"]["original_text"] = sf(
        "0.2 equiv of malonate 2a was pre-stirred with base and transferred to a flask containing the catalyst and the remaining 0.8 equiv of 2a",
        "split 2a into 0.2-equivalent pre-stirred portion and 0.8-equivalent catalyst-flask portion", e)
    r["procedure"]["addition_order"] = sf("pre-stirred portion transferred to catalyst plus remaining substrate",
                                              ["pre-stir 0.2 equiv 2a with base", "transfer into catalyst plus 0.8 equiv 2a"], e)
    b.outcome(r, "yield_5", "yield", ["m_5"], e, "77% yield by NMR", {"number": 77, "unit": "%"},
              "NMR; internal standard not restated for this experiment", product_role="cycloisomerization product", yield_type="NMR")
    b.outcome(r, "recovery_2a", "recovery", ["m_2a_racemic"], e, "16% of the remaining 2a", {"number": 16, "unit": "%"},
              "NMR context in the sentence; separate recovery method not explicitly restated", product_role="recovered starting material")
    b.classify_changes(r, e, prev.cycloisomerization_changes("split_test", "5"), [
        "Net C-C ring formation and C-H redistribution are visible for product 5, but exact atom/proton mapping is unresolved."
    ])
    r["review"]["unresolved"] += [
        "The prose does not restate the base identity/amount, catalyst identity/loading, ligand, solvent, temperature, total time, scale, concentration, workup, or purification for this experiment.",
        "The 16% recovery is described after the NMR yield, but a separate recovery measurement method and denominator are not explicitly stated.",
    ]
    return r["reaction_id"]


def add_glovebox_control(b: Builder):
    e = [ev(8), ev(8, "text")]
    r = b.base_reaction("page8_glovebox_dry_dmf_control", e, "specific_experiment", ["mechanistic_study", "control"],
                        "Page 8, opening paragraph",
                        "Glovebox control using freshly distilled and degassed DMF",
                        "Cycloisomerization control testing exogenous water as proton source")
    r["general_scheme_ids"] = [f"{DOC}:scheme1_current_cycloisomerization_generic"]
    r["participants"] = [
        b.participant("substrate", "m_2a_racemic", "reactant", e),
        b.participant("dmf", "m_dmf", "solvent", e),
        b.participant("product", "m_5", "product", e, product_role="cycloisomerization product"),
    ]
    r["solvent_conditions"]["composition"] = sf("freshly distilled and degassed DMF", "DMF, freshly distilled and degassed", e)
    r["procedure"]["pretreatment"] = sf("freshly distilled and degassed DMF", ["distill DMF", "degas DMF"], e)
    step = b.add_step(r, e, "carry out the cycloisomerization in a glovebox using freshly distilled and degassed DMF",
                      participants=["substrate", "dmf", "product"])
    step["atmosphere"] = sf("in a glovebox", "glovebox atmosphere; gas identity unreported", e)
    b.outcome(r, "reaction_result", "other", ["m_5"], e, "proceeded as normal", "proceeded as normal", None,
              product_role="cycloisomerization product", scope="qualitative author statement; no numeric yield or detection method")
    b.classify_changes(r, e, prev.cycloisomerization_changes("glovebox", "5"), [
        "The author states that the reaction proceeded as normal, but no numeric result or exact atom/proton mapping is provided."
    ])
    r["notes"].append(sf("proceeded as normal, ruling out this possibility",
                         "Author uses this control to argue against exogenous water as the required proton source", e))
    r["review"]["unresolved"] += [
        "Catalyst, ligand, base, loadings, temperature, time, scale, concentration, workup, purification, numeric yield, and measurement method are not restated for this control.",
        "The glovebox gas identity is not reported.",
    ]
    return r["reaction_id"]


def add_scheme5(b: Builder):
    e = [ev(8)]
    data = [
        ("6_to_7", "m_scheme5_6", "m_scheme5_7", 59, "26% D incorporation", {"number": 26, "unit": "%", "site": "C4"}),
        ("8_to_9", "m_scheme5_8", "m_scheme5_9", None, "no D noted at C4", {"qualitative": "no D noted", "site": "C4"}),
        ("10_to_11", "m_scheme5_10", "m_scheme5_11", None, "no D noted at C4", {"qualitative": "no D noted", "site": "C4"}),
        ("12_to_13", "m_scheme5_12", "m_scheme5_13", None, "no D noted at C4", {"qualitative": "no D noted", "site": "C4"}),
    ]
    ids = []
    for suffix, substrate, product, yield_n, isotope_raw, isotope_norm in data:
        r = b.base_reaction(f"scheme5_{suffix}", e, "specific_experiment", ["mechanistic_study"],
                            f"Scheme 5, arrow {suffix.replace('_to_', '→')}",
                            "Deuterium-labeling cycloisomerization",
                            "isotope-tracing cycloisomerization of allenyl malonate")
        r["general_scheme_ids"] = [f"{DOC}:scheme1_current_cycloisomerization_generic"]
        r["participants"] = [
            b.participant("substrate", substrate, "reactant", e),
            b.participant("product", product, "product", e, product_role="deuterium-labeling product"),
        ]
        b.add_step(r, e, "cycloisomerization under standard conditions denoted 'a' in Scheme 5",
                   participants=["substrate", "product"])
        r["procedure"]["original_text"] = sf("a", "standard conditions label 'a'; not expanded on page 8", e)
        if yield_n is not None:
            b.outcome(r, "yield", "yield", [product], e, f"{yield_n}%", {"number": yield_n, "unit": "%"}, None,
                      product_role="deuterium-labeling product")
        b.outcome(r, "deuterium_result", "other", [product], e, isotope_raw, isotope_norm, None,
                  product_role="deuterium-labeling product",
                  scope="deuterium incorporation at product C4; measurement method unreported")
        b.classify_changes(r, e, prev.cycloisomerization_changes(f"scheme5_{suffix}", product), [
            "The net ring-forming C-C connection is visible, but H/D transfer and exact atom mapping remain partial.",
            "Deuterium is treated as a hydrogen isotope for bond-class bookkeeping; isotope-source mapping is not claimed complete.",
        ])
        r["review"]["unresolved"] += [
            "Condition label 'a' is not expanded on page 8; Table 3 conditions were not silently inherited without an explicit page-8 association.",
            "Reaction scale, catalyst/base identities and amounts, ligand, solvent, temperature, time, workup, purification, and isotope measurement method are not reported in Scheme 5.",
            "Absolute/relative product stereochemistry is not assigned; no unchecked SMILES was created from the drawing.",
        ]
        if yield_n is None:
            r["review"]["unresolved"].append("No product yield is shown for this Scheme 5 arrow.")
        else:
            r["review"]["unresolved"].append("The 59% yield measurement basis is not reported on page 8.")
        ids.append(r["reaction_id"])
    return ids


def add_scheme6(b: Builder):
    e = [ev(8)]
    labels = ["2a", "5.1", "5.2", "5.3", "5.4", "5.5"]
    for label in labels:
        mid = "m_scheme6_" + label.replace(".", "_")
        b.molecule(mid, label, f"generic Scheme 6 entity {label}", e, identity="unresolved",
                   chirality="unresolved", geometry="unresolved", representation="unresolved",
                   note="Generic mechanistic entity containing R1/E/RO and/or Pd variables; not a concrete source-verified molecule.")
    r = b.base_reaction("scheme6_proposed_cycloisomerization_mechanism", e, "proposed_mechanism", ["mechanistic_study"],
                        "Scheme 6",
                        "Proposed mechanism for cycloisomerization of allenyl malonates",
                        "author-proposed oxidative-addition/nucleopalladation/reductive-elimination cycle with alternative hydropalladation")
    r["participants"] = [
        b.participant(f"entity_{i}", "m_scheme6_" + label.replace(".", "_"), "unknown", e)
        for i, label in enumerate(labels, 1)
    ]
    r["notes"] = [
        sf("Pd = Pd[P(OPh)3]", "Pd species represented generically as Pd[P(OPh)3]", e),
        sf("oxidative addition into the malonate C-H bond; KH deprotonation; nucleopalladation; enolate regeneration; reductive elimination",
           ["malonate C-H oxidative addition", "KH deprotonation", "nucleopalladation", "enolate regeneration", "reductive elimination"], e),
        sf("Alternative hydropalladation/isomerization through 5.5 cannot be ruled out",
           "alternative Pd(eta3-allyl) pathway remains possible", [ev(8), ev(8, "text")]),
        sf("DFT: nucleopalladation ΔG‡ = 1.7 kcal/mol, ΔGR = −26.8 kcal/mol; reductive elimination ΔG‡ = 2.2 kcal/mol, ΔGR = −30.6 kcal/mol",
           {"nucleopalladation": {"barrier": 1.7, "reaction_free_energy": -26.8, "unit": "kcal/mol"},
            "reductive_elimination": {"barrier": 2.2, "reaction_free_energy": -30.6, "unit": "kcal/mol"}},
           [ev(8, "text")]),
    ]
    b.classify_none(r, e, "Scheme 6 is an author-proposed catalytic mechanism with generic intermediates, not a separate observed net experiment.")
    r["review"]["unresolved"] = [
        "Generic R1/E/RO variables and Pd coordination prevent concrete SMILES encoding.",
        "The main-text DFT numbers cite Figure S1, but no Supporting Information file is associated in this workspace; the cited profile was not inspected.",
        "The alternative hydropalladation/Pd(eta3-allyl) path is explicitly not ruled out by the authors.",
    ]
    return r["reaction_id"]


def add_source_context(bundle):
    new = [
        {"context_id": f"{DOC}:page7_mass_balance_discussion", "document_id": DOC,
         "kind": "author_discussion", "evidence_ids": [ev(7), ev(7, "text")],
         "description": "The authors state a general 60-80% mass balance and uncharacterized byproducts for the cycloisomerization optimization; this is not assigned to a specific Table 3 row."},
        {"context_id": f"{DOC}:page9_conclusion", "document_id": DOC,
         "kind": "conclusion_context", "evidence_ids": [ev(9), ev(9, "text")],
         "description": "Conclusion summarizes coupling/cyclization, chirality transfer and ligand-switched cycloisomerization; it adds no distinct reaction arrow or row."},
        {"context_id": f"{DOC}:page9_si_availability_statement", "document_id": DOC,
         "kind": "supporting_information_statement", "evidence_ids": [ev(9), ev(9, "text")],
         "description": "The article states that SI contains characterization data, optimization tables and additional substrates/catalysts, but no SI file was explicitly supplied or associated in this workspace."},
        {"context_id": f"{DOC}:pages9_10_references", "document_id": DOC,
         "kind": "bibliography", "evidence_ids": [ev(9), ev(10), ev(9, "text"), ev(10, "text")],
         "description": "Pages 9-10 contain acknowledgments and bibliography. Cited papers are not current-paper experiments and no reaction records are created from the reference list."},
    ]
    ids = {x.get("context_id") for x in new}
    bundle["source_context"] = [x for x in bundle["source_context"] if x.get("context_id") not in ids] + new


def main():
    bundle = json.loads((RUN / "submitted.json").read_text())
    template = json.loads((RUN / "raw_model/pages_0007_0010_template.json").read_text())
    # Refresh the document/evidence registry from the required standard-init
    # template.  The prior bounded merge retained only evidence used on pages
    # 1-6, whereas this pass must legitimately cite pages 7-10.
    bundle["documents"] = template["documents"]
    bundle["documents"][0]["title"] = TITLE
    # Idempotent reruns: remove only records owned by this bounded pass.
    owned_prefixes = (
        f"{DOC}:table3_entry_", f"{DOC}:page7_split_substrate_test",
        f"{DOC}:page8_glovebox_dry_dmf_control", f"{DOC}:scheme5_",
        f"{DOC}:scheme6_proposed_cycloisomerization_mechanism",
    )
    bundle["reactions"] = [r for r in bundle["reactions"] if not r["reaction_id"].startswith(owned_prefixes)]
    b = Builder(bundle)
    b.molecules = {m["molecule_id"]: m for m in bundle["molecules"]}
    b.reactions = bundle["reactions"]

    add_named_molecules(b)
    page7_ids = add_table3(b)
    page7_ids.append(add_split_substrate_experiment(b))
    page8_ids = [add_glovebox_control(b)]
    page8_ids.extend(add_scheme5(b))
    page8_ids.append(add_scheme6(b))

    bundle["molecules"] = list(b.molecules.values())
    bundle["reactions"] = b.reactions
    add_source_context(bundle)
    (RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")

    raw = {
        "document_id": DOC,
        "page_range": [7, 10],
        "review_mode": "all four original rendered page images inspected; crops used for Table 3 and Schemes 5-6",
        "article_title": TITLE,
        "semantic_accuracy_verified": False,
        "page_inventory": {
            "7": ["Table 3 entries 1-16", "Table 3 footnotes a-g", "separate 0.2/0.8-equivalent split-substrate experiment", "cycloisomerization discussion"],
            "8": ["glovebox/fresh DMF control", "Scheme 5 four deuterium-labeling arrows", "Scheme 6 proposed mechanism", "opening of conclusions"],
            "9": ["conclusion continuation", "SI availability statement", "author/funding/acknowledgments", "references 1-4"],
            "10": ["references 5-12 only"],
        },
        "new_record_counts": {
            "total": len(page7_ids) + len(page8_ids),
            "specific_experiments": len(page7_ids) + len(page8_ids) - 1,
            "optimization_rows": 16,
            "mechanistic_or_control_experiments": 6,
            "proposed_mechanisms": 1,
        },
        "limitations": [
            "No SI file was explicitly supplied or associated.",
            "Scheme 5 condition label a is not expanded on page 8 and was not guessed from Table 3.",
            "Scheme 5 structures and isotope positions are retained from the drawing without unchecked concrete SMILES.",
            "Table 3 entry 10 has an unresolved 45 min versus 60 min source conflict.",
        ],
    }
    (RUN / "raw_model/pages_0007_0010_extraction.json").write_text(json.dumps(raw, ensure_ascii=False, indent=2) + "\n")

    coverage = {"document_id": DOC, "semantic_accuracy_verified": False, "pages": [
        {"page": 7, "state": "attempted", "regions": [
            {"region": "Table 3 reaction header and structures", "finding": "General cycloisomerization layout for racemic 2a to product 5; common Pd loading, ligand/base variables, DMF, temperature and overnight condition mapped to each row.", "evidence_ids": [ev(7), ev(7, "text")]},
            {"region": "Table 3 rows 1-16", "finding": "Sixteen distinct optimization/control experiments recorded row-by-row, including zero, trace, product yield and starting-material recovery outcomes.", "evidence_ids": [ev(7), ev(7, "text")]},
            {"region": "Table 3 footnotes a-g", "finding": "NMR/isolated-yield bases, KHMDS addition modes, KH formulation and pre-stirring footnote applied to their rows; entry 10 footnote reports 60 min.", "evidence_ids": [ev(7)]},
            {"region": "discussion columns below Table 3", "finding": "Separate 0.2/0.8-equivalent split-substrate experiment gives 77% NMR yield of 5 and 16% remaining 2a; prose reports 45 min for entry 10, conflicting with the 60 min footnote.", "evidence_ids": [ev(7), ev(7, "text")]},
        ], "record_ids": page7_ids, "unresolved": [
            "Table 3 entry 10 pretreatment duration conflicts: 45 min in prose versus 60 min in footnote g.",
            "Table 3 omits reaction scale, concentration, workup, purification and characterization; overnight is not numeric.",
            "The split-substrate experiment does not restate catalyst/base identities or most conditions.",
        ]},
        {"page": 8, "state": "attempted", "regions": [
            {"region": "opening paragraph", "finding": "Real glovebox control using freshly distilled and degassed DMF, reported qualitatively as proceeding normally; no numeric yield.", "evidence_ids": [ev(8), ev(8, "text")]},
            {"region": "Scheme 5 upper-left", "finding": "Four distinct deuterium-labeling cycloisomerization arrows (6→7, 8→9, 10→11, 12→13); 6→7 reports 59% yield and 26% product D incorporation, while the other three report no D at C4.", "evidence_ids": [ev(8)]},
            {"region": "Scheme 5 discussion", "finding": "Authors interpret partial D transfer and lack of C4 labeling; statements retained as mechanistic context/outcomes rather than absolute proof.", "evidence_ids": [ev(8), ev(8, "text")]},
            {"region": "Scheme 6 and adjacent text", "finding": "Author-proposed oxidative-addition/nucleopalladation/reductive-elimination cycle plus an alternative hydropalladation/Pd(eta3-allyl) path; recorded as proposed mechanism, not observed experiment.", "evidence_ids": [ev(8), ev(8, "text")]},
            {"region": "Conclusions opening", "finding": "Summary only; no new reaction row or arrow beyond already recorded experiments.", "evidence_ids": [ev(8), ev(8, "text")]},
        ], "record_ids": page8_ids, "unresolved": [
            "Scheme 5 condition label a is not expanded on page 8; catalyst/base/ligand/solvent/temperature/time were not silently inherited.",
            "Scheme 5 isotope measurement methods and three of four product yields are unreported.",
            "Scheme 6 cites Figure S1 in unavailable/unassociated SI; only main-text DFT values were retained.",
        ]},
        {"page": 9, "state": "attempted", "regions": [
            {"region": "conclusion continuation", "finding": "Narrative summary of previously reported coupling/cyclization, chirality transfer and cycloisomerization; no new experiment.", "evidence_ids": [ev(9), ev(9, "text")]},
            {"region": "Associated Content", "finding": "Article states SI availability and contents, but no SI file is supplied/associated in the workspace, so none is imported.", "evidence_ids": [ev(9), ev(9, "text")]},
            {"region": "author, funding, notes and acknowledgments", "finding": "Administrative/context material; no synthetic or characterization experiment record.", "evidence_ids": [ev(9), ev(9, "text")]},
            {"region": "references 1-4", "finding": "Bibliographic citations only; cited-paper experiments are not current-paper records.", "evidence_ids": [ev(9), ev(9, "text")]},
        ], "record_ids": [], "unresolved": ["No SI association is established despite the article's availability statement."]},
        {"page": 10, "state": "attempted", "regions": [
            {"region": "full two-column page", "finding": "References 5-12 only; no current-paper figure, table, reaction arrow, experimental paragraph or standalone spectrum.", "evidence_ids": [ev(10), ev(10, "text")]},
        ], "record_ids": [], "unresolved": []},
    ]}
    (RUN / "coverage/pages_0007_0010.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")
    dataset = json.loads((RUN / "dataset.json").read_text())
    mentions = st.source_mentions(dataset, ["2a", "5", "6", "7", "8", "9", "10", "11", "12", "13"], context=180)
    (RUN / "raw_model/pages_0007_0010_mentions.json").write_text(json.dumps(mentions, ensure_ascii=False, indent=2) + "\n")
    stereo = {
        "m_2a_racemic_reused_on_page7": st.stereo_inspect("CCCCCC=C=CCC(C(=O)OC)C(=O)OC"),
        "new_scheme5_structures": {
            "status": "no concrete SMILES assigned",
            "reason": "Isotope positions are visible, but complete source-checked connectivity and allene/product stereochemistry were not safely encoded from Scheme 5.",
            "tool_limitation": "The toolkit does not establish allene axial or ring relative stereochemistry from an unresolved drawing."
        },
    }
    (RUN / "raw_model/pages_0007_0010_stereo_diagnostics.json").write_text(json.dumps(stereo, ensure_ascii=False, indent=2) + "\n")
    (RUN / "audit_pages_0007_0010.json").write_text(json.dumps(dataset_audit(dataset), ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(raw["new_record_counts"], ensure_ascii=False))


if __name__ == "__main__":
    main()
