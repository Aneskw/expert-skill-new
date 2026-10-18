#!/usr/bin/env python3
"""Extend the v3.1 bundle from source PDF pages 7--12.

The page images were inspected for this batch.  Page 9 contains three new
preparative experiments.  Methods already represented by page-1--6 records are
deduplicated and gain page-9 image evidence instead of being recreated.
"""

from __future__ import annotations

import copy
import importlib.util
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))

import standard_records as st  # noqa: E402


PREVIOUS_PATH = ROOT / "runs/current/raw_model/build_pages_0001_0006.py"
SPEC = importlib.util.spec_from_file_location("pages_0001_0006_helpers", PREVIOUS_PATH)
PREVIOUS = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(PREVIOUS)

q = PREVIOUS.q
f = PREVIOUS.f
missing = PREVIOUS.missing
not_applicable = PREVIOUS.not_applicable
concrete_molecule = PREVIOUS.concrete_molecule
unresolved_molecule = PREVIOUS.unresolved_molecule
participant = PREVIOUS.participant
step = PREVIOUS.step
outcome = PREVIOUS.outcome
reaction_base = PREVIOUS.reaction_base
set_manual_classification = PREVIOUS.set_manual_classification


DOC = "doc_ca1481599a4b0f91"
P = {n: f"{DOC}:page{n}:text" for n in range(1, 13)}
I = {n: f"{DOC}:page{n}:image" for n in range(1, 13)}
PAGE9 = [P[9], I[9]]


def explicit_z_molecule(mid, label, name, smiles, evidence, note):
    """Named oleate/oleic-acid identity with encoded Z geometry."""
    m = concrete_molecule(mid, label, name, smiles, evidence, identity_note=note)
    diagnostic = st.stereo_inspect(smiles)
    canonical = diagnostic["canonical_isomeric_smiles"]
    m["stereochemistry"]["double_bonds"] = f(
        {"source_identity": name},
        [{"bond_index": 8, "assignment": "Z", "basis": "named oleate/oleic-acid identity"}],
        evidence,
        method="structure_derivation",
        derivation="The candidate encodes the cis (Z) C9=C10 geometry conventionally specified by the named oleate/oleic-acid identity; the graph remains machine-transcribed.",
    )
    m["stereochemistry"]["assignment"] = f(
        {"source_identity": name},
        "Z geometry at the olefin encoded in the candidate isomeric SMILES",
        evidence,
        method="structure_derivation",
        derivation="Common-name structure derivation only; no independent standard or spectrum was used to certify the molecular identity.",
    )
    m["stereochemistry"]["geometry_status"] = f(
        {"source_identity": name}, "specified", evidence,
        method="structure_derivation",
        derivation="The candidate graph contains one specified Z alkene; the carbonyl is not E/Z-applicable.",
    )
    m["stereochemistry"]["representation_status"] = f(
        {"candidate_smiles": canonical}, "encoded", evidence,
        method="structure_derivation",
        derivation="The candidate isomeric SMILES retains the named olefin geometry.",
    )
    return m


def extend_page9_image_evidence(bundle):
    """Add inspected page-9 image lineage wherever page-9 text was already cited."""
    def walk(value):
        if isinstance(value, dict):
            evidence = value.get("evidence_ids")
            if isinstance(evidence, list) and P[9] in evidence and I[9] not in evidence:
                evidence.append(I[9])
            for child in value.values():
                walk(child)
        elif isinstance(value, list):
            for child in value:
                walk(child)

    # Prepared document metadata and evidence hashes are immutable.  Only the
    # authored record payload may gain the inspected image citation.
    walk(bundle["molecules"])
    walk(bundle["reactions"])
    walk(bundle["source_context"])
    for reaction in bundle["reactions"]:
        reaction["review"]["unresolved"] = [
            item for item in reaction["review"]["unresolved"]
            if "Page 9 was used as methods context but is outside" not in item
        ]


def set_inferred_ramp_time(reaction_step):
    reaction_step["time"] = f(
        "heating rate 3.3 °C/min from 110 °C to 300 °C",
        q(57.6, "min", raw_text="~57.6 min calculated from (300-110)/3.3", qualifier="approx"),
        PAGE9,
        method="inference",
        derivation="Calculated ramp duration: (300 °C - 110 °C) / (3.3 °C/min) = 57.6 min. The paper reports the rate, not a rounded duration.",
    )


def add_molecules(bundle):
    existing = {m["molecule_id"] for m in bundle["molecules"]}
    molecules = [
        concrete_molecule(
            f"{DOC}:mol_3_butyn_1_ol", "3-butyn-1-ol", "3-butyn-1-ol",
            "C#CCCO", PAGE9,
        ),
        concrete_molecule(
            f"{DOC}:mol_nitric_acid", "70% HNO3", "nitric acid",
            "O=[N+]([O-])O", PAGE9,
            identity_note="The source reports 70% HNO3. Candidate SMILES represents nitric acid only; solution composition beyond the stated 70% is not expanded.",
        ),
        unresolved_molecule(
            f"{DOC}:mol_sodium_dichromate_dihydrate", "Na2Cr2O7·2H2O",
            "sodium dichromate dihydrate", PAGE9,
            chirality="not_applicable", geometry="not_applicable", representation="unresolved",
            note="The hydrated inorganic salt is source-identified, but a unique explicit ionic/hydration SMILES was not forced from the formula.",
        ),
        unresolved_molecule(
            f"{DOC}:mol_sodium_periodate", "NaIO4", "sodium periodate", PAGE9,
            chirality="not_applicable", geometry="not_applicable", representation="unresolved",
            note="The inorganic salt is source-identified by formula; a unique explicit ionic SMILES was not forced.",
        ),
        concrete_molecule(
            f"{DOC}:mol_diethyl_ether", "diethyl ether", "diethyl ether",
            "CCOCC", PAGE9,
        ),
        unresolved_molecule(
            f"{DOC}:mol_magnesium_sulfate", "anhydrous magnesium sulfate",
            "anhydrous magnesium sulfate", PAGE9,
            chirality="not_applicable", geometry="not_applicable", representation="unresolved",
            note="The drying salt is explicitly named, but a disconnected ionic SMILES was not needed for this extraction and was left unresolved.",
        ),
        concrete_molecule(
            f"{DOC}:mol_pentane", "pentane", "pentane", "CCCCC", PAGE9,
        ),
        unresolved_molecule(
            f"{DOC}:mol_iron_iii_chloride_hexahydrate", "iron(III) chloride hexahydrate",
            "iron(III) chloride hexahydrate", PAGE9,
            chirality="not_applicable", geometry="not_applicable", representation="unresolved",
            note="The hydrated inorganic salt is explicitly named, but its coordination and hydration state are not represented as one concrete SMILES.",
        ),
        explicit_z_molecule(
            f"{DOC}:mol_sodium_oleate", "sodium oleate", "sodium oleate",
            "CCCCCCCC/C=C\\CCCCCCCC(=O)[O-].[Na+]", PAGE9,
            "Candidate isomeric SMILES encodes the named oleate cis double bond; source identity and machine-transcribed graph remain pending review.",
        ),
        concrete_molecule(
            f"{DOC}:mol_ethanol", "ethanol", "ethanol", "CCO", PAGE9,
        ),
        unresolved_molecule(
            f"{DOC}:mol_iron_oleate", "iron oleate", "iron oleate precursor", PAGE9,
            chirality="not_applicable", geometry="unresolved", representation="unresolved",
            note="The source names an iron-oleate precursor but does not give a unique stoichiometry, coordination structure, molecular formula, or stereochemically complete molecular graph.",
        ),
        explicit_z_molecule(
            f"{DOC}:mol_oleic_acid", "oleic acid", "oleic acid",
            "CCCCCCCC/C=C\\CCCCCCCC(=O)O", PAGE9,
            "Candidate isomeric SMILES encodes the named oleic-acid cis double bond; source identity and machine-transcribed graph remain pending review.",
        ),
        concrete_molecule(
            f"{DOC}:mol_1_octadecene", "1-octadecene", "1-octadecene",
            "C=CCCCCCCCCCCCCCCCC", PAGE9,
        ),
        concrete_molecule(
            f"{DOC}:mol_1_tetradecene", "1-tetradecene", "1-tetradecene",
            "C=CCCCCCCCCCCCC", PAGE9,
        ),
        unresolved_molecule(
            f"{DOC}:mol_fe3o4_nanocrystals", "8 nm magnetite (Fe3O4) nanocrystals",
            "8 nm magnetite (Fe3O4) nanocrystals", [P[2], P[3], I[3], P[9], I[9]],
            chirality="not_applicable", geometry="not_applicable", representation="unresolved",
            note="The nanocrystal product has a bulk Fe3O4 identity and reported particle size, but no unique finite molecular graph or stereochemical identity.",
        ),
    ]
    bundle["molecules"].extend(m for m in molecules if m["molecule_id"] not in existing)


def reaction_3_butynoic_acid():
    ev = PAGE9
    alcohol = f"{DOC}:mol_3_butyn_1_ol"
    nitric = f"{DOC}:mol_nitric_acid"
    dichromate = f"{DOC}:mol_sodium_dichromate_dihydrate"
    periodate = f"{DOC}:mol_sodium_periodate"
    water = f"{DOC}:mol_water"
    ether = f"{DOC}:mol_diethyl_ether"
    mgso4 = f"{DOC}:mol_magnesium_sulfate"
    pentane = f"{DOC}:mol_pentane"
    acid = f"{DOC}:mol_3_butynoic_acid"

    r = reaction_base(
        "rxn_3_butynoic_acid_synthesis", ev,
        {"page": 9, "section": "Methods—Synthesis of 3-butynoic acid"},
        "Synthesis of 3-butynoic acid",
        "oxidation of 3-butyn-1-ol to 3-butynoic acid",
    )
    r["participants"] = [
        participant("water_initial", water, "solvent", ev, volume=("deionized water (135 mL)", q(135, "mL"))),
        participant("nitric_acid", nitric, "reagent", ev,
                    volume=("70% HNO3 (0.170 mL)", q(0.170, "mL")),
                    amount=("2.5 mmol", q(2.5, "mmol")),
                    equivalents=("0.05 eq.", q(0.05, "equiv")),
                    reference=("relative to 3-butyn-1-ol", alcohol)),
        participant("sodium_dichromate", dichromate, "reagent", ev,
                    mass=("0.15 g", q(0.15, "g")), amount=("0.5 mmol", q(0.5, "mmol")),
                    equivalents=("0.01 eq.", q(0.01, "equiv")),
                    reference=("relative to 3-butyn-1-ol", alcohol)),
        participant("sodium_periodate", periodate, "reagent", ev,
                    mass=("34.2 g", q(34.2, "g")), amount=("160 mmol", q(160, "mmol")),
                    equivalents=("3.2 eq.", q(3.2, "equiv")),
                    reference=("relative to 3-butyn-1-ol", alcohol)),
        participant("butynol", alcohol, "reactant", ev,
                    mass=("3.50 g", q(3.50, "g")), amount=("50 mmol", q(50, "mmol")),
                    equivalents=("1 eq.", q(1, "equiv")),
                    reference=("3-butyn-1-ol reference", alcohol)),
        participant("water_alcohol_solution", water, "solvent", ev,
                    volume=("chilled deionized water (135 mL)", q(135, "mL"))),
        participant("diethyl_ether", ether, "workup_agent", ev, stage="workup",
                    volume=("diethyl ether (80 mL × 5)", q(400, "mL", raw_text="80 mL × 5"))),
        participant("magnesium_sulfate", mgso4, "workup_agent", ev, stage="workup"),
        participant("pentane_wash", pentane, "purification_agent", ev, stage="purification",
                    volume=("pentane (200 mL × 3)", q(600, "mL", raw_text="200 mL × 3"))),
        participant("acid_product", acid, "product", ev, product_role="main product"),
    ]
    r["reaction_scale"] = f(
        "3-butyn-1-ol: 50 mmol",
        {"reference_molecule_id": alcohol, "amount": q(50, "mmol")}, ev,
    )
    r["solvent_conditions"]["composition"] = f(
        "135 mL deionized water in the flask; 3-butyn-1-ol diluted with a separate 135 mL chilled deionized water",
        [{"molecule_id": water, "volume": q(135, "mL"), "use": "initial charge"},
         {"molecule_id": water, "volume": q(135, "mL"), "use": "alcohol solution"}], ev,
    )
    r["solvent_conditions"]["total_volume"] = missing(ev)
    r["steps"] = [
        step("acid_s1_charge_water_oxidants", ev,
             ["water_initial", "nitric_acid", "sodium_dichromate", "sodium_periodate"],
             "Charge water, 70% nitric acid, sodium dichromate dihydrate and sodium periodate; stir vigorously on an ice bath.",
             temperature=("ice bath", "ice bath"), time=("15 min", q(15, "min"))),
        step("acid_s2_prepare_alcohol_solution", ev, ["butynol", "water_alcohol_solution"],
             "Dilute 3-butyn-1-ol with chilled deionized water."),
        step("acid_s3_slow_addition", ev,
             ["butynol", "water_alcohol_solution", "nitric_acid", "sodium_dichromate", "sodium_periodate"],
             "Add the chilled 3-butyn-1-ol solution slowly to the oxidant mixture."),
        step("acid_s4_overnight_stir", ev,
             ["butynol", "nitric_acid", "sodium_dichromate", "sodium_periodate", "acid_product"],
             "Stir the reaction mixture overnight.", time=("overnight", "overnight")),
        step("acid_s5_filter", ev, ["acid_product"],
             "Filter the reaction mixture to remove insoluble salts.", stage="workup"),
        step("acid_s6_extract", ev, ["acid_product", "diethyl_ether"],
             "Extract the product from the aqueous phase with diethyl ether five times (80 mL each).", stage="workup"),
        step("acid_s7_dry_concentrate", ev, ["acid_product", "magnesium_sulfate"],
             "Combine the organic fractions, dry over anhydrous magnesium sulfate and concentrate by rotary evaporation to an orange viscous liquid.", stage="workup"),
        step("acid_s8_pentane_wash", ev, ["acid_product", "pentane_wash"],
             "Wash with pentane three times (200 mL each) and concentrate again by rotary evaporation to obtain an off-white solid.", stage="purification"),
    ]
    r["procedure"]["original_text"] = f(
        "Deionized water (135 mL) ... 70% HNO3 (0.170 mL, 0.05 eq., 2.5 mmol), Na2Cr2O7·2H2O (0.15 g, 0.01 eq., 0.5 mmol), and NaIO4 (34.2 g, 3.2 eq., 160 mmol) ... ice bath for 15 min. 3-Butyn-1-ol (3.50 g, 1 eq., 50 mmol) ... chilled deionized water (135 mL), added slowly, and stirred overnight ... diethyl ether (80 mL × 5) ... magnesium sulfate ... pentane (200 mL × 3) ... 3-butynoic acid (2.19 g, 52% yield).",
        "See ordered charges, reaction, workup and purification steps.", ev,
    )
    r["procedure"]["general_procedure_reference"] = f(
        "following a previously reported procedure, with modification (reference 36)",
        "modified literature procedure; exact modifications not enumerated", ev,
    )
    r["procedure"]["addition_order"] = f(
        "Water; HNO3; Na2Cr2O7·2H2O; NaIO4; then slow addition of chilled aqueous 3-butyn-1-ol.",
        ["water", "HNO3", "Na2Cr2O7·2H2O", "NaIO4", "aqueous 3-butyn-1-ol"], ev,
    )
    r["procedure"]["pretreatment"] = f(
        "500 mL single-neck round-bottom flask fitted with magnetic stirrer bar",
        "500 mL single-neck round-bottom flask with stir bar", ev,
    )
    r["outcomes"] = [
        outcome("acid_product_mass", "product_amount", [acid], ev,
                ("2.19 g", q(2.19, "g")), product_role="main product",
                measurement_basis="mass of off-white solid after pentane washing and rotary evaporation"),
        outcome("acid_reported_yield", "yield", [acid], ev,
                ("52% yield", q(52, "%")), product_role="main product",
                measurement_basis="reported alongside 2.19 g purified product; calculation method not stated"),
    ]
    set_manual_classification(r, ev, [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}], "before_order": 1, "after_order": 2, "change_type": "order_increased"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}], "before_order": 1, "after_order": 0, "change_type": "hydrogen_loss", "hydrogen_count": 2},
    ], "Net primary-alcohol oxidation is structure-derived from the named substrate and product. Oxygen sources, oxidant reduction products and complete atom mapping are unresolved.")
    r["notes"] = [
        f("The source does not state whether the ice bath was maintained during slow addition or the overnight hold.",
          "Temperature after the initial 15 min ice-bath period is unreported.", ev),
        f("The source reports '52% yield' after purification but does not explicitly label the yield as isolated or specify a measurement method.",
          "Yield type remains unreported; product mass and operational purification basis are retained separately.", ev),
    ]
    r["review"]["unresolved"] = [
        "Temperature during the slow addition and overnight reaction is not reported.",
        "Yield type and calculation basis are not explicitly reported.",
        "Complete atom mapping and oxidant-derived oxygen identity are unresolved; bond classification is partial.",
        "Reference 36 is cited as the prior procedure, but its full experimental record is not imported as a current-paper experiment.",
    ]
    return r


def reaction_iron_oleate():
    ev = PAGE9
    fecl3 = f"{DOC}:mol_iron_iii_chloride_hexahydrate"
    sodium_oleate = f"{DOC}:mol_sodium_oleate"
    ethanol = f"{DOC}:mol_ethanol"
    water = f"{DOC}:mol_water"
    hexane = f"{DOC}:mol_n_hexane"
    product = f"{DOC}:mol_iron_oleate"
    r = reaction_base(
        "rxn_iron_oleate_precursor_synthesis", ev,
        {"page": 9, "section": "Methods—Synthesis of iron-oleate precursor"},
        "Synthesis of iron-oleate precursor",
        "formation and organic-phase isolation of an iron-oleate precursor",
    )
    r["participants"] = [
        participant("iron_chloride", fecl3, "reactant", ev,
                    mass=("10.8 g", q(10.8, "g")), amount=("40 mmol", q(40, "mmol"))),
        participant("ethanol_reaction", ethanol, "solvent", ev,
                    volume=("ethanol (80 mL)", q(80, "mL"))),
        participant("water_reaction", water, "solvent", ev,
                    volume=("nanopure water (60 mL)", q(60, "mL"))),
        participant("sodium_oleate", sodium_oleate, "reactant", ev,
                    mass=("36.5 g", q(36.5, "g")), amount=("120 mmol", q(120, "mmol"))),
        participant("hexane_reaction", hexane, "solvent", ev,
                    volume=("n-hexane (140 mL)", q(140, "mL"))),
        participant("water_wash", water, "workup_agent", ev, stage="workup",
                    volume=("nanopure water (30 mL × 3)", q(90, "mL", raw_text="30 mL × 3"))),
        participant("iron_oleate_product", product, "product", ev, product_role="precursor product"),
    ]
    r["reaction_scale"] = f(
        "iron(III) chloride hexahydrate: 40 mmol",
        {"reference_molecule_id": fecl3, "amount": q(40, "mmol")}, ev,
    )
    r["solvent_conditions"]["composition"] = f(
        "80 mL ethanol, 60 mL nanopure water and 140 mL n-hexane",
        [{"molecule_id": ethanol, "volume": q(80, "mL")},
         {"molecule_id": water, "volume": q(60, "mL")},
         {"molecule_id": hexane, "volume": q(140, "mL")}], ev,
    )
    r["solvent_conditions"]["mixture_ratio"] = f(
        "ethanol:water:n-hexane = 80:60:140 by reported volume",
        {"ethanol": 80, "water": 60, "n_hexane": 140, "basis": "mL charges"}, ev,
    )
    r["solvent_conditions"]["total_volume"] = missing(ev)
    r["steps"] = [
        step("iron_oleate_s1_dissolve_iron", ev,
             ["iron_chloride", "ethanol_reaction", "water_reaction"],
             "Dissolve iron(III) chloride hexahydrate in ethanol and nanopure water in a 500 mL three-neck round-bottom flask."),
        step("iron_oleate_s2_add_oleate_hexane", ev, ["sodium_oleate", "hexane_reaction"],
             "Quickly add sodium oleate together with n-hexane and stir until the sodium oleate is completely dissolved."),
        step("iron_oleate_s3_reflux", ev,
             ["iron_chloride", "sodium_oleate", "iron_oleate_product"],
             "Heat the reaction solution at 70 °C for a 4 h reflux under a nitrogen blanket.",
             temperature=("70 °C", q(70, "°C")), time=("4 h", q(4, "h")), atmosphere="N2 blanket"),
        step("iron_oleate_s4_cool", ev, ["iron_oleate_product"],
             "Cool the reaction solution to room temperature.", temperature=("room temperature", "room temperature")),
        step("iron_oleate_s5_wash", ev, ["iron_oleate_product", "water_wash"],
             "Separate the upper organic layer containing iron oleate and wash it three times with 30 mL nanopure water in a separatory funnel.", stage="workup"),
        step("iron_oleate_s6_concentrate", ev, ["iron_oleate_product"],
             "Concentrate the iron oleate by rotary evaporation.", stage="workup"),
        step("iron_oleate_s7_vacuum", ev, ["iron_oleate_product"],
             "Transfer to a 100 mL round-bottom flask connected to a Schlenk line and place under vacuum overnight.",
             stage="workup", time=("overnight", "overnight"), pressure=("under vacuum", "under vacuum")),
    ]
    r["procedure"]["original_text"] = f(
        "iron (III) chloride hexahydrate (10.8 g, 40 mmol) ... 80 mL ethanol and 60 mL nanopure water ... Sodium oleate (36.5 g, 120 mmol) ... 140 mL n-hexane ... heated to 70 °C for a 4-h reflux under a N2 blanket ... upper organic layer ... washed three times with 30 mL nanopure water ... concentrated ... under vacuum overnight.",
        "See ordered charges, reflux and isolation steps.", ev,
    )
    r["procedure"]["general_procedure_reference"] = f(
        "previously established method (reference 34)",
        "literature method; no modifications stated", ev,
    )
    r["procedure"]["addition_order"] = f(
        "Dissolve iron(III) chloride in ethanol/water; quickly add sodium oleate with n-hexane.",
        ["iron(III) chloride in ethanol/water", "sodium oleate and n-hexane"], ev,
    )
    r["procedure"]["pretreatment"] = f(
        "500 mL three-neck round-bottom flask containing a Teflon-coated egg-shaped magnetic stir bar",
        "500 mL three-neck flask with stir bar", ev,
    )
    r["outcomes"] = [
        outcome("iron_oleate_obtained", "other", [product], ev,
                ("resulting iron oleate", "iron-oleate precursor obtained"),
                product_role="precursor product",
                measurement_basis="operational isolation of the upper organic layer followed by rotary evaporation and overnight vacuum treatment"),
    ]
    set_manual_classification(r, ev, [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "Fe"}, {"atom_map": None, "element": "O"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
    ], "Iron-carboxylate Fe-O formation is inferred from the named iron-oleate product. Product stoichiometry, coordination mode, counterion fate and atom mapping are not reported.")
    r["notes"] = [
        f("The source reports reflux at 70 °C in a mixed ethanol/water/n-hexane system.",
          {"temperature": q(70, "°C"), "mode": "reflux"}, ev),
        f("No mass, molar amount, chemical yield or compositional analysis is reported for the isolated iron-oleate precursor.",
          "Product amount and yield unreported.", ev),
    ]
    r["review"]["unresolved"] = [
        "Iron-oleate stoichiometry and coordination structure are not reported; product SMILES remains unresolved.",
        "No product mass, yield or compositional measurement basis is reported.",
        "Fe-O classification is a partial inference from the named product; complete atom mapping is unavailable.",
    ]
    return r


def reaction_magnetite():
    ev = PAGE9
    iron_oleate = f"{DOC}:mol_iron_oleate"
    oleic_acid = f"{DOC}:mol_oleic_acid"
    octadecene = f"{DOC}:mol_1_octadecene"
    tetradecene = f"{DOC}:mol_1_tetradecene"
    ethanol = f"{DOC}:mol_ethanol"
    hexane = f"{DOC}:mol_n_hexane"
    product = f"{DOC}:mol_fe3o4_nanocrystals"
    r = reaction_base(
        "rxn_8nm_magnetite_nanocrystals_synthesis", ev,
        {"page": 9, "section": "Methods—Synthesis of 8 nm magnetite (Fe3O4) nanocrystals"},
        "Synthesis of 8 nm magnetite (Fe3O4) nanocrystals",
        "thermal decomposition of iron oleate to magnetite nanocrystals",
        categories=("scale_up", "other"),
    )
    r["participants"] = [
        participant("iron_oleate", iron_oleate, "reactant", ev,
                    mass=("9.0 g", q(9.0, "g")), amount=("10 mmol", q(10, "mmol"))),
        participant("oleic_acid", oleic_acid, "reagent", ev,
                    mass=("1.56 g", q(1.56, "g")), amount=("5.5 mmol", q(5.5, "mmol"))),
        participant("octadecene", octadecene, "solvent", ev,
                    mass=("1-octadecene (40 g)", q(40, "g"))),
        participant("tetradecene", tetradecene, "solvent", ev,
                    mass=("1-tetradecene (10 g)", q(10, "g"))),
        participant("ethanol_precipitation", ethanol, "workup_agent", ev, stage="workup"),
        participant("hexane_dispersion", hexane, "workup_agent", ev, stage="workup"),
        participant("magnetite_product", product, "product", ev, product_role="nanocrystal product"),
    ]
    r["reaction_scale"] = f(
        "iron oleate: 10 mmol (9.0 g)",
        {"reference_molecule_id": iron_oleate, "amount": q(10, "mmol"), "mass": q(9.0, "g")}, ev,
    )
    r["solvent_conditions"]["composition"] = f(
        "1-octadecene:1-tetradecene (40:10 g)",
        [{"molecule_id": octadecene, "mass": q(40, "g")},
         {"molecule_id": tetradecene, "mass": q(10, "g")}], ev,
    )
    r["solvent_conditions"]["mixture_ratio"] = f(
        "1-octadecene:1-tetradecene = 40:10 by mass",
        {"1_octadecene": 40, "1_tetradecene": 10, "basis": "g charges"}, ev,
    )
    r["solvent_conditions"]["total_volume"] = missing(ev)
    ramp = step("magnetite_s2_heat_ramp", ev,
                ["iron_oleate", "oleic_acid", "octadecene", "tetradecene"],
                "Place under nitrogen flow and heat from 110 °C to 300 °C at 3.3 °C/min using a temperature controller.",
                temperature=("110 °C to 300 °C", {"minimum": 110, "maximum": 300, "unit": "°C", "raw_text": "110 °C to 300 °C"}),
                atmosphere="nitrogen flow")
    set_inferred_ramp_time(ramp)
    r["steps"] = [
        step("magnetite_s1_degas", ev,
             ["iron_oleate", "oleic_acid", "octadecene", "tetradecene"],
             "Stir and degas the mixture on a Schlenk line at 110 °C.",
             temperature=("110 °C", q(110, "°C")), time=("3 h", q(3, "h")), pressure=("degassed on a Schlenk line", "reduced pressure / Schlenk degassing")),
        ramp,
        step("magnetite_s3_hold_300", ev,
             ["iron_oleate", "oleic_acid", "magnetite_product"],
             "Count reaction time from reaching 300 °C and hold for 1 h under nitrogen flow.",
             temperature=("300 °C", q(300, "°C")), time=("1 h", q(1, "h")), atmosphere="nitrogen flow"),
        step("magnetite_s4_cool", ev, ["magnetite_product"],
             "Quickly cool the reaction solution to room temperature by blowing air across the flask.",
             temperature=("room temperature", "room temperature")),
        step("magnetite_s5_precipitate_isolate", ev, ["magnetite_product", "ethanol_precipitation"],
             "Precipitate the iron-oxide nanocrystals with ethanol and isolate by centrifugation at 5500 × g.",
             stage="workup", time=("20 min", q(20, "min"))),
        step("magnetite_s6_wash_cycles", ev, ["magnetite_product", "ethanol_precipitation"],
             "Wash with three ethanol precipitation/redispersion cycles, each followed by centrifugation at 5500 × g.",
             stage="workup", time=("20 min per centrifugation", q(20, "min"))),
        step("magnetite_s7_dry", ev, ["magnetite_product"],
             "Dry the nanocrystals under vacuum and weigh them.", stage="workup", pressure=("under vacuum", "under vacuum")),
        step("magnetite_s8_redisperse", ev, ["magnetite_product", "hexane_dispersion"],
             "Redisperse the nanocrystals in n-hexane at 20 mg/mL for further use.", stage="workup"),
    ]
    r["special_conditions"]["reactor"] = f(
        "250 mL three-neck round-bottom flask; Schlenk line; temperature controller",
        ["250 mL three-neck flask", "Schlenk line", "temperature controller"], ev,
    )
    r["procedure"]["original_text"] = f(
        "Iron oleate (10 mmol, 9.0 g) and oleic acid (5.5 mmol, 1.56 g) ... 1-octadecene:1-tetradecene (40:10 g) ... degassed ... 110 °C for 3 h ... under nitrogen flow ... heated to 300 °C at 3.3 °C/min ... 1 h ... quickly cooled ... precipitated using ethanol ... centrifugation (5500 × g, 20 min) ... three precipitation/redispersion cycles ... dried under vacuum ... redispersed in n-hexane at 20 mg/mL.",
        "See ordered degassing, ramp, hold, cooling and isolation steps.", ev,
    )
    r["procedure"]["addition_order"] = f(
        "Charge iron oleate and oleic acid into the 1-octadecene/1-tetradecene solvent mixture before degassing.",
        ["iron oleate", "oleic acid", "1-octadecene/1-tetradecene"], ev,
    )
    r["procedure"]["pretreatment"] = f(
        "mixture stirred and degassed on a Schlenk line at 110 °C for 3 h",
        "Schlenk degassing at 110 °C for 3 h", ev,
    )
    r["outcomes"] = [
        outcome("magnetite_nanocrystals_obtained", "other", [product], ev,
                ("resulting iron-oxide nanocrystals", "8 nm magnetite (Fe3O4) nanocrystals obtained"),
                product_role="nanocrystal product",
                measurement_basis="method heading identifies 8 nm magnetite; final material was isolated, dried, weighed and redispersed, but product mass/yield is not reported"),
    ]
    set_manual_classification(r, ev, [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "Fe"}, {"atom_map": None, "element": "O"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
    ], "Fe-O lattice bond formation is inferred from conversion of the named iron-oleate precursor to Fe3O4 nanocrystals. The heterogeneous solid and decomposition coproducts prevent complete molecular atom mapping.")
    r["notes"] = [
        f("The source reports a heating rate of 3.3 °C/min from the 110 °C degassing temperature to 300 °C; the corresponding ~57.6 min ramp duration is calculated and marked inferred.",
          {"temperature_start": q(110, "°C"), "temperature_end": q(300, "°C"), "heating_rate": "3.3 °C/min", "calculated_transition_time": q(57.6, "min", qualifier="approx")},
          ev, method="inference", derivation="(300-110)/3.3 = 57.6 min"),
        f("The final product was redispersed at 20 mg/mL, but the source does not report the isolated product mass or chemical yield.",
          {"final_dispersion_concentration": "20 mg/mL", "isolated_mass": None, "yield": None}, ev),
    ]
    r["review"]["unresolved"] = [
        "No isolated nanocrystal mass or chemical yield is reported.",
        "Iron-oleate precursor composition and decomposition coproducts are unresolved.",
        "The Fe3O4 nanocrystal is heterogeneous and has no unique finite SMILES; Fe-O classification is partial.",
        "The ~57.6 min ramp duration is calculated from the reported rate and endpoints, not directly reported by the authors.",
    ]
    return r


def add_context(bundle):
    contexts = [
        {
            "context_id": f"{DOC}:context_page7_serum_results",
            "kind": "bioanalytical_results",
            "evidence_ids": [P[7], I[7]],
            "description": "Page 7 reports cTnI/HSA enrichment performance, top-down-MS limits of detection (0.06 ng/mL for direct cTnI detection and 0.75 ng/mL for the integrated serum workflow), an NP-Pep enrichment factor of 115-fold, and serum cTnI recoveries of about 51% versus about 17% for Agarose-mAb. These are protein-assay results, not small-molecule synthesis yields or conversions.",
        },
        {
            "context_id": f"{DOC}:context_page8_fig4_proteoforms",
            "kind": "bioanalytical_workflow_and_mass_spectra",
            "evidence_ids": [P[8], I[8]],
            "description": "Fig. 4 shows serum spike-in, NP-Pep binding, magnetic isolation, elution, top-down LC/MS, and six deconvoluted cTnI proteoform spectra. Protein phosphorylation/truncation masses are analytical observations and do not define additional covalent synthesis experiments.",
        },
        {
            "context_id": f"{DOC}:context_page9_characterization_methods",
            "kind": "characterization_methods",
            "evidence_ids": [P[9], I[9]],
            "description": "Page 9 includes NMR and ESI-MS values for 3-butynoic acid and BAPTES plus TEM/FTIR/TGA/zeta-potential and FTICR-MS methods. These support the associated preparation records but are not separate synthesis experiments or peak-table reactions.",
        },
        {
            "context_id": f"{DOC}:context_page10_protein_workflows",
            "kind": "biochemical_sample_preparation_and_enrichment",
            "evidence_ids": [P[10], I[10]],
            "description": "Page 10 gives cardiac-tissue extraction, NP-Pep and agarose cTnI enrichment, serum spike-in, washing, elution, desalting and SDS-PAGE/LC-MS sample preparation. These are noncovalent biochemical workflows; no new synthesized molecular product is claimed.",
        },
        {
            "context_id": f"{DOC}:context_page11_lcms_statistics",
            "kind": "analytical_and_statistical_methods",
            "evidence_ids": [P[11], I[11]],
            "description": "Page 11 reports reverse-phase chromatography, top-down MS/MS, database-search and statistical-analysis parameters, followed by data availability and references. Instrument voltages, collision energies, masses and abundance calculations are analytical conditions, not synthetic experiments.",
        },
        {
            "context_id": f"{DOC}:context_page12_references_declarations",
            "kind": "references_and_article_metadata",
            "evidence_ids": [P[12], I[12]],
            "description": "Page 12 contains references 34–62, acknowledgements, author contributions, competing interests, additional information and the license. Reference 36 is prior allene/thiol chemistry and is not recast as an experiment performed in this paper.",
        },
    ]
    existing = {c.get("context_id") for c in bundle["source_context"]}
    bundle["source_context"].extend(c for c in contexts if c["context_id"] not in existing)


def coverage():
    new_ids = [
        f"{DOC}:rxn_3_butynoic_acid_synthesis",
        f"{DOC}:rxn_iron_oleate_precursor_synthesis",
        f"{DOC}:rxn_8nm_magnetite_nanocrystals_synthesis",
    ]
    existing_page9 = [
        f"{DOC}:rxn_baptes_synthesis",
        f"{DOC}:rxn_np_baptes_silanization",
        f"{DOC}:rxn_soluble_baptes_peptide_test",
        f"{DOC}:rxn_np_pep_coupling",
    ]
    return {
        "document_id": DOC,
        "semantic_accuracy_verified": False,
        "pages": [
            {
                "page": 7,
                "state": "attempted",
                "regions": [
                    {"region": "left column continuation of Fig. 3 results", "finding": "cTnI/HSA enrichment, LOD, enrichment-factor and recovery results; bioanalytical observations, not covalent syntheses.", "evidence_ids": [P[7], I[7]]},
                    {"region": "right column Proteoform-resolved analysis and Discussion", "finding": "Serum spike-in proteoform results and interpretation; no reaction arrow, small-molecule preparation, or optimization table row.", "evidence_ids": [P[7], I[7]]},
                ],
                "record_ids": [],
                "unresolved": ["Cited Supplementary Figs. 21–31 and Tables 7–8 are not supplied as explicitly associated SI in this workspace."],
            },
            {
                "page": 8,
                "state": "attempted",
                "regions": [
                    {"region": "Fig. 4a workflow", "finding": "Noncovalent serum spike-in/NP-Pep enrichment and top-down LC/MS workflow; no new covalent product.", "evidence_ids": [P[8], I[8]]},
                    {"region": "Fig. 4b panels i–vi and caption", "finding": "Six cTnI proteoform mass spectra from donor, diseased and postmortem samples; independent protein masses are characterization, not synthesis experiments.", "evidence_ids": [P[8], I[8]]},
                    {"region": "bottom Discussion continuation", "finding": "Prospective clinical and broader-protein applications; no additional experiment or reaction scheme.", "evidence_ids": [P[8], I[8]]},
                ],
                "record_ids": [],
                "unresolved": ["Supplementary Table 8 and source-data files cited by Fig. 4 are not available as associated source files."],
            },
            {
                "page": 9,
                "state": "attempted",
                "regions": [
                    {"region": "left column Methods—Synthesis of 3-butynoic acid", "finding": "One specific oxidation/preparation experiment with charges, overnight reaction, extraction, purification, 2.19 g product and 52% reported yield.", "evidence_ids": [P[9], I[9]]},
                    {"region": "left column Methods—Synthesis of N-(3-(triethoxysilyl)propyl)buta-2,3-dienamide", "finding": "Detailed method for the existing BAPTES synthesis record; merged as duplicate evidence, not a new experiment. NMR/MS values are characterization.", "evidence_ids": [P[9], I[9]]},
                    {"region": "right column Methods—Synthesis of iron-oleate precursor", "finding": "One specific precursor preparation with FeCl3·6H2O, sodium oleate, mixed solvents, 70 °C/4 h N2 reflux, aqueous washes and overnight vacuum treatment.", "evidence_ids": [P[9], I[9]]},
                    {"region": "right column Methods—Synthesis of 8 nm magnetite nanocrystals", "finding": "One specific thermal nanocrystal synthesis with 110 °C/3 h degassing, 3.3 °C/min ramp to 300 °C, 1 h hold, ethanol isolation/washes and hexane redispersion.", "evidence_ids": [P[9], I[9]]},
                    {"region": "right column Methods—Synthesis of NP-BAPTES and NP-Pep", "finding": "Detailed methods for existing silanization and peptide-coupling records; merged as duplicate evidence, not additional experiments.", "evidence_ids": [P[9], I[9]]},
                    {"region": "right column Material characterization and Small molecule analysis by FTICR-MS", "finding": "TEM/FTIR/TGA/zeta-potential and MS acquisition methods; no standalone synthesis or peak-table reaction.", "evidence_ids": [P[9], I[9]]},
                ],
                "record_ids": new_ids + existing_page9,
                "unresolved": [
                    "References 34 and 36 identify prior procedures but are not imported as current-paper experiments.",
                    "Iron-oleate coordination/stoichiometry and Fe3O4 finite molecular structure are not uniquely specified.",
                    "Cited supplementary characterization figures are not available as explicitly associated SI.",
                ],
            },
            {
                "page": 10,
                "state": "attempted",
                "regions": [
                    {"region": "top-left FTICR-MS continuation", "finding": "Mass-spectrometer settings and MS/MS analysis; analytical method only.", "evidence_ids": [P[10], I[10]]},
                    {"region": "left column cardiac tissue extraction and sarcomere enrichment", "finding": "Biochemical extraction, noncovalent affinity capture, washing, acidic elution, desalting and SDS-PAGE sample preparation; no covalent synthesis.", "evidence_ids": [P[10], I[10]]},
                    {"region": "right column serum enrichment and comparison methods", "finding": "NP-Pep/agarose serum enrichment workflows and controls; no new reaction product and no no-reaction synthesis record invented.", "evidence_ids": [P[10], I[10]]},
                    {"region": "bottom-right start of RPC procedure", "finding": "Chromatographic analytical method begins; no reaction arrow or preparative product.", "evidence_ids": [P[10], I[10]]},
                ],
                "record_ids": [],
                "unresolved": ["Agarose ligand-conjugation follows manufacturer recommendations, but the conjugation procedure and product characterization are not reported; no synthetic record was inferred."],
            },
            {
                "page": 11,
                "state": "attempted",
                "regions": [
                    {"region": "left column RPC and Top-down MS analysis", "finding": "Chromatography gradient, ion-source and mass-spectrometer parameters, peak processing and database search; characterization only.", "evidence_ids": [P[11], I[11]]},
                    {"region": "left column Statistical analysis/Data availability", "finding": "Statistical and repository metadata; no reaction experiment.", "evidence_ids": [P[11], I[11]]},
                    {"region": "references 1–33", "finding": "Cited literature only; not experiments performed in this article.", "evidence_ids": [P[11], I[11]]},
                ],
                "record_ids": [],
                "unresolved": ["Raw proteomics files are referenced through PRIDE PXD019712 but were not used as reaction-source records in this page batch."],
            },
            {
                "page": 12,
                "state": "attempted",
                "regions": [
                    {"region": "references 34–62", "finding": "Prior literature list; reference 36 concerns allenamide/cysteine chemistry but is not a current-paper experiment.", "evidence_ids": [P[12], I[12]]},
                    {"region": "acknowledgements, contributions, competing interests and additional information", "finding": "Article metadata and declarations; no experimental procedure or reaction result.", "evidence_ids": [P[12], I[12]]},
                    {"region": "license block", "finding": "Publication license; no chemistry experiment.", "evidence_ids": [P[12], I[12]]},
                ],
                "record_ids": [],
                "unresolved": ["Supplementary information is linked by DOI but was not explicitly provided or associated in the workspace and was not inferred as SI."],
            },
        ],
    }


def main():
    submitted_path = ROOT / "runs/current/submitted.json"
    bundle = json.loads(submitted_path.read_text())
    dataset = json.loads((ROOT / "runs/current/dataset.json").read_text())
    bundle["documents"] = copy.deepcopy(dataset["documents"])
    extend_page9_image_evidence(bundle)
    add_molecules(bundle)
    existing_reactions = {r["reaction_id"] for r in bundle["reactions"]}
    new_reactions = [reaction_3_butynoic_acid(), reaction_iron_oleate(), reaction_magnetite()]
    bundle["reactions"].extend(r for r in new_reactions if r["reaction_id"] not in existing_reactions)
    add_context(bundle)

    raw = {
        "document_id": DOC,
        "scope_pages": [7, 8, 9, 10, 11, 12],
        "semantic_accuracy_verified": False,
        "reviewed_page_images": [I[n] for n in range(7, 13)],
        "new_molecule_ids": [
            f"{DOC}:mol_3_butyn_1_ol",
            f"{DOC}:mol_nitric_acid",
            f"{DOC}:mol_sodium_dichromate_dihydrate",
            f"{DOC}:mol_sodium_periodate",
            f"{DOC}:mol_diethyl_ether",
            f"{DOC}:mol_magnesium_sulfate",
            f"{DOC}:mol_pentane",
            f"{DOC}:mol_iron_iii_chloride_hexahydrate",
            f"{DOC}:mol_sodium_oleate",
            f"{DOC}:mol_ethanol",
            f"{DOC}:mol_iron_oleate",
            f"{DOC}:mol_oleic_acid",
            f"{DOC}:mol_1_octadecene",
            f"{DOC}:mol_1_tetradecene",
            f"{DOC}:mol_fe3o4_nanocrystals",
        ],
        "new_reaction_ids": [r["reaction_id"] for r in new_reactions],
        "deduplicated_page9_record_ids": [
            f"{DOC}:rxn_baptes_synthesis",
            f"{DOC}:rxn_np_baptes_silanization",
            f"{DOC}:rxn_soluble_baptes_peptide_test",
            f"{DOC}:rxn_np_pep_coupling",
        ],
        "candidate_bundle": copy.deepcopy(bundle),
        "decision": "Three new preparative records from page 9; page-9 methods for existing records deduplicated; pages 7-8 and 10-12 retained as non-synthetic source context.",
    }
    st.write(ROOT / "runs/current/raw_model/pages_0007_0012_extraction.json", raw)
    st.write(submitted_path, bundle)
    st.write(ROOT / "runs/current/coverage/pages_0007_0012.json", coverage())


if __name__ == "__main__":
    main()
