#!/usr/bin/env python3
"""Build the v3.1 source-derived bundle for JACS 2021, 143, 3323-3329.

This file is a reproducible transcription aid.  It does not certify chemistry and
does not associate the unavailable Supporting Information.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DATASET = json.loads((RUN / "dataset.json").read_text())
DOC = "doc_83f760dc954f4ef1"


def ev(page):
    return [f"{DOC}:page{page}:text", f"{DOC}:page{page}:image"]


def f(raw=None, normalized=None, evidence=(), method="original_extraction", derivation=None,
      status=None):
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation, status=status)


def inferred(value, evidence, derivation, raw=None):
    return f(raw, value, evidence, method="inference", derivation=derivation)


bundle = st.empty(DATASET["documents"])
molecules = {}


def add_molecule(key, label, pages, name=None, smiles=None, chirality="unresolved",
                 geometry="unresolved", rep_status="unresolved", assignment=None,
                 mixture=None, axial=None, identity_note=None):
    evidence = []
    for page in pages:
        evidence += ev(page)
    mid = f"{DOC}:mol_{key}"
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    m["names"]["original"] = f(name or label, name or label, evidence)
    m["names"]["english"] = f(name or label, name or label, evidence)
    m["names"]["abbreviations"] = f(label, [label], evidence)
    m["cas_numbers"] = f()
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = f(
            "structure drawing at cited page",
            smiles,
            evidence,
            method="structure_derivation",
            derivation="Candidate connectivity was transcribed from the source drawing; RDKit parsing does not verify the drawing or stereochemical assignment.",
        )
        if rep_status == "encoded":
            diag = st.stereo_inspect(smiles)
            m["stereochemistry"]["isomeric_smiles"] = f(
                None, diag["canonical_isomeric_smiles"], evidence,
                method="structure_derivation",
                derivation="Canonical isomeric SMILES calculated from the transcribed candidate graph."
            )
            m["stereochemistry"]["identity_key"] = f(
                None, diag["stereo_graph_key"], evidence,
                method="structure_derivation",
                derivation="Stereo graph key calculated from the candidate SMILES; it is not a sample identity."
            )
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = f()
    m["stereochemistry"]["chirality_status"] = (
        f(chirality, chirality, evidence) if chirality in {"specified", "racemic", "mixture"}
        else inferred(chirality, evidence, "Conservative stereochemical status after inspection of the cited drawing; no absent tag was interpreted as achiral or racemic.")
    )
    m["stereochemistry"]["geometry_status"] = inferred(
        geometry, evidence, "Geometry status assessed from the cited drawing and labels; unassigned geometry remains unresolved/unspecified."
    )
    m["stereochemistry"]["representation_status"] = inferred(
        rep_status, evidence, "Comparison of the source stereochemical information with the stored molecular representation."
    )
    if assignment is not None:
        m["stereochemistry"]["assignment"] = f(assignment, assignment, evidence)
    if mixture is not None:
        m["stereochemistry"]["mixture"] = f(mixture, mixture, evidence)
    if axial is not None:
        m["stereochemistry"]["axial"] = f(axial, axial, evidence)
    m["identity_note"] = inferred(
        identity_note or "Source label/drawing retained; exact identity and stereochemical encoding remain machine-extracted and need expert review.",
        evidence,
        "Identity note records the limits of the visual transcription."
    )
    molecules[key] = m
    bundle["molecules"].append(m)
    return mid


def mid(key):
    return molecules[key]["molecule_id"]


# Core structures.  Only a small number of confidently legible connectivities are
# encoded; dense scope structures retain source labels rather than unsafe SMILES.
add_molecule("1a", "(±)-1a", [1, 2, 3, 4], "racemic tertiary α-allenyl alcohol 1a",
             "CC(O)(C=C=C)c1ccc2ccccc2c1", "racemic", "not_applicable", "partial",
             mixture="racemic (±)", assignment="racemic tertiary alcohol; no enantiomer selected")
add_molecule("2a", "2a", [1, 2, 3], "TBS ketene acetal 2a",
             "CC(C)=C(OCC)O[Si](C)(C)C(C)(C)C", "achiral", "not_applicable", "encoded")
add_molecule("2b", "2b", [2, 3], "TMS ketene acetal 2b",
             "CC(C)=C(OCC)O[Si](C)(C)C", "achiral", "not_applicable", "encoded")
add_molecule("2c", "2c", [2], "TES ketene acetal 2c",
             "CC(C)=C(OCC)O[Si](CC)(CC)CC", "achiral", "not_applicable", "encoded")
add_molecule("3a", "(R)-3a", [1, 2, 3, 4], "β-allenyl ester (R)-3a",
             "CCOC(=O)C(C)(C)C(C)(C=C=C)c1ccc2ccccc2c1", "specified", "not_applicable", "partial",
             assignment="(R)-3a; absolute configuration stated in the article and drawn with a stereobond; candidate SMILES stores connectivity only")

# Scope molecules.
for letter in "bcdefghijklmnopqrst":
    add_molecule(f"sub_{letter}", f"Table 2 electrophile leading to 3{letter}", [3],
                 f"racemic tertiary α-allenyl alcohol substrate for 3{letter}",
                 chirality="racemic", geometry="not_applicable", rep_status="not_encoded",
                 mixture="racemic (±) substrate series")
for letter in "bcdefghijklmnopqrstuvwxzy":
    # All product panels show a stereobond.  The article does not print an R/S
    # descriptor beside every panel, so the assignment remains textual.
    assignment = "stereobond shown in scope panel; absolute descriptor is not printed beside this individual product"
    if letter == "t":
        assignment = "article states 3t has absolute configuration opposite to parent product 3a"
    if letter == "v":
        assignment = "major product shown; footnote 79 states the same diastereomer is major with Z-2e, but individual stereocenter descriptors are not printed"
    product_pages = [3] if letter <= "t" else ([3, 4, 7] if letter == "v" else [4])
    add_molecule(f"3{letter}", f"3{letter}", product_pages,
                 f"scope product 3{letter}", chirality="partially_specified" if letter == "v" else "specified",
                 geometry="not_applicable", rep_status="not_encoded", assignment=assignment,
                 mixture="diastereomeric product; 2.4:1 dr with E-enriched 2e and 1.7:1 dr with Z-2e" if letter == "v" else None)

# Table 3 nucleophiles, including the explicitly discussed negative tert-butyl derivative.
for key, label, name, geom, mix in [
    ("2d", "2d", "methyl ester-derived TBS ketene acetal 2d", "not_applicable", None),
    ("2e", "(E)-2e", "propionate-derived silyl ketene acetal 2e", "mixture", "E:Z = 6.7:1, E-enriched"),
    ("2f", "2f", "Danishefsky's diene 2f", "unspecified", None),
    ("2g", "2g", "cyclohexane carboxylate-derived silyl ketene acetal 2g", "unspecified", None),
    ("2h", "2h", "4,4-difluorocyclohexane carboxylate-derived silyl ketene acetal 2h", "unspecified", None),
    ("2i", "2i", "tetrahydropyran carboxylate-derived silyl ketene acetal 2i", "unspecified", None),
    ("2tbu", "tert-butyl derivative", "tert-butyl ester-derived silyl ketene acetal", "unresolved", None),
]:
    add_molecule(key, label, [3, 4] if key == "2tbu" else [4], name,
                 chirality="unresolved", geometry=geom, rep_status="not_encoded", mixture=mix)

# Footnote 79 reports the Z isomer as a distinct nucleophile experiment.  It is
# registered separately from the E-enriched 2e mixture; the complete drawing is
# not re-encoded as an unsafe SMILES.
add_molecule(
    "2e_z", "(Z)-2e", [4, 7], "(Z)-isomer of propionate-derived silyl ketene acetal 2e",
    chirality="unresolved", geometry="specified", rep_status="not_encoded",
    assignment="(Z)-isomer explicitly stated in footnote 79",
    identity_note="Separate source identity for the Z geometric isomer used in footnote 79; complete connectivity and isomeric SMILES remain unencoded and need expert review.",
)
molecules["2e_z"]["stereochemistry"]["double_bonds"] = f(
    "(Z)-isomer of silyl ketene acetal",
    [{"assignment": "Z", "location": "silyl ketene acetal C=C; atom indices unavailable because no complete SMILES is asserted"}],
    ev(7),
)
molecules["2e_z"]["stereochemistry"]["geometry_status"] = f(
    "(Z)-isomer", "specified", ev(7)
)

# Scheme 2 products/intermediate.
add_molecule("5", "(R,Z)-5", [3, 4], "allylsilane (R,Z)-5", chirality="specified", geometry="specified",
             rep_status="not_encoded", assignment="(R,Z)-5; source reports 99% ee and >20:1 Z:E")
add_molecule("6", "(R,Z)-6", [3, 4], "amine (R,Z)-6", chirality="specified", geometry="specified",
             rep_status="not_encoded", assignment="(R,Z)-6; source reports 99% ee and >20:1 Z:E")
add_molecule("7", "(R)-7", [3, 4], "carboxylic acid (R)-7", chirality="specified", geometry="not_applicable",
             rep_status="not_encoded", assignment="(R)-7; 98% ee shown")
add_molecule("8", "(R)-8", [3, 4], "carbamate (R)-8", chirality="specified", geometry="not_applicable",
             rep_status="not_encoded", assignment="(R)-8; 98% ee shown")
add_molecule("9", "(R,S)-9", [4], "lactone (R,S)-9", chirality="specified", geometry="unspecified",
             rep_status="not_encoded", assignment="(R,S)-9 assigned by X-ray; 5.2:1 dr and 98% ee shown")

# Generic/current and cited prior-work structures.
for key, label, name, pages in [
    ("gen1", "racemic 1", "generic racemic tertiary α-allenyl alcohol", [2]),
    ("gen2", "2", "generic tetrasubstituted silyl ketene acetal", [2, 4]),
    ("gen3", "3", "generic β-allenyl ester product", [2, 4]),
    ("stoltz_nuc", "Stoltz malonodinitrile", "substituted malonodinitrile (prior work)", [2]),
    ("stoltz_el", "Stoltz allylic carbonate", "3,3-disubstituted allylic carbonate (prior work)", [2]),
    ("stoltz_prod", "Stoltz product", "acyclic vicinal-quaternary product (prior work)", [2]),
    ("jorg_aldehyde", "Jørgensen aldehyde", "α-branched aldehyde (prior work)", [2]),
    ("jorg_prod", "Jørgensen product", "vicinal-quaternary 1,4-dialdehyde (prior work)", [2]),
    ("ir_bound", "Ir-bound tertiary carbocation", "proposed Ir-stabilized allenylic carbocation intermediate", [1, 2]),
]:
    add_molecule(key, label, pages, name, chirality="unresolved", geometry="unresolved", rep_status="unresolved")

# Reagents, catalysts, ligands, and solvents.  Their source identities are retained;
# CAS values and unshown complete structures are deliberately not invented.
agents = {
    "ir": ("[Ir(cod)Cl]2", "bis(1,5-cyclooctadiene)diiridium(I) dichloride", [1, 2, 3, 4]),
    "rh": ("[Rh(cod)Cl]2", "bis(1,5-cyclooctadiene)dirhodium(I) dichloride", [2]),
    "L1": ("(S)-L1", "chiral phosphoramidite-olefin ligand (S)-L1", [1, 2, 3, 4]),
    "L2": ("(S)-L2", "reduced ligand analog (S)-L2", [2]),
    "La": ("La(OTf)3", "lanthanum triflate", [1, 2, 3, 4]),
    "Zn": ("Zn(OTf)2", "zinc triflate", [2]),
    "dioxane": ("1,4-dioxane", "1,4-dioxane", [1, 2, 3, 4]),
    "base4": ("4", "2,6-di-tert-butyl-4-methylpyridine", [2]),
    "moco6": ("Mo(CO)6", "molybdenum hexacarbonyl", [4]),
    "ph2sih2": ("Ph2SiH2", "diphenylsilane", [4]),
    "phme": ("PhMe", "toluene", [4]),
    "pdoac2": ("Pd(OAc)2", "palladium(II) acetate", [4]),
    "dppf": ("dppf", "1,1'-bis(diphenylphosphino)ferrocene", [4]),
    "pyrrolidine": ("pyrrolidine", "pyrrolidine", [4]),
    "phi": ("PhI", "iodobenzene", [4]),
    "mecn": ("MeCN", "acetonitrile", [4]),
    "dibal": ("DIBAL-H", "diisobutylaluminum hydride", [4]),
    "dcm": ("DCM", "dichloromethane", [4]),
    "dmp": ("DMP", "Dess-Martin periodinane", [4]),
    "naclo2": ("NaClO2", "sodium chlorite", [4]),
    "kh2po4": ("KH2PO4", "potassium dihydrogen phosphate", [4]),
    "methylbutene": ("2-methyl-2-butene", "2-methyl-2-butene", [4]),
    "tbuoh": ("tBuOH", "tert-butanol", [4]),
    "water": ("H2O", "water", [4]),
    "dppa": ("DPPA", "diphenyl phosphoryl azide", [4]),
    "et3n": ("Et3N", "triethylamine", [4]),
    "bnoh": ("BnOH", "benzyl alcohol", [4]),
    "dpephos": ("DPEphos", "bis[(2-diphenylphosphino)phenyl] ether", [4]),
    "dce": ("1,2-DCE", "1,2-dichloroethane", [4]),
    "ag2co3": ("Ag2CO3", "silver carbonate", [2]),
    "pnbzoh": ("4-NO2PhCO2H", "4-nitrobenzoic acid", [2]),
    "dabco": ("DABCO", "1,4-diazabicyclo[2.2.2]octane", [2]),
    "tbd": ("TBD", "1,5,7-triazabicyclo[4.4.0]dec-5-ene", [2]),
}
for key, (label, name, pages) in agents.items():
    add_molecule(key, label, pages, name, chirality="unresolved", geometry="unresolved", rep_status="unresolved")


def set_charge(p, kind, raw, value, evidence, ref=None, derived=False):
    if derived:
        p["charges"][kind] = f(raw, value, evidence, method="inference",
                                derivation="Loading/equivalent denominator is inferred as the limiting substrate shown in the same reaction scheme.")
    else:
        p["charges"][kind] = f(raw, value, evidence)
    if kind in ("equivalents", "mol_percent"):
        p["equivalence_reference"] = inferred(
            ref, evidence, "Catalyst/additive loading or reagent equivalents are referenced to the limiting substrate depicted for the experiment.", raw="relative to substrate"
        )


def participant(r, pid, key, role, evidence, stage="reaction", charges=None):
    p = st.new_participant(pid, mid(key), role, evidence, stage=stage)
    for kind, raw, value, ref in charges or []:
        set_charge(p, kind, raw, value, evidence, ref)
    r["participants"].append(p)
    return p


def set_meta(r, record_type, categories, location, evidence, reaction_name, transformation, unresolved):
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, evidence)
    r["location"] = f(location, location, evidence)
    r["reaction_name"] = f(reaction_name, reaction_name, evidence)
    r["transformation_type"] = inferred(transformation, evidence, "Transformation class assigned from the source arrow and prose.")
    r["review"]["unresolved"] = list(unresolved)


def set_solvent(r, key, evidence, concentration=None, basis=None, mixture=None):
    r["solvent_conditions"]["composition"] = f(molecules[key]["names"]["original"]["raw_value"], [mid(key)], evidence)
    if concentration is not None:
        r["solvent_conditions"]["concentration"] = f(concentration[0], concentration[1], evidence)
    if basis is not None:
        r["solvent_conditions"]["concentration_basis"] = f(basis, basis, evidence)
    if mixture is not None:
        r["solvent_conditions"]["mixture_ratio"] = f(mixture, mixture, evidence)


def step(r, sid, evidence, participant_ids, operation, temperature=None, time=None, stage="reaction"):
    s = st.new_step(sid, evidence, stage=stage)
    s["participant_ids"] = participant_ids
    s["operation"] = f(operation, operation, evidence)
    if temperature is not None:
        s["temperature"] = f(temperature[0], temperature[1], evidence)
    if time is not None:
        s["time"] = f(time[0], time[1], evidence)
    r["steps"].append(s)
    return s


def outcome(r, oid, kind, targets, evidence, raw, value, basis_raw=None, basis=None,
            yield_type=None, scope=None, product_role="product"):
    o = st.new_outcome(oid, kind, [mid(t) for t in targets], evidence)
    o["value"] = f(raw, value, evidence)
    o["product_role"] = f(product_role, product_role, evidence)
    if basis_raw is not None or basis is not None:
        o["measurement_basis"] = f(basis_raw, basis, evidence)
    if yield_type is not None:
        o["yield_type"] = f(yield_type, yield_type, evidence)
    if scope is not None:
        o["scope"] = f(scope, scope, evidence)
    r["outcomes"].append(o)
    return o


def classification(r, evidence, labels, changes, coverage="partial", limitations=None):
    r["classification"]["bond_changes"] = [
        f(None, change, evidence, method="structure_derivation",
          derivation="Partial net bond change inferred from the source reactant/product drawings; atom mapping is not claimed complete.")
        for change in changes
    ]
    r["classification"]["labels"] = f(
        None, labels, evidence, method="structure_derivation",
        derivation="Taxonomy labels include all and only the formed/promoted bonds in the stored partial bond-change list."
    )
    r["classification"]["coverage"] = inferred(coverage, evidence, "Coverage reflects whether a complete atom map was available; all source-drawing classifications here are partial or not applicable.")
    r["classification"]["limitations"] = limitations or ["No complete atom mapping was asserted.", "Implicit atom sources and coproducts may be omitted."]


CC_FORM = [{"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1"]}]
CC_CHANGES = [
    {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
    {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}], "before_order": 1, "after_order": 0, "change_type": "broken"},
]


def no_change_classification(r, evidence):
    classification(r, evidence, [], [], coverage="not_applicable", limitations=["No desired transformation was observed; no expected bond label was assigned."])


# Scheme 1: two cited examples, one current generic transformation, and the
# explicitly proposed η2-bound carbocation model.  Panel A is context only.
e2 = ev(2)
for rid, reactants, product, name in [
    ("scheme1B_stoltz_prior", ["stoltz_nuc", "stoltz_el"], "stoltz_prod", "Stoltz (2018) prior-work enantioselective allylic alkylation"),
    ("scheme1B_jorgensen_prior", ["jorg_aldehyde"], "jorg_prod", "Jørgensen (2018) oxidative aldehyde homocoupling"),
]:
    r = st.new_reaction(f"{DOC}:{rid}", DOC, e2)
    set_meta(r, "literature_example", [], "PDF page 2, Scheme 1B", e2, name, "prior-work C-C bond formation",
             ["Cited literature example only; no current-paper experimental procedure or yield is supplied in Scheme 1B."])
    for i, key in enumerate(reactants): participant(r, f"p{i+1}", key, "reactant", e2)
    participant(r, "product", product, "product", e2)
    classification(r, e2, CC_FORM, [CC_CHANGES[0]])
    bundle["reactions"].append(r)

r = st.new_reaction(f"{DOC}:scheme1C_generic", DOC, e2)
set_meta(r, "generic_scheme", [], "PDF page 2, Scheme 1C", e2,
         "Enantioselective allenylic alkylation with tetrasubstituted silyl ketene acetals",
         "allenylic C-C alkylation", ["Generic R-group transformation; individual conditions and outcomes are recorded separately."])
participant(r, "substrate", "gen1", "reactant", e2)
participant(r, "nucleophile", "gen2", "reactant", e2)
participant(r, "product", "gen3", "product", e2)
classification(r, e2, CC_FORM, CC_CHANGES)
bundle["reactions"].append(r)

r = st.new_reaction(f"{DOC}:scheme1C_proposed_intermediate", DOC, ev(1) + e2)
set_meta(r, "proposed_mechanism", [], "PDF pages 1-2, mechanistic rationale and Scheme 1C inset", ev(1) + e2,
         "η2-coordination-induced SN1-type ionization model", "proposed Ir-stabilized tertiary allenylic carbocation",
         ["Author-proposed intermediate/model; not treated as a demonstrated preparative experiment."])
participant(r, "substrate", "gen1", "reactant", ev(1) + e2)
participant(r, "intermediate", "ir_bound", "product", ev(1) + e2)
classification(r, ev(1) + e2, [], [], coverage="unresolved", limitations=["Coordination/ionization drawing is mechanistic and not a complete mapped reaction."])
bundle["reactions"].append(r)


def add_standard_alkylation_participants(r, evidence, nuc="2a", catalyst="ir", ligand="L1", lewis="La",
                                          scale=("0.4 mmol", {"number": 0.4, "unit": "mmol"}),
                                          nuc_equiv=2.5, la_loading=7.5, ligand_loading=20.0,
                                          catalyst_loading=5.0, product="3a"):
    participant(r, "substrate", "1a", "reactant", evidence,
                charges=[] if scale is None else [("amount", scale[0], scale[1], None)])
    participant(r, "nucleophile", nuc, "reactant", evidence,
                charges=[("equivalents", f"{nuc_equiv} equiv", {"number": nuc_equiv, "unit": "equiv"}, mid("1a"))])
    if product:
        participant(r, "product", product, "product", evidence)
    if catalyst:
        participant(r, "catalyst", catalyst, "catalyst", evidence,
                    charges=[("mol_percent", f"{catalyst_loading} mol %", {"number": catalyst_loading, "unit": "mol%"}, mid("1a"))])
    if ligand:
        charges = [] if ligand_loading is None else [("mol_percent", f"{ligand_loading} mol %", {"number": ligand_loading, "unit": "mol%"}, mid("1a"))]
        participant(r, "ligand", ligand, "ligand", evidence, charges=charges)
    if lewis:
        participant(r, "lewis_acid", lewis, "additive", evidence,
                    charges=[("mol_percent", f"{la_loading} mol %", {"number": la_loading, "unit": "mol%"}, mid("1a"))])
    participant(r, "solvent", "dioxane", "solvent", evidence)
    if scale is not None:
        r["reaction_scale"] = f(scale[0], {"reference_molecule_id": mid("1a"), "amount": scale[1]}, evidence)
    set_solvent(r, "dioxane", evidence, ("0.1 M", {"number": 0.1, "unit": "M"}), "1a")


# Table 1 optimization entries.
table1 = [
    (1, "none", "2a", "ir", "L1", "La", 7.5, 45, 83, 75, 99),
    (2, "2b instead of 2a", "2b", "ir", "L1", "La", 7.5, 45, 92, 88, 89),
    (3, "2c instead of 2a", "2c", "ir", "L1", "La", 7.5, 45, 80, None, 96),
    (4, "23 °C instead of 45 °C", "2a", "ir", "L1", "La", 7.5, 23, 85, None, 97),
    (5, "60 °C instead of 45 °C", "2a", "ir", "L1", "La", 7.5, 60, 58, None, 99),
    (6, "no La(OTf)3", "2a", "ir", "L1", None, 0, 45, 0, None, None),
    (7, "Zn(OTf)2 instead of La(OTf)3", "2a", "ir", "L1", "Zn", 7.5, 45, 73, None, 99),
    (8, "(S)-L2 instead of (S)-L1", "2a", "ir", "L2", "La", 7.5, 45, 0, None, None),
    (9, "no [Ir] nor (S)-L1", "2a", None, None, "La", 7.5, 45, 3, None, None),
    (10, "no [Ir]", "2a", None, "L1", "La", 7.5, 45, 1, None, None),
    (11, "no (S)-L1", "2a", "ir", None, "La", 7.5, 45, 1, None, None),
    (12, "1:1 [Ir]:L1", "2a", "ir", "L1", "La", 7.5, 45, 39, None, 78),
    (13, "[Rh] instead of [Ir]", "2a", "rh", "L1", "La", 7.5, 45, 3, None, None),
    (14, "add 4 (10 mol %)", "2a", "ir", "L1", "La", 7.5, 45, 73, None, 99),
]
for entry, variation, nuc, catalyst, ligand, lewis, la_load, temp, nmr_yield, isolated, ee in table1:
    evidence = ev(2)
    if entry == 1:
        evidence = ev(1) + ev(2)
    r = st.new_reaction(f"{DOC}:table1_entry_{entry:02d}", DOC, evidence)
    set_meta(r, "specific_experiment", ["optimization"], f"PDF page 2, Table 1 entry {entry}", evidence,
             "Effect of reaction parameters", "Ir-catalyzed allenylic C-C alkylation",
             ["Supporting Information was not supplied; workup, purification, atmosphere, and detailed addition order are unavailable."])
    ligand_loading = None if entry == 12 else 20.0
    add_standard_alkylation_participants(r, evidence, nuc, catalyst, ligand, lewis,
                                         ligand_loading=ligand_loading, la_loading=la_load,
                                         catalyst_loading=5.0)
    if entry == 14:
        participant(r, "acid_scavenger", "base4", "additive", evidence,
                    charges=[("mol_percent", "10 mol %", {"number": 10, "unit": "mol%"}, mid("1a"))])
    step(r, "reaction", evidence, [p["participant_id"] for p in r["participants"]],
         f"variation from standard conditions: {variation}",
         (f"{temp} °C", {"number": temp, "unit": "°C"}), ("17 h", {"number": 17, "unit": "h"}))
    r["procedure"]["original_text"] = f(variation, variation, evidence)
    outcome(r, "yield_nmr", "yield", ["3a"], evidence, f"{nmr_yield}%", {"number": nmr_yield, "unit": "%"},
            "Determined by 1H NMR analysis of the unpurified reaction mixture using (CHCl2)2 as an internal standard",
            "1H NMR of unpurified reaction mixture; (CHCl2)2 internal standard", "NMR yield")
    if isolated is not None:
        outcome(r, "yield_isolated", "yield", ["3a"], evidence, f"({isolated})", {"number": isolated, "unit": "%"},
                "isolated yield (parenthesized value, footnote d)", "isolated", "isolated yield")
    if ee is not None:
        outcome(r, "ee", "selectivity", ["3a"], evidence, f"{ee}% ee", {"number": ee, "unit": "%"},
                "Determined by HPLC with chiral stationary phase", "HPLC with chiral stationary phase", scope="enantiomeric excess")
    if entry == 1:
        outcome(r, "rr", "selectivity", ["3a"], evidence, ">20:1 selectivity over the corresponding 1,3-diene regioisomer",
                ">20:1", scope="regioselectivity, β-allenyl ester:1,3-diene")
    if entry == 6:
        outcome(r, "starting_material_recovery", "recovery", ["1a"], evidence,
                "only starting material was recovered", "only starting material recovered", product_role="recovered starting material")
    if entry in (6, 8):
        no_change_classification(r, evidence)
    else:
        classification(r, evidence, CC_FORM, CC_CHANGES)
    if entry == 12:
        r["notes"].append(f("1:1 [Ir]:L1", "Catalyst/ligand ratio 1:1; exact ligand mol% is not printed in the row and was not calculated.", evidence))
        r["review"]["unresolved"].append("Table 1 entry 12 gives a 1:1 [Ir]:L1 ratio but not a standalone ligand loading; no mol% was inferred.")
    bundle["reactions"].append(r)


# Table 2 electrophile scope, including the distinct 1 g scale experiment.
table2 = {
    "a": (76, 99, ">20:1", 17, []), "b": (77, 92, ">20:1", 20, ["b"]),
    "c": (73, 97, ">20:1", 7, []), "d": (65, 96, ">20:1", 20, []),
    "e": (76, 92, ">20:1", 14, ["b"]), "f": (68, 91, ">20:1", 27, ["b"]),
    "g": (52, 96, ">20:1", 20, ["b", "c"]), "h": (66, 93, ">20:1", 67, ["b", "d"]),
    "i": (75, 98, ">20:1", 23, ["e"]), "j": (82, 98, ">20:1", 17, ["e"]),
    "k": (82, 93, ">20:1", 24, ["b"]), "l": (76, 96, "10:1", 3, ["f"]),
    "m": (71, 94, ">20:1", 8, ["b"]), "n": (52, 90, ">20:1", 3, ["b"]),
    "o": (65, 93, "20:1", 27, []), "p": (70, 93, "7:1", 2, []),
    "q": (54, 94, ">20:1", 22, ["f"]), "r": (0, None, None, 14, []),
    "s": (52, 90, ">20:1", 53, ["b"]), "t": (55, 94, "20:1", 44, ["b"]),
}
for letter, (yield_value, ee, rr, hours, footnotes) in table2.items():
    evidence = ev(2) + ev(3) if letter == "a" else ev(3)
    r = st.new_reaction(f"{DOC}:table2_3{letter}", DOC, evidence)
    cats = ["substrate_scope"] + (["scale_up"] if letter == "a" else [])
    set_meta(r, "specific_experiment", cats, f"PDF page 3, Table 2 product 3{letter}", evidence,
             "Ir-catalyzed enantioselective alkylation: electrophile scope", "allenylic C-C alkylation",
             ["No SI was supplied; workup, purification, atmosphere, and addition order are unreported in the article."])
    substrate_key = "1a" if letter == "a" else f"sub_{letter}"
    nuc = "2b" if "b" in footnotes else "2a"
    scale = None if letter == "a" else ("0.4 mmol", {"number": 0.4, "unit": "mmol"})
    # Build manually so each scope entry uses its own electrophile molecule.
    participant(r, "substrate", substrate_key, "reactant", evidence,
                charges=[] if scale is None else [("amount", "0.4 mmol", {"number": 0.4, "unit": "mmol"}, None)])
    nuc_eq = 4.0 if "e" in footnotes else 2.5
    participant(r, "nucleophile", nuc, "reactant", evidence,
                charges=[("equivalents", f"{nuc_eq} equiv", {"number": nuc_eq, "unit": "equiv"}, mid(substrate_key))])
    participant(r, "product", f"3{letter}", "product", evidence)
    participant(r, "catalyst", "ir", "catalyst", evidence,
                charges=[("mol_percent", "5.0 mol %", {"number": 5.0, "unit": "mol%"}, mid(substrate_key))])
    participant(r, "ligand", "L1", "ligand", evidence,
                charges=[("mol_percent", "20 mol %", {"number": 20, "unit": "mol%"}, mid(substrate_key))])
    la = 20.0 if "c" in footnotes else 10.0 if "d" in footnotes else 7.5
    participant(r, "lewis_acid", "La", "additive", evidence,
                charges=[("mol_percent", f"{la} mol %", {"number": la, "unit": "mol%"}, mid(substrate_key))])
    participant(r, "solvent", "dioxane", "solvent", evidence)
    if letter == "a":
        r["reaction_scale"] = f("1 gram scale using substrate 1a", {"reference_molecule_id": mid("1a"), "mass": {"number": 1, "unit": "g"}}, evidence)
    else:
        r["reaction_scale"] = f("0.4 mmol", {"reference_molecule_id": mid(substrate_key), "amount": {"number": 0.4, "unit": "mmol"}}, evidence)
    set_solvent(r, "dioxane", evidence, ("0.1 M", {"number": 0.1, "unit": "M"}), substrate_key)
    temp = 23 if "f" in footnotes else 45
    step(r, "reaction", evidence, [p["participant_id"] for p in r["participants"]],
         f"Table 2 scope reaction for 3{letter}", (f"{temp} °C", {"number": temp, "unit": "°C"}),
         (f"{hours} h", {"number": hours, "unit": "h"}))
    outcome(r, "yield", "yield", [f"3{letter}"], evidence, f"{yield_value}% yield", {"number": yield_value, "unit": "%"},
            "Isolated yields shown", "isolated", "isolated yield")
    if ee is not None:
        outcome(r, "ee", "selectivity", [f"3{letter}"], evidence, f"{ee}% ee", {"number": ee, "unit": "%"},
                "HPLC, SFC, or GC analysis using a chiral stationary phase", "HPLC/SFC/GC with chiral stationary phase",
                scope="enantiomeric excess; individual analytical platform not assigned to this panel")
    if rr is not None:
        outcome(r, "rr", "selectivity", [f"3{letter}"], evidence, f"{rr} rr", rr,
                "1H NMR analysis of the unpurified reaction mixture", "1H NMR of unpurified reaction mixture",
                scope="regiochemical ratio; major β-allenyl ester:minor regioisomer")
    if letter == "r":
        no_change_classification(r, evidence)
        r["notes"].append(f("3r; 0% (14 h)", "Desired product 3r is drawn by the source but not observed; no expected bond label assigned.", evidence))
    else:
        classification(r, evidence, CC_FORM, CC_CHANGES)
    bundle["reactions"].append(r)


# Prose-only negative result for the tert-butyl ester-derived nucleophile.
e3 = ev(3)
r = st.new_reaction(f"{DOC}:table3_tbutyl_negative", DOC, e3)
set_meta(r, "specific_experiment", ["substrate_scope", "control"], "PDF page 3, Table 3 discussion", e3,
         "tert-Butyl silyl ketene acetal negative scope result", "attempted allenylic C-C alkylation",
         ["The article states only that the tert-butyl derivative gave none of the desired product; exact time, temperature, loadings, and detection method are not separately reported for this trial."])
participant(r, "substrate", "1a", "reactant", e3)
participant(r, "nucleophile", "2tbu", "reactant", e3)
step(r, "attempt", e3, ["substrate", "nucleophile"], "tert-butyl derivative was employed; none of the desired β-allenyl ester was observed")
outcome(r, "no_product", "no_reaction", [], e3, "none of the desired β-allenyl ester was observed", "not observed",
        basis_raw=None, basis=None, product_role="desired product not identified as a molecule")
no_change_classification(r, e3)
bundle["reactions"].append(r)


# Table 3 nucleophile scope.
table3 = {
    "u": ("2d", 70, 95, ">20:1", None, 36, 10.0),
    "v": ("2e", 72, 85, ">20:1", "2.4:1", 36, 7.5),
    "w": ("2f", 54, 88, ">20:1", None, 62, 7.5),
    "x": ("2g", 82, 91, ">20:1", None, 5, 7.5),
    "y": ("2h", 57, 96, ">20:1", None, 14, 7.5),
    "z": ("2i", 68, 92, ">20:1", None, 14, 7.5),
}
for letter, (nuc, yld, ee, rr, dr, hours, la) in table3.items():
    evidence = ev(3) + ev(4) if letter in ("u", "v", "w") else ev(4)
    r = st.new_reaction(f"{DOC}:table3_3{letter}", DOC, evidence)
    set_meta(r, "specific_experiment", ["substrate_scope"], f"PDF page 4, Table 3 product 3{letter}", evidence,
             "Ir-catalyzed enantioselective alkylation: nucleophile scope", "allenylic C-C alkylation",
             ["No SI was supplied; workup, purification, atmosphere, and addition order are unavailable."])
    add_standard_alkylation_participants(r, evidence, nuc=nuc, product=f"3{letter}", la_loading=la)
    step(r, "reaction", evidence, [p["participant_id"] for p in r["participants"]],
         f"Table 3 scope reaction for 3{letter}", ("45 °C", {"number": 45, "unit": "°C"}),
         (f"{hours} h", {"number": hours, "unit": "h"}))
    outcome(r, "yield", "yield", [f"3{letter}"], evidence, f"{yld}% yield", {"number": yld, "unit": "%"},
            "Isolated yields shown", "isolated", "isolated yield")
    outcome(r, "ee", "selectivity", [f"3{letter}"], evidence, f"{ee}% ee", {"number": ee, "unit": "%"},
            "HPLC or SFC analysis using a chiral stationary phase", "HPLC or SFC with chiral stationary phase",
            scope="enantiomeric excess; individual analytical platform not assigned to this panel")
    outcome(r, "rr", "selectivity", [f"3{letter}"], evidence, f"{rr} rr", rr,
            "1H NMR analysis of the unpurified reaction mixture", "1H NMR of unpurified reaction mixture",
            scope="regiochemical ratio; major β-allenyl product:minor regioisomer")
    if dr:
        outcome(r, "dr", "selectivity", [f"3{letter}"], evidence, f"{dr} dr", dr,
                scope="diastereomeric ratio; major:minor configuration not assigned; measurement method unreported")
        r["review"]["unresolved"].append("The 2.4:1 dr measurement method and full stereoisomer ordering for 3v are not reported on pages 1-6.")
    classification(r, evidence, CC_FORM, CC_CHANGES)
    bundle["reactions"].append(r)


# Page-7 footnote 79 is a distinct substrate-scope experiment, not a correction
# to the E-enriched 2e row.  The footnote reports outcomes but not the detailed
# conditions or analytical methods, so those facts remain explicitly missing.
e7 = ev(7)
e47 = ev(4) + e7
r = st.new_reaction(f"{DOC}:footnote79_z_2e_3v", DOC, e47)
set_meta(
    r, "specific_experiment", ["substrate_scope"], "PDF page 7, footnote 79 (context: page 4, Table 3)", e47,
    "Ir-catalyzed enantioselective alkylation using (Z)-2e",
    "allenylic C-C alkylation",
    [
        "Footnote 79 reports use of the Z isomer and four outcomes but does not separately state reaction time, temperature, scale, catalyst/ligand/additive loadings, solvent, atmosphere, workup, purification, or addition order.",
        "Yield type/measurement basis and the analytical methods for ee, rr, and dr are not stated in footnote 79 and were not inherited from the adjacent Table 3 row.",
        "The major product is stated to be the same diastereomer as in the E-enriched 2e experiment, but its absolute stereocenter descriptors and the minor diastereomer configuration are not printed in footnote 79.",
    ],
)
participant(r, "substrate", "1a", "reactant", e47)
participant(r, "nucleophile", "2e_z", "reactant", e47)
participant(r, "product", "3v", "product", e47)
step(
    r, "reported_trial", e7, ["substrate", "nucleophile", "product"],
    "When the (Z)-isomer of silyl ketene acetal was used, the same diastereomer was obtained as a major product",
)
r["procedure"]["original_text"] = f(
    "When the (Z)-isomer of silyl ketene acetal was used, the same diastereomer was obtained as a major product in 73% yield with 83% ee, >20:1 rr, and 1.7:1 dr.",
    "When the (Z)-isomer of silyl ketene acetal was used, the same diastereomer was obtained as a major product in 73% yield with 83% ee, >20:1 rr, and 1.7:1 dr.",
    e7,
)
outcome(r, "yield", "yield", ["3v"], e7, "73% yield", {"number": 73, "unit": "%"},
        yield_type="reported yield")
outcome(r, "ee", "selectivity", ["3v"], e7, "83% ee", {"number": 83, "unit": "%"},
        scope="enantiomeric excess; major enantiomer/absolute configuration and measurement method unreported")
outcome(r, "rr", "selectivity", ["3v"], e7, ">20:1 rr", ">20:1",
        scope="regiochemical ratio; ordering inferred from the linked 3v experiment context, measurement method unreported")
outcome(r, "dr", "selectivity", ["3v"], e7, "1.7:1 dr", "1.7:1",
        scope="diastereomeric ratio; same diastereomer as the E-enriched 2e experiment is major, minor configuration and measurement method unreported")
classification(r, e47, CC_FORM, CC_CHANGES)
bundle["reactions"].append(r)


def add_scheme2_common(r, evidence, reactant, product, name, transformation, unresolved):
    set_meta(r, "specific_experiment", ["derivatization"], "PDF page 4, Scheme 2", evidence, name, transformation, unresolved)
    participant(r, "substrate", reactant, "reactant", evidence)
    participant(r, "product", product, "product", evidence)


e4 = ev(3) + ev(4)
# 3a -> 5 hydrosilylation.
r = st.new_reaction(f"{DOC}:scheme2_hydrosilylation_5", DOC, e4)
add_scheme2_common(r, e4, "3a", "5", "Mo-catalyzed regioselective hydrosilylation", "allene hydrosilylation",
                   ["Reaction time, atmosphere, workup, purification, and yield measurement basis are not reported in the main article."])
participant(r, "catalyst", "moco6", "catalyst", e4, charges=[("mol_percent", "10 mol %", {"number": 10, "unit": "mol%"}, mid("3a"))])
participant(r, "hydrosilane", "ph2sih2", "reagent", e4, charges=[("equivalents", "1.2 equiv", {"number": 1.2, "unit": "equiv"}, mid("3a"))])
participant(r, "solvent", "phme", "solvent", e4)
set_solvent(r, "phme", e4)
step(r, "reaction", e4, [p["participant_id"] for p in r["participants"]], "Mo-catalyzed hydrosilylation", ("110 °C", {"number": 110, "unit": "°C"}))
outcome(r, "yield", "yield", ["5"], e4, "73% yield", {"number": 73, "unit": "%"}, yield_type="reported yield")
outcome(r, "rr", "selectivity", ["5"], e4, ">20:1 rr", ">20:1", scope="regiochemical ratio")
outcome(r, "ze", "selectivity", ["5"], e4, ">20:1 Z/E", ">20:1", scope="Z:E; product 5")
outcome(r, "ee", "selectivity", ["5"], e4, "99% ee", {"number": 99, "unit": "%"}, scope="enantiomeric excess")
changes = [
    {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "Si"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
    {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
    {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
]
labels = [{"category_number": 12, "code": "C-Si", "supporting_change_ids": ["b1"]}, {"category_number": 1, "code": "C-H", "supporting_change_ids": ["b2"]}]
classification(r, e4, labels, changes)
bundle["reactions"].append(r)

# 3a -> 6 arylamination.
r = st.new_reaction(f"{DOC}:scheme2_arylamination_6", DOC, e4)
add_scheme2_common(r, e4, "3a", "6", "Pd-catalyzed arylamination", "allene arylamination",
                   ["Reaction time, atmosphere, workup, purification, and yield measurement basis are not reported in the main article."])
for pid, key, role, charge in [
    ("catalyst", "pdoac2", "catalyst", ("mol_percent", "2.5 mol %", {"number": 2.5, "unit": "mol%"})),
    ("ligand", "dppf", "ligand", ("mol_percent", "5 mol %", {"number": 5, "unit": "mol%"})),
    ("amine", "pyrrolidine", "reagent", ("equivalents", "6 equiv", {"number": 6, "unit": "equiv"})),
    ("aryl", "phi", "reagent", ("equivalents", "1.25 equiv", {"number": 1.25, "unit": "equiv"})),
]:
    participant(r, pid, key, role, e4, charges=[(charge[0], charge[1], charge[2], mid("3a"))])
participant(r, "solvent", "mecn", "solvent", e4)
set_solvent(r, "mecn", e4)
step(r, "reaction", e4, [p["participant_id"] for p in r["participants"]], "Pd-catalyzed arylamination", ("90 °C", {"number": 90, "unit": "°C"}))
outcome(r, "yield", "yield", ["6"], e4, "70% yield", {"number": 70, "unit": "%"}, yield_type="reported yield")
outcome(r, "rr", "selectivity", ["6"], e4, ">20:1 rr", ">20:1", scope="regiochemical ratio")
outcome(r, "ze", "selectivity", ["6"], e4, ">20:1 Z/E", ">20:1", scope="Z:E; product 6")
outcome(r, "ee", "selectivity", ["6"], e4, "99% ee", {"number": 99, "unit": "%"}, scope="enantiomeric excess")
changes = [
    {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "N"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
    {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
]
labels = [{"category_number": 5, "code": "C-N", "supporting_change_ids": ["b1"]}, {"category_number": 2, "code": "C-C", "supporting_change_ids": ["b2"]}]
classification(r, e4, labels, changes)
bundle["reactions"].append(r)

# 3a -> 7, three-operation ester-to-acid sequence.  The two warming programs are
# represented as separate start/warm steps so temperature and duration remain paired.
r = st.new_reaction(f"{DOC}:scheme2_acid_7", DOC, e4)
add_scheme2_common(r, e4, "3a", "7", "Ester-to-carboxylic-acid sequence", "ester reduction/oxidation to carboxylic acid",
                   ["No product yield, atmosphere, workup, purification, or detailed addition order is reported for acid 7 in the main article."])
for pid, key, role, eq in [
    ("dibal", "dibal", "reagent", 3.0), ("dmp", "dmp", "reagent", 1.5),
    ("naclo2", "naclo2", "reagent", 3.0), ("kh2po4", "kh2po4", "reagent", 3.0),
    ("alkene", "methylbutene", "reagent", 6.6),
]: participant(r, pid, key, role, e4, charges=[("equivalents", f"{eq} equiv", {"number": eq, "unit": "equiv"}, mid("3a"))])
for pid, key in [("dcm", "dcm"), ("tbuoh", "tbuoh"), ("water", "water")]:
    participant(r, pid, key, "solvent", e4)
r["solvent_conditions"]["composition"] = f("DCM; tBuOH/H2O", [mid("dcm"), mid("tbuoh"), mid("water")], e4)
r["solvent_conditions"]["mixture_ratio"] = f("tBuOH/H2O", "ratio not reported", e4)
allp = [p["participant_id"] for p in r["participants"]]
step(r, "a_low_temperature", e4, allp, "DIBAL-H treatment begins at -78 °C", ("-78 °C", {"number": -78, "unit": "°C"}))
step(r, "a_warm", e4, allp, "warm from -78 °C to room temperature", ("rt", "r.t."), ("30 min", {"number": 30, "unit": "min"}))
step(r, "b_oxidation", e4, allp, "DMP oxidation in DCM", ("0 °C", {"number": 0, "unit": "°C"}), ("30 min", {"number": 30, "unit": "min"}))
step(r, "c_low_temperature", e4, allp, "NaClO2/KH2PO4/2-methyl-2-butene treatment begins at 0 °C", ("0 °C", {"number": 0, "unit": "°C"}))
step(r, "c_warm", e4, allp, "warm from 0 °C to room temperature", ("rt", "r.t."), ("11 h", {"number": 11, "unit": "h"}))
outcome(r, "ee", "selectivity", ["7"], e4, "98% ee", {"number": 98, "unit": "%"}, scope="enantiomeric excess; measurement method not reported")
classification(r, e4, [], [], coverage="unresolved", limitations=["Multistep atom correspondence is not established; no formed-bond label is asserted."])
bundle["reactions"].append(r)

# 7 -> 8 Curtius rearrangement/Cbz capture.
r = st.new_reaction(f"{DOC}:scheme2_curtius_8", DOC, e4)
add_scheme2_common(r, e4, "7", "8", "Curtius rearrangement followed by BnOH quench", "Curtius rearrangement/carbamate formation",
                   ["Both stage durations, atmosphere, workup, purification, and yield measurement basis are unreported."])
participant(r, "dppa", "dppa", "reagent", e4, charges=[("equivalents", "1.6 equiv", {"number": 1.6, "unit": "equiv"}, mid("7"))])
participant(r, "et3n_stage1", "et3n", "reagent", e4, charges=[("equivalents", "3 equiv", {"number": 3, "unit": "equiv"}, mid("7"))])
participant(r, "bnoh", "bnoh", "reagent", e4, charges=[("equivalents", "7 equiv", {"number": 7, "unit": "equiv"}, mid("7"))])
participant(r, "et3n_stage2", "et3n", "reagent", e4, charges=[("equivalents", "7 equiv", {"number": 7, "unit": "equiv"}, mid("7"))])
participant(r, "solvent", "phme", "solvent", e4)
set_solvent(r, "phme", e4)
step(r, "curtius", e4, ["substrate", "dppa", "et3n_stage1", "solvent"], "Curtius rearrangement", ("80 °C", {"number": 80, "unit": "°C"}))
step(r, "benzyl_alcohol_capture", e4, ["product", "bnoh", "et3n_stage2", "solvent"], "quench/capture with BnOH", ("80 °C", {"number": 80, "unit": "°C"}))
outcome(r, "yield", "yield", ["8"], e4, "51% yield", {"number": 51, "unit": "%"}, yield_type="reported yield")
outcome(r, "ee", "selectivity", ["8"], e4, "98% ee", {"number": 98, "unit": "%"}, scope="enantiomeric excess")
changes = [{"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "N"}], "before_order": 0, "after_order": 1, "change_type": "formed"}]
classification(r, e4, [{"category_number": 5, "code": "C-N", "supporting_change_ids": ["b1"]}], changes)
bundle["reactions"].append(r)

# 7 -> 9 Rh-catalyzed lactonization.
r = st.new_reaction(f"{DOC}:scheme2_lactonization_9", DOC, e4)
add_scheme2_common(r, e4, "7", "9", "Rh-catalyzed intramolecular cyclization", "allene-carboxylic-acid lactonization",
                   ["Reaction time, atmosphere, workup, purification, and yield/dr/ee measurement methods are unreported."])
participant(r, "catalyst", "rh", "catalyst", e4, charges=[("mol_percent", "5 mol %", {"number": 5, "unit": "mol%"}, mid("7"))])
participant(r, "ligand", "dpephos", "ligand", e4, charges=[("mol_percent", "2.5 mol %", {"number": 2.5, "unit": "mol%"}, mid("7"))])
participant(r, "solvent", "dce", "solvent", e4)
set_solvent(r, "dce", e4)
step(r, "reaction", e4, [p["participant_id"] for p in r["participants"]], "Rh-catalyzed intramolecular cyclization", ("70 °C", {"number": 70, "unit": "°C"}))
outcome(r, "yield", "yield", ["9"], e4, "94% yield", {"number": 94, "unit": "%"}, yield_type="reported yield")
outcome(r, "dr", "selectivity", ["9"], e4, "5.2:1 dr", "5.2:1", scope="diastereomeric ratio; (R,S)-9 shown, minor configuration unreported")
outcome(r, "ee", "selectivity", ["9"], e4, "98% ee", {"number": 98, "unit": "%"}, scope="enantiomeric excess")
changes = [
    {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
    {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
]
labels = [{"category_number": 3, "code": "C-O", "supporting_change_ids": ["b1"]}, {"category_number": 1, "code": "C-H", "supporting_change_ids": ["b2"]}]
classification(r, e4, labels, changes)
bundle["reactions"].append(r)


# Non-synthetic source context from the reviewed pages.
bundle["source_context"] = [
    {"context_id": f"{DOC}:context_title", "document_id": DOC, "evidence_ids": ev(1),
     "kind": "bibliographic_metadata", "value": {"title": "Construction of Vicinal Quaternary Centers via Iridium-Catalyzed Asymmetric Allenylic Alkylation of Racemic Tertiary Alcohols", "doi": "10.1021/jacs.1c00609"}},
    {"context_id": f"{DOC}:context_scheme1A", "document_id": DOC, "evidence_ids": ev(2),
     "kind": "conceptual_figure", "value": "Scheme 1A compares quaternary-center motifs; it is not a reaction experiment."},
    {"context_id": f"{DOC}:context_si_notice", "document_id": DOC, "evidence_ids": ev(4),
     "kind": "supplement_notice", "value": "The article states Supporting Information exists, but no SI file was supplied or associated in this workspace."},
    {"context_id": f"{DOC}:context_pages5_6", "document_id": DOC, "evidence_ids": ev(5) + ev(6),
     "kind": "non_experimental_pages", "value": "Pages 5-6 contain author information, acknowledgments, and references; no current-paper experimental arrow/table/spectrum was identified."},
    {"context_id": f"{DOC}:context_page7_references", "document_id": DOC, "evidence_ids": ev(7),
     "kind": "reference_page", "value": "Page 7 contains references 65-82 and no reaction drawing, table, or standalone characterization table. Footnote 79 is the only current-paper experimental result on the page and is represented by reaction footnote79_z_2e_3v."},
]

(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"])}, ensure_ascii=False))
