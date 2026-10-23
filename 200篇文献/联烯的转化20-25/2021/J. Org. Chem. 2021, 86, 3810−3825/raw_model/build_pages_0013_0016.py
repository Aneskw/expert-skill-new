#!/usr/bin/env python3
"""Build evidence-linked v3.1 records for original PDF pages 13--16."""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as sr

RUN = ROOT / "runs/current"
DOC = "doc_0341d857f4f70e0f"
DATASET = RUN / "dataset.json"
OLD_STANDARD = RUN / "standard_records.json"
SUBMITTED = RUN / "submitted.json"


def dump(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


dataset = json.loads(DATASET.read_text())
old = json.loads(OLD_STANDARD.read_text())
bundle = sr.empty(dataset["documents"])
molecules = {m["molecule_id"]: copy.deepcopy(m) for m in old["molecules"]}
old_reactions = {r["reaction_id"]: r for r in old["reactions"]}


def ev(page: int) -> list[str]:
    return [f"{DOC}:page{page}:text", f"{DOC}:page{page}:image"]


def f(raw=None, normalized=None, evidence=(), method="original_extraction",
      derivation=None, status=None):
    if normalized is None and raw is not None:
        normalized = raw
    return sr.fact(raw=raw, normalized=normalized, evidence=evidence, method=method,
                   derivation=derivation, status=status)


def mid(key: str) -> str:
    return f"{DOC}:mol:{key}"


def append_evidence(record: dict, page: int):
    for evidence_id in ev(page):
        if evidence_id not in record["evidence_ids"]:
            record["evidence_ids"].append(evidence_id)


def add_molecule(key: str, label: str, name: str, *, chirality="unresolved",
                 geometry="unresolved", representation="unresolved",
                 identity_note: str | None = None) -> str:
    molecule_id = mid(key)
    if molecule_id in molecules:
        return molecule_id
    evidence = ev(13)
    m = sr.new_molecule(molecule_id, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    m["names"]["original"] = f(name, name, evidence)
    if identity_note:
        m["identity_note"] = f(identity_note, identity_note, evidence)
    reason = (
        "Assigned after inspection of the page-13 name/label and shown chemical role; "
        "no unreported stereochemistry was converted to achiral or racemic."
    )
    for field, value in (("chirality_status", chirality),
                         ("geometry_status", geometry),
                         ("representation_status", representation)):
        m["stereochemistry"][field] = f(
            None, value, evidence, "inference", reason
        )
    molecules[molecule_id] = m
    return molecule_id


# Page-13 workup materials and generic labels. No guessed concrete structure is assigned.
add_molecule(
    "sodium_thiosulfate_aqueous", "aqueous sodium thiosulfate",
    "aqueous sodium thiosulfate", chirality="not_applicable",
    geometry="not_applicable", representation="unresolved",
    identity_note="Concentration, hydration state, and amount are unreported."
)
add_molecule(
    "Celite", "Celite", "Celite filtration aid", chirality="not_applicable",
    geometry="not_applicable", representation="unresolved",
    identity_note="Heterogeneous filtration material; no molecular SMILES applies."
)
add_molecule(
    "10_general", "10", "generic 1-aroyl-beta-carboline 10",
    identity_note="Generic source label spanning 10a and 10b; no concrete SMILES assigned."
)
add_molecule(
    "11_general", "11", "generic debenzylated 1-aroyl-beta-carboline 11",
    identity_note="Generic source label spanning 11a and 11b; no concrete SMILES assigned."
)


def set_named_achiral(key: str, name: str):
    m = molecules[mid(key)]
    evidence = ev(13)
    m["names"]["original"] = f(name, name, evidence)
    reason = (
        "The complete page-13 product name and the source structure identify no stereogenic "
        "center, stereogenic axis, E/Z alkene, or ring cis/trans assignment; no SMILES was "
        "encoded because source-level structure transcription remains machine-unverified."
    )
    for field, value in (("chirality_status", "achiral"),
                         ("geometry_status", "not_applicable"),
                         ("representation_status", "not_encoded")):
        m["stereochemistry"][field] = f(
            None, value, evidence, "inference", reason
        )


set_named_achiral(
    "10a", "(9-Benzyl-9H-pyrido[3,4-b]indol-1-yl)(phenyl)methanone"
)
set_named_achiral(
    "10b", "(9-Benzyl-9H-pyrido[3,4-b]indol-1-yl)(4-methoxyphenyl)methanone"
)
set_named_achiral(
    "11a", "Phenyl(9H-pyrido[3,4-b]indol-1-yl)methanone"
)
set_named_achiral(
    "11b", "(4-Methoxyphenyl)(9H-pyrido[3,4-b]indol-1-yl)methanone"
)


def participant(r: dict, molecule_id: str, role: str, pid: str, *, stage="reaction",
                charges=None, reference=None, method="original_extraction"):
    p = sr.new_participant(pid, molecule_id, role, ev(13), stage)
    for field, (raw, normalized) in (charges or {}).items():
        p["charges"][field] = f(raw, normalized, ev(13), method)
    if reference:
        p["equivalence_reference"] = f(reference, reference, ev(13), method)
    r["participants"].append(p)
    return p


def step(r: dict, step_id: str, stage: str, participant_ids: list[str], operation: str,
         *, temperature=None, time=None):
    s = sr.new_step(step_id, ev(13), stage)
    s["participant_ids"] = participant_ids
    s["operation"] = f(operation, operation, ev(13))
    if temperature:
        s["temperature"] = f(temperature[0], temperature[1], ev(13))
    if time:
        s["time"] = f(time[0], time[1], ev(13))
    r["steps"].append(s)
    return s


def product_amount(r: dict, oid: str, target: str, raw: str, number: float):
    o = sr.new_outcome(oid, "product_amount", [target], ev(13))
    o["product_role"] = f("product", "product", ev(13))
    o["value"] = f(raw, {"number": number, "unit": "mg"}, ev(13))
    o["measurement_basis"] = f("isolated product", "isolated product", ev(13))
    r["outcomes"].append(o)


def reaction(suffix: str, name: str, transformation: str) -> dict:
    r = sr.new_reaction(f"{DOC}:{suffix}", DOC, ev(13))
    r["record_type"] = "generic_scheme"
    r["reaction_name"] = f(name, name, ev(13))
    r["transformation_type"] = f(transformation, transformation, ev(13))
    r["location"] = f("PDF page 13", {"pdf_page": 13}, ev(13))
    r["review"]["unresolved"] = [
        "Generic procedure remains needs_review; schema validity does not establish chemical correctness."
    ]
    r["completeness"] = "partial"
    return r


def classify(r: dict, mode: str):
    if mode == "carbonyl":
        change = {
            "change_id": "b1", "atoms": [{"atom_map": None, "element": "C"},
                                             {"atom_map": None, "element": "O"}],
            "before_order": 0, "after_order": 2, "change_type": "formed"
        }
        labels = [{"category_number": 14, "code": "DB-C-O",
                   "supporting_change_ids": ["b1"]}]
    else:
        change = {
            "change_id": "b1", "atoms": [{"atom_map": None, "element": "N"},
                                             {"atom_map": None, "element": "H"}],
            "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain",
            "hydrogen_count": 1
        }
        labels = [{"category_number": 29, "code": "N-H",
                   "supporting_change_ids": ["b1"]}]
    r["classification"]["bond_changes"] = [f(
        change, change, ev(13), "structure_derivation",
        "Partial net bond change inferred from the page-13 named/drawn transformation without atom mapping."
    )]
    r["classification"]["labels"] = f(
        labels, labels, ev(13), "structure_derivation",
        "Taxonomy label inferred from the net source transformation; mapping remains incomplete."
    )
    r["classification"]["coverage"] = f(
        None, "partial", ev(13), "inference",
        "Only the clearly indicated new bond is classified; no atom-mapped reaction SMILES was supplied."
    )
    r["classification"]["limitations"] = [
        "No atom mapping; bond cleavage, leaving groups, and other implicit-H changes may be incomplete."
    ]


# Generic page-13 oxidative-cleavage procedure.
gen10 = reaction(
    "general_synthesis_aroyl_beta_carbolines_10",
    "General procedure for 1-aroyl-beta-carbolines 10a,b",
    "oxidative cleavage of exocyclic alkene to aroyl beta-carboline"
)
ref2 = mid("2_general")
participant(gen10, ref2, "reactant", "p_2", charges={
    "amount": ("0.05 mmol", {"number": 0.05, "unit": "mmol"})
})
participant(gen10, mid("RuCl3"), "catalyst", "p_RuCl3", charges={
    "mass": ("0.5 mg", {"number": 0.5, "unit": "mg"}),
    "amount": ("0.0025 mmol", {"number": 0.0025, "unit": "mmol"}),
    "mol_percent": ("5 mol%", {"number": 5, "unit": "mol%"})
}, reference=ref2)
participant(gen10, mid("NaIO4"), "reagent", "p_NaIO4", charges={
    "mass": ("64.2 mg", {"number": 64.2, "unit": "mg"}),
    "amount": ("0.3 mmol", {"number": 0.3, "unit": "mmol"}),
    "equivalents": ("6 equiv", {"number": 6, "unit": "equiv"})
}, reference=ref2)
for key, pid in (("MeCN", "p_MeCN"), ("CCl4", "p_CCl4"), ("water", "p_water_solvent")):
    participant(gen10, mid(key), "solvent", pid)
participant(gen10, mid("10_general"), "product", "p_10")
gen10["reaction_scale"] = f(
    "2e-f (0.05 mmol)",
    {"reference_molecule_id": ref2, "amount": {"number": 0.05, "unit": "mmol"}}, ev(13)
)
gen10["solvent_conditions"]["composition"] = f(
    "CH3CN/CCl4/H2O", ["CH3CN", "CCl4", "H2O"], ev(13)
)
gen10["solvent_conditions"]["mixture_ratio"] = f("2:1:1", [2, 1, 1], ev(13))
gen10["special_conditions"]["reactor"] = f("sealed tube", "sealed tube", ev(13))
step(gen10, "reaction", "reaction", ["p_2", "p_RuCl3", "p_NaIO4", "p_MeCN", "p_CCl4", "p_water_solvent", "p_10"],
     "stirred in a sealed tube and monitored by TLC",
     temperature=("rt", "room temperature"), time=("1 h", {"number": 1, "unit": "h"}))
participant(gen10, mid("ethyl_acetate"), "workup_agent", "p_EtOAc_workup", stage="workup",
            charges={"volume": ("10 mL", {"number": 10, "unit": "mL"})})
participant(gen10, mid("sodium_thiosulfate_aqueous"), "workup_agent", "p_thiosulfate", stage="workup")
participant(gen10, mid("Celite"), "workup_agent", "p_Celite", stage="workup")
participant(gen10, mid("water"), "workup_agent", "p_water_wash", stage="workup")
participant(gen10, mid("brine"), "workup_agent", "p_brine", stage="workup")
participant(gen10, mid("sodium_sulfate"), "workup_agent", "p_Na2SO4", stage="workup")
participant(gen10, mid("silica_gel"), "purification_agent", "p_silica", stage="purification")
participant(gen10, mid("ethyl_acetate"), "purification_agent", "p_EtOAc_eluent", stage="purification")
participant(gen10, mid("petroleum_ether"), "purification_agent", "p_petroleum_ether", stage="purification")
step(gen10, "dilute_and_quench", "workup", ["p_EtOAc_workup", "p_thiosulfate"],
     "diluted with ethyl acetate and quenched with aqueous sodium thiosulfate")
step(gen10, "filter", "workup", ["p_Celite"], "filtered through a plug of Celite")
step(gen10, "wash_dry_concentrate", "workup", ["p_water_wash", "p_brine", "p_Na2SO4"],
     "filtrate was washed with water and brine; ethyl acetate extracts were dried and concentrated")
step(gen10, "purification", "purification", ["p_silica", "p_EtOAc_eluent", "p_petroleum_ether"],
     "silica-gel column chromatography (100-200 mesh) using 15% ethyl acetate-petroleum ether (v/v)")
gen10["procedure"]["original_text"] = f(
    "2e-f (0.05 mmol), RuCl3 (0.5 mg, 0.0025 mmol), and NaIO4 (64.2 mg, 0.3 mmol) in CH3CN/CCl4/H2O (2:1:1) were stirred in a sealed tube at rt for 1 h; dilute with EtOAc (10 mL), quench with aqueous sodium thiosulfate, filter through Celite, wash, dry, concentrate, and chromatograph.",
    "Page-13 general procedure for 10a,b", ev(13)
)
gen10["notes"].append(f(
    "The source reports a 40-45% range for products 10 collectively; the range is not assigned as a generic-experiment yield because product-specific yields are recorded on 10a and 10b.",
    "General range only; entry-specific yields remain on 10a/10b", ev(13)
))
gen10["review"]["unresolved"].extend([
    "Individual volumes of CH3CN, CCl4, and H2O are unreported; only the 2:1:1 ratio is stated.",
    "Amounts of aqueous sodium thiosulfate, Celite, wash water, brine, drying agent, and chromatography solvents are unreported except for the 10 mL dilution EtOAc."
])
classify(gen10, "carbonyl")
bundle["reactions"].append(gen10)


# Generic page-13 debenzylation procedure.
gen11 = reaction(
    "general_synthesis_debenzylated_beta_carbolines_11",
    "General procedure for debenzylated 1-aroyl-beta-carbolines 11a,b",
    "AlCl3-mediated N-debenzylation"
)
ref10 = mid("10_general")
participant(gen11, ref10, "reactant", "p_10", charges={
    "amount": ("0.02 mmol", {"number": 0.02, "unit": "mmol"})
})
participant(gen11, mid("AlCl3"), "reagent", "p_AlCl3", charges={
    "amount": ("0.06 mmol", {"number": 0.06, "unit": "mmol"}),
    "equivalents": ("3.0 equiv", {"number": 3.0, "unit": "equiv"})
}, reference=ref10)
participant(gen11, mid("toluene"), "solvent", "p_toluene", charges={
    "volume": ("1 mL", {"number": 1, "unit": "mL"})
})
participant(gen11, mid("11_general"), "product", "p_11")
gen11["reaction_scale"] = f(
    "1-aroyl-beta-carboline 10 (0.02 mmol)",
    {"reference_molecule_id": ref10, "amount": {"number": 0.02, "unit": "mmol"}}, ev(13)
)
gen11["solvent_conditions"]["composition"] = f("dry toluene", ["toluene"], ev(13))
step(gen11, "reaction", "reaction", ["p_10", "p_AlCl3", "p_toluene", "p_11"],
     "stirred in dry toluene", temperature=("rt", "room temperature"),
     time=("5 h", {"number": 5, "unit": "h"}))
gen11["procedure"]["original_text"] = f(
    "1-Aroyl-beta-carboline 10 (0.02 mmol) and anhydrous AlCl3 (0.06 mmol) in dry toluene (1 mL) were stirred at rt for 5 h.",
    "Page-13 general procedure for 11a,b", ev(13)
)
gen11["notes"].append(f(
    "The source reports an 80-85% range for 11a,b collectively; the range is not assigned as a generic-experiment yield because product-specific yields are recorded on 11a and 11b.",
    "General range only; entry-specific yields remain on 11a/11b", ev(13)
))
gen11["review"]["unresolved"].append(
    "The page-13 general paragraph does not report workup details for the AlCl3 step; no workup was inferred."
)
classify(gen11, "nh")
bundle["reactions"].append(gen11)


# Existing Scheme 3 scope experiments: attach characterization-derived isolated masses and bases.
scheme3_data = {
    "9a": (27.0, "combined yield of isolated inseparable diastereomers", "30% ethyl acetate-petroleum ether (v/v)"),
    "9b": (31.6, "isolated product", "50% ethyl acetate-petroleum ether (v/v)"),
    "9i": (28.9, "isolated product", "35% ethyl acetate-petroleum ether (v/v)"),
    "9g": (24.8, "combined yield of isolated inseparable diastereomers", "30% ethyl acetate-petroleum ether (v/v)"),
}
for label, (mass, basis, eluent) in scheme3_data.items():
    r = copy.deepcopy(old_reactions[f"{DOC}:scheme3_{label}"])
    append_evidence(r, 13)
    r["location"] = f("PDF pages 6 and 13", {"pdf_pages": [6, 13]}, ev(6) + ev(13))
    y = next(o for o in r["outcomes"] if o["outcome_id"] == "yield")
    y["yield_type"] = f("isolated", "isolated", ev(13))
    y["measurement_basis"] = f(basis, basis, ev(13))
    product_amount(r, f"amount_{label}", mid(label), f"{mass} mg", mass)
    r["notes"].append(f(
        f"Product {label} was isolated using {eluent}; NMR/HRMS peak lists are characterization context, not separate experiments.",
        {"product": label, "isolation_eluent": eluent}, ev(13)
    ))
    r["review"]["unresolved"] = [
        u for u in r["review"]["unresolved"]
        if "Measurement basis" not in u and "measurement basis" not in u
    ]
    bundle["reactions"].append(r)


def set_charge_by_molecule(r: dict, molecule_key: str, field: str, raw, normalized):
    p = next(p for p in r["participants"] if p["molecule_id"] == mid(molecule_key))
    p["charges"][field] = f(raw, normalized, ev(13))
    return p["participant_id"]


# Existing Scheme 5 oxidative-cleavage experiments: exact scale, workup, and product mass.
for substrate, product, mass in (("2e", "10a", 8.2), ("2f", "10b", 7.8)):
    r = copy.deepcopy(old_reactions[f"{DOC}:scheme5_{substrate}_to_{product}"])
    append_evidence(r, 13)
    r["location"] = f("PDF pages 6 and 13", {"pdf_pages": [6, 13]}, ev(6) + ev(13))
    r["general_scheme_ids"] = list(dict.fromkeys(
        r["general_scheme_ids"] + [gen10["reaction_id"]]
    ))
    set_charge_by_molecule(r, substrate, "amount", "0.05 mmol", {"number": 0.05, "unit": "mmol"})
    set_charge_by_molecule(r, "RuCl3", "mass", "0.5 mg", {"number": 0.5, "unit": "mg"})
    set_charge_by_molecule(r, "RuCl3", "amount", "0.0025 mmol", {"number": 0.0025, "unit": "mmol"})
    set_charge_by_molecule(r, "NaIO4", "mass", "64.2 mg", {"number": 64.2, "unit": "mg"})
    set_charge_by_molecule(r, "NaIO4", "amount", "0.3 mmol", {"number": 0.3, "unit": "mmol"})
    r["reaction_scale"] = f(
        f"{substrate} (0.05 mmol)",
        {"reference_molecule_id": mid(substrate), "amount": {"number": 0.05, "unit": "mmol"}}, ev(13)
    )
    r["special_conditions"]["reactor"] = f("sealed tube", "sealed tube", ev(13))
    r["solvent_conditions"]["composition"] = f(
        "CH3CN/CCl4/H2O", ["CH3CN", "CCl4", "H2O"], ev(13)
    )
    r["solvent_conditions"]["mixture_ratio"] = f("2:1:1", [2, 1, 1], ev(13))
    # Workup/purification participants are inherited from the page-13 general procedure.
    for molecule_key, role, pid, stage_name, charges in [
        ("ethyl_acetate", "workup_agent", "p_EtOAc_workup", "workup", {"volume": ("10 mL", {"number": 10, "unit": "mL"})}),
        ("sodium_thiosulfate_aqueous", "workup_agent", "p_thiosulfate", "workup", None),
        ("Celite", "workup_agent", "p_Celite", "workup", None),
        ("water", "workup_agent", "p_water_wash", "workup", None),
        ("brine", "workup_agent", "p_brine", "workup", None),
        ("sodium_sulfate", "workup_agent", "p_Na2SO4", "workup", None),
        ("silica_gel", "purification_agent", "p_silica", "purification", None),
        ("ethyl_acetate", "purification_agent", "p_EtOAc_eluent", "purification", None),
        ("petroleum_ether", "purification_agent", "p_petroleum_ether", "purification", None),
    ]:
        participant(r, mid(molecule_key), role, pid, stage=stage_name,
                    charges=charges, method="inherited_general_condition")
    reaction_pids = [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"]
    r["steps"] = []
    step(r, "reaction", "reaction", reaction_pids,
         "stirred in a sealed tube and monitored by TLC",
         temperature=("rt", "room temperature"), time=("1 h", {"number": 1, "unit": "h"}))
    step(r, "dilute_and_quench", "workup", ["p_EtOAc_workup", "p_thiosulfate"],
         "diluted with ethyl acetate and quenched with aqueous sodium thiosulfate")
    step(r, "filter", "workup", ["p_Celite"], "filtered through a plug of Celite")
    step(r, "wash_dry_concentrate", "workup", ["p_water_wash", "p_brine", "p_Na2SO4"],
         "filtrate was washed with water and brine; ethyl acetate extracts were dried and concentrated")
    step(r, "purification", "purification", ["p_silica", "p_EtOAc_eluent", "p_petroleum_ether"],
         "silica-gel column chromatography (100-200 mesh) using 15% ethyl acetate-petroleum ether (v/v)")
    r["procedure"]["general_procedure_reference"] = f(
        "General Procedure to synthesize 10a,b", gen10["reaction_id"], ev(13),
        "inherited_general_condition"
    )
    r["procedure"]["original_text"] = f(
        "2e-f (0.05 mmol), RuCl3 (0.5 mg, 0.0025 mmol), and NaIO4 (64.2 mg, 0.3 mmol) in CH3CN/CCl4/H2O (2:1:1), sealed tube, rt, 1 h; EtOAc dilution, thiosulfate quench, Celite filtration, washing, drying, concentration, and chromatography.",
        "Page-13 detailed procedure inherited by the product-specific experiment", ev(13),
        "inherited_general_condition"
    )
    product_amount(r, f"amount_{product}", mid(product), f"{mass} mg", mass)
    r["review"]["unresolved"].extend([
        "Individual volumes of the 2:1:1 solvent system are unreported.",
        "Workup-material amounts are unreported except for the 10 mL ethyl acetate dilution."
    ])
    bundle["reactions"].append(r)


# Existing Scheme 5 N-debenzylations: exact scale/solvent volume and isolated masses.
for substrate, product, mass in (("10a", "11a", 4.6), ("10b", "11b", 4.8)):
    r = copy.deepcopy(old_reactions[f"{DOC}:scheme5_{substrate}_to_{product}"])
    append_evidence(r, 13)
    r["location"] = f("PDF pages 6 and 13", {"pdf_pages": [6, 13]}, ev(6) + ev(13))
    r["general_scheme_ids"] = list(dict.fromkeys(
        r["general_scheme_ids"] + [gen11["reaction_id"]]
    ))
    set_charge_by_molecule(r, substrate, "amount", "0.02 mmol", {"number": 0.02, "unit": "mmol"})
    set_charge_by_molecule(r, "AlCl3", "amount", "0.06 mmol", {"number": 0.06, "unit": "mmol"})
    set_charge_by_molecule(r, "toluene", "volume", "1 mL", {"number": 1, "unit": "mL"})
    r["reaction_scale"] = f(
        f"{substrate} (0.02 mmol)",
        {"reference_molecule_id": mid(substrate), "amount": {"number": 0.02, "unit": "mmol"}}, ev(13)
    )
    r["solvent_conditions"]["composition"] = f("dry toluene", ["toluene"], ev(13))
    reaction_pids = [p["participant_id"] for p in r["participants"]]
    r["steps"] = []
    step(r, "reaction", "reaction", reaction_pids, "stirred in dry toluene",
         temperature=("rt", "room temperature"), time=("5 h", {"number": 5, "unit": "h"}))
    r["procedure"]["general_procedure_reference"] = f(
        "General Procedure to synthesize 11a,b", gen11["reaction_id"], ev(13),
        "inherited_general_condition"
    )
    r["procedure"]["original_text"] = f(
        "1-Aroyl-beta-carboline 10 (0.02 mmol) and anhydrous AlCl3 (0.06 mmol) in dry toluene (1 mL) were stirred at rt for 5 h.",
        "Page-13 detailed procedure inherited by the product-specific experiment", ev(13),
        "inherited_general_condition"
    )
    product_amount(r, f"amount_{product}", mid(product), f"{mass} mg", mass)
    r["notes"].append(f(
        f"Product {product} was reported as a known useful product and isolated using 20% ethyl acetate-petroleum ether (v/v); the cited literature supports prior identity/use, not a separate experiment of this paper.",
        {"product": product, "isolation_eluent": "20% ethyl acetate-petroleum ether (v/v)"}, ev(13)
    ))
    r["review"]["unresolved"].append(
        "The AlCl3-step workup and full purification operation are not reported in the general paragraph; only the characterization block's isolation eluent is stated."
    )
    bundle["reactions"].append(r)


# Include the already-existing generic Scheme 3 procedure only as a validation dependency
# for the four updated scope records. It is neither new nor semantically covered in this batch.
bundle["reactions"].append(copy.deepcopy(
    old_reactions[f"{DOC}:general_synthesis_spiroindolenines_9"]
))
bundle["molecules"] = list(molecules.values())
bundle["source_context"] = [
    {
        "context_id": "page13_characterization_blocks",
        "document_id": DOC,
        "evidence_ids": ev(13),
        "type": "characterization_context",
        "text": "NMR, HRMS, melting-point, and mass-spectral peak lists for 9a/9b/9i/9g/10a/10b/11a/11b are characterization rather than separate synthesis experiments; isolated masses and yield bases were attached to their existing experiment identities."
    },
    {
        "context_id": "page13_associated_content",
        "document_id": DOC,
        "evidence_ids": ev(13),
        "type": "source_availability_context",
        "text": "The article advertises Supporting Information and CCDC deposits, but no SI file was explicitly provided or associated in source_mapping.json; no SI association was inferred."
    },
    {
        "context_id": "page14_back_matter",
        "document_id": DOC,
        "evidence_ids": ev(14),
        "type": "nonreaction_context",
        "text": "Author information, declarations, acknowledgments, and references only; no experiment of this paper appears on the page."
    },
    {
        "context_id": "page15_references",
        "document_id": DOC,
        "evidence_ids": ev(15),
        "type": "literature_context",
        "text": "Reference list and citation annotations only. Cited prior-work statements, including reference 25's control description, were not promoted to experiments of this paper."
    },
    {
        "context_id": "page16_reference_continuation",
        "document_id": DOC,
        "evidence_ids": ev(16),
        "type": "literature_context",
        "text": "Continuation/end of reference 28 only; no reaction arrow, table row, or experiment of this paper is present."
    },
]

page13_ids = [
    r["reaction_id"] for r in bundle["reactions"]
    if any(evidence_id.startswith(f"{DOC}:page13:") for evidence_id in r["evidence_ids"])
]
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 13,
            "state": "attempted",
            "regions": [
                {
                    "region": "upper half, continuation of Characterization of 9a,b/9i/9g",
                    "finding": "Product masses and isolated/combined-yield bases for 9a, 9b, 9i, and 9g were merged into the existing Scheme 3 experiment identities; NMR/HRMS peak lists remain characterization context.",
                    "evidence_ids": ev(13)
                },
                {
                    "region": "middle, General Procedure to synthesize 10a,b and 11a,b",
                    "finding": "Two generic procedures were created for the oxidative-cleavage and AlCl3-debenzylation stages; exact scale, reagent doses, solvent ratio/volume, ordered workup, and reported missing details were retained.",
                    "evidence_ids": ev(13)
                },
                {
                    "region": "middle-lower, characterization of 10a/10b/11a/11b",
                    "finding": "Product masses and isolated yields were merged into the four existing Scheme 5 experiments. Known-product citations 27/28 were treated as prior-identity references, not new experiments.",
                    "evidence_ids": ev(13)
                },
                {
                    "region": "lower right, Associated Content and Accession Codes",
                    "finding": "SI and crystallographic deposits are advertised, but no SI was supplied in the source mapping; no SI association or SI-derived fact was inferred.",
                    "evidence_ids": ev(13)
                }
            ],
            "record_ids": page13_ids,
            "unresolved": [
                "Individual volumes for CH3CN/CCl4/H2O are unreported; only the 2:1:1 ratio is given.",
                "The AlCl3-step workup is not described; the product blocks state only the isolation eluent.",
                "Complex product SMILES remain unencoded rather than guessed, and all machine records remain needs_review.",
                "Advertised Supporting Information was not explicitly provided or associated."
            ]
        },
        {
            "page": 14,
            "state": "attempted",
            "regions": [
                {
                    "region": "full page",
                    "finding": "Author information, competing-interest note, acknowledgments, and reference list; no reaction table, arrow, or experiment of this paper was found.",
                    "evidence_ids": ev(14)
                }
            ],
            "record_ids": [],
            "unresolved": []
        },
        {
            "page": 15,
            "state": "attempted",
            "regions": [
                {
                    "region": "both columns, references 7-28",
                    "finding": "Bibliographic references and cited prior-work annotations only. These descriptions were not treated as experiments in the current paper.",
                    "evidence_ids": ev(15)
                }
            ],
            "record_ids": [],
            "unresolved": [
                "External cited articles were not opened or reconstructed because this batch covers only the supplied source pages."
            ]
        },
        {
            "page": 16,
            "state": "attempted",
            "regions": [
                {
                    "region": "top-left reference continuation and otherwise blank page",
                    "finding": "The final line of reference 28 continues from page 15; no experimental content, reaction arrow, or table is present.",
                    "evidence_ids": ev(16)
                }
            ],
            "record_ids": [],
            "unresolved": []
        }
    ]
}

raw_inventory = {
    "document_id": DOC,
    "page_range": [13, 16],
    "semantic_accuracy_verified": False,
    "new_reaction_ids": [gen10["reaction_id"], gen11["reaction_id"]],
    "updated_reaction_ids": [
        f"{DOC}:scheme3_9a", f"{DOC}:scheme3_9b", f"{DOC}:scheme3_9i", f"{DOC}:scheme3_9g",
        f"{DOC}:scheme5_2e_to_10a", f"{DOC}:scheme5_2f_to_10b",
        f"{DOC}:scheme5_10a_to_11a", f"{DOC}:scheme5_10b_to_11b"
    ],
    "notes": [
        "All four supplied page images were visually inspected at original detail.",
        "Pages 14-16 contain back matter/references only and therefore create no reaction records.",
        "Spectral peak lists remain source_context, not synthesis experiments.",
        "No SI association was inferred and no complex SMILES was guessed."
    ]
}

dump(RUN / "raw_model/pages_0013_0016_extraction.json", raw_inventory)
dump(SUBMITTED, bundle)
dump(RUN / "coverage/pages_0013_0016.json", coverage)
print(json.dumps({
    "submitted_reactions": len(bundle["reactions"]),
    "new_reactions": len(raw_inventory["new_reaction_ids"]),
    "updated_reactions": len(raw_inventory["updated_reaction_ids"]),
    "bundle_molecules": len(bundle["molecules"]),
    "coverage_pages": [p["page"] for p in coverage["pages"]]
}, ensure_ascii=False, indent=2))
