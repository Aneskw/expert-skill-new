#!/usr/bin/env python3
"""Extend the existing v3.1 bundle with source facts from PDF pages 7-11."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_484c4de81c12e684"
RUN = ROOT / "runs/current"
BUNDLE_PATH = RUN / "submitted.json"
BUNDLE = json.loads(BUNDLE_PATH.read_text())
RAW = json.loads((RUN / "raw_model/page7_11_extraction.json").read_text())
TEMPLATE = json.loads((RUN / "raw_model/page7_11_template.json").read_text())

# The previous page-bounded bundle contains all page text but only images 1-6.
# Bring in the page 7-11 image evidence objects from the freshly initialized
# template before any new facts cite them.
bundle_document = next(item for item in BUNDLE["documents"] if item["document_id"] == DOC)
template_document = next(item for item in TEMPLATE["documents"] if item["document_id"] == DOC)
known_evidence = {item["evidence_id"] for item in bundle_document["evidence"]}
for evidence_item in template_document["evidence"]:
    if evidence_item["evidence_id"] not in known_evidence:
        bundle_document["evidence"].append(evidence_item)
        known_evidence.add(evidence_item["evidence_id"])


def ev(*pages: int) -> list[str]:
    return [item for page in pages for item in (f"{DOC}:page{page}:text", f"{DOC}:page{page}:image")]


def source(raw, normalized=None, evidence=()):
    return st.fact(raw, raw if normalized is None else normalized, evidence)


def inferred(normalized, evidence, derivation, raw=None):
    return st.fact(raw, normalized, evidence, method="inference", derivation=derivation)


def not_applicable(evidence, derivation):
    return st.fact(
        evidence=evidence,
        method="not_applicable",
        status="not_applicable",
        derivation=derivation,
    )


def append_unique(values: list[str], additions: list[str]) -> None:
    for value in additions:
        if value not in values:
            values.append(value)


def extend_fact_evidence(fact: dict, evidence: list[str]) -> None:
    append_unique(fact["evidence_ids"], evidence)


def molecule(mid: str) -> dict:
    return next(item for item in BUNDLE["molecules"] if item["molecule_id"] == mid)


def reaction(rid: str) -> dict:
    return next(item for item in BUNDLE["reactions"] if item["reaction_id"] == rid)


def participant(r: dict, pid: str) -> dict:
    return next(item for item in r["participants"] if item["participant_id"] == pid)


def add_participant_once(r: dict, pid: str, mid: str, role: str, evidence: list[str], stage="reaction") -> dict:
    found = [item for item in r["participants"] if item["participant_id"] == pid]
    if found:
        return found[0]
    item = st.new_participant(pid, mid, role, evidence, stage)
    r["participants"].append(item)
    return item


def set_charge(p: dict, key: str, raw, normalized, evidence: list[str], method="original_extraction", derivation=None):
    p["charges"][key] = st.fact(raw, normalized, evidence, method=method, derivation=derivation)


def set_reference(p: dict, raw, mid: str, evidence: list[str], method="original_extraction", derivation=None):
    p["equivalence_reference"] = st.fact(raw, mid, evidence, method=method, derivation=derivation)


def make_step(step_id, stage, participant_ids, operation, evidence, temperature=None, time=None,
              operation_method="original_extraction", condition_method="original_extraction",
              operation_derivation=None, condition_derivation=None):
    item = st.new_step(step_id, evidence, stage)
    item["participant_ids"] = participant_ids
    item["operation"] = st.fact(
        operation,
        operation,
        evidence,
        method=operation_method,
        derivation=operation_derivation,
    )
    if temperature is not None:
        item["temperature"] = st.fact(
            temperature[0],
            temperature[1],
            evidence,
            method=condition_method,
            derivation=condition_derivation,
        )
    if time is not None:
        item["time"] = st.fact(
            time[0],
            time[1],
            evidence,
            method=condition_method,
            derivation=condition_derivation,
        )
    return item


PRODUCT_NAMES = {
    "3c": "Benzyl (E)-2-(7-methoxy-5-oxo-4-phenyl-1,10b-dihydro-2H,5H-pyrano[3,4-c]chromen-2-ylidene)acetate",
    "3d": "Benzyl (E)-2-(8-methoxy-5-oxo-4-phenyl-1,10b-dihydro-2H,5H-pyrano[3,4-c]chromen-2-ylidene)acetate",
    "3e": "Benzyl (E)-2-(9-methyl-5-oxo-4-phenyl-1,10b-dihydro-2H,5H-pyrano[3,4-c]chromen-2-ylidene)acetate",
    "3f": "Benzyl (E)-2-(7,9-di-tert-butyl-5-oxo-4-phenyl-1,10b-dihydro-2H,5H-pyrano[3,4-c]chromen-2-ylidene)acetate",
    "3g": "Benzyl (E)-2-(9-fluoro-5-oxo-4-phenyl-1,10b-dihydro-2H,5H-pyrano[3,4-c]chromen-2-ylidene)acetate",
    "3h": "Benzyl (E)-2-(9-chloro-5-oxo-4-phenyl-1,10b-dihydro-2H,5H-pyrano[3,4-c]chromen-2-ylidene)acetate",
    "3i": "Benzyl (E)-2-(9-bromo-5-oxo-4-phenyl-1,10b-dihydro-2H,5H-pyrano[3,4-c]chromen-2-ylidene)acetate",
    "3j": "Benzyl (E)-2-(8-bromo-5-oxo-4-phenyl-1,10b-dihydro-2H,5H-pyrano[3,4-c]chromen-2-ylidene)acetate",
    "3k": "Benzyl (E)-2-(5-oxo-4-(p-tolyl)-1,10b-dihydro-2H,5H-pyrano[3,4-c]chromen-2-ylidene)acetate",
    "3l": "Benzyl (E)-2-(5-oxo-4-(m-tolyl)-1,10b-dihydro-2H,5H-pyrano[3,4-c]chromen-2-ylidene)acetate",
    "3m": "Benzyl (E)-2-(4-(4-chlorophenyl)-5-oxo-1,10b-dihydro-2H,5H-pyrano[3,4-c]chromen-2-ylidene)acetate",
    "3n": "Benzyl (E)-2-(4-(3-chlorophenyl)-5-oxo-1,10b-dihydro-2H,5H-pyrano[3,4-c]chromen-2-ylidene)acetate",
    "3o": "Benzyl (E)-2-(4-([1,1'-biphenyl]-4-yl)-5-oxo-1,10b-dihydro-2H,5H-pyrano[3,4-c]chromen-2-ylidene)acetate",
    "3p": "Benzyl (E)-2-(5-oxo-4-(thiophen-3-yl)-1,10b-dihydro-2H,5H-pyrano[3,4-c]chromen-2-ylidene)acetate",
    "5": "Ethyl (E)-2-(4-(4-chlorophenyl)-5-oxo-1,10b-dihydro-2H,5H-pyrano[3,4-c]chromen-2-ylidene)acetate",
}

NAME_PAGES = {
    "3c": (7,), "3d": (7,), "3e": (7,), "3f": (7,), "3g": (7,),
    "3h": (8,), "3i": (8,), "3j": (8,), "3k": (8,), "3l": (8,),
    "3m": (9,), "3n": (9,), "3o": (9,), "3p": (9,), "5": (10,),
}

for label, name in PRODUCT_NAMES.items():
    mid = f"m_{label}"
    item = molecule(mid)
    evidence = ev(*NAME_PAGES[label])
    item["names"]["original"] = source(name, name, evidence)
    item["names"]["english"] = source(name, name, evidence)
    extend_fact_evidence(item["names"]["abbreviations"], evidence)
    append_unique(item["occurrences"][0]["evidence_ids"], evidence)
    extend_fact_evidence(item["occurrences"][0]["compound_label"], evidence)
    # The compound heading itself explicitly carries (E); absolute configuration remains absent.
    extend_fact_evidence(item["stereochemistry"]["double_bonds"], evidence)
    extend_fact_evidence(item["stereochemistry"]["assignment"], evidence)
    extend_fact_evidence(item["stereochemistry"]["geometry_status"], evidence)


def make_free_acid() -> dict:
    evidence = ev(9, 10)
    item = st.new_molecule("m_free_acid_3m", DOC, evidence)
    item["identity_status"] = "source_identified"
    item["occurrences"][0]["compound_label"] = source("free acid intermediate", "free acid intermediate", evidence)
    item["names"]["original"] = source("free acid intermediate", "free acid intermediate", evidence)
    item["names"]["english"] = inferred(
        "free-acid intermediate obtained from 3m by Pd/C-H2 treatment",
        evidence,
        "Normalized from the section 3.3 prose; no systematic name or isolated characterization is reported.",
    )
    item["names"]["abbreviations"] = st.fact(evidence=evidence)
    item["names"]["aliases"] = st.fact(evidence=evidence)
    item["cas_numbers"] = st.fact(evidence=evidence)
    item["representation_type"] = "unresolved"
    item["smiles"] = st.fact(evidence=evidence)
    item["generic_template"] = None
    item["stereochemistry"]["chirality_status"] = inferred(
        "unresolved", evidence,
        "The intermediate is not structurally depicted or separately assigned; stereochemistry is not inferred solely from its precursor."
    )
    item["stereochemistry"]["geometry_status"] = inferred(
        "unresolved", evidence,
        "The source calls it only a free-acid intermediate and gives no separate E/Z assignment."
    )
    item["stereochemistry"]["representation_status"] = inferred(
        "unresolved", evidence,
        "No source-verified structure representation is available for this intermediate."
    )
    item["identity_note"] = source(
        "free acid intermediate",
        "Intermediate isolated by the first silica-gel column; amount, yield, and full identity are unreported.",
        evidence,
    )
    item["external_properties"] = []
    return item


def make_celite() -> dict:
    evidence = ev(9)
    item = st.new_molecule("m_celite", DOC, evidence)
    item["identity_status"] = "source_identified"
    item["occurrences"][0]["compound_label"] = source("Celite", "Celite", evidence)
    item["names"]["original"] = source("Celite", "Celite", evidence)
    item["names"]["english"] = source("Celite", "Celite filtration aid", evidence)
    item["names"]["abbreviations"] = st.fact(evidence=evidence)
    item["names"]["aliases"] = st.fact(evidence=evidence)
    item["cas_numbers"] = st.fact(evidence=evidence)
    item["representation_type"] = "unresolved"
    item["smiles"] = st.fact(evidence=evidence)
    item["generic_template"] = None
    item["stereochemistry"]["mixture"] = source(
        "Celite", {"material": "Celite", "composition": "not specified"}, evidence
    )
    item["stereochemistry"]["chirality_status"] = inferred(
        "not_applicable", evidence,
        "Celite is used as a filtration material rather than a discrete stereochemical molecular entity."
    )
    item["stereochemistry"]["geometry_status"] = inferred(
        "not_applicable", evidence,
        "Molecular E/Z or ring geometry is not applicable to the unspecified filtration material."
    )
    item["stereochemistry"]["representation_status"] = inferred(
        "not_applicable", evidence,
        "A single-molecule stereochemical representation is not applicable to this filtration material."
    )
    item["identity_note"] = source(
        "a pad of Celite", "Filtration aid; composition and amount are not reported.", evidence
    )
    item["external_properties"] = []
    return item


for new_item in (make_free_acid(), make_celite()):
    existing_indexes = [
        index for index, item in enumerate(BUNDLE["molecules"])
        if item["molecule_id"] == new_item["molecule_id"]
    ]
    if existing_indexes:
        BUNDLE["molecules"][existing_indexes[0]] = new_item
    else:
        BUNDLE["molecules"].append(new_item)


# Treat 10 wt% as the source specification of Pd/C, not as a catalyst charge.
pdc = molecule("m_pdc")
pdc["stereochemistry"]["mixture"] = source(
    "Pd/C (wt. 10%)",
    {"material": "palladium on carbon", "palladium_mass_fraction": {"number": 10, "unit": "%"}},
    ev(9),
)
pdc["stereochemistry"]["chirality_status"] = inferred(
    "not_applicable", ev(9),
    "A supported-metal catalyst material is not represented as a discrete stereochemical molecule."
)
pdc["stereochemistry"]["geometry_status"] = inferred(
    "not_applicable", ev(9),
    "Molecular E/Z or ring geometry is not applicable to the supported catalyst material."
)
pdc["stereochemistry"]["representation_status"] = inferred(
    "not_applicable", ev(9),
    "A single-molecule stereochemical representation is not applicable to Pd/C."
)
pdc["identity_note"] = source(
    "Pd/C (wt. 10%)",
    "10 wt% palladium-on-carbon catalyst grade; the amount charged is not reported.",
    ev(9),
)


# Attach compound-specific HPLC methods to the existing selectivity outcomes.
for label, method in RAW["hplc_methods"].items():
    if label == "5":
        target_reactions = [
            reaction(f"{DOC}:figure3_direct_to_5"),
            reaction(f"{DOC}:figure3_three_step_to_5"),
        ]
    else:
        target_reactions = [reaction(f"{DOC}:scope_{label}")]
    evidence = ev(*method["pages"])
    raw_method = (
        f"{method['column']}; {method['mobile_phase']}; flow rate = "
        f"{method['flow_rate']['number']} {method['flow_rate']['unit']}; retention time: "
        f"{method['retention_times_min']['major']} min (major) and "
        f"{method['retention_times_min']['minor']} min (minor)"
    )
    normalized_method = {
        "method": "chiral HPLC",
        "column": method["column"],
        "mobile_phase": method["mobile_phase"],
        "flow_rate": method["flow_rate"],
        "retention_times": {
            "major": {"number": method["retention_times_min"]["major"], "unit": "min"},
            "minor": {"number": method["retention_times_min"]["minor"], "unit": "min"},
        },
        "configuration_assignment": "major and minor absolute configurations unassigned",
    }
    for target_reaction in target_reactions:
        ee = next(item for item in target_reaction["outcomes"] if item["outcome_id"] == "ee")
        ee["measurement_basis"] = source(raw_method, normalized_method, evidence)
        ee["scope"] = source(
            f"retention time: {method['retention_times_min']['major']} min (major) and {method['retention_times_min']['minor']} min (minor)",
            {
                "metric": "ee",
                "major_enantiomer": "unassigned",
                "minor_enantiomer": "unassigned",
                "reported_peak_order": ["major", "minor"],
            },
            evidence,
        )
        append_unique(target_reaction["evidence_ids"], evidence)


# Replace the abbreviated lower Figure 3 route with the complete section 3.3 procedure.
route = reaction(f"{DOC}:figure3_three_step_to_5")
append_unique(route["evidence_ids"], ev(9, 10))
route["location"] = source(
    "Figure 3 lower route and section 3.3 (pages 5, 9-10)",
    {"pages": [4, 5, 9, 10], "regions": ["Figure 3 lower route", "section 3.3"]},
    ev(4, 5, 9, 10),
)

p_1m = participant(route, "p_1m")
set_charge(p_1m, "amount", "1m (1 mmol)", {"number": 1, "unit": "mmol"}, ev(9))

p_2b = participant(route, "p_2b")
set_charge(p_2b, "amount", "2b (1.25 mmol)", {"number": 1.25, "unit": "mmol"}, ev(9))
alternative_standard = st.fact(
    "Figure 3: standard conditions; Figure 2 defines 2b as 1.5 equiv.",
    {"number": 1.5, "unit": "equiv"},
    ev(5),
    method="inherited_general_condition",
    derivation="Figure 3 invokes the standard conditions defined in Figure 2.",
)
alternative_charges = st.fact(
    "section 3.3: 1m (1 mmol) and 2b (1.25 mmol)",
    {"number": 1.25, "unit": "equiv"},
    ev(9),
    method="inference",
    derivation="Calculated as 1.25 mmol 2b divided by 1.00 mmol 1m.",
)
p_2b["charges"]["equivalents"] = st.fact(
    [alternative_standard["raw_value"], alternative_charges["raw_value"]],
    None,
    ev(5, 9),
    method="inference",
    status="conflict",
    derivation="The explicit section 3.3 charges imply 1.25 equiv, whereas the invoked standard conditions specify 1.5 equiv.",
    alternatives=[alternative_standard, alternative_charges],
)
set_reference(p_2b, "relative to 1m", "m_1m", ev(5, 9), method="inference",
              derivation="Both the standard-condition equivalent and the charge-derived equivalent use 1m as denominator.")
route["reaction_scale"] = source(
    "1m (1 mmol)",
    {"reference_molecule_id": "m_1m", "amount": {"number": 1, "unit": "mmol"}},
    ev(9),
)

p_dcc = participant(route, "p_dcc")
set_charge(p_dcc, "equivalents", "DCC (2.0 equiv.)", {"number": 2.0, "unit": "equiv"}, ev(10))
set_reference(
    p_dcc, None, "m_free_acid_3m", ev(10), method="inference",
    derivation="The source prints 2.0 equiv after isolating the free-acid intermediate but does not explicitly restate the denominator."
)
p_dmap = participant(route, "p_dmap")
set_charge(p_dmap, "equivalents", "DMAP (2.0 equiv.)", {"number": 2.0, "unit": "equiv"}, ev(10))
set_reference(
    p_dmap, None, "m_free_acid_3m", ev(10), method="inference",
    derivation="The source prints 2.0 equiv after isolating the free-acid intermediate but does not explicitly restate the denominator."
)
p_etoh = participant(route, "p_etoh")
set_charge(p_etoh, "volume", "EtOH (2 mL)", {"number": 2, "unit": "mL"}, ev(10))
p_ch2cl2 = participant(route, "p_ch2cl2")
set_charge(p_ch2cl2, "volume", "CH2Cl2 (5 mL)", {"number": 5, "unit": "mL"}, ev(10))

p_free_product = add_participant_once(route, "p_free_acid_product", "m_free_acid_3m", "product", ev(9, 10))
p_free_product["product_role"] = source("free acid intermediate", "intermediate", ev(9, 10))
add_participant_once(route, "p_free_acid_reactant", "m_free_acid_3m", "reactant", ev(10))
add_participant_once(route, "p_celite", "m_celite", "workup_agent", ev(9), "workup")
add_participant_once(route, "p_etoac_elution", "m_etoac", "workup_agent", ev(9), "workup")
add_participant_once(route, "p_pe_first_column", "m_pe", "purification_agent", ev(10), "purification")
add_participant_once(route, "p_etoac_first_column", "m_etoac", "purification_agent", ev(10), "purification")
add_participant_once(route, "p_pe_final_column", "m_pe", "purification_agent", ev(10), "purification")
add_participant_once(route, "p_etoac_final_column", "m_etoac", "purification_agent", ev(10), "purification")

route["solvent_conditions"]["composition"] = st.fact(
    ["standard conditions", "dissolved in CH2Cl2 (5 mL)", "EtOH (2 mL)"],
    {
        "annulation": ["ethyl acetate"],
        "esterification": ["dichloromethane", "ethanol"],
        "stage_specific_volumes": {"CH2Cl2": {"number": 5, "unit": "mL"}, "EtOH": {"number": 2, "unit": "mL"}},
    },
    ev(5, 9, 10),
    method="inference",
    status="pending_review",
    derivation="Stage-specific solvent composition is assembled from the invoked standard conditions and the explicit esterification additions.",
)
route["solvent_conditions"]["total_volume"] = st.fact(
    None,
    None,
    ev(10),
    method="inference",
    status="unreported",
    derivation="The source reports stage-specific additions, not a measured total solution volume.",
)

route["steps"] = [
    make_step(
        "annulation_pre_stir", "reaction", ["p_1m", "p_4l", "p_etoac"],
        "standard-condition pre-stir before charging 2b", ev(5, 6, 9),
        ("0 °C", {"number": 0, "unit": "°C"}), ("15 min", {"number": 15, "unit": "min"}),
        "inherited_general_condition", "inherited_general_condition",
        "Section 3.3 says the reaction mixture was obtained under standard conditions; section 3.2 defines the pre-stir.",
        "Inherited from the standard-condition general procedure invoked by section 3.3.",
    ),
    make_step(
        "annulation", "reaction", ["p_1m", "p_2b", "p_4l", "p_etoac", "p_3m"],
        "obtain the reaction mixture under the standard condition using 1m and 2b", ev(5, 9),
        ("0 °C", {"number": 0, "unit": "°C"}), None,
        "original_extraction", "inherited_general_condition", None,
        "The 0 °C condition is inherited from the standard conditions invoked in section 3.3; annulation time is not reported.",
    ),
    make_step(
        "hydrogenolysis", "reaction", ["p_3m", "p_pdc", "p_h2", "p_free_acid_product"],
        "add Pd/C (wt. 10%) and stir under H2 atmosphere (H2 balloon)", ev(9),
        ("room temperature", "room temperature"), ("48 h", {"number": 48, "unit": "h"}),
    ),
    make_step(
        "celite_filtration", "workup", ["p_celite", "p_etoac_elution"],
        "filter through a pad of Celite and elute with EtOAc", ev(9),
    ),
    make_step(
        "first_concentration", "workup", [],
        "concentrate the filtrate under reduced pressure", ev(9, 10),
    ),
    make_step(
        "free_acid_purification", "purification", ["p_pe_first_column", "p_etoac_first_column", "p_free_acid_product"],
        "silica gel column chromatography (PE:EtOAc = 5:1-2:1) to afford the free acid intermediate", ev(10),
    ),
    make_step(
        "esterification_dissolution", "reaction", ["p_free_acid_reactant", "p_ch2cl2"],
        "dissolve the free acid intermediate in CH2Cl2 (5 mL)", ev(10),
    ),
    make_step(
        "esterification_cold_addition", "reaction", ["p_free_acid_reactant", "p_dcc", "p_dmap", "p_etoh", "p_ch2cl2"],
        "cool the solution to 0 °C before adding DCC, DMAP, and EtOH", ev(10),
        ("0 °C", {"number": 0, "unit": "°C"}), None,
    ),
    make_step(
        "esterification_room_temperature_hold", "reaction", ["p_free_acid_reactant", "p_dcc", "p_dmap", "p_etoh", "p_ch2cl2", "p_5"],
        "move the reaction mixture to room temperature and stir; monitor completion by TLC", ev(10),
        ("rt", "rt"), ("overnight", "overnight"),
    ),
    make_step(
        "final_evaporation", "workup", [],
        "remove the solvent by rotary evaporation", ev(10),
    ),
    make_step(
        "final_purification", "purification", ["p_pe_final_column", "p_etoac_final_column", "p_5"],
        "column chromatography using PE/EtOAc (15:1-9:1) as eluent to afford product 5", ev(10),
    ),
]
route["steps"][2]["atmosphere"] = source("H2 atmosphere (H2 balloon)", "hydrogen balloon", ev(9))

route["procedure"]["original_text"] = source(
    "To the reaction mixture obtained under standard condition using 1m (1 mmol) and 2b (1.25 mmol) was added Pd/C (wt. 10%), then the mixture was stirred under H2 atmosphere (H2 balloon) at room temperature for 48 h. The resulting reaction mixture was filtered through a pad of Celite and eluted with EtOAc. The filtration was concentrated under reduced pressure and the residue was purified by silica gel column chromatography (PE:EtOAc = 5:1-2:1) to afford the free acid intermediate, which was then dissolved in CH2Cl2 (5 mL). This solution was cooled to 0 °C before DCC (2.0 equiv.), DMAP (2.0 equiv.) and EtOH (2 mL) were added. After that, the reaction mixture was moved to rt and stirred overnight. After completion of the reaction by TLC analysis, the solvent was removed by rotary evaporation and the residue was directly subjected to column chromatography using PE/EtOAc (15:1-9:1) as the eluent to afford product 5.",
    "section 3.3 procedure spanning pages 9-10",
    ev(9, 10),
)
route["procedure"]["general_procedure_reference"] = source(
    "reaction mixture obtained under standard condition",
    f"{DOC}:scheme_current_general",
    ev(5, 6, 9),
)
route["procedure"]["modifications"] = source(
    "1m (1 mmol), 2b (1.25 mmol), followed by Pd/C-H2 hydrogenolysis and DCC/DMAP/EtOH esterification",
    {
        "annulation_scale": {"1m": {"number": 1, "unit": "mmol"}, "2b": {"number": 1.25, "unit": "mmol"}},
        "post_annulation_sequence": ["Pd/C-H2 hydrogenolysis", "DCC/DMAP/EtOH esterification"],
    },
    ev(9, 10),
)
route["procedure"]["addition_order"] = source(
    "standard-condition reaction mixture; add Pd/C; H2 balloon; filter/concentrate/column; dissolve free acid in CH2Cl2; cool to 0 °C; add DCC, DMAP, EtOH; move to rt; evaporate; column",
    [
        "standard-condition annulation",
        "Pd/C addition and H2 hold",
        "Celite filtration and first purification",
        "CH2Cl2 dissolution",
        "cool to 0 °C",
        "DCC, DMAP, EtOH addition",
        "room-temperature overnight hold",
        "evaporation and final purification",
    ],
    ev(9, 10),
)
route["procedure"]["pretreatment"] = source(
    "free acid intermediate solution was cooled to 0 °C before DCC, DMAP and EtOH were added",
    "cool free-acid solution in CH2Cl2 to 0 °C before reagent addition",
    ev(10),
)

# Add the observed ester C-O cleavage to the already partial multistep change inventory.
changes = route["classification"]["bond_changes"]
if not any((item["normalized_value"] or {}).get("change_id") == "c4" for item in changes):
    changes.append(st.fact(
        None,
        {
            "change_id": "c4",
            "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}],
            "before_order": 1,
            "after_order": 0,
            "change_type": "broken",
            "source_basis": "benzyl ester hydrogenolysis to the reported free acid intermediate",
        },
        ev(9, 10),
        method="structure_derivation",
        derivation="The section 3.3 sequence reports Pd/C-H2 treatment of benzyl ester 3m followed by isolation of a free acid; atom mapping remains unresolved.",
    ))
append_unique(route["classification"]["labels"]["evidence_ids"], ev(9, 10))
append_unique(route["classification"]["coverage"]["evidence_ids"], ev(9, 10))
route["classification"]["limitations"] = [
    "No complete atom-mapped reaction SMILES was asserted.",
    "The multistep record retains annulation C-C/C-O formation, benzyl-ester C-O cleavage, and ethyl-ester C-O formation; proton transfers and coproducts remain unresolved.",
    "The free-acid intermediate is source-mentioned but not fully identified or quantified.",
    "Classification is a machine derivation from the supplied drawings and procedure, not source verification.",
]

route["notes"] = [
    st.fact(
        [
            "Figure 3 lower route: standard conditions",
            "section 3.3: 1m (1 mmol) and 2b (1.25 mmol)",
        ],
        None,
        ev(5, 9),
        method="inference",
        status="conflict",
        derivation="The invoked standard conditions use 1.5 equiv 2b, while the explicit section 3.3 charges imply 1.25 equiv; neither is selected as an authorial correction.",
        alternatives=[alternative_standard, alternative_charges],
    ),
    inferred(
        "The 10 wt% value describes the Pd/C catalyst grade, not the unreported amount charged.",
        ev(9),
        "The source provides a mass fraction for Pd/C but no charge amount.",
        raw="Pd/C (wt. 10%)",
    ),
]
route["review"] = {
    "status": "needs_review",
    "method": "model",
    "reviewer": None,
    "checked_fields": [],
    "unresolved": [
        "Machine-extracted record; chemical identity and full atom mapping require expert review.",
        "Figure 3 invokes standard conditions with 2b at 1.5 equiv, but section 3.3 explicitly gives 1m 1 mmol and 2b 1.25 mmol (1.25 equiv by calculation); the conflict is retained.",
        "The amount of 10 wt% Pd/C charged is not reported.",
        "The free-acid intermediate amount, isolated yield, full identity, and stereochemical assignment are not reported.",
        "The standard-condition annulation time is not reported for this experiment.",
        "The 48% value remains a three-step total yield from 1m to 5, not a last-step yield.",
        "Absolute configuration remains unassigned; the HPLC peaks are only identified as major/minor.",
    ],
}
route["completeness"] = "partial"


# Product 5 characterization supports both routes but is not a separate synthesis.
for rid in (f"{DOC}:figure3_direct_to_5", f"{DOC}:figure3_three_step_to_5"):
    current = reaction(rid)
    append_unique(current["evidence_ids"], ev(10))


contexts = {item["context_id"]: item for item in BUNDLE["source_context"]}
new_contexts = [
    {
        "context_id": f"{DOC}:characterization_3c_3p",
        "kind": "characterization",
        "evidence_ids": ev(7, 8, 9),
        "note": "Pages 7-9 contain compound headings, melting points, NMR, optical rotations, compound-specific chiral-HPLC conditions, and HRMS for 3c-3p. These observations support identity and selectivity measurement context but are not additional synthesis experiments.",
    },
    {
        "context_id": f"{DOC}:characterization_5",
        "kind": "characterization",
        "evidence_ids": ev(10),
        "note": "Page 10 reports the name, NMR, optical rotation, and compound-specific chiral-HPLC method for 5; no separate yield is reported in this characterization block.",
    },
    {
        "context_id": f"{DOC}:supplementary_statement_page10",
        "kind": "source_limitation",
        "evidence_ids": ev(10),
        "note": "The paper says spectra and HPLC data are available online, but no supplementary file was explicitly provided or associated in this workspace.",
    },
    {
        "context_id": f"{DOC}:references_pages10_11",
        "kind": "bibliography",
        "evidence_ids": ev(10, 11),
        "note": "The reference list contains cited literature titles only; it does not provide additional current-paper experimental procedures or results.",
    },
]
for item in new_contexts:
    if item["context_id"] not in contexts:
        BUNDLE["source_context"].append(item)


BUNDLE_PATH.write_text(json.dumps(BUNDLE, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({
    "molecules": len(BUNDLE["molecules"]),
    "reactions": len(BUNDLE["reactions"]),
    "new_molecule_ids": ["m_free_acid_3m", "m_celite"],
    "new_reaction_ids": [],
    "updated_reaction_ids": [
        f"{DOC}:scope_3c", f"{DOC}:scope_3d", f"{DOC}:scope_3e", f"{DOC}:scope_3f",
        f"{DOC}:scope_3g", f"{DOC}:scope_3h", f"{DOC}:scope_3i", f"{DOC}:scope_3j",
        f"{DOC}:scope_3k", f"{DOC}:scope_3l", f"{DOC}:scope_3m", f"{DOC}:scope_3n",
        f"{DOC}:scope_3o", f"{DOC}:scope_3p", f"{DOC}:figure3_direct_to_5",
        f"{DOC}:figure3_three_step_to_5",
    ],
}, ensure_ascii=False))
