#!/usr/bin/env python3
"""Build the source-derived v3.1 delta for original PDF pages 7--9.

The script is intentionally incremental: it starts from the current v3.1 records,
adds the single scale-up experiment on page 7, enriches the six existing Table 4
records whose characterization continues on page 7, and records pages 8--9 as
bibliography context. It does not certify chemical correctness.
"""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st


RUN = ROOT / "runs/current"
DOC = "doc_64f4a881f331038c"
E = {p: f"{DOC}:page{p}:text" for p in range(1, 10)}
I = {p: f"{DOC}:page{p}:image" for p in range(1, 10)}
SCALE_ID = f"{DOC}:scale_up_5_mmol_1a_2a_to_4a"


def dump(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def uniq(values):
    return list(dict.fromkeys(values))


_UNSET = object()


def f(raw=None, normalized=_UNSET, evidence=(), method="original_extraction",
      status=None, derivation=None, alternatives=()):
    if normalized is _UNSET:
        normalized = raw
    return st.fact(raw, normalized, tuple(dict.fromkeys(evidence)), method=method,
                   status=status, derivation=derivation,
                   alternatives=alternatives)


def missing(evidence=(), status="unreported"):
    return f(evidence=evidence, status=status)


def not_applicable(evidence=()):
    return f(evidence=evidence, method="not_applicable", status="not_applicable")


def derived(raw, normalized, evidence, derivation):
    return f(raw, normalized, evidence, method="structure_derivation",
             derivation=derivation)


def q(number, unit, raw_text=None, **extra):
    value = {"number": number, "unit": unit}
    if raw_text is not None:
        value["raw_text"] = raw_text
    value.update(extra)
    return value


def participant(pid, molecule_id, role, evidence, stage="reaction", mass=None,
                amount=None, volume=None, product_role=None):
    p = st.new_participant(pid, molecule_id, role, evidence, stage=stage)
    if mass is not None:
        p["charges"]["mass"] = f(mass.get("raw_text", str(mass)), mass, evidence)
    if amount is not None:
        p["charges"]["amount"] = f(amount.get("raw_text", str(amount)), amount, evidence)
    if volume is not None:
        raw = volume.get("raw_text", str(volume)) if isinstance(volume, dict) else str(volume)
        p["charges"]["volume"] = f(raw, volume, evidence)
    if product_role is not None:
        p["product_role"] = f(product_role, product_role, evidence)
    return p


def step(step_id, stage, participant_ids, evidence, operation, temperature=None,
         time=None, pressure=None):
    s = st.new_step(step_id, evidence, stage=stage)
    s["participant_ids"] = participant_ids
    s["operation"] = f(operation, operation, evidence)
    s["temperature"] = (f(str(temperature), temperature, evidence)
                        if temperature is not None else missing(evidence))
    s["time"] = (f(str(time), time, evidence)
                 if time is not None else missing(evidence))
    s["pressure"] = (f(pressure, pressure, evidence)
                     if pressure is not None else missing(evidence))
    return s


def outcome(outcome_id, kind, targets, evidence, value, product_role=None,
            yield_type=None, basis=None, scope=None):
    o = st.new_outcome(outcome_id, kind, targets, evidence)
    raw = value.get("raw_text", str(value)) if isinstance(value, dict) else str(value)
    o["value"] = f(raw, value, evidence)
    o["product_role"] = (f(product_role, product_role, evidence)
                         if product_role else not_applicable(evidence))
    o["yield_type"] = (f(yield_type, yield_type, evidence)
                       if yield_type else missing(evidence))
    o["measurement_basis"] = (f(basis, basis, evidence)
                              if basis else missing(evidence))
    o["scope"] = f(scope, scope, evidence) if scope is not None else missing(evidence)
    return o


def add_page7_image_where_text_is_cited(value):
    """Add page-7 image lineage only to facts already citing page-7 text."""
    if isinstance(value, dict):
        if "evidence_ids" in value and isinstance(value["evidence_ids"], list):
            if E[7] in value["evidence_ids"]:
                value["evidence_ids"] = uniq(value["evidence_ids"] + [I[7]])
        for nested in value.values():
            add_page7_image_where_text_is_cited(nested)
    elif isinstance(value, list):
        for nested in value:
            add_page7_image_where_text_is_cited(nested)


init = json.loads((RUN / "raw_model/pages_0007_0009_init.json").read_text())
dataset = json.loads((RUN / "dataset.json").read_text())
bundle = copy.deepcopy(dataset["standard"])

# The init output contains the authoritative prepared document evidence, including
# rendered page images 7--9. Keep the established document identity/title.
bundle["documents"] = init["documents"]


# Add page-7 occurrence lineage to reused scale-up molecules without pretending
# that the prose alone re-verifies their complete structures/stereochemistry.
for molecule in bundle["molecules"]:
    mid = molecule["molecule_id"]
    if mid in {"m_1a", "m_2a", "m_4a"}:
        occurrence = molecule["occurrences"][0]
        occurrence["evidence_ids"] = uniq(occurrence["evidence_ids"] + [E[7], I[7]])
        occurrence["compound_label"]["evidence_ids"] = uniq(
            occurrence["compound_label"]["evidence_ids"] + [E[7], I[7]])
    if mid in {"m_4p", "m_4q", "m_4r", "m_4s", "m_4t", "m_4u"}:
        add_page7_image_where_text_is_cited(molecule)


# Page 7 finishes 4p and gives the characterization paragraphs for 4q--4u.
# These are the same Table 4 experiments, so enrich rather than duplicate them.
hplc = {
    4: {"product": "m_4p", "column": "Chiralpak OJ-H", "wavelength_nm": 254,
        "hexane_isopropanol": "98:2", "flow_mL_min": 1.0,
        "retention_times_min": {"minor": 8.08, "major": 8.75}, "mass_mg": None},
    5: {"product": "m_4q", "column": "Chiralpak IC", "wavelength_nm": 254,
        "hexane_isopropanol": "97:3", "flow_mL_min": 1.0,
        "retention_times_min": {"minor": 7.12, "major": 8.45}, "mass_mg": 40.8},
    6: {"product": "m_4r", "column": "Chiralpak OJ-H", "wavelength_nm": 254,
        "hexane_isopropanol": "95:5", "flow_mL_min": 1.0,
        "retention_times_min": {"major": 9.73, "minor": 11.42}, "mass_mg": 46.9},
    7: {"product": "m_4s", "column": "Chiralpak OJ-H", "wavelength_nm": 220,
        "hexane_isopropanol": "95:5", "flow_mL_min": 1.0,
        "retention_times_min": {"minor": 9.33, "major": 12.20}, "mass_mg": 49.5},
    8: {"product": "m_4t", "column": "Chiralpak OJ-H", "wavelength_nm": 220,
        "hexane_isopropanol": "90:10", "flow_mL_min": 1.0,
        "retention_times_min": {"major": 29.02, "minor": 32.16}, "mass_mg": 47.8},
    9: {"product": "m_4u", "column": "Chiralpak OJ-H", "wavelength_nm": 254,
        "hexane_isopropanol": "95:5", "flow_mL_min": 1.0,
        "retention_times_min": {"minor": 9.28, "major": 9.94}, "mass_mg": 50.2},
}
updated_reactions = []
for entry, method in hplc.items():
    rid = f"{DOC}:table4_allenone_scope_entry_{entry:02d}"
    reaction = next(r for r in bundle["reactions"] if r["reaction_id"] == rid)
    add_page7_image_where_text_is_cited(reaction)
    reaction["evidence_ids"] = uniq(reaction["evidence_ids"] + [E[7], I[7]])
    for o in reaction["outcomes"]:
        if o["outcome_id"] == "product_ee":
            o["measurement_basis"] = f(
                "HPLC analysis on the stated chiral stationary phase",
                {"method": "HPLC on a chiral stationary phase", **{k: v for k, v in method.items()
                 if k not in {"product", "mass_mg"}}},
                (E[7], I[7]))
            o["scope"] = f(
                {"metric": "ee", "peak_assignment": method["retention_times_min"]},
                {"metric": "ee", "peak_assignment": method["retention_times_min"]},
                (E[7], I[7]))
    reaction["outcomes"] = [o for o in reaction["outcomes"]
                            if o["outcome_id"] != "product_mass_page7"]
    if method["mass_mg"] is not None:
        reaction["outcomes"].append(outcome(
            "product_mass_page7", "product_amount", [method["product"]],
            (E[7], I[7]), q(method["mass_mg"], "mg", f"{method['mass_mg']} mg"),
            product_role="isolated main product", basis="isolated material",
            scope={"table": 4, "entry": entry, "characterization_paragraph": True}))
    updated_reactions.append(copy.deepcopy(reaction))


# Build the single new experiment on the bounded pages.
ev = (E[7], I[7])
r = st.new_reaction(SCALE_ID, DOC, ev)
r["record_type"] = "specific_experiment"
r["experiment_categories"] = f(["scale_up"], ["scale_up"], ev)
r["location"] = f("page 7, A 5 mmol Scale Reaction", "page 7, A 5 mmol Scale Reaction", ev)
r["general_scheme_ids"] = [f"{DOC}:generic_aza_claisen_hydrolysis"]
r["reaction_name"] = f("A 5 mmol Scale Reaction", "aza-Claisen rearrangement followed by hydrolysis", ev)
r["transformation_type"] = f(
    "allylamine 1a and allenone 2a to δ-chiral β-diketone 4a",
    "enantioenriched allylamine + allenone to delta-chiral beta-diketone",
    ev, method="inference",
    derivation="Transformation label is summarized from the named participants and isolated product.")
r["participants"] = [
    participant("allylamine_1a", "m_1a", "reactant", ev,
                mass=q(1.01, "g", "1.01 g"), amount=q(5.0, "mmol", "5.0 mmol")),
    participant("allenone_2a", "m_2a", "reactant", ev,
                mass=q(865, "mg", "865 mg"), amount=q(6.0, "mmol", "6.0 mmol")),
    participant("methanol", "m_meoh", "solvent", ev,
                volume=q(25, "mL", "25 mL")),
    participant("wet_dcm", "m_dcm", "solvent", ev, stage="workup",
                volume=q(50, "mL", "50 mL")),
    participant("water_in_dcm", "m_water", "workup_agent", ev, stage="workup"),
    participant("silica_hydrolysis", "m_silica", "workup_agent", ev, stage="workup",
                mass=q(12.5, "g", "12.5 g")),
    participant("dcm_wash", "m_dcm", "workup_agent", ev, stage="workup",
                volume=q(75, "mL", "3 × 25 mL")),
    participant("petroleum_ether", "m_petroleum_ether", "purification_agent", ev,
                stage="purification"),
    participant("ethyl_acetate", "m_etoac", "purification_agent", ev,
                stage="purification"),
    participant("product_4a", "m_4a", "product", ev,
                product_role="isolated main product"),
]
r["reaction_scale"] = f(
    "allylamine 1a (5.0 mmol)",
    {"reference_molecule_id": "m_1a", "amount": q(5.0, "mmol")}, ev)
r["solvent_conditions"]["composition"] = f(
    "methanol; then dichloromethane containing 1.0% (v/v) water",
    ["m_meoh", "m_dcm", "m_water"], ev)
r["solvent_conditions"]["mixture_ratio"] = f(
    "dichloromethane containing 1.0% (v/v) water",
    {"water": q(1.0, "%"), "basis": "v/v", "stage": "wet-silica hydrolysis"}, ev)
r["solvent_conditions"]["total_volume"] = f(
    "methanol (25 mL); wet dichloromethane (50 mL); separate DCM wash 3 × 25 mL",
    {"reaction_solvent": q(25, "mL"), "hydrolysis_solvent": q(50, "mL"),
     "wash_solvent": q(75, "mL", "3 × 25 mL")}, ev)
r["steps"] = [
    step("addition", "reaction", ["allylamine_1a", "methanol", "allenone_2a"], ev,
         "allenone 2a was added to a stirred solution of allylamine 1a in methanol"),
    step("rearrangement_hold", "reaction", ["allylamine_1a", "allenone_2a", "methanol"], ev,
         "the mixture was stirred", "room temperature", q(12, "h")),
    step("remove_methanol", "workup", ["methanol"], ev,
         "concentrated under reduced pressure to remove methanol", pressure="reduced pressure"),
    step("prepare_hydrolysis", "workup", ["wet_dcm", "water_in_dcm", "silica_hydrolysis"], ev,
         "the residue was diluted with wet dichloromethane and silica gel was added"),
    step("hydrolysis_hold", "workup", ["wet_dcm", "water_in_dcm", "silica_hydrolysis", "product_4a"], ev,
         "the mixture was stirred", "room temperature", q(12, "h")),
    step("filter_and_wash", "workup", ["silica_hydrolysis", "dcm_wash", "product_4a"], ev,
         "filtered; silica gel was washed with dichloromethane (3 × 25 mL)"),
    step("concentrate_filtrates", "workup", ["product_4a"], ev,
         "combined filtrates were concentrated under reduced pressure", pressure="reduced pressure"),
    step("chromatography", "purification", ["petroleum_ether", "ethyl_acetate", "product_4a"], ev,
         "silica gel chromatography, petroleum ether/ethyl acetate = 10:1"),
]
r["procedure"]["original_text"] = f(
    "To a stirred solution of allylamine 1a (1.01 g, 5.0 mmol) in methanol (25 mL) was added allenone 2a (865 mg, 6.0 mmol). The mixture was stirred at room temperature for 12 h and concentrated under reduced pressure to remove the methanol. The residue was diluted with dichloromethane (50 mL, containing 1.0% (v/v) water), and silica gel (12.5 g) was added. The mixture was stirred at room temperature for 12 h and then filtered, and the silica gel was washed with dichloromethane (3 × 25 mL). The combined filtrates were concentrated under reduced pressure, and the residue was purified by silica gel chromatography (petroleum ether:ethyl acetate = 10:1) to give compound 4a as a yellow solid (1.27 g, 87%).",
    "5 mmol scale reaction procedure as printed", ev)
r["procedure"]["general_procedure_reference"] = f(
    "A 5 mmol Scale Reaction", "A 5 mmol Scale Reaction", ev)
r["procedure"]["addition_order"] = f(
    "allenone 2a was added to a stirred solution of allylamine 1a in methanol",
    ["allylamine 1a in methanol", "add allenone 2a"], ev)
r["outcomes"] = [
    outcome("isolated_product_mass", "product_amount", ["m_4a"], ev,
            q(1.27, "g", "1.27 g"), product_role="isolated main product",
            basis="isolated after silica gel chromatography",
            scope={"start": "1a (5.0 mmol)", "end": "4a"}),
    outcome("isolated_yield", "yield", ["m_4a"], ev,
            q(87, "%", "87%"), product_role="isolated main product",
            yield_type="isolated yield", basis="isolated after silica gel chromatography",
            scope={"start": "1a (5.0 mmol)", "end": "4a", "scale_up": True}),
]
classification_derivation = (
    "Partial net-bond interpretation for the same 1a/2a to 4a transformation shown in the paper. "
    "Atom mapping, hydrogen/coproduct balance, and elementary mechanism are not asserted.")
changes = [
    {"change_id": "b1", "atoms": [{"element": "C"}, {"element": "C"}],
     "before_order": 0, "after_order": 1, "change_type": "formed"},
    {"change_id": "b2", "atoms": [{"element": "C"}, {"element": "O"}],
     "before_order": 0, "after_order": 2, "change_type": "formed"},
    {"change_id": "b3", "atoms": [{"element": "C"}, {"element": "N"}],
     "before_order": 1, "after_order": 0, "change_type": "broken"},
    {"change_id": "b4", "atoms": [{"element": "C"}, {"element": "C"}],
     "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
]
labels = [
    {"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1"]},
    {"category_number": 14, "code": "DB-C-O", "supporting_change_ids": ["b2"]},
]
r["classification"] = {
    "standard_id": st.STANDARD,
    "mapped_reaction_smiles": missing(ev),
    "bond_changes": [derived("source transformation", change, ev, classification_derivation)
                     for change in changes],
    "labels": derived("source transformation", labels, ev, classification_derivation),
    "coverage": f("partial net-bond interpretation", "partial", ev, method="inference",
                  derivation=classification_derivation),
    "limitations": [classification_derivation, "No three-model votes were performed."],
    "model_votes": [],
}
r["notes"] = [
    f("The scale-up paragraph reports no ee, er, dr, or E/Z measurement for this isolated batch.",
      "No scale-specific stereoselectivity measurement reported.", ev, method="inference",
      derivation="The bounded paragraph ends with product mass and isolated yield only."),
]
r["review"] = {
    "status": "needs_review", "method": "model", "reviewer": None,
    "checked_fields": [],
    "unresolved": [
        "Machine record; source structure transcription and partial bond classification require expert review.",
        "The scale-up paragraph does not report a scale-specific ee, er, dr, E/Z ratio, atmosphere, concentration, or water volume as a separate charge.",
        "The named diketone graph is reused from the existing molecule registry and remains a candidate representation; the authors state products 4 are predominantly enol in NMR spectra.",
    ],
}
r["completeness"] = "partial"

# Idempotently replace any prior version of this bounded record.
bundle["reactions"] = [x for x in bundle["reactions"] if x["reaction_id"] != SCALE_ID]
bundle["reactions"].append(r)


contexts = [
    {
        "context_id": f"{DOC}:page7_characterization_4p_4u",
        "document_id": DOC,
        "kind": "characterization",
        "evidence_ids": [E[7], I[7]],
        "description": "Continuation of 4p and product-specific physical form, mass/yield, NMR, HRMS, optical rotation and chiral-HPLC data for 4q-4u. These observations enrich existing Table 4 experiments and are not separate syntheses.",
        "linked_reaction_ids": [f"{DOC}:table4_allenone_scope_entry_{entry:02d}" for entry in range(4, 10)],
    },
    {
        "context_id": f"{DOC}:page7_supporting_information_statement",
        "document_id": DOC,
        "kind": "source_context",
        "evidence_ids": [E[7], I[7]],
        "description": "The article states that X-ray crystallographic data and copies of NMR spectra/LC traces are available as Supporting Information. No SI file was supplied or explicitly associated in this workspace, so no SI evidence was imported.",
    },
    {
        "context_id": f"{DOC}:page8_references",
        "document_id": DOC,
        "kind": "bibliography",
        "evidence_ids": [E[8], I[8]],
        "description": "Reference list entries 1-9: cited prior literature, not experimental records performed in this article.",
    },
    {
        "context_id": f"{DOC}:page9_references",
        "document_id": DOC,
        "kind": "bibliography",
        "evidence_ids": [E[9], I[9]],
        "description": "Continuation of reference list entries 9-19: cited prior literature, not experimental records performed in this article.",
    },
]
context_ids = {c["context_id"] for c in contexts}
bundle["source_context"] = [c for c in bundle["source_context"]
                            if c.get("context_id") not in context_ids] + contexts


coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 7,
            "state": "attempted",
            "regions": [
                {
                    "region": "upper-left continuation of the 4p paragraph",
                    "finding": "Completes 4p NMR/HRMS and reports 97% ee with the Chiralpak OJ-H method; linked to existing Table 4 entry 4, not a new synthesis.",
                    "evidence_ids": [E[7], I[7]],
                },
                {
                    "region": "left column, characterization paragraphs for 4q-4u",
                    "finding": "Physical form, isolated mass/yield, spectra, optical rotation and 97% ee/HPLC methods for 4q-4u. The five product masses are added to existing Table 4 entries 5-9; spectra remain source_context.",
                    "evidence_ids": [E[7], I[7]],
                },
                {
                    "region": "bottom-left heading and upper-right continuation, A 5 mmol Scale Reaction",
                    "finding": "One new scale-up experiment: 1a (1.01 g, 5.0 mmol) and 2a (865 mg, 6.0 mmol) in MeOH (25 mL), followed by wet-silica hydrolysis/workup and chromatography, giving 4a (1.27 g, 87%).",
                    "evidence_ids": [E[7], I[7]],
                },
                {
                    "region": "right column, Associated Content through Acknowledgments",
                    "finding": "Supporting-information availability, CCDC accession, author metadata and acknowledgments; no supplied SI experiment is associated and no additional synthesis record is created.",
                    "evidence_ids": [E[7], I[7]],
                },
            ],
            "record_ids": [
                f"{DOC}:table4_allenone_scope_entry_{entry:02d}" for entry in range(4, 10)
            ] + [SCALE_ID],
            "unresolved": [
                "The 4p characterization paragraph begins on page 6; page 7 only supplies its continuation, and this coverage file does not claim page 6 was re-covered.",
                "No scale-specific ee/er/dr/E/Z measurement is reported for the 5 mmol batch.",
                "The article mentions SI, but no SI file was supplied or explicitly associated in this session.",
            ],
        },
        {
            "page": 8,
            "state": "attempted",
            "regions": [
                {
                    "region": "full page, References 1-9",
                    "finding": "Bibliography only. The cited reviews and prior reaction papers are literature references, not current-paper experiments, schemes, tables or reaction arrows.",
                    "evidence_ids": [E[8], I[8]],
                }
            ],
            "record_ids": [],
            "unresolved": [
                "Individual cited-paper experiments were not extracted because the cited papers themselves were not supplied as source documents for this bounded task."
            ],
        },
        {
            "page": 9,
            "state": "attempted",
            "regions": [
                {
                    "region": "left column, continuation of References 9-19; remaining page blank",
                    "finding": "Bibliography only. No current-paper synthesis experiment, no-reaction control, mechanism diagram, table or reaction arrow is present.",
                    "evidence_ids": [E[9], I[9]],
                }
            ],
            "record_ids": [],
            "unresolved": [
                "Individual cited-paper experiments were not extracted because the cited papers themselves were not supplied as source documents for this bounded task."
            ],
        },
    ],
}


diagnostic_ids = {"m_1a", "m_2a", "m_4a", "m_4p", "m_4q", "m_4r", "m_4s", "m_4t", "m_4u"}
diagnostics = {}
for molecule in bundle["molecules"]:
    if molecule["molecule_id"] in diagnostic_ids:
        smiles = st.value(molecule["smiles"])
        diagnostics[molecule["molecule_id"]] = st.stereo_inspect(smiles) if smiles else {
            "source_verified": False,
            "limitations": ["No candidate SMILES is available; source stereochemistry remains in registry facts."],
        }

raw_delta = {
    "batch_id": "allene-stereo-v31-200-20261017",
    "document_id": DOC,
    "bounded_pages": [7, 8, 9],
    "semantic_accuracy_verified": False,
    "new_reactions": [r],
    "updated_reactions": updated_reactions,
    "source_context_additions": contexts,
    "molecule_registry_ids_reused": sorted(diagnostic_ids),
    "note": "Machine-authored candidates after direct visual inspection of rendered source pages 7-9; all formal reaction records remain needs_review.",
}

dump(RUN / "raw_model/pages_0007_0009_extraction.json", raw_delta)
dump(RUN / "raw_model/stereo_diagnostics_pages_0007_0009.json", diagnostics)
dump(RUN / "coverage/pages_0007_0009.json", coverage)
dump(RUN / "submitted.json", bundle)
print(json.dumps({
    "new_reactions": 1,
    "updated_existing_reactions": len(updated_reactions),
    "source_context_additions": len(contexts),
    "total_molecules": len(bundle["molecules"]),
    "total_reactions": len(bundle["reactions"]),
    "total_source_context": len(bundle["source_context"]),
    "pages": [7, 8, 9],
}, ensure_ascii=False))
