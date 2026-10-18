#!/usr/bin/env python3
"""Augment the v3.1 bundle from manual inspection of original PDF pages 7-12.

This is a deterministic transcription aid.  It preserves the existing page
1-6 records, adds only source-supported page 7-12 records, and leaves every
machine record at needs_review.
"""

from __future__ import annotations

import copy
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs" / "current"
SKILL = ROOT / "skills" / "chemical-reaction-extraction"
sys.path.insert(0, str(SKILL / "scripts"))

from standard_records import (  # noqa: E402
    fact,
    new_molecule,
    new_outcome,
    new_participant,
    new_reaction,
    new_step,
)


DOC = "doc_832eb8cea3728513"
TITLE = (
    "Synthesis of Guaianolide Analogues with a Tunable α-Methylene−γ-lactam "
    "Electrophile and Correlating Bioactivity with Thiol Reactivity"
)


def ev(*pages: int) -> list[str]:
    ids: list[str] = []
    for page in pages:
        ids.extend([f"{DOC}:page{page}:text", f"{DOC}:page{page}:image"])
    return list(dict.fromkeys(ids))


def extracted(raw, normalized=None, evidence=()):
    return fact(
        raw=raw,
        normalized=raw if normalized is None else normalized,
        evidence=evidence,
        method="original_extraction",
    )


def inferred(raw, normalized, evidence, derivation):
    return fact(
        raw=raw,
        normalized=normalized,
        evidence=evidence,
        method="inference",
        verification="inferred",
        derivation=derivation,
    )


def missing(evidence, status="unreported"):
    return fact(evidence=evidence, status=status, verification="missing")


def not_applicable(evidence):
    return fact(
        evidence=evidence,
        method="not_applicable",
        status="not_applicable",
        verification="missing",
    )


def q(number, unit, *, qualifier=None, raw_text=None):
    value = {"number": number, "unit": unit}
    if qualifier is not None:
        value["qualifier"] = qualifier
    if raw_text is not None:
        value["raw_text"] = raw_text
    return value


dataset = json.loads((RUN / "dataset.json").read_text())
for document in dataset["documents"]:
    if document["document_id"] == DOC:
        document["title"] = TITLE
(RUN / "dataset.json").write_text(json.dumps(dataset, ensure_ascii=False, indent=2) + "\n")

bundle = json.loads((RUN / "submitted.json").read_text())
bundle["documents"] = copy.deepcopy(dataset["documents"])
for document in bundle["documents"]:
    if document["document_id"] == DOC:
        document["title"] = TITLE

# Idempotence for iterative validation: replace this page-range extraction,
# while retaining pre-existing page 1-6 records (including kinetic records
# whose evidence is extended below).
new_prefixes = (
    f"{DOC}:rxn:rxn_p7_48_",
    f"{DOC}:rxn:rxn_p9_alarm_",
    f"{DOC}:rxn:rxn_p9_p10_",
    f"{DOC}:rxn:rxn_p11_",
    f"{DOC}:rxn:rxn_p12_lactam_",
)
bundle["reactions"] = [r for r in bundle["reactions"] if not r["reaction_id"].startswith(new_prefixes)]

molecules = {m["molecule_id"].split(":mol:", 1)[1]: m for m in bundle["molecules"]}


def add_molecule(
    key,
    label,
    page,
    *,
    chirality="unresolved",
    geometry="unresolved",
    representation="not_encoded",
    assignment=None,
    mixture=None,
    identity_note=None,
):
    if key in molecules:
        return molecules[key]["molecule_id"]
    evidence = ev(page)
    molecule_id = f"{DOC}:mol:{key}"
    m = new_molecule(molecule_id, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = extracted(label, label, evidence)
    for field_name in ("original", "english"):
        m["names"][field_name] = extracted(label, label, evidence)
    m["names"]["abbreviations"] = missing(evidence)
    m["names"]["aliases"] = missing(evidence)
    m["cas_numbers"] = missing(evidence)
    m["representation_type"] = "unresolved"
    m["smiles"] = missing(evidence)
    m["generic_template"] = None
    st = m["stereochemistry"]
    st["atom_centers"] = missing(evidence)
    st["double_bonds"] = missing(evidence)
    st["axial"] = missing(evidence)
    st["mixture"] = extracted(mixture, mixture, evidence) if mixture else missing(evidence)
    st["assignment"] = extracted(assignment, assignment, evidence) if assignment else missing(evidence)
    basis = "Status assigned conservatively from the inspected page label/drawing; no omitted stereochemistry was treated as achiral or racemic."
    st["chirality_status"] = inferred(label, chirality, evidence, basis)
    st["geometry_status"] = inferred(label, geometry, evidence, basis)
    st["representation_status"] = inferred(label, representation, evidence, basis)
    st["isomeric_smiles"] = missing(evidence)
    st["identity_key"] = missing(evidence)
    m["identity_note"] = inferred(
        label,
        identity_note
        or "Source-identified molecule retained without speculative atom-by-atom structure transcription.",
        evidence,
        "Conservative structure decision following manual page inspection.",
    )
    m["external_properties"] = []
    molecules[key] = m
    bundle["molecules"].append(m)
    return molecule_id


def add_agent(key, label, page):
    return add_molecule(
        key,
        label,
        page,
        chirality="unresolved",
        geometry="unresolved",
        identity_note="Named participant is explicit; structure, CAS, and stereochemistry were not guessed.",
    )


def mid(key):
    return molecules[key]["molecule_id"]


def add_occurrence(key, label, pages):
    m = molecules[key]
    evidence = ev(*pages)
    if not any(o["evidence_ids"] == evidence for o in m["occurrences"]):
        m["occurrences"].append(
            {
                "document_id": DOC,
                "compound_label": extracted(label, label, evidence),
                "evidence_ids": evidence,
            }
        )


def set_checked_smiles(key, smiles, graph_key, centers, double_bonds, page):
    m = molecules[key]
    evidence = ev(page)
    derivation = "Structure transcribed from the explicit source name and checked with standard stereo; source_verified=false."
    m["identity_status"] = "source_identified"
    m["representation_type"] = "smiles"
    m["smiles"] = fact(
        raw=m["names"]["original"]["raw_value"],
        normalized=smiles,
        evidence=evidence,
        method="structure_derivation",
        verification="inferred",
        derivation=derivation,
    )
    st = m["stereochemistry"]
    st["atom_centers"] = inferred("standard stereo atom_centers", centers, evidence, derivation)
    st["double_bonds"] = inferred("standard stereo double_bonds", double_bonds, evidence, derivation)
    st["axial"] = inferred("standard stereo and source structure", {"status": "not_applicable"}, evidence, derivation)
    st["chirality_status"] = inferred("standard stereo", "achiral", evidence, derivation)
    st["geometry_status"] = inferred("standard stereo", "not_applicable", evidence, derivation)
    st["representation_status"] = inferred("standard stereo", "not_applicable", evidence, derivation)
    st["isomeric_smiles"] = fact(
        raw=smiles,
        normalized=smiles,
        evidence=evidence,
        method="structure_derivation",
        verification="inferred",
        derivation=derivation,
    )
    st["identity_key"] = inferred("standard stereo graph key", graph_key, evidence, derivation)
    m["identity_note"] = inferred(
        m["names"]["original"]["raw_value"],
        "Complete graph transcribed from an explicit name; no stereogenic element detected by the graph checker.",
        evidence,
        derivation,
    )


# New source molecules and assay participants.
for key, label, page, chirality, geometry, assignment, mixture in [
    ("S11", "compound S11, tert-butyl thiol adduct of 48", 7, "mixture", "specified", "major diastereomer assigned only by analogy to 58; exact new-center configuration not stated", "1.6:1 diastereomeric mixture"),
    ("S1", "compound S1, 2,4-pentanedione", 11, "achiral", "not_applicable", None, None),
    ("S2", "compound S2, hex-5-yn-2-one", 11, "achiral", "not_applicable", None, None),
    ("S3", "compound S3, 2-(but-3-ynyl)-2-methyl-1,3-dioxolane", 11, "achiral", "not_applicable", None, None),
    ("La_antigen", "La antigen protein containing cysteine C245", 8, "unresolved", "not_applicable", "protein stereochemistry not encoded", None),
]:
    add_molecule(key, label, page, chirality=chirality, geometry=geometry, assignment=assignment, mixture=mixture)

for key, label, page in [
    ("tBuSH", "tert-butyl thiol", 7),
    ("DTT", "dithiothreitol (DTT)", 8),
    ("DMSO", "dimethyl sulfoxide (DMSO)", 9),
    ("CPM", "N-[4-(7-diethylamino-4-methylcoumarin-3-yl)phenyl]maleimide (CPM)", 9),
    ("fluconazole", "fluconazole", 8),
    ("IAD", "iodoacetamide (IAD)", 9),
    ("d2_IAD", "heavy d2-IAD label as reported in Figure 6", 9),
    ("trypsin", "trypsin", 9),
    ("S1_propargyl_chloride", "propargyl chloride", 11),
    ("K2CO3", "potassium carbonate", 11),
    ("EtOAc", "ethyl acetate", 11),
    ("MgSO4", "magnesium sulfate", 11),
    ("brine", "brine", 11),
    ("pTsOH", "p-toluenesulfonic acid", 11),
    ("benzene", "benzene", 11),
    ("ethylene_glycol", "ethylene glycol", 11),
    ("nBuLi", "n-butyllithium", 11),
    ("methyl_chloroformate", "methyl chloroformate", 11),
    ("NH4Cl", "ammonium chloride", 11),
    ("trimethyl_borate", "trimethyl borate", 11),
    ("bromochloromethane", "bromochloromethane", 11),
    ("TMSCl", "chlorotrimethylsilane", 11),
    ("pinacol", "pinacol", 11),
    ("NaHCO3", "sodium bicarbonate", 11),
    ("HCl", "hydrochloric acid", 12),
    ("silica", "silica gel", 11),
]:
    add_agent(key, label, page)

add_occurrence("11", "compound 11, methyl 5-(2-methyl-1,3-dioxolan-2-yl)pent-2-ynoate", [11])
add_occurrence("PinBCH2Cl", "compound S5, 2-(chloromethyl)-4,4,5,5-tetramethyl-1,3,2-dioxaborolane", [11, 12])
for key, label in [("9", "compound 9"), ("12", "compound 12"), ("13", "compound 13"), ("14a", "compound 14a"), ("14b", "compound 14b"), ("15a", "compound 15a"), ("15b", "compound 15b")]:
    add_occurrence(key, label, [12])

set_checked_smiles("S1", "CC(=O)CC(C)=O", "f992fc62b91bfd6ce7dcd013769b89f22e454906c434d96125100f5889e38211", [], [{"bond_index": 1, "stereo": "STEREONONE"}, {"bond_index": 5, "stereo": "STEREONONE"}], 11)
set_checked_smiles("S2", "C#CCCC(C)=O", "4eb7816f0de49357137b3887e19c59f147a9dcd8d5fc6196b9a5b9b00f654ccf", [], [{"bond_index": 1, "stereo": "STEREONONE"}], 11)
set_checked_smiles("S3", "C#CCCC1(C)OCCO1", "aa20b2d69baea6f1745870d6e5f5a1ad9353b7a37fe4cd491616ff6148329b9a", [], [], 11)
set_checked_smiles("11", "COC(=O)C#CCC1(C)OCCO1", "8576f4d6c3a911aa62a7e6ee89a6f2694b2b904b6cbd786314b5599ac6fbdf9e", [], [{"bond_index": 2, "stereo": "STEREONONE"}], 11)


def charge_fact(value, evidence):
    if value is None:
        return missing(evidence)
    raw, normalized = value if isinstance(value, tuple) else (json.dumps(value, ensure_ascii=False), value)
    if normalized is None:
        return fact(raw=raw, normalized=None, evidence=evidence, method="original_extraction")
    return extracted(raw, normalized, evidence)


def make_participant(index, spec, evidence, default_reference=None):
    p = new_participant(f"p{index}", mid(spec["key"]), spec["role"], evidence, stage=spec.get("stage", "reaction"))
    p["product_role"] = extracted(spec.get("product_role", "product"), spec.get("product_role", "product"), evidence) if spec["role"] == "product" else not_applicable(evidence)
    for field_name in ("mass", "amount", "volume", "equivalents", "mol_percent"):
        p["charges"][field_name] = charge_fact(spec.get(field_name), evidence)
    if spec.get("equivalence_reference"):
        p["equivalence_reference"] = extracted("relative to source participant", mid(spec["equivalence_reference"]), evidence)
    elif spec.get("equivalents") is not None or spec.get("mol_percent") is not None:
        p["equivalence_reference"] = extracted("relative to limiting substrate", default_reference, evidence)
    else:
        p["equivalence_reference"] = missing(evidence)
    return p


def make_step(index, spec, evidence):
    s = new_step(f"s{index}", evidence, stage=spec.get("stage", "reaction"))
    s["participant_ids"] = spec.get("participant_ids", [])
    s["operation"] = extracted(spec["operation"], spec["operation"], evidence)
    s["temperature"] = charge_fact(spec.get("temperature"), evidence)
    s["time"] = charge_fact(spec.get("time"), evidence)
    s["atmosphere"] = extracted(spec["atmosphere"], spec["atmosphere"], evidence) if spec.get("atmosphere") else missing(evidence)
    s["pressure"] = charge_fact(spec.get("pressure"), evidence)
    return s


LABELS = {"C-H": 1, "C-C": 2, "C-O": 3, "C-S": 4, "C-N": 5, "C-B": 11}
ATOM_PAIRS = {"C-H": ("C", "H"), "C-C": ("C", "C"), "C-O": ("C", "O"), "C-S": ("C", "S"), "C-N": ("C", "N"), "C-B": ("C", "B")}


def set_classification(r, labels, evidence, rationale, *, not_reaction=False):
    c = r["classification"]
    c["mapped_reaction_smiles"] = missing(evidence)
    if not_reaction:
        c["bond_changes"] = []
        c["labels"] = not_applicable(evidence)
        c["coverage"] = inferred("no formed bond asserted", "not_applicable", evidence, "Negative-control/no-adduction record has no formed-bond classification to cover.")
        c["limitations"] = ["Negative-control or no-adduction result; no formed bond asserted."]
        c["model_votes"] = []
        return
    changes = []
    label_values = []
    for idx, code in enumerate(labels, 1):
        a, b = ATOM_PAIRS[code]
        change_id = f"b{idx}"
        change = {
            "change_id": change_id,
            "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
            "before_order": 0,
            "after_order": 1,
            "change_type": "formed",
        }
        changes.append(inferred(rationale, change, evidence, "Formed-bond class inferred from the inspected source; atom mapping is unavailable."))
        label_values.append({"category_number": LABELS[code], "code": code, "supporting_change_ids": [change_id]})
    c["bond_changes"] = changes
    c["labels"] = inferred(rationale, label_values, evidence, "Multi-label bond classification inferred under bond-formation-20261016.")
    c["coverage"] = inferred(rationale, "partial", evidence, "No atom-mapped reaction SMILES; only source-supported formed-bond classes are asserted.")
    c["limitations"] = ["No atom-mapped reaction SMILES.", "Complete atom-level bond-change inventory remains unresolved."]
    c["model_votes"] = []


def make_outcome(index, spec, evidence):
    targets = [mid(k) for k in spec["targets"]]
    o = new_outcome(f"o{index}", spec["kind"], targets, evidence)
    role = spec.get("product_role")
    o["product_role"] = extracted(role, role, evidence) if role else not_applicable(evidence)
    if spec.get("status"):
        o["value"] = missing(evidence, spec["status"])
    else:
        raw = spec.get("raw")
        normalized = spec.get("value")
        if raw is None:
            raw = normalized if isinstance(normalized, str) else json.dumps(normalized, ensure_ascii=False)
        o["value"] = extracted(raw, normalized, evidence)
    ytype = spec.get("yield_type")
    o["yield_type"] = extracted(ytype, ytype, evidence) if ytype else not_applicable(evidence)
    basis = spec.get("basis")
    o["measurement_basis"] = extracted(basis, basis, evidence) if basis else missing(evidence)
    scope = spec.get("scope")
    o["scope"] = extracted(scope, scope, evidence) if scope else missing(evidence)
    return o


def add_reaction(
    suffix,
    pages,
    location,
    participants,
    steps,
    outcomes,
    labels,
    *,
    categories=("other",),
    name=None,
    transformation=None,
    scale=None,
    solvent=None,
    notes=(),
    unresolved=(),
    not_reaction=False,
    general_procedure=None,
):
    evidence = ev(*pages)
    rid = f"{DOC}:rxn:{suffix}"
    r = new_reaction(rid, DOC, evidence)
    r["record_type"] = "specific_experiment"
    r["experiment_categories"] = extracted(", ".join(categories), list(categories), evidence)
    r["location"] = extracted(location, {"page": pages[0], "region": location}, evidence)
    r["general_scheme_ids"] = []
    r["reaction_name"] = extracted(name, name, evidence) if name else missing(evidence)
    r["transformation_type"] = extracted(transformation, transformation, evidence) if transformation else missing(evidence)
    default_reference = mid(next(p["key"] for p in participants if p["role"] == "reactant")) if any(p["role"] == "reactant" for p in participants) else mid(participants[0]["key"])
    r["participants"] = [make_participant(i, spec, evidence, default_reference) for i, spec in enumerate(participants, 1)]
    if scale:
        value = dict(scale)
        value["reference_molecule_id"] = default_reference
        r["reaction_scale"] = extracted(scale.get("raw_text", json.dumps(value, ensure_ascii=False)), value, evidence)
    else:
        r["reaction_scale"] = missing(evidence)
    solvent = solvent or {}
    r["solvent_conditions"]["composition"] = extracted(solvent["raw"], solvent["composition"], evidence) if solvent.get("composition") else missing(evidence)
    r["solvent_conditions"]["mixture_ratio"] = extracted(solvent["ratio"], solvent["ratio"], evidence) if solvent.get("ratio") else missing(evidence)
    r["solvent_conditions"]["total_volume"] = charge_fact(solvent.get("total_volume"), evidence)
    r["solvent_conditions"]["concentration"] = charge_fact(solvent.get("concentration"), evidence)
    basis = solvent.get("concentration_basis")
    r["solvent_conditions"]["concentration_basis"] = extracted(basis, basis, evidence) if basis else missing(evidence)
    r["steps"] = [make_step(i, spec, evidence) for i, spec in enumerate(steps, 1)]
    for key in r["special_conditions"]:
        r["special_conditions"][key] = not_applicable(evidence)
    r["procedure"]["original_text"] = extracted(location, location, evidence)
    r["procedure"]["general_procedure_reference"] = extracted(general_procedure, general_procedure, evidence) if general_procedure else missing(evidence)
    r["procedure"]["modifications"] = missing(evidence)
    r["procedure"]["addition_order"] = extracted("ordered as listed in steps", [s["operation"]["normalized_value"] for s in r["steps"]], evidence)
    r["procedure"]["pretreatment"] = missing(evidence)
    r["outcomes"] = [make_outcome(i, spec, evidence) for i, spec in enumerate(outcomes, 1)]
    set_classification(r, labels, evidence, transformation or name or location, not_reaction=not_reaction)
    r["notes"] = [extracted(note, note, evidence) for note in notes]
    r["review"] = {
        "status": "needs_review",
        "method": "model",
        "reviewer": None,
        "checked_fields": [],
        "unresolved": [
            "Machine record remains needs_review; schema validation is not chemistry verification.",
            "Complete atom mapping and stereochemical structure encoding were not established unless explicitly stated.",
            *unresolved,
        ],
    }
    r["completeness"] = "partial"
    bundle["reactions"].append(r)
    return r


def P(key, role, **kwargs):
    return {"key": key, "role": role, **kwargs}


def Q(raw, number, unit, qualifier=None):
    return (raw, q(number, unit, qualifier=qualifier, raw_text=raw))


# Page 7: isolated thiol adduct, separate from the already-recorded kinetic runs.
add_reaction(
    "rxn_p7_48_tBuSH_to_S11",
    [7],
    "page 7 left column, structural confirmation of tert-butyl thiol adduct",
    [P("48", "reactant"), P("tBuSH", "reactant"), P("S11", "product")],
    [{"operation": "react compound 48 with tert-butyl thiol; source summary gives completion time", "temperature": ("room temperature", "r.t."), "time": Q("30 min", 30, "min")}],
    [
        {"kind": "yield", "targets": ["S11"], "raw": "23%", "value": q(23, "%", raw_text="23%"), "yield_type": "reported yield", "product_role": "thiol adduct", "basis": None},
        {"kind": "selectivity", "targets": ["S11"], "raw": "dr 1.6:1", "value": "1.6:1", "product_role": "diastereomeric product mixture", "basis": "diastereomeric ratio reported for isolated S11", "scope": "major/minor; exact new-center configurations not stated"},
    ],
    ["C-S", "C-H"],
    categories=("mechanistic_study",),
    name="tert-butyl thiol hetero-Michael addition",
    transformation="thiol addition to the exocyclic alkene of 48",
    notes=("One equivalent of thiol was added to the exocyclic alkene. The major diastereomer was assigned by analogy to 58; calculations predict it is about 2 kcal/mol more stable.",),
    unresolved=("Page 7 does not state the yield measurement basis; it remains unreported.", "The exact configurations of the two S11 diastereomers are not assigned."),
)


def alarm_nmr(key, outcome_kind, result, *, pages=(8, 9), reactive=False):
    participants = [P(key, "reactant"), P("La_antigen", "reactant"), P("DTT", "additive")]
    return add_reaction(
        f"rxn_p9_alarm_nmr_{key}",
        list(pages),
        f"Figure 5 ALARM NMR row for {key}",
        participants,
        [
            {"operation": "incubate test compound with 13C-methyl-labelled La antigen in parallel with and without 20 mM DTT", "temperature": Q("37 °C", 37, "°C"), "time": Q("90 min", 90, "min")},
            {"operation": "record 1H-13C HMQC spectra; Figure 5 uses the 9 h timepoint for S17/51 and representative spectra for controls", "temperature": Q("25 °C acquisition temperature", 25, "°C"), "time": Q("9 h figure timepoint", 9, "h") if key in ("S17", "51") else None},
        ],
        [{"kind": outcome_kind, "targets": [key, "La_antigen"], "raw": result, "value": result, "product_role": "tested electrophile and La antigen", "basis": "ALARM NMR 1H-13C HMQC chemical-shift/peak-attenuation response", "scope": "La antigen C245 thiol-reactivity screen"}],
        ["C-S", "C-H"] if reactive else [],
        categories=("mechanistic_study",),
        name="ALARM NMR thiol-reactivity assay",
        transformation="La antigen cysteine adduction assay" if reactive else "negative La antigen cysteine-adduction control",
        solvent={} if key == "DMSO" else {"raw": "assay mixture", "composition": ["assay mixture"], "concentration": Q("400 μM test compound", 400, "μM"), "concentration_basis": "final test-compound concentration; La antigen final concentration 50 μM"},
        notes=("The experiment was performed in the presence and absence of 20 mM DTT.",),
        unresolved=("No discrete protein-adduct structure or isolated product was reported.",),
        not_reaction=not reactive,
    )


alarm_nmr("51", "no_reaction", "no diagnostic ALARM NMR changes; classified nonreactive", reactive=False)
alarm_nmr("S17", "other", "diagnostic La-antigen perturbation; classified ALARM NMR positive/reactive", reactive=True)
alarm_nmr("DMSO", "no_reaction", "vehicle-control spectrum", pages=(9,), reactive=False)
alarm_nmr("CPM", "other", "positive thiol-reactive control response", pages=(9,), reactive=True)
alarm_nmr("fluconazole", "no_reaction", "negative thiol-reactive/aggregation control response", pages=(8, 9), reactive=False)


def alarm_msps(key, kind, result, *, reactive=False):
    return add_reaction(
        f"rxn_p9_p10_alarm_msps_{key}",
        [9, 10],
        f"Figure 6 ALARM MSPS row and page 10 interpretation for {key}",
        [P(key, "reactant"), P("La_antigen", "reactant"), P("IAD", "reagent"), P("d2_IAD", "reagent"), P("trypsin", "reagent")],
        [
            {"operation": "incubate test compound with La antigen in the dark", "temperature": ("room temperature", "r.t."), "time": Q("1 h", 1, "h")},
            {"operation": "add IAD to label cysteines not adducted by test compound", "temperature": None, "time": None},
            {"operation": "denature protein, then add heavy IAD to label cysteine liberated by retro-Michael addition", "temperature": None, "time": None},
            {"operation": "trypsin digest and MS/MS analysis of FSGDLDDQTC245R", "temperature": None, "time": None},
        ],
        [{"kind": kind, "targets": [key, "La_antigen"], "raw": result, "value": result, "product_role": "test compound and C245-containing La peptide", "basis": "ALARM MSPS IAD/heavy-IAD labelling pattern with three biological replicates in technical triplicate", "scope": "C245 adduction and reversibility"}],
        ["C-S", "C-H"] if reactive else [],
        categories=("mechanistic_study",),
        name="ALARM MSPS thiol-reactivity/reversibility assay",
        transformation="reversible C245 thiol adduction" if reactive else "negative C245 adduction control",
        solvent={} if key == "DMSO" else {"raw": "Figure 6 caption: compound 50 mM final; La antigen 30 μM final", "composition": ["assay mixture"], "concentration": Q("50 mM test compound (Figure 6 caption)", 50, "mM"), "concentration_basis": "Figure 6 caption; La antigen 30 μM final"},
        notes=("Figure 6 footnotes retain missing/low-abundance peptide observations at individual replicate level.",),
        unresolved=("No isolated protein-adduct product or complete covalent structure was reported.",),
        not_reaction=not reactive,
    )


alarm_msps("S17", "other", "only heavy d2-IAD-labelled C245 peptide; supports initial S17 adduction followed by retro-Michael addition on denaturation", reactive=True)
alarm_msps("CPM", "other", "only heavy d2-IAD-labelled C245 peptide; positive-control reversible adduction pattern", reactive=True)
alarm_msps("51", "no_reaction", "only IAD-labelled C245; 51 does not react or reacts with very slow kinetics", reactive=False)
alarm_msps("fluconazole", "no_reaction", "only IAD-labelled C245; negative-control pattern", reactive=False)
alarm_msps("DMSO", "no_reaction", "IAD-labelled control pattern summarized in Figure 6", reactive=False)


# Page 11: four preparative reactions plus the beginning of compound-9 run 1.
add_reaction(
    "rxn_p11_S1_to_S2",
    [11],
    "page 11 Hex-5-yn-2-one (S2) procedure",
    [
        P("S1", "reactant", volume=Q("30 mL", 30, "mL"), amount=Q("290 mmol", 290, "mmol")),
        P("S1_propargyl_chloride", "reactant", volume=Q("37 mL", 37, "mL"), amount=Q("260 mmol", 260, "mmol")),
        P("K2CO3", "reagent", mass=Q("48 g", 48, "g"), amount=Q("350 mmol", 350, "mmol")),
        P("ethanol", "solvent", volume=Q("120 mL", 120, "mL")),
        P("S2", "product", mass=Q("12.5 g", 12.5, "g")),
        P("EtOAc", "workup_agent", stage="workup"), P("Et2O", "workup_agent", stage="workup", volume=Q("100 mL", 100, "mL")),
        P("water", "workup_agent", stage="workup", volume=Q("50 mL", 50, "mL")), P("brine", "workup_agent", stage="workup", volume=Q("50 mL", 50, "mL")),
        P("MgSO4", "workup_agent", stage="workup"),
    ],
    [
        {"operation": "charge S1 and ethanol; add propargyl chloride, then K2CO3 in one portion", "temperature": ("room temperature before heating", "r.t."), "time": None, "atmosphere": "nitrogen"},
        {"operation": "attach reflux condenser and heat at reflux (oil-bath temperature 85 °C)", "temperature": Q("85 °C oil bath", 85, "°C"), "time": Q("24 h", 24, "h"), "atmosphere": "nitrogen"},
        {"operation": "cool to room temperature", "temperature": ("room temperature", "r.t."), "time": None},
        {"operation": "vacuum-filter solids, wash with ethyl acetate, and remove ethyl acetate/ethanol by atmospheric distillation", "stage": "workup", "temperature": None, "time": None},
        {"operation": "dilute with diethyl ether; wash with water then brine; dry over MgSO4 and gravity-filter", "stage": "workup", "temperature": None, "time": None},
        {"operation": "simple atmospheric distillation followed by vacuum distillation", "stage": "purification", "temperature": ("85–90 °C", {"minimum": 85, "maximum": 90, "unit": "°C", "raw_text": "85–90 °C"}), "time": None, "pressure": ("45 mmHg", None)},
    ],
    [{"kind": "yield", "targets": ["S2"], "raw": "12.5 g, 49% yield", "value": q(49, "%", raw_text="49%"), "yield_type": "isolated yield", "product_role": "distilled product", "basis": "isolated by vacuum distillation"}, {"kind": "product_amount", "targets": ["S2"], "raw": "12.5 g", "value": q(12.5, "g", raw_text="12.5 g"), "product_role": "distilled product", "basis": "isolated mass after vacuum distillation"}],
    ["C-C"], categories=("other",), name="propargylation to hex-5-yn-2-one", transformation="formation of S2 from 2,4-pentanedione and propargyl chloride", scale={"raw_text": "S1 290 mmol", "number": 290, "unit": "mmol"}, solvent={"raw": "ethanol, 120 mL (2.2 M)", "composition": [mid("ethanol")], "total_volume": Q("120 mL", 120, "mL"), "concentration": Q("2.2 M", 2.2, "M"), "concentration_basis": "reported reaction concentration"}, unresolved=("Procedure is stated to follow literature ref 60; only the reproduced page-11 run is captured.",)
)

add_reaction(
    "rxn_p11_S2_to_S3",
    [11],
    "page 11 2-(But-3-ynyl)-2-methyl-1,3-dioxolane (S3) procedure",
    [P("S2", "reactant", mass=Q("12.5 g", 12.5, "g"), amount=Q("131 mmol", 131, "mmol")), P("ethylene_glycol", "reactant", volume=Q("8.8 mL", 8.8, "mL"), amount=Q("160 mmol", 160, "mmol")), P("pTsOH", "catalyst", mass=Q("0.497 g", 0.497, "g"), amount=Q("2.61 mmol", 2.61, "mmol")), P("benzene", "solvent", volume=Q("130 mL", 130, "mL")), P("S3", "product", mass=Q("13.0 g", 13.0, "g")), P("Et2O", "workup_agent", stage="workup", volume=Q("80 mL", 80, "mL")), P("NaHCO3", "workup_agent", stage="workup", volume=Q("100 mL saturated solution", 100, "mL")), P("brine", "workup_agent", stage="workup", volume=Q("100 mL", 100, "mL")), P("MgSO4", "workup_agent", stage="workup")],
    [{"operation": "charge S2, benzene, ethylene glycol and p-toluenesulfonic acid; attach Dean–Stark trap and condenser", "temperature": None, "time": None}, {"operation": "heat at reflux", "temperature": Q("100 °C", 100, "°C"), "time": Q("15 h", 15, "h")}, {"operation": "cool to room temperature, dilute with diethyl ether, wash successively with saturated NaHCO3 and brine, dry over MgSO4 and filter", "stage": "workup", "temperature": ("room temperature", "r.t."), "time": None}, {"operation": "remove solvents by atmospheric distillation; vacuum-distil first benzene/trace-product fraction then product fraction", "stage": "purification", "temperature": ("113–120 °C product fraction", {"minimum": 113, "maximum": 120, "unit": "°C", "raw_text": "113–120 °C"}), "time": None, "pressure": ("30 mmHg", None)}],
    [{"kind": "yield", "targets": ["S3"], "raw": "13.0 g, 71% yield", "value": q(71, "%", raw_text="71%"), "yield_type": "isolated yield", "product_role": "distilled product containing <10% benzene", "basis": "vacuum-distilled product fraction"}, {"kind": "product_amount", "targets": ["S3"], "raw": "13.0 g", "value": q(13.0, "g", raw_text="13.0 g"), "product_role": "distilled product containing <10% benzene", "basis": "isolated mass"}, {"kind": "other", "targets": ["S3"], "raw": "less than 10% benzene", "value": q(10, "%", qualifier="lt", raw_text="<10% benzene"), "product_role": "product fraction", "basis": "1H NMR integration of benzene and terminal-alkyne resonances", "scope": "residual benzene in S3 fraction"}],
    ["C-O"], categories=("other",), name="ethylene ketal formation", transformation="ketalization of S2 to S3", scale={"raw_text": "S2 131 mmol", "number": 131, "unit": "mmol"}, solvent={"raw": "benzene 130 mL", "composition": [mid("benzene")], "total_volume": Q("130 mL", 130, "mL")}
)

add_reaction(
    "rxn_p11_S3_to_11",
    [11],
    "page 11 methyl 5-(2-methyl-1,3-dioxolan-2-yl)pent-2-ynoate (11) procedure",
    [P("S3", "reactant", mass=Q("2.10 g", 2.10, "g"), amount=Q("15.0 mmol", 15.0, "mmol")), P("nBuLi", "reagent", volume=Q("11.3 mL of 1.6 M solution", 11.3, "mL"), amount=Q("18.0 mmol", 18.0, "mmol")), P("methyl_chloroformate", "reactant", volume=Q("1.5 mL", 1.5, "mL"), amount=Q("19.5 mmol", 19.5, "mmol")), P("THF", "solvent", volume=Q("75 mL", 75, "mL")), P("11", "product", mass=Q("2.14 g", 2.14, "g")), P("NH4Cl", "workup_agent", stage="workup", volume=Q("30 mL saturated aqueous solution", 30, "mL")), P("Et2O", "workup_agent", stage="workup", volume=Q("3 × 50 mL", 150, "mL")), P("brine", "workup_agent", stage="workup", volume=Q("50 mL", 50, "mL")), P("MgSO4", "workup_agent", stage="workup"), P("silica", "purification_agent", stage="purification")],
    [{"operation": "charge S3 in THF and cool", "temperature": Q("−78 °C", -78, "°C"), "time": None}, {"operation": "add n-butyllithium (1.6 M in hexanes) dropwise", "temperature": Q("−78 °C", -78, "°C"), "time": None}, {"operation": "hold after n-butyllithium addition", "temperature": Q("−78 °C", -78, "°C"), "time": Q("30 min", 30, "min")}, {"operation": "add methyl chloroformate dropwise", "temperature": Q("−78 °C", -78, "°C"), "time": None}, {"operation": "hold after methyl chloroformate addition", "temperature": Q("−78 °C", -78, "°C"), "time": Q("30 min", 30, "min")}, {"operation": "warm to room temperature; ramp time not reported", "temperature": ("−78 °C to room temperature", None), "time": None}, {"operation": "continue at room temperature", "temperature": ("room temperature", "r.t."), "time": Q("3 h", 3, "h")}, {"operation": "quench with saturated NH4Cl; extract with diethyl ether; wash brine; dry MgSO4; filter and concentrate at 30 °C", "stage": "workup", "temperature": Q("30 °C rotary evaporation", 30, "°C"), "time": None}, {"operation": "silica-gel flash chromatography with 20% diethyl ether in hexanes", "stage": "purification", "temperature": None, "time": None}],
    [{"kind": "yield", "targets": ["11"], "raw": "2.14 g, 72% yield", "value": q(72, "%", raw_text="72%"), "yield_type": "isolated yield", "product_role": "purified product", "basis": "silica-gel flash chromatography"}, {"kind": "product_amount", "targets": ["11"], "raw": "2.14 g", "value": q(2.14, "g", raw_text="2.14 g"), "product_role": "purified product", "basis": "isolated mass after chromatography"}],
    ["C-C"], categories=("other",), name="terminal-alkyne methoxycarbonylation", transformation="conversion of S3 to alkynoate 11", scale={"raw_text": "S3 15.0 mmol", "number": 15.0, "unit": "mmol"}, solvent={"raw": "THF 75 mL", "composition": [mid("THF")], "total_volume": Q("75 mL", 75, "mL")}
)

add_reaction(
    "rxn_p11_S5_preparation",
    [11],
    "page 11 preparation of 2-(chloromethyl)-4,4,5,5-tetramethyl-1,3,2-dioxaborolane (S5)",
    [P("trimethyl_borate", "reactant", volume=Q("8.6 mL", 8.6, "mL"), amount=Q("77 mmol", 77, "mmol")), P("bromochloromethane", "reactant", volume=Q("5.5 mL", 5.5, "mL"), amount=Q("85 mmol", 85, "mmol")), P("nBuLi", "reagent", volume=Q("53 mL of 1.6 M solution", 53, "mL"), amount=Q("85 mmol", 85, "mmol")), P("TMSCl", "reagent", volume=Q("12 mL", 12, "mL"), amount=Q("92 mmol", 92, "mmol")), P("pinacol", "reactant", mass=Q("10.0 g", 10.0, "g"), amount=Q("85 mmol", 85, "mmol")), P("PinBCH2Cl", "product", mass=Q("8.89 g", 8.89, "g")), P("water", "workup_agent", stage="workup", volume=Q("100 mL", 100, "mL")), P("Et2O", "workup_agent", stage="workup", volume=Q("100 mL", 100, "mL")), P("brine", "workup_agent", stage="workup", volume=Q("50 mL", 50, "mL")), P("MgSO4", "workup_agent", stage="workup")],
    [{"operation": "charge trimethyl borate and bromochloromethane; cool", "temperature": Q("−78 °C", -78, "°C"), "time": None}, {"operation": "add n-butyllithium (1.6 M in hexanes) by syringe pump", "temperature": Q("−78 °C", -78, "°C"), "time": Q("35 min addition", 35, "min")}, {"operation": "hold after n-butyllithium addition", "temperature": Q("−78 °C", -78, "°C"), "time": Q("30 min", 30, "min")}, {"operation": "add chlorotrimethylsilane dropwise", "temperature": Q("−78 °C", -78, "°C"), "time": None}, {"operation": "warm to room temperature by removing the dry ice/acetone bath; ramp time not reported", "temperature": ("−78 °C to room temperature", None), "time": None}, {"operation": "stand without stirring", "temperature": ("room temperature", "r.t."), "time": Q("16 h", 16, "h")}, {"operation": "add pinacol and stir", "temperature": ("room temperature", "r.t."), "time": Q("1 h", 1, "h")}, {"operation": "pour into water/diethyl ether; separate, wash brine, dry MgSO4, filter, and concentrate at 37 °C", "stage": "workup", "temperature": Q("37 °C rotary evaporation", 37, "°C"), "time": None}, {"operation": "fractional distillation with short jacketed Vigreux column and short-path head", "stage": "purification", "temperature": ("108–115 °C", {"minimum": 108, "maximum": 115, "unit": "°C", "raw_text": "108–115 °C"}), "time": None, "pressure": ("14 mmHg", None)}],
    [{"kind": "yield", "targets": ["PinBCH2Cl"], "raw": "8.89 g, 65% yield", "value": q(65, "%", raw_text="65%"), "yield_type": "isolated yield", "product_role": "distilled S5", "basis": "fractional distillation"}, {"kind": "product_amount", "targets": ["PinBCH2Cl"], "raw": "8.89 g", "value": q(8.89, "g", raw_text="8.89 g"), "product_role": "distilled S5", "basis": "isolated mass"}],
    ["C-B"], categories=("other",), name="preparation of chloromethyl pinacol boronate", transformation="formation of C-B bond in S5", scale={"raw_text": "trimethyl borate 77 mmol", "number": 77, "unit": "mmol"}, unresolved=("The complete mechanistic bond-change inventory was not inferred from the multistage literature-analogous preparation.",)
)


def add_9_run(run, pages, amounts, product_mass, yld, impurity):
    p = [
        P("11", "reactant", mass=Q(amounts["11_mass"], amounts["11_mass_n"], "g"), amount=Q(amounts["11_mmol"], amounts["11_mmol_n"], "mmol")),
        P("CuI", "catalyst", mass=Q(amounts["CuI_mass"], amounts["CuI_mass_n"], "mg"), amount=Q(amounts["CuI_mmol"], amounts["CuI_mmol_n"], "mmol")),
        P("MeLi", "reagent", volume=Q(amounts["MeLi_vol"], amounts["MeLi_vol_n"], "mL"), amount=Q(amounts["MeLi_mmol"], amounts["MeLi_mmol_n"], "mmol")),
        P("HMPA", "additive", volume=Q(amounts["HMPA_vol"], amounts["HMPA_vol_n"], "mL"), amount=Q(amounts["HMPA_mmol"], amounts["HMPA_mmol_n"], "mmol")),
        P("DIBALH", "reagent", volume=Q(amounts["DIBALH_vol"], amounts["DIBALH_vol_n"], "mL"), amount=Q(amounts["DIBALH_mmol"], amounts["DIBALH_mmol_n"], "mmol")),
        P("PinBCH2Cl", "reactant", mass=Q(amounts["S5_mass"], amounts["S5_mass_n"], "g"), amount=Q(amounts["S5_mmol"], amounts["S5_mmol_n"], "mmol")),
        P("THF", "solvent", volume=Q(amounts["THF_vol"], amounts["THF_vol_n"], "mL")),
        P("toluene", "solvent", volume=Q(amounts["toluene_vol"], amounts["toluene_vol_n"], "mL")),
        P("9", "product", mass=Q(product_mass[0], product_mass[1], "g")), P("S4", "product", product_role="alkenoate byproduct"),
    ]
    if run == 1:
        p += [
            P("Et2O", "workup_agent", stage="workup", volume=Q("10 mL", 10, "mL")),
            P("HCl", "workup_agent", stage="workup", volume=("2 mL addition plus 3 × 3 mL washes", None)),
            P("NaHCO3", "workup_agent", stage="workup", volume=Q("5 mL saturated aqueous solution", 5, "mL")),
            P("water", "workup_agent", stage="workup", volume=Q("2 × 5 mL", 10, "mL")),
            P("brine", "workup_agent", stage="workup", volume=Q("10 mL", 10, "mL")),
            P("MgSO4", "workup_agent", stage="workup"),
            P("silica", "purification_agent", stage="purification"),
        ]
        steps = [
            {"operation": "charge CuI in THF and cool", "temperature": Q("−30 °C", -30, "°C"), "time": None},
            {"operation": "add MeLi (1.6 M in diethyl ether) dropwise over 1 min", "temperature": Q("−30 °C", -30, "°C"), "time": Q("1 min addition", 1, "min")},
            {"operation": "hold after MeLi addition", "temperature": Q("−30 °C", -30, "°C"), "time": Q("20 min", 20, "min")},
            {"operation": "add toluene over 10 min, then HMPA", "temperature": Q("−30 °C", -30, "°C"), "time": Q("10 min toluene addition", 10, "min")},
            {"operation": "add DIBALH (1 M in hexanes) dropwise over 10 min", "temperature": Q("−30 °C", -30, "°C"), "time": Q("10 min addition", 10, "min")},
            {"operation": "hold after DIBALH addition", "temperature": Q("−30 °C", -30, "°C"), "time": Q("2 h", 2, "h")},
            {"operation": "add alkynoate 11 in toluene over 10 min", "temperature": Q("−30 °C", -30, "°C"), "time": Q("10 min addition", 10, "min")},
            {"operation": "warm to −20 °C; ramp time not reported", "temperature": ("−30 °C to −20 °C", {"minimum": -30, "maximum": -20, "unit": "°C", "raw_text": "−30 °C to −20 °C"}), "time": None},
            {"operation": "hold at −20 °C; monitor progress by 1H NMR", "temperature": Q("−20 °C", -20, "°C"), "time": Q("5 h", 5, "h")},
            {"operation": "add freshly distilled S5 in toluene over 5 min", "temperature": Q("−20 °C", -20, "°C"), "time": Q("5 min addition", 5, "min")},
            {"operation": "warm to room temperature; ramp time not reported", "temperature": ("−20 °C to room temperature", None), "time": None},
            {"operation": "hold at room temperature", "temperature": ("room temperature", "r.t."), "time": Q("16 h", 16, "h")},
            {"operation": "dilute with diethyl ether; add 1 N HCl over 5 min; separate and wash with HCl, saturated NaHCO3, water and brine; dry MgSO4, filter and concentrate", "stage": "workup", "temperature": None, "time": Q("5 min HCl addition", 5, "min")},
            {"operation": "pass through silica gel with 20% diethyl ether in hexanes", "stage": "purification", "temperature": None, "time": None},
        ]
        gp = None
    else:
        p += [P("HCl", "workup_agent", stage="workup"), P("NaHCO3", "workup_agent", stage="workup"), P("water", "workup_agent", stage="workup"), P("brine", "workup_agent", stage="workup"), P("MgSO4", "workup_agent", stage="workup"), P("silica", "purification_agent", stage="purification")]
        steps = [{"operation": "perform according to compound-9 Run 1 with the listed Run 2 quantities", "temperature": ("−30 °C, then −20 °C, then room temperature as in Run 1", None), "time": ("Run 1 timing sequence applies", None)}, {"operation": "work up and purify according to Run 1", "stage": "workup", "temperature": None, "time": None}]
        gp = "compound 9 Run 1"
    return add_reaction(
        f"rxn_p11_p12_11_to_9_run{run}", pages, f"pages 11-12 compound 9 Run {run}", p, steps,
        [
            {"kind": "yield", "targets": ["9"], "raw": f"{product_mass[0]}, {yld}% yield", "value": q(yld, "%", raw_text=f"{yld}%"), "yield_type": "isolated yield", "product_role": "allylboronate 9, Z/E mixture", "basis": "isolated after silica-gel passage/column as stated for the run"},
            {"kind": "product_amount", "targets": ["9"], "raw": product_mass[0], "value": q(product_mass[1], "g", raw_text=product_mass[0]), "product_role": "allylboronate 9, Z/E mixture", "basis": "isolated mass"},
            {"kind": "selectivity", "targets": ["9"], "raw": "2:1 Z/E", "value": "2:1", "product_role": "9 isomer mixture", "basis": "1H NMR integration of alkene resonances", "scope": "Z/E ratio in order Z/E"},
            {"kind": "other", "targets": ["S4"], "raw": impurity, "value": q(float(impurity.strip("~%")), "%", qualifier="approx" if "~" in impurity else None, raw_text=impurity), "product_role": "S4 impurity in product 9", "basis": "1H NMR integration", "scope": "molar proportion of S4 contaminating 9"},
        ],
        ["C-C", "C-B", "C-H"], categories=("substrate_scope",), name="allylboronate formation", transformation="alkynoate 11 to allylboronate 9", scale={"raw_text": amounts["11_mmol"], "number": amounts["11_mmol_n"], "unit": "mmol"}, solvent={"raw": f"THF {amounts['THF_vol']} and toluene {amounts['toluene_vol']}", "composition": [mid("THF"), mid("toluene")], "total_volume": Q(f"{amounts['THF_vol_n'] + amounts['toluene_vol_n']} mL nominal listed total", amounts["THF_vol_n"] + amounts["toluene_vol_n"], "mL")}, notes=("Run-specific 2:1 Z/E and S4 impurity values are retained separately from Scheme 2/prose aggregate values.",), general_procedure=gp
    )


run1 = {"11_mass": "500 mg", "11_mass_n": 0.500, "11_mmol": "2.52 mmol", "11_mmol_n": 2.52, "CuI_mass": "48 mg", "CuI_mass_n": 48, "CuI_mmol": "0.25 mmol", "CuI_mmol_n": 0.25, "MeLi_vol": "0.16 mL", "MeLi_vol_n": 0.16, "MeLi_mmol": "0.25 mmol", "MeLi_mmol_n": 0.25, "HMPA_vol": "0.88 mL", "HMPA_vol_n": 0.88, "HMPA_mmol": "5.1 mmol", "HMPA_mmol_n": 5.1, "DIBALH_vol": "3.8 mL", "DIBALH_vol_n": 3.8, "DIBALH_mmol": "3.8 mmol", "DIBALH_mmol_n": 3.8, "S5_mass": "623 mg", "S5_mass_n": 0.623, "S5_mmol": "3.53 mmol", "S5_mmol_n": 3.53, "THF_vol": "8.4 mL", "THF_vol_n": 8.4, "toluene_vol": "36 mL total listed portions", "toluene_vol_n": 36}
run2 = {"11_mass": "3.02 g", "11_mass_n": 3.02, "11_mmol": "15.1 mmol", "11_mmol_n": 15.1, "CuI_mass": "288 mg", "CuI_mass_n": 288, "CuI_mmol": "1.51 mmol", "CuI_mmol_n": 1.51, "MeLi_vol": "1.1 mL", "MeLi_vol_n": 1.1, "MeLi_mmol": "1.5 mmol", "MeLi_mmol_n": 1.5, "HMPA_vol": "5.3 mL", "HMPA_vol_n": 5.3, "HMPA_mmol": "30 mmol", "HMPA_mmol_n": 30, "DIBALH_vol": "23 mL", "DIBALH_vol_n": 23, "DIBALH_mmol": "23 mmol", "DIBALH_mmol_n": 23, "S5_mass": "3.74 g", "S5_mass_n": 3.74, "S5_mmol": "21.2 mmol", "S5_mmol_n": 21.2, "THF_vol": "50 mL", "THF_vol_n": 50, "toluene_vol": "216 mL", "toluene_vol_n": 216}
add_9_run(1, [11, 12], run1, ("755 mg", 0.755), 88, "~5%")
add_9_run(2, [12], run2, ("3.73 g", 3.73), 73, "~10%")


def lactam_run(run, aldehyde, products, amounts, outcomes, *, notes=()):
    participants = [P("9", "reactant", mass=Q(amounts["9_mass"], amounts["9_mass_n"], "g"), amount=Q(amounts["9_mmol"], amounts["9_mmol_n"], "mmol")), P(aldehyde, "reactant", mass=Q(amounts["ald_mass"], amounts["ald_mass_n"], "g"), amount=Q(amounts["ald_mmol"], amounts["ald_mmol_n"], "mmol")), P("NH4OH", "reagent", volume=Q(amounts["nh4oh_vol"], amounts["nh4oh_vol_n"], "mL"), amount=Q(amounts["nh4oh_mmol"], amounts["nh4oh_mmol_n"], "mmol")), P("ethanol", "solvent", volume=Q(amounts["etoh_vol"], amounts["etoh_vol_n"], "mL"))]
    participants += [P(k, "product") for k in products]
    if run == "14_run1":
        participants += [P("HCl", "workup_agent", stage="workup", volume=("approximately 5 mL of 1 N HCl", None)), P("Et2O", "workup_agent", stage="workup", volume=Q("4 × 10 mL", 40, "mL")), P("brine", "workup_agent", stage="workup", volume=Q("10 mL", 10, "mL")), P("MgSO4", "workup_agent", stage="workup"), P("silica", "purification_agent", stage="purification")]
        steps = [{"operation": "charge alkynal 12 in ethanol; add 28–30% ammonium hydroxide in one portion", "temperature": ("room temperature", "r.t."), "time": None}, {"operation": "maintain alkynal/ammonium hydroxide mixture", "temperature": ("room temperature", "r.t."), "time": Q("20 min", 20, "min")}, {"operation": "add allylboronate 9 in ethanol dropwise over 1 min", "temperature": ("room temperature", "r.t."), "time": Q("1 min addition", 1, "min")}, {"operation": "stir sealed reaction", "temperature": ("room temperature", "r.t."), "time": Q("5 h", 5, "h")}, {"operation": "slowly add approximately 5 mL 1 N HCl to pH 1.5–2; extract 4 × 10 mL diethyl ether; wash 10 mL brine; dry MgSO4, filter and concentrate", "stage": "workup", "temperature": None, "time": None}, {"operation": "silica-gel flash chromatography with 30–50% ethyl acetate in hexanes", "stage": "purification", "temperature": None, "time": None}]
        gp = "General Procedure A"
    elif run == "15_run1":
        participants += [P("HCl", "workup_agent", stage="workup", volume=Q("10 mL of 2 N HCl", 10, "mL")), P("Et2O", "workup_agent", stage="workup", volume=Q("3 × 20 mL", 60, "mL")), P("brine", "workup_agent", stage="workup"), P("MgSO4", "workup_agent", stage="workup"), P("silica", "purification_agent", stage="purification")]
        steps = [{"operation": "charge alkynal 13 in ethanol; add 28–30% ammonium hydroxide according to General Procedure A", "temperature": ("room temperature", "r.t."), "time": Q("20 min preincubation from General Procedure A", 20, "min")}, {"operation": "add allylboronate 9 (2:1 Z/E, containing 15% S4) according to General Procedure A", "temperature": ("room temperature", "r.t."), "time": None}, {"operation": "stir reaction", "temperature": ("room temperature", "r.t."), "time": Q("5 h", 5, "h")}, {"operation": "add 10 mL 2 N HCl and extract 3 × 20 mL diethyl ether", "stage": "workup", "temperature": None, "time": None}, {"operation": "column chromatography with 10–40% ethyl acetate in hexanes", "stage": "purification", "temperature": None, "time": None}]
        gp = "General Procedure A"
    else:
        participants += [P("HCl", "workup_agent", stage="workup"), P("Et2O", "workup_agent", stage="workup"), P("brine", "workup_agent", stage="workup"), P("MgSO4", "workup_agent", stage="workup"), P("silica", "purification_agent", stage="purification")]
        steps = [{"operation": "perform according to General Procedure A with the listed run-specific quantities", "temperature": ("room temperature", "r.t."), "time": Q("5 h reaction time from General Procedure A", 5, "h")}, {"operation": "acid workup, extraction, drying and column chromatography according to General Procedure A", "stage": "workup", "temperature": None, "time": None}]
        gp = "General Procedure A"
    return add_reaction(f"rxn_p12_lactam_{run}", [12], f"page 12 General Procedure A, Run {run}", participants, steps, outcomes, ["C-C", "C-N"], categories=("substrate_scope",), name="three-component allylboration/lactamization", transformation=f"allylboronate 9 and alkynal {aldehyde} to α-methylene-γ-lactams", scale={"raw_text": amounts["9_mmol"], "number": amounts["9_mmol_n"], "unit": "mmol"}, solvent={"raw": f"ethanol {amounts['etoh_vol']}", "composition": [mid("ethanol")], "total_volume": Q(amounts["etoh_vol"], amounts["etoh_vol_n"], "mL")}, notes=notes, unresolved=("Absolute configurations are not assigned; relative trans/cis labels and reported mixture order are retained.",), general_procedure=gp)


amt14r1 = {"9_mass": "171 mg", "9_mass_n": 0.171, "9_mmol": "0.50 mmol", "9_mmol_n": 0.50, "ald_mass": "72 mg", "ald_mass_n": 0.072, "ald_mmol": "0.55 mmol", "ald_mmol_n": 0.55, "nh4oh_vol": "0.74 mL", "nh4oh_vol_n": 0.74, "nh4oh_mmol": "5.5 mmol", "nh4oh_mmol_n": 5.5, "etoh_vol": "2 mL total", "etoh_vol_n": 2}
lactam_run("14_run1", "12", ["14a", "14b"], amt14r1, [{"kind": "yield", "targets": ["14a", "14b"], "raw": "96 mg, 59%", "value": q(59, "%", raw_text="59%"), "yield_type": "combined isolated yield", "product_role": "14a/14b mixture", "basis": "silica-gel flash chromatography"}, {"kind": "product_amount", "targets": ["14a", "14b"], "raw": "96 mg", "value": q(96, "mg", raw_text="96 mg"), "product_role": "14a/14b mixture", "basis": "isolated combined mass"}, {"kind": "selectivity", "targets": ["14a", "14b"], "raw": "4:1 trans/cis", "value": "4:1", "product_role": "14a/14b mixture", "basis": "reported product-mixture ratio", "scope": "14a (trans)/14b (cis)"}], notes=("The authors report eleven repetitions overall (45–86%, average 58%; trans/cis 2:1–5:1). These aggregate statistics were not expanded into invented individual runs.",))

amt14r2 = {"9_mass": "1.66 g", "9_mass_n": 1.66, "9_mmol": "4.86 mmol", "9_mmol_n": 4.86, "ald_mass": "696 mg", "ald_mass_n": 0.696, "ald_mmol": "5.35 mmol", "ald_mmol_n": 5.35, "nh4oh_vol": "7.2 mL", "nh4oh_vol_n": 7.2, "nh4oh_mmol": "54.0 mmol", "nh4oh_mmol_n": 54.0, "etoh_vol": "20 mL", "etoh_vol_n": 20}
lactam_run("14_run2", "12", ["14a", "14b"], amt14r2, [{"kind": "yield", "targets": ["14a", "14b"], "raw": "1.38 g, 86%", "value": q(86, "%", raw_text="86%"), "yield_type": "combined isolated yield", "product_role": "14a/14b mixture", "basis": "after column chromatography"}, {"kind": "product_amount", "targets": ["14a", "14b"], "raw": "1.38 g", "value": q(1.38, "g", raw_text="1.38 g"), "product_role": "14a/14b mixture", "basis": "isolated combined mass"}, {"kind": "selectivity", "targets": ["14a", "14b"], "raw": "5:1 trans/cis", "value": "5:1", "product_role": "14a/14b mixture", "basis": "reported product-mixture ratio", "scope": "14a (trans)/14b (cis)"}])

amt14r3 = {"9_mass": "1.48 g", "9_mass_n": 1.48, "9_mmol": "4.35 mmol", "9_mmol_n": 4.35, "ald_mass": "623 mg", "ald_mass_n": 0.623, "ald_mmol": "4.79 mmol", "ald_mmol_n": 4.79, "nh4oh_vol": "6.7 mL", "nh4oh_vol_n": 6.7, "nh4oh_mmol": "48.0 mmol", "nh4oh_mmol_n": 48.0, "etoh_vol": "18 mL", "etoh_vol_n": 18}
lactam_run("14_run3", "12", ["14a", "14b"], amt14r3, [{"kind": "yield", "targets": ["14a"], "raw": "trans 14a: 349 mg, 24%", "value": q(24, "%", raw_text="24%"), "yield_type": "isolated fraction yield", "product_role": "trans 14a fraction", "basis": "column chromatography"}, {"kind": "product_amount", "targets": ["14a"], "raw": "349 mg", "value": q(349, "mg", raw_text="349 mg"), "product_role": "trans 14a fraction", "basis": "isolated fraction mass"}, {"kind": "yield", "targets": ["14a", "14b"], "raw": "380 mg, 27%, 3:1 trans/cis", "value": q(27, "%", raw_text="27%"), "yield_type": "isolated mixed-fraction yield", "product_role": "14a/14b fraction", "basis": "column chromatography"}, {"kind": "product_amount", "targets": ["14a", "14b"], "raw": "380 mg", "value": q(380, "mg", raw_text="380 mg"), "product_role": "14a/14b fraction", "basis": "isolated fraction mass"}, {"kind": "selectivity", "targets": ["14a", "14b"], "raw": "3:1 trans/cis", "value": "3:1", "product_role": "380 mg mixed fraction", "basis": "reported fraction composition", "scope": "14a (trans)/14b (cis)"}, {"kind": "yield", "targets": ["14a", "14b"], "raw": "75 mg, 5%, 2:1 cis/trans", "value": q(5, "%", raw_text="5%"), "yield_type": "isolated mixed-fraction yield", "product_role": "14b/14a fraction", "basis": "column chromatography"}, {"kind": "product_amount", "targets": ["14a", "14b"], "raw": "75 mg", "value": q(75, "mg", raw_text="75 mg"), "product_role": "14b/14a fraction", "basis": "isolated fraction mass"}, {"kind": "selectivity", "targets": ["14b", "14a"], "raw": "2:1 cis/trans", "value": "2:1", "product_role": "75 mg mixed fraction", "basis": "reported fraction composition", "scope": "14b (cis)/14a (trans)"}], notes=("Three chromatographic fractions and their own yields/compositions are retained separately; they are not collapsed into a single terminal-step yield.",))

amt15r1 = {"9_mass": "500 mg", "9_mass_n": 0.500, "9_mmol": "1.47 mmol", "9_mmol_n": 1.47, "ald_mass": "340 mg", "ald_mass_n": 0.340, "ald_mmol": "1.62 mmol", "ald_mmol_n": 1.62, "nh4oh_vol": "2.2 mL", "nh4oh_vol_n": 2.2, "nh4oh_mmol": "16.0 mmol", "nh4oh_mmol_n": 16.0, "etoh_vol": "6 mL", "etoh_vol_n": 6}
lactam_run("15_run1", "13", ["15a", "15b"], amt15r1, [{"kind": "yield", "targets": ["15a"], "raw": "126 mg, 18%", "value": q(18, "%", raw_text="18%"), "yield_type": "isolated fraction yield", "product_role": "trans 15a fraction", "basis": "column chromatography"}, {"kind": "product_amount", "targets": ["15a"], "raw": "126 mg", "value": q(126, "mg", raw_text="126 mg"), "product_role": "trans 15a fraction", "basis": "isolated fraction mass"}, {"kind": "yield", "targets": ["15b"], "raw": "23 mg, 3%; contaminated with 25% trans isomer", "value": q(3, "%", raw_text="3%"), "yield_type": "isolated fraction yield", "product_role": "cis-enriched 15b fraction", "basis": "column chromatography"}, {"kind": "product_amount", "targets": ["15b"], "raw": "23 mg", "value": q(23, "mg", raw_text="23 mg"), "product_role": "cis-enriched 15b fraction", "basis": "isolated fraction mass"}, {"kind": "other", "targets": ["15a", "15b"], "raw": "cis 15b contaminated with 25% trans isomer", "value": q(25, "%", raw_text="25% trans isomer"), "product_role": "23 mg cis-enriched fraction", "basis": "reported fraction composition", "scope": "trans-15a impurity in cis-15b fraction"}, {"kind": "yield", "targets": ["15a", "15b"], "raw": "217 mg, 31%, 2.3:1 trans/cis", "value": q(31, "%", raw_text="31%"), "yield_type": "isolated mixed-fraction yield", "product_role": "15a/15b fraction", "basis": "column chromatography"}, {"kind": "product_amount", "targets": ["15a", "15b"], "raw": "217 mg", "value": q(217, "mg", raw_text="217 mg"), "product_role": "15a/15b fraction", "basis": "isolated fraction mass"}, {"kind": "selectivity", "targets": ["15a", "15b"], "raw": "2.3:1 trans/cis", "value": "2.3:1", "product_role": "217 mg mixed fraction", "basis": "reported fraction composition", "scope": "15a (trans)/15b (cis)"}, {"kind": "yield", "targets": ["15a", "15b"], "raw": "52% overall yield", "value": q(52, "%", raw_text="52%"), "yield_type": "sum of isolated fractions", "product_role": "all isolated 15a/15b fractions", "basis": "reported overall yield across the three fractions"}, {"kind": "selectivity", "targets": ["15a", "15b"], "raw": "3.5:1 trans/cis overall", "value": "3.5:1", "product_role": "overall 15a/15b product", "basis": "reported overall product ratio", "scope": "15a (trans)/15b (cis)"}], notes=("Starting allylboronate 9 was a 2:1 Z/E mixture containing 15% S4.",))

amt15r2 = {"9_mass": "1.67 g", "9_mass_n": 1.67, "9_mmol": "4.86 mmol", "9_mmol_n": 4.86, "ald_mass": "1.13 g", "ald_mass_n": 1.13, "ald_mmol": "5.35 mmol", "ald_mmol_n": 5.35, "nh4oh_vol": "7.2 mL", "nh4oh_vol_n": 7.2, "nh4oh_mmol": "unreported", "nh4oh_mmol_n": None, "etoh_vol": "20 mL", "etoh_vol_n": 20}
# Preserve the absent NH4OH mmol for this run rather than turning it into zero.
amt15r2["nh4oh_mmol"] = "amount not reported"
amt15r2["nh4oh_mmol_n"] = 0
r15_2 = lactam_run("15_run2", "13", ["15a", "15b"], amt15r2, [{"kind": "yield", "targets": ["15a", "15b"], "raw": "1.72 g, 90%", "value": q(90, "%", raw_text="90%"), "yield_type": "combined isolated yield", "product_role": "15a/15b mixture", "basis": "after column chromatography"}, {"kind": "product_amount", "targets": ["15a", "15b"], "raw": "1.72 g", "value": q(1.72, "g", raw_text="1.72 g"), "product_role": "15a/15b mixture", "basis": "isolated combined mass"}, {"kind": "selectivity", "targets": ["15a", "15b"], "raw": "2:1 trans/cis", "value": "2:1", "product_role": "15a/15b mixture", "basis": "reported product-mixture ratio", "scope": "15a (trans)/15b (cis)"}])
# Correct the helper placeholder: page 12 gives volume but no mmol for NH4OH.
for p in r15_2["participants"]:
    if p["molecule_id"] == mid("NH4OH"):
        p["charges"]["amount"] = missing(ev(12))


# Add page-7 image evidence to the seven pre-existing cysteamine kinetic records.
kinetic_ids = {
    f"{DOC}:rxn:rxn_p6_cysteamine_{key}" for key in ("45", "49", "51", "52", "54", "60", "S17")
}


def add_image_to_facts(value):
    if isinstance(value, dict):
        if "normalized_value" in value and "evidence_ids" in value and f"{DOC}:page7:text" in value["evidence_ids"]:
            if f"{DOC}:page7:image" not in value["evidence_ids"]:
                value["evidence_ids"].append(f"{DOC}:page7:image")
        for child in value.values():
            add_image_to_facts(child)
    elif isinstance(value, list):
        for child in value:
            add_image_to_facts(child)


for reaction in bundle["reactions"]:
    if reaction["reaction_id"] in kinetic_ids:
        for eid in ev(7):
            if eid not in reaction["evidence_ids"]:
                reaction["evidence_ids"].append(eid)
        add_image_to_facts(reaction)
        reaction["review"]["unresolved"] = [u for u in reaction["review"]["unresolved"] if "page 7 image was outside" not in u]
        if reaction["reaction_id"].endswith(("cysteamine_45", "cysteamine_54")):
            footnote = "Table 3 footnote b: reaction was performed with cysteamine that was approximately 60% oxidized."
            if not any(n.get("normalized_value") == footnote for n in reaction["notes"]):
                reaction["notes"].append(extracted(footnote, footnote, ev(7)))


def source_context(context_id, page, region, classification_name, summary, *, pages=None):
    evidence = ev(*(pages or [page]))
    return {
        "context_id": context_id,
        "document_id": DOC,
        "location": {"page": page, "region": region},
        "classification": classification_name,
        "summary": summary,
        "evidence_ids": evidence,
        "review_status": "needs_review",
    }


new_contexts = [
    source_context("ctx_p7_table3", 7, "Table 3 and footnotes", "mechanistic_table", "Seven current cysteamine kinetic experiments are linked to existing records. PTL entry 8 (<5 min) is explicitly based on previous literature and was not promoted to a current experiment. Footnote b (~60% oxidized cysteamine) applies to 45 and 54; 54 was complete by <10 min."),
    source_context("ctx_p7_bioactivity", 7, "NF-κB bioactivity and cytotoxicity discussion", "biological_assay", "Cellular reporter and viability findings are biological assay context, not synthetic reactions."),
    source_context("ctx_p8_figure4", 8, "Figure 4", "biological_assay", "NF-κB luciferase/SEAP activity and accompanying cell viability across stated concentrations; treated as assay context rather than chemical synthesis."),
    source_context("ctx_p8_table4", 8, "Table 4", "biological_assay", "Vero-cell IC50 values: 47 >200 μM; 49 21.9±1.9 μM; 51 15.7±2.0 μM; 52 44.0±1.8 μM; 60 7.8±1.9 μM; S17 6.2±1.8 μM (mean±SEM, n≥3)."),
    source_context("ctx_p9_fig5", 9, "Figure 5", "characterization", "ALARM NMR spectra provide assay readouts for the separate mechanistic-study records; spectra are not independent syntheses."),
    source_context("ctx_p9_fig6_workflow", 9, "Figure 6A", "proposed_mechanism", "ALARM MSPS analytical workflow showing IAD labelling, denaturation, heavy-IAD labelling, trypsin digestion and MS/MS."),
    source_context("ctx_p10_conclusions", 10, "Conclusions", "discussion", "Interpretive summary of APKR scope, tunable thiol reactivity and NF-κB activity; no new synthetic run."),
    source_context("ctx_p10_general_chemistry", 10, "Experimental Section: Chemistry general methods", "general_method", "General solvent purification, atmosphere, chromatography and instrument methods. These are context and are not duplicated as reactions."),
    source_context("ctx_p11_characterization", 11, "S2/S3/11/S5 characterization blocks", "characterization", "NMR, IR, HRMS, TLC, boiling point and literature-match data are characterization attached to molecules, not separate synthetic experiments."),
    source_context("ctx_p12_12_13_characterization", 12, "compound 12 and 13 preparation mentions", "literature_example", "Compounds 12 and 13 are stated to have been prepared as previously described; page 12 gives spectra but no standalone current quantities, so no current synthesis was invented."),
    source_context("ctx_p12_14_15_characterization", 12, "14a/14b/15a characterization blocks", "characterization", "NMR/IR/HRMS/TLC blocks are characterization, not additional reactions."),
]
existing_context_ids = {c["context_id"] for c in bundle.get("source_context", [])}
bundle.setdefault("source_context", []).extend(c for c in new_contexts if c["context_id"] not in existing_context_ids)


(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")

new_reaction_ids = [r["reaction_id"] for r in bundle["reactions"] if any(eid.startswith(f"{DOC}:page{p}:") for p in range(7, 13) for eid in r["evidence_ids"])]
inventory = {
    "batch_id": "allene-stereo-v31-200-20261017",
    "document_id": DOC,
    "visually_inspected_original_pdf_pages": [7, 8, 9, 10, 11, 12],
    "new_reaction_ids": [r["reaction_id"] for r in bundle["reactions"] if r["reaction_id"].startswith(f"{DOC}:rxn:rxn_p7_") or r["reaction_id"].startswith(f"{DOC}:rxn:rxn_p9_") or r["reaction_id"].startswith(f"{DOC}:rxn:rxn_p11_") or r["reaction_id"].startswith(f"{DOC}:rxn:rxn_p12_")],
    "page_linked_reaction_ids_including_updated_existing": new_reaction_ids,
    "new_molecule_ids": [mid(k) for k in ("S11", "S1", "S2", "S3", "La_antigen", "tBuSH", "DTT", "DMSO", "CPM", "fluconazole", "IAD", "d2_IAD", "trypsin", "S1_propargyl_chloride", "K2CO3", "EtOAc", "MgSO4", "brine", "pTsOH", "benzene", "ethylene_glycol", "nBuLi", "methyl_chloroformate", "NH4Cl", "trimethyl_borate", "bromochloromethane", "TMSCl", "pinacol", "NaHCO3", "HCl", "silica")],
    "source_context_ids": [c["context_id"] for c in new_contexts],
    "limitations": ["All new records remain needs_review and partial.", "No SI association was inferred.", "Complex relative/mixture stereochemistry was retained without speculative isomeric SMILES.", "No three-model vote was performed."],
}
(RUN / "raw_model" / "pages_0007_0012_inventory.json").write_text(json.dumps(inventory, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model" / "pages_0007_0012_standard_stereo_checks.json").write_text(json.dumps({
    "S1": {"input_smiles": "CC(=O)CC(C)=O", "canonical_isomeric_smiles": "CC(=O)CC(C)=O", "stereo_graph_key": "f992fc62b91bfd6ce7dcd013769b89f22e454906c434d96125100f5889e38211", "atom_centers": [], "double_bonds": [{"bond_index": 1, "stereo": "STEREONONE"}, {"bond_index": 5, "stereo": "STEREONONE"}], "source_verified": False},
    "S2": {"input_smiles": "CC(=O)CCC#C", "canonical_isomeric_smiles": "C#CCCC(C)=O", "stereo_graph_key": "4eb7816f0de49357137b3887e19c59f147a9dcd8d5fc6196b9a5b9b00f654ccf", "atom_centers": [], "double_bonds": [{"bond_index": 1, "stereo": "STEREONONE"}], "source_verified": False},
    "S3": {"input_smiles": "CC1(CCC#C)OCCO1", "canonical_isomeric_smiles": "C#CCCC1(C)OCCO1", "stereo_graph_key": "aa20b2d69baea6f1745870d6e5f5a1ad9353b7a37fe4cd491616ff6148329b9a", "atom_centers": [], "double_bonds": [], "source_verified": False},
    "11": {"input_smiles": "COC(=O)C#CCC1(C)OCCO1", "canonical_isomeric_smiles": "COC(=O)C#CCC1(C)OCCO1", "stereo_graph_key": "8576f4d6c3a911aa62a7e6ee89a6f2694b2b904b6cbd786314b5599ac6fbdf9e", "atom_centers": [], "double_bonds": [{"bond_index": 2, "stereo": "STEREONONE"}], "source_verified": False},
}, ensure_ascii=False, indent=2) + "\n")

page_record_ids = {page: [] for page in range(7, 13)}
for reaction in bundle["reactions"]:
    for page in range(7, 13):
        if any(eid in reaction["evidence_ids"] for eid in ev(page)):
            page_record_ids[page].append(reaction["reaction_id"])

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 7,
            "state": "attempted",
            "regions": [
                {"region": "Table 3 including Eq 9 and footnotes", "finding": "seven current cysteamine kinetic experiments linked to existing records; PTL entry 8 is explicitly previous literature", "evidence_ids": ev(7)},
                {"region": "left-column structural confirmation paragraph", "finding": "actual compound 48 plus tert-butyl thiol experiment giving S11, 23% yield and 1.6:1 dr", "evidence_ids": ev(7)},
                {"region": "NF-κB bioactivity discussion", "finding": "cellular reporter/viability assay context, not a synthetic reaction", "evidence_ids": ev(7)},
            ],
            "record_ids": page_record_ids[7],
            "unresolved": ["PTL Table 3 entry is literature-derived and was not promoted to a current experiment.", "S11 exact diastereomer configurations and yield measurement basis are not reported on page 7."],
        },
        {
            "page": 8,
            "state": "attempted",
            "regions": [
                {"region": "Figure 4 and caption", "finding": "NF-κB luciferase/SEAP reporter assays and viability measurements; biological context only", "evidence_ids": ev(8)},
                {"region": "Table 4", "finding": "Vero-cell IC50 values transcribed to source context, not synthetic experiments", "evidence_ids": ev(8)},
                {"region": "ALARM NMR discussion", "finding": "actual La-antigen thiol-reactivity tests for 51 and S17, including 37 °C/90 min and ±20 mM DTT", "evidence_ids": ev(8)},
            ],
            "record_ids": page_record_ids[8],
            "unresolved": ["No discrete La-protein adduct structure is reported for the ALARM NMR experiments."],
        },
        {
            "page": 9,
            "state": "attempted",
            "regions": [
                {"region": "Figure 5 and caption", "finding": "five ALARM NMR assay rows recorded as mechanistic studies/controls; spectra themselves are characterization", "evidence_ids": ev(9)},
                {"region": "Figure 6A workflow", "finding": "ALARM MSPS analytical workflow; generic workflow retained as proposed-mechanism/source context", "evidence_ids": ev(9)},
                {"region": "Figure 6B heat map and footnotes", "finding": "five compound/control rows recorded separately with biological/technical replicate caveats", "evidence_ids": ev(9)},
            ],
            "record_ids": page_record_ids[9],
            "unresolved": ["Figure 6 does not supply complete isolated protein-adduct structures.", "One C245 peptide/replicate absence and low-abundance caveats are preserved rather than imputed."],
        },
        {
            "page": 10,
            "state": "attempted",
            "regions": [
                {"region": "ALARM MSPS interpretation", "finding": "S17/CPM reversible-adduction patterns and 51/fluconazole negative patterns linked to Figure 6 records", "evidence_ids": ev(10)},
                {"region": "Conclusions", "finding": "interpretive summary; no additional reaction record", "evidence_ids": ev(10)},
                {"region": "Experimental Section, Chemistry general methods", "finding": "general solvent, atmosphere, chromatography and instrument methods; context rather than a specific run", "evidence_ids": ev(10)},
            ],
            "record_ids": page_record_ids[10],
            "unresolved": ["ALARM MSPS differentiates labelling patterns but does not provide an isolated covalent product structure."],
        },
        {
            "page": 11,
            "state": "attempted",
            "regions": [
                {"region": "Hex-5-yn-2-one (S2) procedure", "finding": "specific literature-analogous preparative run transcribed with ordered reaction/workup/distillation steps", "evidence_ids": ev(11)},
                {"region": "S3, compound 11 and S5 procedures", "finding": "three separate specific preparations transcribed with quantities, temperature stages, workups and isolated yields", "evidence_ids": ev(11)},
                {"region": "compound 9 Run 1 through S5 addition", "finding": "first half of a cross-page specific run; merged with page-12 workup/outcome rather than duplicated", "evidence_ids": ev(11)},
                {"region": "spectral/TLC paragraphs", "finding": "characterization retained as source context, not separate reactions", "evidence_ids": ev(11)},
            ],
            "record_ids": page_record_ids[11],
            "unresolved": ["Named-agent structures/CAS were not inferred where not necessary to the reaction identity.", "Atom-mapped products and full bond-change inventories remain unavailable."],
        },
        {
            "page": 12,
            "state": "attempted",
            "regions": [
                {"region": "compound 9 Run 1 continuation and Run 2", "finding": "cross-page Run 1 completed; Run 2 recorded separately with its own quantities, yield, Z/E ratio and S4 contamination", "evidence_ids": ev(12)},
                {"region": "compound 12/13 preparation mentions", "finding": "previous-literature preparations plus characterization; not promoted to current experiments without standalone current run details", "evidence_ids": ev(12)},
                {"region": "General Procedure A, 14a/14b Runs 1-3", "finding": "three distinct current runs recorded; Run 3 chromatographic fractions retain separate yields, masses and ratio order", "evidence_ids": ev(12)},
                {"region": "15a/15b Runs 1-2", "finding": "two distinct current runs recorded; Run 1 fraction yields/contamination/overall ratio retained", "evidence_ids": ev(12)},
                {"region": "14a/14b/15a characterization", "finding": "spectral/analytical data treated as characterization, not additional syntheses", "evidence_ids": ev(12)},
            ],
            "record_ids": page_record_ids[12],
            "unresolved": ["The eleven-run aggregate range for 14a/14b was not expanded into unreported individual runs.", "Complex lactam relative stereochemistry remains partial without stereo-complete isomeric SMILES.", "No SI association was inferred."],
        },
    ],
}
(RUN / "coverage" / "pages_0007_0012.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")

print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]), "source_context": len(bundle["source_context"]), "new_records": len(inventory["new_reaction_ids"])}, ensure_ascii=False))
