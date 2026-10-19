#!/usr/bin/env python3
"""Extend the v3.1 source bundle after visual review of PDF pages 7–12.

The previous pages 1–6 bundle is treated only as the current workspace state.
This script does not recreate those records.  It adds page-specific evidence,
characterization amounts/names, and the one newly encountered synthesis of 6.
All machine-produced records remain needs_review.
"""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st


DOC = "doc_1d33e0721b132f96"
EVT = {page: f"{DOC}:page{page}:text" for page in range(1, 13)}
EVI = {page: f"{DOC}:page{page}:image" for page in range(1, 13)}


def ev(*pages: int) -> list[str]:
    values: list[str] = []
    for page in pages:
        values.extend((EVT[page], EVI[page]))
    return list(dict.fromkeys(values))


def extend_unique(values: list, additions: list) -> None:
    for value in additions:
        if value not in values:
            values.append(value)


def add_fact_evidence(fact: dict, evidence: list[str]) -> None:
    extend_unique(fact["evidence_ids"], evidence)


def add_page_image_where_text(value, page: int) -> None:
    """Attach an inspected image only to facts already sourced to that page text."""
    if isinstance(value, dict):
        evidence = value.get("evidence_ids")
        if isinstance(evidence, list) and EVT[page] in evidence:
            extend_unique(evidence, [EVI[page]])
        for child in value.values():
            add_page_image_where_text(child, page)
    elif isinstance(value, list):
        for child in value:
            add_page_image_where_text(child, page)


def src(raw, normalized=None, evidence=()):
    return st.fact(raw, raw if normalized is None else normalized, evidence)


def inf(normalized, evidence, derivation):
    return st.fact(None, normalized, evidence, method="inference", derivation=derivation)


def derived(normalized, evidence, derivation, raw):
    return st.fact(raw, normalized, evidence, method="structure_derivation", derivation=derivation)


def quantity(number, unit, raw_text):
    return {"number": number, "unit": unit, "raw_text": raw_text}


def participant(pid, mid, role, evidence, stage="reaction", product_role=None,
                mass=None, amount=None, volume=None, molpct=None, ref=None):
    item = st.new_participant(pid, mid, role, evidence, stage)
    if product_role is not None:
        item["product_role"] = src(product_role, product_role, evidence)
    for field, value in {
        "mass": mass,
        "amount": amount,
        "volume": volume,
        "mol_percent": molpct,
    }.items():
        if value is not None:
            item["charges"][field] = src(value["raw_text"], value, evidence)
    if ref is not None:
        item["equivalence_reference"] = src(ref, ref, evidence)
    return item


def step(step_id, stage, participant_ids, evidence, operation, temperature=None,
         time=None, atmosphere=None):
    item = st.new_step(step_id, evidence, stage)
    item["participant_ids"] = participant_ids
    item["operation"] = src(operation, operation, evidence)
    if temperature is not None:
        raw = temperature.get("raw_text", temperature) if isinstance(temperature, dict) else temperature
        item["temperature"] = src(raw, temperature, evidence)
    if time is not None:
        raw = time.get("raw_text", time) if isinstance(time, dict) else time
        item["time"] = src(raw, time, evidence)
    if atmosphere is not None:
        item["atmosphere"] = src(atmosphere, atmosphere, evidence)
    return item


def outcome(oid, kind, targets, evidence, value, basis, role="main product", yield_type=None):
    item = st.new_outcome(oid, kind, targets, evidence)
    raw = value.get("raw_text", value) if isinstance(value, dict) else value
    item["value"] = src(raw, value, evidence)
    item["product_role"] = src(role, role, evidence)
    item["measurement_basis"] = src(basis, basis, evidence)
    if yield_type is not None:
        item["yield_type"] = src(yield_type, yield_type, evidence)
    return item


prepared_dataset = json.loads((RUN / "dataset.json").read_text())
bundle = copy.deepcopy(json.loads((RUN / "standard_records.json").read_text()))
# The previous merge intentionally contained only the images then claimed as
# visually covered.  Reuse the prepared dataset's immutable document/evidence
# block now that pages 7–12 have also been inspected.
bundle["documents"] = copy.deepcopy(prepared_dataset["documents"])
# Make reruns deterministic even after this continuation has already been
# merged once into standard_records.json.
bundle["reactions"] = [
    item for item in bundle["reactions"]
    if item["reaction_id"] != f"{DOC}:preparation_substrate_6"
]
for item in bundle["reactions"]:
    item["outcomes"] = [
        outcome_item for outcome_item in item["outcomes"]
        if outcome_item["outcome_id"] != "product_amount_characterization"
    ]
bundle["molecules"] = [
    item for item in bundle["molecules"]
    if item["molecule_id"] != f"{DOC}:mol:methanol"
]
bundle["source_context"] = [
    item for item in bundle["source_context"]
    if item["context_id"] not in {
        f"{DOC}:pages7_9_characterization",
        f"{DOC}:pages10_12_references",
    }
]
reaction_by_id = {item["reaction_id"]: item for item in bundle["reactions"]}
molecule_by_id = {item["molecule_id"]: item for item in bundle["molecules"]}


# Characterization blocks repeat the isolated yields and newly provide masses.
# They are merged into the existing scope experiment records, never counted as
# separate synthesis reactions.
characterization = {
    "3e": (37.6, (7,)),
    "3f": (66.9, (7,)),
    "3g": (70.5, (7,)),
    "3h": (54.8, (7,)),
    "3i": (62.8, (7,)),
    "3j": (28.7, (7,)),
    "3k": (56.6, (7, 8)),
    "3l": (53.9, (8,)),
    "3m": (61.9, (8,)),
    "3n": (53.0, (8,)),
    "3o": (43.3, (8,)),
    "3p": (57.9, (8,)),
    "3q": (21.7, (8,)),
    "3r": (29.7, (8,)),
    "4a": (67.6, (8,)),
    "4b": (76.4, (8,)),
    "4c": (61.0, (8,)),
    "4d": (56.6, (8, 9)),
    "4e": (53.9, (9,)),
    "4f": (66.8, (9,)),
    "4g": (64.3, (9,)),
    "4h": (37.0, (9,)),
    "4i": (19.2, (9,)),
    "4j": (23.8, (9,)),
    "4k": (28.1, (9,)),
}

product_names = {
    "3e": "4-(diphenylphosphoryl)-5-(4-fluorophenyl)-3-methylfuran-2(5H)-one",
    "3f": "5-(4-chlorophenyl)-4-(diphenylphosphoryl)-3-methylfuran-2(5H)-one",
    "3g": "5-(4-bromophenyl)-4-(diphenylphosphoryl)-3-methylfuran-2(5H)-one",
    "3h": "4-(diphenylphosphoryl)-3-methyl-5-(4-(trifluoromethyl)phenyl)furan-2(5H)-one",
    "3i": "5-(3-chlorophenyl)-4-(diphenylphosphoryl)-3-methylfuran-2(5H)-one",
    "3j": "4-(diphenylphosphoryl)-3-methyl-5-(o-tolyl)furan-2(5H)-one",
    "3k": "4-(diphenylphosphoryl)-5-(2-methoxyphenyl)-3-methylfuran-2(5H)-one",
    "3l": "5-(2-chlorophenyl)-4-(diphenylphosphoryl)-3-methylfuran-2(5H)-one",
    "3m": "5-(2,5-dimethylphenyl)-4-(diphenylphosphoryl)-3-methylfuran-2(5H)-one",
    "3n": "5-(3,5-bis(trifluoromethyl)phenyl)-4-(diphenylphosphoryl)-3-methylfuran-2(5H)-one",
    "3o": "4-(diphenylphosphoryl)-3-methyl-5-(naphthalen-2-yl)furan-2(5H)-one",
    "3p": "4-(diphenylphosphoryl)-5-phenyl-3-propylfuran-2(5H)-one",
    "3q": "4-(diphenylphosphoryl)-3,5-dimethyl-5-phenylfuran-2(5H)-one",
    "3r": "4-(diphenylphosphoryl)-3-methyl-5,5-diphenylfuran-2(5H)-one",
    "4a": "4-(di-p-tolylphosphoryl)-3-methyl-5-phenylfuran-2(5H)-one",
    "4b": "4-(bis(4-methoxyphenyl)phosphoryl)-3-methyl-5-phenylfuran-2(5H)-one",
    "4c": "4-(di([1,1'-biphenyl]-4-yl)phosphoryl)-3-methyl-5-phenylfuran-2(5H)-one",
    "4d": "4-(bis(4-fluorophenyl)phosphoryl)-3-methyl-5-phenylfuran-2(5H)-one",
    "4e": "4-(bis(4-chlorophenyl)phosphoryl)-3-methyl-5-phenylfuran-2(5H)-one",
    "4f": "4-(bis(4-bromophenyl)phosphoryl)-3-methyl-5-phenylfuran-2(5H)-one",
    "4g": "4-(di-o-tolylphosphoryl)-3-methyl-5-phenylfuran-2(5H)-one",
    "4h": "4-(bis(3,5-dimethylphenyl)phosphoryl)-3-methyl-5-phenylfuran-2(5H)-one",
    "4i": "Diethyl (4-methyl-5-oxo-2-phenyl-2,5-dihydrofuran-3-yl)phosphonate",
    "4j": "Diisopropyl (4-methyl-5-oxo-2-phenyl-2,5-dihydrofuran-3-yl)phosphonate",
    "4k": "Ethyl (4-methyl-5-oxo-2-phenyl-2,5-dihydrofuran-3-yl)(phenyl)phosphinate",
}

for label, (mass, pages) in characterization.items():
    evidence = ev(*pages)
    prefix = "scheme2_" if label.startswith("3") else "scheme3_"
    reaction = reaction_by_id[f"{DOC}:{prefix}{label}"]
    extend_unique(reaction["evidence_ids"], evidence)
    add_fact_evidence(reaction["location"], evidence)
    for existing in reaction["outcomes"]:
        if existing["kind"] == "yield":
            for field in ("product_role", "value", "yield_type", "measurement_basis", "scope"):
                add_fact_evidence(existing[field], evidence)
    product_id = f"{DOC}:mol:{label}"
    reaction["outcomes"].append(outcome(
        "product_amount_characterization",
        "product_amount",
        [product_id],
        evidence,
        quantity(mass, "mg", f"{mass:.1f} mg"),
        "reported isolated product amount in characterization block",
    ))

    molecule = molecule_by_id[product_id]
    extend_unique(molecule["occurrences"][0]["evidence_ids"], evidence)
    add_fact_evidence(molecule["occurrences"][0]["compound_label"], evidence)
    molecule["names"]["original"] = src(product_names[label], product_names[label], evidence)
    molecule["names"]["english"] = src(product_names[label], product_names[label], evidence)


# Scheme 6 is a proposed pathway.  The drawing itself is now visually covered,
# but it remains one mechanism record rather than a set of isolated experiments.
mechanism = reaction_by_id[f"{DOC}:scheme6_proposed_radical_pathway"]
extend_unique(mechanism["evidence_ids"], ev(7))
mechanism["location"] = src(
    "PDF page 7, Scheme 6; mechanism discussion on PDF page 3",
    "PDF page 7, Scheme 6; mechanism discussion on PDF page 3",
    ev(3, 7),
)
for field in mechanism["classification"]["bond_changes"]:
    add_fact_evidence(field, ev(7))
for field_name in ("mapped_reaction_smiles", "labels", "coverage"):
    add_fact_evidence(mechanism["classification"][field_name], ev(7))
add_fact_evidence(mechanism["procedure"]["original_text"], ev(7))
mechanism["review"]["unresolved"] = [
    "This is an author-proposed mechanism, not a set of isolated intermediates; intermediate A/B structures and both radical-generation paths remain proposals."
]


# The page-9 original image confirms that the later prose really prints
# t-BuOK 0.4 mL/overnight.  Preserve it alongside the contradictory page-5
# Scheme 4 values; do not choose or repair either source statement.
ring_opening = reaction_by_id[f"{DOC}:scheme4_ring_opening_3a_to_5"]
add_page_image_where_text(ring_opening, 9)
ring_opening["location"] = src(
    "PDF page 5 Scheme 4 and PDF page 9 experimental procedure",
    "PDF page 5 Scheme 4 and PDF page 9 experimental procedure",
    ev(5, 9),
)
ring_opening["review"]["unresolved"] = [
    "Source conflict confirmed by both original images: Scheme 4 reports t-BuOH 4 mL and 12 h, whereas the page-9 procedure prints t-BuOK 0.4 mL and overnight; no value was selected or corrected.",
    "Complete bond changes for ring opening were not atom mapped.",
]
for molecule_id in (f"{DOC}:mol:5", f"{DOC}:mol:tbuok"):
    molecule = molecule_by_id[molecule_id]
    add_page_image_where_text(molecule, 9)
    extend_unique(molecule["occurrences"][0]["evidence_ids"], ev(9))


# Methanol is newly mentioned as the reaction solvent for preparing 6.
evidence9 = ev(9)
methanol_id = f"{DOC}:mol:methanol"
methanol_diag = st.stereo_inspect("CO")
methanol = st.new_molecule(methanol_id, DOC, evidence9)
methanol["identity_status"] = "source_identified"
methanol["representation_type"] = "smiles"
methanol["occurrences"][0]["compound_label"] = src("methanol", "methanol", evidence9)
methanol["names"]["original"] = src("methanol", "methanol", evidence9)
methanol["names"]["english"] = src("methanol", "methanol", evidence9)
methanol["smiles"] = derived(
    methanol_diag["canonical_isomeric_smiles"], evidence9,
    "Connectivity derived from the source chemical name and canonicalized by RDKit; parsing is not source verification.",
    "methanol",
)
methanol["stereochemistry"]["assignment"] = inf(
    "no stereogenic center/axis identified in the machine connectivity candidate",
    evidence9,
    "Graph inspection of the methanol candidate; absence of stereochemical labels alone was not used.",
)
methanol["stereochemistry"]["chirality_status"] = inf(
    "achiral", evidence9, "Standard stereo inspection of the methanol connectivity candidate found no stereogenic element."
)
methanol["stereochemistry"]["geometry_status"] = inf(
    "not_applicable", evidence9, "The methanol connectivity candidate has no E/Z or ring cis/trans element."
)
methanol["stereochemistry"]["representation_status"] = inf(
    "encoded", evidence9, "The achiral machine connectivity candidate is fully represented by the isomeric SMILES."
)
methanol["stereochemistry"]["isomeric_smiles"] = derived(
    methanol_diag["canonical_isomeric_smiles"], evidence9,
    "Canonical isomeric SMILES from the same machine connectivity candidate.",
    "methanol",
)
methanol["stereochemistry"]["identity_key"] = inf(
    methanol_diag["stereo_graph_key"], evidence9,
    "SHA-256 stereo graph key calculated from the canonical isomeric SMILES; not a sample identity or source verification.",
)
methanol["identity_note"] = src(
    "Machine connectivity candidate from the source name; no CAS number was reported.",
    "Machine connectivity candidate from the source name; no CAS number was reported.",
    evidence9,
)
bundle["molecules"].append(methanol)
molecule_by_id[methanol_id] = methanol


# Add page-9 evidence to molecules directly participating in the newly found
# substrate-preparation procedure.  Their pre-existing machine structures and
# stereo statuses are retained unchanged.
for molecule_id in (
    f"{DOC}:mol:1a", f"{DOC}:mol:cucl", f"{DOC}:mol:6",
    f"{DOC}:mol:pe", f"{DOC}:mol:ea",
):
    molecule = molecule_by_id[molecule_id]
    extend_unique(molecule["occurrences"][0]["evidence_ids"], evidence9)
    add_fact_evidence(molecule["occurrences"][0]["compound_label"], evidence9)
    add_fact_evidence(molecule["names"]["original"], evidence9)
    add_fact_evidence(molecule["names"]["english"], evidence9)
    if molecule_id in (f"{DOC}:mol:1a", f"{DOC}:mol:6"):
        for stereo_field in ("atom_centers", "double_bonds", "axial", "mixture", "assignment",
                             "chirality_status", "geometry_status", "representation_status",
                             "isomeric_smiles", "identity_key"):
            add_fact_evidence(molecule["stereochemistry"][stereo_field], evidence9)
        add_fact_evidence(molecule["smiles"], evidence9)


# New specific experiment: preparation of substrate 6 from allenoic acid 1a.
preparation_id = f"{DOC}:preparation_substrate_6"
preparation = st.new_reaction(preparation_id, DOC, evidence9)
preparation["record_type"] = "specific_experiment"
preparation["experiment_categories"] = src(["other"], ["other"], evidence9)
preparation["location"] = src(
    "PDF page 9, General Procedure for the Synthesis of Substrates 6",
    "PDF page 9, General Procedure for the Synthesis of Substrates 6",
    evidence9,
)
preparation["reaction_name"] = src(
    "CuCl-catalyzed cycloisomerization of 1a to substrate 6",
    "CuCl-catalyzed cycloisomerization of 1a to substrate 6",
    evidence9,
)
preparation["transformation_type"] = inf(
    "allenoic-acid cycloisomerization to a butenolide",
    evidence9,
    "Machine description of the source procedure and substrate/product identities.",
)
preparation["participants"] = [
    participant(
        "substrate_1a", f"{DOC}:mol:1a", "reactant", evidence9,
        mass=quantity(380, "mg", "380 mg"), amount=quantity(2.0, "mmol", "2.0 mmol"),
    ),
    participant(
        "cucl", f"{DOC}:mol:cucl", "catalyst", evidence9,
        amount=quantity(0.08, "mmol", "0.08 mmol"),
        molpct=quantity(4, "mol%", "4 mol %"), ref=f"{DOC}:mol:1a",
    ),
    participant(
        "methanol", methanol_id, "solvent", evidence9,
        volume=quantity(20, "mL", "20 mL"),
    ),
    participant(
        "product_6", f"{DOC}:mol:6", "product", evidence9,
        product_role="isolated substrate 6",
    ),
    participant("pe", f"{DOC}:mol:pe", "purification_agent", evidence9, stage="purification"),
    participant("ea", f"{DOC}:mol:ea", "purification_agent", evidence9, stage="purification"),
]
preparation["reaction_scale"] = src(
    "1a: 2.0 mmol",
    {
        "reference_molecule_id": f"{DOC}:mol:1a",
        "amount": {"number": 2.0, "unit": "mmol", "raw_text": "2.0 mmol"},
    },
    evidence9,
)
preparation["solvent_conditions"]["composition"] = src("methanol", ["methanol"], evidence9)
preparation["solvent_conditions"]["total_volume"] = src(
    "20 mL", quantity(20, "mL", "20 mL"), evidence9
)
preparation["steps"] = [
    step(
        "charge_under_n2", "reaction", ["substrate_1a", "cucl", "methanol"], evidence9,
        "prepare a solution of 1a and CuCl in methanol under N2", atmosphere="N2",
    ),
    step(
        "reflux", "reaction", ["substrate_1a", "cucl", "methanol"], evidence9,
        "stir under reflux", temperature="reflux (methanol; numerical temperature unreported)",
        time=quantity(2, "h", "2 h"), atmosphere="N2",
    ),
    step(
        "evaporate", "workup", ["product_6"], evidence9,
        "evaporate the solvent",
    ),
    step(
        "flash_chromatography", "purification", ["product_6", "pe", "ea"], evidence9,
        "flash chromatography on silica gel with PE/EA = 10:1",
    ),
]
preparation["procedure"]["original_text"] = src(
    "Under N2, 1a (380 mg, 2.0 mmol) and CuCl (0.08 mmol, 4 mol %) in methanol (20 mL) were stirred under reflux for 2 h; after evaporation, flash chromatography on silica gel (PE/EA 10:1) gave 6 (226 mg, 65%) as a colorless liquid.",
    "source procedure retained with ordered reaction, evaporation, and purification operations",
    evidence9,
)
preparation["procedure"]["addition_order"] = src(
    "1a and CuCl were dissolved in methanol under N2 before reflux",
    ["1a", "CuCl", "methanol", "reflux"],
    evidence9,
)
preparation["outcomes"] = [
    outcome(
        "isolated_yield", "yield", [f"{DOC}:mol:6"], evidence9,
        quantity(65, "%", "65%"), "isolated after flash chromatography", yield_type="isolated",
    ),
    outcome(
        "isolated_product_amount", "product_amount", [f"{DOC}:mol:6"], evidence9,
        quantity(226, "mg", "226 mg"), "isolated after flash chromatography",
    ),
]
preparation["classification"]["mapped_reaction_smiles"] = st.fact(evidence=evidence9)
preparation["classification"]["bond_changes"] = [
    inf(
        {
            "change_id": "b1",
            "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}],
            "before_order": 0,
            "after_order": 1,
            "change_type": "formed",
        },
        evidence9,
        "New intramolecular C–O connection inferred from the source identities and the retained structure candidates; atom indices and complete mapping are unresolved.",
    )
]
preparation["classification"]["labels"] = inf(
    [{"category_number": 3, "code": "C-O", "supporting_change_ids": ["b1"]}],
    evidence9,
    "C–O label generated from the partial inferred net bond change.",
)
preparation["classification"]["coverage"] = inf(
    "partial",
    evidence9,
    "The visually supported ring-closing C–O bond is represented; complete atom mapping and all bond-order/proton-transfer changes were not established.",
)
preparation["classification"]["limitations"] = [
    "No atom-mapped reaction SMILES was asserted.",
    "Allene-to-alkene bond-order changes and proton transfer were not fully mapped.",
]
preparation["classification"]["model_votes"] = []
preparation["notes"] = [
    src(
        "No R/S, enantiomeric ratio, diastereomeric ratio, or other stereoselectivity is reported for product 6.",
        "Product 6 configuration remains unspecified; no racemate or single stereoisomer was inferred.",
        evidence9,
    )
]
preparation["review"] = {
    "status": "needs_review",
    "method": "model",
    "reviewer": None,
    "checked_fields": [],
    "unresolved": [
        "The numerical reflux temperature is unreported.",
        "The tetrahedral configuration of product 6 is unreported.",
        "Complete atom mapping and all net bond-order/proton changes remain unresolved.",
    ],
}
preparation["completeness"] = "partial"
bundle["reactions"].append(preparation)
reaction_by_id[preparation_id] = preparation


# Source-context coverage: characterization/spectra and references are retained
# as context but never converted into synthesis records.
source_context_by_id = {item["context_id"]: item for item in bundle["source_context"]}
si_context = source_context_by_id[f"{DOC}:si_notice"]
extend_unique(si_context["evidence_ids"], evidence9)
si_context["note"] = (
    "Page 9 states that Supporting Information is available online, but no SI file was explicitly provided or associated in this workspace."
)
page9_context = source_context_by_id[f"{DOC}:page9_ring_opening_context"]
page9_context["kind"] = "procedure_cross_reference"
page9_context["evidence_ids"] = evidence9
page9_context["note"] = (
    "The page-9 original procedure is direct evidence for the existing ring-opening record and confirms, rather than resolves, its conflict with Scheme 4."
)
bundle["source_context"].append({
    "context_id": f"{DOC}:pages7_9_characterization",
    "document_id": DOC,
    "kind": "characterization",
    "evidence_ids": ev(7, 8, 9),
    "note": (
        "Compound names, isolated masses/yields, melting points, NMR and HRMS blocks for 3e–4k; names/masses support existing reactions, while spectral peaks are characterization rather than synthesis experiments."
    ),
})
bundle["source_context"].append({
    "context_id": f"{DOC}:pages10_12_references",
    "document_id": DOC,
    "kind": "bibliography",
    "evidence_ids": ev(10, 11, 12),
    "note": (
        "Author information, acknowledgments, and cited references; no current-paper reaction procedures or result rows were identified on these pages."
    ),
})


# Raw extraction inventory and exact requested page-coverage file.
rid = lambda suffix: f"{DOC}:{suffix}"
page7_scope = [rid(f"scheme2_3{label}") for label in "efghijk"]
page8_scope = [rid(f"scheme2_3{label}") for label in "klmnopqr"] + [rid(f"scheme3_4{label}") for label in "abcd"]
page9_scope = [rid(f"scheme3_4{label}") for label in "defghijk"]
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 7,
            "state": "attempted",
            "regions": [
                {
                    "region": "upper half, Scheme 6 Proposed Reaction Mechanism",
                    "finding": "Author-proposed P-centered-radical generation (Paths A/B), addition, oxidation, and cyclization; proposed mechanism only, not isolated intermediate experiments.",
                    "evidence_ids": ev(7),
                },
                {
                    "region": "lower columns, characterization of 3e–3k",
                    "finding": "Names, isolated masses/yields, melting points and spectra supplement existing Scheme 2 records; spectral peak lists are not synthesis experiments.",
                    "evidence_ids": ev(7),
                },
            ],
            "record_ids": [rid("scheme6_proposed_radical_pathway")] + page7_scope,
            "unresolved": [
                "Mechanistic intermediates A/B and the alternative radical-generation paths are author proposals, not isolated species.",
                "No absolute/relative product configurations or stereoselectivity metrics are reported in the characterization blocks.",
            ],
        },
        {
            "page": 8,
            "state": "attempted",
            "regions": [
                {
                    "region": "full two-column page, characterization continuation for 3k–3r and 4a–4d",
                    "finding": "Compound names, isolated masses/yields and spectral data were merged into existing scope reactions; no independent reaction arrow or new procedure appears.",
                    "evidence_ids": ev(8),
                }
            ],
            "record_ids": page8_scope,
            "unresolved": [
                "Characterization does not assign the tetrahedral configurations of the products; no racemate or single stereoisomer was inferred.",
                "NMR/HRMS values are retained as source context rather than synthetic outcomes.",
            ],
        },
        {
            "page": 9,
            "state": "attempted",
            "regions": [
                {
                    "region": "upper and middle columns, characterization of 4d–4k",
                    "finding": "Names and isolated masses/yields supplement existing Scheme 3 records; spectra are characterization, not separate syntheses.",
                    "evidence_ids": ev(9),
                },
                {
                    "region": "lower left, Open Ring of 3a procedure",
                    "finding": "Direct procedure for the existing derivatization record; the image confirms printed t-BuOK 0.4 mL/overnight, conflicting with Scheme 4 t-BuOH 4 mL/12 h.",
                    "evidence_ids": ev(9),
                },
                {
                    "region": "lower left to right, General Procedure for the Synthesis of Substrates 6",
                    "finding": "Distinct CuCl-catalyzed cycloisomerization of 1a in methanol under N2/reflux for 2 h, followed by evaporation and PE/EA 10:1 chromatography, affording 6 (226 mg, 65%).",
                    "evidence_ids": ev(9),
                },
                {
                    "region": "lower right, Associated Content",
                    "finding": "The article states that SI is available online, but no SI file is provided or associated in this workspace.",
                    "evidence_ids": ev(9),
                },
            ],
            "record_ids": page9_scope + [rid("scheme4_ring_opening_3a_to_5"), preparation_id],
            "unresolved": [
                "The t-BuOH/t-BuOK and 4/0.4 mL plus 12 h/overnight conflicts are preserved without choosing a value.",
                "The exact reflux temperature and product-6 configuration are unreported.",
                "Online SI was mentioned but not supplied; no SI association was inferred.",
            ],
        },
        {
            "page": 10,
            "state": "attempted",
            "regions": [
                {
                    "region": "author information, acknowledgments, and references 1–10",
                    "finding": "Administrative and bibliographic content only; cited works are not experiments of the current paper.",
                    "evidence_ids": ev(10),
                }
            ],
            "record_ids": [],
            "unresolved": ["Cited primary experiments were not supplied and were not reconstructed from bibliography entries."],
        },
        {
            "page": 11,
            "state": "attempted",
            "regions": [
                {
                    "region": "references 10–16",
                    "finding": "Bibliographic references only; no current-paper reaction procedure, table row, or reaction arrow.",
                    "evidence_ids": ev(11),
                }
            ],
            "record_ids": [],
            "unresolved": ["Cited primary experiments were not supplied and were not reconstructed from reference text."],
        },
        {
            "page": 12,
            "state": "attempted",
            "regions": [
                {
                    "region": "references 16–20 and remaining blank page area",
                    "finding": "Bibliographic references only; no current-paper experiment or characterization block.",
                    "evidence_ids": ev(12),
                }
            ],
            "record_ids": [],
            "unresolved": ["Cited primary experiments were not supplied and were not reconstructed from reference text."],
        },
    ],
}

inventory = {
    "document_id": DOC,
    "batch_id": "allene-stereo-v31-200-20261017",
    "claimed_visual_pages": [7, 8, 9, 10, 11, 12],
    "new_reaction_ids": [preparation_id],
    "updated_reaction_ids": sorted(
        [rid(f"scheme2_{label}") for label in characterization if label.startswith("3")]
        + [rid(f"scheme3_{label}") for label in characterization if label.startswith("4")]
        + [rid("scheme4_ring_opening_3a_to_5"), rid("scheme6_proposed_radical_pathway")]
    ),
    "characterization_product_amounts_mg": {
        label: mass for label, (mass, _) in characterization.items()
    },
    "nonreaction_pages": [10, 11, 12],
    "decisions": [
        "Characterization blocks were merged into existing reactions and not counted as experiments.",
        "Scheme 6 remained a proposed-mechanism record.",
        "The page-9 preparation of 6 was added as one distinct specific experiment.",
        "Page-5/page-9 ring-opening conflicts were preserved without resolution.",
        "No Supporting Information association was inferred.",
    ],
}

(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
(RUN / "coverage/pages_0007_0012.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/pages_0007_0012_inventory.json").write_text(json.dumps(inventory, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/pages_0007_0012_stereo_diagnostics.json").write_text(json.dumps({
    "methanol": methanol_diag,
    "1a": st.stereo_inspect("CC(=C=Cc1ccccc1)C(=O)O"),
    "6": st.stereo_inspect("CC1=CC(c2ccccc2)OC1=O"),
}, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/pages_0007_0012_mentions.json").write_text(json.dumps(
    st.source_mentions(
        prepared_dataset,
        ["General Procedure for the Synthesis of Substrates 6", "3-methyl-5-phenyl-2(5H)-furanone", "1a"],
        context=240,
    ),
    ensure_ascii=False,
    indent=2,
) + "\n")

tasks_path = RUN / "tasks.json"
tasks = json.loads(tasks_path.read_text())
for task in tasks:
    if task.get("document_id") == DOC and task.get("evidence_id") in {EVT[page] for page in range(7, 13)}:
        task["state"] = "attempted"
        task["reason"] = "Original page image and extracted text inspected; detailed region inventory is in coverage/pages_0007_0012.json."
tasks_path.write_text(json.dumps(tasks, ensure_ascii=False, indent=2) + "\n")

print(json.dumps({
    "total_reactions": len(bundle["reactions"]),
    "total_molecules": len(bundle["molecules"]),
    "new_reactions": 1,
    "updated_scope_reactions": len(characterization),
}, ensure_ascii=False))
