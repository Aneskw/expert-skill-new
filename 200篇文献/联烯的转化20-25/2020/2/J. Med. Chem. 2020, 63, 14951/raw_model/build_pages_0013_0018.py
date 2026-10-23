#!/usr/bin/env python3
"""Manual v3.1 transcription for original PDF pages 13-18.

The script is intentionally deterministic and idempotent.  It augments the
existing submission, preserves earlier-page records, and leaves all machine
records as needs_review/partial.
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

from standard_records import fact, new_molecule, new_outcome, new_participant, new_reaction, new_step  # noqa: E402


DOC = "doc_832eb8cea3728513"
TITLE = (
    "Synthesis of Guaianolide Analogues with a Tunable α-Methylene−γ-lactam "
    "Electrophile and Correlating Bioactivity with Thiol Reactivity"
)


def ev(*pages: int) -> list[str]:
    out = []
    for page in pages:
        out.extend([f"{DOC}:page{page}:text", f"{DOC}:page{page}:image"])
    return list(dict.fromkeys(out))


def extracted(raw, normalized=None, evidence=()):
    return fact(raw=raw, normalized=raw if normalized is None else normalized, evidence=evidence, method="original_extraction")


def inferred(raw, normalized, evidence, derivation):
    return fact(raw=raw, normalized=normalized, evidence=evidence, method="inference", verification="inferred", derivation=derivation)


def missing(evidence, status="unreported"):
    return fact(evidence=evidence, status=status, verification="missing")


def not_applicable(evidence):
    return fact(evidence=evidence, method="not_applicable", status="not_applicable", verification="missing")


def q(number, unit, *, qualifier=None, raw_text=None):
    value = {"number": number, "unit": unit}
    if qualifier is not None:
        value["qualifier"] = qualifier
    if raw_text is not None:
        value["raw_text"] = raw_text
    return value


def Q(raw, number, unit, qualifier=None):
    return (raw, q(number, unit, qualifier=qualifier, raw_text=raw))


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

# Replace only this batch's deterministic output on reruns.
prefix = f"{DOC}:rxn:rxn_p13_"
batch_prefixes = tuple(f"{DOC}:rxn:rxn_p{p}_" for p in range(13, 19))
bundle["reactions"] = [r for r in bundle["reactions"] if not r["reaction_id"].startswith(batch_prefixes)]
batch_context_prefixes = tuple(f"ctx_p{p}_" for p in range(13, 19))
bundle["source_context"] = [c for c in bundle.get("source_context", []) if not c["context_id"].startswith(batch_context_prefixes)]

molecules = {m["molecule_id"].split(":mol:", 1)[1]: m for m in bundle["molecules"]}


def add_molecule(key, label, page, *, chirality="unresolved", geometry="unresolved", representation="not_encoded", assignment=None, mixture=None, note=None):
    if key in molecules:
        return molecules[key]["molecule_id"]
    evidence = ev(page)
    m = new_molecule(f"{DOC}:mol:{key}", DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = extracted(label, label, evidence)
    m["names"]["original"] = extracted(label, label, evidence)
    m["names"]["english"] = extracted(label, label, evidence)
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
    basis = "Conservative manual review of the source label/drawing; unspecified stereochemistry was not treated as achiral or racemic."
    st["chirality_status"] = inferred(label, chirality, evidence, basis)
    st["geometry_status"] = inferred(label, geometry, evidence, basis)
    st["representation_status"] = inferred(label, representation, evidence, basis)
    st["isomeric_smiles"] = missing(evidence)
    st["identity_key"] = missing(evidence)
    m["identity_note"] = inferred(label, note or "Source-identified entity retained without speculative full-structure encoding.", evidence, basis)
    m["external_properties"] = []
    molecules[key] = m
    bundle["molecules"].append(m)
    return m["molecule_id"]


def mid(key):
    return molecules[key]["molecule_id"]


def add_occurrence(key, label, pages):
    evidence = ev(*pages)
    m = molecules[key]
    if not any(o.get("evidence_ids") == evidence for o in m["occurrences"]):
        m["occurrences"].append({"document_id": DOC, "compound_label": extracted(label, label, evidence), "evidence_ids": evidence})


add_molecule(
    "27cis",
    "cis ketone fraction from S8 hydrolysis (unnumbered; corresponding cis isomer of 27a)",
    13,
    chirality="partially_specified",
    geometry="specified",
    assignment="source reports cis/trans fraction order but no absolute configuration",
    mixture="7:1 cis/trans fraction",
)
add_molecule(
    "trans_DMCHDA",
    "trans-N,N′-dimethylcyclohexane-1,2-diamine",
    17,
    chirality="partially_specified",
    geometry="not_applicable",
    assignment="trans relative relationship is named; absolute configuration/enantiomeric composition is unreported",
    note="Named ligand retained without a guessed isomeric SMILES or CAS number.",
)

for key in (
    "14a", "14b", "15a", "15b", "18a", "18b", "19a", "19b", "20a", "20b", "21a", "21b", "22a", "22b", "23a", "23b",
    "24a", "25a", "26", "S8", "27a", "28a", "29a", "30a", "31a", "32a", "37", "38", "39", "40", "42", "43", "44", "56", "57",
):
    pages = [p for p in range(13, 19) if key in {
        13: {"14a", "14b", "15a", "15b", "18a", "18b", "19a", "S8", "27a"},
        14: {"18a", "18b", "19a", "19b", "20a", "20b", "21a", "27a", "28a"},
        15: {"20a", "20b", "21a", "21b", "22a", "22b", "23a", "28a", "29a"},
        16: {"21a", "21b", "22a", "22b", "23a", "23b", "24a", "30a", "31a", "57"},
        17: {"22b", "23a", "25a", "26", "32a", "37", "38", "56", "57"},
        18: {"24a", "37", "38", "39", "40", "42", "43", "44"},
    }[p]]
    if pages:
        add_occurrence(key, f"compound {key}", pages)


def charge_fact(value, evidence):
    if value is None:
        return missing(evidence)
    raw, normalized = value if isinstance(value, tuple) else (json.dumps(value, ensure_ascii=False), value)
    return fact(raw=raw, normalized=normalized, evidence=evidence, method="original_extraction")


def P(key, role, **kwargs):
    return {"key": key, "role": role, **kwargs}


def make_participant(index, spec, evidence, default_reference):
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


LABELS = {"C-H": 1, "C-C": 2, "C-O": 3, "C-N": 5, "DB-C-C": 13, "DB-C-O": 14, "O-H": 24, "N-S": 31}
ATOM_PAIRS = {"C-H": ("C", "H"), "C-C": ("C", "C"), "C-O": ("C", "O"), "C-N": ("C", "N"), "DB-C-C": ("C", "C"), "DB-C-O": ("C", "O"), "O-H": ("O", "H"), "N-S": ("N", "S")}


def set_classification(r, labels, evidence, rationale):
    c = r["classification"]
    c["mapped_reaction_smiles"] = missing(evidence)
    changes = []
    label_values = []
    for idx, code in enumerate(labels, 1):
        a, b = ATOM_PAIRS[code]
        change_id = f"b{idx}"
        after_order = 2 if code.startswith("DB-") else 1
        change = {"change_id": change_id, "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}], "before_order": 0, "after_order": after_order, "change_type": "formed"}
        changes.append(inferred(rationale, change, evidence, "Bond class inferred from the source transformation; atom mapping is unavailable."))
        label_values.append({"category_number": LABELS[code], "code": code, "supporting_change_ids": [change_id]})
    c["bond_changes"] = changes
    c["labels"] = inferred(rationale, label_values, evidence, "Multi-label classification inferred under bond-formation-20261016.")
    c["coverage"] = inferred(rationale, "partial", evidence, "Only source-supported newly formed bond classes are asserted; full atom mapping is absent.")
    c["limitations"] = ["No atom-mapped reaction SMILES.", "Complete atom-level bond-change inventory remains unresolved."]
    c["model_votes"] = []


def make_outcome(index, spec, evidence):
    o = new_outcome(f"o{index}", spec["kind"], [mid(k) for k in spec["targets"]], evidence)
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


def add_reaction(suffix, pages, location, participants, steps, outcomes, labels, *, name, transformation, category="other", scale=None, solvent=None, notes=(), unresolved=(), general_procedure=None, evidence_override=None):
    evidence = evidence_override or ev(*pages)
    r = new_reaction(f"{DOC}:rxn:{suffix}", DOC, evidence)
    r["record_type"] = "specific_experiment"
    r["experiment_categories"] = extracted(category, [category], evidence)
    r["location"] = extracted(location, {"page": pages[0], "region": location}, evidence)
    r["general_scheme_ids"] = []
    r["reaction_name"] = extracted(name, name, evidence)
    r["transformation_type"] = extracted(transformation, transformation, evidence)
    first = next((p for p in participants if p["role"] == "reactant"), participants[0])
    default_reference = mid(first["key"])
    r["participants"] = [make_participant(i, p, evidence, default_reference) for i, p in enumerate(participants, 1)]
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
    r["steps"] = [make_step(i, s, evidence) for i, s in enumerate(steps, 1)]
    for key in r["special_conditions"]:
        r["special_conditions"][key] = not_applicable(evidence)
    r["procedure"]["original_text"] = extracted(location, location, evidence)
    r["procedure"]["general_procedure_reference"] = extracted(general_procedure, general_procedure, evidence) if general_procedure else missing(evidence)
    r["procedure"]["modifications"] = missing(evidence)
    r["procedure"]["addition_order"] = extracted("ordered as transcribed", [s["operation"]["normalized_value"] for s in r["steps"]], evidence)
    r["procedure"]["pretreatment"] = missing(evidence)
    r["outcomes"] = [make_outcome(i, o, evidence) for i, o in enumerate(outcomes, 1)]
    set_classification(r, labels, evidence, transformation)
    r["notes"] = [extracted(n, n, evidence) for n in notes]
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [], "unresolved": ["Machine record remains needs_review; schema validation is not chemistry verification.", "Full stereochemical structure encoding and atom mapping remain unresolved.", *unresolved]}
    r["completeness"] = "partial"
    bundle["reactions"].append(r)
    return r


def Y(targets, pct, *, raw=None, ytype="isolated yield", role="product", basis=None, scope=None):
    return {"kind": "yield", "targets": targets, "raw": raw or f"{pct}%", "value": q(pct, "%", raw_text=raw or f"{pct}%"), "yield_type": ytype, "product_role": role, "basis": basis, "scope": scope}


def A(targets, number, unit, *, raw=None, role="product", basis="reported isolated mass"):
    raw = raw or f"{number} {unit}"
    return {"kind": "product_amount", "targets": targets, "raw": raw, "value": q(number, unit, raw_text=raw), "product_role": role, "basis": basis}


def S(targets, raw, *, basis, scope, role="product mixture"):
    return {"kind": "selectivity", "targets": targets, "raw": raw, "value": raw, "product_role": role, "basis": basis, "scope": scope}


def mass_amt(key, mass_raw, mass_n, mass_unit, mmol_raw, mmol_n, role="reactant", **kw):
    return P(key, role, mass=Q(mass_raw, mass_n, mass_unit), amount=Q(mmol_raw, mmol_n, "mmol"), **kw)


def vol_amt(key, vol_raw, vol_n, mmol_raw, mmol_n, role="reagent", **kw):
    return P(key, role, volume=Q(vol_raw, vol_n, "mL"), amount=Q(mmol_raw, mmol_n, "mmol"), **kw)


def hydrolysis_run(suffix, page, substrate_keys, product_keys, substrate_mass, substrate_mmol, ratio, ppts_mass, ppts_mmol, acetone_ml, water_ml, outcomes, *, notes=(), unresolved=()):
    participants = [P(k, "reactant") for k in substrate_keys]
    participants[0]["mass"] = Q(substrate_mass[0], substrate_mass[1], substrate_mass[2])
    participants[0]["amount"] = Q(substrate_mmol[0], substrate_mmol[1], "mmol")
    participants += [mass_amt("PPTS", ppts_mass[0], ppts_mass[1], ppts_mass[2], ppts_mmol[0], ppts_mmol[1], role="catalyst"), P("acetone", "solvent", volume=Q(f"{acetone_ml} mL", acetone_ml, "mL")), P("water", "solvent", volume=Q(f"{water_ml} mL", water_ml, "mL"))]
    participants += [P(k, "product") for k in product_keys]
    if suffix == "rxn_p13_14ab_to_18ab_run1":
        participants += [P("EtOAc", "workup_agent", stage="workup", volume=Q("20 mL dilution", 20, "mL")), P("water", "workup_agent", stage="workup", volume=Q("2 × 5 mL wash", 10, "mL")), P("brine", "workup_agent", stage="workup", volume=Q("5 mL", 5, "mL")), P("MgSO4", "workup_agent", stage="workup")]
    steps = [
        {"operation": "charge substrate, acetone and water; add PPTS in one portion", "temperature": ("room temperature before heating", "r.t."), "time": None, "atmosphere": "nitrogen"},
        {"operation": "heat at reflux (oil-bath temperature 70 °C)", "temperature": Q("70 °C oil bath", 70, "°C"), "time": Q("16 h", 16, "h"), "atmosphere": "nitrogen"},
        {"operation": "cool to room temperature; work up and purify as General Procedure B", "stage": "workup", "temperature": ("room temperature", "r.t."), "time": None},
    ]
    return add_reaction(suffix, [page], f"page {page} General Procedure B / {suffix}", participants, steps, outcomes, ["DB-C-O"], name="acid-catalysed ketal hydrolysis", transformation="ketal deprotection to methyl ketone", scale={"raw_text": substrate_mmol[0], "number": substrate_mmol[1], "unit": "mmol"}, solvent={"raw": f"acetone {acetone_ml} mL / water {water_ml} mL", "composition": [mid("acetone"), mid("water")], "total_volume": Q(f"{acetone_ml + water_ml:g} mL", acetone_ml + water_ml, "mL")}, notes=(f"Starting material ratio: {ratio}.", *notes), unresolved=unresolved, general_procedure="General Procedure B")


# Page 13: S8 methylation, ketal hydrolyses, and the beginning of 19a Run 1.
def s8_run(run, mass, mmol, ratio, nah_mass, nah_mmol, mei_ml, mei_mmol, dmf_ml, product_mass, yld, product_ratio):
    participants = [mass_amt("14a", mass[0], mass[1], mass[2], mmol[0], mmol[1]), P("14b", "reactant"), mass_amt("NaH", nah_mass[0], nah_mass[1], "mg", nah_mmol[0], nah_mmol[1], role="reagent"), vol_amt("MeI", mei_ml[0], mei_ml[1], mei_mmol[0], mei_mmol[1]), P("DMF", "solvent", volume=Q(f"{dmf_ml} mL", dmf_ml, "mL")), P("S8", "product", mass=Q(product_mass[0], product_mass[1], "mg"))]
    if run == 1:
        participants += [P("NH4Cl", "workup_agent", stage="workup", volume=Q("20 mL saturated aqueous solution", 20, "mL")), P("Et2O", "workup_agent", stage="workup", volume=Q("3 × 25 mL", 75, "mL")), P("brine", "workup_agent", stage="workup", volume=Q("10 mL", 10, "mL")), P("MgSO4", "workup_agent", stage="workup")]
    return add_reaction(
        f"rxn_p13_S8_run{run}", [13], f"page 13 compound S8 Run {run}",
        participants,
        [
            {"operation": "dissolve 14a/b in DMF and cool", "temperature": Q("0 °C", 0, "°C"), "time": None, "atmosphere": "nitrogen"},
            {"operation": "add sodium hydride (60% in mineral oil) in one portion", "temperature": Q("0 °C", 0, "°C"), "time": Q("15 min", 15, "min")},
            {"operation": "add iodomethane dropwise over 1 min", "temperature": Q("0 °C", 0, "°C"), "time": Q("1 min addition", 1, "min")},
            {"operation": "hold at 0 °C", "temperature": Q("0 °C", 0, "°C"), "time": Q("15 min", 15, "min")},
            {"operation": "remove ice bath and warm to room temperature", "temperature": ("0 °C to room temperature", None), "time": Q("15 min", 15, "min")},
            {"operation": "quench into saturated ammonium chloride, extract with diethyl ether, wash brine, dry, filter, concentrate and chromatograph", "stage": "workup", "temperature": None, "time": None},
        ],
        [Y(["S8"], yld, raw=f"{product_mass[0]}, {yld}%", basis="isolated after silica-gel purification"), A(["S8"], product_mass[1], "mg", raw=product_mass[0]), S(["S8"], product_ratio, basis="reported product-mixture composition", scope="trans/cis")],
        ["C-N"], name="N-methylation of lactam 14a/b", transformation="formation of N-methyl lactam S8", scale={"raw_text": mmol[0], "number": mmol[1], "unit": "mmol"}, solvent={"raw": f"DMF {dmf_ml} mL", "composition": [mid("DMF")], "total_volume": Q(f"{dmf_ml} mL", dmf_ml, "mL")}, notes=(f"Starting 14a/b ratio: {ratio}.", "Run 1 was repeated three times in total with 67–74% yields; unreported individual repetitions were not expanded.") if run == 1 else (f"Starting 14a/b ratio: {ratio}.",), unresolved=("Absolute stereochemistry is not encoded; source uses relative trans/cis labels.",)
    )


s8_run(1, ("337 mg", 337, "mg"), ("1.03 mmol", 1.03), "approximately 10:1 trans/cis", ("62 mg", 62), ("1.5 mmol", 1.5), ("0.13 mL", 0.13), ("2.1 mmol", 2.1), 5.2, ("225 mg", 225), 67, "contains 6% cis isomer")
s8_run(2, ("430 mg", 430, "mg"), ("1.3 mmol", 1.3), "4:1 trans/cis", ("105 mg", 105), ("2.6 mmol", 2.6), ("0.21 mL", 0.21), ("3.3 mmol", 3.3), 7.0, ("321 mg", 321), 71, "4:1 trans/cis")

hydrolysis_run("rxn_p13_14ab_to_18ab_run1", 13, ["14a", "14b"], ["18a", "18b"], ("100 mg", 100, "mg"), ("0.304 mmol", 0.304), "5:1 trans/cis", ("38 mg", 38, "mg"), ("0.152 mmol", 0.152), 4.0, 0.3, [Y(["18a", "18b"], 73, raw="63 mg, 73%", ytype="combined isolated yield", basis="silica-gel plug (75% ethyl acetate/hexanes)"), A(["18a", "18b"], 63, "mg"), S(["18a", "18b"], "5:1 trans/cis", basis="integrated lactam methine resonances", scope="18a (trans)/18b (cis)")], notes=("The source reports nine total repetitions with yields 71–91% and average 83%; unreported repetitions were not expanded.",))
hydrolysis_run("rxn_p13_14ab_to_18ab_run2", 13, ["14a", "14b"], ["18a", "18b"], ("1.66 g", 1.66, "g"), ("5.33 mmol", 5.33), "3:1 trans/cis", ("670 mg", 670, "mg"), ("2.7 mmol", 2.7), 72.0, 4.0, [Y(["18a", "18b"], 91, raw="1.3 g, 91%", ytype="combined isolated yield", basis="column chromatography"), A(["18a", "18b"], 1.3, "g"), S(["18a", "18b"], "3:1 trans/cis", basis="reported product-mixture ratio", scope="18a (trans)/18b (cis)")])
hydrolysis_run("rxn_p13_18b_characterization_batch", 13, ["14a", "14b"], ["18a", "18b"], ("426 mg", 426, "mg"), ("1.30 mmol", 1.30), "3:1 trans/cis", ("163 mg", 163, "mg"), ("0.648 mmol", 0.648), 18.0, 1.0, [Y(["18a", "18b"], 88, raw="307 mg, 88%", ytype="combined isolated yield", basis="silica-gel plug, then repurification"), A(["18a", "18b"], 307, "mg"), S(["18a", "18b"], "3:1 trans/cis", basis="reported combined-product ratio", scope="18a (trans)/18b (cis)")], notes=("Repurification gave cis-enriched characterization fractions, but no independent cis-isomer mass or yield was reported.",), unresolved=("No independent yield can be assigned to 18b from the unquantified cis-enriched characterization fractions.",))
hydrolysis_run("rxn_p13_S8_to_27a_run1", 13, ["S8"], ["27a"], ("290 mg", 290, "mg"), ("0.882 mmol", 0.882), "contained 6% cis isomer", ("111 mg", 111, "mg"), ("0.441 mmol", 0.441), 12.0, 0.6, [Y(["27a"], 66, raw="164 mg, 66%", basis="column chromatography"), A(["27a"], 164, "mg")])
hydrolysis_run("rxn_p13_S8_to_27a_run2", 13, ["S8"], ["27a", "27cis"], ("490 mg", 490, "mg"), ("1.4 mmol", 1.4), "5:1 trans/cis", ("180 mg", 180, "mg"), ("0.72 mmol", 0.72), 19.0, 1.0, [Y(["27a"], 73, raw="310 mg, 73%", ytype="isolated trans-fraction yield", role="trans ketone 27a fraction", basis="column chromatography"), A(["27a"], 310, "mg", role="trans ketone fraction"), S(["27a", "27cis"], "19:1 trans/cis", basis="reported trans-fraction composition", scope="27a trans/cis composition", role="310 mg trans-enriched fraction"), Y(["27cis"], 4, raw="18 mg, 4%", ytype="isolated cis-fraction yield", role="unnumbered cis ketone fraction", basis="column chromatography"), A(["27cis"], 18, "mg", role="cis ketone fraction"), S(["27cis", "27a"], "7:1 cis/trans", basis="reported cis-fraction composition", scope="cis ketone/27a", role="18 mg cis-enriched fraction")], unresolved=("The source does not assign a compound number or complete absolute configuration to the isolated cis ketone fraction.",))

# Page 14: hydrolysis of 15a/b and General Procedure C additions.
hydrolysis_run("rxn_p13_p14_15a_to_19a_run1", 13, ["15a"], ["19a"], ("120 mg", 120, "mg"), ("0.306 mmol", 0.306), "pure trans", ("39 mg", 39, "mg"), ("0.15 mmol", 0.15), 4.1, 0.3, [Y(["19a"], 70, raw="74 mg, 70%", basis="silica plug (50% ethyl acetate/hexanes)"), A(["19a"], 74, "mg")], notes=("Procedure continues from page 13 onto page 14. Six total repetitions gave 58–85%, average 73%; unreported repetitions were not expanded.",))
# Extend evidence for the cross-page record to both inspected pages.
bundle["reactions"][-1]["evidence_ids"] = ev(13, 14)
def add_evidence_recursive(value, eids):
    if isinstance(value, dict):
        if "normalized_value" in value and "evidence_ids" in value:
            value["evidence_ids"] = list(dict.fromkeys(value["evidence_ids"] + eids))
        for child in value.values(): add_evidence_recursive(child, eids)
    elif isinstance(value, list):
        for child in value: add_evidence_recursive(child, eids)
add_evidence_recursive(bundle["reactions"][-1], ev(14))

hydrolysis_run("rxn_p14_15ab_to_19ab_run2", 14, ["15a", "15b"], ["19a", "19b"], ("1.74 g", 1.74, "g"), ("4.4 mmol", 4.4), "2:1 trans/cis", ("552 mg", 552, "mg"), ("2.2 mmol", 2.2), 60.0, 3.0, [Y(["19a", "19b"], 68, raw="1.05 g, 68%", ytype="combined isolated yield", basis="column chromatography"), A(["19a", "19b"], 1.05, "g"), S(["19a", "19b"], "2:1 trans/cis", basis="reported product-mixture ratio", scope="19a (trans)/19b (cis)")])
hydrolysis_run("rxn_p14_15ab_to_19ab_run3", 14, ["15a", "15b"], ["19a", "19b"], ("350 mg", 350, "mg"), ("0.89 mmol", 0.89), "2:1 trans/cis", ("112 mg", 112, "mg"), ("0.45 mmol", 0.45), 13.0, 1.0, [Y(["19a", "19b"], 85, raw="265 mg, 85%", ytype="combined isolated yield", basis="column chromatography"), A(["19a", "19b"], 265, "mg"), S(["19a", "19b"], "1.5:1 trans/cis", basis="reported product-mixture ratio", scope="19a (trans)/19b (cis)")])
hydrolysis_run("rxn_p14_19b_characterization_batch", 14, ["15a", "15b"], ["19a", "19b"], ("1.72 g", 1.72, "g"), ("4.39 mmol", 4.39), "4:1 trans/cis", ("552 mg", 552, "mg"), ("2.20 mmol", 2.20), 60.0, 3.0, [Y(["19a", "19b"], 80, raw="1.22 g, 80%", ytype="isolated mixed-fraction yield", role="first 19a/b fraction", basis="column chromatography"), A(["19a", "19b"], 1.22, "g", role="first 19a/b fraction"), S(["19a", "19b"], "4:1 trans/cis", basis="reported first-fraction composition", scope="19a/19b", role="first fraction"), Y(["19b", "19a"], 2, raw="30 mg, 2%", ytype="isolated mixed-fraction yield", role="second cis-enriched fraction", basis="column chromatography"), A(["19b", "19a"], 30, "mg", role="second cis-enriched fraction"), S(["19b", "19a"], "5:1 cis/trans", basis="reported second-fraction composition", scope="19b/19a", role="second fraction")])


def grignard_run(suffix, pages, substrate_keys, product_keys, substrate_mass, substrate_mmol, ratio, reagent_ml, reagent_mmol, thf_ml, outcomes, *, notes=(), unresolved=()):
    participants = [P(k, "reactant") for k in substrate_keys]
    participants[0]["mass"] = Q(substrate_mass[0], substrate_mass[1], substrate_mass[2]); participants[0]["amount"] = Q(substrate_mmol[0], substrate_mmol[1], "mmol")
    participants += [P("ethynylMgBr", "reagent", volume=Q(f"{reagent_ml} mL of 0.5 M solution", reagent_ml, "mL"), amount=Q(f"{reagent_mmol} mmol", reagent_mmol, "mmol")), P("THF", "solvent", volume=Q(f"{thf_ml} mL", thf_ml, "mL"))]
    participants += [P(k, "product") for k in product_keys]
    if suffix == "rxn_p14_18ab_to_20ab_run1":
        participants += [P("HCl", "workup_agent", stage="workup", volume=Q("25 mL of 1 N HCl", 25, "mL")), P("Et2O", "workup_agent", stage="workup", volume=Q("4 × 40 mL", 160, "mL")), P("brine", "workup_agent", stage="workup", volume=Q("20 mL", 20, "mL")), P("MgSO4", "workup_agent", stage="workup")]
    steps = [
        {"operation": "dissolve ketone in THF and cool", "temperature": Q("0 °C", 0, "°C"), "time": None, "atmosphere": "nitrogen"},
        {"operation": "add ethynylmagnesium bromide (0.5 M in THF) dropwise", "temperature": Q("0 °C", 0, "°C"), "time": Q("15 min addition", 15, "min")},
        {"operation": "maintain at 0 °C", "temperature": Q("0 °C", 0, "°C"), "time": Q("3 h", 3, "h")},
        {"operation": "quench at 0 °C with 1 N HCl; extract, wash, dry, filter, concentrate and chromatograph as General Procedure C", "stage": "workup", "temperature": Q("0 °C quench", 0, "°C"), "time": None},
    ]
    return add_reaction(suffix, pages, f"page {pages[0]} General Procedure C / {suffix}", participants, steps, outcomes, ["C-C", "O-H"], name="ethynylmagnesium bromide addition", transformation="methyl-ketone addition to form a tertiary propargyl alcohol", scale={"raw_text": substrate_mmol[0], "number": substrate_mmol[1], "unit": "mmol"}, solvent={"raw": f"THF {thf_ml} mL plus ethynylmagnesium bromide solution", "composition": [mid("THF")], "total_volume": Q(f"{thf_ml + reagent_ml:g} mL before quench", thf_ml + reagent_ml, "mL")}, notes=(f"Starting-material ratio: {ratio}.", *notes), unresolved=unresolved, general_procedure="General Procedure C")


grignard_run("rxn_p14_18ab_to_20ab_run1", [14], ["18a", "18b"], ["20a", "20b"], ("558 mg", 558, "mg"), ("2.08 mmol", 2.08), "approximately 3:1 trans/cis", 17.0, 8.5, 25.0, [Y(["20a"], 7, raw="41 mg, 7%", ytype="isolated trans-fraction yield", role="trans 20a fraction", basis="column chromatography"), A(["20a"], 41, "mg", role="trans fraction"), Y(["20b"], 1, raw="1 mg, 1%", ytype="isolated cis-fraction yield", role="cis 20b fraction", basis="column chromatography"), A(["20b"], 1, "mg", role="cis fraction"), Y(["20a", "20b"], 63, raw="385 mg, 63%", ytype="isolated mixed-fraction yield", role="4:1 trans/cis fraction", basis="column chromatography"), A(["20a", "20b"], 385, "mg", role="mixed fraction"), S(["20a", "20b"], "4:1 trans/cis", basis="reported mixed-fraction composition", scope="20a/20b"), S(["20a", "20b"], "1:1 at the newly formed stereocenter", basis="13C NMR", scope="diastereomers at the tertiary propargyl alcohol center")], notes=("Six reported repetitions gave 65–81%; unreported repetitions were not expanded.",))
grignard_run("rxn_p14_18ab_to_20ab_run2", [14], ["18a", "18b"], ["20a", "20b"], ("1.30 g", 1.30, "g"), ("4.9 mmol", 4.9), "3:1 trans/cis", 38.0, 19.0, 60.0, [Y(["20a", "20b"], 78, raw="1.13 g, 78%", ytype="combined isolated yield", basis="column chromatography"), A(["20a", "20b"], 1.13, "g"), S(["20a", "20b"], "3:1 trans/cis", basis="reported product-mixture ratio", scope="20a/20b")])
grignard_run("rxn_p14_18ab_to_20ab_run3", [14], ["18a", "18b"], ["20a", "20b"], ("640 mg", 640, "mg"), ("2.4 mmol", 2.4), "4:1 trans/cis", 19.0, 9.5, 29.0, [Y(["20a", "20b"], 75, raw="528 mg, 75%", ytype="isolated mixed-fraction yield", role="trans-enriched fraction", basis="column chromatography"), A(["20a", "20b"], 528, "mg", role="trans-enriched fraction"), S(["20a", "20b"], "6:1 trans/cis", basis="reported first-fraction composition", scope="20a/20b", role="528 mg fraction"), Y(["20b", "20a"], 6, raw="39 mg, 6%", ytype="isolated mixed-fraction yield", role="cis-enriched fraction", basis="column chromatography"), A(["20b", "20a"], 39, "mg", role="cis-enriched fraction"), S(["20b", "20a"], "2:1 cis/trans", basis="reported second-fraction composition", scope="20b/20a", role="39 mg fraction")])
grignard_run("rxn_p14_27a_to_28a", [14], ["27a"], ["28a"], ("310 mg", 310, "mg"), ("1.04 mmol", 1.04), "7:1 trans/cis", 6.3, 3.1, 14.0, [Y(["28a"], 62, raw="198 mg, 62%", basis="column chromatography (20–60% ethyl acetate/hexanes)"), A(["28a"], 198, "mg"), S(["28a"], "7:1 trans/cis", basis="reported product-mixture ratio", scope="lactam trans/cis"), S(["28a"], "1:1 diastereomers", basis="13C NMR", scope="new tertiary alcohol stereocenter"), {"kind": "other", "targets": ["27a"], "raw": "NMR contains approximately 7% ketone starting material", "value": q(7, "%", qualifier="approx", raw_text="~7% ketone starting material"), "product_role": "residual starting material in characterized 28a", "basis": "NMR", "scope": "27a impurity in 28a sample"}], unresolved=("The source does not assign absolute configuration at the newly formed center.",))

# Page 15: 21a/21b, 29a and propargyl pivalates 22a/b.
grignard_run("rxn_p14_p15_19a_to_21a_run1", [14, 15], ["19a"], ["21a"], ("270 mg", 270, "mg"), ("0.78 mmol", 0.78), "pure trans", 6.2, 3.1, 9.3, [Y(["21a"], 91, raw="265 mg, 91%", basis="isolated product"), A(["21a"], 265, "mg"), S(["21a"], "1:1 diastereomers", basis="13C NMR", scope="new tertiary alcohol stereocenter")], notes=("Five reported repetitions gave 64–91%, average 82%; unreported repetitions were not expanded.",))
grignard_run("rxn_p15_19ab_to_21ab_cis_characterization", [15], ["19a", "19b"], ["21a", "21b"], ("1.22 g", 1.22, "g"), ("3.51 mmol", 3.51), "4:1 trans/cis", 28.0, 14.0, 42.0, [Y(["21a", "21b"], 82, raw="1.07 g, 82%", ytype="isolated mixed-fraction yield", role="predominantly trans fraction", basis="column chromatography"), A(["21a", "21b"], 1.07, "g", role="predominantly trans fraction"), S(["21a", "21b"], "4:1 trans/cis", basis="reported first-fraction composition", scope="21a/21b", role="1.07 g fraction"), Y(["21b", "21a"], 4, raw="52 mg, 4%", ytype="isolated mixed-fraction yield", role="predominantly cis fraction", basis="column chromatography"), A(["21b", "21a"], 52, "mg", role="cis-enriched fraction"), S(["21b", "21a"], "2:1 cis/trans", basis="reported second-fraction composition", scope="21b/21a", role="52 mg fraction"), S(["21b"], "approximately 1.5:1 diastereomers", basis="estimated from 13C NMR peaks; not separately integrable in 1H NMR", scope="new tertiary alcohol stereocenter", role="characterized cis-enriched material"), {"kind": "other", "targets": ["21a", "21b"], "raw": "NMR contains approximately 56% trans isomer", "value": q(56, "%", qualifier="approx", raw_text="~56% trans isomer"), "product_role": "characterized cis-enriched material", "basis": "NMR", "scope": "trans-isomer content"}], unresolved=("The reported 2:1 cis/trans fraction and approximately 56% trans content in the characterization sample refer to different descriptions; both are preserved without reconciliation.",))

add_reaction(
    "rxn_p15_28a_to_29a_acetylation", [15], "page 15 compound 29a direct acetylation",
    [mass_amt("28a", "86 mg", 86, "mg", "0.27 mmol", 0.27), mass_amt("DMAP", "4 mg", 4, "mg", "0.03 mmol", 0.03, role="catalyst"), vol_amt("Et3N", "0.37 mL", 0.37, "2.7 mmol", 2.7, role="reagent"), vol_amt("Ac2O", "0.12 mL", 0.12, "1.3 mmol", 1.3), P("DCM", "solvent", volume=Q("1.3 mL", 1.3, "mL")), P("29a", "product", mass=Q("50 mg", 50, "mg")), P("DCM", "workup_agent", stage="workup", volume=Q("25 mL dilution", 25, "mL")), P("NH4Cl", "workup_agent", stage="workup", volume=Q("5 mL saturated aqueous solution", 5, "mL")), P("brine", "workup_agent", stage="workup", volume=Q("5 mL", 5, "mL")), P("MgSO4", "workup_agent", stage="workup")],
    [{"operation": "charge DMAP and 28a in dichloromethane; cool", "temperature": Q("0 °C", 0, "°C"), "time": None, "atmosphere": "nitrogen"}, {"operation": "add triethylamine dropwise", "temperature": Q("0 °C", 0, "°C"), "time": Q("2 min addition", 2, "min")}, {"operation": "add acetic anhydride dropwise", "temperature": Q("0 °C", 0, "°C"), "time": Q("1 min addition", 1, "min")}, {"operation": "remove ice bath, warm to room temperature and stir", "temperature": ("0 °C to room temperature", None), "time": Q("3 h", 3, "h")}, {"operation": "dilute, wash with saturated ammonium chloride and brine, dry, filter, concentrate and chromatograph", "stage": "workup", "temperature": None, "time": None}],
    [Y(["29a"], 54, raw="50 mg, 54%", basis="flash chromatography (40–50% ethyl acetate/hexanes)"), A(["29a"], 50, "mg"), S(["29a"], "dr 1:1", basis="reported product dr", scope="diastereomers at the propargylic center")], ["C-O"], name="acetylation of propargyl alcohol 28a", transformation="formation of propargyl acetate 29a", scale={"raw_text": "28a 0.27 mmol", "number": 0.27, "unit": "mmol"}, solvent={"raw": "CH2Cl2 1.3 mL", "composition": [mid("DCM")], "total_volume": Q("1.3 mL", 1.3, "mL")}
)
add_reaction(
    "rxn_p15_28a_to_29a_gpD_alt", [15], "page 15 one-sentence alternative: following General Procedure D afforded 29a",
    [P("28a", "reactant"), P("29a", "product")],
    [{"operation": "follow General Procedure D as stated; individual charge and any modification were not reported for this run", "temperature": None, "time": None}],
    [Y(["29a"], 30, raw="30% yield of 29a", basis="after column chromatography")], ["C-O"], name="source-stated General Procedure D preparation of 29a", transformation="source identifies 29a as the product, despite General Procedure D being the pivaloylation protocol", notes=("The source sentence is preserved verbatim in substance and is chemically inconsistent with General Procedure D, which specifies pivalic anhydride and pivalate products.",), unresolved=("Reagents, quantities, and the reason General Procedure D is said to afford acetate 29a are unresolved; no correction or reagent inference was made.",), general_procedure="General Procedure D (source-stated, potentially inconsistent)"
)


def pivalate_run(suffix, pages, substrate_keys, product_keys, substrate_mass, substrate_mmol, ratio, piv_ml, piv_mmol, sc_mass, sc_mmol, mecn_ml, outcomes, *, notes=(), unresolved=()):
    participants = [P(k, "reactant") for k in substrate_keys]
    participants[0]["mass"] = Q(substrate_mass[0], substrate_mass[1], substrate_mass[2]); participants[0]["amount"] = Q(substrate_mmol[0], substrate_mmol[1], "mmol")
    participants += [vol_amt("piv2O", piv_ml[0], piv_ml[1], piv_mmol[0], piv_mmol[1]), mass_amt("ScOTf3", sc_mass[0], sc_mass[1], sc_mass[2], sc_mmol[0], sc_mmol[1], role="catalyst"), P("MeCN", "solvent", volume=Q(f"{mecn_ml} mL", mecn_ml, "mL"))]
    participants += [P(k, "product") for k in product_keys]
    if suffix == "rxn_p15_20a_to_22a_run1":
        participants += [P("Et2O", "workup_agent", stage="workup", volume=Q("20 mL dilution plus 2 × 10 mL extraction", 40, "mL")), P("NaHCO3", "workup_agent", stage="workup", volume=Q("10 mL saturated solution", 10, "mL")), P("brine", "workup_agent", stage="workup"), P("MgSO4", "workup_agent", stage="workup")]
    steps = [{"operation": "dissolve propargyl alcohol in acetonitrile; add pivalic anhydride then scandium(III) triflate", "temperature": ("room temperature", "r.t."), "time": None, "atmosphere": "nitrogen"}, {"operation": "stir at room temperature", "temperature": ("room temperature", "r.t."), "time": Q("16 h", 16, "h")}, {"operation": "dilute with diethyl ether, wash saturated sodium bicarbonate, extract, wash brine, dry, filter, concentrate and chromatograph as General Procedure D", "stage": "workup", "temperature": None, "time": None}]
    return add_reaction(suffix, pages, f"page {pages[0]} General Procedure D / {suffix}", participants, steps, outcomes, ["C-O"], name="scandium-catalysed pivaloylation", transformation="conversion of propargyl alcohol to propargyl pivalate", scale={"raw_text": substrate_mmol[0], "number": substrate_mmol[1], "unit": "mmol"}, solvent={"raw": f"acetonitrile {mecn_ml} mL", "composition": [mid("MeCN")], "total_volume": Q(f"{mecn_ml} mL", mecn_ml, "mL")}, notes=(f"Starting-material ratio: {ratio}.", *notes), unresolved=unresolved, general_procedure="General Procedure D")


pivalate_run("rxn_p15_20a_to_22a_run1", [15], ["20a"], ["22a"], ("56 mg", 56, "mg"), ("0.18 mmol", 0.18), "pure trans", ("0.05 mL", 0.05), ("0.2 mmol", 0.2), ("36 mg", 36, "mg"), ("0.073 mmol", 0.073), 0.72, [Y(["22a"], 76, raw="54 mg, 76%", basis="silica chromatography (15–30% ethyl acetate/hexanes)"), A(["22a"], 54, "mg"), S(["22a"], "dr 1:1", basis="reported product dr", scope="propargylic-center diastereomers")])
pivalate_run("rxn_p15_20ab_to_22ab_run2", [15], ["20a", "20b"], ["22a", "22b"], ("368 mg", 368, "mg"), ("1.25 mmol", 1.25), "3.5:1 trans/cis", ("0.36 mL", 0.36), ("1.8 mmol", 1.8), ("246 mg", 246, "mg"), ("0.50 mmol", 0.50), 5.0, [Y(["22a"], 61, raw="299 mg, 61%", ytype="isolated trans-fraction yield", role="trans 22a", basis="column chromatography"), A(["22a"], 299, "mg", role="trans 22a"), Y(["22b"], 16, raw="79 mg, 16%", ytype="isolated cis-fraction yield", role="cis 22b", basis="column chromatography"), A(["22b"], 79, "mg", role="cis 22b")])
pivalate_run("rxn_p15_20ab_to_22ab_run3", [15], ["20a", "20b"], ["22a", "22b"], ("1.13 g", 1.13, "g"), ("3.85 mmol", 3.85), "3:1 trans/cis", ("0.94 mL", 0.94), ("4.62 mmol", 4.62), ("758 mg", 758, "mg"), ("1.54 mmol", 1.54), 15.0, [Y(["22a"], 34, raw="490 mg, 34%", ytype="isolated trans-fraction yield", role="trans 22a", basis="column chromatography"), A(["22a"], 490, "mg", role="trans 22a"), Y(["22b"], 11, raw="163 mg, 11%", ytype="isolated cis-fraction yield", role="cis 22b", basis="column chromatography"), A(["22b"], 163, "mg", role="cis 22b"), Y(["22a", "22b"], 12, raw="176 mg, 12%", ytype="isolated mixed-fraction yield", role="4:1 trans/cis fraction", basis="column chromatography"), A(["22a", "22b"], 176, "mg", role="mixed fraction"), S(["22a", "22b"], "4:1 trans/cis", basis="reported mixed-fraction composition", scope="22a/22b")])
pivalate_run("rxn_p15_22b_characterization_batch", [15], ["20a", "20b"], ["22a", "22b"], ("140 mg", 140, "mg"), ("0.48 mmol", 0.48), "5:1 trans/cis", ("0.12 mL", 0.12), ("0.57 mmol", 0.57), ("94 mg", 94, "mg"), ("0.19 mmol", 0.19), 1.9, [Y(["22a"], 55, raw="99 mg, 55%", ytype="isolated trans-fraction yield", role="trans 22a", basis="column chromatography"), A(["22a"], 99, "mg", role="trans fraction"), Y(["22b"], 11, raw="20 mg, 11%", ytype="isolated cis-fraction yield", role="cis 22b", basis="column chromatography"), A(["22b"], 20, "mg", role="cis fraction"), Y(["22a", "22b"], 4, raw="7 mg, 4%", ytype="isolated mixed-fraction yield", role="2:1 trans/cis mixture", basis="column chromatography"), A(["22a", "22b"], 7, "mg", role="mixed fraction"), S(["22a", "22b"], "2:1 trans/cis", basis="reported mixed-fraction composition", scope="22a/22b"), S(["22b"], "1.5:1 dr", basis="1H NMR integration of terminal-alkyne resonances at 2.48 (major) and 2.43 ppm (minor)", scope="22b propargylic-center diastereomers")])

# Page 16: pivalates 23a/b and Stryker-reagent allene formation.
pivalate_run("rxn_p15_p16_21a_to_23a_run1", [15, 16], ["21a"], ["23a"], ("100 mg", 100, "mg"), ("0.27 mmol", 0.27), "pure trans", ("0.07 mL", 0.07), ("0.35 mmol", 0.35), ("53 mg", 53, "mg"), ("0.11 mmol", 0.11), 1.1, [Y(["23a"], 96, raw="118 mg, 96%", basis="column chromatography (20–50% ethyl acetate/hexanes)"), A(["23a"], 118, "mg"), S(["23a"], "1:1 diastereomers", basis="13C NMR", scope="propargylic-center diastereomers")])
pivalate_run("rxn_p16_21ab_to_23ab_run2", [16], ["21a", "21b"], ["23a", "23b"], ("230 mg", 230, "mg"), ("0.62 mmol", 0.62), "2.5:1 trans/cis", ("0.16 mL", 0.16), ("0.80 mmol", 0.80), ("121 mg", 121, "mg"), ("0.25 mmol", 0.25), 2.5, [Y(["23a"], 39, raw="110 mg, 39%", ytype="isolated trans-fraction yield", role="trans 23a", basis="column chromatography"), A(["23a"], 110, "mg", role="trans 23a"), Y(["23b"], 16, raw="44 mg, 16%", ytype="isolated cis-fraction yield", role="cis 23b", basis="column chromatography"), A(["23b"], 44, "mg", role="cis 23b"), Y(["23a", "23b"], 55, raw="55% total yield", ytype="sum of isolated fraction yields", role="all isolated 23a/23b fractions", basis="reported total yield")])
pivalate_run("rxn_p16_23b_run3", [16], ["21a", "21b"], ["23a", "23b"], ("1.07 g", 1.07, "g"), ("2.87 mmol", 2.87), "4:1 trans/cis", ("0.70 mL", 0.70), ("3.5 mmol", 3.5), ("565 mg", 565, "mg"), ("1.15 mmol", 1.15), 12.0, [Y(["23a"], 43, raw="565 mg, 43%", ytype="isolated trans-fraction yield", role="trans 23a", basis="column chromatography"), A(["23a"], 565, "mg", role="trans 23a"), Y(["23b"], 11, raw="139 mg, 11%", ytype="isolated cis-fraction yield", role="cis 23b", basis="column chromatography"), A(["23b"], 139, "mg", role="cis 23b"), Y(["23a", "23b"], 5, raw="72 mg, 5%", ytype="isolated mixed-fraction yield", role="2:1 trans/cis mixture", basis="column chromatography"), A(["23a", "23b"], 72, "mg", role="mixed fraction"), S(["23a", "23b"], "2:1 trans/cis", basis="reported mixed-fraction composition; ratios based on methylene resonances at 5.38 and 5.40 ppm", scope="23a/23b")])


def allene_run(suffix, pages, substrate, product_keys, sm_mass, sm_mmol, reagent_mass, reagent_mmol, toluene_ml, water_ml, temperature, time_h, outcomes, *, notes=(), unresolved=()):
    participants = [mass_amt(substrate, sm_mass[0], sm_mass[1], sm_mass[2], sm_mmol[0], sm_mmol[1]), mass_amt("Stryker", reagent_mass[0], reagent_mass[1], "mg", reagent_mmol[0], reagent_mmol[1], role="reagent"), P("toluene", "solvent", volume=Q(f"{toluene_ml} mL", toluene_ml, "mL"))]
    if water_ml is not None: participants.append(P("water", "additive", volume=Q(f"{water_ml} mL", water_ml, "mL")))
    participants += [P(k, "product") for k in product_keys]
    if suffix == "rxn_p16_22a_to_24a_gpE":
        participants += [P("NH4Cl", "workup_agent", stage="workup"), P("Et2O", "workup_agent", stage="workup", volume=Q("10 mL dilution plus 3 × 10 mL extraction", 40, "mL")), P("MgSO4", "workup_agent", stage="workup"), P("silica", "purification_agent", stage="purification")]
    steps = [{"operation": "charge substrate, toluene and water; degas by bubbling nitrogen", "temperature": ("room temperature before cooling", "r.t."), "time": Q("10 min degassing", 10, "min"), "atmosphere": "nitrogen"}, {"operation": "cool reaction mixture", "temperature": Q(temperature, float(temperature.split()[0].replace("−", "-")), "°C"), "time": None, "atmosphere": "nitrogen"}, {"operation": "add Stryker reagent in one portion after approximately 1 min air exposure; evacuate/refill with nitrogen three times", "temperature": Q(temperature, float(temperature.split()[0].replace("−", "-")), "°C"), "time": None, "atmosphere": "nitrogen"}, {"operation": "stir under nitrogen", "temperature": Q(temperature, float(temperature.split()[0].replace("−", "-")), "°C"), "time": Q(f"{time_h} h", time_h, "h"), "atmosphere": "nitrogen"}, {"operation": "quench into cooled saturated ammonium chloride, dilute with diethyl ether, stir open to air 30 min, extract, dry, filter and chromatograph as General Procedure E", "stage": "workup", "temperature": Q("0 °C initial quench", 0, "°C"), "time": Q("30 min open-air stir", 30, "min"), "atmosphere": "air"}]
    return add_reaction(suffix, pages, f"page {pages[0]} General Procedure E / {suffix}", participants, steps, outcomes, ["DB-C-C", "C-H"], name="copper-hydride formation of 3,3-disubstituted allene", transformation="reductive rearrangement of propargyl carboxylate to allene", scale={"raw_text": sm_mmol[0], "number": sm_mmol[1], "unit": "mmol"}, solvent={"raw": f"toluene {toluene_ml} mL", "composition": [mid("toluene")], "total_volume": Q(f"{toluene_ml} mL", toluene_ml, "mL")}, notes=notes, unresolved=unresolved, general_procedure="General Procedure E")


allene_run("rxn_p16_22a_to_24a_gpE", [16], "22a", ["24a"], ("107 mg", 107, "mg"), ("0.272 mmol", 0.272), ("533 mg", 533), ("0.272 mmol", 0.272), 7.0, 0.01, "−10 °C", 2, [Y(["24a"], 33, raw="24 mg, 33%", basis="silica-gel flash chromatography (15–30% ethyl acetate/hexanes)"), A(["24a"], 24, "mg")], notes=("The source reports reagent-source/temperature observations without complete individual run quantities: aged Sigma reagent gave no observed α-methylene reduction; new Acros reagent at room temperature gave 24a/31a 1:1; lowering to −10 °C gave 10:1; other batches/fresh reagent at −10 °C gave little over-reduction.",), unresolved=("The reagent-source observations lack individual yields and complete charges and were retained as context rather than fabricated separate runs.",))
allene_run("rxn_p16_29a_to_30a", [16], "29a", ["30a"], ("69 mg", 69, "mg"), ("0.20 mmol", 0.20), ("385 mg", 385), ("0.196 mmol", 0.196), 5.0, 0.01, "−10 °C", 2, [Y(["30a"], 53, raw="30 mg, 53%", basis="column chromatography (10–20% ethyl acetate/hexanes)"), A(["30a"], 30, "mg")], notes=("Baseline impurities from triphenylphosphine byproducts were present in the characterization spectra.",))
allene_run("rxn_p16_22b_to_57", [16], "22b", ["57"], ("75 mg", 75, "mg"), ("0.19 mmol", 0.19), ("373 mg", 373), ("0.19 mmol", 0.19), 5.0, 0.01, "−10 °C", 2, [Y(["57"], 36, raw="19 mg, 36%", basis="column chromatography (15–25% ethyl acetate/hexanes)"), A(["57"], 19, "mg"), {"kind": "other", "targets": ["57"], "raw": "the only product", "value": "only product reported", "product_role": "57", "basis": "source statement after chromatography", "scope": "observed product identity"}], notes=("This is the cis-lactam substrate run; product 57 is the α-methyl lactam produced by reduction of the exocyclic methylene during allene formation.",))

# Page 17: 25a, desilylation runs, and N-arylation of 24a.
allene_run("rxn_p17_23a_to_25a", [17], "23a", ["25a", "32a"], ("110 mg", 110, "mg"), ("0.24 mmol", 0.24), ("424 mg", 424), ("0.22 mmol", 0.22), 6.0, None, "−20 °C", 2, [Y(["25a"], 46, raw="39 mg, 46%", ytype="isolated product yield", role="allene 25a", basis="column chromatography (5–25% ethyl acetate/hexanes)"), A(["25a"], 39, "mg", role="allene 25a"), Y(["32a"], 11, raw="9 mg, 11%", ytype="isolated byproduct yield", role="over-reduced 32a", basis="column chromatography"), A(["32a"], 9, "mg", role="over-reduced 32a"), S(["25a", "32a"], "4:1", basis="reported isolated-product ratio", scope="25a/32a")], notes=("Other runs gave 25a in 31–53% and 32a in 0–31%; these aggregate ranges were not expanded into unreported runs.",))


def tbaf_run(suffix, substrate, sm_mass, sm_mmol, tbaf_ml, tbaf_mmol, thf_ml, product, product_mass, yld, *, ratio=None):
    outs = [Y([product], yld, raw=f"{product_mass[0]}, {yld}%", basis="silica-gel flash chromatography"), A([product], product_mass[1], "mg", raw=product_mass[0])]
    if ratio: outs.append(S([product], ratio, basis="integrated lactam methine resonances", scope="diastereomeric ratio"))
    participants = [mass_amt(substrate, sm_mass[0], sm_mass[1], "mg", sm_mmol[0], sm_mmol[1]), P("TBAF", "reagent", volume=Q(f"{tbaf_ml} mL of 1 M solution", tbaf_ml, "mL"), amount=Q(f"{tbaf_mmol} mmol", tbaf_mmol, "mmol")), P("THF", "solvent", volume=Q(f"{thf_ml} mL", thf_ml, "mL")), P(product, "product", mass=Q(product_mass[0], product_mass[1], "mg"))]
    if suffix == "rxn_p17_25a_to_26_run1":
        participants += [P("NH4Cl", "workup_agent", stage="workup", volume=Q("2 mL saturated aqueous solution", 2, "mL")), P("EtOAc", "workup_agent", stage="workup", volume=Q("4 × 5 mL", 20, "mL")), P("brine", "workup_agent", stage="workup", volume=Q("5 mL", 5, "mL")), P("MgSO4", "workup_agent", stage="workup")]
    return add_reaction(suffix, [17], f"page 17 General Procedure F / {suffix}", participants, [{"operation": "dissolve TIPS allene in THF and cool", "temperature": Q("0 °C", 0, "°C"), "time": None, "atmosphere": "nitrogen"}, {"operation": "add TBAF (1 M in THF) dropwise and stir", "temperature": Q("0 °C", 0, "°C"), "time": Q("45 min", 45, "min")}, {"operation": "add saturated ammonium chloride while cooled; extract with ethyl acetate, wash brine, dry, filter, concentrate and chromatograph", "stage": "workup", "temperature": Q("0 °C initial quench", 0, "°C"), "time": None}], outs, ["C-H"], name="TBAF desilylation", transformation="replacement of terminal TIPS substituent by hydrogen", scale={"raw_text": sm_mmol[0], "number": sm_mmol[1], "unit": "mmol"}, solvent={"raw": f"THF {thf_ml} mL plus TBAF solution", "composition": [mid("THF")], "total_volume": Q(f"{thf_ml + tbaf_ml:g} mL before quench", thf_ml + tbaf_ml, "mL")}, general_procedure="General Procedure F")


tbaf_run("rxn_p17_25a_to_26_run1", "25a", ("38 mg", 38), ("0.11 mmol", 0.11), 0.15, 0.15, 1.1, "26", ("12 mg", 12), 57)
bundle["reactions"][-1]["notes"].append(extracted("Four total repetitions gave 57–93%, average 79%; unreported repetitions were not expanded.", "Four total repetitions gave 57–93%, average 79%; unreported repetitions were not expanded.", ev(17)))
tbaf_run("rxn_p17_25a_to_26_run2", "25a", ("40 mg", 40), ("0.11 mmol", 0.11), 0.17, 0.17, 1.1, "26", ("18 mg", 18), 82)
tbaf_run("rxn_p17_32a_to_56", "32a", ("26 mg", 26), ("0.072 mmol", 0.072), 0.11, 0.10, 0.8, "56", ("14 mg", 14), 86, ratio="9:1 dr")


def narylation_run(suffix, pages, aryl_key, product, sm_mass, sm_mmol, aryl_charge, cui_mass, cui_mmol, ligand_key, ligand_ul, ligand_mmol, base_key, base_mass, base_mmol, toluene_ml, product_mass, yld, *, notes=(), unresolved=()):
    aryl = P(aryl_key, "reactant")
    if aryl_charge[2] == "mg": aryl["mass"] = Q(aryl_charge[0], aryl_charge[1], "mg")
    else: aryl["volume"] = Q(aryl_charge[0], aryl_charge[1], aryl_charge[2])
    aryl["amount"] = Q(aryl_charge[3], aryl_charge[4], "mmol")
    participants = [mass_amt("24a", sm_mass[0], sm_mass[1], "mg", sm_mmol[0], sm_mmol[1]), aryl, mass_amt("CuI", cui_mass[0], cui_mass[1], "mg", cui_mmol[0], cui_mmol[1], role="catalyst"), P(ligand_key, "ligand", volume=Q(f"{ligand_ul} μL", ligand_ul, "μL"), amount=Q(f"{ligand_mmol} mmol", ligand_mmol, "mmol")), mass_amt(base_key, base_mass[0], base_mass[1], "mg", base_mmol[0], base_mmol[1], role="reagent"), P("toluene", "solvent", volume=Q(f"{toluene_ml} mL", toluene_ml, "mL")), P(product, "product", mass=Q(f"{product_mass} mg", product_mass, "mg"))]
    steps = [{"operation": "charge CuI and base in a flame-dried tube; cap, evacuate and fill with nitrogen three times", "temperature": ("room temperature", "r.t."), "time": None, "atmosphere": "nitrogen"}, {"operation": "add aryl iodide and ligand, then add 24a in toluene; evacuate/fill with nitrogen three times", "temperature": ("room temperature", "r.t."), "time": None, "atmosphere": "nitrogen"}, {"operation": "place in preheated oil bath", "temperature": Q("80 °C", 80, "°C"), "time": Q("20 h", 20, "h"), "atmosphere": "nitrogen"}, {"operation": "cool, filter through silica, concentrate and chromatograph", "stage": "purification", "temperature": ("room temperature", "r.t."), "time": None}]
    return add_reaction(suffix, pages, f"page {pages[0]} General Procedure G / {suffix}", participants, steps, [Y([product], yld, raw=f"{product_mass} mg, {yld}%", basis="silica-gel column chromatography"), A([product], product_mass, "mg")], ["C-N"], name="copper-catalysed N-arylation", transformation=f"N-arylation of 24a to {product}", scale={"raw_text": sm_mmol[0], "number": sm_mmol[1], "unit": "mmol"}, solvent={"raw": f"toluene {toluene_ml} mL", "composition": [mid("toluene")], "total_volume": Q(f"{toluene_ml} mL", toluene_ml, "mL")}, notes=notes, unresolved=unresolved, general_procedure="General Procedure G")


narylation_run("rxn_p17_24a_33_to_37", [17], "33", "37", ("21 mg", 21), ("0.076 mmol", 0.076), ("13 μL", 13, "μL", "0.11 mmol", 0.11), ("4 mg", 4), ("0.02 mmol", 0.02), "trans_DMCHDA", 8, 0.05, "Cs2CO3", ("50 mg", 50), ("0.15 mmol", 0.15), 0.8, 12, 44)
# Alternative 37 run: detailed charge conflicts with the preceding prose ligand statement.
narylation_run("rxn_p17_24a_33_to_37_K2CO3", [17], "33", "37", ("10 mg", 10), ("0.033 mmol", 0.033), ("5 μL", 5, "μL", "0.04 mmol", 0.04), ("5 mg", 5), ("0.03 mmol", 0.03), "trans_DMCHDA", 4, 0.03, "K2CO3", ("10 mg", 10), ("0.065 mmol", 0.065), 0.4, 3, 27, notes=("The preceding sentence says the ligand was changed to N,N′-dimethylethylenediamine, but the detailed charge lists trans-N,N′-dimethylcyclohexane-1,2-diamine (4 μL, 0.03 mmol). The detailed charge is represented and the conflict is retained.", "Starting material was recovered, but no amount or recovery percentage was reported."), unresolved=("Ligand identity conflicts between prose and the detailed charge; no authorial correction was made.", "Starting-material recovery amount is unreported."))

# Page 18 N-arylation scope entries.
narylation_run("rxn_p17_p18_24a_34_to_38", [17, 18], "34", "38", ("33 mg", 33), ("0.12 mmol", 0.12), ("38 mg", 38, "mg", "0.17 mmol", 0.17), ("5 mg", 5), ("0.02 mmol", 0.02), "DMEDA", 5, 0.05, "Cs2CO3", ("78 mg", 78), ("0.24 mmol", 0.24), 1.2, 19, 42)
narylation_run("rxn_p18_24a_35_to_39", [18], "35", "39", ("24 mg", 24), ("0.087 mmol", 0.087), ("0.02 mL", 0.02, "mL", "0.1 mmol", 0.1), ("2 mg", 2), ("0.02 mmol", 0.02), "DMEDA", 4, 0.04, "Cs2CO3", ("57 mg", 57), ("0.17 mmol", 0.17), 0.9, 17, 47)
add_reaction("rxn_p18_24a_35_to_39_K2CO3_aggregate", [18], "page 18 aggregate K2CO3 variant for compound 39", [P("24a", "reactant"), P("35", "reactant"), P("CuI", "catalyst"), P("DMEDA", "ligand"), P("K2CO3", "reagent"), P("toluene", "solvent"), P("39", "product")], [{"operation": "use potassium carbonate under the N-arylation conditions; individual charges and run count are not stated", "temperature": Q("80 °C (General Procedure G)", 80, "°C"), "time": Q("20 h (General Procedure G)", 20, "h"), "atmosphere": "nitrogen"}], [Y(["39"], 14, raw="14–16% isolated 39", ytype="isolated-yield range", basis="source reports an aggregate range", scope="minimum of 14–16% range"), {"kind": "yield", "targets": ["39"], "raw": "14–16% isolated 39", "value": {"minimum": 14, "maximum": 16, "unit": "%", "raw_text": "14–16%"}, "yield_type": "isolated-yield range", "product_role": "39", "basis": "source reports aggregate range", "scope": "range across unspecified run(s)"}, {"kind": "recovery", "targets": ["24a"], "raw": "20–25% starting material recovered", "value": {"minimum": 20, "maximum": 25, "unit": "%", "raw_text": "20–25%"}, "yield_type": "recovery", "product_role": "recovered 24a", "basis": "source reports aggregate range", "scope": "starting-material recovery"}], ["C-N"], name="potassium-carbonate N-arylation variant", transformation="N-arylation of 24a to 39 under lower-yield K2CO3 conditions", notes=("The source gives only aggregate yield/recovery ranges; exact run count and charges are unreported.",), unresolved=("The aggregate 14–16% yield and 20–25% recovery cannot be assigned to individual runs.",), general_procedure="General Procedure G with K2CO3")
# Drop the redundant point-yield helper outcome; retain the proper range only.
bundle["reactions"][-1]["outcomes"] = bundle["reactions"][-1]["outcomes"][1:]
narylation_run("rxn_p18_24a_36_to_40_Cs2CO3", [18], "36", "40", ("29 mg", 29), ("0.11 mmol", 0.11), ("37 mg", 37, "mg", "0.16 mmol", 0.16), ("4 mg", 4), ("0.02 mmol", 0.02), "DMEDA", 5, 0.04, "Cs2CO3", ("68 mg", 68), ("0.21 mmol", 0.21), 1.0, 25, 63)
narylation_run("rxn_p18_24a_36_to_40_K2CO3", [18], "36", "40", ("20 mg", 20), ("0.072 mmol", 0.072), ("20 mg", 20, "mg", "0.087 mmol", 0.087), ("3 mg", 3), ("0.01 mmol", 0.01), "DMEDA", 4, 0.03, "K2CO3", ("20 mg", 20), ("0.14 mmol", 0.14), 0.72, 5, 19, notes=("The run also recovered 5 mg of 24a (25%).",))
bundle["reactions"][-1]["outcomes"].extend([make_outcome(3, {"kind": "recovery", "targets": ["24a"], "raw": "5 mg, 25% starting material recovered", "value": q(25, "%", raw_text="25%"), "yield_type": "recovery", "product_role": "recovered 24a", "basis": "isolated after column chromatography"}, ev(18)), make_outcome(4, {"kind": "product_amount", "targets": ["24a"], "raw": "5 mg recovered", "value": q(5, "mg", raw_text="5 mg"), "product_role": "recovered 24a", "basis": "isolated recovered mass"}, ev(18))])
narylation_run("rxn_p18_24a_41_to_42_run1", [18], "41", "42", ("21 mg", 21), ("0.076 mmol", 0.076), ("12 μL", 12, "μL", "0.10 mmol", 0.10), ("3 mg", 3), ("0.02 mmol", 0.02), "trans_DMCHDA", 5, 0.03, "Cs2CO3", ("50 mg", 50), ("0.15 mmol", 0.15), 0.8, 15, 56)
narylation_run("rxn_p18_24a_41_to_42_run2", [18], "41", "42", ("23 mg", 23), ("0.083 mmol", 0.083), ("10 μL", 10, "μL", "0.11 mmol", 0.11), ("3 mg", 3), ("0.02 mmol", 0.02), "trans_DMCHDA", 4, 0.03, "Cs2CO3", ("52 mg", 52), ("0.16 mmol", 0.16), 0.8, 19, 63)


def tosyl_run(run, sm_mass, sm_mmol, nah_mass, nah_mmol, tscl_mass, tscl_mmol, dmf_ml, stages, outcomes, notes=()):
    steps = [{"operation": "dissolve 24a in DMF and cool", "temperature": Q("0 °C", 0, "°C"), "time": None, "atmosphere": "nitrogen"}, {"operation": "add sodium hydride (60% in mineral oil) in one portion and stir", "temperature": Q("0 °C", 0, "°C"), "time": Q("15 min", 15, "min")}, {"operation": "add p-toluenesulfonyl chloride in one portion", "temperature": Q("0 °C", 0, "°C"), "time": None}]
    steps += stages
    steps.append({"operation": "quench into saturated ammonium chloride, extract with ethyl acetate, dry, filter, concentrate and chromatograph", "stage": "workup", "temperature": None, "time": None})
    participants = [mass_amt("24a", sm_mass[0], sm_mass[1], "mg", sm_mmol[0], sm_mmol[1]), mass_amt("NaH", nah_mass[0], nah_mass[1], "mg", nah_mmol[0], nah_mmol[1], role="reagent"), mass_amt("pTsCl", tscl_mass[0], tscl_mass[1], "mg", tscl_mmol[0], tscl_mmol[1]), P("DMF", "solvent", volume=Q(f"{dmf_ml} mL", dmf_ml, "mL")), P("43", "product")]
    if run == 1:
        participants += [P("NH4Cl", "workup_agent", stage="workup", volume=Q("5 mL saturated aqueous solution", 5, "mL")), P("EtOAc", "workup_agent", stage="workup", volume=Q("3 × 10 mL", 30, "mL")), P("MgSO4", "workup_agent", stage="workup")]
    return add_reaction(f"rxn_p18_24a_to_43_run{run}", [18], f"page 18 compound 43 Run {run}", participants, steps, outcomes, ["N-S"], name="N-tosylation of lactam 24a", transformation="formation of N-tosyl lactam 43", scale={"raw_text": sm_mmol[0], "number": sm_mmol[1], "unit": "mmol"}, solvent={"raw": f"DMF {dmf_ml} mL", "composition": [mid("DMF")], "total_volume": Q(f"{dmf_ml} mL", dmf_ml, "mL")}, notes=notes)


tosyl_run(1, ("20 mg", 20), ("0.07 mmol", 0.07), ("6 mg", 6), ("0.14 mmol", 0.14), ("28 mg", 28), ("0.14 mmol", 0.14), 0.75, [{"operation": "remove ice bath, warm to room temperature and stir", "temperature": ("0 °C to room temperature", None), "time": Q("2.5 h", 2.5, "h")}], [Y(["43"], 29, raw="9 mg, 29%", basis="column chromatography (10–25% ethyl acetate/hexanes)"), A(["43"], 9, "mg")])
tosyl_run(2, ("18 mg", 18), ("0.069 mmol", 0.069), ("5 mg", 5), ("0.1 mmol", 0.1), ("25 mg", 25), ("0.13 mmol", 0.13), 0.7, [{"operation": "hold after p-toluenesulfonyl chloride addition", "temperature": Q("0 °C", 0, "°C"), "time": Q("1 h", 1, "h")}, {"operation": "remove ice bath, warm to room temperature and stir", "temperature": ("0 °C to room temperature", None), "time": Q("1 h", 1, "h")}], [Y(["43"], 48, raw="14 mg, 48%", ytype="initial isolated yield", role="43 containing unidentified contaminants", basis="first column chromatography"), A(["43"], 14, "mg", role="initial contaminated fraction"), Y(["43"], 18, raw="4 mg, 18% pure 43", ytype="repurified isolated yield", role="pure 43 after additional column", basis="additional column chromatography"), A(["43"], 4, "mg", role="pure repurified 43")], notes=("The initial 14 mg fraction contained unidentified contaminants by 1H NMR; additional chromatography afforded 4 mg of pure 43.",))

# Compound 44 Run 1 begins on page 18 and its workup/outcome continues in page-19 text.
e44 = ev(18) + [f"{DOC}:page19:text"]
add_reaction(
    "rxn_p18_p19_24a_to_44_run1", [18], "page 18 compound 44 Run 1, continued in page-19 text",
    [mass_amt("24a", "15 mg", 15, "mg", "0.054 mmol", 0.054), mass_amt("DMAP", "1 mg", 1, "mg", "0.008 mmol", 0.008, role="catalyst"), vol_amt("Et3N", "0.08 mL", 0.08, "0.5 mmol", 0.5, role="reagent"), mass_amt("Boc2O", "59 mg", 59, "mg", "0.27 mmol", 0.27), P("DCM", "solvent", volume=Q("0.3 mL", 0.3, "mL")), P("44", "product", mass=Q("8 mg", 8, "mg")), P("DCM", "workup_agent", stage="workup", volume=Q("15 mL dilution", 15, "mL")), P("NH4Cl", "workup_agent", stage="workup", volume=Q("5 mL saturated aqueous solution", 5, "mL")), P("brine", "workup_agent", stage="workup", volume=Q("5 mL", 5, "mL")), P("MgSO4", "workup_agent", stage="workup")],
    [{"operation": "dissolve 24a in dichloromethane, add DMAP and cool", "temperature": Q("0 °C", 0, "°C"), "time": None, "atmosphere": "nitrogen"}, {"operation": "add triethylamine dropwise", "temperature": Q("0 °C", 0, "°C"), "time": None}, {"operation": "add di-tert-butyl dicarbonate in one portion", "temperature": Q("0 °C", 0, "°C"), "time": None}, {"operation": "remove ice bath, warm to room temperature and stir", "temperature": ("0 °C to room temperature", None), "time": Q("2 h", 2, "h")}, {"operation": "dilute with dichloromethane; wash saturated ammonium chloride then brine; dry, filter, concentrate and chromatograph", "stage": "workup", "temperature": None, "time": None}],
    [Y(["44"], 40, raw="8 mg, 40%", basis="column chromatography (5–25% ethyl acetate/hexanes)"), A(["44"], 8, "mg")], ["C-N"], name="N-Boc protection of lactam 24a", transformation="formation of N-Boc lactam 44", scale={"raw_text": "24a 0.054 mmol", "number": 0.054, "unit": "mmol"}, solvent={"raw": "CH2Cl2 0.3 mL", "composition": [mid("DCM")], "total_volume": Q("0.3 mL", 0.3, "mL")}, notes=("The reaction setup is on visually inspected page 18; workup and 8 mg/40% outcome were taken from the extracted text continuation on page 19. Page 19 image is outside this batch and is not claimed as visually covered.",), unresolved=("Page 19 image was not part of this six-page visual batch; the cross-page outcome relies on page19:text only.",), general_procedure="compound 44 specific procedure", evidence_override=e44
)


def source_context(context_id, page, region, classification_name, summary, *, pages=None, evidence=None):
    return {"context_id": context_id, "document_id": DOC, "location": {"page": page, "region": region}, "classification": classification_name, "summary": summary, "evidence_ids": evidence or ev(*(pages or [page])), "review_status": "needs_review"}


contexts = [
    source_context("ctx_p13_characterization", 13, "15b, S8, 18a/18b and 27a spectral blocks", "characterization", "NMR/IR/HRMS/TLC data were retained as molecule/source context and not made into synthesis records."),
    source_context("ctx_p14_characterization", 14, "19a/19b, 20a/20b and 28a spectral blocks", "characterization", "Spectral and chromatographic characterization is not an independent synthetic experiment."),
    source_context("ctx_p15_characterization", 15, "21a/21b, 29a, 22a/22b spectral blocks", "characterization", "Characterization blocks and peak assignments were not promoted to reactions."),
    source_context("ctx_p16_stryker_observations", 16, "General Procedure E reagent-source discussion", "optimization_context", "Aged/new commercial and freshly prepared Stryker reagent observations and 24a/31a ratios are reported without complete per-run quantities/yields; they remain context attached to the General Procedure E record."),
    source_context("ctx_p16_31a_characterization", 16, "compound 31a paragraph", "characterization", "The source says the structure of 31a was not confirmed because it could not be separated from 24a; only matching NMR resonances are discussed. No standalone synthesis was invented."),
    source_context("ctx_p17_characterization", 17, "57, 25a, 32a, 26, 56 and 37 spectral blocks", "characterization", "NMR/IR/HRMS/TLC blocks are characterization rather than additional experiments."),
    source_context("ctx_p18_characterization", 18, "38, 39, 40, 42 and 43 spectral blocks", "characterization", "Spectral/analytical data were not treated as synthesis records."),
]
bundle["source_context"].extend(contexts)

(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")

batch_reactions = [r for r in bundle["reactions"] if r["reaction_id"].startswith(batch_prefixes)]
inventory = {
    "batch_id": "allene-stereo-v31-200-20261017",
    "document_id": DOC,
    "visually_inspected_original_pdf_pages": [13, 14, 15, 16, 17, 18],
    "new_reaction_ids": [r["reaction_id"] for r in batch_reactions],
    "new_molecule_ids": [mid("27cis"), mid("trans_DMCHDA")],
    "source_context_ids": [c["context_id"] for c in contexts],
    "limitations": ["All new records remain needs_review and partial.", "No SI association was inferred.", "Complex relative/mixture stereochemistry was retained without speculative isomeric SMILES.", "Aggregate repeated-run yield ranges were not expanded into fictitious runs.", "Page 19 text was used only to complete the cross-page compound-44 Run 1; page 19 image was not claimed as inspected."],
}
(RUN / "raw_model" / "pages_0013_0018_inventory.json").write_text(json.dumps(inventory, ensure_ascii=False, indent=2) + "\n")

page_record_ids = {p: [] for p in range(13, 19)}
for r in batch_reactions:
    for p in range(13, 19):
        if any(e in r["evidence_ids"] for e in ev(p)):
            page_record_ids[p].append(r["reaction_id"])

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 13, "state": "attempted", "regions": [
            {"region": "top-left continuation of compound 15b", "finding": "spectral characterization continuation only; no new synthesis", "evidence_ids": ev(13)},
            {"region": "compound S8 Runs 1–2", "finding": "two distinct N-methylation runs transcribed with quantities, staged 0 °C additions/warm-up, yields and trans/cis compositions", "evidence_ids": ev(13)},
            {"region": "General Procedure B and 18a/18b/27a preparations", "finding": "five ketal-hydrolysis runs recorded; cis-enriched fractions and unreported independent cis yield kept distinct", "evidence_ids": ev(13)},
            {"region": "bottom-right start of compound 19a Run 1", "finding": "cross-page run merged with page 14 rather than duplicated", "evidence_ids": ev(13)},
        ], "record_ids": page_record_ids[13], "unresolved": ["Aggregate repeated-run ranges for S8 and 18a/b were not expanded into unreported runs.", "The isolated cis ketone corresponding to 27a is unnumbered and lacks complete absolute configuration."]},
        {"page": 14, "state": "attempted", "regions": [
            {"region": "19a/19b General Procedure B runs", "finding": "four hydrolysis runs/fractionation outcomes linked, including page-13 continuation", "evidence_ids": ev(14)},
            {"region": "General Procedure C and 20a/20b Runs 1–3", "finding": "three ethynylmagnesium bromide addition runs recorded with separate chromatographic fractions and dr/trans-cis results", "evidence_ids": ev(14)},
            {"region": "compound 28a procedure", "finding": "specific 27a-to-28a run recorded with 62% yield, 7:1 trans/cis, 1:1 dr and ~7% residual ketone by NMR", "evidence_ids": ev(14)},
            {"region": "spectral blocks", "finding": "characterization only, not independent syntheses", "evidence_ids": ev(14)},
        ], "record_ids": page_record_ids[14], "unresolved": ["Exact absolute configurations at newly formed propargylic centers are not assigned.", "The cis-enriched 21b characterization sample has both a reported 2:1 cis/trans fraction description and ~56% trans by NMR; both are preserved."]},
        {"page": 15, "state": "attempted", "regions": [
            {"region": "21a/21b procedures", "finding": "cross-page 21a run and cis-enriched fractionation run recorded with yields and NMR-based dr", "evidence_ids": ev(15)},
            {"region": "compound 29a direct acetylation and one-sentence General Procedure D alternative", "finding": "direct 54% acetylation and separate source-stated 30% alternative retained; the latter is flagged because General Procedure D is a pivaloylation protocol", "evidence_ids": ev(15)},
            {"region": "General Procedure D, compounds 22a/22b", "finding": "four pivaloylation runs recorded with every isolated fraction and product-specific yield/ratio", "evidence_ids": ev(15)},
            {"region": "bottom-right start of 23a Run 1", "finding": "cross-page run merged with page 16", "evidence_ids": ev(15)},
        ], "record_ids": page_record_ids[15], "unresolved": ["The source-stated 'Following General Procedure D' preparation of acetate 29a is chemically inconsistent with the printed pivaloylation procedure; no reagent was guessed or corrected.", "Several complex stereoisomer mixtures remain without isomeric SMILES."]},
        {"page": 16, "state": "attempted", "regions": [
            {"region": "23a/23b General Procedure D runs", "finding": "cross-page Run 1 plus Runs 2–3 recorded with trans, cis and mixed fraction yields", "evidence_ids": ev(16)},
            {"region": "General Procedure E and compound 24a", "finding": "specific 22a-to-24a run recorded with ordered degassing, cooling, reagent addition, quench and 33% yield", "evidence_ids": ev(16)},
            {"region": "Stryker-reagent source/temperature discussion", "finding": "optimization observations retained as source context because individual charges/yields are incomplete", "evidence_ids": ev(16)},
            {"region": "compound 30a and 57 procedures", "finding": "two specific allene-formation/reduction runs recorded", "evidence_ids": ev(16)},
            {"region": "compound 31a paragraph", "finding": "unconfirmed inseparable byproduct characterization only; no standalone synthesis invented", "evidence_ids": ev(16)},
        ], "record_ids": page_record_ids[16], "unresolved": ["Individual yields for reagent-source observations are unreported.", "Compound 31a structure is explicitly unconfirmed by the authors."]},
        {"page": 17, "state": "attempted", "regions": [
            {"region": "compound 25a preparation", "finding": "23a reaction recorded with 25a and 32a as separately targeted outcomes, including 4:1 ratio", "evidence_ids": ev(17)},
            {"region": "General Procedure F, compounds 26 and 56", "finding": "two 25a-to-26 runs and one 32a-to-56 run recorded with 0 °C/45 min stages and dr", "evidence_ids": ev(17)},
            {"region": "General Procedure G and compound 37", "finding": "Cs2CO3 run and lower-yield K2CO3 variant recorded separately; ligand prose/detail conflict retained", "evidence_ids": ev(17)},
            {"region": "bottom-right start of compound 38 procedure", "finding": "cross-page run merged with page 18", "evidence_ids": ev(17)},
            {"region": "spectral blocks", "finding": "characterization only", "evidence_ids": ev(17)},
        ], "record_ids": page_record_ids[17], "unresolved": ["The compound-37 alternative run names DMEDA in prose but trans-N,N′-dimethylcyclohexane-1,2-diamine in the explicit charge; no correction was made.", "Aggregate 25a/32a yield ranges from other runs were not expanded."]},
        {"page": 18, "state": "attempted", "regions": [
            {"region": "compounds 38–40 N-arylation procedures", "finding": "38 cross-page run, 39 Cs/K2CO3 variants, and 40 Cs/K2CO3 runs recorded separately with yield/recovery targets", "evidence_ids": ev(18)},
            {"region": "compound 42 Runs 1–2", "finding": "two independent 2-iodothiophene coupling runs recorded", "evidence_ids": ev(18)},
            {"region": "compound 43 Runs 1–2", "finding": "two N-tosylation runs recorded; Run 2 preserves initial contaminated and repurified yields", "evidence_ids": ev(18)},
            {"region": "compound 44 Run 1 start", "finding": "reaction setup read from page-18 image; page-19 text used only to complete workup and 40% outcome", "evidence_ids": ev(18)},
            {"region": "spectral blocks", "finding": "characterization only, not additional reactions", "evidence_ids": ev(18)},
        ], "record_ids": page_record_ids[18], "unresolved": ["K2CO3 compound-39 data are aggregate ranges with unspecified individual runs.", "Compound 44 Run 1 outcome continues on page 19 text; page 19 image was not inspected or claimed covered."]},
    ],
}
(RUN / "coverage" / "pages_0013_0018.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")

print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]), "source_context": len(bundle["source_context"]), "new_records": len(batch_reactions)}, ensure_ascii=False))
