#!/usr/bin/env python3
"""Deterministic v3.1 transcription for original PDF pages 19--24.

This augments the current document bundle, replaces only records created by
this page batch on reruns, and keeps every machine-created reaction partial and
needs_review.  It intentionally leaves complex structures unresolved instead
of inventing stereo-complete SMILES.
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


def na(evidence):
    return fact(evidence=evidence, method="not_applicable", status="not_applicable", verification="missing")


def q(number, unit, *, qualifier=None, raw_text=None):
    v = {"number": number, "unit": unit}
    qualifier = {"approximately": "approx", "less_than": "lt"}.get(qualifier, qualifier)
    if qualifier in {"within", "start of qualitative ramp"}:
        qualifier = None
    if qualifier is not None:
        v["qualifier"] = qualifier
    if raw_text is not None:
        v["raw_text"] = raw_text
    return v


def Q(raw, number, unit, qualifier=None):
    return (raw, q(number, unit, qualifier=qualifier, raw_text=raw))


bundle = json.loads((RUN / "submitted.json").read_text())
dataset = json.loads((RUN / "dataset.json").read_text())
bundle["documents"] = copy.deepcopy(dataset["documents"])

# Idempotence: replace only newly created page-19--24 records and context.
batch_prefixes = tuple(f"{DOC}:rxn:rxn_p{p}_" for p in range(19, 25))
bundle["reactions"] = [r for r in bundle["reactions"] if not r["reaction_id"].startswith(batch_prefixes)]
context_prefixes = tuple(f"ctx_p{p}_" for p in range(19, 25))
bundle["source_context"] = [c for c in bundle.get("source_context", []) if not c["context_id"].startswith(context_prefixes)]

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
    basis = "Source label/drawing reviewed conservatively; unreported configuration was not converted to achiral or racemic."
    st["chirality_status"] = inferred(label, chirality, evidence, basis)
    st["geometry_status"] = inferred(label, geometry, evidence, basis)
    st["representation_status"] = inferred(label, representation, evidence, basis)
    st["isomeric_smiles"] = missing(evidence)
    st["identity_key"] = missing(evidence)
    m["identity_note"] = inferred(label, note or "Source entity retained without speculative complete structure encoding.", evidence, basis)
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


# Source products/intermediates newly encountered in this batch.
add_molecule("S9", "compound S9", 21, chirality="partially_specified", geometry="specified", assignment="(3aS*,9bS*) relative stereochemistry; absolute configuration unreported")
add_molecule("S10", "compound S10", 21, chirality="mixture", geometry="specified", assignment="relative fused-ring stereochemistry retained; carbon adjacent to lactam carbonyl unresolved", mixture="prepared from 58 reported as a 9:1 diastereomeric mixture")
add_molecule("S12", "compound S12", 21, chirality="unspecified", geometry="mixture", assignment="E/Z order is source-reported; complete structure is not encoded", mixture="E/Z = 2.5:7.5")
add_molecule("S13", "compound S13", 21, chirality="mixture", geometry="mixture", assignment="trans/cis relationship source-reported; pure trans fraction only used for characterization", mixture="trans/cis 8.5:2.5")
add_molecule("S14", "compound S14 (pure trans)", 22, chirality="partially_specified", geometry="specified", assignment="pure trans relationship reported; absolute configuration unreported")
add_molecule("S15_alcohol", "crude propargyl alcohol intermediate preceding acetate S15", 22, chirality="unspecified", geometry="specified", assignment="trans lactone precursor retained; new propargylic configuration unreported")
add_molecule("S15", "compound S15 (propargyl acetate)", 22, chirality="partially_specified", geometry="specified", assignment="relative lactone stereochemistry retained; propargylic configuration unreported")
add_molecule("S16_michael", "crude dimethylamine Michael adduct en route to S16", 22, chirality="unresolved", geometry="unresolved", assignment="temporary crude intermediate named only by role")
add_molecule("S16_crude_allene", "crude allene-yne intermediate en route to S16", 22, chirality="partially_specified", geometry="specified", assignment="crude intermediate identity is source-described but not fully characterized")
add_molecule("S16", "compound S16", 22, chirality="partially_specified", geometry="specified", assignment="relative lactone stereochemistry retained; terminal-CH2 allene is not axially stereogenic")

# Named reagents/materials not present in the registry.  No CAS values are guessed.
for key, label, page in [
    ("polyPPh3", "polymer-bound triphenylphosphine (3 mmol/g)", 19),
    ("pentane", "pentane/pentanes", 19),
    ("TfOH", "trifluoromethanesulfonic acid", 21),
    ("Na2SO4", "sodium sulfate", 22),
    ("DMP", "Dess–Martin periodinane", 22),
    ("Celite", "Celite", 22),
    ("CeCl3", "anhydrous cerium(III) chloride beads", 22),
    ("dimethylamine", "dimethylamine (2 M in THF)", 22),
    ("PBS", "phosphate-buffered saline", 23),
    ("KOH", "potassium hydroxide", 23),
    ("hexamethylbenzene", "hexamethylbenzene", 23),
]:
    add_molecule(key, label, page, chirality="unresolved", geometry="unresolved", representation="unresolved", note="Named reagent/material; structure and CAS were not inferred for this extraction.")

for key in ("24a", "25a", "26", "30a", "37", "38", "39", "40", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58", "59", "60", "S11", "S17"):
    pages = [p for p in range(19, 25) if key in {
        19: {"24a", "25a", "26", "30a", "44", "45", "46", "47", "48"},
        20: {"37", "38", "39", "40", "42", "43", "44", "49", "50", "51", "52", "53", "54", "55"},
        21: {"44", "45", "47", "55", "56", "57", "58", "59", "60"},
        22: {"S17"},
        23: {"45", "48", "49", "51", "52", "54", "60", "S11", "S17"},
        24: {"51", "S17"},
    }[p]]
    if pages:
        add_occurrence(key, f"compound {key}", pages)


def charge_fact(value, evidence):
    if value is None:
        return missing(evidence)
    raw, normalized = value if isinstance(value, tuple) else (json.dumps(value, ensure_ascii=False), value)
    return fact(raw=raw, normalized=normalized, evidence=evidence, method="original_extraction")


def P(key, role, **kwargs):
    if role == "purification_solvent":
        role = "purification_agent"
    return {"key": key, "role": role, **kwargs}


def make_participant(index, spec, evidence, reference):
    p = new_participant(f"p{index}", mid(spec["key"]), spec["role"], evidence, stage=spec.get("stage", "reaction"))
    p["product_role"] = extracted(spec.get("product_role", "product"), spec.get("product_role", "product"), evidence) if spec["role"] == "product" else na(evidence)
    for field_name in ("mass", "amount", "volume", "equivalents", "mol_percent"):
        p["charges"][field_name] = charge_fact(spec.get(field_name), evidence)
        if field_name == "mol_percent" and isinstance(p["charges"][field_name].get("normalized_value"), dict):
            p["charges"][field_name]["normalized_value"]["unit"] = "mol%"
    p["equivalence_reference"] = extracted("relative to limiting substrate", reference, evidence) if spec.get("equivalents") is not None or spec.get("mol_percent") is not None else missing(evidence)
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


LABELS = {"C-H": 1, "C-C": 2, "C-O": 3, "C-S": 4, "C-N": 5, "DB-C-C": 13, "DB-C-O": 14}
ATOM_PAIRS = {"C-H": ("C", "H"), "C-C": ("C", "C"), "C-O": ("C", "O"), "C-N": ("C", "N"), "C-S": ("C", "S"), "DB-C-C": ("C", "C"), "DB-C-O": ("C", "O")}


def set_classification(r, labels, evidence, rationale):
    c = r["classification"]
    c["mapped_reaction_smiles"] = missing(evidence)
    changes, values = [], []
    for i, code in enumerate(labels, 1):
        a, b = ATOM_PAIRS[code]
        change_id = f"b{i}"
        change = {"change_id": change_id, "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}], "before_order": 0, "after_order": 2 if code.startswith("DB-") else 1, "change_type": "formed"}
        changes.append(inferred(rationale, change, evidence, "Bond class inferred from the source transformation; atom mapping is unavailable."))
        values.append({"category_number": LABELS[code], "code": code, "supporting_change_ids": [change_id]})
    c["bond_changes"] = changes
    c["labels"] = inferred(rationale, values, evidence, "Multi-label classification inferred under bond-formation-20261016.")
    c["coverage"] = inferred(rationale, "partial", evidence, "Only source-supported newly formed bond classes are asserted.")
    c["limitations"] = ["No atom-mapped reaction SMILES.", "Complete atom-level bond-change inventory remains unresolved."]
    c["model_votes"] = []


def make_outcome(index, spec, evidence):
    o = new_outcome(f"o{index}", spec["kind"], [mid(k) for k in spec["targets"]], evidence)
    o["product_role"] = extracted(spec["role"], spec["role"], evidence) if spec.get("role") else na(evidence)
    if spec.get("status"):
        o["value"] = missing(evidence, spec["status"])
    else:
        o["value"] = extracted(spec.get("raw", json.dumps(spec["value"], ensure_ascii=False)), spec["value"], evidence)
    o["yield_type"] = extracted(spec["yield_type"], spec["yield_type"], evidence) if spec.get("yield_type") else na(evidence)
    o["measurement_basis"] = extracted(spec["basis"], spec["basis"], evidence) if spec.get("basis") else missing(evidence)
    o["scope"] = extracted(spec["scope"], spec["scope"], evidence) if spec.get("scope") else missing(evidence)
    return o


def add_reaction(suffix, pages, location, participants, steps, outcomes, labels, *, name, transformation, category="other", scale=None, solvent=None, notes=(), unresolved=(), general=None, evidence_override=None):
    evidence = evidence_override or ev(*pages)
    r = new_reaction(f"{DOC}:rxn:{suffix}", DOC, evidence)
    r["record_type"] = "specific_experiment"
    r["experiment_categories"] = extracted(category, [category], evidence)
    r["location"] = extracted(location, {"page": pages[0], "region": location}, evidence)
    r["general_scheme_ids"] = []
    r["reaction_name"] = extracted(name, name, evidence)
    r["transformation_type"] = extracted(transformation, transformation, evidence)
    first = next((p for p in participants if p["role"] == "reactant"), participants[0])
    reference = mid(first["key"])
    r["participants"] = [make_participant(i, p, evidence, reference) for i, p in enumerate(participants, 1)]
    if scale:
        sv = dict(scale); sv["reference_molecule_id"] = reference
        r["reaction_scale"] = extracted(scale.get("raw_text", json.dumps(sv, ensure_ascii=False)), sv, evidence)
    else:
        r["reaction_scale"] = missing(evidence)
    solvent = solvent or {}
    r["solvent_conditions"]["composition"] = extracted(solvent["raw"], solvent["composition"], evidence) if solvent.get("composition") else missing(evidence)
    r["solvent_conditions"]["mixture_ratio"] = extracted(solvent["ratio"], solvent["ratio"], evidence) if solvent.get("ratio") else missing(evidence)
    r["solvent_conditions"]["total_volume"] = charge_fact(solvent.get("total_volume"), evidence)
    r["solvent_conditions"]["concentration"] = charge_fact(solvent.get("concentration"), evidence)
    r["solvent_conditions"]["concentration_basis"] = extracted(solvent["basis"], solvent["basis"], evidence) if solvent.get("basis") else missing(evidence)
    r["steps"] = [make_step(i, s, evidence) for i, s in enumerate(steps, 1)]
    for key in r["special_conditions"]:
        r["special_conditions"][key] = na(evidence)
    r["procedure"]["original_text"] = extracted(location, location, evidence)
    r["procedure"]["general_procedure_reference"] = extracted(general, general, evidence) if general else missing(evidence)
    r["procedure"]["modifications"] = missing(evidence)
    r["procedure"]["addition_order"] = extracted("ordered as transcribed", [s["operation"]["normalized_value"] for s in r["steps"]], evidence)
    r["procedure"]["pretreatment"] = missing(evidence)
    r["outcomes"] = [make_outcome(i, o, evidence) for i, o in enumerate(outcomes, 1)]
    set_classification(r, labels, evidence, transformation)
    r["notes"] = [extracted(n, n, evidence) for n in notes]
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [], "unresolved": ["Machine record remains needs_review; schema validation is not chemistry verification.", "Stereo-complete structures and atom mapping remain unresolved.", *unresolved]}
    r["completeness"] = "partial"
    bundle["reactions"].append(r)
    return r


def Y(target, pct, *, raw=None, basis="isolated after column chromatography", ytype="isolated yield", scope=None):
    raw = raw or f"{pct}%"
    return {"kind": "yield", "targets": [target], "raw": raw, "value": q(pct, "%", raw_text=raw), "yield_type": ytype, "role": target, "basis": basis, "scope": scope}


def A(target, number, unit, *, raw=None, role=None, basis="reported isolated mass"):
    raw = raw or f"{number} {unit}"
    return {"kind": "product_amount", "targets": [target], "raw": raw, "value": q(number, unit, raw_text=raw), "role": role or target, "basis": basis}


def S(targets, raw, *, basis, scope, value=None, role="product mixture"):
    return {"kind": "selectivity", "targets": targets, "raw": raw, "value": raw if value is None else value, "role": role, "basis": basis, "scope": scope}


def mass_amt(key, mass_raw, mass, mass_unit, amt_raw, amt, role="reactant", **kw):
    return P(key, role, mass=Q(mass_raw, mass, mass_unit), amount=Q(amt_raw, amt, "mmol"), **kw)


def vol_amt(key, vol_raw, vol, vol_unit, amt_raw, amt, role="reagent", **kw):
    return P(key, role, volume=Q(vol_raw, vol, vol_unit), amount=Q(amt_raw, amt, "mmol"), **kw)


def apkr_run(suffix, page, substrate, product, sm_mass, sm_mmol, cat, solvent_ml, product_mass, yld, *, procedure="I", time=None, color="solid", note=(), unresolved=()):
    cat_kind, cat_v, cat_amt, cat_mmol = cat
    ps = [mass_amt(substrate, f"{sm_mass} mg", sm_mass, "mg", sm_mmol[0], sm_mmol[1]), P("Rhcat", "catalyst"), P("toluene", "solvent", volume=Q(f"{solvent_ml} mL", solvent_ml, "mL")), P("CO", "reagent"), P(product, "product", mass=Q(f"{product_mass} mg", product_mass, "mg"))]
    rh = ps[1]
    if cat_kind == "mass":
        rh["mass"] = Q(f"{cat_v} mg", cat_v, "mg"); rh["amount"] = Q(cat_amt, cat_mmol, "mmol")
    else:
        rh["volume"] = Q(f"{cat_v} mL catalyst solution", cat_v, "mL")
        if cat_mmol is not None:
            rh["amount"] = Q(cat_amt, cat_mmol, "mmol")
    if procedure == "H":
        ps.append(P("polyPPh3", "workup_agent", stage="workup"))
        steps = [
            {"operation": "place rhodium catalyst/toluene under vacuum for 2–3 s and refill with CO three times; preheat", "temperature": Q("110 °C", 110, "°C"), "time": {"raw_text": "2–3 s per vacuum cycle", "minimum": 2, "maximum": 3, "unit": "s"}, "atmosphere": "carbon monoxide"},
            {"operation": "add allene-yne solution dropwise by syringe pump", "temperature": Q("110 °C", 110, "°C"), "time": Q("1 h addition", 1, "h"), "atmosphere": "carbon monoxide"},
            {"operation": "hold after addition", "temperature": Q("110 °C", 110, "°C"), "time": Q(time or "30 min", 30 if not time else time[1], "min" if not time else time[2]), "atmosphere": "carbon monoxide"},
            {"operation": "cool; treat with polymer-bound triphenylphosphine, filter, concentrate and chromatograph", "stage": "workup", "temperature": ("room temperature", "r.t."), "time": Q("4 h polymer treatment (General Procedure H)", 4, "h")},
        ]
    else:
        steps = [
            {"operation": "place rhodium catalyst/toluene under vacuum for 2–3 s and refill with CO three times; preheat", "temperature": Q("110 °C", 110, "°C"), "time": {"raw_text": "2–3 s per vacuum cycle", "minimum": 2, "maximum": 3, "unit": "s"}, "atmosphere": "carbon monoxide"},
            {"operation": "add allene-yne solution in a single portion and heat", "temperature": Q("110 °C", 110, "°C"), "time": Q(time[0], time[1], time[2]) if time else Q("1 h (General Procedure I)", 1, "h"), "atmosphere": "carbon monoxide"},
            {"operation": "concentrate and purify by silica-gel column chromatography", "stage": "purification", "temperature": ("room temperature", "r.t."), "time": None},
        ]
    return add_reaction(suffix, [page], f"page {page} compound {product} preparation", ps, steps, [Y(product, yld, raw=f"{product_mass} mg, {yld}%"), A(product, product_mass, "mg")], ["C-C", "DB-C-O"], name="Rh(I)-catalysed allenic Pauson–Khand reaction", transformation=f"intramolecular APKR of {substrate} to fused lactam {product}", category="substrate_scope", scale={"raw_text": sm_mmol[0], "number": sm_mmol[1], "unit": "mmol"}, solvent={"raw": f"toluene {solvent_ml} mL", "composition": [mid("toluene")], "total_volume": Q(f"{solvent_ml} mL", solvent_ml, "mL")}, notes=note, unresolved=unresolved, general=f"General Procedure {procedure}")


# Page 19: finish the existing cross-page 44 Run 1 with visually inspected evidence.
r44 = next(r for r in bundle["reactions"] if r["reaction_id"].endswith("rxn_p18_p19_24a_to_44_run1"))
r44["evidence_ids"] = ev(18, 19)
for obj in [r44, *r44["participants"], *r44["steps"], *r44["outcomes"], *r44["notes"]]:
    for node in ([obj] if isinstance(obj, dict) else []):
        if "evidence_ids" in node:
            node["evidence_ids"] = ev(18, 19)
r44["notes"] = [extracted("The setup is on page 18 and the workup/yield continue on visually inspected page 19.", evidence=ev(18, 19))]
r44["review"]["unresolved"] = [x for x in r44["review"]["unresolved"] if "Page 19 image" not in x]

add_reaction("rxn_p19_24a_to_44_run2", [19], "page 19 compound 44 Run 2", [mass_amt("24a", "23 mg", 23, "mg", "0.083 mmol", 0.083), mass_amt("DMAP", "1 mg", 1, "mg", "0.008 mmol", 0.008, role="catalyst"), vol_amt("Et3N", "0.12 mL", 0.12, "mL", "0.83 mmol", 0.83), mass_amt("Boc2O", "90 mg", 90, "mg", "0.41 mmol", 0.41), P("DCM", "solvent", volume=Q("0.4 mL", 0.4, "mL")), P("44", "product", mass=Q("17 mg", 17, "mg"))], [{"operation": "perform N-Boc protection according to the compound-44 procedure", "temperature": Q("0 °C to room temperature", 0, "°C", "start of qualitative ramp"), "time": Q("2 h (preceding specific procedure)", 2, "h"), "atmosphere": "nitrogen"}, {"operation": "work up and purify by column chromatography", "stage": "workup", "temperature": None, "time": None}], [Y("44", 54, raw="17 mg, 54%"), A("44", 17, "mg")], ["C-N"], name="N-Boc protection of 24a", transformation="formation of N-Boc allene-yne 44", scale={"raw_text": "0.083 mmol 24a", "number": 0.083, "unit": "mmol"}, solvent={"raw": "CH2Cl2 0.4 mL", "composition": [mid("DCM")], "total_volume": Q("0.4 mL", 0.4, "mL")}, general="compound 44 specific procedure")

apkr_run("rxn_p19_24a_to_45_run1", 19, "24a", "45", 13, ("0.044 mmol", 0.044), ("mass", 2, "0.004 mmol", 0.004), 4.5, 9, 64, procedure="H", note=("Catalyst solution used 3.4 mL toluene and substrate solution used 1.1 mL; total 4.5 mL is retained.",))
apkr_run("rxn_p19_24a_to_45_run2", 19, "24a", "45", 21, ("0.076 mmol", 0.076), ("mass", 3, "0.008 mmol", 0.008), 7.6, 19, 79, procedure="H", note=("Compound 45 was crystallized by EtOAc/pentane vapor diffusion and confirmed by X-ray crystallography; crystallization is characterization, not a separate synthesis.",))
apkr_run("rxn_p19_24a_to_45_run3", 19, "24a", "45", 30, ("0.11 mmol", 0.11), ("mass", 4, "0.01 mmol", 0.01), 11.0, 25, 75, procedure="I", time=("1 h", 1, "h"), note=("Catalyst solution used 7 mL toluene and substrate solution used 4 mL; total 11 mL is retained.",))
apkr_run("rxn_p19_25a_to_46_run1", 19, "25a", "46", 20, ("0.056 mmol", 0.056), ("mass", 2, "0.006 mmol", 0.006), 5.7, 11, 50, procedure="H", time=("3 h", 180, "min"), note=("TLC showed consumption after 3 h; whether this is total or post-addition time is not explicitly resolved.",), unresolved=("The 3 h timing basis relative to the 1 h slow addition is ambiguous." ,))
apkr_run("rxn_p19_25a_to_46_run2", 19, "25a", "46", 25, ("0.078 mmol", 0.078), ("mass", 3, "0.008 mmol", 0.008), 8.0, 11, 41, procedure="H")
apkr_run("rxn_p19_26_to_47_run1", 19, "26", "47", 12, ("0.060 mmol", 0.060), ("mass", 3, "0.006 mmol", 0.006), 6.0, 9, 64, procedure="H")
apkr_run("rxn_p19_26_to_47_run2", 19, "26", "47", 17, ("0.085 mmol", 0.085), ("mass", 4, "0.008 mmol", 0.008), 8.6, 8, 42, procedure="H")
apkr_run("rxn_p19_30a_to_48", 19, "30a", "48", 30, ("0.1 mmol", 0.1), ("mass", 4, "0.009 mmol", 0.009), 8.6, 24, 75, procedure="H")

# Page 20: individual APKR runs.  Catalyst solution volumes are retained as volumes where printed.
apkr_run("rxn_p20_37_to_49_run1", 20, "37", "49", 11, ("0.031 mmol", 0.031), ("solution", 1.2, "0.003 mmol", 0.003), 3.1, 8, 67, procedure="I")
apkr_run("rxn_p20_37_to_49_run2", 20, "37", "49", 3, ("0.008 mmol", 0.008), ("solution", 0.5, "amount unreported", None), 0.5, 2, 66, procedure="I", unresolved=("Catalyst mmol is not reported for Run 2.",))
apkr_run("rxn_p20_38_to_50", 20, "38", "50", 17, ("0.045 mmol", 0.045), ("mass", 2, "0.0045 mmol", 0.0045), 4.5, 13, 72, procedure="H")
apkr_run("rxn_p20_39_to_51_run1", 20, "39", "51", 17, ("0.04 mmol", 0.04), ("mass", 2, "0.004 mmol", 0.004), 4.0, 13, 72, procedure="H")
apkr_run("rxn_p20_39_to_51_run2", 20, "39", "51", 6, ("0.01 mmol", 0.01), ("solution", 0.5, "amount unreported", None), 0.9, 4, 67, procedure="H", unresolved=("Catalyst mmol is not reported for Run 2.",))
apkr_run("rxn_p20_40_to_52_run1", 20, "40", "52", 6, ("0.2 mmol", 0.2), ("solution", 0.5, "0.001 mmol", 0.001), 1.3, 5, 76, procedure="I", note=("The printed 6 mg / 0.2 mmol substrate charge is preserved without correction despite its internal inconsistency.",), unresolved=("Compound 40 Run 1 reports 6 mg and 0.2 mmol; the discrepancy requires expert/source review." ,))
apkr_run("rxn_p20_40_to_52_run2", 20, "40", "52", 25, ("0.065 mmol", 0.065), ("mass", 3, "0.006 mmol", 0.006), 6.6, 20, 77, procedure="I")
apkr_run("rxn_p20_42_to_53_run1", 20, "42", "53", 12, ("0.033 mmol", 0.033), ("solution", 1.3, "0.003 mmol", 0.003), 3.3, 7, 54, procedure="I")
apkr_run("rxn_p20_42_to_53_run2", 20, "42", "53", 18, ("0.050 mmol", 0.050), ("mass", 2, "0.005 mmol", 0.005), 5.0, 15, 79, procedure="I")
apkr_run("rxn_p20_43_to_54_run1", 20, "43", "54", 4, ("0.009 mmol", 0.009), ("solution", 0.5, "0.001 mmol", 0.001), 1.0, 3, 75, procedure="I")
apkr_run("rxn_p20_43_to_54_run2", 20, "43", "54", 9, ("0.02 mmol", 0.02), ("mass", 1, "0.002 mmol", 0.002), 2.1, 6, 67, procedure="I")

# Page 21: final APKR examples and N-acetylation.
apkr_run("rxn_p21_44_to_55_run1", 21, "44", "55", 7, ("0.018 mmol", 0.018), ("mass", 1, "0.002 mmol", 0.002), 1.9, 3, 40, procedure="H")
apkr_run("rxn_p21_44_to_55_run2", 21, "44", "55", 16, ("0.042 mmol", 0.042), ("mass", 2, "0.004 mmol", 0.004), 4.2, 7, 41, procedure="I")
r58 = apkr_run("rxn_p21_56_to_58", 21, "56", "58", 12, ("0.059 mmol", 0.059), ("mass", 2, "0.005 mmol", 0.005), 5.0, 11, 78, procedure="H", note=("Starting allene-yne 56 was a 9:1 dr mixture; isolated 58 was 7.3:1 dr.", "The stereochemistry adjacent to the lactam carbonyl was not assigned."))
r58["outcomes"].extend([make_outcome(3, S(["56"], "9:1 dr", basis="source description of charged starting material", scope="starting-material diastereomeric composition"), ev(21)), make_outcome(4, S(["58"], "7.3:1 dr", basis="1H NMR integration of lactam-methine resonances at 4.56 and 4.51 ppm", scope="isolated-product diastereomeric ratio"), ev(21))])
apkr_run("rxn_p21_57_to_59", 21, "57", "59", 17, ("0.061 mmol", 0.061), ("mass", 3, "0.0061 mmol", 0.0061), 6.0, 14, 77, procedure="I", note=("The stereochemistry adjacent to the lactam carbonyl was not assigned.",))


def acetyl_run(suffix, substrate, product, sm_mass, sm_mmol, et3n_vol, et3n_mmol, aco_vol, aco_mmol, dcm, hours, out_mass, yld, *, notes=(), unresolved=()):
    ps = [mass_amt(substrate, f"{sm_mass} mg", sm_mass, "mg", f"{sm_mmol} mmol", sm_mmol), mass_amt("DMAP", "1 mg", 1, "mg", "0.008 mmol", 0.008, role="catalyst"), vol_amt("Et3N", f"{et3n_vol} mL", et3n_vol, "mL", f"{et3n_mmol} mmol", et3n_mmol), P("Ac2O", "reactant", volume=Q(f"{aco_vol} mL", aco_vol, "mL")), P("DCM", "solvent", volume=Q(f"{dcm} mL", dcm, "mL")), P(product, "product", mass=Q(f"{out_mass} mg", out_mass, "mg"))]
    if aco_mmol is not None:
        ps[3]["amount"] = Q(f"{aco_mmol} mmol", aco_mmol, "mmol")
    return add_reaction(suffix, [21], f"page 21 compound {product} preparation", ps, [{"operation": "dissolve lactam in DCM, add DMAP and cool", "temperature": Q("0 °C", 0, "°C"), "time": None, "atmosphere": "nitrogen"}, {"operation": "add triethylamine then acetic anhydride dropwise", "temperature": Q("0 °C", 0, "°C"), "time": None, "atmosphere": "nitrogen"}, {"operation": "allow to warm to room temperature and stir", "temperature": ("0 °C to room temperature", None), "time": Q(f"{hours} h", hours, "h"), "atmosphere": "nitrogen"}, {"operation": "work up and purify by silica-gel column chromatography", "stage": "workup", "temperature": None, "time": None}], [Y(product, yld, raw=f"{out_mass} mg, {yld}%"), A(product, out_mass, "mg")], ["C-N"], name="N-acetylation of secondary lactam", transformation=f"formation of N-acetyl lactam {product} from {substrate}", scale={"raw_text": f"{sm_mmol} mmol {substrate}", "number": sm_mmol, "unit": "mmol"}, solvent={"raw": f"DCM {dcm} mL", "composition": [mid("DCM")], "total_volume": Q(f"{dcm} mL", dcm, "mL")}, notes=notes, unresolved=unresolved, general="General Procedure J")


acetyl_run("rxn_p21_47_to_60", "47", "60", 8, 0.03, 0.05, 0.3, 0.02, 0.2, 0.4, 2, 7, 78)
acetyl_run("rxn_p21_45_to_S9", "45", "S9", 9, 0.03, 0.04, 0.3, 0.02, 0.2, 0.2, 20, 6, 60)
acetyl_run("rxn_p21_58_to_S10", "58", "S10", 11, 0.05, 0.07, 0.5, 0.02, None, 0.5, 6.5, 6, 46, notes=("Starting 58 was reported as a 9:1 diastereomeric mixture; product dr was assessed by lactam-methine integration but no numeric product ratio is stated.",), unresolved=("Acetic anhydride amount in mmol and numeric S10 product dr are unreported.",))

# Pages 21--22: lactone sequence.  S13 spans the page boundary.
e13 = ev(21, 22)
r13 = add_reaction("rxn_p21_p22_S12_to_S13", [21, 22], "pages 21–22 compound S13 preparation", [mass_amt("S12", "1.2 g", 1.2, "g", "2.9 mmol", 2.9), mass_amt("12", "757 mg", 757, "mg", "5.8 mmol", 5.8), P("toluene", "solvent", volume=Q("7.8 mL", 7.8, "mL")), mass_amt("TfOH", "43.6 mg", 43.6, "mg", "0.29 mmol", 0.29, role="catalyst"), P("NH4Cl", "workup_agent", stage="workup"), P("NH4OH", "workup_agent", stage="workup"), P("Et2O", "workup_agent", stage="workup", volume=Q("3 × 60 mL", 180, "mL")), P("brine", "workup_agent", stage="workup", volume=Q("2 × 30 mL", 60, "mL")), P("Na2SO4", "workup_agent", stage="workup"), P("silica", "purification_agent", stage="purification", mass=Q("100 g SNAP column", 100, "g")), P("pentane", "purification_solvent", stage="purification"), P("S13", "product", mass=Q("0.375 g", 0.375, "g"))], [{"operation": "charge S12, aldehyde 12 and toluene; cool", "temperature": Q("0 °C", 0, "°C"), "time": None, "atmosphere": "nitrogen"}, {"operation": "add trifluoromethanesulfonic acid and maintain", "temperature": Q("0 °C", 0, "°C"), "time": Q("12 h", 12, "h"), "atmosphere": "nitrogen"}, {"operation": "dilute with 9:1 aqueous NH4Cl/NH4OH (58 mL); extract, wash, dry, filter and concentrate", "stage": "workup", "temperature": None, "time": None}, {"operation": "purify on 100 g SNAP column with 15–90% diethyl ether/pentane gradient", "stage": "purification", "temperature": None, "time": None}], [Y("S13", 48, raw="0.375 g, 48%"), A("S13", 0.375, "g"), S(["S12"], "E/Z = 2.5:7.5", basis="source description of charged S12", scope="starting-material geometry composition"), S(["S13"], "trans/cis 8.5:2.5", basis="source-reported isolated product mixture", scope="isolated-product relative-isomer ratio")], ["C-C", "C-O"], name="acid-promoted formation of α-methylene-γ-lactone S13", transformation="S12 and 3-phenyl-2-propynal to hydroxy lactone S13", scale={"raw_text": "S12 2.9 mmol", "number": 2.9, "unit": "mmol"}, solvent={"raw": "toluene 7.8 mL", "composition": [mid("toluene")], "total_volume": Q("7.8 mL", 7.8, "mL")}, notes=("A pure trans-lactone fraction was collected only for characterization; the reported preparative product was a trans/cis 8.5:2.5 mixture.",), evidence_override=e13)

add_reaction("rxn_p22_S13_to_S14", [22], "page 22 compound S14 preparation", [mass_amt("S13", "0.11 g", 0.11, "g", "0.41 mmol", 0.41), mass_amt("DMP", "0.21 g", 0.21, "g", "0.49 mmol", 0.49), P("DCM", "solvent", volume=Q("5.8 mL", 5.8, "mL")), P("silica", "purification_agent", stage="purification", mass=Q("25 g SNAP column", 25, "g")), P("pentane", "purification_solvent", stage="purification"), P("S14", "product", mass=Q("96 mg", 96, "mg"))], [{"operation": "charge pure-trans S13 in DCM and add Dess–Martin periodinane", "temperature": ("room temperature", "r.t."), "time": Q("3 h", 3, "h")}, {"operation": "concentrate and purify on 25 g SNAP column with 50–100% diethyl ether/pentane", "stage": "purification", "temperature": None, "time": None}], [Y("S14", 88, raw="96 mg, 88%"), A("S14", 96, "mg")], ["DB-C-O"], name="Dess–Martin oxidation", transformation="oxidation of trans alcohol S13 to ketone S14", scale={"raw_text": "S13 0.41 mmol", "number": 0.41, "unit": "mmol"}, solvent={"raw": "CH2Cl2 5.8 mL", "composition": [mid("DCM")], "total_volume": Q("5.8 mL", 5.8, "mL")})

r15 = add_reaction("rxn_p22_S14_to_S15_two_step", [22], "page 22 two-step compound S15 preparation", [mass_amt("S14", "0.129 g", 0.129, "g", "0.48 mmol", 0.48), mass_amt("CeCl3", "1.13 g", 1.13, "g", "4.56 mmol", 4.56, role="additive"), P("THF", "solvent", volume=Q("8.7 mL plus 11.7 mL substrate solution", 20.4, "mL")), vol_amt("ethynylMgBr", "9.62 mL of 0.5 M in THF", 9.62, "mL", "4.80 mmol", 4.80), P("S15_alcohol", "product", mass=Q("0.178 g crude", 0.178, "g"), amount=Q("0.60 mmol crude", 0.60, "mmol"), product_role="intermediate"), P("DCM", "solvent", volume=Q("2.6 mL", 2.6, "mL")), mass_amt("DMAP", "7.42 mg", 7.42, "mg", "0.061 mmol", 0.061, role="catalyst"), vol_amt("Et3N", "0.85 mL", 0.85, "mL", "6.07 mmol", 6.07), vol_amt("Ac2O", "0.29 mL", 0.29, "mL", "3.03 mmol", 3.03, role="reactant"), P("NH4Cl", "workup_agent", stage="workup", volume=Q("5 mL at each workup", 10, "mL")), P("Et2O", "workup_agent", stage="workup"), P("Celite", "workup_agent", stage="workup"), P("brine", "workup_agent", stage="workup"), P("Na2SO4", "workup_agent", stage="workup"), P("pentane", "purification_solvent", stage="purification"), P("S15", "product", mass=Q("0.132 g", 0.132, "g"))], [{"operation": "charge anhydrous CeCl3 in a nitrogen glovebox; add THF and sonicate", "temperature": ("room temperature", "r.t."), "time": Q("4 h", 4, "h"), "atmosphere": "nitrogen"}, {"operation": "stir CeCl3 suspension", "temperature": ("room temperature", "r.t."), "time": Q("12 h", 12, "h"), "atmosphere": "nitrogen"}, {"operation": "cool, add ethynylmagnesium bromide and hold", "temperature": Q("0 °C", 0, "°C"), "time": Q("45 min", 45, "min"), "atmosphere": "nitrogen"}, {"operation": "cool reagent solution to −78 °C and cannulate into S14/THF at −78 °C; hold", "temperature": Q("−78 °C", -78, "°C"), "time": Q("30 min", 30, "min"), "atmosphere": "nitrogen"}, {"operation": "quench with saturated NH4Cl/ether, warm and stir, then Celite-filter/extract/wash/dry/concentrate", "stage": "workup", "temperature": ("−78 °C to room temperature", None), "time": Q("1 h at room temperature", 1, "h")}, {"operation": "take unstable crude propargyl alcohol directly into DCM; add DMAP, triethylamine and acetic anhydride", "temperature": ("room temperature", "r.t."), "time": Q("12 h", 12, "h"), "atmosphere": "nitrogen"}, {"operation": "dilute with NH4Cl, Celite-filter, extract, wash, dry and purify with 25% diethyl ether/pentane", "stage": "workup", "temperature": None, "time": None}], [A("S15_alcohol", 0.178, "g", raw="0.178 g crude", role="crude propargyl alcohol intermediate", basis="crude mass after first-step workup"), {"kind": "product_amount", "targets": ["S15_alcohol"], "raw": "0.60 mmol crude", "value": q(0.60, "mmol", raw_text="0.60 mmol"), "role": "crude propargyl alcohol intermediate", "basis": "source-reported crude amount"}, Y("S15", 82, raw="0.132 g, 82%", basis="isolated after two-step sequence and flash chromatography", ytype="two-step overall isolated yield", scope="S14 to final acetate S15"), A("S15", 0.132, "g")], ["C-C", "C-O"], name="two-step ethynyl addition and acetylation to S15", transformation="S14 to propargyl alcohol intermediate then acetate S15", scale={"raw_text": "S14 0.48 mmol", "number": 0.48, "unit": "mmol"}, solvent={"raw": "THF 20.4 mL across addition solutions; DCM 2.6 mL for acetylation", "composition": [mid("THF"), mid("DCM")]}, notes=("The crude alcohol was reported as 0.178 g, 0.60 mmol from 0.48 mmol S14; both printed values are preserved without correction.", "The 82% yield spans S14 to isolated acetate S15 and is not assigned to acetylation alone."), unresolved=("The crude-intermediate mass/amount exceeds the charged S14 amount and may include impurities; no correction was made.",))

add_reaction("rxn_p22_S15_to_S16", [22], "page 22 compound S16 preparation", [mass_amt("S15", "115 mg", 115, "mg", "0.33 mmol", 0.33), vol_amt("dimethylamine", "0.45 mL of 2 M in THF", 0.45, "mL", "0.83 mmol", 0.83), P("THF", "solvent", volume=Q("3.3 mL", 3.3, "mL")), P("S16_michael", "product", mass=Q("139 mg crude", 139, "mg"), amount=Q("0.37 mmol crude", 0.37, "mmol"), product_role="intermediate"), mass_amt("Stryker", "716 mg", 716, "mg", "0.37 mmol", 0.37, role="reagent"), P("toluene", "solvent", volume=Q("5.7 mL", 5.7, "mL")), P("S16_crude_allene", "product", mass=Q("153 mg crude", 153, "mg"), amount=Q("0.47 mmol crude", 0.47, "mmol"), product_role="intermediate"), P("NH4Cl", "workup_agent", stage="workup", volume=Q("5 mL saturated solution", 5, "mL")), P("Celite", "workup_agent", stage="workup"), P("Et2O", "workup_agent", stage="workup"), P("brine", "workup_agent", stage="workup", volume=Q("10 mL", 10, "mL")), P("Na2SO4", "workup_agent", stage="workup"), P("DCM", "solvent", volume=Q("18 mL", 18, "mL")), P("silica", "reagent", mass=Q("4 g", 4, "g"), amount=Q("66 mmol", 66, "mmol")), P("pentane", "purification_solvent", stage="purification"), P("S16", "product", mass=Q("36.1 mg", 36.1, "mg"))], [{"operation": "treat S15 in THF with dimethylamine", "temperature": ("room temperature", "r.t."), "time": Q("11 h", 11, "h"), "atmosphere": "nitrogen"}, {"operation": "concentrate to crude Michael adduct", "stage": "workup", "temperature": None, "time": None}, {"operation": "dilute crude adduct with degassed toluene, add Stryker reagent and stir", "temperature": ("room temperature", "r.t."), "time": Q("3 h", 3, "h"), "atmosphere": "nitrogen"}, {"operation": "quench with saturated NH4Cl and stir open to air", "stage": "workup", "temperature": ("room temperature", "r.t."), "time": Q("12 h", 12, "h"), "atmosphere": "air"}, {"operation": "Celite-filter, extract, wash, dry and concentrate to crude allene-yne", "stage": "workup", "temperature": None, "time": None}, {"operation": "dilute crude material with DCM, add silica and stir capped", "temperature": ("room temperature", "r.t."), "time": Q("11 h", 11, "h")}, {"operation": "filter, concentrate and flash chromatograph with 10% diethyl ether/pentane", "stage": "purification", "temperature": None, "time": None}], [A("S16_michael", 139, "mg", raw="139 mg crude", role="crude Michael-adduct intermediate", basis="crude mass after concentration"), {"kind": "product_amount", "targets": ["S16_michael"], "raw": "0.37 mmol crude", "value": q(0.37, "mmol", raw_text="0.37 mmol"), "role": "crude Michael-adduct intermediate", "basis": "source-reported crude amount"}, A("S16_crude_allene", 153, "mg", raw="153 mg crude", role="crude allene-yne intermediate", basis="crude mass after Stryker/workup stage"), {"kind": "product_amount", "targets": ["S16_crude_allene"], "raw": "0.47 mmol crude", "value": q(0.47, "mmol", raw_text="0.47 mmol"), "role": "crude allene-yne intermediate", "basis": "source-reported crude amount"}, Y("S16", 39, raw="36.1 mg, 39%"), A("S16", 36.1, "mg")], ["DB-C-C", "C-H"], name="dimethylamine-assisted reductive deoxyallenylation to S16", transformation="propargyl acetate S15 to allene-yne S16 through crude intermediates", scale={"raw_text": "S15 0.33 mmol", "number": 0.33, "unit": "mmol"}, solvent={"raw": "THF 3.3 mL; degassed toluene 5.7 mL; DCM 18 mL", "composition": [mid("THF"), mid("toluene"), mid("DCM")]}, notes=("Crude intermediate amounts exceed the starting amount and were not treated as yields or conversions.",), unresolved=("The crude Michael-adduct and crude allene-yne structures are not independently characterized in the source.",))

add_reaction("rxn_p22_S16_to_S17", [22], "page 22 compound S17 preparation", [mass_amt("S16", "15 mg", 15, "mg", "0.054 mmol", 0.054), mass_amt("Rhcat", "1 mg", 1, "mg", "5.6 × 10−3 mmol", 0.0056, role="catalyst", mol_percent=Q("0.10 equiv", 10, "%")), P("toluene", "solvent", volume=Q("1.6 mL", 1.6, "mL")), P("CO", "reagent"), P("Celite", "workup_agent", stage="workup"), P("Et2O", "workup_agent", stage="workup"), P("pentane", "purification_solvent", stage="purification"), P("S17", "product", mass=Q("8.8 mg", 8.8, "mg"))], [{"operation": "degass S16/toluene by bubbling nitrogen", "temperature": ("room temperature", "r.t."), "time": Q("approximately 5 min", 5, "min", "approximately"), "atmosphere": "nitrogen"}, {"operation": "place under vacuum for 3–5 s and refill with CO three times", "temperature": ("room temperature", "r.t."), "time": {"raw_text": "3–5 s", "minimum": 3, "maximum": 5, "unit": "s"}, "atmosphere": "carbon monoxide"}, {"operation": "add rhodium catalyst in one portion; repeat vacuum/CO cycle three times", "temperature": ("room temperature", "r.t."), "time": {"raw_text": "3–5 s", "minimum": 3, "maximum": 5, "unit": "s"}, "atmosphere": "carbon monoxide"}, {"operation": "place in preheated oil bath and stir", "temperature": Q("90 °C", 90, "°C"), "time": Q("25 min", 25, "min"), "atmosphere": "carbon monoxide"}, {"operation": "cool, pass through Celite with ether, concentrate and flash chromatograph with 75% ether/pentane", "stage": "workup", "temperature": ("room temperature", "r.t."), "time": None}], [Y("S17", 54, raw="8.8 mg, 54%"), A("S17", 8.8, "mg")], ["C-C", "DB-C-O"], name="Rh(I)-catalysed allenic Pauson–Khand reaction to S17", transformation="intramolecular APKR of lactone-tethered allene-yne S16 to fused lactone S17", category="substrate_scope", scale={"raw_text": "S16 0.054 mmol", "number": 0.054, "unit": "mmol"}, solvent={"raw": "toluene 1.6 mL", "composition": [mid("toluene")], "total_volume": Q("1.6 mL", 1.6, "mL")})


def update_kinetic(suffix, substrate, sm_mass, sm_mmol, cs_mass, cs_mmol, cdcl3_ml, *, special=False):
    r = next(x for x in bundle["reactions"] if x["reaction_id"].endswith(suffix))
    evidence = list(dict.fromkeys(r["evidence_ids"] + ev(23)))
    r["evidence_ids"] = evidence
    specs = [mass_amt(substrate, f"{sm_mass} mg", sm_mass, "mg", f"{sm_mmol} mmol", sm_mmol), mass_amt("cysteamine", f"{cs_mass} mg", cs_mass, "mg", f"{cs_mmol} mmol", cs_mmol, role="reagent", equivalents=Q("15 equiv", 15, "equiv")), P("CDCl3", "solvent", volume=Q(f"{cdcl3_ml} mL", cdcl3_ml, "mL"))]
    r["participants"] = [make_participant(i, p, evidence, mid(substrate)) for i, p in enumerate(specs, 1)]
    if special:
        steps = [{"operation": "dissolve lactam 54 in CDCl3; separately dissolve cysteamine in an unreported CDCl3 volume; combine, cap and shake vigorously", "temperature": ("room temperature", "r.t."), "time": Q("2 min", 2, "min")}, {"operation": "monitor disappearance of methylene doublets at δ 6.17 and 5.44 by 1H NMR", "temperature": ("room temperature", "r.t."), "time": {"raw_text": "less than 10 min", "number": 10, "unit": "min", "qualifier": "lt"}}]
    else:
        steps = [{"operation": "prepare separate 11.5 mM substrate and 15-equiv cysteamine solutions in equal CDCl3 volumes; combine to 5.7 mM substrate", "temperature": Q("22 °C", 22, "°C"), "time": Q("10 s vortex mix", 10, "s")}, {"operation": "transfer to 5 mm NMR tube; acquire 1H NMR immediately and at regular intervals while held at 22 °C", "temperature": Q("22 °C", 22, "°C"), "time": None}, {"operation": "derive fraction remaining from Ha/Hb/Hc/Hd integrations; fit ln(fraction remaining) versus time and calculate pseudo-first-order half-life", "stage": "characterization", "temperature": Q("22 °C", 22, "°C"), "time": None}]
    r["steps"] = [make_step(i, s, evidence) for i, s in enumerate(steps, 1)]
    r["procedure"]["general_procedure_reference"] = extracted("General Procedure K", "General Procedure K", ev(23))
    r["notes"] = [n for n in r["notes"] if n.get("raw_value") != "Page 23 supplies the specific masses, amounts and CDCl3 volume for this already-recorded kinetic experiment."]
    r["notes"].append(extracted("Page 23 supplies the specific masses, amounts and CDCl3 volume for this already-recorded kinetic experiment.", evidence=ev(23)))
    return r


update_kinetic("rxn_p6_cysteamine_45", "45", 8, 0.026, 30, 0.39, 1.0)
update_kinetic("rxn_p6_cysteamine_49", "49", 1.65, 0.0043, 5.0, 0.0649, 0.76)
update_kinetic("rxn_p6_cysteamine_51", "51", 0.98, 0.0022, 2.8, 0.0334, 0.4)
update_kinetic("rxn_p6_cysteamine_52", "52", 1.2, 0.0029, 3.4, 0.044, 0.5)
update_kinetic("rxn_p6_cysteamine_54", "54", 2, 0.004, 5, 0.07, 0.7, special=True)
update_kinetic("rxn_p6_cysteamine_60", "60", 1.1, 0.0041, 5.0, 0.065, 0.80)
update_kinetic("rxn_p6_cysteamine_S17", "S17", 1.0, 0.0033, 3.8, 0.049, 0.66)


def internal_standard_run(suffix, substrate, sm_mass, sm_mmol, cs_mass, cs_mmol, hmb_raw, hmb_vol, hmb_mmol, cdcl3, ratio):
    ps = [mass_amt(substrate, f"{sm_mass} mg", sm_mass, "mg", f"{sm_mmol} mmol", sm_mmol), mass_amt("cysteamine", f"{cs_mass} mg", cs_mass, "mg", f"{cs_mmol} mmol", cs_mmol, role="reagent", equivalents=Q("15 equiv", 15, "equiv")), P("hexamethylbenzene", "additive"), P("CDCl3", "solvent", volume=Q(f"{cdcl3} mL", cdcl3, "mL"))]
    if hmb_vol is not None:
        ps[2]["volume"] = Q(hmb_raw, hmb_vol, "μL")
    if hmb_mmol is not None:
        ps[2]["amount"] = Q(f"{hmb_mmol} mmol", hmb_mmol, "mmol")
    return add_reaction(suffix, [23], f"page 23 {substrate}/cysteamine kinetic run with internal standard", ps, [{"operation": "prepare and combine substrate/cysteamine solutions according to General Procedure K with one equivalent hexamethylbenzene", "temperature": Q("22 °C", 22, "°C"), "time": Q("10 s vortex mix", 10, "s")}, {"operation": f"acquire 1H NMR at regular intervals; initial hexamethylbenzene/substrate integration was {ratio}", "temperature": Q("22 °C", 22, "°C"), "time": None}, {"operation": "calculate fraction remaining using the source-referenced internal-standard equation", "stage": "characterization", "temperature": None, "time": None}], [{"kind": "other", "targets": [substrate], "raw": "half-life nearly identical to run without internal standard", "value": "nearly identical to no-internal-standard result", "role": substrate, "basis": "source qualitative comparison of fitted pseudo-first-order half-lives", "scope": "internal-standard validation run"}], ["C-S", "C-H"], name="cysteamine pseudo-first-order kinetic run with internal standard", transformation=f"NMR monitoring of thiol addition to {substrate}", category="mechanistic_study", scale={"raw_text": f"{sm_mmol} mmol {substrate}", "number": sm_mmol, "unit": "mmol"}, solvent={"raw": f"CDCl3 {cdcl3} mL", "composition": [mid("CDCl3")], "total_volume": Q(f"{cdcl3} mL", cdcl3, "mL")}, notes=(f"Initial hexamethylbenzene/{substrate} integration ratio: {ratio}.",), unresolved=(("The source prints '0.7 μL mmol' for hexamethylbenzene without a separable mmol value; no amount was inferred." ,) if hmb_mmol is None else ()), general="General Procedure K with internal standard")


internal_standard_run("rxn_p23_cysteamine_49_internal_standard", "49", 1.7, 0.0043, 5, 0.065, "0.7 μL mmol", 0.7, None, 0.76, "2:1")
internal_standard_run("rxn_p23_cysteamine_52_internal_standard", "52", 1.2, 0.0029, 3.4, 0.044, "0.4 μL", 0.4, 0.0029, 0.05, "1:1")

# Enrich the existing S11 structural-confirmation experiment rather than duplicate it.
rs11 = next(r for r in bundle["reactions"] if r["reaction_id"].endswith("rxn_p7_48_tBuSH_to_S11"))
e11 = list(dict.fromkeys(rs11["evidence_ids"] + ev(23)))
rs11["evidence_ids"] = e11
s11_specs = [mass_amt("48", "10 mg", 10, "mg", "0.03 mmol", 0.03), vol_amt("tBuSH", "0.02 mL", 0.02, "mL", "0.2 mmol", 0.2, role="reactant"), P("acetone", "solvent", volume=Q("0.5 mL", 0.5, "mL")), P("PBS", "solvent", volume=Q("0.5 mL", 0.5, "mL")), mass_amt("KOH", "2 mg", 2, "mg", "0.04 mmol", 0.04, role="reagent"), P("water", "workup_agent", stage="workup", volume=Q("3 mL", 3, "mL")), P("DCM", "workup_agent", stage="workup", volume=Q("10 mL plus 2 × 5 mL extraction", 20, "mL")), P("brine", "workup_agent", stage="workup", volume=Q("5 mL", 5, "mL")), P("MgSO4", "workup_agent", stage="workup"), P("S11", "product", mass=Q("3 mg", 3, "mg"))]
rs11["participants"] = [make_participant(i, p, e11, mid("48")) for i, p in enumerate(s11_specs, 1)]
rs11["steps"] = [make_step(i, s, e11) for i, s in enumerate([{ "operation": "stir 48 and tert-butyl thiol in acetone/PBS at pH 8; TLC shows only starting material", "temperature": ("room temperature", "r.t."), "time": Q("24 h", 24, "h")}, {"operation": "add KOH to raise pH to 12 and continue until starting material is consumed by TLC", "temperature": ("room temperature", "r.t."), "time": Q("within 30 min", 30, "min", "within")}, {"operation": "dilute with water/DCM, extract, brine wash, MgSO4 dry, filter, concentrate and flash chromatograph (25–50% EtOAc/hexanes)", "stage": "workup", "temperature": None, "time": None}], 1)]
rs11["outcomes"][0]["measurement_basis"] = extracted("isolated after flash column chromatography", "isolated after flash column chromatography", ev(23))
rs11["outcomes"][1]["measurement_basis"] = extracted("1H NMR integration of Ha at δ 4.66 (minor) and 4.50 (major)", "1H NMR integration of Ha at δ 4.66 (minor) and 4.50 (major)", ev(23))
rs11["outcomes"] = [o for o in rs11["outcomes"] if o["outcome_id"] != "o3"]
rs11["outcomes"].append(make_outcome(3, A("S11", 3, "mg"), ev(23)))
rs11["notes"] = [n for n in rs11["notes"] if n.get("raw_value") != "The pH 8 stage was a true no-reaction observation after 24 h; reaction occurred only after adjustment to pH 12."]
rs11["notes"].append(extracted("The pH 8 stage was a true no-reaction observation after 24 h; reaction occurred only after adjustment to pH 12.", evidence=ev(23)))

# Page 24 adds detailed assay methods to the existing mechanistic records.
for suffix in ("rxn_p9_alarm_nmr_51", "rxn_p9_alarm_nmr_S17"):
    r = next(x for x in bundle["reactions"] if x["reaction_id"].endswith(suffix))
    r["evidence_ids"] = list(dict.fromkeys(r["evidence_ids"] + ev(24)))
    r["steps"] = [s for s in r["steps"] if not (s["operation"].get("raw_value") or "").startswith(("load 160 μL samples", "record data at 25 °C"))]
    r["steps"].extend([make_step(len(r["steps"]) + 1, {"operation": "load 160 μL samples into 3 mm SampleJet tubes; store at ambient temperature while queued", "stage": "characterization", "temperature": ("ambient temperature", "ambient"), "time": None}, ev(24)), make_step(len(r["steps"]) + 2, {"operation": "record data at 25 °C on a 700 MHz NMR cryoprobe every hour", "stage": "characterization", "temperature": Q("25 °C", 25, "°C"), "time": Q("10 h total", 10, "h")}, ev(24))])
    r["notes"] = [n for n in r["notes"] if not (n.get("raw_value") or "").startswith("Page 24 method:")]
    r["notes"].append(extracted("Page 24 method: test compound 400 μM, 13C-methyl-labelled La antigen 50 μM, ±20 mM DTT, 37 °C for 90 min; 160 μL sample volume.", evidence=ev(24)))

for suffix in ("rxn_p9_p10_alarm_msps_S17", "rxn_p9_p10_alarm_msps_CPM", "rxn_p9_p10_alarm_msps_51", "rxn_p9_p10_alarm_msps_fluconazole", "rxn_p9_p10_alarm_msps_DMSO"):
    r = next(x for x in bundle["reactions"] if x["reaction_id"].endswith(suffix))
    r["evidence_ids"] = list(dict.fromkeys(r["evidence_ids"] + ev(24)))
    r["notes"] = [n for n in r["notes"] if not (n.get("raw_value") or "").startswith(("Page 24 sample preparation", "Figure 6 caption reports compounds"))]
    r["notes"].append(extracted("Page 24 sample preparation used 2.5 μL La antigen (2 μg/μL), 15.5 μL of 10 mM test-compound stock, dilution to 310 μL with 20 mM Tris pH 8.0, and 1 h at 37 °C before 50 mM IAD treatment.", evidence=ev(24)))
    r["notes"].append(extracted("Figure 6 caption reports compounds at 50 mM final concentration, whereas the page-24 charge/dilution corresponds arithmetically to 0.5 mM; both source statements are retained without choosing one.", evidence=ev(9, 24)))
    r["review"]["unresolved"] = [
        item for item in r["review"]["unresolved"]
        if not item.startswith("Figure 6 concentration conflict:")
    ]
    r["review"]["unresolved"].append("Figure 6 concentration conflict: caption says 50 mM final, but page-24 15.5 μL of 10 mM diluted to 310 μL corresponds to 0.5 mM; authorial clarification is required.")


def context(context_id, page, region, classification, summary, *, pages=None):
    return {"context_id": context_id, "document_id": DOC, "location": {"page": page, "region": region}, "classification": classification, "summary": summary, "evidence_ids": ev(*(pages or [page])), "review_status": "needs_review"}


contexts = [
    context("ctx_p19_characterization", 19, "44–48 spectral blocks and compound-45 crystallography statement", "characterization", "NMR/IR/HRMS/TLC and X-ray confirmation were retained as characterization, not separate syntheses."),
    context("ctx_p20_characterization", 20, "49–54 spectral blocks", "characterization", "Spectral peak tables and analytical values were not promoted to reactions."),
    context("ctx_p21_characterization", 21, "55, 58, 59, 60, S9 and S10 spectral blocks", "characterization", "Spectral data remain source context; the reported dr measurement for 58 is linked to the specific reaction outcome."),
    context("ctx_p22_characterization", 22, "S13–S17 spectral and analytical blocks", "characterization", "NMR, IR, HRMS and TLC blocks are characterization, not additional preparations."),
    context("ctx_p23_general_procedure_K", 23, "General Procedure K and kinetic calculation text", "analytical_method", "Pseudo-first-order NMR preparation, integration and fitting method was attached to the specific kinetic records; equations/peak definitions were not treated as reactions."),
    context("ctx_p23_S11_characterization", 23, "S11 spectral block", "characterization", "S11 peak lists and relative-abundance-style spectral data were not made into synthesis records."),
    context("ctx_p24_HPLC", 24, "General HPLC protocol", "characterization_method", "Analytical purity/HPLC gradient, wavelength and column details are characterization rather than synthesis."),
    context("ctx_p24_biological_assays", 24, "cell culture, NF-κB reporter and Vero assay methods", "biological_method", "Cell culture and reporter/viability assay methods are biological context, not chemical synthesis experiments."),
    context("ctx_p24_alarm_methods", 24, "ALARM NMR and La-antigen MS sample-preparation methods", "mechanistic_method", "Detailed ALARM conditions were linked to existing mechanistic records. The MS workflow continues onto page 25, which is outside this visual batch."),
]
bundle["source_context"].extend(contexts)

# Coverage record links: include new records plus same experiments enriched from prior pages.
new_batch = [r for r in bundle["reactions"] if r["reaction_id"].startswith(batch_prefixes)]
linked_existing = {
    19: [r44["reaction_id"]],
    23: [next(r for r in bundle["reactions"] if r["reaction_id"].endswith(s))["reaction_id"] for s in ("rxn_p6_cysteamine_45", "rxn_p6_cysteamine_49", "rxn_p6_cysteamine_51", "rxn_p6_cysteamine_52", "rxn_p6_cysteamine_54", "rxn_p6_cysteamine_60", "rxn_p6_cysteamine_S17", "rxn_p7_48_tBuSH_to_S11")],
    24: [next(r for r in bundle["reactions"] if r["reaction_id"].endswith(s))["reaction_id"] for s in ("rxn_p9_alarm_nmr_51", "rxn_p9_alarm_nmr_S17", "rxn_p9_p10_alarm_msps_S17", "rxn_p9_p10_alarm_msps_CPM", "rxn_p9_p10_alarm_msps_51", "rxn_p9_p10_alarm_msps_fluconazole", "rxn_p9_p10_alarm_msps_DMSO")],
}
page_ids = {p: [] for p in range(19, 25)}
for r in new_batch:
    for p in range(19, 25):
        if any(e in r["evidence_ids"] for e in ev(p)):
            page_ids[p].append(r["reaction_id"])
for p, ids in linked_existing.items():
    page_ids[p].extend(ids)
for p in page_ids:
    page_ids[p] = list(dict.fromkeys(page_ids[p]))

coverage = {"document_id": DOC, "semantic_accuracy_verified": False, "pages": [
    {"page": 19, "state": "attempted", "regions": [
        {"region": "top-left continuation of compound 44 Run 1 and Run 2", "finding": "cross-page Run 1 completed with page-19 image evidence; separate Run 2 recorded with 17 mg/54% outcome", "evidence_ids": ev(19)},
        {"region": "General Procedures H/I and compound 45 Runs 1–3", "finding": "slow-addition and single-portion APKR operations transcribed separately, including CO cycles, 110 °C stages and workup", "evidence_ids": ev(19)},
        {"region": "compounds 46–48 preparations", "finding": "five distinct current APKR runs recorded with their own quantities and yields", "evidence_ids": ev(19)},
        {"region": "spectral blocks and 45 crystallography", "finding": "characterization only; no extra synthesis invented", "evidence_ids": ev(19)},
    ], "record_ids": page_ids[19], "unresolved": ["Compound 46 Run 1 reports consumption after 3 h but does not unambiguously state whether this is total or post-addition time.", "Complex fused products remain without stereo-complete isomeric SMILES."]},
    {"page": 20, "state": "attempted", "regions": [
        {"region": "compounds 49–54 APKR preparations", "finding": "eleven individual runs/conditions recorded rather than merged by product", "evidence_ids": ev(20)},
        {"region": "compound 55 heading at bottom", "finding": "heading continues to the page-21 preparation and is linked there", "evidence_ids": ev(20)},
        {"region": "spectral and analytical blocks", "finding": "characterization retained as context, not reactions", "evidence_ids": ev(20)},
    ], "record_ids": page_ids[20], "unresolved": ["Compound 40 Run 1 prints 6 mg and 0.2 mmol; both are preserved without correction.", "Catalyst mmol is not printed for two catalyst-solution runs."]},
    {"page": 21, "state": "attempted", "regions": [
        {"region": "compound 55 Runs 1–2 and compounds 58/59", "finding": "four APKR experiments recorded, including 56 starting dr, 58 isolated dr/method and unassigned new-center stereochemistry", "evidence_ids": ev(21)},
        {"region": "General Procedure J and compounds 60, S9 and S10", "finding": "three N-acetylations recorded with separate quantities, times and isolated yields", "evidence_ids": ev(21)},
        {"region": "start of S13 preparation", "finding": "cross-page reaction merged with page 22, preserving S12 E/Z and ordered low-temperature stages", "evidence_ids": ev(21)},
        {"region": "spectral blocks", "finding": "characterization only", "evidence_ids": ev(21)},
    ], "record_ids": page_ids[21], "unresolved": ["S10 acetic-anhydride mmol and numeric product dr are unreported.", "S12 complete structure is unavailable on these pages."]},
    {"page": 22, "state": "attempted", "regions": [
        {"region": "completion of S13 and preparations of S14–S17", "finding": "five distinct preparations/overall sequences recorded with ordered temperature, atmosphere, workup and purification stages", "evidence_ids": ev(22)},
        {"region": "S15 two-step sequence", "finding": "crude propargyl alcohol and final acetate kept distinct; 82% retained as S14-to-S15 overall yield", "evidence_ids": ev(22)},
        {"region": "S16 crude-intermediate sequence", "finding": "crude Michael-adduct/allene amounts retained as amounts rather than yields", "evidence_ids": ev(22)},
        {"region": "spectral blocks", "finding": "NMR/IR/HRMS/TLC characterization not promoted to synthesis", "evidence_ids": ev(22)},
    ], "record_ids": page_ids[22], "unresolved": ["Crude S15/S16 intermediate amounts exceed charged precursor amounts and likely include impurities; source values are preserved.", "Temporary crude intermediates lack independently verified complete structures."]},
    {"page": 23, "state": "attempted", "regions": [
        {"region": "General Procedure K and cysteamine runs", "finding": "seven existing kinetic experiments enriched with exact masses/amounts/volumes; two internal-standard runs added separately", "evidence_ids": ev(23)},
        {"region": "compound 54 special kinetic run", "finding": "separate dissolution, 2 min shaking and <10 min disappearance retained; cysteamine-solution CDCl3 volume remains unreported", "evidence_ids": ev(23)},
        {"region": "48 plus tert-butyl thiol to S11", "finding": "existing mechanistic experiment enriched with pH-8 no-reaction stage, pH-12 trigger, workup, 23% yield and NMR-based 1.6:1 dr", "evidence_ids": ev(23)},
        {"region": "S11 spectral table", "finding": "characterization only; no independent mass-spectrum experiment invented", "evidence_ids": ev(23)},
    ], "record_ids": page_ids[23], "unresolved": ["The printed hexamethylbenzene entry '0.7 μL mmol' lacks a separable amount value.", "Exact configurations of the S11 diastereomers remain unresolved."]},
    {"page": 24, "state": "attempted", "regions": [
        {"region": "HPLC protocol", "finding": "analytical characterization method only", "evidence_ids": ev(24)},
        {"region": "cell culture, NF-κB reporter and Vero methods", "finding": "biological assay context rather than chemical synthesis", "evidence_ids": ev(24)},
        {"region": "ALARM NMR method", "finding": "detailed 51/S17 assay conditions linked to existing mechanistic records", "evidence_ids": ev(24)},
        {"region": "La antigen expression and MS sample preparation", "finding": "method details linked to existing ALARM MSPS rows; analysis continues on page 25 outside this batch", "evidence_ids": ev(24)},
    ], "record_ids": page_ids[24], "unresolved": ["Figure 6 says 50 mM final compound, while page-24 stock/dilution arithmetic gives 0.5 mM; both source statements are retained as an unresolved conflict.", "The ALARM MS analysis continues onto page 25, which was not visually inspected or claimed covered in this batch."]},
]}

inventory = {"batch_id": "allene-stereo-v31-200-20261017", "document_id": DOC, "visually_inspected_original_pdf_pages": [19, 20, 21, 22, 23, 24], "new_reaction_ids": [r["reaction_id"] for r in new_batch], "updated_existing_reaction_ids": list(dict.fromkeys(linked_existing[19] + linked_existing[23] + linked_existing[24])), "new_molecule_ids": [mid(k) for k in ("S9", "S10", "S12", "S13", "S14", "S15_alcohol", "S15", "S16_michael", "S16_crude_allene", "S16", "polyPPh3", "pentane", "TfOH", "Na2SO4", "DMP", "Celite", "CeCl3", "dimethylamine", "PBS", "KOH", "hexamethylbenzene")], "source_context_ids": [c["context_id"] for c in contexts], "limitations": ["All new records remain needs_review and partial.", "No SI association was inferred.", "Complex relative/mixture stereochemistry was retained without speculative isomeric SMILES.", "Page 25 continuation was not inspected or claimed covered."]}

(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
(RUN / "coverage" / "pages_0019_0024.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model" / "pages_0019_0024_inventory.json").write_text(json.dumps(inventory, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]), "source_context": len(bundle["source_context"]), "new_records": len(new_batch), "updated_records": len(inventory["updated_existing_reaction_ids"])}, ensure_ascii=False))
