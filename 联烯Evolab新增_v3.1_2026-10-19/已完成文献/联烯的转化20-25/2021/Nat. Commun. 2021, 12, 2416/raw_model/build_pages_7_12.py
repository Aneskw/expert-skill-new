#!/usr/bin/env python3
"""Add source-derived records and evidence for original PDF pages 7--12."""

from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st


DOC = "doc_fb3fead4a044bff7"
SUBMITTED = RUN / "submitted.json"


def ev(page: int, kind: str = "image") -> str:
    return f"{DOC}:page{page}:{kind}"


def q(number, unit, **extra):
    value = {"number": number, "unit": unit}
    value.update(extra)
    return value


def f(raw=None, normalized=None, evidence=(), method="original_extraction", status=None,
      derivation=None):
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   derivation=derivation)


def missing(status="unreported", evidence=()):
    method = "not_applicable" if status == "not_applicable" else "original_extraction"
    return st.fact(status=status, method=method)


bundle = json.loads(SUBMITTED.read_text())
# The earlier page-1--6 bundle predates rendering pages 7--12.  Replace only
# the immutable document/evidence registry with the fresh `standard init`
# output so the submitted evidence set exactly matches the prepared dataset.
fresh_template = json.loads((RUN / "raw_model/blank_template_pages_0007_0012.json").read_text())
bundle["documents"] = fresh_template["documents"]
molecules = {m["molecule_id"]: m for m in bundle["molecules"]}
reactions = {r["reaction_id"]: r for r in bundle["reactions"]}


def add_molecule(mid, label, page, *, name=None, chirality="unresolved",
                 geometry="unresolved", representation="not_encoded", note=None):
    if mid in molecules:
        return mid
    evidence = [ev(page, "image"), ev(page, "text")]
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    m["names"]["original"] = f(name or label, name or label, evidence)
    m["names"]["english"] = f(name or label, name or label, evidence)
    m["cas_numbers"] = missing(evidence=evidence)
    m["representation_type"] = "unresolved"
    m["smiles"] = missing(evidence=evidence)
    m["stereochemistry"]["chirality_status"] = f(
        chirality, chirality, evidence, method="inference",
        derivation="No complete structure was encoded; stereochemical status remains unresolved unless explicitly labelled by the source.")
    m["stereochemistry"]["geometry_status"] = f(
        geometry, geometry, evidence, method="inference",
        derivation="No complete structure was encoded; geometric stereochemistry remains unresolved unless explicitly labelled by the source.")
    m["stereochemistry"]["representation_status"] = f(
        representation, representation, evidence, method="inference",
        derivation="The source label is retained without asserting a complete stereo-capable graph.")
    if note:
        m["identity_note"] = f(note, note, evidence)
    molecules[mid] = m
    bundle["molecules"].append(m)
    return mid


def add_occurrence(mid, label, page):
    m = molecules[mid]
    evidence = [ev(page, "image"), ev(page, "text")]
    if not any(o["compound_label"]["raw_value"] == label and ev(page, "image") in o["evidence_ids"]
               for o in m["occurrences"]):
        m["occurrences"].append({
            "document_id": DOC,
            "compound_label": f(label, label, evidence),
            "evidence_ids": evidence,
        })


def add_reaction_evidence(r, page, region):
    evidence = [ev(page, "image"), ev(page, "text")]
    for eid in evidence:
        if eid not in r["evidence_ids"]:
            r["evidence_ids"].append(eid)
    alt = f(region, {"page": page, "region": region}, evidence)
    if not any(a.get("raw_value") == region for a in r["location"]["alternatives"]):
        r["location"]["alternatives"].append(alt)


def participant(pid, mid, role, evidence, stage="reaction", product_role=None):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    if product_role is not None:
        p["product_role"] = f(product_role, product_role, evidence)
    return p


def set_charge(p, key, raw, normalized, evidence, reference=None,
               method="original_extraction", derivation=None):
    p["charges"][key] = f(raw, normalized, evidence, method=method, derivation=derivation)
    if key in {"equivalents", "mol_percent"} and reference:
        p["equivalence_reference"] = f(
            reference, reference, evidence, method="inference",
            derivation="Loading is interpreted relative to the indicated substrate used as the reaction-scale reference.")


def add_participant(r, p):
    if not any(x["participant_id"] == p["participant_id"] for x in r["participants"]):
        r["participants"].append(p)


def outcome(oid, kind, target, raw, normalized, evidence, *, product_role,
            yield_type=None, basis_raw=None, basis_norm=None, scope=None):
    o = st.new_outcome(oid, kind, [target] if target else [], evidence)
    o["product_role"] = f(product_role, product_role, evidence)
    o["value"] = f(raw, normalized, evidence)
    if yield_type is not None:
        o["yield_type"] = f(yield_type, yield_type, evidence)
    if basis_raw is not None:
        o["measurement_basis"] = f(basis_raw, basis_norm or basis_raw, evidence)
    if scope is not None:
        o["scope"] = f(scope, scope, evidence)
    return o


TAXONOMY = {"C-C": 2, "C-O": 3, "DB-C-C": 13, "TB-C-C": 21}


def set_classification(r, specs, evidence, limitation):
    changes = []
    labels = []
    for i, (code, before, after, change_type) in enumerate(specs, 1):
        change_id = f"b{i}"
        changes.append(f(
            "source-drawing bond change",
            {"change_id": change_id,
             "atoms": [{"atom_map": None, "element": "C"},
                       {"atom_map": None, "element": "O" if code == "C-O" else "C"}],
             "before_order": before, "after_order": after, "change_type": change_type},
            evidence, method="structure_derivation",
            derivation="Conservative change type from the inspected reaction drawing; complete atom mapping is not asserted."))
        row = next((x for x in labels if x["code"] == code), None)
        if row:
            row["supporting_change_ids"].append(change_id)
        else:
            labels.append({"category_number": TAXONOMY[code], "code": code,
                           "supporting_change_ids": [change_id]})
    r["classification"]["bond_changes"] = changes
    r["classification"]["labels"] = f(
        "source-drawing bond types", labels, evidence, method="structure_derivation",
        derivation="Only defensible formed connections or bond-order increases are indexed; other changes remain unresolved.")
    r["classification"]["coverage"] = f(
        "partial", "partial", evidence, method="inference",
        derivation="The figure supports the listed transformation class, but no reviewed full atom mapping is available.")
    r["classification"]["limitations"] = [limitation]


def new_reaction(rid, page, record_type, categories, name, transformation,
                 parts, *, outcomes=None, temperature=None, time=None,
                 atmosphere=None, pressure=None, solvent=None,
                 concentration=None, scale=None, classification=(),
                 location=None, unresolved=(), notes=()):
    evidence = [ev(page, "image"), ev(page, "text")]
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, evidence)
    region = location or f"PDF page {page}"
    r["location"] = f(region, {"page": page, "region": region}, evidence)
    r["reaction_name"] = f(name, name, evidence)
    r["transformation_type"] = f(
        transformation, transformation, evidence, method="inference",
        derivation="Normalized from the source caption and reaction drawing.")
    r["participants"] = parts
    if scale:
        r["reaction_scale"] = f(scale[0], {"reference_molecule_id": scale[1],
                                             "amount": scale[2]}, evidence)
    if solvent:
        r["solvent_conditions"]["composition"] = f(solvent, [solvent], evidence)
    if concentration:
        r["solvent_conditions"]["concentration"] = f(
            concentration[0], concentration[1], evidence)
        r["solvent_conditions"]["concentration_basis"] = f(
            concentration[2], concentration[2], evidence, method="inference",
            derivation="The source scale identifies the allene substrate as the concentration reference.")
    step = st.new_step("reaction", evidence)
    step["participant_ids"] = [p["participant_id"] for p in parts if p["stage"] == "reaction"]
    step["operation"] = f(name, name, evidence)
    if temperature:
        step["temperature"] = f(temperature[0], temperature[1], evidence)
    if time:
        step["time"] = f(time[0], time[1], evidence)
    if atmosphere:
        step["atmosphere"] = f(atmosphere, atmosphere, evidence)
    if pressure:
        step["pressure"] = f(pressure, pressure, evidence)
    r["steps"] = [step]
    r["procedure"]["original_text"] = f(name, name, evidence)
    r["outcomes"] = outcomes or []
    if classification:
        set_classification(r, classification, evidence,
                           "Partial source-drawing classification; complete atom mapping, broken bonds, and implicit-H changes remain unresolved.")
    r["notes"] = [f(x, x, evidence) for x in notes]
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": list(unresolved) or
                   ["Machine extraction requires expert chemistry review."]}
    r["completeness"] = "partial"
    bundle["reactions"].append(r)
    reactions[r["reaction_id"]] = r
    return r


# New identities used on pages 7--11. Complete structures are not asserted.
new_ids = {
    "rh": add_molecule("m_rh_co_2_cl_dimer", "[Rh(CO)2Cl]2", 7),
    "agsbf6": add_molecule("m_agsbf6", "AgSbF6", 7),
    "co": add_molecule("m_carbon_monoxide", "CO", 7, name="carbon monoxide"),
    "bipy": add_molecule("m_2_2_bipyridine", "2,2'-bipy", 7, name="2,2'-bipyridine"),
    "cuoac": add_molecule("m_cu_oac_2_h2o", "Cu(OAc)2·H2O", 7),
    "et3n": add_molecule("m_et3n", "Et3N", 8, name="triethylamine"),
    "nah": add_molecule("m_nah", "NaH", 8, name="sodium hydride"),
    "etoac": add_molecule("m_ethyl_acetate", "ethyl acetate", 11),
    "silica": add_molecule("m_silica_gel", "silica gel", 11,
                           note="Chromatography/filter medium; no unique molecular structure is asserted."),
    "petether": add_molecule("m_petroleum_ether_60_90", "petroleum ether (60–90 °C)", 11,
                             note="Named solvent fraction with unspecified composition."),
    "iproh": add_molecule("m_isopropanol", "i-PrOH", 11, name="2-propanol"),
}


# Fig. 4b, reaction i: distinct 4.0 mmol scale-up experiment.
e7 = [ev(7, "image"), ev(7, "text")]
parts = []
p = participant("allene_1a", "m_1a", "reactant", e7)
set_charge(p, "amount", "1a (4.0 mmol)", q(4.0, "mmol"), e7)
parts.append(p)
p = participant("pyrazolone_2m", "m_2m", "reactant", e7)
set_charge(p, "equivalents", "2m (1.2 equiv)", q(1.2, "equiv"), e7, "m_1a")
parts.append(p)
p = participant("pd_catalyst", "m_pd2(dba)3_chcl3", "catalyst", e7)
set_charge(p, "mol_percent", "Pd2(dba)3·CHCl3 (2.5 mol%)", q(2.5, "mol%"), e7, "m_1a")
parts.append(p)
p = participant("ligand", "m_DACH_ZYC_Phos_C1", "ligand", e7)
set_charge(p, "mol_percent", "DACH-ZYC-Phos-C1 (5.0 mol%)", q(5.0, "mol%"), e7, "m_1a")
parts.append(p)
parts.append(participant("toluene", "m_toluene", "solvent", e7))
parts.append(participant("product_3am", "m_3am_S", "product", e7,
                         product_role="main C-allenylation product"))
scale_outcomes = [
    outcome("product_amount", "product_amount", "m_3am_S", "1.20 g", q(1.20, "g"), e7,
            product_role="main C-allenylation product",
            basis_raw="isolated product mass", basis_norm="isolated product mass"),
    outcome("isolated_yield", "yield", "m_3am_S", "93% y.", q(93, "%"), e7,
            product_role="main C-allenylation product", yield_type="isolated yield",
            basis_raw="Fig. 4 caption: y. isolated yield", basis_norm="isolated yield"),
    outcome("ee", "selectivity", "m_3am_S", "93% ee", q(93, "%"), e7,
            product_role="main C-allenylation product",
            scope={"metric": "ee", "major_isomer": "(S)-3am"}),
]
scale_rx = new_reaction(
    "fig4_scaleup_3am", 7, "specific_experiment", ["scale_up"],
    "gram-scale asymmetric allenylation to (S)-3am", "C-allenylation", parts,
    outcomes=scale_outcomes, temperature=("60 °C", q(60, "°C")),
    time=("13 h", q(13, "h")), solvent="toluene",
    scale=("1a (4.0 mmol)", "m_1a", q(4.0, "mmol")),
    classification=(("C-C", 0, 1, "formed"),), location="Fig. 4b, reaction i",
    unresolved=("The Fig. 4 caption does not state the ee measurement method or toluene volume.",
                "Full atom mapping, leaving-group changes, and implicit-H changes remain unresolved."))
scale_rx["general_scheme_ids"] = [f"{DOC}:fig1c_generic_allenylation"]


# Fig. 4b, reactions ii--iv: enrich the page-5 prose records with exact conditions.
deriv_specs = {
    "prose_3am_to_5am": {
        "region": "Fig. 4b, reaction ii",
        "substrate": "m_3am_S", "temperature": ("50 °C", q(50, "°C")),
        "time": ("13 h", q(13, "h")), "solvent": "toluene",
        "add": [("rh_catalyst", new_ids["rh"], "catalyst", "mol_percent", "5.0 mol%", q(5.0, "mol%")),
                ("agsbf6", new_ids["agsbf6"], "additive", "mol_percent", "12.5 mol%", q(12.5, "mol%")),
                ("carbon_monoxide", new_ids["co"], "reagent", None, None, None),
                ("toluene", "m_toluene", "solvent", None, None, None)],
        "classification": (("C-C", 0, 1, "formed"), ("C-C", 0, 1, "formed")),
        "yield": "yield", "yield_type": "isolated yield",
        "extra_unresolved": "The source does not report the ee measurement method, substrate scale, or solvent volume in Fig. 4.",
    },
    "prose_3da_to_6da": {
        "region": "Fig. 4b, reaction iii",
        "substrate": "m_3da", "temperature": ("110 °C", q(110, "°C")),
        "time": ("18 h", q(18, "h")), "solvent": "toluene",
        "add": [("ophenol", "m_o_(propen_2_yl)phenol", "reactant", "equivalents", "3.0 equiv", q(3.0, "equiv")),
                ("pd_oac_2", "m_pd(oac)2", "catalyst", "mol_percent", "7.5 mol%", q(7.5, "mol%")),
                ("bipy", new_ids["bipy"], "ligand", "mol_percent", "7.5 mol%", q(7.5, "mol%")),
                ("cu_oac_2_h2o", new_ids["cuoac"], "additive", "equivalents", "0.5 equiv", q(0.5, "equiv")),
                ("toluene", "m_toluene", "solvent", None, None, None)],
        "classification": (("C-C", 0, 1, "formed"), ("C-O", 0, 1, "formed")),
        "yield": "yield", "yield_type": "isolated yield",
        "extra_unresolved": "The source does not report the ee/Z:E measurement methods, substrate scale, or solvent volume in Fig. 4.",
    },
    "prose_3da_to_7da": {
        "region": "Fig. 4b, reaction iv",
        "substrate": "m_3da", "temperature": None, "time": None,
        "solvent": "dichloromethane",
        "add": [("ticl4", "m_ticl4", "reagent", "equivalents", "1.5 equiv", q(1.5, "equiv")),
                ("dcm", "m_dichloromethane", "solvent", None, None, None)],
        "classification": (("TB-C-C", 2, 3, "order_increased"),),
        "yield": None, "yield_type": "isolated yield",
        "extra_unresolved": "The source does not report the ee measurement method, substrate scale, or solvent volume in Fig. 4; hydrogen migration is not atom-mapped.",
    },
}

# Existing molecule aliases differ from the printed formulas; add missing named identities when needed.
add_molecule("m_pd(oac)2", "Pd(OAc)2", 7)
add_molecule("m_dichloromethane", "dichloromethane", 7)

for suffix, spec in deriv_specs.items():
    rid = f"{DOC}:{suffix}"
    r = reactions[rid]
    add_reaction_evidence(r, 7, spec["region"])
    for pid, mid, role, key, raw, normalized in spec["add"]:
        existing = next(
            (x for x in r["participants"]
             if x["participant_id"] == pid
             or (x["molecule_id"] == mid and x["role"] == role)),
            None,
        )
        if existing is not None:
            p = existing
        else:
            p = participant(pid, mid, role, e7)
            add_participant(r, p)
        if key:
            set_charge(p, key, raw, normalized, e7, spec["substrate"])
    r["solvent_conditions"]["composition"] = f(spec["solvent"], [spec["solvent"]], e7)
    if suffix != "prose_3da_to_7da":
        step = st.new_step("fig4_reaction", e7)
        step["participant_ids"] = [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"]
        step["operation"] = f(spec["region"], spec["region"], e7)
        step["temperature"] = f(spec["temperature"][0], spec["temperature"][1], e7)
        step["time"] = f(spec["time"][0], spec["time"][1], e7)
        if suffix == "prose_3am_to_5am":
            step["pressure"] = f("CO balloon", "CO balloon", e7)
        r["steps"] = [step]
    else:
        r["steps"] = []
        for sid, temp_raw, temp_norm, time_raw, time_norm in [
            ("low_temperature_hold", "−78 °C", q(-78, "°C"), "2.5 h", q(2.5, "h")),
            ("room_temperature_hold", "rt", "rt", "5.5 h", q(5.5, "h")),
        ]:
            step = st.new_step(sid, e7)
            step["participant_ids"] = [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"]
            step["operation"] = f(f"{temp_raw} for {time_raw}", f"{temp_raw} for {time_raw}", e7)
            step["temperature"] = f(temp_raw, temp_norm, e7)
            step["time"] = f(time_raw, time_norm, e7)
            r["steps"].append(step)
    r["procedure"]["original_text"] = f(spec["region"], spec["region"], e7)
    if spec["yield"]:
        y = next(o for o in r["outcomes"] if o["outcome_id"] == spec["yield"])
        y["yield_type"] = f("isolated yield", "isolated yield", e7)
        y["measurement_basis"] = f(
            "Fig. 4 caption: y. isolated yield", "isolated yield", e7)
    if suffix == "prose_3da_to_7da":
        r["outcomes"] = [
            outcome("yield", "yield", "m_7da", "86% y.", q(86, "%"), e7,
                    product_role="main product", yield_type="isolated yield",
                    basis_raw="Fig. 4 caption: y. isolated yield", basis_norm="isolated yield"),
            outcome("ee", "selectivity", "m_7da", "95% ee", q(95, "%"), e7,
                    product_role="main product",
                    scope={"metric": "ee", "major_isomer": "(S)-7da"}),
        ]
    set_classification(r, spec["classification"], e7,
                       "Partial classification from Fig. 4; full atom mapping, other bond-order changes, and implicit-H changes remain unresolved.")
    r["review"]["unresolved"] = [spec["extra_unresolved"],
                                   "Complete candidate structures and reviewed atom mapping remain unresolved."]
    r["notes"].append(f(
        f"Exact Fig. 4 conditions supersede the earlier page-5 prose limitation for {spec['region']}; the record is not duplicated.",
        f"Exact Fig. 4 conditions supersede the earlier page-5 prose limitation for {spec['region']}; the record is not duplicated.", e7))


# Fig. 4a ORTEP is characterization evidence for (S)-3ak, not a synthesis.
add_occurrence("m_3ak", "(S)-3ak; CCDC 1971578", 7)


# Fig. 5a proposed catalytic cycle is kept separate from experimental records.
e8 = [ev(8, "image"), ev(8, "text")]
mech_parts = [
    participant("allene_family", "m_generic_1", "reactant", e8),
    participant("pyrazolone_family", "m_generic_2", "reactant", e8),
    participant("product_family", "m_generic_3", "product", e8,
                product_role="proposed C-allenylation product"),
    participant("pd_catalyst", "m_[pd]_generic_catalyst", "catalyst", e8),
    participant("ligand", "m_generic_stretchable_ligand", "ligand", e8),
]
new_reaction(
    "fig5a_proposed_mechanism", 8, "proposed_mechanism", ["mechanistic_study"],
    "proposed mechanism for the allenylation of pyrazolone", "proposed Pd catalytic cycle",
    mech_parts, classification=(("C-C", 0, 1, "formed"),), location="Fig. 5a",
    unresolved=("This is an author-proposed catalytic cycle, not a standalone preparative experiment.",
                "Intermediates Int-I/Int-II and Pd oxidation/coordination states are not encoded as verified molecules."))


# Fig. 5b: replace the page-5 generic 'excess bases' description with the exact two stages.
complex_prep = reactions[f"{DOC}:prose_complex_III_prep"]
add_reaction_evidence(complex_prep, 8, "Fig. 5b")
complex_prep["participants"] = [p for p in complex_prep["participants"] if p["participant_id"] != "bases"]
metal = next(p for p in complex_prep["participants"] if p["participant_id"] == "metal")
set_charge(metal, "amount", "PdCl2 (0.1 mmol)", q(0.1, "mmol"), e8)
ligand = next(p for p in complex_prep["participants"] if p["participant_id"] == "ligand")
set_charge(ligand, "equivalents", "DACH-ZYC-Phos-C1 (1.0 equiv)", q(1.0, "equiv"), e8,
           "m_pdcl2")
et3n = participant("et3n", new_ids["et3n"], "reagent", e8)
set_charge(et3n, "equivalents", "Et3N (2.0 equiv)", q(2.0, "equiv"), e8, "m_pdcl2")
add_participant(complex_prep, et3n)
nah = participant("nah", new_ids["nah"], "reagent", e8)
set_charge(nah, "equivalents", "NaH (4.5 equiv)", q(4.5, "equiv"), e8, "m_pdcl2")
add_participant(complex_prep, nah)
complex_prep["reaction_scale"] = f(
    "PdCl2 (0.1 mmol)", {"reference_molecule_id": "m_pdcl2", "amount": q(0.1, "mmol")}, e8)
complex_prep["steps"] = []
for sid, operation, pids, temp, duration in [
    ("et3n_treatment", "Et3N treatment in toluene", ["metal", "ligand", "et3n", "solvent"],
     ("60 °C", q(60, "°C")), ("28 h", q(28, "h"))),
    ("nah_treatment", "NaH treatment", ["metal", "ligand", "et3n", "nah", "solvent"],
     None, ("27.5 h", q(27.5, "h"))),
]:
    step = st.new_step(sid, e8)
    step["participant_ids"] = pids
    step["operation"] = f(operation, operation, e8)
    if temp:
        step["temperature"] = f(temp[0], temp[1], e8)
    step["time"] = f(duration[0], duration[1], e8)
    step["atmosphere"] = f("N2", "nitrogen", e8)
    complex_prep["steps"].append(step)
recryst = st.new_step("recrystallization", [ev(5, "text"), ev(8, "image")], stage="purification")
recryst["participant_ids"] = ["complex", "recrystallization_chcl3", "recrystallization_hexane"]
recryst["operation"] = f(
    "recrystallized from CHCl3/n-hexane", "recrystallization from CHCl3/n-hexane",
    [ev(5, "text"), ev(8, "image")])
complex_prep["steps"].append(recryst)
complex_prep["procedure"]["original_text"] = f(
    "1) Et3N (2.0 equiv), toluene, 60 °C, 28 h; 2) NaH (4.5 equiv), 27.5 h, N2",
    "Et3N treatment followed by NaH treatment under N2", e8)
complex_prep["outcomes"] = [
    outcome("yield", "yield", "m_complex_III", "64%", q(64, "%"), e8,
            product_role="Pd(II)-DACH-ZYC-Phos-C1 complex III")
]
complex_prep["review"]["unresolved"] = [
    "The temperature of the NaH stage, solvent volume, and yield measurement basis are not reported in Fig. 5b.",
    "The coordination graph and complete stereochemical representation remain unencoded.",
]
complex_prep["notes"].append(f(
    "Page 5 says 'excess amount of bases'; Fig. 5b identifies Et3N (2.0 equiv) then NaH (4.5 equiv). Both sources are retained without treating the prose as a separate experiment.",
    "Page 5 says 'excess amount of bases'; Fig. 5b identifies Et3N (2.0 equiv) then NaH (4.5 equiv). Both sources are retained without treating the prose as a separate experiment.",
    [ev(5, "text"), ev(8, "image")]))


# Fig. 5c supplies crystal-structure identity evidence for complexes III and IV.
add_occurrence("m_complex_III", "Pd(II)-DACH-ZYC-Phos-C1·CHCl3; CCDC 2023084", 8)
add_occurrence("m_complex_IV", "Pd(II)-DACH-Phenyl Trost ligand complex IV; CCDC 2031404", 8)
complex_iv = reactions[f"{DOC}:prose_complex_IV_prep"]
add_reaction_evidence(complex_iv, 8, "Fig. 5c ORTEP of complex IV")
complex_iv["review"]["unresolved"] = [
    "Fig. 5c is crystal characterization; the preparation amounts, base identity, time, temperature, and yield remain unreported in the supplied main text.",
    "The complete coordination graph and stereochemical representation remain unencoded.",
]


# Fig. 5d: two true low-activity controls report trace product and recovered 1a.
control_specs = {
    "prose_complex_III_control_standard": (0.2, None, 22, 96),
    "prose_complex_III_control_acoh": (0.1, 20, 25.5, 93),
}
for suffix, (scale, ac_pct, hours, recovery) in control_specs.items():
    r = reactions[f"{DOC}:{suffix}"]
    add_reaction_evidence(r, 8, "Fig. 5d")
    r["reaction_name"] = f(
        "complex-III low-activity allenylation control",
        "complex-III low-activity allenylation control", e8)
    r["transformation_type"] = f(
        "trace C-allenylation under control conditions",
        "trace C-allenylation under control conditions", e8,
        method="inference",
        derivation="Fig. 5d explicitly shows trace (S)-3aa rather than a strictly zero-product result.")
    allene = next(p for p in r["participants"] if p["participant_id"] == "allene")
    set_charge(allene, "amount", f"1a ({scale} mmol)", q(scale, "mmol"), e8)
    pyrazolone = next(p for p in r["participants"] if p["participant_id"] == "pyrazolone")
    set_charge(pyrazolone, "equivalents", "2a (1.2 equiv)", q(1.2, "equiv"), e8, "m_1a")
    complex_p = next(p for p in r["participants"] if p["participant_id"] == "complex")
    set_charge(complex_p, "mol_percent", "III (5.0 mol%)", q(5.0, "mol%"), e8, "m_1a")
    if ac_pct is not None:
        ac = next(p for p in r["participants"] if p["participant_id"] == "acoh")
        set_charge(ac, "mol_percent", "AcOH (20 mol%)", q(20, "mol%"), e8, "m_1a")
    add_participant(r, participant("product_3aa", "m_3aa", "product", e8,
                                   product_role="trace C-allenylation product"))
    r["reaction_scale"] = f(
        f"1a ({scale} mmol)", {"reference_molecule_id": "m_1a", "amount": q(scale, "mmol")}, e8)
    r["solvent_conditions"]["composition"] = f("toluene", ["toluene"], e8)
    r["solvent_conditions"]["concentration"] = f("0.02 M", q(0.02, "M"), e8)
    r["solvent_conditions"]["concentration_basis"] = f(
        "based on 1a", "m_1a", e8, method="inference",
        derivation="The allene 1a is the stated reaction-scale substrate in Fig. 5d.")
    step = st.new_step("reaction", e8)
    step["participant_ids"] = [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"]
    step["operation"] = f("reaction with isolated complex III", "reaction with isolated complex III", e8)
    step["temperature"] = f("60 °C", q(60, "°C"), e8)
    step["time"] = f(f"{hours} h", q(hours, "h"), e8)
    step["atmosphere"] = f("N2", "nitrogen", e8)
    r["steps"] = [step]
    r["outcomes"] = [
        outcome("3aa_trace", "trace", "m_3aa", "(S)-3aa trace", "trace", e8,
                product_role="C-allenylation product"),
        outcome("1a_recovery", "recovery", "m_1a", f"recovery of 1a = {recovery}%",
                q(recovery, "%"), e8, product_role="recovered starting material"),
    ]
    set_classification(r, (("C-C", 0, 1, "formed"),), e8,
                       "A trace C-allenylation product is shown, but no full atom mapping or detection method is reported.")
    r["review"]["unresolved"] = [
        "The trace-product and recovery measurement methods/detection limit are not reported in Fig. 5d.",
        "Full atom mapping, leaving-group changes, and implicit-H changes remain unresolved.",
    ]
    r["notes"].append(f(
        "The page-5 prose says complex III 'failed to catalyze'; Fig. 5d more specifically reports trace (S)-3aa and recovered 1a, so the quantitative figure result is retained rather than a zero-product no_reaction outcome.",
        "The page-5 prose says complex III 'failed to catalyze'; Fig. 5d more specifically reports trace (S)-3aa and recovered 1a, so the quantitative figure result is retained rather than a zero-product no_reaction outcome.",
        [ev(5, "text"), ev(8, "image")]))


# Fig. 7b: a distinct DACH-Phenyl-Trost-ligated experimental comparison.
add_occurrence("m_L11", "DACH-Phenyl-Trost ligand (L11)", 10)
molecules["m_L11"]["names"]["aliases"] = f(
    "DACH-Phenyl-Trost ligand", ["DACH-Phenyl-Trost ligand", "L11"],
    [ev(3, "image"), ev(10, "image")])
add_occurrence("m_3ba", "(S)-3ba", 10)
e10 = [ev(10, "image"), ev(10, "text")]
trost_parts = []
p = participant("allene_1b", "m_1b", "reactant", e10)
set_charge(p, "amount", "1b (0.2 mmol)", q(0.2, "mmol"), e10)
trost_parts.append(p)
p = participant("pyrazolone_2a", "m_2a", "reactant", e10)
set_charge(p, "equivalents", "2a (1.2 equiv)", q(1.2, "equiv"), e10, "m_1b")
trost_parts.append(p)
p = participant("pd_catalyst", "m_pd2(dba)3_chcl3", "catalyst", e10)
set_charge(p, "mol_percent", "Pd2(dba)3·CHCl3 (2.5 mol%)", q(2.5, "mol%"), e10, "m_1b")
trost_parts.append(p)
p = participant("ligand", "m_L11", "ligand", e10)
set_charge(p, "mol_percent", "DACH-Phenyl Trost ligand (5.0 mol%)", q(5.0, "mol%"), e10, "m_1b")
trost_parts.append(p)
trost_parts.append(participant("toluene", "m_toluene", "solvent", e10))
trost_parts.append(participant("product_3ba", "m_3ba", "product", e10,
                               product_role="C-allenylation product"))
trost_outcomes = [
    outcome("nmr_yield", "yield", "m_3ba", "96% y.", q(96, "%"), e10,
            product_role="C-allenylation product", yield_type="NMR yield",
            basis_raw="Fig. 7 caption: y. NMR yield", basis_norm="NMR yield"),
    outcome("ee", "selectivity", "m_3ba", "71% ee", q(71, "%"), e10,
            product_role="C-allenylation product",
            scope={"metric": "ee", "major_isomer": "(S)-3ba"}),
]
trost_rx = new_reaction(
    "fig7b_trost_control", 10, "specific_experiment", ["control", "mechanistic_study"],
    "allenylation of 1b and 2a with DACH-Phenyl-Trost ligand", "C-allenylation",
    trost_parts, outcomes=trost_outcomes, temperature=("60 °C", q(60, "°C")),
    time=("12 h", q(12, "h")), solvent="toluene",
    concentration=("0.02 M", q(0.02, "M"), "m_1b"),
    scale=("1b (0.2 mmol)", "m_1b", q(0.2, "mmol")),
    classification=(("C-C", 0, 1, "formed"),), location="Fig. 7b",
    unresolved=("The ee measurement method is not stated in Fig. 7b.",
                "Full atom mapping, leaving-group changes, and implicit-H changes remain unresolved."))
trost_rx["general_scheme_ids"] = [f"{DOC}:fig1c_generic_allenylation"]


# Page 11 Methods: same 0.5 mmol experiment as Table 2 entry 1, merged rather than duplicated.
e11 = [ev(11, "image"), ev(11, "text")]
method_rx = reactions[f"{DOC}:table2_entry_01"]
add_reaction_evidence(method_rx, 11, "Methods: general procedure and characterization of (S)-3aa")
method_rx["record_type"] = "specific_experiment"
method_rx["experiment_categories"] = f(
    ["substrate_scope"], ["substrate_scope"], [ev(6, "image"), ev(11, "text")])
charges = {
    "allene": {"mass": ("1a (131.5 mg)", q(131.5, "mg")),
               "amount": ("1a (0.5 mmol)", q(0.5, "mmol"))},
    "pyrazolone": {"mass": ("2a (158.7 mg)", q(158.7, "mg")),
                   "amount": ("2a (0.6 mmol)", q(0.6, "mmol")),
                   "equivalents": ("2a (1.2 equiv)", q(1.2, "equiv"))},
    "pd_catalyst": {"mass": ("Pd2(dba)3·CHCl3 (12.9 mg)", q(12.9, "mg")),
                    "amount": ("Pd2(dba)3·CHCl3 (0.0125 mmol)", q(0.0125, "mmol")),
                    "mol_percent": ("2.5 mol%", q(2.5, "mol%"))},
    "ligand": {"mass": ("(R,R)-DACH-ZYC-Phos-C1 (14.3 mg)", q(14.3, "mg")),
               "amount": ("(R,R)-DACH-ZYC-Phos-C1 (0.025 mmol)", q(0.025, "mmol")),
               "mol_percent": ("5.0 mol%", q(5.0, "mol%"))},
    "solvent": {"volume": ("toluene: 5 mL + 11.7 mL + 8.3 mL", q(25.0, "mL"))},
}
for pid, values in charges.items():
    p = next(x for x in method_rx["participants"] if x["participant_id"] == pid)
    for key, (raw, norm) in values.items():
        set_charge(p, key, raw, norm, e11, "m_1a" if key in {"equivalents", "mol_percent"} else None)
method_rx["reaction_scale"] = f(
    "1a (131.5 mg, 0.5 mmol)",
    {"reference_molecule_id": "m_1a", "amount": q(0.5, "mmol")}, e11)
method_rx["solvent_conditions"]["composition"] = f("toluene", ["toluene"], e11)
method_rx["solvent_conditions"]["total_volume"] = f(
    {"Table 2": "25 mL", "Methods": "5 mL + 11.7 mL + 8.3 mL"}, q(25.0, "mL"),
    [ev(6, "image"), ev(11, "text")], method="inference",
    derivation="The three explicitly reported Methods portions sum to the 25 mL total stated in Table 2.")

workup = participant("workup_ethyl_acetate", new_ids["etoac"], "workup_agent", e11, stage="workup")
set_charge(workup, "volume", "ethyl acetate (15 mL × 3)",
           q(15, "mL", repetition_count=3, raw_text="15 mL × 3"), e11,
           method="structure_derivation",
           derivation="The normalized quantity is the explicit 15 mL aliquot; raw_text and repetition_count preserve ×3 without inventing a separately reported total.")
add_participant(method_rx, workup)
add_participant(method_rx, participant("filter_silica", new_ids["silica"], "workup_agent", e11, stage="workup"))
add_participant(method_rx, participant("purification_silica", new_ids["silica"], "purification_agent", e11, stage="purification"))
add_participant(method_rx, participant("petroleum_ether", new_ids["petether"], "purification_agent", e11, stage="purification"))
add_participant(method_rx, participant("purification_ethyl_acetate", new_ids["etoac"], "purification_agent", e11, stage="purification"))
add_participant(method_rx, participant("hplc_hexane", "m_n_hexane", "solvent", e11, stage="characterization"))
add_participant(method_rx, participant("hplc_iproh", new_ids["iproh"], "solvent", e11, stage="characterization"))

method_rx["steps"] = []
for sid, stage, operation, pids, temp, duration, atmosphere in [
    ("sequential_addition", "reaction",
     "Sequential addition to a flame-dried Schlenk flask at room temperature",
     ["pd_catalyst", "ligand", "solvent", "allene", "pyrazolone"], "rt", None, "nitrogen"),
    ("heated_reaction", "reaction",
     "Stir in a pre-heated oil bath; monitor by TLC (petroleum ether/ethyl acetate = 20/1)",
     ["pd_catalyst", "ligand", "solvent", "allene", "pyrazolone", "product"],
     q(60, "°C"), q(12, "h"), "nitrogen"),
    ("silica_filtration", "workup",
     "Filter through a short silica-gel column and elute with ethyl acetate (15 mL × 3)",
     ["filter_silica", "workup_ethyl_acetate"], None, None, None),
    ("evaporation", "workup", "Evaporate the eluate", [], None, None, None),
    ("chromatography_round_1", "purification",
     "Silica-gel chromatography, petroleum ether (60–90 °C)/ethyl acetate = 60/1",
     ["purification_silica", "petroleum_ether", "purification_ethyl_acetate", "product"],
     None, None, None),
    ("chromatography_round_2", "purification",
     "Further silica-gel chromatography of the impure part, petroleum ether (60–90 °C)/ethyl acetate = 40/1",
     ["purification_silica", "petroleum_ether", "purification_ethyl_acetate", "product"],
     None, None, None),
]:
    step = st.new_step(sid, e11, stage=stage)
    step["participant_ids"] = pids
    step["operation"] = f(operation, operation, e11)
    if temp is not None:
        step["temperature"] = f("rt" if temp == "rt" else "60 °C", temp, e11)
    if duration is not None:
        step["time"] = f("12 h", duration, e11)
    if atmosphere is not None:
        step["atmosphere"] = f("under nitrogen atmosphere", atmosphere, e11)
    method_rx["steps"].append(step)

method_text = (
    "Flame-dried Schlenk flask; sequential addition of Pd2(dba)3·CHCl3 (12.9 mg, 0.0125 mmol), "
    "(R,R)-DACH-ZYC-Phos-C1 (14.3 mg, 0.025 mmol), toluene (5 mL), 1a (131.5 mg, 0.5 mmol) "
    "in toluene (11.7 mL), and 2a (158.7 mg, 0.6 mmol) in toluene (8.3 mL) under N2 at rt; "
    "60 °C for 12 h; silica filtration/EtOAc elution; evaporation; two rounds of silica chromatography."
)
method_rx["procedure"]["original_text"] = f(method_text, method_text, e11)
method_rx["procedure"]["addition_order"] = f(
    "Pd catalyst, ligand, toluene, 1a/toluene, then 2a/toluene sequentially",
    ["Pd2(dba)3·CHCl3", "(R,R)-DACH-ZYC-Phos-C1", "toluene", "1a/toluene", "2a/toluene"], e11)
method_rx["procedure"]["pretreatment"] = f(
    "flame-dried Schlenk flask; pre-heated oil bath", ["flame-dried Schlenk flask", "pre-heated oil bath"], e11)
for o in method_rx["outcomes"]:
    if o["outcome_id"] == "yield":
        o["value"] = f("(S)-3aa (179.5 mg, 96%)", q(96, "%"), e11)
        o["yield_type"] = f("isolated yield", "isolated yield", e11)
        o["measurement_basis"] = f(
            "Two-round chromatography afforded (S)-3aa (179.5 mg, 96%)",
            "isolated after two rounds of silica-gel chromatography", e11)
    elif o["outcome_id"] == "ee":
        o["measurement_basis"] = f(
            "Chiralcel IA, n-hexane/i-PrOH = 90/10, 1.0 mL/min, λ = 254 nm; tR major 5.7 min, minor 8.9 min",
            {"method": "chiral HPLC", "column": "Chiralcel IA", "mobile_phase": "n-hexane/i-PrOH 90/10",
             "flow_rate": q(1.0, "mL/min"), "wavelength_nm": 254,
             "retention_times_min": {"major": 5.7, "minor": 8.9}}, e11)
method_rx["outcomes"].append(
    outcome("product_amount", "product_amount", "m_3aa", "179.5 mg", q(179.5, "mg"), e11,
            product_role="C-allenylation product",
            basis_raw="isolated after two rounds of silica-gel chromatography",
            basis_norm="isolated product mass"))
method_rx["review"]["unresolved"] = [
    "Complete product SMILES, reviewed atom mapping, leaving-group changes, and implicit-H changes remain unresolved.",
    "The main-text Methods provide one representative procedure only; unavailable SI details are not inferred for other Table 2 entries.",
]
method_rx["notes"].append(f(
    "The page-11 Methods procedure matches Table 2 entry 1 (1a/2a, 0.5 mmol, 60 °C, 12 h, 96% yield, 95% ee) and is merged into that record rather than counted as a new experiment.",
    "The page-11 Methods procedure matches Table 2 entry 1 (1a/2a, 0.5 mmol, 60 °C, 12 h, 96% yield, 95% ee) and is merged into that record rather than counted as a new experiment.",
    [ev(6, "image"), ev(11, "text")]))
add_occurrence("m_3aa", "(S)-3aa; Methods characterization", 11)


# Non-synthetic source context for figures, characterization, references, and administrative text.
contexts = [
    {"context_id": f"{DOC}:page7_3ak_ortep", "document_id": DOC, "kind": "crystal_structure_characterization",
     "evidence_ids": [ev(7, "image"), ev(7, "text")],
     "content": "Fig. 4a ORTEP representation of (S)-3ak, CCDC 1971578; not a separate synthesis."},
    {"context_id": f"{DOC}:page8_complex_ortep", "document_id": DOC, "kind": "crystal_structure_characterization",
     "evidence_ids": [ev(8, "image"), ev(8, "text")],
     "content": "Fig. 5c ORTEP structures of complexes III·CHCl3 and IV (CCDC 2023084 and 2031404); characterization, not separate reactions."},
    {"context_id": f"{DOC}:pages7_10_dft", "document_id": DOC, "kind": "computational_mechanistic_study",
     "evidence_ids": [ev(7, "text"), ev(9, "image"), ev(9, "text"), ev(10, "image"), ev(10, "text")],
     "content": "DFT transition-state comparisons for DACH-ZYC-Phos-C1 and DACH-Phenyl-Trost catalysts; calculated structures/energies are not preparative experiments."},
    {"context_id": f"{DOC}:page11_3aa_characterization", "document_id": DOC, "kind": "compound_characterization",
     "evidence_ids": [ev(11, "image"), ev(11, "text")],
     "content": "(S)-3aa HPLC, optical rotation, 1H/13C NMR, IR, EI-MS, and HRMS data. Spectral peaks are linked context, not standalone syntheses."},
    {"context_id": f"{DOC}:pages11_12_back_matter", "document_id": DOC, "kind": "data_availability_references_and_back_matter",
     "evidence_ids": [ev(11, "image"), ev(11, "text"), ev(12, "image"), ev(12, "text")],
     "content": "Data-availability statement, references, acknowledgements, author contributions, competing interests, and licensing. SI is mentioned online but was not supplied/associated in this workspace."},
]
known_contexts = {x.get("context_id") for x in bundle["source_context"]}
bundle["source_context"].extend(x for x in contexts if x["context_id"] not in known_contexts)


# Coverage and visual inventory for exactly the requested six-page segment.
page_records = {
    7: [
        f"{DOC}:fig4_scaleup_3am",
        f"{DOC}:prose_3am_to_5am",
        f"{DOC}:prose_3da_to_6da",
        f"{DOC}:prose_3da_to_7da",
    ],
    8: [
        f"{DOC}:fig5a_proposed_mechanism",
        f"{DOC}:prose_complex_III_prep",
        f"{DOC}:prose_complex_IV_prep",
        f"{DOC}:prose_complex_III_control_standard",
        f"{DOC}:prose_complex_III_control_acoh",
    ],
    9: [],
    10: [f"{DOC}:fig7b_trost_control"],
    11: [f"{DOC}:table2_entry_01"],
    12: [],
}
coverage_regions = {
    7: [
        {"region": "Fig. 4a", "finding": "ORTEP representation of (S)-3ak (CCDC 1971578); crystal characterization only.",
         "evidence_ids": [ev(7, "image"), ev(7, "text")]},
        {"region": "Fig. 4b reaction i", "finding": "Distinct 4.0 mmol gram-scale allenylation to (S)-3am; recorded as scale-up experiment.",
         "evidence_ids": [ev(7, "image"), ev(7, "text")]},
        {"region": "Fig. 4b reactions ii-iv and caption", "finding": "Three synthetic applications with exact catalyst/reagent loads, temperatures, times, isolated yields, and selectivities; merged with page-5 prose records.",
         "evidence_ids": [ev(7, "image"), ev(7, "text")]},
        {"region": "lower two-column prose", "finding": "DFT comparison/summary; no additional preparative experiment beyond Fig. 7b.",
         "evidence_ids": [ev(7, "image"), ev(7, "text")]},
    ],
    8: [
        {"region": "Fig. 5a", "finding": "Author-proposed Pd catalytic cycle; stored as proposed_mechanism, not experimental proof.",
         "evidence_ids": [ev(8, "image"), ev(8, "text")]},
        {"region": "Fig. 5b", "finding": "Two-stage synthesis of Pd(II)-DACH-ZYC-Phos-C1 complex III with 64% reported yield; merged with page-5 prose record.",
         "evidence_ids": [ev(8, "image"), ev(8, "text")]},
        {"region": "Fig. 5c", "finding": "ORTEP structures of complexes III and IV; identity/characterization evidence, not separate syntheses.",
         "evidence_ids": [ev(8, "image"), ev(8, "text")]},
        {"region": "Fig. 5d", "finding": "Two distinct low-activity mechanistic controls; both show trace (S)-3aa and separately report 1a recovery.",
         "evidence_ids": [ev(8, "image"), ev(8, "text")]},
    ],
    9: [
        {"region": "full page Fig. 6", "finding": "DFT-optimized transition states, energies, hydrogen-bond distances/angles, and bond energies; computational context only.",
         "evidence_ids": [ev(9, "image"), ev(9, "text")]},
    ],
    10: [
        {"region": "Fig. 7a", "finding": "DFT transition states for DACH-Phenyl-Trost catalyst; computational context only.",
         "evidence_ids": [ev(10, "image"), ev(10, "text")]},
        {"region": "Fig. 7b and caption", "finding": "Distinct experimental comparison using DACH-Phenyl-Trost ligand: 96% NMR yield and 71% ee (S)-3ba.",
         "evidence_ids": [ev(10, "image"), ev(10, "text")]},
    ],
    11: [
        {"region": "Methods, left column", "finding": "General procedure and workup/purification for (S)-3aa; matches and enriches Table 2 entry 1 rather than creating a duplicate.",
         "evidence_ids": [ev(11, "image"), ev(11, "text")]},
        {"region": "(S)-3aa characterization block", "finding": "HPLC/rotation/NMR/IR/MS/HRMS characterization; retained as source context, not separate synthesis records.",
         "evidence_ids": [ev(11, "image"), ev(11, "text")]},
        {"region": "Data availability and references 1-34", "finding": "Back matter and cited literature only; online SI is mentioned but was not supplied or associated.",
         "evidence_ids": [ev(11, "image"), ev(11, "text")]},
    ],
    12: [
        {"region": "references 35-57", "finding": "Cited literature list only; no current-paper experimental record.",
         "evidence_ids": [ev(12, "image"), ev(12, "text")]},
        {"region": "acknowledgements through licence", "finding": "Administrative/back-matter content only; no reaction or standalone characterization experiment.",
         "evidence_ids": [ev(12, "image"), ev(12, "text")]},
    ],
}
page_unresolved = {
    7: ["Full structures/atom mappings for 5am, 6da, and 7da are not safely encoded; selectivity measurement methods are not stated in Fig. 4."],
    8: ["Complex III NaH-stage temperature, solvent volume, yield basis, and trace/recovery measurement methods are unreported; coordination graphs remain unencoded."],
    9: ["No preparative experiment is present; computational structures are not converted into molecule-registry identities."],
    10: ["Fig. 7b states NMR yield but not the ee measurement method; full atom mapping remains unresolved."],
    11: ["Only the representative 3aa procedure is in the supplied main text; SI is unavailable and not inferred for other products."],
    12: ["No synthetic experiment is present on this page."],
}
coverage = {"document_id": DOC, "semantic_accuracy_verified": False, "pages": []}
for page in range(7, 13):
    coverage["pages"].append({
        "page": page, "state": "attempted", "regions": coverage_regions[page],
        "record_ids": page_records[page], "unresolved": page_unresolved[page],
    })
(RUN / "coverage/pages_0007_0012.json").write_text(
    json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")

inventory = {
    "document_id": DOC,
    "source_review_mode": "manual visual inspection of original-detail rendered page images",
    "pages": {
        "7": ["Fig. 4a ORTEP", "Fig. 4b scale-up and three synthetic applications", "DFT/summary prose"],
        "8": ["Fig. 5a proposed mechanism", "Fig. 5b complex III synthesis", "Fig. 5c ORTEPs", "Fig. 5d two controls"],
        "9": ["full-page Fig. 6 DFT transition states"],
        "10": ["Fig. 7a DFT transition states", "Fig. 7b experimental comparison"],
        "11": ["Methods/general procedure", "3aa characterization", "data availability", "references 1-34"],
        "12": ["references 35-57", "acknowledgements and article back matter"],
    },
    "reaction_record_ids_by_page": {str(k): v for k, v in page_records.items()},
    "semantic_accuracy_verified": False,
}
(RUN / "raw_model/visual_inventory_pages_0007_0012.json").write_text(
    json.dumps(inventory, ensure_ascii=False, indent=2) + "\n")


# Machine records must remain unreviewed/needs_review.
for r in bundle["reactions"]:
    r["review"]["status"] = "needs_review"
    r["review"]["method"] = "model"
    r["review"]["reviewer"] = None
    r["review"]["checked_fields"] = []

st.write(SUBMITTED, bundle)
print(json.dumps({
    "molecules": len(bundle["molecules"]),
    "reactions": len(bundle["reactions"]),
    "new_reactions": [f"{DOC}:fig4_scaleup_3am", f"{DOC}:fig5a_proposed_mechanism",
                      f"{DOC}:fig7b_trost_control"],
    "pages": {str(k): len(v) for k, v in page_records.items()},
}, ensure_ascii=False, indent=2))
