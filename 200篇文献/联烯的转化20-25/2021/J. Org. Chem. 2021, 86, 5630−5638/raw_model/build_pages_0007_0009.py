#!/usr/bin/env python3
"""Extend the v3.1 source bundle with the page 7--9 extraction.

The input is the existing machine bundle for pages 1--6.  Existing reaction
identities are retained and receive the additional page evidence; only the
preparation of Intermediate F is a newly created experiment.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_b27cb890048f595e"
E = {p: f"{DOC}:page{p}:text" for p in range(1, 10)}
I = {p: f"{DOC}:page{p}:image" for p in range(1, 10)}

# This stage consumes the pages-1--6 source bundle produced by the preceding
# extraction stage.  Re-running that stage first makes this transform
# reproducible without accumulating duplicate workup participants or steps.
bundle = json.loads((RUN / "submitted.json").read_text())
# ``standard init`` is the authority for the complete document evidence
# registry.  The previously merged pages-1--6 bundle intentionally contained
# only the evidence that had then been processed.
initialized = json.loads((RUN / "raw_model/pages_0007_0009_init.json").read_text())
bundle["documents"] = initialized["documents"]
molecules = {m["molecule_id"]: m for m in bundle["molecules"]}
reactions = {r["reaction_id"]: r for r in bundle["reactions"]}


def f(raw=None, normalized=None, ev=(), *, method="original_extraction",
      derivation=None, status=None, verification=None, alternatives=()):
    if normalized is None and raw is not None:
        normalized = raw
    return st.fact(
        raw,
        normalized,
        ev,
        method=method,
        derivation=derivation,
        status=status,
        verification=verification,
        alternatives=alternatives,
    )


def inferred(normalized, ev, derivation):
    return f(None, normalized, ev, method="inference", derivation=derivation)


def derived(normalized, ev, derivation):
    return f(None, normalized, ev, method="structure_derivation", derivation=derivation)


def add_evidence_recursive(obj, evidence_ids):
    """Attach a page image only where that page's text already supports a fact."""
    if isinstance(obj, dict):
        if isinstance(obj.get("evidence_ids"), list) and obj["evidence_ids"]:
            for eid in evidence_ids:
                page = next((p for p, image_id in I.items() if image_id == eid), None)
                if page is not None and E[page] in obj["evidence_ids"] and eid not in obj["evidence_ids"]:
                    obj["evidence_ids"].append(eid)
        for value in obj.values():
            add_evidence_recursive(value, evidence_ids)
    elif isinstance(obj, list):
        for value in obj:
            add_evidence_recursive(value, evidence_ids)


def drop_evidence_recursive(obj, evidence_ids):
    if isinstance(obj, dict):
        if isinstance(obj.get("evidence_ids"), list):
            obj["evidence_ids"] = [x for x in obj["evidence_ids"] if x not in evidence_ids]
        for value in obj.values():
            drop_evidence_recursive(value, evidence_ids)
    elif isinstance(obj, list):
        for value in obj:
            drop_evidence_recursive(value, evidence_ids)


def add_record_evidence(reaction_id, evidence_ids, drop=()):
    r = reactions[reaction_id]
    if drop:
        drop_evidence_recursive(r, set(drop))
    for eid in evidence_ids:
        if eid not in r["evidence_ids"]:
            r["evidence_ids"].append(eid)
    add_evidence_recursive(r, tuple(eid for eid in evidence_ids if eid in I.values()))


def molecule_evidence(molecule_id, evidence_ids):
    m = molecules[molecule_id]
    occurrence = m["occurrences"][0]
    for eid in evidence_ids:
        if eid not in occurrence["evidence_ids"]:
            occurrence["evidence_ids"].append(eid)
    add_evidence_recursive(m, tuple(eid for eid in evidence_ids if eid in I.values()))


def charge_fact(raw, normalized, ev, *, method="original_extraction", derivation=None):
    return f(raw, normalized, ev, method=method, derivation=derivation)


def participant(pid, mid, role, ev, *, stage="reaction", product_role=None,
                mass=None, amount=None, volume=None, equivalents=None,
                mol_percent=None, reference=None):
    p = st.new_participant(pid, mid, role, ev, stage)
    if product_role:
        p["product_role"] = f(product_role, product_role, ev)
    for key, pair in {
        "mass": mass,
        "amount": amount,
        "volume": volume,
        "equivalents": equivalents,
        "mol_percent": mol_percent,
    }.items():
        if pair is not None:
            p["charges"][key] = charge_fact(pair[0], pair[1], ev)
    if reference:
        p["equivalence_reference"] = f(reference, reference, ev)
    return p


def new_step(sid, ev, pids, operation, *, stage="reaction", temperature=None,
             time=None, atmosphere=None, operation_method="original_extraction",
             operation_derivation=None):
    s = st.new_step(sid, ev, stage)
    s["participant_ids"] = list(pids)
    s["operation"] = f(
        operation,
        operation,
        ev,
        method=operation_method,
        derivation=operation_derivation,
    )
    if temperature:
        s["temperature"] = f(temperature[0], temperature[1], ev)
    if time:
        s["time"] = f(time[0], time[1], ev)
    if atmosphere:
        raw, normalized, method, derivation, at_ev = atmosphere
        s["atmosphere"] = f(
            raw,
            normalized,
            at_ev,
            method=method,
            derivation=derivation,
        )
    return s


def outcome(oid, kind, target, raw, normalized, ev, *, basis=None,
            yield_type=None, scope=None):
    o = st.new_outcome(oid, kind, [target] if target else [], ev)
    o["value"] = f(raw, normalized, ev)
    if kind in {"yield", "product_amount"}:
        o["product_role"] = f("main product", "main_product", ev)
    if yield_type:
        o["yield_type"] = f(yield_type, yield_type, ev)
    if basis:
        o["measurement_basis"] = f(basis, basis, ev)
    if scope:
        o["scope"] = f(scope[0], scope[1], ev)
    return o


def set_reactant_charge(r, pid, mass_raw, mass, amount_raw, amount, ev):
    p = next(x for x in r["participants"] if x["participant_id"] == pid)
    p["charges"]["mass"] = charge_fact(mass_raw, {"number": mass, "unit": "mg"}, ev)
    p["charges"]["amount"] = charge_fact(amount_raw, {"number": amount, "unit": "mmol"}, ev)
    r["reaction_scale"] = f(
        amount_raw + f" of {p['molecule_id']}",
        {"reference_molecule_id": p["molecule_id"], "amount": {"number": amount, "unit": "mmol"}},
        ev,
    )


def append_participant(r, p):
    if all(existing["participant_id"] != p["participant_id"] for existing in r["participants"]):
        r["participants"].append(p)


def add_workup(r, ev, ratio):
    entries = [
        participant("extract_etoac", "ethyl_acetate", "workup_agent", ev, stage="workup"),
        participant("workup_water", "water", "workup_agent", ev, stage="workup"),
        participant("dry_na2so4", "sodium_sulfate", "workup_agent", ev, stage="workup"),
        participant("silica_column", "silica", "purification_agent", ev, stage="purification"),
        participant("etoac_eluent", "ethyl_acetate", "purification_agent", ev, stage="purification"),
        participant("hexane_eluent", "hexane", "purification_agent", ev, stage="purification"),
    ]
    for p in entries:
        append_participant(r, p)
    r["steps"] = [s for s in r["steps"] if s["step_id"] not in {"workup", "purification"}]
    r["steps"].append(new_step(
        "workup",
        ev,
        ["extract_etoac", "workup_water", "dry_na2so4"],
        "the solvent was evaporated off, followed by extraction using ethyl acetate and water; the organic layer was dried over Na2SO4 and evaporated under reduced pressure",
        stage="workup",
    ))
    r["steps"].append(new_step(
        "purification",
        ev,
        ["silica_column", "etoac_eluent", "hexane_eluent"],
        f"column chromatography (eluent: ethyl acetate/hexane = {ratio})",
        stage="purification",
    ))


def set_solvent(r, raw, mid, volume, ev):
    r["solvent_conditions"]["composition"] = f(raw, raw, ev)
    r["solvent_conditions"]["total_volume"] = f(
        f"{volume:g} mL",
        {"number": volume, "unit": "mL"},
        ev,
    )
    solvent = next((p for p in r["participants"] if p["molecule_id"] == mid and p["stage"] == "reaction"), None)
    if solvent:
        solvent["charges"]["volume"] = f(
            f"{volume:g} mL",
            {"number": volume, "unit": "mL"},
            ev,
        )


def set_procedure(r, text, ev):
    r["procedure"]["original_text"] = f(text, text, ev)


def set_candidate_structure(mid, smiles, ev, *, chirality, geometry,
                            representation="partial", note=None, mixture=None,
                            assignment=None):
    m = molecules[mid]
    diag = st.stereo_inspect(smiles)
    m["representation_type"] = "smiles"
    m["smiles"] = derived(
        smiles,
        ev,
        "Connectivity candidate derived from the evidence-linked systematic name and Scheme 4 drawing; RDKit parsing/formula agreement is diagnostic, not source verification.",
    )
    m["stereochemistry"]["atom_centers"] = derived(
        diag["atom_centers"],
        ev,
        "RDKit potential/assigned-center diagnostic on the machine structure candidate; '?' retains unassigned source stereocenters.",
    ) if diag["atom_centers"] else f()
    m["stereochemistry"]["double_bonds"] = derived(
        diag["double_bonds"],
        ev,
        "RDKit double-bond diagnostic on the machine structure candidate; STEREONONE is not an author assignment.",
    ) if diag["double_bonds"] else f()
    m["stereochemistry"]["axial"] = inferred(
        "not_applicable",
        ev,
        "No axial/allene stereogenic element is present in this product connectivity candidate.",
    )
    m["stereochemistry"]["chirality_status"] = inferred(
        chirality,
        ev,
        "The source reports no absolute configurations; mixture status is retained where a dr is reported.",
    )
    m["stereochemistry"]["geometry_status"] = inferred(
        geometry,
        ev,
        "Geometry applicability/status assessed from the source name and drawing.",
    )
    m["stereochemistry"]["representation_status"] = inferred(
        representation,
        ev,
        "The connectivity is encoded, but unassigned stereocenters and/or mixture composition prevent a fully specified stereoisomer representation.",
    )
    m["stereochemistry"]["isomeric_smiles"] = derived(
        diag["canonical_isomeric_smiles"],
        ev,
        "Canonical isomeric SMILES generated from the machine candidate; unspecified centers remain unspecified.",
    )
    m["stereochemistry"]["identity_key"] = derived(
        diag["stereo_graph_key"],
        ev,
        "Stereo graph key generated from the machine candidate; it does not uniquely identify an unassigned mixture/sample.",
    )
    if mixture:
        m["stereochemistry"]["mixture"] = f(mixture[0], mixture[1], ev)
    if assignment:
        m["stereochemistry"]["assignment"] = f(assignment[0], assignment[1], ev)
    if note:
        m["identity_note"] = f(note, note, ev)


# Add the newly reviewed page images to existing cross-page records rather
# than duplicating the same experiments.
scope_ids = [f"{DOC}:scheme3_scope_4{x}" for x in "cdefgh"]
for rid in scope_ids:
    add_record_evidence(rid, (E[7], I[7]))

for suffix in ("scheme4_4b_to_5", "scheme4_5_to_6", "scheme4_6_to_7"):
    add_record_evidence(f"{DOC}:{suffix}", (E[7], I[7]), drop=(E[8], I[8]))
add_record_evidence(f"{DOC}:scheme4_2a_to_8", (E[7], I[7], E[8], I[8]), drop=(E[8], I[8]))
add_record_evidence(f"{DOC}:scheme4_2a_to_9", (E[8], I[8]), drop=(E[7], I[7]))

for mid in ("4b", "4c", "4d", "4e", "4f", "4g", "4h", "5", "6", "7"):
    molecule_evidence(mid, (E[7], I[7]))
for mid in ("8", "9", "intermediate_f"):
    molecule_evidence(mid, (E[8], I[8]))

# Source-derived structure candidates for the molecules whose complete names
# and drawings are available on pages 4, 7, and 8.  Configurations are not
# invented: unassigned centers remain unspecified and dr products remain
# mixture identities.
ev_4b = (I[4], E[6], I[6], E[7], I[7])
set_candidate_structure(
    "4b",
    "CCOC(=O)C1=C(c2ccccc2)C(=O)NC1(C)c1c[nH]c2ccccc12",
    ev_4b,
    chirality="unspecified",
    geometry="not_applicable",
    note="The C2 tetrahedral center is drawn without absolute configuration; no racemic assignment is inferred.",
)

ev_567 = (I[4], E[7], I[7])
full_names = {
    "5": "1-(tert-Butyl)-3-ethyl-2-(1-(tert-butoxycarbonyl)-1H-indol-3-yl)-2-methyl-5-oxo-4-phenyl-2,5-dihydro-1H-pyrrole-1,3-dicarboxylate",
    "6": "1-(tert-Butyl)-3-ethyl-2-(1-(tert-butoxycarbonyl)-1H-indol-3-yl)-5-hydroxy-2-methyl-4-phenyl-2,5-dihydro-1H-pyrrole-1,3-dicarboxylate",
    "7": "1-(tert-Butyl)-3-ethyl-2-(1-(tert-butoxycarbonyl)-1H-indol-3-yl)-2-methyl-5-(1-methyl-1H-indol-3-yl)-4-phenyl-2,5-dihydro-1H-pyrrole-1,3-dicarboxylate",
}
for mid, name in full_names.items():
    molecules[mid]["names"]["original"] = f(name, name, ev_567)
    molecules[mid]["names"]["english"] = f(name, name, ev_567)

set_candidate_structure(
    "5",
    "CCOC(=O)C1=C(c2ccccc2)C(=O)N(C(=O)OC(C)(C)C)C1(C)c1cn(C(=O)OC(C)(C)C)c2ccccc12",
    ev_567,
    chirality="unspecified",
    geometry="not_applicable",
    note="Connectivity candidate agrees with the reported formula; the C2 center is unassigned.",
)
set_candidate_structure(
    "6",
    "CCOC(=O)C1=C(c2ccccc2)C(O)N(C(=O)OC(C)(C)C)C1(C)c1cn(C(=O)OC(C)(C)C)c2ccccc12",
    ev_567,
    chirality="mixture",
    geometry="not_applicable",
    mixture=("inseparable mixture of diastereomers, 5:1", {"type": "diastereomeric_mixture", "metric": "dr", "reported_ratio": "5:1", "order": "major/minor; configurations unassigned"}),
    note="Two tetrahedral centers are present in the connectivity candidate; configurations and the dr measurement method are unreported.",
)
set_candidate_structure(
    "7",
    "CCOC(=O)C1=C(c2ccccc2)C(c2cn(C)c3ccccc23)N(C(=O)OC(C)(C)C)C1(C)c1cn(C(=O)OC(C)(C)C)c2ccccc12",
    ev_567,
    chirality="mixture",
    geometry="not_applicable",
    mixture=("diastereomeric ratio: 3:1", {"type": "diastereomeric_mixture", "metric": "dr", "reported_ratio": "3:1", "order": "major/minor; configurations unassigned"}),
    note="Two tetrahedral centers are present in the connectivity candidate; configurations and the dr measurement method are unreported.",
)

ev_89 = (I[4], E[7], I[7], E[8], I[8])
set_candidate_structure(
    "8",
    "CCOC(=O)C1=C(c2ccccc2)C(=O)NC1C",
    ev_89,
    chirality="unspecified",
    geometry="not_applicable",
    note="The newly hydrogen-bearing C2 center is stereogenic but no configuration or mixture composition is reported.",
)
set_candidate_structure(
    "9",
    "CCOC(=O)C1=C(c2ccccc2)C(=O)NC1(C)OC",
    (I[4], E[8], I[8]),
    chirality="unspecified",
    geometry="not_applicable",
    note="The methoxy/methyl-bearing C2 center is stereogenic but no configuration or mixture composition is reported.",
)

f_smiles = "CCOC(=O)C1=C(C)O/C(=N\\O)/C1c1ccccc1"
set_candidate_structure(
    "intermediate_f",
    f_smiles,
    (I[4], E[8], I[8]),
    chirality="unspecified",
    geometry="specified",
    representation="partial",
    assignment=("(Z)-Ethyl-5-(hydroxyimino)...", {"type": "geometric", "element": "C=N-OH", "configuration": "Z", "assignment_method": "reported name"}),
    note="The reported Z oxime geometry is encoded; the C4 tetrahedral center remains unassigned.",
)
f_diag = st.stereo_inspect(f_smiles)
molecules["intermediate_f"]["stereochemistry"]["double_bonds"] = derived(
    f_diag["double_bonds"],
    (I[4], E[8], I[8]),
    "The reported (Z) oxime was encoded with directional bonds and checked by RDKit (C=N STEREOZ); other double bonds have no applicable E/Z assignment.",
)
molecules["intermediate_f"]["stereochemistry"]["geometry_status"] = f(
    "(Z)-hydroxyimino",
    "specified",
    (E[8], I[8]),
)

# Complete the source charges, ordered operations, workup, and purification for
# the five existing Scheme 4 experiments described on pages 7--8.
r = reactions[f"{DOC}:scheme4_4b_to_5"]
set_reactant_charge(r, "s", "144 mg", 144, "0.4 mmol", 0.4, (E[7], I[7]))
set_solvent(r, "THF", "thf", 3, (E[7], I[7]))
add_workup(r, (E[7], I[7]), "20/80")
set_procedure(r, "To a mixture of 4b (0.4 mmol, 144 mg), triethylamine (0.4 mmol, 40 mg) and DMAP (0.2 mmol, 24 mg) in THF (3 mL) was added (Boc)2O (1.0 mmol, 218 mg) slowly at 0 °C. The reaction mixture was warmed up to 25 °C and stirred for 8 h; solvent evaporation, ethyl acetate/water extraction, Na2SO4 drying, and ethyl acetate/hexane 20/80 column chromatography afforded 5 (122 mg, 55%).", (E[7], I[7]))

r = reactions[f"{DOC}:scheme4_5_to_6"]
set_reactant_charge(r, "s", "112 mg", 112, "0.20 mmol", 0.20, (E[7], I[7]))
set_solvent(r, "THF", "thf", 3, (E[7], I[7]))
append_participant(r, participant("dibal_carrier", "hexane", "solvent", (E[7], I[7])))
r["steps"] = [
    new_step("low_temperature_addition", (E[7], I[7]), ["s", "red", "solv", "dibal_carrier"], "add a solution of DIBAL-H (0.8 mL, 1.0 M in hexane) to the solution of 5", temperature=("−78 °C", {"number": -78, "unit": "°C"}), atmosphere=("N2", "N2", "original_extraction", None, (E[7], I[7]))),
    new_step("low_temperature_hold", (E[7], I[7]), ["s", "red", "solv", "dibal_carrier", "p"], "stir at the same temperature", temperature=("−78 °C", {"number": -78, "unit": "°C"}), time=("30 min", {"number": 30, "unit": "min"}), atmosphere=("N2", "N2", "original_extraction", None, (E[7], I[7]))),
]
add_workup(r, (E[7], I[7]), "15/85")
set_procedure(r, "A solution of DIBAL-H (0.8 mL, 1.0 M in hexane) was added to 5 (112 mg, 0.20 mmol) in THF (3 mL) under N2 at −78 °C and stirred for 30 min at the same temperature; solvent evaporation, ethyl acetate/water extraction, Na2SO4 drying, and ethyl acetate/hexane 15/85 column chromatography afforded 6 (82 mg, 73%; dr 5:1).", (E[7], I[7]))
dr6 = next(o for o in r["outcomes"] if o["outcome_id"] == "dr_6")
dr6["value"] = f("Diastereomeric ratio: 5:1", {"metric": "dr", "ratio": "5:1", "order": "major/minor; configurations unassigned"}, (E[7], I[7]))
dr6["scope"] = f("Diastereomeric ratio: 5:1", {"metric": "dr", "order": "major/minor", "major_configuration": None, "minor_configuration": None}, (E[7], I[7]))

r = reactions[f"{DOC}:scheme4_6_to_7"]
set_reactant_charge(r, "s", "112 mg", 112, "0.20 mmol", 0.20, (E[7], I[7]))
set_solvent(r, "DCM", "dcm", 2, (E[7], I[7]))
add_workup(r, (E[7], I[7]), "20/80")
set_procedure(r, "DCM (2.0 mL) was added to 6 (112 mg, 0.20 mmol), 1-methylindole (39 mg, 0.30 mmol), and In(OTf)3 (11 mg, 0.02 mmol) under air. The mixture was stirred at 25 °C for 15 min; solvent evaporation, ethyl acetate/water extraction, Na2SO4 drying, and ethyl acetate/hexane 20/80 column chromatography afforded 7 (90 mg, 67%; dr 3:1).", (E[7], I[7]))
dr7 = next(o for o in r["outcomes"] if o["outcome_id"] == "dr_7")
dr7["value"] = f("Diastereomeric ratio: 3:1", {"metric": "dr", "ratio": "3:1", "order": "major/minor; configurations unassigned"}, (E[7], I[7]))
dr7["scope"] = f("Diastereomeric ratio: 3:1", {"metric": "dr", "order": "major/minor", "major_configuration": None, "minor_configuration": None}, (E[7], I[7]))

r = reactions[f"{DOC}:scheme4_2a_to_8"]
set_reactant_charge(r, "s", "52 mg", 52, "0.20 mmol", 0.20, (E[7], I[7]))
set_solvent(r, "THF", "thf", 2, (E[7], I[7]))
r["steps"] = [
    new_step("reductant_addition", (E[7], I[7]), ["s", "red", "solv"], "add NaBH3CN in one portion", temperature=("25 °C", {"number": 25, "unit": "°C"})),
    new_step("pre_acidification_hold", (E[7], I[7]), ["s", "red", "solv"], "hold before acidification", temperature=("25 °C", {"number": 25, "unit": "°C"}), time=("5 min", {"number": 5, "unit": "min"})),
    new_step("acidification_and_hold", (E[7], I[7]), ["s", "red", "acid", "solv", "p"], "acidify with 1.5 equiv of 35% HCl and stir", temperature=("25 °C", {"number": 25, "unit": "°C"}), time=("5 h", {"number": 5, "unit": "h"})),
]
add_workup(r, (E[7], I[7]), "30/70")
set_procedure(r, "NaBH3CN (38 mg, 0.60 mmol) was added in one portion to 2a (52 mg, 0.20 mmol) in THF (2 mL) at 25 °C. After 5 min, the mixture was acidified with 1.5 equiv of 35% HCl and stirred for 5 h; solvent evaporation, ethyl acetate/water extraction, Na2SO4 drying, and ethyl acetate/hexane 30/70 column chromatography afforded 8 (27 mg, 56%).", (E[7], I[7]))

r = reactions[f"{DOC}:scheme4_2a_to_9"]
set_reactant_charge(r, "s", "52 mg", 52, "0.20 mmol", 0.20, (E[8], I[8]))
set_solvent(r, "MeOH", "methanol", 2, (E[8], I[8]))
add_workup(r, (E[8], I[8]), "30/70")
set_procedure(r, "2a (52 mg, 0.20 mmol) in MeOH (2 mL) was acidified with 1 equiv of 35% HCl and stirred for 1 h at 25 °C; solvent evaporation, ethyl acetate/water extraction, Na2SO4 drying, and ethyl acetate/hexane 30/70 column chromatography afforded 9 (39 mg, 70%).", (E[8], I[8]))


# New molecule explicitly introduced on page 8.
if "pbu3" not in molecules:
    m = st.new_molecule("pbu3", DOC, (E[8], I[8]))
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f("PBu3", "PBu3", (E[8], I[8]))
    m["names"]["original"] = f("PBu3", "tri-n-butylphosphine", (E[8], I[8]))
    m["names"]["english"] = f("tributylphosphine", "tributylphosphine", (E[8], I[8]))
    m["names"]["abbreviations"] = f("PBu3", ["PBu3"], (E[8], I[8]))
    m["names"]["aliases"] = f()
    m["cas_numbers"] = f()
    m["representation_type"] = "smiles"
    diag = st.stereo_inspect("CCCCP(CCCC)CCCC")
    m["smiles"] = derived("CCCCP(CCCC)CCCC", (E[8], I[8]), "Machine connectivity candidate derived from the explicit PBu3 identity; not an author-provided SMILES.")
    m["stereochemistry"]["atom_centers"] = f()
    m["stereochemistry"]["double_bonds"] = f()
    m["stereochemistry"]["axial"] = inferred("not_applicable", (E[8], I[8]), "No axial stereogenic element in the connectivity candidate.")
    m["stereochemistry"]["mixture"] = f()
    m["stereochemistry"]["assignment"] = f()
    m["stereochemistry"]["chirality_status"] = inferred("achiral", (E[8], I[8]), "Three identical n-butyl substituents preclude a configurationally stereogenic phosphorus center in this candidate.")
    m["stereochemistry"]["geometry_status"] = inferred("not_applicable", (E[8], I[8]), "No applicable E/Z element.")
    m["stereochemistry"]["representation_status"] = inferred("encoded", (E[8], I[8]), "Connectivity has no source-specified stereochemical element to encode.")
    m["stereochemistry"]["isomeric_smiles"] = derived(diag["canonical_isomeric_smiles"], (E[8], I[8]), "Canonical isomeric SMILES generated from the machine candidate.")
    m["stereochemistry"]["identity_key"] = derived(diag["stereo_graph_key"], (E[8], I[8]), "Stereo graph key generated from the machine candidate.")
    m["identity_note"] = f("The article reports PBu3; CAS is not reported.", "The article reports PBu3; CAS is not reported.", (E[8], I[8]))
    m["external_properties"] = []
    molecules["pbu3"] = m
    bundle["molecules"].append(m)


# Page 8 contains a distinct preparation of Intermediate F that was absent
# from the prior-page bundle.
new_id = f"{DOC}:preparation_intermediate_f"
if new_id not in reactions:
    ev = (E[8], I[8], E[5])
    r = st.new_reaction(new_id, DOC, ev)
    r["record_type"] = "specific_experiment"
    r["experiment_categories"] = f(["mechanistic_study"], ["mechanistic_study"], (E[8], I[8]))
    r["location"] = f({"page": 8, "section": "Experimental Section", "compound": "Intermediate F"}, {"page": 8, "section": "Experimental Section", "compound": "Intermediate F"}, (E[8], I[8]))
    r["general_scheme_ids"] = []
    title = "Synthesis of (Z)-Ethyl-5-(hydroxyimino)-2-methyl-4-phenyl-4,5-dihydrofuran-3-carboxylate (Intermediate F)"
    r["reaction_name"] = f(title, title, (E[8], I[8]))
    r["transformation_type"] = inferred("phosphine-mediated intramolecular cyclization to a (Z)-furanone oxime", (E[8], I[8]), "Transformation interpretation from the named allenoate starting material, named product, and the cited Kwon procedure.")
    r["participants"] = [
        participant("substrate", "1a", "reactant", (E[8], I[8]), mass=("52 mg", {"number": 52, "unit": "mg"}), amount=("0.20 mmol", {"number": 0.20, "unit": "mmol"})),
        participant("pbu3", "pbu3", "catalyst", (E[8], I[8]), mass=("2 mg", {"number": 2, "unit": "mg"}), amount=("0.01 mmol", {"number": 0.01, "unit": "mmol"}), equivalents=("0.05 equiv (calculated)", {"number": 0.05, "unit": "equiv"}), mol_percent=("5 mol% (calculated)", {"number": 5, "unit": "mol%"}), reference="1a"),
        participant("solvent", "acetonitrile", "solvent", (E[8], I[8]), volume=("2 mL", {"number": 2, "unit": "mL"})),
        participant("product", "intermediate_f", "product", (E[8], I[8]), product_role="main_product"),
        participant("silica_column", "silica", "purification_agent", (E[8], I[8], E[5]), stage="purification"),
        participant("etoac_eluent", "ethyl_acetate", "purification_agent", (E[8], I[8], E[5]), stage="purification"),
        participant("hexane_eluent", "hexane", "purification_agent", (E[8], I[8], E[5]), stage="purification"),
    ]
    # Mark the calculated loading as inferred rather than source-reported.
    pbu3 = next(p for p in r["participants"] if p["participant_id"] == "pbu3")
    pbu3["charges"]["equivalents"] = f(None, {"number": 0.05, "unit": "equiv"}, (E[8], I[8]), method="inference", derivation="0.01 mmol PBu3 / 0.20 mmol 1a = 0.05 equiv.")
    pbu3["charges"]["mol_percent"] = f(None, {"number": 5, "unit": "mol%"}, (E[8], I[8]), method="inference", derivation="0.01 mmol PBu3 / 0.20 mmol 1a × 100 = 5 mol%.")
    r["reaction_scale"] = f("1a (52 mg, 0.20 mmol)", {"reference_molecule_id": "1a", "amount": {"number": 0.20, "unit": "mmol"}}, (E[8], I[8]))
    r["solvent_conditions"]["composition"] = f("CH3CN", "CH3CN", (E[8], I[8]))
    r["solvent_conditions"]["total_volume"] = f("2 mL", {"number": 2, "unit": "mL"}, (E[8], I[8]))
    air = ("Unless otherwise specified, all reactions were performed under air atmosphere", "air", "inherited_general_condition", "Applied from the article's General Experimental Information; no exception is stated for this experiment.", (E[5],))
    r["steps"] = [
        new_step("reaction", (E[8], I[8], E[5]), ["substrate", "pbu3", "solvent", "product"], "treat 1a with PBu3 in CH3CN", temperature=("25 °C", {"number": 25, "unit": "°C"}), time=("40 min", {"number": 40, "unit": "min"}), atmosphere=air),
        new_step("purification", (E[8], I[8], E[5]), ["silica_column", "etoac_eluent", "hexane_eluent"], "followed by column chromatography; general experimental information specifies 100–200 mesh silica gel with an ethyl acetate–hexane gradient", stage="purification", operation_method="inherited_general_condition", operation_derivation="The specific paragraph states column chromatography; stationary phase and solvent system are inherited from the General Experimental Information. The reported Rf 3/7 ratio is not treated as the isolation eluent ratio."),
    ]
    proc = "Following the Kwons procedure, treatment of ethyl-2-(2-nitro-1-phenylethyl)buta-2,3-dienoate 1a (52 mg, 0.20 mmol) with PBu3 (2 mg, 0.01 mmol) in CH3CN (2 mL) at 25 °C for 40 min followed by column chromatography afforded Intermediate F as colorless oil (36 mg, 70%)."
    r["procedure"]["original_text"] = f(proc, proc, (E[8], I[8]))
    r["procedure"]["general_procedure_reference"] = f("Following the Kwons procedure, superscript 14", "Kwon procedure; source prints superscript 14", (E[8], I[8]))
    r["outcomes"] = [
        outcome("yield_intermediate_f", "yield", "intermediate_f", "70%", {"number": 70, "unit": "%"}, (E[8], I[8]), basis="isolated after column chromatography", yield_type="isolated", scope=("Intermediate F", {"start": "1a", "end": "intermediate_f", "single_step": True})),
        outcome("amount_intermediate_f", "product_amount", "intermediate_f", "36 mg", {"number": 36, "unit": "mg"}, (E[8], I[8]), basis="isolated product", scope=("Intermediate F", {"start": "1a", "end": "intermediate_f", "single_step": True})),
    ]
    changes = [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "N"}], "before_order": 1, "after_order": 2, "change_type": "order_increased"},
    ]
    derivation = "Partial net interpretation of the evidence-linked 1a and Intermediate F structures: furan-ring C–O closure and conversion of the pre-existing nitro C–N bond to the oxime C=N bond; no atom map is asserted."
    r["classification"] = {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": f(),
        "bond_changes": [derived(x, (I[4], E[8], I[8]), derivation) for x in changes],
        "labels": derived([
            {"category_number": 3, "code": "C-O", "supporting_change_ids": ["b1"]},
            {"category_number": 16, "code": "DB-C-N", "supporting_change_ids": ["b2"]},
        ], (I[4], E[8], I[8]), derivation),
        "coverage": inferred("partial", (I[4], E[8], I[8]), "Hydrogen transfers and the complete allene bond-order rearrangement are not atom-mapped."),
        "limitations": ["No atom-mapped reaction SMILES; hydrogen transfers and remaining bond-order changes are unresolved."],
        "model_votes": [],
    }
    r["notes"] = [
        f("Following the Kwons procedure, superscript 14", "The source paragraph visually prints citation 14, but reference 14 is Yao et al.; the Kwon furanone-oxime article appears as reference 16. The mismatch is retained without correction.", (E[8], I[8]), method="inference", derivation="Comparison of the page-8 citation marker with the page-9 reference list."),
        f("Rf (EtOAc/hexane: 3/7) = 0.15", "Rf solvent ratio is characterization context and is not assumed to be the preparative column eluent ratio.", (E[8], I[8]), method="inference", derivation="The paragraph separately states column chromatography without an eluent ratio."),
    ]
    r["review"] = {
        "status": "needs_review",
        "method": "model",
        "reviewer": None,
        "checked_fields": [],
        "unresolved": [
            "The printed Kwon-procedure citation is superscript 14, while the relevant Kwon article is reference 16; no correction was imposed.",
            "The C4 tetrahedral configuration is unreported; only the oxime Z geometry is specified.",
            "Workup and the specific isolation-column eluent ratio are unreported.",
            "Bond classification is partial because no atom mapping establishes all bond-order and hydrogen changes.",
        ],
    }
    r["completeness"] = "partial"
    bundle["reactions"].append(r)
    reactions[new_id] = r


# Non-synthetic material is retained as context rather than converted into
# spurious reactions.
context_updates = {
    f"{DOC}:page7_characterization": {
        "context_id": f"{DOC}:page7_characterization",
        "document_id": DOC,
        "evidence_ids": [E[7], I[7]],
        "kind": "characterization_blocks",
        "description": "Page 7 continues analytical characterization for 4c and gives characterization for 4d–4h and 5–8; these blocks support existing syntheses but are not independent reactions.",
    },
    f"{DOC}:page8_characterization": {
        "context_id": f"{DOC}:page8_characterization",
        "document_id": DOC,
        "evidence_ids": [E[8], I[8]],
        "kind": "characterization_blocks",
        "description": "Page 8 continues characterization for 8 and reports characterization for 9 and Intermediate F; peak lists and HRMS values are not independent syntheses.",
    },
    f"{DOC}:page8_si_notice_unassociated": {
        "context_id": f"{DOC}:page8_si_notice_unassociated",
        "document_id": DOC,
        "evidence_ids": [E[8], I[8]],
        "kind": "source_availability",
        "description": "The article states that Supporting Information and FAIR NMR files are available online, but no SI file is explicitly associated in this workspace.",
    },
    f"{DOC}:page8_9_bibliographic_matter": {
        "context_id": f"{DOC}:page8_9_bibliographic_matter",
        "document_id": DOC,
        "evidence_ids": [E[8], I[8], E[9], I[9]],
        "kind": "nonreaction_bibliographic_content",
        "description": "Author information, acknowledgments, accession information, and references are bibliographic/contextual material, not experiments in the present article.",
    },
}
existing_context = {x.get("context_id"): x for x in bundle["source_context"]}
existing_context.update(context_updates)
bundle["source_context"] = list(existing_context.values())


coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 7,
            "state": "attempted",
            "regions": [
                {"region": "left column, continuation of 4c and characterization of 4d–4h", "finding": "Analytical characterization only; linked to the existing Scheme 3 synthesis records and not counted as new reactions.", "evidence_ids": [E[7], I[7]]},
                {"region": "lower left to upper right, synthesis of 5", "finding": "Detailed 4b→5 bis-Boc protection; low-temperature slow addition, warm hold, workup, purification, charges, isolated mass and yield linked to the existing record.", "evidence_ids": [E[7], I[7], I[4]]},
                {"region": "upper-right middle, synthesis of 6", "finding": "Detailed DIBAL-H reduction 5→6 with separate −78 °C addition and 30 min hold, N2 atmosphere, workup/purification, 73% isolated yield and unassigned dr 5:1 linked to the existing record.", "evidence_ids": [E[7], I[7], I[4]]},
                {"region": "right column middle, synthesis of 7", "finding": "Detailed In(OTf)3-mediated 6→7 indolylation with workup/purification, 67% isolated yield and unassigned dr 3:1 linked to the existing record.", "evidence_ids": [E[7], I[7], I[4]]},
                {"region": "right column bottom, synthesis of 8", "finding": "Detailed 2a→8 reduction starts on this page; NaBH3CN addition, 5 min pre-acidification hold, 35% HCl addition and 5 h hold are represented as ordered steps. Characterization continues on page 8.", "evidence_ids": [E[7], I[7], E[8], I[8], I[4]]},
            ],
            "record_ids": scope_ids + [
                f"{DOC}:scheme4_4b_to_5",
                f"{DOC}:scheme4_5_to_6",
                f"{DOC}:scheme4_6_to_7",
                f"{DOC}:scheme4_2a_to_8",
            ],
            "unresolved": [
                "Absolute configurations of 4c–4h, 5, and 8 are not reported.",
                "Compounds 6 and 7 have reported dr values, but major/minor configurations and measurement methods are not reported.",
                "No reagent-addition duration is reported for the DIBAL-H addition at −78 °C.",
            ],
        },
        {
            "page": 8,
            "state": "attempted",
            "regions": [
                {"region": "upper-left continuation", "finding": "Analytical characterization of product 8; not a separate synthesis.", "evidence_ids": [E[8], I[8]]},
                {"region": "left column, synthesis and characterization of 9", "finding": "Detailed 2a→9 methoxylation with charges, 35% HCl, 25 °C/1 h, workup, purification, isolated mass and 70% yield linked to the existing record.", "evidence_ids": [E[8], I[8], I[4]]},
                {"region": "left column, Intermediate F paragraph", "finding": "Distinct performed preparation 1a→(Z)-Intermediate F using PBu3, CH3CN, 25 °C/40 min, column chromatography, 36 mg and 70% isolated yield; created as a mechanistic-study experiment.", "evidence_ids": [E[8], I[8], I[4]]},
                {"region": "left column lower half, analytical data for 9 and Intermediate F", "finding": "Rf/NMR/HRMS characterization retained as source context, not additional reactions.", "evidence_ids": [E[8], I[8]]},
                {"region": "associated content, accession codes, author information, acknowledgments and references 1–6", "finding": "Nonreaction material. SI availability is stated, but no SI file is associated with this workspace.", "evidence_ids": [E[8], I[8]]},
            ],
            "record_ids": [
                f"{DOC}:scheme4_2a_to_8",
                f"{DOC}:scheme4_2a_to_9",
                new_id,
            ],
            "unresolved": [
                "The paragraph prints a Kwon-procedure citation marker 14 although the matching Kwon paper is reference 16; retained as an internal source mismatch.",
                "Intermediate F has reported Z oxime geometry, but its tetrahedral C4 configuration is unreported.",
                "Intermediate F workup and specific column-eluent ratio are not reported; the Rf solvent ratio is not substituted for the isolation eluent.",
                "The article advertises SI, but no SI was explicitly provided or associated.",
            ],
        },
        {
            "page": 9,
            "state": "attempted",
            "regions": [
                {"region": "full page, references 7–18", "finding": "Reference list only. These cited publications are not experiments of the present article, and no reaction records are created from this page.", "evidence_ids": [E[9], I[9]]},
            ],
            "record_ids": [],
            "unresolved": [],
        },
    ],
}

inventory = {
    "document_id": DOC,
    "page_range": [7, 9],
    "chemical_accuracy_verified": False,
    "visual_pages_reviewed": [7, 8, 9],
    "new_reaction_ids": [new_id],
    "existing_reaction_ids_updated": [
        *scope_ids,
        f"{DOC}:scheme4_4b_to_5",
        f"{DOC}:scheme4_5_to_6",
        f"{DOC}:scheme4_6_to_7",
        f"{DOC}:scheme4_2a_to_8",
        f"{DOC}:scheme4_2a_to_9",
    ],
    "nonreaction_findings": [
        "Characterization peak lists on pages 7–8 were retained as source context.",
        "Associated-content notices do not establish a linked SI file.",
        "Page 9 is entirely references and produced no reaction records.",
    ],
}

stereo_diagnostics = {
    mid: st.stereo_inspect(st.value(molecules[mid]["smiles"]))
    for mid in ("4b", "5", "6", "7", "8", "9", "intermediate_f", "pbu3")
}

(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
(RUN / "coverage/pages_0007_0009.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/pages_0007_0009_inventory.json").write_text(json.dumps(inventory, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/stereo_diagnostics_pages_0007_0009.json").write_text(json.dumps(stereo_diagnostics, ensure_ascii=False, indent=2) + "\n")

print(json.dumps({
    "molecules": len(bundle["molecules"]),
    "reactions": len(bundle["reactions"]),
    "new_reactions": 1,
    "updated_existing_reactions": len(inventory["existing_reaction_ids_updated"]),
    "coverage_pages": len(coverage["pages"]),
}, ensure_ascii=False))
