#!/usr/bin/env python3
"""Augment the v3.1 submission with visually reviewed PDF pages 7--12."""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
TOOLS = ROOT / "skills" / "chemical-reaction-extraction" / "scripts"
sys.path.insert(0, str(TOOLS))
from standard_records import (  # noqa: E402
    fact,
    new_molecule,
    new_outcome,
    new_participant,
    new_reaction,
    new_step,
    stereo_inspect,
    write,
)

DOC = "doc_dfa009e2c8318d59"
RUN = ROOT / "runs" / "current"
SUBMITTED = RUN / "submitted.json"


def ev(*pages):
    out = []
    for page in pages:
        out.extend([f"{DOC}:page{page}:image", f"{DOC}:page{page}:text"])
    return out


def src(raw, normalized=None, evidence=(), method="original_extraction", *, status=None,
        verification=None, derivation=None, alternatives=()):
    return fact(raw=raw, normalized=raw if normalized is None else normalized,
                evidence=evidence, method=method, status=status,
                verification=verification, derivation=derivation,
                alternatives=alternatives)


def missing(evidence=(), status="unreported"):
    return fact(evidence=evidence, status=status)


def raw_only(raw, evidence=(), method="original_extraction", *, status=None, alternatives=()):
    return fact(raw=raw, normalized=None, evidence=evidence, method=method,
                status=status, alternatives=alternatives)


def na(evidence=()):
    return fact(evidence=evidence, method="not_applicable", status="not_applicable")


def qty(number, unit, raw_text, qualifier=None):
    v = {"number": number, "unit": unit, "raw_text": raw_text}
    if qualifier:
        v["qualifier"] = qualifier
    return v


def qrange(minimum, maximum, unit, raw_text, qualifier=None):
    v = {"minimum": minimum, "maximum": maximum, "unit": unit, "raw_text": raw_text}
    if qualifier:
        v["qualifier"] = qualifier
    return v


def participant(pid, mid, role, evidence, stage="reaction", product_role=None):
    p = new_participant(pid, mid, role, evidence, stage=stage)
    p["product_role"] = src(product_role, product_role, evidence) if product_role else na(evidence)
    return p


def charge(p, field, raw, normalized, evidence, method="original_extraction"):
    p["charges"][field] = src(raw, normalized, evidence, method=method)


def equivalence(p, raw, number, ref_mid, evidence):
    charge(p, "equivalents", raw, qty(number, "equiv", raw), evidence)
    p["equivalence_reference"] = src(
        ref_mid, ref_mid, evidence, method="inference",
        derivation="The source gives the equivalent next to this reagent; the named stoichiometric substrate is retained as the denominator.",
    )


def mol_percent(p, raw, number, ref_mid, evidence):
    charge(p, "mol_percent", raw, qty(number, "mol%", raw), evidence)
    p["equivalence_reference"] = src(
        ref_mid, ref_mid, evidence, method="inference",
        derivation="The catalyst loading is normalized against the charged 1-aza-1,3,4-triene amount; this denominator is inferred from the experimental paragraph.",
    )


def step(sid, evidence, pids, operation, temperature=None, duration=None,
         stage="reaction", atmosphere=None, inherited=False):
    s = new_step(sid, evidence, stage=stage)
    s["participant_ids"] = list(pids)
    method = "inherited_general_condition" if inherited else "original_extraction"
    s["operation"] = src(operation, operation, evidence, method=method)
    s["temperature"] = (raw_only(temperature[0], evidence, method=method)
                        if temperature and temperature[1] is None
                        else src(temperature[0], temperature[1], evidence, method=method)
                        if temperature else missing(evidence))
    s["time"] = (raw_only(duration[0], evidence, method=method)
                 if duration and duration[1] is None
                 else src(duration[0], duration[1], evidence, method=method)
                 if duration else missing(evidence))
    s["atmosphere"] = src(atmosphere, atmosphere, evidence, method=method) if atmosphere else missing(evidence)
    s["pressure"] = missing(evidence)
    return s


def outcome(oid, kind, targets, raw, normalized, evidence, *, product_role=None,
            yield_type=None, basis=None, scope=None):
    o = new_outcome(oid, kind, targets, evidence)
    o["product_role"] = src(product_role, product_role, evidence) if product_role else na(evidence)
    o["value"] = src(raw, normalized, evidence)
    o["yield_type"] = src(yield_type, yield_type, evidence) if yield_type else na(evidence)
    o["measurement_basis"] = src(basis, basis, evidence) if basis else missing(evidence)
    o["scope"] = src(scope, scope, evidence) if scope else missing(evidence)
    return o


def add_named_molecule(bundle, mid, label, name, evidence, smiles=None, *, chirality="achiral",
                       geometry="not_applicable", representation=None, note=None):
    if any(m["molecule_id"] == mid for m in bundle["molecules"]):
        return
    m = new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = src(label, label, evidence)
    m["names"]["original"] = src(name, name, evidence)
    m["names"]["english"] = src(name, name, evidence)
    m["names"]["abbreviations"] = missing(evidence)
    m["names"]["aliases"] = missing(evidence)
    m["cas_numbers"] = missing(evidence)
    m["generic_template"] = None
    st = m["stereochemistry"]
    if smiles:
        d = stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = src(
            {"source_label": label, "candidate_smiles": smiles}, d["canonical_isomeric_smiles"], evidence,
            method="structure_derivation",
            derivation="Candidate connectivity was derived from the source reagent name and normalized with RDKit; parsing is not source verification.",
        )
        st["atom_centers"] = src(label, d["atom_centers"], evidence, method="structure_derivation", derivation="RDKit graph diagnostic.")
        st["double_bonds"] = src(label, d["double_bonds"], evidence, method="structure_derivation", derivation="RDKit graph diagnostic.")
        st["axial"] = src(label, [], evidence, method="structure_derivation", derivation="No explicitly assigned axis is present in the candidate graph.")
        st["mixture"] = src("No stereoisomeric mixture is reported for this reagent.", "No stereoisomeric mixture is reported for this reagent.", evidence)
        st["assignment"] = src(label, "No configurational assignment is applicable to the candidate graph.", evidence, method="structure_derivation", derivation="Source name and graph review.")
        st["chirality_status"] = src(label, chirality, evidence, method="structure_derivation", derivation="Status is based on the candidate graph, not merely absence of stereobonds.")
        st["geometry_status"] = src(label, geometry, evidence, method="structure_derivation", derivation="No E/Z-capable stereogenic double bond is present unless separately stated.")
        has_unspecified = any(x["specified"] == "Unspecified" for x in d["potential_stereo"])
        rep = representation or ("partial" if has_unspecified else "encoded")
        st["representation_status"] = src(label, rep, evidence, method="structure_derivation", derivation="Candidate graph stereo diagnostic.")
        st["isomeric_smiles"] = src(smiles, d["canonical_isomeric_smiles"], evidence, method="structure_derivation", derivation="Canonical isomeric SMILES from the candidate graph.")
        st["identity_key"] = src(smiles, d["stereo_graph_key"], evidence, method="structure_derivation", derivation="RDKit stereo graph key; not a sample identity.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing(evidence)
        st["atom_centers"] = missing(evidence)
        st["double_bonds"] = missing(evidence)
        st["axial"] = missing(evidence)
        st["mixture"] = missing(evidence)
        st["assignment"] = src("No wedge/dash, R/S, Ra/Sa, E/Z, or cis/trans assignment is shown.", "unspecified", evidence)
        st["chirality_status"] = src(label, chirality, evidence, method="inference", derivation="Unspecified is not treated as achiral or racemic.")
        st["geometry_status"] = src(label, geometry, evidence, method="inference", derivation="Geometry status follows the explicit source description; no configuration is invented.")
        st["representation_status"] = src(label, representation or "not_encoded", evidence, method="inference", derivation="No complete reliable SMILES is asserted.")
        st["isomeric_smiles"] = missing(evidence)
        st["identity_key"] = missing(evidence)
    m["identity_note"] = src(note, note, evidence) if note else missing(evidence)
    m["external_properties"] = []
    bundle["molecules"].append(m)


def add_occurrence(bundle, mid, label, evidence):
    m = next(x for x in bundle["molecules"] if x["molecule_id"] == mid)
    if not any(o["evidence_ids"] == list(evidence) for o in m["occurrences"]):
        m["occurrences"].append({"document_id": DOC, "compound_label": src(label, label, evidence), "evidence_ids": list(evidence)})


def by_id(bundle, rid):
    return next(r for r in bundle["reactions"] if r["reaction_id"] == f"{DOC}:{rid}")


def set_common_reaction_fields(r, evidence, location, procedure_raw, procedure_norm, unresolved):
    r["evidence_ids"] = list(dict.fromkeys(r.get("evidence_ids", []) + list(evidence)))
    r["location"] = src(location, location, evidence)
    r["procedure"]["original_text"] = src(procedure_raw, procedure_norm, evidence)
    r["procedure"]["general_procedure_reference"] = src("Section 4.3 general procedure", "Section 4.3 general procedure", ev(7))
    r["procedure"]["modifications"] = missing(evidence)
    r["procedure"]["addition_order"] = src(
        "allene/THF; n-BuLi; hold; cool; isothiocyanate/THF; hold; cool; alkylating agent; warm; Cu(I); heat; quench/workup; chromatography",
        "allene/THF → n-BuLi → isothiocyanate/THF → alkylating agent → Cu(I) → workup → chromatography",
        evidence,
    )
    r["procedure"]["pretreatment"] = src("dry THF; anhydrous argon conditions", "dry THF under argon", ev(7))
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [], "unresolved": list(unresolved)}
    r["completeness"] = "partial"


def standard_workup_participants(evidence, alumina_layer="1.5 cm"):
    ps = []
    p = participant("p_nh4cl", "m_nh4cl", "workup_agent", evidence, stage="workup")
    charge(p, "volume", "30 mL saturated aqueous NH4Cl solution", qty(30, "mL", "30 mL"), evidence)
    ps.append(p)
    p = participant("p_nacn", "m_nacn", "workup_agent", evidence, stage="workup")
    charge(p, "mass", "0.50 g", qty(0.50, "g", "0.50 g"), evidence)
    ps.append(p)
    p = participant("p_ether_workup", "m_ether", "workup_agent", evidence, stage="workup")
    charge(p, "volume", "3 × 30 mL", {"number": 30, "unit": "mL", "raw_text": "3 × 30 mL"}, evidence)
    ps.append(p)
    p = participant("p_water_wash", "m_water", "workup_agent", evidence, stage="workup")
    charge(p, "volume", "50 mL", qty(50, "mL", "50 mL"), evidence)
    ps.append(p)
    ps.append(participant("p_mgso4", "m_mgso4", "workup_agent", evidence, stage="workup"))
    p = participant("p_alumina", "m_alumina", "purification_agent", evidence, stage="purification")
    p["charges"]["amount"] = src(alumina_layer, None, evidence)
    ps.append(p)
    ps.append(participant("p_hexane", "m_hexane", "solvent", evidence, stage="purification"))
    ps.append(participant("p_ether_eluent", "m_ether", "solvent", evidence, stage="purification"))
    return ps


def onepot_steps(evidence, ncs_pids, ncs_temp, ncs_time, cat_pid):
    inherited_ev = list(dict.fromkeys(ev(7) + list(evidence)))
    return [
        step("s1_nBuLi_addition", inherited_ev, ["p_allene", "p_thf_main", "p_nBuLi"],
             "add n-BuLi to the vigorously stirred allene solution",
             ("−100 to −90 °C", qrange(-100, -90, "°C", "−100 to −90 °C")),
             atmosphere="argon", inherited=True),
        step("s2_lithiation_hold", inherited_ev, ["p_allene", "p_nBuLi"],
             "stir after n-BuLi addition",
             ("ca. −65 to −40 °C", qrange(-65, -40, "°C", "ca. −65 to −40 °C", "approx")),
             ("5–15 min", qrange(5, 15, "min", "5–15 min")), atmosphere="argon", inherited=True),
        step("s3_cool_for_isothiocyanate", inherited_ev, ["p_allene", "p_nBuLi"],
             "cool before isothiocyanate addition",
             ("ca. −90 to −70 °C", qrange(-90, -70, "°C", "ca. −90 to −70 °C", "approx")),
             atmosphere="argon", inherited=True),
        step("s4_isothiocyanate_addition_hold", inherited_ev, ncs_pids,
             "add isothiocyanate/THF in one portion and stir",
             (ncs_temp[0], ncs_temp[1]), (ncs_time[0], ncs_time[1]), atmosphere="argon", inherited=True),
        step("s5_alkylator_addition", inherited_ev, ["p_alkylator"],
             "cool again and add the alkylating agent in one portion",
             ("ca. −80 °C", qty(-80, "°C", "ca. −80 °C", "approx")), atmosphere="argon", inherited=True),
        step("s6_warm_and_add_copper", inherited_ev, [cat_pid],
             "remove cooling bath; add Cu(I) catalyst when the mixture reaches −10 to −5 °C",
             ("−10 to −5 °C", qrange(-10, -5, "°C", "−10 to −5 °C")), atmosphere="argon", inherited=True),
        step("s7_cyclization", inherited_ev, ["p_allene", "p_ncs", "p_alkylator", cat_pid, "p_product"],
             "stir to effect intramolecular cyclization",
             ("32–50 °C", qrange(32, 50, "°C", "32–50 °C")),
             ("15–50 min", qrange(15, 50, "min", "15–50 min")), atmosphere="argon", inherited=True),
        step("s8_quench", inherited_ev, ["p_nh4cl", "p_nacn"],
             "cool to room temperature and treat with saturated aqueous NH4Cl containing NaCN",
             ("r.t.", "r.t."), stage="workup", inherited=True),
        step("s9_extract_wash_dry", inherited_ev, ["p_ether_workup", "p_water_wash", "p_mgso4"],
             "separate layers; extract aqueous phase with Et2O; combine organics, wash with water, dry over MgSO4, and concentrate",
             stage="workup", inherited=True),
        step("s10_purification", evidence, ["p_alumina", "p_hexane", "p_ether_eluent", "p_product"],
             "pass the ethereal solution through neutral alumina and purify by neutral-alumina column chromatography",
             stage="purification"),
    ]


bundle = json.loads(SUBMITTED.read_text())
bundle["documents"] = json.loads((RUN / "dataset.json").read_text())["documents"]

# Reagents/workup agents first seen explicitly in pages 7--12.
add_named_molecule(bundle, "m_cui", "CuI", "copper(I) iodide", ev(7), "[Cu+].[I-]")
add_named_molecule(bundle, "m_nh4cl", "NH4Cl", "ammonium chloride", ev(7), "[NH4+].[Cl-]")
add_named_molecule(bundle, "m_nacn", "NaCN", "sodium cyanide", ev(7), "[Na+].[C-]#N")
add_named_molecule(bundle, "m_ether", "Et2O", "diethyl ether", ev(7), "CCOCC")
add_named_molecule(bundle, "m_mgso4", "MgSO4", "magnesium sulfate", ev(7), "[Mg+2].[O-]S(=O)(=O)[O-]")
add_named_molecule(bundle, "m_alumina", "neutral Al2O3", "neutral alumina", ev(7), representation="not_encoded", chirality="not_applicable", geometry="not_applicable", note="Stationary phase; no molecular SMILES is asserted.")
add_named_molecule(bundle, "m_acetone", "acetone", "acetone", ev(11), "CC(C)=O")
add_named_molecule(bundle, "m_5j", "5j", "methyl {[(Z)-2-(ferrocenylmethoxy)-N-methylbuta-2,3-dienimidoyl]sulfanyl}acetate (1-aza-1,3,4-triene 5j)", ev(10), chirality="unspecified", geometry="mixture", representation="not_encoded", note="The heading names a Z form, while the isolated crude material is explicitly an E/Z mixture (~70:30 major/minor); no assignment of the major isomer is made.")
add_named_molecule(bundle, "m_5k", "5k", "tert-butyl {[(Z)-2-(ferrocenylmethoxy)-N-methylbuta-2,3-dienimidoyl]sulfanyl}acetate (1-aza-1,3,4-triene 5k)", ev(10, 11), chirality="unspecified", geometry="mixture", representation="not_encoded", note="The heading names a Z form, while the isolated crude material is explicitly an E/Z mixture (~64:36 major/minor); no assignment of the major isomer is made.")
add_named_molecule(bundle, "m_5l", "5l", "bromoacetonitrile-derived 1-aza-1,3,4-triene 5l", ev(11), chirality="unspecified", geometry="unresolved", representation="not_encoded", note="Only a small amount in a crude mixture is reported; structure and E/Z composition are not fully established in this passage.")

# Append page-specific occurrences without replacing the earlier table/scheme occurrences.
for key, pages in {"a": (7,), "b": (7,), "c": (8,), "d": (8,), "e": (8,), "f": (8, 9), "g": (9,), "h": (9,), "i": (10,), "j": (10,), "k": (11,), "l": (11,)}.items():
    add_occurrence(bundle, f"m_2{key}", f"2{key}", ev(*pages))
for mid, label, pages in [("m_7f", "7f", (9,)), ("m_12a", "12a", (11, 12)), ("m_12b", "12b", (12,)), ("m_17", "17", (12,))]:
    add_occurrence(bundle, mid, label, ev(*pages))

# Page 7 general procedure enriches the existing generic scheme, not a duplicate experiment.
generic = by_id(bundle, "scheme_table1_onepot")
generic["evidence_ids"] = list(dict.fromkeys(generic["evidence_ids"] + ev(7)))
generic["location"] = src("Pages 2–3 generic route and page 7 section 4.3 general procedure", "Pages 2–3 generic route; page 7 section 4.3", ev(2, 3, 7))
generic["procedure"]["original_text"] = src(
    {"section": "4.3", "allene": "10.0–15.0 mmol", "THF": "20–40 mL", "n-BuLi": "15.0–17.6 mmol", "isothiocyanate": "10.0–15.0 mmol", "alkylator": "MeI 1.5–2.4 equiv / propargyl bromide 1.2 equiv / benzyl bromide 1.0 equiv", "CuI": "2.10–3.15 mmol", "cyclization": "32–50 °C, 15–50 min"},
    "General lithiation → isothiocyanate addition → S-alkylation → Cu(I)-cyclization procedure with substrate-dependent low-temperature holds",
    ev(7),
)
generic["review"]["unresolved"] = ["The procedure gives ranges spanning multiple substrates; row records retain exact charges but inherit the applicable temperature/time branch."]

# Exact one-pot examples 2a--2i. Values are transcribed, not recalculated.
onepot = {
    "a": dict(pages=(7,), allene=("3.78 g", 3.78, 14.9, 40), nbu=(10, 1.6, 16.0), ncs=("m_me_ncs", "1.09 g", 1.09, 14.9, 3), alk=("m_mei", "3.25 g", 3.25, 22.9, 1.5), cat=("m_cui", "0.60 g", 0.60), ntemp=("−50 to −40 °C", qrange(-50, -40, "°C", "−50 to −40 °C")), ntime=("5–10 min", qrange(5, 10, "min", "5–10 min")), product=("2.81 g", 2.81, 55), eluents="hexane/Et2O 10:1, 3:1, 1:1"),
    "b": dict(pages=(7,), allene=("3.78 g", 3.78, 14.9, 40), nbu=(10, 1.6, 16.0), ncs=("m_ipr_ncs", "1.51 g", 1.51, 14.9, 5), alk=("m_mei", "3.75 g", 3.75, 26.4, 1.8), cat=("m_cui", "0.60 g", 0.60), ntemp=("ca. −30 °C", qty(-30, "°C", "ca. −30 °C", "approx")), ntime=("20 min", qty(20, "min", "20 min")), product=("2.41 g", 2.41, 44), eluents="hexane/Et2O 10:1, 3:1, 1:1"),
    "c": dict(pages=(8,), allene=("3.81 g", 3.81, 15.0, 25), nbu=(6, 2.5, 15.0), ncs=("m_voe_ncs", "1.94 g", 1.94, 15.0, 5), alk=("m_mei", "4.20 g", 4.20, 29.6, 2.0), cat=("m_cui", "0.60 g", 0.60), ntemp=("−45 to −40 °C", qrange(-45, -40, "°C", "−45 to −40 °C")), ntime=("10 min", qty(10, "min", "10 min")), product=("3.97 g", 3.97, 67), eluents="hexane/Et2O 10:1"),
    "d": dict(pages=(8,), allene=("3.78 g", 3.78, 14.9, 40), nbu=(11, 1.6, 17.6), ncs=("m_cpent_ncs", "1.89 g", 1.89, 14.9, 5), alk=("m_mei", "4.07 g", 4.07, 28.7, 1.9), cat=("m_cui", "0.60 g", 0.60), ntemp=("ca. −30 °C", qty(-30, "°C", "ca. −30 °C", "approx")), ntime=("20 min", qty(20, "min", "20 min")), product=("1.31 g", 1.31, 22), eluents="hexane/Et2O 10:1"),
    "e": dict(pages=(8,), allene=("3.78 g", 3.78, 14.9, 40), nbu=(11, 1.6, 17.6), ncs=("m_chept_ncs", "2.31 g", 2.31, 14.9, 5), alk=("m_mei", "4.37 g", 4.37, 30.8, 2.1), cat=("m_cui", "0.60 g", 0.60), ntemp=("ca. −30 °C", qty(-30, "°C", "ca. −30 °C", "approx")), ntime=("20 min", qty(20, "min", "20 min")), product=("2.32 g", 2.32, 37), eluents="hexane/Et2O 10:1"),
    "f": dict(pages=(8, 9), allene=("3.78 g", 3.78, 14.9, 40), nbu=(11, 1.6, 17.6), ncs=("m_ph_ncs", "2.00 g", 2.00, 14.8, 5), alk=("m_mei", "5.07 g", 5.07, 35.7, 2.4), cat=("m_cui", "0.60 g", 0.60), ntemp=("ca. −65 °C", qty(-65, "°C", "ca. −65 °C", "approx")), ntime=("15 min", qty(15, "min", "15 min")), product=("2.46 g", 2.46, 41), eluents="hexane/Et2O 10:1"),
    "g": dict(pages=(9,), allene=("3.81 g", 3.81, 15.0, 25), nbu=(6, 2.5, 15.0), ncs=("m_voe_ncs", "1.94 g", 1.94, 15.0, 5), alk=("m_eti", "4.79 g", 4.79, 30.7, 2.0), cat=("m_cubr", "0.60 g", 0.60), ntemp=("−45 to −40 °C", qrange(-45, -40, "°C", "−45 to −40 °C")), ntime=("10 min", qty(10, "min", "10 min")), product=("2.80 g", 2.80, 45), eluents="hexane/Et2O 10:1"),
    "h": dict(pages=(9,), allene=("3.78 g", 3.78, 14.9, 40), nbu=(6, 2.5, 15.0), ncs=("m_me_ncs", "1.09 g", 1.09, 14.9, 3), alk=("m_propargyl_br", "2.20 g", 2.20, 18.5, 1.2), cat=("m_cui", "0.60 g", 0.60), ntemp=("−50 to −40 °C", qrange(-50, -40, "°C", "−50 to −40 °C")), ntime=("5–10 min", qrange(5, 10, "min", "5–10 min")), product=("0.953 g", 0.953, 18), eluents="hexane/Et2O 10:1"),
    "i": dict(pages=(10,), allene=("2.54 g", 2.54, 10.0, 20), nbu=(4.5, 2.5, 11.2), ncs=("m_me_ncs", "0.73 g", 0.73, 10.0, 3), alk=("m_benzyl_br", "1.71 g", 1.71, 10.0, 1.0), cat=("m_cui", "0.40 g", 0.40), ntemp=("−50 to −40 °C", qrange(-50, -40, "°C", "−50 to −40 °C")), ntime=("5–10 min", qrange(5, 10, "min", "5–10 min")), product=("2.07 g", 2.07, 50), eluents="hexane/Et2O 10:1"),
}

for key, c in onepot.items():
    evidence = ev(*c["pages"])
    r = by_id(bundle, f"table1_2{key}")
    allene_mass_raw, allene_mass, scale_mmol, thf_main_ml = c["allene"]
    nbu_ml, nbu_m, nbu_mmol = c["nbu"]
    ncs_mid, ncs_mass_raw, ncs_mass, ncs_mmol, thf_ncs_ml = c["ncs"]
    alk_mid, alk_mass_raw, alk_mass, alk_mmol, alk_eq = c["alk"]
    cat_mid, cat_mass_raw, cat_mass = c["cat"]
    prod_mass_raw, prod_mass, percent = c["product"]
    p_allene = participant("p_allene", "m_allene1", "reactant", evidence)
    charge(p_allene, "mass", allene_mass_raw, qty(allene_mass, "g", allene_mass_raw), evidence)
    charge(p_allene, "amount", f"{scale_mmol} mmol", qty(scale_mmol, "mmol", f"{scale_mmol} mmol"), evidence)
    p_thf_main = participant("p_thf_main", "m_thf", "solvent", evidence)
    charge(p_thf_main, "volume", f"{thf_main_ml} mL", qty(thf_main_ml, "mL", f"{thf_main_ml} mL"), evidence)
    p_nbu = participant("p_nBuLi", "m_nBuLi", "reagent", evidence)
    charge(p_nbu, "volume", f"{nbu_ml} mL of {nbu_m} M solution in hexane", qty(nbu_ml, "mL", f"{nbu_ml} mL of {nbu_m} M solution in hexane"), evidence)
    charge(p_nbu, "amount", f"{nbu_mmol} mmol", qty(nbu_mmol, "mmol", f"{nbu_mmol} mmol"), evidence)
    p_ncs = participant("p_ncs", ncs_mid, "reactant", evidence)
    charge(p_ncs, "mass", ncs_mass_raw, qty(ncs_mass, "g", ncs_mass_raw), evidence)
    charge(p_ncs, "amount", f"{ncs_mmol} mmol", qty(ncs_mmol, "mmol", f"{ncs_mmol} mmol"), evidence)
    p_thf_ncs = participant("p_thf_ncs", "m_thf", "solvent", evidence)
    charge(p_thf_ncs, "volume", f"{thf_ncs_ml} mL", qty(thf_ncs_ml, "mL", f"{thf_ncs_ml} mL"), evidence)
    p_alk = participant("p_alkylator", alk_mid, "reagent", evidence)
    charge(p_alk, "mass", alk_mass_raw, qty(alk_mass, "g", alk_mass_raw), evidence)
    charge(p_alk, "amount", f"{alk_mmol} mmol", qty(alk_mmol, "mmol", f"{alk_mmol} mmol"), evidence)
    equivalence(p_alk, f"{alk_eq} equiv", alk_eq, ncs_mid, evidence)
    p_cat = participant("p_copper", cat_mid, "catalyst", evidence)
    charge(p_cat, "mass", cat_mass_raw, qty(cat_mass, "g", cat_mass_raw), evidence)
    p_prod = participant("p_product", f"m_2{key}", "product", evidence, product_role="main product")
    charge(p_prod, "mass", prod_mass_raw, qty(prod_mass, "g", prod_mass_raw), evidence)
    r["participants"] = [p_allene, p_thf_main, p_nbu, p_ncs, p_thf_ncs, p_alk, p_cat, p_prod] + standard_workup_participants(evidence)
    r["reaction_scale"] = src(f"{ncs_mmol} mmol isothiocyanate", {"reference_molecule_id": ncs_mid, "amount": qty(ncs_mmol, "mmol", f"{ncs_mmol} mmol")}, evidence)
    r["solvent_conditions"]["composition"] = src("dry THF; n-BuLi solution in hexane", ["THF", "hexane"], evidence)
    r["solvent_conditions"]["mixture_ratio"] = missing(evidence)
    r["solvent_conditions"]["total_volume"] = raw_only(f"THF {thf_main_ml} mL + {thf_ncs_ml} mL", evidence)
    r["solvent_conditions"]["concentration"] = missing(evidence)
    r["solvent_conditions"]["concentration_basis"] = missing(evidence)
    r["steps"] = onepot_steps(evidence, ["p_ncs", "p_thf_ncs"], c["ntemp"], c["ntime"], "p_copper")
    r["steps"][-1]["operation"] = src(f"pass through a 1.5 cm neutral-alumina layer; column chromatography on neutral alumina using hexane then {c['eluents']}", f"neutral-alumina cleanup and chromatography ({c['eluents']})", evidence)
    r["outcomes"] = [
        outcome("o_product_mass", "product_amount", [f"m_2{key}"], prod_mass_raw, qty(prod_mass, "g", prod_mass_raw), evidence, product_role="main product", basis="isolated product"),
        outcome("o_yield", "yield", [f"m_2{key}"], f"{percent}%", qty(percent, "%", f"{percent}%"), evidence, product_role="main product", yield_type="isolated yield", basis="isolated product; yield based on isothiocyanate", scope={"starting_material": "isothiocyanate", "product": f"2{key}"}),
    ]
    if key == "g":
        p7g = participant("p_trace_7g", "m_7g", "product", evidence, product_role="rearrangement contaminant")
        r["participants"].append(p7g)
        r["outcomes"].append(outcome("o_7g_contamination", "trace", ["m_7g"], "~6%", qty(6, "%", "~6%", "approx"), evidence, product_role="rearrangement contaminant", basis="composition of isolated 2g; analytical method not reported", scope={"matrix": "isolated 2g", "component": "7g"}))
    set_common_reaction_fields(r, evidence, f"Section 4.3.{1 + ord(key) - ord('a')}, pages {c['pages'][0]}–{c['pages'][-1]}" if len(c["pages"]) > 1 else f"Section 4.3.{1 + ord(key) - ord('a')}, page {c['pages'][0]}",
                               {"section": f"4.3.{1 + ord(key) - ord('a')}", "product": f"2{key}", "exact_charges": True, "general_procedure": "4.3"},
                               f"Exact preparation of 2{key} under the section 4.3 general procedure",
                               ["The row-specific cyclization temperature/time is reported only as the section 4.3 range (32–50 °C, 15–50 min).", "Complete organometallic product SMILES and atom mapping remain unresolved."] + (["The isolated 2g contains approximately 6% 7g; the analytical basis for this composition is not stated."] if key == "g" else []))


def two_step_record(key, pages, scale, thf_main, nbu, ncs, water, alk_mid, alk, stage1_temp,
                    stage1_crude, stage1_ratio, isomer_ratio, precursor_charge, cubr,
                    cyclization, crude2, isolated, yield_percent, eluent, extra_unresolved):
    evidence = ev(*pages)
    r = by_id(bundle, f"table1_2{key}")
    interm = f"m_5{key}"
    ps = []
    p = participant("p_allene", "m_allene1", "reactant", evidence); charge(p, "mass", scale[0], qty(scale[1], "g", scale[0]), evidence); charge(p, "amount", scale[2], qty(scale[3], "mmol", scale[2]), evidence); ps.append(p)
    p = participant("p_thf_main", "m_thf", "solvent", evidence); charge(p, "volume", thf_main, qty(float(thf_main.split()[0]), "mL", thf_main), evidence); ps.append(p)
    p = participant("p_nBuLi", "m_nBuLi", "reagent", evidence); charge(p, "volume", nbu[0], qty(nbu[1], "mL", nbu[0]), evidence); charge(p, "amount", nbu[2], qty(nbu[3], "mmol", nbu[2]), evidence); ps.append(p)
    p = participant("p_ncs", "m_me_ncs", "reactant", evidence); charge(p, "mass", ncs[0], qty(ncs[1], "g", ncs[0]), evidence); charge(p, "amount", ncs[2], qty(ncs[3], "mmol", ncs[2]), evidence); ps.append(p)
    p = participant("p_thf_ncs", "m_thf", "solvent", evidence); charge(p, "volume", "3 mL", qty(3, "mL", "3 mL"), evidence); ps.append(p)
    p = participant("p_water", "m_water", "additive", evidence); charge(p, "mass", water[0], qty(water[1], "g", water[0]), evidence); charge(p, "amount", water[2], qty(water[3], "mmol", water[2]), evidence); ps.append(p)
    p = participant("p_thf_water", "m_thf", "solvent", evidence); charge(p, "volume", "2 mL", qty(2, "mL", "2 mL"), evidence); ps.append(p)
    p = participant("p_alkylator", alk_mid, "reagent", evidence); charge(p, "mass", alk[0], qty(alk[1], "g", alk[0]), evidence); charge(p, "amount", alk[2], qty(alk[3], "mmol", alk[2]), evidence); equivalence(p, alk[4], alk[5], "m_me_ncs", evidence); ps.append(p)
    ps.append(participant("p_intermediate", interm, "product", evidence, product_role="crude isolated intermediate"))
    p = participant("p_intermediate_charge", interm, "reactant", evidence); charge(p, "mass", precursor_charge[0], qty(precursor_charge[1], "g", precursor_charge[0]), evidence); charge(p, "amount", precursor_charge[2], qty(precursor_charge[3], "mmol", precursor_charge[2]), evidence); ps.append(p)
    p = participant("p_thf_cyclization", "m_thf", "solvent", evidence); charge(p, "volume", "20 mL", qty(20, "mL", "20 mL"), evidence); ps.append(p)
    p = participant("p_cubr", "m_cubr", "catalyst", evidence); charge(p, "mass", cubr[0], qty(cubr[1], "g", cubr[0]), evidence); mol_percent(p, cubr[2], cubr[3], interm, evidence); ps.append(p)
    p = participant("p_product", f"m_2{key}", "product", evidence, product_role="main product"); charge(p, "mass", isolated[0], qty(isolated[1], "g", isolated[0]), evidence); ps.append(p)
    # Explicit stage-1 quench and stage-2 workup doses.
    for suffix in ("stage1", "stage2"):
        p = participant(f"p_nh4cl_{suffix}", "m_nh4cl", "workup_agent", evidence, stage="workup"); charge(p, "volume", "10 mL saturated aqueous NH4Cl solution", qty(10, "mL", "10 mL"), evidence); ps.append(p)
    p = participant("p_nacn_stage2", "m_nacn", "workup_agent", evidence, stage="workup"); charge(p, "mass", "0.20 g", qty(0.20, "g", "0.20 g"), evidence); ps.append(p)
    p = participant("p_ether_extract", "m_ether", "workup_agent", evidence, stage="workup"); charge(p, "volume", "2 × 15 mL", {"number": 15, "unit": "mL", "raw_text": "2 × 15 mL"}, evidence); ps.append(p)
    p = participant("p_water_wash", "m_water", "workup_agent", evidence, stage="workup"); charge(p, "volume", "2 × 20 mL", {"number": 20, "unit": "mL", "raw_text": "2 × 20 mL"}, evidence); ps.append(p)
    ps += [participant("p_mgso4", "m_mgso4", "workup_agent", evidence, stage="workup"), participant("p_alumina", "m_alumina", "purification_agent", evidence, stage="purification"), participant("p_hexane", "m_hexane", "solvent", evidence, stage="purification"), participant("p_ether_eluent", "m_ether", "solvent", evidence, stage="purification")]
    r["participants"] = ps
    r["reaction_scale"] = src(scale[2] + " allene 1", {"reference_molecule_id": "m_allene1", "amount": qty(scale[3], "mmol", scale[2])}, evidence)
    r["solvent_conditions"]["composition"] = src("dry THF; n-BuLi solution in hexane", ["THF", "hexane"], evidence)
    r["solvent_conditions"]["mixture_ratio"] = missing(evidence)
    r["solvent_conditions"]["total_volume"] = raw_only(f"THF {thf_main} + 3 mL + 2 mL in stage 1; THF 20 mL in stage 2", evidence)
    r["solvent_conditions"]["concentration"] = missing(evidence)
    r["solvent_conditions"]["concentration_basis"] = missing(evidence)
    r["steps"] = [
        step("s1_nBuLi_addition", evidence, ["p_allene", "p_thf_main", "p_nBuLi"], "add n-BuLi to allene 1 in dry THF", ("−90 °C", qty(-90, "°C", "−90 °C")), atmosphere="argon"),
        step("s2_lithiation_hold", evidence, ["p_allene", "p_nBuLi"], "allow temperature to rise and stir", ("−45 to −40 °C", qrange(-45, -40, "°C", "−45 to −40 °C")), ("10 min", qty(10, "min", "10 min")), atmosphere="argon"),
        step("s3_cool_add_ncs", evidence, ["p_ncs", "p_thf_ncs"], "cool to −90 °C, add MeNCS/THF in one portion, and stir", ("−48 to −40 °C", qrange(-48, -40, "°C", "−48 to −40 °C")), ("10 min", qty(10, "min", "10 min")), atmosphere="argon"),
        step("s4_water_then_alkylator", evidence, ["p_water", "p_thf_water", "p_alkylator"], f"cool to −95 °C; add water/THF, then add alkylator at {stage1_temp}", (f"−95 °C then {stage1_temp}", None), atmosphere="argon"),
        step("s5_warm_hold", evidence, ["p_alkylator"], "remove cooling bath and allow temperature to rise", ("−30 °C", qty(-30, "°C", "−30 °C")), ("10 min", qty(10, "min", "10 min")), atmosphere="argon"),
        step("s6_cool_quench_stage1", evidence, ["p_nh4cl_stage1", "p_intermediate"], "cool again and quench with saturated aqueous NH4Cl; perform usual workup", ("−80 °C", qty(-80, "°C", "−80 °C")), stage="workup"),
        step("s7_add_cubr", evidence, ["p_intermediate_charge", "p_thf_cyclization", "p_cubr"], "add CuBr to the stirred crude 1-aza-1,3,4-triene solution at room temperature", ("r.t.; exotherm to 43 °C", None)),
        step("s8_cyclization_hold", evidence, ["p_intermediate_charge", "p_cubr", "p_product"], "stir for an additional 30 min at room temperature", ("r.t.", "r.t."), ("30 min", qty(30, "min", "30 min"))),
        step("s9_stage2_workup", evidence, ["p_nh4cl_stage2", "p_nacn_stage2", "p_ether_extract", "p_water_wash", "p_mgso4"], "treat with saturated NH4Cl containing NaCN; separate, extract with Et2O, wash with water, dry over MgSO4", stage="workup"),
        step("s10_stage2_purification", evidence, ["p_alumina", "p_hexane", "p_ether_eluent", "p_product"], f"pass through a 2 cm neutral-alumina layer, concentrate, and purify by neutral-alumina chromatography ({eluent})", stage="purification"),
    ]
    r["outcomes"] = [
        outcome("o_stage1_crude_mass", "product_amount", [interm, "m_allene1"], stage1_crude[0], qty(stage1_crude[1], "g", stage1_crude[0]), evidence, product_role="crude intermediate mixture", basis="total crude residue"),
        outcome("o_stage1_composition", "other", [interm, "m_allene1"], stage1_ratio, stage1_ratio, evidence, product_role="crude intermediate mixture", basis="1H NMR", scope={"components": [f"5{key}", "unreacted allene 1"]}),
        outcome("o_stage1_geometry", "selectivity", [interm], isomer_ratio, isomer_ratio, evidence, product_role="E/Z intermediate mixture", basis="1H NMR", scope={"metric": "major:minor", "configuration_order": "unknown; source does not assign major E or Z"}),
    ]
    if crude2:
        r["outcomes"].append(outcome("o_stage2_crude", "other", [f"m_2{key}"], crude2[0], qty(crude2[1], "g", crude2[0]), evidence, product_role="crude product", basis=crude2[2], scope={"stage": "CuBr cyclization before column chromatography"}))
    r["outcomes"] += [
        outcome("o_product_mass", "product_amount", [f"m_2{key}"], isolated[0], qty(isolated[1], "g", isolated[0]), evidence, product_role="main product", basis="isolated product"),
        outcome("o_yield", "yield", [f"m_2{key}"], f"{yield_percent}%", qty(yield_percent, "%", f"{yield_percent}%"), evidence, product_role="main product", yield_type="isolated yield", basis="experimental paragraph reports the percentage after charging the corrected crude intermediate amount", scope={"starting_material": f"crude 5{key} charge", "product": f"2{key}", "table_scope": "Table 1 presents the row as a two-preparative-step synthesis based on isothiocyanate"}),
    ]
    set_common_reaction_fields(r, evidence, f"Section 4.3.10{'–4.3.11' if key == 'k' else ''}, pages {pages[0]}–{pages[-1]}" if pages[0] != pages[-1] else f"Section 4.3.{10 if key == 'j' else 11}, page {pages[0]}",
                               {"section": f"4.3.{10 if key == 'j' else 11}", "stages": [f"allene 1 to crude 5{key}", f"crude 5{key} to 2{key}"]},
                               f"Two-preparative-step synthesis of 2{key} via crude 1-aza-1,3,4-triene 5{key}",
                               ["The heading assigns (Z), but the isolated intermediate is explicitly an E/Z mixture and the major isomer is not assigned.", "Table 1 states yields are based on isothiocyanate, whereas the detailed cyclization paragraph reports the same final percentage after charging a corrected mmol amount of crude intermediate; both scopes are retained without choosing one."] + list(extra_unresolved))


two_step_record("j", (10,), ("3.81 g", 3.81, "15.0 mmol", 15.0), "25 mL", ("7.5 mL of 2.37 M solution in hexane", 7.5, "17.8 mmol", 17.8), ("1.10 g", 1.10, "15.0 mmol", 15.0), ("0.28 g", 0.28, "15.5 mmol", 15.5), "m_methyl_bromoacetate", ("2.36 g", 2.36, "15.4 mmol", 15.4, "1.0 equiv", 1.0), "−80 °C", ("6.17 g", 6.17), "~73:27", "~70:30", ("6.02 g", 6.02, "11.0 mmol (accounting for 27% allene impurity)", 11.0), ("0.15 g", 0.15, "8.7 mol%", 8.7), None, None, ("1.37 g", 1.37), 31, "hexane/Et2O 10:1", [])
two_step_record("k", (10, 11), ("3.81 g", 3.81, "15.0 mmol", 15.0), "25 mL", ("7.5 mL of 2.37 M solution in hexane", 7.5, "17.8 mmol", 17.8), ("1.10 g", 1.10, "15.0 mmol", 15.0), ("0.28 g", 0.28, "15.5 mmol", 15.5), "m_tbutyl_bromoacetate", ("3.00 g", 3.00, "15.3 mmol", 15.3, "1.0 equiv", 1.0), "−80 °C", ("6.86 g", 6.86), "~70:30", "~64:36", ("6.80 g", 6.80, "10.8 mmol (accounting for 30% allene impurity)", 10.8), ("0.15 g", 0.15, "8.8 mol%", 8.8), None, ("5.90 g (87%)", 5.90, "source reports 87% for crude residue consisting of 2k"), ("0.96 g", 0.96), 20, "hexane/Et2O 10:1 then 3:1", ["The 87% value is attached to the crude residue before column chromatography and is not substituted for the 20% isolated yield."])

# Two-stage low-yield bromoacetonitrile example 2l.
evidence = ev(11)
r = by_id(bundle, "table1_2l")
ps = []
for pid, mid, role, stage in [("p_allene", "m_allene1", "reactant", "reaction"), ("p_thf_main", "m_thf", "solvent", "reaction"), ("p_nBuLi", "m_nBuLi", "reagent", "reaction"), ("p_ncs", "m_me_ncs", "reactant", "reaction"), ("p_thf_ncs", "m_thf", "solvent", "reaction"), ("p_water", "m_water", "additive", "reaction"), ("p_thf_water", "m_thf", "solvent", "reaction"), ("p_alkylator", "m_bromoacetonitrile", "reagent", "reaction"), ("p_intermediate", "m_5l", "product", "reaction"), ("p_cubr", "m_cubr", "catalyst", "reaction"), ("p_thf_cyclization", "m_thf", "solvent", "reaction"), ("p_product", "m_2l", "product", "reaction")]:
    ps.append(participant(pid, mid, role, evidence, stage=stage, product_role="main product" if pid == "p_product" else "crude intermediate" if pid == "p_intermediate" else None))
charges = {
    "p_allene": [("mass", "2.54 g", qty(2.54, "g", "2.54 g")), ("amount", "10.0 mmol", qty(10.0, "mmol", "10.0 mmol"))],
    "p_thf_main": [("volume", "20 mL", qty(20, "mL", "20 mL"))],
    "p_nBuLi": [("volume", "4.5 mL of 2.5 M solution in hexane", qty(4.5, "mL", "4.5 mL of 2.5 M solution in hexane")), ("amount", "11 mmol", qty(11, "mmol", "11 mmol"))],
    "p_ncs": [("mass", "0.73 g", qty(0.73, "g", "0.73 g")), ("amount", "10.0 mmol", qty(10.0, "mmol", "10.0 mmol"))],
    "p_thf_ncs": [("volume", "3 mL", qty(3, "mL", "3 mL"))],
    "p_water": [("mass", "0.20 g", qty(0.20, "g", "0.20 g")), ("amount", "11 mmol", qty(11, "mmol", "11 mmol"))],
    "p_thf_water": [("volume", "2 mL", qty(2, "mL", "2 mL"))],
    "p_alkylator": [("mass", "1.36 g", qty(1.36, "g", "1.36 g")), ("amount", "11 mmol", qty(11, "mmol", "11 mmol")), ("equivalents", "1.1 equiv", qty(1.1, "equiv", "1.1 equiv"))],
    "p_cubr": [("mass", "0.40 g", qty(0.40, "g", "0.40 g"))],
    "p_thf_cyclization": [("volume", "20 mL", qty(20, "mL", "20 mL"))],
    "p_product": [("mass", "0.11 g", qty(0.11, "g", "0.11 g"))],
}
for p in ps:
    for fld, raw, norm in charges.get(p["participant_id"], []):
        charge(p, fld, raw, norm, evidence)
    if p["participant_id"] == "p_alkylator":
        p["equivalence_reference"] = src("m_me_ncs", "m_me_ncs", evidence, method="inference", derivation="The stated 1.1 equiv is normalized against the stoichiometric MeNCS charge.")
ps += [participant("p_nh4cl", "m_nh4cl", "workup_agent", evidence, stage="workup"), participant("p_alumina", "m_alumina", "purification_agent", evidence, stage="purification"), participant("p_hexane", "m_hexane", "solvent", evidence, stage="purification"), participant("p_ether", "m_ether", "solvent", evidence, stage="purification"), participant("p_acetone", "m_acetone", "solvent", evidence, stage="purification")]
charge(ps[-5], "volume", "10 mL saturated aqueous NH4Cl solution", qty(10, "mL", "10 mL"), evidence)
r["participants"] = ps
r["reaction_scale"] = src("10.0 mmol allene 1", {"reference_molecule_id": "m_allene1", "amount": qty(10.0, "mmol", "10.0 mmol")}, evidence)
r["solvent_conditions"]["composition"] = src("dry THF; n-BuLi solution in hexane", ["THF", "hexane"], evidence)
r["solvent_conditions"]["mixture_ratio"] = missing(evidence); r["solvent_conditions"]["total_volume"] = raw_only("THF 20 + 3 + 2 mL in stage 1; THF 20 mL in stage 2", evidence); r["solvent_conditions"]["concentration"] = missing(evidence); r["solvent_conditions"]["concentration_basis"] = missing(evidence)
r["steps"] = [
    step("s1_nBuLi_addition", evidence, ["p_allene", "p_thf_main", "p_nBuLi"], "add n-BuLi to allene 1 in dry THF", ("−90 °C", qty(-90, "°C", "−90 °C")), atmosphere="argon"),
    step("s2_lithiation_hold", evidence, ["p_allene", "p_nBuLi"], "stir after addition", ("−46 to −38 °C", qrange(-46, -38, "°C", "−46 to −38 °C")), ("13 min", qty(13, "min", "13 min")), atmosphere="argon"),
    step("s3_ncs_addition_hold", evidence, ["p_ncs", "p_thf_ncs"], "cool to −80 °C, add MeNCS/THF in one portion, then stir", ("−43 to −40 °C", qrange(-43, -40, "°C", "−43 to −40 °C")), ("10 min", qty(10, "min", "10 min")), atmosphere="argon"),
    step("s4_water_bromoacetonitrile", evidence, ["p_water", "p_thf_water", "p_alkylator"], "cool to −95 °C; add water/THF, then add BrCH2CN at −83 °C", ("−95 °C then −83 °C", None), atmosphere="argon"),
    step("s5_warm_recool_quench", evidence, ["p_nh4cl", "p_intermediate"], "raise to −40 °C, cool again to −60 °C, and quench with saturated aqueous NH4Cl; usual workup", ("−40 °C then −60 °C", None), stage="workup"),
    step("s6_cubr_cyclization", evidence, ["p_intermediate", "p_thf_cyclization", "p_cubr", "p_product"], "add CuBr to the crude product in THF and stir", ("40–50 °C", qrange(40, 50, "°C", "40–50 °C")), ("40 min", qty(40, "min", "40 min"))),
    step("s7_cleanup", evidence, ["p_alumina", "p_product"], "perform usual workup, pass through a 1.5 cm neutral-alumina layer, and concentrate", stage="purification"),
    step("s8_final_purification", evidence, ["p_alumina", "p_hexane", "p_ether", "p_acetone", "p_product"], "column chromatography on neutral alumina (hexane/Et2O 1:1), then recrystallize from acetone/Et2O", stage="purification"),
]
r["outcomes"] = [
    outcome("o_stage1_crude", "other", ["m_5l", "m_allene1"], "4.19 g crude residue; small amount 5l, unreacted allene 1, and unidentified tar-like products", {"mass": qty(4.19, "g", "4.19 g"), "composition": "small amount 5l + unreacted allene 1 + unidentified tar-like products"}, evidence, product_role="crude intermediate mixture", basis="1H NMR spectroscopic analysis"),
    outcome("o_stage2_composition", "other", ["m_2l", "m_allene1"], "1.72 g residue containing ~22% 2l and ~78% allene 1", {"mass": qty(1.72, "g", "1.72 g"), "composition": "2l:allene 1 ~22:78"}, evidence, product_role="crude product mixture", basis="source-reported composition; method not explicitly restated"),
    outcome("o_product_mass", "product_amount", ["m_2l"], "0.11 g", qty(0.11, "g", "0.11 g"), evidence, product_role="main product", basis="isolated after chromatography and recrystallization"),
    outcome("o_yield", "yield", ["m_2l"], "3%", qty(3, "%", "3%"), evidence, product_role="main product", yield_type="isolated yield", basis="isolated product", scope={"starting_material": "allene 1/MeNCS sequence", "product": "2l"}),
]
set_common_reaction_fields(r, evidence, "Section 4.3.12, page 11", {"section": "4.3.12", "stages": ["allene 1 to crude 5l", "crude 5l to 2l"]}, "Two-stage bromoacetonitrile route to 2l", ["The amount and E/Z composition of 5l are not quantified.", "The method used for the reported ~22:78 stage-2 residue composition is not explicitly restated.", "Unidentified tar-like products remain unresolved."])


def comparator_record(key, pages, start_mid, start_mass, start_mmol, thf_main, nbu_vol, nbu_m, nbu_mmol, nbu_add_temp, lith_temp, lith_time, phncs_mass, phncs_mmol, thf_ncs, ncs_temp, ncs_time, mei_mass, mei_mmol, cui_mass, cu_add_temp, heat_temp, heat_time, nh4cl_vol, nacn_mass, extraction, wash, residue_mass, product_mass, percent, eluent):
    evidence = ev(*pages)
    r = by_id(bundle, f"scheme3_synthesis_12{key}")
    specs = [
        ("p_start", start_mid, "reactant", "reaction", [("mass", start_mass, qty(float(start_mass.split()[0]), "g", start_mass)), ("amount", start_mmol, qty(float(start_mmol.split()[0]), "mmol", start_mmol))]),
        ("p_thf_main", "m_thf", "solvent", "reaction", [("volume", thf_main, qty(float(thf_main.split()[0]), "mL", thf_main))]),
        ("p_nBuLi", "m_nBuLi", "reagent", "reaction", [("volume", f"{nbu_vol} of {nbu_m} M solution in hexane", qty(float(nbu_vol.split()[0]), "mL", f"{nbu_vol} of {nbu_m} M solution in hexane")), ("amount", nbu_mmol, qty(float(nbu_mmol.split()[0]), "mmol", nbu_mmol))]),
        ("p_phncs", "m_ph_ncs", "reactant", "reaction", [("mass", phncs_mass, qty(float(phncs_mass.split()[0]), "g", phncs_mass)), ("amount", phncs_mmol, qty(float(phncs_mmol.split()[0]), "mmol", phncs_mmol))]),
        ("p_thf_ncs", "m_thf", "solvent", "reaction", [("volume", thf_ncs, qty(float(thf_ncs.split()[0]), "mL", thf_ncs))]),
        ("p_mei", "m_mei", "reagent", "reaction", [("mass", mei_mass, qty(float(mei_mass.split()[0]), "g", mei_mass)), ("amount", mei_mmol, qty(float(mei_mmol.split()[0]), "mmol", mei_mmol)), ("equivalents", "2.0 equiv" if key == "b" else "2 equiv", qty(2.0, "equiv", "2.0 equiv" if key == "b" else "2 equiv"))]),
        ("p_cui", "m_cui", "catalyst", "reaction", [("mass", cui_mass, qty(float(cui_mass.split()[0]), "g", cui_mass))]),
        ("p_product", f"m_12{key}", "product", "reaction", [("mass", product_mass, qty(float(product_mass.split()[0]), "g", product_mass))]),
        ("p_nh4cl", "m_nh4cl", "workup_agent", "workup", [("volume", nh4cl_vol, qty(float(nh4cl_vol.split()[0]), "mL", nh4cl_vol))]),
        ("p_nacn", "m_nacn", "workup_agent", "workup", [("mass", nacn_mass, qty(float(nacn_mass.replace('~','').split()[0]), "g", nacn_mass, "approx" if "~" in nacn_mass else None))]),
        ("p_ether", "m_ether", "workup_agent", "workup", [("volume", extraction, {"raw_text": extraction, "unit": "mL", "number": float(extraction.split('×')[-1].split()[0])})]),
        ("p_wash", "m_nh4cl" if "NH4Cl" in wash else "m_water", "workup_agent", "workup", [("volume", wash, qty(float(wash.split()[0]), "mL", wash))]),
        ("p_mgso4", "m_mgso4", "workup_agent", "workup", []),
        ("p_alumina", "m_alumina", "purification_agent", "purification", []),
        ("p_hexane", "m_hexane", "solvent", "purification", []),
        ("p_ether_eluent", "m_ether", "solvent", "purification", []),
    ]
    ps = []
    for pid, mid, role, stage, cs in specs:
        p = participant(pid, mid, role, evidence, stage=stage, product_role="main product" if pid == "p_product" else None)
        for fld, raw, norm in cs: charge(p, fld, raw, norm, evidence)
        if pid == "p_mei":
            p["equivalence_reference"] = src("m_ph_ncs", "m_ph_ncs", evidence, method="inference", derivation="The stated 2 equiv is normalized against the stoichiometric PhNCS charge.")
        ps.append(p)
    r["participants"] = ps
    r["reaction_scale"] = src(phncs_mmol + " PhNCS", {"reference_molecule_id": "m_ph_ncs", "amount": qty(float(phncs_mmol.split()[0]), "mmol", phncs_mmol)}, evidence)
    r["solvent_conditions"]["composition"] = src("dry THF; n-BuLi solution in hexane", ["THF", "hexane"], evidence)
    r["solvent_conditions"]["mixture_ratio"] = missing(evidence); r["solvent_conditions"]["total_volume"] = raw_only(f"THF {thf_main} + {thf_ncs}", evidence); r["solvent_conditions"]["concentration"] = missing(evidence); r["solvent_conditions"]["concentration_basis"] = missing(evidence)
    r["steps"] = [
        step("s1_nBuLi_addition", evidence, ["p_start", "p_thf_main", "p_nBuLi"], "add n-BuLi to the vigorously stirred allene solution", (nbu_add_temp, None), atmosphere="argon"),
        step("s2_lithiation_hold", evidence, ["p_start", "p_nBuLi"], "stir after addition", (lith_temp, None), (lith_time, qty(float(lith_time.split()[0]), "min", lith_time)), atmosphere="argon"),
        step("s3_phncs_addition_hold", evidence, ["p_phncs", "p_thf_ncs"], "cool, add PhNCS/THF, and stir", (ncs_temp, None), (ncs_time, qty(float(ncs_time.split()[0]), "min", ncs_time)), atmosphere="argon"),
        step("s4_mei_addition", evidence, ["p_mei"], "cool to about −80 °C and add MeI in one portion", ("ca. −80 °C", qty(-80, "°C", "ca. −80 °C", "approx")), atmosphere="argon"),
        step("s5_cui_addition", evidence, ["p_cui"], "remove cooling bath and add CuI when the stated temperature is reached", (cu_add_temp, None), atmosphere="argon"),
        step("s6_cyclization", evidence, ["p_start", "p_phncs", "p_mei", "p_cui", "p_product"], "stir to cyclize", (heat_temp, None), (heat_time, qty(float(heat_time.split()[0]), "min", heat_time)), atmosphere="argon"),
        step("s7_workup", evidence, ["p_nh4cl", "p_nacn", "p_ether", "p_wash", "p_mgso4"], "cool to room temperature; treat with saturated NH4Cl containing NaCN; separate, extract, wash, dry, and concentrate", ("r.t.", "r.t."), stage="workup"),
        step("s8_purification", evidence, ["p_alumina", "p_hexane", "p_ether_eluent", "p_product"], f"purify on neutral alumina using hexane and {eluent}", stage="purification"),
    ]
    r["outcomes"] = [
        outcome("o_crude_residue", "product_amount", [f"m_12{key}"], residue_mass, qty(float(residue_mass.split()[0]), "g", residue_mass), evidence, product_role="crude product", basis="total residue described as consisting of the product by 1H NMR" if key == "a" else "total residue before chromatography"),
        outcome("o_product_mass", "product_amount", [f"m_12{key}"], product_mass, qty(float(product_mass.split()[0]), "g", product_mass), evidence, product_role="main product", basis="isolated product"),
        outcome("o_yield", "yield", [f"m_12{key}"], f"{percent}%", qty(percent, "%", f"{percent}%"), evidence, product_role="main product", yield_type="isolated yield", basis="isolated product", scope={"starting_material": "PhNCS", "product": f"12{key}"}),
    ]
    # Preserve the Scheme 3 CuBr versus experimental CuI discrepancy explicitly.
    alt_scheme = src("CuBr", "CuBr", ev(4))
    alt_exp = src(f"CuI ({cui_mass})", "CuI", evidence)
    r["notes"] = [fact(raw={"Scheme 3": "CuBr", "experimental section": f"CuI ({cui_mass})"}, normalized=None, evidence=list(dict.fromkeys(ev(4) + evidence)), status="conflict", verification="unreviewed", alternatives=[alt_scheme, alt_exp])]
    set_common_reaction_fields(r, evidence, f"Section 4.{4 if key == 'a' else 5}, pages {pages[0]}–{pages[-1]}" if pages[0] != pages[-1] else f"Section 4.{4 if key == 'a' else 5}, page {pages[0]}", {"section": f"4.{4 if key == 'a' else 5}", "exact_charges": True}, f"Exact preparation of comparator pyrrole 12{key}", ["Scheme 3 labels the catalyst CuBr, while the detailed experimental paragraph uses CuI; both source values are retained as a conflict and the exact procedure participant is CuI.", "Complete product SMILES and atom mapping remain unresolved."])


comparator_record("a", (11, 12), "m_allyloxyallene", "5.00 g", "52.0 mmol", "70 mL", "32 mL", 1.6, "51.2 mmol", "−90 °C", "ca. −70 to −65 °C", "10 min", "6.76 g", "50.0 mmol", "10 mL", "−70 to −65 °C", "10 min", "14.00 g", "98.6 mmol", "2.00 g", "−10 °C", "40–45 °C", "30 min", "10 mL saturated aqueous NH4Cl solution", "~2.00 g", "4 × 30 mL", "60 mL saturated aqueous NH4Cl solution", "12.10 g", "5.49 g", 45, "hexane/Et2O 10:1")
comparator_record("b", (12,), "m_benzyloxyallene", "2.47 g", "16.9 mmol", "25 mL", "6 mL", 2.5, "15.0 mmol", "−93 °C", "ca. −52 to −45 °C", "5 min", "2.03 g", "15.0 mmol", "5 mL", "−52 to −42 °C", "5 min", "4.30 g", "30.3 mmol", "0.60 g", "−5 °C", "33–45 °C", "35 min", "15 mL saturated aqueous NH4Cl solution", "0.50 g", "15 mL", "20 mL water", "4.70 g", "2.22 g", 50, "hexane/Et2O 10:1")

# Preparative CDCl3 rearrangement of 2f is distinct from the Table 2 monitoring sample.
bundle["reactions"] = [x for x in bundle["reactions"] if x["reaction_id"] != f"{DOC}:prep_7f_cdcl3_8d"]
base = copy.deepcopy(by_id(bundle, "thermal_rearrangement_2f"))
base["reaction_id"] = f"{DOC}:prep_7f_cdcl3_8d"
base["record_type"] = "specific_experiment"
base["experiment_categories"] = src(["mechanistic_study"], ["mechanistic_study"], ev(9))
base["location"] = src("Section 4.3.6.1(a), page 9", "Section 4.3.6.1(a), page 9", ev(9))
base["reaction_name"] = src("Preparative rearrangement of 2f to 7f in CDCl3", "Preparative rearrangement of 2f to 7f in CDCl3", ev(9))
base["transformation_type"] = src("formal [1,3]-O-to-C rearrangement", "formal [1,3]-O-to-C rearrangement", ev(9), method="inference", derivation="Inferred from the named starting material and isolated product.")
p_start = participant("p_start", "m_2f", "reactant", ev(9)); charge(p_start, "mass", "0.11 g", qty(0.11, "g", "0.11 g"), ev(9)); charge(p_start, "amount", "0.27 mmol", qty(0.27, "mmol", "0.27 mmol"), ev(9))
p_solvent = participant("p_cdcl3", "m_cdcl3", "solvent", ev(9)); charge(p_solvent, "volume", "1 mL", qty(1, "mL", "1 mL"), ev(9))
p_product = participant("p_product", "m_7f", "product", ev(9), product_role="main product"); charge(p_product, "mass", "0.05 g", qty(0.05, "g", "0.05 g"), ev(9))
p_ether = participant("p_ether", "m_ether", "solvent", ev(9), stage="purification")
p_alumina = participant("p_alumina", "m_alumina", "purification_agent", ev(9), stage="purification"); p_alumina["charges"]["amount"] = src("1 cm layer", None, ev(9))
base["participants"] = [p_start, p_solvent, p_product, p_ether, p_alumina]
base["reaction_scale"] = src("0.27 mmol 2f", {"reference_molecule_id": "m_2f", "amount": qty(0.27, "mmol", "0.27 mmol")}, ev(9))
base["solvent_conditions"]["composition"] = src("CDCl3", ["CDCl3"], ev(9)); base["solvent_conditions"]["mixture_ratio"] = missing(ev(9)); base["solvent_conditions"]["total_volume"] = src("1 mL", qty(1, "mL", "1 mL"), ev(9)); base["solvent_conditions"]["concentration"] = missing(ev(9)); base["solvent_conditions"]["concentration_basis"] = missing(ev(9))
base["steps"] = [step("s1_storage", ev(9), ["p_start", "p_cdcl3", "p_product"], "keep the CDCl3 solution at room temperature", ("r.t.", "r.t."), ("8 days", qty(8, "d", "8 days"))), step("s2_evaporate", ev(9), ["p_start", "p_product"], "evaporate CDCl3", stage="workup"), step("s3_alumina_cleanup", ev(9), ["p_ether", "p_alumina", "p_product"], "dissolve residue in Et2O, pass through a 1 cm neutral-alumina layer, and remove Et2O", stage="purification")]
base["outcomes"] = [outcome("o_product_mass", "product_amount", ["m_7f"], "0.05 g", qty(0.05, "g", "0.05 g"), ev(9), product_role="main product", basis="isolated after alumina cleanup"), outcome("o_yield", "yield", ["m_7f"], "45%", qty(45, "%", "45%"), ev(9), product_role="main product", yield_type="isolated yield", basis="isolated product", scope={"starting_material": "2f", "product": "7f"})]
base["procedure"]["original_text"] = src({"section": "4.3.6.1(a)", "2f": "0.11 g, 0.27 mmol", "CDCl3": "1 mL", "time": "8 days", "product": "0.05 g, 45%"}, "CDCl3 storage followed by alumina cleanup", ev(9)); base["procedure"]["general_procedure_reference"] = missing(ev(9)); base["procedure"]["modifications"] = missing(ev(9)); base["procedure"]["addition_order"] = src("2f in CDCl3 → 8 days at r.t. → evaporate → Et2O → neutral alumina → evaporate", "2f/CDCl3 storage → alumina cleanup", ev(9)); base["procedure"]["pretreatment"] = missing(ev(9))
base["evidence_ids"] = ev(9)
base["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [], "unresolved": ["The source does not report conversion or recovery for this isolated 8-day sample; it is not merged with the separate Table 2 time-course sample.", "Complete organometallic SMILES and atom mapping remain unresolved."]}
base["completeness"] = "partial"
bundle["reactions"].append(base)

# Enrich the thermal time course with its page-9 preparative quantities and isolated outcome.
r = by_id(bundle, "thermal_rearrangement_2f")
evidence = ev(5, 9)
p_start = participant("p_start", "m_2f", "reactant", evidence); charge(p_start, "mass", "0.21 g", qty(0.21, "g", "0.21 g"), ev(9)); charge(p_start, "amount", "0.52 mmol", qty(0.52, "mmol", "0.52 mmol"), ev(9))
p_tol = participant("p_toluene", "m_toluene", "solvent", evidence); charge(p_tol, "volume", "2 mL", qty(2, "mL", "2 mL"), ev(9))
p_prod = participant("p_product", "m_7f", "product", evidence, product_role="main product"); charge(p_prod, "mass", "0.18 g", qty(0.18, "g", "0.18 g"), ev(9))
p_ether = participant("p_ether", "m_ether", "solvent", ev(9), stage="purification")
p_alumina = participant("p_alumina", "m_alumina", "purification_agent", ev(9), stage="purification"); p_alumina["charges"]["amount"] = src("1 cm layer", None, ev(9))
r["participants"] = [p_start, p_tol, p_prod, p_ether, p_alumina]
r["reaction_scale"] = src("0.52 mmol 2f", {"reference_molecule_id": "m_2f", "amount": qty(0.52, "mmol", "0.52 mmol")}, ev(9))
r["solvent_conditions"]["composition"] = src("toluene", ["toluene"], evidence); r["solvent_conditions"]["mixture_ratio"] = missing(evidence); r["solvent_conditions"]["total_volume"] = src("2 mL", qty(2, "mL", "2 mL"), ev(9)); r["solvent_conditions"]["concentration"] = missing(evidence); r["solvent_conditions"]["concentration_basis"] = missing(evidence)
r["steps"] = [step("s1_3h_monitor", evidence, ["p_start", "p_toluene", "p_product"], "warm and monitor by 1H NMR", ("98–100 °C", qrange(98, 100, "°C", "98–100 °C")), ("3 h", qty(3, "h", "3 h"))), step("s2_8h_total", evidence, ["p_start", "p_toluene", "p_product"], "continue warming to a total of 8 h", ("98–100 °C", qrange(98, 100, "°C", "98–100 °C")), ("8 h total", qty(8, "h", "8 h total"))), step("s3_cleanup", ev(9), ["p_ether", "p_alumina", "p_product"], "add Et2O, pass through a 1 cm neutral-alumina layer, and remove solvents", stage="purification")]
r["outcomes"] = [outcome("o_conversion_3h", "conversion", ["m_2f"], "45%", qty(45, "%", "45%"), ev(5), basis="1H NMR", scope={"time": "3 h", "product": "7f"}), outcome("o_conversion_8h", "conversion", ["m_2f"], "96%", qty(96, "%", "96%"), ev(5), basis="1H NMR", scope={"time": "8 h", "product": "7f"}), outcome("o_product_mass", "product_amount", ["m_7f"], "0.18 g", qty(0.18, "g", "0.18 g"), ev(9), product_role="main product", basis="solid after alumina cleanup"), outcome("o_yield", "yield", ["m_7f"], "79%", qty(79, "%", "79%"), ev(9), product_role="main product", yield_type="isolated yield", basis="source-reported yield after alumina cleanup", scope={"starting_material": "2f", "product": "7f"}), outcome("o_residual_2f", "other", ["m_2f"], "containing 8% pyrrole 2f", qty(8, "%", "8%"), ev(9), basis="composition of isolated solid; method not reported", scope={"matrix": "7f-rich solid", "component": "2f"})]
r["procedure"]["original_text"] = src({"section": "4.3.6.1(b)", "2f": "0.21 g, 0.52 mmol", "toluene": "2 mL", "temperature": "98–100 °C", "time": "8 h", "product": "0.18 g, 79%, containing 8% 2f"}, "thermal time-course and preparative isolation", evidence); r["procedure"]["general_procedure_reference"] = missing(evidence); r["procedure"]["modifications"] = missing(evidence); r["procedure"]["addition_order"] = src("2f/toluene → warm 3 h/8 h → Et2O → neutral alumina → evaporate", "thermal monitoring → alumina cleanup", evidence); r["procedure"]["pretreatment"] = missing(evidence)
r["evidence_ids"] = evidence
r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [], "unresolved": ["The 96% 1H-NMR conversion at 8 h and the isolated solid containing 8% 2f use different measurement bases; both are retained rather than reconciled arithmetically.", "The analytical method for the 8% residual-2f composition is not reported.", "Complete organometallic SMILES and atom mapping remain unresolved."]}

# Exact page-12 derivatization 2c -> 17, including the distinct 30 min monitoring statement.
r = by_id(bundle, "scheme7_2c_to_17")
evidence = ev(6, 12)
specs = [("p_start", "m_2c", "reactant", "reaction", [("mass", "0.40 g", qty(0.40, "g", "0.40 g")), ("amount", "1.0 mmol", qty(1.0, "mmol", "1.0 mmol"))]), ("p_tbuok", "m_tbuok", "reagent", "reaction", [("mass", "0.29 g", qty(0.29, "g", "0.29 g")), ("amount", "2.6 mmol", qty(2.6, "mmol", "2.6 mmol"))]), ("p_dmso", "m_dmso", "solvent", "reaction", [("volume", "3 mL", qty(3, "mL", "3 mL"))]), ("p_water", "m_water", "workup_agent", "workup", [("volume", "10 mL", qty(10, "mL", "10 mL"))]), ("p_ether", "m_ether", "workup_agent", "workup", [("volume", "3 × 5 mL", {"number": 5, "unit": "mL", "raw_text": "3 × 5 mL"})]), ("p_water_wash", "m_water", "workup_agent", "workup", [("volume", "3 × 5 mL", {"number": 5, "unit": "mL", "raw_text": "3 × 5 mL"})]), ("p_mgso4", "m_mgso4", "workup_agent", "workup", []), ("p_hexane", "m_hexane", "solvent", "purification", []), ("p_product", "m_17", "product", "reaction", [("mass", "0.338 g", qty(0.338, "g", "0.338 g"))])]
ps = []
for pid, mid, role, stage, cs in specs:
    p = participant(pid, mid, role, evidence, stage=stage, product_role="main product" if pid == "p_product" else None)
    for fld, raw, norm in cs: charge(p, fld, raw, norm, ev(12))
    ps.append(p)
r["participants"] = ps
r["reaction_scale"] = src("1.0 mmol 2c", {"reference_molecule_id": "m_2c", "amount": qty(1.0, "mmol", "1.0 mmol")}, ev(12))
r["solvent_conditions"]["composition"] = src("DMSO", ["DMSO"], ev(12)); r["solvent_conditions"]["mixture_ratio"] = missing(evidence); r["solvent_conditions"]["total_volume"] = src("3 mL", qty(3, "mL", "3 mL"), ev(12)); r["solvent_conditions"]["concentration"] = missing(evidence); r["solvent_conditions"]["concentration_basis"] = missing(evidence)
r["steps"] = [step("s1_add_2c", ev(12), ["p_tbuok", "p_dmso", "p_start"], "add 2c to a stirred solution of t-BuOK in DMSO", ("r.t.", "r.t.")), step("s2_reaction_hold", ev(12), ["p_start", "p_tbuok", "p_dmso", "p_product"], "stir at room temperature", ("r.t.", "r.t."), ("50 min", qty(50, "min", "50 min"))), step("s3_workup", ev(12), ["p_water", "p_ether", "p_water_wash", "p_mgso4"], "dilute with water; extract with Et2O; wash combined organic solution with water; dry over MgSO4 and concentrate", stage="workup"), step("s4_recrystallization", ev(12), ["p_hexane", "p_product"], "isolate 17 by recrystallization from hexane", stage="purification")]
r["outcomes"] = [outcome("o_monitoring_conversion", "conversion", ["m_2c"], "quantitative transformation after 30 min", {"number": 100, "unit": "%", "raw_text": "quantitative transformation after 30 min"}, ev(12), basis="monitoring; analytical method not reported", scope={"time": "30 min", "product": "17"}), outcome("o_crude_mass", "product_amount", ["m_17"], "0.35 g residue", qty(0.35, "g", "0.35 g residue"), ev(12), product_role="crude product", basis="solid residue after concentration"), outcome("o_product_mass", "product_amount", ["m_17"], "0.338 g", qty(0.338, "g", "0.338 g"), ev(12), product_role="main product", basis="isolated after recrystallization"), outcome("o_yield", "yield", ["m_17"], "96%", qty(96, "%", "96%"), ev(12), product_role="main product", yield_type="isolated yield", basis="isolated recrystallized product", scope={"starting_material": "2c", "product": "17"})]
r["procedure"]["original_text"] = src({"section": "4.6", "t-BuOK": "0.29 g, 2.6 mmol", "DMSO": "3 mL", "2c": "0.40 g, 1.0 mmol", "reaction": "r.t., 50 min", "workup": "water 10 mL; Et2O 3×5 mL; water wash 3×5 mL; MgSO4", "product": "0.338 g, 96%"}, "Exact t-BuOK/DMSO preparation of 17", ev(12)); r["procedure"]["general_procedure_reference"] = missing(evidence); r["procedure"]["modifications"] = missing(evidence); r["procedure"]["addition_order"] = src("t-BuOK/DMSO → 2c → 50 min at r.t. → water → Et2O extraction/wash → MgSO4 → hexane recrystallization", "t-BuOK/DMSO → 2c → workup → recrystallization", ev(12)); r["procedure"]["pretreatment"] = missing(evidence)
r["evidence_ids"] = evidence
r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [], "unresolved": ["The monitoring method for the quantitative 30 min conversion is not reported.", "Scheme 7 depicts loss as CH2=CHOK whereas the prose calls it elimination of vinyl alcohol; the coproduct identity is not asserted.", "Complete organometallic SMILES and atom mapping remain unresolved."]}

# Characterization passages remain source context, never synthetic experiments.
for item in [
    {"context_id": "ctx_page7_characterization", "type": "characterization", "evidence_ids": ev(7), "summary": "Elemental analysis, IR, NMR, and EI-MS data for 2a/2b and NMR assignments for 7a/7b; no separate synthesis record."},
    {"context_id": "ctx_page8_characterization", "type": "characterization", "evidence_ids": ev(8), "summary": "Elemental analysis, IR, NMR, and EI-MS blocks for 2c–2f and NMR assignments for 7c–7e; characterization only."},
    {"context_id": "ctx_page9_characterization", "type": "characterization", "evidence_ids": ev(9), "summary": "IR/NMR/EI-MS characterization for 2f/7f/2g/7g/2h; peak lists, molecular ions, and fragment abundances are not synthesis experiments."},
    {"context_id": "ctx_page10_characterization", "type": "characterization", "evidence_ids": ev(10), "summary": "NMR/IR/MS characterization for 7h, 2i, 7i, 5j, 2j, and 7j; analytical peak lists are retained as context only."},
    {"context_id": "ctx_page11_characterization", "type": "characterization", "evidence_ids": ev(11), "summary": "NMR/IR/elemental-analysis data for 5k, 2k, 7k, and 2l; no standalone mass-spectrum experiment is created."},
    {"context_id": "ctx_page12_characterization", "type": "characterization_and_nonexperimental_endmatter", "evidence_ids": ev(12), "summary": "NMR/IR/MS characterization for 12a, 12b, 17, and 18, followed by declarations, acknowledgments, supplementary-material notice, and references; these are not additional synthesis experiments."},
]:
    if not any(x["context_id"] == item["context_id"] for x in bundle["source_context"]):
        bundle["source_context"].append(item)

# Raw page inventory and required coverage artifact.
inventory = {
    "document_id": DOC,
    "pages_reviewed": [7, 8, 9, 10, 11, 12],
    "visual_review": "original-resolution JPEG page render inspected for every page",
    "pages": {
        "7": ["section 4.3 general procedure", "exact preparations 2a and 2b", "characterization blocks for 2a/7a/2b/7b"],
        "8": ["exact preparations 2c, 2d, 2e and start of 2f", "characterization blocks for products and rearrangement products"],
        "9": ["completion of 2f preparation", "two preparative 2f→7f procedures", "exact preparations 2g and 2h", "characterization blocks"],
        "10": ["exact preparation 2i", "5j formation and 5j→2j cyclization", "characterization blocks", "start of 5k formation"],
        "11": ["completion of 5k formation and 5k→2k cyclization", "two-stage 2l preparation", "start of 12a preparation", "characterization blocks"],
        "12": ["completion of 12a preparation", "exact 12b and 17 preparations", "18 characterization", "nonexperimental endmatter/references"],
    },
}
write(RUN / "raw_model" / "page_inventory_0007_0012.json", inventory)

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 7, "state": "attempted", "regions": [{"region": "upper half, section 4.3", "finding": "General one-pot lithiation/isothiocyanate/S-alkylation/Cu(I)-cyclization procedure, workup, chromatography, two-step exception, CDCl3 storage, and thermal 2f note.", "evidence_ids": ev(7)}, {"region": "middle and lower columns, sections 4.3.1–4.3.2", "finding": "Exact preparations of 2a and 2b; 7a/7b paragraphs are characterization assignments only.", "evidence_ids": ev(7)}], "record_ids": [f"{DOC}:scheme_table1_onepot", f"{DOC}:table1_2a", f"{DOC}:table1_2b"], "unresolved": ["Row-specific cyclization time/temperature remains the reported general range."]},
        {"page": 8, "state": "attempted", "regions": [{"region": "full page, sections 4.3.3–4.3.6", "finding": "Exact preparations of 2c, 2d, and 2e and the beginning of 2f; intervening 7c–7e and product spectra are characterization only.", "evidence_ids": ev(8)}], "record_ids": [f"{DOC}:table1_2c", f"{DOC}:table1_2d", f"{DOC}:table1_2e", f"{DOC}:table1_2f"], "unresolved": ["The 2f procedure continues onto page 9."]},
        {"page": 9, "state": "attempted", "regions": [{"region": "top-left continuation of 4.3.6", "finding": "Completion, isolation, and 41% yield of 2f.", "evidence_ids": ev(9)}, {"region": "section 4.3.6.1", "finding": "Two concrete 2f→7f preparations: CDCl3/r.t./8 d (45% isolated) and toluene/98–100 °C/8 h (79% isolated, 8% residual 2f); the thermal experiment is merged with the page-5 NMR time course.", "evidence_ids": ev(9)}, {"region": "lower half, sections 4.3.7–4.3.8", "finding": "Exact preparations of 2g (45%, ~6% 7g contaminant) and 2h (18%); spectra are characterization.", "evidence_ids": ev(9)}], "record_ids": [f"{DOC}:table1_2f", f"{DOC}:prep_7f_cdcl3_8d", f"{DOC}:thermal_rearrangement_2f", f"{DOC}:table1_2g", f"{DOC}:table1_2h"], "unresolved": ["The analytical method for the 8% residual-2f and ~6% 7g composition statements is not reported."]},
        {"page": 10, "state": "attempted", "regions": [{"region": "upper third, section 4.3.9", "finding": "Exact preparation of 2i (50%); 7h/2i/7i spectra are characterization.", "evidence_ids": ev(10)}, {"region": "middle and lower section 4.3.10", "finding": "Two-preparative-step 2j experiment: crude E/Z 5j formation and separate CuBr cyclization to 2j; crude composition and major:minor ratio retained.", "evidence_ids": ev(10)}, {"region": "bottom start of section 4.3.11", "finding": "Beginning of the analogous 5k formation; continued on page 11.", "evidence_ids": ev(10)}], "record_ids": [f"{DOC}:table1_2i", f"{DOC}:table1_2j", f"{DOC}:table1_2k"], "unresolved": ["For 5j/5k, the heading says Z but the crude material is an E/Z mixture and the major configuration is not assigned.", "The Table 1 versus detailed-paragraph basis of the 31% 2j value remains scope-ambiguous."]},
        {"page": 11, "state": "attempted", "regions": [{"region": "upper half, section 4.3.11", "finding": "Completion of crude E/Z 5k formation and CuBr cyclization to 2k; 87% crude and 20% isolated outcomes retained separately.", "evidence_ids": ev(11)}, {"region": "middle-lower section 4.3.12", "finding": "Two-stage low-yield bromoacetonitrile route to 2l, including tar/unreacted-allene crude composition and 3% isolated product.", "evidence_ids": ev(11)}, {"region": "bottom start of section 4.4", "finding": "Beginning of exact comparator synthesis 12a; continued on page 12.", "evidence_ids": ev(11)}], "record_ids": [f"{DOC}:table1_2k", f"{DOC}:table1_2l", f"{DOC}:scheme3_synthesis_12a"], "unresolved": ["5l amount and geometry are not quantified.", "The Table 1 versus detailed-paragraph basis of the 20% 2k value remains scope-ambiguous."]},
        {"page": 12, "state": "attempted", "regions": [{"region": "top continuation of section 4.4", "finding": "Completion, workup, and 45% isolated yield for 12a.", "evidence_ids": ev(12)}, {"region": "middle section 4.5", "finding": "Exact synthesis of 12b in 50% isolated yield.", "evidence_ids": ev(12)}, {"region": "lower-middle section 4.6", "finding": "Exact t-BuOK/DMSO conversion of 2c to 17: quantitative monitoring at 30 min and 96% isolated after a 50 min preparative hold.", "evidence_ids": ev(12)}, {"region": "lower section 4.6.1 and endmatter", "finding": "Compound 18 contains characterization data only on this page; declarations, SI notice, and references are nonexperimental context.", "evidence_ids": ev(12)}], "record_ids": [f"{DOC}:scheme3_synthesis_12a", f"{DOC}:scheme3_synthesis_12b", f"{DOC}:scheme7_2c_to_17"], "unresolved": ["Scheme 3 shows CuBr whereas both detailed comparator procedures use CuI; conflict retained.", "The analytical method for the 30 min quantitative conversion of 2c is not stated."]},
    ],
}
write(RUN / "coverage" / "pages_0007_0012.json", coverage)

write(SUBMITTED, bundle)
write(RUN / "raw_model" / "submitted_pages_0007_0012.json", bundle)
print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]), "source_context": len(bundle["source_context"]), "updated_records": 17, "new_records": 1}, ensure_ascii=False))
