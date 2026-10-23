#!/usr/bin/env python3
"""Build the v3.1 incremental submission for source PDF pages 7--12.

The script reuses the current v3 registry, enriches existing scope records
instead of duplicating them, and adds the seven specific 5a--5g acetylations.
It does not modify validators, schemas, tools, or scientific rules.
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

import standard_records as sr  # noqa: E402


DOC = "doc_ee191f89470cf89e"


def ev(page: int, kind: str = "image") -> str:
    return f"{DOC}:page{page}:{kind}"


def evidence(*pages: int) -> list[str]:
    return [item for page in pages for item in (ev(page, "image"), ev(page, "text"))]


def q(number, unit, qualifier=None, raw_text=None):
    value = {"number": number, "unit": unit}
    if qualifier is not None:
        value["qualifier"] = qualifier
    if raw_text is not None:
        value["raw_text"] = raw_text
    return value


def rng(low, high, unit, raw_text):
    return {"minimum": low, "maximum": high, "unit": unit, "raw_text": raw_text}


def qraw(value):
    """Render a quantity object as its source-like text, not a Python dict."""
    if not isinstance(value, dict):
        return str(value)
    if value.get("raw_text"):
        return value["raw_text"]
    if "minimum" in value and "maximum" in value:
        return f"{value['minimum']}-{value['maximum']} {value['unit']}"
    qualifier = {"lt": "<", "le": "≤", "gt": ">", "ge": "≥", "approx": "~"}.get(value.get("qualifier"), "")
    return f"{qualifier}{value['number']} {value['unit']}"


def sf(raw=None, normalized=None, evs=(), method="original_extraction",
       status=None, derivation=None, alternatives=()):
    return sr.fact(raw=raw, normalized=normalized, evidence=evs, method=method,
                   status=status, derivation=derivation, alternatives=alternatives)


def unique(values):
    return list(dict.fromkeys(values))


def merge_evidence(fact, evs):
    for eid in evs:
        if eid not in fact["evidence_ids"]:
            fact["evidence_ids"].append(eid)


def replace_or_append(items, key, value):
    for index, item in enumerate(items):
        if item.get(key) == value.get(key):
            items[index] = value
            return
    items.append(value)


template = json.loads((RUN / "raw_model/pages_0007_0012_init.json").read_text())
current = json.loads((RUN / "standard_records.json").read_text())
bundle = copy.deepcopy(template)
bundle["molecules"] = copy.deepcopy(current["molecules"])
bundle["reactions"] = []
bundle["source_context"] = []

molecules = {m["molecule_id"]: m for m in bundle["molecules"]}
reactions = {r["reaction_id"]: copy.deepcopy(r) for r in current["reactions"]}


def mid(label: str) -> str:
    return f"{DOC}:mol:{label}"


def touch_molecule(label, evs):
    molecule = molecules[mid(label)]
    occurrence = molecule["occurrences"][0]
    occurrence["evidence_ids"] = unique(occurrence["evidence_ids"] + list(evs))
    merge_evidence(occurrence["compound_label"], evs)
    merge_evidence(molecule["names"]["original"], evs)
    return molecule["molecule_id"]


def new_molecule(label, evs, name=None, smiles=None, profile="unresolved"):
    molecule_id = mid(label)
    if molecule_id in molecules:
        touch_molecule(label, evs)
        return molecule_id
    m = sr.new_molecule(molecule_id, DOC, evs)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = sf(label, label, evs)
    m["names"]["original"] = sf(name or label, name or label, evs)
    m["identity_note"] = sf(
        "Identity is registered from the source name/label; machine structure assignments remain subject to expert review.",
        "Identity is registered from the source name/label; machine structure assignments remain subject to expert review.",
        evs,
    )
    st = m["stereochemistry"]
    if smiles:
        diag = sr.stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = sf(
            smiles,
            diag["canonical_isomeric_smiles"],
            evs,
            method="structure_derivation",
            derivation="Candidate SMILES for an unambiguous common reagent was checked with standard stereo; graph diagnostics are not source verification.",
        )
        st["atom_centers"] = sf(
            None, diag["atom_centers"], evs, method="structure_derivation",
            derivation="RDKit standard stereo inspection of the candidate SMILES.")
        st["double_bonds"] = sf(
            None, diag["double_bonds"], evs, method="structure_derivation",
            derivation="RDKit standard stereo inspection of the candidate SMILES.")
        st["chirality_status"] = sf(
            None, "achiral", evs, method="structure_derivation",
            derivation="The candidate graph contains no tetrahedral or axial stereogenic element.")
        st["geometry_status"] = sf(
            None, "not_applicable", evs, method="structure_derivation",
            derivation="The candidate graph has no source-relevant configurable E/Z or ring cis/trans element.")
        st["representation_status"] = sf(
            None, "not_applicable", evs, method="structure_derivation",
            derivation="No stereochemical encoding is applicable to this candidate graph.")
        st["isomeric_smiles"] = sf(
            None, diag["canonical_isomeric_smiles"], evs, method="structure_derivation",
            derivation="Canonical isomeric SMILES from standard stereo inspection.")
        st["identity_key"] = sf(
            None, diag["stereo_graph_key"], evs, method="structure_derivation",
            derivation="Stereo graph key from standard stereo inspection; not a sample identity.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = sf()
        st["representation_status"] = sf(
            None, "not_encoded", evs, method="inference",
            derivation="No complete source-verified SMILES was assigned from the inspected region.")
        if profile == "allenol_center_unspecified":
            st["atom_centers"] = sf(
                "The hydroxy-bearing carbon is structurally stereogenic, but no R/S, wedge, racemate, or composition assignment is reported.",
                [{"location": "hydroxy-bearing carbon", "assignment": "unspecified"}],
                evs,
            )
            st["axial"] = sf(
                "The source name/drawing contains a terminal CH2 allene; no axial stereogenicity is assigned.",
                {"axis_type": "allene", "status": "not_stereogenic_terminal_CH2"},
                evs,
                method="inference",
                derivation="A terminal CH2 end gives two identical substituents at one allene terminus.",
            )
            st["chirality_status"] = sf(
                None, "unspecified", evs, method="inference",
                derivation="The hydroxy-bearing stereocenter is unassigned; absence of a descriptor is not treated as achiral or racemic.")
            st["geometry_status"] = sf(
                None, "not_applicable", evs, method="inference",
                derivation="The cumulated allene is assessed under axial stereochemistry, not conventional E/Z geometry.")
        elif profile == "acetoxy_center_unspecified":
            st["atom_centers"] = sf(
                "The acetoxy-bearing carbon is structurally stereogenic, but no R/S, wedge, racemate, or composition assignment is reported.",
                [{"location": "acetoxy-bearing carbon", "assignment": "unspecified"}],
                evs,
            )
            st["axial"] = sf(
                "The source name/drawing contains a terminal CH2 allene; no axial stereogenicity is assigned.",
                {"axis_type": "allene", "status": "not_stereogenic_terminal_CH2"},
                evs,
                method="inference",
                derivation="A terminal CH2 end gives two identical substituents at one allene terminus.",
            )
            st["chirality_status"] = sf(
                None, "unspecified", evs, method="inference",
                derivation="The acetoxy-bearing stereocenter is unassigned; absence of a descriptor is not treated as achiral or racemic.")
            st["geometry_status"] = sf(
                None, "not_applicable", evs, method="inference",
                derivation="The cumulated allene is assessed under axial stereochemistry, not conventional E/Z geometry.")
        else:
            st["chirality_status"] = sf(
                None, "unresolved", evs, method="inference",
                derivation="The source item is a mixture/material or an unresolved structure for which molecular chirality cannot be assigned.")
            st["geometry_status"] = sf(
                None, "unresolved", evs, method="inference",
                derivation="The source item is a mixture/material or an unresolved structure for which E/Z and ring geometry cannot be assigned.")
    molecules[molecule_id] = m
    bundle["molecules"].append(m)
    return molecule_id


def set_acetoxy_stereo(label, evs):
    m = molecules[mid(label)]
    st = m["stereochemistry"]
    st["atom_centers"] = sf(
        "Acetoxy-bearing stereogenic carbon is present; no R/S, wedge, racemate, or composition assignment is reported.",
        [{"location": "acetoxy-bearing carbon", "assignment": "unspecified"}], evs)
    st["axial"] = sf(
        "Terminal CH2 allene is present; the allene axis is not stereogenic.",
        {"axis_type": "allene", "status": "not_stereogenic_terminal_CH2"}, evs,
        method="inference",
        derivation="The source name/drawing gives a terminal CH2 allene with two identical substituents at one terminus.")
    st["chirality_status"] = sf(
        None, "unspecified", evs, method="inference",
        derivation="The acetoxy-bearing stereocenter is unassigned; no achiral or racemic claim is made.")
    st["geometry_status"] = sf(
        None, "not_applicable", evs, method="inference",
        derivation="Allene axial assessment is separate from conventional E/Z geometry.")
    st["representation_status"] = sf(
        None, "not_encoded", evs, method="inference",
        derivation="No complete source-verified isomeric SMILES was assigned.")
    m["representation_type"] = "unresolved"
    m["smiles"] = sf()


def set_z_geometry(label, evs):
    m = molecules[mid(label)]
    st = m["stereochemistry"]
    st["double_bonds"] = sf(
        f"Experimental heading explicitly names ({'Z'})-{label}.",
        [{"location": "2'-tosylimino C=N", "assignment": "Z"}], evs)
    st["geometry_status"] = sf(
        None, "specified", evs, method="inference",
        derivation="Normalized from the explicit (Z) stereodescriptor in the experimental heading.")
    st["representation_status"] = sf(
        None, "not_encoded", evs, method="inference",
        derivation="The explicit Z assignment is retained textually; no complete source-verified isomeric SMILES was assigned.")


def get_reaction(suffix):
    return reactions[f"{DOC}:{suffix}"]


def include_reaction(r):
    replace_or_append(bundle["reactions"], "reaction_id", r)


def participant(pid, molecule_id, role, evs, stage="reaction", mass=None,
                amount=None, volume=None, equiv=None, molpct=None, ref=None,
                product_role=None):
    p = sr.new_participant(pid, molecule_id, role, evs, stage=stage)
    if mass is not None:
        p["charges"]["mass"] = sf(qraw(mass), mass, evs)
    if amount is not None:
        p["charges"]["amount"] = sf(qraw(amount), amount, evs)
    if volume is not None:
        p["charges"]["volume"] = sf(qraw(volume), volume, evs)
    if equiv is not None:
        p["charges"]["equivalents"] = sf(qraw(equiv), equiv, evs)
    if molpct is not None:
        p["charges"]["mol_percent"] = sf(qraw(molpct), molpct, evs)
    if ref is not None:
        p["equivalence_reference"] = sf(ref, ref, evs)
    if product_role is not None:
        p["product_role"] = sf(product_role, product_role, evs)
    return p


def add_step(r, sid, evs, pids, operation, temperature=None, time=None,
             stage="reaction"):
    s = sr.new_step(sid, evs, stage=stage)
    s["participant_ids"] = list(pids)
    s["operation"] = sf(operation, operation, evs)
    if temperature is not None:
        s["temperature"] = sf(qraw(temperature), temperature, evs)
    if time is not None:
        s["time"] = sf(qraw(time), time, evs)
    r["steps"].append(s)


def add_outcome(r, oid, kind, target, evs, value, raw, basis=None,
                yield_type=None, product_role="desired product", scope=None):
    o = sr.new_outcome(oid, kind, [target], evs)
    o["value"] = sf(raw, value, evs)
    if basis is not None:
        o["measurement_basis"] = sf(basis, basis, evs)
    if yield_type is not None:
        o["yield_type"] = sf(yield_type, yield_type, evs)
    if product_role is not None:
        o["product_role"] = sf(product_role, product_role, evs)
    if scope is not None:
        o["scope"] = sf(scope, scope, evs)
    r["outcomes"].append(o)


def add_note(r, text, evs, method="original_extraction", derivation=None):
    for note in r["notes"]:
        if note.get("normalized_value") == text:
            merge_evidence(note, evs)
            return
    r["notes"].append(sf(text if method == "original_extraction" else None,
                         text, evs, method=method, derivation=derivation))


def classification(r, evs, bonds):
    c = r["classification"]
    c["mapped_reaction_smiles"] = sf()
    changes = []
    labels = {}
    for index, (a, b, order) in enumerate(bonds, 1):
        cid = f"b{index}"
        change = {
            "change_id": cid,
            "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
            "before_order": 0,
            "after_order": order,
            "change_type": "formed",
        }
        changes.append(sf(
            None, change, evs, method="inference",
            derivation="Formed-bond type inferred from the inspected transformation; atom mapping is not claimed."))
        label = sr.bond_label(a, b, order)
        labels.setdefault(label["code"], dict(label, supporting_change_ids=[]))["supporting_change_ids"].append(cid)
    c["bond_changes"] = changes
    c["labels"] = sf(
        None, list(labels.values()), evs, method="inference",
        derivation="Multi-label taxonomy inferred from the source transformation without full atom mapping.")
    c["coverage"] = sf(
        None, "partial", evs, method="inference",
        derivation="Only source-supportable newly formed heavy-atom bonds are listed; implicit-H changes and full atom balance remain unresolved.")
    c["limitations"] = [
        "No source-verified atom mapping was available.",
        "Implicit-hydrogen changes and complete atom balance remain unresolved.",
    ]
    c["model_votes"] = []


def find_participant(r, pid):
    return next(p for p in r["participants"] if p["participant_id"] == pid)


def update_charge(r, pid, field, value, evs, raw=None):
    find_participant(r, pid)["charges"][field] = sf(raw or qraw(value), value, evs)


def merge_record_evidence(r, evs):
    r["evidence_ids"] = unique(r["evidence_ids"] + list(evs))
    for key in ("location", "reaction_name", "transformation_type"):
        merge_evidence(r[key], evs)


def conflict_fact(alternatives):
    return sf(None, None, (), status="conflict", alternatives=[
        sf(raw, normalized, evs) for raw, normalized, evs in alternatives
    ])


E6, E7, E8, E9, E10, E11, E12 = (evidence(page) for page in range(6, 13))

# Standard stereo checks are applied to the only newly encoded concrete reagent.
hexane = new_molecule("hexane", E7 + E8 + E10 + E11, name="hexane", smiles="CCCCCC")
brine = new_molecule("brine solution", E7 + E8 + E10 + E11, name="brine solution")

# Add page-specific occurrences to existing common agents.
for label, evs in [
    ("DCM", E7 + E8), ("pyridine", E7 + E8), ("acetyl_chloride", E7 + E8),
    ("EtOAc", E7 + E8 + E10 + E11), ("Na2SO4", E7 + E8 + E10 + E11),
    ("H2O", E8 + E10 + E11), ("PhMe", E8 + E10 + E11),
    ("PPh3", E8 + E9 + E10 + E11), ("Cs2CO3", E10 + E11),
    ("pTol3P", E11 + E12),
]:
    touch_molecule(label, evs)


# Complete the cross-page general acetylation procedure begun on page 6.
general5 = get_reaction("page6_general_synthesis_5a_5g")
merge_record_evidence(general5, E7)
general5["location"] = sf(
    "Experimental section iii, pages 6-7",
    "Experimental section iii, pages 6-7", E6 + E7)
general5["procedure"]["original_text"] = sf(
    "General acetylation: allenol in dry dichloromethane at -5 °C; add pyridine, hold 5 min; add acetyl chloride over 5 min; hold 45 min; extract with dichloromethane, wash with brine, dry over Na2SO4, concentrate, and chromatograph (EtOAc/hexane 1:9).",
    "General acetylation with staged -5 °C additions, 45 min hold, aqueous workup, and silica chromatography.",
    E6 + E7)
general5["procedure"]["addition_order"] = sf(
    "allenol/DCM; pyridine; acetyl chloride",
    ["allenol in DCM", "pyridine", "acetyl chloride"], E6 + E7)
general5["review"]["unresolved"] = [
    "Generic aryl/heteroaryl variables are not expanded into a concrete structure or atom mapping.",
    "The amount of drying agent and chromatography solvent volumes are not reported.",
]
general_product = new_molecule(
    "benzyl 2-(acetoxy-4-aryl/heteroaryl)buta-2,3-dienoate generic",
    E6 + E7,
    profile="acetoxy_center_unspecified",
)
if not any(p["participant_id"] == "p_product_family" for p in general5["participants"]):
    general5["participants"].extend([
        participant("p_extract_DCM", mid("DCM"), "workup_agent", E7, stage="workup",
                    volume=q(25, "mL", raw_text="3 × 25 mL")),
        participant("p_brine", brine, "workup_agent", E7, stage="workup",
                    volume=q(20, "mL", raw_text="2 × 20 mL")),
        participant("p_drying", mid("Na2SO4"), "workup_agent", E7, stage="workup"),
        participant("p_column_EtOAc", mid("EtOAc"), "purification_agent", E7, stage="purification"),
        participant("p_column_hexane", hexane, "purification_agent", E7, stage="purification"),
        participant("p_product_family", general_product, "product", E6 + E7,
                    product_role="generic 5a-g product family"),
    ])
# The first four reaction steps came from page 6; explicitly preserve -5 °C for the 45 min hold.
for step in general5["steps"]:
    if step["step_id"] == "post_addition_hold":
        step["temperature"] = sf("same temperature (-5 °C)", q(-5, "°C"), E6)
general5["steps"] = [s for s in general5["steps"] if s["step_id"] not in {"extraction", "wash_dry", "purification"}]
add_step(general5, "extraction", E7, ["p_extract_DCM"],
         "extract aqueous layer with dichloromethane (3 × 25 mL)", stage="workup")
add_step(general5, "wash_dry", E7, ["p_brine", "p_drying"],
         "wash combined organic layer with brine (2 × 20 mL), dry over anhydrous Na2SO4, and concentrate", stage="workup")
add_step(general5, "purification", E7, ["p_column_EtOAc", "p_column_hexane", "p_product_family"],
         "silica gel chromatography using ethyl acetate/hexane 1:9", stage="purification")
classification(general5, E6 + E7, (("C", "O", 1),))
include_reaction(general5)


# Seven independent 5a--5g preparations.  Characterization peak lists are not
# converted to synthetic experiments; the physical form and melting point are
# retained as notes only.
five_series = [
    ("5a", "benzyl 2-(hydroxy(phenyl)methyl)buta-2,3-dienoate", 7, 2.80, 1.84, 57, "gummy liquid", None),
    ("5b", "benzyl 2-(hydroxy(4-methoxyphenyl)methyl)buta-2,3-dienoate", 8, 3.10, 1.86, 53, "gummy liquid", None),
    ("5c", "benzyl 2-((4-chlorophenyl)(hydroxy)methyl)buta-2,3-dienoate", 8, 3.14, 1.99, 56, "gummy liquid", None),
    ("5d", "benzyl 2-(hydroxy(4-(trifluoromethyl)phenyl)methyl)buta-2,3-dienoate", 8, 3.48, 2.30, 60, "gummy liquid", None),
    ("5e", "benzyl 2-(hydroxy(pyren-1-yl)methyl)buta-2,3-dienoate", 8, 4.04, 2.00, 45, "brown solid", (114, 116)),
    ("5f", "benzyl 2-(hydroxy(thiophene-2-yl)methyl)buta-2,3-dienoate", 8, 2.86, 1.87, 57, "gummy liquid", None),
    ("5g", "benzyl 3-hydroxy-5-phenyl-2-vinylidenepentanoate", 8, 3.08, 1.72, 49, "gummy liquid", None),
]
new_five_ids = []
for product, start_name, page, start_mass, product_mass, yld, physical, mp in five_series:
    specific = E7 + E8 if product == "5a" else E8
    product_mid = touch_molecule(product, specific)
    set_acetoxy_stereo(product, specific)
    start_mid = new_molecule(f"{product}_allenol", specific, name=start_name,
                             profile="allenol_center_unspecified")
    rid = f"{DOC}:synthesis_{product}"
    r = sr.new_reaction(rid, DOC, unique(E6 + E7 + specific))
    r["record_type"] = "specific_experiment"
    r["experiment_categories"] = sf(["other"], ["other"], specific)
    r["location"] = sf(f"Experimental section iii, compound {product}",
                       f"Experimental section iii, compound {product}", specific)
    r["reaction_name"] = sf(f"Synthesis of {product}", f"Synthesis of {product}", specific)
    r["transformation_type"] = sf("acetylation of an allenic alcohol",
                                   "acetylation of an allenic alcohol", E6 + E7 + specific)
    r["general_scheme_ids"] = [general5["reaction_id"]]
    r["review"] = {
        "status": "needs_review", "method": "model", "reviewer": None,
        "checked_fields": [],
        "unresolved": [
            "Complete source-verified structures and atom mapping remain unresolved.",
            "No enantiomeric composition or absolute configuration is reported for the hydroxy or acetoxy stereocenter.",
            "Amounts of Na2SO4 and chromatography solvents are unreported.",
        ],
    }
    r["completeness"] = "partial"
    r["participants"] = [
        participant("p_allenol", start_mid, "reactant", specific,
                    mass=q(start_mass, "g"), amount=q(10.0, "mmol"),
                    equiv=q(1.0, "equiv"), ref=start_mid),
        participant("p_pyridine", mid("pyridine"), "reagent", E6 + specific,
                    volume=q(1.20, "mL"), amount=q(15.0, "mmol"),
                    equiv=q(1.5, "equiv"), ref=start_mid),
        participant("p_AcCl", mid("acetyl_chloride"), "reagent", E6 + specific,
                    volume=q(1.10, "mL"), amount=q(15.0, "mmol"),
                    equiv=q(1.5, "equiv"), ref=start_mid),
        participant("p_DCM", mid("DCM"), "solvent", E6, volume=q(20, "mL")),
        participant("p_extract_DCM", mid("DCM"), "workup_agent", E7, stage="workup",
                    volume=q(25, "mL", raw_text="3 × 25 mL")),
        participant("p_brine", brine, "workup_agent", E7, stage="workup",
                    volume=q(20, "mL", raw_text="2 × 20 mL")),
        participant("p_drying", mid("Na2SO4"), "workup_agent", E7, stage="workup"),
        participant("p_column_EtOAc", mid("EtOAc"), "purification_agent", specific, stage="purification"),
        participant("p_column_hexane", hexane, "purification_agent", specific, stage="purification"),
        participant("p_product", product_mid, "product", specific,
                    mass=q(product_mass, "g"), product_role=f"{product}"),
    ]
    r["reaction_scale"] = sf(
        "allenic alcohol (10.0 mmol)",
        {"reference_molecule_id": start_mid, "amount": q(10.0, "mmol")}, E6 + specific)
    r["solvent_conditions"]["composition"] = sf("dry dichloromethane", [mid("DCM")], E6)
    r["solvent_conditions"]["total_volume"] = sf("20 mL", q(20, "mL"), E6)
    r["procedure"]["original_text"] = sf(
        "Following the general procedure: -5 °C pyridine/acetyl chloride acetylation, 45 min hold, aqueous extraction, brine wash, Na2SO4 drying, and silica chromatography (EtOAc/hexane 1:9).",
        "General acetylation procedure with compound-specific masses and isolated product amount.", E6 + E7 + specific)
    r["procedure"]["general_procedure_reference"] = sf(
        "iii. General Procedure for the Synthesis of Compounds 5a-g",
        general5["reaction_id"], E6 + E7 + specific)
    r["procedure"]["addition_order"] = sf(
        "allenol/DCM; pyridine; acetyl chloride",
        ["allenol in DCM", "pyridine", "acetyl chloride"], E6 + E7)
    add_step(r, "charge", E6, ["p_allenol", "p_DCM"],
             "prepare allenol solution in dry dichloromethane", q(-5, "°C"))
    add_step(r, "pyridine_hold", E6, ["p_pyridine"],
             "add pyridine and stir", q(-5, "°C"), q(5.0, "min"))
    add_step(r, "acetyl_chloride_addition", E6, ["p_AcCl"],
             "slowly add acetyl chloride", q(-5, "°C"), q(5.0, "min"))
    add_step(r, "post_addition_hold", E6 + E7, ["p_allenol", "p_pyridine", "p_AcCl"],
             "stir after addition", q(-5, "°C"), q(45, "min"))
    add_step(r, "extraction", E7, ["p_extract_DCM"],
             "extract aqueous layer with dichloromethane (3 × 25 mL)", stage="workup")
    add_step(r, "wash_dry", E7, ["p_brine", "p_drying"],
             "wash with brine (2 × 20 mL), dry over anhydrous Na2SO4, and concentrate", stage="workup")
    add_step(r, "purification", specific, ["p_column_EtOAc", "p_column_hexane", "p_product"],
             "silica gel chromatography using ethyl acetate/hexane 1:9", stage="purification")
    add_outcome(r, "yield", "yield", product_mid, specific, q(yld, "%"), f"{yld}%",
                basis="isolated after silica gel chromatography", yield_type="isolated",
                product_role=product)
    add_outcome(r, "product_amount", "product_amount", product_mid, specific,
                q(product_mass, "g"), f"{product_mass:g} g",
                basis="reported isolated mass", product_role=product)
    physical_note = f"Product {product} was reported as a {physical}."
    if mp:
        physical_note += f" Melting point {mp[0]}-{mp[1]} °C."
    add_note(r, physical_note, specific)
    classification(r, E6 + E7 + specific, (("C", "O", 1),))
    bundle["reactions"].append(r)
    new_five_ids.append(rid)


# Shared workup/purification molecules for the three annulation procedures.
def workup_participants(extract_volume, evs, ratio):
    return [
        participant("p_water", mid("H2O"), "workup_agent", evs, stage="workup", volume=q(10, "mL")),
        participant("p_extract_EtOAc", mid("EtOAc"), "workup_agent", evs, stage="workup",
                    volume=q(extract_volume, "mL", raw_text=f"3 × {extract_volume} mL")),
        participant("p_brine", brine, "workup_agent", evs, stage="workup", volume=q(20, "mL")),
        participant("p_drying", mid("Na2SO4"), "workup_agent", evs, stage="workup"),
        participant("p_column_EtOAc", mid("EtOAc"), "purification_agent", evs, stage="purification"),
        participant("p_column_hexane", hexane, "purification_agent", evs, stage="purification"),
    ]


def replace_product_amount(r, target, mass, evs):
    r["outcomes"] = [o for o in r["outcomes"] if o["outcome_id"] != "product_amount"]
    add_outcome(r, "product_amount", "product_amount", target, evs, mass,
                mass.get("raw_text", f"{mass['number']:g} {mass['unit']}"),
                basis="reported isolated mass")


def enrich_yield(r, evs):
    for outcome in r["outcomes"]:
        if outcome["kind"] == "yield":
            merge_evidence(outcome["value"], evs)
            outcome["measurement_basis"] = sf("isolated yield", "isolated yield", evs)
            outcome["yield_type"] = sf("isolated", "isolated", evs)


def append_workup_steps(r, evs, extract_volume, ratio):
    add_step(r, "quench", evs, ["p_water"], "quench by addition of water (10 mL)", stage="workup")
    add_step(r, "extraction", evs, ["p_extract_EtOAc"],
             f"extract aqueous layer with ethyl acetate (3 × {extract_volume} mL)", stage="workup")
    add_step(r, "wash_dry", evs, ["p_brine", "p_drying"],
             "wash combined organic layer with brine (20 mL), dry over anhydrous Na2SO4, and concentrate", stage="workup")
    add_step(r, "purification", evs, ["p_column_EtOAc", "p_column_hexane", "p_product"],
             f"silica gel chromatography using ethyl acetate/hexane {ratio}", stage="purification")


# Scheme 3 records: merge the explicit experimental paragraphs on pages 8--10.
scheme3_details = [
    ("3aa", 8, 63.0, 62.4, 85.4, 83, (148, 150)),
    ("3ab", 8, 63.0, 70.0, 92.5, 85, (155, 157)),
    ("3ac", 8, 63.0, 70.7, 85.5, 78, (142, 144)),
    ("3ad", 9, 63.0, 79.0, 90.0, 77, (156, 158)),
    ("3ae", 9, 63.0, 70.0, 82.7, 76, (153, 155)),
    ("3af", 9, 63.0, 95.1, 97.5, 75, (132, 134)),
    ("3ag", 9, 63.0, 74.48, 81.2, 72, (141, 143)),
    ("3ah", 9, 63.0, 63.91, 76.0, 73, (135, 137)),
    ("3ai", 9, 63.0, 95.5, 98.0, 75, (114, 116)),
    ("3aj", 9, 63.0, 47.5, 64.2, 71, None),
    ("3ba", 9, 60.0, 62.4, 84.0, 84, (125, 127)),
    ("3ca", 10, 65.3, 62.4, 82.0, 78, (124, 126)),
    ("3da", 10, 75.3, 62.4, 84.1, 73, (124, 126)),
    ("3ea", 10, 78.6, 62.4, 87.6, 74, (160, 162)),
    ("3fa", 10, 78.6, 62.4, 84.0, 71, (149, 151)),
]
scheme3_ids = []
for product, page, m1_mass, m2_mass, product_mass, yld, mp in scheme3_details:
    specific = evidence(page)
    if product == "3ac":
        specific = unique(specific + E9)
    if product == "3ba":
        specific = unique(specific + E10)
    r = get_reaction(f"scheme3_{product}")
    merge_record_evidence(r, E8 + specific)
    r["location"] = sf(f"Scheme 3 and experimental section iv, {product}",
                       f"Scheme 3 and experimental section iv, {product}", E8 + specific)
    update_charge(r, "p_sub1", "mass", q(m1_mass, "mg"), specific)
    update_charge(r, "p_sub2", "mass", q(m2_mass, "mg"), specific)
    update_charge(r, "p_cat", "mass", q(10.5, "mg"), specific)
    update_charge(r, "p_solvent", "volume", q(2.0, "mL"), E8)
    update_charge(r, "p_product", "mass", q(product_mass, "mg"), specific)
    r["participants"] = [p for p in r["participants"] if p["participant_id"] not in {
        "p_water", "p_extract_EtOAc", "p_brine", "p_drying", "p_column_EtOAc", "p_column_hexane"}]
    r["participants"].extend(workup_participants(10, E8, "1:9"))
    r["solvent_conditions"]["total_volume"] = sf("2.0 mL", q(2.0, "mL"), E8)
    r["procedure"]["original_text"] = sf(
        "Charge 2-sulfonamidoindole, delta-acetoxy allenoate, and toluene; add PPh3 at 25 °C; stir 36 h; quench with water; extract with EtOAc; wash with brine; dry; concentrate; chromatograph (EtOAc/hexane 1:9).",
        "Room-temperature annulation, aqueous workup, and silica purification.", E8 + specific)
    r["procedure"]["general_procedure_reference"] = sf(
        "iv. General Procedure for the Synthesis of Compounds 3aa-3aj, 3ba-3da, and 3ea-3fa",
        "experimental section iv", E8)
    r["procedure"]["addition_order"] = sf(
        "2-sulfonamidoindole + delta-acetoxy allenoate + toluene; then PPh3",
        ["2-sulfonamidoindole", "delta-acetoxy allenoate", "toluene", "PPh3"], E8)
    r["steps"] = []
    add_step(r, "charge", E8, ["p_sub1", "p_sub2", "p_solvent"],
             "charge substrates and toluene", q(25, "°C"))
    add_step(r, "catalyst_addition", E8, ["p_cat"], "add PPh3", q(25, "°C"))
    add_step(r, "reaction_hold", E8, ["p_sub1", "p_sub2", "p_cat", "p_solvent"],
             "stir and monitor by TLC", q(25, "°C"), q(36, "h"))
    append_workup_steps(r, E8, 10, "1:9")
    replace_product_amount(r, find_participant(r, "p_product")["molecule_id"], q(product_mass, "mg"), specific)
    enrich_yield(r, specific)
    add_note(r, f"Specific paragraph reports {product} as " +
             (f"a solid, mp {mp[0]}-{mp[1]} °C." if mp else "a gummy liquid; no melting point is reported."), specific)
    r["review"]["unresolved"] = unique(r["review"]["unresolved"] + [
        "Complete source-verified structures, atom mapping, and absolute configuration remain unresolved.",
        "Amounts of drying agent and chromatography solvents are unreported.",
    ])
    include_reaction(r)
    scheme3_ids.append(r["reaction_id"])


# Scheme 4 records: preserve the three printed identifier conflicts explicitly.
scheme4_details = [
    ("4aa", 10, 63.0, 62.4, 77.0, 75, (162, 164)),
    ("4ab", 10, 63.0, 70.0, 85.0, 78, (164, 166)),
    ("4ac", 10, 63.0, 71.0, 80.0, 73, (189, 191)),
    ("4ad", 10, 63.0, 79.0, 81.0, 70, (138, 140)),
    ("4af", 11, 63.0, 95.1, 88.0, 68, (176, 178)),
    ("4ah", 11, 63.0, 64.0, 66.3, 64, (148, 150)),
    ("4aj", 11, 63.0, 47.5, 61.2, 68, (133, 135)),
    ("4ak", 11, 63.0, 100.0, 96.2, 72, (115, 117)),
    ("4ba", 11, 60.0, 62.4, 76.0, 76, (130, 132)),
    ("4ea", 11, 79.0, 62.4, 77.0, 65, (202, 204)),
]
scheme4_ids = []
for product, page, m1_mass, m2_mass, product_mass, yld, mp in scheme4_details:
    specific = evidence(page)
    if product == "4ad":
        specific = unique(specific + E11)
    r = get_reaction(f"scheme4_{product}")
    merge_record_evidence(r, E10 + specific)
    r["location"] = sf(f"Scheme 4 and experimental section v, {product}",
                       f"Scheme 4 and experimental section v, {product}", E10 + specific)
    update_charge(r, "p_sub1", "mass", q(m1_mass, "mg", raw_text="60. mg") if product == "4ba" else q(m1_mass, "mg"), specific)
    update_charge(r, "p_sub2", "mass", q(m2_mass, "mg"), specific)
    update_charge(r, "p_cat", "mass", q(10.5, "mg"), specific)
    update_charge(r, "p_base", "mass", q(130.0, "mg"), specific)
    # The source explicitly divides toluene into charge and feed portions.
    if any(p["participant_id"] == "p_solvent" for p in r["participants"]):
        find_participant(r, "p_solvent")["participant_id"] = "p_solvent_charge"
    update_charge(r, "p_solvent_charge", "volume", q(1.0, "mL"), E10)
    r["participants"] = [p for p in r["participants"] if p["participant_id"] not in {
        "p_solvent_feed", "p_water", "p_extract_EtOAc", "p_brine", "p_drying", "p_column_EtOAc", "p_column_hexane"}]
    r["participants"].append(participant("p_solvent_feed", mid("PhMe"), "solvent", E10, volume=q(1.0, "mL")))
    update_charge(r, "p_product", "mass", q(product_mass, "mg"), specific)
    r["participants"].extend(workup_participants(15, E10, "1:9"))
    r["solvent_conditions"]["total_volume"] = sf("2.0 mL total (1.0 mL charge + 1.0 mL feed)", q(2.0, "mL"), E10)
    r["procedure"]["original_text"] = sf(
        "Charge sulfonamidoindole, PPh3, Cs2CO3, and 1.0 mL toluene; add allenoate in 1.0 mL toluene over 30 min at 80 °C; hold 36 h; aqueous workup; silica chromatography (EtOAc/hexane 1:9).",
        "Staged 80 °C allenoate addition, 36 h hold, aqueous workup, and silica purification.", E10 + specific)
    r["procedure"]["general_procedure_reference"] = sf(
        "v. General Procedure for the Synthesis of Compounds 4aa-4ad, 4af, 4ah, 4aj, 4ak, 4ba, and 4ea",
        "experimental section v", E10)
    r["procedure"]["addition_order"] = sf(
        "sulfonamidoindole + PPh3 + Cs2CO3 + toluene; then allenoate/toluene feed",
        ["sulfonamidoindole", "PPh3", "Cs2CO3", "toluene", "allenoate solution"], E10)
    conflict_note = None
    if product == "4aj":
        r["procedure"]["modifications"] = conflict_fact([
            ("Scheme 4 assigns 4aj to allenoate 2j", {"product": mid("4aj"), "substrate": mid("2j")}, evidence(4)),
            ("Specific paragraph prints '2a (47.5 mg, 0.24 mmol)'", {"product": mid("4aj"), "printed_substrate_label": "2a", "mass": q(47.5, "mg")}, E11),
        ])
        conflict_note = "Source conflict: Scheme 4/product identity assigns 2j, while the specific 4aj paragraph prints 2a with 47.5 mg; no silent correction was made."
    elif product == "4ak":
        r["procedure"]["modifications"] = conflict_fact([
            ("Section heading and Scheme 4 identify product 4ak", {"product": mid("4ak")}, evidence(4) + E10),
            ("Specific compound heading ends with '(4aj)' but the following sentence says afforded 4ak", {"printed_heading_label": "4aj", "body_product_label": "4ak"}, E11),
        ])
        conflict_note = "Source conflict: the compound heading is printed as 4aj, but Scheme 4 and the procedure body identify the product as 4ak."
    elif product == "4ba":
        r["procedure"]["modifications"] = conflict_fact([
            ("Scheme 4 assigns 4ba to allenoate 2a", {"product": mid("4ba"), "substrate": mid("2a")}, evidence(4)),
            ("Specific paragraph prints 'with 4a (62.4 mg, 0.24 mmol)'", {"product": mid("4ba"), "printed_substrate_label": "4a", "mass": q(62.4, "mg")}, E11),
        ])
        conflict_note = "Source conflict: Scheme 4 assigns substrate 2a, while the specific 4ba paragraph prints 4a; both statements are retained."
    r["steps"] = []
    add_step(r, "initial_charge", E10, ["p_sub1", "p_cat", "p_base", "p_solvent_charge"],
             "charge sulfonamidoindole, PPh3, Cs2CO3, and toluene")
    add_step(r, "allenoate_addition", E10, ["p_sub2", "p_solvent_feed"],
             "add allenoate solution gradually", q(80, "°C"), q(30, "min"))
    add_step(r, "reaction_hold", E10, ["p_sub1", "p_sub2", "p_cat", "p_base", "p_solvent_charge", "p_solvent_feed"],
             "hold after addition and monitor by TLC", q(80, "°C"), q(36, "h"))
    append_workup_steps(r, E10, 15, "1:9")
    replace_product_amount(r, find_participant(r, "p_product")["molecule_id"], q(product_mass, "mg"), specific)
    enrich_yield(r, specific)
    add_note(r, f"Specific paragraph reports {product} as a solid, mp {mp[0]}-{mp[1]} °C.", specific)
    if conflict_note:
        add_note(r, conflict_note, specific)
        r["review"]["unresolved"].append(conflict_note)
    r["review"]["unresolved"] = unique(r["review"]["unresolved"] + [
        "Full atom mapping, implicit-H changes, and complete structure encoding remain unresolved.",
        "Amounts of drying agent and chromatography solvents are unreported.",
    ])
    include_reaction(r)
    scheme4_ids.append(r["reaction_id"])


# Scheme 5 records with specific procedures visible on pages 11--12.
scheme5_details = [
    ("6aa", 11, 63.0, 77.4, 91.0, 79, (148, 150), True),
    ("6ab", 12, 63.0, 84.5, 100.6, 83, (175, 177), True),
    ("6ac", 12, 63.0, 85.4, 92.7, 76, (170, 172), True),
    ("6ad", 12, 63.0, 94.0, 91.5, 71, (123, 125), True),
    ("6ae", 12, 63.0, 107.1, 91.0, 65, (194, 196), True),
    ("6af", 12, 63.0, 79.0, 81.5, 70, (77, 79), True),
    ("6ag", 12, 63.0, 84.0, 83.4, 69, None, False),
    ("6ba", 12, 60.0, 77.4, 92.2, 82, (105, 107), True),
    ("6ca", 12, 65.3, 77.4, 93.0, 79, (96, 98), True),
]
scheme5_ids = []
for product, page, m1_mass, m2_mass, product_mass, yld, mp, explicit_z in scheme5_details:
    specific = evidence(page)
    r = get_reaction(f"scheme5_{product}")
    merge_record_evidence(r, E11 + specific)
    r["location"] = sf(f"Scheme 5 and experimental section vi, {product}",
                       f"Scheme 5 and experimental section vi, {product}", E11 + specific)
    if explicit_z:
        set_z_geometry(product, specific)
    update_charge(r, "p_sub1", "mass", q(m1_mass, "mg"), specific)
    update_charge(r, "p_sub2", "mass", q(m2_mass, "mg"), specific)
    update_charge(r, "p_cat", "mass", q(12.2, "mg"), specific)
    if any(p["participant_id"] == "p_solvent" for p in r["participants"]):
        find_participant(r, "p_solvent")["participant_id"] = "p_solvent_charge"
    update_charge(r, "p_solvent_charge", "volume", q(1.0, "mL"), E11)
    r["participants"] = [p for p in r["participants"] if p["participant_id"] not in {
        "p_solvent_feed", "p_water", "p_extract_EtOAc", "p_brine", "p_drying", "p_column_EtOAc", "p_column_hexane"}]
    r["participants"].append(participant("p_solvent_feed", mid("PhMe"), "solvent", E11, volume=q(1.0, "mL")))
    update_charge(r, "p_product", "mass", q(product_mass, "mg"), specific)
    r["participants"].extend(workup_participants(15, E11, "1:4"))
    r["solvent_conditions"]["total_volume"] = sf("2.0 mL total (1.0 mL charge + 1.0 mL feed)", q(2.0, "mL"), E11)
    r["procedure"]["original_text"] = sf(
        "Charge sulfonamidoindole, tri(p-tolyl)phosphine, and 1.0 mL toluene; add beta-prime-acetoxyallenoate in 1.0 mL toluene over 30 min at 80 °C; hold 36 h; aqueous workup; silica chromatography (EtOAc/hexane 1:4).",
        "Staged 80 °C allenoate addition, 36 h hold, aqueous workup, and silica purification.", E11 + specific)
    r["procedure"]["general_procedure_reference"] = sf(
        "vi. General Procedure for the Synthesis of Compounds 6aa-6ag, 6ba-6da, 6bf, 6ga, and 6hb",
        "experimental section vi", E11)
    r["procedure"]["addition_order"] = sf(
        "sulfonamidoindole + tri(p-tolyl)phosphine + toluene; then beta-prime-acetoxyallenoate/toluene feed",
        ["sulfonamidoindole", "tri(p-tolyl)phosphine", "toluene", "beta-prime-acetoxyallenoate solution"], E11)
    r["steps"] = []
    add_step(r, "initial_charge", E11, ["p_sub1", "p_cat", "p_solvent_charge"],
             "charge sulfonamidoindole, tri(p-tolyl)phosphine, and toluene")
    add_step(r, "allenoate_addition", E11, ["p_sub2", "p_solvent_feed"],
             "add beta-prime-acetoxyallenoate solution gradually", q(80, "°C"), q(30, "min"))
    add_step(r, "reaction_hold", E11, ["p_sub1", "p_sub2", "p_cat", "p_solvent_charge", "p_solvent_feed"],
             "hold after addition and monitor by TLC", q(80, "°C"), q(36, "h"))
    append_workup_steps(r, E11, 15, "1:4")
    replace_product_amount(r, find_participant(r, "p_product")["molecule_id"], q(product_mass, "mg"), specific)
    enrich_yield(r, specific)
    add_note(r, f"Specific paragraph reports {product} as " +
             (f"a solid, mp {mp[0]}-{mp[1]} °C." if mp else "a gummy liquid; no melting point is reported."), specific)
    if not explicit_z:
        r["review"]["unresolved"].append(
            "The specific 6ag heading does not give a Z/E descriptor; no geometry is inferred from omission.")
    r["review"]["unresolved"] = unique(r["review"]["unresolved"] + [
        "Absolute configuration remains unassigned; relative drawing/dr are not converted into an absolute stereochemical claim.",
        "Diastereomeric-ratio measurement method is not reported in Scheme 5.",
        "Amounts of drying agent and chromatography solvents are unreported.",
    ])
    include_reaction(r)
    scheme5_ids.append(r["reaction_id"])


# The generic family records are required in the same submission because the
# enriched scope records reference them through general_scheme_ids.
for suffix in ("generic_3_series", "generic_4_series", "generic_6_series"):
    include_reaction(get_reaction(suffix))


# Scheme 9 on page 7 is a second source region for the three existing scale-up
# experiments.  Add only facts visible in the graphic; masses/workup remain
# supported by their original extracted-text region.
scaleup_ids = []
for product in ("3aa", "4aa", "6aa"):
    r = get_reaction(f"scaleup_{product}")
    merge_record_evidence(r, E7)
    merge_evidence(r["reaction_scale"], E7)
    r["solvent_conditions"]["concentration"] = sf("PhMe (0.1 M)", q(0.1, "M"), E7)
    r["solvent_conditions"]["concentration_basis"] = sf(
        "Scheme 9 reports PhMe (0.1 M); denominator is not explicitly stated in the graphic.",
        "reported reaction concentration; denominator not explicitly defined", E7)
    for p in r["participants"]:
        if p["participant_id"] in {"p_1a", "p_sub2"}:
            merge_evidence(p["charges"]["amount"], E7)
        if p["participant_id"] == "p_cat":
            merge_evidence(p["charges"]["mol_percent"], E7)
    for s in r["steps"]:
        if s["stage"] == "reaction":
            merge_evidence(s["temperature"], E7)
            merge_evidence(s["time"], E7)
    enrich_yield(r, E7)
    r["location"] = sf(
        f"Scheme 9 (page 7) and scale-up procedure text (page 13), {product}",
        f"Scheme 9 and scale-up procedure, {product}", unique(E7 + [ev(13, "text")]))
    r["review"]["unresolved"] = unique(r["review"]["unresolved"] + [
        "PDF page 13 procedure text remains contextual evidence; page 13 is not included in the visual coverage claim for this batch."
    ])
    include_reaction(r)
    scaleup_ids.append(r["reaction_id"])


bundle["source_context"].extend([
    {
        "context_id": f"{DOC}:scheme8_proposed_pathways_page7",
        "kind": "proposed_mechanism",
        "description": "Scheme 8 depicts proposed intermediates and electron-flow pathways for the [3+3], tosyl-migration, and [4+1] annulations. These are mechanistic hypotheses, not independent experiments.",
        "evidence_ids": E7,
    },
    {
        "context_id": f"{DOC}:analytical_characterization_pages7_12",
        "kind": "characterization_context",
        "description": "NMR, IR, HRMS, physical form, and melting-point text associated with compounds 5a-g, 3-series, 4-series, and 6-series is characterization material. Peak lists and calculated/found masses are not counted as synthetic experiments.",
        "evidence_ids": unique(E7 + E8 + E9 + E10 + E11 + E12),
    },
])


# Coverage must list exactly the six visually inspected pages.  Record IDs are
# assigned by actual page occurrence; analytical-only continuations are stated
# as such rather than inflated into new experiments.
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 7,
            "state": "attempted",
            "regions": [
                {"region": "Scheme 8, upper and middle page", "finding": "Proposed [3+3], tosyl-migration, and [4+1] pathways; retained as proposed-mechanism source context, not experimental records.", "evidence_ids": E7},
                {"region": "Scheme 9, lower left", "finding": "Three existing 1.0 mmol scale-up experiments for 3aa, 4aa, and 6aa; page-7 conditions/yields were merged into those records.", "evidence_ids": E7},
                {"region": "general procedure iii continuation, lower columns", "finding": "45 min hold, DCM extraction, brine wash, Na2SO4 drying, and chromatography complete the cross-page general acetylation procedure.", "evidence_ids": E7},
                {"region": "specific 5a paragraph", "finding": "Independent synthesis of 5a: stated reagent volumes, starting/product masses, and 57% isolated yield; analytical peak list is characterization only.", "evidence_ids": E7},
            ],
            "record_ids": unique([general5["reaction_id"], new_five_ids[0]] + scaleup_ids),
            "unresolved": ["Scheme 8 is a proposed mechanism and does not establish intermediate isolation.", "Specific 5a stereocenter configuration/enantiomeric composition and complete atom mapping are not reported."],
        },
        {
            "page": 8,
            "state": "attempted",
            "regions": [
                {"region": "left column and upper right", "finding": "Specific syntheses 5b-5g with starting/product masses and isolated yields; the top continuation of 5a is analytical characterization, not a new reaction.", "evidence_ids": E8},
                {"region": "general procedure iv", "finding": "Room-temperature PPh3 annulation procedure, including 36 h hold, workup, extraction volumes, and chromatography, merged into each corresponding Scheme 3 experiment.", "evidence_ids": E8},
                {"region": "specific 3aa, 3ab, and start of 3ac paragraphs", "finding": "Compound-specific charges, isolated masses, yields, and physical data were merged into existing Scheme 3 records.", "evidence_ids": E8},
            ],
            "record_ids": unique(new_five_ids + scheme3_ids[:3]),
            "unresolved": ["NMR/IR/HRMS peak lists are retained only as characterization context.", "Absolute configurations of the unassigned stereogenic products remain unresolved."],
        },
        {
            "page": 9,
            "state": "attempted",
            "regions": [
                {"region": "top left continuation", "finding": "Characterization continuation for 3ac; no separate synthetic experiment created.", "evidence_ids": E9},
                {"region": "specific 3ad-3aj, 3ba paragraphs", "finding": "Compound-specific charges, isolated product masses/yields, and physical data merged into existing Scheme 3 records.", "evidence_ids": E9},
            ],
            "record_ids": unique([f"{DOC}:scheme3_3ac"] + scheme3_ids[3:11]),
            "unresolved": ["Peak tables and calculated/found HRMS values are characterization, not reaction outcomes.", "No source-verified complete structures or absolute configurations were encoded."],
        },
        {
            "page": 10,
            "state": "attempted",
            "regions": [
                {"region": "top left continuation", "finding": "Characterization continuation for 3ba; no duplicate experiment created.", "evidence_ids": E10},
                {"region": "specific 3ca-3fa paragraphs", "finding": "Specific charges, isolated masses/yields, workup, and chromatography merged into existing Scheme 3 records.", "evidence_ids": E10},
                {"region": "general procedure v", "finding": "Initial charge, 30 min allenoate feed at 80 °C, 36 h hold, workup, and purification merged into Scheme 4 records as separate ordered steps.", "evidence_ids": E10},
                {"region": "specific 4aa-4ad paragraphs", "finding": "Compound-specific charges and isolated outcomes merged into existing Scheme 4 records; 4ad characterization continues on page 11.", "evidence_ids": E10},
            ],
            "record_ids": unique([f"{DOC}:scheme3_3ba"] + scheme3_ids[11:] + scheme4_ids[:4]),
            "unresolved": ["Complete atom mapping and implicit-hydrogen changes remain unresolved.", "Drying-agent and chromatography-solvent quantities are unreported."],
        },
        {
            "page": 11,
            "state": "attempted",
            "regions": [
                {"region": "top left continuation", "finding": "Characterization continuation for 4ad; no duplicate reaction created.", "evidence_ids": E11},
                {"region": "specific 4af, 4ah, 4aj, 4ak, 4ba, 4ea paragraphs", "finding": "Specific charges and isolated outcomes merged into existing Scheme 4 records; printed identifier conflicts for 4aj/4ak/4ba were retained explicitly.", "evidence_ids": E11},
                {"region": "general procedure vi", "finding": "Initial charge, 30 min beta-prime-acetoxyallenoate feed at 80 °C, 36 h hold, workup, and 1:4 chromatography merged into Scheme 5 records.", "evidence_ids": E11},
                {"region": "specific 6aa paragraph", "finding": "Specific charges, isolated mass/yield, explicit Z heading, and physical data merged into the existing 6aa scope record.", "evidence_ids": E11},
            ],
            "record_ids": unique([f"{DOC}:scheme4_4ad"] + scheme4_ids[4:] + [scheme5_ids[0]]),
            "unresolved": ["4aj paragraph prints 2a although Scheme 4 assigns 2j; 4ak heading prints 4aj; 4ba paragraph prints 4a although Scheme 4 assigns 2a.", "The source does not report a dr measurement method."],
        },
        {
            "page": 12,
            "state": "attempted",
            "regions": [
                {"region": "specific 6ab-6ca paragraphs", "finding": "Eight visible product sections were inventoried; 6ab-6af, 6ag, 6ba, and 6ca procedures/outcomes were merged into existing Scheme 5 records. 6ag lacks an explicit Z/E descriptor and was not assigned one.", "evidence_ids": E12},
                {"region": "NMR/IR/HRMS text throughout", "finding": "Analytical peak lists and calculated/found masses were retained as characterization context, not synthetic experiments.", "evidence_ids": E12},
            ],
            "record_ids": scheme5_ids[1:],
            "unresolved": ["6ag geometry is not explicitly assigned in its heading; omission is not interpreted as achiral, racemic, or a particular E/Z isomer.", "Absolute configurations remain unassigned; relative stereochemistry and dr stay product-specific."],
        },
    ],
}


inventory = {
    "document_id": DOC,
    "batch_id": "allene-stereo-v31-200-20261017",
    "visually_inspected_pages": [7, 8, 9, 10, 11, 12],
    "new_reactions": new_five_ids,
    "enriched_reactions": unique(
        [general5["reaction_id"]] + scheme3_ids + scheme4_ids + scheme5_ids + scaleup_ids),
    "supporting_generic_records_in_submission": [
        f"{DOC}:generic_3_series", f"{DOC}:generic_4_series", f"{DOC}:generic_6_series"],
    "submitted_reaction_count": len(bundle["reactions"]),
    "registered_molecule_count_in_bundle": len(bundle["molecules"]),
    "source_context_added": [c["context_id"] for c in bundle["source_context"]],
    "chemistry_accuracy_verified": False,
    "notes": [
        "All machine records remain needs_review after standard merge.",
        "No supporting-information document was associated.",
        "Characterization peak lists were not converted into reaction experiments.",
        "Printed identifier conflicts were retained rather than silently corrected.",
    ],
}


# Earlier records used Python dict rendering for some quantity raw values.
# Normalize only those mechanical strings in this submission to source-like
# quantity text; retain authored prose, conflicts, and qualitative values.
for reaction in bundle["reactions"]:
    for fact in sr.facts(reaction):
        fact["evidence_ids"] = unique(fact["evidence_ids"])
        normalized = fact.get("normalized_value")
        raw = fact.get("raw_value")
        if (isinstance(normalized, dict) and "number" in normalized and "unit" in normalized
                and isinstance(raw, str) and raw.startswith("{'number':")):
            fact["raw_value"] = qraw(normalized)


(RUN / "coverage").mkdir(parents=True, exist_ok=True)
(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
(RUN / "coverage/pages_0007_0012.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/pages_0007_0012_extraction.json").write_text(json.dumps(inventory, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({
    "new_reactions": len(new_five_ids),
    "enriched_reactions": len(inventory["enriched_reactions"]),
    "submitted_reactions": len(bundle["reactions"]),
    "molecules": len(bundle["molecules"]),
}, ensure_ascii=False))
