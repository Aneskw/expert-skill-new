#!/usr/bin/env python3
"""Build the source-derived v3.1 bundle for DOI 10.1039/d1qo00459j.

This is a deterministic transcription of the page-1--6 visual inventory.  It
does not certify structures or chemistry; the merge step quarantines all
records as needs_review.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st  # noqa: E402


DOC = "doc_d31888e9faa3be67"
DATASET = ROOT / "runs/current/dataset.json"
OUT = ROOT / "runs/current/submitted.json"


def ev(page: int, *kinds: str) -> list[str]:
    return [f"{DOC}:page{page}:{kind}" for kind in kinds]


PAGE = {p: ev(p, "text", "image") for p in range(1, 7)}


def f(raw=None, normalized=None, evidence=(), method="original_extraction",
      status=None, derivation=None):
    if normalized is None and raw is not None and method == "original_extraction":
        normalized = raw
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   derivation=derivation)


def missing(evidence=(), status="unreported"):
    return st.fact(evidence=evidence, status=status)


bundle = st.empty(json.loads(DATASET.read_text())["documents"])
molecules: dict[str, dict] = {}


def molecule_id(label: str) -> str:
    safe = (label.replace("′", "prime").replace("'", "prime")
            .replace(" ", "_").replace("[", "").replace("]", "")
            .replace("+", "plus").replace("/", "_"))
    return f"{DOC}:mol_{safe}"


def ensure_molecule(label: str, evidence: list[str], *, name: str | None = None,
                    smiles: str | None = None, allene: bool = False,
                    identity_note: str | None = None) -> str:
    mid = molecule_id(label)
    if mid in molecules:
        occ = molecules[mid]["occurrences"][0]
        occ["evidence_ids"] = list(dict.fromkeys(occ["evidence_ids"] + evidence))
        occ["compound_label"]["evidence_ids"] = occ["evidence_ids"]
        return mid
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    m["names"]["original"] = f(name or label, name or label, evidence)
    m["names"]["english"] = missing(evidence)
    m["names"]["abbreviations"] = f(label, [label], evidence)
    m["names"]["aliases"] = missing(evidence)
    m["cas_numbers"] = missing(evidence)
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = f(
            {"source": "name/formula shown on page", "candidate_smiles": smiles},
            smiles,
            evidence,
            method="structure_derivation",
            derivation="Candidate connectivity encoded from the explicitly named simple molecule; parsing is not source-structure certification.",
        )
        m["stereochemistry"]["chirality_status"] = f(
            None, "not_applicable", evidence, method="structure_derivation",
            derivation="The candidate graph has no tetrahedral or axial stereogenic element."
        )
        m["stereochemistry"]["geometry_status"] = f(
            None, "not_applicable", evidence, method="structure_derivation",
            derivation="No source-relevant E/Z element occurs in the candidate graph."
        )
        m["stereochemistry"]["representation_status"] = f(
            None, "not_applicable", evidence, method="structure_derivation",
            derivation="No stereochemical encoding is applicable to this candidate graph."
        )
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing(evidence)
        if allene:
            m["stereochemistry"]["axial"] = f(
                "Allene axis shown in the source drawing without an R_a/S_a or P/M descriptor.",
                {"axis_type": "allene", "configuration": "unspecified"}, evidence)
        m["stereochemistry"]["chirality_status"] = f(
            "No absolute/relative stereochemical descriptor is printed for this drawing.",
            "unspecified", evidence)
        m["stereochemistry"]["geometry_status"] = f(
            "No E/Z or cis/trans descriptor is printed for this drawing.",
            "unspecified", evidence)
        m["stereochemistry"]["representation_status"] = f(
            "Complex graphical structure not safely converted to an isomeric SMILES from the supplied main-text image.",
            "not_encoded", evidence)
    m["stereochemistry"]["isomeric_smiles"] = missing(evidence)
    m["stereochemistry"]["identity_key"] = missing(evidence)
    m["identity_note"] = f(
        identity_note or "Identity retained by paper-local compound label and source drawing; exact structure remains pending expert review.",
        identity_note or "Identity retained by paper-local compound label and source drawing; exact structure remains pending expert review.",
        evidence)
    molecules[mid] = m
    return mid


def template_molecule(label: str, evidence: list[str], fixed: str,
                      variables: list[tuple[str, str]]) -> str:
    mid = molecule_id(label)
    if mid in molecules:
        return mid
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    m["names"]["original"] = f(label, label, evidence)
    m["names"]["abbreviations"] = f(label, [label], evidence)
    m["representation_type"] = "generic_template"
    m["smiles"] = missing(evidence)
    m["generic_template"] = {
        "fixed_structure": f({"graphical_region": fixed}, fixed, evidence),
        "variables": [
            {
                "label": var,
                "attachment_points": f("Source graphical attachment point", "source graphical attachment point", evidence),
                "definition": f(definition, definition, evidence),
                "allowed_values": missing(evidence),
                "constraints": missing(evidence),
            }
            for var, definition in variables
        ],
        "instance_mapping": [],
        "reason": f(
            "The scheme uses R variables and/or a curved tether; it is a generic source template, not a concrete molecule.",
            "generic source drawing", evidence),
    }
    m["stereochemistry"]["chirality_status"] = f(
        "Generic allene/template drawing has no printed absolute configuration.",
        "unresolved", evidence)
    m["stereochemistry"]["geometry_status"] = f(
        "Generic template has no printed E/Z or cis/trans descriptor.",
        "unresolved", evidence)
    m["stereochemistry"]["representation_status"] = f(
        "Stereochemistry cannot be fully encoded for an R-group template.",
        "unresolved", evidence)
    m["identity_note"] = f(
        "Generic graphical structure; variables are preserved without expansion.",
        "Generic graphical structure; variables are preserved without expansion.", evidence)
    molecules[mid] = m
    return mid


TOLUENE = ensure_molecule("toluene", sum(PAGE.values(), []), name="toluene", smiles="Cc1ccccc1")
AGSBF6 = ensure_molecule("AgSbF6", PAGE[3] + PAGE[4], name="silver hexafluoroantimonate")


def participant(pid: str, mid: str, role: str, evidence: list[str], *,
                product_role: str | None = None, mol_percent: float | None = None,
                reference: str | None = None):
    p = st.new_participant(pid, mid, role, evidence)
    if product_role is not None:
        p["product_role"] = f(product_role, product_role, evidence)
    if mol_percent is not None:
        p["charges"]["mol_percent"] = f(
            f"{mol_percent:g} mol%", {"number": mol_percent, "unit": "mol%"}, evidence)
        p["equivalence_reference"] = f(reference, reference, evidence)
    return p


def set_classification(r: dict, mode: str, evidence: list[str]):
    c = r["classification"]
    if mode == "none":
        c["labels"] = f("No observed product-forming transformation.", [], evidence,
                         method="inference", derivation="The source reports decomposition/no desired product; expected-product bonds are not classified.")
        c["coverage"] = f("No observed transformation", "not_applicable", evidence,
                           method="inference", derivation="Negative/decomposition experiment.")
        c["limitations"] = ["No desired product was observed; no intended bond formation is asserted."]
        return
    changes = []
    if mode == "alder_ene":
        changes = [
            {"change_id": "c1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
             "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "c2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}],
             "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
        ]
    elif mode in {"two_plus_two", "diels_alder"}:
        changes = [
            {"change_id": "c1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
             "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "c2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
             "before_order": 0, "after_order": 1, "change_type": "formed"},
        ]
    else:
        raise ValueError(mode)
    labels = {}
    for change in changes:
        lab = st.bond_label(change["atoms"][0]["element"], change["atoms"][1]["element"], change["after_order"])
        labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(change["change_id"])
    derivation = (
        "Partial net-bond annotation from the source reaction class and drawings. Exact atom correspondence and the complete HDDA/aromatization bond inventory are unresolved; no atom-mapped reaction is claimed."
    )
    c["bond_changes"] = [f(None, x, evidence, method="inference", derivation=derivation) for x in changes]
    c["labels"] = f(None, list(labels.values()), evidence, method="inference", derivation=derivation)
    c["coverage"] = f(None, "partial", evidence, method="inference", derivation=derivation)
    c["limitations"] = [
        "Partial classification only; paper-local structures are not atom mapped.",
        "HDDA/aromatization and all bond-order changes are not fully enumerated.",
        "Implicit-H transfer is a net annotation and does not identify the proton-transfer mechanism.",
    ]


def make_specific(rid: str, page_evidence: list[str], substrate: str, *,
                  product: str | None, yield_percent: float | None,
                  reaction_name: str, mode: str, temperature: float,
                  time_value: dict | None, location: str,
                  isolated: bool, category: str = "substrate_scope",
                  silver: bool = False, source_note: str | None = None,
                  product_allene: bool = False):
    smid = ensure_molecule(substrate, page_evidence, allene=(substrate != "5f"))
    pmid = ensure_molecule(product, page_evidence, allene=product_allene) if product else None
    r = st.new_reaction(f"{DOC}:{rid}", DOC, page_evidence)
    r["record_type"] = "specific_experiment"
    r["experiment_categories"] = f(category, [category], page_evidence)
    r["location"] = f(location, {"source_location": location}, page_evidence)
    r["reaction_name"] = f(reaction_name, reaction_name, page_evidence)
    r["transformation_type"] = f(reaction_name, reaction_name, page_evidence)
    r["participants"] = [participant("substrate", smid, "reactant", page_evidence)]
    if pmid:
        r["participants"].append(participant("product", pmid, "product", page_evidence, product_role="main_product"))
    r["participants"].append(participant("solvent", TOLUENE, "solvent", page_evidence))
    if silver:
        r["participants"].append(participant("AgSbF6", AGSBF6, "additive", page_evidence,
                                             mol_percent=5.0, reference=smid))
    r["solvent_conditions"]["composition"] = f("toluene", [TOLUENE], page_evidence)
    step = st.new_step("reaction", page_evidence)
    step["participant_ids"] = [p["participant_id"] for p in r["participants"] if p["role"] != "product"]
    raw_operation = "toluene, 55 °C, 5 mol% AgSbF6" if silver else (
        f"toluene, {temperature:g} °C" + (f", {time_value.get('raw_text')}" if time_value else ""))
    step["operation"] = f(raw_operation, "heated in toluene", page_evidence)
    step["temperature"] = f(f"{temperature:g} °C", {"number": temperature, "unit": "°C"}, page_evidence)
    step["time"] = f(time_value["raw_text"], time_value, page_evidence) if time_value else missing(page_evidence)
    r["steps"] = [step]
    if pmid and yield_percent is not None:
        out = st.new_outcome("yield", "yield", [pmid], page_evidence)
        out["product_role"] = f("reported product", "main_product", page_evidence)
        out["value"] = f(f"{yield_percent:g}%", {"number": yield_percent, "unit": "%"}, page_evidence)
        if isolated:
            out["yield_type"] = f("Isolated yield", "isolated", page_evidence)
            out["measurement_basis"] = f("Isolated yield", "isolated", page_evidence)
        else:
            out["yield_type"] = missing(page_evidence)
            out["measurement_basis"] = missing(page_evidence)
        out["scope"] = f(f"{substrate} to {product}", {"from": substrate, "to": product}, page_evidence)
        r["outcomes"] = [out]
    else:
        out = st.new_outcome("decomposition", "decomposition", [], page_evidence)
        out["value"] = f(
            {"table_or_scheme_result": "0%", "prose_observation": "decomposition of the starting material"},
            "decomposition", page_evidence)
        out["scope"] = f(substrate, {"starting_material": substrate}, page_evidence)
        r["outcomes"] = [out]
    if source_note:
        r["notes"].append(f(source_note, source_note, page_evidence))
    unresolved = [
        "Paper-local complex substrate/product structures remain graphical and require expert structure encoding.",
        "Exact atom mapping and complete net bond-change inventory are unresolved.",
        "Substrate charge, reaction scale, solvent volume/concentration, atmosphere, pressure, workup and purification are unreported in the supplied main text.",
        "No SI file was supplied; the ESI footnote is not treated as an SI association.",
    ]
    if silver and time_value is None:
        unresolved.append("The AgSbF6/55 °C footnote does not report a reaction time.")
    if pmid and not isolated:
        unresolved.append("Yield measurement basis is not reported for Scheme 2.")
    if not pmid:
        unresolved.append("No product participant is created because the source reports decomposition/no desired product.")
    r["review"]["unresolved"] = unresolved
    set_classification(r, mode if pmid else "none", page_evidence)
    bundle["reactions"].append(r)
    return r["reaction_id"]


record_pages = {p: [] for p in range(1, 7)}


# Page 1: the three black reactions are prior literature summaries.  The blue
# box is a conceptual branch map rather than a reported experiment.
g1_ev = PAGE[1]
g1_r = template_molecule("Scheme1_eq4_generic_allene-aryne", g1_ev,
                         "allene tethered to an aryne-forming aryl framework", [("R1", "source-defined substituent")])
g1_p = template_molecule("Scheme1_eq4_generic_outcomes", g1_ev,
                         "three possible products: [2+2], allenic-C-H ene, allylic-C-H ene", [("R1", "source-defined substituent")])
r = st.new_reaction(f"{DOC}:scheme1_eq4_competing_modes", DOC, g1_ev)
r["record_type"] = "generic_scheme"
r["experiment_categories"] = f("Conceptual current-work reaction map", ["other"], g1_ev)
r["location"] = f("PDF page 1, Scheme 1, blue boxed eqn (4)", {"page": 1, "scheme": "1", "equation": 4}, g1_ev)
r["reaction_name"] = f("Alder–ene reaction or [2 + 2] cycloaddition", "competing Alder–ene/[2+2] modes", g1_ev)
r["transformation_type"] = f("Three alternative conceptual branches", "alternative generic transformations", g1_ev)
r["participants"] = [participant("generic_input", g1_r, "reactant", g1_ev), participant("generic_outputs", g1_p, "product", g1_ev, product_role="alternative_products")]
r["review"]["unresolved"] = ["Generic conceptual map only; it is not a concrete experiment and has no unique atom mapping."]
r["classification"]["coverage"] = f("Multiple alternative branches", "unresolved", g1_ev, method="inference", derivation="A single generic record represents three alternative branches; no one bond inventory applies to all branches.")
r["classification"]["labels"] = f("Alternative branches not collapsed into one label set", [], g1_ev, method="inference", derivation="Avoids merging mutually exclusive transformations.")
r["classification"]["limitations"] = ["Alternative outcomes are not a single mapped reaction."]
bundle["reactions"].append(r)
record_pages[1].append(r["reaction_id"])


# Page 2, Scheme 2.
for args in [
    ("scheme2_1a_to_2a", "1a", "2a", 73, "Alder–ene reaction", "alder_ene", "PDF page 2, Scheme 2, arrow 1"),
    ("scheme2_1b_decomposition", "1b", None, None, "decomposition", "none", "PDF page 2, Scheme 2, arrow 2"),
    ("scheme2_1c_to_2c", "1c", "2c", 58, "[2 + 2] cycloaddition", "two_plus_two", "PDF page 2, Scheme 2, arrow 3"),
]:
    rid, sub, prod, yld, name, mode, loc = args
    note = None
    if sub == "1b":
        note = "The drawing prints 2b, 0%; the prose states decomposition of 1b, so 2b is not entered as an observed product."
    if sub == "1c":
        note = "The prose states that only cycloadduct 2c was formed and the Alder–ene product was not observed."
    record_pages[2].append(make_specific(rid, PAGE[2], sub, product=prod, yield_percent=yld,
        reaction_name=name, mode=mode, temperature=90,
        time_value={"number": 8, "unit": "h", "raw_text": "8 h"}, location=loc,
        isolated=False, category="other", source_note=note))


# Page 3, Table 1.  Footnote c is a distinct AgSbF6/55 °C experiment for rows 2 and 3.
t1 = [
    ("1", "3a", "4a", 63, "Alder–ene reaction", "alder_ene"),
    ("2", "3b", None, None, "decomposition", "none"),
    ("3", "3c", None, None, "decomposition", "none"),
    ("4", "3d", "4d", 64, "Alder–ene reaction of an allenic C–H bond", "alder_ene"),
    ("5", "3d′", "4d′", 60, "Alder–ene reaction of an allenic C–H bond", "alder_ene"),
    ("6", "3e", "4e", 55, "Alder–ene reaction of an allenic C–H bond", "alder_ene"),
    ("7", "3f", "4f", 62, "[2 + 2] cycloaddition", "two_plus_two"),
    ("8", "3g", "4g", 34, "Alder–ene reaction of an allylic C–H bond", "alder_ene"),
    ("9", "3h", "4h", 58, "Alder–ene reaction of an allylic C–H bond", "alder_ene"),
]
t1_ev = list(dict.fromkeys(PAGE[3] + ev(2, "text")))
for entry, sub, prod, yld, name, mode in t1:
    note = None
    if not prod:
        note = f"Table 1 entry {entry} reports 0% and decomposition; intended products shown in the table are not entered as observed products."
    rid = make_specific(f"table1_entry{entry}_thermal", t1_ev, sub, product=prod, yield_percent=yld,
        reaction_name=name, mode=mode, temperature=90,
        time_value={"minimum": 6, "maximum": 8, "unit": "h", "raw_text": "6–8 h"},
        location=f"PDF page 3, Table 1, entry {entry}", isolated=True, source_note=note)
    record_pages[3].append(rid)
    if entry in {"2", "3"}:
        rid = make_specific(f"table1_entry{entry}_AgSbF6", t1_ev, sub, product=None, yield_percent=None,
            reaction_name="decomposition under AgSbF6 conditions", mode="none", temperature=55,
            time_value=None, location=f"PDF page 3, Table 1, entry {entry}, footnote c",
            isolated=False, silver=True,
            source_note="Footnote c reports decomposition in toluene at 55 °C with 5 mol% AgSbF6; this is retained as a condition-distinct experiment.")
        record_pages[3].append(rid)


# Page 4, Table 2.  Entries 1 and 2 also have distinct footnote-c silver runs;
# entries 7 and 8 use footnote-d conditions instead of the table header.
t2 = [
    ("1", "5a", None, None, "decomposition", "none", 90, {"minimum": 6, "maximum": 7, "unit": "h", "raw_text": "6–7 h"}),
    ("2", "5b", None, None, "decomposition", "none", 90, {"minimum": 6, "maximum": 7, "unit": "h", "raw_text": "6–7 h"}),
    ("3", "5c", "6c", 72, "Alder–ene reaction of an allenic C–H bond", "alder_ene", 90, {"minimum": 6, "maximum": 7, "unit": "h", "raw_text": "6–7 h"}),
    ("4", "5d", "6d", 50, "Alder–ene reaction", "alder_ene", 90, {"minimum": 6, "maximum": 7, "unit": "h", "raw_text": "6–7 h"}),
    ("5", "5e", "6e", 34, "Diels–Alder reaction to benzobarrelene", "diels_alder", 90, {"minimum": 6, "maximum": 7, "unit": "h", "raw_text": "6–7 h"}),
    ("6", "5f", "6f", 58, "Diels–Alder reaction to benzobarrelene", "diels_alder", 90, {"minimum": 6, "maximum": 7, "unit": "h", "raw_text": "6–7 h"}),
    ("7", "5g", "6g", 32, "Alder–ene reaction", "alder_ene", 150, {"number": 12, "unit": "h", "raw_text": "12 h"}),
    ("8", "5h", "6h", 62, "Alder–ene reaction", "alder_ene", 150, {"number": 12, "unit": "h", "raw_text": "12 h"}),
]
t2_ev = list(dict.fromkeys(PAGE[4] + ev(2, "text") + ev(3, "text")))
for entry, sub, prod, yld, name, mode, temp, timev in t2:
    note = None
    if entry in {"1", "2"}:
        note = f"Table 2 entry {entry} reports 0% and decomposition; the drawn expected product is not treated as observed."
    elif entry == "5":
        note = "The prose states that the NTs toluene moiety, rather than the allene, participates in the Diels–Alder reaction; the allene remains intact in 6e."
    elif entry == "7":
        note = "The authors attribute the marginal yield to instability of the propiolate ester at elevated temperature."
    rid = make_specific(f"table2_entry{entry}_thermal", t2_ev, sub, product=prod, yield_percent=yld,
        reaction_name=name, mode=mode, temperature=temp, time_value=timev,
        location=f"PDF page 4, Table 2, entry {entry}", isolated=True,
        source_note=note, product_allene=(entry == "5"))
    record_pages[4].append(rid)
    if entry in {"1", "2"}:
        rid = make_specific(f"table2_entry{entry}_AgSbF6", t2_ev, sub, product=None, yield_percent=None,
            reaction_name="decomposition under AgSbF6 conditions", mode="none", temperature=55,
            time_value=None, location=f"PDF page 4, Table 2, entry {entry}, footnote c",
            isolated=False, silver=True,
            source_note="Footnote c reports decomposition in toluene at 55 °C with 5 mol% AgSbF6; this is a condition-distinct experiment.")
        record_pages[4].append(rid)


# Page 4, Table 3.
t3 = [
    ("1", "7a", "8a", 82), ("2", "7b", "8b", 85), ("3", "7c", "8c", 80),
    ("4", "7d", "8d", 78), ("5", "7e", "8e", 65), ("6", "7f", "8f", 78),
    ("7", "7ga", None, None), ("8", "7gb", None, None),
    ("9", "7gc", None, None), ("10", "7gd", None, None),
]
t3_ev = PAGE[4]
for entry, sub, prod, yld in t3:
    note = None
    if not prod:
        note = "The monoyne substrate decomposed and no expected [2 + 2] cycloadduct was observed."
    if entry == "10":
        note += " Table 3 prints expected-product label 8gc in row 10, while the prose refers collectively to 8ga–8gd; no observed product is created and the label discrepancy remains unresolved."
    rid = make_specific(f"table3_entry{entry}", t3_ev, sub, product=prod, yield_percent=yld,
        reaction_name="[2 + 2] cycloaddition" if prod else "decomposition",
        mode="two_plus_two" if prod else "none", temperature=150,
        time_value={"number": 12, "unit": "h", "raw_text": "12 h"},
        location=f"PDF page 4, Table 3, entry {entry}", isolated=True, source_note=note)
    record_pages[4].append(rid)


# Page 5, Scheme 3: five generic trend arrows, not concrete experiments.
scheme3 = [
    ("5", "Type-I Alder–ene reaction of gem-dimethyl allenes", "alder_ene", [("R", "H or carbon substituent")]),
    ("6", "Type-I Alder–ene reaction at an allenic C–H bond", "alder_ene", [("R", "H")]),
    ("7", "Type-I Alder–ene reaction at an allylic C–H bond", "alder_ene", [("R", "non-hydrogen substituent")]),
    ("8", "[2 + 2] cycloaddition at the distal allene π-bond", "two_plus_two", [("R", "source-defined substituent")]),
    ("9", "[2 + 2] cycloaddition of a 1,3-diyne-tethered allene", "two_plus_two", [("R", "source-defined alkyne substituent"), ("R1/R2/R3", "allene substituents")]),
]
for eqn, name, mode, variables in scheme3:
    ge = PAGE[5]
    inp = template_molecule(f"Scheme3_eq{eqn}_input", ge, f"Scheme 3 equation ({eqn}) generic input drawing", variables)
    prod = template_molecule(f"Scheme3_eq{eqn}_product", ge, f"Scheme 3 equation ({eqn}) generic product drawing", variables)
    r = st.new_reaction(f"{DOC}:scheme3_eq{eqn}", DOC, ge)
    r["record_type"] = "generic_scheme"
    r["experiment_categories"] = f("General selectivity trend", ["other"], ge)
    r["location"] = f(f"PDF page 5, Scheme 3, equation ({eqn})", {"page": 5, "scheme": 3, "equation": int(eqn)}, ge)
    r["reaction_name"] = f(name, name, ge)
    r["transformation_type"] = f(name, name, ge)
    r["participants"] = [participant("generic_input", inp, "reactant", ge), participant("generic_product", prod, "product", ge, product_role="generic_product")]
    r["review"]["unresolved"] = ["Generic trend only; R variables and exact atom mapping are unresolved.", "No concrete reaction conditions or yield belong to this summary arrow."]
    set_classification(r, mode, ge)
    bundle["reactions"].append(r)
    record_pages[5].append(r["reaction_id"])


bundle["molecules"] = list(molecules.values())
bundle["source_context"] = [
    {
        "context_id": f"{DOC}:metadata",
        "kind": "bibliographic_metadata",
        "title": "Selectivity between an Alder–ene reaction and a [2 + 2] cycloaddition in the intramolecular reactions of allene-tethered arynes",
        "authors": ["Anh Le", "Daesung Lee"],
        "doi": "10.1039/d1qo00459j",
        "evidence_ids": PAGE[1],
    },
    {
        "context_id": f"{DOC}:scheme1_prior_work",
        "kind": "literature_context",
        "description": "Scheme 1 equations (1)–(3) summarize Cheng (2006), Yin (2013), and Lee (2015); they are cited prior work, not current-paper experiments.",
        "evidence_ids": PAGE[1],
    },
    {
        "context_id": f"{DOC}:characterization_inventory_pages1_6",
        "kind": "coverage_note",
        "description": "No standalone mass-spectrum peak table, molecular-weight table, fragment-relative-abundance dataset, or other characterization-only reaction candidate was found on inspected pages 1–6.",
        "evidence_ids": sum(PAGE.values(), []),
    },
    {
        "context_id": f"{DOC}:si_status",
        "kind": "source_availability",
        "description": "Page 1 states that ESI is available, but no SI file was supplied or explicitly associated in this session.",
        "evidence_ids": PAGE[1],
    },
]


OUT.write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({
    "molecules": len(bundle["molecules"]),
    "reactions": len(bundle["reactions"]),
    "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in bundle["reactions"]),
    "generic_schemes": sum(r["record_type"] == "generic_scheme" for r in bundle["reactions"]),
    "record_pages": record_pages,
}, ensure_ascii=False, indent=2))
