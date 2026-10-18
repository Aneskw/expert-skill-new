#!/usr/bin/env python3
"""Incrementally add source-backed PDF pages 7--10 to the v3.1 bundle."""
from __future__ import annotations

import json
import sys
from copy import deepcopy
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
import standard_records as st

DOC = "doc_a413938d5770a9e8"
P = {n: f"{DOC}:page{n}:text" for n in range(7, 11)}
PI = {n: f"{DOC}:page{n}:image" for n in range(7, 11)}


def F(raw=None, norm=None, ev=(), method="original_extraction", derivation=None, status=None):
    return st.fact(raw, norm, list(ev), method=method, derivation=derivation, status=status)


def q(number, unit):
    return {"number": number, "unit": unit}


def add_evidence(fact, evidence_ids):
    if not isinstance(fact, dict) or "evidence_ids" not in fact:
        return
    fact["evidence_ids"] = list(dict.fromkeys(fact["evidence_ids"] + list(evidence_ids)))


template = json.loads((RUN / "raw_model/init_pages_0007_0010.json").read_text())
base = json.loads((RUN / "standard_records.json").read_text())
bundle = deepcopy(template)
# Keep the fresh init document block because it contains all ten page-image
# evidence objects; the previous incremental standard bundle only retained
# images for pages 1--6.
bundle["documents"] = deepcopy(template["documents"])
bundle["molecules"] = deepcopy(base["molecules"])
bundle["reactions"] = deepcopy(base["reactions"])
bundle["source_context"] = deepcopy(base["source_context"])

molecules = {m["molecule_id"]: m for m in bundle["molecules"]}
reactions = {r["reaction_id"]: r for r in bundle["reactions"]}


PRODUCTS = {
    "3c": {
        "page": 7,
        "name": "Dimethyl-(Z)-3-(cyclohex-1-en-1-yl)-4-((Z)-3-(2-ethoxy-2-oxoethyl)-2-(prop-1-en-2-yl)hexa-2,5-dien-1-ylidene)cyclopent-2-ene-1,1-dicarboxylate",
        "assignments": ["cyclopentene exocyclic ylidene: Z", "hexa-2-ene: Z"],
        "methods": {"a": (83, 77)}, "eluent": "ethyl acetate/pentane, 5:95",
    },
    "3d": {
        "page": 7,
        "name": "Dimethyl-(Z)-3-(3,3-dimethylbut-1-en-2-yl)-4-((Z)-3-(2-ethoxy-2-oxoethyl)-2-(prop-1-en-2-yl)hexa-2,5-dien-1-ylidene)cyclopent-2-ene-1,1-dicarboxylate",
        "assignments": ["cyclopentene exocyclic ylidene: Z", "hexa-2-ene: Z"],
        "methods": {"a": (82, 78), "b": (78, 74)}, "eluent": "ethyl acetate/pentane, 5:95",
    },
    "3e": {
        "page": 7,
        "name": "Ethyl-(Z)-3-allyl-4-((Z)-(4,4-bis(methoxymethyl)-2-(prop-1-en-2-yl)cyclopent-2-en-1-ylidene)methyl)-5-methylhexa-3,5-dienoate",
        "assignments": ["ethyl enoate alkene: Z", "ylidene-methyl linkage: Z"],
        "methods": {"a": (80, 64)}, "eluent": "ethyl acetate/pentane, 4:96",
    },
    "3f": {
        "page": 7,
        "name": "Ethyl-(Z)-3-allyl-4-((Z)-(2-(cyclohex-1-en-1-yl)-4,4-bis(methoxymethyl)cyclopent-2-en-1-ylidene)methyl)-5-methylhexa-3,5-dienoate",
        "assignments": ["ethyl enoate alkene: Z", "ylidene-methyl linkage: Z"],
        "methods": {"a": (77, 68)}, "eluent": "ethyl acetate/pentane, 5:95",
    },
    "3g": {
        "page": 7,
        "name": "Ethyl-(Z)-3-allyl-4-((Z)-(2-(3,3-dimethylbut-1-en-2-yl)-4,4-bis(methoxymethyl)cyclopent-2-en-1-ylidene)methyl)-5-methylhexa-3,5-dienoate",
        "assignments": ["ethyl enoate alkene: Z", "ylidene-methyl linkage: Z"],
        "methods": {"a": (81, 71)}, "eluent": "ethyl acetate/pentane, 4:96",
    },
    "3h": {
        "page": 7,
        "name": "Ethyl-(Z)-3-allyl-4-((Z)-(8,8-dimethyl-3-(prop-1-en-2-yl)-7,9-dioxaspiro[4.5]dec-3-en-2-ylidene)methyl)-5-methylhexa-3,5-dienoate",
        "assignments": ["ethyl enoate alkene: Z", "ylidene-methyl linkage: Z"],
        "methods": {"a": (75, 62)}, "eluent": "ethyl acetate/pentane, 5:95",
    },
    "3i": {
        "page": 7,
        "name": "Ethyl-(Z)-3-allyl-4-((Z)-(3-(3,3-dimethylbut-1-en-2-yl)-8,8-dimethyl-7,9-dioxaspiro[4.5]dec-3-en-2-ylidene)methyl)-5-methylhexa-3,5-dienoate",
        "assignments": ["ethyl enoate alkene: Z", "ylidene-methyl linkage: Z"],
        "methods": {"a": (34, 31)}, "eluent": "ethyl acetate/pentane, 5:95",
    },
    "3j": {
        "page": 7,
        "name": "Ethyl-(Z)-3-allyl-4-((Z)-(4,4-bis((benzyloxy)methyl)-2-(prop-1-en-2-yl)cyclopent-2-en-1-ylidene)methyl)-5-methylhexa-3,5-dienoate",
        "assignments": ["ethyl enoate alkene: Z", "ylidene-methyl linkage: Z"],
        "methods": {"a": (83, 91)}, "eluent": "ethyl acetate/pentane, 3:97",
    },
    "3k": {
        "page": 7,
        "name": "Dimethyl-(Z)-4-((Z)-3-(ethoxycarbonyl)-2-(prop-1-en-2-yl)hexa-2,5-dien-1-ylidene)-3-(prop-1-en-2-yl)cyclopent-2-ene-1,1-dicarboxylate",
        "assignments": ["cyclopentene exocyclic ylidene: Z", "hexa-2-ene: Z"],
        "methods": {"a": (85, 70), "b": (87, 72)}, "eluent": "ethyl acetate/pentane, 5:95",
    },
    "3l": {
        "page": 7,
        "name": "Dimethyl-(Z)-4-((Z)-3-(2-ethoxy-2-oxoethyl)-2-(prop-1-en-2-yl-d5)hexa-2,5-dien-1-ylidene)-3-(prop-1-en-2-yl)cyclopent-2-ene-1,1-dicarboxylate",
        "assignments": ["cyclopentene exocyclic ylidene: Z", "hexa-2-ene: Z"],
        "methods": {"a": (72, 62)}, "eluent": "ethyl acetate/pentane, 5:95",
    },
    "3m": {
        "page": 7,
        "name": "Dimethyl-(Z)-4-((Z)-2-(cyclohex-1-en-1-yl)-3-(2-ethoxy-2-oxoethyl)hexa-2,5-dien-1-ylidene)-3-(prop-1-en-2-yl)cyclopent-2-ene-1,1-dicarboxylate",
        "assignments": ["cyclopentene exocyclic ylidene: Z", "hexa-2-ene: Z"],
        "methods": {"a": (78, 74), "b": (80, 76)}, "eluent": "ethyl acetate/pentane, 5:95",
    },
    "3n": {
        "page": 8,
        "name": "Dimethyl-(Z)-4-((Z)-2-(cyclopent-1-en-1-yl)-3-(2-ethoxy-2-oxoethyl)hexa-2,5-dien-1-ylidene)-3-(prop-1-en-2-yl)cyclopent-2-ene-1,1-dicarboxylate",
        "assignments": ["cyclopentene exocyclic ylidene: Z", "hexa-2-ene: Z"],
        "methods": {"a": (75, 69), "b": (77, 71)}, "eluent": "ethyl acetate/pentane, 5:95",
    },
    "3o": {
        "page": 8,
        "name": "Dimethyl-(Z)-4-((Z)-3-(2-ethoxy-2-oxoethyl)-2-vinylhexa-2,5-dien-1-ylidene)-3-(prop-1-en-2-yl)cyclopent-2-ene-1,1-dicarboxylate",
        "assignments": ["cyclopentene exocyclic ylidene: Z", "hexa-2-ene: Z"],
        "methods": {"a": (61, 50)}, "eluent": "ethyl acetate/pentane, 5:95",
    },
    "3p": {
        "page": 8,
        "name": "Dimethyl-(Z)-4-((2Z,5E)-3-(ethoxycarbonyl)-2-(prop-1-en-2-yl)octa-2,5-dien-1-ylidene)-3-(prop-1-en-2-yl)cyclopent-2-ene-1,1-dicarboxylate",
        "assignments": ["cyclopentene exocyclic ylidene: Z", "octa-2-ene: Z", "octa-5-ene: E"],
        "methods": {"a": (76, 67)}, "eluent": "ethyl acetate/pentane, 7:93",
    },
    "3q": {
        "page": 8,
        "name": "Dimethyl-(Z)-4-((2Z,5E)-3-(ethoxycarbonyl)-6-phenyl-2-(prop-1-en-2-yl)hexa-2,5-dien-1-ylidene)-3-(prop-1-en-2-yl)cyclopent-2-ene-1,1-dicarboxylate",
        "assignments": ["cyclopentene exocyclic ylidene: Z", "hexa-2-ene: Z", "hexa-5-ene: E"],
        "methods": {"a": (57, 56)}, "eluent": "ethyl acetate/pentane, 5:95",
    },
    "3r": {
        "page": 8,
        "name": "Dimethyl-(Z)-4-((Z)-3-(ethoxycarbonyl)-6-methyl-2-(prop-1-en-2-yl)hepta-2,5-dien-1-ylidene)-3-(prop-1-en-2-yl)cyclopent-2-ene-1,1-dicarboxylate",
        "assignments": ["cyclopentene exocyclic ylidene: Z", "hepta-2-ene: Z"],
        "methods": {"a": (71, 63)}, "eluent": "ethyl acetate/pentane, 5:95",
    },
    "3s": {
        "page": 8,
        "name": "Dimethyl-(Z)-4-((Z)-3-(2-ethoxy-2-oxoethyl)-5-methyl-2-(prop-1-en-2-yl)hexa-2,5-dien-1-ylidene)-3-(prop-1-en-2-yl)cyclopent-2-ene-1,1-dicarboxylate",
        "assignments": ["cyclopentene exocyclic ylidene: Z", "hexa-2-ene: Z"],
        "methods": {"a": (80, 71)}, "eluent": "ethyl acetate/pentane, 5:95",
    },
    "3t": {
        "page": 8,
        "name": "Dimethyl-(Z)-4-((E)-3-allyl-2-(prop-1-en-2-yl)hept-2-en-1-ylidene)-3-(prop-1-en-2-yl)cyclopent-2-ene-1,1-dicarboxylate",
        "assignments": ["cyclopentene exocyclic ylidene: Z", "hept-2-ene: E"],
        "methods": {"a": (80, 64)}, "eluent": "ethyl acetate/pentane, 3:97",
    },
    "3u": {
        "page": 8,
        "name": "Dimethyl-(Z)-4-((Z)-3-phenyl-2-(prop-1-en-2-yl)hexa-2,5-dien-1-ylidene)-3-(prop-1-en-2-yl)cyclopent-2-ene-1,1-dicarboxylate",
        "assignments": ["characterization heading names Z/Z component", "Table 3 reports E/Z mixture without ratio"],
        "methods": {"a": (60, 50)}, "eluent": "ethyl acetate/pentane, 3:97", "mixture": True,
    },
    "3v": {
        "page": 8,
        "name": "Dimethyl-(Z)-4-((Z)-3-(2-hydroxyethyl)-2-(prop-1-en-2-yl)hexa-2,5-dien-1-ylidene)-3-(prop-1-en-2-yl)cyclopent-2-ene-1,1-dicarboxylate",
        "assignments": ["cyclopentene exocyclic ylidene: Z", "hexa-2-ene: Z"],
        "methods": {"a": (63, 49)}, "eluent": "ethyl acetate/pentane, 10:90",
    },
}


def reaction_id(product, method):
    table = 2 if product[-1] in "cdefghij" else 3
    return f"{DOC}:table{table}_{product}_method_{method}"


for label, info in PRODUCTS.items():
    ev = [P[info["page"]], PI[info["page"]]]
    mid = "m_" + label
    m = molecules[mid]
    for occurrence in m["occurrences"]:
        if occurrence["document_id"] == DOC:
            occurrence["evidence_ids"] = list(dict.fromkeys(occurrence["evidence_ids"] + ev))
            add_evidence(occurrence["compound_label"], ev)
    m["names"]["original"] = F(info["name"], info["name"], ev)
    m["names"]["english"] = F(info["name"], info["name"], ev)
    m["identity_note"] = F(
        f"Characterization heading identifies {label}; no independently checked complete product graph was encoded.",
        f"Source systematic name retained for {label}; full structure remains unresolved.", ev)
    m["stereochemistry"]["double_bonds"] = F(info["name"], info["assignments"], ev)
    m["stereochemistry"]["chirality_status"] = F(
        "No center/axis assignment is reported in the characterization heading", "unresolved", ev,
        method="inference", derivation="Absence of an R/S or axial descriptor is not treated as achiral or racemic.")
    if info.get("mixture"):
        mix_ev = [f"{DOC}:page4:image", f"{DOC}:page4:text"] + ev
        m["stereochemistry"]["mixture"] = F(
            "Table 3: (E:Z); characterization heading names a Z/Z component",
            {"type": "E/Z mixture", "ratio": None, "named_component": "Z/Z"}, mix_ev)
        m["stereochemistry"]["geometry_status"] = F(
            "E/Z mixture; ratio not reported", "mixture", mix_ev,
            method="inference", derivation="The scope table explicitly marks E:Z, while the characterization heading names one Z/Z component.")
    else:
        m["stereochemistry"]["geometry_status"] = F(
            info["name"], "partially_specified", ev, method="inference",
            derivation="Only the locants/configurations explicit in the systematic heading are retained; the complete graph was not encoded.")
    m["stereochemistry"]["representation_status"] = F(
        "Source E/Z descriptors retained textually; no isomeric SMILES encoded", "not_encoded", ev,
        method="inference", derivation="Complex product structure was not independently transcribed from the supplied source.")

    for method, (yield_pct, mass_mg) in info["methods"].items():
        rid = reaction_id(label, method)
        r = reactions[rid]
        r["evidence_ids"] = list(dict.fromkeys(r["evidence_ids"] + ev))
        add_evidence(r["location"], ev)
        for outcome in r["outcomes"]:
            if outcome["kind"] == "yield":
                add_evidence(outcome["value"], ev)
                add_evidence(outcome["yield_type"], ev)
                add_evidence(outcome["measurement_basis"], ev)
        mass_outcome = st.new_outcome("product_mass_characterization", "product_amount", [mid], ev)
        mass_outcome["product_role"] = F("product", "product", ev)
        mass_outcome["value"] = F(f"{mass_mg} mg", q(mass_mg, "mg"), ev)
        mass_outcome["measurement_basis"] = F(
            f"Isolated yield: Method {method.upper()}, {yield_pct}% ({mass_mg} mg); column chromatography",
            "isolated product mass after column chromatography", ev)
        r["outcomes"].append(mass_outcome)

        pids = {p["participant_id"] for p in r["participants"]}
        for pid, pmid in (("purif_silica", "m_silica"), ("purif_etoac", "m_etoac"), ("purif_pentane", "m_pentane")):
            if pid not in pids:
                r["participants"].append(st.new_participant(pid, pmid, "purification_agent", ev, stage="purification"))
        for s in r["steps"]:
            if s["stage"] == "purification":
                s["operation"] = F(
                    f"Column chromatography on silica gel (eluent: {info['eluent']})",
                    f"silica-gel column chromatography; {info['eluent']}", ev)
                s["participant_ids"] = ["purif_silica", "purif_etoac", "purif_pentane"]
                add_evidence(s["temperature"], ev)
                add_evidence(s["time"], ev)
        r["notes"].append(F("as a liquid", f"{label} was reported as a liquid", ev))


# Resolve the missing method field for the already-recorded 3u E/Z observation without inventing a ratio.
r3u = reactions[f"{DOC}:table3_3u_method_a"]
for outcome in r3u["outcomes"]:
    if outcome["outcome_id"] == "geometry":
        outcome["measurement_basis"] = F(
            "Table 3 annotation '(E:Z)'; characterization NMR/heading on page 8",
            "source annotation and characterization; numerical ratio and analytical method not reported", [f"{DOC}:page4:image", P[8], PI[8]])
        outcome["scope"] = F(
            "E:Z", {"metric": "E/Z", "order": "E:Z", "ratio": None, "named_component": "Z/Z"},
            [f"{DOC}:page4:image", P[8], PI[8]])
r3u["notes"].append(F(
    "Table 3 labels 3u '(E:Z)', while the page-8 characterization heading names the (Z)/(Z) form.",
    "E/Z mixture retained; the source does not report a ratio or fully assign both components.",
    [f"{DOC}:page4:image", P[8], PI[8]]))
r3u["review"]["unresolved"] = list(dict.fromkeys(r3u["review"]["unresolved"] + [
    "3u is marked E:Z in Table 3, but no numerical ratio or analytical method is reported; the characterization heading names only a Z/Z component."
]))


def add_unresolved_molecule(mid, label, note):
    ev = [P[10], PI[10]]
    m = st.new_molecule(mid, DOC, ev)
    m["occurrences"][0]["compound_label"] = F(label, label, ev)
    m["names"]["original"] = F(label, label, ev)
    m["names"]["english"] = F(label, label, ev)
    m["identity_note"] = F(note, note, ev)
    m["stereochemistry"]["chirality_status"] = F(
        "structure not specified", "unresolved", ev, method="inference",
        derivation="Footnote 19 provides a class label only; absence of stereo descriptors is not achirality.")
    m["stereochemistry"]["geometry_status"] = F(
        "structure not specified", "unresolved", ev, method="inference",
        derivation="Footnote 19 does not identify the substrate structure or alkene/allene configuration.")
    m["stereochemistry"]["representation_status"] = F(
        "structure not specified", "unresolved", ev, method="inference",
        derivation="No concrete structure can be encoded from the footnote alone.")
    bundle["molecules"].append(m)
    molecules[mid] = m
    return m


add_unresolved_molecule(
    "m_footnote19_enallene", "enallene in footnote 19 competition experiment",
    "Identity and substitution pattern are not stated on page 10; details are referred to unavailable SI.")
add_unresolved_molecule(
    "m_footnote19_allenyne", "allenyne in footnote 19 competition experiment",
    "Identity and substitution pattern are not stated on page 10; details are referred to unavailable SI.")
add_unresolved_molecule(
    "m_footnote19_pd_catalyst", "Pd catalyst in footnote 19 competition experiment",
    "The main-text footnote states Pd-catalyzed coupling but does not identify the palladium compound.")

# Reuse the existing B2pin2 identity but add the page-10 occurrence/evidence; do not infer a structure.
b2 = molecules["m_b2pin2"]
for occurrence in b2["occurrences"]:
    if occurrence["document_id"] == DOC:
        occurrence["evidence_ids"] = list(dict.fromkeys(occurrence["evidence_ids"] + [P[10], PI[10]]))
        add_evidence(occurrence["compound_label"], [P[10], PI[10]])
add_evidence(b2["names"]["original"], [P[10], PI[10]])
add_evidence(b2["names"]["english"], [P[10], PI[10]])
b2["identity_note"] = F(
    "B2pin2 is named in footnote 19; no source structure is drawn on page 10.",
    "Named reagent retained without a source-verified structure.", [P[10], PI[10]])

# Phenylboronic acid is explicitly named by formula; the connectivity is a machine candidate.
ev10 = [P[10], PI[10]]
phb = st.new_molecule("m_phboronic_acid", DOC, ev10)
phb["occurrences"][0]["compound_label"] = F("PhB(OH)2", "PhB(OH)2", ev10)
phb["names"]["original"] = F("PhB(OH)2", "phenylboronic acid", ev10)
phb["names"]["english"] = F("PhB(OH)2", "phenylboronic acid", ev10)
phb["identity_status"] = "source_identified"
phb["representation_type"] = "smiles"
phb["smiles"] = F(
    "PhB(OH)2", "OB(O)c1ccccc1", ev10, method="structure_derivation",
    derivation="Connectivity candidate expanded from the explicit formula/name; RDKit parsing is not source verification.")
stereo = st.stereo_inspect("OB(O)c1ccccc1")
phb["stereochemistry"]["chirality_status"] = F(
    "No stereogenic element in the derived graph", "achiral", ev10, method="structure_derivation",
    derivation="RDKit standard stereo inspection found no atom, axial, or other potential stereogenic element in the candidate graph.")
phb["stereochemistry"]["geometry_status"] = F(
    "No applicable E/Z element in the derived graph", "not_applicable", ev10, method="structure_derivation",
    derivation="RDKit standard stereo inspection found no stereogenic double bond in the candidate graph.")
phb["stereochemistry"]["representation_status"] = F(
    "candidate graph encoded", "encoded", ev10, method="structure_derivation",
    derivation="The candidate connectivity has no applicable stereochemical descriptor.")
phb["stereochemistry"]["isomeric_smiles"] = F(
    "PhB(OH)2", stereo["canonical_isomeric_smiles"], ev10, method="structure_derivation",
    derivation="Canonical isomeric SMILES from standard stereo inspection of the candidate graph.")
phb["stereochemistry"]["identity_key"] = F(
    None, stereo["stereo_graph_key"], ev10, method="structure_derivation",
    derivation="Stereo graph key from standard stereo inspection; not a sample identity claim.")
phb["identity_note"] = F(
    "PhB(OH)2 is explicitly named; candidate connectivity is not source-verified.",
    "Phenylboronic acid name is source-derived; graph remains a machine candidate.", ev10)
bundle["molecules"].append(phb)
molecules[phb["molecule_id"]] = phb


def competition_record(suffix, reagent_mid, reagent_label):
    ev = [P[10], PI[10]]
    r = st.new_reaction(f"{DOC}:{suffix}", DOC, ev)
    r["record_type"] = "specific_experiment"
    r["experiment_categories"] = F(["mechanistic_study"], ["mechanistic_study"], ev)
    r["location"] = F("reference footnote 19 on PDF page 10", "reference footnote 19", ev)
    r["reaction_name"] = F(
        "Pd-catalyzed competition in activation of an enallene and an allenyne",
        "Pd-catalyzed enallene/allenyne competition experiment", ev)
    r["transformation_type"] = F(
        f"coupling of competing substrates with {reagent_label}",
        "mechanistic competition coupling", ev)
    r["participants"] = [
        st.new_participant("enallene", "m_footnote19_enallene", "reactant", ev),
        st.new_participant("allenyne", "m_footnote19_allenyne", "reactant", ev),
        st.new_participant("coupling_reagent", reagent_mid, "reagent", ev),
        st.new_participant("pd_catalyst", "m_footnote19_pd_catalyst", "catalyst", ev),
    ]
    for participant in r["participants"][:2]:
        participant["charges"]["equivalents"] = F("1:1 ratio", q(1, "equiv"), ev)
    r["participants"][0]["equivalence_reference"] = F(
        "enallene:allenyne = 1:1", "m_footnote19_allenyne", ev)
    r["participants"][1]["equivalence_reference"] = F(
        "enallene:allenyne = 1:1", "m_footnote19_enallene", ev)
    s = st.new_step("competition_coupling", ev, stage="reaction")
    s["participant_ids"] = [p["participant_id"] for p in r["participants"]]
    s["operation"] = F(
        f"Pd-catalyzed coupling of competing enallene and allenyne with {reagent_label}",
        "allow a 1:1 enallene/allenyne mixture to compete under Pd-catalyzed coupling", ev)
    r["steps"] = [s]
    o = st.new_outcome("relative_reactivity", "selectivity", ["m_footnote19_enallene", "m_footnote19_allenyne"], ev)
    o["value"] = F(
        "the enallene reacts faster than the allenyne",
        {"order": ["m_footnote19_enallene", "m_footnote19_allenyne"], "relation": "faster_than", "ratio": None}, ev)
    o["measurement_basis"] = F(
        f"Pd-catalyzed coupling of the competing substrates with {reagent_label}",
        f"competition coupling with {reagent_label}; numerical method/result not reported in the main text", ev)
    o["scope"] = F(
        "activation by palladium(II)",
        {"metric": "relative activation/reactivity", "order": ["enallene", "allenyne"], "result": "enallene faster"}, ev)
    r["outcomes"] = [o]
    r["classification"]["bond_changes"] = []
    r["classification"]["labels"] = F(
        [], [], ev, method="inference",
        derivation="The footnote does not identify products or net bond changes, so no formation label is assigned.")
    r["classification"]["coverage"] = F(
        "products/net bond changes not reported", "unresolved", ev, method="inference",
        derivation="The mechanistic outcome is relative substrate reactivity, not an identified product transformation.")
    r["classification"]["limitations"] = [
        "Products and net bond changes are not reported in the supplied main-text page.",
        "Detailed conditions/results are referred to Supporting Information that was not supplied or associated.",
    ]
    r["classification"]["model_votes"] = []
    r["procedure"]["original_text"] = F(
        "These studies involved Pd-catalyzed coupling of the competing substrates with B2pin2 or PhB(OH)2 (see Supporting Information).",
        f"Footnote-level description of the {reagent_label} competition experiment", ev)
    r["review"] = {
        "status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
        "unresolved": [
            "Substrate identities, palladium catalyst identity, product identities, temperature, time, scale, solvent, and numerical relative-rate/selectivity result are not reported on page 10.",
            "Supporting Information is referenced by the paper but was not supplied or explicitly associated in this workspace.",
        ],
    }
    r["completeness"] = "partial"
    bundle["reactions"].append(r)
    reactions[r["reaction_id"]] = r
    return r["reaction_id"]


FOOTNOTE_B2 = competition_record("footnote19_competition_b2pin2", "m_b2pin2", "B2pin2")
FOOTNOTE_PHB = competition_record("footnote19_competition_phboronic_acid", "m_phboronic_acid", "PhB(OH)2")


def scope_ids(page):
    ids = []
    for label, info in PRODUCTS.items():
        if info["page"] != page:
            continue
        for method in info["methods"]:
            ids.append(reaction_id(label, method))
    return ids


bundle["source_context"].extend([
    {
        "context_id": f"{DOC}:page7_product_characterization",
        "kind": "characterization",
        "evidence_ids": [P[7], PI[7]],
        "description": "Continuation of 3c and NMR/HRMS/physical-form characterization blocks for 3d-3m. The isolated masses/yields and chromatography details were merged into the existing Table 2/3 experiment records; spectral peaks were not converted into synthesis records.",
        "molecule_ids": ["m_" + x for x in PRODUCTS if PRODUCTS[x]["page"] == 7],
    },
    {
        "context_id": f"{DOC}:page8_product_characterization",
        "kind": "characterization",
        "evidence_ids": [P[8], PI[8]],
        "description": "NMR/HRMS/physical-form characterization blocks for 3n-3v. The isolated masses/yields and chromatography details were merged into existing Table 3 experiment records; spectral peaks were not converted into synthesis records.",
        "molecule_ids": ["m_" + x for x in PRODUCTS if PRODUCTS[x]["page"] == 8],
    },
    {
        "context_id": f"{DOC}:page8_si_notice",
        "kind": "supplement_notice",
        "evidence_ids": [P[8], PI[8]],
        "description": "The article states that Supporting Information contains starting-material preparation, mechanistic/KIE details, and NMR spectra. No SI file was supplied or explicitly associated; no SI evidence was imported.",
    },
    {
        "context_id": f"{DOC}:page9_back_matter_references",
        "kind": "bibliography_and_back_matter",
        "evidence_ids": [P[9], PI[9]],
        "description": "Author contact continuation, competing-interest note, acknowledgments, and references 1-17. Bibliographic citations were not treated as current-paper synthesis experiments.",
    },
    {
        "context_id": f"{DOC}:page10_references",
        "kind": "bibliography",
        "evidence_ids": [P[10], PI[10]],
        "description": "References 17-22. Footnote 19 contains a current-work competition-study statement and is represented separately by two partial mechanistic-study records; the remaining citations are bibliography, not current experiments.",
    },
])


coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 7,
            "state": "attempted",
            "regions": [
                {
                    "region": "two-column continuation of product 3c and characterization blocks 3d-3m",
                    "finding": "No new reaction arrow. NMR/HRMS blocks are characterization; isolated masses/yields, physical form, and silica-gel eluent ratios repeat/enrich existing Table 2/3 Method A/B experiments.",
                    "evidence_ids": [P[7], PI[7]],
                }
            ],
            "record_ids": scope_ids(7),
            "unresolved": [
                "Complete product graphs and isomeric SMILES were not independently established from the characterization headings.",
                "Spectral assignments were retained as source context rather than interpreted as new synthesis experiments.",
            ],
        },
        {
            "page": 8,
            "state": "attempted",
            "regions": [
                {
                    "region": "product characterization blocks 3n-3v",
                    "finding": "No new reaction arrow. NMR/HRMS blocks are characterization; isolated masses/yields and chromatography details enrich existing Table 3 records. Product 3u remains an E/Z mixture with no numerical ratio, while its heading names a Z/Z component.",
                    "evidence_ids": [P[8], PI[8]],
                },
                {
                    "region": "Associated Content / Supporting Information notice",
                    "finding": "The page advertises an SI PDF, but no SI file was supplied or explicitly associated in this workspace; no SI facts were inferred.",
                    "evidence_ids": [P[8], PI[8]],
                },
                {
                    "region": "Author Information",
                    "finding": "Administrative metadata, not a reaction experiment.",
                    "evidence_ids": [P[8], PI[8]],
                },
            ],
            "record_ids": scope_ids(8),
            "unresolved": [
                "3u E/Z ratio and analytical method are unreported on the supplied page.",
                "Complete product graphs and isomeric SMILES were not independently established.",
            ],
        },
        {
            "page": 9,
            "state": "attempted",
            "regions": [
                {
                    "region": "author contact, Notes, Acknowledgments, and References 1-17",
                    "finding": "No current-paper reaction arrow, table, or experimental procedure. References are bibliographic citations and were not promoted to experiments.",
                    "evidence_ids": [P[9], PI[9]],
                }
            ],
            "record_ids": [],
            "unresolved": [],
        },
        {
            "page": 10,
            "state": "attempted",
            "regions": [
                {
                    "region": "References 17-22",
                    "finding": "Bibliographic citations are prior work, not current-paper experiments.",
                    "evidence_ids": [P[10], PI[10]],
                },
                {
                    "region": "reference footnote 19",
                    "finding": "Current-work mechanistic competition statement: a 1:1 enallene/allenyne mixture was tested in Pd-catalyzed coupling with B2pin2 or PhB(OH)2, and the enallene reacted faster. Two reagent-distinct partial mechanistic-study records were created.",
                    "evidence_ids": [P[10], PI[10]],
                },
            ],
            "record_ids": [FOOTNOTE_B2, FOOTNOTE_PHB],
            "unresolved": [
                "Footnote 19 omits substrate identities, catalyst identity, products, conditions, and numerical relative-rate/selectivity data; the cited SI was not supplied.",
            ],
        },
    ],
}

raw_extract = {
    "document_id": DOC,
    "bounded_pages": [7, 8, 9, 10],
    "semantic_accuracy_verified": False,
    "visual_review": "All four supplied page images were opened at original detail.",
    "product_characterization": PRODUCTS,
    "new_reaction_ids": [FOOTNOTE_B2, FOOTNOTE_PHB],
    "existing_reaction_ids_enriched": scope_ids(7) + scope_ids(8),
    "decisions": [
        "NMR/HRMS peak lists are source context, not synthesis experiments.",
        "Repeated isolated yield/mass/purification statements were merged into existing scope experiments.",
        "Footnote 19 was split by the explicitly alternative coupling reagent; missing SI details were not inferred.",
        "No bibliography entry was promoted to a current-paper experiment.",
    ],
}

(RUN / "coverage/pages_0007_0010.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/raw_extraction_pages_0007_0010.json").write_text(json.dumps(raw_extract, ensure_ascii=False, indent=2) + "\n")
(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({
    "molecules": len(bundle["molecules"]),
    "reactions": len(bundle["reactions"]),
    "source_context": len(bundle["source_context"]),
    "new_reactions": [FOOTNOTE_B2, FOOTNOTE_PHB],
    "enriched_existing_reactions": len(scope_ids(7) + scope_ids(8)),
}, ensure_ascii=False))
