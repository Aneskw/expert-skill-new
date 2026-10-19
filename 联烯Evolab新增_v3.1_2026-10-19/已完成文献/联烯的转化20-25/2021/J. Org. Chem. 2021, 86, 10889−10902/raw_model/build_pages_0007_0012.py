#!/usr/bin/env python3
"""Extend the v3.1 extraction with visually reviewed PDF pages 7--12.

The input is the already merged page-1--6 standard bundle. Existing scope and
derivatization experiments are patched with characterization evidence rather
than duplicated. New records are limited to preparations first encountered on
pages 7--9. All records remain machine-produced and needs_review.
"""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as sr  # noqa: E402

RUN = ROOT / "runs/current"
DOC = "doc_036d1c8f4ae80aff"
bundle = json.loads((RUN / "standard_records.json").read_text())
# The earlier page-1--6 merge predates rendered image evidence for later pages.
# Refresh document/evidence metadata from this batch's required standard-init
# template while preserving all prior chemistry records.
init_bundle = json.loads((RUN / "raw_model/pages_0007_0012_init.json").read_text())
bundle["documents"] = copy.deepcopy(init_bundle["documents"])


def ev(page: int) -> list[str]:
    return [f"{DOC}:page{page}:image", f"{DOC}:page{page}:text"]


def f(raw=None, normalized=None, evidence=(), *, method="original_extraction",
      status=None, verification=None, derivation=None, alternatives=()):
    if normalized is None and raw is not None:
        normalized = raw
    return sr.fact(raw=raw, normalized=normalized, evidence=evidence, method=method,
                   status=status, verification=verification, derivation=derivation,
                   alternatives=alternatives)


def inferred(normalized, evidence, derivation, raw=None):
    return f(raw, normalized, evidence, method="inference", status="pending_review",
             verification="inferred", derivation=derivation)


def q(number, unit, raw=None, qualifier=None):
    out = {"number": number, "unit": unit, "raw_text": raw or f"{number} {unit}"}
    if qualifier:
        out["qualifier"] = qualifier
    return out


def add_evidence_to_fact(fact, evidence):
    if not isinstance(fact, dict) or "evidence_ids" not in fact:
        return
    fact["evidence_ids"] = list(dict.fromkeys(fact["evidence_ids"] + list(evidence)))


def add_evidence_recursive(obj, evidence, only_if_page_text=None):
    if isinstance(obj, dict):
        if "evidence_ids" in obj:
            if only_if_page_text is None or only_if_page_text in obj["evidence_ids"]:
                add_evidence_to_fact(obj, evidence)
        for value in obj.values():
            add_evidence_recursive(value, evidence, only_if_page_text)
    elif isinstance(obj, list):
        for value in obj:
            add_evidence_recursive(value, evidence, only_if_page_text)


molecules_by_id = {m["molecule_id"]: m for m in bundle["molecules"]}
reactions_by_id = {r["reaction_id"]: r for r in bundle["reactions"]}


def label_of(m):
    return sr.value(m["occurrences"][0]["compound_label"]) if m.get("occurrences") else None


labels = {}
for _m in bundle["molecules"]:
    lab = label_of(_m)
    if lab and lab not in labels:
        labels[lab] = _m["molecule_id"]


SMILES = {
    # Page-7 propargyl anthranilate derivatives and their literature precursors.
    "S1h": "COC(=O)c1c(NS(=O)(=O)c2ccc(C)cc2)cc(C(F)(F)F)cc1",
    "S1i": "COC(=O)c1c(NS(=O)(=O)c2ccc(C)cc2)ccc(OC)c1",
    "S1j": "COC(=O)c1c(NS(=O)(=O)c2ccc(C)cc2)ccc(Br)c1",
    "S1n": "COC(=O)c1c(NS(=O)(=O)c2ccc(C)cc2)cc(OC)c(OC)c1",
    "S10h": "C#CCN(c1c(C(=O)OC)cc(C(F)(F)F)cc1)S(=O)(=O)c1ccc(C)cc1",
    "S10i": "C#CCN(c1c(C(=O)OC)ccc(OC)c1)S(=O)(=O)c1ccc(C)cc1",
    "S10j": "C#CCN(c1c(C(=O)OC)ccc(Br)c1)S(=O)(=O)c1ccc(C)cc1",
    "S10n": "C#CCN(c1c(C(=O)OC)cc(OC)c(OC)c1)S(=O)(=O)c1ccc(C)cc1",
    # Homologated allenes.
    "S11d": "C=C=CCN(c1ccccc1C(=O)OCC)C(=O)OC(C)(C)C",
    "S11e": "C=C=CCN(c1ccccc1C(=O)OCC)S(=O)(=O)c1ccc(C)cc1",
    "S11f": "C=C=CCN(c1c(C(=O)OC)cc(C)cc1)S(=O)(=O)c1ccc(C)cc1",
    "S11h": "C=C=CCN(c1c(C(=O)OC)cc(C(F)(F)F)cc1)S(=O)(=O)c1ccc(C)cc1",
    "S11i": "C=C=CCN(c1c(C(=O)OC)ccc(OC)c1)S(=O)(=O)c1ccc(C)cc1",
    "S11j": "C=C=CCN(c1c(C(=O)OC)ccc(Br)c1)S(=O)(=O)c1ccc(C)cc1",
    "S11n": "C=C=CCN(c1c(C(=O)OC)cc(OC)c(OC)c1)S(=O)(=O)c1ccc(C)cc1",
    # Ester reductions.
    "S12d": "C=C=CCN(c1ccccc1CO)C(=O)OC(C)(C)C",
    "S12f": "C=C=CCN(c1c(CO)ccc(C)c1)S(=O)(=O)c1ccc(C)cc1",
    "S12h": "C=C=CCN(c1c(CO)ccc(C(F)(F)F)c1)S(=O)(=O)c1ccc(C)cc1",
    "S12i": "C=C=CCN(c1c(CO)cc(OC)cc1)S(=O)(=O)c1ccc(C)cc1",
    "S12j": "C=C=CCN(c1c(CO)cc(Br)cc1)S(=O)(=O)c1ccc(C)cc1",
    "S12n": "C=C=CCN(c1c(CO)cc(OC)c(OC)c1)S(=O)(=O)c1ccc(C)cc1",
    # Scope substrates now named explicitly in the Experimental Section.
    "3f": "C=C=CCN(c1c(CNc2ccc(OC)cc2)ccc(C)c1)S(=O)(=O)c1ccc(C)cc1",
    "3g": "C=C=CCN(c1c(CNc2ccc(OC)cc2)ccc(Cl)c1)S(=O)(=O)c1ccc(C)cc1",
    "3h": "C=C=CCN(c1c(CNc2ccc(OC)cc2)ccc(C(F)(F)F)c1)S(=O)(=O)c1ccc(C)cc1",
    "3i": "C=C=CCN(c1c(CNc2ccc(OC)cc2)cc(OC)cc1)S(=O)(=O)c1ccc(C)cc1",
    "3j": "C=C=CCN(c1c(CNc2ccc(OC)cc2)cc(Br)cc1)S(=O)(=O)c1ccc(C)cc1",
    "3k": "C=C=CCN(c1c(CNc2ccc(OC)cc2)cc(Cl)cc1)S(=O)(=O)c1ccc(C)cc1",
    "3l": "C=C=CCN(c1c(CNc2ccc(OC)cc2)c(F)ccc1)S(=O)(=O)c1ccc(C)cc1",
    "3m": "C=C=CCN(c1c(CNc2ccc(OC)cc2)c(C)ccc1)S(=O)(=O)c1ccc(C)cc1",
    "3n": "C=C=CCN(c1c(CNc2ccc(OC)cc2)cc(OC)c(OC)c1)S(=O)(=O)c1ccc(C)cc1",
    # Alternative S13/S14/S15 route for g/k/l/m.
    "S13g": "O=S(=O)(Nc1c(CO)ccc(Cl)c1)c1ccc(C)cc1",
    "S13k": "O=S(=O)(Nc1c(CO)cc(Cl)cc1)c1ccc(C)cc1",
    "S13l": "O=S(=O)(Nc1c(CO)c(F)ccc1)c1ccc(C)cc1",
    "S13m": "O=S(=O)(Nc1c(CO)c(C)ccc1)c1ccc(C)cc1",
    "S14g": "C#CCN(c1c(CO)ccc(Cl)c1)S(=O)(=O)c1ccc(C)cc1",
    "S14k": "C#CCN(c1c(CO)cc(Cl)cc1)S(=O)(=O)c1ccc(C)cc1",
    "S14l": "C#CCN(c1c(CO)c(F)ccc1)S(=O)(=O)c1ccc(C)cc1",
    "S14m": "C#CCN(c1c(CO)c(C)ccc1)S(=O)(=O)c1ccc(C)cc1",
    "S15g": "C=C=CCN(c1c(CO)ccc(Cl)c1)S(=O)(=O)c1ccc(C)cc1",
    "S15k": "C=C=CCN(c1c(CO)cc(Cl)cc1)S(=O)(=O)c1ccc(C)cc1",
    "S15l": "C=C=CCN(c1c(CO)c(F)ccc1)S(=O)(=O)c1ccc(C)cc1",
    "S15m": "C=C=CCN(c1c(CO)c(C)ccc1)S(=O)(=O)c1ccc(C)cc1",
    # Reagents and work-up materials first encountered in this page range.
    "CuI": "[Cu]I",
    "dicyclohexylamine": "C1CCC(CC1)NC1CCCCC1",
    "1,4-dioxane": "O1CCOCC1",
    "chloroform": "ClC(Cl)Cl",
    "ammonium hydroxide": "[NH4+].[OH-]",
    "Rochelle salt": "[Na+].[K+].[O-]C(=O)[C@@H](O)[C@H](O)C(=O)[O-]",
    "NaOH": "[Na+].[OH-]",
    "Na2SO4": "[Na+].[Na+].[O-]S(=O)(=O)[O-]",
    "brine": "[Na+].[Cl-]",
    "EtOAc": "CCOC(C)=O",
}

FULL_NAMES = {
    "S10h": "Methyl 2-((4-methyl-N-(prop-2-yn-1-yl)phenyl)sulfonamido)-4-(trifluoromethyl)benzoate",
    "S10i": "Methyl 5-methoxy-2-((4-methyl-N-(prop-2-yn-1-yl)phenyl)sulfonamido)benzoate",
    "S10j": "Methyl 5-bromo-2-((4-methyl-N-(prop-2-yn-1-yl)phenyl)sulfonamido)benzoate",
    "S10n": "Methyl 4,5-dimethoxy-2-((4-methyl-N-(prop-2-yn-1-yl)phenyl)sulfonamido)benzoate",
    "S11d": "Ethyl 2-(buta-2,3-dien-1-yl(tert-butoxycarbonyl)amino)benzoate",
    "S11e": "Ethyl 2-((N-(buta-2,3-dien-1-yl)-4-methylphenyl)sulfonamido)benzoate",
    "S11f": "Methyl 2-((N-(buta-2,3-dien-1-yl)-4-methylphenyl)sulfonamido)-4-methylbenzoate",
    "S11h": "Methyl 2-((N-(buta-2,3-dien-1-yl)-4-methylphenyl)sulfonamido)-4-(trifluoromethyl)benzoate",
    "S11i": "Methyl 2-((N-(buta-2,3-dien-1-yl)-4-methylphenyl)sulfonamido)-5-methoxybenzoate",
    "S11j": "Methyl 5-bromo-2-((N-(buta-2,3-dien-1-yl)-4-methylphenyl)sulfonamido)benzoate",
    "S11n": "Methyl 2-((N-(buta-2,3-dien-1-yl)-4-methylphenyl)sulfonamido)-4,5-dimethoxybenzoate",
}


def molecule(label, page, *, full_name=None, smiles=None, note=None):
    evidence = ev(page)
    if label in labels:
        m = molecules_by_id[labels[label]]
        if not any(sr.value(o["compound_label"]) == label and set(evidence) <= set(o["evidence_ids"])
                   for o in m["occurrences"]):
            m["occurrences"].append({"document_id": DOC,
                                     "compound_label": f(label, label, evidence),
                                     "evidence_ids": evidence})
        if full_name:
            m["names"]["english"] = f(full_name, full_name, evidence)
        smi = smiles or SMILES.get(label)
        if smi and sr.value(m["smiles"]) is None:
            apply_smiles(m, label, smi, evidence)
        return m["molecule_id"]
    mid = f"{DOC}:mol_{label.replace(' ', '_').replace('/', '_').replace(',', '')}"
    m = sr.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    m["names"]["original"] = f(label, label, evidence)
    if full_name:
        m["names"]["english"] = f(full_name, full_name, evidence)
    smi = smiles or SMILES.get(label)
    if smi:
        apply_smiles(m, label, smi, evidence)
    else:
        m["identity_note"] = f(note or "Source identity retained without a concrete structure.",
                               note or "Source identity retained without a concrete structure.", evidence)
        m["stereochemistry"]["chirality_status"] = inferred(
            "unresolved", evidence, "No concrete source-derived graph is available for a stereo audit.")
        m["stereochemistry"]["geometry_status"] = inferred(
            "unresolved", evidence, "No concrete source-derived graph is available for an E/Z or ring-cis/trans audit.")
        m["stereochemistry"]["representation_status"] = inferred(
            "unresolved", evidence, "No concrete structure was encoded.")
    bundle["molecules"].append(m)
    molecules_by_id[mid] = m
    labels[label] = mid
    return mid


def apply_smiles(m, label, smiles, evidence):
    diag = sr.stereo_inspect(smiles)
    canon = diag["canonical_isomeric_smiles"]
    potential = diag["potential_stereo"]
    if any(p["specified"] == "Unspecified" for p in potential):
        chirality = "unspecified" if any(p["type"].startswith("Atom_") for p in potential) else "achiral"
        geometry = "unspecified" if any(p["type"].startswith("Bond_") for p in potential) else "not_applicable"
        representation = "partial"
    else:
        chirality, geometry, representation = "achiral", "not_applicable", "encoded"
    m["representation_type"] = "smiles"
    m["smiles"] = f(label, canon, evidence, method="structure_derivation",
                      verification="inferred",
                      derivation="Candidate transcribed from the source systematic name; parsing is not source verification.")
    m["stereochemistry"]["isomeric_smiles"] = f(
        smiles, canon, evidence, method="structure_derivation", verification="inferred",
        derivation="Canonical isomeric SMILES generated by standard stereo diagnostics.")
    m["stereochemistry"]["identity_key"] = f(
        canon, diag["stereo_graph_key"], evidence, method="structure_derivation",
        verification="inferred", derivation="Standard stereo graph key; not a sample identity.")
    m["stereochemistry"]["chirality_status"] = inferred(
        chirality, evidence, "Source name and standard stereo graph were inspected; terminal allenes have a CH2 end and therefore no allene axis assignment.")
    m["stereochemistry"]["geometry_status"] = inferred(
        geometry, evidence, "No E/Z-eligible alkene or ring cis/trans assignment is present in the source-derived graph.")
    m["stereochemistry"]["representation_status"] = inferred(
        representation, evidence, "Compared with standard stereo diagnostics; source verification remains pending.")
    if "C=C=C" in canon:
        m["stereochemistry"]["axial"] = inferred(
            {"status": "not_applicable", "reason": "terminal allene end is CH2"}, evidence,
            "The systematic name and graph show a terminal CH2=C= unit, whose terminal carbon bears two identical hydrogens.")
    m["identity_note"] = inferred(
        "Candidate structure derived from the printed systematic name; expert source verification pending.",
        evidence, "No database identity or CAS number was inferred.")


def agent(label, page=7, smiles=None, note=None):
    return molecule(label, page, smiles=smiles, note=note)


def participant(pid, mid, role, evidence, *, eq=None, molpct=None, mass=None,
                amount=None, volume=None, ref=None, stage="reaction"):
    p = sr.new_participant(pid, mid, role, evidence, stage)
    for key, value in (("equivalents", eq), ("mol_percent", molpct), ("mass", mass),
                       ("amount", amount), ("volume", volume)):
        if value is not None:
            raw = value.get("raw_text") if isinstance(value, dict) else str(value)
            p["charges"][key] = f(raw, value, evidence)
    if ref:
        p["equivalence_reference"] = f(ref, ref, evidence)
    return p


def step(sid, evidence, operation, *, temperature=None, time=None, atmosphere=None,
         pids=(), stage="reaction"):
    s = sr.new_step(sid, evidence, stage)
    s["participant_ids"] = list(pids)
    s["operation"] = f(operation, operation, evidence)
    if temperature is not None:
        raw = temperature.get("raw_text") if isinstance(temperature, dict) else str(temperature)
        s["temperature"] = f(raw, temperature, evidence)
    if time is not None:
        raw = time.get("raw_text") if isinstance(time, dict) else str(time)
        s["time"] = f(raw, time, evidence)
    if atmosphere is not None:
        s["atmosphere"] = f(atmosphere, atmosphere, evidence)
    return s


def classify(r, changes, evidence, limitation):
    labels_out = {}
    for idx, (a, b, before, after, ctype, hcount) in enumerate(changes, 1):
        cid = f"b{idx}"
        ch = {"change_id": cid, "atoms": [{"atom_map": None, "element": a},
                                            {"atom_map": None, "element": b}],
              "before_order": before, "after_order": after, "change_type": ctype}
        if hcount is not None:
            ch["hydrogen_count"] = hcount
        r["classification"]["bond_changes"].append(inferred(
            ch, evidence, "Net change inferred from the named source and transformation; atom mapping was not asserted.",
            raw=f"{ctype} {a}-{b}"))
        if ctype in {"formed", "order_increased", "hydrogen_gain"}:
            lab = sr.bond_label(a, b, after)
            labels_out.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(cid)
    r["classification"]["labels"] = inferred(
        list(labels_out.values()), evidence, "Calculated from retained net changes under bond-formation-20261016.")
    r["classification"]["coverage"] = inferred(
        "partial", evidence, "No complete atom mapping or coproduct accounting was asserted.")
    r["classification"]["limitations"] = [limitation, "No atom-mapped reaction SMILES was asserted."]


def yield_outcome(target, raw, number, evidence, basis="isolated after silica gel chromatography"):
    o = sr.new_outcome("yield", "yield", [target], evidence)
    o["product_role"] = f("product", "product", evidence)
    o["value"] = f(raw, q(number, "%", raw), evidence)
    o["yield_type"] = f("isolated yield", "isolated yield", evidence)
    o["measurement_basis"] = f(basis, basis, evidence)
    return o


def make_reaction(suffix, page, name, transformation, participants, steps, product_mid,
                  yield_raw, yield_number, general_ids, changes, *, extra_evidence=(), notes=(),
                  unresolved=()):
    evidence = list(dict.fromkeys(ev(page) + list(extra_evidence)))
    r = sr.new_reaction(f"{DOC}:{suffix}", DOC, evidence)
    r["record_type"] = "specific_experiment"
    r["experiment_categories"] = f(["other"], ["other"], evidence)
    r["location"] = f(f"Page {page}, Experimental Section", f"Page {page}, Experimental Section", evidence)
    r["reaction_name"] = f(name, name, evidence)
    r["transformation_type"] = f(transformation, transformation, evidence)
    r["general_scheme_ids"] = list(general_ids)
    r["participants"] = participants
    r["steps"] = steps
    r["outcomes"] = [yield_outcome(product_mid, yield_raw, yield_number, evidence)]
    r["notes"] = [f(n, n, evidence) for n in notes]
    r["review"]["unresolved"] = list(unresolved) or [
        "Machine structure candidates and inferred bond changes require expert review."]
    classify(r, changes, evidence, "Only the principal named transformation is classified; work-up species and unshown coproducts are excluded.")
    bundle["reactions"].append(r)
    reactions_by_id[r["reaction_id"]] = r
    return r


def generic(suffix, page, name, transformation, participants, steps, notes=()):
    evidence = ev(page)
    r = sr.new_reaction(f"{DOC}:{suffix}", DOC, evidence)
    r["record_type"] = "generic_scheme"
    r["experiment_categories"] = f([], [], evidence)
    r["location"] = f(f"Page {page}, Experimental Section", f"Page {page}, Experimental Section", evidence)
    r["reaction_name"] = f(name, name, evidence)
    r["transformation_type"] = f(transformation, transformation, evidence)
    r["participants"] = participants
    r["steps"] = steps
    r["notes"] = [f(n, n, evidence) for n in notes]
    r["review"]["unresolved"] = ["Generic procedure retained separately; not counted as a specific experiment."]
    bundle["reactions"].append(r)
    reactions_by_id[r["reaction_id"]] = r
    return r["reaction_id"]


# Shared/reused agents. No CAS values are invented.
for label, smiles in {
    "NaH": "[Na+].[H-]", "propargyl bromide derivative": "C#CCBr", "DMF": "CN(C)C=O",
    "CuI": SMILES["CuI"], "(CHO)n": None, "dicyclohexylamine": SMILES["dicyclohexylamine"],
    "1,4-dioxane": SMILES["1,4-dioxane"], "chloroform": SMILES["chloroform"],
    "ammonium hydroxide": SMILES["ammonium hydroxide"], "water": "O", "Na2SO4": SMILES["Na2SO4"],
    "DIBAL-H": "CC(C)[AlH]C(C)C", "DCM": "ClCCl", "MeOH": "CO",
    "Rochelle salt": SMILES["Rochelle salt"], "brine": SMILES["brine"], "DMP": "CC(=O)OI1(OC(C)=O)OC(=O)c2ccccc21",
    "NaOH": SMILES["NaOH"], "p-anisidine": "COc1ccc(N)cc1", "NaBH4": "[Na+].[BH4-]",
    "K2CO3": "[K+].[K+].[O-]C([O-])=O", "NH4Cl": "[NH4+].[Cl-]", "EtOAc": SMILES["EtOAc"],
}.items():
    agent(label, 7 if label not in {"K2CO3", "NH4Cl", "EtOAc"} else 8, smiles=smiles,
          note="Paraformaldehyde is retained as the printed polymeric formula; no discrete-chain SMILES is asserted." if label == "(CHO)n" else None)


# New generic procedures, with all explicitly quantified reagents preserved.
hom_ev = ev(7)
generic_s10 = molecule("generic S10", 7, note="Generic terminal alkyne member of the S10 series; no concrete substituent is asserted.")
HOM = generic(
    "gp_homologation_alkynes_to_allenes", 7, "Homologation of alkynes to allenes",
    "terminal alkyne homologation to allene",
    [participant("generic_substrate", generic_s10, "reactant", hom_ev, eq=q(1, "equiv", "1 equiv"), ref=generic_s10),
     participant("CuI", labels["CuI"], "catalyst", hom_ev, eq=q(0.5, "equiv", "0.5 equiv"), ref=generic_s10),
     participant("paraformaldehyde", labels["(CHO)n"], "reagent", hom_ev, eq=q(2.5, "equiv", "2.5 equiv"), ref=generic_s10),
     participant("Cy2NH", labels["dicyclohexylamine"], "reagent", hom_ev, eq=q(1.8, "equiv", "1.8 equiv"), ref=generic_s10),
     participant("dioxane", labels["1,4-dioxane"], "solvent", hom_ev)],
    [step("reflux", hom_ev, "Add CuI, (CHO)n and Cy2NH to alkyne S10 in dioxane (0.2 M); heat to reflux", temperature="reflux", time=q(6, "h", "6 h")),
     step("cool_concentrate", hom_ev, "Cool to room temperature and remove solvent in vacuo", temperature="rt", stage="workup"),
     step("wash", hom_ev, "Dissolve in CHCl3 (50 mL); wash with 10% NH4OH (2 × 50 mL) and H2O (2 × 50 mL); dry Na2SO4", stage="workup"),
     step("purification", hom_ev, "Silica gel chromatography, hexanes/EtOAc 8:2", stage="purification")])

nalk2_ev = ev(8)
generic_s13 = molecule("generic S13", 8, note="Generic corresponding tosylamide S13; substituent identity varies by substrate.")
NALK2 = generic(
    "gp_S13_N_alkylation", 8, "N-alkylation of tosylamides S13", "N-propargylation",
    [participant("generic_substrate", generic_s13, "reactant", nalk2_ev, eq=q(1, "equiv", "1 equiv"), ref=generic_s13),
     participant("K2CO3", labels["K2CO3"], "reagent", nalk2_ev, eq=q(2, "equiv", "2 equiv"), ref=generic_s13),
     participant("propargyl_bromide", labels["propargyl bromide derivative"], "reagent", nalk2_ev, eq=q(1.5, "equiv", "1.5 equiv"), ref=generic_s13),
     participant("DMF", labels["DMF"], "solvent", nalk2_ev)],
    [step("premix", nalk2_ev, "Stir corresponding tosylamide S13 with K2CO3 in DMF (0.25 M)", temperature="rt", time=q(30, "min", "30 min"), atmosphere="argon"),
     step("alkylation", nalk2_ev, "Add propargyl bromide derivative and warm to 80 °C", temperature=q(80, "°C", "80 °C"), time=q(16, "h", "16 h"), atmosphere="argon"),
     step("workup", nalk2_ev, "Quench with saturated NH4Cl; extract EtOAc; wash saturated NH4Cl (3 × 100 mL); dry Na2SO4; concentrate", stage="workup"),
     step("purification", nalk2_ev, "Silica gel chromatography, hexanes/EtOAc 7:3", stage="purification")],
    notes=("S13 starting tosylamides are stated to be literature preparations; no preparation of S13 was created.",))

asym_ev = ev(9)
generic_allene3 = molecule("generic allene 3", 9, note="Generic member of the allene-3 substrate series; substituent identity varies by scope entry.")
ASYM_ALLENE = generic(
    "gp_asymmetric_allene_cyclization", 9, "General procedure for asymmetric cyclization of allenes 3",
    "enantioselective intramolecular allene hydrofunctionalization",
    [participant("generic_substrate", generic_allene3, "reactant", asym_ev, eq=q(1, "equiv", "1 equiv"), amount=q(0.15, "mmol", "0.15 mmol"), ref=generic_allene3),
     participant("Rh", agent("[Rh(cod)Cl]2", 9), "catalyst", asym_ev, eq=q(0.04, "equiv", "0.04 equiv"), mass=q(3, "mg", "3 mg"), amount=q(6, "μmol", "6 μmol"), ref=generic_allene3),
     participant("PPTS", agent("PPTS", 9), "additive", asym_ev, eq=q(0.1, "equiv", "0.1 equiv"), mass=q(4, "mg", "4 mg"), amount=q(15, "μmol", "15 μmol"), ref=generic_allene3),
     participant("chloroacetic_acid", agent("chloroacetic acid", 9), "additive", asym_ev, eq=q(0.1, "equiv", "0.1 equiv"), mass=q(1.4, "mg", "1.4 mg"), amount=q(15, "μmol", "15 μmol"), ref=generic_allene3),
     participant("ligand", agent("(R)-DTBM-Garphos", 9), "ligand", asym_ev, eq=q(0.1, "equiv", "0.1 equiv"), mass=q(19, "mg", "19 mg"), amount=q(15, "μmol", "15 μmol"), ref=generic_allene3),
     participant("DCE", agent("DCE", 9), "solvent", asym_ev, volume=q(0.4, "mL", "0.4 mL"))],
    [step("catalyst_setup", asym_ev, "Charge catalyst, one Brønsted acid and ligand; evacuate/backfill argon three times", atmosphere="argon"),
     step("preactivation", asym_ev, "Add DCE and stir", temperature="rt", time=q(10, "min", "10 min"), atmosphere="argon"),
     step("cyclization", asym_ev, "Add allene 3 (0.15 mmol, 1 equiv) under argon and heat", temperature=q(50, "°C", "50 °C"), time=q(24, "h", "24 h"), atmosphere="argon"),
     step("purification", asym_ev, "Cool; remove solvent; silica gel chromatography, hexanes/EtOAc 9:1", stage="purification")],
    notes=("PPTS and chloroacetic acid are alternatives selected per substrate, not co-reagents in one experiment.",))


# Page-7 N-alkylations S1 -> S10, continuing the page-4 general procedure.
old_nalk = f"{DOC}:gp_n_alkylation"
for label, page, y, mass, mmol in [
    ("S10h", 7, 84, 1.8, 4.5), ("S10i", 7, 83, 0.62, 1.7),
    ("S10j", 7, 70, 2.95, 7.0), ("S10n", 7, 99, 0.8, 1.98),
]:
    precursor = "S1" + label[-1]
    smid = molecule(precursor, page, smiles=SMILES[precursor], note="Literature-derived starting tosylamide; preparation not extracted as a current experiment.")
    pmid = molecule(label, page, full_name=FULL_NAMES[label])
    evidence = ev(page) + ev(4)
    parts = [participant("substrate", smid, "reactant", evidence, eq=q(1, "equiv", "1 equiv"), ref=smid),
             participant("product", pmid, "product", evidence, mass=q(mass, "g", f"{mass:g} g"), amount=q(mmol, "mmol", f"{mmol:g} mmol")),
             participant("NaH", labels["NaH"], "reagent", evidence, eq=q(1.2, "equiv", "1.2 equiv"), ref=smid),
             participant("propargyl_bromide", labels["propargyl bromide derivative"], "reagent", evidence, eq=q(1.3, "equiv", "1.3 equiv"), ref=smid),
             participant("DMF", labels["DMF"], "solvent", evidence)]
    steps = [step("deprotonation", evidence, "Add substrate solution slowly to NaH in DMF (0.33 M) and hold", temperature=q(0, "°C", "0 °C"), time=q(2, "h", "2 h"), atmosphere="argon"),
             step("alkylation", evidence, "Add propargyl bromide derivative; warm slowly to room temperature and stir", temperature="rt", time=q(16, "h", "16 h"), atmosphere="argon"),
             step("workup", evidence, "Quench with saturated NH4Cl; extract EtOAc; wash saturated NH4Cl (3 × 100 mL); dry Na2SO4; concentrate", stage="workup"),
             step("purification", evidence, "Silica gel chromatography, hexanes/EtOAc 8:2", stage="purification")]
    make_reaction(f"prep_{label}", page, f"N-alkylation to {label}", "N-propargylation", parts, steps, pmid,
                  f"{y}%", y, [old_nalk], [("C", "N", 0, 1, "formed", None)],
                  extra_evidence=ev(4), unresolved=("Starting tosylamide preparation is literature-derived.",))


# Page-7 homologations S10 -> S11.
for suffix, page, y, mass, mmol in [
    ("d", 7, 82, 2.44, 7.7), ("e", 7, 90, 3.0, 8.1), ("f", 7, 85, 2.84, 7.65),
    ("h", 7, 70, 1.3, 3.0), ("i", 7, 70, 0.45, 1.16), ("j", 7, 80, 2.44, 5.6),
    ("n", 7, 56, 0.46, 1.1),
]:
    smid = molecule("S10" + suffix, page)
    pmid = molecule("S11" + suffix, page, full_name=FULL_NAMES["S11" + suffix])
    evidence = ev(page)
    parts = [participant("substrate", smid, "reactant", evidence, eq=q(1, "equiv", "1 equiv"), ref=smid),
             participant("CuI", labels["CuI"], "catalyst", evidence, eq=q(0.5, "equiv", "0.5 equiv"), ref=smid),
             participant("paraformaldehyde", labels["(CHO)n"], "reagent", evidence, eq=q(2.5, "equiv", "2.5 equiv"), ref=smid),
             participant("Cy2NH", labels["dicyclohexylamine"], "reagent", evidence, eq=q(1.8, "equiv", "1.8 equiv"), ref=smid),
             participant("dioxane", labels["1,4-dioxane"], "solvent", evidence),
             participant("CHCl3", labels["chloroform"], "workup_agent", evidence, volume=q(50, "mL", "50 mL"), stage="workup"),
             participant("NH4OH", labels["ammonium hydroxide"], "workup_agent", evidence, volume=q(100, "mL", "2 × 50 mL"), stage="workup"),
             participant("water", labels["water"], "workup_agent", evidence, volume=q(100, "mL", "2 × 50 mL"), stage="workup"),
             participant("Na2SO4", labels["Na2SO4"], "workup_agent", evidence, stage="workup"),
             participant("product", pmid, "product", evidence, mass=q(mass, "g", f"{mass:g} g"), amount=q(mmol, "mmol", f"{mmol:g} mmol"))]
    steps = copy.deepcopy(reactions_by_id[HOM]["steps"])
    make_reaction(f"prep_S11{suffix}", page, f"Homologation of S10{suffix} to S11{suffix}",
                  "terminal alkyne homologation to allene", parts, steps, pmid, f"{y}%", y, [HOM],
                  [("C", "C", 0, 2, "formed", None)],
                  unresolved=("New allene carbon is inferred to derive from paraformaldehyde; no atom mapping was reported.",))


# Ester reductions S11 -> 3a/S12d/f/h/i/j/n (pages 7--8).
reduce_products = {"e": "3a", "d": "S12d", "f": "S12f", "h": "S12h", "i": "S12i", "j": "S12j", "n": "S12n"}
reduce_data = {
    "e": (7, 94, 2.5, 7.6), "d": (7, 88, None, None), "f": (7, 90, 2.37, 6.9),
    "h": (8, 64, 0.7, 1.86), "i": (8, 54, 0.23, 0.663), "j": (8, 90, 2.0, 5.0),
    "n": (8, 57, 0.24, 0.61),
}
REDUCE = f"{DOC}:gp_ester_reduction"
for suffix, product_label in reduce_products.items():
    page, y, mass, mmol = reduce_data[suffix]
    smid = molecule("S11" + suffix, page)
    pmid = molecule(product_label, page)
    evidence = list(dict.fromkeys(ev(page) + ev(4)))
    parts = [participant("substrate", smid, "reactant", evidence, eq=q(1, "equiv", "1 equiv"), ref=smid),
             participant("DIBAL-H", labels["DIBAL-H"], "reagent", evidence, eq=q(2.2, "equiv", "2.2 equiv"), ref=smid),
             participant("DCM_reaction", labels["DCM"], "solvent", evidence),
             participant("MeOH", labels["MeOH"], "workup_agent", evidence, volume=q(5, "mL", "5 mL"), stage="workup"),
             participant("Rochelle_salt", labels["Rochelle salt"], "workup_agent", evidence, stage="workup"),
             participant("DCM_extract", labels["DCM"], "workup_agent", evidence, volume=q(90, "mL", "3 × 30 mL"), stage="workup"),
             participant("brine", labels["brine"], "workup_agent", evidence, stage="workup"),
             participant("Na2SO4", labels["Na2SO4"], "workup_agent", evidence, stage="workup")]
    prod_kwargs = {}
    if mass is not None: prod_kwargs["mass"] = q(mass, "g", f"{mass:g} g")
    if mmol is not None: prod_kwargs["amount"] = q(mmol, "mmol", f"{mmol:g} mmol")
    parts.append(participant("product", pmid, "product", evidence, **prod_kwargs))
    steps = [step("low_temperature_addition", evidence, "Add 1 M DIBAL-H in DCM dropwise to ester in DCM (0.3 M)", temperature=q(-78, "°C", "−78 °C"), atmosphere="argon"),
             step("low_temperature_hold", evidence, "Stir at −78 °C", temperature=q(-78, "°C", "−78 °C"), time=q(3, "h", "3 h"), atmosphere="argon"),
             step("quench", evidence, "Add MeOH (5 mL), then saturated Rochelle salt at −78 °C", temperature=q(-78, "°C", "−78 °C"), stage="workup"),
             step("warm_workup", evidence, "Warm to room temperature and stir", temperature="rt", time=q(1, "h", "1 h"), stage="workup"),
             step("extraction", evidence, "Extract DCM (3 × 30 mL); wash saturated NaCl; dry Na2SO4; concentrate", stage="workup"),
             step("purification", evidence, "Silica gel chromatography, hexanes/EtOAc 9:1 to 7:3", stage="purification")]
    make_reaction(f"prep_{product_label}", page, f"Ester reduction of S11{suffix} to {product_label}",
                  "ester reduction to benzyl alcohol", parts, steps, pmid, f"{y}%", y, [REDUCE],
                  [("C", "H", 0, 1, "hydrogen_gain", 2)], extra_evidence=ev(4),
                  unresolved=(("Product mass and amount are not reported." if mass is None else "Complete atom mapping and removed ethoxy coproduct are unreported."),))


# Oxidation/reductive-amination sequences to 3d/e/f/h/i/j/n and 3g/k/l/m.
OXRA = f"{DOC}:gp_oxidation_reductive_amination"
ox_data = {
    "d": ("S12d", 8, 70, 1.8, 4.8), "e": ("3a", 8, 70, 1.5, 3.5),
    "f": ("S12f", 8, 80, 0.716, 1.6), "h": ("S12h", 8, 73, 0.7, 1.4),
    "i": ("S12i", 8, 50, 0.26, 0.56), "j": ("S12j", 8, 65, 1.28, 2.5),
    "n": ("S12n", 8, 64, 0.18, 0.36),
}


def oxra_record(suffix, source_label, page, y, mass_g, mmol):
    smid = molecule(source_label, page)
    pmid = molecule("3" + suffix, page)
    evidence = list(dict.fromkeys(ev(page) + ev(5)))
    parts = [participant("alcohol", smid, "reactant", evidence, eq=q(1, "equiv", "1 equiv"), ref=smid),
             participant("DMP", labels["DMP"], "reagent", evidence, eq=q(1.1, "equiv", "1.1 equiv"), ref=smid),
             participant("DCM_oxidation", labels["DCM"], "solvent", evidence),
             participant("NaOH", labels["NaOH"], "workup_agent", evidence, stage="workup"),
             participant("p_anisidine", labels["p-anisidine"], "reactant", evidence, eq=q(1.5, "equiv", "1.5 equiv"), ref=smid),
             participant("MeOH", labels["MeOH"], "solvent", evidence),
             participant("DCM_reductive_amination", labels["DCM"], "solvent", evidence),
             participant("NaBH4", labels["NaBH4"], "reagent", evidence, eq=q(1.1, "equiv", "1.1 equiv"), ref=smid),
             participant("product", pmid, "product", evidence, mass=q(mass_g, "g", f"{mass_g:g} g" if mass_g >= 1 else f"{mass_g*1000:g} mg"), amount=q(mmol, "mmol", f"{mmol:g} mmol"))]
    steps = [step("oxidation", evidence, "Add DMP to alcohol in DCM (0.25 M) and stir", temperature="rt", time=q(30, "min", "30 min"), atmosphere="argon"),
             step("oxidation_workup", evidence, "Quench with 1 M NaOH; extract DCM (2 × 30 mL); wash brine; dry Na2SO4; concentrate", stage="workup"),
             step("imine_formation", evidence, "Add p-anisidine in MeOH (0.25 M) or 1:1 MeOH/DCM; stir", temperature="rt", time=q(18, "h", "18 h"), atmosphere="argon"),
             step("reductant_addition", evidence, "Cool to 0 °C and add NaBH4 portionwise", temperature=q(0, "°C", "0 °C"), atmosphere="argon"),
             step("reduction_hold", evidence, "Warm to room temperature and stir", temperature="rt", time=q(2, "h", "2 h"), atmosphere="argon"),
             step("final_workup", evidence, "Quench water; extract DCM; wash brine; dry Na2SO4; concentrate", stage="workup"),
             step("purification", evidence, "Silica gel chromatography, hexanes/EtOAc 9:1 to 8:2", stage="purification")]
    make_reaction(f"prep_3{suffix}", page, f"Oxidation/reductive amination to 3{suffix}",
                  "benzyl alcohol oxidation and reductive amination", parts, steps, pmid, f"{y}%", y, [OXRA],
                  [("C", "N", 0, 1, "formed", None)], extra_evidence=ev(5),
                  notes=("The general procedure lists MeOH (0.25 M) or 1:1 MeOH/DCM; the substrate-specific solvent variant is unreported.",),
                  unresolved=("The reported yield covers the oxidation/reductive-amination sequence; the oxidation intermediate was used without purification.",))


for suffix, args in ox_data.items():
    oxra_record(suffix, *args)


# Page-8/9 S13 -> S14 N-alkylations and page-9 S14 -> S15 homologations.
for suffix, page, y, mass, mmol in [
    ("g", 8, 60, 1.0, 3.0), ("k", 9, 60, 1.0, 2.8),
    ("l", 9, 53, 0.600, 1.8), ("m", 9, 53, 0.600, 1.82),
]:
    smid = molecule("S13" + suffix, page, smiles=SMILES["S13" + suffix],
                    note="Corresponding S13 tosylamide is literature-derived; candidate graph inferred from S14 and the stated N-alkylation.")
    pmid = molecule("S14" + suffix, page)
    evidence = list(dict.fromkeys(ev(page) + ev(8)))
    parts = [participant("substrate", smid, "reactant", evidence, eq=q(1, "equiv", "1 equiv"), ref=smid),
             participant("K2CO3", labels["K2CO3"], "reagent", evidence, eq=q(2, "equiv", "2 equiv"), ref=smid),
             participant("propargyl_bromide", labels["propargyl bromide derivative"], "reagent", evidence, eq=q(1.5, "equiv", "1.5 equiv"), ref=smid),
             participant("DMF", labels["DMF"], "solvent", evidence),
             participant("product", pmid, "product", evidence, mass=q(mass, "g", f"{mass:g} g" if mass >= 1 else f"{mass*1000:g} mg"), amount=q(mmol, "mmol", f"{mmol:g} mmol"))]
    steps = copy.deepcopy(reactions_by_id[NALK2]["steps"])
    make_reaction(f"prep_S14{suffix}", page, f"N-alkylation of S13{suffix} to S14{suffix}", "N-propargylation",
                  parts, steps, pmid, f"{y}%", y, [NALK2], [("C", "N", 0, 1, "formed", None)],
                  extra_evidence=ev(8), unresolved=("S13 starting-material preparation is literature-derived.",))

for suffix, y, mass_g, mmol in [("g", 54, 0.450, 1.24), ("k", 72, 0.750, 2.0),
                                ("l", 50, 0.310, 0.9), ("m", 70, 0.435, 1.27)]:
    smid = molecule("S14" + suffix, 9)
    pmid = molecule("S15" + suffix, 9)
    evidence = ev(9)
    parts = [participant("substrate", smid, "reactant", evidence, eq=q(1, "equiv", "1 equiv"), ref=smid),
             participant("CuI", labels["CuI"], "catalyst", evidence, eq=q(0.5, "equiv", "0.5 equiv"), ref=smid),
             participant("paraformaldehyde", labels["(CHO)n"], "reagent", evidence, eq=q(2.5, "equiv", "2.5 equiv"), ref=smid),
             participant("Cy2NH", labels["dicyclohexylamine"], "reagent", evidence, eq=q(1.8, "equiv", "1.8 equiv"), ref=smid),
             participant("dioxane", labels["1,4-dioxane"], "solvent", evidence),
             participant("product", pmid, "product", evidence, mass=q(mass_g, "g", f"{mass_g*1000:g} mg"), amount=q(mmol, "mmol", f"{mmol:g} mmol"))]
    make_reaction(f"prep_S15{suffix}", 9, f"Homologation of S14{suffix} to S15{suffix}",
                  "terminal alkyne homologation to allene", parts, copy.deepcopy(reactions_by_id[HOM]["steps"]), pmid,
                  f"{y}%", y, [HOM], [("C", "C", 0, 2, "formed", None)], extra_evidence=ev(7),
                  unresolved=("New allene carbon is inferred to derive from paraformaldehyde; no atom mapping was reported.",))

for suffix, y, mass_g, mmol in [("g", 57, 0.330, 0.71), ("k", 56, 0.530, 1.13),
                                ("l", 52, 0.210, 0.47), ("m", 58, 0.330, 0.74)]:
    oxra_record(suffix, "S15" + suffix, 9, y, mass_g, mmol)


# Patch existing scope experiments with the actual general procedure, isolated-yield
# basis, product headings and SFC ee results. These are the same experiments because
# substrate/product/yield pairs match the earlier Scheme-3/Table-2 arrows.
sfc = {
    "2a": ("scheme3_3a_2a", 9, 56, "unknown", "30% MeOH; Phenomenex Amylose-1; 40 °C; CO2/MeOH 70:30; 1 mL/min; 210 nm; tR major 5.98 min, minor 6.96 min"),
    "2e": ("table2_entry_5", 10, 90, "R", "30% MeOH; Phenomenex Amylose-1; 40 °C; CO2/MeOH 70:30; 1 mL/min; 210 nm; tR major 19.26 min, minor 22.02 min"),
    "2f": ("scheme3_3f_2f", 10, 86, "R", "30% MeOH; Phenomenex Amylose-1; 40 °C; CO2/MeOH 70:30; 1 mL/min; 210 nm; tR major 33.15 min, minor 30.60 min"),
    "2g": ("scheme3_3g_2g", 10, 92, "R", "30% MeOH; Phenomenex Amylose-1; 40 °C; CO2/MeOH 70:30; 1 mL/min; 210 nm; tR major 19.78 min, minor 16.24 min"),
    "2h": ("scheme3_3h_2h", 10, 90, "R", "20% MeOH; Phenomenex Amylose-1; 40 °C; CO2/MeOH 80:20; 1 mL/min; 210 nm; tR major 15.05 min, minor 12.15 min"),
    "2i": ("scheme3_3i_2i", 10, 90, "R", "20% MeOH; Phenomenex Amylose-1; 40 °C; CO2/MeOH 80:20; 1 mL/min; 210 nm; tR major 15.71 min, minor 16.97 min"),
    "2j": ("scheme3_3j_2j", 10, 88, "R", "30% MeOH; Phenomenex Amylose-1; 40 °C; CO2/MeOH 70:30; 1 mL/min; 210 nm; tR major 24.37 min, minor 21.15 min"),
    "2k": ("scheme3_3k_2k", 10, 94, "R", "30% MeOH; Phenomenex Amylose-1; 40 °C; CO2/MeOH 70:30; 1 mL/min; 210 nm; tR major 20.09 min, minor 17.33 min"),
    "2l": ("scheme3_3l_2l", 10, 90, "R", "20% MeOH; Phenomenex Amylose-1; 40 °C; CO2/MeOH 80:20; 1 mL/min; 210 nm; tR major 19.08 min, minor 17.66 min"),
    "2m": ("scheme3_3m_2m", 10, 94, "R", "20% MeOH; Phenomenex Amylose-1; 40 °C; CO2/MeOH 80:20; 1 mL/min; 210 nm; tR major 20.98 min, minor 19.36 min"),
    "2n": ("scheme3_3n_2n", 11, 96, "R", "20% MeOH; Phenomenex Amylose-1; 40 °C; CO2/MeOH 80:20; 2 mL/min; 210 nm; tR major 11.50 min, minor 12.90 min"),
}
patched_scope_ids = []
for product_label, (suffix, page, ee_value, major_config, basis) in sfc.items():
    rid = f"{DOC}:{suffix}"
    r = reactions_by_id[rid]
    evidence = ev(page) + (ev(9) if product_label == "2e" else [])
    r["evidence_ids"] = list(dict.fromkeys(r["evidence_ids"] + evidence))
    if ASYM_ALLENE not in r["general_scheme_ids"]:
        r["general_scheme_ids"].append(ASYM_ALLENE)
    target = r["outcomes"][0]["target_molecule_ids"][0]
    product_m = molecules_by_id[target]
    product_m["occurrences"].append({"document_id": DOC,
                                     "compound_label": f(product_label, product_label, evidence),
                                     "evidence_ids": evidence})
    product_m["stereochemistry"]["chirality_status"] = inferred(
        "mixture", evidence, "The isolated sample has a reported nonzero ee and therefore contains both enantiomers.")
    product_m["stereochemistry"]["geometry_status"] = inferred(
        "not_applicable", evidence, "The product contains a terminal vinyl group; no E/Z or ring cis/trans descriptor is applicable in the printed name.")
    product_m["stereochemistry"]["representation_status"] = inferred(
        "not_encoded", evidence, "Absolute configuration is preserved textually; no source-verified stereochemical SMILES is asserted.")
    assignment = "major enantiomer configuration unassigned" if major_config == "unknown" else f"major enantiomer assigned ({major_config}) in the product heading"
    product_m["stereochemistry"]["assignment"] = f(assignment, assignment, evidence)
    product_m["stereochemistry"]["mixture"] = f(f"enantioenriched sample; {ee_value}% ee", f"enantioenriched sample; {ee_value}% ee", evidence)
    for o in r["outcomes"]:
        if o["kind"] == "yield":
            o["yield_type"] = f("isolated yield", "isolated yield", evidence)
            o["measurement_basis"] = f("isolated after silica gel chromatography", "isolated after silica gel chromatography", evidence)
        if o["kind"] == "selectivity" and sr.value(o["scope"]) and sr.value(o["scope"]).get("metric") == "er":
            o["measurement_basis"] = f("SFC; " + basis, "SFC; " + basis, evidence)
            scope = copy.deepcopy(sr.value(o["scope"]))
            scope["major_configuration"] = major_config
            o["scope"] = f(scope, scope, list(dict.fromkeys(o["scope"]["evidence_ids"] + evidence)))
    ee = sr.new_outcome("ee_characterization", "selectivity", [target], evidence)
    ee["product_role"] = f("product", "product", evidence)
    ee["value"] = f(f"{ee_value}% ee", q(ee_value, "%", f"{ee_value}% ee"), evidence)
    ee["measurement_basis"] = f("SFC; " + basis, "SFC; " + basis, evidence)
    ee_scope = {"metric": "ee", "major_configuration": major_config,
                "order": "major enantiomer:opposite enantiomer"}
    ee["scope"] = f(ee_scope, ee_scope, evidence)
    if not any(o["outcome_id"] == "ee_characterization" for o in r["outcomes"]):
        r["outcomes"].append(ee)
    r["review"]["unresolved"] = [u for u in r["review"]["unresolved"]
                                      if "Yield measurement basis" not in u and "Absolute configuration" not in u]
    if major_config == "unknown":
        r["review"]["unresolved"].append("Absolute configuration of the 2a major enantiomer remains unassigned.")
    if "Source-verified stereochemical SMILES remains unencoded." not in r["review"]["unresolved"]:
        r["review"]["unresolved"].append("Source-verified stereochemical SMILES remains unencoded.")
    patched_scope_ids.append(rid)


# Add explicit page-11 image evidence and detailed charges/steps to the four
# derivatization records already created from Scheme 4.
scheme4_ids = [f"{DOC}:scheme4_2e_to_4", f"{DOC}:scheme4_2e_to_5",
               f"{DOC}:scheme4_4_to_6", f"{DOC}:scheme4_6_to_7"]
for rid in scheme4_ids:
    r = reactions_by_id[rid]
    r["evidence_ids"] = list(dict.fromkeys(r["evidence_ids"] + ev(11)))
    add_evidence_recursive(r, [f"{DOC}:page11:image"], f"{DOC}:page11:text")


def set_charge(rid, pid, key, raw, normalized, evidence=ev(11)):
    r = reactions_by_id[rid]
    p = next(p for p in r["participants"] if p["participant_id"] == pid)
    p["charges"][key] = f(raw, normalized, evidence)


rid4 = f"{DOC}:scheme4_2e_to_4"
set_charge(rid4, "substrate", "mass", "110 mg", q(110, "mg", "110 mg"))
set_charge(rid4, "substrate", "amount", "0.25 mmol", q(0.25, "mmol", "0.25 mmol"))
set_charge(rid4, "CAN", "mass", "275 mg", q(275, "mg", "275 mg"))
set_charge(rid4, "CAN", "amount", "0.5 mmol", q(0.5, "mmol", "0.5 mmol"))
set_charge(rid4, "product", "mass", "77 mg", q(77, "mg", "77 mg"))
r = reactions_by_id[rid4]
me_cn = agent("MeCN", 11, smiles="CC#N")
water = labels["water"]
r["participants"].insert(-1, participant("MeCN", me_cn, "solvent", ev(11), volume=q(8, "mL", "8.0 mL")))
r["participants"].insert(-1, participant("water", water, "solvent", ev(11), volume=q(8, "mL", "8 mL")))
r["notes"].append(f("CAN charge is printed as 275 mg, 0.5 mmol, 2.5 equiv for 0.25 mmol 2e; the mmol ratio is 2.0 equiv, so the stated equivalents conflict is retained.",
                    "CAN charge is printed as 275 mg, 0.5 mmol, 2.5 equiv for 0.25 mmol 2e; the mmol ratio is 2.0 equiv, so the stated equivalents conflict is retained.", ev(11)))

rid5 = f"{DOC}:scheme4_2e_to_5"
for pid, vals in {"substrate": (("mass", "62 mg", q(62, "mg", "62 mg")), ("amount", "0.143 mmol", q(0.143, "mmol", "0.143 mmol"))),
                  "sodium": (("mass", "21.2 mg", q(21.2, "mg", "21.2 mg")), ("amount", "0.885 mmol", q(0.885, "mmol", "0.885 mmol"))),
                  "naphthalene": (("mass", "4 mg", q(4, "mg", "4 mg")), ("amount", "0.02 mmol", q(0.02, "mmol", "0.02 mmol"))),
                  "THF": (("volume", "5 mL total", q(5, "mL", "1 mL + 4 mL")),),
                  "product": (("mass", "34 mg", q(34, "mg", "34 mg")),)}.items():
    for key, raw, normalized in vals:
        set_charge(rid5, pid, key, raw, normalized)
r = reactions_by_id[rid5]
r["steps"] = [step("sodium_naphthalenide", ev(11), "Sonicate sodium sheets in naphthalene/THF until a green color persists", temperature="rt", atmosphere="argon"),
              step("substrate_addition", ev(11), "Add 2e solution in THF; rapid loss of green color observed", temperature="rt", atmosphere="argon"),
              step("reaction_hold", ev(11), "Remove from sonicator and stir", temperature="rt", time=q(15, "h", "15 h"), atmosphere="argon"),
              step("cold_quench", ev(11), "Cool to 0 °C; slowly add MeOH (10 mL), then saturated NaHCO3 after sodium consumption", temperature=q(0, "°C", "0 °C"), stage="workup"),
              step("purification", ev(11), "Dilute Et2O; wash NaHCO3 and H2O; dry Na2SO4; silica gel hexanes/EtOAc 8:2", stage="purification")]
for o in r["outcomes"]:
    if o["kind"] == "yield":
        o["yield_type"] = f("isolated yield", "isolated yield", ev(11))
        o["measurement_basis"] = f("isolated after silica gel chromatography", "isolated after silica gel chromatography", ev(11))
r["review"]["unresolved"] = [u for u in r["review"]["unresolved"] if "Yield measurement basis" not in u]

rid6 = f"{DOC}:scheme4_4_to_6"
for pid, vals in {"substrate": (("mass", "45 mg", q(45, "mg", "45 mg")), ("amount", "0.123 mmol", q(0.123, "mmol", "0.123 mmol"))),
                  "acryloyl_chloride": (("volume", "20 μL", q(20, "μL", "20 μL")), ("amount", "0.247 mmol", q(0.247, "mmol", "0.247 mmol"))),
                  "DMAP": (("mass", "1.5 mg", q(1.5, "mg", "1.5 mg")), ("amount", "0.012 mmol", q(0.012, "mmol", "0.012 mmol"))),
                  "Et3N": (("volume", "60 μL", q(60, "μL", "60 μL")), ("amount", "0.5 mmol", q(0.5, "mmol", "0.5 mmol"))),
                  "product": (("mass", "30 mg", q(30, "mg", "30 mg")),)}.items():
    for key, raw, normalized in vals:
        set_charge(rid6, pid, key, raw, normalized)

rid7 = f"{DOC}:scheme4_6_to_7"
for pid, vals in {"prose_substrate_5": (("mass", "15 mg", q(15, "mg", "15 mg")), ("amount", "0.04 mmol", q(0.04, "mmol", "0.04 mmol"))),
                  "catalyst": (("mass", "2.5 mg", q(2.5, "mg", "2.5 mg")), ("amount", "0.004 mmol", q(0.004, "mmol", "0.004 mmol"))),
                  "DCM": (("volume", "1.5 mL", q(1.5, "mL", "1.5 mL")),),
                  "product": (("mass", "12 mg", q(12, "mg", "12 mg")),)}.items():
    for key, raw, normalized in vals:
        set_charge(rid7, pid, key, raw, normalized)
r = reactions_by_id[rid7]
target7 = r["outcomes"][0]["target_molecule_ids"][0]
ee7 = sr.new_outcome("ee_characterization", "selectivity", [target7], ev(11))
ee7["product_role"] = f("product", "product", ev(11))
ee7["value"] = f("88% ee", q(88, "%", "88% ee"), ev(11))
ee7["measurement_basis"] = f("SFC; 40% MeOH; Phenomenex Amylose-1; 40 °C; CO2/MeOH 60:40; 1 mL/min; 210 nm; tR major 9.12 min, minor 11.06 min",
                               "SFC; 40% MeOH; Phenomenex Amylose-1; 40 °C; CO2/MeOH 60:40; 1 mL/min; 210 nm; tR major 9.12 min, minor 11.06 min", ev(11))
ee7["scope"] = f({"metric": "ee", "major_configuration": "R", "order": "major enantiomer:opposite enantiomer"},
                 {"metric": "ee", "major_configuration": "R", "order": "major enantiomer:opposite enantiomer"}, ev(11))
if not any(o["outcome_id"] == "ee_characterization" for o in r["outcomes"]):
    r["outcomes"].append(ee7)
mol7 = molecules_by_id[target7]
mol7["stereochemistry"]["chirality_status"] = inferred("mixture", ev(11), "The isolated product has 88% ee.")
mol7["stereochemistry"]["geometry_status"] = inferred("not_applicable", ev(11), "No E/Z or ring cis/trans descriptor is applicable in the printed product name.")
mol7["stereochemistry"]["representation_status"] = inferred("not_encoded", ev(11), "The (R) assignment is textual and no source-verified stereochemical SMILES is asserted.")
mol7["stereochemistry"]["assignment"] = f("major enantiomer assigned (R) in product heading", "major enantiomer assigned (R) in product heading", ev(11))
mol7["stereochemistry"]["mixture"] = f("enantioenriched sample; 88% ee", "enantioenriched sample; 88% ee", ev(11))


# Update labels 3f--3n with page-8/9 names and candidates; standard stereo
# diagnostics are re-run for every concrete graph below.
for lab in ["3f", "3h", "3i", "3j", "3n"]:
    molecule(lab, 8, smiles=SMILES[lab])
for lab in ["3g", "3k", "3l", "3m"]:
    molecule(lab, 9, smiles=SMILES[lab])


# Source context: spectral peaks and bibliographic matter are not synthetic records.
bundle["source_context"].extend([
    {"context_id": f"{DOC}:context_characterization_pages7_11", "document_id": DOC,
     "pages": [7, 8, 9, 10, 11], "kind": "characterization",
     "evidence_ids": ev(7) + ev(8) + ev(9) + ev(10) + ev(11),
     "note": "NMR, MS/HRMS, optical rotations and SFC peak lists support molecule/sample facts; mass-spectrum peaks and fragment abundances are not synthetic experiments."},
    {"context_id": f"{DOC}:context_literature_S13", "document_id": DOC, "pages": [8, 9],
     "kind": "literature_preparations", "evidence_ids": ev(8) + ev(9),
     "note": "Tosylamides S13 are explicitly stated to be synthesized according to literature; no current-paper S13 preparation was created."},
    {"context_id": f"{DOC}:context_SI_availability", "document_id": DOC, "pages": [11],
     "kind": "associated_content_notice", "evidence_ids": ev(11),
     "note": "The main paper prints an SI availability link, but no SI file was supplied in this workspace; no SI document was associated or used."},
    {"context_id": f"{DOC}:context_authors_acknowledgments_references", "document_id": DOC, "pages": [11, 12],
     "kind": "bibliographic_context", "evidence_ids": ev(11) + ev(12),
     "note": "Author information, acknowledgments and cited references; no current-paper experiments on page 12."},
])


# Ensure all records remain machine review only and all newly registered concrete
# molecules carry a standard stereo diagnostic artifact.
for r in bundle["reactions"]:
    r["review"]["status"] = "needs_review"
    r["review"]["method"] = "model"
    r["review"]["reviewer"] = None

diagnostics = {}
for m in bundle["molecules"]:
    smi = sr.value(m["smiles"])
    if smi:
        diagnostics[m["molecule_id"]] = sr.stereo_inspect(smi)

# Bounded source inventory retained separately from the formal bundle.
new_ids = [r["reaction_id"] for r in bundle["reactions"] if any(f":page{p}:" in eid for p in range(7, 13) for eid in r["evidence_ids"])]
raw_inventory = {
    "document_id": DOC,
    "pages": [7, 8, 9, 10, 11, 12],
    "visual_inspection": "Each supplied full-page JPEG was viewed at original detail.",
    "new_specific_record_ids": [rid for rid in new_ids if rid not in patched_scope_ids + scheme4_ids and reactions_by_id[rid]["record_type"] == "specific_experiment"],
    "patched_existing_record_ids": patched_scope_ids + scheme4_ids,
    "generic_record_ids": [HOM, NALK2, ASYM_ALLENE],
    "non_experiment_regions": ["spectral peak lists pages 7-11", "SI availability notice page 11", "authors/acknowledgments/references page 12"],
    "source_conflicts": [
        "CAN is printed as 0.5 mmol and 2.5 equiv for 0.25 mmol 2e (arithmetic gives 2.0 equiv).",
        "Scheme 4 states Et3N 2 equiv for 4->6, while the detailed procedure states 4 equiv.",
        "Scheme 4 depicts 6->7, while the detailed RCM prose says a solution of 5 was used.",
    ],
}


# Exact six-page coverage. Existing scope/derivatization record IDs are included
# because these pages add evidence to the same experiments.
by_page = {p: [] for p in range(7, 13)}
for r in bundle["reactions"]:
    for p in by_page:
        if any(f":page{p}:" in eid for eid in r["evidence_ids"]):
            by_page[p].append(r["reaction_id"])
for p in by_page:
    by_page[p] = list(dict.fromkeys(by_page[p]))

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 7, "state": "attempted", "regions": [
            {"region": "upper-left compound paragraphs", "finding": "four specific N-alkylation products S10h/S10i/S10j/S10n with yields and isolated amounts; spectral peak lists retained as characterization", "evidence_ids": ev(7)},
            {"region": "lower-left Homologation of Alkynes to Allenes", "finding": "generic quantified homologation procedure and S11d/S11e/S11f preparations", "evidence_ids": ev(7)},
            {"region": "right column S11/S12 paragraphs", "finding": "S11h/S11i/S11j/S11n homologations and ester reductions to 3a/S12d/S12f", "evidence_ids": ev(7)},
        ], "record_ids": by_page[7], "unresolved": ["S12d isolated mass and mmol are not reported.", "Literature preparation of S1h/S1i/S1j/S1n was not inferred as current work."]},
        {"page": 8, "state": "attempted", "regions": [
            {"region": "left column S12h/S12i/S12j/S12n", "finding": "four ester-reduction products with yields and isolated amounts; spectral data are characterization", "evidence_ids": ev(8)},
            {"region": "middle/right compound paragraphs 3d/3e/3f/3h/3i/3j/3n", "finding": "seven oxidation/reductive-amination sequence products", "evidence_ids": ev(8)},
            {"region": "lower-right Preparation of Allenes 3 / N-Alkylation", "finding": "literature S13 statement plus a distinct K2CO3/DMF N-alkylation procedure; S14g paragraph starts", "evidence_ids": ev(8)},
        ], "record_ids": by_page[8], "unresolved": ["Specific reactions do not identify whether MeOH alone or 1:1 MeOH/DCM was used in reductive amination.", "S13 preparations are literature-derived and not counted as experiments in this paper."]},
        {"page": 9, "state": "attempted", "regions": [
            {"region": "left column S14k/S14l/S14m and S15g/S15k/S15l/S15m", "finding": "three continued N-alkylations and four alkyne-to-allene homologations; S14g continuation is merged with page 8", "evidence_ids": ev(9)},
            {"region": "upper-right 3g/3k/3l/3m", "finding": "four oxidation/reductive-amination sequence products", "evidence_ids": ev(9)},
            {"region": "lower-right asymmetric cyclization procedure", "finding": "generic allene cyclization procedure with alternative Brønsted acids, explicit loadings, temperature, time and purification", "evidence_ids": ev(9)},
            {"region": "bottom product 2a and start of 2e", "finding": "duplicate characterization evidence for existing scope experiments; isolated-yield and SFC ee facts merged, not duplicated", "evidence_ids": ev(9)},
        ], "record_ids": by_page[9], "unresolved": ["Full stereochemical SMILES of cyclization products remains unencoded.", "PPTS versus chloroacetic acid is selected per substrate and is not a combined reagent charge."]},
        {"page": 10, "state": "attempted", "regions": [
            {"region": "both columns product paragraphs 2e-2n", "finding": "isolated-product characterization for existing asymmetric cyclization scope experiments, including (R) headings, ee and SFC methods; merged into existing records", "evidence_ids": ev(10)},
        ], "record_ids": by_page[10], "unresolved": ["Product connectivity is retained from earlier scope labels/names; source-verified stereochemical SMILES is not asserted.", "2n characterization continues on page 11."]},
        {"page": 11, "state": "attempted", "regions": [
            {"region": "top-left continuation of 2n", "finding": "completion of existing 2n characterization and SFC method", "evidence_ids": ev(11)},
            {"region": "Derivatization paragraphs 2e->4, 2e->5, 4->6, RCM->7", "finding": "detailed operations and reagent quantities merged into the four existing Scheme-4 records", "evidence_ids": ev(11)},
            {"region": "Associated Content and accession code", "finding": "SI availability and crystallographic context only; no SI file was supplied or associated", "evidence_ids": ev(11)},
            {"region": "Author Information", "finding": "bibliographic context, not an experiment", "evidence_ids": ev(11)},
        ], "record_ids": by_page[11], "unresolved": ["CAN 0.5 mmol versus stated 2.5 equiv conflict retained.", "Et3N 2 equiv in Scheme 4 versus 4 equiv in detailed prose retained.", "RCM substrate is 6 in Scheme 4 but printed as 5 in detailed prose; no correction was imposed."]},
        {"page": 12, "state": "attempted", "regions": [
            {"region": "author continuation, notes and acknowledgments", "finding": "bibliographic/administrative content; no current synthetic experiment", "evidence_ids": ev(12)},
            {"region": "references 1-18", "finding": "cited literature only; no current-paper reaction records created", "evidence_ids": ev(12)},
        ], "record_ids": by_page[12], "unresolved": []},
    ],
}

(RUN / "coverage").mkdir(parents=True, exist_ok=True)
sr.write(RUN / "raw_model/pages_0007_0012_extraction.json", raw_inventory)
sr.write(RUN / "raw_model/pages_0007_0012_stereo_diagnostics.json", diagnostics)
sr.write(RUN / "coverage/pages_0007_0012.json", coverage)
sr.write(RUN / "submitted.json", bundle)

print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
                  "new_specific": len(raw_inventory["new_specific_record_ids"]),
                  "patched_existing": len(raw_inventory["patched_existing_record_ids"]),
                  "generic_new": len(raw_inventory["generic_record_ids"])}, ensure_ascii=False))
