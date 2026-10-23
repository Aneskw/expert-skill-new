#!/usr/bin/env python3
"""Build the page 1--6 evidence-linked v3.1 submission for this paper.

This is a source-derived machine extraction.  It intentionally leaves uncertain
structures and unreported measurement bases unresolved; merge will quarantine
all records as needs_review.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as sr  # noqa: E402

RUN = ROOT / "runs/current"
DATASET = json.loads((RUN / "dataset.json").read_text())
DOC = "doc_036d1c8f4ae80aff"


def ev(page: int, include_text: bool = True) -> list[str]:
    result = [f"{DOC}:page{page}:image"]
    if include_text:
        result.append(f"{DOC}:page{page}:text")
    return result


def f(raw, normalized=None, evidence=(), method="original_extraction", status="pending_review", derivation=None):
    if normalized is None:
        normalized = raw
    return sr.fact(raw=raw, normalized=normalized, evidence=evidence, method=method,
                   status=status, derivation=derivation)


def missing(evidence=(), status="unreported"):
    return sr.fact(evidence=evidence, status=status)


def inferred(normalized, evidence, derivation, raw=None):
    return sr.fact(raw=raw, normalized=normalized, evidence=evidence, method="inference",
                   status="pending_review", verification="inferred", derivation=derivation)


def conflict(alternatives, evidence):
    return sr.fact(evidence=evidence, status="conflict", alternatives=alternatives)


bundle = sr.empty(DATASET["documents"])
molecules: dict[str, dict] = {}
reactions: list[dict] = []


SMILES = {
    "1a": "CC#CN(c1ccccc1CO)S(=O)(=O)c1ccc(C)cc1",
    "1b": "CC#CCOc1ccccc1CNc1ccc(OC)cc1",
    "1c": "CC#CCCc1ccccc1CNc1ccc(OC)cc1",
    "1d": "CC#CN(c1ccccc1CNc1ccc(OC)cc1)C(=O)OC(C)(C)C",
    "1e": "CC#CN(c1ccccc1CNc1ccc(OC)cc1)S(=O)(=O)c1ccc(C)cc1",
    "3a": "C=C=CCN(c1ccccc1CO)S(=O)(=O)c1ccc(C)cc1",
    "3d": "C=C=CCN(c1ccccc1CNc1ccc(OC)cc1)C(=O)OC(C)(C)C",
    "3e": "C=C=CCN(c1ccccc1CNc1ccc(OC)cc1)S(=O)(=O)c1ccc(C)cc1",
    "S1a": "CCOC(=O)c1ccccc1NC(=O)OC(C)(C)C",
    "S2d": "CC#CN(c1ccccc1C(=O)OCC)C(=O)OC(C)(C)C",
    "S2e": "CC#CN(c1ccccc1C(=O)OCC)S(=O)(=O)c1ccc(C)cc1",
    "S3d": "CC#CN(c1ccccc1CO)C(=O)OC(C)(C)C",
    "S4": "CC#CCOc1ccccc1C=O",
    "S6": "CCOC(OCC)c1ccccc1CCC=O",
    "S7": "C#CCCc1ccccc1C(OCC)OCC",
    "S8": "CC#CCCc1ccccc1C(OCC)OCC",
    "S9": "CC#CCCc1ccccc1C=O",
    "S10d": "C#CN(c1ccccc1C(=O)OCC)C(=O)OC(C)(C)C",
    "S10e": "C#CN(c1ccccc1C(=O)OCC)S(=O)(=O)c1ccc(C)cc1",
    "S10f": "C#CN(c1c(C(=O)OC)cc(C)cc1)S(=O)(=O)c1ccc(C)cc1",
    "ethyl 2-aminobenzoate": "CCOC(=O)c1ccccc1N",
    "Boc2O": "CC(C)(C)OC(=O)OC(=O)OC(C)(C)C",
    "DMAP": "CN(C)c1ccncc1",
    "Et3N": "CCN(CC)CC",
    "NaH": "[Na+].[H-]",
    "propargyl bromide derivative": "C#CCBr",
    "but-2-yn-1-yl bromide": "CC#CCBr",
    "DIBAL-H": "CC(C)[AlH]C(C)C",
    "DMP": "CC(=O)OI1(OC(C)=O)OC(=O)c2ccccc21",
    "p-anisidine": "COc1ccc(N)cc1",
    "NaBH4": "[Na+].[BH4-]",
    "K2CO3": "[K+].[K+].[O-]C([O-])=O",
    "salicylaldehyde": "O=Cc1ccccc1O",
    "1-bromo-2-butyne": "CC#CCBr",
    "Pd(OAc)2": "CC(=O)[O-].CC(=O)[O-].[Pd+2]",
    "NaHCO3": "[Na+].O=C([O-])O",
    "Bu4NCl": "CCCC[N+](CCCC)(CCCC)CCCC.[Cl-]",
    "DIPA": "CC(C)NC(C)C",
    "trimethylsilyldiazomethane": "C[Si](C)(C)C=[N+]=[N-]",
    "nBuLi": "[Li]CCCC",
    "MeI": "CI",
    "PPTS": "Cc1ccc(S(=O)(=O)O)cc1.c1cc[nH+]cc1",
    "TFA": "O=C(O)C(F)(F)F",
    "chloroacetic acid": "O=C(O)CCl",
    "DCE": "ClCCCl",
    "THF": "C1CCOC1",
    "DMF": "CN(C)C=O",
    "DCM": "ClCCl",
    "MeOH": "CO",
    "acetone": "CC(C)=O",
    "water": "O",
    "Et2O": "CCOCC",
    "CAN": "[Ce+4]",
    "HCl": "Cl",
    "sodium metal": "[Na]",
    "naphthalene": "c1ccc2ccccc2c1",
    "acryloyl chloride": "C=CC(=O)Cl",
}


def mol(key: str, label: str | None = None, page: int = 2, *, smiles: str | None = None,
        chirality="unresolved", geometry="unresolved", representation=None,
        assignment=None, mixture=None, note=None) -> str:
    if key in molecules:
        return molecules[key]["molecule_id"]
    label = label or key
    eid = ev(page)
    mid = f"{DOC}:mol_{key.replace(' ', '_').replace('/', '_')}"
    m = sr.new_molecule(mid, DOC, eid)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, eid)
    m["names"]["original"] = f(label, label, eid)
    m["identity_note"] = f(note, note, eid) if note else inferred(
        "Identity retained at the source-label level.", eid,
        "The source label/drawing was inspected, but no external identity lookup was performed.")
    smi = smiles if smiles is not None else SMILES.get(label)
    if smi:
        diag = sr.stereo_inspect(smi)
        m["representation_type"] = "smiles"
        m["smiles"] = sr.fact(raw=label, normalized=diag["canonical_isomeric_smiles"], evidence=eid,
                              method="structure_derivation", status="pending_review",
                              verification="inferred",
                              derivation="Candidate transcribed from the source name/drawing; RDKit parsing is not source verification.")
        has_unspecified_potential = any(p["specified"] == "Unspecified" for p in diag["potential_stereo"])
        rep = representation or ("partial" if has_unspecified_potential else "encoded")
        if chirality == "unresolved":
            has_atom = any(p["type"].startswith("Atom_") for p in diag["potential_stereo"])
            chirality = "unspecified" if has_atom else "achiral"
        if geometry == "unresolved":
            has_bond = any(p["type"].startswith("Bond_") for p in diag["potential_stereo"])
            geometry = "unspecified" if has_bond else "not_applicable"
        m["stereochemistry"]["isomeric_smiles"] = sr.fact(
            raw=smi, normalized=diag["canonical_isomeric_smiles"], evidence=eid,
            method="structure_derivation", status="pending_review", verification="inferred",
            derivation="Canonical isomeric SMILES generated by standard stereo diagnostics.")
        m["stereochemistry"]["identity_key"] = sr.fact(
            raw=diag["canonical_isomeric_smiles"], normalized=diag["stereo_graph_key"], evidence=eid,
            method="structure_derivation", status="pending_review", verification="inferred",
            derivation="SHA-256 stereo graph key from standard stereo diagnostics; not a sample identity.")
    else:
        m["representation_type"] = "unresolved"
        rep = representation or "unresolved"
    stereoev = eid
    for field_name, value in (("chirality_status", chirality), ("geometry_status", geometry),
                              ("representation_status", rep)):
        m["stereochemistry"][field_name] = inferred(
            value, stereoev,
            "Status assigned from the inspected source drawing/name plus standard stereo graph diagnostics where a candidate SMILES exists.")
    if assignment is not None:
        m["stereochemistry"]["assignment"] = f(assignment, assignment, stereoev)
    if mixture is not None:
        m["stereochemistry"]["mixture"] = f(mixture, mixture, stereoev)
    molecules[key] = m
    return mid


def q(number, unit, raw=None, qualifier=None):
    out = {"number": number, "unit": unit}
    if raw is not None:
        out["raw_text"] = raw
    if qualifier:
        out["qualifier"] = qualifier
    return out


def participant(pid, mid, role, evidence, *, eq=None, molpct=None, mass=None, amount=None,
                volume=None, ref=None, stage="reaction"):
    p = sr.new_participant(pid, mid, role, evidence, stage)
    for name, value in (("equivalents", eq), ("mol_percent", molpct), ("mass", mass),
                        ("amount", amount), ("volume", volume)):
        if value is not None:
            if isinstance(value, dict) and "normalized_value" in value and "acquisition_method" in value:
                p["charges"][name] = value
                continue
            raw = value.get("raw_text") if isinstance(value, dict) else str(value)
            p["charges"][name] = f(raw or value, value, evidence)
    if ref:
        p["equivalence_reference"] = f(ref, ref, evidence)
    return p


def step(step_id, evidence, operation, *, temperature=None, time=None, atmosphere=None,
         participant_ids=(), stage="reaction"):
    s = sr.new_step(step_id, evidence, stage)
    s["participant_ids"] = list(participant_ids)
    s["operation"] = f(operation, operation, evidence)
    if temperature is not None:
        raw = temperature.get("raw_text") if isinstance(temperature, dict) else str(temperature)
        s["temperature"] = f(raw or temperature, temperature, evidence)
    if time is not None:
        raw = time.get("raw_text") if isinstance(time, dict) else str(time)
        s["time"] = f(raw or time, time, evidence)
    if atmosphere is not None:
        s["atmosphere"] = f(atmosphere, atmosphere, evidence)
    return s


CHANGESETS = {
    "C-N+H": [("C", "N", 0, 1, "formed"), ("C", "H", 0, 1, "hydrogen_gain")],
    "C-O+H": [("C", "O", 0, 1, "formed"), ("C", "H", 0, 1, "hydrogen_gain")],
    "C-N": [("C", "N", 0, 1, "formed")],
    "C-O": [("C", "O", 0, 1, "formed")],
    "C-C": [("C", "C", 0, 1, "formed")],
    "C-C+TB-C-C": [("C", "C", 0, 1, "formed"), ("C", "C", 1, 3, "order_increased")],
    "C-H": [("C", "H", 0, 1, "hydrogen_gain")],
    "N-H": [("N", "H", 0, 1, "hydrogen_gain")],
    "DB-C-O": [("C", "O", 1, 2, "order_increased")],
    "C-C+DB-C-C": [("C", "C", 0, 1, "formed"), ("C", "C", 1, 2, "order_increased")],
}


def set_classification(r, kind, evidence):
    changes = CHANGESETS.get(kind)
    if not changes:
        return
    labels = {}
    for idx, (a, b, before, after, ctype) in enumerate(changes, 1):
        cid = f"b{idx}"
        change = {"change_id": cid, "atoms": [{"atom_map": None, "element": a},
                                                {"atom_map": None, "element": b}],
                  "before_order": before, "after_order": after, "change_type": ctype}
        if "hydrogen" in ctype:
            change["hydrogen_count"] = 1
        r["classification"]["bond_changes"].append(inferred(
            change, evidence, "Net bond change inferred from source reactant/product drawings; atom mapping remains unresolved.",
            raw=f"source drawings imply {ctype} {a}-{b}"))
        if ctype in {"formed", "order_increased", "hydrogen_gain"}:
            lab = sr.bond_label(a, b, after)
            labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(cid)
    r["classification"]["labels"] = inferred(
        list(labels.values()), evidence,
        "Labels calculated from the explicitly retained net bond changes under bond-formation-20261016.",
        raw="formation labels inferred from inspected reaction arrow")
    r["classification"]["coverage"] = inferred(
        "partial", evidence,
        "Only net changes clear from the source drawings are represented; no complete atom mapping or coproduct accounting was performed.")
    r["classification"]["limitations"] = [
        "No atom-mapped reaction SMILES was asserted.",
        "Hydrogen gains are net changes and do not identify the proton source or elementary mechanism.",
        "Unshown coproducts and all bond-order rearrangements outside the reaction center remain unresolved.",
    ]


def reaction(suffix, page, *, categories, name, transformation, location, participants=(),
             steps=(), outcomes=(), general=(), class_kind=None, notes=(), unresolved=(),
             evidence_extra=()):
    evidence = list(dict.fromkeys(ev(page) + list(evidence_extra)))
    rid = f"{DOC}:{suffix}"
    r = sr.new_reaction(rid, DOC, evidence)
    r["record_type"] = "specific_experiment"
    r["experiment_categories"] = f(categories, categories, evidence)
    r["location"] = f(location, location, evidence)
    r["reaction_name"] = f(name, name, evidence)
    r["transformation_type"] = f(transformation, transformation, evidence)
    r["general_scheme_ids"] = list(general)
    r["participants"] = list(participants)
    r["steps"] = list(steps)
    r["outcomes"] = list(outcomes)
    r["notes"] = [f(n, n, evidence) for n in notes]
    r["review"]["unresolved"] = list(unresolved) or [
        "Machine extraction requires expert comparison of every structure and stereochemical depiction."]
    if class_kind:
        set_classification(r, class_kind, evidence)
    reactions.append(r)
    return r


def outcome(oid, kind, target, raw, normalized, evidence, *, product_role="product",
            yield_type=None, basis=None, scope=None):
    o = sr.new_outcome(oid, kind, [target], evidence)
    o["product_role"] = f(product_role, product_role, evidence)
    o["value"] = f(raw, normalized, evidence)
    if yield_type is not None:
        o["yield_type"] = f(yield_type, yield_type, evidence)
    if basis is not None:
        o["measurement_basis"] = f(basis, basis, evidence)
    if scope is not None:
        o["scope"] = f(scope, scope, evidence)
    return o


def generic(suffix, page, name, location, *, record_type="generic_scheme", notes=()):
    evidence = ev(page)
    r = sr.new_reaction(f"{DOC}:{suffix}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = f([], [], evidence)
    r["location"] = f(location, location, evidence)
    r["reaction_name"] = f(name, name, evidence)
    r["transformation_type"] = f(name, name, evidence)
    r["notes"] = [f(n, n, evidence) for n in notes]
    r["review"]["unresolved"] = [
        "Generic or proposed content is retained separately and is not counted as a specific experiment."
    ]
    reactions.append(r)
    return r["reaction_id"]


# Generic/current-work arrows, procedures, and the proposed mechanism are not experiments.
G_THIS = generic("scheme1_this_work", 2, "Rh-catalyzed hydrofunctionalization of alkynes/allene",
                 "Page 2, Scheme 1, This work",
                 notes=("The three upper panels are cited prior work and remain source context.",))
G_NALK = generic("gp_n_alkylation", 4, "General procedure for N-alkylation",
                 "Page 4, Experimental Section")
G_REDUCE = generic("gp_ester_reduction", 4, "General procedure for ester reduction",
                   "Pages 4-5, Experimental Section")
G_OXRA = generic("gp_oxidation_reductive_amination", 5,
                 "Alcohol oxidation followed by reductive amination",
                 "Page 5, Experimental Section")
G_RAC = generic("gp_racemic_alkyne_cyclization", 6,
                "General procedure for racemic cyclization of alkynes 1",
                "Page 6, Experimental Section")
G_ASYM = generic("gp_asymmetric_alkyne_cyclization", 6,
                 "General procedure for asymmetric cyclization of alkynes 1",
                 "Page 6, Experimental Section")
generic("proposed_inner_outer_mechanism", 3,
        "Proposed inner-sphere versus outer-sphere hydrofunctionalization pathways",
        "Page 3, mechanistic discussion", record_type="proposed_mechanism",
        notes=("No dedicated mechanistic experiment is reported in the reviewed region.",))


# Shared molecules for the hydrofunctionalization experiments.
for label, page in [("1a", 2), ("1b", 2), ("1c", 2), ("1d", 2), ("1e", 2),
                    ("3a", 3), ("3d", 3), ("3e", 3)]:
    mol(label, label, page)


def agent(label, page=2):
    return mol("agent_" + label, label, page)


RH = agent("[Rh(cod)Cl]2")
RAC_BINAP = agent("rac-BINAP")
RAC_BNP = agent("rac-BNP acid")
DCE = agent("DCE")


def hydro_participants(substrate, product_mid, page, ligand, acid, ligand_pct, acid_pct,
                       *, substrate_amount=None, concentration=None):
    evidence = ev(page)
    smid = mol(substrate, substrate, page)
    lmid = agent(ligand, page)
    amid = agent(acid, page)
    ps = [participant("substrate", smid, "reactant", evidence,
                      amount=substrate_amount, ref=None),
          participant("rh_catalyst", RH, "catalyst", evidence,
                      molpct=q(4, "mol%", "4 mol%"), ref=smid),
          participant("ligand", lmid, "ligand", evidence, molpct=ligand_pct, ref=smid),
          participant("bronsted_acid", amid, "additive", evidence, molpct=acid_pct, ref=smid),
          participant("dce", DCE, "solvent", evidence),
          participant("product", product_mid, "product", evidence)]
    return ps


# Scheme 2: four current-work experiments. Page 6 repeats the procedures/product data.
scheme2 = [
    ("scheme2_1a_2a", "1a", "2a_rac", "2a", 85, "C-O+H", "intramolecular hydroalkoxylation"),
    ("scheme2_1b_2b", "1b", "2b_rac", "2b", 41, "C-N+H", "intramolecular hydroamination"),
    ("scheme2_1c_2c", "1c", "2c_rac", "2c (with 2c′ byproduct)", 57, "C-N+H", "intramolecular hydroamination"),
    ("scheme2_1d_2d", "1d", "2d_rac", "2d", 70, "C-N+H", "intramolecular hydroamination"),
]
for suffix, substrate, pkey, plabel, yld, ckind, ttype in scheme2:
    assignment = "racemic cyclization product; absolute configuration not assigned"
    pmid = mol(pkey, plabel, 2, chirality="racemic", geometry="unspecified",
               representation="not_encoded", assignment=assignment)
    parts = hydro_participants(substrate, pmid, 2, "rac-BINAP", "rac-BNP acid",
                               q(8, "mol%", "8 mol%"), q(8, "mol%", "8 mol%"),
                               substrate_amount=q(0.2, "mmol", "0.2 mmol"))
    st = [step("catalyst_preformation", ev(6), "DCE added to Rh/acid/ligand; stir",
               temperature="rt", time=q(10, "min", "10 min"), atmosphere="argon",
               participant_ids=[p["participant_id"] for p in parts[:5]]),
          step("cyclization", ev(2) + ev(6), "Add alkyne under argon and heat",
               temperature=q(70, "°C", "70 °C"), time=q(24, "h", "24 h"), atmosphere="argon",
               participant_ids=["substrate", "rh_catalyst", "ligand", "bronsted_acid", "dce"]),
          step("purification", ev(6), "Evaporate solvent; silica gel, EtOAc/hexanes 1:9",
               participant_ids=["product"], stage="purification")]
    outs = [outcome("yield", "yield", pmid, f"{yld}%", q(yld, "%", f"{yld}%"), ev(2) + ev(6),
                    yield_type="reported product yield")]
    notes = []
    unresolved = ["Yield measurement basis is not explicitly stated for the Scheme 2 value.",
                  "Product stereostructure candidate was not encoded from the drawing."]
    if substrate == "1c":
        bymid = mol("2c_prime_mix", "2c′", 6, chirality="mixture", geometry="mixture",
                    representation="not_encoded", mixture="mixture of E/Z isomers; composition not reported")
        parts.append(participant("byproduct", bymid, "product", ev(2) + ev(6)))
        outs.append(outcome("byproduct_yield", "yield", bymid, "20%", q(20, "%", "20%"), ev(2) + ev(6),
                            product_role="byproduct", yield_type="reported byproduct yield",
                            scope={"identity": "2c′", "mixture": "E/Z isomers"}))
        notes.append("The prose and page-6 characterization report 2c′ as a 20% E/Z-isomer mixture; it is not merged with 2c.")
    reaction(suffix, 2, categories=["other"], name=f"Cyclization of {substrate} to {plabel}",
             transformation=ttype, location="Page 2 Scheme 2; page 6 experimental repetition",
             participants=parts, steps=st, outcomes=outs, general=[G_RAC], class_kind=ckind,
             notes=notes, unresolved=unresolved, evidence_extra=ev(6))


# Table 1 optimization: each row is a distinct experiment.  The table says 10 mol%
# ligand/acid, whereas the page-6 general procedure charges 8 mol%; preserve conflict.
table1 = [
    (1, "1d", "Josiphos-SL-J002-1", "rac-BNP acid", "<5%", None, 50, 24),
    (2, "1d", "(R)-BINAP", "PPTS", "traces", None, 50, 24),
    (3, "1d", "(R)-BINAP", "rac-BNP acid", "19%", "57:43", 50, 24),
    (4, "1d", "(R)-BINAP", "rac-BNP acid", "58%", "57:43", 50, "5 days"),
    (5, "1d", "(R)-DTBM-Segphos", "rac-BNP acid", "73%", "57:43", 50, 24),
    (6, "1d", "(S)-DTBM-Garphos", "rac-BNP acid", "81%", "67:33", 50, 24),
    (7, "1e", "(R)-DTBM-Garphos", "TFA", "50%", "80:20", 50, 24),
    (8, "1e", "(R)-DTBM-Garphos", "TFA", "60%", "60:40", 70, 24),
]
for entry, substrate, ligand, acid, yraw, er, temp, duration in table1:
    pkey = f"t1e{entry}_product"
    assignment = "major enantiomer unassigned"
    if entry == 6:
        assignment = "(S)-2d is the major enantiomer (Table 1 footnote c)"
    pmid = mol(pkey, "2d" if substrate == "1d" else "2e", 2,
               chirality="mixture", geometry="unspecified", representation="not_encoded",
               assignment=assignment, mixture=f"enantiomeric mixture; er {er}" if er else None)
    conflict_pct = conflict([
        f("10 mol% (Table 1 footnote)", q(10, "mol%", "10 mol%"), ev(2)),
        f("16 μmol per 0.2 mmol = 8 mol% (general procedure)", q(8, "mol%", "8 mol%"), ev(6),
          method="inference", derivation="16 μmol / 0.2 mmol × 100"),
    ], ev(2) + ev(6))
    parts = hydro_participants(substrate, pmid, 2, ligand, acid, conflict_pct, conflict_pct,
                               substrate_amount=q(0.2, "mmol", "0.2 mmol"), concentration=q(0.4, "M", "0.4 M"))
    stime = q(5, "d", "5 days") if duration == "5 days" else q(duration, "h", f"{duration} h")
    st = [step("cyclization", ev(2) + ev(6), "Asymmetric cyclization under Table 1 conditions",
               temperature=q(temp, "°C", f"{temp} °C"), time=stime, atmosphere="argon",
               participant_ids=[p["participant_id"] for p in parts[:-1]])]
    if yraw == "traces":
        outs = [outcome("trace_product", "trace", pmid, "traces", "traces", ev(2),
                        yield_type="qualitative table result")]
    else:
        number = int(yraw.strip("<%"))
        value = q(number, "%", yraw, "lt" if yraw.startswith("<") else None)
        outs = [outcome("yield", "yield", pmid, yraw, value, ev(2), yield_type="reported product yield")]
    if er:
        a, b = map(int, er.split(":"))
        outs.append(outcome("er", "selectivity", pmid, er,
                            {"metric": "er", "major": a, "minor": b, "order": "major:minor"}, ev(2),
                            scope={"metric": "er", "order": "major:minor",
                                   "major_configuration": "S" if entry == 6 else "unassigned"}))
    reaction(f"table1_entry_{entry}", 2, categories=["optimization"],
             name=f"Table 1 entry {entry}: asymmetric hydroamination of {substrate}",
             transformation="Rh-catalyzed asymmetric intramolecular hydroamination",
             location=f"Page 2, Table 1, entry {entry}", participants=parts, steps=st, outcomes=outs,
             general=[G_ASYM], class_kind="C-N+H",
             notes=("Table footnote and page-6 general procedure disagree on ligand/acid loading: 10 versus 8 mol%.",),
             unresolved=("Ligand and Brønsted-acid loadings conflict between Table 1 and the general procedure.",
                         "Yield measurement basis is not stated.",
                         "Absolute configuration is unassigned except for the entry-6 major enantiomer."),
             evidence_extra=ev(6))


# Table 2 optimization.  Entry 5 is also the 2e scope arrow in Scheme 3, so its
# scope evidence is merged into the same record rather than counted twice.
table2 = [
    (1, "3d", "(R)-DTBM-Segphos", "rac-BNP acid", 50, "74:26", 50, "0.4 M"),
    (2, "3d", "(R)-DTBM-Garphos", "rac-BNP acid", 76, "75:25", 50, "0.4 M"),
    (3, "3d", "(R)-DTBM-Garphos", "rac-BNP acid", 70, "78:22", 50, "0.2 M"),
    (4, "3e", "(R)-DTBM-Segphos", "PPTS", 80, "90:10", 50, "0.4 M"),
    (5, "3e", "(R)-DTBM-Garphos", "PPTS", 70, "95:5", 50, "0.4 M"),
    (6, "3e", "(R)-DTBM-Garphos", "chloroacetic acid", 90, "91:9", 50, "0.4 M"),
    (7, "3e", "(R)-DTBM-Garphos", "PPTS", 60, "95:5", 70, "0.4 M"),
]
table2_records = {}
for entry, substrate, ligand, acid, yld, er, temp, conc in table2:
    assignment = "major enantiomer depicted by a solid wedge; absolute label not printed in the reviewed Table 2"
    pmid = mol(f"t2e{entry}_product", "2d" if substrate == "3d" else "2e", 3,
               chirality="mixture", geometry="unspecified", representation="not_encoded",
               assignment=assignment, mixture=f"enantiomeric mixture; er {er}")
    parts = hydro_participants(substrate, pmid, 3, ligand, acid,
                               q(10, "mol%", "10 mol%"), q(10, "mol%", "10 mol%"),
                               concentration=q(float(conc.split()[0]), "M", conc))
    st = [step("cyclization", ev(3), "Asymmetric allene cyclization under Table 2 conditions",
               temperature=q(temp, "°C", f"{temp} °C"), time=q(24, "h", "24 h"),
               participant_ids=[p["participant_id"] for p in parts[:-1]])]
    a, b = map(int, er.split(":"))
    outs = [outcome("yield", "yield", pmid, f"{yld}%", q(yld, "%", f"{yld}%"), ev(3),
                    yield_type="reported product yield"),
            outcome("er", "selectivity", pmid, er,
                    {"metric": "er", "major": a, "minor": b, "order": "major:minor"}, ev(3),
                    scope={"metric": "er", "order": "major:minor", "major_configuration": "unassigned"})]
    cats = ["optimization", "substrate_scope"] if entry == 5 else ["optimization"]
    r = reaction(f"table2_entry_{entry}", 3, categories=cats,
                 name=f"Table 2 entry {entry}: asymmetric hydroamination of {substrate}",
                 transformation="Rh-catalyzed asymmetric intramolecular hydroamination of allene",
                 location=f"Page 3, Table 2, entry {entry}" + (" and Scheme 3 2e arrow" if entry == 5 else ""),
                 participants=parts, steps=st, outcomes=outs, class_kind="C-N+H",
                 unresolved=("Yield measurement basis is not stated.",
                             "Absolute configuration of the major enantiomer is not printed in the reviewed region."))
    table2_records[entry] = r["reaction_id"]


# Scheme 3 scope entries not already represented by Table 2 entry 5.
scope = [
    ("3f", "2f", 86, "93:7", "chloroacetic acid"),
    ("3g", "2g", 55, "96:4", "chloroacetic acid"),
    ("3h", "2h", 99, "95:5", "chloroacetic acid"),
    ("3i", "2i", 86, "95:5", "chloroacetic acid"),
    ("3j", "2j", 70, "94:6", "chloroacetic acid"),
    ("3k", "2k", 72, "97:3", "chloroacetic acid"),
    ("3l", "2l", 60, "95:5", "chloroacetic acid"),
    ("3m", "2m", 86, "97:3", "chloroacetic acid"),
    ("3n", "2n", 78, "98:2", "chloroacetic acid"),
    ("3a", "2a", 90, "78:22", "PPTS"),
]
for substrate, product_label, yld, er, acid in scope:
    smid = mol("scope_" + substrate, substrate, 3,
               chirality="achiral", geometry="not_applicable", representation="unresolved",
               note="Source allene drawing inspected; full substituted structure was not encoded from the page image.")
    pmid = mol("scope_product_" + product_label, product_label, 3,
               chirality="mixture", geometry="unspecified", representation="not_encoded",
               assignment="major enantiomer depicted by solid wedge; absolute R/S label not printed in Scheme 3",
               mixture=f"enantiomeric mixture; er {er}")
    lmid, amid = agent("(R)-DTBM-Garphos", 3), agent(acid, 3)
    parts = [participant("substrate", smid, "reactant", ev(3)),
             participant("rh_catalyst", RH, "catalyst", ev(3), molpct=q(4, "mol%", "4 mol%"), ref=smid),
             participant("ligand", lmid, "ligand", ev(3), molpct=q(10, "mol%", "10 mol%"), ref=smid),
             participant("bronsted_acid", amid, "additive", ev(3), molpct=q(10, "mol%", "10 mol%"), ref=smid),
             participant("dce", DCE, "solvent", ev(3)),
             participant("product", pmid, "product", ev(3))]
    a, b = map(int, er.split(":"))
    outs = [outcome("yield", "yield", pmid, f"{yld}%", q(yld, "%", f"{yld}%"), ev(3),
                    yield_type="reported product yield"),
            outcome("er", "selectivity", pmid, er,
                    {"metric": "er", "major": a, "minor": b, "order": "major:minor"}, ev(3),
                    scope={"metric": "er", "order": "major:minor", "major_configuration": "unassigned"})]
    ckind = "C-O+H" if substrate == "3a" else "C-N+H"
    reaction(f"scheme3_{substrate}_{product_label}", 3, categories=["substrate_scope"],
             name=f"Scheme 3 cyclization of {substrate} to {product_label}",
             transformation="asymmetric hydroalkoxylation" if substrate == "3a" else "asymmetric hydroamination",
             location=f"Page 3, Scheme 3, {product_label} arrow", participants=parts,
             steps=[step("cyclization", ev(3), "Cyclize under Scheme 3 conditions",
                         temperature=q(50, "°C", "50 °C"), time=q(24, "h", "24 h"),
                         participant_ids=[p["participant_id"] for p in parts[:-1]])],
             outcomes=outs, class_kind=ckind,
             unresolved=("Full substrate/product SMILES were not encoded from the compact scope drawings.",
                         "Yield measurement basis is not stated.",
                         "Absolute configuration of the major enantiomer is not printed in Scheme 3."))


# Helpers for preparative experiments on pages 4--6.
def prep_mol(label, page):
    return mol(label, label, page)


def prep_reaction(suffix, page, start_label, product_label, yraw, *, categories=("other",),
                  name=None, transformation="preparative transformation", class_kind=None,
                  agents=(), steps=(), general=(), amount_raw=None, evidence_extra=(), notes=(), unresolved=()):
    evidence = ev(page)
    smid, pmid = prep_mol(start_label, page), prep_mol(product_label, page)
    parts = [participant("substrate", smid, "reactant", evidence),
             participant("product", pmid, "product", evidence)]
    for idx, item in enumerate(agents, 1):
        label, role, kwargs = item
        parts.append(participant(f"agent_{idx}", prep_mol(label, page), role, evidence, **kwargs))
    number = float(yraw.rstrip("%"))
    outs = [outcome("yield", "yield", pmid, yraw, q(number, "%", yraw), evidence + list(evidence_extra),
                    yield_type="reported isolated material after stated purification",
                    basis="isolated after the stated workup/purification")]
    if amount_raw:
        raw_mass, mass_num, raw_amount, amount_num = amount_raw
        outs.append(outcome("product_mass", "product_amount", pmid, raw_mass, q(mass_num, "g" if " g" in raw_mass else "mg", raw_mass), evidence,
                            basis="isolated material", scope={"reported_amount": raw_amount,
                                                               "amount": q(amount_num, "mmol", raw_amount)}))
    return reaction(suffix, page, categories=list(categories), name=name or f"Preparation of {product_label}",
                    transformation=transformation, location=f"Page {page}, Experimental Section, {product_label}",
                    participants=parts, steps=steps, outcomes=outs, general=general,
                    class_kind=class_kind, notes=notes,
                    unresolved=unresolved or ("Candidate structures remain machine-derived and need expert checking.",),
                    evidence_extra=evidence_extra)


# S1a protection.
s1a_agents = [
    ("Boc2O", "reagent", {"mass": q(9.3, "g", "9.3 g"), "amount": q(42, "mmol", "42 mmol"),
                            "eq": q(1.7, "equiv", "1.7 equiv"), "ref": prep_mol("ethyl 2-aminobenzoate", 4)}),
    ("DMAP", "catalyst", {"mass": q(0.92, "g", "0.92 g"), "amount": q(7.5, "mmol", "7.5 mmol"),
                           "eq": q(0.3, "equiv", "0.3 equiv"), "ref": prep_mol("ethyl 2-aminobenzoate", 4)}),
    ("Et3N", "reagent", {"volume": q(3.5, "mL", "3.5 mL"), "amount": q(25, "mmol", "25 mmol"),
                          "eq": q(1, "equiv", "1 equiv"), "ref": prep_mol("ethyl 2-aminobenzoate", 4)}),
    ("THF", "solvent", {"volume": q(250, "mL", "250 mL")}),
]
prep_reaction("prep_S1a", 4, "ethyl 2-aminobenzoate", "S1a", "70%",
              name="PG-amine protection: preparation of S1a", transformation="Boc protection", class_kind="C-N",
              agents=s1a_agents,
              steps=[step("reaction", ev(4), "Add Boc2O, DMAP, and Et3N to ester in dry THF; heat",
                          temperature=q(60, "°C", "60 °C"), time=q(24, "h", "24 h")),
                     step("workup", ev(4), "Quench with water; separate; extract aqueous layer with Et2O (3 × 50 mL); wash brine; dry Na2SO4; concentrate", stage="workup"),
                     step("purification", ev(4), "Silica gel, hexane/EtOAc 39:1", stage="purification")],
              amount_raw=("5.57 g", 5.57, "21 mmol", 21))


# N-alkylation products S2d, S2e and page-6 S10 series.
for suffix, page, start, product, y, mass, mmol, bromide in [
    ("prep_S2d", 4, "S1a", "S2d", 93, 2.95, 9.3, "but-2-yn-1-yl bromide"),
    ("prep_S2e", 4, "S1b", "S2e", 90, 3.34, 9.0, "but-2-yn-1-yl bromide"),
    ("prep_S10d", 6, "S1a", "S10d", 94, 2.85, 9.4, "propargyl bromide derivative"),
    ("prep_S10e", 6, "S1b", "S10e", 90, 3.2, 9.0, "propargyl bromide derivative"),
    ("prep_S10f", 6, "S1f", "S10f", 90, 3.21, 9.0, "propargyl bromide derivative"),
]:
    ref = prep_mol(start, page)
    agents = [("NaH", "reagent", {"eq": q(1.2, "equiv", "1.2 equiv"), "ref": ref}),
              (bromide, "reagent", {"eq": q(1.3, "equiv", "1.3 equiv"), "ref": ref}),
              ("DMF", "solvent", {})]
    prep_reaction(suffix, page, start, product, f"{y}%", name=f"N-alkylation to {product}",
                  transformation="N-alkylation", class_kind="C-N", agents=agents, general=[G_NALK],
                  steps=[step("deprotonation", ev(4), "Add substrate solution slowly to NaH in DMF and hold",
                              temperature=q(0, "°C", "0 °C"), time=q(2, "h", "2 h"), atmosphere="argon"),
                         step("alkylation", ev(4), f"Add {bromide}; warm slowly to room temperature and stir",
                              temperature="rt", time=q(16, "h", "16 h"), atmosphere="argon"),
                         step("workup", ev(4), "Quench with saturated NH4Cl; extract EtOAc; wash saturated NH4Cl (3 × 100 mL); dry Na2SO4; concentrate", stage="workup"),
                         step("purification", ev(4), "Silica gel, hexanes/EtOAc 8:2", stage="purification")],
                  amount_raw=(f"{mass} g", mass, f"{mmol} mmol", mmol),
                  evidence_extra=ev(4) if page == 6 else (),
                  unresolved=("Starting materials S1b/S1f are stated to be literature-derived; their preparation is not counted here.",
                              "Named product structure remains a machine candidate pending expert review."))


# Ester reductions to S3d and 1a.
for suffix, start, product, y, mass, mmol in [
    ("prep_S3d", "S2d", "S3d", 99, 2.53, 9.2),
    ("prep_1a", "S2e", "1a", 91, 2.8, 8.5),
]:
    ref = prep_mol(start, 5)
    prep_reaction(suffix, 5, start, product, f"{y}%", transformation="ester reduction to benzyl alcohol",
                  class_kind="C-H", general=[G_REDUCE],
                  agents=[("DIBAL-H", "reagent", {"eq": q(2.2, "equiv", "2.2 equiv"), "ref": ref}),
                          ("DCM", "solvent", {}), ("MeOH", "workup_agent", {"volume": q(5, "mL", "5 mL")})],
                  steps=[step("addition", ev(4), "Add 1 M DIBAL-H in DCM dropwise to ester in DCM (0.3 M)",
                              temperature=q(-78, "°C", "−78 °C")),
                         step("low_temperature_hold", ev(4), "Stir at −78 °C",
                              temperature=q(-78, "°C", "−78 °C"), time=q(3, "h", "3 h")),
                         step("quench", ev(4), "Add MeOH then saturated Rochelle salt at −78 °C",
                              temperature=q(-78, "°C", "−78 °C"), stage="workup"),
                         step("warm_workup", ev(4), "Warm to room temperature and stir",
                              temperature="rt", time=q(1, "h", "1 h"), stage="workup"),
                         step("purification", ev(4), "Extract DCM (3 × 30 mL), wash saturated NaCl, dry Na2SO4; silica hexanes/EtOAc 9:1 to 7:3", stage="purification")],
                  amount_raw=(f"{mass} g", mass, f"{mmol} mmol", mmol), evidence_extra=ev(4))


# Oxidation/reductive-amination sequences to 1d/1e; yield spans the stated start and final product.
for suffix, start, product, y, mass, mmol in [
    ("prep_1d", "S3d", "1d", 80, None, None),
    ("prep_1e", "1a", "1e", 70, 2.0, 4.8),
]:
    ref = prep_mol(start, 5)
    amount = None if mass is None else (f"{mass:g} g", mass, f"{mmol} mmol", mmol)
    prep_reaction(suffix, 5, start, product, f"{y}%", transformation="alcohol oxidation then reductive amination",
                  class_kind="C-N", general=[G_OXRA],
                  agents=[("DMP", "reagent", {"eq": q(1.1, "equiv", "1.1 equiv"), "ref": ref}),
                          ("p-anisidine", "reagent", {"eq": q(1.5, "equiv", "1.5 equiv"), "ref": ref}),
                          ("NaBH4", "reagent", {"eq": q(1.1, "equiv", "1.1 equiv"), "ref": ref}),
                          ("DCM", "solvent", {}), ("MeOH", "solvent", {})],
                  steps=[step("oxidation", ev(5), "Add DMP to alcohol in DCM (0.25 M); stir",
                              temperature="rt", time=q(30, "min", "30 min")),
                         step("oxidation_workup", ev(5), "Quench 1 M NaOH; extract DCM (2 × 30 mL), wash brine, dry Na2SO4, concentrate; use crude aldehyde", stage="workup"),
                         step("imine_formation", ev(5), "Add p-anisidine in MeOH (0.25 M) or 1:1 MeOH/DCM; stir",
                              temperature="rt", time=q(18, "h", "18 h"), atmosphere="argon"),
                         step("reduction_addition", ev(5), "Cool and add NaBH4 portionwise",
                              temperature=q(0, "°C", "0 °C")),
                         step("reduction_hold", ev(5), "Warm to room temperature and stir",
                              temperature="rt", time=q(2, "h", "2 h")),
                         step("purification", ev(5), "Quench water; extract DCM, wash brine, dry Na2SO4; silica hexanes/EtOAc 9:1 to 8:2", stage="purification")],
                  amount_raw=amount,
                  notes=("Reported yield is retained for the stated alcohol-to-amine sequence and is not relabeled as a last-step-only yield.",))


# Alkyne 1b route.
ref_sal = prep_mol("salicylaldehyde", 5)
prep_reaction("prep_S4", 5, "salicylaldehyde", "S4", "83%", transformation="O-alkylation",
              class_kind="C-O",
              agents=[("K2CO3", "reagent", {"mass": q(1, "g", "1 g"), "amount": q(7.2, "mmol", "7.2 mmol"),
                                               "eq": q(1.2, "equiv", "1.2 equiv"), "ref": ref_sal}),
                      ("1-bromo-2-butyne", "reagent", {"volume": q(0.58, "mL", "0.58 mL"),
                                                       "amount": q(6.6, "mmol", "6.6 mmol"),
                                                       "eq": q(1.1, "equiv", "1.1 equiv"), "ref": ref_sal}),
                      ("DMF", "solvent", {"volume": q(3, "mL", "3 mL")})],
              steps=[step("reaction", ev(5), "Add salicylaldehyde then 1-bromo-2-butyne to K2CO3 in DMF; stir",
                          temperature="rt", time=q(16, "h", "16 h")),
                     step("workup", ev(5), "Quench water (10 mL); extract EtOAc (3 × 10 mL); wash water (3 × 10 mL), brine (2 × 10 mL); dry Na2SO4", stage="workup"),
                     step("purification", ev(5), "Silica hexanes/EtOAc 19:1", stage="purification")],
              amount_raw=("870 mg", 870, "5 mmol", 5),
              notes=("Heading names the product S4, but the procedure sentence says ‘give S6’; S4 is retained as the intended label and the source conflict is unresolved.",),
              unresolved=("Source label conflict: S4 in heading versus S6 in the procedure sentence.",))

for suffix, start, product, y, mass, mmol, page in [
    ("prep_1b", "S4", "1b", 70, 0.983, 3.5, 5),
    ("prep_1c", "S9", "1c", 85, 1.5, 5.4, 6),
]:
    ref = prep_mol(start, page)
    prep_reaction(suffix, page, start, product, f"{y}%", transformation="reductive amination",
                  class_kind="C-N", general=[G_OXRA],
                  agents=[("p-anisidine", "reagent", {"eq": q(1.5, "equiv", "1.5 equiv"), "ref": ref}),
                          ("NaBH4", "reagent", {"eq": q(1.1, "equiv", "1.1 equiv"), "ref": ref}),
                          ("MeOH", "solvent", {})],
                  steps=[step("imine_formation", ev(5), "Use the general reductive-amination procedure; stir with p-anisidine",
                              temperature="rt", time=q(18, "h", "18 h"), atmosphere="argon"),
                         step("reduction_addition", ev(5), "Cool and add NaBH4 portionwise", temperature=q(0, "°C", "0 °C")),
                         step("reduction_hold", ev(5), "Warm to room temperature and stir", temperature="rt", time=q(2, "h", "2 h")),
                         step("purification", ev(5), "Aqueous workup; silica hexanes/EtOAc 9:1 to 8:2", stage="purification")],
                  amount_raw=(f"{mass} g", mass, f"{mmol} mmol", mmol), evidence_extra=ev(5) if page == 6 else ())


# Alkyne 1c route: S6--S9.
ref_s5 = prep_mol("S5", 5)
prep_reaction("prep_S6", 5, "S5", "S6", "86%", transformation="Pd-catalyzed aldehyde-chain construction",
              class_kind="C-C",
              agents=[("allylic alcohol", "reactant", {"volume": q(3.18, "mL", "3.18 mL"),
                                                        "amount": q(46.8, "mmol", "46.8 mmol"),
                                                        "eq": q(2.5, "equiv", "2.5 equiv"), "ref": ref_s5}),
                      ("Pd(OAc)2", "catalyst", {"mass": q(168, "mg", "168 mg"),
                                                "amount": q(0.75, "mmol", "0.75 mmol"),
                                                "molpct": q(4, "mol%", "4 mol%"), "ref": ref_s5}),
                      ("NaHCO3", "reagent", {"mass": q(7.54, "g", "7.54 g"),
                                             "amount": q(89.8, "mmol", "89.8 mmol"),
                                             "eq": q(4.8, "equiv", "4.8 equiv"), "ref": ref_s5}),
                      ("Bu4NCl", "additive", {"mass": q(5.2, "g", "5.2 g"),
                                              "amount": q(18.7, "mmol", "18.7 mmol"),
                                              "eq": q(1, "equiv", "1 equiv"), "ref": ref_s5}),
                      ("DMF", "solvent", {"volume": q(35, "mL", "35 mL")})],
              steps=[step("reaction", ev(5), "Stir S5 and allylic alcohol with Pd(OAc)2, NaHCO3, Bu4NCl in DMF",
                          temperature=q(50, "°C", "50 °C"), time=q(4, "h", "4 h")),
                     step("workup", ev(5), "Filter through silica; add H2O (30 mL); extract EtOAc (2 × 50 mL); wash H2O (2 × 50 mL) and brine (50 mL); dry Na2SO4", stage="workup"),
                     step("purification", ev(5), "Silica hexanes/EtOAc 19:1", stage="purification")],
              amount_raw=("3.8 g", 3.8, "16.1 mmol", 16.1),
              unresolved=("S5 is literature-derived and its preparation is not counted as a current experiment.",
                          "Reaction-center mapping is incomplete."))

ref_s6 = prep_mol("S6", 5)
prep_reaction("prep_S7", 5, "S6", "S7", "51%", transformation="aldehyde-to-terminal-alkyne homologation",
              class_kind="C-C+TB-C-C",
              agents=[("nBuLi", "reagent", {"volume": q(7.73, "mL", "7.73 mL"),
                                            "amount": q(19.32, "mmol", "19.32 mmol"),
                                            "eq": q(1.2, "equiv", "1.2 equiv"), "ref": ref_s6}),
                      ("DIPA", "reagent", {"volume": q(2.71, "mL", "2.71 mL"),
                                           "amount": q(19.32, "mmol", "19.32 mmol"),
                                           "eq": q(1.2, "equiv", "1.2 equiv"), "ref": ref_s6}),
                      ("trimethylsilyldiazomethane", "reagent", {"volume": q(8.05, "mL", "8.05 mL"),
                                                                 "amount": q(16.1, "mmol", "16.1 mmol"),
                                                                 "eq": q(1, "equiv", "1 equiv"), "ref": ref_s6}),
                      ("THF", "solvent", {"volume": q(130, "mL", "130 mL")})],
              steps=[step("base_addition", ev(5), "Add nBuLi dropwise to DIPA in THF", temperature=q(-78, "°C", "−78 °C")),
                     step("base_warm", ev(5), "Warm to 0 °C and stir", temperature=q(0, "°C", "0 °C"), time=q(15, "min", "15 min")),
                     step("diazo_addition", ev(5), "Cool to −78 °C, add trimethylsilyldiazomethane, and stir", temperature=q(-78, "°C", "−78 °C"), time=q(30, "min", "30 min")),
                     step("aldehyde_addition", ev(5), "Add S6 solution in THF (33 mL) and stir", temperature=q(-78, "°C", "−78 °C"), time=q(1, "h", "1 h")),
                     step("reflux", ev(5), "Heat to reflux", temperature="reflux", time=q(3, "h", "3 h")),
                     step("purification", ev(5), "Quench H2O (100 mL); extract Et2O (2 × 60 mL); wash H2O (100 mL); dry; silica 99:1", stage="purification")],
              amount_raw=("1.92 g", 1.92, "8.26 mmol", 8.26))

ref_s7 = prep_mol("S7", 5)
prep_reaction("prep_S8", 5, "S7", "S8", "92%", transformation="terminal alkyne methylation", class_kind="C-C",
              agents=[("nBuLi", "reagent", {"volume": q(3.2, "mL", "3.2 mL"), "amount": q(8, "mmol", "8 mmol"),
                                            "eq": q(1.1, "equiv", "1.1 equiv"), "ref": ref_s7}),
                      ("MeI", "reagent", {"volume": q(2.27, "mL", "2.27 mL"), "amount": q(36.5, "mmol", "36.5 mmol"),
                                          "eq": q(5, "equiv", "5 equiv"), "ref": ref_s7}),
                      ("THF", "solvent", {})],
              steps=[step("deprotonation", ev(5), "Add nBuLi to S7 in THF (0.3 M) and stir", temperature=q(-78, "°C", "−78 °C"), time=q(50, "min", "50 min")),
                     step("methylation", ev(5), "Add MeI and stir at room temperature", temperature="rt", time=q(16, "h", "16 h")),
                     step("purification", ev(5), "Quench saturated NH4Cl; extract EtOAc; wash brine; dry; silica 99:1", stage="purification")],
              amount_raw=("1.65 g", 1.65, "6.7 mmol", 6.7))

ref_s8 = prep_mol("S8", 5)
prep_reaction("prep_S9", 5, "S8", "S9", "94%", transformation="acetal deprotection to aldehyde", class_kind="DB-C-O",
              agents=[("PPTS", "reagent", {"mass": q(505, "mg", "505 mg"), "amount": q(2.01, "mmol", "2.01 mmol"),
                                            "eq": q(0.3, "equiv", "0.3 equiv"), "ref": ref_s8}),
                      ("acetone", "solvent", {"volume": q(260, "mL", "260 mL")}),
                      ("water", "reagent", {"volume": q(7, "mL", "7 mL")})],
              steps=[step("deprotection", ev(5), "Heat S8 and PPTS in acetone/H2O to reflux", temperature="reflux", time=q(15, "h", "15 h")),
                     step("workup", ev(6), "Remove volatiles; dissolve in DCM (30 mL), wash H2O (30 mL), extract aqueous layer DCM (3 × 20 mL); brine; dry", stage="workup"),
                     step("purification", ev(6), "Silica hexanes/EtOAc 98:2", stage="purification")],
              amount_raw=("1.08 g", 1.08, "6.3 mmol", 6.3), evidence_extra=ev(6))


# Scheme 4 derivatizations.  Page-11 text is used only as necessary context to
# expose two internal source conflicts; page 11 is not included in bounded coverage.
m2e = mol("2e_deriv_start", "2e", 4, chirality="specified", geometry="unspecified",
          representation="not_encoded", assignment="enantioenriched 2e; Scheme 4 retains configuration")
m4 = mol("4", "4", 4, chirality="specified", geometry="unspecified", representation="not_encoded",
         assignment="configuration retained from 2e")
m5 = mol("5", "5", 4, chirality="specified", geometry="unspecified", representation="not_encoded",
         assignment="configuration retained from 2e")
m6 = mol("6", "6", 4, chirality="specified", geometry="unspecified", representation="not_encoded",
         assignment="configuration retained from 4")
m7 = mol("7", "7", 4, chirality="mixture", geometry="unspecified", representation="not_encoded",
         assignment="major configuration retained through derivatization", mixture="94:6 er in Scheme 4")

can_mid, hcl_mid = prep_mol("CAN", 4), prep_mol("HCl", 4)
can_conf = conflict([
    f("2.5 equiv (Scheme 4 and prose parenthetical)", q(2.5, "equiv", "2.5 equiv"), ev(4) + [f"{DOC}:page11:text"]),
    f("275 mg, 0.5 mmol for 0.25 mmol 2e = 2.0 equiv", q(2.0, "equiv", "2.0 equiv"), [f"{DOC}:page11:text"],
      method="inference", derivation="0.5 mmol / 0.25 mmol"),
], ev(4) + [f"{DOC}:page11:text"])
p4parts = [participant("substrate", m2e, "reactant", ev(4)),
           participant("CAN", can_mid, "reagent", ev(4) + [f"{DOC}:page11:text"]),
           participant("HCl", hcl_mid, "workup_agent", ev(4)),
           participant("product", m4, "product", ev(4))]
p4parts[1]["charges"]["equivalents"] = can_conf
p4parts[1]["equivalence_reference"] = f(m2e, m2e, ev(4))
reaction("scheme4_2e_to_4", 4, categories=["derivatization"], name="Oxidative PMP deprotection of 2e",
         transformation="oxidative N-deprotection", location="Page 4 Scheme 4; page 11 detailed procedure used as context",
         participants=p4parts,
         steps=[step("CAN_addition", ev(4) + [f"{DOC}:page11:text"], "Add aqueous CAN dropwise to 2e in MeCN in an ice bath", temperature=q(0, "°C", "ice bath")),
                step("oxidation", [f"{DOC}:page11:text"], "Warm to 25 °C and stir", temperature=q(25, "°C", "25 °C"), time=q(3, "h", "3 h")),
                step("salt_formation", ev(4), "Work up, then treat with 1 M HCl in Et2O", stage="workup")],
         outcomes=[outcome("yield", "yield", m4, "85%", q(85, "%", "85%"), ev(4) + [f"{DOC}:page11:text"],
                           yield_type="reported isolated hydrochloride salt", basis="isolated material")],
         class_kind="N-H",
         notes=("CAN loading conflicts internally; both numerical alternatives are retained on the CAN participant charge rather than selecting one.",),
         unresolved=("CAN equivalents conflict within the source.", "Product stereostructure is not encoded."),
         evidence_extra=[f"{DOC}:page11:text"])

ref2e = m2e
reaction("scheme4_2e_to_5", 4, categories=["derivatization"], name="Reductive tosyl deprotection of 2e",
         transformation="N-detosylation", location="Page 4 Scheme 4",
         participants=[participant("substrate", ref2e, "reactant", ev(4)),
                       participant("sodium", prep_mol("sodium metal", 4), "reagent", ev(4), eq=q(6, "equiv", "6 equiv"), ref=ref2e),
                       participant("naphthalene", prep_mol("naphthalene", 4), "additive", ev(4), eq=q(0.2, "equiv", "0.2 equiv"), ref=ref2e),
                       participant("THF", prep_mol("THF", 4), "solvent", ev(4)),
                       participant("product", m5, "product", ev(4))],
         steps=[step("detosylation", ev(4), "Treat with sodium and naphthalene in THF", temperature="rt", time=q(16, "h", "16 h"))],
         outcomes=[outcome("yield", "yield", m5, "85%", q(85, "%", "85%"), ev(4), yield_type="reported product yield")],
         class_kind="N-H", unresolved=("Yield measurement basis is not stated in Scheme 4.", "Product stereostructure is not encoded."))

et3n_conf = conflict([
    f("2 equiv (Scheme 4)", q(2, "equiv", "2 equiv"), ev(4)),
    f("4 equiv (60 μL, 0.5 mmol for 0.123 mmol 4; detailed procedure)", q(4, "equiv", "4 equiv"), [f"{DOC}:page11:text"]),
], ev(4) + [f"{DOC}:page11:text"])
p6parts = [participant("substrate", m4, "reactant", ev(4)),
           participant("acryloyl_chloride", prep_mol("acryloyl chloride", 4), "reagent", ev(4), eq=q(2, "equiv", "2 equiv"), ref=m4),
           participant("Et3N", prep_mol("Et3N", 4), "reagent", ev(4) + [f"{DOC}:page11:text"]),
           participant("DMAP", prep_mol("DMAP", 4), "catalyst", ev(4), eq=q(0.1, "equiv", "0.1 equiv"), ref=m4),
           participant("DCM", prep_mol("DCM", 4), "solvent", ev(4)),
           participant("product", m6, "product", ev(4))]
p6parts[2]["charges"]["equivalents"] = et3n_conf
p6parts[2]["equivalence_reference"] = f(m4, m4, ev(4))
reaction("scheme4_4_to_6", 4, categories=["derivatization"], name="Acryloylation of 4",
         transformation="amide formation", location="Page 4 Scheme 4; page 11 detailed procedure used as context",
         participants=p6parts,
         steps=[step("addition", ev(4), "Add acryloyl chloride to 4, Et3N and DMAP in DCM at 0 °C", temperature=q(0, "°C", "0 °C")),
                step("warm_hold", ev(4), "Warm to room temperature and stir", temperature="rt", time=q(2, "h", "2 h"))],
         outcomes=[outcome("yield", "yield", m6, "65%", q(65, "%", "65%"), ev(4) + [f"{DOC}:page11:text"],
                           yield_type="reported isolated product", basis="isolated after silica chromatography")],
         class_kind="C-N",
         notes=("Et3N loading conflicts between Scheme 4 and the detailed procedure; both numerical alternatives are retained on the Et3N participant charge.",),
         unresolved=("Et3N loading conflict between Scheme 4 and detailed procedure.", "Product stereostructure is not encoded."),
         evidence_extra=[f"{DOC}:page11:text"])

hgc = prep_mol("Hoveyda-Grubbs second-generation catalyst", 4)
reaction("scheme4_6_to_7", 4, categories=["derivatization"], name="Ring-closing metathesis to 7",
         transformation="ring-closing metathesis", location="Page 4 Scheme 4; page 11 detailed procedure used as context",
         participants=[participant("scheme_substrate_6", m6, "reactant", ev(4)),
                       participant("prose_substrate_5", m5, "unknown", [f"{DOC}:page11:text"]),
                       participant("catalyst", hgc, "catalyst", ev(4), molpct=q(10, "mol%", "10 mol%"), ref=m6),
                       participant("DCM", prep_mol("DCM", 4), "solvent", ev(4)),
                       participant("product", m7, "product", ev(4))],
         steps=[step("metathesis", ev(4), "Reflux in DCM with Hoveyda-Grubbs second-generation catalyst",
                     temperature="reflux", time=q(36, "h", "36 h"))],
         outcomes=[outcome("yield", "yield", m7, "87%", q(87, "%", "87%"), ev(4) + [f"{DOC}:page11:text"],
                           yield_type="reported isolated product", basis="isolated after silica chromatography"),
                   outcome("er", "selectivity", m7, "94:6 er", {"metric": "er", "major": 94, "minor": 6,
                                                                   "order": "major:minor"}, ev(4),
                           scope={"metric": "er", "order": "major:minor", "major_configuration": "retained/unassigned"})],
         class_kind="C-C+DB-C-C",
         notes=("Scheme 4 shows 6→7, while the page-11 detailed RCM prose says a solution of 5 was added; the substrate identity conflict is retained.",),
         unresolved=("RCM substrate identity conflicts: 6 in Scheme 4 versus 5 in the detailed prose.",
                     "Product stereostructure is not encoded."), evidence_extra=[f"{DOC}:page11:text"])


# General context and characterization are deliberately not synthetic records.
bundle["source_context"] = [
    {"context_id": f"{DOC}:context_page1", "document_id": DOC, "pages": [1],
     "kind": "background_and_graphical_abstract",
     "evidence_ids": ev(1),
     "note": "Title, abstract, graphical abstract, Figure 1 bioactive examples, and cited background; no current specific experiment."},
    {"context_id": f"{DOC}:context_scheme1_prior_work", "document_id": DOC, "pages": [2],
     "kind": "literature_examples",
     "evidence_ids": ev(2),
     "note": "Scheme 1 panels a-c are cited prior work and were not counted as experiments in this paper."},
    {"context_id": f"{DOC}:context_characterization_pages4_6", "document_id": DOC, "pages": [4, 5, 6],
     "kind": "characterization",
     "evidence_ids": ev(4) + ev(5) + ev(6),
     "note": "NMR, MS/HRMS, optical/instrumental general information and spectral peak lists support compound descriptions but are not separate syntheses."},
    {"context_id": f"{DOC}:context_literature_starting_materials", "document_id": DOC, "pages": [4, 5, 6],
     "kind": "literature_preparations",
     "evidence_ids": ev(4) + ev(5) + ev(6),
     "note": "S1b, S1f and S5 are stated to be synthesized according to literature; no current-paper preparation record was invented."},
]

bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions

# Save raw stereo diagnostics for every concrete SMILES candidate.
diagnostics = {}
for m in bundle["molecules"]:
    smi = sr.value(m["smiles"])
    if smi:
        diagnostics[m["molecule_id"]] = sr.stereo_inspect(smi)

(RUN / "raw_model").mkdir(parents=True, exist_ok=True)
sr.write(RUN / "raw_model/pages_0001_0006_stereo_diagnostics.json", diagnostics)
sr.write(RUN / "submitted.json", bundle)

# Bounded per-page coverage.  Record IDs are drawn from this submission only.
by_page = {i: [] for i in range(1, 7)}
for r in reactions:
    for eid in r["evidence_ids"]:
        for page in range(1, 7):
            if f":page{page}:" in eid:
                by_page[page].append(r["reaction_id"])
for page in by_page:
    by_page[page] = list(dict.fromkeys(by_page[page]))

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title, abstract and graphical abstract", "finding": "current-work generic transformation only; no specific experiment", "evidence_ids": ev(1)},
            {"region": "Figure 1 and introductory text", "finding": "bioactive examples and cited background, not experiments in this paper", "evidence_ids": ev(1)},
        ], "record_ids": by_page[1], "unresolved": ["Small structures in Figure 1 were not registered because they are contextual examples, not current experiments."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Scheme 1", "finding": "three cited literature panels plus one current-work generic arrow; kept separate", "evidence_ids": ev(2)},
            {"region": "Scheme 2", "finding": "four racemic alkyne cyclizations; 1c also reports 2c′ as a 20% isomer mixture", "evidence_ids": ev(2)},
            {"region": "Table 1 and footnotes", "finding": "eight independent optimization entries", "evidence_ids": ev(2)},
            {"region": "optimization prose", "finding": "discussion supports table rows; no extra experiment beyond the table", "evidence_ids": ev(2)},
        ], "record_ids": by_page[2], "unresolved": ["Yield basis is not stated.", "Table 1 loadings conflict with the page-6 general procedure (10 versus 8 mol%)."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Table 2 and footnotes", "finding": "seven independent allene optimization entries", "evidence_ids": ev(3)},
            {"region": "Scheme 3", "finding": "ten additional scope arrows; the 2e arrow duplicates Table 2 entry 5 and was merged as evidence", "evidence_ids": ev(3)},
            {"region": "mechanistic discussion", "finding": "author-proposed inner/outer-sphere alternatives, not a mechanistic experiment", "evidence_ids": ev(3)},
        ], "record_ids": by_page[3], "unresolved": ["Compact scope structures were not fully encoded.", "Absolute R/S labels are not printed in the reviewed Scheme 3 region."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "Scheme 4", "finding": "four derivatization arrows with yields and conditions", "evidence_ids": ev(4)},
            {"region": "General Information", "finding": "instrumental/analytical context, not separate syntheses", "evidence_ids": ev(4)},
            {"region": "PG protection and N-alkylation", "finding": "S1a preparation, general N-alkylation, and S2d/S2e product records", "evidence_ids": ev(4)},
            {"region": "ester reduction procedure", "finding": "general multistage procedure continued on page 5", "evidence_ids": ev(4)},
        ], "record_ids": by_page[4], "unresolved": ["Scheme 4 CAN and Et3N loadings conflict with detailed page-11 prose; both values retained."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "S3d/1a and oxidation/reductive amination", "finding": "two ester reductions and two oxidation/reductive-amination sequences", "evidence_ids": ev(5)},
            {"region": "S4/1b route", "finding": "O-alkylation and reductive amination; S4/S6 source-label conflict retained", "evidence_ids": ev(5)},
            {"region": "S6-S9 route", "finding": "S6, S7, S8 and start of S9 preparations", "evidence_ids": ev(5)},
            {"region": "compound characterization paragraphs", "finding": "spectral data are source context, not separate synthetic experiments", "evidence_ids": ev(5)},
        ], "record_ids": by_page[5], "unresolved": ["S9 procedure continues onto page 6.", "S5 is literature-derived; no current preparation was inferred."]},
        {"page": 6, "state": "attempted", "regions": [
            {"region": "S9 continuation and 1c", "finding": "completion of S9 workup/purification and one reductive amination to 1c", "evidence_ids": ev(6)},
            {"region": "racemic/asymmetric alkyne cyclization procedures", "finding": "general procedures; Scheme 2 products are duplicate evidence, not new records", "evidence_ids": ev(6)},
            {"region": "S10d/S10e/S10f", "finding": "three N-alkylation product preparations", "evidence_ids": ev(6)},
            {"region": "product characterization", "finding": "spectral peak lists retained as context only", "evidence_ids": ev(6)},
        ], "record_ids": by_page[6], "unresolved": ["Page ends during S10f characterization; later pages were not claimed as covered."]},
    ],
}
(RUN / "coverage").mkdir(parents=True, exist_ok=True)
sr.write(RUN / "coverage/pages_0001_0006.json", coverage)

# Literal cross-region mention hints are raw review aids, not merged facts.
mentions = sr.source_mentions(DATASET, ["1a", "1d", "1e", "2a", "2e", "3d", "3e", "S4", "S9", "Scheme 4"], 140)
sr.write(RUN / "raw_model/pages_0001_0006_mentions.json", mentions)

print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
                  "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in reactions),
                  "generic_or_mechanistic": sum(r["record_type"] != "specific_experiment" for r in reactions)},
                 ensure_ascii=False, indent=2))
