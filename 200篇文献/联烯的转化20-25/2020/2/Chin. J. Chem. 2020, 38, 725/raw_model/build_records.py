#!/usr/bin/env python3
"""Build page-1--6 source-derived v3.1 records for doc_bd193f10bf2e3fd4.

This file is provenance/build logic, not a source-reviewed chemistry assertion.
Complex drawing-defined structures are deliberately left unresolved rather than
transcribed into unsafe exact SMILES.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_bd193f10bf2e3fd4"
RUN = ROOT / "runs/current"
BUNDLE = st.read(RUN / "submitted.json")


def ev(page: int, image: bool = True) -> list[str]:
    ids = [f"{DOC}:page{page}:text"]
    if image:
        ids.append(f"{DOC}:page{page}:image")
    return ids


def evs(*pages: int) -> list[str]:
    out = []
    for page in pages:
        out.extend(ev(page))
    return list(dict.fromkeys(out))


def F(raw=None, normalized=None, pages=(), method="original_extraction",
      status=None, derivation=None, alternatives=()):
    evidence = evs(*pages) if pages else []
    if normalized is None and raw is not None and not isinstance(raw, str):
        normalized = raw
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   derivation=derivation, alternatives=alternatives)


def Q(number, unit, *, qualifier=None, raw_text=None):
    value = {"number": number, "unit": unit}
    if qualifier:
        value["qualifier"] = qualifier
    if raw_text:
        value["raw_text"] = raw_text
    return value


def QR(minimum, maximum, unit, raw_text):
    return {"minimum": minimum, "maximum": maximum, "unit": unit, "raw_text": raw_text}


MOLECULES = {}
STEREO_DIAGNOSTICS = {}


def add_molecule(mid, label, pages, name=None, *, smiles=None,
                 chirality="unresolved", geometry="unresolved",
                 representation="unresolved", assignment=None,
                 atom_centers=None, double_bonds=None, axial=None, mixture=None,
                 note=None):
    evidence = evs(*pages)
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"] = [{
        "document_id": DOC,
        "compound_label": st.fact(label, label, evidence),
        "evidence_ids": evidence,
    }]
    original = name or label
    m["names"]["original"] = st.fact(original, original, evidence)
    m["names"]["english"] = st.fact(original, original, evidence)
    m["names"]["abbreviations"] = st.fact()
    m["names"]["aliases"] = st.fact()
    m["cas_numbers"] = st.fact()
    m["generic_template"] = None
    if smiles:
        diag = st.stereo_inspect(smiles)
        STEREO_DIAGNOSTICS[mid] = diag
        m["representation_type"] = "smiles"
        m["smiles"] = st.fact(
            f"structure transcribed from the source label/drawing: {label}", smiles,
            evidence, method="structure_derivation",
            derivation="Candidate connectivity was transcribed from the supplied page; RDKit parsing is not source verification.")
        iso = diag["canonical_isomeric_smiles"]
        m["stereochemistry"]["isomeric_smiles"] = st.fact(
            None, iso, evidence, method="structure_derivation",
            derivation="Canonical isomeric SMILES calculated from the unreviewed candidate SMILES.")
        m["stereochemistry"]["identity_key"] = st.fact(
            None, diag["stereo_graph_key"], evidence, method="structure_derivation",
            derivation="RDKit stereo graph key; not a sample identity or source verification.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = st.fact()
    def sf(key, value, raw=None, derived=False):
        if value is None:
            m["stereochemistry"][key] = st.fact()
        elif derived:
            m["stereochemistry"][key] = st.fact(
                raw or value, value, evidence, method="structure_derivation",
                derivation="Status inferred conservatively from the source drawing/name; exact structure remains unverified.")
        else:
            m["stereochemistry"][key] = st.fact(raw or value, value, evidence)
    sf("atom_centers", atom_centers)
    sf("double_bonds", double_bonds)
    sf("axial", axial)
    sf("mixture", mixture)
    sf("assignment", assignment)
    sf("chirality_status", chirality, derived=chirality in {"achiral", "not_applicable"})
    sf("geometry_status", geometry, derived=geometry == "not_applicable")
    sf("representation_status", representation, derived=representation == "not_applicable")
    m["identity_note"] = st.fact(
        note or "Source-labelled identity retained; exact drawing transcription was not independently verified.",
        note or "Source-labelled identity retained; exact drawing transcription was not independently verified.",
        evidence)
    m["external_properties"] = []
    MOLECULES[mid] = m
    return mid


# Source compounds and drawing-defined intermediates.
complex_specs = [
    ("tmpl_allenic_alcohol_1", "1", (1,), "generic allenic alcohol 1", "unspecified", "not_applicable", "not_encoded"),
    ("tmpl_diols_2", "2", (1,), "stereodefined various 1,3-diols 2", "specified", "unspecified", "not_encoded"),
    ("tmpl_thps_3", "3", (1,), "stereodefined various tetrahydropyrans 3", "specified", "unspecified", "not_encoded"),
    ("cmp6", "6", (2, 3), "phosphonium salt 6", "unspecified", "specified", "not_encoded"),
    ("cmp7", "7", (2, 3), "A-ring primary alcohol 7", "specified", "not_applicable", "not_encoded"),
    ("cmp10", "10", (2, 5, 6), "B-ring diester 10", "specified", "specified", "not_encoded"),
    ("cmp11_R", "(R)-11", (2, 3), "(R)-allenic alcohol 11 (>99% ee)", "specified", "not_applicable", "not_encoded"),
    ("cmp12_S", "(S)-12", (2, 5), "(S)-allenic acetate 12 (96% ee)", "specified", "not_applicable", "not_encoded"),
    ("cmp14", "14", (2, 3), "anti,anti-1,3-diol 14", "specified", "not_applicable", "not_encoded"),
    ("cmp15", "15", (2, 3), "primary alcohol 15", "specified", "not_applicable", "not_encoded"),
    ("cmp16", "16", (2, 3), "E-nitro olefin 16", "specified", "specified", "not_encoded"),
    ("cmp17", "17", (3,), "tetrahydropyran acetal 17", "specified", "not_applicable", "not_encoded"),
    ("cmp18", "18", (3,), "A-ring acetal 18", "specified", "not_applicable", "not_encoded"),
    ("cmp19", "19", (3,), "silyl-migrated product 19", "specified", "not_applicable", "not_encoded"),
    ("cmp20", "20", (3,), "borinate intermediate 20", "specified", "not_applicable", "not_encoded"),
    ("cmp21", "21", (3, 6), "top fragment 21", "specified", "specified", "not_encoded"),
    ("cmp22_rac_Z", "rac-(Z)-22", (4,), "racemic (Z)-1,4-diol 22", "racemic", "specified", "not_encoded"),
    ("cmp22_rac_E", "rac-(E)-22", (4,), "racemic (E)-diol 22", "racemic", "specified", "not_encoded"),
    ("cmp22_S_Z", "(S)-(Z)-22", (4, 5), "(S)-(Z)-1,4-diol 22", "specified", "specified", "not_encoded"),
    ("cmp22_S_E", "(S)-(E)-22", (4,), "(S)-(E)-diol 22", "specified", "specified", "not_encoded"),
    ("cmp23a_rac", "rac-23a", (4,), "racemic syn,syn-triol 23a", "racemic", "not_applicable", "not_encoded"),
    ("cmp23b_rac", "rac-23b", (4,), "racemic anti,syn-triol 23b", "racemic", "not_applicable", "not_encoded"),
    ("cmp23c_rac", "rac-23c", (4,), "racemic diol 23c", "racemic", "not_applicable", "not_encoded"),
    ("cmp23a_S", "23a", (4, 5), "enantioenriched syn,syn-triol 23a", "specified", "not_applicable", "not_encoded"),
    ("cmp23b_S", "23b", (4,), "enantioenriched anti,syn-triol 23b", "specified", "not_applicable", "not_encoded"),
    ("cmp23c_S", "23c", (4,), "enantioenriched diol 23c", "specified", "not_applicable", "not_encoded"),
    ("cmp24", "24", (5,), "diol 24", "specified", "not_applicable", "not_encoded"),
    ("cmp25", "25", (5,), "Michael acceptor 25", "specified", "specified", "not_encoded"),
    ("cmp26", "26", (5, 6), "cis-tetrahydropyran 26", "specified", "not_applicable", "not_encoded"),
    ("cmp28", "28", (5, 6), "aldehyde intermediate 28", "specified", "not_applicable", "not_encoded"),
    ("cmp29_11S_Z", "29 (11S,Z)", (5, 6), "desired Z-olefin 29, 11S", "specified", "specified", "not_encoded"),
    ("cmp29_11R_Z", "29 (11R,Z)", (5,), "C11 epimer of Z-olefin 29", "specified", "specified", "not_encoded"),
    ("cmp29_11S_E", "29 (11S,E)", (6,), "E-olefin isomer of 29", "specified", "specified", "not_encoded"),
    ("cmp30", "30", (6,), "aldehyde/allylic alcohol 30", "specified", "specified", "not_encoded"),
    ("cmp31a_EE", "31a (E,E major)", (6,), "diene ethyl ester 31a, E,E major isomer", "specified", "specified", "not_encoded"),
    ("cmp31a_other", "31a (minor geometrical isomer)", (6,), "minor geometrical isomer of 31a", "specified", "unspecified", "not_encoded"),
    ("cmp31b_EE", "31b (E,E major)", (6,), "diene methyl ester 31b, E,E major isomer", "specified", "specified", "not_encoded"),
    ("cmp31b_other", "31b (minor geometrical isomer)", (6,), "minor geometrical isomer of 31b", "specified", "unspecified", "not_encoded"),
    ("cmp32a", "32a", (6,), "enal 32a (R = Et)", "specified", "specified", "not_encoded"),
    ("cmp32b", "32b", (6,), "enal 32b (R = Me)", "specified", "specified", "not_encoded"),
    ("cmp33", "33", (6,), "aldehyde 33", "specified", "specified", "not_encoded"),
    ("cmp34a", "34a", (6,), "Kang's sulfone-thioether linker 34a", "unresolved", "not_applicable", "not_encoded"),
    ("cmp35_E", "35 (E major)", (6,), "E-alkene 35", "specified", "specified", "not_encoded"),
    ("cmp35_Z", "35 (Z minor)", (6,), "Z-alkene isomer of 35", "specified", "specified", "not_encoded"),
    ("cmp36a", "36a", (6,), "BT-sulfone 36a (R = TMS)", "specified", "specified", "not_encoded"),
    ("cmp36b", "36b", (6,), "BT-sulfone 36b (R = H)", "specified", "specified", "not_encoded"),
]
for mid, label, pages, name, chir, geom, rep in complex_specs:
    assignment = None
    mixture = None
    atom_centers = None
    double_bonds = None
    if chir in {"specified", "racemic"}:
        atom_centers = f"Source drawing/name specifies stereochemical identity for {label}; exact R/S transcription not attempted."
    if chir == "racemic":
        mixture = "Source explicitly identifies the substrate/product series as racemic."
    if geom == "specified":
        double_bonds = f"Source drawing or label specifies alkene geometry for {label}."
    if label == "(R)-11": assignment = "R configuration explicitly stated; >99% ee."
    if label == "(S)-12": assignment = "S configuration explicitly stated; 96% ee."
    if label == "17": assignment = "C19/C23 cis relationship confirmed by author-reported NOE study."
    if label == "26": assignment = "Source identifies the tetrahydropyran as cis; dr reported per experiment."
    add_molecule(mid, label, pages, name, chirality=chir, geometry=geom,
                 representation=rep, assignment=assignment,
                 atom_centers=atom_centers, double_bonds=double_bonds,
                 mixture=mixture,
                 note="Complex drawing-defined structure retained by paper label; exact SMILES withheld to avoid unsafe stereochemical transcription.")

# A source-drawn, non-stereogenic aldehyde and simple, safely transcribed agents/solvents.
add_molecule("cmp13", "13", (2, 3), "3-(benzyloxy)propanal", smiles="O=CCCOCc1ccccc1",
             chirality="achiral", geometry="not_applicable", representation="not_applicable",
             note="Candidate connectivity transcribed from Scheme 4; parsing does not verify the drawing.")

simple = [
    ("solv_thf", "THF", (3, 4, 5, 6), "tetrahydrofuran", "C1CCOC1"),
    ("solv_dcm", "CH2Cl2", (3, 5, 6), "dichloromethane", "ClCCl"),
    ("solv_acetone", "acetone", (3,), "acetone", "CC(=O)C"),
    ("solv_toluene", "PhCH3", (3,), "toluene", "Cc1ccccc1"),
    ("solv_meoh", "MeOH", (3,), "methanol", "CO"),
    ("solv_ethanol", "EtOH", (5,), "ethanol", "CCO"),
    ("reag_nitromethane", "CH3NO2", (3,), "nitromethane", "C[N+](=O)[O-]"),
    ("reag_tea", "TEA", (3,), "triethylamine", "CCN(CC)CC"),
    ("reag_mscl", "MsCl", (3,), "methanesulfonyl chloride", "CS(=O)(=O)Cl"),
    ("reag_naoh", "NaOH", (4, 5), "sodium hydroxide", "[Na+].[OH-]"),
    ("reag_h2o2", "H2O2", (4, 5, 6), "hydrogen peroxide", "OO"),
    ("reag_h2", "H2", (3,), "hydrogen", "[H][H]"),
    ("reag_nh4cl", "NH4Cl", (6,), "ammonium chloride", "[NH4+].[Cl-]"),
    ("reag_nahco3", "NaHCO3", (6,), "sodium bicarbonate", "[Na+].O=C([O-])O"),
    ("reag_26lutidine", "2,6-lutidine", (6,), "2,6-lutidine", "Cc1cccc(C)n1"),
    ("reag_imidazole", "imid.", (3,), "imidazole", "c1ncc[nH]1"),
    ("reag_ozone", "O3", (3,), "ozone", "[O-][O+]=O"),
    ("reag_nabh4", "NaBH4", (3,), "sodium borohydride", "[Na+].[BH4-]"),
    ("reag_dmp_acetonide", "DMP (4.0 eq)", (5,), "2,2-dimethoxypropane", "COC(C)(C)OC"),
    ("reag_tbscl", "TBSCl", (3,), "tert-butyldimethylsilyl chloride", "CC(C)(C)[Si](C)(C)Cl"),
    ("reag_ptsa", "PTSA", (3, 5), "p-toluenesulfonic acid", "Cc1ccc(S(=O)(=O)O)cc1"),
    ("reag_pcsa", "PCSA", (5,), "p-chlorobenzenesulfonic acid", "O=S(=O)(O)c1ccc(Cl)cc1"),
    ("reag_pnsa", "PNSA", (5,), "p-nitrobenzenesulfonic acid", "O=S(=O)(O)c1ccc([N+](=O)[O-])cc1"),
]
for mid, label, pages, name, smiles in simple:
    add_molecule(mid, label, pages, name, smiles=smiles, chirality="achiral",
                 geometry="not_applicable", representation="not_applicable",
                 note="Simple source-named/drawn agent; candidate SMILES retained for graph diagnostics only.")

# Agents for which an exact formulation, counterion, stereochemical constitution, or source drawing
# was not safe to encode from these six pages.
unresolved_agents = [
    ("reag_sia2bh", "Sia2BH", (1, 2, 3, 4, 5), "disiamylborane"),
    ("reag_csa", "CSA", (2, 3, 5), "camphorsulfonic acid"),
    ("reag_tbaf", "TBAF", (2, 3, 5), "tetrabutylammonium fluoride"),
    ("reag_baib", "BAIB", (2, 3, 5, 6), "bis(acetoxy)iodobenzene"),
    ("reag_tempo", "TEMPO", (2, 3, 5), "TEMPO"),
    ("reag_tmg", "TMG", (2, 3), "tetramethylguanidine"),
    ("reag_pdc", "10% Pd/C", (3,), "10% palladium on carbon"),
    ("reag_li_naph", "lithium naphthalenide", (3,), "lithium naphthalenide"),
    ("reag_9bbn", "9-BBN", (3, 4, 5), "9-borabicyclo[3.3.1]nonane"),
    ("reag_khmds", "KHMDS", (3, 6), "potassium bis(trimethylsilyl)amide"),
    ("reag_thexbh2", "ThexBH2", (3, 4), "thexylborane"),
    ("reag_ipcbh2_plus", "(+)-IpcBH2", (4, 5), "(+)-isopinocampheylborane"),
    ("reag_ipcbh2_minus", "(–)-IpcBH2", (4,), "(–)-isopinocampheylborane"),
    ("reag_ylide6", "phosphonium salt 6 / KHMDS", (3,), "Wittig reagent generated from phosphonium salt 6"),
    ("reag_ylide_ethyl", "Ph3P=CHCO2Et", (5,), "ethyl(triphenylphosphoranylidene)acetate"),
    ("reag_tcsa", "TCSA", (5,), "trichlorobenzenesulfonic acid shown in Table 2"),
    ("reag_dmp_oxidant", "DMP", (5, 6), "Dess–Martin periodinane"),
    ("reag_tempobf4", "TEMPOBF4", (5,), "TEMPO oxoammonium tetrafluoroborate"),
    ("reag_tcca", "TCCA", (5,), "trichloroisocyanuric acid"),
    ("reag_27a", "27a", (5,), "4-substituted TEMPO derivative 27a"),
    ("reag_27b", "27b", (5,), "4-substituted TEMPO derivative 27b"),
    ("reag_27c", "27c", (5, 6), "4-benzyloxy-TEMPO"),
    ("reag_hwe_p5", "HWE reagent (identity not shown in Table 3)", (5,), "unspecified HWE reagent used after oxidation"),
    ("reag_still_gennari", "Still–Gennari phosphonate", (6,), "methyl bis(2,2,2-trifluoroethyl)phosphonoacetate"),
    ("reag_18c6", "18-c-6", (6,), "18-crown-6"),
    ("reag_tbsotf", "TBSOTf", (6,), "tert-butyldimethylsilyl triflate"),
    ("reag_dibalh", "DIBAL-H", (6,), "diisobutylaluminium hydride"),
    ("cmp9a", "9a", (2, 6), "HWE phosphonate 9a (R = Et)"),
    ("cmp9b", "9b", (2, 6), "HWE phosphonate 9b (R = Me)"),
    ("reag_lihmds", "LiHMDS", (6,), "lithium bis(trimethylsilyl)amide"),
    ("reag_tmsotf", "TMSOTf", (6,), "trimethylsilyl triflate"),
    ("reag_ammonium_molybdate", "(NH4)6Mo7O24·4H2O", (6,), "ammonium heptamolybdate tetrahydrate"),
]
for mid, label, pages, name in unresolved_agents:
    chir = "specified" if mid in {"reag_csa", "reag_ipcbh2_plus", "reag_ipcbh2_minus"} else "unresolved"
    assignment = "Source reagent name/sign identifies a chiral reagent; absolute encoding was not established." if chir == "specified" else None
    add_molecule(mid, label, pages, name, chirality=chir, geometry="not_applicable",
                 representation="not_encoded", assignment=assignment,
                 note="Agent identity retained from the source; exact formulation/structure was not safely encoded from pages 1–6.")


REACTIONS = []


def new_reaction(suffix, pages, *, record_type="specific_experiment", categories=("other",),
                 name=None, transformation=None, location=None):
    evidence = evs(*pages)
    r = st.new_reaction(f"{DOC}:{suffix}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = st.fact(list(categories), list(categories), evidence)
    loc = location or {"pages": list(pages)}
    r["location"] = st.fact(str(loc), loc, evidence)
    r["reaction_name"] = st.fact(name, name, evidence) if name else st.fact()
    if transformation:
        r["transformation_type"] = st.fact(
            None, transformation, evidence, method="inference",
            derivation="Transformation description inferred from the source arrow and adjacent prose.")
    r["general_scheme_ids"] = []
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    REACTIONS.append(r)
    return r


def add_participant(r, pid, mid, role, pages, *, stage="reaction", mass=None,
                    amount=None, volume=None, equiv=None, mol_percent=None,
                    ref=None, product_role=None):
    evidence = evs(*pages)
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    if product_role:
        p["product_role"] = st.fact(product_role, product_role, evidence)
    values = {"mass": mass, "amount": amount, "volume": volume,
              "equivalents": equiv, "mol_percent": mol_percent}
    for key, val in values.items():
        if val is not None:
            raw = val.get("raw_text") if isinstance(val, dict) and val.get("raw_text") else str(val)
            p["charges"][key] = st.fact(raw, val, evidence)
    if equiv is not None or mol_percent is not None:
        p["equivalence_reference"] = st.fact(ref, ref, evidence)
    r["participants"].append(p)
    return pid


def add_step(r, sid, pages, participant_ids, operation, *, temperature=None, time=None,
             atmosphere=None, pressure=None, stage="reaction"):
    evidence = evs(*pages)
    s = st.new_step(sid, evidence, stage=stage)
    s["participant_ids"] = list(participant_ids)
    s["operation"] = st.fact(operation, operation, evidence)
    if temperature is not None:
        s["temperature"] = st.fact(str(temperature), temperature, evidence)
    if time is not None:
        s["time"] = st.fact(str(time), time, evidence)
    if atmosphere is not None:
        s["atmosphere"] = st.fact(atmosphere, atmosphere, evidence)
    if pressure is not None:
        s["pressure"] = st.fact(pressure, pressure, evidence)
    r["steps"].append(s)
    return s


def set_solvent(r, pages, composition=None, *, volume=None, concentration=None, basis=None):
    evidence = evs(*pages)
    if composition is not None:
        r["solvent_conditions"]["composition"] = st.fact(str(composition), composition, evidence)
    if volume is not None:
        r["solvent_conditions"]["total_volume"] = st.fact(str(volume), volume, evidence)
    if concentration is not None:
        r["solvent_conditions"]["concentration"] = st.fact(str(concentration), concentration, evidence)
    if basis is not None:
        r["solvent_conditions"]["concentration_basis"] = st.fact(str(basis), basis, evidence)


def add_outcome(r, oid, kind, targets, pages, value, *, yield_type=None,
                basis=None, scope=None, raw=None, alternatives=()):
    evidence = evs(*pages)
    o = st.new_outcome(oid, kind, targets, evidence)
    raw_value = raw if raw is not None else str(value)
    if alternatives:
        o["value"] = st.fact(raw_value, None, evidence, status="conflict", alternatives=alternatives)
    else:
        o["value"] = st.fact(raw_value, value, evidence)
    if yield_type:
        o["yield_type"] = st.fact(yield_type, yield_type, evidence)
    if basis:
        o["measurement_basis"] = st.fact(basis, basis, evidence)
    if scope:
        o["scope"] = st.fact(str(scope), scope, evidence)
    r["outcomes"].append(o)
    return o


def add_note(r, text, pages):
    r["notes"].append(st.fact(text, text, evs(*pages)))


LABEL_META = {
    "C-H": (1, "C", "H", 1, "hydrogen_gain"),
    "C-C": (2, "C", "C", 1, "formed"),
    "C-O": (3, "C", "O", 1, "formed"),
    "DB-C-C": (13, "C", "C", 2, "formed"),
    "DB-C-O": (14, "C", "O", 2, "order_increased"),
    "DB-O-S": (18, "O", "S", 2, "formed"),
    "O-H": (24, "O", "H", 1, "hydrogen_gain"),
    "O-Si": (28, "O", "Si", 1, "formed"),
}


def classify_partial(r, labels, pages, limitation=None):
    evidence = evs(*pages)
    changes = []
    normalized_labels = []
    for idx, code in enumerate(labels, 1):
        number, a, b, order, kind = LABEL_META[code]
        before = 0 if kind in {"formed", "hydrogen_gain"} else 1
        change = {"change_id": f"b{idx}", "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
                  "before_order": before, "after_order": order, "change_type": kind}
        if kind == "hydrogen_gain": change["hydrogen_count"] = 1
        changes.append(st.fact(None, change, evidence, method="structure_derivation",
                               derivation="Partial bond-change class inferred from the source arrow; atoms are not fully mapped."))
        normalized_labels.append({"category_number": number, "code": code, "supporting_change_ids": [f"b{idx}"]})
    r["classification"]["bond_changes"] = changes
    r["classification"]["labels"] = st.fact(
        None, normalized_labels, evidence, method="structure_derivation",
        derivation="Taxonomy labels inferred from the visibly reported net transformation; no complete atom mapping claimed.")
    r["classification"]["coverage"] = st.fact(
        None, "partial", evidence, method="inference",
        derivation="Complex structures and/or atom correspondence remain unresolved.")
    r["classification"]["limitations"] = [limitation or "Partial source-drawing classification; complete atom mapping and all implicit-H changes are unresolved."]


def classify_na(r, pages, reason):
    evidence = evs(*pages)
    r["classification"]["labels"] = st.fact(None, [], evidence, method="inference", derivation=reason)
    r["classification"]["coverage"] = st.fact(None, "not_applicable", evidence, method="inference", derivation=reason)
    r["classification"]["limitations"] = [reason]


def classify_unresolved(r, pages, reason):
    evidence = evs(*pages)
    r["classification"]["coverage"] = st.fact(None, "unresolved", evidence, method="inference", derivation=reason)
    r["classification"]["limitations"] = [reason]


# Page 1 previous-work generic scheme.
r = new_reaction("scheme1_previous_generic", (1,), record_type="generic_scheme", name="hydroboration–allylation / oxa-Michael approach", transformation="generic stereodivergent tetrahydropyran construction", location={"page": 1, "region": "Scheme 1"})
for pid, mid, role in [("s1_allene", "tmpl_allenic_alcohol_1", "reactant"), ("s1_diols", "tmpl_diols_2", "product"), ("s1_thps", "tmpl_thps_3", "product"), ("s1_sia", "reag_sia2bh", "reagent")]:
    add_participant(r, pid, mid, role, (1,))
add_step(r, "generic_hydroboration_allylation", (1,), ["s1_allene", "s1_sia", "s1_diols"], "generic hydroboration–allylation; aldehyde and detailed conditions not specified in Scheme 1")
add_step(r, "generic_oxa_michael", (1,), ["s1_diols", "s1_thps"], "generic oxa-Michael sequence, 3–4 steps in total")
add_note(r, "Authors label this as their previous approach; it is not a specific current-paper experiment.", (1,))
classify_unresolved(r, (1,), "Generic R-group scheme lacks concrete atom identities and complete conditions.")

# Page 2 conceptual enantio-integrative approach.
r = new_reaction("scheme3_enantio_integrative_concept", (2,), record_type="generic_scheme", name="chirality-economic enantio-integrative approach", transformation="conceptual routing of both resolved allene enantiomers", location={"page": 2, "region": "Scheme 3"})
for pid, mid, role in [("s3_r11", "cmp11_R", "reactant"), ("s3_7", "cmp7", "product"), ("s3_s12", "cmp12_S", "reactant"), ("s3_10", "cmp10", "product")]:
    add_participant(r, pid, mid, role, (2,))
add_note(r, "Retrosynthetic/conceptual arrows only; individual reported experiments are recorded from Schemes 4–11 and Tables 1–3.", (2,))
classify_unresolved(r, (2,), "Conceptual retrosynthetic scheme is not an atom-mapped experiment.")

# Scheme 4, current A-ring synthesis.
r = new_reaction("scheme4_11_13_to_14", (2, 3), name="hydroboration–allylation", transformation="allene hydroboration and aldehyde allylation", location={"pages": [2, 3], "region": "Scheme 4 and adjacent prose"})
for args in [("a11", "cmp11_R", "reactant"), ("a13", "cmp13", "reactant"), ("a14", "cmp14", "product"), ("asia", "reag_sia2bh", "reagent"), ("athf", "solv_thf", "solvent")]: add_participant(r, *args, (2, 3))
set_solvent(r, (3,), "THF")
add_step(r, "hydroboration_allylation", (2, 3), ["a11", "a13", "asia", "athf", "a14"], "hydroboration of (R)-11 with Sia2BH followed by kinetic allylation of aldehyde 13")
add_outcome(r, "yield14", "yield", ["cmp14"], (3,), Q(98, "%"), basis="isolated", raw="98%")
add_outcome(r, "stereo14", "selectivity", ["cmp14"], (2,), "single stereoisomer", scope={"metric": "diastereoselectivity", "assignment": "anti,anti"})
r["review"]["unresolved"] = ["Exact quantities, temperature, time, workup and purification are not reported on pages 2–3."]
classify_partial(r, ["C-C", "C-H"], (2, 3))

r = new_reaction("scheme4_14_to_15", (2, 3), name="acetonide protection and BPS deprotection", transformation="diol protection followed by desilylation", location={"page": 3, "region": "Scheme 4"})
for args in [("b14", "cmp14", "reactant"), ("b15", "cmp15", "product"), ("bcsa", "reag_csa", "reagent"), ("btbaf", "reag_tbaf", "reagent"), ("bac", "solv_acetone", "solvent"), ("bdcm", "solv_dcm", "solvent")]:
    kw = {"equiv": Q(0.25, "equiv"), "ref": "cmp14"} if args[1] == "reag_csa" else {}
    add_participant(r, *args, (3,), **kw)
set_solvent(r, (3,), ["acetone", "CH2Cl2"])
add_step(r, "acetonide_protection", (3,), ["b14", "bcsa", "bac", "bdcm"], "treat with CSA (0.25 equiv) in acetone/CH2Cl2")
add_step(r, "desilylation", (3,), ["b14", "btbaf", "b15"], "then TBAF", temperature="r.t.")
add_outcome(r, "yield15", "yield", ["cmp15"], (2, 3), Q(100, "%"), basis="quantitative", raw="quant.")
r["review"]["unresolved"] = ["Step durations, reagent quantities other than CSA loading, workup and purification are unreported on pages 2–3."]
classify_partial(r, ["C-O"], (3,))

r = new_reaction("scheme4_15_to_16", (2, 3), name="one-pot oxidation/Henry reaction/elimination", transformation="primary alcohol to E-nitro olefin", location={"pages": [2, 3], "region": "Scheme 4"})
for args in [("c15", "cmp15", "reactant"), ("c16", "cmp16", "product"), ("cbaib", "reag_baib", "reagent"), ("ctempo", "reag_tempo", "reagent"), ("ctmg", "reag_tmg", "reagent"), ("cnitro", "reag_nitromethane", "reactant"), ("cmscl", "reag_mscl", "reagent"), ("ctea", "reag_tea", "reagent"), ("cdcm", "solv_dcm", "solvent"), ("ctol", "solv_toluene", "solvent")]: add_participant(r, *args, (2, 3))
set_solvent(r, (3,), ["CH2Cl2", "PhCH3"])
add_step(r, "oxidation", (2, 3), ["c15", "cbaib", "ctempo", "cdcm"], "BAIB/TEMPO oxidation in CH2Cl2")
add_step(r, "henry", (2, 3), ["ctmg", "cnitro", "ctol"], "TMG-promoted Henry reaction with CH3NO2 in PhCH3")
add_step(r, "mesylation_elimination", (2, 3), ["cmscl", "ctea", "c16"], "then MsCl and TEA; elimination affords E-nitro olefin 16")
add_outcome(r, "yield16", "yield", ["cmp16"], (2, 3), Q(81, "%"), raw="81% (2 steps)", scope={"span": "15 to 16", "reported_steps": 2})
r["review"]["unresolved"] = ["The source calls this an 81% two-step yield although the arrow lists several operations; exact charges, temperatures and times are unreported.", "Yield measurement basis is not stated on pages 2–3."]
classify_partial(r, ["DB-C-C"], (2, 3))

r = new_reaction("scheme4_16_to_17", (3,), name="oxa-Michael/Nef/acetalization cascade", transformation="nitroalkene cyclization to tetrahydropyran acetal", location={"page": 3, "region": "Scheme 4"})
for args in [("d16", "cmp16", "reactant"), ("d17", "cmp17", "product"), ("dptsa", "reag_ptsa", "reagent"), ("dmeoh", "solv_meoh", "solvent")]: add_participant(r, *args, (3,))
set_solvent(r, (3,), "MeOH")
add_step(r, "cascade", (3,), ["d16", "dptsa", "dmeoh", "d17"], "PTSA-promoted deketalization, oxa-Michael addition, in-situ Nef reaction and acetalization", temperature=Q(50, "°C"))
add_outcome(r, "yield17", "yield", ["cmp17"], (3,), Q(82, "%"), basis="isolated", raw="82%")
add_outcome(r, "dr17", "selectivity", ["cmp17"], (3,), {"ratio": ">20/1"}, raw="dr >20/1", scope={"metric": "dr", "major": "cis C19/C23"})
r["review"]["unresolved"] = ["Reaction time, reagent charge, workup, purification and dr measurement method are not reported on page 3."]
classify_partial(r, ["C-O"], (3,))

r = new_reaction("scheme4_17_to_18", (3,), name="ozonolysis/reduction/selective silylation", transformation="vinyl group to protected hydroxymethyl group", location={"page": 3, "region": "Scheme 4"})
for args in [("e17", "cmp17", "reactant"), ("e18", "cmp18", "product"), ("eo3", "reag_ozone", "reagent"), ("enabh4", "reag_nabh4", "reagent"), ("etbscl", "reag_tbscl", "reagent"), ("eimid", "reag_imidazole", "reagent")]: add_participant(r, *args, (3,))
add_step(r, "ozonolysis", (3,), ["e17", "eo3"], "ozonolysis", temperature=Q(-78, "°C"))
add_step(r, "reduction", (3,), ["enabh4"], "then NaBH4 reduction", temperature=Q(0, "°C"))
add_step(r, "silylation", (3,), ["etbscl", "eimid", "e18"], "then TBSCl/imidazole selective silylation", temperature="r.t.")
add_outcome(r, "yield18", "yield", ["cmp18"], (3,), Q(90, "%"), raw="90%")
r["review"]["unresolved"] = ["Yield measurement basis, durations, solvent, quantities, workup and purification are not reported on page 3."]
classify_partial(r, ["C-H", "O-H", "O-Si"], (3,))

# Schemes 5 and 6.
r = new_reaction("scheme5_hydrogenolysis_18_to_7", (3,), name="heterogeneous hydrogenation debenzylation", transformation="benzyl ether cleavage", location={"page": 3, "region": "Scheme 5 upper arrow"})
for args in [("f18", "cmp18", "reactant"), ("f7", "cmp7", "product"), ("fpdc", "reag_pdc", "catalyst"), ("fh2", "reag_h2", "reagent")]:
    kw = {"mol_percent": Q(10, "mol%"), "ref": "cmp18"} if args[1] == "reag_pdc" else {}
    add_participant(r, *args, (3,), **kw)
add_step(r, "hydrogenolysis", (3,), ["f18", "fpdc", "fh2", "f7"], "debenzylation with 10% Pd/C under H2", pressure="1 atm")
add_outcome(r, "yield7_quant", "yield", ["cmp7"], (3,), Q(100, "%"), basis="quantitative", raw="quant.")
r["review"]["unresolved"] = ["Temperature, time, solvent and workup are not reported on page 3."]
classify_partial(r, ["O-H", "C-H"], (3,))

r = new_reaction("scheme5_li_naph_18_to_7_19", (3,), name="lithium-naphthalenide debenzylation", transformation="benzyl ether cleavage with silyl migration", location={"page": 3, "region": "Scheme 5 lower arrow"})
for args in [("g18", "cmp18", "reactant"), ("g7", "cmp7", "product"), ("g19", "cmp19", "product"), ("glinaph", "reag_li_naph", "reagent")]: add_participant(r, *args, (3,))
add_step(r, "reduction", (3,), ["g18", "glinaph", "g7", "g19"], "lithium-naphthalenide treatment; desired alcohol and silyl-migrated product formed")
add_outcome(r, "yield7_22", "yield", ["cmp7"], (3,), Q(22, "%"), raw="7 (22%)")
add_outcome(r, "yield19_60", "yield", ["cmp19"], (3,), Q(60, "%"), raw="19 (60%)")
r["review"]["unresolved"] = ["Yield measurement bases, temperature, time, solvent, workup and purification are unreported on page 3."]
classify_partial(r, ["O-H", "C-H"], (3,))

r = new_reaction("scheme6_9bbn_li_naph_18_to_7", (3,), name="9-BBN-protected lithium-naphthalenide debenzylation", transformation="traceless borinate protection and benzyl ether cleavage", location={"page": 3, "region": "Scheme 6, 18 via 20 to 7"})
for args in [("h18", "cmp18", "reactant"), ("h20", "cmp20", "reagent"), ("h7", "cmp7", "product"), ("h9bbn", "reag_9bbn", "reagent"), ("hlinaph", "reag_li_naph", "reagent")]:
    kw = {"equiv": Q(1, "equiv"), "ref": "cmp18"} if args[1] == "reag_9bbn" else {}
    add_participant(r, *args, (3,), **kw)
add_step(r, "borinate_formation", (3,), ["h18", "h9bbn", "h20"], "treat 18 with 9-BBN to form borinate 20")
add_step(r, "debenzylation", (3,), ["h20", "hlinaph"], "lithium naphthalenide treatment")
add_step(r, "acid_workup", (3,), ["h7"], "then HCl workup", stage="workup")
add_outcome(r, "yield7_91", "yield", ["cmp7"], (3,), Q(91, "%"), basis="isolated", raw="91%")
r["review"]["unresolved"] = ["Temperature, time, solvent, lithium-naphthalenide charge and workup details are unreported on page 3."]
classify_partial(r, ["O-H", "C-H"], (3,))

r = new_reaction("scheme6_7_to_21", (3,), name="oxidative Wittig reaction", transformation="primary alcohol oxidation and olefination", location={"page": 3, "region": "Scheme 6, 7 to 21"})
for args in [("i7", "cmp7", "reactant"), ("i6", "cmp6", "reactant"), ("i21", "cmp21", "product"), ("itempo", "reag_tempo", "reagent"), ("ibaib", "reag_baib", "reagent"), ("ikhmds", "reag_khmds", "reagent")]: add_participant(r, *args, (3,))
add_step(r, "oxidation", (3,), ["i7", "itempo", "ibaib"], "selective TEMPO/BAIB oxidation")
add_step(r, "wittig", (3,), ["i6", "ikhmds", "i21"], "then phosphonium salt 6 and KHMDS")
add_outcome(r, "yield21", "yield", ["cmp21"], (3,), Q(70, "%"), basis="isolated", raw="70%")
r["review"]["unresolved"] = ["Exact charges, temperature, time, solvent, workup and purification are unreported on page 3."]
classify_partial(r, ["DB-C-C"], (3,))

# Figure 4 proposed stereochemical models.
r = new_reaction("figure4_proposed_hydroboration_models", (4,), record_type="proposed_mechanism", name="Houk, modified Houk and chelation models", transformation="proposed facial-selectivity models for hydroboration", location={"page": 4, "region": "Figure 4 A–C"})
for args in [("m23a", "cmp23a_S", "product"), ("m23b", "cmp23b_S", "product"), ("m9bbn", "reag_9bbn", "reagent"), ("mipcbh", "reag_ipcbh2_plus", "reagent")]: add_participant(r, *args, (4,))
add_note(r, "Models A and B predict anti,syn-23b; model C is the authors' proposed chelation model for syn,syn-23a. These are proposed pathways, not separate observed experiments.", (4,))
classify_partial(r, ["C-H", "C-O"], (4,), "Proposed mechanism only; product bond classes are not experimental atom mapping.")


# Table 1 optimization rows.
table1 = [
    # first reagent, eq, temp, time, second reagent, eq, temp, time, conversion, Z/E, ds, yield
    (1, "reag_thexbh2", 1.1, QR(-78, 0, "°C", "–78 to 0 °C"), Q(1, "h"), "reag_thexbh2", 3.5, Q(0, "°C"), Q(2, "h"), ">95", "1/6", "1.4/1", None),
    (2, "reag_sia2bh", 2.0, QR(-78, -30, "°C", "–78 to –30 °C"), Q(1, "h"), "reag_thexbh2", 3.5, QR(-30, 0, "°C", "–30 to 0 °C"), Q(4, "h"), None, "1/2", None, None),
    (3, "reag_sia2bh", 1.0, QR(-78, -30, "°C", "–78 to –30 °C"), Q(1, "h"), "reag_thexbh2", 3.5, QR(-30, 0, "°C", "–30 to 0 °C"), Q(4, "h"), None, "1/7", None, None),
    (4, "reag_9bbn", 2.1, QR(-30, -10, "°C", "–30 to –10 °C"), Q(1, "h"), "reag_thexbh2", 3.5, QR(-10, 0, "°C", "–10 to 0 °C"), Q(4, "h"), ">95", None, "2.5/1", None),
    (5, "reag_9bbn", 1.1, QR(-30, -10, "°C", "–30 to –10 °C"), Q(1, "h"), "reag_thexbh2", 3.5, QR(-10, 0, "°C", "–10 to 0 °C"), Q(4, "h"), ">95", None, "4.0/1", "conflict"),
    (6, "reag_9bbn", 1.1, QR(-30, 0, "°C", "–30 to 0 °C"), Q(1, "h"), "reag_ipcbh2_plus", 3.5, Q(-10, "°C"), Q(2, "h"), ">95", None, "5.0/1", 70),
    (7, "reag_9bbn", 1.1, QR(-30, -10, "°C", "–30 to –10 °C"), Q(1, "h"), "reag_ipcbh2_plus", 3.5, Q(-20, "°C"), Q(2, "h"), ">95", None, "6.4/1", 76),
    (8, "reag_9bbn", 1.1, QR(-30, -10, "°C", "–30 to –10 °C"), Q(1, "h"), "reag_ipcbh2_minus", 3.5, Q(0, "°C"), Q(5, "h"), "76", None, "6.2/1", None),
    (9, "reag_ipcbh2_plus", 4.4, QR(-78, -30, "°C", "–78 to –30 °C"), Q(3, "h"), None, None, None, None, ">95", None, "4.1/1", None),
]
for entry, first, feq, ftemp, ftime, second, seq, stemp, stime, conv, ze, ds, yld in table1:
    rac = entry <= 5
    substrate = "cmp22_rac_Z" if rac else "cmp22_S_Z"
    eisomer = "cmp22_rac_E" if rac else "cmp22_S_E"
    products = ["cmp23a_rac", "cmp23b_rac", "cmp23c_rac"] if rac else ["cmp23a_S", "cmp23b_S", "cmp23c_S"]
    r = new_reaction(f"table1_entry{entry}", (4,), categories=("optimization",), name="hydroboration–oxidation optimization", transformation="stereoselective hydroboration–oxidation of (Z)-diol 22", location={"page": 4, "region": f"Table 1 entry {entry}"})
    add_participant(r, "t1_sub", substrate, "reactant", (4,), mass=QR(150, 200, "mg", "150–200 mg"))
    add_participant(r, "t1_e", eisomer, "product", (4,), product_role="alkene-isomerization product")
    add_participant(r, "t1_a", products[0], "product", (4,), product_role="syn,syn major triol")
    add_participant(r, "t1_b", products[1], "product", (4,), product_role="anti,syn triol")
    add_participant(r, "t1_c", products[2], "product", (4,), product_role="isomerization-derived diol")
    add_participant(r, "t1_first", first, "reagent", (4,), equiv=Q(feq, "equiv"), ref=substrate)
    if second:
        add_participant(r, "t1_second", second, "reagent", (4,), equiv=Q(seq, "equiv"), ref=substrate)
    add_participant(r, "t1_naoh", "reag_naoh", "reagent", (4,))
    add_participant(r, "t1_h2o2", "reag_h2o2", "reagent", (4,))
    add_participant(r, "t1_thf", "solv_thf", "solvent", (4,), volume=Q(1.0, "mL"))
    set_solvent(r, (4,), "THF", volume=Q(1.0, "mL"))
    add_step(r, "hydroxy_protection_or_first_borane", (4,), ["t1_sub", "t1_first"], "first borane treatment / hydroxy-protection stage", temperature=ftemp, time=ftime)
    if second:
        add_step(r, "hydroboration", (4,), ["t1_second"], "second borane hydroboration stage", temperature=stemp, time=stime)
    add_step(r, "oxidation", (4,), ["t1_naoh", "t1_h2o2", "t1_a", "t1_b", "t1_c"], "then NaOH/H2O2 oxidation", temperature="r.t.", time=Q(1.5, "h"))
    if conv is not None:
        val = Q(float(conv.strip(">")), "%", qualifier="gt", raw_text=conv + "%") if conv.startswith(">") else Q(float(conv), "%")
        add_outcome(r, f"conv{entry}", "conversion", [substrate], (4,), val, basis="1H NMR (400 MHz, CDCl3)", raw=conv + "%")
    if ze:
        add_outcome(r, f"ze{entry}", "selectivity", [substrate, eisomer], (4,), {"ratio": ze}, raw=f"Z/E {ze}", scope={"metric": "Z/E", "order": [substrate, eisomer]})
    if ds:
        add_outcome(r, f"ds{entry}", "selectivity", products[:2], (4,), {"ratio": ds}, raw=f"ds {ds}", scope={"metric": "diastereoselectivity", "order": products[:2]})
    if yld == "conflict":
        a1 = st.fact("63", Q(63, "%"), evs(4))
        a2 = st.fact("61% isolated yield", Q(61, "%"), evs(4))
        add_outcome(r, "yield23a5", "yield", [products[0]], (4,), None, yield_type="isolated", raw="Table 1: 63%; prose: 61%", alternatives=[a1, a2])
        add_note(r, "Source conflict retained: Table 1 reports 63%, while adjacent prose reports 61% isolated yield for entry 5.", (4,))
        r["review"]["unresolved"].append("Entry 5 isolated yield conflicts between Table 1 (63%) and prose (61%).")
    elif yld is not None:
        add_outcome(r, f"yield23a{entry}", "yield", [products[0]], (4,), Q(yld, "%"), yield_type="isolated", basis="isolated", raw=f"{yld}%")
    if entry == 1:
        add_outcome(r, "ratio_sum_to_23c", "selectivity", products, (4,), {"ratio": "2.4/1"}, raw="Σ(23a+23b)/23c = 2.4/1", scope={"metric": "product ratio", "order": [[products[0], products[1]], products[2]]})
        add_note(r, "Footnote says 23c <5% in all entries but also gives Σ(23a+23b)/23c = 2.4/1 for entry 1; the apparent internal inconsistency is retained.", (4,))
        r["review"]["unresolved"].append("Table 1 footnote internally conflicts for entry 1 regarding 23c (<5% versus 2.4/1 product ratio).")
    else:
        add_outcome(r, f"minor23c_{entry}", "selectivity", [products[2]], (4,), "trace", basis="1H NMR (400 MHz, CDCl3)", raw="23c detected with a trace amount")
    if entry == 9:
        add_note(r, "Entry 9 table cell gives –78 to –30 °C, 3 h; footnote g additionally states the reaction was performed at –30 °C for 3 h.", (4,))
        r["review"]["unresolved"].append("Entry 9 temperature wording is retained without forcing the range and –30 °C footnote into one unsupported program.")
    r["review"]["unresolved"] += ["Workup and purification are not reported on page 4.", "Z/E and ds measurement methods are not explicitly stated in Table 1."]
    classify_partial(r, ["C-H", "C-O"], (4,))

# Scheme 7: preparation steps and a distinct gram-scale hydroboration.
r = new_reaction("scheme7_12_to_22", (5,), name="hydroboration–oxidation/deacetylation", transformation="allenic acetate to (Z)-1,4-diol", location={"page": 5, "region": "Scheme 7, 12 to 22"})
for args in [("j12", "cmp12_S", "reactant"), ("j22", "cmp22_S_Z", "product"), ("jsia", "reag_sia2bh", "reagent"), ("jtbaoh", "reag_tbaf", "reagent"), ("jh2o2", "reag_h2o2", "reagent"), ("jthf", "solv_thf", "solvent")]: add_participant(r, *args, (5,))
set_solvent(r, (5,), "THF")
add_step(r, "hydroboration", (5,), ["j12", "jsia", "jthf"], "Sia2BH in THF")
add_step(r, "oxidation_deacetylation", (5,), ["jtbaoh", "jh2o2", "j22"], "then TBAOH/H2O2")
add_outcome(r, "yield22", "yield", ["cmp22_S_Z"], (5,), Q(86, "%"), raw="86%")
r["review"]["unresolved"] = ["The scheme text appears as TBAOH; exact charge, temperature, time, workup, purification and yield basis are unreported."]
classify_partial(r, ["C-H", "C-O", "O-H"], (5,))

r = new_reaction("scheme7_scaleup_22_to_23a", (5,), categories=("scale_up",), name="gram-scale hydroboration–oxidation", transformation="stereoselective hydroboration–oxidation of (S)-(Z)-22", location={"page": 5, "region": "Scheme 7 and prose"})
for args in [("k22", "cmp22_S_Z", "reactant"), ("k23a", "cmp23a_S", "product"), ("k23b", "cmp23b_S", "product"), ("k9bbn", "reag_9bbn", "reagent"), ("kipc", "reag_ipcbh2_plus", "reagent"), ("knaoh", "reag_naoh", "reagent"), ("kh2o2", "reag_h2o2", "reagent"), ("kthf", "solv_thf", "solvent")]:
    eq = Q(1, "equiv") if args[1] == "reag_9bbn" else Q(3.3, "equiv") if args[1] == "reag_ipcbh2_plus" else None
    add_participant(r, *args, (5,), equiv=eq, ref="cmp22_S_Z" if eq else None)
set_solvent(r, (5,), "THF")
add_step(r, "borinate", (5,), ["k22", "k9bbn"], "22 first reacted with 9-BBN (1 equiv)")
add_step(r, "hydroboration", (5,), ["kipc", "kthf"], "hydroboration with (+)-IpcBH2 (3.3 equiv in THF)")
add_step(r, "oxidation", (5,), ["knaoh", "kh2o2", "k23a", "k23b"], "NaOH/H2O2 oxidation")
add_outcome(r, "yield23a_scale", "yield", ["cmp23a_S"], (5,), Q(73, "%"), raw="73%")
add_outcome(r, "ds23a_scale", "selectivity", ["cmp23a_S", "cmp23b_S"], (5,), {"ratio": "6.4/1"}, raw="ds 6.4/1", scope={"metric": "ds", "order": ["cmp23a_S", "cmp23b_S"]})
r["review"]["unresolved"] = ["Temperature/time/workup/purification and yield/selectivity measurement bases are unreported on page 5."]
classify_partial(r, ["C-H", "C-O"], (5,))

r = new_reaction("scheme7_23a_to_24", (5,), name="vicinal-diol protection and BPS deprotection", transformation="acetonide formation followed by desilylation", location={"page": 5, "region": "Scheme 7"})
for args in [("l23", "cmp23a_S", "reactant"), ("l24", "cmp24", "product"), ("ldmp", "reag_dmp_acetonide", "reagent"), ("lcsa", "reag_csa", "catalyst"), ("ltbaf", "reag_tbaf", "reagent"), ("ldcm", "solv_dcm", "solvent")]:
    eq = Q(4.0, "equiv") if args[1] == "reag_dmp_acetonide" else None
    mp = Q(0.5, "mol%") if args[1] == "reag_csa" else None
    add_participant(r, *args, (5,), equiv=eq, mol_percent=mp, ref="cmp23a_S" if eq or mp else None)
set_solvent(r, (5,), "DCM")
add_step(r, "protection", (5,), ["l23", "ldmp", "lcsa", "ldcm"], "2,2-dimethoxypropane (4.0 equiv), CSA (0.5 mol%) in DCM", temperature="r.t.", time=Q(30, "min"))
add_step(r, "desilylation", (5,), ["ltbaf", "l24"], "then TBAF", time=Q(2, "h"))
a_scheme = st.fact("Scheme 7: 80%", Q(80, "%"), evs(5))
a_prose = st.fact("prose: 82%", Q(82, "%"), evs(5))
add_outcome(r, "yield24", "yield", ["cmp24"], (5,), None, raw="Scheme 7: 80%; prose: 82%", alternatives=[a_scheme, a_prose])
add_note(r, "Source conflict retained: Scheme 7 gives 80%, while adjacent prose gives 82% for 23a→24.", (5,))
r["review"]["unresolved"] = ["Yield conflicts between Scheme 7 (80%) and prose (82%); yield basis, workup and purification are unreported."]
classify_partial(r, ["C-O"], (5,))

r = new_reaction("scheme7_24_to_25", (5,), name="oxidative Wittig reaction", transformation="diol oxidation/olefination to Michael acceptor 25", location={"page": 5, "region": "Scheme 7"})
for args in [("n24", "cmp24", "reactant"), ("n25", "cmp25", "product"), ("ntempo", "reag_tempo", "reagent"), ("nbaib", "reag_baib", "reagent"), ("nylide", "reag_ylide_ethyl", "reagent"), ("nthf", "solv_thf", "solvent")]: add_participant(r, *args, (5,))
set_solvent(r, (5,), "THF")
add_step(r, "oxidative_wittig", (5,), ["n24", "ntempo", "nbaib", "nylide", "nthf", "n25"], "TEMPO/BAIB oxidation followed by Ph3P=CHCO2Et olefination", temperature="r.t.")
add_outcome(r, "yield25", "yield", ["cmp25"], (5,), Q(81, "%"), raw="81%")
add_outcome(r, "ez25", "selectivity", ["cmp25"], (5,), {"ratio": ">20/1"}, raw="E/Z >20/1", scope={"metric": "E/Z"})
r["review"]["unresolved"] = ["Time, reagent charges, workup, purification and measurement bases are unreported on page 5."]
classify_partial(r, ["DB-C-C"], (5,))

# Table 2 cyclization rows; entry 6 also cites Scheme 7 and prose.
table2 = [
    (1, "reag_ptsa", 1.2, Q(50, "°C"), Q(12, "h"), "NR", None),
    (2, "reag_pcsa", 1.2, Q(50, "°C"), Q(10, "h"), "lt10", None),
    (3, "reag_pcsa", 2.0, Q(80, "°C"), Q(60, "h"), 69, "8/1"),
    (4, "reag_pnsa", 2.0, Q(80, "°C"), Q(85, "h"), 67, "7/1"),
    (5, "reag_tcsa", 2.0, Q(80, "°C"), Q(85, "h"), 59, "8/1"),
    (6, "reag_pcsa", 2.0, "reflux", Q(70, "h"), 82, "11/1"),
    (7, "reag_pcsa", 4.0, "reflux", Q(36, "h"), 59, "11/1"),
]
for entry, promoter, eq, temp, tm, yld, ds in table2:
    pages = (5,)
    r = new_reaction(f"table2_entry{entry}", pages, categories=("optimization",), name="oxa-1,4-addition optimization", transformation="acid-promoted cyclization to cis-tetrahydropyran 26", location={"page": 5, "region": f"Table 2 entry {entry}"})
    add_participant(r, "t2_25", "cmp25", "reactant", pages)
    if yld != "NR": add_participant(r, "t2_26", "cmp26", "product", pages)
    add_participant(r, "t2_prom", promoter, "reagent", pages, equiv=Q(eq, "equiv"), ref="cmp25")
    add_participant(r, "t2_etoh", "solv_ethanol", "solvent", pages)
    set_solvent(r, pages, "EtOH")
    pids = ["t2_25", "t2_prom", "t2_etoh"] + ([] if yld == "NR" else ["t2_26"])
    add_step(r, "cyclization", pages, pids, "acid-promoted cyclization", temperature=temp, time=tm)
    if yld == "NR":
        add_outcome(r, "no_reaction", "no_reaction", [], pages, True, raw="N.R.")
        classify_na(r, pages, "No reaction was reported; no expected-product bond label assigned.")
    else:
        value = Q(10, "%", qualifier="lt", raw_text="<10%") if yld == "lt10" else Q(yld, "%")
        add_outcome(r, f"yield26_{entry}", "yield", ["cmp26"], pages, value, yield_type="isolated", basis="isolated", raw="<10%" if yld == "lt10" else f"{yld}%")
        if ds:
            add_outcome(r, f"ds26_{entry}", "selectivity", ["cmp26"], pages, {"ratio": ds}, raw=f"ds {ds}", scope={"metric": "diastereoselectivity"})
        classify_partial(r, ["C-O"], pages)
    if entry == 6:
        add_note(r, "Same optimized 25→26 experiment is shown in Scheme 7 and discussed in prose; evidence is merged rather than duplicated.", pages)
        r["evidence_ids"] = evs(5)
    add_note(r, "The reaction sketch above Table 2 says EtOH, 90 °C generically; row-specific temperatures are retained and the generic 90 °C is not inherited into each row.", pages)
    r["review"]["unresolved"] += ["Workup/purification and ds measurement method are unreported on page 5."]

# Table 3 oxidation/HWE rows.
table3 = [
    (1, [("reag_tempo", 0.1), ("reag_baib", 1.2)], Q(4, "h"), ">95", "1.7/1", None, None),
    (2, [("reag_dmp_oxidant", 1.2)], Q(5, "h"), ">95", "complex", None, None),
    (3, [("reag_tempobf4", 2.6)], Q(29, "h"), ">95", "3/1", None, None),
    (4, [("reag_tempo", 0.1), ("reag_tcca", 0.75)], Q(39, "min"), ">95", "4.4/1", 35, None),
    (5, [("reag_27a", 0.75), ("reag_baib", 1.3)], Q(26, "h"), "<7", "8/1", None, "THF"),
    (6, [("reag_27b", 0.75), ("reag_baib", 1.3)], Q(17, "h"), "15", ">10/1", None, "THF"),
    (7, [("reag_27c", 0.75), ("reag_baib", 1.3)], Q(26, "h"), "NR", None, None, "THF"),
    (8, [("reag_27a", 0.25), ("reag_baib", 1.3)], Q(8, "h"), "80", "5/1", None, "CH2Cl2"),
    (9, [("reag_27b", 0.25), ("reag_baib", 1.3)], Q(8, "h"), "42", "10/1", None, "CH2Cl2"),
    (10, [("reag_27c", 0.25), ("reag_baib", 1.3)], Q(8, "h"), ">90", "10/1", None, "CH2Cl2"),
]
for entry, reagents, tm, conv, ratio, yld, solvent in table3:
    r = new_reaction(f"table3_entry{entry}", (5,), categories=("optimization",), name="oxidation/HWE optimization", transformation="oxidation of 26 followed by HWE olefination", location={"page": 5, "region": f"Table 3 entry {entry}"})
    add_participant(r, "t3_26", "cmp26", "reactant", (5,))
    successful = conv != "NR" and entry != 2
    if successful:
        add_participant(r, "t3_29s", "cmp29_11S_Z", "product", (5,), product_role="11S product")
        add_participant(r, "t3_29r", "cmp29_11R_Z", "product", (5,), product_role="11R epimer")
    for i, (mid, eq) in enumerate(reagents, 1):
        add_participant(r, f"t3_reag{i}", mid, "reagent", (5,), equiv=Q(eq, "equiv"), ref="cmp26")
    add_participant(r, "t3_hwe", "reag_hwe_p5", "reagent", (5,))
    if solvent:
        smid = "solv_thf" if solvent == "THF" else "solv_dcm"
        add_participant(r, "t3_solvent", smid, "solvent", (5,))
        set_solvent(r, (5,), solvent)
    pids = ["t3_26"] + [f"t3_reag{i}" for i in range(1, len(reagents)+1)] + (["t3_solvent"] if solvent else [])
    temp = "r.t." if solvent else None
    add_step(r, "oxidation", (5,), pids, "oxidation conditions listed in Table 3", temperature=temp, time=tm)
    add_step(r, "hwe", (5,), ["t3_hwe"] + (["t3_29s", "t3_29r"] if successful else []), "subsequent HWE reaction; reagent identity and conditions not shown in Table 3")
    if conv == "NR":
        add_outcome(r, "no_reaction", "no_reaction", [], (5,), True, raw="N.R.")
        classify_na(r, (5,), "No reaction reported in entry 7.")
    else:
        if conv.startswith(">"):
            cv = Q(float(conv[1:]), "%", qualifier="gt", raw_text=conv + "%")
        elif conv.startswith("<"):
            cv = Q(float(conv[1:]), "%", qualifier="lt", raw_text=conv + "%")
        else:
            cv = Q(float(conv), "%")
        add_outcome(r, f"conv_t3_{entry}", "conversion", ["cmp26"], (5,), cv, raw=conv + "%")
        if entry == 2:
            add_outcome(r, "complex_mixture", "other", [], (5,), "complex mixture", raw="complex")
            classify_unresolved(r, (5,), "DMP gave a complex mixture and no desired product identity; no bond formation class asserted.")
        else:
            add_outcome(r, f"ratio_t3_{entry}", "selectivity", ["cmp29_11S_Z", "cmp29_11R_Z"], (5,), {"ratio": ratio}, raw=f"29 (11S/11R) {ratio}", scope={"metric": "11S/11R", "order": ["cmp29_11S_Z", "cmp29_11R_Z"]})
            classify_partial(r, ["DB-C-C"], (5,))
        if yld is not None:
            add_outcome(r, "combined_yield29", "yield", ["cmp29_11S_Z", "cmp29_11R_Z"], (5,), Q(yld, "%"), yield_type="combined isolated", basis="isolated (prose); combined yield (footnote)", raw="combined yield of 29 was 35%")
    r["review"]["unresolved"] += ["HWE reagent identity/conditions, workup and purification are not reported in Table 3.", "Conversion and 11S/11R measurement methods are not stated on page 5."]

# Page 6 Scheme 8 scale-up and protection.
r = new_reaction("scheme8_scaleup_26_to_29", (5, 6), categories=("scale_up",), name="optimized oxidation/Still–Gennari olefination", transformation="diol oxidation and Z-selective HWE olefination", location={"pages": [5, 6], "region": "Scheme 8 and preceding prose"})
for args in [("o26", "cmp26", "reactant"), ("o28", "cmp28", "reagent"), ("o29z", "cmp29_11S_Z", "product"), ("o29e", "cmp29_11S_E", "product"), ("o27c", "reag_27c", "reagent"), ("obaib", "reag_baib", "reagent"), ("ohwe", "reag_still_gennari", "reagent"), ("ocrown", "reag_18c6", "additive"), ("okhmds", "reag_khmds", "reagent"), ("onahco3", "reag_nahco3", "additive"), ("odcm", "solv_dcm", "solvent"), ("othf", "solv_thf", "solvent")]: add_participant(r, *args, (5, 6))
set_solvent(r, (6,), ["CH2Cl2", "THF"])
add_step(r, "oxidation", (6,), ["o26", "o27c", "obaib", "onahco3", "odcm", "o28"], "4-benzyloxy-TEMPO/BAIB oxidation with NaHCO3 in CH2Cl2", temperature="r.t.")
add_step(r, "hwe", (6,), ["o28", "ohwe", "ocrown", "okhmds", "othf", "o29z", "o29e"], "then Still–Gennari phosphonate, 18-crown-6 and KHMDS in THF")
add_outcome(r, "yield29_scale", "yield", ["cmp29_11S_Z", "cmp29_11S_E"], (6,), Q(74, "%"), basis="isolated", raw="74%")
add_outcome(r, "ze29_scale", "selectivity", ["cmp29_11S_Z", "cmp29_11S_E"], (6,), {"ratio": "10/1"}, raw="Z/E 10/1", scope={"metric": "Z/E", "order": ["cmp29_11S_Z", "cmp29_11S_E"]})
r["review"]["unresolved"] = ["Exact charges, time, workup, purification and Z/E measurement method are not reported on page 6.", "The 74% yield appears beside the Z/E ratio; whether it is mixture or major-isomer yield is not explicit, so both geometrical product identities are targeted."]
classify_partial(r, ["DB-C-C"], (5, 6))

r = new_reaction("scheme8_29_to_10", (6,), name="TBS protection", transformation="alcohol silylation", location={"page": 6, "region": "Scheme 8"})
for args in [("p29", "cmp29_11S_Z", "reactant"), ("p10", "cmp10", "product"), ("ptbs", "reag_tbsotf", "reagent"), ("plut", "reag_26lutidine", "reagent"), ("pdcm", "solv_dcm", "solvent")]: add_participant(r, *args, (6,))
set_solvent(r, (6,), "CH2Cl2")
add_step(r, "silylation", (6,), ["p29", "ptbs", "plut", "pdcm", "p10"], "TBSOTf/2,6-lutidine in CH2Cl2")
add_outcome(r, "yield10", "yield", ["cmp10"], (6,), Q(87, "%"), raw="87%")
r["review"]["unresolved"] = ["Temperature, time, charges, workup, purification and yield basis are unreported on page 6."]
classify_partial(r, ["O-Si"], (6,))

# Scheme 9 standalone selective reduction.
r = new_reaction("scheme9_10_to_30", (6,), name="chemoselective DIBAL-H reduction", transformation="differential reduction of two esters", location={"page": 6, "region": "Scheme 9"})
for args in [("q10", "cmp10", "reactant"), ("q30", "cmp30", "product"), ("qdibal", "reag_dibalh", "reagent"), ("qdcm", "solv_dcm", "solvent")]:
    eq = Q(2.9, "equiv") if args[1] == "reag_dibalh" else None
    add_participant(r, *args, (6,), equiv=eq, ref="cmp10" if eq else None)
set_solvent(r, (6,), "CH2Cl2")
add_step(r, "reduction", (6,), ["q10", "qdibal", "qdcm", "q30"], "DIBAL-H (2.9 equiv) in CH2Cl2", temperature=Q(-78, "°C"))
add_outcome(r, "yield30", "yield", ["cmp30"], (6,), Q(91, "%"), basis="isolated", raw="91%")
r["review"]["unresolved"] = ["Reaction time, workup and purification are not reported in Scheme 9."]
classify_partial(r, ["C-H", "O-H"], (6,))

# Scheme 10 two distinct telescoped HWE variants and two DMP oxidations.
for suffix, phosphonate, product, minor, yld, ratio in [
    ("a", "cmp9a", "cmp31a_EE", "cmp31a_other", 64, "14/1"),
    ("b", "cmp9b", "cmp31b_EE", "cmp31b_other", 59, "13/1"),
]:
    r = new_reaction(f"scheme10_10_to_31{suffix}", (6,), name="telescoped DIBAL-H reduction/HWE elongation", transformation="ester reduction followed by HWE diene formation", location={"page": 6, "region": f"Scheme 10, 9{suffix} branch"})
    for args in [("r10", "cmp10", "reactant"), ("r30", "cmp30", "reagent"), ("rphos", phosphonate, "reactant"), ("rmajor", product, "product"), ("rminor", minor, "product"), ("rdibal", "reag_dibalh", "reagent"), ("rnh4cl", "reag_nh4cl", "workup_agent"), ("rlihmds", "reag_lihmds", "reagent"), ("rdcm", "solv_dcm", "solvent"), ("rthf", "solv_thf", "solvent")]:
        eq = Q(2.9, "equiv") if args[1] == "reag_dibalh" else None
        stage = "workup" if args[1] == "reag_nh4cl" else "reaction"
        add_participant(r, *args, (6,), equiv=eq, ref="cmp10" if eq else None, stage=stage)
    set_solvent(r, (6,), ["CH2Cl2", "THF"])
    add_step(r, "reduction", (6,), ["r10", "rdibal", "rdcm", "r30"], "DIBAL-H (2.9 equiv) in CH2Cl2; crude 30 not isolated", temperature=Q(-78, "°C"), time=Q(20, "min"))
    add_step(r, "quench", (6,), ["rnh4cl"], "then NH4Cl", stage="workup")
    add_step(r, "hwe", (6,), ["r30", "rphos", "rlihmds", "rthf", "rmajor", "rminor"], "LiHMDS-promoted HWE reaction in THF")
    add_outcome(r, f"yield31{suffix}", "yield", [product, minor], (6,), Q(yld, "%"), raw=f"{yld}% (2 steps)", scope={"span": f"10 to 31{suffix}", "reported_steps": 2})
    add_outcome(r, f"eez31{suffix}", "selectivity", [product, minor], (6,), {"ratio": ratio}, raw=f"E/E selectivity {ratio}", scope={"metric": "E,E/other geometrical isomer", "order": [product, minor]})
    r["review"]["unresolved"] = ["The source labels the ratio E/E in the product caption/prose as excellent E,E stereoselectivity; minor-isomer geometry and measurement method are not assigned.", "HWE temperature/time/charges, workup, purification and yield basis are unreported."]
    classify_partial(r, ["C-H", "O-H", "DB-C-C"], (6,))

for suffix, substrate, product, yld in [("a", "cmp31a_EE", "cmp32a", 94), ("b", "cmp31b_EE", "cmp32b", 93)]:
    r = new_reaction(f"scheme10_31{suffix}_to_32{suffix}", (6,), name="Dess–Martin oxidation", transformation="allylic alcohol oxidation to enal", location={"page": 6, "region": f"Scheme 10, 31{suffix} to 32{suffix}"})
    for args in [("s31", substrate, "reactant"), ("s32", product, "product"), ("sdmp", "reag_dmp_oxidant", "reagent")]: add_participant(r, *args, (6,))
    add_step(r, "oxidation", (6,), ["s31", "sdmp", "s32"], "Dess–Martin periodinane oxidation")
    add_outcome(r, f"yield32{suffix}", "yield", [product], (6,), Q(yld, "%"), raw=f"{yld}%")
    add_note(r, "Adjacent prose states that the Z-C12 alkene configuration and C11 stereochemistry were retained.", (6,))
    r["review"]["unresolved"] = ["Solvent, temperature, time, charges, workup, purification and yield basis are unreported on page 6."]
    classify_partial(r, ["DB-C-O"], (6,))

# Scheme 11 fragment transformations visible on page 6.
r = new_reaction("scheme11_21_to_33", (6,), name="modified Kita acetal conversion", transformation="acetal to aldehyde with concomitant TMS protection", location={"page": 6, "region": "Scheme 11, 21 to 33"})
for args in [("u21", "cmp21", "reactant"), ("u33", "cmp33", "product"), ("utms", "reag_tmsotf", "reagent"), ("ulut", "reag_26lutidine", "reagent"), ("udcm", "solv_dcm", "solvent")]: add_participant(r, *args, (6,))
set_solvent(r, (6,), "CH2Cl2")
add_step(r, "acetal_conversion", (6,), ["u21", "utms", "ulut", "udcm", "u33"], "TMSOTf/2,6-lutidine in CH2Cl2; aldehyde 33 not isolated")
add_outcome(r, "not_isolated33", "other", ["cmp33"], (6,), "not isolated", raw="not isolated")
r["review"]["unresolved"] = ["No standalone yield, temperature, time, charges, workup or purification are reported for 21→33."]
classify_partial(r, ["DB-C-O", "O-Si"], (6,))

r = new_reaction("scheme11_33_34a_to_35", (6,), name="Julia–Kocienski olefination", transformation="aldehyde–sulfone olefination", location={"page": 6, "region": "Scheme 11, 33 and 34a to 35"})
for args in [("v33", "cmp33", "reactant"), ("v34", "cmp34a", "reactant"), ("v35e", "cmp35_E", "product"), ("v35z", "cmp35_Z", "product"), ("vkhmds", "reag_khmds", "reagent"), ("vthf", "solv_thf", "solvent")]: add_participant(r, *args, (6,))
set_solvent(r, (6,), ["THF", "HMPA"])
add_step(r, "julia", (6,), ["v33", "v34", "vkhmds", "vthf", "v35e", "v35z"], "KHMDS-promoted Julia–Kocienski olefination with Kang's linker 34a in THF/HMPA")
add_outcome(r, "yield35", "yield", ["cmp35_E"], (6,), Q(82, "%"), raw="82%")
add_outcome(r, "ez35", "selectivity", ["cmp35_E", "cmp35_Z"], (6,), {"ratio": ">20/1"}, raw="E/Z >20/1", scope={"metric": "E/Z", "order": ["cmp35_E", "cmp35_Z"]})
r["review"]["unresolved"] = ["Exact charges, temperature, time, workup, purification and measurement bases are unreported on page 6."]
classify_partial(r, ["DB-C-C"], (6,))

r = new_reaction("scheme11_35_to_36a_36b", (6,), name="thioether oxidation", transformation="BT-thioether to BT-sulfones", location={"page": 6, "region": "Scheme 11, 35 to 36a/36b"})
for args in [("w35", "cmp35_E", "reactant"), ("w36a", "cmp36a", "product"), ("w36b", "cmp36b", "product"), ("wmoly", "reag_ammonium_molybdate", "reagent"), ("wh2o2", "reag_h2o2", "reagent")]: add_participant(r, *args, (6,))
add_step(r, "oxidation", (6,), ["w35", "wmoly", "wh2o2", "w36a", "w36b"], "oxidation with ammonium heptamolybdate tetrahydrate/H2O2")
add_outcome(r, "yield36a", "yield", ["cmp36a"], (6,), Q(54, "%"), raw="36a, 54%")
add_outcome(r, "yield36b", "yield", ["cmp36b"], (6,), Q(35, "%"), raw="36b, 35%")
r["review"]["unresolved"] = ["The page does not clarify whether 36a and 36b were isolated from one oxidation or distinct workup variants; both product-targeted yields are retained in one arrow-level record.", "Charges, solvent, temperature, time, workup, purification and yield bases are unreported."]
classify_partial(r, ["DB-O-S"], (6,))


# Source context for non-experiment regions and characterization triage.
BUNDLE["source_context"] = [
    {"context_id": f"{DOC}:page1_figure1", "evidence_ids": evs(1), "kind": "structural_background",
     "description": "Figure 1 compares originally proposed and revised lasonolide A structures; it is not a current synthesis experiment."},
    {"context_id": f"{DOC}:page2_scheme2", "evidence_ids": evs(2), "kind": "retrosynthetic_analysis",
     "description": "Scheme 2 is retrosynthetic analysis with disconnections/building blocks 6–10, not a set of executed experiments."},
    {"context_id": f"{DOC}:page2_figure2", "evidence_ids": evs(2), "kind": "structural_comparison",
     "description": "Figure 2 illustrates hidden symmetry between fragments 7 and 10."},
    {"context_id": f"{DOC}:page3_figure3", "evidence_ids": evs(3), "kind": "method_background",
     "description": "Figure 3 contrasts the authors' target syn,syn triol with the cited Evans approach; it is not a current specific experiment."},
    {"context_id": f"{DOC}:page6_known_linkers", "evidence_ids": evs(6), "kind": "literature_context",
     "description": "The Lee (2002), Kang (2003), and Ghosh (2007) linker drawings in Scheme 11 are cited prior-work structures, not current-paper experimental arrows."},
    {"context_id": f"{DOC}:pages1_6_characterization_triage", "evidence_ids": evs(1, 2, 3, 4, 5, 6), "kind": "characterization_triage",
     "description": "No standalone mass-spectrum peak table, molecular-weight table, or fragment-relative-abundance table appears on pages 1–6."},
]

BUNDLE["molecules"] = list(MOLECULES.values())
BUNDLE["reactions"] = REACTIONS
st.write(RUN / "submitted.json", BUNDLE)
st.write(RUN / "raw_model/stereo_checks.json", STEREO_DIAGNOSTICS)

# Raw page inventory and exact coverage artifact.
by_page = {page: [] for page in range(1, 7)}
for r in REACTIONS:
    for eid in r["evidence_ids"]:
        marker = ":page"
        if marker in eid:
            page = int(eid.split(marker, 1)[1].split(":", 1)[0])
            if 1 <= page <= 6 and r["reaction_id"] not in by_page[page]:
                by_page[page].append(r["reaction_id"])

regions = {
    1: [
        {"region": "title/summary and Background text", "finding": "Paper metadata and overview; no specific current experiment.", "evidence_ids": evs(1)},
        {"region": "Scheme 1", "finding": "Authors' previous generic hydroboration–allylation/oxa-Michael approach; recorded as a generic scheme.", "evidence_ids": evs(1)},
        {"region": "Figure 1", "finding": "Structural-revision background, not a synthesis experiment.", "evidence_ids": evs(1)},
    ],
    2: [
        {"region": "Scheme 2", "finding": "Retrosynthetic analysis only; retained in source_context.", "evidence_ids": evs(2)},
        {"region": "Figure 2 and Scheme 3", "finding": "Fragment symmetry and conceptual enantio-integrative routing; Scheme 3 recorded as generic, not specific experiments.", "evidence_ids": evs(2)},
        {"region": "A-ring synthesis prose", "finding": "Identifies the current 11+13→14, 14→15, and 15→16 operations continued in Scheme 4 on page 3.", "evidence_ids": evs(2)},
    ],
    3: [
        {"region": "Scheme 4", "finding": "Five current A-ring transformations 11+13→14→15→16→17→18, including yields and selectivity.", "evidence_ids": evs(3)},
        {"region": "Scheme 5", "finding": "Two distinct debenzylation experiments: Pd/C hydrogenation and lithium-naphthalenide reaction with silyl-migrated byproduct 19.", "evidence_ids": evs(3)},
        {"region": "Scheme 6", "finding": "9-BBN-protected debenzylation and oxidative Wittig preparation of 21.", "evidence_ids": evs(3)},
        {"region": "Figure 3", "finding": "Method/background comparison, not a current specific experiment.", "evidence_ids": evs(3)},
    ],
    4: [
        {"region": "Figure 4 A–C", "finding": "Three author-discussed stereochemical models retained in one proposed-mechanism record.", "evidence_ids": evs(4)},
        {"region": "Table 1 entries 1–9 and footnotes", "finding": "Nine independent optimization experiments inventoried row-by-row, including conversion, Z/E, ds, yields and 23c statements.", "evidence_ids": evs(4)},
        {"region": "Table 1 surrounding prose", "finding": "Entry-5 61% prose versus 63% table conflict and mechanistic discussion retained.", "evidence_ids": evs(4)},
    ],
    5: [
        {"region": "Scheme 7 and prose", "finding": "Four preparative records plus the 25→26 optimized entry shared with Table 2; scale-up 22→23a kept distinct from Table 1 entry 7.", "evidence_ids": evs(5)},
        {"region": "Table 2 entries 1–7", "finding": "Seven cyclization optimization rows, including a true no-reaction control, inventoried separately.", "evidence_ids": evs(5)},
        {"region": "Table 3 entries 1–10", "finding": "Ten oxidation/HWE optimization rows, including a complex-mixture result and a no-reaction control, inventoried separately.", "evidence_ids": evs(5)},
    ],
    6: [
        {"region": "Schemes 8–10", "finding": "Scale-up 26→29, silylation, standalone reduction, two telescoped HWE variants, and two DMP oxidations recorded.", "evidence_ids": evs(6)},
        {"region": "Scheme 11 lower sequence", "finding": "Current 21→33, 33+34a→35, and 35→36a/36b arrows recorded.", "evidence_ids": evs(6)},
        {"region": "Scheme 11 known linkers", "finding": "Lee/Kang/Ghosh structures are cited prior-work context; not imported as current experiments.", "evidence_ids": evs(6)},
    ],
}
unresolved = {
    1: ["Generic R groups and prior-work conditions prevent concrete molecule/experiment reconstruction."],
    2: ["Several exact operations continue onto page 3; quantities/times are absent in the covered prose."],
    3: ["Complex exact structures, many reagent charges, times, workups and purification details are absent without supplied SI."],
    4: ["Table 1 entry 5 yield conflict; entry 1 23c footnote inconsistency; entry 9 temperature wording; exact complex structures remain unresolved."],
    5: ["23a→24 yield conflict; HWE identity/conditions in Table 3 and several measurement methods are not supplied."],
    6: ["36a/36b experiment relationship and many exact charges/workups are not explicit; complex exact structures remain unresolved."],
}
coverage = {"document_id": DOC, "semantic_accuracy_verified": False, "pages": []}
for page in range(1, 7):
    coverage["pages"].append({"page": page, "state": "attempted", "regions": regions[page],
                              "record_ids": by_page[page], "unresolved": unresolved[page]})
st.write(RUN / "coverage/pages_0001_0006.json", coverage)
st.write(RUN / "raw_model/page_inventory.json", {"document_id": DOC, "pages": regions, "record_ids_by_page": by_page})
st.write(RUN / "raw_model/record_ids_by_page.json", by_page)
print(json.dumps({"molecules": len(MOLECULES), "reactions": len(REACTIONS),
                  "records_by_page": {str(k): len(v) for k, v in by_page.items()}}, ensure_ascii=False, indent=2))
