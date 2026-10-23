#!/usr/bin/env python3
"""Build the evidence-linked v3.1 submission for Org. Lett. 2020, 22, 3979–3984.

This is a source transcription generator.  It intentionally leaves unavailable SI
details and uncertain structures unresolved; it does not certify chemistry.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills" / "chemical-reaction-extraction" / "scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st


DOC = "doc_1b75b17cc3d81272"
RUN = ROOT / "runs" / "current"
DATASET = RUN / "dataset.json"
TITLE = "Organoseleno-Catalyzed Synthesis of α,β-Unsaturated α′-Alkoxy Ketones from Allenes Enabled by Se···O Interactions"


def ev(page: int, *kinds: str) -> list[str]:
    kinds = kinds or ("image", "text")
    return [f"{DOC}:page{page}:{kind}" for kind in kinds]


def F(raw=None, normalized=None, evidence=(), method="original_extraction", derivation=None,
      status=None, verification=None, alternatives=()):
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation,
                   status=status, verification=verification, alternatives=alternatives)


def inferred(normalized, evidence, derivation):
    return F(None, normalized, evidence, method="inference", derivation=derivation)


def derived(normalized, evidence, derivation):
    return F(None, normalized, evidence, method="structure_derivation", derivation=derivation)


def na(evidence=()):
    return F(None, None, evidence, method="not_applicable", status="not_applicable", verification="missing")


dataset = st.read(DATASET)
assert [d["document_id"] for d in dataset["documents"]] == [DOC]
dataset["documents"][0]["title"] = TITLE
st.write(DATASET, dataset)
bundle = st.empty(dataset["documents"])

molecules: dict[str, dict] = {}
stereo_checks: dict[str, dict] = {}


def molecule(mid: str, label: str, evidence: list[str], *, name=None, smiles=None,
             chirality="unresolved", geometry="unresolved", representation=None,
             double_bonds=None, axial=None, mixture=None, assignment=None,
             identity_note=None):
    if mid in molecules:
        occ = molecules[mid]["occurrences"][0]
        occ["evidence_ids"] = list(dict.fromkeys(occ["evidence_ids"] + evidence))
        return molecules[mid]
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, evidence)
    m["names"]["original"] = F(label, name or label, evidence)
    m["names"]["english"] = F(name or label, name or label, evidence)
    m["names"]["abbreviations"] = F()
    m["names"]["aliases"] = F()
    m["cas_numbers"] = F()
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = derived(
            smiles, evidence,
            "Candidate SMILES transcribed from the source drawing; parsing checks syntax only and the structure remains unverified.",
        )
        diag = st.stereo_inspect(smiles)
        stereo_checks[mid] = diag
        m["stereochemistry"]["isomeric_smiles"] = derived(
            diag["canonical_isomeric_smiles"], evidence,
            "RDKit canonical isomeric SMILES calculated from the unverified candidate graph.",
        )
        m["stereochemistry"]["identity_key"] = derived(
            diag["stereo_graph_key"], evidence,
            "RDKit stereo graph key calculated from the unverified candidate graph; not a sample identity.",
        )
        representation = representation or "partial"
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = F()
        representation = representation or "not_encoded"
    if chirality in {"achiral", "unresolved", "not_applicable"}:
        m["stereochemistry"]["chirality_status"] = inferred(
            chirality, evidence,
            "Status assigned from the candidate graph/identity review; absence of source stereo marks alone was not treated as proof of achirality.",
        )
    else:
        m["stereochemistry"]["chirality_status"] = F(chirality, chirality, evidence)
    if geometry in {"not_applicable", "unresolved"}:
        m["stereochemistry"]["geometry_status"] = inferred(
            geometry, evidence,
            "E/Z applicability/status assessed from the displayed structure; allene axial stereochemistry is tracked separately.",
        )
    else:
        m["stereochemistry"]["geometry_status"] = F(geometry, geometry, evidence)
    m["stereochemistry"]["representation_status"] = inferred(
        representation, evidence,
        "Encoding coverage assessed by comparing the source drawing with the candidate or unresolved representation.",
    )
    m["stereochemistry"]["atom_centers"] = F()
    m["stereochemistry"]["double_bonds"] = (
        F(double_bonds, double_bonds, evidence) if double_bonds is not None else F()
    )
    m["stereochemistry"]["axial"] = F(axial, axial, evidence) if axial is not None else F()
    m["stereochemistry"]["mixture"] = F(mixture, mixture, evidence) if mixture is not None else F()
    m["stereochemistry"]["assignment"] = F(assignment, assignment, evidence) if assignment is not None else F()
    m["identity_note"] = F(identity_note, identity_note, evidence) if identity_note else F()
    m["external_properties"] = []
    molecules[mid] = m
    return m


# Common reagents and solvents.  CAS values are deliberately unreported.
molecule(f"{DOC}:diphenyl_diselenide", "(PhSe)2", ev(2), name="diphenyl diselenide",
         smiles="c1ccc([Se][Se]c2ccccc2)cc1", chirality="achiral", geometry="not_applicable")
molecule(f"{DOC}:pyfotf", "[pyF][OTf]", ev(2), name="1-fluoropyridinium triflate",
         chirality="achiral", geometry="not_applicable", representation="not_encoded",
         identity_note="The ionic structure is source-drawn; no concrete registry SMILES is asserted.")
molecule(f"{DOC}:mecn", "MeCN", ev(2), name="acetonitrile", smiles="CC#N",
         chirality="achiral", geometry="not_applicable")
molecule(f"{DOC}:thf", "THF", ev(2), name="tetrahydrofuran", smiles="C1CCOC1",
         chirality="achiral", geometry="not_applicable")
molecule(f"{DOC}:water", "H2O", ev(1, "image") + ev(3, "text"), name="water", smiles="O",
         chirality="achiral", geometry="not_applicable")
molecule(f"{DOC}:phsecl", "PhSeCl", ev(2), name="phenylselenyl chloride", smiles="Cl[Se]c1ccccc1",
         chirality="achiral", geometry="not_applicable")
molecule(f"{DOC}:h218o", "H2(18)O", ev(3), name="oxygen-18 water", smiles="[18OH2]",
         chirality="achiral", geometry="not_applicable")
molecule(f"{DOC}:oxygen", "O2", ev(3), name="oxygen", smiles="O=O",
         chirality="achiral", geometry="not_applicable")


# Candidate source-drawing transcriptions.  Products with E/Z ratios are registered
# as the depicted/author-assigned major E structure; the ratio stays experiment-local.
product_smiles = {
    "2a": "COCC(=O)/C(C)=C/c1ccc(OC)cc1",
    "2b": "COCC(=O)/C(C)=C/c1ccccc1",
    "2c": "COCC(=O)/C(C)=C/c1ccc(C)cc1",
    "2d": "COCC(=O)/C(C)=C/c1ccc(Cl)cc1",
    "2e": "COCC(=O)/C(C)=C/c1cccc(Cl)c1",
    "2f": "COCC(=O)/C(C)=C/c1ccccc1Cl",
    "2g": "COCC(=O)/C(C)=C/c1ccccc1Br",
    "2h": "COCC(=O)/C=C/c1ccc(OC)cc1",
    "2i": "COCC(=O)/C=C/c1cccc(Cl)c1",
    "2j": "COCC(=O)/C(c1ccccc1)=C/c1ccccc1",
    "2k": "COCC(=O)/C(c1ccccc1)=C/c1ccc(Cl)cc1",
    "2l": "COCC(=O)/C(C)=C/c1ccsc1",
    "2m": "COCC(=O)C=C(C)C#Cc1ccccc1",
    "2n": "COC(C)C(=O)/C=C/c1ccccc1",
    "2r": "O=C(COCc1ccccc1)/C(C)=C/c1ccccc1",
    "2s": "O=C(COCc1ccccc1)/C=C/c1ccc(C)cc1",
    "2t": "C=CCOCC(=O)/C(C)=C/c1ccccc1Br",
    "2u": "C#CCOCC(=O)/C(C)=C/c1ccccc1",
    "2v": "C#CCOCC(=O)/C=C/c1ccc(C)cc1",
    "2w": "C#CCOCC(=O)/C(C)=C/c1ccccc1Br",
    "2x": "C=C=CCOCC(=O)/C(C)=C/c1ccccc1",
}


def product_molecule(label: str, page: int):
    chirality = "unspecified" if label == "2n" else "achiral"
    geometry = "specified"
    rep = "not_encoded" if label == "2m" else ("partial" if label in {"2n", "2x"} else "encoded")
    return molecule(
        f"{DOC}:{label}", label, ev(page),
        name=f"source-drawn α,β-unsaturated α′-alkoxy ketone {label}",
        smiles=product_smiles[label], chirality=chirality, geometry=geometry, representation=rep,
        double_bonds=[{"assignment": "E major isomer", "basis": "Scheme drawing and E/Z order"}],
        assignment={"source_statement": "(E)-stereochemistry assigned by NOE on 2c and applied by the authors to products 2"},
        identity_note="Candidate connectivity and E encoding require expert verification against the source drawing; reported E/Z composition is retained per experiment.",
    )


for label in [f"2{x}" for x in "abcdefghijklmn"]:
    product_molecule(label, 2)
for label in [f"2{x}" for x in "rstuvwx"]:
    product_molecule(label, 3)


# Only the clearly readable starting-material structures needed to avoid false
# synthesis/characterization ambiguity receive candidate SMILES.
molecule(f"{DOC}:1a", "1a", ev(1) + ev(2), name="allenyl methyl ether 1a",
         smiles="COC(c1ccc(OC)cc1)C(C)=C=C", chirality="unspecified", geometry="not_applicable",
         representation="partial",
         identity_note="The benzylic tetrahedral stereocenter is not configured in the source drawing; the terminal CH2 allene is not assigned an axis.")
for letter in "bcdefghijklmno":
    is_1o = letter == "o"
    molecule(f"{DOC}:1{letter}", f"1{letter}", ev(2), name=f"source-drawn allene 1{letter}",
             chirality="specified" if is_1o else "unspecified", geometry="not_applicable", representation="not_encoded",
             axial=[{"axis": "C=C=C", "assignment": "unspecified or not safely encodable from the supplied drawing"}],
             assignment={"source": "wedge/dash stereochemistry is drawn in 1o; exact descriptors are not safely encoded"} if is_1o else None,
             identity_note="Dense scope drawing retained by label; no unsafe complete SMILES was asserted.")
molecule(f"{DOC}:1p", "1p", ev(3), name="allenol 1p", smiles="OC(c1ccccc1)C(C)=C=C",
         chirality="unspecified", geometry="not_applicable", representation="partial",
         identity_note="Candidate connectivity from Scheme 4; tetrahedral stereochemistry is unassigned.")
molecule(f"{DOC}:1q", "1q", ev(3), name="allenyl acetate 1q",
         smiles="CC(=O)OC(c1ccc(C)cc1)C(C)=C=C", chirality="unspecified", geometry="not_applicable",
         representation="partial", identity_note="Candidate connectivity from Scheme 4; tetrahedral stereochemistry is unassigned.")
for letter in "rstuvwx":
    molecule(f"{DOC}:1{letter}", f"1{letter}", ev(3), name=f"source-drawn allenyl ether 1{letter}",
             chirality="unspecified", geometry="not_applicable", representation="not_encoded",
             axial=[{"axis": "C=C=C", "assignment": "unspecified"}],
             identity_note="Scope substrate retained by source label; product connectivity is represented separately.")

molecule(f"{DOC}:3a", "3a", ev(2), name="Meyer–Schuster-type enone 3a",
         smiles="CC(=O)/C(C)=C/c1ccc(OC)cc1", chirality="achiral", geometry="unspecified",
         representation="partial", identity_note="Candidate connectivity from Scheme 2; alkene geometry was not separately reported.")
molecule(f"{DOC}:3b", "3b", ev(3), name="enone 3b",
         smiles="CC(=O)/C(C)=C/c1ccccc1", chirality="achiral", geometry="unspecified",
         representation="partial", identity_note="Candidate connectivity from Scheme 4; alkene geometry was not separately reported.")
molecule(f"{DOC}:selenated_mixture", "mixture of selenated products", ev(2),
         name="mixture of selenated products", chirality="unresolved", geometry="unresolved",
         representation="not_encoded", mixture={"composition": "unreported"})
molecule(f"{DOC}:complicated_mixture", "complicated mixture", ev(2),
         name="complicated reaction mixture", chirality="unresolved", geometry="unresolved",
         representation="not_encoded", mixture={"composition": "unreported"})
molecule(f"{DOC}:18O-2e", "18O-2e", ev(3), name="oxygen-18-labeled product 2e",
         smiles="COCC(=[18O])/C(C)=C/c1cccc(Cl)c1", chirality="achiral", geometry="specified",
         representation="encoded", double_bonds=[{"assignment": "E major isomer"}],
         identity_note="The carbonyl oxygen is shown as 18O in Scheme 5; candidate SMILES is a structure derivation.")


def participant(r, pid, mid, role, evidence, *, mol_percent=None, equivalents=None,
                stage="reaction", product_role=None):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    if product_role is not None:
        p["product_role"] = F(product_role, product_role, evidence)
    if mol_percent is not None:
        p["charges"]["mol_percent"] = F(
            f"{mol_percent} mol%", {"number": mol_percent, "unit": "mol%"}, evidence
        )
        p["equivalence_reference"] = inferred(
            r["participants"][0]["molecule_id"] if r["participants"] else mid,
            evidence, "Loading denominator inferred as the displayed allene substrate.",
        )
    if equivalents is not None:
        p["charges"]["equivalents"] = F(
            f"{equivalents} equiv", {"number": equivalents, "unit": "equiv"}, evidence
        )
        p["equivalence_reference"] = inferred(
            r["participants"][0]["molecule_id"] if r["participants"] else mid,
            evidence, "Equivalent denominator inferred as the displayed allene substrate.",
        )
    r["participants"].append(p)
    return p


def base_reaction(suffix, page, record_type="specific_experiment", categories=None, location=None):
    evidence = ev(page)
    r = st.new_reaction(f"{DOC}:{suffix}", DOC, evidence)
    r["record_type"] = record_type
    if categories is None:
        r["experiment_categories"] = na(evidence)
    else:
        r["experiment_categories"] = F(categories, categories, evidence)
    loc = location or f"PDF page {page}"
    r["location"] = F(loc, {"pdf_page": page, "region": loc}, evidence)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    return r


def set_solvents_and_step(r, page, time_hours=None, *, anhydrous=False, oxygen=False):
    evidence = ev(page)
    substrate_mid = r["participants"][0]["molecule_id"]
    participant(r, f"{r['reaction_id']}:mecn", f"{DOC}:mecn", "solvent", evidence)
    participant(r, f"{r['reaction_id']}:thf", f"{DOC}:thf", "solvent", evidence)
    r["solvent_conditions"]["composition"] = F(
        "anh. MeCN/anh. THF" if anhydrous else "MeCN/THF",
        ["anhydrous acetonitrile", "anhydrous THF"] if anhydrous else ["acetonitrile", "THF"], evidence,
    )
    r["solvent_conditions"]["mixture_ratio"] = F("1:1", "1:1", evidence)
    s = st.new_step(f"{r['reaction_id']}:reaction", evidence)
    s["participant_ids"] = [p["participant_id"] for p in r["participants"]]
    s["operation"] = F("reaction under the displayed arrow conditions", "reaction", evidence)
    s["temperature"] = F("RT", "room temperature", evidence)
    if time_hours is not None:
        s["time"] = F(f"{time_hours} h", {"number": time_hours, "unit": "h"}, evidence)
    if oxygen:
        s["atmosphere"] = F("O2 (1 bar)", "oxygen atmosphere", evidence)
        s["pressure"] = F("1 bar", "1 bar", evidence)
    r["steps"] = [s]
    r["reaction_scale"] = F()
    return substrate_mid


def outcome(r, oid, kind, target, raw, normalized, evidence, *, basis=None, yield_type=None, scope=None):
    o = st.new_outcome(oid, kind, [target] if target else [], evidence)
    o["value"] = F(raw, normalized, evidence)
    if yield_type is not None:
        o["yield_type"] = F(yield_type, yield_type, evidence)
    if basis is not None:
        o["measurement_basis"] = F(basis, basis, evidence)
    if scope is not None:
        o["scope"] = F(scope, scope, evidence)
    r["outcomes"].append(o)
    return o


def classify_oxygenative(r, evidence, limitation=None):
    changes = [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}],
         "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}],
         "before_order": 0, "after_order": 2, "change_type": "formed"},
        {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
         "before_order": 1, "after_order": 2, "change_type": "order_increased"},
    ]
    labels = [
        {"category_number": 3, "code": "C-O", "supporting_change_ids": ["b1"]},
        {"category_number": 14, "code": "DB-C-O", "supporting_change_ids": ["b2"]},
        {"category_number": 13, "code": "DB-C-C", "supporting_change_ids": ["b3"]},
    ]
    deriv = "Partial abstract bond-change transcription from the source-drawn allene-to-enone rearrangement; no atom mapping is asserted."
    r["classification"]["mapped_reaction_smiles"] = F()
    r["classification"]["bond_changes"] = [derived(c, evidence, deriv) for c in changes]
    r["classification"]["labels"] = derived(labels, evidence, deriv)
    r["classification"]["coverage"] = inferred("partial", evidence, deriv)
    r["classification"]["limitations"] = [
        limitation or "Exact atom correspondence, broken bonds, internal proton transfers, and coproducts remain unresolved without a verified atom map."
    ]
    r["classification"]["model_votes"] = []


def classify_none(r, evidence, reason):
    r["classification"]["mapped_reaction_smiles"] = F()
    r["classification"]["bond_changes"] = []
    r["classification"]["labels"] = inferred([], evidence, reason)
    r["classification"]["coverage"] = inferred("not_applicable", evidence, reason)
    r["classification"]["limitations"] = [reason]
    r["classification"]["model_votes"] = []


def classify_unresolved(r, evidence, reason):
    r["classification"]["mapped_reaction_smiles"] = F()
    r["classification"]["bond_changes"] = []
    r["classification"]["labels"] = inferred([], evidence, reason)
    r["classification"]["coverage"] = inferred("unresolved", evidence, reason)
    r["classification"]["limitations"] = [reason]
    r["classification"]["model_votes"] = []


reactions = []


def add_generic(suffix, page, record_type, location, name, note):
    r = base_reaction(suffix, page, record_type=record_type, categories=None, location=location)
    r["reaction_name"] = F(name, name, ev(page))
    r["transformation_type"] = F(note, note, ev(page))
    classify_unresolved(r, ev(page), "Generic or proposed drawing lacks concrete verified atom mapping and is not a specific experiment.")
    r["notes"] = [F(note, note, ev(page))]
    r["review"]["unresolved"] = ["Generic/literature/mechanistic drawing is contextual and must not be counted as a newly reported experiment."]
    reactions.append(r)
    return r


add_generic("scheme1a_prior_alkene", 2, "literature_example", "Scheme 1a", "prior alkene functionalization",
            "Cited previous selenium-catalyzed allylic or vicinal functionalization; ref. 3.")
add_generic("scheme1b_prior_alkyne", 2, "literature_example", "Scheme 1b", "prior alkyne functionalization",
            "Cited previous selenium-catalyzed access to ynones or oxazoles; ref. 6b.")
add_generic("scheme1c_prior_diene", 2, "literature_example", "Scheme 1c", "prior diene pyridination",
            "Cited previous selenium-catalyzed pyridination; ref. 6a.")
scheme1d = add_generic("scheme1d_current_overview", 2, "generic_scheme", "Scheme 1d", "current-study general transformation",
                       "Generic allene oxygenative rearrangement to α,β-unsaturated α′-alkoxy ketones.")
scheme3_generic = add_generic("scheme3_scope_generic", 2, "generic_scheme", "Scheme 3 general arrow",
                              "controlled synthesis of α,β-unsaturated α′-alkoxy ketones 2a–n",
                              "Scope-level general reaction drawing; individual substrate rows are separate records.")
scheme4_generic = add_generic("scheme4_scope_generic", 3, "generic_scheme", "Scheme 4 general arrow",
                              "controlled synthesis of α,β-unsaturated α′-alkoxy ketones 2r–x",
                              "Scope-level general reaction drawing; individual substrate rows are separate records.")


def scheme2_record(suffix, loading, time_h, yield_2a, *, yield_3a=None, extra_product=None, scope_2a=False):
    evidence = ev(2) + (ev(1, "text") if loading in (100, 10, 20, 50) else [])
    r = base_reaction(suffix, 2, categories=["optimization"] + (["substrate_scope"] if scope_2a else []),
                      location=f"Scheme 2, {suffix.split('_')[0]}")
    r["general_scheme_ids"] = [scheme1d["reaction_id"]] + ([scheme3_generic["reaction_id"]] if scope_2a else [])
    r["reaction_name"] = F("selenium-catalyzed rearrangement–functionalization", "selenium-catalyzed rearrangement–functionalization", evidence)
    r["transformation_type"] = F("oxygenative rearrangement", "allene-to-enone oxygenative rearrangement", evidence)
    participant(r, f"{r['reaction_id']}:substrate", f"{DOC}:1a", "reactant", evidence)
    participant(r, f"{r['reaction_id']}:product2a", f"{DOC}:2a", "product", evidence, product_role="major or reported product")
    if yield_3a is not None:
        participant(r, f"{r['reaction_id']}:product3a", f"{DOC}:3a", "product", evidence, product_role="major Meyer–Schuster-type product")
    if extra_product:
        participant(r, f"{r['reaction_id']}:sideproducts", extra_product, "product", evidence, product_role="unresolved side-product mixture")
    participant(r, f"{r['reaction_id']}:catalyst", f"{DOC}:diphenyl_diselenide", "catalyst", evidence, mol_percent=loading)
    participant(r, f"{r['reaction_id']}:oxidant", f"{DOC}:pyfotf", "reagent", evidence, equivalents=1.2)
    participant(r, f"{r['reaction_id']}:water", f"{DOC}:water", "reagent", evidence)
    set_solvents_and_step(r, 2, time_h)
    outcome(r, f"{r['reaction_id']}:yield2a", "yield", f"{DOC}:2a", f"{yield_2a}%",
            {"number": yield_2a, "unit": "%"}, evidence,
            basis="isolated yield" if scope_2a else None,
            yield_type="isolated" if scope_2a else None, scope="1a to 2a")
    if scope_2a:
        outcome(r, f"{r['reaction_id']}:ez", "selectivity", f"{DOC}:2a", ">99:1 E/Z", ">99:1",
                evidence, scope={"metric": "E/Z ratio", "order": ["E", "Z"]})
    if yield_3a is not None:
        outcome(r, f"{r['reaction_id']}:yield3a", "yield", f"{DOC}:3a", f"{yield_3a}%",
                {"number": yield_3a, "unit": "%"}, evidence, scope="1a to 3a")
    if extra_product:
        outcome(r, f"{r['reaction_id']}:sideproducts", "other", extra_product,
                "mixture of selenated products", "mixture of selenated products", evidence)
    classify_oxygenative(r, evidence, "Labels cover the observed 2a branch only; structures and bond changes for selenated products or 3a are unresolved.")
    r["notes"] = [F("All quantities and arrow conditions are transcribed from Scheme 2.",
                     "All quantities and arrow conditions are transcribed from Scheme 2.", evidence)]
    r["review"]["unresolved"] = ["Reaction scale, workup, purification, atmosphere, and most yield measurement bases are unavailable without SI."]
    reactions.append(r)
    return r


scheme2_record("scheme2a_100mol", 100, 10, 25, extra_product=f"{DOC}:selenated_mixture")
scheme2_record("scheme2b_10mol", 10, 24, 12, yield_3a=82)
scheme2_record("scheme2c_20mol_2a", 20, 6, 69, scope_2a=True)
scheme2_record("scheme2d_50mol", 50, 6, 72)


# Prose-only catalyst replacement control.
r = base_reaction("phsecl_20mol_control", 2, categories=["optimization", "control"], location="page 2 prose below Scheme 2")
r["general_scheme_ids"] = [scheme1d["reaction_id"]]
r["reaction_name"] = F("replacement of (PhSe)2 by PhSeCl", "selenium catalyst replacement control", ev(2))
r["transformation_type"] = F("product 2a formation", "allene-to-enone transformation", ev(2))
participant(r, f"{r['reaction_id']}:substrate", f"{DOC}:1a", "reactant", ev(2))
participant(r, f"{r['reaction_id']}:product", f"{DOC}:2a", "product", ev(2), product_role="reported product")
participant(r, f"{r['reaction_id']}:mixture", f"{DOC}:complicated_mixture", "product", ev(2), product_role="unresolved mixture")
participant(r, f"{r['reaction_id']}:phsecl", f"{DOC}:phsecl", "catalyst", ev(2), mol_percent=20)
outcome(r, f"{r['reaction_id']}:yield", "yield", f"{DOC}:2a", "14%", {"number": 14, "unit": "%"}, ev(2), scope="1a to 2a")
outcome(r, f"{r['reaction_id']}:mixture", "other", f"{DOC}:complicated_mixture", "a complicated mixture was obtained", "complicated mixture", ev(2))
classify_oxygenative(r, ev(2), "Classification covers only reported 2a; the complicated mixture is chemically unresolved.")
r["review"]["unresolved"] = ["The prose gives only catalyst replacement/loading and 2a yield; time, temperature, solvents, oxidant charge, and yield basis were not imported from another experiment."]
reactions.append(r)


scope3 = {
    "2b": (61, "92:8"), "2c": (72, "90:10"), "2d": (76, "95:5"),
    "2e": (76, "93:7"), "2f": (61, ">99:1"), "2g": (53, "96:4"),
    "2h": (76, ">99:1"), "2i": (49, "96:4"), "2j": (46, "91:9"),
    "2k": (40, "84:16"), "2l": (52, "78:22"), "2m": (43, "89:11"),
    "2n": (65, ">99:1"),
}


def scope_record(label, yield_percent, ez_ratio, page, generic_id):
    substrate = "1" + label[1:]
    evidence = ev(page)
    r = base_reaction(f"scheme{page + 1}_{label}" if page == 2 else f"scheme4_{label}", page,
                      categories=["substrate_scope"], location=f"Scheme {3 if page == 2 else 4}, {label}")
    r["general_scheme_ids"] = [generic_id]
    r["reaction_name"] = F("controlled synthesis of α,β-unsaturated α′-alkoxy ketones",
                            "controlled synthesis of α,β-unsaturated α′-alkoxy ketones", evidence)
    r["transformation_type"] = F("oxygenative rearrangement", "allene-to-enone oxygenative rearrangement", evidence)
    participant(r, f"{r['reaction_id']}:substrate", f"{DOC}:{substrate}", "reactant", evidence)
    participant(r, f"{r['reaction_id']}:product", f"{DOC}:{label}", "product", evidence, product_role="reported scope product")
    participant(r, f"{r['reaction_id']}:catalyst", f"{DOC}:diphenyl_diselenide", "catalyst", evidence, mol_percent=20)
    participant(r, f"{r['reaction_id']}:oxidant", f"{DOC}:pyfotf", "reagent", evidence)
    participant(r, f"{r['reaction_id']}:water", f"{DOC}:water", "reagent", ev(3, "text"))
    set_solvents_and_step(r, page, 6)
    outcome(r, f"{r['reaction_id']}:yield", "yield", f"{DOC}:{label}", f"{yield_percent}%",
            {"number": yield_percent, "unit": "%"}, evidence, basis="isolated yield", yield_type="isolated",
            scope=f"{substrate} to {label}")
    outcome(r, f"{r['reaction_id']}:ez", "selectivity", f"{DOC}:{label}", f"E/Z = {ez_ratio}", ez_ratio,
            evidence, scope={"metric": "E/Z ratio", "order": ["E", "Z"]})
    classify_oxygenative(r, evidence)
    r["review"]["unresolved"] = ["Substrate charge/scale, oxidant amount, E/Z analytical method, workup, purification, and atmosphere are unreported on the supplied page."]
    reactions.append(r)
    return r


for label, values in scope3.items():
    scope_record(label, *values, 2, scheme3_generic["reaction_id"])


# Nonproductive alkyl-substituted precursor 1o.
r = base_reaction("scheme3_1o_no_enone", 2, categories=["substrate_scope", "control"], location="Scheme 3, 1o")
r["general_scheme_ids"] = [scheme3_generic["reaction_id"]]
r["reaction_name"] = F("non-productive precursor", "failed scope experiment", ev(2))
r["transformation_type"] = F("enone formation could not be achieved", "no observed enone formation", ev(2))
participant(r, f"{r['reaction_id']}:substrate", f"{DOC}:1o", "reactant", ev(2))
participant(r, f"{r['reaction_id']}:catalyst", f"{DOC}:diphenyl_diselenide", "catalyst", ev(2), mol_percent=20)
participant(r, f"{r['reaction_id']}:oxidant", f"{DOC}:pyfotf", "reagent", ev(2))
set_solvents_and_step(r, 2, 6)
outcome(r, f"{r['reaction_id']}:no_reaction", "no_reaction", f"{DOC}:1o", "non-productive precursor", "non-productive precursor", ev(2))
classify_none(r, ev(2), "No product transformation is reported; expected bond formation is not assigned.")
r["review"]["unresolved"] = ["The exact 1o structure is retained by source label without a safe complete SMILES; conversion, recovery, and detection basis are unreported."]
reactions.append(r)


# Scheme 4 special cases.
r = base_reaction("scheme4_1p_3b", 3, categories=["substrate_scope"], location="Scheme 4, 1p to 3b")
r["general_scheme_ids"] = [scheme4_generic["reaction_id"]]
r["reaction_name"] = F("allenol reaction", "allenol-to-enone reaction", ev(3))
r["transformation_type"] = F("allenol 1p yielded enone 3b", "allene rearrangement to enone 3b", ev(3))
participant(r, f"{r['reaction_id']}:substrate", f"{DOC}:1p", "reactant", ev(3))
participant(r, f"{r['reaction_id']}:product", f"{DOC}:3b", "product", ev(3), product_role="reported product")
participant(r, f"{r['reaction_id']}:catalyst", f"{DOC}:diphenyl_diselenide", "catalyst", ev(3), mol_percent=20)
participant(r, f"{r['reaction_id']}:oxidant", f"{DOC}:pyfotf", "reagent", ev(3))
set_solvents_and_step(r, 3, 24)
outcome(r, f"{r['reaction_id']}:yield", "yield", f"{DOC}:3b", "44%", {"number": 44, "unit": "%"}, ev(3), scope="1p to 3b")
classify_unresolved(r, ev(3), "The 1p-to-3b branch lacks a verified atom map and does not retain the α′-alkoxy group; bond classification remains unresolved.")
r["review"]["unresolved"] = ["Yield basis, scale, workup, purification, atmosphere, and exact atom correspondence are unreported."]
reactions.append(r)

r = base_reaction("scheme4_1q_recovery", 3, categories=["substrate_scope", "control"], location="Scheme 4, 1q")
r["general_scheme_ids"] = [scheme4_generic["reaction_id"]]
r["reaction_name"] = F("allenyl acetate control", "unreactive allenyl acetate control", ev(3))
r["transformation_type"] = F("starting material 1q was recovered", "no observed transformation", ev(3))
participant(r, f"{r['reaction_id']}:substrate", f"{DOC}:1q", "reactant", ev(3))
participant(r, f"{r['reaction_id']}:catalyst", f"{DOC}:diphenyl_diselenide", "catalyst", ev(3), mol_percent=20)
participant(r, f"{r['reaction_id']}:oxidant", f"{DOC}:pyfotf", "reagent", ev(3))
set_solvents_and_step(r, 3, 6)
outcome(r, f"{r['reaction_id']}:recovery", "recovery", f"{DOC}:1q", "starting material 1q was recovered", "starting material recovered", ev(3))
classify_none(r, ev(3), "Recovered starting material is not a product transformation; no expected labels are assigned.")
r["review"]["unresolved"] = ["Recovery percentage and measurement basis are unreported; no product is invented."]
reactions.append(r)

scope4 = {
    "2r": (57, "92:8"), "2s": (41, "91:9"), "2t": (45, "97:3"),
    "2u": (67, "96:4"), "2v": (31, "88:12"), "2w": (60, "94:6"), "2x": (53, "96:4"),
}
for label, values in scope4.items():
    scope_record(label, *values, 3, scheme4_generic["reaction_id"])


# Scheme 5 isotope-labeling experiment.
r = base_reaction("scheme5_h218o", 3, categories=["mechanistic_study"], location="Scheme 5, 18O-labeling experiment")
r["general_scheme_ids"] = [scheme3_generic["reaction_id"]]
r["reaction_name"] = F("18O-labeling experiment", "oxygen-18 labeling experiment", ev(3))
r["transformation_type"] = F("formation of 18O-labeled 2e", "oxygenative allene rearrangement with isotope labeling", ev(3))
participant(r, f"{r['reaction_id']}:substrate", f"{DOC}:1e", "reactant", ev(3))
participant(r, f"{r['reaction_id']}:product", f"{DOC}:18O-2e", "product", ev(3), product_role="isotope-labeled product")
participant(r, f"{r['reaction_id']}:catalyst", f"{DOC}:diphenyl_diselenide", "catalyst", ev(3), mol_percent=20)
participant(r, f"{r['reaction_id']}:oxidant", f"{DOC}:pyfotf", "reagent", ev(3), equivalents=1.2)
participant(r, f"{r['reaction_id']}:h218o", f"{DOC}:h218o", "reagent", ev(3), equivalents=3)
set_solvents_and_step(r, 3, 6, anhydrous=True)
outcome(r, f"{r['reaction_id']}:incorporation", "other", f"{DOC}:18O-2e", "18O 90%", {"number": 90, "unit": "%"}, ev(3),
        basis="high-resolution mass spectrometry", scope={"metric": "18O incorporation", "site": "carbonyl oxygen"})
classify_oxygenative(r, ev(3), "The isotope source and carbonyl label are source-shown, but exact full atom mapping and proton transfers remain unresolved.")
r["notes"] = [F("18O-labeled 2e was observed by high-resolution mass spectrometry.",
                  "18O-labeled 2e was observed by high-resolution mass spectrometry.", ev(3, "text"))]
r["review"]["unresolved"] = ["No product yield, reaction scale, workup, purification, or atmosphere is reported on the supplied page."]
reactions.append(r)


# Scheme 5 oxygen-balloon control.
r = base_reaction("scheme5_o2_control", 3, categories=["control", "mechanistic_study"], location="Scheme 5, O2-balloon control")
r["general_scheme_ids"] = [scheme3_generic["reaction_id"]]
r["reaction_name"] = F("O2-balloon control", "oxygen-atmosphere control", ev(3))
r["transformation_type"] = F("formation of 2a in dehydrated solvents under O2", "allene-to-enone control", ev(3))
participant(r, f"{r['reaction_id']}:substrate", f"{DOC}:1a", "reactant", ev(3))
participant(r, f"{r['reaction_id']}:product", f"{DOC}:2a", "product", ev(3), product_role="reported product")
participant(r, f"{r['reaction_id']}:catalyst", f"{DOC}:diphenyl_diselenide", "catalyst", ev(3), mol_percent=20)
participant(r, f"{r['reaction_id']}:oxidant", f"{DOC}:pyfotf", "reagent", ev(3), equivalents=1.2)
participant(r, f"{r['reaction_id']}:oxygen", f"{DOC}:oxygen", "reagent", ev(3))
set_solvents_and_step(r, 3, 6, anhydrous=True, oxygen=True)
outcome(r, f"{r['reaction_id']}:yield", "yield", f"{DOC}:2a", "10%", {"number": 10, "unit": "%"}, ev(3), scope="1a to 2a")
classify_oxygenative(r, ev(3))
r["review"]["unresolved"] = ["Yield measurement basis, reaction scale, workup, and purification are unreported."]
reactions.append(r)


mechanism = add_generic("scheme6_proposed_mechanism", 3, "proposed_mechanism", "Scheme 6",
                        "rationalization for formation of α,β-unsaturated α′-alkoxy ketones 2",
                        "Author-proposed cycle via PhSe(OTf), seleniranium INT-1, INT-2 to INT-5, and regeneration of selenium π-acid catalyst.")
mechanism["evidence_ids"] = list(dict.fromkeys(mechanism["evidence_ids"] + ev(6)))
mechanism["notes"].append(F(
    "Footnote 14 describes an alternative mechanism in Scheme S1 but states the authors do not think it is operative for 2a–n and 2r–x.",
    "Alternative SI mechanism mentioned but not supplied or adopted as operative by the authors.", ev(6),
))


bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
bundle["source_context"] = [
    {
        "context_id": f"{DOC}:figure1_computational_profile",
        "document_id": DOC,
        "kind": "computational_mechanistic_context",
        "evidence_ids": ev(4),
        "description": "Figure 1 is a computed 1b + PhSeOTf -> INT-1b profile, not a synthetic experiment. Reported values include ΔG‡ 4.6 kcal/mol, ΔGR -15.6 kcal/mol, and Se···O distances 3.199 Å (TS) and 3.029 Å (INT-1b).",
    },
    {
        "context_id": f"{DOC}:si_notice",
        "document_id": DOC,
        "kind": "supplementary_information_notice",
        "evidence_ids": ev(4),
        "description": "The article states SI is available, but no SI file was explicitly provided or associated in this session.",
    },
    {
        "context_id": f"{DOC}:references_pages5_6",
        "document_id": DOC,
        "kind": "bibliography_and_mechanistic_footnote",
        "evidence_ids": ev(5) + ev(6),
        "description": "Pages 5–6 are references except for footnotes 9–14; no current-paper synthetic experiment is created from bibliography entries.",
    },
]


# Save raw transcription separately from the formal bundle.
raw = {
    "document_id": DOC,
    "title": TITLE,
    "pages_attempted": [1, 2, 3, 4, 5, 6],
    "record_count": len(reactions),
    "specific_experiment_count": sum(r["record_type"] == "specific_experiment" for r in reactions),
    "records": [
        {
            "reaction_id": r["reaction_id"],
            "record_type": r["record_type"],
            "location": r["location"]["raw_value"],
            "evidence_ids": r["evidence_ids"],
            "status": "raw_model_transcription",
        }
        for r in reactions
    ],
    "limitations": [
        "No supporting-information file was supplied.",
        "Candidate SMILES are unverified structure transcriptions.",
        "Formal records remain needs_review after machine merge.",
    ],
}
st.write(RUN / "raw_model" / "pages_0001_0006_extraction.json", raw)
st.write(RUN / "raw_model" / "stereo_checks.json", stereo_checks)
mentions = st.source_mentions(dataset, ["1a", "2a", "1o", "1p", "1q", "2e", "3a", "3b"], context=220)
st.write(RUN / "raw_model" / "source_mentions.json", mentions)
st.write(RUN / "submitted.json", bundle)
print(json.dumps({"molecules": len(molecules), "reactions": len(reactions),
                  "specific_experiments": raw["specific_experiment_count"]}, ensure_ascii=False, indent=2))
