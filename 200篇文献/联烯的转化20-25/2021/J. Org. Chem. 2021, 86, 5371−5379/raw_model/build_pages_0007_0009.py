#!/usr/bin/env python3
"""Enrich the v3.1 extraction with source evidence from PDF pages 7--9.

This script deliberately updates already-indexed Scheme 2 / eq 2--4 records
instead of creating duplicate experiments.  Existing records provide IDs and
schema scaffolding only; all new facts below come from visual/text inspection
of pages 7--9.  Page 3 is used only as explicitly identified cross-page context
for the printed reaction arrows and the conflicting "standard conditions".
"""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_2c129c777a77980b"
RUN = ROOT / "runs/current"
DATASET = RUN / "dataset.json"
STANDARD = RUN / "standard_records.json"
OUT = RUN / "submitted_pages_0007_0009.json"
CANONICAL_SUBMITTED = RUN / "submitted.json"
PRIOR_SUBMITTED = RUN / "submitted_pages_0001_0006.json"
RAW = RUN / "raw_model"
COVERAGE = RUN / "coverage/pages_0007_0009.json"


def ev(page: int, kind: str = "text") -> str:
    return f"{DOC}:page{page}:{kind}"


def evidence(*pages: int) -> list[str]:
    return [item for page in pages for item in (ev(page), ev(page, "image"))]


def F(raw=None, normalized=None, evs=(), method="original_extraction",
      status=None, verification=None, derivation=None, alternatives=()):
    if normalized is None and raw is not None and status != "conflict":
        normalized = raw
    return st.fact(raw, normalized, evs, method=method, status=status,
                   verification=verification, derivation=derivation,
                   alternatives=alternatives)


def q(number, unit, **extra):
    return {"number": number, "unit": unit, **extra}


def add_evidence(obj, *ids):
    for item in ids:
        if item not in obj:
            obj.append(item)


dataset = json.loads(DATASET.read_text())
current = json.loads(STANDARD.read_text())
bundle = json.loads(OUT.read_text())
bundle["molecules"] = []
bundle["reactions"] = []
bundle["source_context"] = []

old_molecules = {m["molecule_id"]: m for m in current["molecules"]}
old_reactions = {r["reaction_id"]: r for r in current["reactions"]}


def mid(label: str) -> str:
    return f"{DOC}:mol_{label}"


def rid(suffix: str) -> str:
    return f"{DOC}:{suffix}"


def get_reaction(suffix: str):
    return copy.deepcopy(old_reactions[rid(suffix)])


def make_participant(pid, molecule_id, role, stage="reaction", evs=()):
    return st.new_participant(pid, molecule_id, role, evs, stage)


def set_charge(participant, key, raw, normalized, evs, method="original_extraction",
               status=None, alternatives=()):
    participant["charges"][key] = F(raw, normalized, evs, method=method,
                                      status=status, alternatives=alternatives)


def make_step(sid, stage, pids, operation, evs, temperature=None, time=None,
              atmosphere=None):
    s = st.new_step(sid, evs, stage)
    s["participant_ids"] = list(pids)
    s["operation"] = F(operation, operation, evs)
    if temperature is not None:
        raw, normalized = temperature
        s["temperature"] = F(raw, normalized, evs)
    if time is not None:
        raw, normalized = time
        s["time"] = F(raw, normalized, evs)
    if atmosphere is not None:
        s["atmosphere"] = F(atmosphere, atmosphere, evs)
    return s


def make_outcome(oid, kind, targets, raw, normalized, evs, basis=None,
                 yield_type=None, scope=None):
    o = st.new_outcome(oid, kind, targets, evs)
    o["value"] = F(raw, normalized, evs)
    if basis is not None:
        o["measurement_basis"] = F(basis, basis, evs)
    if yield_type is not None:
        o["yield_type"] = F(yield_type, yield_type, evs)
    if scope is not None:
        o["scope"] = F(scope, scope, evs)
    return o


def note(text, evs, method="original_extraction", derivation=None):
    return F(text if method == "original_extraction" else None, text, evs,
             method=method, derivation=derivation)


def update_location(r, pages, evs):
    r["location"] = F({"pdf_pages": pages}, {"pdf_pages": pages}, evs)
    r["evidence_ids"] = list(evs)
    r["review"]["status"] = "needs_review"
    r["review"]["method"] = "model"
    r["review"]["reviewer"] = None
    r["review"]["checked_fields"] = []
    r["completeness"] = "partial"


def append_unresolved(r, text):
    if text not in r["review"]["unresolved"]:
        r["review"]["unresolved"].append(text)


def new_named_molecule(label, name, evs, smiles=None, chirality="unresolved",
                       geometry="not_applicable", representation="not_encoded",
                       identity_note=None):
    molecule_id = mid(label)
    m = st.new_molecule(molecule_id, DOC, evs)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, evs)
    m["names"]["original"] = F(name, name, evs)
    m["names"]["english"] = F(name, name, evs)
    m["names"]["abbreviations"] = F(label, [label], evs)
    m["names"]["aliases"] = F()
    m["cas_numbers"] = F()
    m["external_properties"] = []
    m["generic_template"] = None
    m["identity_note"] = F(identity_note, identity_note, evs) if identity_note else F()
    if smiles:
        diag = st.stereo_inspect(smiles)
        canonical = diag["canonical_isomeric_smiles"]
        m["representation_type"] = "smiles"
        m["smiles"] = F(smiles, canonical, evs, method="structure_derivation",
                         derivation="Candidate graph encoded from the source name and drawing; graph parsing does not establish source identity.")
        m["stereochemistry"]["isomeric_smiles"] = F(smiles, canonical, evs,
            method="structure_derivation", derivation="Canonical isomeric SMILES from the candidate graph.")
        m["stereochemistry"]["identity_key"] = F(None, diag["stereo_graph_key"], evs,
            method="structure_derivation", derivation="standard stereo graph key for the candidate isomeric SMILES.")
        m["stereochemistry"]["atom_centers"] = F(None, diag["atom_centers"], evs,
            method="structure_derivation", derivation="Potential/assigned centers in the candidate graph; '?' remains unassigned.")
        m["stereochemistry"]["double_bonds"] = F("source name and candidate graph", diag["double_bonds"], evs)
        m["stereochemistry"]["representation_status"] = F(None, representation, evs,
            method="inference", derivation="The E alkene is encoded; the newly formed oxirane stereocenter is unassigned in the source and candidate graph.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = F()
        m["stereochemistry"]["atom_centers"] = F()
        m["stereochemistry"]["double_bonds"] = F()
        m["stereochemistry"]["isomeric_smiles"] = F()
        m["stereochemistry"]["identity_key"] = F()
        m["stereochemistry"]["representation_status"] = F(None, representation, evs,
            method="inference", derivation="Named material is retained without inventing a unique molecular graph.")
    m["stereochemistry"]["axial"] = F()
    m["stereochemistry"]["mixture"] = F()
    m["stereochemistry"]["assignment"] = F()
    m["stereochemistry"]["chirality_status"] = F(None, chirality, evs,
        method="inference", derivation="Unspecified is not treated as achiral or racemic.")
    m["stereochemistry"]["geometry_status"] = F(None, geometry, evs,
        method="inference", derivation="Geometry status follows the source name/drawing.")
    return m


# ---------------------------------------------------------------------------
# Page 7: enrich the twelve existing Scheme 2 experiments, one per product.
# ---------------------------------------------------------------------------
page7_evs = evidence(3, 7)
page7_product_data = {
    "3fa": ("Benzyl (E)-4-methyl-2-methylene-3-styrylpent-3-enoate", 33.3, 35),
    "3ga": ("Benzyl (E)-6-(dimethylamino)-2-methylene-6-oxo-3-(propan-2-ylidene)hex-4-enoate", 62.8, 67),
    "3gi": ("Ethyl (E)-6-(dimethylamino)-2-methylene-6-oxo-3-(propan-2-ylidene)hex-4-enoate", 67.1, 54),
    "3ha": ("Benzyl (E)-3-(3-(dimethylamino)-3-oxoprop-1-en-1-yl)-4-ethyl-2-methylenehex-3-enoate", 64.6, 63),
    "3ia": ("Benzyl (E)-3-cyclohexylidene-6-(dimethylamino)-2-methylene-6-oxohex-4-enoate", 79.9, 75),
    "3ib": ("[1,1′-Biphenyl]-4-ylmethyl (E)-3-cyclohexylidene-6-(dimethylamino)-2-methylene-6-oxohex-4-enoate", 87.2, 68),
    "3ie": ("4-(Trifluoromethyl)benzyl (E)-3-cyclohexylidene-6-(dimethylamino)-2-methylene-6-oxohex-4-enoate", 93.3, 74),
    "3ii": ("Ethyl (E)-3-cyclohexylidene-6-(dimethylamino)-2-methylene-6-oxohex-4-enoate", 103.8, 71),
    "3ja": ("6-Benzyl 1-ethyl (2E,4E)-4-ethylidene-5-methylenehex-2-enedioate", 71.9, 80),
    "3ik": ("6-(((3r,5r,7r)-Adamantan-1-yl)methyl) 1-ethyl (E)-4-cyclohexylidene-5-methylenehex-2-enedioate", 84.5, 69),
    "3im": ("Naphthalen-1-ylmethyl (E)-3-cyclohexylidene-6-(dimethylamino)-2-methylene-6-oxohex-4-enoate", 79.1, 65),
    "3in": ("Ethyl (4E)-3-cyclohexylidene-6-(dimethylamino)-2-ethylidene-6-oxohex-4-enoate", 33.3, 36),
}

modified_molecules = {}
updated_reactions = []

for label, (source_name, mass_mg, yield_pct) in page7_product_data.items():
    molecule_id = mid(label)
    m = copy.deepcopy(old_molecules[molecule_id])
    occ = m["occurrences"][0]
    add_evidence(occ["evidence_ids"], ev(7), ev(7, "image"))
    occ["compound_label"] = F(label, label, evidence(7))
    m["names"]["original"] = F(source_name, source_name, evidence(7))
    m["names"]["english"] = F(source_name, source_name, evidence(7))
    m["names"]["abbreviations"] = F(label, [label], evidence(7))
    m["identity_status"] = "source_identified"
    m["representation_type"] = "unresolved"
    m["smiles"] = F()
    m["stereochemistry"]["atom_centers"] = F()
    m["stereochemistry"]["axial"] = F()
    m["stereochemistry"]["isomeric_smiles"] = F()
    m["stereochemistry"]["identity_key"] = F()
    if label == "3ja":
        m["stereochemistry"]["double_bonds"] = F("(2E,4E) product name; Scheme 2 E/Z = 2:1", [{"assignment": "2E,4E", "applies_to": "named major component"}], page7_evs)
        m["stereochemistry"]["mixture"] = F("E/Z = 2:1", {"metric": "E/Z", "ratio": "2:1", "order": ["E", "Z"], "method": None}, evidence(3, 7))
        m["stereochemistry"]["geometry_status"] = F(None, "mixture", page7_evs, method="inference", derivation="Scheme 2 reports E/Z = 2:1; page-7 name labels the named component (2E,4E).")
        identity_note = "The page-7 product block contains fractional NMR integrations consistent with the reported E/Z mixture; no separate isomer structures were assigned."
    elif label == "3in":
        m["stereochemistry"]["double_bonds"] = F("(4E) product name; Scheme 2 Z/E = 2.2/1", [{"assignment": "4E", "element": "hex-4-ene"}], page7_evs)
        m["stereochemistry"]["mixture"] = F("Z/E = 2.2/1", {"metric": "Z/E", "ratio": "2.2/1", "order": ["Z", "E"], "method": None}, evidence(3, 7))
        m["stereochemistry"]["geometry_status"] = F(None, "mixture", page7_evs, method="inference", derivation="Scheme 2 reports the ethylidene Z/E mixture while the product name fixes the separate 4E alkene.")
        identity_note = "The page-7 NMR block has fractional integrations for the mixture; the ethylidene isomer structures were not independently encoded."
    else:
        m["stereochemistry"]["double_bonds"] = F("(E) in product name", [{"assignment": "E", "basis": "page-7 product name"}], evidence(7))
        m["stereochemistry"]["mixture"] = F()
        m["stereochemistry"]["geometry_status"] = F(None, "specified", evidence(7), method="inference", derivation="The page-7 product name explicitly assigns E geometry.")
        identity_note = "The full structure is visible in Scheme 2 and named on page 7, but no machine SMILES is asserted pending expert structure review."
    m["stereochemistry"]["chirality_status"] = F(None, "partially_specified" if label == "3ik" else "unresolved", evidence(7), method="inference", derivation="No unreported absolute/relative configuration is inferred; the printed (3r,5r,7r) adamantyl descriptor is retained only for 3ik.")
    m["stereochemistry"]["representation_status"] = F(None, "not_encoded", evidence(7), method="inference", derivation="Source geometry is preserved textually but no complete graph is encoded.")
    m["stereochemistry"]["assignment"] = F("(3r,5r,7r) in printed 3ik name", {"descriptor": "(3r,5r,7r)", "scope": "adamantyl substituent in printed name"}, evidence(7)) if label == "3ik" else F()
    if label == "3ik":
        identity_note += " The printed name omits the dimethylamide motif visible in Scheme 2 and indicated by NMR/HRMS; the source conflict is not silently repaired."
    m["identity_note"] = F(identity_note, identity_note, evidence(7))
    modified_molecules[molecule_id] = m

    r = get_reaction(f"scheme2_scope_{label}")
    update_location(r, [3, 7], page7_evs)
    yield_outcome = next(o for o in r["outcomes"] if o["kind"] == "yield")
    yield_outcome["value"] = F(f"{yield_pct}%", q(yield_pct, "%"), page7_evs)
    yield_outcome["measurement_basis"] = F("isolated product", "isolated product", page7_evs)
    yield_outcome["yield_type"] = F("isolated yield", "isolated yield", page7_evs)
    r["outcomes"] = [o for o in r["outcomes"] if o["kind"] != "product_amount"]
    r["outcomes"].append(make_outcome(f"amount_{label}", "product_amount", [molecule_id],
        f"{mass_mg} mg", q(mass_mg, "mg"), evidence(7), basis="isolated product after chromatography"))
    append_unresolved(r, "Page-7 NMR/HRMS/IR/Rf data are characterization evidence, not additional synthesis experiments.")
    if label in {"3ja", "3in"}:
        append_unresolved(r, "The selectivity ratio is source-reported but its analytical measurement method is not stated.")
    if label == "3ik":
        append_unresolved(r, "The printed page-7 name conflicts with the Scheme-2/NMR/HRMS indication of a dimethylamide-containing product; no corrected structure was invented.")
    updated_reactions.append(r)


# ---------------------------------------------------------------------------
# Page 8: product 3' characterization (context only), epoxidation, and D study.
# ---------------------------------------------------------------------------
page8_evs = evidence(8)

# Update the dimer registry entry, but do not invent a new synthesis record from
# its analytical block.  Table 1 already carries every experiment that reports it.
m3p = copy.deepcopy(old_molecules[mid("3prime")])
add_evidence(m3p["occurrences"][0]["evidence_ids"], *page8_evs)
m3p["names"]["original"] = F("Diethyl (2E,6E)-4,5-di(propan-2-ylidene)octa-2,6-dienedioate", "Diethyl (2E,6E)-4,5-di(propan-2-ylidene)octa-2,6-dienedioate", page8_evs)
m3p["names"]["english"] = copy.deepcopy(m3p["names"]["original"])
m3p["stereochemistry"]["double_bonds"] = F("(2E,6E) in product name", [{"assignment": "2E,6E"}], page8_evs)
m3p["stereochemistry"]["geometry_status"] = F(None, "specified", page8_evs, method="inference", derivation="Both terminal enoate geometries are named (2E,6E).")
m3p["stereochemistry"]["chirality_status"] = F(None, "unresolved", page8_evs, method="inference", derivation="No chirality assignment is made; unknown is not converted to achiral.")
m3p["stereochemistry"]["representation_status"] = F(None, "not_encoded", page8_evs, method="inference", derivation="The named geometry is retained without a complete candidate graph.")
m3p["identity_note"] = F("Page 8 provides characterization only; the experiment(s) producing 3′ remain the distinct Table 1 rows on page 2.", "Page 8 provides characterization only; the experiment(s) producing 3′ remain the distinct Table 1 rows on page 2.", page8_evs)
modified_molecules[mid("3prime")] = m3p

# Product 4 receives an evidence-linked candidate graph.  Standard stereo
# diagnostics identify the oxirane center as potential but unassigned.
product4_smiles = "CCOC(=O)/C=C/C1(OC1(C)C)C(=C)C(=O)OCc1ccccc1"
m4 = new_named_molecule("4", "Benzyl (E)-2-(2-(3-ethoxy-3-oxoprop-1-en-1-yl)-3,3-dimethyloxiran-2-yl)acrylate",
    evidence(3, 8), smiles=product4_smiles, chirality="unspecified",
    geometry="specified", representation="partial",
    identity_note="The source assigns E geometry but no absolute configuration at the newly formed oxirane stereocenter; the candidate graph leaves that center unspecified.")
m4["stereochemistry"]["assignment"] = F("(E) in product name", {"double_bond": "3-ethoxy-3-oxoprop-1-en-1-yl", "assignment": "E"}, evidence(8))
modified_molecules[mid("4")] = m4

quench_mid = mid("sat_NaHSO3_Na2S2O5")
water_mid = mid("DI_water")
modified_molecules[quench_mid] = new_named_molecule("sat_NaHSO3_Na2S2O5", "saturated NaHSO3/Na2S2O5 solution", page8_evs,
    chirality="not_applicable", identity_note="The slash notation is preserved as a source-named workup solution; its exact composition is not uniquely specified.")
modified_molecules[water_mid] = new_named_molecule("DI_water", "deionized water", page8_evs, smiles="O",
    chirality="achiral", geometry="not_applicable", representation="encoded")

epox = get_reaction("eq2_epoxidation_3aa_to_4")
epox_evs = evidence(3, 8)
update_location(epox, [3, 8], epox_evs)
epox["participants"] = []
p = make_participant("substrate", mid("3aa"), "reactant", evs=page8_evs)
set_charge(p, "mass", "62.8 mg", q(62.8, "mg"), page8_evs)
set_charge(p, "amount", "0.2 mmol", q(0.2, "mmol"), page8_evs)
set_charge(p, "equivalents", "1.0 equiv", q(1.0, "equiv"), page8_evs)
p["equivalence_reference"] = F(mid("3aa"), mid("3aa"), page8_evs)
epox["participants"].append(p)
p = make_participant("mcpba", mid("m_CPBA"), "reagent", evs=page8_evs)
set_charge(p, "equivalents", "1.8 equiv", q(1.8, "equiv"), page8_evs)
p["equivalence_reference"] = F(mid("3aa"), mid("3aa"), page8_evs)
epox["participants"].append(p)
p = make_participant("dcm_reaction", mid("DCM"), "solvent", evs=page8_evs)
set_charge(p, "volume", "2 mL", q(2, "mL"), page8_evs)
epox["participants"].append(p)
epox["participants"].append(make_participant("product", mid("4"), "product", evs=page8_evs))
epox["participants"].append(make_participant("quench", quench_mid, "workup_agent", "workup", page8_evs))
p = make_participant("dcm_extract", mid("DCM"), "workup_agent", "workup", page8_evs)
set_charge(p, "volume", "5 mL × 3", {"number": 5, "unit": "mL", "repeat_count": 3}, page8_evs)
epox["participants"].append(p)
epox["participants"].append(make_participant("water_wash", water_mid, "workup_agent", "workup", page8_evs))
epox["participants"].append(make_participant("na2so4", mid("Na2SO4"), "workup_agent", "workup", page8_evs))
epox["participants"].append(make_participant("silica", mid("silica"), "purification_agent", "purification", page8_evs))
epox["participants"].append(make_participant("ea", mid("EA"), "purification_agent", "purification", page8_evs))
epox["participants"].append(make_participant("pe", mid("PE"), "purification_agent", "purification", page8_evs))
epox["reaction_scale"] = F("3aa (0.2 mmol, 62.8 mg)", {"reference_molecule_id": mid("3aa"), "amount": q(0.2, "mmol")}, page8_evs)
epox["solvent_conditions"]["composition"] = F("dry DCM", [mid("DCM")], page8_evs)
epox["solvent_conditions"]["total_volume"] = F("2 mL", q(2, "mL"), page8_evs)
epox["procedure"]["original_text"] = F("Transformation of [3]Dendralene 3aa with m-CPBA", "page-8 experimental procedure", page8_evs)
epox["procedure"]["addition_order"] = F("3aa in dry DCM; cool to 0 °C; add m-CPBA in one portion", ["3aa in dry DCM", "cool to 0 °C", "m-CPBA in one portion"], page8_evs)
epox["procedure"]["pretreatment"] = F("10 mL oven-dried reaction flask", "10 mL oven-dried reaction flask", page8_evs)
epox["steps"] = [
    make_step("charge_and_cool", "reaction", ["substrate", "dcm_reaction"], "charge 3aa in dry DCM and cool with stirring", page8_evs, ("0 °C", q(0, "°C"))),
    make_step("mcpba_addition", "reaction", ["mcpba"], "add m-CPBA (85% pure) in one portion", page8_evs, ("0 °C", q(0, "°C"))),
    make_step("cold_hold", "reaction", ["substrate", "mcpba", "dcm_reaction"], "stir", page8_evs, ("0 °C", q(0, "°C")), ("20 min", q(20, "min"))),
    make_step("warm_hold", "reaction", ["substrate", "mcpba", "dcm_reaction"], "allow to warm from 0 °C to room temperature; transition time not reported", page8_evs, ("rt", "room temperature"), ("overnight", "overnight")),
    make_step("workup", "workup", ["quench", "dcm_extract", "water_wash", "na2so4"], "dilute with DCM and saturated NaHSO3/Na2S2O5; transfer, extract with DCM (5 mL × 3), wash combined organics with D.I. water, dry over anhydrous Na2SO4, filter salts, and rotary evaporate", page8_evs, ("room temperature", "room temperature")),
    make_step("purification", "purification", ["silica", "ea", "pe", "product"], "silica-gel column chromatography, EA:PE = 1:10", page8_evs, ("room temperature", "room temperature")),
]
epox["outcomes"] = [
    make_outcome("yield_4", "yield", [mid("4")], "50%", q(50, "%"), epox_evs,
                 basis="isolated product after silica-gel chromatography", yield_type="isolated yield"),
    make_outcome("amount_4", "product_amount", [mid("4")], "32.9 mg", q(32.9, "mg"), page8_evs,
                 basis="isolated colorless oil after silica-gel chromatography"),
]
epox["notes"] = [note("m-CPBA was reported as 85% pure; the source gives 1.8 equivalents but no mass or amount.", page8_evs)]
epox["review"]["unresolved"] = [
    "m-CPBA mass/amount is unreported; only 1.8 equivalents and 85% purity are stated.",
    "The transition time from 0 °C to room temperature is unreported; only the overnight room-temperature stage is retained.",
    "The exact composition of the saturated NaHSO3/Na2S2O5 workup solution is not specified.",
    "The source assigns E geometry but no absolute configuration at the newly formed oxirane stereocenter.",
]
updated_reactions.append(epox)


# The no-AcOH run keeps its page-3 standard-condition arrow.  The detailed
# page-8 procedure explicitly contains AcOH and therefore is not copied into it.
no_acoh = get_reaction("eq3_deuterium_no_acoh")
no_acoh_evs = evidence(3, 8)
update_location(no_acoh, [3, 8], no_acoh_evs)
no_acoh["outcomes"] = [
    make_outcome("yield_d6_product", "yield", [mid("3aa_d6")], "78% without HOAc", q(78, "%"), no_acoh_evs,
                 basis="isolated product reported in the combined page-8 product block", yield_type="isolated yield"),
    make_outcome("amount_d6_product", "product_amount", [mid("3aa_d6")], "50.0 mg without HOAc", q(50.0, "mg"), page8_evs,
                 basis="isolated colorless oil reported in the combined product block"),
    make_outcome("deuterium_scrambling", "other", [mid("3aa_d6")], "without any scrambling of the deuterium content", "without any scrambling of the deuterium content", no_acoh_evs,
                 basis="source statement", scope={"metric": "deuterium retention", "method": "not reported"}),
]
no_acoh["notes"] = [note("The page-8 detailed procedure explicitly adds AcOH; it was not silently transferred to this no-AcOH experiment.", page8_evs)]
no_acoh["review"]["unresolved"] = [
    "The analytical method used to establish absence of deuterium scrambling is not reported in the procedure.",
    "The page-8 product block reports the no-HOAc mass/yield, but its preceding detailed procedure explicitly contains AcOH; full no-HOAc isolation details remain ambiguous.",
]
updated_reactions.append(no_acoh)


# Rebuild the AcOH-containing run from the explicit page-8 procedure while
# retaining the page-3 "standard conditions" conflict as a conflict fact.
with_acoh = get_reaction("eq4_deuterium_with_acoh")
with_acoh_evs = evidence(3, 8)
update_location(with_acoh, [3, 8], with_acoh_evs)
with_acoh["participants"] = []
p = make_participant("allene_d6", mid("1a_d6"), "reactant", evs=page8_evs)
set_charge(p, "amount", "0.3 mmol", q(0.3, "mmol"), page8_evs)
set_charge(p, "equivalents", "1.0 equiv", q(1.0, "equiv"), page8_evs)
p["equivalence_reference"] = F(mid("1a_d6"), mid("1a_d6"), page8_evs)
with_acoh["participants"].append(p)
with_acoh["participants"].append(make_participant("diazo", mid("2a"), "reactant", evs=page8_evs))
p = make_participant("pd", mid("PdOAc2"), "catalyst", evs=with_acoh_evs)
set_charge(p, "mass", "4.4 mg", q(4.4, "mg"), page8_evs)
set_charge(p, "amount", "0.02 mmol", q(0.02, "mmol"), page8_evs)
p["charges"]["mol_percent"] = F(
    {"page3": "standard conditions (Scheme 2: Pd(OAc)2 10 mol%)", "page8": "0.02 mmol Pd(OAc)2 with 0.3 mmol 1a-d6"},
    None, with_acoh_evs, status="conflict", verification="unreviewed",
    alternatives=[
        F("standard conditions (Scheme 2: Pd(OAc)2 10 mol%)", q(10, "mol%"), evidence(3), method="inherited_general_condition"),
        F("0.02 mmol Pd(OAc)2 with 0.3 mmol 1a-d6", {"number": 0.02, "unit": "mmol", "basis": "0.3 mmol 1a-d6"}, page8_evs),
    ])
p["equivalence_reference"] = F(mid("1a_d6"), mid("1a_d6"), with_acoh_evs)
with_acoh["participants"].append(p)
p = make_participant("dmbq", mid("2,6_DMBQ"), "reagent", evs=page8_evs)
set_charge(p, "mass", "30.0 mg", q(30.0, "mg"), page8_evs)
set_charge(p, "amount", "0.22 mmol", q(0.22, "mmol"), page8_evs)
set_charge(p, "equivalents", "1.1 equiv", q(1.1, "equiv"), page8_evs)
p["equivalence_reference"] = F(mid("1a_d6"), mid("1a_d6"), page8_evs)
with_acoh["participants"].append(p)
p = make_participant("thf", mid("THF"), "solvent", evs=page8_evs)
set_charge(p, "volume", "2 mL", q(2, "mL"), page8_evs)
with_acoh["participants"].append(p)
p = make_participant("acoh", mid("AcOH"), "additive", evs=page8_evs)
set_charge(p, "mass", "12.0 mg", q(12.0, "mg"), page8_evs)
set_charge(p, "equivalents", "1.0 equiv", q(1.0, "equiv"), page8_evs)
p["equivalence_reference"] = F(mid("1a_d6"), mid("1a_d6"), page8_evs)
with_acoh["participants"].append(p)
with_acoh["participants"].append(make_participant("product", mid("3aa_d6"), "product", evs=page8_evs))
with_acoh["participants"].append(make_participant("silica", mid("silica"), "purification_agent", "purification", page8_evs))
with_acoh["participants"].append(make_participant("ea", mid("EA"), "purification_agent", "purification", page8_evs))
with_acoh["participants"].append(make_participant("pe", mid("PE"), "purification_agent", "purification", page8_evs))
with_acoh["reaction_scale"] = F("1a-d6 (0.3 mmol)", {"reference_molecule_id": mid("1a_d6"), "amount": q(0.3, "mmol")}, page8_evs)
with_acoh["solvent_conditions"]["composition"] = F("dry degassed THF", [mid("THF")], page8_evs)
with_acoh["solvent_conditions"]["total_volume"] = F("2 mL", q(2, "mL"), page8_evs)
with_acoh["procedure"]["original_text"] = F("Deuterium-Labeling Experiment", "page-8 experimental procedure", page8_evs)
with_acoh["procedure"]["addition_order"] = F("Pd(OAc)2 and 2,6-DMBQ; vacuum/N2 cycles; THF; 1a-d6; 2a; AcOH", ["Pd(OAc)2", "2,6-DMBQ", "vacuum/N2 cycles", "THF", "1a-d6", "2a", "AcOH"], page8_evs)
with_acoh["procedure"]["pretreatment"] = F("10 mL oven-dried flask; vacuum and back-filled with N2 three times; dry degassed THF", "10 mL oven-dried flask and three vacuum/N2 cycles", page8_evs)
with_acoh["steps"] = [
    make_step("pretreat", "reaction", ["pd", "dmbq"], "weigh into 10 mL oven-dried flask; evacuate and back-fill with N2 three times; add dry degassed THF", page8_evs, atmosphere="N2"),
    make_step("addition", "reaction", ["allene_d6", "diazo", "acoh"], "add 1a-d6 and 2a successively by microsyringe, then add AcOH", page8_evs, atmosphere="N2"),
    make_step("reaction", "reaction", ["allene_d6", "diazo", "pd", "dmbq", "thf", "acoh"], "stir in an oil bath", page8_evs, ("90 °C", q(90, "°C")), ("10 h", q(10, "h")), "N2"),
    make_step("post_reaction", "workup", ["silica", "ea", "pe"], "source says the mixture was allowed to warm to room temperature; filter through silica with 1:1 EtOAc/petroleum ether, evaporate, then purify by silica-gel chromatography (EtOAc/petroleum ether; final ratio unreported)", page8_evs, ("room temperature", "room temperature")),
]
with_acoh["outcomes"] = [
    make_outcome("yield_d6_product", "yield", [mid("3aa_d6")], "8% with AcOH", q(8, "%"), with_acoh_evs,
                 basis="isolated product reported in the combined page-8 product block", yield_type="isolated yield"),
    make_outcome("amount_d6_product", "product_amount", [mid("3aa_d6")], "4.9 mg with AcOH", q(4.9, "mg"), page8_evs,
                 basis="isolated colorless oil after silica-gel chromatography"),
    make_outcome("deuterium_scrambling", "other", [mid("3aa_d6")], "without any scrambling of the deuterium content", "without any scrambling of the deuterium content", page8_evs,
                 basis="source statement", scope={"metric": "deuterium retention", "method": "not reported"}),
]
with_acoh["notes"] = [
    note("The page-3 arrow says AcOH (1.0 equiv) under standard conditions, whose Scheme-2 definition uses Pd(OAc)2 (10 mol%); page 8 instead reports Pd(OAc)2 4.4 mg, 0.02 mmol with 1a-d6 0.3 mmol. Both are retained as a source conflict.", with_acoh_evs),
    note("The page-8 source explicitly reports DMBQ as 30.0 mg, 0.22 mmol, 1.1 equiv and 1a-d6 as 0.3 mmol; the numbers are internally inconsistent and were not corrected.", page8_evs),
    note("The phrase 'allowed to warm to room temperature' after a 90 °C hold is preserved verbatim and not changed to 'cool'.", page8_evs),
]
with_acoh["review"]["unresolved"] = [
    "Pd(OAc)2 loading conflicts between page-3 standard conditions (10 mol%) and page-8 quantities (0.02 mmol Pd vs 0.3 mmol 1a-d6); no value was chosen.",
    "Page 8 reports DMBQ as both 0.22 mmol and 1.1 equivalents against 0.3 mmol 1a-d6; these values are internally inconsistent.",
    "The amount of diazo ester 2a is not reported in the page-8 procedure.",
    "The final chromatography EtOAc/petroleum ether ratio is not reported.",
    "The analytical method used to establish absence of deuterium scrambling is not reported.",
]
updated_reactions.append(with_acoh)


# Molecule 3aa-d6: page 8 names the E product and reports retention but does not
# assign any additional stereochemistry or provide an isotopic SMILES.
d6_product = copy.deepcopy(old_molecules[mid("3aa_d6")])
add_evidence(d6_product["occurrences"][0]["evidence_ids"], *page8_evs)
d6_name = "6-Benzyl 1-ethyl (E)-5-methylene-4-(propan-2-ylidene-d6)hex-2-enedioate"
d6_product["names"]["original"] = F(d6_name, d6_name, page8_evs)
d6_product["names"]["english"] = copy.deepcopy(d6_product["names"]["original"])
d6_product["stereochemistry"]["double_bonds"] = F("(E) in product name", [{"assignment": "E", "element": "hex-2-ene"}], page8_evs)
d6_product["stereochemistry"]["chirality_status"] = F(None, "unresolved", page8_evs, method="inference", derivation="No chirality assignment is reported; no isotopic stereochemical claim is invented.")
d6_product["stereochemistry"]["geometry_status"] = F(None, "specified", page8_evs, method="inference", derivation="The product name assigns E geometry.")
d6_product["stereochemistry"]["representation_status"] = F(None, "not_encoded", page8_evs, method="inference", derivation="Exact isotopic structure is retained by name but no candidate isotopic SMILES is asserted.")
d6_product["identity_note"] = F("The product is reported without deuterium scrambling; the measurement method is not stated.", "The product is reported without deuterium scrambling; the measurement method is not stated.", page8_evs)
modified_molecules[mid("3aa_d6")] = d6_product


# Standalone incoming validation requires each referenced generic scheme and
# molecule.  Include the existing generic coupling record as a dependency only.
generic = get_reaction("general_allenic_CH_diazo_coupling")
updated_reactions.insert(0, generic)

all_used_mids = {
    p["molecule_id"]
    for r in updated_reactions
    for p in r["participants"]
}
all_used_mids.update(modified_molecules)
for molecule_id in sorted(all_used_mids):
    source_molecule = modified_molecules[molecule_id] if molecule_id in modified_molecules else old_molecules[molecule_id]
    bundle["molecules"].append(copy.deepcopy(source_molecule))
bundle["reactions"] = updated_reactions

bundle["source_context"] = [
    {
        "context_id": f"{DOC}:page7_characterization",
        "document_id": DOC,
        "evidence_ids": evidence(7),
        "kind": "characterization",
        "note": "Page 7 contains twelve product characterization blocks (3fa, 3ga, 3gi, 3ha, 3ia, 3ib, 3ie, 3ii, 3ja, 3ik, 3im, 3in). NMR/HRMS/IR/Rf data support identity and isolated amounts but are not separate synthesis experiments."
    },
    {
        "context_id": f"{DOC}:page8_3prime_characterization",
        "document_id": DOC,
        "evidence_ids": evidence(8),
        "kind": "characterization",
        "note": "The 3′ NMR/HRMS/IR/Rf block is characterization only; it was not converted into a new synthesis experiment because Table 1 already contains the distinct conditions that report 3′."
    },
    {
        "context_id": f"{DOC}:page8_mechanistic_interpretation",
        "document_id": DOC,
        "evidence_ids": evidence(8),
        "kind": "author_interpretation",
        "note": "The authors state that retained D6 labeling rules out an alternative methyl C−H cleavage/rearrangement pathway. This interpretation is not an isolated reaction outcome beyond the two mechanistic-study records."
    },
    {
        "context_id": f"{DOC}:page8_si_notice_unassociated",
        "document_id": DOC,
        "evidence_ids": evidence(8),
        "kind": "source_availability",
        "note": "Page 8 states that Supporting Information is available, but source_mapping.json marks SI association unestablished and no SI file was associated."
    },
    {
        "context_id": f"{DOC}:page9_references",
        "document_id": DOC,
        "evidence_ids": evidence(9),
        "kind": "references",
        "note": "Page 9 is bibliography only. Cited literature is context, not a present-paper experimental record."
    },
]

if CANONICAL_SUBMITTED.exists() and not PRIOR_SUBMITTED.exists():
    previous_submission = json.loads(CANONICAL_SUBMITTED.read_text())
    if previous_submission.get("reactions"):
        st.write(PRIOR_SUBMITTED, previous_submission)
st.write(OUT, bundle)
st.write(CANONICAL_SUBMITTED, bundle)

page7_ids = [rid(f"scheme2_scope_{label}") for label in page7_product_data]
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 7,
            "state": "attempted",
            "regions": [
                {"region": "双栏全部产品表征区（3fa、3ga、3gi、3ha、3ia、3ib、3ie、3ii、3ja、3ik、3im、3in）", "finding": "12个既有Scheme 2合成实验的产物名称、分离质量和收率；NMR/HRMS/IR/Rf为表征，不另建反应。", "evidence_ids": evidence(7)},
                {"region": "3ja与3in表征块", "finding": "分数积分与第3页E/Z或Z/E混合物标注一致；比例测定方法未报告，未拆成无依据的单一异构体记录。", "evidence_ids": evidence(3, 7)},
            ],
            "record_ids": page7_ids,
            "unresolved": ["3ik印刷名称与Scheme 2/NMR/HRMS所示含二甲基酰胺结构不一致，未擅自改名或补结构。", "12个复杂产物的完整SMILES尚待专家结构复核；来源E/Z信息已文字保留。"]
        },
        {
            "page": 8,
            "state": "attempted",
            "regions": [
                {"region": "左上3′表征块", "finding": "仅为副产物3′的表征，未说明独立制备条件；不新增合成实验。", "evidence_ids": evidence(8)},
                {"region": "左栏Transformation of [3]Dendralene 3aa with m-CPBA", "finding": "完整衍生化实验：低温加料、0 °C保持、升至室温过夜、后处理、柱层析及4的32.9 mg/50%结果。", "evidence_ids": evidence(3, 8)},
                {"region": "左下Deuterium-Labeling Experiment及3aa-d6表征块", "finding": "含AcOH的详述实验与不含/含AcOH两个产率结果；按eq 3与eq 4保留为两个机理研究记录，并记录条件冲突。", "evidence_ids": evidence(3, 8)},
                {"region": "Associated Content", "finding": "仅提示SI可获得；当前来源映射未建立SI关联，因此不导入SI事实。", "evidence_ids": evidence(8)},
            ],
            "record_ids": [rid("eq2_epoxidation_3aa_to_4"), rid("eq3_deuterium_no_acoh"), rid("eq4_deuterium_with_acoh")],
            "unresolved": ["eq 4的Pd载量在第3页标准条件与第8页具体用量之间冲突。", "第8页DMBQ的0.22 mmol、1.1 equiv与1a-d6 0.3 mmol不自洽。", "无AcOH实验的完整操作未单独叙述；氘保留的测定方法未报告。"]
        },
        {
            "page": 9,
            "state": "attempted",
            "regions": [
                {"region": "整页References (6)–(26)", "finding": "仅为引用文献列表，无本篇新增实验箭头、表格行或操作步骤。", "evidence_ids": evidence(9)}
            ],
            "record_ids": [],
            "unresolved": []
        },
    ]
}
st.write(COVERAGE, coverage)

inventory = {
    "document_id": DOC,
    "inspected_original_images": [str(RUN / f"assets/{DOC}/page{p}.jpg") for p in (7, 8, 9)],
    "page_inventory": {
        "7": {"reaction_arrows": 0, "table_rows": 0, "product_characterization_blocks": 12, "records_updated": page7_ids},
        "8": {"reaction_arrows": 0, "complete_procedures": 2, "paired_mechanistic_results": 2, "characterization_only_blocks": ["3′"], "records_updated": coverage["pages"][1]["record_ids"]},
        "9": {"reaction_arrows": 0, "table_rows": 0, "content": "references only", "records_updated": []},
    },
    "deduplication": "Page-7 product blocks enrich existing Scheme-2 records; page-8 eq 2--4 enrich existing records; no duplicate experiment IDs were created.",
    "limitations": ["Image review is model-assisted and all records remain needs_review.", "Spectral peak lists are source context, not synthesis records."]
}
st.write(RAW / "pages_0007_0009_inventory.json", inventory)

mentions = st.source_mentions(dataset, list(page7_product_data) + ["3′", "m-CPBA", "3aa-d6"], context=180)
st.write(RAW / "pages_0007_0009_mentions.json", mentions)

stereo = {
    "product_4": st.stereo_inspect(product4_smiles),
    "limitations": "standard stereo provides graph diagnostics only; it does not source-verify the candidate or assign the unspecified oxirane center."
}
st.write(RAW / "stereo_diagnostics_pages_0007_0009.json", stereo)

print(json.dumps({
    "bundle": str(OUT),
    "molecules": len(bundle["molecules"]),
    "reactions": len(bundle["reactions"]),
    "experimental_records_updated": len(updated_reactions) - 1,
    "generic_dependency_records": 1,
    "coverage": str(COVERAGE),
}, ensure_ascii=False, indent=2))
