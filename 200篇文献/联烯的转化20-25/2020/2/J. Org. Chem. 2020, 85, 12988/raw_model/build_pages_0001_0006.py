#!/usr/bin/env python3
"""Build the page-1--6 v3.1 extraction bundle from the inspected article images."""
from __future__ import annotations

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DOC = "doc_2a371ced54a752e8"


def ev(page: int, kind: str = "image") -> str:
    return f"{DOC}:page{page}:{kind}"


def src(raw, normalized=None, evidence=()):
    return st.fact(raw, raw if normalized is None else normalized, evidence)


def inf(normalized, evidence, why):
    return st.fact(None, normalized, evidence, method="inference", derivation=why)


def derived(normalized, evidence, why):
    return st.fact(None, normalized, evidence, method="structure_derivation", derivation=why)


def missing(evidence=(), status="unreported"):
    return st.fact(evidence=evidence, status=status)


dataset_path = RUN / "dataset.json"
dataset = st.read(dataset_path)
dataset["documents"][0]["title"] = (
    "Bi(cyclopentyl)diol-Derived Boronates in Highly Enantioselective "
    "Chiral Phosphoric Acid-Catalyzed Allylation, Propargylation, and Crotylation of Aldehydes"
)
st.write(dataset_path, dataset)
bundle = st.empty(dataset["documents"])

molecules = {}
stereo_diagnostics = {}


def add_molecule(mid, label, evidence, *, name=None, smiles=None, chirality=None,
                 geometry=None, representation=None, assignment=None, axial=None,
                 note=None):
    if mid in molecules:
        existing = molecules[mid]
        old = {o["document_id"]: o for o in existing["occurrences"]}
        if DOC not in old:
            existing["occurrences"].append({
                "document_id": DOC,
                "compound_label": src(label, label, evidence),
                "evidence_ids": list(evidence),
            })
        return existing
    m = st.new_molecule(mid, DOC, evidence)
    m["occurrences"][0]["compound_label"] = src(label, label, evidence)
    m["names"]["original"] = src(name or label, name or label, evidence)
    m["names"]["english"] = src(name or label, name or label, evidence)
    if label != (name or label):
        m["names"]["abbreviations"] = src(label, [label], evidence)
    if smiles:
        m["identity_status"] = "source_identified"
        m["representation_type"] = "smiles"
        m["smiles"] = derived(smiles, evidence, "SMILES transcribed from the displayed structure; source identity still needs expert review.")
        diag = st.stereo_inspect(smiles)
        stereo_diagnostics[mid] = diag
        chirality = chirality or ("specified" if diag["atom_centers"] else "achiral")
        if geometry is None:
            alkene = [p for p in diag["potential_stereo"] if str(p["type"]).startswith("Bond_")]
            geometry = "specified" if alkene and all(p["specified"] == "Specified" for p in alkene) else "unspecified" if alkene else "not_applicable"
        representation = representation or "encoded"
        m["stereochemistry"]["isomeric_smiles"] = derived(diag["canonical_isomeric_smiles"], evidence, "RDKit canonicalization of the candidate isomeric SMILES.")
        m["stereochemistry"]["identity_key"] = derived(diag["stereo_graph_key"], evidence, "SHA-256 key of the canonical isomeric graph.")
        if diag["atom_centers"]:
            m["stereochemistry"]["atom_centers"] = derived(diag["atom_centers"], evidence, "RDKit diagnostics of the source-transcribed candidate structure.")
        db = [d for d in diag["double_bonds"] if d["stereo"] != "STEREONONE"]
        if db:
            m["stereochemistry"]["double_bonds"] = derived(db, evidence, "Double-bond geometry encoded in the candidate isomeric SMILES.")
    else:
        m["identity_status"] = "source_identified"
        m["representation_type"] = "unresolved"
        chirality = chirality or "unresolved"
        geometry = geometry or "unresolved"
        representation = representation or "not_encoded"
    reason = "Stereochemical status assigned from the inspected drawing/name; not certified by schema validation."
    m["stereochemistry"]["chirality_status"] = inf(chirality, evidence, reason)
    m["stereochemistry"]["geometry_status"] = inf(geometry, evidence, reason)
    m["stereochemistry"]["representation_status"] = inf(representation, evidence, reason)
    if assignment is not None:
        m["stereochemistry"]["assignment"] = src(assignment, assignment, evidence)
    if axial is not None:
        m["stereochemistry"]["axial"] = src(axial, axial, evidence)
    if note:
        m["identity_note"] = src(note, note, evidence)
    molecules[mid] = m
    return m


# Shared small molecules and materials.
add_molecule("m_aldehyde_1", "1", [ev(3), ev(3, "text")], name="benzaldehyde", smiles="O=Cc1ccccc1")
aldehydes = {
    2: ("benzaldehyde", "O=Cc1ccccc1"),
    3: ("4-bromobenzaldehyde", "O=Cc1ccc(Br)cc1"),
    4: ("4-nitrobenzaldehyde", "O=Cc1ccc([N+](=O)[O-])cc1"),
    5: ("4-methoxybenzaldehyde", "COc1ccc(C=O)cc1"),
    6: ("3-methoxybenzaldehyde", "COc1cccc(C=O)c1"),
    7: ("o-tolualdehyde", "Cc1ccccc1C=O"),
    8: ("methyl 4-formylbenzoate", "COC(=O)c1ccc(C=O)cc1"),
    9: ("piperonal", "O=Cc1ccc2c(c1)OCO2"),
    10: ("thiophene-2-carbaldehyde", "O=Cc1cccs1"),
    11: ("(E)-cinnamaldehyde", "O=C/C=C/c1ccccc1"),
    12: ("phenylpropiolaldehyde", "O=CC#Cc1ccccc1"),
    13: ("3-phenylpropanal", "O=CCCc1ccccc1"),
    14: ("phenylacetaldehyde", "O=CCc1ccccc1"),
    15: ("benzyloxyacetaldehyde", "O=CCOCc1ccccc1"),
    16: ("cyclohexanecarboxaldehyde", "O=CC1CCCCC1"),
    18: ("4-chlorobenzaldehyde", "O=Cc1ccc(Cl)cc1"),
    20: ("4-formylbenzonitrile", "N#Cc1ccc(C=O)cc1"),
    23: ("2-naphthaldehyde", "O=Cc1ccc2ccccc2c1"),
    24: ("1-naphthaldehyde", "O=Cc1cccc2ccccc12"),
    25: ("anthracene-9-carbaldehyde", "O=Cc1c2ccccc2cc2ccccc12"),
}
for number, (name, smi) in aldehydes.items():
    mid = "m_aldehyde_1" if number == 2 else f"m_aldehyde_for_{number}"
    add_molecule(mid, name, [ev(4) if number <= 12 or number >= 17 else ev(5)], name=name, smiles=smi,
                 geometry="specified" if number == 11 else None)

products = {
 2:("(R)-1-phenylbut-3-en-1-ol","C=CC[C@@H](O)c1ccccc1"),
 3:("(R)-1-(4-bromophenyl)but-3-en-1-ol","C=CC[C@@H](O)c1ccc(Br)cc1"),
 4:("(R)-1-(4-nitrophenyl)but-3-en-1-ol","C=CC[C@@H](O)c1ccc([N+](=O)[O-])cc1"),
 5:("(R)-1-(4-methoxyphenyl)but-3-en-1-ol","C=CC[C@@H](O)c1ccc(OC)cc1"),
 6:("(R)-1-(3-methoxyphenyl)but-3-en-1-ol","C=CC[C@@H](O)c1cccc(OC)c1"),
 7:("(R)-1-o-tolylbut-3-en-1-ol","C=CC[C@@H](O)c1ccccc1C"),
 8:("(R)-methyl 4-(1-hydroxybut-3-enyl)benzoate","C=CC[C@@H](O)c1ccc(C(=O)OC)cc1"),
 9:("(R)-1-(benzo[d][1,3]dioxol-5-yl)but-3-en-1-ol","C=CC[C@@H](O)c1ccc2c(c1)OCO2"),
 10:("(R)-1-(thiophen-2-yl)but-3-en-1-ol","C=CC[C@@H](O)c1cccs1"),
 11:("(R,E)-1-phenylhexa-1,5-dien-3-ol","C=CC[C@@H](O)/C=C/c1ccccc1"),
 12:("(R)-1-phenylhex-5-en-1-yn-3-ol","C=CC[C@@H](O)C#Cc1ccccc1"),
 13:("(R)-1-phenylhex-5-en-3-ol","C=CC[C@H](O)CCc1ccccc1"),
 14:("(R)-1-phenylpent-4-en-2-ol","C=CC[C@@H](O)Cc1ccccc1"),
 15:("(R)-1-benzyloxypent-4-en-2-ol","C=CC[C@@H](O)COCc1ccccc1"),
 16:("(R)-1-cyclohexylbut-3-en-1-ol","C=CC[C@@H](O)C1CCCCC1"),
 17:("(1R,2R)-2-methyl-1-phenylbut-3-en-1-ol","C=C[C@@H](C)[C@@H](O)c1ccccc1"),
 18:("(1R,2R)-1-(4-chlorophenyl)-2-methylbut-3-en-1-ol","C=C[C@@H](C)[C@@H](O)c1ccc(Cl)cc1"),
 19:("(1R,2R)-1-(4-bromophenyl)-2-methylbut-3-en-1-ol","C=C[C@@H](C)[C@@H](O)c1ccc(Br)cc1"),
 20:("4-((1R,2R)-1-hydroxy-2-methylbut-3-en-1-yl)benzonitrile","C=C[C@@H](C)[C@@H](O)c1ccc(C#N)cc1"),
 21:("(1R,2R)-1-(4-methoxyphenyl)-2-methylbut-3-en-1-ol","C=C[C@@H](C)[C@@H](O)c1ccc(OC)cc1"),
 22:("methyl 4-((1R,2R)-1-hydroxy-2-methylbut-3-en-1-yl)benzoate","C=C[C@@H](C)[C@@H](O)c1ccc(C(=O)OC)cc1"),
 23:("(1R,2R)-2-methyl-1-(naphthalen-2-yl)but-3-en-1-ol","C=C[C@@H](C)[C@@H](O)c1ccc2ccccc2c1"),
 24:("(1R,2R)-2-methyl-1-(naphthalen-1-yl)but-3-en-1-ol","C=C[C@@H](C)[C@@H](O)c1cccc2ccccc12"),
 25:("(1R,2R)-1-(anthracen-9-yl)-2-methylbut-3-en-1-ol","C=C[C@@H](C)[C@@H](O)c1c2ccccc2cc2ccccc12"),
 26:("(1R,2R)-2-methyl-1-(thiophen-2-yl)but-3-en-1-ol","C=C[C@@H](C)[C@@H](O)c1cccs1"),
 27:("(3S,4R,E)-4-methyl-1-phenylhexa-1,5-dien-3-ol","C=C[C@@H](C)[C@H](O)/C=C/c1ccccc1"),
 28:("(3R,4R)-4-methyl-1-phenylhex-5-en-1-yn-3-ol","C=C[C@@H](C)[C@@H](O)C#Cc1ccccc1"),
 29:("(1S,2R)-1-cyclohexyl-2-methylbut-3-en-1-ol","C=C[C@@H](C)[C@@H](O)C1CCCCC1"),
 30:("(1S,2R)-1-cyclohexyl-2-methylbut-3-en-1-yl 3,5-dinitrobenzoate","C=C[C@@H](C)[C@@H](OC(=O)c1cc([N+](=O)[O-])cc([N+](=O)[O-])c1)C1CCCCC1"),
 31:("(R)-2,2-dimethyl-1-phenylbut-3-en-1-ol","C=CC(C)(C)[C@@H](O)c1ccccc1"),
 32:("(1R,2R)-2,6-dimethyl-1-phenyl-2-vinylhept-5-en-1-ol","C=C[C@@](C)(CCC=C(C)C)[C@@H](O)c1ccccc1"),
 33:("(1R,2S)-2,6-dimethyl-1-phenyl-2-vinylhept-5-en-1-ol","C=C[C@](C)(CCC=C(C)C)[C@@H](O)c1ccccc1"),
 34:("(1R,2R)-1-phenyl-2-vinylpentan-1-ol","CCC[C@H](C=C)[C@@H](O)c1ccccc1"),
 35:("(R)-((R)-cyclohex-2-en-1-yl)(phenyl)methanol","O[C@@H](c1ccccc1)[C@H]1C=CCCC1"),
 36:("(R)-1-phenylbut-3-yn-1-ol","C#CC[C@@H](O)c1ccccc1"),
 37:("(R)-1-(4-bromophenyl)but-3-yn-1-ol","C#CC[C@@H](O)c1ccc(Br)cc1"),
 38:("(R)-1-(4-nitrophenyl)but-3-yn-1-ol","C#CC[C@@H](O)c1ccc([N+](=O)[O-])cc1"),
 39:("(R)-1-(4-methoxyphenyl)but-3-yn-1-ol","C#CC[C@@H](O)c1ccc(OC)cc1"),
 40:("(R)-1-(3-methoxyphenyl)but-3-yn-1-ol","C#CC[C@@H](O)c1cccc(OC)c1"),
 41:("(R)-1-o-tolylbut-3-yn-1-ol","C#CC[C@@H](O)c1ccccc1C"),
 42:("(R)-methyl 4-(1-hydroxybut-3-ynyl)benzoate","C#CC[C@@H](O)c1ccc(C(=O)OC)cc1"),
 43:("(R)-1-(benzo[d][1,3]dioxol-5-yl)but-3-yn-1-ol","C#CC[C@@H](O)c1ccc2c(c1)OCO2"),
 44:("(R)-1-(thiophen-2-yl)but-3-yn-1-ol","C#CC[C@@H](O)c1cccs1"),
 45:("(R,E)-1-phenylhex-1-en-5-yn-3-ol","C#CC[C@@H](O)/C=C/c1ccccc1"),
}
for n, (name, smi) in products.items():
    page = 4 if n <= 35 else 6
    evidence = [ev(page), ev(page, "text")]
    add_molecule(f"m_product_{n}", str(n), evidence, name=name, smiles=smi,
                 chirality="specified", geometry="specified" if n in (11, 27, 45) else None,
                 assignment=name.split(")", 1)[0] + ")" if name.startswith("(") else None)

add_molecule("m_B1", "B1", [ev(3)], name="pinacol-derived allylboronate B1", smiles="C=CCB1OC(C)(C)C(C)(C)O1")
add_molecule("m_B2", "B2", [ev(3), ev(4)], name="bi(cyclopentyl)diol-derived allylboronate B2", smiles="C=CCB1OC2(CCCC2)C2(CCCC2)O1")
add_molecule("m_B3", "B3", [ev(4), ev(5)], name="trans-crotyl boronate B3", smiles="C/C=C/CB1OC2(CCCC2)C2(CCCC2)O1", geometry="specified")
for n in range(4, 9):
    add_molecule(f"m_B{n}", f"B{n}", [ev(4), ev(5)], name=f"substituted allylic boronate B{n}",
                 note="The displayed reagent label is retained; complete structure transcription is deferred.")
add_molecule("m_B9", "B9", [ev(5), ev(6)], name="allenyl pinacol boronate B9", smiles="C=C=CB1OC(C)(C)C(C)(C)O1", chirality="achiral", representation="partial",
             note="Terminal-CH2 allene is not assigned axial chirality; toolkit allene stereo diagnostics remain incomplete.")
add_molecule("m_B10", "B10", [ev(5), ev(6)], name="bi(cyclopentyl)diol-derived allenylboronate B10", smiles="C=C=CB1OC2(CCCC2)C2(CCCC2)O1", chirality="achiral", representation="partial",
             note="Terminal-CH2 allene is not assigned axial chirality; toolkit allene stereo diagnostics remain incomplete.")
for mid, label, name, evidence, chir, geom, rep, axis in [
    ("m_PA1", "PA1", "(R)-TRIP phosphoric acid", [ev(2), ev(3), ev(4), ev(5), ev(6)], "specified", "not_applicable", "not_encoded", [{"axis":"BINOL","configuration":"R"}]),
    ("m_PA2", "PA2", "2,4,6-cyclohexylphenyl-substituted BINOL phosphoric acid PA2", [ev(3), ev(5)], "unresolved", "not_applicable", "not_encoded", None),
    ("m_PA3", "PA3", "(R)-H8-TRIP phosphoric acid", [ev(3), ev(4), ev(5)], "specified", "not_applicable", "not_encoded", [{"axis":"BINOL","configuration":"R"}]),
    ("m_MS4A", "4 Å MS", "4 Å molecular sieves", [ev(3), ev(4), ev(5), ev(6)], "not_applicable", "not_applicable", "not_encoded", None),
]:
    add_molecule(mid, label, evidence, name=name, chirality=chir, geometry=geom, representation=rep, axial=axis)
add_molecule("m_toluene", "toluene", [ev(3), ev(4), ev(5), ev(6)], name="toluene", smiles="Cc1ccccc1")
add_molecule("m_cyclohexane", "cyclohexane", [ev(3), ev(4), ev(5), ev(6)], name="cyclohexane", smiles="C1CCCCC1")

# Generic/source-context entities.
for mid, label, name, page in [
    ("m_generic_aldehyde", "RCHO", "generic aldehyde", 2),
    ("m_generic_allylboronate", "allyl/allenyl boronate", "generic allyl or allenyl boronate", 2),
    ("m_generic_alcohol", "homoallylic/propargylic alcohol", "generic chiral alcohol product", 2),
    ("m_valine_46", "46", "N-Boc-protected L-valine 46", 6),
    ("m_product_47", "47", "nannocystin intermediate 47", 6),
    ("m_product_49", "49", "acrylate intermediate 49", 6),
    ("m_product_50", "50", "pyranone 50", 6),
    ("m_acryloyl_chloride", "acryloyl chloride", "acryloyl chloride", 6),
    ("m_grubbs_II", "Grubbs second-generation catalyst", "Grubbs second-generation ruthenium catalyst", 6),
]:
    add_molecule(mid, label, [ev(page), ev(page, "text")], name=name)


def participant(pid, mid, role, evidence, *, amount=None, equiv=None, molpct=None,
                mass=None, volume=None, ref=None, stage="reaction"):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    if amount is not None:
        raw, val = amount
        p["charges"]["amount"] = src(raw, val, evidence)
    if equiv is not None:
        raw, val = equiv
        p["charges"]["equivalents"] = src(raw, val, evidence)
    if molpct is not None:
        raw, val = molpct
        p["charges"]["mol_percent"] = src(raw, val, evidence)
    if mass is not None:
        raw, val = mass
        p["charges"]["mass"] = src(raw, val, evidence)
    if volume is not None:
        raw, val = volume
        p["charges"]["volume"] = src(raw, val, evidence)
    if ref:
        p["equivalence_reference"] = src(ref, ref, evidence)
    return p


def step(step_id, evidence, pids, operation, temp=None, time=None, stage="reaction"):
    s = st.new_step(step_id, evidence, stage=stage)
    s["participant_ids"] = list(pids)
    s["operation"] = src(operation, operation, evidence)
    if temp is not None:
        s["temperature"] = src(temp[0], temp[1], evidence)
    if time is not None:
        s["time"] = src(time[0], time[1], evidence)
    return s


def outcome(oid, kind, target, raw, normalized, evidence, *, basis=None, yield_type=None, scope=None):
    o = st.new_outcome(oid, kind, [target], evidence)
    o["value"] = src(raw, normalized, evidence)
    if basis is not None:
        o["measurement_basis"] = src(basis, basis, evidence)
    if yield_type is not None:
        o["yield_type"] = src(yield_type, yield_type, evidence)
    if scope is not None:
        o["scope"] = src(scope, scope, evidence)
    return o


def classify_addition(r, evidence):
    changes = [
        {"change_id":"b1","atoms":[{"atom_map":None,"element":"C"},{"atom_map":None,"element":"C"}],"before_order":0,"after_order":1,"change_type":"formed"},
        {"change_id":"b2","atoms":[{"atom_map":None,"element":"C"},{"atom_map":None,"element":"O"}],"before_order":2,"after_order":1,"change_type":"order_decreased"},
        {"change_id":"b3","atoms":[{"atom_map":None,"element":"O"},{"atom_map":None,"element":"H"}],"before_order":0,"after_order":1,"change_type":"hydrogen_gain","hydrogen_count":1},
    ]
    r["classification"]["bond_changes"] = [derived(x, evidence, "Inferred net transformation from source reactant/product drawings; atom mapping is not asserted.") for x in changes]
    labels = [
        {"category_number":2,"code":"C-C","supporting_change_ids":["b1"]},
        {"category_number":24,"code":"O-H","supporting_change_ids":["b3"]},
    ]
    r["classification"]["labels"] = derived(labels, evidence, "Taxonomy labels computed from the explicitly listed inferred bond changes.")
    r["classification"]["coverage"] = inf("partial", evidence, "No atom-mapped reaction SMILES; proton source and omitted hydrolysis partners remain unresolved.")
    r["classification"]["limitations"] = ["Inferred net changes only; no complete atom mapping.", "Hydrolysis/proton-transfer partners are not specified on pages 1–6."]


def reaction(rid, evidence, record_type, categories, location, transformation, unresolved=()):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = src(", ".join(categories), categories, evidence)
    r["location"] = src(location, {"pages": sorted({int(e.split(":page")[1].split(":")[0]) for e in evidence})}, evidence)
    r["transformation_type"] = src(transformation, transformation, evidence)
    r["review"]["unresolved"] = list(unresolved) or ["Machine-extracted record requires expert source review."]
    return r


reactions = []

# Page 1/2 overview and background: generic/current work kept separate from performed experiments.
for rid, page, rtype, title in [
    ("scheme1_generic_aldehyde_allylation", 2, "generic_scheme", "Background generic aldehyde allylation/propargylation"),
    ("scheme2a_prior_pinacol_cpa", 2, "literature_example", "Prior CPA-catalyzed pinacol-boronate additions"),
    ("scheme2b_current_generic", 2, "generic_scheme", "Current bi(cyclopentyl)diol-boronate method overview"),
]:
    evidence = [ev(page), ev(page, "text")]
    r = reaction(rid, evidence, rtype, ["other"], "Scheme 1" if "scheme1" in rid else "Scheme 2", "carbonyl allylation/propargylation")
    r["reaction_name"] = src(title, title, evidence)
    if rid == "scheme1_generic_aldehyde_allylation":
        article_title = dataset["documents"][0]["title"]
        r["notes"].append(src(f"Article title: {article_title}", f"Article title: {article_title}", [ev(1), ev(1,"text")]))
    r["participants"] = [
        participant("generic_aldehyde", "m_generic_aldehyde", "reactant", evidence),
        participant("generic_boronate", "m_generic_allylboronate", "reactant", evidence),
        participant("generic_product", "m_generic_alcohol", "product", evidence),
    ]
    classify_addition(r, evidence)
    reactions.append(r)


def aldehyde_mid(product_no):
    mapping = {
        2:"m_aldehyde_1",17:"m_aldehyde_1",31:"m_aldehyde_1",32:"m_aldehyde_1",33:"m_aldehyde_1",34:"m_aldehyde_1",35:"m_aldehyde_1",36:"m_aldehyde_1",
        3:"m_aldehyde_for_3",19:"m_aldehyde_for_3",37:"m_aldehyde_for_3",
        4:"m_aldehyde_for_4",38:"m_aldehyde_for_4",
        5:"m_aldehyde_for_5",21:"m_aldehyde_for_5",39:"m_aldehyde_for_5",
        6:"m_aldehyde_for_6",40:"m_aldehyde_for_6",
        7:"m_aldehyde_for_7",41:"m_aldehyde_for_7",
        8:"m_aldehyde_for_8",22:"m_aldehyde_for_8",42:"m_aldehyde_for_8",
        9:"m_aldehyde_for_9",43:"m_aldehyde_for_9",
        10:"m_aldehyde_for_10",26:"m_aldehyde_for_10",44:"m_aldehyde_for_10",
        11:"m_aldehyde_for_11",27:"m_aldehyde_for_11",45:"m_aldehyde_for_11",
        12:"m_aldehyde_for_12",28:"m_aldehyde_for_12",
        13:"m_aldehyde_for_13",14:"m_aldehyde_for_14",15:"m_aldehyde_for_15",
        16:"m_aldehyde_for_16",29:"m_aldehyde_for_16",
        18:"m_aldehyde_for_18",20:"m_aldehyde_for_20",23:"m_aldehyde_for_23",24:"m_aldehyde_for_24",25:"m_aldehyde_for_25",
    }
    return mapping[product_no]


def add_catalytic_record(rid, evidence, rtype, categories, location, product_no, boronate,
                         catalyst, catalyst_pct, solvent, temp_raw, temp_norm, time_raw, time_norm,
                         *, yield_pct=None, ee=None, dr=None, conversion=None, isolated=False,
                         scale_mmol=0.10, boronate_mmol=0.12, note=None):
    r = reaction(rid, evidence, rtype, categories, location,
                 "enantioselective aldehyde allylation" if product_no <= 35 else "enantioselective aldehyde propargylation",
                 ["Yield/selectivity measurement methods remain unreported unless explicitly supplied.",
                  "Bond classification is partial and inferred from source drawings."])
    am = aldehyde_mid(product_no)
    pm = f"m_product_{product_no}"
    pids = ["aldehyde", "boronate", "catalyst", "ms", "solvent", "product"]
    r["participants"] = [
        participant("aldehyde", am, "reactant", evidence, amount=(f"{scale_mmol:.2f} mmol", {"number":scale_mmol,"unit":"mmol"})),
        participant("boronate", f"m_{boronate}", "reactant", evidence, amount=(f"{boronate_mmol:.2f} mmol", {"number":boronate_mmol,"unit":"mmol"}),
                    equiv=("1.2 equiv", {"number":1.2,"unit":"equiv"}), ref=am),
        participant("catalyst", f"m_{catalyst}", "catalyst", evidence,
                    molpct=(f"{catalyst_pct:g} mol%", {"number":catalyst_pct,"unit":"mol%"}), ref=am),
        participant("ms", "m_MS4A", "additive", evidence, mass=("50 mg", {"number":50,"unit":"mg"})),
        participant("solvent", f"m_{solvent}", "solvent", evidence, volume=("1 mL", {"number":1,"unit":"mL"})),
        participant("product", pm, "product", evidence),
    ]
    r["reaction_scale"] = src(f"{scale_mmol:.2f} mmol aldehyde", {"number":scale_mmol,"unit":"mmol","reference_molecule_id":am}, evidence)
    r["solvent_conditions"]["composition"] = src(solvent, solvent, evidence)
    r["solvent_conditions"]["total_volume"] = src("1 mL", {"number":1,"unit":"mL"}, evidence)
    r["steps"] = [step("reaction", evidence, pids[:-1], "combine reagents and stir", (temp_raw,temp_norm), (time_raw,time_norm))]
    if yield_pct is not None:
        r["outcomes"].append(outcome("yield", "yield", pm, f"{yield_pct}%", {"number":yield_pct,"unit":"%"}, evidence,
                                            basis="isolated" if isolated else None, yield_type="product yield"))
    if conversion is not None:
        r["outcomes"].append(outcome("conversion", "conversion", am, f"{conversion}%", {"number":conversion,"unit":"%"}, evidence))
    if ee is not None:
        raw = f"{ee}% ee" if isinstance(ee,(int,float)) else f"{ee} ee"
        val = {"number":ee,"unit":"%"} if isinstance(ee,(int,float)) else {"number":float(str(ee).lstrip(">")),"unit":"%","qualifier":"gt","raw_text":raw}
        r["outcomes"].append(outcome("ee", "selectivity", pm, raw, {"metric":"ee","value":val,"configuration":"source product configuration retained"}, evidence,
                                            basis="chiral HPLC and comparison of optical rotation" if location.startswith("Scheme 3") or location.startswith("Scheme 4") else None))
    if dr is not None:
        r["outcomes"].append(outcome("dr", "selectivity", pm, f"dr {dr}", {"metric":"dr","ratio":dr,"order":"major:minor","major":"depicted major stereoisomer"}, evidence,
                                            basis="1H NMR analysis of crude product" if product_no == 17 else None))
    if note:
        r["notes"].append(src(note, note, evidence))
    classify_addition(r, evidence)
    reactions.append(r)
    return r


# Table 1: every printed row plus the separately described acidified-PA1 rescreen.
t1 = [
 (1,"B1","PA1",5,"toluene","room temperature","room temperature","6 h",{"number":6,"unit":"h"},93,"literature_example"),
 (2,"B1","PA2",5,"toluene","room temperature","room temperature","6 h",{"number":6,"unit":"h"},93,"specific_experiment"),
 (3,"B1","PA2",5,"cyclohexane","room temperature","room temperature","6 h",{"number":6,"unit":"h"},94,"specific_experiment"),
 (4,"B2","PA1",5,"cyclohexane","room temperature","room temperature","6 h",{"number":6,"unit":"h"},99,"specific_experiment"),
 (5,"B2","PA2",5,"cyclohexane","room temperature","room temperature","6 h",{"number":6,"unit":"h"},94,"specific_experiment"),
 (6,"B2","PA3",5,"cyclohexane","room temperature","room temperature","6 h",{"number":6,"unit":"h"},77,"specific_experiment"),
 (7,"B2","PA1",2,"cyclohexane","room temperature","room temperature","15 min",{"number":15,"unit":"min"},">99","specific_experiment"),
 (8,"B2","PA1",2,"toluene","−30 °C",{"number":-30,"unit":"°C"},"6 h",{"number":6,"unit":"h"},">99","specific_experiment"),
]
for entry,b,pa,pct,sol,traw,tval,timeraw,timeval,eev,rtype in t1:
    evidence=[ev(3),ev(3,"text")]
    if entry in (7,8): evidence += [ev(4),ev(4,"text")]
    add_catalytic_record(f"table1_entry_{entry}", evidence, rtype, ["optimization"] if rtype=="specific_experiment" else ["other"],
                         f"Table 1 entry {entry}" + ("; Scheme 3 product 2, method A" if entry==7 else ""),
                         2,b,pa,pct,sol,traw,tval,timeraw,timeval,conversion=100,ee=eev,
                         note="Entry 1 is explicitly identified as an earlier report." if entry==1 else None)
add_catalytic_record("table1_acidified_pa1_rescreen", [ev(3),ev(3,"text")], "specific_experiment", ["optimization"],
                     "Page 3 prose rescreen", 2,"B1","PA1",5,"toluene","room temperature","room temperature","6 h",{"number":6,"unit":"h"},
                     conversion=100,ee=94,note="Current acidified-PA1 rescreen is distinct from the prior 93% ee Table 1 entry 1.")

# Scheme 3 scope. Product 2 method A is Table 1 entry 7 and is not duplicated.
scope = {
 2:(94,">99",None,"B2","PA1",0.5), 3:(96,">99",None,"B2","PA1",2), 4:(96,">99",None,"B2","PA1",2),
 5:(97,99,None,"B2","PA1",2), 6:(98,99,None,"B2","PA1",2), 7:(94,99,None,"B2","PA1",2), 8:(93,97,None,"B2","PA1",2),
 9:(97,99,None,"B2","PA1",2), 10:(91,97,None,"B2","PA1",2), 11:(94,99,None,"B2","PA1",2), 12:(92,99,None,"B2","PA1",2),
 13:(95,95,None,"B2","PA3",2), 14:(93,95,None,"B2","PA3",2), 15:(96,93,None,"B2","PA3",2), 16:(92,87,None,"B2","PA3",2),
 17:(97,95,"16:1","B3","PA1",2),18:(98,99,">20:1","B3","PA1",2),19:(93,99,">20:1","B3","PA1",2),20:(97,99,"17:1","B3","PA1",2),
 21:(91,99,"12:1","B3","PA1",2),22:(94,99,"15:1","B3","PA1",2),23:(94,99,"16:1","B3","PA1",2),24:(93,97,"14:1","B3","PA1",2),
 25:(92,93,">20:1","B3","PA1",2),26:(93,99,"17:1","B3","PA1",2),27:(92,95,"9:1","B3","PA1",2),28:(91,99,">20:1","B3","PA1",2),29:(90,91,"13:1","B3","PA1",2),
 31:(97,99,None,"B4","PA1",2),32:(92,77,"10:1","B5","PA1",2),33:(87,95,"11:1","B6","PA1",2),34:(97,91,"13:1","B7","PA1",2),35:(94,23,">20:1","B8","PA1",2),
}
for n,(yp,eev,drv,b,pa,pct) in scope.items():
    if n == 2:
        method="B"; sol="toluene"; traw="−30 °C"; tval={"number":-30,"unit":"°C"}; timeraw="6 h"; timeval={"number":6,"unit":"h"}
    elif n <= 12:
        method="A"; sol="cyclohexane"; traw="room temperature"; tval="room temperature"; timeraw="15 min"; timeval={"number":15,"unit":"min"}
    else:
        method="B"; sol="toluene"; traw="−30 °C"; tval={"number":-30,"unit":"°C"}; timeraw="6 h"; timeval={"number":6,"unit":"h"}
    evidence=[ev(4),ev(4,"text")]
    if n>=13: evidence += [ev(5),ev(5,"text")]
    add_catalytic_record(f"scheme3_product_{n}_method_{method.lower()}", evidence, "specific_experiment", ["substrate_scope"],
                         f"Scheme 3 product {n}, method {method}", n,b,pa,pct,sol,traw,tval,timeraw,timeval,
                         yield_pct=yp,ee=eev,dr=drv,isolated=True)

# Scheme 3 arrow 29 -> 30: conditions and yield are not present on pages 1-6.
evidence=[ev(4),ev(5),ev(5,"text")]
r=reaction("scheme3_29_to_30_derivatization",evidence,"specific_experiment",["derivatization"],"Scheme 3, 29 → 30","ester derivatization",
           ["Reagents, temperature, time and yield are not reported in the inspected pages.","The product is used to establish absolute stereochemistry by X-ray analysis."])
r["participants"]=[participant("alcohol29","m_product_29","reactant",evidence),participant("ester30","m_product_30","product",evidence)]
r["steps"]=[step("derivatization",evidence,["alcohol29"],"convert 29 to ester derivative 30; conditions unreported")]
r["classification"]["bond_changes"]=[derived({"change_id":"b1","atoms":[{"atom_map":None,"element":"C"},{"atom_map":None,"element":"O"}],"before_order":0,"after_order":1,"change_type":"formed"},evidence,"Inferred ester C-O bond formation from source structures.")]
r["classification"]["labels"]=derived([{"category_number":3,"code":"C-O","supporting_change_ids":["b1"]}],evidence,"Taxonomy label computed from the inferred formed bond.")
r["classification"]["coverage"]=inf("partial",evidence,"Reagent and complete atom mapping are unavailable on pages 1-6.")
r["classification"]["limitations"]=["Derivatization reagent and coproducts unreported in inspected pages."]
reactions.append(r)

# Table 2 optimization.
t2=[
 (1,"B9","PA1",20,"toluene","−20 °C",{"number":-20,"unit":"°C"},91,"literature_example"),
 (2,"B9","PA1",5,"toluene","room temperature","room temperature",77,"specific_experiment"),
 (3,"B9","PA1",5,"cyclohexane","room temperature","room temperature",88,"specific_experiment"),
 (4,"B9","PA1",5,"cyclohexane","10 °C",{"number":10,"unit":"°C"},92,"literature_example"),
 (5,"B9","PA2",5,"cyclohexane","room temperature","room temperature",87,"specific_experiment"),
 (6,"B10","PA1",5,"cyclohexane","room temperature","room temperature",98,"specific_experiment"),
 (7,"B10","PA2",5,"cyclohexane","room temperature","room temperature",95,"specific_experiment"),
 (8,"B10","PA3",5,"cyclohexane","room temperature","room temperature",74,"specific_experiment"),
]
for entry,b,pa,pct,sol,traw,tval,eev,rtype in t2:
    note=None
    if entry==1: note="Prior reported experiment; page 2 prose also describes a 96 h prior propargylation, so the duration is unresolved rather than silently assigned."
    if entry==4: note="Prose attributes the identical 92% ee condition to Reddy's report."
    add_catalytic_record(f"table2_entry_{entry}",[ev(5),ev(6),ev(6,"text")],rtype,["optimization"] if rtype=="specific_experiment" else ["other"],
                         f"Table 2 entry {entry}",36,b,pa,pct,sol,traw,tval,"8 h",{"number":8,"unit":"h"},
                         conversion=100,ee=eev,note=note)

# Scheme 4 scope.
scheme4={36:(94,97),37:(94,99),38:(93,97),39:(91,95),40:(96,96),41:(94,93),42:(96,97),43:(94,95),44:(95,91),45:(94,91)}
for n,(yp,eev) in scheme4.items():
    add_catalytic_record(f"scheme4_product_{n}",[ev(6),ev(6,"text")],"specific_experiment",["substrate_scope"],f"Scheme 4 product {n}",
                         n,"B10","PA1",5,"cyclohexane","room temperature","room temperature","8 h",{"number":8,"unit":"h"},
                         yield_pct=yp,ee=eev,isolated=True,
                         note="This scope run is retained separately from Table 2 entry 6; the printed ee values differ (97% vs 98%)." if n==36 else None)

# Page-6 explicitly described applications; details absent from pages 1-6 stay missing.
def simple_application(rid, title, reactants, product_mid, operation, outcomes, category, transformation):
    evidence=[ev(6),ev(6,"text")]
    r=reaction(rid,evidence,"specific_experiment",[category],"Page 6 application discussion",transformation,
               ["Only the facts explicitly stated on page 6 are structured; full conditions are not visible in the inspected page range."])
    r["reaction_name"]=src(title,title,evidence)
    r["participants"]=[participant(pid,mid,role,evidence) for pid,mid,role in reactants]+[participant("product",product_mid,"product",evidence)]
    r["steps"]=[step("application",evidence,[x[0] for x in reactants],operation)]
    for o in outcomes: r["outcomes"].append(o)
    r["classification"]["coverage"]=inf("unresolved",evidence,"Page-6 prose does not provide enough structural detail for a reliable bond-change inventory.")
    r["classification"]["labels"]=missing(evidence)
    reactions.append(r)

simple_application("page6_19_to_47","Mitsunobu preparation of 47",[("alcohol19","m_product_19","reactant"),("valine46","m_valine_46","reactant")],"m_product_47",
                   "Mitsunobu reaction; reagents and conditions not stated on page 6",
                   [outcome("ee","selectivity","m_product_47","99% ee",{"metric":"ee","value":{"number":99,"unit":"%"},"configuration":"unknown"},[ev(6),ev(6,"text")])],"derivatization","Mitsunobu esterification")
simple_application("page6_19_to_49","Acrylate preparation 49",[("alcohol19","m_product_19","reactant"),("acryloyl_chloride","m_acryloyl_chloride","reagent")],"m_product_49",
                   "treat 19 with acryloyl chloride; other conditions unreported",
                   [outcome("ee_retention","selectivity","m_product_49","little ee degradation",{"metric":"ee change","value":"little degradation","configuration":"unknown"},[ev(6),ev(6,"text")])],"derivatization","acrylate formation")
simple_application("page6_49_to_50","Ring-closing metathesis to 50",[("diene49","m_product_49","reactant"),("grubbsII","m_grubbs_II","catalyst")],"m_product_50",
                   "ring-closing metathesis with Grubbs second-generation catalyst",
                   [outcome("ee","selectivity","m_product_50","83% ee",{"metric":"ee","value":{"number":83,"unit":"%"},"configuration":"unknown"},[ev(6),ev(6,"text")]),
                    outcome("dr","selectivity","m_product_50","dr >20:1",{"metric":"dr","ratio":">20:1","order":"major:minor","major":"unknown"},[ev(6),ev(6,"text")])],"derivatization","ring-closing metathesis")

# Gram-scale examples reuse product registry entries and report a 10 mmol benzaldehyde scale.
for kind,prod,b,yp,eev,drv,trans in [
    ("allylation",2,"B2",97,99,None,"enantioselective aldehyde allylation"),
    ("crotylation",17,"B3",97,93,">20:1","enantioselective aldehyde crotylation"),
    ("propargylation",36,"B10",97,98,None,"enantioselective aldehyde propargylation"),
]:
    evidence=[ev(6),ev(6,"text")]
    r=reaction(f"page6_gram_{kind}",evidence,"specific_experiment",["scale_up"],"Page 6 gram-scale discussion",trans,
               ["Detailed temperature, time, reagent charges and yield measurement basis are not stated on page 6."])
    r["participants"]=[participant("aldehyde","m_aldehyde_1","reactant",evidence,amount=("10 mmol",{"number":10,"unit":"mmol"})),
                       participant("boronate",f"m_{b}","reactant",evidence),participant("product",f"m_product_{prod}","product",evidence)]
    r["reaction_scale"]=src("10 mmol benzaldehyde",{"number":10,"unit":"mmol","reference_molecule_id":"m_aldehyde_1"},evidence)
    r["steps"]=[step("scale_up",evidence,["aldehyde","boronate"],f"conduct gram-scale {kind}; detailed conditions unreported")]
    r["outcomes"]=[outcome("yield","yield",f"m_product_{prod}",f"{yp}%",{"number":yp,"unit":"%"},evidence),
                   outcome("ee","selectivity",f"m_product_{prod}",f"{eev}% ee",{"metric":"ee","value":{"number":eev,"unit":"%"},"configuration":"source product configuration retained"},evidence)]
    if drv: r["outcomes"].append(outcome("dr","selectivity",f"m_product_{prod}",f"dr {drv}",{"metric":"dr","ratio":drv,"order":"major:minor","major":"depicted major stereoisomer"},evidence))
    classify_addition(r,evidence)
    reactions.append(r)

bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id":f"{DOC}:page1_title_abstract","document_id":DOC,"kind":"article_metadata","evidence_ids":[ev(1),ev(1,"text")],
     "title":dataset["documents"][0]["title"],"finding":"Title, abstract and graphical overview; not a distinct performed experiment."},
    {"context_id":f"{DOC}:page2_background_natural_products","document_id":DOC,"kind":"background","evidence_ids":[ev(2)],
     "finding":"Natural-product structures and cited historical methods are background, not experiments in this paper."},
    {"context_id":f"{DOC}:page5_prior_b1_comparison","document_id":DOC,"kind":"literature_comparison","evidence_ids":[ev(5),ev(5,"text")],
     "finding":"Earlier B1 benzyloxyacetaldehyde result (79% ee) retained as comparison, not recast as this paper's run."},
    {"context_id":f"{DOC}:page6_catalyst_recovery","document_id":DOC,"kind":"workup_observation","evidence_ids":[ev(6),ev(6,"text")],
     "finding":"Catalyst recovery 84% after chromatography with 5% MeOH/DCM and reacidification with 6 M HCl; source does not assign it to one gram-scale row."},
    {"context_id":f"{DOC}:page6_mechanism_prose","document_id":DOC,"kind":"proposed_mechanism","evidence_ids":[ev(6),ev(6,"text")],
     "finding":"Proposed bifunctional, chair-like six-membered transition state; Scheme 6 itself is on an uninspected later page."},
]

st.write(RUN / "submitted.json", bundle)
st.write(RUN / "raw_model/pages_0001_0006_standard_stereo_checks.json", stereo_diagnostics)
raw = {
    "document_id": DOC,
    "visual_pages_inspected": [1,2,3,4,5,6],
    "record_count": len(reactions),
    "molecule_count": len(molecules),
    "decisions": [
        "Table rows and individual scope arrows are separate records.",
        "Table 1 entry 7 and Scheme 3 product 2 method A are the same experiment and share one record.",
        "Prior work and generic schemes are typed separately from current specific experiments.",
        "Table 2 entry 6 and Scheme 4 product 36 remain distinct because their printed ee values differ.",
        "Scope yields are isolated per the Scheme 3/4 footnotes; table conversions retain an unreported measurement basis.",
        "Molecule structures are candidate transcriptions and remain unreviewed; unresolved boronate/application structures are not guessed.",
    ],
    "records": [{"reaction_id":r["reaction_id"],"location":r["location"]["raw_value"],"record_type":r["record_type"],"evidence_ids":r["evidence_ids"]} for r in reactions],
}
st.write(RUN / "raw_model/pages_0001_0006_raw_extraction.json", raw)

by_page = {p: [] for p in range(1,7)}
for r in reactions:
    for p in by_page:
        if ev(p) in r["evidence_ids"] or ev(p,"text") in r["evidence_ids"]:
            by_page[p].append(r["reaction_id"])
coverage = {
    "document_id":DOC,
    "semantic_accuracy_verified":False,
    "pages":[
        {"page":1,"state":"attempted","regions":[
            {"region":"title, graphical abstract and abstract","finding":"Article metadata and a current-work overview; no distinct condition-specific experiment.","evidence_ids":[ev(1),ev(1,"text")]},
            {"region":"introduction text","finding":"Background and cited work; not recast as this paper's experiments.","evidence_ids":[ev(1),ev(1,"text")]}
        ],"record_ids":by_page[1],"unresolved":["Graphical overview is generic and does not supply row-level conditions."]},
        {"page":2,"state":"attempted","regions":[
            {"region":"Scheme 1","finding":"Background natural products and generic open/closed-transition-state allylation framework.","evidence_ids":[ev(2),ev(2,"text")]},
            {"region":"Scheme 2a","finding":"Prior pinacol-boronate literature method, classified as literature_example.","evidence_ids":[ev(2),ev(2,"text")]},
            {"region":"Scheme 2b","finding":"Generic overview of the present bi(cyclopentyl)diol-boronate method.","evidence_ids":[ev(2),ev(2,"text")]}
        ],"record_ids":by_page[2],"unresolved":["Individual substrate identities and conditions are intentionally absent from the generic schemes."]},
        {"page":3,"state":"attempted","regions":[
            {"region":"Table 1 rows 1-8 and footnotes","finding":"Eight optimization rows transcribed; entry 1 is prior work, entries 2-8 are current optimization runs.","evidence_ids":[ev(3),ev(3,"text")]},
            {"region":"results prose above/beside Table 1","finding":"Separate acidified-PA1 rescreen (94% ee) and catalyst/boronate comparison retained.","evidence_ids":[ev(3),ev(3,"text")]}
        ],"record_ids":by_page[3],"unresolved":["Conversion measurement basis and most ee analytical methods are not stated on this page."]},
        {"page":4,"state":"attempted","regions":[
            {"region":"Scheme 3A-B and method headers","finding":"Products 2-16, yields, ee values, Method A/B conditions and footnotes transcribed.","evidence_ids":[ev(4),ev(4,"text")]},
            {"region":"Scheme 3C-D","finding":"Products 17-35, dr/ee/yields and boronate identities transcribed; 29→30 retained as a partial derivatization.","evidence_ids":[ev(4),ev(4,"text")]}
        ],"record_ids":by_page[4],"unresolved":["B4-B8 full structures and some dr-method assignments require expert review; only product 17 is clearly tied to the crude-1H-NMR footnote."]},
        {"page":5,"state":"attempted","regions":[
            {"region":"Scheme 3 continuation discussion","finding":"Aliphatic/crotyl/substituted-boronate scope and stereochemical discussion linked to Scheme 3 records.","evidence_ids":[ev(5),ev(5,"text")]},
            {"region":"Table 2 rows 1-8 and footnotes","finding":"Eight propargylation optimization rows transcribed; entries 1 and 4 are explicitly prior literature conditions.","evidence_ids":[ev(5),ev(5,"text")]}
        ],"record_ids":by_page[5],"unresolved":["Table 2 entry 1 arrow shows 8 h while page-2 prose mentions a 96 h prior propargylation; identity of the duration is unresolved."]},
        {"page":6,"state":"attempted","regions":[
            {"region":"Table 2 discussion","finding":"Entry interpretations and literature attribution checked against the printed table.","evidence_ids":[ev(6),ev(6,"text")]},
            {"region":"Scheme 4 products 36-45","finding":"Ten substrate-scope arrows with isolated yield and ee transcribed.","evidence_ids":[ev(6),ev(6,"text")]},
            {"region":"synthetic applications and gram-scale prose","finding":"19→47, 19→49, 49→50 and three 10 mmol scale experiments retained as partial records.","evidence_ids":[ev(6),ev(6,"text")]},
            {"region":"mechanism discussion","finding":"Proposed transition-state model retained as source context because Scheme 6 lies outside the visually covered pages.","evidence_ids":[ev(6),ev(6,"text")]}
        ],"record_ids":by_page[6],"unresolved":["Application conditions and structures 47/49/50 are not visible on page 6; the catalyst-recovery statement is not assigned to one scale-up record."]},
    ]
}
st.write(RUN / "coverage/pages_0001_0006.json", coverage)
print(json.dumps({"records":len(reactions),"molecules":len(molecules),"coverage_pages":len(coverage["pages"])},ensure_ascii=False))
