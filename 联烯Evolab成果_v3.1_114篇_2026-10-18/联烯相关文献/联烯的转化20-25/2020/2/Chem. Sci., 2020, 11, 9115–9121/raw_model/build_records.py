#!/usr/bin/env python3
"""Build source-bounded v3.1 records for doc_3fdb45b913df167c, PDF pages 1-6."""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_3fdb45b913df167c"
RUN = ROOT / "runs/current"


def ev(page: int) -> list[str]:
    return [f"{DOC}:page{page}:text", f"{DOC}:page{page}:image"]


def evs(*pages: int) -> list[str]:
    return [item for page in pages for item in ev(page)]


def ef(raw, normalized=None, evidence=()):
    return st.fact(raw, raw if normalized is None else normalized, evidence)


def df(normalized, evidence, derivation):
    return st.fact(None, normalized, evidence, method="structure_derivation", derivation=derivation)


def inf(normalized, evidence, derivation):
    return st.fact(None, normalized, evidence, method="inference", derivation=derivation)


def missing(status="unreported", evidence=()):
    return st.fact(None, None, evidence, status=status)


def q(number, unit, raw=None, **extra):
    value = {"number": number, "unit": unit}
    value.update(extra)
    if raw is not None:
        value["raw_text"] = raw
    return value


dataset = json.loads((RUN / "dataset.json").read_text())
doc = dataset["documents"][0]
doc["title"] = "Catalytic enantioselective allene–anhydride approach to β,γ-unsaturated enones bearing an α-all-carbon-quaternary center"
dataset["documents"][0] = doc
st.write(RUN / "dataset.json", dataset)
bundle = st.empty(dataset["documents"])
molecules: dict[str, dict] = {}


def add_molecule(mid, label, page, *, name=None, smiles=None, chirality="unresolved",
                 geometry="unresolved", representation=None, assignment=None,
                 mixture=None, identity_note=None):
    evidence = ev(page)
    m = st.new_molecule(mid, DOC, evidence)
    m["occurrences"][0]["compound_label"] = ef(label, label, evidence)
    m["names"]["original"] = ef(name or label, name or label, evidence)
    m["names"]["english"] = ef(name, name, evidence) if name else missing()
    m["names"]["abbreviations"] = ef(label, [label], evidence)
    m["names"]["aliases"] = missing()
    m["cas_numbers"] = missing()
    if smiles:
        m["identity_status"] = "source_identified"
        m["representation_type"] = "smiles"
        m["smiles"] = df(smiles, evidence, "Candidate connectivity transcribed from the named/source-drawn molecule; parsing does not establish source correctness.")
        representation = representation or "encoded"
        diag = st.stereo_inspect(smiles)
        m["stereochemistry"]["isomeric_smiles"] = df(diag["canonical_isomeric_smiles"], evidence, "RDKit canonicalization of the candidate source-derived SMILES.")
        m["stereochemistry"]["identity_key"] = df(diag["stereo_graph_key"], evidence, "Stereo graph key calculated from the candidate isomeric SMILES.")
    else:
        m["identity_status"] = "source_identified"
        m["representation_type"] = "unresolved"
        m["smiles"] = missing()
        representation = representation or "unresolved"
    m["stereochemistry"]["chirality_status"] = inf(chirality, evidence, "Status preserves the source label/drawing; absence of a tag is not treated as achirality.")
    m["stereochemistry"]["geometry_status"] = inf(geometry, evidence, "Status preserves source alkene/oxime geometry information and unresolved cases.")
    m["stereochemistry"]["representation_status"] = inf(representation, evidence, "Reports whether the retained structure candidate encodes all source-specified stereochemistry.")
    if assignment is not None:
        m["stereochemistry"]["assignment"] = ef(assignment, assignment, evidence)
    if mixture is not None:
        m["stereochemistry"]["mixture"] = ef(mixture, mixture, evidence)
    note = identity_note or ("Exact structure retained only by source label/drawing; no unsupported complete SMILES was asserted." if not smiles else "Candidate structure requires expert comparison with the source drawing.")
    m["identity_note"] = ef(note, note, evidence)
    molecules[mid] = m
    return m


def add_template(mid, label, page, fixed, variables):
    m = add_molecule(mid, label, page, chirality="unresolved", geometry="unresolved", representation="unresolved")
    evidence = ev(page)
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": ef(fixed, {"source_description": fixed}, evidence),
        "variables": [{
            "label": v,
            "attachment_points": ef(f"source-marked attachment for {v}", [f"attachment:{v}"], evidence),
            "definition": ef(desc, desc, evidence),
            "allowed_values": missing(),
            "constraints": missing(),
        } for v, desc in variables],
        "instance_mapping": [],
        "reason": ef("generic R-group drawing", "Source depicts variable substituents; a concrete molecule is not asserted.", evidence),
    }
    return m


# Generic identities used only by the page-1 overview/literature arrows.
add_template("m_gen_anhydride", "carboxylic anhydride (generic)", 1, "R-C(O)-O-C(O)-R", [("R", "aryl or alkyl acyl substituent")])
add_template("m_gen_organometallic", "R1-M (generic)", 1, "R1-M", [("R1", "organic substituent"), ("M", "metal")])
add_template("m_gen_styrene", "substituted styrene (generic)", 1, "Ar2-CH=CH2", [("Ar2", "aryl group")])
add_template("m_gen_monoallene", "monosubstituted terminal allene (generic)", 1, "R-CH=C=CH2", [("R", "substituent")])
add_template("m_gen_disuballene", "1,1-disubstituted terminal allene 1 (generic)", 1, "R1(Ar)-C=C=CH2", [("R1", "alkyl group"), ("Ar", "aryl group")])
add_template("m_gen_ketone", "ketone product (generic)", 1, "generic ketone product shown in Scheme 1", [("R", "acyl substituent"), ("R1", "substrate substituent")])
add_template("m_gen_current_product", "β,γ-unsaturated enone 3 (generic)", 1, "R-C(O)-C*(Ar)(R1)-CH=CH2", [("R", "aryl or alkyl"), ("Ar", "aryl"), ("R1", "alkyl")])

# Core substrates, agents and solvents. Exact complex ligands remain drawing-linked.
add_molecule("m_1a", "1a", 2, name="3-phenyl-nona-1,2-diene", smiles="C=C=C(c1ccccc1)CCCCCC", chirality="achiral", geometry="not_applicable", representation="partial")
add_molecule("m_2a", "2a", 2, name="benzoic anhydride", smiles="O=C(OC(=O)c1ccccc1)c1ccccc1", chirality="achiral", geometry="not_applicable")
add_molecule("m_cu_oac2", "Cu(OAc)2", 2, name="copper(II) acetate", smiles="[Cu+2].CC(=O)[O-].CC(=O)[O-]", chirality="achiral", geometry="not_applicable")
add_molecule("m_silane", "Me(MeO)2SiH", 2, name="methyldimethoxysilane", smiles="CO[SiH](C)OC", chirality="achiral", geometry="not_applicable")
add_molecule("m_thf", "THF", 2, name="tetrahydrofuran", smiles="C1CCOC1", chirality="achiral", geometry="not_applicable")
add_molecule("m_cyclohexane", "cyclohexane", 2, smiles="C1CCCCC1", chirality="achiral", geometry="not_applicable")
add_molecule("m_dioxane", "dioxane", 2, name="1,4-dioxane", smiles="O1CCOCC1", chirality="achiral", geometry="not_applicable")
add_molecule("m_dcm", "DCM", 2, name="dichloromethane", smiles="ClCCl", chirality="achiral", geometry="not_applicable")
add_molecule("m_toluene", "toluene", 2, smiles="Cc1ccccc1", chirality="achiral", geometry="not_applicable")

ligands = ["DPPM", "DPPE", "DPPP", "DPPF", "BINAP", "PCy3", "(R)-BINAP", "(R)-MeO-BIPHEP", "(R)-C3-Tunephos", "(R)-Garphos", "(S)-Synphos", "(R)-Segphos", "(R,R)-BDPP", "(R)-SDP", "(R,R)-BPE", "(S,S)-BPE"]
for label in ligands:
    mid = "m_lig_" + label.lower().replace("(", "").replace(")", "").replace(",", "").replace("-", "_")
    chiral = "specified" if label.startswith("(") else ("unspecified" if label == "BINAP" else "unresolved")
    add_molecule(mid, label, 2 if label != "(S,S)-BPE" else 3, chirality=chiral, geometry="not_applicable", representation="not_encoded" if chiral == "specified" else "unresolved", assignment=label.split("-")[0] if chiral == "specified" else None)

def lig_id(label):
    return "m_lig_" + label.lower().replace("(", "").replace(")", "").replace(",", "").replace("-", "_")


# Scope substrates/anhydrides and stereochemically distinct products.
for label in [f"1{x}" for x in "bcdefghijk"]:
    add_molecule(f"m_{label}", label, 3, chirality="unspecified", geometry="not_applicable", representation="unresolved", identity_note="Source drawing identifies a terminal 1,1-disubstituted allene; exact connectivity is not promoted without independent structure transcription.")
for label in [f"2{x}" for x in "bcdefghij"]:
    add_molecule(f"m_{label}", label, 3, chirality="achiral", geometry="not_applicable", representation="unresolved", identity_note="Source-drawn symmetric carboxylic anhydride; exact connectivity retained by label/drawing only.")

rac_products = ["3aa", "3ab", "3ac", "3ad", "3ae", "3ia", "3ja", "3ka"]
for label in rac_products:
    add_molecule(f"m_rac_{label}", f"rac-{label}", 2, chirality="racemic", geometry="not_applicable", representation="not_encoded", mixture="racemic product", identity_note="Racemic β,γ-unsaturated enone shown in Table 2; terminal alkene has no E/Z assignment.")

r_products = ["3aa", "3ba", "3ca", "3da", "3ea", "3fa", "3ga", "3ha", "3ia", "3ja", "3ka", "3ab", "3ac", "3ad", "3ae", "3af", "3ag", "3ah", "3ai", "3aj"]
for label in r_products:
    assignment = {"configuration": "R", "basis": "source product label"}
    if label == "3ea":
        assignment = {"configuration": "R", "basis": "single-crystal X-ray crystallography", "CCDC": "1939837"}
    add_molecule(f"m_R_{label}", f"(R)-{label}", 3, chirality="specified", geometry="not_applicable", representation="not_encoded", assignment=assignment, identity_note="Source specifies the R product; exact isomeric SMILES was not asserted from the wedge drawing.")
for label in ["3af", "3ai", "3aj"]:
    add_molecule(f"m_S_{label}", f"(S)-{label}", 3, chirality="specified", geometry="not_applicable", representation="not_encoded", assignment={"configuration": "S", "basis": "source product label"}, identity_note="Source specifies the S enantiomer; exact isomeric SMILES was not asserted from the wedge drawing.")

# Scheme 2 products and reagents.
add_molecule("m_RZ_4", "(R,Z)-4", 4, chirality="specified", geometry="specified", representation="not_encoded", assignment={"atom_center": "R", "oxime_geometry": "Z"})
for label in ["5", "6", "7"]:
    add_molecule(f"m_R_{label}", f"(R)-{label}", 4, chirality="specified", geometry="not_applicable", representation="not_encoded", assignment={"configuration": "R", "basis": "source product label"})
scheme2_agents = {
    "m_hydroxylamine_hcl": "NH2OH·HCl", "m_naoac": "NaOAc", "m_etoh": "EtOH", "m_water": "H2O",
    "m_nitroboronic": "3-O2NC6H4B(OH)2", "m_pd_dppf": "Pd(dppf)Cl2", "m_k2co3": "K2CO3", "m_dmso": "DMSO",
    "m_pd_phcn": "Ph(PhCN)2Cl2", "m_cucl2_h2o": "CuCl2·2H2O", "m_agno2": "AgNO2", "m_tbuoh": "t-BuOH", "m_me_no2": "CH3NO2",
    "m_pdc": "Pd/C", "m_etoac": "EtOAc", "m_o2": "O2", "m_h2": "H2", "m_ar": "Ar",
}
for mid, label in scheme2_agents.items():
    add_molecule(mid, label, 4, chirality="unresolved", geometry="not_applicable", representation="unresolved")

reactions = []


def base_reaction(rid, evidence, record_type, categories, location, name, transform):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = ef(categories, categories, evidence)
    r["location"] = ef(location, location, evidence)
    r["reaction_name"] = ef(name, name, evidence)
    r["transformation_type"] = ef(transform, transform, evidence) if transform else missing()
    r["review"]["unresolved"] = []
    return r


def participant(pid, mid, role, evidence, *, ref=None, amount=None, equiv=None, molpct=None, volume=None):
    p = st.new_participant(pid, mid, role, evidence)
    if amount:
        p["charges"]["amount"] = ef(amount[0], q(amount[1], amount[2], amount[0]), evidence)
    if equiv is not None:
        raw = f"{equiv} equiv."
        p["charges"]["equivalents"] = ef(raw, q(equiv, "equiv", raw), evidence)
    if molpct is not None:
        raw = f"{molpct} mol%"
        p["charges"]["mol_percent"] = ef(raw, q(molpct, "mol%", raw), evidence)
    if volume:
        p["charges"]["volume"] = ef(volume[0], q(volume[1], volume[2], volume[0]), evidence)
    if ref and (equiv is not None or molpct is not None):
        p["equivalence_reference"] = inf(ref, evidence, "Denominator assigned from the explicitly stated limiting reagent/amount relationship.")
    return p


def set_classification(r, kind="hydroacylation"):
    evidence = r["evidence_ids"]
    if kind == "not_applicable":
        r["classification"]["coverage"] = inf("not_applicable", evidence, "No observed chemical transformation is reported.")
        r["classification"]["labels"] = inf([], evidence, "No bond-formation labels are assigned to an unobserved transformation.")
        return
    if kind == "unresolved":
        r["classification"]["coverage"] = inf("unresolved", evidence, "Observed transformation is insufficiently specified for bond assignment.")
        r["classification"]["labels"] = inf([], evidence, "No bond labels asserted.")
        return
    changes = []
    labels = []
    if kind == "hydroacylation":
        changes = [
            {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
            {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
        ]
        labels = [{"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1"]}, {"category_number": 1, "code": "C-H", "supporting_change_ids": ["b2"]}]
    elif kind == "oxime":
        changes = [{"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "N"}], "before_order": 0, "after_order": 2, "change_type": "formed"}]
        labels = [{"category_number": 16, "code": "DB-C-N", "supporting_change_ids": ["b1"]}]
    elif kind == "suzuki":
        changes = [{"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"}]
        labels = [{"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1"]}]
    elif kind == "wacker":
        changes = [{"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}], "before_order": 0, "after_order": 2, "change_type": "formed"}]
        labels = [{"category_number": 14, "code": "DB-C-O", "supporting_change_ids": ["b1"]}]
    elif kind == "hydrogenation":
        changes = [
            {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 2},
            {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
        ]
        labels = [{"category_number": 1, "code": "C-H", "supporting_change_ids": ["b1"]}]
    r["classification"]["bond_changes"] = [df(c, evidence, "Partial net bond change derived from the source-drawn reactant/product structures; atom mapping remains unresolved.") for c in changes]
    r["classification"]["labels"] = inf(labels, evidence, "Taxonomy labels include all and only the asserted formed/promoted bonds in the partial source-derived change list.")
    r["classification"]["coverage"] = inf("partial", evidence, "No complete atom mapping; coproducts and exact hydrogen provenance are unresolved.")
    r["classification"]["limitations"] = ["Partial drawing-derived net changes only; no complete atom mapping.", "Implicit hydrogen gain does not identify the elementary proton source."]


def add_outcome(r, oid, kind, target, raw, normalized, *, basis=None, scope=None, yield_type=None):
    o = st.new_outcome(oid, kind, [target] if target else [], r["evidence_ids"])
    o["value"] = ef(raw, normalized, r["evidence_ids"])
    if basis is not None:
        o["measurement_basis"] = ef(basis, basis, r["evidence_ids"])
    if scope is not None:
        o["scope"] = ef(scope, scope, r["evidence_ids"])
    if yield_type is not None:
        o["yield_type"] = ef(yield_type, yield_type, r["evidence_ids"])
    r["outcomes"].append(o)


def hydro(rid, pages, *, categories, location, allene="m_1a", anhydride="m_2a", product="m_rac_3aa",
          ligand="DPPE", solvent="m_thf", allene_amount=None, allene_equiv=None,
          anhydride_amount=None, anhydride_equiv=None, catalyst=2.5, ligand_pct=2.5,
          silane_equiv=3.0, solvent_volume=None, concentration=None, temperature="r.t.", time_h=6,
          yield_pct=None, yield_basis=None, recovery_pct=None, recovery_basis=None, ee=None,
          ee_raw=None, ee_scope=None, classification="hydroacylation", notes=()):
    evidence = evs(*pages)
    r = base_reaction(rid, evidence, "specific_experiment", categories, location, "Cu-catalyzed hydroacylation of an allene with a carboxylic anhydride", "allene–anhydride hydroacylation")
    ref = anhydride if anhydride_amount and allene_equiv is not None else allene
    r["participants"] = [
        participant("allene", allene, "reactant", evidence, ref=ref, amount=allene_amount, equiv=allene_equiv),
        participant("anhydride", anhydride, "reactant", evidence, ref=ref, amount=anhydride_amount, equiv=anhydride_equiv),
        participant("product", product, "product", evidence),
        participant("copper", "m_cu_oac2", "catalyst", evidence, ref=ref, molpct=catalyst),
        participant("ligand", lig_id(ligand), "ligand", evidence, ref=ref, molpct=ligand_pct),
        participant("silane", "m_silane", "reagent", evidence, ref=ref, equiv=silane_equiv),
        participant("solvent", solvent, "solvent", evidence, volume=solvent_volume),
    ]
    if ref == allene and allene_amount:
        r["reaction_scale"] = ef(allene_amount[0], {"reference_molecule_id": allene, "amount": q(allene_amount[1], allene_amount[2], allene_amount[0])}, evidence)
    elif ref == anhydride and anhydride_amount:
        r["reaction_scale"] = ef(anhydride_amount[0], {"reference_molecule_id": anhydride, "amount": q(anhydride_amount[1], anhydride_amount[2], anhydride_amount[0])}, evidence)
    r["solvent_conditions"]["composition"] = ef(molecules[solvent]["occurrences"][0]["compound_label"]["raw_value"], [solvent], evidence)
    if solvent_volume:
        r["solvent_conditions"]["total_volume"] = ef(solvent_volume[0], q(solvent_volume[1], solvent_volume[2], solvent_volume[0]), evidence)
    if concentration:
        r["solvent_conditions"]["concentration"] = ef(concentration, q(float(concentration.split()[0]), "M", concentration), evidence)
        r["solvent_conditions"]["concentration_basis"] = missing(evidence=evidence)
    step = st.new_step("reaction", evidence)
    step["participant_ids"] = [p["participant_id"] for p in r["participants"]]
    step["operation"] = inf("react under the stated conditions", evidence, "Source supplies conditions but not an addition sequence in the main text.")
    step["temperature"] = ef(temperature, temperature if temperature == "r.t." else q(temperature, "°C", f"{temperature} °C"), evidence)
    step["time"] = ef(f"{time_h} h", q(time_h, "h", f"{time_h} h"), evidence)
    r["steps"] = [step]
    r["procedure"]["original_text"] = ef(location, location, evidence)
    if yield_pct is not None:
        add_outcome(r, "yield", "yield", product, f"{yield_pct}%", q(yield_pct, "%", f"{yield_pct}%"), basis=yield_basis, yield_type="NMR yield" if yield_basis and "NMR" in yield_basis else None)
    if recovery_pct is not None:
        add_outcome(r, "recovery", "recovery", allene, f"{recovery_pct}%", q(recovery_pct, "%", f"{recovery_pct}%"), basis=recovery_basis)
    if ee is not None:
        raw = ee_raw or f"{ee}% ee"
        norm = q(abs(ee), "%", raw)
        if ee < 0:
            norm["sign"] = "negative"
        add_outcome(r, "ee", "selectivity", product, raw, norm, scope=ee_scope or {"metric": "ee", "major_configuration": "source product label"})
    set_classification(r, classification)
    r["notes"] = [ef(n, n, evidence) for n in notes]
    if yield_pct is not None and yield_basis is None:
        r["review"]["unresolved"].append("Yield measurement basis is not reported on the supplied pages.")
    if recovery_pct is not None and recovery_basis is None:
        r["review"]["unresolved"].append("Recovery measurement basis is not explicitly reported.")
    if ee is not None:
        r["review"]["unresolved"].append("ee measurement method is not reported on the supplied pages; Table 3 signs are retained from the figure where visible.")
    if any("entry 12" in n for n in notes):
        r["review"]["unresolved"].append("Table 1 entry 12 retains the printed 1.0 mL common footnote and 0.1 M header although the altered reagent scale makes the concentration basis ambiguous.")
    r["review"]["unresolved"].append("Exact procedure, workup and purification are unavailable because no SI was supplied.")
    reactions.append(r)
    return r


# Page 1: five prior-work arrows and one generic current-work arrow.
literature_specs = [
    ("scheme1_lit_a", "Scheme 1A-a", "acylation of organometallic reagents", ["m_gen_anhydride", "m_gen_organometallic"], "m_gen_ketone"),
    ("scheme1_lit_b", "Scheme 1A-b", "Rh-catalyzed hydrometalation of styrenes and reaction with anhydrides", ["m_gen_anhydride", "m_gen_styrene"], "m_gen_ketone"),
    ("scheme1_lit_c", "Scheme 1A-c", "CuH-catalyzed hydroacylation of styrenes", ["m_gen_anhydride", "m_gen_styrene"], "m_gen_ketone"),
    ("scheme1_lit_d", "Scheme 1A-d", "Pd-catalyzed hydroacylation of monosubstituted terminal allenes", ["m_gen_anhydride", "m_gen_monoallene"], "m_gen_ketone"),
    ("scheme1_lit_e", "Scheme 1A-e", "Cu-catalyzed boroacylation of 1,1-disubstituted allenes", ["m_gen_anhydride", "m_gen_disuballene"], "m_gen_ketone"),
]
for rid, location, name, reactants, product in literature_specs:
    r = base_reaction(rid, evs(1), "literature_example", ["other"], location, name, name)
    r["participants"] = [participant(f"r{i+1}", mid, "reactant", evs(1)) for i, mid in enumerate(reactants)] + [participant("p", product, "product", evs(1))]
    r["review"]["unresolved"] = ["Prior-work summary only; cited papers were not supplied, so no experimental conditions or outcomes were imported."]
    set_classification(r, "unresolved")
    reactions.append(r)

r = base_reaction("scheme1_current_generic", evs(1), "generic_scheme", ["other"], "Scheme 1B", "asymmetric construction of all-carbon quaternary centers", "allene–anhydride hydroacylation")
r["participants"] = [participant("allene", "m_gen_disuballene", "reactant", evs(1)), participant("anhydride", "m_gen_anhydride", "reactant", evs(1)), participant("product", "m_gen_current_product", "product", evs(1))]
r["notes"] = [ef("8 racemic examples up to 99% yield; 23 chiral examples up to 96% yield, 98% ee; three regioisomer classes not observed", {"racemic_examples": 8, "chiral_examples": 23, "not_observed": ["three depicted regioisomer classes"]}, evs(1))]
r["review"]["unresolved"] = ["Generic R/Ar/R1 identities remain templates; no concrete experiment is created from the overview counts."]
set_classification(r, "hydroacylation")
reactions.append(r)

# Table 1, twelve distinct optimization rows.
table1 = [
    (1,"DPPM","m_thf",0,100,2.5),(2,"DPPE","m_thf",31,63,2.5),(3,"DPPP","m_thf",0,100,2.5),
    (4,"DPPF","m_thf",8,84,2.5),(5,"BINAP","m_thf",21,75,2.5),(6,"PCy3","m_thf",0,100,5.0),
    (7,"DPPE","m_cyclohexane",0,100,2.5),(8,"DPPE","m_dioxane",15,78,2.5),(9,"DPPE","m_dcm",10,79,2.5),
    (10,"DPPE","m_toluene",24,62,2.5),(11,"DPPE","m_thf",49,4,2.5),(12,"DPPE","m_thf",66,None,2.5),
]
for entry, ligand, solvent, yld, rec, lpct in table1:
    kwargs = dict(allene_amount=("0.1 mmol",0.1,"mmol"), anhydride_equiv=1.5)
    notes=[]
    if entry == 11:
        kwargs = dict(allene_amount=("0.5 mmol",0.5,"mmol"), anhydride_equiv=1.5)
    if entry == 12:
        kwargs = dict(allene_equiv=2.0, anhydride_amount=("0.5 mmol",0.5,"mmol"))
        notes=["Table 1 entry 12: the common footnote still prints solvent (1.0 mL) and the scheme header prints 0.1 M after the scale change."]
    hydro(f"table1_entry{entry:02d}", [2], categories=["optimization"], location=f"Table 1 entry {entry}", ligand=ligand, solvent=solvent,
          ligand_pct=lpct, solvent_volume=("1.0 mL",1.0,"mL"), concentration="0.1 M", yield_pct=yld, yield_basis="NMR yield (internal standard not specified)",
          recovery_pct=rec, **kwargs, notes=notes)

# Table 2, eight independent racemic scope experiments.
table2 = [
    ("3aa",67,"m_1a","m_2a",6,2.0),("3ab",99,"m_1a","m_2b",6,2.0),("3ac",89,"m_1a","m_2c",6,2.0),("3ad",74,"m_1a","m_2d",9,2.0),
    ("3ae",95,"m_1a","m_2e",6,2.0),("3ia",90,"m_1i","m_2a",6,1.5),("3ja",74,"m_1j","m_2a",6,1.5),("3ka",78,"m_1k","m_2a",6,1.5),
]
for label,yld,a,b,t,aeq in table2:
    hydro(f"table2_rac_{label}",[2],categories=["substrate_scope"],location=f"Table 2, rac-{label}",allene=a,anhydride=b,product=f"m_rac_{label}",ligand="DPPE",
          allene_equiv=aeq,anhydride_amount=("0.5 mmol",0.5,"mmol"),solvent_volume=("10.0 mL",10.0,"mL"),concentration="0.1 M",time_h=t,yield_pct=yld)

# Table 3 ligand screening. Negative ee signs are visible in the source figure and retained.
table3 = [
    ("R_BINAP","(R)-BINAP",74,-79),("R_MeO_BIPHEP","(R)-MeO-BIPHEP",65,-84),("R_C3_Tunephos","(R)-C3-Tunephos",64,-86),
    ("R_Garphos","(R)-Garphos",62,-78),("S_Synphos","(S)-Synphos",63,80),("R_Segphos","(R)-Segphos",68,-87),
    ("RR_BDPP","(R,R)-BDPP",60,-65),("R_SDP","(R)-SDP",None,None),("RR_BPE","(R,R)-BPE",34,91),
]
for code, ligand, yld, ee_val in table3:
    hydro(f"table3_{code}",[2],categories=["optimization"],location=f"Table 3, {ligand}",product="m_R_3aa",ligand=ligand,
          allene_amount=("0.1 mmol",0.1,"mmol"),anhydride_equiv=1.5,catalyst=5,ligand_pct=5,silane_equiv=4,
          solvent_volume=("1.0 mL",1.0,"mL"),concentration="0.1 M",temperature=40,time_h=12,yield_pct=yld,
          yield_basis="1H NMR using dibromomethane as internal standard" if yld is not None else None,
          recovery_pct=84 if ligand=="(R)-SDP" else None,ee=ee_val,ee_raw=(f"−{abs(ee_val)}% ee" if ee_val and ee_val<0 else None),
          classification="unresolved" if ligand=="(R)-SDP" else "hydroacylation")

# Table 4, six distinct optimization rows.
table4 = [(1,4.0,12,None,68,93,5,5,0.1),(2,4.0,5,None,92,96,5,5,0.1),(3,3.0,5,None,89,94,5,5,0.1),
          (4,2.0,5,15,82,94,5,5,0.1),(5,1.0,5,50,50,94,5,5,0.1),(6,3.0,6,None,86,95,2.5,2.5,1.0)]
for entry,x,t,rec,yld,ee_val,cp,lp,scale in table4:
    hydro(f"table4_entry{entry:02d}",[3],categories=["optimization"],location=f"Table 4 entry {entry}",product="m_R_3aa",ligand="(R,R)-BPE",
          allene_amount=(f"{scale} mmol",scale,"mmol"),anhydride_equiv=1.5,catalyst=cp,ligand_pct=lp,silane_equiv=x,
          solvent_volume=("1.0 mL",1.0,"mL") if entry<6 else ("10.0 mL",10.0,"mL"),concentration="0.1 M",time_h=t,
          yield_pct=yld,yield_basis="1H NMR using MeNO2 as internal standard",recovery_pct=rec,
          recovery_basis="1H NMR using MeNO2 as internal standard" if rec is not None else None,ee=ee_val)

# Table 5, 23 independent enantioselective scope experiments.
table5_allene = [("3aa",87,95,"1a",None),("3ba",83,97,"1b",None),("3ca",81,95,"1c",None),("3da",86,95,"1d","b"),("3ea",80,94,"1e",None),
                 ("3fa",96,99,"1f",None),("3ga",75,92,"1g",None),("3ha",93,94,"1h","b"),("3ia",87,95,"1i","b"),("3ja",76,93,"1j","b"),("3ka",78,96,"1k","b")]
for label,yld,ee_val,allene_label,foot in table5_allene:
    if foot == "b":
        kwargs=dict(allene_equiv=1.5,anhydride_amount=("1.0 mmol",1.0,"mmol"),silane_equiv=1.5)
    else:
        kwargs=dict(allene_amount=("1.0 mmol",1.0,"mmol"),anhydride_equiv=1.5,silane_equiv=3.0)
    hydro(f"table5_R_{label}",[3],categories=["substrate_scope"],location=f"Table 5, (R)-{label}",allene=f"m_{allene_label}",product=f"m_R_{label}",ligand="(R,R)-BPE",
          solvent_volume=("10.0 mL",10.0,"mL"),concentration="0.1 M",yield_pct=yld,ee=ee_val,**kwargs)

table5_anh = [("R","3ab",94,96,"2b","c"),("R","3ac",90,96,"2c","c"),("R","3ad",81,95,"2d","c"),("R","3ae",98,97,"2e","c"),
              ("R","3af",92,99,"2f","c"),("S","3af",91,99,"2f","d"),("R","3ag",96,98,"2g","c"),("R","3ah",92,91,"2h","c"),
              ("R","3ai",85,93,"2i","c"),("S","3ai",70,92,"2i","d"),("R","3aj",93,95,"2j","c"),("S","3aj",65,97,"2j","d")]
for conf,label,yld,ee_val,anh_label,foot in table5_anh:
    ligand="(S,S)-BPE" if conf=="S" else "(R,R)-BPE"
    aeq=1.5 if foot=="d" else 2.0
    amount=0.5 if foot=="d" else 1.0
    hydro(f"table5_{conf}_{label}",[3],categories=["substrate_scope"],location=f"Table 5, ({conf})-{label}",anhydride=f"m_{anh_label}",product=f"m_{conf}_{label}",ligand=ligand,
          allene_equiv=aeq,anhydride_amount=(f"{amount} mmol",amount,"mmol"),solvent_volume=("10.0 mL",10.0,"mL"),concentration="0.1 M",yield_pct=yld,ee=ee_val)

# Scheme 2a: current-paper gram-scale experiment, linked to the page-3 prose and page-4 arrow.
r = hydro("scheme2a_gram_scale",[3,4],categories=["scale_up"],location="Scheme 2a gram-scale synthesis",allene="m_1k",product="m_R_3ka",ligand="(R,R)-BPE",
          allene_equiv=1.5,anhydride_amount=("10 mmol",10,"mmol"),catalyst=2.5,ligand_pct=2.5,silane_equiv=1.5,
          solvent_volume=None,concentration="0.1 M",yield_pct=79,ee=96,notes=["Allene 1k is printed as 15 mmol; product amount is 2.4974 g."])
r["participants"][0]["charges"]["amount"] = ef("15 mmol",q(15,"mmol","15 mmol"),r["evidence_ids"])
r["reaction_scale"] = ef("2a, 10 mmol",{"reference_molecule_id":"m_2a","amount":q(10,"mmol","10 mmol")},r["evidence_ids"])
add_outcome(r,"product_amount","product_amount","m_R_3ka","2.4974 g",q(2.4974,"g","2.4974 g"))


def derivatization(rid, location, product, kind, yield_pct, ee_val, participants_data, solvent_ids, temp, time_h, atmosphere=None):
    evidence=evs(3,4)
    r=base_reaction(rid,evidence,"specific_experiment",["derivatization"],location,"synthetic transformation of (R)-3ka",kind)
    r["participants"]=[participant("substrate","m_R_3ka","reactant",evidence)]
    for i,(mid,role,equiv,molpct) in enumerate(participants_data):
        r["participants"].append(participant(f"agent{i+1}",mid,role,evidence,ref="m_R_3ka",equiv=equiv,molpct=molpct))
    for i,mid in enumerate(solvent_ids):
        r["participants"].append(participant(f"solvent{i+1}",mid,"solvent",evidence))
    if atmosphere:
        r["participants"].append(participant("atmosphere_agent",atmosphere,"reagent",evidence))
    r["participants"].append(participant("product",product,"product",evidence))
    r["solvent_conditions"]["composition"]=ef([molecules[x]["names"]["original"]["raw_value"] for x in solvent_ids],solvent_ids,evidence)
    if len(solvent_ids)>1:
        ratio="4:1" if rid.endswith("oxime") else "15:1"
        r["solvent_conditions"]["mixture_ratio"]=ef(ratio,ratio,evidence)
    step=st.new_step("reaction",evidence);step["participant_ids"]=[p["participant_id"] for p in r["participants"]]
    step["operation"]=inf("react under Scheme 2 conditions",evidence,"The scheme gives conditions but not an addition sequence.")
    step["temperature"]=ef(temp,temp if temp in ("r.t.","reflux") else q(temp,"°C",f"{temp} °C"),evidence)
    step["time"]=ef(f"{time_h} h",q(time_h,"h",f"{time_h} h"),evidence)
    if atmosphere:
        step["atmosphere"]=ef(molecules[atmosphere]["names"]["original"]["raw_value"]+" balloon",molecules[atmosphere]["names"]["original"]["raw_value"]+" balloon",evidence)
    r["steps"]=[step]
    r["procedure"]["original_text"]=ef(location,location,evidence)
    add_outcome(r,"yield","yield",product,f"{yield_pct}%",q(yield_pct,"%",f"{yield_pct}%"))
    add_outcome(r,"ee","selectivity",product,f"{ee_val}% ee",q(ee_val,"%",f"{ee_val}% ee"),scope={"metric":"ee","major_configuration":"source product label"})
    set_classification(r,kind)
    r["review"]["unresolved"]=["Yield and ee measurement methods are not reported on the supplied pages.","Exact procedure, workup and purification are unavailable because no SI was supplied."]
    reactions.append(r)


derivatization("scheme2b_oxime","Scheme 2b, arrow a","m_RZ_4","oxime",63,95,
               [("m_hydroxylamine_hcl","reagent",5.0,None),("m_naoac","additive",7.0,None)], ["m_etoh","m_water"],"reflux",72)
derivatization("scheme2b_suzuki","Scheme 2b, arrow b","m_R_5","suzuki",90,96,
               [("m_nitroboronic","reactant",1.1,None),("m_pd_dppf","catalyst",None,10),("m_k2co3","reagent",2.0,None)], ["m_dmso"],80,12,"m_ar")
derivatization("scheme2b_wacker","Scheme 2b, arrow c","m_R_6","wacker",63,93,
               [("m_pd_phcn","catalyst",None,12),("m_cucl2_h2o","catalyst",None,12),("m_agno2","additive",None,6)], ["m_tbuoh","m_me_no2"],"r.t.",36,"m_o2")
derivatization("scheme2b_hydrogenation","Scheme 2b, arrow d","m_R_7","hydrogenation",98,96,
               [("m_pdc","catalyst",None,2)], ["m_etoac"],"r.t.",6,"m_h2")

# Scheme 3 is a proposed catalytic cycle, not a proved experiment.
r=base_reaction("scheme3_proposed_mechanism",evs(4),"proposed_mechanism",["mechanistic_study"],"Scheme 3","proposed catalytic cycle","CuH hydrocupration, allylcopper acylation, β-O elimination, CuH regeneration")
r["participants"]=[participant("allene","m_1d","reactant",evs(4)),participant("anhydride","m_2a","reactant",evs(4)),participant("product","m_R_3da","product",evs(4)),participant("silane","m_silane","reagent",evs(4))]
r["notes"]=[ef("hydrocupration gives allylcopper E-II; six-membered chair-like TS; β-O elimination; Cu-OBz reacts with Me(MeO)2SiH",{"status":"author-proposed mechanism","intermediates":["I","E-II","TS","III","IV"]},evs(4))]
r["review"]["unresolved"]=["Author-proposed mechanism and computed intermediates are not treated as experimentally observed species."]
set_classification(r,"hydroacylation")
reactions.append(r)

bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": f"{DOC}:si_availability", "document_id": DOC, "evidence_ids": evs(1), "kind": "unresolved_si_reference", "text": "Page 1 states that ESI is available, but no SI file was explicitly supplied or associated in this workspace."},
    {"context_id": f"{DOC}:dft_schemes4_6", "document_id": DOC, "evidence_ids": evs(4,5), "kind": "computational_mechanistic_analysis", "text": "Schemes 4–6 report DFT transition structures, allylcopper isomerization and the proposed origin of enantioselectivity. These computed structures are not independent synthetic experiments."},
    {"context_id": f"{DOC}:page6_nonreaction", "document_id": DOC, "evidence_ids": evs(6), "kind": "nonreaction_content", "text": "Page 6 contains the conclusion, conflicts statement, acknowledgements and references; it does not add a current-paper reaction arrow, table, procedure or characterization-only peak list."},
    {"context_id": f"{DOC}:characterization_inventory", "document_id": DOC, "evidence_ids": evs(1,2,3,4,5,6), "kind": "characterization_inventory", "text": "No standalone mass-spectrum peak table, molecular-weight table or fragment-relative-abundance table occurs on pages 1–6."},
]

st.write(RUN / "submitted.json", bundle)

# Raw inventories are provenance, not formal records.
page_inventory = {
    "document_id": DOC,
    "pages": {
        "1": {"regions": ["title/abstract", "Scheme 1A-a–e", "Scheme 1B", "ESI notice"], "reaction_rows": 6},
        "2": {"regions": ["Table 1 entries 1–12", "Table 2 eight products", "Table 3 nine ligand screens"], "reaction_rows": 29},
        "3": {"regions": ["Table 4 entries 1–6", "Table 5 twenty-three products", "Scheme 2 prose cross-reference"], "reaction_rows": 29},
        "4": {"regions": ["Scheme 2a", "Scheme 2b arrows a–d", "Scheme 3", "Schemes 4–5 computational"], "reaction_rows": 6},
        "5": {"regions": ["Scheme 6 computational transition-state analysis", "conclusion continuation"], "reaction_rows": 0},
        "6": {"regions": ["conclusion", "conflicts", "acknowledgements", "references 1–12"], "reaction_rows": 0},
    },
    "limitations": ["Inventory counts source arrows/rows, not chemical accuracy.", "Same Scheme 2 experiments are cross-referenced in page-3 prose and counted once formally."],
}
st.write(RUN / "raw_model/page_inventory.json", page_inventory)

stereo_checks = {}
for mid in ["m_1a", "m_2a", "m_cu_oac2", "m_silane", "m_thf", "m_cyclohexane", "m_dioxane", "m_dcm", "m_toluene"]:
    stereo_checks[mid] = st.stereo_inspect(st.value(molecules[mid]["smiles"]))
st.write(RUN / "raw_model/stereo_checks.json", stereo_checks)

mentions = st.source_mentions(dataset, ["1a", "2a", "3aa", "3ka"], context=180)
st.write(RUN / "raw_model/source_mentions.json", mentions)

record_ids_by_page = {str(i): [] for i in range(1, 7)}
for reaction in reactions:
    pages = sorted({int(e.split(":page",1)[1].split(":",1)[0]) for e in reaction["evidence_ids"] if ":page" in e})
    for page in pages:
        if 1 <= page <= 6:
            record_ids_by_page[str(page)].append(reaction["reaction_id"])
st.write(RUN / "raw_model/record_ids_by_page.json", record_ids_by_page)

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 1,
            "state": "attempted",
            "regions": [
                {"region": "title, abstract and article metadata", "finding": "Current paper identified as an enantioselective copper-catalyzed allene–anhydride hydroacylation; DOI and ESI availability statement observed.", "evidence_ids": evs(1)},
                {"region": "Scheme 1A-a through 1A-e", "finding": "Five prior-work/general chemistry arrows classified as literature examples, not experiments of this paper.", "evidence_ids": evs(1)},
                {"region": "Scheme 1B", "finding": "One current-work generic reaction scheme; overview counts and three depicted regioisomer classes marked not observed, without creating concrete experiments.", "evidence_ids": evs(1)},
            ],
            "record_ids": record_ids_by_page["1"],
            "unresolved": ["The page says ESI is available, but no SI file was explicitly provided; no SI association was inferred.", "Generic R/Ar/R1 structures remain templates."],
        },
        {
            "page": 2,
            "state": "attempted",
            "regions": [
                {"region": "Table 1 entries 1–12 and footnotes a–d", "finding": "Twelve independent optimization experiments recorded, including zero-yield rows, starting-material recovery and entry-specific scale/ligand changes.", "evidence_ids": evs(2)},
                {"region": "Table 2, eight drawn racemic products and footnotes a–c", "finding": "Eight independent substrate-scope experiments recorded; racemic product identities remain distinct from enantioenriched products.", "evidence_ids": evs(2)},
                {"region": "Table 3, nine ligand drawings/results and footnote", "finding": "Nine independent ligand-screening experiments recorded; signed ee values visible in the figure are retained and the NMR-yield basis is linked to dibromomethane.", "evidence_ids": evs(2)},
            ],
            "record_ids": record_ids_by_page["2"],
            "unresolved": ["Table 1 entry 12 prints the common 1.0 mL footnote and 0.1 M header after a scale change; the concentration basis is ambiguous.", "Table 1 recovery methods, Table 2 yield methods and Table 3 ee methods are not reported on this page.", "No SI procedures/workups/purifications were supplied."],
        },
        {
            "page": 3,
            "state": "attempted",
            "regions": [
                {"region": "Table 4 entries 1–6 and footnotes a–d", "finding": "Six independent optimization experiments recorded with silane loading, time, recovery, NMR yield and ee; entry 6 retains the 1.0 mmol scale and reduced catalyst/ligand loadings.", "evidence_ids": evs(3)},
                {"region": "Table 5, allene-substrate and anhydride-substrate panels plus footnotes a–d", "finding": "Twenty-three independent scope experiments recorded, including separately registered R/S enantiomers and source-assigned (R)-3ea X-ray configuration (CCDC 1939837).", "evidence_ids": evs(3)},
                {"region": "prose introducing Scheme 2", "finding": "The gram-scale result and four downstream transformations are cross-references to the arrows on page 4 and are linked to the same five records, not duplicated.", "evidence_ids": evs(3)},
            ],
            "record_ids": record_ids_by_page["3"],
            "unresolved": ["Table 5 yield and ee measurement methods are not stated on the supplied pages.", "Most complex source drawings were not promoted to exact isomeric SMILES; R/S assignments are preserved as source labels.", "No SI procedures/workups/purifications were supplied."],
        },
        {
            "page": 4,
            "state": "attempted",
            "regions": [
                {"region": "Scheme 2a gram-scale arrow", "finding": "One 10 mmol benzoic-anhydride/15 mmol allene experiment recorded with all printed loadings, 2.4974 g product, 79% yield and 96% ee.", "evidence_ids": evs(4)},
                {"region": "Scheme 2b arrows a–d and reaction-condition line", "finding": "Four independent derivatization experiments recorded: oxime formation, Suzuki coupling, Wacker-type oxidation and hydrogenation; every printed reagent loading, solvent ratio, atmosphere, time and temperature is structured.", "evidence_ids": evs(4)},
                {"region": "Scheme 3", "finding": "Author-proposed catalytic cycle recorded as a proposed mechanism, not as a proved experiment.", "evidence_ids": evs(4)},
                {"region": "Schemes 4 and 5 with adjacent DFT prose", "finding": "Computed hydrocupration transition structures and allylcopper isomerizations retained in source_context, not promoted to experimental reactions.", "evidence_ids": evs(4)},
            ],
            "record_ids": record_ids_by_page["4"],
            "unresolved": ["Scheme 2 yield/ee measurement methods and exact procedures/workups are not reported in the supplied main text.", "Complex reactant/product identities in the four derivatizations retain drawing-linked unresolved structures, causing documented structure_gap advisories."],
        },
        {
            "page": 5,
            "state": "attempted",
            "regions": [
                {"region": "Scheme 6a–c and explanatory prose", "finding": "DFT transition-state energy comparison and quadrant analysis only; no independent synthetic experiment, condition table, procedure or observed product measurement.", "evidence_ids": evs(5)},
                {"region": "conclusion opening", "finding": "Narrative summary of reactions already represented; no new reaction arrow or row.", "evidence_ids": evs(5)},
            ],
            "record_ids": record_ids_by_page["5"],
            "unresolved": ["Other less-favorable computed structures are referred to the unavailable ESI; they were not reconstructed."],
        },
        {
            "page": 6,
            "state": "attempted",
            "regions": [
                {"region": "conclusion, conflicts and acknowledgements", "finding": "Narrative/non-experimental content; reproduction acknowledgements do not create additional experiments.", "evidence_ids": evs(6)},
                {"region": "references 1–12", "finding": "Citations only; cited-paper chemistry was not imported because the cited papers were not supplied.", "evidence_ids": evs(6)},
                {"region": "full-page characterization inventory", "finding": "No standalone mass-spectrum peak list, molecular-weight table or fragment-relative-abundance table was present.", "evidence_ids": evs(6)},
            ],
            "record_ids": record_ids_by_page["6"],
            "unresolved": ["References continue beyond the required page range; no coverage is claimed for those pages in this session."],
        },
    ],
}
st.write(RUN / "coverage/pages_0001_0006.json", coverage)

# Mark only pages 1-6 attempted; page 7 remains outside this session's claimed coverage.
tasks=json.loads((RUN/"tasks.json").read_text())
for task in tasks:
    if task["document_id"]==DOC and task["evidence_id"].startswith(f"{DOC}:page"):
        page=int(task["evidence_id"].split(":page",1)[1].split(":",1)[0])
        if 1 <= page <= 6:
            task["state"]="attempted"
            task["reason"]="Page image and extracted text inspected; formal mappings are in coverage/pages_0001_0006.json."
st.write(RUN/"tasks.json",tasks)

print(json.dumps({"molecules":len(bundle["molecules"]),"reactions":len(bundle["reactions"]),"specific_experiments":sum(r["record_type"]=="specific_experiment" for r in reactions)},ensure_ascii=False))
