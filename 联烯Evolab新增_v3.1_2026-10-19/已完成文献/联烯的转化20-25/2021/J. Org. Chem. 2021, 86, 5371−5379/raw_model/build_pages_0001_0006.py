#!/usr/bin/env python3
"""Build the source-derived v3.1 bundle for PDF pages 1--6.

This is a reproducible extraction artifact, not a chemistry validator.  It uses
only evidence already present in runs/current/dataset.json and deliberately
leaves uncertain structures unresolved.
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

DATASET = ROOT / "runs/current/dataset.json"
SUBMITTED = ROOT / "runs/current/submitted.json"
RAW = ROOT / "runs/current/raw_model"
COVERAGE = ROOT / "runs/current/coverage/pages_0001_0006.json"
DOC = "doc_2c129c777a77980b"
TITLE = "Palladium-Catalyzed Oxidative Coupling of the Allenic C−H Bond with α-Diazo Esters: Synthesis of [3]Dendralenes"


def ev(page: int, image: bool = False) -> str:
    return f"{DOC}:page{page}:{'image' if image else 'text'}"


def F(raw=None, normalized=None, evidence=(), method="original_extraction", derivation=None,
      status=None, verification=None, alternatives=()):
    if normalized is None and raw is not None:
        normalized = raw
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation,
                   status=status, verification=verification, alternatives=alternatives)


def q(number, unit, **extra):
    return {"number": number, "unit": unit, **extra}


def rng(minimum, maximum, unit, raw_text):
    return {"minimum": minimum, "maximum": maximum, "unit": unit, "raw_text": raw_text}


dataset = json.loads(DATASET.read_text())
assert not dataset.get("experiments"), "This batch must not reuse old experimental records"
assert not dataset.get("candidates"), "Unexpected mechanically generated candidates"
dataset["documents"][0]["title"] = TITLE
DATASET.write_text(json.dumps(dataset, ensure_ascii=False, indent=2) + "\n")

bundle = st.empty(copy.deepcopy(dataset["documents"]))
molecules = {}


def make_template(variables, evidence, fixed, reason):
    return {
        "fixed_structure": F(fixed, fixed, evidence),
        "variables": [
            {
                "label": label,
                "attachment_points": F(f"Source drawing attachment for {label}", {"label": label}, evidence),
                "definition": F(definition, definition, evidence),
                "allowed_values": F(),
                "constraints": F(),
            }
            for label, definition in variables
        ],
        "instance_mapping": [],
        "reason": F(reason, reason, evidence),
    }


def mol(label, name=None, pages=(3,), smiles=None, chirality="unresolved",
        geometry="unresolved", representation=None, double_bonds=None,
        mixture=None, assignment=None, template=None, note=None):
    mid = f"{DOC}:mol_{label.replace('′', 'prime').replace('-', '_').replace('/', '_')}"
    if mid in molecules:
        return mid
    evidence = [ev(p) for p in pages]
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, evidence)
    source_name = name or label
    m["names"]["original"] = F(source_name, source_name, evidence)
    m["names"]["english"] = F(source_name, source_name, evidence)
    m["names"]["abbreviations"] = F(label, [label], evidence)
    m["names"]["aliases"] = F()
    m["cas_numbers"] = F()
    m["external_properties"] = []
    m["generic_template"] = None
    m["identity_note"] = F(note, note, evidence) if note else F()
    if template:
        m["representation_type"] = "generic_template"
        m["generic_template"] = template
        m["smiles"] = F()
        rep_status = "not_encoded"
    elif smiles:
        diag = st.stereo_inspect(smiles)
        canonical = diag["canonical_isomeric_smiles"]
        m["representation_type"] = "smiles"
        m["smiles"] = F(smiles, canonical, evidence, method="structure_derivation",
                         derivation="Candidate structure encoded from the source name/drawing; RDKit canonicalization does not verify source identity.")
        rep_status = representation or "encoded"
        m["stereochemistry"]["isomeric_smiles"] = F(smiles, canonical, evidence,
            method="structure_derivation", derivation="Canonical isomeric SMILES from the candidate graph.")
        m["stereochemistry"]["identity_key"] = F(None, diag["stereo_graph_key"], evidence,
            method="structure_derivation", derivation="standard stereo graph key for the candidate isomeric SMILES.")
        if double_bonds is None and diag["double_bonds"]:
            double_bonds = diag["double_bonds"]
        if chirality == "achiral":
            m["stereochemistry"]["atom_centers"] = F(None, diag["atom_centers"], evidence,
                method="structure_derivation", derivation="No assigned tetrahedral center in the candidate graph; source drawing/name also gives none.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = F()
        rep_status = representation or "not_encoded"
    if not smiles:
        m["stereochemistry"]["isomeric_smiles"] = F()
        m["stereochemistry"]["identity_key"] = F()
    m["stereochemistry"]["chirality_status"] = F(None, chirality, evidence, method="inference",
        derivation="Status preserves the source assignment or its absence; unresolved is not treated as achiral or racemic.")
    m["stereochemistry"]["geometry_status"] = F(None, geometry, evidence, method="inference",
        derivation="Geometry status follows the product label/drawing; unassigned elements remain unresolved.")
    m["stereochemistry"]["representation_status"] = F(None, rep_status, evidence, method="inference",
        derivation="Reports whether the retained source stereochemistry is encoded in the candidate representation.")
    m["stereochemistry"]["double_bonds"] = F("source name/drawing", double_bonds, evidence) if double_bonds else F()
    m["stereochemistry"]["mixture"] = F("source-reported mixture", mixture, evidence) if mixture else F()
    m["stereochemistry"]["assignment"] = F("source assignment", assignment, evidence) if assignment else F()
    m["stereochemistry"]["axial"] = F()
    if m["stereochemistry"]["atom_centers"]["normalized_value"] is None:
        m["stereochemistry"]["atom_centers"] = F()
    molecules[mid] = m
    return mid


def generic(label, name, pages, fixed, variables, reason):
    evidence = [ev(p) for p in pages]
    return mol(label, name, pages, template=make_template(variables, evidence, fixed, reason),
               chirality="unresolved", geometry="unresolved", representation="not_encoded", note=reason)


# Generic structures from source drawings.
M_ALLENE = generic("allene_1", "generic allene 1", (1, 3), "substituted allene bearing an allylic C−H and R1", [("R", "allene terminal substituents"), ("R1", "side-chain substituent")], "R-group drawing in Schemes 1e/2; not a single concrete molecule.")
M_DIAZO = generic("diazo_2", "generic α-diazo ester 2", (1, 3), "α-diazocarboxylate with ester substituent R2", [("R2", "ester substituent")], "R-group drawing in Schemes 1e/2.")
M_DEND = generic("dendralene_3", "generic [3]dendralene 3", (1, 3), "cross-conjugated triene bearing CO2R1/CO2R2 substituents", [("R", "allene-derived substituents"), ("R1", "side chain"), ("R2", "ester substituent")], "Generic product drawing; stereochemical identity varies by instance.")
M_PROPARGYL = generic("propargyl_alcohol_generic", "substituted propargyl alcohol", (4,), "substituted propargyl alcohol R1R2C(OH)C≡CH", [("R1", "substituent"), ("R2", "substituent")], "Generic reactant in the allene preparation scheme.")
M_ALCOHOL = generic("alcohol_generic", "alcohol used for transesterification", (4,), "R2OH", [("R2", "alcohol residue")], "Generic alcohol in the α-diazo ester preparation scheme.")
M_BETA = mol("beta_ketoester", "ethyl 2-methyl-3-oxobutanoate", (4,), smiles="CCOC(=O)C(C)C(C)=O", chirality="unspecified", geometry="not_applicable", representation="partial")
M_2PRIME = generic("2prime_generic", "generic 2-methyl-3-oxobutanoate ester 2′", (4,), "R2O2C-CH(CH3)-COCH3", [("R2", "alcohol-derived ester group")], "Generic intermediate 2′ in the preparation procedure.")

# Core concrete and repeatedly used substances.
M_1A = mol("1a", "ethyl 5-methylhexa-3,4-dienoate", (2, 3, 5), smiles="CCOC(=O)CC=C=C(C)C", chirality="achiral", geometry="not_applicable", representation="partial")
M_2A = mol("2a", "benzyl 2-diazopropanoate", (2, 3, 5), smiles="CC(=[N+]=[N-])C(=O)OCc1ccccc1", chirality="achiral", geometry="not_applicable", representation="partial")
M_3AA = mol("3aa", "6-benzyl 1-ethyl (E)-5-methylene-4-(propan-2-ylidene)hex-2-enedioate", (2, 3, 5), smiles="CCOC(=O)/C=C/C(=C(C)C)C(=C)C(=O)OCc1ccccc1", chirality="achiral", geometry="specified", double_bonds=[{"assignment": "E", "element": "hex-2-ene"}])
M_3P = mol("3prime", "diethyl (2E,6E)-4,5-di(propan-2-ylidene)octa-2,6-dienedioate", (2,), chirality="unresolved", geometry="specified", double_bonds=[{"assignment": "2E,6E", "basis": "Table 1 drawing and later name"}])
M_1AD6 = mol("1a_d6", "ethyl 5-(methyl-d3)hexa-3,4-dienoate-6,6,6-d3", (3, 4), chirality="achiral", geometry="not_applicable")
M_3AAD6 = mol("3aa_d6", "6-benzyl 1-ethyl (E)-5-methylene-4-(propan-2-ylidene-d6)hex-2-enedioate", (3,), chirality="achiral", geometry="specified", double_bonds=[{"assignment": "E", "element": "hex-2-ene"}])
M_4 = mol("4", "benzyl (E)-2-(2-(3-ethoxy-3-oxoprop-1-en-1-yl)-3,3-dimethyloxiran-2-yl)acrylate", (3,), chirality="unspecified", geometry="specified", double_bonds=[{"assignment": "E", "basis": "product label"}], representation="not_encoded")


def simple(label, name, pages=(2,), smiles=None, chirality="achiral", geometry="not_applicable"):
    return mol(label, name, pages, smiles=smiles, chirality=chirality, geometry=geometry,
               representation="encoded" if smiles else "not_encoded")


M_PDOAC = simple("PdOAc2", "palladium(II) acetate", (2, 3, 5), "CC(=O)[O-].CC(=O)[O-].[Pd+2]")
M_PDCL2 = simple("PdCl2", "palladium(II) chloride", (2,), "Cl[Pd]Cl")
M_PDL = mol("PdOAc2_L", "Pd(OAc)2/1,2-bis(phenylmethanesulfinyl)ethane complex", (2,), chirality="unresolved", geometry="not_applicable", note="Source calls this the white catalyst; complete coordination structure is not encoded.")
M_BQ = simple("BQ", "p-benzoquinone", (2,), "O=C1C=CC(=O)C=C1")
M_DMBQ = simple("2,6_DMBQ", "2,6-dimethylbenzoquinone", (2, 3, 5), "CC1=CC(=O)C=CC1=O")
M_DMSO = simple("DMSO", "dimethyl sulfoxide", (2,), "CS(C)=O")
M_SALEN = mol("salenCrCl", "(salen)CrCl", (2,), chirality="unresolved", geometry="unresolved", note="Ligand substitution pattern/coordination structure is not specified in the table.")
M_THF = simple("THF", "tetrahydrofuran", (2, 3, 5), "C1CCOC1")
M_DIOX = simple("dioxane", "dioxane", (2,), "O1CCOCC1")
M_DME = simple("DME", "1,2-dimethoxyethane", (2,), "COCCOC")
M_MENO2 = simple("MeNO2", "nitromethane", (2,), "C[N+](=O)[O-]")
M_ACOH = simple("AcOH", "acetic acid", (3,), "CC(=O)O")
M_MCPBA = mol("m_CPBA", "m-chloroperoxybenzoic acid", (3,), chirality="achiral", geometry="not_applicable")
M_DCM = simple("DCM", "dichloromethane", (4, 5), "ClCCl")
M_EA = simple("EA", "ethyl acetate", (4, 5), "CCOC(C)=O")
M_PE = mol("PE", "petroleum ether", (4, 5), chirality="unresolved", geometry="not_applicable", note="Named solvent mixture; composition is not uniquely specified.")
M_SILICA = mol("silica", "silica gel", (4, 5), chirality="not_applicable", geometry="not_applicable", representation="not_applicable")
M_TOL = simple("toluene", "toluene", (4,), "Cc1ccccc1")
M_HCL = simple("HCl", "hydrochloric acid", (4,), "Cl")
M_BRINE = mol("brine", "brine", (4,), chirality="not_applicable", geometry="not_applicable", representation="not_applicable")
M_NA2SO4 = simple("Na2SO4", "sodium sulfate", (4, 5), "[Na+].[Na+].[O-]S(=O)(=O)[O-]")
M_AMINAL = mol("amino_acetal", "1,1-dimethoxy-N,N-dimethylethan-1-amine", (4,), chirality="achiral", geometry="not_applicable")
M_TMSAC = simple("ethynyl_TMS", "ethynyltrimethylsilane", (4,), "C#C[Si](C)(C)C")
M_NBULI = simple("n_BuLi", "n-butyllithium", (4,), "[Li]CCCC")
M_D6AC = mol("d6_acetone", "acetone-d6", (4,), chirality="achiral", geometry="not_applicable")
M_K2CO3 = simple("K2CO3", "potassium carbonate", (4,), "[K+].[K+].[O-]C([O-])=O")
M_MEOH = simple("MeOH", "methanol", (4,), "CO")
M_D5 = mol("d6_5", "[d6]-5", (4,), chirality="unresolved", geometry="not_applicable", note="Intermediate structure is shown; exact isotopic SMILES is not encoded.")
M_ETCO2H = simple("EtCO2H", "propionic acid", (4,), "CCC(=O)O")
M_ORTHO = simple("triethyl_orthoacetate", "triethyl orthoacetate", (4,), "CCOC(C)(OCC)OCC")
M_PABSA = mol("p_ABSA", "p-acetamidobenzenesulfonyl azide", (4,), chirality="achiral", geometry="not_applicable")
M_DBU = mol("DBU", "1,8-diazabicyclo[5.4.0]undec-7-ene", (4,), chirality="achiral", geometry="not_applicable")
M_MECN = simple("MeCN", "acetonitrile", (4,), "CC#N")
M_NH4CL = simple("NH4Cl", "ammonium chloride", (4,), "[Cl-].[NH4+]")


def change(change_id, a, b, before, after, kind, evidence, derivation):
    return F(None, {"change_id": change_id, "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
                    "before_order": before, "after_order": after, "change_type": kind}, evidence,
             method="structure_derivation", derivation=derivation)


def set_classification(r, mode, evidence):
    c = r["classification"]
    c["mapped_reaction_smiles"] = F()
    c["model_votes"] = []
    if mode == "coupling":
        c["bond_changes"] = [
            change("b1", "C", "C", 0, 1, "formed", evidence, "Net carbon-carbon connection inferred from allene/diazo-to-dendralene drawings."),
            change("b2", "C", "C", 1, 2, "order_increased", evidence, "New alkene in the cross-conjugated product inferred from source drawings."),
        ]
        labels = [
            {"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1"]},
            {"category_number": 13, "code": "DB-C-C", "supporting_change_ids": ["b2"]},
        ]
        c["labels"] = F(None, labels, evidence, method="structure_derivation", derivation="Partial net-change classification from source drawings; no atom mapping asserted.")
        c["coverage"] = F(None, "partial", evidence, method="inference", derivation="Hydrogen loss, N2 extrusion, and full atom correspondence are not completely mapped.")
        c["limitations"] = ["No atom-mapped reaction SMILES was asserted.", "C−H cleavage/N2 loss and all bond-order changes are not completely mapped."]
    elif mode == "epoxidation":
        c["bond_changes"] = [
            change("b1", "C", "O", 0, 1, "formed", evidence, "First epoxide C−O bond inferred from substrate/product drawings."),
            change("b2", "C", "O", 0, 1, "formed", evidence, "Second epoxide C−O bond inferred from substrate/product drawings."),
        ]
        c["labels"] = F(None, [{"category_number": 3, "code": "C-O", "supporting_change_ids": ["b1", "b2"]}], evidence, method="structure_derivation", derivation="Partial epoxidation classification from source drawings.")
        c["coverage"] = F(None, "partial", evidence, method="inference", derivation="Alkene bond-order decrease and oxidant coproducts are not fully mapped.")
        c["limitations"] = ["No atom mapping; alkene bond-order decrease and oxidant products are unresolved."]
    elif mode == "transesterification":
        c["bond_changes"] = [change("b1", "C", "O", 0, 1, "formed", evidence, "Alcohol-derived ester C−O connection inferred from the generic scheme.")]
        c["labels"] = F(None, [{"category_number": 3, "code": "C-O", "supporting_change_ids": ["b1"]}], evidence, method="structure_derivation", derivation="Partial transesterification classification.")
        c["coverage"] = F(None, "partial", evidence, method="inference", derivation="Displaced ethoxy group and proton transfers are not mapped.")
        c["limitations"] = ["No atom mapping; cleavage and proton transfer omitted."]
    elif mode == "diazo_transfer":
        c["bond_changes"] = [change("b1", "C", "N", 0, 2, "formed", evidence, "Diazo carbon–nitrogen multiple bond inferred from precursor/product identities.")]
        c["labels"] = F(None, [{"category_number": 16, "code": "DB-C-N", "supporting_change_ids": ["b1"]}], evidence, method="structure_derivation", derivation="Partial diazo-transfer classification.")
        c["coverage"] = F(None, "partial", evidence, method="inference", derivation="Sulfonamide coproduct and proton transfer are not mapped.")
        c["limitations"] = ["No atom mapping; only the principal product bond is classified."]
    elif mode == "none":
        c["bond_changes"] = []
        c["labels"] = F(None, [], evidence, method="inference", derivation="No observed transformation; no expected-product label assigned.")
        c["coverage"] = F(None, "not_applicable", evidence, method="inference", derivation="Negative control/no product.")
        c["limitations"] = []
    else:
        c["bond_changes"] = []
        c["labels"] = F()
        c["coverage"] = F(None, "unresolved", evidence, method="inference", derivation="Insufficient atom correspondence for safe bond classification.")
        c["limitations"] = ["No atom-mapped reaction SMILES; net changes remain unresolved."]


reactions = []


def reaction(suffix, pages, record_type="specific_experiment", categories=None, name=None,
             transform=None, classification="unresolved", unresolved=None):
    evidence = [ev(p) for p in pages]
    rid = f"{DOC}:{suffix}"
    r = st.new_reaction(rid, DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = F(str(categories or []), categories or [], evidence)
    r["location"] = F({"pdf_pages": list(pages)}, {"pdf_pages": list(pages)}, evidence)
    r["reaction_name"] = F(name, name, evidence) if name else F()
    r["transformation_type"] = F(None, transform, evidence, method="inference", derivation="Machine description from the source transformation.") if transform else F()
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
                   "unresolved": list(unresolved or [])}
    r["completeness"] = "partial"
    set_classification(r, classification, evidence)
    reactions.append(r)
    return r


def participant(r, pid, mid, role, stage="reaction", amount=None, mass=None, volume=None,
                equiv=None, molpct=None, ref=None, evidence=None, method="original_extraction"):
    evidence = evidence or r["evidence_ids"]
    p = st.new_participant(pid, mid, role, evidence, stage)
    vals = {"amount": amount, "mass": mass, "volume": volume, "equivalents": equiv, "mol_percent": molpct}
    for key, item in vals.items():
        if item is not None:
            raw = item if isinstance(item, str) else item.get("raw_text", f"{item.get('number')} {item.get('unit')}")
            p["charges"][key] = F(raw, item, evidence, method=method)
    if ref:
        p["equivalence_reference"] = F(ref, ref, evidence, method=method)
    r["participants"].append(p)
    return p


def step(r, sid, stage, pids, operation, temperature=None, time=None, atmosphere=None, evidence=None, method="original_extraction"):
    evidence = evidence or r["evidence_ids"]
    s = st.new_step(sid, evidence, stage)
    s["participant_ids"] = list(pids)
    s["operation"] = F(operation, operation, evidence, method=method)
    if temperature is not None:
        raw = temperature if isinstance(temperature, str) else temperature.get("raw_text", f"{temperature.get('number')} {temperature.get('unit')}")
        s["temperature"] = F(raw, temperature, evidence, method=method)
    if time is not None:
        raw = time if isinstance(time, str) else time.get("raw_text", f"{time.get('number')} {time.get('unit')}")
        s["time"] = F(raw, time, evidence, method=method)
    if atmosphere is not None:
        s["atmosphere"] = F(atmosphere, atmosphere, evidence, method=method)
    r["steps"].append(s)
    return s


def outcome(r, oid, kind, targets, value, basis=None, yield_type=None, scope=None, evidence=None, method="original_extraction"):
    evidence = evidence or r["evidence_ids"]
    o = st.new_outcome(oid, kind, targets, evidence)
    raw = value if isinstance(value, str) else value.get("raw_text", f"{value.get('number')} {value.get('unit')}")
    o["value"] = F(raw, value, evidence, method=method)
    if basis is not None:
        o["measurement_basis"] = F(basis, basis, evidence, method=method)
    if yield_type is not None:
        o["yield_type"] = F(yield_type, yield_type, evidence, method=method)
    if scope is not None:
        o["scope"] = F(scope, scope, evidence, method=method)
    r["outcomes"].append(o)
    return o


def add_note(r, text, evidence=None, method="original_extraction", derivation=None):
    r["notes"].append(F(text if method == "original_extraction" else None, text, evidence or r["evidence_ids"], method=method, derivation=derivation))


# Page 1 prior-work panels: retained as literature examples, not present-paper experiments.
prior_specs = [
    ("scheme1a_sherburn_literature", "Sherburn (2009) dendralene coupling"),
    ("scheme1b_glorius_literature", "Glorius (2013) Rh-catalyzed route"),
    ("scheme1c_wu_literature", "Wu (2017 & 2019) phosphinyl dendralene route"),
    ("scheme1d_tanaka_literature", "Tanaka (2016) Rh-catalyzed allene/alkyne coupling"),
    ("scheme1d_backvall_literature", "Bäckvall (2017) Pd-catalyzed allene oxidation"),
]
for suffix, title in prior_specs:
    r = reaction(suffix, (1,), "literature_example", [], title, classification="unresolved",
                 unresolved=["Generic structures and detailed experimental conditions belong to cited work and were not reconstructed from this review scheme."])
    participant(r, "literature_inputs", M_ALLENE, "reactant")
    participant(r, "literature_product", M_DEND, "product")

# Present-work generic coupling, merging graphical abstract, Scheme 1e, Scheme 2, and page-5 procedure.
GEN = reaction("general_allenic_CH_diazo_coupling", (1, 3, 5), "generic_scheme", [],
               "General Procedure for the Palladium-Catalyzed Oxidative Coupling of the Allenic C−H Bond with α-Diazo Esters",
               "allene C−H oxidative coupling with α-diazo esters", "coupling",
               ["Generic R-group structures are templates; individual product identities are recorded separately."])
participant(GEN, "allene", M_ALLENE, "reactant", amount=q(0.3, "mmol"), equiv=q(1.0, "equiv"), ref=M_ALLENE, evidence=[ev(3), ev(5)])
participant(GEN, "diazo", M_DIAZO, "reactant", amount=q(0.6, "mmol"), equiv=q(2.0, "equiv"), ref=M_ALLENE, evidence=[ev(3), ev(5)])
participant(GEN, "pd", M_PDOAC, "catalyst", mass=q(6.7, "mg"), amount=q(0.03, "mmol"), molpct=q(10, "mol%"), ref=M_ALLENE, evidence=[ev(3), ev(5)])
participant(GEN, "dmbq", M_DMBQ, "reagent", mass=q(44.9, "mg"), amount=q(0.33, "mmol"), equiv=q(1.1, "equiv"), ref=M_ALLENE, evidence=[ev(3), ev(5)])
participant(GEN, "thf", M_THF, "solvent", volume=q(3, "mL"), evidence=[ev(3), ev(5)])
participant(GEN, "product", M_DEND, "product")
participant(GEN, "silica_workup", M_SILICA, "purification_agent", "purification")
participant(GEN, "ea_workup", M_EA, "purification_agent", "purification")
participant(GEN, "pe_workup", M_PE, "purification_agent", "purification")
GEN["reaction_scale"] = F("allene 1 (0.3 mmol)", {"reference_molecule_id": M_ALLENE, "amount": q(0.3, "mmol")}, [ev(3), ev(5)])
for key, raw, norm in [
    ("composition", "dry degassed THF", [M_THF]),
    ("total_volume", "3 mL", q(3, "mL")),
]: GEN["solvent_conditions"][key] = F(raw, norm, [ev(3), ev(5)])
GEN["procedure"]["original_text"] = F("General Procedure for the Palladium-Catalyzed Oxidative Coupling...", "page-5 general procedure", [ev(5)])
GEN["procedure"]["addition_order"] = F("Pd(OAc)2 and 2,6-DMBQ were weighed; flask degassed/backfilled; THF, allene, then diazo ester were added successively.", ["Pd(OAc)2", "2,6-DMBQ", "degassing/backfill", "THF", "allene", "diazo ester"], [ev(5)])
GEN["procedure"]["pretreatment"] = F("10 mL oven-dried flask; evacuated and back-filled with N2 three times; dry degassed THF", "oven-dried flask and three vacuum/N2 cycles", [ev(5)])
step(GEN, "pretreat", "reaction", ["pd", "dmbq"], "weigh into oven-dried flask; evacuate and back-fill with N2 three times", atmosphere="N2", evidence=[ev(5)])
step(GEN, "reaction", "reaction", ["allene", "diazo", "pd", "dmbq", "thf"], "add THF, allene, and α-diazo ester successively; stir in oil bath", q(90, "°C"), q(10, "h"), "N2", [ev(3), ev(5)])
step(GEN, "post_reaction", "workup", ["silica_workup", "ea_workup", "pe_workup"], "source says mixture was allowed to warm to room temperature, then filtered through silica (EA:PE = 1:1); solvent evaporated; residue purified by silica chromatography (EA-PE)", "room temperature", evidence=[ev(5)])
add_note(GEN, "The page-5 phrase 'allowed to warm to room temperature' is preserved verbatim although the preceding hold is 90 °C; no correction to 'cool' was made.", [ev(5)])

# Table 1 optimization entries.
table = [
    (1, M_PDCL2, 5, M_BQ, 1.5, M_THF, 80, 36, 5, 27, None),
    (2, M_PDL, 5, M_BQ, 1.5, M_THF, 80, 35, 4, 26, None),
    (3, M_PDOAC, 5, M_BQ, 1.5, M_THF, 80, 42, 6, 20, None),
    (4, M_PDOAC, 5, M_BQ, 1.5, M_THF, 80, 39, 14, None, (M_DMSO, 1.0)),
    (5, M_PDOAC, 5, M_BQ, 1.5, M_THF, 80, 43, 15, None, (M_SALEN, 1.0)),
    (6, M_PDOAC, 5, M_BQ, 1.1, M_THF, 80, 42, 4, 20, None),
    (7, M_PDOAC, 5, M_BQ, 2.0, M_THF, 80, 37, 3, 38, None),
    (8, M_PDOAC, 5, M_BQ, 2.5, M_THF, 80, 33, 3, 50, None),
    (9, M_PDOAC, 5, M_BQ, 1.1, M_DIOX, 80, 32, "dash", 6, None),
    (10, M_PDOAC, 5, M_BQ, 1.1, M_DME, 80, 39, "dash", 37, None),
    (11, M_PDOAC, 10, M_BQ, 1.1, M_THF, 80, 48, 8, 7, None),
    (12, M_PDOAC, 10, M_BQ, 1.1, M_THF, 90, 59, 3, 6, None),
    (13, M_PDOAC, 10, M_DMBQ, 1.1, M_THF, 90, 73, None, 6, None),
]
for entry, catalyst, loading, oxidant, oxeq, solvent, temp, ymain, ydimer, recovery, additive in table:
    r = reaction(f"table1_entry_{entry:02d}", (2,), "specific_experiment", ["optimization"],
                 f"Table 1 optimization entry {entry}", "allene C−H oxidative coupling", "coupling",
                 ["Workup and purification are not reported for the optimization microscale experiment."])
    participant(r, "allene_1a", M_1A, "reactant", amount=q(0.1, "mmol"), equiv=q(1.0, "equiv"), ref=M_1A)
    participant(r, "diazo_2a", M_2A, "reactant", amount=q(0.2, "mmol"), equiv=q(2.0, "equiv"), ref=M_1A)
    participant(r, "catalyst", catalyst, "catalyst", molpct=q(loading, "mol%"), ref=M_1A)
    participant(r, "oxidant", oxidant, "reagent", equiv=q(oxeq, "equiv"), ref=M_1A)
    participant(r, "solvent", solvent, "solvent", volume=q(1, "mL"))
    participant(r, "product_3aa", M_3AA, "product")
    participant(r, "dimer_3prime", M_3P, "product")
    participant(r, "internal_standard", M_MENO2, "reagent", "characterization")
    if additive:
        participant(r, "additive", additive[0], "additive", equiv=q(additive[1], "equiv"), ref=M_1A)
    r["reaction_scale"] = F("1a (0.1 mmol)", {"reference_molecule_id": M_1A, "amount": q(0.1, "mmol")}, [ev(2)])
    r["solvent_conditions"]["composition"] = F(molecules[solvent]["names"]["original"]["normalized_value"], [solvent], [ev(2)])
    r["solvent_conditions"]["total_volume"] = F("1 mL", q(1, "mL"), [ev(2)])
    step(r, "reaction", "reaction", [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"], "stir under N2", q(temp, "°C"), q(10, "h"), "N2")
    main_basis = "isolated product" if entry == 13 else "1H NMR (400 MHz) using MeNO2 as internal standard"
    outcome(r, "yield_3aa", "yield", [M_3AA], q(ymain, "%"), main_basis, "isolated yield" if entry == 13 else "NMR yield")
    if ydimer == "dash":
        outcome(r, "table_marker_3prime", "other", [M_3P], "table dash (–)", "1H NMR (400 MHz) using MeNO2 as internal standard", scope={"metric": "reported dimer result", "source_marker": "–"})
    elif ydimer is not None:
        outcome(r, "yield_3prime", "yield", [M_3P], q(ydimer, "%"), "1H NMR (400 MHz) using MeNO2 as internal standard", "NMR yield")
    if recovery is not None:
        outcome(r, "recovered_1a", "recovery", [M_1A], q(recovery, "%"), "1H NMR (400 MHz) using MeNO2 as internal standard", scope={"metric": "recovered starting material"})
    if ydimer is None:
        add_note(r, "The dimer 3′ cell is blank for this entry and remains unreported.")

# Scheme 2 concrete scope molecules and records.
scope_rows = [
    ("3aa", "1a", "2a", 76, 71.7, 5), ("3ab", "1a", "2b", 71, 82.8, 5),
    ("3ac", "1a", "2c", 87, 89.7, 5), ("3ad", "1a", "2d", 87, 89.7, 5),
    ("3ae", "1a", "2e", 83, 94.9, 6), ("3af", "1a", "2f", 81, 95.2, 6),
    ("3ag", "1a", "2g", 74, 80.4, 6), ("3ah", "1a", "2h", 67, 61.0, 6),
    ("3ai", "1a", "2i", 81, 61.5, 6), ("3aj", "1a", "2j", 66, 55.8, 6),
    ("3ak", "1a", "2k", 70, 78.4, 6), ("3al", "1a", "2l", 68, 57.7, 6),
    ("3ba", "1b", "2a", 85, 86.6, 6), ("3ca", "1c", "2a", 82, 87.6, 6),
    ("3da", "1d", "2a", 62, 68.3, 6), ("3ea", "1e", "2a", 60, 64.5, 6),
    ("3fa", "1f", "2a", 35, None, None), ("3ga", "1g", "2a", 67, None, None),
    ("3gi", "1g", "2i", 54, None, None), ("3ha", "1h", "2a", 63, None, None),
    ("3ia", "1i", "2a", 75, None, None), ("3ib", "1i", "2b", 68, None, None),
    ("3ie", "1i", "2e", 74, None, None), ("3ii", "1i", "2i", 71, None, None),
    ("3ja", "1j", "2a", 80, None, None), ("3ik", "1i", "2k", 69, None, None),
    ("3im", "1i", "2m", 65, None, None), ("3in", "1i", "2n", 36, None, None),
]
allene_mids = {"1a": M_1A}
for lab in ["1b", "1c", "1d", "1e", "1f", "1g", "1h", "1i", "1j", "1k"]:
    allene_mids[lab] = mol(lab, f"allene {lab}", (3, 4) if lab in {"1e", "1g", "1h", "1i"} else (3,), chirality="unspecified", geometry="not_applicable", note="Structure shown in Scheme 2; full candidate SMILES not encoded in this machine extraction.")
diazo_mids = {"2a": M_2A}
for lab in ["2b", "2c", "2d", "2e", "2f", "2g", "2h", "2i", "2j", "2k", "2l", "2m", "2n"]:
    diazo_mids[lab] = mol(lab, f"α-diazo ester {lab}", (3, 4, 5), chirality="unresolved", geometry="not_applicable", note="Compound identity follows the source label/drawing; complete candidate SMILES is not encoded.")
product_mids = {"3aa": M_3AA}
for prod, _, _, _, _, _ in scope_rows:
    if prod == "3aa": continue
    if prod == "3ja":
        product_mids[prod] = mol(prod, "6-benzyl 1-ethyl 4-ethylidene-5-methylenehex-2-enedioate (E/Z mixture)", (3, 7), chirality="unresolved", geometry="mixture", mixture={"metric": "E/Z", "ratio": "2:1", "configurations": "major/minor not independently assigned beyond source order"})
    elif prod == "3in":
        product_mids[prod] = mol(prod, "ethyl 3-cyclohexylidene-6-(dimethylamino)-2-ethylidene-6-oxohex-4-enoate (Z/E mixture)", (3, 7), chirality="unresolved", geometry="mixture", mixture={"metric": "Z/E", "ratio": "2.2/1", "configurations": "source-reported order"})
    else:
        pgs = (3, 5) if prod in {"3ab", "3ac", "3ad"} else (3, 6) if prod in {"3ae", "3af", "3ag", "3ah", "3ai", "3aj", "3ak", "3al", "3ba", "3ca", "3da", "3ea"} else (3, 7)
        product_mids[prod] = mol(prod, f"[3]dendralene product {prod}", pgs, chirality="unresolved", geometry="specified", double_bonds=[{"assignment": "E", "basis": "Scheme 2 drawing/product label"}], note="Full structure is visible in Scheme 2 but remains unencoded pending expert structure review.")


def add_standard_coupling_participants(r, allene_mid, diazo_mid, product_mid=None):
    inherited = [ev(3), ev(5)]
    participant(r, "allene", allene_mid, "reactant", amount=q(0.3, "mmol"), equiv=q(1, "equiv"), ref=allene_mid, evidence=inherited, method="inherited_general_condition")
    participant(r, "diazo", diazo_mid, "reactant", amount=q(0.6, "mmol"), equiv=q(2, "equiv"), ref=allene_mid, evidence=inherited, method="inherited_general_condition")
    participant(r, "pd", M_PDOAC, "catalyst", mass=q(6.7, "mg"), amount=q(0.03, "mmol"), molpct=q(10, "mol%"), ref=allene_mid, evidence=inherited, method="inherited_general_condition")
    participant(r, "dmbq", M_DMBQ, "reagent", mass=q(44.9, "mg"), amount=q(0.33, "mmol"), equiv=q(1.1, "equiv"), ref=allene_mid, evidence=inherited, method="inherited_general_condition")
    participant(r, "thf", M_THF, "solvent", volume=q(3, "mL"), evidence=inherited, method="inherited_general_condition")
    if product_mid: participant(r, "product", product_mid, "product")
    participant(r, "silica", M_SILICA, "purification_agent", "purification", evidence=[ev(5)], method="inherited_general_condition")
    participant(r, "ea", M_EA, "purification_agent", "purification", evidence=[ev(5)], method="inherited_general_condition")
    participant(r, "pe", M_PE, "purification_agent", "purification", evidence=[ev(5)], method="inherited_general_condition")
    r["reaction_scale"] = F("allene 1 (0.3 mmol)", {"reference_molecule_id": allene_mid, "amount": q(0.3, "mmol")}, inherited, method="inherited_general_condition")
    r["solvent_conditions"]["composition"] = F("dry degassed THF", [M_THF], inherited, method="inherited_general_condition")
    r["solvent_conditions"]["total_volume"] = F("3 mL", q(3, "mL"), inherited, method="inherited_general_condition")
    r["procedure"]["general_procedure_reference"] = F("General Procedure for the Palladium-Catalyzed Oxidative Coupling", GEN["reaction_id"], [ev(5)], method="inherited_general_condition")
    r["procedure"]["addition_order"] = F("THF, allene, then α-diazo ester", ["THF", "allene", "α-diazo ester"], [ev(5)], method="inherited_general_condition")
    step(r, "reaction", "reaction", ["allene", "diazo", "pd", "dmbq", "thf"], "stir in oil bath", q(90, "°C"), q(10, "h"), "N2", inherited, "inherited_general_condition")
    step(r, "workup", "workup", ["silica", "ea", "pe"], "silica plug (EA:PE = 1:1), evaporation, then silica-gel chromatography (EA-PE)", "room temperature", evidence=[ev(5)], method="inherited_general_condition")

for prod, allene, diazo, yld, mass, mass_page in scope_rows:
    pages = [3]
    if mass_page: pages.append(mass_page)
    r = reaction(f"scheme2_scope_{prod}", tuple(pages), "specific_experiment", ["substrate_scope"],
                 f"Scheme 2 synthesis of {prod}", "allene C−H oxidative coupling", "coupling",
                 ["Complete structure/atom mapping remains pending expert review."])
    r["general_scheme_ids"] = [GEN["reaction_id"]]
    add_standard_coupling_participants(r, allene_mids[allene], diazo_mids[diazo], product_mids[prod])
    outcome(r, f"yield_{prod}", "yield", [product_mids[prod]], q(yld, "%"), "isolated product", "isolated yield", evidence=[ev(3)])
    if mass is not None:
        outcome(r, f"amount_{prod}", "product_amount", [product_mids[prod]], q(mass, "mg"), "isolated product", evidence=[ev(mass_page)])
    if prod == "3ja":
        outcome(r, "selectivity_EZ", "selectivity", [product_mids[prod]], "E/Z = 2:1", None, scope={"metric": "E/Z ratio", "order": ["E", "Z"], "configuration_assignment": "source order; analytical method unreported"}, evidence=[ev(3)])
        r["review"]["unresolved"].append("E/Z-ratio measurement method is not reported on page 3.")
    if prod == "3in":
        outcome(r, "selectivity_ZE", "selectivity", [product_mids[prod]], "Z/E = 2.2/1", None, scope={"metric": "Z/E ratio", "order": ["Z", "E"], "configuration_assignment": "source order; analytical method unreported"}, evidence=[ev(3)])
        r["review"]["unresolved"].append("Z/E-ratio measurement method is not reported on page 3.")

# Explicit no-product scope branch.
NEG = reaction("scheme2_scope_1k_no_product", (3,), "specific_experiment", ["substrate_scope", "control"],
               "Scheme 2 reaction of monosubstituted allene 1k", "attempted allene C−H oxidative coupling", "none",
               ["No product identity or analytical detection method is reported; none was invented."])
NEG["general_scheme_ids"] = [GEN["reaction_id"]]
add_standard_coupling_participants(NEG, allene_mids["1k"], M_2A, None)
outcome(NEG, "no_product", "no_reaction", [], "no product observed", None, scope={"reported statement": "3ka, not observed"}, evidence=[ev(3)])

# Scale-up, epoxidation, and isotope experiments.
SCALE = reaction("eq1_scale_up_3aa", (2, 3, 5), "specific_experiment", ["scale_up"],
                 "Procedure for the Scale-up Experiment", "allene C−H oxidative coupling", "coupling",
                 ["The page-5 procedure gives allene 1a amount but no mass; diazo 2a amount is given only as 2.0 equivalents."])
participant(SCALE, "allene", M_1A, "reactant", amount=q(3.0, "mmol"), equiv=q(1.0, "equiv"), ref=M_1A)
participant(SCALE, "diazo", M_2A, "reactant", equiv=q(2.0, "equiv"), ref=M_1A)
participant(SCALE, "pd", M_PDOAC, "catalyst", mass=q(67.4, "mg"), amount=q(0.3, "mmol"), molpct=q(10, "mol%"), ref=M_1A)
participant(SCALE, "dmbq", M_DMBQ, "reagent", mass=q(449.3, "mg"), amount=q(3.3, "mmol"), equiv=q(1.1, "equiv"), ref=M_1A)
participant(SCALE, "thf", M_THF, "solvent", volume=q(30, "mL"))
participant(SCALE, "product", M_3AA, "product")
participant(SCALE, "silica", M_SILICA, "purification_agent", "purification")
participant(SCALE, "ea", M_EA, "purification_agent", "purification")
participant(SCALE, "pe", M_PE, "purification_agent", "purification")
SCALE["reaction_scale"] = F("1a (3.0 mmol)", {"reference_molecule_id": M_1A, "amount": q(3.0, "mmol")}, [ev(2), ev(5)])
SCALE["solvent_conditions"]["composition"] = F("dry degassed THF: 20 mL initially plus 10 mL substrate solution", [M_THF], [ev(5)])
SCALE["solvent_conditions"]["total_volume"] = F("30 mL THF", q(30, "mL"), [ev(5)])
SCALE["procedure"]["addition_order"] = F("Pd(OAc)2 and DMBQ; degas/backfill; 20 mL THF; solution of 1a and 2a in 10 mL THF", ["Pd(OAc)2", "DMBQ", "20 mL THF", "substrate solution in 10 mL THF"], [ev(5)])
step(SCALE, "pretreat", "reaction", ["pd", "dmbq"], "weigh into 100 mL oven-dried flask; evacuate/back-fill with N2 three times; add 20 mL THF", atmosphere="N2", evidence=[ev(5)])
step(SCALE, "reaction", "reaction", ["allene", "diazo", "pd", "dmbq", "thf"], "add substrate solution and stir in oil bath", q(90, "°C"), q(10, "h"), "N2", [ev(5)])
step(SCALE, "workup", "workup", ["silica", "ea", "pe"], "silica plug (EA:PE 1:1), evaporation, silica chromatography (EA:PE 1:20)", "room temperature", evidence=[ev(5)])
outcome(SCALE, "yield_3aa", "yield", [M_3AA], q(79, "%"), "isolated product", "isolated yield", evidence=[ev(2), ev(3), ev(5)])
outcome(SCALE, "amount_3aa", "product_amount", [M_3AA], q(0.748, "g"), "isolated product", evidence=[ev(2), ev(3), ev(5)])

EPOX = reaction("eq2_epoxidation_3aa_to_4", (3,), "specific_experiment", ["derivatization"],
                "Transformation of [3]dendralene 3aa with m-CPBA", "site-selective epoxidation", "epoxidation",
                ["Detailed workup appears on page 8, outside this batch; it was not claimed as page-1–6 coverage."])
participant(EPOX, "substrate", M_3AA, "reactant", amount=q(0.2, "mmol"), equiv=q(1, "equiv"), ref=M_3AA)
participant(EPOX, "mcpba", M_MCPBA, "reagent", equiv=q(1.8, "equiv"), ref=M_3AA)
participant(EPOX, "dcm", M_DCM, "solvent")
participant(EPOX, "product", M_4, "product")
EPOX["reaction_scale"] = F("3aa, 0.2 mmol", {"reference_molecule_id": M_3AA, "amount": q(0.2, "mmol")}, [ev(3)])
EPOX["solvent_conditions"]["composition"] = F("DCM", [M_DCM], [ev(3)])
step(EPOX, "cold_stage", "reaction", ["substrate", "mcpba", "dcm"], "treat 3aa with m-CPBA in DCM", q(0, "°C"), evidence=[ev(3)])
step(EPOX, "warm_stage", "reaction", ["substrate", "mcpba", "dcm"], "warm from 0 °C to room temperature; ramp time and hold time not stated on page 3", "room temperature", evidence=[ev(3)])
outcome(EPOX, "yield_4", "yield", [M_4], q(50, "%"), "reported product yield; isolation basis not stated on page 3", "reported yield", evidence=[ev(3)])

for suffix, additive, yld in [("eq3_deuterium_no_acoh", False, 78), ("eq4_deuterium_with_acoh", True, 8)]:
    r = reaction(suffix, (3,), "specific_experiment", ["mechanistic_study"],
                 "Deuterium-labeling experiment", "allene C−H oxidative coupling", "coupling",
                 ["Page 3 does not report full product-isolation details for this mechanistic experiment."])
    participant(r, "allene_d6", M_1AD6, "reactant", equiv=q(1, "equiv"), ref=M_1AD6)
    participant(r, "diazo", M_2A, "reactant")
    participant(r, "pd", M_PDOAC, "catalyst", molpct=q(10, "mol%"), ref=M_1AD6)
    participant(r, "dmbq", M_DMBQ, "reagent", equiv=q(1.1, "equiv"), ref=M_1AD6)
    participant(r, "thf", M_THF, "solvent")
    if additive: participant(r, "acoh", M_ACOH, "additive", equiv=q(1.0, "equiv"), ref=M_1AD6)
    participant(r, "product", M_3AAD6, "product")
    step(r, "reaction", "reaction", [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"], "standard conditions", q(90, "°C"), q(10, "h"), "N2", evidence=[ev(3)], method="inherited_general_condition")
    outcome(r, "yield_d6_product", "yield", [M_3AAD6], q(yld, "%"), "reported product yield; detailed basis not stated on page 3", "reported yield")
    outcome(r, "deuterium_scrambling", "other", [M_3AAD6], "without scrambling of deuterium content" if not additive else "product shown with retained D6 labeling", "source statement/drawing", scope={"metric": "deuterium retention"})

MECH = reaction("scheme3_proposed_mechanism", (3,), "proposed_mechanism", [], "Scheme 3. Proposed Reaction Mechanism", "proposed Pd catalytic cycle", "unresolved",
                ["Intermediates A–E are proposed and were not treated as isolated experiments."])
participant(MECH, "allene", M_ALLENE, "reactant")
participant(MECH, "diazo", M_DIAZO, "reactant")
participant(MECH, "pd", M_PDOAC, "catalyst")
participant(MECH, "oxidant", M_DMBQ, "reagent")
participant(MECH, "product", M_DEND, "product")
add_note(MECH, "Author-proposed sequence: allene coordination; selective allenic C−H cleavage to π-allylpalladium B; vinyl-Pd C; palladium carbene D; migratory insertion E; β-H elimination; oxidation regenerates Pd(II).", [ev(3)])

# Page 4 allene preparations.
ALLENE_GP = reaction("general_preparation_allenes_1g_1i", (4,), "generic_scheme", [], "General procedure for the preparation of allenes 1g−1i", "propargyl alcohol to allene rearrangement", "unresolved",
                      ["Generic substituted structures are not expanded to a single molecular graph."])
participant(ALLENE_GP, "propargyl", M_PROPARGYL, "reactant", amount=q(8, "mmol"), equiv=q(1, "equiv"), ref=M_PROPARGYL)
participant(ALLENE_GP, "aminal", M_AMINAL, "reagent", amount=q(40, "mmol"), equiv=q(5, "equiv"), ref=M_PROPARGYL)
participant(ALLENE_GP, "toluene", M_TOL, "solvent", volume=q(20, "mL"))
participant(ALLENE_GP, "product", M_ALLENE, "product")
participant(ALLENE_GP, "ea", M_EA, "workup_agent", "workup", volume=q(10, "mL"))
participant(ALLENE_GP, "hcl", M_HCL, "workup_agent", "workup", volume=q(20, "mL"))
participant(ALLENE_GP, "ea_extract", M_EA, "workup_agent", "workup", volume=q(40, "mL"))
participant(ALLENE_GP, "brine", M_BRINE, "workup_agent", "workup")
participant(ALLENE_GP, "na2so4", M_NA2SO4, "workup_agent", "workup")
participant(ALLENE_GP, "silica", M_SILICA, "purification_agent", "purification")
participant(ALLENE_GP, "pe", M_PE, "purification_agent", "purification")
step(ALLENE_GP, "reaction", "reaction", ["propargyl", "aminal", "toluene", "product"], "dissolve in dry toluene and stir in an oil bath", q(110, "°C"), q(5, "h"), "N2")
step(ALLENE_GP, "cool_quench", "workup", ["ea", "hcl"], "cool to room temperature; add ethyl acetate and 4 M HCl", "room temperature")
step(ALLENE_GP, "extract", "workup", ["ea_extract", "brine", "na2so4"], "extract with ethyl acetate (2 × 20 mL), wash combined organics with brine, dry over Na2SO4, evaporate")
step(ALLENE_GP, "purify", "purification", ["silica", "pe", "ea"], "silica-gel chromatography, PE:EA = 2:1")

prep_allenes = [
    ("1e", 10, 227, 12, False), ("1g", 8, 129, 11, True),
    ("1h", 7, 111, 9, True), ("1i", 8, 241, 16, True),
]
for lab, scale, mass, yld, use_gp in prep_allenes:
    pages = (4,)
    precursor = mol(f"precursor_{lab}", f"source/literature precursor for {lab}", pages, chirality="unresolved", geometry="unresolved", note="Specific precursor identity is not supplied in the product characterization block.")
    r = reaction(f"preparation_{lab}", pages, "specific_experiment", ["other"], f"Preparation of allene {lab}", "allene preparation", "unresolved",
                 ["Specific precursor structure is unresolved." if not use_gp else "Product structure remains unencoded pending expert review."])
    if use_gp:
        r["general_scheme_ids"] = [ALLENE_GP["reaction_id"]]
        participant(r, "precursor", precursor, "reactant", amount=q(scale, "mmol"), equiv=q(1, "equiv"), ref=precursor)
        participant(r, "aminal", M_AMINAL, "reagent", equiv=q(5, "equiv"), ref=precursor, method="inherited_general_condition")
        participant(r, "toluene", M_TOL, "solvent", volume=q(20, "mL"), method="inherited_general_condition")
        participant(r, "product", allene_mids[lab], "product")
        step(r, "reaction", "reaction", ["precursor", "aminal", "toluene", "product"], "general allene procedure", q(110, "°C"), q(5, "h"), "N2", method="inherited_general_condition")
    else:
        participant(r, "precursor", precursor, "reactant")
        participant(r, "product", allene_mids[lab], "product")
        step(r, "reported_literature_method", "reaction", ["precursor", "product"], "synthesized according to cited literature; conditions not reproduced on page 4")
    r["reaction_scale"] = F(f"{scale} mmol scale", {"reference_molecule_id": precursor, "amount": q(scale, "mmol")}, [ev(4)])
    outcome(r, "yield", "yield", [allene_mids[lab]], q(yld, "%"), "reported isolated material; exact measurement basis not explicitly stated", "reported yield")
    outcome(r, "product_amount", "product_amount", [allene_mids[lab]], q(mass, "mg"), "isolated material")

# Three-step deuterated allene synthesis as one reported overall experiment.
D6 = reaction("three_step_synthesis_1a_d6", (4,), "specific_experiment", ["mechanistic_study", "other"], "Synthesis of allene 1a-d6", "multistep synthesis of deuterated allene", "unresolved",
              ["The 31% is an overall yield over three steps; it is not assigned to the final Johnson–Claisen step alone.", "The amount of d6-acetone and triethyl orthoacetate is not stated on page 4."])
participant(D6, "tmsacetylene", M_TMSAC, "reactant", amount=q(20, "mmol"), equiv=q(1, "equiv"), ref=M_TMSAC)
participant(D6, "nbuli", M_NBULI, "reagent", equiv=q(1.2, "equiv"), ref=M_TMSAC)
participant(D6, "d6acetone", M_D6AC, "reagent")
participant(D6, "k2co3", M_K2CO3, "reagent", molpct=q(25, "mol%"), ref=M_TMSAC)
participant(D6, "meoh", M_MEOH, "solvent")
participant(D6, "intermediate", M_D5, "product")
participant(D6, "etco2h", M_ETCO2H, "catalyst", molpct=q(15, "mol%"), ref=M_D5)
participant(D6, "orthoacetate", M_ORTHO, "reagent")
participant(D6, "product", M_1AD6, "product")
D6["reaction_scale"] = F("20 mmol ethynyltrimethylsilane", {"reference_molecule_id": M_TMSAC, "amount": q(20, "mmol")}, [ev(4)])
step(D6, "lithiation", "reaction", ["tmsacetylene", "nbuli"], "add 1.2 equiv n-BuLi in THF", q(-78, "°C"))
step(D6, "d6_acetone", "reaction", ["d6acetone"], "add d6-acetone and hold", "room temperature", q(1, "h"))
step(D6, "desilylation", "reaction", ["k2co3", "meoh", "intermediate"], "treat with K2CO3 in MeOH", "room temperature", q(2.5, "h"))
step(D6, "johnson_claisen", "reaction", ["intermediate", "etco2h", "orthoacetate", "product"], "use [d6]-5 without further purification; react with triethyl orthoacetate under EtCO2H catalysis in an oil bath", q(140, "°C"), q(4, "h"))
outcome(D6, "overall_yield", "yield", [M_1AD6], q(31, "%"), "reported overall yield in three steps", "three-step overall yield", scope={"start": M_TMSAC, "end": M_1AD6})
outcome(D6, "product_amount", "product_amount", [M_1AD6], q(986, "mg"), "isolated product")

# Generic two-stage preparation of diazo esters.
DIAZO_GP = reaction("general_preparation_diazo_2b_2e_2g_2m", (4,), "generic_scheme", [], "Procedure for the Preparation of α-Diazo Esters 2b, 2c, 2d, 2e, 2g, and 2m", "transesterification followed by diazo transfer", "unresolved",
                    ["This generic record spans two operations; product-specific yields are stored separately."])
participant(DIAZO_GP, "beta", M_BETA, "reactant", amount=q(10, "mmol"), equiv=q(1, "equiv"), ref=M_BETA)
participant(DIAZO_GP, "alcohol", M_ALCOHOL, "reactant", amount=q(13, "mmol"), equiv=q(1.3, "equiv"), ref=M_BETA)
participant(DIAZO_GP, "toluene", M_TOL, "solvent", volume=q(3, "mL"))
participant(DIAZO_GP, "intermediate", M_2PRIME, "product")
participant(DIAZO_GP, "pabsa", M_PABSA, "reagent", equiv=q(1.5, "equiv"), ref=M_2PRIME)
participant(DIAZO_GP, "dbu", M_DBU, "reagent", equiv=q(1.5, "equiv"), ref=M_2PRIME)
participant(DIAZO_GP, "mecn", M_MECN, "solvent")
participant(DIAZO_GP, "product", M_DIAZO, "product")
participant(DIAZO_GP, "nh4cl", M_NH4CL, "workup_agent", "workup", volume=q(20, "mL"))
participant(DIAZO_GP, "dcm", M_DCM, "workup_agent", "workup", volume=q(60, "mL"))
participant(DIAZO_GP, "na2so4", M_NA2SO4, "workup_agent", "workup")
participant(DIAZO_GP, "silica", M_SILICA, "purification_agent", "purification")
participant(DIAZO_GP, "ea", M_EA, "purification_agent", "purification")
participant(DIAZO_GP, "pe", M_PE, "purification_agent", "purification")
step(DIAZO_GP, "transesterification", "reaction", ["beta", "alcohol", "toluene", "intermediate"], "heat with distillation of the toluene/ethanol azeotrope until no azeotrope distills", rng(160, 180, "°C", "160−180 °C"), rng(2, 4, "h", "2−4 h"))
step(DIAZO_GP, "intermediate_purification", "purification", ["silica", "ea", "pe"], "flash silica filtration (1:1 EtOAc/petroleum ether), evaporation, then silica chromatography (1:10 EtOAc:petroleum ether)")
step(DIAZO_GP, "dbu_addition", "reaction", ["intermediate", "pabsa", "dbu", "mecn"], "slowly add DBU to 2′/p-ABSA in MeCN (0.3 M)", q(0, "°C"), atmosphere="N2")
step(DIAZO_GP, "diazo_transfer", "reaction", ["intermediate", "pabsa", "dbu", "mecn", "product"], "allow to warm to room temperature and stir overnight", "room temperature", "overnight", "N2")
step(DIAZO_GP, "workup", "workup", ["nh4cl", "dcm", "na2so4"], "quench with saturated NH4Cl; extract with DCM (3 × 20 mL); dry over Na2SO4; filter and concentrate")
step(DIAZO_GP, "final_purification", "purification", ["silica", "ea", "pe"], "silica chromatography (1:20 EtOAc:petroleum ether)")

intermediate_data = [("2bprime", "[1,1′-biphenyl]-4-ylmethyl 2-methyl-3-oxobutanoate", 25, 707, 10),
                     ("2cprime", "4-methoxybenzyl 2-methyl-3-oxobutanoate", 46, None, 10),
                     ("2dprime", "4-chlorobenzyl 2-methyl-3-oxobutanoate", 36, 2144, 25),
                     ("2eprime", "4-(trifluoromethyl)benzyl 2-methyl-3-oxobutanoate", 69, 2826, 15)]
precursors = {}
for lab, name, yld, mass, scale in intermediate_data:
    mid = mol(lab, name, (4, 5), chirality="unspecified", geometry="not_applicable")
    precursors[lab] = mid
    alcohol = mol(f"alcohol_{lab}", f"alcohol precursor for {lab}", (4,), chirality="unresolved", geometry="not_applicable")
    r = reaction(f"preparation_{lab}", (4, 5), "specific_experiment", ["other"], f"Preparation of {lab}", "transesterification", "transesterification",
                 ["Specific alcohol structure is retained by product identity but not encoded as SMILES."])
    r["general_scheme_ids"] = [DIAZO_GP["reaction_id"]]
    participant(r, "beta", M_BETA, "reactant", amount=q(scale, "mmol"), equiv=q(1, "equiv"), ref=M_BETA)
    participant(r, "alcohol", alcohol, "reactant", equiv=q(1.3, "equiv"), ref=M_BETA, method="inherited_general_condition")
    participant(r, "toluene", M_TOL, "solvent", volume=q(3, "mL"), method="inherited_general_condition")
    participant(r, "product", mid, "product")
    step(r, "reaction", "reaction", ["beta", "alcohol", "toluene", "product"], "transesterification with azeotropic distillation", rng(160, 180, "°C", "160−180 °C"), rng(2, 4, "h", "2−4 h"), evidence=[ev(4)], method="inherited_general_condition")
    outcome(r, "yield", "yield", [mid], q(yld, "%"), "reported isolated material; exact basis not explicitly stated", "reported yield")
    if mass is not None: outcome(r, "product_amount", "product_amount", [mid], q(mass, "mg"), "isolated material")

# Precursors for 2g/2m are identified only by labels and cited literature.
for lab in ["2gprime", "2mprime"]:
    precursors[lab] = mol(lab, f"precursor {lab}", (5,), chirality="unresolved", geometry="not_applicable", note="Prepared according to cited literature; no preparation yield is assigned here.")

diazo_final = [
    ("2b", "2bprime", 44, 601), ("2c", "2cprime", 35, 350),
    ("2d", "2dprime", 39, 786), ("2e", "2eprime", 19, 500),
    ("2g", "2gprime", 40, 592), ("2m", "2mprime", 43, 730),
]
for lab, precursor_lab, yld, mass in diazo_final:
    r = reaction(f"preparation_{lab}", (4, 5), "specific_experiment", ["other"], f"Preparation of α-diazo ester {lab}", "diazo transfer", "diazo_transfer",
                 ["Precursor structure remains unencoded pending expert review."])
    r["general_scheme_ids"] = [DIAZO_GP["reaction_id"]]
    participant(r, "precursor", precursors[precursor_lab], "reactant", equiv=q(1, "equiv"), ref=precursors[precursor_lab])
    participant(r, "pabsa", M_PABSA, "reagent", equiv=q(1.5, "equiv"), ref=precursors[precursor_lab], method="inherited_general_condition")
    participant(r, "dbu", M_DBU, "reagent", equiv=q(1.5, "equiv"), ref=precursors[precursor_lab], method="inherited_general_condition")
    participant(r, "mecn", M_MECN, "solvent", method="inherited_general_condition")
    participant(r, "product", diazo_mids[lab], "product")
    step(r, "cold_addition", "reaction", ["precursor", "pabsa", "dbu", "mecn"], "slow DBU addition", q(0, "°C"), atmosphere="N2", evidence=[ev(4)], method="inherited_general_condition")
    step(r, "room_temperature_hold", "reaction", ["precursor", "pabsa", "dbu", "mecn", "product"], "allow to warm and stir", "room temperature", "overnight", "N2", [ev(4)], "inherited_general_condition")
    outcome(r, "yield", "yield", [diazo_mids[lab]], q(yld, "%"), "reported isolated material; exact basis not explicitly stated", "reported yield")
    outcome(r, "product_amount", "product_amount", [diazo_mids[lab]], q(mass, "mg"), "isolated material")

# 2h and 2k have product-specific yields but only cited-literature preparation routes on pages 4--5.
for lab, yld, mass in [("2h", 17, 246), ("2k", 47, 932)]:
    precursor = mol(f"literature_precursor_{lab}", f"literature precursor(s) for {lab}", (4, 5), chirality="unresolved", geometry="unresolved")
    r = reaction(f"preparation_{lab}_literature_route", (4, 5), "specific_experiment", ["other"], f"Preparation of α-diazo ester {lab}", "literature-route α-diazo ester preparation", "unresolved",
                 ["The cited literature route is summarized but complete conditions are not reproduced on pages 4–5."])
    participant(r, "precursor", precursor, "reactant")
    participant(r, "product", diazo_mids[lab], "product")
    step(r, "reported_route", "reaction", ["precursor", "product"], "prepared following the cited literature procedure; page-specific conditions unreported")
    outcome(r, "yield", "yield", [diazo_mids[lab]], q(yld, "%"), "reported isolated material; exact basis not explicitly stated", "reported yield")
    outcome(r, "product_amount", "product_amount", [diazo_mids[lab]], q(mass, "mg"), "isolated material")

bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": f"{DOC}:page1_background", "document_id": DOC, "evidence_ids": [ev(1)], "kind": "background_and_prior_art", "note": "Narrative review and Scheme 1 prior work; literature arrows are separately classified as literature_example."},
    {"context_id": f"{DOC}:pages4_6_characterization", "document_id": DOC, "evidence_ids": [ev(4), ev(5), ev(6)], "kind": "characterization", "note": "NMR, IR, HRMS, Rf, color/state blocks support compound identity and isolated amount but are not separate synthesis experiments."},
    {"context_id": f"{DOC}:si_notice", "document_id": DOC, "evidence_ids": [ev(1)], "kind": "source_availability", "note": "The article indicates Supporting Information exists, but no SI file was supplied or associated in this workspace."},
]

SUBMITTED.write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")

# Raw page inventory and exact-SMILES stereo diagnostics.
diagnostics = {}
for key in [M_1A, M_2A, M_3AA, M_PDOAC, M_PDCL2, M_BQ, M_DMBQ, M_THF, M_DIOX, M_DME, M_MENO2]:
    smi = molecules[key]["smiles"]["normalized_value"]
    if smi: diagnostics[key] = st.stereo_inspect(smi)
(RAW / "stereo_diagnostics_pages_0001_0006.json").write_text(json.dumps(diagnostics, ensure_ascii=False, indent=2) + "\n")

page_records = {p: [] for p in range(1, 7)}
for r in reactions:
    for p in range(1, 7):
        if ev(p) in r["evidence_ids"]:
            page_records[p].append(r["reaction_id"])

inventory = {
    "document_id": DOC,
    "batch_id": "allene-stereo-v31-200-20261017",
    "page_range": [1, 6],
    "source_experiments_before_batch": 0,
    "reaction_count": len(reactions),
    "molecule_count": len(molecules),
    "decisions": {
        "scope_successes": 28,
        "scope_negative_branches": 1,
        "optimization_rows": 13,
        "same_experiment_deduplication": "Scheme 2 records merge page-3 yields with page-5/6 product amounts; scale-up page-2 prose, page-3 eq 1, and page-5 detailed procedure are one record.",
        "si_association": "unestablished",
        "spectral_blocks": "source_context, not reactions",
    },
    "page_record_ids": page_records,
}
(RAW / "pages_0001_0006_inventory.json").write_text(json.dumps(inventory, ensure_ascii=False, indent=2) + "\n")

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title/abstract/graphical abstract", "finding": "Article title and present-work generic allene C−H/α-diazo ester coupling; graphical abstract merged with the same generic record.", "evidence_ids": [ev(1), ev(1, True)]},
            {"region": "Scheme 1a–d", "finding": "Five prior-work transformations retained as literature examples, not experiments of this paper.", "evidence_ids": [ev(1), ev(1, True)]},
            {"region": "Scheme 1e", "finding": "Present-work generic transformation; merged with Scheme 2/general procedure.", "evidence_ids": [ev(1), ev(1, True)]}],
         "record_ids": page_records[1], "unresolved": ["Scheme 1 R-group structures are generic; cited-work experimental details were not reconstructed."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Table 1 rows 1–13 and footnotes a–g", "finding": "Thirteen independent optimization experiments, including additives, catalyst-loading changes, two table dashes, isolated entry 13, and recovered 1a outcomes.", "evidence_ids": [ev(2), ev(2, True)]},
            {"region": "right-column scope/scale-up prose", "finding": "Narrative scope summary and scale-up result linked to Scheme 2 and eq 1 records; not duplicated.", "evidence_ids": [ev(2), ev(2, True)]}],
         "record_ids": page_records[2], "unresolved": ["Blank Table 1 cells remain unreported; table dashes were not converted to numeric zero."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Scheme 2 product grid and footnote", "finding": "Twenty-eight successful substrate-scope experiments plus one explicit no-product branch; all share the stated standard conditions and isolated-yield basis.", "evidence_ids": [ev(3), ev(3, True)]},
            {"region": "eq 1 and eq 2", "finding": "Distinct scale-up synthesis of 3aa and derivatization/epoxidation to 4.", "evidence_ids": [ev(3), ev(3, True)]},
            {"region": "eq 3 and eq 4", "finding": "Two deuterium mechanistic experiments, without and with 1.0 equiv AcOH; 78% and 8% product yields.", "evidence_ids": [ev(3), ev(3, True)]},
            {"region": "Scheme 3 and mechanism prose", "finding": "Author-proposed catalytic cycle, stored as proposed_mechanism rather than isolated experiments.", "evidence_ids": [ev(3), ev(3, True)]}],
         "record_ids": page_records[3], "unresolved": ["E/Z and Z/E measurement methods for 3ja/3in are unreported on this page; several complete product SMILES remain pending expert structure review."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "General Methods and allene 1g–1i procedure", "finding": "General analytical text retained as context; generic allene procedure and product-specific preparations 1e/1g/1h/1i recorded.", "evidence_ids": [ev(4), ev(4, True)]},
            {"region": "Synthesis of allene 1a-d6 scheme/text", "finding": "Three-step deuterated-allene preparation recorded with ordered −78 °C, room-temperature, and 140 °C stages and overall 31% yield.", "evidence_ids": [ev(4), ev(4, True)]},
            {"region": "α-diazo ester preparation scheme/procedure", "finding": "Generic transesterification/diazo-transfer sequence plus specific precursor/product preparations begun on this page.", "evidence_ids": [ev(4), ev(4, True)]}],
         "record_ids": page_records[4], "unresolved": ["Cited-literature precursor details and some reagent quantities are not reproduced; [d6]-5 exact isotopic SMILES is not encoded."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "α-diazo ester product blocks", "finding": "Product-specific preparations/yields for 2b–2e, 2g, 2h, 2k, and 2m; analytical peak lists remain source_context.", "evidence_ids": [ev(5), ev(5, True)]},
            {"region": "General oxidative-coupling procedure", "finding": "Full standard catalyst/oxidant/substrate/solvent charges, addition order, 90 °C/10 h hold, and purification captured in the generic record and inherited by scope entries.", "evidence_ids": [ev(5), ev(5, True)]},
            {"region": "Scale-up procedure and 3aa–3ad characterization", "finding": "Scale-up merged with eq 1; isolated masses for Scheme 2 products 3aa–3ad linked to existing scope records.", "evidence_ids": [ev(5), ev(5, True)]}],
         "record_ids": page_records[5], "unresolved": ["Source wording 'allowed to warm to room temperature' after 90 °C was preserved and flagged rather than corrected."]},
        {"page": 6, "state": "attempted", "regions": [
            {"region": "3ae–3ea characterization blocks", "finding": "Isolated masses and yields linked to the same Scheme 2 experiments; NMR/HRMS/IR/Rf lists are characterization context, not new synthesis records.", "evidence_ids": [ev(6), ev(6, True)]}],
         "record_ids": page_records[6], "unresolved": ["Complete machine SMILES for complex dendralenes are intentionally unencoded pending expert review."]},
    ],
}
COVERAGE.write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")

print(json.dumps({"reactions": len(reactions), "molecules": len(molecules), "coverage_pages": len(coverage["pages"])}, ensure_ascii=False))
