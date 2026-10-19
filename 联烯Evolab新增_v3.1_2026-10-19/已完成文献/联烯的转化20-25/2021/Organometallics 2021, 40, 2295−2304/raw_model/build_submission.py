#!/usr/bin/env python3
"""Build the source-derived v3.1 submission for PDF pages 1--6.

This script contains machine extraction candidates.  It does not certify source
accuracy, complete structures, or experimental reproducibility.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DATASET = json.loads((RUN / "dataset.json").read_text())
DOC = "doc_23c0fdb9606e2366"


def ev(page: int, kind: str = "text") -> str:
    return f"{DOC}:page{page}:{kind}"


def F(raw=None, normalized=None, evidence=(), method="original_extraction", derivation=None,
      status=None):
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation, status=status)


def q(number, unit, qualifier=None, raw_text=None):
    value = {"number": number, "unit": unit}
    if qualifier:
        value["qualifier"] = qualifier
    if raw_text:
        value["raw_text"] = raw_text
    return value


bundle = st.empty(DATASET["documents"])
diagnostics = {}
molecules = {}


def molecule(mid, name, label, evidence, smiles=None, chirality="unresolved",
             geometry="unresolved", representation="not_encoded", axial=None,
             identity_note=None, aliases=None):
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, evidence)
    m["names"]["original"] = F(name, name, evidence)
    m["names"]["english"] = F(name, name, evidence)
    m["names"]["abbreviations"] = F(label, [label], evidence) if label else F()
    m["names"]["aliases"] = F(aliases, aliases, evidence) if aliases else F()
    m["cas_numbers"] = F()
    m["external_properties"] = []
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = F(None, smiles, evidence, method="structure_derivation",
                         derivation="Candidate connectivity encoded from the legible source name/drawing; parsing is not source verification.")
        diag = st.stereo_inspect(smiles)
        diagnostics[mid] = diag
        m["stereochemistry"]["atom_centers"] = F(
            None, diag["atom_centers"], evidence, method="structure_derivation",
            derivation="RDKit potential/assigned atom-center diagnostic for the candidate SMILES; source assignment remains separately stated.")
        m["stereochemistry"]["double_bonds"] = F(
            None, diag["double_bonds"], evidence, method="structure_derivation",
            derivation="RDKit double-bond diagnostic for the candidate SMILES; this does not infer missing source E/Z assignments.")
        m["stereochemistry"]["isomeric_smiles"] = F(
            None, diag["canonical_isomeric_smiles"], evidence, method="structure_derivation",
            derivation="Canonical isomeric graph generated from the candidate SMILES without inventing unreported stereochemistry.")
        m["stereochemistry"]["identity_key"] = F(
            None, diag["stereo_graph_key"], evidence, method="structure_derivation",
            derivation="Stereo-preserving graph key generated from the candidate SMILES; not a sample identity.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = F()
        m["stereochemistry"]["atom_centers"] = F()
        m["stereochemistry"]["double_bonds"] = F()
        m["stereochemistry"]["isomeric_smiles"] = F()
        m["stereochemistry"]["identity_key"] = F()
    if axial is not None:
        m["stereochemistry"]["axial"] = F(
            axial.get("raw"), axial.get("normalized"), evidence,
            method=axial.get("method", "structure_derivation"), derivation=axial.get("derivation"))
    else:
        m["stereochemistry"]["axial"] = F()
    m["stereochemistry"]["mixture"] = F()
    if chirality == "unspecified":
        m["stereochemistry"]["assignment"] = F(
            None, "unassigned", evidence, method="inference",
            derivation="No absolute/relative configuration label was found in the inspected source region; this records an unassigned state, not achirality or racemic composition.")
    else:
        m["stereochemistry"]["assignment"] = F()
    for key, value in (("chirality_status", chirality), ("geometry_status", geometry),
                       ("representation_status", representation)):
        m["stereochemistry"][key] = F(
            None, value, evidence, method="inference",
            derivation="Status records what is and is not specified/encoded in the inspected source and candidate representation.")
    note = identity_note or (
        "Exact organometallic bonding/coordination was not encoded; identity is retained by source label and name."
        if not smiles else
        "Candidate structure is machine-derived from the source name/drawing and remains needs_review.")
    m["identity_note"] = F(None, note, evidence, method="inference",
                           derivation="Extraction note describing identity/representation limits; not a verbatim source quotation.")
    molecules[mid] = m
    return mid


# Molecules and materials on the inspected pages. Dense organometallic structures
# remain label-resolved rather than receiving speculative exact SMILES.
molecule("m_1", "1-phenyl-1-propyne", "1", [ev(2), ev(3), ev(6)], "CC#Cc1ccccc1",
         "achiral", "not_applicable", "not_applicable")
molecule("m_2", "ethyl trifluoropyruvate", "2", [ev(2), ev(3), ev(6)],
         "CCOC(=O)C(=O)C(F)(F)F", "achiral", "not_applicable", "not_applicable")
molecule("m_1a", "eta1-allenyl pentamethylcyclopentadienyliron dicarbonyl complex", "1a",
         [ev(2), ev(3), ev(4), ev(5), ev(6)], chirality="unresolved", geometry="unresolved")
molecule("m_1ap", "eta1-propargyliron complex", "1a'", [ev(4), ev(5)],
         chirality="unresolved", geometry="unresolved")
molecule("m_1app", "Cp*Fe(CO)(eta3-1-phenylallen-1-yl)", "1a''", [ev(5)],
         chirality="unresolved", geometry="unresolved")
molecule("m_3", "vinyliron dihydrofuran complex", "3", [ev(2), ev(3), ev(6)],
         chirality="unspecified", geometry="unspecified")
molecule("m_3a", "ethyl 4-bromo-5-phenyl-2-(trifluoromethyl)-2,3-dihydrofuran-2-carboxylate", "3a",
         [ev(3)], "CCOC(=O)C1(C(F)(F)F)CC(Br)=C(c2ccccc2)O1",
         "unspecified", "unspecified", "not_encoded")
molecule("m_3b", "diethyl 5-phenyl-2-(trifluoromethyl)-2,3-dihydrofuran-2,4-dicarboxylate", "3b",
         [ev(3)], "CCOC(=O)C1(C(F)(F)F)CC(C(=O)OCC)=C(c2ccccc2)O1",
         "unspecified", "unspecified", "not_encoded")
molecule("m_3c", "ethyl 4-iodo-5-phenyl-2-(trifluoromethyl)-2,3-dihydrofuran-2-carboxylate", "3c",
         [ev(3)], "CCOC(=O)C1(C(F)(F)F)CC(I)=C(c2ccccc2)O1",
         "unspecified", "unspecified", "not_encoded")
molecule("m_3d", "ethyl 2-hydroxy-5-oxo-5-phenyl-2-(trifluoromethyl)pentanoate", "3d",
         [ev(3)], "CCOC(=O)C(O)(C(F)(F)F)CCC(=O)c1ccccc1",
         "unspecified", "not_applicable", "not_encoded")
molecule("m_4", "1-phenyl-1-butyne", "4", [ev(3)], "CCC#Cc1ccccc1",
         "achiral", "not_applicable", "not_applicable")
conflict_mid = molecule("m_onepot6_substrate_conflict", "conflicting one-pot-6 alkyne identity", "4 / 1-phenyl-1-propyne", [ev(3), ev(6)],
                        chirality="unresolved", geometry="unresolved", representation="not_encoded",
                        identity_note="Page 3 identifies substrate 4 as 1-phenyl-1-butyne, while the corresponding page-6 one-pot procedure says 1-phenyl-1-propyne; no identity is selected by the machine record.")
name_alternatives = [
    F("1-phenyl-1-butyne (4)", "1-phenyl-1-butyne", [ev(3)]),
    F("1-phenyl-1-propyne", "1-phenyl-1-propyne", [ev(6)])
]
for key in ("original", "english"):
    molecules[conflict_mid]["names"][key] = st.fact(
        ["1-phenyl-1-butyne (4)", "1-phenyl-1-propyne"], None, [ev(3), ev(6)],
        status="conflict", alternatives=name_alternatives)
molecules[conflict_mid]["occurrences"][0]["compound_label"] = st.fact(
    ["4", "1-phenyl-1-propyne"], None, [ev(3), ev(6)], status="conflict",
    alternatives=[F("4", "4", [ev(3)]), F("1-phenyl-1-propyne", "1-phenyl-1-propyne", [ev(6)])])
molecule("m_5", "p-toluenesulfonyl isocyanate", "5", [ev(3), ev(6)],
         "Cc1ccc(S(=O)(=O)N=C=O)cc1", "achiral", "not_applicable", "not_applicable",
         aliases=["tosyl isocyanate", "TsNCO"])
molecule("m_6", "1,5-dihydropyrrolone-derived vinyliron complex", "6", [ev(3), ev(6)],
         chirality="unspecified", geometry="unspecified")
molecule("m_6p", "1,3-dihydro-2H-pyrrol-2-one-derived vinyliron complex", "6'", [ev(3)],
         chirality="unspecified", geometry="unspecified")
axial_unassigned = {
    "raw": None,
    "normalized": {"type": "allene_axis", "assignment": "unassigned"},
    "method": "structure_derivation",
    "derivation": "Both allene termini are differently substituted in the source drawing; no axial descriptor is supplied."
}
molecule("m_7a", "eta1-allenyl pentamethylcyclopentadienyliron dicarbonyl complex", "7a", [ev(3), ev(4), ev(6)],
         chirality="unspecified", geometry="unresolved", representation="not_encoded", axial=axial_unassigned)
molecule("m_7b", "eta1-allenyl cyclopentadienyliron dicarbonyl complex", "7b", [ev(4)],
         chirality="unspecified", geometry="unresolved", representation="not_encoded", axial=axial_unassigned)
molecule("m_11", "dialkyl-substituted eta1-allenyliron complex derived from 3-hexyne", "11", [ev(4)],
         chirality="unspecified", geometry="unresolved", representation="not_encoded", axial=axial_unassigned)
molecule("m_6a", "4-iodo-3-methyl-5-phenyl-1-tosyl-1,5-dihydro-2H-pyrrol-2-one", "6a", [ev(3)],
         "Cc1ccc(S(=O)(=O)N2C(=O)C(I)=C(C)C2c2ccccc2)cc1",
         "unspecified", "unspecified", "not_encoded")

for mid, name, label, smi in [
    ("m_9a", "phenyldiazonium tetrafluoroborate", "9a", "N#[N+]c1ccccc1.[B-](F)(F)(F)F"),
    ("m_9b", "4-methoxybenzenediazonium tetrafluoroborate", "9b", "COc1ccc([N+]#N)cc1.[B-](F)(F)(F)F"),
    ("m_9c", "4-nitrobenzenediazonium tetrafluoroborate", "9c", "O=[N+]([O-])c1ccc([N+]#N)cc1.[B-](F)(F)(F)F"),
    ("m_9d", "4-(methoxycarbonyl)benzenediazonium tetrafluoroborate", "9d", "COC(=O)c1ccc([N+]#N)cc1.[B-](F)(F)(F)F"),
    ("m_10a", "5-methyl-1,3-diphenyl-1H-pyrazole", "10a", "Cc1cc(-c2ccccc2)nn1-c1ccccc1"),
    ("m_10b", "1-(4-methoxyphenyl)-5-methyl-3-phenyl-1H-pyrazole", "10b", "COc1ccc(-n2nc(-c3ccccc3)cc2C)cc1"),
    ("m_10c", "5-methyl-1-(4-nitrophenyl)-3-phenyl-1H-pyrazole", "10c", "Cc1cc(-c2ccccc2)nn1-c1ccc([N+](=O)[O-])cc1"),
    ("m_10cp", "minor regioisomer of 10c", "10c'", "Cc1nn(-c2ccc([N+](=O)[O-])cc2)cc1-c1ccccc1"),
    ("m_10d", "methyl 4-(5-methyl-3-phenyl-1H-pyrazol-1-yl)benzoate", "10d", "COC(=O)c1ccc(-n2nc(-c3ccccc3)cc2C)cc1"),
    ("m_12", "3-ethyl-5-methyl-1-(4-nitrophenyl)-1H-pyrazole", "12", "CCc1cc(C)n(-c2ccc([N+](=O)[O-])cc2)n1"),
    ("m_13", "diisopropyl azodicarboxylate", "13", "CC(C)OC(=O)N=NC(=O)OC(C)C"),
    ("m_14", "isopropyl 5-methyl-3-phenyl-1H-pyrazole-1-carboxylate", "14", "CC(C)OC(=O)n1nc(-c2ccccc2)cc1C"),
]:
    molecule(mid, name, label, [ev(4)], smi, "achiral", "not_applicable", "not_applicable")

simple = [
    ("m_fpthf", "[Fp*(thf)]+BF4-", "[Fp*(thf)]+BF4-", None),
    ("m_fpstar_i", "Fp*I", "Fp*I", None),
    ("m_fpstar2", "Fp*2", "Fp*2", None),
    ("m_agbf4", "silver tetrafluoroborate", "AgBF4", None),
    ("m_bf3", "boron trifluoride diethyl etherate", "BF3·Et2O", None),
    ("m_tmph", "2,2,6,6-tetramethylpiperidine", "TMPH", "CC1(C)CCCC(C)(C)N1"),
    ("m_dce", "1,2-dichloroethane", "DCE", "ClCCCl"),
    ("m_ch2cl2", "dichloromethane", "CH2Cl2", "ClCCl"),
    ("m_chcl3", "chloroform", "CHCl3", "ClC(Cl)Cl"),
    ("m_nbs", "N-bromosuccinimide", "NBS", "O=C1CCC(=O)N1Br"),
    ("m_can", "ceric ammonium nitrate", "CAN", None),
    ("m_k2co3", "potassium carbonate", "K2CO3", "O=C([O-])[O-].[K+].[K+]"),
    ("m_etoh", "ethanol", "EtOH", "CCO"),
    ("m_i2", "iodine", "I2", "II"),
    ("m_tfa", "trifluoroacetic acid", "TFA", "O=C(O)C(F)(F)F"),
    ("m_water", "water in wet dichloromethane", "water", "O"),
    ("m_et3n", "triethylamine", "Et3N", "CCN(CC)CC"),
    ("m_et3nhbf4", "triethylammonium tetrafluoroborate", "[Et3NH]+[BF4]-", "CC[NH+](CC)CC.[B-](F)(F)(F)F"),
    ("m_c6d6", "benzene-d6", "C6D6", "[2H]c1c([2H])c([2H])c([2H])c([2H])c1[2H]"),
    ("m_o2", "oxygen", "O2", "[O][O]"),
    ("m_bpo", "benzoyl peroxide", "BzOOBz", "O=C(OOC(=O)c1ccccc1)c1ccccc1"),
    ("m_tempo", "2,2,6,6-tetramethylpiperidine-1-oxyl", "TEMPO", "CC1(C)CCCC(C)(C)N1[O]"),
    ("m_dha", "9,10-dihydroanthracene", "DHA", "c1ccc2c(c1)Cc1ccccc1C2"),
    ("m_co", "carbon monoxide", "CO", "[C-]#[O+]"),
    ("m_silica", "silica gel", "silica gel", None),
    ("m_alumina", "neutral aluminum oxide", "neutral aluminum oxide", None),
    ("m_hexanes", "hexanes", "hexanes", None),
    ("m_etoac", "ethyl acetate", "EtOAc", "CCOC(C)=O"),
]
simple_pages = {
    "m_fpthf":[3,5,6], "m_fpstar_i":[3,4], "m_fpstar2":[5], "m_agbf4":[4],
    "m_bf3":[3,4,6], "m_tmph":[3,6], "m_dce":[3,6], "m_ch2cl2":[3,4,5,6],
    "m_chcl3":[3], "m_nbs":[3], "m_can":[3], "m_k2co3":[3], "m_etoh":[3],
    "m_i2":[3], "m_tfa":[3], "m_water":[3], "m_et3n":[3,4,5,6],
    "m_et3nhbf4":[5], "m_c6d6":[4,5], "m_o2":[5], "m_bpo":[5],
    "m_tempo":[5], "m_dha":[5], "m_co":[5], "m_silica":[6], "m_alumina":[6],
    "m_hexanes":[3,6], "m_etoac":[6]
}
for mid, name, label, smi in simple:
    source_evidence = [ev(page) for page in simple_pages[mid]]
    molecule(mid, name, label, source_evidence, smi,
             "achiral" if smi else "unresolved",
             "not_applicable" if smi else "unresolved",
             "not_applicable" if smi else "not_encoded")


def charge(raw, normalized, evidence, method="original_extraction", derivation=None):
    return F(raw, normalized, evidence, method=method, derivation=derivation)


def participant(pid, mid, role, evidence, stage="reaction", product_role=None, charges=None,
                reference=None):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    if product_role:
        p["product_role"] = F(product_role, product_role, evidence)
    charges = charges or {}
    for key, value in charges.items():
        if len(value) == 2:
            raw, normalized = value
            p["charges"][key] = charge(raw, normalized, evidence)
        else:
            raw, normalized, method, derivation = value
            p["charges"][key] = charge(raw, normalized, evidence, method, derivation)
    if reference:
        p["equivalence_reference"] = F(reference, reference, evidence)
    return p


def step(step_id, evidence, pids, operation, temperature=None, time=None,
         atmosphere=None, pressure=None, stage="reaction"):
    s = st.new_step(step_id, evidence, stage=stage)
    s["participant_ids"] = list(pids)
    s["operation"] = F(None, operation, evidence, method="inference",
                         derivation="Normalized operation reconstructed from the cited source region; exact procedure prose is retained separately when supplied.") if operation else F()
    s["temperature"] = F(temperature[0], temperature[1], evidence) if temperature else F()
    s["time"] = F(time[0], time[1], evidence) if time else F()
    s["atmosphere"] = F(atmosphere, atmosphere, evidence) if atmosphere else F()
    s["pressure"] = F(pressure[0], pressure[1], evidence) if pressure else F()
    return s


def outcome(oid, kind, targets, evidence, raw, normalized, product_role=None,
            yield_type=None, basis=None, scope=None):
    o = st.new_outcome(oid, kind, targets, evidence)
    o["value"] = F(raw, normalized, evidence)
    if product_role:
        o["product_role"] = F(product_role, product_role, evidence)
    if yield_type:
        o["yield_type"] = F(yield_type, yield_type, evidence)
    if basis:
        braw, bnorm, bev = basis
        o["measurement_basis"] = F(braw, bnorm, bev)
    if scope:
        o["scope"] = F(scope, scope, evidence)
    return o


ISO_BASIS = ("Yields refer to isolated compounds, unless otherwise indicated.", "isolated", [ev(6)])


def classification(r, changes=None, coverage="partial", unresolved=False, not_applicable=False,
                   limitations=None):
    c = r["classification"]
    c["mapped_reaction_smiles"] = F()
    c["model_votes"] = []
    if unresolved:
        c["bond_changes"] = []
        c["labels"] = F()
        c["coverage"] = F(None, "unresolved", r["evidence_ids"], method="inference",
                            derivation="Source drawing uses generic or coordination-rich structures that were not atom-mapped.")
    elif not_applicable:
        c["bond_changes"] = []
        c["labels"] = F(None, [], r["evidence_ids"], method="inference",
                          derivation="No observed bond-forming transformation is classified for this record.")
        c["coverage"] = F(None, "not_applicable", r["evidence_ids"], method="inference",
                            derivation="No observed transformation (negative control or contextual/proposed record).")
    else:
        changes = changes or []
        change_facts = []
        labels = {}
        for i, (a, b, before, after, kind) in enumerate(changes, 1):
            cid = f"b{i}"
            obj = {"change_id": cid,
                   "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
                   "before_order": before, "after_order": after, "change_type": kind}
            change_facts.append(F(None, obj, r["evidence_ids"], method="structure_derivation",
                                  derivation="Partial atom-type change inferred from the source structures; exact atom mapping was not asserted."))
            if kind in ("formed", "order_increased", "hydrogen_gain"):
                lab = st.bond_label(a, b, after)
                labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(cid)
        c["bond_changes"] = change_facts
        c["labels"] = F(None, list(labels.values()), r["evidence_ids"], method="structure_derivation",
                          derivation="Taxonomy labels calculated from the explicitly retained partial bond changes.")
        c["coverage"] = F(None, coverage, r["evidence_ids"], method="inference",
                            derivation="Only visually supportable net changes are retained; no complete atom/coproduct mapping is claimed.")
    c["limitations"] = limitations or [
        "No complete atom mapping was asserted for the dense organometallic transformation.",
        "Hydrogen and coproduct sources may be unresolved; classification is partial unless marked not applicable."
    ]


reactions = []


def reaction(suffix, evidence, record_type="specific_experiment", categories=None, location=None,
             name=None, transform=None, participants=None, steps=None, outcomes=None, notes=None,
             general_scheme_ids=None, class_changes=None, class_unresolved=False,
             class_not_applicable=False, class_limitations=None, solvents=None, procedure=None,
             addition_order=None, reactor=None, special=None, unresolved=None):
    rid = f"{DOC}:{suffix}"
    r = st.new_reaction(rid, DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = F(categories, categories, evidence) if categories is not None else F()
    r["location"] = F(location, location, evidence) if location is not None else F()
    r["reaction_name"] = F(name, name, evidence) if name else F()
    if transform:
        r["transformation_type"] = F(None, transform, evidence, method="inference",
                                       derivation="Concise transformation description inferred from the source arrow and prose.")
    r["general_scheme_ids"] = list(general_scheme_ids or [])
    r["participants"] = list(participants or [])
    r["steps"] = list(steps or [])
    r["outcomes"] = list(outcomes or [])
    r["notes"] = [F(None, n, evidence, method="inference",
                    derivation="Extraction/review note explaining record boundaries, ambiguity, or classification limits; not represented as a verbatim source quotation.")
                  for n in (notes or [])]
    if solvents:
        raw, norm, total, concentration = solvents
        r["solvent_conditions"]["composition"] = F(raw, norm, evidence)
        if total:
            r["solvent_conditions"]["total_volume"] = F(total[0], total[1], evidence)
        if concentration:
            r["solvent_conditions"]["concentration"] = F(concentration[0], concentration[1], evidence)
    if procedure:
        r["procedure"]["original_text"] = F(procedure, procedure, evidence)
    if addition_order:
        r["procedure"]["addition_order"] = F(addition_order, addition_order, evidence)
    if reactor:
        r["special_conditions"]["reactor"] = F(reactor, reactor, evidence)
    for key, val in (special or {}).items():
        raw, norm = val
        r["special_conditions"][key] = F(raw, norm, evidence)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": list(unresolved or [])}
    r["completeness"] = "partial"
    classification(r, class_changes, unresolved=class_unresolved,
                   not_applicable=class_not_applicable, limitations=class_limitations)
    reactions.append(r)
    return rid


# Page 1 and Scheme 2 contextual transformations: correctly typed, not counted as
# current preparative experiments.
rids = {}
rids["p1_generic"] = reaction(
    "p1_graphical_abstract_generic", [ev(1)], "generic_scheme", None,
    {"page": 1, "region": "graphical abstract"},
    "graphical abstract: iron-mediated alkyne/dipolar-reagent coupling", "generic overview",
    notes=["Graphical summary only; it does not supply one uniquely parameterized experiment."],
    class_unresolved=True, unresolved=["Generic R/R' structures preclude a concrete molecule registry entry and atom mapping."])
for suffix, region, title in [
    ("p2_scheme2a_literature", "Scheme 2a", "prior propargyliron [3 + 2] cycloadditions"),
    ("p2_scheme2b_literature", "Scheme 2b", "prior allenyliron [3 + 2] cycloadditions"),
    ("p2_scheme2c_prior_work", "Scheme 2c", "previously reported iron-catalyzed alpha-functionalization"),
]:
    rids[suffix] = reaction(
        suffix, [ev(2)], "literature_example", None, {"page": 2, "region": region}, title,
        "generic prior-work transformation", notes=["This panel summarizes cited prior work and is not a new experiment in the present article."],
        class_unresolved=True, unresolved=["R-group/electrophile classes are generic; no concrete experiment or complete bond map is asserted."])


def ids(*parts):
    return list(parts)


# Scheme 3 and page-6 procedures.
e36 = [ev(3), ev(6)]
p_alt3 = [
    participant("p_1a", "m_1a", "reactant", e36, charges={"mass": ("182 mg", q(182, "mg")), "amount": ("0.50 mmol", q(0.50, "mmol")), "equivalents": ("1.0 equiv", q(1.0, "equiv"))}, reference="m_1a"),
    participant("p_2", "m_2", "reactant", e36, charges={"volume": ("99 uL", q(99, "uL")), "amount": ("0.75 mmol", q(0.75, "mmol")), "equivalents": ("1.5 equiv", q(1.5, "equiv"))}, reference="m_1a"),
    participant("p_ch2cl2", "m_ch2cl2", "solvent", e36, charges={"volume": ("0.5 mL", q(0.5, "mL"))}),
    participant("p_3", "m_3", "product", e36, product_role="main product"),
    participant("p_alumina", "m_alumina", "purification_agent", e36, stage="purification"),
    participant("p_hexanes", "m_hexanes", "purification_agent", e36, stage="purification"),
]
rids["s3_alt3"] = reaction(
    "scheme3_1a_2_to_3", e36, categories=["other"], location={"page": 3, "scheme": "3", "arrow": "top"},
    name="alternative synthesis of iron complex 3 from 1a", transform="formal [3 + 2] cycloaddition",
    participants=p_alt3,
    steps=[step("reaction", e36, ids("p_1a", "p_2", "p_ch2cl2", "p_3"), "stirred at room temperature", ("24 °C", q(24, "°C")), ("2.5 h", q(2.5, "h")), "argon/glovebox"),
           step("concentration", [ev(6)], ids("p_3"), "concentrated in vacuo", stage="workup"),
           step("chromatography", [ev(6)], ids("p_3", "p_alumina", "p_hexanes", "p_ch2cl2"), "flash chromatography on neutral aluminum oxide; hexanes/CH2Cl2 1/0 to 10/1 to 0/1, followed by recrystallization", stage="purification")],
    outcomes=[outcome("yield_3", "yield", ["m_3"], e36, "48% yield", q(48, "%"), "main product", "reported yield", ISO_BASIS, {"product": "3"}),
              outcome("amount_3", "product_amount", ["m_3"], [ev(6)], "128 mg", q(128, "mg"), "main product", basis=("isolated solid", "isolated", [ev(6)]), scope={"product": "3"})],
    notes=["Scheme 3 reports >=20:1 rr for complex 3 in the accompanying prose; exact isomer identities are not separately registered."],
    solvents=("CH2Cl2 (0.5 mL)", ["CH2Cl2"], ("0.5 mL", q(0.5, "mL")), None),
    procedure="In an argon-filled glovebox, a threaded reaction tube (5 mL) equipped with a magnetic stir bar was charged with ethyl trifluoropyruvate (99 μL, 0.75 mmol, 1.5 equiv). A solution of freshly prepared allenyliron complex 1a (182 mg, 0.50 mmol, 1.0 equiv) in CH2Cl2 (0.5 mL) was then added, and the reaction tube was fitted with a septum cap and removed from the glovebox. The reaction mixture was then stirred at room temperature (24 °C) for 2.5 h. The crude mixture was concentrated in vacuo and purified by flash column chromatography on neutral aluminum oxide (hexanes/CH2Cl2 1/0 to 10/1 to 0/1) and subsequent recrystallization to provide the desired product 3 as a bright yellow solid. Yield: 128 mg, 48%.",
    addition_order=["ethyl trifluoropyruvate", "solution of 1a in CH2Cl2"], reactor="threaded reaction tube (5 mL) with stir bar",
    class_changes=[("C", "C", 0, 1, "formed"), ("C", "O", 0, 1, "formed")],
    unresolved=["Exact organometallic atom mapping and stereochemical identity of complex 3 remain unresolved."])

p_onepot3 = [
    participant("p_fpthf", "m_fpthf", "reagent", e36, charges={"mass": ("203 mg", q(203, "mg")), "amount": ("0.50 mmol", q(0.50, "mmol")), "equivalents": ("1.0 equiv", q(1.0, "equiv"))}, reference="m_2"),
    participant("p_1", "m_1", "reactant", e36, charges={"volume": ("125 uL", q(125, "uL")), "amount": ("1.0 mmol", q(1.0, "mmol")), "equivalents": ("2.0 equiv", q(2.0, "equiv"))}, reference="m_2"),
    participant("p_2", "m_2", "reactant", e36, charges={"volume": ("66 uL", q(66, "uL")), "amount": ("0.50 mmol", q(0.50, "mmol")), "equivalents": ("1.0 equiv", q(1.0, "equiv"))}, reference="m_2"),
    participant("p_dce", "m_dce", "solvent", e36, charges={"volume": ("1 mL", q(1, "mL"))}),
    participant("p_bf3", "m_bf3", "additive", e36, charges={"volume": ("123 uL", q(123, "uL")), "amount": ("1.0 mmol", q(1.0, "mmol")), "equivalents": ("2.0 equiv", q(2.0, "equiv"))}, reference="m_2"),
    participant("p_tmph", "m_tmph", "reagent", e36, charges={"volume": ("235 uL", q(235, "uL")), "amount": ("1.5 mmol", q(1.5, "mmol")), "equivalents": ("3.0 equiv", q(3.0, "equiv"))}, reference="m_2"),
    participant("p_3", "m_3", "product", e36, product_role="main product"),
    participant("p_silica", "m_silica", "purification_agent", e36, stage="purification"),
    participant("p_hexanes", "m_hexanes", "purification_agent", e36, stage="purification"),
    participant("p_etoac", "m_etoac", "purification_agent", e36, stage="purification"),
]
rids["s3_onepot3"] = reaction(
    "scheme3_onepot_1_2_to_3", e36, categories=["other"], location={"page": 3, "scheme": "3", "arrow": "middle"},
    name="one-pot synthesis of iron complex 3", transform="in situ allenyliron formation and formal [3 + 2] cycloaddition",
    participants=p_onepot3,
    steps=[step("reaction", e36, [p["participant_id"] for p in p_onepot3 if p["stage"] == "reaction"], "BF3·OEt2 and distilled TMPH were added in rapid succession; stirred", ("80 °C", q(80, "°C")), ("12 h", q(12, "h")), "argon/glovebox"),
           step("chromatography", [ev(6)], ids("p_3", "p_silica", "p_hexanes", "p_etoac"), "loaded directly onto silica gel packed with hexanes; chromatography with hexanes/EtOAc 9/1 to 4/1", stage="purification")],
    outcomes=[outcome("yield_3", "yield", ["m_3"], e36, "51%", q(51, "%"), "main product", "reported yield", ISO_BASIS, {"product": "3"}),
              outcome("amount_3", "product_amount", ["m_3"], [ev(6)], "135 mg", q(135, "mg"), "main product", basis=("isolated product", "isolated", [ev(6)]), scope={"product": "3"})],
    solvents=("DCE (1 mL)", ["DCE"], ("1 mL", q(1, "mL")), None),
    procedure="In an argon-filled glovebox, a threaded reaction tube (5 mL) equipped with a magnetic stir bar was charged with [Fp*(thf)]+BF4- (203 mg, 0.50 mmol, 1.0 equiv), 1-phenyl-1-propyne (125 μL, 1.0 mmol, 2.0 equiv), ethyl trifluoropyruvate (66 μL, 0.50 mmol, 1.0 equiv), and DCE (1 mL). Subsequently, BF3·OEt2 (123 μL, 1.0 mmol, 2.0 equiv) and distilled TMPH (235 μL, 1.5 mmol, 3.0 equiv) were added in rapid succession and the reaction tube was fitted with a septum cap and removed from the glovebox. The reaction mixture was stirred for 12 h at 80 °C, and the crude material was loaded directly onto a silica gel column packed with hexanes and purified by chromatography (hexanes/EtOAc 9/1 to 4/1) to give the product 3 as an oily yellow solid that crystallized upon storage under nitrogen. Yield: 135 mg, 51%.",
    addition_order=["[Fp*(thf)]+BF4-, 1, 2, DCE", "BF3·OEt2", "TMPH"], reactor="threaded reaction tube (5 mL) with stir bar",
    class_changes=[("C", "C", 0, 1, "formed"), ("C", "O", 0, 1, "formed")],
    unresolved=["The in situ 1a intermediate is not assigned a quantitative charge or yield.", "Exact organometallic atom mapping remains unresolved."])

def scheme3_deriv(suffix, product_mid, reagent_mid, reagent_pid, reagent_eq, yield_value, changes,
                  solvents, extra_parts=None, extra_outcomes=None, notes=None):
    evidence = [ev(3)]
    parts = [participant("p_3", "m_3", "reactant", evidence, charges={"equivalents": ("1.0 equiv", q(1.0, "equiv"))}, reference="m_3"),
             participant(reagent_pid, reagent_mid, "reagent", evidence, charges={"equivalents": (f"{reagent_eq} equiv", q(reagent_eq, "equiv"))}, reference="m_3"),
             participant("p_product", product_mid, "product", evidence, product_role="main product")]
    parts += extra_parts or []
    solvent_mid = {"CH2Cl2":"m_ch2cl2", "EtOH":"m_etoh", "water":"m_water"}
    represented = {p["molecule_id"] for p in parts}
    for i, solvent_name in enumerate(solvents[1]):
        mid = solvent_mid.get(solvent_name)
        if mid and mid not in represented:
            role = "additive" if solvent_name == "water" else "solvent"
            parts.append(participant(f"p_solvent_{i}", mid, role, evidence))
            represented.add(mid)
    outs = [outcome("yield_main", "yield", [product_mid], evidence, f"{yield_value}%", q(yield_value, "%"), "main product", "reported yield", ISO_BASIS, {"product": product_mid})]
    outs += extra_outcomes or []
    if product_mid == "m_3c":
        reaction_steps = [
            step("low_temperature", evidence, [p["participant_id"] for p in parts], "reaction initiated at 0 °C", ("0 °C", q(0,"°C"))),
            step("warming", evidence, [p["participant_id"] for p in parts], "warmed to 24 °C; transition duration unreported", ("24 °C", q(24,"°C")))
        ]
    else:
        reaction_steps = [step("reaction", evidence, [p["participant_id"] for p in parts], "reaction under the Scheme 3 footnote conditions", ("24 °C", q(24,"°C")))]
    return reaction(suffix, evidence, categories=["derivatization"], location={"page": 3, "scheme": "3", "arrow": suffix},
                    name="divergent transformation of vinyliron complex 3", transform="demetalative functionalization",
                    participants=parts, steps=reaction_steps, outcomes=outs,
                    notes=notes, solvents=solvents, class_changes=changes,
                    unresolved=["Scale, detailed addition order, workup, and purification are unavailable in the supplied page range."])

rids["s3_3a"] = scheme3_deriv("scheme3_3_to_3a", "m_3a", "m_nbs", "p_nbs", 1.0, 54,
    [("C", "Br", 0, 1, "formed")], ("CH2Cl2, 24 °C", ["CH2Cl2"], None, None))
rids["s3_3b"] = scheme3_deriv("scheme3_3_to_3b", "m_3b", "m_can", "p_can", 1.9, 43,
    [("C", "C", 0, 1, "formed")], ("EtOH/CH2Cl2, 24 °C", ["EtOH", "CH2Cl2"], None, None),
    extra_parts=[participant("p_k2co3", "m_k2co3", "reagent", [ev(3)], charges={"equivalents": ("2.9 equiv", q(2.9, "equiv"))}, reference="m_3"), participant("p_etoh", "m_etoh", "solvent", [ev(3)]), participant("p_ch2cl2", "m_ch2cl2", "solvent", [ev(3)])])
rids["s3_3c"] = scheme3_deriv("scheme3_3_to_3c", "m_3c", "m_i2", "p_i2", 1.0, 42,
    [("C", "I", 0, 1, "formed")], ("CH2Cl2, 0-24 °C", ["CH2Cl2"], None, None),
    extra_parts=[participant("p_fpstar_i", "m_fpstar_i", "product", [ev(3)], product_role="recovered iron fragment")],
    extra_outcomes=[outcome("yield_fpstar_i", "yield", ["m_fpstar_i"], [ev(3)], "63% rec. Fp*I", q(63, "%"), "recovered iron fragment", "reported recovery-like yield", ISO_BASIS, {"reported wording": "rec. Fp*I"})],
    notes=["The source labels Fp*I as recovered; because it is formed/released from complex 3 rather than charged as a separate starting molecule, the machine record retains it as a product-specific yield-like outcome, not starting-material recovery."])
rids["s3_3d"] = scheme3_deriv("scheme3_3_to_3d", "m_3d", "m_tfa", "p_tfa", 2.0, 98,
    [("C", "H", 0, 1, "hydrogen_gain"), ("C", "O", 1, 2, "order_increased"), ("O", "H", 0, 1, "hydrogen_gain")],
    ("wet CH2Cl2, 24 °C", ["CH2Cl2", "water"], None, None),
    extra_parts=[participant("p_ch2cl2", "m_ch2cl2", "solvent", [ev(3)]), participant("p_water", "m_water", "additive", [ev(3)])],
    notes=["C-H/O-H and carbonyl formation are partial net classifications; exact proton/water atom origins are not mapped."])


# Scheme 4.
e46 = [ev(3), ev(6)]
p_onepot6 = [
    participant("p_fpthf", "m_fpthf", "reagent", e46, charges={"mass": ("203 mg", q(203, "mg")), "amount": ("0.50 mmol", q(0.50, "mmol")), "equivalents": ("1.0 equiv", q(1.0, "equiv"))}, reference="m_fpthf"),
    participant("p_alkyne_conflict", "m_onepot6_substrate_conflict", "reactant", e46, charges={"volume": ("125 uL", q(125, "uL")), "amount": ("1.0 mmol", q(1.0, "mmol")), "equivalents": ("2.0 equiv", q(2.0, "equiv"))}, reference="m_fpthf"),
    participant("p_5", "m_5", "reactant", e46, charges={"volume": ("82 uL", q(82, "uL")), "amount": ("0.54 mmol", q(0.54, "mmol")), "equivalents": ("1.1 equiv", q(1.1, "equiv"))}, reference="m_fpthf"),
    participant("p_dce", "m_dce", "solvent", e46, charges={"volume": ("1 mL", q(1, "mL"))}),
    participant("p_bf3", "m_bf3", "additive", e46, charges={"volume": ("123 uL", q(123, "uL")), "amount": ("1.0 mmol", q(1.0, "mmol")), "equivalents": ("2.0 equiv", q(2.0, "equiv"))}, reference="m_fpthf"),
    participant("p_tmph", "m_tmph", "reagent", e46, charges={"volume": ("235 uL", q(235, "uL")), "amount": ("1.5 mmol", q(1.5, "mmol")), "equivalents": ("3.0 equiv", q(3.0, "equiv"))}, reference="m_fpthf"),
    participant("p_6", "m_6", "product", e46, product_role="main product"),
    participant("p_silica", "m_silica", "purification_agent", e46, stage="purification"),
    participant("p_etoac", "m_etoac", "purification_agent", e46, stage="purification"),
    participant("p_ch2cl2_pur", "m_ch2cl2", "purification_agent", e46, stage="purification"),
]
rids["s4_onepot6"] = reaction(
    "scheme4_onepot_4_5_to_6", e46, categories=["other"], location={"page": 3, "scheme": "4", "arrow": "top"},
    name="one-pot synthesis of iron complex 6", transform="in situ allenyliron formation and formal [3 + 2] cycloaddition",
    participants=p_onepot6,
    steps=[step("reaction", e46, [p["participant_id"] for p in p_onepot6 if p["stage"] == "reaction"], "BF3·OEt2 and TMPH added in rapid succession; stirred", ("80 °C", q(80, "°C")), ("12 h", q(12, "h")), "argon"),
           step("chromatography", [ev(6)], ids("p_6", "p_silica", "p_etoac", "p_ch2cl2_pur"), "direct chromatography with CH2Cl2/EtOAc 1/0 to 49/1", stage="purification")],
    outcomes=[outcome("yield_6", "yield", ["m_6"], e46, "44%", q(44, "%"), "main product", "reported yield", ISO_BASIS, {"product": "6"}),
              outcome("amount_6", "product_amount", ["m_6"], [ev(6)], "126 mg", q(126, "mg"), "main product", basis=("isolated product", "isolated", [ev(6)]), scope={"product": "6"})],
    solvents=("DCE (1 mL)", ["DCE"], ("1 mL", q(1, "mL")), None),
    procedure="In a septum-capped, argon-filled reaction tube (5 mL) equipped with a magnetic stir bar and charged with [Fp*(thf)]+BF4- (203 mg, 0.50 mmol, 1.0 equiv) were placed 1-phenyl-1-propyne (125 μL, 1.0 mmol, 2.0 equiv), tosyl isocyanate (82 μL, 0.54 mmol, 1.1 equiv), and DCE (1 mL). Then, in rapid succession, BF3·OEt2 (123 μL, 1.0 mmol, 2.0 equiv) and distilled TMPH (235 μL, 1.5 mmol, 3.0 equiv) were added. After it was stirred for 12 h at 80 °C, the crude reaction mixture was purified directly by chromatography (CH2Cl2/EtOAc 1/0 to 49/1) to give the product 6 as an oily yellow solid that crystallized upon storage under nitrogen. Yield: 126 mg, 44%.",
    addition_order=["[Fp*(thf)]+BF4-, 4, 5, DCE", "BF3·OEt2", "TMPH"], reactor="septum-capped reaction tube (5 mL) with stir bar",
    class_changes=[("C", "C", 0, 1, "formed"), ("C", "N", 0, 1, "formed")],
    notes=["Source conflict: page 3 Scheme 4/prose identifies substrate 4 as 1-phenyl-1-butyne, while the corresponding page-6 one-pot procedure states 1-phenyl-1-propyne. Both alternatives are retained without selecting one."],
    unresolved=["Alkyne substrate identity conflicts between pages 3 and 6; the record uses an explicit conflict molecule rather than choosing an identity.","Exact organometallic atom mapping and configuration at the potential stereogenic ring carbon remain unresolved."])

def simple_scheme4(suffix, product, temp, time_raw, time_norm, yld, solvent="CH2Cl2", notes=None):
    evidence=[ev(3)]
    parts=[participant("p_7a", "m_7a", "reactant", evidence, charges={"equivalents": ("1.0 equiv", q(1.0, "equiv"))}, reference="m_7a"),
           participant("p_5", "m_5", "reactant", evidence, charges={"equivalents": ("1.1 equiv", q(1.1, "equiv"))}, reference="m_7a"),
           participant("p_solvent", "m_dce" if solvent=="DCE" else "m_ch2cl2", "solvent", evidence),
           participant("p_product", product, "product", evidence, product_role="main product")]
    return reaction(suffix, evidence, categories=["other"], location={"page":3,"scheme":"4","arrow":suffix},
        name="cycloaddition of isolated allenyliron 7a with tosyl isocyanate", transform="formal [3 + 2] cycloaddition",
        participants=parts, steps=[step("reaction", evidence, [p["participant_id"] for p in parts], "stirred under Scheme 4 conditions", (f"{temp} °C", q(temp,"°C")), (time_raw,time_norm))],
        outcomes=[outcome("yield", "yield", [product], evidence, f"{yld}%", q(yld,"%"), "main product", "reported yield", ISO_BASIS, {"product":product})],
        solvents=(solvent,[solvent],None,None), notes=notes,
        class_changes=[("C","C",0,1,"formed"),("C","N",0,1,"formed")],
        unresolved=["Scale, workup, purification, exact mapping, and absolute configuration are not supplied on the inspected pages."])
rids["s4_6p"] = simple_scheme4("scheme4_7a_5_to_6prime", "m_6p", 24, "1 h", q(1,"h"), 43, "DCE")

e3=[ev(3)]
rids["s4_iso"] = reaction(
    "scheme4_6prime_to_6", e3, categories=["other"], location={"page":3,"scheme":"4","arrow":"6' to 6"},
    name="Et3N-promoted isomerization of 6' to 6", transform="dihydropyrrolone tautomerization",
    participants=[participant("p_6p","m_6p","reactant",e3), participant("p_et3n","m_et3n","reagent",e3), participant("p_6","m_6","product",e3,product_role="main product")],
    steps=[step("reaction",e3,ids("p_6p","p_et3n","p_6"),"treatment with Et3N; temperature and time not reported in Scheme 4")],
    outcomes=[outcome("conversion","conversion",["m_6p"],e3,">95% conversion",q(95,"%","gt",">95%"),scope={"conversion":"6' to 6"})],
    class_changes=[("C","C",1,2,"order_increased"),("C","H",0,1,"hydrogen_gain")],
    unresolved=["Et3N amount, solvent, temperature, time, conversion method, and exact hydrogen/bond migration mapping are unreported in the inspected source."])

# Alternative 7a -> 6 merges the detailed page-6 setup into the Scheme-4 arrow.
p_alt6=[
    participant("p_7a","m_7a","reactant",e46,charges={"mass":("75 mg",q(75,"mg")),"amount":("0.20 mmol",q(0.20,"mmol")),"equivalents":("1.0 equiv",q(1.0,"equiv"))},reference="m_7a"),
    participant("p_5","m_5","reactant",e46,charges={"volume":("37 uL",q(37,"uL")),"amount":("0.22 mmol",q(0.22,"mmol")),"equivalents":("1.1 equiv",q(1.1,"equiv"))},reference="m_7a"),
    participant("p_ch2cl2","m_ch2cl2","solvent",e46,charges={"volume":("0.6 mL",q(0.6,"mL"))}),
    participant("p_6","m_6","product",e46,product_role="main product"),
    participant("p_alumina","m_alumina","purification_agent",e46,stage="purification"),
    participant("p_hexanes","m_hexanes","purification_agent",e46,stage="purification"),
    participant("p_ch2cl2_pur","m_ch2cl2","purification_agent",e46,stage="purification")]
rids["s4_alt6"] = reaction(
    "scheme4_7a_5_to_6",e46,categories=["other"],location={"page":3,"scheme":"4","arrow":"middle right"},
    name="alternative synthesis of iron complex 6 from 7a",transform="formal [3 + 2] cycloaddition with in situ isomerization",
    participants=p_alt6,
    steps=[step("reaction",e46,ids("p_7a","p_5","p_ch2cl2","p_6"),"stirred",("40 °C",q(40,"°C")),("16 h",q(16,"h")),"argon/glovebox"),
           step("concentration",[ev(6)],ids("p_6"),"concentrated in vacuo",stage="workup"),
           step("purification",[ev(6)],ids("p_6","p_alumina","p_hexanes","p_ch2cl2_pur"),"neutral aluminum oxide chromatography (hexanes/CH2Cl2 1/0 to 10/1 to 0/1), then recrystallization from hexanes/CH2Cl2 (10/1 v/v)",stage="purification")],
    outcomes=[outcome("yield","yield",["m_6"],e46,"75%",q(75,"%"),"main product","reported yield",ISO_BASIS,{"product":"6"})],
    solvents=("CH2Cl2 (0.6 mL)",["CH2Cl2"],("0.6 mL",q(0.6,"mL")),None),
    procedure="In an argon-filled glovebox, a threaded reaction tube (5 mL) equipped with a magnetic stir bar was charged with TsNCO (37 μL, 0.22 mmol, 1.1 equiv). A solution of freshly prepared allenyliron complex 7a (75 mg, 0.20 mmol, 1.0 equiv) in CH2Cl2 (0.6 mL) was then added. The reaction tube was fitted with a septum cap and removed from the glovebox. The reaction mixture was then stirred at 40 °C for 16 h. Subsequently, the reaction mixture was concentrated in vacuo, and the crude material was purified by flash column chromatography on neutral aluminum oxide (hexanes/CH2Cl2 1/0 to 10/1 to 0/1) and subsequent recrystallization from hexanes/CH2Cl2 (10/1 v/v) to provide the desired product 6 as a bright yellow solid.",
    addition_order=["TsNCO","solution of 7a in CH2Cl2"],reactor="threaded reaction tube (5 mL) with stir bar",
    class_changes=[("C","C",0,1,"formed"),("C","N",0,1,"formed")],
    unresolved=["The product mass appears on page 7 outside the inspected range and is not imported here; exact mapping/configuration remains unresolved."])

p_6a=[participant("p_6","m_6","reactant",e3,charges={"equivalents":("1.3 equiv",q(1.3,"equiv"))},reference="m_i2"),
      participant("p_i2","m_i2","reagent",e3,charges={"equivalents":("1.0 equiv",q(1.0,"equiv"))},reference="m_i2"),
      participant("p_ch2cl2","m_ch2cl2","solvent",e3), participant("p_chcl3","m_chcl3","solvent",e3),
      participant("p_6a","m_6a","product",e3,product_role="main product"), participant("p_fpstar_i","m_fpstar_i","product",e3,product_role="recovered iron fragment")]
rids["s4_6a"] = reaction(
    "scheme4_6_to_6a",e3,categories=["derivatization"],location={"page":3,"scheme":"4","arrow":"bottom"},
    name="iododemetalation of complex 6",transform="iododemetalation",participants=p_6a,
    steps=[step("low_temperature",e3,[p["participant_id"] for p in p_6a],"reaction initiated at 0 °C",("0 °C",q(0,"°C"))),
           step("warming",e3,[p["participant_id"] for p in p_6a],"warmed to 24 °C; transition duration unreported",("24 °C",q(24,"°C")))],
    outcomes=[outcome("yield_6a","yield",["m_6a"],e3,"59%",q(59,"%"),"main product","reported yield",ISO_BASIS,{"product":"6a"}),
              outcome("yield_fpstar_i","yield",["m_fpstar_i"],e3,"Fp*I, 53%",q(53,"%"),"recovered iron fragment","reported recovery-like yield",ISO_BASIS,{"reported wording":"Fp*I, 53%"})],
    solvents=("CH2Cl2/CHCl3",["CH2Cl2","CHCl3"],None,None),class_changes=[("C","I",0,1,"formed")],
    notes=["The 0 to 24 °C program is retained as two stages; Scheme 4 does not report either stage duration."],
    unresolved=["Stage times, solvent ratio/volume, workup, purification, and absolute configuration at C5 are unreported in the inspected pages."])


# Schemes 5 and 6.
def pyrazole_scope(suffix,diazo,product,yld,rr,temp=60,minor=None,basis_rr="1H NMR spectrum of crude material"):
    evidence=[ev(4)]
    parts=[participant("p_7b","m_7b","reactant",evidence,charges={"equivalents":("1.0 equiv",q(1.0,"equiv"))},reference="m_7b"),
           participant("p_9",diazo,"reactant",evidence,charges={"equivalents":("1.5 equiv",q(1.5,"equiv"))},reference="m_7b"),
           participant("p_ch2cl2","m_ch2cl2","solvent",evidence),
           participant("p_product",product,"product",evidence,product_role="major/main regioisomer")]
    targets=[product]
    if minor:
        parts.append(participant("p_minor",minor,"product",evidence,product_role="minor regioisomer"));targets.append(minor)
    outs=[outcome("yield","yield",[product],evidence,f"{yld}%",q(yld,"%"),"main isolated product","reported yield",ISO_BASIS,{"product":product}),
          outcome("rr","selectivity",targets,evidence,rr,rr,"regioisomeric ratio",basis=(basis_rr,basis_rr,evidence),scope={"metric":"rr","order":[product,minor] if minor else [product,"unobserved regioisomer"],"sample":"crude material"})]
    return reaction(suffix,evidence,categories=["substrate_scope"],location={"page":4,"scheme":"5","entry":suffix},
        name="pyrazole synthesis using cyclopentadienyl-allenyliron cycloaddition to aryldiazonium salt",transform="formal [3 + 2] cycloaddition",
        participants=parts,steps=[step("reaction",evidence,[p["participant_id"] for p in parts],"heated in CH2Cl2",(f"{temp} °C",q(temp,"°C")))],outcomes=outs,
        solvents=("CH2Cl2",["CH2Cl2"],None,None),class_changes=[("C","N",0,1,"formed"),("C","N",0,1,"formed")],
        unresolved=["Reaction time, scale, workup and purification are not supplied on page 4.","Exact organometallic atom mapping is unresolved."])

for key,diazo,prod,yld,rr in [
    ("s5_10a","m_9a","m_10a",27,">100:1 rr"),("s5_10b","m_9b","m_10b",36,">100:1 rr"),
    ("s5_10c","m_9c","m_10c",35,">100:1 rr"),("s5_10d","m_9d","m_10d",35,">100:1 rr")]:
    rids[key]=pyrazole_scope(f"scheme5_{key}",diazo,prod,yld,rr)

e4=[ev(4)]
p12=[participant("p_11","m_11","reactant",e4,charges={"equivalents":("1.0 equiv",q(1.0,"equiv"))},reference="m_11"),
     participant("p_9c","m_9c","reactant",e4,charges={"equivalents":("1.5 equiv",q(1.5,"equiv"))},reference="m_11"),
     participant("p_ch2cl2","m_ch2cl2","solvent",e4),participant("p_12","m_12","product",e4,product_role="major/main regioisomer")]
rids["s5_12"]=reaction("scheme5_11_9c_to_12",e4,categories=["substrate_scope"],location={"page":4,"scheme":"5","entry":"middle"},
    name="pyrazole synthesis from dialkyl allenyliron 11",transform="formal [3 + 2] cycloaddition",participants=p12,
    steps=[step("reaction",e4,[p["participant_id"] for p in p12],"heated in CH2Cl2",("60 °C",q(60,"°C")))],
    outcomes=[outcome("yield","yield",["m_12"],e4,"40%",q(40,"%"),"main isolated product","reported yield",ISO_BASIS,{"product":"12"}),
              outcome("rr","selectivity",["m_12"],e4,"9.0:1 rr","9.0:1","regioisomeric ratio",basis=("Regioisomeric ratios reported are those of the crude material.","crude material; analytical method not stated on page 4",e4),scope={"metric":"rr","order":["12","unregistered minor regioisomer"],"sample":"crude material"})],
    solvents=("CH2Cl2",["CH2Cl2"],None,None),notes=["NOE is shown as connectivity support in Scheme 5; no separate synthesis record is created for the spectroscopy."],
    class_changes=[("C","N",0,1,"formed"),("C","N",0,1,"formed")],
    unresolved=["The minor regioisomer identity and rr analytical method are not specified on page 4.","Time, scale, workup and purification are not supplied in the inspected range."])

for temp,yld,rr in [(4,29,"5.4:1 rr"),(40,38,"4.4:1 rr"),(60,42,"1.5:1 rr")]:
    key=f"s5_7a_{temp}c"
    evidence=e4
    parts=[participant("p_7a","m_7a","reactant",evidence,charges={"equivalents":("1.0 equiv",q(1.0,"equiv"))},reference="m_7a"),
           participant("p_9c","m_9c","reactant",evidence,charges={"equivalents":("1.5 equiv",q(1.5,"equiv"))},reference="m_7a"),
           participant("p_ch2cl2","m_ch2cl2","solvent",evidence),participant("p_10c","m_10c","product",evidence,product_role="major regioisomer"),
           participant("p_10cp","m_10cp","product",evidence,product_role="minor regioisomer")]
    rids[key]=reaction(f"scheme5_7a_9c_{temp}C",evidence,categories=["mechanistic_study"],location={"page":4,"scheme":"5","entry":f"bottom, {temp} °C"},
        name="temperature-dependent regioselectivity with Fp*-allenyl complex 7a",transform="formal [3 + 2] cycloaddition",participants=parts,
        steps=[step("reaction",evidence,[p["participant_id"] for p in parts],"reaction in CH2Cl2",(f"{temp} °C",q(temp,"°C")))],
        outcomes=[outcome("yield","yield",["m_10c","m_10cp"],evidence,f"{yld}%",q(yld,"%"),"combined regioisomer product yield as displayed","reported yield",ISO_BASIS,{"displayed product pair":["10c","10c'"]}),
                  outcome("rr","selectivity",["m_10c","m_10cp"],evidence,rr,rr,"regioisomeric ratio",basis=("Regioisomeric ratios reported are those of the crude material.","crude material; analytical method not stated on page 4",evidence),scope={"metric":"rr","order":["10c","10c'"],"sample":"crude material"})],
        solvents=("CH2Cl2",["CH2Cl2"],None,None),class_changes=[("C","N",0,1,"formed"),("C","N",0,1,"formed")],
        unresolved=["Reaction time, scale, purification and whether displayed yield is combined or assigned solely to 10c are not explicitly stated on page 4."])

p14=[participant("p_7b","m_7b","reactant",e4,charges={"equivalents":("1.0 equiv",q(1.0,"equiv"))},reference="m_7b"),
     participant("p_13","m_13","reactant",e4,charges={"equivalents":("1.5 equiv",q(1.5,"equiv"))},reference="m_7b"),
     participant("p_bf3","m_bf3","additive",e4,charges={"equivalents":("1.5 equiv",q(1.5,"equiv"))},reference="m_7b"),
     participant("p_ch2cl2","m_ch2cl2","solvent",e4),participant("p_14","m_14","product",e4,product_role="main product")]
rids["s6_14"]=reaction("scheme6_7b_13_to_14",e4,categories=["substrate_scope"],location={"page":4,"scheme":"6"},
    name="pyrazole synthesis using diisopropyl azodicarboxylate",transform="formal [3 + 2] cycloaddition with dealkoxycarbonylation",participants=p14,
    steps=[step("reaction",e4,[p["participant_id"] for p in p14],"heated in CH2Cl2 with BF3·Et2O",("60 °C",q(60,"°C")))],
    outcomes=[outcome("yield","yield",["m_14"],e4,"36%",q(36,"%"),"main product","reported yield",ISO_BASIS,{"product":"14"}),
              outcome("rr","selectivity",["m_14"],e4,">20:1 rr",">20:1","regioisomeric ratio",scope={"metric":"rr","order":["14","unregistered minor regioisomer"]})],
    solvents=("CH2Cl2",["CH2Cl2"],None,None),class_changes=[("C","N",0,1,"formed"),("C","N",0,1,"formed")],
    unresolved=["The rr analytical method, reaction time, scale, workup and purification are not stated on page 4."])


# Page-4/5 isomerization observations, Figure 3 and Table 1.
p_crude=[participant("p_1","m_1","reactant",e4),participant("p_fpstar_i","m_fpstar_i","reagent",e4),
         participant("p_agbf4","m_agbf4","reagent",e4),participant("p_et3n","m_et3n","reagent",e4),
         participant("p_1a","m_1a","product",e4,product_role="expected main complex"),participant("p_1ap","m_1ap","product",e4,product_role="variable side product")]
rids["crude_1a"]=reaction("p4_crude_synthesis_1a_with_1aprime",e4,categories=["mechanistic_study"],location={"page":4,"region":"Allenyliron to Propargyliron Isomerization prose"},
    name="crude synthesis of 1a with variable 1a' formation",transform="allenyliron preparation with variable propargyliron side product",participants=p_crude,
    steps=[step("complexation",e4,ids("p_1","p_fpstar_i","p_agbf4"),"reaction of 1-phenyl-1-propyne with Fp*I/AgBF4"),step("deprotonation",e4,ids("p_et3n","p_1a","p_1ap"),"followed by deprotonation with Et3N")],
    outcomes=[outcome("side_product_fraction","selectivity",["m_1a","m_1ap"],e4,"trace amount (<1%) to ca. 70% of a second iron complex depending on batch",{"minimum":1,"maximum":70,"unit":"%","minimum_qualifier":"lt","maximum_qualifier":"approx"},"batch-dependent product composition",scope={"metric":"1a' fraction","range":"<1% to ca.70%","order":["1a'","1a"]})],
    class_changes=[("C","Fe",0,1,"formed")],notes=["The exact batch conditions and individual values are not reported in the main text."],
    unresolved=["Amounts, solvent, temperature, time, workup, purification and composition measurement basis require the unavailable SI.","The partial C-Fe label does not claim a complete preparation mapping."])

rids["standing_eq"]=reaction("p4_standing_1a_to_1aprime",e4,categories=["mechanistic_study"],location={"page":4,"region":"Allenyliron to Propargyliron Isomerization prose"},
    name="standing isomerization/equilibration of crude 1a",transform="allenyliron-to-propargyliron isomerization",participants=[participant("p_1a","m_1a","reactant",e4),participant("p_1ap","m_1ap","product",e4,product_role="predominant equilibrium species"),participant("p_ch2cl2","m_ch2cl2","solvent",e4),participant("p_c6d6","m_c6d6","solvent",e4)],
    steps=[step("standing",e4,ids("p_1a","p_1ap","p_ch2cl2","p_c6d6"),"standing as a solid or as a solution in CD2Cl2 or C6D6; time and temperature not reported")],
    outcomes=[outcome("ratio","selectivity",["m_1ap","m_1a"],e4,"propargyliron:allenyliron approximately 92:8","92:8","equilibrium-like composition",scope={"metric":"propargyliron:allenyliron ratio","order":["1a'","1a"]})],
    solvents=("solid or solution in CD2Cl2 or C6D6",["solid","CD2Cl2","C6D6"],None,None),
    class_changes=[("C","Fe",0,1,"formed"),("C","C",2,3,"order_increased")],
    unresolved=["Individual solvent results, time, temperature, sample history and ratio measurement method are not separately reported."])

e45=[ev(4),ev(5)]
kin_common=[participant("p_1a","m_1a","reactant",e45),participant("p_1ap","m_1ap","product",e45,product_role="isomerization product"),participant("p_ch2cl2","m_ch2cl2","solvent",e45),participant("p_c6d6","m_c6d6","solvent",e45)]
rids["kin_bg"]=reaction("figure3_background_isomerization",e45,categories=["mechanistic_study"],location={"page":5,"figure":"3","series":"background"},
    name="background isomerization of purified 1a",transform="allenyliron-to-propargyliron isomerization",participants=kin_common,
    steps=[step("reaction",e45,[p["participant_id"] for p in kin_common],"monitored with rigorous exclusion of light and oxygen",("24 °C",q(24,"°C")),("24 h",q(24,"h")),"light and oxygen rigorously excluded")],
    outcomes=[outcome("conversion","conversion",["m_1a"],e45,"<10% conversion after 24 h",q(10,"%","lt","<10%"),scope={"product":"1a'"})],
    solvents=("CD2Cl2 or C6D6",["CD2Cl2","C6D6"],None,None),class_changes=[("C","Fe",0,1,"formed"),("C","C",2,3,"order_increased")],
    unresolved=["The exact solvent-specific curve and conversion measurement basis are not separately stated; Figure 3 is not digitized into point values."])

def kinetic_additive(suffix, additive, loading_raw, loading_norm, raw_result, norm_result, evidence, extra_unresolved=None):
    parts=[participant("p_1a","m_1a","reactant",evidence),participant("p_add",additive,"additive",evidence,charges={"mol_percent":(loading_raw,loading_norm)},reference="m_1a"),
           participant("p_ch2cl2","m_ch2cl2","solvent",evidence),participant("p_1ap","m_1ap","product",evidence,product_role="isomerization product")]
    return reaction(suffix,evidence,categories=["mechanistic_study"],location={"page":5,"figure":"3","series":suffix},
        name="additive effect on 1a isomerization",transform="allenyliron-to-propargyliron isomerization",participants=parts,
        steps=[step("reaction",evidence,[p["participant_id"] for p in parts],"reaction progress monitored by 1H NMR spectroscopy",("24 °C",q(24,"°C")))],
        outcomes=[outcome("rate_effect","other",["m_1a","m_1ap"],evidence,raw_result,norm_result,scope={"metric":"rate effect versus background"})],
        solvents=("CD2Cl2",["CD2Cl2"],None,None),class_changes=[("C","Fe",0,1,"formed"),("C","C",2,3,"order_increased")],
        unresolved=["Exact concentration, time-course point values and kinetic fit are not tabulated on the inspected pages."]+list(extra_unresolved or []))
rids["kin_fp"]=kinetic_additive("figure3_fp_cation_15mol", "m_fpthf", "15 mol%", q(15,"mol%"), "ca. 30 times faster initial rate than uncatalyzed", "ca. 30-fold", [ev(5)])
rids["kin_acid"]=kinetic_additive("figure3_et3nhbf4_100mol", "m_et3nhbf4", "100 mol%", q(100,"mol%"), "significant rate acceleration, less dramatic than [Fp*(thf)]+[BF4]-", "significant; less than Fp* cation", [ev(5)])
rids["kin_et3n"]=kinetic_additive("figure3_et3n_control", "m_et3n", "Et3N (loading not reported)", None, "Et3N did not promote the allenyl-to-propargyl isomerization", "no promotion relative to background", [ev(5),ev(6)], ["Et3N loading is not stated in the main text/figure."])

table_rows=[
    (1,None,None,51.0,0.0125,85),(2,"m_o2","O2 purged",51.0,0.0125,64),(3,"m_bpo","BzOOBz 20%",51.0,0.0125,65),
    (4,None,None,45.0,0.010,60),(5,"m_tempo","TEMPO 5%",45.0,0.010,7),(6,"m_dha","DHA 50%",45.0,0.010,44)]
for entry,scav,scav_raw,temp,amount,conv in table_rows:
    evidence=[ev(5)]
    parts=[participant("p_1a","m_1a","reactant",evidence,charges={"amount":(f"{amount} mmol",q(amount,"mmol"))}),participant("p_c6d6","m_c6d6","solvent",evidence,charges={"volume":("0.5 mL",q(0.5,"mL"))}),participant("p_1ap","m_1ap","product",evidence,product_role="isomerization product")]
    unresolved=["Conversion measurement basis is not stated in Table 1."]
    if scav:
        if entry==2:
            parts.append(participant("p_scav",scav,"additive",evidence))
        else:
            percent=float(scav_raw.split()[-1].rstrip("%"))
            parts.append(participant("p_scav",scav,"additive",evidence,charges={"mol_percent":(scav_raw,q(percent,"mol%"),"inference","The table gives a bare percentage in the radical-scavenger column; normalized as mol% relative to 1a for machine structure, but the source does not explicitly state the percentage basis.")},reference="m_1a"))
            unresolved.append("The scavenger percentage basis/type is not explicitly stated; mol% normalization is inferred and needs review.")
    rids[f"table1_{entry}"]=reaction(f"table1_entry_{entry}",evidence,categories=["mechanistic_study"],location={"page":5,"table":"1","entry":entry},
        name="radical scavenger effect on 1a isomerization",transform="allenyliron-to-propargyliron isomerization",participants=parts,
        steps=[step("reaction",evidence,[p["participant_id"] for p in parts],"heated in C6D6 with the listed radical scavenger",(f"{temp} °C",q(temp,"°C")),("1 h",q(1,"h")),"O2 purged" if entry==2 else None)],
        outcomes=[outcome("conversion","conversion",["m_1a"],evidence,f"{conv}% conversion",q(conv,"%"),scope={"product":"1a'"})],
        solvents=("C6D6 (0.5 mL)",["C6D6"],("0.5 mL",q(0.5,"mL")),None),class_changes=[("C","Fe",0,1,"formed"),("C","C",2,3,"order_increased")],
        unresolved=unresolved)

rids["scheme7"]=reaction("scheme7_proposed_isomerization_mechanisms",[ev(5)],"proposed_mechanism",None,{"page":5,"scheme":"7"},
    "plausible mechanisms for allenyl-to-propargyl isomerization","proposed radical, eta3-intermediate, and electrophilic pathways",
    notes=["Scheme 7 contains proposed pathways I-III; it is not itself an observed experiment."],class_unresolved=True,
    unresolved=["Elementary-step bond changes are proposals and are not promoted to observed reaction classifications."])

p_photo=[participant("p_1a","m_1a","reactant",[ev(5)],charges={"amount":("0.03 mmol",q(0.03,"mmol"))}),participant("p_c6d6","m_c6d6","solvent",[ev(5)],charges={"volume":("0.5 mL",q(0.5,"mL"))}),participant("p_1app","m_1app","product",[ev(5)],product_role="photoproduct")]
rids["photo"]=reaction("scheme8_photolysis_1a_to_1app",[ev(5)],categories=["mechanistic_study"],location={"page":5,"scheme":"8","arrow":"top"},
    name="irradiation of 1a",transform="photodecarbonylation to eta3 complex",participants=p_photo,
    steps=[step("irradiation",[ev(5)],[p["participant_id"] for p in p_photo],"irradiated with visible light",time=("15 min",q(15,"min")))],
    outcomes=[outcome("yield","yield",["m_1app"],[ev(5)],"quantitative","quantitative","photoproduct",scope={"product":"1a''"})],
    solvents=("C6D6 (0.5 mL)",["C6D6"],("0.5 mL",q(0.5,"mL")),None),special={"light_source":("85 W CFL","85 W CFL"),"wavelength":("visible","visible"),"power":("85 W",q(85,"W"))},
    class_unresolved=True,unresolved=["Yield measurement basis and exact eta3 metal-carbon bond representation are not specified; photoproduct spectroscopy is retained as source context."])

p_photo2=[participant("p_1a","m_1a","reactant",[ev(5)],charges={"amount":("0.03 mmol",q(0.03,"mmol"))}),participant("p_fpstar2","m_fpstar2","additive",[ev(5)],charges={"mol_percent":("1 mol%",q(1,"mol%"))},reference="m_1a"),participant("p_c6d6","m_c6d6","solvent",[ev(5)],charges={"volume":("0.5 mL",q(0.5,"mL"))}),participant("p_1app","m_1app","product",[ev(5)],product_role="major product"),participant("p_1ap","m_1ap","product",[ev(5)],product_role="minor product")]
rids["photo_fp2"]=reaction("scheme8_photolysis_with_fpstar2",[ev(5)],categories=["mechanistic_study"],location={"page":5,"scheme":"8","arrow":"middle"},
    name="irradiation of 1a in the presence of Fp*2",transform="photochemical formation of eta3 and propargyliron products",participants=p_photo2,
    steps=[step("irradiation",[ev(5)],[p["participant_id"] for p in p_photo2],"irradiated with visible light",time=("15 min",q(15,"min")))],
    outcomes=[outcome("ratio","selectivity",["m_1app","m_1ap"],[ev(5)],"4.8:1","4.8:1","product ratio",scope={"metric":"product ratio","order":["1a''","1a'"]})],
    solvents=("C6D6 (0.5 mL)",["C6D6"],("0.5 mL",q(0.5,"mL")),None),special={"light_source":("85 W CFL","85 W CFL"),"wavelength":("visible","visible"),"power":("85 W",q(85,"W"))},
    class_unresolved=True,unresolved=["The ratio measurement method and combined yield are not reported in Scheme 8; exact metal coordination prevents reliable bond mapping."])

p_co=[participant("p_1app","m_1app","reactant",[ev(5)]),participant("p_co","m_co","reagent",[ev(5)]),participant("p_c6d6","m_c6d6","solvent",[ev(5)])]
rids["co_control"]=reaction("scheme8_co_recarbonylation_control",[ev(5)],categories=["control","mechanistic_study"],location={"page":5,"scheme":"8","arrow":"bottom crossed arrow"},
    name="CO recarbonylation control of 1a''",transform="no observed conversion to 1a'",participants=p_co,
    steps=[step("co_exposure",[ev(5)],[p["participant_id"] for p in p_co],"CO bubbled in an NMR tube or stirred under CO; dark or CFL irradiation",atmosphere="CO",pressure=("1 atm","1 atm"))],
    outcomes=[outcome("no_conversion","no_reaction",["m_1app"],[ev(5)],"no conversion to propargyliron 1a' under all conditions examined","no observed conversion",scope={"tested product":"1a'","variants":["dark","CFL irradiation","CO bubbling","CO atmosphere"]})],
    solvents=("C6D6",["C6D6"],None,None),special={"light_source":("dark or CFL","dark or CFL")},class_not_applicable=True,
    unresolved=["Individual time, temperature, scale, and analytical results for each dark/light and vessel variant are not separately reported."])


bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id":f"{DOC}:p1_scheme1","document_id":DOC,"evidence_ids":[ev(1)],"kind":"nonreaction_structures","description":"Scheme 1 contains therapeutic/natural-product structures and is not a current reaction experiment."},
    {"context_id":f"{DOC}:p3_xrd","document_id":DOC,"evidence_ids":[ev(3)],"kind":"structure_characterization","description":"Figures 1 and 2 are ORTEP structures of complexes 3 and 6; no new synthesis record is created for the diffraction plots."},
    {"context_id":f"{DOC}:p5_1app_spectra","document_id":DOC,"evidence_ids":[ev(5)],"kind":"spectroscopy","description":"Scheme 8 lists selected NMR/IR properties of 1a''; spectral peaks remain characterization context, not a standalone synthesis."},
    {"context_id":f"{DOC}:p6_general","document_id":DOC,"evidence_ids":[ev(6)],"kind":"general_information","description":"General experimental/characterization methods and the statement that yields are isolated unless otherwise indicated."},
    {"context_id":f"{DOC}:p6_unresolved_controls","document_id":DOC,"evidence_ids":[ev(6)],"kind":"insufficiently_parameterized_experiments","description":"BF4- sources, Ag powder, and AgI are said not to promote isomerization, and no significant isomerization is said for complexes derived from 1-phenyl-1-butyne or 3-hexyne under a variety of conditions; individual conditions/results are in unavailable SI and are not split into invented records."}
]

# Page coverage: exactly pages 1--6. Record IDs are based on evidence membership,
# with page-6 detailed procedures linked to their existing Scheme records.
page_record_ids = {p: [] for p in range(1,7)}
for r in reactions:
    for p in range(1,7):
        if ev(p) in r["evidence_ids"]:
            page_record_ids[p].append(r["reaction_id"])

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page":1,"state":"attempted","regions":[
            {"region":"title/abstract and graphical abstract","finding":"Title/metadata verified; graphical abstract retained as one generic scheme, not a specific experiment.","evidence_ids":[ev(1),ev(1,"image")]},
            {"region":"Scheme 1 and introduction","finding":"Drug/natural-product structures and introductory discussion; no current preparative experiment.","evidence_ids":[ev(1),ev(1,"image")]}
        ],"record_ids":page_record_ids[1],"unresolved":["Graphical R/R' structures are generic and not atom-mapped."]},
        {"page":2,"state":"attempted","regions":[
            {"region":"Scheme 2a","finding":"Cited prior propargyliron cycloaddition classes; literature example only.","evidence_ids":[ev(2),ev(2,"image")]},
            {"region":"Scheme 2b","finding":"Cited prior allenyliron cycloaddition classes; literature example only.","evidence_ids":[ev(2),ev(2,"image")]},
            {"region":"Scheme 2c and surrounding prose","finding":"Authors' previously reported alpha-functionalization concept; prior work, not a present experiment.","evidence_ids":[ev(2),ev(2,"image")]},
            {"region":"Results opening","finding":"Narrative introduces current 1a + 2 -> 3 experiment; concrete arrow appears on page 3.","evidence_ids":[ev(2),ev(2,"image")]}
        ],"record_ids":page_record_ids[2],"unresolved":["Scheme 2 uses generic R groups/electrophile classes and cannot support concrete structures."]},
        {"page":3,"state":"attempted","regions":[
            {"region":"Scheme 3 and footnotes a-f","finding":"Six independent current reaction arrows recorded: two syntheses of 3 and four divergent transformations 3a-d, including Fp*I outcome for 3c.","evidence_ids":[ev(3),ev(3,"image")]},
            {"region":"Scheme 4 and footnotes a-d","finding":"Five independent current arrows recorded: one-pot 6, 7a->6', 6'->6, 7a->6, and iododemetalation to 6a/Fp*I.","evidence_ids":[ev(3),ev(3,"image")]},
            {"region":"Figures 1 and 2","finding":"ORTEP characterization only; retained in source_context, not additional syntheses.","evidence_ids":[ev(3),ev(3,"image")]}
        ],"record_ids":page_record_ids[3],"unresolved":["Dense organometallic structures and exact stereochemical assignments are not encoded; detailed SI procedures are unavailable."]},
        {"page":4,"state":"attempted","regions":[
            {"region":"Scheme 5 top/middle/bottom and footnote","finding":"Eight distinct experiments recorded: four 7b/diazonium scope entries, one 11/9c entry, and three temperature variants for 7a/9c.","evidence_ids":[ev(4),ev(4,"image")]},
            {"region":"Scheme 6 and footnote","finding":"One 7b + DIAD/BF3·Et2O experiment recorded with yield and rr.","evidence_ids":[ev(4),ev(4,"image")]},
            {"region":"Allenyliron-to-propargyliron prose","finding":"Variable crude 1a/1a' formation and standing 92:8 composition recorded separately; purified background experiment continues on page 5.","evidence_ids":[ev(4),ev(4,"image")]}
        ],"record_ids":page_record_ids[4],"unresolved":["Reaction times/workups for Schemes 5/6 and several rr methods are not supplied on page 4; SI is not associated."]},
        {"page":5,"state":"attempted","regions":[
            {"region":"Figure 3 and kinetics prose","finding":"Background, 15 mol% Fp* cation, 100 mol% Et3NHBF4, and Et3N comparison series recorded without digitizing graph points.","evidence_ids":[ev(5),ev(5,"image")]},
            {"region":"Scheme 7","finding":"Three plausible mechanistic pathways retained as one proposed_mechanism record, not experimental proof.","evidence_ids":[ev(5),ev(5,"image")]},
            {"region":"Table 1 rows 1-6 and footnotes","finding":"Six separate radical-scavenger experiments recorded with row-specific scale, temperature, time, additive, and conversion.","evidence_ids":[ev(5),ev(5,"image")]},
            {"region":"Scheme 8 top/middle/bottom","finding":"Two photochemical experiments and one CO no-conversion control recorded; displayed spectral values retained as characterization context.","evidence_ids":[ev(5),ev(5,"image")]}
        ],"record_ids":page_record_ids[5],"unresolved":["Figure 3 point values were not digitized; Table 1 conversion method and scavenger percentage basis are not explicit."]},
        {"page":6,"state":"attempted","regions":[
            {"region":"isomerization controls and conclusion","finding":"Et3N control linked to Figure 3 record; other SI-only negative controls and vague 'variety of conditions' observations retained in source_context rather than invented individual experiments.","evidence_ids":[ev(6),ev(6,"image")]},
            {"region":"General Information","finding":"Global isolated-yield basis and analytical/purification context captured; standalone spectral methods are not reactions.","evidence_ids":[ev(6),ev(6,"image")]},
            {"region":"One-pot synthesis of 3","finding":"Detailed quantities, order, 80 °C/12 h, chromatography, mass and 51% yield merged into existing Scheme 3 record.","evidence_ids":[ev(6),ev(6,"image")]},
            {"region":"Alternative synthesis of 3 from 1a","finding":"Detailed quantities, 24 °C/2.5 h, concentration, chromatography/recrystallization, mass and 48% yield merged into existing Scheme 3 record.","evidence_ids":[ev(6),ev(6,"image")]},
            {"region":"One-pot synthesis of 6 and start of alternative 7a route","finding":"Detailed one-pot 6 procedure and the 7a/5 setup through 40 °C/16 h and purification merged into their existing Scheme 4 records.","evidence_ids":[ev(6),ev(6,"image")]}
        ],"record_ids":page_record_ids[6],"unresolved":["The alternative 7a->6 product mass is continued on page 7 and was not imported; SI-only controls lack parameterized conditions."]}
    ]
}

visual = {
    "document_id": DOC,
    "inspected_original_pdf_pages": [1,2,3,4,5,6],
    "chemical_accuracy_verified": False,
    "scheme_arrow_counts": {"Scheme 3":6,"Scheme 4":5,"Scheme 5":8,"Scheme 6":1,"Scheme 8 observed/control arrows":3},
    "table_rows": {"Table 1":6},
    "nonexperimental_regions": ["Scheme 1 structures","Figures 1/2 ORTEP","Scheme 7 proposed mechanism","Scheme 8 spectral annotations"],
    "notes": ["Counts describe visually inspected arrows/rows, not a claim of whole-paper semantic completeness.","Current experiments repeated in page-6 prose were merged with Scheme records rather than duplicated."]
}

(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/stereo_diagnostics.json").write_text(json.dumps(diagnostics, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/visual_transcription.json").write_text(json.dumps(visual, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/source_mentions.json").write_text(json.dumps(
    st.source_mentions(DATASET, ["1a", "1a'", "3a", "3b", "3c", "3d", "6a", "10c", "12", "14"], context=180),
    ensure_ascii=False, indent=2) + "\n")
(RUN / "coverage/pages_0001_0006.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"molecules":len(bundle["molecules"]),"reactions":len(bundle["reactions"]),"specific":sum(r["record_type"]=="specific_experiment" for r in bundle["reactions"]),"stereo_diagnostics":len(diagnostics)}, ensure_ascii=False))
