#!/usr/bin/env python3
"""Build the source-derived v3.1 submission for doc_b1bd28de50b7c405.

This is a machine-authored extraction artifact. It does not promote any record
to source-reviewed status and does not alter the extraction tools or schemas.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))

import standard_records as st


DOC = "doc_b1bd28de50b7c405"
E = {page: f"{DOC}:page{page}:text" for page in range(1, 6)}
EI = {page: f"{DOC}:page{page}:image" for page in range(1, 6)}


def f(raw=None, normalized=None, evidence=(), method="original_extraction", derivation=None, status=None):
    return st.fact(
        raw,
        normalized,
        evidence,
        method=method,
        derivation=derivation,
        status=status,
    )


def inferred(normalized, evidence, derivation):
    return f(None, normalized, evidence, method="inference", derivation=derivation)


def structural(normalized, evidence, derivation):
    return f(None, normalized, evidence, method="structure_derivation", derivation=derivation)


def missing():
    return f()


def molecule(
    mid,
    label,
    name,
    evidence,
    *,
    smiles=None,
    chirality="unresolved",
    geometry="unresolved",
    representation="unresolved",
    atom_centers=None,
    double_bonds=None,
    axial=None,
    mixture=None,
    assignment=None,
    identity_note=None,
    allene=False,
):
    ev = list(evidence)
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    m["names"]["original"] = f(name, name, ev)
    m["names"]["english"] = f(name, name, ev)
    m["names"]["abbreviations"] = f(label, [label], ev)
    m["names"]["aliases"] = missing()
    m["cas_numbers"] = missing()
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = structural(
            smiles,
            ev,
            "Candidate graph transcribed from the source drawing/name; RDKit parsing is a graph diagnostic, not source identity verification.",
        )
        diag = st.stereo_inspect(smiles)
        m["stereochemistry"]["isomeric_smiles"] = structural(
            diag["canonical_isomeric_smiles"], ev, "Canonical isomeric SMILES generated from the candidate source-derived graph."
        )
        m["stereochemistry"]["identity_key"] = structural(
            diag["stereo_graph_key"], ev, "Stereo graph key calculated from the candidate isomeric SMILES."
        )
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing()
        m["stereochemistry"]["isomeric_smiles"] = missing()
        m["stereochemistry"]["identity_key"] = missing()
    m["generic_template"] = None
    m["stereochemistry"]["atom_centers"] = (
        f(atom_centers, atom_centers, ev) if atom_centers is not None else missing()
    )
    m["stereochemistry"]["double_bonds"] = (
        f(double_bonds, double_bonds, ev) if double_bonds is not None else missing()
    )
    m["stereochemistry"]["axial"] = f(axial, axial, ev) if axial is not None else missing()
    m["stereochemistry"]["mixture"] = f(mixture, mixture, ev) if mixture is not None else missing()
    m["stereochemistry"]["assignment"] = f(assignment, assignment, ev) if assignment is not None else missing()
    m["stereochemistry"]["chirality_status"] = inferred(
        chirality, ev, "Status assigned from the depicted/name-resolved candidate while preserving any unassigned stereogenic elements."
    )
    m["stereochemistry"]["geometry_status"] = inferred(
        geometry, ev, "Status assigned from the depicted/name-resolved candidate; absence of a source descriptor is not converted to a single isomer."
    )
    m["stereochemistry"]["representation_status"] = inferred(
        representation, ev, "Compared source stereochemical claims with the candidate machine representation."
    )
    m["identity_note"] = (
        f(identity_note, identity_note, ev) if identity_note else missing()
    )
    m["external_properties"] = []
    if allene:
        m["allene_features"] = [
            {
                "fragment_id": f"{mid}:allene1",
                "carbon_atom_maps": structural(
                    {"candidate_smiles_atom_indices": [0, 1, 2]},
                    ev,
                    "Indices refer to the candidate SMILES, not an atom mapping asserted by the paper.",
                ),
                "substitution": structural(
                    "terminal CH2=C=C(H)-N motif",
                    ev,
                    "Derived from the terminal allenamide drawing.",
                ),
                "axial_stereochemistry": inferred(
                    "not stereogenic at the terminal CH2 end",
                    ev,
                    "A terminal CH2 end has two identical hydrogen substituents; no axial descriptor was reported.",
                ),
                "local_environment": structural(
                    "N-substituted allenamide",
                    ev,
                    "Derived from the source drawing and allenamide label.",
                ),
            }
        ]
    return m


def generic_material(mid, label, name, evidence, fixed, reason, *, chirality="unresolved", axial=None, assignment=None):
    ev = list(evidence)
    m = molecule(
        mid,
        label,
        name,
        ev,
        chirality=chirality,
        geometry="not_applicable",
        representation="not_encoded",
        axial=axial,
        assignment=assignment,
        identity_note=reason,
    )
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": f(fixed, fixed, ev),
        "variables": [],
        "instance_mapping": [],
        "reason": f(reason, reason, ev),
    }
    return m


def participant(pid, mid, role, evidence, *, mol_percent=None, equivalents=None, reference=None, stage="reaction"):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    if mol_percent is not None:
        p["charges"]["mol_percent"] = f(
            f"{mol_percent:g} mol%",
            {"number": mol_percent, "unit": "mol%"},
            evidence,
        )
        p["equivalence_reference"] = f(reference, reference, evidence)
    if equivalents is not None:
        p["charges"]["equivalents"] = f(
            f"{equivalents:g} equiv",
            {"number": equivalents, "unit": "equiv"},
            evidence,
        )
        p["equivalence_reference"] = f(reference, reference, evidence)
    return p


def base_reaction(rid, record_type, categories, evidence, location, name, transformation, unresolved):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, evidence)
    r["location"] = f(location, location, evidence)
    r["reaction_name"] = f(name, name, evidence)
    r["transformation_type"] = inferred(
        transformation,
        evidence,
        "Normalized transformation description derived from the source reaction name and drawings.",
    )
    r["review"] = {
        "status": "needs_review",
        "method": "model",
        "reviewer": None,
        "checked_fields": [],
        "unresolved": list(unresolved),
    }
    r["completeness"] = "partial"
    return r


def step(step_id, evidence, participants, *, operation=None, temperature=None, time=None, atmosphere=None, stage="reaction"):
    s = st.new_step(step_id, evidence, stage=stage)
    s["participant_ids"] = list(participants)
    if operation is not None:
        s["operation"] = f(operation, operation, evidence)
    if temperature is not None:
        normalized = temperature
        if isinstance(temperature, tuple):
            raw, normalized = temperature
        else:
            raw = temperature
        s["temperature"] = f(raw, normalized, evidence)
    if time is not None:
        normalized = time
        if isinstance(time, tuple):
            raw, normalized = time
        else:
            raw = time
        s["time"] = f(raw, normalized, evidence)
    if atmosphere is not None:
        s["atmosphere"] = f(atmosphere, atmosphere, evidence)
    return s


def outcome(outcome_id, kind, targets, evidence, value, *, basis=None, yield_type=None, scope=None, role="main product"):
    o = st.new_outcome(outcome_id, kind, targets, evidence)
    o["product_role"] = f(role, role, evidence)
    raw_value, normalized_value = value if isinstance(value, tuple) else (value, value)
    o["value"] = f(raw_value, normalized_value, evidence)
    if basis is not None:
        o["measurement_basis"] = f(basis, basis, evidence)
    if yield_type is not None:
        o["yield_type"] = f(yield_type, yield_type, evidence)
    if scope is not None:
        o["scope"] = f(scope, scope, evidence)
    return o


def set_partial_classification(r, evidence, changes, limitation):
    facts = []
    expected = {}
    for item in changes:
        facts.append(
            structural(
                item,
                evidence,
                "Net bond change inferred from the source reactant/product drawings; no complete atom-mapped reaction is asserted.",
            )
        )
        if item["change_type"] in ("formed", "order_increased", "hydrogen_gain"):
            label = st.bond_label(item["atoms"][0]["element"], item["atoms"][1]["element"], item["after_order"])
            entry = expected.setdefault(label["code"], {**label, "supporting_change_ids": []})
            entry["supporting_change_ids"].append(item["change_id"])
    r["classification"] = {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": missing(),
        "bond_changes": facts,
        "labels": structural(
            list(expected.values()),
            evidence,
            "Taxonomy labels computed from the explicitly listed partial bond changes.",
        ),
        "coverage": inferred(
            "partial",
            evidence,
            "Only unambiguous net changes visible in the scheme are listed; aromatic reorganization, atom provenance, and coproducts are not fully mapped.",
        ),
        "limitations": [limitation, "No three-model voting was performed."],
        "model_votes": [],
    }


def set_unresolved_classification(r, evidence, limitation):
    r["classification"] = {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": missing(),
        "bond_changes": [],
        "labels": missing(),
        "coverage": inferred("unresolved", evidence, limitation),
        "limitations": [limitation, "No three-model voting was performed."],
        "model_votes": [],
    }


DEAROM_CHANGES = [
    {
        "change_id": "b1",
        "atoms": [
            {"atom_map": None, "element": "C", "source_position": "naphthol C1"},
            {"atom_map": None, "element": "C", "source_position": "terminal allenamide carbon"},
        ],
        "before_order": 0,
        "after_order": 1,
        "change_type": "formed",
    },
    {
        "change_id": "b2",
        "atoms": [
            {"atom_map": None, "element": "C", "source_position": "allenamide central carbon"},
            {"atom_map": None, "element": "H", "source_position": "transferred proton; precise donor mapping unresolved"},
        ],
        "before_order": 0,
        "after_order": 1,
        "change_type": "hydrogen_gain",
        "hydrogen_count": 1,
    },
    {
        "change_id": "b3",
        "atoms": [
            {"atom_map": None, "element": "C", "source_position": "terminal allenamide carbon"},
            {"atom_map": None, "element": "C", "source_position": "allenamide central carbon"},
        ],
        "before_order": 2,
        "after_order": 1,
        "change_type": "order_decreased",
    },
    {
        "change_id": "b4",
        "atoms": [
            {"atom_map": None, "element": "C", "source_position": "naphthol carbonyl carbon"},
            {"atom_map": None, "element": "O", "source_position": "naphthol oxygen"},
        ],
        "before_order": 1,
        "after_order": 2,
        "change_type": "order_increased",
    },
]


def add_dearom_conditions(r, substrate_mid, catalyst_mid, evidence, *, loading=5, temperature="r.t.", hours=None):
    r["participants"].append(participant("naphthol", substrate_mid, "reactant", evidence))
    r["participants"].append(participant("allenamide", r.pop("_allenamide"), "reactant", evidence))
    r["participants"].append(participant("catalyst_system", catalyst_mid, "catalyst", evidence, mol_percent=loading, reference=substrate_mid))
    r["participants"].append(participant("solvent", "m_toluene", "solvent", evidence))
    r["participants"].append(participant("product", r.pop("_product"), "product", evidence))
    r["reaction_scale"] = missing()
    r["solvent_conditions"]["composition"] = f("toluene", ["m_toluene"], evidence)
    r["steps"] = [
        step(
            "reaction",
            evidence,
            [p["participant_id"] for p in r["participants"]],
            operation="combine the depicted substrates with the in situ catalyst system",
            temperature=temperature if isinstance(temperature, tuple) else temperature,
            time=((f"{hours:g} h", {"number": hours, "unit": "h"}) if hours is not None else None),
            atmosphere="nitrogen",
        )
    ]
    r["procedure"]["original_text"] = f(
        {"source_region": r["location"]["raw_value"], "conditions": "table heading, arrow, and footnotes"},
        "table/scheme conditions only; detailed operation is not in the supplied main text",
        evidence,
    )
    set_partial_classification(
        r,
        evidence,
        DEAROM_CHANGES,
        "Dearomatization includes aromatic bond reorganization and proton transfer whose complete atom mapping cannot be established from the main-text drawings alone.",
    )


def add_yield_ee(r, product_mid, evidence, yield_percent, ee_percent):
    r["outcomes"] = [
        outcome(
            "yield",
            "yield",
            [product_mid],
            evidence,
            (f"{yield_percent}%", {"number": yield_percent, "unit": "%"}),
            basis="isolated after flash chromatography",
            yield_type="isolated",
            scope={"table_entry": r["reaction_id"].split(":")[-1]},
        ),
        outcome(
            "ee",
            "selectivity",
            [product_mid],
            evidence,
            (f"{ee_percent}% ee", {"number": ee_percent, "unit": "%"}),
            basis="chiral HPLC",
            scope={"metric": "ee", "major_configuration": "S", "alkene_geometry": "E"},
            role="enantioenriched product",
        ),
    ]


def make_dearom_entry(rid, table, row, substrate, allenamide, product, catalyst, loading, yld, ee, *, categories, temp, hours, evidence):
    unresolved = [
        "Associated SI is mentioned by the article but was not supplied or explicitly linked; substrate quantities, scale, workup, and purification details are unavailable here.",
        "Candidate reactant/product structures and stereochemical encoding require expert chemical review.",
        "Bond classification is partial; aromatic bond reorganization and full proton provenance are not completely atom-mapped.",
    ]
    if table == 1:
        unresolved.append("Reaction duration is not reported in Table 1 or the supplied main-text prose.")
        unresolved.append("The slash catalyst notation gives a catalyst-system loading; individual component charges/speciation are not separately reported.")
    elif row != 10 and catalyst.startswith("m_sys_"):
        unresolved.append("The slash catalyst notation gives a catalyst-system loading; individual component charges/speciation are not separately reported.")
    r = base_reaction(
        rid,
        "specific_experiment",
        categories,
        evidence,
        {"pdf_page": 2 if table == 1 else 3, "table": table, "entry": row},
        "enantioselective dearomatization of β-naphthols with allenamides",
        "C1-selective allylating dearomatization",
        unresolved,
    )
    r["general_scheme_ids"] = [f"{DOC}:general_dearomatization"]
    r["_allenamide"] = allenamide
    r["_product"] = product
    add_dearom_conditions(r, substrate, catalyst, evidence, loading=loading, temperature=temp, hours=hours)
    add_yield_ee(r, product, evidence, yld, ee)
    return r


def main():
    dataset_path = RUN / "dataset.json"
    data = st.read(dataset_path)
    doc = data["documents"][0]
    doc["title"] = "New Chiral BINOL-Based Phosphates for Enantioselective [Au(I)]-Catalyzed Dearomatization of β-Naphthols with Allenamides"
    doc["supplementary_material"] = {
        "status": "mentioned_not_supplied",
        "note": "Page 1 explicitly states that supporting information is available online, but no SI file was supplied or explicitly associated in this session.",
    }
    st.write(dataset_path, data)

    bundle = st.read(RUN / "submitted.json")
    bundle["documents"] = data["documents"]
    molecules = []

    # Generic present-work reaction families.
    molecules += [
        generic_material("m_naphthol_template", "1", "substituted 2-naphthol family 1", [E[1], EI[1], E[3], EI[3]], "2-naphthol core bearing R and R1 substituents", "Generic R/R1 family in Figure 1c and Table 2; no single concrete SMILES represents all members."),
        generic_material("m_allenamide_template", "2", "N-tosyl allenamide family 2", [E[1], EI[1], E[3], EI[3]], "CH2=C=CH-N(R2)Ts", "Generic R2 allenamide family in Figure 1c and Table 2; the terminal allene is not axially stereogenic at its CH2 end."),
        generic_material("m_product_template", "3", "C1-allylated naphthalenone family 3", [E[1], EI[1], E[3], EI[3]], "dearomatized naphthalenone bearing an E-enamide side chain and a stereogenic C1 center", "Generic product family; the paper assigns major S and E products, but substituent variables prevent one concrete SMILES.", chirality="specified", axial=None, assignment={"major_absolute_configuration": "S", "alkene_geometry": "E", "basis": "3ab by optical-rotation comparison; other 3 products by analogy"}),
    ]

    # Concrete substrates and products. Product candidates encode E geometry but leave the S center untagged; the source assignment is retained separately.
    substrate_smiles = {
        "m_1a": ("1a", "1,3-dimethyl-2-naphthol", "Cc1c(O)c(C)c2ccccc2c1"),
        "m_1b": ("1b", "1-methyl-2-naphthol", "Cc1c(O)cc2ccccc2c1"),
        "m_1c": ("1c", "1-benzyl-2-naphthol", "c1ccccc1Cc1c(O)cc2ccccc2c1"),
        "m_1d": ("1d", "1-methyl-3-ethyl-2-naphthol", "Cc1c(O)c(CC)c2ccccc2c1"),
        "m_1e": ("1e", "1-methyl-3-bromo-2-naphthol", "Cc1c(O)c(Br)c2ccccc2c1"),
        "m_1f": ("1f", "1-methyl-3-benzyl-2-naphthol", "Cc1c(O)c(Cc2ccccc2)c2ccccc2c1"),
    }
    for mid, (label, name, smi) in substrate_smiles.items():
        molecules.append(molecule(mid, label, name, [E[2], EI[2], E[3], EI[3]], smiles=smi, chirality="achiral", geometry="not_applicable", representation="not_applicable", identity_note="Position assignment follows the Table 2 R/R1 drawing; candidate structure is not independently source-certified."))

    allenamides = {
        "m_2a": ("2a", "N-phenyl-N-tosyl allenamide", "C=C=CN(c1ccccc1)S(=O)(=O)c1ccc(C)cc1"),
        "m_2b": ("2b", "N-(4-trifluoromethylphenyl)-N-tosyl allenamide", "C=C=CN(c1ccc(C(F)(F)F)cc1)S(=O)(=O)c1ccc(C)cc1"),
        "m_2c": ("2c", "N-methyl-N-tosyl allenamide", "C=C=CN(C)S(=O)(=O)c1ccc(C)cc1"),
    }
    for mid, (label, name, smi) in allenamides.items():
        molecules.append(molecule(mid, label, name, [E[2], EI[2], E[3], EI[3]], smiles=smi, chirality="achiral", geometry="not_applicable", representation="partial", axial=[{"axis": "C=C=C", "assignment": "not stereogenic because one terminus is CH2"}], identity_note="Candidate connectivity follows the allenamide drawings; no source stereodescriptor is reported.", allene=True))

    product_smiles = {
        "m_3aa": ("3aa", "(S,E)-3aa", "O=C1C(C)=Cc2ccccc2C1(C)C/C=C/N(c1ccccc1)S(=O)(=O)c1ccc(C)cc1"),
        "m_3ba": ("3ba", "(S,E)-3ba", "O=C1C=Cc2ccccc2C1(C)C/C=C/N(c1ccccc1)S(=O)(=O)c1ccc(C)cc1"),
        "m_3ca": ("3ca", "(S,E)-3ca", "O=C1C=Cc2ccccc2C1(Cc1ccccc1)C/C=C/N(c1ccccc1)S(=O)(=O)c1ccc(C)cc1"),
        "m_3da": ("3da", "(S,E)-3da", "O=C1C(CC)=Cc2ccccc2C1(C)C/C=C/N(c1ccccc1)S(=O)(=O)c1ccc(C)cc1"),
        "m_3ea": ("3ea", "(S,E)-3ea", "O=C1C(Br)=Cc2ccccc2C1(C)C/C=C/N(c1ccccc1)S(=O)(=O)c1ccc(C)cc1"),
        "m_3fa": ("3fa", "(S,E)-3fa", "O=C1C(Cc3ccccc3)=Cc2ccccc2C1(C)C/C=C/N(c1ccccc1)S(=O)(=O)c1ccc(C)cc1"),
        "m_3ab": ("3ab", "(S,E)-3ab", "O=C1C(C)=Cc2ccccc2C1(C)C/C=C/N(c1ccc(C(F)(F)F)cc1)S(=O)(=O)c1ccc(C)cc1"),
        "m_3ac": ("3ac", "(S,E)-3ac", "O=C1C(C)=Cc2ccccc2C1(C)C/C=C/N(C)S(=O)(=O)c1ccc(C)cc1"),
    }
    for mid, (label, name, smi) in product_smiles.items():
        if mid == "m_3ab":
            assignment = {"absolute_configuration": "S", "basis": "optical rotation comparison with literature", "alkene_geometry": "E"}
        else:
            assignment = {"absolute_configuration": "S", "basis": "assigned by analogy by the authors", "alkene_geometry": "E"}
        molecules.append(molecule(mid, label, name, [E[2], EI[2], E[3], EI[3]], smiles=smi, chirality="specified", geometry="specified", representation="partial", atom_centers=[{"center": "new C1 stereocenter", "configuration": "S", "encoding": "not encoded in candidate SMILES"}], double_bonds=[{"bond": "enamide C=C", "configuration": "E", "encoding": "encoded by slash notation in candidate SMILES"}], assignment=assignment, mixture={"type": "enantioenriched", "composition": "experiment-specific ee; not a fixed registry composition"}, identity_note="The candidate graph encodes the reported E enamide geometry. The S stereocenter is deliberately left untagged because local @/@@ encoding was not independently established from the drawing; representation is partial."))

    molecules.append(molecule("m_toluene", "toluene", "toluene", [E[2], E[3]], smiles="Cc1ccccc1", chirality="achiral", geometry="not_applicable", representation="not_applicable"))

    # Catalyst systems retain source slash notation as named in-situ mixtures; individual speciation/loadings are unresolved.
    table1_systems = {
        "m_sys_l1_c1": ("L1/C1", "PPh3AuCl/AgC1 in situ catalyst system"),
        "m_sys_l1_c2": ("L1/C2", "PPh3AuCl/AgC2 in situ catalyst system"),
        "m_sys_l1_c3": ("L1/C3", "PPh3AuCl/AgC3 in situ catalyst system"),
        "m_sys_l2_c1": ("L2/C1", "L2(AuCl)2/AgC1 in situ catalyst system"),
        "m_sys_l3_c1": ("L3/C1", "L3AuCl/AgC1 in situ catalyst system"),
        "m_sys_l4_c1": ("L4/C1", "JohnPhosAuCl/AgC1 in situ catalyst system"),
        "m_sys_l5_c1": ("L5/C1", "IPrAuCl/AgC1 in situ catalyst system"),
        "m_sys_l6_c1": ("L6/C1", "Ad-NHCAuCl/AgC1 in situ catalyst system"),
    }
    for mid, (label, name) in table1_systems.items():
        molecules.append(generic_material(mid, label, name, [E[2], EI[2]], name, "Named two-component in situ catalyst system; the main text does not separately report component charge/speciation.", chirality="partially_specified", axial=[{"component": label.split("/")[-1], "type": "chiral BINOL phosphate", "descriptor": "not printed"}]))

    for code, axial_desc in [("c4", "Ra inherited from precursor 4"), ("c5", "Ra inherited from precursor 4"), ("c6", "Ra inherited from precursor 4")]:
        molecules.append(generic_material(f"m_sys_johnphos_ag{code}", f"JohnPhosAuCl/Ag{code.upper()}", f"JohnPhosAuCl/Ag{code.upper()} in situ catalyst system", [E[3], EI[3]], f"JohnPhosAuCl plus Ag{code.upper()}", "Named in situ catalyst system; individual charge/speciation is not distinguished by the common 5 mol% arrow label.", chirality="partially_specified", axial=[{"component": f"Ag{code.upper()}", "configuration": axial_desc}]))

    # Catalyst synthesis molecules and common reagents.
    molecules += [
        molecule("m_4", "4", "(Ra)-(2,2'-dimethoxy-[1,1'-binaphthalene]-3,3'-diyl)diboronic acid", [E[2], E[3], EI[3]], chirality="specified", geometry="not_applicable", representation="not_encoded", axial=[{"axis": "binaphthyl", "configuration": "Ra"}], assignment={"source": "explicit (Ra) name in prose"}, identity_note="Large source-drawn structure retained without a concrete SMILES to avoid unsafe axial-stereo encoding."),
        molecule("m_5", "5", "9,10-dibromoanthracene", [E[2], E[3], EI[3]], smiles="Brc1c2ccccc2c(Br)c2ccccc12", chirality="achiral", geometry="not_applicable", representation="not_applicable"),
        molecule("m_6", "6", "dibromo anthracenyl-substituted dimethoxy BINOL adduct 6", [E[2], E[3], EI[3]], chirality="specified", geometry="not_applicable", representation="not_encoded", axial=[{"axis": "binaphthyl", "configuration": "Ra retained from 4"}], identity_note="Exact large structure is depicted in Scheme 2 but not transcribed to a concrete SMILES; axial retention is a structural-continuity inference."),
        generic_material("m_7", "7", "arylboronic acid family 7", [E[2], E[3], EI[3]], "Ar-B(OH)2 where Ar bears X = 3,5-(CF3)2, H, or 3,5-Me2", "Scheme 2b uses a variable arylboronic-acid family."),
    ]
    for suffix, x in [("a", "3,5-(CF3)2"), ("b", "H"), ("c", "3,5-Me2")]:
        molecules.append(molecule(f"m_8{suffix}", f"8{suffix}", f"diaryl-anthracenyl BINOL diol 8{suffix} (X={x})", [E[2], E[3], EI[3]], chirality="specified", geometry="not_applicable", representation="not_encoded", axial=[{"axis": "binaphthyl", "configuration": "Ra retained from 4"}], identity_note=f"Scheme 2 defines X={x}; large structure retained without concrete SMILES."))
        molecules.append(molecule(f"m_9{suffix}", f"9{suffix}", f"chiral BINOL phosphoric acid 9{suffix} (X={x})", [E[2], E[3], EI[3]], chirality="specified", geometry="not_applicable", representation="not_encoded", axial=[{"axis": "binaphthyl", "configuration": "Ra retained from 4"}], identity_note=f"Scheme 3 phosphoric acid with X={x}; large structure retained without concrete SMILES."))
        molecules.append(molecule(f"m_agc{4 + 'abc'.index(suffix)}", f"AgC{4 + 'abc'.index(suffix)}", f"silver phosphate AgC{4 + 'abc'.index(suffix)} (X={x})", [E[2], E[3], EI[3]], chirality="specified", geometry="not_applicable", representation="not_encoded", axial=[{"axis": "binaphthyl", "configuration": "Ra retained from 4"}], identity_note=f"Scheme 3 silver salt with X={x}; large structure retained without concrete SMILES."))

    simple = [
        ("m_pdpph3_4", "Pd(PPh3)4", "tetrakis(triphenylphosphine)palladium(0)", None),
        ("m_k2co3", "K2CO3", "potassium carbonate", "[K+].[K+].[O-]C([O-])=O"),
        ("m_thf", "THF", "tetrahydrofuran", "C1CCOC1"),
        ("m_water", "H2O", "water", "O"),
        ("m_pdoac2", "Pd(OAc)2", "palladium(II) acetate", "CC(=O)[O-].CC(=O)[O-].[Pd+2]"),
        ("m_otol3p", "(o-Tol)3P", "tri(o-tolyl)phosphine", "Cc1ccccc1P(c1ccccc1C)c1ccccc1C"),
        ("m_dme", "DME", "1,2-dimethoxyethane", "COCCOC"),
        ("m_bbr3", "BBr3", "boron tribromide", "B(Br)(Br)Br"),
        ("m_ch2cl2", "CH2Cl2", "dichloromethane", "ClCCl"),
        ("m_pocl3", "POCl3", "phosphoryl chloride", "O=P(Cl)(Cl)Cl"),
        ("m_pyridine", "py", "pyridine", "n1ccccc1"),
        ("m_hcl", "HCl", "hydrogen chloride", "Cl"),
        ("m_ag2co3", "Ag2CO3", "silver carbonate", "[Ag+].[Ag+].[O-]C([O-])=O"),
        ("m_johnphosaucl", "JohnPhosAuCl", "JohnPhos gold(I) chloride", None),
        ("m_agc2_r", "(R)-AgC2", "(R)-silver phosphate AgC2", None),
        ("m_d8_toluene", "d8-toluene", "toluene-d8", "[2H]C([2H])([2H])c1c([2H])c([2H])c([2H])c([2H])c1[2H]"),
        ("m_johnphosauc4", "JohnPhosAuC4", "chiral JohnPhosAuC4 adduct", None),
    ]
    for mid, label, name, smi in simple:
        kwargs = dict(chirality="achiral", geometry="not_applicable", representation="not_applicable")
        if mid in {"m_agc2_r"}:
            kwargs = dict(chirality="specified", geometry="not_applicable", representation="not_encoded", axial=[{"configuration": "R as printed in Table 2 footnote e"}])
        elif mid in {"m_johnphosaucl", "m_johnphosauc4", "m_pdpph3_4"}:
            kwargs = dict(chirality="unresolved" if mid == "m_johnphosauc4" else "achiral", geometry="not_applicable", representation="unresolved")
        molecules.append(molecule(mid, label, name, [E[2], E[3], EI[3]], smiles=smi, identity_note="Concrete SMILES omitted for complex coordination species." if smi is None else None, **kwargs))

    bundle["molecules"] = molecules
    reactions = []

    # Present-work generic transformation, Figure 1c / Table 2 heading.
    rg = base_reaction(
        "general_dearomatization",
        "generic_scheme",
        ["other"],
        [E[1], EI[1], E[3], EI[3]],
        {"pdf_pages": [1, 3], "regions": ["Figure 1c", "Table 2 reaction heading"]},
        "enantioselective [Au(I)]-catalyzed dearomatization of β-naphthols with allenamides",
        "C1-selective allylating dearomatization",
        ["Generic R/R1/R2 scheme; no single experiment scale or full procedure is asserted.", "Complete atom mapping is unresolved."],
    )
    rg["participants"] = [
        participant("naphthol_family", "m_naphthol_template", "reactant", rg["evidence_ids"]),
        participant("allenamide_family", "m_allenamide_template", "reactant", rg["evidence_ids"]),
        participant("product_family", "m_product_template", "product", rg["evidence_ids"]),
    ]
    rg["steps"] = [step("generic_reaction", rg["evidence_ids"], ["naphthol_family", "allenamide_family", "product_family"], operation="gold/chiral-counterion-catalyzed dearomatization")]
    set_partial_classification(rg, rg["evidence_ids"], DEAROM_CHANGES, "Generic source scheme; substituent instances and full aromatic bond changes are not completely mapped.")
    reactions.append(rg)

    # Table 1, eight independent optimization rows.
    t1 = [
        (1, "m_sys_l1_c1", 5, 79, 75),
        (2, "m_sys_l1_c2", 5, 53, 72),
        (3, "m_sys_l1_c3", 5, 70, 62),
        (4, "m_sys_l2_c1", 2.5, 56, 80),
        (5, "m_sys_l3_c1", 5, 56, 80),
        (6, "m_sys_l4_c1", 5, 92, 83),
        (7, "m_sys_l5_c1", 5, 61, 81),
        (8, "m_sys_l6_c1", 5, 65, 80),
    ]
    for row, cat, loading, yld, ee in t1:
        reactions.append(make_dearom_entry(f"table1_entry{row}", 1, row, "m_1a", "m_2a", "m_3aa", cat, loading, yld, ee, categories=["optimization"], temp="r.t.", hours=None, evidence=[E[2], EI[2]]))

    # Scheme 2a: one concrete, source-reported Suzuki coupling.
    rs2a = base_reaction(
        "scheme2a_compound6",
        "specific_experiment",
        ["other"],
        [E[2], E[3], EI[3]],
        {"pdf_page": 3, "scheme": "2a"},
        "chemoselective Pd(0)-mediated coupling",
        "double Suzuki-Miyaura cross-coupling",
        ["Reaction time, atmosphere, workup, purification, and yield measurement basis are not reported in the supplied main text.", "Exact large structures and complete atom mapping require review; associated SI was not supplied."],
    )
    rs2a["participants"] = [
        participant("compound4", "m_4", "reactant", rs2a["evidence_ids"]),
        participant("compound5", "m_5", "reactant", rs2a["evidence_ids"]),
        participant("pd_catalyst", "m_pdpph3_4", "catalyst", rs2a["evidence_ids"], mol_percent=5, reference="m_4"),
        participant("base", "m_k2co3", "reagent", rs2a["evidence_ids"], equivalents=3, reference="m_4"),
        participant("thf", "m_thf", "solvent", rs2a["evidence_ids"]),
        participant("water", "m_water", "solvent", rs2a["evidence_ids"]),
        participant("compound6", "m_6", "product", rs2a["evidence_ids"]),
    ]
    rs2a["solvent_conditions"]["composition"] = f("THF/H2O", ["m_thf", "m_water"], rs2a["evidence_ids"])
    rs2a["solvent_conditions"]["mixture_ratio"] = f("20:1", {"components": ["THF", "H2O"], "ratio": "20:1"}, rs2a["evidence_ids"])
    rs2a["steps"] = [step("coupling", rs2a["evidence_ids"], [p["participant_id"] for p in rs2a["participants"]], operation="Pd-catalyzed coupling", temperature="reflux")]
    rs2a["outcomes"] = [outcome("yield6", "yield", ["m_6"], rs2a["evidence_ids"], ("65%", {"number": 65, "unit": "%"}), scope={"from": "4 and 5", "to": "6"})]
    set_partial_classification(rs2a, rs2a["evidence_ids"], [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C", "source_position": "BINOL C3"}, {"atom_map": None, "element": "C", "source_position": "anthracene C9"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C", "source_position": "BINOL C3'"}, {"atom_map": None, "element": "C", "source_position": "second anthracene C9"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
    ], "Only the two source-depicted C-C couplings are classified; leaving groups/coproducts are not fully mapped.")
    reactions.append(rs2a)

    # Scheme 2b family: range is across 8a-c, so retain one generic scheme rather than inventing per-compound yields.
    rs2b = base_reaction(
        "scheme2b_family8",
        "generic_scheme",
        ["other"],
        [E[2], E[3], EI[3]],
        {"pdf_page": 3, "scheme": "2b"},
        "coupling/deprotection sequence to 8a-c",
        "Suzuki-Miyaura arylation followed by BBr3 demethylation",
        ["The source reports only a 38-52% overall range across 8a-c; individual yields cannot be assigned.", "Reaction times, atmosphere, workup, and purification are not reported in the supplied main text/Scheme 2b."],
    )
    rs2b["participants"] = [
        participant("compound6", "m_6", "reactant", rs2b["evidence_ids"]),
        participant("arylboronic_acids", "m_7", "reactant", rs2b["evidence_ids"]),
        participant("pd_catalyst", "m_pdoac2", "catalyst", rs2b["evidence_ids"], mol_percent=10, reference="m_6"),
        participant("ligand", "m_otol3p", "ligand", rs2b["evidence_ids"], mol_percent=40, reference="m_6"),
        participant("base", "m_k2co3", "reagent", rs2b["evidence_ids"]),
        participant("dme", "m_dme", "solvent", rs2b["evidence_ids"]),
        participant("bbr3", "m_bbr3", "reagent", rs2b["evidence_ids"]),
        participant("dcm", "m_ch2cl2", "solvent", rs2b["evidence_ids"]),
        participant("product8a", "m_8a", "product", rs2b["evidence_ids"]),
        participant("product8b", "m_8b", "product", rs2b["evidence_ids"]),
        participant("product8c", "m_8c", "product", rs2b["evidence_ids"]),
    ]
    rs2b["steps"] = [
        step("arylation", rs2b["evidence_ids"], ["compound6", "arylboronic_acids", "pd_catalyst", "ligand", "base", "dme"], operation="Pd(OAc)2/(o-Tol)3P coupling with K2CO3 (2 M) in DME", temperature=("80 °C", {"number": 80, "unit": "°C"})),
        step("demethylation", rs2b["evidence_ids"], ["bbr3", "dcm", "product8a", "product8b", "product8c"], operation="BBr3-promoted cleavage in CH2Cl2", temperature=("0 °C", {"number": 0, "unit": "°C"})),
    ]
    rs2b["outcomes"] = [outcome("overall_yield_range", "yield", ["m_8a", "m_8b", "m_8c"], rs2b["evidence_ids"], ("38-52% overall", {"minimum": 38, "maximum": 52, "unit": "%", "raw_text": "38-52%"}), scope={"span": "6/7 to 8a-c", "assignment": "range across family; individual values unreported"})]
    set_partial_classification(rs2b, rs2b["evidence_ids"], [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C", "source_position": "anthracene C-Br site"}, {"atom_map": None, "element": "C", "source_position": "arylboronic acid ipso C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C", "source_position": "second anthracene C-Br site"}, {"atom_map": None, "element": "C", "source_position": "second arylboronic acid ipso C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
    ], "Family scheme; only the two new aryl C-C connections are classified. Demethylation bond cleavages and coproducts are not completely mapped.")
    reactions.append(rs2b)

    # Scheme 3 family: phosphorylation/hydrolysis/metal exchange.
    rs3 = base_reaction(
        "scheme3_family_agc4_6",
        "generic_scheme",
        ["other"],
        [E[2], E[3], EI[3]],
        {"pdf_page": 3, "scheme": "3"},
        "phosphorylation/hydrolysis/metal exchange sequence",
        "BINOL phosphoric acid formation and silver-salt exchange",
        ["The source reports only a 91-94% range for AgC4-6 and does not assign individual silver-salt yields.", "Intermediate phosphoric-acid yields are described only as acceptable overall yield; exact values, workup, and purification are unavailable without SI."],
    )
    rs3["participants"] = [
        *[participant(f"start8{x}", f"m_8{x}", "reactant", rs3["evidence_ids"]) for x in "abc"],
        participant("pocl3", "m_pocl3", "reagent", rs3["evidence_ids"]),
        participant("pyridine", "m_pyridine", "solvent", rs3["evidence_ids"]),
        participant("water", "m_water", "reagent", rs3["evidence_ids"]),
        participant("hcl", "m_hcl", "reagent", rs3["evidence_ids"]),
        participant("ag2co3", "m_ag2co3", "reagent", rs3["evidence_ids"]),
        participant("dcm", "m_ch2cl2", "solvent", rs3["evidence_ids"]),
        *[participant(f"acid9{x}", f"m_9{x}", "product", rs3["evidence_ids"]) for x in "abc"],
        *[participant(f"agc{n}", f"m_agc{n}", "product", rs3["evidence_ids"]) for n in (4, 5, 6)],
    ]
    rs3["steps"] = [
        step("phosphorylation", rs3["evidence_ids"], ["start8a", "start8b", "start8c", "pocl3", "pyridine"], operation="POCl3 in pyridine", temperature=("100 °C", {"number": 100, "unit": "°C"})),
        step("hydrolysis", rs3["evidence_ids"], ["water", "acid9a", "acid9b", "acid9c"], operation="H2O treatment", temperature=("100 °C", {"number": 100, "unit": "°C"})),
        step("acid_treatment", rs3["evidence_ids"], ["hcl", "acid9a", "acid9b", "acid9c"], operation="HCl treatment", time=("6 h", {"number": 6, "unit": "h"})),
        step("silver_exchange", rs3["evidence_ids"], ["acid9a", "acid9b", "acid9c", "ag2co3", "dcm", "water", "agc4", "agc5", "agc6"], operation="Ag2CO3 treatment in DCM/H2O"),
    ]
    rs3["outcomes"] = [outcome("silver_salt_yield_range", "yield", ["m_agc4", "m_agc5", "m_agc6"], rs3["evidence_ids"], ("91-94%", {"minimum": 91, "maximum": 94, "unit": "%", "raw_text": "91-94%"}), scope={"targets": "AgC4-6", "assignment": "range across family; individual values unreported"})]
    set_partial_classification(rs3, rs3["evidence_ids"], [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "O", "source_position": "first BINOL oxygen"}, {"atom_map": None, "element": "P", "source_position": "phosphate phosphorus"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "O", "source_position": "second BINOL oxygen"}, {"atom_map": None, "element": "P", "source_position": "phosphate phosphorus"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
    ], "Only the cyclic phosphate O-P bonds are classified; hydrolysis and ionic Ag exchange are not fully atom-mapped.")
    reactions.append(rs3)

    # Table 2, eleven independent rows. Row 10 has explicitly separate 5 mol% catalyst components.
    t2 = [
        (1, "m_1a", "m_2a", "m_3aa", "m_sys_johnphos_agc4", 50, 87, ["optimization"]),
        (2, "m_1a", "m_2a", "m_3aa", "m_sys_johnphos_agc5", 75, 86, ["optimization"]),
        (3, "m_1a", "m_2a", "m_3aa", "m_sys_johnphos_agc6", 98, 90, ["optimization"]),
        (4, "m_1a", "m_2a", "m_3aa", "m_9b", 53, 81, ["optimization", "control"]),
        (5, "m_1b", "m_2a", "m_3ba", "m_sys_johnphos_agc6", 60, 59, ["substrate_scope"]),
        (6, "m_1c", "m_2a", "m_3ca", "m_sys_johnphos_agc6", 47, 55, ["substrate_scope"]),
        (7, "m_1d", "m_2a", "m_3da", "m_sys_johnphos_agc6", 76, 83, ["substrate_scope"]),
        (8, "m_1e", "m_2a", "m_3ea", "m_sys_johnphos_agc6", 24, 95, ["substrate_scope"]),
        (9, "m_1f", "m_2a", "m_3fa", "m_sys_johnphos_agc6", 49, 67, ["substrate_scope"]),
        (10, "m_1a", "m_2b", "m_3ab", "explicit_row10", 40, 43, ["substrate_scope"]),
        (11, "m_1a", "m_2c", "m_3ac", "m_sys_johnphos_agc6", 50, 5, ["substrate_scope"]),
    ]
    for row, sub, allene, prod, cat, yld, ee, cats in t2:
        if row == 10:
            r = base_reaction(
                f"table2_entry{row}", "specific_experiment", cats, [E[3], EI[3]], {"pdf_page": 3, "table": 2, "entry": row}, "enantioselective dearomatization of β-naphthols with allenamides", "C1-selective allylating dearomatization", ["Associated SI was not supplied; substrate quantities, scale, workup, and purification details are unavailable.", "Candidate structures and partial bond classification require expert review."],
            )
            r["general_scheme_ids"] = [f"{DOC}:general_dearomatization"]
            r["participants"] = [
                participant("naphthol", sub, "reactant", r["evidence_ids"]),
                participant("allenamide", allene, "reactant", r["evidence_ids"]),
                participant("johnphosaucl", "m_johnphosaucl", "catalyst", r["evidence_ids"], mol_percent=5, reference=sub),
                participant("agc2", "m_agc2_r", "reagent", r["evidence_ids"], mol_percent=5, reference=sub),
                participant("solvent", "m_toluene", "solvent", r["evidence_ids"]),
                participant("product", prod, "product", r["evidence_ids"]),
            ]
            r["solvent_conditions"]["composition"] = f("toluene", ["m_toluene"], r["evidence_ids"])
            r["steps"] = [step("reaction", r["evidence_ids"], [p["participant_id"] for p in r["participants"]], operation="in situ catalyst formation and dearomatization", temperature=("0 °C", {"number": 0, "unit": "°C"}), time=("16 h", {"number": 16, "unit": "h"}), atmosphere="nitrogen")]
            set_partial_classification(r, r["evidence_ids"], DEAROM_CHANGES, "Dearomatization aromatic reorganization and complete proton provenance are not fully mapped.")
            add_yield_ee(r, prod, r["evidence_ids"], yld, ee)
        else:
            r = make_dearom_entry(f"table2_entry{row}", 2, row, sub, allene, prod, cat, 5, yld, ee, categories=cats, temp=("0 °C", {"number": 0, "unit": "°C"}), hours=16, evidence=[E[3], EI[3]])
        reactions.append(r)

    # Actual 31P-NMR mixing experiment, not a synthesis/yield record.
    rnmr = base_reaction(
        "mechanistic_31p_nmr_mixing",
        "specific_experiment",
        ["mechanistic_study"],
        [E[3], EI[3]],
        {"pdf_page": 3, "region": "right-column final paragraph"},
        "halide-scavenger performance monitored by 31P NMR",
        "in situ catalyst/adduct formation study",
        ["Amounts, mixing time, temperature, and full speciation are not reported.", "The new signals are attributed by the authors to JohnPhosAuC4; no covalent bond formation is established, so reaction classification remains unresolved."],
    )
    rnmr["participants"] = [
        participant("johnphosaucl", "m_johnphosaucl", "reactant", rnmr["evidence_ids"]),
        participant("agc4", "m_agc4", "reactant", rnmr["evidence_ids"]),
        participant("solvent", "m_d8_toluene", "solvent", rnmr["evidence_ids"]),
        participant("adduct", "m_johnphosauc4", "product", rnmr["evidence_ids"]),
    ]
    rnmr["solvent_conditions"]["composition"] = f("d8-toluene", ["m_d8_toluene"], rnmr["evidence_ids"])
    rnmr["steps"] = [step("mixing_and_nmr", rnmr["evidence_ids"], [p["participant_id"] for p in rnmr["participants"]], operation="mix JohnPhosAuCl and AgC4 and monitor by 31P NMR")]
    rnmr["outcomes"] = [outcome("new_31p_signals", "other", ["m_johnphosauc4"], rnmr["evidence_ids"], ("54.69 ppm and 6.90 ppm", {"signals_ppm": [54.69, 6.90]}), basis="31P NMR in d8-toluene", scope={"author_attribution": "catalytically active chiral JohnPhosAuC4 adduct"}, role="author-attributed in situ species")]
    set_unresolved_classification(rnmr, rnmr["evidence_ids"], "NMR signal changes establish a new species but not a source-verified covalent bond graph.")
    reactions.append(rnmr)

    # Eq. 1 is a proposed mechanistic model, not a mechanistic experiment.
    rmech = base_reaction(
        "eq1_proposed_outer_sphere_mechanism",
        "proposed_mechanism",
        ["mechanistic_study"],
        [E[3], EI[3], E[4], EI[4]],
        {"pdf_pages": [3, 4], "regions": ["page 3 mechanistic discussion", "page 4 Eq. 1"]},
        "tentative enantiodiscriminating outer-sphere approach",
        "proposed ion-pair/bifunctional gold-catalysis mechanism",
        ["This is an author-proposed mechanism, not a directly validated experiment.", "Elementary atom transfers and catalyst coordination are incomplete; overall net classification is partial."],
    )
    rmech["general_scheme_ids"] = [f"{DOC}:general_dearomatization"]
    rmech["participants"] = [
        participant("naphthol_family", "m_naphthol_template", "reactant", rmech["evidence_ids"]),
        participant("allenamide_family", "m_allenamide_template", "reactant", rmech["evidence_ids"]),
        participant("product_family", "m_product_template", "product", rmech["evidence_ids"]),
        participant("catalyst", "m_sys_johnphos_agc6", "catalyst", rmech["evidence_ids"]),
    ]
    rmech["steps"] = [step("proposed_outer_sphere_step", rmech["evidence_ids"], [p["participant_id"] for p in rmech["participants"]], operation="proposed chiral-phosphate activation of naphthol and attack on an in situ cationic allyl-gold species")]
    set_partial_classification(rmech, rmech["evidence_ids"], DEAROM_CHANGES, "Overall proposed net transformation only; elementary coordination, proton transfer, and aromatic bond changes are not fully mapped.")
    reactions.append(rmech)

    bundle["reactions"] = reactions
    bundle["source_context"] = [
        {"context_id": f"{DOC}:figure1a_context", "evidence_ids": [E[1], EI[1]], "kind": "non_experimental_generic_context", "description": "Figure 1a surveys generic oxidative and C-C-forming naphthol dearomatization strategies; it is not a current-paper experiment."},
        {"context_id": f"{DOC}:figure1b_bioactive_examples", "evidence_ids": [E[1], EI[1]], "kind": "non_synthetic_figure", "description": "Bioactive structures are examples, not reaction experiments."},
        {"context_id": f"{DOC}:scheme1_counterion_structures", "evidence_ids": [E[2], EI[2]], "kind": "compound_structure_panel", "description": "Scheme 1 depicts second-generation counterions C4-C6 without a reaction arrow."},
        {"context_id": f"{DOC}:starting_31p_signals", "evidence_ids": [E[3], EI[3]], "kind": "characterization", "description": "Before mixing, JohnPhosAuCl and AgC4 are reported at 59.18 and 11.41 ppm, respectively; these are characterization observations, not separate synthesis records."},
        {"context_id": f"{DOC}:si_statement", "evidence_ids": [E[1], EI[1]], "kind": "source_availability", "description": "Supporting information is stated to be available online, but no SI file was supplied or explicitly associated in this session."},
        {"context_id": f"{DOC}:pages4_5_nonreaction", "evidence_ids": [E[4], EI[4], E[5], EI[5]], "kind": "coverage_context", "description": "After Eq. 1, page 4 contains conclusion/acknowledgements/references; page 5 contains references and manuscript dates only."},
    ]

    st.write(RUN / "submitted.json", bundle)

    # Raw model artifacts: page inventory, source-mention hints, and graph diagnostics.
    inventory = {
        "document_id": DOC,
        "semantic_accuracy_verified": False,
        "pages": {
            "1": ["Figure 1a generic literature context", "Figure 1b non-synthetic bioactive examples", "Figure 1c present-work generic transformation/mechanistic concept", "SI availability statement"],
            "2": ["Table 1 rows 1-8", "Scheme 1 counterion structures without reaction arrow", "prose describing catalyst synthesis and Table 2 setup"],
            "3": ["Scheme 2a", "Scheme 2b", "Scheme 3", "Table 2 rows 1-11", "31P-NMR mixing study", "mechanistic discussion leading to Eq. 1"],
            "4": ["Eq. 1 proposed mechanism", "conclusion and references"],
            "5": ["references and manuscript dates only"],
        },
        "record_accounting": {
            "table1_rows": 8,
            "table2_rows": 11,
            "other_specific_experiments": 2,
            "generic_schemes": 3,
            "proposed_mechanisms": 1,
        },
    }
    st.write(RUN / "raw_model/page1_5_visual_inventory.json", inventory)
    st.write(
        RUN / "raw_model/source_mentions.json",
        st.source_mentions(data, ["1 a", "2 a", "3 aa", "3 ab", "JohnPhosAuCl", "AgC4"], context=180),
    )
    diagnostics = {}
    for m in molecules:
        smi = st.value(m["smiles"])
        if smi:
            diagnostics[m["molecule_id"]] = st.stereo_inspect(smi)
    st.write(RUN / "raw_model/stereo_diagnostics.json", diagnostics)

    by_id = {r["reaction_id"]: r for r in reactions}
    t1_ids = [f"{DOC}:table1_entry{i}" for i in range(1, 9)]
    t2_ids = [f"{DOC}:table2_entry{i}" for i in range(1, 12)]
    coverage = {
        "document_id": DOC,
        "semantic_accuracy_verified": False,
        "pages": [
            {
                "page": 1,
                "state": "attempted",
                "regions": [
                    {"region": "Figure 1a", "finding": "Generic literature-context dearomatization arrows; retained as source_context, not current-paper experiments.", "evidence_ids": [E[1], EI[1]]},
                    {"region": "Figure 1b", "finding": "Bioactive structures only; no synthesis experiment.", "evidence_ids": [E[1], EI[1]]},
                    {"region": "Figure 1c", "finding": "Present-work generic dearomatization and proposed activation concept; mapped to the generic scheme.", "evidence_ids": [E[1], EI[1]]},
                    {"region": "lower-left SI statement", "finding": "SI is stated to be available online but was not supplied/associated.", "evidence_ids": [E[1], EI[1]]},
                ],
                "record_ids": [f"{DOC}:general_dearomatization"],
                "unresolved": ["No SI file was supplied; Figure 1c uses generic R/R1/R2 groups."],
            },
            {
                "page": 2,
                "state": "attempted",
                "regions": [
                    {"region": "Table 1 including reaction arrow, rows 1-8, and footnotes a-d", "finding": "Eight independent optimization experiments with isolated yields and chiral-HPLC ee values.", "evidence_ids": [E[2], EI[2]]},
                    {"region": "Scheme 1", "finding": "Structures C4-C6; no reaction arrow, retained as structural context.", "evidence_ids": [E[2], EI[2]]},
                    {"region": "right-column catalyst-synthesis prose", "finding": "Textual conditions/yields for Scheme 2a, 2b and Scheme 3, linked to the corresponding page-3 scheme records.", "evidence_ids": [E[2], EI[2]]},
                ],
                "record_ids": t1_ids + [f"{DOC}:scheme2a_compound6", f"{DOC}:scheme2b_family8", f"{DOC}:scheme3_family_agc4_6"],
                "unresolved": ["Table 1 duration and individual slash-catalyst component charges are not reported.", "Catalyst-synthesis experimental details refer to SI, which was not supplied."],
            },
            {
                "page": 3,
                "state": "attempted",
                "regions": [
                    {"region": "Scheme 2a", "finding": "Specific synthesis of 6, 65%.", "evidence_ids": [E[3], EI[3]]},
                    {"region": "Scheme 2b", "finding": "Family coupling/deprotection sequence to 8a-c, overall range 38-52%; not separable into individual yields.", "evidence_ids": [E[3], EI[3]]},
                    {"region": "Scheme 3", "finding": "Family phosphorylation/hydrolysis/Ag exchange to AgC4-6; final range 91-94% not individually assigned.", "evidence_ids": [E[3], EI[3]]},
                    {"region": "Table 2 rows 1-11 and footnotes a-e", "finding": "Eleven independent screening/scope experiments with isolated yield and chiral-HPLC ee.", "evidence_ids": [E[3], EI[3]]},
                    {"region": "right-column final paragraph", "finding": "Actual 31P-NMR mixing study and author attribution of new signals to JohnPhosAuC4.", "evidence_ids": [E[3], EI[3]]},
                    {"region": "mechanistic discussion", "finding": "Author-proposed outer-sphere mechanism continued in Eq. 1 on page 4.", "evidence_ids": [E[3], EI[3]]},
                ],
                "record_ids": [f"{DOC}:scheme2a_compound6", f"{DOC}:scheme2b_family8", f"{DOC}:scheme3_family_agc4_6"] + t2_ids + [f"{DOC}:mechanistic_31p_nmr_mixing", f"{DOC}:eq1_proposed_outer_sphere_mechanism"],
                "unresolved": ["Individual yields within the Scheme 2b and Scheme 3 ranges are not reported.", "NMR mixing amounts, temperature, time, and complete adduct structure/speciation are not reported."],
            },
            {
                "page": 4,
                "state": "attempted",
                "regions": [
                    {"region": "upper-left Eq. 1", "finding": "Tentative ion-pair/outer-sphere mechanism; classified as proposed_mechanism, not an experiment.", "evidence_ids": [E[4], EI[4]]},
                    {"region": "conclusion, acknowledgements, conflict statement, references", "finding": "No additional current-paper experiment or reaction table.", "evidence_ids": [E[4], EI[4]]},
                ],
                "record_ids": [f"{DOC}:eq1_proposed_outer_sphere_mechanism"],
                "unresolved": ["Eq. 1 is qualitative and does not establish complete elementary atom mapping."],
            },
            {
                "page": 5,
                "state": "attempted",
                "regions": [
                    {"region": "full page", "finding": "References [15]-[21] and manuscript received/revised dates only; no current-paper experiment, table, reaction arrow, or characterization table.", "evidence_ids": [E[5], EI[5]]},
                ],
                "record_ids": [],
                "unresolved": [],
            },
        ],
    }
    assert set(by_id) == set(sum([p["record_ids"] for p in coverage["pages"]], [])) | {
        f"{DOC}:general_dearomatization"
    }
    st.write(RUN / "coverage/pages_0001_0005.json", coverage)
    print(json.dumps({"molecules": len(molecules), "reactions": len(reactions), "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in reactions)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
