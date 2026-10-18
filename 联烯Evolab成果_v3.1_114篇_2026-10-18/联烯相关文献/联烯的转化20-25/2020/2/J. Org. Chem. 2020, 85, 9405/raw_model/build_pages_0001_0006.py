#!/usr/bin/env python3
"""Build source-linked v3.1 records for PDF pages 1-6 of doc_106b2f72bb5b02d2."""
from __future__ import annotations

import copy
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
import standard_records as st

DOC = "doc_106b2f72bb5b02d2"
DATASET = ROOT / "runs/current/dataset.json"
OUT_RAW = ROOT / "runs/current/raw_model/extracted_bundle_pages_0001_0006.json"
OUT_SUBMITTED = ROOT / "runs/current/submitted.json"
OUT_COVERAGE = ROOT / "runs/current/coverage/pages_0001_0006.json"


def ev(page: int) -> list[str]:
    return [f"{DOC}:page{page}:text", f"{DOC}:page{page}:image"]


def evs(*pages: int) -> list[str]:
    return [item for page in pages for item in ev(page)]


def f(raw=None, normalized=None, evidence=(), method="original_extraction", derivation=None,
      status=None):
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation, status=status)


def q(number, unit, raw=None, qualifier=None):
    value = {"number": number, "unit": unit}
    if raw is not None:
        value["raw_text"] = raw
    if qualifier:
        value["qualifier"] = qualifier
    return value


dataset = json.loads(DATASET.read_text())
bundle = st.empty(dataset["documents"])
molecules: dict[str, dict] = {}
reactions: list[dict] = []


def add_molecule(key: str, label: str, name: str, evidence: list[str], smiles: str | None = None,
                 aliases=None, stereo_note: str | None = None, mixture=None) -> str:
    mid = f"{DOC}:mol_{key}"
    if mid in molecules:
        occurrence = molecules[mid]["occurrences"][0]
        occurrence["evidence_ids"] = list(dict.fromkeys(occurrence["evidence_ids"] + evidence))
        return mid
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    m["names"]["original"] = f(name, name, evidence)
    m["names"]["english"] = f(name, name, evidence)
    m["names"]["abbreviations"] = f(aliases or [], aliases or [], evidence) if aliases else f()
    m["names"]["aliases"] = f(aliases or [], aliases or [], evidence) if aliases else f()
    m["cas_numbers"] = f()
    if smiles:
        diag = st.stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = f(
            f"Source structure/name for {label}", smiles, evidence,
            method="structure_derivation",
            derivation="Candidate SMILES transcribed from the source drawing or explicit systematic name; RDKit parsing is diagnostic only.",
        )
        atom_potential = [x for x in diag["potential_stereo"] if x["type"] == "Atom_Tetrahedral"]
        bond_potential = [x for x in diag["potential_stereo"] if x["type"] == "Bond_Double"]
        m["stereochemistry"]["atom_centers"] = f(
            "standard stereo diagnostic", diag["atom_centers"], evidence,
            method="structure_derivation", derivation="RDKit zero-based atom indices for the candidate SMILES.")
        m["stereochemistry"]["double_bonds"] = f(
            "standard stereo diagnostic", bond_potential, evidence,
            method="structure_derivation", derivation="Potential double-bond stereochemistry from the candidate graph; cumulated allene bonds require manual interpretation.")
        axial = []
        if "=C=" in smiles:
            axial = [{"type": "allene", "assignment": "not stereogenic in the drawn compound because one terminal allene carbon is CH2"}]
        m["stereochemistry"]["axial"] = f(
            "source drawing and graph", axial, evidence, method="structure_derivation",
            derivation="Manual allene-axis check supplementing limited toolkit axial support.")
        if mixture is not None:
            m["stereochemistry"]["mixture"] = f(mixture, mixture, evidence)
        if stereo_note:
            m["stereochemistry"]["assignment"] = f(stereo_note, stereo_note, evidence)
        unspecified_atoms = [x for x in atom_potential if x["specified"] == "Unspecified"]
        specified_atoms = [x for x in atom_potential if x["specified"] == "Specified"]
        if mixture is not None:
            chirality = "mixture"
        elif unspecified_atoms and specified_atoms:
            chirality = "partially_specified"
        elif unspecified_atoms:
            chirality = "unspecified"
        elif specified_atoms:
            chirality = "specified"
        else:
            chirality = "achiral"
        # Terminal O-allenes here are not E/Z elements; toolkit Bond_Double hints on
        # cumulated bonds are not promoted to source geometry assignments.
        geometry = "not_applicable" if not bond_potential or "=C=" in smiles else "unspecified"
        if unspecified_atoms:
            representation = "partial"
        elif specified_atoms:
            representation = "encoded"
        else:
            representation = "not_applicable"
        m["stereochemistry"]["chirality_status"] = f(
            "source drawing plus standard stereo diagnostic", chirality, evidence,
            method="structure_derivation", derivation="Potential tetrahedral centers were checked; absence of a source assignment remains unspecified rather than racemic.")
        m["stereochemistry"]["geometry_status"] = f(
            "source drawing plus standard stereo diagnostic", geometry, evidence,
            method="structure_derivation", derivation="No source-assigned E/Z geometry was found; terminal alkenes and nonstereogenic cumulated bonds are not assigned E/Z.")
        m["stereochemistry"]["representation_status"] = f(
            "candidate SMILES stereo audit", representation, evidence,
            method="structure_derivation", derivation="Partial is used whenever a potential tetrahedral center is left unassigned.")
        if specified_atoms:
            m["stereochemistry"]["isomeric_smiles"] = f(
                "candidate SMILES stereo audit", diag["canonical_isomeric_smiles"], evidence,
                method="structure_derivation", derivation="Canonical isomeric SMILES retains only source-specified stereochemistry.")
        m["stereochemistry"]["identity_key"] = f(
            "standard stereo diagnostic", diag["stereo_graph_key"], evidence,
            method="structure_derivation", derivation="Hash of the canonical isomeric candidate graph; not a sample identity.")
    else:
        if mixture is not None:
            m["stereochemistry"]["mixture"] = f(mixture, mixture, evidence)
            m["stereochemistry"]["chirality_status"] = f(
                mixture, "mixture", evidence, method="inference",
                derivation="The source explicitly reports an isomer mixture, but the compact structure is not encoded.")
            m["stereochemistry"]["representation_status"] = f(
                "structure not encoded", "not_encoded", evidence, method="inference",
                derivation="The source drawing contains unassigned stereochemistry that was not safely encoded.")
        m["identity_note"] = f(
            "Structure not encoded from the bounded pages", "Structure identity requires expert review", evidence,
            method="inference", derivation="The label/name is retained without inventing a compact structure.")
    molecules[mid] = m
    return mid


def add_generic(key: str, label: str, name: str, evidence: list[str], fixed: str,
                variables: list[tuple[str, str]]) -> str:
    mid = f"{DOC}:mol_{key}"
    if mid in molecules:
        return mid
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    m["names"]["original"] = f(name, name, evidence)
    m["names"]["english"] = f(name, name, evidence)
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": f(fixed, {"description": fixed}, evidence),
        "variables": [
            {
                "label": var,
                "attachment_points": f(desc, desc, evidence),
                "definition": f(desc, desc, evidence),
                "allowed_values": f(),
                "constraints": f(),
            }
            for var, desc in variables
        ],
        "instance_mapping": [],
        "reason": f("Source uses R-group notation", "Generic source drawing cannot be represented as a concrete SMILES", evidence),
    }
    for key2, val in (("chirality_status", "unresolved"), ("geometry_status", "unresolved"),
                      ("representation_status", "not_encoded")):
        m["stereochemistry"][key2] = f(
            "generic R-group drawing", val, evidence, method="inference",
            derivation="Stereogenicity depends on source variables and is not assumed.")
    molecules[mid] = m
    return mid


SMILES = {
    "1a": "CC(C)(C)OC(=O)N(C)C(OCC#C)CCc1ccccc1",
    "2a": "CC(C)(C)OC(=O)N(C)C(OC=C=C)CCc1ccccc1",
    "3a": "C=C(C=O)C(CCc1ccccc1)N(C)C(=O)OC(C)(C)C",
    "1b": "CC(C)(C)OC(=O)N(Cc1ccccc1)C(OCC#C)CCc1ccccc1",
    "2b": "CC(C)(C)OC(=O)N(Cc1ccccc1)C(OC=C=C)CCc1ccccc1",
    "3b": "C=C(C=O)C(CCc1ccccc1)N(Cc1ccccc1)C(=O)OC(C)(C)C",
    "1c": "C=CCN(C(OCC#C)CCc1ccccc1)C(=O)OC(C)(C)C",
    "2c": "C=CCN(C(OC=C=C)CCc1ccccc1)C(=O)OC(C)(C)C",
    "3c": "C=CCN(C(CCc1ccccc1)C(=C)C=O)C(=O)OC(C)(C)C",
    "1d": "CC(C)(C)OC(=O)N(C)C(OCC#C)CCCC",
    "2d": "CC(C)(C)OC(=O)N(C)C(OC=C=C)CCCC",
    "3d": "C=C(C=O)C(CCCC)N(C)C(=O)OC(C)(C)C",
    "1e": "CC(C)(C)OC(=O)N(C)C(OCC#C)C",
    "2e": "CC(C)(C)OC(=O)N(C)C(OC=C=C)C",
    "3e": "C=C(C=O)C(C)N(C)C(=O)OC(C)(C)C",
    "1f": "CC(C)(C)OC(=O)N(C)C(OCC#C)Cc1ccccc1",
    "2f": "CC(C)(C)OC(=O)N(C)C(OC=C=C)Cc1ccccc1",
    "3f": "C=C(C=O)C(Cc1ccccc1)N(C)C(=O)OC(C)(C)C",
    "1g": "CC(C)(C)OC(=O)N(C)C(OCC#C)c1ccccc1",
    "2g": "CC(C)(C)OC(=O)N(C)C(OC=C=C)c1ccccc1",
    "3g": "C=C(C=O)C(c1ccccc1)N(C)C(=O)OC(C)(C)C",
    "1h": "CC(C)(C)OC(=O)N(C)C(OCC#C)C(C)C",
    "2h": "CC(C)(C)OC(=O)N(C)C(OC=C=C)C(C)C",
    "3h": "C=C(C=O)C(C(C)C)N(C)C(=O)OC(C)(C)C",
    "1i": "CC(C)(C)OC(=O)N(C)COCC#C",
    "2i": "CC(C)(C)OC(=O)N(C)COC=C=C",
    "3i": "C=C(C=O)CN(C)C(=O)OC(C)(C)C",
    "1j": "CC(C)(C)OC(=O)N(C)C(OCC#C)C[C@H](C)CCC=C(C)C",
    "2j": "CC(C)(C)OC(=O)N(C)C(OC=C=C)C[C@H](C)CCC=C(C)C",
    "3j": "C=C(C=O)C(N(C)C(=O)OC(C)(C)C)C[C@H](C)CCC=C(C)C",
    "1k": "CC(C)(C)OC(=O)N1CCCC1OCC#C",
    "2k": "CC(C)(C)OC(=O)N1CCCC1OC=C=C",
    "3k": "CC(C)(C)OC(=O)N1CCCC1C(=C)C=O",
    "1l": "CC(C)(C)OC(=O)N1CCCCC1OCC#C",
    "2l": "CC(C)(C)OC(=O)N1CCCCC1OC=C=C",
    "3l": "CC(C)(C)OC(=O)N1CCCCC1C(=C)C=O",
    "4": "C=CC(=C)C(C)N(C)C(=O)OC(C)(C)C",
    "6": "CC(C)(C)OC(=O)N1CCCCC1C(=C)CO",
    "8": "CCOC(=O)C(CC(C=O)CN(C)C(=O)OC(C)(C)C)C(=O)OCC",
    "9a": "CC(C)(C)OC(=O)NC(S(=O)(=O)c1ccccc1)CCc1ccccc1",
    "10a": "CC(C)(C)OC(=O)NC(OCC#C)CCc1ccccc1",
    "9i": "CC(C)(C)OC(=O)NCS(=O)(=O)c1ccccc1",
    "11k": "CC(C)(C)OC(=O)N1CCCC1=O",
    "12k": "CC(C)(C)OC(=O)N1CCCC1O",
}

NAMES = {
    "1a": "tert-Butyl methyl(3-phenyl-1-(prop-2-yn-1-yloxy)propyl)carbamate",
    "2a": "tert-Butyl methyl(3-phenyl-1-(propa-1,2-dien-1-yloxy)propyl)carbamate",
    "3a": "tert-Butyl (2-formyl-5-phenylpent-1-en-3-yl)(methyl)carbamate",
    "1b": "tert-Butyl benzyl(3-phenyl-1-(prop-2-yn-1-yloxy)propyl)carbamate",
    "2b": "tert-Butyl benzyl(3-phenyl-1-(propa-1,2-dien-1-yloxy)propyl)carbamate",
    "3b": "tert-Butyl benzyl(2-formyl-5-phenylpent-1-en-3-yl)carbamate",
    "1c": "tert-Butyl allyl(3-phenyl-1-(prop-2-yn-1-yloxy)propyl)carbamate",
    "2c": "tert-Butyl allyl(3-phenyl-1-(propa-1,2-dien-1-yloxy)propyl)carbamate",
    "3c": "tert-Butyl allyl(2-formyl-5-phenylpent-1-en-3-yl)carbamate",
    "1d": "tert-Butyl methyl(1-(prop-2-yn-1-yloxy)pentyl)carbamate",
    "2d": "tert-Butyl methyl(1-(propa-1,2-dien-1-yloxy)pentyl)carbamate",
    "3d": "tert-Butyl (2-formylhept-1-en-3-yl)(methyl)carbamate",
    "1e": "tert-Butyl methyl(1-(prop-2-yn-1-yloxy)ethyl)carbamate",
    "2e": "tert-Butyl methyl(1-(propa-1,2-dien-1-yloxy)ethyl)carbamate",
    "3e": "tert-Butyl (3-formylbut-3-en-2-yl)(methyl)carbamate",
    "1f": "tert-Butyl methyl(2-phenyl-1-(prop-2-yn-1-yloxy)ethyl)carbamate",
    "2f": "tert-Butyl methyl(2-phenyl-1-(propa-1,2-dien-1-yloxy)ethyl)carbamate",
    "3f": "tert-Butyl (3-formyl-1-phenylbut-3-en-2-yl)(methyl)carbamate",
    "1g": "tert-Butyl methyl(phenyl(prop-2-yn-1-yloxy)methyl)carbamate",
    "2g": "tert-Butyl methyl(phenyl(propa-1,2-dien-1-yloxy)methyl)carbamate",
    "3g": "tert-Butyl (2-formyl-1-phenylallyl)(methyl)carbamate",
    "1h": "tert-Butyl methyl(2-methyl-1-(prop-2-yn-1-yloxy)propyl)carbamate",
    "2h": "tert-Butyl methyl(2-methyl-1-(propa-1,2-dien-1-yloxy)propyl)carbamate",
    "3h": "tert-Butyl (2-formyl-4-methylpent-1-en-3-yl)(methyl)carbamate",
    "1i": "tert-Butyl methyl((prop-2-yn-1-yloxy)methyl)carbamate",
    "2i": "tert-Butyl methyl((propa-1,2-dien-1-yloxy)methyl)carbamate",
    "3i": "tert-Butyl (2-formylallyl)(methyl)carbamate",
    "1j": "tert-Butyl ((3R)-3,7-dimethyl-1-(prop-2-yn-1-yloxy)oct-6-en-1-yl)(methyl)carbamate",
    "2j": "tert-Butyl ((3R)-3,7-dimethyl-1-(propa-1,2-dien-1-yloxy)oct-6-en-1-yl)(methyl)carbamate",
    "3j": "tert-Butyl ((5R)-2-formyl-5,9-dimethyldeca-1,8-dien-3-yl)(methyl)carbamate",
    "1k": "tert-Butyl 2-(prop-2-yn-1-yloxy)pyrrolidine-1-carboxylate",
    "2k": "tert-Butyl 2-(propa-1,2-dien-1-yloxy)pyrrolidine-1-carboxylate",
    "3k": "tert-Butyl 2-(3-oxoprop-1-en-2-yl)pyrrolidine-1-carboxylate",
    "1l": "tert-Butyl 2-(prop-2-yn-1-yloxy)piperidine-1-carboxylate",
    "2l": "tert-Butyl 2-(propa-1,2-dien-1-yloxy)piperidine-1-carboxylate",
    "3l": "tert-Butyl 2-(3-oxoprop-1-en-2-yl)piperidine-1-carboxylate",
    "4": "tert-Butyl methyl(3-methylenepent-4-en-2-yl)carbamate",
    "6": "tert-Butyl 2-(3-hydroxyprop-1-en-2-yl)piperidine-1-carboxylate",
    "8": "Diethyl 2-(3-((tert-butoxycarbonyl)(methyl)amino)-2-formylpropyl)malonate",
    "9a": "tert-Butyl (3-phenyl-1-(phenylsulfonyl)propyl)carbamate",
    "10a": "tert-Butyl (3-phenyl-1-(prop-2-yn-1-yloxy)propyl)carbamate",
    "9i": "tert-Butyl ((phenylsulfonyl)methyl)carbamate",
    "11k": "tert-Butyl 2-oxopyrrolidine-1-carboxylate",
    "12k": "tert-Butyl 2-hydroxypyrrolidine-1-carboxylate",
}


def compound(label: str, pages=(2,)) -> str:
    note = None
    mixture = None
    if label in {"1j", "2j", "3j"}:
        note = "The source name/drawing specifies the remote R center; the N,O-acetal or product amino-bearing center is unassigned."
    if label in {"1j", "2j"}:
        mixture = "mixture of diastereomers and rotamers; proportions not assigned"
    if label == "3j":
        mixture = "5/5 dr; configurations/order unassigned"
    return add_molecule(label, label, NAMES[label], evs(*pages), SMILES[label], stereo_note=note, mixture=mixture)


AGENTS = {
    "triton_b": ("Triton B", "benzyltrimethylammonium hydroxide", "[OH-].C[N+](C)(C)Cc1ccccc1"),
    "dmso": ("DMSO", "dimethyl sulfoxide", "CS(C)=O"),
    "tf2nh": ("Tf2NH", "trifluoromethanesulfonimide", "O=S(=O)(NS(=O)(=O)C(F)(F)F)C(F)(F)F"),
    "ch2cl2": ("CH2Cl2", "dichloromethane", "ClCCl"),
    "bf3_ether": ("BF3·OEt2", "boron trifluoride diethyl etherate", None),
    "sc_otf3": ("Sc(OTf)3", "scandium triflate", None),
    "cu_otf2": ("Cu(OTf)2", "copper(II) triflate", None),
    "ptsoh": ("p-TsOH·H2O", "p-toluenesulfonic acid monohydrate", "Cc1ccc(S(=O)(=O)O)cc1.O"),
    "ppts": ("PPTS", "pyridinium p-toluenesulfonate", "[nH+]1ccccc1.Cc1ccc(S(=O)(=O)[O-])cc1"),
    "tbuk": ("tBuOK", "potassium tert-butoxide", "[K+].[O-]C(C)(C)C"),
    "ylide_salt": ("Ph3PCH3I", "methyltriphenylphosphonium iodide", "[I-].C[P+](c1ccccc1)(c1ccccc1)c1ccccc1"),
    "ether": ("Et2O", "diethyl ether", "CCOCC"),
    "naphthoquinone": ("1,4-naphthoquinone", "1,4-naphthoquinone", "O=C1C=CC(=O)c2ccccc12"),
    "toluene": ("toluene", "toluene", "Cc1ccccc1"),
    "nabh4": ("NaBH4", "sodium borohydride", "[Na+].[BH4-]"),
    "cecl3": ("CeCl3·7H2O", "cerium(III) chloride heptahydrate", None),
    "meoh": ("MeOH", "methanol", "CO"),
    "nah": ("NaH", "sodium hydride", "[Na+].[H-]"),
    "thf": ("THF", "tetrahydrofuran", "C1CCOC1"),
    "diethyl_malonate": ("CH2(CO2Et)2", "diethyl malonate", "CCOC(=O)CC(=O)OCC"),
    "etona": ("EtONa", "sodium ethoxide", "CC[O-].[Na+]"),
    "etoh": ("EtOH", "ethanol", "CCO"),
    "bocnh2": ("tert-butyl carbamate", "tert-butyl carbamate", "CC(C)(C)OC(=O)N"),
    "hydrocinnamaldehyde": ("3-phenylpropionaldehyde", "3-phenylpropionaldehyde", "O=CCCc1ccccc1"),
    "sodium_benzenesulfinate": ("sodium benzenesulfinate dihydrate", "sodium benzenesulfinate dihydrate", None),
    "formic_acid": ("HCO2H", "formic acid", "O=CO"),
    "propargyl_alcohol": ("2-propyn-1-ol", "propargyl alcohol", "C#CCO"),
    "mei": ("MeI", "methyl iodide", "CI"),
    "dmf": ("DMF", "N,N-dimethylformamide", "CN(C)C=O"),
    "benzyl_bromide": ("benzyl bromide", "benzyl bromide", "BrCc1ccccc1"),
    "allyl_bromide": ("allyl bromide", "allyl bromide", "C=CCBr"),
    "pentanal": ("n-valeraldehyde", "pentanal", "CCCCC=O"),
    "acetaldehyde": ("acetaldehyde", "acetaldehyde", "CC=O"),
    "phenylacetaldehyde": ("phenylacetaldehyde", "phenylacetaldehyde", "O=CCc1ccccc1"),
    "benzaldehyde": ("benzaldehyde", "benzaldehyde", "O=Cc1ccccc1"),
    "isobutyraldehyde": ("isobutyraldehyde", "2-methylpropanal", "CC(C)C=O"),
    "formaldehyde": ("HCHO", "formaldehyde", "C=O"),
    "citronellal": ("(+)-citronellal", "(R)-citronellal", "O=CC[C@H](C)CCC=C(C)C"),
    "pyrrolidinone": ("2-pyrrolidinone", "2-pyrrolidinone", "O=C1CCCN1"),
    "nbuli": ("nBuLi", "n-butyllithium", "[Li]CCCC"),
    "boc2o": ("(Boc)2O", "di-tert-butyl dicarbonate", "CC(C)(C)OC(=O)OC(=O)OC(C)(C)C"),
    "redal": ("sodium bis(2-methoxyethoxy)aluminum hydride", "sodium bis(2-methoxyethoxy)aluminum hydride", None),
    "molecular_sieves": ("3 Å molecular sieves", "3 Å molecular sieves", None),
    "nahco3": ("saturated aqueous NaHCO3", "sodium bicarbonate", "[Na+].[O-]C(=O)O"),
    "nh4cl": ("saturated aqueous NH4Cl", "ammonium chloride", "[Cl-].[NH4+]"),
    "brine": ("brine", "aqueous sodium chloride", "[Na+].[Cl-]"),
    "na2so4": ("Na2SO4", "sodium sulfate", "[Na+].[Na+].[O-]S(=O)(=O)[O-]"),
    "water": ("H2O", "water", "O"),
    "etoac": ("EtOAc", "ethyl acetate", "CCOC(C)=O"),
    "hexane": ("n-hexane", "n-hexane", "CCCCCC"),
    "silica": ("silica gel", "silica gel", None),
    "potassium_sodium_tartrate": ("saturated aqueous potassium sodium tartrate", "potassium sodium tartrate solution", None),
}


def agent(key: str, evidence: list[str]) -> str:
    label, name, smiles = AGENTS[key]
    return add_molecule(key, label, name, evidence, smiles, aliases=[label] if label != name else None)


def participant(pid: str, mid: str, role: str, evidence: list[str], stage="reaction", charges=None,
                reference=None):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    for field, (raw, norm) in (charges or {}).items():
        p["charges"][field] = f(raw, norm, evidence)
    if reference:
        p["equivalence_reference"] = f(reference, reference, evidence)
    return p


def step(step_id: str, evidence: list[str], operation: str, pids=(), stage="reaction",
         temperature=None, time=None, atmosphere=None, pressure=None, inherited=False):
    s = st.new_step(step_id, evidence, stage=stage)
    s["participant_ids"] = list(pids)
    method = "inherited_general_condition" if inherited else "original_extraction"
    s["operation"] = f(operation, operation, evidence, method=method)
    if temperature is not None:
        raw, norm = temperature
        s["temperature"] = f(raw, norm, evidence, method=method)
    if time is not None:
        raw, norm = time
        s["time"] = f(raw, norm, evidence, method=method)
    if atmosphere is not None:
        s["atmosphere"] = f(atmosphere, atmosphere, evidence, method=method)
    if pressure is not None:
        s["pressure"] = f(pressure, pressure, evidence, method=method)
    return s


def outcome(oid: str, kind: str, target: list[str], evidence: list[str], raw, norm,
            yield_type=None, basis=None, scope=None):
    o = st.new_outcome(oid, kind, target, evidence)
    o["value"] = f(raw, norm, evidence)
    if yield_type is not None:
        o["yield_type"] = f(yield_type, yield_type, evidence)
    if basis is not None:
        o["measurement_basis"] = f(basis, basis, evidence)
    if scope is not None:
        o["scope"] = f(scope, scope, evidence)
    return o


def classification(labels: list[str], evidence: list[str], multiplicity=None, no_reaction=False):
    c = copy.deepcopy(st.new_reaction(f"{DOC}:tmp", DOC, evidence)["classification"])
    if no_reaction:
        c["labels"] = f("no observed transformation", [], evidence, method="inference",
                         derivation="No product-forming bond is assigned for the failed/no-reaction experiment.")
        c["coverage"] = f("no observed transformation", "not_applicable", evidence, method="inference",
                           derivation="Bond-formation classification is not applicable without an observed product.")
        c["limitations"] = ["No expected-product bonds were invented for a negative result."]
        return c
    changes = []
    counts = multiplicity or {}
    for code in labels:
        count = counts.get(code, 1)
        for _ in range(count):
            order = 2 if code.startswith("DB-") else 3 if code.startswith("TB-") else 1
            pair = code.removeprefix("DB-").removeprefix("TB-").split("-")
            if code == "C-H":
                pair = ["C", "H"]
                kind, before = "hydrogen_gain", 0
            elif code == "O-H":
                pair = ["O", "H"]
                kind, before = "hydrogen_gain", 0
            elif code.startswith("DB-"):
                kind, before = "order_increased", 1
            else:
                kind, before = "formed", 0
            changes.append({
                "change_id": f"b{len(changes)+1}",
                "atoms": [{"atom_map": None, "element": pair[0]}, {"atom_map": None, "element": pair[1]}],
                "before_order": before,
                "after_order": order,
                "change_type": kind,
                **({"hydrogen_count": 1} if kind == "hydrogen_gain" else {}),
            })
    # Preserve directly visible companion losses/order decreases for common
    # transformations in the bounded pages.  These changes do not create
    # taxonomy labels, but prevent the partial comparison from looking like a
    # mechanism or mass-balanced atom mapping.
    if labels == ["DB-C-C", "C-H"]:
        changes.extend([
            {"change_id": f"b{len(changes)+1}", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
             "before_order": 3, "after_order": 2, "change_type": "order_decreased"},
            {"change_id": f"b{len(changes)+2}", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}],
             "before_order": 1, "after_order": 0, "change_type": "hydrogen_loss", "hydrogen_count": 1},
        ])
    elif set(labels) == {"C-H", "O-H"}:
        changes.append(
            {"change_id": f"b{len(changes)+1}", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}],
             "before_order": 2, "after_order": 1, "change_type": "order_decreased"}
        )
    elif set(labels) == {"C-C", "C-H"}:
        changes.append(
            {"change_id": f"b{len(changes)+1}", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
             "before_order": 2, "after_order": 1, "change_type": "order_decreased"}
        )
    label_map = {}
    for change in changes:
        if change["change_type"] not in ("formed", "order_increased", "hydrogen_gain"):
            continue
        lab = st.bond_label(change["atoms"][0]["element"], change["atoms"][1]["element"], change["after_order"])
        label_map.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(change["change_id"])
    c["bond_changes"] = [
        f("structure comparison", change, evidence, method="structure_derivation",
          derivation="Partial net bond-change assignment from the source reactant/product drawings; no complete atom mapping asserted.")
        for change in changes
    ]
    c["labels"] = f(
        "source drawing comparison", list(label_map.values()), evidence, method="structure_derivation",
        derivation="Taxonomy labels cover the explicitly represented formed/promoted bonds in this partial comparison.")
    c["coverage"] = f(
        "manual partial comparison", "partial", evidence, method="inference",
        derivation="No complete atom mapping was supplied; proton transfers, leaving groups and omitted coproducts remain outside coverage.")
    c["limitations"] = ["Partial source-drawing comparison only; no complete atom mapping.",
                         "Implicit-H and omitted coproduct changes may be incomplete."]
    return c


def reaction(rid: str, evidence: list[str], record_type: str, categories: list[str], location: str,
             name: str, transformation: str, participants, steps, outcomes, labels=None,
             solvent=None, scale=None, procedure=None, notes=None, unresolved=None,
             no_reaction=False, general_scheme_ids=None):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, evidence)
    r["location"] = f(location, {"pages": sorted({int(x.split(":page")[1].split(":")[0]) for x in evidence})}, evidence)
    r["reaction_name"] = f(name, name, evidence)
    r["transformation_type"] = f(
        transformation, transformation, evidence, method="inference",
        derivation="Normalized transformation description based on the source title/scheme caption.")
    r["participants"] = participants
    r["steps"] = steps
    r["outcomes"] = outcomes
    r["general_scheme_ids"] = general_scheme_ids or []
    if scale is not None:
        r["reaction_scale"] = f(scale[0], scale[1], evidence)
    if solvent is not None:
        r["solvent_conditions"]["composition"] = f(solvent, solvent, evidence)
    if procedure is not None:
        r["procedure"]["original_text"] = f(procedure, procedure, evidence)
    r["classification"] = classification(labels or [], evidence, no_reaction=no_reaction)
    r["notes"] = [f(x, x, evidence) for x in (notes or [])]
    r["review"] = {
        "status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
        "unresolved": (unresolved or []) + ["Machine-extracted record; chemical identity and structure candidates require expert review."],
    }
    r["completeness"] = "partial"
    reactions.append(r)
    return r["reaction_id"]


# Generic and prior-work arrows on page 1.
g_alkenyl = add_generic("generic_no_alkenyl_acetal", "generic N,O-alkenyl acetal", "generic N,O-alkenyl acetal", ev(1),
                        "N,O-acetal bearing an alkenyl ether; R1/R2/R3/PG variable", [("R1", "N substituent"), ("R2", "acetal-carbon substituent"), ("R3", "alkenyl substituent"), ("PG", "N protecting group")])
g_beta_amino = add_generic("generic_beta_amino_aldehyde", "generic β-amino aldehyde", "generic β-amino aldehyde", ev(1),
                           "β-amino aldehyde product with R1/R2/R3/PG variable", [("R1", "N substituent"), ("R2", "amino-carbon substituent"), ("R3", "alpha substituent"), ("PG", "N protecting group")])
g_1 = add_generic("generic_1", "1", "generic N,O-propargyl acetal 1", ev(1), "N,O-propargyl acetal with R1/R2 variables", [("R1", "N substituent"), ("R2", "acetal-carbon substituent")])
g_2 = add_generic("generic_2", "2", "generic N,O-allenyl acetal 2", evs(1, 2), "N,O-allenyl acetal with R1/R2 variables", [("R1", "N substituent"), ("R2", "acetal-carbon substituent")])
g_3 = add_generic("generic_3", "3", "generic β-amino-α-methylene aldehyde 3", evs(1, 3), "β-amino-α-methylene aldehyde with R1/R2 variables", [("R1", "N substituent"), ("R2", "β substituent")])
acid_generic = add_molecule("acid_generic", "Acid", "unspecified acid catalyst", ev(1))
lit_b_sub = add_molecule("lit_scheme1b_substrate", "Scheme 1b substrate", "prior-work O,O-allenyl acetal", ev(1))
lit_b_prod = add_molecule("lit_scheme1b_product", "Scheme 1b product", "prior-work β-alkoxy-α-methylene aldehyde", ev(1))
lit_c_sub = add_molecule("lit_scheme1c_substrate", "Scheme 1c substrate", "prior-work O,O-propargyl acetal", ev(1))
lit_c_prod = add_molecule("lit_scheme1c_product", "Scheme 1c product", "prior-work O,O-allenyl acetal", ev(1))

reaction("scheme1a_prior_generic_aza_ferrier", ev(1), "literature_example", [], "page 1 Scheme 1a",
         "aza-Ferrier reaction of N,O-alkenyl acetal", "prior-work generic aza-Ferrier rearrangement",
         [participant("substrate", g_alkenyl, "reactant", ev(1)), participant("acid", acid_generic, "catalyst", ev(1)), participant("product", g_beta_amino, "product", ev(1))],
         [step("reaction", ev(1), "acid-catalyzed rearrangement", ["substrate", "acid"])], [], ["C-C"],
         procedure="Scheme 1a: acid-catalyzed cleavage followed by sequential C-C bond formation.",
         unresolved=["This is cited prior work/general chemistry, not a new experiment in this paper."])
reaction("scheme1b_prior_ferrier_46pct", ev(1), "literature_example", [], "page 1 Scheme 1b",
         "Ferrier reaction of O,O-allenyl acetal", "prior-work Lewis-acid Ferrier rearrangement",
         [participant("substrate", lit_b_sub, "reactant", ev(1)), participant("bf3", agent("bf3_ether", ev(1)), "catalyst", ev(1)), participant("product", lit_b_prod, "product", ev(1))],
         [step("reaction", ev(1), "BF3·OEt2-catalyzed rearrangement", ["substrate", "bf3"])],
         [outcome("yield", "yield", [lit_b_prod], ev(1), "46%", q(46, "%", "46%"), basis=None)], ["C-C"],
         procedure="Scheme 1b: BF3·OEt2, 46%.", unresolved=["Prior-work example; yield basis is not reported in the bounded source region."])
reaction("scheme1c_prior_isomerization_84pct", ev(1), "literature_example", [], "page 1 Scheme 1c",
         "isomerization of O,O-propargyl acetal", "prior-work propargyl-to-allenyl isomerization",
         [participant("substrate", lit_c_sub, "reactant", ev(1)), participant("triton", agent("triton_b", ev(1)), "catalyst", ev(1), charges={"mol_percent": ("10 mol%", q(10, "mol%", "10 mol%"))}, reference=lit_c_sub), participant("product", lit_c_prod, "product", ev(1))],
         [step("reaction", ev(1), "Triton B-catalyzed isomerization", ["substrate", "triton"])],
         [outcome("yield", "yield", [lit_c_prod], ev(1), "84%", q(84, "%", "84%"))], ["DB-C-C", "C-H"],
         procedure="Scheme 1c: 10 mol% Triton B, 84%.", unresolved=["Prior-work example; time, temperature and yield basis are not reported in the bounded region."])
generic_scheme_id = reaction("scheme1d_current_generic", ev(1), "generic_scheme", [], "page 1 Scheme 1d",
         "aza-Ferrier reaction of N,O-allenyl acetal", "current-work generic two-stage sequence",
         [participant("propargyl", g_1, "reactant", ev(1)), participant("allene", g_2, "product", ev(1)), participant("aldehyde", g_3, "product", ev(1))],
         [step("isomerization", ev(1), "Triton B isomerization of 1 to 2", ["propargyl", "allene"]), step("aza_ferrier", ev(1), "acid-catalyzed aza-Ferrier rearrangement of 2 to 3", ["allene", "aldehyde"])],
         [], ["DB-C-C", "C-H", "C-C"], procedure="Scheme 1d: 1 → 2 by Triton B; 2 → 3 by acid-catalyzed aza-Ferrier reaction.",
         unresolved=["Generic R-group scheme; individual experiments are recorded separately."])

# Source compounds used throughout pages 2-6.
for letter in "abcdefghijkl":
    compound("1" + letter, pages=(2, 4, 5, 6) if letter <= "k" else (2,))
    compound("2" + letter, pages=(2, 4, 5, 6) if letter <= "k" else (2,))
for letter in "abcdefghijkl":
    compound("3" + letter, pages=(2, 3))
for label in ("4", "6", "8"):
    compound(label, pages=(3,))
mol5 = add_molecule("5", "5", "Diels-Alder adduct 5", ev(3), mixture="5/5 dr; configurations/order unassigned")
mol7 = add_molecule("7", "7", "1,3-oxazine derivative 7", ev(3))

# Unsuccessful tBuOK isomerization preceding Scheme 2.
reaction("control_tbuok_isomerization_1a", ev(2), "specific_experiment", ["control"], "page 2 prose before Scheme 2",
         "traditional base-induced isomerization", "failed/low-yield propargyl-to-allenyl isomerization",
         [participant("p1a", compound("1a"), "reactant", ev(2)), participant("tbuk", agent("tbuk", ev(2)), "reagent", ev(2)), participant("p2a", compound("2a"), "product", ev(2))],
         [step("reaction", ev(2), "isomerization attempt with tBuOK", ["p1a", "tbuk"])],
         [outcome("low_yield", "trace", [compound("2a")], ev(2), "yield of desired 2a was low", "low"),
          outcome("recovery", "recovery", [compound("1a")], ev(2), "recovery of 1a", "reported qualitatively"),
          outcome("decomposition", "decomposition", [], ev(2), "small amounts of decomposition materials", "small amounts")],
         [], no_reaction=True, procedure="The authors tried the traditional base-induced isomerization of 1a with tBuOK; desired 2a was low-yielding with recovery and decomposition.",
         unresolved=["Exact charges, solvent, temperature, time and quantitative yield/recovery are unreported."])

# Scheme 2: 1a-l -> 2a-l. Product masses from pages 4-6 are merged into the same records.
iso_yields = dict(zip("abcdefghijkl", [84, 75, 71, 58, 76, 31, 67, 65, 62, 82, 76, 55]))
iso_masses = {"a": (254, 4), "b": (208, 4), "c": (173, 4), "d": (404, 5), "e": (390, 5),
              "f": (84, 5), "g": (92, 5), "h": (91, 5), "i": (228, 6), "j": (265, 6), "k": (141, 6)}
iso_record_ids = {}
for letter, yld in iso_yields.items():
    pages = [2] + ([iso_masses[letter][1]] if letter in iso_masses else [])
    evidence = evs(*pages)
    m1, m2 = compound("1" + letter), compound("2" + letter)
    parts = [participant("substrate", m1, "reactant", evidence),
             participant("triton", agent("triton_b", evidence), "catalyst", evidence,
                         charges={"mol_percent": ("30 mol%", q(30, "mol%", "30 mol%"))}, reference=m1),
             participant("dmso", agent("dmso", evidence), "solvent", evidence),
             participant("product", m2, "product", evidence)]
    if letter == "a":
        parts[0] = participant("substrate", m1, "reactant", evidence,
                               charges={"mass": ("303 mg", q(303, "mg", "303 mg")), "amount": ("1.00 mmol", q(1.00, "mmol", "1.00 mmol"))})
        parts[1] = participant("triton", agent("triton_b", evidence), "catalyst", evidence,
                               charges={"mass": ("125 mg of 40% MeOH solution", q(125, "mg", "125 mg")), "amount": ("0.30 mmol", q(0.30, "mmol", "0.30 mmol")), "mol_percent": ("30 mol%", q(30, "mol%", "30 mol%"))}, reference=m1)
        parts[2] = participant("dmso", agent("dmso", evidence), "solvent", evidence,
                               charges={"volume": ("0.6 mL + 0.6 mL", q(1.2, "mL", "0.6 mL + 0.6 mL"))})
    outs = [outcome("yield", "yield", [m2], evidence, f"{yld}%", q(yld, "%", f"{yld}%"), "isolated yield", "isolated after chromatography")]
    if letter in iso_masses:
        mass, _ = iso_masses[letter]
        outs.append(outcome("product_mass", "product_amount", [m2], evidence, f"{mass} mg", q(mass, "mg", f"{mass} mg"), scope="isolated product"))
    rid = reaction(f"scheme2_isomerization_1{letter}_to_2{letter}", evidence, "specific_experiment", ["substrate_scope"], f"page 2 Scheme 2; product paragraph page {pages[-1] if len(pages)>1 else 2}",
         "Preparation of N,O-allenyl acetal 2 by isomerization", "propargyl-to-allenyl isomerization",
         parts,
         [step("pretreat_triton", evidence, "remove MeOH from Triton B and dry under reduced pressure", ["triton"], time=("30 min", q(30, "min", "30 min")) if letter == "a" else None, inherited=letter != "a"),
          step("reaction", evidence, "stir substrate with dried Triton B in DMSO", ["substrate", "triton", "dmso"], temperature=("room temperature", "room temperature"), time=("16 h", q(16, "h", "16 h")), atmosphere="Ar", inherited=letter != "a")],
         outs, ["DB-C-C", "C-H"], solvent="DMSO", scale=("1.00 mmol of 1a", {"reference_molecule_id": m1, "amount": q(1.00, "mmol", "1.00 mmol")}) if letter == "a" else None,
         procedure="Prepared by the same procedure with 2a" if letter != "a" else "Triton B MeOH solution was evaporated, dried, dissolved in DMSO, then 1a in DMSO was added and stirred 16 h at room temperature under Ar; aqueous workup and silica chromatography followed.",
         unresolved=["For derivatives other than 2a, exact substrate/catalyst charges and workup volumes are not restated."], general_scheme_ids=[generic_scheme_id])
    iso_record_ids[letter] = rid

# Table 1: nine distinct optimization entries.
table = [
    (1, "bf3_ether", 10, 53, None, None, "approximately 0.3 mmol"),
    (2, "sc_otf3", 10, 32, None, None, "approximately 0.3 mmol"),
    (3, "cu_otf2", 10, 0, 95, None, "approximately 0.3 mmol"),
    (4, "ptsoh", 10, 6, None, 19, "approximately 0.3 mmol"),
    (5, "ppts", 10, 0, 88, None, "approximately 0.3 mmol"),
    (6, "tf2nh", 10, 93, None, None, "approximately 0.3 mmol"),
    (7, "tf2nh", 10, 93, None, None, "1.25 mmol"),
    (8, "tf2nh", 5, 83, None, None, "approximately 0.3 mmol"),
    (9, "tf2nh", 1, 70, None, None, "approximately 0.3 mmol"),
]
aldehyde_side = agent("hydrocinnamaldehyde", ev(2))
table_ids = []
for entry, catalyst_key, loading, yld, recovery, side_yield, scale_raw in table:
    m2a, m3a = compound("2a"), compound("3a")
    cat = agent(catalyst_key, ev(2))
    parts = [participant("substrate", m2a, "reactant", ev(2)),
             participant("catalyst", cat, "catalyst", ev(2), charges={"mol_percent": (f"{loading} mol%", q(loading, "mol%", f"{loading} mol%"))}, reference=m2a),
             participant("solvent", agent("ch2cl2", ev(2)), "solvent", ev(2)),
             participant("product", m3a, "product", ev(2))]
    if side_yield is not None:
        parts.append(participant("side_product", aldehyde_side, "product", ev(2)))
    outs = [outcome("yield_3a", "yield", [m3a], ev(2), f"{yld}%", q(yld, "%", f"{yld}%"), "isolated yield", "isolated")]
    if recovery is not None:
        outs.append(outcome("recovery_2a", "recovery", [m2a], ev(2), f"{recovery}% recovery of 2a", q(recovery, "%", f"{recovery}%"), basis=None))
    if side_yield is not None:
        outs.append(outcome("yield_hydrocinnamaldehyde", "yield", [aldehyde_side], ev(2), f"{side_yield}%", q(side_yield, "%", f"{side_yield}%"), basis=None, scope="side product"))
    table_ids.append(reaction(f"table1_entry_{entry}", ev(2), "specific_experiment", ["optimization"] + (["scale_up"] if entry == 7 else []), f"page 2 Table 1 entry {entry}",
         "acid-catalyzed aza-Ferrier reaction", "aza-Ferrier rearrangement",
         parts, [step("reaction", ev(2), "stir 2a with acid catalyst in CH2Cl2", ["substrate", "catalyst", "solvent"], temperature=("−78 °C", q(-78, "°C", "−78 °C")), time=("1 h", q(1, "h", "1 h")))],
         outs, ["C-C"], solvent="CH2Cl2", scale=(scale_raw + " of 2a", {"reference_molecule_id": m2a, "amount": q(1.25, "mmol", "1.25 mmol") if entry == 7 else {"number": 0.3, "unit": "mmol", "qualifier": "approx", "raw_text": "approximately 0.3 mmol"}}),
         procedure="Table 1: 2a, acid catalyst, CH2Cl2, −78 °C, 1 h.",
         unresolved=(["Recovery measurement basis is not reported."] if recovery is not None else []) + (["Side-product yield measurement basis is not reported."] if side_yield is not None else []),
         general_scheme_ids=[generic_scheme_id]))

# Scheme 3 scope entries, including separate >1 mmol runs for 2e and 2i.
scope_yields = {"b": 91, "c": 86, "d": 91, "e": 89, "f": 75, "g": 69, "h": 89, "i": 81, "j": 2, "k": 94, "l": 93}
scope_ids = []
for letter, yld in scope_yields.items():
    m2, m3 = compound("2" + letter), compound("3" + letter)
    outs = [outcome("yield", "yield", [m3], ev(3), f"{yld}%", q(yld, "%", f"{yld}%"), "isolated yield", "isolated")]
    if letter == "j":
        outs.append(outcome("dr", "selectivity", [m3], ev(3), "5/5 dr", {"ratio": [5, 5], "metric": "dr", "order": "unassigned"}, basis=None, scope="diastereomeric ratio; configurations unassigned"))
        outs.append(outcome("side_products", "decomposition", [], ev(3), "various side products observed by TLC analysis", "various side products"))
    rid = reaction(f"scheme3_scope_2{letter}_to_3{letter}", ev(3), "specific_experiment", ["substrate_scope"], f"page 3 Scheme 3, 2{letter} → 3{letter}",
         "Tf2NH-catalyzed aza-Ferrier reaction", "aza-Ferrier rearrangement",
         [participant("substrate", m2, "reactant", ev(3)), participant("catalyst", agent("tf2nh", ev(3)), "catalyst", ev(3), charges={"mol_percent": ("10 mol%", q(10, "mol%", "10 mol%"))}, reference=m2), participant("solvent", agent("ch2cl2", ev(3)), "solvent", ev(3)), participant("product", m3, "product", ev(3))],
         [step("reaction", ev(3), "stir N,O-allenyl acetal with Tf2NH in CH2Cl2", ["substrate", "catalyst", "solvent"], temperature=("−78 °C", q(-78, "°C", "−78 °C")), time=("1 h", q(1, "h", "1 h")))],
         outs, ["C-C"], solvent="CH2Cl2", scale=("approximately 0.2–0.8 mmol of 2", {"reference_molecule_id": m2, "amount": {"minimum": 0.2, "maximum": 0.8, "unit": "mmol", "raw_text": "approximately 0.2–0.8 mmol"}}),
         procedure="Scheme 3: 10 mol% Tf2NH, CH2Cl2, −78 °C, 1 h; isolated yield shown.",
         unresolved=(["The dr method and individual configurations are not reported."] if letter == "j" else []) + (["Source conflict: Scheme 3 labels the 93% product as 3l, while the adjacent prose says 3j; no value is silently selected for the prose typo."] if letter == "l" else []), general_scheme_ids=[generic_scheme_id])
    scope_ids.append(rid)
    if letter == "l":
        conflict = f(["Scheme 3: 2l → 3l (93%)", "prose: 2k and 2l gave 3k (94%) and 3j (93%)"], None, ev(3), status="conflict")
        conflict["alternatives"] = [
            f("Scheme 3: 2l → 3l (93%)", "2l → 3l, 93%", ev(3)),
            f("prose: 2k and 2l gave 3k (94%) and 3j (93%)", "prose identifies the 93% product as 3j", ev(3)),
        ]
        reactions[-1]["notes"].append(conflict)

for letter, scale_num, yld in [("e", 1.49, 88), ("i", 1.14, 86)]:
    m2, m3 = compound("2" + letter), compound("3" + letter)
    scope_ids.append(reaction(f"scheme3_scaleup_2{letter}_to_3{letter}", ev(3), "specific_experiment", ["substrate_scope", "scale_up"], f"page 3 Scheme 3 parenthetical scale-up for 2{letter}",
         "Tf2NH-catalyzed aza-Ferrier reaction", "aza-Ferrier rearrangement",
         [participant("substrate", m2, "reactant", ev(3), charges={"amount": (f"{scale_num} mmol", q(scale_num, "mmol", f"{scale_num} mmol"))}), participant("catalyst", agent("tf2nh", ev(3)), "catalyst", ev(3), charges={"mol_percent": ("10 mol%", q(10, "mol%", "10 mol%"))}, reference=m2), participant("solvent", agent("ch2cl2", ev(3)), "solvent", ev(3)), participant("product", m3, "product", ev(3))],
         [step("reaction", ev(3), "stir N,O-allenyl acetal with Tf2NH in CH2Cl2", ["substrate", "catalyst", "solvent"], temperature=("−78 °C", q(-78, "°C", "−78 °C")), time=("1 h", q(1, "h", "1 h")))],
         [outcome("yield", "yield", [m3], ev(3), f"{yld}%", q(yld, "%", f"{yld}%"), "isolated yield", "isolated")], ["C-C"], solvent="CH2Cl2", scale=(f"{scale_num} mmol of 2{letter}", {"reference_molecule_id": m2, "amount": q(scale_num, "mmol", f"{scale_num} mmol")}),
         procedure="Scheme 3 footnote: yield in parentheses obtained from an over-1 mmol scale experiment.", general_scheme_ids=[generic_scheme_id]))

# Scheme 4 derivatizations.
scheme4_ids = []
m3e, m4 = compound("3e"), compound("4")
scheme4_ids.append(reaction("scheme4a_wittig_3e_to_4", ev(3), "specific_experiment", ["derivatization"], "page 3 Scheme 4a first arrow",
    "Wittig olefination of 3e", "carbonyl methylenation",
    [participant("substrate", m3e, "reactant", ev(3)), participant("phosphonium", agent("ylide_salt", ev(3)), "reagent", ev(3)), participant("base", agent("tbuk", ev(3)), "reagent", ev(3)), participant("solvent", agent("ether", ev(3)), "solvent", ev(3)), participant("product", m4, "product", ev(3))],
    [step("addition_at_0c", ev(3), "perform Wittig olefination at 0 °C", ["substrate", "phosphonium", "base", "solvent"], temperature=("0 °C", q(0, "°C", "0 °C"))), step("warm_to_rt", ev(3), "allow mixture to warm from 0 °C to room temperature", ["substrate", "phosphonium", "base", "solvent"], temperature=("0 °C to rt", "warm to room temperature"))],
    [outcome("yield", "yield", [m4], ev(3), "78%", q(78, "%", "78%"), basis=None)], ["DB-C-C"], solvent="Et2O", procedure="Scheme 4a: Ph3PCH3I, tBuOK, Et2O, 0 °C to rt; 4, 78%.",
    unresolved=["Reaction/addition times and yield measurement basis are not stated on page 3."]))
scheme4_ids.append(reaction("scheme4a_diels_alder_4_to_5", ev(3), "specific_experiment", ["derivatization"], "page 3 Scheme 4a second arrow",
    "Diels-Alder reaction of 4 with 1,4-naphthoquinone", "Diels-Alder cycloaddition",
    [participant("diene", m4, "reactant", ev(3)), participant("dienophile", agent("naphthoquinone", ev(3)), "reactant", ev(3)), participant("solvent", agent("toluene", ev(3)), "solvent", ev(3)), participant("product", mol5, "product", ev(3))],
    [step("reflux", ev(3), "heat in toluene at reflux", ["diene", "dienophile", "solvent"], temperature=("reflux", "reflux"))],
    [outcome("yield", "yield", [mol5], ev(3), "85%", q(85, "%", "85%"), basis=None), outcome("dr", "selectivity", [mol5], ev(3), "5/5 dr", {"ratio": [5, 5], "metric": "dr", "order": "unassigned"}, basis=None, scope="diastereomeric ratio; configurations unassigned")],
    ["C-C"], solvent="toluene", procedure="Scheme 4a: 1,4-naphthoquinone, toluene, reflux; 5, 85%, 5/5 dr.", unresolved=["Reaction time, yield basis, dr method and configurations are not reported on page 3."]))
m3l, m6 = compound("3l"), compound("6")
scheme4_ids.append(reaction("scheme4b_reduction_3l_to_6", ev(3), "specific_experiment", ["derivatization"], "page 3 Scheme 4b first arrow",
    "1,2-reduction of 3l", "aldehyde reduction",
    [participant("substrate", m3l, "reactant", ev(3)), participant("nabh4", agent("nabh4", ev(3)), "reagent", ev(3)), participant("cecl3", agent("cecl3", ev(3)), "reagent", ev(3)), participant("solvent", agent("meoh", ev(3)), "solvent", ev(3)), participant("product", m6, "product", ev(3))],
    [step("reaction", ev(3), "reduce 3l with NaBH4/CeCl3·7H2O in MeOH", ["substrate", "nabh4", "cecl3", "solvent"], temperature=("0 °C", q(0, "°C", "0 °C")))],
    [outcome("yield", "yield", [m6], ev(3), "88%", q(88, "%", "88%"), basis=None)], ["C-H", "O-H"], solvent="MeOH", procedure="Scheme 4b: NaBH4, CeCl3·7H2O, MeOH, 0 °C; 6, 88%.", unresolved=["Reaction time and yield basis are not reported on page 3."]))
scheme4_ids.append(reaction("scheme4b_cyclization_6_to_7", ev(3), "specific_experiment", ["derivatization"], "page 3 Scheme 4b second arrow",
    "base-promoted cyclization of 6", "intramolecular cyclization",
    [participant("substrate", m6, "reactant", ev(3)), participant("nah", agent("nah", ev(3)), "reagent", ev(3)), participant("solvent", agent("thf", ev(3)), "solvent", ev(3)), participant("product", mol7, "product", ev(3))],
    [step("at_0c", ev(3), "treat 6 with NaH in THF at 0 °C", ["substrate", "nah", "solvent"], temperature=("0 °C", q(0, "°C", "0 °C"))), step("warm_to_rt", ev(3), "warm from 0 °C to room temperature", ["substrate", "nah", "solvent"], temperature=("0 °C to rt", "warm to room temperature"))],
    [outcome("yield", "yield", [mol7], ev(3), "74%", q(74, "%", "74%"), basis=None)], ["C-O"], solvent="THF", procedure="Scheme 4b: NaH, THF, 0 °C to rt; 7, 74%.", unresolved=["Stage times and yield basis are not reported on page 3; product 7 compact structure is not encoded."]))
m3i, m8 = compound("3i"), compound("8")
scheme4_ids.append(reaction("scheme4c_michael_3i_to_8", ev(3), "specific_experiment", ["derivatization"], "page 3 Scheme 4c",
    "Michael addition of diethyl malonate to 3i", "conjugate addition",
    [participant("substrate", m3i, "reactant", ev(3)), participant("malonate", agent("diethyl_malonate", ev(3)), "reactant", ev(3)), participant("catalyst", agent("etona", ev(3)), "catalyst", ev(3), charges={"mol_percent": ("10 mol%", q(10, "mol%", "10 mol%"))}, reference=m3i), participant("solvent", agent("etoh", ev(3)), "solvent", ev(3)), participant("product", m8, "product", ev(3))],
    [step("reaction", ev(3), "stir 3i and diethyl malonate with 10 mol% EtONa in EtOH", ["substrate", "malonate", "catalyst", "solvent"], temperature=("rt", "room temperature"))],
    [outcome("yield", "yield", [m8], ev(3), "79%", q(79, "%", "79%"), basis=None)], ["C-C", "C-H"], solvent="EtOH", procedure="Scheme 4c: 10 mol% EtONa, diethyl malonate, EtOH, rt; 8, 79%.", unresolved=["Reaction time and yield basis are not reported on page 3."]))

# Failed analogous cyclization mentioned in prose.
reaction("control_acyclic_3a_cyclization_failure", ev(3), "specific_experiment", ["control", "derivatization"], "page 3 prose below Scheme 3",
         "attempted analogous transformation using acyclic 3a", "failed cyclization control",
         [participant("substrate", compound("3a"), "reactant", ev(3))],
         [step("attempt", ev(3), "attempt the same reduction/cyclization transformation using acyclic 3a", ["substrate"])],
         [outcome("no_product", "no_reaction", [], ev(3), "the cyclization resulted in failure", "failed cyclization")],
         [], no_reaction=True, procedure="The authors attempted the same transformation using acyclic 3a; the cyclization resulted in failure.",
         unresolved=["The exact intermediate entering cyclization, reagent charges, conditions and analytical basis are not restated."])

# Detailed precursor preparations on pages 3-6.
for label, pages in (("9a", (3, 4)), ("10a", (4,)), ("9i", (6,)), ("11k", (6,)), ("12k", (6,))):
    add_molecule(label, label, NAMES[label], evs(*pages), SMILES[label])

# 9a preparation.
m9a = add_molecule("9a", "9a", NAMES["9a"], evs(3, 4), SMILES["9a"])
reaction("prep_9a", evs(3, 4), "specific_experiment", ["other"], "pages 3-4, Representative Procedure for Preparation of 1a, Step 1",
    "Preparation of 9a", "α-amidosulfone formation",
    [participant("bocnh2", agent("bocnh2", evs(3, 4)), "reactant", evs(3, 4), charges={"mass": ("469 mg", q(469, "mg", "469 mg")), "amount": ("4.00 mmol", q(4.00, "mmol", "4.00 mmol"))}),
     participant("aldehyde", agent("hydrocinnamaldehyde", evs(3, 4)), "reactant", evs(3, 4), charges={"volume": ("0.58 mL", q(0.58, "mL", "0.58 mL")), "amount": ("4.4 mmol", q(4.4, "mmol", "4.4 mmol"))}),
     participant("sulfinate", agent("sodium_benzenesulfinate", evs(3, 4)), "reagent", evs(3, 4), charges={"mass": ("0.80 g", q(0.80, "g", "0.80 g")), "amount": ("6.0 mmol", q(6.0, "mmol", "6.0 mmol"))}),
     participant("formic", agent("formic_acid", evs(3, 4)), "reagent", evs(3, 4), charges={"volume": ("1.0 mL", q(1.0, "mL", "1.0 mL")), "amount": ("27 mmol", q(27, "mmol", "27 mmol"))}),
     participant("thf", agent("thf", evs(3, 4)), "solvent", evs(3, 4), charges={"volume": ("1.6 mL", q(1.6, "mL", "1.6 mL"))}),
     participant("water", agent("water", evs(3, 4)), "solvent", evs(3, 4), charges={"volume": ("4.0 mL", q(4.0, "mL", "4.0 mL"))}),
     participant("product", m9a, "product", evs(3, 4))],
    [step("reaction", evs(3, 4), "stir mixture in THF/H2O", ["bocnh2", "aldehyde", "sulfinate", "formic", "thf", "water"], temperature=("room temperature", "room temperature"), time=("18 h", q(18, "h", "18 h"))),
     step("quench", ev(4), "quench with saturated aqueous NaHCO3 at 0 °C", [], stage="workup", temperature=("0 °C", q(0, "°C", "0 °C"))),
     step("isolate", ev(4), "extract with CH2Cl2; wash with brine; dry over Na2SO4; evaporate", [], stage="workup")],
    [outcome("mass", "product_amount", [m9a], ev(4), "1.48 g", q(1.48, "g", "1.48 g"), scope="crude white solid"), outcome("yield", "yield", [m9a], ev(4), "99% yield", q(99, "%", "99%"), "yield", "isolated crude product used without purification")],
    ["C-N", "C-S"], solvent="THF/H2O", scale=("4.00 mmol tert-butyl carbamate", {"reference_molecule_id": agent("bocnh2", evs(3,4)), "amount": q(4.00, "mmol", "4.00 mmol")}),
    procedure="A mixture of tert-butyl carbamate, 3-phenylpropionaldehyde, sodium benzenesulfinate dihydrate and HCO2H in THF/H2O was stirred 18 h at room temperature, quenched at 0 °C and worked up to give 9a.")

# 10a preparation.
m10a = add_molecule("10a", "10a", NAMES["10a"], ev(4), SMILES["10a"])
reaction("prep_10a", ev(4), "specific_experiment", ["other"], "page 4, Representative Procedure for Preparation of 1a, Step 2",
    "Preparation of 10a", "substitution of α-amidosulfone by propargyl alcohol",
    [participant("propargyl_alcohol", agent("propargyl_alcohol", ev(4)), "reactant", ev(4), charges={"volume": ("201 μL", q(201, "uL", "201 μL")), "amount": ("3.41 mmol", q(3.41, "mmol", "3.41 mmol"))}),
     participant("nah", agent("nah", ev(4)), "reagent", ev(4), charges={"mass": ("0.16 g, 60 wt % in oil", q(0.16, "g", "0.16 g")), "amount": ("4.0 mmol", q(4.0, "mmol", "4.0 mmol"))}),
     participant("thf", agent("thf", ev(4)), "solvent", ev(4), charges={"volume": ("34 mL", q(34, "mL", "34 mL"))}),
     participant("9a", m9a, "reactant", ev(4), charges={"mass": ("1.40 g", q(1.40, "g", "1.40 g")), "amount": ("3.73 mmol", q(3.73, "mmol", "3.73 mmol"))}),
     participant("product", m10a, "product", ev(4))],
    [step("alkoxide", ev(4), "treat propargyl alcohol in THF with NaH at 0 °C", ["propargyl_alcohol", "nah", "thf"], temperature=("0 °C", q(0, "°C", "0 °C"))),
     step("pre_stir", ev(4), "stir after NaH addition", ["propargyl_alcohol", "nah", "thf"], temperature=("room temperature", "room temperature"), time=("30 min", q(30, "min", "30 min")), atmosphere="Ar"),
     step("add_9a", ev(4), "cool to 0 °C and add solid 9a", ["9a"], temperature=("0 °C", q(0, "°C", "0 °C"))),
     step("reaction", ev(4), "stir after 9a addition", ["9a", "propargyl_alcohol", "nah"], temperature=("room temperature", "room temperature"), time=("1 h", q(1, "h", "1 h"))),
     step("workup", ev(4), "quench with saturated NH4Cl at 0 °C; extract EtOAc, wash brine, dry Na2SO4, evaporate", [], stage="workup", temperature=("0 °C", q(0, "°C", "0 °C"))),
     step("purify", ev(4), "silica chromatography, n-hexane/EtOAc 6:1 to 5:1", [], stage="purification")],
    [outcome("mass", "product_amount", [m10a], ev(4), "0.92 g", q(0.92, "g", "0.92 g"), scope="isolated product"), outcome("yield", "yield", [m10a], ev(4), "93% yield", q(93, "%", "93%"), "isolated yield", "silica chromatography")],
    ["C-O"], solvent="THF", scale=("3.41 mmol propargyl alcohol", {"reference_molecule_id": agent("propargyl_alcohol", ev(4)), "amount": q(3.41, "mmol", "3.41 mmol")}),
    procedure="Propargyl alcohol was deprotonated with NaH; 9a was added at 0 °C and the mixture was stirred 1 h at room temperature under Ar, worked up, and chromatographed.")

# 1a methylation.
m1a = compound("1a")
reaction("prep_1a_methylation", ev(4), "specific_experiment", ["other"], "page 4, Representative Procedure for Preparation of 1a, Step 3",
    "Preparation of 1a", "N-methylation",
    [participant("10a", m10a, "reactant", ev(4), charges={"mass": ("0.92 g", q(0.92, "g", "0.92 g")), "amount": ("3.2 mmol", q(3.2, "mmol", "3.2 mmol"))}),
     participant("mei", agent("mei", ev(4)), "reagent", ev(4), charges={"volume": ("0.30 mL", q(0.30, "mL", "0.30 mL")), "amount": ("4.8 mmol", q(4.8, "mmol", "4.8 mmol"))}),
     participant("nah", agent("nah", ev(4)), "reagent", ev(4), charges={"mass": ("0.19 g, 60 wt % in oil", q(0.19, "g", "0.19 g")), "amount": ("4.8 mmol", q(4.8, "mmol", "4.8 mmol"))}),
     participant("dmf", agent("dmf", ev(4)), "solvent", ev(4), charges={"volume": ("11 mL", q(11, "mL", "11 mL"))}), participant("product", m1a, "product", ev(4))],
    [step("addition", ev(4), "treat solution of 10a and MeI in DMF with NaH at 0 °C", ["10a", "mei", "nah", "dmf"], temperature=("0 °C", q(0, "°C", "0 °C"))), step("reaction", ev(4), "stir under Ar", ["10a", "mei", "nah", "dmf"], temperature=("room temperature", "room temperature"), time=("1 h", q(1, "h", "1 h")), atmosphere="Ar"), step("workup", ev(4), "cool and quench with saturated NH4Cl at 0 °C; extract EtOAc, wash brine, dry Na2SO4, evaporate", [], stage="workup", temperature=("0 °C", q(0, "°C", "0 °C"))), step("purify", ev(4), "silica chromatography, n-hexane/EtOAc 6:1 to 4:1", [], stage="purification")],
    [outcome("mass", "product_amount", [m1a], ev(4), "840 mg", q(840, "mg", "840 mg"), scope="isolated product"), outcome("yield", "yield", [m1a], ev(4), "87% yield", q(87, "%", "87%"), "isolated yield", "silica chromatography")], ["C-N"], solvent="DMF", scale=("3.2 mmol 10a", {"reference_molecule_id": m10a, "amount": q(3.2, "mmol", "3.2 mmol")}), procedure="10a and MeI in DMF were treated with NaH at 0 °C, stirred 1 h at room temperature under Ar, worked up, and chromatographed.")

# 1b and 1c: Step-3 variants explicitly described on page 4.
for letter, reagent_key, mass, yld, ratio in [("b", "benzyl_bromide", 290, 75, "15:1"), ("c", "allyl_bromide", 248, 62, "15:1")]:
    product_mid = compound("1" + letter)
    reaction(f"prep_1{letter}_alkylation", ev(4), "specific_experiment", ["other"], f"page 4, product 1{letter} paragraph",
        f"Preparation of 1{letter}", "N-alkylation",
        [participant("10a", m10a, "reactant", ev(4)), participant("alkyl_halide", agent(reagent_key, ev(4)), "reagent", ev(4)), participant("nah", agent("nah", ev(4)), "reagent", ev(4)), participant("dmf", agent("dmf", ev(4)), "solvent", ev(4)), participant("product", product_mid, "product", ev(4))],
        [step("addition", ev(4), f"perform Step 3 as for 1a using {AGENTS[reagent_key][0]} instead of MeI", ["10a", "alkyl_halide", "nah", "dmf"], temperature=("0 °C", q(0, "°C", "0 °C")), inherited=True), step("reaction", ev(4), "stir at room temperature under Ar", ["10a", "alkyl_halide", "nah", "dmf"], temperature=("room temperature", "room temperature"), time=("1 h", q(1, "h", "1 h")), atmosphere="Ar", inherited=True), step("purify", ev(4), f"silica chromatography, n-hexane/EtOAc {ratio}", [], stage="purification")],
        [outcome("mass", "product_amount", [product_mid], ev(4), f"{mass} mg", q(mass, "mg", f"{mass} mg"), scope="isolated product"), outcome("yield", "yield", [product_mid], ev(4), f"{yld}% yield from 10a", q(yld, "%", f"{yld}%"), "isolated yield from 10a", "silica chromatography")], ["C-N"], solvent="DMF", procedure=f"Prepared by the same procedure with 1a using {AGENTS[reagent_key][0]} instead of MeI in Step 3.", unresolved=["Exact charges of 10a, alkyl halide, NaH and DMF are not restated for this run."])

# Overall three-step preparations of 1d-h and 1j; exact charges are not restated.
overall = {
    "d": ("pentanal", 175, 54, (4, 5)), "e": ("acetaldehyde", 711, 62, (5,)), "f": ("phenylacetaldehyde", 141, 22, (5,)),
    "g": ("benzaldehyde", 136, 58, (5,)), "h": ("isobutyraldehyde", 198, 78, (5,)), "j": ("citronellal", 877, 72, (6,)),
}
overall_ids = []
for letter, (aldehyde_key, mass, yld, pages) in overall.items():
    evidence = evs(*pages)
    product_mid = compound("1" + letter)
    overall_ids.append(reaction(f"prep_1{letter}_overall_from_aldehyde", evidence, "specific_experiment", ["other"], f"page {'-'.join(map(str, pages))}, product 1{letter} paragraph",
        f"Preparation of 1{letter}", "three-step N,O-propargyl acetal synthesis",
        [participant("bocnh2", agent("bocnh2", evidence), "reactant", evidence), participant("aldehyde", agent(aldehyde_key, evidence), "reactant", evidence), participant("sulfinate", agent("sodium_benzenesulfinate", evidence), "reagent", evidence), participant("formic", agent("formic_acid", evidence), "reagent", evidence), participant("propargyl_alcohol", agent("propargyl_alcohol", evidence), "reactant", evidence), participant("nah", agent("nah", evidence), "reagent", evidence), participant("mei", agent("mei", evidence), "reagent", evidence), participant("product", product_mid, "product", evidence)],
        [step("step1", evidence, f"perform Step 1 as for 1a using {AGENTS[aldehyde_key][0]} instead of 3-phenylpropionaldehyde", ["bocnh2", "aldehyde", "sulfinate", "formic"], temperature=("room temperature", "room temperature"), time=("18 h", q(18, "h", "18 h")), inherited=True), step("step2", evidence, "perform propargyl substitution as in the 1a procedure", ["propargyl_alcohol", "nah"], inherited=True), step("step3", evidence, "perform N-methylation as in the 1a procedure", ["mei", "nah"], inherited=True)],
        [outcome("mass", "product_amount", [product_mid], evidence, f"{mass} mg", q(mass, "mg", f"{mass} mg"), scope="isolated product"), outcome("overall_yield", "yield", [product_mid], evidence, f"{yld}% overall yield", q(yld, "%", f"{yld}%"), "overall isolated yield", "silica chromatography", scope=f"overall from {AGENTS[aldehyde_key][0]} to 1{letter}")],
        ["C-N", "C-O"], procedure=f"Prepared by the same procedure with 1a using {AGENTS[aldehyde_key][0]} instead of 3-phenylpropionaldehyde in Step 1.",
        unresolved=["This is an overall multistep result; intermediate yields and run-specific charges/volumes are not restated and the overall yield is not assigned to only the last step."]))

# 9i and subsequent overall conversion to 1i.
m9i = add_molecule("9i", "9i", NAMES["9i"], ev(6), SMILES["9i"])
reaction("prep_9i", ev(6), "specific_experiment", ["other"], "page 6, product 1i paragraph, Step 1",
    "Preparation of 9i", "aminosulfone formation from formaldehyde",
    [participant("bocnh2", agent("bocnh2", ev(6)), "reactant", ev(6), charges={"mass": ("0.75 g", q(0.75, "g", "0.75 g")), "amount": ("6.4 mmol", q(6.4, "mmol", "6.4 mmol"))}), participant("hcho", agent("formaldehyde", ev(6)), "reactant", ev(6), charges={"volume": ("1.0 mL of 37 wt% aqueous solution", q(1.0, "mL", "1.0 mL")), "amount": ("ca. 13 mmol", {"number": 13, "unit": "mmol", "qualifier": "approx", "raw_text": "ca. 13 mmol"})}), participant("sulfinate", agent("sodium_benzenesulfinate", ev(6)), "reagent", ev(6), charges={"mass": ("2.56 g", q(2.56, "g", "2.56 g")), "amount": ("12.8 mmol", q(12.8, "mmol", "12.8 mmol"))}), participant("formic", agent("formic_acid", ev(6)), "reagent", ev(6), charges={"volume": ("0.75 mL", q(0.75, "mL", "0.75 mL")), "amount": ("20 mmol", q(20, "mmol", "20 mmol"))}), participant("meoh", agent("meoh", ev(6)), "solvent", ev(6), charges={"volume": ("8 mL", q(8, "mL", "8 mL"))}), participant("water", agent("water", ev(6)), "solvent", ev(6), charges={"volume": ("16 mL", q(16, "mL", "16 mL"))}), participant("product", m9i, "product", ev(6))],
    [step("reaction", ev(6), "stir in MeOH/H2O", ["bocnh2", "hcho", "sulfinate", "formic", "meoh", "water"], temperature=("room temperature", "room temperature"), time=("5 days", q(5, "d", "5 days"))), step("isolate", ev(6), "filter crystals, wash with Et2O, dissolve in CH2Cl2, dry over Na2SO4 and evaporate", [], stage="workup")],
    [outcome("mass", "product_amount", [m9i], ev(6), "826 mg", q(826, "mg", "826 mg"), scope="crude white solid"), outcome("yield", "yield", [m9i], ev(6), "48% yield", q(48, "%", "48%"), "yield", "isolated by filtration; crude used without purification")], ["C-N", "C-S"], solvent="MeOH/H2O", scale=("6.4 mmol tert-butyl carbamate", {"reference_molecule_id": agent("bocnh2", ev(6)), "amount": q(6.4, "mmol", "6.4 mmol")}), procedure="The mixture was stirred 5 days at room temperature; crystals were filtered, washed, dissolved, dried and concentrated to give 9i.")

m1i = compound("1i")
reaction("prep_1i_from_9i_overall", ev(6), "specific_experiment", ["other"], "page 6, product 1i paragraph after Step 1",
    "Preparation of 1i from 9i", "two-step propargyl substitution and N-methylation",
    [participant("9i", m9i, "reactant", ev(6)), participant("propargyl_alcohol", agent("propargyl_alcohol", ev(6)), "reactant", ev(6)), participant("nah", agent("nah", ev(6)), "reagent", ev(6)), participant("mei", agent("mei", ev(6)), "reagent", ev(6)), participant("product", m1i, "product", ev(6))],
    [step("substitution", ev(6), "follow the procedure for 1a to replace phenylsulfonyl by propargyloxy", ["9i", "propargyl_alcohol", "nah"], inherited=True), step("methylation", ev(6), "follow the procedure for 1a for N-methylation", ["mei", "nah"], inherited=True)],
    [outcome("mass", "product_amount", [m1i], ev(6), "368 mg", q(368, "mg", "368 mg"), scope="isolated product"), outcome("overall_yield", "yield", [m1i], ev(6), "67% overall yield from 9i", q(67, "%", "67%"), "overall isolated yield", "silica chromatography", scope="9i to 1i over the following steps")], ["C-O", "C-N"], procedure="The following procedure was the same as that with 1a; 368 mg, 67% overall yield from 9i.", unresolved=["Run-specific charges and individual step yields are not restated; 67% remains an overall yield from 9i."])

# 1k three-step sequence.
m11k = add_molecule("11k", "11k", NAMES["11k"], ev(6), SMILES["11k"])
m12k = add_molecule("12k", "12k", NAMES["12k"], ev(6), SMILES["12k"])
m1k = compound("1k")
reaction("prep_11k_boc_protection", ev(6), "specific_experiment", ["other"], "page 6, Representative Procedure for Preparation of 1k, Step 1",
    "Preparation of 11k", "N-Boc protection of 2-pyrrolidinone",
    [participant("lactam", agent("pyrrolidinone", ev(6)), "reactant", ev(6), charges={"volume": ("0.38 mL", q(0.38, "mL", "0.38 mL")), "amount": ("5.0 mmol", q(5.0, "mmol", "5.0 mmol"))}), participant("nbuli", agent("nbuli", ev(6)), "reagent", ev(6), charges={"volume": ("3.1 mL of 1.6 M hexane solution", q(3.1, "mL", "3.1 mL")), "amount": ("5.0 mmol", q(5.0, "mmol", "5.0 mmol"))}), participant("boc2o", agent("boc2o", ev(6)), "reagent", ev(6), charges={"volume": ("1.26 mL", q(1.26, "mL", "1.26 mL")), "amount": ("5.5 mmol", q(5.5, "mmol", "5.5 mmol"))}), participant("thf", agent("thf", ev(6)), "solvent", ev(6), charges={"volume": ("12.5 mL", q(12.5, "mL", "12.5 mL"))}), participant("product", m11k, "product", ev(6))],
    [step("deprotonation", ev(6), "treat 2-pyrrolidinone in THF with nBuLi at 0 °C", ["lactam", "nbuli", "thf"], temperature=("0 °C", q(0, "°C", "0 °C")), time=("30 min", q(30, "min", "30 min")), atmosphere="Ar"), step("boc_addition", ev(6), "add (Boc)2O at 0 °C", ["boc2o"], temperature=("0 °C", q(0, "°C", "0 °C"))), step("hold_rt", ev(6), "stir after (Boc)2O addition", ["lactam", "nbuli", "boc2o"], temperature=("room temperature", "room temperature"), time=("2 h", q(2, "h", "2 h"))), step("workup", ev(6), "quench with saturated NH4Cl; extract EtOAc, wash brine, dry Na2SO4 and evaporate", [], stage="workup")],
    [outcome("crude_product", "other", [m11k], ev(6), "crude product used without purification", "crude product")], ["C-N"], solvent="THF", scale=("5.0 mmol 2-pyrrolidinone", {"reference_molecule_id": agent("pyrrolidinone", ev(6)), "amount": q(5.0, "mmol", "5.0 mmol")}), procedure="2-Pyrrolidinone was treated with nBuLi at 0 °C, then (Boc)2O, stirred 2 h at room temperature and worked up to crude 11k.", unresolved=["No isolated mass or step yield is reported for crude 11k."])

reaction("prep_12k_reduction", ev(6), "specific_experiment", ["other"], "page 6, Representative Procedure for Preparation of 1k, Step 2",
    "Preparation of 12k", "partial reduction of lactam",
    [participant("11k", m11k, "reactant", ev(6)), participant("redal", agent("redal", ev(6)), "reagent", ev(6), charges={"volume": ("6.0 mL of 1.0 M toluene solution", q(6.0, "mL", "6.0 mL")), "amount": ("6.0 mmol", q(6.0, "mmol", "6.0 mmol"))}), participant("ch2cl2", agent("ch2cl2", ev(6)), "solvent", ev(6), charges={"volume": ("25 mL", q(25, "mL", "25 mL"))}), participant("product", m12k, "product", ev(6))],
    [step("reaction", ev(6), "add Red-Al solution to 11k in CH2Cl2 at −78 °C", ["11k", "redal", "ch2cl2"], temperature=("−78 °C", q(-78, "°C", "−78 °C")), time=("2 h", q(2, "h", "2 h")), atmosphere="Ar"), step("methanol_quench", ev(6), "quench with MeOH", [], stage="workup", temperature=("−78 °C", q(-78, "°C", "−78 °C"))), step("tartrate", ev(6), "add saturated potassium sodium tartrate and stir until a clear two-phase solution", [], stage="workup", temperature=("room temperature", "room temperature"), time=("a few hours", "a few hours")), step("isolate", ev(6), "dilute with H2O, extract CH2Cl2, dry Na2SO4 and evaporate", [], stage="workup"), step("purify", ev(6), "silica chromatography, n-hexane/EtOAc 3:1 to 1.5:1", [], stage="purification")],
    [outcome("mass", "product_amount", [m12k], ev(6), "633 mg", q(633, "mg", "633 mg"), scope="isolated product"), outcome("yield", "yield", [m12k], ev(6), "68% yield from Step 1", q(68, "%", "68%"), "yield from Step 1", "silica chromatography", scope="combined from Step 1 starting material to 12k")], ["C-H", "O-H"], solvent="CH2Cl2", procedure="Red-Al was added to 11k in CH2Cl2 at −78 °C under Ar and stirred 2 h; methanol/tartrate workup and chromatography gave 12k.", unresolved=["The source reports 68% yield from Step 1; it is retained as a multistep-span yield rather than a last-step-only yield."])

reaction("prep_1k_acetalization", ev(6), "specific_experiment", ["other"], "page 6, Representative Procedure for Preparation of 1k, Step 3",
    "Preparation of 1k", "formation of N,O-propargyl acetal",
    [participant("12k", m12k, "reactant", ev(6), charges={"mass": ("633 mg", q(633, "mg", "633 mg")), "amount": ("3.38 mmol", q(3.38, "mmol", "3.38 mmol"))}), participant("propargyl", agent("propargyl_alcohol", ev(6)), "reactant", ev(6), charges={"volume": ("0.24 mL", q(0.24, "mL", "0.24 mL")), "amount": ("4.1 mmol", q(4.1, "mmol", "4.1 mmol"))}), participant("sc", agent("sc_otf3", ev(6)), "catalyst", ev(6), charges={"mass": ("33 mg", q(33, "mg", "33 mg")), "amount": ("0.067 mmol", q(0.067, "mmol", "0.067 mmol"))}), participant("sieves", agent("molecular_sieves", ev(6)), "additive", ev(6), charges={"mass": ("1.4 g", q(1.4, "g", "1.4 g"))}), participant("ch2cl2", agent("ch2cl2", ev(6)), "solvent", ev(6), charges={"volume": ("17 mL", q(17, "mL", "17 mL"))}), participant("product", m1k, "product", ev(6))],
    [step("reaction", ev(6), "add Sc(OTf)3 to 12k, propargyl alcohol and 3 Å molecular sieves in CH2Cl2", ["12k", "propargyl", "sc", "sieves", "ch2cl2"], temperature=("room temperature", "room temperature"), time=("18 h", q(18, "h", "18 h")), atmosphere="Ar"), step("workup", ev(6), "treat with saturated NaHCO3, extract CH2Cl2, wash brine, dry Na2SO4 and evaporate", [], stage="workup"), step("purify", ev(6), "silica chromatography, n-hexane/EtOAc 6:1 to 5:1", [], stage="purification")],
    [outcome("mass", "product_amount", [m1k], ev(6), "530 mg", q(530, "mg", "530 mg"), scope="isolated crystals"), outcome("yield", "yield", [m1k], ev(6), "70% yield", q(70, "%", "70%"), "isolated yield", "silica chromatography")], ["C-O"], solvent="CH2Cl2", scale=("3.38 mmol 12k", {"reference_molecule_id": m12k, "amount": q(3.38, "mmol", "3.38 mmol")}), procedure="Sc(OTf)3 was added to 12k, propargyl alcohol and 3 Å molecular sieves in CH2Cl2; the mixture was stirred 18 h at room temperature under Ar and worked up.")

# Source context: characterization and non-record observations are retained but not promoted.
bundle["source_context"] = [
    {"context_id": f"{DOC}:context_title", "kind": "bibliographic_metadata", "evidence_ids": ev(1),
     "title": "Brønsted Acid-Catalyzed Aza-Ferrier Reaction of N,O-Allenyl Acetals: Synthesis of β-Amino-α-methylene Aldehydes",
     "doi": "10.1021/acs.joc.0c01047"},
    {"context_id": f"{DOC}:context_tbu_substrate_failure", "kind": "unresolved_negative_observation", "evidence_ids": ev(2),
     "text": "Preparation of the tert-butyl-substituted N,O-allenyl substrate (R2=tBu) was unsuccessful; no concrete procedure is supplied on page 2."},
    {"context_id": f"{DOC}:context_general_methods", "kind": "general_methods", "evidence_ids": ev(3),
     "text": "General spectroscopy, inert-atmosphere, low-temperature apparatus, reagent handling, TLC and chromatography statements."},
    {"context_id": f"{DOC}:context_characterization_pages4_6", "kind": "characterization", "evidence_ids": evs(4, 5, 6),
     "text": "NMR, IR, HRMS, melting point, physical appearance and rotamer/diastereomer descriptions for compounds 1a-k and 2a-k; these data support compound identity but are not separate synthetic experiments."},
]

bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
OUT_RAW.parent.mkdir(parents=True, exist_ok=True)
OUT_RAW.write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
OUT_SUBMITTED.write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")

regions = {
    1: [
        {"region": "title/abstract and graphical abstract", "finding": "Current paper title, DOI and two-stage 1→2→3 summary; bibliographic/source context, not a standalone optimization row.", "evidence_ids": ev(1)},
        {"region": "Scheme 1a", "finding": "Prior-work generic aza-Ferrier transformation; classified as literature_example.", "evidence_ids": ev(1)},
        {"region": "Scheme 1b-c", "finding": "Two prior-work examples (46% Ferrier rearrangement; 84% Triton B isomerization), kept as literature_example records.", "evidence_ids": ev(1)},
        {"region": "Scheme 1d", "finding": "Current-work generic sequence 1→2→3; generic_scheme only, with concrete experiments elsewhere.", "evidence_ids": ev(1)},
    ],
    2: [
        {"region": "left-column prose above Scheme 2", "finding": "Low-yield tBuOK isomerization control with qualitative recovery and decomposition.", "evidence_ids": ev(2)},
        {"region": "Scheme 2", "finding": "Twelve distinct 1a-l→2a-l Triton B isomerizations, one record per arrow/product instance.", "evidence_ids": ev(2)},
        {"region": "Table 1 and footnotes", "finding": "Nine distinct acid-screening experiments, including isolated 3a yields, 2a recoveries, side-product yield and the 1.25 mmol scale-up.", "evidence_ids": ev(2)},
        {"region": "right-column scope discussion", "finding": "Narrative cross-check of Scheme 3 outcomes; unsuccessful tBu-substrate preparation retained as source_context because no concrete procedure is supplied.", "evidence_ids": ev(2)},
    ],
    3: [
        {"region": "Scheme 3 and footnotes", "finding": "Eleven standard scope reactions plus separate 2e and 2i over-1 mmol runs; 3j retains 2% yield, 5/5 dr and TLC-observed side products.", "evidence_ids": ev(3)},
        {"region": "Scheme 4a-c", "finding": "Five derivatization experiments: Wittig, Diels-Alder, reduction, cyclization and Michael addition.", "evidence_ids": ev(3)},
        {"region": "left-column discussion below Scheme 3", "finding": "Failed analogous cyclization using acyclic 3a retained as a partial negative-control record.", "evidence_ids": ev(3)},
        {"region": "Experimental Section general methods", "finding": "Instrument and general handling/chromatography information retained as source_context; not synthetic experiments.", "evidence_ids": ev(3)},
        {"region": "bottom-right start of Representative Procedure for 1a", "finding": "Start of Step 1 for 9a; merged with the continuation on page 4.", "evidence_ids": ev(3)},
    ],
    4: [
        {"region": "left column continuation of 1a representative procedure", "finding": "Completed 9a formation, 10a preparation, 1a N-methylation and full 1a→2a isomerization procedure with ordered workup/purification.", "evidence_ids": ev(4)},
        {"region": "right column product paragraphs", "finding": "Concrete 1b/1c N-alkylations and 2b/2c isomerizations; masses/yields merged with Scheme 2 records.", "evidence_ids": ev(4)},
        {"region": "bottom-right start of 1d paragraph", "finding": "Start of overall 1d preparation; merged with continuation on page 5.", "evidence_ids": ev(4)},
        {"region": "NMR/IR/HRMS blocks", "finding": "Characterization only; retained in source_context, not promoted to reactions.", "evidence_ids": ev(4)},
    ],
    5: [
        {"region": "two-column product paragraphs 1d-2h", "finding": "Overall preparations of 1d-h and isomerizations 2d-h; isolated masses/yields captured, exact inherited run charges left unreported.", "evidence_ids": ev(5)},
        {"region": "NMR/IR/HRMS and melting-point lines", "finding": "Characterization only; no separate synthetic records created.", "evidence_ids": ev(5)},
    ],
    6: [
        {"region": "left column 1i/2i and 1j/2j paragraphs", "finding": "9i formation, overall 9i→1i, 1i→2i, overall 1j preparation and 1j→2j isomerization; source-reported mixtures retained.", "evidence_ids": ev(6)},
        {"region": "right column Representative Procedure for 1k", "finding": "Three distinct concrete steps: 11k Boc protection, lactam reduction to 12k, and Sc(OTf)3-catalyzed acetalization to 1k.", "evidence_ids": ev(6)},
        {"region": "bottom-right 2k paragraph", "finding": "1k→2k Triton B isomerization with product mass/yield merged into the Scheme 2 record.", "evidence_ids": ev(6)},
        {"region": "NMR/IR/HRMS blocks and page-5 continuation", "finding": "Characterization only; retained in source_context.", "evidence_ids": ev(6)},
    ],
}
coverage_pages = []
for page in range(1, 7):
    page_eid = f"{DOC}:page{page}:text"
    record_ids = [r["reaction_id"] for r in reactions if page_eid in r["evidence_ids"]]
    unresolved = {
        1: ["Prior-work Scheme 1 structures were not converted into definitive compact structures from the small drawings."],
        2: ["Exact conditions of the tBuOK attempt and unsuccessful R2=tBu substrate preparation are not reported.", "Run-specific charges for Scheme 2 derivatives 2b-l are not shown on this page."],
        3: ["Several Scheme 4 times and measurement bases are absent; dr methods/configurations are not reported.", "The failed acyclic cyclization does not identify the exact cyclization intermediate or charges."],
        4: ["Compound structures are machine candidates and remain needs_review; characterization was not promoted to reaction facts."],
        5: ["The phrase 'same procedure with 1a/2a' does not restate run-specific charges; these remain unreported."],
        6: ["Overall/multistep yields for 1i and 12k retain their stated spans; per-step yields are unavailable.", "SI is mentioned by the paper but no SI file was provided or associated."],
    }[page]
    coverage_pages.append({"page": page, "state": "attempted", "regions": regions[page], "record_ids": record_ids, "unresolved": unresolved})
coverage = {"document_id": DOC, "semantic_accuracy_verified": False, "pages": coverage_pages}
OUT_COVERAGE.parent.mkdir(parents=True, exist_ok=True)
OUT_COVERAGE.write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")

print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]), "raw": str(OUT_RAW), "submitted": str(OUT_SUBMITTED), "coverage": str(OUT_COVERAGE)}, ensure_ascii=False))
