#!/usr/bin/env python3
"""Build the bounded v3.1 extraction for original PDF pages 1--6.

This is a deterministic transcription of the page-image review.  It does not
claim human review or chemical accuracy and deliberately preserves unresolved
structures, source conflicts, and unavailable SI details.
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

RUN = ROOT / "runs/current"
DOC = "doc_975f54c401a61631"
MISS = object()


def ev(*pages: int) -> list[str]:
    ids: list[str] = []
    for page in pages:
        ids.append(f"{DOC}:page{page}:text")
        if 1 <= page <= 6:
            ids.append(f"{DOC}:page{page}:image")
    return list(dict.fromkeys(ids))


def F(raw=None, norm=MISS, pages=(), method="original_extraction", status=None, derivation=None, alternatives=()):
    if norm is MISS:
        norm = raw
    return st.fact(raw, norm, ev(*pages), method=method, status=status,
                   derivation=derivation, alternatives=alternatives)


def inferred(norm, pages, why, raw=None):
    return F(raw, norm, pages, "structure_derivation", derivation=why)


def q(number, unit, qualifier=None, raw_text=None):
    out = {"number": number, "unit": unit}
    if qualifier:
        out["qualifier"] = qualifier
    if raw_text:
        out["raw_text"] = raw_text
    return out


def generic_template(description, variables, pages):
    return {
        "fixed_structure": F(description, {"description": description}, pages),
        "variables": [
            {
                "label": label,
                "attachment_points": F(text, text, pages),
                "definition": F(text, text, pages),
                "allowed_values": F(None, None, pages),
                "constraints": F(None, None, pages),
            }
            for label, text in variables
        ],
        "instance_mapping": [],
        "reason": F("Generic R-group structure in the source drawing.",
                    "Generic R-group structure in the source drawing.", pages),
    }


diagnostics = {}
molecules = []


def add_molecule(mid, label, name, pages, smiles=None, stereo="achiral", template=None,
                 identity_note=None, original_name=None):
    evidence = ev(*pages)
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, pages)
    m["names"]["original"] = F(original_name or name, original_name or name, pages)
    m["names"]["english"] = F(name, name, pages)
    m["names"]["abbreviations"] = F(None, None, pages)
    m["names"]["aliases"] = F(None, None, pages)
    m["cas_numbers"] = F(None, None, pages)
    m["external_properties"] = []
    if template is not None:
        m["representation_type"] = "generic_template"
        m["smiles"] = F(None, None, pages)
        m["generic_template"] = template
        for key in ("atom_centers", "double_bonds", "axial", "mixture", "assignment",
                    "isomeric_smiles", "identity_key"):
            m["stereochemistry"][key] = F(None, None, pages)
        m["stereochemistry"]["chirality_status"] = inferred(
            "unresolved", pages, "R-group template does not specify a single molecular stereostructure.")
        m["stereochemistry"]["geometry_status"] = inferred(
            "unresolved", pages, "R-group template does not specify a single molecular stereostructure.")
        m["stereochemistry"]["representation_status"] = inferred(
            "unresolved", pages, "No single isomeric SMILES can encode the R-group template.")
    elif smiles is None:
        m["representation_type"] = "unresolved"
        m["smiles"] = F(None, None, pages)
        m["generic_template"] = None
        for key in ("atom_centers", "double_bonds", "axial", "mixture", "assignment",
                    "isomeric_smiles", "identity_key"):
            m["stereochemistry"][key] = F(None, None, pages)
        m["stereochemistry"]["chirality_status"] = inferred(
            "unresolved", pages, "The source identity is named but a source-checked complete graph was not established.")
        m["stereochemistry"]["geometry_status"] = inferred(
            "unresolved", pages, "The source identity is named but a source-checked complete graph was not established.")
        m["stereochemistry"]["representation_status"] = inferred(
            "not_encoded", pages, "No concrete source-checked SMILES was assigned.")
    else:
        diag = st.stereo_inspect(smiles)
        diagnostics[mid] = {"input_smiles": smiles, **diag}
        m["representation_type"] = "smiles"
        m["smiles"] = inferred(smiles, pages,
            "Candidate SMILES was transcribed/derived from the inspected source drawing and name; RDKit parsing is diagnostic, not source verification.",
            raw={"source_label": label, "source_name": original_name or name})
        m["generic_template"] = None
        m["stereochemistry"]["atom_centers"] = inferred(
            diag["atom_centers"], pages, "Computed from the candidate SMILES with the bundled standard stereo diagnostic.")
        m["stereochemistry"]["identity_key"] = inferred(
            diag["stereo_graph_key"], pages, "SHA-256 graph key from the canonical isomeric SMILES; not a sample identity.")
        if stereo == "E":
            m["stereochemistry"]["double_bonds"] = F(
                "(E) in the compound name and E geometry in the source drawing",
                [{"configuration": "E", "location": "new allylic-sulfone alkene"}], pages)
            m["stereochemistry"]["axial"] = inferred([], pages, "The product graph contains no allene axis.")
            m["stereochemistry"]["mixture"] = F(None, None, pages)
            m["stereochemistry"]["assignment"] = F("(E)", "E alkene", pages)
            m["stereochemistry"]["chirality_status"] = inferred(
                "achiral", pages, "No tetrahedral or axial stereocenter is present in the candidate product graph.")
            m["stereochemistry"]["geometry_status"] = F("(E)", "specified", pages)
            m["stereochemistry"]["representation_status"] = inferred(
                "encoded", pages, "The candidate SMILES contains explicit alkene direction markers.")
            m["stereochemistry"]["isomeric_smiles"] = inferred(
                diag["canonical_isomeric_smiles"], pages, "Canonicalized from the E-configured candidate SMILES.")
        elif stereo == "EE":
            m["stereochemistry"]["double_bonds"] = F(
                "1,3-Bis((E)-3-(phenylsulfonyl)prop-1-en-1-yl)",
                [{"configuration": "E", "count": 2, "location": "both allylic-sulfone alkenes"}], pages)
            m["stereochemistry"]["axial"] = inferred([], pages, "The product graph contains no allene axis.")
            m["stereochemistry"]["mixture"] = F(None, None, pages)
            m["stereochemistry"]["assignment"] = F("both (E)", "both alkenes E", pages)
            m["stereochemistry"]["chirality_status"] = inferred(
                "achiral", pages, "No tetrahedral or axial stereocenter is present in the candidate product graph.")
            m["stereochemistry"]["geometry_status"] = F("both (E)", "specified", pages)
            m["stereochemistry"]["representation_status"] = inferred(
                "encoded", pages, "Both candidate alkene configurations are encoded.")
            m["stereochemistry"]["isomeric_smiles"] = inferred(
                diag["canonical_isomeric_smiles"], pages, "Canonicalized from the doubly E-configured candidate SMILES.")
        elif stereo == "allene_terminal":
            m["stereochemistry"]["double_bonds"] = inferred(
                [{"type": "cumulated allene", "E_Z_applicable": False}], pages,
                "E/Z descriptors do not apply independently to the cumulated allene bonds.")
            m["stereochemistry"]["axial"] = inferred(
                [{"type": "allene", "stereogenic": False, "reason": "terminal CH2 has two identical H substituents"}],
                pages, "Assessed from the terminal N-allenyl candidate graph.")
            m["stereochemistry"]["mixture"] = F(None, None, pages)
            m["stereochemistry"]["assignment"] = F(None, None, pages)
            m["stereochemistry"]["chirality_status"] = inferred(
                "achiral", pages, "The terminal allene has two identical substituents at one end and no other stereocenter.")
            m["stereochemistry"]["geometry_status"] = inferred(
                "not_applicable", pages, "E/Z is not applicable to the terminal cumulated allene.")
            m["stereochemistry"]["representation_status"] = inferred(
                "not_applicable", pages, "No source-specified stereochemical element requires encoding.")
            m["stereochemistry"]["isomeric_smiles"] = inferred(
                diag["canonical_isomeric_smiles"], pages, "Canonicalized from the candidate connectivity.")
        elif stereo == "racemic":
            m["stereochemistry"]["double_bonds"] = inferred([], pages, "The product contains no alkene.")
            m["stereochemistry"]["axial"] = inferred([], pages, "The product contains no allene or other assigned axis.")
            m["stereochemistry"]["mixture"] = F(
                "(R,S) in the source compound name", {"type": "R/S mixture", "ratio": "unreported"}, pages)
            m["stereochemistry"]["assignment"] = F("(R,S)", "R and S configurations reported", pages)
            m["stereochemistry"]["chirality_status"] = F("(R,S)", "racemic", pages)
            m["stereochemistry"]["geometry_status"] = inferred(
                "not_applicable", pages, "The saturated product has no applicable E/Z element.")
            m["stereochemistry"]["representation_status"] = inferred(
                "partial", pages, "Connectivity SMILES cannot encode a racemic sample as a single tagged stereoisomer.")
            m["stereochemistry"]["isomeric_smiles"] = F(None, None, pages)
        elif stereo == "sulfur_unspecified":
            m["stereochemistry"]["double_bonds"] = inferred(
                diag["double_bonds"], pages, "Computed from the candidate connectivity; no E/Z assignment applies.")
            m["stereochemistry"]["axial"] = inferred([], pages, "No axial assignment is present in the candidate graph.")
            m["stereochemistry"]["mixture"] = F(None, None, pages)
            m["stereochemistry"]["assignment"] = inferred(
                {"center": "sulfinyl sulfur", "configuration": "unassigned"}, pages,
                "The bundled stereo diagnostic flags an unspecified potential sulfur center; the paper gives no sulfur configuration.")
            m["stereochemistry"]["chirality_status"] = inferred(
                "unspecified", pages, "Potential sulfinyl-sulfur stereochemistry is not specified in the source or candidate SMILES.")
            m["stereochemistry"]["geometry_status"] = inferred(
                "not_applicable", pages, "No source-assigned alkene geometry applies to this participant.")
            m["stereochemistry"]["representation_status"] = inferred(
                "partial", pages, "Connectivity is encoded but the potential sulfur configuration remains unspecified.")
            m["stereochemistry"]["isomeric_smiles"] = inferred(
                diag["canonical_isomeric_smiles"], pages, "Canonical connectivity retains the source's unspecified sulfur configuration.")
        else:
            m["stereochemistry"]["double_bonds"] = inferred(
                diag["double_bonds"], pages, "Computed from the candidate SMILES; no source stereochemical assignment was added.")
            m["stereochemistry"]["axial"] = inferred([], pages, "No axial assignment is present in the candidate graph.")
            m["stereochemistry"]["mixture"] = F(None, None, pages)
            m["stereochemistry"]["assignment"] = F(None, None, pages)
            m["stereochemistry"]["chirality_status"] = inferred(
                "achiral", pages, "The candidate graph contains no assigned or potential tetrahedral/axial stereogenic element.")
            m["stereochemistry"]["geometry_status"] = inferred(
                "not_applicable", pages, "No source-assigned alkene geometry applies to this participant.")
            m["stereochemistry"]["representation_status"] = inferred(
                "not_applicable", pages, "No source-specified stereochemical element requires encoding.")
            m["stereochemistry"]["isomeric_smiles"] = inferred(
                diag["canonical_isomeric_smiles"], pages, "Canonicalized from the candidate connectivity.")
    m["identity_note"] = F(identity_note, identity_note, pages) if identity_note else F(None, None, pages)
    molecules.append(m)
    return mid


# Core substrates.
M = {}
M["1a"] = add_molecule(f"{DOC}:mol:1a", "1a", "N-allenyl-2-pyrrolidinone", (1, 2, 3, 4, 5, 6),
                       "O=C1CCCN1C=C=C", "allene_terminal")
M["1b"] = add_molecule(f"{DOC}:mol:1b", "1b", "N-allenyl-2-oxazolidinone", (3, 5, 6),
                       "O=C1OCCN1C=C=C", "allene_terminal")
M["1c"] = add_molecule(f"{DOC}:mol:1c", "1c", "N,N-diallenyl-2-imidazolidinone", (3, 6),
                       "O=C1N(C=C=C)CCN1C=C=C", "allene_terminal")
M["1d"] = add_molecule(f"{DOC}:mol:1d", "1d", "N-allenyl-N-tosylaniline", (3, 6),
                       "C=C=CN(c1ccccc1)S(=O)(=O)c1ccc(C)cc1", "allene_terminal")

# Sulfinate salts 2a--2n. Candidate graphs are retained as derived, not source-verified.
sulfinates = {
    "2a": ("benzenesulfinic acid sodium salt", "[Na+].[O-]S(=O)c1ccccc1"),
    "2b": ("p-toluenesulfinic acid sodium salt", "[Na+].[O-]S(=O)c1ccc(C)cc1"),
    "2c": ("4-fluorobenzenesulfinic acid sodium salt", "[Na+].[O-]S(=O)c1ccc(F)cc1"),
    "2d": ("4-chlorobenzenesulfinic acid sodium salt", "[Na+].[O-]S(=O)c1ccc(Cl)cc1"),
    "2e": ("3-chloro-4-methylbenzenesulfinic acid sodium salt", "[Na+].[O-]S(=O)c1cc(Cl)c(C)cc1"),
    "2f": ("3-(trifluoromethyl)benzenesulfinic acid sodium salt", "[Na+].[O-]S(=O)c1cccc(C(F)(F)F)c1"),
    "2g": ("pyridine-2-sulfinic acid sodium salt", "[Na+].[O-]S(=O)c1ccccn1"),
    "2h": ("2,5-dichlorothiophene-3-sulfinic acid sodium salt", None),
    "2i": ("methanesulfinic acid sodium salt", "[Na+].CS(=O)[O-]"),
    "2j": ("isopropylsulfinic acid sodium salt", "[Na+].CC(C)S(=O)[O-]"),
    "2k": ("cyclopropylsulfinic acid sodium salt", "[Na+].[O-]S(=O)C1CC1"),
    "2l": ("cyclohexylsulfinic acid sodium salt", "[Na+].[O-]S(=O)C1CCCCC1"),
    "2m": ("sodium triflinate", "[Na+].[O-]S(=O)C(F)(F)F"),
    "2n": ("sodium 1-methyl 3-sulfinopropanoate", "[Na+].[O-]S(=O)CCC(=O)OC"),
}
for label, (name, smi) in sulfinates.items():
    M[label] = add_molecule(f"{DOC}:mol:{label}", label, name, (2, 3, 4, 5, 6), smi,
                            "sulfur_unspecified",
                            identity_note=("The heteroaryl connectivity was not forced into a guessed SMILES."
                                           if smi is None else None))

# Product candidates. Every named (E) product remains a separate molecule identity.
core_a = "O=C1CCCN1/C=C/CS(=O)(=O){}"
prod_a = {
    "3aa": ("(E)-1-(3-(phenylsulfonyl)prop-1-en-1-yl)pyrrolidin-2-one", core_a.format("c1ccccc1")),
    "3ab": ("(E)-1-(3-tosylprop-1-en-1-yl)pyrrolidin-2-one", core_a.format("c1ccc(C)cc1")),
    "3ac": ("(E)-1-(3-((4-fluorophenyl)sulfonyl)prop-1-en-1-yl)pyrrolidin-2-one", core_a.format("c1ccc(F)cc1")),
    "3ad": ("(E)-1-(3-((4-chlorophenyl)sulfonyl)prop-1-en-1-yl)pyrrolidin-2-one", core_a.format("c1ccc(Cl)cc1")),
    "3ae": ("(E)-1-(3-((3-chloro-4-methylphenyl)sulfonyl)prop-1-en-1-yl)pyrrolidin-2-one", core_a.format("c1cc(Cl)c(C)cc1")),
    "3af": ("(E)-1-(3-((3-(trifluoromethyl)phenyl)sulfonyl)prop-1-en-1-yl)pyrrolidin-2-one", core_a.format("c1cccc(C(F)(F)F)c1")),
    "3ag": ("(E)-1-(3-(pyridin-2-ylsulfonyl)prop-1-en-1-yl)pyrrolidin-2-one", core_a.format("c1ccccn1")),
    "3ah": ("(E)-1-(3-((2,5-dichlorothiophen-3-yl)sulfonyl)prop-1-en-1-yl)pyrrolidin-2-one", None),
    "3ai": ("(E)-1-(3-(methylsulfonyl)prop-1-en-1-yl)pyrrolidin-2-one", core_a.format("C")),
    "3aj": ("(E)-1-(3-(isopropylsulfonyl)prop-1-en-1-yl)pyrrolidin-2-one", core_a.format("C(C)C")),
    "3ak": ("(E)-1-(3-(cyclopropylsulfonyl)prop-1-en-1-yl)pyrrolidin-2-one", core_a.format("C1CC1")),
    "3al": ("(E)-1-(3-(cyclohexylsulfonyl)prop-1-en-1-yl)pyrrolidin-2-one", core_a.format("C1CCCCC1")),
    "3am": ("(E)-1-(3-((trifluoromethyl)sulfonyl)prop-1-en-1-yl)pyrrolidin-2-one", core_a.format("C(F)(F)F")),
    "3an": ("methyl (E)-3-((3-(2-oxopyrrolidin-1-yl)allyl)sulfonyl)propanoate", core_a.format("CCC(=O)OC")),
}
product_pages = {**{x: (2, 4) for x in ("3aa", "3ab", "3ac", "3ad", "3ae", "3af")},
                 **{x: (2, 5) for x in ("3ag", "3ah", "3ai", "3aj", "3ak", "3al", "3am", "3an")}}
for label, (name, smi) in prod_a.items():
    M[label] = add_molecule(f"{DOC}:mol:{label}", label, name, product_pages[label], smi, "E",
                            identity_note=("The heteroaryl product graph was left unresolved rather than guessing the thiophene substitution order."
                                           if smi is None else None))

core_b = "O=C1OCCN1/C=C/CS(=O)(=O){}"
prod_b = {
    "3ba": ("(E)-3-(3-(phenylsulfonyl)prop-1-en-1-yl)oxazolidin-2-one", core_b.format("c1ccccc1")),
    "3bb": ("(E)-3-(3-tosylprop-1-en-1-yl)oxazolidin-2-one", core_b.format("c1ccc(C)cc1")),
    "3bi": ("source-conflicted 3bi methyl/isopropylsulfonyl oxazolidinone", core_b.format("C")),
    "3bj": ("(E)-3-(3-(isopropylsulfonyl)prop-1-en-1-yl)oxazolidin-2-one", core_b.format("C(C)C")),
    "3bm": ("(E)-3-(3-((trifluoromethyl)sulfonyl)prop-1-en-1-yl)oxazolidin-2-one", core_b.format("C(F)(F)F")),
}
for label, (name, smi) in prod_b.items():
    note = None
    original = None
    if label == "3bi":
        original = "(E)-3-(3-(Isopropylsulfonyl)prop-1-en-1-yl)oxazolidin-2-one"
        note = ("Source conflict: the page-5 heading says isopropylsulfonyl, while the Scheme 3 drawing, "
                "reactant 2i (methanesulfinate), and formula C7H11NO4S support a methylsulfonyl graph. "
                "The candidate SMILES follows the drawing/formula and remains needs_review.")
    M[label] = add_molecule(f"{DOC}:mol:{label}", label, name, (3, 5, 6), smi, "E",
                            identity_note=note, original_name=original)

M["3ca"] = add_molecule(f"{DOC}:mol:3ca", "3ca",
    "1,3-bis((E)-3-(phenylsulfonyl)prop-1-en-1-yl)imidazolidin-2-one", (3, 6),
    "O=C1N(/C=C/CS(=O)(=O)c2ccccc2)CCN1/C=C/CS(=O)(=O)c1ccccc1", "EE")
M["3da"] = add_molecule(f"{DOC}:mol:3da", "3da",
    "(E)-4-methyl-N-phenyl-N-(3-(phenylsulfonyl)prop-1-en-1-yl)benzenesulfonamide", (3, 6),
    "Cc1ccc(S(=O)(=O)N(c2ccccc2)/C=C/CS(=O)(=O)c2ccccc2)cc1", "E")
M["3db"] = add_molecule(f"{DOC}:mol:3db", "3db",
    "(E)-4-methyl-N-phenyl-N-(3-tosylprop-1-en-1-yl)benzenesulfonamide", (3, 6),
    "Cc1ccc(S(=O)(=O)N(c2ccccc2)/C=C/CS(=O)(=O)c2ccc(C)cc2)cc1", "E")
M["4aa"] = add_molecule(f"{DOC}:mol:4aa", "4aa",
    "(R,S)-1-(1,3-bis(phenylsulfonyl)propyl)pyrrolidin-2-one", (2, 3, 6),
    "O=C1CCCN1C(S(=O)(=O)c1ccccc1)CCS(=O)(=O)c1ccccc1", "racemic")
M["d3aa"] = add_molecule(f"{DOC}:mol:d-3aa", "d-3aa",
    "(E)-1-(3-(phenylsulfonyl)prop-1-en-1-yl-2-d)pyrrolidin-2-one", (3, 6), None,
    identity_note="The E geometry and central deuterium are explicit, but no concrete isotope/isomeric SMILES was forced.")

# Named reagents, solvents, and workup agents.
simple = {
    "TFA": ("trifluoroacetic acid", "O=C(O)C(F)(F)F", (1, 2, 3, 4, 5, 6)),
    "TFA-d": ("deuterated trifluoroacetic acid", "O=C(O[2H])C(F)(F)F", (3, 6)),
    "H2O": ("water", "O", (1, 2, 3, 4)),
    "D2O": ("deuterium oxide", "[2H]O[2H]", (3, 6)),
    "THF": ("tetrahydrofuran", "C1CCOC1", (2, 3, 4, 6)),
    "DCM": ("dichloromethane", "ClCCl", (4, 5, 6)),
    "MgSO4": ("anhydrous magnesium sulfate", "[Mg+2].[O-]S(=O)(=O)[O-]", (4,)),
    "EtOAc": ("ethyl acetate", "CCOC(=O)C", (5,)),
    "hexane": ("n-hexane", "CCCCCC", (5,)),
    "PhSO2H": ("benzenesulfinic acid", "O=S(O)c1ccccc1", (3,)),
    "p-toluic": ("p-toluic acid", "Cc1ccc(C(=O)O)cc1", (2,)),
    "TMBA": ("2,4,6-trimethylbenzoic acid", "Cc1cc(C)c(C(=O)O)c(C)c1", (2,)),
    "PTSA": ("p-toluenesulfonic acid", "Cc1ccc(S(=O)(=O)O)cc1", (2,)),
    "Bpy": ("2,2'-bipyridine", "n1ccccc1-c1ccccn1", (2,)),
    "Phen": ("1,10-phenanthroline", "c1cnc2c(c1)ccc1cccnc12", (2,)),
    "TMEDA": ("N,N,N',N'-tetramethylethylenediamine", "CN(C)CCN(C)C", (2,)),
}
for key, (name, smi, pages) in simple.items():
    M[key] = add_molecule(f"{DOC}:mol:{key}", key, name, pages, smi,
                          "sulfur_unspecified" if key == "PhSO2H" else "achiral")
for key, name, pages in [
    ("Cu", "tetrakis(acetonitrile)copper(I) hexafluorophosphate", (2,)),
    ("silica", "triethylamine-treated silica gel", (4, 5)),
    ("trichloroethylene", "trichloroethylene internal standard", (2, 4)),
]:
    M[key] = add_molecule(f"{DOC}:mol:{key}", key, name, pages, None)

# R-group entities used only for the generic/current and literature arrows.
templates = {
    "gen_allene": ("generic N-allenyl derivative 1x", "N-substituted terminal allene", [("R1", "substituent on N"), ("R2", "substituent on N")]),
    "gen_sulfinate": ("generic sodium sulfinate 2y", "R3-SO2Na", [("R3", "(hetero)aryl or (cyclo)alkyl")]),
    "gen_product": ("generic (E)-allylic sulfone 3xy", "N-substituted allylic sulfone", [("R1", "substituent on N"), ("R2", "substituent on N"), ("R3", "sulfone substituent")]),
    "lit_allene": ("generic allene", "R1-substituted allene", [("R1", "allene substituent")]),
    "lit_hydrazine": ("generic sulfonylhydrazine", "R2SO2NHNH2", [("R2", "sulfonyl substituent")]),
    "lit_thiol": ("generic thiol", "R2SH", [("R2", "thiol substituent")]),
    "lit_linear": ("generic linear allylic sulfone", "linear allylic sulfone", [("R1", "allene-derived substituent"), ("R2", "sulfone substituent")]),
    "lit_branched": ("generic branched allylic sulfone", "branched allylic sulfone", [("R1", "allene-derived substituent"), ("R2", "sulfone substituent")]),
    "lit_chiral": ("generic chiral branched allylic sulfone", "source-drawn chiral branched allylic sulfone", [("R1", "allene-derived substituent"), ("R2", "sulfone substituent")]),
    "lit_cyclopropyl": ("generic 3-cyclopropylideneprop-2-en-1-one", "cyclopropylidene allene ketone", [("R1", "ketone substituent"), ("R2", "cyclopropyl substituent")]),
    "lit_cyclopropyl_prod": ("generic Miao-Ren allylic sulfone", "cyclopropyl allylic sulfone", [("R1", "ketone substituent"), ("R2", "cyclopropyl substituent"), ("R3", "sulfone substituent")]),
    "lit_ketoallene": ("generic alpha-ketoallene", "alpha-ketoallene", [("R1", "ketone substituent")]),
    "lit_sulfinic": ("generic sulfinic acid", "R2SO2H", [("R2", "sulfonyl substituent")]),
    "lit_vinyl": ("generic vinylic sulfone", "alpha-keto vinylic sulfone", [("R1", "ketone substituent"), ("R2", "sulfone substituent")]),
}
for key, (name, desc, variables) in templates.items():
    M[key] = add_molecule(f"{DOC}:mol:{key}", key, name, (1, 2, 3),
                          template=generic_template(desc, variables, (1, 2, 3)))

for key, name, pages in [
    ("PdOAc2", "palladium(II) acetate", (1,)), ("DPPF", "DPPF", (1,)),
    ("AcOH", "acetic acid", (1,)), ("RhCODCl", "[Rh(COD)Cl]2", (1,)),
    ("DPEphos", "DPEphos", (1,)), ("PhCO2H", "benzoic acid", (1,)),
    ("Difluorphos", "(R)-Difluorphos", (1,)), ("EtOH", "ethanol", (1,)),
    ("mCPBA", "m-CPBA", (1,)), ("DMSO", "dimethyl sulfoxide", (1,)),
]:
    M[key] = add_molecule(f"{DOC}:mol:{key}", key, name, pages, None)


def participant(pid, mid, role, pages, stage="reaction", amount=None, mass=None, volume=None,
                equiv=None, mol_percent=None, reference=None, product_role=None):
    p = st.new_participant(pid, mid, role, ev(*pages), stage)
    if product_role is not None:
        p["product_role"] = F(product_role, product_role, pages)
    if amount is not None:
        p["charges"]["amount"] = F(amount[0], amount[1], pages)
    if mass is not None:
        p["charges"]["mass"] = F(mass[0], mass[1], pages)
    if volume is not None:
        p["charges"]["volume"] = F(volume[0], volume[1], pages)
    if equiv is not None:
        p["charges"]["equivalents"] = equiv if isinstance(equiv, dict) and "raw_value" in equiv else F(equiv[0], equiv[1], pages)
    if mol_percent is not None:
        p["charges"]["mol_percent"] = F(mol_percent[0], mol_percent[1], pages)
    if reference is not None:
        p["equivalence_reference"] = F(reference, reference, pages)
    return p


def step(step_id, stage, pids, pages, operation=None, temperature=None, time=None,
         atmosphere=None, pressure=None):
    s = st.new_step(step_id, ev(*pages), stage)
    s["participant_ids"] = list(pids)
    if operation is not None:
        s["operation"] = F(operation, operation, pages)
    if temperature is not None:
        s["temperature"] = F(temperature[0], temperature[1], pages)
    if time is not None:
        s["time"] = F(time[0], time[1], pages)
    if atmosphere is not None:
        s["atmosphere"] = F(atmosphere, atmosphere, pages)
    if pressure is not None:
        s["pressure"] = F(pressure, pressure, pages)
    return s


def reaction(rid, pages, record_type="specific_experiment", categories=(), location=None,
             reaction_name=None, transformation=None):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev(*pages))
    r["record_type"] = record_type
    r["experiment_categories"] = F(list(categories), list(categories), pages)
    r["location"] = F(location or f"PDF pages {', '.join(map(str, pages))}",
                      {"pdf_pages": list(pages), "source_location": location}, pages)
    if reaction_name:
        r["reaction_name"] = F(reaction_name, reaction_name, pages)
    if transformation:
        r["transformation_type"] = inferred(transformation, pages,
            "Transformation description derived from the source arrow and named reaction.")
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    return r


def set_solvent(r, pages, composition=None, ratio=None, volume=None, concentration=None, basis=None):
    if composition is not None:
        r["solvent_conditions"]["composition"] = F(composition, composition, pages)
    if ratio is not None:
        r["solvent_conditions"]["mixture_ratio"] = F(ratio, ratio, pages)
    if volume is not None:
        r["solvent_conditions"]["total_volume"] = F(volume[0], volume[1], pages)
    if concentration is not None:
        r["solvent_conditions"]["concentration"] = F(concentration[0], concentration[1], pages)
    if basis is not None:
        r["solvent_conditions"]["concentration_basis"] = F(basis, basis, pages)


def set_procedure(r, pages, text, general=None, modifications=None, addition_order=None):
    r["procedure"]["original_text"] = F(text, text, pages)
    if general:
        r["procedure"]["general_procedure_reference"] = F(general, general, pages)
    if modifications:
        r["procedure"]["modifications"] = F(modifications, modifications, pages)
    if addition_order:
        r["procedure"]["addition_order"] = F(addition_order, addition_order, pages)


def add_outcome(r, oid, kind, targets, pages, raw, normalized, basis=None,
                yield_type=None, role=None, scope=None):
    o = st.new_outcome(oid, kind, targets, ev(*pages))
    o["value"] = F(raw, normalized, pages)
    if basis is not None:
        o["measurement_basis"] = F(basis, basis, pages)
    if yield_type is not None:
        o["yield_type"] = F(yield_type, yield_type, pages)
    if role is not None:
        o["product_role"] = F(role, role, pages)
    if scope is not None:
        o["scope"] = F(scope, scope, pages)
    r["outcomes"].append(o)


def set_classification(r, pages, cs=1, ch=1, applicable=True):
    if not applicable:
        r["classification"]["mapped_reaction_smiles"] = F(None, None, pages)
        r["classification"]["bond_changes"] = []
        r["classification"]["labels"] = inferred([], pages, "No product-forming transformation was observed.")
        r["classification"]["coverage"] = inferred("not_applicable", pages,
            "No product-forming bond changes are assigned to the negative/no-product experiment.")
        r["classification"]["limitations"] = ["No bond-formation label inferred from an intended but unobserved reaction."]
        return
    changes = []
    n = 1
    for _ in range(cs):
        changes.append({"change_id": f"b{n}", "atoms": [{"atom_map": None, "element": "C"},
                        {"atom_map": None, "element": "S"}], "before_order": 0,
                        "after_order": 1, "change_type": "formed"})
        n += 1
    for _ in range(ch):
        changes.append({"change_id": f"b{n}", "atoms": [{"atom_map": None, "element": "C"},
                        {"atom_map": None, "element": "H"}], "before_order": 0,
                        "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1})
        n += 1
    labels = {}
    for c in changes:
        lab = st.bond_label(c["atoms"][0]["element"], c["atoms"][1]["element"], c["after_order"])
        labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(c["change_id"])
    why = ("Partial transformation-level differences transcribed from the hydrosulfonylation arrow. "
           "No atom-mapped reaction SMILES or complete atom correspondence is asserted.")
    r["classification"]["mapped_reaction_smiles"] = F(None, None, pages)
    r["classification"]["bond_changes"] = [inferred(c, pages, why) for c in changes]
    r["classification"]["labels"] = inferred(list(labels.values()), pages, why)
    r["classification"]["coverage"] = inferred("partial", pages, why)
    r["classification"]["limitations"] = [
        "Transformation-level C-S/C-H changes only; atom mapping is absent.",
        "Allene bond-order changes and atom/proton provenance are not claimed complete.",
    ]


reactions = []

# Scheme 1, equations 1--5: explicitly prior literature, never current-paper experiments.
lit_rows = [
    ("scheme1_eq1_yamamoto", "Yamamoto (1998), Scheme 1 eq 1", "lit_hydrazine", "lit_linear",
     [("PdOAc2", "catalyst", ("10 mol%", q(10, "mol%"))), ("DPPF", "ligand", ("10 mol%", q(10, "mol%"))),
      ("AcOH", "reagent", ("1 equiv", q(1, "equiv"))), ("THF", "solvent", None)],
     ("70 °C", q(70, "°C")), ("1 h", q(1, "h")), "THF"),
    ("scheme1_eq2_breit", "Breit (2017), Scheme 1 eq 2", "lit_hydrazine", "lit_branched",
     [("RhCODCl", "catalyst", ("2.5 mol%", q(2.5, "mol%"))), ("DPEphos", "ligand", ("10 mol%", q(10, "mol%"))),
      ("PhCO2H", "reagent", ("5 equiv", q(5, "equiv")))],
     ("80 °C", q(80, "°C")), ("18 h", q(18, "h")), "toluene (0.4 M)"),
    ("scheme1_eq3_breit", "Breit (2015), Scheme 1 eq 3", "lit_thiol", "lit_chiral",
     [("RhCODCl", "catalyst", ("1 mol%", q(1, "mol%"))), ("Difluorphos", "ligand", ("2 mol%", q(2, "mol%"))),
      ("EtOH", "solvent", None), ("mCPBA", "reagent", ("2.5 equiv", q(2.5, "equiv"))), ("DCM", "solvent", None)],
     ("-15 °C then 0 °C", None), ("16 h then 1 h", None), "EtOH then DCM"),
    ("scheme1_eq4_miao_ren", "Miao-Ren (2016), Scheme 1 eq 4", "gen_sulfinate", "lit_cyclopropyl_prod",
     [("AcOH", "reagent", ("2 equiv", q(2, "equiv"))), ("DMSO", "solvent", None)],
     ("RT", "RT"), ("8-72 h", {"minimum": 8, "maximum": 72, "unit": "h", "raw_text": "8-72 h"}), "DMSO (0.1 M)"),
    ("scheme1_eq5_loh", "Loh (2021), Scheme 1 eq 5", "lit_sulfinic", "lit_vinyl",
     [("H2O", "solvent", None), ("EtOH", "solvent", None)],
     ("RT", "RT"), ("1 h", q(1, "h")), "H2O/EtOH 1:1 (0.15 M)"),
]
for idx, (rid, loc, partner, product, extras, temp, tm, solvent_text) in enumerate(lit_rows, 1):
    r = reaction(rid, (1,), "literature_example", ("other",), loc,
                 "literature hydrosulfonylation/hydrothiolation example", "allene hydrofunctionalization")
    substrate = "lit_cyclopropyl" if idx == 4 else "lit_ketoallene" if idx == 5 else "lit_allene"
    r["participants"] = [
        participant("p_allene", M[substrate], "reactant", (1,)),
        participant("p_partner", M[partner], "reactant", (1,)),
        participant("p_product", M[product], "product", (1,), product_role="reported literature product"),
    ]
    for n, (key, role, loading) in enumerate(extras, 1):
        kwargs = {}
        if loading:
            if loading[1]["unit"] == "mol%": kwargs.update(mol_percent=loading, reference=M[substrate])
            elif loading[1]["unit"] == "equiv": kwargs.update(equiv=loading, reference=M[substrate])
        r["participants"].append(participant(f"p_extra{n}", M[key], role, (1,), **kwargs))
    set_solvent(r, (1,), solvent_text)
    if idx == 3:
        r["steps"] = [
            step("s1_hydrothiolation", "reaction", [p["participant_id"] for p in r["participants"][:5]], (1,),
                 "rhodium-catalyzed hydrothiolation", ("-15 °C", q(-15, "°C")), ("16 h", q(16, "h"))),
            step("s2_oxidation", "reaction", ["p_product", "p_extra4", "p_extra5"], (1,),
                 "m-CPBA oxidation", ("0 °C", q(0, "°C")), ("1 h", q(1, "h"))),
        ]
    else:
        r["steps"] = [step("s1", "reaction", [p["participant_id"] for p in r["participants"]], (1,),
                           "conditions printed above Scheme 1 arrow", temp, tm)]
    set_procedure(r, (1,), "Scheme 1 literature-example arrow and printed conditions.")
    set_classification(r, (1,))
    r["review"]["unresolved"] = ["Generic R groups prevent concrete structures and complete atom mapping.",
                                     "No yield is printed in Scheme 1 for this literature example."]
    reactions.append(r)

# One current-paper generic scheme consolidates repeated generic arrows on pages 1--3.
r = reaction("generic_hydrosulfonylation", (1, 2, 3), "generic_scheme", ("substrate_scope",),
             "Scheme 1 eq 6; Scheme 2/3 generic arrows", "TFA-mediated hydrosulfonylation of N-allenyl derivatives",
             "allene hydrosulfonylation")
r["participants"] = [
    participant("p_allene", M["gen_allene"], "reactant", (1, 2, 3), equiv=("1 equiv", q(1, "equiv")), reference=M["gen_allene"]),
    participant("p_sulfinate", M["gen_sulfinate"], "reactant", (1, 2, 3), equiv=("2 equiv", q(2, "equiv")), reference=M["gen_allene"]),
    participant("p_tfa", M["TFA"], "reagent", (1, 2, 3), equiv=("1 equiv", q(1, "equiv")), reference=M["gen_allene"]),
    participant("p_water", M["H2O"], "solvent", (1, 2, 3)),
    participant("p_product", M["gen_product"], "product", (1, 2, 3), product_role="generic (E)-allylic sulfone"),
]
set_solvent(r, (1, 2, 3), "H2O", volume=("0.5 mL", q(0.5, "mL")),
            concentration=("0.5 M", q(0.5, "M")), basis="allene 1")
r["steps"] = [step("s1", "reaction", [p["participant_id"] for p in r["participants"]], (1, 2, 3),
                   "stir under open air", ("25 °C", q(25, "°C")), ("4 h", q(4, "h")), "open air")]
set_procedure(r, (1, 2, 3), "Generic conditions: allene 1, sodium sulfinate 2, TFA, H2O, open air, 25 °C, 4 h.")
set_classification(r, (1, 2, 3))
r["review"]["unresolved"] = ["Generic R1/R2/R3 identities are intentionally retained as variables."]
reactions.append(r)
GENERIC_ID = r["reaction_id"]

# Table 1, nine independent optimization experiments.
table1 = [
    (1, 0, 0, None, None, 0, 0),
    (2, 10, 0, None, "p-toluic", 63, 50),
    (3, 10, 10, "Bpy", "p-toluic", 86, 57),
    (4, 10, 10, "Phen", "p-toluic", 74, 33),
    (5, 10, 10, "TMEDA", "p-toluic", 95, 51),
    (6, 10, 10, "Bpy", "TMBA", 79, 39),
    (7, 10, 10, "Bpy", "PTSA", 76, 30),
    (8, 10, 10, "Bpy", "TFA", ">99", 58),
    (9, 0, 0, None, "TFA", ">99", 90),
]
for entry, x, y, ligand, acid, conversion, yld in table1:
    r = reaction(f"table1_entry{entry}", (2,), "specific_experiment", ("optimization",),
                 f"Table 1 entry {entry}", "Copper-catalyzed hydrosulfonylation of 1a with 2a",
                 "allene hydrosulfonylation")
    r["general_scheme_ids"] = [GENERIC_ID]
    r["participants"] = [
        participant("p_1a", M["1a"], "reactant", (2,), amount=("0.25 mmol", q(0.25, "mmol")),
                    equiv=("1 equiv", q(1, "equiv")), reference=M["1a"]),
        participant("p_2a", M["2a"], "reactant", (2,), amount=("0.5 mmol", q(0.5, "mmol")),
                    equiv=("2 equiv", q(2, "equiv")), reference=M["1a"]),
        participant("p_thf", M["THF"], "solvent", (2,), volume=("0.5 mL", q(0.5, "mL"))),
        participant("p_3aa", M["3aa"], "product", (2,), product_role="target product"),
    ]
    if x:
        r["participants"].append(participant("p_cu", M["Cu"], "catalyst", (2,),
                                             mol_percent=(f"{x} mol%", q(x, "mol%")), reference=M["1a"]))
    if ligand:
        r["participants"].append(participant("p_ligand", M[ligand], "ligand", (2,),
                                             mol_percent=(f"{y} mol%", q(y, "mol%")), reference=M["1a"]))
    if acid:
        r["participants"].append(participant("p_acid", M[acid], "reagent", (2,),
                                             amount=("0.25 mmol", q(0.25, "mmol")),
                                             equiv=("1 equiv", q(1, "equiv")), reference=M["1a"]))
    r["reaction_scale"] = F("1a (0.25 mmol)", {"reference_molecule_id": M["1a"], "amount": q(0.25, "mmol")}, (2,))
    set_solvent(r, (2,), "THF", volume=("0.5 mL", q(0.5, "mL")), concentration=("0.5 M", q(0.5, "M")), basis="1a")
    r["steps"] = [step("s1", "reaction", [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"],
                       (2,), "catalyst and ligand placed in a Schlenk tube under argon; reaction mixture held",
                       ("55 °C", q(55, "°C")), ("18 h", q(18, "h")), "argon")]
    set_procedure(r, (2,), "Table 1 footnote a conditions and the selected entry row.")
    conv_norm = q(99, "%", "gt", ">99%") if conversion == ">99" else q(conversion, "%")
    add_outcome(r, "o_conversion", "conversion", [M["1a"]], (2,), f"{conversion}%", conv_norm,
                basis="based on recovered starting material 1a", scope="conversion of 1a")
    add_outcome(r, "o_yield", "yield", [M["3aa"]], (2,), f"{yld}%", q(yld, "%"),
                basis="1H NMR using trichloroethylene (0.25 mmol, 1 equiv) as internal standard",
                yield_type="NMR yield", role="target product")
    set_classification(r, (2,), applicable=(yld > 0))
    r["review"]["unresolved"] = ([] if yld > 0 else
        ["No product-forming bond changes are assigned to the zero-yield/no-conversion row."])
    reactions.append(r)

# Table 2, five independent parametric-study experiments.
table2 = [
    (1, "THF instead of H2O", "THF", "open air", 25, True, ">99", 80),
    (2, "under argon instead of open air", "H2O", "argon", 25, True, ">99", 90),
    (3, "55 °C instead of 25 °C", "H2O", "open air", 55, True, ">99", 96),
    (4, "no TFA", "H2O", "open air", 25, False, 66, 0),
    (5, "none (standard conditions)", "H2O", "open air", 25, True, ">99", ">99"),
]
for entry, variation, solvent, atmosphere, temp_num, use_tfa, conversion, yld in table2:
    r = reaction(f"table2_entry{entry}", (2,), "specific_experiment", ("optimization",),
                 f"Table 2 entry {entry}: {variation}", "Acid-mediated hydrosulfonylation of 1a with 2a",
                 "allene hydrosulfonylation")
    r["general_scheme_ids"] = [GENERIC_ID]
    r["participants"] = [
        participant("p_1a", M["1a"], "reactant", (2,), amount=("0.25 mmol", q(0.25, "mmol")),
                    equiv=("1 equiv", q(1, "equiv")), reference=M["1a"]),
        participant("p_2a", M["2a"], "reactant", (2,), amount=("0.5 mmol", q(0.5, "mmol")),
                    equiv=("2 equiv", q(2, "equiv")), reference=M["1a"]),
        participant("p_solvent", M[solvent], "solvent", (2,), volume=("0.5 mL", q(0.5, "mL"))),
        participant("p_3aa", M["3aa"], "product", (2,), product_role="target product"),
    ]
    if use_tfa:
        r["participants"].append(participant("p_tfa", M["TFA"], "reagent", (2,),
                                             amount=("0.25 mmol", q(0.25, "mmol")),
                                             equiv=("1 equiv", q(1, "equiv")), reference=M["1a"]))
    r["reaction_scale"] = F("1a (0.25 mmol)", {"reference_molecule_id": M["1a"], "amount": q(0.25, "mmol")}, (2,))
    set_solvent(r, (2,), solvent, volume=("0.5 mL", q(0.5, "mL")), concentration=("0.5 M", q(0.5, "M")), basis="1a")
    r["steps"] = [step("s1", "reaction", [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"],
                       (2,), "placed in a tube and held under the entry conditions",
                       (f"{temp_num} °C", q(temp_num, "°C")), ("4 h", q(4, "h")), atmosphere)]
    set_procedure(r, (2,), "Table 2 standard conditions modified only as stated in the selected row.")
    conv_norm = q(99, "%", "gt", ">99%") if conversion == ">99" else q(conversion, "%")
    y_norm = q(99, "%", "gt", ">99%") if yld == ">99" else q(yld, "%")
    add_outcome(r, "o_conversion", "conversion", [M["1a"]], (2,), f"{conversion}%", conv_norm,
                basis="based on recovered starting material", scope="conversion of 1a")
    add_outcome(r, "o_yield", "yield", [M["3aa"]], (2,), f"{yld}%", y_norm,
                basis="1H NMR using trichloroethylene (0.25 mmol, 1 equiv) as internal standard",
                yield_type="NMR yield", role="target product")
    set_classification(r, (2,), applicable=(yld != 0))
    if yld == 0:
        r["review"]["unresolved"] = ["66% conversion with 0% target-product yield; identities/basis of consumed material are unreported."]
    reactions.append(r)


def general_scope_record(label, substrate, sulfinate, product, product_page, mass_mg, yield_pct,
                         tfa_equiv=1, extra=None):
    pages = (2 if substrate == "1a" else 3, 4, product_page)
    r = reaction(f"scope_{label}", tuple(dict.fromkeys(pages)), "specific_experiment", ("substrate_scope",),
                 f"Scheme {'2' if substrate == '1a' else '3'} product {label}; experimental paragraph page {product_page}",
                 "TFA-mediated hydrosulfonylation", "allene hydrosulfonylation")
    r["general_scheme_ids"] = [GENERIC_ID]
    r["participants"] = [
        participant(f"p_{substrate}", M[substrate], "reactant", pages, amount=("0.25 mmol", q(0.25, "mmol")),
                    equiv=("1 equiv", q(1, "equiv")), reference=M[substrate]),
        participant(f"p_{sulfinate}", M[sulfinate], "reactant", pages, amount=("0.5 mmol", q(0.5, "mmol")),
                    equiv=("2 equiv", q(2, "equiv")), reference=M[substrate]),
        participant("p_tfa", M["TFA"], "reagent", pages,
                    amount=((f"{0.25*tfa_equiv:g} mmol", q(0.25*tfa_equiv, "mmol")) if tfa_equiv == 1 else None),
                    equiv=(f"{tfa_equiv} equiv", q(tfa_equiv, "equiv")), reference=M[substrate]),
        participant("p_water", M["H2O"], "solvent", pages, volume=("0.5 mL", q(0.5, "mL"))),
        participant("p_dcm", M["DCM"], "workup_agent", (4,), stage="workup"),
        participant("p_mgso4", M["MgSO4"], "workup_agent", (4,), stage="workup"),
        participant(f"p_{product}", M[product], "product", pages, product_role="main product"),
    ]
    r["reaction_scale"] = F(f"{substrate} (0.25 mmol)",
        {"reference_molecule_id": M[substrate], "amount": q(0.25, "mmol")}, pages)
    set_solvent(r, pages, "H2O", volume=("0.5 mL", q(0.5, "mL")),
                concentration=("0.5 M", q(0.5, "M")), basis=substrate)
    r["steps"] = [
        step("s1_reaction", "reaction", [f"p_{substrate}", f"p_{sulfinate}", "p_tfa", "p_water", f"p_{product}"],
             pages, "stirred under open air", ("25 °C", q(25, "°C")), ("4 h", q(4, "h")), "open air"),
        step("s2_workup", "workup", ["p_dcm", "p_mgso4", f"p_{product}"], (4,),
             "aqueous workup; separate organic phase; extract aqueous layer with DCM; combine, dry over anhydrous MgSO4, and concentrate under vacuum"),
    ]
    set_procedure(r, pages,
        f"The reaction of {substrate} (0.25 mmol) with {sulfinate} (2 equiv, 0.5 mmol) afforded {product} ({mass_mg} mg, {yield_pct}%) according to the general procedure.",
        general="General Procedure on PDF page 4",
        modifications=(extra if extra else (f"{tfa_equiv} equiv TFA" if tfa_equiv != 1 else None)),
        addition_order="sodium sulfinate, solvent, allene, then TFA")
    add_outcome(r, "o_mass", "product_amount", [M[product]], (product_page,), f"{mass_mg} mg", q(mass_mg, "mg"), role="main product")
    add_outcome(r, "o_yield", "yield", [M[product]], pages, f"{yield_pct}%", q(yield_pct, "%"),
                basis="isolated yield (Scheme footnote)", yield_type="isolated yield", role="main product")
    set_classification(r, pages)
    return r


scope_a = [
    ("3aa", "2a", 4, 63.0, 95, 1, None), ("3ab", "2b", 4, 63.6, 91, 1, None),
    ("3ac", "2c", 4, 63.0, 89, 1, None), ("3ad", "2d", 4, 63.7, 85, 1, None),
    ("3ae", "2e", 4, 62.8, 80, 1, None), ("3af", "2f", 4, 74.2, 89, 1, None),
    ("3ag", "2g", 5, 39.9, 60, 2, "2 equiv TFA instead of 1 equiv"),
    ("3ah", "2h", 5, 70.6, 83, 1, None), ("3ai", "2i", 5, 32.0, 63, 1, None),
    ("3aj", "2j", 5, 35.3, 61, 1, None), ("3ak", "2k", 5, 49.3, 90, 1, None),
    ("3al", "2l", 5, 54.3, 80, 1, "silica gel chromatography, AcOEt/n-hexane 20:80"),
    ("3am", "2m", 5, 46.9, 73, 1, None), ("3an", "2n", 5, 48.2, 70, 1, None),
]
for product, sulfinate, page, mass, yld, tfaeq, extra in scope_a:
    r = general_scope_record(product, "1a", sulfinate, product, page, mass, yld, tfaeq, extra)
    if product == "3al":
        r["participants"].extend([
            participant("p_etoac", M["EtOAc"], "purification_agent", (5,), stage="purification"),
            participant("p_hexane", M["hexane"], "purification_agent", (5,), stage="purification"),
            participant("p_silica", M["silica"], "purification_agent", (5,), stage="purification"),
        ])
        r["steps"].append(step("s3_purification", "purification", ["p_etoac", "p_hexane", "p_silica", "p_3al"],
                                    (5,), "silica gel column chromatography (AcOEt/n-hexane 20:80)"))
    if product == "3ah":
        r["review"]["unresolved"].append("Reactant/product thiophene connectivity was not forced into a guessed concrete SMILES.")
    reactions.append(r)

# The page-4 large-scale procedure is a separate experiment from the 0.25 mmol scope entry.
r = reaction("large_scale_3aa", (2, 3, 4), "specific_experiment", ("scale_up",),
             "Scheme 2 footnote a; scale-up discussion page 3; General Procedure for Large Scale Synthesis page 4",
             "large-scale TFA-mediated hydrosulfonylation", "allene hydrosulfonylation")
r["general_scheme_ids"] = [GENERIC_ID]
r["participants"] = [
    participant("p_1a", M["1a"], "reactant", (2, 3, 4), amount=("2 mmol", q(2, "mmol")), equiv=("1 equiv", q(1, "equiv")), reference=M["1a"]),
    participant("p_2a", M["2a"], "reactant", (2, 3, 4), amount=("4 mmol", q(4, "mmol")), equiv=("2 equiv", q(2, "equiv")), reference=M["1a"]),
    participant("p_tfa", M["TFA"], "reagent", (2, 3, 4), amount=("2 mmol", q(2, "mmol")), equiv=("1 equiv", q(1, "equiv")), reference=M["1a"]),
    participant("p_water", M["H2O"], "solvent", (2, 3, 4), volume=("4 mL", q(4, "mL"))),
    participant("p_dcm", M["DCM"], "workup_agent", (4,), stage="workup"),
    participant("p_mgso4", M["MgSO4"], "workup_agent", (4,), stage="workup"),
    participant("p_3aa", M["3aa"], "product", (2, 3, 4), product_role="main product"),
]
r["reaction_scale"] = F("1a (2 mmol)", {"reference_molecule_id": M["1a"], "amount": q(2, "mmol")}, (2, 3, 4))
set_solvent(r, (2, 3, 4), "H2O (identified by the Scheme 2 context; page-4 large-scale paragraph says only 'solvent')",
            volume=("4 mL", q(4, "mL")), concentration=("0.5 M", q(0.5, "M")), basis="1a")
r["steps"] = [
    step("s1_reaction", "reaction", ["p_1a", "p_2a", "p_tfa", "p_water", "p_3aa"], (2, 3, 4),
         "stirred under open air", ("25 °C", q(25, "°C")), ("4 h", q(4, "h")), "open air"),
    step("s2_workup", "workup", ["p_dcm", "p_mgso4", "p_3aa"], (4,),
         "aqueous workup; DCM extraction; dry over anhydrous magnesium sulfate; concentrate under vacuum"),
]
set_procedure(r, (2, 3, 4), "General Procedure for Large Scale Synthesis: sodium sulfinate, solvent, allene, and TFA; stir 25 °C for 4 h under open air; aqueous workup.")
add_outcome(r, "o_mass", "product_amount", [M["3aa"]], (4,), "514.7 mg", q(514.7, "mg"), role="main product")
add_outcome(r, "o_yield", "yield", [M["3aa"]], (2, 3, 4), "97%", q(97, "%"), basis="isolated", yield_type="isolated yield", role="main product")
set_classification(r, (2, 3, 4))
r["review"]["unresolved"] = ["The large-scale paragraph uses the generic word 'solvent'; H2O is inherited from the associated Scheme 2 context, not repeated in that paragraph."]
reactions.append(r)

# Double hydrosulfonylation product 4aa.
r = reaction("double_hydrosulfonylation_4aa", (2, 3, 6), "specific_experiment", ("derivatization",),
             "Scheme 2 footnote b; product 4aa paragraph page 6", "double hydrosulfonylation", "double allene hydrosulfonylation")
r["general_scheme_ids"] = [GENERIC_ID]
r["participants"] = [
    participant("p_1a", M["1a"], "reactant", (2, 3, 6), amount=("0.25 mmol", q(0.25, "mmol")), equiv=("1 equiv", q(1, "equiv")), reference=M["1a"]),
    participant("p_2a", M["2a"], "reactant", (2, 3, 6), amount=("1 mmol", q(1, "mmol")), equiv=("4 equiv", q(4, "equiv")), reference=M["1a"]),
    participant("p_tfa", M["TFA"], "reagent", (2, 3, 6), equiv=("4 equiv", q(4, "equiv")), reference=M["1a"]),
    participant("p_water", M["H2O"], "solvent", (2, 3, 6)),
    participant("p_thf", M["THF"], "solvent", (2, 3, 6)),
    participant("p_4aa", M["4aa"], "product", (2, 3, 6), product_role="1,3-disulfone product"),
]
r["reaction_scale"] = F("1a (0.25 mmol)", {"reference_molecule_id": M["1a"], "amount": q(0.25, "mmol")}, (2, 3, 6))
set_solvent(r, (2, 3, 6), "H2O/THF", ratio="1:1", concentration=("0.25 M", q(0.25, "M")), basis="1a")
r["steps"] = [step("s1", "reaction", [p["participant_id"] for p in r["participants"]], (2, 3, 6),
                   "reaction performed under open air", ("40 °C", q(40, "°C")), ("2 days", q(2, "d")), "open air")]
set_procedure(r, (2, 3, 6), "Modified procedure: 40 °C over two days under open air with 4 equiv TFA, 4 equiv PhSO2Na, and H2O/THF 1:1 (0.25 M).")
add_outcome(r, "o_mass", "product_amount", [M["4aa"]], (6,), "83.54 mg", q(83.54, "mg"), role="1,3-disulfone product")
add_outcome(r, "o_yield", "yield", [M["4aa"]], (2, 6), "82%", q(82, "%"), basis="isolated", yield_type="isolated yield", role="1,3-disulfone product")
set_classification(r, (2, 3, 6), cs=2, ch=2)
r["review"]["unresolved"] = ["The separate statement that 4aa can also be obtained from 3aa lacks a detailed independent procedure/yield and is retained in source_context, not counted as another experiment."]
reactions.append(r)

# Scheme 3 products following the page-4 general procedure, including explicit modifications.
for product, substrate, sulfinate, page, mass, yld in [
    ("3ba", "1b", "2a", 5, 60.1, 86), ("3bb", "1b", "2b", 5, 56.3, 80),
    ("3bi", "1b", "2i", 5, 26.2, 51), ("3bj", "1b", "2j", 6, 33.2, 57),
    ("3bm", "1b", "2m", 6, 40.8, 63), ("3da", "1d", "2a", 6, 64.1, 60),
    ("3db", "1d", "2b", 6, 38.6, 35),
]:
    r = general_scope_record(product, substrate, sulfinate, product, page, mass, yld)
    if product == "3bi":
        r["review"]["unresolved"].append("Source name/formula/drawing conflict for 3bi is preserved in the molecule identity note.")
    reactions.append(r)

# Symmetrical diallene 1c has a source conflict between the paragraph charge and modified-equivalent statement.
r = reaction("scope_3ca", (3, 4, 6), "specific_experiment", ("substrate_scope",),
             "Scheme 3 product 3ca; experimental paragraph page 6", "double hydrosulfonylation of N,N-diallenyl-2-imidazolidinone",
             "double allene hydrosulfonylation")
r["general_scheme_ids"] = [GENERIC_ID]
eq_alt1 = F("2 equiv", q(2, "equiv"), (6,))
eq_alt2 = F("4 equiv PhSO2Na", q(4, "equiv"), (3, 6))
eq_conflict = F(["2 equiv", "4 equiv"], None, (3, 6), status="conflict", alternatives=(eq_alt1, eq_alt2))
r["participants"] = [
    participant("p_1c", M["1c"], "reactant", (3, 4, 6), amount=("0.25 mmol", q(0.25, "mmol")), equiv=("1 equiv", q(1, "equiv")), reference=M["1c"]),
    participant("p_2a", M["2a"], "reactant", (3, 6), amount=("0.5 mmol in the opening sentence", q(0.5, "mmol")), equiv=eq_conflict, reference=M["1c"]),
    participant("p_tfa", M["TFA"], "reagent", (3, 6), equiv=("2 equiv", q(2, "equiv")), reference=M["1c"]),
    participant("p_water", M["H2O"], "solvent", (3, 6)),
    participant("p_thf", M["THF"], "solvent", (3, 6)),
    participant("p_dcm", M["DCM"], "workup_agent", (4,), stage="workup"),
    participant("p_mgso4", M["MgSO4"], "workup_agent", (4,), stage="workup"),
    participant("p_3ca", M["3ca"], "product", (3, 6), product_role="double hydrosulfonylation product"),
]
r["reaction_scale"] = F("1c (0.25 mmol)", {"reference_molecule_id": M["1c"], "amount": q(0.25, "mmol")}, (3, 6))
set_solvent(r, (3, 6), "H2O/THF", ratio="1:1", concentration=("0.5 M", q(0.5, "M")), basis="1c")
r["steps"] = [
    step("s1_reaction", "reaction", ["p_1c", "p_2a", "p_tfa", "p_water", "p_thf", "p_3ca"], (3, 6),
         "modified general procedure", ("25 °C", q(25, "°C")), ("4 h", q(4, "h")), "open air"),
    step("s2_workup", "workup", ["p_dcm", "p_mgso4", "p_3ca"], (4,),
         "aqueous workup, DCM extraction, MgSO4 drying, concentration"),
]
set_procedure(r, (3, 4, 6), "The 3ca paragraph states 2 equiv (0.5 mmol) 2a, then says the modified procedure used 4 equiv PhSO2Na, 2 equiv TFA, and H2O/THF 1:1.",
              general="General Procedure on page 4", modifications="2 equiv TFA; 4 equiv PhSO2Na; H2O/THF 1:1")
add_outcome(r, "o_mass", "product_amount", [M["3ca"]], (6,), "70 mg", q(70, "mg"), role="double hydrosulfonylation product")
add_outcome(r, "o_yield", "yield", [M["3ca"]], (3, 6), "60%", q(60, "%"), basis="isolated yield (Scheme footnote)", yield_type="isolated yield", role="double hydrosulfonylation product")
set_classification(r, (3, 6), cs=2, ch=2)
r["notes"] = [F("Source conflict: 3ca paragraph first gives 2a as 2 equiv/0.5 mmol, but its modified-procedure sentence and Scheme 3 footnote specify 4 equiv PhSO2Na.",
                  "Preserve both 2 equiv/0.5 mmol and 4 equiv; no preferred charge selected.", (3, 6))]
r["review"]["unresolved"] = ["2a loading conflict: 2 equiv/0.5 mmol versus 4 equiv; exact 4-equiv mmol charge is not printed."]
reactions.append(r)

# Mechanistic control experiments. They are independent from the scope/optimization rows.
def control_record(rid, pages, acid_key, solvent_key, product_keys, outcome_kind, outcome_raw, outcome_norm,
                   basis=None, note=None):
    r = reaction(rid, pages, "specific_experiment", ("control", "mechanistic_study"),
                 f"Scheme 4B control: {rid}", "control experiment for proposed hydrosulfonylation mechanism",
                 "allene hydrosulfonylation control")
    r["participants"] = [participant("p_1a", M["1a"], "reactant", pages)]
    partner_key = "PhSO2H" if rid == "control_phso2h_no_product" else "2a"
    r["participants"].append(participant("p_partner", M[partner_key], "reactant", pages))
    if acid_key:
        r["participants"].append(participant("p_acid", M[acid_key], "reagent", pages))
    r["participants"].append(participant("p_solvent", M[solvent_key], "solvent", pages))
    for i, product in enumerate(product_keys, 1):
        r["participants"].append(participant(f"p_product{i}", M[product], "product", pages, product_role="observed control product"))
    r["steps"] = [step("s1", "reaction", [p["participant_id"] for p in r["participants"]], pages,
                       "control conditions shown in Scheme 4B; unprinted time/temperature remain missing")]
    set_solvent(r, pages, solvent_key)
    set_procedure(r, pages, "Scheme 4B arrow and accompanying mechanistic-control prose.")
    targets = [M[x] for x in product_keys]
    add_outcome(r, "o1", outcome_kind, targets, pages, outcome_raw, outcome_norm, basis=basis,
                scope=("H/D incorporation at the central allene-derived position" if "D" in outcome_raw else None))
    set_classification(r, pages, applicable=bool(product_keys))
    r["review"]["unresolved"] = ["Control-experiment time, temperature, scale, workup, and purification are unreported in pages 1-6."]
    if basis is None and outcome_kind == "yield":
        r["review"]["unresolved"].append("The 100% control outcome has no explicitly stated measurement basis.")
    if note:
        r["notes"].append(F(note, note, pages))
    return r

reactions.append(control_record("control_h2o_standard", (3,), "TFA", "H2O", ["3aa"], "yield", "100%", q(100, "%"),
                                note="This Scheme 4B control is retained separately because its 100% outcome differs from the 95% isolated scope experiment and >99% NMR optimization row."))
reactions.append(control_record("control_d2o_tfa", (3, 6), "TFA", "D2O", ["d3aa", "3aa"], "selectivity",
                                "90:10 mixture of d-3aa and 3aa", {"ratio": [90, 10], "order": [M["d3aa"], M["3aa"]]},
                                basis="1H NMR", note="Prose states D2O with either TFA or deuterated TFA gave 90% deuterium incorporation."))
reactions.append(control_record("control_d2o_tfa_d", (3,), "TFA-d", "D2O", ["d3aa", "3aa"], "selectivity",
                                "H/D = 10:90", {"ratio": [10, 90], "order": ["H", "D"]},
                                basis="Scheme 4B label; 1H NMR attribution from the d-3aa paragraph",
                                note="The drawing specifies deuterated TFA/D2O; the exact scale and time are not printed."))
r = control_record("control_phso2h_no_product", (3,), None, "H2O", [], "no_reaction", "not observed", "not observed",
                   note="Benzenesulfinic acid was used instead of sodium benzenesulfinate; no product was formed.")
set_classification(r, (3,), applicable=False)
reactions.append(r)

# Proposed mechanism is kept distinct from mechanistic experiments.
r = reaction("scheme4_proposed_mechanism", (3,), "proposed_mechanism", ("other",),
             "Scheme 4A", "Proposed mechanism", "proposed hydrocarboxylation/iminium formation and sulfinate addition")
r["participants"] = [
    participant("p_1a", M["1a"], "reactant", (3,)), participant("p_tfa", M["TFA"], "reagent", (3,)),
    participant("p_2", M["gen_sulfinate"], "reactant", (3,)), participant("p_3", M["gen_product"], "product", (3,), product_role="proposed net product"),
]
r["steps"] = [step("s1_proposed_hydrocarboxylation", "reaction", ["p_1a", "p_tfa"], (3,),
                   "proposed protonation/hydrocarboxylation and equilibrium with iminium/carbocation forms"),
              step("s2_proposed_sulfinate_addition", "reaction", ["p_2", "p_3"], (3,),
                   "proposed nucleophilic substitution or conjugated addition")]
set_procedure(r, (3,), "Author-proposed mechanism in Scheme 4A; not an experimental procedure.")
r["classification"]["mapped_reaction_smiles"] = F(None, None, (3,))
r["classification"]["bond_changes"] = []
r["classification"]["labels"] = F(None, None, (3,))
r["classification"]["coverage"] = inferred("unresolved", (3,), "Mechanistic alternatives are proposed and not atom-mapped.")
r["classification"]["limitations"] = ["Proposed alternatives are not experimental proof and are not given definitive bond-change labels."]
r["review"]["unresolved"] = ["Nucleophilic substitution versus conjugated addition is explicitly unresolved by the authors."]
reactions.append(r)

# Source-context items are observations or characterization, not extra syntheses.
source_context = [
    {"context_id": f"{DOC}:context:page1_graphical_abstract", "document_id": DOC,
     "kind": "duplicate_summary_graphic", "evidence_ids": ev(1),
     "note": "Graphical abstract duplicates the current generic hydrosulfonylation and is not counted as another experiment."},
    {"context_id": f"{DOC}:context:intermediates_A_B", "document_id": DOC,
     "kind": "mechanistic_observation", "evidence_ids": ev(3, 6),
     "note": "Intermediate A was observed by 1H/19F NMR and intermediate B by HRMS; their characterization does not create synthetic experiments."},
    {"context_id": f"{DOC}:context:3aa_to_4aa_statement", "document_id": DOC,
     "kind": "underdetailed_reaction_statement", "evidence_ids": ev(3),
     "note": "The text says 4aa can also be obtained by hydrosulfonylation of 3aa, but pages 1-6 give no separate conditions or yield for that route."},
    {"context_id": f"{DOC}:context:characterization_pages4_6", "document_id": DOC,
     "kind": "characterization", "evidence_ids": ev(4, 5, 6),
     "note": "NMR, HRMS, melting points, formulas, appearances, and instrument settings are characterization context, not independent synthesis records."},
    {"context_id": f"{DOC}:context:nmr_internal_standard", "document_id": DOC,
     "kind": "measurement_context", "evidence_ids": ev(2, 4),
     "note": "Trichloroethylene was added after concentration as an NMR internal standard; it is not represented as a reaction-time reagent."},
    {"context_id": f"{DOC}:context:si_announcement", "document_id": DOC,
     "kind": "supplement_announcement", "evidence_ids": ev(4, 6),
     "note": "The article announces Supporting Information, but no SI file is supplied or explicitly mapped in this workspace; no association is inferred."},
]

# Every generated record remains partial/needs_review and gets an explicit unresolved note when none arose.
for r in reactions:
    r["completeness"] = "partial"
    r["review"]["status"] = "needs_review"
    r["review"]["method"] = "model"
    r["review"]["reviewer"] = None
    r["review"]["checked_fields"] = []
    if not r["review"]["unresolved"]:
        r["review"]["unresolved"] = ["Machine-extracted record requires expert source review; no SI was associated."]

dataset = json.loads((RUN / "dataset.json").read_text())
bundle = st.empty(dataset["documents"])
bundle["molecules"] = molecules
bundle["reactions"] = reactions
bundle["source_context"] = source_context

# Bounded mention scan for cross-region checks; literal hits are only review leads.
mentions = st.source_mentions(dataset, ["1a", "2a", "3aa", "3ca", "4aa", "d-3aa"], context=180)

# Exact per-page coverage, based on the six page-image inspections.
ids = lambda prefix: [r["reaction_id"] for r in reactions if r["reaction_id"].split(":", 1)[1].startswith(prefix)]
table1_ids = ids("table1_")
table2_ids = ids("table2_")
scope_a_ids = [f"{DOC}:scope_{x[0]}" for x in scope_a]
scheme3_ids = [f"{DOC}:scope_{x}" for x in ("3ba", "3bb", "3bi", "3bj", "3bm", "3ca", "3da", "3db")]
lit_ids = [f"{DOC}:{row[0]}" for row in lit_rows]
control_ids = [f"{DOC}:{x}" for x in ("control_h2o_standard", "control_d2o_tfa", "control_d2o_tfa_d", "control_phso2h_no_product")]
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title/abstract/graphical abstract", "finding": "Title restored; graphical abstract is a duplicate summary of the current generic reaction.", "evidence_ids": ev(1)},
            {"region": "Scheme 1 equations 1-5", "finding": "Five cited prior-work arrows recorded as literature_example records.", "evidence_ids": ev(1)},
            {"region": "Scheme 1 equation 6", "finding": "Current-paper generic TFA-mediated hydrosulfonylation; shared with repeated generic arrows on pages 2-3.", "evidence_ids": ev(1)},
        ], "record_ids": lit_ids + [GENERIC_ID], "unresolved": ["R-group structures in Scheme 1 remain generic; no yields are printed for prior examples."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Table 1 rows 1-9 and footnote", "finding": "Nine independent optimization experiments, including zero-reaction row 1, transcribed with NMR yield and conversion bases.", "evidence_ids": ev(2)},
            {"region": "Table 2 rows 1-5 and footnote", "finding": "Five independent parameter-study experiments; entry 4 retains 66% conversion and 0% product yield.", "evidence_ids": ev(2)},
            {"region": "Scheme 2 and footnotes a-d", "finding": "Fourteen 1a scope products, separate 2 mmol 3aa experiment, and modified double hydrosulfonylation 4aa linked to detailed paragraphs.", "evidence_ids": ev(2)},
        ], "record_ids": table1_ids + table2_ids + [GENERIC_ID] + scope_a_ids + [f"{DOC}:large_scale_3aa", f"{DOC}:double_hydrosulfonylation_4aa"],
         "unresolved": ["Scheme drawings were used for candidate structures; source verification remains pending."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Scheme 2 discussion", "finding": "Cross-page prose supports scope yields, scale-up 3aa, and 4aa; under-detailed 3aa-to-4aa route retained only in source_context.", "evidence_ids": ev(3)},
            {"region": "Scheme 3 and footnotes", "finding": "Eight independent products/experiments linked to pages 5-6; 3ca special conditions retained with charge conflict.", "evidence_ids": ev(3)},
            {"region": "Scheme 4A", "finding": "Author-proposed mechanism recorded separately from experiments; intermediate observations retained as source_context.", "evidence_ids": ev(3)},
            {"region": "Scheme 4B and prose", "finding": "H2O positive control, TFA/D2O and TFA-d/D2O isotope controls, and PhSO2H no-product control recorded independently.", "evidence_ids": ev(3)},
        ], "record_ids": [GENERIC_ID] + scope_a_ids + [f"{DOC}:large_scale_3aa", f"{DOC}:double_hydrosulfonylation_4aa"] + scheme3_ids + control_ids + [f"{DOC}:scheme4_proposed_mechanism"],
         "unresolved": ["Scheme 4 does not report control scale, time, temperature, workup, or purification; the proposed second step has two alternatives."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "General Considerations and General Procedure", "finding": "General addition order, 25 °C/4 h/open-air conditions, aqueous workup, DCM extraction, MgSO4 drying, and conditional purification linked to applicable scope records.", "evidence_ids": ev(4)},
            {"region": "General Procedure for Large Scale Synthesis", "finding": "2 mmol scale experiment linked to the separate large_scale_3aa record.", "evidence_ids": ev(4)},
            {"region": "3aa-3af preparation paragraphs", "finding": "Product masses and yields linked to six scope records; 3aa large-scale mass/yield linked separately.", "evidence_ids": ev(4)},
            {"region": "spectroscopy/instrument and product characterization blocks", "finding": "Characterization context only; no standalone spectra or mass peaks counted as synthesis.", "evidence_ids": ev(4)},
        ], "record_ids": scope_a_ids + scheme3_ids + [f"{DOC}:large_scale_3aa"],
         "unresolved": ["Conditional triethylamine-treated silica purification is not assigned to products unless explicitly stated; SI is not provided."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "3ag-3an preparation paragraphs", "finding": "Eight product-specific masses/yields linked to Scheme 2 scope records; 3ag 2-equiv TFA and 3al chromatography retained.", "evidence_ids": ev(5)},
            {"region": "3ba, 3bb, and 3bi preparation paragraphs", "finding": "Three Scheme 3 experiments linked; 3bi heading versus drawing/formula/reactant conflict preserved.", "evidence_ids": ev(5)},
            {"region": "NMR/HRMS/melting-point blocks", "finding": "Characterization only, not additional synthesis records.", "evidence_ids": ev(5)},
        ], "record_ids": [f"{DOC}:scope_{x}" for x in ("3ag", "3ah", "3ai", "3aj", "3ak", "3al", "3am", "3an", "3ba", "3bb", "3bi")],
         "unresolved": ["3bi source name conflicts with its formula, reactant, and drawn methylsulfonyl structure; candidate graph remains needs_review."]},
        {"page": 6, "state": "attempted", "regions": [
            {"region": "3bj, 3bm, 3ca, 3da, and 3db preparation paragraphs", "finding": "Five Scheme 3 products linked; 3ca 2-equiv/0.5 mmol versus 4-equiv PhSO2Na conflict retained.", "evidence_ids": ev(6)},
            {"region": "4aa modified procedure", "finding": "40 °C/two-day double hydrosulfonylation, mass, isolated yield, and racemic R/S assignment linked to one record.", "evidence_ids": ev(6)},
            {"region": "Intermediate A/B and d-3aa characterization", "finding": "Intermediates are source_context; d-3aa paragraph supports the isotope-control record rather than a duplicate synthesis record.", "evidence_ids": ev(6)},
            {"region": "Associated Content", "finding": "SI is announced but no SI file is supplied or mapped; association remains unestablished.", "evidence_ids": ev(6)},
        ], "record_ids": [f"{DOC}:scope_{x}" for x in ("3bj", "3bm", "3ca", "3da", "3db")] + [f"{DOC}:double_hydrosulfonylation_4aa", f"{DOC}:control_d2o_tfa"],
         "unresolved": ["No SI details are imported; isotope-control scale/workup is missing; 3ca loading conflict remains unresolved."]},
    ],
}

# Workflow is intentionally replaced because the pre-existing file described a different document/session.
workflow = f"""# Extraction workflow

## Batch and bounded source

- Batch: `allene-stereo-v31-200-20261017`
- Document: `{DOC}`
- Source: *Easy Access to Allylic Sulfones Through Transition-Metal-Free Hydrosulfonylation of Allenes*
- DOI: `10.1021/acs.joc.1c01345`
- This session claims visual coverage of exactly original PDF pages 1–6. All six supplied JPG renderings were opened at original detail.
- `runs/current/source_mapping.json` marks this as a main-text PDF with `presentation_conversion=false`; it is not a PPTX-derived presentation.
- SI association remains unestablished. The article announces SI, but no SI file is supplied or explicitly mapped, so none was inferred.
- The prepared dataset contained zero candidates and zero experiments. The prior `runs/workflow.md` described another document and was not used as a source; this log replaces it.

## Instructions and method

- Read `AGENTS.md`, the complete extraction skill, `extraction_standard.md`, `stereochemistry.md`, `structured_fields.md`, `quality_rules.md`, `literature.schema.json`, and `bond_taxonomy.json` before building records.
- Restored the article title in prepared document metadata without altering source evidence or hashes.
- Used `standard init` for the v3.1 submission template and generated records with the deterministic builder `runs/current/raw_model/build_pages_0001_0006.py`.
- Candidate concrete SMILES were checked with the bundled standard stereo diagnostic. Diagnostics are stored in `runs/current/raw_model/stereo_diagnostics_pages_0001_0006.json`; parsing is not treated as source verification.

## Page decisions

### Page 1

- Scheme 1 equations 1–5 are five `literature_example` records, not current-paper experiments. Generic R groups remain templates and the absent yields remain missing.
- Equation 6 and the graphical abstract describe the current generic reaction. Repeated generic arrows on pages 2–3 are linked to the same `generic_scheme` record rather than counted again.

### Page 2

- Table 1 rows 1–9 are nine independent optimization records. Entry 1 is a true zero-conversion/zero-yield experiment; no product-forming classification is invented.
- Table 2 rows 1–5 are five independent records. Entry 4 preserves 66% conversion of 1a and 0% yield of 3aa as distinct outcomes.
- Scheme 2 contributes fourteen 0.25 mmol scope experiments (3aa–3an), one separate 2 mmol 3aa scale-up, and one modified 4aa double hydrosulfonylation.
- Table yields/conversions retain their 1H-NMR/internal-standard bases; Scheme 2 yields retain the footnote's isolated-yield basis.

### Page 3

- Scheme 3 contributes eight experiments (3ba, 3bb, 3bi, 3bj, 3bm, 3ca, 3da, 3db). They are linked to page-4 general conditions and product paragraphs on pages 5–6.
- Source conflict for 3ca is explicit: the paragraph first states 2a = 2 equiv/0.5 mmol, but the modified-procedure sentence and Scheme 3 footnote state 4 equiv PhSO2Na. Both are retained; no preferred loading or inferred 1.0 mmol charge was selected.
- Scheme 4A is one `proposed_mechanism`, not proof. Intermediate A (NMR) and B (HRMS) remain `source_context`.
- Scheme 4B/prose yields four mechanistic controls: an H2O positive control, TFA/D2O and TFA-d/D2O isotope controls, and a PhSO2H no-product control. Missing scales/times/temperatures are preserved.

### Pages 4–6

- The page-4 general procedure supplies addition order, 25 °C, 4 h, open air, aqueous workup, DCM extraction, MgSO4 drying, and conditional purification. Conditional chromatography is not assigned unless product-specific text says so.
- Product-specific masses/yields were linked without duplicating Scheme entries. Large-scale 3aa and 4aa remain separate experiments.
- Trichloroethylene is post-reaction NMR measurement context, not a reaction-time reagent.
- Page 5 has an internal source conflict for 3bi: the heading says isopropylsulfonyl, while reactant 2i, formula C7H11NO4S, and the Scheme 3 drawing indicate methylsulfonyl. The candidate graph follows the drawing/formula but remains `needs_review`; the conflict is retained in the molecule note.
- Page 6 product 4aa retains the source `(R,S)` assignment as a racemic R/S identity with an unreported ratio; a single stereotagged isomer was not invented.
- NMR/HRMS/melting points and standalone intermediate characterizations are `source_context`, not extra syntheses.

## Stereochemistry and structures

- Every molecule has explicit `chirality_status`, `geometry_status`, and `representation_status` facts.
- N-allenyl substrates are terminal allenes; the candidate graphs indicate a non-stereogenic terminal allene, but this is recorded as a structure derivation rather than source wording.
- Products 3aa–3an, 3ba, 3bb, 3bi, 3bj, 3bm, 3ca, 3da, and 3db retain source E geometry. The E directions are encoded in candidate isomeric SMILES where a concrete graph was assigned.
- 4aa remains an R/S mixture/racemate without a guessed single-isomer SMILES. d-3aa retains E and central-deuterium wording but no forced isotope/isomeric SMILES.
- Generic R-group and unclear thiophene structures remain templates or unresolved. CAS numbers were not guessed. No different configurations were merged.

## Classification

- Positive hydrosulfonylation records carry partial transformation-level C–S and C–H formation labels; double addition records carry two of each. All are inferred and lack atom-mapped reaction SMILES.
- Zero-yield/no-product controls have `not_applicable` classification and no intended-product labels.
- Unchanged product bonds and allene bond-order decreases are not mislabeled as new sigma bonds. No three-model votes were fabricated.

## Inventory and files

- Formal records: {len(reactions)} total — {sum(r['record_type']=='specific_experiment' for r in reactions)} specific experiments, {sum(r['record_type']=='literature_example' for r in reactions)} literature examples, 1 generic scheme, and 1 proposed mechanism.
- Molecule registry: {len(molecules)} entries. Source-context items: {len(source_context)}.
- Raw bundle: `runs/current/raw_model/pages_0001_0006_extraction.json`
- Stereo diagnostics: `runs/current/raw_model/stereo_diagnostics_pages_0001_0006.json`
- Mention scan: `runs/current/raw_model/pages_0001_0006_mentions.json`
- Submission: `runs/current/submitted.json`
- Coverage: `runs/current/coverage/pages_0001_0006.json`
- Validation/audit outputs are persisted after the builder is run and reviewed.

## Known unresolved source limitations

- No SI is linked, so procedures that refer to SI or omit quantities remain partial.
- Generic R groups, 2h/3ah thiophene connectivity, and d-3aa isotope encoding remain unresolved rather than guessed.
- 3ca has a 2-equiv/0.5 mmol versus 4-equiv PhSO2Na conflict.
- 3bi has a heading-versus-formula/reactant/drawing identity conflict.
- Mechanistic controls omit scale, time, temperature, workup, and purification; the 100% H2O-control measurement basis is not stated.
- The statement that 3aa can also furnish 4aa lacks a separate procedure/yield and remains source context.
- All records remain `needs_review`; schema validation and RDKit parsing do not establish chemical correctness, source accuracy, completeness, or reproducibility.

## Validation and audit review

- `standard validate --bundle runs/current/submitted.json`: valid, zero errors.
- `standard merge --dataset runs/current/dataset.json --bundle runs/current/submitted.json`: 49 reactions merged; 85 molecules registered.
- `validate --dataset runs/current/dataset.json --strict`: valid, zero errors, v3 present; chemical accuracy is not asserted.
- The advisory audit reports zero `reaction_identity`, zero `outcome_target`, and zero `structure_gap` warnings.
- One `measurement_basis` warning remains for the Scheme 4B H2O positive control's 100% yield. The page gives `100%` but does not state isolated/NMR/other basis, so the field remains null/missing and the warning is intentionally unresolved rather than guessed.
- `documents_without_linked_si` contains this document because no SI was explicitly supplied or associated.
- Persisted checks: `runs/current/quality_checks.json` and `runs/current/raw_model/pages_0001_0006_audit.json`.
"""

assert len(coverage["pages"]) == 6 and [p["page"] for p in coverage["pages"]] == list(range(1, 7))
assert all(p["regions"] for p in coverage["pages"])
assert len(reactions) == 49, len(reactions)

st.write(RUN / "submitted.json", bundle)
st.write(RUN / "raw_model/pages_0001_0006_extraction.json", bundle)
st.write(RUN / "raw_model/stereo_diagnostics_pages_0001_0006.json", diagnostics)
st.write(RUN / "raw_model/pages_0001_0006_mentions.json", mentions)
st.write(RUN / "coverage/pages_0001_0006.json", coverage)
(ROOT / "runs/workflow.md").write_text(workflow)
print(json.dumps({"molecules": len(molecules), "reactions": len(reactions),
                  "source_context": len(source_context), "coverage_pages": len(coverage["pages"])}, ensure_ascii=False))
