#!/usr/bin/env python3
"""Build the bounded pages 1--6 v3.1 extraction for doc_4e4afb29c8afe53f.

This is a retained construction artifact.  It does not claim source review; the
standard merge deliberately quarantines every molecule and reaction.
"""

from __future__ import annotations

import json
import sys
from collections import OrderedDict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills" / "chemical-reaction-extraction" / "scripts"
sys.path.insert(0, str(SCRIPTS))

import standard_records as st


DOC = "doc_4e4afb29c8afe53f"
DATASET = ROOT / "runs" / "current" / "dataset.json"
RAW_DIR = ROOT / "runs" / "current" / "raw_model"
OUT = ROOT / "runs" / "current" / "submitted.json"


def tev(page: int) -> str:
    return f"{DOC}:page{page}:text"


def iev(page: int) -> str:
    return f"{DOC}:page{page}:image"


def ev(page: int) -> list[str]:
    return [iev(page), tev(page)]


def f(raw=None, normalized=None, evidence=(), method="original_extraction", derivation=None,
      status=None):
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation,
                   status=status)


def inferred(normalized, evidence, derivation, raw=None):
    return f(raw, normalized, evidence, method="structure_derivation", derivation=derivation)


def missing(evidence=(), status="unreported"):
    return f(evidence=evidence, status=status)


def q(number, unit, raw_text=None, qualifier=None):
    result = {"number": number, "unit": unit}
    if raw_text is not None:
        result["raw_text"] = raw_text
    if qualifier is not None:
        result["qualifier"] = qualifier
    return result


def mid(label: str) -> str:
    return f"{DOC}:mol_{label}"


dataset = json.loads(DATASET.read_text())
bundle = st.empty(dataset["documents"])
stereo_checks: dict[str, dict] = {}
molecules: OrderedDict[str, dict] = OrderedDict()


def molecule(label: str, name: str, evidence: list[str], smiles: str | None = None,
             chirality: str | None = None, geometry: str | None = None,
             representation_status: str | None = None, mixture=None,
             assignment=None, axial=None, note: str | None = None):
    identity = mid(label)
    if identity in molecules:
        m = molecules[identity]
        occurrence = {
            "document_id": DOC,
            "compound_label": f(label, label, evidence),
            "evidence_ids": list(evidence),
        }
        if occurrence not in m["occurrences"]:
            m["occurrences"].append(occurrence)
        return identity

    m = st.new_molecule(identity, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"] = [{
        "document_id": DOC,
        "compound_label": f(label, label, evidence),
        "evidence_ids": list(evidence),
    }]
    m["names"]["original"] = f(name, name, evidence)
    m["names"]["english"] = f(name, name, evidence)
    m["names"]["abbreviations"] = missing(evidence)
    m["names"]["aliases"] = f(label, [label], evidence)
    m["cas_numbers"] = missing(evidence)
    m["identity_note"] = f(
        note or "Source identity retained from the bounded article pages; no CAS or external lookup was used.",
        note or "Source identity retained from the bounded article pages; no CAS or external lookup was used.",
        evidence,
    )
    m["external_properties"] = []

    if smiles:
        diag = st.stereo_inspect(smiles)
        stereo_checks[identity] = diag
        m["representation_type"] = "smiles"
        m["smiles"] = inferred(
            smiles, evidence,
            "Candidate graph transcribed from the source drawing/name and checked only for RDKit parsing; parsing is not source verification.",
            raw="source drawing/name candidate",
        )
        atom_potential = any(x["type"].startswith("Atom_") for x in diag["potential_stereo"])
        bond_potential = any(x["type"].startswith("Bond_") for x in diag["potential_stereo"])
        chirality = chirality or ("unspecified" if atom_potential else "achiral")
        geometry = geometry or ("unspecified" if bond_potential else "not_applicable")
        representation_status = representation_status or (
            "partial" if chirality in {"unspecified", "partially_specified", "mixture", "racemic"}
            or geometry in {"unspecified", "partially_specified", "mixture"}
            else "not_applicable"
        )
        m["stereochemistry"]["isomeric_smiles"] = inferred(
            diag["canonical_isomeric_smiles"], evidence,
            "RDKit canonical isomeric graph for the candidate; unassigned potential stereo remains unassigned.",
            raw=smiles,
        )
        m["stereochemistry"]["identity_key"] = inferred(
            diag["stereo_graph_key"], evidence,
            "SHA-256 stereo graph key calculated by standard stereo_inspect; not a sample identity.",
            raw=diag["canonical_isomeric_smiles"],
        )
        m["stereochemistry"]["atom_centers"] = inferred(
            diag["atom_centers"], evidence,
            "Potential/assigned tetrahedral centers calculated for the candidate graph; the article supplies no assignment unless separately stated.",
            raw="candidate graph diagnostic",
        )
        m["stereochemistry"]["double_bonds"] = inferred(
            diag["double_bonds"], evidence,
            "Double-bond stereo diagnostic calculated for the candidate graph; aromatic bonds are not E/Z assignments.",
            raw="candidate graph diagnostic",
        )
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing(evidence)
        m["stereochemistry"]["atom_centers"] = missing(evidence)
        m["stereochemistry"]["double_bonds"] = missing(evidence)
        m["stereochemistry"]["isomeric_smiles"] = missing(evidence)
        m["stereochemistry"]["identity_key"] = missing(evidence)
        chirality = chirality or "unresolved"
        geometry = geometry or "unresolved"
        representation_status = representation_status or "unresolved"

    m["generic_template"] = None
    m["stereochemistry"]["axial"] = (
        f(axial, axial, evidence) if axial is not None else missing(evidence)
    )
    m["stereochemistry"]["mixture"] = (
        f(mixture, mixture, evidence) if mixture is not None else missing(evidence)
    )
    m["stereochemistry"]["assignment"] = (
        f(assignment, assignment, evidence) if assignment is not None else missing(evidence)
    )
    for key, val in (
        ("chirality_status", chirality),
        ("geometry_status", geometry),
        ("representation_status", representation_status),
    ):
        m["stereochemistry"][key] = inferred(
            val, evidence,
            "Stereo status distinguishes unassigned potential stereo from achiral/racemic claims and follows the bounded source plus candidate-graph diagnostics.",
            raw=val,
        )
    molecules[identity] = m
    return identity


def generic_molecule(label: str, name: str, fixed: str, variables: list[tuple[str, str]],
                     evidence: list[str]):
    identity = mid(label)
    m = st.new_molecule(identity, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"] = [{"document_id": DOC, "compound_label": f(label, label, evidence),
                         "evidence_ids": list(evidence)}]
    m["names"]["original"] = f(name, name, evidence)
    m["names"]["english"] = f(name, name, evidence)
    m["names"]["abbreviations"] = missing(evidence)
    m["names"]["aliases"] = f(label, [label], evidence)
    m["cas_numbers"] = missing(evidence)
    m["representation_type"] = "generic_template"
    m["smiles"] = missing(evidence)
    m["generic_template"] = {
        "fixed_structure": f(fixed, fixed, evidence),
        "variables": [{
            "label": var,
            "attachment_points": f("source drawing attachment", attachment, evidence),
            "definition": f(definition, definition, evidence),
            "allowed_values": missing(evidence),
            "constraints": missing(evidence),
        } for var, definition, attachment in variables],
        "instance_mapping": [],
        "reason": f("R-group structure in the source", "Generic source drawing; no wildcard concrete SMILES asserted.", evidence),
    }
    for key, val in (("chirality_status", "unresolved"), ("geometry_status", "unresolved"),
                     ("representation_status", "not_encoded")):
        m["stereochemistry"][key] = inferred(
            val, evidence, "Generic R-group drawing does not define one stereoisomer.", raw=val)
    for key in ("atom_centers", "double_bonds", "axial", "mixture", "assignment",
                "isomeric_smiles", "identity_key"):
        m["stereochemistry"][key] = missing(evidence)
    m["identity_note"] = f(
        "Generic template retained because R groups vary across the source scope.",
        "Generic template retained because R groups vary across the source scope.", evidence)
    m["external_properties"] = []
    molecules[identity] = m
    return identity


# Generic source drawings.
generic_allenol = generic_molecule(
    "generic_allenol_1", "generic 2,3-allenol 1",
    "HO-C(R1)(R2)-C(R3)=C=C(R4)(R5)",
    [("R1", "Ar, alkyl, alkenyl, or alkynyl", "carbinol carbon"),
     ("R2", "H or Me", "carbinol carbon"),
     ("R3", "H or Me", "proximal allene carbon"),
     ("R4", "H or Ph", "distal allene carbon"),
     ("R5", "H, Ph, or nBu", "distal allene carbon")], ev(1))
generic_phosphine = generic_molecule(
    "generic_poh_2", "generic P(O)H compound 2", "H-P(=O)(R6)2",
    [("R6", "aryl or n-butyl", "phosphorus")], ev(1))
generic_product = generic_molecule(
    "generic_ketophosphine_3", "generic gamma-ketophosphine oxide 3",
    "rearranged ketone bearing a beta C-P(=O)(R6)2 bond",
    [("R1-R6", "substituents inherited from source templates", "source-defined positions")], ev(1))


# Common reagents and solvents.  Salts/complexes remain unresolved rather than
# forcing fragile ionic SMILES.
simple = {
    "1a": ("1-phenyl-2,3-butadien-1-ol", "C=C=CC(O)c1ccccc1"),
    "2a": ("diphenylphosphine oxide", "[PH](=O)(c1ccccc1)c1ccccc1"),
    "2b": ("bis(4-methylphenyl)phosphine oxide", "[PH](=O)(c1ccc(C)cc1)c1ccc(C)cc1"),
    "2c": ("bis(4-chlorophenyl)phosphine oxide", "[PH](=O)(c1ccc(Cl)cc1)c1ccc(Cl)cc1"),
    "2d": ("bis(4-methoxyphenyl)phosphine oxide", "[PH](=O)(c1ccc(OC)cc1)c1ccc(OC)cc1"),
    "2e": ("dibutylphosphine oxide", "[PH](=O)(CCCC)CCCC"),
    "2f": ("ethyl phenylphosphinate", "[PH](=O)(OCC)c1ccccc1"),
    "2g": ("diethyl phosphonate", "[PH](=O)(OCC)OCC"),
    "TfOH": ("trifluoromethanesulfonic acid", "OS(=O)(=O)C(F)(F)F"),
    "TFA": ("trifluoroacetic acid", "O=C(O)C(F)(F)F"),
    "H2SO4": ("sulfuric acid", "OS(=O)(=O)O"),
    "THF": ("tetrahydrofuran", "C1CCOC1"),
    "DCE": ("1,2-dichloroethane", "ClCCCl"),
    "DCM": ("dichloromethane", "ClCCl"),
    "toluene": ("toluene", "Cc1ccccc1"),
    "dioxane": ("1,4-dioxane", "O1CCOCC1"),
    "EtOAc": ("ethyl acetate", "CCOC(=O)C"),
    "DMSO": ("dimethyl sulfoxide", "CS(=O)C"),
    "HBr": ("hydrobromic acid", "Br"),
    "petroleum_ether": ("petroleum ether", None),
    "silica_gel": ("silica gel 300-400 mesh", None),
}
for label, (name, smi) in simple.items():
    page = 4 if label in {"petroleum_ether", "silica_gel"} else 3 if label in {"EtOAc", "DMSO", "HBr"} else 2
    molecule(label, name, ev(page), smi)

for label, name in {
    "ZnOTf2": "zinc trifluoromethanesulfonate",
    "AlOTf3": "aluminum trifluoromethanesulfonate",
    "NHTf2": "bis(trifluoromethanesulfonyl)amine",
    "BPhF5_3": "tris(pentafluorophenyl)borane",
    "BF3Et2O": "boron trifluoride diethyl etherate",
}.items():
    molecule(label, name, ev(2), None, chirality="unresolved", geometry="not_applicable",
             representation_status="unresolved")

# The P-stereogenic starting material is explicitly identified but is not forced
# into a toolkit stereochemical encoding.
molecule("2h", "(S)-tert-butyl(phenyl)phosphine oxide", ev(3), None,
         chirality="specified", geometry="not_applicable", representation_status="not_encoded",
         assignment="(S), 99% ee (source statement)",
         note="Absolute P configuration and 99% ee are source-reported; no P-stereo SMILES was asserted.")


PDP = "P(=O)(c1ccccc1)c1ccccc1"


def ketone(aryl: str) -> str:
    return f"CC(=O)CC({aryl}){PDP}"


aryl = {
    "a": "c1ccccc1", "b": "c1ccc(Cl)cc1", "c": "c1cccc(Cl)c1",
    "d": "c1ccccc1Cl", "e": "c1ccc(Br)cc1", "f": "c1ccc(F)cc1",
    "g": "c1ccc(C)cc1", "h": "c1ccc(C(F)(F)F)cc1", "i": "c1ccc(SC)cc1",
    "j": "c1ccc(C(=O)OC)cc1", "k": "c1ccc([N+](=O)[O-])cc1",
    "l": "c1ccc(C#N)cc1", "m": "c1ccc(-c2ccccc2)cc1", "n": "c1ccc(OC)cc1",
    "o": "c1cccc2ccccc12", "p": "c1ccoc1", "q": "c1ccsc1",
    "r": "c1cccnc1", "s": "c1csc2ccccc12",
}
product_smiles = {letter: ketone(group) for letter, group in aryl.items()}
product_smiles.update({
    "t": f"CC(=O)CC(CCCCC){PDP}",
    "u": f"CC(=O)CC(C1CCCCC1){PDP}",
    "v": f"CC(=O)CC(Cc1ccccc1){PDP}",
    "w": f"CC(=O)CC(C=Cc1ccccc1){PDP}",
    "x": f"CC(=O)CC(C#Cc1ccccc1){PDP}",
    "y": f"CC(=O)CC(C)(CCCCCC){PDP}",
    "z": f"CC(=O)CC(C)(c1ccccc1){PDP}",
    "aa": f"CC(=O)CC1({PDP})CCCCC1",
    "ac": f"O=C(CC(c1ccccc1){PDP})C(c1ccccc1)c1ccccc1",
    "ad": f"O=C(CC(c1ccccc1){PDP})Cc1ccccc1",
    "ae": f"CCCCC(=O)CC(CCCCC){PDP}",
    "af": f"CCCCC(=O)CC(c1ccccc1){PDP}",
    "ag": f"CCCCC(=O)CC(c1ccc(Cl)cc1){PDP}",
    "ah": f"CCCCC(=O)CC(c1ccc(C)cc1){PDP}",
    "ai": f"CCCCC(=O)CC(c1ccc(C(F)(F)F)cc1){PDP}",
})

reactant_smiles = {letter: f"C=C=CC(O){group}" for letter, group in aryl.items()}
reactant_smiles.update({
    "t": "C=C=CC(O)CCCCC", "u": "C=C=CC(O)C1CCCCC1",
    "v": "C=C=CC(O)Cc1ccccc1", "w": "C=C=CC(O)C=Cc1ccccc1",
    "x": "C=C=CC(O)C#Cc1ccccc1", "y": "C=C=CC(O)(C)CCCCCC",
    "z": "C=C=CC(O)(C)c1ccccc1", "aa": "C=C=CC1(O)CCCCC1",
    "ab": "C=C=C(C)C(O)c1ccccc1",
})

product_names = {
    "a": "4-(diphenylphosphinyl)-4-phenylbutan-2-one",
    "b": "4-(diphenylphosphinyl)-4-(4-chlorophenyl)butan-2-one",
    "c": "4-(diphenylphosphinyl)-4-(3-chlorophenyl)butan-2-one",
    "d": "4-(diphenylphosphinyl)-4-(2-chlorophenyl)butan-2-one",
    "e": "4-(diphenylphosphinyl)-4-(4-bromophenyl)butan-2-one",
    "f": "4-(diphenylphosphinyl)-4-(4-fluorophenyl)butan-2-one",
    "g": "4-(diphenylphosphinyl)-4-(4-tolyl)butan-2-one",
    "h": "4-(diphenylphosphinyl)-4-(4-(trifluoromethyl)phenyl)butan-2-one",
    "i": "4-(diphenylphosphinyl)-4-(4-methylthiophenyl)butan-2-one",
    "j": "methyl 4-(1-(diphenylphosphinyl)-3-oxo-3-methylpropyl)benzoate",
    "k": "4-(diphenylphosphinyl)-4-(4-nitrophenyl)butan-2-one",
    "l": "4-(1-(diphenylphosphinyl)-3-oxo-3-methylpropyl)benzonitrile",
    "m": "4-(diphenylphosphinyl)-4-([1,1'-biphenyl]-4-yl)butan-2-one",
    "n": "4-(diphenylphosphinyl)-4-(4-methoxyphenyl)butan-2-one",
    "o": "4-(diphenylphosphinyl)-4-(naphthalen-1-yl)butan-2-one",
    "p": "4-(diphenylphosphinyl)-4-(furan-2-yl)butan-2-one",
    "q": "4-(diphenylphosphinyl)-4-(thiophen-2-yl)butan-2-one",
    "r": "4-(diphenylphosphinyl)-4-(pyridin-3-yl)butan-2-one",
    "s": "4-(diphenylphosphinyl)-4-(benzothiophen-3-yl)butan-2-one",
    "t": "4-(diphenylphosphinyl)nonan-2-one",
    "u": "4-(diphenylphosphinyl)-4-cyclohexylbutan-2-one",
    "v": "4-(diphenylphosphinyl)-5-phenylpentan-2-one",
    "w": "4-(diphenylphosphinyl)-6-phenylhex-5-en-2-one",
    "x": "4-(diphenylphosphinyl)-6-phenylhex-5-yn-2-one",
    "y": "4-(diphenylphosphinyl)-4-methyldecan-2-one",
    "z": "4-methyl-4-(diphenylphosphinyl)-4-phenylbutan-2-one",
    "aa": "1-(1-diphenylphosphinylcyclohexyl)propan-2-one",
    "ab": "3-methyl-4-(diphenylphosphinyl)-4-phenylbutan-2-one (diastereomer mixture)",
    "ac": "4-(diphenylphosphinyl)-1,1,4-triphenylbutan-2-one",
    "ad": "4-(diphenylphosphinyl)-1,4-diphenylbutan-2-one",
    "ae": "7-(diphenylphosphinyl)dodecan-5-one",
    "af": "1-(diphenylphosphinyl)-1-phenylheptan-3-one",
    "ag": "1-(diphenylphosphinyl)-1-(4-chlorophenyl)heptan-3-one",
    "ah": "1-(diphenylphosphinyl)-1-(4-methylphenyl)heptan-3-one",
    "ai": "1-(diphenylphosphinyl)-1-(4-(trifluoromethyl)phenyl)heptan-3-one",
}

letters = list(product_names)
for letter in letters:
    p_ev = ev(2) + (ev(3) if letter in {"r", "ab", "ac", "ad", "ae", "af", "ag", "ah", "ai"} else [])
    r_smi = reactant_smiles.get(letter)
    molecule(f"1{letter}", f"2,3-allenol 1{letter}", p_ev, r_smi,
             chirality="unspecified" if r_smi else "unresolved",
             geometry="unspecified" if letter in {"w", "ad", "ae", "af", "ag", "ah", "ai"} else "not_applicable",
             representation_status="partial" if r_smi else "unresolved",
             axial="allene axial configuration not assigned by the source" if letter in {"ad", "ae", "af", "ag", "ah", "ai"} else None,
             note="Scope substrate identity follows Table 2. Candidate graph is omitted where the compact drawing is not sufficiently certain.")
    if letter == "ab":
        molecule(f"3{letter}", product_names[letter], p_ev, None, chirality="mixture",
                 geometry="not_applicable", representation_status="not_encoded",
                 mixture="diastereomeric mixture, dr 65:35; configurations unassigned")
    else:
        molecule(f"3{letter}", product_names[letter], p_ev,
                 product_smiles.get(letter), chirality="unspecified",
                 geometry="unspecified" if letter == "w" else "not_applicable",
                 representation_status="partial",
                 note="Product identity follows Table 2 and, where present on pages 4-6, the experimental heading; absolute configuration is not assigned.")


# Table 3 products and special product-mixture identities.
table3_products = {
    "aj": ("4-(bis(4-methylphenyl)phosphinyl)-4-phenylbutan-2-one",
           "CC(=O)CC(c1ccccc1)P(=O)(c1ccc(C)cc1)c1ccc(C)cc1"),
    "ak": ("4-(bis(4-chlorophenyl)phosphinyl)-4-phenylbutan-2-one",
           "CC(=O)CC(c1ccccc1)P(=O)(c1ccc(Cl)cc1)c1ccc(Cl)cc1"),
    "al": ("4-(bis(4-methoxyphenyl)phosphinyl)-4-phenylbutan-2-one",
           "CC(=O)CC(c1ccccc1)P(=O)(c1ccc(OC)cc1)c1ccc(OC)cc1"),
    "am": ("4-(di-n-butylphosphinyl)-4-phenylbutan-2-one",
           "CC(=O)CC(c1ccccc1)P(=O)(CCCC)CCCC"),
    "ao": ("diethyl (3-oxo-1-phenylbutyl)phosphonate (target 3ao)",
           "CC(=O)CC(c1ccccc1)P(=O)(OCC)OCC"),
}
for label, (name, smi) in table3_products.items():
    molecule(f"3{label}", name, ev(3), smi, chirality="unspecified",
             geometry="not_applicable", representation_status="partial")
molecule("3an", "ethyl (3-oxobutyl)(phenyl)phosphinate (reported diastereomer mixture 3an)",
         ev(3), None, chirality="mixture", geometry="not_applicable",
         representation_status="not_encoded", mixture="dr 1.6:1; configurations unassigned")
molecule("3ap", "4-(tert-butyl(phenyl)phosphinyl)-4-phenylbutan-2-one (3ap mixture)",
         ev(3), None, chirality="racemic", geometry="not_applicable",
         representation_status="not_encoded",
         mixture="source reports a racemic form and dr 3.4:1; individual configurations unassigned",
         assignment="product reported as racemic form")

molecule("4", "3-bromo-4-(diphenylphosphinyl)-4-phenylbutan-2-one", ev(3),
         "CC(=O)C(Br)C(c1ccccc1)P(=O)(c1ccccc1)c1ccccc1",
         chirality="unspecified", geometry="not_applicable", representation_status="partial")
molecule("5", "(E)-4-phenyl-3-buten-2-one", ev(3), "CC(=O)/C=C/c1ccccc1",
         chirality="achiral", geometry="specified", representation_status="encoded",
         assignment="(E) assignment shown/named for enone 5")


reactions: list[dict] = []


def participant(pid: str, molecule_id: str, role: str, evidence: list[str], *,
                amount=None, mass=None, volume=None, equivalents=None, mol_percent=None,
                eq_ref=None, stage="reaction", product_role=None):
    p = st.new_participant(pid, molecule_id, role, evidence, stage=stage)
    p["product_role"] = f(product_role, product_role, evidence) if product_role else missing(evidence)
    for key, value_ in (("amount", amount), ("mass", mass), ("volume", volume),
                        ("equivalents", equivalents), ("mol_percent", mol_percent)):
        p["charges"][key] = f(value_, value_, evidence) if value_ is not None else missing(evidence)
    p["equivalence_reference"] = (
        inferred(eq_ref, evidence, "Loading/equivalent denominator explicitly or contextually refers to the named substrate.", raw=eq_ref)
        if eq_ref else missing(evidence)
    )
    return p


def step(step_id: str, evidence: list[str], pids: list[str], operation: str,
         temperature=None, time=None, atmosphere=None, stage="reaction"):
    s = st.new_step(step_id, evidence, stage=stage)
    s["participant_ids"] = list(pids)
    s["operation"] = f(operation, operation, evidence)
    s["temperature"] = f(temperature, temperature, evidence) if temperature is not None else missing(evidence)
    s["time"] = f(time, time, evidence) if time is not None else missing(evidence)
    s["atmosphere"] = f(atmosphere, atmosphere, evidence) if atmosphere is not None else missing(evidence)
    s["pressure"] = missing(evidence)
    return s


def outcome(oid: str, kind: str, targets: list[str], evidence: list[str], value,
            yield_type=None, basis=None, scope=None, product_role=None):
    o = st.new_outcome(oid, kind, targets, evidence)
    o["product_role"] = f(product_role, product_role, evidence) if product_role else missing(evidence)
    o["value"] = f(value, value, evidence)
    o["yield_type"] = f(yield_type, yield_type, evidence) if yield_type else missing(evidence)
    o["measurement_basis"] = f(basis, basis, evidence) if basis else missing(evidence)
    o["scope"] = f(scope, scope, evidence) if scope else missing(evidence)
    return o


def set_classification(r: dict, evidence: list[str], changes: list[dict] | None,
                       coverage="partial", limitations=None):
    c = r["classification"]
    c["mapped_reaction_smiles"] = missing(evidence)
    if changes is None:
        c["bond_changes"] = []
        c["labels"] = inferred([], evidence, "No observed transformation (or classification not applicable).", raw=[])
        c["coverage"] = inferred(coverage, evidence, "No bond formation is asserted for this record.", raw=coverage)
    else:
        labels: OrderedDict[str, dict] = OrderedDict()
        for i, change in enumerate(changes, 1):
            change = dict(change)
            change.setdefault("change_id", f"b{i}")
            c["bond_changes"].append(inferred(
                change, evidence,
                "Partial net bond-change model derived from the source structures/text; no complete atom mapping is claimed.",
                raw="source transformation drawing",
            ))
            if change["change_type"] in {"formed", "order_increased", "hydrogen_gain"}:
                lab = st.bond_label(change["atoms"][0]["element"], change["atoms"][1]["element"],
                                    change["after_order"])
                labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})[
                    "supporting_change_ids"].append(change["change_id"])
        c["labels"] = inferred(
            list(labels.values()), evidence,
            "Taxonomy labels calculated from all listed formed/promoted bonds under bond-formation-20261016.",
            raw="source transformation drawing",
        )
        c["coverage"] = inferred(coverage, evidence, "Partial source-derived net changes; atom correspondence and coproducts remain unresolved.", raw=coverage)
    c["limitations"] = limitations or [
        "No complete atom-mapped reaction SMILES is asserted.",
        "Oxygen/proton provenance, hydrogen transfers, leaving groups, and coproduct balance are unresolved.",
    ]
    c["model_votes"] = []


MAIN_CHANGES = [
    {"atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "P"}],
     "before_order": 0, "after_order": 1, "change_type": "formed"},
    {"atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}],
     "before_order": 1, "after_order": 2, "change_type": "order_increased"},
]


def reaction(suffix: str, evidence: list[str], record_type: str, categories: list[str],
             name: str, transformation: str, parts=None, steps=None, outcomes=None,
             notes=None, changes=MAIN_CHANGES, coverage="partial", unresolved=None,
             schemes=None, location=None):
    rid = f"{DOC}:{suffix}"
    r = st.new_reaction(rid, DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, evidence)
    r["location"] = f(location or suffix, location or suffix, evidence)
    r["general_scheme_ids"] = list(schemes or [])
    r["reaction_name"] = f(name, name, evidence)
    r["transformation_type"] = inferred(
        transformation, evidence, "Transformation description derived from the source reaction arrow/prose.", raw=name)
    r["participants"] = list(parts or [])
    r["reaction_scale"] = missing(evidence)
    r["solvent_conditions"] = {key: missing(evidence) for key in (
        "composition", "mixture_ratio", "total_volume", "concentration", "concentration_basis")}
    r["steps"] = list(steps or [])
    r["special_conditions"] = {key: missing(evidence) for key in (
        "light_source", "wavelength", "power", "current", "voltage", "electrodes",
        "flow_parameters", "reactor")}
    r["procedure"] = {key: missing(evidence) for key in (
        "original_text", "general_procedure_reference", "modifications", "addition_order", "pretreatment")}
    r["outcomes"] = list(outcomes or [])
    r["notes"] = [f(n, n, evidence) for n in (notes or [])]
    r["evidence_ids"] = list(dict.fromkeys(evidence))
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": list(unresolved or [])}
    r["completeness"] = "partial"
    set_classification(r, evidence, changes, coverage)
    reactions.append(r)
    return r


# Page 1 prior-work arrows are classified explicitly as literature examples and
# are not mixed with current-paper experiments.
for suffix, name in [
    ("literature_scheme1a_base", "base-promoted phospha-Michael addition (prior work)"),
    ("literature_scheme1a_nhc", "NHC-catalyzed phospha-Michael addition (prior work)"),
    ("literature_scheme1b_ag", "Ag-promoted radical phosphinylation of diaryl allylic alcohols (prior work)"),
    ("literature_scheme1b_eosin", "Eosin Y/blue-LED phosphinylation of allylic alcohols (prior work)"),
    ("literature_scheme1c_znotf", "Zn(OTf)2-catalyzed phosphinylation of propargylic alcohols (prior work)"),
]:
    reaction(suffix, ev(1), "literature_example", [], name, "prior-work context",
             changes=None, coverage="unresolved",
             notes=["Displayed in Scheme 1 as cited background, not an experiment newly reported by this paper."],
             unresolved=["Generic prior-work structures and conditions were not expanded into concrete molecules in this bounded current-paper extraction."],
             location="PDF page 1, Scheme 1 prior-work arrow")


current_scheme = reaction(
    "generic_current_transformation", ev(1) + ev(2) + ev(3), "generic_scheme", ["other"],
    "TfOH-catalyzed phosphinylation of 2,3-allenols", "allenol rearrangement followed by phospha-Michael addition",
    parts=[participant("p_allenol", generic_allenol, "reactant", ev(1)),
           participant("p_poh", generic_phosphine, "reactant", ev(1)),
           participant("p_tfoh", mid("TfOH"), "catalyst", ev(1), mol_percent=q(30, "mol%"), eq_ref=generic_allenol),
           participant("p_thf", mid("THF"), "solvent", ev(1)),
           participant("p_product", generic_product, "product", ev(1), product_role="target")],
    steps=[step("generic_reaction", ev(1), ["p_allenol", "p_poh", "p_tfoh", "p_thf"],
                "couple under the displayed current-work conditions", q(100, "°C"), q(24, "h"), "Ar")],
    notes=["Scheme 1d/current-work abstract graphic is a generic transformation, not a separate concrete experiment."],
    unresolved=["Generic R-group template does not define one atom mapping or stereoisomer."],
    location="PDF pages 1-3, Scheme 1d and generic arrows in Tables 2-3")
current_scheme_id = current_scheme["reaction_id"]


# Table 1: every optimization row is an independent experiment.
opt_rows = [
    (1, "ZnOTf2", 20, "DCE", 100, 81), (2, "ZnOTf2", 20, "toluene", 100, 69),
    (3, "ZnOTf2", 20, "dioxane", 100, 87), (4, "ZnOTf2", 20, "THF", 100, 90),
    (5, "AlOTf3", 20, "THF", 100, 88), (6, "NHTf2", 20, "THF", 100, 76),
    (7, "BPhF5_3", 20, "THF", 100, 35), (8, "BF3Et2O", 20, "THF", 100, 64),
    (9, "TFA", 20, "THF", 100, 24), (10, "H2SO4", 20, "THF", 100, 82),
    (11, "TfOH", 20, "THF", 100, 92), (12, None, None, "THF", 100, 0),
    (13, "TfOH", 30, "THF", 100, 95), (14, "TfOH", 40, "THF", 100, 93),
    (15, "TfOH", 10, "toluene", 100, 58), (16, "TfOH", 10, "DCM", 100, 27),
    (17, "TfOH", 10, "DCE", 100, 89), (18, "TfOH", 30, "THF", 90, 91),
    (19, "TfOH", 30, "THF", 110, 88),
]
for entry, catalyst, loading, solvent, temp, yld in opt_rows:
    e = ev(2)
    parts = [
        participant("p_1a", mid("1a"), "reactant", e, amount=q(0.4, "mmol")),
        participant("p_2a", mid("2a"), "reactant", e, amount=q(0.48, "mmol")),
        participant("p_solvent", mid(solvent), "solvent", e, volume=q(2, "mL")),
        participant("p_3a", mid("3a"), "product", e, product_role="target"),
    ]
    if catalyst:
        parts.insert(2, participant("p_catalyst", mid(catalyst), "catalyst", e,
                                    mol_percent=q(loading, "mol%"), eq_ref=mid("1a")))
    y = outcome("yield_3a", "yield", [mid("3a")], e, q(yld, "%"),
                basis=None, product_role="target")
    r = reaction(
        f"table1_entry_{entry:02d}", e, "specific_experiment", ["optimization"],
        "optimization of the phosphinylation of 1a with 2a", "TfOH/Lewis-acid-catalyzed allenol phosphinylation",
        parts=parts,
        steps=[step("reaction", e, [p["participant_id"] for p in parts if p["role"] != "product"],
                    "heat the charged reaction mixture", q(temp, "°C"), q(24, "h"), "Ar")],
        outcomes=[y], changes=None if yld == 0 else MAIN_CHANGES,
        coverage="not_applicable" if yld == 0 else "partial",
        schemes=[current_scheme_id], location=f"PDF page 2, Table 1 entry {entry}",
        unresolved=["Yield measurement basis is not reported in Table 1; it was not guessed."],
        notes=["Entry 12 explicitly omits a catalyst and reports only recovered starting materials in the adjoining prose."] if entry == 12 else [],
    )
    r["reaction_scale"] = f("1a (0.4 mmol)", {"reference_molecule_id": mid("1a"), "amount": q(0.4, "mmol")}, e)
    r["solvent_conditions"]["composition"] = f(solvent, [solvent], e)
    r["solvent_conditions"]["total_volume"] = f("2 mL", q(2, "mL"), e)
    r["procedure"]["original_text"] = f("Table 1 footnote a and row values", "Table 1 footnote a and row values", e)


scope_yields = {
    "a": 95, "b": 97, "c": 99, "d": 94, "e": 92, "f": 88, "g": 99,
    "h": 71, "i": 50, "j": 91, "k": 94, "l": 74, "m": 98, "n": 99,
    "o": 95, "p": 73, "q": 99, "r": 42, "s": 92, "t": 91, "u": 85,
    "v": 77, "w": 86, "x": 55, "y": 79, "z": 99, "aa": 89, "ab": 87,
    "ac": 56, "ad": 25, "ae": 42, "af": 90, "ag": 85, "ah": 88, "ai": 57,
}
scope_masses = {
    "a": 132, "b": 149, "c": 152, "d": 144, "e": 157, "f": 129, "g": 144,
    "h": 118, "i": 79, "j": 147, "k": 148, "l": 116, "m": 166, "n": 150,
    "o": 151, "p": 99, "q": 140, "r": 58, "s": 149, "t": 125, "u": 121,
    "v": 112, "w": 129,
}


def char_page(letter: str) -> int | None:
    if letter in {"a", "b", "c", "d"}: return 4
    if letter in {"e", "f", "g", "h", "i", "j", "k", "l", "m", "n"}: return 5
    if letter in {"o", "p", "q", "r", "s", "t", "u", "v", "w"}: return 6
    return None


for letter, yld in scope_yields.items():
    pages = [2]
    if letter in {"r", "ab", "ac", "ad", "ae", "af", "ag", "ah", "ai"}:
        pages.append(3)
    cp = char_page(letter)
    if cp:
        pages.append(cp)
    e = []
    for page in pages:
        e += ev(page)
    e = list(dict.fromkeys(e))
    tfoh_loading = 130 if letter == "r" else 30
    temp = 90 if letter == "r" else 100
    parts = [
        participant("p_allenol", mid(f"1{letter}"), "reactant", e, amount=q(0.4, "mmol")),
        participant("p_2a", mid("2a"), "reactant", e, amount=q(0.48, "mmol")),
        participant("p_tfoh", mid("TfOH"), "catalyst", e, mol_percent=q(tfoh_loading, "mol%"), eq_ref=mid(f"1{letter}")),
        participant("p_thf", mid("THF"), "solvent", e, volume=q(2, "mL")),
        participant("p_product", mid(f"3{letter}"), "product", e, product_role="target"),
    ]
    outs = [outcome("yield", "yield", [mid(f"3{letter}")], e, q(yld, "%"),
                       yield_type="isolated" if cp else None,
                       basis="isolated after flash chromatography" if cp else None,
                       product_role="target")]
    if letter in scope_masses:
        outs.append(outcome("product_amount", "product_amount", [mid(f"3{letter}")], e,
                            q(scope_masses[letter], "mg"), basis="isolated product mass", product_role="target"))
    if letter == "ab":
        outs.append(outcome("dr", "selectivity", [mid("3ab")], ev(2) + ev(3), "65:35",
                            basis="31P NMR analysis of the crude products",
                            scope="diastereomer ratio; configurations/order unassigned"))
    unresolved = ["Absolute configuration of the product's potential stereocenter(s) is not assigned."]
    if not cp:
        unresolved.append("Yield measurement basis is not stated on bounded pages 1-6.")
    r = reaction(
        f"table2_scope_3{letter}", e, "specific_experiment", ["substrate_scope"],
        "TfOH-catalyzed phosphinylation of a 2,3-allenol", "allenol rearrangement/phospha-Michael coupling",
        parts=parts,
        steps=[step("reaction", e, ["p_allenol", "p_2a", "p_tfoh", "p_thf"],
                    "stir in an oil bath", q(temp, "°C"), q(24, "h"), "Ar")],
        outcomes=outs, schemes=[current_scheme_id], unresolved=unresolved,
        location=f"PDF pages {','.join(map(str, pages))}; Table 2 product 3{letter}" if len(pages) > 1 else f"PDF page 2, Table 2 product 3{letter}",
    )
    r["reaction_scale"] = f("2,3-allenol (0.4 mmol)", {"reference_molecule_id": mid(f"1{letter}"), "amount": q(0.4, "mmol")}, e)
    r["solvent_conditions"]["composition"] = f("THF", ["THF"], e)
    r["solvent_conditions"]["total_volume"] = f("2 mL", q(2, "mL"), e)
    r["procedure"]["general_procedure_reference"] = f(
        "General Procedure for the Synthesis of gamma-Ketophosphine Oxides", "page-4 general procedure", e,
        method="inherited_general_condition")


table3 = [
    ("aj", "2b", 99, 100, "THF", None), ("ak", "2c", 95, 100, "THF", None),
    ("al", "2d", 94, 100, "THF", None), ("am", "2e", 67, 140, "THF", None),
    ("an", "2f", 12, 100, "EtOAc", "1.6:1"), ("ao", "2g", 0, 100, "THF", None),
    ("ap", "2h", 81, 100, "THF", "3.4:1"),
]
for prod, reagent, yld, temp, solvent, dr in table3:
    e = ev(3)
    parts = [
        participant("p_1a", mid("1a"), "reactant", e, amount=q(0.4, "mmol")),
        participant("p_poh", mid(reagent), "reactant", e,
                    amount=q(0.48, "mmol")),
        participant("p_tfoh", mid("TfOH"), "catalyst", e, mol_percent=q(30, "mol%"), eq_ref=mid("1a")),
        participant("p_solvent", mid(solvent), "solvent", e, volume=q(2, "mL")),
        participant("p_product", mid(f"3{prod}"), "product", e, product_role="target"),
    ]
    outs = [outcome("yield", "yield", [mid(f"3{prod}")], e, q(yld, "%"), product_role="target")]
    if dr:
        outs.append(outcome("dr", "selectivity", [mid(f"3{prod}")], e, dr,
                            basis="1H NMR analysis of the mixture" if prod == "ap" else None,
                            scope="diastereomer ratio; configurations/order unassigned"))
    r = reaction(
        f"table3_scope_3{prod}", e, "specific_experiment", ["substrate_scope"],
        "TfOH-catalyzed coupling with a P(O)H compound", "allenol rearrangement/phospha-Michael coupling",
        parts=parts,
        steps=[step("reaction", e, ["p_1a", "p_poh", "p_tfoh", "p_solvent"],
                    "heat under the Table 3 conditions", q(temp, "°C"), q(24, "h"), "Ar")],
        outcomes=outs, changes=None if yld == 0 else MAIN_CHANGES,
        coverage="not_applicable" if yld == 0 else "partial", schemes=[current_scheme_id],
        unresolved=[
            "Table 3 caption literally states 2a (0.48 mmol) while the row varies 2b-2h; 0.48 mmol is retained as an inherited row charge but needs expert confirmation.",
            "Yield measurement basis is not reported on bounded pages 1-6.",
        ] + (["The 3an dr measurement method is not stated; superscript c denotes EtOAc, not the method."] if prod == "an" else []),
        notes=["Table 3 footnote b: 140 °C."] if prod == "am" else (["Table 3 footnote c: in EtOAc."] if prod == "an" else []),
        location=f"PDF page 3, Table 3 product 3{prod}",
    )
    r["reaction_scale"] = f("1a (0.4 mmol)", {"reference_molecule_id": mid("1a"), "amount": q(0.4, "mmol")}, e)
    r["solvent_conditions"]["composition"] = f(solvent, [solvent], e)
    r["solvent_conditions"]["total_volume"] = f("2 mL", q(2, "mL"), e)


# Application and mechanistic/control experiments on page 3.
e3 = ev(3)
scale_parts = [
    participant("p_1a", mid("1a"), "reactant", e3, mass=q(1.46, "g"), amount=q(10, "mmol")),
    participant("p_2a", mid("2a"), "reactant", e3, mass=q(2.43, "g"), amount=q(12, "mmol")),
    participant("p_tfoh", mid("TfOH"), "catalyst", e3, mol_percent=q(30, "mol%"), eq_ref=mid("1a")),
    participant("p_thf", mid("THF"), "solvent", e3),
    participant("p_3a", mid("3a"), "product", e3, product_role="target"),
]
r = reaction(
    "scheme2a_gram_scale_3a", e3, "specific_experiment", ["scale_up"],
    "gram-scale phosphinylation of 1a", "allenol rearrangement/phospha-Michael coupling",
    parts=scale_parts,
    steps=[step("reaction", e3, ["p_1a", "p_2a", "p_tfoh", "p_thf"],
                "conduct the gram-scale reaction", q(100, "°C"), q(24, "h"), "Ar")],
    outcomes=[outcome("yield", "yield", [mid("3a")], e3, q(99, "%"), product_role="target"),
              outcome("product_amount", "product_amount", [mid("3a")], e3, q(3.46, "g"),
                      basis="reported product mass", product_role="target")],
    schemes=[current_scheme_id],
    unresolved=["Yield measurement/isolation basis is not explicitly stated in the bounded Scheme 2a/prose."],
    location="PDF page 3, Scheme 2a")
r["reaction_scale"] = f("1a (10 mmol)", {"reference_molecule_id": mid("1a"), "amount": q(10, "mmol")}, e3)

deriv_parts = [
    participant("p_3a", mid("3a"), "reactant", e3, amount=q(1.0, "mmol")),
    participant("p_hbr", mid("HBr"), "reagent", e3, equivalents=q(3, "equiv"), eq_ref=mid("3a")),
    participant("p_dmso", mid("DMSO"), "solvent", e3),
    participant("p_4", mid("4"), "product", e3, product_role="target"),
]
reaction(
    "scheme2b_bromination_4", e3, "specific_experiment", ["derivatization"],
    "conversion of 3a to alpha-bromocarbonyl compound 4", "alpha bromination",
    parts=deriv_parts,
    steps=[step("reaction", e3, ["p_3a", "p_hbr", "p_dmso"], "react with HBr in DMSO", q(55, "°C"), q(4, "h"))],
    outcomes=[outcome("yield", "yield", [mid("4")], e3, q(56, "%"), product_role="target")],
    changes=[{"atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "Br"}],
              "before_order": 0, "after_order": 1, "change_type": "formed"}],
    unresolved=["Yield measurement basis and workup are not stated until outside the bounded pages."],
    location="PDF page 3, Scheme 2b")

def control_1a_to_5(suffix, loading, yld):
    parts = [participant("p_1a", mid("1a"), "reactant", e3, amount=q(1.0, "mmol")),
             participant("p_thf", mid("THF"), "solvent", e3, volume=q(5, "mL")),
             participant("p_5", mid("5"), "product", e3, product_role="target")]
    if loading:
        parts.insert(1, participant("p_tfoh", mid("TfOH"), "catalyst", e3,
                                    mol_percent=q(loading, "mol%"), eq_ref=mid("1a")))
    return reaction(
        suffix, e3, "specific_experiment", ["control", "mechanistic_study"],
        "allenol-to-enone control experiment", "acid-catalyzed allenol rearrangement",
        parts=parts,
        steps=[step("reaction", e3, [p["participant_id"] for p in parts if p["role"] != "product"],
                    "heat 1a in THF without phosphine oxide", q(100, "°C"), q(24, "h"), "Ar")],
        outcomes=[outcome("yield_5", "yield", [mid("5")], e3, q(yld, "%"), product_role="target")],
        changes=None if yld == 0 else [{"atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}],
                                       "before_order": 1, "after_order": 2, "change_type": "order_increased"}],
        coverage="not_applicable" if yld == 0 else "partial",
        unresolved=["Yield measurement basis is not reported on page 3."],
        notes=["No product 5 was observed and only starting material was recovered."] if yld == 0 else [],
        location=f"PDF page 3, Scheme 3a, TfOH {loading or 0} mol%")

control_1a_to_5("scheme3a_1a_to_5_with_tfoh", 30, 51)
control_1a_to_5("scheme3a_1a_to_5_without_tfoh", 0, 0)

addition_parts = [
    participant("p_5", mid("5"), "reactant", e3, amount=q(0.4, "mmol")),
    participant("p_2a", mid("2a"), "reactant", e3, amount=q(0.48, "mmol")),
    participant("p_tfoh", mid("TfOH"), "catalyst", e3, mol_percent=q(30, "mol%"), eq_ref=mid("5")),
    participant("p_thf", mid("THF"), "solvent", e3, volume=q(2, "mL")),
    participant("p_3a", mid("3a"), "product", e3, product_role="target"),
]
reaction(
    "scheme3b_enone_5_to_3a", e3, "specific_experiment", ["control", "mechanistic_study"],
    "addition of diphenylphosphine oxide 2a to enone 5", "phospha-Michael addition",
    parts=addition_parts,
    steps=[step("reaction", e3, ["p_5", "p_2a", "p_tfoh", "p_thf"], "heat the enone/phosphine-oxide mixture", q(100, "°C"), q(24, "h"), "Ar")],
    outcomes=[outcome("yield", "yield", [mid("3a")], e3, q(96, "%"), product_role="target")],
    changes=[
        {"atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "P"}],
         "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}],
         "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
        {"atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
         "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
    ],
    unresolved=["Yield measurement basis and complete proton transfer/coproduct balance are unreported."],
    location="PDF page 3, Scheme 3b")


# Page 4 proposed mechanism and operational general procedure.
reaction(
    "scheme4_proposed_mechanism", ev(3) + ev(4), "proposed_mechanism", [],
    "plausible mechanism", "proposed allenol rearrangement and phospha-Michael addition",
    notes=[
        "Authors propose dehydration to carbenium ion A, resonance form B, water attack to enol C, isomerization to enone D, phospha-Michael addition to E, and enol-keto tautomerism to 3.",
        "Intermediates A-E are proposed, not reported as isolated compounds.",
    ], changes=None, coverage="unresolved",
    unresolved=["Proposed intermediates are generic and lack experimental atom mapping; no mechanistic step is promoted to an observed isolated reaction."],
    location="PDF pages 3-4, Scheme 4")

gp_ev = ev(4)
gp_parts = [
    participant("p_allenol", generic_allenol, "reactant", gp_ev, amount=q(0.4, "mmol")),
    participant("p_poh", generic_phosphine, "reactant", gp_ev, amount=q(0.48, "mmol")),
    participant("p_tfoh", mid("TfOH"), "catalyst", gp_ev, amount=q(0.12, "mmol"), mol_percent=q(30, "mol%"), eq_ref=generic_allenol),
    participant("p_thf", mid("THF"), "solvent", gp_ev, volume=q(2, "mL")),
    participant("p_pe", mid("petroleum_ether"), "purification_agent", gp_ev, stage="purification"),
    participant("p_ea", mid("EtOAc"), "purification_agent", gp_ev, stage="purification"),
    participant("p_product", generic_product, "product", gp_ev, product_role="target"),
]
r = reaction(
    "general_procedure_gamma_ketophosphine_oxides", gp_ev, "generic_scheme", ["other"],
    "General Procedure for the Synthesis of gamma-Ketophosphine Oxides (3aa-3ap)",
    "allenol rearrangement/phospha-Michael coupling",
    parts=gp_parts,
    steps=[
        step("charge", gp_ev, ["p_allenol", "p_poh", "p_tfoh", "p_thf"],
             "add 2,3-allenol, H-phosphine oxide, and TfOH in THF to a Schlenk tube with stir bar under argon",
             atmosphere="Ar"),
        step("heat", gp_ev, ["p_allenol", "p_poh", "p_tfoh", "p_thf"],
             "stir in an oil bath", q(100, "°C"), q(24, "h"), "Ar"),
        step("cool", gp_ev, ["p_product"], "cool the resulting mixture to room temperature", "room temperature", stage="workup"),
        step("purify", gp_ev, ["p_product", "p_pe", "p_ea"],
             "purify by flash chromatography (petroleum ether/ethyl acetate)", stage="purification"),
    ], schemes=[current_scheme_id],
    unresolved=["The heading states 3aa-3ap although the same operation is referenced across the broader scope; no unavailable SI detail was inferred.",
                "Petroleum ether/ethyl acetate gradient composition is not specified."],
    location="PDF page 4, Experimental Section general procedure")
r["reaction_scale"] = f("2,3-allenol (0.4 mmol)", {"reference_molecule_id": generic_allenol, "amount": q(0.4, "mmol")}, gp_ev)
r["solvent_conditions"]["composition"] = f("THF", ["THF"], gp_ev)
r["solvent_conditions"]["total_volume"] = f("2 mL", q(2, "mL"), gp_ev)
r["procedure"]["original_text"] = f(
    "To a Schlenk tube equipped with a magnetic stir bar was added 2,3-allenols (0.4 mmol), H-phosphine oxides (0.48 mmol), and TfOH (0.12 mmol) in THF (2 mL) under argon ... 100 °C ... 24 h ... cooled ... flash chromatography (petroleum ether/ethyl acetate).",
    "page-4 general procedure transcription", gp_ev)
r["procedure"]["addition_order"] = f(
    "2,3-allenol + H-phosphine oxide + TfOH in THF; then heat", ["2,3-allenol", "H-phosphine oxide", "TfOH", "THF", "heat"], gp_ev)


bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
bundle["source_context"] = [
    {
        "context_id": f"{DOC}:bibliographic_metadata",
        "document_id": DOC,
        "kind": "bibliographic_metadata",
        "title": "TfOH-Catalyzed Phosphinylation of 2,3-Allenols into gamma-Ketophosphine Oxides",
        "citation": "J. Org. Chem. 2020, 85, 8185-8195",
        "doi": "10.1021/acs.joc.0c00328",
        "evidence_ids": ev(1),
    },
    {
        "context_id": f"{DOC}:characterization_methods",
        "document_id": DOC,
        "kind": "source_context",
        "description": "General instrumentation, chromatography media, and spectral methods on page 4 are characterization/purification context, not standalone reactions.",
        "evidence_ids": ev(4),
    },
    {
        "context_id": f"{DOC}:product_characterization_pages_4_6",
        "document_id": DOC,
        "kind": "characterization_only",
        "description": "NMR, IR, HRMS, melting-point, and physical-form data for 3a-3w were associated with their scope records where relevant; individual spectral peaks were not promoted to reactions.",
        "evidence_ids": ev(4) + ev(5) + ev(6),
    },
    {
        "context_id": f"{DOC}:si_unassociated",
        "document_id": DOC,
        "kind": "source_limitation",
        "description": "The article advertises Supporting Information, but source_mapping.json does not explicitly associate an SI file; no SI content was inferred.",
        "evidence_ids": ev(1),
    },
]


coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 1, "state": "attempted",
            "regions": [
                {"region": "title/abstract/current-work graphical abstract", "finding": "bibliographic metadata and one generic current transformation; not a concrete experiment", "evidence_ids": ev(1)},
                {"region": "Scheme 1a-c", "finding": "five displayed prior-work arrows retained as literature examples, not current-paper experiments", "evidence_ids": ev(1)},
                {"region": "Scheme 1d", "finding": "generic current TfOH-catalyzed allenol phosphinylation", "evidence_ids": ev(1)},
            ],
            "record_ids": [f"{DOC}:literature_scheme1a_base", f"{DOC}:literature_scheme1a_nhc", f"{DOC}:literature_scheme1b_ag", f"{DOC}:literature_scheme1b_eosin", f"{DOC}:literature_scheme1c_znotf", current_scheme_id],
            "unresolved": ["Prior-work generic structures were not expanded into concrete molecules; cited papers were not supplied."]
        },
        {
            "page": 2, "state": "attempted",
            "regions": [
                {"region": "Table 1 including footnotes", "finding": "19 independent optimization experiments, including a 0% catalyst-free row", "evidence_ids": ev(2)},
                {"region": "Table 2 full graphical scope", "finding": "35 concrete substrate-scope products 3a-3ai; table footnotes include 3r modification and 3ab dr method", "evidence_ids": ev(2)},
                {"region": "discussion around Tables 1-2", "finding": "optimization prose merged with rows; scope discussion continues on page 3", "evidence_ids": ev(2)},
            ],
            "record_ids": [f"{DOC}:table1_entry_{i:02d}" for i in range(1, 20)] + [f"{DOC}:table2_scope_3{x}" for x in letters] + [current_scheme_id],
            "unresolved": ["Table 1 and many compact scope yields do not state measurement basis; no basis was guessed.", "Absolute configurations of potential product stereocenters are unassigned."]
        },
        {
            "page": 3, "state": "attempted",
            "regions": [
                {"region": "continuing Table 2 discussion", "finding": "product-specific modifications, 3ab dr, X-ray statement for 3e, and one failed bulky tertiary substrate discussed", "evidence_ids": ev(3)},
                {"region": "Table 3 including footnotes", "finding": "seven P(O)H scope experiments 3aj-3ap, including inert 3ao and dr results for 3an/3ap", "evidence_ids": ev(3)},
                {"region": "Scheme 2a-b", "finding": "one gram-scale experiment and one bromination derivatization", "evidence_ids": ev(3)},
                {"region": "Scheme 3a-b", "finding": "two allenol-to-enone controls and one enone-to-3a mechanistic experiment", "evidence_ids": ev(3)},
                {"region": "mechanistic prose", "finding": "author interpretation linked to Scheme 4; not converted into observed intermediates", "evidence_ids": ev(3)},
            ],
            "record_ids": [f"{DOC}:table2_scope_3{x}" for x in ["r", "ab", "ac", "ad", "ae", "af", "ag", "ah", "ai"]] + [f"{DOC}:table3_scope_3{x}" for x, *_ in table3] + [f"{DOC}:scheme2a_gram_scale_3a", f"{DOC}:scheme2b_bromination_4", f"{DOC}:scheme3a_1a_to_5_with_tfoh", f"{DOC}:scheme3a_1a_to_5_without_tfoh", f"{DOC}:scheme3b_enone_5_to_3a", f"{DOC}:scheme4_proposed_mechanism"],
            "unresolved": ["Table 3 caption names 2a despite varying 2b-2h; row-specific 0.48 mmol assignments need expert confirmation.", "The failed bulky tertiary allenol is discussed without a compound label/product and remains narrative context rather than a fabricated reaction identity."]
        },
        {
            "page": 4, "state": "attempted",
            "regions": [
                {"region": "Scheme 4", "finding": "one proposed mechanism with non-isolated generic intermediates A-E", "evidence_ids": ev(4)},
                {"region": "Experimental Section general information", "finding": "instrumentation, chromatography, and characterization context; not standalone reactions", "evidence_ids": ev(4)},
                {"region": "general procedure", "finding": "one generic operational procedure with explicit charges, argon, 100 °C/24 h, cooling, and flash chromatography", "evidence_ids": ev(4)},
                {"region": "3a-3d experimental/characterization paragraphs", "finding": "isolated masses/yields merged into corresponding Table 2 scope records; spectra remain characterization context", "evidence_ids": ev(4)},
            ],
            "record_ids": [f"{DOC}:scheme4_proposed_mechanism", f"{DOC}:general_procedure_gamma_ketophosphine_oxides"] + [f"{DOC}:table2_scope_3{x}" for x in ["a", "b", "c", "d"]],
            "unresolved": ["The general-procedure heading says 3aa-3ap even though analogous conditions are referenced across the scope; unavailable SI was not used to resolve this."]
        },
        {
            "page": 5, "state": "attempted",
            "regions": [
                {"region": "3d continuation and 3e-3n experimental/characterization paragraphs", "finding": "isolated masses/yields merged into Table 2 records; spectra, HRMS, IR, and melting points not made into reactions", "evidence_ids": ev(5)},
            ],
            "record_ids": [f"{DOC}:table2_scope_3{x}" for x in ["d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n"]],
            "unresolved": ["No standalone synthesis arrows occur on this page; structural candidates remain needs_review despite successful parsing."]
        },
        {
            "page": 6, "state": "attempted",
            "regions": [
                {"region": "3n continuation and 3o-3w experimental/characterization paragraphs", "finding": "isolated masses/yields merged into Table 2 records; spectral data retained as context", "evidence_ids": ev(6)},
            ],
            "record_ids": [f"{DOC}:table2_scope_3{x}" for x in ["n", "o", "p", "q", "r", "s", "t", "u", "v", "w"]],
            "unresolved": ["The 3w characterization paragraph continues onto page 7; only the page-6 product name, mass, yield, and visible characterization portion were used."]
        },
    ]
}


RAW_DIR.mkdir(parents=True, exist_ok=True)
st.write(RAW_DIR / "extracted_bundle_pages_0001_0006.json", bundle)
st.write(RAW_DIR / "stereo_checks_pages_0001_0006.json", stereo_checks)
st.write(OUT, bundle)
st.write(ROOT / "runs" / "current" / "coverage" / "pages_0001_0006.json", coverage)

print(json.dumps({
    "molecules": len(bundle["molecules"]),
    "reactions": len(bundle["reactions"]),
    "source_context": len(bundle["source_context"]),
    "stereo_checks": len(stereo_checks),
}, indent=2))
