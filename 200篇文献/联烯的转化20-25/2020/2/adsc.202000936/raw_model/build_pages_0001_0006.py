#!/usr/bin/env python3
"""Build evidence-linked v3.1 records for adsc.202000936, PDF pages 1-6.

This is a source transcription helper, not a chemistry-certification program.
All structures inferred from drawings remain machine candidates and every
reaction remains needs_review after standard merge.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

from rdkit import Chem

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_2abd3d7449bbf0b6"
RUN = ROOT / "runs/current"
E = {n: [f"{DOC}:page{n}:image", f"{DOC}:page{n}:text"] for n in range(1, 7)}


def fact(raw=None, normalized=None, ev=(), method="original_extraction", derivation=None,
         status=None):
    return st.fact(raw, normalized, ev, method=method, derivation=derivation, status=status)


def missing(ev=()):
    return st.fact(evidence=ev)


def not_applicable(ev=()):
    return st.fact(evidence=ev, method="not_applicable", status="not_applicable")


bundle = json.loads((RUN / "submitted.json").read_text())
bundle["documents"] = json.loads((RUN / "dataset.json").read_text())["documents"]
bundle["molecules"] = []
bundle["reactions"] = []
bundle["source_context"] = []
molecules = {}
stereo_diagnostics = {}


def chiral_smiles(base: str, desired: str) -> str:
    """Assign the single alcohol stereocenter to desired CIP R/S."""
    mol = Chem.MolFromSmiles(base)
    candidates = []
    for atom in mol.GetAtoms():
        if (atom.GetAtomicNum() != 6 or atom.GetTotalNumHs() != 1 or
                atom.GetHybridization() != Chem.HybridizationType.SP3):
            continue
        if any(n.GetAtomicNum() == 8 for n in atom.GetNeighbors()):
            candidates.append(atom.GetIdx())
    if len(candidates) != 1:
        raise ValueError(f"Expected one alcohol stereocenter in {base}, found {candidates}")
    idx = candidates[0]
    for tag in (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW):
        trial = Chem.Mol(mol)
        trial.GetAtomWithIdx(idx).SetChiralTag(tag)
        Chem.AssignStereochemistry(trial, cleanIt=True, force=True)
        if trial.GetAtomWithIdx(idx).GetProp("_CIPCode") == desired:
            return Chem.MolToSmiles(trial, isomericSmiles=True)
    raise ValueError(f"Could not encode {desired} for {base}")


def add_molecule(mid, label, name, ev, smiles=None, chirality="achiral",
                 geometry="not_applicable", assignment=None, axial=None,
                 representation=None):
    if mid in molecules:
        return molecules[mid]
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = fact(label, label, ev)
    m["names"]["original"] = fact(name, name, ev)
    m["names"]["english"] = fact(name, name, ev)
    m["names"]["abbreviations"] = missing(ev)
    m["names"]["aliases"] = missing(ev)
    m["cas_numbers"] = missing(ev)
    m["external_properties"] = []
    m["stereochemistry"]["mixture"] = missing(ev)
    m["stereochemistry"]["axial"] = missing(ev) if axial is None else fact(
        axial, axial, ev)
    m["stereochemistry"]["assignment"] = missing(ev) if assignment is None else fact(
        assignment, assignment, ev)
    m["stereochemistry"]["chirality_status"] = fact(
        {"source_label": label}, chirality, ev, "inference",
        "Status assigned from the displayed structure and any explicit source configuration; it is not a sample-purity claim.")
    m["stereochemistry"]["geometry_status"] = fact(
        {"source_label": label}, geometry, ev, "inference",
        "Status assigned from the displayed structure; unlabelled E/Z geometry is retained as unspecified.")
    if smiles:
        diag = st.stereo_inspect(smiles)
        stereo_diagnostics[mid] = diag
        m["representation_type"] = "smiles"
        m["smiles"] = fact(
            {"source": "structure drawing", "compound_label": label}, smiles, ev,
            "structure_derivation",
            "Candidate SMILES transcribed from the source drawing and checked with standard stereo; parsing is not source verification.")
        m["generic_template"] = None
        m["stereochemistry"]["atom_centers"] = fact(
            {"source": "structure drawing"}, diag["atom_centers"], ev,
            "structure_derivation", "RDKit standard stereo diagnostics for the candidate SMILES.")
        m["stereochemistry"]["double_bonds"] = fact(
            {"source": "structure drawing"}, diag["double_bonds"], ev,
            "structure_derivation", "RDKit standard stereo diagnostics for the candidate SMILES.")
        rep = representation or ("partial" if chirality in ("unspecified", "partially_specified") or geometry == "unspecified" else "encoded")
        m["stereochemistry"]["representation_status"] = fact(
            {"source_label": label}, rep, ev, "inference",
            "Compared explicit source stereochemistry with stereochemical information encoded in the candidate SMILES.")
        if rep in ("encoded", "partial") and ("@" in smiles or rep == "encoded"):
            m["stereochemistry"]["isomeric_smiles"] = fact(
                {"candidate_smiles": smiles}, diag["canonical_isomeric_smiles"], ev,
                "structure_derivation", "Canonical isomeric SMILES from standard stereo.")
        else:
            m["stereochemistry"]["isomeric_smiles"] = missing(ev)
        m["stereochemistry"]["identity_key"] = fact(
            {"candidate_smiles": smiles}, diag["stereo_graph_key"], ev,
            "structure_derivation", "Stereo graph key from standard stereo; not a sample identity.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing(ev)
        m["generic_template"] = None
        m["stereochemistry"]["atom_centers"] = missing(ev)
        m["stereochemistry"]["double_bonds"] = missing(ev)
        m["stereochemistry"]["representation_status"] = fact(
            {"source_label": label}, representation or "unresolved", ev, "inference",
            "The compact or generic source structure was not converted to a complete candidate SMILES.")
        m["stereochemistry"]["isomeric_smiles"] = missing(ev)
        m["stereochemistry"]["identity_key"] = missing(ev)
    m["identity_note"] = fact(
        {"compound_label": label},
        "Machine transcription from pages 1-6; identity and any candidate structure require expert review.",
        ev)
    molecules[mid] = m
    bundle["molecules"].append(m)
    return m


def add_axial_catalyst(code, substituent, config="S_a"):
    suffix = code + ("_Ra" if config == "R_a" else "")
    ev = E[2] + E[3] + E[4] + E[5]
    return add_molecule(
        f"m_cat_{suffix}", f"({config})-7{code}",
        f"atropisomeric bipyridine N,N'-dioxide 7{code} (R = {substituent})",
        ev, chirality="specified", geometry="not_applicable",
        axial=[{"axis": "biaryl atropisomeric axis", "configuration": config, "source_label": f"7{code}"}],
        assignment={"absolute_axial_configuration": config}, representation="not_encoded")


# Common reagents and solvents.
add_molecule("m_allenylsilane_2", "2", "allenyltrichlorosilane", E[1] + E[2] + E[5],
             "C=C=C[Si](Cl)(Cl)Cl", chirality="unspecified", representation="partial")
add_molecule("m_dipea", "Et-i-Pr2N", "N,N-diisopropylethylamine", E[3] + E[5], "CCN(C(C)C)C(C)C")
solvents = {
    "DCM": ("m_dcm", "dichloromethane", "ClCCl"),
    "EtCN": ("m_etcn", "propionitrile", "CCC#N"),
    "THF": ("m_thf", "tetrahydrofuran", "C1CCOC1"),
    "CHCl3": ("m_chcl3", "chloroform", "ClC(Cl)Cl"),
    "EtOAc": ("m_etoac", "ethyl acetate", "CCOC(C)=O"),
    "Ether": ("m_ether", "diethyl ether", "CCOCC"),
    "Toluene": ("m_toluene", "toluene", "Cc1ccccc1"),
}
for label, (mid, name, smi) in solvents.items():
    add_molecule(mid, label, name, E[3] + E[5], smi)

catalysts = {
    "a": "Ph", "b": "CF3", "c": "CN", "d": "tBu",
    "e": "3,5-(CF3)2C6H3", "f": "4-ClC6H4", "g": "4-MeOC6H4",
    "h": "2-naphthyl", "i": "2-furyl", "j": "2-thienyl",
}
for c, r_group in catalysts.items():
    add_axial_catalyst(c, r_group)
add_axial_catalyst("a", "Ph", "R_a")
add_molecule("m_cat_7_family", "(S_a)-7a-j", "atropisomeric bipyridine N,N'-dioxide family 7a-j",
             E[2], chirality="specified", geometry="not_applicable",
             axial=[{"axis": "biaryl atropisomeric axis", "configuration": "S_a", "source_label": "7a-j"}],
             assignment={"absolute_axial_configuration": "S_a"}, representation="not_encoded")

# Scheme 2 labels retained without guessed complete structures.
for mid, label, name in [
    ("m_ketone_8_family", "8a,b,d-j", "methyl ketone family 8a,b,d-j"),
    ("m_pyridine", "Py", "pyridine"), ("m_iodine", "I2", "iodine"),
    ("m_pyridinium_9_family", "9a,b,d-j", "pyridinium salt family 9a,b,d-j"),
    ("m_myrtenal_10", "10", "R-(-)-myrtenal"),
    ("m_ammonium_acetate", "NH4OAc", "ammonium acetate"),
    ("m_formamide", "formamide", "formamide"),
    ("m_pyridine_11_family", "11a,b,d-j", "terpene-derived bipyridine precursor family 11a,b,d-j"),
    ("m_11i", "11i", "2-furyl terpene-derived bipyridine precursor 11i"),
    ("m_11c", "11c", "nitrile terpene-derived bipyridine precursor 11c"),
    ("m_mcpba", "mCPBA", "meta-chloroperoxybenzoic acid"),
    ("m_noxide_12_family", "12a-j", "pyridine N-oxide family 12a-j"),
    ("m_oxidant_unspecified", "[O]", "unspecified oxidant"),
]:
    add_molecule(mid, label, name, E[2], representation="unresolved",
                 chirality="unresolved", geometry="unresolved")

# Substrate/product structures from Table 4.  E/Z remains unspecified where shown without a descriptor.
substrates = {
    "a": ("benzaldehyde", "O=Cc1ccccc1", "not_applicable"),
    "b": ("2-nitrobenzaldehyde", "O=Cc1ccccc1[N+](=O)[O-]", "not_applicable"),
    "c": ("4-nitrobenzaldehyde", "O=Cc1ccc([N+](=O)[O-])cc1", "not_applicable"),
    "d": ("4-(trifluoromethyl)benzaldehyde", "O=Cc1ccc(C(F)(F)F)cc1", "not_applicable"),
    "e": ("2-methoxybenzaldehyde", "COc1ccccc1C=O", "not_applicable"),
    "f": ("4-methoxybenzaldehyde", "COc1ccc(C=O)cc1", "not_applicable"),
    "g": ("2,4-dimethoxybenzaldehyde", "O=Cc1c(OC)cc(OC)cc1", "not_applicable"),
    "h": ("2,3-dimethoxybenzaldehyde", "O=Cc1c(OC)c(OC)ccc1", "not_applicable"),
    "i": ("2-fluorobenzaldehyde", "O=Cc1ccccc1F", "not_applicable"),
    "j": ("4-fluorobenzaldehyde", "O=Cc1ccc(F)cc1", "not_applicable"),
    "k": ("2-chlorobenzaldehyde", "O=Cc1ccccc1Cl", "not_applicable"),
    "l": ("3-chlorobenzaldehyde", "O=Cc1cccc(Cl)c1", "not_applicable"),
    "m": ("4-chlorobenzaldehyde", "O=Cc1ccc(Cl)cc1", "not_applicable"),
    "n": ("2-bromobenzaldehyde", "O=Cc1ccccc1Br", "not_applicable"),
    "o": ("4-bromobenzaldehyde", "O=Cc1ccc(Br)cc1", "not_applicable"),
    "p": ("2-methylbenzaldehyde", "Cc1ccccc1C=O", "not_applicable"),
    "q": ("3-methylbenzaldehyde", "Cc1cccc(C=O)c1", "not_applicable"),
    "r": ("4-methylbenzaldehyde", "Cc1ccc(C=O)cc1", "not_applicable"),
    "s": ("2-naphthaldehyde", "O=Cc1ccc2ccccc2c1", "not_applicable"),
    "t": ("furan-2-carbaldehyde (furfural)", "O=Cc1ccco1", "not_applicable"),
    "u": ("furan-3-carbaldehyde", "O=Cc1ccoc1", "not_applicable"),
    "v": ("thiophene-2-carbaldehyde", "O=Cc1cccs1", "not_applicable"),
    "w": ("5-methylthiophene-2-carbaldehyde", "Cc1ccc(C=O)s1", "not_applicable"),
    "x": ("cinnamaldehyde", "O=CC=Cc1ccccc1", "unspecified"),
    "y": ("alpha-methylcinnamaldehyde", "O=CC(C)=Cc1ccccc1", "unspecified"),
    "z": ("1-cyclohexene-1-carboxaldehyde", "O=CC1=CCCCC1", "not_applicable"),
    "aa": ("3-phenylpropanal", "O=CCCc1ccccc1", "not_applicable"),
    "ab": ("cyclohexanecarboxaldehyde", "O=CC1CCCCC1", "not_applicable"),
}
product_bases = {
    "a": "C#CCC(O)c1ccccc1", "b": "C#CCC(O)c1ccccc1[N+](=O)[O-]",
    "c": "C#CCC(O)c1ccc([N+](=O)[O-])cc1", "d": "C#CCC(O)c1ccc(C(F)(F)F)cc1",
    "e": "COc1ccccc1C(O)CC#C", "f": "COc1ccc(C(O)CC#C)cc1",
    "g": "C#CCC(O)c1c(OC)cc(OC)cc1", "h": "C#CCC(O)c1c(OC)c(OC)ccc1",
    "i": "C#CCC(O)c1ccccc1F", "j": "C#CCC(O)c1ccc(F)cc1",
    "k": "C#CCC(O)c1ccccc1Cl", "l": "C#CCC(O)c1cccc(Cl)c1",
    "m": "C#CCC(O)c1ccc(Cl)cc1", "n": "C#CCC(O)c1ccccc1Br",
    "o": "C#CCC(O)c1ccc(Br)cc1", "p": "Cc1ccccc1C(O)CC#C",
    "q": "Cc1cccc(C(O)CC#C)c1", "r": "Cc1ccc(C(O)CC#C)cc1",
    "s": "C#CCC(O)c1ccc2ccccc2c1", "t": "C#CCC(O)c1ccco1",
    "u": "C#CCC(O)c1ccoc1", "v": "C#CCC(O)c1cccs1",
    "w": "C#CCC(O)c1ccc(C)s1", "x": "C#CCC(O)C=Cc1ccccc1",
    "y": "C#CCC(O)C(C)=Cc1ccccc1", "z": "C#CCC(O)C1=CCCCC1",
    "aa": "C#CCC(O)CCc1ccccc1", "ab": "C#CCC(O)C1CCCCC1",
}
for code, (name, smi, geom) in substrates.items():
    ev = E[5] if code != "a" else E[1] + E[2] + E[3] + E[4] + E[5]
    add_molecule(f"m_1{code}", f"1{code}", name, ev, smi,
                 chirality="achiral", geometry=geom,
                 representation="partial" if geom == "unspecified" else "encoded")

# Optimization product with unassigned major-enantiomer identity.
add_molecule("m_3a_unspecified", "3a", "1-phenylpent-4-yn-1-ol (configuration unspecified)",
             E[2] + E[3] + E[4], product_bases["a"], chirality="unspecified",
             representation="partial")
add_molecule("m_3a_S", "(S)-3a", "(S)-1-phenylpent-4-yn-1-ol", E[3],
             chiral_smiles(product_bases["a"], "S"), chirality="specified",
             assignment={"absolute_configuration": "S", "source": "Table 1 footnote d"})

for code, base in product_bases.items():
    desired = "S" if code in ("b", "y") else "R"
    assignment_source = "Table 4 label" if code == "y" else (
        "opposite absolute configuration stated in prose and opposite wedge in Table 4" if code == "b" else
        "Table 4 footnote: R unless stated otherwise")
    geom = "unspecified" if code in ("x", "y") else "not_applicable"
    product_ev = E[3] + E[5] if code == "a" else E[5]
    add_molecule(f"m_3{code}_{desired}", f"3{code}", f"homopropargylic alcohol 3{code}", product_ev,
                 chiral_smiles(base, desired), chirality="specified", geometry=geom,
                 assignment={"absolute_configuration": desired, "source": assignment_source},
                 representation="partial" if geom == "unspecified" else "encoded")


def add_participant(r, pid, mid, role, ev, eq=None, molp=None, amount=None, mass=None, stage="reaction"):
    p = st.new_participant(pid, mid, role, ev, stage=stage)
    if role == "product":
        p["product_role"] = fact("main product", "main product", ev)
    else:
        p["product_role"] = not_applicable(ev)
    if eq is not None:
        p["charges"]["equivalents"] = fact(f"{eq} equiv", {"number": eq, "unit": "equiv"}, ev)
        p["equivalence_reference"] = fact("relative to aldehyde", r["_substrate_mid"], ev)
    if molp is not None:
        p["charges"]["mol_percent"] = fact(f"{molp} mol%", {"number": molp, "unit": "mol%"}, ev)
        p["equivalence_reference"] = fact("relative to aldehyde", r["_substrate_mid"], ev)
    if amount is not None:
        p["charges"]["amount"] = fact(f"{amount} mmol", {"number": amount, "unit": "mmol"}, ev)
    if mass is not None:
        p["charges"]["mass"] = fact(f"{mass} g", {"number": mass, "unit": "g"}, ev)
    r["participants"].append(p)
    return p


def propargyl_changes(ev):
    raw = "net bond changes inferred from Scheme 1/Table 4 aldehyde + allenyltrichlorosilane -> homopropargylic alcohol"
    changes = [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}], "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
        {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 2, "after_order": 3, "change_type": "order_increased"},
        {"change_id": "b4", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
        {"change_id": "b5", "atoms": [{"atom_map": None, "element": "O"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
        {"change_id": "b6", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "Si"}], "before_order": 1, "after_order": 0, "change_type": "broken"},
    ]
    labels = [
        {"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1"]},
        {"category_number": 21, "code": "TB-C-C", "supporting_change_ids": ["b3"]},
        {"category_number": 24, "code": "O-H", "supporting_change_ids": ["b5"]},
    ]
    return raw, changes, labels


def base_reaction(rid, ev, record_type="specific_experiment", categories=None,
                  substrate="m_1a", product="m_3a_unspecified"):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["_substrate_mid"] = substrate
    r["experiment_categories"] = fact(categories or [], categories or [], ev)
    r["location"] = fact(rid, {"pages": sorted({int(x.split(":page")[1].split(":")[0]) for x in ev if ":page" in x})}, ev)
    r["reaction_name"] = fact("asymmetric propargylation", "asymmetric propargylation", ev)
    r["transformation_type"] = fact(
        "aldehyde + allenyltrichlorosilane -> homopropargylic alcohol",
        "enantioselective nucleophilic propargylation of an aldehyde", ev)
    r["general_scheme_ids"] = [] if rid == "scheme1_general_propargylation" else [f"{DOC}:scheme1_general_propargylation"]
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    raw, changes, labels = propargyl_changes(ev)
    r["classification"]["mapped_reaction_smiles"] = missing(ev)
    r["classification"]["bond_changes"] = [fact(raw, c, ev, "structure_derivation",
        "Net graph change inferred from the displayed generic transformation; atom correspondence is not fully mapped.") for c in changes]
    r["classification"]["labels"] = fact(raw, labels, ev, "structure_derivation",
        "Labels cover the inferred new C-C connection, alkyne bond-order promotion, and implicit O-H formation.")
    r["classification"]["coverage"] = fact(raw, "partial", ev, "inference",
        "Allene carbon correspondence, proton source, and omitted silicon/chloride coproducts prevent complete mapping.")
    r["classification"]["limitations"] = [
        "No complete atom-mapped reaction was asserted.",
        "The O-H label is a net hydrogen gain; the proton source is unreported on pages 1-6.",
        "C-Si cleavage and allene-to-alkyne bond-order changes are retained, but coproduct identities are unavailable.",
    ]
    r["classification"]["model_votes"] = []
    bundle["reactions"].append(r)
    return r


def set_conditions(r, ev, solvent_label, catalyst_mid, loading, time_h, scale_mmol,
                   temperature=-60, inherited=False, include_reagent_doses=True, mass_g=None):
    method = "inherited_general_condition" if inherited else "original_extraction"
    raw_basis = "Table 1 model-reaction conditions carried into the stated optimization series" if inherited else "table row and footnote"
    add_participant(r, "substrate", r["_substrate_mid"], "reactant", ev,
                    amount=scale_mmol, mass=mass_g)
    add_participant(r, "allenylsilane", "m_allenylsilane_2", "reactant", ev,
                    eq=1.5 if include_reagent_doses else None)
    if include_reagent_doses:
        add_participant(r, "base", "m_dipea", "reagent", ev, eq=1.5)
    add_participant(r, "catalyst", catalyst_mid, "catalyst", ev, molp=loading)
    solvent_mid = solvents[solvent_label][0]
    add_participant(r, "solvent", solvent_mid, "solvent", ev)
    r["reaction_scale"] = fact(
        f"{scale_mmol} mmol aldehyde scale" + (f" ({mass_g} g)" if mass_g is not None else ""),
        {"reference_molecule_id": r["_substrate_mid"], "amount": {"number": scale_mmol, "unit": "mmol"},
         **({"mass": {"number": mass_g, "unit": "g"}} if mass_g is not None else {})},
        ev, method, "Inherited only where the prose explicitly continues the same model reaction." if inherited else None)
    r["solvent_conditions"]["composition"] = fact(solvent_label, [{"molecule_id": solvent_mid}], ev, method,
        "Solvent series belongs to the same model reaction." if inherited else None)
    step = st.new_step("reaction", ev)
    step["participant_ids"] = [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"]
    step["operation"] = fact(raw_basis, "stir/hold under stated conditions", ev, method,
                             "Inherited from Table 1 model reaction." if inherited else None)
    step["temperature"] = fact(f"{temperature} °C", {"number": temperature, "unit": "°C"}, ev, method,
                               "Inherited from Table 1 model reaction." if inherited else None)
    step["time"] = fact(f"{time_h} h", {"number": time_h, "unit": "h"}, ev)
    step["atmosphere"] = fact("under argon", "argon", ev)
    step["pressure"] = missing(ev)
    r["steps"] = [step]
    r["procedure"]["original_text"] = fact(raw_basis, raw_basis, ev, method,
        "Table 1 supplies the shared model-reaction doses." if inherited else None)
    r["procedure"]["general_procedure_reference"] = fact(
        "Scheme 1 / displayed table conditions", f"{DOC}:scheme1_general_propargylation", ev)


def add_conversion(r, value, ev):
    o = st.new_outcome("conversion", "conversion", [r["_substrate_mid"]], ev)
    o["product_role"] = not_applicable(ev)
    normalized = ({"number": value, "unit": "%"} if isinstance(value, (int, float)) else value)
    o["value"] = fact(f"{value}%" if isinstance(value, (int, float)) else str(value), normalized, ev)
    o["yield_type"] = not_applicable(ev)
    o["measurement_basis"] = fact("Conversion was calculated from GC data.", "GC", ev)
    o["scope"] = fact("aldehyde conversion", {"metric": "conversion", "target": r["_substrate_mid"]}, ev)
    r["outcomes"].append(o)


def add_ee(r, value, product_mid, ev, config_order="major enantiomer not assigned"):
    o = st.new_outcome("ee", "selectivity", [product_mid], ev)
    o["product_role"] = fact("main product", "main product", ev)
    o["value"] = fact(f"{value}% ee" if isinstance(value, (int, float)) else str(value),
                      {"number": value, "unit": "%"} if isinstance(value, (int, float)) else value, ev)
    o["yield_type"] = not_applicable(ev)
    o["measurement_basis"] = fact("Determined by GC on a chiral column.", "chiral GC", ev)
    o["scope"] = fact("ee", {"metric": "ee", "configuration_order": config_order}, ev)
    r["outcomes"].append(o)


def add_yield_and_scope_ee(r, product_mid, yld, ee, ev):
    oy = st.new_outcome("isolated_yield", "yield", [product_mid], ev)
    oy["product_role"] = fact("main product", "main product", ev)
    oy["value"] = fact(f"{yld}%", {"number": yld, "unit": "%"}, ev)
    oy["yield_type"] = fact("The yield is shown for isolated yield.", "isolated yield", ev)
    oy["measurement_basis"] = fact("The yield is shown for isolated yield.", "isolated", ev)
    oy["scope"] = fact("product 3", {"metric": "isolated yield"}, ev)
    r["outcomes"].append(oy)
    oe = st.new_outcome("ee", "selectivity", [product_mid], ev)
    oe["product_role"] = fact("main product", "main product", ev)
    oe["value"] = fact(f"ee {ee}%", {"number": ee, "unit": "%"}, ev)
    oe["yield_type"] = not_applicable(ev)
    oe["measurement_basis"] = fact(
        "Enantiopurity was determined by chiral GC or chiral HPLC, see SI for details.",
        "chiral GC or chiral HPLC; compound-specific method unavailable without SI", ev)
    desired = "S" if product_mid.endswith("_S") else "R"
    oe["scope"] = fact("ee", {"metric": "ee", "major_absolute_configuration": desired,
                               "configuration_order": f"{desired} major; minor enantiomer not separately identified"}, ev)
    r["outcomes"].append(oe)


# Page 2 generic propargylation scheme.
r = base_reaction("scheme1_general_propargylation", E[2], "generic_scheme", ["other"])
add_participant(r, "aldehyde_1", "m_1a", "reactant", E[2])
add_participant(r, "allenylsilane_2", "m_allenylsilane_2", "reactant", E[2])
add_participant(r, "alcohol_3", "m_3a_unspecified", "product", E[2])
r["procedure"]["original_text"] = fact("Scheme 1 reaction arrow", "generic aldehyde propargylation", E[2])
r["review"]["unresolved"] = ["Scheme 1 is a generic transformation and includes prior-work catalysts 4-6; no current-paper yield is assigned to it."]
r.pop("_substrate_mid")


def generic_scheme2(rid, reactant_ids, product_id, reagent_ids, arrow_text):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, E[2])
    r["record_type"] = "generic_scheme"
    r["experiment_categories"] = fact(["other"], ["other"], E[2])
    r["location"] = fact("Scheme 2", {"page": 2, "arrow": arrow_text}, E[2])
    r["reaction_name"] = fact("Synthesis of atropisomeric bipyridine N,N'-dioxides 7", arrow_text, E[2])
    r["transformation_type"] = fact(arrow_text, arrow_text, E[2])
    for i, mid in enumerate(reactant_ids):
        add_participant(r, f"reactant_{i+1}", mid, "reactant", E[2])
    for i, mid in enumerate(reagent_ids):
        add_participant(r, f"reagent_{i+1}", mid, "reagent", E[2])
    add_participant(r, "product", product_id, "product", E[2])
    step = st.new_step("reaction", E[2]); step["participant_ids"] = [p["participant_id"] for p in r["participants"]]
    step["operation"] = fact(arrow_text, arrow_text, E[2]); r["steps"] = [step]
    r["procedure"]["original_text"] = fact("Scheme 2 arrow", arrow_text, E[2])
    r["classification"]["coverage"] = fact("compact generic Scheme 2 structures", "unresolved", E[2], "inference",
        "No complete atom mapping was attempted for the compact generic catalyst synthesis.")
    r["classification"]["limitations"] = ["Generic family structures and missing SI prevent defensible atom-level classification."]
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
                   "unresolved": ["Scheme 2 gives a generic arrow without yield or experimental conditions; SI was not provided.",
                                  "Compact catalyst-family structures were not forced into guessed SMILES."]}
    r["completeness"] = "partial"
    bundle["reactions"].append(r)


generic_scheme2("scheme2_step1_8_to_9", ["m_ketone_8_family", "m_pyridine"], "m_pyridinium_9_family", ["m_iodine"], "8a,b,d-j to pyridinium salts 9a,b,d-j using I2 and pyridine")
generic_scheme2("scheme2_step2_9_to_11", ["m_pyridinium_9_family", "m_myrtenal_10"], "m_pyridine_11_family", ["m_ammonium_acetate", "m_formamide"], "Kroenke annulation 9 to 11 with R-(-)-myrtenal, NH4OAc and formamide")
generic_scheme2("scheme2_step3_11_to_12", ["m_pyridine_11_family"], "m_noxide_12_family", ["m_mcpba"], "oxidation of 11a,b,d-j with mCPBA to 12a-j")
generic_scheme2("scheme2_step4_12_to_7", ["m_noxide_12_family"], "m_cat_7_family", ["m_oxidant_unspecified"], "oxidative assembly of 12a-j to atropisomeric N,N'-dioxides (S_a)-7a-j")

# Page 3, Table 1: 11 catalyst-screening entries.
table1 = [
    (1, "m_cat_a", 99, 96, "m_3a_R", "R major inferred from opposite-series statement and Table 4 assignment"),
    (2, "m_cat_a_Ra", 99, 96, "m_3a_S", "S major; Table 1 footnote d"),
    (3, "m_cat_b", "trace", "not determined", "m_3a_unspecified", "major enantiomer not assigned"),
    (4, "m_cat_c", {"number": 5, "unit": "%", "qualifier": "le", "raw_text": "≤5%"}, 69, "m_3a_unspecified", "major enantiomer not assigned"),
    (5, "m_cat_d", {"number": 5, "unit": "%", "qualifier": "le", "raw_text": "≤5%"}, 50, "m_3a_unspecified", "major enantiomer not assigned"),
    (6, "m_cat_e", 99, 93, "m_3a_unspecified", "major enantiomer not assigned"),
    (7, "m_cat_f", 99, 97, "m_3a_unspecified", "major enantiomer not assigned"),
    (8, "m_cat_g", 70, 93, "m_3a_unspecified", "major enantiomer not assigned"),
    (9, "m_cat_h", 99, 95, "m_3a_unspecified", "major enantiomer not assigned"),
    (10, "m_cat_i", 60, 60, "m_3a_unspecified", "major enantiomer not assigned"),
    (11, "m_cat_j", 99, 76, "m_3a_unspecified", "major enantiomer not assigned"),
]
for entry, cat, conv, ee, product, order in table1:
    r = base_reaction(f"table1_entry_{entry:02d}", E[3], categories=["optimization"], product=product)
    set_conditions(r, E[3], "DCM", cat, 10, 18, 0.5)
    add_participant(r, "product", product, "product", E[3])
    add_conversion(r, conv, E[3]); add_ee(r, ee, product, E[3], order)
    r["review"]["unresolved"] = ["No isolated yield is reported for this optimization row.", "SI-only experimental detail is unavailable."]
    r.pop("_substrate_mid")

# Page 3, Table 2: 11 solvent-screening entries.
table2 = [
    (1, "m_cat_a", "DCM", 99, 96), (2, "m_cat_a", "EtCN", {"number": 5, "unit": "%", "qualifier": "le", "raw_text": "≤5%"}, 76),
    (3, "m_cat_e", "DCM", 99, 93), (4, "m_cat_e", "EtCN", 14, 93),
    (5, "m_cat_g", "DCM", 70, 93), (6, "m_cat_g", "EtCN", 99, 72),
    (7, "m_cat_g", "THF", 85, 37), (8, "m_cat_g", "CHCl3", {"number": 5, "unit": "%", "qualifier": "le", "raw_text": "≤5%"}, 80),
    (9, "m_cat_g", "EtOAc", 43, 74), (10, "m_cat_g", "Ether", 7, 40),
    (11, "m_cat_g", "Toluene", 9, 40),
]
for entry, cat, solv, conv, ee in table2:
    r = base_reaction(f"table2_entry_{entry:02d}", E[3], categories=["optimization"])
    set_conditions(r, E[3], solv, cat, 10, 18, 0.5, inherited=True)
    add_participant(r, "product", "m_3a_unspecified", "product", E[3])
    add_conversion(r, conv, E[3]); add_ee(r, ee, "m_3a_unspecified", E[3])
    r["review"]["unresolved"] = ["Allenylsilane/base doses are explicitly inherited from the adjacent Table 1 model reaction; compound-specific SI detail is unavailable.", "Absolute configuration of the major enantiomer is not assigned in Table 2."]
    r.pop("_substrate_mid")

# Page 4, Table 3: 14 loading/time optimization entries.
table3 = [
    (1, "m_cat_e", 10, 0.25, 99, 91), (2, "m_cat_e", 5, 0.25, 99, 96),
    (3, "m_cat_e", 2, 0.5, 98, None), (4, "m_cat_e", 2, 1, 99, 97),
    (5, "m_cat_e", 1, 2, 66, None), (6, "m_cat_e", 1, 3, 93, None),
    (7, "m_cat_e", 1, 3.5, 99, 98), (8, "m_cat_f", 10, 0.5, 99, 98),
    (9, "m_cat_f", 5, 0.5, 90, None), (10, "m_cat_f", 5, 1, 99, 98),
    (11, "m_cat_f", 2, 3, 69, None), (12, "m_cat_f", 2, 8, 99, 97),
    (13, "m_cat_f", 1, 3, 35, None), (14, "m_cat_f", 1, 15, 58, 94),
]
for entry, cat, loading, hours, conv, ee in table3:
    r = base_reaction(f"table3_entry_{entry:02d}", E[4], categories=["optimization"])
    set_conditions(r, E[3] + E[4], "DCM", cat, loading, hours, 0.5, inherited=True)
    add_participant(r, "product", "m_3a_unspecified", "product", E[4])
    add_conversion(r, conv, E[4])
    if ee is not None:
        add_ee(r, ee, "m_3a_unspecified", E[4])
    r["review"]["unresolved"] = ["Table 3 omits ee for rows where no ee value is printed; it was not imputed.", "Model-reaction doses/solvent/scale are inherited from Table 1 context; SI is unavailable."]
    r.pop("_substrate_mid")

# Page 5, Table 4: 29 distinct scope experiments, including duplicate 3a conditions.
scope_rows = [
    ("3a_small", "a", 92, 98, None), ("3a_scale", "a", 81, 92, "scale"),
    ("3b", "b", 84, 82, None), ("3c", "c", 86, 82, None), ("3d", "d", 96, 84, None),
    ("3e", "e", 90, 96, None), ("3f", "f", 88, 96, None), ("3g", "g", 83, 96, None),
    ("3h", "h", 81, 93, None), ("3i", "i", 90, 94, None), ("3j", "j", 93, 96, None),
    ("3k", "k", 88, 92, None), ("3l", "l", 94, 94, None), ("3m", "m", 87, 95, None),
    ("3n", "n", 94, 91, None), ("3o", "o", 92, 95, None), ("3p", "p", 95, 84, None),
    ("3q", "q", 91, 96, None), ("3r", "r", 93, 98, None), ("3s", "s", 93, 96, None),
    ("3t", "t", 75, 80, None), ("3u", "u", 79, 90, None), ("3v", "v", 87, 95, None),
    ("3w", "w", 79, 93, None), ("3x", "x", 79, 90, None), ("3y", "y", 82, 95, None),
    ("3z", "z", 89, 98, None), ("3aa", "aa", 68, 57, "modified"), ("3ab", "ab", 72, 57, None),
]
for row_id, code, yld, ee, special in scope_rows:
    desired = "S" if code in ("b", "y") else "R"
    product_mid = f"m_3{code}_{desired}"
    r = base_reaction(f"table4_{row_id}", E[5], categories=["substrate_scope"] + (["scale_up"] if special == "scale" else []), substrate=f"m_1{code}", product=product_mid)
    if special == "modified":
        set_conditions(r, E[5], "EtCN", "m_cat_g", 10, 18, 0.5)
        r["procedure"]["modifications"] = fact(
            "The reaction was carried out with catalyst 7g (10 mol%) in EtCN at -60 °C for 18 h.",
            {"catalyst": "7g", "loading": {"number": 10, "unit": "mol%"}, "solvent": "EtCN", "temperature": {"number": -60, "unit": "°C"}, "time": {"number": 18, "unit": "h"}}, E[5])
    elif special == "scale":
        set_conditions(r, E[5], "DCM", "m_cat_e", 1, 8, 9.4, mass_g=1.0)
        r["procedure"]["modifications"] = fact("The reaction was performed on a 1.0 g (9.4 mmol) scale.", {"scale_up": True}, E[5])
    else:
        set_conditions(r, E[5], "DCM", "m_cat_e", 1, 8, 0.5)
    add_participant(r, "product", product_mid, "product", E[5])
    add_yield_and_scope_ee(r, product_mid, yld, ee, E[5])
    r["review"]["unresolved"] = [
        "SI was advertised but not supplied or explicitly associated; compound-specific operations and the choice of chiral GC versus HPLC remain unavailable.",
        "Candidate substrate/product SMILES and stereochemical encoding are machine-derived from the Table 4 drawing and require expert review.",
    ]
    if code in ("x", "y"):
        r["review"]["unresolved"].append("The source drawing does not label alkene E/Z geometry; it remains unspecified despite encoding the alcohol-center assignment.")
    r.pop("_substrate_mid")

# Page 2 prose reports a concrete chemical modification but defers all details to unavailable SI.
r = st.new_reaction(f"{DOC}:prose_11i_to_11c", DOC, E[2])
r["record_type"] = "specific_experiment"
r["experiment_categories"] = fact(["derivatization"], ["derivatization"], E[2])
r["location"] = fact("Results and Discussion, page 2", {"page": 2, "region": "right column prose above Scheme 2"}, E[2])
r["reaction_name"] = fact("chemical modification of 2-furyl derivative 11i to nitrile 11c", "11i to 11c derivatization", E[2])
r["transformation_type"] = fact("2-furyl derivative 11i was transformed into nitrile 11c", "furyl-to-nitrile scaffold modification", E[2])
add_participant(r, "reactant", "m_11i", "reactant", E[2])
add_participant(r, "product", "m_11c", "product", E[2])
step = st.new_step("reported_transformation", E[2]); step["participant_ids"] = ["reactant", "product"]
step["operation"] = fact("2-furyl derivative 11i was transformed into nitrile 11c, see SI for details", "chemical modification; details unavailable", E[2])
r["steps"] = [step]
r["procedure"]["original_text"] = fact("2-furyl derivative 11i was transformed into nitrile 11c, see SI for details", "reported transformation with SI-only details", E[2])
r["classification"]["coverage"] = fact("structures/procedure deferred to SI", "unresolved", E[2], "inference",
    "No defensible atom mapping or bond-change assignment is possible from the compact prose reference alone.")
r["classification"]["limitations"] = ["SI not supplied; exact structures, reagents, conditions, outcome and atom mapping remain unresolved."]
r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
               "unresolved": ["The main text reports the transformation, but every experimental detail and any yield are deferred to unavailable SI.",
                              "Reactant/product compact structures were not converted into guessed SMILES."]}
r["completeness"] = "partial"
bundle["reactions"].append(r)

# Source context: prose/figures that are not current-paper synthetic experiments.
bundle["source_context"] = [
    {"context_id": f"{DOC}:page1_prior_work", "document_id": DOC, "kind": "literature_context",
     "evidence_ids": E[1], "note": "Introduction cites earlier asymmetric propargylations and ee values; cited work was not converted into current-paper experiments."},
    {"context_id": f"{DOC}:scheme1_prior_catalysts", "document_id": DOC, "kind": "literature_and_prediction_context",
     "evidence_ids": E[2], "note": "Catalysts 4 and 5 are prior literature; structures 6a-d and printed ee values are computational predictions, not experiments in this paper."},
    {"context_id": f"{DOC}:figure1_computation", "document_id": DOC, "kind": "computational_mechanism",
     "evidence_ids": E[4], "note": "BP1-BP5 and TS-6a/TS-7a are computed stereochemical models; no synthetic experiment record created."},
    {"context_id": f"{DOC}:figure2_computation", "document_id": DOC, "kind": "computational_mechanism",
     "evidence_ids": E[6], "note": "Computed TS structures and relative energies for 1a-c with 7e; retained as mechanism context, not experimental proof."},
    {"context_id": f"{DOC}:mass_spectra_inventory", "document_id": DOC, "kind": "coverage_note",
     "evidence_ids": E[1] + E[2] + E[3] + E[4] + E[5] + E[6], "note": "No standalone mass-spectrum peak table, molecular-weight list, or fragment-abundance observation appears on pages 1-6."},
]

# Save diagnostics and raw candidate bundle before formal submission.
(RUN / "raw_model").mkdir(parents=True, exist_ok=True)
(RUN / "raw_model/stereo_diagnostics_pages_0001_0006.json").write_text(
    json.dumps(stereo_diagnostics, ensure_ascii=False, indent=2) + "\n")
mentions = st.source_mentions(json.loads((RUN / "dataset.json").read_text()),
                              ["1a", "3a", "1b", "3b", "7e", "3aa"], context=180)
(RUN / "raw_model/source_mentions_pages_0001_0006.json").write_text(
    json.dumps(mentions, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/extracted_bundle_pages_0001_0006.json").write_text(
    json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
                  "source_context": len(bundle["source_context"]),
                  "stereo_diagnostics": len(stereo_diagnostics)}, ensure_ascii=False))
