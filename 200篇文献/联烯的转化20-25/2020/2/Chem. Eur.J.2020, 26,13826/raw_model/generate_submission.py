#!/usr/bin/env python3
"""Build the source-derived v3.1 bundle for Chem. Eur. J. 2020, 26, 13826.

This is an extraction artifact, not a modification of the shared tools/schema.
All machine-created records remain needs_review after standard merge.
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
DATA = json.loads((RUN / "dataset.json").read_text())
DOC = "doc_9e0daaaad7571125"


def ev(page):
    return [f"{DOC}:page{page}:text", f"{DOC}:page{page}:image"]


def evs(*pages):
    out = []
    for p in pages:
        out.extend(ev(p))
    return out


def f(raw=None, normalized=None, evidence=(), method="original_extraction", derivation=None,
      status=None, verification=None, alternatives=()):
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation,
                   status=status, verification=verification, alternatives=alternatives)


def na(evidence=()):
    return f(evidence=evidence, method="not_applicable", status="not_applicable", verification="missing")


bundle = st.empty(DATA["documents"])
molecules = {}
reactions = []


def _add_occurrence(m, label, evidence):
    occ = m["occurrences"][0]
    occ["evidence_ids"] = list(dict.fromkeys(occ["evidence_ids"] + list(evidence)))
    occ["compound_label"]["evidence_ids"] = list(dict.fromkeys(occ["compound_label"]["evidence_ids"] + list(evidence)))


def add_mol(mid, label, name=None, smiles=None, evidence=None, chirality=None,
            geometry=None, stereo_note=None):
    evidence = list(evidence or ev(3))
    if mid in molecules:
        _add_occurrence(molecules[mid], label, evidence)
        return mid
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    source_name = name or label
    m["names"]["original"] = f(source_name, source_name, evidence)
    m["names"]["english"] = f(source_name, source_name, evidence)
    m["names"]["abbreviations"] = f()
    m["names"]["aliases"] = f()
    m["cas_numbers"] = f()
    m["external_properties"] = []
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = f({"source": "source structure/name", "candidate_smiles": smiles}, smiles,
                         evidence, method="structure_derivation",
                         derivation="Candidate graph transcribed from the source drawing/name; RDKit parsing does not verify identity.")
        diag = st.stereo_inspect(smiles)
        centers = diag["atom_centers"]
        potential_db = [x for x in diag["potential_stereo"] if "Bond_Double" in x["type"]]
        if chirality is None:
            chirality = "unspecified" if any(x["CIP"] == "?" for x in centers) else "specified" if centers else "achiral"
        if geometry is None:
            if any(x["specified"] == "Specified" for x in potential_db):
                geometry = "specified"
            elif potential_db:
                geometry = "unspecified"
            else:
                geometry = "not_applicable"
        rep = "partial" if chirality in ("unspecified", "partially_specified", "unresolved") or geometry in ("unspecified", "partially_specified", "unresolved") else "encoded"
        m["stereochemistry"]["atom_centers"] = f(None, centers, evidence, method="structure_derivation", derivation="RDKit graph diagnostic on candidate SMILES; source assignment remains machine-unverified.")
        m["stereochemistry"]["double_bonds"] = f(None, potential_db, evidence, method="structure_derivation", derivation="RDKit potential-stereo diagnostic on candidate SMILES.")
        m["stereochemistry"]["axial"] = f(None, [], evidence, method="structure_derivation", derivation="No axial assignment encoded; source allene substitution was reviewed separately where applicable.")
        m["stereochemistry"]["mixture"] = f()
        m["stereochemistry"]["assignment"] = f(stereo_note, stereo_note, evidence) if stereo_note else f()
        m["stereochemistry"]["chirality_status"] = f(None, chirality, evidence, method="inference", derivation="Status from source drawing plus candidate-graph stereo diagnostic.")
        m["stereochemistry"]["geometry_status"] = f(None, geometry, evidence, method="inference", derivation="Status from source drawing plus candidate-graph stereo diagnostic.")
        m["stereochemistry"]["representation_status"] = f(None, rep, evidence, method="inference", derivation="Encoded only when all relevant stereo in the candidate graph is specified or not applicable.")
        m["stereochemistry"]["isomeric_smiles"] = f(None, diag["canonical_isomeric_smiles"], evidence, method="structure_derivation", derivation="RDKit canonical isomeric SMILES from the candidate graph.")
        m["stereochemistry"]["identity_key"] = f(None, diag["stereo_graph_key"], evidence, method="structure_derivation", derivation="SHA-256 stereo graph key from the candidate isomeric SMILES; not a sample identity.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = f()
        for k in ("atom_centers", "double_bonds", "axial", "mixture", "assignment", "isomeric_smiles", "identity_key"):
            m["stereochemistry"][k] = f()
        if stereo_note:
            m["stereochemistry"]["assignment"] = f(stereo_note, stereo_note, evidence)
        m["stereochemistry"]["chirality_status"] = f(None, chirality or "unresolved", evidence, method="inference", derivation="Source drawing was inspected; the machine record does not encode a complete stereochemical graph.")
        m["stereochemistry"]["geometry_status"] = f(None, geometry or "unresolved", evidence, method="inference", derivation="Source drawing was inspected; unspecified is not treated as achiral or racemic.")
        m["stereochemistry"]["representation_status"] = f(None, "not_encoded", evidence, method="inference", derivation="No complete source-verified SMILES was asserted.")
    m["identity_note"] = f("Machine-extracted identity; structure and stereochemistry require expert review.",
                            "Machine-extracted identity; structure and stereochemistry require expert review.", evidence)
    molecules[mid] = m
    return mid


def add_template(mid, label, fixed, variables, evidence):
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    m["names"]["original"] = f(label, label, evidence)
    m["names"]["english"] = f(label, label, evidence)
    m["names"]["abbreviations"] = f()
    m["names"]["aliases"] = f()
    m["cas_numbers"] = f()
    m["representation_type"] = "generic_template"
    m["smiles"] = f()
    m["generic_template"] = {
        "fixed_structure": f(fixed, fixed, evidence),
        "variables": [{"label": v, "attachment_points": f("source drawing", "source drawing", evidence),
                       "definition": f("R-group variable", "R-group variable", evidence),
                       "allowed_values": f(), "constraints": f()} for v in variables],
        "instance_mapping": [],
        "reason": f("Generic R-group structure in Scheme 1/Table 2 header.", "Generic R-group structure in Scheme 1/Table 2 header.", evidence),
    }
    for k in ("atom_centers", "double_bonds", "axial", "mixture", "assignment", "isomeric_smiles", "identity_key"):
        m["stereochemistry"][k] = f()
    m["stereochemistry"]["chirality_status"] = f(None, "unresolved", evidence, method="inference", derivation="Variable substituents prevent a single chirality assignment.")
    m["stereochemistry"]["geometry_status"] = f(None, "unresolved", evidence, method="inference", derivation="Generic scheme does not define a single molecular identity.")
    m["stereochemistry"]["representation_status"] = f(None, "not_encoded", evidence, method="inference", derivation="Generic template is not a concrete stereochemical SMILES.")
    m["identity_note"] = f("Generic source template; not a concrete compound.", "Generic source template; not a concrete compound.", evidence)
    m["external_properties"] = []
    molecules[mid] = m
    return mid


def ester_product(acid_smiles):
    """Create a candidate 3-series ester from candidate acid + source-derived E alcohol graph."""
    from rdkit import Chem
    alcohol = Chem.MolFromSmiles("Cc1ccc(S(=O)(=O)N(Cc2ccccc2)/C=C/CO)cc1")
    acid = Chem.MolFromSmiles(acid_smiles)
    if alcohol is None or acid is None:
        return None
    matches = acid.GetSubstructMatches(Chem.MolFromSmarts("[C:1](=[O:2])[O;H1:3]"))
    if len(matches) != 1:
        return None
    c_idx, _, acid_o = matches[0]
    combo = Chem.CombineMols(acid, alcohol)
    rw = Chem.RWMol(combo)
    offset = acid.GetNumAtoms()
    alcohol_o = next(a.GetIdx() for a in alcohol.GetAtoms()
                     if a.GetSymbol() == "O" and a.GetDegree() == 1 and a.GetTotalNumHs() == 1)
    rw.AddBond(c_idx, offset + alcohol_o, Chem.BondType.SINGLE)
    rw.RemoveAtom(acid_o)
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def set_fact(obj, key, raw, norm, evidence, method="original_extraction", derivation=None):
    obj[key] = f(raw, norm, evidence, method=method, derivation=derivation)


def participant(pid, mid, role, evidence, stage="reaction", amount=None, mass=None,
                volume=None, equiv=None, mol_percent=None, ref=None, product_role=None):
    p = st.new_participant(pid, mid, role, evidence, stage)
    if product_role is not None:
        p["product_role"] = f(product_role, product_role, evidence)
    if amount is not None:
        p["charges"]["amount"] = f(str(amount), amount, evidence)
    if mass is not None:
        p["charges"]["mass"] = f(str(mass), mass, evidence)
    if volume is not None:
        p["charges"]["volume"] = f(str(volume), volume, evidence)
    if equiv is not None:
        p["charges"]["equivalents"] = f(str(equiv), equiv, evidence)
    if mol_percent is not None:
        p["charges"]["mol_percent"] = f(str(mol_percent), mol_percent, evidence)
    if ref is not None:
        p["equivalence_reference"] = f(ref, ref, evidence)
    return p


def base_reaction(suffix, evidence, record_type="specific_experiment", categories=None,
                  location=None, general=True):
    rid = f"{DOC}:{suffix}"
    r = st.new_reaction(rid, DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories or ["other"], categories or ["other"], evidence)
    r["location"] = f(location or suffix, location or suffix, evidence)
    r["general_scheme_ids"] = [f"{DOC}:generic_hydrocarboxylation"] if general else []
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    return r


def set_hydro_class(r, evidence):
    changes = [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
        {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
    ]
    r["classification"]["bond_changes"] = [f(None, x, evidence, method="structure_derivation", derivation="Net bond change inferred from the source allene-to-(E)-enamide transformation drawing; atom mapping remains unresolved.") for x in changes]
    labels = [{"category_number": 3, "code": "C-O", "supporting_change_ids": ["b1"]},
              {"category_number": 1, "code": "C-H", "supporting_change_ids": ["b3"]}]
    r["classification"]["labels"] = f(None, labels, evidence, method="structure_derivation", derivation="Taxonomy labels derived from the explicit net changes; the C=C to C-C decrease is not a newly formed C-C sigma bond.")
    r["classification"]["coverage"] = f(None, "partial", evidence, method="inference", derivation="No complete atom-mapped reaction or explicit proton-source atom mapping is asserted.")
    r["classification"]["limitations"] = ["No complete atom mapping; hydrogen source is not atom-mapped.", "Source structures remain machine-extracted and need expert review."]
    r["classification"]["model_votes"] = []


def set_changes(r, evidence, specs):
    changes = []
    labels = {}
    for i, (a, b, before, after, kind) in enumerate(specs, 1):
        cid = f"b{i}"
        change = {"change_id": cid, "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
                  "before_order": before, "after_order": after, "change_type": kind}
        changes.append(change)
        if kind in ("formed", "hydrogen_gain", "order_increased"):
            lab = st.bond_label(a, b, after)
            labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(cid)
    r["classification"]["bond_changes"] = [f(None, x, evidence, method="structure_derivation", derivation="Net change inferred from the source transformation drawing; no complete atom mapping asserted.") for x in changes]
    r["classification"]["labels"] = f(None, list(labels.values()), evidence, method="structure_derivation", derivation="Labels are all and only the inferred formed/promoted bonds.")
    r["classification"]["coverage"] = f(None, "partial", evidence, method="inference", derivation="Graph-level net changes are partial because full atom mapping/coproduct accounting is unavailable.")
    r["classification"]["limitations"] = ["No complete atom mapping or coproduct accounting."]
    r["classification"]["model_votes"] = []


def no_class(r, evidence, coverage="unresolved"):
    r["classification"]["bond_changes"] = []
    r["classification"]["labels"] = f(None, [] if coverage == "not_applicable" else None, evidence,
                                            method="inference", derivation="No observed, structurally resolved transformation supports a bond-formation label.")
    r["classification"]["coverage"] = f(None, coverage, evidence, method="inference", derivation="Failed or unresolved experiment is not classified by its intended product.")
    r["classification"]["limitations"] = ["Do not classify a failed experiment from the intended product."]
    r["classification"]["model_votes"] = []


def step(sid, evidence, stage="reaction", pids=(), operation=None, temp=None, time=None):
    s = st.new_step(sid, evidence, stage)
    s["participant_ids"] = list(pids)
    if operation is not None:
        s["operation"] = f(operation, operation, evidence)
    if temp is not None:
        s["temperature"] = f(str(temp), temp, evidence)
    if time is not None:
        s["time"] = f(str(time), time, evidence)
    return s


def outcome(oid, kind, targets, evidence, raw=None, norm=None, yield_type=None, basis=None, scope=None):
    o = st.new_outcome(oid, kind, targets, evidence)
    if raw is not None or norm is not None:
        o["value"] = f(raw, norm, evidence)
    if yield_type is not None:
        o["yield_type"] = f(yield_type, yield_type, evidence)
    if basis is not None:
        o["measurement_basis"] = f(basis, basis, evidence)
    if scope is not None:
        o["scope"] = f(scope, scope, evidence)
    return o


# Generic R-group structures and general transformation.
add_template(f"{DOC}:mol:generic_allenamide", "generic allenamide 1", "R1R2N-C=C=CH2", ["R1", "R2"], evs(2, 3))
add_template(f"{DOC}:mol:generic_acid", "generic carboxylic acid 2", "R3-C(=O)OH", ["R3"], evs(2, 3))
add_template(f"{DOC}:mol:generic_product", "generic gamma-acyloxy (E)-enamide 3", "R1R2N-CH=CH-CH2-O-C(=O)-R3", ["R1", "R2", "R3"], evs(2, 3))
g = base_reaction("generic_hydrocarboxylation", evs(2, 3, 5, 6), "generic_scheme", ["other"], "Scheme 1B, Table 2 header, representative procedure", general=False)
g["reaction_name"] = f("metal- and reagent-free hydrocarboxylation of allenamides", "metal- and reagent-free hydrocarboxylation of allenamides", evs(1, 2))
g["transformation_type"] = f(None, "O-aminoallylation / hydrocarboxylation", evs(1, 2), method="inference", derivation="Normalized from the source title and scheme caption.")
g["participants"] = [participant("generic_allenamide", f"{DOC}:mol:generic_allenamide", "reactant", evs(2, 3)),
                     participant("generic_acid", f"{DOC}:mol:generic_acid", "reactant", evs(2, 3)),
                     participant("generic_product", f"{DOC}:mol:generic_product", "product", evs(2, 3), product_role="desired generic product")]
g["solvent_conditions"]["composition"] = f("DCE or DCE:ACN = 3:1", ["DCE", "DCE:ACN"], ev(3))
g["solvent_conditions"]["mixture_ratio"] = f("DCE:ACN = 3:1 when footnote d applies", "3:1", ev(3))
g["solvent_conditions"]["concentration"] = f("0.3 M", {"number": 0.3, "unit": "M"}, ev(3))
g["steps"] = [step("reaction", ev(3), operation="stir", temp={"number": 55, "unit": "°C"}, time={"minimum": 12, "maximum": 24, "unit": "h", "raw_text": "12-24 h"})]
g["procedure"]["original_text"] = f("Representative procedure spans pp. 5-6: 1 (0.3 mmol), 2 (0.39 mmol), DCE (1 mL, 0.3 M), 55 °C for 12 h; TLC; bicarbonate workup; DCM extraction; silica chromatography.",
                                             "Representative procedure spans pp. 5-6", evs(5, 6))
set_hydro_class(g, evs(2, 3))
g["review"]["unresolved"] = ["Generic scheme is not a single concrete experiment; product-specific exceptions are recorded separately."]
reactions.append(g)


# Core concrete molecules and common agents.
core_ev = evs(2, 3)
add_mol(f"{DOC}:mol:1a", "1a", "N-benzyl-N-allenyl-p-toluenesulfonamide", "Cc1ccc(S(=O)(=O)N(Cc2ccccc2)C=C=C)cc1", core_ev, chirality="achiral", geometry="not_applicable", stereo_note="Terminal CH2 makes the allene axis non-stereogenic.")
molecules[f"{DOC}:mol:1a"]["stereochemistry"]["representation_status"] = f(
    None, "partial", core_ev, method="inference",
    derivation="RDKit flags potential allene stereochemistry, while the source structure has a terminal CH2; axial status is retained in the assignment note rather than claimed fully encoded.")
add_mol(f"{DOC}:mol:2a", "2a", "benzoic acid", "O=C(O)c1ccccc1", core_ev)
add_mol(f"{DOC}:mol:3a", "3a", "gamma-benzoyloxy (E)-enamide 3a", "Cc1ccc(S(=O)(=O)N(Cc2ccccc2)/C=C/COC(=O)c2ccccc2)cc1", core_ev, geometry="specified", stereo_note="Source explicitly assigns the product alkene as E.")
add_mol(f"{DOC}:mol:2l", "2l", "acetic acid", "CC(=O)O", evs(3, 6))
add_mol(f"{DOC}:mol:3l", "3l", "gamma-acetoxy (E)-enamide 3l", "CC(=O)OC/C=C/N(Cc1ccccc1)S(=O)(=O)c1ccc(C)cc1", evs(3, 5, 6), geometry="specified", stereo_note="Source Scheme 2 depicts the E-enamide.")

agent_smiles = {
    "DCE": ("1,2-dichloroethane", "ClCCCl", evs(2, 3, 5, 6)), "DCM": ("dichloromethane", "ClCCl", evs(2, 5, 6)),
    "CHCl3": ("chloroform", "ClC(Cl)Cl", ev(2)), "toluene": ("toluene", "Cc1ccccc1", ev(2)),
    "THF": ("tetrahydrofuran", "C1CCOC1", ev(2)), "dioxane": ("1,4-dioxane", "O1CCOCC1", ev(2)),
    "DME": ("1,2-dimethoxyethane", "COCCOC", ev(2)), "ACN": ("acetonitrile", "CC#N", evs(2, 3, 5)),
    "MeOH": ("methanol", "CO", evs(2, 5)), "HFIP": ("hexafluoroisopropanol", "OC(C(F)(F)F)C(F)(F)F", ev(2)),
    "iPrOH": ("isopropanol", "CC(C)O", ev(2)), "EtOAc": ("ethyl acetate", "CCOC(C)=O", ev(2)),
    "NaHCO3": ("sodium bicarbonate", "[Na+].O=C([O-])O", evs(5, 6)), "silica": ("silica gel", None, evs(5, 6)),
    "K2CO3": ("potassium carbonate", "[K+].[K+].[O-]C([O-])=O", ev(5)), "phenylacetylene": ("phenylacetylene", "C#Cc1ccccc1", ev(5)),
    "ZnBr2": ("zinc bromide", "Br[Zn]Br", ev(5)), "DIPEA": ("diisopropylethylamine", "CCN(C(C)C)C(C)C", ev(5)),
    "BF3OEt2": ("boron trifluoride diethyl etherate", "CCOCC.F[B](F)F", ev(5)), "2phenylindole": ("2-phenylindole", "c1ccc(-c2[nH]c3ccccc3c2)cc1", ev(5)),
    "NIS": ("N-iodosuccinimide", "O=C1CCC(=O)N1I", ev(5)), "NBS": ("N-bromosuccinimide", "O=C1CCC(=O)N1Br", ev(5)),
    "mCPBA": ("m-chloroperbenzoic acid", "OOC(=O)c1cccc(Cl)c1", ev(5)),
    "CuBase": ("Cu/base catalyst system reported in ref. 9g", None, ev(5)),
    "dibromomethane": ("1,2-dibromomethane (as printed)", None, ev(2)),
}
for key, (name, smi, evidence) in agent_smiles.items():
    add_mol(f"{DOC}:mol:{key}", key, name, smi, evidence)


# Table 1: every optimization row is a distinct experiment.
opt_rows = [
    (1, "DCM", 25, 24, 74), (2, "DCM", 55, 12, 81), (3, "DCE", 55, 12, 97),
    (4, "CHCl3", 55, 12, 93), (5, "toluene", 55, 16, 71), (6, "THF", 55, 24, 36),
    (7, "dioxane", 55, 24, 72), (8, "DME", 55, 24, 57), (9, "ACN", 55, 24, 96),
    (10, "MeOH", 55, 24, 38), (11, "HFIP", 55, 1, "trace"), (12, "iPrOH", 55, 24, 41),
    (13, "EtOAc", 55, 48, 67),
]
recoveries = {1: 23, 6: 46, 7: 11, 8: 31, 13: 31}
for n, solv, tempv, timev, yld in opt_rows:
    e = ev(2)
    r = base_reaction(f"table1_entry_{n:02d}", e, categories=["optimization"], location=f"Table 1 entry {n}")
    r["reaction_name"] = f("metal-free hydrocarboxylation", "metal-free hydrocarboxylation", e)
    r["transformation_type"] = f(None, "allenamide hydrocarboxylation", e, method="inference", derivation="Normalized from Table 1 title.")
    r["participants"] = [participant("p_1a", f"{DOC}:mol:1a", "reactant", e, amount={"number": 0.1, "unit": "mmol"}),
                         participant("p_2a", f"{DOC}:mol:2a", "reactant", e, amount={"number": 0.13, "unit": "mmol"}, equiv={"number": 1.3, "unit": "equiv"}, ref=f"{DOC}:mol:1a"),
                         participant("p_solvent", f"{DOC}:mol:{solv}", "solvent", e),
                         participant("p_3a", f"{DOC}:mol:3a", "product", e, product_role="desired product")]
    r["reaction_scale"] = f("0.1 mmol of 1a", {"reference_molecule_id": f"{DOC}:mol:1a", "amount": {"number": 0.1, "unit": "mmol"}}, e)
    r["solvent_conditions"]["composition"] = f(solv, [solv], e)
    r["solvent_conditions"]["concentration"] = f("0.3 M", {"number": 0.3, "unit": "M"}, e)
    r["solvent_conditions"]["concentration_basis"] = f("1a", f"{DOC}:mol:1a", e)
    r["steps"] = [step("reaction", e, pids=["p_1a", "p_2a", "p_solvent", "p_3a"], operation="stir", temp={"number": tempv, "unit": "°C"}, time={"number": timev, "unit": "h"})]
    r["procedure"]["original_text"] = f(f"Table 1 entry {n}: {solv}, {tempv} °C, {timev} h.", f"Table 1 entry {n}", e)
    if yld == "trace":
        r["outcomes"] = [outcome("desired_trace", "trace", [f"{DOC}:mol:3a"], e, "trace", "trace", basis="internal-standard yield measurement; technique not stated")]
    else:
        r["outcomes"] = [outcome("yield_internal_standard", "yield", [f"{DOC}:mol:3a"], e, f"{yld}%", {"number": yld, "unit": "%"}, yield_type="reported yield", basis="measured using 1,2-dibromomethane as internal standard")]
    if n == 3:
        r["outcomes"].append(outcome("yield_isolated", "yield", [f"{DOC}:mol:3a"], e, "(96%)", {"number": 96, "unit": "%"}, yield_type="isolated yield", basis="after column chromatography"))
    if n in recoveries:
        r["outcomes"].append(outcome("recovery_1a", "recovery", [f"{DOC}:mol:1a"], e, f"{recoveries[n]}% of 1a remained unreacted", {"number": recoveries[n], "unit": "%"}, basis="reported as unreacted 1a; method not stated"))
    if n == 6:
        r["outcomes"].append(outcome("unidentified_products", "other", [], e, "unidentifiable compounds formed", "unidentifiable compounds formed"))
    if n == 11:
        r["outcomes"].append(outcome("decomposition_1a", "decomposition", [f"{DOC}:mol:1a"], e, "Complete decomposition of 1a was observed", "complete decomposition"))
    set_hydro_class(r, e)
    r["review"]["unresolved"] = ["The internal-standard technique is not specified on the supplied pages."]
    reactions.append(r)


# Table 2 carboxylic-acid scope. Known candidate SMILES are used only where the drawing/name is clear.
acid_smiles = {
    "2a": "O=C(O)c1ccccc1", "2b": "O=C(O)c1c(C)cccc1C", "2c": "COc1ccc(C(=O)O)cc1",
    "2d": "O=C(O)c1ccccc1Cl", "2e": "O=C(O)c1ccccc1I", "2f": "O=C(O)c1cccc([N+](=O)[O-])c1",
    "2g": "O=C(O)c1ccccc1-c1ccccc1", "2i": "O=C(O)Cc1ccccc1[N+](=O)[O-]",
    "2j": "CC(C(=O)O)c1ccc(Cl)cc1", "2k": "CC(C(=O)O)c1ccc([N+](=O)[O-])cc1", "2l": "CC(=O)O",
    "2m": "O=C(O)C1CCCCC1", "2n": "CC(Oc1ccccc1)C(=O)O", "2o": "CC(C)(C)C(=O)O",
    "2p": "C=CC(=O)O", "2v": "O=C(O)c1ccccc1O", "2w": "O=C(O)CNC(=O)c1ccccc1",
    "2aa": "O=C(O)c1cc2ccccc2s1", "2ab": "O=C(O)c1cc2ccccc2o1", "2ag": "O=C(O)CCCCC1CCSS1",
    "2ah": "CC(=O)Oc1ccccc1C(=O)O", "2aj": "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
}
scope_rows = [
    ("3a", "2a", 96, "benzoic acid"), ("3b", "2b", 95, "2,6-dimethylbenzoic acid"),
    ("3c", "2c", 90, "4-methoxybenzoic acid"), ("3d", "2d", 92, "2-chlorobenzoic acid"),
    ("3e", "2e", 93, "2-iodobenzoic acid"), ("3f", "2f", 96, "nitrobenzoic acid 2f (position from drawing)"),
    ("3g", "2g", 93, "2-phenylbenzoic acid"), ("3h", "2h", 96, "bromo-methylbenzoic acid 2h (structure shown)"),
    ("3i", "2i", 91, "2-(2-nitrophenyl)acetic acid"), ("3j", "2j", 90, "2-(4-chlorophenyl)propionic acid; configuration not named"),
    ("3k", "2k", 90, "2-(4-nitrophenyl)propionic acid; configuration not named"), ("3l", "2l", 95, "acetic acid"),
    ("3m", "2m", 91, "cyclohexanecarboxylic acid"), ("3n", "2n", 87, "2-phenoxypropionic acid; configuration not named"),
    ("3o", "2o", 89, "pivalic acid"), ("3p", "2p", 88, "acrylic acid"),
    ("3q", "2q", 90, "norbornene carboxylic acid 2q (structure shown)"), ("3r", "2r", 90, "10-undecenoic acid"),
    ("3s", "2s", 91, "linoleic acid"), ("3t", "2t", 93, "bromoalkanoic acid 2t (chain length from drawing unresolved)"),
    ("3u", "2u", 86, "hydroxy acid 2u (structure shown)"), ("3v", "2v", 94, "salicylic acid"),
    ("3w", "2w", 82, "hippuric acid"), ("3x", "2x", 86, "N-acetyl-L-phenylalanine"),
    ("3y", "2y", 85, "N-Boc-L-leucine"), ("3z", "2z", 88, "N-Boc-L-proline"),
    ("3aa", "2aa", 84, "benzothiophene carboxylic acid 2aa"), ("3ab", "2ab", 92, "benzofuran carboxylic acid 2ab"),
    ("3ac", "2ac", 81, "quinoline carboxylic acid 2ac"), ("3ad", "2ad", 89, "chloro fused N-heterocyclic carboxylic acid 2ad"),
    ("3ae", "2ae", 90, "fluorophenyl methyl isoxazole carboxylic acid 2ae"), ("3af", "2af", 84, "indole-tethered carboxylic acid 2af"),
    ("3ag", "2ag", 86, "alpha-lipoic acid"), ("3ah", "2ah", 93, "aspirin"),
    ("3ai", "2ai", 92, "naproxen"), ("3aj", "2aj", 91, "ibuprofen"),
    ("3ak", "2ak", 71, "mycophenolic acid"), ("3al", "2al", 84, "indomethacin"),
    ("3am", "2am", 71, "chenodeoxycholic acid"), ("3an", "2an", 67, "enoxolone"),
]
chiral_acids = {"2j", "2k", "2n", "2x", "2y", "2z", "2ai", "2aj", "2am", "2an"}
poor_sol = {"3ad", "3ae", "3af", "3ag", "3ak", "3am", "3an"}
fast = {"3v", "3ab"}
for prod, acid, yld, acid_name in scope_rows:
    e = evs(3, 4) if prod in {"3aa", "3ab", "3ac", "3ad", "3ae", "3af", "3ag", "3ah", "3ai", "3aj", "3ak", "3al", "3am", "3an"} else ev(3)
    amid = f"{DOC}:mol:{acid}"
    if acid not in ("2a", "2l"):
        add_mol(amid, acid, acid_name, acid_smiles.get(acid), e,
                chirality="partially_specified" if acid in chiral_acids else None,
                stereo_note="Source drawing/name indicates stereochemical content; absolute/mixture assignment is not fully encoded." if acid in chiral_acids else None)
    pmid = f"{DOC}:mol:{prod}"
    if prod not in ("3a", "3l"):
        psmi = ester_product(acid_smiles[acid]) if acid in acid_smiles else None
        add_mol(pmid, prod, f"gamma-acyloxy (E)-enamide {prod}", psmi, e,
                chirality="partially_specified" if acid in chiral_acids else None,
                geometry="specified", stereo_note="New enamide alkene is assigned E; inherited stereocenters, where present, require expert review.")
    r = base_reaction(f"table2_{prod}", e, categories=["substrate_scope"], location=f"Table 2A, {prod}")
    r["reaction_name"] = f("hydrocarboxylation", "hydrocarboxylation", e)
    r["transformation_type"] = f(None, "linear-selective O-aminoallylation to gamma-acyloxy (E)-enamide", e, method="inference", derivation="Normalized from the general scheme and product class.")
    solvent_mid = f"{DOC}:mol:ACN" if False else f"{DOC}:mol:DCE"
    r["participants"] = [participant("p_1a", f"{DOC}:mol:1a", "reactant", e, amount={"number": 0.3, "unit": "mmol"}),
                         participant("p_acid", amid, "reactant", e, amount={"number": 0.39, "unit": "mmol"}, equiv={"number": 1.3, "unit": "equiv"}, ref=f"{DOC}:mol:1a"),
                         participant("p_DCE", f"{DOC}:mol:DCE", "solvent", e),
                         participant("p_product", pmid, "product", e, product_role="desired E-regioisomer")]
    if prod in poor_sol:
        r["participants"].append(participant("p_ACN", f"{DOC}:mol:ACN", "solvent", e))
    r["reaction_scale"] = f("0.3 mmol of 1", {"reference_molecule_id": f"{DOC}:mol:1a", "amount": {"number": 0.3, "unit": "mmol"}}, e)
    comp = ["DCE", "ACN"] if prod in poor_sol else ["DCE"]
    r["solvent_conditions"]["composition"] = f("DCE:ACN = 3:1" if prod in poor_sol else "DCE", comp, e)
    if prod in poor_sol:
        r["solvent_conditions"]["mixture_ratio"] = f("3:1 DCE:ACN", "3:1", e)
    r["solvent_conditions"]["total_volume"] = f("1 mL", {"number": 1, "unit": "mL"}, evs(5, 6), method="inherited_general_condition")
    r["solvent_conditions"]["concentration"] = f("0.3 M", {"number": 0.3, "unit": "M"}, e)
    r["solvent_conditions"]["concentration_basis"] = f("allenamide 1", f"{DOC}:mol:1a", e)
    t = 30 if prod in fast else 24 if prod == "3an" else 12
    tu = "min" if prod in fast else "h"
    r["steps"] = [step("reaction", e, pids=[p["participant_id"] for p in r["participants"]], operation="stir", temp={"number": 55, "unit": "°C"}, time={"number": t, "unit": tu})]
    r["procedure"]["general_procedure_reference"] = f("Representative procedure for hydrocarboxylation, pp. 5-6", "Representative procedure for hydrocarboxylation", evs(5, 6))
    r["outcomes"] = [outcome("isolated_yield", "yield", [pmid], e, f"{yld}%", {"number": yld, "unit": "%"}, yield_type="isolated yield", basis="isolated yield given below structure")]
    r["outcomes"].append(outcome("stereo_regioselectivity", "selectivity", [pmid], evs(1, 3, 4), ">49:1; exclusive E-selectivity / single E-regioisomer", {"qualifier": "gt", "ratio": [49, 1], "descriptor": "E and linear regioselectivity"}, basis=None, scope={"metric": "E/regioselectivity", "major": "drawn E linear product", "minor": "not identified"}))
    set_hydro_class(r, e)
    r["review"]["unresolved"] = ["Selectivity measurement method is not reported on supplied pages.", "Supporting Information was not provided; product characterization and exceptions beyond shown footnotes remain unavailable."]
    reactions.append(r)


# Table 2B: other allenamides, including the two limitations.
other_rows = [("3ao", "1b", 91), ("3ap", "1c", 93), ("3aq", "1d", 90), ("3ar", "1e", 88),
              ("3as", "1f", 90), ("3at", "1g", 92), ("3au", "1h", 94)]
for prod, allene, yld in other_rows:
    e = evs(3, 4)
    add_mol(f"{DOC}:mol:{allene}", allene, f"allenamide {allene} (structure corresponding to {prod})", None, e, chirality="unresolved", geometry="not_applicable", stereo_note="Terminal allene substrate; substituent identity is retained from the source drawing/text but not encoded.")
    add_mol(f"{DOC}:mol:{prod}", prod, f"gamma-benzoyloxy (E)-enamide {prod}", None, e, geometry="specified", stereo_note="Source reports similar selectivity; product drawing is E.")
    r = base_reaction(f"table2_{prod}", e, categories=["substrate_scope"], location=f"Table 2B, {prod}")
    r["reaction_name"] = f("hydrocarboxylation", "hydrocarboxylation", e)
    r["transformation_type"] = f(None, "allenamide hydrocarboxylation", e, method="inference", derivation="Normalized from Table 2B heading.")
    r["participants"] = [participant("p_allene", f"{DOC}:mol:{allene}", "reactant", e, amount={"number": 0.3, "unit": "mmol"}),
                         participant("p_2a", f"{DOC}:mol:2a", "reactant", e, amount={"number": 0.39, "unit": "mmol"}, equiv={"number": 1.3, "unit": "equiv"}, ref=f"{DOC}:mol:{allene}"),
                         participant("p_DCE", f"{DOC}:mol:DCE", "solvent", e), participant("p_product", f"{DOC}:mol:{prod}", "product", e, product_role="desired E-regioisomer")]
    r["reaction_scale"] = f("0.3 mmol of 1", {"reference_molecule_id": f"{DOC}:mol:{allene}", "amount": {"number": 0.3, "unit": "mmol"}}, e)
    r["solvent_conditions"]["composition"] = f("DCE", ["DCE"], e)
    r["solvent_conditions"]["concentration"] = f("0.3 M", {"number": 0.3, "unit": "M"}, e)
    r["steps"] = [step("reaction", e, pids=["p_allene", "p_2a", "p_DCE", "p_product"], operation="stir", temp={"number": 55, "unit": "°C"}, time={"number": 12, "unit": "h"})]
    r["outcomes"] = [outcome("isolated_yield", "yield", [f"{DOC}:mol:{prod}"], e, f"{yld}%", {"number": yld, "unit": "%"}, yield_type="isolated yield", basis="isolated yield given below structure")]
    r["procedure"]["general_procedure_reference"] = f("Table 2 optimal conditions; benzoic acid 2a", "Table 2 footnote f", e)
    set_hydro_class(r, e)
    r["review"]["unresolved"] = ["Complete concrete substrate/product SMILES require expert structure transcription."]
    reactions.append(r)

for suffix, name, finding in [("3av", "allenyl phosphine oxide associated with 3av", "allenyl phosphine oxide remained intact"),
                              ("3aw", "OBn allenyl ether associated with 3aw", "OBn ether gave an unidentifiable mixture of products")]:
    e = evs(3, 4)
    mid = f"{DOC}:mol:input_{suffix}"
    add_mol(mid, f"input for {suffix}", name, None, e, geometry="not_applicable")
    r = base_reaction(f"table2_{suffix}_limitation", e, categories=["substrate_scope", "control"], location=f"Table 2B limitation, {suffix}")
    r["reaction_name"] = f("attempted hydrocarboxylation", "attempted hydrocarboxylation", e)
    r["participants"] = [participant("p_allene", mid, "reactant", e, amount={"number": 0.3, "unit": "mmol"}),
                         participant("p_2a", f"{DOC}:mol:2a", "reactant", e, amount={"number": 0.39, "unit": "mmol"}, equiv={"number": 1.3, "unit": "equiv"}, ref=mid),
                         participant("p_DCE", f"{DOC}:mol:DCE", "solvent", e)]
    r["reaction_scale"] = f("0.3 mmol of allene derivative", {"reference_molecule_id": mid, "amount": {"number": 0.3, "unit": "mmol"}}, e)
    r["solvent_conditions"]["composition"] = f("DCE", ["DCE"], e)
    r["solvent_conditions"]["concentration"] = f("0.3 M", {"number": 0.3, "unit": "M"}, e)
    r["steps"] = [step("reaction", e, pids=["p_allene", "p_2a", "p_DCE"], operation="stir under Table 2 conditions", temp={"number": 55, "unit": "°C"}, time={"number": 12, "unit": "h"})]
    kind = "no_reaction" if suffix == "3av" else "other"
    norm = True if suffix == "3av" else "unidentifiable mixture"
    r["outcomes"] = [outcome("limitation_result", kind, [mid] if suffix == "3av" else [], e, finding, norm, basis="visual/source statement; analytical method not reported")]
    no_class(r, e, "not_applicable" if suffix == "3av" else "unresolved")
    r["review"]["unresolved"] = ["No desired product structure/yield was reported; unidentified mixture composition is unavailable." if suffix == "3aw" else "Extent/method for 'remained intact' is not reported."]
    reactions.append(r)


# Table 3: paired experiments for entries 1-7; entries 8-11 are cited reported yields.
table3 = [(1, "3h", "2h", 97, 66), (2, "3o", "2o", 90, 0), (3, "3p", "2p", 92, 32),
          (4, "3af", "2af", 87, 32), (5, "3ag", "2ag", 88, "complex mixture"),
          (6, "3ak", "2ak", 74, "<5"), (7, "3an", "2an", 70, 0)]
add_mol(f"{DOC}:mol:NHTsBn", "NHTs(Bn)", "N-benzyl-p-toluenesulfonamide", "Cc1ccc(S(=O)(=O)NCc2ccccc2)cc1", ev(5))
for n, prod, acid, ya, yb in table3:
    e = evs(4, 5)
    for cond, val in (("A", ya), ("B", yb)):
        r = base_reaction(f"table3_entry_{n:02d}_{cond}", e, categories=["control"], location=f"Table 3 entry {n}, condition {cond}")
        r["reaction_name"] = f("metal-free hydrocarboxylation" if cond == "A" else "Cu/base-catalyzed hydrocarboxylation comparison", "metal-free hydrocarboxylation" if cond == "A" else "Cu/base-catalyzed hydrocarboxylation comparison", e)
        r["transformation_type"] = f(None, "allenamide hydrocarboxylation comparison", e, method="inference", derivation="Normalized from Table 3 heading.")
        r["participants"] = [participant("p_1a", f"{DOC}:mol:1a", "reactant", e, amount={"number": 0.1, "unit": "mmol"}),
                             participant("p_acid", f"{DOC}:mol:{acid}", "reactant", e),
                             participant("p_intended", f"{DOC}:mol:{prod}", "product", e, product_role="desired product (may be undetected)")]
        if cond == "A":
            r["participants"].append(participant("p_DCE", f"{DOC}:mol:DCE", "solvent", e))
            r["solvent_conditions"]["composition"] = f("similar conditions to Table 2", ["DCE"], e, method="inherited_general_condition")
            r["solvent_conditions"]["concentration"] = f("0.3 M (Table 2 conditions)", {"number": 0.3, "unit": "M"}, e, method="inherited_general_condition")
            hold = 24 if prod == "3an" else 12
            r["steps"] = [step("reaction", e, pids=[p["participant_id"] for p in r["participants"]], operation="condition A: similar to Table 2", temp={"number": 55, "unit": "°C"}, time={"number": hold, "unit": "h"})]
            r["outcomes"] = [outcome("nmr_yield", "yield", [f"{DOC}:mol:{prod}"], e, f"{val}%", {"number": val, "unit": "%"}, yield_type="1H NMR yield", basis="crude 1H NMR analysis using an internal standard")]
            set_hydro_class(r, e)
        else:
            r["participants"].append(participant("p_CuBase", f"{DOC}:mol:CuBase", "catalyst", e))
            t = 36 if n in (2, 7) else None
            r["steps"] = [step("reaction", e, pids=[p["participant_id"] for p in r["participants"]], operation="condition B: Cu-/base-catalyzed conditions reported in ref. 9g were employed", time={"number": t, "unit": "h"} if t else None)]
            if isinstance(val, int) and val > 0:
                r["outcomes"] = [outcome("nmr_yield", "yield", [f"{DOC}:mol:{prod}"], e, f"{val}%", {"number": val, "unit": "%"}, yield_type="1H NMR yield", basis="crude 1H NMR analysis using an internal standard")]
                set_hydro_class(r, e)
            elif val == "<5":
                r["outcomes"] = [outcome("nmr_yield", "yield", [f"{DOC}:mol:{prod}"], e, "<5%", {"number": 5, "unit": "%", "qualifier": "lt", "raw_text": "<5%"}, yield_type="1H NMR yield", basis="crude 1H NMR analysis using an internal standard")]
                set_hydro_class(r, e)
            elif val == 0:
                side = 91 if n == 2 else 93
                r["participants"].append(participant("p_NHTsBn", f"{DOC}:mol:NHTsBn", "product", e, product_role="observed decomposition/cleavage product"))
                r["outcomes"] = [outcome("desired_not_detected", "no_reaction", [f"{DOC}:mol:{prod}"], e, "Even trace amounts of the desired product were not detected", True, basis="analytical detection limit not reported"),
                                     outcome("NHTsBn_yield", "yield", [f"{DOC}:mol:NHTsBn"], e, f"{side}%", {"number": side, "unit": "%"}, yield_type="1H NMR yield", basis="crude 1H NMR analysis using an internal standard")]
                no_class(r, e, "unresolved")
            else:
                r["outcomes"] = [outcome("mixture", "other", [], e, "complex mixture", "complex mixture", basis="crude 1H NMR analysis; individual products not assigned")]
                no_class(r, e, "unresolved")
        r["reaction_scale"] = f("0.1 mmol of 1a", {"reference_molecule_id": f"{DOC}:mol:1a", "amount": {"number": 0.1, "unit": "mmol"}}, e)
        r["review"]["unresolved"] = ["Condition B reagent loadings, solvent and temperature are only referenced to ref. 9g and are not provided on supplied pages."] if cond == "B" else ["Acid charge and internal-standard identity are not reported in Table 3; condition A only says similar to Table 2."]
        reactions.append(r)

for n, prod, allene, yld in [(8, "3ao", "1b", 43), (9, "3aq", "1d", 47), (10, "3at", "1g", 65), (11, "3au", "1h", 80)]:
    e = ev(5)
    r = base_reaction(f"table3_entry_{n:02d}_reported_ref9g", e, record_type="literature_example", categories=["other"], location=f"Table 3 entry {n}, parenthesized reported yield", general=False)
    r["reaction_name"] = f("reported Cu-catalyzed hydrocarboxylation from ref. 9g", "reported Cu-catalyzed hydrocarboxylation from ref. 9g", e)
    r["participants"] = [participant("p_allene", f"{DOC}:mol:{allene}", "reactant", e), participant("p_2a", f"{DOC}:mol:2a", "reactant", e), participant("p_product", f"{DOC}:mol:{prod}", "product", e, product_role="reported literature product")]
    r["outcomes"] = [outcome("reported_yield", "yield", [f"{DOC}:mol:{prod}"], e, f"({yld})", {"number": yld, "unit": "%"}, yield_type="reported isolated yield", basis="reported yield from ref. 9g, reproduced in parentheses for direct comparison")]
    set_hydro_class(r, e)
    r["review"]["unresolved"] = ["This is prior literature, not a new experiment in the present paper; detailed conditions are not reproduced on the supplied pages."]
    reactions.append(r)


# Scheme 2 derivatization experiments.
add_mol(f"{DOC}:mol:4a", "4a", "gamma-hydroxy (E)-enamide 4a", "OC/C=C/N(Cc1ccccc1)S(=O)(=O)c1ccc(C)cc1", ev(5), geometry="specified", stereo_note="Source drawing depicts the E-enamide.")
for label, note in [("4b", "ynenamide 4b"), ("4c", "enamide-tethered indole 4c"),
                    ("4d", "alpha-benzoyloxy-beta-iodo-N,O-acetal 4d"), ("4e", "bromomethoxylation product 4e"),
                    ("4f", "hydroxycarboxylation product 4f"), ("4g", "iodoamination product 4g")]:
    ch = "partially_specified" if label in {"4d", "4e", "4f", "4g"} else "unresolved"
    geom = "not_applicable" if label in {"4d", "4e", "4f", "4g"} else "specified"
    add_mol(f"{DOC}:mol:{label}", label, note, None, ev(5), chirality=ch, geometry=geom,
            stereo_note="Source uses wedge/dash bonds for relative stereochemistry; absolute configuration and mixture composition are not reported/encoded." if ch == "partially_specified" else "Source product drawing retains an E-enamide; complete graph not encoded.")

deriv = [
    ("a", "4a", 91, [("K2CO3", "reagent", {"number": 3, "unit": "equiv"}), ("MeOH", "solvent", None)], "rt", {"number": 12, "unit": "h"}, "hydrolysis", [("C", "O", 1, 0, "broken"), ("O", "H", 0, 1, "hydrogen_gain")]),
    ("b", "4b", 39, [("phenylacetylene", "reactant", {"number": 1.5, "unit": "equiv"}), ("ZnBr2", "catalyst", {"number": 20, "unit": "mol%"}), ("DIPEA", "reagent", {"number": 1.5, "unit": "equiv"}), ("BF3OEt2", "reagent", {"number": 1.1, "unit": "equiv"}), ("DCM", "solvent", None)], "rt", {"number": 2, "unit": "h"}, "alkynylation", [("C", "O", 1, 0, "broken"), ("C", "C", 0, 1, "formed")]),
    ("c", "4c", 87, [("2phenylindole", "reactant", {"number": 1.1, "unit": "equiv"}), ("BF3OEt2", "reagent", {"number": 1.1, "unit": "equiv"}), ("DCM", "solvent", None)], {"number": 0, "unit": "°C"}, {"number": 5, "unit": "min"}, "indolation", [("C", "O", 1, 0, "broken"), ("C", "C", 0, 1, "formed")]),
    ("d", "4d", 89, [("NIS", "reagent", {"number": 1.1, "unit": "equiv"}), ("2a", "reactant", {"number": 1.3, "unit": "equiv"}), ("DCM", "solvent", None)], "rt", {"number": 1, "unit": "h"}, "iodocarboxylation", [("C", "C", 2, 1, "order_decreased"), ("C", "I", 0, 1, "formed"), ("C", "O", 0, 1, "formed")]),
    ("e", "4e", 94, [("NBS", "reagent", {"number": 1.1, "unit": "equiv"}), ("MeOH", "solvent", None)], "rt", {"number": 30, "unit": "min"}, "bromomethoxylation", [("C", "C", 2, 1, "order_decreased"), ("C", "Br", 0, 1, "formed"), ("C", "O", 0, 1, "formed")]),
    ("f", "4f", 75, [("mCPBA", "reagent", {"number": 1.3, "unit": "equiv"}), ("DCM", "solvent", None)], "rt", {"number": 6, "unit": "h"}, "hydroxycarboxylation", [("C", "C", 2, 1, "order_decreased"), ("C", "O", 0, 1, "formed"), ("C", "O", 0, 1, "formed")]),
    ("g", "4g", 94, [("NIS", "reagent", {"number": 1.1, "unit": "equiv"}), ("ACN", "solvent", None)], "rt", {"number": 1, "unit": "h"}, "iodoamination", [("C", "C", 2, 1, "order_decreased"), ("C", "I", 0, 1, "formed"), ("C", "N", 0, 1, "formed")]),
]
for letter, prod, yld, agents, tempv, timev, name, changes in deriv:
    e = ev(5)
    r = base_reaction(f"scheme2_{letter}_{prod}", e, categories=["derivatization"], location=f"Scheme 2, arrow {letter}", general=False)
    r["reaction_name"] = f(name, name, e)
    r["transformation_type"] = f(None, name, e, method="inference", derivation="Normalized from Scheme 2 arrow label/caption.")
    r["participants"] = [participant("p_3l", f"{DOC}:mol:3l", "reactant", e, amount={"number": 0.2, "unit": "mmol"}), participant("p_product", f"{DOC}:mol:{prod}", "product", e, product_role="desired product")]
    for key, role, dose in agents:
        mid = f"{DOC}:mol:{key}"
        pid = f"p_{key}"
        kwargs = {}
        if dose:
            if dose["unit"] == "equiv": kwargs.update(equiv=dose, ref=f"{DOC}:mol:3l")
            elif dose["unit"] == "mol%": kwargs.update(mol_percent=dose, ref=f"{DOC}:mol:3l")
        r["participants"].append(participant(pid, mid, role, e, **kwargs))
    r["reaction_scale"] = f("0.2 mmol of 3l", {"reference_molecule_id": f"{DOC}:mol:3l", "amount": {"number": 0.2, "unit": "mmol"}}, e)
    solvents = [k for k, role, _ in agents if role == "solvent"]
    if solvents:
        r["solvent_conditions"]["composition"] = f(", ".join(solvents), solvents, e)
    if letter == "e":
        r["solvent_conditions"]["concentration"] = f("MeOH (0.3 M)", {"number": 0.3, "unit": "M"}, e)
        r["solvent_conditions"]["concentration_basis"] = f("3l", f"{DOC}:mol:3l", e)
    r["steps"] = [step("reaction", e, pids=[p["participant_id"] for p in r["participants"]], operation=name, temp=tempv, time=timev)]
    r["procedure"]["original_text"] = f(f"Scheme 2({letter}) reagent list and conditions", f"Scheme 2({letter})", e)
    r["outcomes"] = [outcome("yield", "yield", [f"{DOC}:mol:{prod}"], e, f"{yld}%", {"number": yld, "unit": "%"}, yield_type="reported yield", basis=None)]
    set_changes(r, e, changes)
    r["review"]["unresolved"] = ["Yield measurement/isolation basis is not stated in the supplied Scheme 2 caption; SI was not provided.", "Absolute stereochemical assignments for new centers are not reported in the supplied pages." if prod in {"4d", "4e", "4f", "4g"} else "Complete product graph needs expert review."]
    reactions.append(r)


# Gram-scale preparation of 3l; Scheme 2 and Experimental Section disagree (96% vs 95%).
e = evs(4, 5, 6)
r = base_reaction("gram_scale_3l", e, categories=["scale_up"], location="Scheme 2 footnote a and Experimental Section gram-scale synthesis", general=True)
r["reaction_name"] = f("Gram scale synthesis of hydrocarboxylation product 3l", "Gram scale synthesis of hydrocarboxylation product 3l", e)
r["transformation_type"] = f(None, "allenamide hydrocarboxylation scale-up", e, method="inference", derivation="Normalized from the experimental heading.")
r["participants"] = [participant("p_1a", f"{DOC}:mol:1a", "reactant", e, amount={"number": 5, "unit": "mmol"}, mass={"number": 1.495, "unit": "g"}),
                     participant("p_2l", f"{DOC}:mol:2l", "reactant", e, amount={"number": 6.5, "unit": "mmol"}, volume={"number": 0.37, "unit": "mL"}, equiv={"number": 1.3, "unit": "equiv"}, ref=f"{DOC}:mol:1a"),
                     participant("p_DCE", f"{DOC}:mol:DCE", "solvent", e, volume={"number": 16.7, "unit": "mL"}),
                     participant("p_3l", f"{DOC}:mol:3l", "product", e, product_role="desired product"),
                     participant("p_NaHCO3", f"{DOC}:mol:NaHCO3", "workup_agent", e, stage="workup", volume={"number": 15, "unit": "mL"}),
                     participant("p_DCM", f"{DOC}:mol:DCM", "workup_agent", e, stage="workup", volume={"number": 40, "unit": "mL"})]
r["reaction_scale"] = f("5 mmol of 1a", {"reference_molecule_id": f"{DOC}:mol:1a", "amount": {"number": 5, "unit": "mmol"}}, e)
r["solvent_conditions"]["composition"] = f("dichloroethane", ["DCE"], ev(6))
r["solvent_conditions"]["total_volume"] = f("16.7 mL", {"number": 16.7, "unit": "mL"}, ev(6))
r["solvent_conditions"]["concentration"] = f("0.4 M", {"number": 0.4, "unit": "M"}, ev(6))
r["solvent_conditions"]["concentration_basis"] = f("1a", f"{DOC}:mol:1a", ev(6))
r["steps"] = [
    step("addition", ev(6), pids=["p_1a", "p_2l", "p_DCE"], operation="add acetic acid 2l to 1a in DCE at room temperature", temp="room temperature"),
    step("reaction", ev(6), pids=["p_1a", "p_2l", "p_DCE", "p_3l"], operation="heat", temp={"number": 55, "unit": "°C"}, time={"number": 12, "unit": "h"}),
    step("cool", ev(6), pids=["p_3l"], operation="cool to room temperature", temp="room temperature"),
    step("quench", ev(6), stage="workup", pids=["p_NaHCO3", "p_3l"], operation="add saturated aqueous NaHCO3 and stir to quench unreacted acetic acid", time={"number": 1, "unit": "h"}),
    step("extract", ev(6), stage="workup", pids=["p_DCM", "p_3l"], operation="separate organic layer; extract aqueous layer with DCM (2 × 20 mL); combine and concentrate; dry under vacuum"),
]
r["procedure"]["original_text"] = f("To a 25 mL round-bottom flask containing 1a (1.495 g, 5 mmol) in DCE (0.4 M, 16.7 mL) was added acetic acid 2l (0.37 mL, 6.5 mmol) at room temperature; heat 55 °C 12 h; cool; NaHCO3 workup; DCM extraction; concentrate and vacuum-dry.", "Gram scale synthesis paragraph", ev(6))
r["procedure"]["addition_order"] = f("1a in DCE, then acetic acid 2l at room temperature", ["1a in DCE", "2l"], ev(6))
alt96 = f("Scheme 2: 5 mmol scale; 96%", {"number": 96, "unit": "%"}, ev(5))
alt95 = f("Experimental Section: 1.7 g ... (95% yield)", {"number": 95, "unit": "%"}, ev(6))
conflict = f(["5 mmol scale; 96%", "1.7 g; 95% yield"], None, e, status="conflict", verification="unreviewed", alternatives=[alt96, alt95])
oy = st.new_outcome("yield_conflict", "yield", [f"{DOC}:mol:3l"], e)
oy["value"] = conflict
oy["yield_type"] = f("isolated/dried product yield", "isolated/dried product yield", ev(6))
oy["measurement_basis"] = f("off-white solid obtained after concentration and vacuum drying; no column chromatography", "isolated dried product", ev(6))
r["outcomes"] = [oy, outcome("product_mass", "product_amount", [f"{DOC}:mol:3l"], evs(4, 6), "1.7 g", {"number": 1.7, "unit": "g"}, basis="dried under vacuum")]
set_hydro_class(r, e)
r["notes"] = [f("Source conflict: Scheme 2 reports 96% for the 5 mmol scale, while the detailed experimental paragraph reports 95%; both are retained without choosing.", "Source conflict retained", e),
              f("Table 2 footnote b says the yield was similar for a 5 mmol scale reaction.", "Table 2 footnote b", ev(3))]
r["review"]["unresolved"] = ["Conflicting 96% (Scheme 2) and 95% (experimental paragraph) yields require author/SI clarification.", "Atmosphere and pressure are unreported."]
reactions.append(r)


# Proposed DFT mechanism: not an experimental synthesis record.
e = ev(4)
r = base_reaction("figure1_2_proposed_mechanism", e, record_type="proposed_mechanism", categories=["mechanistic_study"], location="Figures 1-2 and DFT discussion", general=False)
r["reaction_name"] = f("Protonation First pathway", "Protonation First pathway", e)
r["transformation_type"] = f("stepwise protonation then 1,4-carboxylate addition", "stepwise protonation then 1,4-carboxylate addition", e)
r["notes"] = [f("PBE/6-31G* with SMD solvation corrections for DCE; protonation barrier 14.5 kcal mol-1, iminium-carboxylate complex ΔG = 3.6 kcal mol-1, product ΔG = -20.2 kcal mol-1.", "DFT energy profile", e),
              f("Concerted-1,2 was not a stationary point; Concerted-1,4 was geometrically impossible; carboxylation-first fragmented back to starting materials.", "Author DFT interpretation", e)]
no_class(r, e, "not_applicable")
r["review"]["unresolved"] = ["This is a proposed/computed mechanism, not a preparative experiment; no mechanistic experiment is reported on the supplied pages."]
reactions.append(r)


# Non-reaction source context: prior-work schemes and characterization/general remarks.
bundle["source_context"] = [
    {"context_id": f"{DOC}:context:scheme1_prior_work", "document_id": DOC, "kind": "cited_prior_work",
     "evidence_ids": ev(2), "note": "Scheme 1A summarizes prior transition-metal and Cu-catalyzed literature; it is not counted as a present-paper experiment."},
    {"context_id": f"{DOC}:context:general_characterization", "document_id": DOC, "kind": "characterization_methods",
     "evidence_ids": ev(5), "note": "General 1H/13C NMR, HRMS and melting-point instrumentation is context, not an independent synthesis experiment."},
    {"context_id": f"{DOC}:context:si_unavailable", "document_id": DOC, "kind": "availability_limit",
     "evidence_ids": evs(1, 5), "note": "The paper points to Supporting Information, but no SI was explicitly provided or associated in this session."},
]

bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions

raw_manifest = {
    "document_id": DOC,
    "pages_visually_inspected": [1, 2, 3, 4, 5, 6],
    "counts": {"molecules": len(bundle["molecules"]), "reactions": len(reactions)},
    "reaction_ids": [r["reaction_id"] for r in reactions],
    "table1_rows": 13,
    "table2_successes": 47,
    "table2_limitations": 2,
    "table3_present_paper_experiments": 14,
    "table3_prior_literature_rows": 4,
    "scheme2_derivatizations": 7,
    "scale_up_records": 1,
    "proposed_mechanism_records": 1,
    "semantic_accuracy_verified": False,
}
(RUN / "raw_model/model_extraction.json").write_text(json.dumps(raw_manifest, ensure_ascii=False, indent=2) + "\n")
(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
print(json.dumps(raw_manifest["counts"], ensure_ascii=False))
