#!/usr/bin/env python3
"""Build source-derived v3.1 records for original PDF pages 1--6.

This script is intentionally paper-specific.  It preserves visual transcription,
unknown structures, source conflicts, and missing analytical/procedural details.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))

import standard_records as st
from reaction_tools import dataset_audit

DOC = "doc_dfedfa74f4eea49f"
TITLE = (
    "Regio- and Stereoselective Synthesis of Multi-Alkylated Allylic "
    "Boronates through Three-Component Coupling Reactions between "
    "Allenes, Alkyl Halides, and a Diboron Reagent"
)

E = {p: f"{DOC}:page{p}:image" for p in range(1, 7)}
T = {p: f"{DOC}:page{p}:text" for p in range(1, 7)}


def f(raw=None, normalized=None, evidence=(), method="original_extraction",
      status=None, derivation=None):
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   derivation=derivation)


def q(number, unit, raw=None, qualifier=None):
    value = {"number": number, "unit": unit}
    if qualifier:
        value["qualifier"] = qualifier
    if raw:
        value["raw_text"] = raw
    return value


def ratio(raw):
    """Retain qualitative bounds and n.d. verbatim; ratios are not percentages."""
    return raw


dataset = json.loads((RUN / "dataset.json").read_text())
doc = next(d for d in dataset["documents"] if d["document_id"] == DOC)
doc["title"] = TITLE
dataset["documents"] = [doc]
st.write(RUN / "dataset.json", dataset)

bundle = st.empty(dataset["documents"])
molecules: dict[str, dict] = {}
reactions: list[dict] = []


def add_occurrence(m, label, evidence):
    ids = list(dict.fromkeys(evidence))
    if not any(o["evidence_ids"] == ids and
               o["compound_label"]["normalized_value"] == label
               for o in m["occurrences"]):
        m["occurrences"].append({
            "document_id": DOC,
            "compound_label": f(label, label, ids),
            "evidence_ids": ids,
        })


def mol(mid, label, evidence, *, smiles=None, english=None,
        chirality="unspecified", geometry="unresolved",
        representation=None, stereo_note=None, identity_note=None):
    evidence = list(dict.fromkeys(evidence))
    if mid in molecules:
        add_occurrence(molecules[mid], label, evidence)
        return mid
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    m["names"]["original"] = f(label, label, evidence)
    if english:
        m["names"]["english"] = f(english, english, evidence)
    if label and label[0].isdigit() or label.startswith("("):
        m["names"]["abbreviations"] = f(label, [label], evidence)
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = f(
            {"source_label": label, "candidate_smiles": smiles}, smiles, evidence,
            method="structure_derivation",
            derivation=("Candidate connectivity transcribed from the displayed structure or "
                        "unambiguous source chemical name; RDKit parsing is diagnostic only."),
        )
        diag = st.stereo_inspect(smiles)
        m["stereochemistry"]["isomeric_smiles"] = f(
            smiles, diag["canonical_isomeric_smiles"], evidence,
            method="structure_derivation",
            derivation="Canonical isomeric SMILES computed from the candidate graph; source identity is not thereby verified.",
        )
        m["stereochemistry"]["identity_key"] = f(
            None, diag["stereo_graph_key"], evidence, method="structure_derivation",
            derivation="RDKit stereo graph key from the candidate SMILES.",
        )
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = f(evidence=evidence)
    rep = representation or ("not_applicable" if smiles else "not_encoded")
    m["stereochemistry"]["chirality_status"] = f(
        {"source_structure": label}, chirality, evidence, method="structure_derivation",
        derivation="Status is limited to the stereochemical information visible in the cited source region.",
    )
    m["stereochemistry"]["geometry_status"] = f(
        {"source_structure": label}, geometry, evidence, method="structure_derivation",
        derivation="Status is limited to the stereochemical information visible in the cited source region.",
    )
    m["stereochemistry"]["representation_status"] = f(
        {"source_structure": label}, rep, evidence, method="structure_derivation",
        derivation="Complete source stereochemistry is not claimed beyond the candidate representation.",
    )
    if stereo_note:
        m["stereochemistry"]["assignment"] = f(stereo_note, stereo_note, evidence)
        if "E" in stereo_note or "Z" in stereo_note:
            m["stereochemistry"]["double_bonds"] = f(stereo_note, stereo_note, evidence)
    note = identity_note or (
        "Source-labelled molecule; complete connectivity/absolute stereochemistry was not safely encoded."
        if not smiles else
        "Source-labelled molecule with a machine candidate graph; structure remains unreviewed."
    )
    m["identity_note"] = f(note, note, evidence)
    molecules[mid] = m
    return mid


def product_pair(code, evidence, major="E"):
    """Register E and Z separately; reported yields may target both isomers."""
    e_id = mol(f"{DOC}:mol_E_{code}", f"(E)-{code}", evidence,
               chirality="unspecified", geometry="specified",
               representation="not_encoded", stereo_note="E alkene geometry explicitly reported/drawn")
    z_id = mol(f"{DOC}:mol_Z_{code}", f"(Z)-{code}", evidence,
               chirality="unspecified", geometry="specified",
               representation="not_encoded",
               stereo_note=("Z alkene geometry is the second component explicitly implied by the reported E/Z metric; "
                            f"{major} is the displayed/major order where stated"))
    return [e_id, z_id]


def regio_product(code, evidence):
    return mol(f"{DOC}:mol_{code}_regioisomer", f"{code} regioisomer (E and/or Z)", evidence,
               chirality="unspecified", geometry="unresolved", representation="not_encoded",
               stereo_note="Source reports the regioisomer fraction without a complete E/Z assignment")


# Concrete recurring inputs/reagents.  Complex catalysts and dense scope structures
# remain source-labelled without guessed SMILES.
M = {}
M["1a"] = mol(f"{DOC}:mol_1a", "1a", [E[4], E[5]], smiles="C=C=C(C)C1CCCCC1",
              english="cyclohexyl methyl allene", chirality="achiral", geometry="not_applicable",
              identity_note="gem-Dialkylallene 1a; terminal CH2 makes the displayed allene axis nonstereogenic.")
M["1b"] = mol(f"{DOC}:mol_1b", "1b", [E[4], E[5]], smiles="C=C=C(C)C(C)(C)C",
              english="tert-butyl methyl allene", chirality="achiral", geometry="not_applicable",
              identity_note="gem-Dialkylallene 1b; terminal CH2 makes the displayed allene axis nonstereogenic.")
for x in "cdefghijk":
    M[f"1{x}"] = mol(f"{DOC}:mol_1{x}", f"1{x}", [E[5]],
                       chirality="unspecified", geometry="not_applicable",
                       identity_note="Full substituent connectivity is visible in Table 2 but was not converted to a complete candidate SMILES; any central stereochemistry is not assigned by the source.")

halides = {
    "2a": ("1-iodo-4,4,4-trifluorobutane (2a)", "ICCCC(F)(F)F"),
    "2b": ("ethyl iodide (2b)", "CCI"),
    "2c": ("1-iodopentane (2c)", "CCCCCI"),
    "2d": ("phenylalkyl iodide 2d", None),
    "2e": ("1-iododecane (2e)", "CCCCCCCCCCI"),
    "2f": ("isobutyl iodide (2f)", "CC(C)CI"),
    "2g": ("benzyl bromide (2g)", "BrCc1ccccc1"),
    "2h": ("methoxy-functionalized alkyl iodide 2h", "COCCCI"),
    "2i": ("silyl-ether-functionalized alkyl iodide 2i", None),
    "2j": ("chloro-functionalized alkyl iodide 2j", "ClCCCI"),
    "2k": ("methyl iodide (2k)", "CI"),
    "2l": ("phenethyl iodide (2l)", "ICCc1ccccc1"),
    "2m": ("iodocyclohexane (2m)", "IC1CCCCC1"),
    "2n": ("2-bromo-2-methylbutane (2n)", "CCC(C)(Br)C"),
}
for key, (name, smi) in halides.items():
    M[key] = mol(f"{DOC}:mol_{key}", key, [E[5]] if key != "2a" else [E[4], E[5]],
                 smiles=smi, english=name, chirality="achiral" if smi else "unspecified",
                 geometry="not_applicable")

named = {
    "B2pin2": ("bis(pinacolato)diboron (3)", None),
    "KOtBu": ("potassium tert-butoxide", "[K+].CC(C)(C)[O-]"),
    "KOMe": ("potassium methoxide", "[K+].C[O-]"),
    "CuCl": ("copper(I) chloride", "[Cu]Cl"),
    "Xantphos": ("Xantphos", None),
    "dppp": ("dppp", None), "dppf": ("dppf", None), "PCy3": ("PCy3", None),
    "SIMesCuCl": ("SIMesCuCl", None), "IMesCuCl": ("IMesCuCl", None),
    "IAdCuCl": ("IAdCuCl", None),
    "DMF": ("N,N-dimethylformamide", "CN(C)C=O"),
    "THF": ("tetrahydrofuran", "C1CCOC1"),
    "MeCN": ("acetonitrile", "CC#N"),
    "DMA": ("N,N-dimethylacetamide", "CC(=O)N(C)C"),
    "toluene": ("toluene", "Cc1ccccc1"),
    "n-hexane": ("n-hexane", "CCCCCC"),
    "NaBO3.4H2O": ("sodium perborate tetrahydrate", None),
    "H2NOMe": ("methoxyamine", "CON"),
    "n-BuLi": ("n-butyllithium", "[Li]CCCC"),
    "Boc2O": ("di-tert-butyl dicarbonate", "CC(C)(C)OC(=O)OC(=O)OC(C)(C)C"),
    "CH2ClBr": ("bromochloromethane", "ClCBr"),
    "paraformaldehyde": ("paraformaldehyde", None),
    "BF3OEt2": ("boron trifluoride diethyl etherate", None),
    "DCM": ("dichloromethane", "ClCCl"),
    "benzaldehyde": ("benzaldehyde", "O=Cc1ccccc1"),
    "3bromo4fluorobenzaldehyde": ("3-bromo-4-fluorobenzaldehyde", None),
    "2naphthaldehyde": ("2-naphthaldehyde", "O=Cc1ccc2ccccc2c1"),
}
for key, (name, smi) in named.items():
    pages = [E[4], E[5]] if key in {"B2pin2", "KOtBu", "CuCl", "Xantphos", "SIMesCuCl", "DMF"} else [E[6]] if key in {"NaBO3.4H2O", "H2NOMe", "n-BuLi", "Boc2O", "CH2ClBr", "paraformaldehyde", "BF3OEt2", "DCM", "benzaldehyde", "3bromo4fluorobenzaldehyde", "2naphthaldehyde"} else [E[4]]
    M[key] = mol(f"{DOC}:mol_{key.replace('.', '_')}", name, pages, smiles=smi,
                 chirality="achiral" if smi else "unspecified", geometry="not_applicable")


def new_rx(suffix, evidence, *, record_type="specific_experiment", categories=None,
           location=None, name=None, transformation=None, unresolved=None):
    r = st.new_reaction(f"{DOC}:{suffix}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, evidence) if categories is not None else f(evidence=evidence)
    r["location"] = f(location, location, evidence)
    if name:
        r["reaction_name"] = f(name, name, evidence)
    if transformation:
        r["transformation_type"] = f(
            None, transformation, evidence, method="inference",
            derivation="Transformation class derived from the displayed reactant/product connectivity.",
        )
    r["review"] = {
        "status": "needs_review", "method": "model", "reviewer": None,
        "checked_fields": [],
        "unresolved": list(unresolved or ["Complete structures and unavailable SI procedures remain unresolved."]),
    }
    r["completeness"] = "partial"
    reactions.append(r)
    return r


def addp(r, pid, mid, role, evidence, *, stage="reaction", amount=None,
         equivalents=None, mol_percent=None, volume=None, reference=None,
         product_role=None):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    if product_role:
        p["product_role"] = f(product_role, product_role, evidence)
    for key, value in (("amount", amount), ("equivalents", equivalents),
                       ("mol_percent", mol_percent), ("volume", volume)):
        if value is not None:
            raw, norm = value
            p["charges"][key] = f(raw, norm, evidence)
    if any(v is not None for v in (equivalents, mol_percent)):
        p["equivalence_reference"] = f(reference, reference, evidence)
    r["participants"].append(p)
    return p


def add_step(r, sid, evidence, *, operation=None, temperature=None, time=None,
             participant_ids=None, stage="reaction"):
    s = st.new_step(sid, evidence, stage=stage)
    s["participant_ids"] = list(participant_ids or [])
    if operation is not None:
        s["operation"] = f(operation, operation, evidence)
    if temperature is not None:
        s["temperature"] = f(temperature[0], temperature[1], evidence)
    if time is not None:
        s["time"] = f(time[0], time[1], evidence)
    r["steps"].append(s)
    return s


def outcome(r, oid, kind, targets, evidence, raw, normalized, *, role=None,
            yield_type=None, basis=None, scope=None):
    o = st.new_outcome(oid, kind, targets, evidence)
    o["value"] = f(raw, normalized, evidence)
    if role:
        o["product_role"] = f(role, role, evidence)
    if yield_type:
        o["yield_type"] = f(yield_type, yield_type, evidence)
    if basis and basis != "not reported":
        o["measurement_basis"] = f(basis, basis, evidence)
    if scope:
        o["scope"] = f(scope, scope, evidence)
    r["outcomes"].append(o)
    return o


def note(r, text, evidence):
    r["notes"].append(f(text, text, evidence))


def classify(r, evidence, specs=None, coverage="partial", limitation=None):
    c = r["classification"]
    if coverage in {"unresolved", "not_applicable"}:
        c["coverage"] = f(None, coverage, evidence, method="inference",
                          derivation="No observed/source-resolved bond transformation is asserted for this record.")
        if coverage == "not_applicable":
            c["labels"] = f(None, [], evidence, method="inference",
                            derivation="Failed/no-reaction outcome; no expected-product label assigned.")
        c["limitations"] = [limitation or "No source-resolved atom mapping."]
        return
    changes = []
    grouped = {}
    for i, spec in enumerate(specs or [], 1):
        a, b, before, after, kind = spec
        ch = {
            "change_id": f"b{i}",
            "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
            "before_order": before, "after_order": after, "change_type": kind,
        }
        changes.append(ch)
        if kind in {"formed", "order_increased", "hydrogen_gain"}:
            lab = st.bond_label(a, b, after)
            grouped.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(ch["change_id"])
    derivation = "Partial bond changes derived from the displayed transformation; atom identities are not mapped."
    c["bond_changes"] = [f(None, ch, evidence, method="structure_derivation", derivation=derivation) for ch in changes]
    c["labels"] = f(None, list(grouped.values()), evidence, method="structure_derivation", derivation=derivation)
    c["coverage"] = f(None, coverage, evidence, method="inference",
                      derivation="Only source-visible connectivity classes are represented; full atom mapping is unavailable.")
    c["limitations"] = [limitation or "No complete atom mapping; leaving-group, coproduct, and implicit-H balance remain unresolved."]


CARBO = [("C", "C", 0, 1, "formed"), ("C", "B", 0, 1, "formed"),
         ("C", "C", 2, 1, "order_decreased")]
ALLYLBORATION = [("C", "C", 0, 1, "formed"), ("C", "C", 1, 2, "order_increased"),
                 ("C", "C", 2, 1, "order_decreased"), ("C", "O", 2, 1, "order_decreased"),
                 ("O", "H", 0, 1, "formed"), ("C", "B", 1, 0, "broken")]


def common_carboboration(r, allene_key, halide_key, evidence, condition, *,
                         catalyst_key=None, ligand_key=None, base_key="KOtBu",
                         solvent_key="DMF", temperature=None, time_h=24,
                         catalyst_mol=None, catalyst_amount=None, ligand_mol=None,
                         ligand_amount=None, concentration=None, solvent_volume=None):
    ref = M[allene_key]
    addp(r, "allene", ref, "reactant", evidence,
         amount=("0.5 mmol", q(0.5, "mmol")), equivalents=("1.0 equiv", q(1.0, "equiv")), reference=ref)
    addp(r, "halide", M[halide_key], "reactant", evidence,
         amount=("1.0 mmol", q(1.0, "mmol")), equivalents=("2.0 equiv", q(2.0, "equiv")), reference=ref)
    addp(r, "diboron", M["B2pin2"], "reagent", evidence,
         amount=("0.6 mmol", q(0.6, "mmol")), equivalents=("1.2 equiv", q(1.2, "equiv")), reference=ref)
    if catalyst_key:
        addp(r, "catalyst", M[catalyst_key], "catalyst", evidence,
             amount=(f"{catalyst_amount} mmol", q(catalyst_amount, "mmol")) if catalyst_amount is not None else None,
             mol_percent=(f"{catalyst_mol} mol%", q(catalyst_mol, "mol%")) if catalyst_mol is not None else None,
             reference=ref)
    if ligand_key:
        addp(r, "ligand", M[ligand_key], "ligand", evidence,
             amount=(f"{ligand_amount} mmol", q(ligand_amount, "mmol")) if ligand_amount is not None else None,
             mol_percent=(f"{ligand_mol} mol%", q(ligand_mol, "mol%")) if ligand_mol is not None else None,
             reference=ref)
    addp(r, "base", M[base_key], "reagent", evidence,
         amount=("0.6 mmol", q(0.6, "mmol")), equivalents=("1.2 equiv", q(1.2, "equiv")), reference=ref)
    addp(r, "solvent", M[solvent_key], "solvent", evidence,
         volume=(f"{solvent_volume} mL", q(solvent_volume, "mL")) if solvent_volume is not None else None)
    r["reaction_scale"] = f(
        f"{allene_key} (0.5 mmol)", {"reference_molecule_id": ref, "amount": q(0.5, "mmol")}, evidence)
    r["solvent_conditions"]["composition"] = f(solvent_key, [solvent_key], evidence)
    if solvent_volume is not None:
        r["solvent_conditions"]["total_volume"] = f(f"{solvent_volume} mL", q(solvent_volume, "mL"), evidence)
    if concentration is not None:
        r["solvent_conditions"]["concentration"] = f(f"{concentration} M", q(concentration, "M"), evidence)
        r["solvent_conditions"]["concentration_basis"] = f(allene_key, ref, evidence)
    temp = temperature if temperature is not None else (30 if condition == "B" else -5)
    add_step(r, "reaction", evidence,
             operation=f"carboboration under conditions {condition}",
             temperature=(f"{temp} °C", q(temp, "°C")),
             time=(f"{time_h} h", q(time_h, "h")),
             participant_ids=[p["participant_id"] for p in r["participants"]])


# One generic current-work reaction, shown in the graphical abstract and Figure 2C.
generic = new_rx(
    "generic_allylic_boronate_carboboration", [E[1], E[3]], record_type="generic_scheme",
    categories=["other"], location="page 1 graphical abstract and Figure 2C",
    name="three-component coupling reaction / allylic boronate-selective carboboration",
    transformation="copper-catalyzed carboboration of gem-dialkylallenes",
    unresolved=["Generic R groups and exact conditions are intentionally unresolved."],
)
g_allene = mol(f"{DOC}:mol_generic_gem_dialkylallene", "generic gem-dialkylallene", [E[1], E[3]], chirality="unspecified", geometry="not_applicable")
g_halide = mol(f"{DOC}:mol_generic_alkyl_halide", "generic alkyl halide", [E[1], E[3]], chirality="unspecified", geometry="not_applicable")
g_product = mol(f"{DOC}:mol_generic_allylic_boronate", "generic 2,3,3-trialkyl-substituted allylic boronate", [E[1], E[3]], chirality="unspecified", geometry="specified", representation="not_encoded", stereo_note="stereodefined tetrasubstituted alkene; E/Z up to >95:5")
for pid, mid, role in (("allene", g_allene, "reactant"), ("halide", g_halide, "reactant"), ("diboron", M["B2pin2"], "reagent"), ("product", g_product, "product")):
    addp(generic, pid, mid, role, [E[1], E[3]], product_role="main product" if role == "product" else None)
add_step(generic, "generic_reaction", [E[1], E[3]], operation="Cu-catalyzed three-component coupling with base", participant_ids=[p["participant_id"] for p in generic["participants"]])
classify(generic, [E[1], E[3]], CARBO)


# Figure 1B: explicitly cited Hall literature precedent, not current-work data.
hall = new_rx(
    "figure1B_hall_literature_example", [E[2]], record_type="literature_example",
    location="Figure 1B", name="Hall's synthesis of differentially 2,3,3-trisubstituted allylic boronates",
    transformation="conjugate organocuprate addition followed by electrophilic borylation",
    unresolved=["Generic R groups, EWG identity, quantities, temperature, time, and yield are not specified in Figure 1B."],
)
for pid, label, role in (("alkyne", "generic alkynyl EWG substrate", "reactant"),
                         ("cuprate", "generic dialkylcuprate", "reagent"),
                         ("boryl", "electrophilic B(pin) reagent", "reagent"),
                         ("product", "generic Hall allylic boronate", "product")):
    mid = mol(f"{DOC}:mol_hall_{pid}", label, [E[2]], chirality="unspecified", geometry="specified" if pid == "product" else "unresolved", representation="not_encoded",
              stereo_note="Figure 1B reports the literature product as typically E/Z >20:1" if pid == "product" else None)
    addp(hall, pid, mid, role, [E[2]], product_role="literature-example product" if role == "product" else None)
addp(hall, "solvent", M["THF"], "solvent", [E[2]])
add_step(hall, "literature_scheme", [E[2]], operation="two-stage in situ literature scheme", participant_ids=[p["participant_id"] for p in hall["participants"]])
classify(hall, [E[2]], CARBO, limitation="Generic cited literature transformation; no complete atom mapping or current-work experiment is claimed.")


# Table 1 optimization rows.
table1 = [
    (1, "1a", "CuCl", "Xantphos", 5, 30, None, 88, 88, "60:40", ">99:1"),
    (2, "1a", "CuCl", "dppp", 5, 30, None, 30, 68, "70:30", "44:56"),
    (3, "1a", "CuCl", "dppf", 5, 30, None, 5, 80, "n.d.", "6:94"),
    (4, "1a", "CuCl", "PCy3", 5, 30, None, 2, 3, "n.d.", "81:19"),
    (5, "1a", "SIMesCuCl", None, 2, 0, None, 73, 79, ">95:5", "92:8"),
    (6, "1a", "IMesCuCl", None, 2, 0, None, 55, 58, ">95:5", "94:6"),
    (7, "1a", "IAdCuCl", None, 2, 0, None, 71, 72, "32:68", "98:2"),
    (8, "1a", "SIMesCuCl", None, 2, 0, "THF", 19, 22, "n.d.", "87:13"),
    (9, "1a", "SIMesCuCl", None, 2, 0, "MeCN", 13, 16, "n.d.", "84:16"),
    (10, "1a", "SIMesCuCl", None, 2, 0, "DMA", 48, 56, ">95:5", "85:15"),
    (11, "1a", "SIMesCuCl", None, 2, 0, "toluene", "<1", "<1", "n.d.", "n.d."),
    (12, "1a", "SIMesCuCl", None, 2, 0, "n-hexane", 5, 8, "n.d.", "60:40"),
    (13, "1a", "SIMesCuCl", None, 2, 0, "base:KOMe", 5, 6, "n.d.", "79:21"),
    (14, "1a", "SIMesCuCl", None, 2, -5, None, 77, 83, ">95:5", "93:7"),
    (15, "1a", "SIMesCuCl", None, 2, -10, None, 51, 54, ">95:5", "94:6"),
    (16, "1a", "SIMesCuCl", None, 2, 0, "concentration:0.25M", 68, 72, ">95:5", "95:5"),
]


def percent_norm(v):
    if isinstance(v, (int, float)):
        return q(v, "%")
    if v == "<1":
        return q(1, "%", "<1%", "lt")
    return v


for entry, allene, catalyst, ligand, loading, temp, variation, y4, ytotal, ez, rr in table1:
    evidence = [E[4], T[4]] + ([T[2]] if entry <= 9 else [T[3]])
    r = new_rx(f"table1_entry_{entry:02d}", evidence, categories=["optimization"],
               location=f"Table 1 entry {entry}", name="Reaction optimization",
               transformation="copper-catalyzed alkylboration of a gem-dialkylallene",
               unresolved=["Complete product structures and SI procedure are unavailable."])
    solvent = variation if variation in {"THF", "MeCN", "DMA", "toluene", "n-hexane"} else "DMF"
    base = "KOMe" if variation == "base:KOMe" else "KOtBu"
    conc = 0.25 if variation == "concentration:0.25M" else 0.5
    vol = None if variation == "concentration:0.25M" else 1.0
    common_carboboration(r, allene, "2a", [E[4]], "optimization",
                         catalyst_key=catalyst, ligand_key=ligand, base_key=base,
                         solvent_key=solvent, temperature=temp, catalyst_mol=loading,
                         catalyst_amount=0.025, ligand_mol=loading if ligand else None,
                         concentration=conc, solvent_volume=vol)
    code = "4aa"
    products = product_pair(code, [E[4]])
    reg = regio_product("5aa", [E[4]])
    for j, mid in enumerate(products):
        addp(r, f"product_{j+1}", mid, "product", [E[4]], product_role="major E / minor Z product fraction")
    addp(r, "regioisomer", reg, "product", [E[4]], product_role="regioisomer fraction")
    outcome(r, "yield_4", "yield", products, [E[4]], f"{y4}%" if y4 != "<1" else "<1%", percent_norm(y4),
            role="4aa product fraction", yield_type="GC yield", basis="GC analysis with an internal standard",
            scope="yield of 4xa")
    outcome(r, "combined_yield", "yield", products + [reg], [E[4]], f"({ytotal}%)" if ytotal != "<1" else "(<1%)", percent_norm(ytotal),
            role="4aa + 5aa", yield_type="combined GC yield", basis="GC analysis with an internal standard",
            scope="parenthesized yield of 4xa+5xa")
    if ez != "n.d.":
        outcome(r, "E_Z", "selectivity", products, [E[4]], ez, ratio(ez),
                basis="13C NMR analysis of roughly purified material", scope="E:Z ratio of 4aa; order E:Z")
    else:
        outcome(r, "E_Z_not_determined", "other", products, [E[4]], "n.d.", "not determined",
                scope="E/Z of 4aa")
    if rr != "n.d.":
        outcome(r, "4_5", "selectivity", products + [reg], [E[4]], rr, ratio(rr),
                basis="GC analysis of the reaction mixture using an internal standard", scope="4aa:5aa regioisomer ratio")
    else:
        outcome(r, "4_5_not_determined", "other", products + [reg], [E[4]], "n.d.", "not determined", scope="4aa:5aa")
    if loading == 2:
        msg = ("Table 1 footnote a gives Cu(I) catalyst 0.025 mmol, while the row gives 2 mol% "
               "for 0.5 mmol allene; both source values are retained without choosing one.")
        note(r, msg, [E[4]])
        r["review"]["unresolved"].append("Catalyst charge conflict: 0.025 mmol versus row loading 2 mol% (0.5 mmol basis).")
    if entry == 16:
        note(r, "The 0.25 M variation supersedes the 0.5 M standard; exact solvent volume is not printed and was not inferred.", [E[4]])
    r["general_scheme_ids"] = [generic["reaction_id"]]
    classify(r, [E[4]], CARBO)


# Table 1 entries 17--18 (1b); entry 18 also appears as 4ba in Table 2.
for entry, catalyst, ligand, loading, temp, y, ez, rr in [
    (17, "SIMesCuCl", None, 2, -5, 58, "77:23", None),
    (18, "CuCl", "Xantphos", 5, 30, 83, ">95:5", ">99:1"),
]:
    evidence = [E[4], T[4]] + ([E[5], T[5]] if entry == 18 else [])
    r = new_rx(f"table1_entry_{entry:02d}", evidence, categories=["optimization"],
               location=f"Table 1 entry {entry}" + ("; repeated as (E)-4ba in Table 2" if entry == 18 else ""),
               name="Reaction optimization", transformation="copper-catalyzed alkylboration of a gem-dialkylallene",
               unresolved=["Complete product structures and SI procedure are unavailable."])
    common_carboboration(r, "1b", "2a", [E[4]] + ([E[5]] if entry == 18 else []), "B" if entry == 18 else "A",
                         catalyst_key=catalyst, ligand_key=ligand, temperature=temp,
                         catalyst_mol=loading, catalyst_amount=0.025,
                         ligand_mol=loading if ligand else None, ligand_amount=0.025 if ligand else None,
                         concentration=0.5, solvent_volume=1.0)
    products = product_pair("4ba", evidence)
    for j, mid in enumerate(products):
        addp(r, f"product_{j+1}", mid, "product", evidence, product_role="4ba product fraction")
    if entry == 17:
        minor = mol(f"{DOC}:mol_entry17_unidentified_minor", "unidentified minor isomer in Table 1 entry 17", [E[4]], chirality="unresolved", geometry="unresolved", representation="not_encoded")
        addp(r, "unidentified_minor", minor, "product", [E[4]], product_role="unidentified minor isomer")
        outcome(r, "yield", "yield", products, [E[4]], "58%", q(58, "%"), role="4ba-containing product", yield_type="isolated yield", basis="isolated yield")
        outcome(r, "isomer_ratio", "selectivity", [products[0], minor], [E[4]], "77:23 mixture of isomers", "77:23",
                basis="1H NMR analysis of purified material", scope="major product:unidentified minor isomer; not assigned as a definitive E/Z or 4/5 ratio")
        note(r, "Product contained small amounts of protoboration products of 1b and a boryl-substitution product of 2a.", [E[4]])
        note(r, "The structure of the minor isomer could not be identified.", [E[4]])
        r["review"]["unresolved"].extend(["Minor isomer structure is unidentified.", "Catalyst charge conflict: 0.025 mmol versus row loading 2 mol% (0.5 mmol basis)."])
    else:
        reg = regio_product("5ba", evidence)
        addp(r, "regioisomer", reg, "product", evidence, product_role="regioisomer fraction")
        outcome(r, "yield", "yield", products, evidence, "83%", q(83, "%"), role="4ba product fraction", yield_type="isolated yield", basis="isolated yield")
        outcome(r, "E_Z", "selectivity", products, evidence, ">95:5", ">95:5",
                basis=["13C NMR of roughly purified material (Table 1)", "1H NMR after column chromatography (Table 2)"],
                scope="E:Z ratio of 4ba")
        outcome(r, "4_5", "selectivity", products + [reg], evidence, ">99:1", ">99:1",
                basis=["GC reaction-mixture analysis (Table 1)", "1H NMR after column chromatography (Table 2)"],
                scope="4ba:5ba; Table 2 rounds/generalizes this as >95:5")
        note(r, "The Table 2 4ba panel repeats the Table 1 entry 18 conditions, 83% yield, and high selectivity; it is treated as duplicate evidence, not a second experiment.", evidence)
    r["general_scheme_ids"] = [generic["reaction_id"]]
    classify(r, evidence, CARBO)


# Table 2 successful entries.  4ba is deliberately omitted here because it is the
# duplicate presentation already merged into Table 1 entry 18.
scope_rows = [
    ("4aa", "1a", "2a", "A", 66, ">95:5", ">95:5", "bd"),
    ("4ab", "1a", "2b", "A", 62, ">95:5", ">95:5", "bd"),
    ("4ac", "1a", "2c", "A", 49, ">95:5", ">95:5", "bde"),
    ("4ad", "1a", "2d", "A", 63, ">95:5", ">95:5", "bd"),
    ("4cb", "1c", "2b", "B", 88, "88:12", ">95:5", "cf"),
    ("4cc", "1c", "2c", "B", 89, "88:12", ">95:5", "cf"),
    ("4db", "1d", "2b", "B", 87, "89:11", ">95:5", "cf"),
    ("4dc", "1d", "2c", "B", 87, "89:11", ">95:5", "cf"),
    ("4eb", "1e", "2b", "B", 82, "89:11", ">95:5", "cf"),
    ("4ec", "1e", "2c", "B", 50, "87:13", ">95:5", "cd"),
    ("4fb", "1f", "2b", "B", 79, "89:11", ">95:5", "cf"),
    ("4fc", "1f", "2c", "B", 76, "91:9", ">95:5", "cf"),
    ("4ga", "1g", "2a", "B", 64, ">95:5", ">95:5", "cg"),
    ("4ha", "1h", "2a", "B", 69, ">95:5", ">95:5", "ch"),
    ("4hc", "1h", "2c", "B", 92, ">95:5", ">95:5", "c"),
    ("4ia", "1i", "2a", "A", 59, "89:11", ">95:5", "bdi"),
    ("4ic", "1i", "2c", "A", 52, "90:10", ">95:5", "bd"),
    ("4ja", "1j", "2a", "A", 64, "92:8", ">95:5", "bd"),
    ("4jc", "1j", "2c", "A", 46, "94:6", "95:5", "bd"),
    ("4kc", "1k", "2c", "B", 67, "50:50", ">95:5", "cd"),
    ("4bb", "1b", "2b", "B", 75, ">95:5", ">95:5", "c"),
    ("4bc", "1b", "2c", "B", 83, ">95:5", ">95:5", "c"),
    ("4be", "1b", "2e", "B", 88, ">95:5", ">95:5", "c"),
    ("4bf", "1b", "2f", "B", 81, "91:9", ">95:5", "c"),
    ("4bg", "1b", "2g", "B", 55, ">95:5", ">95:5", "cj"),
    ("4bh", "1b", "2h", "B", 74, ">95:5", ">95:5", "c"),
    ("4bi", "1b", "2i", "B", 44, ">95:5", ">95:5", "cd"),
    ("4bj", "1b", "2j", "B", 51, ">95:5", ">95:5", "cdk"),
]

oxidation_derived = set("4aa 4ab 4ac 4ad 4ec 4ia 4ic 4ja 4jc 4kc 4bi 4bj".split())
for code, allene, halide, cond, y, ez, rr, flags in scope_rows:
    evidence = [E[5], T[5]]
    r = new_rx(f"table2_{code}", evidence, categories=["substrate_scope"],
               location=f"Table 2, {code}", name="Screening the substrate scope",
               transformation="copper-catalyzed alkylboration of a gem-dialkylallene",
               unresolved=["Full dense substrate/product structures and SI procedure are unavailable."])
    temp = 30 if "e" in flags else 0 if "i" in flags else (30 if cond == "B" else -5)
    time_h = 48 if "h" in flags else 24
    if cond == "A":
        common_carboboration(r, allene, halide, evidence, "A", catalyst_key="SIMesCuCl",
                             catalyst_mol=2, catalyst_amount=0.01, temperature=temp, time_h=time_h,
                             concentration=0.5, solvent_volume=1.0)
    else:
        common_carboboration(r, allene, halide, evidence, "B", catalyst_key="CuCl", ligand_key="Xantphos",
                             catalyst_mol=5, catalyst_amount=0.025, ligand_mol=5, ligand_amount=0.025,
                             temperature=temp, time_h=time_h, concentration=0.5, solvent_volume=1.0)
    products = product_pair(code, evidence)
    for jx, mid in enumerate(products):
        addp(r, f"product_{jx+1}", mid, "product", evidence, product_role="E/Z borylation-product fraction")
    reg = regio_product("5" + code[1:], evidence)
    addp(r, "regioisomer", reg, "product", evidence, product_role="regioisomer fraction")
    targets = products
    class_specs = list(CARBO)
    if code in oxidation_derived:
        derivative = mol(f"{DOC}:mol_oxidized_{code}", f"oxidized derivative used to determine the isolated yield of {code}", evidence,
                         chirality="unspecified", geometry="specified", representation="not_encoded",
                         stereo_note=f"Derived from the E/Z {code} product fraction; exact derivative stereostructure not encoded")
        addp(r, "oxidized_isolation_product", derivative, "product", evidence, stage="workup", product_role="actual isolated oxidation derivative")
        add_step(r, "oxidation_before_isolation", evidence, stage="workup",
                 operation="oxidation of the boryl group before isolation (reagents/temperature/time not reported here)",
                 participant_ids=["product_1", "product_2", "oxidized_isolation_product"])
        targets = [derivative]
        class_specs += [("C", "B", 1, 0, "broken"), ("C", "O", 0, 1, "formed"), ("O", "H", 0, 1, "formed")]
        basis = "isolated yield after oxidation of the boryl group"
        r["review"]["unresolved"].append("Oxidation reagents and conditions are not printed in Table 2 and no SI was supplied.")
    else:
        basis = "isolated yield of the borylation product"
    outcome(r, "yield", "yield", targets, evidence, f"{y}%", q(y, "%"),
            role="isolated product", yield_type=basis, basis=basis,
            scope="Table 2 isolated yield")
    selectivity_basis = "GC analysis after silica gel column chromatography" if "f" in flags else "1H NMR analysis after column chromatography"
    outcome(r, "E_Z", "selectivity", products, evidence, ez, ratio(ez), basis=selectivity_basis,
            scope=f"E:Z ratio of {code}; order E:Z")
    outcome(r, "4_5", "selectivity", products + [reg], evidence, rr, ratio(rr), basis=selectivity_basis,
            scope=f"{code}:5{code[1:]} regioisomer ratio")
    if "g" in flags:
        note(r, "A small amount of byproduct, possibly the borylation product of the alkyl iodide, was removed by selective oxidation.", evidence)
    if code == "4kc":
        note(r, "The source explicitly calls this an unsuccessful allene case because E/Z = 50:50 (no stereoselectivity), despite 67% yield.", evidence)
    if code in {"4bi", "4bj"}:
        note(r, "The reaction produced the alkylboration product together with unidentified side products.", evidence)
        r["review"]["unresolved"].append("Unidentified side products are reported but not structurally resolved.")
    if "k" in flags:
        alt = mol(f"{DOC}:mol_{code}_chloride_site_product", f"minor product from C-C bond formation at chloride in {code} experiment", evidence,
                  chirality="unspecified", geometry="unresolved", representation="not_encoded")
        addp(r, "chloride_site_product", alt, "product", evidence, product_role="minor chemoselectivity product")
        outcome(r, "iodide_chloride_selectivity", "selectivity", products + [alt], evidence, "86:14", "86:14",
                basis="not reported", scope="C-C bond formation at iodide:at chloride")
    r["general_scheme_ids"] = [generic["reaction_id"]]
    classify(r, evidence, class_specs)


# Four unsuccessful alkyl-halide scope entries under Conditions B.
for halide, raw_result, kind in [
    ("2k", "decomposition of halide", "decomposition"),
    ("2l", "decomposition of halide", "decomposition"),
    ("2m", "complex mixture", "other"),
    ("2n", "no reaction", "no_reaction"),
]:
    evidence = [E[5], T[5]]
    r = new_rx(f"table2_unsuitable_{halide}", evidence, categories=["substrate_scope", "control"],
               location=f"Table 2, unsuitable halide {halide}", name="Unsuccessful alkyl-halide scope experiment",
               transformation="attempted copper-catalyzed alkylboration",
               unresolved=["No product identity or analytical detection limit is reported."])
    common_carboboration(r, "1b", halide, evidence, "B", catalyst_key="CuCl", ligand_key="Xantphos",
                         catalyst_mol=5, catalyst_amount=0.025, ligand_mol=5, ligand_amount=0.025,
                         temperature=30, time_h=24, concentration=0.5, solvent_volume=1.0)
    outcome(r, "reported_result", kind, [M[halide]], evidence, raw_result, raw_result,
            basis="not reported", scope=f"reported outcome for {halide}")
    if halide == "2l":
        styrene = mol(f"{DOC}:mol_styrene", "styrene side product", evidence, smiles="C=Cc1ccccc1", chirality="achiral", geometry="not_applicable")
        addp(r, "styrene", styrene, "product", evidence, product_role="reported side product")
        outcome(r, "styrene_side_product", "other", [styrene], evidence, "styrene as a side product", "styrene observed",
                basis="not reported", scope="side product attributed by authors to beta-elimination")
    r["general_scheme_ids"] = [generic["reaction_id"]]
    classify(r, evidence, coverage="not_applicable", limitation="Failed/complex/decomposition outcome; no expected-product bond labels are assigned.")


# Page 6 Figure 3A gram-scale carboboration.
ev6 = [E[6], T[6]]
scale = new_rx("figure3A_gram_scale_4hc", ev6, categories=["scale_up"], location="Figure 3A",
               name="Gram-scale synthesis of allylic boronate", transformation="copper-catalyzed alkylboration of gem-dialkylallene 1h",
               unresolved=["Atmosphere, addition order, workup, purification, and yield measurement basis are not printed on page 6."])
ref = M["1h"]
addp(scale, "allene", ref, "reactant", ev6, amount=("6.0 mmol", q(6.0, "mmol")), equivalents=("1.0 equiv", q(1.0, "equiv")), reference=ref)
addp(scale, "halide", M["2c"], "reactant", ev6, equivalents=("2.0 equiv", q(2.0, "equiv")), reference=ref)
addp(scale, "diboron", M["B2pin2"], "reagent", ev6, equivalents=("1.2 equiv", q(1.2, "equiv")), reference=ref)
addp(scale, "catalyst", M["CuCl"], "catalyst", ev6, mol_percent=("5 mol%", q(5, "mol%")), reference=ref)
addp(scale, "ligand", M["Xantphos"], "ligand", ev6, mol_percent=("5 mol%", q(5, "mol%")), reference=ref)
addp(scale, "base", M["KOtBu"], "reagent", ev6, equivalents=("1.2 equiv", q(1.2, "equiv")), reference=ref)
addp(scale, "solvent", M["DMF"], "solvent", ev6)
prod4hc = product_pair("4hc", [E[5], E[6]])
for jx, mid in enumerate(prod4hc): addp(scale, f"product_{jx+1}", mid, "product", ev6, product_role="E/Z 4hc product fraction")
scale["reaction_scale"] = f("1h (6.0 mmol)", {"reference_molecule_id": ref, "amount": q(6.0, "mmol")}, ev6)
scale["solvent_conditions"]["composition"] = f("DMF", ["DMF"], ev6)
scale["solvent_conditions"]["concentration"] = f("0.5 M", q(0.5, "M"), ev6)
scale["solvent_conditions"]["concentration_basis"] = f("1h", ref, ev6)
add_step(scale, "reaction", ev6, operation="gram-scale carboboration", temperature=("30 °C", q(30, "°C")), time=("24 h", q(24, "h")), participant_ids=[p["participant_id"] for p in scale["participants"]])
outcome(scale, "yield", "yield", prod4hc, ev6, "88%", q(88, "%"), role="4hc product fraction", yield_type="reported yield", scope="gram-scale product yield")
outcome(scale, "product_amount", "product_amount", prod4hc, ev6, "1.86 g", q(1.86, "g"), role="4hc product fraction", scope="isolated/reported product mass; isolation basis not explicitly stated")
outcome(scale, "E_Z", "selectivity", prod4hc, ev6, ">95:5", ">95:5", basis="not reported", scope="E:Z of 4hc")
reg5hc = regio_product("5hc", ev6); addp(scale, "regioisomer", reg5hc, "product", ev6, product_role="regioisomer fraction")
outcome(scale, "4_5", "selectivity", prod4hc + [reg5hc], ev6, ">95:5", ">95:5", basis="not reported", scope="4hc:5hc")
scale["general_scheme_ids"] = [generic["reaction_id"]]
classify(scale, ev6, CARBO)


# Figure 3B derivatizations from (E)-4hc.
def add_4hc_start(r):
    addp(r, "starting_4hc", prod4hc[0], "reactant", ev6)


for suffix, label, major_geom, y, reagents, steps, specs in [
    ("figure3B_oxidation_6hc", "6hc", "Z", 89, ["NaBO3.4H2O", "THF"],
     [("oxidation", "NaBO3·4H2O in THF/H2O", ("rt", "rt"), ("21 h", q(21, "h")))],
     [("C", "B", 1, 0, "broken"), ("C", "O", 0, 1, "formed"), ("O", "H", 0, 1, "formed")]),
    ("figure3B_amination_7hc", "7hc", "Z", 46, ["H2NOMe", "n-BuLi", "Boc2O"],
     [("amination", "H2NOMe and n-BuLi; −78 to 60 °C", ("−78 to 60 °C", {"minimum": -78, "maximum": 60, "unit": "°C", "raw_text": "−78 to 60 °C"}), ("48 h", q(48, "h"))),
      ("Boc_protection", "then (Boc)2O", ("rt", "rt"), ("2 h", q(2, "h")))],
     [("C", "B", 1, 0, "broken"), ("C", "N", 0, 1, "formed")]),
    ("figure3B_homologation_8hc", "8hc", "E", 73, ["CH2ClBr", "n-BuLi"],
     [("low_temperature_addition", "CH2ClBr and n-BuLi at −78 °C; stage duration not separately reported", ("−78 °C", q(-78, "°C")), None),
      ("warm_to_room_temperature", "warm from −78 °C to rt; the source reports 21 h for the complete program", ("rt", "rt"), ("21 h overall", q(21, "h", "21 h overall")))],
     [("C", "B", 1, 0, "broken"), ("C", "C", 0, 1, "formed"), ("C", "B", 0, 1, "formed")]),
]:
    r = new_rx(suffix, ev6, categories=["derivatization"], location="Figure 3B",
               name="alpha-selective transformation of the boryl group", transformation=label,
               unresolved=["Reagent quantities, atmosphere, workup, purification, and measurement methods are not printed on page 6."])
    add_4hc_start(r)
    for ix, reagent in enumerate(reagents): addp(r, f"reagent_{ix+1}", M[reagent], "solvent" if reagent == "THF" else "reagent", ev6)
    product_ids = product_pair(label, ev6, major=major_geom)
    major_index = 0 if major_geom == "E" else 1
    for ix, mid in enumerate(product_ids): addp(r, f"product_{ix+1}", mid, "product", ev6, product_role=f"{major_geom}-major product fraction")
    for sid, op, temp, tm in steps:
        add_step(r, sid, ev6, operation=op, temperature=temp, time=tm, participant_ids=[p["participant_id"] for p in r["participants"]])
    outcome(r, "yield", "yield", product_ids, ev6, f"{y}%", q(y, "%"), role=f"{label} product fraction", yield_type="reported yield", scope=f"yield of {major_geom}-{label}")
    outcome(r, "E_Z", "selectivity", product_ids, ev6, ">95:5", ">95:5", basis="not reported", scope=f"E:Z ratio; source labels {major_geom}-{label} as major")
    classify(r, ev6, specs)


# Figure 4A formaldehyde allylboration.
f4a = new_rx("figure4A_formaldehyde_9a", ev6, categories=["derivatization"], location="Figure 4A",
             name="Allylboration of formaldehyde", transformation="allylboration of formaldehyde",
             unresolved=["Formaldehyde polymer composition, workup, purification, and yield measurement basis are not reported on page 6."])
addp(f4a, "starting_4hc", prod4hc[0], "reactant", ev6)
addp(f4a, "formaldehyde", M["paraformaldehyde"], "reactant", ev6, equivalents=("5 equiv", q(5, "equiv")), reference=prod4hc[0])
addp(f4a, "solvent", M["THF"], "solvent", ev6)
p9a = mol(f"{DOC}:mol_9a", "9a", ev6, chirality="unspecified", geometry="unspecified", representation="not_encoded",
          identity_note="Homoallylic alcohol 9a; the page does not provide an absolute/relative stereochemical descriptor.")
addp(f4a, "product", p9a, "product", ev6, product_role="homoallylic alcohol")
f4a["solvent_conditions"]["composition"] = f("THF", ["THF"], ev6)
add_step(f4a, "reaction", ev6, operation="heat allylic boronate with paraformaldehyde", temperature=("65 °C", q(65, "°C")), time=("48 h", q(48, "h")), participant_ids=[p["participant_id"] for p in f4a["participants"]])
outcome(f4a, "yield", "yield", [p9a], ev6, "90%", q(90, "%"), role="9a", yield_type="reported yield", scope="yield of 9a")
classify(f4a, ev6, ALLYLBORATION)


# Figure 4B generic scheme and three specific aryl-aldehyde products.
f4generic = new_rx("figure4B_generic_aryl_aldehyde_allylboration", ev6, record_type="generic_scheme",
                   categories=["derivatization"], location="Figure 4B generic arrow",
                   name="Allylboration of aryl aldehydes", transformation="Lewis-acid-mediated allylboration",
                   unresolved=["R groups are generic; product absolute configurations are not named on page 6."])
g4 = mol(f"{DOC}:mol_generic_E4", "generic (E)-4 allylic boronate", ev6, chirality="unspecified", geometry="specified", representation="not_encoded", stereo_note="E/Z >95:5")
gald = mol(f"{DOC}:mol_generic_aryl_aldehyde", "generic aryl aldehyde RCHO", ev6, chirality="achiral", geometry="not_applicable", representation="not_encoded")
g9 = mol(f"{DOC}:mol_generic_9", "generic homoallylic alcohol 9", ev6, chirality="partially_specified", geometry="unspecified", representation="not_encoded", stereo_note="relative stereochemistry drawn but no absolute descriptors are printed")
addp(f4generic, "allylic_boronate", g4, "reactant", ev6)
addp(f4generic, "aldehyde", gald, "reactant", ev6, equivalents=("1.5 equiv", q(1.5, "equiv")), reference=g4)
addp(f4generic, "BF3OEt2", M["BF3OEt2"], "reagent", ev6, equivalents=("1.0 equiv", q(1.0, "equiv")), reference=g4)
addp(f4generic, "solvent", M["DCM"], "solvent", ev6)
addp(f4generic, "product", g9, "product", ev6, product_role="generic homoallylic alcohol")
f4generic["solvent_conditions"]["composition"] = f("DCM", ["DCM"], ev6)
add_step(f4generic, "reaction", ev6, operation="BF3·OEt2-mediated allylboration", temperature=("−78 °C", q(-78, "°C")), time=("3 h", q(3, "h")), participant_ids=[p["participant_id"] for p in f4generic["participants"]])
classify(f4generic, ev6, ALLYLBORATION)

for code, aldehyde, y in [("9b", "benzaldehyde", 75), ("9c", "3bromo4fluorobenzaldehyde", 60), ("9d", "2naphthaldehyde", 52)]:
    r = new_rx(f"figure4B_{code}", ev6, categories=["derivatization"], location=f"Figure 4B, {code}",
               name="Allylboration of aryl aldehydes", transformation="Lewis-acid-mediated allylboration",
               unresolved=["Exact allylic-boronate substrate code, absolute configurations, workup, purification, and analytical methods are not printed on page 6."])
    addp(r, "allylic_boronate", g4, "reactant", ev6)
    addp(r, "aldehyde", M[aldehyde], "reactant", ev6, equivalents=("1.5 equiv", q(1.5, "equiv")), reference=g4)
    addp(r, "BF3OEt2", M["BF3OEt2"], "reagent", ev6, equivalents=("1.0 equiv", q(1.0, "equiv")), reference=g4)
    addp(r, "solvent", M["DCM"], "solvent", ev6)
    major = mol(f"{DOC}:mol_{code}_major", f"{code} major diastereomer", ev6, chirality="partially_specified", geometry="unspecified", representation="not_encoded",
                stereo_note="relative configuration drawn; no absolute R/S descriptors printed" + ("; single-crystal X-ray diffraction confirmed the depicted stereochemistry" if code == "9d" else ""))
    minor = mol(f"{DOC}:mol_{code}_minor", f"{code} minor diastereomer", ev6, chirality="unresolved", geometry="unspecified", representation="not_encoded",
                identity_note="Minor diastereomer implied by dr; its structure/configuration is not individually assigned on page 6.")
    addp(r, "major_product", major, "product", ev6, product_role="major diastereomer")
    addp(r, "minor_product", minor, "product", ev6, product_role="minor diastereomer")
    r["solvent_conditions"]["composition"] = f("DCM", ["DCM"], ev6)
    add_step(r, "reaction", ev6, operation="BF3·OEt2-mediated allylboration", temperature=("−78 °C", q(-78, "°C")), time=("3 h", q(3, "h")), participant_ids=[p["participant_id"] for p in r["participants"]])
    outcome(r, "yield", "yield", [major, minor], ev6, f"{y}%", q(y, "%"), role=code, yield_type="reported yield", scope=f"yield of {code} diastereomeric product fraction")
    outcome(r, "dr", "selectivity", [major, minor], ev6, ">95:5", ">95:5", basis="not reported", scope=f"diastereomeric ratio for {code}; depicted major:minor")
    r["general_scheme_ids"] = [f4generic["reaction_id"]]
    classify(r, ev6, ALLYLBORATION)


# Figure 4D is a proposed transition-state model, not an observed synthesis.
mech = new_rx("figure4D_proposed_diastereoselectivity_model", ev6, record_type="proposed_mechanism",
              categories=["mechanistic_study"], location="Figure 4D",
              name="Proposed diastereoselectivity-determining mechanism",
              transformation="proposed six-membered allylboration transition state",
              unresolved=["Proposed model is not represented as an experimentally verified elementary reaction."])
addp(mech, "generic_allylic_boronate", g4, "reactant", ev6)
addp(mech, "generic_aryl_aldehyde", gald, "reactant", ev6)
addp(mech, "generic_product", g9, "product", ev6, product_role="proposed product relationship")
add_step(mech, "proposed_transition_state", ev6, operation="author-proposed low-steric-repulsion transition-state arrangement", participant_ids=[p["participant_id"] for p in mech["participants"]])
classify(mech, ev6, coverage="unresolved", limitation="Author-proposed stereochemical model; no observed bond-change experiment is asserted here.")


# Context that is not converted into current-paper experiments.
bundle["source_context"] = [
    {"context_id": f"{DOC}:bibliographic_title", "document_id": DOC, "kind": "bibliographic_metadata",
     "title": TITLE, "doi": "10.1021/jacs.1c06538", "evidence_ids": [E[1], T[1]]},
    {"context_id": f"{DOC}:figure1A_conceptual_utility", "document_id": DOC, "kind": "conceptual_scheme",
     "description": "Potential downstream homoallylic-alcohol and allylic-amine utility; no specific conditions/yields on page 2.", "evidence_ids": [E[2]]},
    {"context_id": f"{DOC}:figure2_literature_overview", "document_id": DOC, "kind": "literature_overview",
     "description": "Figure 2A/B summarizes regio-divergence and cited prior carboboration classes; individual cited papers are not current-work experiments.", "evidence_ids": [E[3]]},
    {"context_id": f"{DOC}:figure4C_xray", "document_id": DOC, "kind": "structure_assignment",
     "description": "Single-crystal X-ray structure of allylboration product 9d confirms the depicted stereochemistry.", "evidence_ids": [E[6], T[6]]},
]


# Raw visual transcription and exact page coverage.
table2_ids = [f"{DOC}:table2_{row[0]}" for row in scope_rows]
failed_ids = [f"{DOC}:table2_unsuitable_{x}" for x in ("2k", "2l", "2m", "2n")]
table1_ids = [f"{DOC}:table1_entry_{i:02d}" for i in range(1, 19)]
page6_ids = [
    f"{DOC}:figure3A_gram_scale_4hc", f"{DOC}:figure3B_oxidation_6hc",
    f"{DOC}:figure3B_amination_7hc", f"{DOC}:figure3B_homologation_8hc",
    f"{DOC}:figure4A_formaldehyde_9a", f"{DOC}:figure4B_generic_aryl_aldehyde_allylboration",
    f"{DOC}:figure4B_9b", f"{DOC}:figure4B_9c", f"{DOC}:figure4B_9d",
    f"{DOC}:figure4D_proposed_diastereoselectivity_model",
]
coverage = {
    "document_id": DOC, "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title, abstract, and graphical abstract", "finding": "Current-work generic three-component carboboration; bibliographic title/DOI captured; no specific experiment row.", "evidence_ids": [E[1], T[1]]},
            {"region": "introduction columns", "finding": "Background prose and cited prior work; no page-1 specific experiment.", "evidence_ids": [E[1], T[1]]}],
         "record_ids": [generic["reaction_id"]],
         "unresolved": ["Graphical-abstract R groups and exact conditions are generic."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Figure 1A", "finding": "Conceptual downstream utility, retained in source_context rather than counted as experiments.", "evidence_ids": [E[2]]},
            {"region": "Figure 1B", "finding": "Explicitly cited Hall literature example; recorded as literature_example.", "evidence_ids": [E[2]]},
            {"region": "Results and Discussion prose", "finding": "Narrative reports Table 1 entries 1-9; linked to those row records without duplicate experiments.", "evidence_ids": [E[2], T[2]]}],
         "record_ids": [hall["reaction_id"]] + table1_ids[:9],
         "unresolved": ["Figure 1A has no experimental conditions or outcomes."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Figure 2A/B", "finding": "Conceptual/literature overview, not current-paper specific experiments.", "evidence_ids": [E[3]]},
            {"region": "Figure 2C", "finding": "Current-work generic mechanistic/regioselectivity design; merged with graphical abstract generic record.", "evidence_ids": [E[3]]},
            {"region": "bottom prose", "finding": "Narrative reports Table 1 entries 10-16; linked to row records.", "evidence_ids": [E[3], T[3]]}],
         "record_ids": [generic["reaction_id"]] + table1_ids[9:16],
         "unresolved": ["Figure 2 contains generic structures and cited examples without row-level conditions."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "Table 1 scheme, entries 1-18, and footnotes a-f", "finding": "All 18 optimization rows recorded independently; entry 17 minor isomer remains unidentified; entry 18 links duplicate Table 2 presentation.", "evidence_ids": [E[4], T[4]]},
            {"region": "lower substrate-scope prose", "finding": "Narrative summary used as context; row-specific Table 2 values are taken from page 5.", "evidence_ids": [E[4], T[4]]}],
         "record_ids": table1_ids,
         "unresolved": ["2 mol% versus 0.025 mmol catalyst conflict in applicable rows.", "Entry 17 minor-isomer structure is not identified."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "Table 2 conditions A/B and scope of allenes", "finding": "Twenty-one allene-scope panels inventoried; 4ba is duplicate evidence for Table 1 entry 18; 4kc is a yielded but stereochemically unsuccessful case.", "evidence_ids": [E[5], T[5]]},
            {"region": "Table 2 scope of alkyl halides", "finding": "Eight successful halide-scope arrows and four unsuccessful/decomposition entries recorded independently.", "evidence_ids": [E[5], T[5]]},
            {"region": "Table 2 footnotes a-k", "finding": "Conditions, oxidation-derived yields, modified temperature/time/electrophile, GC/NMR bases, selective oxidation, and 86:14 I/Cl selectivity structured.", "evidence_ids": [E[5], T[5]]}],
         "record_ids": table2_ids + failed_ids + [f"{DOC}:table1_entry_18"],
         "unresolved": ["Dense scope structures lack safe complete SMILES.", "Oxidation reagents/conditions for footnote-d yields are absent without SI.", "Unidentified side products for 4bi/4bj remain unresolved."]},
        {"page": 6, "state": "attempted", "regions": [
            {"region": "Figure 3A", "finding": "One gram-scale carboboration experiment recorded with 6.0 mmol scale, 1.86 g/88% outcome and selectivities.", "evidence_ids": [E[6], T[6]]},
            {"region": "Figure 3B", "finding": "Three separate derivatization arrows recorded, including ordered two-stage amination/Boc protection.", "evidence_ids": [E[6], T[6]]},
            {"region": "Figure 4A/B", "finding": "One formaldehyde and three aryl-aldehyde allylborations recorded; Figure 4B generic arrow retained separately.", "evidence_ids": [E[6], T[6]]},
            {"region": "Figure 4C/D", "finding": "9d X-ray assignment retained as source_context; proposed transition-state model recorded as proposed_mechanism, not experiment.", "evidence_ids": [E[6], T[6]]},
            {"region": "bottom prose referencing Figure 5", "finding": "Figure 5 experiments occur on page 7 and are outside this session's claimed coverage; no page-7 record created.", "evidence_ids": [E[6], T[6]]}],
         "record_ids": page6_ids,
         "unresolved": ["Figure 3/4 yield/selectivity analytical bases are not reported on page 6.", "Absolute configurations of 9b-9d are not printed; only depicted relative stereochemistry is retained."]},
    ],
}

visual = {
    "document_id": DOC, "title": TITLE, "pages_inspected": [1, 2, 3, 4, 5, 6],
    "page_images": [f"assets/{DOC}/page{p}.jpg" for p in range(1, 7)],
    "table1_rows": [list(row) for row in table1] + [[17, "1b", "SIMesCuCl", None, 2, -5, None, 58, None, "77:23 unidentified mixture", None], [18, "1b", "CuCl", "Xantphos", 5, 30, None, 83, None, ">95:5", ">99:1"]],
    "table2_success_rows": [list(row) for row in scope_rows] + [["4ba", "1b", "2a", "B", 83, ">95:5", ">95:5", "duplicate of Table 1 entry 18"]],
    "table2_unsuccessful_rows": [["2k", "decomposition of halide"], ["2l", "decomposition of halide; styrene side product in prose"], ["2m", "complex mixture"], ["2n", "no reaction"]],
    "page6_arrows": ["Figure 3A gram-scale 4hc", "Figure 3B 6hc", "Figure 3B 7hc", "Figure 3B 8hc", "Figure 4A 9a", "Figure 4B 9b", "Figure 4B 9c", "Figure 4B 9d"],
    "non_experiment_regions": ["Figure 1A conceptual utility", "Figure 2A/B literature overview", "Figure 4C X-ray image", "Figure 4D proposed model"],
}

bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
st.write(RUN / "submitted.json", bundle)
st.write(RUN / "coverage/pages_0001_0006.json", coverage)
st.write(RUN / "raw_model/visual_transcription.json", visual)
st.write(RUN / "raw_model/source_mentions.json", st.source_mentions(dataset, ["1a", "2a", "4aa", "4hc", "9d"]))
st.write(RUN / "raw_model/stereo_diagnostics.json", {
    "1a": st.stereo_inspect("C=C=C(C)C1CCCCC1"),
    "1b": st.stereo_inspect("C=C=C(C)C(C)(C)C"),
    "ethyl_iodide": st.stereo_inspect("CCI"),
    "benzyl_bromide": st.stereo_inspect("BrCc1ccccc1"),
    "limitations": "Graph diagnostics only; unresolved product E/Z and relative/central stereochemistry remain source-text facts, not encoded structures.",
})
st.write(RUN / "raw_model/premerge_counts.json", {
    "molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
    "reaction_types": dict(Counter(r["record_type"] for r in bundle["reactions"])),
    "needs_review": sum(r["review"]["status"] == "needs_review" for r in bundle["reactions"]),
})

workflow = f"""# Extraction workflow: J. Am. Chem. Soc. 2021, 143, 13865–13877, PDF pages 1–6

- Batch: `allene-stereo-v31-200-20261017`
- Document: `{DOC}`
- Title: *{TITLE}*
- DOI: `10.1021/jacs.1c06538`
- Claimed visual scope: original PDF pages 1–6 exactly. All six supplied full-page JPEGs were opened at original detail. No page outside 1–6 is claimed as covered.
- Source association: main article only. Supporting Information is explicitly mentioned by the article but no SI file was supplied or associated; `source_mapping.json` remains unestablished and no SI was inferred.
- Review state: all imported molecules/reactions remain machine records; every reaction is `needs_review`, with semantic/chemical accuracy flags false.

## Source handling and visual inventory

1. Read `AGENTS.md`, the complete extraction skill, `extraction_standard.md`, `literature.schema.json`, `stereochemistry.md`, `structured_fields.md`, `quality_rules.md`, and `bond_taxonomy.json` before authoring records.
2. Confirmed the mechanically prepared dataset had zero candidates and zero experiments. The pre-existing `runs/workflow.md` described another JACS paper and was treated as stale output, never as source data; this file replaces it.
3. Opened page images 1–6 at original detail. Lineage-tracked crops and receipts for Table 1, Table 2, and Figures 3–4 are under `runs/current/raw_model/`.
4. Page 1: title/abstract and a current-work generic carboboration graphical abstract. Page 2: Figure 1A conceptual utility, Figure 1B one cited Hall literature example, and prose for Table 1 entries 1–9. Page 3: Figure 2 conceptual/literature/current-work generic scheme and prose for entries 10–16. Page 4: all 18 Table 1 rows and footnotes. Page 5: 29 successful Table 2 product panels plus four unsuccessful halide panels; the 4ba panel exactly repeats Table 1 entry 18 and is merged as duplicate evidence. Page 6: one gram-scale arrow, three Figure 3B derivatizations, one formaldehyde and three aryl-aldehyde allylborations, an X-ray structure assignment, and one proposed stereochemical model.
5. Figure 5 is only referenced at the bottom of page 6; its image/arrows occur on page 7. No Figure 5 experiment is claimed or created in this page-1–6 session.

## Record construction decisions

- Table 1 optimization rows are independent experiments. Parenthesized `(4xa+5xa)` yield is a second outcome, not folded into the 4xa yield. E/Z and 4/5 are distinct targeted selectivity outcomes. GC, 13C NMR, 1H NMR, isolated-yield, and missing methods remain separate.
- Table 1 footnote a states Cu(I) catalyst `0.025 mmol`, while 2 mol% rows on a 0.5 mmol allene basis imply a different charge. Both printed values are retained, with a note and `review.unresolved`; no value is silently corrected. For entry 16, 0.25 M supersedes standard 0.5 M, but the unprinted adjusted solvent volume is not inferred.
- Table 1 entry 17 has a 77:23 mixture whose minor structure is explicitly unidentified. It is not forced into a definitive E/Z or 4/5 assignment. Reported protoboration and boryl-substitution contaminants are retained in notes.
- Table 2 conditions A/B include all printed mmol, equivalents, mol%, DMF volume/concentration, temperature, and time. Footnotes e/h/i are applied only to the indicated entries. Footnote j is represented as benzyl bromide. Footnote k is a separate 86:14 iodide:chloride chemoselectivity outcome.
- Footnote-d Table 2 yields target an explicit unresolved oxidized isolation derivative, with a separate oxidation step whose reagents/temperature/time remain unreported. This avoids pretending that a directly isolated boryl product was measured. Footnote-f selectivities use GC after silica; other Table 2 selectivities use 1H NMR after column chromatography.
- The four unsuitable halides retain decomposition/complex-mixture/no-reaction outcomes and no invented product. Phenethyl iodide retains styrene as a reported side product. Failed records are classification `not_applicable`, not assigned expected-product bond labels.
- Figure 3B amination is two ordered stages: H2NOMe/n-BuLi at −78 to 60 °C for 48 h, then Boc2O at rt for 2 h. The source does not resolve ramp/hold timing within the first 48 h. Other page-6 temperature programs are retained verbatim.
- Figure 4B's generic arrow is distinct from the three specific products. 9b–9d have separate major/minor diastereomer identities; dr is not converted to ee/er or absolute configuration. Figure 4C X-ray evidence is source context for 9d. Figure 4D is a `proposed_mechanism`, not an observed experiment.

## Identity, stereochemistry, and classification

- Every molecule has explicit `chirality_status`, `geometry_status`, and `representation_status`. E and Z products are separate molecule IDs where the metric explicitly establishes both components. Unidentified/complex scope structures remain without complete SMILES; missing stereochemical marks are never converted into racemic/achiral claims.
- Candidate SMILES are limited to unambiguous small reagents and the two model allenes. `standard stereo` graph diagnostics are saved in `raw_model/stereo_diagnostics.json`; parsing does not elevate identity to `structure_checked`. Product alkene geometry and 9b–9d relative stereochemistry remain textual/`not_encoded`.
- Successful carboborations have partial C–C and C–B formation labels plus the source-visible allene bond-order decrease. Oxidation-isolation sequences also retain C–O/O–H formation and C–B cleavage. Allylborations retain C–C, promoted C=C, carbonyl reduction, O–H formation, and C–B cleavage. No full atom mapping or three-model votes are fabricated.
- CAS numbers were not reported and none were invented. No external property lookup was performed.

## Expected unresolved audit items

- Page-6 Figure 3/4 yields and selectivities do not print isolated/NMR/GC bases; these remain missing rather than inherited from Table 2.
- Dense allene/product drawings and catalyst complexes are not converted into unsafe complete SMILES. The labels, roles, quantities, stereochemical text, outcomes, and evidence remain available.
- Table 2 oxidation reagents/conditions and unidentified side products are unavailable without SI.
- Validation establishes schema integrity only, not chemical correctness, semantic completeness, or reproducibility. `semantic_accuracy_verified` and `chemical_accuracy_verified` remain false.

## Outputs

- Submission/merged v3 data: `runs/current/submitted.json`, `dataset.json`, `standard_records.json`, `molecules.json`, `records/`, `classification_index.json`, and `规范记录预览.md`.
- Required coverage: `runs/current/coverage/pages_0001_0006.json`.
- Raw artifacts: visual transcription, source-mention hints, stereo diagnostics, lineage crops/receipts, build script, and count summaries under `runs/current/raw_model/`.
- Deterministic audit: `runs/current/audit.json`.

All chemistry-tool commands use `/EvoLab/student-workspaces/xllin/OpenEvo-chem-integration/.venv/bin/python`. Machine records remain `needs_review`.
"""
(ROOT / "runs/workflow.md").write_text(workflow)


def finalize():
    merged = json.loads((RUN / "dataset.json").read_text())
    audit = dataset_audit(merged)
    st.write(RUN / "audit.json", audit)
    standard = merged.get("standard", {})
    counts = {
        "documents": len(merged.get("documents", [])),
        "molecules": len(standard.get("molecules", [])),
        "reactions": len(standard.get("reactions", [])),
        "reaction_types": dict(Counter(r["record_type"] for r in standard.get("reactions", []))),
        "needs_review": sum(r["review"]["status"] == "needs_review" for r in standard.get("reactions", [])),
        "warning_records": len(audit.get("warnings", [])),
        "warning_categories": dict(Counter(w.split(":", 1)[0] for row in audit.get("warnings", []) for w in row["warnings"])),
        "semantic_accuracy_verified": False,
        "chemical_accuracy_verified": False,
    }
    st.write(RUN / "raw_model/final_counts.json", counts)
    print(json.dumps(counts, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--finalize", action="store_true")
    args = parser.parse_args()
    if args.finalize:
        finalize()
    else:
        print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"])}, ensure_ascii=False))
