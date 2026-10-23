#!/usr/bin/env python3
"""Extend the source-derived v3.1 bundle with original PDF pages 7--12.

The script is deliberately evidence preserving.  It updates page-3 summary
records when pages 10--12 provide the full procedure, creates independent
records for distinct arrows/table experiments, and keeps mechanisms and cited
literature separate from experiments performed in the current paper.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SKILL = ROOT / "skills" / "chemical-reaction-extraction"
sys.path.insert(0, str(SKILL / "scripts"))

from standard_records import (  # noqa: E402
    bond_label,
    fact,
    new_molecule,
    new_outcome,
    new_participant,
    new_reaction,
    new_step,
    stereo_inspect,
)

DOC = "doc_2401da5b5fbe4542"
T = {p: f"{DOC}:page{p}:text" for p in range(7, 13)}
I = {p: f"{DOC}:page{p}:image" for p in range(7, 13)}


def ev(*pages: int) -> list[str]:
    out: list[str] = []
    for page in pages:
        out.extend([I[page], T[page]])
    return out


def F(raw=None, normalized=None, evidence=(), method="original_extraction", *,
      status=None, verification=None, derivation=None, alternatives=()):
    if normalized is None and raw is not None:
        normalized = raw
    return fact(raw, normalized, evidence, method, status, verification,
                derivation, alternatives)


def missing(evidence=(), reason=None, *, status="unreported"):
    return fact(evidence=evidence, derivation=reason, status=status)


def q(number, unit, **extra):
    result = {"number": number, "unit": unit}
    result.update(extra)
    return result


bundle_path = ROOT / "runs/current/submitted.json"
bundle = json.loads(bundle_path.read_text())
dataset = json.loads((ROOT / "runs/current/dataset.json").read_text())
# Standard init/prepare may add image evidence after an earlier page batch.
# Always carry the authoritative document/evidence registry into the bundle.
bundle["documents"] = dataset["documents"]
molecules = {m["molecule_id"]: m for m in bundle["molecules"]}
reactions = {r["reaction_id"]: r for r in bundle["reactions"]}


def add_occurrence(m, label, evidence):
    if any(o["compound_label"]["normalized_value"] == label and
           set(o["evidence_ids"]) == set(evidence) for o in m["occurrences"]):
        return
    m["occurrences"].append({
        "document_id": DOC,
        "compound_label": F(label, label, evidence),
        "evidence_ids": list(evidence),
    })


def add_molecule(mid, label, name, evidence, *, smiles=None,
                 chirality="unresolved", geometry="unresolved",
                 representation=None, assignment=None, mixture=None,
                 axial=None, note=None):
    if mid in molecules:
        add_occurrence(molecules[mid], label, evidence)
        return mid
    m = new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, evidence)
    m["names"]["original"] = F(label, label, evidence)
    m["names"]["english"] = F(name, name, evidence)
    m["names"]["abbreviations"] = F([label], [label], evidence)
    m["names"]["aliases"] = missing(evidence, "No additional source alias is reported on pages 7--12.")
    m["cas_numbers"] = missing(evidence, "CAS number is not reported on pages 7--12 and was not guessed.")
    m["generic_template"] = None
    if smiles:
        diagnostic = stereo_inspect(smiles)
        canonical = diagnostic["canonical_isomeric_smiles"]
        m["representation_type"] = "smiles"
        m["smiles"] = F(smiles, canonical, evidence, "structure_derivation",
                          derivation="Connectivity was conservatively encoded from the displayed structure or unambiguous source name; toolkit parsing is not identity verification.")
        if representation is None:
            representation = "partial" if any(
                item["specified"] == "Unspecified"
                for item in diagnostic["potential_stereo"]
            ) else "encoded"
        m["stereochemistry"]["isomeric_smiles"] = F(
            canonical, canonical, evidence, "structure_derivation",
            derivation="Canonical isomeric SMILES generated from the conservatively transcribed source graph; source identity remains needs_review.")
        m["stereochemistry"]["identity_key"] = F(
            diagnostic["stereo_graph_key"], diagnostic["stereo_graph_key"], evidence,
            "structure_derivation", derivation="Toolkit graph key preserves the encoded isomeric graph and does not establish sample identity.")
        m["stereochemistry"]["atom_centers"] = F(
            diagnostic["atom_centers"], diagnostic["atom_centers"], evidence,
            "structure_derivation", derivation="Standard stereo graph inspection; source configuration is not inferred from an untagged graph.")
        m["stereochemistry"]["double_bonds"] = F(
            diagnostic["double_bonds"], diagnostic["double_bonds"], evidence,
            "structure_derivation", derivation="Standard stereo graph inspection; E/Z is not inferred where the source is silent.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing(evidence, "Exact source structure or stereochemical state was not safely encoded; label and drawing evidence are retained.")
        m["stereochemistry"]["isomeric_smiles"] = missing(evidence)
        m["stereochemistry"]["identity_key"] = missing(evidence)
        m["stereochemistry"]["atom_centers"] = F(
            chirality, chirality, evidence, "structure_derivation",
            derivation="Original page drawing inspected; no absolute configuration was invented.")
        m["stereochemistry"]["double_bonds"] = F(
            geometry, geometry, evidence, "structure_derivation",
            derivation="Original page drawing inspected; no E/Z or ring relationship was invented.")
        representation = representation or "not_encoded"
    m["stereochemistry"]["chirality_status"] = F(
        chirality, chirality, evidence, "structure_derivation",
        derivation="Assigned only to the level supported by the original page drawing/name; unspecified is not treated as achiral or racemic.")
    m["stereochemistry"]["geometry_status"] = F(
        geometry, geometry, evidence, "structure_derivation",
        derivation="Assigned only to the level supported by the original page drawing/name.")
    m["stereochemistry"]["representation_status"] = F(
        representation, representation, evidence, "inference",
        derivation="Compares the machine representation with the stereochemical content visible in the source.")
    m["stereochemistry"]["assignment"] = (
        F(assignment, assignment, evidence) if assignment else
        missing(evidence, "No absolute or relative stereochemical assignment is reported for this identity on pages 7--12.")
    )
    m["stereochemistry"]["mixture"] = (
        F(mixture, mixture, evidence) if mixture else
        missing(evidence, "No stereoisomer mixture statement is reported for this identity on pages 7--12.")
    )
    m["stereochemistry"]["axial"] = (
        F(axial, axial, evidence, "structure_derivation",
          derivation="Allene/axis state was inspected in the source drawing.") if axial else missing(evidence)
    )
    m["identity_note"] = F(note or name, note or name, evidence)
    molecules[mid] = m
    bundle["molecules"].append(m)
    return mid


def charge(p, field, raw, number, unit, evidence):
    p["charges"][field] = F(raw, q(number, unit), evidence)


def participant(r, mid, role, evidence, *, stage="reaction", mass=None,
                amount=None, volume=None, equiv=None, molpct=None,
                product_role=None):
    pid = f"{r['reaction_id']}:{mid}:{len(r['participants']) + 1}"
    p = new_participant(pid, mid, role, evidence, stage=stage)
    if mass:
        charge(p, "mass", mass[0], mass[1], mass[2], evidence)
    if amount:
        charge(p, "amount", amount[0], amount[1], amount[2], evidence)
    if volume:
        charge(p, "volume", volume[0], volume[1], volume[2], evidence)
    if equiv:
        charge(p, "equivalents", equiv[0], equiv[1], "equiv", evidence)
    if molpct:
        charge(p, "mol_percent", molpct[0], molpct[1], "mol%", evidence)
    if equiv or molpct:
        ref = next((x["molecule_id"] for x in r["participants"] if x["role"] == "reactant"), mid)
        p["equivalence_reference"] = F(ref, ref, evidence)
    if product_role:
        p["product_role"] = F(product_role, product_role, evidence)
    r["participants"].append(p)
    return pid


def step(r, operation, evidence, *, temperature=None, time=None,
         atmosphere=None, pressure=None, stage="reaction", participants=()):
    s = new_step(f"{r['reaction_id']}:step{len(r['steps']) + 1}", evidence, stage=stage)
    s["participant_ids"] = list(participants)
    s["operation"] = F(operation, operation, evidence)
    if temperature is not None:
        if isinstance(temperature, tuple):
            raw, value = temperature
            norm = q(value, "°C") if isinstance(value, (int, float)) else value
            s["temperature"] = F(raw, norm, evidence)
        else:
            norm = q(temperature, "°C") if isinstance(temperature, (int, float)) else temperature
            s["temperature"] = F(str(temperature), norm, evidence)
    if time is not None:
        if isinstance(time, tuple) and len(time) == 3 and time[0] == "range":
            s["time"] = F(f"{time[1]}-{time[2]} h", {"minimum": time[1], "maximum": time[2], "unit": "h"}, evidence)
        elif isinstance(time, tuple):
            raw, value, unit = time
            s["time"] = F(raw, q(value, unit), evidence)
        else:
            s["time"] = F(str(time), time, evidence)
    if atmosphere:
        s["atmosphere"] = F(atmosphere, atmosphere, evidence)
    if pressure:
        s["pressure"] = F(pressure, pressure, evidence)
    r["steps"].append(s)
    return s["step_id"]


def outcome(r, kind, targets, raw, normalized, evidence, *, yield_type=None,
            basis=None, scope=None, product_role=None):
    o = new_outcome(f"{r['reaction_id']}:out{len(r['outcomes']) + 1}", kind, targets, evidence)
    o["value"] = F(raw, normalized, evidence)
    if yield_type:
        o["yield_type"] = F(yield_type, yield_type, evidence)
    else:
        o["yield_type"] = missing(evidence, "Yield type is not applicable or not reported for this outcome.")
    if basis:
        o["measurement_basis"] = F(basis, basis, evidence)
    else:
        o["measurement_basis"] = missing(evidence, "Measurement method/basis is not reported in the source region.")
    o["scope"] = F(scope or "reaction outcome", scope or "reaction outcome", evidence)
    if product_role:
        o["product_role"] = F(product_role, product_role, evidence)
    r["outcomes"].append(o)
    return o["outcome_id"]


def set_classification(r, labels, evidence, limitation):
    changes = []
    normalized = []
    for i, code in enumerate(labels, 1):
        parts = code.split("-")
        if code.startswith("DB-") or code.startswith("TB-"):
            parts = parts[1:]
        a, b = parts[-2], parts[-1]
        order = 2 if code.startswith("DB-") else 3 if code.startswith("TB-") else 1
        item = bond_label(a, b, order)
        change_id = f"b{i}"
        changes.append(F(
            {"change_id": change_id, "atoms": [{"atom_map": i * 2 - 1, "element": a}, {"atom_map": i * 2, "element": b}],
             "before_order": 0, "after_order": order, "change_type": "formed"},
            {"change_id": change_id, "atoms": [{"atom_map": i * 2 - 1, "element": a}, {"atom_map": i * 2, "element": b}],
             "before_order": 0, "after_order": order, "change_type": "formed"},
            evidence, "inference",
            derivation="Placeholder atom maps encode only the source-inferred net newly formed bond class; no complete atom mapping is asserted."))
        normalized.append({"category_number": item["category_number"], "code": item["code"], "supporting_change_ids": [change_id]})
    r["classification"]["bond_changes"] = changes
    r["classification"]["labels"] = F(labels, normalized, evidence, "inference",
                                         derivation="Net newly formed bond classes inferred from the original source drawing; unchanged bonds and bond-order decreases are not classified as new bonds.")
    r["classification"]["coverage"] = F("partial", "partial", evidence, "inference",
                                           derivation="Only confidently visible net newly formed bond classes are retained; full atom mapping is not claimed.")
    r["classification"]["limitations"] = [limitation]


def reaction(short_id, evidence, record_type, categories, name, transformation,
             location, labels=(), unresolved=()):
    rid = f"{DOC}:{short_id}"
    if rid in reactions:
        raise ValueError(f"Refusing to duplicate existing reaction {rid}")
    r = new_reaction(rid, DOC, evidence)
    r["record_type"] = record_type
    category_map = {
        "proposed_mechanism": "other",
        "cited_literature": "other",
        "substrate_preparation": "other",
        "application_synthesis": "derivatization",
        "scope": "substrate_scope",
    }
    normalized_categories = [category_map.get(category, category) for category in categories]
    r["experiment_categories"] = F(categories, normalized_categories, evidence)
    r["location"] = F(location, location, evidence)
    r["reaction_name"] = F(name, name, evidence)
    r["transformation_type"] = F(transformation, transformation, evidence)
    r["review"] = {
        "status": "needs_review", "method": "assistant_visual", "reviewer": None,
        "checked_fields": ["original page image", "compound labels", "arrow/table row", "conditions", "yield/footnotes"],
        "unresolved": list(unresolved) or ["Complete atom mapping and chemical identity remain unverified."],
    }
    if labels:
        set_classification(r, labels, evidence, "No complete atom-mapped reaction SMILES; classification is an evidence-linked inference from the displayed arrow.")
    reactions[rid] = r
    bundle["reactions"].append(r)
    return r


def add_simple(mid, label, name, pages, smiles, *, chirality="achiral", geometry="not_applicable", note=None):
    return add_molecule(mid, label, name, ev(*pages), smiles=smiles,
                        chirality=chirality, geometry=geometry, note=note)


# Reagents, solvents, substrates and products newly visible on pages 7--12.
simple = [
    ("mol_benzoyl_acetonitrile", "19", "benzoyl acetonitrile", (7,), "N#CCC(=O)c1ccccc1"),
    ("mol_tfa", "TFA", "trifluoroacetic acid", (8,), "O=C(O)C(F)(F)F"),
    ("mol_ethyl_chloroformate", "ClCO2Et", "ethyl chloroformate", (8,), "CCOC(=O)Cl"),
    ("mol_nabh4", "NaBH4", "sodium borohydride", (8,), "[Na+].[BH4-]"),
    ("mol_mg", "Mg", "magnesium", (8,), "[Mg]"),
    ("mol_meoh", "MeOH", "methanol", (8, 9), "CO"),
    ("mol_propargyl_alcohol", "S1", "propargyl alcohol", (9,), "OCC#C"),
    ("mol_dhp", "3,4-dihydro-2H-pyran", "3,4-dihydro-2H-pyran", (9,), "C1=COCCC1"),
    ("mol_ptsoh", "p-TsOH", "p-toluenesulfonic acid", (9,), "Cc1ccc(S(=O)(=O)O)cc1"),
    ("mol_thf", "THF", "tetrahydrofuran", (9, 10, 11), "C1CCOC1"),
    ("mol_nbuli", "n-BuLi", "n-butyllithium", (9,), "[Li]CCCC"),
    ("mol_boc2o", "Boc2O", "di-tert-butyl dicarbonate", (9,), "CC(C)(C)OC(=O)OC(=O)OC(C)(C)C"),
    ("mol_nh4cl", "NH4Cl", "ammonium chloride", (9, 10, 11), "[NH4+].[Cl-]"),
    ("mol_etoac", "EtOAc", "ethyl acetate", (9, 10, 11, 12), "CCOC(C)=O"),
    ("mol_nacl", "NaCl", "sodium chloride", (9, 10, 11), "[Na+].[Cl-]"),
    ("mol_na2so4", "Na2SO4", "sodium sulfate", (9, 10, 11), "[Na+].[Na+].[O-]S(=O)(=O)[O-]"),
    ("mol_etoh", "EtOH", "ethanol", (9,), "CCO"),
    ("mol_ppts", "PPTS", "pyridinium p-toluenesulfonate", (9,), "Cc1ccc(S(=O)(=O)[O-])cc1.c1cc[nH+]cc1"),
    ("mol_nahco3", "NaHCO3", "sodium bicarbonate", (9, 10), "[Na+].O=C([O-])O"),
    ("mol_et2o", "Et2O", "diethyl ether", (9, 11), "CCOCC"),
    ("mol_s6", "S6", "tert-butyl 4-hydroxy-2-butynoate", (9, 10), "CC(C)(C)OC(=O)C#CCO"),
    ("mol_mscl", "MsCl", "methanesulfonyl chloride", (9, 10), "CS(=O)(=O)Cl"),
    ("mol_et3n", "Et3N", "triethylamine", (8, 9, 10), "CCN(CC)CC"),
    ("mol_hcl", "HCl", "hydrochloric acid", (10,), "Cl"),
    ("mol_libr", "LiBr", "lithium bromide", (10,), "[Li+].[Br-]"),
    ("mol_acetone", "acetone", "acetone", (10,), "CC(C)=O"),
    ("mol_s7", "S7", "tert-butyl 4-bromo-2-butynoate", (10, 11), "CC(C)(C)OC(=O)C#CCBr"),
    ("mol_bu3snh", "Bu3SnH", "tri-n-butyltin hydride", (10,), "[H][Sn](CCCC)(CCCC)CCCC"),
    ("mol_ethyl_succinyl_chloride", "ethyl succinyl chloride", "ethyl succinyl chloride", (10,), "CCOC(=O)CCC(=O)Cl"),
    ("mol_nmp", "NMP", "N-methylpyrrolidinone", (10,), "CN1CCCC1=O"),
    ("mol_s8", "S8", "ethyl 4-oxobutanoate", (10, 11), "CCOC(=O)CCC=O"),
    ("mol_triphosgene", "triphosgene", "bis(trichloromethyl) carbonate", (10,), "O=C(OC(Cl)(Cl)Cl)OC(Cl)(Cl)Cl"),
    ("mol_pmosulfonamide", "p-methoxybenzenesulfonamide", "p-methoxybenzenesulfonamide", (10,), "COc1ccc(S(N)(=O)=O)cc1"),
    ("mol_iprnco", "i-PrNCO", "isopropyl isocyanate", (10,), "CC(C)N=C=O"),
    ("mol_odcb", "1,2-dichlorobenzene", "1,2-dichlorobenzene", (10,), "Clc1ccccc1Cl"),
    ("mol_s9", "S9", "p-methoxybenzenesulfonyl isocyanate", (10,), "COc1ccc(S(=O)(=O)N=C=O)cc1"),
    ("mol_pnosulfonamide", "p-nitrobenzenesulfonamide", "p-nitrobenzenesulfonamide", (10,), "NS(=O)(=O)c1ccc([N+](=O)[O-])cc1"),
    ("mol_nbunco", "n-BuNCO", "n-butyl isocyanate", (10,), "CCCCN=C=O"),
    ("mol_s10", "S10", "p-nitrobenzenesulfonyl isocyanate", (10,), "O=[N+]([O-])c1ccc(S(=O)(=O)N=C=O)cc1"),
    ("mol_sncl2", "SnCl2", "tin(II) chloride", (11,), "[Cl-].[Cl-].[Sn+2]"),
    ("mol_nai", "NaI", "sodium iodide", (11,), "[Na+].[I-]"),
    ("mol_dmpu", "DMPU", "N,N'-dimethylpropyleneurea", (11,), "CN1CCCN(C)C1=O"),
    ("mol_propionaldehyde", "propionaldehyde", "propionaldehyde", (11,), "CCC=O"),
    ("mol_butyraldehyde", "butyraldehyde", "butyraldehyde", (11,), "CCCC=O"),
    ("mol_cyclopropanecarboxaldehyde", "cyclopropanecarboxaldehyde", "cyclopropanecarboxaldehyde", (11,), "O=CC1CC1"),
    ("mol_cyclopentanecarboxaldehyde", "cyclopentanecarboxaldehyde", "cyclopentanecarboxaldehyde", (11,), "O=CC1CCCC1"),
    ("mol_cyclohexanecarboxaldehyde", "cyclohexanecarboxaldehyde", "cyclohexanecarboxaldehyde", (11,), "O=CC1CCCCC1"),
    ("mol_4pentynal", "4-pentynal", "4-pentynal", (10, 11), "C#CCCC=O"),
    ("mol_hydrocinnamaldehyde", "hydrocinnamaldehyde", "3-phenylpropanal", (11,), "O=CCCc1ccccc1"),
    ("mol_tbdmso_butanal", "4-(tert-butyldimethylsilyloxy)butanal", "4-(tert-butyldimethylsilyloxy)butanal", (10, 11), "CC(C)(C)[Si](C)(C)OCCCC=O"),
    ("mol_s4", "S4", "ethyl 4-bromo-2-butynoate", (9, 11), "CCOC(=O)C#CCBr"),
]
for mid, label, name, pages, smiles in simple:
    add_simple(mid, label, name, pages, smiles)

# Drawing-defined identities that are kept unresolved or only partially encoded.
add_molecule("mol_30", "30", "ethyl 4-cyano-5-oxo-5-phenyl-2-vinylidenepentanoate", ev(7),
             smiles="CCOC(=O)C(=C=C)CC(C#N)C(=O)c1ccccc1",
             chirality="unspecified", geometry="not_applicable", representation="partial",
             note="The carbon bearing cyano/benzoyl substituents is a potential stereocenter; no configuration is reported.")
add_molecule("mol_31", "31", "ethyl 5-cyano-2-methyl-6-phenyl-4H-pyran-3-carboxylate", ev(7),
             chirality="unresolved", geometry="unresolved", representation="not_encoded")
add_molecule("mol_32", "32", "ethyl 5-cyano-4-methyl-6-phenyl-2H-pyran-3-carboxylate", ev(7),
             chirality="unresolved", geometry="unresolved", representation="not_encoded")
add_molecule("mol_39", "39", "benzyl tosyl carbamate", ev(8),
             smiles="Cc1ccc(S(=O)(=O)NC(=O)OCc2ccccc2)cc1",
             chirality="achiral", geometry="not_applicable")
add_molecule("mol_40", "40", "bicyclic lactam 40", ev(8), chirality="partially_specified",
             geometry="not_applicable", representation="not_encoded",
             assignment="relative wedge/dash stereochemistry is shown in Scheme 10; absolute configuration is not assigned")
add_molecule("mol_trachelanthamidine", "(±)-trachelanthamidine", "racemic trachelanthamidine", ev(8),
             chirality="racemic", geometry="not_applicable", representation="not_encoded",
             mixture="racemate denoted by (±)", assignment="relative stereochemistry shown; absolute enantiomer not selected")
add_molecule("mol_s12", "S12", "free carboxylic acid S12 from 3m", ev(8), chirality="unspecified",
             geometry="not_applicable", representation="not_encoded",
             assignment="pyrroline stereocenter is drawn without an absolute assignment")
add_molecule("mol_41", "41", "allylic alcohol 41", ev(8), chirality="unspecified",
             geometry="not_applicable", representation="not_encoded",
             assignment="pyrroline stereocenter is drawn without an absolute assignment")
add_molecule("mol_42", "42", "bicyclic lactam 42", ev(8), chirality="unspecified",
             geometry="not_applicable", representation="not_encoded")
add_molecule("mol_supinidine", "(±)-supinidine", "racemic supinidine", ev(8, 9),
             chirality="racemic", geometry="specified", representation="not_encoded",
             mixture="racemate denoted by (±)", assignment="alkene geometry is fixed by the ring drawing; absolute enantiomer not selected")


# Page 7: two proposed alternatives for the cited Tong reaction and one real experiment.
for num, title, description in [
    (5, "Alternative gamma-umpolung/SN2 mechanism for Tong's (4+1) reaction", "proposed route through intermediates 20 and 25--27"),
    (6, "Alternative allenoate pathway for Tong's (4+1) reaction", "proposed SN2 formation of allenoate 28 followed by intramolecular gamma-umpolung cyclization"),
]:
    r = reaction(f"mechanism_scheme{num}", ev(7), "proposed_mechanism", ["proposed_mechanism"],
                 title, description, f"Page 7, Scheme {num}", labels=["C-C"],
                 unresolved=["Intermediates are proposed rather than directly observed; this is interpretation of cited Tong chemistry, not a current-paper experiment."])
    for mid, role in [("mol_4", "reactant"), ("mol_19", "reactant"), ("mol_pph3", "catalyst"), ("mol_cs2co3", "reagent"), ("mol_6", "product")]:
        participant(r, mid, role, ev(7), product_role="cited product" if role == "product" else None)

r = reaction("rxn_scheme7_9_19", ev(7), "specific_experiment", ["mechanistic_study"],
             "Reaction of allenylic acetate 9 with benzoyl acetonitrile 19",
             "competition between SN2 and SN2-prime pathways followed by pyran formation",
             "Page 7, Scheme 7", labels=["C-C"],
             unresolved=["The page reports isolated product yields but does not state analytical detection limits or stereochemical composition of 30--32."])
p9 = participant(r, "mol_9", "reactant", ev(7), equiv=("1.0 equiv", 1.0))
p19 = participant(r, "mol_benzoyl_acetonitrile", "reactant", ev(7), equiv=("1.2 equiv", 1.2))
pcs = participant(r, "mol_cs2co3", "reagent", ev(7), equiv=("1.3 equiv", 1.3))
pph = participant(r, "mol_benzene", "solvent", ev(7))
for mid in ("mol_30", "mol_31", "mol_32"):
    participant(r, mid, "product", ev(7), product_role="isolated product")
step(r, "stir allenylic acetate 9, benzoyl acetonitrile 19 and cesium carbonate in benzene", ev(7),
     temperature="room temperature", time=("20 h", 20, "h"), participants=[p9, p19, pcs, pph])
for mid, value in (("mol_30", 19), ("mol_31", 27), ("mol_32", 12)):
    outcome(r, "yield", [mid], f"{value}%", q(value, "%"), ev(7), yield_type="isolated yield",
            basis="isolated product", product_role="isolated product")
outcome(r, "other", ["mol_30", "mol_31"], "46% total SN2-pathway products", q(46, "%"), ev(7),
        basis="author-summed isolated yields of 30 and 31", scope="author-assigned SN2 pathway total")
outcome(r, "other", ["mol_32"], "12% SN2′-pathway product", q(12, "%"), ev(7),
        basis="isolated yield interpreted by authors", scope="author-assigned SN2-prime pathway")

r = reaction("mechanism_scheme8", ev(7, 8), "proposed_mechanism", ["proposed_mechanism"],
             "Proposed mechanism for the current phosphine-catalyzed annulative rearrangement",
             "phosphonium diene formation, beta-prime addition and intramolecular gamma-umpolung cyclization",
             "Pages 7--8, Scheme 8", labels=["C-N"],
             unresolved=["Attack at the beta-prime carbon is the authors' favored proposal rather than a directly observed elementary step; intermediates 33--38 are proposed."])
for mid, role in [("mol_generic_1", "reactant"), ("mol_pph3", "catalyst"), ("mol_generic_2", "unknown"), ("mol_8", "unknown"), ("mol_generic_3", "product")]:
    participant(r, mid, role, ev(7, 8), product_role="proposed product" if role == "product" else None)


# Page 8: two distinct no-reaction controls and the current/cited synthesis arrows.
for level, molpct, label in (("catalytic", 20, "20 mol%"), ("stoichiometric", 100, "100 mol%")):
    r = reaction(f"rxn_scheme9_39_pph3_{level}", ev(8), "specific_experiment", ["mechanistic_study", "control"],
                 f"Benzyl tosyl carbamate 39 with {level} PPh3", "control for direct carbamate displacement/carbonyl addition",
                 f"Page 8, Scheme 9, {label} PPh3 condition", unresolved=["No analytical detection limit is reported for the no-reaction conclusion."])
    p39 = participant(r, "mol_39", "reactant", ev(8))
    pcat = participant(r, "mol_pph3", "catalyst" if molpct == 20 else "reagent", ev(8), molpct=(label, molpct))
    psol = participant(r, "mol_mecn", "solvent", ev(8))
    step(r, "stir benzyl tosyl carbamate 39 with PPh3 in acetonitrile", ev(8), temperature="room temperature",
         time=("24 h", 24, "h"), participants=[p39, pcat, psol])
    outcome(r, "no_reaction", ["mol_39"], "no reaction", "no reaction", ev(8),
            basis="author observation; analytical method not reported", scope=f"{label} PPh3 control")
    r["classification"]["labels"] = F([], [], ev(8), "not_applicable")
    r["classification"]["coverage"] = F("no observed transformation", "not_applicable", ev(8), "inference",
                                           derivation="The source explicitly reports no reaction for this control.")
    r["classification"]["limitations"] = ["No observed transformation; no product bond class assigned."]

r = reaction("rxn_scheme10_3l_40", ev(8), "specific_experiment", ["application_synthesis"],
             "Conversion of pyrroline 3l to lactam 40", "detosylation, lactamization, 1,4-reduction and transesterification",
             "Page 8, Scheme 10, first arrow", labels=["C-N", "C-H"],
             unresolved=["The one-pot arrow combines several transformations; a complete atom map is not supplied."])
p3l = participant(r, "mol_3l", "reactant", ev(8))
pmg = participant(r, "mol_mg", "reagent", ev(8), equiv=("25 equiv", 25))
pmeoh = participant(r, "mol_meoh", "solvent", ev(8))
participant(r, "mol_40", "product", ev(8), product_role="desired product")
step(r, "sonicate pyrroline 3l with magnesium in methanol", ev(8), time=("4.5 h", 4.5, "h"), participants=[p3l, pmg, pmeoh])
outcome(r, "yield", ["mol_40"], "86%", q(86, "%"), ev(8), yield_type="isolated yield", basis="isolated product")

r = reaction("lit_scheme10_40_trachelanthamidine", ev(8), "literature_example", ["cited_literature"],
             "Literature conversion of lactam 40 to (±)-trachelanthamidine", "formal-synthesis completion reported in reference 35",
             "Page 8, Scheme 10, arrow labeled [35]", labels=["C-H"],
             unresolved=["The current page gives no reagents, conditions or yield for the cited literature conversion."])
participant(r, "mol_40", "reactant", ev(8)); participant(r, "mol_trachelanthamidine", "product", ev(8), product_role="cited literature product")

r = reaction("rxn_scheme11_3m_s12", ev(8), "specific_experiment", ["application_synthesis"],
             "TFA-mediated deprotection of pyrroline 3m", "tert-butyl ester deprotection to carboxylic acid S12",
             "Page 8, Scheme 11, first arrow", labels=["O-H"], unresolved=["Coproducts from tert-butyl deprotection are not specified."])
p3m = participant(r, "mol_3m", "reactant", ev(8)); ptfa = participant(r, "mol_tfa", "reagent", ev(8)); pdcm = participant(r, "mol_dcm", "solvent", ev(8))
participant(r, "mol_s12", "product", ev(8), product_role="desired product")
step(r, "stir 3m in 15% trifluoroacetic acid in dichloromethane", ev(8), temperature="room temperature", time=("3 h", 3, "h"), participants=[p3m, ptfa, pdcm])
outcome(r, "yield", ["mol_s12"], "95%", q(95, "%"), ev(8), yield_type="isolated yield", basis="isolated product")

r = reaction("rxn_scheme11_s12_41", ev(8), "specific_experiment", ["application_synthesis"],
             "Conversion of carboxylic acid S12 to allylic alcohol 41", "mixed-anhydride formation followed by sodium borohydride reduction",
             "Page 8, Scheme 11, second arrow", labels=["C-H", "O-H"],
             unresolved=["The reported 73% is explicitly over two operations from S12 to 41 and is not assigned solely to the NaBH4 step."])
ps12 = participant(r, "mol_s12", "reactant", ev(8)); pcl = participant(r, "mol_ethyl_chloroformate", "reagent", ev(8)); pet3n = participant(r, "mol_et3n", "reagent", ev(8)); pthf = participant(r, "mol_thf", "solvent", ev(8)); pbh = participant(r, "mol_nabh4", "reagent", ev(8)); pme = participant(r, "mol_meoh", "solvent", ev(8))
participant(r, "mol_41", "product", ev(8), product_role="desired product")
step(r, "form the mixed anhydride from S12 with ethyl chloroformate and triethylamine in THF", ev(8), temperature=0, time=("1 h", 1, "h"), participants=[ps12, pcl, pet3n, pthf])
step(r, "cool and reduce the mixed anhydride with sodium borohydride in methanol", ev(8), temperature=-78, time=("1.5 h", 1.5, "h"), participants=[pbh, pme])
outcome(r, "yield", ["mol_41"], "73% over two steps", q(73, "%"), ev(8), yield_type="isolated yield over two steps", basis="isolated endpoint 41 from starting acid S12", scope="two-step sequence S12 to 41")

r = reaction("rxn_scheme11_41_42", ev(8, 9), "specific_experiment", ["application_synthesis"],
             "Detosylation/lactamization of allylic alcohol 41", "magnesium-mediated detosylation with in situ lactamization",
             "Pages 8--9, Scheme 11 and continuation text", labels=["C-N"], unresolved=["The source does not report separate yields for detosylation and lactamization."])
p41 = participant(r, "mol_41", "reactant", ev(8, 9)); pmg = participant(r, "mol_mg", "reagent", ev(8, 9), equiv=("50 equiv", 50)); pmeoh = participant(r, "mol_meoh", "solvent", ev(8, 9))
participant(r, "mol_42", "product", ev(8, 9), product_role="desired product")
step(r, "sonicate allylic alcohol 41 with magnesium in methanol", ev(8, 9), time=("3.5 h", 3.5, "h"), participants=[p41, pmg, pmeoh])
outcome(r, "yield", ["mol_42"], "72%", q(72, "%"), ev(8, 9), yield_type="isolated yield", basis="isolated product")
outcome(r, "yield", ["mol_42"], "50% overall from 3m", q(50, "%"), ev(9), yield_type="overall yield over four steps",
        basis="author-reported overall yield from 3m to endpoint lactam 42", scope="four-step sequence beginning at 3m and ending at 42")

r = reaction("lit_scheme11_42_supinidine", ev(8, 9), "literature_example", ["cited_literature"],
             "Literature conversion of lactam 42 to (±)-supinidine", "formal-synthesis completion reported in reference 37",
             "Page 8, Scheme 11, arrow labeled [37]", labels=["C-H"], unresolved=["The current page gives no reagents, conditions or yield for the cited literature conversion."])
participant(r, "mol_42", "reactant", ev(8, 9)); participant(r, "mol_supinidine", "product", ev(8, 9), product_role="cited literature product")


# Pages 9--10: multistage preparations of S6 and S7.
r = reaction("rxn_s6_from_s1", ev(9), "specific_experiment", ["substrate_preparation"],
             "Preparation of tert-butyl 4-hydroxy-2-butynoate S6", "three-step protection, alkynylation/carboxylation and deprotection sequence",
             "Page 9, Experimental Section, S6", labels=["C-C"],
             unresolved=["The 82% yield is explicitly over three steps; crude intermediate S5 is not assigned a final-step yield."])
parts = {}
for mid, role, kwargs in [
    ("mol_dhp", "reactant", {"volume": ("14.9 mL", 14.9, "mL"), "amount": ("160 mmol", 160, "mmol")}),
    ("mol_propargyl_alcohol", "reactant", {"volume": ("8.65 mL", 8.65, "mL"), "amount": ("145 mmol", 145, "mmol")}),
    ("mol_ptsoh", "catalyst", {"mass": ("0.276 g", 0.276, "g"), "amount": ("1.45 mmol", 1.45, "mmol")}),
    ("mol_thf", "solvent", {"volume": ("230 mL", 230, "mL")}),
    ("mol_nbuli", "reagent", {"volume": ("99.7 mL", 99.7, "mL"), "amount": ("160 mmol", 160, "mmol")}),
    ("mol_boc2o", "reactant", {"mass": ("39.5 g", 39.5, "g"), "amount": ("181 mmol", 181, "mmol")}),
    ("mol_nh4cl", "workup_agent", {"volume": ("450 mL saturated aqueous NH4Cl", 450, "mL"), "stage": "workup"}),
    ("mol_etoac", "workup_agent", {"volume": ("3 × 150 mL", 450, "mL"), "stage": "workup"}),
    ("mol_nacl", "workup_agent", {"volume": ("100 mL brine", 100, "mL"), "stage": "workup"}),
    ("mol_na2so4", "workup_agent", {"stage": "workup"}),
    ("mol_etoh", "solvent", {"volume": ("570 mL", 570, "mL")}),
    ("mol_ppts", "reagent", {"mass": ("3.72 g", 3.72, "g"), "amount": ("14.5 mmol", 14.5, "mmol")}),
    ("mol_nahco3", "workup_agent", {"volume": ("500 mL saturated aqueous NaHCO3", 500, "mL"), "stage": "workup"}),
    ("mol_et2o", "workup_agent", {"volume": ("5 × 100 mL", 500, "mL"), "stage": "workup"}),
]:
    parts[mid] = participant(r, mid, role, ev(9), **kwargs)
participant(r, "mol_s6", "product", ev(9), mass=("18.7 g", 18.7, "g"), product_role="desired product")
step(r, "cool 3,4-dihydro-2H-pyran and propargyl alcohol; add p-toluenesulfonic acid", ev(9), temperature=0, participants=[parts["mol_dhp"], parts["mol_propargyl_alcohol"], parts["mol_ptsoh"]])
step(r, "stir protected-alcohol formation mixture at 0 °C", ev(9), temperature=0, time=("10 min", 10, "min"))
step(r, "warm to room temperature and stir", ev(9), temperature="room temperature", time=("2 h", 2, "h"))
step(r, "add dry THF and cool", ev(9), temperature=-78, participants=[parts["mol_thf"]])
step(r, "slowly add 1.6 M n-BuLi in hexanes and hold", ev(9), temperature=-78, time=("1.5 h", 1.5, "h"), participants=[parts["mol_nbuli"]])
step(r, "add Boc2O solution in dry THF in 20 mL portions over 20 min", ev(9), temperature=-78, time=("20 min addition", 20, "min"), participants=[parts["mol_boc2o"]])
step(r, "stir after Boc2O addition", ev(9), temperature=-78, time=("1 h", 1, "h"))
step(r, "move to ice-water bath and stir", ev(9), temperature=0, time=("30 min", 30, "min"))
step(r, "quench, extract, wash, dry, filter and concentrate to crude S5 (43.3 g)", ev(9), stage="workup", participants=[parts["mol_nh4cl"], parts["mol_etoac"], parts["mol_nacl"], parts["mol_na2so4"]])
step(r, "dissolve crude S5 in absolute ethanol, add PPTS and heat", ev(9), temperature=55, time=("2 h", 2, "h"), participants=[parts["mol_etoh"], parts["mol_ppts"]])
step(r, "quench, remove ethanol, extract, wash, dry and concentrate", ev(9), stage="workup", participants=[parts["mol_nahco3"], parts["mol_et2o"], parts["mol_nacl"], parts["mol_na2so4"]])
step(r, "purify by flash chromatography (hexanes to 5%, 10%, then 15% EtOAc/hexanes)", ev(9), stage="purification")
outcome(r, "yield", ["mol_s6"], "82% over the three steps", q(82, "%"), ev(9), yield_type="isolated yield over three steps", basis="isolated S6", scope="three-step sequence from propargyl alcohol S1 to S6")
outcome(r, "product_amount", ["mol_s6"], "18.7 g", q(18.7, "g"), ev(9), basis="isolated product mass")

r = reaction("rxn_s7_from_s6", ev(9, 10), "specific_experiment", ["substrate_preparation"],
             "Preparation of tert-butyl 4-bromo-2-butynoate S7", "mesylation followed by lithium bromide displacement",
             "Pages 9--10, Experimental Section, S7", labels=["C-Br"],
             unresolved=["The source reports LiBr as 23.0 g and 26.5 mmol; these values are preserved without correcting their apparent inconsistency.", "The 80% yield is explicitly over two steps."])
parts = {}
for mid, role, kwargs in [
    ("mol_s6", "reactant", {"mass": ("11.3 g", 11.3, "g"), "amount": ("72.6 mmol", 72.6, "mmol")}),
    ("mol_mscl", "reagent", {"volume": ("6.20 mL", 6.20, "mL"), "amount": ("79.8 mmol", 79.8, "mmol")}),
    ("mol_et3n", "reagent", {"volume": ("11.0 mL", 11.0, "mL"), "amount": ("76.2 mmol", 76.2, "mmol")}),
    ("mol_dcm", "solvent", {"volume": ("400 mL", 400, "mL")}),
    ("mol_hcl", "workup_agent", {"volume": ("3 × 250 mL of 2 M HCl", 750, "mL"), "stage": "workup"}),
    ("mol_nahco3", "workup_agent", {"volume": ("2 × 250 mL saturated aqueous NaHCO3", 500, "mL"), "stage": "workup"}),
    ("mol_nacl", "workup_agent", {"volume": ("250 mL brine", 250, "mL"), "stage": "workup"}),
    ("mol_na2so4", "workup_agent", {"stage": "workup"}),
    ("mol_libr", "reactant", {"mass": ("23.0 g", 23.0, "g"), "amount": ("26.5 mmol", 26.5, "mmol")}),
    ("mol_acetone", "solvent", {"volume": ("20 mL + 460 mL", 480, "mL")}),
]:
    parts[mid] = participant(r, mid, role, ev(9, 10), **kwargs)
participant(r, "mol_s7", "product", ev(10), mass=("12.6 g", 12.6, "g"), product_role="desired product")
step(r, "add methanesulfonyl chloride dropwise to S6 and triethylamine in dry dichloromethane", ev(9), temperature=0, participants=[parts["mol_s6"], parts["mol_mscl"], parts["mol_et3n"], parts["mol_dcm"]])
step(r, "stir mesylation mixture", ev(9), temperature=0, time=("1 h", 1, "h"))
step(r, "add water, separate, wash, dry, filter and concentrate to crude propargylic mesylate (16.3 g)", ev(9, 10), stage="workup", participants=[parts["mol_hcl"], parts["mol_nahco3"], parts["mol_nacl"], parts["mol_na2so4"]])
step(r, "add crude mesylate solution in acetone to anhydrous LiBr in acetone", ev(10), participants=[parts["mol_libr"], parts["mol_acetone"]])
step(r, "heat displacement mixture under reflux", ev(10), temperature="reflux", time=("1.75 h", 1.75, "h"))
step(r, "cool, filter through Celite, concentrate, partition, extract, dry and concentrate", ev(10), stage="workup")
step(r, "purify by flash chromatography (hexanes to 10% EtOAc/hexanes)", ev(10), stage="purification")
outcome(r, "yield", ["mol_s7"], "80% over two steps", q(80, "%"), ev(10), yield_type="isolated yield over two steps", basis="isolated S7", scope="two-step sequence from S6 to S7")
outcome(r, "product_amount", ["mol_s7"], "12.6 g", q(12.6, "g"), ev(10), basis="isolated product mass")


# Page 10: S8--S10.
r = reaction("rxn_s8", ev(10), "specific_experiment", ["substrate_preparation"],
             "Preparation of ethyl 4-oxobutanoate S8", "selective reduction of ethyl succinyl chloride to aldehyde",
             "Page 10, Experimental Section, S8", labels=["C-H"])
p1 = participant(r, "mol_ethyl_succinyl_chloride", "reactant", ev(10), volume=("1.14 mL", 1.14, "mL"), amount=("8.00 mmol", 8.00, "mmol"))
p2 = participant(r, "mol_bu3snh", "reagent", ev(10), volume=("2.30 mL", 2.30, "mL"), amount=("8.40 mmol", 8.40, "mmol"))
p3 = participant(r, "mol_nmp", "solvent", ev(10), volume=("8 mL", 8, "mL"))
participant(r, "mol_s8", "product", ev(10), mass=("0.602 g", 0.602, "g"), product_role="desired product")
step(r, "add tri-n-butyltin hydride dropwise over 5 min to ethyl succinyl chloride in NMP", ev(10), temperature=0, time=("5 min addition", 5, "min"), participants=[p1, p2, p3])
step(r, "stir at 0 °C", ev(10), temperature=0, time=("30 min", 30, "min"))
step(r, "warm to room temperature and stir", ev(10), temperature="room temperature", time=("45 min", 45, "min"))
step(r, "add water and EtOAc; extract, wash with saturated NH4Cl and brine, dry, filter and concentrate", ev(10), stage="workup")
step(r, "distill under reduced pressure (head 100-110 °C; bath 150 °C)", ev(10), temperature=("head 100-110 °C; bath 150 °C", {"head_minimum": 100, "head_maximum": 110, "bath": 150, "unit": "°C"}), pressure="reduced pressure", stage="purification")
outcome(r, "yield", ["mol_s8"], "69%", q(69, "%"), ev(10), yield_type="isolated yield", basis="distilled isolated product")

r = reaction("rxn_s9", ev(10), "specific_experiment", ["substrate_preparation"],
             "Preparation of p-methoxybenzenesulfonyl isocyanate S9", "conversion of sulfonamide to sulfonyl isocyanate with triphosgene",
             "Page 10, Experimental Section, S9", labels=["C-N"])
p1 = participant(r, "mol_triphosgene", "reagent", ev(10), mass=("5.88 g", 5.88, "g"), amount=("19.8 mmol", 19.8, "mmol"))
p2 = participant(r, "mol_pmosulfonamide", "reactant", ev(10), mass=("11.0 g", 11.0, "g"), amount=("58.8 mmol", 58.8, "mmol"))
p3 = participant(r, "mol_iprnco", "reagent", ev(10), volume=("0.550 mL", 0.550, "mL"), amount=("5.60 mmol", 5.60, "mmol"))
p4 = participant(r, "mol_odcb", "solvent", ev(10), volume=("30 mL + 30 mL", 60, "mL"))
participant(r, "mol_s9", "product", ev(10), mass=("7.09 g", 7.09, "g"), product_role="desired product")
step(r, "add triphosgene solution in dry 1,2-dichlorobenzene by syringe pump to sulfonamide and isopropyl isocyanate solution", ev(10), temperature=150, time=("4.5 h addition", 4.5, "h"), participants=[p1, p2, p3, p4])
step(r, "heat reaction mixture", ev(10), temperature=180, time=("2 h", 2, "h"))
step(r, "remove from oil bath, cool to room temperature and stir", ev(10), temperature="room temperature", time=("12 h", 12, "h"))
step(r, "distill off solvent and isopropyl isocyanate (3 mmHg; bath 155 °C; head 120 °C)", ev(10), pressure="3 mmHg", temperature=("bath 155 °C; head 120 °C", {"bath": 155, "head": 120, "unit": "°C"}), stage="purification")
step(r, "further distill product (3 mmHg; bath 200 °C; head 160-165 °C)", ev(10), pressure="3 mmHg", temperature=("bath 200 °C; head 160-165 °C", {"bath": 200, "head_minimum": 160, "head_maximum": 165, "unit": "°C"}), stage="purification")
outcome(r, "yield", ["mol_s9"], "57%", q(57, "%"), ev(10), yield_type="isolated yield", basis="distilled isolated product")

r = reaction("rxn_s10", ev(10), "specific_experiment", ["substrate_preparation"],
             "Preparation of p-nitrobenzenesulfonyl isocyanate S10", "conversion of sulfonamide to sulfonyl isocyanate with triphosgene",
             "Page 10, Experimental Section, S10", labels=["C-N"])
p1 = participant(r, "mol_pnosulfonamide", "reactant", ev(10), mass=("5.50 g", 5.50, "g"), amount=("27.2 mmol", 27.2, "mmol"))
p2 = participant(r, "mol_nbunco", "reagent", ev(10), volume=("1.21 mL", 1.21, "mL"), amount=("10.9 mmol", 10.9, "mmol"))
p3 = participant(r, "mol_triphosgene", "reagent", ev(10), mass=("3.30 g", 3.30, "g"), amount=("11.1 mmol", 11.1, "mmol"))
p4 = participant(r, "mol_odcb", "solvent", ev(10), volume=("25 mL + 20 mL", 45, "mL"))
participant(r, "mol_s10", "product", ev(10), mass=("0.316 g", 0.316, "g"), product_role="desired product")
step(r, "heat p-nitrobenzenesulfonamide and n-butyl isocyanate in dry 1,2-dichlorobenzene under reflux until clear-yellow", ev(10), temperature="reflux", participants=[p1, p2, p4])
step(r, "add triphosgene solution by syringe pump", ev(10), temperature="reflux", time=("3.5 h addition", 3.5, "h"), participants=[p3])
step(r, "continue heating under reflux", ev(10), temperature="reflux", time=("2 h", 2, "h"))
step(r, "cool to room temperature and stir", ev(10), temperature="room temperature", time=("12 h", 12, "h"))
step(r, "distill off solvent and n-butyl isocyanate, then Kugelrohr-distill product", ev(10), pressure="reduced pressure", stage="purification")
outcome(r, "yield", ["mol_s10"], "5%", q(5, "%"), ev(10), yield_type="isolated yield", basis="Kugelrohr-distilled isolated product")


# Replace the page-3 summaries for 7a/7b with full page-10/11 procedures.
def reset_record(r, evidence, location, unresolved):
    r["evidence_ids"] = list(dict.fromkeys(r["evidence_ids"] + evidence))
    r["location"] = F(location, location, r["evidence_ids"])
    r["participants"] = []
    r["steps"] = []
    r["outcomes"] = []
    r["procedure"]["original_text"] = missing(r["evidence_ids"], "Procedure is represented as ordered structured steps from the source text.")
    r["review"]["status"] = "needs_review"
    r["review"]["method"] = "assistant_visual"
    r["review"]["checked_fields"] = ["original page image", "compound labels", "all stated quantities", "ordered conditions", "workup", "purification", "yield basis"]
    r["review"]["unresolved"] = unresolved


for suffix in ("a", "b"):
    rid = f"{DOC}:rxn_mbh_7{suffix}"
    page = 10 if suffix == "a" else 11
    r = reactions[rid]
    reset_record(r, ev(page), f"Page {page}, Experimental Section, preparation of 7{suffix}; summary also on page 3",
                 ["No external identity verification was performed; machine record remains needs_review."])
    if suffix == "a":
        pform = participant(r, "mol_paraformaldehyde", "reactant", ev(10), mass=("2.82 g", 2.82, "g"), amount=("89.2 mmol", 89.2, "mmol"))
        pcat = participant(r, "mol_hydroxyquinuclidine", "catalyst", ev(10), mass=("472 mg", 472, "mg"), amount=("3.57 mmol", 3.57, "mmol"))
        pall = participant(r, "mol_ethyl_allenoate", "reactant", ev(10), mass=("2.01 g", 2.01, "g"), amount=("17.90 mmol", 17.90, "mmol"))
        pthf = participant(r, "mol_thf", "solvent", ev(10), volume=("50 + 25 + 25 mL", 100, "mL"))
        participant(r, "mol_7a", "product", ev(10), mass=("1.67 g", 1.67, "g"), product_role="desired product")
        step(r, "hold paraformaldehyde under vacuum to remove water", ev(10), temperature=55, time=("1 h", 1, "h"), participants=[pform])
        step(r, "cool paraformaldehyde suspension in THF", ev(10), temperature=-15, participants=[pform, pthf])
        step(r, "add 3-hydroxyquinuclidine solution dropwise", ev(10), temperature=-15, time=("5 min addition", 5, "min"), participants=[pcat])
        step(r, "add ethyl 2,3-butadienoate solution dropwise", ev(10), temperature=-15, time=("5 min addition", 5, "min"), participants=[pall])
        step(r, "stir at -15 °C", ev(10), temperature=-15, time=("1 h", 1, "h"))
        step(r, "remove cooling bath, warm to room temperature and stir", ev(10), temperature="room temperature", time=("3.75 h", 3.75, "h"))
        step(r, "quench with saturated NH4Cl; extract with EtOAc; wash with brine; dry, filter and concentrate", ev(10), stage="workup")
        step(r, "purify by flash chromatography (15% EtOAc/hexanes)", ev(10), stage="purification")
        outcome(r, "yield", ["mol_7a"], "66% yield", q(66, "%"), ev(10), yield_type="isolated yield", basis="isolated chromatographed product")
        outcome(r, "product_amount", ["mol_7a"], "1.67 g", q(1.67, "g"), ev(10), basis="isolated product mass")
    else:
        pcat = participant(r, "mol_hydroxyquinuclidine", "catalyst", ev(10, 11), mass=("460 mg", 460, "mg"), amount=("3.55 mmol", 3.55, "mmol"))
        pall = participant(r, "mol_ethyl_allenoate", "reactant", ev(10, 11), mass=("2.03 g", 2.03, "g"), amount=("18.1 mmol", 18.1, "mmol"))
        pace = participant(r, "mol_acetaldehyde", "reactant", ev(10, 11), volume=("2.00 mL", 2.00, "mL"), amount=("35.8 mmol", 35.8, "mmol"))
        pether = participant(r, "mol_et2o", "solvent", ev(10, 11), volume=("10 mL", 10, "mL"))
        participant(r, "mol_7b", "product", ev(11), mass=("2.20 g", 2.20, "g"), product_role="desired product")
        step(r, "add ethyl 2,3-butadienoate dropwise to 3-hydroxyquinuclidine in dry ether", ev(10, 11), temperature=-15, time=("5 min addition", 5, "min"), participants=[pall, pcat, pether])
        step(r, "add pre-cooled acetaldehyde dropwise", ev(10, 11), temperature=-15, time=("1 min addition", 1, "min"), atmosphere="argon", participants=[pace])
        step(r, "stir at -15 °C", ev(10, 11), temperature=-15, time=("50 min", 50, "min"))
        step(r, "quench with saturated NH4Cl; extract with EtOAc; wash with water and brine; dry, filter and concentrate", ev(11), stage="workup")
        step(r, "purify by flash chromatography (10% EtOAc/hexanes)", ev(11), stage="purification")
        outcome(r, "yield", ["mol_7b"], "79% yield", q(79, "%"), ev(11), yield_type="isolated yield", basis="isolated chromatographed product")
        outcome(r, "product_amount", ["mol_7b"], "2.20 g", q(2.20, "g"), ev(11), basis="isolated product mass")


# Page 11: the explicit 7c general procedure and every stated analogue are
# separate experiments because the reactant/product and yield differ.
allenylation = {
    "c": ("mol_propionaldehyde", 52, ("373 mg", 373, "mg")),
    "d": ("mol_butyraldehyde", 55, None),
    "e": ("mol_cyclopropanecarboxaldehyde", 65, None),
    "f": ("mol_cyclopentanecarboxaldehyde", 62, None),
    "g": ("mol_cyclohexanecarboxaldehyde", 80, None),
    "i": ("mol_4pentynal", 29, None),
    "j": ("mol_hydrocinnamaldehyde", 67, None),
    "k": ("mol_tbdmso_butanal", 44, None),
    "l": ("mol_s8", 61, None),
    "m": ("mol_s8", 72, None),
}
for suffix, (aldehyde, yld, product_mass) in allenylation.items():
    r = reaction(f"rxn_allenylation_7{suffix}", ev(11), "specific_experiment", ["substrate_preparation", "substrate_scope"],
                 f"Preparation of allenylic alcohol 7{suffix}", "tin(II)-mediated allenylation of propargyl bromide to aldehyde",
                 f"Page 11, Experimental Section, 7{suffix}", labels=["C-C"],
                 unresolved=["The exact aldehyde quantity/scale is stated only for 7c; analogues inherit the reported general procedure but their individual starting quantities are unreported on this page."] if suffix != "c" else ["No external identity verification was performed; machine record remains needs_review."])
    bromide = "mol_s7" if suffix == "m" else "mol_s4"
    pb = participant(r, bromide, "reactant", ev(11), **({"mass": ("1.00 g", 1.00, "g"), "amount": ("5.24 mmol", 5.24, "mmol")} if suffix == "c" else {}))
    psn = participant(r, "mol_sncl2", "reagent", ev(11), **({"mass": ("1.19 g", 1.19, "g"), "amount": ("6.28 mmol", 6.28, "mmol")} if suffix == "c" else {}))
    pnai = participant(r, "mol_nai", "reagent", ev(11), **({"mass": ("0.942 g", 0.942, "g"), "amount": ("6.28 mmol", 6.28, "mmol")} if suffix == "c" else {}))
    pdmpu = participant(r, "mol_dmpu", "solvent", ev(11), **({"volume": ("8 mL + 4 mL", 12, "mL")} if suffix == "c" else {}))
    pa = participant(r, aldehyde, "reactant", ev(11), **({"volume": ("0.300 mL", 0.300, "mL"), "amount": ("4.19 mmol", 4.19, "mmol")} if suffix == "c" else {}))
    participant(r, f"mol_7{suffix}", "product", ev(11), **({"mass": product_mass} if product_mass else {}), product_role="desired product")
    step(r, "add SnCl2 and NaI sequentially in one portion to propargyl bromide solution in DMPU", ev(11), temperature="room temperature", participants=[pb, psn, pnai, pdmpu])
    step(r, "stir protected from light", ev(11), temperature="room temperature", time=("7-8 h", {"minimum": 7, "maximum": 8, "unit": "h"}, "object"))
    # Normalize the range without losing its printed form.
    r["steps"][-1]["time"] = F("7-8 h", {"minimum": 7, "maximum": 8, "unit": "h"}, ev(11))
    step(r, "cool reaction mixture", ev(11), temperature=0)
    step(r, "add aldehyde solution over 15 min", ev(11), temperature=0, time=("15 min addition", 15, "min"), participants=[pa])
    step(r, "stir after aldehyde addition", ev(11), temperature=0, time=("2 h", 2, "h"))
    step(r, "warm to room temperature and stir protected from light", ev(11), temperature="room temperature")
    r["steps"][-1]["time"] = F("18-24 h", {"minimum": 18, "maximum": 24, "unit": "h"}, ev(11))
    step(r, "cool to 0 °C, dilute with ether, quench with saturated NH4Cl, extract, wash with brine, dry, filter and concentrate", ev(11), temperature=0, stage="workup")
    step(r, "purify by flash chromatography (5-10% EtOAc/hexanes)", ev(11), stage="purification")
    outcome(r, "yield", [f"mol_7{suffix}"], f"{yld}% yield", q(yld, "%"), ev(11), yield_type="isolated yield", basis="isolated chromatographed product")
    if product_mass:
        outcome(r, "product_amount", [f"mol_7{suffix}"], product_mass[0], q(product_mass[1], product_mass[2]), ev(11), basis="isolated product mass")


# Page 12: enrich existing carbamate records 1a--1k; do not duplicate them.
for suffix in list("abcdefghijk"):
    rid = f"{DOC}:rxn_carb_1{suffix}"
    r = reactions[rid]
    reset_record(r, ev(12), f"Page 12, Experimental Section, preparation of 1{suffix}; scope summary also on page 3",
                 ["For 1c--1k the page states that the procedure for 1b was followed; analogue-specific reagent amounts are not reported on page 12."] if suffix not in "ab" else ["No external identity verification was performed; machine record remains needs_review."])
    alcohol = f"mol_7{suffix}"
    if suffix == "a":
        pa = participant(r, alcohol, "reactant", ev(12), mass=("0.500 g", 0.500, "g"), amount=("3.52 mmol", 3.52, "mmol"))
        pi = participant(r, "mol_tsnco", "reactant", ev(12), volume=("0.590 mL", 0.590, "mL"), amount=("3.71 mmol", 3.71, "mmol"))
        ps = participant(r, "mol_dcm", "solvent", ev(12), volume=("8.8 mL", 8.8, "mL"))
        participant(r, "mol_1a", "product", ev(12), mass=("1.10 g", 1.10, "g"), product_role="desired product")
        step(r, "add p-toluenesulfonyl isocyanate dropwise quickly to allenylic alcohol 7a in dry dichloromethane", ev(12), temperature="room temperature", participants=[pa, pi, ps])
        step(r, "stir at room temperature", ev(12), temperature="room temperature", time=("1.5 h", 1.5, "h"))
        step(r, "concentrate under reduced pressure", ev(12), stage="workup")
        step(r, "purify by flash chromatography (10-30% EtOAc/hexanes); optional recrystallization from EtOAc/hexanes is also reported", ev(12), stage="purification")
        outcome(r, "yield", ["mol_1a"], "92% yield", q(92, "%"), ev(12), yield_type="isolated yield", basis="isolated chromatographed product")
        outcome(r, "product_amount", ["mol_1a"], "1.10 g", q(1.10, "g"), ev(12), basis="isolated product mass")
    else:
        kwargs_alcohol = {"mass": ("0.400 g", 0.400, "g"), "amount": ("2.56 mmol", 2.56, "mmol")} if suffix == "b" else {}
        kwargs_iso = {"volume": ("0.430 mL", 0.430, "mL"), "amount": ("2.70 mmol", 2.70, "mmol")} if suffix == "b" else {}
        kwargs_sol = {"volume": ("6.4 mL", 6.4, "mL")} if suffix == "b" else {}
        pa = participant(r, alcohol, "reactant", ev(12), **kwargs_alcohol)
        pi = participant(r, "mol_tsnco", "reactant", ev(12), **kwargs_iso)
        ps = participant(r, "mol_dcm", "solvent", ev(12), **kwargs_sol)
        prod_kwargs = {"mass": ("0.846 g", 0.846, "g")} if suffix == "b" else {}
        participant(r, f"mol_1{suffix}", "product", ev(12), **prod_kwargs, product_role="desired product")
        step(r, "add p-toluenesulfonyl isocyanate dropwise to allenylic alcohol in dry dichloromethane", ev(12), temperature=0, participants=[pa, pi, ps])
        step(r, "stir at 0 °C", ev(12), temperature=0, time=("1 h", 1, "h"))
        step(r, "remove ice bath, warm to room temperature and stir", ev(12), temperature="room temperature", time=("30 min", 30, "min"))
        step(r, "concentrate under reduced pressure", ev(12), stage="workup")
        step(r, "purify by flash chromatography (5-20% EtOAc/hexanes)", ev(12), stage="purification")
        yld = {"b": 94, "c": 84, "d": 69, "e": 88, "f": 86, "g": 88, "h": 88, "i": 89, "j": 72, "k": 72}[suffix]
        outcome(r, "yield", [f"mol_1{suffix}"], f"{yld}% yield", q(yld, "%"), ev(12), yield_type="isolated yield", basis="isolated chromatographed product")
        if suffix == "b":
            outcome(r, "product_amount", ["mol_1b"], "0.846 g", q(0.846, "g"), ev(12), basis="isolated product mass")


# Source context: characterization and literature-only preparation statements
# are retained without becoming synthetic experiment records.
contexts = [
    {"context_id": "pages7_8_mechanistic_interpretation", "document_id": DOC, "evidence_ids": ev(7, 8), "kind": "author_interpretation",
     "note": "The prose compares Schemes 4--8 and favors an SN2/allenoate pathway for Tong's cited chemistry while proposing beta-prime addition for the current annulation; these interpretations are not promoted to observed elementary steps."},
    {"context_id": "pages9_12_characterization_inventory", "document_id": DOC, "evidence_ids": ev(9, 10, 11, 12), "kind": "characterization_inventory",
     "note": "IR, NMR and MS characterization paragraphs identify isolated compounds but are not independent synthetic experiments; no standalone mass-spectrum peak/fragment table appears on pages 7--12."},
    {"context_id": "page9_general_information", "document_id": DOC, "evidence_ids": ev(9), "kind": "experimental_context",
     "note": "General dry-solvent, argon, chromatography and analytical-instrument information is retained as document context and not copied as an unreported condition into every reaction."},
    {"context_id": "page9_starting_material_literature", "document_id": DOC, "evidence_ids": ev(9, 10), "kind": "cited_literature",
     "note": "S4 and 4-(tert-butyldimethylsilyloxy)butanal/4-pentynal are stated to follow literature procedures; no current-paper preparation was invented from those citations."},
    {"context_id": "page10_7h_literature_preparation", "document_id": DOC, "evidence_ids": ev(10), "kind": "cited_literature",
     "note": "Allenylic alcohol 7h is stated to have been prepared according to a literature procedure; the present page gives no current experimental details, so no new current-paper experiment record was created."},
]
existing_contexts = {c["context_id"] for c in bundle.get("source_context", [])}
bundle.setdefault("source_context", []).extend(c for c in contexts if c["context_id"] not in existing_contexts)

# Machine records must remain needs_review.
for r in bundle["reactions"]:
    r["review"]["status"] = "needs_review"

bundle_path.write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")

# Raw inventory retained separately from formal records.
inventory = {
    "document_id": DOC,
    "source_name": "ChemCatChem 2020, 12, 4352.pdf",
    "inspected_range": [7, 12],
    "inspection_mode": "original-detail page image checked against extracted page text",
    "pages": [
        {"page": 7, "visual_evidence_id": I[7], "text_evidence_id": T[7], "visually_inspected": True,
         "inventory": ["Scheme 5 proposed mechanism", "Scheme 6 proposed mechanism", "Scheme 7 experiment", "start of Scheme 8 discussion"]},
        {"page": 8, "visual_evidence_id": I[8], "text_evidence_id": T[8], "visually_inspected": True,
         "inventory": ["Scheme 8 proposed mechanism", "Scheme 9 two no-reaction controls", "Scheme 10 current and cited arrows", "Scheme 11 current and cited arrows"]},
        {"page": 9, "visual_evidence_id": I[9], "text_evidence_id": T[9], "visually_inspected": True,
         "inventory": ["Scheme 11 continuation and 50% overall yield", "conclusions", "general information", "starting-material context", "S6 procedure", "start of S7"]},
        {"page": 10, "visual_evidence_id": I[10], "text_evidence_id": T[10], "visually_inspected": True,
         "inventory": ["S7 continuation", "S8", "S9", "S10", "general alcohol-series statement", "full 7a", "start of 7b"]},
        {"page": 11, "visual_evidence_id": I[11], "text_evidence_id": T[11], "visually_inspected": True,
         "inventory": ["7b continuation", "general 7c procedure", "7c--7g", "7i--7m", "spectral characterization"]},
        {"page": 12, "visual_evidence_id": I[12], "text_evidence_id": T[12], "visually_inspected": True,
         "inventory": ["general carbamate statement", "full 1a", "full 1b", "1c--1k yields/characterization"]},
    ],
    "outside_range_claimed": False,
}
(ROOT / "runs/current/raw_model/page_inventory_0007_0012.json").write_text(json.dumps(inventory, ensure_ascii=False, indent=2) + "\n")

new_ids = [r["reaction_id"] for r in bundle["reactions"] if any(e in r["evidence_ids"] for p in range(7, 13) for e in ev(p))]
(ROOT / "runs/current/raw_model/record_ids_pages_0007_0012.json").write_text(json.dumps(new_ids, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({
    "molecules": len(bundle["molecules"]),
    "reactions": len(bundle["reactions"]),
    "source_context": len(bundle.get("source_context", [])),
    "records_with_pages_7_12_evidence": len(new_ids),
}, ensure_ascii=False))
