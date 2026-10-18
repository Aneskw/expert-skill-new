#!/usr/bin/env python3
"""Build the incremental v3.1 extraction for visually inspected PDF pages 7-12.

Existing pages 1-6 reactions are only copied when the detailed experimental
paragraph on pages 7-12 is the same experiment and needs an added evidence
link/measurement basis.  New reaction records are limited to experiments first
reported in the bounded pages.
"""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = ROOT / "skills" / "chemical-reaction-extraction" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from standard_records import (  # noqa: E402
    bond_label,
    fact,
    new_molecule,
    new_outcome,
    new_participant,
    new_reaction,
    new_step,
    read,
    stereo_inspect,
    write,
)

DOC = "doc_66a75af3fbad9e13"
OUT = ROOT / "runs" / "current"
BASE = read(OUT / "standard_records.json")
SUBMITTED = read(OUT / "submitted_pages_0007_0012.json")
P = {i: f"{DOC}:page{i}:image" for i in range(7, 13)}


def F(raw=None, norm=None, ev=(), method="original_extraction", status=None,
      verification=None, derivation=None, alternatives=()):
    return fact(raw, norm, ev, method, status, verification, derivation, alternatives)


def NA(ev=()):
    return fact(evidence=ev, method="not_applicable", status="not_applicable", verification="missing")


def q(number, unit, raw=None, qualifier=None):
    v = {"number": number, "unit": unit, "raw_text": raw or f"{number} {unit}"}
    if qualifier:
        v["qualifier"] = qualifier
    return v


def qr(minimum, maximum, unit, raw):
    return {"minimum": minimum, "maximum": maximum, "unit": unit, "raw_text": raw}


old_molecules = {m["molecule_id"]: m for m in BASE["molecules"]}
old_reactions = {r["reaction_id"]: r for r in BASE["reactions"]}
new_molecules = {}
updated_reactions = []
new_reactions = []
stereo_checks = {}


def mid(label):
    return f"{DOC}:mol_{label}"


def ensure_old(label):
    molecule_id = mid(label)
    if molecule_id not in old_molecules:
        raise KeyError(f"missing existing molecule {molecule_id}")
    return molecule_id


def add_molecule(label, name, ev, smiles=None, aliases=None, allene=False,
                 unresolved_reason=None):
    molecule_id = mid(label)
    m = new_molecule(molecule_id, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, ev)
    m["names"]["original"] = F(name, name, ev)
    m["names"]["english"] = F(name, name, ev)
    m["names"]["abbreviations"] = F(aliases, aliases, ev) if aliases else fact(evidence=ev)
    m["names"]["aliases"] = fact(evidence=ev)
    m["cas_numbers"] = fact(evidence=ev)
    m["external_properties"] = []
    m["generic_template"] = None
    if smiles:
        diag = stereo_inspect(smiles)
        stereo_checks[molecule_id] = diag
        m["representation_type"] = "smiles"
        derivation = (
            "Connectivity candidate was transcribed from the source drawing or explicit "
            "compound heading and normalized with the bundled stereo inspector; parsing "
            "does not establish chemical correctness."
        )
        m["smiles"] = F(name, diag["canonical_isomeric_smiles"], ev,
                        "structure_derivation", derivation=derivation)
        m["stereochemistry"]["isomeric_smiles"] = F(
            name, diag["canonical_isomeric_smiles"], ev, "structure_derivation",
            derivation="Canonical isomeric SMILES returned for the machine-transcribed candidate graph.")
        m["stereochemistry"]["identity_key"] = F(
            name, diag["stereo_graph_key"], ev, "structure_derivation",
            derivation="Stereo graph key returned by the bundled inspector; not a sample identity.")
        m["stereochemistry"]["atom_centers"] = F(
            name, diag["atom_centers"], ev, "structure_derivation",
            derivation="Potential/assigned tetrahedral centers from the candidate graph.")
        m["stereochemistry"]["double_bonds"] = F(
            name, diag["double_bonds"], ev, "structure_derivation",
            derivation="Double-bond diagnostics from the candidate graph; no source-unreported E/Z is inferred.")
        unspecified_atoms = [x for x in diag["potential_stereo"]
                             if x["specified"] == "Unspecified" and str(x["type"]).startswith("Atom_")]
        unspecified_bonds = [x for x in diag["potential_stereo"]
                             if x["specified"] == "Unspecified" and str(x["type"]).startswith("Bond_")]
        chirality = "unspecified" if allene or unspecified_atoms else "achiral"
        geometry = "unspecified" if unspecified_bonds else "not_applicable"
        representation = "partial" if allene or unspecified_atoms or unspecified_bonds else "encoded"
        m["stereochemistry"]["axial"] = (
            F("Allene axis is present; no R_a/S_a or P/M descriptor is stated", {
                "axis": "allene", "assignment": "unspecified"
            }, ev, "structure_derivation",
              derivation="The allene functionality is source-supported, but the toolkit does not certify its axial configuration.")
            if allene else NA(ev)
        )
        m["stereochemistry"]["mixture"] = fact(evidence=ev)
        m["stereochemistry"]["assignment"] = F(
            name, {"tool_diagnostic": diag["potential_stereo"]}, ev, "structure_derivation",
            derivation="Candidate-graph diagnostic only; absent source descriptors remain unspecified.")
        m["stereochemistry"]["chirality_status"] = F(
            name, chirality, ev, "structure_derivation",
            derivation="Status follows the inspected drawing/name and candidate-graph diagnostics; unspecified is not treated as racemic or achiral.")
        m["stereochemistry"]["geometry_status"] = F(
            name, geometry, ev, "structure_derivation",
            derivation="No source E/Z or ring cis/trans descriptor is invented.")
        m["stereochemistry"]["representation_status"] = F(
            name, representation, ev, "structure_derivation",
            derivation="Partial is used when an allene axis or potential stereo element is not encoded.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = fact(evidence=ev)
        m["stereochemistry"]["atom_centers"] = fact(evidence=ev)
        m["stereochemistry"]["double_bonds"] = fact(evidence=ev)
        m["stereochemistry"]["axial"] = NA(ev)
        m["stereochemistry"]["mixture"] = fact(evidence=ev)
        reason = unresolved_reason or "The bounded source does not establish a unique machine-readable structure."
        m["stereochemistry"]["assignment"] = F(
            reason, None, ev, "inference", "unreadable", "inferred",
            "Structure/stereochemistry remains unresolved.")
        for key in ("chirality_status", "geometry_status", "representation_status"):
            status = "not_encoded" if key == "representation_status" else "unresolved"
            m["stereochemistry"][key] = F(name, status, ev, "inference", derivation=reason)
        m["stereochemistry"]["isomeric_smiles"] = fact(evidence=ev)
        m["stereochemistry"]["identity_key"] = fact(evidence=ev)
    note = unresolved_reason or (
        "Machine registry entry derived from the cited source region; structure and stereochemical "
        "identity remain subject to expert review."
    )
    m["identity_note"] = F(note, note, ev, "inference", derivation="Machine extraction decision.")
    new_molecules[molecule_id] = m
    return molecule_id


def add_named(label, name, ev, smiles=None, aliases=None):
    return add_molecule(label, name, ev, smiles=smiles, aliases=aliases,
                        unresolved_reason=None if smiles else
                        "Named reagent is source-supported, but no unique structural representation was asserted in this extraction.")


def append_unique(seq, value):
    if value not in seq:
        seq.append(value)


def update_fact_evidence(f, ev):
    for eid in ev:
        append_unique(f["evidence_ids"], eid)


def annotate_existing(short_id, pages, heading, basis=True, solvent_fix=None):
    rid = f"{DOC}:{short_id}"
    r = copy.deepcopy(old_reactions[rid])
    ev = [P[p] for p in pages]
    for eid in ev:
        append_unique(r["evidence_ids"], eid)
    r["procedure"]["original_text"] = F(
        heading, heading, ev, "original_extraction")
    r["notes"].append(F(
        f"The detailed experimental paragraph for this already-recorded scheme/table experiment occurs on PDF page(s) {', '.join(map(str, pages))}; it is merged as corroborating evidence rather than counted as a new experiment.",
        f"Detailed-procedure duplicate merged from pages {pages}", ev, "inference",
        derivation="Deduplication by matching compound labels, substrate/product identities, conditions, and yield context."))
    if basis:
        for o in r["outcomes"]:
            if o["kind"] == "yield" and o["measurement_basis"]["normalized_value"] is None:
                o["measurement_basis"] = F(
                    "afforded after the stated chromatographic or crystallization purification",
                    "isolated", ev)
    if solvent_fix:
        participant_id, molecule_label, raw = solvent_fix
        p = next(p for p in r["participants"] if p["participant_id"] == participant_id)
        p["molecule_id"] = mid(molecule_label)
        r["solvent_conditions"]["composition"] = F(raw, [raw], ev)
        r["review"]["unresolved"].append(
            f"Earlier scheme extraction listed a different solvent; the detailed paragraph on page {pages[0]} explicitly states {raw}, which is retained here.")
    updated_reactions.append(r)
    return r


def charge(r, participant_id, ev, *, mass=None, amount=None, volume=None, equivalents=None, mol_percent=None):
    p = next(p for p in r["participants"] if p["participant_id"] == participant_id)
    for key, value in (("mass", mass), ("amount", amount), ("volume", volume),
                       ("equivalents", equivalents), ("mol_percent", mol_percent)):
        if value is not None:
            raw, norm = value
            p["charges"][key] = F(raw, norm, ev)
    return p


def add_detail_outcome(r, outcome_id, target, raw, norm, ev, scope):
    if any(o["outcome_id"] == outcome_id for o in r["outcomes"]):
        return
    o = new_outcome(outcome_id, "yield", [mid(target)], ev)
    o["product_role"] = F("main product", "main product", ev)
    o["value"] = F(raw, norm, ev)
    o["yield_type"] = F("isolated yield", "isolated", ev)
    o["measurement_basis"] = F("isolated after column chromatography", "isolated", ev)
    o["scope"] = F(scope, scope, ev)
    r["outcomes"].append(o)


# Same-experiment detailed paragraphs: evidence and product-specific measurement
# basis are merged. Selected explicit charges are added where they materially
# improve the earlier scheme/table records.
r = annotate_existing("scheme3_1_to_2", [7], "1-(Cyclohex-1-en-1-yl)-5-phenylpent-1-yn-3-ol (2): detailed experimental procedure")
charge(r, "s3_1", [P[7]], mass=("0.584 g", q(0.584, "g")), amount=("5.50 mmol", q(5.50, "mmol")), equivalents=("1.2 equiv", q(1.2, "equiv")))
charge(r, "s3_nbuli", [P[7]], volume=("1.81 mL", q(1.81, "mL")), amount=("4.58 mmol", q(4.58, "mmol")), equivalents=("1.0 equiv", q(1.0, "equiv")))
charge(r, "s3_thf", [P[7]], volume=("40 mL", q(40, "mL")))
charge(r, "s3_aldehyde", [P[7]], mass=("0.615 g", q(0.615, "g")), amount=("4.58 mmol", q(4.58, "mmol")), equivalents=("1.0 equiv", q(1.0, "equiv")))
charge(r, "s3_2", [P[7]], mass=("1.073 g", q(1.073, "g")))
for participant_id in ("s3_1", "s3_nbuli", "s3_aldehyde"):
    participant = next(x for x in r["participants"] if x["participant_id"] == participant_id)
    participant["equivalence_reference"] = F(
        "3-phenylpropionaldehyde (1.0 equiv)", mid("phenylpropanal"), [P[7]])
add_detail_outcome(r, "yield_2_detailed_97", "2", "97%", q(97, "%"), [P[7]], "detailed page-7 preparation; the earlier Scheme 3 range 94-99% is retained separately")

r = annotate_existing("scheme3_2_to_3", [7, 8], "5-(Cyclohex-1-en-1-yl)hexa-3,4-dien-1-ylbenzene (3): detailed experimental procedure")
charge(r, "s3b_thf", [P[7], P[8]], volume=("0.5 mL in mesylate preparation; 3 mL in cuprate preparation; 0.25 mL rinse", None))
charge(r, "s3b_mscl", [P[7]], volume=("27 uL", q(27, "uL")), mass=("0.039 g", q(0.039, "g")), amount=("0.35 mmol", q(0.35, "mmol")))
charge(r, "s3b_cubr", [P[7]], mass=("0.119 g", q(0.119, "g")), amount=("0.790 mmol", q(0.790, "mmol")))
charge(r, "s3b_libr", [P[7]], mass=("0.068 g", q(0.068, "g")), amount=("0.790 mmol", q(0.790, "mmol")))
charge(r, "s3b_3", [P[8]], mass=("0.058 g", q(0.058, "g")))

for short_id, pages, heading in [
    ("table2_entry_2", [8], "3-Ethyl-4-methyl-2-tosyl-5,6,7,8-tetrahydroquinoline (7): detailed procedure"),
    ("table2_entry_1", [9], "4-Methyl-3-(3-phenylpropyl)-2-tosyl-5,6,7,8-tetrahydroquinoline (13): detailed procedure"),
    ("table2_entry_3", [9], "3-Ethyl-2-tosyl-5,6,7,8-tetrahydroquinoline (14): detailed procedure"),
    ("table2_entry_4", [9], "3-Methyl-4-phenyl-2-tosyl-5,6,7,8-tetrahydroquinoline (15): detailed procedure"),
    ("table2_entry_5", [9], "4-Ethyl-3-methyl-6-phenyl-2-tosylpyridine (16): detailed procedure"),
    ("table2_entry_6", [9], "Regioisomeric pyridines 17 and 18: detailed procedure"),
    ("scheme5_23_to_24", [9, 10], "1-(1-Methoxybut-2-yn-1-yl)cyclohex-1-ene (24): detailed procedure"),
    ("scheme5_24_to_26", [10], "3-Ethyl-4-methoxy-2-tosyl-5,6,7,8-tetrahydroquinoline (26): detailed procedure"),
    ("scheme6_27_to_28", [10], "1-(3-Cyclohexyl-1-methoxyprop-2-yn-1-yl)cyclohex-1-ene (28): detailed procedure"),
    ("scheme6_28_to_30", [10], "3-(Cyclohexylmethyl)-4-methoxy-2-tosyl-5,6,7,8-tetrahydroquinoline (30): detailed procedure"),
    ("scheme7_31_to_32", [10, 11], "Ethyl 3-(cyclohex-1-en-1-yl)hexa-3,4-dienoate (32): detailed procedure"),
    ("scheme7_31_to_33", [11], "3-(Cyclohex-1-en-1-yl)-N,N-dimethylhexa-3,4-dienamide (33): detailed procedure"),
    ("scheme7_33_to_34", [11], "Pyridine amide 34: detailed procedure"),
    ("scheme7_32_to_35", [11], "Pyridine ester 35: detailed procedure"),
    ("scheme8_36_to_37", [11], "Propargylic methyl ether 37: detailed procedure"),
    ("scheme8_37_to_39", [11, 12], "7-Azaindole 39: detailed procedure"),
    ("scheme8_37_to_40", [12], "7-Azaindole 40: detailed procedure"),
    ("scheme9_41_to_42", [12], "Multistage preparation of propargylic alcohol 42: detailed procedure"),
    ("scheme9_42_to_43", [12], "Allenylpyrrole 43: detailed procedure"),
    ("scheme9_43_to_44", [12], "7-Azaindole 44: detailed procedure"),
]:
    solvent_fix = None
    if short_id == "scheme7_33_to_34":
        solvent_fix = (f"{DOC}:scheme7_33_to_34:solvent", "toluene", "toluene")
    if short_id == "scheme7_32_to_35":
        solvent_fix = (f"{DOC}:scheme7_32_to_35:solvent", "toluene", "toluene")
    annotate_existing(short_id, pages, heading, solvent_fix=solvent_fix)


# Add exact core charges for the six Table 2 detailed preparations.
detail_specs = {
    "table2_entry_2": dict(page=8, substrate=("0.150 g; 1.01 mmol", None), tscn=("0.201 g; 1.11 mmol; 1.1 equiv", None), solvent=("10 mL", q(10, "mL")), ms=("0.150 g", q(0.150, "g")), dbu=("0.17 mL; 0.17 g; 1.1 mmol; 1.1 equiv", None), product=("0.226 g", q(0.226, "g"))),
    "table2_entry_1": dict(page=9, substrate=("0.040 g; 0.17 mmol; 1.0 equiv", None), tscn=("0.033 g; 0.19 mmol; 1.1 equiv", None), solvent=("1.7 mL", q(1.7, "mL")), ms=("0.045 g", q(0.045, "g")), dbu=("28 uL; 0.028 g; 0.19 mmol; 1.1 equiv", None), product=("0.048 g", q(0.048, "g"))),
    "table2_entry_3": dict(page=9, substrate=("0.092 g; 0.685 mmol; 1.0 equiv", None), tscn=("0.137 g; 0.754 mmol; 1.1 equiv", None), solvent=("7 mL", q(7, "mL")), ms=("0.100 g", q(0.100, "g")), dbu=("0.11 mL; 0.114 g; 0.754 mmol; 1.1 equiv", None), product=("0.085 g", q(0.085, "g"))),
    "table2_entry_4": dict(page=9, substrate=("0.119 g; 0.60 mmol; 1.0 equiv in 3 mL toluene", None), tscn=("0.120 g; 0.66 mmol; 1.1 equiv", None), solvent=("3 mL initial plus 3 mL substrate solution", None), ms=("0.170 g", q(0.170, "g")), dbu=("0.1 mL; 0.100 g; 0.66 mmol; 1.1 equiv", None), product=("0.093 g", q(0.093, "g"))),
    "table2_entry_5": dict(page=9, substrate=("0.763 g; 4.480 mmol; 1.0 equiv", None), tscn=("0.889 g; 4.93 mmol; 1.1 equiv", None), solvent=("45 mL", q(45, "mL")), ms=("0.780 g", q(0.780, "g")), product=("1.136 g", q(1.136, "g"))),
    "table2_entry_6": dict(page=9, substrate=("0.135 g; 1.0 mmol; 1.0 equiv", None), tscn=("0.198 g; 1.09 mmol; 1.1 equiv", None), solvent=("10 mL", q(10, "mL")), ms=("0.140 g", q(0.140, "g")), dbu=("0.16 mL; 1.09 mmol; 1.1 equiv", None), product17=("0.078 g", q(0.078, "g")), product18=("0.042 g", q(0.042, "g"))),
}
updates_by_short = {r["reaction_id"].split(":", 1)[1]: r for r in updated_reactions}
for short_id, spec in detail_specs.items():
    r = updates_by_short[short_id]
    ev = [P[spec["page"]]]
    base = f"{DOC}:{short_id}"
    charge(r, base + ":substrate", ev, mass=spec["substrate"])
    charge(r, base + ":tscn", ev, mass=spec["tscn"])
    charge(r, base + ":solvent", ev, volume=spec["solvent"])
    charge(r, base + ":ms", ev, mass=spec["ms"])
    if "dbu" in spec:
        charge(r, base + ":dbu", ev, volume=spec["dbu"])
    if short_id == "table2_entry_6":
        charge(r, "t2e6_product_17", ev, mass=spec["product17"])
        charge(r, "t2e6_product_18", ev, mass=spec["product18"])
    else:
        product_pid = next(p["participant_id"] for p in r["participants"] if p["role"] == "product")
        charge(r, product_pid, ev, mass=spec["product"])


# New molecules first encountered in pages 7-12.
add_molecule("S1", "propargylic ether S1", [P[8]], unresolved_reason="Compound S1 is only referenced by label in the bounded page; its drawing is in unavailable cited material/SI.")
add_molecule("S2", "propargylic substrate S2", [P[8]], unresolved_reason="Compound S2 is only referenced by label in the bounded page; its complete structure is not shown here.")
add_molecule("S3", "propargylic alcohol S3", [P[8]], unresolved_reason="Compound S3 is only referenced by label in the bounded page; its complete structure is not shown here.")
add_molecule("S4", "6-methylhept-6-en-4-yn-3-yl benzoate", [P[8]],
             smiles="CCC(OC(=O)c1ccccc1)C#CC(=C)C")

product_names = {
    "45": ("4-ethyl-3-methyl-6-phenyl-2-methylpyridine", "n1c(C)c(C)c(CC)cc1-c1ccccc1"),
    "46": ("4-ethyl-3-methyl-6-phenyl-2-propylpyridine", "n1c(CCC)c(C)c(CC)cc1-c1ccccc1"),
    "47": ("4-ethyl-3-methyl-2,6-diphenylpyridine", "n1c(-c2ccccc2)c(C)c(CC)cc1-c1ccccc1"),
    "48": ("4-ethyl-3-methyl-6-phenylpicolinonitrile", "n1c(C#N)c(C)c(CC)cc1-c1ccccc1"),
    "49": ("2-ethoxy-4-ethyl-3-methyl-6-phenylpyridine", "n1c(OCC)c(C)c(CC)cc1-c1ccccc1"),
    "50": ("3-ethyl-2,4-dimethyl-5,6,7,8-tetrahydroquinoline", "CCc1c(C)nc2c(c1C)CCCC2"),
    "51": ("3-ethyl-4-methoxy-2-methyl-5,6,7,8-tetrahydroquinoline", "COc1c(CC)c(C)nc2CCCCc12"),
    "52": ("4-methyl-3-(3-phenylpropyl)-5,6,7,8-tetrahydroquinolin-2(1H)-one", "Cc1c(CCCc2ccccc2)c2CCCCc2[nH]c1=O"),
}
for label, (name, smiles) in product_names.items():
    ev = [P[7], P[12]] if label == "45" else [P[7]]
    add_molecule(label, name, ev, smiles=smiles)

new_reagent_specs = {
    "phmgbr": ("phenylmagnesium bromide", None, [P[8]], ["PhMgBr"]),
    "etmgbr": ("ethylmagnesium bromide", None, [P[8]], ["EtMgBr"]),
    "ch2cl2": ("dichloromethane", "ClCCl", [P[8]], ["CH2Cl2"]),
    "benzoyl_chloride": ("benzoyl chloride", "O=C(Cl)c1ccccc1", [P[8]], None),
    "prmgbr": ("propylmagnesium bromide", None, [P[7]], ["PrMgBr"]),
    "phmgcl": ("phenylmagnesium chloride", None, [P[7]], ["PhMgCl"]),
    "nacn": ("sodium cyanide", "[Na+].[C-]#N", [P[7]], ["NaCN"]),
    "18c6": ("18-crown-6", "C1COCCOCCOCCOCCOCCO1", [P[7]], ["18-crown-6"]),
    "dmf": ("N,N-dimethylformamide", "CN(C)C=O", [P[7]], ["DMF"]),
    "ethanol": ("ethanol", "CCO", [P[7]], ["EtOH"]),
    "memgcl": ("methylmagnesium chloride", None, [P[7]], ["CH3MgCl"]),
    "water": ("water", "O", [P[7], P[12]], ["H2O"]),
    "nh4cl": ("ammonium chloride", "[Cl-].[NH4+]", [P[12]], ["NH4Cl"]),
    "etoac": ("ethyl acetate", "CCOC(C)=O", [P[12]], ["EtOAc"]),
    "nacl": ("sodium chloride", "[Na+].[Cl-]", [P[12]], ["NaCl"]),
    "mgso4": ("magnesium sulfate", "[Mg+2].[O-]S(=O)(=O)[O-]", [P[12]], ["MgSO4"]),
    "silica": ("silica gel", None, [P[12]], None),
    "hexanes": ("hexanes (composition unreported)", None, [P[12]], None),
}
for label, (name, smiles, ev, aliases) in new_reagent_specs.items():
    add_named(label, name, ev, smiles, aliases)


def set_classification(r, ev, formed=(), unresolved=False, extra_changes=()):
    c = r["classification"]
    c["mapped_reaction_smiles"] = fact(evidence=ev)
    c["model_votes"] = []
    c["limitations"] = [
        "No complete atom mapping was asserted; element-pair changes are partial structure derivations from the source drawing."
    ]
    if unresolved:
        c["bond_changes"] = []
        c["labels"] = F("Insufficient unique structures for bond classification", [], ev, "inference",
                          derivation="One or more reactant structures are only identified by source label.")
        c["coverage"] = F("unresolved", "unresolved", ev, "inference",
                            derivation="Source does not support a complete mapped transformation in the bounded pages.")
        return
    changes = []
    labels = {}
    for a, b, before, after, change_type in list(formed) + list(extra_changes):
        cid = f"b{len(changes)+1}"
        change = {"change_id": cid, "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
                  "before_order": before, "after_order": after, "change_type": change_type}
        if change_type.startswith("hydrogen_"):
            change["hydrogen_count"] = 1
        changes.append(F("net bond change inferred from the depicted transformation", change, ev,
                         "structure_derivation", derivation="Element-pair change inferred from the source drawing; atom correspondence remains unresolved."))
        if change_type in ("formed", "order_increased", "hydrogen_gain"):
            label = bond_label(a, b, after)
            labels.setdefault(label["code"], {**label, "supporting_change_ids": []})["supporting_change_ids"].append(cid)
    c["bond_changes"] = changes
    c["labels"] = F("formed-bond labels inferred from source drawings", list(labels.values()), ev,
                      "structure_derivation", derivation="Labels follow bond-formation-20261016; no complete atom mapping is claimed.")
    c["coverage"] = F("partial", "partial", ev, "inference",
                        derivation="Only visually defensible net bond changes are recorded.")


def reaction(short_id, ev, categories, name, transformation, location, formed=(), unresolved=False):
    r = new_reaction(f"{DOC}:{short_id}", DOC, ev)
    r["record_type"] = "specific_experiment"
    r["experiment_categories"] = F(categories, categories, ev)
    r["location"] = F(location, location, ev)
    r["reaction_name"] = F(name, name, ev)
    r["transformation_type"] = F(transformation, transformation, ev, "inference",
                                  derivation="Plain-language transformation inferred from source heading/drawing.")
    r["general_scheme_ids"] = []
    r["reaction_scale"] = fact(evidence=ev)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": [
                       "Machine record: structures, roles, and partial bond classification require expert review."
                   ]}
    r["completeness"] = "partial"
    set_classification(r, ev, formed, unresolved)
    new_reactions.append(r)
    return r


def participant(r, pid, molecule, role, ev, stage="reaction", *, mass=None,
                amount=None, volume=None, equivalents=None, mol_percent=None,
                product_role=None, reference=None):
    molecule_id = mid(molecule)
    if molecule_id not in old_molecules and molecule_id not in new_molecules:
        raise KeyError(f"unknown molecule {molecule_id}")
    p = new_participant(pid, molecule_id, role, ev, stage)
    if product_role:
        p["product_role"] = F(product_role, product_role, ev)
    for key, value in (("mass", mass), ("amount", amount), ("volume", volume),
                       ("equivalents", equivalents), ("mol_percent", mol_percent)):
        if value:
            raw, norm = value
            p["charges"][key] = F(raw, norm, ev)
    if reference:
        p["equivalence_reference"] = F(reference, mid(reference), ev)
    r["participants"].append(p)
    return p


def step(r, sid, ev, pids, operation, temp=None, time=None, stage="reaction"):
    s = new_step(sid, ev, stage)
    s["participant_ids"] = list(pids)
    s["operation"] = F(operation, operation, ev)
    if temp:
        s["temperature"] = F(temp[0], temp[1], ev)
    if time:
        s["time"] = F(time[0], time[1], ev)
    r["steps"].append(s)
    return s


def outcome(r, oid, target, raw, norm, ev, basis, scope=None):
    o = new_outcome(oid, "yield", [mid(target)], ev)
    o["product_role"] = F("main product", "main product", ev)
    o["value"] = F(raw, norm, ev)
    o["yield_type"] = F("isolated yield", "isolated", ev)
    o["measurement_basis"] = F(basis, "isolated", ev)
    o["scope"] = F(scope, scope, ev) if scope else fact(evidence=ev)
    r["outcomes"].append(o)
    return o


CC = ("C", "C", 0, 1, "formed")
CO = ("C", "O", 0, 1, "formed")
DBCO = ("C", "O", 0, 2, "formed")
NH = ("N", "H", 0, 1, "hydrogen_gain")


# Four preparations first reported in the bounded experimental section.
r = reaction("prep_10", [P[8]], ["other"], "Preparation of vinylallene 10",
             "copper-mediated propargylic substitution to an allene",
             "PDF page 8, (1-(cyclohex-1-en-1-yl)propa-1,2-dien-1-yl)benzene (10)", [CC], unresolved=True)
participant(r, "prep10_s1", "S1", "reactant", [P[8]], mass=("0.200 g", q(0.200, "g")), amount=("1.33 mmol", q(1.33, "mmol")), equivalents=("1.0 equiv", q(1.0, "equiv")), reference="S1")
participant(r, "prep10_cubr", "cubr", "reagent", [P[8]], mass=("0.039 g", q(0.039, "g")), amount=("0.26 mmol", q(0.26, "mmol")), equivalents=("0.2 equiv", q(0.2, "equiv")), reference="S1")
participant(r, "prep10_phmgbr", "phmgbr", "reagent", [P[8]], volume=("1.15 mL of 3.0 M solution", q(1.15, "mL")), amount=("3.32 mmol", q(3.32, "mmol")), equivalents=("2.5 equiv", q(2.5, "equiv")), reference="S1")
participant(r, "prep10_ether", "ether", "solvent", [P[8]], volume=("8 mL", q(8, "mL")))
participant(r, "prep10_product", "10", "product", [P[8]], mass=("0.189 g", q(0.189, "g")), product_role="main product")
step(r, "cool_and_add", [P[8]], ["prep10_s1", "prep10_cubr", "prep10_phmgbr", "prep10_ether"], "add phenylmagnesium bromide dropwise", ("-30 °C", q(-30, "°C")))
step(r, "warm_and_hold", [P[8]], [], "allow to warm to room temperature and stir", ("rt", "rt"), ("15 h", q(15, "h")))
step(r, "purification", [P[8]], [], "silica-gel column chromatography with hexanes", stage="purification")
outcome(r, "yield_10", "10", "72%", q(72, "%"), [P[8]], "isolated after column chromatography")
r["procedure"]["original_text"] = F("CuBr (0.2 equiv), S1 (1.0 equiv), and Et2O; PhMgBr (2.5 equiv) at -30 °C; warm to rt for 15 h; column chromatography", "detailed experimental procedure", [P[8]])

r = reaction("prep_11", [P[8]], ["other"], "Preparation of vinylallene 11",
             "copper-mediated propargylic substitution to an allene",
             "PDF page 8, (E)-(3-ethylpenta-1,3,4-trien-1-yl)benzene (11)", [CC], unresolved=True)
participant(r, "prep11_s2", "S2", "reactant", [P[8]], mass=("0.400 g", q(0.400, "g")), amount=("2.32 mmol", q(2.32, "mmol")), equivalents=("1.0 equiv", q(1.0, "equiv")), reference="S2")
participant(r, "prep11_cubr", "cubr", "reagent", [P[8]], mass=("0.066 g", q(0.066, "g")), amount=("0.464 mmol", q(0.464, "mmol")), equivalents=("0.2 equiv", q(0.2, "equiv")), reference="S2")
participant(r, "prep11_etmgbr", "etmgbr", "reagent", [P[8]], volume=("1.54 mL of 3.0 M solution", q(1.54, "mL")), amount=("4.64 mmol", q(4.64, "mmol")), equivalents=("2.0 equiv", q(2.0, "equiv")), reference="S2")
participant(r, "prep11_thf", "thf", "solvent", [P[8]], volume=("20 mL", q(20, "mL")))
participant(r, "prep11_product", "11", "product", [P[8]], mass=("0.317 g", q(0.317, "g")), product_role="main product")
step(r, "addition", [P[8]], ["prep11_s2", "prep11_cubr", "prep11_etmgbr", "prep11_thf"], "add ethylmagnesium bromide dropwise", ("0 °C", q(0, "°C")), ("1 min", q(1, "min")))
step(r, "warm_and_hold", [P[8]], [], "warm to room temperature and stir", ("rt", "rt"), ("45 min", q(45, "min")))
step(r, "purification", [P[8]], [], "silica-gel column chromatography with hexanes", stage="purification")
outcome(r, "yield_11", "11", "80%", q(80, "%"), [P[8]], "isolated after column chromatography")

r = reaction("prep_s4", [P[8]], ["other"], "Preparation of propargylic benzoate S4",
             "benzoylation of a propargylic alcohol",
             "PDF page 8, 6-methylhept-6-en-4-yn-3-yl benzoate (S4)", [CO], unresolved=True)
participant(r, "preps4_s3", "S3", "reactant", [P[8]], mass=("0.314 g", q(0.314, "g")), amount=("2.53 mmol", q(2.53, "mmol")), equivalents=("1.0 equiv", q(1.0, "equiv")), reference="S3")
participant(r, "preps4_dmap", "dmap", "reagent", [P[8]], mass=("0.031 g", q(0.031, "g")), amount=("0.253 mmol", q(0.253, "mmol")), equivalents=("0.1 equiv", q(0.1, "equiv")), reference="S3")
participant(r, "preps4_et3n", "et3n", "reagent", [P[8]], volume=("0.49 mL", q(0.49, "mL")), mass=("0.36 g", q(0.36, "g")), amount=("3.5 mmol", q(3.5, "mmol")), equivalents=("1.4 equiv", q(1.4, "equiv")), reference="S3")
participant(r, "preps4_bzcl", "benzoyl_chloride", "reagent", [P[8]], volume=("0.41 mL", q(0.41, "mL")), mass=("0.50 g", q(0.50, "g")), amount=("3.5 mmol", q(3.5, "mmol")), equivalents=("1.4 equiv", q(1.4, "equiv")), reference="S3")
participant(r, "preps4_ch2cl2", "ch2cl2", "solvent", [P[8]], volume=("25 mL", q(25, "mL")))
participant(r, "preps4_product", "S4", "product", [P[8]], mass=("0.570 g", q(0.570, "g")), product_role="main product")
step(r, "benzoylation", [P[8]], ["preps4_s3", "preps4_dmap", "preps4_et3n", "preps4_bzcl", "preps4_ch2cl2"], "add triethylamine then benzoyl chloride and stir", ("rt", "rt"), ("5 h", q(5, "h")))
step(r, "purification", [P[8]], [], "silica-gel column chromatography, 0-50% CH2Cl2/hexanes", stage="purification")
outcome(r, "yield_s4", "S4", "99%", q(99, "%"), [P[8]], "isolated after column chromatography")

r = reaction("prep_12", [P[8]], ["other"], "Preparation of vinylallene 12",
             "copper-mediated SN2'-type substitution of a propargylic benzoate",
             "PDF page 8, 3-ethyl-2-methylhepta-1,3,4-triene (12)", [CC])
participant(r, "prep12_s4", "S4", "reactant", [P[8]], mass=("0.300 g", q(0.300, "g")), amount=("1.31 mmol", q(1.31, "mmol")), equivalents=("1.0 equiv", q(1.0, "equiv")), reference="S4")
participant(r, "prep12_libr", "libr", "reagent", [P[8]], mass=("0.684 g", q(0.684, "g")), amount=("7.88 mmol", q(7.88, "mmol")), equivalents=("6.0 equiv", q(6.0, "equiv")), reference="S4")
participant(r, "prep12_cubr", "cubr", "reagent", [P[8]], mass=("1.13 g", q(1.13, "g")), amount=("7.88 mmol", q(7.88, "mmol")), equivalents=("6.0 equiv", q(6.0, "equiv")), reference="S4")
participant(r, "prep12_etmgbr", "etmgbr", "reagent", [P[8]], volume=("2.62 mL of 3.0 M solution", q(2.62, "mL")), amount=("7.88 mmol", q(7.88, "mmol")), equivalents=("6.0 equiv", q(6.0, "equiv")), reference="S4")
participant(r, "prep12_thf", "thf", "solvent", [P[8]], volume=("10 mL plus 2 mL substrate solution and 2 mL rinse", None))
participant(r, "prep12_product", "12", "product", [P[8]], mass=("0.135 g", q(0.135, "g")), product_role="main product")
step(r, "organocopper_formation", [P[8]], ["prep12_libr", "prep12_cubr", "prep12_etmgbr", "prep12_thf"], "add EtMgBr to LiBr/CuBr and stir", ("-60 °C", q(-60, "°C")), ("10 min", q(10, "min")))
step(r, "substrate_addition", [P[8]], ["prep12_s4"], "add S4 solution dropwise", ("-60 °C", q(-60, "°C")), ("1 min", q(1, "min")))
step(r, "warm_and_hold", [P[8]], [], "allow to warm to room temperature and stir", ("rt", "rt"), ("1 h", q(1, "h")))
step(r, "purification", [P[8]], [], "silica-gel column chromatography with pentane", stage="purification")
outcome(r, "yield_12", "12", "76%", q(76, "%"), [P[8]], "isolated after column chromatography")


# Table 3: every row is a distinct experiment. The table footnote is the
# source for the isolated-yield/column-purification measurement basis.
table3 = {
    1: dict(substrate="16", product="45", reagent="memgbr", eq=3.0, solvent="thf", temp="reflux", time=(15, "h"), y=(90, 90), ev=[P[7], P[12]], formed=[CC]),
    2: dict(substrate="16", product="46", reagent="prmgbr", eq=3.0, solvent="thf", temp="reflux", time=(13, "h"), y=(77, 82), ev=[P[7]], formed=[CC]),
    3: dict(substrate="16", product="47", reagent="phmgcl", eq=6.0, solvent="thf", temp="reflux", time=(16, "h"), y=(87, 87), ev=[P[7]], formed=[CC]),
    4: dict(substrate="16", product="48", reagent="nacn", eq=5.0, extra=("18c6", 5.0), solvent="dmf", temp="reflux", time=(48, "h"), y=(86, 86), ev=[P[7]], formed=[CC]),
    5: dict(substrate="16", product="49", reagent="kotbu", eq=2.4, extra=("ethanol", 2.4), solvent="thf", temp="rt", time=(40, "min"), y=(91, 91), ev=[P[7]], formed=[CO]),
    6: dict(substrate="7", product="50", reagent="memgbr", eq=6.0, solvent="thf", temp=70, time=(1.5, "h"), y=(53, 53), ev=[P[7]], formed=[CC]),
    7: dict(substrate="7", product="50", reagent="memgbr", eq=12.0, solvent="thf", temp="rt", time=(24, "h"), y=(54, 54), ev=[P[7]], formed=[CC]),
    8: dict(substrate="26", product="51", reagent="memgcl", eq=12.0, solvent="thf", temp="rt", time=(16, "h"), y=(78, 78), ev=[P[7]], formed=[CC]),
    9: dict(substrate="13", product="52", reagent="kotbu", eq=6.0, extra=("water", 2.0), solvent="thf", temp="reflux", time=(16, "h"), y=(75, 75), ev=[P[7]], formed=[DBCO, NH]),
}

for entry, spec in table3.items():
    ev = spec["ev"]
    r = reaction(f"table3_entry_{entry}", ev, ["derivatization", "substrate_scope"],
                 "Ipso substitution reaction of a 2-sulfonylpyridine",
                 "nucleophilic ipso substitution of the 2-sulfonyl group",
                 f"PDF page 7, Table 3, entry {entry}" + ("; detailed procedure begins on page 12" if entry == 1 else ""),
                 spec["formed"])
    participant(r, f"t3e{entry}_substrate", spec["substrate"], "reactant", ev,
                equivalents=("1.0 equiv", q(1.0, "equiv")), reference=spec["substrate"])
    participant(r, f"t3e{entry}_reagent", spec["reagent"], "reagent", ev,
                equivalents=(f"{spec['eq']:g} equiv", q(spec["eq"], "equiv")), reference=spec["substrate"])
    if spec.get("extra"):
        extra, eq = spec["extra"]
        participant(r, f"t3e{entry}_extra", extra, "reagent", ev,
                    equivalents=(f"{eq:g} equiv", q(eq, "equiv")), reference=spec["substrate"])
    participant(r, f"t3e{entry}_solvent", spec["solvent"], "solvent", ev)
    participant(r, f"t3e{entry}_product", spec["product"], "product", ev, product_role="main product")
    r["solvent_conditions"]["composition"] = F(spec["solvent"], [spec["solvent"]], ev)
    temp = (spec["temp"], spec["temp"]) if isinstance(spec["temp"], str) else (f"{spec['temp']} °C", q(spec["temp"], "°C"))
    time_num, time_unit = spec["time"]
    step(r, "reaction", ev, [p["participant_id"] for p in r["participants"] if p["role"] != "product"],
         "stir or heat under the Table 3 conditions", temp, (f"{time_num:g} {time_unit}", q(time_num, time_unit)))
    step(r, "purification", [P[7]], [], "column chromatography (Table 3 footnote)", stage="purification")
    ymin, ymax = spec["y"]
    if ymin == ymax:
        yraw, ynorm = f"{ymin}%", q(ymin, "%")
    else:
        yraw, ynorm = f"{ymin}-{ymax}%", qr(ymin, ymax, "%", f"{ymin}-{ymax}%")
    outcome(r, f"yield_{spec['product']}", spec["product"], yraw, ynorm, [P[7]],
            "isolated yield of product purified by column chromatography (Table 3 footnote)")
    r["procedure"]["original_text"] = F(
        f"Table 3 entry {entry}: {spec['reagent']} ({spec['eq']:g} equiv), {spec['solvent']}, {spec['temp']}, {time_num:g} {time_unit}",
        f"Table 3 entry {entry} conditions", [P[7]])

    if entry == 1:
        # Page 12 supplies exact charges and the beginning of the detailed workup.
        charge(r, "t3e1_substrate", [P[12]], mass=("0.150 g", q(0.150, "g")), amount=("0.43 mmol", q(0.43, "mmol")))
        charge(r, "t3e1_reagent", [P[12]], volume=("0.50 mL of 2.54 M solution", q(0.50, "mL")), amount=("1.28 mmol", q(1.28, "mmol")))
        charge(r, "t3e1_solvent", [P[12]], volume=("7 mL", q(7, "mL")))
        charge(r, "t3e1_product", [P[12]], mass=("0.081 g", q(0.081, "g")))
        participant(r, "t3e1_nh4cl", "nh4cl", "workup_agent", [P[12]], "workup", volume=("4 mL saturated aqueous solution", q(4, "mL")))
        participant(r, "t3e1_water", "water", "workup_agent", [P[12]], "workup", volume=("3 mL", q(3, "mL")))
        participant(r, "t3e1_etoac", "etoac", "workup_agent", [P[12]], "workup", volume=("three 10 mL portions", q(30, "mL")))
        participant(r, "t3e1_nacl", "nacl", "workup_agent", [P[12]], "workup", volume=("30 mL saturated aqueous solution", q(30, "mL")))
        participant(r, "t3e1_mgso4", "mgso4", "workup_agent", [P[12]], "workup")
        participant(r, "t3e1_silica", "silica", "purification_agent", [P[12]], "purification", mass=("9 g", q(9, "g")))
        participant(r, "t3e1_hexanes", "hexanes", "purification_agent", [P[12]], "purification")
        step(r, "detailed_workup", [P[12]], ["t3e1_nh4cl", "t3e1_water", "t3e1_etoac", "t3e1_nacl", "t3e1_mgso4"],
             "cool, quench, extract, wash, dry, filter, and concentrate", stage="workup")
        step(r, "detailed_purification", [P[12]], ["t3e1_silica", "t3e1_hexanes"],
             "column chromatography, 0-5% EtOAc/hexanes", stage="purification")


# Source context keeps independent instrumentation/characterization material out
# of the reaction count while documenting that it was inspected.
source_context = [
    {"context_id": f"{DOC}:pages7_12_characterization", "document_id": DOC,
     "kind": "characterization_data",
     "evidence_ids": [P[i] for i in range(7, 13)],
     "note": "IR, 1H/13C NMR, HRMS, melting points, nOe values, crude weights and spectral peak lists accompany the preparations. They corroborate product characterization but are not independent synthesis experiments."},
    {"context_id": f"{DOC}:page7_instrumentation", "document_id": DOC,
     "kind": "general_instrumentation",
     "evidence_ids": [P[7]],
     "note": "Instrumentation and chemical-shift reference statements were inspected as general analytical context and were not inherited as reaction conditions."},
]


# Assemble an incremental bundle. Existing referenced molecules are included
# unchanged so that the standalone incoming bundle validates atomically.
all_reactions = updated_reactions + new_reactions
referenced = {p["molecule_id"] for r in all_reactions for p in r["participants"]}
all_molecules = []
for molecule_id in sorted(referenced):
    if molecule_id in new_molecules:
        all_molecules.append(new_molecules[molecule_id])
    else:
        all_molecules.append(copy.deepcopy(old_molecules[molecule_id]))

SUBMITTED["molecules"] = all_molecules
SUBMITTED["reactions"] = all_reactions
SUBMITTED["source_context"] = source_context
write(OUT / "submitted_pages_0007_0012.json", SUBMITTED)
write(OUT / "submitted.json", SUBMITTED)
write(OUT / "raw_model" / "stereo_checks_pages_0007_0012.json", stereo_checks)

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 7, "state": "attempted", "regions": [
            {"region": "Table 3 (entries 1-9 and footnote a)", "finding": "Nine distinct ipso-substitution experiments were transcribed row-by-row; footnote a assigns isolated yields after column chromatography.", "evidence_ids": [P[7]]},
            {"region": "Instrumentation and compound 2 procedure", "finding": "General analytical instrumentation was retained as source context. The detailed preparation of compound 2 was merged into the existing Scheme 3 record; spectral data were not made into reactions.", "evidence_ids": [P[7]]},
            {"region": "compound 3 procedure beginning at lower right", "finding": "The procedure continues onto page 8 and was merged with the existing 2-to-3 Scheme 3 record.", "evidence_ids": [P[7]]}
        ], "record_ids": [f"{DOC}:table3_entry_{i}" for i in range(1, 10)] + [f"{DOC}:scheme3_1_to_2", f"{DOC}:scheme3_2_to_3"],
         "unresolved": ["Table 3 entries 2-9 have no detailed procedures within pages 7-12; only the table conditions and footnote measurement basis are claimed here."]},
        {"page": 8, "state": "attempted", "regions": [
            {"region": "continuation of compound 3 and procedures for 7, 10, 11, S4 and 12", "finding": "Compound 3 and product 7 procedures were merged into existing records. Preparations of 10, 11, S4 and 12 were first encountered here and became four new experiments.", "evidence_ids": [P[8]]},
            {"region": "spectral/HRMS paragraphs", "finding": "Analytical peak lists were classified as characterization context, not synthesis experiments.", "evidence_ids": [P[8]]}
        ], "record_ids": [f"{DOC}:scheme3_2_to_3", f"{DOC}:table2_entry_2", f"{DOC}:prep_10", f"{DOC}:prep_11", f"{DOC}:prep_s4", f"{DOC}:prep_12"],
         "unresolved": ["S1, S2 and S3 are referenced by label without complete structures in the bounded pages, so their SMILES and stereochemical identities remain unresolved."]},
        {"page": 9, "state": "attempted", "regions": [
            {"region": "procedures for products 13-18", "finding": "These are detailed counterparts of Table 2 entries 1 and 3-6; evidence and explicit charges were merged without creating duplicate experiments.", "evidence_ids": [P[9]]},
            {"region": "compound 24 procedure beginning at lower right", "finding": "The procedure continues onto page 10 and maps to the existing Scheme 5 record.", "evidence_ids": [P[9]]}
        ], "record_ids": [f"{DOC}:table2_entry_{i}" for i in (1, 3, 4, 5, 6)] + [f"{DOC}:scheme5_23_to_24"],
         "unresolved": ["The previously recorded Table 2 versus detailed-paragraph yield conflicts for products 14 and 16 remain unresolved; neither value was selected."]},
        {"page": 10, "state": "attempted", "regions": [
            {"region": "continuation of 24 and procedures for 26, 28 and 30", "finding": "All are detailed counterparts of Schemes 5 and 6 and were merged into those existing records.", "evidence_ids": [P[10]]},
            {"region": "compound 32 procedure beginning at lower right", "finding": "The Johnson rearrangement procedure continues onto page 11 and maps to the existing Scheme 7 record.", "evidence_ids": [P[10]]}
        ], "record_ids": [f"{DOC}:scheme5_23_to_24", f"{DOC}:scheme5_24_to_26", f"{DOC}:scheme6_27_to_28", f"{DOC}:scheme6_28_to_30", f"{DOC}:scheme7_31_to_32"],
         "unresolved": ["Transient allenyl ethers 25 and 29 were not isolated; their axial configuration remains unspecified and not fully encodable."]},
        {"page": 11, "state": "attempted", "regions": [
            {"region": "continuation of 32 and procedures for 33-35 and 37", "finding": "Detailed procedures were merged with the existing Scheme 7 and Scheme 8 records. The explicit toluene solvent for 34 and 35 corrects the earlier scheme-only extraction while retaining an unresolved note.", "evidence_ids": [P[11]]},
            {"region": "compound 39 procedure beginning at lower right", "finding": "The procedure continues onto page 12 and maps to the existing Scheme 8 record.", "evidence_ids": [P[11]]}
        ], "record_ids": [f"{DOC}:scheme7_31_to_32", f"{DOC}:scheme7_31_to_33", f"{DOC}:scheme7_33_to_34", f"{DOC}:scheme7_32_to_35", f"{DOC}:scheme8_36_to_37", f"{DOC}:scheme8_37_to_39"],
         "unresolved": ["The Scheme 8 versus detailed-paragraph loading conflict for compound 36 remains explicit in the existing record."]},
        {"page": 12, "state": "attempted", "regions": [
            {"region": "continuation of 39 and procedures for 40, 42, 43 and 44", "finding": "These are detailed counterparts of Schemes 8 and 9 and were merged into the existing experiment records; characterization-only spectra remain source context.", "evidence_ids": [P[12]]},
            {"region": "compound 45 procedure at lower right", "finding": "This is the detailed counterpart of Table 3 entry 1; exact charges, ordered reaction/workup/purification stages and product mass were merged into the same new record.", "evidence_ids": [P[12]]}
        ], "record_ids": [f"{DOC}:scheme8_37_to_39", f"{DOC}:scheme8_37_to_40", f"{DOC}:scheme9_41_to_42", f"{DOC}:scheme9_42_to_43", f"{DOC}:scheme9_43_to_44", f"{DOC}:table3_entry_1"],
         "unresolved": ["The Scheme 9 versus detailed-paragraph yield conflict for product 44 (35% versus 34%) remains unresolved; neither value was selected.", "The compound 45 paragraph continues on page 13; only the portion visible on page 12 plus Table 3 is claimed in this batch."]},
    ],
}
write(OUT / "coverage" / "pages_0007_0012.json", coverage)

raw_summary = {
    "batch_id": "allene-stereo-v31-200-20261017",
    "document_id": DOC,
    "page_range": [7, 12],
    "visual_pages_inspected": [7, 8, 9, 10, 11, 12],
    "semantic_accuracy_verified": False,
    "new_reaction_count": len(new_reactions),
    "updated_existing_reaction_count": len(updated_reactions),
    "new_molecule_count": len(new_molecules),
    "decisions": {
        "table3": "nine independent row-level derivatization/substrate-scope experiments",
        "preparations_10_11_S4_12": "four new preparative experiments first described in page-8 text",
        "detailed_duplicates": "merged evidence into existing Scheme/Table records; no duplicate experiments created",
        "characterization": "spectral peak lists and instrumentation retained in source_context",
        "si": "not associated; source mapping says unestablished",
    },
    "unresolved": [
        "S1, S2 and S3 structures are not shown in the bounded pages.",
        "Table 3 entries 2-9 lack detailed experimental paragraphs within pages 7-12.",
        "All machine records remain needs_review; validation is not chemical verification.",
    ],
}
write(OUT / "raw_model" / "raw_extraction_pages_0007_0012.json", raw_summary)
print(json.dumps({
    "new_reactions": len(new_reactions),
    "updated_reactions": len(updated_reactions),
    "bundle_reactions": len(all_reactions),
    "bundle_molecules": len(all_molecules),
}, ensure_ascii=False))
