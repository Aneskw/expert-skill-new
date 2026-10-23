#!/usr/bin/env python3
"""Extend the v3.1 source bundle with original PDF pages 13--18.

This script merges detailed procedures into experiments already inventoried on
earlier pages, creates only the newly encountered preparation of allenylic
acetate 9, and registers the separated E/Z identities of 13a--c.  Spectral
characterization remains source context rather than a synthetic experiment.
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
T = {p: f"{DOC}:page{p}:text" for p in range(13, 19)}
I = {p: f"{DOC}:page{p}:image" for p in range(13, 19)}


def ev(*pages):
    out = []
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
bundle["documents"] = dataset["documents"]
molecules = {m["molecule_id"]: m for m in bundle["molecules"]}
reactions = {r["reaction_id"]: r for r in bundle["reactions"]}


def rid(short):
    return f"{DOC}:{short}"


def getr(short):
    return reactions[rid(short)]


def union_evidence(obj, evidence):
    obj["evidence_ids"] = list(dict.fromkeys(obj.get("evidence_ids", []) + list(evidence)))


def add_occurrence(mid, label, evidence):
    m = molecules[mid]
    if not any(set(o["evidence_ids"]) == set(evidence) for o in m["occurrences"]):
        m["occurrences"].append({
            "document_id": DOC,
            "compound_label": F(label, label, evidence),
            "evidence_ids": list(evidence),
        })


def add_molecule(mid, label, name, evidence, *, smiles=None,
                 chirality="achiral", geometry="not_applicable",
                 representation=None, assignment=None, mixture=None,
                 note=None):
    if mid in molecules:
        add_occurrence(mid, label, evidence)
        return mid
    m = new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, evidence)
    m["names"]["original"] = F(label, label, evidence)
    m["names"]["english"] = F(name, name, evidence)
    m["names"]["abbreviations"] = F([label], [label], evidence)
    m["names"]["aliases"] = missing(evidence, "No additional alias is reported on pages 13--18.")
    m["cas_numbers"] = missing(evidence, "CAS number is not reported on pages 13--18 and was not guessed.")
    m["generic_template"] = None
    if smiles:
        diag = stereo_inspect(smiles)
        canonical = diag["canonical_isomeric_smiles"]
        m["representation_type"] = "smiles"
        m["smiles"] = F(smiles, canonical, evidence, "structure_derivation",
                          derivation="Connectivity transcribed from an unambiguous source name; toolkit parsing is not identity verification.")
        if representation is None:
            representation = "partial" if any(x["specified"] == "Unspecified" for x in diag["potential_stereo"]) else "encoded"
        m["stereochemistry"]["isomeric_smiles"] = F(canonical, canonical, evidence, "structure_derivation", derivation="Canonical isomeric SMILES generated from the conservatively transcribed graph.")
        m["stereochemistry"]["identity_key"] = F(diag["stereo_graph_key"], diag["stereo_graph_key"], evidence, "structure_derivation", derivation="Toolkit graph key preserves the encoded isomeric graph; it is not source identity verification.")
        m["stereochemistry"]["atom_centers"] = F(diag["atom_centers"], diag["atom_centers"], evidence, "structure_derivation", derivation="Standard stereo graph inspection; no source-silent configuration is inferred.")
        m["stereochemistry"]["double_bonds"] = F(diag["double_bonds"], diag["double_bonds"], evidence, "structure_derivation", derivation="Standard stereo graph inspection; no source-silent E/Z assignment is inferred.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing(evidence, "The source structure/stereochemical identity was retained by label; exact SMILES was not safely encoded.")
        m["stereochemistry"]["isomeric_smiles"] = missing(evidence)
        m["stereochemistry"]["identity_key"] = missing(evidence)
        m["stereochemistry"]["atom_centers"] = F(chirality, chirality, evidence, "structure_derivation",
                                                     derivation="Original page drawing and name inspected; no configuration was invented.")
        m["stereochemistry"]["double_bonds"] = F(geometry, geometry, evidence, "structure_derivation",
                                                     derivation="Original page assignment inspected; E/Z identities are kept separate.")
        representation = representation or "not_encoded"
    st = m["stereochemistry"]
    st["chirality_status"] = F(chirality, chirality, evidence, "structure_derivation",
                                  derivation="Unspecified configuration is not treated as achiral or racemic.")
    st["geometry_status"] = F(geometry, geometry, evidence, "structure_derivation", derivation="Geometry status follows only the original page name/drawing and explicit assignment text.")
    st["representation_status"] = F(representation, representation, evidence, "inference",
                                       derivation="Compares the machine representation with the source stereochemical content.")
    st["assignment"] = F(assignment, assignment, evidence) if assignment else missing(evidence)
    st["mixture"] = F(mixture, mixture, evidence) if mixture else missing(evidence)
    st["axial"] = missing(evidence, "No axial stereochemistry is applicable or reported for this identity.", status="not_applicable")
    m["identity_note"] = F(note or name, note or name, evidence)
    molecules[mid] = m
    bundle["molecules"].append(m)
    return mid


def simple(mid, label, name, pages, smiles):
    return add_molecule(mid, label, name, ev(*pages), smiles=smiles)


for args in [
    ("mol_cf3ph_snbu3", "(3-CF3C6H4)SnBu3", "(3-trifluoromethylphenyl)tributyltin", (13,), "FC(F)(F)c1cccc([Sn](CCCC)(CCCC)CCCC)c1"),
    ("mol_chlorosulfonyl_isocyanate", "ClSO2NCO", "chlorosulfonyl isocyanate", (13,), "O=C=NS(=O)(=O)Cl"),
    ("mol_chlorobenzene", "PhCl", "chlorobenzene", (13,), "Clc1ccccc1"),
    ("mol_26lutidine", "2,6-lutidine", "2,6-lutidine", (15,), "Cc1cccc(C)n1"),
    ("mol_acetyl_chloride", "AcCl", "acetyl chloride", (15,), "CC(=O)Cl"),
    ("mol_acetic_acid", "AcOH", "acetic acid", (17,), "CC(=O)O"),
]:
    simple(*args)

# Separated E/Z identities reported on pages 15--17.  The isolated reaction
# products remain mixture registry entries; these identities support explicit
# selectivity targets and prevent configured isomers from being merged.
isomers = [
    ("mol_13aa", "13aa", "(E)-ethyl 4-(N-benzyl-4-methylphenylsulfonamido)-2-(ethanoyloxymethyl)but-2-enoate", "achiral", "E assigned; by analogy to NOE-assigned 13ba/13bb"),
    ("mol_13ab", "13ab", "(Z)-ethyl 4-(N-benzyl-4-methylphenylsulfonamido)-2-(ethanoyloxymethyl)but-2-enoate", "achiral", "Z assigned; by analogy to NOE-assigned 13ba/13bb"),
    ("mol_13ba", "13ba", "(E)-ethyl 5-ethanoyl-2-(ethanoyloxymethyl)-5-methyl-6-oxohept-2-enoate", "achiral", "E assigned by NOE experiments"),
    ("mol_13bb", "13bb", "(Z)-ethyl 5-ethanoyl-2-(ethanoyloxymethyl)-5-methyl-6-oxohept-2-enoate", "achiral", "Z assigned by NOE experiments"),
    ("mol_13ca", "13ca", "(E)-ethyl 5-cyano-2-(ethanoyloxymethyl)-5-methyl-6-oxo-6-phenylhex-2-enoate", "unspecified", "E assigned by analogy to NOE-assigned 13ba/13bb; configuration at the carbon bearing cyano, methyl, benzoyl and chain substituents is unreported"),
    ("mol_13cb", "13cb", "(Z)-ethyl 5-cyano-2-(ethanoyloxymethyl)-5-methyl-6-oxo-6-phenylhex-2-enoate", "unspecified", "Z assigned by analogy to NOE-assigned 13ba/13bb; configuration at the carbon bearing cyano, methyl, benzoyl and chain substituents is unreported"),
]
for mid, label, name, chiral, assignment in isomers:
    pages = (15, 16) if label.startswith("13a") or label.startswith("13b") else (15, 16, 17)
    add_molecule(mid, label, name, ev(*pages), chirality=chiral,
                 geometry="specified", representation="not_encoded",
                 assignment=assignment,
                 note="The E/Z identity is source-assigned and separately registered; exact stereochemical SMILES is withheld pending structure review.")


def charge(p, field, raw, number, unit, evidence):
    p["charges"][field] = F(raw, q(number, unit), evidence)


def participant(r, mid, role, evidence, *, stage="reaction", mass=None,
                amount=None, volume=None, equiv=None, molpct=None,
                product_role=None):
    pid = f"{r['reaction_id']}:{mid}:{len(r['participants']) + 1}"
    p = new_participant(pid, mid, role, evidence, stage=stage)
    for field, value in (("mass", mass), ("amount", amount), ("volume", volume)):
        if value:
            charge(p, field, value[0], value[1], value[2], evidence)
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
         atmosphere=None, stage="reaction", participants=()):
    s = new_step(f"{r['reaction_id']}:step{len(r['steps']) + 1}", evidence, stage=stage)
    s["participant_ids"] = list(participants)
    s["operation"] = F(operation, operation, evidence)
    if temperature is not None:
        norm = q(temperature, "°C") if isinstance(temperature, (int, float)) else temperature
        s["temperature"] = F(str(temperature), norm, evidence)
    if time is not None:
        if isinstance(time, tuple) and time[0] == "range":
            s["time"] = F(f"{time[1]}-{time[2]} h", {"minimum": time[1], "maximum": time[2], "unit": "h"}, evidence)
        elif isinstance(time, tuple):
            s["time"] = F(time[0], q(time[1], time[2]), evidence)
        else:
            s["time"] = F(str(time), time, evidence)
    if atmosphere:
        s["atmosphere"] = F(atmosphere, atmosphere, evidence)
    r["steps"].append(s)
    return s["step_id"]


def outcome(r, kind, targets, raw, normalized, evidence, *, yield_type=None,
            basis=None, scope="reaction outcome"):
    o = new_outcome(f"{r['reaction_id']}:out{len(r['outcomes']) + 1}", kind, targets, evidence)
    o["value"] = F(raw, normalized, evidence)
    o["yield_type"] = F(yield_type, yield_type, evidence) if yield_type else missing(evidence)
    o["measurement_basis"] = F(basis, basis, evidence) if basis else missing(evidence, "Measurement basis is not reported.")
    o["scope"] = F(scope, scope, evidence)
    r["outcomes"].append(o)
    return o


def reset_detail(r, evidence):
    union_evidence(r, evidence)
    r["participants"] = []
    r["steps"] = []
    r["review"]["status"] = "needs_review"
    r["review"]["method"] = "assistant_visual"
    r["review"]["reviewer"] = None
    r["review"]["checked_fields"] = ["original page image", "compound labels", "reagent quantities", "ordered operations", "yield and analytical basis"]


def replace_yield_evidence(r, target, evidence, basis="isolated product"):
    for o in r["outcomes"]:
        if o["kind"] == "yield" and target in o["target_molecule_ids"]:
            o["evidence_ids"] = list(dict.fromkeys(o.get("evidence_ids", []) + evidence))
            o["measurement_basis"] = F(basis, basis, evidence)


def set_classification(r, labels, evidence, limitation):
    changes, normalized = [], []
    for i, code in enumerate(labels, 1):
        a, b = code.split("-")[-2:]
        item = bond_label(a, b, 1)
        cid = f"b{i}"
        change = {"change_id": cid, "atoms": [{"atom_map": i * 2 - 1, "element": a}, {"atom_map": i * 2, "element": b}],
                  "before_order": 0, "after_order": 1, "change_type": "formed"}
        changes.append(F(change, change, evidence, "inference", derivation="Placeholder maps encode only the inferred net formed-bond class."))
        normalized.append({"category_number": item["category_number"], "code": item["code"], "supporting_change_ids": [cid]})
    r["classification"]["bond_changes"] = changes
    r["classification"]["labels"] = F(labels, normalized, evidence, "inference", derivation="Net formed-bond class inferred from the source arrow.")
    r["classification"]["coverage"] = F("partial", "partial", evidence, "inference", derivation="No complete atom mapping is asserted.")
    r["classification"]["limitations"] = [limitation]


def annotate_existing(short, pages, note):
    r = getr(short)
    evidence = ev(*pages)
    union_evidence(r, evidence)
    r["notes"].append(F(note, note, evidence))
    r["review"]["status"] = "needs_review"


# Page 13 carbamates: analogue details merge into the page-3 records.
for suffix in ("l", "m", "ba", "bc", "bd"):
    annotate_existing(f"rxn_carb_1{suffix}", (13,),
                      "Page 13 states this preparation follows the 1b procedure and reports the product-specific isolated yield; characterization is not a separate experiment.")
    add_occurrence(f"mol_1{suffix}", f"1{suffix}", ev(13))

# 1bb is a distinct in-situ sulfonyl-isocyanate procedure and is rewritten from
# the detailed source instead of inheriting the generic Table-2 footnote.
r = getr("rxn_carb_1bb")
reset_detail(r, ev(13))
r["location"] = F("Page 13, detailed preparation of 1bb", "Page 13, detailed preparation of 1bb", ev(13))
r["review"]["unresolved"] = ["The in-situ sulfonyl isocyanate is not isolated; no yield for its formation is assigned.",
                              "The 79% value is based on recovered 7b, whereas 22% is the isolated 1bb yield; both are retained separately."]
p_sn = participant(r, "mol_cf3ph_snbu3", "reactant", ev(13), mass=("0.587 g", 0.587, "g"), amount=("1.35 mmol", 1.35, "mmol"))
p_csi = participant(r, "mol_chlorosulfonyl_isocyanate", "reactant", ev(13), volume=("0.120 mL", 0.120, "mL"), amount=("1.35 mmol", 1.35, "mmol"))
p_cb = participant(r, "mol_chlorobenzene", "solvent", ev(13), volume=("2 mL", 2, "mL"))
p_7b = participant(r, "mol_7b", "reactant", ev(13), mass=("0.200 g", 0.200, "g"), amount=("1.28 mmol", 1.28, "mmol"))
p_dcm = participant(r, "mol_dcm", "solvent", ev(13), volume=("2 mL", 2, "mL"))
participant(r, "mol_1bb", "product", ev(13), mass=("114 mg", 114, "mg"), product_role="desired product")
participant(r, "mol_etoac", "purification_agent", ev(13), stage="purification")
step(r, "heat aryltributyltin and chlorosulfonyl isocyanate in dry chlorobenzene under reflux", ev(13), time=("15.5 h", 15.5, "h"), participants=[p_sn, p_csi, p_cb])
step(r, "cool to room temperature and then to 0 °C", ev(13), temperature=0, participants=[p_sn, p_csi, p_cb])
step(r, "add allenylic alcohol 7b dropwise", ev(13), temperature=0, participants=[p_7b])
step(r, "after 10 min add dry dichloromethane because of a solubility problem", ev(13), temperature=0, time=("10 min", 10, "min"), participants=[p_dcm])
step(r, "stir at 0 °C", ev(13), temperature=0, time=("30 min", 30, "min"), participants=[p_7b, p_dcm])
step(r, "remove the ice bath and stir at room temperature", ev(13), temperature="room temperature", time=("30 min", 30, "min"), participants=[p_7b, p_dcm])
step(r, "quench with water, extract/wash, dry over sodium sulfate and concentrate", ev(13), stage="workup")
step(r, "flash chromatography: hexanes wash, then 5-20% ethyl acetate/hexanes", ev(13), stage="purification")
replace_yield_evidence(r, "mol_1bb", ev(13), "based on recovered starting material")
for o in r["outcomes"]:
    if o["kind"] == "yield" and o["value"]["raw_value"] == "22":
        o["measurement_basis"] = F("isolated product", "isolated product", ev(13))
        o["evidence_ids"] = list(dict.fromkeys(o["evidence_ids"] + ev(13)))
outcome(r, "product_amount", ["mol_7b"], "145 mg recovered 7b", q(145, "mg"), ev(13), yield_type="recovered starting material amount", basis="isolated recovered 7b", scope="recovered starting material from the 1bb preparation")
add_occurrence("mol_1bb", "1bb", ev(13))


# Detailed annulations 3a and 3b.
r = getr("rxn_ann_3a")
reset_detail(r, ev(13))
r["review"]["unresolved"] = ["No reaction temperature is explicitly stated in the detailed 3a paragraph; it was not copied from another procedure."]
p1 = participant(r, "mol_1a", "reactant", ev(13), mass=("50.0 mg", 50.0, "mg"), amount=("0.147 mmol", 0.147, "mmol"), equiv=("1.0 equiv", 1.0))
pp = participant(r, "mol_pph3", "catalyst", ev(13), mass=("38.6 mg", 38.6, "mg"), amount=("0.147 mmol", 0.147, "mmol"), equiv=("1.0 equiv", 1.0))
ps = participant(r, "mol_mecn", "solvent", ev(13), volume=("2 mL + 2 mL", 4, "mL"))
participant(r, "mol_3a", "product", ev(13), mass=("29.9 mg", 29.9, "mg"), product_role="desired product")
step(r, "add the 1a solution to the PPh3 solution by syringe pump; fasten the syringe/needle joint to limit evaporation", ev(13), time=("24 h", 24, "h"), participants=[p1, pp, ps])
step(r, "stir after completion of the addition", ev(13), time=("4 h", 4, "h"), participants=[p1, pp, ps])
step(r, "evaporate solvent and purify by flash chromatography (10% ethyl acetate/hexanes)", ev(13), stage="purification")
replace_yield_evidence(r, "mol_3a", ev(13))
add_occurrence("mol_3a", "3a", ev(13))

r = getr("rxn_ann_3b")
reset_detail(r, ev(13, 14))
r["review"]["unresolved"] = ["The detailed paragraph first describes 1b in 5-6 mL MeCN but later says rinsing continued until the syringe contained 4 mL of solution; both source statements are retained and no volume was silently selected."]
p1 = participant(r, "mol_1b", "reactant", ev(13, 14), mass=("50.0 mg", 50.0, "mg"), amount=("0.142 mmol", 0.142, "mmol"), equiv=("1.0 equiv", 1.0))
pp = participant(r, "mol_pbu3", "catalyst", ev(14), volume=("0.035 mL", 0.035, "mL"), amount=("0.1415 mmol", 0.1415, "mmol"), equiv=("1.0 equiv", 1.0))
ps = participant(r, "mol_mecn", "solvent", ev(13, 14), volume=("5-6 mL for initial 1b solution; 4 mL stock in catalyst flask; syringe finally stated as 4 mL", 9, "mL"))
participant(r, "mol_3b", "product", ev(14), mass=("32.3 mg", 32.3, "mg"), product_role="desired product")
step(r, "bubble argon through dry MeCN stock and separately through the 1b solution", ev(13), time=("1.5 h", 1.5, "h"), atmosphere="argon", participants=[p1, ps])
step(r, "evacuate the receiving flask for 5 min and backfill with argon; repeat twice", ev(13, 14), time=("5 min per evacuation", 5, "min"), atmosphere="argon")
step(r, "add deaerated MeCN and tri-n-butylphosphine to the receiving flask", ev(14), temperature="room temperature", atmosphere="argon", participants=[pp, ps])
step(r, "transfer/rinse the carbamate solution and add it to the phosphine solution by syringe pump", ev(14), temperature="room temperature", time=("24 h", 24, "h"), atmosphere="argon", participants=[p1, pp, ps])
step(r, "stir after completion of the addition", ev(14), temperature="room temperature", time=("4 h", 4, "h"), atmosphere="argon", participants=[p1, pp, ps])
step(r, "evaporate solvent and purify by flash chromatography (10% ethyl acetate/hexanes)", ev(14), stage="purification")
replace_yield_evidence(r, "mol_3b", ev(14))
add_occurrence("mol_3b", "3b", ev(13, 14))

# The remaining products explicitly follow the 3b procedure; product-specific
# doses are not printed and remain unreported.  3k retains its no-deaeration exception.
page14 = ["c", "d", "e", "f", "g", "h", "i", "j", "k", "l"]
page15 = ["m", "ba", "bb", "bc", "bd"]
for suffix, pages in [(s, (14,)) for s in page14] + [(s, (15,)) for s in page15]:
    annotate_existing(f"rxn_ann_3{suffix}", pages,
                      "The product paragraph explicitly inherits the detailed 3b procedure; product-specific starting amounts are not reported. Characterization is not a separate experiment." + (" This 3k experiment explicitly omits solvent/carbamate deaeration." if suffix == "k" else ""))
    add_occurrence(f"mol_3{suffix}", f"3{suffix}", ev(*pages))


# New current-paper experiment: preparation of allenylic acetate 9.
new_id = rid("rxn_acetate_9")
if new_id in reactions:
    raise ValueError("Refusing to duplicate the page-15 preparation of 9")
r = new_reaction(new_id, DOC, ev(15))
r["record_type"] = "specific_experiment"
r["experiment_categories"] = F(["substrate_preparation", "mechanistic_study"], ["other", "mechanistic_study"], ev(15))
r["location"] = F("Page 15, Synthesis Performed for the Mechanistic Investigation, Allenylic Acetate 9", "Page 15, preparation of 9", ev(15))
r["reaction_name"] = F("Preparation of allenylic acetate 9", "Preparation of allenylic acetate 9", ev(15))
r["transformation_type"] = F("acetylation of allenylic alcohol 7a", "acetylation of allenylic alcohol 7a", ev(15))
r["review"] = {"status": "needs_review", "method": "assistant_visual", "reviewer": None,
                 "checked_fields": ["original page image", "compound labels", "all reagent quantities", "ordered temperature/time stages", "workup", "isolated yield"],
                 "unresolved": ["Complete atom mapping is not source-verified; the net C-O classification is inferred from the named product and source structure."]}
set_classification(r, ["C-O"], ev(15), "No complete atom-mapped reaction SMILES; net ester bond formation is inferred from the source.")
p7 = participant(r, "mol_7a", "reactant", ev(15), mass=("0.200 g", 0.200, "g"), amount=("1.41 mmol", 1.41, "mmol"))
pl = participant(r, "mol_26lutidine", "reagent", ev(15), volume=("265 μL", 265, "μL"), amount=("2.28 mmol", 2.28, "mmol"))
pa = participant(r, "mol_acetyl_chloride", "reactant", ev(15), volume=("225 μL", 225, "μL"), amount=("3.10 mmol", 3.10, "mmol"))
pd = participant(r, "mol_dcm", "solvent", ev(15), volume=("10 mL", 10, "mL"))
participant(r, "mol_hcl", "workup_agent", ev(15), stage="workup", volume=("10 mL of 2 M HCl", 10, "mL"))
participant(r, "mol_et2o", "workup_agent", ev(15), stage="workup", volume=("3 × 5 mL", 15, "mL"))
participant(r, "mol_nacl", "workup_agent", ev(15), stage="workup", volume=("10 mL brine for each organic set", 10, "mL"))
participant(r, "mol_na2so4", "workup_agent", ev(15), stage="workup")
participant(r, "mol_etoac", "purification_agent", ev(15), stage="purification")
participant(r, "mol_9", "product", ev(15), mass=("219 mg", 219, "mg"), product_role="desired product")
step(r, "add 2,6-lutidine dropwise to 7a in dry dichloromethane", ev(15), temperature=0, participants=[p7, pl, pd])
step(r, "stir before adding acetyl chloride", ev(15), temperature=0, time=("5 min", 5, "min"), participants=[p7, pl, pd])
step(r, "add acetyl chloride dropwise", ev(15), temperature=0, participants=[pa])
step(r, "stir at 0 °C", ev(15), temperature=0, time=("1 h", 1, "h"), participants=[p7, pl, pa, pd])
step(r, "remove the ice bath, warm to room temperature and stir", ev(15), temperature="room temperature", time=("17 h", 17, "h"), participants=[p7, pl, pa, pd])
step(r, "quench with 2 M HCl, separate, extract with dichloromethane then ether, wash each organic set with brine, combine, dry and concentrate", ev(15), stage="workup")
step(r, "flash chromatography (0-10% ethyl acetate/hexanes)", ev(15), stage="purification")
outcome(r, "yield", ["mol_9"], "85%", q(85, "%"), ev(15), yield_type="isolated yield", basis="isolated product 9", scope="preparation of allenylic acetate 9")
bundle["reactions"].append(r)
reactions[new_id] = r
add_occurrence("mol_9", "9", ev(15))


# Tables 4--6 receive the detailed general procedures and distinct E/Z targets.
for table, count, pages in [(4, 6, (15, 16)), (5, 7, (16,)), (6, 4, (16,))]:
    for entry in range(1, count + 1):
        short = f"rxn_mech_t{table}_e{entry:02d}"
        r = getr(short)
        union_evidence(r, ev(*pages))
        r["notes"].append(F(
            "Detailed experimental-section general procedure: 9 (0.190 mmol, 1.0 equiv), pronucleophile (1.2 equiv), and the row-specific PPh3/Cs2CO3/NaOAc charges; Table 4 uses dry MeCN with 6 h addition and 1-1.5 h post-stir, while Tables 5/6 use dry benzene with 1 h addition and the row-specific post-stir time. Product fractions were isolated by flash chromatography.",
            None, ev(*pages)))
        r["review"]["status"] = "needs_review"
        if table == 5:
            r["review"]["unresolved"].append("The detailed experimental procedure on page 16 and Table 5 specify dry benzene; earlier page-4 prose was recorded as a conflicting CH3CN statement. The conflict is preserved, with the detailed procedure supporting benzene.")

# Update selectivity targets and assignment notes without pretending the ratio
# measurement method is reported.
target_map = {
    "a": (["mol_13ab", "mol_13aa"], "Z (13ab):E (13aa); individual geometry assigned by analogy to NOE-assigned 13b isomers"),
    "b": (["mol_13ba", "mol_13bb"], "E (13ba):Z (13bb); individual geometry assigned by NOE experiments"),
    "c": (["mol_13ca", "mol_13cb"], "E (13ca):Z (13cb); individual geometry assigned by analogy to NOE-assigned 13b isomers"),
}
for table, letter, entries in [(4, "a", range(1, 7)), (5, "b", range(1, 8)), (6, "c", range(1, 5))]:
    for entry in entries:
        r = getr(f"rxn_mech_t{table}_e{entry:02d}")
        for o in r["outcomes"]:
            if o["kind"] == "selectivity":
                o["target_molecule_ids"] = target_map[letter][0]
                o["scope"] = F(target_map[letter][1], target_map[letter][1], ev(15, 16))
                o["measurement_basis"] = missing(ev(15, 16), "The pages report isolated-mixture E/Z ratios but do not state the analytical method used to measure the ratios; NOE/analogy concerns assignment of isomer identity, not ratio quantitation.")
                o["evidence_ids"] = list(dict.fromkeys(o.get("evidence_ids", []) + ev(15, 16)))
                for target in target_map[letter][0]:
                    if not any(p["molecule_id"] == target for p in r["participants"]):
                        participant(r, target, "product", ev(15, 16), product_role="component identity of isolated E/Z mixture")


# Page 17 detailed mechanistic experiment (Scheme 7).
r = getr("rxn_scheme7_9_19")
reset_detail(r, ev(17))
r["review"]["unresolved"] = ["The source attributes 30/31 to the SN2 pathway and 32 to SN2-prime; these pathway assignments remain author interpretation, not directly observed elementary steps."]
p9 = participant(r, "mol_9", "reactant", ev(17), mass=("35.0 mg", 35.0, "mg"), amount=("0.190 mmol", 0.190, "mmol"), equiv=("1.0 equiv", 1.0))
pcs = participant(r, "mol_cs2co3", "reagent", ev(17), mass=("80.5 mg", 80.5, "mg"), amount=("0.247 mmol", 0.247, "mmol"))
p19 = participant(r, "mol_19", "reactant", ev(17), mass=("33.1 mg", 33.1, "mg"), amount=("0.228 mmol", 0.228, "mmol"), equiv=("1.2 equiv", 1.2))
pb = participant(r, "mol_benzene", "solvent", ev(17), volume=("2.8 mL + 2.8 mL", 5.6, "mL"))
participant(r, "mol_nh4cl", "workup_agent", ev(17), stage="workup", volume=("5 mL saturated aqueous NH4Cl", 5, "mL"))
participant(r, "mol_et2o", "workup_agent", ev(17), stage="workup", volume=("3 × 5 mL", 15, "mL"))
participant(r, "mol_nacl", "workup_agent", ev(17), stage="workup", volume=("10 mL brine", 10, "mL"))
participant(r, "mol_na2so4", "workup_agent", ev(17), stage="workup")
participant(r, "mol_30", "product", ev(17), mass=("9.8 mg", 9.8, "mg"), product_role="isolated product")
participant(r, "mol_31", "product", ev(17), mass=("13.8 mg", 13.8, "mg"), product_role="isolated product")
participant(r, "mol_32", "product", ev(17), mass=("5.9 mg", 5.9, "mg"), product_role="isolated product")
step(r, "add the solution of 9 to the Cs2CO3/benzoyl acetonitrile suspension by syringe pump", ev(17), temperature="room temperature", time=("1 h", 1, "h"), participants=[p9, pcs, p19, pb])
step(r, "stir after completion of the addition", ev(17), temperature="room temperature", time=("19 h", 19, "h"), participants=[p9, pcs, p19, pb])
step(r, "quench with saturated aqueous NH4Cl, extract with ether, wash with brine, dry and concentrate", ev(17), stage="workup")
step(r, "flash chromatography (0-15% ethyl acetate/hexanes) gives 30 and co-eluted 31/32; preparative TLC separates 31 and 32", ev(17), stage="purification")
for target in ("mol_30", "mol_31", "mol_32"):
    replace_yield_evidence(r, target, ev(17))


# Scheme 9 no-reaction controls: now the analytical basis is explicit.
for short, loading in [("rxn_scheme9_39_pph3_catalytic", "0.2 equiv PPh3"),
                       ("rxn_scheme9_39_pph3_stoichiometric", "1.0 equiv PPh3")]:
    r = getr(short)
    union_evidence(r, ev(17))
    for o in r["outcomes"]:
        if o["kind"] == "no_reaction":
            o["evidence_ids"] = list(dict.fromkeys(o.get("evidence_ids", []) + ev(17)))
            o["measurement_basis"] = F(
                "TLC and 1H NMR monitoring; crude 1H NMR showed only 39/PPh3 and 31P NMR showed only PPh3 after 24 h",
                ["TLC", "1H NMR", "31P NMR"], ev(17))
    r["notes"].append(F("The detailed page-17 control quantities and the catalytic/stoichiometric PPh3 distinction are structured in participants; TLC, 1H NMR and 31P NMR support the author-reported no-reaction outcome.", None, ev(17)))
    r["review"]["unresolved"] = ["No instrumental detection limit is reported; the record preserves the authors' no-reaction conclusion and stated monitoring methods."]


# 3l -> 40: preserve the 4.5 h scheme label and 4.25 h detailed procedure as
# explicit alternatives rather than silently choosing one.
r = getr("rxn_scheme10_3l_40")
reset_detail(r, ev(17))
r["review"]["unresolved"] = ["Time conflict: Scheme 10 on page 8 says 4.5 h, while the detailed page-17 procedure says 4.25 h; neither was overwritten.",
                              "Compound 40 has relative wedge/dash depiction without a reported absolute assignment; structure remains partially specified/not encoded."]
p3 = participant(r, "mol_3l", "reactant", ev(17), mass=("25.0 mg", 25.0, "mg"), amount=("0.0630 mmol", 0.0630, "mmol"))
pmg = participant(r, "mol_mg", "reagent", ev(17), mass=("38.4 mg", 38.4, "mg"), amount=("1.58 mmol", 1.58, "mmol"), equiv=("25 equiv", 25))
pme = participant(r, "mol_meoh", "solvent", ev(17), volume=("3.2 mL", 3.2, "mL"))
participant(r, "mol_40", "product", ev(17), mass=("10.0 mg", 10.0, "mg"), product_role="desired product")
step_id = step(r, "sonicate 3l with freshly crushed magnesium in dry methanol under an argon balloon", ev(17), atmosphere="argon", time=("4.25 h (detailed procedure); Scheme 10: 4.5 h", 4.25, "h"), participants=[p3, pmg, pme])
page8_evidence = [f"{DOC}:page8:image", f"{DOC}:page8:text"]
r["steps"][-1]["time"] = fact(status="conflict", evidence=page8_evidence + ev(17), derivation="Scheme arrow and detailed procedure report different sonication times.", alternatives=[F("4.5 h (Scheme 10)", q(4.5, "h"), page8_evidence), F("4.25 h (detailed procedure)", q(4.25, "h"), ev(17))])
step(r, "dilute with dichloromethane/saturated NH4Cl, perform the stated washes/extractions, dry and concentrate", ev(17), stage="workup")
step(r, "flash chromatography with the stated hexanes/ethyl acetate gradient; visualize with iodine", ev(17), stage="purification")
replace_yield_evidence(r, "mol_40", ev(17))
add_occurrence("mol_40", "40", ev(17))


# 3m -> S12 detailed deprotection.
r = getr("rxn_scheme11_3m_s12")
reset_detail(r, ev(17))
r["review"]["unresolved"] = ["No configuration is assigned to the stereogenic pyrroline center; unspecified is not treated as racemic."]
p3 = participant(r, "mol_3m", "reactant", ev(17), mass=("53 mg", 53, "mg"), amount=("0.13 mmol", 0.13, "mmol"))
pt = participant(r, "mol_tfa", "reagent", ev(17), volume=("0.63 mL", 0.63, "mL"))
pd = participant(r, "mol_dcm", "solvent", ev(17), volume=("3.6 mL", 3.6, "mL"))
participant(r, "mol_acetic_acid", "purification_agent", ev(17), stage="purification")
participant(r, "mol_s12", "product", ev(17), mass=("43.7 mg", 43.7, "mg"), product_role="desired product")
step(r, "add TFA to 3m in dry dichloromethane and stir", ev(17), temperature="room temperature", time=("3 h", 3, "h"), participants=[p3, pt, pd])
step(r, "evaporate solvent and purify by flash chromatography (10-30% ethyl acetate/hexanes, then 30% ethyl acetate/hexanes with 1% acetic acid)", ev(17), stage="purification")
replace_yield_evidence(r, "mol_s12", ev(17))
add_occurrence("mol_s12", "S12", ev(17))


# S12 -> 41: two operations and two separate NaBH4 charges; the 73% remains
# scoped over the complete sequence.
r = getr("rxn_scheme11_s12_41")
reset_detail(r, ev(17, 18))
r["review"]["unresolved"] = ["The 73% isolated yield is explicitly over the two operations from S12 to 41 and is not a yield for the second NaBH4 addition alone."]
ps = participant(r, "mol_s12", "reactant", ev(17), mass=("21 mg", 21, "mg"), amount=("0.057 mmol", 0.057, "mmol"))
pe = participant(r, "mol_et3n", "reagent", ev(17), volume=("16 μL", 16, "μL"), amount=("0.11 mmol", 0.11, "mmol"))
pc = participant(r, "mol_ethyl_chloroformate", "reagent", ev(17), volume=("8.0 μL of 97% material", 8.0, "μL"), amount=("0.081 mmol", 0.081, "mmol"))
pthf = participant(r, "mol_thf", "solvent", ev(17), volume=("500 μL", 500, "μL"))
pme = participant(r, "mol_meoh", "solvent", ev(17), volume=("1 mL", 1, "mL"))
pbh1 = participant(r, "mol_nabh4", "reagent", ev(17), mass=("6 mg first portion", 6, "mg"), amount=("0.14 mmol first portion", 0.14, "mmol"))
pbh2 = participant(r, "mol_nabh4", "reagent", ev(17), mass=("6 mg second portion", 6, "mg"), amount=("0.14 mmol second portion", 0.14, "mmol"))
participant(r, "mol_41", "product", ev(17, 18), mass=("14.7 mg", 14.7, "mg"), product_role="desired product")
step(r, "cool S12 in dry THF and add triethylamine and ethyl chloroformate sequentially dropwise", ev(17), temperature=0, participants=[ps, pe, pc, pthf])
step(r, "remove the bath and stir to form the mixed anhydride", ev(17), temperature="room temperature", time=("1.2 h", 1.2, "h"), participants=[ps, pe, pc, pthf])
step(r, "filter through Celite, wash with ethyl acetate, concentrate and azeotropically remove water with benzene", ev(17), stage="workup")
step(r, "dissolve the crude material in dry methanol, cool, add the first sodium borohydride portion and stir", ev(17), temperature=-78, time=("30 min", 30, "min"), participants=[pme, pbh1])
step(r, "add the second equal sodium borohydride portion and continue stirring", ev(17), temperature=-78, time=("1 h", 1, "h"), participants=[pme, pbh2])
step(r, "quench with saturated NH4Cl, warm slowly to room temperature, extract/wash/dry and concentrate", ev(17), temperature="room temperature", stage="workup")
step(r, "flash chromatography (0-60% ethyl acetate/hexanes)", ev(17), stage="purification")
replace_yield_evidence(r, "mol_41", ev(17, 18), "isolated endpoint 41 from starting acid S12 over two operations")
add_occurrence("mol_41", "41", ev(17, 18))


# 41 -> 42: preserve two separate magnesium additions and sonication periods.
r = getr("rxn_scheme11_41_42")
reset_detail(r, ev(18))
r["review"]["unresolved"] = ["The source reports one isolated yield for the combined detosylation/in situ lactamization; no stage-specific yields are available."]
p41 = participant(r, "mol_41", "reactant", ev(18), mass=("14 mg", 14, "mg"), amount=("0.040 mmol", 0.040, "mmol"))
pm1 = participant(r, "mol_mg", "reagent", ev(18), mass=("25 mg first portion", 25, "mg"), amount=("1.03 mmol first portion", 1.03, "mmol"))
pm2 = participant(r, "mol_mg", "reagent", ev(18), mass=("25 mg second portion", 25, "mg"), amount=("1.03 mmol second portion", 1.03, "mmol"))
pme = participant(r, "mol_meoh", "solvent", ev(18), volume=("2 mL", 2, "mL"))
participant(r, "mol_42", "product", ev(18), mass=("4.4 mg", 4.4, "mg"), product_role="desired product")
step(r, "sonicate 41 with the first portion of freshly crushed magnesium in dry methanol under an argon balloon", ev(18), atmosphere="argon", time=("1.5 h", 1.5, "h"), participants=[p41, pm1, pme])
step(r, "add the second portion of magnesium and sonicate until the magnesium dissolves", ev(18), atmosphere="argon", time=("2 h", 2, "h"), participants=[p41, pm2, pme])
step(r, "dilute with dichloromethane/saturated NH4Cl, perform the stated washes/extractions, dry and concentrate", ev(18), stage="workup")
step(r, "flash chromatography in a Pasteur pipette with the stated hexanes/ethyl acetate/methanol gradient; visualize with iodine or permanganate", ev(18), stage="purification")
replace_yield_evidence(r, "mol_42", ev(18))
add_occurrence("mol_42", "42", ev(18))


# Source context: characterization and the unlinked SI availability statement.
bundle["source_context"].extend([
    {"context_id": "pages13_18_characterization_inventory", "document_id": DOC,
     "evidence_ids": ev(13, 14, 15, 16, 17, 18), "kind": "characterization_inventory",
     "note": "IR/NMR/MS characterization paragraphs and individual spectra descriptions support product identities but are not independent synthesis experiments or mass-spectrum peak-table reactions."},
    {"context_id": "page18_si_availability_unlinked", "document_id": DOC,
     "evidence_ids": ev(18), "kind": "supplementary_information_statement",
     "note": "Page 18 states that Supporting Information is available online, but no SI file was explicitly supplied in this workspace; no SI document association was inferred."},
    {"context_id": "page18_references_context", "document_id": DOC,
     "evidence_ids": ev(18), "kind": "bibliographic_context",
     "note": "Acknowledgements and references are not current-paper experiments and were not converted into reaction records."},
])

# Exact page coverage required for this batch.
page_records = {
    13: ["rxn_carb_1l", "rxn_carb_1m", "rxn_carb_1ba", "rxn_carb_1bb", "rxn_carb_1bc", "rxn_carb_1bd", "rxn_ann_3a", "rxn_ann_3b"],
    14: ["rxn_ann_3b"] + [f"rxn_ann_3{x}" for x in page14],
    15: [f"rxn_ann_3{x}" for x in page15] + ["rxn_acetate_9"] + [f"rxn_mech_t4_e{x:02d}" for x in range(1, 7)],
    16: [f"rxn_mech_t4_e{x:02d}" for x in range(1, 7)] + [f"rxn_mech_t5_e{x:02d}" for x in range(1, 8)] + [f"rxn_mech_t6_e{x:02d}" for x in range(1, 5)],
    17: ["rxn_scheme7_9_19", "rxn_scheme9_39_pph3_catalytic", "rxn_scheme9_39_pph3_stoichiometric", "rxn_scheme10_3l_40", "rxn_scheme11_3m_s12", "rxn_scheme11_s12_41"],
    18: ["rxn_scheme11_s12_41", "rxn_scheme11_41_42"],
}
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 13, "state": "attempted", "regions": [
            {"region": "upper-left continuation and carbamate paragraphs 1l/1m/1ba/1bc/1bd", "finding": "These analogue preparations explicitly inherit the 1b procedure and report isolated yields; they update existing scope records rather than create duplicate experiments. Spectral text is characterization only.", "evidence_ids": ev(13)},
            {"region": "middle-right detailed 1bb procedure", "finding": "The in-situ sulfonyl-isocyanate sequence, 22% isolated 1bb yield, 79% yield based on recovered starting material and 145 mg recovered 7b were transcribed with separate outcomes.", "evidence_ids": ev(13)},
            {"region": "lower section, detailed 3a and start of 3b", "finding": "The 3a syringe-pump annulation was fully detailed; the 3b deaeration setup starts here and continues on page 14. These merge into existing Table-1/Table-3 experiments.", "evidence_ids": ev(13)}],
         "record_ids": [rid(x) for x in page_records[13]],
         "unresolved": ["Analogue-specific starting doses for 1l/1m/1ba/1bc/1bd are not printed.", "No reaction temperature is stated in the detailed 3a paragraph.", "The 1bb in-situ sulfonyl isocyanate is not isolated and receives no independent yield."]},
        {"page": 14, "state": "attempted", "regions": [
            {"region": "upper-left continuation of 3b", "finding": "The PBu3 charge, 24 h room-temperature addition, 4 h hold, purification, 32.3 mg product and 74% isolated yield complete the existing 3b experiment.", "evidence_ids": ev(14)},
            {"region": "3c-3l analogue product paragraphs", "finding": "The source says these products follow the 3b procedure and provides product-specific yields/characterization; existing scope entries were updated without creating duplicates. 3k explicitly omits deaeration.", "evidence_ids": ev(14)}],
         "record_ids": [rid(x) for x in page_records[14]],
         "unresolved": ["The 3b paragraph describes an initial 5-6 mL carbamate solution but later a final 4 mL syringe solution; both statements are preserved as a source inconsistency.", "Analogue-specific starting doses are not reported for 3c-3l."]},
        {"page": 15, "state": "attempted", "regions": [
            {"region": "top continuation and 3m/3ba-3bd paragraphs", "finding": "Product-specific isolated yields and characterization update the existing annulation scope records; no duplicate experiments were made.", "evidence_ids": ev(15)},
            {"region": "middle, Allenylic Acetate 9", "finding": "A new current-paper preparation of 9 from 7a, 2,6-lutidine and acetyl chloride was recorded with staged 0 °C/room-temperature operations, workup, purification and 85% isolated yield.", "evidence_ids": ev(15)},
            {"region": "lower, Table-4 general procedure and E/Z discussion", "finding": "The detailed Table-4 procedure was attached to all six existing rows. The isolated 13a-c mixtures remain reaction products, while E/Z isomers 13aa/ab, 13ba/bb and 13ca/cb were registered separately; NOE assignment for 13b and analogy for 13a/c were retained.", "evidence_ids": ev(15)}],
         "record_ids": [rid(x) for x in page_records[15]],
         "unresolved": ["The analytical method used to quantify the reported E/Z ratios is not stated; NOE/analogy supplies identity assignment, not ratio measurement.", "No SI spectra file was supplied, so the referenced NOE data were not independently inspected."]},
        {"page": 16, "state": "attempted", "regions": [
            {"region": "upper characterization of 12a/13aa/13ab", "finding": "Characterization confirms separate named E and Z identities but is not a new synthesis experiment.", "evidence_ids": ev(16)},
            {"region": "middle, General Procedure for 9 with 10b/10c", "finding": "The dry-benzene general procedure, optional PPh3/Cs2CO3/NaOAc charges, addition times, row-specific post-stir times and chromatography fractions were linked to the existing Table-5 and Table-6 records.", "evidence_ids": ev(16)},
            {"region": "lower characterization of 11b/12b/13ba/13bb/11c/12c/13ca", "finding": "These are analytical characterization blocks; E/Z identities remain separate molecule records and are not duplicate reactions.", "evidence_ids": ev(16)}],
         "record_ids": [rid(x) for x in page_records[16]],
         "unresolved": ["Earlier page-4 prose was recorded as conflicting with Table 5 on solvent; page 16 explicitly supports dry benzene, but the conflict was not silently erased.", "Ratio quantitation method remains unreported."]},
        {"page": 17, "state": "attempted", "regions": [
            {"region": "top continuation of 13cb characterization and Scheme-7 experiment", "finding": "Characterization was kept as context. The detailed 9/19/Cs2CO3 reaction was merged into the existing Scheme-7 experiment with isolated masses/yields for 30-32 and separate flash/preparative-TLC operations.", "evidence_ids": ev(17)},
            {"region": "middle, reaction of 39 with PPh3", "finding": "Both catalytic and stoichiometric PPh3 no-reaction controls remain independent; TLC, crude 1H NMR and 31P NMR are now recorded as the no-reaction measurement basis.", "evidence_ids": ev(17)},
            {"region": "lower, preparations of 40, S12 and start of 41", "finding": "Detailed application-synthesis procedures update the existing Scheme-10/11 arrows. The 4.25 h detailed time versus 4.5 h scheme time is preserved as a conflict, and the 73% S12-to-41 yield remains scoped over two operations.", "evidence_ids": ev(17)}],
         "record_ids": [rid(x) for x in page_records[17]],
         "unresolved": ["No detection limit is reported for the Scheme-9 no-reaction controls.", "The Scheme-10/detailed-procedure time conflict is unresolved.", "Absolute configuration of 40 and the Scheme-11 intermediates is not assigned."]},
        {"page": 18, "state": "attempted", "regions": [
            {"region": "upper continuation of 41 characterization", "finding": "Spectral characterization supports 41 but is not an independent reaction.", "evidence_ids": ev(18)},
            {"region": "middle, conversion of 41 to 42", "finding": "Two separate magnesium additions and 1.5 h/2 h sonication periods, full workup/purification and 72% isolated yield were merged into the existing 41-to-42 record.", "evidence_ids": ev(18)},
            {"region": "Supporting Information, acknowledgements and references", "finding": "The paper states that SI is available online, but no SI file was supplied and no association was inferred. Acknowledgements/references are context, not current experiments.", "evidence_ids": ev(18)}],
         "record_ids": [rid(x) for x in page_records[18]],
         "unresolved": ["No stage-specific yields are reported for detosylation versus in-situ lactamization.", "The online SI statement does not establish a supplied SI document in this workspace."]},
    ],
}

# Mark only this batch's tasks attempted.
tasks_path = ROOT / "runs/current/tasks.json"
tasks = json.loads(tasks_path.read_text())
for task in tasks:
    if task.get("document_id") == DOC and task.get("evidence_id") in {x for p in range(13, 19) for x in ev(p)}:
        task["state"] = "attempted"
        task["reason"] = "Original page image inspected; regions and records mapped in pages_0013_0018.json. Machine records remain needs_review."

# Outcome evidence belongs on each evidence-bearing fact in the v3.1 schema;
# outcomes themselves intentionally have no top-level evidence_ids property.
for reaction_record in bundle["reactions"]:
    for outcome_record in reaction_record["outcomes"]:
        outcome_record.pop("evidence_ids", None)

bundle_path.write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
tasks_path.write_text(json.dumps(tasks, ensure_ascii=False, indent=2) + "\n")
(ROOT / "runs/current/coverage/pages_0013_0018.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")

# Standard stereo diagnostics for all newly registered pages-13--18 identities.
diagnostics = []
for mid in [x[0] for x in isomers] + ["mol_cf3ph_snbu3", "mol_chlorosulfonyl_isocyanate", "mol_chlorobenzene", "mol_26lutidine", "mol_acetyl_chloride", "mol_acetic_acid"]:
    m = molecules[mid]
    smi = m["smiles"]["normalized_value"]
    diagnostics.append({"molecule_id": mid, "source_verified": False,
                        "diagnostic": stereo_inspect(smi) if smi else None,
                        "recorded_status": {k: m["stereochemistry"][k]["normalized_value"] for k in ("chirality_status", "geometry_status", "representation_status")}})
(ROOT / "runs/current/raw_model/stereo_checks_pages_0013_0018.json").write_text(json.dumps(diagnostics, ensure_ascii=False, indent=2) + "\n")

print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
                  "new_reactions": 1, "new_isomer_identities": 6}, ensure_ascii=False))
