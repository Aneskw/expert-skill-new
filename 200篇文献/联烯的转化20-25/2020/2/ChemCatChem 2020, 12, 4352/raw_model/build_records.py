#!/usr/bin/env python3
"""Build the source-derived v3.1 submission for PDF pages 1--6.

This is an evidence-preserving transcription recipe.  It deliberately leaves
complex drawn structures unresolved and never promotes schema validity to a
claim of chemical correctness.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SKILL = ROOT / "skills" / "chemical-reaction-extraction"
sys.path.insert(0, str(SKILL / "scripts"))
from standard_records import (  # noqa: E402
    fact,
    new_molecule,
    new_outcome,
    new_participant,
    new_reaction,
    new_step,
    source_mentions,
)

DOC = "doc_2401da5b5fbe4542"
E = {p: f"{DOC}:page{p}:text" for p in range(1, 7)}


def F(raw=None, normalized=None, evidence=(), method="original_extraction", *,
      status=None, verification=None, derivation=None, alternatives=()):
    if normalized is None and raw is not None:
        normalized = raw
    return fact(raw, normalized, evidence, method, status, verification,
                derivation, alternatives)


def missing(evidence=(), reason=None):
    return fact(evidence=evidence, derivation=reason)


def qty(value, unit, **extra):
    out = {"number": value, "unit": unit}
    out.update(extra)
    return out


bundle = json.loads((ROOT / "runs/current/submitted.json").read_text())
dataset_path = ROOT / "runs/current/dataset.json"
dataset = json.loads(dataset_path.read_text())
document = dataset["documents"][0]
document["title"] = (
    "Phosphine-Catalyzed (4 + 1) Annulation: Rearrangement of Allenylic "
    "Carbamates to 3-Pyrrolines through Phosphonium Diene Intermediates"
)
document["source_name"] = "ChemCatChem 2020, 12, 4352.pdf"
document["supplementary_material"] = {
    "status": "unestablished",
    "note": "The main text links supporting information, but no SI file was explicitly supplied in this workspace; no association was inferred.",
}
dataset_path.write_text(json.dumps(dataset, ensure_ascii=False, indent=2) + "\n")
bundle["documents"] = dataset["documents"]

molecules = {}


def add_molecule(mid, label, evidence, *, description=None, smiles=None,
                 chirality="unresolved", geometry="unresolved",
                 representation="partial", assignment=None, mixture=None,
                 axial=None, identity="source_identified"):
    ev = list(dict.fromkeys(evidence))
    m = new_molecule(mid, DOC, ev)
    m["identity_status"] = identity
    m["occurrences"][0]["compound_label"] = F(label, label, ev)
    m["names"]["original"] = F(label, label, ev)
    m["names"]["english"] = F(description or label, description or label, ev)
    m["names"]["abbreviations"] = F([label], [label], ev)
    m["names"]["aliases"] = missing(ev, "No additional alias is reported on pages 1--6.")
    m["cas_numbers"] = missing(ev, "CAS number is not reported on pages 1--6 and was not guessed.")
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = F(smiles, smiles, ev, "structure_derivation",
                          derivation="Conservative transcription of the source drawing/name; parsing is not identity verification.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing(ev, "Exact SMILES withheld where the page drawing or stereochemical state was not safely encoded.")
    m["generic_template"] = None
    st = m["stereochemistry"]
    ctext = chirality
    if "no stereogenic center" in chirality:
        chirality = "achiral"
    elif "not applicable" in chirality:
        chirality = "not_applicable"
    elif any(x in chirality for x in ("unspecified", "possible", "no source-specified", "no assigned", "prochiral")):
        chirality = "unspecified"
    else:
        chirality = "unresolved"
    gtext = geometry
    if "mixture" in geometry:
        geometry = "mixture"
    elif any(x in geometry for x in ("not applicable", "endocyclic alkene", "allene geometry", "terminal allene")):
        geometry = "not_applicable"
    elif "unspecified" in geometry or "retained only in drawing" in geometry or "not source-specified" in geometry:
        geometry = "unspecified"
    else:
        geometry = "unresolved"
    rtext = representation
    if smiles and axial:
        representation = "partial"
    elif smiles and chirality in ("achiral", "not_applicable") and geometry == "not_applicable":
        representation = "encoded"
    elif smiles:
        representation = "partial"
    elif "not applicable" in representation:
        representation = "not_applicable"
    elif "exact" in representation or "mixture retained" in representation or "generic template" in representation:
        representation = "not_encoded"
    else:
        representation = "unresolved"
    st["chirality_status"] = F(ctext, chirality, ev, "structure_derivation",
                                 derivation="Assessed from the displayed 2D structure and explicit source labels; machine record remains needs_review.")
    st["geometry_status"] = F(gtext, geometry, ev, "structure_derivation",
                                derivation="Assessed from the displayed 2D structure and explicit source labels; absence of a label was not converted to racemic or achiral.")
    st["representation_status"] = F(rtext, representation, ev, "inference",
                                      derivation="Compares encoded representation with the source depiction.")
    st["assignment"] = F(assignment, assignment, ev) if assignment else missing(ev, "No absolute/relative stereochemical assignment reported for this identity on pages 1--6.")
    st["mixture"] = F(mixture, mixture, ev) if mixture else missing(ev, "No stereoisomer mixture statement reported for this identity on pages 1--6.")
    st["axial"] = F(axial, axial, ev, "structure_derivation",
                      derivation="Terminal allene substitution was inspected in the source drawing.") if axial else missing(ev)
    st["atom_centers"] = F(ctext, ctext, ev, "structure_derivation",
                             derivation="Source drawing inspection; no configuration was invented.")
    st["double_bonds"] = F(gtext, gtext, ev, "structure_derivation",
                             derivation="Source drawing inspection; no E/Z assignment was invented.")
    st["isomeric_smiles"] = F(smiles, smiles, ev, "structure_derivation",
                                derivation="No stereochemical token is claimed beyond the source.") if smiles else missing(ev)
    st["identity_key"] = missing(ev, "No structure-derived identity key generated for unresolved or partially encoded source identity.")
    m["identity_note"] = F(description or label, description or label, ev)
    molecules[mid] = m
    return mid


def add_generic(mid, label, evidence, fixed, variables, mapping):
    ev = list(dict.fromkeys(evidence))
    add_molecule(mid, label, ev, description=f"Generic source template {label}",
                 chirality="variable or unspecified", geometry="variable or unspecified",
                 representation="generic template; not a concrete molecule")
    m = molecules[mid]
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": F(fixed, fixed, ev),
        "variables": [{
            "label": v,
            "attachment_points": F("source drawing attachment", "source drawing attachment", ev),
            "definition": F(desc, desc, ev),
            "allowed_values": F(vals, vals, ev),
            "constraints": F("Only values explicitly shown on pages 1--6", "Only values explicitly shown on pages 1--6", ev),
        } for v, desc, vals in variables],
        "instance_mapping": [],
        "reason": F("The source uses R groups and a numbered series.", "Preserve the source template without inventing a single concrete SMILES.", ev),
    }
    return mid


# Common named reagents/solvents/catalysts used in the six-page range.
simple = {
    "mol_pph3": ("PPh3", "triphenylphosphine", [E[2], E[4], E[5], E[6]], "P(c1ccccc1)(c1ccccc1)c1ccccc1"),
    "mol_pmeph2": ("PMePh2", "methyldiphenylphosphine", [E[2]], None),
    "mol_pme3": ("PMe3", "trimethylphosphine", [E[2]], "CP(C)C"),
    "mol_pet3": ("PEt3", "triethylphosphine", [E[2]], None),
    "mol_pbu3": ("PBu3", "tri-n-butylphosphine", [E[2], E[3]], None),
    "mol_pcy3": ("PCy3", "tricyclohexylphosphine", [E[2]], None),
    "mol_dcm": ("CH2Cl2", "dichloromethane", [E[2], E[3]], "ClCCl"),
    "mol_mecn": ("CH3CN", "acetonitrile", [E[2], E[3], E[4]], "CC#N"),
    "mol_benzene": ("PhH", "benzene", [E[4], E[5]], "c1ccccc1"),
    "mol_cs2co3": ("Cs2CO3", "cesium carbonate", [E[4], E[5]], None),
    "mol_naoac": ("NaOAc", "sodium acetate", [E[4], E[6]], None),
    "mol_tsnco": ("TsNCO", "p-toluenesulfonyl isocyanate", [E[3]], None),
    "mol_pnosnco": ("p-NO2C6H4SO2NCO", "p-nitrobenzenesulfonyl isocyanate", [E[3]], None),
    "mol_mcf3snco": ("m-CF3C6H4SO2NCO", "m-trifluoromethylbenzenesulfonyl isocyanate", [E[3]], None),
    "mol_bsnco": ("PhSO2NCO", "benzenesulfonyl isocyanate", [E[3]], None),
    "mol_pmosnco": ("p-MeOC6H4SO2NCO", "p-methoxybenzenesulfonyl isocyanate", [E[3]], None),
    "mol_ethyl_allenoate": ("ethyl 2,3-butadienoate", "ethyl 2,3-butadienoate", [E[3]], "CCOC(=O)C=C=C"),
    "mol_paraformaldehyde": ("paraformaldehyde", "paraformaldehyde", [E[3]], None),
    "mol_acetaldehyde": ("acetaldehyde", "acetaldehyde", [E[3]], "CC=O"),
    "mol_hydroxyquinuclidine": ("3-hydroxyquinuclidine", "3-hydroxyquinuclidine", [E[3]], None),
    "mol_co2": ("CO2", "carbon dioxide", [E[2], E[4]], "O=C=O"),
}
for mid, (lab, name, ev, smi) in simple.items():
    add_molecule(mid, lab, ev, description=name, smiles=smi,
                 chirality="no source-specified stereochemistry for the named reagent",
                 geometry="not applicable or not source-specified",
                 representation="complete for source label" if smi else "text identity only; structure not encoded",
                 identity="source_identified")

# Drawn substrate/product series.  Terminal allene carbon is CH2; most b--m
# alcohols/carbamates and 3b--m products contain an unassigned stereocenter.
descriptors = {
    "a": "R2 = H; ethyl ester", "b": "R2 = Me; ethyl ester",
    "c": "R2 = Et; ethyl ester", "d": "R2 = n-Pr; ethyl ester",
    "e": "R2 = cyclopropyl; ethyl ester", "f": "R2 = cyclopentyl; ethyl ester",
    "g": "R2 = cyclohexyl; ethyl ester", "h": "R2 = phenyl; ethyl ester",
    "i": "R2 = but-3-ynyl; ethyl ester", "j": "R2 = phenethyl; ethyl ester",
    "k": "R2 = 3-(tert-butyldimethylsilyloxy)propyl; ethyl ester",
    "l": "R2 = 3-ethoxy-3-oxopropyl; ethyl ester",
    "m": "R2 = 3-ethoxy-3-oxopropyl; tert-butyl allenoate/pyrroline ester",
}
for s, desc in descriptors.items():
    ch = "no stereogenic center shown" if s == "a" else "stereogenic center present; configuration unspecified"
    add_molecule(f"mol_7{s}", f"7{s}", [E[3]], description=f"allenylic alcohol 7{s}; {desc}",
                 chirality=ch, geometry="allene geometry shown without an E/Z descriptor",
                 representation="source drawing retained; exact SMILES not encoded",
                 axial="terminal allene carbon is CH2; allene is not axially stereogenic as drawn")
    add_molecule(f"mol_1{s}", f"1{s}", [E[2], E[3]] if s in "ab" else [E[3]],
                 description=f"allenylic tosyl carbamate 1{s}; {desc}", chirality=ch,
                 geometry="allene geometry shown without an E/Z descriptor",
                 representation="source drawing retained; exact SMILES not encoded",
                 axial="terminal allene carbon is CH2; allene is not axially stereogenic as drawn")
    pch = "no stereogenic center shown" if s == "a" else "pyrroline stereogenic center present; configuration unspecified"
    add_molecule(f"mol_3{s}", f"3{s}", [E[2], E[3]] if s in "ab" else [E[3]],
                 description=f"3-pyrroline 3{s}; {desc}", chirality=pch,
                 geometry="endocyclic alkene; no E/Z assignment applicable to the depicted five-membered ring",
                 representation="source drawing retained; exact stereochemical SMILES not encoded")

extra = {
    "ba": "p-nitrophenylsulfonyl", "bb": "m-trifluoromethylphenylsulfonyl",
    "bc": "phenylsulfonyl", "bd": "p-methoxyphenylsulfonyl",
}
for s, desc in extra.items():
    add_molecule(f"mol_1{s}", f"1{s}", [E[2], E[3]],
                 description=f"allenylic carbamate 1{s}; {desc}, methyl-substituted scaffold",
                 chirality="stereogenic center present; configuration unspecified",
                 geometry="allene geometry shown without an E/Z descriptor",
                 representation="source drawing retained; exact SMILES not encoded",
                 axial="terminal allene carbon is CH2; allene is not axially stereogenic as drawn")
    add_molecule(f"mol_3{s}", f"3{s}", [E[2], E[3]],
                 description=f"3-pyrroline 3{s}; {desc}, methyl-substituted scaffold",
                 chirality="pyrroline stereogenic center present; configuration unspecified",
                 geometry="endocyclic alkene; no E/Z assignment applicable to the depicted five-membered ring",
                 representation="source drawing retained; exact stereochemical SMILES not encoded")

# Mechanistic-study identities from Tables 4--6.
add_molecule("mol_9", "9", [E[4], E[5], E[6]], description="ethyl 2-(acetoxymethyl)-2,3-butadienoate",
             smiles="CCOC(=O)C(=C=C)COC(C)=O", chirality="no stereogenic center shown",
             geometry="terminal allene; no E/Z descriptor", representation="connectivity encoded; source identity remains needs_review",
             axial="terminal allene carbon is CH2; allene is not axially stereogenic as drawn")
add_molecule("mol_10a", "10a", [E[4], E[6]], description="benzyl toluenesulfonamide",
             smiles="Cc1ccc(S(=O)(=O)NCc2ccccc2)cc1", chirality="no stereogenic center shown",
             geometry="not applicable", representation="connectivity encoded; source identity remains needs_review")
add_molecule("mol_10b", "10b", [E[4], E[5], E[6]], description="3-methyl-2,4-pentanedione",
             smiles="CC(=O)C(C)C(C)=O", chirality="prochiral/enolizable center; no stable source assignment",
             geometry="tautomer/geometry not specified", representation="one neutral connectivity form encoded; tautomeric state not asserted")
add_molecule("mol_10c", "10c", [E[4], E[5]], description="benzoyl propionitrile",
             smiles="CC(C#N)C(=O)c1ccccc1", chirality="stereogenic center possible; configuration unspecified",
             geometry="not applicable", representation="connectivity encoded; stereochemistry not encoded")
for letter in "abc":
    for num, label in ((11, f"11{letter}"), (12, f"12{letter}")):
        add_molecule(f"mol_{label}", label, [E[4], E[5]], description=f"mechanistic-study product {label}",
                     chirality="source drawing has no assigned absolute configuration",
                     geometry="source alkene/allene geometry retained only in drawing",
                     representation="source drawing retained; exact SMILES not encoded")
ratios = {"a": "Z:E = 1.11:1 to 1.19:1", "b": "E:Z = 1.05:1 to 1.18:1", "c": "E:Z = 1.08:1"}
for letter in "abc":
    ev = [E[4], E[5]] if letter != "c" else [E[5]]
    add_molecule(f"mol_13{letter}", f"13{letter}", ev, description=f"allylic acetate 13{letter}, E/Z mixture",
                 chirality="no absolute configuration assigned",
                 geometry="E/Z mixture", representation="mixture retained; individual isomers not merged into a single configured SMILES",
                 mixture=ratios[letter], assignment="E/Z ratio order reported; individual structures shown with a wavy bond")

# Generic/source-mechanism identities and prior-work structures.
add_generic("mol_generic_1", "1 (generic)", [E[1], E[2], E[4]], "allenylic carbamate scaffold",
            [("R1", "ester alkyl group", ["Et", "t-Bu"]), ("R2", "substituent at the alcohol-derived carbon", ["H", "alkyl", "cycloalkyl", "phenyl", "functionalized alkyl"])],
            [(f"mol_1{s}", f"1{s}", descriptors[s]) for s in descriptors])
add_generic("mol_generic_2", "2 (generic)", [E[1], E[2], E[4]], "phosphonium diene intermediate",
            [("R1", "ester alkyl group", ["Et", "t-Bu"]), ("R2", "diene substituent", ["H", "substituted carbon group"]), ("R3", "phosphine substituent", ["Ph", "Bu"])], [])
add_generic("mol_generic_3", "3 (generic)", [E[1], E[2], E[4]], "3-pyrroline scaffold",
            [("R1", "ester alkyl group", ["Et", "t-Bu"]), ("R2", "2-substituent", ["H", "alkyl", "cycloalkyl", "phenyl", "functionalized alkyl"])],
            [(f"mol_3{s}", f"3{s}", descriptors[s]) for s in descriptors])
for mid, lab, desc, ev in [
    ("mol_4", "4", "Tong prior-work allenylic acetate", [E[2], E[4], E[6]]),
    ("mol_5", "5", "Tong proposed phosphonium diene intermediate", [E[2], E[4]]),
    ("mol_6", "6", "Tong prior-work cyclopentene product", [E[2], E[6]]),
    ("mol_8", "8", "toluenesulfonamide anion", [E[4]]),
    ("mol_19", "19", "benzoyl acetonitrile bisnucleophile", [E[2], E[6]]),
]:
    add_molecule(mid, lab, ev, description=desc, chirality="unresolved from source context",
                 geometry="unresolved from source context", representation="source drawing retained; exact SMILES not encoded")


def part(r, mid, role, ev, *, amount=None, equiv=None, molpct=None, volume=None, product_role=None):
    p = new_participant(f"{r['reaction_id']}:{mid}:{len(r['participants'])+1}", mid, role, ev)
    if amount is not None:
        p["charges"]["amount"] = F(str(amount), qty(amount, "mmol"), ev)
    if equiv is not None:
        p["charges"]["equivalents"] = F(str(equiv), qty(equiv, "equiv"), ev)
    if molpct is not None:
        p["charges"]["mol_percent"] = F(str(molpct), qty(molpct, "mol%"), ev)
    if volume is not None:
        p["charges"]["volume"] = F(str(volume), qty(volume, "mL"), ev)
    if product_role:
        p["product_role"] = F(product_role, product_role, ev)
    if equiv is not None or molpct is not None:
        reference = next((x["molecule_id"] for x in r["participants"] if x["role"] == "reactant"), mid)
        p["equivalence_reference"] = F(reference, reference, ev)
    r["participants"].append(p)
    return p["participant_id"]


def step(r, operation, ev, *, time=None, temperature=None, atmosphere=None, participants=()):
    s = new_step(f"{r['reaction_id']}:step{len(r['steps'])+1}", ev)
    s["participant_ids"] = list(participants)
    s["operation"] = F(operation, operation, ev)
    if time is not None:
        if isinstance(time, tuple):
            normalized_time = {"minimum": time[0], "maximum": time[1], "unit": "h"}
            raw_time = f"{time[0]}-{time[1]} h"
        else:
            normalized_time = qty(time, "h")
            raw_time = str(time)
        s["time"] = F(raw_time, normalized_time, ev)
    if temperature is not None:
        norm = qty(temperature, "°C") if isinstance(temperature, (int, float)) else temperature
        s["temperature"] = F(str(temperature), norm, ev)
    if atmosphere is not None:
        s["atmosphere"] = F(atmosphere, atmosphere, ev)
    r["steps"].append(s)


def outcome(r, kind, targets, value, ev, *, yield_type=None, basis=None, scope="reaction outcome"):
    o = new_outcome(f"{r['reaction_id']}:out{len(r['outcomes'])+1}", kind, targets, ev)
    if value is not None:
        norm = qty(value, "%") if isinstance(value, (int, float)) else value
        o["value"] = F(str(value), norm, ev)
    if yield_type is not None:
        o["yield_type"] = F(yield_type, yield_type, ev)
    if basis is not None:
        o["measurement_basis"] = F(basis, basis, ev)
    o["scope"] = F(scope, scope, ev)
    r["outcomes"].append(o)


def reaction(rid, ev, record_type, category, name, transformation, *, location, labels=(), unresolved=()):
    rid = f"{DOC}:{rid}"
    r = new_reaction(rid, DOC, ev)
    r["record_type"] = record_type
    category_map = {
        "generic_scheme": "other", "cited_literature": "other",
        "substrate_preparation": "other", "scope": "substrate_scope",
        "proposed_mechanism": "other",
    }
    categories = [category_map.get(x, x) for x in category]
    r["experiment_categories"] = F(category, categories, ev)
    r["location"] = F(location, location, ev)
    r["reaction_name"] = F(name, name, ev)
    r["transformation_type"] = F(transformation, transformation, ev)
    r["review"] = {
        "status": "needs_review", "method": "assistant_visual", "reviewer": None,
        "checked_fields": ["page region", "compound labels", "table values", "footnotes"],
        "unresolved": list(unresolved) or ["Exact atom mapping and complete structures are not source-verified."],
    }
    category_numbers = {"C-C": 2, "C-O": 3, "C-N": 5}
    changes = []
    normalized_labels = []
    for i, label in enumerate(labels, 1):
        other = label.split("-")[1]
        change_id = f"b{i}"
        change = {"change_id": change_id,
                  "atoms": [{"atom_map": i * 2 - 1, "element": "C"}, {"atom_map": i * 2, "element": other}],
                  "before_order": 0, "after_order": 1, "change_type": "formed"}
        changes.append(F(change, change, ev, "inference",
                         derivation="Placeholder atom-map indices identify the inferred net formed-bond class only; no full reaction mapping is asserted."))
        normalized_labels.append({"category_number": category_numbers[label], "code": label,
                                  "supporting_change_ids": [change_id]})
    r["classification"]["bond_changes"] = changes
    r["classification"]["labels"] = F(list(labels), normalized_labels, ev, "inference",
                                         derivation="Net formed-bond class inferred from source drawings; no complete atom mapping asserted.")
    r["classification"]["coverage"] = F("partial", "partial", ev, "inference",
                                               derivation="Only confidently visible net formed-bond classes are retained.")
    r["classification"]["limitations"] = ["No complete atom-mapped reaction SMILES; bond classification is inferred from the source drawing."]
    return r


reactions = []

# Scheme 1: current generic transformation and the displayed Tong precedent.
r = reaction("scheme_current_annulation", [E[1], E[2]], "generic_scheme", ["generic_scheme"],
             "Phosphine-catalyzed (4+1) annulative rearrangement", "allenylic carbamate to 3-pyrroline",
             location="Page 2, Scheme 1, top ('This work')", labels=["C-N"],
             unresolved=["R1/R2/R3 are source variables; this is not a concrete experiment."])
part(r, "mol_generic_1", "reactant", [E[2]]); part(r, "mol_pph3", "catalyst", [E[2]])
part(r, "mol_generic_2", "unknown", [E[2]]); part(r, "mol_generic_3", "product", [E[2]], product_role="desired product")
part(r, "mol_co2", "product", [E[2]], product_role="byproduct")
r["notes"].append(F("R1 = Ph or Bu; R2 = H, alkyl, phenyl; PR3 shown as a variable phosphine.",
                    "Generic source scope; no conditions or yield asserted.", [E[2]]))
reactions.append(r)

r = reaction("lit_tong_scheme1", [E[2]], "literature_example", ["cited_literature"],
             "Tong phosphine-catalyzed (4+1) annulation", "allenylic acetate and carbon bisnucleophile to cyclopentene",
             location="Page 2, Scheme 1, bottom ('Previous work by Tong')", labels=["C-C"],
             unresolved=["Cited-paper arrow only; no current-paper experiment or full literature procedure is asserted."])
for mid, role in [("mol_4", "reactant"), ("mol_19", "reactant"), ("mol_pph3", "catalyst"), ("mol_cs2co3", "reagent"), ("mol_5", "unknown"), ("mol_6", "product")]:
    part(r, mid, role, [E[2]])
reactions.append(r)

# Two MBH preparations explicitly summarized in page-3 prose.
for s, aldehyde, y in [("a", "mol_paraformaldehyde", 66), ("b", "mol_acetaldehyde", 79)]:
    rid = f"rxn_mbh_7{s}"
    r = reaction(rid, [E[3]], "specific_experiment", ["substrate_preparation"],
                 f"Preparation of allenylic alcohol 7{s}", "Morita-Baylis-Hillman reaction",
                 location=f"Page 3, prose below Tables 2/3, product 7{s}", labels=["C-C"],
                 unresolved=["Pages 1--6 do not report scale, solvent, temperature, time, workup, purification, or yield measurement basis."])
    part(r, "mol_ethyl_allenoate", "reactant", [E[3]])
    part(r, aldehyde, "reactant", [E[3]])
    pcat = part(r, "mol_hydroxyquinuclidine", "catalyst", [E[3]])
    part(r, f"mol_7{s}", "product", [E[3]], product_role="desired product")
    step(r, "perform Morita-Baylis-Hillman reaction in the presence of 3-hydroxyquinuclidine", [E[3]],
         participants=[r["participants"][0]["participant_id"], r["participants"][1]["participant_id"], pcat])
    outcome(r, "yield", [f"mol_7{s}"], y, [E[3]], yield_type="yield type unreported", basis=None)
    r["procedure"]["original_text"] = F("MBH reaction between ethyl 2,3-butadienoate and paraformaldehyde or acetaldehyde in the presence of 3-hydroxyquinuclidine", None, [E[3]])
    reactions.append(r)

# Carbamate preparations: 1a--m plus the four non-tosyl variants.
carb_yields = {"a": 92, "b": 94, "c": 84, "d": 69, "e": 88, "f": 86, "g": 88,
               "h": 88, "i": 89, "j": 72, "k": 72, "l": 92, "m": 95,
               "ba": 99, "bb": 79, "bc": 99, "bd": 96}
iso_for = {"ba": "mol_pnosnco", "bb": "mol_mcf3snco", "bc": "mol_bsnco", "bd": "mol_pmosnco"}
for s, y in carb_yields.items():
    ev = [E[3]]
    rid = f"rxn_carb_1{s}"
    r = reaction(rid, ev, "specific_experiment", ["substrate_preparation", "scope"],
                 f"Preparation of allenylic carbamate 1{s}", "alcohol carbamoylation with sulfonyl isocyanate",
                 location=f"Page 3, {'Table 2' if len(s)>1 else 'Table 3'}, entry {s}", labels=["C-O"],
                 unresolved=["Scale, workup and purification are not reported on pages 1--6; exact atom mapping is unresolved."])
    alcohol = "mol_7b" if len(s) > 1 else f"mol_7{s}"
    p1 = part(r, alcohol, "reactant", ev, equiv=1.0)
    p2 = part(r, iso_for.get(s, "mol_tsnco"), "reactant", ev, equiv=1.05)
    ps = part(r, "mol_dcm", "solvent", ev)
    part(r, f"mol_1{s}", "product", ev, product_role="desired product")
    r["solvent_conditions"]["composition"] = F("dry CH2Cl2", ["CH2Cl2"], ev)
    r["solvent_conditions"]["concentration"] = F("0.4 M", qty(0.4, "M"), ev)
    r["solvent_conditions"]["concentration_basis"] = missing(ev, "The source does not explicitly state the concentration denominator.")
    step(r, "add sulfonyl isocyanate to allenylic alcohol and stir", ev, time=1, temperature=0, participants=[p1, p2, ps])
    step(r, "allow solution to warm/stir at room temperature", ev, time=0.5, temperature="room temperature", participants=[p1, p2, ps])
    outcome(r, "yield", [f"mol_1{s}"], y, ev,
            yield_type="based on recovered starting material" if s == "bb" else "isolated yield",
            basis="isolated product" if s != "bb" else "based on recovered starting material")
    if s == "bb":
        outcome(r, "yield", ["mol_1bb"], 22, ev, yield_type="isolated yield", basis="isolated product")
        r["notes"].append(F("79% based on recovered starting material; isolated yield 22%", None, ev))
    reactions.append(r)

# Table 1 optimization. Entries 3 and 12 are the same displayed experiments
# reused as Table 3 entry 1 and Tables 2/3 entry 2, so those IDs are merged.
t1 = [
    (1, "a", "mol_pph3", "mol_dcm", 0.035, 19, False, True),
    (2, "a", "mol_pph3", "mol_dcm", 0.035, 47, False, False),
    (3, "a", "mol_pph3", "mol_mecn", 0.035, 69, False, False),
    (4, "b", "mol_pph3", "mol_mecn", 0.035, 17, False, False),
    (5, "b", "mol_pmeph2", "mol_mecn", 0.035, 33, False, False),
    (6, "b", "mol_pme3", "mol_mecn", 0.035, "trace", False, False),
    (7, "b", "mol_pet3", "mol_mecn", 0.035, 55, False, False),
    (8, "b", "mol_pbu3", "mol_mecn", 0.035, 55, False, False),
    (9, "b", "mol_pcy3", "mol_mecn", 0.035, 14, False, False),
    (10, "b", "mol_pbu3", "mol_mecn", 0.018, 66, False, False),
    (11, "b", "mol_pbu3", "mol_mecn", 0.007, 51, False, False),
    (12, "b", "mol_pbu3", "mol_mecn", 0.018, 74, True, False),
]
ann_by_suffix = {}
for entry, s, catalyst, solvent, conc, y, deaerated, one_portion in t1:
    rid = f"rxn_ann_3{s}" if entry in (3, 12) else f"rxn_opt_t1_e{entry:02d}"
    ev = [E[2], E[3]] if entry in (3, 12) else [E[2]]
    r = reaction(rid, ev, "specific_experiment", ["optimization"] + (["scope"] if entry in (3, 12) else []),
                 f"Phosphine-catalyzed annulation to 3{s}", "(4+1) annulative rearrangement",
                 location=f"Page 2, Table 1, entry {entry}" + ("; repeated on page 3" if entry in (3, 12) else ""),
                 labels=["C-N"], unresolved=["No workup/purification details on pages 1--6; complete atom mapping unresolved.",
                 "The printed concentration basis is not defined; 0.14 mmol with 4 mL + 4 mL also permits a different total-mixture concentration calculation."])
    psub = part(r, f"mol_1{s}", "reactant", ev, amount=0.14, equiv=1.0)
    pcat = part(r, catalyst, "catalyst", ev, amount=0.14, equiv=1.0)
    psol = part(r, solvent, "solvent", ev, volume=8.0)
    part(r, f"mol_3{s}", "product", ev, product_role="desired product")
    r["reaction_scale"] = F("1 (0.14 mmol)", qty(0.14, "mmol", reference_molecule_id=f"mol_1{s}"), ev)
    r["solvent_conditions"]["composition"] = F(molecules[solvent]["names"]["original"]["normalized_value"], [molecules[solvent]["names"]["original"]["normalized_value"]], ev)
    r["solvent_conditions"]["total_volume"] = F("4 mL + 4 mL", qty(8.0, "mL"), ev)
    r["solvent_conditions"]["concentration"] = F(f"{conc} M", qty(conc, "M"), ev)
    r["solvent_conditions"]["concentration_basis"] = missing(ev, "The table does not define whether concentration refers to each solution or the combined mixture.")
    if deaerated:
        step(r, "deaerate solvent by bubbling with argon", ev, time=1.5, atmosphere="argon", participants=[psol])
    if one_portion:
        step(r, "add allenylic carbamate to phosphine in one portion", ev, participants=[psub, pcat, psol])
    else:
        step(r, "add allenylic carbamate solution to phosphine solution", ev, time=24, participants=[psub, pcat, psol])
    step(r, "stir after addition", ev, time=4, participants=[psub, pcat, psol])
    if y == "trace":
        outcome(r, "trace", [f"mol_3{s}"], "trace", ev, yield_type="isolated yield", basis="isolated product")
        outcome(r, "decomposition", [f"mol_1{s}"], "mainly decomposition", [E[2]], basis="qualitative author statement")
    else:
        outcome(r, "yield", [f"mol_3{s}"], y, ev, yield_type="isolated yield", basis="isolated product")
    reactions.append(r)
    if entry in (3, 12):
        ann_by_suffix[s] = r

# Remaining Table 3 annulations plus Table 2 sulfonyl variants.  3a and 3b
# have already been created above and receive the page-3 source facts there.
ann_yields = {"a": 69, "b": 74, "c": 62, "d": 55, "e": 28, "f": 37, "g": 30,
              "h": 66, "i": 52, "j": 52, "k": 48, "l": 56, "m": 55,
              "ba": 37, "bb": 45, "bc": 58, "bd": 72}
for s, y in ann_yields.items():
    if s in ann_by_suffix:
        r = ann_by_suffix[s]
        r["experiment_categories"] = F(["optimization", "scope"], ["optimization", "substrate_scope"], [E[2], E[3]])
        r["notes"].append(F(f"The same {s} experiment/yield is repeated in page-3 scope tables and is not counted as a second physical experiment.", None, [E[2], E[3]], "inference",
                            derivation="Identical substrate/product, catalyst, solvent, dose, timing, deaeration state and isolated yield."))
        if s == "a":
            # Table 3 explicitly gives room temperature and 28 h.
            for st in r["steps"]:
                st["temperature"] = F("rt", "room temperature", [E[3]])
        else:
            for st in r["steps"]:
                st["temperature"] = F("rt", "room temperature", [E[3]])
        continue
    ev = [E[3]]
    rid = f"rxn_ann_3{s}"
    r = reaction(rid, ev, "specific_experiment", ["scope"],
                 f"Phosphine-catalyzed annulation to 3{s}", "(4+1) annulative rearrangement",
                 location=f"Page 3, {'Table 2' if len(s)>1 else 'Table 3'}, product 3{s}", labels=["C-N"],
                 unresolved=["No workup/purification details on pages 1--6; complete atom mapping unresolved."])
    psub = part(r, f"mol_1{s}", "reactant", ev, amount=0.14, equiv=1.0)
    pcat = part(r, "mol_pbu3", "catalyst", ev, amount=0.14, equiv=1.0)
    psol = part(r, "mol_mecn", "solvent", ev, volume=8.0)
    part(r, f"mol_3{s}", "product", ev, product_role="desired product")
    r["reaction_scale"] = F("1 (0.14 mmol)", qty(0.14, "mmol", reference_molecule_id=f"mol_1{s}"), ev)
    r["solvent_conditions"]["composition"] = F("CH3CN", ["CH3CN"], ev)
    r["solvent_conditions"]["total_volume"] = F("4 mL + 4 mL", qty(8.0, "mL"), ev)
    r["solvent_conditions"]["concentration"] = F("0.018 M", qty(0.018, "M"), ev)
    r["solvent_conditions"]["concentration_basis"] = missing(ev, "Printed basis not explicitly defined.")
    if s != "k":
        step(r, "deaerate solvent by bubbling with argon", ev, time=1.5, atmosphere="argon", participants=[psol])
    step(r, "add carbamate solution to phosphine solution at room temperature", ev, time=24, temperature="room temperature", participants=[psub, pcat, psol])
    step(r, "stir after addition at room temperature", ev, time=4, temperature="room temperature", participants=[psub, pcat, psol])
    outcome(r, "yield", [f"mol_3{s}"], y, ev, yield_type="isolated yield", basis="isolated product")
    reactions.append(r)

# Scheme 2 is a proposed mechanism, not an experiment.
r = reaction("mechanism_scheme2", [E[4]], "proposed_mechanism", ["proposed_mechanism"],
             "Proposed conversion through phosphonium diene 2", "mechanistic proposal for the annulative rearrangement",
             location="Page 4, Scheme 2", labels=["C-N"],
             unresolved=["The initial attack site (beta-prime versus gamma) is explicitly unresolved by the authors here."])
for mid, role in [("mol_generic_1", "reactant"), ("mol_pph3", "catalyst"), ("mol_generic_2", "unknown"), ("mol_8", "unknown"), ("mol_co2", "product"), ("mol_generic_3", "product")]:
    part(r, mid, role, [E[4]])
r["notes"].append(F("Authors propose displacement/elimination to phosphonium diene 2, CO2 and sulfonamide anion 8, followed by cyclization.", None, [E[4]]))
reactions.append(r)

# Tables 4--6: every row is an independent mechanistic experiment.
table_defs = {
    4: {"nu": "a", "solvent": "mol_mecn", "addition": 6, "post": (1, 1.5),
        "rows": [(20,130,None,77,0,0,0),(None,130,None,81,0,0,0),(20,None,None,0,20,15,6),(110,None,None,0,0,14,0),(20,None,130,0,19,14,6),(None,None,130,0,0,0,100)]},
    5: {"nu": "b", "solvent": "mol_benzene", "addition": 1, "post": (5, 5.5),
        "rows": [(20,130,None,0,42,19,0),(None,130,None,8,74,0,0),(20,None,None,0,34,19,9),(110,None,None,0,0,23,0),(20,None,130,0,34,18,7),(110,None,130,0,0,24,0),(None,None,130,0,0,0,100)]},
    6: {"nu": "c", "solvent": "mol_benzene", "addition": 1, "post": [8.5,29,72,9],
        "rows": [(20,130,None,18,26,0,5),(None,130,None,28,44,0,2),(20,None,None,0,7,0,2),(110,None,None,0,0,20,0)]},
}
for table, spec in table_defs.items():
    letter = spec["nu"]
    for idx, row in enumerate(spec["rows"], 1):
        pph, cs, nao, y11, y12, y13, rec9 = row
        if table in (4, 5):
            ev = [E[4]]
            if idx <= 4:
                ev.append(E[5])
            if idx >= 3:
                ev.append(E[6])
        else:
            ev = [E[5]]
        rid = f"rxn_mech_t{table}_e{idx:02d}"
        r = reaction(rid, ev, "specific_experiment", ["mechanistic_study", "control"],
                     f"Mechanistic reaction of 9 with 10{letter}, Table {table} entry {idx}",
                     "competition among SN2/SN2-prime/gamma-umpolung pathways",
                     location=f"Page {4 if table in (4,5) else 5}, Table {table}, entry {idx}",
                     labels=["C-N"] if letter == "a" else ["C-C"],
                     unresolved=["Product-specific isolation is reported, but complete atom mapping and any detection limits for zero entries are not reported."])
        p9 = part(r, "mol_9", "reactant", ev, amount=0.190, equiv=1.0)
        p10 = part(r, f"mol_10{letter}", "reactant", ev, equiv=1.2)
        if pph is not None:
            part(r, "mol_pph3", "catalyst" if pph == 20 else "reagent", ev, molpct=pph)
        if cs is not None:
            part(r, "mol_cs2co3", "reagent", ev, molpct=cs, equiv=1.3)
        if nao is not None:
            part(r, "mol_naoac", "additive", ev, molpct=nao, equiv=1.3)
        psol = part(r, spec["solvent"], "solvent", ev, volume=5.6)
        alt_psol = None
        if table == 5:
            alt_psol = part(r, "mol_mecn", "solvent", ev)
            r["review"]["unresolved"].append("Page-4 prose says the standard 10a/10b reactions used CH3CN, while Table 5 arrow and footnote specify dry benzene; neither value was selected as authoritative.")
        for prod, val in ((f"mol_11{letter}", y11), (f"mol_12{letter}", y12), (f"mol_13{letter}", y13)):
            part(r, prod, "product", ev, product_role="observed product" if val > 0 else "assayed but not observed (0% reported)")
        r["reaction_scale"] = F("9 (0.190 mmol)", qty(0.190, "mmol", reference_molecule_id="mol_9"), ev)
        solvent_label = "dry CH3CN" if table == 4 else "dry benzene"
        if table == 5:
            r["solvent_conditions"]["composition"] = fact(
                evidence=ev, status="conflict", verification="unreviewed",
                derivation="Page-4 prose and Table 5 disagree on solvent.",
                alternatives=[F("dry benzene (Table 5 arrow/footnote)", ["benzene"], ev),
                              F("CH3CN (page-4 prose describing 10a and 10b)", ["CH3CN"], ev)])
            r["notes"].append(F("Solvent conflict: Table 5 specifies dry benzene/PhH, whereas the adjacent prose says the standard 10a and 10b reactions were in CH3CN.", None, ev))
        else:
            r["solvent_conditions"]["composition"] = F(solvent_label, [solvent_label.replace("dry ", "")], ev)
        r["solvent_conditions"]["total_volume"] = F("2.8 mL + 2.8 mL", qty(5.6, "mL"), ev)
        r["solvent_conditions"]["concentration"] = missing(ev, "No concentration value is printed for Tables 4--6.")
        r["solvent_conditions"]["concentration_basis"] = missing(ev)
        step(r, "add solution of 9 by syringe pump to solution/suspension of pronucleophile and listed reagents", ev,
             time=spec["addition"], temperature="room temperature", participants=[p9, p10, psol] + ([alt_psol] if alt_psol else []))
        post = spec["post"][idx-1] if isinstance(spec["post"], list) else spec["post"]
        step(r, "stir after addition", ev, time=post, temperature="room temperature", participants=[p9, p10, psol])
        for num, val in ((11,y11),(12,y12),(13,y13)):
            outcome(r, "yield", [f"mol_{num}{letter}"], val, ev, yield_type="isolated yield", basis="isolated product")
        outcome(r, "recovery", ["mol_9"], rec9, ev, yield_type="recovered starting material", basis="isolated starting material")
        if (table, idx) in ((4,6),(5,7)):
            r["classification"]["bond_changes"] = []
            r["classification"]["labels"] = F([], [], ev, "not_applicable")
            r["classification"]["coverage"] = F("no observed transformation", "not_applicable", ev, "inference",
                                                  derivation="All three product columns are 0% and starting material 9 is recovered in 100% yield.")
            r["classification"]["limitations"] = ["No observed transformation; expected-product bond classes were not assigned."]
            outcome(r, "no_reaction", ["mol_9"], "no reaction", [E[6]] if table in (4,5) else ev,
                    basis="100% isolated recovery of 9 and zero reported product yields")
        if y13 > 0:
            ratio = ratios[letter]
            out = new_outcome(f"{rid}:out{len(r['outcomes'])+1}", "selectivity", [f"mol_13{letter}"], ev)
            out["value"] = F(ratio, ratio, ev)
            out["measurement_basis"] = missing(ev, "E/Z analytical method is not reported in the table footnote.")
            out["scope"] = F("isolated allylic acetate 13 mixture", "isolated allylic acetate 13 mixture", ev)
            r["outcomes"].append(out)
        reactions.append(r)

# Page 6 contains two distinct proposed mechanisms, not additional experiments.
r = reaction("mechanism_scheme3", [E[6]], "proposed_mechanism", ["proposed_mechanism"],
             "Mechanistic rationale for allenoate 12 without external base", "phosphonium-dienolate general-base/SN2 pathway",
             location="Page 6, Scheme 3", labels=["C-N", "C-C"],
             unresolved=["Nu is a generic nitrogen- or carbon-centered pronucleophile; intermediates 14--17 are proposed, not observed here."])
for mid, role in [("mol_9", "reactant"), ("mol_10a", "reactant"), ("mol_pph3", "catalyst"), ("mol_12a", "product")]:
    part(r, mid, role, [E[6]])
reactions.append(r)

r = reaction("mechanism_scheme4", [E[6]], "proposed_mechanism", ["proposed_mechanism"],
             "Plausible mechanism for Tong's (4+1) reaction", "phosphonium-dienolate and intramolecular gamma-umpolung pathway",
             location="Page 6, Scheme 4", labels=["C-C"],
             unresolved=["Intermediates 18 and 20--24 are proposed; the authors state multiple mechanisms may operate."])
for mid, role in [("mol_4", "reactant"), ("mol_19", "reactant"), ("mol_pph3", "catalyst"), ("mol_cs2co3", "reagent"), ("mol_6", "product")]:
    part(r, mid, role, [E[6]])
reactions.append(r)

bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": "page1_figure1", "document_id": DOC, "evidence_ids": [E[1]], "kind": "background_figure",
     "note": "Figure 1 shows representative pyrrolizidine alkaloids; it is not a current-paper synthetic experiment."},
    {"context_id": "page3_alcohol_scope_context", "document_id": DOC, "evidence_ids": [E[3]], "kind": "current_work_summary",
     "note": "Prose says 7c--m were prepared by tin(II)- or indium-mediated allenylation, but pages 1--6 provide no discrete row-wise conditions/yields; no fabricated specific experiments were created from that grouped statement."},
    {"context_id": "page4_prior_allenoate", "document_id": DOC, "evidence_ids": [E[4]], "kind": "cited_literature",
     "note": "A cited Tong reaction of 4 with butyl toluenesulfonamide and 20 mol% PPh3 is said to give an allenoate analogue of 12a in 76%; it is prior work, not a current-paper experiment."},
    {"context_id": "pages4_6_mechanistic_discussion", "document_id": DOC, "evidence_ids": [E[4], E[5], E[6]], "kind": "author_interpretation",
     "note": "Mechanistic prose interpreting Tables 4--6 is retained separately from the row-wise experiments and from proposed Schemes 2--4."},
    {"context_id": "pages1_6_characterization_check", "document_id": DOC, "evidence_ids": [E[p] for p in range(1,7)], "kind": "characterization_inventory",
     "note": "No standalone mass-spectrum peak table, molecular-weight table, or fragment-relative-abundance table appears on pages 1--6."},
]

out = ROOT / "runs/current/submitted.json"
out.write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
(ROOT / "runs/current/raw_model/source_mentions.json").write_text(
    json.dumps(source_mentions(dataset, ["1a", "3a", "1b", "3b", "9", "10a", "10b", "10c", "13a", "13b", "13c"]),
               ensure_ascii=False, indent=2) + "\n"
)
page_map = {}
for page in range(1, 7):
    eid = E[page]
    page_map[str(page)] = [r["reaction_id"] for r in reactions if eid in r["evidence_ids"]]
(ROOT / "runs/current/raw_model/record_ids_by_page.json").write_text(
    json.dumps(page_map, ensure_ascii=False, indent=2) + "\n"
)
tasks_path = ROOT / "runs/current/tasks.json"
tasks = json.loads(tasks_path.read_text())
for task in tasks:
    if task.get("document_id") == DOC and task.get("evidence_id") in {E[p] for p in range(1, 7)}:
        task["state"] = "attempted"
        task["reason"] = "Original page image inspected; reactions/contexts mapped in pages_0001_0006.json. Machine records remain needs_review."
tasks_path.write_text(json.dumps(tasks, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(reactions),
                  "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in reactions)}, ensure_ascii=False))
