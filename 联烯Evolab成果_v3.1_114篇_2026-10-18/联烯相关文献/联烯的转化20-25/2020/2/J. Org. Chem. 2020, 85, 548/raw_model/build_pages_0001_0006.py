#!/usr/bin/env python3
"""Build the bounded pages 1-6 v3.1 extraction from visually inspected source pages."""
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
SOURCE = OUT / "submitted.json"
# The bounded pages were visually inspected, so formal figure/table facts cite
# their rendered page images. Later pages are used only as extracted-text
# context for explicitly noted same-document conflicts.
P = {i: f"{DOC}:page{i}:{'image' if i <= 6 else 'text'}" for i in range(1, 17)}


def F(raw=None, norm=None, ev=(), method="original_extraction", status=None, derivation=None):
    verification = None
    if method in ("inference", "structure_derivation") and (raw is not None or norm is not None):
        verification = "inferred"
    return fact(raw, norm, ev, method, status, verification, derivation)


def NA(ev=()):
    return fact(evidence=ev, method="not_applicable", status="not_applicable", verification="missing")


def q(number, unit, raw=None, qualifier=None):
    value = {"number": number, "unit": unit, "raw_text": raw or f"{number} {unit}"}
    if qualifier:
        value["qualifier"] = qualifier
    return value


def qr(minimum, maximum, unit, raw):
    return {"minimum": minimum, "maximum": maximum, "unit": unit, "raw_text": raw}


def alt(raw, norm, ev):
    return F(raw, norm, ev)


def conflict(alternatives, ev):
    return fact(None, None, ev, "original_extraction", "conflict", "missing", None, alternatives)


bundle = read(SOURCE)
bundle["molecules"] = []
bundle["reactions"] = []
bundle["source_context"] = []
molecules = {}
stereo_checks = {}


def add_molecule(mid, label, ev, name=None, smiles=None, allene=False, geometry=None,
                 identity_note=None, aliases=None):
    ev = list(dict.fromkeys(ev))
    m = new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, ev)
    m["names"]["original"] = F(name or label, name or label, ev)
    m["names"]["english"] = F(name, name, ev) if name else fact(evidence=ev)
    m["names"]["abbreviations"] = F(aliases, aliases, ev) if aliases else fact(evidence=ev)
    m["names"]["aliases"] = fact(evidence=ev)
    m["cas_numbers"] = fact(evidence=ev)
    m["external_properties"] = []
    if smiles:
        diag = stereo_inspect(smiles)
        stereo_checks[mid] = diag
        m["representation_type"] = "smiles"
        m["smiles"] = F(label, diag["canonical_isomeric_smiles"], ev, "structure_derivation",
                         derivation="Connectivity candidate was transcribed from the source drawing or explicit chemical name and normalized with the bundled stereochemistry inspector; this is not source verification.")
        m["stereochemistry"]["isomeric_smiles"] = F(
            label, diag["canonical_isomeric_smiles"], ev, "structure_derivation",
            derivation="Bundled standard stereo inspection of the machine-transcribed candidate SMILES.")
        m["stereochemistry"]["identity_key"] = F(
            label, diag["stereo_graph_key"], ev, "structure_derivation",
            derivation="SHA-256 stereo graph key returned by the bundled standard stereo inspector.")
        unspecified = [x for x in diag["potential_stereo"] if x["specified"] == "Unspecified"]
        atom_unspecified = [x for x in unspecified if str(x["type"]).startswith("Atom_")]
        bond_unspecified = [x for x in unspecified if str(x["type"]).startswith("Bond_")]
        # RDKit's ordinary potential-stereo list does not certify allene axial
        # stereochemistry. Never let an otherwise tag-free allene candidate be
        # collapsed to achiral/fully encoded on that basis.
        chirality = "unspecified" if allene or atom_unspecified else "achiral"
        geom = "unspecified" if bond_unspecified else "not_applicable"
        rep = "partial" if allene or unspecified else "encoded"
        m["stereochemistry"]["atom_centers"] = F(
            label, diag["atom_centers"], ev, "structure_derivation",
            derivation="Potential and assigned tetrahedral centers returned by standard stereo inspection; missing tags remain unspecified.")
        m["stereochemistry"]["double_bonds"] = F(
            label, diag["double_bonds"], ev, "structure_derivation",
            derivation="Double-bond diagnostics returned by standard stereo inspection; this does not infer source-unreported E/Z.")
        m["stereochemistry"]["axial"] = F(
            "No axial descriptor stated", None, ev, "inference", "unreadable",
            "The toolkit does not certify allene/axial stereochemistry; no source descriptor is asserted.") if allene else NA(ev)
        m["stereochemistry"]["mixture"] = fact(evidence=ev)
        m["stereochemistry"]["assignment"] = F(
            label, {"tool_diagnostic": diag["potential_stereo"]}, ev, "structure_derivation",
            derivation="Graph diagnostic only; the source drawing/name remains authoritative.")
        m["stereochemistry"]["chirality_status"] = F(label, chirality, ev, "structure_derivation",
            derivation="Status follows standard stereo diagnostics for the candidate graph and does not convert absent source stereochemistry into a racemate.")
        m["stereochemistry"]["geometry_status"] = F(label, geometry or geom, ev, "structure_derivation",
            derivation="Status follows the explicit source descriptor when supplied, otherwise the candidate graph diagnostic.")
        m["stereochemistry"]["representation_status"] = F(label, rep, ev, "structure_derivation",
            derivation="Encoded only when the candidate graph has no unspecified potential stereogenic element; otherwise partial.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = fact(evidence=ev)
        m["generic_template"] = None
        if allene:
            m["stereochemistry"]["axial"] = F(
                "Allene axis shown without R_a/S_a or P/M descriptor", {"axis": "allene", "assignment": "unspecified"}, ev,
                "structure_derivation", derivation="Visual inspection identifies an allene axis, but the source gives no absolute axial descriptor.")
            chirality = "unspecified"
        else:
            m["stereochemistry"]["axial"] = NA(ev)
            chirality = "unresolved"
        m["stereochemistry"]["atom_centers"] = fact(evidence=ev)
        if geometry == "specified_E":
            m["stereochemistry"]["double_bonds"] = F(
                "(E) in the source compound heading", [{"assignment": "E", "location": "source name"}], ev)
            geom_status = "specified"
        else:
            m["stereochemistry"]["double_bonds"] = fact(evidence=ev)
            geom_status = "unresolved" if geometry is None else geometry
        m["stereochemistry"]["mixture"] = fact(evidence=ev)
        m["stereochemistry"]["assignment"] = F(
            "No additional stereochemical assignment is stated in the inspected region", None, ev,
            "inference", "unreadable", "The drawing/name was checked for wedges, E/Z, cis/trans and axial descriptors; none beyond separately recorded fields was assigned.")
        m["stereochemistry"]["chirality_status"] = F(
            label, chirality, ev, "inference", derivation="Complete stereochemical identity cannot be established from the bounded drawing/name; absence of a descriptor is not treated as achiral or racemic.")
        m["stereochemistry"]["geometry_status"] = F(
            label, geom_status, ev, "inference", derivation="Complete E/Z or ring-relative stereochemistry is not established from the bounded source unless explicitly recorded.")
        m["stereochemistry"]["representation_status"] = F(
            label, "not_encoded", ev, "inference", derivation="No complete source-verified isomeric SMILES was assigned.")
        m["stereochemistry"]["isomeric_smiles"] = fact(evidence=ev)
        m["stereochemistry"]["identity_key"] = fact(evidence=ev)
    note = identity_note or "Machine registry entry derived from the cited paper region; structure and stereochemical identity remain subject to expert review."
    m["identity_note"] = F(note, note, ev, "inference", derivation="Extraction decision for a machine-generated record.")
    molecules[mid] = m
    bundle["molecules"].append(m)
    return mid


def add_generic(mid, label, ev, fixed, variables):
    m = new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, ev)
    m["names"]["original"] = F(label, label, ev)
    m["representation_type"] = "generic_template"
    m["smiles"] = fact(evidence=ev)
    m["generic_template"] = {
        "fixed_structure": F(fixed, fixed, ev),
        "variables": [{
            "label": v,
            "attachment_points": F(f"attachment labeled {v}", f"source attachment labeled {v}", ev),
            "definition": F(f"{v} is an unspecified substituent in Scheme 2D", "unspecified substituent", ev),
            "allowed_values": fact(evidence=ev),
            "constraints": fact(evidence=ev),
        } for v in variables],
        "instance_mapping": [],
        "reason": F("Scheme 2D is an R-group drawing", "Generic template retained because the source does not define concrete substituents", ev),
    }
    m["cas_numbers"] = fact(evidence=ev)
    m["stereochemistry"]["atom_centers"] = fact(evidence=ev)
    m["stereochemistry"]["double_bonds"] = fact(evidence=ev)
    m["stereochemistry"]["axial"] = F(
        "Generic vinylallene axis is drawn without assignment", {"assignment": "unspecified"}, ev,
        "structure_derivation", derivation="R-group template contains an allene; substituent identities needed for axial chirality are unspecified.") if "vinylallene" in label else NA(ev)
    m["stereochemistry"]["mixture"] = fact(evidence=ev)
    m["stereochemistry"]["assignment"] = fact(evidence=ev)
    m["stereochemistry"]["chirality_status"] = F(label, "unresolved", ev, "inference", derivation="Generic substituents prevent a unique chirality assessment.")
    m["stereochemistry"]["geometry_status"] = F(label, "unresolved", ev, "inference", derivation="Generic substituents prevent a unique geometry assessment.")
    m["stereochemistry"]["representation_status"] = F(label, "not_encoded", ev, "inference", derivation="Generic template is intentionally not represented as a concrete isomeric SMILES.")
    m["stereochemistry"]["isomeric_smiles"] = fact(evidence=ev)
    m["stereochemistry"]["identity_key"] = fact(evidence=ev)
    m["identity_note"] = F(label, "Generic Scheme 2D identity; not a unique chemical substance", ev, "inference", derivation="Source uses variable substituents.")
    m["external_properties"] = []
    molecules[mid] = m
    bundle["molecules"].append(m)
    return mid


def add_named(mid, name, ev, smiles=None, aliases=None):
    return add_molecule(mid, name, ev, name=name, smiles=smiles, aliases=aliases,
                        identity_note="Named reagent/solvent identity is source-supported; any SMILES is a machine-derived connectivity candidate and not independently source-verified.")


numbered_names = {
    1: "ethynylcyclohex-1-ene (compound 1)",
    2: "1-(cyclohex-1-en-1-yl)-5-phenylpent-1-yn-3-ol",
    3: "5-(cyclohex-1-en-1-yl)hexa-3,4-dien-1-ylbenzene",
    4: "direct substitution product 4",
    5: "vinylallene 5",
    6: "isopyridine cycloadduct 6",
    7: "3-ethyl-4-methyl-2-tosyl-5,6,7,8-tetrahydroquinoline",
    8: "lactam byproduct 8",
    9: "vinylallene 9",
    10: "(1-(cyclohex-1-en-1-yl)propa-1,2-dien-1-yl)benzene",
    11: "(E)-(3-ethylpenta-1,3,4-trien-1-yl)benzene",
    12: "3-ethyl-2-methylhepta-1,3,4-triene",
    13: "4-methyl-3-(3-phenylpropyl)-2-tosyl-5,6,7,8-tetrahydroquinoline",
    14: "3-ethyl-2-tosyl-5,6,7,8-tetrahydroquinoline",
    15: "3-methyl-4-phenyl-2-tosyl-5,6,7,8-tetrahydroquinoline",
    16: "4-ethyl-3-methyl-6-phenyl-2-tosylpyridine",
    17: "3-ethyl-4-methyl-2-propyl-6-tosylpyridine",
    18: "4-ethyl-5-methyl-3-propyl-2-tosylpyridine",
    23: "1-bromoprop-1-ene reagent 23 (source geometry unresolved)",
    24: "1-(1-methoxybut-2-yn-1-yl)cyclohex-1-ene",
    25: "vinylallenyl ether 25",
    26: "3-ethyl-4-methoxy-2-tosyl-5,6,7,8-tetrahydroquinoline",
    27: "cyclohexylacetylene (compound 27)",
    28: "1-(3-cyclohexyl-1-methoxyprop-2-yn-1-yl)cyclohex-1-ene",
    29: "vinylallenyl ether 29",
    30: "3-(cyclohexylmethyl)-4-methoxy-2-tosyl-5,6,7,8-tetrahydroquinoline",
    31: "enynyl alcohol 31",
    32: "ethyl 3-(cyclohex-1-en-1-yl)hexa-3,4-dienoate",
    33: "3-(cyclohex-1-en-1-yl)-N,N-dimethylhexa-3,4-dienamide",
    34: "2-(3-ethyl-2-tosyl-5,6,7,8-tetrahydroquinolin-4-yl)-N,N-dimethylacetamide",
    35: "ethyl 2-(3-ethyl-2-tosyl-5,6,7,8-tetrahydroquinolin-4-yl)acetate",
    36: "tert-butyl 3-formyl-1H-pyrrole-1-carboxylate",
    37: "tert-butyl 3-(1-methoxy-3-phenylprop-2-yn-1-yl)-1H-pyrrole-1-carboxylate",
    38: "vinylallenyl pyrrole 38",
    39: "tert-butyl 5-benzyl-4-methoxy-6-(toluenesulfonyl)-1H-pyrrolo[2,3-b]pyridine-1-carboxylate",
    40: "tert-butyl 5-benzyl-4-methoxy-6-(phenylsulfonyl)-1H-pyrrolo[2,3-b]pyridine-1-carboxylate",
    41: "N-silyl-3-iodopyrrole 41",
    42: "tert-butyl 3-(3-hydroxyprop-1-ynyl)-1H-pyrrole-1-carboxylate",
    43: "tert-butyl 3-(5-(dimethylamino)-5-oxopenta-1,2-dien-3-yl)-1H-pyrrole-1-carboxylate",
    44: "tert-butyl 4-(2-(dimethylamino)-2-oxoethyl)-5-methyl-6-tosyl-1H-pyrrolo[2,3-b]pyridine-1-carboxylate",
}

page_for = {}
for n in range(1, 19):
    page_for[n] = [P[3], P[4]]
for n in range(23, 31):
    page_for[n] = [P[5]]
for n in range(31, 45):
    page_for[n] = [P[6]]
page_for[11] += [P[8]]

smiles_candidates = {
    1: "C#CC1=CCCCC1",
    24: "CC#CC(OC)C1=CCCCC1",
    27: "C#CC1CCCCC1",
    28: "COC(C#CC1CCCCC1)C1=CCCCC1",
    31: "CC(O)C#CC1=CCCCC1",
    36: "CC(C)(C)OC(=O)n1cc(C=O)cc1",
    37: "COC(C#Cc1ccccc1)c1ccn(C(=O)OC(C)(C)C)c1",
    42: "CC(C)(C)OC(=O)n1cc(C#CCO)cc1",
}
allene_numbers = {3, 5, 9, 10, 11, 12, 25, 29, 32, 33, 38, 43}
for n, name in numbered_names.items():
    add_molecule(f"{DOC}:mol_{n}", str(n), page_for[n], name=name,
                 smiles=smiles_candidates.get(n), allene=n in allene_numbers,
                 geometry="specified_E" if n == 11 else None)

add_molecule(f"{DOC}:mol_2_benzoate", "benzoate derivative of 2", [P[3]],
             name="corresponding benzoate derivative of propargylic alcohol 2")

all_ev = [P[i] for i in range(2, 7)]
reagents = {
    "tscn": ("tosyl cyanide", "Cc1ccc(S(=O)(=O)C#N)cc1", [P[2], P[3], P[4], P[5], P[6]], ["TsCN"]),
    "toluene": ("toluene", "Cc1ccccc1", [P[3], P[4], P[6]], None),
    "thf": ("tetrahydrofuran", "C1CCOC1", all_ev, ["THF"]),
    "ether": ("diethyl ether", "CCOCC", [P[5]], ["Et2O"]),
    "ms4a": ("4 Å molecular sieves", None, [P[3], P[4], P[5], P[6]], ["4 Å MS"]),
    "dbu": ("1,8-diazabicyclo[5.4.0]undec-7-ene", "C1CCC2=NCCCN2CC1", [P[3], P[4], P[5], P[6]], ["DBU"]),
    "kotbu": ("potassium tert-butoxide", "CC(C)(C)[O-].[K+]", [P[3], P[6]], ["KOt-Bu"]),
    "caco3": ("calcium carbonate", "O=C([O-])[O-].[Ca+2]", [P[3]], None),
    "dipea": ("N-ethyl-N-isopropylpropan-2-amine", "CCN(C(C)C)C(C)C", [P[3]], ["(i-Pr)2EtN"]),
    "dtbp": ("2,6-di-tert-butylpyridine", None, [P[3]], None),
    "nbuli": ("n-butyllithium", None, [P[3], P[5], P[6]], ["n-BuLi"]),
    "phenylpropanal": ("3-phenylpropionaldehyde", "O=CCCc1ccccc1", [P[3]], None),
    "memgbr": ("methylmagnesium bromide", None, [P[3]], ["CH3MgBr"]),
    "mscl": ("methanesulfonyl chloride", "CS(=O)(=O)Cl", [P[3]], ["MsCl"]),
    "cubr": ("copper(I) bromide", "Br[Cu]", [P[3]], ["CuBr"]),
    "libr": ("lithium bromide", "[Li+].[Br-]", [P[3]], ["LiBr"]),
    "cyclohexenal": ("cyclohex-1-ene-1-carbaldehyde", "O=CC1=CCCCC1", [P[5]], ["cyclohexenecarbaldehyde"]),
    "mei": ("methyl iodide", "CI", [P[5], P[6]], ["MeI"]),
    "tbuli": ("tert-butyllithium", None, [P[5]], ["t-BuLi"]),
    "tmeda": ("N,N,N',N'-tetramethylethylenediamine", "CN(C)CCN(C)C", [P[5]], ["TMEDA"]),
    "meoh": ("methanol", "CO", [P[5]], ["MeOH"]),
    "acetaldehyde": ("acetaldehyde", "CC=O", [P[6]], None),
    "orthoacetate": ("triethyl orthoacetate", "CC(OCC)(OCC)OCC", [P[6]], ["CH3C(OEt)3"]),
    "propionic_acid": ("propionic acid", "CCC(=O)O", [P[6]], ["EtCO2H"]),
    "xylenes": ("xylenes (isomer composition unreported)", None, [P[6]], None),
    "amide_acetal": ("N,N-dimethylacetamide dimethyl acetal", "CC(OC)(OC)N(C)C", [P[6]], ["CH3C(OMe)2NMe2"]),
    "li_phenylacetylide": ("lithium phenylacetylide", None, [P[6]], None),
    "phso2cn": ("benzenesulfonyl cyanide", "N#CS(=O)(=O)c1ccccc1", [P[6]], ["PhSO2CN"]),
    "tbaf": ("tetrabutylammonium fluoride", "CCCC[N+](CCCC)(CCCC)CCCC.[F-]", [P[6]], ["TBAF"]),
    "et3n": ("triethylamine", "CCN(CC)CC", [P[6]], ["Et3N"]),
    "boc2o": ("di-tert-butyl dicarbonate", "CC(C)(C)OC(=O)OC(=O)OC(C)(C)C", [P[6]], ["Boc2O"]),
    "dmap": ("4-(dimethylamino)pyridine", "CN(C)c1ccncc1", [P[6]], ["DMAP"]),
    "propargyl_alcohol": ("propargyl alcohol", "C#CCO", [P[6]], None),
    "pdpph3cl2": ("bis(triphenylphosphine)palladium(II) dichloride", None, [P[6]], ["(PPh3)2PdCl2"]),
    "cui": ("copper(I) iodide", "I[Cu]", [P[6]], ["CuI"]),
    "piperidine": ("piperidine", "C1CCNCC1", [P[6]], None),
}
for key, (name, smiles, ev, aliases) in reagents.items():
    add_named(f"{DOC}:mol_{key}", name, ev, smiles, aliases)

add_generic(f"{DOC}:mol_generic_vinylallene", "generic preformed vinylallene", [P[2]],
            "R1/R2-substituted allene conjugated to an R3/R4-substituted alkene", ["R1", "R2", "R3", "R4"])
add_generic(f"{DOC}:mol_generic_sulfonyl_cyanide", "generic arylsulfonyl cyanide", [P[2]],
            "ArSO2CN", ["Ar"])
add_generic(f"{DOC}:mol_generic_sulfonylpyridine", "generic 2-sulfonylpyridine product", [P[2]],
            "2-(ArSO2)-pyridine bearing substituents derived from R1-R4", ["Ar", "R1", "R2", "R3", "R4"])
add_generic(f"{DOC}:mol_generic_substituted_pyridine", "generic nucleophilic substitution product", [P[2]],
            "2-Nuc-pyridine bearing substituents derived from R1-R4", ["Nuc", "R1", "R2", "R3", "R4"])


def mid(x):
    return f"{DOC}:mol_{x}"


def set_classification(r, ev, formed=(), unresolved=False, limitation=None):
    c = r["classification"]
    c["mapped_reaction_smiles"] = fact(evidence=ev)
    c["model_votes"] = []
    c["limitations"] = [limitation or "No atom-mapped reaction was asserted; formed-bond element pairs are a partial visual derivation from the source drawing."]
    if unresolved:
        c["bond_changes"] = []
        c["labels"] = F("No unique product graph could be assigned", [], ev, "inference",
                          derivation="Complex mixture, generic variables, or insufficient structural identity prevents defensible bond classification.")
        c["coverage"] = F("unresolved", "unresolved", ev, "inference",
                            derivation="Source does not support a unique mapped transformation.")
        return
    changes = []
    labels = {}
    for i, (a, b, order, change_type) in enumerate(formed, 1):
        change = {"change_id": f"b{i}", "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
                  "before_order": 0, "after_order": order, "change_type": change_type}
        if change_type.startswith("hydrogen_"):
            change["hydrogen_count"] = 1
        changes.append(F("formed bond inferred from the depicted net transformation", change, ev, "structure_derivation",
                         derivation="Element-pair change was inferred from the source drawing; atom correspondence is unresolved."))
        lab = bond_label(a, b, order)
        labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(change["change_id"])
    c["bond_changes"] = changes
    c["labels"] = F("formed-bond labels inferred from source drawing", list(labels.values()), ev,
                      "structure_derivation", derivation="Labels follow bond-formation-20261016 from partial element-pair changes; no complete atom mapping is claimed.")
    c["coverage"] = F("partial", "partial", ev, "inference",
                        derivation="Only visually defensible net formed bonds are recorded; broken and bond-order changes remain unresolved.")


def reaction(rid, ev, category=None, record_type="specific_experiment", name=None,
             transformation=None, location=None, formed=(), unresolved_class=False):
    ev = list(dict.fromkeys(ev))
    r = new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = F(category, category, ev) if category is not None else NA(ev)
    r["location"] = F(location or rid, location or rid, ev)
    r["reaction_name"] = F(name, name, ev) if name else fact(evidence=ev)
    r["transformation_type"] = F(transformation, transformation, ev, "inference",
                                  derivation="Plain-language transformation class inferred from the source arrow/title.") if transformation else fact(evidence=ev)
    r["general_scheme_ids"] = []
    r["reaction_scale"] = fact(evidence=ev)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
                   "unresolved": ["Complete source-verified structures and atom mapping are not established; record remains machine-generated."]}
    r["completeness"] = "partial"
    set_classification(r, ev, formed, unresolved_class)
    bundle["reactions"].append(r)
    return r


def addp(r, pid, molecule, role, ev, stage="reaction", eq=None, mass=None, amount=None,
         volume=None, molpct=None, ref=None, product_role=None):
    p = new_participant(pid, mid(molecule), role, ev, stage)
    if product_role:
        p["product_role"] = F(product_role, product_role, ev)
    if eq is not None:
        p["charges"]["equivalents"] = F(eq[0], eq[1], ev)
    if mass is not None:
        p["charges"]["mass"] = F(mass[0], mass[1], ev)
    if amount is not None:
        p["charges"]["amount"] = F(amount[0], amount[1], ev)
    if volume is not None:
        p["charges"]["volume"] = F(volume[0], volume[1], ev)
    if molpct is not None:
        p["charges"]["mol_percent"] = F(molpct[0], molpct[1], ev)
    if ref:
        p["equivalence_reference"] = F(ref, mid(ref), ev)
    r["participants"].append(p)
    return p


def add_step(r, sid, ev, pids=(), op=None, temp=None, time=None, stage="reaction"):
    s = new_step(sid, ev, stage)
    s["participant_ids"] = list(pids)
    if op is not None:
        s["operation"] = F(op, op, ev)
    if temp is not None:
        s["temperature"] = F(temp[0], temp[1], ev)
    if time is not None:
        s["time"] = F(time[0], time[1], ev)
    r["steps"].append(s)
    return s


def add_yield(r, oid, targets, raw, norm, ev, basis=None, scope=None, kind="yield", role="product"):
    o = new_outcome(oid, kind, [mid(x) for x in targets], ev)
    o["product_role"] = F(role, role, ev)
    o["value"] = F(raw, norm, ev)
    o["yield_type"] = F("reported yield", "reported", ev) if kind == "yield" else NA(ev)
    if basis:
        o["measurement_basis"] = F(basis, basis, ev)
    if scope:
        o["scope"] = F(scope, scope, ev)
    r["outcomes"].append(o)
    return o


def add_note(r, text, ev, method="inference", derivation="Extraction decision or unresolved source relation."):
    r["notes"].append(F(text, text, ev, method, derivation=derivation if method != "original_extraction" else None))


CC = ("C", "C", 1, "formed")
CN = ("C", "N", 1, "formed")
CO = ("C", "O", 1, "formed")
OH = ("O", "H", 1, "hydrogen_gain")

# Scheme 2D: current-work generic sequence, intentionally separate from specific experiments.
r = reaction("scheme2d_generic", [P[2]], record_type="generic_scheme", name="[4 + 2] cycloadditions of preformed vinylallenes and sulfonyl cyanides",
             transformation="generic cycloaddition, isomerization, and nucleophilic substitution sequence",
             location="PDF page 2, Scheme 2D", unresolved_class=True)
for pid, molecule, role in [("generic_vinylallene", "generic_vinylallene", "reactant"),
                            ("generic_sulfonyl_cyanide", "generic_sulfonyl_cyanide", "reactant"),
                            ("generic_sulfonylpyridine", "generic_sulfonylpyridine", "product"),
                            ("generic_substituted_pyridine", "generic_substituted_pyridine", "product")]:
    addp(r, pid, molecule, role, [P[2]], product_role="intermediate" if pid == "generic_sulfonylpyridine" else "final product" if role == "product" else None)
add_step(r, "cycloaddition", [P[2]], ["generic_vinylallene", "generic_sulfonyl_cyanide"], "[4 + 2] cycloaddition", None, None)
add_step(r, "isomerization", [P[2]], ["generic_sulfonylpyridine"], "isomerization", None, None)
add_step(r, "nucleophilic_substitution", [P[2]], ["generic_sulfonylpyridine", "generic_substituted_pyridine"], "nucleophilic substitution", None, None)

# Scheme 3.
r = reaction("scheme3_1_to_2", [P[3]], ["other"], name="propargylic alcohol synthesis",
             transformation="acetylide addition to aldehyde", location="PDF page 3, Scheme 3, first arrow", formed=[CC, OH])
addp(r, "s3_1", "1", "reactant", [P[3]])
addp(r, "s3_nbuli", "nbuli", "reagent", [P[3]], eq=("1.0 equiv", q(1.0, "equiv", "1.0 equiv")), ref="1")
addp(r, "s3_thf", "thf", "solvent", [P[3]])
addp(r, "s3_aldehyde", "phenylpropanal", "reactant", [P[3]], eq=("1.0 equiv", q(1.0, "equiv", "1.0 equiv")), ref="1")
addp(r, "s3_2", "2", "product", [P[3]], product_role="main product")
add_step(r, "deprotonation", [P[3]], ["s3_1", "s3_nbuli", "s3_thf"], "treat alkyne with n-BuLi", ("−40 °C", q(-40, "°C", "−40 °C")), ("30 min", q(30, "min", "30 min")))
add_step(r, "aldehyde_addition", [P[3]], ["s3_aldehyde"], "add 3-phenylpropionaldehyde and warm", ("−40 to 0 °C", "warming program"), None)
add_yield(r, "yield_2", ["2"], "94-99%", qr(94, 99, "%", "94-99%"), [P[3]], scope="compound 1 to compound 2; measurement basis not stated on page 3")

r = reaction("scheme3_2_to_3", [P[3]], ["other"], name="organocopper substitution synthesis of vinylallene 3",
             transformation="SN2'-type organocopper substitution", location="PDF page 3, Scheme 3, second arrow", formed=[CC])
addp(r, "s3b_2", "2", "reactant", [P[3]])
addp(r, "s3b_megbr1", "memgbr", "reagent", [P[3]], eq=("1.0 equiv", q(1.0, "equiv", "1.0 equiv")), ref="2")
addp(r, "s3b_mscl", "mscl", "reagent", [P[3]], eq=("1.0 equiv", q(1.0, "equiv", "1.0 equiv")), ref="2")
for pid, mol in [("s3b_megbr2", "memgbr"), ("s3b_cubr", "cubr"), ("s3b_libr", "libr")]:
    addp(r, pid, mol, "reagent", [P[3]], eq=("2.3 equiv", q(2.3, "equiv", "2.3 equiv")), ref="2")
addp(r, "s3b_thf", "thf", "solvent", [P[3]])
addp(r, "s3b_3", "3", "product", [P[3]], product_role="main product")
add_step(r, "mesylate_formation_1", [P[3]], ["s3b_2", "s3b_megbr1", "s3b_thf"], "treat alcohol with CH3MgBr", ("0 °C", q(0, "°C", "0 °C")), ("30 min", q(30, "min", "30 min")))
add_step(r, "mesylate_formation_2", [P[3]], ["s3b_mscl"], "add MsCl then warm", ("−78 °C to rt", "warming program"), None)
add_step(r, "organocopper_substitution", [P[3]], ["s3b_megbr2", "s3b_cubr", "s3b_libr"], "organocopper substitution", ("−78 °C to rt", "warming program"), ("1 h", q(1, "h", "1 h")))
add_yield(r, "yield_3", ["3"], "70%", q(70, "%", "70%"), [P[3]], scope="compound 2 to compound 3; measurement basis not stated in Scheme 3")

r = reaction("scheme3_benzoate_mixture", [P[3]], ["other"], name="organocopper reaction of the corresponding benzoate",
             transformation="competing SN2' and direct substitution", location="PDF page 3 prose below Scheme 3", formed=[CC])
addp(r, "s3c_benzoate", "2_benzoate", "reactant", [P[3]])
for pid, mol in [("s3c_megbr", "memgbr"), ("s3c_cubr", "cubr"), ("s3c_libr", "libr")]:
    addp(r, pid, mol, "reagent", [P[3]])
addp(r, "s3c_3", "3", "product", [P[3]], product_role="SN2' product")
addp(r, "s3c_4", "4", "product", [P[3]], product_role="direct substitution product")
add_step(r, "same_conditions", [P[3]], ["s3c_benzoate", "s3c_megbr", "s3c_cubr", "s3c_libr"], "reaction under the same organocopper conditions described for Scheme 3", None, None)
add_yield(r, "yield_3", ["3"], "35%", q(35, "%", "35%"), [P[3]])
add_yield(r, "yield_4", ["4"], "25%", q(25, "%", "25%"), [P[3]])
add_note(r, "The exact benzoate charge and whether every Scheme 3 operation was repeated are not enumerated in this prose sentence.", [P[3]])


def setup_cycloaddition(r, substrate, ev, tscn_eq=1.1, solvent="toluene", ms=False, dbu=None, base=None):
    addp(r, f"{r['reaction_id']}:substrate", substrate, "reactant", ev)
    addp(r, f"{r['reaction_id']}:tscn", "tscn", "reactant", ev,
         eq=(f"{tscn_eq} equiv", q(tscn_eq, "equiv", f"{tscn_eq} equiv")), ref=substrate)
    addp(r, f"{r['reaction_id']}:solvent", solvent, "solvent", ev)
    if ms:
        # The schema has no ratio-valued charge slot. Preserve the source ratio as
        # the raw mass field and leave its normalized quantity unresolved instead
        # of misrepresenting the ratio as a mass.
        addp(r, f"{r['reaction_id']}:ms", "ms4a", "additive", ev,
             mass=("1:1 w/w sieves/vinylallene", None), ref=substrate)
    if dbu:
        raw, norm = dbu
        addp(r, f"{r['reaction_id']}:dbu", "dbu", "reagent", ev, eq=(raw, norm), ref=substrate)
    if base:
        mol, raw, norm = base
        addp(r, f"{r['reaction_id']}:base", mol, "reagent", ev, eq=(raw, norm), ref=substrate)
    r["solvent_conditions"]["composition"] = F(reagents[solvent][0], [reagents[solvent][0]], ev)
    if solvent == "toluene":
        r["solvent_conditions"]["concentration"] = F("0.1 M", q(0.1, "M", "0.1 M"), ev)
        r["solvent_conditions"]["concentration_basis"] = F("vinylallene", mid(substrate), ev)


# Table 1, every row retained independently.
table1 = {
    1: {"temp": ("reflux", "reflux"), "time": ("1 h", q(1, "h", "1 h")), "other": "complex mixture"},
    2: {"temp": ("rt", "rt"), "time": ("72 h", q(72, "h", "72 h")), "other": "complex mixture"},
    3: {"temp": ("reflux", "reflux"), "time": ("1 h", q(1, "h", "1 h")), "base": ("caco3", "1.0-1.1 equiv", qr(1.0, 1.1, "equiv", "1.0-1.1 equiv")), "other": "complex mixture"},
    4: {"temp": ("50 °C", q(50, "°C", "50 °C")), "time": ("15 h", q(15, "h", "15 h")), "base": ("caco3", "1.0-1.1 equiv", qr(1.0, 1.1, "equiv", "1.0-1.1 equiv")), "products": [("6", "trace", "trace"), ("8", "51%", q(51, "%", "51%"))]},
    5: {"temp": ("50 °C", q(50, "°C", "50 °C")), "time": ("15 h", q(15, "h", "15 h")), "base": ("dipea", "2 equiv", q(2, "equiv", "2 equiv")), "other": "complex mixture"},
    6: {"temp": ("50 °C", q(50, "°C", "50 °C")), "time": ("15 h", q(15, "h", "15 h")), "base": ("dtbp", "2 equiv", q(2, "equiv", "2 equiv")), "products": [("7", "mixture of 7 and 8", "mixture"), ("8", "mixture of 7 and 8", "mixture")]},
    7: {"temp": ("50 °C", q(50, "°C", "50 °C")), "time": ("15 h", q(15, "h", "15 h")), "ms": True, "products": [("6", "56%", q(56, "%", "56%")), ("8", "trace", "trace")]},
    8: {"temp": ("reflux", "reflux"), "time": ("5 h", q(5, "h", "5 h")), "ms": True, "other": "complex mixture"},
    9: {"temp": ("50 °C", q(50, "°C", "50 °C")), "time": ("20 h", q(20, "h", "20 h")), "ms": True, "base": ("dtbp", "2 equiv", q(2, "equiv", "2 equiv")), "products": [("6", "38%", q(38, "%", "38%")), ("7", "15%", q(15, "%", "15%"))]},
    10: {"temp": ("50 °C", q(50, "°C", "50 °C")), "time": ("15 h", q(15, "h", "15 h")), "ms": True, "post": ("dbu", "1.0-1.1 equiv", qr(1.0, 1.1, "equiv", "1.0-1.1 equiv")), "products": [("7", "55%", q(55, "%", "55%")), ("7", "68% using freshly purified vinylallene", q(68, "%", "68%"))]},
    11: {"temp": ("50 °C", q(50, "°C", "50 °C")), "time": ("15 h", q(15, "h", "15 h")), "ms": True, "post": ("kotbu", "1.0-1.1 equiv", qr(1.0, 1.1, "equiv", "1.0-1.1 equiv")), "products": [("7", "55%", q(55, "%", "55%"))]},
}
for entry, spec in table1.items():
    identified = bool(spec.get("products"))
    formed = [CC, CN]
    if entry in (4, 7):
        formed.append(CO)
    r = reaction(f"table1_entry_{entry}", [P[3]], ["optimization"], name="vinylallene-tosyl cyanide cycloaddition optimization",
                 transformation="hetero-Diels-Alder cycloaddition/isomerization optimization",
                 location=f"PDF page 3, Table 1, entry {entry}", formed=formed if identified else (), unresolved_class=not identified)
    setup_cycloaddition(r, "5", [P[3]], ms=spec.get("ms", False), base=spec.get("base"))
    if spec.get("post"):
        mol, raw, norm = spec["post"]
        addp(r, f"t1e{entry}_postbase", mol, "reagent", [P[3]], eq=(raw, norm), ref="5")
    add_step(r, f"entry_{entry}_cycloaddition", [P[3]], [p["participant_id"] for p in r["participants"] if p["role"] != "product"], "heat or stir under Table 1 conditions", spec["temp"], spec["time"])
    if spec.get("post"):
        add_step(r, f"entry_{entry}_base_isomerization", [P[3]], [f"t1e{entry}_postbase"], "add base after cooling", ("rt", "rt"), ("1 h", q(1, "h", "1 h")))
    if spec.get("other"):
        o = new_outcome(f"entry_{entry}_observation", "other", [], [P[3]])
        o["product_role"] = NA([P[3]])
        o["value"] = F(spec["other"], spec["other"], [P[3]])
        o["yield_type"] = NA([P[3]])
        o["scope"] = F("Table 1 product mixture", "unresolved mixture", [P[3]])
        r["outcomes"].append(o)
    for idx, (prod, raw, norm) in enumerate(spec.get("products", []), 1):
        if not any(p["molecule_id"] == mid(prod) for p in r["participants"]):
            addp(r, f"t1e{entry}_product_{prod}", prod, "product", [P[3]], product_role="reported product or byproduct")
        if raw == "trace":
            o = new_outcome(f"entry_{entry}_trace_{prod}", "trace", [mid(prod)], [P[3]])
            o["product_role"] = F("trace byproduct", "byproduct", [P[3]])
            o["value"] = F("trace", "trace", [P[3]])
            o["yield_type"] = NA([P[3]])
            r["outcomes"].append(o)
        elif norm == "mixture":
            o = new_outcome(f"entry_{entry}_mixture_{prod}", "other", [mid(prod)], [P[3]])
            o["product_role"] = F("component of mixture", "mixture component", [P[3]])
            o["value"] = F(raw, raw, [P[3]])
            o["yield_type"] = NA([P[3]])
            r["outcomes"].append(o)
        else:
            basis = "yield estimated by 1H NMR analysis; product approximately 80% pure" if entry == 7 and prod == "6" else None
            scope = "freshly purified vinylallene 5" if entry == 10 and raw.startswith("68%") else None
            add_yield(r, f"entry_{entry}_yield_{prod}_{idx}", [prod], raw, norm, [P[3]], basis=basis, scope=scope)
    if entry == 10:
        r["outcomes"][-1]["measurement_basis"] = F("isolated after chromatographic purification (prose)", "isolated", [P[3]])


# Table 2 substrate scope.
table2 = {
    1: ("3", ["13"], [("13", "68%", q(68, "%", "68%"))], ("90 °C", q(90, "°C", "90 °C")), ("3 h", q(3, "h", "3 h")), ("1 h", q(1, "h", "1 h"))),
    2: ("5", ["7"], [("7", "68%", q(68, "%", "68%"))], ("50 °C", q(50, "°C", "50 °C")), ("12 h", q(12, "h", "12 h")), ("1 h", q(1, "h", "1 h"))),
    3: ("9", ["14"], [("14", "45% in Table 2; 39% in experimental paragraph", None)], ("80 °C", q(80, "°C", "80 °C")), ("8 h", q(8, "h", "8 h")), ("1 h", q(1, "h", "1 h"))),
    4: ("10", ["15"], [("15", "42%", q(42, "%", "42%"))], ("95 °C", q(95, "°C", "95 °C")), ("8 h", q(8, "h", "8 h")), ("1 h", q(1, "h", "1 h"))),
    5: ("11", ["16"], [("16", "73% in Table 2; 72% in experimental paragraph", None)], ("50 °C", q(50, "°C", "50 °C")), ("24 h", q(24, "h", "24 h")), None),
    6: ("12", ["17", "18"], [("17", "25%", q(25, "%", "25%")), ("18", "13%", q(13, "%", "13%"))], ("80 °C", q(80, "°C", "80 °C")), ("8 h", q(8, "h", "8 h")), ("15 h", q(15, "h", "15 h"))),
}
for entry, (substrate, products, yields, temp, tm, dbutm) in table2.items():
    ev = [P[4]] + ([P[8]] if entry == 2 else [P[9]] if entry in (1, 3, 4, 5, 6) else [])
    r = reaction(f"table2_entry_{entry}", ev, ["substrate_scope"], name="cycloaddition of vinylallene and tosyl cyanide",
                 transformation="hetero-Diels-Alder cycloaddition followed by pyridine isomerization",
                 location=f"PDF page 4, Table 2, entry {entry}", formed=[CC, CN])
    setup_cycloaddition(r, substrate, ev, ms=True, dbu=("1.1 equiv", q(1.1, "equiv", "1.1 equiv")) if dbutm else None)
    for prod in products:
        addp(r, f"t2e{entry}_product_{prod}", prod, "product", ev, product_role="main product" if prod == products[0] else "regioisomeric product")
    if entry == 4:
        add_step(r, "slow_addition", ev, [f"{r['reaction_id']}:substrate"], "slowly add vinylallene to TsCN and sieves", temp, ("5 h", q(5, "h", "5 h")))
        add_step(r, "heated_hold", ev, [p["participant_id"] for p in r["participants"] if p["role"] != "product"], "continue heating after addition", temp, ("3 h", q(3, "h", "3 h")))
    else:
        add_step(r, "cycloaddition", ev, [p["participant_id"] for p in r["participants"] if p["role"] != "product"], "cycloaddition", temp, tm)
    if dbutm:
        add_step(r, "dbu_isomerization", ev, [f"{r['reaction_id']}:dbu"], "cool to room temperature, add DBU, and stir", ("rt", "rt"), dbutm)
    add_step(r, "purification", [P[4]], [], "column chromatography", None, None, "purification")
    for prod, raw, norm in yields:
        if norm is None and entry == 3:
            o = add_yield(r, f"yield_{prod}", [prod], raw, None, ev, basis="isolated after column chromatography")
            o["value"] = conflict([alt("45% (Table 2)", q(45, "%", "45%"), [P[4]]), alt("39% (experimental paragraph)", q(39, "%", "39%"), [P[9]])], ev)
            r["review"]["unresolved"].append("Yield conflict: Table 2 reports 45%, whereas the experimental paragraph reports 39%; neither value was selected.")
        elif norm is None and entry == 5:
            o = add_yield(r, f"yield_{prod}", [prod], raw, None, ev, basis="isolated after column chromatography")
            o["value"] = conflict([alt("73% (Table 2)", q(73, "%", "73%"), [P[4]]), alt("72% (experimental paragraph)", q(72, "%", "72%"), [P[9]])], ev)
            r["review"]["unresolved"].append("Yield conflict: Table 2 reports 73%, whereas the experimental paragraph reports 72%; neither value was selected.")
        else:
            add_yield(r, f"yield_{prod}", [prod], raw, norm, [P[4]], basis="isolated after column chromatography")
    if entry == 1:
        add_note(r, "Table 2 also gives 60 °C for 12 h as an alternative to 90 °C for 3 h; the experimental paragraph describes 90 °C for 3 h, so the structured steps follow that explicit experiment without discarding the table alternative.", [P[4], P[9]])
    if entry == 6:
        o = new_outcome("regioisomer_ratio", "selectivity", [mid("17"), mid("18")], [P[4]])
        o["product_role"] = F("17:18", "expected product:regioisomer", [P[4]])
        o["value"] = F("66:34", {"ratio": [66, 34], "order": [mid("17"), mid("18")]}, [P[4]])
        o["yield_type"] = NA([P[4]])
        o["measurement_basis"] = fact(evidence=[P[4]])
        o["scope"] = F("regioisomeric product ratio 17:18", "17:18", [P[4]])
        r["outcomes"].append(o)
        r["review"]["unresolved"].append("The analytical method for the 66:34 regioisomer ratio is not reported on page 4.")


def sequence_reaction(rid, ev, substrate, product, category, name, transformation, location, formed,
                      reagents_with_eq, solvent, steps, yield_data, intermediates=()):
    r = reaction(rid, ev, category, name=name, transformation=transformation, location=location, formed=formed)
    addp(r, rid + "_substrate", substrate, "reactant", ev)
    for idx, (mol, raw, norm) in enumerate(reagents_with_eq, 1):
        addp(r, f"{rid}_reagent_{idx}", mol, "reagent", ev, eq=(raw, norm), ref=substrate)
    addp(r, rid + "_solvent", solvent, "solvent", ev)
    for inter in intermediates:
        addp(r, f"{rid}_intermediate_{inter}", inter, "product", ev, product_role="transient intermediate")
    addp(r, rid + "_product", product, "product", ev, product_role="main product")
    for sid, op, temp, tm in steps:
        add_step(r, sid, ev, [p["participant_id"] for p in r["participants"] if p["role"] != "product"], op, temp, tm)
    raw, norm, basis, scope = yield_data
    add_yield(r, "yield_" + product, [product], raw, norm, ev, basis=basis, scope=scope)
    return r


# Schemes 5 and 6.
sequence_reaction("scheme5_23_to_24", [P[5]], "23", "24", ["other"], "synthesis of propargylic ether 24",
                  "acetylide addition and O-methylation", "PDF page 5, Scheme 5, 23 to 24", [CC, CO],
                  [("nbuli", "2.2 equiv", q(2.2, "equiv", "2.2 equiv")), ("cyclohexenal", "1 equiv", q(1, "equiv", "1 equiv")), ("mei", "excess", "excess")],
                  "thf", [("deprotonation", "treat bromopropene-derived acetylide precursor with n-BuLi", ("−78 °C", q(-78, "°C", "−78 °C")), ("2 h", q(2, "h", "2 h"))),
                          ("aldehyde_addition", "add cyclohexenecarbaldehyde", ("−78 °C", q(-78, "°C", "−78 °C")), ("45 min", q(45, "min", "45 min"))),
                          ("methylation", "add excess MeI and warm", ("−78 °C to rt", "warming program"), ("18 h", q(18, "h", "18 h")))],
                  ("92%", q(92, "%", "92%"), None, "compound 23/cyclohexenecarbaldehyde to 24"))
r = sequence_reaction("scheme5_24_to_26", [P[5]], "24", "26", ["other"], "isomerization and cycloaddition of vinylallenyl ether 25",
                  "propargyl-to-allene isomerization, hetero-Diels-Alder cycloaddition, and pyridine isomerization", "PDF page 5, Scheme 5, 24 via 25 to 26", [CC, CN],
                  [("tbuli", "1.1 equiv", q(1.1, "equiv", "1.1 equiv")), ("tmeda", "1.1 equiv", q(1.1, "equiv", "1.1 equiv")),
                   ("meoh", "3.5 equiv", q(3.5, "equiv", "3.5 equiv")), ("tscn", "1 equiv", q(1, "equiv", "1 equiv")),
                   ("ms4a", "dose not stated in Scheme 5", None), ("dbu", "1 equiv", q(1, "equiv", "1 equiv"))],
                  "ether", [("allene_isomerization_low_temp", "treat 24 with t-BuLi/TMEDA", ("−78 °C", q(-78, "°C", "−78 °C")), None),
                            ("methanol_quench_and_warm", "add MeOH and warm", ("−78 °C to rt", "warming program"), ("1 h", q(1, "h", "1 h"))),
                            ("cycloaddition", "add TsCN and molecular sieves in Et2O/THF", ("rt", "rt"), ("1 h", q(1, "h", "1 h"))),
                            ("dbu_isomerization", "add DBU", ("rt", "rt"), ("1 h", q(1, "h", "1 h")))],
                  ("51% overall", q(51, "%", "51%"), None, "overall from propargylic ether 24 to pyridine 26"), intermediates=("25",))
r["solvent_conditions"]["composition"] = F("Et2O-THF", ["diethyl ether", "tetrahydrofuran"], [P[5]])

sequence_reaction("scheme6_27_to_28", [P[5]], "27", "28", ["other"], "synthesis of propargylic ether 28",
                  "acetylide addition and O-methylation", "PDF page 5, Scheme 6, 27 to 28", [CC, CO],
                  [("nbuli", "1.1 equiv", q(1.1, "equiv", "1.1 equiv")), ("cyclohexenal", "1 equiv", q(1, "equiv", "1 equiv")), ("mei", "excess", "excess")],
                  "thf", [("deprotonation", "treat cyclohexylacetylene with n-BuLi", ("−78 °C", q(-78, "°C", "−78 °C")), ("30 min", q(30, "min", "30 min"))),
                          ("aldehyde_addition", "add cyclohexenecarbaldehyde", ("−78 °C", q(-78, "°C", "−78 °C")), ("30 min", q(30, "min", "30 min"))),
                          ("methylation", "add excess MeI and warm", ("−78 °C to rt", "warming program"), ("16 h", q(16, "h", "16 h")))],
                  ("94%", q(94, "%", "94%"), None, "compound 27/cyclohexenecarbaldehyde to 28"))
r = sequence_reaction("scheme6_28_to_30", [P[5]], "28", "30", ["other"], "isomerization and cycloaddition of vinylallenyl ether 29",
                  "propargyl-to-allene isomerization, hetero-Diels-Alder cycloaddition, and pyridine isomerization", "PDF page 5, Scheme 6, 28 via 29 to 30", [CC, CN],
                  [("tbuli", "1.1 equiv", q(1.1, "equiv", "1.1 equiv")), ("tmeda", "1.1 equiv", q(1.1, "equiv", "1.1 equiv")),
                   ("meoh", "3.5 equiv", q(3.5, "equiv", "3.5 equiv")), ("tscn", "1.1 equiv", q(1.1, "equiv", "1.1 equiv")),
                   ("ms4a", "dose not stated in Scheme 6", None), ("dbu", "1.1 equiv", q(1.1, "equiv", "1.1 equiv"))],
                  "ether", [("allene_isomerization_low_temp", "treat 28 with t-BuLi/TMEDA", ("−78 °C", q(-78, "°C", "−78 °C")), None),
                            ("methanol_quench_and_warm", "add MeOH and warm", ("−78 °C to rt", "warming program"), ("1 h", q(1, "h", "1 h"))),
                            ("cycloaddition", "add TsCN and molecular sieves in Et2O/THF", ("rt", "rt"), ("2.5 h", q(2.5, "h", "2.5 h"))),
                            ("dbu_isomerization", "add DBU", ("rt", "rt"), ("14 h", q(14, "h", "14 h")))],
                  ("61% overall", q(61, "%", "61%"), None, "overall from propargylic ether 28 to pyridine 30"), intermediates=("29",))
r["solvent_conditions"]["composition"] = F("Et2O-THF", ["diethyl ether", "tetrahydrofuran"], [P[5]])

# Additional explicit page-5 condition for 33, distinct from Scheme 7's 70 °C experiment.
r = reaction("text_33_50c_24h", [P[5]], ["optimization"], name="lower-temperature cycloaddition of vinylallene 33",
             transformation="hetero-Diels-Alder cycloaddition followed by DBU isomerization", location="PDF page 5 prose above Synthesis of 7-Azaindoles", formed=[CC, CN])
setup_cycloaddition(r, "33", [P[5]], ms=True, dbu=("DBU amount unreported in this sentence", None))
addp(r, "text33_product34", "34", "product", [P[5]], product_role="desired pyridine")
add_step(r, "cycloaddition", [P[5]], [p["participant_id"] for p in r["participants"] if p["role"] != "product"], "cycloaddition at 50 °C", ("50 °C", q(50, "°C", "50 °C")), ("24 h", q(24, "h", "24 h")))
add_step(r, "dbu_treatment", [P[5]], [f"{r['reaction_id']}:dbu"], "treat with DBU", None, None)
add_yield(r, "yield_34", ["34"], "31%", q(31, "%", "31%"), [P[5]], scope="distinct 50 °C / 24 h condition; measurement basis not stated in this sentence")

# Scheme 7.
sequence_reaction("scheme7_1_to_31", [P[6]], "1", "31", ["other"], "synthesis of enynyl alcohol 31",
                  "acetylide addition to acetaldehyde", "PDF page 6, Scheme 7, compound 1 to 31", [CC, OH],
                  [("nbuli", "1.0 equiv", q(1.0, "equiv", "1.0 equiv")), ("acetaldehyde", "amount unreported in Scheme 7", None)],
                  "thf", [("deprotonation", "treat alkyne 1 with n-BuLi", ("−78 °C", q(-78, "°C", "−78 °C")), ("30 min", q(30, "min", "30 min"))),
                          ("acetaldehyde_addition", "add acetaldehyde and warm", ("−78 to 0 °C", "warming program"), None)],
                  ("99%", q(99, "%", "99%"), None, "compound 1 to enynyl alcohol 31"))
sequence_reaction("scheme7_31_to_32", [P[6]], "31", "32", ["other"], "Johnson orthoester rearrangement",
                  "sigmatropic rearrangement to vinylallene ester", "PDF page 6, Scheme 7, 31 to 32", [CC],
                  [("orthoacetate", "10 equiv", q(10, "equiv", "10 equiv")), ("propionic_acid", "cat.", "catalytic")],
                  "xylenes", [("rearrangement", "heat with triethyl orthoacetate and catalytic propionic acid", ("reflux", "reflux"), ("4 h", q(4, "h", "4 h")))],
                  ("57-64%", qr(57, 64, "%", "57-64%"), None, "compound 31 to vinylallene 32"))
sequence_reaction("scheme7_31_to_33", [P[6]], "31", "33", ["other"], "Eschenmoser rearrangement",
                  "sigmatropic rearrangement to vinylallenyl amide", "PDF page 6, Scheme 7, 31 to 33", [CC],
                  [("amide_acetal", "5 equiv", q(5, "equiv", "5 equiv"))], "toluene",
                  [("rearrangement", "heat with dimethylacetamide dimethyl acetal", ("reflux", "reflux"), ("3 h", q(3, "h", "3 h")))],
                  ("76-90%", qr(76, 90, "%", "76-90%"), None, "compound 31 to vinylallene 33"))

for substrate, product, temp, tm, dbutm, yld in [("32", "35", 90, 8, 15, 50), ("33", "34", 70, 5, 4 * 60, 53)]:
    r = reaction(f"scheme7_{substrate}_to_{product}", [P[6]], ["substrate_scope"], name="vinylallene-tosyl cyanide cycloaddition",
                 transformation="hetero-Diels-Alder cycloaddition followed by DBU isomerization",
                 location=f"PDF page 6, Scheme 7, {substrate} to {product}", formed=[CC, CN])
    setup_cycloaddition(r, substrate, [P[6]], solvent="thf", ms=True, dbu=("1.1 equiv", q(1.1, "equiv", "1.1 equiv")))
    addp(r, f"s7_product_{product}", product, "product", [P[6]], product_role="main product")
    add_step(r, "cycloaddition", [P[6]], [p["participant_id"] for p in r["participants"] if p["role"] != "product"], "cycloaddition", (f"{temp} °C", q(temp, "°C", f"{temp} °C")), (f"{tm} h", q(tm, "h", f"{tm} h")))
    add_step(r, "dbu_isomerization", [P[6]], [f"{r['reaction_id']}:dbu"], "add DBU at room temperature", ("rt", "rt"), ("15 min", q(15, "min", "15 min")) if substrate == "32" else ("4 h", q(4, "h", "4 h")))
    add_yield(r, f"yield_{product}", [product], f"{yld}%", q(yld, "%", f"{yld}%"), [P[6]])

# Scheme 8.
r = sequence_reaction("scheme8_36_to_37", [P[6], P[11]], "36", "37", ["other"], "synthesis of propargylic methyl ether 37",
                  "alkynyl addition and O-methylation", "PDF page 6, Scheme 8, 36 to 37", [CC, CO],
                  [("li_phenylacetylide", "1.1 equiv", q(1.1, "equiv", "1.1 equiv")), ("mei", "excess", "excess")], "thf",
                  [("alkynyl_addition", "add aldehyde 36 to lithium phenylacetylide", ("−78 °C", q(-78, "°C", "−78 °C")), ("30 min", q(30, "min", "30 min"))),
                   ("methylation", "add excess MeI and warm", ("−78 °C to rt", "warming program"), ("14 h", q(14, "h", "14 h")))],
                  ("92%", q(92, "%", "92%"), None, "aldehyde 36 to propargylic ether 37"))
p36 = next(p for p in r["participants"] if p["molecule_id"] == mid("36"))
p36["charges"]["equivalents"] = conflict([alt("1.5 equiv (Scheme 8)", q(1.5, "equiv", "1.5 equiv"), [P[6]]), alt("1.0 equiv (experimental paragraph)", q(1.0, "equiv", "1.0 equiv"), [P[11]])], [P[6], P[11]])
p36["equivalence_reference"] = F("lithium phenylacetylide / reaction basis", mid("li_phenylacetylide"), [P[6]])
r["review"]["unresolved"].append("Charge conflict: Scheme 8 labels compound 36 as 1.5 equiv, whereas the experimental paragraph uses aldehyde 36 as 1.0 equiv; both are retained.")
add_note(r, "The scheme condenses lithium phenylacetylide as a single 1.1-equiv reagent, while the experimental paragraph separately reports phenylacetylene and n-BuLi charges.", [P[6], P[11]])

for product, cyanide, baseeq, ev, steps, yld in [
    ("39", "tscn", 0.1, [P[6], P[11]], [("isomerization", "catalytic KOt-Bu at −78 °C", ("−78 °C", q(-78, "°C", "−78 °C")), ("30 min", q(30, "min", "30 min"))), ("cycloaddition", "add tosyl cyanide then warm", ("−78 °C to rt", "warming program"), ("1 h", q(1, "h", "1 h")))], 44),
    ("40", "phso2cn", 0.2, [P[6], P[12]], [("isomerization", "catalytic KOt-Bu at −78 °C", ("−78 °C", q(-78, "°C", "−78 °C")), ("30 min", q(30, "min", "30 min"))), ("cycloaddition_zero", "add benzenesulfonyl cyanide and warm to 0 °C", ("0 °C", q(0, "°C", "0 °C")), ("8 h", q(8, "h", "8 h"))), ("cycloaddition_rt", "warm to room temperature", ("rt", "rt"), ("14 h", q(14, "h", "14 h")))], 35),
]:
    r = reaction(f"scheme8_37_to_{product}", ev, ["substrate_scope"], name="allenylpyrrole cycloaddition",
                 transformation="base-promoted propargyl-to-allene isomerization and hetero-Diels-Alder cycloaddition",
                 location=f"PDF page 6, Scheme 8, 37 via 38 to {product}", formed=[CC, CN])
    addp(r, f"s8_{product}_37", "37", "reactant", ev)
    addp(r, f"s8_{product}_base", "kotbu", "reagent", ev, eq=(f"{baseeq} equiv", q(baseeq, "equiv", f"{baseeq} equiv")), ref="37")
    addp(r, f"s8_{product}_cyanide", cyanide, "reactant", ev, eq=("1.1 equiv", q(1.1, "equiv", "1.1 equiv")), ref="37")
    addp(r, f"s8_{product}_ms", "ms4a", "additive", ev)
    addp(r, f"s8_{product}_thf", "thf", "solvent", ev)
    addp(r, f"s8_{product}_toluene", "toluene", "solvent", ev)
    addp(r, f"s8_{product}_38", "38", "product", ev, product_role="transient allenyl ether")
    addp(r, f"s8_{product}_product", product, "product", ev, product_role="main product")
    r["solvent_conditions"]["composition"] = F("THF/toluene", ["tetrahydrofuran", "toluene"], ev)
    for sid, op, temp, tm in steps:
        add_step(r, sid, ev, [p["participant_id"] for p in r["participants"] if p["role"] != "product"], op, temp, tm)
    add_yield(r, f"yield_{product}", [product], f"{yld}%", q(yld, "%", f"{yld}%"), [P[6]])

# Scheme 9.
r = reaction("scheme9_41_to_42", [P[6]], ["other"], name="deprotection/Boc protection/Sonogashira sequence to 42",
             transformation="multistage protecting-group exchange and C-C coupling", location="PDF page 6, Scheme 9, 41 to 42", formed=[CC, CN])
addp(r, "s9_41", "41", "reactant", [P[6]])
for pid, mol, raw, norm in [
    ("s9_tbaf", "tbaf", "1.2 equiv", q(1.2, "equiv", "1.2 equiv")),
    ("s9_et3n", "et3n", "amount unreported in Scheme 9", None),
    ("s9_boc2o", "boc2o", "amount unreported in Scheme 9", None),
    ("s9_dmap", "dmap", "cat.", "catalytic"),
    ("s9_propargyl", "propargyl_alcohol", "2 equiv", q(2, "equiv", "2 equiv")),
    ("s9_pd", "pdpph3cl2", "0.07 equiv", q(0.07, "equiv", "0.07 equiv")),
    ("s9_cui", "cui", "0.15 equiv", q(0.15, "equiv", "0.15 equiv")),
    ("s9_piperidine", "piperidine", "5.0 equiv", q(5.0, "equiv", "5.0 equiv")),
]:
    addp(r, pid, mol, "catalyst" if pid in ("s9_pd", "s9_cui") else "reagent", [P[6]], eq=(raw, norm) if norm is not None else None, ref="41")
addp(r, "s9_thf", "thf", "solvent", [P[6]])
addp(r, "s9_42", "42", "product", [P[6]], product_role="main product")
add_step(r, "desilylation", [P[6]], ["s9_41", "s9_tbaf", "s9_thf"], "TBAF treatment", ("rt", "rt"), ("2 h", q(2, "h", "2 h")))
add_step(r, "boc_protection", [P[6]], ["s9_et3n", "s9_boc2o", "s9_dmap"], "add Et3N, Boc2O, and catalytic DMAP", ("rt", "rt"), ("1 h", q(1, "h", "1 h")))
add_step(r, "sonogashira_coupling", [P[6]], ["s9_propargyl", "s9_pd", "s9_cui", "s9_piperidine"], "couple with propargyl alcohol", ("rt", "rt"), ("1 h", q(1, "h", "1 h")))
add_yield(r, "yield_42", ["42"], "87% overall", q(87, "%", "87%"), [P[6]], scope="overall from 41 to 42 across the two numbered operations")

sequence_reaction("scheme9_42_to_43", [P[6]], "42", "43", ["other"], "Eschenmoser rearrangement to allenylpyrrole 43",
                  "sigmatropic rearrangement", "PDF page 6, Scheme 9, 42 to 43", [CC],
                  [("amide_acetal", "5 equiv", q(5, "equiv", "5 equiv"))], "toluene",
                  [("rearrangement", "heat with dimethylacetamide dimethyl acetal", ("reflux", "reflux"), ("4 h", q(4, "h", "4 h")))],
                  ("48%", q(48, "%", "48%"), None, "compound 42 to allenylpyrrole 43"))
r = reaction("scheme9_43_to_44", [P[6], P[12]], ["substrate_scope"], name="allenylpyrrole 43-tosyl cyanide cycloaddition",
             transformation="hetero-Diels-Alder cycloaddition and spontaneous isomerization", location="PDF page 6, Scheme 9, 43 to 44", formed=[CC, CN])
setup_cycloaddition(r, "43", [P[6], P[12]], solvent="toluene", ms=True)
addp(r, "s9_44", "44", "product", [P[6], P[12]], product_role="main product")
add_step(r, "cycloaddition", [P[6], P[12]], [p["participant_id"] for p in r["participants"] if p["role"] != "product"], "heat with TsCN and molecular sieves", ("60 °C", q(60, "°C", "60 °C")), ("72 h", q(72, "h", "72 h")))
o = add_yield(r, "yield_44", ["44"], "35% in Scheme 9; 34% in experimental paragraph", None, [P[6], P[12]])
o["value"] = conflict([alt("35% (Scheme 9)", q(35, "%", "35%"), [P[6]]), alt("34% (experimental paragraph)", q(34, "%", "34%"), [P[12]])], [P[6], P[12]])
r["review"]["unresolved"].append("Yield conflict: Scheme 9 reports 35%, whereas the experimental paragraph reports 34%; neither value was selected.")


bundle["source_context"] = [
    {"context_id": f"{DOC}:page1_metadata", "document_id": DOC, "kind": "paper_metadata",
     "evidence_ids": [P[1]], "title": "Synthesis of Highly Substituted Pyridines via [4 + 2] Cycloadditions of Vinylallenes and Sulfonyl Cyanides",
     "doi": "10.1021/acs.joc.9b02628", "note": "Title and abstract were visually inspected; not a reaction experiment."},
    {"context_id": f"{DOC}:scheme1_generic_strategies", "document_id": DOC, "kind": "generic_strategy_overview",
     "evidence_ids": [P[1]], "note": "Scheme 1 contains four conceptual cyano Diels-Alder aromatization pathways (elimination, cycloreversion, and two isomerization modes). No concrete experiment, amounts, or yield is reported, so no specific experiment was created."},
    {"context_id": f"{DOC}:scheme2_prior_work", "document_id": DOC, "kind": "prior_work_overview",
     "evidence_ids": [P[2]], "note": "Scheme 2A-C summarizes cited earlier intramolecular propargylic ene/cyano or imino Diels-Alder cascades and prior bimolecular variants. These are classified as cited prior work, not current-paper experiments."},
    {"context_id": f"{DOC}:scheme4_regiochemical_model", "document_id": DOC, "kind": "author_mechanistic_interpretation",
     "evidence_ids": [P[4]], "note": "Scheme 4 depicts proposed normal and alternate transition-state/regiochemical courses (19-22); it is an explanatory model, not a separately executed mechanistic experiment."},
    {"context_id": f"{DOC}:page5_unsuccessful_one_flask", "document_id": DOC, "kind": "qualitative_negative_observation",
     "evidence_ids": [P[5]], "note": "Attempts to combine sigmatropic rearrangement and Diels-Alder reaction in one flask were unsuccessful because of side reactions, but substrate identity, exact conditions, product and measurement basis are not specified; retained as context rather than fabricating a control experiment."},
    {"context_id": f"{DOC}:page5_ipso_discussion", "document_id": DOC, "kind": "forward_reference",
     "evidence_ids": [P[5]], "note": "Page 5 discusses ipso substitution products 45-52 and points to Table 3 on page 7. No Table 3 row or arrow is visible in the bounded page range, so those experiments are not claimed as covered here."},
    {"context_id": f"{DOC}:page6_general_procedures", "document_id": DOC, "kind": "general_procedure_context",
     "evidence_ids": [P[6]], "note": "General glassware, argon, transfer, concentration, TLC/column and material-purification statements begin on page 6. They were not silently inherited into every scheme record unless the scheme/table explicitly supplied the condition."},
]

raw = {
    "batch_id": "allene-stereo-v31-200-20261017",
    "document_id": DOC,
    "page_range": [1, 6],
    "visual_pages_inspected": [1, 2, 3, 4, 5, 6],
    "semantic_accuracy_verified": False,
    "record_count": len(bundle["reactions"]),
    "molecule_count": len(bundle["molecules"]),
    "decisions": {
        "scheme1": "generic conceptual strategies only; source_context",
        "scheme2_panels_a_to_c": "cited prior work; source_context",
        "scheme2_panel_d": "current-work generic scheme; formal generic_scheme record",
        "table1": "eleven independent optimization records",
        "table2": "six independent substrate-scope records; regioisomers 17/18 remain one experiment with product-specific outcomes",
        "scheme4": "author regiochemical model; source_context",
        "schemes5_to_9": "each preparative arrow or telescoped one-pot sequence recorded separately",
        "si": "not associated; explicit paper notice does not establish supplied SI",
    },
    "known_conflicts": [
        "Table 2 entry 3: 45% on page 4 versus 39% in page 9 experimental text",
        "Table 2 entry 5: 73% on page 4 versus 72% in page 9 experimental text",
        "Scheme 8 compound 36: 1.5 equiv on page 6 versus 1.0 equiv in page 11 experimental text",
        "Scheme 9 product 44: 35% on page 6 versus 34% in page 12 experimental text",
    ],
}

write(OUT / "submitted.json", bundle)
write(OUT / "raw_model" / "raw_extraction_pages_0001_0006.json", raw)
write(OUT / "raw_model" / "stereo_checks.json", stereo_checks)
print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]), "stereo_checks": len(stereo_checks)}, ensure_ascii=False))
