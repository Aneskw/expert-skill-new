#!/usr/bin/env python3
"""Build the source-specific v3.1 bundle for doc_3a22feeb08de83f6.

This script is deliberately source-specific.  It does not read historical reaction
records and does not infer supporting-information content.
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

DOC = "doc_3a22feeb08de83f6"
RUN = ROOT / "runs/current"


def ev(page):
    return [f"{DOC}:page{page}:text", f"{DOC}:page{page}:image"]


def uniq(items):
    return list(dict.fromkeys(items))


def f(raw=None, normalized=None, evidence=(), method="original_extraction", derivation=None,
      status=None):
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation,
                   status=status)


def fq(number, unit, evidence, raw=None, qualifier=None):
    value = {"number": number, "unit": unit}
    if qualifier:
        value["qualifier"] = qualifier
    if raw is not None:
        value["raw_text"] = raw
    return f(raw if raw is not None else f"{number} {unit}", value, evidence)


def inferred(raw, normalized, evidence, derivation):
    return f(raw, normalized, evidence, method="inference", derivation=derivation)


bundle = json.loads((RUN / "raw_model/standard_init.json").read_text())
bundle["molecules"] = []
bundle["reactions"] = []
bundle["source_context"] = []
molecules = {}
stereo_checks = {}


def add_molecule(key, label, evidence, *, name=None, smiles=None, kind="ordinary",
                 chirality="unspecified", geometry="unspecified", assignment=None,
                 mixture=None, axial=None):
    mid = f"{DOC}:mol_{key}"
    if mid in molecules:
        return mid
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    m["occurrences"][0]["evidence_ids"] = list(evidence)
    m["names"]["original"] = f(name or label, name or label, evidence)
    m["names"]["english"] = f(name or label, name or label, evidence)
    m["names"]["abbreviations"] = f(label, [label], evidence)
    m["names"]["aliases"] = f()
    m["cas_numbers"] = f()
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = f(
            f"Source drawing/name for {label}", smiles, evidence,
            method="structure_derivation",
            derivation="Candidate SMILES transcribed conservatively from the supplied source drawing or explicit chemical name; parsing is not source identity verification.",
        )
        diag = st.stereo_inspect(smiles)
        stereo_checks[mid] = diag
        m["stereochemistry"]["isomeric_smiles"] = f(
            None, diag["canonical_isomeric_smiles"], evidence,
            method="structure_derivation",
            derivation="Canonical isomeric SMILES computed from the candidate registry graph with the standard stereo routine.",
        )
        m["stereochemistry"]["identity_key"] = f(
            None, diag["stereo_graph_key"], evidence,
            method="structure_derivation",
            derivation="Stereo graph key computed from the candidate graph; it is not a sample identity or source verification.",
        )
        rep_status = "partial" if any(x["specified"] == "Unspecified" for x in diag["potential_stereo"]) else "encoded"
        if kind == "simple" and not diag["potential_stereo"]:
            chirality = "achiral"
            geometry = "not_applicable"
        m["stereochemistry"]["representation_status"] = inferred(
            "candidate graph inspected", rep_status, evidence,
            "The standard stereo routine checked the candidate graph; source identity and axial/relative assignments remain outside parsing verification.",
        )
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = f()
        m["stereochemistry"]["representation_status"] = inferred(
            "source structure not encoded", "not_encoded", evidence,
            "The compact drawing, generic structure, mixture, or stereochemical detail was not converted to a complete safe SMILES.",
        )
    m["stereochemistry"]["chirality_status"] = inferred(
        f"source stereochemistry for {label}", chirality, evidence,
        "Status assigned from the inspected drawing and explicit descriptors; absence of a descriptor is retained as unspecified, not racemic or achiral.",
    )
    m["stereochemistry"]["geometry_status"] = inferred(
        f"source geometry for {label}", geometry, evidence,
        "Status assigned from the inspected drawing and explicit descriptors; unassigned alkene/relative geometry remains unspecified.",
    )
    if assignment is not None:
        m["stereochemistry"]["assignment"] = f(assignment, assignment, evidence)
    if mixture is not None:
        m["stereochemistry"]["mixture"] = f(mixture, mixture, evidence)
    if axial is not None:
        m["stereochemistry"]["axial"] = f(axial, axial, evidence)
    if kind == "allene":
        m["stereochemistry"]["axial"] = f(
            "Substituted allene drawn without R_a/S_a or P/M descriptor",
            {"element": "allene axis", "assignment": "unassigned"}, evidence,
        )
        m["identity_note"] = f(
            "The allene axis was inspected; no axial configuration is assigned in the supplied pages.",
            "Potential allene axial stereochemistry is retained as unassigned and not encoded.", evidence,
        )
    elif kind == "simple":
        m["stereochemistry"]["assignment"] = inferred(
            "candidate graph has no potential stereo element", {"basis": "standard stereo graph diagnostic"}, evidence,
            "The candidate graph contains no potential tetrahedral, double-bond, or encoded axial stereogenic element; this applies only to the candidate graph.",
        )
    bundle["molecules"].append(m)
    molecules[mid] = m
    return mid


def add_template(key, label, evidence, fixed, variables):
    mid = f"{DOC}:mol_{key}"
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    m["names"]["original"] = f(label, label, evidence)
    m["names"]["english"] = f(label, label, evidence)
    m["names"]["abbreviations"] = f(label, [label], evidence)
    m["representation_type"] = "generic_template"
    m["smiles"] = f()
    tv = []
    for item in variables:
        tv.append({
            "label": item,
            "attachment_points": f(f"{item} attachment shown in source generic drawing", {"label": item, "attachment": "source-defined"}, evidence),
            "definition": f(f"source variable {item}", item, evidence),
            "allowed_values": f(),
            "constraints": f(),
        })
    m["generic_template"] = {
        "fixed_structure": f(fixed, {"source_fixed_structure": fixed}, evidence),
        "variables": tv,
        "instance_mapping": [],
        "reason": f("R-group structure in the source", "Generic source drawing cannot be represented as one concrete molecule without inventing substituents.", evidence),
    }
    for k, value in [("chirality_status", "unresolved"), ("geometry_status", "unresolved"), ("representation_status", "not_encoded")]:
        m["stereochemistry"][k] = inferred(
            "generic R-group family", value, evidence,
            "Stereochemistry is not globally assignable for a variable family; no concrete stereoisomer is asserted.",
        )
    m["identity_note"] = f("generic source family", "Template retained without wildcard SMILES.", evidence)
    bundle["molecules"].append(m)
    molecules[mid] = m
    return mid


# Recurrent concrete reagents with conservative candidate graphs.
simple = {
    "2a": ("2a", "alpha-diazo methyl ketone 2a", "CC(=O)C(=[N+]=[N-])c1ccccc1", ev(2)),
    "2b": ("2b", "p-tolyl alpha-diazo methyl ketone 2b", "CC(=O)C(=[N+]=[N-])c1ccc(C)cc1", ev(3)),
    "2c": ("2c", "p-bromophenyl alpha-diazo methyl ketone 2c", "CC(=O)C(=[N+]=[N-])c1ccc(Br)cc1", ev(3)),
    "2d": ("2d", "p-chlorophenyl alpha-diazo methyl ketone 2d", "CC(=O)C(=[N+]=[N-])c1ccc(Cl)cc1", ev(3)),
    "2e": ("2e", "p-fluorophenyl alpha-diazo methyl ketone 2e", "CC(=O)C(=[N+]=[N-])c1ccc(F)cc1", ev(3)),
    "2f": ("2f", "m-tolyl alpha-diazo methyl ketone 2f", "CC(=O)C(=[N+]=[N-])c1cccc(C)c1", ev(3)),
    "2g": ("2g", "o-fluorophenyl alpha-diazo methyl ketone 2g", "CC(=O)C(=[N+]=[N-])c1ccccc1F", ev(3)),
    "2h": ("2h", "3,5-dimethylphenyl alpha-diazo methyl ketone 2h", "CC(=O)C(=[N+]=[N-])c1cc(C)cc(C)c1", ev(3)),
    "2i": ("2i", "2-naphthyl alpha-diazo methyl ketone 2i", "CC(=O)C(=[N+]=[N-])c1ccc2ccccc2c1", ev(3)),
    "2j": ("2j", "alpha-diazo isopropyl ketone 2j", "CC(C)C(=O)C(=[N+]=[N-])c1ccccc1", ev(3)),
    "2k": ("2k", "alpha-diazo diphenyl ketone 2k", "O=C(c1ccccc1)C(=[N+]=[N-])c1ccccc1", uniq(ev(3)+ev(4))),
    "2l": ("2l", "alpha-diazo p-tolyl phenyl ketone 2l", "O=C(c1ccc(C)cc1)C(=[N+]=[N-])c1ccccc1", uniq(ev(3)+ev(4))),
    "2m": ("2m", "alpha-diazo p-chlorophenyl phenyl ketone 2m", "O=C(c1ccc(Cl)cc1)C(=[N+]=[N-])c1ccccc1", uniq(ev(3)+ev(4))),
    "2n": ("2n", "alpha-diazo p-bromophenyl phenyl ketone 2n", "O=C(c1ccc(Br)cc1)C(=[N+]=[N-])c1ccccc1", uniq(ev(3)+ev(4))),
    "2o": ("2o", "p-chlorophenyl-substituted alpha-diazo phenyl ketone 2o", "O=C(c1ccccc1)C(=[N+]=[N-])c1ccc(Cl)cc1", ev(3)),
    "2p": ("2p", "p-bromophenyl-substituted alpha-diazo phenyl ketone 2p", "O=C(c1ccccc1)C(=[N+]=[N-])c1ccc(Br)cc1", ev(3)),
    "2q": ("2q", "alpha-diazo p-methoxyphenyl phenyl ketone 2q", "COc1ccc(C(=O)C(=[N+]=[N-])c2ccccc2)cc1", ev(4)),
    "2r": ("2r", "alpha-diazo m-tolyl phenyl ketone 2r", "Cc1cccc(C(=O)C(=[N+]=[N-])c2ccccc2)c1", ev(4)),
    "2s": ("2s", "alpha-diazo 2-thienyl phenyl ketone 2s", "O=C(c1cccs1)C(=[N+]=[N-])c1ccccc1", ev(4)),
    "diazo_ester": ("PhC(N2)CO2Et", "ethyl 2-diazo-2-phenylacetate", "CCOC(=O)C(=[N+]=[N-])c1ccccc1", ev(3)),
    "water": ("H2O", "water", "O", uniq(ev(3)+ev(4))),
    "h2o2": ("H2O2", "hydrogen peroxide", "OO", ev(3)),
    "h2": ("H2", "hydrogen", "[H][H]", ev(3)),
    "phenylboronic_acid": ("PhB(OH)2", "phenylboronic acid", "OB(O)c1ccccc1", ev(3)),
    "phenylacetylene": ("PhCCH", "phenylacetylene", "C#Cc1ccccc1", ev(3)),
    "mcpba": ("m-CPBA", "meta-chloroperoxybenzoic acid", "O=C(OO)c1cccc(Cl)c1", ev(3)),
    "bh3": ("BH3", "borane", "[BH3]", ev(3)),
    "dce": ("DCE", "1,2-dichloroethane", "ClCCCl", uniq(ev(2)+ev(3)+ev(4)+ev(6))),
    "dcm": ("DCM", "dichloromethane", "ClCCl", uniq(ev(2)+ev(3))),
    "thf": ("THF", "tetrahydrofuran", "C1CCOC1", uniq(ev(2)+ev(3))),
    "toluene": ("toluene", "toluene", "Cc1ccccc1", uniq(ev(2)+ev(3)+ev(4))),
    "ethanol": ("EtOH", "ethanol", "CCO", ev(3)),
    "dioxane": ("dioxane", "1,4-dioxane", "O1CCOCC1", ev(3)),
    "triethylamine": ("Et3N", "triethylamine", "CCN(CC)CC", ev(3)),
    "diketone": ("2k-O", "1,2-diphenylethane-1,2-dione", "O=C(c1ccccc1)C(=O)c1ccccc1", ev(4)),
}
mid = {}
for key, (label, name, smiles, evidence) in simple.items():
    mid[key] = add_molecule(key, label, evidence, name=name, smiles=smiles, kind="simple")

# Vinylallenes, cyclopentadienes, products, and compact catalysts are retained without unsafe SMILES.
for label in ["1a","1b","1c","1d","1e","1f","1g","1h","1i","1j","1k","1l","1m","1n","1o","1p"]:
    evidence = uniq((ev(2) if label != "1d" else ev(2)+ev(3)+ev(4)))
    mid[label] = add_molecule(label, label, evidence, name=f"vinylallene {label}", kind="allene",
                              chirality="unspecified", geometry="unspecified")
for label in ["8b","8f","8g","8h","8i","8j","8k"]:
    mid[label] = add_molecule(label, label, ev(4), name=f"vinylallene {label}", kind="allene")
mid["1a_prime"] = add_molecule("1a_prime", "1a'", uniq(ev(2)+ev(4)), name="cyclopentadiene 1a'",
                                chirality="unspecified", geometry="not_applicable")
mid["1d_prime"] = add_molecule("1d_prime", "1d'", ev(3), name="cyclopentadiene 1d'",
                                chirality="unspecified", geometry="not_applicable")

for label in ["3a","3a_prime"] + [f"3{x}" for x in "bcdefghijklmnop"] + [f"4{x}" for x in "bcdefghijklmnop"]:
    evidence = ev(2) if label.startswith("3") else ev(3)
    mid[label] = add_molecule(label, label.replace("_prime", "'"), evidence,
                              name=f"bicyclic 2,3-dihydrofuran {label.replace('_prime', chr(39))}",
                              chirality="unspecified", geometry="not_applicable")

# Distinct enantioenriched products are not merged with unassigned/racemic-condition products.
mid["6d_unassigned"] = add_molecule("6d_unassigned", "6d", ev(3), name="bicyclic furan 6d from achiral conditions",
                                     chirality="unspecified", geometry="not_applicable")
mid["6d_enantioenriched"] = add_molecule(
    "6d_enantioenriched", "(+)-6d", uniq(ev(3)+ev(4)), name="enantioenriched (+)-6d",
    chirality="partially_specified", geometry="not_applicable",
    assignment="Positive rotation and wedge/dash relative drawing reported; absolute configuration inferred by the authors from (+)-6q, without an R/S name on the supplied pages.",
    mixture="Enantiomeric mixture; experiment-specific (+)/(−) ratios are stored as reaction outcomes.",
)
for label in ["4k_plus","6a","6b","6c","6e","6f","6g","6h","6i","6j","6k","6l","6m","6n","6o","6p","6q"]:
    source_label = "(+)-4k" if label == "4k_plus" else f"(+)-{label}"
    assignment = "Absolute configuration characterized by X-ray diffraction." if label == "6q" else "Positive rotation and source wedge/dash drawing reported; no explicit R/S name transcribed."
    mid[label] = add_molecule(label, source_label, ev(4), name=f"enantioenriched product {source_label}",
                              chirality="partially_specified", geometry="not_applicable",
                              assignment=assignment,
                              mixture="Enantiomeric mixture; experiment-specific er is stored on the reaction record.")

for label, assignment, mixture in [
    ("5a", None, None), ("5b", None, None),
    ("5c", "Relative stereochemistry shown; absolute configuration not assigned.", "Diastereomeric mixture, dr >20:1."),
    ("5d", "Relative stereochemistry shown; absolute configuration not assigned.", "Diastereomeric mixture, dr >20:1."),
    ("5e", "Structure confirmed by X-ray diffraction; no R/S name transcribed.", None),
]:
    mid[label] = add_molecule(label, label, ev(3), name=f"functionalization product {label}",
                              chirality="partially_specified" if assignment else "unspecified",
                              geometry="not_applicable", assignment=assignment, mixture=mixture)

mid["eq2_printed_3a"] = add_molecule(
    "eq2_printed_3a", "3a (printed in Eq. 2)", ev(4),
    name="annulation product printed as 3a in control Equation 2",
    chirality="unspecified", geometry="not_applicable",
)
for label in ["9a","9b"]:
    mid[label] = add_molecule(label, label, ev(6), name=f"acyclic diene {label}",
                              chirality="unspecified", geometry="unspecified")
mid["9q"] = add_molecule("9q", "9q", ev(6), name="1,3-disubstituted cyclopentadiene 9q",
                          chirality="unspecified", geometry="not_applicable")
mid["10a"] = add_molecule("10a", "10a", ev(6), name="annulation product 10a",
                           chirality="unspecified", geometry="unspecified")
mid["complicated_mixture"] = add_molecule(
    "complicated_mixture", "complicated mixture", ev(6), name="several unidentified minor products",
    chirality="unresolved", geometry="unresolved", mixture="Composition and individual product identities are not reported."
)
mid["three_inseparable_products"] = add_molecule(
    "three_inseparable_products", "three inseparable products", ev(3),
    name="three unidentified inseparable products from the alpha-diazo ester",
    chirality="unresolved", geometry="unresolved",
    mixture="Three products were reported as inseparable; individual structures and composition are not given."
)

# Generic families on page 1.
mid["generic_alkene"] = add_template("generic_alkene", "mono-, di-, or trisubstituted alkene", ev(1),
                                      "generic substituted alkene", ["R1","R2","R3","R4"])
mid["generic_diazo"] = add_template("generic_diazo", "generic alpha-diazo carbonyl", ev(1),
                                     "N2=C(X)-C(=O)Y family", ["X","Y"])
mid["generic_cyclopropane"] = add_template("generic_cyclopropane", "cyclopropane I", ev(1),
                                            "generic cyclopropane carbonyl product", ["R1","R2","R3","X","Y"])
mid["generic_dhf"] = add_template("generic_dhf", "2,3-dihydrofuran II", ev(1),
                                   "generic 2,3-dihydrofuran family", ["R1","R2","R3","X","Y"])
mid["generic_cp"] = add_template("generic_cp", "substituted cyclopentadiene", ev(1),
                                  "cyclopentadiene bearing Ar1, R1 and R2", ["Ar1","R1","R2"])
mid["generic_product"] = add_template("generic_product", "enantioenriched bicyclic 2,3-dihydrofuran", ev(1),
                                       "bicyclic [3+2]-annulation product", ["Ar1","Ar2","Ar3","R1","R2"])

# Reagents/catalysts/materials represented by source names only.
agent_names = {
    "silica": "silica gel", "LAuCl": "LAuCl (L=P(t-Bu)2(o-biphenyl))", "AgOTf": "AgOTf",
    "PhO3PAuCl": "(PhO)3PAuCl", "PPh3AuCl": "PPh3AuCl", "IPrAuCl": "IPrAuCl",
    "AgNTf2": "AgNTf2", "AgSbF6": "AgSbF6", "NaBARF": "NaBARF",
    "dpppAuCl2": "dppp(AuCl)2", "Rh2OAc4": "Rh2(OAc)4", "CuOTf2": "Cu(OTf)2",
    "chiral_gold": "L(AuCl)2, L=(R)-SEGPHOS", "MS4A": "4 Å molecular sieves",
    "CPA7a": "chiral phosphoric acid 7a", "CPA7b": "chiral phosphoric acid 7b",
    "CPA7c": "chiral phosphoric acid 7c", "CPA7d": "chiral phosphoric acid 7d",
    "CPA7e": "chiral phosphoric acid 7e", "CPA7f": "chiral phosphoric acid 7f",
    "CPA7g": "chiral phosphoric acid 7g", "CPA7aAg": "chiral silver phosphate 7a-Ag",
    "PdPPh34": "Pd(PPh3)4", "BaOH2_8H2O": "Ba(OH)2·8H2O",
    "PdCl2PPh32": "PdCl2(PPh3)2", "CuI": "CuI", "PdC": "Pd/C",
    "IPrAuSbF6": "IPrAuSbF6", "Int1": "gold carbene Int-1",
}
for key, name in agent_names.items():
    pages = uniq(ev(2)+ev(3)+ev(4)+ev(6))
    mid[key] = add_molecule(key, name, pages, name=name, chirality="unresolved", geometry="unresolved")


def participant(r, key, role, evidence, *, stage="reaction", amount=None, equivalents=None,
                mol_percent=None, volume=None, mass=None, reference=None, pid=None):
    pid = pid or f"p{len(r['participants'])+1}_{key}"
    p = st.new_participant(pid, mid[key], role, evidence, stage=stage)
    if amount is not None:
        p["charges"]["amount"] = amount
    if equivalents is not None:
        p["charges"]["equivalents"] = equivalents
    if mol_percent is not None:
        p["charges"]["mol_percent"] = mol_percent
    if volume is not None:
        p["charges"]["volume"] = volume
    if mass is not None:
        p["charges"]["mass"] = mass
    if reference is not None:
        p["equivalence_reference"] = f(reference, mid[reference], evidence)
    r["participants"].append(p)
    return pid


def step(r, sid, evidence, participant_ids, *, operation=None, temperature=None, time=None,
         atmosphere=None, stage="reaction"):
    s = st.new_step(sid, evidence, stage=stage)
    s["participant_ids"] = list(participant_ids)
    if operation is not None:
        s["operation"] = f(operation, operation, evidence)
    if temperature is not None:
        s["temperature"] = temperature
    if time is not None:
        s["time"] = time
    if atmosphere is not None:
        s["atmosphere"] = f(atmosphere, atmosphere, evidence)
    r["steps"].append(s)


def outcome(r, oid, kind, targets, evidence, raw, normalized, *, basis=None, scope=None,
            yield_type=None):
    o = st.new_outcome(oid, kind, [mid[t] for t in targets], evidence)
    o["value"] = f(raw, normalized, evidence)
    if basis is not None:
        o["measurement_basis"] = f(basis, basis, evidence)
    if scope is not None:
        o["scope"] = f(scope, scope, evidence)
    if yield_type is not None:
        o["yield_type"] = f(yield_type, yield_type, evidence)
    r["outcomes"].append(o)


def add_note(r, text, evidence, *, method="original_extraction", derivation=None):
    r["notes"].append(f(text, text, evidence, method=method, derivation=derivation))


def set_class(r, evidence, changes=None, *, coverage="partial", limitation=None):
    c = r["classification"]
    if not changes:
        c["coverage"] = inferred("No observed classifiable transformation", coverage, evidence,
                                 "No bond-formation label is asserted for no-reaction, unresolved-mixture, literature-context, or untranscribed generic records.")
        c["labels"] = inferred("No labels asserted", [], evidence,
                               "Classification is not applicable or remains unresolved from the supplied outcome.")
        c["limitations"] = [limitation or "No observed product graph was safely mapped."]
        return
    vals = []
    labels = {}
    for i, (a, b, before, after, change_type) in enumerate(changes, 1):
        cid = f"b{i}"
        val = {"change_id": cid, "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
               "before_order": before, "after_order": after, "change_type": change_type}
        vals.append(val)
        if change_type in {"formed", "order_increased", "hydrogen_gain"}:
            label = st.bond_label(a, b, after)
            labels.setdefault(label["code"], {**label, "supporting_change_ids": []})["supporting_change_ids"].append(cid)
    c["bond_changes"] = [f("class-level source transformation", x, evidence, method="inference",
                                derivation="Class-level bond change inferred from the inspected source reactant/product drawing; exact atom identities and full mapping are unresolved.") for x in vals]
    c["labels"] = f("class-level source transformation", list(labels.values()), evidence,
                    method="inference", derivation="All and only formed/promoted class-level bonds in the retained partial change list are indexed.")
    c["coverage"] = inferred("partial class-level classification", coverage, evidence,
                             "No complete atom mapping was asserted; implicit hydrogen and coproduct accounting may remain incomplete.")
    c["limitations"] = [limitation or "Exact atom identities, full atom mapping, coproducts, and complete hydrogen accounting remain unresolved."]


def reaction(suffix, evidence, *, record_type="specific_experiment", categories=None,
             location=None, name=None, transformation=None, scheme_ids=None):
    rid = f"{DOC}:{suffix}"
    r = st.new_reaction(rid, DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories or [], categories or [], evidence)
    r["location"] = f(location or suffix, {"source": location or suffix}, evidence)
    if name:
        r["reaction_name"] = f(name, name, evidence)
    if transformation:
        r["transformation_type"] = f(transformation, transformation, evidence)
    r["general_scheme_ids"] = list(scheme_ids or [])
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    bundle["reactions"].append(r)
    return r


ANNULATION = [("C","C",0,1,"formed"), ("C","O",0,1,"formed")]

# Page 1: prior-work families and this-work generic reaction.
r = reaction("scheme1_prior_cyclopropanation", ev(1), record_type="literature_example",
             location="Page 1 Scheme 1, previous work, first arrow", name="metal-catalyzed cyclopropanation",
             transformation="generic alkene cyclopropanation")
participant(r, "generic_alkene", "reactant", ev(1)); participant(r, "generic_diazo", "reactant", ev(1)); participant(r, "generic_cyclopropane", "product", ev(1))
set_class(r, ev(1), None, coverage="unresolved", limitation="Generic cited family; exact structures and mappings are not supplied as one experiment.")

r = reaction("scheme1_prior_cloke_wilson", ev(1), record_type="literature_example",
             location="Page 1 Scheme 1, previous work, I to II", name="Cloke-Wilson rearrangement",
             transformation="cyclopropane to 2,3-dihydrofuran rearrangement")
participant(r, "generic_cyclopropane", "reactant", ev(1)); participant(r, "generic_dhf", "product", ev(1))
set_class(r, ev(1), None, coverage="unresolved", limitation="Cited generic rearrangement; no single concrete experiment is reported here.")

generic_id = f"{DOC}:scheme1_current_generic"
r = reaction("scheme1_current_generic", ev(1), record_type="generic_scheme", location="Page 1 Scheme 1, This work",
             name="gold-catalyzed [3+2]-annulation", transformation="enantioselective [3+2]-annulation")
participant(r, "generic_cp", "reactant", ev(1)); participant(r, "generic_diazo", "reactant", ev(1)); participant(r, "generic_product", "product", ev(1))
set_class(r, ev(1), ANNULATION)


def add_silica_purification(r, evidence):
    pid = participant(r, "silica", "purification_agent", evidence, stage="purification")
    step(r, f"s{len(r['steps'])+1}_purification", evidence, [pid], operation="purification from a silica column", stage="purification")


def add_table1(entry, catalyst_keys, minutes, solvent_key, vals, *, starting_prime=False, prose_lt2=False):
    evidence = ev(2)
    r = reaction(f"table1_entry{entry}", evidence, categories=["optimization"], location=f"Page 2 Table 1 entry {entry}",
                 name="one-pot catalytic reaction", transformation="gold-catalyzed [3+2]-annulation",
                 scheme_ids=[generic_id])
    start_key = "1a_prime" if starting_prime else "1a"
    ps = participant(r, start_key, "reactant", evidence)
    p2 = participant(r, "2a", "reactant", evidence, equivalents=fq(2.0,"equiv",evidence,"2.0 equiv"), reference=start_key)
    catpids = []
    for i, ck in enumerate(catalyst_keys):
        catpids.append(participant(r, ck, "catalyst" if i == 0 else "additive", evidence,
                                   mol_percent=fq(5,"mol%",evidence,"5 mol% catalyst") if i == 0 else None,
                                   reference=start_key if i == 0 else None))
    solv = participant(r, solvent_key, "solvent", evidence)
    pids = [ps,p2,*catpids,solv]
    if not starting_prime:
        step(r,"s1_cyclization",evidence,[ps,*catpids,solv],operation="treat vinylallene with catalyst before diazo addition",
             temperature=fq(25,"°C",evidence,"25 °C"),time=fq(3,"min",evidence,"3 min"))
        r["solvent_conditions"]["concentration"] = fq(0.15,"M",evidence,"[1a] = 0.15 M")
        r["solvent_conditions"]["concentration_basis"] = f("[1a]",mid["1a"],evidence)
    else:
        add_note(r,"Entries 13–15 use cyclopentadiene 1a' as the reactant; the table footnote reports [1a] = 0.15 M, so no 1a' concentration is assigned.",evidence)
        r["review"]["unresolved"].append("Starting cyclopentadiene 1a' concentration is not explicitly reported; the footnote names 1a instead.")
    step(r,"s2_annulation",evidence,pids,operation="add diazo species 2a and continue the one-pot reaction",
         time=fq(minutes,"min",evidence,f"{minutes} min"))
    r["solvent_conditions"]["composition"] = f(agent_names.get(solvent_key, simple.get(solvent_key,(solvent_key,))[0]), mid[solvent_key], evidence)
    basis = "after purification from a silica column"
    if vals.get("1a_prime") is not None:
        value = vals["1a_prime"]
        if starting_prime:
            outcome(r,"recovery_1a_prime","recovery",["1a_prime"],evidence,f"{value}%",{"number":value,"unit":"%"},basis=basis,scope="recovered starting cyclopentadiene 1a'",yield_type="recovery")
        else:
            if value > 0 or entry == 8:
                if not any(p["molecule_id"] == mid["1a_prime"] for p in r["participants"]): participant(r,"1a_prime","product",evidence)
                outcome(r,"yield_1a_prime","yield",["1a_prime"],evidence,f"{value}%",{"number":value,"unit":"%"},basis=basis,scope="cyclopentadiene intermediate isolated after run",yield_type="product yield")
    if vals.get("3a") is not None:
        participant(r,"3a","product",evidence)
        outcome(r,"yield_3a","yield",["3a"],evidence,f"{vals['3a']}%",{"number":vals["3a"],"unit":"%"},basis=basis,yield_type="product yield")
    if vals.get("3a_prime") is not None:
        participant(r,"3a_prime","product",evidence)
        outcome(r,"yield_3a_prime","yield",["3a_prime"],evidence,f"{vals['3a_prime']}%",{"number":vals["3a_prime"],"unit":"%"},basis=basis,yield_type="product yield")
    if prose_lt2:
        if not any(p["molecule_id"] == mid["3a_prime"] for p in r["participants"]): participant(r,"3a_prime","product",evidence)
        outcome(r,"yield_3a_prime_lt2","yield",["3a_prime"],evidence,"negligible yields (<2%)",{"number":2,"unit":"%","qualifier":"lt","raw_text":"<2%"},basis=basis,yield_type="product yield")
        add_note(r,"Table cell for 3a' is a dash, while adjacent prose reports a negligible yield below 2%; both source expressions are retained.",evidence)
    dash_cols = [k for k in ("1a_prime","3a","3a_prime") if k not in vals]
    if dash_cols:
        add_note(r,"Table 1 dash (–) in columns: "+", ".join(dash_cols),evidence)
    if entry in (8,14):
        outcome(r,"no_desired_annulation_product","no_reaction",[],evidence,"failed to give noticeable annulation product","failed to give noticeable annulation product",scope="desired annulation product")
    add_silica_purification(r,evidence)
    set_class(r,evidence,ANNULATION if vals.get("3a",0)>0 or vals.get("3a_prime",0)>0 else None,
              coverage="not_applicable" if entry in (8,14) else "partial")
    return r


table1 = [
    (1,["LAuCl","AgOTf"],20,"dce",{"1a_prime":28,"3a":36}),
    (2,["PhO3PAuCl","AgOTf"],20,"dce",{"3a":37,"3a_prime":18}),
    (3,["PPh3AuCl","AgOTf"],30,"dce",{"3a":48,"3a_prime":24}),
    (4,["IPrAuCl","AgOTf"],10,"dce",{"1a_prime":5,"3a":68}),
    (5,["IPrAuCl","AgNTf2"],10,"dce",{"3a":74}),
    (6,["IPrAuCl","AgSbF6"],10,"dce",{"3a":78}),
    (7,["IPrAuCl","NaBARF"],60,"dce",{"3a":71}),
    (8,["AgSbF6"],60,"dce",{"1a_prime":0}),
    (9,["IPrAuCl","AgSbF6"],15,"dcm",{"3a":56}),
    (10,["IPrAuCl","AgSbF6"],20,"thf",{"1a_prime":27,"3a":41}),
    (11,["IPrAuCl","AgSbF6"],30,"toluene",{"1a_prime":35,"3a":29}),
    (12,["dpppAuCl2","AgSbF6"],60,"dce",{"3a":68,"3a_prime":6}),
    (13,["IPrAuCl","AgSbF6"],20,"dce",{"3a":75}),
    (14,["Rh2OAc4"],30,"dce",{"1a_prime":36}),
    (15,["CuOTf2"],180,"dce",{"1a_prime":5,"3a":16,"3a_prime":8}),
]
for row in table1:
    entry = row[0]
    add_table1(*row, starting_prime=entry >= 13, prose_lt2=entry == 4)


def add_scope_annulation(table, entry, substrate, diazo, product, yield_pct, page, *, concentration, pre_time=3,
                         catalyst="IPrAuCl", activator="AgSbF6"):
    evidence = ev(page)
    r = reaction(f"{table}_entry{entry}_{product}",evidence,categories=["substrate_scope"],
                 location=f"Page {page} {table.replace('_',' ').title()} entry {entry}",name="[3+2]-annulation",
                 transformation="gold-catalyzed [3+2]-annulation",scheme_ids=[generic_id])
    ps=participant(r,substrate,"reactant",evidence)
    pd=participant(r,diazo,"reactant",evidence,equivalents=fq(2.0,"equiv",evidence,"2.0 equiv"),reference=substrate)
    pc=participant(r,catalyst,"catalyst",evidence,mol_percent=fq(5,"mol%",evidence,"5 mol%"),reference=substrate)
    pa=participant(r,activator,"additive",evidence,mol_percent=fq(5,"mol%",evidence,"5 mol%"),reference=substrate)
    psol=participant(r,"dce","solvent",evidence)
    pp=participant(r,product,"product",evidence)
    step(r,"s1_cyclization",evidence,[ps,pc,pa,psol],operation="gold-catalyzed formation of substituted cyclopentadiene",
         temperature=f("rt","rt",evidence),time=fq(pre_time,"min",evidence,f"{pre_time} min"))
    step(r,"s2_annulation",evidence,[ps,pd,pc,pa,psol,pp],operation="add diazo ketone and continue one pot")
    r["solvent_conditions"]["composition"] = f("DCE",mid["dce"],evidence)
    r["solvent_conditions"]["concentration"] = fq(concentration,"M",evidence,f"[{substrate[0]}] = {concentration} M")
    r["solvent_conditions"]["concentration_basis"] = f(substrate,mid[substrate],evidence)
    outcome(r,"yield","yield",[product],evidence,f"{yield_pct}%",{"number":yield_pct,"unit":"%"},basis="after purification from a silica column",yield_type="product yield")
    add_silica_purification(r,evidence)
    set_class(r,evidence,ANNULATION)
    r["review"]["unresolved"].append("The time and temperature after diazo-ketone addition are not separately reported in the table header.")
    return r


table2 = [
    (1,"1b","3b",52),(2,"1c","3c",65),(3,"1d","3d",62),(4,"1e","3e",67),(5,"1f","3f",59),
    (6,"1g","3g",65),(7,"1h","3h",63),(8,"1i","3i",62),(9,"1j","3j",68),(10,"1k","3k",59),
    (11,"1l","3l",62),(12,"1m","3m",42),(13,"1n","3n",48),(14,"1o","3o",46),(15,"1p","3p",49),
]
for e,s,p,y in table2: add_scope_annulation("table2",e,s,"2a",p,y,2,concentration=0.15)

table3 = [
    (1,"2b","4b",71),(2,"2c","4c",72),(3,"2d","4d",67),(4,"2e","4e",66),(5,"2f","4f",66),
    (6,"2g","4g",38),(7,"2h","4h",68),(8,"2i","4i",71),(9,"2j","4j",63),(10,"2k","4k",75),
    (11,"2l","4l",51),(12,"2m","4m",67),(13,"2n","4n",62),(14,"2o","4o",68),(15,"2p","4p",70),
]
for e,d,p,y in table3: add_scope_annulation("table3",e,"1a",d,p,y,3,concentration=0.15)

# Prose-only Table 3 boundary substrate: a real experiment without identifiable products/yields.
r=reaction("table3_diazo_ester_inapplicable",ev(3),categories=["substrate_scope"],
           location="Page 3 prose below Table 3, alpha-diazo ester",
           name="attempted [3+2]-annulation with an alpha-diazo ester",
           transformation="unresolved reaction mixture",scheme_ids=[generic_id])
p1=participant(r,"1a","reactant",ev(3)); p2=participant(r,"diazo_ester","reactant",ev(3),equivalents=fq(2.0,"equiv",ev(3),"2.0 equiv"),reference="1a")
pc=participant(r,"IPrAuCl","catalyst",ev(3),mol_percent=fq(5,"mol%",ev(3),"5 mol%"),reference="1a"); pa=participant(r,"AgSbF6","additive",ev(3),mol_percent=fq(5,"mol%",ev(3),"5 mol%"),reference="1a")
ps=participant(r,"dce","solvent",ev(3)); pp=participant(r,"three_inseparable_products","product",ev(3))
step(r,"s1_cyclization",ev(3),[p1,pc,pa,ps],operation="gold-catalyzed formation of cyclopentadiene 1a'",temperature=f("rt","rt",ev(3)),time=fq(3,"min",ev(3),"3 min"))
step(r,"s2_attempt",ev(3),[p1,p2,pc,pa,ps,pp],operation="add alpha-diazo ester under Table 3 optimized conditions; post-addition duration is not reported")
r["solvent_conditions"]["composition"] = f("DCE",mid["dce"],ev(3)); r["solvent_conditions"]["concentration"] = fq(0.15,"M",ev(3),"[1a] = 0.15 M"); r["solvent_conditions"]["concentration_basis"] = f("1a",mid["1a"],ev(3))
outcome(r,"inseparable_products","other",["three_inseparable_products"],ev(3),"inapplicable substrate; three inseparable products","three inseparable products",scope="unidentified mixture")
set_class(r,ev(3),None,coverage="unresolved")
r["review"]["unresolved"] += ["Individual product identities, yields, conversion, and measurement methods are unreported.","Post-diazo-addition duration and temperature are not separately reported."]

# Page 3 Scheme 2: scale-up plus five independent functionalizations.
r = reaction("scheme2_gram_scale_6d",ev(3),categories=["scale_up"],location="Page 3 Scheme 2, gram-scale 6d",
             name="gram-scale [3+2]-annulation",transformation="gold-catalyzed [3+2]-annulation",scheme_ids=[generic_id])
p1=participant(r,"1d","reactant",ev(3)); p2=participant(r,"2k","reactant",ev(3),equivalents=fq(2.0,"equiv",ev(3),"2.0 equiv"),reference="1d")
pc=participant(r,"IPrAuCl","catalyst",ev(3),mol_percent=fq(5,"mol%",ev(3),"5 mol%"),reference="1d")
pa=participant(r,"AgSbF6","additive",ev(3),mol_percent=fq(5,"mol%",ev(3),"5 mol%"),reference="1d")
ps=participant(r,"dce","solvent",ev(3)); pp=participant(r,"6d_unassigned","product",ev(3))
step(r,"s1_reaction",ev(3),[p1,p2,pc,pa,ps,pp],operation="gram-scale one-pot annulation",temperature=f("rt","rt",ev(3)),time=fq(30,"min",ev(3),"30 min"))
outcome(r,"yield","yield",["6d_unassigned"],ev(3),"70% gram scale",{"number":70,"unit":"%"},scope="1d to 6d gram-scale reaction")
outcome(r,"product_mass","product_amount",["6d_unassigned"],ev(3),"1.05 g",{"number":1.05,"unit":"g"})
set_class(r,ev(3),ANNULATION)
r["review"]["unresolved"].append("Yield measurement basis is not stated in Scheme 2.")

def deriv(suffix, product, reactants, conditions, yield_pct, *, changes=None, selectivity=None, review=None):
    evidence=ev(3)
    r=reaction(suffix,evidence,categories=["derivatization"],location=f"Page 3 Scheme 2, {suffix}",name="chemical functionalization of 6d",
               transformation=conditions["transformation"])
    pids=[]
    for key,role in reactants:
        kwargs={}
        if key=="mcpba": kwargs={"equivalents":fq(2.2,"equiv",evidence,"2.2 equiv"),"reference":"5d"}
        pids.append(participant(r,key,role,evidence,**kwargs))
    pids.append(participant(r,product,"product",evidence))
    for i,spec in enumerate(conditions["steps"],1):
        step(r,f"s{i}_{spec['id']}",evidence,pids,operation=spec.get("operation"),temperature=spec.get("temperature"),time=spec.get("time"))
    r["solvent_conditions"]["composition"] = f(conditions.get("solvent_raw"),conditions.get("solvent_norm"),evidence) if conditions.get("solvent_raw") else f()
    outcome(r,"yield","yield",[product],evidence,f"{yield_pct}%",{"number":yield_pct,"unit":"%"},yield_type="product yield")
    if selectivity:
        outcome(r,"selectivity","selectivity",[product],evidence,selectivity,{"metric":"dr","qualifier":"gt","ratio":[20,1]},scope={"metric":"dr","order":"unassigned major/minor"})
    set_class(r,evidence,changes,coverage="unresolved" if changes is None else "partial")
    r["review"]["unresolved"].append("Yield measurement basis is not stated in Scheme 2.")
    if selectivity:r["review"]["unresolved"].append("Diastereomer identities/order and selectivity measurement method are not reported.")
    if review:r["review"]["unresolved"].append(review)
    return r

deriv("scheme2_suzuki_5a","5a",[("6d_unassigned","reactant"),("phenylboronic_acid","reactant"),("PdPPh34","catalyst"),("BaOH2_8H2O","reagent"),("dioxane","solvent"),("water","solvent")],
      {"transformation":"Suzuki coupling","solvent_raw":"dioxane/H2O=1:1","solvent_norm":{"components":[mid['dioxane'],mid['water']],"ratio":"1:1"},"steps":[{"id":"coupling","operation":"Suzuki coupling","temperature":fq(90,"°C",ev(3),"90 °C"),"time":fq(12,"h",ev(3),"12 h")}]},74,
      changes=[("C","C",0,1,"formed"),("C","Br",1,0,"broken"),("C","B",1,0,"broken")])
deriv("scheme2_sonogashira_5b","5b",[("6d_unassigned","reactant"),("phenylacetylene","reactant"),("PdCl2PPh32","catalyst"),("CuI","catalyst"),("triethylamine","reagent")],
      {"transformation":"Sonogashira coupling","steps":[{"id":"coupling","operation":"Sonogashira coupling","temperature":fq(60,"°C",ev(3),"60 °C"),"time":fq(3,"h",ev(3),"3 h")}]},86,
      changes=[("C","C",0,1,"formed"),("C","Br",1,0,"broken"),("C","H",1,0,"hydrogen_loss")])
deriv("scheme2_hydroboration_5c","5c",[("6d_unassigned","reactant"),("bh3","reactant"),("h2o2","reagent"),("thf","solvent")],
      {"transformation":"hydroboration followed by oxidation","solvent_raw":"THF","solvent_norm":mid["thf"],"steps":[{"id":"hydroboration_oxidation","operation":"BH3-THF hydroboration followed by H2O2 oxidation","temperature":fq(60,"°C",ev(3),"60 °C"),"time":fq(3,"h",ev(3),"3 h")}]},72,
      changes=[("C","C",2,1,"order_decreased"),("C","O",0,1,"formed"),("C","H",0,1,"hydrogen_gain")],selectivity="dr >20:1")
deriv("scheme2_hydrogenation_5d","5d",[("6d_unassigned","reactant"),("h2","reactant"),("PdC","catalyst"),("ethanol","solvent"),("dcm","solvent")],
      {"transformation":"selective hydrogenation","solvent_raw":"EtOH/DCM=1:1","solvent_norm":{"components":[mid['ethanol'],mid['dcm']],"ratio":"1:1"},"steps":[{"id":"hydrogenation","operation":"selective hydrogenation","temperature":fq(0,"°C",ev(3),"0 °C"),"time":fq(3,"h",ev(3),"3 h")}]},88,
      changes=[("C","C",2,1,"order_decreased"),("C","H",0,1,"hydrogen_gain")],selectivity="dr >20:1")
deriv("scheme2_oxidative_cleavage_5e","5e",[("5d","reactant"),("mcpba","reactant"),("dcm","solvent")],
      {"transformation":"oxidative cleavage","solvent_raw":"DCM","solvent_norm":mid["dcm"],"steps":[
          {"id":"cold_addition","operation":"begin oxidative cleavage at 0 °C; addition/hold time not reported","temperature":fq(0,"°C",ev(3),"0 °C")},
          {"id":"warm_to_rt","operation":"warm to room temperature; warming duration is not separately reported","temperature":f("rt","rt",ev(3)),"time":fq(2,"h",ev(3),"2 h")},
      ]},73,
      changes=None,review="Exact net bond changes for the oxidative cleavage were not asserted without a complete safe graph mapping.")


def add_chiral_annulation(suffix, entry, substrate, diazo, product, yld, er, *, table4=False, gold=2.2, nabarf=2.0, cpa="CPA7d", cpa_load=4.0, starting_prime=False):
    page=3 if table4 else 4; evidence=ev(page)
    table="Table 4" if table4 else "Table 5"
    r=reaction(suffix,evidence,categories=["optimization"] if table4 else ["substrate_scope"],location=f"Page {page} {table} entry {entry}",
               name="one-pot chiral catalytic [3+2]-annulation",transformation="enantioselective [3+2]-annulation",scheme_ids=[generic_id])
    ps=participant(r,substrate,"reactant",evidence)
    pd=participant(r,diazo,"reactant",evidence,equivalents=fq(1.2,"equiv",evidence,"1.2 equiv"),reference=substrate)
    pg=participant(r,"chiral_gold","catalyst",evidence,mol_percent=fq(gold,"mol%",evidence,f"{gold} mol%"),reference=substrate)
    pids=[ps,pd,pg]
    if nabarf is not None:pids.append(participant(r,"NaBARF","additive",evidence,mol_percent=fq(nabarf,"mol%",evidence,f"{nabarf} mol%"),reference=substrate))
    if cpa is not None:pids.append(participant(r,cpa,"catalyst",evidence,mol_percent=fq(cpa_load,"mol%",evidence,f"{cpa_load} mol%"),reference=substrate))
    pids.append(participant(r,"MS4A","additive",evidence)); pids.append(participant(r,"toluene","solvent",evidence)); pids.append(participant(r,product,"product",evidence))
    if not starting_prime:
        step(r,"s1_cyclization",evidence,pids,operation="gold-catalyzed vinylallene cyclization in the presence of 4 Å molecular sieves",temperature=f("rt","rt",evidence),time=fq(3,"min",evidence,"3 min"))
        r["solvent_conditions"]["concentration"] = fq(0.125,"M",evidence,f"[{substrate}] = 0.125 M")
        r["solvent_conditions"]["concentration_basis"] = f(substrate,mid[substrate],evidence)
    else:
        add_note(r,"Table 4 entry 11 uses cyclopentadiene 1d' as starting material; footnote [1d] = 0.125 M is not assigned to 1d'.",evidence)
        r["review"]["unresolved"].append("Concentration of starting 1d' is not explicitly reported.")
    step(r,"s2_cpa_complexation",evidence,pids,operation="cool and hold with the listed CPA or phosphate before diazo addition",temperature=fq(0,"°C",evidence,"under 0 °C"),time=fq(1,"min",evidence,"1 min"))
    step(r,"s3_cold_annulation",evidence,pids,operation="add diazo ketone and hold at 0 °C",temperature=fq(0,"°C",evidence,"0 °C"),time=fq(3,"h",evidence,"3 h"))
    step(r,"s4_room_temperature_hold",evidence,pids,operation="warm to room temperature and hold; warming duration not reported",temperature=f("rt","rt",evidence),time=fq(1,"h",evidence,"1 h"))
    r["solvent_conditions"]["composition"] = f("toluene",mid["toluene"],evidence)
    outcome(r,"yield","yield",[product],evidence,f"{yld}%",{"number":yld,"unit":"%"},basis="after purification from a silica column",yield_type="product yield")
    a,b=map(int,er.split(":"))
    outcome(r,"er","selectivity",[product],evidence,er,{"metric":"er","ratio":[a,b],"order":["(+)","(-)"]},basis="HPLC analysis on a chiral stationary phase (see Supporting Information)",scope={"metric":"er","order":["(+)","(-)"]})
    add_silica_purification(r,evidence); set_class(r,evidence,ANNULATION)
    r["review"]["unresolved"].append("The warming transition time from 0 °C to room temperature is not reported.")
    return r


table4 = [
    (1,5.5,5.0,None,None,78,"78:22"),(2,5.5,5.0,"CPA7a",5.0,75,"87:13"),(3,2.2,2.0,"CPA7a",2.0,72,"87:13"),
    (4,2.2,2.0,"CPA7a",4.0,70,"88:12"),(5,2.2,2.0,"CPA7b",4.0,70,"88:12"),(6,2.2,2.0,"CPA7c",4.0,74,"91:9"),
    (7,2.2,2.0,"CPA7d",4.0,73,"93:7"),(8,2.2,2.0,"CPA7e",4.0,75,"89:11"),(9,2.2,2.0,"CPA7f",4.0,72,"89:11"),
    (10,2.2,2.0,"CPA7g",4.0,68,"87:13"),(11,2.2,None,"CPA7aAg",2.0,21,"82:18"),
]
for entry,gold,nb,cpa,cload,y,eratio in table4:
    add_chiral_annulation(f"table4_entry{entry}",entry,"1d_prime" if entry==11 else "1d","2k","6d_enantioenriched",y,eratio,
                          table4=True,gold=gold,nabarf=nb,cpa=cpa,cpa_load=cload or 0,starting_prime=entry==11)

def table4_qualitative_attempt(suffix, *, cpa_loading, modified_temperature=False):
    evidence=ev(4)
    r=reaction(suffix,evidence,categories=["optimization"],location="Page 4 prose following Table 4",
               name="qualitative modification of the optimized chiral annulation",
               transformation="enantioselective [3+2]-annulation",scheme_ids=[generic_id])
    p1=participant(r,"1d","reactant",evidence); p2=participant(r,"2k","reactant",evidence,equivalents=fq(1.2,"equiv",evidence,"1.2 equiv"),reference="1d")
    pg=participant(r,"chiral_gold","catalyst",evidence,mol_percent=fq(2.2,"mol%",evidence,"2.2 mol%"),reference="1d")
    pn=participant(r,"NaBARF","additive",evidence,mol_percent=fq(2.0,"mol%",evidence,"2.0 mol%"),reference="1d")
    pc=participant(r,"CPA7d","catalyst",evidence,mol_percent=fq(cpa_loading,"mol%",evidence,f"{cpa_loading} mol%"),reference="1d")
    pm=participant(r,"MS4A","additive",evidence); ps=participant(r,"toluene","solvent",evidence); pp=participant(r,"6d_enantioenriched","product",evidence)
    pids=[p1,p2,pg,pn,pc,pm,ps,pp]
    step(r,"s1_cyclization",evidence,pids,operation="optimized room-temperature vinylallene cyclization",temperature=f("rt","rt",evidence),time=fq(3,"min",evidence,"3 min"))
    if modified_temperature:
        step(r,"s2_lower_temperature_attempt",evidence,pids,operation="repeat annulation at a decreased reaction temperature; exact temperature and duration are not reported",temperature=f("a decrease in reaction temperature","lower temperature (exact value unreported)",evidence))
        r["review"]["unresolved"].append("Exact decreased temperature, time profile, numerical yield, er, and yield method are not reported.")
    else:
        step(r,"s2_cpa_hold",evidence,pids,operation="hold with 6 mol% CPA 7d before diazo addition",temperature=fq(0,"°C",evidence,"under 0 °C"),time=fq(1,"min",evidence,"1 min"))
        step(r,"s3_cold_annulation",evidence,pids,operation="annulation at 0 °C",temperature=fq(0,"°C",evidence,"0 °C"),time=fq(3,"h",evidence,"3 h"))
        step(r,"s4_rt_hold",evidence,pids,operation="warm to room temperature and hold; warming duration unreported",temperature=f("rt","rt",evidence),time=fq(1,"h",evidence,"1 h"))
        r["review"]["unresolved"].append("Numerical decreased yield, er, yield method, and warming transition time are not reported.")
    r["solvent_conditions"]["composition"] = f("toluene",mid["toluene"],evidence); r["solvent_conditions"]["concentration"] = fq(0.125,"M",evidence,"[1d] = 0.125 M"); r["solvent_conditions"]["concentration_basis"] = f("1d",mid["1d"],evidence)
    outcome(r,"decreased_yield","other",["6d_enantioenriched"],evidence,"led to decreased yields of (+)-6d","decreased yield",scope="qualitative product-yield comparison")
    set_class(r,evidence,ANNULATION)
    return r

table4_qualitative_attempt("table4_prose_lower_temperature",cpa_loading=4.0,modified_temperature=True)
table4_qualitative_attempt("table4_prose_cpa7d_6mol",cpa_loading=6.0)

table5 = [
    (1,"1a","2k","4k_plus",78,"85:15"),(2,"1b","2k","6a",64,"77:23"),(3,"8b","2k","6b",78,"89:11"),
    (4,"1e","2k","6c",75,"91:9"),(5,"1f","2k","6e",66,"96:4"),(6,"8f","2k","6f",62,"91:9"),
    (7,"8g","2k","6g",71,"91:9"),(8,"8h","2k","6h",65,"93:7"),(9,"8i","2k","6i",55,"93:7"),
    (10,"8j","2k","6j",51,"96:4"),(11,"8k","2k","6k",64,"93:7"),(12,"1d","2l","6l",72,"91:9"),
    (13,"1d","2m","6m",76,"87:13"),(14,"1d","2n","6n",71,"87:13"),(15,"1d","2q","6o",54,"87:13"),
    (16,"1d","2r","6p",62,"88:12"),(17,"1d","2s","6q",70,"92:8"),
]
for entry,s,d,p,y,eratio in table5:
    add_chiral_annulation(f"table5_entry{entry}_{p}",entry,s,d,p,y,eratio)

# Page 4 water-trapping mechanistic controls.
r=reaction("eq1_water_trap",ev(4),categories=["mechanistic_study","control"],location="Page 4 Equation (1)",name="water trapping of gold carbene",
           transformation="oxidation of alpha-diazo ketone to 1,2-diketone")
p2=participant(r,"2k","reactant",ev(4)); pw=participant(r,"water","reactant",ev(4),equivalents=fq(5,"equiv",ev(4),"5 equiv"),reference="2k")
pc=participant(r,"IPrAuSbF6","catalyst",ev(4),mol_percent=fq(5,"mol%",ev(4),"5 mol%"),reference="2k"); ps=participant(r,"dce","solvent",ev(4)); pp=participant(r,"diketone","product",ev(4))
step(r,"s1_water_trap",ev(4),[p2,pw,pc,ps,pp],operation="treat 2k with water under IPrAuSbF6 catalysis",temperature=fq(25,"°C",ev(4),"25 °C"),time=fq(3,"min",ev(4),"3 min"))
outcome(r,"yield","yield",["diketone"],ev(4),"74%",{"number":74,"unit":"%"},yield_type="product yield")
add_note(r,"Gold carbene Int-1 is drawn as a proposed intermediate and was not treated as an isolated product.",ev(4))
set_class(r,ev(4),[("C","O",0,2,"formed"),("C","N",2,0,"broken")])
r["review"]["unresolved"].append("Yield measurement basis is not stated for Equation (1).")

r=reaction("eq2_water_competition",ev(4),categories=["mechanistic_study","control"],location="Page 4 Equation (2)",name="water competition during [3+2]-annulation",
           transformation="competing annulation and diketone formation")
p1=participant(r,"1a_prime","reactant",ev(4)); p2=participant(r,"2k","reactant",ev(4),equivalents=fq(1.1,"equiv",ev(4),"1.1 equiv"),reference="1a_prime")
pw=participant(r,"water","reactant",ev(4),equivalents=fq(5,"equiv",ev(4),"5 equiv"),reference="1a_prime"); pc=participant(r,"IPrAuSbF6","catalyst",ev(4),mol_percent=fq(5,"mol%",ev(4),"5 mol%"),reference="1a_prime")
ps=participant(r,"dce","solvent",ev(4)); ppa=participant(r,"eq2_printed_3a","product",ev(4)); ppd=participant(r,"diketone","product",ev(4))
step(r,"s1_competition",ev(4),[p1,p2,pw,pc,ps,ppa,ppd],operation="run annulation in the presence of water",temperature=fq(25,"°C",ev(4),"25 °C"),time=fq(3,"min",ev(4),"3 min"))
outcome(r,"yield_printed_3a","yield",["eq2_printed_3a"],ev(4),"3a (8%)",{"number":8,"unit":"%"},yield_type="product yield")
outcome(r,"yield_diketone","yield",["diketone"],ev(4),"2k-O (62%)",{"number":62,"unit":"%"},yield_type="product yield")
add_note(r,"Equation (2) prints the annulation product label as 3a although the diazo substrate is 2k; this identity is kept separate from Table 1 compound 3a.",ev(4))
r["review"]["unresolved"] += ["Printed product label 3a in Equation (2) conflicts with the 2k substrate/other numbering; identity not silently merged.","Yield measurement bases are not stated for Equation (2)."]
set_class(r,ev(4),ANNULATION+[("C","O",0,2,"formed"),("C","N",2,0,"broken")])

# Scheme 3 is a proposal, not an executed series of intermediates.
r=reaction("scheme3_proposed_mechanism",uniq(ev(4)+ev(5)),record_type="proposed_mechanism",location="Pages 4–5 Scheme 3 and discussion",
           name="DFT-supported proposed [3+2]-annulation mechanism",transformation="proposed gold-carbene/gold-enolate annulation")
participant(r,"generic_cp","reactant",ev(4)); participant(r,"generic_diazo","reactant",ev(4)); participant(r,"generic_product","product",ev(4))
add_note(r,"Author-proposed sequence: cyclopentadiene attack on gold carbene, gold-enolate/allylic-cation intermediates, phosphoric-acid ion pairing, and ring closure; DFT energies are mechanistic context, not isolated experiments.",uniq(ev(4)+ev(5)))
set_class(r,uniq(ev(4)+ev(5)),ANNULATION)
r["review"]["unresolved"].append("Proposed intermediates A–E/Int-1 are not experimentally isolated on the supplied pages.")

# Page 6 footnote 13: three separately reported boundary-substrate experiments.
r=reaction("note13_diene9a",ev(6),categories=["substrate_scope"],location="Page 6 footnote 13, diene 9a",name="annulation of acyclic trisubstituted diene",
           transformation="[3+2]-annulation")
p1=participant(r,"9a","reactant",ev(6)); p2=participant(r,"2a","reactant",ev(6),equivalents=fq(2.0,"equiv",ev(6),"2.0 equiv"),reference="9a")
pc=participant(r,"IPrAuCl","catalyst",ev(6),mol_percent=fq(5,"mol%",ev(6),"5 mol%"),reference="9a"); pa=participant(r,"AgSbF6","additive",ev(6),mol_percent=fq(5,"mol%",ev(6),"5 mol%"),reference="9a"); ps=participant(r,"dce","solvent",ev(6)); pp=participant(r,"10a","product",ev(6))
step(r,"s1_annulation",ev(6),[p1,p2,pc,pa,ps,pp],operation="annulation of diene 9a with diazo ketone 2a",temperature=f("rt","rt",ev(6)),time=fq(2,"h",ev(6),"2 h"))
outcome(r,"yield","yield",["10a"],ev(6),"25%",{"number":25,"unit":"%"},yield_type="product yield"); set_class(r,ev(6),ANNULATION)
r["review"]["unresolved"].append("Yield measurement basis is not stated in footnote 13.")

r=reaction("note13_diene9b_recovery",ev(6),categories=["substrate_scope","control"],location="Page 6 footnote 13, diene 9b",name="unsuccessful annulation of tetrasubstituted acyclic diene",
           transformation="no applicable annulation")
p1=participant(r,"9b","reactant",ev(6)); p2=participant(r,"2a","reactant",ev(6),equivalents=fq(2.0,"equiv",ev(6),"2.0 equiv"),reference="9b")
pc=participant(r,"IPrAuCl","catalyst",ev(6),mol_percent=fq(5,"mol%",ev(6),"5 mol%"),reference="9b"); pa=participant(r,"AgSbF6","additive",ev(6),mol_percent=fq(5,"mol%",ev(6),"5 mol%"),reference="9b"); ps=participant(r,"dce","solvent",ev(6))
step(r,"s1_attempt",ev(6),[p1,p2,pc,pa,ps],operation="attempted annulation; source states reaction was inapplicable",temperature=f("rt","rt",ev(6)),time=fq(2,"h",ev(6),"2 h"))
outcome(r,"recovery","recovery",["9b"],ev(6),"82% recovery",{"number":82,"unit":"%"},yield_type="starting-material recovery")
outcome(r,"no_annulation","no_reaction",[],ev(6),"annulation was inapplicable","annulation was inapplicable",scope="desired annulation")
set_class(r,ev(6),None,coverage="not_applicable"); r["review"]["unresolved"].append("Recovery measurement basis is not stated in footnote 13.")

r=reaction("note13_cyclopentadiene9q",ev(6),categories=["substrate_scope"],location="Page 6 footnote 13, cyclopentadiene 9q",name="annulation attempt with 1,3-disubstituted cyclopentadiene",
           transformation="unresolved reaction mixture")
p1=participant(r,"9q","reactant",ev(6)); p2=participant(r,"2a","reactant",ev(6)); pc=participant(r,"IPrAuCl","catalyst",ev(6),mol_percent=fq(5,"mol%",ev(6),"5 mol%"),reference="9q"); pa=participant(r,"AgSbF6","additive",ev(6),mol_percent=fq(5,"mol%",ev(6),"5 mol%"),reference="9q"); ps=participant(r,"dce","solvent",ev(6)); pp=participant(r,"complicated_mixture","product",ev(6))
step(r,"s1_attempt",ev(6),[p1,p2,pc,pa,ps,pp],operation="reaction of 9q with 2a",temperature=f("rt","rt",ev(6)),time=fq(2,"h",ev(6),"2 h"))
outcome(r,"mixture","other",["complicated_mixture"],ev(6),"complicated mixture / several minor products","complicated mixture / several minor products",scope="unidentified products")
set_class(r,ev(6),None,coverage="unresolved"); r["review"]["unresolved"].append("Individual products, yields, conversion, and measurement method are not reported.")

# Non-synthetic source context.
bundle["source_context"] = [
    {"context_id":f"{DOC}:context_page1_bioactive","document_id":DOC,"kind":"biological_scaffold_context","evidence_ids":ev(1),"note":"Figure 1 shows selected bioactive molecules containing the core; it is not a current-paper synthesis experiment."},
    {"context_id":f"{DOC}:context_xray","document_id":DOC,"kind":"structure_assignment","evidence_ids":uniq(ev(2)+ev(3)+ev(4)+ev(5)),"note":"X-ray assignments are reported for 5e and (+)-6q (CCDC 2026964 and 2026956); external crystallographic files were not supplied or associated."},
    {"context_id":f"{DOC}:context_dft","document_id":DOC,"kind":"computational_mechanism","evidence_ids":uniq(ev(4)+ev(5)),"note":"DFT energy profiles and ion-pair discussion support the proposed mechanism; they are not additional executed synthesis experiments."},
    {"context_id":f"{DOC}:context_characterization_inventory","document_id":DOC,"kind":"negative_characterization_inventory","evidence_ids":uniq(ev(1)+ev(2)+ev(3)+ev(4)+ev(5)+ev(6)),"note":"No standalone mass-spectrum peak table, molecular-weight table, or fragment-relative-abundance block occurs on pages 1–6."},
]

# Coverage and raw visual inventory.
ids = [r["reaction_id"] for r in bundle["reactions"]]
by_prefix = lambda prefix: [x for x in ids if x.split(":",1)[1].startswith(prefix)]
page_records = {
    1: by_prefix("scheme1_"),
    2: by_prefix("table1_")+by_prefix("table2_"),
    3: by_prefix("table3_")+by_prefix("scheme2_")+by_prefix("table4_entry"),
    4: by_prefix("table4_prose")+by_prefix("table5_")+by_prefix("eq1_")+by_prefix("eq2_")+[f"{DOC}:scheme3_proposed_mechanism"],
    5: [f"{DOC}:scheme3_proposed_mechanism"],
    6: by_prefix("note13_"),
}
regions = {
    1: [
        {"region":"Scheme 1, right-column middle","finding":"Two cited prior-work transformation families and one current generic [3+2]-annulation; not concrete current experiments.","evidence_ids":ev(1)},
        {"region":"Figure 1, lower right","finding":"Bioactive scaffold context only; no current synthesis experiment.","evidence_ids":ev(1)},
        {"region":"title/abstract/introduction","finding":"Title, DOI, scope claim, and explicit online SI notice; no SI file supplied.","evidence_ids":ev(1)},
    ],
    2: [
        {"region":"Table 1, lower left","finding":"All 15 optimization rows inventoried separately; row 4 prose <2% 3a' retained; rows 13–15 use 1a' as starting material.","evidence_ids":ev(2)},
        {"region":"Table 2, lower right","finding":"All 15 scope arrows/products 3b–3p inventoried independently with isolated yields.","evidence_ids":ev(2)},
        {"region":"prose above tables","finding":"Cross-check of catalyst/solvent yields, negligible 3a' in entry 4, and no-product statements for entries 8 and 14.","evidence_ids":ev(2)},
    ],
    3: [
        {"region":"Table 3, upper left","finding":"All 15 diazoketone scope outcomes 4b–4p inventoried separately.","evidence_ids":ev(3)},
        {"region":"Scheme 2, upper right","finding":"One gram-scale annulation and five independent functionalization arrows; dr outcomes attached to 5c and 5d.","evidence_ids":ev(3)},
        {"region":"Table 4, lower right","finding":"All 11 chiral-catalyst optimization rows inventoried with yield and (+)/(−) er.","evidence_ids":ev(3)},
    ],
    4: [
        {"region":"Table 5, lower left","finding":"All 17 enantioselective scope outcomes inventoried with isolated yield and HPLC er.","evidence_ids":ev(4)},
        {"region":"prose following Table 4, upper left","finding":"Two qualitative optimization attempts retained separately: decreased reaction temperature and increased CPA 7d loading to 6 mol%, each reported only as decreased product yield.","evidence_ids":ev(4)},
        {"region":"Equations (1) and (2), upper right","finding":"Two real water-trapping/competition mechanistic controls; Eq. 2 printed-product label conflict retained.","evidence_ids":ev(4)},
        {"region":"Scheme 3 and adjacent prose","finding":"DFT-supported proposed mechanism, not a series of executed intermediate syntheses.","evidence_ids":ev(4)},
    ],
    5: [
        {"region":"left column mechanism continuation","finding":"Continuation of the Scheme 3 ion-pair/DFT discussion; no new executed reaction arrow.","evidence_ids":ev(5)},
        {"region":"summary/acknowledgements/references","finding":"Conclusion and bibliography; no new current-paper synthesis experiment or standalone spectral peak table.","evidence_ids":ev(5)},
    ],
    6: [
        {"region":"footnote 13 and reaction drawings, upper half","finding":"Three distinct boundary-substrate experiments: 9a gives 10a (25%), 9b gives 82% recovery, and 9q gives a complicated mixture.","evidence_ids":ev(6)},
        {"region":"footnote 14 and manuscript dates","finding":"SI-only additional studies are mentioned but SI is not supplied; no records imported from unavailable SI.","evidence_ids":ev(6)},
        {"region":"reference 12","finding":"Review citations only; not imported as current experiments.","evidence_ids":ev(6)},
    ],
}
unresolved = {
    1:["Generic R-group structures and prior-work examples are not concrete experiments."],
    2:["Table 2 does not separately report the post-diazo-addition time/temperature."],
    3:["Scheme 2 yield bases and dr methods are unreported; Table 4 warming transition time is unreported."],
    4:["Equation yields lack measurement bases; Eq. 2 printed label 3a conflicts with substrate numbering."],
    5:["Computational discussion does not establish isolated intermediates or experimental mechanism proof."],
    6:["Footnote 13 yield/recovery methods and 9q product identities are unreported; referenced SI is unavailable."],
}
coverage={"document_id":DOC,"semantic_accuracy_verified":False,"pages":[
    {"page":p,"state":"attempted","regions":regions[p],"record_ids":page_records[p],"unresolved":unresolved[p]} for p in range(1,7)
]}
visual_inventory={"document_id":DOC,"pages":[{"page":p,"image_evidence_id":f"{DOC}:page{p}:image","attempted":True,"regions":regions[p],"record_ids":page_records[p]} for p in range(1,7)],
                  "deduplication_decisions":[
                      "Table 1 entry 6 and Table 2 scope records are different experiments; Table 2 does not repeat 3a.",
                      "Scheme 2 gram-scale 6d is distinct from Table 4 (+)-6d optimization because catalyst chirality, scale, yield, and stereochemical sample differ.",
                      "Table 3 4k and Table 5 (+)-4k are distinct unassigned/enantioenriched samples and experiments.",
                      "Equation 2 printed 3a was not merged with Table 1 3a because the equation uses diazoketone 2k and presents a numbering conflict.",
                  ]}

# Update mechanical task state without claiming semantic verification.
tasks=json.loads((RUN/"tasks.json").read_text())
for t in tasks:
    if t["document_id"]==DOC and t["evidence_id"] in {f"{DOC}:page{p}:text" for p in range(1,7)}:
        t["state"]="attempted"
        t["reason"]="Original page image and extracted text inspected; findings and unresolved items are in coverage/pages_0001_0006.json."

st.write(RUN/"submitted.json",bundle)
st.write(RUN/"coverage/pages_0001_0006.json",coverage)
st.write(RUN/"raw_model/page_extraction.json",visual_inventory)
st.write(RUN/"raw_model/stereo_checks.json",stereo_checks)
st.write(RUN/"tasks.json",tasks)
st.write(RUN/"raw_model/build_summary.json",{
    "document_id":DOC,"molecules":len(bundle["molecules"]),"reactions":len(bundle["reactions"]),
    "specific_experiments":sum(r["record_type"]=="specific_experiment" for r in bundle["reactions"]),
    "machine_review_status":"needs_review","semantic_accuracy_verified":False,
})
print(json.dumps(json.loads((RUN/"raw_model/build_summary.json").read_text()),ensure_ascii=False,indent=2))
