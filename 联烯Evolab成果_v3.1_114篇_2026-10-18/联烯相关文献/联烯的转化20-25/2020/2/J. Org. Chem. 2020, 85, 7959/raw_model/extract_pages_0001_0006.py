#!/usr/bin/env python3
"""Build the page-1--6 v3.1 machine submission for doc_f3d393b861fb548d.

The construction below is deliberately source bounded.  It records machine
interpretation of the six inspected page images and leaves all records in
needs_review state; it is not a source-review promotion.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SKILL = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SKILL))

from standard_records import (  # noqa: E402
    fact,
    new_molecule,
    new_outcome,
    new_participant,
    new_reaction,
    new_step,
    source_mentions,
    stereo_inspect,
)

DOC = "doc_f3d393b861fb548d"
RUN = ROOT / "runs/current"
RAW = RUN / "raw_model"
OUT = RUN / "submitted.json"


def tids(page: int) -> list[str]:
    return [f"{DOC}:page{page}:text", f"{DOC}:page{page}:image"]


def F(raw, normalized=None, ev=(), method="original_extraction", derivation=None,
      status=None):
    if normalized is None:
        normalized = raw
    return fact(raw=raw, normalized=normalized, evidence=ev, method=method,
                derivation=derivation, status=status)


def Q(number, unit, **extra):
    result = {"number": number, "unit": unit}
    result.update(extra)
    return result


def inferred(raw, normalized, ev, derivation):
    return F(raw, normalized, ev, method="inference", derivation=derivation)


bundle = json.loads((RAW / "init_template.json").read_text())


def add_molecule(mid, label, name, pages, smiles=None, *, axial_note=None,
                 source_stereo=None, unresolved=False):
    ev = []
    for p in pages:
        ev.extend(tids(p))
    m = new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, ev)
    m["names"]["original"] = F(name, name, ev)
    m["names"]["english"] = F(name, name, ev)
    m["identity_note"] = F(
        "Machine identity candidate from the inspected source page(s); CAS was not reported and was not guessed.",
        "machine candidate; CAS unreported", ev)
    if unresolved or not smiles:
        m["representation_type"] = "unresolved"
        chirality_state = ("mixture" if source_stereo and "mixture" in source_stereo.lower()
                           else "specified" if source_stereo else "unresolved")
        m["stereochemistry"]["chirality_status"] = inferred(
            "source structure not fully encoded", chirality_state, ev,
            "The inspected drawing/name does not support a complete stereochemical molecular graph.")
        m["stereochemistry"]["geometry_status"] = inferred(
            "source structure not fully encoded", "unresolved", ev,
            "The inspected drawing/name does not support a complete geometry assignment.")
        m["stereochemistry"]["representation_status"] = inferred(
            "no complete candidate encoding", "not_encoded", ev,
            "No complete candidate SMILES is asserted for this source entity.")
        if source_stereo:
            m["stereochemistry"]["assignment"] = F(source_stereo, source_stereo, ev)
        return m

    diag = stereo_inspect(smiles)
    m["representation_type"] = "smiles"
    m["smiles"] = F(
        "candidate transcribed/derived from source structure", smiles, ev,
        method="structure_derivation",
        derivation="Connectivity was derived from the inspected drawing and source name; RDKit parsing is only a graph diagnostic.")
    m["stereochemistry"]["identity_key"] = F(
        "RDKit stereo graph diagnostic", diag["stereo_graph_key"], ev,
        method="structure_derivation", derivation="Computed by standard_records.stereo_inspect from the candidate SMILES.")
    m["stereochemistry"]["isomeric_smiles"] = F(
        "RDKit canonical isomeric candidate", diag["canonical_isomeric_smiles"], ev,
        method="structure_derivation", derivation="Computed by standard_records.stereo_inspect; no source verification is implied.")
    m["stereochemistry"]["atom_centers"] = F(
        "RDKit atom-center diagnostic", diag["atom_centers"], ev,
        method="structure_derivation", derivation="Zero-based candidate-SMILES atom indices from stereo_inspect.")
    m["stereochemistry"]["double_bonds"] = F(
        "RDKit double-bond diagnostic", diag["double_bonds"], ev,
        method="structure_derivation", derivation="Candidate-graph double bonds from stereo_inspect; carbonyls are not E/Z assignments.")
    potential_atom = any(x["type"].startswith("Atom_") and x["specified"] == "Unspecified"
                         for x in diag["potential_stereo"])
    potential_bond = any(x["type"].startswith("Bond_") and x["specified"] == "Unspecified"
                         for x in diag["potential_stereo"])
    if source_stereo:
        m["stereochemistry"]["assignment"] = F(source_stereo, source_stereo, ev)
        chir = "specified"
    elif potential_atom:
        chir = "unspecified"
        m["stereochemistry"]["assignment"] = F(
            "drawing has no wedge/dash assignment", "potential tetrahedral stereochemistry unassigned", ev)
    else:
        chir = "unresolved" if axial_note and "potential" in axial_note else "achiral"
        m["stereochemistry"]["assignment"] = inferred(
            "no assigned tetrahedral center in candidate graph", "no assigned tetrahedral center", ev,
            "stereo_inspect found no potential tetrahedral center; axial stereochemistry is assessed separately.")
    if axial_note:
        m["stereochemistry"]["axial"] = F(axial_note, axial_note, ev,
            method="structure_derivation",
            derivation="Allene terminal substitution was inspected; toolkit axial support is limited.")
    else:
        m["stereochemistry"]["axial"] = inferred(
            "not an allene/atropisomeric axis in the candidate", "not_applicable", ev,
            "No allene or other stereogenic axis is present in the candidate graph.")
    geom = "unspecified" if potential_bond else "not_applicable"
    rep = "partial" if potential_atom or potential_bond or (axial_note and "potential" in axial_note) else "not_applicable"
    m["stereochemistry"]["chirality_status"] = inferred(
        "stereochemistry review of inspected drawing and candidate graph", chir, ev,
        "Status combines the source drawing with stereo_inspect; absence of source tags is not treated as racemic.")
    m["stereochemistry"]["geometry_status"] = inferred(
        "alkene geometry review", geom, ev,
        "Potential stereogenic double bonds were checked; terminal methylenes and carbonyls are not assigned E/Z.")
    m["stereochemistry"]["representation_status"] = inferred(
        "candidate encoding review", rep, ev,
        "Potential unassigned elements remain partial; otherwise no stereochemical encoding is applicable.")
    return m


def add_template(mid, label, name, page, fixed, variable, allowed, instances=()):
    ev = tids(page)
    m = new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, ev)
    m["names"]["original"] = F(name, name, ev)
    m["names"]["english"] = F(name, name, ev)
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": F(fixed, fixed, ev),
        "variables": [{
            "label": variable,
            "attachment_points": F("source drawing attachment point", "source-defined attachment point", ev),
            "definition": F(f"{variable} as drawn/defined in the source", f"source variable {variable}", ev),
            "allowed_values": F(allowed, allowed, ev),
            "constraints": fact(raw="No additional constraints stated on the inspected page",
                                normalized=None, evidence=ev, status="unreadable"),
        }],
        "instance_mapping": [],
        "reason": F("Source uses an R-group/general structure", "generic R-group structure", ev),
    }
    m["identity_note"] = F("Generic source structure; not a concrete molecule", "generic template", ev)
    for key, val in (("chirality_status", "unresolved"), ("geometry_status", "unresolved"),
                     ("representation_status", "not_encoded")):
        m["stereochemistry"][key] = inferred(
            "generic source structure", val, ev,
            "R-group variation prevents a complete concrete stereochemical graph.")
    return m


# Concrete source structures used in the six-page extraction.
mol_specs = [
    ("m_1a", "1a", "N-(pyridin-2-ylmethyl)benzamide", [1, 2, 4, 6], "O=C(NCc1ccccn1)c1ccccc1"),
    ("m_1b", "1b", "2-methyl-N-(pyridin-2-ylmethyl)benzamide", [3, 6], "Cc1ccccc1C(=O)NCc1ccccn1"),
    ("m_1c", "1c", "2-fluoro-N-(pyridin-2-ylmethyl)benzamide", [3, 6], "O=C(NCc1ccccn1)c1ccccc1F"),
    ("m_1d", "1d", "2-methoxy-N-(pyridin-2-ylmethyl)benzamide", [3, 6], "COc1ccccc1C(=O)NCc1ccccn1"),
    ("m_1e", "1e", "N-(pyridin-2-ylmethyl)-1-naphthamide", [3], "O=C(NCc1ccccn1)c1cccc2ccccc12"),
    ("m_1f", "1f", "N-(pyridin-2-ylmethyl)picolinamide", [3], "O=C(NCc1ccccn1)c1ccccn1"),
    ("m_1g", "1g", "N-(pyridin-2-ylmethyl)thiophene-2-carboxamide", [3], "O=C(NCc1ccccn1)c1cccs1"),
    ("m_1h", "1h", "N-(pyridin-2-ylmethyl)acetamide", [3], "CC(=O)NCc1ccccn1"),
    ("m_1i", "1i", "2-phenyl-N-(pyridin-2-ylmethyl)acetamide", [3], "O=C(Cc1ccccc1)NCc1ccccn1"),
    ("m_1j", "1j", "N-(pyridin-2-ylmethyl)acrylamide", [3], "C=CC(=O)NCc1ccccn1"),
    ("m_1k", "1k", "tert-butyl (pyridin-2-ylmethyl)carbamate", [3], "CC(C)(C)OC(=O)NCc1ccccn1"),
    ("m_1l", "1l", "benzyl (2-oxo-2-((pyridin-2-ylmethyl)amino)ethyl)carbamate", [3], "O=C(CNC(=O)OCc1ccccc1)NCc1ccccn1"),
    ("m_1m", "1m", "N-(pyridin-2-ylmethyl)pivalamide", [3], "CC(C)(C)C(=O)NCc1ccccn1"),
    ("m_1n", "1n", "tertiary 2-picolyl benzamide 1n", [2, 3], "O=C(c1ccccc1)N(C)Cc1ccccn1"),
    ("m_1p", "1p", "N-(pyridin-4-ylmethyl)benzamide", [2, 3], "O=C(NCc1ccncc1)c1ccccc1"),
    ("m_11", "11", "N-benzylbenzamide", [4, 5], "O=C(NCc1ccccc1)c1ccccc1"),
    ("m_2a", "2a", "1-methyl-4-(propa-1,2-dien-1-ylsulfonyl)benzene", [1, 2, 3, 4, 5], "C=C=CS(=O)(=O)c1ccc(C)cc1"),
    ("m_2b", "2b", "1-(buta-1,2-dien-1-ylsulfonyl)-4-methylbenzene", [3], "CC=C=CS(=O)(=O)c1ccc(C)cc1"),
    ("m_2c", "2c", "1-(buta-2,3-dien-2-ylsulfonyl)-4-methylbenzene", [3], "C=C=C(C)S(=O)(=O)c1ccc(C)cc1"),
    ("m_3a", "3a", "N-(1-(pyridin-2-yl)-3-tosylbut-3-en-1-yl)benzamide", [1, 2, 4], "C=C(S(=O)(=O)c1ccc(C)cc1)CC(NC(=O)c1ccccc1)c1ccccn1"),
    ("m_4a", "4a", "N-(1-(pyridin-2-yl)-2-(tosylmethyl)allyl)benzamide", [1, 2, 4], "C=C(CS(=O)(=O)c1ccc(C)cc1)C(NC(=O)c1ccccc1)c1ccccn1"),
    ("m_5a", "5a", "4,4'-((prop-2-ene-1,2-diyl)disulfonyl)bis(methylbenzene)", [2, 3, 4, 5], "C=C(S(=O)(=O)c1ccc(C)cc1)CS(=O)(=O)c1ccc(C)cc1"),
    ("m_5b", "5b", "1-methyl-4-((2-(phenylsulfonyl)allyl)sulfonyl)benzene", [3, 4], "C=C(S(=O)(=O)c1ccccc1)CS(=O)(=O)c1ccc(C)cc1"),
    ("m_3q", "3q", "N-(3-(phenylsulfonyl)-1-(pyridin-2-yl)but-3-en-1-yl)benzamide", [3, 4], "C=C(S(=O)(=O)c1ccccc1)CC(NC(=O)c1ccccc1)c1ccccn1"),
]

for mid, label, name, pages, smi in mol_specs:
    axial = None
    if mid == "m_2a":
        axial = "allene has a terminal CH2 end; axial chirality is not stereogenic for this substitution pattern"
    elif mid == "m_2b":
        axial = "potential allene axial stereochemistry; the source page does not assign an axial configuration"
    elif mid == "m_2c":
        axial = "allene has a terminal CH2 end; axial chirality is not stereogenic for this substitution pattern"
    bundle["molecules"].append(add_molecule(mid, label, name, pages, smi, axial_note=axial))

# Scope products 3b-m and minor products 4b-m are concrete connectivity candidates.
acyl = {
    "b": "Cc1ccccc1C(=O)", "c": "O=C(c1ccccc1F)", "d": "COc1ccccc1C(=O)",
    "e": "O=C(c1cccc2ccccc12)", "f": "O=C(c1ccccn1)", "g": "O=C(c1cccs1)",
    "h": "CC(=O)", "i": "O=C(Cc1ccccc1)", "j": "C=CC(=O)",
    "k": "CC(C)(C)OC(=O)", "l": "O=C(CNC(=O)OCc1ccccc1)", "m": "CC(C)(C)C(=O)",
}
for suffix, acyl_part in acyl.items():
    major = f"{acyl_part}NC(c1ccccn1)CC(=C)S(=O)(=O)c1ccc(C)cc1"
    minor = f"{acyl_part}NC(c1ccccn1)C(=C)CS(=O)(=O)c1ccc(C)cc1"
    bundle["molecules"].append(add_molecule(
        f"m_3{suffix}", f"3{suffix}", f"vinyl sulfone product 3{suffix}", [3], major))
    bundle["molecules"].append(add_molecule(
        f"m_4{suffix}", f"4{suffix}", f"allyl sulfone product 4{suffix}", [3], minor))

# Unresolved/mixture/intermediate and operational molecules.
unresolved_specs = [
    ("m_1o", "1o", "chiral oxazolidinone/carbamate substrate 1o", [2, 3], "specified R center shown by wedge in Figure 2"),
    ("m_2a_d2", "2a-d2", "dideuterated allenyl sulfone 2a-d2", [3, 4], None),
    ("m_3a_d", "3a-d", "deuterated vinyl-sulfone product mixture 3a-d", [3, 4], "mixture of positional D isotopologues; stereochemical assignment not encoded"),
    ("m_4a_d", "4a-d", "deuterated allyl-sulfone product mixture 4a-d", [4], "mixture of positional D isotopologues; stereochemical assignment not encoded"),
    ("m_unknown_products", "uncharacterizable products", "uncharacterizable products", [3], None),
]
for mid, label, name, pages, stereo in unresolved_specs:
    bundle["molecules"].append(add_molecule(mid, label, name, pages, unresolved=True, source_stereo=stereo))

simple = {
    "m_thf": ("THF", "tetrahydrofuran", "C1CCOC1", [2, 3, 4, 6]),
    "m_toluene": ("toluene", "toluene", "Cc1ccccc1", [2]),
    "m_toluene_d8": ("toluene-d8", "toluene-d8", None, [4]),
    "m_dmso": ("DMSO", "dimethyl sulfoxide", "CS(C)=O", [2]),
    "m_acetonitrile": ("CH3CN", "acetonitrile", "CC#N", [2]),
    "m_ethanol": ("EtOH", "ethanol", "CCO", [2]),
    "m_ni_otf2": ("Ni(OTf)2", "nickel(II) triflate", None, [1, 2]),
    "m_nats": ("NaTs", "sodium p-toluenesulfinate", None, [2]),
    "m_bht": ("BHT", "butylated hydroxytoluene", "CC(C)(C)c1cc(C)c(O)c(C(C)(C)C)c1", [3, 4]),
    "m_aibn": ("AIBN", "azobisisobutyronitrile", "CC(C#N)(C)N=NC(C)(C)C#N", [3, 4]),
    "m_k2co3": ("K2CO3", "potassium carbonate", None, [3]),
    "m_pyridine": ("pyridine", "pyridine", "n1ccccc1", [4, 5]),
    "m_benzoic_acid": ("benzoic acid", "benzoic acid", "O=C(O)c1ccccc1", [6]),
    "m_2methylbenzoic_acid": ("2-methylbenzoic acid", "2-methylbenzoic acid", "Cc1ccccc1C(=O)O", [6]),
    "m_2fluorobenzoic_acid": ("2-fluorobenzoic acid", "2-fluorobenzoic acid", "O=C(O)c1ccccc1F", [6]),
    "m_2methoxybenzoic_acid": ("2-methoxybenzoic acid", "2-methoxybenzoic acid", "COc1ccccc1C(=O)O", [6]),
    "m_ch2cl2": ("CH2Cl2", "dichloromethane", "ClCCl", [5, 6]),
    "m_dmf": ("DMF", "N,N-dimethylformamide", "CN(C)C=O", [6]),
    "m_oxalyl_chloride": ("oxalyl chloride", "oxalyl chloride", "O=C(Cl)C(=O)Cl", [6]),
    "m_2picolylamine": ("2-picolylamine", "pyridin-2-ylmethanamine", "NCc1ccccn1", [6]),
    "m_nahco3": ("NaHCO3", "sodium bicarbonate", None, [6]),
    "m_naoh": ("NaOH", "sodium hydroxide", None, [6]),
    "m_water": ("water", "water", "O", [6]),
    "m_mgso4": ("MgSO4", "magnesium sulfate", None, [6]),
    "m_silica": ("silica gel", "silica gel", None, [5, 6]),
    "m_etoac": ("EtOAc", "ethyl acetate", "CCOC(C)=O", [3, 6]),
    "m_acetone": ("acetone", "acetone", "CC(C)=O", [3]),
    "m_hexane": ("hexane", "hexane", "CCCCCC", [6]),
}
for mid, (label, name, smi, pages) in simple.items():
    bundle["molecules"].append(add_molecule(mid, label, name, pages, smi, unresolved=smi is None))

# Generic source structures and mechanism intermediates.
bundle["molecules"].extend([
    add_template("m_gen_1", "1", "generic 2-picolyl amide 1", 2,
                 "R-C(=O)-NH-CH2-(2-pyridyl)", "R", "aryl, heteroaryl, alkyl, protected amino-acyl groups shown for 1b-m"),
    add_template("m_gen_3", "3", "generic vinyl sulfone product 3", 2,
                 "R-C(=O)-NH-CH(2-pyridyl)-CH2-C(=CH2)-SO2Ts", "R", "inherits the acyl group of generic substrate 1"),
    add_template("m_gen_4", "4", "generic allyl sulfone product 4", 3,
                 "R-C(=O)-NH-CH(2-pyridyl)-C(=CH2)-CH2-SO2Ts", "R", "inherits the acyl group of generic substrate 1"),
    add_template("m_gen_aryl_acid", "aryl carboxylic acid", "generic aryl carboxylic acid", 6,
                 "Ar-C(=O)OH", "Ar", "aryl group"),
    add_template("m_gen_alkyl_acid", "alkyl carboxylic acid", "generic alkyl carboxylic acid", 6,
                 "R-C(=O)OH", "R", "alkyl group"),
])
for label in ("6", "7", "8", "9", "10", "12"):
    bundle["molecules"].append(add_molecule(
        f"m_int_{label}", label, f"proposed mechanistic intermediate {label}", [4, 5], unresolved=True))


def set_charge(p, key, raw, normalized, ev, reference=None, inherited=False):
    method = "inherited_general_condition" if inherited else "original_extraction"
    p["charges"][key] = F(raw, normalized, ev, method=method)
    if key in ("equivalents", "mol_percent"):
        p["equivalence_reference"] = F(
            f"relative to {reference}", reference, ev, method=method)


def participant(r, mid, role, ev, *, eq=None, molpct=None, mass=None,
                amount=None, volume=None, ref=None, stage="reaction", product_role=None,
                inherited=False):
    p = new_participant(f"{r['reaction_id']}:p{len(r['participants'])+1}", mid, role, ev, stage=stage)
    if product_role:
        p["product_role"] = F(product_role, product_role, ev)
    if eq is not None:
        set_charge(p, "equivalents", f"{eq} equiv", Q(eq, "equiv"), ev, ref, inherited)
    if molpct is not None:
        set_charge(p, "mol_percent", f"{molpct} mol%", Q(molpct, "mol%"), ev, ref, inherited)
    if mass is not None:
        set_charge(p, "mass", mass[0], Q(mass[1], mass[2]), ev, inherited=inherited)
    if amount is not None:
        set_charge(p, "amount", amount[0], Q(amount[1], amount[2]), ev, inherited=inherited)
    if volume is not None:
        set_charge(p, "volume", volume[0], Q(volume[1], volume[2]), ev, inherited=inherited)
    r["participants"].append(p)
    return p


def classify_main(r, ev):
    changes = [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}], "before_order": 1, "after_order": 0, "change_type": "hydrogen_loss"},
        {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain"},
        {"change_id": "b4", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
    ]
    r["classification"]["bond_changes"] = [
        F("partial source-derived net-change candidate", c, ev, method="inference",
          derivation="Derived from the substrate/product drawings; atom correspondence and coproduct mapping remain unresolved.")
        for c in changes
    ]
    labels = [
        {"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1"]},
        {"category_number": 1, "code": "C-H", "supporting_change_ids": ["b3"]},
    ]
    r["classification"]["labels"] = inferred(
        "formed-bond labels from partial net-change model", labels, ev,
        "Taxonomy labels include only formed C-C and net gained C-H; the decreased allene C=C is retained only as a change.")
    r["classification"]["coverage"] = inferred(
        "partial source-derived classification", "partial", ev,
        "No complete atom mapping, proton-source mapping, or coproduct balance is asserted.")
    r["classification"]["limitations"] = [
        "No complete atom-mapped reaction SMILES is asserted.",
        "Internal hydrogen transfer and allene atom correspondence remain unresolved.",
    ]


def classify_cn(r, ev):
    change = {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "N"}], "before_order": 0, "after_order": 1, "change_type": "formed"}
    r["classification"]["bond_changes"] = [inferred(
        "amide C-N formation candidate", change, ev,
        "Derived from carboxylic-acid/amine starting identities and the named amide product; activation coproducts are unmapped.")]
    r["classification"]["labels"] = inferred(
        "C-N formed bond", [{"category_number": 5, "code": "C-N", "supporting_change_ids": ["b1"]}], ev,
        "The amide C-N connection is the only asserted formed-bond label.")
    r["classification"]["coverage"] = inferred(
        "partial source-derived classification", "partial", ev,
        "Acid activation, leaving groups, proton transfers, and coproducts are not completely mapped.")
    r["classification"]["limitations"] = ["No complete atom mapping or coproduct balance is asserted."]


def classify_none(r, ev, reason="No observed transformation was reported."):
    r["classification"]["labels"] = inferred(reason, [], ev, reason)
    r["classification"]["coverage"] = inferred(reason, "not_applicable", ev, reason)
    r["classification"]["limitations"] = [reason]


def classify_unresolved(r, ev, reason):
    r["classification"]["labels"] = inferred(reason, [], ev, reason)
    r["classification"]["coverage"] = inferred(reason, "unresolved", ev, reason)
    r["classification"]["limitations"] = [reason]


def base_reaction(tag, page, rtype, categories, location, name, transformation,
                  unresolved=None, extra_evidence=()):
    ev = tids(page) + list(extra_evidence)
    rid = f"{DOC}:{tag}"
    r = new_reaction(rid, DOC, ev)
    r["record_type"] = rtype
    r["experiment_categories"] = F(categories, categories, ev)
    r["location"] = F(location, location, ev)
    r["reaction_name"] = F(name, name, ev)
    r["transformation_type"] = F(transformation, transformation, ev)
    r["review"] = {
        "status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
        "unresolved": unresolved or ["Machine extraction requires expert source review."],
    }
    bundle["reactions"].append(r)
    return r, ev


def reaction_step(r, ev, operation, temperature=None, time=None, atmosphere=None,
                  stage="reaction", participant_ids=None):
    s = new_step(f"{r['reaction_id']}:s{len(r['steps'])+1}", ev, stage=stage)
    s["participant_ids"] = participant_ids if participant_ids is not None else [p["participant_id"] for p in r["participants"] if p["stage"] == stage]
    s["operation"] = F(operation, operation, ev)
    if temperature is not None:
        s["temperature"] = F(temperature[0], temperature[1], ev)
    if time is not None:
        s["time"] = F(time[0], time[1], ev)
    if atmosphere is not None:
        s["atmosphere"] = F(atmosphere, atmosphere, ev)
    r["steps"].append(s)
    return s


def outcome(r, ev, kind, targets, raw, normalized, *, product_role=None,
            yield_type=None, basis=None, scope=None):
    o = new_outcome(f"{r['reaction_id']}:o{len(r['outcomes'])+1}", kind, targets, ev)
    o["value"] = F(raw, normalized, ev)
    if product_role:
        o["product_role"] = F(product_role, product_role, ev)
    if yield_type:
        o["yield_type"] = F(yield_type, yield_type, ev)
    if basis:
        o["measurement_basis"] = F(basis, basis, ev)
    if scope:
        o["scope"] = F(scope, scope, ev)
    r["outcomes"].append(o)
    return o


def set_solvent(r, ev, mid, name, volume=None, concentration=None, inherited=False):
    method = "inherited_general_condition" if inherited else "original_extraction"
    r["solvent_conditions"]["composition"] = F(name, [mid], ev, method=method)
    if volume:
        r["solvent_conditions"]["total_volume"] = F(volume[0], Q(volume[1], volume[2]), ev, method=method)
    if concentration:
        r["solvent_conditions"]["concentration"] = F(concentration[0], Q(concentration[1], concentration[2]), ev, method=method)


# Scheme 1b generic current transformation.
r, ev = base_reaction("scheme1b_generic", 2, "generic_scheme", ["other"],
    "Scheme 1b", "metal-free coupling of generic 2-picolyl amide 1 with allenyl sulfone 2a",
    "C(sp3)-H functionalization giving vinyl sulfone 3", extra_evidence=tids(1))
participant(r, "m_gen_1", "reactant", ev)
participant(r, "m_2a", "reactant", ev)
participant(r, "m_gen_3", "product", ev, product_role="main")
participant(r, "m_thf", "solvent", ev)
set_solvent(r, ev, "m_thf", "THF")
reaction_step(r, ev, "generic transformation", ("50 °C", Q(50, "°C")))
classify_main(r, ev)


# Table 1 optimization rows.
table_rows = [
    (1, 3.0, "m_toluene", "toluene", 64, "75:25", 120, 1, "m_ni_otf2", 10, None, None),
    (2, 3.0, "m_toluene", "toluene", 67, "75:25", 120, 1, None, None, None, None),
    (3, 2.0, "m_toluene", "toluene", 64, "94:6", 50, 24, None, None, None, None),
    (4, 2.0, "m_dmso", "DMSO", 15, None, 50, 24, None, None, None, None),
    (5, 2.0, "m_acetonitrile", "CH3CN", 68, "67:33", 50, 24, None, None, None, None),
    (6, 2.0, "m_ethanol", "EtOH", 53, "63:37", 50, 24, None, None, None, 9),
    (7, 2.0, "m_thf", "THF", 67, "98:2", 50, 24, None, None, None, None),
    (8, 1.5, "m_thf", "THF", 56, "98:2", 50, 24, None, None, None, 17),
    (9, 2.0, "m_thf", "THF", 0, None, 50, 24, None, None, 2.0, None),
    (10, 2.0, "m_thf", "THF", 0, None, 50, 2, None, None, 1.0, None),
    (11, 2.0, "m_thf", "THF", 0, None, 50, 24, None, None, 0.1, None),
]
for num, eq2a, solv_mid, solv_name, yld, sel, temp, hours, cat_mid, cat_pct, nats_eq, recovery in table_rows:
    unresolved = ["Yield measurement basis is not stated in Table 1."]
    if sel is None:
        unresolved.append("Selectivity cell is blank and is retained as unreported.")
    if num in (9, 10, 11):
        unresolved.append("A zero table yield is preserved; no product was fabricated.")
    r, ev = base_reaction(f"table1_entry_{num:02d}", 2, "specific_experiment", ["optimization"],
        f"Table 1, entry {num}", "optimization of 1a with allenyl sulfone 2a", "C(sp3)-H functionalization", unresolved,
        extra_evidence=tids(1) if num == 1 else ())
    r["general_scheme_ids"] = [f"{DOC}:scheme1b_generic"]
    scale = 0.07 if num == 11 else 0.19
    vol = 0.37 if num == 11 else 1.0
    participant(r, "m_1a", "reactant", ev, eq=1.0, amount=(f"{scale} mmol", scale, "mmol"), ref="m_1a")
    participant(r, "m_2a", "reactant", ev, eq=eq2a, ref="m_1a")
    participant(r, "m_3a", "product", ev, product_role="main")
    participant(r, "m_4a", "product", ev, product_role="minor")
    participant(r, solv_mid, "solvent", ev, volume=(f"{vol} mL", vol, "mL"))
    if cat_mid:
        participant(r, cat_mid, "catalyst", ev, molpct=cat_pct, ref="m_1a")
    if nats_eq is not None:
        participant(r, "m_nats", "additive", ev, eq=nats_eq, ref="m_1a")
    r["reaction_scale"] = F(f"{scale} mmol of 1a", {"reference_molecule_id": "m_1a", "amount": Q(scale, "mmol")}, ev)
    set_solvent(r, ev, solv_mid, solv_name, (f"{vol} mL", vol, "mL"), ("0.19 M", 0.19, "M"))
    reaction_step(r, ev, "stir", (f"{temp} °C", Q(temp, "°C")), (f"{hours} h", Q(hours, "h")))
    basis = "NMR yield" if num == 4 else None
    outcome(r, ev, "yield", ["m_3a", "m_4a"], f"{yld}%", Q(yld, "%"), product_role="combined 3a/4a", basis=basis)
    if sel:
        a, b = map(int, sel.split(":"))
        outcome(r, ev, "selectivity", ["m_3a", "m_4a"], sel,
                {"ratio": [a, b], "order": ["m_3a", "m_4a"]},
                scope="product ratio 3a/4a; measurement method not reported in the table")
    if recovery:
        outcome(r, ev, "recovery", ["m_1a"], f"{recovery}% recovered 1a", Q(recovery, "%"),
                basis="recovered starting material; entry 8 prose states after column chromatography" if num == 8 else "recovered starting material")
    if num in (9, 10, 11):
        outcome(r, ev, "decomposition", ["m_2a"], "complete decomposition of 2a", "complete decomposition",
                basis="reaction observation described in the Table 1 discussion")
        classify_none(r, ev, "No vinyl/allyl sulfone product was reported for this NaTs-additive entry.")
    else:
        classify_main(r, ev)


# Equation 1: explicitly isolated disulfone byproduct, retained separately because its exact Table-1 row is unclear.
r, ev = base_reaction("eq1_disulfone_5a", 2, "specific_experiment", ["mechanistic_study"],
    "Equation 1 and adjacent prose", "formation and isolation of disulfone byproduct 5a",
    "allenyl sulfone consumption forming disulfone", [
        "The displayed equation gives solvent/50 °C/24 h but does not unambiguously identify the corresponding Table 1 entry.",
        "Solvent identity and measurement basis are unreported in this bounded evidence.",
    ])
participant(r, "m_1a", "reactant", ev)
participant(r, "m_2a", "reactant", ev)
participant(r, "m_3a", "product", ev, product_role="main")
participant(r, "m_4a", "product", ev, product_role="minor")
participant(r, "m_5a", "product", ev, product_role="byproduct")
reaction_step(r, ev, "reaction shown in Equation 1", ("50 °C", Q(50, "°C")), ("24 h", Q(24, "h")))
outcome(r, ev, "yield", ["m_5a"], "10% based upon allenyl sulfone equivalents", Q(10, "%"),
        product_role="isolated byproduct", scope="calculated based upon allenyl sulfone equivalents")
classify_unresolved(r, ev, "Disulfone atom mapping and stoichiometric incorporation of two sulfonyl units were not resolved from the page drawing.")


# Figure 1 generic scope and twelve scope rows.
r, ev = base_reaction("figure1_generic", 3, "generic_scheme", ["substrate_scope"],
    "Figure 1 generic arrow", "substrate scope of 2-picolyl amides 1b-m with 2a",
    "C(sp3)-H functionalization to vinyl/allyl sulfones")
participant(r, "m_gen_1", "reactant", ev, eq=1.0, ref="m_gen_1")
participant(r, "m_2a", "reactant", ev, eq=2.0, ref="m_gen_1")
participant(r, "m_gen_3", "product", ev, product_role="main")
participant(r, "m_gen_4", "product", ev, product_role="minor")
participant(r, "m_thf", "solvent", ev, volume=("1 mL", 1, "mL"))
set_solvent(r, ev, "m_thf", "THF", ("1 mL", 1, "mL"), ("0.2 M", 0.2, "M"))
reaction_step(r, ev, "scope reaction; generic arrow shows 48 h", ("50 °C", Q(50, "°C")), ("48 h", Q(48, "h")))
classify_main(r, ev)

scope_rows = {
    "b": (64, "96:4", None), "c": (64, "96:4", None), "d": (54, "97:3", 66),
    "e": (60, "95:5", None), "f": (64, "99:1", None), "g": (57, "94:6", None),
    "h": (59, "94:6", None), "i": (61, "92:8", None), "j": (70, "93:7", None),
    "k": (57, "96:4", None), "l": (60, "95:5", None), "m": (59, "91:9", None),
}
for suffix, (yld, sel, brsm) in scope_rows.items():
    unresolved = [
        "Figure 1 reports averaged yield/selectivity, but the bounded page does not state the yield measurement basis.",
        "Generic arrow says 48 h; caption says modified protocols for solids/oils are in unavailable SI, so product-specific duration is not generalized.",
        "The newly generated stereocenter is drawn without a stereochemical assignment; racemic composition was not inferred.",
    ]
    r, ev = base_reaction(f"figure1_scope_3{suffix}", 3, "specific_experiment", ["substrate_scope"],
        f"Figure 1, products 3{suffix}/4{suffix}", f"reaction of 1{suffix} with 2a", "C(sp3)-H functionalization", unresolved)
    r["general_scheme_ids"] = [f"{DOC}:figure1_generic"]
    participant(r, f"m_1{suffix}", "reactant", ev, eq=1.0, amount=("0.2 mmol", 0.2, "mmol"), ref=f"m_1{suffix}")
    participant(r, "m_2a", "reactant", ev, eq=2.0, ref=f"m_1{suffix}")
    participant(r, f"m_3{suffix}", "product", ev, product_role="main")
    participant(r, f"m_4{suffix}", "product", ev, product_role="minor")
    participant(r, "m_thf", "solvent", ev, volume=("1 mL", 1, "mL"))
    r["reaction_scale"] = F("0.2 mmol of picolyl amide", {"reference_molecule_id": f"m_1{suffix}", "amount": Q(0.2, "mmol")}, ev)
    set_solvent(r, ev, "m_thf", "THF", ("1 mL", 1, "mL"), ("0.2 M", 0.2, "M"), inherited=True)
    time = 72 if suffix == "k" else 24 if suffix == "l" else None
    reaction_step(r, ev, "scope reaction under Figure 1 conditions",
                  ("50 °C", Q(50, "°C")), (f"{time} h", Q(time, "h")) if time else None)
    outcome(r, ev, "yield", [f"m_3{suffix}", f"m_4{suffix}"], f"{yld}%", Q(yld, "%"), product_role="combined product mixture")
    a, b = map(int, sel.split(":"))
    outcome(r, ev, "selectivity", [f"m_3{suffix}", f"m_4{suffix}"], sel,
            {"ratio": [a, b], "order": [f"m_3{suffix}", f"m_4{suffix}"]},
            scope=f"3{suffix}/4{suffix} product ratio; average of reported runs")
    if brsm:
        outcome(r, ev, "other", [f"m_3{suffix}", f"m_4{suffix}"], f"{brsm}% BRSM", Q(brsm, "%"),
                scope="yield based on recovered starting material")
    if suffix == "l":
        r["notes"].append(F(
            "Four runs; one 85:15 run was excluded from the reported average selectivity.",
            {"runs": 4, "excluded_ratio": "85:15"}, ev))
    classify_main(r, ev)


# Figure 2 and adjacent prose: six distinct failed/unreactive experiments.
failed = [
    ("unreactive_1n", "m_1n", "m_2a", 50, "24 h", "unreacted 1n and 2a; 5a observed", "m_5a"),
    ("unreactive_1o", "m_1o", "m_2a", 50, "24 h", "unreacted 1o and 2a; 5a observed", "m_5a"),
    ("unreactive_1p_50c", "m_1p", "m_2a", 50, "24 h", "uncharacterizable products and 5a", "m_unknown_products"),
    ("unreactive_1p_rt", "m_1p", "m_2a", None, None, "decomposition of 1p and formation of 5a", "m_5a"),
    ("unreactive_2b", "m_1a", "m_2b", 50, None, "unreacted 1a; 2b disappeared; no disulfone; uncharacterizable products", "m_unknown_products"),
    ("unreactive_2c", "m_1a", "m_2c", 50, None, "unreacted 1a and 2c; no disulfone", None),
]
for tag, sub, allene, temp, time_raw, obs, product in failed:
    r, ev = base_reaction(tag, 3, "specific_experiment", ["control", "mechanistic_study"],
        "Figure 2 and adjacent prose", f"unreactive-substrate experiment {tag}", "reactivity control", [
            "Exact solvent, charge, and/or duration are missing where the prose only invokes optimized conditions.",
            "Uncharacterizable material is not assigned a concrete structure.",
        ], extra_evidence=tids(2) if tag in ("unreactive_1n", "unreactive_1o", "unreactive_1p_50c") else ())
    participant(r, sub, "reactant", ev, eq=1.0, ref=sub)
    participant(r, allene, "reactant", ev, eq=1.0 if allene in ("m_2b", "m_2c") else 2.0, ref=sub)
    if product:
        participant(r, product, "product", ev, product_role="observed byproduct/uncharacterizable material")
    if tag not in ("unreactive_2b", "unreactive_2c"):
        participant(r, "m_thf", "solvent", ev)
    reaction_step(r, ev, "reactivity test", ("room temperature", "r.t.") if temp is None else (f"{temp} °C", Q(temp, "°C")),
                  (time_raw, Q(24, "h")) if time_raw else None)
    if tag == "unreactive_2c":
        outcome(r, ev, "no_reaction", [], obs, obs, basis="crude 1H NMR")
        classify_none(r, ev)
    elif product == "m_5a":
        outcome(r, ev, "trace", ["m_5a"], obs, "5a observed", basis="crude 1H NMR")
        classify_unresolved(r, ev, "Formation of 5a is qualitative and its complete atom mapping was not resolved.")
    elif product:
        outcome(r, ev, "decomposition", [product], obs, obs, basis="crude 1H NMR")
        classify_unresolved(r, ev, "Uncharacterizable/decomposition products prevent a formed-bond assignment.")


# Mechanistic experiments on pages 3-5.
r, ev = base_reaction("mass_balance_1c_2a", 3, "specific_experiment", ["mechanistic_study"],
    "page 3 mass-balance paragraph", "mass-balance experiment of 1c with 2a", "mechanistic mass-balance study", [
        "Reaction charges, temperature, time, and product yields are not specified in this paragraph.",
        "Baseline material was uncharacterizable and is not assigned a product identity.",
    ])
participant(r, "m_1c", "reactant", ev)
participant(r, "m_2a", "reactant", ev)
participant(r, "m_etoac", "workup_agent", ev, stage="workup")
participant(r, "m_acetone", "workup_agent", ev, stage="workup")
reaction_step(r, ev, "flush chromatography silica several times with ethyl acetate and acetone", stage="workup")
outcome(r, ev, "other", ["m_1c", "m_2a"], "total isolated fractions accounted for the initial mass of the starting materials", "mass balance accounted", basis="weighed isolated fractions and crude/baseline 1H NMR")
classify_unresolved(r, ev, "The paragraph reports mass balance rather than a uniquely assigned transformation.")

for tag, additive, eq_add, yld, sel in (
    ("scheme2a_bht", "m_bht", 5.0, 67, "96:4"),
    ("scheme2b_aibn", "m_aibn", 1.0, 60, "85:15"),
):
    r, ev = base_reaction(tag, 4, "specific_experiment", ["mechanistic_study", "control"],
        f"Scheme 2{'a' if 'bht' in tag else 'b'}", f"radical {'inhibitor' if 'bht' in tag else 'initiator'} experiment",
        "C(sp3)-H functionalization under radical-probe conditions", ["Yield measurement basis is not reported on the inspected page."],
        extra_evidence=tids(3))
    participant(r, "m_1a", "reactant", ev, eq=1.0, ref="m_1a")
    participant(r, "m_2a", "reactant", ev, eq=2.0, ref="m_1a")
    participant(r, additive, "additive", ev, eq=eq_add, ref="m_1a")
    participant(r, "m_3a", "product", ev, product_role="main")
    participant(r, "m_4a", "product", ev, product_role="minor")
    participant(r, "m_thf", "solvent", ev)
    set_solvent(r, ev, "m_thf", "THF")
    reaction_step(r, ev, "stir", ("50 °C", Q(50, "°C")), ("24 h", Q(24, "h")))
    outcome(r, ev, "yield", ["m_3a", "m_4a"], f"{yld}%", Q(yld, "%"), product_role="combined 3a/4a")
    a, b = map(int, sel.split(":"))
    outcome(r, ev, "selectivity", ["m_3a", "m_4a"], sel, {"ratio": [a, b], "order": ["m_3a", "m_4a"]}, scope="3a/4a ratio")
    classify_main(r, ev)

r, ev = base_reaction("scheme2c_crossover", 4, "specific_experiment", ["mechanistic_study"],
    "Scheme 2c", "crossover experiment with disulfone 5b", "crossover vinyl-sulfone formation", ["Combined yield measurement basis is not reported on the inspected page."],
    extra_evidence=tids(3))
for mid, eq in (("m_1a", 1.0), ("m_2a", 1.0), ("m_5b", 2.0)):
    participant(r, mid, "reactant", ev, eq=eq, ref="m_1a")
participant(r, "m_3a", "product", ev, product_role="minor")
participant(r, "m_3q", "product", ev, product_role="main crossover product")
participant(r, "m_thf", "solvent", ev)
set_solvent(r, ev, "m_thf", "THF")
reaction_step(r, ev, "stir", ("50 °C", Q(50, "°C")), ("24 h", Q(24, "h")))
outcome(r, ev, "yield", ["m_3a", "m_3q"], "41%", Q(41, "%"), product_role="combined 3q/3a")
outcome(r, ev, "selectivity", ["m_3q", "m_3a"], "3q/3a = 93:7", {"ratio": [93, 7], "order": ["m_3q", "m_3a"]}, scope="crossover product ratio")
classify_main(r, ev)

r, ev = base_reaction("control_1a_5b_k2co3", 3, "specific_experiment", ["mechanistic_study", "control"],
    "page 3 prose; detailed conditions referred to unavailable SI", "reaction of 1a with 5b in the presence of potassium carbonate",
    "negative mechanistic control", ["Detailed conditions are in unavailable SI and were not inferred."])
participant(r, "m_1a", "reactant", ev)
participant(r, "m_5b", "reactant", ev)
participant(r, "m_k2co3", "reagent", ev)
reaction_step(r, ev, "combine 1a and 5b in the presence of potassium carbonate; detailed conditions are referred to unavailable SI")
outcome(r, ev, "no_reaction", [], "recovered starting materials", "recovered starting materials", basis="source prose; detailed method unavailable")
classify_none(r, ev)

r, ev = base_reaction("scheme2d_deuterium", 4, "specific_experiment", ["mechanistic_study"],
    "Scheme 2d", "deuterium-labeling experiment", "isotope-labeling C(sp3)-H functionalization", [
        "The products are positional-isotopologue mixtures and are not forced into one concrete stereochemical SMILES.",
        "Yield measurement basis is not reported on the inspected page.",
    ], extra_evidence=tids(3))
participant(r, "m_1a", "reactant", ev, eq=1.0, ref="m_1a")
participant(r, "m_2a_d2", "reactant", ev, eq=2.0, ref="m_1a")
participant(r, "m_3a_d", "product", ev, product_role="main isotope mixture")
participant(r, "m_4a_d", "product", ev, product_role="minor isotope mixture")
participant(r, "m_thf", "solvent", ev)
reaction_step(r, ev, "stir", ("50 °C", Q(50, "°C")), ("24 h", Q(24, "h")))
outcome(r, ev, "yield", ["m_3a_d", "m_4a_d"], "61%", Q(61, "%"), product_role="combined isotope-labeled products")
outcome(r, ev, "selectivity", ["m_3a_d", "m_4a_d"], "96:4", {"ratio": [96, 4], "order": ["m_3a_d", "m_4a_d"]}, scope="3a-d/4a-d")
outcome(r, ev, "selectivity", ["m_3a_d"], "allylic C-D ~68.5%; alkenyl C-D ~36.5%", {"allylic_d_percent": 68.5, "alkenyl_d_percent": 36.5}, scope="deuterium content from 1H/2H NMR; at least five labeled components")
outcome(r, ev, "selectivity", ["m_3a_d"], "3a1-d2 44%, 3a2-d3 20%, 3a3-d3 27%, 3a4-d2 5%", {"components_percent": {"3a1-d2": 44, "3a2-d3": 20, "3a3-d3": 27, "3a4-d2": 5}}, scope="positional isotopologue composition shown in Scheme 2d")
classify_main(r, ev)

r, ev = base_reaction("interconversion_probe", 4, "specific_experiment", ["mechanistic_study", "control"],
    "page 4 prose", "probe for conversion of allyl sulfone 4a to vinyl sulfone 3a", "product-interconversion control", [
        "Charges are not stated in the bounded prose.", "The observed final ratio contains newly formed material plus the charged 59:41 mixture."])
for mid, role in (("m_1a", "reactant"), ("m_2a", "reactant"), ("m_3a", "product"), ("m_4a", "product"), ("m_thf", "solvent")):
    participant(r, mid, role, ev, product_role="charged and observed product" if mid in ("m_3a", "m_4a") else None)
reaction_step(r, ev, "react in the presence of an initial 3a/4a mixture", ("50 °C", Q(50, "°C")), ("24 h", Q(24, "h")))
outcome(r, ev, "selectivity", ["m_3a", "m_4a"], "initial 59:41; crude after 24 h 78:22", {"initial_ratio": [59, 41], "final_ratio": [78, 22], "order": ["m_3a", "m_4a"]}, basis="crude 1H NMR", scope="interconversion probe")
classify_main(r, ev)

r, ev = base_reaction("in_situ_nmr_monitoring", 4, "specific_experiment", ["mechanistic_study"],
    "page 4 prose", "in situ 1H NMR monitoring of 1a with 2a", "reaction time-course monitoring", [
        "Charges are not stated in the bounded prose.", "Signals for 4a were below observation; this is not a quantified zero yield."])
participant(r, "m_1a", "reactant", ev)
participant(r, "m_2a", "reactant", ev)
participant(r, "m_3a", "product", ev, product_role="observed product")
participant(r, "m_4a", "product", ev, product_role="not observed at low abundance")
participant(r, "m_toluene_d8", "solvent", ev)
set_solvent(r, ev, "m_toluene_d8", "toluene-d8")
reaction_step(r, ev, "monitor every hour by 1H NMR", ("50 °C", Q(50, "°C")), ("25.5 h", Q(25.5, "h")))
outcome(r, ev, "other", ["m_1a", "m_3a", "m_4a"], "steady increase of 3a correlated with decrease of 1a; 4a not observed due to low abundance", "time-course qualitative observation", basis="hourly 1H NMR")
classify_main(r, ev)

for tag, pyr_eq, obs in (
    ("nbenzyl_pyridine_1eq", 1.0, "complete decomposition of 2a and unreacted N-benzylbenzamide"),
    ("nbenzyl_pyridine_0p2eq", 0.2, "same result, with 5a approximately 17% based on 2 equiv of 2a"),
):
    r, ev = base_reaction(tag, 5, "specific_experiment", ["mechanistic_study", "control"],
        "page 5 prose", "N-benzylbenzamide/2a experiment with pyridine", "intramolecular-pyridine control", [
            "Temperature, time, solvent, charges, and detailed procedure are not stated in the bounded paragraph.",
            "The approximately 17% 5a observation is explicitly only reported for the 0.2 equiv pyridine run.",
        ])
    participant(r, "m_11", "reactant", ev, eq=1.0, ref="m_11")
    participant(r, "m_2a", "reactant", ev, eq=2.0, ref="m_11")
    participant(r, "m_pyridine", "additive", ev, eq=pyr_eq, ref="m_11")
    reaction_step(r, ev, "react N-benzylbenzamide 11 with 2a and pyridine; temperature, solvent, and duration are not stated in this bounded paragraph")
    if pyr_eq == 0.2:
        participant(r, "m_5a", "product", ev, product_role="byproduct")
        outcome(r, ev, "yield", ["m_5a"], "~17% based on 2 equiv of 2a", Q(17, "%", qualifier="approx", raw_text="~17%"), product_role="byproduct", basis="crude 1H NMR", scope="based on 2 equivalents of 2a")
        classify_unresolved(r, ev, "Qualitative/approximate 5a formation is retained without complete atom mapping.")
    else:
        outcome(r, ev, "decomposition", ["m_2a"], obs, "complete decomposition of 2a", basis="source observation")
        classify_none(r, ev, "No coupling product formed; 2a decomposed.")


# Proposed mechanisms/computations are not promoted to observed experiments.
r, ev = base_reaction("scheme3_proposed_mechanism", 4, "proposed_mechanism", ["mechanistic_study"],
    "Scheme 3a-c", "proposed mechanisms for vinyl and allyl sulfone formation", "proposed ionic mechanism", [
        "Arrow-pushing intermediates 6-10 are proposed, not isolated experimental products.",
        "Alternative mechanisms are stated to be in unavailable SI.",
    ], extra_evidence=tids(5))
for mid, role in (("m_1a", "reactant"), ("m_2a_d2", "reactant"), ("m_int_6", "product"),
                  ("m_int_7", "product"), ("m_int_8", "product"), ("m_int_9", "product"),
                  ("m_int_10", "product"), ("m_5a", "product"), ("m_3a_d", "product"), ("m_4a_d", "product")):
    participant(r, mid, role, ev, product_role="proposed intermediate/product" if role == "product" else None)
reaction_step(r, ev, "Scheme 3a: proposed pyridine attack, proton transfer, addition-elimination, and Ts- generation")
reaction_step(r, ev, "Scheme 3b: proposed Ts- addition, deprotonation, and SN2'/SN2 vinyl-sulfone pathways")
reaction_step(r, ev, "Scheme 3c: proposed direct addition pathway to allyl sulfone")
classify_unresolved(r, ev, "Proposed elementary steps are not treated as one experimentally observed net reaction.")

r, ev = base_reaction("scheme4_dft_acidity", 5, "proposed_mechanism", ["mechanistic_study"],
    "Scheme 4a-b", "DFT relative-acidity model for N-H versus C(sp3)-H bonds", "computational acidity comparison", [
        "DFT free energies are computational results, not preparative reaction yields.",
        "Only the displayed model reaction is recorded; no experimental conversion is asserted.",
    ])
for mid, role in (("m_1a", "reactant"), ("m_int_9", "reactant"), ("m_int_10", "product"), ("m_5a", "product"), ("m_11", "reactant"), ("m_int_12", "reactant")):
    participant(r, mid, role, ev, product_role="computed species" if role == "product" else None)
reaction_step(r, ev, "M06-2X/6-311+G(d,p)/SMD(toluene)//B3LYP-D3/6-31G(d) relative-free-energy calculations")
outcome(r, ev, "other", ["m_1a", "m_11", "m_int_12"],
        "1a: N-H 6.4, C-H 17.5; 11: C-H 26.0; 12: N-H -32.7, C-H -43.3 kcal/mol",
        {"1a_NH": Q(6.4, "kcal/mol"), "1a_CH": Q(17.5, "kcal/mol"),
         "11_CH": Q(26.0, "kcal/mol"), "12_NH": Q(-32.7, "kcal/mol"),
         "12_CH": Q(-43.3, "kcal/mol")},
        basis="DFT level stated in Scheme 4 footnote", scope="relative acidity/free-energy values")
classify_unresolved(r, ev, "Computational acid-base arrows are not classified as experimentally observed bond formation.")


# Page-6 general procedures A-D.  They are generic schemes, not duplicated concrete experiments.
for letter, acid_mid, base_mid, base_name in (
    ("A", "m_gen_aryl_acid", "m_nahco3", "sat. aq. sodium bicarbonate"),
    ("B", "m_gen_alkyl_acid", "m_naoh", "sat. aq. sodium hydroxide"),
):
    r, ev = base_reaction(f"general_procedure_{letter.lower()}", 6, "generic_scheme", ["other"],
        f"General Procedure {letter}", f"conversion of {'aryl' if letter == 'A' else 'alkyl'} carboxylic acids to 2-picolyl amides",
        "amide formation via acid chloride activation", ["Substrate identity, exact 2-4 h hold, overnight duration, and chromatography gradient are instance-specific or qualitative."])
    participant(r, acid_mid, "reactant", ev, eq=1.0, ref=acid_mid)
    participant(r, "m_ch2cl2", "solvent", ev)
    participant(r, "m_dmf", "reagent", ev, eq=0.05, ref=acid_mid)
    participant(r, "m_oxalyl_chloride", "reagent", ev, eq=1.1, ref=acid_mid)
    participant(r, "m_2picolylamine", "reactant", ev, eq=1.1, ref=acid_mid)
    participant(r, base_mid, "workup_agent", ev, stage="workup")
    participant(r, "m_water", "workup_agent", ev, stage="workup")
    participant(r, "m_mgso4", "workup_agent", ev, stage="workup")
    participant(r, "m_silica", "purification_agent", ev, stage="purification")
    set_solvent(r, ev, "m_ch2cl2", "dichloromethane", concentration=("0.5 M", 0.5, "M"))
    reaction_step(r, ev, "charge acid, CH2Cl2, and DMF; cool in ice/water bath")
    reaction_step(r, ev, "add oxalyl chloride dropwise", ("ice/water bath", "ice/water bath"), ("over 5 min", Q(5, "min")))
    reaction_step(r, ev, "after 5 min remove bath and warm to room temperature", ("room temperature", "r.t."), ("5 min before warming", Q(5, "min")))
    reaction_step(r, ev, "maintain and monitor by TLC and 1H NMR", ("room temperature", "r.t."), ("2-4 h", {"minimum": 2, "maximum": 4, "unit": "h", "raw_text": "2-4 h"}))
    reaction_step(r, ev, "cool in ice/water bath and add 2-picolylamine dropwise")
    reaction_step(r, ev, "remove bath, warm to room temperature, and stir vigorously", ("room temperature", "r.t."), ("overnight", "overnight"))
    reaction_step(r, ev, f"add {base_name}; separate and wash/extract as described; dry MgSO4; filter; concentrate", stage="workup")
    reaction_step(r, ev, "silica-gel flash chromatography", stage="purification")
    r["procedure"]["addition_order"] = F("acid/CH2Cl2/DMF; oxalyl chloride; then 2-picolylamine", [acid_mid, "m_ch2cl2", "m_dmf", "m_oxalyl_chloride", "m_2picolylamine"], ev)
    classify_cn(r, ev)

for letter, physical, duration in (("C", "solid", None), ("D", "liquid", ("48-72 h", {"minimum": 48, "maximum": 72, "unit": "h", "raw_text": "48-72 h"}))):
    r, ev = base_reaction(f"general_procedure_{letter.lower()}", 6, "generic_scheme", ["other"],
        f"General Procedure {letter}", f"reaction of {physical} 2-picolyl amide with 2a", "C(sp3)-H functionalization", [
            "Exact substrate/allene charges and product-specific completion time are not supplied by the general procedure paragraph.",
            "The modified scope protocol referenced elsewhere is in unavailable SI.",
        ])
    participant(r, "m_gen_1", "reactant", ev)
    participant(r, "m_2a", "reactant", ev)
    participant(r, "m_gen_3", "product", ev, product_role="main generic product")
    participant(r, "m_thf", "solvent", ev)
    participant(r, "m_ch2cl2", "workup_agent", ev, stage="workup")
    participant(r, "m_silica", "purification_agent", ev, stage="purification")
    set_solvent(r, ev, "m_thf", "THF")
    reaction_step(r, ev, f"charge {physical} picolyl amide and 2a{' with THF' if letter == 'D' else ''} in air; seal tube")
    reaction_step(r, ev, "evacuate for 1-2 s and fill with nitrogen, repeated three times", None, ("1-2 s evacuation per cycle", {"minimum": 1, "maximum": 2, "unit": "s", "raw_text": "1-2 s"}), "nitrogen")
    if letter == "C":
        reaction_step(r, ev, "add THF by syringe; wrap cap with parafilm; lower into preheated oil bath and monitor by 1H NMR to disappearance", ("50 °C", Q(50, "°C")))
    else:
        reaction_step(r, ev, "wrap cap with parafilm; lower into preheated oil bath and monitor by 1H NMR", ("50 °C", Q(50, "°C")), duration)
    reaction_step(r, ev, "dilute with dichloromethane, transfer, and concentrate under reduced pressure", stage="workup")
    reaction_step(r, ev, "silica-gel flash chromatography", stage="purification")
    classify_main(r, ev)


# Concrete precursor preparations beginning on page 6.
precursors = [
    ("1a", "m_benzoic_acid", "m_1a", ("1.0 g", 1.0, "g"), ("8.2 mmol", 8.2, "mmol"), ("16 mL", 16, "mL"), ("0.03 mL", 0.03, "mL"), ("0.4 mmol", 0.4, "mmol"), ("0.77 mL", 0.77, "mL"), ("9.0 mmol", 9.0, "mmol"), ("0.92 mL", 0.92, "mL"), ("1.6 g", 1.6, "g"), 94, "40-100% ethyl acetate/hexane", False),
    ("1b", "m_2methylbenzoic_acid", "m_1b", ("300 mg", 300, "mg"), ("2.2 mmol", 2.2, "mmol"), ("4 mL", 4, "mL"), ("0.01 mL", 0.01, "mL"), ("0.1 mmol", 0.1, "mmol"), ("0.21 mL", 0.21, "mL"), ("2.5 mmol", 2.5, "mmol"), ("0.25 mL", 0.25, "mL"), ("353 mg", 353, "mg"), 71, "10-100% ethyl acetate/hexane", False),
    ("1c", "m_2fluorobenzoic_acid", "m_1c", ("350 mg", 350, "mg"), ("2.5 mmol", 2.5, "mmol"), ("5 mL", 5, "mL"), ("0.01 mL", 0.01, "mL"), ("0.1 mmol", 0.1, "mmol"), ("0.24 mL", 0.24, "mL"), ("2.8 mmol", 2.8, "mmol"), ("0.29 mL", 0.29, "mL"), ("481 mg", 481, "mg"), 84, "20-100% ethyl acetate/hexane", False),
    ("1d", "m_2methoxybenzoic_acid", "m_1d", ("350 mg", 350, "mg"), ("2.3 mmol", 2.3, "mmol"), ("4.6 mL", 4.6, "mL"), ("0.01 mL", 0.01, "mL"), ("0.1 mmol", 0.1, "mmol"), ("0.22 mL", 0.22, "mL"), ("2.6 mmol", 2.6, "mmol"), ("0.27 mL", 0.27, "mL"), None, None, None, True),
]
for label, acid, prod, acid_mass, acid_amt, solv_vol, dmf_vol, dmf_amt, ox_vol, ox_amt, amine_vol, prod_mass, yld, gradient, crosses in precursors:
    unresolved = ["The experiment follows General Procedure A; product-specific timing is inherited from that general procedure."]
    if crosses:
        unresolved.append("The paragraph continues on page 7; purification, product mass, and yield were not taken into this page-1--6 pass.")
    r, ev = base_reaction(f"synthesis_{label}", 6, "specific_experiment", ["other"],
        f"Compound synthesis paragraph for {label}", f"synthesis of {label}", "amide formation", unresolved)
    r["general_scheme_ids"] = [f"{DOC}:general_procedure_a"]
    participant(r, acid, "reactant", ev, eq=1.0, mass=acid_mass, amount=acid_amt, ref=acid)
    participant(r, "m_ch2cl2", "solvent", ev, volume=solv_vol)
    p = participant(r, "m_dmf", "reagent", ev, volume=dmf_vol, amount=dmf_amt)
    set_charge(p, "equivalents", "0.05 equiv (General Procedure A)", Q(0.05, "equiv"), ev, acid, inherited=True)
    participant(r, "m_oxalyl_chloride", "reagent", ev, eq=1.1, volume=ox_vol, amount=ox_amt, ref=acid)
    # The page gives amine volume and mmol; preserve both.
    amine_amt = ox_amt
    participant(r, "m_2picolylamine", "reactant", ev, eq=1.1, volume=amine_vol, amount=amine_amt, ref=acid)
    participant(r, prod, "product", ev, product_role="main")
    participant(r, "m_nahco3", "workup_agent", ev, stage="workup")
    participant(r, "m_water", "workup_agent", ev, stage="workup")
    participant(r, "m_mgso4", "workup_agent", ev, stage="workup")
    participant(r, "m_silica", "purification_agent", ev, stage="purification")
    participant(r, "m_etoac", "purification_agent", ev, stage="purification")
    participant(r, "m_hexane", "purification_agent", ev, stage="purification")
    r["reaction_scale"] = F(acid_amt[0] + f" of acid precursor for {label}", {"reference_molecule_id": acid, "amount": Q(acid_amt[1], acid_amt[2])}, ev)
    set_solvent(r, ev, "m_ch2cl2", "dichloromethane", solv_vol, ("0.5 M per General Procedure A", 0.5, "M"), inherited=True)
    reaction_step(r, ev, "perform General Procedure A: acid activation at ice-bath then room temperature, followed by amine addition and overnight stirring", ("ice bath then room temperature", "ice bath then r.t."), ("2-4 h activation plus overnight", "multistage timing inherited"))
    reaction_step(r, ev, "aqueous bicarbonate workup, water wash, MgSO4 drying, filtration, concentration", stage="workup")
    if gradient:
        reaction_step(r, ev, f"silica-gel flash chromatography ({gradient})", stage="purification")
    if prod_mass:
        outcome(r, ev, "product_amount", [prod], prod_mass[0], Q(prod_mass[1], prod_mass[2]), product_role="isolated product")
    if yld is not None:
        outcome(r, ev, "yield", [prod], f"{yld}%", Q(yld, "%"), product_role="isolated product", basis="isolated after silica-gel flash chromatography")
    classify_cn(r, ev)


# Non-synthetic observations retained as source context, not reaction records.
bundle["source_context"].extend([
    {
        "context_id": f"{DOC}:context_bibliography",
        "document_id": DOC,
        "evidence_ids": tids(1),
        "kind": "bibliographic_metadata",
        "description": "Metal-Free C-C Coupling of an Allenyl Sulfone with Picolyl Amides to Access Vinyl Sulfones via Pyridine-Initiated In Situ Generation of Sulfinate Anion; J. Org. Chem. 2020, 85, 7959-7975; DOI 10.1021/acs.joc.0c00788.",
    },
    {
        "context_id": f"{DOC}:context_page2_scheme1a",
        "document_id": DOC,
        "evidence_ids": tids(2),
        "kind": "biological_examples_not_experiments",
        "description": "Scheme 1a depicts selected biologically active vinyl sulfones without reaction arrows; no current-paper synthesis record was created.",
    },
    {
        "context_id": f"{DOC}:context_icp_ms",
        "document_id": DOC,
        "evidence_ids": tids(3),
        "kind": "characterization_or_analysis",
        "description": "ICP-MS of a 10 mg sample of 2a found Ag 4.5 ± 0.4 ug/L; this analytical measurement is not a synthesis experiment.",
    },
    {
        "context_id": f"{DOC}:context_deuterium_ms",
        "document_id": DOC,
        "evidence_ids": tids(4),
        "kind": "characterization_or_analysis",
        "description": "ESI masses 409.15 and 410.16 m/z support di- and tri-deuterated products; peak data were not converted into separate synthesis records.",
    },
    {
        "context_id": f"{DOC}:context_general_methods",
        "document_id": DOC,
        "evidence_ids": tids(5) + tids(6),
        "kind": "general_characterization_methods",
        "description": "NMR, IR, HRMS, melting-point, chromatography, solvent-purification and ICP-MS instrument descriptions are source context, not standalone reactions.",
    },
])


# Retain literal mention hints separately; they are not evidence reconciliation.
dataset = json.loads((RUN / "dataset.json").read_text())
mentions = source_mentions(dataset, ["1a", "2a", "3a", "4a", "5a", "1n", "1o", "1p"], context=220)
(RAW / "source_mentions_key_labels.json").write_text(json.dumps(mentions, ensure_ascii=False, indent=2) + "\n")

# Save all concrete stereo diagnostics separately for reviewer inspection.
stereo_checks = {}
for m in bundle["molecules"]:
    smi = m["smiles"]["normalized_value"]
    if smi:
        stereo_checks[m["molecule_id"]] = stereo_inspect(smi)
(RAW / "stereo_checks.json").write_text(json.dumps(stereo_checks, ensure_ascii=False, indent=2) + "\n")

OUT.write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]), "source_context": len(bundle["source_context"])}, ensure_ascii=False))
