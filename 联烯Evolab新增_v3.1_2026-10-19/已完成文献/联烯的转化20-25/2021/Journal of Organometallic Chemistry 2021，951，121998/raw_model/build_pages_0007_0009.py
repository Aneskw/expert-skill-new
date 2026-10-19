#!/usr/bin/env python3
"""Build the page 7--9 machine bundle from visually inspected source images."""

from copy import deepcopy
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
SKILL = ROOT / "skills" / "chemical-reaction-extraction" / "scripts"
sys.path.insert(0, str(SKILL))
import standard_records as st  # noqa: E402

DATASET = ROOT / "runs" / "current" / "raw_model" / "dataset_before_pages_0007_0009.json"
PREPARED = ROOT / "runs" / "current" / "dataset.json"
OUT = ROOT / "runs" / "current" / "submitted.json"
RAW = ROOT / "runs" / "current" / "raw_model"
DOC = "doc_e6209d110a05e1c9"
P5 = f"{DOC}:page5:image"
P6 = f"{DOC}:page6:image"
P7 = f"{DOC}:page7:image"
P8 = f"{DOC}:page8:image"
P9 = f"{DOC}:page9:image"


def add_unique(seq, value):
    if value not in seq:
        seq.append(value)


def index(items, key, value):
    return next(x for x in items if x[key] == value)


def set_name(mol, name, evidence):
    mol["names"]["original"] = st.fact(name, name, evidence)
    mol["names"]["english"] = st.fact(name, name, evidence)


def stereo_candidate(mol, smiles, chirality, atom_note, evidence):
    first = st.stereo_inspect(smiles)
    canonical = first["canonical_isomeric_smiles"]
    diag = st.stereo_inspect(canonical)
    deriv = (
        "Candidate connectivity transcribed from the source compound name/drawing and "
        "canonicalized by RDKit; parsing is not source verification."
    )
    mol["representation_type"] = "smiles"
    mol["smiles"] = st.fact(None, canonical, evidence, "structure_derivation", derivation=deriv)
    stereo = mol["stereochemistry"]
    stereo["atom_centers"] = st.fact(atom_note, atom_note, evidence)
    stereo["double_bonds"] = st.fact(
        "terminal CH2-substituted allene; no E/Z assignment reported",
        "terminal CH2-substituted allene; E/Z not applicable",
        evidence,
        "inference",
        derivation="The terminal allene carbon is CH2 in the source structure/name; no E/Z or axial assignment is introduced.",
    )
    stereo["axial"] = st.fact(
        "No allene axial configuration is assigned; one terminus is CH2.",
        "not applicable for the depicted terminal allene",
        evidence,
        "inference",
        derivation="Axial chirality requires two different substituents at each allene terminus; the depicted terminal carbon is CH2.",
    )
    stereo["mixture"] = st.fact(evidence=evidence, status="not_applicable")
    stereo["assignment"] = st.fact(
        "No absolute configuration is assigned." if chirality == "unspecified" else "No stereogenic atom center in the candidate connectivity.",
        "configuration unassigned" if chirality == "unspecified" else "achiral candidate connectivity",
        evidence,
    )
    stereo["chirality_status"] = st.fact(None, chirality, evidence, "inference", derivation="Status follows the source structure/name without treating unassigned configuration as achiral or racemic.")
    stereo["geometry_status"] = st.fact(None, "not_applicable", evidence, "inference", derivation="Terminal allene does not define an E/Z alkene in the depicted molecule.")
    stereo["representation_status"] = st.fact(None, "partial", evidence, "inference", derivation="Connectivity is encoded; source stereochemical assignment is absent or toolkit allene diagnostics are not treated as an axial assignment.")
    stereo["isomeric_smiles"] = st.fact(None, canonical, evidence, "structure_derivation", derivation="Canonical isomeric SMILES from the candidate graph; unassigned elements remain unassigned.")
    stereo["identity_key"] = st.fact(None, diag["stereo_graph_key"], evidence, "structure_derivation", derivation="RDKit stereo_graph_key from the candidate graph; not a sample identity.")
    return {"molecule_id": mol["molecule_id"], "input_smiles": smiles, "diagnostic": diag}


def conflict_fact(raw_values, alternatives, evidence):
    return st.fact(raw_values, None, evidence, status="conflict", alternatives=alternatives)


dataset = st.read(DATASET)
bundle = deepcopy(dataset["standard"])
# The official prepare pass registered page-image evidence for all nine pages.
# Only the document/evidence registry is replaced; the prior machine records are
# retained as the starting point for this bounded page update.
bundle["documents"] = deepcopy(st.read(PREPARED)["documents"])
molecules = {m["molecule_id"]: m for m in bundle["molecules"]}
reactions = {r["reaction_id"]: r for r in bundle["reactions"]}

# Resolve the audit's reactant structure gaps with source-derived candidates. These
# remain inferred, partial, and needs-review; successful parsing is not verification.
candidate_smiles = {
    "4a": "COC(=O)C(N(C)Cc1cn(COC(=O)C(C)(C)C)nn1)(C(F)(F)F)C=C=C",
    "4b": "COC(=O)C(N(C)Cc1cn(Cc2ccccc2)nn1)(C(F)(F)F)C=C=C",
    "4c": "CCOC(=O)C(N(C)Cc1cn(COC(=O)C(C)(C)C)nn1)(C(=O)OCC)C=C=C",
    "4d": "CCOC(=O)C(N(C)Cc1cn(Cc2ccccc2)nn1)(C(=O)OCC)C=C=C",
    "4e": "CC(C)(C)C(=O)OCn1cc(CN(C)C(C(F)(F)F)(C(F)(F)F)C=C=C)nn1",
    "4f": "c1ccccc1Cn2cc(CN(C)C(C(F)(F)F)(C(F)(F)F)C=C=C)nn2",
}
diagnostics = []
for label, smiles in candidate_smiles.items():
    mol = molecules[f"{DOC}:mol:{label}"]
    chirality = "unspecified" if label in {"4a", "4b"} else "achiral"
    atom_note = "potential tetrahedral center; configuration not reported" if chirality == "unspecified" else "no stereogenic atom center in the candidate connectivity"
    diagnostics.append(stereo_candidate(mol, smiles, chirality, atom_note, tuple(mol["occurrences"][0]["evidence_ids"])))

# Page 7 supplies complete product names and explicitly describes mixtures. Keep
# product SMILES blank because a registry entry represents multiple stereoisomers.
product_names = {
    "5a": "Dimethyl 2,2′-(cyclobutane-1,3-diylidenebis(methanylylidene))bis(3,3,3-trifluoro-2-(methyl((1-((pivaloyloxy)methyl)-1H-1,2,3-triazol-4-yl)methyl)amino)propanoate)",
    "5b": "Dimethyl 2,2′-(cyclobutane-1,3-diylidenebis(methanylylidene))bis(2-(((1-benzyl-1H-1,2,3-triazol-4-yl)methyl)(methyl)amino)-3,3,3-trifluoropropanoate)",
    "5c": "Tetraethyl 2,2′-(cyclobutane-1,3-diylidenebis(methanylylidene))bis(2-(methyl((1-((pivaloyloxy)methyl)-1H-1,2,3-triazol-4-yl)methyl)amino)malonate)",
    "5d": "Tetraethyl 2,2′-(cyclobutane-1,3-diylidenebis(methanylylidene))bis(2-(((1-benzyl-1H-1,2,3-triazol-4-yl)methyl)(methyl)amino)malonate)",
    "5e": "(4,4′-((((Cyclobutane-1,3-diylidenebis(methanylylidene))bis(1,1,1,3,3,3-hexafluoropropane-2,2-diyl))bis(methylazanediyl))bis(methylene))bis(1H-1,2,3-triazole-4,1-diyl))bis(methylene) bis(2,2-dimethylpropanoate)",
    "5f": "2,2′-(Cyclobutane-1,3-diylidenebis(methanylylidene))bis(N-((1-benzyl-1H-1,2,3-triazol-4-yl)methyl)-1,1,1,3,3,3-hexafluoro-N-methylpropan-2-amine)",
}
page7_ratios = {"5a": "~1:1", "5b": "~1:1", "5c": "~1:2", "5d": "~1:1", "5e": "3:1", "5f": "4:1"}
for label, name in product_names.items():
    mol = molecules[f"{DOC}:mol:{label}"]
    set_name(mol, name, (P7,))
    add_unique(mol["occurrences"][0]["evidence_ids"], P7)
    mol["identity_note"] = st.fact(
        f"{label} is reported as a Z/E mixture" + (" and as diastereomers A/B" if label in {"5a", "5b"} else "") + "; no unique single-isomer graph is asserted.",
        "reported stereoisomer mixture; unique isomeric SMILES intentionally unresolved",
        (P7,),
    )
    stereo = mol["stereochemistry"]
    stereo["double_bonds"] = st.fact(
        f"Z/E mixture, page-7 experimental description {page7_ratios[label]}",
        f"Z/E mixture ({page7_ratios[label]})",
        (P7,),
    )
    stereo["axial"] = st.fact(
        "No allene remains in the named cyclobutane dimer product.",
        "not applicable",
        (P7,),
        "inference",
        derivation="The source product name/drawing is a dimethylenecyclobutane rather than an axially chiral allene.",
    )
    if label in {"5a", "5b"}:
        stereo["atom_centers"] = st.fact("diastereomers A and B; absolute configurations not assigned", "unassigned stereogenic centers represented by A/B fractions", (P7,))
        stereo["chirality_status"] = st.fact(None, "mixture", (P7,), "inference", derivation="The source explicitly reports a 1:1 mixture of diastereomers A and B without absolute configurations.")
        stereo["mixture"] = st.fact(f"diastereomers A/B = 1:1; Z/E = {page7_ratios[label]}", {"diastereomer_ratio": "1:1", "geometry_ratio": page7_ratios[label]}, (P7,))
    else:
        duplicate_group = "two ester substituents at each candidate center" if label in {"5c", "5d"} else "two CF3 substituents at each candidate center"
        stereo["atom_centers"] = st.fact(
            None,
            f"no stereogenic atom center: {duplicate_group}",
            (P7,),
            "inference",
            derivation="The complete source name/drawing gives two identical substituents at each otherwise candidate tetrahedral center.",
        )
        stereo["assignment"] = st.fact(
            "Z/E mixture; no atom-centered configuration is assigned.",
            "achiral with respect to atom centers; geometry remains a mixture",
            (P7,),
        )
        stereo["chirality_status"] = st.fact(
            None,
            "achiral",
            (P7,),
            "inference",
            derivation="Each otherwise candidate tetrahedral center bears duplicate ester groups (5c/5d) or duplicate CF3 groups (5e/5f).",
        )
        stereo["mixture"] = st.fact(f"Z/E = {page7_ratios[label]}", {"geometry_ratio": page7_ratios[label]}, (P7,))
    stereo["geometry_status"] = st.fact(None, "mixture", (P7,), "inference", derivation="The source explicitly reports Z/E mixtures.")
    stereo["representation_status"] = st.fact(None, "not_encoded", (P7,), "inference", derivation="The registry entry is a reported stereoisomer mixture, so no unique isomeric SMILES is asserted.")
    stereo["isomeric_smiles"] = st.fact(evidence=(P7,))
    stereo["identity_key"] = st.fact(evidence=(P7,))

# Page 7 continues the same 2d -> 3d synthesis and reports an isolated B/Z fraction.
fraction_id = f"{DOC}:mol:3d_B_Z"
fraction = st.new_molecule(fraction_id, DOC, (P7,))
fraction["identity_status"] = "source_identified"
fraction["representation_type"] = "unresolved"
fraction["occurrences"][0]["compound_label"] = st.fact("Diastereomer B, Z-isomer", "3d Diastereomer B, Z-isomer", (P7,))
set_name(fraction, "Diastereomer B, Z-isomer of 3d", (P7,))
fraction["names"]["abbreviations"] = st.fact("3d, Diastereomer B, Z-isomer", ["3d-B-Z"], (P7,))
fraction["stereochemistry"]["atom_centers"] = st.fact("Diastereomer B; absolute configuration not reported", "relative fraction label B; absolute configuration unresolved", (P7,))
fraction["stereochemistry"]["double_bonds"] = st.fact("Z-isomer", "Z (source assignment; atom mapping unavailable)", (P7,))
fraction["stereochemistry"]["axial"] = st.fact(evidence=(P7,), status="not_applicable")
fraction["stereochemistry"]["mixture"] = st.fact(evidence=(P7,), status="not_applicable")
fraction["stereochemistry"]["assignment"] = st.fact("Diastereomer B, Z-isomer", "B/Z source label; absolute configuration unresolved", (P7,))
fraction["stereochemistry"]["chirality_status"] = st.fact(None, "unspecified", (P7,), "inference", derivation="B is a relative fraction label, not an absolute stereochemical assignment.")
fraction["stereochemistry"]["geometry_status"] = st.fact(None, "specified", (P7,), "inference", derivation="The source explicitly calls the isolated fraction the Z-isomer.")
fraction["stereochemistry"]["representation_status"] = st.fact(None, "not_encoded", (P7,), "inference", derivation="The source lacks atom-mapped stereochemical connectivity for a unique isomeric SMILES.")
fraction["identity_note"] = st.fact("Separate chromatographic fraction; no absolute configuration supplied.", "source-labelled B/Z fraction of 3d", (P7,))
bundle["molecules"].append(fraction)
molecules[fraction_id] = fraction

r3d = reactions[f"{DOC}:scheme2_2d_to_3d"]
add_unique(r3d["evidence_ids"], P7)
r3d["location"] = st.fact({"pdf_pages": [3, 6, 7]}, {"pdf_pages": [3, 6, 7]}, (f"{DOC}:page3:image", P6, P7))
part = st.new_participant("product_B_Z", fraction_id, "product", (P7,), stage="purification")
part["product_role"] = st.fact("separated Diastereomer B, Z-isomer fraction", "separated stereoisomer fraction", (P7,))
r3d["participants"].append(part)
index(r3d["steps"], "step_id", "purification")["participant_ids"].append("product_B_Z")
out = st.new_outcome("yield_B_Z_fraction", "yield", [fraction_id], (P7,))
out["product_role"] = st.fact("Diastereomer B, Z-isomer", "isolated stereoisomer fraction", (P7,))
out["value"] = st.fact("yield 21%", {"number": 21, "unit": "%"}, (P7,))
out["yield_type"] = st.fact("isolated fraction yield", "isolated yield", (P7,))
out["measurement_basis"] = st.fact("reported after column chromatographic separation", "isolated after column chromatography", (P7,))
out["scope"] = st.fact({"stereoisomer": "Diastereomer B, Z-isomer"}, {"stereoisomer": "B/Z fraction"}, (P7,))
r3d["outcomes"].append(out)
r3d["review"]["unresolved"].append("Page 7 labels the 21% fraction only as Diastereomer B, Z-isomer; absolute configuration and atom-mapped Z bond remain unresolved.")

# Merge page-7 compound-specific workup/results into the six existing Scheme 4
# reactions. Different source values remain explicit conflicts.
scheme_yields = {"5a": 78, "5b": 76, "5c": 76, "5d": 70, "5e": 78, "5f": 70}
page7_yields = {"5a": 78, "5b": 73, "5c": 76, "5d": 70, "5e": 78, "5f": 70}
scheme_ratios = {"5a": "1:1", "5b": "3:1", "5c": "1:2", "5d": "1:1", "5e": "3:1", "5f": "4:1"}
eluents = {"5a": "8/1", "5b": "5/1", "5c": "8/1", "5d": "6/1", "5e": "3/1", "5f": "3/1"}
for label in product_names:
    rid = f"{DOC}:scheme4_4{label[-1]}_to_{label}"
    r = reactions[rid]
    add_unique(r["evidence_ids"], P7)
    r["location"] = st.fact({"pdf_pages": [5, 6, 7]}, {"pdf_pages": [5, 6, 7]}, (P5, P6, P7))
    purification = index(r["steps"], "step_id", "purification")
    purification["operation"] = st.fact(
        f"column chromatography on silica gel, petroleum ether/ethyl acetate = {eluents[label]}",
        f"column chromatography on silica gel with petroleum ether/ethyl acetate = {eluents[label]}",
        (P7,),
    )
    y = index(r["outcomes"], "outcome_id", "yield")
    y["target_molecule_ids"] = [f"{DOC}:mol:{label}"]
    y["yield_type"] = st.fact("Yield reported after synthesis and column chromatography", "isolated yield", (P7,))
    y["measurement_basis"] = st.fact(
        "Section 4.7 compound entry follows the general synthesis and compound-specific column chromatography",
        "isolated after column chromatography",
        (P6, P7),
    )
    if label == "5b":
        alt_scheme = st.fact("76% in Scheme 4", {"number": 76, "unit": "%"}, (P5,))
        alt_exp = st.fact("Yield 73% in Section 4.7", {"number": 73, "unit": "%"}, (P7,))
        y["value"] = conflict_fact(["76%", "73%"], [alt_scheme, alt_exp], (P5, P7))
    else:
        y["value"] = st.fact(f"Yield {page7_yields[label]}%", {"number": page7_yields[label], "unit": "%"}, (P7,))
    sel = index(r["outcomes"], "outcome_id", "selectivity")
    sel["target_molecule_ids"] = [f"{DOC}:mol:{label}"]
    sel["measurement_basis"] = st.fact(evidence=(P5, P7))
    sel["scope"] = st.fact({"metric": "Z/E ratio", "order": ["Z", "E"]}, {"metric": "Z/E ratio", "order": ["Z", "E"]}, (P5, P7))
    if label == "5b":
        a = st.fact("Z/E = 3:1 in Scheme 4", "Z/E = 3:1", (P5,))
        b = st.fact("Mixture of Z/E isomers in ratio of ~1:1 in Section 4.7", "Z/E ≈ 1:1", (P7,))
        sel["value"] = conflict_fact(["Z/E = 3:1", "Z/E ~1:1"], [a, b], (P5, P7))
    else:
        sel["value"] = st.fact(f"Mixture of Z/E isomers in ratio of {page7_ratios[label]}", f"Z/E = {page7_ratios[label]}", (P7,))
    if label in {"5a", "5b"}:
        dr = st.new_outcome("diastereomer_ratio", "selectivity", [f"{DOC}:mol:{label}"], (P7,))
        dr["product_role"] = st.fact("target product mixture", "target product mixture", (P7,))
        dr["value"] = st.fact("diastereomers A and B at the ratio 1:1", "A/B = 1:1", (P7,))
        dr["measurement_basis"] = st.fact(evidence=(P7,))
        dr["scope"] = st.fact({"metric": "diastereomer ratio", "order": ["A", "B"]}, {"metric": "diastereomer ratio", "order": ["A", "B"]}, (P7,))
        r["outcomes"].append(dr)
    r["notes"].append(st.fact(
        f"Section 4.7 reports {label}: {page7_yields[label]}% and Z/E {page7_ratios[label]}; eluent petroleum ether/ethyl acetate {eluents[label]}.",
        "page-7 compound-specific yield, stereoisomer ratio, and purification retained",
        (P7,),
    ))
    r["review"]["unresolved"] = [u for u in r["review"]["unresolved"] if "Page 7 was context only" not in u and "Yield measurement basis is not stated in Scheme 4" not in u]
    r["review"]["unresolved"].append("The source does not explicitly assign a measurement method to the reported Z/E ratio; measurement_basis remains missing.")
    if label == "5b":
        r["review"]["unresolved"].append("Source conflict retained: Scheme 4 gives 76% and Z/E 3:1, while the page-7 Section 4.7 entry gives 73% and Z/E approximately 1:1.")

bundle["source_context"].extend([
    {
        "context_id": f"{DOC}:context:characterization_page7",
        "kind": "characterization",
        "evidence_ids": [P7],
        "description": "Continuation of 3d and characterization blocks for 5a-f (NMR and elemental analysis), plus Section 4.8 X-ray instrumentation/refinement and CCDC identifiers; not counted as separate syntheses.",
    },
    {
        "context_id": f"{DOC}:context:backmatter_pages8_9",
        "kind": "general_information",
        "evidence_ids": [P8, P9],
        "description": "Acknowledgments, an online supplementary-material notice, and references. No current-paper reaction arrow/table row appears; the notice does not establish an associated SI file in this workspace.",
    },
])

# Every machine record remains explicitly unreviewed.
for reaction in bundle["reactions"]:
    reaction["review"]["status"] = "needs_review"
    reaction["review"]["method"] = "model"
    reaction["review"]["reviewer"] = None
    reaction["review"]["checked_fields"] = []

st.write(OUT, bundle)
st.write(RAW / "pages_0007_0009_stereo_diagnostics.json", {
    "semantic_accuracy_verified": False,
    "limitations": "RDKit diagnostics validate candidate syntax/graph only; they do not verify source identity or stereochemical correctness.",
    "candidates": diagnostics,
})
st.write(RAW / "pages_0007_0009_mentions.json", st.source_mentions(dataset, ["Diastereomer B", "Yield 73%", "Supplementary materials"], context=220))
st.write(RAW / "pages_0007_0009_visual_inventory.json", {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "inspected_images": [P7, P8, P9],
    "page7": {
        "synthesis_related": ["3d Diastereomer B, Z-isomer, 21%", "5a-f compound-specific yield, stereoisomer ratio, and chromatography eluent"],
        "characterization_only": ["NMR and elemental analysis for 3d/5a-f", "Section 4.8 X-ray method and CCDC data"],
    },
    "page8": {"backmatter": ["acknowledgments", "online supplementary-material notice", "references 1-15"]},
    "page9": {"backmatter": ["references 15-21"]},
})
print(f"wrote {OUT}")
