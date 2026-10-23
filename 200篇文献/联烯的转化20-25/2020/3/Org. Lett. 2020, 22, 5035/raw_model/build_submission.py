#!/usr/bin/env python3
"""Build the evidence-linked v3.1 transcription for Org. Lett. 2020, 22, 5035–5040."""

from __future__ import annotations

import copy
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills" / "chemical-reaction-extraction" / "scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st  # noqa: E402
import reaction_tools as rt  # noqa: E402


RUN = ROOT / "runs" / "current"
DATASET = RUN / "dataset.json"
DOC = "doc_cab8818902b8627f"
TITLE = "Dearomative Cycloadditions Utilizing an Organic Photosensitizer: An Alternative to Iridium Catalysis"


def write(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def img(page: int) -> str:
    return f"{DOC}:page{page}:image"


def txt(page: int) -> str:
    return f"{DOC}:page{page}:text"


def F(raw=None, normalized=None, evidence=(), **kwargs):
    if normalized is None and raw is not None:
        normalized = raw
    return st.fact(raw, normalized, evidence, **kwargs)


def D(normalized, evidence, derivation, raw=None):
    return st.fact(raw, normalized, evidence, method="structure_derivation", derivation=derivation)


def I(normalized, evidence, derivation, raw=None):
    return st.fact(raw, normalized, evidence, method="inference", derivation=derivation)


def NA():
    return st.fact(method="not_applicable", status="not_applicable")


dataset = json.loads(DATASET.read_text())
dataset["documents"][0]["title"] = TITLE
dataset["documents"][0]["supplementary_material"] = {
    "status": "unknown",
    "note": "The main text advertises Supporting Information, but no SI file was explicitly supplied or associated in this session.",
}
write(DATASET, dataset)
bundle = st.empty(dataset["documents"])


stereo_checks = {}


def molecule(
    mid,
    label,
    page,
    *,
    name=None,
    smiles=None,
    chirality="unresolved",
    geometry="unresolved",
    representation_status=None,
    mixture=None,
    assignment=None,
    axial=None,
    note=None,
    extra_evidence=(),
):
    ev = [img(page), *extra_evidence]
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, ev)
    m["names"]["original"] = F(name or label, name or label, ev)
    m["names"]["english"] = F(name, name, ev) if name else F()
    m["names"]["abbreviations"] = F(label, [label], ev) if label != (name or label) else F()
    m["names"]["aliases"] = F()
    m["cas_numbers"] = F()
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = D(
            smiles,
            ev,
            "Candidate connectivity was transcribed from the source-drawn structure or explicit chemical name; parsing does not verify source identity.",
            raw={"source": "source drawing/name", "candidate_smiles": smiles},
        )
        diag = st.stereo_inspect(smiles)
        stereo_checks[mid] = diag
        m["stereochemistry"]["isomeric_smiles"] = D(
            diag["canonical_isomeric_smiles"], ev, "RDKit canonical isomeric SMILES from the candidate graph.", raw=smiles
        )
        m["stereochemistry"]["identity_key"] = D(
            diag["stereo_graph_key"], ev, "RDKit stereo graph key from the candidate graph.", raw=smiles
        )
        representation_status = representation_status or "encoded"
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = F()
        representation_status = representation_status or "unresolved"
    if chirality == "racemic":
        m["stereochemistry"]["chirality_status"] = F("rac-" if str(label).startswith("rac-") else "racemic", "racemic", ev)
        m["stereochemistry"]["mixture"] = F(
            "rac- label in source", {"type": "racemic", "composition": "not numerically reported"}, ev
        )
    else:
        m["stereochemistry"]["chirality_status"] = I(
            chirality,
            ev,
            "Status assigned from the inspected source drawing/name; unspecified or unresolved is retained when the drawing does not establish configuration.",
            raw=label,
        )
        m["stereochemistry"]["mixture"] = F(mixture, mixture, ev) if mixture is not None else F()
    m["stereochemistry"]["geometry_status"] = I(
        geometry,
        ev,
        "E/Z and ring-geometry status was assessed from the inspected source drawing without inventing absent descriptors.",
        raw=label,
    )
    m["stereochemistry"]["representation_status"] = I(
        representation_status,
        ev,
        "Representation status compares the candidate registry representation with source-depicted stereochemistry.",
        raw=label,
    )
    m["stereochemistry"]["assignment"] = F(assignment, assignment, ev) if assignment else F()
    m["stereochemistry"]["axial"] = I(
        axial,
        ev,
        "Allene axial-stereogenicity assessment derived from terminal substitution visible in the source drawing.",
        raw=label,
    ) if axial else F()
    m["stereochemistry"]["atom_centers"] = F()
    m["stereochemistry"]["double_bonds"] = F()
    m["identity_note"] = I(
        note or "Source-identified compound; complete structure was not reconstructed unless a candidate SMILES is present.",
        ev,
        "Identity note records the boundary of the transcription.",
        raw=label,
    )
    m["external_properties"] = []
    bundle["molecules"].append(m)
    return mid


def generic_molecule(mid, label, page, fixed, variables, note=None):
    ev = [img(page)]
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, ev)
    m["names"]["original"] = F(label, label, ev)
    m["names"]["english"] = F()
    m["names"]["abbreviations"] = F()
    m["names"]["aliases"] = F()
    m["cas_numbers"] = F()
    m["representation_type"] = "generic_template"
    m["smiles"] = F()
    m["generic_template"] = {
        "fixed_structure": F(fixed, fixed, ev),
        "variables": [
            {
                "label": v,
                "attachment_points": F(f"source-labeled attachment for {v}", f"source-labeled attachment for {v}", ev),
                "definition": F(f"{v} as drawn in the source", f"{v} as drawn in the source", ev),
                "allowed_values": F(),
                "constraints": F(),
            }
            for v in variables
        ],
        "instance_mapping": [],
        "reason": F(
            "Generic R-group structure in source",
            note or "The source uses variable substituents; no single concrete molecule is asserted.",
            ev,
        ),
    }
    for key in ("chirality_status", "geometry_status", "representation_status"):
        m["stereochemistry"][key] = I(
            "unresolved" if key != "representation_status" else "not_encoded",
            ev,
            "The generic R-group drawing does not establish a unique stereochemical identity.",
            raw=label,
        )
    m["stereochemistry"]["atom_centers"] = F()
    m["stereochemistry"]["double_bonds"] = F()
    m["stereochemistry"]["axial"] = F()
    m["stereochemistry"]["mixture"] = F()
    m["stereochemistry"]["assignment"] = F()
    m["stereochemistry"]["isomeric_smiles"] = F()
    m["stereochemistry"]["identity_key"] = F()
    m["identity_note"] = F(
        "Generic structure retained as a template", "Generic structure retained as a template", ev
    )
    m["external_properties"] = []
    bundle["molecules"].append(m)
    return mid


# Named catalysts, additives, and solvents.
molecule("m_2czpn", "2CzPN", 1, name="1,2-bis(carbazol-9-yl)-4,5-dicyanobenzene", smiles="N#Cc1cc(n2c3ccccc3c3ccccc32)c(n2c3ccccc3c3ccccc32)cc1C#N", chirality="achiral", geometry="not_applicable")
molecule("m_tx", "TX", 2, name="thioxanthone", smiles="O=C1c2ccccc2Sc2ccccc21", chirality="achiral", geometry="not_applicable")
molecule("m_chcl3", "CHCl3", 2, name="chloroform", smiles="ClC(Cl)Cl", chirality="achiral", geometry="not_applicable")
molecule("m_dioxane", "1,4-dioxane", 2, name="1,4-dioxane", smiles="O1CCOCC1", chirality="achiral", geometry="not_applicable")
molecule("m_phme", "PhMe", 2, name="toluene", smiles="Cc1ccccc1", chirality="achiral", geometry="not_applicable")
molecule("m_dcm", "DCM", 2, name="dichloromethane", smiles="ClCCl", chirality="achiral", geometry="not_applicable")
molecule("m_mecn", "CH3CN", 2, name="acetonitrile", smiles="CC#N", chirality="achiral", geometry="not_applicable")
molecule("m_scotf3", "Sc(OTf)3", 2, name="scandium triflate", chirality="unresolved", geometry="not_applicable", representation_status="unresolved")
molecule("m_irf", "[Ir-F]", 2, name="(Ir[dF(CF3)ppy]2(dtbpy))PF6", chirality="unspecified", geometry="not_applicable", representation_status="not_encoded", note="The source drawing gives the complex connectivity but no metal-centered stereochemical descriptor; a complete candidate SMILES was not reconstructed.")

# Concrete compounds in the current paper.
for label in ("1a", "1b", "1c", "1d", "1e"):
    molecule(f"m_{label}", label, 2, chirality="achiral", geometry="not_applicable", representation_status="unresolved")
for label in ("2a", "2b", "2c", "2d", "2e", "3a"):
    molecule(f"m_{label}", f"rac-{label}", 2, chirality="racemic", geometry="not_applicable", representation_status="not_encoded", assignment="Relative stereochemistry is source-depicted; absolute enantiomer is not assigned.")
for label in ("4a", "4b", "4c", "4d", "4e", "4f"):
    molecule(f"m_{label}", label, 3, chirality="achiral", geometry="not_applicable", representation_status="unresolved")
for label in ("5a", "5b", "5c", "5d", "5e", "5f"):
    molecule(f"m_{label}", f"rac-{label}", 3, chirality="racemic", geometry="not_applicable", representation_status="not_encoded", assignment="Source-depicted relative stereochemistry; racemic product, with no absolute enantiomer assigned.")
molecule("m_6", "6", 4, chirality="achiral", geometry="not_applicable", representation_status="unresolved", axial={"allene": "terminal allene; no axial stereogenicity because one terminus is CH2"})
for label in ("7", "8"):
    molecule(f"m_{label}", f"rac-{label}", 4, chirality="racemic", geometry="not_applicable", representation_status="not_encoded", assignment="Source-depicted relative stereochemistry; absolute configuration is unassigned.")
for label in ("9a", "9b"):
    molecule(f"m_{label}", label, 4, chirality="achiral", geometry="not_applicable", representation_status="unresolved", axial={"allene": "terminal allene; no axial stereogenicity because one terminus is CH2"})
for label in ("10a", "11a", "11b"):
    molecule(f"m_{label}", f"rac-{label}", 4, chirality="racemic", geometry="not_applicable", representation_status="not_encoded", assignment="Source-depicted relative stereochemistry; absolute enantiomer is not assigned.")
molecule("m_10_entry2", "10 (entry 2, R1 = OPh; not observed)", 4, chirality="unresolved", geometry="not_applicable", representation_status="not_encoded", note="The table explicitly reports 0% in product column 10 for entry 2; no separate compound number or isolated structure is provided.")
molecule("m_12", "12", 4, chirality="unresolved", geometry="not_applicable", representation_status="not_encoded", note="Proposed dearomatized intermediate, not reported as isolated.")

# Generic structures from graphical abstract and scheme headers/prior work.
generic_molecule("m_graphical_substrate", "graphical-abstract generic allene substrate", 1, "arene/heteroarene tethered to an allene with X, Ar, R1 and R2 variables", ["X", "Ar", "R1", "R2"])
generic_molecule("m_graphical_product", "graphical-abstract generic polycyclic product", 1, "polycyclic product of the generic tethered-allene dearomatization", ["X", "Ar", "R1", "R2"])
generic_molecule("m_litA_substrate", "Scheme 1A generic naphthol derivative", 2, "naphthol-derived tethered alkene with EWG, R1 and R2 variables", ["EWG", "R1", "R2"])
generic_molecule("m_litA_product", "Scheme 1A generic racemic product", 2, "racemic dearomatized polycyclic product from Scheme 1A", ["EWG", "R1", "R2"])
generic_molecule("m_litB_substrate", "Scheme 1B generic indole derivative", 2, "indole-derived tethered alkene with X, R1 and R2 variables", ["X", "R1", "R2"])
generic_molecule("m_litB_product", "Scheme 1B generic racemic product", 2, "racemic spiroindoline product from Scheme 1B", ["X", "R1", "R2"])
generic_molecule("m_scheme2_substrate", "Scheme 2 generic naphthol 1", 2, "allyloxy naphthol derivative with R1/R2 variables", ["R1", "R2"])
generic_molecule("m_scheme2_product", "Scheme 2 generic rac-2", 2, "racemic dearomatized cycloadduct with R1/R2 variables", ["R1", "R2"])
generic_molecule("m_scheme4_substrate", "Scheme 4 generic indole 4", 3, "tethered-alkene indole with X, R1 and R2 variables", ["X", "R1", "R2"])
generic_molecule("m_scheme4_product", "Scheme 4 generic rac-5", 3, "racemic dearomatized spiroindoline with X, R1 and R2 variables", ["X", "R1", "R2"])
generic_molecule("m_scheme6_substrate", "Scheme 6 generic naphthol allene 9", 4, "allenyloxy naphthalene bearing COR1", ["R1"])
generic_molecule("m_scheme6_product10", "Scheme 6 generic rac-10", 4, "racemic aromatic cyclic acetal product family 10", ["R1"])
generic_molecule("m_scheme6_product11", "Scheme 6 generic rac-11", 4, "racemic dearomatized product family 11", ["R1"])


def P(pid, mid, role, evidence, *, amount=None, mol_percent=None, equivalents=None, reference=None, product_role=None):
    p = st.new_participant(pid, mid, role, evidence)
    if amount is not None:
        raw_amount = amount.get("raw_text") or f"{amount['number']} {amount['unit']}"
        p["charges"]["amount"] = F(raw_amount, amount, evidence)
    if mol_percent is not None:
        p["charges"]["mol_percent"] = F(f"{mol_percent} mol %", {"number": mol_percent, "unit": "mol%"}, evidence)
    if equivalents is not None:
        p["charges"]["equivalents"] = F(f"{equivalents} equiv", {"number": equivalents, "unit": "equiv"}, evidence)
    if reference is not None:
        p["equivalence_reference"] = I(reference, evidence, "Denominator assigned to the displayed substrate scale/loading basis.", raw="displayed limiting substrate")
    if product_role is not None:
        p["product_role"] = F(product_role, product_role, evidence)
    return p


def base_reaction(rid, evidence, *, record_type="specific_experiment", categories=None, location=None, reaction_name=None, transformation=None):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = F(categories, categories, evidence) if categories is not None else NA()
    r["location"] = F(location or rid, {"location": location or rid}, evidence)
    r["reaction_name"] = F(reaction_name, reaction_name, evidence) if reaction_name else F()
    r["transformation_type"] = I(transformation, evidence, "Transformation type inferred from the source caption/drawing.", raw=reaction_name or rid) if transformation else F()
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    return r


def set_common_conditions(r, evidence, *, substrate_mid=None, scale=None, solvent_mid=None, concentration=None, time=None, wavelength=None, temperature="rt", atmosphere="N2", light_source=None, operation="irradiation"):
    if scale is not None and substrate_mid:
        raw_scale = scale.get("raw_text") or f"{scale['number']} {scale['unit']}"
        r["reaction_scale"] = F(
            raw_scale + " scale",
            {"reference_molecule_id": substrate_mid, "amount": scale},
            evidence,
        )
    if solvent_mid:
        r["solvent_conditions"]["composition"] = F(solvent_mid, [solvent_mid], evidence)
    if concentration is not None:
        r["solvent_conditions"]["concentration"] = F(
            f"{concentration} M", {"number": concentration, "unit": "M"}, evidence
        )
        if substrate_mid:
            r["solvent_conditions"]["concentration_basis"] = I(
                substrate_mid, evidence, "Concentration basis assigned to the displayed substrate/reaction scale.", raw="scheme/table concentration"
            )
    participant_ids = [p["participant_id"] for p in r["participants"]]
    step = st.new_step("reaction", evidence)
    step["participant_ids"] = participant_ids
    step["operation"] = F(operation, operation, evidence)
    if temperature is not None:
        step["temperature"] = F(temperature, temperature, evidence)
    if time is not None:
        step["time"] = F(f"{time} h", {"number": time, "unit": "h"}, evidence)
    if atmosphere is not None:
        step["atmosphere"] = F(atmosphere, atmosphere, evidence)
    r["steps"] = [step]
    if wavelength is not None:
        r["special_conditions"]["wavelength"] = F(
            f"{wavelength} nm", {"number": wavelength, "unit": "nm"}, evidence
        )
    if light_source:
        r["special_conditions"]["light_source"] = F(light_source, light_source, evidence)
    elif wavelength is not None:
        r["special_conditions"]["light_source"] = F(f"{wavelength} nm light", "light irradiation", evidence)


def outcome(oid, kind, targets, evidence, *, raw, normalized=None, yield_type=None, basis=None, scope=None):
    o = st.new_outcome(oid, kind, targets, evidence)
    o["value"] = F(raw, normalized if normalized is not None else raw, evidence)
    if yield_type is not None:
        o["yield_type"] = F(yield_type, yield_type, evidence)
    if basis is not None:
        o["measurement_basis"] = F(basis, basis, evidence)
    if scope is not None:
        o["scope"] = F(scope, scope, evidence)
    return o


def class_cc(r, evidence, count=2, limitation=None):
    changes = []
    ids = []
    for i in range(count):
        cid = f"b{i+1}"
        ids.append(cid)
        change = {
            "change_id": cid,
            "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
            "before_order": 0,
            "after_order": 1,
            "change_type": "formed",
        }
        changes.append(D(change, evidence, "New C-C connection inferred from the source-drawn intramolecular cycloaddition; exact atom maps are unresolved.", raw="source-drawn cycloaddition connection"))
    label = {"category_number": 2, "code": "C-C", "supporting_change_ids": ids}
    r["classification"]["bond_changes"] = changes
    r["classification"]["labels"] = D([label], evidence, "Taxonomy label derived from the source-drawn new carbon-carbon connections.", raw="new C-C connection(s)")
    r["classification"]["coverage"] = I("partial", evidence, "Only source-visible new sigma connections are represented; exact atom correspondence, bond-order migration, hydrogens and coproducts are unresolved.", raw="source drawing")
    r["classification"]["limitations"] = [limitation or "Exact atom mapping and hydrogen/bond-order changes were not reconstructed from the dense source drawing."]
    r["classification"]["model_votes"] = []


def class_unresolved(r, evidence, message):
    r["classification"]["mapped_reaction_smiles"] = F()
    r["classification"]["bond_changes"] = []
    r["classification"]["labels"] = F()
    r["classification"]["coverage"] = I("unresolved", evidence, message, raw="source drawing")
    r["classification"]["limitations"] = [message]
    r["classification"]["model_votes"] = []


def class_na(r, evidence, message):
    r["classification"]["mapped_reaction_smiles"] = F()
    r["classification"]["bond_changes"] = []
    r["classification"]["labels"] = I([], evidence, message, raw="no observed transformation")
    r["classification"]["coverage"] = I("not_applicable", evidence, message, raw="no observed transformation")
    r["classification"]["limitations"] = [message]
    r["classification"]["model_votes"] = []


def add(r):
    if not r["review"]["unresolved"]:
        r["review"]["unresolved"] = ["No explicitly associated SI was supplied; detailed procedure, workup, purification, and complete atom mapping remain unavailable."]
    bundle["reactions"].append(r)
    return r["reaction_id"]


# Graphical abstract and prior-work arrows.
e = [img(1)]
r = base_reaction("graphical_abstract_generic", e, record_type="generic_scheme", location="page 1 graphical abstract", reaction_name="visible-light dearomative cycloaddition", transformation="generic dearomative cycloaddition of an allene-tethered arene/heteroarene")
r["participants"] = [P("ga_sub", "m_graphical_substrate", "reactant", e), P("ga_prod", "m_graphical_product", "product", e)]
class_unresolved(r, e, "The graphical abstract uses R-group variables and does not support complete bond mapping.")
add(r)

for suffix, smid, pmid, solvent, time_h, loading, title in [
    ("scheme1_literature_A", "m_litA_substrate", "m_litA_product", "m_dioxane", 18, 1, "Glorius (2018) previously reported naphthol dearomatization"),
    ("scheme1_literature_B", "m_litB_substrate", "m_litB_product", "m_chcl3", 48, 4, "You (2019) previously reported indole dearomatization"),
]:
    e = [img(2)]
    r = base_reaction(suffix, e, record_type="literature_example", location=f"page 2 Scheme 1 {'A' if suffix.endswith('A') else 'B'}", reaction_name=title, transformation="previously reported intramolecular dearomative cycloaddition")
    r["participants"] = [P("substrate", smid, "reactant", e), P("product", pmid, "product", e), P("irf", "m_irf", "catalyst", e, mol_percent=loading, reference=smid)]
    if suffix.endswith("B"):
        r["participants"] += [P("solvent_dcm", "m_dcm", "solvent", e), P("solvent_mecn", "m_mecn", "solvent", e)]
        r["solvent_conditions"]["composition"] = F("DCM/CH3CN", ["DCM", "CH3CN"], e)
        set_common_conditions(r, e, substrate_mid=smid, solvent_mid=None, time=time_h, wavelength=None, temperature=None, atmosphere=None, light_source="blue LED")
    else:
        r["participants"].append(P("solvent", solvent, "solvent", e))
        set_common_conditions(r, e, substrate_mid=smid, solvent_mid=solvent, time=time_h, wavelength=None, temperature=None, atmosphere=None, light_source="blue LED")
    class_cc(r, e)
    add(r)


# Table 1: ten independent optimization rows.
table_rows = [
    (1, "m_2czpn", 5, 17, "m_dioxane", "1.2:1", 455, 0.04, None),
    (2, "m_2czpn", 5, 17, "m_dioxane", "1.2:1", 405, 0.1, None),
    (3, "m_2czpn", 5, 17, "m_phme", "1.8:1", 455, 0.1, None),
    (4, "m_2czpn", 5, 17, "m_chcl3", "16:1", 455, 0.1, None),
    (5, "m_2czpn", 1, 14, "m_chcl3", "16:1", 455, 0.1, None),
    (6, "m_2czpn", 0.5, 14, "m_chcl3", "3:1", 455, 0.1, None),
    (7, "m_2czpn", 0.2, 14, "m_chcl3", "1:1.6", 455, 0.1, None),
    (8, "m_2czpn", 1, 17, "m_chcl3", "1:20", 455, 0.1, 0.1),
    (9, "m_tx", 5, 17, "m_chcl3", None, 455, 0.1, None),
    (10, None, None, 17, "m_dioxane", None, 455, 0.1, None),
]
ratio_values = {"1.2:1": [1.2, 1], "1.8:1": [1.8, 1], "16:1": [16, 1], "3:1": [3, 1], "1:1.6": [1, 1.6], "1:20": [1, 20]}
table_ids = {}
for entry, catalyst, loading, time_h, solvent, ratio, wavelength, conc, sc_eq in table_rows:
    ev = [img(2), txt(2)] + ([img(1), txt(1)] if entry == 1 else [])
    rid = f"table1_entry{entry:02d}"
    cats = ["optimization"] + (["substrate_scope"] if entry == 5 else [])
    r = base_reaction(rid, ev, categories=cats, location=f"page 2 Table 1 entry {entry}" + (" and Scheme 2 rac-2a" if entry == 5 else ""), reaction_name="organocatalytic dearomatization of naphthol 1a", transformation="visible-light intramolecular dearomative cycloaddition")
    r["general_scheme_ids"] = [f"{DOC}:scheme2_generic"] if entry == 5 else []
    r["participants"] = [P("substrate_1a", "m_1a", "reactant", ev, amount={"number": 0.2, "unit": "mmol"})]
    if ratio:
        r["participants"] += [P("product_2a", "m_2a", "product", ev, product_role="major/minor product according to reported ratio"), P("product_3a", "m_3a", "product", ev, product_role="major/minor product according to reported ratio")]
    if catalyst:
        r["participants"].append(P("photocatalyst", catalyst, "catalyst", ev, mol_percent=loading, reference="m_1a"))
    r["participants"].append(P("solvent", solvent, "solvent", ev))
    if sc_eq is not None:
        r["participants"].append(P("scandium_triflate", "m_scotf3", "additive", ev, equivalents=sc_eq, reference="m_1a"))
    set_common_conditions(r, ev, substrate_mid="m_1a", scale={"number": 0.2, "unit": "mmol"}, solvent_mid=solvent, concentration=conc, time=time_h, wavelength=wavelength, temperature="rt", atmosphere="N2", operation=f"irradiation at {wavelength} nm")
    if ratio:
        r["outcomes"] = [
            outcome("ratio_2a_3a", "selectivity", ["m_2a", "m_3a"], ev, raw=ratio, normalized={"ratio": ratio_values[ratio], "order": ["m_2a", "m_3a"]}, basis="1H crude NMR", scope={"metric": "product ratio", "order": ["2a", "3a"]}),
            outcome("conversion", "conversion", ["m_1a"], ev, raw="full conversion", normalized="full conversion", basis="no other proton signals in the 1H crude NMR"),
            outcome("combined_isolated_yield", "yield", ["m_2a", "m_3a"], ev, raw="near quantitative yield", normalized="near quantitative", yield_type="isolated mixture", basis="isolation of a mixture of 2a and 3a"),
        ]
        if entry == 5:
            r["outcomes"].append(outcome("isolated_yield_2a", "yield", ["m_2a"], ev, raw="94%", normalized={"number": 94, "unit": "%"}, yield_type="isolated", basis="isolated yield"))
            r["evidence_ids"] = list(dict.fromkeys(ev + [img(2)]))
        class_cc(r, ev)
    else:
        r["outcomes"] = [outcome("conversion", "conversion", ["m_1a"], ev, raw="No conversion of starting material", normalized="no conversion", basis=None)]
        class_na(r, ev, "Control experiment reports no conversion; no expected-product bond class is assigned.")
        r["review"]["unresolved"].append("The source reports no conversion but does not state the analytical method for the control result.")
    if entry == 8:
        r["notes"].append(F("0.1 equiv of Sc(OTf)3; the effect was not general to other substrates", "Sc(OTf)3 directs the model reaction toward 3a but was not general", ev))
    if entry == 1:
        r["notes"].append(F("0.04 M", "Entry 1 concentration differs from the 0.1 M table header condition", ev))
    table_ids[entry] = add(r)


# Proposed sequential triplet-energy mechanism described for the Table 1 system.
e = [img(2), txt(2)]
r = base_reaction("table1_sequential_triplet_mechanism", e, record_type="proposed_mechanism", location="page 2 prose below Table 1", reaction_name="two sequential triplet-energy excitation processes", transformation="1a to 3a by [2+2] cycloaddition, then 3a to 2a by excitation/rearrangement")
r["participants"] = [P("starting_material", "m_1a", "reactant", e), P("intermediate_product", "m_3a", "product", e), P("rearranged_product", "m_2a", "product", e)]
class_unresolved(r, e, "The prose proposes sequential photochemical events but complete atom-level mapping is not supplied in the main text.")
r["review"]["unresolved"].append("The complete mechanism is deferred to Supporting Information 5.1, which was not supplied.")
add(r)


# Scheme 2 generic arrow and four additional concrete scope entries (2a is Table 1 entry 5).
e = [img(2)]
r = base_reaction("scheme2_generic", e, record_type="generic_scheme", location="page 2 Scheme 2 generic arrow", reaction_name="substrate scope of naphthols", transformation="generic visible-light intramolecular dearomative cycloaddition")
r["participants"] = [P("generic_substrate", "m_scheme2_substrate", "reactant", e), P("generic_product", "m_scheme2_product", "product", e), P("photocatalyst", "m_2czpn", "catalyst", e, mol_percent=1, reference="m_scheme2_substrate"), P("solvent", "m_chcl3", "solvent", e)]
set_common_conditions(r, e, substrate_mid="m_scheme2_substrate", scale={"number": 0.2, "unit": "mmol"}, solvent_mid="m_chcl3", concentration=0.1, time=14, wavelength=455, temperature="rt", atmosphere="N2")
class_unresolved(r, e, "R-group scheme does not provide a single complete atom mapping.")
add(r)

scope2 = [("b", 62, 1, 14), ("c", 99, 1, 14), ("d", 76, 1, 14), ("e", 70, 2, 42)]
for suffix, yld, loading, time_h in scope2:
    e = [img(2)]
    r = base_reaction(f"scheme2_{suffix}", e, categories=["substrate_scope"], location=f"page 2 Scheme 2 rac-2{suffix}", reaction_name="dearomative cycloaddition of a naphthol derivative", transformation="visible-light intramolecular dearomative cycloaddition")
    r["general_scheme_ids"] = [f"{DOC}:scheme2_generic"]
    r["participants"] = [P("substrate", f"m_1{suffix}", "reactant", e, amount={"number": 0.2, "unit": "mmol"}), P("product", f"m_2{suffix}", "product", e, product_role="isolated product"), P("photocatalyst", "m_2czpn", "catalyst", e, mol_percent=loading, reference=f"m_1{suffix}"), P("solvent", "m_chcl3", "solvent", e)]
    set_common_conditions(r, e, substrate_mid=f"m_1{suffix}", scale={"number": 0.2, "unit": "mmol"}, solvent_mid="m_chcl3", concentration=0.1, time=time_h, wavelength=455, temperature="rt", atmosphere="N2")
    r["outcomes"] = [outcome("isolated_yield", "yield", [f"m_2{suffix}"], e, raw=f"{yld}%", normalized={"number": yld, "unit": "%"}, yield_type="isolated", basis="isolated yield")]
    class_cc(r, e)
    add(r)


# Scheme 3 scale-up and catalyst-reuse experiments.
e = [img(3), txt(3)]
r = base_reaction("scheme3_scaleup_1c", e, categories=["scale_up"], location="page 3 Scheme 3A", reaction_name="scale-up with recovery of 2CzPN", transformation="visible-light intramolecular dearomative cycloaddition")
r["general_scheme_ids"] = [f"{DOC}:scheme2_generic"]
r["participants"] = [P("substrate_1c", "m_1c", "reactant", e, amount={"number": 7.00, "unit": "mmol", "raw_text": "7.00 mmol"}), P("product_2c", "m_2c", "product", e, product_role="isolated product"), P("photocatalyst", "m_2czpn", "catalyst", e, mol_percent=1, reference="m_1c"), P("solvent", "m_chcl3", "solvent", e)]
set_common_conditions(r, e, substrate_mid="m_1c", scale={"number": 7.00, "unit": "mmol", "raw_text": "7.00 mmol"}, solvent_mid="m_chcl3", time=19, wavelength=455, temperature="rt", atmosphere="N2")
r["outcomes"] = [
    outcome("isolated_yield", "yield", ["m_2c"], e, raw="92%", normalized={"number": 92, "unit": "%"}, yield_type="isolated", basis="isolated product"),
    outcome("product_amount", "product_amount", ["m_2c"], e, raw="2.04 g", normalized={"number": 2.04, "unit": "g"}, basis="isolated product mass"),
    outcome("catalyst_recovery", "recovery", ["m_2czpn"], e, raw="88%", normalized={"number": 88, "unit": "%"}, basis="recovered by column chromatography"),
]
class_cc(r, e)
add(r)

e = [img(3), txt(3)]
r = base_reaction("scheme3_recycled_catalyst_1a", e, categories=["control"], location="page 3 Scheme 3B", reaction_name="reaction with recycled 2CzPN", transformation="visible-light intramolecular dearomative cycloaddition")
r["general_scheme_ids"] = [f"{DOC}:scheme2_generic"]
r["participants"] = [P("substrate_1a", "m_1a", "reactant", e), P("product_2a", "m_2a", "product", e, product_role="product"), P("recycled_photocatalyst", "m_2czpn", "catalyst", e), P("solvent", "m_chcl3", "solvent", e)]
set_common_conditions(r, e, substrate_mid="m_1a", solvent_mid="m_chcl3", time=14, wavelength=455, temperature="rt", atmosphere="N2")
r["outcomes"] = [outcome("yield", "yield", ["m_2a"], e, raw="90%", normalized={"number": 90, "unit": "%"})]
class_cc(r, e)
r["review"]["unresolved"].append("Scheme 3B does not print the recycled-catalyst loading or the yield measurement basis.")
add(r)


# Scheme 4 generic arrow and six current-paper indole scope entries.
e = [img(3)]
r = base_reaction("scheme4_generic", e, record_type="generic_scheme", location="page 3 Scheme 4 generic arrow", reaction_name="substrate scope of indoles", transformation="generic visible-light intramolecular dearomative cycloaddition")
r["participants"] = [P("generic_substrate", "m_scheme4_substrate", "reactant", e), P("generic_product", "m_scheme4_product", "product", e), P("photocatalyst", "m_2czpn", "catalyst", e, mol_percent=4, reference="m_scheme4_substrate"), P("solvent", "m_phme", "solvent", e)]
set_common_conditions(r, e, substrate_mid="m_scheme4_substrate", scale={"number": 0.05, "unit": "mmol"}, solvent_mid="m_phme", concentration=0.0125, time=24, wavelength=455, temperature="rt", atmosphere="N2")
class_unresolved(r, e, "R-group scheme does not provide a single complete atom mapping.")
add(r)

scheme4 = [("a", 95, 24), ("b", 99, 24), ("c", 99, 19), ("d", 76, 24), ("e", 51, 24), ("f", 76, 24)]
for suffix, yld, time_h in scheme4:
    e = [img(3)]
    r = base_reaction(f"scheme4_{suffix}", e, categories=["substrate_scope"], location=f"page 3 Scheme 4 rac-5{suffix}", reaction_name="dearomatization of an indole derivative", transformation="visible-light intramolecular dearomative cycloaddition")
    r["general_scheme_ids"] = [f"{DOC}:scheme4_generic"]
    r["participants"] = [P("substrate", f"m_4{suffix}", "reactant", e, amount={"number": 0.05, "unit": "mmol"}), P("product", f"m_5{suffix}", "product", e, product_role="isolated product"), P("photocatalyst", "m_2czpn", "catalyst", e, mol_percent=4, reference=f"m_4{suffix}"), P("solvent", "m_phme", "solvent", e)]
    set_common_conditions(r, e, substrate_mid=f"m_4{suffix}", scale={"number": 0.05, "unit": "mmol"}, solvent_mid="m_phme", concentration=0.0125, time=time_h, wavelength=455, temperature="rt", atmosphere="N2")
    r["outcomes"] = [outcome("isolated_yield", "yield", [f"m_5{suffix}"], e, raw=f"{yld}%", normalized={"number": yld, "unit": "%"}, yield_type="isolated", basis="isolated yield")]
    class_cc(r, e)
    add(r)


# Scheme 5: allene-substituted indole, duplicated as the lower arrow of Scheme 1D.
e = [img(4), txt(4), img(2)]
r = base_reaction("scheme5_allene_indole6", e, categories=["substrate_scope"], location="page 4 Scheme 5; overview duplicate on page 2 Scheme 1D", reaction_name="formation of methylenecyclobutane products from allene-substituted indole 6", transformation="visible-light intramolecular dearomative [2+2] cycloaddition of an allene")
r["participants"] = [P("substrate_6", "m_6", "reactant", e, amount={"number": 0.05, "unit": "mmol"}), P("product_7", "m_7", "product", e, product_role="major product"), P("product_8", "m_8", "product", e, product_role="minor product"), P("photocatalyst", "m_2czpn", "catalyst", e, mol_percent=4, reference="m_6"), P("solvent", "m_phme", "solvent", e)]
set_common_conditions(r, e, substrate_mid="m_6", scale={"number": 0.05, "unit": "mmol"}, solvent_mid="m_phme", concentration=0.0125, time=19, wavelength=455, temperature="rt", atmosphere="N2", light_source="blue LED")
r["outcomes"] = [
    outcome("combined_isolated_yield", "yield", ["m_7", "m_8"], e, raw="99%", normalized={"number": 99, "unit": "%"}, yield_type="isolated combined products", basis="isolated yield"),
    outcome("ratio_7_8", "selectivity", ["m_7", "m_8"], e, raw="5.3:1", normalized={"ratio": [5.3, 1], "order": ["m_7", "m_8"]}, scope={"metric": "product ratio", "order": ["7", "8"]}),
    outcome("conversion", "conversion", ["m_6"], e, raw="full conversion", normalized="full conversion"),
]
class_cc(r, e)
r["review"]["unresolved"].append("The product-ratio and full-conversion measurement methods are not reported in the supplied main text.")
add(r)


# Scheme 6 generic arrow, two concrete experiments, and proposed intermediate 12.
e = [img(4)]
r = base_reaction("scheme6_generic", e, record_type="generic_scheme", location="page 4 Scheme 6 generic arrow", reaction_name="dearomative cycloaddition of allene-substituted naphthol derivatives", transformation="generic visible-light photochemistry of naphthol allenes")
r["participants"] = [P("generic_substrate", "m_scheme6_substrate", "reactant", e), P("generic_product10", "m_scheme6_product10", "product", e), P("generic_product11", "m_scheme6_product11", "product", e), P("photocatalyst", "m_2czpn", "catalyst", e, mol_percent=2, reference="m_scheme6_substrate"), P("solvent", "m_chcl3", "solvent", e)]
set_common_conditions(r, e, substrate_mid="m_scheme6_substrate", scale={"number": 0.1, "unit": "mmol"}, solvent_mid="m_chcl3", concentration=0.05, time=16, wavelength=455, temperature="rt", atmosphere="N2")
class_unresolved(r, e, "The branched generic scheme has R1 variables and product-dependent rearrangement; complete bond changes are unresolved.")
add(r)

e = [img(4), txt(4)]
r = base_reaction("scheme6_entry1_9a", e, categories=["substrate_scope"], location="page 4 Scheme 6 entry 1", reaction_name="visible-light photochemistry of naphthol ketone 9a", transformation="formation of aromatic cyclic acetal 10a and dearomatized dihydrofuran 11a")
r["general_scheme_ids"] = [f"{DOC}:scheme6_generic"]
r["participants"] = [P("substrate_9a", "m_9a", "reactant", e, amount={"number": 0.1, "unit": "mmol"}), P("product_10a", "m_10a", "product", e, product_role="main aromatic cyclic acetal product"), P("product_11a", "m_11a", "product", e, product_role="minor 1,4-cycloaddition product"), P("photocatalyst", "m_2czpn", "catalyst", e, mol_percent=2, reference="m_9a"), P("solvent", "m_chcl3", "solvent", e)]
set_common_conditions(r, e, substrate_mid="m_9a", scale={"number": 0.1, "unit": "mmol"}, solvent_mid="m_chcl3", concentration=0.05, time=16, wavelength=455, temperature="rt", atmosphere="N2")
r["outcomes"] = [
    outcome("yield_10a", "yield", ["m_10a"], e, raw="61%", normalized={"number": 61, "unit": "%"}, yield_type="isolated", basis="isolated yield"),
    outcome("yield_11a", "yield", ["m_11a"], e, raw="15%", normalized={"number": 15, "unit": "%"}, yield_type="isolated", basis="isolated yield"),
]
class_cc(r, e, limitation="Both products arise from source-described cycloaddition pathways, but product-specific atom mapping and the additional acetal/rearrangement bond changes are unresolved.")
r["notes"].append(F("10a was confirmed by single-crystal X-ray analysis", "10a structure confirmed by single-crystal X-ray analysis", e))
r["review"]["unresolved"].append("The two product pathways have different net rearrangements; only the clearly shared cycloaddition C-C classification is retained at partial coverage.")
add(r)

e = [img(4), txt(4), img(2)]
r = base_reaction("scheme6_entry2_9b", e, categories=["substrate_scope"], location="page 4 Scheme 6 entry 2; overview duplicate on page 2 Scheme 1D", reaction_name="visible-light photochemistry of naphthol phenyl ester 9b", transformation="selective dearomative cycloaddition to 11b")
r["general_scheme_ids"] = [f"{DOC}:scheme6_generic"]
r["participants"] = [P("substrate_9b", "m_9b", "reactant", e, amount={"number": 0.1, "unit": "mmol"}), P("unobserved_product10", "m_10_entry2", "product", e, product_role="product family 10, not observed"), P("product_11b", "m_11b", "product", e, product_role="only observed product"), P("photocatalyst", "m_2czpn", "catalyst", e, mol_percent=2, reference="m_9b"), P("solvent", "m_chcl3", "solvent", e)]
set_common_conditions(r, e, substrate_mid="m_9b", scale={"number": 0.1, "unit": "mmol"}, solvent_mid="m_chcl3", concentration=0.05, time=16, wavelength=455, temperature="rt", atmosphere="N2", light_source="blue LED")
r["outcomes"] = [
    outcome("yield_product10", "yield", ["m_10_entry2"], e, raw="0%", normalized={"number": 0, "unit": "%"}, yield_type="isolated", basis="isolated yield"),
    outcome("yield_11b", "yield", ["m_11b"], e, raw="81%", normalized={"number": 81, "unit": "%"}, yield_type="isolated", basis="isolated yield"),
]
class_cc(r, e)
r["notes"].append(F("only the dearomatized product 11b was observed", "11b was the only observed product", e))
add(r)

e = [img(4), txt(4)]
r = base_reaction("scheme6_proposed_intermediate12", e, record_type="proposed_mechanism", location="page 4 Scheme 6 and accompanying prose", reaction_name="proposed formation of 10a via dearomatized intermediate 12", transformation="proposed cycloaddition followed by 1,3-allylic transposition")
r["participants"] = [P("substrate_9a", "m_9a", "reactant", e), P("intermediate_12", "m_12", "product", e, product_role="proposed intermediate"), P("product_10a", "m_10a", "product", e, product_role="observed product")]
class_unresolved(r, e, "Intermediate 12 and the 1,3-allylic transposition are proposed rather than isolated; exact bond mapping is unresolved.")
r["notes"].append(F("We believe that 10a is formed via dearomatized intermediate 12 that upon 1,3-allylic transposition yields product 10a.", "Author-proposed pathway, not a separately demonstrated experiment", e))
add(r)


# Non-synthetic observations and page-level source context.
bundle["source_context"] = [
    {
        "context_id": f"{DOC}:photosensitizer_properties",
        "kind": "photophysical_property",
        "evidence_ids": [img(1), txt(1), img(2)],
        "observations": [
            {"subject_molecule_id": "m_2czpn", "metric": "triplet energy", "raw_value": "60.6 kcal/mol (T1 = 2.63 eV)"}
        ],
        "note": "Contextual catalyst property, not a synthetic experiment.",
    },
    {
        "context_id": f"{DOC}:table2_solubility",
        "kind": "solubility_table",
        "evidence_ids": [img(3), txt(3)],
        "observations": {
            "units": "M",
            "rows": [
                {"solvent": "PhMe", "[Ir-F]": "7.0e-5", "2CzPN": "1.3e-2"},
                {"solvent": "1,4-dioxane", "[Ir-F]": "2.2e-4", "2CzPN": "1.6e-2"},
                {"solvent": "DCM", "[Ir-F]": "5.6e-3", "2CzPN": "8.3e-2"},
                {"solvent": "methyl tert-butyl ether", "[Ir-F]": "7.0e-5", "2CzPN": "9.2e-4"},
                {"solvent": "DMSO", "[Ir-F]": "1.6e-1", "2CzPN": "1.6e-2"},
            ],
        },
        "note": "Physical-property measurements are not reaction experiments; the 2CzPN method is deferred to unsupplied SI 3.4.",
    },
    {
        "context_id": f"{DOC}:xray_and_accession",
        "kind": "structure_characterization",
        "evidence_ids": [img(4), txt(4)],
        "observations": {"10a": "confirmed by single-crystal X-ray analysis", "ccdc": ["1998653-1998654", "1998663", "2000974"]},
        "note": "Characterization/accession information; no standalone synthesis record created.",
    },
    {
        "context_id": f"{DOC}:pages5_6_bibliography",
        "kind": "bibliography_only",
        "evidence_ids": [img(5), img(6)],
        "note": "Pages 5-6 contain cited literature only and no current-paper experimental arrow/table row.",
    },
]


# Reproducible raw transcription and literal mention search.
page_records = {str(i): [] for i in range(1, 7)}
for rxn in bundle["reactions"]:
    for eid in rxn["evidence_ids"]:
        if ":page" in eid and eid.endswith(":image"):
            page = eid.split(":page", 1)[1].split(":", 1)[0]
            if rxn["reaction_id"] not in page_records[page]:
                page_records[page].append(rxn["reaction_id"])
raw = {
    "document_id": DOC,
    "title": TITLE,
    "attempted_pages": [1, 2, 3, 4, 5, 6],
    "records": [
        {
            "reaction_id": r["reaction_id"],
            "record_type": r["record_type"],
            "location": r["location"]["raw_value"],
            "evidence_ids": r["evidence_ids"],
            "decision": "retained as one source-derived record; repeated overview appearances are linked as evidence rather than duplicated",
        }
        for r in bundle["reactions"]
    ],
    "page_record_ids": page_records,
    "limitations": [
        "No SI file was supplied or associated.",
        "Dense structures were retained as source labels/templates when a safe complete stereochemical SMILES was not reconstructed.",
        "Bracketed [Ir-F] comparison yields in Schemes 2 and 4 are cited prior-work values and were not made into current-paper experiments.",
    ],
}
write(RUN / "raw_model" / "pages_0001_0006_extraction.json", raw)
write(RUN / "raw_model" / "stereo_checks.json", stereo_checks)
write(
    RUN / "raw_model" / "source_mentions.json",
    st.source_mentions(dataset, ["1a", "2a", "3a", "4b", "6", "7", "8", "9a", "9b", "10a", "11a", "11b", "12", "2CzPN"]),
)
write(RUN / "submitted.json", bundle)


# Coverage contract: exactly original PDF pages 1-6.
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 1,
            "state": "attempted",
            "regions": [
                {"region": "title/abstract and graphical abstract", "finding": "title metadata and one generic current-study dearomative-cycloaddition scheme", "evidence_ids": [img(1)]},
                {"region": "lower-right prose", "finding": "Table 1 entry 1 begins here; merged with the page-2 table row rather than duplicated", "evidence_ids": [img(1), txt(1)]},
            ],
            "record_ids": page_records["1"],
            "unresolved": ["Graphical-abstract R-group structures do not support a unique molecule or complete atom map."],
        },
        {
            "page": 2,
            "state": "attempted",
            "regions": [
                {"region": "Scheme 1A-B", "finding": "two cited prior-work generic reaction arrows, classified as literature examples", "evidence_ids": [img(2)]},
                {"region": "Scheme 1C", "finding": "[Ir-F] catalyst structure/property context, not a reaction", "evidence_ids": [img(2)]},
                {"region": "Scheme 1D", "finding": "two current-work overview arrows; linked to Scheme 6 entry 2 and Scheme 5 instead of duplicated", "evidence_ids": [img(2)]},
                {"region": "Table 1", "finding": "ten independent optimization/control rows and all footnotes transcribed", "evidence_ids": [img(2), txt(2)]},
                {"region": "Scheme 2", "finding": "generic arrow plus five scope products; rac-2a merged with identical Table 1 entry 5 and four additional scope experiments retained", "evidence_ids": [img(2)]},
                {"region": "prose around Table 1", "finding": "sequential excitation mechanism and control interpretation retained separately from experiments", "evidence_ids": [img(2), txt(2)]},
            ],
            "record_ids": page_records["2"],
            "unresolved": ["No SI 5.1 was supplied for the complete mechanism.", "Dense product structures and exact atom maps remain unresolved."],
        },
        {
            "page": 3,
            "state": "attempted",
            "regions": [
                {"region": "Scheme 3A-B", "finding": "one 7.00 mmol scale-up/catalyst-recovery experiment and one recycled-catalyst experiment", "evidence_ids": [img(3), txt(3)]},
                {"region": "Table 2", "finding": "solubility measurements retained in source_context, not as reactions", "evidence_ids": [img(3), txt(3)]},
                {"region": "Scheme 4", "finding": "generic indole scope arrow plus six independent current-paper product rows; bracketed [Ir-F] yields kept as cited comparisons", "evidence_ids": [img(3)]},
                {"region": "lower-right prose", "finding": "allene motivation and cited Arai/Ohkuma precedent; no new current-paper arrow on this page", "evidence_ids": [img(3), txt(3)]},
            ],
            "record_ids": page_records["3"],
            "unresolved": ["Scheme 3B does not print recycled-catalyst loading or yield measurement basis.", "The SI-only indole optimization is unavailable."],
        },
        {
            "page": 4,
            "state": "attempted",
            "regions": [
                {"region": "Scheme 5 and preceding prose", "finding": "allene-indole 6 gives products 7 and 8 in one experiment; combined isolated yield, ratio and full conversion retained", "evidence_ids": [img(4), txt(4)]},
                {"region": "Scheme 6 table and arrows", "finding": "generic arrow plus two independent naphthol-allene experiments, including both product columns and the 0% result", "evidence_ids": [img(4), txt(4)]},
                {"region": "Scheme 6 intermediate 12/prose", "finding": "author-proposed intermediate and 1,3-allylic transposition recorded as proposed mechanism", "evidence_ids": [img(4), txt(4)]},
                {"region": "Associated Content and accession codes", "finding": "SI advertised but not supplied; X-ray/accession data retained as context", "evidence_ids": [img(4), txt(4)]},
            ],
            "record_ids": page_records["4"],
            "unresolved": ["Scheme 5 ratio/conversion methods are unreported.", "Scheme 6 product-specific atom maps and rearrangement bond changes remain partial."],
        },
        {
            "page": 5,
            "state": "attempted",
            "regions": [{"region": "both columns", "finding": "bibliography only; no current-paper reaction arrow, table row, or standalone experimental record", "evidence_ids": [img(5)]}],
            "record_ids": [],
            "unresolved": [],
        },
        {
            "page": 6,
            "state": "attempted",
            "regions": [{"region": "both columns", "finding": "bibliography only, including prior allene photochemistry references; no current-paper experiment", "evidence_ids": [img(6)]}],
            "record_ids": [],
            "unresolved": [],
        },
    ],
}
write(RUN / "coverage" / "pages_0001_0006.json", coverage)


# Mark mechanical page tasks attempted after visual and textual review.
tasks = json.loads((RUN / "tasks.json").read_text())
for task in tasks:
    task["state"] = "attempted"
    task["reason"] = "Original rendered page and extracted text inspected; findings and record links are documented in coverage/pages_0001_0006.json."
write(RUN / "tasks.json", tasks)

if "standard" in dataset:
    write(RUN / "raw_model" / "audit.json", rt.dataset_audit(dataset))

print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]), "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in bundle["reactions"]), "source_context": len(bundle["source_context"])}, ensure_ascii=False, indent=2))
