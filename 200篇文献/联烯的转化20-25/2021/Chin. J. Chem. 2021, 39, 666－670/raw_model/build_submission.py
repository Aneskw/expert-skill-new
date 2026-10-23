#!/usr/bin/env python3
"""Build the source-derived v3.1 submission for doc_1149787c4988fff5.

The script contains extraction decisions only.  It imports the unchanged project
helpers and writes the machine submission, visual inventory and required page
coverage files beside the prepared dataset.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DOC = "doc_1149787c4988fff5"
P = {i: f"{DOC}:page{i}:image" for i in range(1, 6)}
PT = {i: f"{DOC}:page{i}:text" for i in range(1, 6)}


def save(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def F(raw=None, normalized=None, ev=(), method="original_extraction", derivation=None,
      status=None):
    return st.fact(raw, normalized, ev, method=method, derivation=derivation,
                   status=status)


def inferred(normalized, ev, why, raw=None):
    return F(raw, normalized, ev, method="inference", derivation=why)


def derived(normalized, ev, why, raw=None):
    return F(raw, normalized, ev, method="structure_derivation", derivation=why)


def conflict(alternatives, ev):
    result = F([a[0] for a in alternatives], None, ev, status="conflict")
    result["alternatives"] = [F(raw, normalized, aev) for raw, normalized, aev in alternatives]
    return result


def mixed_conflict(alternatives, ev):
    """Conflict whose alternatives may mix direct source values and calculations."""
    result = F([a[0] for a in alternatives], None, ev, status="conflict")
    result["alternatives"] = []
    for raw, normalized, aev, method, why in alternatives:
        result["alternatives"].append(F(raw, normalized, aev, method=method, derivation=why))
    return result


def qty(number, unit, raw=None):
    value = {"number": number, "unit": unit}
    if raw is not None:
        value["raw_text"] = raw
    return value


bundle = json.loads((RUN / "submitted.json").read_text())
bundle["molecules"] = []
bundle["reactions"] = []
bundle["source_context"] = []
STEREO_DIAGNOSTICS = {}


def add_molecule(mid, label, name=None, ev=(P[3],), smiles=None,
                 chirality="unresolved", geometry="unresolved",
                 representation="unresolved", assignment=None,
                 axial=None, generic=None, note=None):
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, ev)
    m["occurrences"][0]["evidence_ids"] = list(ev)
    m["names"]["original"] = F(name or label, name or label, ev)
    m["names"]["english"] = F(name or label, name or label, ev)
    m["names"]["abbreviations"] = F(label, [label], ev)
    m["cas_numbers"] = F()
    m["identity_note"] = inferred(
        note or "Identity is retained at the source label/drawing level; no external identity enrichment was used.",
        ev, "Machine identity note based on the supplied article drawing and label.")
    m["stereochemistry"]["chirality_status"] = inferred(
        chirality, ev, "Status assigned conservatively from visible stereochemical markings; absence of a named configuration is not treated as achiral or racemic.")
    m["stereochemistry"]["geometry_status"] = inferred(
        geometry, ev, "Status assigned conservatively from the visible drawing and labels.")
    m["stereochemistry"]["representation_status"] = inferred(
        representation, ev, "Compared the registry representation with the source drawing; unencoded axial/relative information is retained separately.")
    if assignment is not None:
        m["stereochemistry"]["assignment"] = F(assignment, assignment, ev)
    if axial is not None:
        m["stereochemistry"]["axial"] = F(axial, axial, ev)
    if generic is not None:
        m["representation_type"] = "generic_template"
        m["generic_template"] = {
            "fixed_structure": F(generic["fixed"], generic["fixed"], ev),
            "variables": [
                {
                    "label": var,
                    "attachment_points": F(desc, desc, ev),
                    "definition": F(desc, desc, ev),
                    "allowed_values": F(),
                    "constraints": F(),
                }
                for var, desc in generic.get("variables", [])
            ],
            "instance_mapping": [],
            "reason": F("Generic R-group drawing", "Source uses an R-group template that is not a single concrete molecule.", ev),
        }
    elif smiles is not None:
        m["representation_type"] = "smiles"
        m["smiles"] = derived(smiles, ev, "Candidate connectivity transcribed from the source structure/name; parsing is not identity verification.", raw="source structure/name")
        diag = st.stereo_inspect(smiles)
        STEREO_DIAGNOSTICS[mid] = diag
        m["stereochemistry"]["isomeric_smiles"] = derived(
            diag["canonical_isomeric_smiles"], ev, "RDKit canonical isomeric SMILES generated from the machine candidate.")
        m["stereochemistry"]["identity_key"] = derived(
            diag["stereo_graph_key"], ev, "RDKit stereo graph key generated from the machine candidate.")
    else:
        m["representation_type"] = "unresolved"
    bundle["molecules"].append(m)
    return mid


def add_simple(mid, label, smiles, ev):
    if mid == "m_DABCO":
        return add_molecule(mid, label, label, ev, smiles=smiles,
                            chirality="unresolved", geometry="not_applicable",
                            representation="partial",
                            note="Named reagent candidate; toolkit exposes an unspecified potential stereo element, so no achiral/fully encoded claim is made.")
    return add_molecule(mid, label, label, ev, smiles=smiles,
                        chirality="achiral", geometry="not_applicable",
                        representation="encoded",
                        note="Simple named reagent/solvent candidate; graph was checked but source identity remains machine-level.")


# Generic article-level structures.
add_molecule("m_gen_1", "1", "allylamines 1 (generic)", (P[1], P[2], P[4]),
             chirality="unresolved", geometry="unspecified", representation="not_encoded",
             generic={"fixed": "R1-NH-CH(R2)-CH=CH2 as drawn",
                      "variables": [("R1", "N-substituent"), ("R2", "substituent on the allylic carbon")]})
add_molecule("m_gen_2", "2", "electron-deficient allene diesters 2 (generic)", (P[1], P[2], P[4]),
             chirality="unspecified", geometry="not_applicable", representation="not_encoded",
             axial="Allene axis is present; no R_a/S_a assignment is reported.",
             generic={"fixed": "R4-CH=C=CH-CO2R3 as drawn",
                      "variables": [("R3", "ester alkyl group"), ("R4", "terminal allene substituent shown in red")]})
add_molecule("m_gen_3", "3", "3-azabicyclo[3.1.0]hexane products 3 (generic)", (P[1], P[2], P[3], P[4]),
             chirality="partially_specified", geometry="unspecified", representation="not_encoded",
             assignment="Relative bonds are drawn, but no generic absolute configuration is named.",
             generic={"fixed": "3-azabicyclo[3.1.0]hexane product framework with ester-substituted alkene",
                      "variables": [("R1", "N-substituent"), ("R2", "cyclopropane substituent"),
                                    ("R3", "ester alkyl group"), ("R4", "alkene substituent")]})

# Main specific organic inputs with cautious graph candidates.
add_molecule("m_1a", "1a", "diallylamine 1a", (P[2], P[3], P[4]),
             smiles="C=CCNCC=C", chirality="achiral", geometry="not_applicable", representation="encoded")
add_molecule("m_2a", "2a", "dimethyl allene-1,3-dicarboxylate 2a", (P[2], P[3], P[4]),
             smiles="COC(=O)C=C=CC(=O)OC", chirality="unspecified", geometry="not_applicable",
             representation="partial", axial="Potential allene axial stereochemistry is not assigned in the article.")
add_molecule("m_2b", "2b", "diethyl allene-1,3-dicarboxylate 2b", (P[3],),
             smiles="CCOC(=O)C=C=CC(=O)OCC", chirality="unspecified", geometry="not_applicable",
             representation="partial", axial="Potential allene axial stereochemistry is not assigned in the article.")

# Allylamines 1b-1r are retained by label and visible substituent description.
allyl_names = {
    "1b": "N-phenyl allylamine 1b", "1c": "N-(4-methoxybenzyl) allylamine 1c",
    "1d": "N-(cyclohexylmethyl) allylamine 1d", "1e": "N-(4-bromobenzyl) allylamine 1e",
    "1f": "N-(thiophen-2-ylmethyl) allylamine 1f", "1g": "N-(naphthylmethyl) allylamine 1g",
    "1h": "benzyl-substituted allylamine 1h (R2 = methyl as drawn)",
    "1i": "benzyl-substituted allylamine 1i (R2 = phenyl as drawn)",
    "1j": "benzyl-substituted aryl allylamine 1j", "1k": "benzyl-substituted 4-methoxyaryl allylamine 1k",
    "1l": "benzyl-substituted methoxyaryl allylamine 1l", "1m": "benzyl-substituted methoxyaryl allylamine 1m",
    "1n": "benzyl-substituted chloroaryl allylamine 1n", "1o": "benzyl-substituted chloroaryl allylamine 1o",
    "1p": "benzyl-substituted chloroaryl allylamine 1p", "1q": "benzyl-substituted bromoaryl allylamine 1q",
    "1r": "benzyl-substituted bromoaryl allylamine 1r",
}
for label, name in allyl_names.items():
    ev = (P[3],) if label != "1k" else (P[3], PT[3])
    add_molecule(f"m_{label}", label, name, ev, chirality="unspecified",
                 geometry="not_applicable", representation="not_encoded",
                 note="Specific compact drawing was read at label/substituent level; a full stereochemically safe SMILES was not asserted.")

product_yields = {
    "3aa": 80, "3ba": 85, "3ca": 73, "3da": 68, "3ea": 42,
    "3fa": 65, "3ga": 46, "3ha": 67, "3ia": 62, "3ja": 65,
    "3ka": 54, "3la": 61, "3ma": 53, "3na": 76, "3oa": 63,
    "3pa": 60, "3qa": 56, "3ra": 57, "3bb": 75,
}
for label in product_yields:
    ev = [P[3]]
    if label == "3aa": ev += [P[2], P[4]]
    assignment = "Scheme drawings contain relative wedge/bond information; absolute configuration is not named."
    if label == "3ra":
        assignment += " The article states that the structure was established by single-crystal X-ray diffraction."
    add_molecule(f"m_{label}", label, f"3-azabicyclo[3.1.0]hexane product {label}", tuple(ev),
                 chirality="partially_specified", geometry="unspecified", representation="not_encoded",
                 assignment=assignment,
                 note="Compact bicyclic product drawing retained without forcing an isomeric SMILES; ring/relative stereochemistry requires expert review.")
add_molecule("m_4aa", "4aa", "brominated byproduct 4aa", (P[2],),
             chirality="partially_specified", geometry="unspecified", representation="not_encoded",
             assignment="Relative bonds are drawn; no absolute configuration is named.")
add_molecule("m_5", "5", "carboxylic acid derivative 5", (P[3],),
             chirality="partially_specified", geometry="not_applicable", representation="not_encoded",
             assignment="Figure states relative configuration; structure established by single-crystal X-ray diffraction, but absolute configuration is not stated.")

# Prior-work Scheme 1 structures are intentionally not conflated with current-work labels.
for label, name in [("lit_I", "Baskaran ketenimine I"), ("lit_allylamine", "generic allylamine in Baskaran's work"),
                    ("lit_II", "Baskaran intermediate II"), ("lit_product", "Baskaran cyclopropane-fused cyclic amidine product")]:
    add_molecule(f"m_{label}", label, name, (P[2],), chirality="unresolved",
                 geometry="unresolved", representation="not_encoded")

# Named agents, solvents and operation materials.
simple = {
    "CuBr2": "[Cu+2].[Br-].[Br-]", "CuI": "[Cu+].[I-]",
    "FeCl3": "[Fe+3].[Cl-].[Cl-].[Cl-]", "CuOAc2": "[Cu+2].CC(=O)[O-].CC(=O)[O-]",
    "THF": "C1CCOC1", "dioxane": "O1CCOCC1", "DCM": "ClCCl", "DCE": "ClCCCl",
    "diethyl_ether": "CCOCC", "toluene": "Cc1ccccc1", "MeCN": "CC#N",
    "DBU": "C1CCC2=NCCCN2CC1", "DMAP": "CN(C)c1ccncc1", "DABCO": "C1CN2CCN1CC2",
    "dibutylamine": "CCCCNCCCC", "NaOH": "[Na+].[OH-]", "Cs2CO3": "[Cs+].[Cs+].[O-]C([O-])=O",
    "KOH": "[K+].[OH-]", "water": "O", "HCl": "Cl", "ethyl_acetate": "CCOC(C)=O",
}
display = {"CuOAc2": "Cu(OAc)2", "dioxane": "1,4-dioxane", "DCM": "CH2Cl2",
           "DCE": "1,2-dichloroethane (DCE)", "MeCN": "CH3CN", "Cs2CO3": "Cs2CO3"}
for key, smiles in simple.items():
    add_simple(f"m_{key}", display.get(key, key.replace("_", " ")), smiles, (P[2], P[3], P[4]))
add_molecule("m_petroleum_ether", "petroleum ether", "petroleum ether", (P[4],),
             chirality="not_applicable", geometry="not_applicable", representation="unresolved",
             note="Named mixture; composition is not reported and no single-component identity is invented.")
add_molecule("m_silica", "silica gel", "silica gel", (P[4],), chirality="not_applicable",
             geometry="not_applicable", representation="unresolved")


def participant(pid, mid, role, ev, stage="reaction", amount=None, mass=None,
                volume=None, equiv=None, molpct=None, ref=None, product_role=None):
    p = st.new_participant(pid, mid, role, ev, stage)
    if amount is not None:
        p["charges"]["amount"] = F(f"{amount} mmol", qty(amount, "mmol"), ev)
    if mass is not None:
        number, unit = mass
        p["charges"]["mass"] = F(f"{number} {unit}", qty(number, unit), ev)
    if volume is not None:
        p["charges"]["volume"] = F(f"{volume} mL", qty(volume, "mL"), ev)
    if equiv is not None:
        p["charges"]["equivalents"] = F(f"{equiv} equiv", qty(equiv, "equiv"), ev)
    if molpct is not None:
        p["charges"]["mol_percent"] = F(f"{molpct} mol%", qty(molpct, "mol%"), ev)
    if ref is not None:
        p["equivalence_reference"] = F(f"relative to {ref}", ref, ev)
    if product_role is not None:
        p["product_role"] = F(product_role, product_role, ev)
    return p


def base_reaction(rid, ev, record_type, categories, location, name, transformation):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = F(categories, categories, ev)
    r["location"] = F(location, {"description": location}, ev)
    r["reaction_name"] = F(name, name, ev)
    r["transformation_type"] = inferred(transformation, ev, "Transformation description summarized from the source title/scheme.")
    r["review"]["unresolved"] = []
    r["completeness"] = "partial"
    return r


def set_classification(r, ev, observed=True, brominated=False, unresolved=False):
    c = r["classification"]
    if unresolved:
        c["coverage"] = inferred("unresolved", ev, "Exact atom correspondence cannot be assigned safely from the compact drawing.")
        c["limitations"] = ["No atom mapping is asserted; exact bond changes remain unresolved."]
        return
    if not observed:
        c["labels"] = inferred([], ev, "No observed product transformation is represented; expected chemistry is not classified.")
        c["coverage"] = inferred("not_applicable", ev, "Table reports no observed product for this control/condition.")
        c["limitations"] = ["No inferred labels are attached to a non-observed transformation."]
        return
    changes = [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "N"}, {"atom_map": None, "element": "C"}],
         "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
         "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
         "before_order": 0, "after_order": 1, "change_type": "formed"},
    ]
    labels = [
        {"category_number": 5, "code": "C-N", "supporting_change_ids": ["b1"]},
        {"category_number": 2, "code": "C-C", "supporting_change_ids": ["b2", "b3"]},
    ]
    if brominated:
        changes.append({"change_id": "b4", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "Br"}],
                        "before_order": 0, "after_order": 1, "change_type": "formed"})
        labels.append({"category_number": 9, "code": "C-Br", "supporting_change_ids": ["b4"]})
    why = "Net bond types inferred from the aza-Michael/cyclization product framework; atom mapping is not claimed."
    c["bond_changes"] = [derived(x, ev, why) for x in changes]
    c["labels"] = derived(labels, ev, why)
    c["coverage"] = inferred("partial", ev, "Only confidently visible net formed-bond types are represented; bond-order and hydrogen changes remain unresolved.")
    c["limitations"] = ["No atom mapping is asserted.", "Hydrogen and exact bond-order changes are unresolved.",
                        "Compact source drawings require expert identity review."]


def outcome(kind, oid, targets, raw, normalized, ev, basis=None, yield_type=None, scope=None):
    o = st.new_outcome(oid, kind, targets, ev)
    o["value"] = F(raw, normalized, ev)
    if basis is not None:
        o["measurement_basis"] = F(basis, basis, ev)
    if yield_type is not None:
        o["yield_type"] = F(yield_type, yield_type, ev)
    if scope is not None:
        o["scope"] = F(scope, scope, ev)
    return o


RID_GENERAL = f"{DOC}:general_current_work"
general = base_reaction("general_current_work", (P[1], P[2]), "generic_scheme", ["other"],
                        "Page 1 comprehensive graphic and page 2 Scheme 1b", "Construction of 3-azabicyclo[3.1.0]hexane backbone",
                        "tandem aza-Michael addition, oxidative carbanion radical cyclization and cyclopropanation")
general["participants"] = [participant("gen1", "m_gen_1", "reactant", (P[1], P[2])),
                           participant("gen2", "m_gen_2", "reactant", (P[1], P[2])),
                           participant("gen3", "m_gen_3", "product", (P[1], P[2]), product_role="main product"),
                           participant("g_dbu", "m_DBU", "reagent", (P[1], P[2])),
                           participant("g_cui", "m_CuI", "catalyst", (P[1], P[2])),
                           participant("g_thf", "m_THF", "solvent", (P[1],))]
general["solvent_conditions"]["composition"] = F("THF", ["THF"], (P[1],))
general["steps"] = [st.new_step("one_pot", (P[1], P[2]))]
general["steps"][0]["participant_ids"] = [p["participant_id"] for p in general["participants"]]
general["steps"][0]["operation"] = F("3 Steps in One-Pot", "aza-Michael addition; radical cyclization; oxidative cyclopropanation", (P[1],))
general["steps"][0]["temperature"] = F("50 °C", qty(50, "°C"), (P[1],))
set_classification(general, (P[1], P[2]))
general["review"]["unresolved"] = ["Generic R groups and exact atom mapping remain unresolved.", "No generic reaction time is reported on pages 1-2."]
bundle["reactions"].append(general)

lit = base_reaction("scheme1a_baskaran", (P[2],), "literature_example", ["other"],
                    "Page 2 Scheme 1a", "Baskaran's prior one-pot cyclopropane-fused cyclic amidine synthesis",
                    "prior-literature aza-Michael/oxidative cyclization")
lit["participants"] = [participant("lit_i", "m_lit_I", "reactant", (P[2],)),
                       participant("lit_amine", "m_lit_allylamine", "reactant", (P[2],)),
                       participant("lit_ii", "m_lit_II", "product", (P[2],), product_role="drawn intermediate"),
                       participant("lit_p", "m_lit_product", "product", (P[2],), product_role="literature product"),
                       participant("lit_dbu", "m_DBU", "reagent", (P[2],)),
                       participant("lit_cubr", "m_CuBr2", "reagent", (P[2],))]
lit["steps"] = [st.new_step("prior_step1", (P[2],)), st.new_step("prior_step2", (P[2],))]
lit["steps"][0]["participant_ids"] = ["lit_i", "lit_amine", "lit_dbu", "lit_ii"]
lit["steps"][0]["operation"] = F("RNH-allyl, DBU", "aza-Michael addition", (P[2],))
lit["steps"][1]["participant_ids"] = ["lit_ii", "lit_cubr", "lit_dbu", "lit_p"]
lit["steps"][1]["operation"] = F("CuBr2, DBU", "oxidative cyclization", (P[2],))
set_classification(lit, (P[2],), unresolved=True)
lit["review"]["unresolved"] = ["Cited prior work is shown only as a schematic example; quantities, time, temperature and yield are absent."]
bundle["reactions"].append(lit)


# Table 1: every row is a distinct optimization experiment.
t1 = [
    (1, "CuBr2", 2.0, "THF", "DBU", 2.2, 60, 17, 45),
    (2, "CuBr2", 1.5, "THF", "DBU", 2.2, 60, 50, 14),
    (3, "CuI", 1.5, "THF", "DBU", 2.2, 60, 72, None),
    (4, "FeCl3", 1.5, "THF", "DBU", 2.2, 60, None, None),
    (5, "CuOAc2", 1.5, "THF", "DBU", 2.2, 60, "trace", None),
    (6, "CuI", 2.5, "THF", "DBU", 2.2, 60, 68, None),
    (7, "CuI", 2.0, "THF", "DBU", 2.2, 60, 80, None),
    (8, "CuI", 2.0, "dioxane", "DBU", 2.2, 60, 43, None),
    (9, "CuI", 2.0, "DCM", "DBU", 2.2, "reflux", 40, None),
    (10, "CuI", 2.0, "DCE", "DBU", 2.2, 60, 57, None),
    (11, "CuI", 2.0, "diethyl_ether", "DBU", 2.2, "reflux", 56, None),
    (12, "CuI", 2.0, "toluene", "DBU", 2.2, 60, 45, None),
    (13, "CuI", 2.0, "MeCN", "DBU", 2.2, 60, 52, None),
    (14, "CuI", 2.0, "THF", "DMAP", 2.2, 60, 56, None),
    (15, "CuI", 2.0, "THF", "DABCO", 2.2, 60, None, None),
    (16, "CuI", 2.0, "THF", "dibutylamine", 2.2, 60, 55, None),
    (17, "CuI", 2.0, "THF", "NaOH", 2.2, 60, None, None),
    (18, "CuI", 2.0, "THF", "Cs2CO3", 2.2, 60, None, None),
    (19, "CuI", 2.0, "THF", "DBU", 2.2, 30, 78, None),
]
no_product_rows = {4, 15, 17, 18}
for row, cat, cat_eq, solv, base, base_eq, temp, y3, y4 in t1:
    rid = f"table1_entry_{row:02d}"
    r = base_reaction(rid, (P[2],), "specific_experiment", ["optimization"],
                      f"Page 2 Table 1 entry {row}", "Optimization of reaction conditions",
                      "tandem aza-Michael addition and oxidative cyclization")
    r["general_scheme_ids"] = [RID_GENERAL]
    ps = [participant("t1_1a", "m_1a", "reactant", (P[2],), amount=0.2),
          participant("t1_2a", "m_2a", "reactant", (P[2],), amount=0.2),
          participant("t1_cat", f"m_{cat}", "catalyst", (P[2],), equiv=cat_eq, ref="m_1a"),
          participant("t1_base", f"m_{base}", "reagent", (P[2],), equiv=base_eq, ref="m_1a"),
          participant("t1_solv", f"m_{solv}", "solvent", (P[2],), volume=1.5)]
    if row not in no_product_rows:
        ps.append(participant("t1_3aa", "m_3aa", "product", (P[2],), product_role="main product"))
    if y4 is not None:
        ps.append(participant("t1_4aa", "m_4aa", "product", (P[2],), product_role="brominated byproduct"))
    r["participants"] = ps
    r["reaction_scale"] = F("1a (0.2 mmol)", {"reference_molecule_id": "m_1a", "amount": qty(0.2, "mmol")}, (P[2],))
    r["solvent_conditions"]["composition"] = F(display.get(solv, solv.replace("_", " ")), [f"m_{solv}"], (P[2],))
    r["solvent_conditions"]["total_volume"] = F("1.5 mL", qty(1.5, "mL"), (P[2],))
    s = st.new_step("reaction", (P[2],))
    s["participant_ids"] = [p["participant_id"] for p in ps if p["stage"] == "reaction"]
    s["operation"] = F(f"Table 1 entry {row}", "conduct reaction under listed conditions", (P[2],))
    s["temperature"] = F("reflux" if temp == "reflux" else f"{temp} °C",
                             "reflux" if temp == "reflux" else qty(temp, "°C"), (P[2],))
    s["atmosphere"] = F("under air atmosphere", "air", (P[2],))
    r["steps"] = [s]
    if y3 == "trace":
        r["outcomes"].append(outcome("trace", "3aa_trace", ["m_3aa"], "trace", "trace", (P[2],), basis="isolated yield column footnote"))
    elif isinstance(y3, int):
        r["outcomes"].append(outcome("yield", "3aa_yield", ["m_3aa"], f"{y3}%", qty(y3, "%"), (P[2],),
                                             basis="isolated", yield_type="isolated yield"))
    else:
        r["outcomes"].append(outcome("no_reaction", "3aa_dash", [], "–", "table dash (–)", (P[2],),
                                             basis="isolated yield column footnote"))
    if y4 is not None:
        r["outcomes"].append(outcome("yield", "4aa_yield", ["m_4aa"], f"{y4}%", qty(y4, "%"), (P[2],),
                                             basis="isolated", yield_type="isolated yield"))
    set_classification(r, (P[2],), observed=row not in no_product_rows, brominated=y4 is not None)
    r["review"]["unresolved"] = ["Reaction time is not reported in Table 1.", "Product stereochemistry is not named and is not encoded."]
    bundle["reactions"].append(r)


# Table 2: every row is a distinct optimization experiment.
t2 = [
    (1, 1.0, 2.2, 2.0, 30, 34), (2, 1.1, 2.2, 2.0, 30, 40),
    (3, 1.3, 2.2, 2.0, 30, 40), (4, 1.5, 2.2, 2.0, 30, 40),
    (5, 2.0, 2.2, 2.0, 30, 40), (6, 1.1, 2.5, 2.0, 30, 42),
    (7, 1.1, 2.8, 2.0, 30, 48), (8, 1.1, 3.0, 2.0, 30, 48),
    (9, 1.1, 3.5, 2.0, 30, 47), (10, 1.1, 2.8, 2.2, 30, 49),
    (11, 1.1, 2.8, 2.4, 30, 51), (12, 1.1, 2.8, 2.7, 30, 45),
    (13, 1.1, 2.8, 3.0, 30, 42), (14, 1.1, 2.8, 2.4, 40, 50),
    (15, 1.1, 2.8, 2.4, 50, 54), (16, 1.1, 2.8, 2.4, 70, 48),
    (17, 1.1, 2.8, 2.4, 90, 43),
]
for row, ratio2, dbu_eq, cui_eq, temp, yld in t2:
    r = base_reaction(f"table2_entry_{row:02d}", (P[3],), "specific_experiment", ["optimization"],
                      f"Page 3 Table 2 entry {row}", "Further screening of reaction conditions",
                      "tandem aza-Michael addition and oxidative cyclization")
    r["general_scheme_ids"] = [RID_GENERAL]
    ps = [participant("t2_1k", "m_1k", "reactant", (P[3],), amount=0.2, equiv=1.0, ref="m_1k"),
          participant("t2_2a", "m_2a", "reactant", (P[3],), equiv=ratio2, ref="m_1k"),
          participant("t2_dbu", "m_DBU", "reagent", (P[3],), equiv=dbu_eq, ref="m_1k"),
          participant("t2_cui", "m_CuI", "catalyst", (P[3],), equiv=cui_eq, ref="m_1k"),
          participant("t2_thf", "m_THF", "solvent", (P[3],), volume=1.5),
          participant("t2_3ka", "m_3ka", "product", (P[3],), product_role="main product")]
    r["participants"] = ps
    r["reaction_scale"] = F("1k (0.2 mmol)", {"reference_molecule_id": "m_1k", "amount": qty(0.2, "mmol")}, (P[3],))
    r["solvent_conditions"]["composition"] = F("THF", ["m_THF"], (P[3],))
    r["solvent_conditions"]["total_volume"] = F("1.5 mL", qty(1.5, "mL"), (P[3],))
    s = st.new_step("reaction", (P[3],))
    s["participant_ids"] = [p["participant_id"] for p in ps]
    s["operation"] = F(f"Table 2 entry {row}", "conduct reaction under listed conditions", (P[3],))
    s["temperature"] = F(f"{temp} °C", qty(temp, "°C"), (P[3],))
    s["atmosphere"] = F("under air atmosphere", "air", (P[3],))
    r["steps"] = [s]
    r["outcomes"] = [outcome("yield", "3ka_yield", ["m_3ka"], f"{yld}%", qty(yld, "%"), (P[3],),
                                     basis="isolated", yield_type="isolated yield")]
    set_classification(r, (P[3],))
    r["review"]["unresolved"] = ["Reaction time is not reported in Table 2.", "Full stereochemical structure of 1k/3ka is not encoded."]
    bundle["reactions"].append(r)


# Scheme 2 substrate scope.  The scheme arrow, footnote and page-4 procedure are
# intentionally retained as conflicts rather than silently selecting one set.
scope_map = [("1a", "2a", "3aa"), ("1b", "2a", "3ba"), ("1c", "2a", "3ca"),
             ("1d", "2a", "3da"), ("1e", "2a", "3ea"), ("1f", "2a", "3fa"),
             ("1g", "2a", "3ga"), ("1h", "2a", "3ha"), ("1i", "2a", "3ia"),
             ("1j", "2a", "3ja"), ("1k", "2a", "3ka"), ("1l", "2a", "3la"),
             ("1m", "2a", "3ma"), ("1n", "2a", "3na"), ("1o", "2a", "3oa"),
             ("1p", "2a", "3pa"), ("1q", "2a", "3qa"), ("1r", "2a", "3ra"),
             ("1b", "2b", "3bb")]
scope_conflict_note = (
    "Source conflict retained: Scheme 2 arrow gives DBU 2.8 equiv, CuI 2.4 equiv and 50 °C; "
    "Scheme 2 footnote gives DBU 0.44 mmol, CuI 0.4 mmol and 30 °C on 0.2 mmol substrates; "
    "page 4 Experimental gives allene 0.22 mmol, CuI 0.48 mmol, DBU 0.56 mmol and 50 °C."
)
scope_ids = []
for a, allene, prod in scope_map:
    rid = f"scheme2_{prod}"
    scope_ids.append(f"{DOC}:{rid}")
    r = base_reaction(rid, (P[3], P[4]), "specific_experiment", ["substrate_scope"],
                      f"Page 3 Scheme 2 product {prod}; page 4 Experimental general procedure",
                      "Investigation of substrate scope of allylamine and allene",
                      "tandem aza-Michael addition and oxidative cyclization")
    r["general_scheme_ids"] = [RID_GENERAL]
    ps = [participant("sc_amine", f"m_{a}", "reactant", (P[3], P[4]), amount=0.2),
          participant("sc_allene", f"m_{allene}", "reactant", (P[3], P[4])),
          participant("sc_dbu", "m_DBU", "reagent", (P[3], P[4]), ref=f"m_{a}"),
          participant("sc_cui", "m_CuI", "catalyst", (P[3], P[4]), ref=f"m_{a}"),
          participant("sc_thf", "m_THF", "solvent", (P[3], P[4]), volume=1.5),
          participant("sc_product", f"m_{prod}", "product", (P[3], P[4]), product_role="main product"),
          participant("sc_silica", "m_silica", "purification_agent", (P[4],), stage="purification"),
          participant("sc_pe", "m_petroleum_ether", "purification_agent", (P[4],), stage="purification"),
          participant("sc_ea", "m_ethyl_acetate", "purification_agent", (P[4],), stage="purification")]
    ps[1]["charges"]["amount"] = conflict([
        ("2 (0.2 mmol), Scheme 2 footnote", qty(0.2, "mmol"), (P[3],)),
        ("allene 2 (0.22 mmol), Experimental", qty(0.22, "mmol"), (P[4],)),
    ], (P[3], P[4]))
    ps[2]["charges"]["equivalents"] = mixed_conflict([
        ("DBU (2.8 equiv), Scheme 2 arrow", qty(2.8, "equiv"), (P[3],), "original_extraction", None),
        ("0.44 mmol / 0.2 mmol = 2.2 equiv", qty(2.2, "equiv"), (P[3],), "inference", "Calculated from the Scheme 2 footnote charges."),
        ("0.56 mmol / 0.2 mmol = 2.8 equiv", qty(2.8, "equiv"), (P[4],), "inference", "Calculated from the Experimental charges."),
    ], (P[3], P[4]))
    ps[2]["charges"]["amount"] = conflict([
        ("DBU (0.44 mmol), Scheme 2 footnote", qty(0.44, "mmol"), (P[3],)),
        ("DBU (0.56 mmol), Experimental", qty(0.56, "mmol"), (P[4],)),
    ], (P[3], P[4]))
    ps[3]["charges"]["equivalents"] = mixed_conflict([
        ("CuI (2.4 equiv), Scheme 2 arrow", qty(2.4, "equiv"), (P[3],), "original_extraction", None),
        ("0.4 mmol / 0.2 mmol = 2.0 equiv", qty(2.0, "equiv"), (P[3],), "inference", "Calculated from the Scheme 2 footnote charges."),
        ("0.48 mmol / 0.2 mmol = 2.4 equiv", qty(2.4, "equiv"), (P[4],), "inference", "Calculated from the Experimental charges."),
    ], (P[3], P[4]))
    ps[3]["charges"]["amount"] = conflict([
        ("CuI (0.4 mmol), Scheme 2 footnote", qty(0.4, "mmol"), (P[3],)),
        ("CuI (0.48 mmol), Experimental", qty(0.48, "mmol"), (P[4],)),
    ], (P[3], P[4]))
    r["participants"] = ps
    r["reaction_scale"] = F("allylamine 1 (0.2 mmol)", {"reference_molecule_id": f"m_{a}", "amount": qty(0.2, "mmol")}, (P[3], P[4]))
    r["solvent_conditions"]["composition"] = F("THF", ["m_THF"], (P[3], P[4]))
    r["solvent_conditions"]["total_volume"] = F("1.5 mL", qty(1.5, "mL"), (P[3], P[4]))
    s1 = st.new_step("pre_copper_reaction", (P[3], P[4]))
    s1["participant_ids"] = ["sc_amine", "sc_allene", "sc_thf"]
    s1["operation"] = F("Allylamine 1 and allene 2 were added to THF; stirred until nearly completed monitored by TLC.",
                         "combine allylamine and allene in THF; stir until nearly complete by TLC", (P[4],))
    s1["temperature"] = conflict([
        ("30 °C, Scheme 2 footnote", qty(30, "°C"), (P[3],)),
        ("50 °C, Scheme 2 arrow and Experimental", qty(50, "°C"), (P[3], P[4])),
    ], (P[3], P[4]))
    s1["time"] = F("until the reaction was nearly completed monitored by TLC", "until nearly complete by TLC", (P[4],))
    s1["atmosphere"] = F("under air atmosphere", "air", (P[3],))
    s2 = st.new_step("copper_base_stage", (P[4],))
    s2["participant_ids"] = ["sc_amine", "sc_allene", "sc_dbu", "sc_cui", "sc_thf", "sc_product"]
    s2["operation"] = F("Afterwards, CuI and DBU were added; after completion indicated by TLC, concentrate under reduced pressure.",
                         "add CuI and DBU; react to TLC completion; concentrate", (P[4],))
    s2["temperature"] = conflict([
        ("30 °C, Scheme 2 footnote", qty(30, "°C"), (P[3],)),
        ("50 °C, Scheme 2 arrow and Experimental", qty(50, "°C"), (P[3], P[4])),
    ], (P[3], P[4]))
    s2["time"] = F("until completion indicated by TLC", "until TLC completion", (P[4],))
    s2["atmosphere"] = F("under air atmosphere", "air", (P[3],))
    s3 = st.new_step("flash_chromatography", (P[4],), "purification")
    s3["participant_ids"] = ["sc_product", "sc_silica", "sc_pe", "sc_ea"]
    s3["operation"] = F("flash chromatography on silica gel (petroleum ether/ethyl acetate = 2/1)",
                         {"method": "flash chromatography", "stationary_phase": "silica gel", "eluent_ratio": "petroleum ether/ethyl acetate 2:1"}, (P[4],))
    r["steps"] = [s1, s2, s3]
    r["procedure"]["original_text"] = F(
        "Allylamine 1 (0.2 mmol) and allene 2 (0.22 mmol) were added to THF (1.5 mL) at 50 °C ... CuI (0.48 mmol) and DBU (0.56 mmol) were added ... flash chromatography ... petroleum ether/ethyl acetate = 2/1.",
        "page 4 Experimental general procedure", (P[4],))
    r["procedure"]["general_procedure_reference"] = F("Experimental, page 669", "page 4 Experimental", (P[4],))
    r["procedure"]["addition_order"] = F("allylamine + allene in THF first; CuI + DBU afterwards", ["allylamine", "allene", "THF", "CuI", "DBU"], (P[4],))
    r["outcomes"] = [outcome("yield", f"{prod}_yield", [f"m_{prod}"], f"{product_yields[prod]}% yield",
                                     qty(product_yields[prod], "%"), (P[3],), basis=None, yield_type=None)]
    set_classification(r, (P[3], P[4]))
    r["notes"] = [F(scope_conflict_note, scope_conflict_note, (P[3], P[4]))]
    r["review"]["unresolved"] = [scope_conflict_note, "Scheme 2 does not state the yield measurement basis.",
                                  "Exact product stereochemistry and full isomeric SMILES remain unresolved."]
    bundle["reactions"].append(r)


# Scheme 3a gram-scale demonstration.
gram = base_reaction("scheme3a_gram_scale_3aa", (P[3],), "specific_experiment", ["scale_up"],
                     "Page 3 Scheme 3a", "Gram-scale experiment", "scale-up of tandem aza-Michael/oxidative cyclization")
gram["general_scheme_ids"] = [RID_GENERAL]
gram["participants"] = [
    participant("gr_1a", "m_1a", "reactant", (P[3],), amount=7.27, mass=(706, "mg")),
    participant("gr_2a", "m_2a", "reactant", (P[3],), amount=8.0, mass=(1.25, "g")),
    participant("gr_dbu", "m_DBU", "reagent", (P[3],)), participant("gr_cui", "m_CuI", "catalyst", (P[3],)),
    participant("gr_thf", "m_THF", "solvent", (P[3],)),
    participant("gr_3aa", "m_3aa", "product", (P[3],), mass=(1.5, "g"), product_role="main product")]
gram["reaction_scale"] = F("1a, 7.27 mmol", {"reference_molecule_id": "m_1a", "amount": qty(7.27, "mmol")}, (P[3],))
gram["solvent_conditions"]["composition"] = F("THF", ["m_THF"], (P[3],))
sg = st.new_step("reaction", (P[3],)); sg["participant_ids"] = [p["participant_id"] for p in gram["participants"]]
sg["operation"] = F("DBU, THF, CuI", "gram-scale reaction under scheme conditions", (P[3],))
sg["temperature"] = F("50 °C", qty(50, "°C"), (P[3],)); gram["steps"] = [sg]
gram["outcomes"] = [outcome("yield", "gram_yield", ["m_3aa"], "yield 82%", qty(82, "%"), (P[3],)),
                    outcome("product_amount", "gram_product_mass", ["m_3aa"], "1.5 g", qty(1.5, "g"), (P[3],))]
set_classification(gram, (P[3],))
gram["review"]["unresolved"] = ["DBU, CuI and THF amounts, reaction time, atmosphere, workup and yield basis are not reported in Scheme 3a."]
bundle["reactions"].append(gram)


# Scheme 3b derivatization.
deriv = base_reaction("scheme3b_derivatization_5", (P[3],), "specific_experiment", ["derivatization"],
                      "Page 3 Scheme 3b", "Derivatization of product", "alkaline aqueous transformation to carboxylic acid derivative 5")
deriv["participants"] = [participant("d_3aa", "m_3aa", "reactant", (P[3],)),
                         participant("d_koh", "m_KOH", "reagent", (P[3],)),
                         participant("d_water", "m_water", "solvent", (P[3],)),
                         participant("d_hcl", "m_HCl", "workup_agent", (P[3],), stage="workup"),
                         participant("d_5", "m_5", "product", (P[3],), product_role="derivatization product")]
sd1 = st.new_step("alkaline_aqueous_treatment", (P[3],)); sd1["participant_ids"] = ["d_3aa", "d_koh", "d_water", "d_5"]
sd1["operation"] = F("(1) KOH/H2O", "treat with aqueous KOH", (P[3],))
sd2 = st.new_step("acidification", (P[3],), "workup"); sd2["participant_ids"] = ["d_hcl", "d_5"]
sd2["operation"] = F("(2) HCl (pH=1)", "acidify with HCl to pH 1", (P[3],))
deriv["steps"] = [sd1, sd2]
deriv["outcomes"] = [outcome("yield", "product_5_yield", ["m_5"], "92%", qty(92, "%"), (P[3],))]
set_classification(deriv, (P[3],), unresolved=True)
deriv["review"]["unresolved"] = ["Amounts, temperature, time and yield measurement basis are not reported.",
                                  "Exact net bond changes cannot be mapped safely from the compact drawings."]
bundle["reactions"].append(deriv)


# Proposed mechanism, not an observed mechanistic experiment.
mech = base_reaction("scheme4_proposed_mechanism", (P[1], P[4]), "proposed_mechanism", ["other"],
                     "Page 4 Scheme 4 (with overview on page 1)", "Proposed mechanism",
                     "proposed aza-Michael addition, SET oxidation, 5-exo-trig radical cyclization and oxidative cyclopropanation")
mech["participants"] = [participant("me_1a", "m_1a", "reactant", (P[4],)),
                        participant("me_2a", "m_2a", "reactant", (P[4],)),
                        participant("me_dbu", "m_DBU", "reagent", (P[4],)),
                        participant("me_cui", "m_CuI", "reagent", (P[4],)),
                        participant("me_3aa", "m_3aa", "product", (P[4],), product_role="proposed final product")]
for i, (raw, norm) in enumerate([
    ("DBU", "aza-Michael addition/deprotonation to A"),
    ("SET oxidation, CuI", "SET oxidation and radical cyclization A to B/C"),
    ("DBU", "deprotonation C to D"),
    ("SET oxidation, CuI", "oxidative cyclopropanation D to 3aa"),
], 1):
    s = st.new_step(f"proposed_step_{i}", (P[4],))
    s["participant_ids"] = [p["participant_id"] for p in mech["participants"]]
    s["operation"] = F(raw, norm, (P[4],)); mech["steps"].append(s)
set_classification(mech, (P[1], P[4]))
mech["review"]["unresolved"] = ["Scheme 4 is an author-proposed mechanism, not direct experimental proof.",
                                 "Intermediate atom mapping, hydrogen transfers and copper speciation are unresolved."]
bundle["reactions"].append(mech)


# Non-synthetic observations and page-level context.
bundle["source_context"] = [
    {"context_id": f"{DOC}:figure1_bioactive_context", "document_id": DOC,
     "evidence_ids": [P[2]], "kind": "background_figure",
     "note": "Figure 1 contains bioactive example structures, not synthesis experiments in this paper."},
    {"context_id": f"{DOC}:xray_3ra", "document_id": DOC, "evidence_ids": [P[3]],
     "kind": "structure_assignment", "note": "The structure of 3ra was established by single-crystal X-ray diffraction."},
    {"context_id": f"{DOC}:xray_5", "document_id": DOC, "evidence_ids": [P[3]],
     "kind": "structure_assignment", "note": "The structure and relative configuration of 5 were established by single-crystal X-ray diffraction."},
    {"context_id": f"{DOC}:si_statement", "document_id": DOC, "evidence_ids": [P[4]],
     "kind": "availability_statement", "note": "Article states SI is available online, but no SI file was supplied or explicitly associated in this session."},
    {"context_id": f"{DOC}:page5_references", "document_id": DOC, "evidence_ids": [P[5]],
     "kind": "references", "note": "Page 5 contains continuation of references and publication dates; no current-paper synthesis experiment is shown."},
]

save(RUN / "submitted.json", bundle)

# Machine visual inventory is kept separate from the formal records.
visual_inventory = {
    "document_id": DOC,
    "pages_visually_inspected": [1, 2, 3, 4, 5],
    "chemical_accuracy_verified": False,
    "items": [
        {"page": 1, "items": ["title and abstract", "generic current-work reaction", "mechanistic overview A-C"]},
        {"page": 2, "items": ["Figure 1 background molecules", "Scheme 1a prior literature", "Scheme 1b current generic scheme", "Table 1 entries 1-19 and footnotes"]},
        {"page": 3, "items": ["Table 2 entries 1-17 and footnotes", "Scheme 2 products 3aa-3ra and 3bb", "Scheme 3a scale-up", "Scheme 3b derivatization", "X-ray renderings"]},
        {"page": 4, "items": ["Scheme 4 proposed mechanism", "Experimental general procedure", "SI availability statement", "references"]},
        {"page": 5, "items": ["references continuation", "publication dates"]},
    ],
    "source_conflicts": [scope_conflict_note],
}
save(RUN / "raw_model/page1_5_visual_inventory.json", visual_inventory)
save(RUN / "raw_model/stereo_diagnostics.json", STEREO_DIAGNOSTICS)
mentions = st.source_mentions(json.loads((RUN / "dataset.json").read_text()),
                              ["1a", "2a", "3aa", "1k", "3ka", "3ra", "5"])
save(RUN / "raw_model/source_mentions.json", mentions)

all_ids = [r["reaction_id"] for r in bundle["reactions"]]
t1_ids = [f"{DOC}:table1_entry_{i:02d}" for i in range(1, 20)]
t2_ids = [f"{DOC}:table2_entry_{i:02d}" for i in range(1, 18)]
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title/abstract", "finding": "article metadata and reported overall 42%-85% range; not a distinct experiment", "evidence_ids": [P[1]]},
            {"region": "Comprehensive Graphic Content", "finding": "generic current-work reaction plus mechanistic overview; mapped to generic scheme and proposed mechanism", "evidence_ids": [P[1]]}],
         "record_ids": [RID_GENERAL, f"{DOC}:scheme4_proposed_mechanism"],
         "unresolved": ["Generic R groups and exact atom mapping are unresolved."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Figure 1", "finding": "background bioactive structures; non-synthetic context", "evidence_ids": [P[2]]},
            {"region": "Scheme 1a", "finding": "prior Baskaran literature example, not a current-paper experiment", "evidence_ids": [P[2]]},
            {"region": "Scheme 1b", "finding": "current generic transformation, duplicate context for page-1 generic scheme", "evidence_ids": [P[2]]},
            {"region": "Table 1 rows 1-19 and footnotes", "finding": "19 independent optimization experiments; rows 1-2 include brominated byproduct 4aa", "evidence_ids": [P[2]]}],
         "record_ids": [f"{DOC}:scheme1a_baskaran", RID_GENERAL] + t1_ids,
         "unresolved": ["Table 1 does not report reaction time; primary product stereochemistry is not named."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Table 2 rows 1-17 and footnotes", "finding": "17 independent optimization experiments", "evidence_ids": [P[3]]},
            {"region": "Scheme 2 arrow, products and footnote", "finding": "19 substrate-scope experiments; arrow and footnote conditions conflict", "evidence_ids": [P[3]]},
            {"region": "Scheme 3a", "finding": "one gram-scale experiment", "evidence_ids": [P[3]]},
            {"region": "Scheme 3b", "finding": "one derivatization experiment to product 5", "evidence_ids": [P[3]]},
            {"region": "3ra and 5 molecular renderings/prose", "finding": "X-ray structure assignments; not separate syntheses", "evidence_ids": [P[3]]}],
         "record_ids": t2_ids + scope_ids + [f"{DOC}:scheme3a_gram_scale_3aa", f"{DOC}:scheme3b_derivatization_5"],
         "unresolved": [scope_conflict_note, "Scheme 2 and Scheme 3 do not state yield measurement bases."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "Scheme 4", "finding": "author-proposed mechanism, not a mechanistic experiment", "evidence_ids": [P[4]]},
            {"region": "Experimental", "finding": "general two-stage procedure and purification for products 3; linked to Scheme 2 records", "evidence_ids": [P[4]]},
            {"region": "Supporting Information statement", "finding": "SI availability stated, but no SI file supplied or associated", "evidence_ids": [P[4]]},
            {"region": "References", "finding": "citations only; no additional current-paper experiments", "evidence_ids": [P[4]]}],
         "record_ids": [f"{DOC}:scheme4_proposed_mechanism"] + scope_ids,
         "unresolved": ["General procedure conflicts with Scheme 2 footnote for allene/CuI/DBU amounts and temperature.", "SI association remains unestablished."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "full page", "finding": "references continuation and manuscript dates; no current-paper reaction experiment or standalone characterization table", "evidence_ids": [P[5]]}],
         "record_ids": [], "unresolved": []},
    ],
}
save(RUN / "coverage/pages_0001_0005.json", coverage)

tasks = json.loads((RUN / "tasks.json").read_text())
for task in tasks:
    if task["document_id"] == DOC and task["evidence_id"] in PT.values():
        task["state"] = "attempted"
        task["reason"] = "Original page image inspected; findings and record IDs are documented in coverage/pages_0001_0005.json."
save(RUN / "tasks.json", tasks)

print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
                  "scope_records": len(scope_ids)}, ensure_ascii=False))
