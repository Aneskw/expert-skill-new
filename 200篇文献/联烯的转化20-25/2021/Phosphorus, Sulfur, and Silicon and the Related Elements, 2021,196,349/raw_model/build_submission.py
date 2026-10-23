#!/usr/bin/env python3
"""Build the evidence-linked v3.1 machine submission for Liu et al. (2021).

This is a reproducible transcription of the five supplied original PDF pages.
It deliberately leaves source-uncertain identities and procedural details unresolved.
"""

from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills" / "chemical-reaction-extraction" / "scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st


DOC = "doc_a2cb06c9c35e36e6"
P2T = f"{DOC}:page2:text"
P2I = f"{DOC}:page2:image"
P3T = f"{DOC}:page3:text"
P3I = f"{DOC}:page3:image"
P4T = f"{DOC}:page4:text"
P4I = f"{DOC}:page4:image"
EV_SCHEME = [P2I, P2T, P3I, P3T]
EV_TABLE = [P3I, P3T]
EV_SCOPE = [P3I, P3T, P4T]


def source(raw, normalized=None, evidence=(), method="original_extraction"):
    if normalized is None and raw is not None:
        normalized = raw
    return st.fact(raw, normalized, evidence, method=method)


def derived(normalized, evidence, derivation, raw=None, method="structure_derivation"):
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation)


def q(number, unit, **extra):
    value = {"number": number, "unit": unit}
    value.update(extra)
    return value


def mol_id(label):
    return f"{DOC}:mol_{label}"


def rxn_id(label):
    return f"{DOC}:{label}"


def set_status(molecule, key, value, evidence, reason):
    molecule["stereochemistry"][key] = derived(
        value, evidence, reason, raw={"source": "structure drawing and candidate graph"}, method="inference"
    )


def make_molecule(label, name, evidence, smiles=None, chirality="unresolved",
                  geometry="unresolved", representation="unresolved", note=None,
                  abbreviations=None):
    mid = mol_id(label)
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = source(label, label, evidence)
    m["names"]["original"] = source(name, name, evidence)
    m["names"]["english"] = source(name, name, evidence)
    if abbreviations:
        m["names"]["abbreviations"] = source(abbreviations, list(abbreviations), evidence)
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = derived(
            smiles,
            evidence,
            "Machine candidate transcribed from the displayed structure or explicit chemical name; RDKit parsing does not verify source identity.",
            raw={"source_structure_label": label},
        )
        diag = st.stereo_inspect(smiles)
        m["stereochemistry"]["atom_centers"] = derived(
            diag["atom_centers"], evidence, "RDKit diagnostic on the machine candidate SMILES; source configuration remains governed by the drawing.",
            raw={"candidate_smiles": smiles},
        )
        m["stereochemistry"]["double_bonds"] = derived(
            diag["double_bonds"], evidence, "RDKit diagnostic on the machine candidate SMILES; terminal double bonds have no E/Z assignment.",
            raw={"candidate_smiles": smiles},
        )
        m["stereochemistry"]["identity_key"] = derived(
            diag["stereo_graph_key"], evidence, "Stereo-preserving graph key calculated from the machine candidate SMILES.",
            raw={"candidate_smiles": smiles},
        )
    set_status(m, "chirality_status", chirality, evidence,
               "Status follows the displayed structure plus graph diagnostics; silence is not treated as racemic.")
    set_status(m, "geometry_status", geometry, evidence,
               "Status follows the displayed double bonds; terminal alkenes are not assigned E/Z.")
    set_status(m, "representation_status", representation, evidence,
               "The representation status records whether source stereochemistry is encoded, not whether SMILES parses.")
    if note:
        m["identity_note"] = source(note, note, evidence)
    return m


def template_fact(raw, normalized, evidence, derivation=None):
    if derivation:
        return derived(normalized, evidence, derivation, raw=raw)
    return source(raw, normalized, evidence)


def make_template(label, name, evidence, fixed, variables, reason,
                  chirality="unresolved", geometry="unresolved", representation="partial"):
    m = make_molecule(label, name, evidence, chirality=chirality,
                      geometry=geometry, representation=representation, note=reason)
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": template_fact(fixed, fixed, evidence),
        "variables": [
            {
                "label": v,
                "attachment_points": template_fact(f"{v} attached at the labeled allene/product carbon", [v], evidence),
                "definition": template_fact(desc, desc, evidence),
                "allowed_values": st.fact(evidence=evidence),
                "constraints": st.fact(evidence=evidence),
            }
            for v, desc in variables
        ],
        "instance_mapping": [],
        "reason": template_fact(reason, reason, evidence),
    }
    return m


def participant(pid, mid, role, evidence, stage="reaction", charges=None,
                reference=None, product_role=None):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    if product_role:
        p["product_role"] = source(product_role, product_role, evidence)
    for key, (raw, normalized, ev) in (charges or {}).items():
        p["charges"][key] = source(raw, normalized, ev)
    if reference:
        p["equivalence_reference"] = source(reference, reference, evidence)
    return p


def conflict_mol_percent(evidence):
    p3 = source("Pd(acac)2 (2 mol%)", q(2, "mol%"), [P3I, P3T])
    p4 = st.fact("Pd(acac)2 (6 mg, 2 mmol%)", None, [P4I, P4T],
                 status="unreadable", verification="unreviewed",
                 derivation="The printed unit mmol% is not a supported loading unit and is preserved verbatim.")
    return st.fact(["2 mol%", "2 mmol%"], None, evidence, status="conflict",
                   alternatives=[p3, p4],
                   derivation="Scheme 2 reports 2 mol%; the printed general procedure says 2 mmol%. No preferred reading is selected.")


def step(step_id, evidence, participant_ids, operation, temperature=None,
         time=None, atmosphere=None, stage="reaction"):
    s = st.new_step(step_id, evidence, stage=stage)
    s["participant_ids"] = list(participant_ids)
    s["operation"] = source(operation, operation, evidence)
    if temperature:
        s["temperature"] = source(temperature[0], temperature[1], evidence)
    if time:
        s["time"] = source(time[0], time[1], evidence)
    if atmosphere:
        s["atmosphere"] = source(atmosphere, atmosphere, evidence)
    return s


def add_hydrosilylation_classification(r, evidence):
    changes = [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "Si"}],
         "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
         "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
        {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}],
         "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
        {"change_id": "b4", "atoms": [{"atom_map": None, "element": "Si"}, {"atom_map": None, "element": "H"}],
         "before_order": 1, "after_order": 0, "change_type": "hydrogen_loss", "hydrogen_count": 1},
    ]
    labels = [
        {"category_number": 12, "code": "C-Si", "supporting_change_ids": ["b1"]},
        {"category_number": 1, "code": "C-H", "supporting_change_ids": ["b3"]},
    ]
    reason = "Net allene hydrosilylation inferred from the displayed substrate/product connectivity; exact atom mapping and coproduct balance are not asserted."
    r["classification"]["bond_changes"] = [derived(c, evidence, reason) for c in changes]
    r["classification"]["labels"] = derived(labels, evidence, reason)
    r["classification"]["coverage"] = derived("partial", evidence, reason, method="inference")
    r["classification"]["limitations"] = [
        "Exact atom mapping was not asserted from the drawing.",
        "Hydrogen transfer is a net change; the elementary proton/hydride sequence is not treated as experimentally proven.",
        "Unreported side products and isomer structures are outside this partial classification.",
    ]


def common_reaction(rid, evidence, record_type, categories, location, reaction_name):
    r = st.new_reaction(rxn_id(rid), DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = source(categories, categories, evidence)
    r["location"] = source(location, location, evidence)
    r["reaction_name"] = source(reaction_name, reaction_name, evidence)
    r["transformation_type"] = derived(
        "regioselective allene hydrosilylation to a branched allylsilane",
        evidence,
        "Transformation type inferred from the source reaction name and displayed connectivity.",
        raw=reaction_name,
        method="inference",
    )
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    return r


def set_solvent(r, evidence, total_volume=None):
    r["solvent_conditions"]["composition"] = source("THF", [mol_id("thf")], evidence)
    if total_volume:
        r["solvent_conditions"]["total_volume"] = source(total_volume[0], total_volume[1], evidence)


def outcome(oid, kind, targets, evidence, raw_value, normalized_value,
            product_role=None, yield_type=None, basis=None, scope=None):
    o = st.new_outcome(oid, kind, targets, evidence)
    o["value"] = source(raw_value, normalized_value, evidence)
    if product_role:
        o["product_role"] = source(product_role, product_role, evidence)
    if yield_type:
        o["yield_type"] = source(yield_type, yield_type, evidence)
    if basis:
        o["measurement_basis"] = source(basis, basis, evidence)
    if scope:
        o["scope"] = source(scope, scope, evidence)
    return o


dataset = json.loads((ROOT / "runs/current/dataset.json").read_text())
bundle = st.empty(dataset["documents"])

# Shared named materials. Catalyst/ligand structures are intentionally not guessed from names.
molecules = [
    make_molecule("phenylsilane", "phenylsilane", EV_SCHEME, "[SiH3]c1ccccc1", "achiral", "not_applicable", "not_applicable", abbreviations=["PhSiH3"]),
    make_molecule("thf", "tetrahydrofuran", EV_TABLE, "C1CCOC1", "achiral", "not_applicable", "not_applicable", abbreviations=["THF"]),
    make_molecule("pd_acac2", "palladium(II) acetylacetonate", EV_SCHEME, note="Source reagent Pd(acac)2; exact coordination representation not encoded.", abbreviations=["Pd(acac)2"]),
    make_molecule("pd_oac2", "palladium(II) acetate", EV_TABLE, note="Source reagent Pd(OAc)2; exact coordination representation not encoded.", abbreviations=["Pd(OAc)2"]),
    make_molecule("pdcl2", "palladium(II) chloride", EV_TABLE, note="Source reagent PdCl2; exact coordination representation not encoded.", abbreviations=["PdCl2"]),
    make_molecule("pd_me_cn2_cl2", "bis(acetonitrile)palladium(II) dichloride", EV_TABLE, note="Source reagent Pd(MeCN)2Cl2; exact coordination representation not encoded.", abbreviations=["Pd(MeCN)2Cl2"]),
    make_molecule("pph3", "triphenylphosphine", EV_TABLE, note="Ligand identity is source-named; structure was not independently encoded.", abbreviations=["PPh3"]),
    make_molecule("dppm", "bis(diphenylphosphino)methane", EV_TABLE, note="Ligand identity is source-named; structure was not independently encoded.", abbreviations=["dppm"]),
    make_molecule("dppbz", "1,2-bis(diphenylphosphino)benzene", EV_TABLE, note="Ligand identity is source-named; structure was not independently encoded.", abbreviations=["dppbz"]),
    make_molecule("dppf", "1,1'-bis(diphenylphosphino)ferrocene", EV_TABLE, note="Ligand identity is source-named; exact organometallic representation was not encoded.", abbreviations=["dppf"]),
    make_molecule("dpephos", "bis[(2-diphenylphosphino)phenyl] ether", EV_TABLE, note="Ligand identity is source-named; structure was not independently encoded.", abbreviations=["DPEphos"]),
    make_molecule("xantphos", "4,5-bis(diphenylphosphino)-9,9-dimethylxanthene", EV_SCHEME, note="Ligand identity is source-named; structure was not independently encoded.", abbreviations=["Xantphos"]),
    make_molecule("tmob", "1,3,5-trimethoxybenzene", EV_TABLE, "COc1cc(OC)cc(OC)c1", "achiral", "not_applicable", "not_applicable"),
    make_molecule("ethyl_acetate", "ethyl acetate", [P4I, P4T], "CCOC(C)=O", "achiral", "not_applicable", "not_applicable"),
    make_molecule("petroleum_ether", "petroleum ether", [P4I, P4T], note="Named chromatographic eluent mixture; composition is not reported."),
    make_molecule("silica_gel", "silica gel", [P4I, P4T], note="Named chromatography stationary phase; no discrete molecular structure is asserted."),
]
molecules[0]["stereochemistry"]["assignment"] = derived(
    {"assessment": "achiral", "basis": "the silicon bears three equivalent hydrogens"},
    EV_SCHEME,
    "Phenylsilane has three equivalent Si-H substituents; the toolkit's generic silicon potential-center diagnostic is not an assigned stereocenter.",
    raw={"candidate_smiles": "[SiH3]c1ccccc1"},
)

generic_allene = make_template(
    "generic_allene", "substituted terminal allene", EV_SCHEME,
    "R1R2C=C=CH2", [("R1", "source-defined organic substituent"), ("R2", "source-defined organic substituent or H")],
    "Generic R-group drawing in the graphical abstract and Scheme 1c; it is not a concrete molecule.",
    chirality="achiral", geometry="not_applicable", representation="partial",
)
generic_product = make_template(
    "generic_branched_allylsilane", "branched allylsilane", EV_SCHEME,
    "R1R2C(SiH2Ph)-CH=CH2", [("R1", "substituent inherited from the allene"), ("R2", "substituent or H inherited from the allene")],
    "Generic R-group product drawing; stereogenicity depends on R1/R2 and no configuration is assigned.",
    chirality="unresolved", geometry="not_applicable", representation="partial",
)
molecules += [generic_allene, generic_product]

# Concrete Scheme 2 candidates. Chain lengths and connectivities are machine transcriptions of the drawings.
scope_rows = [
    ("2a", "1a", "1-octylallene", "CCCCCCCCC=C=C", "branched allylsilane 2a", "CCCCCCCCC(C=C)[SiH2]c1ccccc1", 93),
    ("2b", "1b", "long-chain terminal allene precursor to 2b", "CCCCCCCCCCCCCCC=C=C", "branched allylsilane 2b", "CCCCCCCCCCCCCCC(C=C)[SiH2]c1ccccc1", 91),
    ("2c", "1c", "chloro-functionalized terminal allene precursor to 2c", "ClCCCCC=C=C", "branched allylsilane 2c", "ClCCCCC(C=C)[SiH2]c1ccccc1", 88),
    ("2d", "1d", "hydroxy-functionalized terminal allene precursor to 2d", "OCCCCC=C=C", "branched allylsilane 2d", "OCCCCC(C=C)[SiH2]c1ccccc1", 75),
    ("2e", "1e", "acrylate-functionalized terminal allene precursor to 2e", "C=CC(=O)OCCCC=C=C", "branched allylsilane 2e", "C=CC(=O)OCCCC(C=C)[SiH2]c1ccccc1", 83),
    ("2f", "1f", "p-methoxybenzoate-functionalized terminal allene precursor to 2f", "COc1ccc(C(=O)OCCCC=C=C)cc1", "branched allylsilane 2f", "COc1ccc(C(=O)OCCCC(C=C)[SiH2]c2ccccc2)cc1", 86),
    ("2g", "1g", "acetate-protected terminal allene precursor to 2g", "CC(=O)OCCCC=C=C", "branched allylsilane 2g", "CC(=O)OCCCC(C=C)[SiH2]c1ccccc1", 94),
    ("2h", "1h", "benzyl ether-functionalized terminal allene precursor to 2h", "c1ccccc1COCCCC=C=C", "branched allylsilane 2h", "c1ccccc1COCCCC(C=C)[SiH2]c2ccccc2", 90),
    ("2i", "1i", "tosylate-functionalized terminal allene precursor to 2i", "Cc1ccc(S(=O)(=O)OCCCC=C=C)cc1", "branched allylsilane 2i", "Cc1ccc(S(=O)(=O)OCCCC(C=C)[SiH2]c2ccccc2)cc1", 92),
    ("2j", "1j", "TBDMS ether-functionalized terminal allene precursor to 2j", "CC(C)(C)[Si](C)(C)OCCCC=C=C", "branched allylsilane 2j", "CC(C)(C)[Si](C)(C)OCCCC(C=C)[SiH2]c1ccccc1", 91),
    ("2k", "1k", "p-methoxyphenyl-substituted terminal allene precursor to 2k", "COc1ccc(C=C=C)cc1", "branched allylsilane 2k", "COc1ccc(C(C=C)[SiH2]c2ccccc2)cc1", 91),
    ("2l", "1l", "p-tert-butylphenyl-substituted terminal allene precursor to 2l", "CC(C)(C)c1ccc(C=C=C)cc1", "branched allylsilane 2l", "CC(C)(C)c1ccc(C(C=C)[SiH2]c2ccccc2)cc1", 84),
    ("2m", "1m", "p-tolyl-substituted terminal allene precursor to 2m", "Cc1ccc(C=C=C)cc1", "branched allylsilane 2m", "Cc1ccc(C(C=C)[SiH2]c2ccccc2)cc1", 80),
]

for product_label, substrate_label, substrate_name, substrate_smiles, product_name, product_smiles, _ in scope_rows:
    molecules.append(make_molecule(substrate_label, substrate_name, [P3I], substrate_smiles,
                                   "achiral", "not_applicable", "not_applicable",
                                   note="Candidate connectivity was transcribed from the Scheme 2 structure; no SI identity was available."))
    molecules.append(make_molecule(product_label, product_name, [P3I], product_smiles,
                                   "unspecified", "not_applicable", "not_encoded",
                                   note="The displayed tetrahedral carbon has no wedge/dash or absolute/relative assignment; no racemic claim is inferred."))
    molecules.append(make_template(
        f"other_{product_label}_isomers", f"other hydrosilylation isomers corresponding to {product_label}", [P3I, P3T],
        "unspecified alternative hydrosilylation regio-/stereoisomers", [("parent", f"same substrate family as {product_label}")],
        "The source reports a combined denominator 'all other isomers' without enumerating a unique composition.",
        chirality="mixture", geometry="mixture", representation="unresolved",
    ))

# Mechanistic intermediates are source labels with unresolved exact coordination/structure encodings.
for label, name in [("mech1", "(L)Pd(0) intermediate 1"),
                    ("mech2", "phenylsilyl palladium hydride intermediate 2"),
                    ("mech3", "allene-coordinated palladium hydride intermediate 3"),
                    ("mech4", "allylpalladium intermediate 4")]:
    molecules.append(make_molecule(label, name, [P4I, P4T], note="Proposed catalytic-cycle intermediate; exact coordination SMILES and stereochemistry are unresolved."))

bundle["molecules"] = molecules

# One generic scheme consolidates the graphical abstract, Scheme 1c, and the Scheme 2 header arrow.
generic = common_reaction(
    "generic_branched_hydrosilylation", EV_SCHEME, "generic_scheme", ["other"],
    {"pages": [2, 3], "regions": ["graphical abstract", "Scheme 1c", "Scheme 2 header arrow"]},
    "Pd(acac)2/Xantphos-catalyzed regioselective hydrosilylation of allenes",
)
generic["participants"] = [
    participant("generic_allene", mol_id("generic_allene"), "reactant", EV_SCHEME),
    participant("generic_phsih3", mol_id("phenylsilane"), "reactant", EV_SCHEME),
    participant("generic_pd", mol_id("pd_acac2"), "catalyst", EV_SCHEME,
                charges={"mol_percent": ("2 mol%", q(2, "mol%"), EV_SCHEME)}, reference=mol_id("generic_allene")),
    participant("generic_xantphos", mol_id("xantphos"), "ligand", EV_SCHEME,
                charges={"mol_percent": ("2 mol%", q(2, "mol%"), EV_SCHEME)}, reference=mol_id("generic_allene")),
    participant("generic_thf", mol_id("thf"), "solvent", EV_SCHEME),
    participant("generic_product", mol_id("generic_branched_allylsilane"), "product", EV_SCHEME,
                product_role="branched allylsilane"),
]
generic["steps"] = [step("generic_reaction", EV_SCHEME,
                         [p["participant_id"] for p in generic["participants"]],
                         "allene hydrosilylation", ("30 °C", q(30, "°C")),
                         ("1 h (shown in Scheme 2 header arrow)", q(1, "h")))]
set_solvent(generic, EV_SCHEME)
generic["notes"] = [source(
    "Graphical abstract states high yield up to 94% and r.r. >15:1; these are scope-level summary claims, not assigned as an outcome of the generic scheme.",
    "Summary claims are not converted into a generic experiment yield.", EV_SCHEME
)]
generic["review"]["unresolved"] = ["Generic R groups are not concrete molecular identities.",
                                      "The 1 h shown over the Scheme 2 header arrow is not assigned to every scope member."]
add_hydrosilylation_classification(generic, EV_SCHEME)

reactions = [generic]

# Table 1: one optimization experiment per row.
optimization = [
    (1, "pd_acac2", "pph3", 2, 2, 46, "1:1"),
    (2, "pd_acac2", "dppm", 2, 2, 25, "1:2"),
    (3, "pd_acac2", "dppbz", 2, 2, 15, "1:3"),
    (4, "pd_acac2", "dppf", 2, 2, "trace", "table dash (–)"),
    (5, "pd_acac2", "dpephos", 2, 2, 90, "30:1"),
    (6, "pd_acac2", "xantphos", 2, 2, 96, "32:1"),
    (7, "pd_acac2", "xantphos", 1, 1, 65, "13:1"),
    (8, "pd_oac2", "xantphos", 2, 2, 44, "3:1"),
    (9, "pdcl2", "xantphos", 2, 2, "trace", "table dash (–)"),
    (10, "pd_me_cn2_cl2", "xantphos", 2, 2, "trace", "table dash (–)"),
]

for entry, catalyst, ligand, cat_loading, ligand_loading, yld, rr in optimization:
    r = common_reaction(
        f"table1_entry_{entry:02d}", EV_TABLE, "specific_experiment", ["optimization"],
        {"page": 3, "table": "Table 1", "entry": entry},
        "optimization of hydrosilylation of 1-octylallene and PhSiH3",
    )
    r["general_scheme_ids"] = [generic["reaction_id"]]
    r["participants"] = [
        participant("allene_1a", mol_id("1a"), "reactant", EV_TABLE,
                    charges={"amount": ("0.25 mmol", q(0.25, "mmol"), EV_TABLE)}),
        participant("phenylsilane", mol_id("phenylsilane"), "reactant", EV_TABLE,
                    charges={"amount": ("0.3 mmol", q(0.3, "mmol"), EV_TABLE)}),
        participant("catalyst", mol_id(catalyst), "catalyst", EV_TABLE,
                    charges={"mol_percent": (f"{cat_loading} mol%", q(cat_loading, "mol%"), EV_TABLE)},
                    reference=mol_id("1a")),
        participant("ligand", mol_id(ligand), "ligand", EV_TABLE,
                    charges={"mol_percent": (f"{ligand_loading} mol%", q(ligand_loading, "mol%"), EV_TABLE)},
                    reference=mol_id("1a")),
        participant("thf", mol_id("thf"), "solvent", EV_TABLE,
                    charges={"volume": ("1 mL", q(1, "mL"), EV_TABLE)}),
        participant("product_2a", mol_id("2a"), "product", EV_TABLE, product_role="branched allylsilane 2a"),
        participant("other_isomers", mol_id("other_2a_isomers"), "product", EV_TABLE,
                    product_role="combined other hydrosilylation isomers"),
        participant("internal_standard", mol_id("tmob"), "additive", EV_TABLE, stage="characterization"),
    ]
    r["reaction_scale"] = source(
        "1-octylallene (0.25 mmol)", {"reference_molecule_id": mol_id("1a"), "amount": q(0.25, "mmol")}, EV_TABLE
    )
    set_solvent(r, EV_TABLE, ("THF (1 mL)", q(1, "mL")))
    r["steps"] = [step("reaction", EV_TABLE, [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"],
                       "hydrosilylation", ("30 °C", q(30, "°C")), ("1 h", q(1, "h")), "N2")]
    r["procedure"]["original_text"] = source(
        "1-octylallene (0.25 mmol), phenylsilane (0.3 mmol), L, Pd source, THF (1 mL), 30 °C, 1 h, N2 atmosphere.",
        "Table 1 footnote a conditions", EV_TABLE
    )
    if isinstance(yld, int):
        r["outcomes"].append(outcome(
            f"entry{entry}_yield", "yield", [mol_id("2a")], EV_TABLE,
            f"{yld}", q(yld, "%"), "branched allylsilane 2a", "1H NMR yield",
            "1H NMR spectroscopy with 1,3,5-trimethoxybenzene as an internal standard",
            {"table_column": "2a/yield%"},
        ))
    else:
        r["outcomes"].append(outcome(
            f"entry{entry}_trace", "trace", [mol_id("2a")], EV_TABLE,
            "trace", "trace", "branched allylsilane 2a", "qualitative entry in the NMR-yield column",
            "1H NMR spectroscopy with 1,3,5-trimethoxybenzene as an internal standard",
            {"table_column": "2a/yield%"},
        ))
    r["outcomes"].append(outcome(
        f"entry{entry}_rr", "selectivity", [mol_id("2a"), mol_id("other_2a_isomers")], EV_TABLE,
        rr, {"ratio": rr, "order": [mol_id("2a"), mol_id("other_2a_isomers")]},
        "regioisomeric outcome", None, "1H NMR spectroscopy",
        {"metric": "r.r.", "numerator": "2a", "denominator": "other isomers"},
    ))
    r["review"]["unresolved"] = [
        "Individual structures and quantities of the combined 'other isomers' are not reported.",
        "The internal-standard amount and characterization timing are not reported.",
        "Workup and purification are not reported for optimization rows.",
    ]
    add_hydrosilylation_classification(r, EV_TABLE)
    reactions.append(r)

# Scheme 2: thirteen distinct substrate-scope experiments sharing the caption and page-4 general procedure.
for product_label, substrate_label, _, _, _, _, yld in scope_rows:
    evidence = [P3I, P3T, P4I, P4T]
    r = common_reaction(
        f"scheme2_{product_label}", evidence, "specific_experiment", ["substrate_scope"],
        {"pages": [3, 4], "scheme": "Scheme 2", "product": product_label, "general_procedure_footnote": 1},
        "Pd(acac)2/Xantphos-catalyzed regioselective hydrosilylation of allenes",
    )
    r["general_scheme_ids"] = [generic["reaction_id"]]
    substrate_mid = mol_id(substrate_label)
    product_mid = mol_id(product_label)
    other_mid = mol_id(f"other_{product_label}_isomers")
    r["participants"] = [
        participant("allene", substrate_mid, "reactant", evidence,
                    charges={"amount": ("allene (1 mmol)", q(1, "mmol"), [P3I, P3T])}),
        participant("phenylsilane", mol_id("phenylsilane"), "reactant", evidence,
                    charges={"mass": ("PhSiH3 (130 mg)", q(130, "mg"), [P4I, P4T]),
                             "amount": ("PhSiH3 (1.2 mmol)", q(1.2, "mmol"), evidence)}),
        participant("pd_acac2", mol_id("pd_acac2"), "catalyst", evidence,
                    charges={"mass": ("Pd(acac)2 (6 mg)", q(6, "mg"), [P4I, P4T])},
                    reference=substrate_mid),
        participant("xantphos", mol_id("xantphos"), "ligand", evidence,
                    charges={"mass": ("Xantphos (12 mg)", q(12, "mg"), [P4I, P4T]),
                             "mol_percent": ("Xantphos (2 mol%)", q(2, "mol%"), evidence)},
                    reference=substrate_mid),
        participant("thf_initial", mol_id("thf"), "solvent", evidence,
                    charges={"volume": ("THF (1 mL) added to Pd(acac)2/Xantphos", q(1, "mL"), [P4I, P4T])}),
        participant("thf_with_substrates", mol_id("thf"), "solvent", evidence,
                    charges={"volume": ("allene and PhSiH3 in THF (1 mL)", q(1, "mL"), [P4I, P4T])}),
        participant("product", product_mid, "product", evidence, product_role=f"branched allylsilane {product_label}"),
        participant("other_isomers", other_mid, "product", evidence, product_role="combined other hydrosilylation isomers"),
        participant("silica", mol_id("silica_gel"), "purification_agent", [P4I, P4T], stage="purification"),
        participant("petroleum_ether", mol_id("petroleum_ether"), "purification_agent", [P4I, P4T], stage="purification"),
        participant("ethyl_acetate", mol_id("ethyl_acetate"), "purification_agent", [P4I, P4T], stage="purification"),
    ]
    r["participants"][2]["charges"]["mol_percent"] = conflict_mol_percent(evidence)
    r["reaction_scale"] = source(
        "allenes (1 mmol)", {"reference_molecule_id": substrate_mid, "amount": q(1, "mmol")}, [P3I, P3T]
    )
    set_solvent(r, evidence, ("THF (2 mL)", q(2, "mL")))
    r["steps"] = [
        step("reaction", evidence, [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"],
             "stirred until completion; monitored by TLC", ("30 °C", q(30, "°C")),
             ("0.5–3 h", {"minimum": 0.5, "maximum": 3, "unit": "h", "raw_text": "0.5–3 h"}), "N2"),
        step("workup_remove_solvent", [P4I, P4T], [], "solvent was removed in vacuum", stage="workup"),
        step("purification_column", [P4I, P4T], ["silica", "petroleum_ether", "ethyl_acetate"],
             "silica gel column chromatography eluting with petroleum ether and ethyl acetate", stage="purification"),
    ]
    r["procedure"]["original_text"] = source(
        "In a nitrogen-filled Schlenk tube, Pd(acac)2 (6 mg, 2 mmol%), Xantphos (12 mg, 2 mol%) were added to THF (1 mL), followed by allene (1 mmol) and PhSiH3 (130 mg, 1.2 mmol) in THF (1 mL); stirred at 30 °C for 0.5–3 h, monitored by TLC; solvent removed in vacuum; silica chromatography (petroleum ether/ethyl acetate).",
        "General procedure footnote 1", [P4I, P4T]
    )
    r["procedure"]["general_procedure_reference"] = source("General procedure for hydrosilylation of allenes (footnote 1)", "footnote 1", [P4I, P4T])
    r["procedure"]["addition_order"] = source(
        "Pd(acac)2 and Xantphos were added to THF, followed by allene and PhSiH3 in THF.",
        ["Pd(acac)2", "Xantphos", "THF", "allene", "PhSiH3 in THF"], [P4I, P4T]
    )
    r["outcomes"] = [
        outcome(f"{product_label}_yield", "yield", [product_mid], [P3I, P3T],
                f"{yld}%", q(yld, "%"), f"branched allylsilane {product_label}", "isolated yield",
                "Yields of isolated products", {"scheme": "Scheme 2", "product": product_label}),
        outcome(f"{product_label}_rr", "selectivity", [product_mid, other_mid], [P3I, P3T],
                ">15:1", {"ratio": ">15:1", "qualifier": "gt", "order": [product_mid, other_mid]},
                "branched allylsilane versus all other isomers", None, "1H NMR spectroscopy",
                {"metric": "r.r.", "numerator": product_label, "denominator": "all other isomers"}),
    ]
    r["notes"] = [source(
        "Scheme 2 caption reports Pd(acac)2 (2 mol%), while the printed page-4 general procedure reads Pd(acac)2 (6 mg, 2 mmol%).",
        "Source loading-unit conflict retained without correction.", evidence
    )]
    r["review"]["unresolved"] = [
        "The exact product-specific time within the inherited 0.5–3 h range is not reported.",
        "The printed Pd(acac)2 loading unit conflicts between Scheme 2 (mol%) and the general procedure (mmol%).",
        "The petroleum ether/ethyl acetate ratio and isolated product mass are unreported.",
        "The structures and quantities of 'all other isomers' are not enumerated.",
        "Candidate SMILES are machine transcriptions of the displayed structures and require expert structure review.",
        "The newly formed tetrahedral product center has no assigned configuration; no racemic composition is inferred.",
    ]
    add_hydrosilylation_classification(r, [P3I, P4I])
    reactions.append(r)

# Scheme 3 is a proposed mechanism, not an experimentally demonstrated sequence.
mechanism = common_reaction(
    "scheme3_proposed_mechanism", [P4I, P4T], "proposed_mechanism", ["other"],
    {"page": 4, "scheme": "Scheme 3"}, "proposed Pd(0)-involved mechanism",
)
mechanism["participants"] = [
    participant("pd_precursor", mol_id("pd_acac2"), "catalyst", [P4I, P4T]),
    participant("ligand", mol_id("xantphos"), "ligand", [P4I, P4T]),
    participant("phenylsilane", mol_id("phenylsilane"), "reactant", [P4I, P4T]),
    participant("allene", mol_id("generic_allene"), "reactant", [P4I, P4T]),
    participant("intermediate_1", mol_id("mech1"), "catalyst", [P4I, P4T]),
    participant("intermediate_2", mol_id("mech2"), "catalyst", [P4I, P4T]),
    participant("intermediate_3", mol_id("mech3"), "catalyst", [P4I, P4T]),
    participant("intermediate_4", mol_id("mech4"), "catalyst", [P4I, P4T]),
    participant("product", mol_id("generic_branched_allylsilane"), "product", [P4I, P4T], product_role="branched allylsilane"),
]
mechanism["steps"] = [
    step("activation", [P4I, P4T], ["pd_precursor", "ligand", "intermediate_1"], "Pd(acac)2 was proposed to be reduced to (L)Pd(0)"),
    step("oxidative_addition", [P4I, P4T], ["intermediate_1", "phenylsilane", "intermediate_2"], "oxidative addition of phenylsilane to Pd(0) was proposed"),
    step("allene_coordination", [P4I, P4T], ["intermediate_2", "allene", "intermediate_3"], "allene coordination was proposed"),
    step("hydrometalation", [P4I, P4T], ["intermediate_3", "intermediate_4"], "hydrometalation to an allylpalladium intermediate was proposed"),
    step("reductive_elimination", [P4I, P4T], ["intermediate_4", "product", "intermediate_1"], "C-Si reductive elimination and catalyst regeneration were proposed"),
]
mechanism["notes"] = [source(
    "Further detailed mechanistic studies are ongoing to establish unambiguously the real mechanistic nature of the reaction.",
    "Authors explicitly describe Scheme 3 as proposed rather than established.", [P4I, P4T]
)]
mechanism["classification"]["coverage"] = derived(
    "unresolved", [P4I, P4T], "The drawing is a proposed elementary cycle with unresolved exact coordination and atom mapping.", method="inference"
)
mechanism["classification"]["limitations"] = ["Proposed mechanism only; no elementary step is promoted to an experimentally verified reaction."]
mechanism["review"]["unresolved"] = [
    "Exact coordination structures, oxidation states, atom mapping, and stereochemical course are unresolved.",
    "No kinetic, isotopic, or intermediate-observation experiment is reported on the supplied pages.",
]
reactions.append(mechanism)

bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": f"{DOC}:front_matter", "document_id": DOC, "context_type": "bibliographic_metadata",
     "evidence_ids": [f"{DOC}:page1:image", f"{DOC}:page1:text"],
     "note": "PDF page 1 is publisher front matter; it contains no reaction experiment."},
    {"context_id": f"{DOC}:scheme1_prior_work", "document_id": DOC, "context_type": "literature_prior_work",
     "evidence_ids": [P3I, P2T],
     "note": "Scheme 1a-b summarizes cited prior hydrosilylation selectivities and ligand systems. It is not counted as a current-paper experiment because individual substrates/conditions/results are not supplied there."},
    {"context_id": f"{DOC}:general_procedure", "document_id": DOC, "context_type": "general_procedure",
     "evidence_ids": [P4I, P4T],
     "note": "The footnote procedure is inherited by the thirteen Scheme 2 scope records; it is not duplicated as a separate concrete experiment."},
    {"context_id": f"{DOC}:si_notice", "document_id": DOC, "context_type": "supplementary_material_notice",
     "evidence_ids": [f"{DOC}:page1:image", P2T, P4T],
     "note": "The article advertises supplemental data and NMR spectra, but no SI file was provided or explicitly associated in this workspace."},
    {"context_id": f"{DOC}:references_page5", "document_id": DOC, "context_type": "references_only",
     "evidence_ids": [f"{DOC}:page5:image", f"{DOC}:page5:text"],
     "note": "PDF page 5 contains references only and contributes no current-paper experiment."},
]

diagnostics = {
    m["molecule_id"]: st.stereo_inspect(m["smiles"]["normalized_value"])
    for m in bundle["molecules"]
    if m["representation_type"] == "smiles"
}
st.write(ROOT / "runs/current/raw_model/stereo_diagnostics.json", {
    "document_id": DOC,
    "source_verified": False,
    "note": "Generated with the v3.1 standard stereo implementation; graph diagnostics do not verify source identity.",
    "molecules": diagnostics,
})
st.write(ROOT / "runs/current/raw_model/source_mentions.json", st.source_mentions(
    dataset, ["1a", "2a", "2e", "2m", "general procedure", "r.r."], context=320
))

out = ROOT / "runs/current/submitted.json"
st.write(out, bundle)
print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
                  "stereo_diagnostics": len(diagnostics), "output": str(out)}, ensure_ascii=False))
