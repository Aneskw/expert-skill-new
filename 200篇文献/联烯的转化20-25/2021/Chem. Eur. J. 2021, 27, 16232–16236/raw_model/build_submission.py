#!/usr/bin/env python3
"""Build the v3.1 machine submission from visually inspected pages 1--5.

This is the retained raw extraction program.  It deliberately leaves structures,
procedural details, and measurement bases unresolved where the main article does
not supply them; no SI was provided or associated.
"""
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_2f1500a108e364c7"
DATASET = ROOT / "runs/current/dataset.json"
OUT = ROOT / "runs/current/submitted.json"
EV = {n: f"{DOC}:page{n}:image" for n in range(1, 6)}


def f(raw=None, norm=None, ev=(), method="original_extraction", derivation=None,
      status=None, verification=None, alternatives=()):
    return st.fact(raw, norm, ev, method=method, derivation=derivation,
                   status=status, verification=verification,
                   alternatives=alternatives)


def inferred(norm, ev, why, raw=None):
    return f(raw, norm, ev, method="inference", derivation=why)


def derived(norm, ev, why, raw=None):
    return f(raw, norm, ev, method="structure_derivation", derivation=why)


def unreadable(raw, ev, why):
    return f(raw, None, ev, status="unreadable", derivation=why)


def q(number, unit, raw=None, qualifier=None):
    x = {"number": number, "unit": unit}
    if raw is not None:
        x["raw_text"] = raw
    if qualifier:
        x["qualifier"] = qualifier
    return x


bundle = st.empty(json.loads(DATASET.read_text())["documents"])
molecules = {}


def base_molecule(mid, label, ev, name=None, chirality="unresolved",
                  geometry="unresolved", representation="unresolved",
                  stereo_note=None):
    ev = list(dict.fromkeys(ev))
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    m["names"]["original"] = f(name or label, name or label, ev)
    m["names"]["english"] = f(name or label, name or label, ev)
    m["cas_numbers"] = f()
    m["identity_note"] = f(
        "Identity is limited to the displayed label/name and inspected main-text drawing; no SI was supplied.",
        "Identity is source-identified but the structure remains a machine candidate.", ev)
    cue = stereo_note or "Stereochemistry was inspected in the cited drawing; no stronger assignment is asserted."
    m["stereochemistry"]["chirality_status"] = inferred(chirality, ev, cue, label)
    m["stereochemistry"]["geometry_status"] = inferred(geometry, ev, cue, label)
    m["stereochemistry"]["representation_status"] = inferred(representation, ev, cue, label)
    molecules[mid] = m
    return m


def unresolved_molecule(mid, label, ev, name=None, kind="ordinary", mixture=None):
    if kind == "allene":
        m = base_molecule(mid, label, ev, name, "unspecified", "not_applicable", "not_encoded",
                          "The allene axis was inspected; the source drawing does not assign R_a/S_a or P/M, and no axial stereo is encoded.")
        m["stereochemistry"]["axial"] = f(
            "Substituted C=C=C axis shown without an absolute assignment.",
            [{"axis": "allene", "configuration": "unspecified"}], ev)
    elif kind == "cis_product":
        m = base_molecule(mid, label, ev, name, "partially_specified", "not_applicable", "not_encoded",
                          "Wedge/dash drawing and text assign the cis relative product; absolute configuration is not assigned, so a concrete stereochemical SMILES is withheld.")
        m["stereochemistry"]["atom_centers"] = f(
            "cis-adduct shown by wedge/dash bonds; absolute centers not named",
            [{"assignment": "cis relative configuration", "absolute_configuration": "unassigned"}], ev)
        m["stereochemistry"]["assignment"] = f(
            "complete diastereoselectivity for the cis-adduct",
            {"type": "relative", "assignment": "cis", "absolute": "unassigned"}, ev)
    elif kind == "enriched_mixture":
        m = base_molecule(mid, label, ev, name, "mixture", "not_applicable", "not_encoded",
                          "The sample is reported only as an unassigned 67:33 enantiomeric mixture; no absolute configuration is inferred.")
        m["stereochemistry"]["mixture"] = f("67:33 er", {"metric": "er", "ratio": "67:33", "configuration_order": "unassigned"}, ev)
        m["stereochemistry"]["assignment"] = f("absolute configuration not assigned", "unassigned", ev)
    else:
        m = base_molecule(mid, label, ev, name, "unresolved", "unresolved", "unresolved")
    if mixture is not None:
        m["stereochemistry"]["mixture"] = f(mixture, mixture, ev)
    return m


def smiles_molecule(mid, label, name, smiles, ev, chirality, geometry, representation,
                    stereo_reason, axial=None):
    m = base_molecule(mid, label, ev, name, chirality, geometry, representation, stereo_reason)
    diag = st.stereo_inspect(smiles)
    m["representation_type"] = "smiles"
    m["smiles"] = derived(smiles, ev, "Candidate connectivity was transcribed from the inspected source drawing; RDKit parsing is not source verification.", smiles)
    m["stereochemistry"]["isomeric_smiles"] = derived(diag["canonical_isomeric_smiles"], ev, "RDKit canonical isomeric SMILES generated from the candidate source-derived graph.")
    m["stereochemistry"]["identity_key"] = derived(diag["stereo_graph_key"], ev, "RDKit stereo_graph_key generated from the candidate graph; it is not a sample identity.")
    if axial is not None:
        m["stereochemistry"]["axial"] = f(axial, axial, ev)
    return m


def generic_molecule(mid, label, fixed, variables, ev):
    m = base_molecule(mid, label, ev, label, "unresolved", "unresolved", "not_encoded",
                      "The figure is an R-group template; individual center/axis/geometry assignments cannot be made without an instance.")
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": f(fixed, {"description": fixed}, ev),
        "variables": [{
            "label": v,
            "attachment_points": f(f"attachment shown for {v}", f"source drawing attachment for {v}", ev),
            "definition": f(f"{v} is an undefined substituent in Figure 1D", None, ev, status="unreadable"),
            "allowed_values": f(),
            "constraints": f(),
        } for v in variables],
        "instance_mapping": [],
        "reason": f("Figure 1D uses EWG, Pg and R variables.", "Generic R-group drawing cannot be represented as a concrete molecule.", ev),
    }
    return m


def ensure(mid, label=None, ev=None, **kw):
    if mid not in molecules:
        unresolved_molecule(mid, label or mid, ev or [EV[1]], **kw)
    return mid


# Key source-derived connectivity candidates and common named agents.
smiles_molecule("m_10a", "10a", "amino crotonate 10a", "CCOC(=O)C=CCNS(=O)(=O)c1ccc(C)cc1", [EV[1], EV[2], EV[3], EV[4]],
                "unspecified", "unspecified", "partial",
                "No absolute chirality applies to the candidate graph; the potentially stereogenic C=C geometry is drawn but not named E/Z, so it remains unspecified.")
smiles_molecule("m_12a", "12a", "allenoate 12a", "O=C(OCc1ccccc1)C=C=Cc1ccccc1", [EV[2], EV[4]],
                "unspecified", "not_applicable", "partial",
                "The substituted allene can possess axial chirality, but Figure 1/Table 1/Scheme 1 gives no R_a/S_a or P/M assignment; toolkit double-bond flags do not encode the axis.",
                {"axis": "C=C=C", "configuration": "unspecified", "representation": "not encoded"})

common_smiles = {
    "m_pme2ph": ("P(Me)2Ph", "dimethylphenylphosphine", "CP(C)c1ccccc1"),
    "m_pmeph2": ("PMe(Ph)2", "methyldiphenylphosphine", "CP(c1ccccc1)c1ccccc1"),
    "m_pbu3": ("PBu3", "tributylphosphine", "CCCCP(CCCC)CCCC"),
    "m_pph3": ("PPh3", "triphenylphosphine", "P(c1ccccc1)(c1ccccc1)c1ccccc1"),
    "m_pcy3": ("P(C6H11)3", "tricyclohexylphosphine", "P(C1CCCCC1)(C1CCCCC1)C1CCCCC1"),
    "m_pnme2_3": ("P(NMe2)3", "hexamethylphosphorous triamide", "CN(C)P(N(C)C)N(C)C"),
    "m_ch2cl2": ("CH2Cl2", "dichloromethane", "ClCCl"),
    "m_dce": ("(CH2Cl)2", "1,2-dichloroethane", "ClCCCl"),
    "m_chlorobenzene": ("ClC6H5", "chlorobenzene", "c1ccccc1Cl"),
    "m_thf": ("THF", "tetrahydrofuran", "C1CCOC1"),
    "m_etoh": ("EtOH", "ethanol", "CCO"),
    "m_h2": ("H2", "hydrogen", "[H][H]"),
    "m_n_buli": ("n-BuLi", "n-butyllithium", "[Li]CCCC"),
    "m_hexanethiol": ("C6H13SH", "1-hexanethiol", "CCCCCCS"),
}
for mid, (label, name, smi) in common_smiles.items():
    ev = [EV[2], EV[3], EV[4]] if mid == "m_pme2ph" else ([EV[4]] if mid in {"m_etoh", "m_h2", "m_n_buli", "m_hexanethiol"} else [EV[2]])
    smiles_molecule(mid, label, name, smi, ev, "achiral", "not_applicable", "encoded",
                    "The candidate molecular graph has no potential tetrahedral, axial, or E/Z stereogenic element in the toolkit check; this is a structure-derived, not source-certified, conclusion.")

# Products and source-labelled materials.
unresolved_molecule("m_13a", "13a", [EV[1], EV[2], EV[3]], "pyrroline 13a", kind="cis_product")
unresolved_molecule("m_15a", "15a", [EV[2], EV[3], EV[4]], "pyrrole/by-product 15a", kind="cis_product")
unresolved_molecule("m_13a_er", "13a (67:33 er sample)", [EV[4]], "enantioenriched pyrroline 13a sample", kind="enriched_mixture")
unresolved_molecule("m_10d", "(D2)-10d", [EV[2]], "dideuterated amino crotonate (D2)-10d")
unresolved_molecule("m_14a", "14a", [EV[2]], "imidazolidine 14a", kind="cis_product")
unresolved_molecule("m_14d", "(D4)-14d", [EV[2]], "tetradeuterated imidazolidine (D4)-14d", kind="cis_product")
unresolved_molecule("m_D3", "D3", [EV[4]], "chiral phosphepine D3", kind="ordinary")
unresolved_molecule("m_pdc", "10% Pd/C", [EV[4]], "10% palladium on carbon")

# Table 2 source products and partners. Concrete SMILES are withheld because the
# compact main-text drawings do not safely preserve every ester/protecting group
# and relative stereochemical detail without the SI.
for label in [*[f"13{x}" for x in "bcdefghijklmnopqrstuvwxy"], *[f"15{x}" for x in "bcdefghijklmnopqrstuvwxy"]]:
    unresolved_molecule(f"m_{label}", label, [EV[3]], label, kind="cis_product")
for label in [*[f"10{x}" for x in "bcdefghijk"], *[f"12{x}" for x in "bcdefghijklmnopqrstuvwxy"]]:
    unresolved_molecule(f"m_{label}", label, [EV[3]], label, kind="allene" if label.startswith("12") else "ordinary")
for label in ("13m", "13t"):
    # Add page-4 occurrences to the already-created registry entries.
    molecules[f"m_{label}"]["occurrences"].append({"document_id": DOC, "compound_label": f(label, label, [EV[4]]), "evidence_ids": [EV[4]]})
unresolved_molecule("m_16m", "16m", [EV[4]], "thiolate addition product 16m", kind="cis_product")
unresolved_molecule("m_17t", "17t", [EV[4]], "hydrogenation product 17t", kind="cis_product")


def participant(pid, mid, role, ev, product_role=None, mol_percent=None, ref=None, stage="reaction"):
    p = st.new_participant(pid, mid, role, ev, stage=stage)
    if product_role:
        p["product_role"] = f(product_role, product_role, ev)
    if mol_percent is not None:
        p["charges"]["mol_percent"] = f(f"{mol_percent} mol%", q(mol_percent, "mol%", f"{mol_percent} mol%"), ev)
        p["equivalence_reference"] = f(ref, ref, ev)
    return p


def step(step_id, ev, pids, temp=None, time=None, operation=None, atmosphere=None):
    s = st.new_step(step_id, ev)
    s["participant_ids"] = list(pids)
    if operation is not None:
        s["operation"] = f(operation, operation, ev)
    if temp is not None:
        raw, norm = temp
        s["temperature"] = f(raw, norm, ev)
    if time is not None:
        raw, norm = time
        s["time"] = f(raw, norm, ev)
    if atmosphere is not None:
        s["atmosphere"] = f(atmosphere, atmosphere, ev)
    return s


def outcome(oid, kind, targets, ev, raw=None, norm=None, role=None, ytype=None, basis=None, scope=None):
    o = st.new_outcome(oid, kind, targets, ev)
    if raw is not None or norm is not None:
        o["value"] = f(raw, norm, ev)
    if role is not None:
        o["product_role"] = f(role, role, ev)
    if ytype is not None:
        o["yield_type"] = f(ytype, ytype, ev)
    if basis is not None:
        o["measurement_basis"] = f(basis, basis, ev)
    if scope is not None:
        o["scope"] = f(scope, scope, ev)
    return o


def reaction(rid, ev, rtype="specific_experiment", categories=None, location=None, name=None, transform=None):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = rtype
    r["experiment_categories"] = f(categories or [], categories or [], ev)
    r["location"] = f(location or rid, {"description": location or rid}, ev)
    if name:
        r["reaction_name"] = f(name, name, ev)
    if transform:
        r["transformation_type"] = inferred(transform, ev, "Transformation family assigned from the source caption and displayed reactant/product connectivity.", transform)
    r["review"]["unresolved"] = [
        "No Supporting Information was supplied or associated; charges, scale, atmosphere, pressure, workup and purification details absent from the main article remain unreported.",
        "Machine-extracted record requires expert confirmation of compact chemical drawings and stereochemical correspondence."
    ]
    return r


def solvent(r, mid, label, ev):
    r["participants"].append(participant("solvent", mid, "solvent", ev))
    r["solvent_conditions"]["composition"] = f(label, [label], ev)


def classification_unresolved(r, ev, why):
    r["classification"]["coverage"] = inferred("unresolved", ev, why)
    r["classification"]["limitations"] = [why]


def classification_annulation(r, ev):
    changes = [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "N"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
        {"change_id": "b4", "atoms": [{"atom_map": None, "element": "N"}, {"atom_map": None, "element": "H"}], "before_order": 1, "after_order": 0, "change_type": "hydrogen_loss", "hydrogen_count": 1},
    ]
    labels = [
        {"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1"]},
        {"category_number": 5, "code": "C-N", "supporting_change_ids": ["b2"]},
        {"category_number": 1, "code": "C-H", "supporting_change_ids": ["b3"]},
    ]
    r["classification"]["bond_changes"] = [derived(x, ev, "Net bond-change class inferred from the displayed redox-isomerization/(3+2) annulation; exact atom mapping is unresolved.") for x in changes]
    r["classification"]["labels"] = derived(labels, ev, "Taxonomy labels calculated from the explicitly listed inferred bond changes.")
    r["classification"]["coverage"] = inferred("partial", ev, "No atom-mapped complete structures are asserted; hydrogen source and exact atom correspondence remain unresolved.")
    r["classification"]["limitations"] = ["Partial net classification only; no atom-mapped reaction SMILES.", "Hydrogen gain/loss is net and does not identify a proton-transfer mechanism."]


def classification_thiol(r, ev):
    changes = [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "S"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
        {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
    ]
    labels = [{"category_number": 4, "code": "C-S", "supporting_change_ids": ["b1"]}, {"category_number": 1, "code": "C-H", "supporting_change_ids": ["b2"]}]
    r["classification"]["bond_changes"] = [derived(x, ev, "Net conjugate thiolate addition inferred from Figure 3A; atom mapping is not asserted.") for x in changes]
    r["classification"]["labels"] = derived(labels, ev, "Taxonomy labels calculated from the listed inferred changes.")
    r["classification"]["coverage"] = inferred("partial", ev, "Exact atom mapping and proton source are unresolved.")
    r["classification"]["limitations"] = ["Partial net graph change only."]


def classification_hydrogenation(r, ev):
    changes = [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 2},
    ]
    labels = [{"category_number": 1, "code": "C-H", "supporting_change_ids": ["b2"]}]
    r["classification"]["bond_changes"] = [derived(x, ev, "Net alkene hydrogenation inferred from Figure 3A; no atom mapping is asserted.") for x in changes]
    r["classification"]["labels"] = derived(labels, ev, "A reduced C-C double bond is not classified as a newly formed C-C sigma bond; only C-H formation is indexed.")
    r["classification"]["coverage"] = inferred("partial", ev, "Implicit hydrogen mapping is net only.")
    r["classification"]["limitations"] = ["Hydrogen atom identities are unresolved; reduced C-C sigma bond is not a new connection."]


reactions = []

# Figure 1: three prior-work examples and this paper's generic design.
for mid, label in [("m_lit_1", "1"), ("m_lit_2", "2"), ("m_lit_3", "3"), ("m_lit_5", "5"), ("m_lit_6", "6"), ("m_lit_7", "7"), ("m_lit_8", "8")]:
    unresolved_molecule(mid, label, [EV[1]], f"literature-example compound {label}")
unresolved_molecule("m_rhco", "Rh[I]/CO", [EV[1]], "Rh(I)/carbon monoxide catalyst system")
unresolved_molecule("m_mecn", "MeCN", [EV[1]], "acetonitrile")
unresolved_molecule("m_nhc", "NHC", [EV[1]], "N-heterocyclic carbene catalyst")
unresolved_molecule("m_pme3", "PMe3", [EV[1]], "trimethylphosphine")

r = reaction("figure1A_literature_example", [EV[1]], "literature_example", [], "Figure 1A", "Rh[I] catalyzed allylation/Pauson-Khand reaction", "assisted auto-tandem allylation/Pauson-Khand")
r["participants"] = [participant("r2", "m_lit_2", "reactant", [EV[1]]), participant("r3", "m_lit_3", "reactant", [EV[1]]), participant("cat", "m_rhco", "catalyst", [EV[1]], mol_percent=10, ref="m_lit_2"), participant("solvent", "m_mecn", "solvent", [EV[1]]), participant("p1", "m_lit_1", "product", [EV[1]], "main product")]
r["solvent_conditions"]["composition"] = f("MeCN", ["MeCN"], [EV[1]])
r["steps"] = [
    step("cycle_A", [EV[1]], ["r2", "r3", "cat", "solvent"], ("30 °C", q(30, "°C", "30 °C")), operation="allylic allylation"),
    step("cycle_B", [EV[1]], ["r2", "r3", "cat", "solvent", "p1"], ("reflux", "reflux"), operation="temperature increased; Pauson-Khand cycle; ramp time unreported"),
]
r["review"]["unresolved"].append("Figure 1A gives 30 °C then reflux, but no duration for either stage or the heating ramp.")
r["outcomes"] = [outcome("yield", "yield", ["m_lit_1"], [EV[1]], "87% y", q(87, "%", "87% y"), "main product", basis=None), outcome("dr", "selectivity", ["m_lit_1"], [EV[1]], "7:1 dr", "7:1", scope="diastereomeric ratio")]
classification_unresolved(r, [EV[1]], "Prior-work illustration is not atom-mapped and is not classified as a new experiment in this paper.")
reactions.append(r)

r = reaction("figure1B_literature_example", [EV[1]], "literature_example", [], "Figure 1B", "NHC catalyzed ring opening/cycloisomerization", "auto-tandem ring opening/cycloisomerization")
r["participants"] = [participant("r5", "m_lit_5", "reactant", [EV[1]]), participant("cat", "m_nhc", "catalyst", [EV[1]], mol_percent=10, ref="m_lit_5"), participant("solvent", "m_ch2cl2", "solvent", [EV[1]]), participant("p6", "m_lit_6", "product", [EV[1]], "main product")]
r["solvent_conditions"]["composition"] = f("CH2Cl2", ["CH2Cl2"], [EV[1]])
r["steps"] = [step("reaction", [EV[1]], ["r5", "cat", "solvent", "p6"], time=("16 h", q(16, "h", "16 h")), operation="electrocyclic opening followed by cycloisomerization")]
r["outcomes"] = [outcome("yield", "yield", ["m_lit_6"], [EV[1]], "70% y", q(70, "%", "70% y"), "main product")]
classification_unresolved(r, [EV[1]], "Prior-work illustration is not atom-mapped.")
reactions.append(r)

r = reaction("figure1C_literature_example", [EV[1]], "literature_example", [], "Figure 1C", "Phosphine redox isomerization", "alkynoate-to-dienoate redox isomerization")
r["participants"] = [participant("r7", "m_lit_7", "reactant", [EV[1]]), participant("cat", "m_pme3", "catalyst", [EV[1]], mol_percent=10, ref="m_lit_7"), participant("solvent", "m_ch2cl2", "solvent", [EV[1]]), participant("p8", "m_lit_8", "product", [EV[1]], "main product")]
r["solvent_conditions"]["composition"] = f("CH2Cl2", ["CH2Cl2"], [EV[1]])
r["steps"] = [step("reaction", [EV[1]], ["r7", "cat", "solvent", "p8"], time=("16 h", q(16, "h", "16 h")), operation="redox isomerization")]
classification_unresolved(r, [EV[1]], "Prior-work illustration is not atom-mapped; hydrogen transfer and double-bond changes remain unresolved.")
reactions.append(r)

generic_molecule("m_10_generic", "10", "amino crotonate framework bearing EWG1 and Pg", ["EWG1", "Pg"], [EV[1]])
generic_molecule("m_12_generic", "12", "allene framework bearing EWG2 and R", ["EWG2", "R"], [EV[1]])
generic_molecule("m_13_generic", "13", "pyrroline framework bearing EWG1, EWG2, Pg and R", ["EWG1", "EWG2", "Pg", "R"], [EV[1]])
r = reaction("figure1D_generic_scheme", [EV[1]], "generic_scheme", [], "Figure 1D", "Redox isomerization/(3+2) annulation by auto-tandem phosphine catalysis", "auto-tandem redox isomerization/(3+2) annulation")
r["participants"] = [participant("r10", "m_10_generic", "reactant", [EV[1]]), participant("r12", "m_12_generic", "reactant", [EV[1]]), participant("cat", "m_pme2ph", "catalyst", [EV[1]], mol_percent=10, ref="m_10_generic"), participant("solvent", "m_ch2cl2", "solvent", [EV[1]]), participant("p13", "m_13_generic", "product", [EV[1]], "main product")]
r["solvent_conditions"]["composition"] = f("CH2Cl2", ["CH2Cl2"], [EV[1]])
r["steps"] = [step("reaction", [EV[1]], ["r10", "r12", "cat", "solvent", "p13"], time=("16 h", q(16, "h", "16 h")), operation="redox isomerization followed by (3+2) annulation")]
r["outcomes"] = [outcome("scope_summary", "other", ["m_13_generic"], [EV[1]], "19 examples", 19, scope="generic-scheme scope summary; count of examples, not a chemical amount"), outcome("dr", "selectivity", ["m_13_generic"], [EV[1]], "all >20:1 dr", {"ratio": ">20:1", "qualifier": "gt"}, scope="diastereomeric ratio")]
classification_annulation(r, [EV[1]])
reactions.append(r)

# Figure 2A: one experiment per catalyst, including three individually represented dash results.
fig2a = [("pme2ph", "m_pme2ph", "27 (12)^b"), ("pmeph2", "m_pmeph2", "45"), ("pbu3", "m_pbu3", "19"), ("pph3", "m_pph3", "–"), ("pcy3", "m_pcy3", "–"), ("pnme2_3", "m_pnme2_3", "–")]
for key, cat_mid, yraw in fig2a:
    r = reaction(f"figure2A_{key}", [EV[2]], "specific_experiment", ["mechanistic_study"], f"Figure 2A catalyst row: {key}", "Redox isomerization/annulation of amino crotonate 10a", "redox isomerization followed by self-annulation")
    r["participants"] = [participant("r10a_1", "m_10a", "reactant", [EV[2]]), participant("r10a_2", "m_10a", "reactant", [EV[2]]), participant("cat", cat_mid, "catalyst", [EV[2]], mol_percent=10, ref="m_10a"), participant("solvent", "m_ch2cl2", "solvent", [EV[2]]), participant("p14a", "m_14a", "product", [EV[2]], "main product")]
    r["solvent_conditions"]["composition"] = f("CH2Cl2", ["CH2Cl2"], [EV[2]])
    r["steps"] = [step("reaction", [EV[2]], ["r10a_1", "r10a_2", "cat", "solvent", "p14a"], ("55 °C", q(55, "°C", "55 °C")), operation="redox isomerization followed by reaction with a second equivalent of 10a")]
    if yraw == "–":
        r["outcomes"] = [outcome("table_result", "other", ["m_14a"], [EV[2]], "table dash (–)", "table dash (–)", scope="yield column source marker; not converted to zero or no reaction")]
    elif yraw == "27 (12)^b":
        o = outcome("yield", "yield", ["m_14a"], [EV[2]], yraw, None, "displayed product 14a")
        o["value"] = unreadable(yraw, [EV[2]], "The meaning of the parenthetical 12 and superscript b is not defined on the supplied main-text page; no value is selected.")
        r["outcomes"] = [o]
        r["review"]["unresolved"].append("Figure 2A entry P(Me)2Ph reports '27 (12)^b'; the supplied page does not expose footnote b, so the yield value/basis remains unresolved.")
    else:
        r["outcomes"] = [outcome("yield", "yield", ["m_14a"], [EV[2]], f"{yraw}%", q(int(yraw), "%", f"{yraw}%"), "displayed product 14a")]
    classification_annulation(r, [EV[2]])
    reactions.append(r)

r = reaction("figure2B_deuterium", [EV[2]], "specific_experiment", ["mechanistic_study"], "Figure 2B", "Deuterium-labelling redox isomerization/annulation", "isotopically labelled redox isomerization followed by self-annulation")
r["participants"] = [participant("r10d_1", "m_10d", "reactant", [EV[2]]), participant("r10d_2", "m_10d", "reactant", [EV[2]]), participant("cat", "m_pme2ph", "catalyst", [EV[2]], mol_percent=10, ref="m_10d"), participant("solvent", "m_ch2cl2", "solvent", [EV[2]]), participant("p14d", "m_14d", "product", [EV[2]], "main isotopologue")]
r["solvent_conditions"]["composition"] = f("CH2Cl2", ["CH2Cl2"], [EV[2]])
r["steps"] = [step("reaction", [EV[2]], ["r10d_1", "r10d_2", "cat", "solvent", "p14d"], ("55 °C", q(55, "°C", "55 °C")), operation="deuterium-labelling experiment")]
r["outcomes"] = [outcome("yield", "yield", ["m_14d"], [EV[2]], "25% y", q(25, "%", "25% y"), "main isotopologue")]
classification_annulation(r, [EV[2]])
r["review"]["unresolved"].append("Reaction duration and yield measurement basis are not stated in Figure 2B.")
reactions.append(r)

# Table 1 optimization, one independent row per entry.
table1 = [
    (1, "m_pme2ph", "m_ch2cl2", "CH2Cl2", 25, "2:1", 40, "m_13a"),
    (2, "m_pmeph2", "m_ch2cl2", "CH2Cl2", 25, "2:1", 37, "m_13a"),
    (3, "m_pbu3", "m_ch2cl2", "CH2Cl2", 25, "1:1", 23, "m_13a"),
    (4, "m_pph3", "m_ch2cl2", "CH2Cl2", 25, "1:6", 6, "m_13a"),
    (5, "m_pcy3", "m_ch2cl2", "CH2Cl2", 25, ">1:20", 44, "m_15a"),
    (6, "m_pnme2_3", "m_ch2cl2", "CH2Cl2", 25, ">1:20", 53, "m_15a"),
    (7, "m_pme2ph", "m_ch2cl2", "CH2Cl2", 55, "6:1", 53, "m_13a"),
    (8, "m_pme2ph", "m_dce", "(CH2Cl)2", 55, ">20:1", 39, "m_13a"),
    (9, "m_pme2ph", "m_dce", "(CH2Cl)2", 70, ">20:1", 31, "m_13a"),
    (10, "m_pme2ph", "m_dce", "(CH2Cl)2", 90, ">20:1", 25, "m_13a"),
    (11, "m_pme2ph", "m_chlorobenzene", "ClC6H5", 55, ">20:1", 30, "m_13a"),
    (12, "m_pme2ph", "m_thf", "THF", 55, ">20:1", 12, "m_13a"),
]
for entry, cat_mid, solv_mid, solv_label, temp, ratio, y, ytarget in table1:
    ev = [EV[2], EV[3]] if entry == 7 else [EV[2]]
    r = reaction(f"table1_entry_{entry:02d}", ev, "specific_experiment", ["optimization"], f"Table 1, entry {entry}", "Pyrroline 13a synthesis optimization", "redox isomerization/(3+2) allenoate annulation")
    r["participants"] = [participant("r10a", "m_10a", "reactant", ev), participant("r12a", "m_12a", "reactant", ev), participant("cat", cat_mid, "catalyst", ev, mol_percent=20, ref="m_10a"), participant("p13a", "m_13a", "product", ev, "desired pyrroline"), participant("p15a", "m_15a", "product", ev, "alternate annulation product")]
    solvent(r, solv_mid, solv_label, ev)
    r["steps"] = [step("reaction", ev, [p["participant_id"] for p in r["participants"]], (f"{temp} °C", q(temp, "°C", f"{temp} °C")), ("16 h", q(16, "h", "16 h")), "auto-tandem redox isomerization/(3+2) annulation")]
    r["outcomes"] = [
        outcome("reported_yield", "yield", [ytarget], ev, f"{y}%", q(y, "%", f"{y}%"), "pyrroline 13a" if ytarget == "m_13a" else "pyrrole 15a", "isolated", "isolated yield; rows 5-6 parenthetical values refer to 15a"),
        outcome("product_ratio", "selectivity", ["m_13a", "m_15a"], ev, ratio, ratio, scope="13a:15a ratio"),
        outcome("ratio_basis", "other", ["m_13a", "m_15a"], ev, "determined by 1H NMR analysis", "1H NMR", scope="measurement basis for 13a:15a ratio"),
        outcome("dr_13a", "selectivity", ["m_13a"], ev, "all >20:1 dr", {"ratio": ">20:1", "qualifier": "gt"}, scope="diastereomeric ratio of 13a")
    ]
    if entry == 7:
        r["outcomes"].append(outcome("yield_15a", "yield", ["m_15a"], [EV[3]], "9% of pyrrole 15a", q(9, "%", "9%"), "alternate annulation product", "isolated", "textual by-product yield"))
    r["notes"].append(f("Table heading assigns the main yield column to 13a while footnote c visibly reads 'Isolated yield of 15a'; prose and footnote d distinguish the intended targets. The contradiction is retained for review.",
                          "Internal source-label conflict retained; target assignment follows row notation plus prose, not silent correction.", ev))
    r["review"]["unresolved"].append("Table 1 footnote c visibly names 15a although the header/prose assign ordinary-row yields to 13a; expert review should confirm this apparent publication typo.")
    classification_annulation(r, ev)
    reactions.append(r)

# Table 2 scope. Each item is one displayed substrate/product row.
table2 = [
    ("b", "4:1", 79, 63, True), ("c", "11:1", 75, 69, True), ("d", "7:1", 66, 62, True),
    ("e", ">20:1", 77, 77, True), ("f", "12:1", 64, 59, True), ("g", "13:1", 54, 47, True),
    ("h", ">20:1", 17, 17, False), ("i", "4:1", 84, 60, True), ("j", "6:1", 53, 45, True),
    ("k", "10:1", 60, 39, False), ("l", "5:1", 70, 55, True), ("m", "16:1", 73, 68, True),
    ("n", ">20:1", 56, 56, True), ("o", "5:1", 78, 53, True), ("a", "6:1", 62, 53, True),
    ("p", "6:1", 57, 49, True), ("y", ">1:20", 27, 27, False), ("q", "15:1", 44, 41, True),
    ("r", "8:1", 57, 46, True), ("s", ">20:1", 40, 40, True), ("t", "18:1", 79, 72, False),
    ("u", "16:1", 58, 53, False), ("v", "9:1", 55, 42, False), ("w", ">20:1", 48, 48, False),
    ("x", "18:1", 52, 49, False),
]
for suffix, ratio, combined_yield, pyrroline_yield, footnote_c in table2:
    main_mid, side_mid = f"m_13{suffix}", f"m_15{suffix}"
    if suffix == "a":
        main_mid, side_mid = "m_13a", "m_15a"
    amine_mid = f"m_10{suffix}" if suffix in "bcdefghijk" else "m_10a"
    allene_mid = f"m_12{suffix}" if suffix not in "bcdefghijk" and suffix != "a" else "m_12a"
    loading = 30 if footnote_c else 20
    r = reaction(f"table2_scope_{suffix}", [EV[3]], "specific_experiment", ["substrate_scope"], f"Table 2, product pair 13{suffix}:15{suffix}", "Redox isomerization/(3+2) annulation scope", "redox isomerization/(3+2) allenoate annulation")
    r["general_scheme_ids"] = [f"{DOC}:figure1D_generic_scheme"]
    r["participants"] = [participant("amine", amine_mid, "reactant", [EV[3]]), participant("allene", allene_mid, "reactant", [EV[3]]), participant("cat", "m_pme2ph", "catalyst", [EV[3]], mol_percent=loading, ref=amine_mid), participant("main", main_mid, "product", [EV[3]], "pyrroline 13"), participant("side", side_mid, "product", [EV[3]], "pyrrole 15")]
    solvent(r, "m_ch2cl2", "CH2Cl2", [EV[3]])
    r["steps"] = [step("reaction", [EV[3]], [p["participant_id"] for p in r["participants"]], ("55 °C", q(55, "°C", "55 °C")), ("16 h", q(16, "h", "16 h")), "scope reaction under Table 2 conditions")]
    r["outcomes"] = [
        outcome("combined_yield", "yield", [main_mid, side_mid], [EV[3]], f"{combined_yield}% y", q(combined_yield, "%", f"{combined_yield}% y"), "combined purified 13 and 15", "combined yield of purified 13 and 15 after column chromatography", "combined yield across both product families"),
        outcome("pyrroline_yield", "yield", [main_mid], [EV[3]], f"({pyrroline_yield}% y)", q(pyrroline_yield, "%", f"({pyrroline_yield}% y)"), "pyrroline 13", "yield of 13 shown in parentheses; products purified by column chromatography", "product-specific yield"),
        outcome("product_ratio", "selectivity", [main_mid, side_mid], [EV[3]], ratio, ratio, scope=f"13{suffix}:15{suffix} ratio"),
        outcome("ratio_basis", "other", [main_mid, side_mid], [EV[3]], "ratio of 13:15 by 1H NMR analysis", "1H NMR", scope="measurement basis for product ratio"),
        outcome("dr", "selectivity", [main_mid, side_mid], [EV[3]], "all >20:1 dr", {"ratio": ">20:1", "qualifier": "gt"}, scope="diastereomeric ratio; text assigns the cis-adduct")
    ]
    if suffix == "y":
        r["notes"].append(f("Text states that the aryl ketone substrate gives pyrrole 15y; the 13y:15y ratio is >1:20.", "15y is the predominant observed product; 13y remains only the ratio-defined counterpart.", [EV[3]]))
    classification_annulation(r, [EV[3]])
    reactions.append(r)

# Enantioselective experiment, Scheme 1.
r = reaction("scheme1_enantioselective", [EV[4]], "specific_experiment", ["other"], "Scheme 1", "Enantioselective synthesis of pyrroline 13a", "enantioselective redox isomerization/(3+2) annulation")
r["participants"] = [participant("r10a", "m_10a", "reactant", [EV[4]]), participant("r12a", "m_12a", "reactant", [EV[4]]), participant("cat", "m_D3", "catalyst", [EV[4]], mol_percent=20, ref="m_10a"), participant("solvent", "m_ch2cl2", "solvent", [EV[4]]), participant("p13a", "m_13a_er", "product", [EV[4]], "enantioenriched pyrroline"), participant("p15a", "m_15a", "product", [EV[4]], "alternate annulation product")]
r["solvent_conditions"]["composition"] = f("CH2Cl2", ["CH2Cl2"], [EV[4]])
r["steps"] = [step("reaction", [EV[4]], [p["participant_id"] for p in r["participants"]], ("50 °C", q(50, "°C", "50 °C")), ("16 h", q(16, "h", "16 h")), "enantioselective auto-tandem annulation")]
r["outcomes"] = [
    outcome("yield_13a", "yield", ["m_13a_er"], [EV[4]], "32% y", q(32, "%", "32% y"), "enantioenriched pyrroline", "isolated", "Scheme 1 product 13a"),
    outcome("yield_15a", "yield", ["m_15a"], [EV[4]], "6% y", q(6, "%", "6% y"), "alternate annulation product", "isolated", "Scheme 1 product 15a"),
    outcome("dr", "selectivity", ["m_13a_er"], [EV[4]], ">20:1 dr", {"ratio": ">20:1", "qualifier": "gt"}, scope="diastereomeric ratio"),
    outcome("er", "selectivity", ["m_13a_er"], [EV[4]], "67:33 er", "67:33", scope="enantiomeric ratio; major/minor absolute configurations unassigned"),
    outcome("er_basis", "other", ["m_13a_er"], [EV[4]], "Determined by HPLC with chiral stationary phases", "chiral-stationary-phase HPLC", scope="measurement basis for er")
]
classification_annulation(r, [EV[4]])
reactions.append(r)

# Figure 3A derivatizations.
r = reaction("figure3A_thiolate_addition", [EV[4]], "specific_experiment", ["derivatization"], "Figure 3A, upper arrow", "Thiolate addition to pyrroline 13m", "thiolate conjugate addition")
r["participants"] = [participant("r13m", "m_13m", "reactant", [EV[4]]), participant("thiol", "m_hexanethiol", "reagent", [EV[4]]), participant("base", "m_n_buli", "reagent", [EV[4]]), participant("solvent", "m_thf", "solvent", [EV[4]]), participant("p16m", "m_16m", "product", [EV[4]], "main product")]
r["solvent_conditions"]["composition"] = f("THF", ["THF"], [EV[4]])
r["steps"] = [
    step("low_temperature_addition", [EV[4]], ["r13m", "thiol", "base", "solvent"], ("-78 °C", q(-78, "°C", "-78 °C")), operation="reagents combined/initiated at -78 °C; exact addition order and duration unreported"),
    step("warming_and_reaction", [EV[4]], ["r13m", "thiol", "base", "solvent", "p16m"], ("rt", "rt"), ("20 h", q(20, "h", "20 h")), "warm from -78 °C to rt; 20 h shown for the overall arrow")
]
r["outcomes"] = [outcome("yield", "yield", ["m_16m"], [EV[4]], "63% y", q(63, "%", "63% y"), "main product", basis=None), outcome("dr", "selectivity", ["m_16m"], [EV[4]], "2:1 dr", "2:1", scope="diastereomeric ratio")]
r["review"]["unresolved"].append("The 20 h duration cannot be apportioned between the -78 °C stage, warming ramp and any room-temperature hold; the ramp time is unreported.")
classification_thiol(r, [EV[4]])
reactions.append(r)

r = reaction("figure3A_hydrogenation", [EV[4]], "specific_experiment", ["derivatization"], "Figure 3A, lower arrow", "Hydrogenation of pyrroline 13t", "alkene hydrogenation")
r["participants"] = [participant("r13t", "m_13t", "reactant", [EV[4]]), participant("cat", "m_pdc", "catalyst", [EV[4]]), participant("hydrogen", "m_h2", "reagent", [EV[4]]), participant("solvent", "m_etoh", "solvent", [EV[4]]), participant("p17t", "m_17t", "product", [EV[4]], "main product")]
r["solvent_conditions"]["composition"] = f("EtOH", ["EtOH"], [EV[4]])
r["steps"] = [step("reaction", [EV[4]], [p["participant_id"] for p in r["participants"]], ("rt", "rt"), ("16 h", q(16, "h", "16 h")), "hydrogenation", "H2")]
r["outcomes"] = [outcome("yield", "yield", ["m_17t"], [EV[4]], "72% y", q(72, "%", "72% y"), "main product"), outcome("dr", "selectivity", ["m_17t"], [EV[4]], "2:1 dr", "2:1", scope="diastereomeric ratio")]
r["notes"].append(f("10% Pd/C denotes the catalyst material composition as printed; it is not converted into a reaction mol% loading.", "catalyst composition only", [EV[4]]))
r["review"]["unresolved"].append("Hydrogen pressure and yield measurement basis are not reported in Figure 3A.")
classification_hydrogenation(r, [EV[4]])
reactions.append(r)

# Figure 3B is an author-proposed mechanism, not a performed experiment.
for label in ["11a", "12b", "13b", "18", "19", "20", "21", "22"]:
    ensure(f"m_{label}", label, [EV[4]], name=f"proposed-mechanism species {label}", kind="allene" if label == "12b" else "ordinary")
r = reaction("figure3B_proposed_mechanism", [EV[4]], "proposed_mechanism", [], "Figure 3B", "Proposed mechanism for phosphine-catalyzed synthesis of pyrroline 13b", "author-proposed dual-cycle mechanism")
r["participants"] = [participant(f"p{i}", f"m_{label}", "unknown", [EV[4]]) for i, label in enumerate(["10a", "12b", "18", "19", "11a", "20", "21", "22", "13b"], 1)]
r["notes"] = [f("Cycle A: phosphine 1,4-addition to 10a, elimination to enamine 19, then imine 11a. Cycle B: addition to allenoate 12b, union with 11a, annulation and catalyst elimination to 13b.",
                  "Author-proposed mechanism; not a separate experimental observation.", [EV[4]])]
classification_unresolved(r, [EV[4]], "Proposed elementary mechanism is not classified as a verified experimental net reaction.")
reactions.append(r)

bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": f"{DOC}:figure2C_computation", "document_id": DOC, "kind": "computational_mechanistic_context", "evidence_ids": [EV[2]], "summary": "M06-2X/6-311+G(d,p), SMD dichloromethane, 298.15 K comparison of phosphine addition barriers; not a synthetic experiment."},
    {"context_id": f"{DOC}:page5_references", "document_id": DOC, "kind": "bibliography", "evidence_ids": [EV[5]], "summary": "Page 5 contains references and manuscript dates; no new reaction experiment is depicted."},
    {"context_id": f"{DOC}:xray_13x", "document_id": DOC, "kind": "structure_assignment_context", "evidence_ids": [EV[3], EV[5]], "summary": "Single-crystal X-ray analysis of 13x was used to confirm the stereochemical assignment; CCDC 2086419 is cited, but the crystallographic file is not supplied here."},
]

OUT.write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]), "source_context": len(bundle["source_context"])}, ensure_ascii=False))
