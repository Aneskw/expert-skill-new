#!/usr/bin/env python3
"""Session-local transcription for DOI 10.1021/acs.orglett.1c00997.

This file records machine transcription decisions from the five supplied page
images.  It intentionally leaves records at needs_review and does not modify
the extraction tools, schemas, or quality rules.
"""

from pathlib import Path
import copy
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))

import standard_records as st


DOC = "doc_5db5615b0ea4506c"
P = {n: f"{DOC}:page{n}:image" for n in range(1, 6)}
T = {n: f"{DOC}:page{n}:text" for n in range(1, 6)}
OUT = ROOT / "runs/current"


def extracted(raw, normalized=None, ev=(), status=None):
    if normalized is None:
        normalized = raw
    return st.fact(raw, normalized, ev, method="original_extraction", status=status)


def inherited(raw, normalized=None, ev=(), derivation=None):
    if normalized is None:
        normalized = raw
    return st.fact(raw, normalized, ev, method="inherited_general_condition",
                   derivation=derivation)


def inferred(normalized, ev, derivation, raw=None):
    return st.fact(raw, normalized, ev, method="inference", derivation=derivation)


def derived(raw, normalized, ev, derivation):
    return st.fact(raw, normalized, ev, method="structure_derivation", derivation=derivation)


def missing(ev=(), status="unreported"):
    return st.fact(evidence=ev, status=status)


stereo_diagnostics = {}


def add_molecule(bundle, mid, label, name, ev, smiles=None, *,
                 chirality=None, geometry=None, assignment=None, mixture=None,
                 identity_note=None, aliases=None):
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = extracted(label, label, ev)
    m["names"]["original"] = extracted(name, name, ev)
    m["names"]["english"] = extracted(name, name, ev)
    m["names"]["abbreviations"] = extracted(aliases, aliases, ev) if aliases else missing(ev)
    m["names"]["aliases"] = extracted(aliases, aliases, ev) if aliases else missing(ev)
    m["cas_numbers"] = missing(ev)
    m["external_properties"] = []
    m["identity_note"] = extracted(identity_note, identity_note, ev) if identity_note else missing(ev)

    if smiles:
        diag = st.stereo_inspect(smiles)
        stereo_diagnostics[mid] = diag
        m["representation_type"] = "smiles"
        m["smiles"] = derived(
            f"Structure drawing/identity for {label}", smiles, ev,
            "Candidate SMILES transcribed from the supplied source drawing; RDKit parsing and stereo diagnostics do not establish source correctness.")
        iso = diag["canonical_isomeric_smiles"]
        m["stereochemistry"]["isomeric_smiles"] = derived(
            f"Candidate structure for {label}", iso, ev,
            "Canonical isomeric SMILES calculated from the transcribed candidate graph with standard stereo diagnostics.")
        m["stereochemistry"]["identity_key"] = derived(
            f"Candidate structure for {label}", diag["stereo_graph_key"], ev,
            "Stereo graph key calculated from the candidate isomeric graph; it is not a sample identity claim.")
        potential_atoms = [x for x in diag["potential_stereo"] if x["type"].startswith("Atom_")]
        if chirality is None:
            chirality = "unspecified" if potential_atoms else "achiral"
        if geometry is None:
            geometry = "unspecified" if any(
                x["type"] == "Bond_Double" and x["specified"] == "Unspecified"
                for x in diag["potential_stereo"]
            ) else "not_applicable"
        if chirality == "achiral":
            m["stereochemistry"]["chirality_status"] = inferred(
                "achiral", ev,
                "Candidate graph has no tetrahedral stereocenter and no source-indicated stereogenic axis; terminal allenes were separately assessed for identical terminal substituents.")
        else:
            m["stereochemistry"]["chirality_status"] = extracted(chirality, chirality, ev)
        if geometry in {"specified", "mixture", "partially_specified"}:
            m["stereochemistry"]["geometry_status"] = extracted(geometry, geometry, ev)
        else:
            m["stereochemistry"]["geometry_status"] = inferred(
                geometry, ev, "Geometry status derived from the candidate graph and the source drawing.")
        unspecified = any(x["specified"] == "Unspecified" for x in diag["potential_stereo"])
        source_has_stereo = geometry in {"specified", "mixture", "partially_specified"} or chirality in {"specified", "partially_specified"}
        representation = "encoded" if source_has_stereo and not unspecified else "partial" if source_has_stereo or unspecified else "not_applicable"
        m["stereochemistry"]["representation_status"] = inferred(
            representation, ev,
            "Compared source-indicated stereo with the candidate graph; partial is retained whenever a potential/source stereo element is not fully encoded.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing(ev)
        m["stereochemistry"]["isomeric_smiles"] = missing(ev)
        m["stereochemistry"]["identity_key"] = missing(ev)
        m["stereochemistry"]["chirality_status"] = extracted(
            chirality or "unresolved", chirality or "unresolved", ev) if chirality else inferred(
                "unresolved", ev, "The source identity was retained, but a safe complete graph was not transcribed.")
        m["stereochemistry"]["geometry_status"] = extracted(
            geometry or "unresolved", geometry or "unresolved", ev) if geometry else inferred(
                "unresolved", ev, "The source identity was retained, but a safe complete graph was not transcribed.")
        rep = "not_encoded" if geometry in {"specified", "mixture", "partially_specified"} or chirality in {"specified", "partially_specified"} else "unresolved"
        m["stereochemistry"]["representation_status"] = inferred(
            rep, ev, "Source stereo description retained textually because no safe complete SMILES was asserted.")

    if assignment is not None:
        m["stereochemistry"]["assignment"] = extracted(assignment, assignment, ev)
    else:
        m["stereochemistry"]["assignment"] = missing(ev)
    if geometry == "specified" and assignment:
        m["stereochemistry"]["double_bonds"] = extracted(assignment, assignment, ev)
    else:
        m["stereochemistry"]["double_bonds"] = missing(ev)
    m["stereochemistry"]["atom_centers"] = missing(ev)
    m["stereochemistry"]["axial"] = missing(ev)
    m["stereochemistry"]["mixture"] = extracted(mixture, mixture, ev) if mixture is not None else missing(ev)
    bundle["molecules"].append(m)
    return m


def add_template(bundle, mid, label, name, fixed, variables, ev, note=None):
    m = add_molecule(bundle, mid, label, name, ev, identity_note=note)
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": extracted(fixed, {"description": fixed}, ev),
        "variables": [
            {
                "label": v,
                "attachment_points": extracted(f"Attachment point for {v} shown in {label}", [v], ev),
                "definition": extracted(desc, desc, ev),
                "allowed_values": missing(ev),
                "constraints": missing(ev),
            }
            for v, desc in variables
        ],
        "instance_mapping": [],
        "reason": extracted(
            f"Generic R/Ar notation in {label}",
            "The source uses variable substituents; a concrete SMILES would invent an instance.", ev),
    }
    return m


def participant(pid, mid, role, ev, *, eq=None, molpct=None, ref=None, stage="reaction", product_role=None):
    p = st.new_participant(pid, mid, role, ev, stage=stage)
    if product_role is not None:
        p["product_role"] = extracted(product_role, product_role, ev)
    if eq is not None:
        p["charges"]["equivalents"] = extracted(eq[0], {"number": eq[1], "unit": "equiv"}, ev)
    if molpct is not None:
        p["charges"]["mol_percent"] = extracted(molpct[0], {"number": molpct[1], "unit": "mol%"}, ev)
    if (eq is not None or molpct is not None) and ref:
        p["equivalence_reference"] = extracted(f"relative to {ref}", ref, ev)
    return p


def step(step_id, ev, pids, operation, temp=None, time=None, *, inherited_fields=False):
    s = st.new_step(step_id, ev)
    s["participant_ids"] = list(pids)
    f = inherited if inherited_fields else extracted
    s["operation"] = f(operation, operation, ev)
    s["temperature"] = f(temp[0], {"number": temp[1], "unit": "°C"}, ev) if temp else missing(ev)
    s["time"] = f(time[0], {"number": time[1], "unit": time[2]}, ev) if time else missing(ev)
    s["atmosphere"] = missing(ev)
    s["pressure"] = missing(ev)
    return s


def outcome(oid, kind, targets, ev, value_raw, value_norm, *, role=None,
            yield_type=None, basis=None, scope=None, inherited_value=False):
    o = st.new_outcome(oid, kind, targets, ev)
    vf = inherited if inherited_value else extracted
    o["value"] = vf(value_raw, value_norm, ev)
    o["product_role"] = extracted(role, role, ev) if role else missing(ev)
    o["yield_type"] = extracted(yield_type, yield_type, ev) if yield_type else missing(ev)
    o["measurement_basis"] = extracted(basis, basis, ev) if basis else missing(ev)
    o["scope"] = extracted(scope, scope, ev) if scope else missing(ev)
    return o


def set_classification(r, kind, ev):
    changes = []
    if kind == "digermylation":
        changes = [
            {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "Ge"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "Ge"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 3, "after_order": 2, "change_type": "order_decreased"},
            {"change_id": "b4", "atoms": [{"atom_map": None, "element": "Ge"}, {"atom_map": None, "element": "Ge"}], "before_order": 1, "after_order": 0, "change_type": "broken"},
        ]
    elif kind == "hydrogermylation":
        changes = [
            {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "Ge"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
            {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
            {"change_id": "b4", "atoms": [{"atom_map": None, "element": "Ge"}, {"atom_map": None, "element": "H"}], "before_order": 1, "after_order": 0, "change_type": "broken"},
        ]
    elif kind == "hydrogermylation_branches":
        changes = [
            {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "Ge"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
            {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 1, "after_order": 2, "change_type": "order_increased"},
        ]
    elif kind == "cross_coupling_double":
        changes = [
            {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "Ge"}], "before_order": 1, "after_order": 0, "change_type": "broken"},
            {"change_id": "b4", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "Ge"}], "before_order": 1, "after_order": 0, "change_type": "broken"},
            {"change_id": "b5", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "I"}], "before_order": 1, "after_order": 0, "change_type": "broken"},
            {"change_id": "b6", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "I"}], "before_order": 1, "after_order": 0, "change_type": "broken"},
        ]
    elif kind == "cross_coupling_single":
        changes = [
            {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "Ge"}], "before_order": 1, "after_order": 0, "change_type": "broken"},
            {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "I"}], "before_order": 1, "after_order": 0, "change_type": "broken"},
        ]
    elif kind == "overgermylation":
        changes = [
            {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "Ge"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "Ge"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b3", "atoms": [{"atom_map": None, "element": "O"}, {"atom_map": None, "element": "Ge"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b4", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 3, "after_order": 2, "change_type": "order_decreased"},
            {"change_id": "b5", "atoms": [{"atom_map": None, "element": "O"}, {"atom_map": None, "element": "H"}], "before_order": 1, "after_order": 0, "change_type": "broken"},
            {"change_id": "b6", "atoms": [{"atom_map": None, "element": "Ge"}, {"atom_map": None, "element": "Ge"}], "before_order": 1, "after_order": 0, "change_type": "broken"},
        ]
    elif kind == "not_applicable":
        r["classification"]["mapped_reaction_smiles"] = missing(ev)
        r["classification"]["bond_changes"] = []
        r["classification"]["labels"] = inferred([], ev, "No observed net transformation is classified for this generic/hypothetical record.")
        r["classification"]["coverage"] = inferred("not_applicable", ev, "Proposed mechanisms and context-only records are not assigned observed bond-formation labels.")
        r["classification"]["limitations"] = ["Hypothetical elementary arrows were not converted into experimentally verified bond changes."]
        r["classification"]["model_votes"] = []
        return
    else:
        raise ValueError(kind)

    labels_by_code = {}
    for ch in changes:
        if ch["change_type"] in {"formed", "hydrogen_gain", "order_increased"}:
            label = st.bond_label(ch["atoms"][0]["element"], ch["atoms"][1]["element"], ch["after_order"])
            labels_by_code.setdefault(label["code"], {**label, "supporting_change_ids": []})["supporting_change_ids"].append(ch["change_id"])
    r["classification"]["mapped_reaction_smiles"] = missing(ev)
    r["classification"]["bond_changes"] = [
        derived(None, ch, ev, "Partial net bond change transcribed from the drawn source transformation; no complete atom mapping was asserted.")
        for ch in changes
    ]
    r["classification"]["labels"] = inferred(
        list(labels_by_code.values()), ev,
        "Taxonomy labels calculated from all recorded formed/promoted bonds under bond-formation-20261016.")
    r["classification"]["coverage"] = inferred(
        "partial", ev,
        "Source-supported net changes are recorded without claiming complete atom mapping or coproduct accounting.")
    r["classification"]["limitations"] = [
        "No complete atom mapping was asserted.",
        "Hydrogen changes are net changes; elementary proton/hydride provenance is not claimed.",
        "For product mixtures, labels are the union of source-supported branches and do not imply one product contains every change.",
    ]
    r["classification"]["model_votes"] = []


def base_reaction(bundle, rid, ev, record_type, categories, location, name, transformation,
                  unresolved):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = extracted(categories, categories, ev) if categories is not None else missing(ev)
    r["location"] = extracted(location, location, ev)
    r["reaction_name"] = extracted(name, name, ev) if name else missing(ev)
    r["transformation_type"] = inferred(transformation, ev, "Transformation type normalized from the source title/scheme.") if transformation else missing(ev)
    r["reaction_scale"] = missing(ev)
    r["solvent_conditions"] = {k: missing(ev) for k in ("composition", "mixture_ratio", "total_volume", "concentration", "concentration_basis")}
    r["special_conditions"] = {k: missing(ev) for k in ("light_source", "wavelength", "power", "current", "voltage", "electrodes", "flow_parameters", "reactor")}
    r["procedure"] = {k: missing(ev) for k in ("original_text", "general_procedure_reference", "modifications", "addition_order", "pretreatment")}
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [], "unresolved": list(unresolved)}
    r["completeness"] = "partial"
    bundle["reactions"].append(r)
    return r


dataset = st.read(OUT / "dataset.json")
dataset["documents"][0]["title"] = "Synthesis of Vinylgermanes via the Au/TiO2-Catalyzed cis-1,2-Digermylation of Alkynes and the Regioselective Hydrogermylation of Allenes"
st.write(OUT / "dataset.json", dataset)
bundle = st.empty(dataset["documents"])


# Shared reagents, catalysts, solvents and generic source structures.
add_template(bundle, "m_generic_1", "1", "generic terminal alkyne 1", "R-C≡CH", [("R", "substituent shown on the terminal alkyne series")], [P[1], P[2]])
add_template(bundle, "m_generic_3", "3", "generic cis-1,2-digermylated alkene 3", "R-C(GeMe3)=CH-GeMe3 with the two GeMe3 groups cis", [("R", "substituent inherited from terminal alkyne 1")], [P[1], P[2]], "The paper describes the series as cis/syn products.")
add_template(bundle, "m_generic_4", "4", "generic terminal allene 4", "R1R2C=C=CH2", [("R1", "first substituent"), ("R2", "second substituent; may be H")], [P[1], P[3]])
add_template(bundle, "m_generic_6", "6", "generic major vinylgermane 6", "R1R2CH-C(GeEt3)=CH2", [("R1", "first allene substituent"), ("R2", "second allene substituent; may be H")], [P[1], P[3]])
add_template(bundle, "m_generic_activated_3", "activated cis adduct", "generic cis adduct from ClMe2Ge-GeMe2Cl", "R-C(GeMe2Cl)=CH-GeMe2Cl", [("R", "alkyne substituent")], [P[1]])
add_template(bundle, "m_generic_7", "7", "generic aryl iodide 7", "Ar'-I", [("Ar'", "aryl group introduced in cross coupling")], [P[3]])
add_template(bundle, "m_generic_8", "8-gem/8-cis/8-trans", "generic diarylated products 8", "Products formed from digermyl alkene 3 and two Ar' groups", [("Ar", "aryl group already present in 3"), ("Ar'", "aryl group from iodide 7")], [P[3], P[4]])

simple_molecules = [
    ("m_digermane_2", "2", "hexamethyldigermane", "C[Ge](C)(C)[Ge](C)(C)C", [P[1], P[2]]),
    ("m_activated_digermane", "ClMe2Ge-GeMe2Cl", "1,2-dichlorotetramethyldigermane", "Cl[Ge](C)(C)[Ge](C)(C)Cl", [P[1]]),
    ("m_et3geh_5", "5", "triethylgermane", "[GeH](CC)(CC)CC", [P[1], P[3]]),
    ("m_dce", "DCE", "1,2-dichloroethane", "ClCCCl", [P[2], P[3]]),
    ("m_dmf", "DMF", "N,N-dimethylformamide", "CN(C)C=O", [P[3]]),
    ("m_p_iodotoluene", "p-iodotoluene", "1-iodo-4-methylbenzene", "Cc1ccc(I)cc1", [P[3], P[4]]),
    ("m_p_fluoroiodobenzene", "p-fluoroiodobenzene", "1-fluoro-4-iodobenzene", "Fc1ccc(I)cc1", [P[3], P[4]]),
]
for mid, label, name, smi, ev in simple_molecules:
    add_molecule(bundle, mid, label, name, ev, smi)

for mid, label, name, ev in [
    ("m_au_tio2", "Au/TiO2", "gold nanoparticles supported on titanium dioxide", [P[1], P[2], P[3]]),
    ("m_pd0", "Pd(0)", "palladium(0) catalyst", [P[1]]),
    ("m_pt_acac2", "Pt(acac)2", "platinum(II) acetylacetonate", [P[1], T[2]]),
    ("m_ph3pau_ntf2", "Ph3PAuNTf2", "triphenylphosphine gold(I) bis(trifluoromethanesulfonyl)imide", [P[2]]),
    ("m_pd_ag_system", "Pd2dba3/AgBF4", "Pd2(dba)3/AgBF4 catalytic system", [P[3], P[4]]),
]:
    add_molecule(bundle, mid, label, name, ev, identity_note="Catalyst/support or combined catalytic system is source-identified; no safe single-molecule SMILES was asserted.")


# Table 1 substrates and products. Dense or visually ambiguous entries stay unresolved.
table1 = [
    ("a", "phenylacetylene", "C#Cc1ccccc1", "C[Ge](C)(C)/C(=C\\[Ge](C)(C)C)c1ccccc1", 95),
    ("b", "4-ethynyltoluene", "C#Cc1ccc(C)cc1", "C[Ge](C)(C)/C(=C\\[Ge](C)(C)C)c1ccc(C)cc1", 93),
    ("c", "4-ethynylanisole", "C#Cc1ccc(OC)cc1", "C[Ge](C)(C)/C(=C\\[Ge](C)(C)C)c1ccc(OC)cc1", 96),
    ("d", "substituted arylacetylene 1d", None, None, 96),
    ("e", "substituted arylacetylene 1e", None, None, 95),
    ("f", "4-bromophenylacetylene", "C#Cc1ccc(Br)cc1", "C[Ge](C)(C)/C(=C\\[Ge](C)(C)C)c1ccc(Br)cc1", 90),
    ("g", "ethynylcyclohexene 1g", "C#CC1=CCCCC1", "C[Ge](C)(C)/C(=C\\[Ge](C)(C)C)C1=CCCCC1", 93),
    ("h", "cyclohexyl-substituted terminal alkyne 1h", None, None, 85),
    ("i", "ethynylcyclopropane", "C#CC1CC1", "C[Ge](C)(C)/C(=C\\[Ge](C)(C)C)C1CC1", 84),
    ("j", "1-hexyne", "C#CCCCC", "C[Ge](C)(C)/C(=C\\[Ge](C)(C)C)CCCC", 89),
    ("k", "3-phenyl-1-propyne", "C#CCc1ccccc1", "C[Ge](C)(C)/C(=C\\[Ge](C)(C)C)Cc1ccccc1", 93),
    ("l", "ethynylcyclohexane 1l", "C#CC1CCCCC1", "C[Ge](C)(C)/C(=C\\[Ge](C)(C)C)C1CCCCC1", 94),
    ("m", "propargyl benzoate", "C#CCOC(=O)c1ccccc1", "C[Ge](C)(C)/C(=C\\[Ge](C)(C)C)COC(=O)c1ccccc1", 93),
    ("n", "pent-4-ynenitrile", "N#CCCC#C", "C[Ge](C)(C)/C(=C\\[Ge](C)(C)C)CCC#N", 96),
    ("o", "trimethylsilylacetylene", "C#C[Si](C)(C)C", None, 76),
    ("p", "5-chloropent-1-yne", "C#CCCCCCl", "C[Ge](C)(C)/C(=C\\[Ge](C)(C)C)CCCCl", 94),
    ("q", "propargyl alcohol", "C#CCO", "C[Ge](C)(C)/C(=C\\[Ge](C)(C)C)CO", 91),
]

for suffix, sname, ssmi, psmi, _yield in table1:
    add_molecule(bundle, f"m_1{suffix}", f"1{suffix}", sname, [P[2]], ssmi,
                 chirality="achiral" if ssmi else None,
                 identity_note="Candidate identity follows the Table 1 drawing; ambiguous dense entries remain unresolved.")
    if suffix != "o":
        add_molecule(bundle, f"m_3{suffix}", f"3{suffix}", f"cis-1,2-digermylated product 3{suffix}", [P[2]], psmi,
                     chirality="achiral" if psmi else "unresolved", geometry="specified",
                     assignment={"source_term": "cis-1,2-digermylated", "normalized_assignment": "Z at the product C=C"},
                     identity_note="The paper assigns syn/cis addition by 2D NMR; candidate graph remains a machine transcription.")

add_molecule(bundle, "m_3o_Z", "3o-Z", "Z isomer of digermylated trimethylsilylacetylene", [P[2], P[3]],
             "C[Ge](C)(C)/C(=C\\[Ge](C)(C)C)[Si](C)(C)C", chirality="achiral", geometry="specified",
             assignment={"configuration": "Z", "source": "Table 1 footnote c and page 3 discussion"})
add_molecule(bundle, "m_3o_E", "3o-E", "E isomer of digermylated trimethylsilylacetylene", [P[2], P[3]],
             "C[Ge](C)(C)/C(=C/[Ge](C)(C)C)[Si](C)(C)C", chirality="achiral", geometry="specified",
             assignment={"configuration": "E", "source": "Table 1 footnote c and page 3 discussion"})
add_molecule(bundle, "m_3q_prime", "3q'", "O-trimethylgermyl-protected derivative of 3q", [P[2]],
             "C[Ge](C)(C)/C(=C\\[Ge](C)(C)C)CO[Ge](C)(C)C", chirality="achiral", geometry="specified",
             assignment={"source_term": "cis digermylated alkene; O-trimethylgermyl protected"})


# Table 2 substrates, major products and explicitly ratio-bearing minor regioisomers.
table2 = [
    ("a", "phenylallene", "C=C=Cc1ccccc1", "C=C([Ge](CC)(CC)CC)Cc1ccccc1", 82, 2, (94, 6)),
    ("b", "4-methylphenylallene", "C=C=Cc1ccc(C)cc1", "C=C([Ge](CC)(CC)CC)Cc1ccc(C)cc1", 83, 2, (92, 8)),
    ("c", "4-methoxyphenylallene", "C=C=Cc1ccc(OC)cc1", "C=C([Ge](CC)(CC)CC)Cc1ccc(OC)cc1", 86, 2, (98, 2)),
    ("d", "4-chlorophenylallene", "C=C=Cc1ccc(Cl)cc1", "C=C([Ge](CC)(CC)CC)Cc1ccc(Cl)cc1", 81, 2, (96, 4)),
    ("e", "cyclohexylallene", "C=C=CC1CCCCC1", "C=C([Ge](CC)(CC)CC)CC1CCCCC1", 84, 2, (99, 1)),
    ("f", "methoxyalkyl-substituted allene 4f", None, None, 75, 2, None),
    ("g", "1,1-dimethylallene", "C=C=C(C)C", "C=C([Ge](CC)(CC)CC)C(C)C", 70, 6, None),
    ("h", "ethyl 2-methylbuta-2,3-dienoate", "C=C=C(C)C(=O)OCC", "C=C([Ge](CC)(CC)CC)C(C)C(=O)OCC", 86, 6, None),
    ("i", "diethyl allene-1,1-dicarboxylate 4i", "C=C=C(C(=O)OCC)C(=O)OCC", "C=C([Ge](CC)(CC)CC)C(C(=O)OCC)C(=O)OCC", 88, 6, None),
]
for suffix, sname, ssmi, psmi, _yield, _time, ratio in table2:
    add_molecule(bundle, f"m_4{suffix}", f"4{suffix}", sname, [P[3]], ssmi,
                 chirality="achiral" if suffix not in {"h"} else "achiral",
                 identity_note="Terminal allene axial stereogenicity was assessed; the CH2 terminus prevents an allene stereogenic axis in the transcribed candidate.")
    add_molecule(bundle, f"m_6{suffix}", f"6{suffix}", f"major internal-double-bond hydrogermylation product 6{suffix}", [P[3]], psmi,
                 chirality="unspecified" if suffix == "h" else ("achiral" if psmi else None),
                 geometry="not_applicable" if psmi else None,
                 identity_note="Major vinylgermane drawn in Table 2; source gives no absolute configuration.")
    if ratio:
        add_molecule(bundle, f"m_6{suffix}_minor", f"6{suffix}-minor", f"minor terminal-double-bond addition regioisomer corresponding to 4{suffix}", [P[3]],
                     identity_note="The generic bracketed minor regioisomer is source-indicated, but its exact E/Z identity is not assigned in the main text.")

for mid, label, name in [
    ("m_4j", "4j", "conjugated cyclohexenylallene 4j"),
    ("m_6j", "6j", "anticipated hydrogermylation product 6j"),
    ("m_6j_prime", "6j'", "conjugated major vinylgermane 6j'"),
    ("m_4k", "4k", "vinylidenecyclohexane"),
    ("m_6k", "6k", "anticipated hydrogermylation product 6k"),
    ("m_6k_prime", "6k'", "dehydrogenated vinylgermane 6k'"),
]:
    smi = "C=C=C1CCCCC1" if mid == "m_4k" else None
    add_molecule(bundle, mid, label, name, [P[3]], smi,
                 chirality="achiral" if smi else None,
                 identity_note="Source drawing retained; dense cyclic/conjugated products were not forced into a complete candidate graph.")


# Cross-coupling products and monoarylated intermediates.
add_molecule(bundle, "m_8a", "8a", "tris(4-methylphenyl)ethylene", [P[3]],
             "Cc1ccc(C(=C(c2ccc(C)cc2)c2ccc(C)cc2))cc1", chirality="achiral", geometry="not_applicable")
for series in ("b", "c", "d"):
    add_molecule(bundle, f"m_8{series}_gem", f"8{series}-gem", f"geminal diarylated product 8{series}-gem", [P[3], P[4]],
                 chirality="unresolved", geometry="not_applicable",
                 identity_note="Constitutional isomer label and drawing retained; no safe complete SMILES asserted.")
    for geo in ("cis", "trans"):
        add_molecule(bundle, f"m_8{series}_{geo}", f"8{series}-{geo}", f"{geo} diarylated product 8{series}-{geo}", [P[3], P[4]],
                     chirality="unresolved", geometry="specified",
                     assignment={"source_geometric_label": geo},
                     identity_note="Source cis/trans label retained without converting it automatically to E/Z.")
add_molecule(bundle, "m_9a", "9a", "E-monoarylated germyl alkene 9a", [P[3], P[4], P[5]],
             chirality="unresolved", geometry="specified",
             assignment={"configuration": "E", "assignment_method": "NOESY and HMBC; MM2 comparison in footnote 21"},
             identity_note="E assignment is explicit, but a complete graph was not safely transcribed.")
add_molecule(bundle, "m_9c", "9c", "E-monoarylated germyl alkene 9c", [P[4]],
             chirality="unresolved", geometry="specified",
             assignment={"configuration": "E", "source_statement": "isolated and characterized as a single E-isomer"},
             identity_note="E assignment is explicit, but a complete graph was not safely transcribed.")


# Generic current-work schemes.
r = base_reaction(bundle, "generic_digermylation", [P[1], P[2]], "generic_scheme", ["substrate_scope"],
                  {"pages": [1, 2], "regions": ["page 1 abstract graphic", "page 2 Table 1 header"]},
                  "cis-1,2-digermylation of terminal alkynes", "1,2-digermylation",
                  ["Generic R group is not a concrete experiment.", "Workup, atmosphere and scale are unreported in the supplied main text."])
r["participants"] = [
    participant("generic1", "m_generic_1", "reactant", [P[1], P[2]]),
    participant("digermane", "m_digermane_2", "reactant", [P[1], P[2]], eq=("1.2 equiv", 1.2), ref="m_generic_1"),
    participant("catalyst", "m_au_tio2", "catalyst", [P[1], P[2]], molpct=("1 mol % in Au", 1), ref="m_generic_1"),
    participant("solvent", "m_dce", "solvent", [P[2]]),
    participant("generic3", "m_generic_3", "product", [P[1], P[2]], product_role="major product"),
]
r["solvent_conditions"]["composition"] = extracted("dry 1,2-dichloroethane", ["1,2-dichloroethane"], [P[2]])
r["steps"] = [step("reaction", [P[2]], [p["participant_id"] for p in r["participants"]], "react in dry DCE", ("70 °C", 70), ("2 h", 2, "h"))]
r["procedure"]["original_text"] = extracted("terminal alkynes with Me3GeGeMe3 (1.2 equiv), Au/TiO2 (1 mol % in Au), dry DCE, 70 °C, 2 h", "terminal alkynes with Me3GeGeMe3 (1.2 equiv), Au/TiO2 (1 mol % in Au), dry DCE, 70 °C, 2 h", [P[2]])
set_classification(r, "digermylation", [P[1], P[2]])

r = base_reaction(bundle, "generic_hydrogermylation", [P[1], P[2], P[3]], "generic_scheme", ["substrate_scope"],
                  {"pages": [1, 2, 3], "regions": ["abstract graphic", "Table 2 header"]},
                  "regioselective hydrogermylation of terminal allenes", "allene hydrogermylation",
                  ["Generic R groups are not a concrete experiment.", "Workup, atmosphere, scale and yield basis are unreported."])
r["participants"] = [
    participant("generic4", "m_generic_4", "reactant", [P[1], P[3]]),
    participant("et3geh", "m_et3geh_5", "reactant", [P[1], P[3]], eq=("1.5 equiv", 1.5), ref="m_generic_4"),
    participant("catalyst", "m_au_tio2", "catalyst", [P[1], P[3]], molpct=("1 mol %", 1), ref="m_generic_4"),
    participant("solvent", "m_dce", "solvent", [P[2], P[3]]),
    participant("generic6", "m_generic_6", "product", [P[1], P[3]], product_role="major or only product"),
]
r["solvent_conditions"]["composition"] = extracted("DCE", ["1,2-dichloroethane"], [P[3]])
r["steps"] = [step("reaction", [P[3]], [p["participant_id"] for p in r["participants"]], "hydrogermylation", ("70 °C", 70), ("2–6 h", 2, "h"))]
r["steps"][0]["time"] = extracted("2–6 h", {"minimum": 2, "maximum": 6, "unit": "h", "raw_text": "2–6 h"}, [P[3]])
set_classification(r, "hydrogermylation", [P[1], P[3]])


# Page 1 literature examples are retained as literature examples, never current experiments.
lit_specs = [
    ("literature_pd_activated", "m_activated_digermane", "m_generic_activated_3", "m_pd0", "moderate to good yields", "moderate to good yield", None, [P[1]]),
    ("literature_pd_unactivated", "m_digermane_2", "m_generic_3", "m_pd0", "trace amounts of products", "trace", None, [P[1]]),
    ("literature_pt_unactivated", "m_digermane_2", "m_generic_3", "m_pt_acac2", "up to 60% conversion yield", {"number": 60, "unit": "%", "qualifier": "le", "raw_text": "up to 60%"}, 5, [P[1], P[2]]),
]
for rid, reagent_mid, product_mid, cat_mid, raw_out, norm_out, alkyne_eq, ev in lit_specs:
    r = base_reaction(bundle, rid, ev, "literature_example", None,
                      {"pages": sorted({1, 2} if rid == "literature_pt_unactivated" else {1}), "region": "Scheme 1 and adjacent prior-work prose"},
                      "prior catalytic 1,2-digermylation of alkynes", "1,2-digermylation",
                      ["This is cited prior work, not a new experiment in the paper.", "Exact substrate, scale, workup and measurement basis are not supplied here."])
    r["participants"] = [
        participant("alkyne", "m_generic_1", "reactant", ev, eq=(("5-fold excess", 5) if alkyne_eq else None), ref=reagent_mid if alkyne_eq else None),
        participant("digermane", reagent_mid, "reactant", ev),
        participant("catalyst", cat_mid, "catalyst", ev, molpct=(("5 mol %", 5) if cat_mid == "m_pt_acac2" else None), ref=reagent_mid if cat_mid == "m_pt_acac2" else None),
        participant("product", product_mid, "product", ev, product_role="reported cis adduct"),
    ]
    r["steps"] = [step("reaction", ev, [p["participant_id"] for p in r["participants"]], "literature digermylation", ("120 °C", 120), None)]
    if rid == "literature_pd_unactivated":
        r["outcomes"] = [outcome("trace", "trace", [product_mid], ev, raw_out, norm_out, role="trace product")]
    elif rid == "literature_pt_unactivated":
        r["outcomes"] = [outcome("conversion", "conversion", [reagent_mid], ev, raw_out, norm_out, scope="digermylation conversion yield as described in the citing paper")]
    else:
        r["outcomes"] = [outcome("yield", "yield", [product_mid], ev, raw_out, norm_out, role="cis adduct")]
    set_classification(r, "digermylation", ev)


# Seventeen Table 1 entries.
for suffix, _sname, _ssmi, _psmi, y in table1:
    rid = f"table1_3{suffix}"
    ev = [P[2]] + ([P[3]] if suffix == "o" else [])
    r = base_reaction(bundle, rid, ev, "specific_experiment", ["substrate_scope"],
                      {"page": 2, "region": f"Table 1 entry 3{suffix}"},
                      "cis-1,2-digermylation of terminal alkyne", "1,2-digermylation",
                      ["Reaction scale, atmosphere, workup and purification are unreported in the supplied main text."])
    product_mids = ["m_3o_Z", "m_3o_E"] if suffix == "o" else [f"m_3{suffix}"]
    dig_eq = 2.0 if suffix == "f" else 1.2
    hours = 4 if suffix == "f" else 2
    r["participants"] = [
        participant("substrate", f"m_1{suffix}", "reactant", ev),
        participant("digermane", "m_digermane_2", "reactant", ev, eq=(f"{dig_eq:g} equiv", dig_eq), ref=f"m_1{suffix}"),
        participant("catalyst", "m_au_tio2", "catalyst", ev, molpct=("1 mol % in Au", 1), ref=f"m_1{suffix}"),
        participant("solvent", "m_dce", "solvent", ev),
    ] + [participant(f"product_{i+1}", mid, "product", ev, product_role="product" if suffix == "o" else "major product") for i, mid in enumerate(product_mids)]
    r["solvent_conditions"]["composition"] = inherited("dry 1,2-dichloroethane", ["1,2-dichloroethane"], [P[2]], "Inherited from the Table 1 header/prose for all entries.")
    r["steps"] = [step("reaction", [P[2]], [p["participant_id"] for p in r["participants"]], "digermylation under Table 1 conditions", ("70 °C", 70), (f"{hours} h", hours, "h"), inherited_fields=True)]
    r["procedure"]["general_procedure_reference"] = inherited("Table 1 general conditions", "Table 1 general conditions", [P[2]], "Applied to this table entry; footnote overrides are structured separately.")
    r["outcomes"] = [outcome("isolated_yield", "yield", product_mids, [P[2]], f"{y}%", {"number": y, "unit": "%"}, role="product mixture" if suffix == "o" else "major product", yield_type="isolated yield", basis="isolated")]
    if suffix == "o":
        r["outcomes"].append(outcome("ez_ratio", "selectivity", product_mids, [P[2], P[3]], "Z/E = 75/25", {"ratio": [75, 25], "order": ["Z", "E"]}, basis=None, scope="geometric-isomer ratio Z/E"))
        r["review"]["unresolved"].append("The method used to determine the Z/E ratio is not reported in the supplied main text.")
    else:
        r["outcomes"].append(outcome("cis_selectivity", "selectivity", product_mids, [P[2]], "1,2-cis-digermylated adduct", "cis", basis="2D NMR experiments", scope="syn/cis addition assignment", inherited_value=True))
    set_classification(r, "digermylation", ev)


# Footnote/prose condition variant for propargyl alcohol and explicit Au(I) low-yield control.
r = base_reaction(bundle, "propargyl_overgermylation_3q_prime", [P[2]], "specific_experiment", ["substrate_scope"],
                  {"page": 2, "region": "Table 1 footnote d and adjacent prose"},
                  "digermylation and O-germylation of propargyl alcohol", "digermylation with O-germylation",
                  ["Exact time, product-isolation basis, workup and purification are unreported.", "The prose gives 2.0–2.5 equiv while the table footnote specifies 2.5 equiv; both statements are preserved in notes."])
r["participants"] = [
    participant("substrate", "m_1q", "reactant", [P[2]]),
    participant("digermane", "m_digermane_2", "reactant", [P[2]], eq=("2.5 equiv", 2.5), ref="m_1q"),
    participant("catalyst", "m_au_tio2", "catalyst", [P[2]], molpct=("1 mol % in Au", 1), ref="m_1q"),
    participant("solvent", "m_dce", "solvent", [P[2]]),
    participant("product_3q", "m_3q", "product", [P[2]], product_role="unprotected product"),
    participant("product_3q_prime", "m_3q_prime", "product", [P[2]], product_role="O-trimethylgermyl-protected product"),
]
r["solvent_conditions"]["composition"] = inherited("dry 1,2-dichloroethane", ["1,2-dichloroethane"], [P[2]], "Inherited from Table 1 conditions.")
r["steps"] = [step("reaction", [P[2]], [p["participant_id"] for p in r["participants"]], "digermylation with competing O-germylation", ("70 °C", 70), None, inherited_fields=True)]
r["outcomes"] = [
    outcome("product_ratio", "selectivity", ["m_3q", "m_3q_prime"], [P[2]], "mixture ... in ~1/1 ratio", {"ratio": [1, 1], "qualifier": "approx", "order": ["3q", "3q'"]}, scope="3q/3q' product ratio"),
    outcome("protected_product_yield", "yield", ["m_3q_prime"], [P[2]], "~50% yield", {"number": 50, "unit": "%", "qualifier": "approx", "raw_text": "~50% yield"}, role="O-protected product"),
]
r["notes"] = [extracted("prose: with 2.0–2.5 equiv; footnote: with 2.5 equiv", {"prose": "2.0–2.5 equiv", "table_footnote": "2.5 equiv"}, [P[2]])]
set_classification(r, "overgermylation", [P[2]])

r = base_reaction(bundle, "auI_ph3pau_ntf2_control", [P[2]], "specific_experiment", ["optimization", "control"],
                  {"page": 2, "region": "right-column Au(I) catalyst discussion"},
                  "Au(I)-catalyzed attempted digermylation", "1,2-digermylation",
                  ["The exact alkyne substrate, scale, solvent, workup and yield basis are not reported in the supplied main text."])
r["participants"] = [
    participant("alkyne", "m_generic_1", "reactant", [P[2]]),
    participant("digermane", "m_digermane_2", "reactant", [P[2]]),
    participant("catalyst", "m_ph3pau_ntf2", "catalyst", [P[2]], molpct=("5 mol %", 5), ref="m_generic_1"),
    participant("product", "m_generic_3", "product", [P[2]], product_role="digermylated product"),
]
r["steps"] = [step("reaction", [P[2]], [p["participant_id"] for p in r["participants"]], "Au(I)-catalyzed digermylation attempt", ("70 °C", 70), ("24 h", 24, "h"))]
r["outcomes"] = [outcome("yield", "yield", ["m_generic_3"], [P[2]], "5% yield", {"number": 5, "unit": "%"}, role="digermylated product")]
set_classification(r, "digermylation", [P[2]])


# Eleven Table 2 entries.
for suffix, _sname, _ssmi, _psmi, y, hours, ratio in table2:
    r = base_reaction(bundle, f"table2_6{suffix}", [P[3]], "specific_experiment", ["substrate_scope"],
                      {"page": 3, "region": f"Table 2 entry 6{suffix}"},
                      "regioselective hydrogermylation of terminal allene", "allene hydrogermylation",
                      ["Reaction scale, atmosphere, workup, purification and yield basis are unreported."])
    product_mids = [f"m_6{suffix}"] + ([f"m_6{suffix}_minor"] if ratio else [])
    r["participants"] = [
        participant("allene", f"m_4{suffix}", "reactant", [P[3]]),
        participant("et3geh", "m_et3geh_5", "reactant", [P[3]], eq=("1.5 equiv", 1.5), ref=f"m_4{suffix}"),
        participant("catalyst", "m_au_tio2", "catalyst", [P[3]], molpct=("1 mol %", 1), ref=f"m_4{suffix}"),
        participant("solvent", "m_dce", "solvent", [P[3]]),
    ] + [participant(f"product_{i+1}", mid, "product", [P[3]], product_role="major regioisomer" if i == 0 else "minor regioisomer") for i, mid in enumerate(product_mids)]
    r["solvent_conditions"]["composition"] = inherited("DCE", ["1,2-dichloroethane"], [P[3]], "Inherited from the Table 2 header.")
    r["steps"] = [step("reaction", [P[3]], [p["participant_id"] for p in r["participants"]], "hydrogermylation under Table 2 conditions", ("70 °C", 70), (f"{hours} h", hours, "h"), inherited_fields=True)]
    r["procedure"]["general_procedure_reference"] = inherited("Table 2 general conditions and footnote a", "Table 2 general conditions", [P[3]], "Entry time assigned from the source class: monosubstituted 2 h; 1,1-disubstituted 6 h.")
    r["outcomes"] = [outcome("yield", "yield", [f"m_6{suffix}"], [P[3]], f"{y}%", {"number": y, "unit": "%"}, role="major product")]
    if ratio:
        r["outcomes"].append(outcome("regioselectivity", "selectivity", product_mids, [P[3]], f"internal/terminal double bond addition: {ratio[0]}/{ratio[1]}", {"ratio": list(ratio), "order": ["internal-double-bond addition", "terminal-double-bond addition"]}, scope="regioisomeric product ratio"))
        r["review"]["unresolved"].append("The selectivity measurement method is not reported.")
    elif suffix in {"g", "h", "i"}:
        r["outcomes"].append(outcome("regioselectivity", "selectivity", [f"m_6{suffix}"], [P[2], P[3]], ">99% selectivity", {"number": 99, "unit": "%", "qualifier": "gt", "raw_text": ">99%"}, scope="addition at the more substituted/internal double bond"))
        r["review"]["unresolved"].append("The exact minor product identity and selectivity measurement method are not reported.")
    set_classification(r, "hydrogermylation", [P[3]])

for suffix, y, ratio, hours, branch_note in [
    ("j", 81, (20, 80), 2, "6j is minor; conjugated 6j' is major and is not formed by isomerization of 6j under the reaction conditions."),
    ("k", 79, (66, 34), 6, "6k' is the dehydrogenated side product."),
]:
    products = [f"m_6{suffix}", f"m_6{suffix}_prime"]
    r = base_reaction(bundle, f"table2_6{suffix}", [P[2], P[3]], "specific_experiment", ["substrate_scope"],
                      {"pages": [2, 3], "region": f"Table 2 special entry 4{suffix}"},
                      "regioselective hydrogermylation of terminal allene", "branched allene hydrogermylation",
                      ["Reaction scale, atmosphere, workup, purification, yield basis and ratio method are unreported."])
    r["participants"] = [
        participant("allene", f"m_4{suffix}", "reactant", [P[3]]),
        participant("et3geh", "m_et3geh_5", "reactant", [P[3]], eq=("1.5 equiv", 1.5), ref=f"m_4{suffix}"),
        participant("catalyst", "m_au_tio2", "catalyst", [P[3]], molpct=("1 mol %", 1), ref=f"m_4{suffix}"),
        participant("solvent", "m_dce", "solvent", [P[2], P[3]]),
        participant("product1", products[0], "product", [P[3]], product_role="anticipated product"),
        participant("product2", products[1], "product", [P[3]], product_role="major conjugated product" if suffix == "j" else "dehydrogenated side product"),
    ]
    r["solvent_conditions"]["composition"] = inherited("DCE", ["1,2-dichloroethane"], [P[2], P[3]], "Inherited from the optimized/Table 2 conditions.")
    r["steps"] = [step("reaction", [P[3]], [p["participant_id"] for p in r["participants"]], "hydrogermylation under Table 2 conditions", ("70 °C", 70), (f"{hours} h", hours, "h"), inherited_fields=True)]
    r["outcomes"] = [
        outcome("combined_yield", "yield", products, [P[3]], f"{y}%", {"number": y, "unit": "%"}, role="combined product mixture"),
        outcome("product_ratio", "selectivity", products, [P[3]], f"6{suffix}/6{suffix}' = {ratio[0]}/{ratio[1]}", {"ratio": list(ratio), "order": [f"6{suffix}", f"6{suffix}'"]}, scope="product distribution"),
    ]
    r["notes"] = [extracted(branch_note, branch_note, [P[2], P[3]])]
    set_classification(r, "hydrogermylation_branches", [P[3]])


# Scheme 2 proposed mechanisms: separated by the four depicted pathways, not experiments.
mechanisms = [
    ("mechanism_scheme2_digermylation", "proposed Au-cluster mechanism for alkyne digermylation", "Scheme 2 top: Ge-Ge insertion on Au cluster, alkyne coordination/addition through I, product release."),
    ("mechanism_scheme2_hydrogermylation", "proposed Au-cluster mechanism for allene hydrogermylation", "Scheme 2 middle: eta-1-like II delivers hydride to the more electrophilic substituted carbon, then releases 6."),
    ("mechanism_scheme2_4j", "proposed charge-delocalized pathway from 4j to 6j'", "Scheme 2 lower middle: resonance structures IV/V rationalize hydride transfer and conjugated product 6j'."),
    ("mechanism_scheme2_4k", "proposed dehydrogenative pathway from 4k to 6k'", "Scheme 2 bottom: VI undergoes intramolecular H2 elimination through a proposed six-membered transition state."),
]
for rid, name, note in mechanisms:
    r = base_reaction(bundle, rid, [P[3]], "proposed_mechanism", ["mechanistic_study"],
                      {"page": 3, "region": "Scheme 2"}, name, "proposed mechanism",
                      ["This is an author-proposed mechanism, not a separately observed experiment."])
    r["notes"] = [extracted(note, note, [P[3]])]
    set_classification(r, "not_applicable", [P[3]])


# Scheme 3 four concrete cross-coupling combinations.
cross_specs = [
    ("a", "m_3b", "m_p_iodotoluene", ["m_8a"], 72, None),
    ("b", "m_3a", "m_p_iodotoluene", ["m_8b_gem", "m_8b_cis", "m_8b_trans"], 75, (55, 11, 34)),
    ("c", "m_3c", "m_p_fluoroiodobenzene", ["m_8c_gem", "m_8c_cis", "m_8c_trans"], 72, (53, 23, 24)),
    ("d", "m_3a", "m_p_fluoroiodobenzene", ["m_8d_gem", "m_8d_cis", "m_8d_trans"], 73, (54, 15, 31)),
]
for suffix, substrate_mid, iodide_mid, products, y, ratio in cross_specs:
    r = base_reaction(bundle, f"scheme3_8{suffix}", [P[3]], "specific_experiment", ["substrate_scope", "derivatization"],
                      {"page": 3, "region": f"Scheme 3 row 8{suffix}"},
                      "Pd-catalyzed cross coupling of cis-1,2-digermylated adduct with aryl iodide", "double aryl cross coupling",
                      ["Reaction time, temperature, scale, atmosphere, workup, purification and yield basis are unreported.", "The combined notation Pd2dba3/AgBF4 (5 mol%) does not apportion a separate AgBF4 loading."])
    r["participants"] = [
        participant("substrate", substrate_mid, "reactant", [P[3]]),
        participant("aryl_iodide", iodide_mid, "reactant", [P[3]], eq=("3 equiv", 3), ref=substrate_mid),
        participant("catalytic_system", "m_pd_ag_system", "catalyst", [P[3]], molpct=("5 mol %", 5), ref=substrate_mid),
        participant("solvent", "m_dmf", "solvent", [P[3]]),
    ] + [participant(f"product_{i+1}", mid, "product", [P[3]], product_role="major product" if i == 0 else "minor isomer") for i, mid in enumerate(products)]
    r["solvent_conditions"]["composition"] = extracted("DMF", ["DMF"], [P[3]])
    r["steps"] = [step("reaction", [P[3]], [p["participant_id"] for p in r["participants"]], "Pd2dba3/AgBF4-catalyzed coupling in DMF")]
    r["outcomes"] = [outcome("yield", "yield", products, [P[3]], f"{y}% yield", {"number": y, "unit": "%"}, role="product or product mixture")]
    if ratio:
        r["outcomes"].append(outcome("isomer_ratio", "selectivity", products, [P[3]], f"8{suffix}-gem/8{suffix}-cis/8{suffix}-trans = {ratio[0]}/{ratio[1]}/{ratio[2]}", {"ratio": list(ratio), "order": [f"8{suffix}-gem", f"8{suffix}-cis", f"8{suffix}-trans"]}, scope="constitutional/geometric product distribution"))
    set_classification(r, "cross_coupling_double", [P[3]])


# Monoarylated intermediates and the second coupling of isolated 9c.
for suffix, substrate_mid, iodide_mid, product_mid, ev in [
    ("a", "m_3b", "m_p_iodotoluene", "m_9a", [P[3], P[4], P[5]]),
    ("c", "m_3c", "m_p_fluoroiodobenzene", "m_9c", [P[3], P[4]]),
]:
    core_ev = [P[3], P[4]]
    r = base_reaction(bundle, f"monoarylation_9{suffix}", ev, "specific_experiment", ["mechanistic_study", "derivatization"],
                      {"pages": sorted({3, 4, 5} if suffix == "a" else {3, 4}), "region": f"prose around Schemes 3–4, intermediate 9{suffix}"},
                      "Pd-catalyzed monoarylation of cis-1,2-digermylated alkene", "single aryl cross coupling",
                      ["Exact time, temperature, scale, workup, purification and isolated yield are unreported.", "The stated 25–30% is a relative yield in the crude reaction mixture, not an isolated-yield value."])
    r["participants"] = [
        participant("substrate", substrate_mid, "reactant", core_ev),
        participant("aryl_iodide", iodide_mid, "reactant", core_ev, eq=("1.5 equiv", 1.5), ref=substrate_mid),
        participant("catalytic_system", "m_pd_ag_system", "catalyst", core_ev),
        participant("product", product_mid, "product", core_ev, product_role="single E monoarylated intermediate"),
    ]
    r["steps"] = [step("reaction", core_ev, [p["participant_id"] for p in r["participants"]], "monoarylation under Pd2dba3/AgBF4 conditions")]
    r["outcomes"] = [
        outcome("relative_crude_yield", "yield", [product_mid], core_ev, "approximately 25–30% relative yield in the crude reaction mixture", {"minimum": 25, "maximum": 30, "unit": "%", "qualifier": "approx", "raw_text": "approximately 25–30%"}, role="monoarylated intermediate", yield_type="relative crude yield", basis="relative yield in the crude reaction mixture"),
        outcome("geometry", "selectivity", [product_mid], core_ev, "single E-isomer", "E only", basis="NOESY and HMBC" if suffix == "a" else None, scope="alkene geometry"),
    ]
    if suffix == "a":
        r["notes"] = [extracted("Simple MM2 calculations show that 9a-E is more stable relative to its Z-isomer by ~6 kcal/mol.", {"calculation": "MM2", "relative_stability": "~6 kcal/mol", "more_stable": "9a-E"}, [P[5]])]
    set_classification(r, "cross_coupling_single", ev)

r = base_reaction(bundle, "second_coupling_9c_to_8c", [P[3], P[4]], "specific_experiment", ["mechanistic_study", "derivatization"],
                  {"pages": [3, 4], "region": "page 4 prose below Scheme 4 with Scheme 3 ratio reference"},
                  "second Pd-catalyzed arylation of isolated 9c", "single aryl cross coupling",
                  ["Reaction time, temperature, scale, solvent, workup, purification, total yield and ratio method are unreported."])
r["participants"] = [
    participant("substrate", "m_9c", "reactant", [P[4]]),
    participant("aryl_iodide", "m_p_fluoroiodobenzene", "reactant", [P[4]], eq=("1.5 equiv", 1.5), ref="m_9c"),
    participant("catalytic_system", "m_pd_ag_system", "catalyst", [P[4]]),
    participant("product_gem", "m_8c_gem", "product", [P[3], P[4]], product_role="major product"),
    participant("product_cis", "m_8c_cis", "product", [P[3], P[4]], product_role="minor geometric isomer"),
    participant("product_trans", "m_8c_trans", "product", [P[3], P[4]], product_role="minor geometric isomer"),
]
r["steps"] = [step("reaction", [P[4]], [p["participant_id"] for p in r["participants"]], "treat 9c with Pd2dba3/AgBF4 and p-fluoroiodobenzene")]
r["outcomes"] = [outcome("product_ratio", "selectivity", ["m_8c_gem", "m_8c_cis", "m_8c_trans"], [P[3], P[4]], "ratio identical to Scheme 3: 53/23/24", {"ratio": [53, 23, 24], "order": ["8c-gem", "8c-cis", "8c-trans"]}, scope="product distribution")]
set_classification(r, "cross_coupling_single", [P[3], P[4]])

r = base_reaction(bundle, "mechanism_scheme4_diarylation", [P[4]], "proposed_mechanism", ["mechanistic_study"],
                  {"page": 4, "region": "Scheme 4"},
                  "proposed mechanism for Pd-catalyzed diarylation of cis-1,2-digermylated alkenes", "proposed mechanism",
                  ["This is an author-proposed mechanism, not a separately observed experiment."])
r["notes"] = [extracted("Scheme 4 proposes electrophilic aryl addition, fluoride capture, E-monoarylated 9, competing second additions through VIII/IX, H shift, and Me3GeF elimination.", "Author-proposed branching pathway to 8-gem and 8-cis/8-trans", [P[4]])]
set_classification(r, "not_applicable", [P[4]])


# Context that must not become synthetic experiments.
bundle["source_context"] = [
    {"context_id": f"{DOC}:context:page1_abstract", "document_id": DOC, "evidence_ids": [P[1]], "kind": "graphical_abstract", "note": "Overview graphics repeat the two current generic transformations and reported yield ranges; they are not duplicate concrete experiments."},
    {"context_id": f"{DOC}:context:page2_optimization", "document_id": DOC, "evidence_ids": [P[2]], "kind": "discussion", "note": "Au/Al2O3 and Au/ZnO were less efficient; 0.1 mol% Au/TiO2 gave about 40–70% lower yields in 2 h. No exact substrate-row results are supplied in the main text."},
    {"context_id": f"{DOC}:context:page2_internal_alkynes", "document_id": DOC, "evidence_ids": [P[2]], "kind": "negative_scope_boundary", "note": "Internal alkynes are reported completely unreactive, but no individual substrate/experiment is enumerated; no fabricated product or formal specific record was created."},
    {"context_id": f"{DOC}:context:page2_stereo", "document_id": DOC, "evidence_ids": [P[2]], "kind": "stereochemical_assignment", "note": "The syn-1,2-addition of digermane was confirmed by 2D NMR experiments."},
    {"context_id": f"{DOC}:context:page3_mechanistic_discussion", "document_id": DOC, "evidence_ids": [P[3]], "kind": "mechanistic_discussion", "note": "Text discusses proposed intermediates I–VI and branch selectivity; these statements are not independent synthesis experiments."},
    {"context_id": f"{DOC}:context:page4_supporting_information_notice", "document_id": DOC, "evidence_ids": [P[4]], "kind": "supplement_notice", "note": "The article announces an SI PDF containing spectra. No SI file was supplied or explicitly associated in this workspace."},
    {"context_id": f"{DOC}:context:page5_references", "document_id": DOC, "evidence_ids": [P[5]], "kind": "bibliography_and_footnotes", "note": "Page 5 is bibliography plus footnotes 21–22. Footnote 21's MM2 comparison is linked to the 9a record; references are not experiments."},
]


# Persist the filled machine bundle and raw diagnostics/inventory.
st.write(OUT / "submitted.json", bundle)
st.write(OUT / "raw_model/stereo_diagnostics.json", stereo_diagnostics)
st.write(OUT / "raw_model/source_inventory.json", {
    "document_id": DOC,
    "scope": "original PDF pages 1-5",
    "full_page_images_opened": [f"assets/{DOC}/page{n}.jpg" for n in range(1, 6)],
    "source_lineage_crops": [
        "page1_abstract_scheme.jpg", "page1_scheme1.jpg", "page2_table1.jpg",
        "page3_table2.jpg", "page3_scheme2.jpg", "page3_scheme3.jpg", "page4_scheme4.jpg"
    ],
    "inventory": {
        "page1": ["title/abstract", "graphical abstract", "Scheme 1 prior work/current generic route"],
        "page2": ["Table 1 rows 3a-3q", "footnotes b-d", "Au(I) control prose", "hydrogermylation setup prose"],
        "page3": ["Table 2 rows 6a-6k", "Scheme 2 proposed mechanisms", "Scheme 3 four coupling rows", "monoarylation discussion starts"],
        "page4": ["Scheme 4 proposed mechanism", "9a/9c/second coupling prose", "SI notice", "bibliography begins"],
        "page5": ["bibliography", "footnote 21 MM2 comparison", "footnote 22 reference"]
    },
    "si_association": "unestablished",
    "presentation_conversion": False,
    "chemical_accuracy_verified": False,
})

# Mark mechanical page tasks completed only after every supplied page image was opened.
tasks = st.read(OUT / "tasks.json")
for task in tasks:
    task["state"] = "completed"
    task["reason"] = "Full original page image inspected; findings mapped in coverage/pages_0001_0005.json."
st.write(OUT / "tasks.json", tasks)

print(json.dumps({
    "molecules": len(bundle["molecules"]),
    "reactions": len(bundle["reactions"]),
    "source_context": len(bundle["source_context"]),
    "stereo_diagnostics": len(stereo_diagnostics),
}, ensure_ascii=False, indent=2))
