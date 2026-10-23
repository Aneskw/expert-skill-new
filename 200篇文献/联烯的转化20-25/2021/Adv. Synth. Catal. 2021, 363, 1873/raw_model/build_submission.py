#!/usr/bin/env python3
"""Build the source-derived v3.1 machine bundle for adsc.202001513.

This is an extraction artifact, not a chemistry-validation program.  The source
images were reviewed before this file was authored.  It deliberately preserves
unreported values and source conflicts.
"""
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DOC = "doc_ad9cc345885d5e84"
P = {i: f"{DOC}:page{i}:image" for i in range(1, 6)}
T = {i: f"{DOC}:page{i}:text" for i in range(1, 6)}


def f(raw=None, norm=None, ev=(), method="original_extraction", status=None,
      derivation=None):
    return st.fact(raw, norm, ev, method=method, status=status,
                   derivation=derivation)


def missing(ev=()):
    return f(ev=ev)


def na(ev=()):
    return f(ev=ev, method="not_applicable", status="not_applicable")


def q(raw, number, unit, ev):
    return f(raw, {"number": number, "unit": unit}, ev)


def conflict(raw, alternatives, ev):
    x = f(raw, None, ev, status="conflict")
    x["alternatives"] = alternatives
    return x


def inferred(raw, norm, ev, why):
    return f(raw, norm, ev, method="inference", derivation=why)


bundle = json.loads((RUN / "submitted.json").read_text())
bundle["molecules"] = []
bundle["reactions"] = []
bundle["source_context"] = []
stereo_checks = {}
molecules = {}


def add_molecule(mid, label, ev, *, name=None, smiles=None,
                 chirality="unresolved", geometry="unresolved",
                 representation="not_encoded", stereo_note=None,
                 axial=None, double_bonds=None, assignment=None,
                 reason=None):
    if mid in molecules:
        return molecules[mid]
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    source_name = name or label
    m["names"]["original"] = f(source_name, source_name, ev)
    m["names"]["english"] = f(source_name, source_name, ev)
    m["names"]["abbreviations"] = missing(ev)
    m["names"]["aliases"] = missing(ev)
    m["cas_numbers"] = missing(ev)
    if smiles:
        diag = st.stereo_inspect(smiles)
        stereo_checks[mid] = diag
        m["representation_type"] = "smiles"
        m["smiles"] = f(
            {"source": "source name/drawing", "candidate": smiles}, smiles, ev,
            method="structure_derivation",
            derivation="Candidate connectivity transcribed from the cited source name/drawing; RDKit parsing is a graph check, not source verification.",
        )
        m["stereochemistry"]["atom_centers"] = inferred(
            "RDKit graph diagnostic", diag["atom_centers"], ev,
            "Computed from the candidate SMILES; the toolkit output does not establish source identity.")
        m["stereochemistry"]["double_bonds"] = (
            f(double_bonds, double_bonds, ev) if double_bonds is not None
            else inferred("RDKit graph diagnostic", diag["double_bonds"], ev,
                          "Computed from the candidate SMILES; cumulated/allene and aromatic handling remains toolkit-limited."))
        m["stereochemistry"]["axial"] = (
            inferred("source drawing and graph", axial, ev,
                     "Axial status derived from the source drawing and candidate connectivity; toolkit axial support is limited.")
            if axial is not None else missing(ev))
        m["stereochemistry"]["mixture"] = missing(ev)
        m["stereochemistry"]["assignment"] = (
            f(assignment, assignment, ev) if assignment is not None else missing(ev))
        m["stereochemistry"]["chirality_status"] = inferred(
            "source drawing/name plus graph diagnostic", chirality, ev,
            stereo_note or "Status preserves the source assignment and candidate-graph limitations; absence of a tag is not treated as racemic.")
        m["stereochemistry"]["geometry_status"] = inferred(
            "source drawing/name plus graph diagnostic", geometry, ev,
            "Status preserves explicit source geometry when present; otherwise it is not invented.")
        m["stereochemistry"]["representation_status"] = inferred(
            "candidate SMILES representation", representation, ev,
            "Compared with the source-level stereo statement; toolkit parsing alone is not verification.")
        m["stereochemistry"]["isomeric_smiles"] = inferred(
            "RDKit canonical isomeric SMILES", diag["canonical_isomeric_smiles"], ev,
            "Canonicalized from the unverified candidate SMILES.")
        m["stereochemistry"]["identity_key"] = inferred(
            "RDKit stereo graph key", diag["stereo_graph_key"], ev,
            "Hash of the candidate canonical isomeric graph; not a sample identity.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing(ev)
        m["stereochemistry"]["atom_centers"] = missing(ev)
        m["stereochemistry"]["double_bonds"] = (
            f(double_bonds, double_bonds, ev) if double_bonds is not None else missing(ev))
        m["stereochemistry"]["axial"] = (
            f(axial, axial, ev) if axial is not None else missing(ev))
        m["stereochemistry"]["mixture"] = missing(ev)
        m["stereochemistry"]["assignment"] = (
            f(assignment, assignment, ev) if assignment is not None else missing(ev))
        for key, value in (("chirality_status", chirality),
                           ("geometry_status", geometry),
                           ("representation_status", representation)):
            m["stereochemistry"][key] = inferred(
                "source drawing/name reviewed; complete stereo graph not encoded", value, ev,
                stereo_note or "The source-level status is retained without manufacturing a complete stereochemical structure.")
        m["stereochemistry"]["isomeric_smiles"] = missing(ev)
        m["stereochemistry"]["identity_key"] = missing(ev)
    m["identity_note"] = f(
        reason or "Machine-extracted literature identity; structure remains subject to expert review.",
        reason or "Machine-extracted literature identity; structure remains subject to expert review.", ev)
    m["external_properties"] = []
    molecules[mid] = m
    bundle["molecules"].append(m)
    return m


# Core substrates/products with candidate structures where the source is sufficiently explicit.
M1A = f"{DOC}:mol:enamine_1a_Z"
M2A = f"{DOC}:mol:allenoate_2a"
M3AA = f"{DOC}:mol:product_3aa"
add_molecule(
    M1A, "1a", [P[2], T[2], P[3]],
    name="ethyl (Z)-2-((4-methylphenyl)sulfonamido)-3-phenylacrylate",
    smiles="CCOC(=O)/C(NS(=O)(=O)c1ccc(C)cc1)=C/c1ccccc1",
    chirality="achiral", geometry="specified", representation="encoded",
    double_bonds=[{"source_assignment": "Z", "alkene": "enamine C=C"}],
    assignment="The prose explicitly names substrate 1a as the (Z) isomer.",
    reason="The Z assignment is explicit in page-2 prose; the candidate isomeric SMILES encodes that assignment.")

allenoate_smiles = {
    "2a": "C=C=C(COC(C)=O)C(=O)OCc1ccccc1",
    "2b": "C=C=C(COC(C)=O)C(=O)OC",
    "2c": "C=C=C(COC(C)=O)C(=O)OCC",
    "2d": "C=C=C(COC(C)=O)C(=O)OC(C)C",
    "2e": "C=C=C(COC(C)=O)C(=O)OC(C)(C)C",
    "2f": "C=C=C(COC(C)=O)C(=O)OCCc1ccccc1",
    "2g": "C=C=C(COC(C)=O)C(=O)OCc1ccccc1C",
    "2h": "C=C=C(COC(C)=O)C(=O)OCc1ccccc1F",
    "2i": "C=C=C(COC(C)=O)C(=O)OCc1cccc2ccccc12",
    "2j": "C=C=C(COC(C)=O)C(=O)OCc1c2ccccc2cc2ccccc12",
}
allenoates = {}
for label, smi in allenoate_smiles.items():
    ev = [P[2], P[3]] if label == "2a" else [P[3]]
    mid = f"{DOC}:mol:allenoate_{label}"
    allenoates[label] = mid
    add_molecule(
        mid, label, ev,
        name=("benzyl 2-(acetoxymethyl)buta-2,3-dienoate" if label == "2a"
              else f"beta-prime-acetoxy allenoate {label} (structure drawn in Table 3)"),
        smiles=smi, chirality="achiral", geometry="not_applicable",
        representation="partial",
        axial=[{"type": "allene", "status": "not stereogenic", "basis": "one terminal carbon is CH2"}],
        stereo_note="The drawn allene has a terminal CH2 group, so an axial stereoisomer is not created; no racemic state is inferred.")
M2A = allenoates["2a"]
add_molecule(
    M3AA, "3aa", [P[2], P[3]],
    name="5-benzyl 2-ethyl [1,1'-biphenyl]-2,5-dicarboxylate",
    smiles="CCOC(=O)c1ccc(C(=O)OCc2ccccc2)c(-c2ccccc2)c1",
    chirality="unspecified", geometry="not_applicable", representation="partial",
    reason="Candidate constitutional SMILES is based on the page-2 name and drawing; possible atropisomerism is not assigned by the source.")


def unresolved(mid_suffix, label, ev, name=None, *, chirality="unresolved",
               geometry="unresolved", representation="not_encoded",
               assignment=None, reason=None):
    return add_molecule(f"{DOC}:mol:{mid_suffix}", label, ev, name=name,
                        chirality=chirality, geometry=geometry,
                        representation=representation, assignment=assignment,
                        reason=reason).get("molecule_id")


# Common agents.  Their exact source labels are preserved; structures are not forced merely to silence audits.
agent_specs = {
    "PPh3": "triphenylphosphine", "Ph2PEt": "ethyldiphenylphosphine",
    "Cs2CO3": "cesium carbonate", "K2CO3": "potassium carbonate",
    "Na2CO3": "sodium carbonate", "LiOH": "lithium hydroxide",
    "K3PO4": "potassium phosphate", "p-xylene": "p-xylene",
    "mesitylene": "mesitylene", "trichlorotoluene": "trichlorotoluene (isomer not specified)",
    "DMSO": "dimethyl sulfoxide", "DMF": "N,N-dimethylformamide",
    "toluene": "toluene", "tri_p_fluorophenyl_phosphine": "tris(4-fluorophenyl)phosphine",
    "tri_p_methoxyphenyl_phosphine": "tris(4-methoxyphenyl)phosphine",
    "PhEt2P": "diethylphenylphosphine", "Bu3P": "tributylphosphine",
    "petroleum_ether": "petroleum ether", "ethyl_acetate": "ethyl acetate",
    "THF": "tetrahydrofuran", "MeOH": "methanol", "KOH": "potassium hydroxide",
    "H2O": "water", "L_methionine_methyl_ester_HCl": "L-methionine methyl ester hydrochloride",
    "BOP": "BOP coupling reagent", "Et3N": "triethylamine", "LiAlH4": "lithium aluminium hydride",
    "Pd_PPh3_4": "tetrakis(triphenylphosphine)palladium(0)", "CsF": "cesium fluoride",
    "p_anisylboronic_acid": "4-methoxyphenylboronic acid", "D2O": "deuterium oxide",
    "ethanol": "ethanol", "ethylene": "ethylene", "Ru_catalyst": "ruthenium catalyst (unspecified)",
    "amine_generic": "amine (generic)", "Ni_catalyst": "nickel catalyst (unspecified)",
}
agents = {}
for key, name in agent_specs.items():
    ev = [P[2]]
    if key in {"petroleum_ether", "ethyl_acetate"}: ev = [P[5]]
    elif key in {"THF", "MeOH", "KOH", "H2O", "L_methionine_methyl_ester_HCl", "BOP", "Et3N", "LiAlH4", "Pd_PPh3_4", "CsF", "p_anisylboronic_acid", "D2O"}: ev = [P[4]]
    chirality = "specified" if key == "L_methionine_methyl_ester_HCl" else "unresolved"
    assignment = "L-methionine-derived stereocenter; absolute R/S is not explicitly stated on pages 1–5." if key == "L_methionine_methyl_ester_HCl" else None
    if key == "p_anisylboronic_acid":
        mid = f"{DOC}:mol:agent_{key}"
        add_molecule(mid, key, ev, name=name, smiles="COc1ccc(B(O)O)cc1",
                     chirality="achiral", geometry="not_applicable", representation="encoded",
                     reason="Candidate connectivity is directly transcribed from the Scheme-2 formula 4-MeOC6H4B(OH)2; no external lookup was used.")
        agents[key] = mid
    else:
        agents[key] = unresolved(f"agent_{key}", key if key not in {"tri_p_fluorophenyl_phosphine", "tri_p_methoxyphenyl_phosphine"} else name,
                                 ev, name=name, chirality=chirality,
                                 geometry="not_applicable", assignment=assignment,
                                 reason="Identity is source-named; no external CAS or structure lookup was performed.")


def participant(pid, mid, role, ev, *, amount=None, equivalents=None,
                mol_percent=None, volume=None, eqref=None, product_role=None,
                stage="reaction"):
    p = st.new_participant(pid, mid, role, ev, stage=stage)
    if product_role:
        p["product_role"] = f(product_role, product_role, ev)
    else:
        p["product_role"] = na(ev) if role != "product" else missing(ev)
    for key, value in (("amount", amount), ("equivalents", equivalents),
                       ("mol_percent", mol_percent), ("volume", volume)):
        if value is not None:
            p["charges"][key] = value
    if eqref:
        p["equivalence_reference"] = f(eqref, eqref, ev)
    return p


def step(step_id, ev, participant_ids, operation, *, temp=None, time=None,
         atmosphere=None, stage="reaction"):
    s = st.new_step(step_id, ev, stage=stage)
    s["participant_ids"] = list(participant_ids)
    s["operation"] = f(operation, operation, ev)
    s["temperature"] = temp if temp is not None else missing(ev)
    s["time"] = time if time is not None else missing(ev)
    s["atmosphere"] = atmosphere if atmosphere is not None else missing(ev)
    s["pressure"] = missing(ev)
    return s


def outcome(oid, kind, targets, ev, raw, norm, *, yield_type=None,
            basis=None, scope=None, role="main product"):
    o = st.new_outcome(oid, kind, targets, ev)
    o["product_role"] = f(role, role, ev) if targets else na(ev)
    o["value"] = f(raw, norm, ev)
    o["yield_type"] = f(yield_type, yield_type, ev) if yield_type else na(ev)
    o["measurement_basis"] = f(basis, basis, ev) if basis else missing(ev)
    o["scope"] = f(scope, scope, ev) if scope else missing(ev)
    return o


def classification_unresolved(ev, why):
    c = st.new_reaction(f"{DOC}:temporary", DOC, ev)["classification"]
    c["coverage"] = inferred(why, "unresolved", ev, why)
    c["labels"] = missing(ev)
    c["limitations"] = [why]
    return c


def classification_na(ev, why):
    c = st.new_reaction(f"{DOC}:temporary", DOC, ev)["classification"]
    c["coverage"] = inferred(why, "not_applicable", ev, why)
    c["labels"] = inferred(why, [], ev, why)
    c["limitations"] = [why]
    return c


def classification_changes(ev, changes, label_specs, why):
    c = st.new_reaction(f"{DOC}:temporary", DOC, ev)["classification"]
    c["bond_changes"] = [inferred("source transformation", x, ev, why) for x in changes]
    labels = [{"category_number": n, "code": code, "supporting_change_ids": ids}
              for n, code, ids in label_specs]
    c["labels"] = inferred("source transformation", labels, ev, why)
    c["coverage"] = inferred("partial source-derived classification", "partial", ev, why)
    c["limitations"] = [why, "No atom-mapped reaction SMILES is asserted; exact atom correspondence and all bond-order/hydrogen changes remain unresolved."]
    return c


def cc2(ev):
    changes = [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
    ]
    return classification_changes(ev, changes, [(2, "C-C", ["b1", "b2"])],
                                  "The source identifies a [4+2]/domino benzannulation, supporting two new C-C ring-closing connections; aromatic bond-order changes are not safely mapped.")


def finish_reaction(rid, ev, *, record_type="specific_experiment", categories=None,
                    location=None, name=None, transformation=None,
                    participants=(), steps=(), outcomes=(), classification=None,
                    notes=(), unresolved_items=(), general=(), solvent=None,
                    procedure=None, scale=None):
    r = st.new_reaction(rid, DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = (f(categories, categories, ev) if categories is not None else na(ev))
    r["location"] = f(location, location, ev)
    r["general_scheme_ids"] = list(general)
    r["reaction_name"] = f(name, name, ev) if name else missing(ev)
    r["transformation_type"] = (inferred("source reaction depiction", transformation, ev,
                                         "Normalized transformation type inferred from the source depiction/name.")
                                     if transformation else missing(ev))
    r["participants"] = list(participants)
    r["steps"] = list(steps)
    r["outcomes"] = list(outcomes)
    r["classification"] = classification or classification_unresolved(ev, "Insufficient source detail for safe bond-change assignment.")
    r["notes"] = [f(x, x, ev) for x in notes]
    r["reaction_scale"] = scale if scale is not None else missing(ev)
    if solvent:
        r["solvent_conditions"]["composition"] = f(solvent, solvent, ev)
    r["procedure"]["original_text"] = f(procedure, procedure, ev) if procedure else missing(ev)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": list(unresolved_items)}
    r["completeness"] = "partial"
    bundle["reactions"].append(r)
    return r


# Scheme 1: three cited examples and the authors' generic proposal.
lit_mols = {}
for key, label in {
    "TTMA": "TTMA (source abbreviation; exact expansion not shown)",
    "lit_terephthalate": "terephthalate product (R groups generic)",
    "dimethyl_terephthalate": "dimethyl terephthalate",
    "amino_terephthalate": "amino-substituted dimethyl terephthalate",
    "ethyl_propiolate": "ethyl propiolate",
    "cyclotrimer_product": "carboxylate-substituted aromatic product",
    "this_work_enamine_template": "enamine 1 (R1/R2 generic)",
    "this_work_allenoate_template": "beta-prime-acetoxy allenoate 2 (R1/R3 generic)",
    "this_work_product_template": "substituted terephthalate 3 (R groups generic)",
}.items():
    lit_mols[key] = unresolved(key, label, [P[2]], name=label,
                               chirality="unresolved", geometry="unresolved",
                               reason="Generic or abbreviated Scheme-1 structure; variables were not expanded into a concrete SMILES.")

lit_specs = [
    ("scheme1_prior_diels_alder", [lit_mols["TTMA"], agents["ethanol"], agents["ethylene"], lit_mols["lit_terephthalate"]],
     "Diels-Alder/esterification/dehydrogenation route", "Scheme 1a, upper prior-work example"),
    ("scheme1_prior_CH_functionalization", [lit_mols["dimethyl_terephthalate"], agents["Ru_catalyst"], agents["amine_generic"], lit_mols["amino_terephthalate"]],
     "C-H functionalization", "Scheme 1a, middle prior-work example"),
    ("scheme1_prior_222_benzannulation", [lit_mols["ethyl_propiolate"], agents["Ni_catalyst"], lit_mols["cyclotrimer_product"]],
     "[2+2+2] benzannulation", "Scheme 1a, lower prior-work example"),
]
for suffix, mids, rname, loc in lit_specs:
    ps = []
    for i, mid in enumerate(mids):
        role = "product" if i == len(mids)-1 else ("catalyst" if "catalyst" in mid else "reactant")
        ps.append(participant(f"p{i+1}", mid, role, [P[2]], product_role="depicted product" if role == "product" else None))
    finish_reaction(f"{DOC}:{suffix}", [P[2]], record_type="literature_example",
                    location=loc, name=rname, participants=ps,
                    classification=classification_unresolved([P[2]], "Cited overview lacks enough atom-level detail for safe mapping."),
                    unresolved_items=["Cited prior work is represented only at Scheme-1 overview resolution; no experimental procedure was imported."],
                    procedure="Overview of cited prior work; not an experiment reported in this paper.")

GENERIC = f"{DOC}:scheme1_this_work_generic"
finish_reaction(
    GENERIC, [P[2]], record_type="generic_scheme", location="Scheme 1b",
    name="De novo construction of substituted terephthalates",
    transformation="phosphine-catalyzed domino benzannulation",
    participants=[
        participant("enamine", lit_mols["this_work_enamine_template"], "reactant", [P[2]]),
        participant("allenoate", lit_mols["this_work_allenoate_template"], "reactant", [P[2]]),
        participant("PPh3", agents["PPh3"], "catalyst", [P[2]]),
        participant("product", lit_mols["this_work_product_template"], "product", [P[2]], product_role="generic product"),
    ], classification=cc2([P[2]]),
    unresolved_items=["R1/R2/R3 are variables; the generic scheme is not a concrete experiment."],
    procedure="Scheme 1b generic reaction; PPh3 shown above the arrow and 21 examples stated.")


# Table 1 optimization rows.
opt_rows = [
    (1, "Ph2PEt", "Cs2CO3", "toluene", 62, 5),
    (2, "Ph2PEt", "K2CO3", "toluene", 43, 5),
    (3, "Ph2PEt", "Na2CO3", "toluene", 68, 5),
    (4, "Ph2PEt", "LiOH", "toluene", 58, 5),
    (5, "Ph2PEt", "K3PO4", "toluene", 66, 5),
    (6, "Ph2PEt", "Na2CO3", "p-xylene", 13, 5),
    (7, "Ph2PEt", "Na2CO3", "mesitylene", 53, 5),
    (8, "Ph2PEt", "Na2CO3", "trichlorotoluene", 22, 5),
    (9, "Ph2PEt", "Na2CO3", "DMSO", 47, 24),
    (10, "Ph2PEt", "Na2CO3", "DMF", None, 24),
    (11, "tri_p_fluorophenyl_phosphine", "Cs2CO3", "toluene", 44, 5),
    (12, "tri_p_methoxyphenyl_phosphine", "Cs2CO3", "toluene", 17, 5),
    (13, "PhEt2P", "Cs2CO3", "toluene", 36, 5),
    (14, "Bu3P", "Cs2CO3", "toluene", 28, 5),
    (15, "PPh3", "Cs2CO3", "toluene", 67, 5),
    (16, "PPh3", "Na2CO3", "toluene", 75, 5),
    (17, "PPh3", None, "toluene", 19, 5),
    (18, None, "Cs2CO3", "toluene", 5, 5),
]
opt_ids = []
for n, cat, base, solv, yld, hours in opt_rows:
    ev = [P[2]]
    ps = [
        participant("enamine_1a", M1A, "reactant", ev, amount=q("0.1 mmol", 0.1, "mmol", ev)),
        participant("allenoate_2a", M2A, "reactant", ev, amount=q("0.12 mmol", 0.12, "mmol", ev)),
    ]
    if cat:
        ps.append(participant("catalyst", agents[cat], "catalyst", ev,
                              amount=q("0.02 mmol", 0.02, "mmol", ev), eqref=M1A))
    if base:
        ps.append(participant("base", agents[base], "reagent", ev,
                              amount=q("0.12 mmol", 0.12, "mmol", ev), eqref=M1A))
    ps += [participant("solvent", agents[solv], "solvent", ev,
                       volume=q("1 mL", 1, "mL", ev)),
           participant("product_3aa", M3AA, "product", ev, product_role="desired product")]
    out = (outcome("yield", "yield", [M3AA], ev, f"{yld}%", {"number": yld, "unit": "%"},
                   yield_type="isolated yield", basis="isolated", scope=f"Table 1 entry {n}")
           if yld is not None else
           outcome("table_dash", "other", [M3AA], ev, "table dash (–)", "table dash (–)",
                   scope=f"Table 1 entry {n}", role="desired product; yield not reported"))
    rid = f"{DOC}:table1_entry_{n:02d}"
    opt_ids.append(rid)
    finish_reaction(
        rid, ev, categories=["optimization"], location=f"Table 1 entry {n}",
        name="screening conditions of the reaction",
        transformation="domino [4+2] benzannulation/aromatization",
        participants=ps,
        steps=[step("reaction", ev, [p["participant_id"] for p in ps if p["role"] != "product"],
                    "stirred", temp=f("room temperature", "room temperature", ev),
                    time=q(f"{hours} h", hours, "h", ev))],
        outcomes=[out], classification=(cc2(ev) if yld is not None else classification_unresolved(ev, "The dash does not establish an observed transformation.")),
        general=[GENERIC], solvent=solv,
        scale=f("1a (0.1 mmol)", {"reference_molecule_id": M1A, "amount": {"number": 0.1, "unit": "mmol"}}, ev),
        procedure="Unless otherwise noted: 1a (0.1 mmol), 2a (0.12 mmol), phosphine (0.02 mmol), base (0.12 mmol), solvent (1 mL), room temperature; entries 9 and 10 use 24 h, others 5 h.",
        unresolved_items=(["Meaning of the Table-1 dash is not specified; it was not converted to zero or no reaction."] if yld is None else []) +
                         ["Optimization-table workup details are not reported on pages 1–5; SI was not supplied."])


# Table 2 enamine scope.
table2 = [
    ("3aa", "1a", 75), ("3ba", "1b", 64), ("3ca", "1c", 38),
    ("3da", "1d", 63), ("3ea", "1e", 55), ("3fa", "1f", 61),
    ("3ga", "1g", 71), ("3ha", "1h", 90), ("3ia", "1i", 79),
    ("3ja", "1j", 71), ("3ka", "1k", 53), ("3la", "1l", 83),
    ("3ma", "1m", 93), ("3na", "1n", 0),
]
enamines = {"1a": M1A}
products = {"3aa": M3AA}
for _, label, _ in table2:
    if label not in enamines:
        enamines[label] = unresolved(f"enamine_{label}_table2", label, [P[3]],
                                     name=f"enamine {label} (structure drawn in Table 2)",
                                     chirality="unspecified", geometry="unspecified",
                                     reason="The alkene structure is visible, but no E/Z descriptor is explicitly assigned for this individual entry on pages 1–5.")
for product_label, _, _ in table2:
    if product_label not in products:
        if product_label == "3ha":
            mid = f"{DOC}:mol:product_{product_label}"
            add_molecule(mid, product_label, [P[3], P[4]],
                         name="methyl benzyl chlorobiphenyl-2,5-dicarboxylate 3ha",
                         smiles="COC(=O)c1ccc(C(=O)OCc2ccccc2)c(-c2cccc(Cl)c2)c1",
                         chirality="unspecified", geometry="not_applicable", representation="partial",
                         reason="Candidate connectivity follows the Table-2/Scheme-2 drawing; no axial configuration is assigned.")
            products[product_label] = mid
        else:
            products[product_label] = unresolved(f"product_{product_label}", product_label, [P[3]],
                                                 name=f"substituted terephthalate {product_label} (structure drawn in Table 2)",
                                                 chirality="unspecified", geometry="not_applicable",
                                                 reason="Constitutional drawing is present; a complete machine SMILES was withheld where substituent position/identity could not be transcribed with high confidence from the main-page image.")

PET = agents["petroleum_ether"]
EA = agents["ethyl_acetate"]


def scope_amount_conflict(label, table_value, gp_value, table_ev, gp_ev):
    return conflict(
        {"table": f"{label} ({table_value} mmol)", "general_procedure": f"{label} ({gp_value} mmol)"},
        [q(f"{table_value} mmol", table_value, "mmol", table_ev),
         f(f"{gp_value} mmol", {"number": gp_value, "unit": "mmol"}, gp_ev,
           method="inherited_general_condition")], table_ev + gp_ev)


def scope_reaction(rid, location, substrate_mid, allene_mid, product_mid, yld, ev_table,
                   *, zero=False):
    ev = ev_table + [P[5], T[5]] if not zero else ev_table
    s1 = participant("substrate_1", substrate_mid, "reactant", ev_table,
                     amount=(q("0.1 mmol", 0.1, "mmol", ev_table) if zero else
                             scope_amount_conflict("enamine 1", 0.1, 0.12, ev_table, [P[5], T[5]])))
    s2 = participant("allenoate_2", allene_mid, "reactant", ev_table,
                     amount=(q("0.12 mmol", 0.12, "mmol", ev_table) if zero else
                             scope_amount_conflict("allenoate 2", 0.12, 0.1, ev_table, [P[5], T[5]])))
    cat = participant("PPh3", agents["PPh3"], "catalyst", ev_table,
                      amount=q("0.02 mmol", 0.02, "mmol", ev_table),
                      mol_percent=q("20 mol%", 20, "mol%", ev_table), eqref=substrate_mid)
    base = participant("Na2CO3", agents["Na2CO3"], "reagent", ev_table,
                       amount=q("0.12 mmol", 0.12, "mmol", ev_table),
                       equivalents=q("1.2 equiv", 1.2, "equiv", ev_table), eqref=substrate_mid)
    solv = participant("toluene", agents["toluene"], "solvent", ev_table,
                       volume=q("1 mL", 1, "mL", ev_table))
    ps = [s1, s2, cat, base, solv]
    outs = []
    if zero:
        outs.append(outcome("no_product", "no_reaction", [], ev_table,
                            "0%; product could not be obtained", True,
                            basis="isolated-yield reporting context", scope=location,
                            role="desired product not obtained"))
        cls = classification_na(ev_table, "No product was obtained, so the intended benzannulation is not assigned a bond-formation label.")
    else:
        ps += [participant("product", product_mid, "product", ev_table, product_role="main product"),
               participant("petroleum_ether", PET, "purification_agent", [P[5]], stage="purification"),
               participant("ethyl_acetate", EA, "purification_agent", [P[5]], stage="purification")]
        outs.append(outcome("yield", "yield", [product_mid], ev_table, f"{yld}%", {"number": yld, "unit": "%"},
                            yield_type="isolated yield", basis="isolated", scope=location))
        cls = cc2(ev_table)
    steps = [step("reaction", ev_table, ["substrate_1", "allenoate_2", "PPh3", "Na2CO3", "toluene"],
                  "stirred", temp=f("room temperature", "room temperature", ev_table),
                  time=q("5 h", 5, "h", ev_table))]
    if not zero:
        steps += [
            step("solvent_removal", [P[5], T[5]], ["toluene"], "removed solvent under reduced pressure", stage="workup"),
            step("flash_chromatography", [P[5], T[5]], ["petroleum_ether", "ethyl_acetate"],
                 "flash column chromatography; petroleum ether/ethyl acetate = 5:1", stage="purification"),
        ]
    notes = []
    unresolved_items = ["SI was not supplied; substrate preparation, full workup details, and characterization remain unavailable."]
    if not zero:
        notes.append("Source conflict: the table footnote charges enamine 1 at 0.1 mmol and allenoate 2 at 0.12 mmol, while the page-5 General Procedure reverses them to 0.12 and 0.1 mmol, respectively.")
        notes.append("The table-local catalyst/base loadings and concentration use substrate 1 as the apparent 0.1 mmol basis; the General Procedure instead declares allenoate 2 as the equivalence basis.")
        unresolved_items.append("Reactant amounts/equivalence basis conflict between the scope-table footnote and General Procedure; neither alternative was selected as the normalized amount.")
    r = finish_reaction(
        rid, ev, categories=["substrate_scope"], location=location,
        name="substrate scope of phosphine-catalyzed domino benzannulation",
        transformation="domino [4+2] benzannulation/aromatization",
        participants=ps, steps=steps, outcomes=outs, classification=cls,
        notes=notes, unresolved_items=unresolved_items, general=[GENERIC], solvent="toluene",
        procedure=("Table scope conditions; successful entries also inherit the page-5 General Procedure workup/purification."
                   if not zero else "Table scope conditions; product could not be obtained."),
        scale=(f("1 (0.1 mmol)", {"reference_molecule_id": substrate_mid, "amount": {"number": 0.1, "unit": "mmol"}}, ev_table) if zero else
               conflict({"table": "1 (0.1 mmol)", "general_procedure": "2 (0.1 mmol)"},
                        [f("1 (0.1 mmol)", {"reference_molecule_id": substrate_mid, "amount": {"number": 0.1, "unit": "mmol"}}, ev_table),
                         f("2 (0.1 mmol)", {"reference_molecule_id": allene_mid, "amount": {"number": 0.1, "unit": "mmol"}}, [P[5], T[5]], method="inherited_general_condition")], ev)))
    if not zero:
        r["procedure"]["general_procedure_reference"] = f("General Procedure for the Synthesis of 3", "page-5 General Procedure", [P[5], T[5]], method="inherited_general_condition")
        r["solvent_conditions"]["total_volume"] = f("1 mL", {"number": 1, "unit": "mL"}, [P[5], T[5]], method="inherited_general_condition")
        r["solvent_conditions"]["concentration"] = f("toluene (0.1 M)", {"number": 0.1, "unit": "M"}, ev_table)
        r["solvent_conditions"]["concentration_basis"] = conflict(
            {"table_basis": substrate_mid, "general_procedure_basis": allene_mid},
            [f("table 0.1 mmol substrate 1 in 1 mL", substrate_mid, ev_table),
             f("General Procedure 0.1 mmol allenoate 2 in 1 mL", allene_mid, [P[5], T[5]], method="inherited_general_condition")], ev)
    return r


table2_ids = []
for product_label, substrate_label, yld in table2:
    rid = f"{DOC}:table2_{product_label}"
    table2_ids.append(rid)
    scope_reaction(rid, f"Table 2, {product_label}", enamines[substrate_label], M2A,
                   products[product_label], yld, [P[3]], zero=(yld == 0))


# Table 3 allenoate scope.
table3 = [("3ab", "2b", 61), ("3ac", "2c", 57), ("3ad", "2d", 29),
          ("3ae", "2e", 22), ("3af", "2f", 46), ("3ag", "2g", 50),
          ("3ah", "2h", 60), ("3ai", "2i", 47), ("3aj", "2j", 53)]
for product_label, _, _ in table3:
    products[product_label] = unresolved(f"product_{product_label}", product_label, [P[3]],
                                         name=f"substituted terephthalate {product_label} (structure drawn in Table 3)",
                                         chirality="unspecified", geometry="not_applicable",
                                         reason="Product identity is source-labelled; complete candidate SMILES was withheld pending expert structure transcription.")
table3_ids = []
for product_label, allene_label, yld in table3:
    rid = f"{DOC}:table3_{product_label}"
    table3_ids.append(rid)
    scope_reaction(rid, f"Table 3, {product_label}", M1A, allenoates[allene_label],
                   products[product_label], yld, [P[3]])


# Equation 1, methanesulfonyl-protected enamine.  Its printed label duplicates Table-2 1n.
MEQ = unresolved("enamine_eq1_1n_methanesulfonyl", "1n (eq 1)", [P[3]],
                  name="methanesulfonyl-protected phenyl enamine shown in eq 1",
                  chirality="unspecified", geometry="unspecified",
                  reason="The source reuses label 1n for a methanesulfonyl substrate even though Table 2 uses 1n for a p-OMe substrate; the identities are kept separate.")
eq_ev = [P[3], P[5], T[5]]
eq_ps = [
    participant("enamine", MEQ, "reactant", [P[3]], amount=f("0.12 mmol", {"number": 0.12, "unit": "mmol"}, [P[5], T[5]], method="inherited_general_condition"), equivalents=f("1.2 equiv", {"number": 1.2, "unit": "equiv"}, [P[5], T[5]], method="inherited_general_condition"), eqref=M2A),
    participant("allenoate", M2A, "reactant", [P[3]], amount=f("0.1 mmol", {"number": 0.1, "unit": "mmol"}, [P[5], T[5]], method="inherited_general_condition"), equivalents=f("1.0 equiv", {"number": 1.0, "unit": "equiv"}, [P[5], T[5]], method="inherited_general_condition"), eqref=M2A),
    participant("PPh3", agents["PPh3"], "catalyst", [P[3]], amount=f("0.02 mmol", {"number": 0.02, "unit": "mmol"}, [P[5], T[5]], method="inherited_general_condition"), mol_percent=q("20 mol%", 20, "mol%", [P[3]]), eqref=M2A),
    participant("Na2CO3", agents["Na2CO3"], "reagent", [P[3]], amount=f("0.12 mmol", {"number": 0.12, "unit": "mmol"}, [P[5], T[5]], method="inherited_general_condition"), equivalents=q("1.2 equiv", 1.2, "equiv", [P[3]]), eqref=M2A),
    participant("toluene", agents["toluene"], "solvent", [P[3]], volume=f("1 mL", {"number": 1, "unit": "mL"}, [P[5], T[5]], method="inherited_general_condition")),
    participant("product", M3AA, "product", [P[3]], product_role="main product"),
    participant("petroleum_ether", PET, "purification_agent", [P[5]], stage="purification"),
    participant("ethyl_acetate", EA, "purification_agent", [P[5]], stage="purification"),
]
eq_rid = f"{DOC}:eq1_methanesulfonyl_enamine"
eq_r = finish_reaction(
    eq_rid, eq_ev, categories=["substrate_scope"], location="Equation 1",
    name="benzannulation with methanesulfonyl-protected enamine",
    transformation="domino [4+2] benzannulation/aromatization",
    participants=eq_ps,
    steps=[step("reaction", [P[3]], ["enamine", "allenoate", "PPh3", "Na2CO3", "toluene"], "stirred", temp=f("rt", "room temperature", [P[3]]), time=q("5 h", 5, "h", [P[3]])),
           step("solvent_removal", [P[5]], ["toluene"], "removed solvent under reduced pressure", stage="workup"),
           step("flash_chromatography", [P[5]], ["petroleum_ether", "ethyl_acetate"], "flash column chromatography; petroleum ether/ethyl acetate = 5:1", stage="purification")],
    outcomes=[outcome("yield", "yield", [M3AA], [P[3]], "68%", {"number": 68, "unit": "%"}, yield_type="isolated yield", basis="isolated", scope="eq 1")],
    classification=cc2([P[3]]), notes=["The printed label 1n conflicts with Table 2, where 1n denotes the p-OMe enamine; separate molecule identities are retained."],
    unresolved_items=["The equation does not give absolute charges; page-5 General Procedure amounts and purification are inherited because the depicted conditions otherwise match.", "SI was not supplied."],
    general=[GENERIC], solvent="toluene", procedure="Eq 1 conditions plus matching General Procedure workup/purification.",
    scale=f("allenoate 2a (0.1 mmol), inherited from General Procedure", {"reference_molecule_id": M2A, "amount": {"number": 0.1, "unit": "mmol"}}, [P[5], T[5]], method="inherited_general_condition"))
eq_r["procedure"]["general_procedure_reference"] = f("General Procedure for the Synthesis of 3", "page-5 General Procedure", [P[5], T[5]], method="inherited_general_condition")


# Scheme 2 products and agents.
M1H = enamines["1h"]
M3HA = products["3ha"]
M3FA = products["3fa"]
MACID = f"{DOC}:mol:product_3ha_acid"
add_molecule(MACID, "hydrolyzed 3ha acid", [P[4]],
             name="selectively hydrolyzed carboxylic acid derived from 3ha",
             smiles="COC(=O)c1ccc(C(=O)O)c(-c2cccc(Cl)c2)c1",
             chirality="unspecified", geometry="not_applicable", representation="partial",
             reason="Candidate acid connectivity is derived from Scheme 2 by replacing the depicted benzyl ester of 3ha with the corresponding carboxylic acid; no axial assignment is made.")
M4 = unresolved("product_4", "4", [P[4]], name="L-methionine methyl ester amide 4",
                chirality="specified", geometry="not_applicable", assignment="L-methionine-derived stereocenter; the main pages do not state an R/S descriptor.",
                reason="The product drawing and L-methionine precursor specify a chiral residue, but absolute configuration was not independently derived into SMILES.")
M5 = unresolved("product_5", "5", [P[4]], name="selectively reduced alcohol 5",
                chirality="unspecified", geometry="not_applicable",
                reason="Product 5 is drawn; complete candidate SMILES was withheld pending expert structure transcription.")
M6 = unresolved("product_6", "6", [P[4]], name="Suzuki-coupling product 6",
                chirality="unspecified", geometry="not_applicable",
                reason="Product 6 is drawn; complete candidate SMILES was withheld pending expert structure transcription.")
MSUZ = unresolved("suzuki_substrate_conflict_3fa_3ha", "3fa in Scheme 2 / 3ha in prose", [P[4], T[4]],
                   name="conflicted Suzuki substrate identity",
                   chirality="unresolved", geometry="not_applicable",
                   reason="Scheme 2 labels and draws brominated 3fa, while the prose names 3ha; the record does not choose between them.")
MD = unresolved("deuterated_benzannulation_product", "deuterated benzannulation product", [P[4], T[4]],
                 name="deuterium-incorporated benzannulation product (positions/extent not quantified)",
                 chirality="unresolved", geometry="not_applicable",
                 reason="Main-text prose reports deuterium incorporation but delegates details to unavailable SI; Scheme 3 marks corresponding mechanistic positions without a quantified final isotopomer.")

scheme2_ids = []
scale_rid = f"{DOC}:scheme2_gram_scale_3ha"
scheme2_ids.append(scale_rid)
scale_ps = [
    participant("enamine_1h", M1H, "reactant", [P[4]], amount=q("4.0 mmol", 4.0, "mmol", [P[4]])),
    participant("allenoate_2a", M2A, "reactant", [P[4]], amount=q("4.0 mmol", 4.0, "mmol", [P[4]])),
    participant("PPh3", agents["PPh3"], "catalyst", [P[4]], mol_percent=q("20 mol%", 20, "mol%", [P[4]]), eqref=M1H),
    participant("Na2CO3", agents["Na2CO3"], "reagent", [P[4]], equivalents=q("1.2 equiv", 1.2, "equiv", [P[4]]), eqref=M1H),
    participant("toluene", agents["toluene"], "solvent", [P[4]]),
    participant("product", M3HA, "product", [P[4]], product_role="main product"),
]
finish_reaction(scale_rid, [P[4]], categories=["scale_up"], location="Scheme 2, gram-scale synthesis",
                name="gram-scale domino benzannulation", transformation="domino [4+2] benzannulation/aromatization",
                participants=scale_ps,
                steps=[step("reaction", [P[4]], ["enamine_1h", "allenoate_2a", "PPh3", "Na2CO3", "toluene"], "stirred", temp=f("RT", "room temperature", [P[4]]), time=q("5 h", 5, "h", [P[4]]))],
                outcomes=[outcome("yield", "yield", [M3HA], [P[4]], "1.48 g, 97%", {"number": 97, "unit": "%"}, scope="1h and 2a to 3ha"),
                          outcome("product_amount", "product_amount", [M3HA], [P[4]], "1.48 g", {"number": 1.48, "unit": "g"}, scope="isolated/obtained product amount")],
                classification=cc2([P[4]]), general=[GENERIC], solvent="toluene",
                scale=f("1h 4.0 mmol and 2a 4.0 mmol", {"reference_molecule_id": M1H, "amount": {"number": 4.0, "unit": "mmol"}}, [P[4]]),
                procedure="Scheme 2: PPh3 (20 mol%), Na2CO3 (1.2 equiv), toluene (0.1 M), RT, 5 h.",
                unresolved_items=["Yield measurement basis and workup/purification are not stated in Scheme 2; SI was not supplied."])

hyd_rid = f"{DOC}:scheme2_hydrolysis_3ha_to_acid"
scheme2_ids.append(hyd_rid)
hyd_ps = [participant("substrate", M3HA, "reactant", [P[4]], amount=q("0.1 mmol", 0.1, "mmol", [P[4]])),
          participant("KOH", agents["KOH"], "reagent", [P[4]]),
          participant("water", agents["H2O"], "reagent", [P[4]]),
          participant("THF", agents["THF"], "solvent", [P[4]]),
          participant("MeOH", agents["MeOH"], "solvent", [P[4]]),
          participant("acid", MACID, "product", [P[4]], product_role="hydrolyzed intermediate")]
hyd_changes = [
    {"change_id": "b1", "atoms": [{"atom_map": None, "element": "O"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain"},
    {"change_id": "b2", "atoms": [{"atom_map": None, "element": "O"}, {"atom_map": None, "element": "C"}], "before_order": 1, "after_order": 0, "change_type": "broken"},
]
finish_reaction(hyd_rid, [P[4], T[4]], categories=["derivatization"], location="Scheme 2, hydrolysis stage",
                name="regioselective ester hydrolysis", transformation="selective ester hydrolysis",
                participants=hyd_ps,
                steps=[step("hydrolysis", [P[4]], ["substrate", "KOH", "water", "THF", "MeOH"], "treated under inorganic basic conditions", temp=q("35 °C", 35, "°C", [P[4]]))],
                outcomes=[outcome("yield", "yield", [MACID], [P[4]], "85%", {"number": 85, "unit": "%"}, scope="3ha to corresponding acid")],
                classification=classification_changes([P[4]], hyd_changes, [(24, "O-H", ["b1"])], "Selective ester hydrolysis implies carboxylic O-H formation and alkoxy C-O cleavage; exact atom mapping is unresolved."),
                solvent="THF/MeOH", scale=f("3ha, 0.1 mmol", {"reference_molecule_id": M3HA, "amount": {"number": 0.1, "unit": "mmol"}}, [P[4]]),
                procedure="1. THF/MeOH, KOH/H2O, 35 °C, 85%.",
                unresolved_items=["Reaction time, reagent quantities, solvent ratio/volume, atmosphere, workup, and yield measurement basis are unreported on pages 1–5."])

couple_rid = f"{DOC}:scheme2_coupling_acid_to_4"
scheme2_ids.append(couple_rid)
couple_ps = [participant("acid", MACID, "reactant", [P[4]]),
             participant("methionine", agents["L_methionine_methyl_ester_HCl"], "reactant", [P[4]]),
             participant("BOP", agents["BOP"], "reagent", [P[4]]),
             participant("Et3N", agents["Et3N"], "reagent", [P[4]]),
             participant("THF", agents["THF"], "solvent", [P[4]]),
             participant("product", M4, "product", [P[4]], product_role="amide product")]
cn_change = [{"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "N"}], "before_order": 0, "after_order": 1, "change_type": "formed"}]
finish_reaction(couple_rid, [P[4], T[4]], categories=["derivatization"], location="Scheme 2, coupling stage",
                name="coupling to L-methionine methyl ester hydrochloride", transformation="amide coupling",
                participants=couple_ps,
                steps=[
                    step("low_temperature_coupling", [P[4]], ["acid", "methionine", "BOP", "Et3N", "THF"], "coupling initiated at 0 °C", temp=q("0 °C", 0, "°C", [P[4]])),
                    step("warm_to_room_temperature", [P[4]], ["acid", "methionine", "BOP", "Et3N", "THF"], "warmed to room temperature; transition/hold time unreported", temp=f("rt", "room temperature", [P[4]])),
                ],
                outcomes=[outcome("yield", "yield", [M4], [P[4]], "90%", {"number": 90, "unit": "%"}, scope="corresponding acid to 4")],
                classification=classification_changes([P[4]], cn_change, [(5, "C-N", ["b1"])], "The depicted amide coupling forms one carbon-nitrogen bond; other activation/dehydration changes are unresolved."),
                solvent="THF", procedure="2. L-methionine methyl ester hydrochloride, BOP, Et3N, THF, 0 °C→rt, 90%.",
                unresolved_items=["The 0 °C→rt warming duration and hold time are not reported; reagent quantities, atmosphere, workup, and yield measurement basis are also unreported."])

red_rid = f"{DOC}:scheme2_reduction_3ha_to_5"
scheme2_ids.append(red_rid)
red_ps = [participant("substrate", M3HA, "reactant", [P[4]], amount=q("0.1 mmol", 0.1, "mmol", [P[4]])),
          participant("LiAlH4", agents["LiAlH4"], "reagent", [P[4]]),
          participant("THF", agents["THF"], "solvent", [P[4]]),
          participant("product", M5, "product", [P[4]], product_role="selectively reduced product")]
red_changes = [
    {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 2},
    {"change_id": "b2", "atoms": [{"atom_map": None, "element": "O"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain"},
]
finish_reaction(red_rid, [P[4], T[4]], categories=["derivatization"], location="Scheme 2, reduction",
                name="LiAlH4 reduction", transformation="selective ester reduction to alcohol",
                participants=red_ps,
                steps=[step("reduction", [P[4]], ["substrate", "LiAlH4", "THF"], "treated with LiAlH4", temp=f("rt", "room temperature", [P[4]]))],
                outcomes=[outcome("yield", "yield", [M5], [P[4]], "75% yield", {"number": 75, "unit": "%"}, scope="3ha to 5")],
                classification=classification_changes([P[4]], red_changes, [(1, "C-H", ["b1"]), (24, "O-H", ["b2"])], "The ester-to-alcohol product drawing requires net C-H and O-H gains; hydride/proton sources and removed atoms are not fully mapped."),
                solvent="THF", scale=f("3ha, 0.1 mmol", {"reference_molecule_id": M3HA, "amount": {"number": 0.1, "unit": "mmol"}}, [P[4]]),
                procedure="LiAlH4, THF, rt.",
                unresolved_items=["Reaction time, LiAlH4 quantity, atmosphere, workup, and yield measurement basis are unreported on pages 1–5."])

suz_rid = f"{DOC}:scheme2_suzuki_to_6"
scheme2_ids.append(suz_rid)
suz_ps = [participant("substrate_conflict", MSUZ, "reactant", [P[4], T[4]], amount=q("0.1 mmol", 0.1, "mmol", [P[4]])),
          participant("boronic_acid", agents["p_anisylboronic_acid"], "reactant", [P[4]]),
          participant("Pd", agents["Pd_PPh3_4"], "catalyst", [P[4]]),
          participant("CsF", agents["CsF"], "reagent", [P[4]]),
          participant("DMF", agents["DMF"], "solvent", [P[4]]),
          participant("product", M6, "product", [P[4]], product_role="Suzuki-coupling product")]
suz_changes = [{"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
               {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "Br"}], "before_order": 1, "after_order": 0, "change_type": "broken"},
               {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "B"}], "before_order": 1, "after_order": 0, "change_type": "broken"}]
finish_reaction(suz_rid, [P[4], T[4]], categories=["derivatization"], location="Scheme 2, Suzuki coupling",
                name="Suzuki coupling", transformation="Suzuki cross-coupling",
                participants=suz_ps,
                steps=[step("coupling", [P[4]], ["substrate_conflict", "boronic_acid", "Pd", "CsF", "DMF"], "Suzuki coupling", temp=q("100 °C", 100, "°C", [P[4]]), atmosphere=f("Ar", "argon", [P[4]]))],
                outcomes=[outcome("yield", "yield", [M6], [P[4]], "64% yield", {"number": 64, "unit": "%"}, scope="conflicted 3fa/3ha substrate to 6")],
                classification=classification_changes([P[4]], suz_changes, [(2, "C-C", ["b1"])], "Suzuki coupling forms one biaryl C-C bond and cleaves the depicted C-Br/C-B connections; exact atom mapping is unresolved."),
                notes=["Source conflict: Scheme 2 labels/draws the brominated starting material as 3fa, while the adjacent prose names 3ha and calls it 4'-bromo. Table 2 depicts 3fa as the bromo entry and 3ha as a chloro entry; no participant identity was silently selected."],
                solvent="DMF", scale=f("Scheme 2 substrate, 0.1 mmol", {"reference_molecule_id": MSUZ, "amount": {"number": 0.1, "unit": "mmol"}}, [P[4]]),
                procedure="Pd(PPh3)4, CsF, DMF, Ar, 100 °C, 4-MeOC6H4B(OH)2.",
                unresolved_items=["Starting-material identity conflicts between Scheme 2 (3fa) and prose (3ha).", "Reaction time, reagent/catalyst quantities, workup, and yield measurement basis are unreported."])

d_rid = f"{DOC}:deuterium_experiment_D2O"
scheme2_ids.append(d_rid)
d_ps = [participant("enamine", M1A, "reactant", [P[4], T[4]], amount=q("0.1 mmol", 0.1, "mmol", [P[2]])),
        participant("allenoate", M2A, "reactant", [P[4], T[4]], amount=q("0.12 mmol", 0.12, "mmol", [P[2]])),
        participant("PPh3", agents["PPh3"], "catalyst", [P[4]], amount=q("0.02 mmol", 0.02, "mmol", [P[2]]), eqref=M1A),
        participant("Na2CO3", agents["Na2CO3"], "reagent", [P[4]], amount=q("0.12 mmol", 0.12, "mmol", [P[2]]), eqref=M1A),
        participant("D2O", agents["D2O"], "additive", [P[4], T[4]], equivalents=q("20 equiv", 20, "equiv", [P[4], T[4]]), eqref=M1A),
        participant("toluene", agents["toluene"], "solvent", [P[4]]),
        participant("deuterated_product", MD, "product", [P[4], T[4]], product_role="deuterium-incorporated product")]
finish_reaction(d_rid, [P[4], T[4]], categories=["mechanistic_study"], location="page 4 deuterium experiment prose; corresponding positions annotated in Scheme 3",
                name="deuterium-labeling experiment", transformation="benzannulation under D2O with deuterium incorporation",
                participants=d_ps,
                steps=[step("reaction", [P[4], T[4], P[2]], ["enamine", "allenoate", "PPh3", "Na2CO3", "D2O", "toluene"], "standard reaction with 20 equiv D2O introduced", temp=f("room temperature", "room temperature", [P[2]]), time=q("5 h", 5, "h", [P[2]]))],
                outcomes=[outcome("deuterium_incorporation", "other", [MD], [P[4], T[4]], "Incorporation of deuterium in the corresponding positions was observed", "deuterium incorporation observed", scope="positions referenced by Scheme 3; degree and analytical method unreported")],
                classification=cc2([P[4], T[4]]),
                notes=["The main text states that further experimental details are in Supporting Information, which was not supplied. Isotope substitution is retained as an outcome and is not forced into an element-only C-H bond label."],
                general=[GENERIC], solvent="toluene",
                procedure="All else being equal, 20 equiv D2O was introduced into the reaction system.",
                unresolved_items=["Deuterium incorporation percentage, exact final-product positions, analytical method, yield, and full procedure are unavailable without SI.", "The 20-equiv denominator is inferred from the standard 0.1 mmol substrate basis and remains subject to review."])

mech_rid = f"{DOC}:scheme3_proposed_mechanism"
finish_reaction(mech_rid, [P[4]], record_type="proposed_mechanism", location="Scheme 3",
                name="plausible mechanism of the reaction of 1a and 2a",
                transformation="proposed phosphine-catalyzed domino benzannulation mechanism",
                participants=[participant("allenoate_structure", M2A, "reactant", [P[4]]),
                              participant("enamine_structure", M1A, "reactant", [P[4]]),
                              participant("PPh3", agents["PPh3"], "catalyst", [P[4]]),
                              participant("product", M3AA, "product", [P[4]], product_role="overall proposed product")],
                classification=cc2([P[4]]),
                notes=["Scheme 3 labels the allenoate drawing as 1a and the enamine drawing as 2a, opposite to Tables 1–3 and the Experimental Section; identities are assigned by structure and the printed-label inversion is retained as a conflict note.",
                       "Proposed sequence: acetate elimination/phosphine attack to A; base-assisted enamine interaction to C; cyclization/elimination/proton transfer through D–J; deprotonation/isomerization to aromatic product with phosphine regeneration."],
                procedure="Author-proposed arrow-pushing mechanism; not an experimentally executed reaction record.",
                unresolved_items=["Intermediates A–J are proposed, not independently observed on pages 1–5.", "Printed substrate labels in Scheme 3 are inverted relative to the rest of the paper."])


# Non-synthetic observations and context, kept out of the reaction count.
bundle["source_context"] = [
    {"context_id": f"{DOC}:context:figure1", "document_id": DOC,
     "kind": "non_synthetic_figure", "location": "Figure 1, page 1",
     "evidence_ids": [P[1], T[1]],
     "content": "Sampling of pharmaceuticals containing terephthalates; structural examples are context, not experiments in this paper."},
    {"context_id": f"{DOC}:context:xray_3ea", "document_id": DOC,
     "kind": "structure_assignment", "location": "page 3 prose and reference 16",
     "evidence_ids": [P[3], T[3], P[5], T[5]],
     "content": "The paper states that a product structure was confirmed by single-crystal X-ray diffraction and cites CCDC 2021955 for 3ea; this is not a separate synthesis experiment."},
    {"context_id": f"{DOC}:context:scheme3_intermediates", "document_id": DOC,
     "kind": "proposed_mechanistic_intermediates", "location": "Scheme 3, page 4",
     "evidence_ids": [P[4], T[4]],
     "content": "Intermediates A–J and red deuterated annotations are author-proposed mechanistic context; they are not independent isolated reactions."},
    {"context_id": f"{DOC}:context:characterization_inventory", "document_id": DOC,
     "kind": "characterization_inventory", "location": "pages 1–5",
     "evidence_ids": [P[i] for i in range(1, 6)],
     "content": "No standalone mass-spectrum peak table, molecular-weight table, or fragment-relative-abundance block appears on pages 1–5."},
]


# Save raw inventory and exact page coverage before validation.
page2_records = [x["reaction_id"] for x in bundle["reactions"] if x["reaction_id"].startswith(f"{DOC}:scheme1_")] + opt_ids
page3_records = table2_ids + table3_ids + [eq_rid]
page4_records = scheme2_ids + [mech_rid]
page5_records = [rid for rid in table2_ids if not rid.endswith("table2_3na")] + table3_ids + [eq_rid]
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title/abstract/keywords and both text columns", "finding": "Article metadata and introductory context; no concrete experiment arrow or table row.", "evidence_ids": [P[1], T[1]]},
            {"region": "Figure 1, lower right", "finding": "Pharmaceutical examples containing terephthalates; non-synthetic context retained in source_context.", "evidence_ids": [P[1]]},
        ], "record_ids": [], "unresolved": ["Supporting information is advertised online but no SI file was supplied or associated."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Scheme 1a, upper left", "finding": "Three cited prior-work transformations recorded as literature_example, not this paper's experiments.", "evidence_ids": [P[2], T[2]]},
            {"region": "Scheme 1b, middle left", "finding": "Authors' generic R-group benzannulation recorded as generic_scheme.", "evidence_ids": [P[2], T[2]]},
            {"region": "Table 1, lower right including footnotes", "finding": "All 18 optimization rows recorded independently; entries 9/10 use 24 h, and entry 10 retains the yield dash without converting it to zero.", "evidence_ids": [P[2], T[2]]},
        ], "record_ids": page2_records, "unresolved": ["Optimization workup details are absent from pages 1–5.", "The exact meaning of Table-1 entry-10 dash is not stated."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Table 2, lower left including footnote", "finding": "Fourteen enamine-scope rows 3aa–3na recorded; 3na is a nonproductive 0% experiment without an invented produced molecule.", "evidence_ids": [P[3], T[3]]},
            {"region": "Equation 1, upper right", "finding": "Methanesulfonyl-protected substrate experiment recorded separately; printed 1n label collision retained.", "evidence_ids": [P[3], T[3]]},
            {"region": "Table 3, lower right including footnote", "finding": "Nine allenoate-scope rows 3ab–3aj recorded independently.", "evidence_ids": [P[3], T[3]]},
            {"region": "prose above tables", "finding": "Scope discussion, 3na failure, and X-ray structure statement linked to records/context rather than duplicated.", "evidence_ids": [P[3], T[3]]},
        ], "record_ids": page3_records, "unresolved": ["Several complex substituent structures remain source-labelled without candidate SMILES pending expert transcription."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "Scheme 2, lower left", "finding": "Gram-scale reaction, hydrolysis, coupling, reduction, and Suzuki arrows recorded as separate experiments; two stage yields remain attached to their own start/end points.", "evidence_ids": [P[4], T[4]]},
            {"region": "deuterium-experiment prose, upper columns", "finding": "Real D2O mechanistic study recorded separately with 20 equiv D2O; incorporation extent/method left missing.", "evidence_ids": [P[4], T[4]]},
            {"region": "Scheme 3, lower right", "finding": "One proposed_mechanism record covers A–J; proposed arrows are not counted as experiments.", "evidence_ids": [P[4], T[4]]},
        ], "record_ids": page4_records, "unresolved": ["Scheme-2 Suzuki substrate is 3fa in the drawing but 3ha in prose.", "Scheme 3 swaps printed 1a/2a labels relative to Tables 1–3.", "Most Scheme-2 reagent amounts/times and yield measurement bases are unreported."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "Experimental Section, upper left", "finding": "General Procedure linked to matching successful scope experiments; reaction, solvent removal, and flash chromatography represented in order.", "evidence_ids": [P[5], T[5]]},
            {"region": "General Procedure quantities", "finding": "Enamine 1/allenoate 2 charges are reversed relative to Tables 2–3; both source alternatives are retained as conflicts.", "evidence_ids": [P[5], T[5], P[3], T[3]]},
            {"region": "Acknowledgements and references", "finding": "No new experiment; cited bibliography not expanded into invented experiments.", "evidence_ids": [P[5], T[5]]},
        ], "record_ids": page5_records, "unresolved": ["No SI experimental or characterization detail was supplied."]},
    ],
}

raw_inventory = {
    "document_id": DOC,
    "title": "De Novo Construction of Substituted Terephthalates via Phosphine Catalyzed Domino Benzannulation Reactions",
    "doi": "10.1002/adsc.202001513",
    "pages_reviewed": [1, 2, 3, 4, 5],
    "visual_review": True,
    "si_association": "unestablished; no SI file supplied",
    "record_inventory": {
        "literature_examples": 3, "generic_schemes": 1, "optimization_rows": 18,
        "table2_scope_rows": 14, "table3_scope_rows": 9, "equation1_experiments": 1,
        "scheme2_concrete_experiments": 5, "deuterium_experiments": 1,
        "proposed_mechanisms": 1,
    },
    "source_conflicts": [
        "Tables 2–3: enamine 1 0.1 mmol/allenoate 2 0.12 mmol; General Procedure: enamine 1 0.12 mmol/allenoate 2 0.1 mmol.",
        "Eq 1 reuses printed label 1n for a methanesulfonyl substrate; Table 2 uses 1n for the p-OMe entry.",
        "Scheme 2 Suzuki substrate is 3fa in the drawing and 3ha in adjacent prose.",
        "Scheme 3 labels allenoate as 1a and enamine as 2a, opposite to Tables 1–3 and the Experimental Section.",
    ],
}

errors = st.validate(bundle, RUN)
if errors:
    raise SystemExit("\n".join(errors))
st.write(RUN / "submitted.json", bundle)
st.write(RUN / "raw_model/page_extraction.json", raw_inventory)
st.write(RUN / "raw_model/stereo_checks.json", stereo_checks)
st.write(RUN / "coverage/pages_0001_0005.json", coverage)
mentions = st.source_mentions(json.loads((RUN / "dataset.json").read_text()),
                              ["1a", "2a", "3aa", "3ha", "3fa", "3na"], 180)
st.write(RUN / "raw_model/source_mentions.json", mentions)
tasks = json.loads((RUN / "tasks.json").read_text())
for task in tasks:
    if task.get("document_id") == DOC and 1 <= int(task["evidence_id"].split(":page")[1].split(":")[0]) <= 5:
        task["state"] = "attempted"
        task["reason"] = "Original page image and extracted text were inspected; page-bounded findings are in coverage/pages_0001_0005.json."
st.write(RUN / "tasks.json", tasks)
import reaction_tools
st.write(RUN / "audit.json", reaction_tools.dataset_audit(json.loads((RUN / "dataset.json").read_text())))
print(json.dumps({"molecules": len(bundle["molecules"]),
                  "reactions": len(bundle["reactions"]),
                  "stereo_checks": len(stereo_checks),
                  "validation_errors": len(errors)}, ensure_ascii=False, indent=2))
