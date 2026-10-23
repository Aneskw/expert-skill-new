#!/usr/bin/env python3
"""Build source-derived v3.1 records for original PDF pages 1--6.

This is a session-local extraction builder, not a modification of the bundled
validation or scientific-quality tools.  All graph assignments remain machine
candidates and all imported reactions remain needs_review.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st


DOC = "doc_2992316b1e64fc22"
RUN = ROOT / "runs/current"


def ev(*pages: int) -> list[str]:
    result = []
    for page in pages:
        if page <= 6:
            result.append(f"{DOC}:page{page}:image")
        result.append(f"{DOC}:page{page}:text")
    return result


def txt(page: int) -> str:
    return f"{DOC}:page{page}:text"


def img(page: int) -> str:
    return f"{DOC}:page{page}:image"


def f(raw=None, normalized=None, evidence=(), method="original_extraction", derivation=None,
      status=None):
    if normalized is None and raw is not None and method == "original_extraction":
        normalized = raw
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation, status=status)


def missing(evidence=()):
    return st.fact(evidence=evidence)


def inferred(normalized, evidence, derivation, raw=None):
    return st.fact(raw, normalized, evidence, method="inference", derivation=derivation)


def derived(normalized, evidence, derivation, raw=None):
    return st.fact(raw, normalized, evidence, method="structure_derivation", derivation=derivation)


def conflict(alternatives, evidence):
    # alternatives: (raw, normalized, evidence_ids[, method, derivation])
    def alt(row):
        method = row[3] if len(row) > 3 else "original_extraction"
        derivation = row[4] if len(row) > 4 else None
        return f(row[0], row[1], row[2], method=method, derivation=derivation)
    return st.fact(
        raw=[row[0] for row in alternatives],
        normalized=None,
        evidence=evidence,
        status="conflict",
        alternatives=[alt(row) for row in alternatives],
    )


def q(number, unit, **extra):
    return {"number": number, "unit": unit, **extra}


def qr(minimum, maximum, unit, raw_text):
    return {"minimum": minimum, "maximum": maximum, "unit": unit, "raw_text": raw_text}


def set_stereo(mol, evidence, *, chirality, geometry, representation,
               atom_centers=None, double_bonds=None, axial=None, mixture=None,
               assignment=None, smiles=None, encode_smiles=True):
    basis = "Status assigned from the source name/drawing and the candidate molecular graph; it is not source verification."
    mol["stereochemistry"]["chirality_status"] = inferred(chirality, evidence, basis, chirality)
    mol["stereochemistry"]["geometry_status"] = inferred(geometry, evidence, basis, geometry)
    mol["stereochemistry"]["representation_status"] = inferred(representation, evidence, basis, representation)
    mol["stereochemistry"]["atom_centers"] = f(atom_centers, atom_centers, evidence) if atom_centers is not None else missing(evidence)
    mol["stereochemistry"]["double_bonds"] = f(double_bonds, double_bonds, evidence) if double_bonds is not None else missing(evidence)
    mol["stereochemistry"]["axial"] = inferred(axial, evidence, "Allene axial status follows the displayed terminal-CH2 allene substitution pattern.", axial) if axial is not None else missing(evidence)
    mol["stereochemistry"]["mixture"] = f(mixture, mixture, evidence) if mixture is not None else missing(evidence)
    mol["stereochemistry"]["assignment"] = f(assignment, assignment, evidence) if assignment is not None else missing(evidence)
    if smiles and encode_smiles:
        diag = st.stereo_inspect(smiles)
        mol["stereochemistry"]["isomeric_smiles"] = derived(
            diag["canonical_isomeric_smiles"], evidence,
            "RDKit canonical isomeric SMILES generated from the source-derived candidate graph.", smiles)
        mol["stereochemistry"]["identity_key"] = derived(
            diag["stereo_graph_key"], evidence,
            "RDKit stereo graph key generated from the candidate graph; this is not a sample identity.")
    else:
        mol["stereochemistry"]["isomeric_smiles"] = missing(evidence)
        mol["stereochemistry"]["identity_key"] = missing(evidence)


def molecule(mid, label, name, evidence, smiles=None, *, chirality="unresolved",
             geometry="unresolved", representation=None, atom_centers=None,
             double_bonds=None, axial=None, mixture=None, assignment=None,
             encode_smiles=True, note=None):
    mol = st.new_molecule(mid, DOC, evidence)
    mol["identity_status"] = "source_identified"
    mol["occurrences"][0]["compound_label"] = f(label, label, evidence)
    mol["names"]["original"] = f(name, name, evidence)
    mol["names"]["english"] = f(name, name, evidence)
    mol["names"]["abbreviations"] = missing(evidence)
    mol["names"]["aliases"] = missing(evidence)
    mol["cas_numbers"] = missing(evidence)
    if smiles:
        mol["representation_type"] = "smiles"
        mol["smiles"] = derived(
            smiles, evidence,
            "Candidate connectivity/stereochemistry transcribed from the source drawing and compound name; parsing does not establish source correctness.",
            {"compound_label": label, "candidate_smiles": smiles})
        representation = representation or "encoded"
    else:
        mol["representation_type"] = "unresolved"
        mol["smiles"] = missing(evidence)
        representation = representation or "not_encoded"
    set_stereo(mol, evidence, chirality=chirality, geometry=geometry,
               representation=representation, atom_centers=atom_centers,
               double_bonds=double_bonds, axial=axial, mixture=mixture,
               assignment=assignment, smiles=smiles, encode_smiles=encode_smiles)
    mol["identity_note"] = f(note, note, evidence) if note else missing(evidence)
    return mol


def named(mid, name, evidence, note=None):
    return molecule(mid, name, name, evidence, note=note,
                    chirality="unresolved", geometry="unresolved", representation="not_encoded")


def charge_fact(raw, normalized, evidence):
    return f(raw, normalized, evidence)


def participant(pid, mid, role, evidence, *, stage="reaction", mass=None,
                amount=None, volume=None, equivalents=None, mol_percent=None,
                reference=None, product_role=None):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    p["product_role"] = f(product_role, product_role, evidence) if product_role else missing(evidence)
    for key, value in (("mass", mass), ("amount", amount), ("volume", volume),
                       ("equivalents", equivalents), ("mol_percent", mol_percent)):
        p["charges"][key] = missing(evidence) if value is None else charge_fact(value[0], value[1], evidence)
    p["equivalence_reference"] = f(reference, reference, evidence) if reference else missing(evidence)
    return p


def step(step_id, evidence, participant_ids, operation, *, temperature=None,
         time=None, atmosphere=None, pressure=None, stage="reaction"):
    s = st.new_step(step_id, evidence, stage=stage)
    s["participant_ids"] = participant_ids
    s["operation"] = f(operation, operation, evidence)
    s["temperature"] = missing(evidence) if temperature is None else (
        temperature if isinstance(temperature, dict) and "normalized_value" in temperature
        else f(temperature[0], temperature[1], evidence))
    s["time"] = missing(evidence) if time is None else (
        time if isinstance(time, dict) and "normalized_value" in time
        else f(time[0], time[1], evidence))
    s["atmosphere"] = missing(evidence) if atmosphere is None else f(atmosphere, atmosphere, evidence)
    s["pressure"] = missing(evidence) if pressure is None else f(pressure[0], pressure[1], evidence)
    return s


def outcome(oid, kind, targets, evidence, raw, normalized, *, product_role=None,
            yield_type=None, basis=None, scope=None):
    o = st.new_outcome(oid, kind, targets, evidence)
    o["product_role"] = f(product_role, product_role, evidence) if product_role else missing(evidence)
    o["value"] = raw if isinstance(raw, dict) and "normalized_value" in raw else f(raw, normalized, evidence)
    o["yield_type"] = f(yield_type, yield_type, evidence) if yield_type else missing(evidence)
    o["measurement_basis"] = f(basis, basis, evidence) if basis else missing(evidence)
    o["scope"] = f(scope, scope, evidence) if scope else missing(evidence)
    return o


def changes_and_labels(specs, evidence):
    facts = []
    labels = {}
    for index, spec in enumerate(specs, 1):
        change = {
            "change_id": f"b{index}",
            "atoms": [{"atom_map": None, "element": spec[0]},
                      {"atom_map": None, "element": spec[1]}],
            "before_order": spec[2], "after_order": spec[3],
            "change_type": spec[4],
        }
        if len(spec) > 5:
            change["hydrogen_count"] = spec[5]
        facts.append(derived(change, evidence,
                             "Partial net bond change inferred from the before/after source drawings; atom mapping is unavailable."))
        if spec[4] in ("formed", "hydrogen_gain", "order_increased"):
            label = st.bond_label(spec[0], spec[1], spec[3])
            labels.setdefault(label["code"], dict(label, supporting_change_ids=[]))["supporting_change_ids"].append(change["change_id"])
    return facts, list(labels.values())


def reaction(rid, evidence, *, record_type="specific_experiment", categories=None,
             location, reaction_name, transformation, participants, steps,
             outcomes=None, change_specs=None, procedure=None, notes=None,
             unresolved=None, solvent=None, scale=None):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories or ["other"], categories or ["other"], evidence)
    r["location"] = f(location, location, evidence)
    r["reaction_name"] = f(reaction_name, reaction_name, evidence)
    r["transformation_type"] = inferred(transformation, evidence,
        "Transformation type inferred from the source before/after structures and reagents.", transformation)
    r["participants"] = participants
    r["steps"] = steps
    r["outcomes"] = outcomes or []
    r["procedure"]["original_text"] = f(procedure, procedure, evidence) if procedure else missing(evidence)
    r["procedure"]["general_procedure_reference"] = missing(evidence)
    r["procedure"]["modifications"] = missing(evidence)
    r["procedure"]["addition_order"] = missing(evidence)
    r["procedure"]["pretreatment"] = missing(evidence)
    if solvent:
        r["solvent_conditions"]["composition"] = f(solvent, solvent, evidence)
    if scale:
        r["reaction_scale"] = f(scale[0], scale[1], evidence)
    r["notes"] = [f(note, note, evidence) for note in (notes or [])]
    r["review"] = {
        "status": "needs_review", "method": "model", "reviewer": None,
        "checked_fields": [], "unresolved": unresolved or [
            "Machine transcription requires expert source and structure review."
        ],
    }
    r["completeness"] = "partial"
    if change_specs is None:
        r["classification"]["coverage"] = inferred(
            "unresolved", evidence, "No atom mapping or defensible bond-change assignment was made.")
        r["classification"]["labels"] = missing(evidence)
    else:
        bond_facts, labels = changes_and_labels(change_specs, evidence)
        r["classification"]["bond_changes"] = bond_facts
        r["classification"]["labels"] = derived(
            labels, evidence, "Taxonomy labels calculated from the partial inferred bond-change list.")
        r["classification"]["coverage"] = inferred(
            "partial", evidence,
            "Source drawings support only a partial net transformation; complete atom mapping and coproduct accounting are unavailable.")
        r["classification"]["limitations"] = [
            "No atom-mapped reaction SMILES was asserted.",
            "Implicit-H changes are net changes; proton sources and complete atom correspondence are unresolved.",
            "Unchanged bonds and bond-order decreases are not counted as newly formed sigma bonds.",
        ]
    return r


def set_conflict_step(s, field_name, alternatives, all_evidence):
    s[field_name] = conflict(alternatives, all_evidence)
    return s


def build():
    dataset = json.loads((RUN / "dataset.json").read_text())
    assert dataset["documents"][0]["document_id"] == DOC
    assert not dataset.get("experiments"), "This increment expects the mechanically prepared zero-record dataset."
    dataset["documents"][0]["title"] = "Total Synthesis of Natural Lembehyne C and Investigation of Its Cytotoxic Properties"
    st.write(RUN / "dataset.json", dataset)
    bundle = st.empty(dataset["documents"])

    z1 = [{"double_bond": "C3=C4", "configuration": "Z"}]
    z5 = [{"double_bond": "C5=C6", "configuration": "Z"}]
    z6 = [{"double_bond": "C6=C7", "configuration": "Z"}]
    z4 = [{"double_bond": "C4=C5", "configuration": "Z"}]
    ztri_17 = [
        {"double_bond": "C11=C12", "configuration": "Z"},
        {"double_bond": "C15=C16", "configuration": "Z"},
        {"double_bond": "C19=C20", "configuration": "Z"},
    ]
    ztri_lem = [
        {"double_bond": "C13=C14", "configuration": "Z"},
        {"double_bond": "C17=C18", "configuration": "Z"},
        {"double_bond": "C21=C22", "configuration": "Z"},
    ]

    # Candidate source-derived graphs.  Long chains are generated explicitly to
    # reduce transcription/counting mistakes while retaining the source labels.
    smi = {}
    smi["2"] = "Br" + "C" * 16
    smi["3"] = "C#C" + "C" * 16
    smi["5"] = "OCCC#C" + "C" * 16
    smi["6"] = "OCC/C=C\\" + "C" * 16
    smi["7"] = "CS(=O)(=O)OCC/C=C\\" + "C" * 16
    smi["8"] = "BrCC/C=C\\" + "C" * 16
    smi["9"] = "C#CCC/C=C\\" + "C" * 16
    smi["10"] = "C=C=CCC/C=C\\" + "C" * 16
    smi["11"] = "OCCCC#C"
    smi["12"] = "OCCCC#C" + "C" * 16
    smi["13"] = "OCCC/C=C\\" + "C" * 16
    smi["14"] = "O=CCC/C=C\\" + "C" * 16
    smi["15a"] = "O1CCOC1" + "C" * 10 + "=C=C"
    smi["15b"] = "O1CCCCC1O" + "C" * 11 + "=C=C"
    tri_tail = "/C=C\\CC/C=C\\CC/C=C\\" + "C" * 16
    smi["17"] = "O=C" + "C" * 9 + tri_tail
    smi["19"] = "O" + "C" * 10 + tri_tail
    smi["18"] = "O1CCCCC1O" + "C" * 10 + tri_tail
    smi["20"] = "C[Si](C)(C)C#CC(O)" + "C" * 9 + tri_tail
    smi["21"] = "C#CC(O)" + "C" * 9 + tri_tail
    smi["22"] = "C#CC(=O)" + "C" * 9 + tri_tail
    smi["1"] = smi["21"]

    mols = []
    add = mols.append
    add(molecule(f"{DOC}:mol_2", "2", "1-bromohexadecane", ev(2), smi["2"], chirality="achiral", geometry="not_applicable"))
    add(molecule(f"{DOC}:mol_3", "3", "octadec-1-yne", ev(2), smi["3"], chirality="achiral", geometry="not_applicable"))
    add(molecule(f"{DOC}:mol_4", "4", "octadec-1-ynylmagnesium bromide", ev(2), note="Transient organomagnesium structure is shown; complete ionic representation was not encoded."))
    add(molecule(f"{DOC}:mol_5", "5", "icos-3-yn-1-ol", ev(2, 3, 9), smi["5"], chirality="achiral", geometry="not_applicable"))
    add(molecule(f"{DOC}:mol_6", "6", "(Z)-icos-3-en-1-ol", ev(2, 3, 9), smi["6"], chirality="achiral", geometry="specified", double_bonds=z1))
    add(molecule(f"{DOC}:mol_7", "7", "(Z)-icos-3-en-1-yl methanesulfonate", ev(2, 3, 9), smi["7"], chirality="achiral", geometry="specified", double_bonds=z1))
    add(molecule(f"{DOC}:mol_8", "8", "(Z)-1-bromoicos-3-ene", ev(2, 3, 9), smi["8"], chirality="achiral", geometry="specified", double_bonds=z1))
    add(molecule(f"{DOC}:mol_9", "9", "(5Z)-docos-5-en-1-yne", ev(2, 3, 9), smi["9"], chirality="achiral", geometry="specified", double_bonds=z5))
    add(molecule(f"{DOC}:mol_10", "10", "(6Z)-tricosa-1,2,6-triene", ev(2, 3, 9), smi["10"], chirality="achiral", geometry="specified", representation="partial", double_bonds=z6,
                 axial=[{"fragment": "C1=C2=C3", "status": "not axially stereogenic", "basis": "terminal CH2 has two H substituents"}]))
    add(molecule(f"{DOC}:mol_11", "11", "pent-4-yn-1-ol", ev(2, 3), smi["11"], chirality="achiral", geometry="not_applicable"))
    add(molecule(f"{DOC}:mol_12", "12", "heneicos-4-yn-1-ol", ev(2, 3), smi["12"], chirality="achiral", geometry="not_applicable"))
    add(molecule(f"{DOC}:mol_13", "13", "(4Z)-heneicos-4-en-1-ol", ev(2, 3), smi["13"], chirality="achiral", geometry="specified", double_bonds=z4))
    add(molecule(f"{DOC}:mol_14", "14", "(4Z)-heneicos-4-enal", ev(2, 3), smi["14"], chirality="achiral", geometry="specified", double_bonds=z4))
    add(molecule(f"{DOC}:mol_15a", "15a", "2-(dodeca-10,11-dien-1-yl)-1,3-dioxolane", ev(1, 2, 3, 9), smi["15a"], chirality="achiral", geometry="not_applicable", representation="partial",
                 axial=[{"fragment": "C10=C11=C12", "status": "not axially stereogenic", "basis": "terminal CH2 has two H substituents"}]))
    add(molecule(f"{DOC}:mol_15b", "15b", "11,12-tridecadien-1-ol tetrahydropyranyl ether", ev(1, 2, 3, 9), smi["15b"], chirality="unspecified", geometry="not_applicable", representation="partial",
                 atom_centers=[{"center": "THP anomeric carbon", "configuration": "unreported"}],
                 axial=[{"fragment": "terminal allene", "status": "not axially stereogenic", "basis": "terminal CH2 has two H substituents"}], encode_smiles=True))
    add(molecule(f"{DOC}:mol_16a", "16a", "magnesacyclopentane intermediate 16a", ev(3), note="Intermediate shown in brackets; substituent-resolved graph not encoded."))
    add(molecule(f"{DOC}:mol_16b", "16b", "magnesacyclopentane intermediate 16b", ev(3), note="Intermediate shown in brackets; substituent-resolved graph not encoded."))
    add(molecule(f"{DOC}:mol_17", "17", "(11Z,15Z,19Z)-hexatriaconta-11,15,19-trienal", ev(2, 3, 9, 10), smi["17"], chirality="achiral", geometry="specified", double_bonds=ztri_17))
    add(molecule(f"{DOC}:mol_18", "18", "2-[(11Z,15Z,19Z)-hexatriaconta-11,15,19-trien-1-yloxy]tetrahydro-2H-pyran", ev(3, 10), smi["18"], chirality="unspecified", geometry="specified", representation="partial", atom_centers=[{"center": "THP anomeric carbon", "configuration": "unreported"}], double_bonds=ztri_17))
    add(molecule(f"{DOC}:mol_19", "19", "(11Z,15Z,19Z)-hexatriaconta-11,15,19-trien-1-ol", ev(3, 10), smi["19"], chirality="achiral", geometry="specified", double_bonds=ztri_17))
    add(molecule(f"{DOC}:mol_20", "20", "trimethylsilyl-protected lembehyne C", ev(2, 3, 10), smi["20"], chirality="unspecified", geometry="specified", representation="partial", atom_centers=[{"center": "propargylic alcohol carbon", "configuration": "unreported"}], double_bonds=ztri_lem))
    add(molecule(f"{DOC}:mol_21", "21", "rac-lembehyne C", ev(2, 3, 4, 10), smi["21"], chirality="racemic", geometry="specified", representation="partial", atom_centers=[{"center": "propargylic alcohol carbon", "configuration": "racemic; individual R/S not assigned in source"}], double_bonds=ztri_lem, mixture={"type": "racemate", "ratio": "unreported"}, encode_smiles=False))
    add(molecule(f"{DOC}:mol_22", "22", "(13Z,17Z,21Z)-octatriaconta-13,17,21-trien-1-yn-3-one", ev(2, 3, 4, 10), smi["22"], chirality="achiral", geometry="specified", double_bonds=ztri_lem))
    add(molecule(f"{DOC}:mol_1", "1", "lembehyne C", ev(1, 2, 3, 4, 5, 6, 10), smi["1"], chirality="specified", geometry="specified", representation="partial",
                 atom_centers=[{"center": "propargylic alcohol carbon", "configuration": "single natural enantiomer depicted; R/S not derived"}], double_bonds=ztri_lem,
                 assignment={"source": "wedge drawing and description as enantiomerically pure natural lembehyne C", "ee": "96%"}, encode_smiles=False))

    reagents = {
        "li_acetylide_eda": "lithium acetylide ethylenediamine complex", "dmso": "dimethyl sulfoxide",
        "etmgbr": "ethylmagnesium bromide", "thf": "tetrahydrofuran", "ethylene_oxide": "ethylene oxide",
        "acid": "acid (H+; identity unreported in scheme)", "h2": "hydrogen", "ni_oac2": "nickel(II) acetate tetrahydrate",
        "etoh": "ethanol", "mscl": "methanesulfonyl chloride", "dcm": "dichloromethane", "libr": "lithium bromide",
        "acetone": "acetone", "paraformaldehyde": "paraformaldehyde", "dihexylamine": "dihexylamine",
        "cui": "copper(I) iodide", "dioxane": "1,4-dioxane", "nbuli": "n-butyllithium",
        "hmpa": "hexamethylphosphoramide", "lindlar": "Lindlar catalyst", "quinoline": "quinoline",
        "dmp": "Dess-Martin periodinane", "ob_reagent": "dimethyl (1-diazo-2-oxopropyl)phosphonate",
        "k2co3": "potassium carbonate", "meoh": "methanol", "cp2ticl2": "titanocene dichloride",
        "mg": "magnesium", "ether": "diethyl ether", "ptsoh": "p-toluenesulfonic acid",
        "chcl3": "chloroform", "tms_acetylene": "trimethylsilylacetylene",
        "li_tms_acetylide": "lithium trimethylsilylacetylide", "tbaf": "tetrabutylammonium fluoride",
        "cbs_r": "(R)-2-methyl-CBS-oxazaborolidine",
    }
    reagent_mid = {}
    for key, name in reagents.items():
        mid = f"{DOC}:reagent_{key}"
        reagent_mid[key] = mid
        evidence = ev(2, 3)
        add(named(mid, name, evidence, note="Named reaction component; structure/CAS was not supplied by the article and was not externally completed."))
    bundle["molecules"] = mols

    M = lambda label: f"{DOC}:mol_{label}"
    R = lambda key: reagent_mid[key]
    reactions = []

    # Scheme 1 is a retrosynthetic plan, not a performed experiment.
    reactions.append(reaction(
        "scheme1_retrosynthetic_analysis", ev(2), record_type="generic_scheme",
        location="PDF page 2, Scheme 1", reaction_name="Retrosynthetic Analysis of Lembehyne C",
        transformation="retrosynthetic plan", participants=[], steps=[], outcomes=[],
        procedure="Scheme 1 presents planned disconnections and forward transformations without a specific run.",
        notes=["The arrows are planning relationships, not independent performed experiments; no scheme yield is assigned."],
        unresolved=["Retrosynthetic structures are not promoted to experimental records."], change_specs=None))

    # Scheme 2.
    e2 = ev(2, 3)
    reactions.append(reaction(
        "scheme2_2_to_3", e2, location="PDF page 2, Scheme 2, 2→3 arrow",
        reaction_name="Alkylation with lithium acetylide ethylenediamine complex", transformation="C-C bond formation by acetylide substitution",
        participants=[participant("s2", M("2"), "reactant", e2), participant("r_liace", R("li_acetylide_eda"), "reagent", e2), participant("sol_dmso", R("dmso"), "solvent", e2), participant("p3", M("3"), "product", e2, product_role="main product")],
        steps=[step("reaction_rt_overnight", e2, ["s2", "r_liace", "sol_dmso"], "stir in DMSO", temperature=("rt", "rt"), time=("overnight", "overnight"))],
        outcomes=[outcome("yield_3", "yield", [M("3")], e2, "93%", q(93, "%"), product_role="main product")],
        change_specs=[("C", "C", 0, 1, "formed"), ("C", "Br", 1, 0, "broken")], solvent="DMSO",
        procedure="Li·(H2NCH2)2, DMSO, rt, overnight, 93%.",
        unresolved=["Yield measurement basis is unreported.", "Citation 11 accompanies the transformation; whether this exact run was newly performed is not explicit on pages 1–6."]))

    reactions.append(reaction(
        "scheme2_3_to_4", e2, location="PDF page 2, Scheme 2, 3→4 arrow",
        reaction_name="Formation of magnesium acetylide 4", transformation="terminal alkyne deprotonation / C-Mg formation",
        participants=[participant("s3", M("3"), "reactant", e2), participant("r_etmgbr", R("etmgbr"), "reagent", e2), participant("sol_thf", R("thf"), "solvent", e2), participant("p4", M("4"), "product", e2, product_role="organomagnesium intermediate")],
        steps=[step("reflux_3h", e2, ["s3", "r_etmgbr", "sol_thf"], "heat under reflux", temperature=("reflux", "reflux"), time=("3 h", q(3, "h")))],
        outcomes=[outcome("yield_4", "yield", [M("4")], e2, "87%", q(87, "%"), product_role="organomagnesium intermediate")],
        change_specs=[("C", "Mg", 0, 1, "formed"), ("C", "H", 1, 0, "hydrogen_loss", 1)], solvent="THF",
        procedure="EtMgBr, THF, reflux, 3 h; 87% is printed under the arrow to 4.",
        unresolved=["It is unclear whether the printed 87% is an isolated yield of transient 4 or an arrow-associated value; measurement basis is unreported."]))

    yield5 = conflict([
        ("78% (Scheme 2)", q(78, "%"), [img(2)]),
        ("16.0 g, 68% (Experimental Section)", q(68, "%"), [txt(9)]),
    ], [img(2), txt(9)])
    reactions.append(reaction(
        "scheme2_4_to_5", [img(2), txt(3), txt(9)], location="PDF page 2, Scheme 2, 4→5 arrow; page 9 extracted experimental text",
        reaction_name="Ethylene oxide homologation and hydrolysis", transformation="acetylide addition to ethylene oxide followed by protonation",
        participants=[participant("s4", M("4"), "reactant", e2), participant("r_ethylene_oxide", R("ethylene_oxide"), "reagent", e2), participant("r_acid", R("acid"), "workup_agent", e2, stage="workup"), participant("sol_thf", R("thf"), "solvent", e2), participant("p5", M("5"), "product", e2, product_role="main product")],
        steps=[step("ethylene_oxide_rt_overnight", e2, ["s4", "r_ethylene_oxide", "sol_thf"], "add ethylene oxide and stir", temperature=("rt", "rt"), time=("overnight", "overnight")), step("acid_hydrolysis", e2, ["r_acid"], "acid hydrolysis", stage="workup")],
        outcomes=[outcome("yield_5", "yield", [M("5")], [img(2), txt(9)], yield5, None, product_role="main product")],
        change_specs=[("C", "C", 0, 1, "formed"), ("O", "H", 0, 1, "formed")], solvent="THF",
        procedure="1. ethylene oxide, THF, rt, overnight; 2. H+.",
        notes=["The later experimental text reports the one-pot preparation from octadec-1-yne/EtMgBr and gives 68%, whereas Scheme 2 prints 78% for the 4→5 arrow."],
        unresolved=["Conflicting 78% versus 68% yield retained without choosing a value.", "Yield measurement basis is unreported."]))

    time6 = conflict([
        ("rt, 24 h (Scheme 2)", q(24, "h"), [img(2)]),
        ("room temperature, 12 h (Experimental Section)", q(12, "h"), [txt(9)]),
    ], [img(2), txt(9)])
    reactions.append(reaction(
        "scheme2_5_to_6", [img(2), txt(3), txt(9)], location="PDF page 2, Scheme 2, 5→6 arrow; page 9 extracted experimental text",
        reaction_name="Z-selective partial hydrogenation of alkynol 5", transformation="alkyne semihydrogenation",
        participants=[participant("s5", M("5"), "reactant", e2), participant("r_h2", R("h2"), "reagent", e2), participant("cat_ni", R("ni_oac2"), "catalyst", e2), participant("sol_etoh", R("etoh"), "solvent", e2), participant("p6", M("6"), "product", e2, product_role="main product")],
        steps=[step("hydrogenation", [img(2), txt(9)], ["s5", "r_h2", "cat_ni", "sol_etoh"], "stir under hydrogen at room temperature", temperature=("rt", "rt"), time=time6)],
        outcomes=[outcome("yield_6", "yield", [M("6")], [img(2), txt(9)], "96%", q(96, "%"), product_role="main product")],
        change_specs=[("C", "C", 3, 2, "order_decreased"), ("C", "H", 0, 1, "hydrogen_gain", 2)], solvent="EtOH",
        procedure="H2, Ni(OAc)2·4H2O, EtOH, rt; Scheme 2 prints 24 h and 96%.",
        unresolved=["Reaction time conflicts: 24 h in Scheme 2 versus 12 h in the later experimental paragraph.", "Yield measurement basis is unreported."]))

    temp7 = conflict([
        ("rt (Scheme 2)", "rt", [img(2)]),
        ("0 °C (Experimental Section)", q(0, "°C"), [txt(9)]),
    ], [img(2), txt(9)])
    time7 = conflict([
        ("4 h (Scheme 2/prose discussion)", q(4, "h"), [img(2), txt(3)]),
        ("1 h (Experimental Section)", q(1, "h"), [txt(9)]),
    ], [img(2), txt(3), txt(9)])
    reactions.append(reaction(
        "scheme2_6_to_7", [img(2), txt(3), txt(9)], location="PDF page 2, Scheme 2, 6→7 arrow; page 9 extracted experimental text",
        reaction_name="Mesylation of alcohol 6", transformation="sulfonate ester formation",
        participants=[participant("s6", M("6"), "reactant", e2), participant("r_mscl", R("mscl"), "reagent", e2), participant("sol_dcm", R("dcm"), "solvent", e2), participant("p7", M("7"), "product", e2, product_role="main product")],
        steps=[step("mesylation", [img(2), txt(9)], ["s6", "r_mscl", "sol_dcm"], "stir with methanesulfonyl chloride", temperature=temp7, time=time7)],
        outcomes=[outcome("yield_7", "yield", [M("7")], [img(2), txt(9)], "94%", q(94, "%"), product_role="main product")],
        change_specs=[("O", "S", 0, 1, "formed"), ("O", "H", 1, 0, "broken")], solvent="CH2Cl2",
        procedure="MsCl, CH2Cl2; Scheme 2 prints rt, 4 h, 94%.",
        unresolved=["Temperature/time conflict retained: rt/4 h in Scheme 2 versus 0 °C/1 h in the later experimental paragraph.", "Yield measurement basis is unreported."]))

    reactions.append(reaction(
        "scheme2_7_to_8", [img(2), txt(3), txt(9)], location="PDF page 2, Scheme 2, 7→8 arrow; page 9 extracted experimental text",
        reaction_name="Bromide substitution of mesylate 7", transformation="C-Br formation by substitution",
        participants=[participant("s7", M("7"), "reactant", e2), participant("r_libr", R("libr"), "reagent", e2), participant("sol_acetone", R("acetone"), "solvent", e2), participant("p8", M("8"), "product", e2, product_role="main product")],
        steps=[step("reflux_3h", e2, ["s7", "r_libr", "sol_acetone"], "reflux with LiBr", temperature=("reflux", "reflux"), time=("3 h", q(3, "h")))],
        outcomes=[outcome("yield_8", "yield", [M("8")], e2, "85%", q(85, "%"), product_role="main product")],
        change_specs=[("C", "Br", 0, 1, "formed"), ("C", "O", 1, 0, "broken")], solvent="acetone",
        procedure="LiBr, (CH2)2CO, reflux, 3 h, 85%.",
        unresolved=["Yield measurement basis is unreported."]))

    reactions.append(reaction(
        "scheme2_8_to_9", [img(2), txt(3), txt(9)], location="PDF page 2, Scheme 2, 8→9 arrow; page 9 extracted experimental text",
        reaction_name="Alkylation with lithium acetylide ethylenediamine complex", transformation="C-C bond formation by acetylide substitution",
        participants=[participant("s8", M("8"), "reactant", e2), participant("r_liace", R("li_acetylide_eda"), "reagent", e2), participant("sol_dmso", R("dmso"), "solvent", e2), participant("p9", M("9"), "product", e2, product_role="main product")],
        steps=[step("reaction_rt_overnight", e2, ["s8", "r_liace", "sol_dmso"], "stir in DMSO", temperature=("rt", "rt"), time=("overnight", "overnight"))],
        outcomes=[outcome("yield_9", "yield", [M("9")], e2, "67%", q(67, "%"), product_role="main product")],
        change_specs=[("C", "C", 0, 1, "formed"), ("C", "Br", 1, 0, "broken")], solvent="DMSO",
        procedure="Li·(H2NCH2)2, DMSO, rt, overnight, 67%.",
        unresolved=["Yield measurement basis is unreported."]))

    reactions.append(reaction(
        "schemes2_3_9_to_10", [img(2), txt(3), txt(9)], location="PDF page 2, Schemes 2 and 3, 9→10 arrows; page 9 extracted experimental text",
        reaction_name="Acetylene–allene rearrangement / homologation", transformation="terminal alkyne to terminal allene homologation",
        participants=[participant("s9", M("9"), "reactant", e2), participant("r_paraform", R("paraformaldehyde"), "reagent", e2), participant("r_dihexylamine", R("dihexylamine"), "reagent", e2), participant("cat_cui", R("cui"), "catalyst", e2), participant("sol_dioxane", R("dioxane"), "solvent", e2), participant("p10", M("10"), "product", e2, product_role="main product")],
        steps=[step("reflux_100c_8h", e2, ["s9", "r_paraform", "r_dihexylamine", "cat_cui", "sol_dioxane"], "heat under reflux", temperature=("100 °C", q(100, "°C")), time=("8 h", q(8, "h")))],
        outcomes=[outcome("yield_10", "yield", [M("10")], e2, "76%", q(76, "%"), product_role="main product")],
        change_specs=[("C", "C", 0, 1, "formed"), ("C", "C", 1, 2, "order_increased"), ("C", "H", 0, 1, "hydrogen_gain", 1)], solvent="1,4-dioxane",
        procedure="(CH2O)n, HN(C6H11)2, CuI, 1,4-dioxane, 100 °C, 8 h, 76%.",
        unresolved=["Partial bond classification does not assert a complete atom map or elementary Mannich/rearrangement mechanism.", "Yield measurement basis is unreported."]))

    # Scheme 3 alternative route.
    reactions.append(reaction(
        "scheme3_11_2_to_12", [img(2), txt(3)], location="PDF page 2, Scheme 3, 11+2→12 arrow",
        reaction_name="Alkylation of pent-4-yn-1-ol-derived acetylide", transformation="C-C bond formation by alkylation",
        participants=[participant("s11", M("11"), "reactant", e2), participant("s2", M("2"), "reactant", e2), participant("r_nbuli", R("nbuli"), "reagent", e2), participant("r_hmpa", R("hmpa"), "additive", e2), participant("sol_thf", R("thf"), "solvent", e2), participant("p12", M("12"), "product", e2, product_role="main product")],
        steps=[step("alkylation_minus15_4h", e2, ["s11", "s2", "r_nbuli", "r_hmpa", "sol_thf"], "successive treatment with n-BuLi and 1-bromohexadecane", temperature=("-15 °C", q(-15, "°C")), time=("4 h", q(4, "h")))],
        outcomes=[outcome("yield_12", "yield", [M("12")], e2, "88%", q(88, "%"), product_role="main product")],
        change_specs=[("C", "C", 0, 1, "formed"), ("C", "Br", 1, 0, "broken")], solvent="THF",
        procedure="n-BuLi, HMPA, THF, -15 °C, 4 h, 88%.",
        unresolved=["Yield measurement basis is unreported."]))

    reactions.append(reaction(
        "scheme3_12_to_13", [img(2), txt(3)], location="PDF page 2, Scheme 3, 12→13 arrow and page 3 prose",
        reaction_name="Z-selective partial hydrogenation of alkynol 12", transformation="alkyne semihydrogenation",
        participants=[participant("s12", M("12"), "reactant", e2), participant("r_h2", R("h2"), "reagent", e2), participant("cat_lindlar", R("lindlar"), "unknown", [img(2)]), participant("cat_ni", R("ni_oac2"), "unknown", [txt(3)]), participant("r_quinoline", R("quinoline"), "additive", [img(2)]), participant("p13", M("13"), "product", e2, product_role="main product")],
        steps=[step("hydrogenation_rt", e2, ["s12", "r_h2", "cat_lindlar", "cat_ni", "r_quinoline"], "hydrogenation; catalyst identity conflicts between scheme and prose", temperature=("rt", "rt"))],
        outcomes=[outcome("yield_13_step", "yield", [M("13")], [img(2)], "96%", q(96, "%"), product_role="main product", scope="12→13 step"), outcome("yield_13_two_step", "yield", [M("13")], [txt(3)], "85% over two steps", q(85, "%"), product_role="main product", scope="11→13 two-step overall yield")],
        change_specs=[("C", "C", 3, 2, "order_decreased"), ("C", "H", 0, 1, "hydrogen_gain", 2)],
        procedure="Scheme 3: H2, Lindlar catalyst, quinoline, rt, 96%. Prose: molecular hydrogen/Ni(OAc)2; 85% over two steps.",
        notes=["Catalyst conflict: Scheme 3 prints Lindlar catalyst/quinoline, while page 3 prose names Ni(OAc)2. Both are retained without choosing."],
        unresolved=["Catalyst identity conflict requires authorial/SI clarification.", "Individual and two-step yield measurement bases are unreported.", "Hydrogenation duration is unreported."]))

    reactions.append(reaction(
        "scheme3_13_to_14", [img(2), txt(3)], location="PDF page 2, Scheme 3, 13→14 arrow",
        reaction_name="Dess-Martin oxidation of alcohol 13", transformation="primary alcohol oxidation to aldehyde",
        participants=[participant("s13", M("13"), "reactant", e2), participant("r_dmp", R("dmp"), "reagent", e2), participant("sol_dcm", R("dcm"), "solvent", e2), participant("p14", M("14"), "product", e2, product_role="main product")],
        steps=[step("oxidation_rt_1h", e2, ["s13", "r_dmp", "sol_dcm"], "stir with Dess-Martin periodinane", temperature=("rt", "rt"), time=("1 h", q(1, "h")))],
        outcomes=[outcome("yield_14", "yield", [M("14")], e2, "95%", q(95, "%"), product_role="main product")],
        change_specs=[("C", "O", 1, 2, "order_increased"), ("C", "H", 1, 0, "hydrogen_loss", 1)], solvent="CH2Cl2",
        procedure="Dess-Martin periodinane, CH2Cl2, rt, 1 h, 95%.",
        unresolved=["Yield measurement basis is unreported."]))

    reactions.append(reaction(
        "scheme3_14_to_9", [img(2), txt(3)], location="PDF page 2, Scheme 3, 14→9 arrow",
        reaction_name="Seyferth-Gilbert-type aldehyde homologation", transformation="aldehyde to terminal alkyne homologation",
        participants=[participant("s14", M("14"), "reactant", e2), participant("r_ob", R("ob_reagent"), "reagent", e2), participant("r_k2co3", R("k2co3"), "reagent", e2), participant("sol_meoh", R("meoh"), "solvent", e2), participant("p9", M("9"), "product", e2, product_role="main product")],
        steps=[step("homologation_rt_24h", e2, ["s14", "r_ob", "r_k2co3", "sol_meoh"], "stir aldehyde with diazophosphonate and K2CO3", temperature=("rt", "rt"), time=("24 h", q(24, "h")))],
        outcomes=[outcome("yield_9", "yield", [M("9")], e2, "87%", q(87, "%"), product_role="main product")],
        change_specs=[("C", "C", 0, 3, "formed")], solvent="MeOH",
        procedure="dimethyl (1-diazo-2-oxopropyl)phosphonate, K2CO3, MeOH, rt, 24 h, 87%.",
        unresolved=["Partial net TB-C-C classification lacks complete atom mapping and coproduct accounting.", "Yield measurement basis is unreported."]))

    # Scheme 4 cross-cyclomagnesiation.  The scheme, prose and later experimental
    # paragraph disagree on loadings/time; all alternatives are retained.
    cross_ev = [img(3), txt(3), txt(9)]
    def cross_reaction(suffix, protected, product, product_label, yield_value):
        reference = M(protected)
        cp = participant("cat_cp2ticl2", R("cp2ticl2"), "catalyst", cross_ev, reference=reference)
        cp["charges"]["mol_percent"] = conflict([
            ("5 mol% (Scheme 4)", q(5, "mol%"), [img(3)]),
            ("0.1 relative to 10 of 15a/15b = 1 mol% (page 3 ratio)", q(1, "mol%"), [txt(3)], "inference", "Calculated as 0.1/10 × 100 mol% using the printed ratio and 15a/15b as denominator."),
            ("0.4 mmol relative to 4.0 mmol 15a/15b = 10 mol% (Experimental Section)", q(10, "mol%"), [txt(9)], "inference", "Calculated as 0.4/4.0 × 100 mol% from the printed charges."),
        ], cross_ev)
        et = participant("r_etmgbr", R("etmgbr"), "reagent", cross_ev, reference=reference)
        et["charges"]["equivalents"] = conflict([
            ("4 equiv (Scheme 4)", q(4, "equiv"), [img(3)]),
            ("36:10 = 3.6 equiv (page 3 ratio)", q(3.6, "equiv"), [txt(3)], "inference", "Calculated as 36/10 using 15a/15b as the one-equivalent denominator."),
            ("40.0 mmol relative to 4.0 mmol 15a/15b = 10 equiv (Experimental Section)", q(10, "equiv"), [txt(9)], "inference", "Calculated as 40.0/4.0 using 15a/15b as the one-equivalent denominator."),
        ], cross_ev)
        mg = participant("r_mg", R("mg"), "reagent", cross_ev, reference=reference)
        mg["charges"]["equivalents"] = conflict([
            ("3 equiv (Scheme 4)", q(3, "equiv"), [img(3)]),
            ("24:10 = 2.4 equiv (page 3 ratio)", q(2.4, "equiv"), [txt(3)], "inference", "Calculated as 24/10 using 15a/15b as the one-equivalent denominator."),
            ("30.0 mmol relative to 4.0 mmol 15a/15b = 7.5 equiv (Experimental Section)", q(7.5, "equiv"), [txt(9)], "inference", "Calculated as 30.0/4.0 using 15a/15b as the one-equivalent denominator."),
        ], cross_ev)
        rx_step = step("cross_cyclomagnesiation", cross_ev, ["s10", "s15", "cat_cp2ticl2", "r_etmgbr", "r_mg", "sol_ether"], "cross-cyclomagnesiation", temperature=("rt / 20–22 °C", qr(20, 22, "°C", "20–22 °C")))
        rx_step["time"] = conflict([
            ("6 h (Scheme 4 and Experimental Section)", q(6, "h"), [img(3), txt(9)]),
            ("10 h (page 3 discussion)", q(10, "h"), [txt(3)]),
        ], cross_ev)
        return reaction(
            f"scheme4_10_{protected}_to_{product}", cross_ev,
            location=f"PDF page 3, Scheme 4, 10+{protected}→{product_label} branch; page 9 extracted experimental text",
            reaction_name="Ti-catalyzed cross-cyclomagnesiation of two 1,2-dienes", transformation="allene cross-coupling followed by hydrolysis",
            participants=[participant("s10", M("10"), "reactant", cross_ev), participant("s15", M(protected), "reactant", cross_ev), cp, et, mg, participant("sol_ether", R("ether"), "solvent", cross_ev), participant("r_acid", R("acid"), "workup_agent", cross_ev, stage="workup"), participant("p_main", M(product), "product", cross_ev, product_role="main product")],
            steps=[rx_step, step("acid_hydrolysis", cross_ev, ["r_acid"], "acid hydrolysis of the in-situ magnesacyclopentane", stage="workup")],
            outcomes=[outcome(f"yield_{product}", "yield", [M(product)], [txt(3), txt(9)], f"{yield_value}%", q(yield_value, "%"), product_role="main product")],
            change_specs=[("C", "C", 0, 1, "formed"), ("C", "C", 1, 2, "order_increased"), ("C", "C", 1, 2, "order_increased"), ("C", "H", 0, 1, "hydrogen_gain", 2)], solvent="Et2O",
            procedure="Scheme 4 prints Cp2TiCl2 (5 mol%), EtMgBr (4 equiv), Mg (3 equiv), Et2O, rt, 6 h, then H+. Page 3 prose and page 9 experimental text report different ratios/loadings/time.",
            notes=["All scheme/prose/experimental loading and time alternatives are retained as conflict facts; no value was selected."],
            unresolved=["Cp2TiCl2, EtMgBr and Mg loadings conflict across Scheme 4, page 3 prose and page 9 experimental text.", "Reaction duration conflicts (6 h versus 10 h).", "Yield measurement basis is unreported.", "Complete atom mapping and exact allene carbon correspondence require expert review."])

    reactions.append(cross_reaction("15a", "15a", "17", "17", 81))
    reactions.append(cross_reaction("15b", "15b", "18", "18", 89))

    reactions.append(reaction(
        "scheme4_18_to_19", [img(3), txt(3), txt(10)], location="PDF page 3, Scheme 4, 18→19 arrow; page 10 extracted text",
        reaction_name="THP deprotection of ether 18", transformation="ether deprotection to alcohol",
        participants=[participant("s18", M("18"), "reactant", ev(3, 10)), participant("r_ptsoh", R("ptsoh"), "reagent", ev(3, 10)), participant("sol_meoh", R("meoh"), "solvent", ev(3, 10)), participant("sol_chcl3", R("chcl3"), "solvent", ev(3, 10)), participant("p19", M("19"), "product", ev(3, 10), product_role="main product")],
        steps=[step("deprotect_55c_2h", ev(3, 10), ["s18", "r_ptsoh", "sol_meoh", "sol_chcl3"], "THP deprotection", temperature=("55 °C", q(55, "°C")), time=("2 h", q(2, "h")))],
        outcomes=[], change_specs=[("O", "H", 0, 1, "formed"), ("C", "O", 1, 0, "broken")], solvent="MeOH/CHCl3",
        procedure="p-TsOH, MeOH/CHCl3, 55 °C, 2 h.",
        notes=["The page 3 prose reports approximately 64% for conversion of 18 to 17 after deprotection and subsequent oxidation; it is not assigned as the isolated yield of 19."],
        unresolved=["Yield of isolated 19 is unreported on pages 1–6.", "The later experimental text cites a known method but gives no quantities in the supplied main article."]))

    reactions.append(reaction(
        "scheme4_19_to_17", [img(3), txt(3), txt(10)], location="PDF page 3, Scheme 4, 19→17 arrow",
        reaction_name="Dess-Martin oxidation of alcohol 19", transformation="primary alcohol oxidation to aldehyde",
        participants=[participant("s19", M("19"), "reactant", ev(3, 10)), participant("r_dmp", R("dmp"), "reagent", ev(3, 10)), participant("sol_dcm", R("dcm"), "solvent", ev(3, 10)), participant("p17", M("17"), "product", ev(3, 10), product_role="main product")],
        steps=[step("oxidation_rt_1h", ev(3, 10), ["s19", "r_dmp", "sol_dcm"], "Dess-Martin oxidation", temperature=("rt", "rt"), time=("1 h", q(1, "h")))],
        outcomes=[outcome("yield_17_step", "yield", [M("17")], [img(3)], "94%", q(94, "%"), product_role="main product", scope="19→17 step"), outcome("yield_17_two_step", "yield", [M("17")], [txt(3)], "~64% from ether 18", q(64, "%", qualifier="approx", raw_text="~64%"), product_role="main product", scope="18→17 two-step overall yield")],
        change_specs=[("C", "O", 1, 2, "order_increased"), ("C", "H", 1, 0, "hydrogen_loss", 1)], solvent="CH2Cl2",
        procedure="Dess-Martin periodinane, CH2Cl2, rt, 1 h, 94%.",
        unresolved=["Measurement bases for the 94% step yield and ~64% two-step yield are unreported."]))

    # Scheme 5.
    reactions.append(reaction(
        "scheme5_17_to_20", [img(3), txt(3), txt(10)], location="PDF page 3, Scheme 5, 17→20 arrow; page 10 extracted experimental text",
        reaction_name="Addition of lithium trimethylsilylacetylide to aldehyde 17", transformation="nucleophilic alkynyl addition to aldehyde",
        participants=[participant("s17", M("17"), "reactant", ev(3, 10)), participant("r_tms_acetylene", R("tms_acetylene"), "reagent", ev(3, 10)), participant("r_nbuli", R("nbuli"), "reagent", ev(3, 10)), participant("r_li_tms", R("li_tms_acetylide"), "reagent", ev(3, 10)), participant("sol_thf", R("thf"), "solvent", ev(3, 10)), participant("p20", M("20"), "product", ev(3, 10), product_role="main product")],
        steps=[step("generate_li_acetylide", [txt(3), txt(10)], ["r_tms_acetylene", "r_nbuli", "sol_thf"], "generate lithium trimethylsilylacetylide", temperature=("-78 to -20 °C", qr(-78, -20, "°C", "-78 to -20 °C"))), step("addition_rt_3d", ev(3, 10), ["s17", "r_li_tms", "sol_thf"], "add aldehyde and warm/stir at room temperature", temperature=("rt / 20–22 °C", qr(20, 22, "°C", "20–22 °C")), time=("3 d", q(3, "d")))],
        outcomes=[outcome("yield_20", "yield", [M("20")], ev(3, 10), "91%", q(91, "%"), product_role="main product")],
        change_specs=[("C", "C", 0, 1, "formed"), ("C", "O", 2, 1, "order_decreased"), ("O", "H", 0, 1, "formed")], solvent="THF",
        procedure="Lithium trimethylsilylacetylide (generated from trimethylsilylacetylene/n-BuLi), THF, rt, 3 d, 91%.",
        unresolved=["Product 20 stereocenter configuration/enantiomer ratio is unreported.", "Yield measurement basis is unreported.", "Ramp duration from -20 °C to room temperature is unreported."]))

    time21 = conflict([
        ("2 h (Scheme 5 and Experimental Section)", q(2, "h"), [img(3), txt(10)]),
        ("4 h (page 3 discussion)", q(4, "h"), [txt(3)]),
    ], [img(3), txt(3), txt(10)])
    reactions.append(reaction(
        "scheme5_20_to_21", [img(3), txt(3), txt(10)], location="PDF page 3, Scheme 5, 20→21 arrow; page 10 extracted experimental text",
        reaction_name="Desilylation to rac-lembehyne C", transformation="terminal alkyne desilylation",
        participants=[participant("s20", M("20"), "reactant", ev(3, 10)), participant("r_tbaf", R("tbaf"), "reagent", ev(3, 10)), participant("sol_thf", R("thf"), "solvent", ev(3, 10)), participant("p21", M("21"), "product", ev(3, 10), product_role="main product")],
        steps=[step("desilylation_rt", [img(3), txt(3), txt(10)], ["s20", "r_tbaf", "sol_thf"], "stir with TBAF at room temperature", temperature=("rt", "rt"), time=time21)],
        outcomes=[outcome("yield_21", "yield", [M("21")], [img(3), txt(10)], "~99% / 99%", q(99, "%", qualifier="approx", raw_text="~99%"), product_role="main product")],
        change_specs=[("C", "Si", 1, 0, "broken"), ("C", "H", 0, 1, "hydrogen_gain", 1)], solvent="THF",
        procedure="TBAF, THF, rt; Scheme 5 prints 2 h and ~99%.",
        unresolved=["Reaction time conflict retained: 2 h in Scheme 5/later experimental paragraph versus 4 h in page 3 discussion.", "Yield measurement basis is unreported."]))

    reactions.append(reaction(
        "scheme5_21_to_22", [img(3), txt(4), txt(10)], location="PDF page 3, Scheme 5, 21→22 arrow and page 4 prose",
        reaction_name="Dess-Martin oxidation of rac-lembehyne C", transformation="secondary alcohol oxidation to ketone",
        participants=[participant("s21", M("21"), "reactant", ev(3, 4, 10)), participant("r_dmp", R("dmp"), "reagent", ev(3, 4, 10)), participant("sol_dcm", R("dcm"), "solvent", ev(3, 4, 10)), participant("p22", M("22"), "product", ev(3, 4, 10), product_role="main product")],
        steps=[step("oxidation_rt_1h", ev(3, 4, 10), ["s21", "r_dmp", "sol_dcm"], "Dess-Martin oxidation", temperature=("rt", "rt"), time=("1 h", q(1, "h")))],
        outcomes=[outcome("yield_22", "yield", [M("22")], ev(3, 4, 10), "87%", q(87, "%"), product_role="main product")],
        change_specs=[("C", "O", 1, 2, "order_increased"), ("C", "H", 1, 0, "hydrogen_loss", 1)], solvent="CH2Cl2",
        procedure="Dess-Martin periodinane, CH2Cl2, rt, 1 h, 87%.",
        unresolved=["Yield measurement basis is unreported."]))

    reactions.append(reaction(
        "scheme5_22_to_1", [img(3), txt(4), txt(10)], location="PDF page 3, Scheme 5, 22→1 arrow and page 4 prose",
        reaction_name="CBS stereoselective reduction of ketone 22", transformation="enantioselective ketone reduction",
        participants=[participant("s22", M("22"), "reactant", ev(3, 4, 10)), participant("cat_cbs", R("cbs_r"), "catalyst", ev(3, 4, 10)), participant("sol_thf", R("thf"), "solvent", ev(3, 4, 10)), participant("p1", M("1"), "product", ev(3, 4, 10), product_role="main product")],
        steps=[step("reduction_minus40_1h", ev(3, 4, 10), ["s22", "cat_cbs", "sol_thf"], "stereoselective reduction", temperature=("-40 °C", q(-40, "°C")), time=("1 h", q(1, "h")))],
        outcomes=[outcome("yield_1", "yield", [M("1")], ev(3, 4, 10), "82%", q(82, "%"), product_role="main product"), outcome("ee_1", "selectivity", [M("1")], ev(3, 4, 10), "96% ee", q(96, "%"), product_role="main product", scope={"metric": "ee", "major_enantiomer": "natural lembehyne C; R/S not derived", "order": "major/minor not reported"})],
        change_specs=[("C", "O", 2, 1, "order_decreased"), ("O", "H", 0, 1, "hydrogen_gain", 1), ("C", "H", 0, 1, "hydrogen_gain", 1)], solvent="THF",
        procedure="(R)-2-methyl-CBS-oxazaborolidine, THF, -40 °C, 1 h, 82% (96% ee).",
        unresolved=["Yield and ee measurement bases are unreported in the supplied article text.", "Absolute R/S configuration was not derived from the wedge drawing/CBS reagent.", "The detailed reduction procedure is cited to reference 14 rather than printed in the supplied article."]))

    bundle["reactions"] = reactions
    bundle["source_context"] = [
        {"context_id": f"{DOC}:ctx_title", "type": "bibliographic_metadata", "evidence_ids": ev(1),
         "content": {"title": dataset["documents"][0]["title"], "doi": "10.1021/acs.jnatprod.0c00261"}},
        {"context_id": f"{DOC}:ctx_page1_prior_work", "type": "literature_background", "evidence_ids": ev(1),
         "content": "Prior Cp2TiCl2-catalyzed lembehyne B synthesis and Figure 1 structures are background, not new current-paper reaction runs."},
        {"context_id": f"{DOC}:ctx_table1_cytotoxicity", "type": "biological_assay", "evidence_ids": ev(4, 5),
         "content": {"analyte": "lembehyne C (1)", "metric": "IC50", "unit": "nM", "values": {"Jurkat": "161 ± 32", "K562": "373 ± 57", "U937": "124 ± 19", "HL60": "93 ± 11"}, "note": "Non-synthetic cell assay retained as source context, not a chemical reaction."}},
        {"context_id": f"{DOC}:ctx_figure2_apoptosis", "type": "biological_mechanistic_assay", "evidence_ids": ev(4, 5, 6),
         "content": {"figure": 2, "system": "Jurkat cells; annexin V/7AAD", "lembehyne_C_concentrations_nM": [400, 200, 100, 50], "incubation": "24 h", "replicates": "n=5", "note": "Flow-cytometry panels and comparator lembehyne B data are not synthetic reactions."}},
        {"context_id": f"{DOC}:ctx_figure3_cell_cycle", "type": "biological_mechanistic_assay", "evidence_ids": ev(5, 6),
         "content": {"figure": 3, "system": "Jurkat cell cycle", "lembehyne_C_concentrations_nM": [400, 200, 100, 50], "incubation": "24 h", "reported_subG1": "18.1% at 200 nM in prose", "note": "Histogram/flow plots retained as assay context."}},
        {"context_id": f"{DOC}:ctx_figure4_h2ax", "type": "biological_mechanistic_assay", "evidence_ids": ev(6),
         "content": {"figure": 4, "incubation": "2 h", "conditions": ["control", "etoposide 0.2 μM", "lembehyne C 0.4 μM", "lembehyne B 2 μM"], "finding": "lembehynes B/C did not affect H2A.X phosphorylation in the authors' analysis", "note": "Mechanistic biological observation, not an isolated chemical product."}},
        {"context_id": f"{DOC}:ctx_si_status", "type": "source_limitation", "evidence_ids": [txt(11)],
         "content": "The article lists Supporting Information, but runs/current/source_mapping.json leaves association unestablished and no SI was imported."},
    ]

    # Save candidate graph diagnostics separately.  These are raw machine checks,
    # not source verification and are not used to promote review status.
    diagnostics = {}
    for label, candidate in smi.items():
        diagnostics[label] = st.stereo_inspect(candidate)
    st.write(RUN / "raw_model/pages_0001_0006_standard_stereo_checks.json", diagnostics)

    raw = {
        "batch_id": "allene-stereo-v31-200-20261017",
        "document_id": DOC,
        "visual_pages": [1, 2, 3, 4, 5, 6],
        "title": dataset["documents"][0]["title"],
        "observations": [
            {"page": 1, "regions": ["article header/abstract", "graphical abstract", "Figure 1", "introduction"], "decision": "No specific current synthetic run; title/background and overview only."},
            {"page": 2, "regions": ["Scheme 1", "Scheme 2", "Scheme 3", "Results and Discussion"], "decision": "Scheme 1 generic retrosynthesis; Scheme 2/3 arrows transcribed as specific source-reported transformations, with duplicate 9→10 merged."},
            {"page": 3, "regions": ["Scheme 4", "Scheme 5", "two-column chemistry discussion"], "decision": "Scheme branches separated; source conflicts preserved; cross-page experimental text used only as supporting context, not visual coverage."},
            {"page": 4, "regions": ["completion of Scheme 5 prose", "Table 1", "Figure 2"], "decision": "Synthetic 22→1 evidence linked; biological table/figure retained as source_context."},
            {"page": 5, "regions": ["Figure 3", "cytotoxicity/apoptosis prose"], "decision": "Non-synthetic biological assay context only."},
            {"page": 6, "regions": ["Figure 4", "apoptosis/cell-cycle/H2A.X prose"], "decision": "Non-synthetic mechanistic assay context only."},
        ],
        "reaction_count": len(reactions),
        "molecule_count": len(mols),
        "limitations": [
            "No SI was explicitly associated.",
            "Pages 9 and 10 extracted text was used to expose conflicts and complete cross-page context; their images were not claimed as reviewed.",
            "Candidate SMILES and bond classifications are inferred and remain needs_review.",
        ],
    }
    st.write(RUN / "raw_model/pages_0001_0006_raw_extraction.json", raw)

    ids_by_page = {
        1: [],
        2: [r["reaction_id"] for r in reactions if "scheme1_" in r["reaction_id"] or "scheme2_" in r["reaction_id"] or "scheme3_" in r["reaction_id"] or "schemes2_3" in r["reaction_id"]],
        3: [r["reaction_id"] for r in reactions if "scheme4_" in r["reaction_id"] or "scheme5_" in r["reaction_id"] or "schemes2_3" in r["reaction_id"]],
        4: [f"{DOC}:scheme5_21_to_22", f"{DOC}:scheme5_22_to_1"],
        5: [],
        6: [],
    }
    coverage = {
        "document_id": DOC,
        "semantic_accuracy_verified": False,
        "pages": [
            {"page": 1, "state": "attempted", "regions": [
                {"region": "title, abstract and graphical abstract", "finding": "Article identity and whole-route summary; not a specific run.", "evidence_ids": [img(1), txt(1)]},
                {"region": "Figure 1 and introduction", "finding": "Lembehynes B/C structures and cited prior lembehyne B work; retained as background, not current experiments.", "evidence_ids": [img(1), txt(1)]}],
             "record_ids": ids_by_page[1], "unresolved": ["Graphical abstract gives a 41% overall synthesis yield but does not define a single experimental run or measurement basis."]},
            {"page": 2, "state": "attempted", "regions": [
                {"region": "Scheme 1", "finding": "Retrosynthetic analysis; one generic-scheme record, not performed experiments.", "evidence_ids": [img(2), txt(2)]},
                {"region": "Scheme 2, every arrow 2→3→4→5→6→7→8→9→10", "finding": "Eight source-reported transformations; arrow-specific conditions/yields transcribed, and later text conflicts retained.", "evidence_ids": [img(2), txt(2), txt(3), txt(9)]},
                {"region": "Scheme 3, every arrow 11+2→12→13→14→9→10", "finding": "Four new alternative-route records; repeated 9→10 linked to the Scheme 2 record; catalyst discrepancy retained.", "evidence_ids": [img(2), txt(2), txt(3)]},
                {"region": "bottom Results and Discussion", "finding": "Chemistry rationale and citations; no additional table row or reaction arrow.", "evidence_ids": [img(2), txt(2)]}],
             "record_ids": ids_by_page[2], "unresolved": ["Several arrows cite prior methods; exact status of a newly executed versus literature-adapted run is not always explicit.", "Yield measurement bases are generally unreported."]},
            {"page": 3, "state": "attempted", "regions": [
                {"region": "Scheme 4, both 15a/15b branches and 18→19→17", "finding": "Four specific transformations; scheme/prose/experimental loading and time conflicts retained.", "evidence_ids": [img(3), txt(3), txt(9), txt(10)]},
                {"region": "Scheme 5, 17→20→21→22→1", "finding": "Four specific transformations; TBAF time conflict and 96% ee retained.", "evidence_ids": [img(3), txt(3), txt(4), txt(10)]},
                {"region": "two-column chemistry discussion", "finding": "Narrative duplicates of schemes plus yields/ratios; linked to records rather than duplicated.", "evidence_ids": [img(3), txt(3)]}],
             "record_ids": ids_by_page[3], "unresolved": ["Cross-cyclomagnesiation loadings differ among Scheme 4, prose ratio and later experimental text.", "Absolute R/S of natural lembehyne C was not derived.", "Several measurement bases remain unreported."]},
            {"page": 4, "state": "attempted", "regions": [
                {"region": "top-left Scheme 5 completion prose", "finding": "Duplicate evidence for 21→22 and 22→1; merged with page-3 records.", "evidence_ids": [img(4), txt(4)]},
                {"region": "Table 1", "finding": "Four-cell-line IC50 assay values; non-synthetic and retained in source_context.", "evidence_ids": [img(4), txt(4)]},
                {"region": "Figure 2 and caption", "finding": "Jurkat apoptosis flow-cytometry dose series; non-synthetic biological assay context.", "evidence_ids": [img(4), txt(4)]}],
             "record_ids": ids_by_page[4], "unresolved": ["IC50 measurement details are outside this chemistry-record increment; no chemical product was invented."]},
            {"page": 5, "state": "attempted", "regions": [
                {"region": "Figure 3 panels A–F and caption", "finding": "Jurkat cell-cycle dose series and comparison with cited lembehyne B data; retained as source_context, not synthetic records.", "evidence_ids": [img(5), txt(5)]},
                {"region": "bottom continuation prose", "finding": "IC50 recap and apoptosis discussion; no reaction arrow/table of chemical synthesis.", "evidence_ids": [img(5), txt(5)]}],
             "record_ids": ids_by_page[5], "unresolved": ["Individual plotted cell-cycle percentages are not all legible enough for exhaustive numerical transcription; figure/caption context is preserved."]},
            {"page": 6, "state": "attempted", "regions": [
                {"region": "Figure 4 panels A–D and caption", "finding": "H2A.X phosphorylation control/test assay; mechanistic biological observation retained as source_context.", "evidence_ids": [img(6), txt(6)]},
                {"region": "apoptosis and cell-cycle prose", "finding": "Reports late/early apoptosis and sub-G1 observations; no synthetic reaction arrow.", "evidence_ids": [img(6), txt(6)]},
                {"region": "DNA-damage discussion", "finding": "Author interpretation and cited biology, not a chemical synthesis experiment.", "evidence_ids": [img(6), txt(6)]}],
             "record_ids": ids_by_page[6], "unresolved": ["The text says 'all three cell lines' despite four named leukemia lines elsewhere; no correction was imposed."]},
        ],
    }
    st.write(RUN / "coverage/pages_0001_0006.json", coverage)
    st.write(RUN / "submitted.json", bundle)
    print(json.dumps({"molecules": len(mols), "reactions": len(reactions), "source_context": len(bundle["source_context"])}, ensure_ascii=False))


if __name__ == "__main__":
    build()
