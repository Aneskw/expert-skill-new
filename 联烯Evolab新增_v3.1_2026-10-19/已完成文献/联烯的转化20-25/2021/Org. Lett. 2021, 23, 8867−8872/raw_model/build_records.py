#!/usr/bin/env python3
"""Build evidence-linked v3.1 records for doc_89290e15b92920b2.

This is a source-specific extraction artifact, not a modification of the shared
tools, schemas, quality rules, or tests.
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
DATASET_PATH = RUN / "dataset.json"
DOC = "doc_89290e15b92920b2"


def ev(page: int, kind: str = "text") -> str:
    return f"{DOC}:page{page}:{kind}"


def write(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n")


def F(raw=None, normalized=None, evidence=(), method="original_extraction", derivation=None,
      status=None):
    if normalized is None and raw is not None and method in ("original_extraction", "inherited_general_condition"):
        normalized = raw
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation, status=status)


def missing(evidence=()):
    return st.fact(evidence=evidence)


MOLECULES = {}


def add_molecule(mid, label, name, evidence, smiles=None, aliases=None, stereo="auto",
                 assignment=None, double_bonds=None, axial=None, mixture=None,
                 chirality=None, geometry=None, representation=None, identity_note=None):
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, evidence)
    m["names"]["original"] = F(name, name, evidence)
    m["names"]["english"] = F(name, name, evidence)
    if aliases:
        m["names"]["aliases"] = F(aliases, aliases, evidence)
    m["cas_numbers"] = missing()
    if smiles:
        diag = st.stereo_inspect(smiles)
        canonical = diag["canonical_isomeric_smiles"]
        basis = f"Candidate connectivity transcribed from the source label/drawing for {label}; RDKit canonicalization is a graph diagnostic, not source verification."
        m["representation_type"] = "smiles"
        m["smiles"] = F({"source_label": label, "candidate_connectivity": smiles}, canonical,
                          evidence, "structure_derivation", basis)
        m["stereochemistry"]["identity_key"] = F(None, diag["stereo_graph_key"], evidence,
                                                     "structure_derivation", basis)
        m["stereochemistry"]["isomeric_smiles"] = F(None, canonical, evidence,
                                                        "structure_derivation", basis)
        if stereo == "auto":
            potential = diag["potential_stereo"]
            unspecified = [x for x in potential if x["specified"] == "Unspecified"]
            chiral_potential = [x for x in potential if x["type"].startswith("Atom_")]
            bond_potential = [x for x in potential if "Bond" in x["type"]]
            chirality = ("specified" if chiral_potential and all(x["specified"] != "Unspecified" for x in chiral_potential)
                         else "unspecified" if chiral_potential else "achiral")
            geometry = ("specified" if bond_potential and all(x["specified"] != "Unspecified" for x in bond_potential)
                        else "unspecified" if bond_potential else "not_applicable")
            representation = "partial" if unspecified else "encoded"
            m["stereochemistry"]["atom_centers"] = F(None, diag["atom_centers"], evidence,
                                                         "structure_derivation", basis)
            m["stereochemistry"]["double_bonds"] = F(None, diag["double_bonds"], evidence,
                                                          "structure_derivation", basis)
        elif stereo == "allene_axis_unassigned":
            chirality = chirality or "unspecified"
            geometry = geometry or "not_applicable"
            representation = representation or "partial"
            axial = axial or [{"axis": "C=C=C", "configuration": "unassigned",
                               "note": "Both allene termini have nonequivalent substituents; no R_a/S_a label is printed."}]
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing()
        chirality = chirality or "unresolved"
        geometry = geometry or "unresolved"
        representation = representation or "not_encoded"
    if assignment is not None:
        m["stereochemistry"]["assignment"] = F(assignment, assignment, evidence)
    if double_bonds is not None:
        m["stereochemistry"]["double_bonds"] = F(double_bonds, double_bonds, evidence)
    if axial is not None:
        m["stereochemistry"]["axial"] = F(axial, axial, evidence)
    if mixture is not None:
        m["stereochemistry"]["mixture"] = F(mixture, mixture, evidence)
    for key, value in (("chirality_status", chirality or "unresolved"),
                       ("geometry_status", geometry or "unresolved"),
                       ("representation_status", representation or "unresolved")):
        m["stereochemistry"][key] = F(None, value, evidence, "inference",
                                         f"Stereo status assigned conservatively from the source drawing/text for {label}; absence of a label is not treated as achirality.")
    if identity_note:
        m["identity_note"] = F(identity_note, identity_note, evidence)
    MOLECULES[mid] = m
    return mid


def participant(pid, mid, role, evidence, product_role=None):
    p = st.new_participant(pid, mid, role, evidence)
    if product_role:
        p["product_role"] = F(product_role, product_role, evidence)
    return p


def base_reaction(suffix, page, record_type="specific_experiment", categories=None,
                  location=None, evidence=None):
    evidence = evidence or [ev(page, "text"), ev(page, "image")]
    r = st.new_reaction(f"{DOC}:{suffix}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = F(categories or ["other"], categories or ["other"], evidence)
    r["location"] = F(location or f"PDF page {page}", {"page": page, "region": location or "source page"}, evidence)
    r["review"]["unresolved"] = [
        "No explicitly associated Supporting Information was supplied; main-text omissions remain unreported.",
        "Dense fullerene structures and atom correspondence are not fully encoded; bond classification is partial or unresolved.",
    ]
    return r


def set_names(r, reaction_name, transformation, evidence):
    r["reaction_name"] = F(reaction_name, reaction_name, evidence)
    r["transformation_type"] = F(None, transformation, evidence, "inference",
                                    "Normalized transformation type derived from the source reaction name and depicted arrow.")


def set_solvent(r, mid, raw, evidence):
    r["solvent_conditions"]["composition"] = F(raw, [mid], evidence)


def add_step(r, sid, participant_ids, evidence, operation, temperature=None, time=None,
             atmosphere=None, temp_method="original_extraction", time_method="original_extraction"):
    s = st.new_step(sid, evidence)
    s["participant_ids"] = participant_ids
    s["operation"] = F(operation, operation, evidence)
    if temperature is not None:
        raw, norm = temperature
        s["temperature"] = F(raw, norm, evidence, temp_method,
                               "Inherited from the scheme arrow applying to this substrate entry." if temp_method == "inherited_general_condition" else None)
    if time is not None:
        raw, norm = time
        s["time"] = F(raw, norm, evidence, time_method,
                        "Inherited from the scheme arrow applying to this substrate entry." if time_method == "inherited_general_condition" else None)
    if atmosphere is not None:
        s["atmosphere"] = F(atmosphere, atmosphere, evidence)
    r["steps"].append(s)


def outcome(r, oid, kind, targets, evidence, raw, normalized, product_role="product",
            yield_type=None, measurement_basis=None, scope=None):
    o = st.new_outcome(oid, kind, targets, evidence)
    o["product_role"] = F(product_role, product_role, evidence)
    o["value"] = F(raw, normalized, evidence)
    if yield_type:
        o["yield_type"] = F(yield_type, yield_type, evidence)
    if measurement_basis:
        o["measurement_basis"] = F(measurement_basis, measurement_basis, evidence)
    if scope:
        o["scope"] = F(scope, scope, evidence)
    r["outcomes"].append(o)


def classify(r, codes, evidence, coverage="partial", limitation=None):
    changes = []
    labels = []
    pairs = {"C-C": ("C", "C", 1), "DB-C-O": ("C", "O", 2), "N-H": ("N", "H", 1)}
    for i, code in enumerate(codes, 1):
        a, b, order = pairs[code]
        cid = f"b{i}"
        change = {"change_id": cid, "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
                  "before_order": 0, "after_order": order, "change_type": "formed",
                  "count_status": "class-level partial representation; exact atom mapping/count unresolved"}
        changes.append(F(None, change, evidence, "structure_derivation",
                         "Bond class inferred from the depicted reactant/product connectivity; dense fullerene atom mapping was not assigned."))
        lab = st.bond_label(a, b, order)
        lab["supporting_change_ids"] = [cid]
        labels.append(lab)
    r["classification"]["bond_changes"] = changes
    r["classification"]["labels"] = F(None, labels, evidence, "inference",
                                          "Taxonomy labels are derived from the partial class-level bond changes, not from a complete atom-mapped reaction.")
    r["classification"]["coverage"] = F(None, coverage, evidence, "inference",
                                            "Coverage is conservative because product atom correspondence is not fully encoded.")
    r["classification"]["limitations"] = [limitation or "Only bond classes visible with confidence are represented; exact atom mapping, bond count, bond-order decreases, and coproduct balance remain unresolved."]


def classify_empty(r, evidence, coverage="unresolved", limitation=None):
    r["classification"]["labels"] = F(None, [], evidence, "inference",
                                          "No bond-formation class is asserted for this record.")
    r["classification"]["coverage"] = F(None, coverage, evidence, "inference",
                                            "No verified product mapping is available for classification.")
    r["classification"]["limitations"] = [limitation or "Observed product connectivity is insufficient for a bond-change assignment."]


def add_note(r, raw, evidence, normalized=None, method="original_extraction", derivation=None):
    r["notes"].append(F(raw, normalized if normalized is not None else raw, evidence, method, derivation))


def main():
    dataset = json.loads(DATASET_PATH.read_text())
    document = dataset["documents"][0]
    document["title"] = "Isocyanide-Induced Annulation Leading to Cyclopento-, Methano-, and Cyclopentano-[60]Fullerene Derivatives"
    DATASET_PATH.write_text(json.dumps(dataset, ensure_ascii=False, indent=2) + "\n")

    e12 = [ev(1, "text"), ev(2, "text"), ev(2, "image")]
    e2 = [ev(2, "text"), ev(2, "image")]
    e23 = [ev(2, "text"), ev(3, "text"), ev(3, "image")]
    e3 = [ev(3, "text"), ev(3, "image")]
    e34 = [ev(3, "text"), ev(4, "text"), ev(4, "image")]
    e4 = [ev(4, "text"), ev(4, "image")]

    # Named/simple inputs and media. Candidate SMILES are graph diagnostics only.
    add_molecule("mol_c60", "C60", "[60]fullerene", e12, identity_note="Source identifies the C60 cage; a complete cage SMILES was not transcribed.")
    add_molecule("mol_tbunc", "tBuNC", "tert-butyl isocyanide", e12, "[C-]#[N+]C(C)(C)C", aliases=["tBuNC"])
    add_molecule("mol_cync", "CyNC", "cyclohexyl isocyanide", e2, "[C-]#[N+]C1CCCCC1", aliases=["CyNC"])
    add_molecule("mol_dmad", "DMAD", "dimethyl acetylenedicarboxylate", e12, "COC(=O)C#CC(=O)OC", aliases=["DMAD"])
    add_molecule("mol_odcb", "oDCB", "1,2-dichlorobenzene", e2, "Clc1ccccc1Cl", aliases=["o-dichlorobenzene", "oDCB"])
    add_molecule("mol_toluene", "toluene", "toluene", e2, "Cc1ccccc1")
    add_molecule("mol_tsoh_h2o", "TsOH·H2O", "p-toluenesulfonic acid monohydrate", e2, "Cc1ccc(S(=O)(=O)O)cc1.O")
    add_molecule("mol_chcl3", "CHCl3", "chloroform", e2, "ClC(Cl)Cl", aliases=["CHCl3"])
    add_molecule("mol_water", "H2O", "water", e4, "O", aliases=["H2O"])
    add_molecule("mol_air", "air", "air", [ev(1, "text")], identity_note="Named atmospheric mixture; composition was not specified and no unique concrete SMILES is asserted.")

    # Scheme 1 products; dense structures are kept source-labelled.
    add_molecule("mol_4a", "4a", "cyclopent-2-en-1-imino-[60]fullerene 4a", e2,
                 assignment={"relative_configuration": "tert-butyl group trans to fullerene cage", "basis": "single-crystal X-ray structure"},
                 chirality="partially_specified", geometry="not_applicable", representation="not_encoded")
    add_molecule("mol_4b", "4b", "cyclopent-2-en-1-imino-[60]fullerene 4b stereoisomer mixture", e2,
                 assignment={"major": "hypothesized trans to fullerene cage", "minor": "other stereoisomer; not assigned"},
                 mixture={"type": "stereoisomer mixture", "minor_13C_peak_intensity": "~1/5 of major"},
                 chirality="mixture", geometry="not_applicable", representation="not_encoded")
    for label, name in [("5a", "ketenimine methano-[60]fullerene 5a"), ("5b", "ketenimine methano-[60]fullerene 5b"),
                        ("6a", "fullerene amide 6a"), ("6b", "fullerene amide 6b")]:
        add_molecule(f"mol_{label}", label, name, e2, chirality="unresolved", geometry="unresolved", representation="not_encoded")

    # Scheme 2 alkyne scope.
    alkynes = {
        "7a": ("methyl 3-(4-nitrophenyl)propiolate", "COC(=O)C#Cc1ccc([N+](=O)[O-])cc1", e23, "auto"),
        "7b": ("dimethyl (E)-hex-2-en-4-ynedioate", "COC(=O)/C=C/C#CC(=O)OC", e3, "auto"),
        "7c": ("methyl propiolate", "COC(=O)C#C", e23, "auto"),
        "7d": ("methyl but-2-ynoate", "COC(=O)C#CC", e23, "auto"),
        "7e": ("methyl 4,4,4-trifluorobut-2-ynoate", "COC(=O)C#CC(F)(F)F", e23, "auto"),
        "7f": ("methyl 3-phenylpropiolate", "COC(=O)C#Cc1ccccc1", e23, "auto"),
        "7g": ("hexyl propiolate", "CCCCCCOC(=O)C#C", e23, "auto"),
        "7h": ("but-3-yn-2-one", "CC(=O)C#C", e23, "auto"),
        "7i": ("dimethyl hexa-2,4-diynedioate", "COC(=O)C#CC#CC(=O)OC", e3, "auto"),
        "7j": ("bis(4-nitrophenyl)acetylene", "O=[N+]([O-])c1ccc(C#Cc2ccc([N+](=O)[O-])cc2)cc1", e3, "auto"),
    }
    for label, (name, smi, evidence, stereo) in alkynes.items():
        add_molecule(f"mol_{label}", label, name, evidence, smi, stereo=stereo)
    for label, name, evidence, extra in [
        ("8a", "cyclopent-2-en-1-imino-[60]fullerene 8a", e3, {}),
        ("9a", "ketenimine methano-[60]fullerene 9a", e3, {}),
        ("9b", "ketenimine methano-[60]fullerene 9b", e3, {"geometry": "specified", "double_bonds": [{"configuration": "shown in source drawing; retained C=C bond"}]}),
        ("10a", "fullerene amide 10a", e3, {}),
        ("10b", "fullerene amide 10b", e3, {"geometry": "specified", "double_bonds": [{"configuration": "shown in source drawing; retained C=C bond"}]}),
        ("11a", "cyano fullerene byproduct 11a", e3, {}),
        ("7i_trace_product", "trace fullerene product from 7i", e3, {}),
    ]:
        add_molecule(f"mol_{label}", label, name, evidence, chirality="unresolved",
                     geometry=extra.get("geometry", "unresolved"), representation="not_encoded",
                     double_bonds=extra.get("double_bonds"))

    # Scheme 4 allene scope and stereochemically discussed products.
    add_molecule("mol_12a", "12a", "dimethyl 2,3-pentadienedioate", e34,
                 "COC(=O)C=C=CC(=O)OC", stereo="allene_axis_unassigned",
                 identity_note="Allene-axis configuration is not assigned in the source; the connectivity candidate does not encode axial stereochemistry.")
    add_molecule("mol_12b", "12b", "ethyl buta-2,3-dienoate", e34, "CCOC(=O)C=C=C",
                 identity_note="Terminal CH2 allene shown in Scheme 4; no axial configuration is applicable to this terminal allene connectivity.")
    add_molecule("mol_13a", "13a", "cyclopentano-[60]fullerene amide 13a", e34,
                 assignment={"relative_stereochemistry": ["Hb trans to Ha", "Hb cis to Hc"],
                             "basis": "coupling constants and HMBC discussed on pages 3-4"},
                 chirality="partially_specified", geometry="not_applicable", representation="not_encoded")
    add_molecule("mol_13b", "13b", "cyclopentano-[60]fullerene amide 13b", e34,
                 assignment={"relative_stereochemistry": ["Hh cis to He", "Hd cis to He", "Hd trans to Hg", "Hd cis to Hh"],
                             "basis": "long-range and vicinal coupling constants plus HMBC on pages 3-4"},
                 chirality="partially_specified", geometry="not_applicable", representation="not_encoded")

    reactions = []
    generic_id = f"{DOC}:graphical_abstract_general"
    g = base_reaction("graphical_abstract_general", 1, "generic_scheme", ["other"], "Page 1 graphical abstract", [ev(1, "image"), ev(1, "text")])
    set_names(g, "Isocyanide-induced annulation of C60 with alkynes or allenoates", "generic multicomponent fullerene annulation", g["evidence_ids"])
    add_note(g, "Graphical abstract summarizes the alkyne annulation/hydration and allenoate four-component annulation families; it is not an additional experiment.", g["evidence_ids"])
    classify_empty(g, g["evidence_ids"])
    reactions.append(g)

    lit1 = base_reaction("literature_c60_dmad", 2, "literature_example", ["other"], "Scheme 1a, products 1 and 2", e2)
    set_names(lit1, "Reported three-component reaction of C60, alkyl isocyanides, and DMAD", "prior-work fullerene cycloaddition", e2)
    add_note(lit1, "Scheme 1a attributes fullerene derivatives 1 and 2 to reference 8; this is prior work, not a current-paper experiment.", e2)
    classify_empty(lit1, e2)
    reactions.append(lit1)
    lit2 = base_reaction("literature_lu3n_c80_dmad", 2, "literature_example", ["other"], "Scheme 1a, product 3", e2)
    set_names(lit2, "Reported Lu3N@Ih-C80/tBuNC/DMAD reaction", "prior-work metallofulleroid annulation", e2)
    add_note(lit2, "Scheme 1a attributes ketenimine metallofulleroid 3 to reference 9; this is prior work.", e2)
    classify_empty(lit2, e2)
    reactions.append(lit2)

    def scheme1_ann(suffix, rnc_mid, p4, p5, y4, y5):
        r = base_reaction(suffix, 2, categories=["substrate_scope"], location="Scheme 1b", evidence=e12)
        r["general_scheme_ids"] = [generic_id]
        set_names(r, "Three-component annulation of C60, alkyl isocyanide, and DMAD", "fullerene annulation", e12)
        r["participants"] = [participant("c60", "mol_c60", "reactant", e12), participant("rnc", rnc_mid, "reactant", e12),
                             participant("dmad", "mol_dmad", "reactant", e12), participant("odcb", "mol_odcb", "solvent", e12),
                             participant("p4", p4, "product", e12, "co-product"), participant("p5", p5, "product", e12, "major product")]
        set_solvent(r, "mol_odcb", "oDCB", e12)
        add_step(r, "annulation", ["c60", "rnc", "dmad", "odcb"], e12,
                 "C60, alkyl isocyanide, and DMAD were reacted in oDCB", ("80 °C", {"number": 80, "unit": "°C"}),
                 ("1 h", {"number": 1, "unit": "h"}))
        outcome(r, "yield_p4", "yield", [p4], e2, f"{y4}%", {"number": y4, "unit": "%"}, "co-product")
        outcome(r, "yield_p5", "yield", [p5], e2, f"{y5}%", {"number": y5, "unit": "%"}, "major product")
        if p4 == "mol_4b":
            outcome(r, "stereoisomer_signal_ratio", "selectivity", [p4], e2,
                    "a set of minor peaks of ~1/5 intensity", "minor 13C NMR peak set ~1/5 intensity of major",
                    product_role="product mixture", scope="major/minor stereoisomer signal intensity; major trans assignment is an author hypothesis", measurement_basis="13C NMR peak intensity")
        classify(r, ["C-C"], e2)
        return r

    reactions += [scheme1_ann("scheme1b_tbunc", "mol_tbunc", "mol_4a", "mol_5a", 7, 36),
                  scheme1_ann("scheme1b_cync", "mol_cync", "mol_4b", "mol_5b", 6, 26)]

    def hydration(suffix, reactant_mid, product_mid, y):
        r = base_reaction(suffix, 2, categories=["derivatization"], location="Scheme 1b hydration", evidence=e12)
        r["general_scheme_ids"] = [generic_id]
        set_names(r, "Acid-promoted hydration of ketenimine methano-[60]fullerene", "ketenimine hydration to amide", e12)
        r["participants"] = [participant("ketenimine", reactant_mid, "reactant", e12),
                             participant("tsoh", "mol_tsoh_h2o", "reagent", e2), participant("chcl3", "mol_chcl3", "solvent", e2),
                             participant("amide", product_mid, "product", e12, "main product")]
        set_solvent(r, "mol_chcl3", "CHCl3", e2)
        add_step(r, "hydration", ["ketenimine", "tsoh", "chcl3"], e2, "Treatment with TsOH·H2O in CHCl3",
                 ("r.t.", "r.t."), ("20 min", {"number": 20, "unit": "min"}))
        outcome(r, "amide_yield", "yield", [product_mid], e2, f"{y}%", {"number": y, "unit": "%"}, "main product")
        classify(r, ["DB-C-O", "N-H"], e2)
        return r

    reactions += [hydration("scheme1b_hydration_5a", "mol_5a", "mol_6a", 95),
                  hydration("scheme1b_hydration_5b", "mol_5b", "mol_6b", 92)]

    tol = base_reaction("toluene_repeat", 2, categories=["control"], location="Page 2 right-column solvent-control prose", evidence=e2)
    tol["general_scheme_ids"] = [generic_id]
    set_names(tol, "Repeat of C60/tBuNC/DMAD reaction in toluene", "solvent-control fullerene annulation", e2)
    tol["participants"] = [participant("c60", "mol_c60", "reactant", e2), participant("tbunc", "mol_tbunc", "reactant", e2),
                           participant("dmad", "mol_dmad", "reactant", e2), participant("toluene", "mol_toluene", "solvent", e2),
                           participant("p4a", "mol_4a", "product", e2, "product"), participant("p5a", "mol_5a", "product", e2, "product")]
    set_solvent(tol, "mol_toluene", "toluene, the same solvent used in the original report", e2)
    add_step(tol, "repeat", ["c60", "tbunc", "dmad", "toluene"], e2, "Reaction repeated in toluene; other conditions not stated")
    outcome(tol, "products_obtained", "other", ["mol_4a", "mol_5a"], e2, "Still, the same products 4a and 5a were obtained.", "4a and 5a obtained", product_role="products")
    classify(tol, ["C-C"], e2)
    reactions.append(tol)

    air = base_reaction("air_conversion_5b", 1, categories=["derivatization"], location="Page 1 right-column prose", evidence=[ev(1, "text"), ev(1, "image")])
    air["general_scheme_ids"] = [generic_id]
    set_names(air, "Gradual conversion of 5b to 6b in air", "ambient conversion of ketenimine fullerene to amide", air["evidence_ids"])
    air["participants"] = [participant("p5b", "mol_5b", "reactant", air["evidence_ids"]), participant("air", "mol_air", "additive", air["evidence_ids"]), participant("p6b", "mol_6b", "product", air["evidence_ids"], "product")]
    add_step(air, "air_exposure", ["p5b", "air"], air["evidence_ids"], "Exposure/conversion in air", time=("gradually", "gradually"), atmosphere="air")
    outcome(air, "conversion_observation", "other", ["mol_6b"], air["evidence_ids"], "5b could gradually be converted to 6b in air", "gradual conversion to 6b")
    classify(air, ["DB-C-O", "N-H"], air["evidence_ids"], limitation="Net amide-forming bond classes are inferred from 5b/6b drawings; the oxygen/hydrogen source and complete atom balance are not reported.")
    reactions.append(air)

    def scheme2_ann(suffix, alkyne, products):
        r = base_reaction(suffix, 3, categories=["substrate_scope"], location="Scheme 2 successful substrate entry", evidence=e23)
        r["general_scheme_ids"] = [generic_id]
        set_names(r, "tBuNC-induced annulation of C60 with unsymmetric alkyne", "fullerene annulation", e23)
        r["participants"] = [participant("c60", "mol_c60", "reactant", e23), participant("tbunc", "mol_tbunc", "reactant", e23),
                             participant("alkyne", alkyne, "reactant", e23), participant("odcb", "mol_odcb", "solvent", e3)]
        for pid, mid, y, role in products:
            r["participants"].append(participant(pid, mid, "product", e3, role))
            outcome(r, f"yield_{pid}", "yield", [mid], e3, f"{y}%", {"number": y, "unit": "%"}, role)
        set_solvent(r, "mol_odcb", "oDCB", e3)
        add_step(r, "annulation", ["c60", "tbunc", "alkyne", "odcb"], e3,
                 "C60, tBuNC, and alkyne reacted in oDCB under N2", ("120 °C", {"number": 120, "unit": "°C"}),
                 ("1 h", {"number": 1, "unit": "h"}), "N2")
        classify(r, ["C-C"], e3)
        return r

    reactions.append(scheme2_ann("scheme2_7a", "mol_7a", [("p8a", "mol_8a", 4, "co-product"), ("p9a", "mol_9a", 19, "major product")]))
    reactions.append(scheme2_ann("scheme2_7b", "mol_7b", [("p9b", "mol_9b", 32, "only isolable product")]))

    for label in ("7c", "7d", "7e", "7f", "7g", "7h", "7j"):
        r = base_reaction(f"scheme2_unsuccessful_{label}", 3, categories=["substrate_scope", "control"],
                          location=f"Scheme 2 unsuccessful substrate {label}", evidence=e23 if label != "7j" else e3)
        r["general_scheme_ids"] = [generic_id]
        set_names(r, "Unsuccessful unsymmetric-alkyne substrate test", "substrate-compatibility test", r["evidence_ids"])
        r["participants"] = [participant("c60", "mol_c60", "reactant", r["evidence_ids"]), participant("tbunc", "mol_tbunc", "reactant", r["evidence_ids"]), participant("alkyne", f"mol_{label}", "reactant", r["evidence_ids"])]
        wording = "failed to give similar results" if label == "7j" else "Only 7a afforded similar products; this substrate was among the tested unsuccessful entries."
        add_step(r, "substrate_test", ["c60", "tbunc", "alkyne"], r["evidence_ids"], "Substrate tested while stoichiometry and temperature were optimized; exact entry conditions are not reported")
        outcome(r, "reported_result", "other", [], r["evidence_ids"], wording, "no similar isolable product reported", product_role="no assigned product")
        classify_empty(r, r["evidence_ids"], "unresolved", "Failure to give similar products does not establish zero conversion or identify any alternative products; no intended bond label is imposed.")
        reactions.append(r)

    trace = base_reaction("scheme2_trace_7i", 3, categories=["substrate_scope"], location="Scheme 2/text trace 7i result", evidence=e3)
    trace["general_scheme_ids"] = [generic_id]
    set_names(trace, "Trace annulation product from diyne 7i", "trace fullerene annulation", e3)
    trace["participants"] = [participant("c60", "mol_c60", "reactant", e3), participant("tbunc", "mol_tbunc", "reactant", e3), participant("diyne", "mol_7i", "reactant", e3), participant("product", "mol_7i_trace_product", "product", e3, "trace product")]
    add_step(trace, "trace_reaction", ["c60", "tbunc", "diyne"], e3, "Reaction at a much higher temperature for a much longer time", ("much higher temperature", "much higher temperature"), ("much longer reaction time", "much longer reaction time"))
    outcome(trace, "trace_product", "trace", ["mol_7i_trace_product"], e3, "trace amount", "trace", product_role="trace product")
    classify(trace, ["C-C"], e3, limitation="The source states a similar trace product but gives no compound label or complete structure; C-C annulation classification is provisional and partial.")
    reactions.append(trace)

    def scheme2_hydration(suffix, substrate, products):
        r = base_reaction(suffix, 3, categories=["derivatization"], location="Scheme 2 hydration arrow", evidence=e23)
        r["general_scheme_ids"] = [generic_id]
        set_names(r, "Acid-promoted hydration of ketenimine methano-[60]fullerene", "ketenimine hydration/rearrangement", e23)
        r["participants"] = [participant("ketenimine", substrate, "reactant", e23), participant("tsoh", "mol_tsoh_h2o", "reagent", e23), participant("chcl3", "mol_chcl3", "solvent", e23)]
        for pid, mid, y, role in products:
            r["participants"].append(participant(pid, mid, "product", e23, role))
            outcome(r, f"yield_{pid}", "yield", [mid], e3, f"{y}%", {"number": y, "unit": "%"}, role)
        set_solvent(r, "mol_chcl3", "CHCl3", e3)
        add_step(r, "hydration", ["ketenimine", "tsoh", "chcl3"], e3, "Treatment with TsOH·H2O in CHCl3", ("r.t.", "r.t."), ("20 min", {"number": 20, "unit": "min"}))
        if len(products) == 1:
            classify(r, ["DB-C-O", "N-H"], e3)
        else:
            classify_empty(r, e3, "unresolved", "The single experiment gives amide 10a and cyano byproduct 11a through divergent ketenimine transformations; complete product-specific atom mapping is unresolved.")
        return r

    reactions.append(scheme2_hydration("scheme2_hydration_9a", "mol_9a", [("p10a", "mol_10a", 74, "main product"), ("p11a", "mol_11a", 21, "byproduct")]))
    reactions.append(scheme2_hydration("scheme2_hydration_9b", "mol_9b", [("p10b", "mol_10b", 96, "main product")]))

    mech3 = base_reaction("scheme3_proposed_mechanism", 3, "proposed_mechanism", ["mechanistic_study"], "Scheme 3", e3)
    set_names(mech3, "Proposed mechanism of methano- and cyclopenteno-fullerene formation", "proposed zwitterionic addition/cyclization mechanism", e3)
    add_note(mech3, "Author-proposed pathways: concerted [3+2], 5-exo-dig, and 3-exo-trig cyclizations via intermediates A-C; not an independently observed experiment.", e3)
    classify_empty(mech3, e3)
    reactions.append(mech3)

    def allene_ann(suffix, allene, product, y):
        r = base_reaction(suffix, 4, categories=["substrate_scope"], location="Scheme 4", evidence=e34)
        r["general_scheme_ids"] = [generic_id]
        set_names(r, "Four-component annulation of C60, tBuNC, allenoate, and water", "allenoate cyclopentanofullerene annulation", e34)
        r["participants"] = [participant("c60", "mol_c60", "reactant", e34), participant("tbunc", "mol_tbunc", "reactant", e34), participant("allene", allene, "reactant", e34), participant("water", "mol_water", "reactant", e4), participant("odcb", "mol_odcb", "solvent", e4), participant("product", product, "product", e34, "only isolable product")]
        set_solvent(r, "mol_odcb", "oDCB", e4)
        add_step(r, "four_component_annulation", ["c60", "tbunc", "allene", "water", "odcb"], e4,
                 "C60, tBuNC, allenoate, and H2O reacted in oDCB under N2", ("100 °C", {"number": 100, "unit": "°C"}), ("1 h", {"number": 1, "unit": "h"}), "N2")
        outcome(r, "product_yield", "yield", [product], e4, f"{y}%", {"number": y, "unit": "%"}, "only isolable product")
        classify(r, ["C-C", "DB-C-O", "N-H"], e4, limitation="C-C annulation plus amide-forming C=O and N-H bonds are inferred from Scheme 4/5; exact cage atom mapping, bond counts, hydrogen source assignment, and all bond-order changes remain unresolved.")
        return r

    reactions += [allene_ann("scheme4_12a", "mol_12a", "mol_13a", 23), allene_ann("scheme4_12b", "mol_12b", "mol_13b", 20)]

    mech5 = base_reaction("scheme5_proposed_mechanism", 4, "proposed_mechanism", ["mechanistic_study"], "Scheme 5", e4)
    set_names(mech5, "Proposed mechanism of cyclopentanofullerene formation exemplified with 13a", "proposed allenoate addition/5-endo-trig/hydration mechanism", e4)
    add_note(mech5, "Author-proposed intermediates D-F and favored 5-endo-trig closure followed by hydration; not an independently observed experiment.", e4)
    classify_empty(mech5, e4)
    reactions.append(mech5)

    bundle = st.empty(dataset["documents"])
    bundle["molecules"] = list(MOLECULES.values())
    bundle["reactions"] = reactions
    bundle["source_context"] = [
        {"context_id": f"{DOC}:context_spectra_4_6", "document_id": DOC, "kind": "characterization",
         "evidence_ids": [ev(1, "text"), ev(2, "text")], "summary": "MALDI-TOF molecular weights and 1D/2D NMR support products 4-6; standalone spectra are not separate syntheses."},
        {"context_id": f"{DOC}:context_xray_4a_5a", "document_id": DOC, "kind": "structure_assignment",
         "evidence_ids": [ev(2, "text"), ev(2, "image")], "summary": "Single-crystal structures confirm trans 4a and ketenimine methanofullerene 5a."},
        {"context_id": f"{DOC}:context_4b_stereoisomers", "document_id": DOC, "kind": "stereochemical_characterization",
         "evidence_ids": [ev(2, "text")], "summary": "4b shows a minor 13C NMR peak set at about one-fifth the major intensity; major trans assignment is hypothesized."},
        {"context_id": f"{DOC}:context_9b_xray", "document_id": DOC, "kind": "structure_assignment",
         "evidence_ids": [ev(3, "text"), ev(3, "image")], "summary": "A crystal structure is reported for 9b; this characterization is not an extra reaction."},
        {"context_id": f"{DOC}:context_13_mass_nmr", "document_id": DOC, "kind": "characterization",
         "evidence_ids": [ev(3, "text"), ev(4, "text"), ev(4, "image")], "summary": "MALDI-TOF and NMR establish water incorporation and relative stereochemistry of 13a/13b; peak/coupling tables are not synthetic experiments."},
        {"context_id": f"{DOC}:context_si_notice", "document_id": DOC, "kind": "availability_notice",
         "evidence_ids": [ev(4, "text"), ev(4, "image")], "summary": "The article advertises Supporting Information, but no SI file was explicitly supplied or associated in this workspace."},
    ]
    write(RUN / "submitted.json", bundle)

    stereo_checks = {}
    for mid, molecule in MOLECULES.items():
        smi = st.value(molecule["smiles"])
        if smi:
            stereo_checks[mid] = st.stereo_inspect(smi)
    write(RUN / "raw_model/stereo_checks.json", stereo_checks)
    write(RUN / "raw_model/source_mentions.json", st.source_mentions(dataset, ["4a", "4b", "5a", "5b", "7a", "7b", "7i", "9a", "9b", "12a", "12b", "13a", "13b"], 180))
    write(RUN / "raw_model/visual_transcription.json", {
        "document_id": DOC,
        "pages_inspected": [1, 2, 3, 4, 5, 6],
        "scheme1": {"conditions": "C60, RNC, DMAD, oDCB, 80 °C, 1 h", "yields": {"4a": "7%", "5a": "36%", "4b": "6%", "5b": "26%", "6a": "95%", "6b": "92%"}, "hydration": "TsOH·H2O, CHCl3, r.t., 20 min"},
        "scheme2": {"conditions": "C60, tBuNC, oDCB, N2, 120 °C, 1 h", "yields": {"8a": "4%", "9a": "19%", "9b": "32%", "10a": "74%", "10b": "96%", "11a": "21%"}, "unsuccessful": ["7c", "7d", "7e", "7f", "7g", "7h", "7j"], "trace": ["7i"], "hydration": "TsOH·H2O, CHCl3, r.t., 20 min"},
        "scheme4": {"conditions": "C60, tBuNC, H2O, oDCB, N2, 100 °C, 1 h", "yields": {"13a": "23%", "13b": "20%"}},
        "non_experiments": ["page 1 graphical abstract overview", "Scheme 3 proposed mechanism", "Scheme 5 proposed mechanism", "Figure 1 and spectral/coupling characterization", "pages 5-6 bibliography"]
    })

    ids = {r["reaction_id"]: r for r in reactions}
    coverage = {"document_id": DOC, "semantic_accuracy_verified": False, "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title/abstract and graphical abstract", "finding": "Generic overview of current alkyne/allenoate annulation families; not an additional experiment.", "evidence_ids": [ev(1, "image"), ev(1, "text")]},
            {"region": "right-column Scheme 1b discussion", "finding": "Current C60/DMAD annulations, acid hydrations, and gradual 5b-to-6b conversion in air are reported; exact charges are absent.", "evidence_ids": [ev(1, "image"), ev(1, "text")]},
        ], "record_ids": [generic_id, f"{DOC}:scheme1b_tbunc", f"{DOC}:scheme1b_cync", f"{DOC}:scheme1b_hydration_5a", f"{DOC}:scheme1b_hydration_5b", f"{DOC}:air_conversion_5b"], "unresolved": ["Graphical abstract contains no independent quantities.", "Main text gives no exact reagent charges or yield measurement basis."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Scheme 1a", "finding": "Two prior-work examples (C60 products 1/2 and Lu3N@C80 product 3); classified as literature examples.", "evidence_ids": [ev(2, "image"), ev(2, "text")]},
            {"region": "Scheme 1b", "finding": "Two annulation experiments (tBuNC and CyNC) and two acid-hydration experiments with explicit conditions/yields.", "evidence_ids": [ev(2, "image"), ev(2, "text")]},
            {"region": "Figure 1 and adjacent prose", "finding": "X-ray/13C NMR structure and stereoisomer evidence; characterization, not extra synthesis. A separate toluene repeat yields 4a and 5a without numeric yields.", "evidence_ids": [ev(2, "image"), ev(2, "text")]},
        ], "record_ids": [f"{DOC}:literature_c60_dmad", f"{DOC}:literature_lu3n_c80_dmad", f"{DOC}:scheme1b_tbunc", f"{DOC}:scheme1b_cync", f"{DOC}:scheme1b_hydration_5a", f"{DOC}:scheme1b_hydration_5b", f"{DOC}:toluene_repeat"], "unresolved": ["Yield measurement basis is not stated.", "The exact major/minor 4b stereoisomer ratio and minor configuration are not assigned.", "The toluene-repeat temperature, time, charges, and yields are unreported."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Scheme 2 successful arrows", "finding": "7a and 7b annulation experiments and 9a/9b acid hydrations with explicit yields and common arrow conditions.", "evidence_ids": [ev(3, "image"), ev(3, "text")]},
            {"region": "Scheme 2 unsuccessful-substrate panel and prose", "finding": "7c-7h and 7j do not afford similar products; 7i gives a similar product in trace amount only at longer time/higher temperature. Each substrate test is separate.", "evidence_ids": [ev(3, "image"), ev(2, "text"), ev(3, "text")]},
            {"region": "Scheme 3", "finding": "Author-proposed mechanism via A-C; no intermediates are treated as isolated experiments.", "evidence_ids": [ev(3, "image"), ev(3, "text")]},
            {"region": "right-column allenoate prose", "finding": "Introduces separate 12a and 12b four-component annulations and MALDI/NMR characterization; the synthesis arrows are on page 4.", "evidence_ids": [ev(3, "image"), ev(3, "text")]},
        ], "record_ids": [f"{DOC}:scheme2_7a", f"{DOC}:scheme2_7b"] + [f"{DOC}:scheme2_unsuccessful_{x}" for x in ("7c", "7d", "7e", "7f", "7g", "7h", "7j")] + [f"{DOC}:scheme2_trace_7i", f"{DOC}:scheme2_hydration_9a", f"{DOC}:scheme2_hydration_9b", f"{DOC}:scheme3_proposed_mechanism", f"{DOC}:scheme4_12a", f"{DOC}:scheme4_12b"], "unresolved": ["Failed substrate entries lack exact individual conditions and do not establish zero conversion.", "7i trace product has no compound number or complete structure.", "Scope yield measurement bases are not stated."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "Scheme 4 reaction arrows", "finding": "Two separate four-component annulations, 12a to 13a (23%) and 12b to 13b (20%), under the printed common conditions.", "evidence_ids": [ev(4, "image"), ev(4, "text")]},
            {"region": "Scheme 4 stereochemical drawings and coupling/HMBC table", "finding": "Relative stereochemistry of 13a/13b is supported by coupling and HMBC observations; these measurements are source_context, not extra reactions.", "evidence_ids": [ev(4, "image"), ev(3, "text"), ev(4, "text")]},
            {"region": "Scheme 5", "finding": "Proposed D-F mechanism and hydration; classified as proposed mechanism, not experiment.", "evidence_ids": [ev(4, "image"), ev(4, "text")]},
            {"region": "Associated Content", "finding": "SI availability is advertised, but no SI file is supplied or explicitly associated in this workspace.", "evidence_ids": [ev(4, "image"), ev(4, "text")]},
        ], "record_ids": [f"{DOC}:scheme4_12a", f"{DOC}:scheme4_12b", f"{DOC}:scheme5_proposed_mechanism"], "unresolved": ["Absolute configurations of 13a/13b are not assigned.", "Dense fullerene products are not encoded as complete isomeric SMILES.", "Yield measurement basis and exact charges are unreported."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "author information, acknowledgments, references 1-9", "finding": "No current-paper reaction arrow, experimental table row, or synthesis procedure; bibliography only after front matter.", "evidence_ids": [ev(5, "image"), ev(5, "text")]},
        ], "record_ids": [], "unresolved": ["Reference descriptions were not converted into current-paper experiments."]},
        {"page": 6, "state": "attempted", "regions": [
            {"region": "references 9-16", "finding": "Bibliography only; no current-paper experimental record.", "evidence_ids": [ev(6, "image"), ev(6, "text")]},
        ], "record_ids": [], "unresolved": ["No reaction extraction required from bibliography entries."]},
    ]}
    assert set(ids) == {rid for p in coverage["pages"] for rid in p["record_ids"]}
    write(RUN / "coverage/pages_0001_0006.json", coverage)

    tasks_path = RUN / "tasks.json"
    tasks = json.loads(tasks_path.read_text())
    for task in tasks:
        if task["document_id"] == DOC:
            task["state"] = "attempted"
            task["reason"] = "Original rendered page inspected; detailed regions and record mappings are in coverage/pages_0001_0006.json."
    write(tasks_path, tasks)

    workflow = f"""# Workflow log — original PDF pages 1–6 of 6

- Batch: `allene-stereo-v31-200-20261017`
- Document: `{DOC}`
- Title: *Isocyanide-Induced Annulation Leading to Cyclopento-, Methano-, and Cyclopentano-[60]Fullerene Derivatives*
- Journal: *Organic Letters* 2021, 23, 8867–8872
- DOI: `10.1021/acs.orglett.1c03371`
- Source SHA-256: `89290e15b92920b2daf5e57736b456b2b1458085488aad10ddd83c6042a7efeb`
- Source mapping: `presentation_conversion=false`; all six inspected pages are original journal-PDF pages.
- Visually inspected: original PDF pages 1, 2, 3, 4, 5, and 6, exactly the supplied range.
- SI association: `unestablished`. Page 4 advertises SI, but no SI file was explicitly supplied or associated.

## Isolation from stale workspace content

The pre-existing workflow described a different article (`doc_0c91dd7389d7f800`, DOI `10.1021/acs.orglett.1c03270`) and 57 unrelated records. It was treated as cross-paper contamination and fully replaced. No molecule, condition, outcome, record ID, or classification from that paper was reused. The prepared current dataset contained zero candidates and zero experiments before this extraction. Shared tools, schemas, taxonomies, scientific quality rules, and tests were not edited.

## Visual inspection and raw artifacts

Every rendered page image was opened at original detail. Crops of Scheme 1, Scheme 2, Scheme 3, Scheme 4 with its coupling/HMBC table, and Scheme 5 are under `runs/current/raw_model/`; each crop has a tool-generated evidence receipt. `visual_transcription.json` preserves the literal visible conditions/yields, `stereo_checks.json` preserves RDKit graph diagnostics for candidate SMILES, and `source_mentions.json` preserves literal cross-page mention hints. Coverage decisions are in `runs/current/coverage/pages_0001_0006.json`; every page is `attempted` and `semantic_accuracy_verified=false`.

## Record-boundary decisions

1. The page-1 graphical abstract is one `generic_scheme`, not another specific experiment.
2. Scheme 1a is prior work: the C60 products 1/2 and Lu3N@C80 product 3 are two `literature_example` records, not current-paper experiments.
3. Scheme 1b contains two annulation experiments (tBuNC and CyNC) and two separate acid-hydration experiments. Each yield targets its own product. The 4b minor stereoisomer NMR signal is a selectivity outcome, not a second synthesis.
4. The prose-reported toluene repeat is a separate solvent-control experiment because its solvent differs; yields and all other conditions are unreported. Gradual 5b-to-6b conversion in air is retained as a separate qualitative derivatization observation.
5. Scheme 2 contains separate 7a and 7b successful annulations, seven distinct unsuccessful substrate tests (7c–h and 7j), one trace 7i experiment, and separate hydration experiments for 9a and 9b. “Failed to give similar results” is not rewritten as zero conversion; no desired product or bond label is invented. Exact failed-entry conditions remain unreported because the prose says a range was optimized.
6. Scheme 3 and Scheme 5 are `proposed_mechanism` records. Intermediates A–F are not treated as isolated compounds or verified experiments.
7. Scheme 4 contains two specific allenoate experiments, 12a→13a and 12b→13b, under the printed common conditions. The coupling/HMBC table is characterization in `source_context`, not extra syntheses.
8. Pages 5–6 are author/bibliography content. No reference entry, molecular-weight statement, MALDI signal, NMR peak, or coupling-constant row was made into a synthesis record.

## Molecules and stereochemistry

All registered molecules explicitly carry `chirality_status`, `geometry_status`, and `representation_status`. Candidate SMILES were generated only for visually/name-explicit small inputs and solvents, then checked with `standard stereo`; parsing/canonicalization remains inferred and does not certify source identity. Dense fullerene products remain source-labelled with null concrete SMILES.

- 4a retains the X-ray-confirmed trans relationship between tert-butyl and the fullerene cage; 4b is a stereoisomer mixture with a hypothesized trans major component and an unassigned minor component.
- 7b preserves the depicted E alkene in the candidate isomeric SMILES; 9b/10b retain source-drawing geometry without a fabricated complete fullerene SMILES.
- 12a is a potentially axially chiral allene and has no printed R_a/S_a assignment; connectivity is recorded with `chirality_status=unspecified` and `representation_status=partial`. 12b has a terminal CH2 allene and no assigned axial configuration.
- 13a/13b retain only the reported relative cis/trans relationships from coupling/HMBC evidence. No absolute configuration, ee, er, or dr is invented; their dense structures remain `not_encoded`.
- No CAS numbers were reported and none were guessed.

## Conditions, outcomes, and classification

The extraction contains 25 formal v3.1 records: 20 specific experiments/observations, one generic scheme, two literature examples, and two proposed mechanisms. All machine-imported records remain `needs_review`; none is source-reviewed. Explicit arrow temperatures/times/atmospheres are paired within ordered reaction steps. The main article reports no exact charges, scale, solvent volumes, concentration, workup, or purification parameters, so these remain unreported rather than inferred from the unavailable SI.

All numeric yields and the 4b selectivity observation explicitly target molecule IDs. The main article does not state whether the displayed yields are isolated, NMR, or another basis; measurement-basis fields remain null. Bond classifications are class-level, inferred, and partial: fullerene annulations receive C-C; clear ketenimine-to-amide hydrations receive DB-C-O and N-H; four-component allenoate products receive C-C, DB-C-O, and N-H. Exact atom maps/counts, bond-order decreases, hydrogen provenance, and coproduct balance remain unresolved. Failed substrate tests and divergent 9a hydration are left unclassified. No model votes were fabricated.

## Validation and unresolved limitations

`standard validate` and strict dataset validation pass. Final audit flags 11 records: 10 have `measurement_basis` warnings because the main article gives percentages without stating isolated/NMR/other basis, and 5 have `structure_gap` warnings because all primary reactant/product structures in those dense fullerene hydration/air-conversion records remain unencoded. The warning sets overlap. There are zero `reaction_identity`, zero `outcome_target`, and zero `missing_steps` warnings. The missing yield bases are not filled from convention, and dense structures are not guessed merely to silence audit.

Remaining scientific limitations are: no associated SI; no exact charges/workup/purification; unreported yield methods; incomplete dense fullerene structures/atom mapping; unassigned absolute stereochemistry of the 4b minor isomer and 13a/13b; unreported exact conditions and alternative products for failed Scheme 2 entries; and unidentified trace 7i product. Generated Markdown was checked for the printed reagent/solvent identities, all available temperatures/times, yields, and the absence of any source quantity conflict. Schema validation is not chemical accuracy certification, and page coverage is attempted rather than verified complete.
"""
    (ROOT / "runs/workflow.md").write_text(workflow)
    print(json.dumps({"molecules": len(MOLECULES), "reactions": len(reactions), "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in reactions)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
