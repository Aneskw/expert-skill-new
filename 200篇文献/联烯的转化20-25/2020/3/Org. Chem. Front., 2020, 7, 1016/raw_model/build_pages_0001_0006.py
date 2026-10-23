#!/usr/bin/env python3
"""Build evidence-linked v3.1 records for DOI 10.1039/d0qo00140f, PDF pages 1--6."""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_673d4d3f5080ad90"
RUN = ROOT / "runs/current"


def ev(page: int, kind: str = "image") -> str:
    return f"{DOC}:page{page}:{kind}"


def sf(raw, normalized=None, page=2, method="original_extraction", derivation=None,
       status=None):
    return st.fact(raw, normalized if normalized is not None else raw, [ev(page)],
                   method=method, derivation=derivation, status=status)


def inferred(value, page, why):
    return st.fact(None, value, [ev(page)], method="inference", derivation=why)


def set_stereo(m, page, chirality, geometry, representation, *, atom_centers=None,
               axial=None, assignment=None, mixture=None, double_bonds=None):
    e = [ev(page)]
    why = "Status assigned from the supplied source drawing; it is not an independent stereochemical determination."
    for key, val in (("chirality_status", chirality), ("geometry_status", geometry),
                     ("representation_status", representation)):
        m["stereochemistry"][key] = st.fact(None, val, e, method="inference", derivation=why)
    details = {
        "atom_centers": atom_centers,
        "axial": axial,
        "assignment": assignment,
        "mixture": mixture,
        "double_bonds": double_bonds,
    }
    for key, val in details.items():
        if val is not None:
            m["stereochemistry"][key] = st.fact(val, val, e)


def molecule(mid, label, page, *, name=None, smiles=None, stereo="unresolved",
             geometry="unresolved", representation=None, assignment=None,
             axial=None, atom_centers=None, double_bonds=None):
    e = [ev(page)]
    m = st.new_molecule(mid, DOC, e)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = st.fact(label, label, e)
    m["names"]["original"] = st.fact(name or label, name or label, e)
    m["identity_note"] = st.fact(
        f"Source identity: {label}; no CAS number is printed on covered pages.",
        f"Source-identified as {label}; CAS unreported.", e)
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = st.fact(
            label, smiles, e, method="structure_derivation",
            derivation="Machine candidate transcribed from the source name/drawing; parsing is not source verification.")
        diag = st.stereo_inspect(smiles)
        rep = representation or "encoded"
        set_stereo(m, page, stereo, geometry, rep, assignment=assignment,
                   axial=axial, atom_centers=atom_centers, double_bonds=double_bonds)
        m["stereochemistry"]["isomeric_smiles"] = st.fact(
            None, diag["canonical_isomeric_smiles"], e, method="structure_derivation",
            derivation="RDKit canonical isomeric SMILES from the machine candidate; source_verified=false.")
        m["stereochemistry"]["identity_key"] = st.fact(
            None, diag["stereo_graph_key"], e, method="structure_derivation",
            derivation="RDKit stereo graph key computed from the machine candidate; not a sample identity.")
    else:
        m["representation_type"] = "unresolved"
        set_stereo(m, page, stereo, geometry, representation or "not_encoded",
                   assignment=assignment, axial=axial, atom_centers=atom_centers,
                   double_bonds=double_bonds)
    return m


def generic_molecule(mid, label, page, fixed, variables, reason):
    m = molecule(mid, label, page, stereo="unresolved", geometry="unspecified",
                 representation="not_encoded")
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": sf(fixed, {"description": fixed}, page),
        "variables": [
            {
                "label": v,
                "attachment_points": sf(f"{v} attachment shown in source drawing", {"label": v}, page),
                "definition": sf(desc, desc, page),
                "allowed_values": st.fact(evidence=[ev(page)]),
                "constraints": st.fact(evidence=[ev(page)]),
            }
            for v, desc in variables
        ],
        "instance_mapping": [],
        "reason": sf(reason, reason, page),
    }
    return m


def charge(p, key, raw, normalized, page=2, reference=None):
    p["charges"][key] = sf(raw, normalized, page)
    if reference:
        p["equivalence_reference"] = sf(reference, reference, page)


def participant(pid, mid, role, page=2):
    return st.new_participant(pid, mid, role, [ev(page)])


def classification_cc(page, generic=False):
    reason = ("The source [3+2] annulation drawing closes a five-membered ring through two new C-C sigma connections; "
              "atom mapping and the complete bond-order inventory remain unresolved.")
    changes = [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
         "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
         "before_order": 0, "after_order": 1, "change_type": "formed"},
    ]
    return {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": st.fact(evidence=[ev(page)]),
        "bond_changes": [st.fact(None, c, [ev(page)], method="structure_derivation", derivation=reason) for c in changes],
        "labels": st.fact(None, [{"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1", "b2"]}],
                          [ev(page)], method="structure_derivation", derivation=reason),
        "coverage": inferred("partial", page, reason),
        "limitations": [
            "No atom-mapped reaction is asserted.",
            "Implicit-H changes and bond-order changes are unresolved.",
            "Classification is derived from the net ring-forming drawing, not a proposed elementary mechanism.",
        ],
        "model_votes": [],
    }


def base_reaction(rid, page, record_type, categories, location):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, [ev(page)])
    r["record_type"] = record_type
    r["experiment_categories"] = sf(categories, categories, page)
    r["location"] = sf(location, {"page": page, "location": location}, page)
    r["reaction_name"] = sf("[3 + 2] annulation", "enantioselective phosphine-catalyzed [3+2] annulation", page)
    r["transformation_type"] = sf("[3 + 2] annulation", "annulation/cycloaddition", page)
    r["review"] = {
        "status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
        "unresolved": ["Machine extraction remains needs_review; schema validity is not chemistry verification."],
    }
    r["classification"] = classification_cc(page)
    return r


def add_standard_participants(r, allene_mid, alkene_mid, product_mid, catalyst_mid,
                              solvent_mid, page, catalyst_loading=True):
    ps = [
        participant("allene", allene_mid, "reactant", page),
        participant("alkene", alkene_mid, "reactant", page),
        participant("product", product_mid, "product", page),
        participant("catalyst", catalyst_mid, "catalyst", page),
        participant("base", "m_k2co3", "additive", page),
        participant("solvent", solvent_mid, "solvent", page),
    ]
    charge(ps[0], "amount", "0.18 mmol", {"number": 0.18, "unit": "mmol"}, page)
    charge(ps[1], "amount", "0.15 mmol", {"number": 0.15, "unit": "mmol"}, page)
    if catalyst_loading:
        charge(ps[3], "mol_percent", "10 mol%", {"number": 10, "unit": "mol%"}, page, alkene_mid)
    charge(ps[4], "mol_percent", "10 mol%", {"number": 10, "unit": "mol%"}, page, alkene_mid)
    r["participants"] = ps
    r["reaction_scale"] = sf(
        "2/2a (0.15 mmol)", {"reference_molecule_id": alkene_mid, "amount": {"number": 0.15, "unit": "mmol"}}, page)


def reaction_step(r, page, temperature, time, atmosphere="N2"):
    s = st.new_step("reaction", [ev(page)])
    s["participant_ids"] = [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"]
    s["operation"] = sf("Reaction conditions", "annulation under the tabulated/scheme conditions", page)
    if temperature is not None:
        norm = temperature if isinstance(temperature, str) else {"number": temperature, "unit": "°C"}
        raw = "r.t." if temperature == "r.t." else f"{temperature} °C"
        s["temperature"] = sf(raw, norm, page)
    if time is not None:
        s["time"] = sf("4–6 h", {"minimum": 4, "maximum": 6, "unit": "h", "raw_text": "4–6 h"}, page)
    if atmosphere is not None:
        s["atmosphere"] = sf("under a N2 atmosphere", "N2", page)
    r["steps"] = [s]


def outcome(r, oid, kind, target, raw, normalized, page, *, basis_raw=None,
            basis_norm=None, scope=None, yield_type=None):
    o = st.new_outcome(oid, kind, [target], [ev(page)])
    o["product_role"] = sf("product", "main product", page)
    o["value"] = sf(raw, normalized, page)
    if basis_raw is not None:
        o["measurement_basis"] = sf(basis_raw, basis_norm or basis_raw, page)
    if scope is not None:
        o["scope"] = sf(scope, scope, page)
    if yield_type is not None:
        o["yield_type"] = sf(yield_type, yield_type, page)
    r["outcomes"].append(o)


def add_common_outcomes(r, product_mid, yld, ee_value, page):
    outcome(r, "yield", "yield", product_mid, f"{yld}%", {"number": yld, "unit": "%"}, page,
            basis_raw="isolated yield", basis_norm="isolated", yield_type="isolated yield")
    outcome(r, "dr", "selectivity", product_mid, "dr > 19 : 1",
            {"metric": "dr", "relation": ">", "major": 19, "minor": 1}, page,
            basis_raw="determined by 1H-NMR analysis of the crude product(s)",
            basis_norm="1H NMR of crude product", scope={"metric": "dr", "isomer_order": "major:minor"})
    if ee_value == "—":
        outcome(r, "ee", "selectivity", product_mid, "—", "table dash (—)", page,
                basis_raw="the ee value was determined by HPLC", basis_norm="HPLC",
                scope={"metric": "ee", "major_enantiomer_configuration": "unassigned"})
    else:
        sign = "negative" if ee_value < 0 else "positive_or_unsigned"
        outcome(r, "ee", "selectivity", product_mid, f"{ee_value}%", {"number": abs(ee_value), "unit": "%", "reported_sign": sign}, page,
                basis_raw="the ee value was determined by HPLC", basis_norm="HPLC",
                scope={"metric": "ee", "major_enantiomer_configuration": "unassigned", "reported_sign": sign})


def main():
    dataset = json.loads((RUN / "dataset.json").read_text())
    bundle = st.empty(dataset["documents"])
    mols = []

    # Generic families in Scheme 1; R groups are retained as templates rather than wildcard SMILES.
    mols += [
        generic_molecule("m_generic_1", "1", 2, "gamma-phenyl allenoate bearing CO2R2", [("R2", "ester alkyl/benzyl group")],
                         "Scheme 1 is an R-group family drawing."),
        generic_molecule("m_generic_2", "2", 2, "pyrazoloneylidene oxindole bearing R1 on oxindole N", [("R1", "oxindole N-substituent")],
                         "Scheme 1 is an R-group family drawing."),
        generic_molecule("m_generic_3", "3", 2, "spiro-cyclopentene-linked oxindole/pyrazolone bearing R1 and CO2R2", [("R1", "oxindole N-substituent"), ("R2", "ester alkyl/benzyl group")],
                         "Product is a generic R-group family with stereochemical wedges/dashes."),
    ]
    set_stereo(mols[-1], 2, "specified", "not_applicable", "not_encoded",
               atom_centers="Three stereogenic centres are shown by wedge/dash bonds; exact CIP descriptors are not printed.",
               assignment="Relative/depicted configuration as drawn; absolute descriptors unreported in Scheme 1.")

    # Concrete allenes. The allene axis is potentially stereogenic and is not assigned by the source.
    allene_specs = {
        "a": ("1a", "ethyl 4-phenylbuta-2,3-dienoate", "CCOC(=O)C=C=Cc1ccccc1", 2),
        "b": ("1b", "methyl 4-phenylbuta-2,3-dienoate", "COC(=O)C=C=Cc1ccccc1", 3),
        "c": ("1c", "benzyl 4-phenylbuta-2,3-dienoate", "O=C(OCc1ccccc1)C=C=Cc1ccccc1", 3),
    }
    for key, (label, name, smi, page) in allene_specs.items():
        mols.append(molecule(f"m_1{key}", label, page, name=name, smiles=smi,
                             stereo="unspecified", geometry="not_applicable", representation="partial",
                             axial="Potential allene axial chirality is not assigned in the source drawing.",
                             assignment="No Ra/Sa or enantiomeric composition reported for the allene substrate."))

    # Electrophiles 2a--2s are source-identified drawings; complete machine SMILES are withheld.
    for letter in "abcdefghijklmnopqrs":
        page = 2 if letter == "a" else 3
        mols.append(molecule(f"m_2{letter}", f"2{letter}", page,
                             name=f"pyrazoloneylidene oxindole {letter}",
                             stereo="unspecified", geometry="unspecified", representation="not_encoded",
                             double_bonds="Tetrasubstituted exocyclic alkene is drawn without an explicit E/Z descriptor.",
                             assignment="No absolute/relative stereochemical assignment is printed for this starting material."))

    # Scope products 3aa--3as, plus 3bk and 3ck.
    product_labels = [f"3a{x}" for x in "abcdefghijklmnopqrs"] + ["3bk", "3ck"]
    for label in product_labels:
        assignment = "Three stereogenic centres are depicted by wedge/dash bonds; exact CIP descriptors are not printed."
        if label == "3ag":
            assignment += " Absolute configuration is stated to be confirmed by X-ray crystallography (CCDC 1961000)."
        mols.append(molecule(f"m_{label}", label, 3, name=f"spiro-cyclopentene product {label}",
                             stereo="specified", geometry="not_applicable", representation="not_encoded",
                             atom_centers="Three source-depicted stereogenic centres, including two vicinal spiroquaternary centres.",
                             assignment=assignment))

    # Catalysts, additive, and solvents. Source drawings preserve catalyst stereochemistry without forcing SMILES.
    mols.append(molecule("m_pph3", "PPh3", 2, name="triphenylphosphine",
                         smiles="P(c1ccccc1)(c1ccccc1)c1ccccc1", stereo="achiral", geometry="not_applicable"))
    for i in range(1, 9):
        mols.append(molecule(f"m_p{i}", f"P{i}", 2, name=f"ferrocene-derived bifunctional phosphine P{i}",
                             stereo="specified", geometry="not_applicable", representation="not_encoded",
                             atom_centers="Source catalyst drawing contains explicit wedge stereochemistry.",
                             assignment="Configuration retained as drawn; complete stereochemical SMILES not transcribed."))
    mols += [
        molecule("m_k2co3", "K2CO3", 2, name="potassium carbonate", smiles="[K+].[K+].[O-]C([O-])=O", stereo="achiral", geometry="not_applicable"),
        molecule("m_toluene", "Toluene", 2, name="toluene", smiles="Cc1ccccc1", stereo="achiral", geometry="not_applicable"),
        molecule("m_et2o", "Et2O", 2, name="diethyl ether", smiles="CCOCC", stereo="achiral", geometry="not_applicable"),
        molecule("m_dcm", "DCM", 2, name="dichloromethane", smiles="ClCCl", stereo="achiral", geometry="not_applicable"),
        molecule("m_dmf", "DMF", 2, name="N,N-dimethylformamide", smiles="CN(C)C=O", stereo="achiral", geometry="not_applicable"),
    ]
    bundle["molecules"] = mols

    reactions = []
    scheme = base_reaction("scheme1_generic", 2, "generic_scheme", [], "Scheme 1")
    scheme["participants"] = [participant("allene_family", "m_generic_1", "reactant"),
                              participant("alkene_family", "m_generic_2", "reactant"),
                              participant("product_family", "m_generic_3", "product")]
    scheme["review"]["unresolved"].append("Generic R1/R2 structures and exact stereochemical descriptors are not expanded.")
    reactions.append(scheme)

    optimization = [
        (1, "pph3", "toluene", "r.t.", 90, "—"),
        (2, "p1", "toluene", "r.t.", 76, -65),
        (3, "p2", "toluene", "r.t.", 78, -68),
        (4, "p3", "toluene", "r.t.", 77, -41),
        (5, "p4", "toluene", "r.t.", 75, -81),
        (6, "p5", "toluene", "r.t.", 76, 62),
        (7, "p6", "toluene", "r.t.", 28, -60),
        (8, "p7", "toluene", "r.t.", 83, 76),
        (9, "p8", "toluene", "r.t.", 85, 87),
        (10, "p8", "toluene", 0, 85, 92),
        (11, "p8", "toluene", -10, 86, 97),
        (12, "p8", "toluene", -20, 88, 98),
        (13, "p8", "et2o", -20, 83, 95),
        (14, "p8", "dcm", -20, 76, 92),
        (15, "p8", "dmf", -20, 20, 75),
    ]
    for entry, catalyst, solvent, temp, yld, ee_value in optimization:
        r = base_reaction(f"table1_entry{entry:02d}", 2, "specific_experiment", ["optimization"], f"Table 1, entry {entry}")
        r["general_scheme_ids"] = [f"{DOC}:scheme1_generic"]
        add_standard_participants(r, "m_1a", "m_2a", "m_3aa", f"m_{catalyst}", f"m_{solvent}", 2)
        charge(r["participants"][-1], "volume", "3.0 mL", {"number": 3.0, "unit": "mL"}, 2)
        r["solvent_conditions"]["composition"] = sf(solvent if solvent != "et2o" else "Et2O", solvent, 2)
        r["solvent_conditions"]["total_volume"] = sf("3.0 mL", {"number": 3.0, "unit": "mL"}, 2)
        reaction_step(r, 2, temp, "4–6 h")
        add_common_outcomes(r, "m_3aa", yld, ee_value, 2)
        r["procedure"]["original_text"] = sf(
            "1a (0.18 mmol), 2a (0.15 mmol), catalyst (10 mol%), K2CO3 (10 mol%) in solvent (3.0 mL), 4–6 h, under a N2 atmosphere.",
            "Table 1 footnote conditions", 2)
        r["review"]["unresolved"] += [
            "Addition order, work-up and purification details are not reported on pages 1–6.",
            "Negative ee signs identify opposite optical sense, but the major enantiomer configuration is unassigned.",
        ]
        reactions.append(r)

    scope_rows = [
        ("3aa", "a", "a", 88, 98), ("3ab", "a", "b", 90, 92), ("3ac", "a", "c", 91, 90),
        ("3ad", "a", "d", 85, 90), ("3ae", "a", "e", 82, 92), ("3af", "a", "f", 87, 90),
        ("3ag", "a", "g", 85, 92), ("3ah", "a", "h", 82, 96), ("3ai", "a", "i", 79, 95),
        ("3aj", "a", "j", 82, 94), ("3ak", "a", "k", 93, 90), ("3bk", "b", "k", 91, 92),
        ("3ck", "c", "k", 90, 80), ("3al", "a", "l", 95, 70), ("3am", "a", "m", 86, 88),
        ("3an", "a", "n", 87, 86), ("3ao", "a", "o", 85, 84), ("3ap", "a", "p", 86, 90),
        ("3aq", "a", "q", 80, 98), ("3ar", "a", "r", 92, 94), ("3as", "a", "s", 85, 88),
    ]
    for label, allene_letter, alkene_letter, yld, ee_value in scope_rows:
        scheme_no = 2 if label in {x[0] for x in scope_rows[:13]} else 3
        r = base_reaction(f"scheme{scheme_no}_{label}", 3, "specific_experiment", ["substrate_scope"], f"Scheme {scheme_no}, product {label}")
        r["general_scheme_ids"] = [f"{DOC}:scheme1_generic"]
        add_standard_participants(r, f"m_1{allene_letter}", f"m_2{alkene_letter}", f"m_{label}", "m_p8", "m_toluene", 3)
        r["solvent_conditions"]["composition"] = sf("toluene", "toluene", 3)
        reaction_step(r, 3, -20, None)
        add_common_outcomes(r, f"m_{label}", yld, ee_value, 3)
        r["procedure"]["original_text"] = sf(
            "1 (0.18 mmol), 2 (0.15 mmol), P8 (10 mol%), K2CO3 (10 mol%) in toluene, under a N2 atmosphere.",
            f"Scheme {scheme_no} caption conditions", 3)
        r["review"]["unresolved"] += [
            "Reaction time, toluene volume/concentration, addition order, work-up and purification are unreported on pages 1–6.",
            "Exact product CIP descriptors are not printed; wedge/dash depiction is retained without invented R/S labels.",
        ]
        if label == "3ag":
            r["notes"].append(sf("absolute configuration of 3ag was confirmed by X-ray crystallography (CCDC 1961000)",
                                  {"method": "X-ray crystallography", "ccdc": "1961000"}, 3))
        reactions.append(r)

    mechanism = base_reaction("scheme4_proposed_mechanism", 3, "proposed_mechanism", [], "Scheme 4 and page 4 discussion")
    mechanism["evidence_ids"] = [ev(3), ev(4)]
    mechanism["participants"] = [participant("allene", "m_1a", "reactant", 3),
                                  participant("alkene", "m_2a", "reactant", 3),
                                  participant("catalyst", "m_p8", "catalyst", 3),
                                  participant("product", "m_3aa", "product", 3)]
    mechanism["procedure"]["original_text"] = st.fact(
        "P8 attacks 1a to form Int-I; gamma-addition to 2 gives II; intramolecular Michael addition gives III; proton transfer and catalyst release gives 3.",
        "Author-proposed catalytic sequence Int-I -> II -> III -> 3", [ev(3), ev(4)])
    mechanism["notes"].append(st.fact(
        "Based on previous reports, a possible mechanism and an activation model were proposed.",
        "Proposed mechanism, not a separately measured experiment.", [ev(4)]))
    mechanism["review"]["unresolved"].append("Intermediates Int-I, II and III are proposed, not isolated or experimentally verified on covered pages.")
    reactions.append(mechanism)

    bundle["reactions"] = reactions
    bundle["source_context"] = [
        {"context_id": f"{DOC}:page1_figure1", "document_id": DOC, "evidence_ids": [ev(1)],
         "kind": "background_figure", "note": "Figure 1 shows cited bioactive SCPO/SCPP examples; not current synthesis experiments."},
        {"context_id": f"{DOC}:page3_ccdc", "document_id": DOC, "evidence_ids": [ev(3), ev(6)],
         "kind": "characterization", "note": "3ag X-ray structure; CCDC 1961000. No separate synthesis record created."},
        {"context_id": f"{DOC}:pages4_6_references", "document_id": DOC, "evidence_ids": [ev(4), ev(5), ev(6)],
         "kind": "references_and_conclusion", "note": "Conclusions and cited prior work; cited transformations are not treated as current experiments."},
    ]

    raw_path = RUN / "raw_model/extracted_bundle_pages_0001_0006.json"
    raw_path.write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
    (RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
    stereo_diagnostics = {
        m["molecule_id"]: st.stereo_inspect(st.value(m["smiles"]))
        for m in mols if m["representation_type"] == "smiles"
    }
    (RUN / "raw_model/stereo_diagnostics.json").write_text(
        json.dumps(stereo_diagnostics, ensure_ascii=False, indent=2) + "\n")
    mention_hints = st.source_mentions(dataset, ["1a", "2a", "3aa", "3ag", "P8"], context=120)
    (RUN / "raw_model/source_mentions.json").write_text(
        json.dumps(mention_hints, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"molecules": len(mols), "reactions": len(reactions),
                      "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in reactions)},
                     ensure_ascii=False))


if __name__ == "__main__":
    main()
