#!/usr/bin/env python3
"""Build the source-bounded v3.1 extraction for PDF pages 1--6."""

from __future__ import annotations

import copy
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st  # noqa: E402


DOC = "doc_93797a8609d47752"
E = {page: f"{DOC}:page{page}:text" for page in range(1, 7)}
I = {page: f"{DOC}:page{page}:image" for page in range(1, 7)}
OUT = ROOT / "runs/current"


def F(raw=None, normalized=None, ev=(), *, method="original_extraction", derivation=None,
      status=None):
    return st.fact(raw, normalized, ev, method=method, derivation=derivation, status=status)


def inferred(raw, normalized, ev, derivation):
    return F(raw, normalized, ev, method="inference", derivation=derivation)


def molecule(mid, label, name, ev, *, representation="unresolved", smiles=None,
             identity_status="source_identified", chirality="unresolved",
             geometry="unresolved", representation_status="unresolved",
             stereo_note=None, template=None, abbreviations=(), axial=None):
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = identity_status
    m["occurrences"][0]["compound_label"] = F(label, label, ev)
    m["names"]["original"] = F(name, name, ev)
    m["names"]["english"] = F(name, name, ev)
    if abbreviations:
        m["names"]["abbreviations"] = F(list(abbreviations), list(abbreviations), ev)
    m["representation_type"] = representation
    if smiles is not None:
        m["smiles"] = F(smiles, smiles, ev, method="structure_derivation",
                         derivation="Candidate connectivity transcribed from the source drawing/name; parsing does not verify source identity.")
        diagnostic = st.stereo_inspect(smiles)
        m["stereochemistry"]["isomeric_smiles"] = F(
            None, diagnostic["canonical_isomeric_smiles"], ev, method="structure_derivation",
            derivation="RDKit canonical isomeric SMILES from the candidate graph; this does not add an unreported stereochemical assignment.")
        m["stereochemistry"]["identity_key"] = F(
            None, diagnostic["stereo_graph_key"], ev, method="structure_derivation",
            derivation="SHA-256 stereo graph key calculated from the candidate isomeric graph; it is not a sample identity.")
        if diagnostic["double_bonds"]:
            m["stereochemistry"]["double_bonds"] = F(
                "double bonds in candidate graph", diagnostic["double_bonds"], ev, method="structure_derivation",
                derivation="RDKit graph diagnostic; STEREONONE preserves that no E/Z assignment is encoded.")
    m["generic_template"] = template
    observation = stereo_note or "No source stereochemical assignment was identified in the inspected pages."
    m["stereochemistry"]["chirality_status"] = inferred(
        observation, chirality, ev,
        "Status preserves the inspected source assignment or lack of assignment; absence of a descriptor is not treated as racemic or achiral.")
    m["stereochemistry"]["geometry_status"] = inferred(
        observation, geometry, ev,
        "Status preserves source-specified geometry separately from connectivity.")
    m["stereochemistry"]["representation_status"] = inferred(
        observation, representation_status, ev,
        "Status compares the retained source stereochemical information with the available molecular representation.")
    if stereo_note:
        m["stereochemistry"]["assignment"] = F(stereo_note, stereo_note, ev)
    if axial is not None:
        axial_raw, axial_normalized, axial_derivation = axial
        m["stereochemistry"]["axial"] = F(
            axial_raw, axial_normalized, ev, method="structure_derivation",
            derivation=axial_derivation)
    m["identity_note"] = F(
        f"Source identity retained as {label}; no CAS number is reported in inspected pages 1–6.",
        f"Source identity retained as {label}; no CAS number is reported in inspected pages 1–6.", ev)
    return m


def generic_template(ev, fixed_raw, fixed_normalized, variables, reason):
    return {
        "fixed_structure": F(fixed_raw, fixed_normalized, ev),
        "variables": [
            {
                "label": label,
                "attachment_points": F(points_raw, points_norm, ev),
                "definition": F(def_raw, def_norm, ev),
                "allowed_values": F(allowed_raw, allowed_norm, ev) if allowed_raw is not None else F(),
                "constraints": F(constraints_raw, constraints_norm, ev) if constraints_raw is not None else F(),
            }
            for label, points_raw, points_norm, def_raw, def_norm, allowed_raw, allowed_norm,
            constraints_raw, constraints_norm in variables
        ],
        "instance_mapping": [],
        "reason": F(reason, reason, ev),
    }


def participant(pid, mid, role, ev, *, product_role=None, mol_percent=None,
                equivalents=None, reference=None, stage="reaction"):
    p = st.new_participant(pid, mid, role, ev, stage=stage)
    if product_role is not None:
        p["product_role"] = F(product_role, product_role, ev)
    if mol_percent is not None:
        raw, number, derivation = mol_percent
        p["charges"]["mol_percent"] = F(
            raw, {"number": number, "unit": "mol%"}, ev, method="inference", derivation=derivation)
    if equivalents is not None:
        raw, number, method, derivation = equivalents
        p["charges"]["equivalents"] = F(
            raw, {"number": number, "unit": "equiv"}, ev, method=method, derivation=derivation)
    if reference is not None:
        p["equivalence_reference"] = inferred(
            f"loading/ratio relative to {reference}", reference, ev,
            "The source ratio or equivalents are expressed relative to this reaction participant.")
    return p


def reaction(suffix, ev, record_type, categories, location, name, transformation, unresolved):
    rid = f"{DOC}:{suffix}"
    r = st.new_reaction(rid, DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = F(categories, categories, ev)
    r["location"] = F(location, location, ev)
    r["reaction_name"] = F(name, name, ev)
    r["transformation_type"] = inferred(
        transformation, transformation, ev,
        "Transformation label is a normalized description of the source-labelled reaction, not an independent source claim.")
    r["review"]["unresolved"] = list(unresolved)
    return r


def set_solvent(r, raw, composition, ev, *, inherited=False, concentration=None):
    method = "inherited_general_condition" if inherited else "original_extraction"
    r["solvent_conditions"]["composition"] = F(raw, composition, ev, method=method)
    if concentration is not None:
        raw_c, number, basis, method_c = concentration
        r["solvent_conditions"]["concentration"] = F(
            raw_c, {"number": number, "unit": "M"}, ev, method=method_c,
            derivation="The comparison is stated as a repetition of the Figure 1 experiment." if method_c == "inference" else None)
        r["solvent_conditions"]["concentration_basis"] = F(basis, basis, ev, method=method_c,
            derivation="The concentration basis follows the explicitly stated monomer concentration." if method_c == "inference" else None)


def add_step(r, sid, ev, participants, operation, *, temperature=None, time=None,
             stage="reaction"):
    s = st.new_step(sid, ev, stage=stage)
    s["participant_ids"] = list(participants)
    s["operation"] = F(operation, operation, ev)
    if temperature is not None:
        raw, normalized, method, derivation = temperature
        s["temperature"] = F(raw, normalized, ev, method=method, derivation=derivation)
    if time is not None:
        raw, normalized, method, derivation = time
        s["time"] = F(raw, normalized, ev, method=method, derivation=derivation)
    r["steps"].append(s)


def add_outcome(r, oid, kind, targets, ev, raw, normalized, *, product_role=None,
                yield_type=None, basis=None, scope=None, method="original_extraction",
                derivation=None):
    o = st.new_outcome(oid, kind, targets, ev)
    o["value"] = F(raw, normalized, ev, method=method, derivation=derivation)
    if product_role is not None:
        o["product_role"] = F(product_role, product_role, ev)
    if yield_type is not None:
        o["yield_type"] = F(yield_type, yield_type, ev)
    if basis is not None:
        o["measurement_basis"] = F(basis, basis, ev)
    if scope is not None:
        o["scope"] = F(scope, scope, ev)
    r["outcomes"].append(o)


def classification(r, ev, mode="unresolved", *, nbe=False):
    c = r["classification"]
    c["mapped_reaction_smiles"] = F()
    c["model_votes"] = []
    if mode == "not_applicable":
        c["bond_changes"] = []
        c["labels"] = inferred("proposed pathway, not a net experiment", [], ev,
                               "Bond-formation taxonomy is not applied to a proposed elementary-step cycle.")
        c["coverage"] = inferred("proposed mechanism", "not_applicable", ev,
                                 "The record is a proposed mechanism rather than a standalone observed net transformation.")
        c["limitations"] = ["No three-model votes were run or fabricated."]
        return
    if mode == "unresolved":
        c["bond_changes"] = []
        c["labels"] = inferred("structure or atom correspondence unresolved", [], ev,
                               "No formation label is asserted without a reliable atom correspondence.")
        c["coverage"] = inferred("unresolved", "unresolved", ev,
                                 "The inspected source does not support a safe complete mapping.")
        c["limitations"] = ["No atom mapping was fabricated.", "No three-model votes were run or fabricated."]
        return
    changes = [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
         "before_order": 0, "after_order": 1, "change_type": "formed"}
    ]
    labels = [{"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1"]}]
    if nbe:
        changes.append(
            {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
             "before_order": 1, "after_order": 2, "change_type": "order_increased"})
        labels.append({"category_number": 13, "code": "DB-C-C", "supporting_change_ids": ["b2"]})
    c["bond_changes"] = [F(None, change, ev, method="structure_derivation",
                            derivation="Partial net bond change inferred from the drawn polymerization scheme; atom identities remain unresolved.")
                           for change in changes]
    c["labels"] = F(None, labels, ev, method="structure_derivation",
                    derivation="Taxonomy labels follow only the retained partial formed/promoted bonds.")
    c["coverage"] = inferred("partial scheme-level classification", "partial", ev,
                             "Polymer repeat-unit drawings do not establish a complete finite atom mapping or end-group balance.")
    c["limitations"] = [
        "Scheme-level polymerization classification only; complete atom mapping and chain-end stoichiometry are unresolved.",
        "No three-model votes were run or fabricated.",
    ]


def table_reaction(index, row, molecules):
    ratio, conc, temp, hours, conv, mn_theo, mn_sec, dispersity, catalyst_mid, catalyst_label = row
    ev = [E[3], I[3]]
    suffix = f"table1_entry{index:02d}"
    r = reaction(
        suffix, ev, "specific_experiment", ["optimization"],
        {"page": 3, "region": "Table 1", "entry": index},
        "Screening of Polymerization Conditions for Poly(1)",
        "ring-opening allene metathesis polymerization",
        ["Reaction scale, atmosphere, pressure, workup, and isolation yield are not reported in Table 1.",
         "THF is inherited from the adjacent statement that it was chosen as solvent of choice; the table itself does not repeat the solvent.",
         "Polymer repeat stereochemistry and chain-end structures are unassigned."],
    )
    molpct = 100.0 / ratio
    r["participants"] = [
        participant(f"{suffix}_monomer", "m_cyclononadiene_1", "reactant", ev),
        participant(f"{suffix}_catalyst", catalyst_mid, "catalyst", ev,
                    mol_percent=(f"[M]/[I] = {ratio}", molpct,
                                 f"Calculated as 100/[M]/[I] = {molpct:g} mol% relative to monomer 1."),
                    reference="m_cyclononadiene_1"),
        participant(f"{suffix}_solvent", "m_thf", "solvent", ev),
        participant(f"{suffix}_product", "m_poly1", "product", ev, product_role="polymer product"),
    ]
    set_solvent(r, "THF was chosen as our solvent of choice", ["THF"], ev, inherited=True,
                concentration=(f"[M]0 = {conc} M", conc, "monomer 1", "original_extraction"))
    add_step(
        r, f"{suffix}_reaction", ev, [p["participant_id"] for p in r["participants"]],
        f"Polymerization with {catalyst_label} under Table 1 conditions.",
        temperature=(f"{temp} °C", {"number": temp, "unit": "°C"}, "original_extraction", None),
        time=(f"{hours} h", {"number": hours, "unit": "h"}, "original_extraction", None),
    )
    add_outcome(r, f"{suffix}_conversion", "conversion", ["m_cyclononadiene_1"], ev,
                f"{conv}%", {"number": conv, "unit": "%"},
                product_role="reactant consumption",
                basis="1H NMR analysis in CDCl3 at 25 °C", scope="monomer conversion")
    add_outcome(r, f"{suffix}_mn_theo", "other", ["m_poly1"], ev,
                f"[M]n,theo. = {mn_theo} kDa", f"{mn_theo} kDa",
                product_role="polymer product",
                basis="([M]/[I] × conversion) × FW1", scope="theoretical number-average molar mass")
    add_outcome(r, f"{suffix}_mn_sec", "other", ["m_poly1"], ev,
                f"Mn = {mn_sec} kDa", f"{mn_sec} kDa",
                product_role="polymer product",
                basis="SEC (RI) in THF; 12-point conventional calibration with narrowly dispersed polystyrene standards",
                scope="number-average molar mass")
    add_outcome(r, f"{suffix}_dispersity", "other", ["m_poly1"], ev,
                f"Đ = {dispersity}", dispersity, product_role="polymer product",
                basis="SEC (RI) in THF; 12-point conventional calibration with narrowly dispersed polystyrene standards",
                scope="dispersity")
    r["procedure"]["original_text"] = F(
        {"table": "Table 1", "entry": index, "values": [ratio, conc, temp, hours, conv, mn_theo, mn_sec, dispersity]},
        {"table": "Table 1", "entry": index}, ev)
    classification(r, ev, "polymer")
    return r


def main():
    dataset_path = OUT / "dataset.json"
    dataset = json.loads(dataset_path.read_text())
    doc = dataset["documents"][0]
    doc["title"] = "Selective Ring-Opening Allene Metathesis: Polymerization or Ruthenium Vinylidene Formation"
    st.write(dataset_path, dataset)

    bundle = json.loads((OUT / "submitted.json").read_text())
    bundle["documents"] = copy.deepcopy(dataset["documents"])

    poly1_template = generic_template(
        [E[2], I[2], E[3], I[3]],
        "poly(1) repeat unit drawn with one allene per nine-carbon backbone repeat and two polymer continuation bonds",
        {"repeat_unit": "ring-opened 1,2-cyclononadiene connectivity", "repeat_mass_g_mol": 122.2,
         "polymer_continuation_points": 2},
        [("n", "polymer repeat count", ["repeat_unit"], "number of constitutional repeating units", "unreported positive integer",
          None, None, None, None)],
        "Polymer has variable chain length/end groups and is therefore represented as a source-defined polymer template rather than a finite SMILES.")
    polynbe_template = generic_template(
        [E[5], I[5]],
        "poly(NBE-iPr) repeat structure drawn in Figure 3",
        {"repeat_unit": "ring-opened NBE-iPr connectivity", "polymer_continuation_points": 2},
        [("n", "polymer repeat count", ["repeat_unit"], "number of constitutional repeating units", "unreported positive integer",
          None, None, None, None)],
        "Variable chain length, unassigned chain ends, and unresolved cis/trans microstructure prevent a finite concrete SMILES.")
    rcam_start_template = generic_template(
        [E[2], I[2]], "terminal allene pair joined by the drawn (CH2)12 tether",
        {"source_drawing": "allenyl-(CH2)12-allenyl precursor"}, [],
        "The cited-literature drawing is retained without assigning an unchecked concrete SMILES.")
    rcam_product_template = generic_template(
        [E[2], I[2]], "macrocyclic allene product drawn in the RCAM previous-work panel",
        {"source_drawing": "macrocyclic allene"}, [],
        "The cited-literature macrocycle is retained as a drawing-defined template.")
    cam_start_template = generic_template(
        [E[2], I[2]], "R-substituted terminal allene", {"fixed_fragment": "R-CH=C=CH2"},
        [("R", "bond to terminal allene carbon", ["terminal carbon"], "substituent R", "unspecified substituent",
          None, None, None, None)],
        "R is not expanded in the previous-work generic scheme.")
    cam_product_template = generic_template(
        [E[2], I[2]], "1,3-disubstituted allene", {"fixed_fragment": "R-CH=C=CH-R"},
        [("R", "bonds to both terminal allene carbons", ["two terminal carbons"], "substituent R", "unspecified substituent",
          None, None, "same R shown on both termini", "same generic R")],
        "R is not expanded in the previous-work generic scheme.")

    mols = [
        molecule("m_cyclononadiene_1", "1", "1,2-cyclononadiene", [E[1], I[1], E[2], I[2], E[3], E[4], I[4], E[5], I[5]],
                 representation="smiles", smiles="C1=C=CCCCCCC1", chirality="unspecified",
                 geometry="not_applicable", representation_status="partial",
                 stereo_note="The cyclic allene axis is drawn without R_a/S_a assignment; the candidate connectivity SMILES does not encode axial configuration.",
                 axial=("Cyclic allene C=C=C axis drawn without R_a/S_a descriptor",
                        {"axis": "cyclic C=C=C", "assignment": "unspecified"},
                        "The source drawing establishes the allene axis, but no absolute axial descriptor is supplied and RDKit does not encode it in the retained candidate SMILES.")),
        molecule("m_poly1", "poly(1)", "poly(1,2-cyclononadiene) / poly(1)", [E[1], I[1], E[2], I[2], E[3], I[3], E[4], I[4], E[5]],
                 representation="generic_template", template=poly1_template, chirality="unspecified",
                 geometry="unresolved", representation_status="not_encoded",
                 stereo_note="Allene-containing repeat units are shown, but axial configuration, any E/Z descriptors, chain ends, and microstructure distributions are not assigned.",
                 axial=("Poly(1) repeat unit is drawn with an allene but without R_a/S_a assignment",
                        {"axis": "repeat-unit C=C=C", "assignment": "unspecified"},
                        "The source repeat-unit drawings contain an allene axis; no repeat-unit axial configuration or distribution is reported.")),
        molecule("m_g1", "G1", "Grubbs first-generation catalyst", [E[1], I[1], E[2], I[2], E[4], E[5], I[5]], abbreviations=["G1"],
                 stereo_note="Catalyst structure is drawn, but coordination stereochemistry is not encoded in this extraction."),
        molecule("m_g2", "G2", "Grubbs second-generation catalyst", [E[1], I[1], E[2], I[2], E[3], E[4]], abbreviations=["G2"],
                 stereo_note="Catalyst structure is drawn, but coordination stereochemistry is not encoded in this extraction."),
        molecule("m_hg2", "HG2", "Hoveyda–Grubbs second-generation catalyst", [E[3], E[4]], abbreviations=["HG2"],
                 stereo_note="Catalyst is named by footnote; coordination stereochemistry is unresolved."),
        molecule("m_g3", "G3", "Grubbs third-generation catalyst", [E[3], E[4]], abbreviations=["G3"],
                 stereo_note="Catalyst is named by footnote; coordination stereochemistry is unresolved."),
        molecule("m_thf_d8", "THF-d8", "tetrahydrofuran-d8", [E[2], E[3]],
                 stereo_note="Isotopic substitution is named, but a complete isotope-resolved structure is not encoded."),
        molecule("m_toluene_d8", "toluene-d8", "toluene-d8", [E[3], E[4], E[5]],
                 stereo_note="Isotopic substitution is named, but a complete isotope-resolved structure is not encoded."),
        molecule("m_thf", "THF", "tetrahydrofuran", [E[3], E[4]], representation="smiles", smiles="C1CCOC1",
                 chirality="achiral", geometry="not_applicable", representation_status="not_applicable",
                 stereo_note="The retained THF connectivity has no stereogenic atom center or E/Z element."),
        molecule("m_methanol", "MeOH", "methanol", [E[3], E[5]], representation="smiles", smiles="CO",
                 chirality="achiral", geometry="not_applicable", representation_status="not_applicable",
                 stereo_note="The retained methanol connectivity has no stereogenic atom center or E/Z element."),
        molecule("m_g1_vinylidene", "G1-derived ruthenium vinylidene", "G1-derived ruthenium vinylidene", [E[1], I[1], E[4], I[4], E[5], I[5]],
                 chirality="unresolved", geometry="mixture", representation_status="not_encoded",
                 stereo_note="Figure 3 labels c and t as cis and trans resonances; their ratio and complete coordination structures are not reported in the inspected pages."),
        molecule("m_g1_unknown_polymer", "unassigned G1-derived polymer", "unassigned polymer from 1 and G1 at 50 °C", [E[4], E[5]],
                 chirality="unresolved", geometry="unresolved", representation_status="unresolved",
                 stereo_note="The authors state that they cannot assign the correct microstructure and observe no allene signal in the isolated material."),
        molecule("m_nbe_ipr", "NBE-iPr", "N-isopropyl norbornene derivative (NBE-iPr)", [E[1], E[5], I[5]], abbreviations=["NBE-iPr"],
                 chirality="unresolved", geometry="unresolved", representation_status="unresolved",
                 stereo_note="The Figure 3 drawing is source-identified, but its complete relative stereochemistry was not assigned to a checked SMILES."),
        molecule("m_poly_nbe_ipr", "poly(NBE-iPr)", "poly(NBE-iPr)", [E[1], E[5], I[5]],
                 representation="generic_template", template=polynbe_template, chirality="unresolved",
                 geometry="mixture", representation_status="not_encoded",
                 stereo_note="The polymer repeat is drawn, while backbone cis/trans distribution and chain ends are not assigned; no ratio is reported."),
        molecule("m_rcam_precursor_generic", "RCAM precursor", "generic tethered terminal-diallene precursor", [E[2], I[2]],
                 representation="generic_template", template=rcam_start_template, chirality="unspecified",
                 geometry="not_applicable", representation_status="not_encoded",
                 stereo_note="Terminal allene configurations are not assigned in the cited previous-work panel."),
        molecule("m_rcam_product_generic", "RCAM macrocyclic allene", "generic macrocyclic allene product", [E[2], I[2]],
                 representation="generic_template", template=rcam_product_template, chirality="unspecified",
                 geometry="not_applicable", representation_status="not_encoded",
                 stereo_note="The macrocyclic allene configuration is not assigned in the cited previous-work panel.",
                 axial=("Macrocyclic allene is drawn without an axial descriptor",
                        {"axis": "macrocyclic C=C=C", "assignment": "unspecified"},
                        "The cited previous-work product drawing contains a cyclic allene axis but provides no absolute descriptor.")),
        molecule("m_cam_terminal_allene_generic", "R-substituted terminal allene", "generic R-substituted terminal allene", [E[2], I[2]],
                 representation="generic_template", template=cam_start_template, chirality="unspecified",
                 geometry="not_applicable", representation_status="not_encoded",
                 stereo_note="R and any allene axial configuration are unspecified."),
        molecule("m_cam_product_generic", "1,3-disubstituted allene", "generic 1,3-disubstituted allene", [E[2], I[2]],
                 representation="generic_template", template=cam_product_template, chirality="unspecified",
                 geometry="not_applicable", representation_status="not_encoded",
                 stereo_note="R and any allene axial configuration are unspecified.",
                 axial=("Generic R-CH=C=CH-R product is drawn without R_a/S_a assignment",
                        {"axis": "C=C=C", "assignment": "unspecified", "depends_on": "R substituent identity"},
                        "The drawn 1,3-disubstituted allene can carry axial stereochemistry, but R and the axial descriptor are not specified.")),
    ]
    bundle["molecules"] = mols

    reactions = []

    r = reaction("scheme1_rcam_prior", [E[2], I[2]], "literature_example", ["other"],
                 {"page": 2, "region": "Scheme 1, Previous Work, RCAM panel"},
                 "Ring Closing Allene Metathesis (RCAM)", "ring-closing allene metathesis",
                 ["This is cited previous work, not a current-paper experiment.",
                  "The panel gives G1 15 mol% but no solvent, temperature, time, yield, or stereochemical assignment."])
    r["participants"] = [
        participant("rcam_start", "m_rcam_precursor_generic", "reactant", [E[2]]),
        participant("rcam_g1", "m_g1", "catalyst", [E[2]], mol_percent=("15 mol%", 15,
                    "Directly read from the Scheme 1 arrow."), reference="m_rcam_precursor_generic"),
        participant("rcam_product", "m_rcam_product_generic", "product", [E[2]], product_role="macrocyclic allene"),
    ]
    add_step(r, "rcam_arrow", [E[2]], [p["participant_id"] for p in r["participants"]], "Cited RCAM arrow with G1.")
    classification(r, [E[2]], "unresolved")
    reactions.append(r)

    r = reaction("scheme1_cam_prior", [E[2], I[2]], "literature_example", ["other"],
                 {"page": 2, "region": "Scheme 1, Previous Work, CAM panel"},
                 "Cross Allene Metathesis (CAM)", "cross-allene metathesis",
                 ["This is cited previous work, not a current-paper experiment.",
                  "R groups, solvent, temperature, time, yield, and stereochemical assignment are not supplied in the panel."])
    r["participants"] = [
        participant("cam_start", "m_cam_terminal_allene_generic", "reactant", [E[2]]),
        participant("cam_g1", "m_g1", "catalyst", [E[2]], mol_percent=("5 mol%", 5,
                    "Directly read from the Scheme 1 arrow."), reference="m_cam_terminal_allene_generic"),
        participant("cam_product", "m_cam_product_generic", "product", [E[2]], product_role="1,3-disubstituted allene"),
    ]
    add_step(r, "cam_arrow", [E[2]], [p["participant_id"] for p in r["participants"]], "Cited CAM arrow with G1.")
    classification(r, [E[2]], "unresolved")
    reactions.append(r)

    r = reaction("scheme1_roalmp_generic", [E[1], I[1], E[2], I[2]], "generic_scheme", ["other"],
                 {"pages": [1, 2], "regions": ["graphical abstract, blue pathway", "Scheme 1, This Work"]},
                 "Ring Opening Allene Metathesis Polymerization (ROAlMP)",
                 "ring-opening allene metathesis polymerization",
                 ["This is a generic current-work transformation, not an additional experiment.",
                  "Scheme 1 does not give solvent, temperature, time, conversion, molecular weight, or polymer stereochemistry."])
    r["participants"] = [
        participant("roalmp_generic_1", "m_cyclononadiene_1", "reactant", [E[1], E[2]]),
        participant("roalmp_generic_g2", "m_g2", "catalyst", [E[1], E[2]]),
        participant("roalmp_generic_poly1", "m_poly1", "product", [E[1], E[2]], product_role="polyallene"),
    ]
    add_step(r, "roalmp_generic_arrow", [E[1], E[2]], [p["participant_id"] for p in r["participants"]],
             "Generic G2-catalyzed ROAlMP arrow.")
    classification(r, [E[1], E[2]], "polymer")
    reactions.append(r)

    r = reaction("graphical_g1_vinylidene_generic", [E[1], I[1]], "generic_scheme", ["mechanistic_study"],
                 {"page": 1, "region": "graphical abstract, red Vinylidene Generation pathway"},
                 "Ruthenium Vinylidene Generation", "G1-selective ruthenium vinylidene formation",
                 ["This graphical summary is a generic transformation and duplicates no concrete experimental row.",
                  "Stoichiometry, time, solvent, and temperature are not stated in the graphical abstract."])
    r["participants"] = [
        participant("g1_generic_1", "m_cyclononadiene_1", "reactant", [E[1]]),
        participant("g1_generic_catalyst", "m_g1", "reactant", [E[1]]),
        participant("g1_generic_vinylidene", "m_g1_vinylidene", "product", [E[1]], product_role="ruthenium vinylidene"),
    ]
    add_step(r, "g1_generic_arrow", [E[1]], [p["participant_id"] for p in r["participants"]],
             "Generic G1-selective vinylidene-generation arrow.")
    classification(r, [E[1]], "unresolved")
    reactions.append(r)

    r = reaction("figure1_g2_thf_kinetics", [E[2], I[2], E[3], E[4]], "specific_experiment",
                 ["mechanistic_study", "optimization"],
                 {"page": 2, "region": "Figure 1 and adjacent text", "additional_discussion_pages": [3, 4]},
                 "Ring Opening Allene Metathesis Polymerization (ROAlMP)",
                 "G2-catalyzed ROAlMP kinetic experiment",
                 ["Reaction scale, atmosphere, pressure, workup, isolated yield, polymer Mn, and stereochemistry are not reported for this kinetic trial.",
                  "Residual-carbene percentages refer to an NMR signal population and are not assigned to a fully characterized discrete molecule."])
    r["general_scheme_ids"] = [f"{DOC}:scheme1_roalmp_generic"]
    r["participants"] = [
        participant("fig1_monomer", "m_cyclononadiene_1", "reactant", [E[2]]),
        participant("fig1_g2", "m_g2", "catalyst", [E[2]], mol_percent=("[M]0/[I]0 = 10", 10,
                    "Calculated as 100/10 = 10 mol% G2 relative to monomer 1."), reference="m_cyclononadiene_1"),
        participant("fig1_thf_d8", "m_thf_d8", "solvent", [E[2]]),
        participant("fig1_poly1", "m_poly1", "product", [E[2]], product_role="polyallene"),
    ]
    set_solvent(r, "THF-d8", ["THF-d8"], [E[2]], concentration=("[M]0 = 1 M", 1, "monomer 1", "original_extraction"))
    add_step(r, "fig1_reaction", [E[2]], [p["participant_id"] for p in r["participants"]],
             "Monitor disappearance of monomer and formation of poly(1) by 1H NMR.",
             temperature=("25 °C", {"number": 25, "unit": "°C"}, "original_extraction", None),
             time=("2 h", {"number": 2, "unit": "h"}, "original_extraction", None))
    add_outcome(r, "fig1_conversion", "conversion", ["m_cyclononadiene_1"], [E[2]],
                "92% conversion after 2 h", {"number": 92, "unit": "%"},
                product_role="reactant consumption", basis="1H NMR", scope="monomer conversion after 2 h")
    add_outcome(r, "fig1_kobs", "other", ["m_cyclononadiene_1"], [E[2]],
                "kobs = 0.026 ± 0.004 min−1", "0.026 ± 0.004 min−1",
                basis="first-order fit of ln([M]/[M]0) versus time", scope="observed propagation rate constant")
    add_outcome(r, "fig1_carbene_2h", "other", ["m_g2"], [E[2]],
                "37% residual carbene after 2 h", {"number": 37, "unit": "%"},
                basis="1H NMR carbene signal", scope="residual carbene signal population; discrete species not fully assigned")
    add_outcome(r, "fig1_carbene_17_5min", "other", ["m_g2"], [E[3]],
                "81% after ca. 17.5 min in THF-d8", {"number": 81, "unit": "%", "qualifier": "approx"},
                basis="1H NMR carbene signal", scope="comparative carbene stability")
    r["procedure"]["original_text"] = F(
        "[M]0/[I]0 = 10 with [M]0 = 1 M in THF-d8 at 25 °C; monitored by 1H NMR for 2 h",
        "[M]0/[I]0 = 10; [M]0 = 1 M; THF-d8; 25 °C; 2 h; 1H NMR monitoring", [E[2]])
    classification(r, [E[2]], "polymer")
    reactions.append(r)

    r = reaction("toluene_g2_comparison", [E[2], E[3]], "specific_experiment",
                 ["optimization", "mechanistic_study"],
                 {"page": 3, "region": "first paragraph, toluene-d8 comparison", "inherits_from": "Figure 1 repetition"},
                 "Ring Opening Allene Metathesis Polymerization (ROAlMP)",
                 "G2-catalyzed ROAlMP solvent-comparison experiment",
                 ["[M]/[I], monomer concentration, and 25 °C are inherited because the source states the Figure 1 polymerization was repeated in toluene-d8.",
                  "Reaction scale, atmosphere, pressure, workup, isolated yield, Mn, and polymer stereochemistry are not reported."])
    r["general_scheme_ids"] = [f"{DOC}:scheme1_roalmp_generic"]
    r["participants"] = [
        participant("tol_monomer", "m_cyclononadiene_1", "reactant", [E[2], E[3]]),
        participant("tol_g2", "m_g2", "catalyst", [E[2], E[3]], mol_percent=("repeated [M]0/[I]0 = 10 trial", 10,
                    "The text says the Figure 1 trial was repeated in toluene-d8; 100/10 = 10 mol%."), reference="m_cyclononadiene_1"),
        participant("tol_solvent", "m_toluene_d8", "solvent", [E[3]]),
        participant("tol_poly1", "m_poly1", "product", [E[3]], product_role="polyallene"),
    ]
    set_solvent(r, "repeated in toluene-d8", ["toluene-d8"], [E[3]],
                concentration=("repeated Figure 1 [M]0 = 1 M", 1, "monomer 1", "inference"))
    add_step(r, "tol_reaction", [E[2], E[3]], [p["participant_id"] for p in r["participants"]],
             "Repeat the Figure 1 polymerization in toluene-d8.",
             temperature=("repeated Figure 1 condition, 25 °C", {"number": 25, "unit": "°C"}, "inference",
                          "The page 3 text explicitly calls this a repetition of the Figure 1 polymerization, changing the solvent."),
             time=("2.5 h", {"number": 2.5, "unit": "h"}, "original_extraction", None))
    add_outcome(r, "tol_conversion", "conversion", ["m_cyclononadiene_1"], [E[3]],
                "60% over 2.5 h", {"number": 60, "unit": "%"},
                product_role="reactant consumption", basis="1H NMR monitoring", scope="monomer conversion")
    add_outcome(r, "tol_kobs", "other", ["m_cyclononadiene_1"], [E[3]],
                "kobs = 0.0069 ± 0.004 min−1", "0.0069 ± 0.004 min−1",
                basis="first-order kinetic analysis", scope="observed propagation rate constant")
    add_outcome(r, "tol_carbene_17_5min", "other", ["m_g2"], [E[3]],
                "75% after ca. 17.5 min in toluene-d8", {"number": 75, "unit": "%", "qualifier": "approx"},
                basis="1H NMR carbene signal", scope="comparative carbene stability")
    classification(r, [E[3]], "polymer")
    reactions.append(r)

    rows = [
        (25, 3, 25, 2, 53, 1.6, 8.8, 2.5, "m_g2", "G2"),
        (50, 3, 25, 2, 44, 2.7, 13.7, 1.9, "m_g2", "G2"),
        (75, 3, 25, 2, 25, 2.3, 11.5, 1.9, "m_g2", "G2"),
        (100, 3, 25, 2, 17, 2.1, 11.8, 1.8, "m_g2", "G2"),
        (150, 3, 25, 2, 14, 2.5, 11, 1.8, "m_g2", "G2"),
        (200, 3, 25, 2, 10, 2.5, 11, 1.8, "m_g2", "G2"),
        (50, 2, 25, 2, 49, 3, 9.9, 2.1, "m_g2", "G2"),
        (50, 1, 25, 2, 52, 3.1, 5.8, 2.7, "m_g2", "G2"),
        (50, 0.5, 25, 2, 48, 3, 3.8, 2.8, "m_g2", "G2"),
        (50, 1, 0, 3, 9.1, 0.6, 8.3, 1.9, "m_g2", "G2"),
        (50, 1, 50, 1, 81, 4.9, 4.6, 4.3, "m_g2", "G2"),
        (100, 1, 50, 1, 39, 4.8, 10.8, 1.8, "m_g2", "G2"),
        (250, 1, 50, 1, 23, 7.1, 10.1, 1.8, "m_g2", "G2"),
        (50, 3, 25, 2, 74, 4.5, 6.8, 2.6, "m_hg2", "HG2"),
        (50, 3, 25, 2, 54, 3.3, 11.8, 2.0, "m_g3", "G3"),
    ]
    reactions.extend(table_reaction(i, row, mols) for i, row in enumerate(rows, 1))

    r = reaction("scheme2_mechanistic_paths", [E[4], I[4]], "proposed_mechanism", ["mechanistic_study"],
                 {"page": 4, "region": "Scheme 2, Mechanistic Pathways of ROAlMP"},
                 "Mechanistic Pathways of ROAlMP", "alternative 2,1- and 1,2-insertion pathways",
                 ["Pathway 1 (vinylidene) and Pathway 2 (alkylidene) are author-proposed catalytic cycles, not separate observed net experiments.",
                  "Generic [Ru] intermediates and elementary-step atom mappings are not fully resolved."])
    r["participants"] = [
        participant("scheme2_monomer", "m_cyclononadiene_1", "reactant", [E[4]]),
        participant("scheme2_poly1", "m_poly1", "product", [E[4]], product_role="shared proposed CRU"),
    ]
    r["procedure"]["original_text"] = F(
        "Scheme 2 shows two possible pathways that lead to the same CRU", "two proposed catalytic pathways", [E[4]])
    classification(r, [E[4]], "not_applicable")
    reactions.append(r)

    r = reaction("g1_rt_vinylidene_no_polymerization", [E[4]], "specific_experiment",
                 ["control", "mechanistic_study"],
                 {"page": 4, "region": "right column, G1 [M]0/[I]0 = 10 room-temperature trial"},
                 "Ruthenium Vinylidene Formation", "G1 reaction with 1 at room temperature",
                 ["The source says no polymerization after reaction times similar to G2 trials but does not state one exact total duration.",
                  "No isolated yield or complete structure of the vinylidene is reported in the inspected pages."])
    r["general_scheme_ids"] = [f"{DOC}:graphical_g1_vinylidene_generic"]
    r["participants"] = [
        participant("g1rt_monomer", "m_cyclononadiene_1", "reactant", [E[4]]),
        participant("g1rt_g1", "m_g1", "reactant", [E[4]], mol_percent=("[M]0/[I]0 = 10", 10,
                    "Calculated as 100/10 = 10 mol% G1 relative to monomer 1."), reference="m_cyclononadiene_1"),
        participant("g1rt_tol", "m_toluene_d8", "solvent", [E[4]]),
        participant("g1rt_vinylidene", "m_g1_vinylidene", "product", [E[4]], product_role="ruthenium vinylidene"),
    ]
    set_solvent(r, "1 M in toluene-d8", ["toluene-d8"], [E[4]],
                concentration=("[M]0 = 1 M", 1, "monomer 1", "original_extraction"))
    add_step(r, "g1rt_formation", [E[4]], [p["participant_id"] for p in r["participants"]],
             "Observe rapid disappearance of the G1 alkylidene signal and formation of a vinylidene-like species.",
             temperature=("room temperature", "room temperature", "original_extraction", None),
             time=("5 min to complete alkylidene-signal disappearance", {"number": 5, "unit": "min"}, "original_extraction", None))
    add_outcome(r, "g1rt_no_polymerization", "no_reaction", [], [E[4]],
                "no polymerization had occurred", "no polymerization observed",
                basis="1H NMR observation", scope="homopolymerization of 1")
    add_outcome(r, "g1rt_vinylidene_observed", "other", ["m_g1_vinylidene"], [E[4]],
                "new 13C NMR peak at δ 336 ppm", "δ 336 ppm",
                product_role="ruthenium vinylidene", basis="13C NMR; comparison to reported first-generation ruthenium vinylidenes",
                scope="vinylidene assignment")
    classification(r, [E[4]], "unresolved")
    reactions.append(r)

    r = reaction("g1_equimolar_vinylidene_hrms", [E[4]], "specific_experiment",
                 ["mechanistic_study"],
                 {"page": 4, "region": "right column, [M]0/[I]0 = 1 ESI-HRMS trial"},
                 "Ruthenium Vinylidene Formation", "equimolar G1/1 vinylidene-generation trial",
                 ["Temperature, time, solvent, scale, and isolation details are not restated for the [M]/[I] = 1 trial.",
                  "The ESI-HRMS peak is characterization evidence and is retained in source_context rather than treated as a synthetic yield."])
    r["general_scheme_ids"] = [f"{DOC}:graphical_g1_vinylidene_generic"]
    r["participants"] = [
        participant("g1eq_monomer", "m_cyclononadiene_1", "reactant", [E[4]]),
        participant("g1eq_g1", "m_g1", "reactant", [E[4]], mol_percent=("[M]0/[I]0 = 1", 100,
                    "Calculated as 100/1 = 100 mol% G1 relative to monomer 1."), reference="m_cyclononadiene_1"),
        participant("g1eq_product", "m_g1_vinylidene", "product", [E[4]], product_role="ruthenium vinylidene"),
    ]
    add_step(r, "g1eq_formation", [E[4]], [p["participant_id"] for p in r["participants"]],
             "Generate G1-derived ruthenium vinylidene in an [M]0/[I]0 = 1 trial and analyze by ESI-HRMS.")
    add_outcome(r, "g1eq_vinylidene_assignment", "other", ["m_g1_vinylidene"], [E[4]],
                "[M − Cl]+• at m/z 909.4908 (909.4937 calcd)",
                {"observed_m_z": 909.4908, "calculated_m_z": 909.4937, "ion": "[M-Cl]+ radical"},
                product_role="ruthenium vinylidene", basis="ESI-HRMS", scope="species assignment")
    classification(r, [E[4]], "unresolved")
    reactions.append(r)

    r = reaction("g1_50c_unassigned_polymer", [E[4], E[5]], "specific_experiment",
                 ["mechanistic_study", "control"],
                 {"pages": [4, 5], "region": "G1 [M]0/[I]0 = 10 trial heated to 50 °C and continuation"},
                 "G1 Trial at Elevated Temperature", "G1/vinylidene-mediated consumption of 1 to an unassigned polymer",
                 ["Room-temperature vinylidene-formation duration is not restated for this repeated trial.",
                  "No yield, molecular weight, or correct polymer microstructure is assigned; classification remains unresolved.",
                  "The source reports new 1H NMR peaks before precipitation and no allene signal in the isolated material."])
    r["participants"] = [
        participant("g150_monomer", "m_cyclononadiene_1", "reactant", [E[4], E[5]]),
        participant("g150_g1", "m_g1", "catalyst", [E[4]], mol_percent=("[M]0/[I]0 = 10", 10,
                    "Calculated as 100/10 = 10 mol% G1 relative to monomer 1."), reference="m_cyclononadiene_1"),
        participant("g150_tol", "m_toluene_d8", "solvent", [E[4]]),
        participant("g150_product", "m_g1_unknown_polymer", "product", [E[4], E[5]], product_role="unassigned polymer"),
        participant("g150_meoh", "m_methanol", "purification_agent", [E[5]], stage="purification"),
    ]
    set_solvent(r, "repeated G1 trial in toluene-d8", ["toluene-d8"], [E[4]],
                concentration=("repeated [M]0 = 1 M G1 trial", 1, "monomer 1", "inference"))
    add_step(r, "g150_vinylidene_formation", [E[4]], ["g150_monomer", "g150_g1", "g150_tol"],
             "Allow complete vinylidene formation at room temperature.",
             temperature=("room temperature", "room temperature", "original_extraction", None))
    add_step(r, "g150_heating", [E[4]], ["g150_monomer", "g150_g1", "g150_tol", "g150_product"],
             "Heat after complete vinylidene formation and monitor by 1H NMR.",
             temperature=("50 °C", {"number": 50, "unit": "°C"}, "original_extraction", None),
             time=("over an hour", {"number": 1, "unit": "h", "qualifier": "approx"}, "original_extraction", None))
    add_step(r, "g150_precipitation", [E[5]], ["g150_product", "g150_meoh"],
             "Triple precipitation into methanol at −78 °C.",
             temperature=("−78 °C", {"number": -78, "unit": "°C"}, "original_extraction", None), stage="purification")
    add_outcome(r, "g150_consumption", "conversion", ["m_cyclononadiene_1"], [E[4]],
                "Disappearance of monomer over an hour", "monomer disappeared qualitatively",
                product_role="reactant consumption", basis="1H NMR", scope="qualitative monomer consumption")
    add_outcome(r, "g150_unassigned_product", "other", ["m_g1_unknown_polymer"], [E[5]],
                "No allene signal was observed in 13C NMR; correct microstructure could not be assigned",
                "isolated polymer microstructure unassigned; no allene signal observed",
                product_role="unassigned polymer", basis="13C NMR of isolated material", scope="microstructure observation")
    classification(r, [E[4], E[5]], "unresolved")
    reactions.append(r)

    r = reaction("figure3_nbe_ipr_polymerization", [E[5], I[5]], "specific_experiment",
                 ["mechanistic_study"],
                 {"page": 5, "region": "Figure 3 and adjacent text"},
                 "Ring-Opening Metathesis Polymerization of NBE-iPr Using a Ruthenium Vinylidene",
                 "vinylidene-mediated ring-opening metathesis polymerization of NBE-iPr",
                 ["The prose gives a slight excess of 1 while Figure 3 displays [M]/[I] = 1; the exact excess and whether the displayed ratio is rounded are unresolved.",
                  "Near-quantitative is a qualitative conversion statement, not an isolated yield.",
                  "The source calls the process uncontrolled and does not assign complete polymer stereochemistry or chain ends."])
    r["participants"] = [
        participant("nbe_stage1_1", "m_cyclononadiene_1", "reactant", [E[5], I[5]], reference="m_g1"),
        participant("nbe_stage1_g1", "m_g1", "reactant", [E[5], I[5]],
                    mol_percent=("Figure 3: [M]/[I] = 1", 100,
                                 "Calculated as 100/1 = 100 mol% G1 relative to monomer 1 from the displayed scheme ratio; the prose separately says a slight excess of 1."),
                    reference="m_cyclononadiene_1"),
        participant("nbe_vinylidene", "m_g1_vinylidene", "product", [E[5]], product_role="in situ catalytic species"),
        participant("nbe_monomer", "m_nbe_ipr", "reactant", [E[5]],
                    equivalents=("10 equiv of NBE-iPr", 10, "original_extraction", None), reference="m_g1"),
        participant("nbe_solvent", "m_toluene_d8", "solvent", [E[5]]),
        participant("nbe_polymer", "m_poly_nbe_ipr", "product", [E[5]], product_role="polymer product"),
    ]
    ratio_note = F(None, None, [E[5], I[5]], status="conflict")
    ratio_note["alternatives"] = [
        F("only slight excess of 1", "slight excess relative to G1", [E[5]]),
        F("[M]/[I] = 1", {"monomer_to_initiator_ratio": 1}, [I[5]]),
    ]
    r["notes"].append(ratio_note)
    set_solvent(r, "Figure 3 spectra in toluene-d8", ["toluene-d8"], [E[5]])
    add_step(r, "nbe_vinylidene_generation", [E[5]], ["nbe_stage1_1", "nbe_stage1_g1", "nbe_solvent", "nbe_vinylidene"],
             "Generate ruthenium vinylidene from G1 with slight excess of 1; scheme labels [M]/[I] = 1.",
             temperature=("room temperature", "room temperature", "inherited_general_condition",
                          None))
    add_step(r, "nbe_room_temperature_hold", [E[5]], ["nbe_vinylidene", "nbe_monomer", "nbe_solvent"],
             "Add 10 equiv NBE-iPr and monitor at 25 °C; no polymerization is observed even after extended time.",
             temperature=("25 °C", {"number": 25, "unit": "°C"}, "original_extraction", None))
    add_step(r, "nbe_heating", [E[5]], ["nbe_vinylidene", "nbe_monomer", "nbe_solvent", "nbe_polymer"],
             "Heat the reaction mixture to form poly(NBE-iPr).",
             temperature=("80 °C", {"number": 80, "unit": "°C"}, "original_extraction", None),
             time=("2.5 h", {"number": 2.5, "unit": "h"}, "original_extraction", None))
    add_outcome(r, "nbe_no_rt_polymerization", "no_reaction", [], [E[5]],
                "no polymerization of NBE-iPr was observed at room temperature", "no polymerization observed at 25 °C",
                basis="1H NMR monitoring", scope="room-temperature hold")
    add_outcome(r, "nbe_conversion", "conversion", ["m_nbe_ipr"], [E[5]],
                "near quantitative conversion", {"number": 100, "unit": "%", "qualifier": "approx"},
                product_role="reactant consumption", basis="1H NMR monitoring", scope="NBE-iPr conversion after 2.5 h at 80 °C")
    add_outcome(r, "nbe_dp", "other", ["m_poly_nbe_ipr"], [E[5]],
                "DPNMR = 47", 47, product_role="polymer product", basis="NMR", scope="degree of polymerization")
    add_outcome(r, "nbe_dispersity", "other", ["m_poly_nbe_ipr"], [E[5]],
                "Đ = 1.6", 1.6, product_role="polymer product", scope="dispersity; main-text measurement method not stated")
    classification(r, [E[5]], "polymer", nbe=True)
    reactions.append(r)

    bundle["reactions"] = reactions
    bundle["source_context"] = [
        {"context_id": f"{DOC}:article_metadata", "document_id": DOC, "evidence_ids": [E[1]],
         "kind": "bibliographic_metadata", "title": doc["title"], "doi": doc["doi"],
         "note": "ACS Macro Letters 2021, 10, 642–648; inspected source pages 1–6 only."},
        {"context_id": f"{DOC}:ring_strain_context", "document_id": DOC, "evidence_ids": [E[2]],
         "kind": "computational_context", "observations": ["1: 12.7 kcal/mol ring strain", "1,2-cyclodecadiene: 15.9 kcal/mol"],
         "note": "Isodesmic calculations are mechanistic/material context, not synthetic reaction records."},
        {"context_id": f"{DOC}:figure2_poly1_characterization", "document_id": DOC, "evidence_ids": [E[3], I[3], E[4]],
         "kind": "characterization", "related_reaction_ids": [f"{DOC}:figure1_g2_thf_kinetics"],
         "observations": ["Poly(1) precipitated into MeOH at −78 °C and vacuum-dried before characterization.",
                          "13C NMR central allene carbon δC 204; monomer 1 δC 206.",
                          "MALDI-TOF repeating mass 122.2 g mol−1."],
         "note": "NMR and MALDI spectra are retained as characterization, not additional syntheses."},
        {"context_id": f"{DOC}:mechanistic_spectroscopy", "document_id": DOC, "evidence_ids": [E[4], I[4]],
         "kind": "mechanistic_characterization",
         "related_reaction_ids": [f"{DOC}:figure1_g2_thf_kinetics", f"{DOC}:g1_rt_vinylidene_no_polymerization",
                                  f"{DOC}:g1_equimolar_vinylidene_hrms"],
         "observations": ["G2 [M]/[I] = 10: alkylidene peak δC 295 plus monomer/polymer allene peaks.",
                          "G1 [M]/[I] = 10: new peak δ 336 ppm.",
                          "G1 [M]/[I] = 1: ESI-HRMS m/z 909.4908, 909.4937 calculated for [M−Cl]+ radical."],
         "note": "Spectral values support species assignments but are not independent synthetic records."},
        {"context_id": f"{DOC}:table1_narrative_precision", "document_id": DOC, "evidence_ids": [E[3], E[4]],
         "kind": "cross_region_table_discussion",
         "related_reaction_ids": [f"{DOC}:table1_entry01", f"{DOC}:table1_entry02", f"{DOC}:table1_entry06",
                                  f"{DOC}:table1_entry10", f"{DOC}:table1_entry11", f"{DOC}:table1_entry14",
                                  f"{DOC}:table1_entry15"],
         "observations": ["The prose gives 52.6% and 10.2% for the Table 1 catalyst-loading endpoints, while the table displays 53% and 10%.",
                          "The prose compares 44% and 81% at [M]/[I] = 50 and describes HG2 74% versus G2 44%."],
         "note": "The endpoint differences are retained as higher-precision narrative values consistent with table rounding, not silently substituted for the table cells and not treated as a chemical conflict."},
        {"context_id": f"{DOC}:poly1_material_properties", "document_id": DOC, "evidence_ids": [E[5], I[5]],
         "kind": "material_characterization", "related_molecule_ids": ["m_poly1"],
         "observations": ["TGA under nitrogen: 5% mass loss at 390 °C.",
                          "DSC −60 to 100 °C at 10 and 20 °C min−1: no Tg observed.",
                          "THF, 1.0 mg mL−1: broad emission λEM 360–660 nm."],
         "note": "Thermal and photophysical measurements are not synthesis experiments."},
        {"context_id": f"{DOC}:page6_inventory", "document_id": DOC, "evidence_ids": [E[6], I[6]],
         "kind": "nonreaction_page", "note": "Page 6 contains acknowledgments and references only; cited papers were not promoted to current-paper experiments."},
        {"context_id": f"{DOC}:si_status", "document_id": DOC, "evidence_ids": [E[5]],
         "kind": "source_limitation", "note": "The article states SI is available, but no SI file is explicitly associated in source_mapping.json; no SI facts were imported."},
    ]

    raw_inventory = {
        "document_id": DOC,
        "page_range": [1, 6],
        "semantic_accuracy_verified": False,
        "note": "Raw model inventory retained separately from normalized records; all items remain needs_review.",
        "pages": [
            {"page": 1, "regions": ["title/abstract", "graphical abstract"],
             "findings": ["Two generic current-work pathways: G2→poly(1), G1→ruthenium vinylidene.",
                          "Introduction/citations are context, not experiments."],
             "record_ids": [f"{DOC}:scheme1_roalmp_generic", f"{DOC}:graphical_g1_vinylidene_generic"]},
            {"page": 2, "regions": ["Scheme 1", "Figure 1", "adjacent text"],
             "findings": ["Two cited previous-work schemes, one current generic ROAlMP scheme, one concrete G2/THF-d8 kinetic trial."],
             "record_ids": [f"{DOC}:scheme1_rcam_prior", f"{DOC}:scheme1_cam_prior",
                            f"{DOC}:scheme1_roalmp_generic", f"{DOC}:figure1_g2_thf_kinetics"]},
            {"page": 3, "regions": ["Figure 2", "Table 1", "solvent-comparison prose"],
             "findings": ["Figure 2 is characterization only.", "Table 1 has 15 independent optimization rows.",
                          "Toluene-d8 repetition is a distinct concrete experiment."],
             "table1_rows": [list(row[:8]) + [row[9]] for row in rows],
             "record_ids": [f"{DOC}:toluene_g2_comparison"] + [f"{DOC}:table1_entry{i:02d}" for i in range(1, 16)]},
            {"page": 4, "regions": ["condition-screen discussion", "Scheme 2", "G1/G2 mechanistic prose"],
             "findings": ["Table 1 discussion adds no new rows.", "Scheme 2 is a proposed mechanism.",
                          "Concrete G1 room-temperature and [M]/[I]=1 trials are distinct; 50 °C trial begins here."],
             "record_ids": [f"{DOC}:table1_entry{i:02d}" for i in range(1, 16)] +
                           [f"{DOC}:scheme2_mechanistic_paths", f"{DOC}:g1_rt_vinylidene_no_polymerization",
                            f"{DOC}:g1_equimolar_vinylidene_hrms", f"{DOC}:g1_50c_unassigned_polymer"]},
            {"page": 5, "regions": ["G1 50 °C continuation", "Figure 3", "material-property prose", "SI statement"],
             "findings": ["G1 50 °C trial yields an unassigned polymer.", "Figure 3 is one sequential NBE-iPr polymerization experiment.",
                          "TGA/DSC/fluorescence are characterization only; SI is mentioned but not supplied."],
             "record_ids": [f"{DOC}:g1_50c_unassigned_polymer", f"{DOC}:figure3_nbe_ipr_polymerization"]},
            {"page": 6, "regions": ["acknowledgments", "references 1–37"],
             "findings": ["No current-paper reaction experiment on this page."], "record_ids": []},
        ],
    }
    st.write(OUT / "raw_model/pages_0001_0006_extraction.json", raw_inventory)

    coverage = {
        "document_id": DOC,
        "semantic_accuracy_verified": False,
        "pages": [
            {"page": 1, "state": "attempted",
             "regions": [
                 {"region": "title and abstract", "finding": "Bibliographic metadata and current-work summary; no additional concrete experiment row.", "evidence_ids": [E[1], I[1]]},
                 {"region": "graphical abstract", "finding": "Generic G2 ROAlMP and G1 vinylidene-generation pathways; mapped to generic records.", "evidence_ids": [E[1], I[1]]}],
             "record_ids": [f"{DOC}:scheme1_roalmp_generic", f"{DOC}:graphical_g1_vinylidene_generic"],
             "unresolved": ["Graphical abstract gives no stoichiometry, time, solvent, or temperature."]},
            {"page": 2, "state": "attempted",
             "regions": [
                 {"region": "Scheme 1", "finding": "Two cited literature examples and one current-work generic ROAlMP arrow.", "evidence_ids": [E[2], I[2]]},
                 {"region": "Figure 1 and right-column text", "finding": "Concrete G2/THF-d8 kinetic polymerization with conversion, carbene signal, and kobs.", "evidence_ids": [E[2], I[2]]}],
             "record_ids": [f"{DOC}:scheme1_rcam_prior", f"{DOC}:scheme1_cam_prior", f"{DOC}:scheme1_roalmp_generic",
                            f"{DOC}:figure1_g2_thf_kinetics"],
             "unresolved": ["Scheme 1 prior-work structures/conditions are generic and incomplete.",
                            "Figure 1 does not report scale, atmosphere, pressure, workup, isolated yield, Mn, or polymer stereochemistry."]},
            {"page": 3, "state": "attempted",
             "regions": [
                 {"region": "Figure 2", "finding": "Poly(1) NMR/MALDI characterization; not a separate synthesis.", "evidence_ids": [E[3], I[3]]},
                 {"region": "Table 1 rows 1–15 and footnotes a–e", "finding": "Fifteen distinct optimization experiments with conversion, theoretical Mn, SEC Mn, and dispersity.", "evidence_ids": [E[3], I[3]]},
                 {"region": "lower-left solvent paragraph", "finding": "Distinct toluene-d8 repetition with 60% conversion, kobs, and residual-carbene comparison.", "evidence_ids": [E[3], I[3]]}],
             "record_ids": [f"{DOC}:toluene_g2_comparison"] + [f"{DOC}:table1_entry{i:02d}" for i in range(1, 16)],
             "unresolved": ["Table 1 omits scale, atmosphere, pressure, workup, isolation yield, and polymer stereochemistry.",
                            "THF is inherited from adjacent prose rather than repeated in each row."]},
            {"page": 4, "state": "attempted",
             "regions": [
                 {"region": "left-column Table 1 discussion", "finding": "Narrative interpretation of already-recorded Table 1 rows; no new experiment row.", "evidence_ids": [E[4], I[4]]},
                 {"region": "Scheme 2", "finding": "Author-proposed vinylidene and alkylidene catalytic pathways; one proposed-mechanism record.", "evidence_ids": [E[4], I[4]]},
                 {"region": "right-column G1/G2 studies", "finding": "G1 room-temperature no-polymerization/vinylidene trial, distinct equimolar HRMS trial, and start of the repeated 50 °C trial.", "evidence_ids": [E[4], I[4]]}],
             "record_ids": [f"{DOC}:figure1_g2_thf_kinetics", f"{DOC}:scheme2_mechanistic_paths",
                            *[f"{DOC}:table1_entry{i:02d}" for i in range(1, 16)],
                            f"{DOC}:g1_rt_vinylidene_no_polymerization", f"{DOC}:g1_equimolar_vinylidene_hrms",
                            f"{DOC}:g1_50c_unassigned_polymer"],
             "unresolved": ["Generic [Ru] intermediates are not safely encoded or atom-mapped.",
                            "The equimolar HRMS trial does not restate solvent, temperature, time, or scale."]},
            {"page": 5, "state": "attempted",
             "regions": [
                 {"region": "upper-left continuation", "finding": "50 °C G1 trial workup and unassigned polymer observation.", "evidence_ids": [E[5], I[5]]},
                 {"region": "Figure 3 and lower-left prose", "finding": "Sequential G1→vinylidene→NBE-iPr experiment; room-temperature negative stage and 80 °C polymerization stage.", "evidence_ids": [E[5], I[5]]},
                 {"region": "right-column material properties", "finding": "TGA, DSC, and fluorescence characterization only; stored in source_context.", "evidence_ids": [E[5], I[5]]},
                 {"region": "Associated Content", "finding": "SI availability is stated, but no SI file is explicitly associated in this workspace.", "evidence_ids": [E[5], I[5]]}],
             "record_ids": [f"{DOC}:g1_50c_unassigned_polymer", f"{DOC}:figure3_nbe_ipr_polymerization"],
             "unresolved": ["G1-derived 50 °C polymer microstructure remains author-unassigned.",
                            "NBE-iPr exact structure/stereochemistry and polymer cis/trans distribution are not encoded.",
                            "Figure 3 prose says slight excess of 1 while the scheme displays [M]/[I] = 1; exact excess/rounding is unresolved, and the duration of the 25 °C hold is not quantified."]},
            {"page": 6, "state": "attempted",
             "regions": [
                 {"region": "acknowledgments", "finding": "Funding and acknowledgments; no reaction experiment.", "evidence_ids": [E[6], I[6]]},
                 {"region": "references 1–37", "finding": "Bibliographic references only; cited publications were not promoted to current-paper experiments.", "evidence_ids": [E[6], I[6]]}],
             "record_ids": [],
             "unresolved": []},
        ],
    }
    st.write(OUT / "coverage/pages_0001_0006.json", coverage)
    st.write(OUT / "submitted.json", bundle)
    print(json.dumps({"molecules": len(mols), "reactions": len(reactions),
                      "page_range": [1, 6]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
