#!/usr/bin/env python3
"""Build the bounded pages 1-6 v3.1 extraction for doc_403ed550996f2a54."""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st


DOC = "doc_403ed550996f2a54"
RUN = ROOT / "runs/current"
DATASET = st.read(RUN / "dataset.json")
EV = {p: f"{DOC}:page{p}:text" for p in range(1, 7)}
IMG = {p: f"{DOC}:page{p}:image" for p in range(1, 7)}


def f(raw=None, normalized=None, evidence=(), method="original_extraction", derivation=None,
      status=None, verification=None, alternatives=()):
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation,
                   status=status, verification=verification, alternatives=alternatives)


def derived(normalized, evidence, derivation):
    return f(None, normalized, evidence, method="structure_derivation", derivation=derivation)


def inferred(normalized, evidence, derivation):
    return f(None, normalized, evidence, method="inference", derivation=derivation)


def set_name(m, name, evidence, label=None):
    m["identity_status"] = "source_identified"
    m["names"]["original"] = f(name, name, evidence)
    m["names"]["english"] = f(name, name, evidence)
    m["occurrences"][0]["compound_label"] = f(label, label, evidence) if label else f()
    m["occurrences"][0]["evidence_ids"] = list(evidence)
    return m


def statuses(m, evidence, chirality, geometry, representation, *, atom=None, double=None,
             axial=None, mixture=None, assignment=None, reason="Source structure/stereo review."):
    for key, value in (("chirality_status", chirality), ("geometry_status", geometry),
                       ("representation_status", representation)):
        m["stereochemistry"][key] = inferred(value, evidence, reason)
    for key, value in (("atom_centers", atom), ("double_bonds", double), ("axial", axial),
                       ("mixture", mixture), ("assignment", assignment)):
        if value is not None:
            m["stereochemistry"][key] = f(value, value, evidence)
    return m


def unresolved_molecule(mid, name, label, evidence, *, chirality="unresolved",
                        geometry="unresolved", representation="unresolved", atom=None,
                        double=None, axial=None, mixture=None, assignment=None, note=None):
    m = set_name(st.new_molecule(mid, DOC, evidence), name, evidence, label)
    statuses(m, evidence, chirality, geometry, representation, atom=atom, double=double,
             axial=axial, mixture=mixture, assignment=assignment,
             reason="The bounded source was inspected, but a complete stereochemical graph was not safely encoded.")
    if note:
        m["identity_note"] = f(note, note, evidence)
    return m


def smiles_molecule(mid, name, label, smiles, evidence, *, chirality="achiral",
                    geometry="not_applicable", representation="encoded", axial=None,
                    note=None, allene=False):
    m = set_name(st.new_molecule(mid, DOC, evidence), name, evidence, label)
    diag = st.stereo_inspect(smiles)
    canonical = diag["canonical_isomeric_smiles"]
    m["representation_type"] = "smiles"
    m["smiles"] = derived(canonical, evidence,
                          "Candidate connectivity transcribed from the source label/drawing and canonicalized with RDKit; parsing does not verify source identity.")
    m["stereochemistry"]["isomeric_smiles"] = derived(canonical, evidence,
                          "RDKit canonical isomeric SMILES for the candidate source graph.")
    m["stereochemistry"]["identity_key"] = derived(diag["stereo_graph_key"], evidence,
                          "SHA-256 stereo graph key returned by standard stereo diagnostics.")
    statuses(m, evidence, chirality, geometry, representation, axial=axial,
             assignment={"tool_diagnostics": diag["potential_stereo"],
                         "source_stereo_assertion": "none beyond the displayed connectivity"},
             reason="Status combines the inspected source drawing with standard stereo graph diagnostics.")
    if note:
        m["identity_note"] = f(note, note, evidence)
    if allene:
        m["allene_features"] = [{
            "fragment_id": "allene_1",
            "carbon_atom_maps": derived([0, 1, 2], evidence,
                "Zero-based atom indices refer to the stored canonical candidate SMILES; verify against the drawing before reuse."),
            "substitution": derived("terminal allene; one terminus is CH2", evidence,
                "Derived from the source generic structure R-CH=C=CH2 and the concrete R assignment."),
            "axial_stereochemistry": inferred("not stereogenic at the CH2 terminus", evidence,
                "A terminal allene bearing two identical hydrogens at one terminus is not configurationally axial-chiral."),
            "local_environment": derived(name, evidence, "Source R-group assignment expanded into a candidate connectivity."),
        }]
    return m


def generic_source_molecule(mid, name, label, evidence, note):
    return unresolved_molecule(mid, name, label, evidence, chirality="unresolved",
                               geometry="unresolved", representation="not_encoded", note=note)


def charge(p, key, raw, normalized, evidence, method="original_extraction", derivation=None):
    p["charges"][key] = f(raw, normalized, evidence, method=method, derivation=derivation)


def participant(pid, mid, role, evidence, stage="reaction", product_role=None):
    p = st.new_participant(pid, mid, role, evidence, stage)
    if product_role is not None:
        p["product_role"] = f(product_role, product_role, evidence)
    return p


def step(step_id, evidence, stage, pids, operation, temperature=None, time=None,
         atmosphere=None, pressure=None, method="inherited_general_condition"):
    s = st.new_step(step_id, evidence, stage)
    s["participant_ids"] = list(pids)
    s["operation"] = f(operation, operation, evidence, method=method)
    if temperature is not None:
        s["temperature"] = f(temperature[0], temperature[1], evidence, method=method)
    if time is not None:
        s["time"] = f(time[0], time[1], evidence, method=method)
    if atmosphere is not None:
        s["atmosphere"] = f(atmosphere, atmosphere, evidence, method=method)
    if pressure is not None:
        s["pressure"] = f(pressure, pressure, evidence, method=method)
    return s


def partial_cycloaddition_classification(r, evidence):
    changes = [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
         "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
         "before_order": 0, "after_order": 1, "change_type": "formed"},
    ]
    r["classification"]["bond_changes"] = [derived(x, evidence,
        "The displayed [6π+2π] ring closure forms two new C-C sigma connections; exact atom mapping and accompanying bond-order changes remain unresolved.") for x in changes]
    r["classification"]["labels"] = derived([
        {"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1", "b2"]}
    ], evidence, "Taxonomy label derived only from the two visually apparent new C-C connections.")
    r["classification"]["coverage"] = inferred("partial", evidence,
        "Reactant/product atom mapping, implicit-H changes, and all bond-order changes were not established safely from the scheme.")
    r["classification"]["limitations"] = [
        "Two new C-C sigma bonds are represented at scheme level; their exact atom maps are unresolved.",
        "No complete atom mapping or hydrogen balance is claimed.",
        "Unchanged product bonds and bond-order decreases are not labeled as newly formed bonds.",
    ]


def set_reaction_header(r, record_type, categories, location, name, transform, evidence):
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, evidence)
    r["location"] = f(location, location, evidence)
    r["reaction_name"] = f(name, name, evidence)
    r["transformation_type"] = inferred(transform, evidence,
        "Transformation type normalized from the source reaction name and depicted bond construction.")


def add_common_specific_participants(r, substrate_mid, product_mid, table_ev):
    e34 = [table_ev, EV[4]]
    ps = [
        participant("p_azepine_1", "m_azepine_1", "reactant", e34),
        participant("p_substrate", substrate_mid, "reactant", e34),
        participant("p_product", product_mid, "product", e34, product_role="target product"),
        participant("p_co", "m_co_acac2_dppe", "catalyst", e34),
        participant("p_zn", "m_zinc", "reagent", e34),
        participant("p_zni2", "m_zni2", "additive", e34),
        participant("p_dce", "m_dce", "solvent", e34),
        participant("p_petether_workup", "m_petroleum_ether", "workup_agent", [EV[4]], "workup"),
        participant("p_silica_filter", "m_silica_gel", "workup_agent", [EV[4]], "workup"),
        participant("p_silica_column", "m_silica_gel", "purification_agent", [EV[4]], "purification"),
        participant("p_petether_eluent", "m_petroleum_ether", "purification_agent", [EV[4]], "purification"),
        participant("p_etoac", "m_ethyl_acetate", "purification_agent", [EV[4]], "purification"),
    ]
    for p in ps:
        r["participants"].append(p)
    byid = {p["participant_id"]: p for p in r["participants"]}
    charge(byid["p_azepine_1"], "amount", "1 (1 mmol)", {"number": 1.0, "unit": "mmol"}, e34)
    charge(byid["p_substrate"], "amount", "2 or 4 (1.5 mmol)", {"number": 1.5, "unit": "mmol"}, e34)
    charge(byid["p_co"], "amount", "Co(acac)2(dppe) (0.10 mmol)", {"number": 0.10, "unit": "mmol"}, e34)
    charge(byid["p_co"], "mol_percent", "10 mol % Co(acac)2(dppe)", {"number": 10.0, "unit": "mol%"}, e34)
    byid["p_co"]["equivalence_reference"] = f("relative to 1", "m_azepine_1", e34)
    charge(byid["p_zn"], "amount", "Zn (0.3 mmol)", {"number": 0.3, "unit": "mmol"}, e34)
    charge(byid["p_zn"], "mol_percent", "30 mol % Zn", {"number": 30.0, "unit": "mol%"}, e34)
    byid["p_zn"]["equivalence_reference"] = f("relative to 1", "m_azepine_1", e34)
    charge(byid["p_zni2"], "amount", "ZnI2 (0.20 mmol)", {"number": 0.20, "unit": "mmol"}, e34)
    charge(byid["p_zni2"], "mol_percent", "20 mol % ZnI2", {"number": 20.0, "unit": "mol%"}, e34)
    byid["p_zni2"]["equivalence_reference"] = f("relative to 1", "m_azepine_1", e34)
    charge(byid["p_dce"], "volume", "DCE (3 mL); procedure: 1.5 mL + 1.5 mL", {"number": 3.0, "unit": "mL"}, e34)
    return byid


def add_common_specific_conditions(r, table_ev):
    e4 = [EV[4]]
    all_reaction = ["p_azepine_1", "p_substrate", "p_product", "p_co", "p_zn", "p_zni2", "p_dce"]
    r["reaction_scale"] = f("N-carbocholesteroxyazepine 1 (1.0 mmol)",
                            {"reference_molecule_id": "m_azepine_1", "amount": {"number": 1.0, "unit": "mmol"}},
                            [table_ev, EV[4]], method="inherited_general_condition")
    r["solvent_conditions"]["composition"] = f("DCE (1,2-dichloroethane)", ["m_dce"], [table_ev, EV[4]], method="inherited_general_condition")
    r["solvent_conditions"]["total_volume"] = f("DCE (3 mL)", {"number": 3.0, "unit": "mL"}, [table_ev, EV[4]], method="inherited_general_condition")
    r["steps"] = [
        step("s1_catalyst_preactivation", e4, "reaction", ["p_co", "p_zn", "p_dce"],
             "Zinc powder was added to Co(acac)2(dppe) in DCE (1.5 mL) and stirred.",
             temperature=("room temperature", "room temperature"),
             time=("2 min", {"number": 2.0, "unit": "min"}), atmosphere="dry argon"),
        step("s2_successive_addition", e4, "reaction", ["p_azepine_1", "p_substrate", "p_zni2", "p_dce"],
             "N-carbocholesteroxyazepine 1, substrate in DCE (1.5 mL), and dry ZnI2 were added successively; ZnI2 was weighed in a dry weighing bottle.",
             atmosphere="dry argon"),
        step("s3_reaction_hold", [table_ev, EV[4]], "reaction", all_reaction,
             "The reaction mixture was heated.", temperature=("60 °C", {"number": 60.0, "unit": "°C"}),
             time=("20 h", {"number": 20.0, "unit": "h"}), atmosphere="dry argon"),
        step("s4_catalyst_deactivation", e4, "workup", ["p_petether_workup"],
             "Petroleum ether was added and the mixture was stirred in air to deactivate the catalyst.",
             time=("10 min", {"number": 10.0, "unit": "min"}), atmosphere="air"),
        step("s5_filtration_evaporation", e4, "workup", ["p_silica_filter"],
             "The mixture was filtered through a short pad of silica and volatiles were removed under vacuum."),
        step("s6_column_chromatography", e4, "purification", ["p_silica_column", "p_petether_eluent", "p_etoac", "p_product"],
             "Chromatographic purification over silica gel with petroleum ether/ethyl acetate 5:1 afforded the target product."),
    ]
    r["procedure"]["original_text"] = f(
        "Zinc powder (30 mol %) was added to Co(acac)2(dppe) (10 mol %) in DCE (1.5 mL) under dry argon and stirred at room temperature for 2 min; 1 (1.0 mmol), 1,2-diene or 1,3-diyne (1.5 mmol) in DCE (1.5 mL), and dry ZnI2 (20 mol %) were added successively; heat 60 °C for 20 h; petroleum ether/air 10 min; silica filtration; evaporation; silica chromatography (petroleum ether/ethyl acetate 5:1).",
        "Cycloaddition general reaction procedure on page 4", e4, method="inherited_general_condition")
    r["procedure"]["general_procedure_reference"] = f("Cycloaddition of 1,2-Dienes and 1,3-Diynes to N-Carbocholesteroxyazepine (General Reaction Procedure)",
                                                        "page 4 general reaction procedure", e4)
    r["procedure"]["addition_order"] = f("Zn to Co complex/DCE; then 1; then substrate in DCE; then dry ZnI2", ["Zn", "Co(acac)2(dppe)/DCE", "1", "substrate/DCE", "ZnI2"], e4)
    r["procedure"]["pretreatment"] = f("Co(acac)2(dppe)/Zn stirred for 2 min at room temperature under dry argon before substrates were added", "2 min catalyst preactivation", e4)


def make_specific(kind, idx, substrate_mid, product_mid, yield_value, table_page, char_pages,
                  product_mass=None, product_mass_page=None, e_assignment=False):
    table_ev = EV[table_page]
    rid = f"{DOC}:{kind}_{idx}"
    evidence = [table_ev, EV[4]] + [EV[p] for p in char_pages]
    if e_assignment:
        evidence.append(EV[2])
    evidence = list(dict.fromkeys(evidence))
    r = st.new_reaction(rid, DOC, evidence)
    set_reaction_header(r, "specific_experiment", ["substrate_scope"],
                        {"table": 1 if kind == "allene" else 2, "entry": idx, "pages": sorted(set([table_page, 4] + char_pages))},
                        "Cobalt(I)-Catalyzed [6π+2π] Cycloaddition of 1,2-Dienes to N-Carbocholesteroxyazepine" if kind == "allene" else
                        "Cobalt(I)-Catalyzed [6π+2π] Cycloaddition of 1,3-Diynes to N-Carbocholesteroxyazepine",
                        "[6π+2π] cycloaddition", evidence)
    r["general_scheme_ids"] = [f"{DOC}:current_{kind}_generic"]
    add_common_specific_participants(r, substrate_mid, product_mid, table_ev)
    add_common_specific_conditions(r, table_ev)
    y = st.new_outcome("o_yield", "yield", [product_mid], [table_ev])
    y["product_role"] = f("target product", "target product", [table_ev])
    y["value"] = f(f"{yield_value}%", {"number": float(yield_value), "unit": "%"}, [table_ev])
    y["yield_type"] = f("isolated", "isolated", [table_ev])
    y["measurement_basis"] = f("Yields of products isolated by column chromatography.", "isolated after column chromatography", [table_ev])
    y["scope"] = f(f"Table {1 if kind == 'allene' else 2}, entry {idx}", {"table": 1 if kind == "allene" else 2, "entry": idx}, [table_ev])
    r["outcomes"].append(y)
    if product_mass is not None:
        amount_ev = [EV[product_mass_page]]
        o = st.new_outcome("o_product_mass", "product_amount", [product_mid], amount_ev)
        o["product_role"] = f("target product", "target product", amount_ev)
        o["value"] = f(f"{product_mass:.3f} g", {"number": product_mass, "unit": "g"}, amount_ev)
        o["measurement_basis"] = f("reported with isolated product characterization", "isolated product mass", amount_ev)
        o["scope"] = f("product characterization paragraph", "single isolated product", amount_ev)
        r["outcomes"].append(o)
    if e_assignment:
        o = st.new_outcome("o_exocyclic_geometry", "selectivity", [product_mid], [EV[2]])
        o["product_role"] = f("target product", "target product", [EV[2]])
        o["value"] = f("(E)-configuration of the exocyclic C(7)-C(10) double bond", "E", [EV[2]])
        o["measurement_basis"] = f("NOESY cross-peaks described for 3a,c-e; phenyl/H2C(8) NOE described for 3b", "NOESY assignment", [EV[2]])
        o["scope"] = f("exocyclic C(7)-C(10) double bond", {"metric": "E/Z assignment", "assigned": "E", "ratio": None}, [EV[2]])
        r["outcomes"].append(o)
    partial_cycloaddition_classification(r, [table_ev])
    r["notes"] = [
        f("Products were observed as two N-(CO)O-cholesteryl rotamers in a 1:1 ratio.",
          "two N-(CO)O-cholesteryl rotamers, 1:1", [EV[2]] if kind == "allene" else [EV[3]]),
        f("Supporting Information is mentioned by the article but no SI file is explicitly supplied or associated in this workspace.",
          "SI not associated", evidence),
    ]
    unresolved = [
        "Complete stereochemical SMILES for compound 1 and the cholesterol-containing product were not safely established from pages 1-6.",
        "Exact atom mapping, all bond-order changes, and implicit-hydrogen balance remain unresolved; classification is partial.",
        "Amounts of petroleum ether, silica gel, and ethyl acetate used in workup/purification are unreported.",
        "No explicitly associated SI was provided; details beyond the main-text general procedure remain unavailable.",
    ]
    if product_mass is None:
        unresolved.append("Product mass is not reported within the bounded pages used for this record.")
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [], "unresolved": unresolved}
    r["completeness"] = "partial"
    return r


def make_generic(kind):
    page = 3 if kind == "allene" else 5
    evidence = [EV[1], EV[3] if kind == "allene" else EV[5]]
    r = st.new_reaction(f"{DOC}:current_{kind}_generic", DOC, evidence)
    name = "[6π+2π] cycloaddition of terminal 1,2-dienes to N-carbocholesteroxyazepine" if kind == "allene" else \
           "[6π+2π] cycloaddition of symmetric 1,3-diynes to N-carbocholesteroxyazepine"
    set_reaction_header(r, "generic_scheme", ["other"],
                        {"regions": ["page 1 graphical abstract", f"page {page} scheme/table header"]},
                        name, "[6π+2π] cycloaddition", evidence)
    substrate_mid = "m_allene_template" if kind == "allene" else "m_diyne_template"
    product_mid = "m_product3_template" if kind == "allene" else "m_product5_template"
    r["participants"] = [
        participant("p_azepine", "m_azepine_1", "reactant", evidence),
        participant("p_unsaturated", substrate_mid, "reactant", evidence),
        participant("p_product", product_mid, "product", evidence, product_role="generic target product family"),
        participant("p_co", "m_co_acac2_dppe", "catalyst", evidence),
        participant("p_zn", "m_zinc", "reagent", evidence),
        participant("p_zni2", "m_zni2", "additive", evidence),
        participant("p_dce", "m_dce", "solvent", evidence),
    ]
    r["solvent_conditions"]["composition"] = f("DCE", ["m_dce"], evidence)
    r["steps"] = [step("s1_generic_hold", evidence, "reaction",
                       [p["participant_id"] for p in r["participants"]], "Generic scheme reaction hold",
                       temperature=("60 °C", {"number": 60.0, "unit": "°C"}),
                       time=("20 h", {"number": 20.0, "unit": "h"}), atmosphere="dry argon",
                       method="original_extraction")]
    y0, y1 = (79, 92) if kind == "allene" else (77, 90)
    o = st.new_outcome("o_scope_yield_range", "yield", [product_mid], evidence)
    o["product_role"] = f("generic target product family", "generic target product family", evidence)
    o["value"] = f(f"{y0}-{y1}%", {"minimum": y0, "maximum": y1, "unit": "%", "raw_text": f"{y0}-{y1}%"}, evidence)
    o["yield_type"] = f("isolated", "isolated", [EV[page]])
    o["measurement_basis"] = f("Yields of products isolated by column chromatography.", "isolated after column chromatography", [EV[page]])
    o["scope"] = f("five examples", {"examples": 5, "aggregate_range": True}, evidence)
    r["outcomes"] = [o]
    partial_cycloaddition_classification(r, evidence)
    r["notes"] = [f("Generic current-work scheme; concrete table entries are stored separately and are not duplicated by this record.",
                    "generic scheme linked to five table entries", evidence)]
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
                   "unresolved": ["Generic R-group structures are not complete concrete molecules.",
                                  "Classification is scheme-level and does not claim complete atom mapping."]}
    r["completeness"] = "partial"
    return r


def make_literature(suffix, location, reactant_mids, product_mid, source_name, condition=None, page=2):
    evidence = [EV[page]]
    r = st.new_reaction(f"{DOC}:literature_{suffix}", DOC, evidence)
    set_reaction_header(r, "literature_example", [], location, source_name,
                        "cited [6π+2π] cycloaddition", evidence)
    for i, mid in enumerate(reactant_mids, 1):
        r["participants"].append(participant(f"p_r{i}", mid, "reactant", evidence))
    r["participants"].append(participant("p_product", product_mid, "product", evidence,
                                        product_role="cited literature product"))
    if condition:
        r["steps"] = [step("s1_scheme_arrow", evidence, "reaction",
                           [p["participant_id"] for p in r["participants"]], condition,
                           method="original_extraction")]
    r["classification"]["coverage"] = inferred("unresolved", evidence,
        "The cited schematic structures were inventoried, but complete atom mapping/classification was not reconstructed from this article.")
    r["classification"]["limitations"] = ["Cited prior work, not a current-paper experiment.",
                                              "Conditions and yields are incomplete in this review scheme."]
    r["notes"] = [f("This arrow is explicitly presented as previously published work.",
                    "cited prior-work arrow", evidence)]
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
                   "unresolved": ["Full identities, conditions, yields, stereochemistry, and atom mapping require the cited paper."]}
    r["completeness"] = "partial"
    return r


# Molecule registry.
molecules = []
molecules.append(unresolved_molecule(
    "m_azepine_1", "N-carbocholesteroxyazepine", "1", [EV[3], EV[4], IMG[3]],
    chirality="partially_specified", geometry="unspecified", representation="not_encoded",
    atom={"cholesteryl_fragment": "multiple wedge/dash stereocenters shown; exact isomeric SMILES not encoded"},
    assignment={"azepine": "N-carbamate of cholesteryl group", "cholesteryl_stereo": "source drawing retained"},
    note="The cholesteryl stereocenters are visibly specified in the source drawing, but a complete checked isomeric SMILES was not established in this bounded pass."))

allene_data = [
    ("2a", "1,2-diene 2a (R = Hex)", "C=C=CCCCCCC"),
    ("2b", "1,2-diene 2b (R = Ph)", "C=C=Cc1ccccc1"),
    ("2c", "1,2-diene 2c (R = Bn)", "C=C=CCc1ccccc1"),
    ("2d", "1,2-diene 2d (R = (CH2)4OH)", "C=C=CCCCCO"),
    ("2e", "1,2-diene 2e (R = (CH2)4Br)", "C=C=CCCCCBr"),
]
for label, name, smi in allene_data:
    molecules.append(smiles_molecule(f"m_{label}", name, label, smi, [EV[3], IMG[3]],
        chirality="achiral", geometry="not_applicable", representation="partial",
        axial={"allene_axis": "not stereogenic because one terminus is CH2"}, allene=True,
        note="Connectivity is a machine candidate expanded from the table R group; RDKit reports an unspecified cumulene double-bond feature, so representation_status remains partial."))

diyne_data = [
    ("4a", "1,3-diyne 4a (R = Hex)", "CCCCCCC#CC#CCCCCCC"),
    ("4b", "1,3-diyne 4b (R = SiMe3)", "C[Si](C)(C)C#CC#C[Si](C)(C)C"),
    ("4c", "1,3-diyne 4c (R = (CH2)2OH)", "OCCC#CC#CCCO"),
    ("4d", "1,3-diyne 4d (R = (CH2)3SBut)", "CC(C)(C)SCCC#CC#CCCSC(C)(C)C"),
    ("4e", "1,3-diyne 4e (R = CH(CH2)2)", "C1CC1C#CC#CC1CC1"),
]
for label, name, smi in diyne_data:
    molecules.append(smiles_molecule(f"m_{label}", name, label, smi, [EV[5], IMG[5]],
        note="Connectivity is a machine candidate expanded from the Table 2 R-group assignment; it remains needs_review."))

product_names = {
    "3a": "Cholesteryl (1S*,6R*)-7-[(E)-heptylidene]-9-azabicyclo[4.2.1]nona-2,4-diene-9-carboxylate equivalent with the (1R*,6S*) form",
    "3b": "Cholesteryl (1S*,6R*)-7-[(E)-phenylmethylidene]-9-azabicyclo[4.2.1]nona-2,4-diene-9-carboxylate equivalent with the (1R*,6S*) form",
    "3c": "Cholesteryl (1S*,6R*)-7-[(E)-2-phenylethylidene]-9-azabicyclo[4.2.1]nona-2,4-diene-9-carboxylate equivalent with the (1R*,6S*) form",
    "3d": "Cholesteryl (1S*,6R*)-7-[(E)-5-hydroxypentylidene]-9-azabicyclo[4.2.1]nona-2,4-diene-9-carboxylate equivalent with the (1R*,6S*) form",
    "3e": "Cholesteryl (1S*,6R*)-7-[(E)-5-bromopentylidene]-9-azabicyclo[4.2.1]nona-2,4-diene-9-carboxylate equivalent with the (1R*,6S*) form",
    "5a": "Cholesteryl (1S*,6R*)-7-hexyl-8-(oct-1-yn-1-yl)-9-azabicyclo[4.2.1]nona-2,4,7-triene-9-carboxylate equivalent with the (1R*,6S*) form",
    "5b": "Cholesteryl (1S*,6R*)-7-(trimethylsilyl)-8-((trimethylsilyl)ethynyl)-9-azabicyclo[4.2.1]nona-2,4,7-triene-9-carboxylate equivalent with the (1R*,6S*) form",
    "5c": "Cholesteryl (1S*,6R*)-7-(4-hydroxybut-1-yn-1-yl)-8-(2-hydroxyethyl)-9-azabicyclo[4.2.1]nona-2,4,7-triene-9-carboxylate equivalent with the (1R*,6S*) form",
    "5d": "Cholesteryl (1S*,6R*)-7-(5-(tert-butylthio)pent-1-yn-1-yl)-8-(3-(tert-butylthio)propyl)-9-azabicyclo[4.2.1]nona-2,4,7-triene-9-carboxylate equivalent with the (1R*,6S*) form",
    "5e": "Table 2 product 5e (R = CH(CH2)2, cyclopropyl)",
}
product_pages = {"3a": [4, 5], "3b": [4, 5], "3c": [5], "3d": [5, 6], "3e": [6],
                 "5a": [6], "5b": [6], "5c": [6], "5d": [6], "5e": [5]}
for label, name in product_names.items():
    is_allene_product = label.startswith("3")
    evidence = list(dict.fromkeys(([EV[2]] if is_allene_product else [EV[3]]) + [EV[p] for p in product_pages[label]]))
    molecules.append(unresolved_molecule(
        f"m_{label}", name, label, evidence, chirality="partially_specified",
        geometry="specified" if is_allene_product else "unspecified", representation="not_encoded",
        atom={"bridgeheads": "relative descriptors (1S*,6R*) equivalent with (1R*,6S*) are reported",
              "cholesteryl_fragment": "stereocenters shown/retained in the source"},
        double={"exocyclic_C7_C10": "E"} if is_allene_product else None,
        mixture={"type": "N-(CO)O-cholesteryl rotamers", "ratio": "1:1"},
        assignment={"method": "NOESY for exocyclic E assignment"} if is_allene_product else
                   {"source_description": "two N-(CO)O-cholesteryl conformational stereoisomers/rotamers"},
        note="A complete isomeric SMILES was not encoded: the source reports/draws cholesterol stereocenters, relative bicyclic descriptors, and a 1:1 rotamer pair."))

simple_reagents = [
    ("m_zinc", "zinc", "Zn", "[Zn]", [EV[3], EV[4], EV[5]], "not_applicable"),
    ("m_zni2", "zinc iodide", "ZnI2", "[I-].[I-].[Zn+2]", [EV[3], EV[4], EV[5]], "not_applicable"),
    ("m_dce", "1,2-dichloroethane", "DCE", "ClCCCl", [EV[3], EV[4], EV[5]], "encoded"),
    ("m_ethyl_acetate", "ethyl acetate", "ethyl acetate", "CCOC(C)=O", [EV[4]], "encoded"),
]
for mid, name, label, smi, evidence, rep in simple_reagents:
    molecules.append(smiles_molecule(mid, name, label, smi, evidence,
        chirality="not_applicable" if mid in ("m_zinc", "m_zni2") else "achiral",
        geometry="not_applicable", representation=rep,
        note="Candidate standard representation for an explicitly named reagent/solvent."))
molecules += [
    unresolved_molecule("m_co_acac2_dppe", "Co(acac)2(dppe)", "Co(acac)2(dppe)", [EV[3], EV[4], EV[5]],
                        chirality="unresolved", geometry="unresolved", representation="not_encoded",
                        note="Coordination structure/oxidation-state representation was not reconstructed from the main-text abbreviation."),
    unresolved_molecule("m_petroleum_ether", "petroleum ether", "petroleum ether", [EV[4]],
                        chirality="unresolved", geometry="not_applicable", representation="not_encoded",
                        mixture={"composition": "unreported hydrocarbon mixture"},
                        note="Named solvent mixture; composition and boiling range are unreported."),
    unresolved_molecule("m_silica_gel", "silica gel", "silica gel", [EV[4]], chirality="not_applicable",
                        geometry="not_applicable", representation="not_encoded", note="Stationary/filter medium; no unique molecular SMILES."),
]

for mid, name, label, evidence, note in [
    ("m_allene_template", "terminal 1,2-diene family R-CH=C=CH2", "2a-e", [EV[1], EV[3]], "Generic R-group family; concrete instances 2a-e are separate registry entries."),
    ("m_diyne_template", "symmetric 1,3-diyne family R-C≡C-C≡C-R", "4a-e", [EV[1], EV[5]], "Generic R-group family; concrete instances 4a-e are separate registry entries."),
    ("m_product3_template", "(E)-9-azabicyclo[4.2.1]nona-2,4-diene product family", "3a-e", [EV[1], EV[3]], "Generic current-work product family."),
    ("m_product5_template", "9-azabicyclo[4.2.1]nona-2,4,7-triene product family", "5a-e", [EV[1], EV[5]], "Generic current-work product family."),
    ("m_lit_azepine_me", "N-carbomethoxyazepine", "R = Me azepine", [EV[2]], "Cited Scheme 1 structure; no current-paper experimental identity claim."),
    ("m_lit_azepine_et", "N-carbethoxyazepine", "R = Et azepine", [EV[2]], "Cited Scheme 1 structure; no current-paper experimental identity claim."),
    ("m_lit_product_me", "cited methyl-carbamate 9-azabicyclo[4.2.1]nonane cycloadduct", "Scheme 1 upper-left product", [EV[2]], "Cited product drawing."),
    ("m_lit_product_et", "cited ethyl-carbamate 9-azabicyclo[4.2.1]nonane cycloadduct", "Scheme 1 lower-left product", [EV[2]], "Cited product drawing."),
    ("m_lit_cr_ncyano_azepine", "tricarbonyl(η6-N-cyanoazepine)chromium(0)", "Cr(CO)3 N-cyanoazepine", [EV[2]], "Cited noncatalytic Scheme 1 complex."),
    ("m_lit_diyne_co2et", "symmetric diyne, R = CO2Et", "R = CO2Et diyne", [EV[2]], "Cited Scheme 1 generic diyne."),
    ("m_lit_diyne_cn", "symmetric diyne, R = CN", "R = CN diyne", [EV[2]], "Cited Scheme 1 generic diyne."),
    ("m_lit_cr_product_co2et", "Cr(CO)3-coordinated diyne cycloadduct, R = CO2Et", "Scheme 1 R=CO2Et product", [EV[2]], "Cited product drawing."),
    ("m_lit_cr_product_cn", "Cr(CO)3-coordinated diyne cycloadduct, R = CN", "Scheme 1 R=CN product", [EV[2]], "Cited product drawing."),
    ("m_lit_cr_acrylate_product", "Cr(CO)3-coordinated ethyl acrylate cycloadduct", "Scheme 1 lower-right product", [EV[2]], "Cited product drawing."),
    ("m_lit_alkyne_template", "alkyne family HC≡C-R", "previous-work alkyne", [EV[3]], "Our previous work panel; generic R group."),
    ("m_lit_previous_product", "previously reported 7-R-9-azabicyclo[4.2.1]nona-2,4,7-triene product family", "previous-work product", [EV[3]], "Cited Molecules 2021 product family."),
]:
    molecules.append(generic_source_molecule(mid, name, label, evidence, note))
molecules.append(smiles_molecule("m_ethyl_acrylate", "ethyl acrylate", "ethyl acrylate", "C=CC(=O)OCC", [EV[2]],
                                 note="Candidate connectivity for the explicitly labeled cited reagent."))


# Current generic schemes and ten independent scope experiments.
reactions = [make_generic("allene"), make_generic("diyne")]
allene_rows = [
    (1, "m_2a", "m_3a", 92, [4, 5], 0.580, 4),
    (2, "m_2b", "m_3b", 86, [4, 5], 0.535, 5),
    (3, "m_2c", "m_3c", 89, [5], 0.566, 5),
    (4, "m_2d", "m_3d", 79, [5, 6], 0.488, 5),
    (5, "m_2e", "m_3e", 84, [6], 0.572, 6),
]
for idx, sub, prod, y, pages, mass, mass_page in allene_rows:
    reactions.append(make_specific("allene", idx, sub, prod, y, 3, pages, mass, mass_page, True))
diyne_rows = [
    (1, "m_4a", "m_5a", 90, [6], 0.652, 6),
    (2, "m_4b", "m_5b", 88, [6], 0.616, 6),
    (3, "m_4c", "m_5c", 78, [6], 0.502, 6),
    (4, "m_4d", "m_5d", 83, [6], 0.678, 6),
    (5, "m_4e", "m_5e", 77, [], None, None),
]
for idx, sub, prod, y, pages, mass, mass_page in diyne_rows:
    reactions.append(make_specific("diyne", idx, sub, prod, y, 5, pages, mass, mass_page, False))

# Five distinct cited arrows in Scheme 1 plus one cited arrow in Scheme 2.
reactions += [
    make_literature("scheme1_cr_catalytic_me", {"page": 2, "scheme": 1, "arrow": "upper-left catalytic"},
                    ["m_lit_azepine_me", "m_ethyl_acrylate"], "m_lit_product_me",
                    "Cr(0)-catalyzed [6π+2π] cycloaddition to ethyl acrylate (Tetrahedron 1993)", "NpCr(CO)3; Np = naphthalene"),
    make_literature("scheme1_cr_catalytic_et", {"page": 2, "scheme": 1, "arrow": "lower-left catalytic"},
                    ["m_lit_azepine_et", "m_ethyl_acrylate"], "m_lit_product_et",
                    "Cr-catalyzed [6π+2π] cycloaddition to ethyl acrylate (Org. Lett. 2000)", "resin-based [Cr] catalyst"),
    make_literature("scheme1_photochemical_diyne_co2et", {"page": 2, "scheme": 1, "arrow": "upper-right noncatalytic"},
                    ["m_lit_cr_ncyano_azepine", "m_lit_diyne_co2et"], "m_lit_cr_product_co2et",
                    "Photochemical noncatalytic [6π+2π] cycloaddition, R = CO2Et", "hν"),
    make_literature("scheme1_photochemical_diyne_cn", {"page": 2, "scheme": 1, "arrow": "middle-right noncatalytic"},
                    ["m_lit_cr_ncyano_azepine", "m_lit_diyne_cn"], "m_lit_cr_product_cn",
                    "Photochemical noncatalytic [6π+2π] cycloaddition, R = CN", "hν"),
    make_literature("scheme1_photochemical_acrylate", {"page": 2, "scheme": 1, "arrow": "lower-right noncatalytic"},
                    ["m_lit_cr_ncyano_azepine", "m_ethyl_acrylate"], "m_lit_cr_acrylate_product",
                    "Photochemical noncatalytic [6π+2π] cycloaddition to ethyl acrylate", "hν"),
    make_literature("scheme2_previous_alkyne", {"page": 3, "scheme": 2, "panel": "Our previous work"},
                    ["m_azepine_1", "m_lit_alkyne_template"], "m_lit_previous_product",
                    "Previously reported cobalt-catalyzed cycloaddition of an alkyne to N-carbocholesteroxyazepine", None, page=3),
]


source_context = [
    {"context_id": f"{DOC}:bibliographic_metadata", "kind": "document_metadata",
     "title": "Cobalt(I)-Catalyzed [6π + 2π] Cycloaddition of 1,2-Dienes and 1,3-Diynes to N-Carbocholesteroxyazepine in the Synthesis of Previously Undescribed Heterofunctional 9-Azabicyclo[4.2.1]nonadi(tri)enes",
     "doi": "10.1021/acsomega.1c03321", "evidence_ids": [EV[1]]},
    {"context_id": f"{DOC}:noesy_figures_2_3", "kind": "characterization",
     "finding": "Figures 2 and 3 are NOESY spectra supporting the E assignment of 3a and 3b; they are not separate syntheses.",
     "evidence_ids": [EV[4], IMG[4]], "related_reaction_ids": [f"{DOC}:allene_1", f"{DOC}:allene_2"]},
    {"context_id": f"{DOC}:nmr_hrms_characterization_pages_4_6", "kind": "characterization",
     "finding": "NMR, HRMS, melting point, optical rotation, and Rf data on pages 4-6 characterize already-recorded products; no standalone synthesis records were created from peak lists.",
     "evidence_ids": [EV[4], EV[5], EV[6]],
     "related_reaction_ids": [r["reaction_id"] for r in reactions if r["record_type"] == "specific_experiment"]},
    {"context_id": f"{DOC}:si_limitation", "kind": "source_limitation",
     "finding": "The article advertises Supporting Information, but source_mapping.json leaves SI association unestablished and no SI document was supplied; no SI facts were imported.",
     "evidence_ids": [EV[1]]},
]

bundle = st.empty(DATASET["documents"])
bundle["molecules"] = molecules
bundle["reactions"] = reactions
bundle["source_context"] = source_context
st.write(RUN / "submitted.json", bundle)

raw_inventory = {
    "document_id": DOC,
    "bounded_pages": [1, 2, 3, 4, 5, 6],
    "semantic_accuracy_verified": False,
    "inventory": {
        "current_generic_schemes": 2,
        "current_specific_experiments": 10,
        "cited_literature_arrows": 6,
        "standalone_characterization_records": 0,
    },
    "table_rows": {
        "Table 1": [{"entry": i, "substrate": f"2{'abcde'[i-1]}", "product": f"3{'abcde'[i-1]}", "yield_percent": y}
                    for i, y in enumerate([92, 86, 89, 79, 84], 1)],
        "Table 2": [{"entry": i, "substrate": f"4{'abcde'[i-1]}", "product": f"5{'abcde'[i-1]}", "yield_percent": y}
                    for i, y in enumerate([90, 88, 78, 83, 77], 1)],
    },
    "limitations": [
        "Only PDF pages 1-6 were claimed as visually covered.",
        "Page 7 characterization of 5e was not used as covered evidence.",
        "No SI document was explicitly associated.",
        "Complete stereochemical SMILES for the cholesterol-containing substrate/products were not safely established.",
    ],
}
st.write(RUN / "raw_model/pages_0001_0006_extraction.json", raw_inventory)

specific_allene = [f"{DOC}:allene_{i}" for i in range(1, 6)]
specific_diyne = [f"{DOC}:diyne_{i}" for i in range(1, 6)]
lit_scheme1 = [r["reaction_id"] for r in reactions if r["reaction_id"].startswith(f"{DOC}:literature_scheme1")]
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title, DOI, and graphical abstract", "finding": "Article metadata plus two current-work generic [6π+2π] pathways; concrete examples are tabulated later.", "evidence_ids": [EV[1], IMG[1]]},
            {"region": "abstract and introduction", "finding": "Aggregate 79-92%/77-90% scope statement and contextual cited literature; no extra concrete experiment.", "evidence_ids": [EV[1], IMG[1]]}],
         "record_ids": [f"{DOC}:current_allene_generic", f"{DOC}:current_diyne_generic"],
         "unresolved": ["Graphical abstract omits concrete R identities and full experimental procedure."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Figure 1", "finding": "Biologically active azacycle examples; contextual structures, not current synthesis experiments.", "evidence_ids": [EV[2], IMG[2]]},
            {"region": "Scheme 1, five arrows", "finding": "Two catalytic and three noncatalytic/photochemical cited prior-work cycloadditions; stored as literature_example records.", "evidence_ids": [EV[2], IMG[2]]},
            {"region": "Results and Discussion prose", "finding": "Current 3a-e E assignments and 1:1 rotamer observation; Table 2 family introduction.", "evidence_ids": [EV[2], IMG[2]]}],
         "record_ids": lit_scheme1 + specific_allene,
         "unresolved": ["Scheme 1 does not reproduce complete procedures, yields, or fully machine-readable structures from the cited papers."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Scheme 2 upper previous-work panel", "finding": "One cited alkyne cycloaddition arrow; separated from this work.", "evidence_ids": [EV[3], IMG[3]]},
            {"region": "Scheme 2 current-work panels", "finding": "Two generic current-work pathways, repeated by the graphical abstract and linked to concrete table rows.", "evidence_ids": [EV[3], IMG[3]]},
            {"region": "Table 1 rows 1-5 and footnotes", "finding": "Five independent 1,2-diene substrate-scope experiments with explicit charges, conditions, isolated yields, and chromatography basis.", "evidence_ids": [EV[3], IMG[3]]},
            {"region": "Experimental general information", "finding": "Analytical instrumentation and dry-argon/general solvent statements; no separate synthesis record.", "evidence_ids": [EV[3], IMG[3]]}],
         "record_ids": [f"{DOC}:literature_scheme2_previous_alkyne", f"{DOC}:current_allene_generic"] + specific_allene,
         "unresolved": ["The previous-work panel lacks a complete procedure and concrete R identity."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "Figures 2 and 3", "finding": "NOESY characterization supporting E geometry for 3a and 3b; associated with existing reactions, not new syntheses.", "evidence_ids": [EV[4], IMG[4]]},
            {"region": "General Reaction Procedure", "finding": "Shared catalyst preactivation, sequential addition, 60 °C/20 h hold, deactivation, filtration, evaporation, and chromatography applied to all ten current experiments.", "evidence_ids": [EV[4], IMG[4]]},
            {"region": "3a characterization beginning", "finding": "Isolated mass/yield and characterization for 3a; peak lists retained as characterization context.", "evidence_ids": [EV[4], IMG[4]]}],
         "record_ids": specific_allene + specific_diyne,
         "unresolved": ["Volumes/amounts for petroleum ether, silica gel, and ethyl acetate are not reported."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "Table 2 rows 1-5 and footnotes", "finding": "Five independent 1,3-diyne substrate-scope experiments with explicit charges, conditions, isolated yields, and chromatography basis.", "evidence_ids": [EV[5], IMG[5]]},
            {"region": "3a continuation and 3b-3d characterization", "finding": "NMR/HRMS and isolated mass data linked to existing 3a-3d reactions; not separate syntheses.", "evidence_ids": [EV[5], IMG[5]]}],
         "record_ids": [f"{DOC}:current_diyne_generic"] + specific_diyne + [f"{DOC}:allene_{i}" for i in range(1, 5)],
         "unresolved": ["5e characterization continues outside the bounded page range and was not imported."]},
        {"page": 6, "state": "attempted", "regions": [
            {"region": "3d continuation, 3e, and 5a-5d characterization", "finding": "Isolated masses and analytical characterization associated with existing table experiments; no additional reaction arrows or table rows.", "evidence_ids": [EV[6], IMG[6]]}],
         "record_ids": [f"{DOC}:allene_4", f"{DOC}:allene_5"] + [f"{DOC}:diyne_{i}" for i in range(1, 5)],
         "unresolved": ["The 5d characterization continues on page 7, which is outside this pass; only page-6 facts were used."]},
    ],
}
st.write(RUN / "coverage/pages_0001_0006.json", coverage)

print(json.dumps({"molecules": len(molecules), "reactions": len(reactions),
                  "specific_experiments": 10, "generic_schemes": 2,
                  "literature_examples": 6}, ensure_ascii=False, indent=2))
