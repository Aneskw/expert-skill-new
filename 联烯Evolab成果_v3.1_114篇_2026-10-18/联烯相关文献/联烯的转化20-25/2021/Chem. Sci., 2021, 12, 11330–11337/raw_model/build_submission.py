#!/usr/bin/env python3
"""Build the page-1--6 v3.1 machine submission from the visual transcription."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DOC = "doc_b620dbe6e0542262"
PREFIX = "d1sc02167b"
raw = json.loads((RUN / "raw_model/visual_transcription.json").read_text())
dataset = json.loads((RUN / "dataset.json").read_text())
dataset["documents"][0]["title"] = raw["title"]


def evid(page: int):
    return [f"{DOC}:page{page}:image", f"{DOC}:page{page}:text"]


def direct(raw_value, normalized=None, page=1):
    if normalized is None:
        normalized = raw_value
    return st.fact(raw_value, normalized, evid(page))


def inferred(normalized, page, derivation):
    return st.fact(None, normalized, evid(page), method="inference", derivation=derivation)


def derived(normalized, page, derivation, raw_value=None):
    return st.fact(raw_value, normalized, evid(page), method="structure_derivation", derivation=derivation)


def mid(label: str):
    return f"{PREFIX}_mol_{label.lower().replace('-', '_')}"


def rid(suffix: str):
    return f"{DOC}:{suffix}"


SMILES = {
    "1a": "c1ccc2c(c1)ccn2-c1ncccn1",
    "2a": "COC(=O)OCC=C=Cc1ccccc1",
    "2b": "COC(=O)OCC=C=Cc1ccc(C)cc1",
    "2c": "COC(=O)OCC=C=Cc1ccc(OC)cc1",
    "2d": "COC(=O)OCC=C=Cc1ccc(F)cc1",
    "2e": "COC(=O)OCC=C=Cc1ccc(Cl)cc1",
    "2f": "COC(=O)OCC=C=Cc1ccc(Br)cc1",
    "2g": "COC(=O)OCC=C=Cc1ccc(I)cc1",
    "2h": "COC(=O)OCC=C=Cc1ccc(C#N)cc1",
    "2i": "COC(=O)OCC=C=Cc1ccc(C(F)(F)F)cc1",
    "2j": "COC(=O)OCC=C=Cc1cccc2ccccc12",
    "2n": "COC(=O)OCC=C=CC1CCCCC1",
    "PivOH": "CC(C)(C)C(=O)O",
    "toluene": "Cc1ccccc1",
    "dioxane": "O1CCOCC1",
    "acetonitrile": "CC#N"
}

NAMES = {
    "1a": "1-(pyrimidin-2-yl)indole",
    "2a": "4-phenyl-2,3-allenyl methyl carbonate",
    "2f": "4-(4-bromophenyl)-2,3-allenyl methyl carbonate",
    "2h": "4-(4-cyanophenyl)-2,3-allenyl methyl carbonate",
    "2n": "4-cyclohexyl-2,3-allenyl methyl carbonate",
    "RhCat": "[Cp*RhCl2]2",
    "AgSbF6": "silver hexafluoroantimonate",
    "AgBF4": "silver tetrafluoroborate",
    "PivOH": "pivalic acid",
    "toluene": "toluene",
    "dioxane": "1,4-dioxane",
    "acetonitrile": "acetonitrile"
}


def source_page(label: str):
    if label in {"1a", "2a", "3aa", "3am", "3ah_ZE98_2", "3ah_ZE95_5", "RhCat", "AgSbF6", "AgBF4", "PivOH", "toluene", "dioxane", "acetonitrile"}:
        return 2
    if label in {"1h", "1i", "1j", "1k", "3ka"}:
        return 4
    if label.startswith("2n") or label in {"2o", "2p", "2q", "2r", "2s", "2t"} or label in {"3an", "3ao", "3ap", "3aq", "3ar", "3as", "3at"}:
        return 6
    return 3


product_labels = {
    row["product"] for row in raw["table2"]
} | {row.get("product") for row in raw["table3"] if row.get("product")} | {
    row["product"] for row in raw["table4"]
} | {"3am", "3ah_ZE98_2", "3ah_ZE95_5"}

labels = {
    "1a", "1b", "1c", "1d", "1e", "1f", "1g", "1h", "1i", "1j", "1k",
    "2a", "2b", "2c", "2d", "2e", "2f", "2g", "2h", "2i", "2j", "2k", "2l", "2m",
    "2n", "2o", "2p", "2q", "2r", "2s", "2t",
    "RhCat", "AgSbF6", "AgBF4", "PivOH", "toluene", "dioxane", "acetonitrile"
} | product_labels


def make_molecule(label: str):
    page = source_page(label)
    m = st.new_molecule(mid(label), DOC, evid(page))
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = direct(label, label, page)
    name = NAMES.get(label)
    if label in product_labels and not name:
        name = f"reported product {label} (geometric-isomer sample)"
    elif not name:
        name = f"source compound {label}"
    m["names"]["original"] = direct(name, name, page)
    m["names"]["english"] = direct(name, name, page)
    m["identity_note"] = direct(
        f"Identity retained from the source label {label}; a concrete structure is only encoded when confidently transcribed from the page drawing.",
        f"source label {label}; machine identity needs review", page)
    smi = SMILES.get(label)
    if smi:
        m["representation_type"] = "smiles"
        m["smiles"] = derived(smi, page, "Candidate connectivity transcribed from the source drawing; RDKit parsing does not verify source identity.", f"structure drawing for {label}")
        m["stereochemistry"]["isomeric_smiles"] = derived(smi, page, "No additional stereotags are asserted beyond the candidate registry SMILES.")
        if label == "1a":
            m["stereochemistry"]["chirality_status"] = derived("achiral", page, "Candidate graph has no tetrahedral or axial stereogenic element.")
            m["stereochemistry"]["geometry_status"] = derived("not_applicable", page, "Candidate graph has no configurable nonaromatic alkene geometry.")
            m["stereochemistry"]["representation_status"] = derived("not_applicable", page, "No stereogenic element is represented in the candidate graph.")
        elif label.startswith("2"):
            m["stereochemistry"]["chirality_status"] = inferred("unspecified", page, "Allenyl substrates can possess axial stereochemistry; the experimental entries do not assign a configuration.")
            m["stereochemistry"]["axial"] = derived(
                [{"axis": "C=C=C", "configuration": "unspecified"}], page,
                "Potential allene axial stereochemistry is retained without inventing R_a/S_a; Fig. 4 uses (R_a)-2a/(R_a)-2n only as calculation examples." if label in {"2a", "2n"} else "Potential allene axial stereochemistry is not assigned in the source entry.")
            m["stereochemistry"]["geometry_status"] = inferred("not_applicable", page, "Allene axial stereochemistry is recorded under chirality rather than alkene E/Z.")
            m["stereochemistry"]["representation_status"] = inferred("not_encoded", page, "The candidate SMILES does not encode allene axial configuration.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = st.fact(evidence=evid(page))
        m["stereochemistry"]["representation_status"] = inferred("unresolved", page, "No concrete stereochemical SMILES was asserted from the small page drawing.")
        m["stereochemistry"]["chirality_status"] = inferred("unresolved", page, "Center/axial chirality was inspected but cannot be fully assigned from the provided page alone.")
        m["stereochemistry"]["geometry_status"] = inferred("unresolved", page, "E/Z or other geometry is unresolved for this isolated registry identity.")
    if label in product_labels:
        m["stereochemistry"]["geometry_status"] = direct(
            "Z/E or E/Z mixture; composition is experiment-specific",
            "mixture", page)
        m["stereochemistry"]["mixture"] = direct(
            "geometric-isomer mixture reported; see per-experiment selectivity outcome",
            {"type": "geometric_isomers", "composition": "experiment_specific"}, page)
        m["stereochemistry"]["representation_status"] = inferred(
            "not_encoded", page, "A single SMILES cannot encode the experiment-specific geometric-isomer mixture retained under one source product label.")
    if label in {"2k", "3ak", "2t", "3at"}:
        m["stereochemistry"]["chirality_status"] = direct(
            "wedge/dash stereochemistry shown in the source drawing; absolute descriptor not transcribed",
            "partially_specified", page)
        m["stereochemistry"]["assignment"] = direct(
            "source drawing contains specified relative/absolute-looking stereobonds; descriptors and full encoding remain unresolved",
            {"source": "wedge/dash drawing", "descriptor": "not assigned"}, page)
        m["stereochemistry"]["representation_status"] = inferred("not_encoded", page, "Source stereobonds are retained descriptively because no verified isomeric SMILES was created.")
    return m


bundle = st.empty(dataset["documents"])
bundle["molecules"] = [make_molecule(label) for label in sorted(labels)]


def set_fact(obj, key, fact):
    obj[key] = fact


def add_charge(p, key, raw_value, normalized, page):
    p["charges"][key] = direct(raw_value, normalized, page)


def participant(pid, label, role, page, ref=None):
    p = st.new_participant(pid, mid(label), role, evid(page))
    if role == "product":
        p["product_role"] = direct("reported product", "main_product", page)
    if ref:
        p["equivalence_reference"] = direct(f"relative to {ref}", mid(ref), page)
    return p


def productive_classification(r, page):
    changes = [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 1, "after_order": 2, "change_type": "order_increased"},
        {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
        {"change_id": "b4", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}], "before_order": 1, "after_order": 0, "change_type": "broken"}
    ]
    why = "Net changes were derived from the source reactant/product drawings without a complete atom mapping; atom identities and hydrogen transfers remain partial."
    r["classification"]["bond_changes"] = [derived(c, page, why) for c in changes]
    r["classification"]["labels"] = derived([
        {"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1"]},
        {"category_number": 13, "code": "DB-C-C", "supporting_change_ids": ["b2"]}
    ], page, why)
    r["classification"]["coverage"] = inferred("partial", page, why)
    r["classification"]["limitations"] = [
        "No complete atom-mapped reaction SMILES was asserted.",
        "Implicit-H changes and eliminated carbonate coproducts are unresolved.",
        "Classification is a structure derivation, not an author-provided label."
    ]


def no_covalent_classification(r, page, reason):
    r["classification"]["bond_changes"] = []
    r["classification"]["labels"] = inferred([], page, reason)
    r["classification"]["coverage"] = inferred("not_applicable", page, reason)
    r["classification"]["limitations"] = [reason]


def add_outcome(r, oid, kind, targets, page, raw_value, normalized, *, yield_type=None, basis=None, scope=None):
    o = st.new_outcome(oid, kind, [mid(t) for t in targets], evid(page))
    o["value"] = direct(raw_value, normalized, page)
    if kind in {"yield", "product_amount"}:
        o["product_role"] = direct("reported product", "main_product", page)
    if yield_type is not None:
        o["yield_type"] = direct(yield_type, yield_type, page)
    if basis is not None:
        o["measurement_basis"] = direct(basis, basis, page)
    if scope is not None:
        o["scope"] = direct(scope, scope, page)
    r["outcomes"].append(o)
    return o


def make_reaction(suffix, page, categories, substrates, product=None, *, scale=0.2, allene_equiv=1.1,
                  solvent="toluene", silver="AgSbF6", temp="r.t.", time_h=24, variation=None,
                  source_name="rhodium-catalyzed 2-C–H 1,3-dienylation", volume_override="auto"):
    r = st.new_reaction(rid(suffix), DOC, evid(page))
    r["record_type"] = "specific_experiment"
    r["experiment_categories"] = direct(", ".join(categories), categories, page)
    r["location"] = direct(f"PDF page {page}, {suffix}", {"page": page, "region": suffix}, page)
    r["reaction_name"] = direct(source_name, source_name, page)
    r["transformation_type"] = inferred("indole C–H dienylation with an allenyl carbonate", page, "Normalized from the paper title and source reaction drawing.")
    ps = []
    p1 = participant("p_substrate_1", substrates[0], "reactant", page)
    add_charge(p1, "amount", f"{scale} mmol", {"number": scale, "unit": "mmol"}, page)
    ps.append(p1)
    if len(substrates) > 1:
        p2 = participant("p_substrate_2", substrates[1], "reactant", page, substrates[0])
        add_charge(p2, "equivalents", f"{allene_equiv} equiv.", {"number": allene_equiv, "unit": "equiv"}, page)
        ps.append(p2)
    if product:
        ps.append(participant("p_product", product, "product", page))
    cat = participant("p_rh_catalyst", "RhCat", "catalyst", page, substrates[0])
    add_charge(cat, "mol_percent", "2.5 mol%", {"number": 2.5, "unit": "mol%"}, page)
    ps.append(cat)
    if silver:
        ag = participant("p_silver", silver, "additive", page, substrates[0])
        add_charge(ag, "mol_percent", "30 mol%", {"number": 30, "unit": "mol%"}, page)
        ps.append(ag)
    piv = participant("p_pivoh", "PivOH", "additive", page, substrates[0])
    add_charge(piv, "equivalents", "1 equiv.", {"number": 1, "unit": "equiv"}, page)
    ps.append(piv)
    solv = participant("p_solvent", solvent, "solvent", page)
    volume = (5.0 if page in {3, 6} else 1.0) if volume_override == "auto" else volume_override
    if volume is not None:
        add_charge(solv, "volume", f"{volume:g} mL", {"number": volume, "unit": "mL"}, page)
    ps.append(solv)
    r["participants"] = ps
    r["reaction_scale"] = direct(f"{scale} mmol of {substrates[0]}", {"reference_molecule_id": mid(substrates[0]), "amount": {"number": scale, "unit": "mmol"}}, page)
    r["solvent_conditions"]["composition"] = direct(solvent, [mid(solvent)], page)
    if volume is not None:
        r["solvent_conditions"]["total_volume"] = direct(f"{volume:g} mL", {"number": volume, "unit": "mL"}, page)
    step = st.new_step("reaction", evid(page))
    step["participant_ids"] = [p["participant_id"] for p in ps if p["role"] != "product"]
    op = "combined under the conditions shown in the table/scheme"
    step["operation"] = direct(op, op, page)
    if isinstance(temp, (int, float)):
        step["temperature"] = direct(f"{temp} °C", {"number": temp, "unit": "°C"}, page)
    else:
        step["temperature"] = direct("room temperature" if temp == "r.t." else temp, temp, page)
    if time_h is not None:
        step["time"] = direct(f"{time_h:g} h", {"number": time_h, "unit": "h"}, page)
    step["atmosphere"] = direct("Ar", "argon", page)
    r["steps"] = [step]
    caption = {
        2: "The reactions were carried out using the quantities and conditions shown in Table 1/Table 3 or the displayed equation.",
        3: "The reactions were carried out using 1 mmol of 1, 1.1 equiv. of 2, 2.5 mol% of [Cp*RhCl2]2, 30 mol% of AgSbF6, and 1 equiv. of PivOH in 5 mL of toluene.",
        6: "The reactions were carried out using 1 mmol of 1a, 1.0 equiv. of 2, 2.5 mol% of [Cp*RhCl2]2, 30 mol% of AgSbF6, and 1 equiv. of PivOH in 5 mL of toluene."
    }.get(page, "conditions shown in the source")
    r["procedure"]["original_text"] = direct(caption, caption, page)
    r["procedure"]["general_procedure_reference"] = direct(f"page {page} table/scheme conditions", f"page {page} table/scheme conditions", page)
    if variation:
        r["procedure"]["modifications"] = direct(variation, variation, page)
    r["review"]["unresolved"] = [
        "No explicitly associated SI was supplied; addition order, workup, purification details, and complete characterization are unreported here.",
        "Concrete structures not encoded for every source-labeled scope compound require expert structure review.",
        "Machine transcription and derived bond classification require review."
    ]
    r["completeness"] = "partial"
    if product:
        productive_classification(r, page)
    return r


reactions = []

# Scheme 1b is a generic transformation, not a concrete experiment.
scheme = st.new_reaction(rid("scheme1b_generic"), DOC, evid(1))
scheme["record_type"] = "generic_scheme"
scheme["experiment_categories"] = direct("generic reaction concept", ["other"], 1)
scheme["location"] = direct("PDF page 1, Scheme 1b", {"page": 1, "region": "Scheme 1b"}, 1)
scheme["reaction_name"] = direct("C–H functionalization of indoles with 4-aryl-2,3-allenyl carbonates", "rhodium-catalyzed indole dienylation", 1)
scheme["transformation_type"] = inferred("generic indole C–H dienylation", 1, "Normalized from Scheme 1b; no concrete substrate instance is asserted.")
productive_classification(scheme, 1)
scheme["notes"] = [direct("Scheme 1b illustrates the proposed general transformation and directing-group π–π-stacking concept; it is not one concrete experiment.", "generic scheme only", 1)]
scheme["review"]["unresolved"] = ["Generic R/Ar/X variables are not expanded into concrete molecules in this record."]
reactions.append(scheme)

# Table 1 optimization rows; page-1 prose is additional evidence for entry 1.
for row in raw["table1"]:
    product = "3aa" if row.get("yield") is not None else None
    r = make_reaction(f"table1_entry{row['entry']}", 2, ["optimization", "control"], ["1a", "2a"], product,
                      solvent=row["solvent"], silver=row["silver"], variation=row["variation"])
    if row["entry"] == 1:
        r["evidence_ids"] = list(dict.fromkeys(r["evidence_ids"] + evid(1)))
        r["notes"].append(direct("Page 1 prose reports the same optimized reaction as Table 1 entry 1.", "duplicate source mention for same experiment", 1))
    if row.get("recovery") is not None:
        add_outcome(r, "o_recovery", "recovery", ["1a"], 2, f"{row['recovery']} under NMR recovery of 1a (%)", {"number": row["recovery"], "unit": "%"}, yield_type="recovery", basis="1H NMR using CH3NO2 as internal standard")
    if row.get("yield") is not None:
        add_outcome(r, "o_yield", "yield", ["3aa"], 2, f"{row['yield']} under NMR yield of 3aa (%)", {"number": row["yield"], "unit": "%"}, yield_type="NMR yield", basis="1H NMR using CH3NO2 as internal standard")
        add_outcome(r, "o_selectivity", "selectivity", ["3aa"], 2, f"Z/E of 3aa: {row['z_e']}", row["z_e"], basis="1H NMR analysis of crude reaction mixture before chromatography", scope={"metric": "Z/E", "order": ["Z", "E"]})
    else:
        add_outcome(r, "o_no_product", "no_reaction", [], 2, "table dash; prose states no expected product", True, basis="expected-product observation; detection limit not reported")
        no_covalent_classification(r, 2, "No expected product was observed; no bond-formation class is assigned to a negative control.")
    reactions.append(r)

# Equation 1.
eq1 = raw["equations_page2"][0]
r = make_reaction("eq1_3am", 2, ["substrate_scope"], eq1["reactants"], eq1["product"], scale=eq1["scale_mmol"], allene_equiv=eq1["allene_equiv"], time_h=eq1["time_h"], variation="secondary 4-aryl-2,3-butadienyl carbonate 2m", volume_override=None)
add_outcome(r, "o_yield", "yield", ["3am"], 2, "3am, 81%", {"number": 81, "unit": "%"}, yield_type="reported yield", basis=None)
for i, sel in enumerate(eq1["selectivity"], 1):
    add_outcome(r, f"o_selectivity_{i}", "selectivity", ["3am"], 2, f"{sel['order']} of {sel['bond']} bond = {sel['ratio']}", sel["ratio"], basis="source equation; analytical method not stated at the equation", scope={"metric": sel["order"], "bond": sel["bond"], "order": sel["order"].split("/")})
reactions.append(r)

# Equation 2 is an isomerization/stability control with no asserted connectivity change.
eq2 = raw["equations_page2"][1]
r = make_reaction("eq2_3ah_isomerization_control", 2, ["control", "mechanistic_study"], ["3ah_ZE98_2"], "3ah_ZE95_5", scale=0.1, allene_equiv=1.0, temp=70, time_h=24, variation="Z-3ah enriched by preparative HPLC was treated under the standard condition", volume_override=0.5)
add_outcome(r, "o_yield", "yield", ["3ah_ZE95_5"], 2, "3ah, NMR yield = 78%", {"number": 78, "unit": "%"}, yield_type="NMR yield", basis="NMR; internal standard not stated for this equation")
add_outcome(r, "o_selectivity", "selectivity", ["3ah_ZE95_5"], 2, "Z/E = 95/5", "95/5", basis="source equation; method not stated", scope={"metric": "Z/E", "order": ["Z", "E"], "input_ratio": "98/2"})
no_covalent_classification(r, 2, "This control reports a change in geometric-isomer ratio, not a new connectivity or bond-order transformation.")
reactions.append(r)

# Table 2 substrate scope, including the separately reported 4 mmol scale experiment.
for row in raw["table2"]:
    suffix = f"table2_{row['product']}" + (f"_{row['record_suffix']}" if row.get("record_suffix") else "")
    scale = row.get("scale_mmol", 1.0)
    r = make_reaction(suffix, 3, ["substrate_scope"] + (["scale_up"] if row.get("record_suffix") else []), row["substrates"], row["product"], scale=scale, temp=row.get("temperature_c", "r.t."), time_h=row.get("time_h", 24))
    add_outcome(r, "o_yield", "yield", [row["product"]], 3, f"{row['product']}, {row['yield']}%", {"number": row["yield"], "unit": "%"}, yield_type="isolated yield", basis="isolated product")
    add_outcome(r, "o_selectivity", "selectivity", [row["product"]], 3, f"Z/E = {row['z_e']}", row["z_e"], basis="1H NMR analysis", scope={"metric": "Z/E", "order": ["Z", "E"]})
    if row.get("product_mass_g") is not None:
        add_outcome(r, "o_product_mass", "product_amount", [row["product"]], 3, f"{row['product_mass_g']} g", {"number": row["product_mass_g"], "unit": "g"}, basis="isolated product mass")
    if row.get("assignment"):
        r["notes"].append(direct(row["assignment"], row["assignment"], 3))
    reactions.append(r)

# Table 3 control rows.
for row in raw["table3"]:
    r = make_reaction(f"table3_entry{row['entry']}", 4, ["control", "mechanistic_study"], row["substrates"], row.get("product"), scale=row.get("scale_mmol", 0.2), temp=row.get("temperature_c", "r.t."), time_h=row.get("time_h", 24), variation=("absence of light" if row.get("absence_of_light") else None))
    if row.get("absence_of_light"):
        r["special_conditions"]["light_source"] = direct("reaction was run in the absence of light", "dark", 4)
    if row.get("recovery") is not None:
        add_outcome(r, "o_recovery", "recovery", [row["substrates"][0]], 4, f"NMR recovery of 1: {row['recovery']}%", {"number": row["recovery"], "unit": "%"}, yield_type="recovery", basis="1H NMR using CH3NO2 as internal standard")
    if row.get("yield") is not None:
        add_outcome(r, "o_nmr_yield", "yield", [row["product"]], 4, f"NMR yield of 3: {row['yield']}%", {"number": row["yield"], "unit": "%"}, yield_type="NMR yield", basis="1H NMR using CH3NO2 as internal standard")
        if row.get("isolated_yield") is not None:
            add_outcome(r, "o_isolated_yield", "yield", [row["product"]], 4, f"isolated yield {row['isolated_yield']}%", {"number": row["isolated_yield"], "unit": "%"}, yield_type="isolated yield", basis="isolated product")
        add_outcome(r, "o_selectivity", "selectivity", [row["product"]], 4, f"Z/E = {row['z_e']}", row["z_e"], basis="1H NMR analysis", scope={"metric": "Z/E", "order": ["Z", "E"]})
    else:
        add_outcome(r, "o_no_product", "no_reaction", [], 4, "table dash; prose states no product was generated/failed to afford product", True, basis="expected-product observation; detection limit not reported")
        no_covalent_classification(r, 4, "No product was reported in this negative control; no intended-reaction bond labels are assigned.")
    reactions.append(r)

# Computational mechanism figures are separate proposed-mechanism records.
for page, fig, description in [
    (4, "fig1_dft_pathway", "DFT free-energy profile for C–H activation, carborhodation and beta-oxygen elimination"),
    (5, "fig2_dft_ze", "DFT comparison of transition states proposed to control Z/E selectivity for 1a + 2a"),
    (5, "fig3_dft_aryl_alkyl", "DFT comparison of aryl- and alkyl-allene pathways leading to opposite major geometry")
]:
    r = st.new_reaction(rid(fig), DOC, evid(page))
    r["record_type"] = "proposed_mechanism"
    r["experiment_categories"] = direct("computational mechanistic analysis", ["other"], page)
    r["location"] = direct(f"PDF page {page}, {fig.split('_')[0].title()}", {"page": page, "region": fig}, page)
    r["reaction_name"] = direct(description, description, page)
    r["transformation_type"] = inferred("proposed computational mechanism", page, "The source labels these regions as calculated free-energy/transition-state analyses, not preparative experiments.")
    r["notes"] = [direct(description, description, page), direct("Not counted as a synthetic experiment.", "proposed mechanism only", page)]
    r["review"]["unresolved"] = ["Computational structures/energies are retained as source context; no experimental reaction yield is assigned."]
    r["classification"]["coverage"] = inferred("unresolved", page, "No single experimentally observed atom-mapped net reaction is represented by this transition-state figure.")
    reactions.append(r)

# Table 4 alkyl-allene scope.
for row in raw["table4"]:
    r = make_reaction(f"table4_{row['product']}", 6, ["substrate_scope"], row["substrates"], row["product"], scale=row.get("scale_mmol", 1.0), allene_equiv=row.get("allene_equiv", 1.0))
    basis = "isolated product"
    ytype = "isolated yield"
    if row.get("yield_note"):
        basis = "1H NMR analysis of E-3 and E/Z mixture after column chromatography"
        ytype = "calculated yield after chromatography"
        r["notes"].append(direct(row["yield_note"], row["yield_note"], 6))
        r["review"]["unresolved"].append("An unidentified unknown product was reported for 3ap; its identity and amount are not given.")
    add_outcome(r, "o_yield", "yield", [row["product"]], 6, f"{row['product']}, {row['yield']}%", {"number": row["yield"], "unit": "%"}, yield_type=ytype, basis=basis)
    o = add_outcome(r, "o_selectivity", "selectivity", [row["product"]], 6, f"E/Z = {row['e_z']}", row["e_z"], basis="1H NMR analysis", scope=None)
    alt_ez = direct("product label E/Z", {"metric": "E/Z", "order": ["E", "Z"]}, 6)
    alt_ze = direct("caption says ratio of Z/E", {"metric": "Z/E", "order": ["Z", "E"]}, 6)
    o["scope"] = st.fact(["product labels: E/Z", "caption: Z/E"], None, evid(6), status="conflict", alternatives=[alt_ez, alt_ze])
    r["notes"].append(direct("Source-order conflict: product labels show E/Z, while the Table 4 caption says the ratio of Z/E was determined.", "selectivity-order conflict preserved", 6))
    r["review"]["unresolved"].append("Resolve Table 4 E/Z versus caption Z/E ordering against the authors' SI or an expert reading; values were not silently reordered.")
    if row.get("assignment"):
        r["notes"].append(direct(row["assignment"], row["assignment"], 6))
    reactions.append(r)

bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": rid("page4_fig1_computational_context"), "document_id": DOC, "kind": "DFT_energy_profile", "evidence_ids": evid(4), "note": "Fig. 1 is a calculated pathway, represented separately as a proposed_mechanism record."},
    {"context_id": rid("page5_dft_context"), "document_id": DOC, "kind": "DFT_selectivity_analysis", "evidence_ids": evid(5), "note": "Figs. 2 and 3 compare calculated transition states and are not synthetic experiments."},
    {"context_id": rid("page6_fig4_lumo"), "document_id": DOC, "kind": "computed_orbitals", "evidence_ids": evid(6), "note": "Fig. 4 shows calculated LUMOs for (R_a)-2a and (R_a)-2n as examples; no experimental enantiopurity is inferred."},
    {"context_id": rid("si_limit"), "document_id": DOC, "kind": "source_limitation", "evidence_ids": evid(1), "note": "The article states ESI exists, but no SI PDF was explicitly supplied or associated in this workspace."}
]

# Update prepared dataset metadata/tasks before merge; no experiment records are injected here.
(RUN / "dataset.json").write_text(json.dumps(dataset, ensure_ascii=False, indent=2) + "\n")
tasks = json.loads((RUN / "tasks.json").read_text())
for task in tasks:
    if task["document_id"] == DOC and 1 <= int(task["evidence_id"].split(":page")[1].split(":")[0]) <= 6:
        task["state"] = "attempted"
        task["reason"] = "Page image visually inspected and represented in coverage; unresolved items remain explicit."
(RUN / "tasks.json").write_text(json.dumps(tasks, ensure_ascii=False, indent=2) + "\n")
(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]), "source_context": len(bundle["source_context"])}, ensure_ascii=False))
