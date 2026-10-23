#!/usr/bin/env python3
"""Build the source-derived v3.1 submission for the four-page ChemComm paper.

This script only assembles facts transcribed into raw_model/extraction_observations.json
and visually checked against the prepared page images. It does not modify the extraction
tools, schema, or quality rules.
"""
from copy import deepcopy
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / "runs/current"
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_b91b082127244b8b"
E1T, E1I = f"{DOC}:page1:text", f"{DOC}:page1:image"
E2T, E2I = f"{DOC}:page2:text", f"{DOC}:page2:image"
E3T, E3I = f"{DOC}:page3:text", f"{DOC}:page3:image"
E4T, E4I = f"{DOC}:page4:text", f"{DOC}:page4:image"


def extracted(raw, normalized, ev):
    return st.fact(raw, normalized, ev, method="original_extraction")


def inherited(raw, normalized, ev):
    return st.fact(raw, normalized, ev, method="inherited_general_condition")


def inferred(normalized, ev, why, raw=None):
    return st.fact(raw, normalized, ev, method="inference", derivation=why)


def derived(normalized, ev, why, raw=None):
    return st.fact(raw, normalized, ev, method="structure_derivation", derivation=why)


def missing(status="unreported", ev=()):
    return st.fact(evidence=ev, status=status)


def conflict(raw, alternatives, ev):
    return st.fact(raw, None, ev, status="conflict", alternatives=alternatives)


def q(number, unit, **extra):
    return {"number": number, "unit": unit, **extra}


dataset = json.loads((RUN / "dataset.json").read_text())
raw = json.loads((RUN / "raw_model/extraction_observations.json").read_text())
bundle = st.empty(dataset["documents"])
molecules = {}
stereo_diagnostics = {}


def add_molecule(label, ev, name=None, smiles=None, kind="unresolved", note=None,
                 chiral="unresolved", geometry="unresolved", representation="unresolved",
                 stereo_assignment=None, axial=None, mixture=None, occurrence_labels=None):
    mid = f"{DOC}:mol:{label.replace('/', '_')}"
    m = st.new_molecule(mid, DOC, ev)
    labels = occurrence_labels or label
    m["occurrences"][0]["compound_label"] = extracted(labels, labels, ev)
    m["occurrences"][0]["evidence_ids"] = list(ev)
    m["names"]["original"] = extracted(name or label, name or label, ev)
    m["names"]["english"] = extracted(name, name, ev) if name else missing()
    m["names"]["abbreviations"] = extracted([label], [label], ev)
    m["names"]["aliases"] = missing()
    m["cas_numbers"] = missing()
    m["identity_note"] = extracted(note, note, ev) if note else missing()
    m["identity_status"] = "source_identified" if name or label else "unresolved"
    m["representation_type"] = kind
    if smiles:
        diag = st.stereo_inspect(smiles)
        stereo_diagnostics[label] = diag
        m["smiles"] = derived(smiles, ev, "SMILES candidate derived from the explicit source name/structure label; RDKit parsing does not verify the source drawing.", raw=name or label)
        m["stereochemistry"]["isomeric_smiles"] = derived(diag["canonical_isomeric_smiles"], ev, "Canonical isomeric SMILES computed by standard stereo inspection from the candidate SMILES.")
        m["stereochemistry"]["identity_key"] = derived(diag["stereo_graph_key"], ev, "Stereo graph key computed by standard stereo inspection; not a sample identity.")
    else:
        m["smiles"] = missing()
    m["stereochemistry"]["chirality_status"] = inferred(chiral, ev, "Status retained from the visible drawing/text; absence of an assignment is not treated as achiral or racemic.")
    m["stereochemistry"]["geometry_status"] = inferred(geometry, ev, "Status retained from the visible drawing/text; unlabelled alkene geometry is not assigned.")
    m["stereochemistry"]["representation_status"] = inferred(representation, ev, "No stereochemically complete concrete structure is asserted unless explicitly encoded in the candidate SMILES.")
    if stereo_assignment is not None:
        m["stereochemistry"]["assignment"] = extracted(stereo_assignment, stereo_assignment, ev)
    if axial is not None:
        m["stereochemistry"]["axial"] = extracted(axial, axial, ev)
    if mixture is not None:
        m["stereochemistry"]["mixture"] = extracted(mixture, mixture, ev)
    molecules[label] = m
    bundle["molecules"].append(m)
    return mid


# Source-labelled substrates. The exact main-text structures are not forced into SMILES.
scope_substitutions = {
    "1a": "aldimine ester 1a (R1=H, R2=Et)", "1b": "aldimine ester 1b (R1=4-Me)",
    "1c": "aldimine ester 1c (R1=4-OMe)", "1d": "aldimine ester 1d (R1=4-OCF3)",
    "1e": "aldimine ester 1e (R1=4-CH2CO2Me)", "1f": "aldimine ester 1f (R1=4-F)",
    "1g": "aldimine ester 1g (R1=4-Cl)", "1h": "aldimine ester 1h (R1=4-Br)",
    "1i": "aldimine ester 1i (R1=5-N(Et)2)", "1j": "aldimine ester 1j (R1=5-OMe)",
    "1k": "aldimine ester 1k (R1=6-F)", "1l": "aldimine ester 1l (R1=6-OMe)",
    "1m": "aldimine ester 1m (R1=6-OEt)", "1n": "aldimine ester 1n (R1=6-Me)",
    "1o": "aryl-substituted aldimine ester 1o", "1p": "methylenedioxy aldimine ester 1p",
    "1q": "dichloro aldimine ester 1q", "1r": "naphthyl aldimine ester 1r",
    "1s": "aldimine ester 1s (para-CO2Me)", "1t": "aldimine ester 1t (para-CF3)",
    "1u": "aldimine ester 1u (meta-CO2Me)", "1v": "Bpin-substituted aldimine ester 1v",
    "1w": "nitro-substituted aldimine ester 1w", "1x": "aldimine ester 1x (R2=Me)",
    "1y": "aldimine ester 1y (R2=iPr)",
}
for label, name in scope_substitutions.items():
    ev = [E2I, E2T] + ([E3I, E3T] if label in {"1a", "1r", "1v", "1w", "1x"} else [])
    add_molecule(label, ev, name=name, note="Concrete stereochemically complete SMILES not asserted from the compact scope drawing.",
                 chiral="unresolved", geometry="unspecified", representation="not_encoded",
                 occurrence_labels="1r (written as 1ra in mechanism discussion)" if label == "1r" else label)

for label, ester in {"2a": "Bn", "2b": "Me", "2c": "Et", "2d": "iPr"}.items():
    ev = [E2I, E2T] + ([E3I, E3T] if label == "2a" else [])
    add_molecule(label, ev, name=f"γ-vinyl allenoate {label} (R3={ester})",
                 note="Allene-containing substrate; the paper does not assign an allene axial configuration.",
                 chiral="unspecified", geometry="unspecified", representation="not_encoded",
                 axial={"type": "allene", "configuration": "not assigned by source"})

product_labels = [row[2] for row in raw["scope_entries"]]
for label in dict.fromkeys(product_labels):
    ev = [E2I, E2T] + ([E3I, E3T] if label in {"3aa", "3ra", "3sa", "3ta", "3ua", "3va", "3wa", "3xa", "3ya", "3ab", "3ac", "3ad"} else []) + ([E4I, E4T] if label in {"3aa", "3ra"} else [])
    assignment = {
        "type": "relative configuration",
        "source": "wedge/dash drawing in Table 1",
        "absolute_configuration": "not assigned",
        "additional_support": "3aa relative configuration confirmed by X-ray crystallographic analysis" if label == "3aa" else None,
    }
    mix = "Depicted major diastereomer; Table 1 footnote reports dr >20:1 for new products." if label != "3wa" else "Expected product drawing only; product was not obtained."
    add_molecule(label, ev, name=f"chromeno[4,3-b]pyrrole product {label}",
                 note="Three contiguous stereogenic centers are reported; absolute configurations are not assigned in the supplied pages.",
                 chiral="partially_specified", geometry="unspecified", representation="not_encoded",
                 stereo_assignment=assignment, mixture=mix)

# Control, labelled, and derivatization compounds.
add_molecule("6", [E3I, E3T], name="ethyl (E)-2-[(2-hydroxybenzylidene)amino]acetate 6",
             note="The text explicitly names the imine as E; no complete SMILES is asserted here.",
             chiral="unresolved", geometry="specified", representation="not_encoded",
             stereo_assignment={"imine_geometry": "E (explicitly named)"})
add_molecule("7a", [E3I, E3T], name="amido analogue 7a (R=H)", chiral="unresolved", geometry="unspecified", representation="not_encoded")
add_molecule("7b", [E3I, E3T], name="4-methylbenzenesulfonamido analogue 7b (R=p-Ts)", chiral="unresolved", geometry="unspecified", representation="not_encoded")
add_molecule("8", [E3I, E3T], name="methoxy analogue 8", chiral="unresolved", geometry="unspecified", representation="not_encoded")
add_molecule("d-3aa", [E3I, E3T], name="deuterium-labelled product d-3aa",
             note="Partial D incorporation at C9, C15 and C16; this is an isotopologue mixture, not a single fully labelled structure.",
             chiral="partially_specified", geometry="unspecified", representation="not_encoded",
             stereo_assignment={"type": "relative configuration depicted; absolute configuration not assigned"},
             mixture={"type": "isotopologue mixture", "D_incorporation": {"C9": "26%", "C15": "82%", "C16": "82%"}})
for label in ("4", "5"):
    add_molecule(label, [E3I, E3T], name=f"derivatization product {label}",
                 note="Relative stereochemical drawing retained; absolute configuration and concrete isomeric SMILES are not asserted.",
                 chiral="partially_specified", geometry="unspecified", representation="not_encoded",
                 stereo_assignment={"type": "relative wedge/dash configuration in Scheme 2", "absolute_configuration": "not assigned"})

# Named reagents and solvents with conservative candidate structures.
simple = {
    "pFPh3P": ("tris(4-fluorophenyl)phosphine", "P(c1ccc(F)cc1)(c1ccc(F)cc1)c1ccc(F)cc1", [E2I, E2T, E3I, E3T]),
    "PPh3": ("triphenylphosphine", "P(c1ccccc1)(c1ccccc1)c1ccccc1", [E2T, E2I]),
    "PhCOOH": ("benzoic acid", "O=C(O)c1ccccc1", [E2I, E2T, E3I, E3T]),
    "PhCOONa": ("sodium benzoate", "O=C([O-])c1ccccc1.[Na+]", [E2I, E2T, E3I, E3T]),
    "toluene": ("toluene", "Cc1ccccc1", [E2I, E2T, E3I, E3T]),
    "D2O": ("deuterium oxide", "[2H]O[2H]", [E3I, E3T]),
    "H2": ("hydrogen", "[H][H]", [E3I, E3T]),
    "MeOH": ("methanol", "CO", [E3I, E3T]),
    "H2SO4": ("sulfuric acid", "OS(=O)(=O)O", [E3I, E3T]),
    "LiCl": ("lithium chloride", "[Li+].[Cl-]", [E3I]),
    "H2O": ("water", "O", [E3I]),
    "DMF": ("N,N-dimethylformamide", "CN(C)C=O", [E3I]),
}
for label, (name, smiles, ev) in simple.items():
    add_molecule(label, ev, name=name, smiles=smiles, kind="smiles",
                 note="Candidate structure derived from the explicit reagent/solvent name; not a claim about purity or grade.",
                 chiral="achiral", geometry="not_applicable", representation="not_applicable")
add_molecule("Pd_C", [E3I, E3T], name="palladium on carbon (Pd/C)",
             note="Heterogeneous named catalyst; support composition and loading are not reported, so no unique concrete structure is asserted.",
             chiral="not_applicable", geometry="not_applicable", representation="not_applicable")


def participant(pid, label, role, ev, reference=None, eq=None, molpct=None, mass=None, amount=None, stage="reaction", product_role=None):
    p = st.new_participant(pid, molecules[label]["molecule_id"], role, ev, stage=stage)
    if product_role:
        p["product_role"] = extracted(product_role, product_role, ev)
    if mass is not None:
        p["charges"]["mass"] = mass if isinstance(mass, dict) and "data_status" in mass else extracted(str(mass), mass, ev)
    if amount is not None:
        p["charges"]["amount"] = extracted(f"{amount['number']} {amount['unit']}", amount, ev)
    if eq is not None:
        p["charges"]["equivalents"] = extracted(f"{eq} equiv.", q(eq, "equiv"), ev)
    if molpct is not None:
        p["charges"]["mol_percent"] = extracted(f"{molpct} mol%", q(molpct, "mol%"), ev)
    if reference:
        p["equivalence_reference"] = inferred(molecules[reference]["molecule_id"], ev, "The source loading is referenced to the principal aldimine/starting substrate for this reaction.", raw=f"relative to {reference}")
    return p


def base_reaction(suffix, record_type, ev, categories=None, location=None, reaction_name=None, transformation=None, unresolved=None):
    rid = f"{DOC}:{suffix}"
    r = st.new_reaction(rid, DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = extracted(categories, categories, ev) if categories is not None else missing()
    r["location"] = extracted(location, location, ev) if location is not None else missing()
    r["reaction_name"] = extracted(reaction_name, reaction_name, ev) if reaction_name else missing()
    r["transformation_type"] = inferred(transformation, ev, "Transformation type normalized from the source reaction name and depicted transformation.") if transformation else missing()
    r["review"]["unresolved"] = list(unresolved or [])
    r["completeness"] = "partial"
    return r


def set_step(step, ev, operation=None, temperature=None, time=None, inherited_values=False):
    if operation is not None:
        step["operation"] = (inherited if inherited_values else extracted)(operation, operation, ev)
    if temperature is not None:
        raw_temp, norm_temp = temperature
        step["temperature"] = (inherited if inherited_values else extracted)(raw_temp, norm_temp, ev)
    if time is not None:
        raw_time, norm_time = time
        step["time"] = (inherited if inherited_values else extracted)(raw_time, norm_time, ev)


def set_solvent(r, label, ev, inherited_value=False):
    maker = inherited if inherited_value else extracted
    r["solvent_conditions"]["composition"] = maker(label, [label], ev)


def set_annulation_classification(r, ev):
    why = "The article states that the cascade constructs three sigma bonds; the depicted mechanism supports two new C-C connections and one new C-O connection. Exact atom mapping is unavailable, so coverage remains partial."
    changes = [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
    ]
    r["classification"]["bond_changes"] = [derived(x, ev, why) for x in changes]
    labels = [
        {"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1", "b2"]},
        {"category_number": 3, "code": "C-O", "supporting_change_ids": ["b3"]},
    ]
    r["classification"]["labels"] = derived(labels, ev, why)
    r["classification"]["coverage"] = inferred("partial", ev, why)
    r["classification"]["limitations"] = ["No complete atom mapping is asserted.", "Hydrogen transfers and catalyst-bound intermediates are not fully mapped."]


def set_no_reaction_classification(r, ev):
    why = "The source reports that the corresponding product was not obtained/no reaction; no expected-product bond formation is assigned."
    r["classification"]["bond_changes"] = []
    r["classification"]["labels"] = inferred([], ev, why)
    r["classification"]["coverage"] = inferred("not_applicable", ev, why)
    r["classification"]["limitations"] = ["Negative control; intended transformation is not classified as observed chemistry."]


def yield_outcome(oid, product_label, value, ev, basis=None, scope=None):
    o = st.new_outcome(oid, "yield", [molecules[product_label]["molecule_id"]], ev)
    o["product_role"] = extracted("target product", "target product", ev)
    o["value"] = extracted(f"{value}%", q(value, "%"), ev)
    if basis:
        o["yield_type"] = extracted(basis, basis, ev)
        o["measurement_basis"] = extracted(basis, basis, ev)
    if scope:
        o["scope"] = extracted(scope, scope, ev)
    return o


def amount_outcome(oid, product_label, amount, ev):
    o = st.new_outcome(oid, "product_amount", [molecules[product_label]["molecule_id"]], ev)
    o["product_role"] = extracted("target product", "target product", ev)
    o["value"] = extracted(f"{amount['number']} {amount['unit']}", amount, ev)
    o["measurement_basis"] = extracted("reported isolated product mass", "reported product mass", ev)
    return o


def dr_outcome(oid, product_label, ev):
    o = st.new_outcome(oid, "selectivity", [molecules[product_label]["molecule_id"]], ev)
    o["product_role"] = extracted("depicted major diastereomer", "major diastereomer", ev)
    o["value"] = extracted("dr >20:1", {"metric": "dr", "qualifier": "gt", "major_minor_ratio": "20:1"}, ev)
    o["scope"] = extracted("diastereomeric ratio; major/minor absolute configurations not named", {"metric": "dr", "order": "depicted major : other diastereomers", "absolute_configuration": "unassigned"}, ev)
    return o


GENERAL_ID = f"{DOC}:scheme1_current_general"

# Scheme 1 prior work and this-work generic transformation.
for idx, name in [(1, "authors' previous phosphine-catalyzed (3+2)/(2+3) sequential annulation"),
                  (2, "authors' previous phosphine-catalyzed (3+2)/(3+2) annulation")]:
    r = base_reaction(f"scheme1_previous_eq{idx}", "literature_example", [E2I, E2T], location=f"page 2, Scheme 1, previous work Eq. {idx}", reaction_name=name,
                      unresolved=["Generic R groups and cited-paper experimental details were not expanded; this is not a current-paper experiment."])
    r["procedure"]["original_text"] = extracted("(p-MeOC6H4)3P, toluene, 110 °C", "(p-MeOC6H4)3P; toluene; 110 °C", [E2I])
    s = st.new_step("depicted_reaction", [E2I])
    set_step(s, [E2I], operation="literature reaction depicted in Scheme 1", temperature=("110 °C", q(110, "°C")))
    r["steps"] = [s]
    bundle["reactions"].append(r)

r = base_reaction("scheme1_current_general", "generic_scheme", [E2I, E2T], location="page 2, Scheme 1, This work",
                  reaction_name="phosphine-catalyzed cascade (2+3)/(2+4) annulation",
                  transformation="cascade annulation to chromeno[4,3-b]pyrroles",
                  unresolved=["Generic R1/R2/R3 scheme is not a concrete experiment; instance conditions are recorded separately."])
set_annulation_classification(r, [E2I, E4T])
bundle["reactions"].append(r)

# Main-text PPh3 screening experiment (SI itself is unavailable).
r = base_reaction("optimization_pph3_3aa_50pct", "specific_experiment", [E2T, E2I], ["optimization"],
                  "page 2 text; Table S1 entry 1 is cited but SI not supplied",
                  "model reaction of 1a and 2a", "phosphine-catalyzed cascade annulation",
                  ["Reaction time, reactant charges, yield measurement basis, workup and purification are not reported in the supplied main text.", "Referenced Table S1 was not supplied and was not inferred."])
r["general_scheme_ids"] = [GENERAL_ID]
r["participants"] = [
    participant("p_1a", "1a", "reactant", [E2T]), participant("p_2a", "2a", "reactant", [E2T]),
    participant("p_3aa", "3aa", "product", [E2T], product_role="target product"),
    participant("p_PPh3", "PPh3", "catalyst", [E2T], reference="1a", molpct=30),
    participant("p_toluene", "toluene", "solvent", [E2T]),
]
set_solvent(r, "toluene", [E2T])
s = st.new_step("reaction", [E2T]); s["participant_ids"] = [p["participant_id"] for p in r["participants"]]
set_step(s, [E2T], operation="1a and 2a reacted with PPh3 in toluene", temperature=("25 °C", q(25, "°C")))
r["steps"] = [s]
r["outcomes"] = [yield_outcome("yield_3aa", "3aa", 50, [E2T], scope="PPh3 screening entry")]
set_annulation_classification(r, [E2T, E2I])
bundle["reactions"].append(r)


def add_scope_reaction(row, index):
    sub, allene, prod, value, descriptor = row
    is_buffer = sub == "1v" and value == 41
    is_no = prod == "3wa"
    suffix = f"table1_{prod}_{'buffer' if is_buffer else 'standard'}"
    categories = ["substrate_scope"] + (["optimization"] if prod == "3aa" else [])
    unresolved = ["Reactant charges, solvent volume/concentration, atmosphere, pressure, workup and purification are not reported in the supplied pages."]
    if prod != "3aa":
        unresolved.append("The table gives a shared 2–4 h range, not a product-specific exact time.")
    r = base_reaction(suffix, "specific_experiment", [E2I, E2T] + ([E3T] if prod in {"3sa", "3ta", "3ua", "3va", "3wa", "3xa", "3ya", "3ab", "3ac", "3ad"} else []),
                      categories, f"page 2, Table 1 entry {index}: {prod}",
                      "phosphine-catalyzed cascade (2+3)/(2+4) annulation",
                      "cascade annulation to chromeno[4,3-b]pyrrole", unresolved)
    r["general_scheme_ids"] = [GENERAL_ID]
    ev = [E2I]
    r["participants"] = [
        participant("p_substrate", sub, "reactant", ev), participant("p_allenoate", allene, "reactant", ev),
        participant("p_product", prod, "product", ev, product_role="expected target product" if is_no else "target product"),
        participant("p_catalyst", "pFPh3P", "catalyst", ev, reference=sub, eq=0.3),
        participant("p_PhCOOH", "PhCOOH", "additive", ev, reference=sub, eq=1.0),
        participant("p_toluene", "toluene", "solvent", ev),
    ]
    if prod == "3aa":
        r["participants"][3]["charges"]["mol_percent"] = extracted("30 mol%", q(30, "mol%"), [E2T])
    if is_buffer:
        r["participants"].append(participant("p_PhCOONa", "PhCOONa", "additive", ev, reference=sub, eq=1.0))
        r["procedure"]["modifications"] = extracted("Added PhCOONa (1.0 equiv.) and PhCOOH (1.0 equiv.) as buffer.", "PhCOONa/PhCOOH buffer, 1.0 equiv/1.0 equiv", ev)
    r["procedure"]["original_text"] = inherited("Unless otherwise noted, all reactions were carried out under the optimal conditions.", raw["scope_conditions"], [E2I])
    r["procedure"]["general_procedure_reference"] = inherited("optimal conditions", "Table 1 optimal conditions", [E2I])
    set_solvent(r, "toluene", ev, inherited_value=True)
    s = st.new_step("reaction", ev); s["participant_ids"] = [p["participant_id"] for p in r["participants"]]
    time = ("2 h", q(2, "h")) if prod == "3aa" else ("2–4 h", {"minimum": 2, "maximum": 4, "unit": "h", "raw_text": "2–4 h"})
    set_step(s, [E2T, E2I] if prod == "3aa" else ev, operation="reaction under optimal conditions", temperature=("25 °C", q(25, "°C")), time=time, inherited_values=True)
    r["steps"] = [s]
    if is_no:
        o = st.new_outcome("no_product_3wa", "no_reaction", [], [E2I, E3T])
        o["value"] = extracted("3wa, 0; expected product could not be obtained", True, [E2I, E3T])
        o["scope"] = extracted("expected nitro-substituted product 3wa", prod, [E2I, E3T])
        r["outcomes"] = [o]
        set_no_reaction_classification(r, [E2I, E3T])
    else:
        r["outcomes"] = [yield_outcome(f"yield_{prod}", prod, value, [E2I], basis="isolated yield", scope=descriptor), dr_outcome(f"dr_{prod}", prod, [E2I])]
        set_annulation_classification(r, [E2I, E4T])
    bundle["reactions"].append(r)


for i, row in enumerate(raw["scope_entries"], 1):
    add_scope_reaction(row, i)


def standard_annulation_participants(sub, prod, ev, mass_fact=None, amount=None):
    return [
        participant("p_substrate", sub, "reactant", ev, mass=mass_fact, amount=amount),
        participant("p_allenoate", "2a", "reactant", ev),
        participant("p_product", prod, "product", ev, product_role="target product"),
        participant("p_catalyst", "pFPh3P", "catalyst", ev, reference=sub, eq=0.3),
        participant("p_PhCOOH", "PhCOOH", "additive", ev, reference=sub, eq=1.0),
        participant("p_toluene", "toluene", "solvent", ev),
    ]


# Scheme 2a: two distinct gram-scale experiments.
mass_1a = conflict({"prose": "1.0 g", "scheme": "1.26 g"},
                   [extracted("1.0 g", q(1.0, "g"), [E3T]), extracted("1.26 g", q(1.26, "g"), [E3I])], [E3T, E3I])
for sub, prod, mass, mmol, pmass, pyield in [
    ("1a", "3aa", mass_1a, 3.6, 1.4, 65), ("1x", "3xa", extracted("1.0 g", q(1.0, "g"), [E3T, E3I]), 4.0, 1.14, 71)
]:
    r = base_reaction(f"scheme2a_gram_scale_{prod}", "specific_experiment", [E3I, E3T], ["scale_up"],
                      f"page 3, Scheme 2a, gram-scale {prod}",
                      "gram-scale phosphine-catalyzed cascade annulation", "cascade annulation to chromeno[4,3-b]pyrrole",
                      ["Amount of 2a, workup, purification, atmosphere, pressure and yield measurement basis are not reported in the supplied pages."] +
                      (["Source conflict: prose reports 1a 1.0 g, while Scheme 2a depicts 1a 1.26 g; both values are retained without selection."] if sub == "1a" else []))
    r["general_scheme_ids"] = [GENERAL_ID]
    r["participants"] = standard_annulation_participants(sub, prod, [E3I, E3T], mass, q(mmol, "mmol"))
    r["reaction_scale"] = extracted(f"{sub}: {mmol} mmol", {"reference_molecule_id": molecules[sub]["molecule_id"], "amount": q(mmol, "mmol")}, [E3T])
    set_solvent(r, "toluene", [E3I])
    s = st.new_step("reaction", [E3I]); s["participant_ids"] = [p["participant_id"] for p in r["participants"]]
    set_step(s, [E3I], operation="gram-scale annulation", temperature=("25 °C", q(25, "°C")), time=("2 h", q(2, "h")))
    r["steps"] = [s]
    r["outcomes"] = [amount_outcome(f"amount_{prod}", prod, q(pmass, "g"), [E3T, E3I]), yield_outcome(f"yield_{prod}", prod, pyield, [E3T, E3I], scope=f"gram-scale {sub} to {prod}")]
    set_annulation_classification(r, [E3I, E4T])
    bundle["reactions"].append(r)

# Scheme 2b, one-pot hydrogenation/debenzylation/Fischer esterification.
r = base_reaction("scheme2b_3xa_to_4", "specific_experiment", [E3I, E3T], ["derivatization"], "page 3, Scheme 2b, 3xa to 4",
                  "heterogeneous catalytic hydrogenation, debenzylation and Fischer esterification", "one-pot reduction/deprotection/esterification",
                  ["Reagent loadings (Pd/C and H2SO4), temperature, time, pressure, workup, purification and yield measurement basis are not reported.", "Complete atom mapping across hydrogenation, debenzylation and esterification is unavailable."])
r["participants"] = [
    participant("p_3xa", "3xa", "reactant", [E3I, E3T]), participant("p_4", "4", "product", [E3I, E3T], product_role="target product"),
    participant("p_PdC", "Pd_C", "catalyst", [E3I, E3T]), participant("p_H2", "H2", "reagent", [E3I, E3T]),
    participant("p_MeOH", "MeOH", "solvent", [E3I, E3T]), participant("p_H2SO4", "H2SO4", "catalyst", [E3I, E3T]),
]
set_solvent(r, "MeOH", [E3I])
s1 = st.new_step("hydrogenation_and_debenzylation", [E3I, E3T]); s1["participant_ids"] = ["p_3xa", "p_PdC", "p_H2", "p_MeOH"]
set_step(s1, [E3I, E3T], operation="The carbon-carbon double bond was reduced by heterogeneous catalytic hydrogenation and the benzyl group was removed at the same time.")
s2 = st.new_step("fischer_esterification", [E3I, E3T]); s2["participant_ids"] = ["p_3xa", "p_MeOH", "p_H2SO4", "p_4"]
set_step(s2, [E3I, E3T], operation="followed by addition of catalytic amounts of sulfuric acid to afford Fischer esterification product 4")
r["steps"] = [s1, s2]
r["procedure"]["addition_order"] = extracted("Pd/C, H2, MeOH, then H2SO4", ["Pd/C", "H2", "MeOH", "then H2SO4"], [E3I])
r["outcomes"] = [yield_outcome("yield_4", "4", 96, [E3I, E3T], scope="one-pot 3xa to 4")]
why = "The text explicitly reports reduction of a C=C bond by hydrogenation; the net C-H gains are classified, while debenzylation/esterification atom origins remain unresolved."
changes = [
    {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
    {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 2},
]
r["classification"]["bond_changes"] = [derived(x, [E3I, E3T], why) for x in changes]
r["classification"]["labels"] = derived([{"category_number": 1, "code": "C-H", "supporting_change_ids": ["b2"]}], [E3I, E3T], why)
r["classification"]["coverage"] = inferred("partial", [E3I, E3T], why)
r["classification"]["limitations"] = ["No full atom mapping; ester oxygen and benzyl-group atom fates are unresolved."]
bundle["reactions"].append(r)

# Scheme 2b, Krapcho decarboxylation.
r = base_reaction("scheme2b_4_to_5", "specific_experiment", [E3I, E3T], ["derivatization"], "page 3, Scheme 2b, 4 to 5",
                  "Krapcho decarboxylation", "decarboxylation",
                  ["Reaction time, solvent volume, workup, purification and yield measurement basis are not reported.", "Exact atom mapping and proton source are unresolved."])
r["participants"] = [
    participant("p_4", "4", "reactant", [E3I, E3T]), participant("p_5", "5", "product", [E3I, E3T], product_role="target product"),
    participant("p_LiCl", "LiCl", "reagent", [E3I], reference="4", eq=2.0), participant("p_H2O", "H2O", "reagent", [E3I], reference="4", eq=2.0),
    participant("p_DMF", "DMF", "solvent", [E3I]),
]
set_solvent(r, "DMF", [E3I])
s = st.new_step("krapcho_decarboxylation", [E3I, E3T]); s["participant_ids"] = [p["participant_id"] for p in r["participants"]]
set_step(s, [E3I], operation="Krapcho decarboxylation with LiCl and H2O in DMF", temperature=("140 °C", q(140, "°C")))
r["steps"] = [s]
r["outcomes"] = [yield_outcome("yield_5", "5", 35, [E3I, E3T], scope="4 to 5")]
why = "Krapcho decarboxylation removes a carbon substituent and gives net replacement by hydrogen; only the supported net C-H formation and C-C cleavage are retained."
changes = [
    {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 1, "after_order": 0, "change_type": "broken"},
    {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain"},
]
r["classification"]["bond_changes"] = [derived(x, [E3I, E3T], why) for x in changes]
r["classification"]["labels"] = derived([{"category_number": 1, "code": "C-H", "supporting_change_ids": ["b2"]}], [E3I, E3T], why)
r["classification"]["coverage"] = inferred("partial", [E3I, E3T], why)
r["classification"]["limitations"] = ["Proton source and leaving fragment are not fully mapped."]
bundle["reactions"].append(r)


def add_control(label, eq_no, wording):
    r = base_reaction(f"scheme3_control_eq{eq_no}_{label}", "specific_experiment", [E3I, E3T], ["control", "mechanistic_study"],
                      f"page 3, Scheme 3, Eq. {eq_no}, substrate {label}", "control experiment under optimal conditions", "negative-control annulation",
                      ["No product identity/yield method, reactant charges, workup or purification are reported; no product is fabricated."])
    r["general_scheme_ids"] = [GENERAL_ID]
    r["participants"] = [
        participant("p_control", label, "reactant", [E3I, E3T]), participant("p_2a", "2a", "reactant", [E3I, E3T]),
        participant("p_catalyst", "pFPh3P", "catalyst", [E3I, E2I], reference=label, eq=0.3),
        participant("p_PhCOOH", "PhCOOH", "additive", [E3I, E2I], reference=label, eq=1.0),
        participant("p_toluene", "toluene", "solvent", [E3I, E2I]),
    ]
    set_solvent(r, "toluene", [E2I], inherited_value=True)
    s = st.new_step("reaction", [E3I, E2I]); s["participant_ids"] = [p["participant_id"] for p in r["participants"]]
    set_step(s, [E3I, E2I], operation="optimal conditions", temperature=("25 °C", q(25, "°C")), time=("2 h", q(2, "h")), inherited_values=True)
    r["steps"] = [s]
    o = st.new_outcome(f"no_reaction_{label}", "no_reaction", [], [E3I, E3T])
    o["value"] = extracted(wording, True, [E3I, E3T])
    o["scope"] = extracted(f"control substrate {label} with 2a", {"substrate": label, "partner": "2a"}, [E3I])
    r["outcomes"] = [o]
    set_no_reaction_classification(r, [E3I, E3T])
    bundle["reactions"].append(r)


add_control("6", 1, "corresponding product not obtained")
add_control("7a", 2, "desired product not obtained (R=H)")
add_control("7b", 2, "desired product not obtained (R=p-Ts)")
add_control("8", 3, "no reaction")

# Scheme 3 deuterium-labelling experiment.
r = base_reaction("scheme3_deuterium_label_d3aa", "specific_experiment", [E3I, E3T], ["mechanistic_study"],
                  "page 3, Scheme 3, Eq. 4", "deuterium-labelling experiment under standard conditions", "H/D exchange during cascade annulation",
                  ["Reactant charges, workup, purification and yield measurement basis are not reported.", "D incorporation is partial and the exact isotopologue distribution is not reported."])
r["general_scheme_ids"] = [GENERAL_ID]
r["participants"] = [
    participant("p_1a", "1a", "reactant", [E3I, E3T]), participant("p_2a", "2a", "reactant", [E3I, E3T]),
    participant("p_D2O", "D2O", "reagent", [E3I, E3T], reference="1a", eq=10.0),
    participant("p_d3aa", "d-3aa", "product", [E3I, E3T], product_role="deuterium-labelled target product"),
    participant("p_catalyst", "pFPh3P", "catalyst", [E2I, E3I], reference="1a", eq=0.3),
    participant("p_PhCOOH", "PhCOOH", "additive", [E2I, E3I], reference="1a", eq=1.0),
    participant("p_toluene", "toluene", "solvent", [E2I, E3I]),
]
set_solvent(r, "toluene", [E2I], inherited_value=True)
s = st.new_step("reaction", [E3I, E2I]); s["participant_ids"] = [p["participant_id"] for p in r["participants"]]
set_step(s, [E3I, E2I], operation="standard conditions with D2O (10 equiv.)", temperature=("25 °C", q(25, "°C")), time=("2 h", q(2, "h")), inherited_values=True)
r["steps"] = [s]
r["outcomes"] = [yield_outcome("yield_d3aa", "d-3aa", 58, [E3I, E3T], scope="deuterium-labelling experiment")]
for pos, pct in [("C9", 26), ("C15", 82), ("C16", 82)]:
    o = st.new_outcome(f"D_incorporation_{pos}", "selectivity", [molecules["d-3aa"]["molecule_id"]], [E3I, E3T])
    o["product_role"] = extracted("deuterated product d-3aa", "target isotopologue mixture", [E3I])
    o["value"] = extracted(f"{pct}% D", q(pct, "%"), [E3I])
    o["measurement_basis"] = extracted("1H NMR analysis", "1H NMR", [E3T])
    o["scope"] = extracted(f"H/D exchange at {pos}", {"metric": "D incorporation", "position": pos}, [E3I, E3T])
    r["outcomes"].append(o)
why = "Scheme 3 reports partial replacement of C-H by C-D at C9, C15 and C16. Deuterium is classified under the C-H taxonomy category, with isotope retained on the H atom; fractional incorporation and proton sources make coverage partial."
changes = []
formed_ids = []
for i, pos in enumerate(("C9", "C15", "C16"), 1):
    fid, bid = f"b{2*i-1}", f"b{2*i}"
    formed_ids.append(fid)
    changes.append({"change_id": fid, "atoms": [{"atom_map": None, "element": "C", "source_position": pos}, {"atom_map": None, "element": "H", "isotope": 2}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain"})
    changes.append({"change_id": bid, "atoms": [{"atom_map": None, "element": "C", "source_position": pos}, {"atom_map": None, "element": "H", "isotope": 1}], "before_order": 1, "after_order": 0, "change_type": "hydrogen_loss"})
r["classification"]["bond_changes"] = [derived(x, [E3I, E3T], why) for x in changes]
r["classification"]["labels"] = derived([{"category_number": 1, "code": "C-H", "supporting_change_ids": formed_ids}], [E3I, E3T], why)
r["classification"]["coverage"] = inferred("partial", [E3I, E3T], why)
r["classification"]["limitations"] = ["Partial isotopic incorporation; individual isotopologue abundances and complete atom mapping are unavailable."]
bundle["reactions"].append(r)


# Figure 2 proposed mechanisms are records, not experimental proof.
for suffix, sub, prod, branch, ev in [
    ("figure2_main_mechanism", "1a", "3aa", "main catalytic cycle via I–VII", [E3I, E3T, E4T]),
    ("figure2_naphthyl_branch", "1r", "3ra", "naphthyl branch via IV-1ra, [1,5] H-shift and o-quinone methide VIII", [E3I, E4I, E4T]),
]:
    r = base_reaction(suffix, "proposed_mechanism", ev, None, f"page 3 Figure 2 and page 4 continuation: {branch}",
                      "putative mechanism of the cascade annulation", "proposed catalytic mechanism",
                      ["Author-proposed mechanism; intermediates and elementary steps are not treated as experimentally isolated reactions.", "No kinetic or intermediate-isolation evidence is reported on the supplied pages."])
    r["general_scheme_ids"] = [GENERAL_ID]
    r["participants"] = [
        participant("p_substrate", sub, "reactant", ev), participant("p_2a", "2a", "reactant", ev),
        participant("p_product", prod, "product", ev, product_role="proposed-cycle product"),
        participant("p_catalyst", "pFPh3P", "catalyst", ev), participant("p_PhCOOH", "PhCOOH", "additive", ev),
    ]
    operations = ([
        "nucleophilic phosphine addition to 2a followed by protonation by PhCOOH",
        "remote 1,7-addition of enolate I to acceptor II to form the first sigma bond",
        "intramolecular carbon-anion addition to vinyl phosphonium",
        "phenolate addition/formal (2+4) cycloaddition",
        "phosphine elimination to product 3aa and catalyst regeneration",
    ] if sub == "1a" else [
        "remote 1,7-addition to resonated intermediate IV-1ra",
        "proposed [1,5] H-shift and dearomatization to o-quinone methide VIII",
        "subsequent annulation and phosphine elimination to 3ra",
    ])
    r["steps"] = []
    for i, op in enumerate(operations, 1):
        s = st.new_step(f"proposed_{i}", ev); s["participant_ids"] = [p["participant_id"] for p in r["participants"]]
        set_step(s, ev, operation=op)
        r["steps"].append(s)
    r["procedure"]["original_text"] = extracted(branch, branch, ev)
    set_annulation_classification(r, ev)
    bundle["reactions"].append(r)


# Non-synthetic observations and document context.
bundle["source_context"] = [
    {"context_id": f"{DOC}:context:figure1", "document_id": DOC, "kind": "non_synthetic_figure",
     "evidence_ids": [E1I, E1T], "description": "Figure 1 presents the chromeno[4,3-b]pyrrole core and biological/natural-product examples; no reaction arrow or current-paper experiment is present."},
    {"context_id": f"{DOC}:context:xray_3aa", "document_id": DOC, "kind": "structural_assignment",
     "evidence_ids": [E1T, E2T], "description": "CCDC 2075585 is associated with 3aa; the text states that X-ray crystallography confirmed its structure and relative configuration."},
    {"context_id": f"{DOC}:context:references", "document_id": DOC, "kind": "bibliography",
     "evidence_ids": [E4I, E4T], "description": "Page 4 bibliography is not extracted as current-paper experiments."},
]

# Save the raw standard-stereo diagnostics used for every concrete candidate SMILES.
st.write(RUN / "raw_model/stereo_diagnostics.json", stereo_diagnostics)
st.write(RUN / "submitted.json", bundle)

ids = [r["reaction_id"] for r in bundle["reactions"]]
page2_ids = [x for x in ids if any(k in x for k in ("scheme1_", "optimization_", "table1_"))]
page3_ids = [x for x in ids if any(k in x for k in ("table1_3s", "table1_3t", "table1_3u", "table1_3v", "table1_3w", "table1_3x", "table1_3y", "table1_3a", "scheme2", "scheme3", "figure2"))]
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title, abstract and bibliographic header", "finding": "Current paper identity and transformation summary; no concrete experimental row.", "evidence_ids": [E1I, E1T]},
            {"region": "Figure 1", "finding": "Core/natural-product structures and biological context only; non-synthetic source_context, not a reaction.", "evidence_ids": [E1I, E1T]},
            {"region": "introductory text and ESI/CCDC footnote", "finding": "Prior literature discussion; CCDC 2075585 identifies 3aa structural evidence. ESI availability is stated, but no SI file was supplied or associated.", "evidence_ids": [E1I, E1T]}
        ], "record_ids": [f"{DOC}:table1_3aa_standard"], "unresolved": ["No SI was supplied; cited supporting procedures and crystallographic file were not inspected."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Scheme 1, previous work Eqs. 1–2", "finding": "Two literature examples, explicitly separated from this paper's experiments.", "evidence_ids": [E2I, E2T]},
            {"region": "Scheme 1, This work", "finding": "Generic cascade-annulation scheme, recorded as generic_scheme rather than a specific experiment.", "evidence_ids": [E2I, E2T]},
            {"region": "left-column model-reaction and optimization prose", "finding": "PPh3/30 mol% model experiment affords 3aa in 50%; optimal p-fluorophenylphosphine/PhCOOH conditions afford 3aa in 82% and >20:1 dr.", "evidence_ids": [E2I, E2T]},
            {"region": "Table 1 header, 29 displayed outcomes/condition variants and footnotes", "finding": "All table entries inventoried individually; 3va standard and buffer are separate records, and 3wa is a negative outcome rather than an invented product experiment.", "evidence_ids": [E2I, E2T]},
            {"region": "right-column scope discussion", "finding": "Text cross-checks substituent groups and yields; continuation occurs on page 3.", "evidence_ids": [E2I, E2T]}
        ], "record_ids": page2_ids, "unresolved": ["Individual substrate charges and exact per-entry times (beyond shared 2–4 h) are not shown.", "Concrete stereochemically complete SMILES are not asserted from compact scope drawings.", "Table S1 and other SI-only optimization rows are unavailable."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "top-left continuation of Table 1 discussion", "finding": "Confirms 3sa/3ta/3ua, buffered 3va, failed 3wa, and 3xa/3ya/3ab–3ad outcomes.", "evidence_ids": [E3I, E3T]},
            {"region": "Scheme 2a gram-scale synthesis", "finding": "Two scale-up experiments; 1a mass conflict retained (text 1.0 g vs scheme 1.26 g).", "evidence_ids": [E3I, E3T]},
            {"region": "Scheme 2b transformation", "finding": "3xa→4 one-pot hydrogenation/debenzylation/esterification and 4→5 Krapcho decarboxylation recorded separately.", "evidence_ids": [E3I, E3T]},
            {"region": "Scheme 3 control experiments Eqs. 1–3", "finding": "Four negative controls recorded separately (6, 7a, 7b, 8); no product fabricated.", "evidence_ids": [E3I, E3T]},
            {"region": "Scheme 3 Eq. 4 deuterium labelling", "finding": "d-3aa, 58%; C9 26% D, C15/C16 82% D; 1H NMR basis retained for D incorporation.", "evidence_ids": [E3I, E3T]},
            {"region": "Figure 2", "finding": "Author-proposed catalytic cycle and naphthyl branch recorded as proposed_mechanism, not experimental proof.", "evidence_ids": [E3I, E3T]}
        ], "record_ids": page3_ids, "unresolved": ["Most procedures, workups, purification, times and yield measurement bases require unavailable SI.", "Scale-up 1a mass conflict cannot be resolved from the supplied paper.", "Mechanism intermediates are proposed and not independently verified."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "mechanism continuation, left column", "finding": "Completes phosphine elimination/catalyst regeneration and describes separate naphthyl [1,5]-H-shift/o-quinone-methide proposal.", "evidence_ids": [E4I, E4T]},
            {"region": "conclusion", "finding": "Author summary: three sigma bonds and three contiguous stereogenic centers; no new concrete experiment.", "evidence_ids": [E4I, E4T]},
            {"region": "references", "finding": "Bibliography only; no current-paper experiment created.", "evidence_ids": [E4I, E4T]}
        ], "record_ids": [f"{DOC}:figure2_main_mechanism", f"{DOC}:figure2_naphthyl_branch", f"{DOC}:table1_3ra_standard"], "unresolved": ["The energetic rationale for the naphthyl pathway is an author proposal, not a measured intermediate or barrier."]}
    ]
}
(RUN / "coverage").mkdir(exist_ok=True)
st.write(RUN / "coverage/pages_0001_0004.json", coverage)
print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]), "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in bundle["reactions"]), "stereo_diagnostics": len(stereo_diagnostics)}, ensure_ascii=False))
