#!/usr/bin/env python3
"""Build source-derived v3.1 records for adsc.202000150, PDF pages 1-6."""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_65d1cdf1cb09f398"
RUN = ROOT / "runs/current"
DATASET = st.read(RUN / "dataset.json")
EVT = lambda p: f"{DOC}:page{p}:text"
EVI = lambda p: f"{DOC}:page{p}:image"


def write(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n")


def f(raw=None, norm=None, ev=(), method="original_extraction", status=None, derivation=None):
    return st.fact(raw, norm, ev, method=method, status=status, derivation=derivation)


def missing(ev=(), status="unreported"):
    return f(ev=ev, status=status)


def inferred(norm, ev, derivation, raw=None):
    return f(raw, norm, ev, method="inference", derivation=derivation)


def derived(norm, ev, derivation, raw=None):
    return f(raw, norm, ev, method="structure_derivation", derivation=derivation)


def set_status(m, chirality, geometry, representation, ev, *, assignment=None,
               axial=None, mixture=None, atom_centers=None, double_bonds=None):
    s = m["stereochemistry"]
    s["chirality_status"] = inferred(chirality, ev, "Status preserves what the cited drawing/text does and does not specify.")
    s["geometry_status"] = inferred(geometry, ev, "Status preserves what the cited drawing/text does and does not specify.")
    s["representation_status"] = inferred(representation, ev, "No complete stereochemical structure is claimed unless an isomeric SMILES is supplied.")
    if assignment is not None:
        s["assignment"] = f(assignment, assignment, ev)
    if axial is not None:
        s["axial"] = f(axial, axial, ev)
    if mixture is not None:
        s["mixture"] = f(mixture, mixture, ev)
    if atom_centers is not None:
        s["atom_centers"] = f(atom_centers, atom_centers, ev)
    if double_bonds is not None:
        s["double_bonds"] = f(double_bonds, double_bonds, ev)


bundle = st.empty(DATASET["documents"])
molecules = {}


def add_unresolved(mid, label, ev, name=None, note=None, *, chirality="unresolved",
                   geometry="unresolved", representation="not_encoded", assignment=None,
                   axial=None, mixture=None):
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["representation_type"] = "unresolved"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    m["names"]["original"] = f(name or label, name or label, ev)
    m["names"]["english"] = f(name, name, ev) if name else missing(ev)
    m["names"]["abbreviations"] = f(label, [label], ev)
    m["identity_note"] = f(note or "Identity is source-linked by label/name; compact structure was not encoded as a complete SMILES.",
                           note or "Identity is source-linked by label/name; compact structure was not encoded as a complete SMILES.", ev)
    set_status(m, chirality, geometry, representation, ev, assignment=assignment, axial=axial, mixture=mixture)
    molecules[mid] = m
    return m


def add_template(mid, label, ev, fixed, variables, note, *, chirality="unresolved",
                 geometry="unresolved", axial=None):
    m = add_unresolved(mid, label, ev, name=label, note=note, chirality=chirality,
                       geometry=geometry, representation="not_encoded", axial=axial)
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": f(fixed, fixed, ev),
        "variables": [{
            "label": v,
            "attachment_points": missing(ev, "unreadable"),
            "definition": missing(ev, "unreadable"),
            "allowed_values": missing(ev),
            "constraints": missing(ev),
        } for v in variables],
        "instance_mapping": [],
        "reason": f(note, note, ev),
    }
    return m


stereo_diagnostics = {}


def add_smiles(mid, label, name, smiles, ev):
    diag = st.stereo_inspect(smiles)
    stereo_diagnostics[mid] = diag
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["representation_type"] = "smiles"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    m["names"]["original"] = f(name, name, ev)
    m["names"]["english"] = f(name, name, ev)
    m["names"]["abbreviations"] = f(label, [label], ev)
    m["smiles"] = derived(smiles, ev, "Candidate SMILES derived from the explicitly named simple reagent; parsing is not source verification.", raw=name)
    s = m["stereochemistry"]
    s["chirality_status"] = derived("achiral", ev, "RDKit graph inspection found no potential stereocenter in the candidate simple-reagent graph.")
    s["geometry_status"] = derived("not_applicable", ev, "No applicable configurable alkene geometry occurs in the candidate simple-reagent graph.")
    s["representation_status"] = derived("encoded", ev, "The candidate simple-reagent graph contains no unresolved stereochemical element.")
    s["isomeric_smiles"] = derived(diag["canonical_isomeric_smiles"], ev, "Computed by standard stereo from the candidate SMILES.")
    s["identity_key"] = derived(diag["stereo_graph_key"], ev, "Computed by standard stereo; graph key is not sample identity.")
    m["identity_note"] = f("Simple named reagent; candidate graph is machine-derived and remains needs_review.",
                           "Simple named reagent; candidate graph is machine-derived and remains needs_review.", ev)
    molecules[mid] = m
    return m


# Generic substrates and products.
add_template("m_isoindigo_1", "isoindigos 1", [EVI(2), EVI(3), EVI(5)],
             "Isoindigo-derived C4 synthon with PG/Ar substituent variables as drawn in Scheme 1 and Tables 2/4.",
             ["PG", "Ar1", "Ar2"], "Generic R/PG structure; no wildcard concrete SMILES asserted.")
add_template("m_allenoate_2", "allenoates 2", [EVI(2), EVI(3), EVI(4)],
             "Electron-deficient allene C2 synthon bearing ester/substituent variables.",
             ["R", "R1", "R2", "FG", "FG'"], "Generic allenoate family; substituent mapping is not fully legible in every compact row.",
             chirality="unspecified", geometry="not_applicable",
             axial={"element": "allene axis", "configuration": "not specified in source"})
add_template("m_pyrazoloneyldiene_4", "pyrazoloneyldiene oxindoles 4", [EVI(4), EVT(4), EVT(5)],
             "Pyrazoloneyldiene oxindole C4 synthon with Ar1/Het/R1 variables.",
             ["Ar1", "Het", "R1"], "Generic scope structure in Table 3; no wildcard concrete SMILES asserted.")
add_template("m_product_3_generic", "products 3", [EVI(2), EVI(3), EVI(5)],
             "Fused/spirocyclic 4H-pyran annulation product with variable PG/aryl/allenoate substituents.",
             ["PG", "Ar1", "Ar2", "R1", "R2"], "Generic product family; stereochemistry varies by experiment.",
             chirality="partially_specified", geometry="unresolved")
add_template("m_product_5_generic", "products 5", [EVI(4), EVT(5)],
             "Pyrano[2,3-c]pyrazole annulation product with variable oxindole/pyrazolone/allenoate substituents.",
             ["Ar1", "Het", "R1", "R2"], "Generic product family; stereochemistry varies by experiment.",
             chirality="partially_specified", geometry="unresolved")
molecules["m_product_3_generic"]["stereochemistry"]["assignment"] = f(
    "Generic product drawings show relative stereochemical bonds, but absolute configuration varies or is unassigned.",
    {"relative_depiction": True, "absolute_configuration": "not generally assigned"}, [EVI(2), EVI(3), EVI(5)])
molecules["m_product_5_generic"]["stereochemistry"]["assignment"] = f(
    "Generic product drawings show relative stereochemical bonds, but absolute configuration varies or is unassigned.",
    {"relative_depiction": True, "absolute_configuration": "not generally assigned"}, [EVI(4), EVT(5)])

# Specifically reused labeled molecules.
add_unresolved("m_1a", "1a", [EVI(3), EVI(4), EVI(5)], "N-Boc isoindigo 1a",
               "Specific source structure is legible, but a complete SMILES was not safely transcribed.",
               chirality="unspecified", geometry="specified", representation="not_encoded",
               assignment={"central_alkene": "configuration depicted; descriptor not printed"})
add_unresolved("m_2a", "2a", [EVI(3), EVI(4), EVI(5)], "allenoate 2a",
               "Specific electron-deficient allene structure is source-depicted; axial configuration is not assigned.",
               chirality="unspecified", geometry="not_applicable", representation="not_encoded",
               axial={"element": "allene axis", "configuration": "unassigned"})
add_unresolved("m_1h", "1h", [EVI(4)], "brominated isoindigo 1h",
               chirality="unspecified", geometry="specified", assignment={"central_alkene": "configuration depicted; descriptor not printed"})
add_unresolved("m_1_pg_me", "1 (PG = Me)", [EVI(5), EVT(5)], "N-methyl isoindigo control substrate",
               chirality="unspecified", geometry="specified", assignment={"central_alkene": "configuration depicted; descriptor not printed"})

product3_labels = [f"3{c}" for c in "abcdefghijkl"]
for label in product3_labels:
    ev = [EVI(3)] + ([EVT(4), EVI(4)] if label == "3h" else [])
    assignment = {"relative_stereochemistry": "wedge/dash product depiction; major/minor order not fully named"}
    chirality = "partially_specified"
    if label == "3h":
        assignment = {"absolute_configuration": "(3'S, 4S)", "method": "X-ray analysis", "CCDC": "1937754"}
        chirality = "specified"
    add_unresolved(f"m_{label}", label, ev, f"annulation product {label}",
                   "Compact product drawing retained by label; no complete stereochemical SMILES was safely transcribed.",
                   chirality=chirality, geometry="unresolved", representation="not_encoded", assignment=assignment)

for label in [f"5{c}" for c in "abcdefghijk"]:
    assignment = {"relative_stereochemistry": "wedge/dash product depiction; major/minor order not fully named"}
    chirality = "partially_specified"
    if label == "5c":
        assignment = {"absolute_configuration": "(3'S, 4S)", "method": "X-ray analysis", "CCDC": "1968335"}
        chirality = "specified"
    add_unresolved(f"m_{label}", label, [EVI(4), EVT(4)], f"annulation product {label}",
                   "Compact product drawing retained by label; no complete stereochemical SMILES was safely transcribed.",
                   chirality=chirality, geometry="unresolved", representation="not_encoded", assignment=assignment)

add_unresolved("m_3_pg_me", "3 (PG = Me)", [EVI(5), EVT(5)], "N-methyl annulation product",
               chirality="racemic", geometry="unresolved", representation="not_encoded",
               mixture={"enantiomeric_excess": "0% ee; racemic version stated"},
               assignment={"relative_stereochemistry": "90:10 dr; absolute configurations not assigned"})
add_unresolved("m_6", "6", [EVI(4), EVT(4)], "hydrolysed compound 6",
               chirality="partially_specified", geometry="unresolved", representation="not_encoded",
               assignment={"relationship": "stereoselectivity retained from 3h; >20:1 dr, 92% ee"})
add_unresolved("m_d3a", "d-3a", [EVI(4), EVT(4)], "deuterium-labelled product d-3a",
               chirality="partially_specified", geometry="unresolved", representation="not_encoded",
               assignment={"label_site": "3'-position", "incorporation": "30% D"})

# Catalysts are independently registered by source label; compact drawings are not retyped.
for label in [f"P{i}" for i in range(1, 11)] + ["P6-1", "P6-2"]:
    ev = [EVI(3)] if label.startswith("P") and label not in ("P6-1", "P6-2") else [EVI(5)]
    if label in ("P6", "P10"):
        ev = sorted(set(ev + [EVI(4), EVI(5)]))
    add_unresolved(f"m_cat_{label.lower().replace('-', '_')}", label, ev,
                   f"chiral phosphonium salt catalyst {label}",
                   "Catalyst drawing/source label retained; complete stereochemical structure not encoded.",
                   chirality="specified", geometry="unresolved", representation="not_encoded",
                   assignment={"source": "chiral amino-acid/dipeptide-derived catalyst drawing"})

# Simple reagents and named process materials.
add_smiles("m_cs2co3", "Cs2CO3", "cesium carbonate", "[Cs+].[Cs+].[O-]C([O-])=O", [EVI(3), EVI(4), EVI(5)])
add_smiles("m_meoh", "MeOH", "methanol", "CO", [EVI(5), EVT(5)])
add_smiles("m_fecl3", "FeCl3", "ferric chloride", "[Fe+3].[Cl-].[Cl-].[Cl-]", [EVI(4), EVT(4)])
add_smiles("m_d2o", "D2O", "deuterium oxide", "[2H]O[2H]", [EVI(4), EVT(4)])
add_smiles("m_dcm", "CH2Cl2", "dichloromethane", "ClCCl", [EVI(4)])
add_smiles("m_etoac", "EtOAc", "ethyl acetate", "CCOC(C)=O", [EVT(5)])
add_template("m_xylene", "xylene", [EVI(3), EVI(4), EVI(5), EVT(5)],
             "Xylene solvent; positional-isomer composition is not stated.", ["isomer composition"],
             "Named solvent of unspecified isomer composition; not forced to one concrete SMILES.",
             chirality="not_applicable", geometry="not_applicable")
add_template("m_petroleum_ether", "petroleum ether", [EVT(5)],
             "Petroleum ether hydrocarbon mixture; composition not reported.", ["mixture composition"],
             "Named mixture of unknown composition; no unique concrete molecule asserted.",
             chirality="not_applicable", geometry="not_applicable")
add_template("m_silica_gel", "silica gel", [EVT(5)],
             "Silica gel stationary phase/material.", ["surface/composition"],
             "Process material rather than a unique discrete molecular identity.",
             chirality="not_applicable", geometry="not_applicable")

reactions = []


def base_reaction(suffix, ev, record_type="specific_experiment", categories=None, location=None,
                  name="[4 + 2] annulation", transformation="enantioselective [4 + 2] annulation"):
    rid = f"{DOC}:{suffix}"
    r = st.new_reaction(rid, DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, ev) if categories is not None else missing(ev, "not_applicable")
    r["location"] = f(location, location, ev)
    r["reaction_name"] = f(name, name, ev)
    r["transformation_type"] = inferred(transformation, ev, "Normalized from the authors' reaction name and depicted net transformation.")
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
                   "unresolved": ["No explicitly associated SI was supplied; compound-specific SI details remain unavailable."]}
    r["completeness"] = "partial"
    reactions.append(r)
    return r


def participant(r, pid, mid, role, ev, stage="reaction", *, amount=None, mass=None, volume=None,
                equivalents=None, mol_percent=None, ref=None, raw_equiv=None):
    p = st.new_participant(pid, mid, role, ev, stage=stage)
    if amount is not None:
        amount_raw = {0.10: "0.10", 0.40: "0.40", 0.20: "0.20"}.get(amount, f"{amount:g}")
        p["charges"]["amount"] = f(f"{amount_raw} mmol", {"number": amount, "unit": "mmol"}, ev)
    if mass is not None:
        p["charges"]["mass"] = f(f"{mass:g} mg", {"number": mass, "unit": "mg"}, ev)
    if volume is not None:
        volume_raw = f"{volume:.1f}" if float(volume).is_integer() else f"{volume:g}"
        p["charges"]["volume"] = f(f"{volume_raw} mL", {"number": volume, "unit": "mL"}, ev)
    if equivalents is not None:
        equiv_raw = f"{equivalents:.1f}" if float(equivalents).is_integer() else f"{equivalents:g}"
        p["charges"]["equivalents"] = f(f"{equiv_raw} equiv.", {"number": equivalents, "unit": "equiv"}, ev)
    if raw_equiv is not None:
        p["charges"]["equivalents"] = f(raw_equiv, None, ev, status="pending_review")
    if mol_percent is not None:
        p["charges"]["mol_percent"] = f(f"{mol_percent:g} mol%", {"number": mol_percent, "unit": "mol%"}, ev)
    if ref is not None:
        p["equivalence_reference"] = f(ref, ref, ev)
    r["participants"].append(p)
    return p


def reaction_step(r, sid, ev, *, temp=None, time=None, operation="stirred", pids=None, method="original_extraction"):
    s = st.new_step(sid, ev)
    s["participant_ids"] = list(pids or [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"])
    s["operation"] = f(operation, operation, ev, method=method)
    def quantity_raw(value):
        if not isinstance(value, dict):
            return value
        if value.get("raw_text"):
            return value["raw_text"]
        if "number" in value:
            return f"{value['number']:g} {value['unit']}"
        if "minimum" in value and "maximum" in value:
            return f"{value['minimum']:g}–{value['maximum']:g} {value['unit']}"
        return value
    if temp is not None:
        norm = temp if isinstance(temp, dict) else temp
        s["temperature"] = f(quantity_raw(temp), norm, ev, method=method)
    if time is not None:
        norm = time if isinstance(time, dict) else time
        s["time"] = f(quantity_raw(time), norm, ev, method=method)
    r["steps"].append(s)
    return s


def outcome(r, oid, kind, targets, raw, norm, ev, *, product_role=None, yield_type=None, basis=None, scope=None, method="original_extraction"):
    o = st.new_outcome(oid, kind, targets, ev)
    o["value"] = f(raw, norm, ev, method=method)
    if product_role is not None:
        o["product_role"] = f(product_role, product_role, ev)
    if yield_type is not None:
        o["yield_type"] = f(yield_type, yield_type, ev)
    if basis is not None:
        o["measurement_basis"] = f(basis, basis, ev)
    if scope is not None:
        o["scope"] = f(scope, scope, ev)
    r["outcomes"].append(o)
    return o


LABEL_NUM = {"C-H": 1, "C-C": 2, "C-O": 3, "O-H": 24}


def set_classification(r, formed, ev, *, broken=(), coverage="partial", rationale=None):
    changes = []
    labels = []
    support = {}
    for a, b in formed:
        cid = f"b{len(changes)+1}"
        changes.append({"change_id": cid, "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
                        "before_order": 0, "after_order": 1, "change_type": "formed"})
        code = "-".join((a, b)) if f"{a}-{b}" in LABEL_NUM else "-".join((b, a))
        support.setdefault(code, []).append(cid)
    for a, b in broken:
        cid = f"b{len(changes)+1}"
        changes.append({"change_id": cid, "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
                        "before_order": 1, "after_order": 0, "change_type": "broken"})
    for code, ids in support.items():
        labels.append({"category_number": LABEL_NUM[code], "code": code, "supporting_change_ids": ids})
    deriv = rationale or "Partial net-bond interpretation from the cited source drawing; atom mapping and atom sources are unresolved."
    r["classification"]["bond_changes"] = [derived(c, ev, deriv) for c in changes]
    r["classification"]["labels"] = derived(labels, ev, deriv)
    r["classification"]["coverage"] = inferred(coverage, ev, deriv)
    r["classification"]["limitations"] = [
        "No complete atom-mapped reaction SMILES was asserted.",
        "Hydrogen/atom sources, coproducts, and all bond-order changes remain unresolved unless explicitly listed.",
    ]


def set_no_reaction_classification(r, ev, rationale):
    r["classification"]["bond_changes"] = []
    r["classification"]["labels"] = inferred([], ev, rationale)
    r["classification"]["coverage"] = inferred("not_applicable", ev, rationale)
    r["classification"]["limitations"] = [rationale]


def add_annulation_class(r, ev):
    set_classification(r, [("C", "C"), ("C", "O")], ev,
                       rationale="The depicted [4+2] annulation closes the 4H-pyran through one new C-C and one new C-O connection; exact atoms and H transfers are unresolved.")


def selectivity(r, prefix, target, metric, raw, ev, basis, *, norm=None, order="unreported"):
    value = norm if norm is not None else raw
    return outcome(r, prefix, "selectivity", [target], raw, value, ev,
                   product_role="desired annulation product",
                   basis=basis,
                   scope={"metric": metric, "configuration_order": order, "major_configuration": "unassigned unless separately stated"})


# Page 2: this-paper generic strategy; bioactive Figure 1 is context, not experiment.
r = base_reaction("scheme1_generic", [EVI(2), EVT(2)], record_type="generic_scheme", categories=[],
                  location={"page": 2, "region": "Scheme 1"},
                  name="bifunctional phosphonium salt-catalyzed [4 + 2] strategy")
participant(r, "p_isoindigo", "m_isoindigo_1", "reactant", [EVI(2)])
participant(r, "p_allene", "m_allenoate_2", "reactant", [EVI(2)])
participant(r, "p_product", "m_product_3_generic", "product", [EVI(2)])
r["procedure"]["original_text"] = f("Scheme 1 generic transformation: unsaturated heterocycles (C4 synthon) + allenes (C2 synthon) under bifunctional phosphonium salt/PTC.",
                                      "Scheme 1 generic transformation", [EVI(2)])
add_annulation_class(r, [EVI(2)])


# Page 3, Table 1 optimization: one row per experiment.
t1 = [
    (1, "P1", 2, 1, 88, "97:3", 2, "r.t.", 20),
    (2, "P2", 2, 2, 91, "94:6", 18, "r.t.", 20),
    (3, "P3", 2, 1, 85, "94:6", 12, "r.t.", 20),
    (4, "P4", 2, 8, 87, "94:6", 45, "r.t.", 20),
    (5, "P5", 2, 1, 90, "96:4", 75, "r.t.", 20),
    (6, "P6", 2, 0.5, 93, "95:5", 88, "r.t.", 20),
    (7, "P7", 2, 1, 89, "96:4", 8, "r.t.", 20),
    (8, "P8", 2, 7, 87, "95:5", 26, "r.t.", 20),
    (9, "P9", 2, 20, 77, "96:4", 44, "r.t.", 20),
    (10, "P10", 2, 5, 84, "90:10", 69, "r.t.", 20),
    (11, "P6", 1, 1, 95, "95:5", 88, "r.t.", 20),
    (12, "P6", 4, 0.3, 95, "94:6", 90, "r.t.", 20),
    (13, "P6", 4, 2, 94, "94:6", 94, {"number": 0, "unit": "°C"}, 20),
    (14, "P6", 4, 6, 96, "96:4", 96, {"number": -20, "unit": "°C"}, 20),
    (15, "P6", 4, 8, 96, "95:5", 96, {"number": -20, "unit": "°C"}, 10),
]
for entry, cat, base_eq, hours, yld, dr, ee, temp, loading in t1:
    ev = [EVI(3), EVT(3)]
    r = base_reaction(f"table1_entry_{entry:02d}", ev, categories=["optimization"],
                      location={"page": 3, "table": 1, "entry": entry})
    participant(r, "p_1a", "m_1a", "reactant", ev, raw_equiv="0.05 equiv.")
    participant(r, "p_2a", "m_2a", "reactant", ev, raw_equiv="0.075 equiv.")
    participant(r, "p_cat", f"m_cat_{cat.lower().replace('-', '_')}", "catalyst", ev,
                mol_percent=loading, ref="m_1a")
    participant(r, "p_base", "m_cs2co3", "reagent", ev, raw_equiv=f"{base_eq} equiv.")
    participant(r, "p_xylene", "m_xylene", "solvent", ev, volume=1.0)
    participant(r, "p_3a", "m_3a", "product", ev)
    r["solvent_conditions"]["composition"] = f("xylene", ["xylene"], ev)
    r["solvent_conditions"]["total_volume"] = f("1.0 mL", {"number": 1.0, "unit": "mL"}, ev)
    temp_raw = "room temperature" if temp == "r.t." else f"{temp['number']:g} °C"
    temp_norm = "r.t." if temp == "r.t." else temp
    reaction_step(r, "reaction", ev, temp=temp_norm, time={"number": hours, "unit": "h"})
    r["steps"][-1]["temperature"]["raw_value"] = temp_raw
    r["procedure"]["original_text"] = f({"table_footnote": "1a (0.05 equiv.), 2a (0.075 equiv.), catalyst, Cs2CO3 (× equiv.) in xylene (1.0 mL)", "entry": entry},
                                          {"source": "Table 1 row and footnotes", "entry": entry}, ev)
    r["notes"].append(f("The printed footnote reports 0.05 equiv. and 0.075 equiv. for 1a/2a; the denominator is not stated, so these were not normalized or silently corrected to mmol.",
                        "Printed reactant units are unresolved and retained verbatim.", ev))
    r["review"]["unresolved"].append("Table 1 reactant quantities are printed as equiv. without a denominator; normalized amounts remain unresolved.")
    outcome(r, "yield", "yield", ["m_3a"], f"{yld}%", {"number": yld, "unit": "%"}, ev,
            product_role="desired product 3a", yield_type="isolated", basis="isolated yield")
    selectivity(r, "dr", "m_3a", "dr", dr, ev, "1H NMR analysis of crude product", order="major:minor; configurations unassigned")
    selectivity(r, "ee", "m_3a", "ee", f"{ee}% ee", ev, "HPLC analysis",
                norm={"number": ee, "unit": "%"})
    add_annulation_class(r, ev)


# Page 3, Table 2 generic and 12 product rows.
r = base_reaction("table2_generic", [EVI(3), EVT(3), EVT(5)], record_type="generic_scheme", categories=[],
                  location={"page": 3, "table": 2, "region": "generic reaction and footnotes"})
participant(r, "p_1", "m_isoindigo_1", "reactant", [EVI(3)], amount=0.10)
participant(r, "p_2", "m_allenoate_2", "reactant", [EVI(3)], amount=0.15)
participant(r, "p_p6", "m_cat_p6", "catalyst", [EVI(3)], amount=0.01, mol_percent=10, ref="m_isoindigo_1")
participant(r, "p_base", "m_cs2co3", "reagent", [EVI(3)], amount=0.40)
participant(r, "p_xylene", "m_xylene", "solvent", [EVI(3)], volume=2.0)
participant(r, "p_3", "m_product_3_generic", "product", [EVI(3)])
r["solvent_conditions"]["composition"] = f("xylene", ["xylene"], [EVI(3)])
r["solvent_conditions"]["total_volume"] = f("2.0 mL", {"number": 2.0, "unit": "mL"}, [EVI(3)])
reaction_step(r, "reaction", [EVI(3)], temp={"number": -20, "unit": "°C"}, time={"minimum": 4, "maximum": 12, "unit": "h", "raw_text": "4–12 h"})
r["procedure"]["original_text"] = f("General Procedure for Asymmetric Synthesis of 3 (page 5) plus Table 2 footnotes.",
                                      "General Procedure for Asymmetric Synthesis of 3", [EVI(3), EVT(5)])
add_annulation_class(r, [EVI(3)])

t2 = {
    "3a": (96, 8, None, "95:5", 96), "3b": (98, 8, None, ">20:1", 96),
    "3c": (90, 12, None, ">20:1", 84), "3d": (95, 6, None, "96:4", 90),
    "3e": (93, 6, None, "91:9", 88), "3f": (92, 12, None, ">20:1", 85),
    "3g": (None, None, None, None, None), "3h": (94, 5, "12.5:1", "93:7", 91),
    "3i": (75, 12, "4:1", ">20:1", 95), "3j": (96, 4, None, "91:9", 95),
    "3k": (95, 5, None, ">20:1", 96), "3l": (95, 5, None, ">20:1", 94),
}


def add_scope_workup(r, ev, product_mid):
    participant(r, "p_etoac_wash", "m_etoac", "workup_agent", ev, stage="workup", volume=2.0)
    participant(r, "p_silica", "m_silica_gel", "purification_agent", ev, stage="purification")
    participant(r, "p_petether", "m_petroleum_ether", "purification_agent", ev, stage="purification")
    participant(r, "p_etoac_eluent", "m_etoac", "purification_agent", ev, stage="purification")
    s = st.new_step("filter_wash_evaporate", ev, stage="workup")
    s["participant_ids"] = ["p_etoac_wash"]
    s["operation"] = f("filtered to remove residual base; solid washed with ethyl acetate (2.0 mL); evaporated under reduced pressure",
                       ["filter", "wash with EtOAc 2.0 mL", "evaporate under reduced pressure"], ev, method="inherited_general_condition")
    r["steps"].append(s)
    s = st.new_step("flash_chromatography", ev, stage="purification")
    s["participant_ids"] = ["p_silica", "p_petether", "p_etoac_eluent"]
    s["operation"] = f("flash column chromatography on silica gel (petroleum ether/ethyl acetate = 8:1 to 5:1)",
                       {"method": "flash silica chromatography", "eluent": "petroleum ether/ethyl acetate", "gradient": "8:1 to 5:1"},
                       ev, method="inherited_general_condition")
    r["steps"].append(s)


for label, (yld, hours, rr, dr, ee) in t2.items():
    ev = [EVI(3), EVT(3), EVT(5)] + ([EVT(4), EVI(4)] if label == "3h" else [])
    r = base_reaction(f"table2_{label}", ev, categories=["substrate_scope"],
                      location={"page": 3, "table": 2, "product_row": label})
    r["general_scheme_ids"] = [f"{DOC}:table2_generic"]
    participant(r, "p_1", "m_isoindigo_1", "reactant", ev, amount=0.10)
    participant(r, "p_2", "m_allenoate_2", "reactant", ev, amount=0.15)
    participant(r, "p_p6", "m_cat_p6", "catalyst", ev, amount=0.01, mol_percent=10, ref="m_isoindigo_1")
    participant(r, "p_base", "m_cs2co3", "reagent", ev, amount=0.40)
    participant(r, "p_xylene", "m_xylene", "solvent", ev, volume=2.0)
    if label != "3g":
        participant(r, "p_product", f"m_{label}", "product", ev)
    r["reaction_scale"] = f("1 (0.10 mmol)", {"reference_molecule_id": "m_isoindigo_1", "amount": {"number": 0.10, "unit": "mmol"}}, ev)
    r["solvent_conditions"]["composition"] = f("xylene", ["xylene"], ev)
    r["solvent_conditions"]["total_volume"] = f("2.0 mL", {"number": 2.0, "unit": "mL"}, ev)
    if hours is not None:
        reaction_step(r, "reaction", ev, temp={"number": -20, "unit": "°C"}, time={"number": hours, "unit": "h"}, method="inherited_general_condition")
    else:
        reaction_step(r, "reaction", ev, temp={"number": -20, "unit": "°C"}, time={"minimum": 4, "maximum": 12, "unit": "h", "raw_text": "4–12 h"}, method="inherited_general_condition")
    r["procedure"]["original_text"] = f("General Procedure for Asymmetric Synthesis of 3: mixture stirred at -20 °C to TLC completion; filtered, washed, evaporated, flash chromatographed.",
                                          "General Procedure for Asymmetric Synthesis of 3", [EVT(5)], method="inherited_general_condition")
    if label == "3g":
        outcome(r, "no_product", "no_reaction", [], "3g, n.r.; annulation product 3g was not observed", "not observed", ev,
                scope={"intended_product_label": "3g", "observation": "not observed"})
        set_no_reaction_classification(r, ev, "No annulation product was observed; no intended-product bond labels were assigned.")
        r["review"]["unresolved"].append("No detection limit or method for the no-reaction observation is reported.")
    else:
        add_scope_workup(r, [EVT(5)], f"m_{label}")
        outcome(r, "yield", "yield", [f"m_{label}"], f"{yld}%", {"number": yld, "unit": "%"}, ev,
                product_role=f"desired product {label}", yield_type="isolated", basis="isolated yield")
        if rr:
            selectivity(r, "rr", f"m_{label}", "rr", rr, ev, "1H NMR or HPLC analysis of crude product", order="major regioisomer:minor regioisomer; structures not separately assigned")
        selectivity(r, "dr", f"m_{label}", "dr", dr, ev, "1H NMR or HPLC analysis of crude product", order="major diastereomer:minor diastereomer; configurations not separately assigned")
        selectivity(r, "ee", f"m_{label}", "ee", f"{ee}% ee", ev, "HPLC analysis", norm={"number": ee, "unit": "%"})
        add_annulation_class(r, ev)
        if label == "3h":
            r["notes"].append(f("Absolute configuration of 3h was determined as (3'S, 4S) by X-ray analysis; CCDC 1937754.",
                                {"assignment": "(3'S, 4S)", "method": "X-ray", "CCDC": "1937754"}, ev))
    r["review"]["unresolved"].append("Specific compact substrate structures for this scope row are represented by the Table 2 family templates rather than guessed concrete SMILES.")


# Page 4, Table 3 generic and 11 scope rows.
r = base_reaction("table3_generic", [EVI(4), EVT(4), EVT(5)], record_type="generic_scheme", categories=[],
                  location={"page": 4, "table": 3, "region": "generic reaction and footnotes"})
participant(r, "p_4", "m_pyrazoloneyldiene_4", "reactant", [EVI(4)], amount=0.10)
participant(r, "p_2", "m_allenoate_2", "reactant", [EVI(4)], amount=0.15)
participant(r, "p_p10", "m_cat_p10", "catalyst", [EVI(4)], amount=0.01, mol_percent=10, ref="m_pyrazoloneyldiene_4")
participant(r, "p_base", "m_cs2co3", "reagent", [EVI(4)], amount=0.40)
participant(r, "p_xylene", "m_xylene", "solvent", [EVI(4)], volume=2.0)
participant(r, "p_5", "m_product_5_generic", "product", [EVI(4)])
r["solvent_conditions"]["composition"] = f("xylene", ["xylene"], [EVI(4)])
r["solvent_conditions"]["total_volume"] = f("2.0 mL", {"number": 2.0, "unit": "mL"}, [EVI(4)])
reaction_step(r, "reaction", [EVI(4)], temp={"number": 0, "unit": "°C"}, time={"minimum": 3.5, "maximum": 24, "unit": "h", "raw_text": "3.5–24 h"})
r["procedure"]["original_text"] = f("General Procedure for Asymmetric Synthesis of 5 (page 5) plus Table 3 footnotes.",
                                      "General Procedure for Asymmetric Synthesis of 5", [EVI(4), EVT(5)])
add_annulation_class(r, [EVI(4)])

t3 = {
    "5a": (96, 3.5, "19:1", "94:6", 92), "5b": (92, 5, ">20:1", "96:4", 90),
    "5c": (93, 6, "19:1", "94:6", 97), "5d": (92, 9, ">20:1", "97:3", 92),
    "5e": (90, 6, "16:1", "93:7", 93), "5f": (86, 6, ">20:1", "94:6", 87),
    "5g": (88, 6, "5.3:1", "94:6", 93), "5h": (82, 6, ">20:1", "97:3", 90),
    "5i": (83, 6, "7:1", "94:6", 82), "5j": (85, 24, ">20:1", "98:2", 90),
    "5k": (88, 9, ">20:1", "97:3", 88),
}
for label, (yld, hours, rr, dr, ee) in t3.items():
    ev = [EVI(4), EVT(4), EVT(5)]
    r = base_reaction(f"table3_{label}", ev, categories=["substrate_scope"],
                      location={"page": 4, "table": 3, "product_row": label})
    r["general_scheme_ids"] = [f"{DOC}:table3_generic"]
    participant(r, "p_4", "m_pyrazoloneyldiene_4", "reactant", ev, amount=0.10)
    participant(r, "p_2", "m_allenoate_2", "reactant", ev, amount=0.15)
    participant(r, "p_p10", "m_cat_p10", "catalyst", ev, amount=0.01, mol_percent=10, ref="m_pyrazoloneyldiene_4")
    participant(r, "p_base", "m_cs2co3", "reagent", ev, amount=0.40)
    participant(r, "p_xylene", "m_xylene", "solvent", ev, volume=2.0)
    participant(r, "p_product", f"m_{label}", "product", ev)
    r["reaction_scale"] = f("4 (0.10 mmol)", {"reference_molecule_id": "m_pyrazoloneyldiene_4", "amount": {"number": 0.10, "unit": "mmol"}}, ev)
    r["solvent_conditions"]["composition"] = f("xylene", ["xylene"], ev)
    r["solvent_conditions"]["total_volume"] = f("2.0 mL", {"number": 2.0, "unit": "mL"}, ev)
    reaction_step(r, "reaction", ev, temp={"number": 0, "unit": "°C"}, time={"number": hours, "unit": "h"}, method="inherited_general_condition")
    r["procedure"]["original_text"] = f("General Procedure for Asymmetric Synthesis of 5: mixture stirred at 0 °C to TLC completion; filtered, washed, evaporated, flash chromatographed.",
                                          "General Procedure for Asymmetric Synthesis of 5", [EVT(5)], method="inherited_general_condition")
    add_scope_workup(r, [EVT(5)], f"m_{label}")
    outcome(r, "yield", "yield", [f"m_{label}"], f"{yld}%", {"number": yld, "unit": "%"}, ev,
            product_role=f"desired product {label}", yield_type="isolated", basis="isolated yield")
    selectivity(r, "rr", f"m_{label}", "rr", rr, ev, "1H NMR or HPLC analysis of crude product", order="major regioisomer:minor regioisomer; structures not separately assigned")
    selectivity(r, "dr", f"m_{label}", "dr", dr, ev, "1H NMR or HPLC analysis of crude product", order="major diastereomer:minor diastereomer; configurations not separately assigned")
    selectivity(r, "ee", f"m_{label}", "ee", f"{ee}% ee", ev, "HPLC analysis", norm={"number": ee, "unit": "%"})
    add_annulation_class(r, ev)
    if label == "5c":
        r["notes"].append(f("Relative and absolute configuration of 5c was determined as (3'S, 4S) by X-ray analysis; CCDC 1968335.",
                            {"assignment": "(3'S, 4S)", "method": "X-ray", "CCDC": "1968335"}, ev))
    r["review"]["unresolved"].append("Specific compact substrate structures for this scope row are represented by the Table 3 family templates rather than guessed concrete SMILES.")


# Page 4 Scheme 2a: preparative-scale run.
ev = [EVI(4), EVT(4)]
r = base_reaction("scheme2a_scaleup_3h", ev, categories=["scale_up"], location={"page": 4, "scheme": "2a"})
participant(r, "p_1h", "m_1h", "reactant", ev, amount=1.0, mass=541)
participant(r, "p_2a", "m_2a", "reactant", ev, equivalents=1.5, ref="m_1h")
participant(r, "p_p6", "m_cat_p6", "catalyst", ev, mol_percent=10, ref="m_1h")
participant(r, "p_base", "m_cs2co3", "reagent", ev, equivalents=4.0, ref="m_1h")
participant(r, "p_xylene", "m_xylene", "solvent", ev)
participant(r, "p_3h", "m_3h", "product", ev)
r["reaction_scale"] = f("1h (541 mg, 1 mmol)", {"reference_molecule_id": "m_1h", "amount": {"number": 1.0, "unit": "mmol"}, "mass": {"number": 541, "unit": "mg"}}, ev)
r["solvent_conditions"]["composition"] = f("xylene", ["xylene"], ev)
reaction_step(r, "reaction", ev, temp={"number": -20, "unit": "°C"}, time={"number": 5, "unit": "h"})
outcome(r, "yield", "yield", ["m_3h"], "96% yield", {"number": 96, "unit": "%"}, ev, product_role="desired product 3h")
selectivity(r, "rr", "m_3h", "rr", "11.5:1", ev, None, order="major:minor; configurations unassigned")
selectivity(r, "dr", "m_3h", "dr", ">20:1", ev, None, order="major:minor; configurations unassigned")
selectivity(r, "ee", "m_3h", "ee", "90% ee", ev, None, norm={"number": 90, "unit": "%"})
add_annulation_class(r, ev)
r["review"]["unresolved"].append("Scheme 2a does not state solvent volume, workup, purification, or measurement methods for yield/rr/dr/ee.")


# Page 4 Scheme 2b: hydrolysis.
r = base_reaction("scheme2b_hydrolysis_6", ev, categories=["derivatization"],
                  location={"page": 4, "scheme": "2b"}, name="selective hydrolysis",
                  transformation="ester hydrolysis of 3h to compound 6")
participant(r, "p_3h", "m_3h", "reactant", ev)
participant(r, "p_fecl3", "m_fecl3", "reagent", ev, mol_percent=15, ref="m_3h")
participant(r, "p_dcm", "m_dcm", "solvent", ev)
participant(r, "p_6", "m_6", "product", ev)
r["solvent_conditions"]["composition"] = f("CH2Cl2", ["CH2Cl2"], ev)
reaction_step(r, "reaction", ev, temp="r.t.", time={"number": 10, "unit": "min"})
r["steps"][-1]["temperature"]["raw_value"] = "RT"
outcome(r, "yield", "yield", ["m_6"], "94% yield", {"number": 94, "unit": "%"}, ev,
        product_role="hydrolysed compound 6", yield_type="isolated", basis="prose states isolated yield")
selectivity(r, "dr", "m_6", "dr", ">20:1", ev, None, order="major:minor; configurations unassigned")
selectivity(r, "ee", "m_6", "ee", "92% ee", ev, None, norm={"number": 92, "unit": "%"})
set_classification(r, [("O", "H")], ev, broken=[("C", "O")],
                   rationale="The depicted tert-butyl ester-to-carboxylic-acid hydrolysis forms O-H and cleaves C-O; full atom mapping and coproducts are unresolved.")
r["review"]["unresolved"].append("Scheme 2b does not state substrate amount, solvent volume, workup, purification, or dr/ee measurement method.")


# Page 4 Scheme 2c: deuterium-labelling mechanistic experiment.
r = base_reaction("scheme2c_deuterium", ev, categories=["mechanistic_study"],
                  location={"page": 4, "scheme": "2c"}, name="deuterium-labelling experiment",
                  transformation="[4+2] annulation with deuterium incorporation")
participant(r, "p_1a", "m_1a", "reactant", ev)
participant(r, "p_2a", "m_2a", "reactant", ev)
participant(r, "p_p6", "m_cat_p6", "catalyst", ev, mol_percent=10, ref="m_1a")
participant(r, "p_base", "m_cs2co3", "reagent", ev, equivalents=4.0, ref="m_1a")
participant(r, "p_d2o", "m_d2o", "additive", ev, equivalents=10.0, ref="m_1a")
participant(r, "p_xylene", "m_xylene", "solvent", ev)
participant(r, "p_d3a", "m_d3a", "product", ev)
r["solvent_conditions"]["composition"] = f("xylene", ["xylene"], ev)
reaction_step(r, "reaction", ev, temp="r.t.", time={"number": 0.3, "unit": "h"})
r["steps"][-1]["temperature"]["raw_value"] = "RT"
outcome(r, "yield", "yield", ["m_d3a"], "91% yield", {"number": 91, "unit": "%"}, ev, product_role="deuterium-labelled d-3a")
outcome(r, "deuterium_incorporation", "other", ["m_d3a"], "30% D at 3'-position", {"number": 30, "unit": "%"}, ev,
        scope={"metric": "deuterium incorporation", "site": "3'-position"})
set_classification(r, [("C", "C"), ("C", "O"), ("C", "H")], ev,
                   rationale="Scheme 2c depicts the annulation C-C/C-O closures plus net C-D incorporation; D is classified under C-H taxonomy while isotope/site are retained in the outcome.")
r["review"]["unresolved"].append("Starting-material amounts, solvent volume, workup/purification, yield basis, and D-incorporation measurement method are unreported.")


# Page 5, Table 4 generic and five control rows.
r = base_reaction("table4_generic", [EVI(5), EVT(5)], record_type="generic_scheme", categories=[],
                  location={"page": 5, "table": 4, "region": "generic reaction above rows"})
participant(r, "p_1", "m_isoindigo_1", "reactant", [EVI(5)], amount=0.05)
participant(r, "p_2a", "m_2a", "reactant", [EVI(5)], amount=0.075)
participant(r, "p_p6", "m_cat_p6", "catalyst", [EVI(5)], amount=0.01, mol_percent=20, ref="m_isoindigo_1")
participant(r, "p_base", "m_cs2co3", "reagent", [EVI(5)], amount=0.20)
participant(r, "p_xylene", "m_xylene", "solvent", [EVI(5)], volume=1.0)
participant(r, "p_3", "m_product_3_generic", "product", [EVI(5)])
r["solvent_conditions"]["composition"] = f("xylene", ["xylene"], [EVI(5)])
r["solvent_conditions"]["total_volume"] = f("1.0 mL", {"number": 1.0, "unit": "mL"}, [EVI(5)])
reaction_step(r, "reaction", [EVI(5)], temp="r.t.", time=None)
r["steps"][-1]["temperature"]["raw_value"] = "room temperature"
add_annulation_class(r, [EVI(5)])

t4 = [
    (1, "P6", "Boc", 0.3, 95, "95:5", 90, "xylene"),
    (2, "P6-1", "Boc", 2, 90, "84:16", 25, "xylene"),
    (3, "P6-2", "Boc", 5, 84, "71:29", 86, "xylene"),
    (4, "P6", "Me", 5, 94, "90:10", 0, "xylene"),
    (5, "P6", "Boc", 1/60, None, None, None, "MeOH"),
]
for entry, cat, pg, hours, yld, dr, ee, solvent in t4:
    ev = [EVI(5), EVT(5)]
    r = base_reaction(f"table4_entry_{entry:02d}", ev, categories=["control", "mechanistic_study"],
                      location={"page": 5, "table": 4, "entry": entry})
    r["general_scheme_ids"] = [f"{DOC}:table4_generic"]
    substrate = "m_1_pg_me" if pg == "Me" else "m_1a"
    product = "m_3_pg_me" if pg == "Me" else "m_3a"
    participant(r, "p_1", substrate, "reactant", ev, amount=0.05)
    participant(r, "p_2a", "m_2a", "reactant", ev, amount=0.075)
    participant(r, "p_cat", f"m_cat_{cat.lower().replace('-', '_')}", "catalyst", ev, amount=0.01, mol_percent=20, ref=substrate)
    participant(r, "p_base", "m_cs2co3", "reagent", ev, amount=0.20)
    solvent_mid = "m_meoh" if solvent == "MeOH" else "m_xylene"
    participant(r, "p_solvent", solvent_mid, "solvent", ev, volume=1.0)
    if entry != 5:
        participant(r, "p_product", product, "product", ev)
    r["reaction_scale"] = f("1 (0.05 mmol)", {"reference_molecule_id": substrate, "amount": {"number": 0.05, "unit": "mmol"}}, ev)
    r["solvent_conditions"]["composition"] = f(solvent, [solvent], ev)
    r["solvent_conditions"]["total_volume"] = f("1.0 mL", {"number": 1.0, "unit": "mL"}, ev, method="inherited_general_condition")
    tval = {"number": 1, "unit": "min"} if entry == 5 else {"number": hours, "unit": "h"}
    reaction_step(r, "reaction", ev, temp="r.t.", time=tval)
    r["steps"][-1]["temperature"]["raw_value"] = "room temperature"
    if entry == 5:
        outcome(r, "messy", "decomposition", [], "messy reaction; no formation of any major products", "messy reaction; no major product", ev,
                scope={"observation": "no major product", "detection_limit": "unreported"})
        set_no_reaction_classification(r, ev, "No major product was observed in MeOH; intended-product bond labels were not assigned.")
        r["review"]["unresolved"].append("No individual decomposition products, detection method, or detection limit are reported.")
    else:
        outcome(r, "yield", "yield", [product], f"{yld}%", {"number": yld, "unit": "%"}, ev,
                product_role="corresponding annulation product", yield_type="isolated", basis="isolated yield")
        selectivity(r, "dr", product, "dr", dr, ev, "1H NMR analysis of crude product", order="major:minor; configurations unassigned")
        selectivity(r, "ee", product, "ee", f"{ee}% ee", ev, "HPLC analysis", norm={"number": ee, "unit": "%"})
        add_annulation_class(r, ev)
        if entry == 4:
            r["notes"].append(f("N-Methyl isoindigo gave a racemic version; authors infer N-Boc may be important for stereoselectivity.",
                                {"observation": "0% ee", "author_interpretation": "N-Boc may be crucial"}, ev))


# Page 5 TS1 is a proposed transition-state model, not a demonstrated experiment.
r = base_reaction("ts1_proposed_stereo_model", [EVI(5), EVT(5)], record_type="proposed_mechanism", categories=[],
                  location={"page": 5, "region": "TS1 (for alpha-addition)"},
                  name="plausible transition-state model TS1", transformation="proposed intermolecular alpha-addition")
participant(r, "p_1a", "m_1a", "reactant", [EVI(5)])
participant(r, "p_2a", "m_2a", "reactant", [EVI(5)])
participant(r, "p_p6", "m_cat_p6", "catalyst", [EVI(5)])
r["procedure"]["original_text"] = f("TS1 (for alpha-addition), Re-attack; weak-bond interactions including H-bonding and ion-pair manners proposed.",
                                      "TS1 proposed stereocontrol model", [EVI(5), EVT(5)])
set_classification(r, [("C", "C")], [EVI(5)],
                   rationale="TS1 depicts a proposed intermolecular alpha-addition C-C bond-forming step; it is a mechanistic proposal, not experimental proof.")
r["review"]["unresolved"] = ["TS1 is a proposed model; no elementary-step kinetics or complete atom mapping is supplied in the bounded pages."]


bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
bundle["source_context"] = [
    {
        "context_id": f"{DOC}:bibliography",
        "type": "bibliographic_metadata",
        "title": "Bifunctional Phosphonium Salt-catalyzed Enantioselective [4 + 2] Annulation of Isoindigos with Allenes: Access to Complex Heterocycles with Centerpiece of 4H-Pyrans",
        "authors": ["Juan Du", "Jia-Hong Wu", "Lixiang Zhu", "Xiaoyu Ren", "Chunhui Jiang", "Tianli Wang"],
        "doi": "10.1002/adsc.202000150",
        "journal": "Advanced Synthesis & Catalysis",
        "year": 2020,
        "volume": 362,
        "pages": "2510-2516",
        "evidence_ids": [EVI(1), EVT(1)],
    },
    {
        "context_id": f"{DOC}:figure1_bioactive_context",
        "type": "non_experimental_context",
        "description": "Figure 1 shows selected bioactive 4H-pyran scaffolds; these are context structures, not experiments reported in this paper.",
        "evidence_ids": [EVI(2)],
    },
    {
        "context_id": f"{DOC}:si_unassociated",
        "type": "source_limitation",
        "description": "The article states Supporting Information is available, but no SI file was supplied or explicitly associated in this workspace.",
        "evidence_ids": [EVI(1), EVT(1)],
    },
    {
        "context_id": f"{DOC}:page6_references",
        "type": "non_experimental_context",
        "description": "PDF page 6 contains acknowledgements continuation and references only; cited literature was not converted into current-paper experiments.",
        "evidence_ids": [EVI(6), EVT(6)],
    },
]

# Persist raw outputs and formal submission.
write(RUN / "raw_model/stereo_diagnostics_pages_0001_0006.json", stereo_diagnostics)
write(RUN / "raw_model/extracted_bundle_pages_0001_0006.json", bundle)
write(RUN / "submitted.json", bundle)

mentions = st.source_mentions(DATASET, ["1a", "2a", "3a", "3h", "5c", "compound 6"], context=220)
mentions["matches"] = [m for m in mentions["matches"] if m["page"] <= 6]
write(RUN / "raw_model/source_mentions_pages_0001_0006.json", mentions)

ids = [r["reaction_id"] for r in reactions]
by_prefix = lambda prefix: [rid for rid in ids if rid.startswith(f"{DOC}:{prefix}")]
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title, abstract, and both body-text columns", "finding": "Bibliographic metadata and study overview; no current-paper reaction arrow/table row on this page.", "evidence_ids": [EVI(1), EVT(1)]},
            {"region": "Supporting Information notice", "finding": "SI availability is advertised, but no SI was supplied or explicitly associated.", "evidence_ids": [EVI(1), EVT(1)]},
        ], "record_ids": [], "unresolved": ["No SI association is established."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Figure 1", "finding": "Bioactive 4H-pyran context only; not a current-paper synthetic experiment.", "evidence_ids": [EVI(2)]},
            {"region": "Scheme 1", "finding": "This-paper generic bifunctional-phosphonium [4+2] strategy; recorded as a generic scheme.", "evidence_ids": [EVI(2), EVT(2)]},
            {"region": "right-column optimization discussion", "finding": "Narrative descriptions of Table 1 entries; merged into the page-3 row records rather than duplicated.", "evidence_ids": [EVI(2), EVT(2)]},
        ], "record_ids": [f"{DOC}:scheme1_generic"], "unresolved": ["Scheme 1 variable attachment details remain template-level."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Table 1 reaction drawing, catalyst drawings P1-P10, rows 1-15, and footnotes a-f", "finding": "Fifteen optimization experiments recorded separately; printed reactant unit 'equiv.' retained as unresolved.", "evidence_ids": [EVI(3), EVT(3)]},
            {"region": "Table 2 generic arrow, product rows 3a-3l, X-ray inset, and footnotes a-c", "finding": "One generic scheme plus twelve scope rows, including 3g no reaction and 3h X-ray assignment.", "evidence_ids": [EVI(3), EVT(3)]},
            {"region": "lower prose", "finding": "Scope narrative cross-checked against product rows; not duplicated as new experiments.", "evidence_ids": [EVI(3), EVT(3)]},
        ], "record_ids": by_prefix("table1_") + by_prefix("table2_"), "unresolved": ["Most compact structures are not encoded as concrete SMILES.", "Table 1 0.05/0.075 'equiv.' denominator is not stated."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "Table 3 generic arrow, product rows 5a-5k, X-ray inset, and footnotes a-c", "finding": "One generic scheme and eleven scope experiments recorded; 5c X-ray assignment retained.", "evidence_ids": [EVI(4), EVT(4)]},
            {"region": "Scheme 2a", "finding": "Preparative 1 mmol-scale synthesis of 3h recorded separately from Table 2.", "evidence_ids": [EVI(4), EVT(4)]},
            {"region": "Scheme 2b", "finding": "FeCl3-mediated hydrolysis of 3h to 6 recorded as derivatization.", "evidence_ids": [EVI(4), EVT(4)]},
            {"region": "Scheme 2c", "finding": "D2O deuterium-labelling annulation recorded as a mechanistic-study experiment.", "evidence_ids": [EVI(4), EVT(4)]},
        ], "record_ids": by_prefix("table3_") + by_prefix("scheme2"), "unresolved": ["Scheme 2 omits several quantities, workups, and measurement methods."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "Table 4 generic drawing, rows 1-5, and footnotes a-e", "finding": "One generic control scheme and five distinct control experiments; MeOH row retained as messy/no-major-product without inventing a product.", "evidence_ids": [EVI(5), EVT(5)]},
            {"region": "TS1 drawing and surrounding prose", "finding": "Proposed stereocontrol transition-state model recorded as proposed mechanism, not experimental proof.", "evidence_ids": [EVI(5), EVT(5)]},
            {"region": "General Procedures for products 3 and 5", "finding": "Amounts, TLC endpoint, filtration/wash/evaporation, and chromatography gradient inherited only by the matching Table 2/3 scope records.", "evidence_ids": [EVI(5), EVT(5)]},
        ], "record_ids": by_prefix("table4_") + [f"{DOC}:ts1_proposed_stereo_model"], "unresolved": ["Compound-specific SI procedures and characterization are unavailable."]},
        {"page": 6, "state": "attempted", "regions": [
            {"region": "acknowledgements continuation and references [1]-[15]", "finding": "No current-paper experiment, reaction arrow, table row, or standalone spectrum; citations were not treated as new experiments.", "evidence_ids": [EVI(6), EVT(6)]},
        ], "record_ids": [], "unresolved": []},
    ],
}
write(RUN / "coverage/pages_0001_0006.json", coverage)

workflow = f"""# Workflow log — allene-stereo-v31-200-20261017

## Scope and source controls

- Document: `{DOC}`; DOI `10.1002/adsc.202000150`.
- True title: *Bifunctional Phosphonium Salt-catalyzed Enantioselective [4 + 2] Annulation of Isoindigos with Allenes: Access to Complex Heterocycles with Centerpiece of 4H-Pyrans*.
- Bounded pass: original PDF pages 1–6. Every supplied JPEG `page1.jpg` through `page6.jpg` was opened at original detail; page-3 tables were additionally cropped for legibility with evidence receipts under `runs/current/raw_model/`.
- The source mapping identifies a journal PDF (`presentation_conversion: false`), not a PPTX-derived presentation.
- The prepared dataset began with zero candidates and zero experiments. The pre-existing workflow described a different DOI and was replaced; none of its experimental claims were used.
- The article advertises SI, but no SI file was supplied or explicitly associated. Detailed SI-only methods remain missing.

## Page inventory and extraction decisions

- Page 1: title/DOI/authors/abstract and SI notice; no reaction arrow or table.
- Page 2: Figure 1 is biological context, not current-paper experiments. Scheme 1 is one generic transformation. Optimization prose is duplicate narrative for Table 1 rows on page 3.
- Page 3: Table 1 yielded 15 separate optimization records. Table 2 yielded one generic scheme plus 12 separate rows (3a–3l), including the genuine no-reaction 3g row. The 3h X-ray assignment and CCDC 1937754 are retained.
- Page 4: Table 3 yielded one generic scheme plus 11 separate rows (5a–5k). The 5c X-ray assignment and CCDC 1968335 are retained. Scheme 2a scale-up, 2b hydrolysis, and 2c D2O labelling are three distinct experiments.
- Page 5: Table 4 yielded one generic control scheme and five separate rows. Entry 5 is retained as a messy/no-major-product experiment without a fabricated product. TS1 is a proposed-mechanism record, not experimental proof. The two General Procedures were inherited only by their matching Table 2/3 scope records.
- Page 6: acknowledgements continuation and references only; cited literature was not turned into this paper's experiments.
- No independent mass-spectrum peak list, molecular-weight table, or fragmentation-abundance observation occurs on pages 1–6; no characterization-only synthetic records were created.

## Records, quantities, and outcomes

- Raw bundle contains {len(reactions)} reactions: 4 generic schemes, 15 optimization rows, 23 substrate-scope rows, 8 other specific experiments/controls, and 1 proposed mechanism.
- All yields, dr, rr, ee, deuterium incorporation, and no-reaction observations target their actual measured molecule when applicable. No ee was converted to er; configuration order remains unassigned unless the source explicitly assigns it.
- Table 1 footnote visibly prints `1a (0.05 equiv.)` and `2a (0.075 equiv.)`, but gives no denominator. These raw values were retained without normalizing or silently changing them to mmol. Cs2CO3 row equivalents are likewise retained without an invented denominator; this ambiguity is explicit in each record.
- Table 2 and Table 3 scope records inherit the printed general doses and page-5 operational sequence: TLC endpoint, filtration, EtOAc 2.0 mL wash, reduced-pressure evaporation, and silica flash chromatography with petroleum ether/EtOAc 8:1→5:1.
- Scheme 2a and 2c do not identify yield measurement basis; Scheme 2a/2b do not identify rr/dr/ee methods. These fields remain missing. Table/scope basis facts are populated only where footnotes or prose state them.
- No numerical source conflict was found within pages 1–6. The Table 1 printed unit anomaly is an unresolved source expression, not silently treated as a conflict with an unavailable SI.

## Structure, stereochemistry, and classification

- Every molecule has explicit chirality, geometry, and representation statuses. Missing stereochemical marks are never interpreted as racemic or achiral.
- Products 3h and 5c preserve the printed `(3'S, 4S)` X-ray assignments. Other product drawings preserve partial relative depiction without invented absolute configurations. N-methyl control product retains the source statement `0% ee`/racemic version.
- Allenoate templates and 2a retain an unassigned allene axis; no axial descriptor was invented. Compact substrates/products/catalysts remain unresolved or template representations rather than guessed SMILES.
- Candidate SMILES are limited to simple named reagents. Each was run through `standard stereo`; diagnostics are saved separately and remain machine-derived, not source verification.
- Annulations are partially classified for the source-supported new C–C and C–O ring connections. Hydrolysis records O–H formation and C–O cleavage. The D2O experiment also retains C–H taxonomy for C–D incorporation. TS1 has only proposed C–C formation. No complete atom mapping or three-model vote was fabricated.

## Files, validation, and audit

- Init template: `runs/current/raw_model/init_template_pages_0001_0006.json`.
- Reproducible builder: `runs/current/raw_model/build_pages_0001_0006.py`.
- Raw bundle: `runs/current/raw_model/extracted_bundle_pages_0001_0006.json`.
- Formal submission: `runs/current/submitted.json`.
- Source-mention hints and stereo diagnostics are stored separately under `runs/current/raw_model/`.
- Coverage: `runs/current/coverage/pages_0001_0006.json`, exactly pages 1–6, with nonempty inspected regions on every page and `semantic_accuracy_verified: false`.
- `standard merge` imported all 51 records and generated `standard_records.json`, `molecules.json`, `records/*.json`, `records/*.md`, `classification_index.json`, `规范记录预览.md`, the legacy compatibility view, and SQLite output.
- `standard validate` and legacy `validate --strict` both pass with zero errors. Snapshots are saved as `raw_model/standard_validation_pages_0001_0006.json` and `raw_model/strict_validation_pages_0001_0006.json`.
- The audit snapshot is `raw_model/audit_pages_0001_0006.json`. There are no `reaction_identity`, `outcome_target`, `missing_steps`, or characterization-only warnings.
- Audit retains 46 `structure_gap` advisories because the compact complex reactant/product structures were deliberately not converted into guessed SMILES. Two of those records (Scheme 2a and 2c) also retain `measurement_basis` advisories because the source does not state the yield method. These are genuine unresolved source gaps, not repaired by assumption.
- Generated Markdown was checked for the complete Table 2/3 general doses and operations (including EtOAc 2.0 mL and the 8:1→5:1 eluent gradient), the Scheme 2a 541 mg/1 mmol/1.5 equiv/10 mol%/4.0 equiv/-20 °C/5 h details, and the Table 1 printed-unit note.

## Unresolved matters

- All records remain `needs_review`; schema validity will not certify chemistry or complete recall.
- SI is unassociated, so compound-specific procedures, full characterization, and many measurement methods remain unavailable.
- Most compact structures and all complete atom mappings remain unresolved. Scope-row reactants are linked to explicit Table 2/3 family templates rather than guessed concrete SMILES.
- Page 7 was not part of the bounded visual pass and is not claimed as covered.
"""
(ROOT / "runs/workflow.md").write_text(workflow)

print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
                  "source_context": len(bundle["source_context"]), "coverage_pages": len(coverage["pages"])}, ensure_ascii=False))
