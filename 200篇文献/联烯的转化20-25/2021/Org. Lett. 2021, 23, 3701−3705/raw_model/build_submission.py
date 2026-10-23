#!/usr/bin/env python3
"""Build the source-derived v3.1 submission for DOI 10.1021/acs.orglett.1c01079.

This is a session-local transcription program.  It does not alter schemas or extraction
tools.  Dense structures remain unresolved when the five-page main text does not permit
a reliable stereochemically complete SMILES transcription.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills" / "chemical-reaction-extraction" / "scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st


DOC = "doc_da7f8b58e55ca30d"
EV = {p: f"{DOC}:page{p}:image" for p in range(1, 6)}
TEV = {p: f"{DOC}:page{p}:text" for p in range(1, 6)}
OUT = ROOT / "runs" / "current"


def f(raw=None, norm=None, ev=(), method="original_extraction", status=None, derivation=None):
    return st.fact(raw, norm, ev, method=method, status=status, derivation=derivation)


def missing(ev=(), status="unreported"):
    return st.fact(evidence=ev, status=status)


def set_fact(container, key, raw, norm, ev, method="original_extraction", derivation=None):
    container[key] = f(raw, norm, ev, method=method, derivation=derivation)


def unresolved_molecule(mid, label, ev, *, name=None, stereo="unresolved", assignment=None):
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    m["names"]["original"] = f(name or label, name or label, ev)
    m["identity_note"] = f(
        "Structure is source-identified by the cited drawing/label, but no concrete SMILES is asserted.",
        "Dense or variable source structure retained without guessed complete SMILES.", ev,
    )
    m["stereochemistry"]["chirality_status"] = f(stereo, stereo, ev)
    m["stereochemistry"]["geometry_status"] = f("not assigned", "unresolved", ev)
    m["stereochemistry"]["representation_status"] = f("not encoded", "not_encoded", ev)
    if assignment is not None:
        m["stereochemistry"]["assignment"] = f(assignment, assignment, ev)
    return m


def generic_molecule(mid, label, ev, fixed, variables, *, stereo="unresolved", reason=None):
    m = unresolved_molecule(mid, label, ev, stereo=stereo)
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": f(fixed, fixed, ev),
        "variables": [
            {
                "label": v,
                "attachment_points": f("shown attachment point in source drawing", [v], ev),
                "definition": f("source-defined variable", "source-defined variable", ev),
                "allowed_values": f("values shown across the cited scope", "scope-enumerated values", ev),
                "constraints": missing(ev),
            }
            for v in variables
        ],
        "instance_mapping": [],
        "reason": f(reason or "Generic R-group drawing in the source.", reason or "Generic R-group drawing in the source.", ev),
    }
    if stereo == "partially_specified":
        m["stereochemistry"]["assignment"] = f(
            "relative bridge stereochemistry is depicted in the generic product drawing; absolute configuration is not assigned",
            {"type": "relative", "absolute_configuration": "unassigned"}, ev)
    return m


stereo_diagnostics = {}


def smiles_molecule(mid, label, smiles, ev, *, name=None):
    """Register only small, unambiguous named reagents/solvents as candidate SMILES."""
    m = unresolved_molecule(mid, label, ev, name=name or label, stereo="achiral")
    diag = st.stereo_inspect(smiles)
    stereo_diagnostics[mid] = {"input_smiles": smiles, **diag}
    m["representation_type"] = "smiles"
    m["smiles"] = f(label, smiles, ev, method="structure_derivation",
                     derivation="Candidate structure expanded from the unambiguous reagent/solvent name or abbreviation in the source; parsing is not source verification.")
    m["stereochemistry"]["chirality_status"] = f(None, "achiral", ev, method="structure_derivation",
        derivation="No potential tetrahedral or axial stereogenic element in the registered small-molecule graph.")
    m["stereochemistry"]["geometry_status"] = f(None, "not_applicable", ev, method="structure_derivation",
        derivation="The registered graph has no independently configurable E/Z alkene.")
    m["stereochemistry"]["representation_status"] = f(None, "encoded", ev, method="structure_derivation",
        derivation="Canonical isomeric graph was computed from the candidate SMILES.")
    m["stereochemistry"]["isomeric_smiles"] = f(None, diag["canonical_isomeric_smiles"], ev, method="structure_derivation",
        derivation="RDKit canonical isomeric SMILES from the registered candidate graph.")
    m["stereochemistry"]["identity_key"] = f(None, diag["stereo_graph_key"], ev, method="structure_derivation",
        derivation="RDKit stereo graph key; not a sample identity or source verification.")
    return m


def base_reaction(rid, record_type, ev, location, categories, name, transformation):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, ev)
    r["location"] = f(location, location, ev)
    r["reaction_name"] = f(name, name, ev)
    r["transformation_type"] = f(None, transformation, ev, method="inference",
        derivation="Normalized from the source reaction title/arrow and product drawing.")
    r["review"]["unresolved"] = []
    return r


def participant(pid, mid, role, ev, *, stage="reaction", eq=None, amount=None, mass=None,
                volume=None, molpct=None, ref=None, product_role=None, charge_method="original_extraction"):
    p = st.new_participant(pid, mid, role, ev, stage=stage)
    if product_role:
        p["product_role"] = f(product_role, product_role, ev)
    for key, value in (("equivalents", eq), ("amount", amount), ("mass", mass),
                       ("volume", volume), ("mol_percent", molpct)):
        if value is not None:
            p["charges"][key] = f(value.get("raw_text", value), value, ev, method=charge_method)
    if ref is not None:
        p["equivalence_reference"] = f("relative to " + ref, ref, ev, method=charge_method)
    return p


def step(step_id, ev, pids, *, temp=None, time=None, atmosphere=None, operation=None, method="original_extraction"):
    s = st.new_step(step_id, ev)
    s["participant_ids"] = list(pids)
    if operation is not None:
        s["operation"] = f(operation, operation, ev, method=method)
    if temp is not None:
        s["temperature"] = f(temp.get("raw_text", temp) if isinstance(temp, dict) else temp, temp, ev, method=method)
    if time is not None:
        s["time"] = f(time.get("raw_text", time) if isinstance(time, dict) else time, time, ev, method=method)
    if atmosphere is not None:
        s["atmosphere"] = f(atmosphere, atmosphere, ev, method=method)
    return s


def outcome(oid, kind, targets, value, ev, *, yield_type=None, basis=None, scope=None):
    o = st.new_outcome(oid, kind, targets, ev)
    if value is not None:
        o["value"] = f(value.get("raw_text", value) if isinstance(value, dict) else value, value, ev)
    if yield_type is not None:
        o["yield_type"] = f(yield_type, yield_type, ev)
    if basis is not None:
        o["measurement_basis"] = f(basis, basis, ev)
    if scope is not None:
        o["scope"] = f(scope, scope, ev)
    return o


TAX = {"C-H": (1, "C", "H"), "C-C": (2, "C", "C"), "C-O": (3, "C", "O"),
       "C-N": (5, "C", "N"), "DB-C-C": (13, "C", "C"), "DB-C-O": (14, "C", "O"),
       "O-H": (24, "O", "H"), "N-H": (29, "N", "H")}


def classification(r, ev, codes, limitation):
    changes = []
    labels = []
    for idx, code in enumerate(codes, 1):
        number, a, b = TAX[code]
        change = {"change_id": f"b{idx}", "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
                  "before_order": 0 if not code.startswith("DB-") else 1,
                  "after_order": 2 if code.startswith("DB-") else 1,
                  "change_type": "order_increased" if code.startswith("DB-") else ("hydrogen_gain" if b == "H" else "formed")}
        changes.append(f(None, change, ev, method="inference", derivation=limitation))
        labels.append({"category_number": number, "code": code, "supporting_change_ids": [f"b{idx}"]})
    r["classification"]["bond_changes"] = changes
    r["classification"]["labels"] = f(None, labels, ev, method="inference", derivation=limitation)
    r["classification"]["coverage"] = f(None, "partial", ev, method="inference", derivation=limitation)
    r["classification"]["limitations"] = [limitation, "No atom-mapped reaction SMILES is asserted for the dense polycyclic source drawings."]


def no_classification(r, ev, reason):
    r["classification"]["labels"] = f(None, [], ev, method="inference", derivation=reason)
    r["classification"]["coverage"] = f(None, "not_applicable", ev, method="inference", derivation=reason)
    r["classification"]["limitations"] = [reason]


def set_solvent(r, ev, composition=None, volume=None, concentration=None, mixture=None, method="original_extraction"):
    sc = r["solvent_conditions"]
    if composition is not None:
        sc["composition"] = f(composition, composition, ev, method=method)
    if volume is not None:
        sc["total_volume"] = f(volume.get("raw_text", volume), volume, ev, method=method)
    if concentration is not None:
        sc["concentration"] = f(concentration.get("raw_text", concentration), concentration, ev, method=method)
        sc["concentration_basis"] = f("source does not identify a separate concentration denominator", "unspecified source basis", ev, method=method)
    if mixture is not None:
        sc["mixture_ratio"] = f(mixture, mixture, ev, method=method)


def pct(n): return {"number": n, "unit": "%", "raw_text": f"{n}%"}
def mmol(n): return {"number": n, "unit": "mmol", "raw_text": f"{n} mmol"}
def gram(n): return {"number": n, "unit": "g", "raw_text": f"{n} g"}
def mg(n): return {"number": n, "unit": "mg", "raw_text": f"{n} mg"}
def ml(n): return {"number": n, "unit": "mL", "raw_text": f"{n} mL"}
def equiv(n): return {"number": n, "unit": "equiv", "raw_text": f"{n} equiv"}
def molpct(n): return {"number": n, "unit": "mol%", "raw_text": f"{n} mol%"}
def deg(n): return {"number": n, "unit": "°C", "raw_text": f"{n} °C"}
def hours(a, b=None):
    return ({"number": a, "unit": "h", "raw_text": f"{a} h"} if b is None else
            {"minimum": a, "maximum": b, "unit": "h", "raw_text": f"{a}–{b} h"})


dataset = json.loads((OUT / "dataset.json").read_text())
bundle = st.empty(dataset["documents"])
molecules = []

# Generic families and central labelled compounds.
molecules += [
    generic_molecule(f"{DOC}:mol:series1", "1", [EV[2], EV[3]],
        "cyclopropyl-tethered allenyl ketone scaffold shown in Schemes 3 and 5", ["R1", "R2", "n"], stereo="unspecified",
        reason="Scope uses a generic allenyl-ketone drawing; individual input labels are not printed for every row."),
    generic_molecule(f"{DOC}:mol:series2", "2", [EV[2]], "benzo[c]isoxazole scaffold", ["R3"],
        reason="Scheme 3 varies benzo[c]isoxazole substitution."),
    generic_molecule(f"{DOC}:mol:series3", "3", [EV[3]], "N-iminoquinazolinium ylide scaffold", ["R3", "R4"],
        reason="Scheme 5 varies quinazolinium and sulfonyl substituents."),
    generic_molecule(f"{DOC}:mol:series4", "4", [EV[2]], "cyclobutane-bridged [4+2] cycloadduct scaffold", ["R1", "R2", "R3", "n"], stereo="partially_specified",
        reason="Scheme 3 gives a generic product with relative stereochemical drawing but no complete generic stereochemical SMILES."),
    generic_molecule(f"{DOC}:mol:series5", "5", [EV[3]], "cyclobutane-bridged [3+2] cycloadduct scaffold", ["R1", "R2", "R3", "R4", "n"], stereo="partially_specified",
        reason="Scheme 5 gives a generic product with relative stereochemical drawing but no complete generic stereochemical SMILES."),
]

for label, pages, name, stereo in [
    ("1a", [2, 3, 4], "cyclopropyl-tethered allenyl ketone 1a", "unspecified"),
    ("2a", [2], "benzo[c]isoxazole 2a", "unresolved"),
    ("3a", [3], "quinazolin-3-ium-3-yl(tosyl)amide 3a", "unresolved"),
]:
    ev = [EV[p] for p in pages]
    m = unresolved_molecule(f"{DOC}:mol:{label}", label, ev, name=name, stereo=stereo)
    if label == "1a":
        m["stereochemistry"]["axial"] = f("allene axis depicted; no R_a/S_a assignment reported", {"axis": "allene", "configuration": "unspecified"}, ev)
    molecules.append(m)

relative_4 = "Relative bridged-ring stereochemistry is depicted; no absolute configuration is assigned in the main text."
relative_5 = "Major diastereomer is depicted with wedge/dash bonds; no absolute configuration is assigned in the main text."
for label in [f"4{chr(97+i)}" for i in range(26)] + ["4aa"]:
    molecules.append(unresolved_molecule(f"{DOC}:mol:{label}", label, [EV[2]], stereo="partially_specified", assignment=relative_4))
for label in [f"5{chr(97+i)}" for i in range(26)] + [f"5a{chr(97+i)}" for i in range(6)]:
    molecules.append(unresolved_molecule(f"{DOC}:mol:{label}", label, [EV[3]], stereo="partially_specified", assignment=relative_5))

for label, page in [(label, 3) for label in ["6", "8", "9"]] + [(str(i), 4) for i in range(10, 22)]:
    molecules.append(unresolved_molecule(f"{DOC}:mol:{label}", label, [EV[page]], stereo="unresolved"))
for label in ["7a", "7b"]:
    molecules.append(unresolved_molecule(f"{DOC}:mol:{label}", label, [EV[3]], stereo="unresolved"))

# Named catalysts/additives and small solvent/reagent candidates.
unresolved_reagents = {
    "tBuXPhosAuCl": ("t-BuXPhosAuCl", [EV[2], EV[4]]), "PPh3AuOTf": ("PPh3AuOTf", [EV[2], EV[3], EV[4]]),
    "PPh3AuCl": ("PPh3AuCl", [EV[3], EV[4]]), "AgSbF6": ("AgSbF6", [EV[2], EV[3], EV[4]]),
    "AgTFA": ("AgTFA", [EV[4]]), "4AMS": ("4 Å molecular sieves", [EV[2], EV[3], EV[4]]),
    "IBX": ("IBX", [EV[3]]), "RaneyNi": ("Raney Ni", [EV[3]]),
}
for key, (name, ev) in unresolved_reagents.items():
    molecules.append(unresolved_molecule(f"{DOC}:mol:{key}", name, ev, name=name, stereo="unresolved"))

smiles_defs = {
    "DCE": ("DCE", "ClCCCl", [EV[2], EV[3], EV[4]], "1,2-dichloroethane"),
    "DCM": ("DCM", "ClCCl", [EV[2], EV[3]], "dichloromethane"),
    "EA": ("EA", "CCOC(C)=O", [EV[3]], "ethyl acetate"),
    "THF": ("THF", "C1CCOC1", [EV[3], EV[4]], "tetrahydrofuran"),
    "dioxane": ("1,4-dioxane", "O1CCOCC1", [EV[3]], "1,4-dioxane"),
    "Zn": ("Zn", "[Zn]", [EV[3]], "zinc"),
    "HOAc": ("HOAc", "CC(=O)O", [EV[3]], "acetic acid"),
    "propargyl_bromide": ("3-bromopropyne", "C#CCBr", [EV[3]], "3-bromopropyne"),
    "K2CO3": ("K2CO3", "[K+].[K+].[O-]C([O-])=O", [EV[3]], "potassium carbonate"),
    "MeCN": ("MeCN", "CC#N", [EV[3]], "acetonitrile"),
    "DMSO": ("DMSO", "CS(C)=O", [EV[3]], "dimethyl sulfoxide"),
    "TfOH": ("TfOH", "OS(=O)(=O)C(F)(F)F", [EV[3]], "trifluoromethanesulfonic acid"),
    "H2": ("H2", "[H][H]", [EV[3]], "hydrogen"),
    "MeOH": ("MeOH", "CO", [EV[3]], "methanol"),
}
for key, (label, smi, ev, name) in smiles_defs.items():
    molecules.append(smiles_molecule(f"{DOC}:mol:{key}", label, smi, ev, name=name))

# Generic literature/background molecules used only to classify Scheme 1 arrows.
for key, label in [("prior_diene", "generic diene"), ("prior_cyclobutenone", "cyclobutenone"),
                   ("prior_da_product", "bicyclo[4.2.0]octene product"), ("prior_dipole", "generic nitrone/dipole"),
                   ("prior_heterocyclobutene", "azetine or thiete"), ("prior_32_product", "fused isoxazoline product"),
                   ("insitu_intermediate", "in situ furan-fused cyclobutene")]:
    molecules.append(unresolved_molecule(f"{DOC}:mol:{key}", label, [EV[2]], stereo="unresolved"))

bundle["molecules"] = molecules
reactions = []

# General schemes for the two present-work platforms.
for rid, page, reactants, product, name, codes in [
    ("general_4plus2", 2, ["series1", "series2"], "series4", "[4 + 2] cycloaddition", ["C-C", "C-O"]),
    ("general_3plus2", 3, ["series1", "series3"], "series5", "[3 + 2] cycloaddition", ["C-C", "C-O", "C-N"]),
]:
    ev = [EV[page]]
    r = base_reaction(rid, "generic_scheme", ev, f"PDF page {page}, Scheme {'3' if page == 2 else '5'} general arrow", [], name, "gold-catalyzed tandem cycloisomerization/cycloaddition")
    for i, mid in enumerate(reactants): r["participants"].append(participant(f"p{i+1}", f"{DOC}:mol:{mid}", "reactant", ev))
    r["participants"].append(participant("product", f"{DOC}:mol:{product}", "product", ev, product_role="generic product"))
    classification(r, ev, codes, "Category-level partial inference from the general source arrow; individual atom identities are unresolved.")
    r["review"]["unresolved"] = ["Generic R-group structures are not atom mapped."]
    if rid == "general_3plus2":
        r["notes"] = [f("Scheme 5 title literally says '[4 + 2] Cycloaddition', while Scheme 4, the prose, ylide reactant, and product numbering describe the current transformation as [3 + 2]. Both source statements are preserved.",
                         {"scheme5_title":"[4 + 2] Cycloaddition","contextual_assignment":"[3 + 2] cycloaddition","status":"conflict"}, ev)]
        r["review"]["unresolved"].append("Source conflict: Scheme 5 title says [4 + 2], whereas Scheme 4/prose/body identify [3 + 2].")
    reactions.append(r)

# Scheme 1 prior-work arrows, explicitly kept as literature examples.
prior_specs = [
    ("prior_scheme1_eq1", ["prior_diene", "prior_cyclobutenone"], "prior_da_product", "Danishefsky prior-work [4 + 2]", ["C-C"]),
    ("prior_scheme1_eq2", ["prior_dipole", "prior_heterocyclobutene"], "prior_32_product", "Didier prior-work [3 + 2]", ["C-C", "C-O"]),
]
for rid, ins, out, name, codes in prior_specs:
    ev = [EV[2]]
    r = base_reaction(rid, "literature_example", ev, "PDF page 2, Scheme 1B", [], name, "cited prior-work cycloaddition")
    for i, mid in enumerate(ins): r["participants"].append(participant(f"p{i+1}", f"{DOC}:mol:{mid}", "reactant", ev))
    r["participants"].append(participant("product", f"{DOC}:mol:{out}", "product", ev))
    classification(r, ev, codes, "Only the generic literature-example drawing is available; classification is partial and not a record of this paper's experiment.")
    r["review"]["unresolved"] = ["Cited prior work only; experimental conditions and exact structures are outside the supplied pages."]
    reactions.append(r)

# Scheme 1C author working hypothesis/in-situ intermediate, not an observed elementary mechanism.
r = base_reaction("proposed_insitu_platform", "proposed_mechanism", [EV[2]], "PDF page 2, Scheme 1C", [],
                  "in situ furan-fused cyclobutene platform", "author working hypothesis")
r["participants"] = [participant("precursor", f"{DOC}:mol:series1", "reactant", [EV[2]]),
                     participant("intermediate", f"{DOC}:mol:insitu_intermediate", "product", [EV[2]])]
no_classification(r, [EV[2]], "Author working hypothesis is not treated as an observed net experiment.")
r["review"]["unresolved"] = ["Intermediate is proposed/in situ and not independently isolated in the supplied main text."]
reactions.append(r)


def standard_scope4(label, yield_n, substituent, reverse=False):
    ev = [EV[2]]
    r = base_reaction(f"scope_{label}", "specific_experiment", ev, f"PDF page 2, Scheme 3, entry {label}",
                      ["substrate_scope"], "[4 + 2] cycloaddition", "gold-catalyzed tandem cycloisomerization/[4+2] cycloaddition")
    r["general_scheme_ids"] = [f"{DOC}:general_4plus2"]
    ref = f"{DOC}:mol:series2" if reverse else f"{DOC}:mol:series1"
    p1_amt, p1_eq = (mmol(0.6), equiv(3.0)) if reverse else (mmol(0.2), equiv(1.0))
    p2_amt, p2_eq = (mmol(0.2), equiv(1.0)) if reverse else (mmol(0.6), equiv(3.0))
    r["participants"] = [
        participant("allene", f"{DOC}:mol:series1", "reactant", ev, amount=p1_amt, eq=p1_eq, ref=ref, charge_method="inherited_general_condition"),
        participant("isoxazole", f"{DOC}:mol:series2", "reactant", ev, amount=p2_amt, eq=p2_eq, ref=ref, charge_method="inherited_general_condition"),
        participant("product", f"{DOC}:mol:{label}", "product", ev, product_role="major cycloadduct"),
        participant("gold", f"{DOC}:mol:tBuXPhosAuCl", "catalyst", ev, molpct=molpct(2), ref=ref, charge_method="inherited_general_condition"),
        participant("silver", f"{DOC}:mol:AgSbF6", "additive", ev, molpct=molpct(15), ref=ref, charge_method="inherited_general_condition"),
        participant("sieves", f"{DOC}:mol:4AMS", "additive", ev, mass=mg(100), charge_method="inherited_general_condition"),
        participant("solvent", f"{DOC}:mol:DCE", "solvent", ev, volume=ml(2.0), charge_method="inherited_general_condition"),
    ]
    r["reaction_scale"] = f("0.2 mmol limiting component", {"reference_molecule_id": ref, "amount": mmol(0.2)}, ev)
    set_solvent(r, ev, "DCE", ml(2.0), method="inherited_general_condition")
    r["steps"] = [step("reaction", ev, [p["participant_id"] for p in r["participants"]], temp="room temperature", time=hours(0.5, 19), atmosphere="nitrogen atmosphere",
                       operation="scope entry under Scheme 3 conditions", method="inherited_general_condition")]
    r["procedure"]["general_procedure_reference"] = f("Scheme 3 footnote a" + ("; footnote c stoichiometry exception" if reverse else ""),
        "Scheme 3 footnote a" + (" with footnote c exception" if reverse else ""), ev)
    r["procedure"]["modifications"] = f("1:2 amounts reversed to 0.6:0.2 mmol" if reverse else None,
        "reverse substrate ratio" if reverse else None, ev) if reverse else missing(ev)
    r["outcomes"] = [outcome("yield", "yield", [f"{DOC}:mol:{label}"], pct(yield_n), ev, scope=f"entry {label}; substituent {substituent}")]
    classification(r, ev, ["C-C", "C-O"], "Partial net classification from the tandem precursor/product drawings; exact atom mapping and multiplicities remain unresolved.")
    r["notes"] = [f(f"Scheme 3 entry {label}: {substituent}; exact row-specific time is not printed, only the 0.5–19 h scope range.",
                    {"entry": label, "substituent": substituent, "time_scope": "0.5–19 h"}, ev)]
    r["review"]["unresolved"] = ["Exact substrate-instance structures and row-specific reaction time are not separately labelled in the main text.",
                                  "Yield measurement basis is not reported in the supplied main text.",
                                  "Dense relative stereochemistry is not encoded as SMILES."]
    reactions.append(r)


scope4 = [
    ("4a",72,"R=H",False),("4b",79,"R=Cl (C6)",False),("4c",70,"R=Br (C6)",False),("4d",53,"R=NO2 (C6)",True),("4e",84,"R=CO2Me (C6)",False),
    ("4f",53,"R=F (C5)",False),("4g",63,"R=Cl (C5)",False),("4h",65,"R=Br (C5)",False),
    ("4i",78,"R1 aryl 4-F",False),("4j",76,"R1 aryl 4-Cl",False),("4k",68,"R1 aryl 4-Br",False),("4l",75,"R1 aryl 4-I",False),
    ("4m",67,"R1 aryl 2-F",False),("4n",54,"R1 aryl 2-Br",False),("4o",69,"R1 aryl 4-Me",False),("4p",76,"R1 aryl 4-OMe",False),("4q",57,"R1 naphthyl",False),
    ("4r",61,"R1 2-thienyl",False),("4s",46,"R1 2-furyl",False),("4t",67,"R2 aryl 3-F",False),("4u",75,"R2 aryl 4-Me",False),
    ("4v",65,"R2 aryl 4-t-Bu",False),("4w",69,"R2 aryl 4-OMe",False),("4x",32,"R2 naphthyl",False),("4y",76,"R2 aryl 2,5-diOMe",False),
    ("4z",74,"mixed aryl substitution shown",False),("4aa",50,"cycloheptyl-fused substrate",False),
]
for row in scope4: standard_scope4(*row)


def standard_scope5(label, yield_n, dr, substituent, footnote=None):
    ev = [EV[3]]
    r = base_reaction(f"scope_{label}", "specific_experiment", ev, f"PDF page 3, Scheme 5, entry {label}",
                      ["substrate_scope"], "[3 + 2] cycloaddition", "gold-catalyzed tandem cycloisomerization/[3+2] cycloaddition")
    r["general_scheme_ids"] = [f"{DOC}:general_3plus2"]
    ref = f"{DOC}:mol:series1"
    r["participants"] = [
        participant("allene", f"{DOC}:mol:series1", "reactant", ev, amount=mmol(0.2), eq=equiv(1.0), ref=ref, charge_method="inherited_general_condition"),
        participant("ylide", f"{DOC}:mol:series3", "reactant", ev, amount=mmol(0.22), eq=equiv(1.1), ref=ref, charge_method="inherited_general_condition"),
        participant("product", f"{DOC}:mol:{label}", "product", ev, product_role="major cycloadduct"),
        participant("gold", f"{DOC}:mol:PPh3AuCl", "catalyst", ev, molpct=molpct(1), ref=ref, charge_method="inherited_general_condition"),
        participant("silver", f"{DOC}:mol:AgSbF6", "additive", ev, molpct=molpct(10), ref=ref, charge_method="inherited_general_condition"),
        participant("sieves", f"{DOC}:mol:4AMS", "additive", ev, mass=mg(100), charge_method="inherited_general_condition"),
    ]
    if footnote == "c":
        r["participants"] += [participant("solvent1", f"{DOC}:mol:dioxane", "solvent", ev, volume=ml(2.0), charge_method="inherited_general_condition"),
                              participant("solvent2", f"{DOC}:mol:DCE", "solvent", ev, volume=ml(2.0), charge_method="inherited_general_condition")]
        set_solvent(r, ev, ["1,4-dioxane", "DCE"], ml(4.0), mixture="2.0 mL/2.0 mL", method="inherited_general_condition")
        temp = deg(80)
    elif footnote == "d":
        r["participants"].append(participant("solvent", f"{DOC}:mol:THF", "solvent", ev, volume=ml(2.0), charge_method="inherited_general_condition"))
        set_solvent(r, ev, "THF", ml(2.0), method="inherited_general_condition")
        temp = "room temperature"
    else:
        r["participants"].append(participant("solvent", f"{DOC}:mol:EA", "solvent", ev, volume=ml(2.0), charge_method="inherited_general_condition"))
        set_solvent(r, ev, "ethyl acetate", ml(2.0), method="inherited_general_condition")
        temp = "room temperature"
    r["reaction_scale"] = f("1 (0.2 mmol)", {"reference_molecule_id": ref, "amount": mmol(0.2)}, ev)
    r["steps"] = [step("reaction", ev, [p["participant_id"] for p in r["participants"]], temp=temp, time=hours(2,12), atmosphere="nitrogen atmosphere",
                       operation="scope entry under Scheme 5 conditions", method="inherited_general_condition")]
    r["procedure"]["general_procedure_reference"] = f("Scheme 5 footnote a" + (f"; footnote {footnote} solvent/temperature exception" if footnote else ""),
        "Scheme 5 footnote a" + (f" with footnote {footnote} exception" if footnote else ""), ev)
    r["outcomes"] = [
        outcome("yield", "yield", [f"{DOC}:mol:{label}"], pct(yield_n), ev, scope=f"entry {label}; {substituent}"),
        outcome("dr", "selectivity", [f"{DOC}:mol:{label}"], dr, ev, scope={"metric":"diastereomeric ratio", "order":"drawn major product : unassigned minor diastereomer"}),
    ]
    classification(r, ev, ["C-C", "C-O", "C-N"], "Partial net classification from tandem precursor/product drawings; exact atom mapping and multiplicities remain unresolved.")
    r["notes"] = [f(f"Scheme 5 entry {label}: {substituent}; the drawn product is the major diastereomer.",
                    {"entry":label,"substituent":substituent,"dr":dr,"footnote":footnote}, ev)]
    r["review"]["unresolved"] = ["Exact substrate-instance structures and exact row-specific reaction time are not separately labelled in the main text.",
                                  "Yield and dr measurement methods are not reported in the supplied main text.",
                                  "Minor diastereomer configuration and dense major-product stereochemistry are not encoded as SMILES."]
    reactions.append(r)


scope5 = [
    ("5a",82,"9:1","R=H",None),("5b",42,"5:1","R=Cl","c"),("5c",61,"5:1","R=Br",None),("5d",65,"6:1","R=CO2Me","c"),
    ("5e",76,"11:1","R'=H, R=F",None),("5f",55,"6:1","R'=Cl, R=F","c"),("5g",36,"6:1","R'=Br, R=F","c"),("5h",66,"4:1","R/R'=-OCH2O-","c"),
    ("5i",74,"8:1","SO2R, R=thienyl",None),("5j",80,"11:1","SO2R, R=4-FC6H4",None),("5k",81,"9:1","SO2R, R=4-ClC6H4",None),
    ("5l",73,"9:1","SO2R, R=4-OMeC6H4",None),("5m",72,"9:1","SO2R, R=beta-naphthyl",None),
    ("5n",79,"9:1","R1 aryl 4-F",None),("5o",81,"11:1","R1 aryl 4-Cl",None),("5p",76,"11:1","R1 aryl 4-Br",None),
    ("5q",55,"8:1","R1 aryl 2-F",None),("5r",64,"5:1","R1 aryl 2-Br",None),("5s",82,"9:1","R1 aryl 4-Me",None),
    ("5t",59,"10:1","R1 aryl 4-OMe",None),("5u",77,"11:1","R1 aryl 4-acetyl",None),("5v",60,"3:1","R1 alkyl (isopropyl shown)",None),
    ("5w",79,"16:1","R1 2-thienyl",None),("5x",77,"12:1","R1 2-furyl",None),
    ("5y",73,"10:1","R2 aryl 2-Me",None),("5z",62,"9:1","R2 aryl 3-F",None),("5aa",80,"9:1","R2 aryl 4-Me",None),
    ("5ab",75,"9:1","R2 aryl 4-OMe",None),("5ac",69,"9:1","R2 aryl 4-Ac",None),("5ad",78,"9:1","R2 aryl 4-t-Bu",None),
    ("5ae",71,"8:1","cycloheptyl-fused, n=2","d"),("5af",40,"4:1","cyclooctyl-fused, n=3","d"),
]
for row in scope5: standard_scope5(*row)


# Explicit optimization and gram-scale experiments not duplicated with scope entries.
def opt4_initial():
    ev=[EV[2],TEV[2]]; ref=f"{DOC}:mol:1a"
    r=base_reaction("opt_4a_initial_36", "specific_experiment", ev, "PDF page 2, Scheme 2 initial line and prose", ["optimization"], "[4 + 2] cycloaddition", "initial-condition tandem cycloaddition")
    r["general_scheme_ids"]=[f"{DOC}:general_4plus2"]
    r["participants"]=[participant("1a",ref,"reactant",ev,amount=mmol(.2),eq=equiv(1),ref=ref),participant("2a",f"{DOC}:mol:2a","reactant",ev,amount=mmol(.6),eq=equiv(3),ref=ref),participant("4a",f"{DOC}:mol:4a","product",ev,product_role="cycloadduct"),participant("gold",f"{DOC}:mol:PPh3AuOTf","catalyst",ev,molpct=molpct(10),ref=ref),participant("sieves",f"{DOC}:mol:4AMS","additive",ev),participant("solvent",f"{DOC}:mol:DCE","solvent",ev)]
    set_solvent(r,ev,"DCE",concentration={"number":.1,"unit":"M","raw_text":"0.1 M"})
    r["reaction_scale"]=f("1a, 0.2 mmol",{"reference_molecule_id":ref,"amount":mmol(.2)},ev)
    r["steps"]=[step("reaction",ev,[p["participant_id"] for p in r["participants"]],temp="room temperature",time=hours(12),atmosphere="N2 atmosphere",operation="treated under initial conditions")]
    r["outcomes"]=[outcome("yield","yield",[f"{DOC}:mol:4a"],pct(36),ev),outcome("stereoselectivity","selectivity",[f"{DOC}:mol:4a"],"excellent stereoselectivity",ev,scope={"metric":"qualitative stereoselectivity","configuration":"not quantified"})]
    classification(r,ev,["C-C","C-O"],"Partial net classification from source drawings; exact atom mapping unresolved.")
    r["review"]["unresolved"]=["4 Å MS mass and yield/selectivity measurement methods are unreported.","Qualitative stereoselectivity is not numerically defined."]
    reactions.append(r)
opt4_initial()


def opt4_silver():
    ev=[TEV[2],EV[2]]; ref=f"{DOC}:mol:1a"
    r=base_reaction("opt_4a_AgSbF6_only_55", "specific_experiment", ev, "PDF page 2 prose: AgSbF6 sole-catalyst check", ["optimization","control"], "[4 + 2] cycloaddition", "silver-only optimization")
    r["participants"]=[participant("1a",ref,"reactant",ev),participant("2a",f"{DOC}:mol:2a","reactant",ev),participant("4a",f"{DOC}:mol:4a","product",ev),participant("silver",f"{DOC}:mol:AgSbF6","catalyst",ev,molpct=molpct(5),ref=ref),participant("solvent",f"{DOC}:mol:DCM","solvent",ev)]
    set_solvent(r,ev,"DCM")
    r["steps"]=[step("reaction",ev,[p["participant_id"] for p in r["participants"]],operation="AgSbF6 used as sole catalyst")]
    r["outcomes"]=[outcome("yield","yield",[f"{DOC}:mol:4a"],pct(55),ev)]
    classification(r,ev,["C-C","C-O"],"Partial net classification from source drawings; exact atom mapping unresolved.")
    r["review"]["unresolved"]=["Reactant charges, temperature, time, atmosphere, solvent volume, and yield basis are not restated in the main-text sentence; Table S1 is unavailable."]
    reactions.append(r)
opt4_silver()


def gram4():
    ev=[EV[2]];ref=f"{DOC}:mol:1a"
    r=base_reaction("scaleup_4a_77", "specific_experiment", ev, "PDF page 2, Scheme 3 footnote b", ["scale_up"], "[4 + 2] cycloaddition", "gram-scale tandem cycloaddition")
    r["general_scheme_ids"]=[f"{DOC}:general_4plus2"]
    r["participants"]=[participant("1a",ref,"reactant",ev,mass=gram(1.130)),participant("2a",f"{DOC}:mol:2a","reactant",ev,mass=gram(1.641)),participant("4a",f"{DOC}:mol:4a","product",ev,product_role="isolated material"),participant("gold",f"{DOC}:mol:tBuXPhosAuCl","catalyst",ev,molpct=molpct(2),ref=ref,charge_method="inherited_general_condition"),participant("silver",f"{DOC}:mol:AgSbF6","additive",ev,molpct=molpct(15),ref=ref,charge_method="inherited_general_condition"),participant("sieves",f"{DOC}:mol:4AMS","additive",ev),participant("solvent",f"{DOC}:mol:DCE","solvent",ev)]
    r["reaction_scale"]=f("1a, 1.130 g",{"reference_molecule_id":ref,"mass":gram(1.130)},ev)
    set_solvent(r,ev,"DCE",method="inherited_general_condition")
    r["steps"]=[step("reaction",ev,[p["participant_id"] for p in r["participants"]],temp="room temperature",atmosphere="nitrogen atmosphere",operation="gram-scale reaction under Scheme 3 conditions",method="inherited_general_condition")]
    r["outcomes"]=[outcome("yield","yield",[f"{DOC}:mol:4a"],pct(77),ev),outcome("product_mass","product_amount",[f"{DOC}:mol:4a"],gram(1.340),ev)]
    classification(r,ev,["C-C","C-O"],"Partial net classification from source drawings; exact atom mapping unresolved.")
    r["review"]["unresolved"]=["Footnote b does not state the exact reaction time or yield measurement basis.","Catalyst loadings, solvent identity, room temperature, and nitrogen are inherited from Scheme 3; scale-up-specific sieves mass and solvent volume are unreported and were not copied from the 0.2 mmol footnote-a experiment."]
    reactions.append(r)
gram4()


def opt5(rid, initial, yield_n, dr):
    ev=[EV[3],TEV[3]];ref=f"{DOC}:mol:1a"
    r=base_reaction(rid,"specific_experiment",ev,"PDF page 3, Scheme 4 " + ("initial" if initial else "optimized") + " line",["optimization"],"[3 + 2] cycloaddition","condition-screening tandem cycloaddition")
    r["general_scheme_ids"]=[f"{DOC}:general_3plus2"]
    r["participants"]=[participant("1a",ref,"reactant",ev,amount=mmol(.2),eq=equiv(1),ref=ref),participant("3a",f"{DOC}:mol:3a","reactant",ev,amount=mmol(.22),eq=equiv(1.1),ref=ref),participant("5a",f"{DOC}:mol:5a","product",ev,product_role="major cycloadduct"),participant("gold",f"{DOC}:mol:PPh3AuOTf" if initial else f"{DOC}:mol:PPh3AuCl","catalyst",ev,molpct=molpct(10 if initial else 1),ref=ref),participant("sieves",f"{DOC}:mol:4AMS","additive",ev)]
    if not initial:r["participants"].append(participant("silver",f"{DOC}:mol:AgSbF6","additive",ev,molpct=molpct(10),ref=ref))
    solvent="DCE" if initial else "EA";r["participants"].append(participant("solvent",f"{DOC}:mol:{solvent}","solvent",ev))
    set_solvent(r,ev,"DCE" if initial else "ethyl acetate",concentration={"number":.1,"unit":"M","raw_text":"0.1 M"})
    r["reaction_scale"]=f("1a, 0.2 mmol",{"reference_molecule_id":ref,"amount":mmol(.2)},ev)
    r["steps"]=[step("reaction",ev,[p["participant_id"] for p in r["participants"]],temp="room temperature",operation="initial attempt" if initial else "optimized condition")]
    r["outcomes"]=[outcome("yield","yield",[f"{DOC}:mol:5a"],pct(yield_n),ev),outcome("dr","selectivity",[f"{DOC}:mol:5a"],dr,ev,scope={"metric":"diastereomeric ratio","order":"major:minor; configurations unassigned"})]
    classification(r,ev,["C-C","C-O","C-N"],"Partial net classification from source drawings; exact atom mapping unresolved.")
    r["notes"]=[f("The page-3 prose calls the model ylide '2a', but Scheme 4 and Scheme 5 label the drawn ylide 3a; the participant is contextually linked to 3a and the literal conflict is retained.",
                  {"prose_label":"2a","scheme_label":"3a","status":"conflict"},ev)]
    r["review"]["unresolved"]=["Source conflict: page-3 prose labels the model ylide 2a, while Schemes 4/5 label it 3a.","Reaction time, 4 Å MS mass, atmosphere, yield basis, and dr method are not reported for Scheme 4 in the supplied main text."]
    reactions.append(r)
opt5("opt_5a_initial_40",True,40,"8:1")
opt5("opt_5a_optimized_83",False,83,"9:1")


def gram5():
    ev=[EV[3]];ref=f"{DOC}:mol:1a"
    r=base_reaction("scaleup_5a_85","specific_experiment",ev,"PDF page 3, Scheme 5 footnote b",["scale_up"],"[3 + 2] cycloaddition","gram-scale tandem cycloaddition")
    r["general_scheme_ids"]=[f"{DOC}:general_3plus2"]
    r["participants"]=[participant("1a",ref,"reactant",ev,mass=gram(.6617)),participant("3a",f"{DOC}:mol:3a","reactant",ev,mass=gram(.8241)),participant("5a_contextual",f"{DOC}:mol:5a","product",ev,product_role="contextually inferred Scheme 5 product"),participant("gold",f"{DOC}:mol:PPh3AuCl","catalyst",ev,molpct=molpct(1),ref=ref,charge_method="inherited_general_condition"),participant("silver",f"{DOC}:mol:AgSbF6","additive",ev,molpct=molpct(10),ref=ref,charge_method="inherited_general_condition"),participant("sieves",f"{DOC}:mol:4AMS","additive",ev),participant("solvent",f"{DOC}:mol:EA","solvent",ev)]
    r["reaction_scale"]=f("1a, 0.6617 g",{"reference_molecule_id":ref,"mass":gram(.6617)},ev)
    set_solvent(r,ev,"ethyl acetate",method="inherited_general_condition")
    r["steps"]=[step("reaction",ev,[p["participant_id"] for p in r["participants"]],temp="room temperature",atmosphere="nitrogen atmosphere",operation="gram-scale reaction under Scheme 5 conditions",method="inherited_general_condition")]
    r["outcomes"]=[outcome("yield","yield",[f"{DOC}:mol:5a"],pct(85),ev,scope={"literal_source_target":"4a","contextual_scheme_target":"5a","status":"conflict"}),outcome("product_mass","product_amount",[f"{DOC}:mol:5a"],gram(1.161),ev,scope={"literal_source_target":"4a","contextual_scheme_target":"5a","status":"conflict"})]
    classification(r,ev,["C-C","C-O","C-N"],"Partial net classification from Scheme 5 contextual product assignment; exact atom mapping unresolved.")
    r["notes"]=[f("Footnote b literally says 'providing 4a (1.161 g) in 85% yield', while the Scheme 5 reaction of 1a + 3a and parenthetical beside 5a contextually indicate 5a. No silent correction was made.",{"literal":"4a","contextual":"5a","product_mass":"1.161 g","yield":"85%"},ev)]
    r["review"]["unresolved"]=["Source conflict: gram-scale footnote names product 4a, while Scheme 5 context indicates 5a.","Exact time and yield/dr measurement methods are not reported.","Catalyst loadings, solvent identity, room temperature, and nitrogen are inherited from Scheme 5; scale-up-specific sieves mass and solvent volume are unreported and were not copied from the 0.2 mmol experiment."]
    reactions.append(r)
gram5()


# Scheme 6 derivatization arrows.
deriv_specs = [
    ("deriv_4a_to_6","4a","6",89,[('Zn','reagent',equiv(30)),('HOAc','reagent',equiv(5)),('DCM','solvent',None)],"room temperature",["N-H","O-H"],"reductive N-O bond cleavage"),
    ("deriv_6_to_7a","6","7a",62,[('propargyl_bromide','reagent',equiv(4)),('K2CO3','reagent',equiv(2)),('MeCN','solvent',None)],deg(60),["C-N"],"N-alkylation"),
    ("deriv_6_to_7b","6","7b",71,[('IBX','reagent',equiv(1.1)),('DMSO','solvent',None)],"room temperature",["DB-C-O"],"oxidation"),
    ("deriv_5a_to_8","5a","8",92,[('TfOH','reagent',equiv(2)),('DCE','solvent',None)],deg(-20),["C-C"],"ring expansion"),
    ("deriv_8_to_9","8","9",48,[('RaneyNi','catalyst',None),('H2','reagent',None),('MeOH','solvent',None)],"room temperature",["C-H","N-H"],"Raney Ni/H2 reduction"),
]
for rid,sm,pm,yld,agents,temp,codes,name in deriv_specs:
    ev=[EV[3]];r=base_reaction(rid,"specific_experiment",ev,f"PDF page 3, Scheme 6, {sm} to {pm}",["derivatization"],name,name)
    r["participants"]=[participant("substrate",f"{DOC}:mol:{sm}","reactant",ev),participant("product",f"{DOC}:mol:{pm}","product",ev,product_role="reported product")]
    for i,(mid,role,charge) in enumerate(agents,1):
        r["participants"].append(participant(f"agent{i}",f"{DOC}:mol:{mid}",role,ev,eq=charge,ref=f"{DOC}:mol:{sm}" if charge else None))
    solvent_name=next((mid for mid,role,_ in agents if role=='solvent'),None)
    if solvent_name:set_solvent(r,ev,solvent_name)
    r["steps"]=[step("reaction",ev,[p["participant_id"] for p in r["participants"]],temp=temp,operation=name)]
    r["outcomes"]=[outcome("yield","yield",[f"{DOC}:mol:{pm}"],pct(yld),ev)]
    classification(r,ev,codes,"Partial category-level inference from the Scheme 6 before/after drawings; hydrogen source or complete atom mapping is not asserted.")
    r["review"]["unresolved"]=["Reaction time, atmosphere, scale, solvent volume, workup/purification, and yield measurement basis are unreported in the supplied main text.","Dense structures are not encoded as complete SMILES; bond classification is partial."]
    reactions.append(r)


# Scheme 7 six individually drawn experiments and footnote-specific conditions.
scheme7 = [
    ("10","11",62,"a",["C-C","C-O"]),("12","13",60,"b",["C-C","C-O"]),("14","15",94,"c",["C-C","C-O"]),
    ("16","17",80,"c",["C-C","C-O"]),("18","19",92,"c",["C-C","C-O","C-N"]),("20","21",48,"d",["C-C","C-O"]),
]
for partner,prod,yld,fn,codes in scheme7:
    ev=[EV[4]];ref=f"{DOC}:mol:1a"
    r=base_reaction(f"scheme7_{partner}_to_{prod}","specific_experiment",ev,f"PDF page 4, Scheme 7, {partner} to {prod}",["substrate_scope"],"cross cycloaddition with other coupling partner","gold/silver-enabled tandem cycloaddition")
    if fn=='a': a_amt,p_amt,p_eq=mmol(.1),mmol(.15),equiv(1.5)
    else: a_amt,p_amt,p_eq=mmol(.2),mmol(.6 if fn in ('b','d') else .2),equiv(3 if fn in ('b','d') else 1)
    r["participants"]=[participant("1a",ref,"reactant",ev,amount=a_amt,eq=equiv(1),ref=ref),participant("partner",f"{DOC}:mol:{partner}","reactant",ev,amount=p_amt,eq=p_eq,ref=ref),participant("product",f"{DOC}:mol:{prod}","product",ev,product_role="reported cycloadduct")]
    if fn=='a':
        r["participants"] += [participant("gold",f"{DOC}:mol:tBuXPhosAuCl","catalyst",ev,molpct=molpct(2),ref=ref),participant("silver",f"{DOC}:mol:AgTFA","additive",ev,molpct=molpct(15),ref=ref),participant("sieves",f"{DOC}:mol:4AMS","additive",ev,mass=mg(100)),participant("solvent",f"{DOC}:mol:DCE","solvent",ev,volume=ml(2))]
        temp="room temperature";atm="nitrogen atmosphere";solv="DCE"
    elif fn=='b':
        r["participants"] += [participant("gold",f"{DOC}:mol:PPh3AuCl","catalyst",ev,molpct=molpct(2),ref=ref),participant("silver",f"{DOC}:mol:AgTFA","additive",ev,molpct=molpct(15),ref=ref),participant("sieves",f"{DOC}:mol:4AMS","additive",ev,mass=mg(100)),participant("solvent",f"{DOC}:mol:DCE","solvent",ev,volume=ml(2))]
        temp=None;atm=None;solv="DCE"
    elif fn=='c':
        r["participants"] += [participant("gold",f"{DOC}:mol:PPh3AuOTf","catalyst",ev,molpct=molpct(10),ref=ref),participant("sieves",f"{DOC}:mol:4AMS","additive",ev,mass=mg(100)),participant("solvent",f"{DOC}:mol:DCE","solvent",ev,volume=ml(2))]
        temp=None;atm=None;solv="DCE"
    else:
        r["participants"] += [participant("silver",f"{DOC}:mol:AgSbF6","catalyst",ev,molpct=molpct(10),ref=ref),participant("sieves",f"{DOC}:mol:4AMS","additive",ev,mass=mg(100)),participant("solvent",f"{DOC}:mol:THF","solvent",ev,volume=ml(2))]
        temp=None;atm=None;solv="THF"
    r["reaction_scale"]=f(f"1a ({a_amt['raw_text']})",{"reference_molecule_id":ref,"amount":a_amt},ev)
    set_solvent(r,ev,solv,ml(2))
    r["steps"]=[step("reaction",ev,[p["participant_id"] for p in r["participants"]],temp=temp,atmosphere=atm,operation=f"Scheme 7 footnote {fn} conditions")]
    r["outcomes"]=[outcome("yield","yield",[f"{DOC}:mol:{prod}"],pct(yld),ev)]
    classification(r,ev,codes,"Partial category-level inference from Scheme 7 net drawings; tandem cycloisomerization and cycloaddition are not fully atom mapped.")
    r["review"]["unresolved"]=["Reaction time, workup/purification, and yield measurement basis are unreported.","Temperature and atmosphere are unreported for this footnote." if fn!='a' else "Exact stereochemical selectivity is not quantified.","Dense structures and exact bond multiplicities are unresolved."]
    reactions.append(r)

bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id":f"{DOC}:context:abstract_graphic","document_id":DOC,"evidence_ids":[EV[1]],"kind":"duplicate_overview","note":"Abstract graphic summarizes Schemes 3, 5, and 7; it is not counted as separate experiments."},
    {"context_id":f"{DOC}:context:natural_products","document_id":DOC,"evidence_ids":[EV[2]],"kind":"background","note":"Scheme 1A natural-product/supramolecular examples are background structures, not experiments in this paper."},
    {"context_id":f"{DOC}:context:optimization_tables","document_id":DOC,"evidence_ids":[TEV[2],TEV[3]],"kind":"unavailable_supporting_information","note":"Main text summarizes Tables S1/S2, but no SI PDF was supplied or associated; unshown rows were not fabricated."},
    {"context_id":f"{DOC}:context:source_label_conflicts","document_id":DOC,"evidence_ids":[TEV[2],TEV[3],EV[3],EV[4]],"kind":"source_conflict","note":"Literal inconsistencies retained: page-2 optimization prose says product 3a where the discussed product is otherwise 4a; page-3 prose calls the model ylide 2a while Schemes 4/5 label it 3a; Scheme 5 title says [4+2] while its body/prose describe [3+2]; Scheme 5 gram-scale footnote says product 4a while context indicates 5a; page-3 Scheme 7 prose says substrates 16 and 17 furnish 17 and 19, while the page-4 drawing shows 16→17 and 18→19."},
    {"context_id":f"{DOC}:context:xray","document_id":DOC,"evidence_ids":[EV[2],EV[3],EV[4]],"kind":"characterization","note":"CCDC 2073826 and 2073827 are shown for 4a and 5a; additional accession codes 1962639 and 2073828 are announced. These are not synthesis records."},
    {"context_id":f"{DOC}:context:mechanistic_discussion","document_id":DOC,"evidence_ids":[TEV[3],TEV[4]],"kind":"author_interpretation","note":"Authors discuss inhibited ring opening to a 1,4-dipole in O-containing solvent; this is interpretation, not a separate mechanistic experiment."},
    {"context_id":f"{DOC}:context:si_notice","document_id":DOC,"evidence_ids":[EV[4]],"kind":"si_notice","note":"The article announces experimental procedures and spectra in SI, but association remains unestablished because no SI file was provided."},
    {"context_id":f"{DOC}:context:bibliography","document_id":DOC,"evidence_ids":[EV[4],EV[5]],"kind":"bibliography","note":"References occupy page 5 (and part of page 4); they are not current-work experiments."},
]

(OUT / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
(OUT / "raw_model" / "stereo_diagnostics.json").write_text(json.dumps(stereo_diagnostics, ensure_ascii=False, indent=2) + "\n")

page2_ids = [r["reaction_id"] for r in reactions if EV[2] in r["evidence_ids"]]
page3_ids = [r["reaction_id"] for r in reactions if EV[3] in r["evidence_ids"]]
page4_ids = [r["reaction_id"] for r in reactions if EV[4] in r["evidence_ids"]]
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page":1,"state":"attempted","regions":[
            {"region":"title, abstract, and abstract graphical summary","finding":"Title/metadata captured; graphical summary duplicates the current-work scope and is not counted as new experiments.","evidence_ids":[EV[1],TEV[1]]},
            {"region":"introduction columns","finding":"Background and citations only; no present-work concrete experiment arrow on this page.","evidence_ids":[EV[1],TEV[1]]}],
         "record_ids":[],"unresolved":["Small structures in the abstract graphic are treated only as duplicate overview; exact structures are handled from Schemes 3, 5, and 7."]},
        {"page":2,"state":"attempted","regions":[
            {"region":"Scheme 1A-C","finding":"Natural-product background, two cited prior-work arrows, and an author working hypothesis/in-situ intermediate were distinguished from current experiments.","evidence_ids":[EV[2]]},
            {"region":"Scheme 2 and adjacent optimization prose","finding":"Initial 36% experiment, optimized/scope 72% result, and explicit AgSbF6-only 55% check were captured without reconstructing unavailable Table S1.","evidence_ids":[EV[2],TEV[2]]},
            {"region":"Scheme 3 including footnotes","finding":"All 27 scope entries 4a-4aa and the separate 4a gram-scale run were mapped; 4d footnote-c stoichiometry reversal retained.","evidence_ids":[EV[2]]}],
         "record_ids":page2_ids,"unresolved":["Individual exact scope times within 0.5–19 h, yield basis, detailed procedures, and stereochemically complete structures require SI/source review."]},
        {"page":3,"state":"attempted","regions":[
            {"region":"Scheme 4 and prose","finding":"Initial 40%/8:1 and optimized 83%/9:1 experiments captured separately.","evidence_ids":[EV[3],TEV[3]]},
            {"region":"Scheme 5 including footnotes","finding":"All 32 scope entries 5a-5af captured with dr; footnote-c and footnote-d condition exceptions retained; gram-scale footnote product-label conflict preserved.","evidence_ids":[EV[3]]},
            {"region":"Scheme 6 five arrows","finding":"Five independent derivatization reactions captured in arrow order with all printed reagent equivalents and temperatures.","evidence_ids":[EV[3]]},
            {"region":"Scheme 7 discussion continuation","finding":"Narrative duplicates yields shown on page 4 and adds author interpretation; it was linked as context, not duplicate experiments.","evidence_ids":[EV[3],TEV[3]]}],
         "record_ids":page3_ids,"unresolved":["Exact row-specific times within 2–12 h, yield/dr methods, workups, and complete stereochemical structures are unavailable; footnote b literally names 4a although Scheme 5 context indicates 5a."]},
        {"page":4,"state":"attempted","regions":[
            {"region":"Scheme 7 table and footnotes a-d","finding":"Six distinct other-partner experiments (10→11, 12→13, 14→15, 16→17, 18→19, 20→21) captured with footnote-specific charges/catalysts/solvents.","evidence_ids":[EV[4]]},
            {"region":"mechanistic discussion and conclusion","finding":"Author interpretation and conclusion retained as source context, not synthetic experiments.","evidence_ids":[EV[4],TEV[4]]},
            {"region":"Supporting Information and accession-code notices","finding":"SI is announced but not provided/associated; CCDC notices are characterization context, not reactions.","evidence_ids":[EV[4],TEV[4]]},
            {"region":"author/acknowledgment/reference column","finding":"Metadata and bibliography only; no additional current experiment.","evidence_ids":[EV[4]]}],
         "record_ids":page4_ids,"unresolved":["Scheme 7 footnotes b-d omit time, temperature/atmosphere (except footnote a), yield basis, and workup; no SI was supplied."]},
        {"page":5,"state":"attempted","regions":[
            {"region":"two-column bibliography","finding":"References 3-15 only; no present-work experiment, standalone spectrum, or mass-peak table.","evidence_ids":[EV[5],TEV[5]]}],
         "record_ids":[],"unresolved":["Cited papers were not independently inspected; bibliography is context only."]},
    ]
}
(OUT / "coverage").mkdir(exist_ok=True)
(OUT / "coverage" / "pages_0001_0005.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"molecules":len(molecules),"reactions":len(reactions),"stereo_diagnostics":len(stereo_diagnostics)},ensure_ascii=False))
