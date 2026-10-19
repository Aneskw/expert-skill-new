#!/usr/bin/env python3
"""Build source-derived v3.1 records for PDF pages 1--6.

This script is an extraction artifact, not a chemistry verifier.  It preserves
printed inconsistencies and leaves unresolved structures/conditions explicit.
"""
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DATASET = st.read(RUN / "dataset.json")
DOC = "doc_e6209d110a05e1c9"


def ev(page, kind="image"):
    return f"{DOC}:page{page}:{kind}"


def f(raw=None, normalized=None, pages=(), method="original_extraction", derivation=None,
      status=None, alternatives=()):
    return st.fact(raw, normalized, [ev(p, "image") if p <= 6 else ev(p, "text") for p in pages],
                   method=method, derivation=derivation, status=status,
                   alternatives=alternatives)


def inferred(normalized, pages, derivation):
    return f(None, normalized, pages, method="inference", derivation=derivation)


def qty(number, unit, raw=None):
    result = {"number": number, "unit": unit}
    if raw is not None:
        result["raw_text"] = raw
    return result


def molecule(mid, label, name, pages, smiles=None, chirality="unresolved",
             geometry="unresolved", representation="unresolved", mixture=None,
             atom_centers=None, assignment=None, identity_note=None):
    evidence = [ev(p, "image") if p <= 6 else ev(p, "text") for p in pages]
    m = st.new_molecule(f"{DOC}:mol:{mid}", DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = st.fact(label, label, evidence)
    m["names"]["original"] = st.fact(name, name, evidence)
    m["names"]["english"] = st.fact(name, name, evidence)
    m["names"]["abbreviations"] = st.fact(label, [label], evidence)
    m["names"]["aliases"] = st.fact()
    m["cas_numbers"] = st.fact()
    if smiles:
        diag = st.stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = st.fact(None, diag["canonical_isomeric_smiles"], evidence,
                              method="structure_derivation",
                              derivation="Candidate connectivity transcribed from the source drawing/name and canonicalized by RDKit; parsing is not source verification.")
        m["generic_template"] = None
        m["stereochemistry"]["isomeric_smiles"] = st.fact(
            None, diag["canonical_isomeric_smiles"], evidence, method="structure_derivation",
            derivation="Canonical isomeric SMILES from the candidate graph; unassigned elements remain unassigned.")
        m["stereochemistry"]["identity_key"] = st.fact(
            None, diag["stereo_graph_key"], evidence, method="structure_derivation",
            derivation="RDKit stereo_graph_key from the candidate graph; not a sample identity.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = st.fact()
        m["generic_template"] = None
    for key, value in (("chirality_status", chirality), ("geometry_status", geometry),
                       ("representation_status", representation)):
        m["stereochemistry"][key] = st.fact(
            None, value, evidence, method="inference",
            derivation="Status assigned after inspecting the source drawing/name; no unreported stereoisomer or composition was inferred.")
    if mixture is not None:
        m["stereochemistry"]["mixture"] = st.fact(mixture, mixture, evidence)
    if atom_centers is not None:
        m["stereochemistry"]["atom_centers"] = st.fact(atom_centers, atom_centers, evidence)
    if assignment is not None:
        m["stereochemistry"]["assignment"] = st.fact(assignment, assignment, evidence)
    m["identity_note"] = st.fact(
        identity_note, identity_note, evidence) if identity_note else st.fact()
    m["external_properties"] = []
    return m


def generic_molecule(mid, label, fixed, variables, pages, reason):
    m = molecule(mid, label, label, pages, chirality="unresolved",
                 geometry="unresolved", representation="not_encoded", identity_note=reason)
    evidence = [ev(p, "image") for p in pages]
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": st.fact(fixed, fixed, evidence),
        "variables": [{
            "label": v,
            "attachment_points": st.fact("shown in Scheme 1A", "shown in Scheme 1A", evidence),
            "definition": st.fact(desc, desc, evidence),
            "allowed_values": st.fact(vals, vals, evidence),
            "constraints": st.fact()
        } for v, desc, vals in variables],
        "instance_mapping": [],
        "reason": st.fact(reason, reason, evidence)
    }
    return m


M = []

# Concrete allene substrates.  Terminal CH2 makes allene E/Z/axial assignment
# inapplicable; where a tetrahedral carbon is potentially stereogenic, it is
# retained as unspecified rather than racemic.
allene_data = {
    "2a": ("methyl 2-(dimethylamino)-2-(trifluoromethyl)penta-3,4-dienoate",
           "COC(=O)C(N(C)C)(C(F)(F)F)C=C=C", "unspecified", [1, 2, 3, 5, 6]),
    "2b": ("methyl 2-(benzyl(methyl)amino)-2-(trifluoromethyl)penta-3,4-dienoate",
           "COC(=O)C(N(C)Cc1ccccc1)(C(F)(F)F)C=C=C", "unspecified", [2, 3, 5]),
    "2c": ("diethyl 2-(dimethylamino)-2-(propa-1,2-dien-1-yl)malonate",
           "CCOC(=O)C(N(C)C)(C(=O)OCC)C=C=C", "achiral", [2, 3, 5]),
    "2d": ("diethyl phosphonate allene 2d",
           "CCOP(=O)(OCC)C(N(C)C)(C(F)(F)F)C=C=C", "unspecified", [2, 3, 5]),
    "2e": ("methyl propargylamino CF3 allenylic ester 2e",
           "COC(=O)C(N(C)CC#C)(C(F)(F)F)C=C=C", "unspecified", [2, 4, 5]),
    "2f": ("diethyl propargylamino allenylmalonate 2f",
           "CCOC(=O)C(N(C)CC#C)(C(=O)OCC)C=C=C", "achiral", [2, 4, 5]),
    "2g": ("bis(trifluoromethyl) propargylamino allene 2g",
           "C#CCN(C)C(C(F)(F)F)(C(F)(F)F)C=C=C", "achiral", [2, 4, 5]),
}
for label, (name, smi, chiral, pages) in allene_data.items():
    M.append(molecule(label, label, name, pages, smi, chiral, "not_applicable",
                      "partial" if chiral == "unspecified" else "not_applicable",
                      atom_centers=("potential tetrahedral center; configuration not reported" if chiral == "unspecified" else None),
                      assignment="No allene axial configuration is assigned; the depicted allene terminus is CH2."))

product_names = {
    "3a": "dimethyl 2,2'-(cyclobutane-1,3-diylidenebis(methanylylidene))bis(2-(dimethylamino)-3,3,3-trifluoropropanoate)",
    "3b": "dimethyl 2,2'-(cyclobutane-1,3-diylidenebis(methanylylidene))bis(2-(benzyl(methyl)amino)-3,3,3-trifluoropropanoate)",
    "3c": "tetraethyl 2,2'-(cyclobutane-1,3-diylidenebis(methanylylidene))bis(2-(dimethylamino)malonate)",
    "3d": "tetraethyl cyclobutane-1,3-diylidenebis(2-(dimethylamino)-1,1,1-trifluoropropan-2-yl-3-ylidene)bis(phosphonate)",
}
ratios3 = {"3a": "Z/E = 1:1; diastereomers A/B = 1:1", "3b": "Z/E = 1:1; diastereomers = 1:1",
           "3c": "Z/E = 2:3", "3d": "Z/E = 1:1; diastereomers A/B = 1:1"}
for label, name in product_names.items():
    ch = "mixture" if label in ("3a", "3b", "3d") else "achiral"
    M.append(molecule(label, label, name, [2, 3, 4, 6], chirality=ch,
                      geometry="mixture", representation="not_encoded",
                      mixture=ratios3[label],
                      identity_note="Registry entry represents the reported product mixture, not any separated E/Z or diastereomer."))

# Separate crystallographic stereoisomers are not merged into the product mixtures.
M.append(molecule("3a_SR_E", "(S,R)-E-3a", "(S,R)-E stereoisomer of 3a", [4],
                  chirality="specified", geometry="specified", representation="not_encoded",
                  atom_centers=[{"assignment": "S,R", "source": "Fig. 2 caption"}],
                  assignment="X-ray structure; CCDC 2081066"))
M.append(molecule("3b_SR_E", "(S,R)-E-3b", "(S,R)-E stereoisomer of 3b", [4],
                  chirality="specified", geometry="specified", representation="not_encoded",
                  atom_centers=[{"assignment": "S,R", "source": "Fig. 2 caption"}],
                  assignment="X-ray structure; CCDC 2081065"))

names4 = {
    "4a": "methyl 2-(methyl((1-((pivaloyloxy)methyl)-1H-1,2,3-triazol-4-yl)methyl)amino)-2-(trifluoromethyl)penta-3,4-dienoate",
    "4b": "methyl 2-(((1-benzyl-1H-1,2,3-triazol-4-yl)methyl)(methyl)amino)-2-(trifluoromethyl)penta-3,4-dienoate",
    "4c": "diethyl 2-(methyl((1-((pivaloyloxy)methyl)-1H-1,2,3-triazol-4-yl)methyl)amino)-2-(propa-1,2-dien-1-yl)malonate",
    "4d": "diethyl 2-(((1-benzyl-1H-1,2,3-triazol-4-yl)methyl)(methyl)amino)-2-(propa-1,2-dien-1-yl)malonate",
    "4e": "(4-((methyl(1,1,1-trifluoro-2-(trifluoromethyl)penta-3,4-dien-2-yl)amino)methyl)-1H-1,2,3-triazol-1-yl)methyl pivalate",
    "4f": "N-((1-benzyl-1H-1,2,3-triazol-4-yl)methyl)-1,1,1-trifluoro-N-methyl-2-(trifluoromethyl)penta-3,4-dien-2-amine",
}
for label, name in names4.items():
    ch = "unspecified" if label in ("4a", "4b") else "achiral"
    M.append(molecule(label, label, name, [4, 5, 6], chirality=ch,
                      geometry="not_applicable", representation="not_encoded",
                      atom_centers=("potential tetrahedral center; configuration not reported" if ch == "unspecified" else None),
                      identity_note="Connectivity is source-identified, but no concrete SMILES is asserted because triazole connectivity/stereochemistry was not independently structure-verified."))

ratios5 = {"5a": "Z/E = 1:1", "5b": "Z/E = 3:1", "5c": "Z/E = 1:2",
           "5d": "Z/E = 1:1", "5e": "Z/E = 3:1", "5f": "Z/E = 4:1"}
for label, ratio in ratios5.items():
    ch = "mixture" if label in ("5a", "5b", "5e", "5f") else "achiral"
    M.append(molecule(label, label, f"1,3-dimethylenecyclobutane dimerization product {label}", [5, 7],
                      chirality=ch, geometry="mixture", representation="not_encoded",
                      mixture=ratio,
                      identity_note="Reported Z/E product mixture; exact isomeric SMILES intentionally not invented."))

# Reagents, solvents, and generic prior-work structures.
simple_specs = [
    ("toluene", "toluene", "Cc1ccccc1"), ("DCE", "1,2-dichloroethane", "ClCCCl"),
    ("water", "water", "O"), ("tBuOH", "tert-butanol", "CC(C)(C)O"),
    ("ethyl_acetate", "ethyl acetate", "CCOC(C)=O"), ("HCl", "hydrochloric acid", "Cl"),
    ("MgSO4", "magnesium sulfate", "[Mg+2].[O-]S(=O)(=O)[O-]"),
    ("benzyl_azide", "benzyl azide", "[N-]=[N+]=NCC1=CC=CC=C1"),
    ("pivoxymethyl_azide", "pivaloyloxymethyl azide", "CC(C)(C)C(=O)OCN=[N+]=[N-]"),
]
for mid, name, smi in simple_specs:
    M.append(molecule(mid, name, name, [2, 4, 5, 6], smi, "achiral", "not_applicable", "not_applicable",
                      assignment="No stereogenic element in the candidate graph."))

for mid, label, name, pages in [
    ("Ru_cod", "Cp*Ru(cod)Cl", "Cp*Ru(cod)Cl", [2]),
    ("Ru_naph", "[Cp*Ru(C10H8)]PF6", "[Cp*Ru(C10H8)]PF6", [1, 2, 3, 5, 6]),
    ("Et3NBnCl", "[Et3NBn]Cl", "benzyltriethylammonium chloride", [1, 2, 3, 5, 6]),
    ("Cu_tfaa", "Cu(F3-acac)2", "copper trifluoroacetylacetonate", [2, 4]),
    ("sodium_ascorbate", "sodium ascorbate", "sodium ascorbate", [5]),
    ("CuSO4_5H2O", "CuSO4·5H2O", "copper sulfate pentahydrate", [5]),
    ("petroleum_ether", "petroleum ether", "petroleum ether mixture", [4, 5, 6]),
    ("silica", "silica gel", "silica gel", [4, 5, 6]),
]:
    M.append(molecule(mid, label, name, pages, chirality="unresolved", geometry="not_applicable",
                      representation="unresolved", identity_note="Named source material; unique concrete molecular representation was not asserted."))

M.append(generic_molecule("amine_template", "substituted propargyl amine",
                          "HC#C-CH2-N(Me)-R",
                          [("R", "amine substituent", ["Me", "Bn", "alkynyl-containing chains shown in Scheme 1"])],
                          [2, 4], "Scheme 1A/general procedure uses a variable substituted amine."))
M.append(generic_molecule("diazo_template", "CF3-containing diazo compound 1",
                          "CF3-C(N2)-X",
                          [("X", "electron-withdrawing substituent", ["CO2Me", "P(O)(OEt)2", "other source instances unresolved"])],
                          [2, 4], "Scheme 1A/general procedure uses a variable diazo compound."))
M.append(generic_molecule("previous_hydroamination_product", "prior hydroamination product",
                          "Scheme 1B left product drawing", [("R1/R2", "amine substituents", None)], [2],
                          "Prior-work product is generic and is not a current-paper experiment."))
M.append(generic_molecule("previous_PK_product", "prior Pauson-Khand product",
                          "Scheme 1B right bicyclic ketone drawing", [("R", "amine substituent", None)], [2],
                          "Prior-work product is generic and is not a current-paper experiment."))
M.append(generic_molecule("expected_dimer_2e_g", "expected dimer of 2e-g",
                          "Scheme 3 crossed-out dimer", [("X/Y", "allene substituents", None)], [4],
                          "Expected product is drawn only as a generic target and was reported at 0%."))

MID = {m["molecule_id"].split(":mol:", 1)[1]: m["molecule_id"] for m in M}


def participant(pid, mid, role, pages, stage="reaction", mass=None, amount=None,
                volume=None, equivalents=None, mol_percent=None, reference=None,
                product_role=None):
    eids = [ev(p, "image") if p <= 6 else ev(p, "text") for p in pages]
    p = st.new_participant(pid, MID[mid], role, eids, stage=stage)
    if product_role:
        p["product_role"] = st.fact(product_role, product_role, eids)
    for key, value in (("mass", mass), ("amount", amount), ("volume", volume),
                       ("equivalents", equivalents), ("mol_percent", mol_percent)):
        if value is not None:
            raw = value.pop("_raw", None) if isinstance(value, dict) else str(value)
            p["charges"][key] = st.fact(raw or value, value, eids)
    if reference:
        p["equivalence_reference"] = st.fact(reference, MID[reference], eids)
    return p


def base_reaction(suffix, pages, record_type="specific_experiment", categories=None,
                  reaction_name=None, transformation=None):
    eids = [ev(p, "image") if p <= 6 else ev(p, "text") for p in pages]
    r = st.new_reaction(f"{DOC}:{suffix}", DOC, eids)
    r["record_type"] = record_type
    r["experiment_categories"] = st.fact(categories, categories, eids) if categories is not None else st.fact()
    r["location"] = st.fact({"pdf_pages": pages}, {"pdf_pages": pages}, eids)
    if reaction_name:
        r["reaction_name"] = st.fact(reaction_name, reaction_name, eids)
    if transformation:
        r["transformation_type"] = st.fact(None, transformation, eids, method="inference",
                                            derivation="Transformation label derived from the source scheme/text.")
    r["review"]["unresolved"] = []
    r["completeness"] = "partial"
    return r


def set_partial_classification(r, codes, pages, bonds_per_code=None, limitation=None):
    if not codes:
        r["classification"]["coverage"] = inferred("unresolved", pages, "Products or atom correspondence are insufficient for bond-change classification.")
        r["classification"]["limitations"] = [limitation or "No bond labels asserted without a secure product/atom correspondence."]
        return
    bonds_per_code = bonds_per_code or {code: 1 for code in codes}
    changes = []
    labels = []
    code_to_num = {row["code"]: row["number"] for row in st.read(ROOT / "skills/chemical-reaction-extraction/references/bond_taxonomy.json")["categories"]}
    for code in codes:
        support = []
        atoms = code.removeprefix("DB-").removeprefix("TB-").split("-")
        if code == "C-C": atoms = ["C", "C"]
        if code == "C-N": atoms = ["C", "N"]
        order = 2 if code.startswith("DB-") else 3 if code.startswith("TB-") else 1
        for _ in range(bonds_per_code.get(code, 1)):
            cid = f"b{len(changes)+1}"
            change = {"change_id": cid, "atoms": [{"atom_map": None, "element": atoms[0]}, {"atom_map": None, "element": atoms[1]}],
                      "before_order": 0, "after_order": order, "change_type": "formed"}
            changes.append(change); support.append(cid)
        labels.append({"category_number": code_to_num[code], "code": code, "supporting_change_ids": support})
    deriv = "Partial net-bond interpretation from the source transformation; atom mapping was not asserted."
    r["classification"]["bond_changes"] = [f(None, c, pages, method="structure_derivation", derivation=deriv) for c in changes]
    r["classification"]["labels"] = f(None, labels, pages, method="structure_derivation", derivation=deriv)
    r["classification"]["coverage"] = inferred("partial", pages, deriv)
    r["classification"]["limitations"] = [limitation or "No atom-mapped reaction SMILES; bond-order changes and broken bonds may be incomplete."]


def add_step(r, sid, stage, pids, pages, operation=None, temperature=None, time=None,
             atmosphere=None, pressure=None, method="original_extraction"):
    s = st.new_step(sid, [ev(p, "image") if p <= 6 else ev(p, "text") for p in pages], stage=stage)
    s["participant_ids"] = list(pids)
    if operation is not None: s["operation"] = f(operation, operation, pages, method=method)
    if temperature is not None: s["temperature"] = f(str(temperature), temperature, pages, method=method)
    if time is not None: s["time"] = f(str(time), time, pages, method=method)
    if atmosphere is not None: s["atmosphere"] = f(atmosphere, atmosphere, pages, method=method)
    if pressure is not None: s["pressure"] = f(pressure, pressure, pages, method=method)
    r["steps"].append(s)


def outcome(oid, kind, target, raw, normalized, pages, basis=None, yield_type=None, scope=None):
    o = st.new_outcome(oid, kind, [MID[x] for x in target], [ev(p, "image") if p <= 6 else ev(p, "text") for p in pages])
    o["value"] = f(raw, normalized, pages)
    if basis is not None: o["measurement_basis"] = f(basis, basis, pages)
    if yield_type is not None: o["yield_type"] = f(yield_type, yield_type, pages)
    if scope is not None: o["scope"] = f(scope, scope, pages)
    if target: o["product_role"] = f("target product", "target product", pages)
    return o


R = []

# Scheme 1 arrows: prior work is kept separate from current-paper generic scope.
r = base_reaction("scheme1A_prior_sigmatropic", [2], "literature_example", None,
                  "[2,3]-sigmatropic rearrangement", "copper-catalyzed ylide rearrangement to functionalized allene")
r["participants"] = [participant("amine", "amine_template", "reactant", [2]),
                     participant("diazo", "diazo_template", "reactant", [2]),
                     participant("cu", "Cu_tfaa", "catalyst", [2], mol_percent={"number":5,"unit":"mol%"}, reference="diazo_template")]
r["notes"] = [f("Scheme 1A is explicitly labeled Previous work and cites ref. 13.", "Scheme 1A is prior work, not a current experiment.", [2])]
set_partial_classification(r, ["C-C", "C-N"], [2], limitation="Generic R/X structures and no atom mapping; only net newly connected atom types are asserted.")
R.append(r)

for suffix, name, product, transform in [
    ("scheme1B_prior_hydroamination", "Cu-catalyzed hydroamination", "previous_hydroamination_product", "allene hydroamination"),
    ("scheme1B_prior_pauson_khand", "Co-catalyzed Pauson-Khand reaction", "previous_PK_product", "Pauson-Khand cyclization"),
]:
    r = base_reaction(suffix, [2], "literature_example", None, name, transform)
    r["participants"] = [participant("allene", "2e", "reactant", [2]), participant("product", product, "product", [2], product_role="reported prior-work product")]
    r["notes"] = [f("Scheme 1B is explicitly labeled Previous work and cites prior references.", "Prior-work illustration only.", [2])]
    set_partial_classification(r, [], [2], limitation="Generic prior-work drawing was not atom mapped or promoted to a current-paper experiment.")
    R.append(r)

r = base_reaction("scheme1C_general_dimerization", [2,3,5], "generic_scheme", None,
                  "Ruthenium-catalyzed dimerization of functional allenes", "tail-to-tail [2+2] allene dimerization")
r["participants"] = [participant("allene", "2a", "reactant", [2]), participant("product", "3a", "product", [2], product_role="generic 1,3-dimethylenecyclobutane")]
r["notes"] = [f("Scheme 1C is the current-work generic transformation; specific Table 1/Schemes 2 and 4 experiments are separate records.",
                  "Generic scheme only.", [2])]
set_partial_classification(r, ["C-C"], [2], {"C-C": 2}, "Two new C-C sigma connections close the cyclobutane ring; atom mapping and full bond-order changes are unresolved.")
R.append(r)

# Table 1 optimization: every row is independent.
entries = [
    (1, "Ru_cod", 3, "DCE", 85, 78), (2, "Ru_cod", 3, "toluene", 90, 82),
    (3, "Ru_cod", 5, "toluene", 95, 98), (4, "dual", 3, "DCE", 85, 61),
    (5, "dual", 5, "toluene", 90, 80), (6, "dual", 5, "toluene", 100, 100),
    (7, None, None, "toluene", 110, 0),
]
for n, cat, loading, solvent, temp, yld in entries:
    r = base_reaction(f"table1_entry_{n}", [2], "specific_experiment",
                      ["optimization"] + (["control"] if n == 7 else []),
                      "Optimization of allene dimerization", "tail-to-tail [2+2] allene dimerization")
    r["general_scheme_ids"] = [f"{DOC}:scheme1C_general_dimerization"]
    r["participants"] = [participant("substrate", "2a", "reactant", [2], amount={"number":0.2,"unit":"mmol","_raw":"0.2 mmol"})]
    if cat == "Ru_cod":
        r["participants"].append(participant("catalyst", "Ru_cod", "catalyst", [2], mol_percent={"number":loading,"unit":"mol%","_raw":f"{loading} mol%"}, reference="2a"))
    elif cat == "dual":
        r["participants"].extend([
            participant("Ru", "Ru_naph", "catalyst", [2], mol_percent={"number":loading,"unit":"mol%","_raw":f"{loading} mol% catalytic system"}, reference="2a"),
            participant("chloride", "Et3NBnCl", "additive", [2], mol_percent={"number":loading,"unit":"mol%","_raw":f"{loading} mol% catalytic system"}, reference="2a")])
    r["participants"].append(participant("solvent", solvent, "solvent", [2], volume={"number":3,"unit":"mL","_raw":"3 mL"}))
    if n != 7:
        r["participants"].append(participant("product", "3a", "product", [2], product_role="target product"))
    r["reaction_scale"] = f("allene 2a (0.2 mmol)", {"reference_molecule_id": MID["2a"], "amount": qty(0.2,"mmol")}, [2])
    r["solvent_conditions"]["composition"] = f(solvent, [solvent], [2])
    r["solvent_conditions"]["total_volume"] = f("3 mL", qty(3,"mL"), [2])
    duration = qty(2,"h") if n in (1,6) else None
    add_step(r, "reaction", "reaction", [p["participant_id"] for p in r["participants"] if p["stage"]=="reaction"], [2],
             "heated; monitored by 19F NMR spectroscopy and TLC", qty(temp,"°C"), duration)
    if n == 7:
        r["outcomes"] = [outcome("no_product", "no_reaction", [], "0%; no product was observed", {"number":0,"unit":"%"}, [2],
                                 basis="19F NMR spectroscopy and TLC analysis", scope={"target":"3a"})]
        set_partial_classification(r, [], [2], "Negative control: no observed product, so no expected-product bond label was inferred.")
    else:
        r["outcomes"] = [outcome("nmr_yield", "yield", ["3a"], f"{yld}%", qty(yld,"%"), [2], basis="19F NMR spectroscopy", yield_type="NMR yield")]
        if n == 6:
            r["outcomes"].append(outcome("isolated_yield", "yield", ["3a"], "71% isolated", qty(71,"%"), [2], basis="isolated", yield_type="isolated yield"))
        set_partial_classification(r, ["C-C"], [2], {"C-C":2}, "Two C-C ring-closing connections inferred; full atom mapping unavailable.")
    if duration is None:
        r["review"]["unresolved"].append("Reaction time is not stated for this table row; it was not inherited from another row.")
    R.append(r)


def dimerization_record(substrate, product, page, yield_value, ratio, label="scheme2", basis="isolated after column chromatography"):
    r = base_reaction(f"{label}_{substrate}_to_{product}", [page,6] if page != 6 else [6], "specific_experiment",
                      ["substrate_scope"], "Ruthenium-catalyzed allene dimerization", "tail-to-tail [2+2] allene dimerization")
    r["general_scheme_ids"] = [f"{DOC}:scheme1C_general_dimerization"]
    r["participants"] = [
        participant("substrate", substrate, "reactant", [page,6], mass={"number":0.25,"unit":"g","_raw":"0.25 g"}, amount={"number":11.2,"unit":"mmol","_raw":"11.2 mmol (as printed)"}),
        participant("Ru", "Ru_naph", "catalyst", [page,6], mass={"number":28.5,"unit":"mg","_raw":"28.5 mg"}, mol_percent={"number":5,"unit":"mol%","_raw":"5 mol% catalytic system"}, reference=substrate),
        participant("chloride", "Et3NBnCl", "additive", [page,6], mass={"number":12.7,"unit":"mg","_raw":"12.7 mg"}, mol_percent={"number":5,"unit":"mol%","_raw":"5 mol% catalytic system"}, reference=substrate),
        participant("solvent", "toluene", "solvent", [page,6], volume={"number":5,"unit":"mL","_raw":"5 mL"}),
        participant("product", product, "product", [page,6], product_role="target product mixture")]
    r["reaction_scale"] = f("corresponding allene (0.25 g, 11.2 mmol)", {"reference_molecule_id": MID[substrate], "mass":qty(0.25,"g"), "amount":qty(11.2,"mmol")}, [6])
    r["solvent_conditions"]["composition"] = f("anhydrous toluene", ["toluene"], [6])
    r["solvent_conditions"]["total_volume"] = f("5 mL", qty(5,"mL"), [6])
    add_step(r, "reaction", "reaction", ["substrate","Ru","chloride","solvent","product"], [page,6],
             "stirred under heating; completion monitored by TLC", qty(100,"°C"), qty(2,"h"), method="inherited_general_condition")
    add_step(r, "workup", "workup", [], [6], "solvent removed under reduced pressure")
    add_step(r, "purification", "purification", [], [6], "column chromatography on silica gel using petroleum ether/ethyl acetate")
    r["procedure"]["general_procedure_reference"] = f("Section 4.6, general procedure for allenes dimerization", "Section 4.6", [6])
    r["outcomes"] = [outcome("yield", "yield", [product], f"{yield_value}%", qty(yield_value,"%"), [page,6], basis=basis, yield_type="reported yield"),
                     outcome("selectivity", "selectivity", [product], ratio, ratio, [page], basis="1H and 19F NMR spectroscopy", scope={"metric":"Z/E ratio","order":["Z","E"]})]
    r["notes"] = [
        f("General procedure prints '0.25 g, 11.2 mmol' for the corresponding allene and '(28.5 mg/12.7 mg, 0.056 mmol)' for the catalytic system.",
          "Printed mass/amount/loading values are internally suspicious and were preserved without correction.", [6]),
        f("Catalyst-system amount 0.056 mmol is not unambiguously assigned to one or both components.",
          "Catalyst-system amount assignment unresolved.", [6])]
    r["review"]["unresolved"] += ["The printed 0.25 g/11.2 mmol allene charge and 0.056 mmol catalyst-system statement are not reconciled.",
                                     "The 90–120 °C general-procedure range differs from the 100 °C scope scheme; 100 °C is retained for these scheme-specific records."]
    set_partial_classification(r, ["C-C"], [page], {"C-C":2}, "Two C-C connections close the cyclobutane; no atom-mapped reaction SMILES asserted.")
    return r


for substrate, product, yld, ratio in [("2a","3a",71,"Z/E = 1:1"),("2b","3b",85,"Z/E = 1:1"),
                                       ("2c","3c",55,"Z/E = 2:3"),("2d","3d",86,"Z/E = 1:1")]:
    R.append(dimerization_record(substrate, product, 3, yld, ratio))

# Scheme 3 failed dimerizations are actual attempted substrate cases, not a
# fabricated no-reaction product.  The source reports complex unidentified mixtures.
for substrate in ("2e","2f","2g"):
    r = base_reaction(f"scheme3_failed_dimerization_{substrate}", [2,4], "specific_experiment",
                      ["substrate_scope"], "Attempted Ru-catalyzed dimerization of allenyne", "unidentified-product formation")
    r["general_scheme_ids"] = [f"{DOC}:scheme1C_general_dimerization"]
    r["participants"] = [participant("substrate", substrate, "reactant", [2,4]), participant("solvent","toluene","solvent",[4])]
    r["solvent_conditions"]["composition"] = f("toluene", ["toluene"], [4])
    add_step(r, "reaction", "reaction", ["substrate","solvent"], [4], "heated with Ru catalytic conditions", {"minimum":90,"maximum":110,"unit":"°C","raw_text":"90–110 °C"})
    r["outcomes"] = [outcome("unidentified_products", "other", [], "complicated mixture of products; expected dimer 0%", "complicated unidentified mixture", [2,4], scope={"expected_product":"crossed-out dimer in Scheme 3"})]
    r["notes"] = [f("Authors attribute the failure tentatively to coordination of the additional triple bond.", "Author interpretation, not a mechanistic experiment.", [2])]
    r["review"]["unresolved"] = ["Individual mixture components, product yields, reaction time, catalyst loading, and measurement basis are unreported."]
    set_partial_classification(r, [], [4], "Products were not identified; no bond-formation class can be assigned safely.")
    R.append(r)

# Click reactions 2e-g + two azides -> 4a-f.
click_map = [("2e","pivoxymethyl_azide","4a",87,"4/1"),("2e","benzyl_azide","4b",98,"2/1"),
             ("2f","pivoxymethyl_azide","4c",88,"4/1"),("2f","benzyl_azide","4d",84,"2/1"),
             ("2g","pivoxymethyl_azide","4e",50,"3/1"),("2g","benzyl_azide","4f",60,"2/1")]
for substrate, azide, product, yld, eluent in click_map:
    r = base_reaction(f"scheme3_click_{substrate}_to_{product}", [4,5,6], "specific_experiment",
                      ["substrate_scope"], "Cu(I)-catalyzed azide-alkyne ligation", "1,4-disubstituted 1,2,3-triazole formation")
    r["participants"] = [
        participant("substrate", substrate, "reactant", [4,5], mass={"number":0.14,"unit":"g","_raw":"0.14 g"}, amount={"number":0.54,"unit":"mmol","_raw":"0.54 mmol"}),
        participant("azide", azide, "reactant", [4,5], mass={"number":0.08,"unit":"g","_raw":"0.08 g"}, amount={"number":0.54,"unit":"mmol","_raw":"0.54 mmol"}),
        participant("ascorbate", "sodium_ascorbate", "reagent", [5], mass={"number":21.5,"unit":"mg","_raw":"21.5 mg"}, amount={"number":0.1,"unit":"mmol","_raw":"0.1 mmol"}),
        participant("CuSO4", "CuSO4_5H2O", "catalyst", [5], amount={"number":0.03,"unit":"mmol","_raw":"0.5 M copper sulfate pentahydrate (0.03 mmol)"}),
        participant("tBuOH", "tBuOH", "solvent", [5]), participant("water", "water", "solvent", [5]),
        participant("product", product, "product", [4,5,6], product_role="target triazoloallene"),
        participant("EtOAc", "ethyl_acetate", "workup_agent", [5], stage="workup", volume={"number":15,"unit":"mL","_raw":"15 mL"}),
        participant("HCl", "HCl", "workup_agent", [5], stage="workup"),
        participant("wash_water", "water", "workup_agent", [5], stage="workup", volume={"number":15,"unit":"mL","_raw":"15 mL"}),
        participant("drying", "MgSO4", "workup_agent", [5], stage="workup"),
        participant("PE", "petroleum_ether", "purification_agent", [5], stage="purification"),
        participant("EA", "ethyl_acetate", "purification_agent", [5], stage="purification")]
    r["reaction_scale"] = f("acetylene-containing allene (0.14 g, 0.54 mmol)", {"reference_molecule_id":MID[substrate],"mass":qty(0.14,"g"),"amount":qty(0.54,"mmol")}, [5])
    r["solvent_conditions"]["composition"] = f("tBuOH/H2O (1:1)", ["tert-butanol","water"], [5])
    r["solvent_conditions"]["mixture_ratio"] = f("1:1", "1:1", [5])
    r["solvent_conditions"]["total_volume"] = f("5 mL", qty(5,"mL"), [5])
    add_step(r, "addition", "reaction", ["substrate","azide","ascorbate","CuSO4","tBuOH","water"], [5],
             "azide, sodium ascorbate, and copper sulfate pentahydrate added sequentially")
    add_step(r, "reaction", "reaction", ["substrate","azide","ascorbate","CuSO4","tBuOH","water","product"], [5],
             "stirred; completion monitored by TLC", "room temperature", qty(7,"h"))
    add_step(r, "evaporation", "workup", ["EtOAc"], [5], "mixed solvent removed under reduced pressure; residue taken up in ethyl acetate")
    add_step(r, "acid_wash", "workup", ["HCl"], [5], "washed with 0.5 M HCl (2 × 10 mL)")
    add_step(r, "water_wash_dry", "workup", ["wash_water","drying"], [5], "ethyl acetate extracts washed with water, dried over MgSO4, filtered, and evaporated")
    add_step(r, "chromatography", "purification", ["PE","EA"], [5,6], f"column chromatography, petroleum ether/ethyl acetate = {eluent}")
    r["procedure"]["original_text"] = f("Section 4.4 general procedure for the synthesis of triazole-containing allenes", "Section 4.4", [5])
    r["procedure"]["addition_order"] = f("benzyl or pivaloylmethyl azide, sodium ascorbate, then CuSO4·5H2O sequentially", ["azide","sodium ascorbate","CuSO4·5H2O"], [5])
    r["outcomes"] = [outcome("yield", "yield", [product], f"Yield {yld}%", qty(yld,"%"), [4,5,6], basis="isolated after column chromatography", yield_type="isolated yield")]
    r["review"]["unresolved"] = ["The printed universal 0.14 g/0.54 mmol and 0.08 g/0.54 mmol charges are retained for each instance; substrate-specific mass consistency was not corrected.",
                                     "Volume of the 0.5 M CuSO4·5H2O solution is not reported."]
    set_partial_classification(r, ["C-N"], [4], {"C-N":2}, "CuAAC ring formation inferred to create two C-N connections; full azide resonance/bond-order mapping is unresolved.")
    R.append(r)

# Scheme 4 dimerizations 4a-f -> 5a-f.  Scheme yields are retained; 5b has
# a later extracted-text conflict (76% in Scheme 4 vs 73% in section 4.7).
scheme5_yields = {"5a":78,"5b":76,"5c":76,"5d":70,"5e":78,"5f":70}
for substrate, product in zip(names4, ratios5):
    r = dimerization_record(substrate, product, 5, scheme5_yields[product], ratios5[product], label="scheme4", basis=None)
    r["outcomes"][0]["measurement_basis"] = st.fact()  # Scheme 4 does not state the basis.
    r["outcomes"][0]["value"]["evidence_ids"] = [ev(5,"image")]
    r["outcomes"][0]["yield_type"] = st.fact()
    if product == "5b":
        alt_scheme = f("76% in Scheme 4", qty(76,"%"), [5])
        alt_text = f("Yield 73% in section 4.7 extracted text", qty(73,"%"), [7])
        conflict = f(["76%", "73%"], None, [5,7], status="conflict", alternatives=[alt_scheme, alt_text])
        r["outcomes"][0]["value"] = conflict
        r["review"]["unresolved"].append("5b yield conflict: Scheme 4 page 5 shows 76%; page-7 experimental extracted text reports 73%. Page 7 was context only, not part of visual coverage.")
    r["review"]["unresolved"].append("Yield measurement basis is not stated in Scheme 4; it was not guessed as isolated/NMR.")
    R.append(r)

# General [2,3]-sigmatropic preparation applied to 2a-g.  Only 2b/2c have
# explicit yields on the covered pages.
prep_yields = {"2b":42,"2c":75}
prep_eluent = {"2b":"10/1","2c":"5/1"}
for product in allene_data:
    pages = [4,5] if product in prep_yields else [4,5]
    r = base_reaction(f"preparation_{product}", pages, "specific_experiment", ["other"],
                      "[2,3]-sigmatropic rearrangement", "copper-catalyzed ylide rearrangement to functionalized allene")
    r["participants"] = [
        participant("amine", "amine_template", "reactant", [4], mass={"number":210.0,"unit":"mg","_raw":"210.0 mg"}, amount={"number":11.3,"unit":"mmol","_raw":"11.3 mmol (as printed)"}),
        participant("diazo", "diazo_template", "reactant", [4], mass={"number":94,"unit":"mg","_raw":"94 mg"}, amount={"number":11.3,"unit":"mmol","_raw":"11.3 mmol (as printed)"}),
        participant("catalyst", "Cu_tfaa", "catalyst", [4], mass={"number":20.9,"unit":"mg","_raw":"20.9 mg"}, amount={"number":0.056,"unit":"mmol","_raw":"0.056 mmol"}),
        participant("solvent", "toluene", "solvent", [4], volume={"number":5,"unit":"mL","_raw":"5 mL"}),
        participant("product", product, "product", pages, product_role="functionalized allene")]
    r["reaction_scale"] = f("corresponding diazocompound (11.3 mmol as printed)", {"reference_molecule_id":MID["diazo_template"],"amount":qty(11.3,"mmol")}, [4])
    r["solvent_conditions"]["composition"] = f("anhydrous toluene", ["toluene"], [4])
    r["solvent_conditions"]["total_volume"] = f("5 mL", qty(5,"mL"), [4])
    add_step(r, "reaction", "reaction", ["amine","diazo","catalyst","solvent","product"], [4],
             "stirred under heating; completion monitored by TLC", {"minimum":90,"maximum":120,"unit":"°C","raw_text":"90–120 °C"},
             {"minimum":1,"maximum":2,"unit":"h","raw_text":"1–2 h"})
    add_step(r, "evaporation", "workup", [], [4,5], "solvent removed under reduced pressure")
    op = "column chromatography on silica gel with petroleum ether/ethyl acetate"
    if product in prep_eluent: op += f" = {prep_eluent[product]}"
    add_step(r, "chromatography", "purification", [], [4,5], op)
    r["procedure"]["general_procedure_reference"] = f("Section 4.2, general procedure for [2,3]-sigmatropic rearrangement", "Section 4.2", [4])
    if product in prep_yields:
        r["outcomes"] = [outcome("yield", "yield", [product], f"Yield {prep_yields[product]}%", qty(prep_yields[product],"%"), [5], basis="isolated after column chromatography", yield_type="isolated yield")]
    r["notes"] = [f("General procedure prints amine 210.0 mg/11.3 mmol, diazocompound 94 mg/11.3 mmol, and catalyst 20.9 mg/0.056 mmol.",
                    "Printed mass/amount values preserved without correction.", [4])]
    r["review"]["unresolved"] = ["Exact substituted amine and diazocompound instance structures are not specified in the covered procedure and remain generic templates.",
                                     "Printed 11.3 mmol charges appear inconsistent with the 0.056 mmol catalyst and 5 mol% Scheme 1 loading; no decimal correction was made."]
    if product not in prep_yields:
        r["review"]["unresolved"].append("No product-specific yield is reported on pages 1–6; the cited prior preparation was not imported as a yield.")
    set_partial_classification(r, ["C-C","C-N"], [2,4], limitation="Net C-C and C-N connections inferred from Scheme 1A; exact atom mapping and broken N-C bond remain unresolved.")
    R.append(r)

# Source context: characterization/observations are not synthetic experiments.
context = [
    {"context_id":f"{DOC}:context:noesy_3d", "kind":"characterization", "evidence_ids":[ev(3,"image")],
     "description":"Fig. 1, 2D 1H NOESY of cis/trans 3d; authors assign more intense olefin/cyclobutane cross-peaks to trans isomer."},
    {"context_id":f"{DOC}:context:xrd_3a_3c", "kind":"characterization", "evidence_ids":[ev(4,"image")],
     "description":"Fig. 2 X-ray structures: (S,R)-E-3a CCDC 2081066, (S,R)-E-3b CCDC 2081065, and 3c CCDC 2081067."},
    {"context_id":f"{DOC}:context:characterization_pages_5_6", "kind":"characterization", "evidence_ids":[ev(5,"image"),ev(6,"image")],
     "description":"NMR and elemental-analysis blocks for 2b, 2c, 4a-f, and 3a-d; preserved as context and not counted as separate syntheses."},
    {"context_id":f"{DOC}:context:general_analytical_methods", "kind":"general_information", "evidence_ids":[ev(4,"image")],
     "description":"General solvent drying, TLC, chromatography, and NMR instrumentation information."},
]

bundle = st.empty(DATASET["documents"])
bundle["molecules"] = M
bundle["reactions"] = R
bundle["source_context"] = context
st.write(RUN / "submitted.json", bundle)

inventory = {
    "document_id": DOC,
    "pages": [1,2,3,4,5,6],
    "molecule_count": len(M),
    "reaction_count": len(R),
    "record_ids": [r["reaction_id"] for r in R],
    "source_context_ids": [c["context_id"] for c in context],
    "notes": ["Machine extraction only; needs_review after merge.",
              "Page 7 text is used only for the 5b yield conflict and is not claimed as visually covered."]
}
st.write(RUN / "raw_model/pages_0001_0006_inventory.json", inventory)
print(json.dumps({"molecules":len(M),"reactions":len(R)}, ensure_ascii=False))
