#!/usr/bin/env python3
"""Build the source-specific v3.1 bundle for Chem. Commun. 2021, 57, 3488.

The script is intentionally paper-specific.  It consumes only the clean template
made by ``standard init`` and facts transcribed from PDF pages 1--4.
"""

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DOC = "doc_dbfa3548c0baadcd"
EV = {
    p: [f"{DOC}:page{p}:text", f"{DOC}:page{p}:image"] for p in range(1, 5)
}


def f(raw=None, normalized=None, evidence=(), method="original_extraction", derivation=None,
      status=None):
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation, status=status)


def q(number, unit, raw=None, qualifier=None):
    d = {"number": number, "unit": unit}
    if qualifier:
        d["qualifier"] = qualifier
    if raw is not None:
        d["raw_text"] = raw
    return d


def missing():
    return f()


def set_name(m, name, ev):
    m["identity_status"] = "source_identified"
    m["names"]["original"] = f(name, name, ev)
    m["names"]["english"] = f(name, name, ev)
    m["names"]["abbreviations"] = f()
    m["names"]["aliases"] = f()
    m["cas_numbers"] = f()
    m["external_properties"] = []


STEREO_CHECKS = {}


def smiles_molecule(mid, label, name, smiles, ev, *, geometry="not_applicable",
                    chirality="achiral", representation="encoded", assignment=None,
                    note=None):
    m = st.new_molecule(mid, DOC, ev)
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    set_name(m, name, ev)
    m["representation_type"] = "smiles"
    m["smiles"] = f(
        f"source drawing/name for {label}", smiles, ev, method="structure_derivation",
        derivation="Candidate SMILES transcribed conservatively from the source drawing or explicit chemical name; parsing does not certify identity."
    )
    diag = st.stereo_inspect(smiles)
    STEREO_CHECKS[mid] = diag
    stc = m["stereochemistry"]
    stc["atom_centers"] = f("candidate graph inspection", diag["atom_centers"], ev,
                              method="structure_derivation",
                              derivation="RDKit graph diagnostic on the candidate SMILES; source assignment takes precedence.")
    stc["double_bonds"] = f("candidate graph inspection", diag["double_bonds"], ev,
                              method="structure_derivation",
                              derivation="RDKit graph diagnostic; absence of a source E/Z label remains unspecified.")
    stc["axial"] = f()
    stc["mixture"] = f()
    if assignment is None and chirality == "achiral":
        assignment = "No atom-centred stereogenic element is present in the asserted candidate graph; this does not assign any unencoded conformational chirality."
    stc["assignment"] = f(
        "source drawing plus candidate-graph inspection", assignment, ev,
        method="structure_derivation",
        derivation="Stereochemical status is limited to what can be established from the visible drawing and candidate graph."
    ) if assignment else f()
    stc["chirality_status"] = f(
        "source drawing/candidate graph", chirality, ev, method="structure_derivation",
        derivation="Status does not infer racemic composition from absent stereomarks."
    )
    stc["geometry_status"] = f(
        "source drawing/candidate graph", geometry, ev, method="structure_derivation",
        derivation="Visible alkene geometry is left unspecified unless the source explicitly assigns it."
    )
    stc["representation_status"] = f(
        "candidate graph encoding", representation, ev, method="structure_derivation",
        derivation="Candidate encoding is partial when source stereochemistry is unassigned or cumulene handling is limited."
    )
    stc["isomeric_smiles"] = f(
        "candidate graph", diag["canonical_isomeric_smiles"], ev,
        method="structure_derivation", derivation="Canonical isomeric SMILES calculated from the candidate graph."
    )
    stc["identity_key"] = f(
        "candidate graph", diag["stereo_graph_key"], ev,
        method="structure_derivation", derivation="Stereo graph key calculated from the candidate graph; not a sample identity."
    )
    m["identity_note"] = f(note, note, ev) if note else f()
    return m


def unresolved_molecule(mid, label, name, ev, *, chirality="unspecified",
                        geometry="unspecified", assignment=None, note=None):
    m = st.new_molecule(mid, DOC, ev)
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    set_name(m, name, ev)
    m["representation_type"] = "unresolved"
    m["smiles"] = f()
    stc = m["stereochemistry"]
    stc["atom_centers"] = f()
    stc["double_bonds"] = f()
    stc["axial"] = f()
    stc["mixture"] = f()
    stc["assignment"] = f(assignment, assignment, ev) if assignment else f()
    stc["chirality_status"] = f(
        "no explicit absolute configuration in the source drawing", chirality, ev,
        method="structure_derivation",
        derivation="Missing stereomarks are retained as unspecified, not converted to achiral or racemic."
    )
    stc["geometry_status"] = f(
        "no explicit E/Z assignment in the source drawing", geometry, ev,
        method="structure_derivation",
        derivation="The compact source drawing was not converted into an unsupported geometric assignment."
    )
    stc["representation_status"] = f(
        "complete stereo-safe graph not asserted", "not_encoded", ev,
        method="structure_derivation",
        derivation="The source label/drawing is retained without forcing an unsafe complete SMILES."
    )
    stc["isomeric_smiles"] = f()
    stc["identity_key"] = f()
    m["identity_note"] = f(note, note, ev) if note else f(
        "Compact source structure retained by label; complete primary structure requires expert review.",
        "Compact source structure retained by label; complete primary structure requires expert review.", ev)
    return m


def template_molecule(mid, label, name, ev, fixed, variables, reason):
    m = st.new_molecule(mid, DOC, ev)
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    set_name(m, name, ev)
    m["representation_type"] = "generic_template"
    m["smiles"] = f()
    m["generic_template"] = {
        "fixed_structure": f(fixed, fixed, ev),
        "variables": [{
            "label": v,
            "attachment_points": f(f"source-labelled attachment for {v}", f"source-labelled attachment for {v}", ev),
            "definition": f(f"variable {v}", f"variable {v}", ev),
            "allowed_values": f(),
            "constraints": f(),
        } for v in variables],
        "instance_mapping": [],
        "reason": f(reason, reason, ev),
    }
    stc = m["stereochemistry"]
    for k, val in (("chirality_status", "unresolved"), ("geometry_status", "unresolved"),
                   ("representation_status", "not_encoded")):
        stc[k] = f("generic R-group drawing", val, ev, method="structure_derivation",
                   derivation="A generic family cannot be assigned one concrete stereochemical identity.")
    m["identity_note"] = f(reason, reason, ev)
    return m


def participant(pid, mid, role, ev, *, stage="reaction", amount=None, equiv=None,
                molpct=None, mass=None, volume=None, reference=None, product_role=None):
    p = st.new_participant(pid, mid, role, ev, stage)
    p["product_role"] = f(product_role, product_role, ev) if product_role else f()
    if amount:
        p["charges"]["amount"] = f(amount[0], amount[1], ev)
    if equiv:
        p["charges"]["equivalents"] = f(equiv[0], equiv[1], ev)
    if molpct:
        p["charges"]["mol_percent"] = f(molpct[0], molpct[1], ev)
    if mass:
        p["charges"]["mass"] = f(mass[0], mass[1], ev)
    if volume:
        p["charges"]["volume"] = f(volume[0], volume[1], ev)
    if reference:
        p["equivalence_reference"] = f(reference, reference, ev)
    return p


def reaction_base(suffix, ev, record_type, categories, location, name, transform):
    r = st.new_reaction(f"{DOC}:{suffix}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, ev)
    r["location"] = f(location, location, ev)
    r["reaction_name"] = f(name, name, ev)
    r["transformation_type"] = f(transform, transform, ev)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    return r


def add_reaction_step(r, ev, participant_ids, temperature="rt", hours=12,
                      operation="combine and react under the displayed conditions"):
    s = st.new_step("reaction", ev)
    s["participant_ids"] = participant_ids
    s["operation"] = f(operation, operation, ev)
    s["temperature"] = f("rt", "rt", ev) if temperature == "rt" else f(
        f"{temperature} °C", q(temperature, "°C", f"{temperature} °C"), ev)
    s["time"] = f(f"{hours} h", q(hours, "h", f"{hours} h"), ev)
    r["steps"].append(s)


def set_solvent(r, ev, mid, name, volume=None):
    r["solvent_conditions"]["composition"] = f(name, [mid], ev)
    if volume:
        r["solvent_conditions"]["total_volume"] = f(volume[0], volume[1], ev)


def outcome(r, oid, kind, target, raw, normalized, ev, *, basis=None, scope=None,
            yield_type=None, product_role="reported product"):
    o = st.new_outcome(oid, kind, [target] if target else [], ev)
    o["product_role"] = f(product_role, product_role, ev) if target else f()
    o["value"] = f(raw, normalized, ev)
    o["yield_type"] = f(yield_type, yield_type, ev) if yield_type else f()
    o["measurement_basis"] = f(basis, basis, ev) if basis else f()
    o["scope"] = f(scope, scope, ev) if scope else f()
    r["outcomes"].append(o)


def set_classification(r, ev, formed_cc=0, coverage="partial", limitation=None):
    changes = []
    for i in range(formed_cc):
        changes.append({
            "change_id": f"b{i+1}",
            "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
            "before_order": 0, "after_order": 1, "change_type": "formed",
        })
    r["classification"]["bond_changes"] = [
        f(None, c, ev, method="inference", derivation="Class-level bond change inferred from the displayed transformation; exact atom mapping is unresolved.")
        for c in changes
    ]
    labels = []
    if changes:
        labels = [{"category_number": 2, "code": "C-C",
                   "supporting_change_ids": [c["change_id"] for c in changes]}]
    r["classification"]["labels"] = f(
        None, labels, ev, method="inference",
        derivation="Only newly formed class-level C-C connections supported by the annulation/addition drawing are indexed."
    )
    r["classification"]["coverage"] = f(
        None, coverage, ev, method="inference",
        derivation="Full atom mapping, hydrogen transfers, leaving-group fate, and coproducts are unresolved."
    )
    r["classification"]["limitations"] = [limitation or
        "Exact atoms, all bond-order changes, hydrogen changes, coproducts, and full atom mapping remain unresolved."]
    r["classification"]["model_votes"] = []


def set_na_classification(r, ev, reason):
    r["classification"]["bond_changes"] = []
    r["classification"]["labels"] = f(None, [], ev, method="inference", derivation=reason)
    r["classification"]["coverage"] = f(None, "not_applicable", ev, method="inference", derivation=reason)
    r["classification"]["limitations"] = [reason]
    r["classification"]["model_votes"] = []


bundle = json.loads((RUN / "raw_model/standard_init.json").read_text())

# Molecule registry ---------------------------------------------------------
molecules = []
molecules += [
    template_molecule(f"{DOC}:prior_dienolate_family", "A", "β,γ-unsaturated carbonyl dienolate family", EV[1],
                      "R1-substituted conjugated dienolate bearing carbonyl R2", ["R1", "R2", "E"],
                      "Scheme 1a is a generic literature reaction-mode diagram."),
    template_molecule(f"{DOC}:prior_addition_product_family", "prior products", "prior γ/α-addition product families", EV[1],
                      "generic γ- or α-functionalized β,γ-unsaturated carbonyl scaffold", ["R1", "R2", "E"],
                      "The source groups cited precedents and does not define a single concrete molecule."),
    template_molecule(f"{DOC}:generic_substrate1", "1", "β,γ-unsaturated carbonyl substrate family", EV[1] + EV[2],
                      "R1-CO-CH2-CH=CH-R2", ["R1", "R2"], "Scheme 1b/Scheme 2 use variable substituents."),
    template_molecule(f"{DOC}:generic_product3", "3", "cyclopentene annulation product family", EV[1] + EV[2],
                      "substituted cyclopentene bearing carbonyl, alkenyl, and CO2Bn substituents", ["R1", "R2"],
                      "The generic product contains variable substituents and unassigned stereochemistry."),
]

substrates = {
    "1a": ("1,4-diphenylbut-3-en-1-one", "O=C(c1ccccc1)CC=Cc1ccccc1"),
    "1b": ("1-(4-methylphenyl)-4-phenylbut-3-en-1-one", "Cc1ccc(C(=O)CC=Cc2ccccc2)cc1"),
    "1c": ("1-(4-methoxyphenyl)-4-phenylbut-3-en-1-one", "COc1ccc(C(=O)CC=Cc2ccccc2)cc1"),
    "1d": ("1-(4-fluorophenyl)-4-phenylbut-3-en-1-one", "Fc1ccc(C(=O)CC=Cc2ccccc2)cc1"),
    "1e": ("1-(4-chlorophenyl)-4-phenylbut-3-en-1-one", "Clc1ccc(C(=O)CC=Cc2ccccc2)cc1"),
    "1f": ("1-(4-bromophenyl)-4-phenylbut-3-en-1-one", "Brc1ccc(C(=O)CC=Cc2ccccc2)cc1"),
    "1g": ("1-(2-bromophenyl)-4-phenylbut-3-en-1-one", "O=C(c1ccccc1Br)CC=Cc1ccccc1"),
    "1h": ("1-(3-methylphenyl)-4-phenylbut-3-en-1-one", "Cc1cccc(C(=O)CC=Cc2ccccc2)c1"),
    "1i": ("1-(3-methoxyphenyl)-4-phenylbut-3-en-1-one", "COc1cccc(C(=O)CC=Cc2ccccc2)c1"),
    "1j": ("1-(furan-2-yl)-4-phenylbut-3-en-1-one", "O=C(c1ccco1)CC=Cc1ccccc1"),
    "1k": ("1-(thiophen-2-yl)-4-phenylbut-3-en-1-one", "O=C(c1ccsc1)CC=Cc1ccccc1"),
    "1l": ("1-(naphthalen-1-yl)-4-phenylbut-3-en-1-one", "O=C(c1cccc2ccccc12)CC=Cc1ccccc1"),
    "1m": ("4-(4-methylphenyl)-1-phenylbut-3-en-1-one", "O=C(c1ccccc1)CC=Cc1ccc(C)cc1"),
    "1n": ("4-(4-methoxyphenyl)-1-phenylbut-3-en-1-one", "O=C(c1ccccc1)CC=Cc1ccc(OC)cc1"),
    "1o": ("4-(4-fluorophenyl)-1-phenylbut-3-en-1-one", "O=C(c1ccccc1)CC=Cc1ccc(F)cc1"),
    "1p": ("1-phenylbut-3-en-1-one", "O=C(c1ccccc1)CC=C"),
    "1q": ("4-(naphthalen-2-yl)-1-phenylbut-3-en-1-one", "O=C(c1ccccc1)CC=Cc1ccc2ccccc2c1"),
    "1r": ("5-phenylpent-4-en-2-one", "CC(=O)CC=Cc1ccccc1"),
    "1s": ("1-phenylhexa-1,5-dien-3-one", "C=CCC(=O)C=Cc1ccccc1"),
    "1t": ("β,γ-unsaturated N-acyl-3,5-dimethylpyrazole 1t", "Cc1nn(C(=O)CC=Cc2ccccc2)c(C)c1"),
}
for label, (name, smi) in substrates.items():
    ev = EV[2] + (EV[3] if label in {"1a", "1g", "1t"} else [])
    molecules.append(smiles_molecule(
        f"{DOC}:{label}", label, name, smi, ev, geometry="unspecified", representation="partial",
        note="Candidate connectivity follows the visible scope drawing; alkene E/Z is not explicitly assigned on pages 1–4."
    ))

molecules.append(smiles_molecule(
    f"{DOC}:2", "2", "benzyl β′-acetoxy allenoate 2",
    "C=C=C(COC(C)=O)C(=O)OCc1ccccc1", EV[1] + EV[2] + EV[3],
    geometry="unspecified", chirality="achiral", representation="partial",
    assignment="The visible allene has a terminal CH2 end, so that local allene axis is not stereogenic; cumulene/toolkit representation remains partial.",
    note="β′-Acetoxy allenoate drawing transcribed as a candidate; no axial descriptor is reported."
))

for label in [f"3{x}" for x in list("abcdefghijklmnopqrs")]:
    ev = EV[2] + (EV[3] if label in {"3a", "3g"} else [])
    molecules.append(unresolved_molecule(
        f"{DOC}:{label}", label, f"cyclopentene annulation product {label}", ev,
        chirality="unspecified", geometry="unspecified",
        note="The source gives a compact cyclopentene drawing but no absolute configuration, dr, or E/Z assignment; a complete stereo-safe SMILES was not asserted."
    ))

molecules += [
    unresolved_molecule(f"{DOC}:4", "4", "umpolung γ-addition product 4", EV[2] + EV[3]),
    unresolved_molecule(f"{DOC}:5", "5", "alternative annulation product 5", EV[3] + EV[4],
                        chirality="unspecified", geometry="unspecified",
                        note="The displayed product contains an allene/cumulene motif without a source axial assignment; complete encoding was withheld."),
    unresolved_molecule(f"{DOC}:6", "6", "spiro[4.4]nonane derivative 6", EV[3]),
]

simple = {
    "PPh3": ("triphenylphosphine", "P(c1ccccc1)(c1ccccc1)c1ccccc1"),
    "PPh2Me": ("methyldiphenylphosphine", "CP(c1ccccc1)c1ccccc1"),
    "PPhMe2": ("dimethylphenylphosphine", "CP(C)c1ccccc1"),
    "PnBu3": ("tri-n-butylphosphine", "CCCCP(CCCC)CCCC"),
    "Panisyl3": ("tris(4-methoxyphenyl)phosphine", "COc1ccc(P(c2ccc(OC)cc2)c2ccc(OC)cc2)cc1"),
    "Ptolyl3": ("tris(4-methylphenyl)phosphine", "Cc1ccc(P(c2ccc(C)cc2)c2ccc(C)cc2)cc1"),
    "Pfluoro3": ("tris(4-fluorophenyl)phosphine", "Fc1ccc(P(c2ccc(F)cc2)c2ccc(F)cc2)cc1"),
    "PtBu3": ("tri-tert-butylphosphine", "CC(C)(C)P(C(C)(C)C)C(C)(C)C"),
    "Na2CO3": ("sodium carbonate", "[Na+].[Na+].[O-]C([O-])=O"),
    "Cs2CO3": ("caesium carbonate", "[Cs+].[Cs+].[O-]C([O-])=O"),
    "NEt3": ("triethylamine", "CCN(CC)CC"),
    "NaOH": ("sodium hydroxide", "[Na+].[OH-]"),
    "toluene": ("toluene", "Cc1ccccc1"),
    "THF": ("tetrahydrofuran", "C1CCOC1"),
    "PdOAc2": ("palladium(II) acetate", "CC(=O)[O-].CC(=O)[O-].[Pd+2]"),
    "K2CO3": ("potassium carbonate", "[K+].[K+].[O-]C([O-])=O"),
    "MeCN": ("acetonitrile", "CC#N"),
}
for key, (name, smi) in simple.items():
    page = 3 if key in {"PtBu3", "NaOH", "PdOAc2", "K2CO3", "MeCN"} else 2
    if key == "PPh3":
        ev = EV[2] + EV[3]
    else:
        ev = EV[page]
    molecules.append(smiles_molecule(f"{DOC}:{key}", key, name, smi, ev))

bundle["molecules"] = molecules

# Generic and cited reaction records ---------------------------------------
reactions = []
for suffix, name, transform in [
    ("literature_gamma_addition", "cited vinylogous γ-addition of β,γ-unsaturated carbonyl compounds", "literature γ-addition family"),
    ("literature_alpha_addition", "cited α-addition/twofold-addition concept for β,γ-unsaturated carbonyl compounds", "literature α-addition family"),
]:
    r = reaction_base(suffix, EV[1], "literature_example", [], "Page 1, Scheme 1a", name, transform)
    r["participants"] = [
        participant("dienolate", f"{DOC}:prior_dienolate_family", "reactant", EV[1]),
        participant("prior_product", f"{DOC}:prior_addition_product_family", "product", EV[1], product_role="generic cited product family"),
    ]
    r["review"]["unresolved"] = ["Scheme 1a groups cited precedents; no single concrete experiment, conditions, or yield is specified."]
    set_classification(r, EV[1], formed_cc=1, limitation="Generic cited family only; exact structures and atom mapping are unavailable.")
    reactions.append(r)

r = reaction_base("generic_scheme", EV[1] + EV[2], "generic_scheme", [],
                  "Pages 1–2, Scheme 1b and Scheme 2 header",
                  "phosphine-catalysed (4+1) annulation", "generic (4+1) annulation")
r["participants"] = [
    participant("generic_1", f"{DOC}:generic_substrate1", "reactant", EV[1] + EV[2]),
    participant("allenoate_2", f"{DOC}:2", "reactant", EV[1] + EV[2]),
    participant("generic_3", f"{DOC}:generic_product3", "product", EV[1] + EV[2], product_role="annulation product family"),
]
r["review"]["unresolved"] = ["Generic R-group scheme; individual experiments are represented separately."]
set_classification(r, EV[1] + EV[2], formed_cc=2,
                   limitation="Generic family supports two class-level C-C formations; exact atom mapping and hydrogen/leaving-group changes are unresolved.")
reactions.append(r)
GENERIC = r["reaction_id"]

# Optimization table -------------------------------------------------------
opt_rows = [
    ("table1_entry1", "1", "PPh3", "Na2CO3", "toluene", 30),
    ("table1_entry2", "2", "PPh2Me", "Na2CO3", "toluene", 55),
    ("table1_entry3", "3", "PPhMe2", "Na2CO3", "toluene", 41),
    ("table1_entry4_nbu", "4 (first occurrence)", "PnBu3", "Na2CO3", "toluene", 28),
    ("table1_entry4_panisyl", "4 (second occurrence; prose calls this entry 5)", "Panisyl3", "Na2CO3", "toluene", 71),
    ("table1_entry5", "5", "Ptolyl3", "Na2CO3", "toluene", 60),
    ("table1_entry6", "6", "Pfluoro3", "Na2CO3", "toluene", 22),
    ("table1_entry7", "7", "Panisyl3", "Cs2CO3", "toluene", 34),
    ("table1_entry8", "8", "Panisyl3", "NEt3", "toluene", 56),
    ("table1_entry9", "9", "Panisyl3", "Na2CO3", "THF", 47),
]
for suffix, printed, catalyst, base, solvent, yld in opt_rows:
    ev = EV[2]
    is_3a = suffix == "table1_entry4_panisyl"
    location = f"Page 2, Table 1, printed entry {printed}"
    if is_3a:
        location += "; Page 2, Scheme 2, product 3a"
    r = reaction_base(suffix, ev, "specific_experiment",
                      ["optimization", "substrate_scope"] if is_3a else ["optimization"],
                      location, "phosphine-catalysed (4+1) annulation", "(4+1) annulation")
    r["general_scheme_ids"] = [GENERIC]
    r["participants"] = [
        participant("substrate_1a", f"{DOC}:1a", "reactant", ev,
                    amount=("0.24 mmol", q(0.24, "mmol", "0.24 mmol")),
                    equiv=("1.2 equiv", q(1.2, "equiv", "1.2 equiv")), reference=f"{DOC}:2"),
        participant("allenoate_2", f"{DOC}:2", "reactant", ev,
                    amount=("0.2 mmol", q(0.2, "mmol", "0.2 mmol"))),
        participant("phosphine", f"{DOC}:{catalyst}", "catalyst", ev,
                    amount=("0.04 mmol", q(0.04, "mmol", "0.04 mmol")),
                    molpct=("20 mol%", q(20, "mol%", "20 mol%")), reference=f"{DOC}:2"),
        participant("base", f"{DOC}:{base}", "additive", ev,
                    amount=("0.24 mmol", q(0.24, "mmol", "0.24 mmol")),
                    equiv=("1.2 equiv", q(1.2, "equiv", "1.2 equiv")), reference=f"{DOC}:2"),
        participant("solvent", f"{DOC}:{solvent}", "solvent", ev,
                    volume=("3 mL", q(3, "mL", "3 mL"))),
        participant("product_3a", f"{DOC}:3a", "product", ev, product_role="desired annulation product"),
    ]
    r["reaction_scale"] = f("allenoate 2: 0.2 mmol", {"reference_molecule_id": f"{DOC}:2", "amount": q(0.2, "mmol")}, ev)
    set_solvent(r, ev, f"{DOC}:{solvent}", solvent, ("3 mL", q(3, "mL", "3 mL")))
    add_reaction_step(r, ev, [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"])
    r["procedure"]["original_text"] = f("Table 1 footnote a conditions", "Table 1 footnote a conditions", ev)
    outcome(r, "yield_3a", "yield", f"{DOC}:3a", f"{yld}%", q(yld, "%", f"{yld}%"), ev,
            basis="isolated yield", yield_type="isolated yield")
    r["review"]["unresolved"] = ["Workup and purification are not reported on pages 1–4; ESI was not supplied."]
    if is_3a:
        r["notes"].append(f(
            "Table 1 prints the anisylphosphine optimum as a second entry 4, while adjacent prose refers to it as entry 5; Scheme 2 product 3a is the same 71% standard-condition experiment.",
            "Table 1 prints the anisylphosphine optimum as a second entry 4, while adjacent prose refers to it as entry 5; Scheme 2 product 3a is the same 71% standard-condition experiment.", ev))
        r["review"]["unresolved"].append("Printed optimization-entry numbering is internally inconsistent (duplicate entry 4 versus prose entry 5).")
    set_classification(r, ev, formed_cc=2)
    reactions.append(r)

# Scheme 2 scope (3a is already merged into the optimized row) ------------
scope_yields = {
    "3b": 77, "3c": 81, "3d": 58, "3e": 58, "3f": 62, "3g": 50,
    "3h": 63, "3i": 74, "3j": 54, "3k": 65, "3l": 61, "3m": 75,
    "3n": 74, "3o": 45, "3p": 58, "3q": 61, "3r": 38, "3s": 66,
}
for product, yld in scope_yields.items():
    substrate = "1" + product[1:]
    ev = EV[2]
    r = reaction_base(f"scheme2_{product}", ev, "specific_experiment", ["substrate_scope"],
                      f"Page 2, Scheme 2, product {product}",
                      "phosphine-catalysed (4+1) annulation", "(4+1) annulation")
    r["general_scheme_ids"] = [GENERIC]
    r["participants"] = [
        participant("substrate", f"{DOC}:{substrate}", "reactant", ev,
                    equiv=("1.2 equiv", q(1.2, "equiv", "1.2 equiv")), reference=f"{DOC}:2"),
        participant("allenoate_2", f"{DOC}:2", "reactant", ev),
        participant("catalyst", f"{DOC}:Panisyl3", "catalyst", ev,
                    molpct=("20 mol%", q(20, "mol%", "20 mol%")), reference=f"{DOC}:2"),
        participant("base", f"{DOC}:Na2CO3", "additive", ev,
                    equiv=("1.2 equiv", q(1.2, "equiv", "1.2 equiv")), reference=f"{DOC}:2"),
        participant("solvent", f"{DOC}:toluene", "solvent", ev),
        participant("product", f"{DOC}:{product}", "product", ev, product_role="desired annulation product"),
    ]
    set_solvent(r, ev, f"{DOC}:toluene", "toluene")
    add_reaction_step(r, ev, [p["participant_id"] for p in r["participants"]])
    r["procedure"]["general_procedure_reference"] = f("Scheme 2 standard conditions", "Scheme 2 standard conditions", ev)
    outcome(r, f"yield_{product}", "yield", f"{DOC}:{product}", f"{yld}%", q(yld, "%", f"{yld}%"), ev)
    r["review"]["unresolved"] = [
        "Scheme 2 does not state the yield measurement basis.",
        "Absolute scale, solvent volume, atmosphere, workup, and purification are not reported on pages 1–4."
    ]
    set_classification(r, ev, formed_cc=2)
    reactions.append(r)

# Scheme 3, equation 1, Scheme 5, and footnote control ---------------------
ev = EV[2] + EV[3]
r = reaction_base("scheme3_product4", ev, "specific_experiment", ["substrate_scope"],
                  "Page 3, Scheme 3", "phosphine-catalysed reaction of N-acylpyrazole 1t with allenoate 2",
                  "umpolung γ-addition")
r["participants"] = [
    participant("substrate_1t", f"{DOC}:1t", "reactant", ev,
                equiv=("1.2 equiv", q(1.2, "equiv", "1.2 equiv")), reference=f"{DOC}:2"),
    participant("allenoate_2", f"{DOC}:2", "reactant", ev),
    participant("catalyst", f"{DOC}:Panisyl3", "catalyst", ev,
                molpct=("20 mol%", q(20, "mol%", "20 mol%")), reference=f"{DOC}:2"),
    participant("base", f"{DOC}:Na2CO3", "additive", ev,
                equiv=("1.2 equiv", q(1.2, "equiv", "1.2 equiv")), reference=f"{DOC}:2"),
    participant("solvent", f"{DOC}:toluene", "solvent", ev),
    participant("product_4", f"{DOC}:4", "product", ev, product_role="reported γ-addition product"),
]
set_solvent(r, ev, f"{DOC}:toluene", "toluene")
add_reaction_step(r, ev, [p["participant_id"] for p in r["participants"]])
outcome(r, "yield_4", "yield", f"{DOC}:4", "30%", q(30, "%", "30%"), ev)
r["review"]["unresolved"] = ["Yield measurement basis, scale, solvent volume, atmosphere, workup, and purification are unreported."]
set_classification(r, ev, formed_cc=1)
reactions.append(r)

for base, yld in [("Na2CO3", 27), ("NaOH", 69)]:
    suffix = "eq1_ptbu3_na2co3" if base == "Na2CO3" else "eq1_ptbu3_naoh"
    r = reaction_base(suffix, ev, "specific_experiment", ["control", "mechanistic_study"],
                      f"Page 3, equation (1), base = {base}",
                      "PtBu3-catalysed alternative annulation", "regiodivergent annulation to product 5")
    r["participants"] = [
        participant("substrate_1a", f"{DOC}:1a", "reactant", ev,
                    equiv=("1.2 equiv", q(1.2, "equiv", "1.2 equiv")), reference=f"{DOC}:2"),
        participant("allenoate_2", f"{DOC}:2", "reactant", ev),
        participant("catalyst", f"{DOC}:PtBu3", "catalyst", ev,
                    molpct=("20 mol%", q(20, "mol%", "20 mol%")), reference=f"{DOC}:2"),
        participant("base", f"{DOC}:{base}", "additive", ev,
                    equiv=("1.2 equiv", q(1.2, "equiv", "1.2 equiv")), reference=f"{DOC}:2"),
        participant("solvent", f"{DOC}:toluene", "solvent", ev),
        participant("product_5", f"{DOC}:5", "product", ev, product_role="reported alternative annulation product"),
    ]
    set_solvent(r, ev, f"{DOC}:toluene", "toluene")
    add_reaction_step(r, ev, [p["participant_id"] for p in r["participants"]])
    outcome(r, "yield_5", "yield", f"{DOC}:5", f"{yld}%", q(yld, "%", f"{yld}%"), ev)
    r["review"]["unresolved"] = ["Yield measurement basis, scale, solvent volume, atmosphere, workup, and purification are unreported."]
    set_classification(r, ev, formed_cc=1)
    reactions.append(r)

r = reaction_base("scheme5_scaleup_3g", ev, "specific_experiment", ["scale_up"],
                  "Page 3, Scheme 5, upper arrow", "scale-up (4+1) annulation to 3g", "(4+1) annulation")
r["general_scheme_ids"] = [GENERIC]
r["participants"] = [
    participant("substrate_1g", f"{DOC}:1g", "reactant", ev,
                equiv=("1.2 equiv", q(1.2, "equiv", "1.2 equiv")), reference=f"{DOC}:2"),
    participant("allenoate_2", f"{DOC}:2", "reactant", ev,
                amount=("1.0 mmol", q(1.0, "mmol", "1.0 mmol")),
                mass=("0.25 g", q(0.25, "g", "0.25 g"))),
    participant("catalyst", f"{DOC}:Panisyl3", "catalyst", ev,
                molpct=("20 mol%", q(20, "mol%", "20 mol%")), reference=f"{DOC}:2"),
    participant("base", f"{DOC}:Na2CO3", "additive", ev,
                equiv=("1.2 equiv", q(1.2, "equiv", "1.2 equiv")), reference=f"{DOC}:2"),
    participant("solvent", f"{DOC}:toluene", "solvent", ev),
    participant("product_3g", f"{DOC}:3g", "product", ev, product_role="desired annulation product"),
]
r["reaction_scale"] = f("allenoate 2: 1.0 mmol", {"reference_molecule_id": f"{DOC}:2", "amount": q(1.0, "mmol")}, ev)
set_solvent(r, ev, f"{DOC}:toluene", "toluene")
add_reaction_step(r, ev, [p["participant_id"] for p in r["participants"]])
outcome(r, "amount_3g", "product_amount", f"{DOC}:3g", "0.28 g", q(0.28, "g", "0.28 g"), ev)
outcome(r, "yield_3g", "yield", f"{DOC}:3g", "57%", q(57, "%", "57%"), ev)
r["review"]["unresolved"] = ["Yield measurement basis, solvent volume, atmosphere, workup, and purification are unreported."]
set_classification(r, ev, formed_cc=2)
reactions.append(r)

r = reaction_base("scheme5_heck_to6", ev, "specific_experiment", ["derivatization"],
                  "Page 3, Scheme 5, lower arrow", "intramolecular Heck reaction of 3g", "Heck cyclization")
r["participants"] = [
    participant("substrate_3g", f"{DOC}:3g", "reactant", ev),
    participant("palladium", f"{DOC}:PdOAc2", "catalyst", ev,
                molpct=("10 mol%", q(10, "mol%", "10 mol%")), reference=f"{DOC}:3g"),
    participant("phosphine", f"{DOC}:PPh3", "ligand", ev,
                molpct=("40 mol%", q(40, "mol%", "40 mol%")), reference=f"{DOC}:3g"),
    participant("base", f"{DOC}:K2CO3", "additive", ev,
                equiv=("2.0 equiv", q(2.0, "equiv", "2.0 equiv")), reference=f"{DOC}:3g"),
    participant("solvent", f"{DOC}:MeCN", "solvent", ev),
    participant("product_6", f"{DOC}:6", "product", ev, product_role="spirocyclic Heck product"),
]
set_solvent(r, ev, f"{DOC}:MeCN", "MeCN")
add_reaction_step(r, ev, [p["participant_id"] for p in r["participants"]], temperature=85, hours=6)
outcome(r, "yield_6", "yield", f"{DOC}:6", "73%", q(73, "%", "73%"), ev)
r["review"]["unresolved"] = ["Yield measurement basis, scale, solvent volume, atmosphere, workup, and purification are unreported."]
set_classification(r, ev, formed_cc=1, limitation="One class-level intramolecular C-C formation is inferred from the Heck cyclization; exact atom mapping and other changes remain unresolved.")
reactions.append(r)

r = reaction_base("scheme4_proposed_mechanism", ev, "proposed_mechanism", [],
                  "Page 3, Scheme 4", "plausible mechanism for (4+1) annulation and product 5 formation",
                  "author-proposed mechanism")
r["participants"] = [
    participant("substrate_1a", f"{DOC}:1a", "reactant", ev),
    participant("allenoate_2", f"{DOC}:2", "reactant", ev),
    participant("product_3a", f"{DOC}:3a", "product", ev, product_role="path-a product"),
    participant("product_5", f"{DOC}:5", "product", ev, product_role="path-b product"),
]
r["notes"] = [f("Intermediates A–D, B, C and C′ are proposed and were not expanded into executed experiments.",
                 "Intermediates A–D, B, C and C′ are proposed and were not expanded into executed experiments.", ev)]
r["review"]["unresolved"] = ["Proposed intermediates and elementary steps are not independently verified experiments."]
set_classification(r, ev, formed_cc=0, coverage="unresolved", limitation="Proposed mechanism only; no executed reaction record or complete atom mapping is asserted.")
reactions.append(r)

ev = EV[4] + EV[3]
r = reaction_base("no_phosphine_product5_control", ev, "specific_experiment", ["control", "mechanistic_study"],
                  "Page 4, footnote 12; discussed on page 3", "reaction of 1a and 2 without phosphine catalyst",
                  "negative control")
r["participants"] = [
    participant("substrate_1a", f"{DOC}:1a", "reactant", ev),
    participant("allenoate_2", f"{DOC}:2", "reactant", ev),
    participant("candidate_product_5", f"{DOC}:5", "product", ev, product_role="product explicitly not detected"),
]
s = st.new_step("reaction", ev)
s["participant_ids"] = ["substrate_1a", "allenoate_2"]
s["operation"] = f("conducted in the absence of phosphine catalyst", "conducted in the absence of phosphine catalyst", ev)
r["steps"] = [s]
outcome(r, "product5_not_detected", "no_reaction", f"{DOC}:5", "no product 5 could be detected", True, ev,
        product_role="candidate product not detected")
r["review"]["unresolved"] = [
    "Footnote 12 does not identify which base variant, scale, temperature, time, atmosphere, or detection method was used; these were not inherited by guesswork."
]
set_na_classification(r, ev, "Negative control: no product 5 was detected, so the intended transformation is not assigned formation labels.")
reactions.append(r)

bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": f"{DOC}:metadata", "document_id": DOC, "kind": "article_metadata",
     "evidence_ids": EV[1], "value": {"title": "Phosphine-catalysed (4+1) annulations of β′-acetoxy allenoate with β,γ-unsaturated carbonyl compounds",
                                               "doi": "10.1039/d1cc00368b"}},
    {"context_id": f"{DOC}:scheme1_context", "document_id": DOC, "kind": "conceptual_reaction_modes",
     "evidence_ids": EV[1], "value": "Scheme 1 contrasts cited γ/α-addition modes with the current design; proposed intermediates are not concrete experiments."},
    {"context_id": f"{DOC}:characterization_inventory", "document_id": DOC, "kind": "negative_characterization_inventory",
     "evidence_ids": EV[1] + EV[2] + EV[3] + EV[4],
     "value": "No standalone mass-spectrum peak table, molecular-weight table, or fragment-relative-abundance block occurs on pages 1–4."},
    {"context_id": f"{DOC}:external_materials", "document_id": DOC, "kind": "external_material_notice",
     "evidence_ids": EV[1], "value": "The article advertises ESI, but no SI file was supplied or associated in this workspace."},
]

(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/stereo_checks.json").write_text(json.dumps(STEREO_CHECKS, ensure_ascii=False, indent=2) + "\n")

inventory = {
    "document_id": DOC,
    "source_pages": [1, 2, 3, 4],
    "visual_inspection": "all four original rendered page images inspected at original detail",
    "pages": [
        {"page": 1, "regions": ["title/abstract/ESI notice", "Scheme 1a cited reaction modes", "Scheme 1b current design"],
         "decisions": ["two generic literature-example families", "one current generic scheme", "no concrete current experiment"]},
        {"page": 2, "regions": ["Table 1, ten printed data rows", "Scheme 2, products 3a–3s", "scope prose"],
         "decisions": ["ten optimization experiments", "3a merged with optimized anisylphosphine row", "eighteen additional scope experiments", "duplicate printed entry 4 retained as source numbering issue"]},
        {"page": 3, "regions": ["Scheme 3", "equation (1)", "Scheme 4", "Scheme 5", "discussion"],
         "decisions": ["one N-acylpyrazole experiment", "two base-specific product-5 experiments", "one proposed mechanism", "one scale-up", "one Heck derivatization"]},
        {"page": 4, "regions": ["references 1–12", "footnote 12"],
         "decisions": ["one no-phosphine negative control from footnote 12", "bibliography not imported as current experiments"]},
    ],
    "counts": {"reactions": len(reactions), "molecules": len(molecules)},
}
(RUN / "raw_model/page_extraction.json").write_text(json.dumps(inventory, ensure_ascii=False, indent=2) + "\n")
print(json.dumps(inventory["counts"], ensure_ascii=False))
