#!/usr/bin/env python3
"""Build the source-bounded v3.1 submission for Chem. Commun. 2021, 57, 2364-2367."""

from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DOC = "doc_923e87cf958fa06a"
P = {i: f"{DOC}:page{i}:text" for i in range(1, 5)}
I = {i: f"{DOC}:page{i}:image" for i in range(1, 5)}
MID = lambda label: f"{DOC}:mol:{label}"
RID = lambda label: f"{DOC}:{label}"


def f(raw=None, normalized=None, ev=(), method="original_extraction", derivation=None,
      status=None, alternatives=()):
    return st.fact(raw, normalized, ev, method=method, derivation=derivation,
                   status=status, alternatives=alternatives)


def q(number, unit, qualifier=None, raw_text=None):
    value = {"number": number, "unit": unit}
    if qualifier:
        value["qualifier"] = qualifier
    if raw_text:
        value["raw_text"] = raw_text
    return value


def inferred(raw, normalized, ev, reason):
    return f(raw, normalized, ev, method="inference", derivation=reason)


bundle = json.loads((RUN / "submitted.json").read_text())
molecules = {}


def add_molecule(label, name, ev, smiles=None, chirality="unspecified",
                 geometry="unspecified", representation=None, assignment=None,
                 mixture=None, axial=None, note=None):
    mid = MID(label)
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    m["names"]["original"] = f(name, name, ev)
    m["names"]["english"] = f(name, name, ev)
    m["names"]["abbreviations"] = f()
    m["names"]["aliases"] = f()
    m["cas_numbers"] = f()
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = f(
            f"source drawing/name for {label}", smiles, ev,
            method="structure_derivation",
            derivation="Candidate SMILES transcribed conservatively from the supplied drawing/name; graph parsing does not verify identity.",
        )
        rep = representation or "partial"
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = f()
        rep = representation or "not_encoded"
    m["generic_template"] = None
    m["stereochemistry"]["atom_centers"] = f()
    m["stereochemistry"]["double_bonds"] = f()
    m["stereochemistry"]["axial"] = f(axial, axial, ev) if axial is not None else f()
    m["stereochemistry"]["mixture"] = f(mixture, mixture, ev) if mixture is not None else f()
    m["stereochemistry"]["assignment"] = f(assignment, assignment, ev) if assignment is not None else f()
    m["stereochemistry"]["chirality_status"] = inferred(
        f"stereochemical review of {label}", chirality, ev,
        "Status records only what the supplied drawing/text specifies; missing stereomarks are not treated as achiral or racemic.",
    )
    m["stereochemistry"]["geometry_status"] = inferred(
        f"geometric-stereochemistry review of {label}", geometry, ev,
        "Status records the supplied drawing/text; unassigned alkene geometry remains unspecified.",
    )
    m["stereochemistry"]["representation_status"] = inferred(
        f"representation review of {label}", rep, ev,
        "Candidate graph coverage and source stereochemical detail were compared; unresolved source products are not forced into SMILES.",
    )
    m["stereochemistry"]["isomeric_smiles"] = f()
    m["stereochemistry"]["identity_key"] = f()
    m["identity_note"] = f(note, note, ev) if note else f()
    m["external_properties"] = []
    molecules[mid] = m
    return mid


def add_template(label, name, ev, fixed, variables, reason):
    mid = add_molecule(label, name, ev, representation="not_encoded")
    m = molecules[mid]
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": f(fixed, fixed, ev),
        "variables": [
            {
                "label": v,
                "attachment_points": f("source drawing", "source-defined attachment point", ev),
                "definition": f("R-group variable", "source-defined substituent variable", ev),
                "allowed_values": f(),
                "constraints": f(),
            }
            for v in variables
        ],
        "instance_mapping": [],
        "reason": f(reason, reason, ev),
    }
    return mid


def product_molecule(label, family, ev, assignment=None):
    if family == "4":
        mixture = "Diastereomeric product sample; experiment-specific dr is stored on the reaction."
        note = "Chiral alpha-butadienyl amine bearing the (1S,2R)-aminoindanol auxiliary; compact source drawing not encoded as a complete SMILES."
    elif family == "5":
        mixture = "Diastereomeric product sample; experiment-specific dr is stored on the reaction."
        note = "Chiral beta-butadienyl/isoprenyl amine bearing the aminoindanol auxiliary; compact source drawing not encoded as a complete SMILES."
    else:
        mixture = "Product stereoisomer composition is represented by reaction-specific dr where reported."
        note = "Complex downstream product retained from the source drawing without forcing a complete SMILES."
    return add_molecule(label, f"product {label}", ev, chirality="mixture",
                        geometry="unspecified", representation="not_encoded",
                        assignment=assignment, mixture=mixture, note=note)


# Generic/prior-work identities.
lit1_s = add_template("lit1_substrate", "chiral arylpropargyl amide family", (I[2],),
                      "Ar-CH(NHAc)-C#CH", ["Ar"], "Cited Scheme 1a uses an R/aryl family drawing.")
lit1_p = add_template("lit1_product", "chiral 2-(N-acetyl-1-arylmethyl)-1,3-butadiene family", (I[2],),
                      "Ar-CH(NHAc)-C(=CH2)-CH=CH2", ["Ar"], "Cited product is shown as a family.")
lit2_s = add_template("lit2_substrate", "N,O-aminol family", (I[2],),
                      "protected alpha-amino ester/iminium precursor", ["PG", "R", "Rprime"], "Cited Scheme 1b is generic.")
lit2_p = add_template("lit2_product", "enantiopure diene amino-acid family", (I[2],),
                      "protected amino ester bearing a 1,3-diene", ["PG", "R", "Rprime"], "Cited Scheme 1b is generic.")
gen_aldehyde = add_template("generic_aldehyde_1", "aldehyde family 1", (I[2], I[3]), "R-CHO", ["R"], "Schemes 1c-3 use an R-group aldehyde.")
gen4 = add_template("generic_product_4", "alpha-butadienyl amine family 4", (I[2], I[3]), "aminoindanol-NH-CH(R)-alpha-butadienyl", ["R"], "Current-work family drawing.")
gen5 = add_template("generic_product_5", "beta-butadienyl amine family 5", (I[2], I[3]), "aminoindanol-NH-CH(R)-CH2-C(Me)=CH2", ["R"], "Current-work family drawing.")

# Recurrent substrates, reagents, solvents, and aminols.
aldehydes = {
    "1a": ("benzaldehyde", "O=Cc1ccccc1"),
    "pMe": ("4-methylbenzaldehyde", "Cc1ccc(C=O)cc1"),
    "pOMe": ("4-methoxybenzaldehyde", "COc1ccc(C=O)cc1"),
    "pF": ("4-fluorobenzaldehyde", "O=Cc1ccc(F)cc1"),
    "pCl": ("4-chlorobenzaldehyde", "O=Cc1ccc(Cl)cc1"),
    "pBr": ("4-bromobenzaldehyde", "O=Cc1ccc(Br)cc1"),
    "pCO2Me": ("methyl 4-formylbenzoate", "COC(=O)c1ccc(C=O)cc1"),
    "pNO2": ("4-nitrobenzaldehyde", "O=Cc1ccc([N+](=O)[O-])cc1"),
    "pIPr": ("4-isopropylbenzaldehyde", "CC(C)c1ccc(C=O)cc1"),
    "pTBu": ("4-tert-butylbenzaldehyde", "CC(C)(C)c1ccc(C=O)cc1"),
    "pCN": ("4-formylbenzonitrile", "N#Cc1ccc(C=O)cc1"),
    "mMe": ("3-methylbenzaldehyde", "Cc1cccc(C=O)c1"),
    "mOMe": ("3-methoxybenzaldehyde", "COc1cccc(C=O)c1"),
    "mF": ("3-fluorobenzaldehyde", "O=Cc1cccc(F)c1"),
    "mCl": ("3-chlorobenzaldehyde", "O=Cc1cccc(Cl)c1"),
    "mBr": ("3-bromobenzaldehyde", "O=Cc1cccc(Br)c1"),
    "mCF3": ("3-(trifluoromethyl)benzaldehyde", "O=Cc1cccc(C(F)(F)F)c1"),
    "mCN": ("3-formylbenzonitrile", "N#Cc1cccc(C=O)c1"),
    "mCO2Me": ("methyl 3-formylbenzoate", "COC(=O)c1cccc(C=O)c1"),
    "oOMe": ("2-methoxybenzaldehyde", "COc1ccccc1C=O"),
    "oF": ("2-fluorobenzaldehyde", "O=Cc1ccccc1F"),
    "oCF3": ("2-(trifluoromethyl)benzaldehyde", "O=Cc1ccccc1C(F)(F)F"),
    "oNO2": ("2-nitrobenzaldehyde", "O=Cc1ccccc1[N+](=O)[O-]"),
    "dimethyl34": ("3,4-dimethylbenzaldehyde", "Cc1ccc(C=O)cc1C"),
    "dibromo35": ("3,5-dibromobenzaldehyde", "O=Cc1cc(Br)cc(Br)c1"),
    "bromo_nitro35": ("3-bromo-5-nitrobenzaldehyde", "O=Cc1cc(Br)cc([N+](=O)[O-])c1"),
    "naph2": ("naphthalene-2-carbaldehyde", "O=Cc1ccc2ccccc2c1"),
    "fur2": ("furan-2-carbaldehyde", "O=Cc1ccco1"),
    "thio2": ("thiophene-2-carbaldehyde", "O=Cc1cccs1"),
    "pyr2": ("pyridine-2-carbaldehyde", "O=Cc1ccccn1"),
    "butanal": ("butanal", "CCCC=O"),
    "isobutyr": ("2-methylpropanal", "CC(C)C=O"),
    "isoval": ("3-methylbutanal", "CC(C)CC=O"),
    "cycho": ("cyclohexanecarbaldehyde", "O=CC1CCCCC1"),
}
AM = {}
for key, (name, smi) in aldehydes.items():
    aldehyde_ev = (I[2], I[3]) if key == "1a" else (I[3],)
    AM[key] = add_molecule(f"aldehyde_{key}", name, aldehyde_ev, smi,
                           chirality="unspecified", geometry="not_applicable")

bor2a = add_molecule("2a", "pinacol homoallenylboronate 2a", (I[2], I[3]),
                     "C=C=CCB1OC(C)(C)C(C)(C)O1", chirality="unspecified",
                     geometry="unspecified", axial="Terminal-CH2 allene; no axial configuration is assigned on the depicted local structure.")
bor2b = add_molecule("2b", "pinacol isoprenylboronate 2b", (I[2], I[3]),
                     "C=C(C)CB1OC(C)(C)C(C)(C)O1", chirality="unspecified", geometry="unspecified")
aminols = {}
for lab, name in {
    "3a": "chiral aminodiol 3a", "3b": "chiral aminodiol 3b", "3c": "chiral amino alcohol 3c",
    "3d": "(1S,2R)-2-aminocyclohexanol 3d", "3e": "(1S,2R)-1-amino-2-indanol 3e",
    "3f": "indanamine 3f",
}.items():
    assignment = "Source identifies (1S,2R) configuration." if lab in ("3d", "3e") else "Source drawing specifies stereochemistry where shown; exact absolute descriptors were not fully transcribed."
    aminols[lab] = add_molecule(lab, name, (I[2],), chirality="specified" if lab in ("3d", "3e") else "partially_specified",
                                geometry="not_applicable", assignment=assignment,
                                note="Compact Table 1 drawing retained without forcing a complete SMILES.")

simple = {
    "MeOH": ("methanol", "CO"), "DMSO": ("dimethyl sulfoxide", "CS(C)=O"),
    "THF": ("tetrahydrofuran", "C1CCOC1"), "toluene": ("toluene", "Cc1ccccc1"),
    "CH2Cl2": ("dichloromethane", "ClCCl"), "EtOAc": ("ethyl acetate", "CCOC(C)=O"),
    "dimethyl_fumarate": ("dimethyl fumarate", "COC(=O)/C=C/C(=O)OC"),
    "diethyl_acetylenedicarboxylate": ("diethyl acetylenedicarboxylate", "CCOC(=O)C#CC(=O)OCC"),
    "MeI": ("methyl iodide", "CI"), "Ag2O": ("silver(I) oxide", "[O-2].[Ag+].[Ag+]"),
    "MnO2": ("manganese dioxide", "O=[Mn]=O"), "anisole": ("anisole", "COc1ccccc1"),
    "dioxane": ("1,4-dioxane", "O1CCOCC1"),
}
simple_ev = {
    "MeOH": (I[2], I[3]), "DMSO": (I[2], I[3]), "THF": (I[2],),
    "toluene": (I[2], I[4]), "CH2Cl2": (I[2], I[4]), "EtOAc": (I[2],),
    "dimethyl_fumarate": (I[4],), "diethyl_acetylenedicarboxylate": (I[4],),
    "MeI": (I[4],), "Ag2O": (I[4],), "MnO2": (I[4],), "anisole": (I[4],),
    "dioxane": (I[4],),
}
SM = {k: add_molecule(k, name, simple_ev[k], smi, chirality="unspecified", geometry="unspecified")
      for k, (name, smi) in simple.items()}
ms4a = add_molecule("4A_MS", "4 angstrom molecular sieves", (I[2],),
                    note="Non-unique material; no concrete molecular structure asserted.")

# Products.
prod4 = {}
for lab in [f"4{c}" for c in "abcdefghijklmnopqrst"]:
    prod4[lab] = product_molecule(lab, "4", (I[2], I[3]))
table_products = {
    "3a": product_molecule("table_product_from_3a", "4", (I[2],)),
    "3b": product_molecule("table_product_from_3b", "4", (I[2],)),
    "3c": product_molecule("table_product_from_3c", "4", (I[2],)),
}
prod5 = {}
for lab in [f"5{c}" for c in "abcdefghijklmnopqrstuvwxy"]:
    assignment = None
    if lab == "5u":
        assignment = "X-ray crystallography (CCDC 2023204) assigns the newly formed stereogenic carbon of the major product as S."
    elif lab[1] <= "w":
        assignment = "For aryl products, the authors assign the newly formed center by analogy to X-ray-characterized 5u; exact source limits are retained in the reaction note."
    prod5[lab] = product_molecule(lab, "5", (I[3],), assignment=assignment)
prod6 = product_molecule("6", "downstream", (I[4],))
prod7 = product_molecule("7", "downstream", (I[4],))
prod8 = product_molecule("8", "downstream", (I[4],))
prod9 = product_molecule("9", "downstream", (I[4],), assignment="Source describes 9 as optically pure; no numerical ee or complete absolute configuration is provided on page 4.")


def participant(pid, mid, role, ev, stage="reaction", mass=None, amount=None,
                volume=None, equivalents=None, mol_percent=None, ref=None):
    p = st.new_participant(pid, mid, role, ev, stage=stage)
    for key, val in (("mass", mass), ("amount", amount), ("volume", volume),
                     ("equivalents", equivalents), ("mol_percent", mol_percent)):
        if val is not None:
            raw = val.get("raw_text") or f"{val.get('number')} {val.get('unit')}"
            p["charges"][key] = f(raw, val, ev)
    if ref:
        p["equivalence_reference"] = f(ref, ref, ev)
    return p


def classification(r, labels=None, changes=None, coverage="partial", limitation=None):
    labels = labels or []
    changes = changes or []
    r["classification"]["mapped_reaction_smiles"] = f()
    r["classification"]["bond_changes"] = [
        inferred("class-level source transformation", ch, r["evidence_ids"],
                 "Class-level bond change inferred from the source reaction drawing; atom identities and complete mapping are unresolved.")
        for ch in changes
    ]
    r["classification"]["labels"] = inferred(
        "class-level source transformation", labels, r["evidence_ids"],
        "Labels include all and only the asserted formed/promoted class-level bonds; no atom mapping is claimed.",
    )
    r["classification"]["coverage"] = inferred(
        "classification review", coverage, r["evidence_ids"],
        "Coverage is partial unless the experiment reports no observable transformation.",
    )
    r["classification"]["limitations"] = [limitation or "Exact atom identities, full atom mapping, hydrogen transfers, and coproducts remain unresolved."]
    r["classification"]["model_votes"] = []


def reaction(label, ev, record_type="specific_experiment", categories=(), location=None,
             name=None, transformation=None, general=()):
    r = st.new_reaction(RID(label), DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = f(list(categories), list(categories), ev)
    r["location"] = f(location, location, ev)
    r["reaction_name"] = f(name, name, ev) if name else f()
    r["transformation_type"] = inferred(
        transformation, transformation, ev,
        "Transformation description is a class-level interpretation of the source drawing/text.",
    ) if transformation else f()
    r["general_scheme_ids"] = list(general)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    return r


def add_step(r, sid, ev, pids, operation, temp=None, time=None, stage="reaction"):
    s = st.new_step(sid, ev, stage=stage)
    s["participant_ids"] = list(pids)
    s["operation"] = f(operation, operation, ev)
    if temp is not None:
        s["temperature"] = f(temp if isinstance(temp, str) else f"{temp['number']} {temp['unit']}", temp, ev)
    if time is not None:
        s["time"] = f(time if isinstance(time, str) else f"{time['number']} {time['unit']}", time, ev)
    r["steps"].append(s)


def add_outcome(r, oid, kind, target, value, ev, basis=None, scope=None, yield_type=None):
    o = st.new_outcome(oid, kind, [target] if target else [], ev)
    raw = value if isinstance(value, str) else value.get("raw_text") or f"{value.get('number')} {value.get('unit')}"
    o["value"] = f(raw, value, ev)
    if basis is not None:
        o["measurement_basis"] = f(basis, basis, ev)
    if scope is not None:
        o["scope"] = f(scope, scope, ev)
    if yield_type is not None:
        o["yield_type"] = f(yield_type, yield_type, ev)
    r["outcomes"].append(o)


def set_solvent(r, composition, volume, ev):
    r["solvent_conditions"]["composition"] = f(composition, composition, ev)
    if volume is not None:
        r["solvent_conditions"]["total_volume"] = f(f"{volume['number']} {volume['unit']}", volume, ev)


def petasis_class(r, successful=True):
    if not successful:
        classification(r, [], [], "not_applicable", "No observable desired transformation; intended bonds are not classified.")
        return
    changes = [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "N"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
    ]
    labels = [
        {"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1"]},
        {"category_number": 5, "code": "C-N", "supporting_change_ids": ["b2"]},
    ]
    classification(r, labels, changes)


reactions = []

# Prior work and current general transformation.
r = reaction("literature_cross_metathesis", (I[2],), "literature_example", (),
             "page 2 Scheme 1a", "ethylene-alkyne cross-metathesis", "cited cross-metathesis")
r["participants"] = [participant("lit1_s", lit1_s, "reactant", (I[2],)), participant("lit1_p", lit1_p, "product", (I[2],))]
classification(r, [], [], "unresolved", "Cited generic family; exact atom mapping and experiment conditions are outside the supplied source detail.")
r["review"]["unresolved"] = ["Cited literature family only; detailed conditions and individual experiments were not imported."]
reactions.append(r)

r = reaction("literature_aminol_dienylation", (I[2],), "literature_example", (),
             "page 2 Scheme 1b", "1,3-butadienylation and enzymatic kinetic resolution", "cited multistep dienylation/resolution")
r["participants"] = [participant("lit2_s", lit2_s, "reactant", (I[2],)), participant("lit2_p", lit2_p, "product", (I[2],))]
classification(r, [], [], "unresolved", "Cited generic multistep sequence; exact conditions are not represented on the supplied page.")
r["review"]["unresolved"] = ["Cited literature family only; individual operations and yields are not stated in Scheme 1b."]
reactions.append(r)

GENERAL = RID("generic_petasis")
r = reaction("generic_petasis", (I[2], I[3]), "generic_scheme", (),
             "pages 2-3 Schemes 1c-3", "asymmetric Petasis reaction", "three-component Petasis dienylation")
r["participants"] = [participant("gen_aldehyde", gen_aldehyde, "reactant", (I[2], I[3])),
                     participant("gen_2a", bor2a, "reactant", (I[2], I[3])),
                     participant("gen_2b", bor2b, "reactant", (I[2], I[3])),
                     participant("gen_3e", aminols["3e"], "reactant", (I[2], I[3])),
                     participant("gen_4", gen4, "product", (I[2], I[3])),
                     participant("gen_5", gen5, "product", (I[2], I[3]))]
petasis_class(r)
r["review"]["unresolved"] = ["Generic R-group scheme; it is not one executed experiment and carries no yield."]
reactions.append(r)


# Table 1: every printed row is an independent optimization record.
table = [
    (1, "3a", "MeOH", True, 48, 68, "2:1"),
    (2, "3b", "MeOH", True, 48, 65, "1.3:1"),
    (3, "3c", "MeOH", True, 72, 86, "2:1"),
    (4, "3d", "MeOH", True, 72, None, None),
    (5, "3e", "MeOH", True, 12, 94, "3:1"),
    (6, "3e", "MeOH", False, 4, 95, "3.4:1"),
    (7, "3e", "THF", False, 48, 52, "3:1"),
    (8, "3e", "toluene", False, 48, 72, "3:1"),
    (9, "3e", "CH2Cl2", False, 48, 75, "3.5:1"),
    (10, "3e", "EtOAc", False, 72, 87, "3.6:1"),
    (11, "3e", "DMSO", False, 20, 95, "11.3:1"),
    (12, "3f", "DMSO", False, 20, None, None),
]
for entry, am, solv, sieve, hours, yld, dr in table:
    ev = (I[2], I[3]) if entry == 11 else (I[2],)
    label = f"table1_entry{entry}"
    r = reaction(label, ev, "specific_experiment", ("optimization",) + (("substrate_scope",) if entry == 11 else ()),
                 f"page 2 Table 1 entry {entry}" + ("; page 3 Scheme 2 product 4a" if entry == 11 else ""),
                 "asymmetric Petasis reaction", "Petasis 1,3-butadienylation", (GENERAL,))
    r["participants"] = [
        participant("aldehyde", AM["1a"], "reactant", ev, amount=q(0.2, "mmol")),
        participant("aminol", aminols[am], "reactant", ev, amount=q(0.24, "mmol")),
        participant("boronate", bor2a, "reactant", ev, amount=q(0.30, "mmol")),
        participant("solvent", SM[solv], "solvent", ev, volume=q(0.4, "mL")),
    ]
    if sieve:
        r["participants"].append(participant("molecular_sieves", ms4a, "reagent", ev, mass=q(50.0, "mg")))
    if entry >= 7:
        r["participants"].append(participant("methanol", SM["MeOH"], "additive", ev,
                                             equivalents=q(5.0, "equiv"), ref=AM["1a"]))
    product = None
    if yld is not None:
        product = table_products.get(am, prod4["4a"])
        r["participants"].append(participant("product", product, "product", ev))
    set_solvent(r, "MeOH" if entry <= 6 else f"{solv} + 5.0 equiv MeOH", q(0.4, "mL"), ev)
    pids = [p["participant_id"] for p in r["participants"]]
    add_step(r, "condensation", ev, pids, "aldehyde, aminol, solvent, and molecular sieves when specified were stirred", "r.t.", q(2, "h"))
    add_step(r, "boronate_addition_and_hold", ev, pids, "2a was added and the reaction mixture was held", "r.t.", q(hours, "h"))
    if yld is None:
        add_outcome(r, "no_reaction", "no_reaction", None, "No reaction", ev,
                    basis=None, scope={"desired_product": "4"})
        petasis_class(r, False)
        r["review"]["unresolved"] = ["No product, conversion, recovery, byproduct, or detection limit is reported."]
    else:
        add_outcome(r, "yield", "yield", product, q(yld, "%"), ev, basis="isolated yield", yield_type="isolated")
        add_outcome(r, "dr", "selectivity", product, dr, ev, basis="1H NMR",
                    scope={"metric": "dr", "order": "major:minor not structurally assigned"})
        petasis_class(r)
        r["review"]["unresolved"] = ["Complete product atom mapping and full workup/purification are unavailable because SI was not supplied."]
    if entry == 11:
        r["notes"].append(f("Table 1 entry 11 and Scheme 2 product 4a describe the same conditions and 95%/11.3:1 result; evidence was merged.",
                            "cross-page duplicate merged", ev))
    reactions.append(r)

# Prose-only elevated-temperature attempt with 3f.
r = reaction("aminol3f_elevated_temperature_control", (I[2],), "specific_experiment", ("optimization", "control"),
             "page 2 prose below Table 1", "asymmetric Petasis reaction", "Petasis control", (GENERAL,))
r["participants"] = [participant("aldehyde", AM["1a"], "reactant", (I[2],)),
                     participant("aminol", aminols["3f"], "reactant", (I[2],)),
                     participant("boronate", bor2a, "reactant", (I[2],)),
                     participant("DMSO", SM["DMSO"], "solvent", (I[2],)),
                     participant("MeOH", SM["MeOH"], "additive", (I[2],))]
set_solvent(r, "DMSO + MeOH", None, (I[2],))
add_step(r, "elevated_temperature_attempt", (I[2],), [p["participant_id"] for p in r["participants"]],
         "3f was tested under otherwise identical conditions at an elevated temperature", "elevated temperature", None)
add_outcome(r, "no_reaction", "no_reaction", None, "No reaction", (I[2],), scope={"desired_product": "4"})
petasis_class(r, False)
r["review"]["unresolved"] = ["The elevated temperature, duration, exact charges, and detection method are not reported separately from Table 1 entry 12."]
reactions.append(r)


# Scheme 2 scope; 4a is already Table 1 entry 11.
scheme2 = {
    "4b": ("pMe", 14, 98, "5.3:1"), "4c": ("pOMe", 14, 93, "6.5:1"),
    "4d": ("pF", 14, 94, ">20:1"), "4e": ("pCl", 14, 90, ">20:1"),
    "4f": ("pBr", 14, 92, ">20:1"), "4g": ("pCO2Me", 18, 94, ">20:1"),
    "4h": ("pNO2", 14, 93, ">20:1"), "4i": ("mOMe", 12, 95, ">20:1"),
    "4j": ("mBr", 12, 98, ">20:1"), "4k": ("mCF3", 16, 90, ">20:1"),
    "4l": ("mCN", 16, 88, ">20:1"), "4m": ("mCO2Me", 24, 92, ">20:1"),
    "4n": ("oF", 24, 86, ">20:1"), "4o": ("dibromo35", 48, 97, ">20:1"),
    "4p": ("bromo_nitro35", 48, 96, ">20:1"), "4q": ("naph2", 48, 92, ">20:1"),
    "4r": ("butanal", 12, 92, ">20:1"), "4s": ("isobutyr", 12, 90, "12:1"),
    "4t": ("cycho", 18, 92, ">20:1"),
}
for lab, (ald, hours, yld, dr) in scheme2.items():
    ev = (I[3],)
    r = reaction(f"scheme2_{lab}", ev, "specific_experiment", ("substrate_scope",),
                 f"page 3 Scheme 2 product {lab}", "asymmetric Petasis reaction", "Petasis 1,3-butadienylation", (GENERAL,))
    r["participants"] = [
        participant("aldehyde", AM[ald], "reactant", ev, amount=q(0.2, "mmol")),
        participant("aminol", aminols["3e"], "reactant", ev, amount=q(0.24, "mmol")),
        participant("MeOH", SM["MeOH"], "additive", ev, equivalents=q(5.0, "equiv"), ref=AM[ald]),
        participant("boronate", bor2a, "reactant", ev, amount=q(0.30, "mmol")),
        participant("DMSO", SM["DMSO"], "solvent", ev, volume=q(0.4, "mL")),
        participant("product", prod4[lab], "product", ev),
    ]
    set_solvent(r, "DMSO + 5.0 equiv MeOH", q(0.4, "mL"), ev)
    pids = [p["participant_id"] for p in r["participants"]]
    add_step(r, "condensation", ev, pids, "aldehyde, 3e, MeOH, and DMSO were stirred", "r.t.", q(2, "h"))
    add_step(r, "boronate_addition_and_hold", ev, pids, "2a was added and the mixture was held", "r.t.", q(hours, "h"))
    add_outcome(r, "yield", "yield", prod4[lab], q(yld, "%"), ev, basis="isolated yield", yield_type="isolated")
    add_outcome(r, "dr", "selectivity", prod4[lab], dr, ev, basis="1H NMR",
                scope={"metric": "dr", "order": "major:minor not structurally assigned"})
    petasis_class(r)
    r["review"]["unresolved"] = ["Full product structure encoding, workup, and purification require unavailable SI."]
    if lab == "4p":
        r["notes"].append(f("The 3-bromo-5-nitro substitution pattern was transcribed from the compact Scheme 2 drawing and remains a machine candidate for expert confirmation.",
                            "candidate 3-bromo-5-nitro assignment", ev))
    reactions.append(r)


# Scheme 3 scope.
scheme3 = {
    "5a": ("1a", 97, "19:1"), "5b": ("pMe", 96, "8.5:1"), "5c": ("pIPr", 90, "10:1"),
    "5d": ("pTBu", 93, ">20:1"), "5e": ("pCl", 95, ">20:1"), "5f": ("pCN", 99, ">20:1"),
    "5g": ("pCO2Me", 95, ">20:1"), "5h": ("mMe", 85, ">20:1"), "5i": ("mF", 87, ">20:1"),
    "5j": ("mCl", 84, ">20:1"), "5k": ("mBr", 90, ">20:1"), "5l": ("mCN", 95, ">20:1"),
    "5m": ("mCO2Me", 98, ">20:1"), "5n": ("oOMe", 85, "6:1"), "5o": ("oF", 90, ">20:1"),
    "5p": ("oCF3", 86, "8:1"), "5q": ("oNO2", 92, "10:1"), "5r": ("dimethyl34", 96, ">20:1"),
    "5s": ("dibromo35", 91, ">20:1"), "5t": ("naph2", 82, ">20:1"), "5u": ("fur2", 99, "8.5:1"),
    "5v": ("thio2", 86, "10:1"), "5w": ("pyr2", 80, ">20:1"), "5x": ("isoval", 82, ">20:1"),
    "5y": ("cycho", 85, ">20:1"),
}
for lab, (ald, yld, dr) in scheme3.items():
    ev = (I[3],)
    r = reaction(f"scheme3_{lab}", ev, "specific_experiment", ("substrate_scope",),
                 f"page 3 Scheme 3 product {lab}", "asymmetric Petasis reaction", "Petasis isoprenylation", (GENERAL,))
    r["participants"] = [
        participant("aldehyde", AM[ald], "reactant", ev, amount=q(0.2, "mmol")),
        participant("aminol", aminols["3e"], "reactant", ev, amount=q(0.24, "mmol")),
        participant("MeOH", SM["MeOH"], "additive", ev, equivalents=q(5.0, "equiv"), ref=AM[ald]),
        participant("boronate", bor2b, "reactant", ev, amount=q(0.30, "mmol")),
        participant("DMSO", SM["DMSO"], "solvent", ev, volume=q(0.2, "mL")),
        participant("product", prod5[lab], "product", ev),
    ]
    set_solvent(r, "DMSO + 5.0 equiv MeOH", q(0.2, "mL"), ev)
    pids = [p["participant_id"] for p in r["participants"]]
    add_step(r, "condensation", ev, pids, "aldehyde, 3e, MeOH, and DMSO were stirred", "r.t.", q(2, "h"))
    add_step(r, "boronate_addition_and_hold", ev, pids, "2b was added and the mixture was held", "r.t.", q(12, "h"))
    add_outcome(r, "yield", "yield", prod5[lab], q(yld, "%"), ev, basis=None)
    add_outcome(r, "dr", "selectivity", prod5[lab], dr, ev, basis=None,
                scope={"metric": "dr", "order": "major:minor; identities not fully assigned"})
    petasis_class(r)
    r["review"]["unresolved"] = ["Scheme 3 does not state yield or dr measurement methods on the supplied page; workup and purification require unavailable SI."]
    if lab == "5u":
        r["notes"].append(f("The newly formed stereogenic carbon of 5u was assigned S by X-ray crystallography (CCDC 2023204).", "S assignment by X-ray", ev))
    elif lab <= "5w":
        r["notes"].append(f("The authors assign other aryl isoprenylated products by analogy to 5u.", "assignment by analogy", ev))
    reactions.append(r)


# Scheme 4 downstream synthetic utility.
def downstream(label, reactants, product, solvent_key, temp, hours, yld, dr, change_kind):
    ev = (I[4],)
    r = reaction(label, ev, "specific_experiment", ("derivatization",),
                 f"page 4 Scheme 4 {label}", None, change_kind)
    r["participants"] = [participant(pid, mid, "reactant", ev) for pid, mid in reactants]
    r["participants"].append(participant("solvent", SM[solvent_key], "solvent", ev))
    r["participants"].append(participant("product", product, "product", ev))
    set_solvent(r, molecules[SM[solvent_key]]["names"]["original"]["normalized_value"], None, ev)
    add_step(r, "reaction", ev, [p["participant_id"] for p in r["participants"]], "reaction under the displayed Scheme 4 conditions", q(temp, "°C") if isinstance(temp, (int, float)) else temp, q(hours, "h"))
    add_outcome(r, "yield", "yield", product, q(yld, "%"), ev, basis=None)
    if dr:
        add_outcome(r, "dr", "selectivity", product, dr, ev, basis=None,
                    scope={"metric": "dr", "order": "major:minor not structurally assigned"})
    r["review"]["unresolved"] = ["Reactant charges, workup, purification, and yield/selectivity measurement methods are not reported on page 4; SI was not supplied."]
    return r

r = downstream("scheme4_diels_alder_6", (("5a", prod5["5a"]), ("dimethyl_fumarate", SM["dimethyl_fumarate"])), prod6, "toluene", 100, 6, 88, "1:1", "Diels-Alder cycloaddition")
changes = [
    {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
    {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
]
classification(r, [{"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1", "b2"]}], changes)
reactions.append(r)

r = downstream("scheme4_diels_alder_7", (("5a", prod5["5a"]), ("diethyl_acetylenedicarboxylate", SM["diethyl_acetylenedicarboxylate"])), prod7, "anisole", 180, 1, 78, ">20:1", "Diels-Alder cycloaddition")
classification(r, [{"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1", "b2"]}], changes)
reactions.append(r)

r = downstream("scheme4_methylation_8", (("7", prod7), ("MeI", SM["MeI"]), ("Ag2O", SM["Ag2O"])), prod8, "CH2Cl2", "r.t.", 48, 85, ">20:1", "methylation")
ch = [{"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}], "before_order": 0, "after_order": 1, "change_type": "formed"}]
classification(r, [{"category_number": 3, "code": "C-O", "supporting_change_ids": ["b1"]}], ch)
reactions.append(r)

r = downstream("scheme4_oxidation_9", (("8", prod8), ("MnO2", SM["MnO2"])), prod9, "dioxane", "r.t.", 48, 82, ">20:1", "oxidative aromatization")
ch = [{"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 1, "after_order": 2, "change_type": "order_increased"}]
classification(r, [{"category_number": 13, "code": "DB-C-C", "supporting_change_ids": ["b1"]}], ch,
               limitation="Class-level aromatization records a C-C bond-order promotion; exact bond identity, multiplicity, and C-H losses are unresolved.")
reactions.append(r)


bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": RID("context_page1"), "document_id": DOC, "kind": "metadata_and_natural_product_context",
     "evidence_ids": [I[1], P[1]], "summary": "Title/metadata/introduction and Figure 1 natural-product structures; no current concrete experiment."},
    {"context_id": RID("context_page3_structure"), "document_id": DOC, "kind": "structure_assignment",
     "evidence_ids": [I[3], P[3]], "summary": "ORTEP/CCDC 2023204 assigns the newly formed center of 5u as S; other aryl products are assigned by analogy."},
    {"context_id": RID("context_page3_dft"), "document_id": DOC, "kind": "computational_mechanistic_context",
     "evidence_ids": [I[3], P[3]], "summary": "Figure 2 proposed TS-Si/TS-Re models and the reported -4.44 kcal/mol barrier difference; not an executed synthesis experiment."},
    {"context_id": RID("context_page4_references"), "document_id": DOC, "kind": "bibliography_and_negative_spectral_inventory",
     "evidence_ids": [I[4], P[4]], "summary": "Summary, acknowledgements, references, and external CCDC notice. No standalone mass-spectrum peak table, molecular-weight table, or fragment-relative-abundance block occurs on pages 1-4."},
]

# Complete document title from the source image.
bundle["documents"][0]["title"] = "Asymmetric Petasis reaction for the synthesis of chiral alpha- and beta-butadienyl amines"

stereo_checks = {}
for mid, m in molecules.items():
    smiles = m["smiles"]["normalized_value"]
    if smiles:
        stereo_checks[mid] = st.stereo_inspect(smiles)
(RUN / "raw_model/stereo_checks.json").write_text(json.dumps(stereo_checks, ensure_ascii=False, indent=2) + "\n")

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 1,
            "state": "attempted",
            "regions": [
                {"region": "title, citation, abstract, and ESI/CCDC notice", "finding": "Article metadata and summary; ESI is advertised but no SI file is supplied or associated.", "evidence_ids": [I[1], P[1]]},
                {"region": "two-column introduction", "finding": "Background and cited prior methods; no current concrete experiment arrow or optimization row.", "evidence_ids": [I[1], P[1]]},
                {"region": "Figure 1", "finding": "Natural-product alpha-butadienyl amine context; non-synthetic source context for this paper.", "evidence_ids": [I[1]]},
            ],
            "record_ids": [],
            "unresolved": ["No SI was provided despite the availability notice; procedures and characterization outside pages 1-4 remain unavailable."],
        },
        {
            "page": 2,
            "state": "attempted",
            "regions": [
                {"region": "Scheme 1a-b", "finding": "Two cited prior-work reaction families, retained as literature examples rather than current experiments.", "evidence_ids": [I[2], P[2]]},
                {"region": "Scheme 1c", "finding": "Current generic Petasis transformation with 2a/2b; not counted as an executed experiment.", "evidence_ids": [I[2], P[2]]},
                {"region": "Table 1 rows 1-12 and footnotes", "finding": "Twelve independent optimization entries; entries 4 and 12 are no-reaction cases. Entry 11 duplicates Scheme 2 product 4a and is merged across pages.", "evidence_ids": [I[2], P[2]]},
                {"region": "optimization discussion below Table 1", "finding": "A distinct prose-only elevated-temperature 3f attempt also gave no reaction; temperature and duration are not reported.", "evidence_ids": [I[2], P[2]]},
                {"region": "right-column scope discussion", "finding": "Narrative cross-check of Scheme 2 yields/dr and Scheme 3 setup; no additional duplicate experiment created.", "evidence_ids": [I[2], P[2]]},
            ],
            "record_ids": [RID("literature_cross_metathesis"), RID("literature_aminol_dienylation"), RID("generic_petasis")]
                          + [RID(f"table1_entry{i}") for i in range(1, 13)]
                          + [RID("aminol3f_elevated_temperature_control")],
            "unresolved": ["Aminol structures 3a-3f are compact; complete source-verified SMILES were not forced.", "The separate elevated-temperature value for 3f is not stated."],
        },
        {
            "page": 3,
            "state": "attempted",
            "regions": [
                {"region": "Scheme 2 and footnote", "finding": "All 20 alpha-butadienyl products 4a-4t inventoried; 4a is the Table 1 entry-11 duplicate and 4b-4t are 19 additional scope experiments.", "evidence_ids": [I[3], P[3]]},
                {"region": "Scheme 3 and footnote", "finding": "All 25 beta-butadienyl products 5a-5y inventoried as independent scope experiments.", "evidence_ids": [I[3], P[3]]},
                {"region": "ORTEP of 5u", "finding": "X-ray/CCDC 2023204 assigns the newly formed center of major 5u as S; other aryl products are assigned by analogy.", "evidence_ids": [I[3], P[3]]},
                {"region": "Figure 2 and DFT discussion", "finding": "Proposed TS-Si/TS-Re computational models, retained as mechanistic source context rather than an executed reaction.", "evidence_ids": [I[3], P[3]]},
            ],
            "record_ids": [RID("generic_petasis"), RID("table1_entry11")]
                          + [RID(f"scheme2_4{c}") for c in "bcdefghijklmnopqrst"]
                          + [RID(f"scheme3_5{c}") for c in "abcdefghijklmnopqrstuvwxy"],
            "unresolved": ["Scheme 3 does not state yield or dr measurement methods on page 3.", "Complete product structures and full stereochemical encoding are unavailable from the compact drawings without SI."],
        },
        {
            "page": 4,
            "state": "attempted",
            "regions": [
                {"region": "Scheme 4 upper arrows", "finding": "Two independent Diels-Alder derivatizations of 5a give 6 and 7.", "evidence_ids": [I[4], P[4]]},
                {"region": "Scheme 4 lower arrows", "finding": "Methylation 7 to 8 and oxidative aromatization 8 to 9 are two additional experiments.", "evidence_ids": [I[4], P[4]]},
                {"region": "discussion and conclusion", "finding": "Narrative confirms Scheme 4 outcomes and summarizes the Petasis study; no extra experiment created.", "evidence_ids": [I[4], P[4]]},
                {"region": "references and CCDC notice", "finding": "Bibliography/external crystallographic availability only. No standalone mass-spectrum, molecular-weight, or fragment-abundance experiment is present.", "evidence_ids": [I[4], P[4]]},
            ],
            "record_ids": [RID("scheme4_diels_alder_6"), RID("scheme4_diels_alder_7"), RID("scheme4_methylation_8"), RID("scheme4_oxidation_9")],
            "unresolved": ["Scheme 4 omits charges, yield/dr measurement methods, workup, and purification; these were not inferred from unavailable SI."],
        },
    ],
}
(RUN / "coverage").mkdir(parents=True, exist_ok=True)
(RUN / "coverage/pages_0001_0004.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")

page_extraction = {
    "document_id": DOC,
    "page_images_visually_inspected": [1, 2, 3, 4],
    "image_evidence_ids": [I[1], I[2], I[3], I[4]],
    "inventory": {
        "page1": {"reaction_arrows": 0, "table_rows": 0, "notes": ["Figure 1 is biological/natural-product context."]},
        "page2": {"cited_scheme_arrows": 2, "generic_current_schemes": 1, "table_rows": 12, "prose_only_controls": 1,
                  "deduplication": ["Table 1 entry 11 equals Scheme 2 product 4a."]},
        "page3": {"scheme2_outcomes": 20, "scheme3_outcomes": 25, "executed_mechanistic_experiments": 0,
                  "context": ["5u ORTEP/X-ray assignment", "DFT proposed transition states"]},
        "page4": {"scheme4_arrows": 4, "table_rows": 0, "notes": ["Conclusion and references are non-experimental."]},
    },
    "negative_inventory": "No standalone mass-spectrum peak table, molecular-weight table, or fragment-relative-abundance block was found on pages 1-4.",
    "limitations": ["Visual inspection and schema validation do not establish expert chemical accuracy or coverage outside pages 1-4.", "No supporting-information file was supplied."],
}
(RUN / "raw_model/page_extraction.json").write_text(json.dumps(page_extraction, ensure_ascii=False, indent=2) + "\n")

mentions = st.source_mentions(json.loads((RUN / "dataset.json").read_text()), ["4a", "5a", "5u", "7", "3f"])
(RUN / "raw_model/source_mentions.json").write_text(json.dumps(mentions, ensure_ascii=False, indent=2) + "\n")

(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
                  "source_context": len(bundle["source_context"])}, ensure_ascii=False))
