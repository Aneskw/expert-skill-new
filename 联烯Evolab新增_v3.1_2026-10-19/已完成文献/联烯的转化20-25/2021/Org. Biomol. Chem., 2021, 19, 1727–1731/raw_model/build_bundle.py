#!/usr/bin/env python3
"""Build the page-1--5 machine extraction bundle for doc_333593389c9636f2.

This is a source-derived candidate generator.  It deliberately leaves uncertain
structures and unavailable SI details unresolved; merge will quarantine all
records as needs_review.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_333593389c9636f2"
P = {i: f"{DOC}:page{i}:text" for i in range(1, 6)}
I = {i: f"{DOC}:page{i}:image" for i in range(1, 6)}
DATASET = ROOT / "runs/current/dataset.json"
OUT = ROOT / "runs/current/submitted.json"
RAW = ROOT / "runs/current/raw_model/pages_0001_0005_inventory.json"
DIAG = ROOT / "runs/current/raw_model/stereo_diagnostics.json"
COVERAGE = ROOT / "runs/current/coverage/pages_0001_0005.json"


def F(raw=None, normalized=..., ev=(), method="original_extraction", derivation=None,
      status=None):
    if normalized is ...:
        normalized = raw
    return st.fact(raw, normalized, ev, method=method, derivation=derivation,
                   status=status)


def Q(number, unit, raw=None):
    return raw or f"{number} {unit}", {"number": number, "unit": unit}


def set_status(m, key, normalized, ev, raw=None, derivation=None):
    method = "structure_derivation" if derivation else "original_extraction"
    m["stereochemistry"][key] = F(raw if raw is not None else normalized, normalized,
                                     ev, method, derivation)


def with_page_images(ev):
    """Add the rendered page evidence for facts transcribed from drawings."""
    expanded = list(ev)
    for page, text_id in P.items():
        if text_id in ev and I[page] not in expanded:
            expanded.append(I[page])
    return expanded


def molecule(mid, label, ev, name=None, smiles=None, stereo="unresolved",
             geometry="unresolved", representation="unresolved", identity_note=None,
             mixture=None, axial=None, assignment=None):
    m = st.new_molecule(f"{DOC}:mol:{mid}", DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, ev)
    m["names"]["original"] = F(label, label, ev)
    if name and name != label:
        m["names"]["english"] = F(None, name, ev, "inference",
                                            "English expansion/structural description normalized from the source label and drawing; not quoted as a source name.")
    else:
        m["names"]["english"] = F(label, label, ev)
    m["names"]["abbreviations"] = F()
    m["names"]["aliases"] = F()
    m["cas_numbers"] = F()
    if smiles:
        structure_ev = with_page_images(ev)
        m["representation_type"] = "smiles"
        m["smiles"] = F(f"structure drawing for {label}", smiles, structure_ev,
                         "structure_derivation",
                         "Candidate connectivity transcribed from the source drawing; RDKit parsing is not source verification.")
        diag = st.stereo_inspect(smiles)
        m["stereochemistry"]["identity_key"] = F(None, diag["stereo_graph_key"], structure_ev,
                                                    "structure_derivation",
                                                    "RDKit stereo graph key for the candidate SMILES.")
        if representation in ("encoded", "not_applicable"):
            m["stereochemistry"]["isomeric_smiles"] = F(None, diag["canonical_isomeric_smiles"], structure_ev,
                                                           "structure_derivation",
                                                           "Canonical isomeric SMILES calculated from the candidate graph.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = F()
    m["generic_template"] = None
    set_status(m, "chirality_status", stereo, ev,
               derivation="Status assigned from the depicted/source-described stereochemical state.")
    set_status(m, "geometry_status", geometry, ev,
               derivation="Status assigned from the depicted/source-described geometric state.")
    set_status(m, "representation_status", representation, ev,
               derivation="Compared source stereochemical information with the candidate representation.")
    m["stereochemistry"]["mixture"] = F(mixture, mixture, ev) if mixture else F()
    m["stereochemistry"]["axial"] = F(axial, axial, ev) if axial else F()
    m["stereochemistry"]["assignment"] = F(assignment, assignment, ev) if assignment else F()
    m["stereochemistry"]["atom_centers"] = F()
    m["stereochemistry"]["double_bonds"] = F()
    m["identity_note"] = F(identity_note, identity_note, ev) if identity_note else F()
    return m


def generic_molecule(mid, label, ev, fixed, variables):
    m = molecule(mid, label, ev, name=label, stereo="unresolved",
                 geometry="unresolved", representation="not_encoded")
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": F(fixed, {"description": fixed}, ev),
        "variables": [
            {"label": v, "attachment_points": F(f"source attachment for {v}", [v], ev),
             "definition": F(f"source variable {v}", None, ev, status="unreadable"),
             "allowed_values": F(), "constraints": F()}
            for v in variables
        ],
        "instance_mapping": [],
        "reason": F("generic R-group drawing", "The source drawing contains variable substituents and is not a single concrete molecule.", ev),
    }
    return m


def participant(pid, mid, role, ev, stage="reaction", amount=None, volume=None,
                equivalents=None, mol_percent=None, ref=None, product_role=None):
    p = st.new_participant(pid, f"{DOC}:mol:{mid}", role, ev, stage)
    p["product_role"] = F(product_role, product_role, ev) if product_role else F()
    for key, val in (("amount", amount), ("volume", volume),
                     ("equivalents", equivalents), ("mol_percent", mol_percent)):
        if val is not None:
            raw, norm = val
            p["charges"][key] = F(raw, norm, ev)
    if ref:
        p["equivalence_reference"] = F(ref, f"{DOC}:mol:{ref}", ev)
    return p


def classification(ev, pairs, note):
    ev = with_page_images(ev)
    changes = []
    labels = {}
    for i, item in enumerate(pairs, 1):
        a, b, kind, *details = item
        if len(details) >= 2:
            before, after = details[:2]
        else:
            before, after = (0, 1) if kind in ("formed", "hydrogen_gain") else (1, 0)
        ch = {"change_id": f"b{i}", "atoms": [{"atom_map": None, "element": a},
              {"atom_map": None, "element": b}], "before_order": before,
              "after_order": after, "change_type": kind}
        if kind.startswith("hydrogen_"):
            ch["hydrogen_count"] = details[2] if len(details) >= 3 else 1
        changes.append(F(None, ch, ev, "structure_derivation", note))
        if kind in ("formed", "hydrogen_gain"):
            lab = st.bond_label(a, b, after)
            labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(f"b{i}")
    return {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": F(),
        "bond_changes": changes,
        "labels": F(None, list(labels.values()), ev, "structure_derivation", note),
        "coverage": F(None, "partial" if pairs else "unresolved", ev, "inference",
                      "Only source-visible net changes were classified; no complete atom mapping or coproduct balance is available."),
        "limitations": ["No complete atom mapping was asserted.", note],
        "model_votes": [],
    }


def base_reaction(suffix, ev, record_type="specific_experiment", categories=None,
                  location=None, name="[2 + 4]-annulation", transformation="annulation",
                  pairs=(("C", "N", "formed"), ("C", "O", "formed"))):
    rid = f"{DOC}:{suffix}"
    r = st.new_reaction(rid, DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = F(location or suffix, categories or [], ev)
    r["location"] = F(location or suffix, {"source": location or suffix}, ev)
    r["reaction_name"] = F(name, name, ev)
    r["transformation_type"] = F(None, transformation, ev, "inference",
                                  "Normalized from the source reaction name and depicted transformation.")
    r["classification"] = classification(ev, pairs, "Bond classes inferred from the source reaction drawings; product atom identities are not fully mapped.")
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [],
                   "unresolved": ["No SI was supplied; experimental workup, purification and characterization details remain unavailable.",
                                  "Candidate structures and bond classes require expert review."]}
    r["completeness"] = "partial"
    return r


def step(r, sid, ev, pids, operation=None, temperature=None, time=None):
    s = st.new_step(sid, ev)
    s["participant_ids"] = pids
    if operation is not None:
        s["operation"] = F(operation, operation, ev)
    if temperature is not None:
        s["temperature"] = F(*temperature, ev)
    if time is not None:
        s["time"] = F(*time, ev)
    r["steps"].append(s)


def outcome(r, oid, kind, target, ev, raw, norm, yield_type=None, basis=None, scope=None):
    o = st.new_outcome(oid, kind, [f"{DOC}:mol:{target}"], ev)
    o["product_role"] = F("product", "product", ev)
    o["value"] = F(raw, norm, ev)
    if yield_type:
        o["yield_type"] = F(yield_type, yield_type, ev)
    if basis:
        o["measurement_basis"] = F(basis, basis, ev)
    if scope:
        o["scope"] = F(scope, scope, ev)
    r["outcomes"].append(o)


data = json.loads(DATASET.read_text())
bundle = st.empty(data["documents"])
bundle["documents"][0]["title"] = (
    "Organocatalytic enantioselective [2 + 4]-annulation of γ-substituted "
    "allenoates with N-acyldiazenes for the synthesis of optically active "
    "1,3,4-oxadiazines"
)
molecules = {}


def add(m):
    molecules[m["molecule_id"].split(":")[-1]] = m


# Core and directly reusable concrete candidates.
add(molecule("1a", "(±)-1a", [P[1], P[2], P[3], P[4]], "racemic γ-benzyl allenoate 1a",
             "CCOC(=O)C=C=CCc1ccccc1", "racemic", "unspecified", "not_encoded",
             "The source explicitly calls 1a racemic; the candidate SMILES does not encode allene axial configuration.",
             "racemic pair; enantiomer ratio not separately reported", "allene axis; individual R_a/S_a assignment not reported"))
add(molecule("2aa", "2aa", [P[2], P[4]], "N-phenyl-N′-benzoyldiazene 2aa",
             "O=C(c1ccccc1)N=Nc1ccccc1", "achiral", "unspecified", "partial",
             "Candidate connectivity from the source drawing; N=N geometry is not specified."))
add(molecule("3aaa", "3aaa", [P[1], P[2], P[3], P[4]], "1,3,4-oxadiazine 3aaa",
             "CCOC(=O)C=C1OC(=NN(c2ccccc2)C1Cc2ccccc2)c2ccccc2",
             stereo="unspecified", geometry="unspecified", representation="not_encoded",
             identity_note="Connectivity candidate follows the legible product drawings on pages 1-4. The ring stereogenic carbon is marked with an asterisk, but no absolute configuration is assigned; the exocyclic alkene has no written E/Z assignment, and ee must not be converted to R/S."))

# Catalysts and common solvents/reagents.
for c in range(1, 7):
    assignment = "Source drawing contains stereobonds; no absolute R/S descriptor was transcribed."
    add(molecule(f"C{c}", f"C{c}", [P[2]], f"chiral DMAP analogue C{c}",
                 stereo="specified", geometry="unresolved",
                 representation="not_encoded", assignment=assignment))

simple = {
    "toluene": ("toluene", "Cc1ccccc1"), "EtOAc": ("ethyl acetate", "CCOC(C)=O"),
    "MeCN": ("acetonitrile", "CC#N"), "CH2Cl2": ("dichloromethane", "ClCCl"),
    "CHCl3": ("chloroform", "ClC(Cl)Cl"), "Et2O": ("diethyl ether", "CCOCC"),
    "THF": ("tetrahydrofuran", "C1CCOC1"), "dioxane": ("1,4-dioxane", "O1CCOCC1"),
    "MTBE": ("methyl tert-butyl ether", "COC(C)(C)C"), "Br2": ("bromine", "BrBr"),
    "Et3N": ("triethylamine", "CCN(CC)CC"), "DCC": ("dicyclohexylcarbodiimide", "C(=NC1CCCCC1)=NC1CCCCC1"),
    "DMAP": ("4-dimethylaminopyridine", "CN(C)c1ccncc1"),
}
for mid, (nm, smi) in simple.items():
    ev = [P[2]] if mid in {"toluene", "EtOAc", "MeCN", "CHCl3", "Et2O", "THF", "dioxane", "MTBE"} else [P[4]]
    if mid == "CH2Cl2": ev = [P[2], P[4]]
    add(molecule(mid, mid, ev, nm, smi, "achiral", "not_applicable", "not_applicable"))
add(molecule("DIBAL-H", "DIBAL-H", [P[4]], "diisobutylaluminium hydride", stereo="unresolved", geometry="unresolved", representation="unresolved"))
add(molecule("acid1", "acid 1", [P[4]], "4-bromo-2-methylbenzoic acid (acid 1, candidate from drawing)",
             "Cc1cc(Br)ccc1C(=O)O", "achiral", "not_applicable", "not_applicable",
             "Candidate regiochemistry transcribed from the Scheme 2 drawing; it remains a machine candidate pending expert review."))
add(molecule("acid2", "acid 2", [P[4]], "acid 2 (polycyclic carboxylic acid as drawn)", stereo="unresolved", geometry="unresolved", representation="unresolved"))

# Scope diazenes and products.  Only the two repeatedly named diazenes receive
# concrete candidate SMILES; the remaining identities retain source labels and
# Ar1/Ar2 definitions without invented full structures.
t2 = [
    ("aa", "Ph/Ph", 48, 82, 71, "3aaa"), ("ab", "Ph/4-FC6H4", 48, 72, 72, "3aab"),
    ("ac", "Ph/4-ClC6H4", 48, 73, 71, "3aac"), ("ad", "Ph/4-BrC6H4", 60, 66, 70, "3aad"),
    ("ae", "Ph/4-MeC6H4", 48, 76, 70, "3aae"), ("af", "Ph/2-MeC6H4", 48, 88, 71, "3aaf"),
    ("ag", "Ph/1-naphthyl", 60, 77, 70, "3aag"), ("ah", "Ph/3,4-(MeO)2C6H3", 60, 71, 74, "3aah"),
    ("ai", "Ph/3,5-(MeO)2C6H3", 72, 72, 72, "3aai"), ("aj", "Ph/3,4,5-(MeO)3C6H2", 60, 77, 78, "3aaj"),
    ("ak", "Ph/1-adamantyl", 60, 87, 75, "3aak"), ("bj", "2-FC6H4/3,4,5-(MeO)3C6H2", 60, 82, 75, "3abj"),
    ("cj", "2-ClC6H4/3,4,5-(MeO)3C6H2", 60, 84, 78, "3acj"), ("dj", "2-BrC6H4/3,4,5-(MeO)3C6H2", 60, 81, 76, "3adj"),
    ("ej", "2-MeC6H4/3,4,5-(MeO)3C6H2", 60, 81, 78, "3aej"), ("fj", "3-ClC6H4/3,4,5-(MeO)3C6H2", 60, 76, 74, "3afj"),
    ("gj", "3-BrC6H4/3,4,5-(MeO)3C6H2", 60, 81, 75, "3agj"), ("hj", "3-MeC6H4/3,4,5-(MeO)3C6H2", 60, 84, 77, "3ahj"),
    ("ij", "4-FC6H4/3,4,5-(MeO)3C6H2", 60, 83, 79, "3aij"), ("jj", "4-BrC6H4/3,4,5-(MeO)3C6H2", 60, 77, 80, "3ajj"),
    ("kj", "4-IC6H4/3,4,5-(MeO)3C6H2", 60, 78, 78, "3akj"), ("lj", "4-CF3C6H4/3,4,5-(MeO)3C6H2", 60, 75, 72, "3alj"),
    ("mj", "4-i-PrC6H4/3,4,5-(MeO)3C6H2", 60, 72, 73, "3amj"), ("nj", "4-t-BuC6H4/3,4,5-(MeO)3C6H2", 60, 76, 72, "3anj"),
    ("oj", "4-BnOC6H4/3,4,5-(MeO)3C6H2", 60, 85, 73, "3aoj"), ("pj", "2,4-F2C6H3/3,4,5-(MeO)3C6H2", 60, 92, 75, "3apj"),
    ("qj", "2,4-Cl2C6H3/3,4,5-(MeO)3C6H2", 60, 93, 75, "3aqj"), ("rj", "2,5-Cl2C6H3/3,4,5-(MeO)3C6H2", 60, 77, 74, "3arj"),
    ("sj", "3,4-Me2C6H3/3,4,5-(MeO)3C6H2", 60, 87, 76, "3asj"), ("tj", "2-naphthyl/3,4,5-(MeO)3C6H2", 60, 76, 87, "3atj"),
    ("uj", "2-pyridyl/3,4,5-(MeO)3C6H2", 60, 74, 82, "3auj"),
]
for code, ar, _, _, _, product in t2:
    did = f"2{code}"
    if did not in molecules:
        smi = "COc1cc(C(=O)N=Nc2ccc3ccccc3c2)cc(OC)c1OC" if did == "2tj" else None
        add(molecule(did, did, [P[3]], f"N-acyldiazene {did}; Ar1/Ar2 = {ar}", smi,
                     "achiral" if smi else "unresolved", "unspecified", "partial" if smi else "unresolved",
                     f"Identity is retained from Table 2 substituent definition: Ar1/Ar2 = {ar}."))
    if product not in molecules:
        product_smi = None
        if product == "3atj":
            product_smi = "CCOC(=O)C=C1OC(=NN(c2ccc3ccccc3c2)C1Cc2ccccc2)c2cc(OC)c(OC)c(OC)c2"
        note = f"Product identity is retained from Table 2 substituent definition Ar1/Ar2 = {ar}; the major enantiomer configuration and exocyclic alkene geometry are not assigned."
        if product in {"3aac", "3aak"}:
            note += " Crystallographic crystals are reported as racemic, which is not taken to redefine the enantioenriched reaction sample."
        add(molecule(product, product, [P[3]], f"1,3,4-oxadiazine {product}; Ar1/Ar2 = {ar}", product_smi,
                     stereo="unspecified", geometry="unspecified", representation="not_encoded", identity_note=note))

# Table 3 allenoates and products.
t3 = [
    ("1b", "Me/Et", "CCOC(=O)C=C=CC", 60, "3btj", 76, 75),
    ("1c", "Et/Et", "CCOC(=O)C=C=CCC", 60, "3ctj", 73, 82),
    ("1d", "Ph/t-Bu", "CC(C)(C)OC(=O)C=C=Cc1ccccc1", 144, "3dtj", 72, 76),
    ("1e", "Ph/Bn", "O=C(OCc1ccccc1)C=C=Cc1ccccc1", 60, "3etj", 81, 81),
    ("1f", "Ph/CHPh2", "O=C(OC(c1ccccc1)c1ccccc1)C=C=Cc1ccccc1", 60, "3ftj", 82, 74),
]
for aid, rr, smi, _, prod, _, _ in t3:
    add(molecule(aid, f"(±)-{aid}", [P[3]], f"racemic γ-substituted allenoate {aid}; R1/R2 = {rr}", smi,
                 "racemic", "unspecified", "not_encoded",
                 "The source explicitly marks the allenoate series as (±); candidate SMILES omit axial configuration.",
                 "racemic pair; individual enantiomer ratio not reported", "allene axis; R_a/S_a not assigned"))
    add(molecule(prod, prod, [P[3]], f"1,3,4-oxadiazine {prod}; R1/R2 = {rr}, Ar1/Ar2 = 2-naphthyl/3,4,5-(MeO)3C6H2",
                 stereo="unspecified", geometry="unspecified", representation="not_encoded",
                 identity_note=f"Product identity is retained from Table 3 substituent definition R1/R2 = {rr} and the fixed 2tj substituents. The product is enantioenriched by ee but the major absolute configuration and exocyclic alkene geometry are not assigned."))

# Derivatization products.  Connectivity is encoded where the complete drawing
# is legible; stereogenic and alkene configurations remain intentionally open.
deriv_products = {
    "4aaa": "CCOC(=O)C=C1OC(=NN(c2ccccc2)C1(Br)Cc2ccccc2)c2ccccc2",
    "5aaa": "OCC=C1OC(=NN(c2ccccc2)C1Cc2ccccc2)c2ccccc2",
    "6aaa": "Cc1cc(Br)ccc1C(=O)OCC=C1OC(=NN(c2ccccc2)C1Cc2ccccc2)c2ccccc2",
    "7aaa": None,
}
for lab, smi in deriv_products.items():
    note = "Source drawing retains the stereogenic scaffold and reports ee, but no absolute configuration or exocyclic alkene E/Z descriptor is assigned."
    if lab == "7aaa":
        note += " The acid-2 polycyclic acyl fragment was not safely identified, so its product SMILES remains unresolved."
    add(molecule(lab, lab, [P[3], P[4]], f"derivatization product {lab}", smi,
                 stereo="unspecified", geometry="unspecified", representation="not_encoded",
                 identity_note=note))

# Generic source drawings used only by generic/literature records.
for mid, label, fixed, vars_, ev in [
    ("gen_allenoate", "generic γ-substituted allenoate 1", "allenoate with CO2R2 and R1 substituents", ["R1", "R2"], [P[1]]),
    ("gen_diazene", "generic N-acyldiazene 2", "N-acyldiazene with Ar1 and Ar2 substituents", ["Ar1", "Ar2"], [P[1]]),
    ("gen_product", "generic 1,3,4-oxadiazine 3", "1,3,4-oxadiazine annulation product", ["R1", "R2", "Ar1", "Ar2"], [P[1]]),
    ("lit_allenoate", "literature allenoate", "ethyl γ-benzyl allenoate drawing", [], [P[1]]),
    ("lit_diazene", "literature N-acyldiazene", "N-acyldiazene with Ar1/Ar2", ["Ar1", "Ar2"], [P[1]]),
    ("lit_product", "literature racemic 1,3,4-oxadiazine", "racemic oxadiazine product", ["Ar1", "Ar2"], [P[1]]),
]:
    add(generic_molecule(mid, label, ev, fixed, vars_))
for im in ("IM-1", "IM-2", "IM-3"):
    add(molecule(im, im, [P[4]], f"proposed intermediate {im}", stereo="unresolved", geometry="unresolved", representation="not_encoded"))

reactions = []

# Scheme 1B: cited prior work, not a current-paper experiment.
r = base_reaction("literature:scheme1B", [P[1]], "literature_example", [], "Scheme 1B",
                  "organocatalytic [2 + 4]-annulation (cited prior work)", "annulation")
r["participants"] = [participant("lit_1", "lit_allenoate", "reactant", [P[1]]),
                     participant("lit_2", "lit_diazene", "reactant", [P[1]]),
                     participant("lit_cat", "DMAP", "catalyst", [P[1]], mol_percent=Q(20, "mol%"), ref="lit_diazene"),
                     participant("lit_sol", "toluene", "solvent", [P[1]]),
                     participant("lit_3", "lit_product", "product", [P[1]], product_role="racemic product")]
r["solvent_conditions"]["composition"] = F("toluene", [f"{DOC}:mol:toluene"], [P[1]])
step(r, "reaction", [P[1]], ["lit_1", "lit_2", "lit_cat", "lit_sol", "lit_3"],
     "DMAP-catalyzed annulation", ("rt", "rt"), Q(48, "h"))
outcome(r, "yield_range", "yield", "lit_product", [P[1]], "0-80% yields",
        {"minimum": 0, "maximum": 80, "unit": "%", "raw_text": "0-80% yields"}, scope="34 examples")
r["review"]["unresolved"].append("Only the cited summary in Scheme 1B was provided; the cited paper's experimental details were not supplied.")
reactions.append(r)

# Scheme 1C current-paper generic scheme.
r = base_reaction("generic:scheme1C", [P[1]], "generic_scheme", [], "Scheme 1C (this work)")
r["participants"] = [participant("g1", "gen_allenoate", "reactant", [P[1]]),
                     participant("g2", "gen_diazene", "reactant", [P[1]]),
                     participant("gcat", "C5", "catalyst", [P[1]], mol_percent=Q(15, "mol%"), ref="gen_diazene"),
                     participant("gsol", "MTBE", "solvent", [P[1]]),
                     participant("g3", "gen_product", "product", [P[1]], product_role="enantioenriched product")]
r["solvent_conditions"]["composition"] = F("MTBE", [f"{DOC}:mol:MTBE"], [P[1]])
step(r, "reaction", [P[1]], ["g1", "g2", "gcat", "gsol", "g3"], "[2 + 4]-annulation",
     Q(30, "°C"), ("48-144 h", {"minimum": 48, "maximum": 144, "unit": "h", "raw_text": "48-144 h"}))
outcome(r, "yield_range", "yield", "gen_product", [P[1]], "66-93% yield",
        {"minimum": 66, "maximum": 93, "unit": "%", "raw_text": "66-93% yield"}, scope="36 examples")
outcome(r, "ee_range", "selectivity", "gen_product", [P[1]], "70-87% ee",
        {"minimum": 70, "maximum": 87, "unit": "%", "raw_text": "70-87% ee"}, scope="ee; major absolute configuration unassigned")
reactions.append(r)

# Table 1 optimization rows.
opt = [
    (1,"C1","toluene",2.0,"rt",48,"trace",None),(2,"C2","toluene",2.0,"rt",48,32,29),
    (3,"C3","toluene",2.0,"rt",48,23,68),(4,"C4","toluene",2.0,"rt",48,6,35),
    (5,"C5","toluene",2.0,"rt",48,32,75),(6,"C6","toluene",2.0,"rt",48,29,68),
    (7,"C5","EtOAc",2.0,"rt",48,66,69),(8,"C5","MeCN",2.0,"rt",48,59,47),
    (9,"C5","CH2Cl2",2.0,"rt",48,35,70),(10,"C5","CHCl3",2.0,"rt",48,21,71),
    (11,"C5","Et2O",2.0,"rt",48,39,70),(12,"C5","THF",2.0,"rt",48,57,69),
    (13,"C5","dioxane",2.0,"rt",48,66,73),(14,"C5","MTBE",2.0,"rt",48,67,74),
    (15,"C5","MTBE",1.5,"rt",48,65,75),(16,"C5","MTBE",1.0,"rt",48,57,76),
    (17,"C5","MTBE",3.0,"rt",48,77,71),(18,"C5","MTBE",1.5,"rt",48,77,71),
    (19,"C5","MTBE",1.5,"rt",48,70,72),(20,"C5","MTBE",1.5,"rt",48,79,69),
    (21,"C5","MTBE",1.5,30,48,70,74),(22,"C5","MTBE",1.5,35,48,70,71),
    (23,"C5","MTBE",1.5,30,42,74,72),(24,"C5","MTBE",1.5,30,48,82,71),
    (25,"C5","MTBE",1.5,30,54,74,71),
]
for ent, cat, solv, vol, temp, hours, yld, ee in opt:
    ev = [P[2], P[3]] if ent == 24 else [P[2]]
    r = base_reaction(f"optimization:table1_entry{ent:02d}", ev, "specific_experiment", ["optimization"], f"Table 1 entry {ent}")
    a1, a2, load = 0.3, 0.2, 10
    if ent == 18: a1 = 0.4
    if ent == 19: a1, a2 = 0.2, 0.3
    if ent == 20: a1, a2 = 0.3, 0.4
    if ent >= 23: load = 15
    pa = participant("allenoate", "1a", "reactant", ev, amount=Q(a1, "mmol"))
    pd = participant("diazene", "2aa", "reactant", ev, amount=Q(a2, "mmol"))
    if ent == 24:
        pa["charges"]["amount"] = F(["Table 1/Table 2: 0.3 mmol", "optimal-conditions prose: 0.2 mmol"], None, ev, status="conflict")
        pa["charges"]["amount"]["alternatives"] = [F("Table 1 footnote and Table 2 footnote: (±)-1a (0.3 mmol)", {"number":0.3,"unit":"mmol"}, [P[2],P[3]]), F("prose: racemic γ-benzyl allenoate 1a (0.2 mmol)", {"number":0.2,"unit":"mmol"}, [P[2]])]
        pd["charges"]["amount"] = F(["Table 1/Table 2: 0.2 mmol", "optimal-conditions prose: 0.3 mmol"], None, ev, status="conflict")
        pd["charges"]["amount"]["alternatives"] = [F("Table 1 footnote and Table 2 footnote: 2aa (0.2 mmol)", {"number":0.2,"unit":"mmol"}, [P[2],P[3]]), F("prose: 2aa (0.3 mmol)", {"number":0.3,"unit":"mmol"}, [P[2]])]
        r["review"]["unresolved"].append("Source conflict: Table 1 entry 24/Table 2 entry 1 use 1a 0.3 mmol and 2aa 0.2 mmol, while the optimal-conditions prose reverses these amounts.")
        r["notes"].append(F(None, "Table/prose reactant-charge conflict retained without choosing a value.", ev, "inference", "Compared the Table 1 and Table 2 footnotes with the adjacent optimal-conditions prose."))
    r["participants"] = [pa, pd,
                         participant("catalyst", cat, "catalyst", ev, mol_percent=Q(load,"mol%"), ref="2aa"),
                         participant("solvent", solv, "solvent", ev, volume=Q(vol,"mL")),
                         participant("product", "3aaa", "product", ev, product_role="major product")]
    r["reaction_scale"] = F(f"2aa {a2} mmol" if ent != 24 else "conflicting 1a/2aa scale statements",
                             {"reference_molecule_id": f"{DOC}:mol:2aa", "amount": {"number": a2, "unit": "mmol"}} if ent != 24 else None,
                             ev, status="conflict" if ent == 24 else None)
    if ent == 24:
        r["reaction_scale"]["alternatives"] = [F("Table basis: 2aa 0.2 mmol", {"reference_molecule_id":f"{DOC}:mol:2aa","amount":{"number":0.2,"unit":"mmol"}}, [P[2],P[3]]), F("prose basis: 1a 0.2 mmol", {"reference_molecule_id":f"{DOC}:mol:1a","amount":{"number":0.2,"unit":"mmol"}}, [P[2]])]
    r["solvent_conditions"]["composition"] = F(solv, [f"{DOC}:mol:{solv}"], ev)
    r["solvent_conditions"]["total_volume"] = F(*Q(vol,"mL"), ev)
    t = ("rt", "rt") if temp == "rt" else Q(temp, "°C")
    step(r, "reaction", ev, ["allenoate","diazene","catalyst","solvent","product"], "annulation", t, Q(hours,"h"))
    if yld == "trace":
        outcome(r, "trace_product", "trace", "3aaa", ev, "trace", "trace", basis="isolated yield column; trace reported")
    else:
        outcome(r, "isolated_yield", "yield", "3aaa", ev, f"{yld}%", {"number":yld,"unit":"%"}, "isolated", "isolated yield")
    if ee is not None:
        outcome(r, "ee", "selectivity", "3aaa", ev, f"{ee}% ee", {"number":ee,"unit":"%"}, basis="chiral HPLC analysis", scope="ee; major enantiomer absolute configuration unassigned")
    if ent == 24:
        r["procedure"]["addition_order"] = F("2aa was added to a solution of 1a with C5", "2aa added to a solution of 1a and C5", [P[2]])
    reactions.append(r)

# Table 2: entry 1 is the same experiment as Table 1 entry 24 and is already
# merged into that record through page-3 evidence.
for idx, (code, ar, hours, yld, ee, product) in enumerate(t2[1:], 2):
    did = f"2{code}"
    r = base_reaction(f"scope:table2_entry{idx:02d}", [P[3]], "specific_experiment", ["substrate_scope"], f"Table 2 entry {idx}")
    r["general_scheme_ids"] = [f"{DOC}:generic:scheme1C"]
    r["participants"] = [participant("allenoate", "1a", "reactant", [P[3]], amount=Q(0.3,"mmol")),
                         participant("diazene", did, "reactant", [P[3]], amount=Q(0.2,"mmol")),
                         participant("catalyst", "C5", "catalyst", [P[3]], mol_percent=Q(15,"mol%"), ref=did),
                         participant("solvent", "MTBE", "solvent", [P[3]], volume=Q(1.5,"mL")),
                         participant("product", product, "product", [P[3]], product_role="major product")]
    r["reaction_scale"] = F("diazene 2 (0.2 mmol)", {"reference_molecule_id":f"{DOC}:mol:{did}","amount":{"number":0.2,"unit":"mmol"}}, [P[3]])
    r["solvent_conditions"]["composition"] = F("MTBE", [f"{DOC}:mol:MTBE"], [P[3]])
    r["solvent_conditions"]["total_volume"] = F(*Q(1.5,"mL"), [P[3]])
    step(r, "reaction", [P[3]], ["allenoate","diazene","catalyst","solvent","product"], "annulation", Q(30,"°C"), Q(hours,"h"))
    outcome(r, "isolated_yield", "yield", product, [P[3]], f"{yld}%", {"number":yld,"unit":"%"}, "isolated", "isolated yield")
    outcome(r, "ee", "selectivity", product, [P[3]], f"{ee}% ee", {"number":ee,"unit":"%"}, basis="chiral HPLC analysis", scope="ee; major absolute configuration unassigned")
    r["notes"].append(F(f"Ar1/Ar2 = {ar}", {"Ar1/Ar2": ar}, [P[3]]))
    reactions.append(r)

# Table 3 scope rows.
for idx, (aid, rr, _, hours, prod, yld, ee) in enumerate(t3, 1):
    r = base_reaction(f"scope:table3_entry{idx:02d}", [P[3]], "specific_experiment", ["substrate_scope"], f"Table 3 entry {idx}")
    r["general_scheme_ids"] = [f"{DOC}:generic:scheme1C"]
    r["participants"] = [participant("allenoate", aid, "reactant", [P[3]], amount=Q(0.3,"mmol")),
                         participant("diazene", "2tj", "reactant", [P[3]], amount=Q(0.2,"mmol")),
                         participant("catalyst", "C5", "catalyst", [P[3]], mol_percent=Q(15,"mol%"), ref="2tj"),
                         participant("solvent", "MTBE", "solvent", [P[3]], volume=Q(1.5,"mL")),
                         participant("product", prod, "product", [P[3]], product_role="major product")]
    r["reaction_scale"] = F("2tj (0.2 mmol)", {"reference_molecule_id":f"{DOC}:mol:2tj","amount":{"number":0.2,"unit":"mmol"}}, [P[3]])
    r["solvent_conditions"]["composition"] = F("MTBE", [f"{DOC}:mol:MTBE"], [P[3]])
    r["solvent_conditions"]["total_volume"] = F(*Q(1.5,"mL"), [P[3]])
    step(r, "reaction", [P[3]], ["allenoate","diazene","catalyst","solvent","product"], "annulation", Q(30,"°C"), Q(hours,"h"))
    outcome(r, "isolated_yield", "yield", prod, [P[3]], f"{yld}%", {"number":yld,"unit":"%"}, "isolated", "isolated yield")
    outcome(r, "ee", "selectivity", prod, [P[3]], f"{ee}% ee", {"number":ee,"unit":"%"}, basis="chiral HPLC analysis", scope="ee; major absolute configuration unassigned")
    r["notes"].append(F(f"R1/R2 = {rr}", {"R1/R2":rr}, [P[3]]))
    reactions.append(r)

# Scale-up experiment.
r = base_reaction("scale_up:scheme2A", [P[3],P[4]], "specific_experiment", ["scale_up"], "Scheme 2A")
r["general_scheme_ids"] = [f"{DOC}:generic:scheme1C"]
r["participants"] = [participant("allenoate","1a","reactant",[P[3],P[4]],amount=Q(1.5,"mmol")),
                     participant("diazene","2tj","reactant",[P[3],P[4]],amount=Q(1.0,"mmol")),
                     participant("catalyst","C5","catalyst",[P[4]],mol_percent=Q(15,"mol%"),ref="2tj"),
                     participant("solvent","MTBE","solvent",[P[4]],volume=Q(7.5,"mL")),
                     participant("product","3atj","product",[P[3],P[4]],product_role="major product")]
r["reaction_scale"] = F("2tj 1.0 mmol", {"reference_molecule_id":f"{DOC}:mol:2tj","amount":{"number":1.0,"unit":"mmol"}}, [P[3],P[4]])
r["solvent_conditions"]["composition"] = F("MTBE", [f"{DOC}:mol:MTBE"], [P[4]])
r["solvent_conditions"]["total_volume"] = F(*Q(7.5,"mL"), [P[4]])
step(r,"reaction",[P[4]],["allenoate","diazene","catalyst","solvent","product"],"scale-up annulation",Q(30,"°C"),Q(60,"h"))
outcome(r,"product_amount","product_amount","3atj",[P[4]],"421 mg",{"number":421,"unit":"mg"})
outcome(r,"yield","yield","3atj",[P[3],P[4]],"76% yield",{"number":76,"unit":"%"})
outcome(r,"ee","selectivity","3atj",[P[3],P[4]],"87% ee",{"number":87,"unit":"%"},scope="ee; major absolute configuration unassigned")
r["review"]["unresolved"].append("Scheme 2A does not state the yield or ee measurement basis; table footnotes were not silently transferred.")
reactions.append(r)

# Scheme 2B derivatizations, one record per arrow.
deriv = [
    ("bromination", "3aaa", "4aaa", [("Br2","reagent"),("Et3N","reagent"),("CH2Cl2","solvent")], "rt", 62, (("C","Br","formed"),("C","H","hydrogen_loss"))),
    ("ester_reduction", "4aaa", "5aaa", [("DIBAL-H","reagent"),("CH2Cl2","solvent")], -78, 88,
     (("C","H","hydrogen_gain",0,1,2),("O","H","hydrogen_gain",0,1,1),
      ("C","O","broken",1,0),("C","O","order_decreased",2,1))),
    ("esterification_acid1", "5aaa", "6aaa", [("acid1","reagent"),("DCC","reagent"),("DMAP","catalyst"),("CH2Cl2","solvent")], None, 78,
     (("C","O","formed"),("O","H","hydrogen_loss"))),
    ("esterification_acid2", "5aaa", "7aaa", [("acid2","reagent"),("DCC","reagent"),("DMAP","catalyst"),("CH2Cl2","solvent")], None, 68,
     (("C","O","formed"),("O","H","hydrogen_loss"))),
]
for name, start, prod, agents, temp, yld, pairs in deriv:
    r = base_reaction(f"derivatization:{name}", [P[3],P[4]], "specific_experiment", ["derivatization"], f"Scheme 2B: {start}→{prod}", name.replace("_"," "), name, pairs)
    parts = [participant("substrate",start,"reactant",[P[4]]), participant("product",prod,"product",[P[4]],product_role="product")]
    for i,(mid,role) in enumerate(agents,1): parts.append(participant(f"agent{i}",mid,role,[P[4]]))
    r["participants"] = parts
    solvent = next((mid for mid,role in agents if role=="solvent"),None)
    if solvent: r["solvent_conditions"]["composition"] = F(solvent,[f"{DOC}:mol:{solvent}"],[P[4]])
    tv = None if temp is None else (("rt","rt") if temp=="rt" else Q(temp,"°C"))
    step(r,"reaction",[P[4]],[p["participant_id"] for p in parts],name.replace("_"," "),tv,None)
    outcome(r,"yield","yield",prod,[P[3],P[4]],f"{yld}% yield",{"number":yld,"unit":"%"})
    outcome(r,"ee","selectivity",prod,[P[3],P[4]],"76% ee",{"number":76,"unit":"%"},scope="ee; major absolute configuration unassigned")
    r["review"]["unresolved"].append("Scheme 2B does not report reagent amounts, reaction time, workup, purification, or measurement basis.")
    if temp is None: r["review"]["unresolved"].append("Reaction temperature is not shown for this esterification arrow.")
    reactions.append(r)

# Scheme 3 is an author-proposed mechanism, not an experimentally verified row.
r = base_reaction("mechanism:scheme3", [P[3],P[4]], "proposed_mechanism", [], "Scheme 3",
                  "proposed mechanism for [2 + 4]-annulation", "proposed catalytic mechanism")
r["participants"] = [participant("allenoate","1a","reactant",[P[4]]), participant("diazene","2aa","reactant",[P[4]]),
                     participant("catalyst","C5","catalyst",[P[4]]),
                     participant("im1","IM-1","unknown",[P[4]]), participant("im2","IM-2","unknown",[P[4]]),
                     participant("im3","IM-3","unknown",[P[4]]), participant("product","3aaa","product",[P[4]])]
r["procedure"]["original_text"] = F("initial addition ... IM-1; conjugate addition ... IM-2; intramolecular 1,4-addition ... IM-3; elimination ... 3aaa and regeneration of C5",
                                      "C5 addition to 1a → IM-1; addition to 2aa → IM-2; intramolecular 1,4-addition → IM-3; elimination → 3aaa + C5", [P[4]])
r["review"]["unresolved"].append("This is explicitly a probable/proposed mechanism; intermediates were not reported as isolated or directly observed here.")
reactions.append(r)

bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id":f"{DOC}:context:scheme1A","document_id":DOC,"evidence_ids":[P[1],I[1]],"kind":"non_synthetic_context","note":"Scheme 1A shows selected bioactive oxadiazines; these are contextual structures, not experiments in this paper."},
    {"context_id":f"{DOC}:context:crystallization","document_id":DOC,"evidence_ids":[P[3],I[3]],"kind":"stereochemical_observation","note":"Crystals of 3aac and 3aak contained racemic pairs; residual 3aac mother solution was 99% ee. This was not converted into a synthetic reaction row or an absolute-configuration assignment."},
    {"context_id":f"{DOC}:context:references","document_id":DOC,"evidence_ids":[P[4],I[4],P[5],I[5]],"kind":"bibliography","note":"Reference lists and acknowledgements contain no additional current-paper synthetic experiment on pages 4-5."},
]

inventory = {
    "document_id": DOC,
    "pages_reviewed": [1,2,3,4,5],
    "visual_assets": [f"assets/{DOC}/page{i}.jpg" for i in range(1,6)],
    "table_rows": {"table1":25,"table2":31,"table3":5},
    "deduplication": {"table2_entry1": f"{DOC}:optimization:table1_entry24"},
    "other_arrows": {"scale_up":1,"derivatization":4,"proposed_mechanism":1,"current_generic_scheme":1,"cited_literature_scheme":1},
    "expected_reaction_records": 68,
    "notes": ["Page 5 contains references only.", "No supplementary file was supplied or associated.", "Table 1 entry 24/prose charge conflict is retained."],
}
RAW.write_text(json.dumps(inventory, ensure_ascii=False, indent=2) + "\n")
diagnostics = {}
for key in ("1a", "1b", "1c", "1d", "1e", "1f", "2aa", "2tj", "3aaa", "3atj", "4aaa", "5aaa", "6aaa", "acid1"):
    smi = st.value(molecules[key]["smiles"])
    diagnostics[key] = st.stereo_inspect(smi)
DIAG.write_text(json.dumps({
    "document_id": DOC,
    "source_verified": False,
    "note": "Graph diagnostics only; allene axial configuration and unreported E/Z cannot be assigned by toolkit output.",
    "molecules": diagnostics,
}, ensure_ascii=False, indent=2) + "\n")
regions = {
    1: [
        {"region":"title, abstract, and introductory prose","finding":"Article metadata and background; no current-paper concrete experiment in this region.","evidence_ids":[P[1]]},
        {"region":"Scheme 1A","finding":"Bioactive oxadiazine context only; retained as non-synthetic source_context.","evidence_ids":[P[1]]},
        {"region":"Scheme 1B","finding":"Cited prior DMAP-catalyzed racemic annulation; classified as literature_example.","evidence_ids":[P[1]]},
        {"region":"Scheme 1C","finding":"This-work generic asymmetric annulation; classified as generic_scheme, not a concrete experiment.","evidence_ids":[P[1]]},
        {"region":"ESI footnote","finding":"The paper says ESI is available, but no SI file was supplied; association remains unestablished.","evidence_ids":[P[1]]},
    ],
    2: [
        {"region":"Table 1 reaction drawing and catalysts C1-C6","finding":"Model substrates/products and catalyst structures inventoried; catalyst absolute descriptors are not written.","evidence_ids":[P[2]]},
        {"region":"Table 1 entries 1-25 and footnotes a-h","finding":"Twenty-five distinct optimization experiments extracted row-by-row, including trace and all loading/ratio/volume exceptions.","evidence_ids":[P[2]]},
        {"region":"optimal-conditions prose","finding":"Same experiment as Table 1 entry 24; reactant amounts conflict with table footnotes and were retained as alternatives.","evidence_ids":[P[2]]},
        {"region":"Table 2 discussion continuation","finding":"Narrative mentions scope products but the definitive rows are on page 3.","evidence_ids":[P[2]]},
    ],
    3: [
        {"region":"Table 2 entries 1-31 and footnotes","finding":"Entry 1 duplicates Table 1 entry 24 and was merged; entries 2-31 are separate substrate-scope records.","evidence_ids":[P[3]]},
        {"region":"Table 3 entries 1-5 and footnotes","finding":"Five allenoate-scope experiments extracted row-by-row.","evidence_ids":[P[3]]},
        {"region":"crystallization paragraph","finding":"Racemic crystals of 3aac/3aak and 99% ee mother liquor for 3aac retained as stereochemical source_context, not synthesis rows.","evidence_ids":[P[3]]},
        {"region":"scale-up and transformation prose","finding":"Scale-up and four Scheme 2 transformations identified; figure conditions are on page 4.","evidence_ids":[P[3]]},
    ],
    4: [
        {"region":"Scheme 2A","finding":"One scale-up annulation (1a + 2tj to 3atj) extracted with charges, solvent volume, temperature, time, mass, yield, and ee.","evidence_ids":[P[4]]},
        {"region":"Scheme 2B","finding":"Four separate derivatization arrows extracted: bromination, ester reduction, and two esterifications.","evidence_ids":[P[4]]},
        {"region":"Scheme 3 and mechanism prose","finding":"Author-proposed catalytic cycle classified as proposed_mechanism; intermediates not treated as isolated experiments.","evidence_ids":[P[4]]},
        {"region":"conclusion, acknowledgements, references 1-5","finding":"No additional current-paper synthetic experiment.","evidence_ids":[P[4]]},
    ],
    5: [
        {"region":"references 6-14","finding":"Bibliography only; no current-paper experiment, table, or reaction arrow.","evidence_ids":[P[5]]},
    ],
}
for page_regions in regions.values():
    for region in page_regions:
        page = int(region["evidence_ids"][0].split(":page", 1)[1].split(":", 1)[0])
        if I[page] not in region["evidence_ids"]:
            region["evidence_ids"].append(I[page])
page_unresolved = {
    1:["Generic R-group drawings and cited-work experimental details are not concrete structures/procedures.","ESI is mentioned but was not supplied."],
    2:["Table 1 entry 24 and optimal-conditions prose reverse the 1a/2aa amounts.","Catalyst C1-C6 full stereochemical identities are not assigned by textual descriptors."],
    3:["Scope-product absolute configurations are unassigned despite ee values.","Crystallization observations do not establish the major reaction-product absolute configuration."],
    4:["Scheme 2B omits reagent quantities and reaction times; esterification temperatures are not shown.","Scheme 2 outcomes do not state measurement bases.","Scheme 3 is proposed rather than experimentally verified."],
    5:["No reaction content; cited papers were not imported as source documents."],
}
COVERAGE.write_text(json.dumps({
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [{
        "page": page,
        "state": "attempted",
        "regions": regions[page],
        "record_ids": [r["reaction_id"] for r in reactions if P[page] in r["evidence_ids"]],
        "unresolved": page_unresolved[page],
    } for page in range(1,6)],
}, ensure_ascii=False, indent=2) + "\n")
OUT.write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"molecules":len(bundle["molecules"]),"reactions":len(reactions),"out":str(OUT)},ensure_ascii=False))
