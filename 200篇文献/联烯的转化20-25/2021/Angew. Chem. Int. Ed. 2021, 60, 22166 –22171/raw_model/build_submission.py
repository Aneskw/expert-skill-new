#!/usr/bin/env python3
"""Build the source-specific v3.1 submission for doc_7e7fe8ed617d0044.

This script reads only the clean standard-init template and source-derived constants
transcribed from PDF pages 1-6.  It does not read earlier reaction records.
"""
from __future__ import annotations

import json
import sys
from copy import deepcopy
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
RAW = RUN / "raw_model"
DOC = "doc_7e7fe8ed617d0044"
EVT = {p: f"{DOC}:page{p}:text" for p in range(1, 7)}
EVI = {p: f"{DOC}:page{p}:image" for p in range(1, 7)}


def dump(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def ev(page):
    return [EVT[page], EVI[page]]


def f(raw=None, norm=None, evidence=(), method="original_extraction", derivation=None,
      status=None, alternatives=()):
    return st.fact(raw, norm, evidence, method=method, derivation=derivation,
                   status=status, alternatives=alternatives)


def q(number, unit, raw=None, qualifier=None):
    d = {"number": number, "unit": unit}
    if raw is not None:
        d["raw_text"] = raw
    if qualifier:
        d["qualifier"] = qualifier
    return d


def qr(minimum, maximum, unit, raw):
    return {"minimum": minimum, "maximum": maximum, "unit": unit, "raw_text": raw}


bundle = json.loads((RAW / "standard_init_titled.json").read_text())
molecules = {}
stereo_checks = {}


def molecule(mid, label, page, name=None, smiles=None, chirality="unresolved",
             geometry="unresolved", rep_status=None, assignment=None, axial=None,
             mixture=None, note=None):
    evidence = ev(page)
    m = st.new_molecule(f"{DOC}:mol:{mid}", DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    m["names"]["original"] = f(name or label, name or label, evidence)
    m["names"]["english"] = f(name or label, name or label, evidence)
    m["cas_numbers"] = f()
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = f({"source": "structure drawing", "label": label}, smiles, evidence,
                          method="structure_derivation",
                          derivation="Candidate SMILES transcribed conservatively from the displayed structure; parsing does not verify identity.")
        try:
            diag = st.stereo_inspect(smiles)
            stereo_checks[m["molecule_id"]] = diag
            if chirality == "achiral" and not any(x["specified"] == "Unspecified" for x in diag["potential_stereo"]):
                m["stereochemistry"]["isomeric_smiles"] = f(
                    None, diag["canonical_isomeric_smiles"], evidence,
                    method="structure_derivation", derivation="RDKit canonical isomeric SMILES of the candidate source graph.")
                m["stereochemistry"]["identity_key"] = f(
                    None, diag["stereo_graph_key"], evidence, method="structure_derivation",
                    derivation="Stereo graph key computed from the candidate SMILES; not a sample identity.")
        except Exception as exc:
            stereo_checks[m["molecule_id"]] = {"error": str(exc), "source_verified": False}
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = f()
    m["stereochemistry"]["chirality_status"] = f(
        None, chirality, evidence, method="inference",
        derivation="Conservative status from the displayed graph and source stereochemical annotations; absence of a descriptor is not treated as achiral or racemic.")
    m["stereochemistry"]["geometry_status"] = f(
        None, geometry, evidence, method="inference",
        derivation="Conservative status from the displayed graph and source stereochemical annotations.")
    m["stereochemistry"]["representation_status"] = f(
        None, rep_status or ("partial" if smiles else "not_encoded"), evidence,
        method="inference", derivation="Status reflects whether the candidate graph encodes all source stereochemical information.")
    if assignment is not None:
        m["stereochemistry"]["assignment"] = f(assignment, assignment, evidence)
    if axial is not None:
        m["stereochemistry"]["axial"] = f(axial, axial, evidence, method="structure_derivation",
                                            derivation="Allene-axis assessment from the displayed/candidate local graph; toolkit support is limited.")
    if mixture is not None:
        m["stereochemistry"]["mixture"] = f(mixture, mixture, evidence)
    m["identity_note"] = f(note, note, evidence) if note else f()
    molecules[mid] = m
    return m["molecule_id"]


# Generic/source-family identities.
for mid, label, name in [
    ("prior_acid", "prior-work beta-keto acid family", "generic beta-keto acid family"),
    ("prior_allene", "prior-work carbon-substituted allene family", "generic carbon-substituted allene family"),
    ("prior_product", "prior-work racemic addition product family", "generic racemic decarboxylative-addition products"),
    ("precedent_allene", "alkoxyallene family", "generic alkoxyallene family"),
    ("precedent_nu", "H-Nu", "generic heteroatom nucleophile"),
    ("precedent_product", "asymmetric hydrofunctionalization product family", "generic heteroatom-substituted allyl product family"),
    ("current_acid", "beta-keto acid family", "generic beta-keto acid family"),
    ("current_allene", "heteroatom-substituted allene family", "generic alkoxy-/amidoallene family"),
    ("current_product", "beta-heteroatom ketone family", "generic branched gamma,delta-unsaturated ketone family"),
]:
    molecule(mid, label, 2, name=name)

# Concrete recurring substrates and products.
molecule("acid1", "1", 2, "3-oxo-3-phenylpropanoic acid", "O=C(O)CC(=O)c1ccccc1", "achiral", "not_applicable")
molecule("allene2a", "2a", 2, "cyclohexyloxyallene", "C=C=COC1CCCCC1", "achiral", "unspecified", "partial",
         axial=[{"axis": "C=C=C", "assignment": "not stereogenic on candidate graph because one terminal carbon is CH2"}])
molecule("allene2b", "2b", 2, "benzyloxyallene", "C=C=COCc1ccccc1", "achiral", "unspecified", "partial",
         axial=[{"axis": "C=C=C", "assignment": "not stereogenic on candidate graph because one terminal carbon is CH2"}])
molecule("prod3a", "3a", 2, "beta-cyclohexyloxy ketone 3a", "C=CC(OC1CCCCC1)CC(=O)c1ccccc1",
         "specified", "unspecified", "partial", assignment={"center": "beta carbon", "configuration": "R", "source": "Table 1 drawing"})
molecule("prod3b", "3b", 2, "beta-benzyloxy ketone 3b", "C=CC(OCc1ccccc1)CC(=O)c1ccccc1",
         "specified", "unspecified", "partial", assignment={"center": "beta carbon", "configuration": "R", "source": "Table 1 drawing; main text refers to SI determination"})

acid_scope = {
    "acid4": ("beta-keto acid precursor to 4", "O=C(O)CC(=O)c1ccc(OC)cc1"),
    "acid5": ("beta-keto acid precursor to 5", "O=C(O)CC(=O)c1ccc(Br)cc1"),
    "acid6": ("beta-keto acid precursor to 6", "O=C(O)CC(=O)c1ccc2ccccc2c1"),
    "acid7": ("beta-keto acid precursor to 7", "O=C(O)CC(=O)c1ccsc1"),
}
for mid, (name, smi) in acid_scope.items():
    molecule(mid, mid.replace("acid", ""), 3, name, smi, "achiral", "not_applicable")
for n in range(4, 8):
    molecule(f"prod{n}", str(n), 3, f"beta-cyclohexyloxy ketone {n}", None,
             "partially_specified", "unspecified", "not_encoded",
             assignment={"center": "beta carbon", "configuration": "source-depicted wedge; absolute R/S not printed beside this product"},
             mixture={"type": "enantioenriched sample", "composition": "experiment-specific ee"})

molecule("allene8", "8", 3, "N-allenyl-2-oxazolidinone", "C=C=CN1CCOC1=O", "achiral", "unspecified", "partial",
         axial=[{"axis": "C=C=C", "assignment": "not stereogenic on candidate graph because one terminal carbon is CH2"}])
for n in range(9, 17):
    assignment = {"center": "beta carbon", "configuration": "R", "source": "explicit label"} if n in (9, 14) else {"center": "beta carbon", "configuration": "source-depicted wedge; absolute R/S not printed beside this product"}
    molecule(f"prod{n}", str(n), 3, f"beta-amido ketone {n}",
             "C=CC(N1CCOC1=O)CC(=O)c1ccccc1" if n == 9 else None,
             "specified" if n in (9, 14) else "partially_specified", "unspecified", "partial" if n == 9 else "not_encoded",
             assignment=assignment, mixture={"type": "enantioenriched sample", "composition": "experiment-specific ee"})
for n in range(10, 17):
    molecule(f"acid{n}", f"acid for {n}", 3, f"beta-keto acid precursor to product {n}")

molecule("allene17", "17", 3, "acyclic N-benzyl-N-Cbz amidoallene", "C=C=CN(Cc1ccccc1)C(=O)OCc1ccccc1",
         "achiral", "unspecified", "partial", axial=[{"axis": "C=C=C", "assignment": "not stereogenic on candidate graph because one terminal carbon is CH2"}])
for n in range(18, 23):
    molecule(f"prod{n}", str(n), 3, f"branched beta-amido ketone {n}", None, "partially_specified", "unspecified", "not_encoded",
             assignment={"center": "beta carbon", "configuration": "source-depicted wedge; absolute R/S not printed"},
             mixture={"type": "branched/linear and enantiomeric composition", "composition": "experiment-specific b:l and ee"})
    molecule(f"acid{n}", f"acid for {n}", 3, f"beta-keto acid precursor to product {n}")

molecule("malonate23", "23", 3, "tert-butyl hydrogen malonate", "CC(C)(C)OC(=O)CC(=O)O", "achiral", "not_applicable")
molecule("desired24", "24", 3, "desired beta-alkoxy tert-butyl ester 24")
molecule("malonate25", "25", 3, "phenyl hydrogen malonate", "O=C(O)CC(=O)Oc1ccccc1", "achiral", "not_applicable")
molecule("prod26", "26", 3, "beta-alkoxy phenyl ester 26", None, "partially_specified", "unspecified", "not_encoded",
         assignment={"center": "beta carbon", "configuration": "source-depicted wedge"})
molecule("prod27", "27", 3, "beta-amido phenyl ester 27", None, "partially_specified", "unspecified", "not_encoded",
         assignment={"center": "beta carbon", "configuration": "source-depicted wedge"}, mixture={"type": "branched/linear mixture", "ratio": "14:1"})
for n in (28, 30):
    molecule(f"allene{n}", str(n), 3, f"ene-alkoxyallene {n}", None, "unspecified", "unspecified", "not_encoded",
             axial=[{"axis": "allene", "assignment": "source gives no axial descriptor"}])
for n, ring in ((29, "dihydrofuran"), (31, "dihydropyran")):
    molecule(f"prod{n}", str(n), 3, f"enantioenriched {ring} {n}", None, "partially_specified", "unspecified", "not_encoded",
             assignment={"center": "ring stereocenter", "configuration": "source-depicted wedge; absolute R/S not printed"},
             mixture={"type": "enantioenriched sample", "composition": "93% ee"})
    molecule(f"intermediate{n}", f"unnamed pre-RCM intermediate to {n}", 3, f"unnamed hydrofunctionalization intermediate leading to {n}")

molecule("ester32", "32", 4, "methyl 3-oxo-3-phenylpropanoate", "COC(=O)CC(=O)c1ccccc1", "achiral", "not_applicable")
for n in (33, 34, 35):
    molecule(f"species{n}", str(n), 4, f"mechanistic species {n}")
molecule("E36", "(E)-36", 4, "(E)-36", None, "unspecified", "specified", "not_encoded",
         assignment={"alkene": "E", "source": "Scheme 5 label"})
molecule("Z36", "(Z)-36", 4, "(Z)-36", None, "unspecified", "specified", "not_encoded",
         assignment={"alkene": "Z", "source": "Scheme 5 label"})
molecule("species37", "37", 4, "compound 37", None, "unresolved", "unresolved", "not_encoded",
         note="Scheme 5B's observed drawing and the crossover drawing labelled 37 are difficult to reconcile confidently; identity is retained unresolved rather than silently merged to a guessed graph.")
for n in (38, 39, 40, 41):
    molecule(f"prod{n}", str(n), 4, f"hydrocarboxylation adduct {n}", None, "unspecified", "unspecified", "not_encoded")
molecule("benzoic", "benzoic acid", 4, "benzoic acid", "O=C(O)c1ccccc1", "achiral", "not_applicable")
molecule("amido_benzoic_products", "unidentified analogous benzoic-acid/amidoallene adduct pattern", 5,
         "unidentified products from benzoic acid and amidoallene 8")

# Catalysts, ligands, solvents and additives. Complex structures are deliberately unresolved.
for mid, label, page in [
    ("pd2dba3", "Pd2(dba)3", 2), ("L1", "(R,R)-L1", 2), ("L2", "(R,R)-L2", 2),
    ("L3", "(R,R)-L3", 2), ("L4", "L4", 3), ("Et3N", "Et3N", 2),
    ("TFA", "CF3CO2H", 3), ("RuCat", "2nd-generation Hoveyda-Grubbs catalyst", 3),
]:
    molecule(mid, label, page, label)
for mid, label, smi, page in [
    ("DCM", "CH2Cl2", "ClCCl", 2), ("THF", "THF", "C1CCOC1", 2),
    ("MeCN", "CH3CN", "CC#N", 2), ("toluene", "toluene", "Cc1ccccc1", 3),
    ("CD2Cl2", "CD2Cl2", "[2H]C([2H])(Cl)Cl", 4),
]:
    molecule(mid, label, page, label, smi, "achiral", "not_applicable", "partial")


def mid(key):
    return molecules[key]["molecule_id"]


def charge_fact(value, evidence, raw=None):
    return f(raw or value, value, evidence)


def add_participant(r, key, role, evidence, *, eq=None, mol=None, stage="reaction"):
    pid = f"p{len(r['participants']) + 1}_{key}"
    p = st.new_participant(pid, mid(key), role, evidence, stage=stage)
    ref_key = next((x["molecule_id"] for x in r["participants"] if x["role"] == "reactant"), mid(key))
    if eq is not None:
        p["charges"]["equivalents"] = charge_fact(q(eq, "equiv", f"{eq} equiv"), evidence, f"{eq} equiv")
        p["equivalence_reference"] = f(ref_key, ref_key, evidence, method="inference",
                                         derivation="The source expresses equivalents relative to the first-listed substrate.")
    if mol is not None:
        p["charges"]["mol_percent"] = charge_fact(q(mol, "mol%", f"{mol} mol%"), evidence, f"{mol} mol%")
        p["equivalence_reference"] = f(ref_key, ref_key, evidence, method="inference",
                                         derivation="Catalyst/ligand loading denominator is the first-listed substrate.")
    r["participants"].append(p)
    return pid


def classification(r, codes, evidence, coverage="partial", limitation=None):
    if coverage == "not_applicable":
        r["classification"]["bond_changes"] = []
        r["classification"]["labels"] = f(None, [], evidence, method="inference",
                                               derivation="No observed transformation; intended bond formation is not classified.")
        r["classification"]["coverage"] = f(None, "not_applicable", evidence, method="inference",
                                                 derivation="The source reports no applicable observed transformation.")
        r["classification"]["limitations"] = [limitation or "No observed transformation was assigned."]
        return
    if coverage == "unresolved":
        r["classification"]["labels"] = f()
        r["classification"]["coverage"] = f(None, "unresolved", evidence, method="inference",
                                                 derivation="The observed species or source drawing is insufficient for safe bond-change assignment.")
        r["classification"]["limitations"] = [limitation or "Exact observed bond changes are unresolved."]
        return
    category_number = {x["code"]: x["number"] for x in json.loads((ROOT / "skills/chemical-reaction-extraction/references/bond_taxonomy.json").read_text())["categories"]}
    changes, labels = [], []
    atoms_for = {"C-C": ("C", "C", 1, "formed"), "C-H": ("C", "H", 1, "hydrogen_gain"),
                 "C-O": ("C", "O", 1, "formed"), "DB-C-C": ("C", "C", 2, "order_increased")}
    for i, code in enumerate(codes, 1):
        a, b, after, kind = atoms_for[code]
        change = {"change_id": f"b{i}", "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
                  "before_order": 0 if kind != "order_increased" else 1, "after_order": after, "change_type": kind}
        if kind == "hydrogen_gain":
            change["hydrogen_count"] = 1
        changes.append(change)
        labels.append({"category_number": category_number[code], "code": code, "supporting_change_ids": [f"b{i}"]})
    r["classification"]["bond_changes"] = [f(None, c, evidence, method="inference",
                                                    derivation="Class-level net change inferred from the displayed reactant/product transformation; exact atom mapping is not asserted.") for c in changes]
    r["classification"]["labels"] = f(None, labels, evidence, method="inference",
                                           derivation="Labels include all and only the class-level formed/promoted bonds listed in bond_changes.")
    r["classification"]["coverage"] = f(None, "partial", evidence, method="inference",
                                             derivation="Exact atom identities, complete hydrogen accounting, coproducts and atom mapping remain unresolved.")
    r["classification"]["limitations"] = [limitation or "Class-level inference only; no complete atom mapping or coproduct accounting."]


def reaction(rid, page, record_type="specific_experiment", categories=None, evidence=None,
             name="Pd-catalyzed asymmetric decarboxylative addition", transformation="decarboxylative allylation"):
    evidence = evidence or ev(page)
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, evidence) if categories is not None else f()
    r["location"] = f(f"PDF page {page}", {"pdf_page": page}, evidence)
    r["reaction_name"] = f(name, name, evidence)
    r["transformation_type"] = f(None, transformation, evidence, method="inference",
                                      derivation="Transformation class inferred conservatively from the source caption/prose and displayed arrow.")
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
                   "unresolved": ["No Supporting Information was supplied; conditions, workup, purification and characterization absent from pages 1-6 remain unreported.",
                                  "Machine record requires expert source review; schema validity is not chemical verification."]}
    r["completeness"] = "partial"
    return r


def configure_method(r, method, evidence, ref, ligand=None, solvent=None, temp=None, time=None,
                     include_tfa=False):
    if method in ("table", "A", "B", "C"):
        add_participant(r, "pd2dba3", "catalyst", evidence, mol=5 if method == "C" else 2.5)
    if method == "table":
        add_participant(r, ligand, "ligand", evidence, mol=5)
        add_participant(r, "Et3N", "additive", evidence, eq=3)
    elif method == "A":
        add_participant(r, "L3", "ligand", evidence, mol=5)
        add_participant(r, "Et3N", "additive", evidence, eq=3)
    elif method == "B":
        add_participant(r, "L2", "ligand", evidence, mol=5)
        if include_tfa:
            add_participant(r, "TFA", "additive", evidence, mol=5)
    elif method == "C":
        add_participant(r, "L2", "ligand", evidence, mol=10)
        add_participant(r, "Et3N", "additive", evidence, eq=0.1)
    solvent = solvent or {"table": "DCM", "A": "DCM", "B": "THF", "C": "toluene"}.get(method)
    if solvent:
        sid = add_participant(r, solvent, "solvent", evidence)
        r["solvent_conditions"]["composition"] = f(molecules[solvent]["names"]["original"]["normalized_value"], [mid(solvent)], evidence)
        conc = {"table": 0.1, "A": 0.1, "B": 0.1, "C": 0.2}.get(method)
        if conc is not None:
            r["solvent_conditions"]["concentration"] = f(f"{conc} M", q(conc, "M", f"{conc} M"), evidence)
    temp = temp if temp is not None else ({"table": -15, "A": -15, "B": -15, "C": -15}.get(method))
    time = time if time is not None else ({"table": q(18, "h", "18 h"), "A": q(18, "h", "18 h"),
                                           "B": qr(3, 24, "h", "3-24 h"), "C": q(18, "h", "18 h")}.get(method))
    s = st.new_step("reaction", evidence)
    s["participant_ids"] = [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"]
    s["operation"] = f(f"{method} conditions", f"{method} conditions", evidence)
    if temp == "r.t.":
        s["temperature"] = f("r.t.", "r.t.", evidence)
    elif temp is not None:
        s["temperature"] = f(f"{temp} °C", q(temp, "°C", f"{temp} °C"), evidence)
    if time is not None:
        s["time"] = f(time.get("raw_text", str(time)), time, evidence)
    r["steps"] = [s]
    return s


def outcome(r, oid, kind, target_keys, value, evidence, *, basis=None, scope=None, yield_type=None):
    o = st.new_outcome(oid, kind, [mid(k) for k in target_keys], evidence)
    o["value"] = f(value if isinstance(value, str) else value.get("raw_text", value), value, evidence)
    if basis is not None:
        o["measurement_basis"] = f(basis, basis, evidence)
    if scope is not None:
        o["scope"] = f(scope, scope, evidence)
    if yield_type is not None:
        o["yield_type"] = f(yield_type, yield_type, evidence)
    r["outcomes"].append(o)
    return o


def add_standard_products(r, reactants, products, evidence, reactant_eq=None):
    for i, key in enumerate(reactants):
        add_participant(r, key, "reactant", evidence, eq=(reactant_eq[i] if reactant_eq else None))
    for key in products:
        add_participant(r, key, "product", evidence)


reactions = []

# Scheme 1: cited families and current generic transformation.
for rid, reactants, product, name in [
    ("literature_previous_decarboxylative_addition", ["prior_acid", "prior_allene"], "prior_product", "previous work: decarboxylative addition of beta-keto acid to allene"),
    ("literature_asymmetric_hydrofunctionalization_precedent", ["precedent_allene", "precedent_nu"], "precedent_product", "precedent: asymmetric hydrofunctionalization of alkoxyallene"),
]:
    r = reaction(rid, 2, "literature_example", None, name=name, transformation="literature reaction family")
    add_standard_products(r, reactants, [product], ev(2))
    if "decarboxylative" in rid:
        classification(r, ["C-C", "C-H"], ev(2))
    else:
        classification(r, [], ev(2), coverage="unresolved",
                       limitation="The generic H-Nu family spans different heteroatoms; no single C-O or C-N class is asserted.")
    reactions.append(r)
r = reaction("generic_current_transformation", 2, "generic_scheme", None,
             name="current work: Pd-catalyzed asymmetric reaction of heteroatom-substituted allenes")
add_standard_products(r, ["current_acid", "current_allene"], ["current_product"], ev(2))
classification(r, ["C-C", "C-H"], ev(2))
reactions.append(r)
generic_id = r["reaction_id"]

# Table 1 optimization; entry 3 is the same 3a experiment repeated as the optimized Scheme 2 reference.
table = [
    (1, "L1", "allene2a", "DCM", "prod3a", 34, 53),
    (2, "L2", "allene2a", "DCM", "prod3a", 61, 91),
    (3, "L3", "allene2a", "DCM", "prod3a", 91, 97),
    (4, "L3", "allene2a", "THF", "prod3a", 55, 89),
    (5, "L3", "allene2a", "MeCN", "prod3a", 37, 74),
    (6, "L3", "allene2b", "DCM", "prod3b", 42, 87),
]
for entry, lig, allene, solv, prod, y, ee in table:
    evidence = ev(2) + (ev(3) if entry == 3 else [])
    r = reaction(f"table1_entry{entry}", 2, categories=["optimization"] + (["substrate_scope"] if entry == 3 else []), evidence=evidence)
    r["general_scheme_ids"] = [generic_id]
    add_standard_products(r, ["acid1", allene], [prod], evidence, [1, 1] if entry == 3 else [None, None])
    configure_method(r, "table", evidence, "acid1", ligand=lig, solvent=solv)
    outcome(r, "yield", "yield", [prod], q(y, "%", f"{y}%"), evidence, basis="isolated", yield_type="isolated yield")
    outcome(r, "ee", "selectivity", [prod], q(ee, "%", f"{ee}% ee"), evidence, basis="HPLC analysis", scope={"metric": "ee", "major_configuration": "R where explicitly assigned"})
    if entry == 3:
        r["notes"].append(f(None, "Table 1 entry 3 and the optimized 3a scope result are the same 91%/97% ee experiment; evidence is merged rather than double-counted.", evidence,
                              method="inference", derivation="Identical substrate pair, catalyst, ligand, solvent, temperature, time, yield and ee."))
    classification(r, ["C-C", "C-H"], evidence)
    reactions.append(r)

# Scheme 2, products 4-7; 3a is already Table 1 entry 3.
for n, acid, y, ee_value in [(4, "acid4", 86, 97), (5, "acid5", 77, 93), (6, "acid6", 91, 95), (7, "acid7", 72, 92)]:
    evidence = ev(3)
    r = reaction(f"scheme2_product{n}", 3, categories=["substrate_scope"], evidence=evidence)
    r["general_scheme_ids"] = [generic_id]
    add_standard_products(r, [acid, "allene2a"], [f"prod{n}"], evidence, [1, 1])
    configure_method(r, "A", evidence, acid)
    outcome(r, "yield", "yield", [f"prod{n}"], q(y, "%", f"{y}%"), evidence)
    outcome(r, "ee", "selectivity", [f"prod{n}"], q(ee_value, "%", f"{ee_value}% ee"), evidence,
            scope={"metric": "ee", "major_enantiomer": "source-depicted wedge; absolute label not printed"})
    classification(r, ["C-C", "C-H"], evidence)
    reactions.append(r)

# Scheme 3, cyclic amidoallene 8: additive and no-additive 9 are separate experiments.
scheme3_b = [("9_without_tfa", "acid1", 9, 95, 87, False), ("9_with_tfa", "acid1", 9, 86, 93, True),
             ("10", "acid10", 10, 94, 94, True), ("11", "acid11", 11, 94, 92, True),
             ("12", "acid12", 12, 95, 94, True), ("13", "acid13", 13, 90, 89, True),
             ("14", "acid14", 14, 86, 94, True), ("15", "acid15", 15, 81, 94, True),
             ("16", "acid16", 16, 86, 89, False)]
for suffix, acid, n, y, ee_value, tfa in scheme3_b:
    evidence = ev(3)
    r = reaction(f"scheme3_product{suffix}", 3, categories=["substrate_scope"])
    r["general_scheme_ids"] = [generic_id]
    add_standard_products(r, [acid, "allene8"], [f"prod{n}"], evidence, [1, 1])
    configure_method(r, "B", evidence, acid, include_tfa=tfa)
    outcome(r, "yield", "yield", [f"prod{n}"], q(y, "%", f"{y}%"), evidence)
    outcome(r, "ee", "selectivity", [f"prod{n}"], q(ee_value, "%", f"{ee_value}% ee"), evidence,
            scope={"metric": "ee", "major_configuration": "R for 9/14; otherwise source-depicted wedge without printed R/S"})
    if not tfa:
        r["notes"].append(f("CF3CO2H was not used", "CF3CO2H was not used", evidence))
    classification(r, ["C-C", "C-H"], evidence)
    reactions.append(r)

# Scheme 3, acyclic amidoallene 17.
for n, y, ratio, ee_value in [(18, 73, (5, 1), 93), (19, 64, (20, 1), 95), (20, 65, (21, 1), 95),
                              (21, 67, (12, 1), 94), (22, 66, (4, 1), 95)]:
    evidence = ev(3)
    acid = f"acid{n}"
    r = reaction(f"scheme3_product{n}", 3, categories=["substrate_scope"])
    r["general_scheme_ids"] = [generic_id]
    add_standard_products(r, [acid, "allene17"], [f"prod{n}"], evidence, [2, 1])
    configure_method(r, "C", evidence, acid)
    outcome(r, "yield", "yield", [f"prod{n}"], q(y, "%", f"{y}%"), evidence)
    outcome(r, "branch_linear", "selectivity", [f"prod{n}"], {"ratio": list(ratio), "order": ["branched", "linear"], "raw_text": f"b:l = {ratio[0]}:1"}, evidence,
            scope={"metric": "branched:linear ratio", "order": ["branched", "linear"]})
    outcome(r, "ee", "selectivity", [f"prod{n}"], q(ee_value, "%", f"{ee_value}% ee"), evidence,
            scope={"metric": "ee", "target": "branched product; absolute configuration not printed"})
    classification(r, ["C-C", "C-H"], evidence)
    reactions.append(r)

# Scheme 4 / Scheme 5B merged experiment for 23 + 2a.
evidence = ev(3) + ev(4)
r = reaction("scheme4_23_plus_2a", 3, categories=["substrate_scope", "mechanistic_study"], evidence=evidence)
r["general_scheme_ids"] = [generic_id]
add_standard_products(r, ["malonate23", "allene2a"], ["desired24", "prod38", "prod39"], evidence, [1, 1])
configure_method(r, "A", evidence, "malonate23", temp="r.t.")
# Method A implies CH2Cl2 on page 3 while Scheme 5B explicitly prints THF.
alt_dcm = f("Method A: CH2Cl2 (0.1 M)", [mid("DCM")], ev(3), method="inherited_general_condition")
alt_thf = f("THF", [mid("THF")], ev(4))
r["solvent_conditions"]["composition"] = f(["Method A", "THF"], None, evidence, status="conflict", alternatives=[alt_dcm, alt_thf])
add_participant(r, "THF", "solvent", ev(4))
outcome(r, "desired_no_reaction", "no_reaction", ["desired24"], "no rxn", ev(3), scope={"desired_product": "24"})
outcome(r, "adduct38_yield", "yield", ["prod38"], q(33, "%", "33%"), ev(4))
outcome(r, "adduct39_yield", "yield", ["prod39"], q(22, "%", "22%"), ev(4))
r["review"]["unresolved"].append("Solvent conflict retained: Scheme 4 calls Method A (defined with CH2Cl2), while Scheme 5B explicitly prints THF for the same 23 + 2a experiment.")
classification(r, ["C-O", "C-H"], evidence, limitation="Only the observed hydrocarboxylation adduct classes are represented; desired decarboxylative product 24 was not formed.")
reactions.append(r)

# Scheme 4, 25 + 2a and 25 + 8.
for rid, reactants, product, y, ee_value, method, temp, time, extra_selectivity in [
    ("scheme4_product26", ["malonate25", "allene2a"], "prod26", 72, 91, "A", "r.t.", q(18, "h", "18 h"), None),
    ("scheme4_product27", ["malonate25", "allene8"], "prod27", 92, 87, None, 40, q(18, "h", "18 h"), (14, 1)),
]:
    evidence = ev(3)
    r = reaction(rid, 3, categories=["substrate_scope"])
    r["general_scheme_ids"] = [generic_id]
    add_standard_products(r, reactants, [product], evidence, [1, 1] if method == "A" else None)
    if method == "A":
        configure_method(r, "A", evidence, reactants[0], temp=temp, time=time)
    else:
        add_participant(r, "pd2dba3", "catalyst", evidence, mol=2.5)
        add_participant(r, "L4", "ligand", evidence, mol=5)
        add_participant(r, "THF", "solvent", evidence)
        r["solvent_conditions"]["composition"] = f("THF", [mid("THF")], evidence)
        r["solvent_conditions"]["concentration"] = f("0.1 M", q(0.1, "M", "0.1 M"), evidence)
        s = st.new_step("reaction", evidence); s["participant_ids"] = [p["participant_id"] for p in r["participants"]]
        s["operation"] = f("reaction", "reaction", evidence); s["temperature"] = f("40 °C", q(40, "°C", "40 °C"), evidence); s["time"] = f("18 h", q(18, "h", "18 h"), evidence); r["steps"] = [s]
    outcome(r, "yield", "yield", [product], q(y, "%", f"{y}%"), evidence)
    outcome(r, "ee", "selectivity", [product], q(ee_value, "%", f"{ee_value}% ee"), evidence, scope={"metric": "ee"})
    if extra_selectivity:
        outcome(r, "branch_linear", "selectivity", [product], {"ratio": list(extra_selectivity), "order": ["branched", "linear"], "raw_text": "approximately 14:1"}, evidence,
                scope={"metric": "branch selectivity", "order": ["branched", "linear"], "qualifier": "approximately"})
    classification(r, ["C-C", "C-H"], evidence)
    reactions.append(r)

# Scheme 4 two-step RCM applications.
for start, prod, inter, stage1, stage2, overall in [("allene28", "prod29", "intermediate29", 63, 92, 58),
                                                    ("allene30", "prod31", "intermediate31", 69, 90, 62)]:
    evidence = ev(3)
    r = reaction(f"scheme4_{prod}", 3, categories=["derivatization"])
    r["general_scheme_ids"] = [generic_id]
    add_standard_products(r, [start, "acid1"], [inter, prod], evidence, [1, 1])
    # Method A agents and stage 1.
    configure_method(r, "A", evidence, start)
    r["steps"][0]["step_id"] = "hydrofunctionalization"
    # Stage 2 catalyst/solvent and room-temperature hold.
    add_participant(r, "RuCat", "catalyst", evidence, mol=10)
    add_participant(r, "DCM", "solvent", evidence)
    s2 = st.new_step("ring_closing_metathesis", evidence)
    s2["participant_ids"] = [p["participant_id"] for p in r["participants"] if p["molecule_id"] in (mid(inter), mid("RuCat"), mid("DCM"))]
    s2["operation"] = f("ring-closing metathesis", "ring-closing metathesis", evidence)
    s2["temperature"] = f("r.t.", "r.t.", evidence); s2["time"] = f("1 h", q(1, "h", "1 h"), evidence)
    r["steps"].append(s2)
    outcome(r, "stage1_yield", "yield", [inter], q(stage1, "%", f"{stage1}%"), evidence, scope={"stage": "Method A hydrofunctionalization"})
    outcome(r, "stage2_yield", "yield", [prod], q(stage2, "%", f"{stage2}%"), evidence, scope={"stage": "ring-closing metathesis", "basis": "stage intermediate"})
    outcome(r, "overall_yield", "yield", [prod], q(overall, "%", f"{overall}% over two steps"), evidence, scope={"span": f"{start} to {prod}", "type": "two-step overall"})
    outcome(r, "ee", "selectivity", [prod], q(93, "%", "93% ee"), evidence, scope={"metric": "ee"})
    classification(r, ["C-C", "C-H", "DB-C-C"], evidence)
    reactions.append(r)

# Scheme 5A competition control.
evidence = ev(4)
r = reaction("scheme5_competition_32_8", 4, categories=["control", "mechanistic_study"])
add_standard_products(r, ["ester32", "allene8"], ["species33", "species34"], evidence)
configure_method(r, "B", evidence, "ester32", include_tfa=True)
outcome(r, "no_33", "no_reaction", ["species33"], "32 showed no reactivity towards formation of 33", evidence)
outcome(r, "no_34", "other", ["species34"], "ion pair 34 could not be formed (author interpretation)", evidence,
        scope={"type": "author mechanistic interpretation"})
classification(r, [], evidence, coverage="not_applicable")
reactions.append(r)

# Scheme 5B NMR monitoring.
r = reaction("scheme5_nmr_monitoring_1_2a", 4, categories=["mechanistic_study"])
add_standard_products(r, ["acid1", "allene2a"], ["species35", "E36", "species37"], evidence)
configure_method(r, "A", evidence, "acid1", solvent="CD2Cl2")
outcome(r, "35_not_observed", "other", ["species35"], "not observed", evidence, scope={"observation": "1H NMR monitoring"})
outcome(r, "36_not_observed", "other", ["E36"], "not observed", evidence, scope={"observation": "1H NMR monitoring"})
outcome(r, "37_observed", "other", ["species37"], "observed", evidence, scope={"observation": "1H NMR monitoring"})
r["review"]["unresolved"].append("The identity/drawing of species 37 in the monitoring panel is not confidently reconcilable with the crossover panel; it remains unresolved.")
classification(r, [], evidence, coverage="unresolved")
reactions.append(r)

# Scheme 5C benzoic acid hydrocarboxylation.
r = reaction("scheme5_benzoic_acid_2a", 4, categories=["mechanistic_study"])
add_standard_products(r, ["benzoic", "allene2a"], ["prod40", "prod41"], evidence)
configure_method(r, "A", evidence, "benzoic", temp="r.t.")
outcome(r, "yield40", "yield", ["prod40"], q(38, "%", "38%"), evidence)
outcome(r, "yield41", "yield", ["prod41"], q(32, "%", "32%"), evidence)
classification(r, ["C-O", "C-H"], evidence)
reactions.append(r)

# Scheme 5D: five distinct mechanistic reactions.
mechanistic_d = [
    ("scheme5_E36", ["E36"], ["prod3a"], [("prod3a", 30, 63)], None),
    ("scheme5_E36_plus_acid1", ["E36", "acid1"], ["prod3a"], [("prod3a", 46, 75)], None),
    ("scheme5_Z36", ["Z36"], ["prod3a"], [("prod3a", 34, 79)], None),
    ("scheme5_41_plus_acid1", ["prod41", "acid1"], ["prod3a"], [("prod3a", 25, 76)], None),
    ("scheme5_crossover_E36_37", ["E36", "species37"], ["prod3a", "prod4"], [("prod3a", 20, 88), ("prod4", 28, 78)], "37 (1.0 eq)"),
]
for rid, reactants, products, results, note in mechanistic_d:
    r = reaction(rid, 4, categories=["mechanistic_study"])
    add_standard_products(r, reactants, products, evidence)
    if rid == "scheme5_E36_plus_acid1":
        acid_p = next(p for p in r["participants"] if p["molecule_id"] == mid("acid1"))
        acid_p["charges"]["equivalents"] = f("1 equiv", q(1, "equiv", "1 equiv"), evidence)
        acid_p["equivalence_reference"] = f(mid("E36"), mid("E36"), evidence, method="inference",
                                                derivation="Scheme 5 explicitly states addition of 1 equivalent of acid 1 relative to (E)-36.")
    if rid == "scheme5_crossover_E36_37":
        p37 = next(p for p in r["participants"] if p["molecule_id"] == mid("species37"))
        p37["charges"]["equivalents"] = f("1.0 eq", q(1.0, "equiv", "1.0 eq"), evidence)
        p37["equivalence_reference"] = f(mid("E36"), mid("E36"), evidence, method="inference",
                                             derivation="Scheme 5 explicitly states 37 (1.0 eq) relative to (E)-36.")
    configure_method(r, "A", evidence, reactants[0])
    for product, y, ee_value in results:
        outcome(r, f"yield_{product}", "yield", [product], q(y, "%", f"{y}%"), evidence)
        outcome(r, f"ee_{product}", "selectivity", [product], q(ee_value, "%", f"{ee_value}% ee"), evidence,
                scope={"metric": "ee", "product": product})
    if note:
        r["notes"].append(f(note, note, evidence))
    classification(r, ["C-C"], evidence, limitation="Partial decarboxylative allylation class only; exact atom mapping and leaving fragments are unresolved.")
    reactions.append(r)

# Footnote 19: a real but detail-limited benzoic-acid/amidoallene observation.
evidence = ev(5)
r = reaction("footnote19_benzoic_acid_amidoallene8", 5, categories=["mechanistic_study"], evidence=evidence,
             name="benzoic acid reaction with amidoallene 8", transformation="hydrocarboxylation observation")
add_standard_products(r, ["benzoic", "allene8"], ["amido_benzoic_products"], evidence)
missing_step = st.new_step("reported_reaction_conditions_unreported", evidence)
missing_step["participant_ids"] = [p["participant_id"] for p in r["participants"]]
missing_step["operation"] = f("reaction with benzoic acid showed a similar pattern", "reaction with benzoic acid showed a similar pattern", evidence)
r["steps"] = [missing_step]
outcome(r, "similar_pattern", "other", ["amido_benzoic_products"], "Amidoallene 8 showed a similar pattern in the reaction with benzoic acid.", evidence,
        scope={"comparison": "similar to branch-adduct formation/migration of 40 to 41"})
r["review"]["unresolved"].append("Footnote 19 supplies no conditions, structures, yields, ratio, or measurement method; no details are inferred from the alkoxyallene experiment.")
classification(r, [], evidence, coverage="unresolved")
reactions.append(r)

# Scheme 6 is a proposed catalytic cycle, not executed intermediates.
evidence = ev(4)
r = reaction("scheme6_proposed_catalytic_cycle", 4, "proposed_mechanism", None, evidence=evidence,
             name="proposed Pd-H / pi-allyl catalytic cycle", transformation="proposed mechanism")
r["notes"].append(f("Pd-H addition to allene; A/B equilibrium; inner-sphere C-C bond formation via C; decarboxylation to product", 
                      "Pd-H addition to allene; A/B equilibrium; inner-sphere C-C bond formation via C; decarboxylation to product", evidence))
r["review"]["unresolved"].append("Intermediates A-D are author-proposed and are not treated as isolated or directly observed experiments.")
classification(r, [], evidence, coverage="unresolved")
reactions.append(r)

bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": f"{DOC}:context:page1_cover", "document_id": DOC, "type": "metadata_and_cover_art",
     "evidence_ids": ev(1), "note": "Title, authors, citation and cover illustration; no synthesis experiment is created from the artwork."},
    {"context_id": f"{DOC}:context:xray", "document_id": DOC, "type": "structure_assignment",
     "evidence_ids": ev(3) + ev(5), "note": "Absolute configuration of 14 was verified by X-ray and described as consistent with 3b; CCDC 2084084 is external and no CIF/SI was supplied."},
    {"context_id": f"{DOC}:context:mechanism_discussion", "document_id": DOC, "type": "mechanistic_discussion",
     "evidence_ids": ev(4), "note": "Author interpretation of branch-adduct kinetics, pi-allyl equilibration, inner-sphere C-C formation and tight ion pairing; distinct from the concrete Scheme 5 experiments."},
    {"context_id": f"{DOC}:context:references_and_negative_spectral_inventory", "document_id": DOC, "type": "reference_list_and_negative_inventory",
     "evidence_ids": ev(5) + ev(6), "note": "Pages 5-6 are references/footnotes/manuscript dates. No standalone mass-spectrum peak table, molecular-weight table or fragment-relative-abundance block occurs on pages 1-6."},
]

# Raw inventories and exact page-bounded coverage.
page_records = {p: [] for p in range(1, 7)}
for rxn in reactions:
    for eid in rxn["evidence_ids"]:
        if ":page" in eid:
            try:
                p = int(eid.split(":page", 1)[1].split(":", 1)[0])
            except ValueError:
                continue
            if 1 <= p <= 6 and rxn["reaction_id"] not in page_records[p]:
                page_records[p].append(rxn["reaction_id"])

coverage = {"document_id": DOC, "semantic_accuracy_verified": False, "pages": [
    {"page": 1, "state": "attempted", "regions": [
        {"region": "full-page cover/title area", "finding": "Title, authors, citation and cover illustration; no experimental reaction arrow or table row.", "evidence_ids": ev(1)}],
     "record_ids": page_records[1], "unresolved": ["Cover artwork is illustrative and is not interpreted as an experimental structure record."]},
    {"page": 2, "state": "attempted", "regions": [
        {"region": "Scheme 1, upper-right", "finding": "Two cited prior-work reaction families and one current generic transformation; conceptual metal-catalyzed AAA panels are not counted as concrete current experiments.", "evidence_ids": ev(2)},
        {"region": "Table 1, lower-right", "finding": "Six independent optimization rows with isolated yields and HPLC ee values; entry 3 is the repeated optimized 3a scope result.", "evidence_ids": ev(2)},
        {"region": "optimization prose", "finding": "Common -15 °C/18 h/Pd2(dba)3/ligand/Et3N conditions and explicit no-linear-isomer statement.", "evidence_ids": ev(2)}],
     "record_ids": page_records[2], "unresolved": ["Exact reactant mmol, workup and purification are not reported on the supplied page."]},
    {"page": 3, "state": "attempted", "regions": [
        {"region": "Scheme 2", "finding": "Four new scope arrows 4-7; 3a is the identical optimized Table 1 entry 3 and is not duplicated.", "evidence_ids": ev(3)},
        {"region": "Scheme 3", "finding": "Nine cyclic-amidoallene experiments (including separate 9 with/without TFA) and five acyclic-amidoallene experiments 18-22.", "evidence_ids": ev(3)},
        {"region": "Scheme 4", "finding": "Five applications: failed desired 24 (cross-linked to observed 38/39), products 26/27, and two two-step RCM sequences to 29/31.", "evidence_ids": ev(3)},
        {"region": "mechanistic prose at lower-right", "finding": "Introduces Scheme 5 controls and NMR monitoring; concrete arrows are inventoried on page 4.", "evidence_ids": ev(3)}],
     "record_ids": page_records[3], "unresolved": ["Scheme yields outside Table 1 do not state measurement basis; individual Method B times are given only as 3-24 h."]},
    {"page": 4, "state": "attempted", "regions": [
        {"region": "Scheme 5A-B", "finding": "Methyl-ester competition control, 1H NMR monitoring, and the 23+2a adduct experiment; the latter is merged with Scheme 4 rather than duplicated.", "evidence_ids": ev(4)},
        {"region": "Scheme 5C-D", "finding": "Benzoic-acid hydrocarboxylation and five separate E/Z-36, 41, extra-acid, and crossover experiments with product-specific yields/ee.", "evidence_ids": ev(4)},
        {"region": "Scheme 6", "finding": "Author-proposed catalytic cycle retained as proposed mechanism, not executed intermediate experiments.", "evidence_ids": ev(4)},
        {"region": "mechanistic discussion", "finding": "Author interpretation of kinetic branch adduct, pi-allyl equilibration and inner-sphere C-C formation.", "evidence_ids": ev(4)}],
     "record_ids": page_records[4], "unresolved": ["Compound 37 drawing/identity cannot be reconciled confidently across Scheme 5 panels; no SI was supplied."]},
    {"page": 5, "state": "attempted", "regions": [
        {"region": "reference list and footnotes 1-22", "finding": "Bibliography plus footnotes; footnote 19 is a detail-limited real benzoic-acid/amidoallene observation and is retained separately.", "evidence_ids": ev(5)},
        {"region": "CCDC 2084084 notice", "finding": "External crystallographic availability notice, not a supplied synthesis experiment or associated SI file.", "evidence_ids": ev(5)}],
     "record_ids": page_records[5], "unresolved": ["External CCDC/CIF and advertised Supporting Information were not supplied and were not associated."]},
    {"page": 6, "state": "attempted", "regions": [
        {"region": "full page", "finding": "Continuation of references 22-23 and manuscript dates only; no reaction arrow, table row, procedure or standalone spectral table.", "evidence_ids": ev(6)}],
     "record_ids": page_records[6], "unresolved": []},
]}

page_extraction = {
    "document_id": DOC, "pages_attempted": [1, 2, 3, 4, 5, 6],
    "source": "Direct visual inspection of original rendered page images; extracted text used only as a transcription aid.",
    "record_count": len(reactions), "molecule_count": len(molecules),
    "deduplication": [
        {"source_regions": ["Table 1 entry 3", "Scheme 2 optimized 3a"], "decision": "one reaction record", "reaction_id": f"{DOC}:table1_entry3"},
        {"source_regions": ["Scheme 4 23+2a no desired 24", "Scheme 5B 23+2a gives 38/39"], "decision": "one reaction record with all outcomes and solvent conflict", "reaction_id": f"{DOC}:scheme4_23_plus_2a"},
    ],
    "negative_inventory": "No standalone mass-spectrum peak table, molecular-weight table or fragment-relative-abundance block on pages 1-6.",
    "page_record_ids": page_records,
}

dump(RUN / "submitted.json", bundle)
dump(RAW / "stereo_checks.json", stereo_checks)
dump(RAW / "page_extraction.json", page_extraction)
dump(RAW / "source_mentions.json", st.source_mentions(bundle, ["2a", "3a", "36", "37", "CF3CO2H"], context=180))
dump(RUN / "coverage/pages_0001_0006.json", coverage)
print(json.dumps({"reactions": len(reactions), "molecules": len(molecules), "pages": [1, 2, 3, 4, 5, 6]}, ensure_ascii=False))
