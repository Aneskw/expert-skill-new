#!/usr/bin/env python3
"""Build the source-derived v3.1 submission for doc_19e82cae555b0ca9.

This file is an extraction artifact, not a change to the shared extraction tools.
All records remain machine records and are imported as needs_review.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DATASET = ROOT / "runs/current/dataset.json"
OUT = ROOT / "runs/current/submitted.json"
RAW = ROOT / "runs/current/raw_model"
DOC = "doc_19e82cae555b0ca9"
E = {p: f"{DOC}:page{p}:image" for p in range(1, 6)}
T = {p: f"{DOC}:page{p}:text" for p in range(1, 6)}


def f(raw=None, norm=None, ev=(), method="original_extraction", derivation=None,
      status=None, verification=None, alternatives=()):
    return st.fact(raw, norm, ev, method=method, derivation=derivation,
                   status=status, verification=verification, alternatives=alternatives)


def inferred(norm, ev, derivation, raw=None):
    return f(raw, norm, ev, method="inference", derivation=derivation)


def derived(norm, ev, derivation, raw=None):
    return f(raw, norm, ev, method="structure_derivation", derivation=derivation)


def na():
    return f(status="not_applicable", method="not_applicable")


def qty(number, unit, raw=None, qualifier=None):
    x = {"number": number, "unit": unit}
    if raw is not None:
        x["raw_text"] = raw
    if qualifier:
        x["qualifier"] = qualifier
    return x


dataset = json.loads(DATASET.read_text())
bundle = st.empty(dataset["documents"])
molecules = {}
stereo_checks = {}


def add_molecule(mid, label, ev, smiles=None, *, stereo_kind="auto", assignment=None,
                 identity_note=None, aliases=None):
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"] = [{
        "document_id": DOC,
        "compound_label": f(label, label, ev),
        "evidence_ids": list(ev),
    }]
    m["names"]["original"] = f(label, label, ev)
    m["names"]["english"] = f(label, label, ev)
    m["names"]["abbreviations"] = f(aliases, aliases, ev) if aliases else f()
    m["names"]["aliases"] = f()
    m["cas_numbers"] = f()
    m["external_properties"] = []
    m["identity_note"] = (f(identity_note, identity_note, ev) if identity_note else
                          f("Structure identity is a machine candidate from the displayed drawing; RDKit parsing is not source verification.",
                            "Structure identity is a machine candidate from the displayed drawing; RDKit parsing is not source verification.", ev))
    if smiles is None:
        m["representation_type"] = "unresolved"
        m["smiles"] = f()
        m["generic_template"] = None
        for k in ("atom_centers", "double_bonds", "axial", "mixture", "assignment",
                  "isomeric_smiles", "identity_key"):
            m["stereochemistry"][k] = f()
        for k in ("chirality_status", "geometry_status", "representation_status"):
            m["stereochemistry"][k] = inferred("unresolved", ev, "The source label or dense structure was retained without a complete machine-readable structure.")
    else:
        diag = st.stereo_inspect(smiles)
        stereo_checks[mid] = diag
        canonical = diag["canonical_isomeric_smiles"]
        m["representation_type"] = "smiles"
        m["smiles"] = derived(canonical, ev, "Connectivity candidate transcribed from the displayed structure and canonicalized with RDKit; parsing does not certify the drawing assignment.", raw={"source": "displayed structure", "label": label})
        m["generic_template"] = None
        m["stereochemistry"]["atom_centers"] = derived(diag["atom_centers"] or [], ev, "RDKit potential/assigned atom-center diagnostics on the candidate connectivity.", raw="displayed structure")
        m["stereochemistry"]["double_bonds"] = derived(diag["double_bonds"] or [], ev, "RDKit double-bond diagnostics on the candidate connectivity.", raw="displayed structure")
        m["stereochemistry"]["axial"] = f()
        m["stereochemistry"]["mixture"] = f()
        if stereo_kind == "relative":
            m["stereochemistry"]["assignment"] = f(assignment, assignment, ev)
            m["stereochemistry"]["chirality_status"] = f(assignment, "partially_specified", ev)
            m["stereochemistry"]["geometry_status"] = inferred("not_applicable", ev, "The depicted alkene is terminal; E/Z is not applicable. The ring relationship is stored as relative stereochemistry, not alkene geometry.")
            m["stereochemistry"]["representation_status"] = inferred("not_encoded", ev, "The source specifies only a relative syn/anti relationship and no absolute enantiomer; no arbitrary @/@@ encoding was introduced.")
            m["stereochemistry"]["isomeric_smiles"] = f()
            m["stereochemistry"]["identity_key"] = f()
        elif stereo_kind == "unspecified":
            m["stereochemistry"]["assignment"] = f("No absolute configuration is assigned in the displayed structure.", "No absolute configuration assigned", ev)
            m["stereochemistry"]["chirality_status"] = f("No absolute configuration is assigned in the displayed structure.", "unspecified", ev)
            m["stereochemistry"]["geometry_status"] = inferred("not_applicable", ev, "The allene is terminal and not axially stereogenic; no E/Z alkene is present.")
            m["stereochemistry"]["representation_status"] = inferred("not_encoded", ev, "Potential tetrahedral stereochemistry is unassigned in the source and remains unencoded.")
            m["stereochemistry"]["isomeric_smiles"] = f()
            m["stereochemistry"]["identity_key"] = derived(diag["stereo_graph_key"], ev, "RDKit stereo graph key for the explicitly unspecified candidate graph; not a sample identity.", raw="candidate connectivity")
        else:
            potential = [x for x in diag["potential_stereo"] if x["specified"] == "Unspecified"]
            if potential:
                cstat, rstat = "unspecified", "not_encoded"
                iso = None
                key = diag["stereo_graph_key"]
            else:
                cstat, rstat = "achiral", "encoded"
                iso = canonical
                key = diag["stereo_graph_key"]
            m["stereochemistry"]["assignment"] = f()
            m["stereochemistry"]["chirality_status"] = derived(cstat, ev, "Status derived from the candidate graph; this does not establish source identity.", raw="candidate connectivity")
            m["stereochemistry"]["geometry_status"] = derived("not_applicable", ev, "No source-assigned E/Z geometry is applicable to this candidate.", raw="candidate connectivity")
            m["stereochemistry"]["representation_status"] = derived(rstat, ev, "Representation status derived from RDKit diagnostics on the candidate graph.", raw="candidate connectivity")
            m["stereochemistry"]["isomeric_smiles"] = (derived(iso, ev, "Canonical isomeric SMILES from RDKit candidate graph.", raw="candidate connectivity") if iso else f())
            m["stereochemistry"]["identity_key"] = derived(key, ev, "RDKit stereo graph key; not a sample identity.", raw="candidate connectivity")
    molecules[mid] = m
    return mid


def add_template(mid, label, ev, fixed, variable, stereo_status):
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"] = [{"document_id": DOC, "compound_label": f(label, label, ev), "evidence_ids": list(ev)}]
    m["names"]["original"] = f(label, label, ev)
    m["names"]["english"] = f(label, label, ev)
    m["cas_numbers"] = f()
    m["representation_type"] = "generic_template"
    m["smiles"] = f()
    m["generic_template"] = {
        "fixed_structure": f(fixed, fixed, ev),
        "variables": [{
            "label": "R",
            "attachment_points": f("R attached at the carbon alpha to urea N", {"site": "alpha carbon"}, ev),
            "definition": f(variable, variable, ev),
            "allowed_values": f("alkyl or aryl; concrete scope shown as 1b-o", "alkyl or aryl", ev),
            "constraints": f(),
        }],
        "instance_mapping": [],
        "reason": f("The source drawing uses R and the range 1b-o/2b-o/3b-o.", "Generic R-group drawing", ev),
    }
    for k in ("atom_centers", "double_bonds", "axial", "mixture", "isomeric_smiles", "identity_key"):
        m["stereochemistry"][k] = f()
    m["stereochemistry"]["assignment"] = f(stereo_status, stereo_status, ev)
    m["stereochemistry"]["chirality_status"] = f(stereo_status, "partially_specified" if stereo_status in ("anti", "syn") else "unspecified", ev)
    m["stereochemistry"]["geometry_status"] = inferred("not_applicable", ev, "The depicted alkene/allene does not carry an E/Z assignment in this generic structure.")
    m["stereochemistry"]["representation_status"] = inferred("not_encoded", ev, "Generic R-group and relative/unassigned stereochemistry are not encoded as a concrete SMILES.")
    m["identity_note"] = f("Generic structure only; no concrete SMILES or CAS was invented.", "Generic structure only; no concrete SMILES or CAS was invented.", ev)
    m["external_properties"] = []
    molecules[mid] = m
    return mid


# Generic structures and the complete b-o scope registry.
add_template("m_generic_1", "tosylurea-allene 1 (R = alkyl or aryl)", [E[1], E[2], E[3]],
             "O=C(N-CH(R)-CH2-CH=C=CH2)-NHTs", "R = alkyl or aryl", "unassigned")
add_template("m_generic_2", "anti-tetrahydropyrimidinone 2", [E[1], E[2]],
             "4-R-6-vinyl-3-tosyl-tetrahydropyrimidin-2-one", "R = alkyl or aryl", "anti")
add_template("m_generic_3", "syn-tetrahydropyrimidinone 3", [E[1], E[2], E[3]],
             "4-R-6-vinyl-3-tosyl-tetrahydropyrimidin-2-one", "R = alkyl or aryl", "syn")

R = {
    "a": "c2ccccc2", "b": "c2ccccc2Br", "c": "c2cccc(Br)c2", "d": "c2ccc(Br)cc2",
    "e": "c2ccc(Cl)cc2", "f": "c2ccc(F)cc2", "g": "c2ccccc2C", "h": "c2ccc(C)cc2",
    "i": "c2ccc(OC)cc2", "j": "c2ccc(C#N)cc2", "k": "c2ccc3ccccc3c2",
    "l": "C2CCCCC2", "m": "CCCCCCCC", "n": "c2ccc([N+](=O)[O-])cc2", "o": "c2ccco2",
}
for x, rgroup in R.items():
    ev = [E[2]] if x != "a" else [E[1], E[2], E[3]]
    s1 = f"O=C(NC({rgroup})CC=C=C)NS(=O)(=O)c1ccc(C)cc1"
    sp = f"O=C1NC({rgroup})CC(C=C)N1S(=O)(=O)c1ccc(C)cc1"
    add_molecule(f"m_1{x}", f"tosylurea-allene 1{x}", ev, s1, stereo_kind="unspecified",
                 identity_note="Candidate connectivity from the displayed scope; substituent position follows the figure. Absolute configuration is not assigned.")
    add_molecule(f"m_2{x}", f"anti-tetrahydropyrimidinone 2{x}", [E[1], E[2], E[3]], sp,
                 stereo_kind="relative", assignment="anti relative configuration shown; absolute configuration/enantiomer not assigned")
    add_molecule(f"m_3{x}", f"syn-tetrahydropyrimidinone 3{x}", [E[1], E[2], E[3]], sp,
                 stereo_kind="relative", assignment="syn relative configuration shown; absolute configuration/enantiomer not assigned")

# Follow-up products.
add_molecule("m_4", "detosylated cyclic urea 4", [E[2], E[3]], "O=C1NC(c2ccccc2)CC(C=C)N1",
             stereo_kind="relative", assignment="relative configuration depicted; absolute configuration not assigned")
add_molecule("m_5", "1,3-diamine 5 dihydrochloride", [E[2], E[3]], "[NH3+]C(c1ccccc1)CC([NH3+])C=C.[Cl-].[Cl-]",
             stereo_kind="relative", assignment="relative configuration depicted; absolute configuration not assigned")
add_molecule("m_6", "hydrogenated tetrahydropyrimidinone 6", [E[3]], "O=C1NC(c2ccccc2)CC(CC)N1S(=O)(=O)c1ccc(C)cc1",
             stereo_kind="relative", assignment="relative configuration depicted; absolute configuration not assigned")
add_molecule("m_7", "aldehyde 7", [E[3]], "O=C1NC(c2ccccc2)CC(C=O)N1S(=O)(=O)c1ccc(C)cc1",
             stereo_kind="relative", assignment="relative configuration depicted; absolute configuration not assigned")
add_molecule("m_8", "linear hydroformylation product 8", [E[3]], "O=C1NC(c2ccccc2)CC(CCC=O)N1S(=O)(=O)c1ccc(C)cc1",
             stereo_kind="relative", assignment="relative configuration depicted; absolute configuration not assigned")
add_molecule("m_9", "N-methylated tetrahydropyrimidinone 9", [E[3]], "O=C1N(C)C(c2ccccc2)CC(C=C)N1S(=O)(=O)c1ccc(C)cc1",
             stereo_kind="relative", assignment="relative configuration depicted; absolute configuration not assigned")

# Reagents, catalysts, solvents, gases and explicitly reported workup agents.
small = {
    "m_dce": ("1,2-dichloroethane (DCE)", "ClCCCl", [E[2], E[3]]),
    "m_benzoic": ("benzoic acid", "O=C(O)c1ccccc1", [E[1], E[2], E[3]]),
    "m_pentafluoro": ("pentafluorobenzoic acid", "O=C(O)c1c(F)c(F)c(F)c(F)c1F", [E[2], E[3]]),
    "m_phenylisobutyric": ("2-methyl-2-phenylpropanoic acid (PhMe2CCO2H)", "CC(C)(C(=O)O)c1ccccc1", [E[1], E[2]]),
    "m_licl": ("lithium chloride", "[Li+].[Cl-]", [E[2]]),
    "m_mg": ("magnesium powder", "[Mg]", [E[3]]),
    "m_meoh": ("methanol", "CO", [E[3]]),
    "m_hcl": ("hydrochloric acid", "Cl", [E[3]]),
    "m_k2co3": ("potassium carbonate", "[K+].[K+].[O-]C([O-])=O", [E[3]]),
    "m_me2so4": ("dimethyl sulfate", "COS(=O)(=O)OC", [E[3]]),
    "m_mecn": ("acetonitrile", "CC#N", [E[3]]),
    "m_h2": ("hydrogen", "[H][H]", [E[3]]),
    "m_o3": ("ozone", "[O-][O+]=O", [E[3]]),
    "m_dcm": ("dichloromethane", "ClCCl", [E[3]]),
    "m_tpp": ("triphenylphosphine (TPP)", "P(c1ccccc1)(c1ccccc1)c1ccccc1", [E[3]]),
    "m_co": ("carbon monoxide", "[C-]#[O+]", [E[3]]),
    "m_toluene": ("toluene", "Cc1ccccc1", [E[3]]),
}
for mid, (name, smi, ev) in small.items():
    add_molecule(mid, name, ev, smi)
for mid, name, ev in [
    ("m_pd_dba2", "Pd(dba)2", [E[1], E[2], E[3]]),
    ("m_rh_cod_cl2", "[{Rh(cod)Cl}2]", [E[1], E[2], E[3]]),
    ("m_dpephos", "DPEphos", [E[1], E[2], E[3]]),
    ("m_ppts", "PPTS", [E[2]]),
    ("m_pdc", "Pd (5% on charcoal)", [E[3]]),
    ("m_rh_co2acac", "[Rh(CO)2acac]", [E[3]]),
    ("m_6dppon", "6-DPPon", [E[3]]),
]:
    add_molecule(mid, name, ev, None, identity_note="Source label retained; no complete organometallic/mixture structure was invented.")

# Source-labelled participants for the three prior-work transformations.
for mid, label in [
    ("m_menche_sub", "Menche (2010) allylcarbonate substrate"), ("m_menche_prod", "Menche (2010) tetrahydropyrimidinone product"),
    ("m_white_sub", "White (2013) olefin substrate"), ("m_white_prod", "White (2013) anti cyclic urea product"),
    ("m_hara_sub", "Hara (2015) olefin substrate"), ("m_hara_prod", "Hara (2015) anti cyclic urea product"),
    ("m_pdpph34", "[Pd(PPh3)4]"), ("m_pdoac2", "Pd(OAc)2"), ("m_pdtfa2", "Pd(TFA)2/bis-sulfoxide system"),
]:
    add_molecule(mid, label, [E[2]], None, identity_note="Prior-work source label/drawing retained for contextual literature classification; not reconstructed as a verified concrete identity.")


def participant(pid, mid, role, ev, amount=None, amount_kind=None, reference=None,
                stage="reaction", product_role=None):
    p = st.new_participant(pid, mid, role, ev, stage)
    if product_role:
        p["product_role"] = f(product_role, product_role, ev)
    if amount is not None:
        p["charges"][amount_kind] = f(amount.get("raw_text", amount), amount, ev)
        p["equivalence_reference"] = inferred(reference, ev, "Loading/equivalents denominator is the displayed substrate.", raw="relative to displayed substrate")
    return p


def change(cid, a, b, before, after, typ, hcount=None):
    x = {"change_id": cid, "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
         "before_order": before, "after_order": after, "change_type": typ}
    if hcount is not None:
        x["hydrogen_count"] = hcount
    return x


def classification(ev, changes=None, coverage="partial", limitation=None):
    if coverage == "not_applicable":
        return {"standard_id": st.STANDARD, "mapped_reaction_smiles": f(), "bond_changes": [],
                "labels": inferred([], ev, "No observed bond-forming transformation is classified for this no-reaction/stereoisomerization/proposed-mechanism record."),
                "coverage": inferred("not_applicable", ev, "Bond-formation classification is not applicable."),
                "limitations": [limitation or "No bond-forming transformation asserted."], "model_votes": []}
    if coverage == "unresolved":
        return {"standard_id": st.STANDARD, "mapped_reaction_smiles": f(), "bond_changes": [],
                "labels": f(), "coverage": inferred("unresolved", ev, "The contextual source does not support a complete bond-difference assignment."),
                "limitations": [limitation or "Bond changes unresolved."], "model_votes": []}
    changes = changes or []
    labels = {}
    for c in changes:
        if c["change_type"] in ("formed", "hydrogen_gain", "order_increased"):
            lab = st.bond_label(c["atoms"][0]["element"], c["atoms"][1]["element"], c["after_order"])
            labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(c["change_id"])
    deriv = "Class-level bond changes inferred from the displayed reactant/product drawings; atom mapping and full balance are unresolved."
    return {"standard_id": st.STANDARD, "mapped_reaction_smiles": f(),
            "bond_changes": [derived(c, ev, deriv) for c in changes],
            "labels": derived(list(labels.values()), ev, deriv),
            "coverage": inferred("partial", ev, deriv),
            "limitations": [limitation or "No atom-mapped reaction; implicit-H source, coproduct balance and exact atom correspondence remain unresolved."],
            "model_votes": []}


HYDROAMINATION = [
    change("b1", "C", "N", 0, 1, "formed"),
    change("b2", "C", "H", 0, 1, "hydrogen_gain", 1),
    change("b3", "C", "C", 2, 1, "order_decreased"),
]


def base_reaction(rid, ev, record_type="specific_experiment", categories=None,
                  location=None, name=None, transform=None, general=()):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories or [], categories or [], ev)
    r["location"] = f(location, location, ev)
    r["general_scheme_ids"] = list(general)
    r["reaction_name"] = f(name, name, ev) if name else f()
    r["transformation_type"] = inferred(transform, ev, "Transformation class inferred from the displayed reactant/product relationship.") if transform else f()
    r["reaction_scale"] = f()
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": ["No explicitly associated SI was supplied; exact charges, workup and purification remain unreported unless printed in the article."]}
    r["completeness"] = "partial"
    return r


def set_solvent(r, mid, raw, concentration, ev, ratio=None):
    pid = f"p_{len(r['participants'])+1}_{mid}"
    r["participants"].append(participant(pid, mid, "solvent", ev))
    r["solvent_conditions"]["composition"] = f(raw, [mid], ev)
    r["solvent_conditions"]["concentration"] = f(concentration.get("raw_text", concentration), concentration, ev)
    r["solvent_conditions"]["concentration_basis"] = f("substrate concentration", "substrate", ev)
    if ratio:
        r["solvent_conditions"]["mixture_ratio"] = f(ratio, ratio, ev)


def add_step(r, sid, ev, pids, temp=None, time=None, operation="reaction"):
    s = st.new_step(sid, ev)
    s["participant_ids"] = list(pids)
    s["operation"] = f(operation, operation, ev)
    if temp is not None:
        s["temperature"] = f(temp.get("raw_text", temp) if isinstance(temp, dict) else temp, temp, ev)
    if time is not None:
        s["time"] = f(time.get("raw_text", time) if isinstance(time, dict) else time, time, ev)
    r["steps"].append(s)


def add_outcome(r, oid, kind, targets, ev, raw, norm, *, basis=None, yield_type=None, scope=None, role="main product"):
    o = st.new_outcome(oid, kind, targets, ev)
    o["product_role"] = f(role, role, ev)
    o["value"] = f(raw, norm, ev)
    o["yield_type"] = f(yield_type, yield_type, ev) if yield_type else f()
    o["measurement_basis"] = f(basis, basis, ev) if basis else f()
    o["scope"] = f(scope, scope, ev) if scope else f()
    r["outcomes"].append(o)


reactions = []
PDGEN = f"{DOC}:generic_pd_anti"
RHGEN = f"{DOC}:generic_rh_syn"

# Two general methods, with the page-1 graphical abstract retained as duplicate overview evidence.
for rid, catalyst, product, selectivity, ev, loc in [
    ("generic_pd_anti", "m_pd_dba2", "m_generic_2", "anti", [E[1], E[2]], "page 1 graphical abstract and page 2 'This work'/Scheme 1"),
    ("generic_rh_syn", "m_rh_cod_cl2", "m_generic_3", "syn", [E[1], E[2], E[3]], "page 1 graphical abstract and page 3 Scheme 2"),
]:
    r = base_reaction(rid, ev, "generic_scheme", [], loc, "intramolecular addition of tosylureas to allenes", "intramolecular hydroamination/cyclization")
    r["participants"] = [participant("substrate", "m_generic_1", "reactant", ev), participant("catalyst", catalyst, "catalyst", ev),
                             participant("ligand", "m_dpephos", "ligand", ev), participant("product", product, "product", ev, product_role=f"{selectivity} product")]
    r["classification"] = classification(ev, HYDROAMINATION)
    r["review"]["unresolved"].append("Generic scheme does not represent one concrete experiment or a complete R-group identity.")
    reactions.append(r)

# Three clearly labelled prior-work examples on page 2.
for rid, sub, prod, cat, title in [
    ("literature_menche_2010", "m_menche_sub", "m_menche_prod", "m_pdpph34", "Menche (2010): Allylic Substitution"),
    ("literature_white_2013", "m_white_sub", "m_white_prod", "m_pdoac2", "White (2013): Olefin Isomerization/Amination"),
    ("literature_hara_2015", "m_hara_sub", "m_hara_prod", "m_pdtfa2", "Hara (2015): C-H Oxidation"),
]:
    ev = [E[2]]
    r = base_reaction(rid, ev, "literature_example", [], "page 2 prior-work comparison", title, None)
    r["participants"] = [participant("substrate", sub, "reactant", ev), participant("catalyst", cat, "catalyst", ev), participant("product", prod, "product", ev)]
    r["classification"] = classification(ev, [change("b1", "C", "N", 0, 1, "formed")], limitation="Contextual prior-work drawing; only the evident C-N ring closure is classified, with no complete mapping.")
    r["review"]["unresolved"].append("Citation-context record only; no current-paper experimental quantities or yields are assigned.")
    reactions.append(r)

# Table 1: every optimization row is a separate experiment.
table = [
    (1, 80, "pd", "m_phenylisobutyric", 100, 85, "1:1"),
    (2, 80, "pd", "m_benzoic", 50, 99, "20:1"),
    (3, 80, "rh", "m_ppts", 30, 22, "1:2"),
    (4, 80, "rh", "m_benzoic", 50, 83, "1:3"),
    (5, 80, "rh", "m_pentafluoro", 20, 99, "1:5"),
    (6, 50, "rh", "m_pentafluoro", 20, 59, "1:8"),
    (7, 50, "rh_hi", "m_pentafluoro", 20, 95, "1:10"),
]
for no, temp, mode, additive, addmol, yld, dr in table:
    ev = [E[2], T[1] if no in (1, 2) else T[2]]
    catalyst = "m_pd_dba2" if mode == "pd" else "m_rh_cod_cl2"
    catload = 2 if mode != "rh_hi" else 4
    ligload = 5 if mode != "rh_hi" else 10
    general = [PDGEN if mode == "pd" else RHGEN]
    r = base_reaction(f"table1_entry_{no}", ev, "specific_experiment", ["optimization"], f"Table 1, entry {no}",
                      "optimization of anti-/syn-selective cyclization", "intramolecular hydroamination/cyclization", general)
    r["participants"] = [
        participant("substrate", "m_1a", "reactant", ev),
        participant("anti_product", "m_2a", "product", ev, product_role="anti diastereomer"),
        participant("syn_product", "m_3a", "product", ev, product_role="syn diastereomer"),
        participant("catalyst", catalyst, "catalyst", ev, qty(catload, "mol%", f"{catload} mol %"), "mol_percent", "m_1a"),
        participant("ligand", "m_dpephos", "ligand", ev, qty(ligload, "mol%", f"{ligload} mol %"), "mol_percent", "m_1a"),
        participant("additive", additive, "additive", ev, qty(addmol, "mol%", f"{addmol} mol %"), "mol_percent", "m_1a"),
    ]
    if no == 3:
        r["participants"].append(participant("licl", "m_licl", "additive", ev, qty(30, "mol%", "30 mol %"), "mol_percent", "m_1a"))
    set_solvent(r, "m_dce", "DCE", qty(0.2, "M", "0.2 M"), ev)
    add_step(r, "reaction", ev, [p["participant_id"] for p in r["participants"]], qty(temp, "°C", f"{temp} °C"), qty(16, "h", "16 h"))
    add_outcome(r, "yield", "yield", ["m_2a", "m_3a"], ev, f"{yld}%", qty(yld, "%", f"{yld}%"), basis="isolated", scope="combined products 2a + 3a", role="combined diastereomers")
    add_outcome(r, "dr", "selectivity", ["m_2a", "m_3a"], ev, dr, dr, basis="1H NMR of the crude product", scope="d.r. anti:syn", role="anti/syn product ratio")
    r["classification"] = classification(ev, HYDROAMINATION)
    r["review"]["unresolved"].append("Exact substrate charge, workup and purification are not reported in the supplied main article.")
    reactions.append(r)

# Scheme 1 Pd scope: each displayed product/yield/dr is a separate experiment.
pd_scope = {
    "b": (65, "12:1"), "c": (79, ">20:1"), "d": (82, ">20:1"), "e": (90, ">20:1"),
    "f": (99, ">20:1"), "g": (97, "20:1"), "h": (83, "20:1"), "i": (90, "20:1"),
    "j": (78, ">20:1"), "k": (86, ">20:1"), "l": (81, "10:1"), "m": (84, "20:1"),
    "n": (91, "13:1"), "o": (81, ">20:1"),
}
for x, (yld, dr) in pd_scope.items():
    ev = [E[2]]
    r = base_reaction(f"scheme1_pd_{x}", ev, "specific_experiment", ["substrate_scope"], f"Scheme 1, 1{x} to 2{x}",
                      "[Pd]-catalyzed anti-selective cyclization of tosylureas", "intramolecular hydroamination/cyclization", [PDGEN])
    r["participants"] = [
        participant("substrate", f"m_1{x}", "reactant", ev),
        participant("anti_product", f"m_2{x}", "product", ev, product_role="major anti product"),
        participant("syn_product", f"m_3{x}", "product", ev, product_role="minor syn product"),
        participant("catalyst", "m_pd_dba2", "catalyst", ev, qty(2, "mol%", "2.0 mol %"), "mol_percent", f"m_1{x}"),
        participant("ligand", "m_dpephos", "ligand", ev, qty(5, "mol%", "5.0 mol %"), "mol_percent", f"m_1{x}"),
        participant("additive", "m_benzoic", "additive", ev, qty(50, "mol%", "50 mol %"), "mol_percent", f"m_1{x}"),
    ]
    set_solvent(r, "m_dce", "DCE", qty(0.2, "M", "0.2 M"), ev)
    add_step(r, "reaction", ev, [p["participant_id"] for p in r["participants"]], qty(80, "°C", "80 °C"), qty(16, "h", "16 h"))
    add_outcome(r, "yield", "yield", [f"m_2{x}"], ev, f"{yld}%", qty(yld, "%", f"{yld}%"), basis="isolated", scope=f"product 2{x}")
    add_outcome(r, "dr", "selectivity", [f"m_2{x}", f"m_3{x}"], ev, dr, dr, basis="1H NMR of the crude product", scope="d.r. anti:syn")
    r["classification"] = classification(ev, HYDROAMINATION)
    reactions.append(r)

# Scheme 2 Rh scope.
rh_scope = {
    "b": (68, "9:1"), "c": (78, "10:1"), "d": (81, "7:1"), "e": (79, "6:1"),
    "f": (81, "7:1"), "g": (84, "7:1"), "h": (91, "7:1"), "i": (98, "6:1"),
    "j": (56, "9:1"), "k": (87, "6:1"), "l": (81, "6:1"), "m": (70, "5:1"),
    "n": (55, "2:1"), "o": (79, "3:1"),
}
for x, (yld, dr) in rh_scope.items():
    ev = [E[3], T[2]]
    r = base_reaction(f"scheme2_rh_{x}", ev, "specific_experiment", ["substrate_scope"], f"Scheme 2, 1{x} to 3{x}",
                      "[Rh]-catalyzed syn-selective cyclization of tosylureas", "intramolecular hydroamination/cyclization", [RHGEN])
    r["participants"] = [
        participant("substrate", f"m_1{x}", "reactant", ev),
        participant("syn_product", f"m_3{x}", "product", ev, product_role="major syn product"),
        participant("anti_product", f"m_2{x}", "product", ev, product_role="minor anti product"),
        participant("catalyst", "m_rh_cod_cl2", "catalyst", ev, qty(2, "mol%", "2.0 mol %"), "mol_percent", f"m_1{x}"),
        participant("ligand", "m_dpephos", "ligand", ev, qty(5, "mol%", "5.0 mol %"), "mol_percent", f"m_1{x}"),
        participant("additive", "m_pentafluoro", "additive", ev, qty(20, "mol%", "20 mol %"), "mol_percent", f"m_1{x}"),
    ]
    set_solvent(r, "m_dce", "DCE", qty(0.2, "M", "0.2 M"), ev)
    add_step(r, "reaction", ev, [p["participant_id"] for p in r["participants"]], qty(50, "°C", "50 °C"), qty(16, "h", "16 h"))
    add_outcome(r, "yield", "yield", [f"m_3{x}"], ev, f"{yld}%", qty(yld, "%", f"{yld}%"), basis="isolated", scope=f"product 3{x}")
    add_outcome(r, "dr", "selectivity", [f"m_3{x}", f"m_2{x}"], ev, dr, dr, basis="1H NMR of the crude products", scope="d.r. syn:anti")
    r["classification"] = classification(ev, HYDROAMINATION)
    reactions.append(r)


def mark_followup_conflict(r):
    a = f("Scheme 3 central substrate is labelled 2a", "m_2a", [E[3]])
    b = f("Prose states follow-up reactions were performed on 3a", "m_3a", [E[2], T[2]])
    conflict = f([a["raw_value"], b["raw_value"]], None, [E[2], E[3], T[2]], status="conflict", alternatives=[a, b])
    r["notes"].append(conflict)
    r["review"]["unresolved"].append("Source conflict: Scheme 3 explicitly labels the starting structure 2a, while the preceding prose says the follow-up set was performed on 3a. The structured participant follows the scheme and the conflict is not silently resolved.")


# Scheme 3 transformations. Scheme label 2a is used structurally, while the prose conflict is preserved.
followups = []
for rid, product, yld, name in [
    ("scheme3_a_detosylation", "m_4", 97, "reductive cleavage of the tosyl group"),
    ("scheme3_ab_diamine", "m_5", 72, "two-step detosylation and hydrolysis to diamine dihydrochloride"),
    ("scheme3_c_methylation", "m_9", 68, "N-methylation"),
    ("scheme3_d_hydrogenation", "m_6", 86, "hydrogenation"),
    ("scheme3_e_ozonolysis", "m_7", 77, "ozonolysis and reductive workup"),
    ("scheme3_f_hydroformylation", "m_8", 76, "hydroformylation"),
]:
    ev = [E[3], T[2], T[3]]
    r = base_reaction(rid, ev, "specific_experiment", ["derivatization"], f"Scheme 3, arrow {rid.split('_')[1]}", name, name)
    r["participants"] = [participant("substrate", "m_2a", "reactant", ev), participant("product", product, "product", ev, product_role="main product")]
    mark_followup_conflict(r)
    followups.append((r, product, yld))

# a: Mg/MeOH.
r, product, yld = followups[0]
r["participants"] += [participant("mg", "m_mg", "reagent", [E[3]], qty(20, "equiv", "20.0 equiv"), "equivalents", "m_2a")]
set_solvent(r, "m_meoh", "MeOH", qty(0.05, "M", "0.05 M"), [E[3]])
add_step(r, "reaction", [E[3]], [p["participant_id"] for p in r["participants"]], "rt", qty(16, "h", "16 h"))
add_outcome(r, "yield", "yield", [product], [E[3]], "97%", qty(97, "%", "97%"), basis=None, scope="2a/3a source-label conflict to 4")
r["classification"] = classification([E[3]], [change("b1", "N", "S", 1, 0, "broken"), change("b2", "N", "H", 0, 1, "hydrogen_gain")])

# a,b: includes both sequential stages; 72% over two steps.
r, product, yld = followups[1]
r["participants"] += [participant("mg", "m_mg", "reagent", [E[3]], qty(20, "equiv", "20.0 equiv"), "equivalents", "m_2a"),
                       participant("hcl", "m_hcl", "reagent", [E[3]], qty(10, "mL", "10 mL/mmol"), "volume", "m_2a")]
set_solvent(r, "m_meoh", "MeOH", qty(0.05, "M", "0.05 M"), [E[3]])
add_step(r, "step_a_detosylation", [E[3]], ["substrate", "mg", "p_5_m_meoh"], "rt", qty(16, "h", "16 h"))
add_step(r, "step_b_acid_hydrolysis", [E[3]], ["hcl", "product"], qty(100, "°C", "100 °C"), qty(3, "h", "3 h"))
r["steps"][-1]["operation"] = f("aq. 2 M HCl (10 mL/mmol)", {"operation": "acid hydrolysis", "HCl_concentration": {"number": 2, "unit": "M"}, "dose": "10 mL/mmol"}, [E[3]])
add_outcome(r, "overall_yield", "yield", [product], [E[3]], "72% over two steps", qty(72, "%", "72%"), basis=None, yield_type="two-step overall", scope="starting cyclic urea to 1,3-diamine 5 dihydrochloride")
r["classification"] = classification([E[3]], [change("b1", "N", "S", 1, 0, "broken"), change("b2", "N", "C", 1, 0, "broken"), change("b3", "N", "H", 0, 1, "hydrogen_gain")])
r["review"]["unresolved"].append("The printed 10 mL/mmol is retained as a volume normalized per substrate amount; no absolute volume or scale is reported.")

# c: methylation.
r, product, yld = followups[2]
r["participants"] += [participant("base", "m_k2co3", "reagent", [E[3]], qty(5, "equiv", "5.0 equiv"), "equivalents", "m_2a"),
                       participant("methylating_agent", "m_me2so4", "reagent", [E[3]], qty(1.5, "equiv", "1.5 equiv"), "equivalents", "m_2a")]
set_solvent(r, "m_mecn", "MeCN", qty(0.05, "M", "0.05 M"), [E[3]])
add_step(r, "reaction", [E[3]], [p["participant_id"] for p in r["participants"]], qty(70, "°C", "70 °C"), qty(16, "h", "16 h"))
add_outcome(r, "yield", "yield", [product], [E[3]], "68%", qty(68, "%", "68%"), basis=None)
r["classification"] = classification([E[3]], [change("b1", "C", "N", 0, 1, "formed")])

# d: hydrogenation.
r, product, yld = followups[3]
r["participants"] += [participant("hydrogen", "m_h2", "reagent", [E[3]]), participant("catalyst", "m_pdc", "catalyst", [E[3]])]
r["participants"][-1]["charges"]["amount"] = f("10 w %", "10 w %", [E[3]])
add_step(r, "reaction", [E[3]], [p["participant_id"] for p in r["participants"]], "rt", qty(18, "h", "18 h"))
r["steps"][0]["pressure"] = f("10 bar", "10 bar", [E[3]])
add_outcome(r, "yield", "yield", [product], [E[3]], "86%", qty(86, "%", "86%"), basis=None)
r["classification"] = classification([E[3]], [change("b1", "C", "C", 2, 1, "order_decreased"), change("b2", "C", "H", 0, 1, "hydrogen_gain", 2)])
r["review"]["unresolved"].append("Pd/C loading is printed as 10 w%; the exact denominator wording is not expanded in the article.")

# e: ozonolysis then TPP workup, separate temperatures and times.
r, product, yld = followups[4]
r["participants"] += [participant("ozone", "m_o3", "reagent", [E[3]]), participant("tpp", "m_tpp", "workup_agent", [E[3]], qty(5, "equiv", "5 equiv"), "equivalents", "m_2a", stage="workup")]
r["participants"] += [participant("meoh", "m_meoh", "solvent", [E[3]]), participant("dcm", "m_dcm", "solvent", [E[3]])]
r["solvent_conditions"]["composition"] = f("MeOH/DCM", ["m_meoh", "m_dcm"], [E[3]])
r["solvent_conditions"]["mixture_ratio"] = f("5:1", "5:1", [E[3]])
r["solvent_conditions"]["concentration"] = f("0.025 M", qty(0.025, "M", "0.025 M"), [E[3]])
r["solvent_conditions"]["concentration_basis"] = f("substrate concentration", "substrate", [E[3]])
add_step(r, "ozonolysis", [E[3]], ["substrate", "ozone", "meoh", "dcm"], qty(-78, "°C", "−78 °C"), qty(10, "min", "10 min"), "ozonolysis")
s = st.new_step("reductive_workup", [E[3]], stage="workup"); s["participant_ids"] = ["tpp", "product"]; s["operation"] = f("then TPP", "reductive workup with TPP", [E[3]]); s["temperature"] = f("rt", "rt", [E[3]]); s["time"] = f("1 h", qty(1, "h", "1 h"), [E[3]]); r["steps"].append(s)
add_outcome(r, "yield", "yield", [product], [E[3]], "77%", qty(77, "%", "77%"), basis=None)
r["classification"] = classification([E[3]], [change("b1", "C", "O", 0, 2, "formed"), change("b2", "C", "C", 2, 0, "broken")])

# f: hydroformylation.
r, product, yld = followups[5]
r["participants"] += [participant("hydrogen", "m_h2", "reagent", [E[3]]), participant("co", "m_co", "reagent", [E[3]]),
                       participant("catalyst", "m_rh_co2acac", "catalyst", [E[3]], qty(0.5, "mol%", "0.5 mol %"), "mol_percent", "m_2a"),
                       participant("ligand", "m_6dppon", "ligand", [E[3]], qty(5, "mol%", "5 mol %"), "mol_percent", "m_2a")]
set_solvent(r, "m_toluene", "toluene", qty(0.025, "M", "0.025 M"), [E[3]])
add_step(r, "reaction", [E[3]], [p["participant_id"] for p in r["participants"]], qty(80, "°C", "80 °C"), qty(18, "h", "18 h"), "hydroformylation under H2/CO (1:1)")
r["steps"][0]["pressure"] = f("20 bar", "20 bar", [E[3]])
add_outcome(r, "yield", "yield", [product], [E[3], T[3]], "76% isolated yield", qty(76, "%", "76%"), basis="isolated", scope="linear hydroformylation product 8")
add_outcome(r, "regioselectivity", "selectivity", [product], [E[3]], "97lin:3br", "97:3", basis=None, scope="linear:branched")
r["classification"] = classification([E[3]], [change("b1", "C", "C", 0, 1, "formed"), change("b2", "C", "O", 0, 2, "formed"), change("b3", "C", "H", 0, 1, "hydrogen_gain"), change("b4", "C", "C", 2, 1, "order_decreased")])

reactions.extend(r for r, _, _ in followups)

# Scheme 4 control experiments.
ev = [E[3], T[3]]
r = base_reaction("scheme4_pd_isomerization", ev, "specific_experiment", ["control", "mechanistic_study"], "Scheme 4, 3a under Pd conditions", "Pd-catalyzed stereoisomerization", "stereochemical isomerization", [PDGEN])
r["participants"] = [participant("substrate", "m_3a", "reactant", ev), participant("product", "m_2a", "product", ev, product_role="anti product"),
                     participant("catalyst", "m_pd_dba2", "catalyst", ev, qty(2, "mol%", "2.0 mol %"), "mol_percent", "m_3a"),
                     participant("ligand", "m_dpephos", "ligand", ev, qty(5, "mol%", "5 mol %"), "mol_percent", "m_3a"),
                     participant("additive", "m_benzoic", "additive", ev, qty(50, "mol%", "50 mol %"), "mol_percent", "m_3a")]
set_solvent(r, "m_dce", "DCE", qty(0.2, "M", "0.2 M"), ev)
add_step(r, "reaction", ev, [p["participant_id"] for p in r["participants"]], qty(80, "°C", "80 °C"), qty(16, "h", "16 h"))
add_outcome(r, "yield", "yield", ["m_2a"], ev, "quant.", "quantitative", basis=None)
add_outcome(r, "dr", "selectivity", ["m_2a", "m_3a"], ev, "d.r. 96:4 (anti:syn)", "96:4", basis=None, scope="d.r. anti:syn")
r["classification"] = classification(ev, coverage="not_applicable", limitation="Observed interconversion changes relative stereochemistry but not molecular connectivity; no new bond-formation label is assigned.")
reactions.append(r)

r = base_reaction("scheme4_rh_no_reaction", ev, "specific_experiment", ["control", "mechanistic_study"], "Scheme 4, 2a under Rh conditions", "Rh-condition control", "no reaction", [RHGEN])
r["participants"] = [participant("substrate", "m_2a", "reactant", ev),
                     participant("catalyst", "m_rh_cod_cl2", "catalyst", ev, qty(4, "mol%", "4 mol %"), "mol_percent", "m_2a"),
                     participant("ligand", "m_dpephos", "ligand", ev, qty(10, "mol%", "10 mol %"), "mol_percent", "m_2a"),
                     participant("additive", "m_pentafluoro", "additive", ev, qty(20, "mol%", "20 mol %"), "mol_percent", "m_2a")]
set_solvent(r, "m_dce", "DCE", qty(0.2, "M", "0.2 M"), ev)
add_step(r, "reaction", ev, [p["participant_id"] for p in r["participants"]], qty(50, "°C", "50 °C"), qty(16, "h", "16 h"))
add_outcome(r, "no_reaction", "no_reaction", [], ev, "no reaction / no change in diastereomeric ratio", "no reaction", basis=None, scope="2a under Rh conditions")
r["classification"] = classification(ev, coverage="not_applicable", limitation="Negative control: no expected product or bond formation is invented.")
reactions.append(r)

# Scheme 5 is a proposal, not an experiment.
ev = [E[3], E[4], T[3], T[4]]
r = base_reaction("scheme5_proposed_mechanism", ev, "proposed_mechanism", [], "Scheme 5 and continuation on page 4", "proposed mechanism for Rh/Pd-catalyzed intramolecular addition of tosylureas to allenes", None)
r["participants"] = [participant("substrate", "m_1a", "reactant", ev), participant("kinetic_product", "m_3a", "product", ev, product_role="kinetic syn product"), participant("thermodynamic_product", "m_2a", "product", ev, product_role="thermodynamic anti product")]
r["classification"] = classification(ev, coverage="not_applicable", limitation="Arrow-pushing intermediates A-C are proposed and are not treated as isolated compounds or an experimentally atom-mapped reaction.")
r["notes"].append(f("Rh pathway: hydrometalation to cyclic pi-allyl B and inner-sphere reductive elimination to syn-3a; Pd pathway: open pi-allyl C and outer-sphere attack to anti-2a.",
                    "Author-proposed mechanistic interpretation", ev))
r["review"]["unresolved"].append("Mechanistic intermediates A-C and elementary steps are proposals, not separately verified experiments.")
reactions.append(r)

bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": f"{DOC}:metadata", "document_id": DOC, "kind": "publication_metadata",
     "title": "Stereodivergent Palladium- and Rhodium-Catalyzed Intramolecular Addition of Tosylureas to Allenes: Diastereoselective Synthesis of Tetrahydropyrimidinones",
     "authors": ["Arne G. A. Geissler", "Jasmin R. Riesterer", "Bernhard Breit"],
     "journal": "Organic Letters", "year": 2021, "volume": 23, "pages": "9168-9172", "doi": "10.1021/acs.orglett.1c03482",
     "evidence_ids": [E[1], T[1]]},
    {"context_id": f"{DOC}:dft", "document_id": DOC, "kind": "computational_interpretation",
     "raw_value": "DFT (BP86/dev2svp) predicts 2a to be 2.4 kcal more stable than 3a; 2a thermodynamic and 3a kinetic.",
     "evidence_ids": [E[3], T[3]]},
    {"context_id": f"{DOC}:si_advertised", "document_id": DOC, "kind": "supplement_availability",
     "raw_value": "The article advertises Supporting Information with synthetic procedures, analytical data and computational data; no SI file was supplied or explicitly associated in this workspace.",
     "evidence_ids": [E[4], T[4]]},
    {"context_id": f"{DOC}:pages4_5", "document_id": DOC, "kind": "non_experiment_pages",
     "raw_value": "Page 4 contains mechanism discussion, SI notice, author information and references; page 5 contains references only. No standalone spectral table or new synthesis arrow was treated as an experiment.",
     "evidence_ids": [E[4], E[5]]},
]

# Machine-readable raw transcription kept separate from the normalized submission.
raw = {
    "document_id": DOC,
    "pages_viewed": [1, 2, 3, 4, 5],
    "table1": [{"entry": n, "temperature_C": t, "catalyst_mode": m, "additive": a, "additive_mol_percent": ap, "isolated_yield_percent": y, "dr_anti_syn": d}
               for n, t, m, a, ap, y, d in table],
    "scheme1_pd_scope": pd_scope,
    "scheme2_rh_scope": rh_scope,
    "scheme3": {"a": {"product": "4", "yield": "97%"}, "a,b": {"product": "5", "yield": "72% over two steps"},
                "c": {"product": "9", "yield": "68%"}, "d": {"product": "6", "yield": "86%"},
                "e": {"product": "7", "yield": "77%"}, "f": {"product": "8", "yield": "76%, 97lin:3br"}},
    "source_conflicts": [{"field": "Scheme 3 starting material", "scheme": "2a", "preceding_prose": "3a"}],
    "record_count": len(reactions),
    "molecule_count": len(molecules),
}
RAW.mkdir(parents=True, exist_ok=True)
st.write(RAW / "visual_transcription.json", raw)
st.write(RAW / "stereo_checks.json", stereo_checks)
st.write(OUT, bundle)
print(json.dumps({"molecules": len(molecules), "reactions": len(reactions)}, ensure_ascii=False))
