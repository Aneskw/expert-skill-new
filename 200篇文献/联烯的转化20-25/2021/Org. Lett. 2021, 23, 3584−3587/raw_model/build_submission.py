#!/usr/bin/env python3
"""Source-bounded v3.1 transcription for DOI 10.1021/acs.orglett.1c00981.

This is a machine-authored submission builder.  It preserves unknowns and is
not a claim of source review or chemical correctness.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st  # noqa: E402

RUN = ROOT / "runs/current"
DATASET = RUN / "dataset.json"
DOC = "doc_275edf2145e01c65"
E = {n: f"{DOC}:page{n}:image" for n in range(1, 5)}
T = {n: f"{DOC}:page{n}:text" for n in range(1, 5)}


def f(raw=None, normalized=None, ev=(), method="original_extraction", derivation=None,
      status=None):
    return st.fact(raw, normalized, ev, method=method, derivation=derivation,
                   status=status)


def derived(normalized, ev, why, raw=None):
    return f(raw, normalized, ev, method="structure_derivation", derivation=why)


def inferred(normalized, ev, why, raw=None):
    return f(raw, normalized, ev, method="inference", derivation=why)


def na(ev, why="Not applicable to this molecule or experiment."):
    return f(ev=ev, method="not_applicable", status="not_applicable",
             derivation=why)


def q(number, unit, **extra):
    return {"number": number, "unit": unit, **extra}


def molecule(mid, label, ev, *, name=None, smiles=None, chirality="unresolved",
             geometry="unresolved", representation=None, stereo_note=None,
             double_note=None, axial_note=None, identity_note=None):
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    m["names"]["original"] = f(name or label, name or label, ev)
    m["names"]["english"] = f(name, name, ev) if name else f(ev=ev)
    m["cas_numbers"] = f(ev=ev)
    m["identity_note"] = f(identity_note, identity_note, ev) if identity_note else f(ev=ev)
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = derived(smiles, ev,
            "Candidate connectivity transcribed from the supplied source drawing; RDKit parsing does not verify source identity.",
            raw={"source_drawing": label})
        diag = st.stereo_inspect(smiles)
        if representation is None:
            representation = "partial" if any(
                p["specified"] == "Unspecified" for p in diag["potential_stereo"]
            ) else "encoded"
        if representation == "encoded":
            m["stereochemistry"]["isomeric_smiles"] = derived(
                diag["canonical_isomeric_smiles"], ev,
                "RDKit canonicalization of the candidate source-derived graph.", raw=smiles)
            m["stereochemistry"]["identity_key"] = derived(
                diag["stereo_graph_key"], ev,
                "SHA-256 key of the canonical isomeric candidate graph; not a sample identity.", raw=smiles)
        else:
            m["stereochemistry"]["isomeric_smiles"] = f(ev=ev)
            m["stereochemistry"]["identity_key"] = f(ev=ev)
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = f(ev=ev)
        m["generic_template"] = None
        representation = representation or "not_encoded"
    m["stereochemistry"]["chirality_status"] = inferred(
        chirality, ev, "Status assigned from explicit wedges/dashes or from the candidate molecular graph; absence of a source label was not treated as racemic.", raw=stereo_note or chirality)
    m["stereochemistry"]["geometry_status"] = inferred(
        geometry, ev, "Status assigned from the source drawing; missing E/Z labels remain unspecified.", raw=double_note or geometry)
    m["stereochemistry"]["representation_status"] = inferred(
        representation, ev, "Compares source stereochemical content with the candidate representation.", raw=representation)
    m["stereochemistry"]["atom_centers"] = (
        f(stereo_note, stereo_note, ev) if stereo_note else f(ev=ev)
    )
    m["stereochemistry"]["double_bonds"] = (
        f(double_note, double_note, ev) if double_note else f(ev=ev)
    )
    m["stereochemistry"]["axial"] = (
        f(axial_note, axial_note, ev) if axial_note else f(ev=ev)
    )
    m["stereochemistry"]["mixture"] = f(ev=ev)
    m["stereochemistry"]["assignment"] = (
        f(stereo_note, stereo_note, ev) if stereo_note else f(ev=ev)
    )
    return m


def template_molecule(mid, label, ev, fixed, variables, reason, chirality, geometry):
    m = molecule(mid, label, ev, chirality=chirality, geometry=geometry,
                 representation="not_encoded", identity_note=reason,
                 stereo_note=("The generic product drawing contains wedge/dash-defined relative stereochemistry; no absolute descriptors are printed."
                              if chirality == "specified" else None))
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": f(fixed, fixed, ev),
        "variables": [], "instance_mapping": [],
        "reason": f(reason, reason, ev),
    }
    for var in variables:
        m["generic_template"]["variables"].append({
            "label": var[0],
            "attachment_points": f(var[1], var[1], ev),
            "definition": f(var[2], var[2], ev),
            "allowed_values": f(var[3], var[3], ev),
            "constraints": f(var[4], var[4], ev),
        })
    return m


def participant(pid, mid, role, ev, *, stage="reaction", eq=None, mol=None,
                reference=None, product_role=None):
    p = st.new_participant(pid, mid, role, ev, stage=stage)
    if eq is not None:
        p["charges"]["equivalents"] = f(f"{eq} equiv", q(eq, "equiv"), ev)
    if mol is not None:
        p["charges"]["mol_percent"] = f(f"{mol} mol%", q(mol, "mol%"), ev)
    if reference:
        p["equivalence_reference"] = f(reference, reference, ev)
    if product_role:
        p["product_role"] = f(product_role, product_role, ev)
    return p


def step(sid, ev, pids, *, operation, temp=None, time=None, stage="reaction"):
    s = st.new_step(sid, ev, stage=stage)
    s["participant_ids"] = list(pids)
    s["operation"] = f(operation, operation, ev)
    if temp is not None:
        s["temperature"] = f(temp[0], temp[1], ev)
    if time is not None:
        s["time"] = f(time[0], time[1], ev)
    return s


def outcome(oid, kind, targets, ev, *, raw, value, product_role=None,
            yield_type=None, basis=None, scope=None):
    o = st.new_outcome(oid, kind, targets, ev)
    o["value"] = f(raw, value, ev)
    if product_role:
        o["product_role"] = f(product_role, product_role, ev)
    if yield_type:
        o["yield_type"] = f(yield_type, yield_type, ev)
    if basis:
        o["measurement_basis"] = f(basis, basis, ev)
    if scope:
        o["scope"] = f(scope, scope, ev)
    return o


def classification(ev, changes=None, coverage="partial", limitation=None):
    changes = changes or []
    c = {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": f(ev=ev),
        "bond_changes": [],
        "labels": f(ev=ev),
        "coverage": inferred(coverage, ev,
            "Coverage reflects only source-supported net bond classes; no complete atom mapping was asserted.", raw=coverage),
        "limitations": [limitation or "No atom-mapped complete structures were asserted; bond changes are partial structural deductions from the source schemes."],
        "model_votes": [],
    }
    if coverage in ("unresolved", "not_applicable"):
        if coverage == "not_applicable":
            c["labels"] = inferred([], ev, "No observed bond-forming transformation is classified for this record.", raw=[])
        return c
    labels = {}
    for i, ch in enumerate(changes, 1):
        change = {"change_id": f"b{i}", **ch}
        c["bond_changes"].append(derived(change, ev,
            "Net change inferred from the source reaction drawing; atom identities are element-level because a complete mapping is unavailable."))
        if change["change_type"] in ("formed", "hydrogen_gain", "order_increased"):
            lab = st.bond_label(change["atoms"][0]["element"], change["atoms"][1]["element"], change["after_order"])
            labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(change["change_id"])
    c["labels"] = derived(list(labels.values()), ev,
        "Taxonomy labels calculated from the listed partial net bond changes.")
    return c


def formed(a, b, order=1):
    return {"atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
            "before_order": 0, "after_order": order, "change_type": "formed"}


def hgain(a, count=1):
    return {"atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": "H"}],
            "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": count}


def order_down(a, b, before, after):
    return {"atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
            "before_order": before, "after_order": after, "change_type": "order_decreased"}


def broken(a, b, before=1):
    return {"atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
            "before_order": before, "after_order": 0, "change_type": "broken"}


def reaction(rid, ev, *, categories, location, name, transformation, parts,
             steps, outcomes, changes=None, coverage="partial", notes=(),
             unresolved=(), general=(), record_type="specific_experiment",
             solvent=None, concentration=None, limitation=None):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, ev)
    r["location"] = f(location, location, ev)
    r["general_scheme_ids"] = list(general)
    r["reaction_name"] = f(name, name, ev)
    r["transformation_type"] = inferred(transformation, ev,
        "Normalized transformation type from the authors' reaction name and depicted net change.", raw=name)
    r["participants"] = parts
    r["steps"] = steps
    r["outcomes"] = outcomes
    if solvent:
        r["solvent_conditions"]["composition"] = f(solvent, solvent, ev)
    if concentration:
        r["solvent_conditions"]["concentration"] = f(
            concentration[0], concentration[1], ev)
        r["solvent_conditions"]["concentration_basis"] = f(
            "scheme concentration; basis not stated", "basis not stated", ev)
    r["procedure"]["original_text"] = f(location, location, ev)
    r["classification"] = classification(ev, changes, coverage, limitation)
    r["notes"] = [f(x, x, ev) for x in notes]
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": list(unresolved)}
    r["completeness"] = "partial"
    return r


def add_solvent(parts, mid, ev):
    parts.append(participant(f"p_{mid}", mid, "solvent", ev))


def add_standard_catalysis(parts, ev, substrate, ligand, *, pd=2.5, lig=5.0,
                           base=0.1, alcohol="m2", alcohol_eq=1.2, solvent="m_ch2cl2"):
    parts.extend([
        participant("p_substrate", substrate, "reactant", ev),
        participant("p_alcohol", alcohol, "reactant", ev, eq=alcohol_eq,
                    reference=substrate) if alcohol_eq is not None else participant("p_alcohol", alcohol, "reactant", ev),
        participant("p_pd", "m_pd2dba3", "catalyst", ev, mol=pd, reference=substrate),
        participant("p_ligand", ligand, "ligand", ev, mol=lig, reference=substrate),
        participant("p_base", "m_et3n", "additive", ev, eq=base, reference=substrate),
        participant("p_solvent", solvent, "solvent", ev),
    ])


data = st.read(DATASET)
data["documents"][0]["title"] = "Catalytic Asymmetric Synthesis of Hexahydro-furofuran-3-ol and Its Pyran Derivatives"
st.write(DATASET, data)
bundle = st.empty(data["documents"])

# Primary numbered compounds.  Candidate SMILES are used only where the
# connectivity can be conservatively transcribed; dense acyclic acetals remain
# unresolved rather than inventing a complete graph.
primary = [
    ("m1", "1", [E[2]], "ene-alkoxyallene 1", "CCC=CCOC=C=C", "achiral", "unspecified", "partial", None, "The side-chain alkene has no printed E/Z assignment; terminal allene is not treated as axially assigned."),
    ("m2", "2", [E[2]], "2-bromo allylic alcohol 2", "OCC(Br)=C", "achiral", "not_applicable", "encoded", None, None),
    ("m3", "3", [E[2]], "adduct 3", None, "specified", "unspecified", "not_encoded", "A dashed stereobond is printed in the product drawing; absolute R/S is not printed.", None),
    ("m4", "4", [E[2]], "cyclic acetal 4", "C=C(Br)OC1OCC=C1", "specified", "not_applicable", "partial", "The cyclic acetal drawing contains a dashed bond; absolute R/S is not printed.", None),
    ("m5", "5", [E[2]], "bis-benzoate 5", None, "specified", "not_applicable", "not_encoded", "Wedge/dash bonds define relative configuration in the drawing; no R/S descriptors are printed.", None),
    ("m6", "6", [E[2]], "bicyclic exo-olefin 6", "C=C1CC2OCC(O2)C1", "specified", "not_applicable", "partial", "Several wedge/dash bonds define the fused-ring relative configuration; not encoded in the candidate SMILES.", None),
    ("m7", "7", [E[2]], "bicyclic ketone 7", "O=C1CC2OCC(O2)C1", "specified", "not_applicable", "partial", "Several wedge/dash bonds define the fused-ring relative configuration; not encoded in the candidate SMILES.", None),
    ("m8", "8", [E[2]], "bis-hexahydro-furofuranol 8", "OC1CC2OCC(O2)C1", "specified", "not_applicable", "partial", "The source depicts a single relative configuration and reports complete diastereoselectivity; absolute R/S descriptors are not printed.", None),
    ("m9", "9", [E[2]], "benzoate 9", "O=C(OC1CC2OCC(O2)C1)c1ccccc1", "specified", "not_applicable", "partial", "The fused-ring and benzoate stereochemistry is shown by wedges/dashes but is not encoded in the candidate SMILES.", None),
    ("m10", "10", [E[2]], "ene-alkoxyallene 10", "CCC=CCCOC=C=C", "achiral", "unspecified", "partial", None, "The side-chain alkene has no printed E/Z assignment; terminal allene is not treated as axially assigned."),
    ("m11", "11", [E[2]], "acetal intermediate 11", None, "specified", "unspecified", "not_encoded", "A dashed stereobond is printed; no absolute descriptor is printed.", None),
    ("m12", "12", [E[2]], "cyclic acetal 12", "C=C(Br)OC1OCCC=C1", "specified", "not_applicable", "partial", "The acetal stereochemistry is drawn but not encoded in the candidate SMILES.", None),
    ("m13", "13", [E[2]], "bicyclic exo-olefin 13", "C=C1CC2OCCC(O2)C1", "specified", "not_applicable", "partial", "Wedge/dash fused-ring stereochemistry is not encoded in the candidate SMILES.", None),
    ("m14", "14", [E[2], E[3]], "bicyclic ketone 14", "O=C1CC2OCCC(O2)C1", "specified", "not_applicable", "partial", "Wedge/dash fused-ring stereochemistry is not encoded in the candidate SMILES.", None),
    ("m15", "15", [E[2], E[3]], "THF-THP alcohol 15", "OC1CC2OCCC(O2)C1", "specified", "not_applicable", "partial", "The source depicts one diastereomer; the relative configuration is not encoded in the candidate SMILES.", None),
    ("m16", "16", [E[2], E[3]], "benzoate 16", "O=C(OC1CC2OCCC(O2)C1)c1ccccc1", "specified", "not_applicable", "partial", "The source stereochemistry is not encoded; the 90% ee does not establish printed R/S descriptors.", None),
    ("m17", "17", [E[3]], "homoallylic alcohol 17", "OCCC(Br)=C", "achiral", "not_applicable", "encoded", None, None),
    ("m18", "18", [E[3]], "acetal intermediate 18", None, "specified", "unspecified", "not_encoded", "A dashed stereobond is printed; no absolute descriptor is printed.", None),
    ("m19", "19", [E[3]], "cyclic acetal vinyl halide 19", "C=C(Br)COC1OCC=C1", "specified", "not_applicable", "partial", "The cyclic acetal stereochemistry is drawn but not encoded in the candidate SMILES.", None),
    ("m20", "20", [E[3]], "diene 20", None, "specified", "not_applicable", "not_encoded", "Wedge/dash fused-ring stereochemistry is printed; the dense diene graph was not safely encoded.", None),
    ("m21", "21", [E[3]], "exo-olefin 21", "C=C1CCC2OCC(O2)C1", "specified", "not_applicable", "partial", "Wedge/dash fused-ring stereochemistry is not encoded in the candidate SMILES.", None),
    ("m22", "22", [E[3]], "THP-THF alcohol 22", "OC1CCC2OCC(O2)C1", "specified", "not_applicable", "partial", "The source depicts one relative configuration; it is not encoded in the candidate SMILES.", None),
    ("m23", "23", [E[3]], "benzoate 23", "O=C(OC1CCC2OCC(O2)C1)c1ccccc1", "specified", "not_applicable", "partial", "The source stereochemistry is not encoded; 95% ee does not supply a printed absolute descriptor.", None),
]
for row in primary:
    mid, label, ev, name, smi, chi, geo, rep, snote, dnote = row
    bundle["molecules"].append(molecule(mid, label, ev, name=name, smiles=smi,
        chirality=chi, geometry=geo, representation=rep, stereo_note=snote,
        double_note=dnote,
        identity_note="Candidate connectivity only; source wedges/dashes remain authoritative." if rep in ("partial", "not_encoded") else None))

# Reagents, catalysts, ligands and solvents.
aux = [
 ("m_pd2dba3", "Pd2(dba)3", "tris(dibenzylideneacetone)dipalladium(0)", None, [E[2],E[3]], "unresolved"),
 ("m_l1", "L1", "chiral ligand L1", None, [E[2]], "specified"),
 ("m_l2", "L2", "chiral ligand L2", None, [E[2],E[3]], "specified"),
 ("m_et3n", "Et3N", "triethylamine", "CCN(CC)CC", [E[2],E[3]], "achiral"),
 ("m_toluene", "toluene", "toluene", "Cc1ccccc1", [E[2]], "achiral"),
 ("m_ch2cl2", "CH2Cl2", "dichloromethane", "ClCCl", [E[2],E[3]], "achiral"),
 ("m_grubbs1", "Grubbs 1st", "first-generation Grubbs catalyst", None, [E[2],E[3]], "unresolved"),
 ("m_oso4", "OsO4", "osmium tetroxide", "O=[Os](=O)(=O)=O", [E[2]], "achiral"),
 ("m_bzcl", "BzCl", "benzoyl chloride", "O=C(Cl)c1ccccc1", [E[2],E[3]], "achiral"),
 ("m_pyridine", "pyridine", "pyridine", "n1ccccc1", [E[2],E[3]], "achiral"),
 ("m_bu3snh", "Bu3SnH", "tributyltin hydride", "[H][Sn](CCCC)(CCCC)CCCC", [E[2],E[3]], "achiral"),
 ("m_o2", "O2", "oxygen", "O=O", [E[2],E[3]], "achiral"),
 ("m_bet3", "BEt3", "triethylborane", "CCB(CC)CC", [E[2],E[3]], "achiral"),
 ("m_o3", "O3", "ozone", "[O-][O+]=O", [E[2],E[3]], "achiral"),
 ("m_dms", "(CH3)2S", "dimethyl sulfide", "CSC", [E[2],E[3]], "achiral"),
 ("m_nabh4", "NaBH4", "sodium borohydride", "[BH4-].[Na+]", [E[2],E[3]], "achiral"),
 ("m_etoh", "EtOH", "ethanol", "CCO", [E[2],E[3]], "achiral"),
 ("m_dppp", "dppp", "1,3-bis(diphenylphosphino)propane", "c1ccc(P(CCCP(c2ccccc2)c2ccccc2)c2ccccc2)cc1", [E[3]], "achiral"),
 ("m_dmf", "DMF", "N,N-dimethylformamide", "CN(C)C=O", [E[3]], "achiral"),
 ("m_pdc", "Pd/C", "palladium on carbon", None, [E[3]], "unresolved"),
 ("m_h2", "H2", "hydrogen", "[H][H]", [E[3]], "achiral"),
 ("m_bipy", "2,2'-bipyridyl", "2,2'-bipyridyl", "n1ccccc1-c1ccccn1", [E[3]], "achiral"),
 ("m_meoh", "MeOH", "methanol", "CO", [E[3]], "achiral"),
 ("m_grubbs2", "Grubbs second", "second-generation Grubbs catalyst", None, [E[4]], "unresolved"),
 ("m_hg2", "Hoveyda Grubbs second", "second-generation Hoveyda-Grubbs catalyst", None, [E[4]], "unresolved"),
]
for mid,label,name,smi,ev,chi in aux:
    bundle["molecules"].append(molecule(mid,label,ev,name=name,smiles=smi,
        chirality=chi,geometry="not_applicable",representation="encoded" if smi else "not_encoded",
        stereo_note=("The ligand drawing is chiral but no R/S descriptors are printed." if chi=="specified" else None),
        identity_note=("Named heterogeneous or organometallic material; no unique complete molecular SMILES asserted." if not smi else None)))

bundle["molecules"].append(molecule("m_dehal19", "dehalogenated product from 19", [E[3],T[4]],
    name="dehalogenated product from 19", chirality="unresolved", geometry="unresolved",
    representation="not_encoded", identity_note="Footnote 16 reports only a small amount; no structure or compound number is supplied."))

# Generic current-work family from graphical abstract and Scheme 1.
bundle["molecules"].extend([
 template_molecule("m_g_allene", "generic ene-alkoxyallene", [E[1],E[2]],
    "R-substituted ene-alkoxyallene with n = 1 or 2", [("R","allene terminus","unspecified substituent R","not enumerated","must match a depicted instance"),("n","ene-alkoxy tether","integer tether length","1 or 2","source-defined")],
    "Generic R/n drawing cannot be represented as a concrete molecule.", "unspecified", "unresolved"),
 template_molecule("m_g_alcohol", "generic bromo-olefinic alcohol", [E[1],E[2]],
    "bromo-olefinic alcohol with X and m = 1 or 2", [("X","vinylic position","halogen label X","not enumerated","source-defined"),("m","hydroxyalkyl tether","integer tether length","1 or 2","source-defined")],
    "Generic X/m drawing cannot be represented as a concrete molecule.", "unspecified", "unresolved"),
 template_molecule("m_g_product", "generic fused bicyclic alcohol", [E[1],E[2]],
    "fused oxygen heterocycle alcohol with n,m = 1 or 2", [("n","left ring","ring-size variable","1 or 2","source-defined"),("m","right ring","ring-size variable","1 or 2","source-defined")],
    "Generic n/m product family with drawn relative stereochemistry.", "specified", "not_applicable"),
])

R = bundle["reactions"]
GEN = f"{DOC}:generic_current_strategy"
R.append(reaction("generic_current_strategy", [E[1],E[2]], categories=["other"],
    location="Graphical abstract and Scheme 1, 'This work: Catalytic Asymmetric Synthesis'",
    name="Pd-catalyzed asymmetric hydroalkoxylation followed by RCM, cyclization, ozonolysis and reduction",
    transformation="generic multistep strategy",
    parts=[participant("p_ga","m_g_allene","reactant",[E[1],E[2]]), participant("p_gb","m_g_alcohol","reactant",[E[1],E[2]]), participant("p_gp","m_g_product","product",[E[1],E[2]],product_role="generic target")],
    steps=[], outcomes=[], changes=None, coverage="unresolved", record_type="generic_scheme",
    notes=("Graphical abstract reports yield >94% and ee >90% generically; these values are not assigned to a concrete experiment.",),
    unresolved=("R, X, n and m are generic variables; no single complete atom mapping applies.",)))

# Table 1: four independent hydroalkoxylation optimization experiments.
table = [
 (1,"m_l1","m_toluene","toluene","rt","5 min",89,62),
 (2,"m_l2","m_toluene","toluene","rt","5 min",94,88),
 (3,"m_l2","m_ch2cl2","CH2Cl2","rt","5 min",95,88),
 (4,"m_l2","m_ch2cl2","CH2Cl2","0 °C","2 h",94,96),
]
for n,lig,solmid,sol,temp,time,yld,ee in table:
    ev=[E[2],T[2]]
    parts=[]; add_standard_catalysis(parts,ev,"m1",lig,solvent=solmid)
    parts.extend([participant("p_product","m3","product",ev,product_role="adduct"),
                  participant("p_ee_derivative","m5","product",ev,stage="characterization",product_role="analytical derivative")])
    tout=(temp, temp if temp=="rt" else q(0,"°C"))
    ttime=(time, q(5,"min") if "min" in time else q(2,"h"))
    changes=[formed("C","O"),hgain("C"),order_down("C","C",2,1)]
    R.append(reaction(f"table1_entry_{n:02d}",ev,categories=["optimization"],
      location=f"Table 1, entry {n}",name="Pd-catalyzed asymmetric coupling / hydroalkoxylation",
      transformation="asymmetric hydroalkoxylation of alkoxyallene",
      parts=parts,steps=[step("reaction",ev,[p["participant_id"] for p in parts if p["stage"]=="reaction"],operation="combine under Table 1 conditions",temp=tout,time=ttime)],
      outcomes=[
       outcome("yield_3","yield",["m3"],ev,raw=f"{yld}%",value=q(yld,"%"),product_role="adduct 3",yield_type="isolated",basis="isolated yield",scope=f"Table 1 entry {n}"),
       outcome("ee_5","selectivity",["m5"],ev,raw=f"{ee}% ee",value=q(ee,"%"),product_role="analytical derivative 5",basis="HPLC measurement",scope="ee of compound 5; used by authors to assess the precursor sequence"),
      ],changes=changes,general=[GEN],solvent=sol,concentration=("0.5 M",q(0.5,"M")),
      unresolved=("Absolute configuration of the major enantiomer is not printed on the supplied pages.","Workup and purification are in unavailable SI.")))

# Table 1 supporting sequence and its explicitly reported negative controls.
ev2=[E[2],T[2]]
parts=[participant("p_3","m3","reactant",ev2),participant("p_g1","m_grubbs1","catalyst",ev2,mol=10,reference="m3"),participant("p_sol","m_ch2cl2","solvent",ev2),participant("p_4","m4","product",ev2,product_role="cyclic acetal")]
R.append(reaction("table1_rcm_3_to_4",ev2,categories=["other"],location="Table 1 reaction diagram, step ii; accompanying prose",
 name="ring-closing metathesis",transformation="ring-closing metathesis",parts=parts,
 steps=[step("reaction",ev2,[x["participant_id"] for x in parts],operation="RCM with first-generation Grubbs catalyst")],
 outcomes=[outcome("yield_4","yield",["m4"],ev2,raw="84% yield",value=q(84,"%"),product_role="cyclic acetal 4",scope="single RCM step")],
 changes=[formed("C","C",2)],general=[GEN],solvent="CH2Cl2",concentration=("0.1 M",q(0.1,"M")),
 unresolved=("Reaction time, temperature and yield measurement basis are not stated on the supplied pages.","Workup and purification are in unavailable SI.")))

parts=[participant("p_4","m4","reactant",ev2),participant("p_os","m_oso4","reagent",ev2,mol=3,reference="m4"),participant("p_bz","m_bzcl","reagent",ev2,eq=2,reference="m4"),participant("p_py","m_pyridine","reagent",ev2),participant("p_5","m5","product",ev2,product_role="bis-benzoate analytical derivative")]
R.append(reaction("table1_derivatization_4_to_5",ev2,categories=["derivatization"],location="Table 1 reaction diagram, steps iii–iv",
 name="dihydroxylation and benzoylation",transformation="two-step preparation of bis-benzoate analytical derivative",parts=parts,
 steps=[step("dihydroxylation",ev2,["p_4","p_os"],operation="OsO4 dihydroxylation"),step("benzoylation",ev2,["p_bz","p_py","p_5"],operation="BzCl/pyridine benzoylation")],
 outcomes=[outcome("overall_yield_5","yield",["m5"],ev2,raw="79% yield (over two steps)",value=q(79,"%"),product_role="bis-benzoate 5",yield_type="overall over two steps",scope="compound 4 to compound 5 across steps iii–iv")],
 changes=[formed("C","O"),formed("C","O"),formed("C","O"),formed("C","O"),hgain("C"),hgain("C"),order_down("C","C",2,1)],general=[GEN],
 unresolved=("Solvent, temperature, time, BzCl basis and yield measurement method are not stated on the supplied pages.","The figure gives pyridine without an amount.")))

for rid,catmid,label in [("rcm_control_grubbs2","m_grubbs2","Grubbs second"),("rcm_control_hoveyda2","m_hg2","Hoveyda Grubbs second")]:
    ev=[E[4],T[4]]
    parts=[participant("p_3","m3","reactant",ev),participant("p_cat",catmid,"catalyst",ev)]
    R.append(reaction(rid,ev,categories=["control"],location="Footnote 11 on original PDF page 4",name="RCM catalyst control",transformation="no conversion in attempted ring-closing metathesis",parts=parts,
      steps=[step("reaction",ev,["p_3","p_cat"],operation=f"reaction using {label}")],
      outcomes=[outcome("conversion_3","conversion",["m3"],ev,raw="showed no conversion",value="no conversion",product_role="starting material 3",scope="attempted RCM")],
      changes=[],coverage="not_applicable",general=[GEN],unresolved=("Catalyst loading, solvent, temperature, time and observation method are not reported in footnote 11.",)))

# Scheme 2 and the preceding thermal decomposition experiment.
parts=[participant("p_4","m4","reactant",ev2),participant("p_tin","m_bu3snh","reagent",ev2)]
R.append(reaction("thermal_bu3snh_decomposition_4",ev2,categories=["control"],location="Page 2 prose before Scheme 2",name="thermal Bu3SnH cyclization attempt",transformation="decomposition during attempted radical cyclization",parts=parts,
 steps=[step("reaction",ev2,["p_4","p_tin"],operation="employ Bu3SnH under thermal conditions",temp=("thermal conditions","thermal conditions"))],
 outcomes=[outcome("decomposition","decomposition",["m4"],ev2,raw="led to the extensive decomposition of 4",value="extensive decomposition",scope="attempted 5-exo cyclization")],
 changes=[],coverage="not_applicable",unresolved=("Exact temperature, time, solvent, loading and analytical method are unreported.",)))

def scheme2_records():
    ev=[E[2],T[2]]
    parts=[participant("p_4","m4","reactant",ev),participant("p_tin","m_bu3snh","reagent",ev,eq=1.1,reference="m4"),participant("p_o2","m_o2","reagent",ev),participant("p_be","m_bet3","reagent",ev,eq=1.0,reference="m4"),participant("p_tol","m_toluene","solvent",ev),participant("p_6","m6","product",ev,product_role="cyclized exo-olefin")]
    R.append(reaction("scheme2_4_to_6",ev,categories=["other"],location="Scheme 2, arrow 4→6",name="radical-mediated 5-exo cyclization",transformation="reductive radical cyclization",parts=parts,
      steps=[step("reaction",ev,[p["participant_id"] for p in parts],operation="Bu3SnH/O2/BEt3 radical cyclization",temp=("rt","rt"),time=("10 min",q(10,"min")))],
      outcomes=[outcome("yield_6","yield",["m6"],ev,raw="90%",value=q(90,"%"),product_role="compound 6")],changes=[formed("C","C"),hgain("C"),broken("C","Br")],general=[GEN],solvent="toluene",concentration=("0.1 M",q(0.1,"M")),unresolved=("Yield measurement basis and workup are unreported on the supplied pages.",)))
    parts=[participant("p_6","m6","reactant",ev),participant("p_o3","m_o3","reagent",ev),participant("p_dms","m_dms","reagent",ev),participant("p_sol","m_ch2cl2","solvent",ev),participant("p_7","m7","product",ev,product_role="bicyclic ketone")]
    R.append(reaction("scheme2_6_to_7",ev,categories=["other"],location="Scheme 2, arrow 6→7",name="ozonolysis",transformation="ozonolysis to ketone",parts=parts,
      steps=[step("ozonolysis",ev,["p_6","p_o3","p_sol"],operation="O3 in CH2Cl2",temp=("−78 °C",q(-78,"°C")),time=("30 min",q(30,"min"))),step("reductive_workup",ev,["p_dms","p_7"],operation="add (CH3)2S",temp=("−78 °C",q(-78,"°C")),time=("30 min",q(30,"min")))],
      outcomes=[outcome("yield_7","yield",["m7"],ev,raw="87%",value=q(87,"%"),product_role="bicyclic ketone 7")],changes=[formed("C","O",2),broken("C","C",2)],general=[GEN],solvent="CH2Cl2",concentration=("0.05 M",q(0.05,"M")),unresolved=("Yield measurement basis and workup beyond dimethyl sulfide are unreported.",)))
    parts=[participant("p_7","m7","reactant",ev),participant("p_bh","m_nabh4","reagent",ev,eq=1.2,reference="m7"),participant("p_sol","m_etoh","solvent",ev),participant("p_8","m8","product",ev,product_role="target alcohol")]
    R.append(reaction("scheme2_7_to_8",ev,categories=["other"],location="Scheme 2, arrow 7→8",name="NaBH4 reduction",transformation="ketone reduction to alcohol",parts=parts,
      steps=[step("reaction",ev,[p["participant_id"] for p in parts],operation="NaBH4 reduction in EtOH",temp=("0 °C",q(0,"°C")),time=("30 min",q(30,"min")))],
      outcomes=[outcome("yield_8","yield",["m8"],ev,raw="80%",value=q(80,"%"),product_role="alcohol 8"),outcome("diastereoselectivity_8","selectivity",["m8"],ev,raw="complete diastereoselectivity",value="complete diastereoselectivity",product_role="alcohol 8",scope="diastereoselectivity; major configuration shown, absolute R/S not printed")],
      changes=[order_down("C","O",2,1),hgain("C"),hgain("O")],general=[GEN],solvent="EtOH",concentration=("0.1 M",q(0.1,"M")),unresolved=("Yield and diastereoselectivity measurement methods are not reported on the supplied pages.",)))
    parts=[participant("p_8","m8","reactant",ev),participant("p_bz","m_bzcl","reagent",ev,eq=1.1,reference="m8"),participant("p_py","m_pyridine","reagent",ev,eq=1.1,reference="m8"),participant("p_sol","m_ch2cl2","solvent",ev),participant("p_9","m9","product",ev,product_role="benzoate analytical derivative")]
    R.append(reaction("scheme2_8_to_9",ev,categories=["derivatization"],location="Scheme 2, arrow 8→9",name="benzoylation",transformation="alcohol benzoylation for ee determination",parts=parts,
      steps=[step("reaction",ev,[p["participant_id"] for p in parts],operation="BzCl/pyridine benzoylation",temp=("rt","rt"),time=("7 h",q(7,"h")))],
      outcomes=[outcome("yield_9","yield",["m9"],ev,raw="97%",value=q(97,"%"),product_role="benzoate 9"),outcome("ee_9","selectivity",["m9"],ev,raw="96% ee",value=q(96,"%"),product_role="benzoate 9",scope="ee of measured benzoate 9; authors use it to assess alcohol 8")],
      changes=[formed("C","O"),broken("O","H")],general=[GEN],solvent="CH2Cl2",concentration=("0.1 M",q(0.1,"M")),unresolved=("Yield and ee measurement methods are not reported for 9 on the supplied pages; no absolute configuration is printed.",)))
scheme2_records()

# Scheme 3: six arrows retained independently; the prose's 72% combined yield is
# attached to the second arrow without relabeling it as a single-step yield.
ev=[E[2],T[2]]
parts=[]; add_standard_catalysis(parts,ev,"m10","m_l2",alcohol="m2",alcohol_eq=1.2,solvent="m_ch2cl2"); parts.append(participant("p_product","m11","product",ev,product_role="acetal 11"))
R.append(reaction("scheme3_10_2_to_11",ev,categories=["other"],location="Scheme 3, arrow 10+2→11",name="Pd-catalyzed asymmetric hydroalkoxylation",transformation="asymmetric hydroalkoxylation of alkoxyallene",parts=parts,
 steps=[step("reaction",ev,[p["participant_id"] for p in parts],operation="Pd/L2-catalyzed hydroalkoxylation",temp=("0 °C",q(0,"°C")),time=("30 min",q(30,"min")))],outcomes=[outcome("yield_11","yield",["m11"],ev,raw="91%",value=q(91,"%"),product_role="acetal 11")],changes=[formed("C","O"),hgain("C"),order_down("C","C",2,1)],general=[GEN],solvent="CH2Cl2",concentration=("0.5 M",q(0.5,"M")),unresolved=("Yield measurement basis and workup are unreported.",)))
parts=[participant("p_11","m11","reactant",ev),participant("p_g1","m_grubbs1","catalyst",ev,mol=10,reference="m11"),participant("p_sol","m_ch2cl2","solvent",ev),participant("p_12","m12","product",ev,product_role="cyclic acetal 12")]
R.append(reaction("scheme3_11_to_12",ev,categories=["other"],location="Scheme 3, arrow 11→12",name="ring-closing metathesis",transformation="ring-closing metathesis",parts=parts,steps=[step("reaction",ev,[p["participant_id"] for p in parts],operation="first-generation Grubbs RCM",temp=("rt","rt"),time=("4 h",q(4,"h")))],outcomes=[outcome("yield_12","yield",["m12"],ev,raw="99%",value=q(99,"%"),product_role="cyclic acetal 12")],changes=[formed("C","C",2)],general=[GEN],solvent="CH2Cl2",concentration=("0.2 M",q(0.2,"M")),unresolved=("Yield measurement basis and workup are unreported.",)))
parts=[participant("p_12","m12","reactant",ev),participant("p_tin","m_bu3snh","reagent",ev,eq=1.1,reference="m12"),participant("p_o2","m_o2","reagent",ev),participant("p_be","m_bet3","reagent",ev,eq=1.0,reference="m12"),participant("p_sol","m_toluene","solvent",ev),participant("p_13","m13","product",ev,product_role="cyclized exo-olefin 13")]
R.append(reaction("scheme3_12_to_13",ev,categories=["other"],location="Scheme 3, arrow 12→13",name="radical-mediated cyclization",transformation="reductive radical cyclization",parts=parts,steps=[step("reaction",ev,[p["participant_id"] for p in parts],operation="Bu3SnH/O2/BEt3 cyclization",temp=("rt","rt"),time=("1 h",q(1,"h")))],outcomes=[outcome("yield_13","yield",["m13"],ev,raw="90%",value=q(90,"%"),product_role="compound 13")],changes=[formed("C","C"),hgain("C"),broken("C","Br")],general=[GEN],solvent="toluene",concentration=("0.1 M",q(0.1,"M")),unresolved=("Yield measurement basis and workup are unreported.",)))
ev23=[E[2],T[2],E[3],T[3]]
parts=[participant("p_13","m13","reactant",ev23),participant("p_o3","m_o3","reagent",ev23),participant("p_dms","m_dms","reagent",ev23),participant("p_sol","m_ch2cl2","solvent",ev23),participant("p_14","m14","product",ev23,product_role="bicyclic ketone 14")]
R.append(reaction("scheme3_13_to_14",ev23,categories=["other"],location="Scheme 3, arrow 13→14; page 3 continuation prose",name="ozonolysis",transformation="ozonolysis to ketone",parts=parts,steps=[step("ozonolysis",ev23,["p_13","p_o3","p_sol"],operation="O3 in CH2Cl2",temp=("−78 °C",q(-78,"°C")),time=("1 h",q(1,"h"))),step("reductive_workup",ev23,["p_dms","p_14"],operation="add (CH3)2S",temp=("−78 °C",q(-78,"°C")),time=("30 min",q(30,"min")))],outcomes=[outcome("yield_14","yield",["m14"],ev23,raw="81%",value=q(81,"%"),product_role="ketone 14",scope="Scheme 3 single arrow 13→14"),outcome("combined_yield_12_to_14","yield",["m14"],ev23,raw="72% combined yield over two steps",value=q(72,"%"),product_role="ketone 14",yield_type="combined over two steps",scope="compound 12 through compound 13 to compound 14; not a relabeling of the ozonolysis yield")],changes=[formed("C","O",2),broken("C","C",2)],general=[GEN],solvent="CH2Cl2",concentration=("0.5 M",q(0.5,"M")),unresolved=("Measurement bases for both 81% and 72% are not reported on the supplied pages.",)))
parts=[participant("p_14","m14","reactant",ev23),participant("p_bh","m_nabh4","reagent",ev23,eq=1.5,reference="m14"),participant("p_sol","m_etoh","solvent",ev23),participant("p_15","m15","product",ev23,product_role="THF-THP alcohol 15")]
R.append(reaction("scheme3_14_to_15",ev23,categories=["other"],location="Scheme 3, arrow 14→15; page 3 prose",name="NaBH4 reduction",transformation="ketone reduction to alcohol",parts=parts,steps=[step("reaction",ev23,[p["participant_id"] for p in parts],operation="NaBH4 reduction in EtOH",temp=("0 °C",q(0,"°C")),time=("30 min",q(30,"min")))],outcomes=[outcome("yield_15","yield",["m15"],ev23,raw="71%",value=q(71,"%"),product_role="alcohol 15"),outcome("diastereoselectivity_15","selectivity",["m15"],ev23,raw="complete diastereoselectivity",value="complete diastereoselectivity",product_role="alcohol 15",scope="diastereoselectivity; source drawing specifies relative configuration")],changes=[order_down("C","O",2,1),hgain("C"),hgain("O")],general=[GEN],solvent="EtOH",concentration=("0.1 M",q(0.1,"M")),unresolved=("Yield and diastereoselectivity measurement methods are unreported.",)))
parts=[participant("p_15","m15","reactant",ev23),participant("p_bz","m_bzcl","reagent",ev23,eq=1.5,reference="m15"),participant("p_py","m_pyridine","reagent",ev23,eq=1.5,reference="m15"),participant("p_sol","m_ch2cl2","solvent",ev23),participant("p_16","m16","product",ev23,product_role="benzoate analytical derivative")]
R.append(reaction("scheme3_15_to_16",ev23,categories=["derivatization"],location="Scheme 3, arrow 15→16; page 3 prose",name="benzoylation",transformation="alcohol benzoylation for ee determination",parts=parts,steps=[step("reaction",ev23,[p["participant_id"] for p in parts],operation="BzCl/pyridine benzoylation",temp=("rt","rt"),time=("12 h",q(12,"h")))],outcomes=[outcome("yield_16","yield",["m16"],ev23,raw="83%",value=q(83,"%"),product_role="benzoate 16"),outcome("ee_16","selectivity",["m16"],ev23,raw="90% ee",value=q(90,"%"),product_role="benzoate 16",scope="ee of measured benzoate 16; authors use it to assess alcohol 15")],changes=[formed("C","O"),broken("O","H")],general=[GEN],solvent="CH2Cl2",concentration=("0.1 M",q(0.1,"M")),unresolved=("Yield and ee measurement methods are not reported for 16; no absolute configuration is printed.",)))

# Scheme 4: seven arrows, including the failed 6-exo experiment.
ev3=[E[3],T[3]]
parts=[]; add_standard_catalysis(parts,ev3,"m1","m_l2",alcohol="m17",alcohol_eq=None,solvent="m_ch2cl2"); parts.append(participant("p_product","m18","product",ev3,product_role="acetal 18"))
R.append(reaction("scheme4_1_17_to_18",ev3,categories=["other"],location="Scheme 4, arrow 1+17→18",name="Pd-catalyzed asymmetric hydroalkoxylation",transformation="asymmetric hydroalkoxylation of alkoxyallene",parts=parts,steps=[step("reaction",ev3,[p["participant_id"] for p in parts],operation="Pd/L2-catalyzed hydroalkoxylation",temp=("0 °C",q(0,"°C")),time=("2 h",q(2,"h")))],outcomes=[outcome("yield_18","yield",["m18"],ev3,raw="99%",value=q(99,"%"),product_role="acetal 18")],changes=[formed("C","O"),hgain("C"),order_down("C","C",2,1)],general=[GEN],solvent="CH2Cl2",concentration=("0.5 M",q(0.5,"M")),unresolved=("The scheme does not state an equivalent for alcohol 17; yield measurement basis and workup are unreported.",)))
parts=[participant("p_18","m18","reactant",ev3),participant("p_g1","m_grubbs1","catalyst",ev3,mol=5,reference="m18"),participant("p_sol","m_ch2cl2","solvent",ev3),participant("p_19","m19","product",ev3,product_role="cyclic acetal 19")]
R.append(reaction("scheme4_18_to_19",ev3,categories=["other"],location="Scheme 4, arrow 18→19",name="ring-closing metathesis",transformation="ring-closing metathesis",parts=parts,steps=[step("reaction",ev3,[p["participant_id"] for p in parts],operation="first-generation Grubbs RCM",temp=("rt","rt"),time=("2 h",q(2,"h")))],outcomes=[outcome("yield_19","yield",["m19"],ev3,raw="98%",value=q(98,"%"),product_role="cyclic acetal 19")],changes=[formed("C","C",2)],general=[GEN],solvent="CH2Cl2",concentration=("0.1 M",q(0.1,"M")),unresolved=("Yield measurement basis and workup are unreported.",)))
parts=[participant("p_19","m19","reactant",ev3),participant("p_tin","m_bu3snh","reagent",ev3,eq=1.1,reference="m19"),participant("p_o2","m_o2","reagent",ev3),participant("p_be","m_bet3","reagent",ev3,eq=1.0,reference="m19"),participant("p_21","m21","product",ev3,product_role="intended 6-exo product"),participant("p_side","m_dehal19","product",ev3,product_role="minor dehalogenated side product")]
R.append(reaction("scheme4_failed_19_to_21",ev3,categories=["control"],location="Scheme 4, crossed-out arrow 19→21; page 3 prose and footnote 16",name="attempted 6-exo radical cyclization",transformation="failed 6-exo radical cyclization with minor dehalogenation",parts=parts,steps=[step("reaction",ev3,[p["participant_id"] for p in parts],operation="Bu3SnH/O2/BEt3 attempted 6-exo-trig cyclization",temp=("rt","rt"),time=("10 min",q(10,"min")))],outcomes=[outcome("no_21","no_reaction",[],ev3,raw="failed to produce product 21",value=True,scope="intended 6-exo cyclization"),outcome("minor_dehalogenation","trace",["m_dehal19"],ev3,raw="a small amount of the dehalogenated product was obtained",value="small amount",product_role="dehalogenated side product",scope="observed side reaction")],changes=[hgain("C"),broken("C","Br")],general=[GEN],unresolved=("No solvent, analytical method, identity, amount or structure is reported for the minor dehalogenated product.","The intended product is shown only as the failed target and is not reported as obtained."),limitation="Classification covers only the explicitly observed dehalogenation side reaction, not the failed intended cyclization."))
parts=[participant("p_19","m19","reactant",ev3),participant("p_pd","m_pd2dba3","catalyst",ev3,mol=5,reference="m19"),participant("p_dppp","m_dppp","ligand",ev3,mol=10,reference="m19"),participant("p_base","m_et3n","additive",ev3,eq=2,reference="m19"),participant("p_sol","m_dmf","solvent",ev3),participant("p_20","m20","product",ev3,product_role="diene 20")]
R.append(reaction("scheme4_19_to_20",ev3,categories=["other"],location="Scheme 4, arrow 19→20",name="Pd-catalyzed cyclization",transformation="Pd-catalyzed ring-forming cyclization",parts=parts,steps=[step("reaction",ev3,[p["participant_id"] for p in parts],operation="Pd2(dba)3/dppp cyclization",temp=("80 °C",q(80,"°C")),time=("12 h",q(12,"h")))],outcomes=[outcome("yield_20","yield",["m20"],ev3,raw="78%",value=q(78,"%"),product_role="diene 20")],changes=[formed("C","C"),broken("C","Br")],general=[GEN],solvent="DMF",concentration=("0.1 M",q(0.1,"M")),unresolved=("The dense product 20 was not safely encoded as SMILES; yield measurement basis and workup are unreported.",)))
parts=[participant("p_20","m20","reactant",ev3),participant("p_pdc","m_pdc","catalyst",ev3),participant("p_h2","m_h2","reagent",ev3),participant("p_bipy","m_bipy","additive",ev3,eq=0.25,reference="m20"),participant("p_sol","m_meoh","solvent",ev3),participant("p_21","m21","product",ev3,product_role="exo-olefin 21")]
R.append(reaction("scheme4_20_to_21",ev3,categories=["other"],location="Scheme 4, arrow 20→21",name="chemoselective hydrogenation",transformation="selective alkene hydrogenation",parts=parts,steps=[step("reaction",ev3,[p["participant_id"] for p in parts],operation="Pd/C hydrogenation with 2,2'-bipyridyl",temp=("rt","rt"),time=("12 h",q(12,"h")))],outcomes=[outcome("yield_21","yield",["m21"],ev3,raw="80%",value=q(80,"%"),product_role="exo-olefin 21")],changes=[order_down("C","C",2,1),hgain("C"),hgain("C")],general=[GEN],solvent="MeOH",concentration=("0.1 M",q(0.1,"M")),unresolved=("Pressure is printed as H2 (1 atm) but the generic pressure field remains unstructured because the schema quantity units do not include atm; it is retained in notes.","Yield measurement basis and workup are unreported."),notes=("H2 (1 atm)",)))
parts=[participant("p_21","m21","reactant",ev3),participant("p_o3","m_o3","reagent",ev3),participant("p_dms","m_dms","reagent",ev3),participant("p_bh","m_nabh4","reagent",ev3,eq=1.2,reference="m21"),participant("p_dc","m_ch2cl2","solvent",ev3),participant("p_et","m_etoh","solvent",ev3),participant("p_22","m22","product",ev3,product_role="THP-THF alcohol 22")]
R.append(reaction("scheme4_21_to_22",ev3,categories=["other"],location="Scheme 4, arrow 21→22, two-step sequence",name="ozonolysis and NaBH4 reduction",transformation="two-step conversion of exo-olefin to alcohol",parts=parts,steps=[step("ozonolysis",ev3,["p_21","p_o3","p_dc"],operation="O3 in CH2Cl2",temp=("−78 °C",q(-78,"°C")),time=("30 min",q(30,"min"))),step("ozonolysis_workup",ev3,["p_dms"],operation="then (CH3)2S",temp=("−78 °C",q(-78,"°C"))),step("reduction",ev3,["p_bh","p_et","p_22"],operation="NaBH4 reduction in EtOH",temp=("0 °C",q(0,"°C")),time=("30 min",q(30,"min")))],outcomes=[outcome("overall_yield_22","yield",["m22"],ev3,raw="70% (over two steps)",value=q(70,"%"),product_role="alcohol 22",yield_type="overall over two steps",scope="compound 21 to compound 22 across ozonolysis and reduction")],changes=[formed("C","O"),hgain("C"),hgain("O"),broken("C","C",2)],general=[GEN],unresolved=("The time after dimethyl sulfide addition is not printed; the two solvents belong to separate stages.","Yield measurement basis and workup are unreported.")))
parts=[participant("p_22","m22","reactant",ev3),participant("p_bz","m_bzcl","reagent",ev3,eq=1.2,reference="m22"),participant("p_py","m_pyridine","reagent",ev3,eq=1.2,reference="m22"),participant("p_sol","m_ch2cl2","solvent",ev3),participant("p_23","m23","product",ev3,product_role="benzoate analytical derivative")]
R.append(reaction("scheme4_22_to_23",ev3,categories=["derivatization"],location="Scheme 4, arrow 22→23",name="benzoylation",transformation="alcohol benzoylation for ee determination",parts=parts,steps=[step("reaction",ev3,[p["participant_id"] for p in parts],operation="BzCl/pyridine benzoylation",temp=("rt","rt"),time=("5 h",q(5,"h")))],outcomes=[outcome("yield_23","yield",["m23"],ev3,raw="90%",value=q(90,"%"),product_role="benzoate 23"),outcome("ee_23","selectivity",["m23"],ev3,raw="95% ee",value=q(95,"%"),product_role="benzoate 23",scope="ee of measured benzoate 23; authors use it to assess alcohol 22")],changes=[formed("C","O"),broken("O","H")],general=[GEN],solvent="CH2Cl2",concentration=("0.1 M",q(0.1,"M")),unresolved=("Yield and ee measurement methods are not reported for 23; no absolute configuration is printed.",)))

bundle["source_context"] = [
 {"context_id":f"{DOC}:figure1_background","document_id":DOC,"kind":"background_bioactive_structures","evidence_ids":[E[1],T[1]],"note":"Figure 1 shows Darunavir and related inhibitors as background, not current-paper synthesis experiments."},
 {"context_id":f"{DOC}:scheme1_prior_work","document_id":DOC,"kind":"literature_background","evidence_ids":[E[2]],"note":"Scheme 1 chiral-pool and enzymatic-resolution routes summarize prior approaches; they are not attributed as current-paper experiments."},
 {"context_id":f"{DOC}:unspecified_failed_variants","document_id":DOC,"kind":"unresolved_experimental_discussion","evidence_ids":[T[3],E[3]],"note":"The prose says varying initiators and hydrogen donors did little to improve the failed 19 cyclization, but identities and conditions are not enumerable; no fabricated rows were created."},
 {"context_id":f"{DOC}:si_notice","document_id":DOC,"kind":"supporting_information_notice","evidence_ids":[E[3],T[3]],"note":"The article announces SI and FAIR NMR data, but no SI was supplied or explicitly associated in this workspace."},
 {"context_id":f"{DOC}:bibliography_page4","document_id":DOC,"kind":"bibliography","evidence_ids":[E[4],T[4]],"note":"Page 4 is bibliography/footnotes. Footnote 11 contains two current-paper no-conversion controls recorded separately; the remaining citations are not current experiments."},
]

# Preserve source-mention hints and graph diagnostics as raw machine artifacts.
mentions = st.source_mentions(data, ["compound 4", "compound 5", "compound 8", "compound 15", "compound 22", "Grubbs second", "no conversion"], context=220)
st.write(RUN / "raw_model/source_mentions.json", mentions)
diagnostics = {}
for m in bundle["molecules"]:
    smi = st.value(m["smiles"])
    if smi:
        diagnostics[m["molecule_id"]] = st.stereo_inspect(smi)
st.write(RUN / "raw_model/stereo_diagnostics.json", diagnostics)
inventory = {
 "document_id": DOC,
 "pages": [1,2,3,4],
 "full_page_images_opened": [E[i] for i in range(1,5)],
 "formal_reaction_records": len(R),
 "molecules": len(bundle["molecules"]),
 "record_ids": [r["reaction_id"] for r in R],
 "notes": ["Raw inventory is machine-authored and remains needs_review.","Crops and their evidence receipts are retained beside this file."],
}
st.write(RUN / "raw_model/inventory.json", inventory)
st.write(RUN / "submitted.json", bundle)
print(json.dumps({"molecules":len(bundle["molecules"]),"reactions":len(R)},ensure_ascii=False))
