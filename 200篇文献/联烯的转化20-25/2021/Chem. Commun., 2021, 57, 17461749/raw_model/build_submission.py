#!/usr/bin/env python3
"""Build the source-bounded v3.1 bundle for doc_1a0071146ff7cf43.

This file is deliberately source-specific.  It does not read old experiment
records and does not promote machine output to source-reviewed status.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DATA = json.loads((RUN / "dataset.json").read_text())
DOC = "doc_1a0071146ff7cf43"
E = {p: f"{DOC}:page{p}:image" for p in range(1, 5)}
T = {p: f"{DOC}:page{p}:text" for p in range(1, 5)}
INIT = RUN / "raw_model/standard_init_final.json"
B = json.loads(INIT.read_text()) if INIT.exists() else st.empty(DATA["documents"])


def f(raw=None, norm=None, ev=(), method="original_extraction", status=None, derivation=None):
    return st.fact(raw, norm, ev, method=method, status=status, derivation=derivation)


def inferred(norm, ev, why):
    return f(None, norm, ev, method="inference", derivation=why)


def add_molecule(mid, label, ev, name=None, smiles=None, chirality="unspecified",
                 geometry="unspecified", representation=None, atom_centers=None,
                 double_bonds=None, axial=None, mixture=None, assignment=None,
                 note=None):
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    m["names"]["original"] = f(name or label, name or label, ev)
    m["names"]["english"] = f(name or label, name or label, ev)
    m["names"]["abbreviations"] = f()
    m["names"]["aliases"] = f()
    m["cas_numbers"] = f()
    m["identity_note"] = f(note, note, ev) if note else f()
    m["representation_type"] = "smiles" if smiles else "unresolved"
    m["smiles"] = f(smiles, smiles, ev, method="structure_derivation",
                     derivation="Candidate graph transcribed from the source drawing/name; parsing does not verify identity.") if smiles else f()
    m["generic_template"] = None
    s = m["stereochemistry"]
    s["atom_centers"] = f(atom_centers, atom_centers, ev) if atom_centers is not None else f()
    s["double_bonds"] = f(double_bonds, double_bonds, ev) if double_bonds is not None else f()
    s["axial"] = f(axial, axial, ev) if axial is not None else f()
    s["mixture"] = f(mixture, mixture, ev) if mixture is not None else f()
    s["assignment"] = f(assignment, assignment, ev) if assignment is not None else f()
    s["chirality_status"] = inferred(chirality, ev, "Stereo status assigned conservatively from the displayed structure and source labels; missing marks are not treated as achiral or racemic.")
    s["geometry_status"] = inferred(geometry, ev, "Geometry status assigned conservatively from the displayed structure and source labels.")
    rep = representation or ("partial" if smiles else "not_encoded")
    s["representation_status"] = inferred(rep, ev, "Candidate SMILES, when present, does not encode unreported allene/center/relative stereochemistry.")
    if smiles:
        d = st.stereo_inspect(smiles)
        s["isomeric_smiles"] = f(None, d["canonical_isomeric_smiles"], ev, method="structure_derivation",
                                 derivation="RDKit canonicalization of the unverified candidate graph.")
        s["identity_key"] = f(None, d["stereo_graph_key"], ev, method="structure_derivation",
                              derivation="RDKit stereo graph key for the unverified candidate graph; not a sample identity.")
    else:
        s["isomeric_smiles"] = f()
        s["identity_key"] = f()
    B["molecules"].append(m)
    return mid


def add_template(mid, label, ev, fixed, variables, reason):
    m = add_molecule(mid, label, ev, name=label, chirality="unresolved", geometry="unresolved", representation="not_encoded")
    obj = B["molecules"][-1]
    obj["representation_type"] = "generic_template"
    obj["smiles"] = f()
    obj["generic_template"] = {
        "fixed_structure": f(fixed, fixed, ev),
        "variables": [{
            "label": v,
            "attachment_points": f("shown in source drawing", "shown in source drawing", ev),
            "definition": f("source-defined substituent", "source-defined substituent", ev),
            "allowed_values": f(),
            "constraints": f(),
        } for v in variables],
        "instance_mapping": [],
        "reason": f(reason, reason, ev),
    }
    return m


def participant(pid, mid, role, ev, eq=None, molpct=None, amount=None, stage="reaction", product_role=None, ref=None):
    p = st.new_participant(pid, mid, role, ev, stage=stage)
    if product_role:
        p["product_role"] = f(product_role, product_role, ev)
    if eq is not None:
        p["charges"]["equivalents"] = f(f"{eq:g} equiv", {"number": eq, "unit": "equiv"}, ev)
    if molpct is not None:
        p["charges"]["mol_percent"] = f(f"{molpct:g} mol%", {"number": molpct, "unit": "mol%"}, ev)
    if amount is not None:
        p["charges"]["amount"] = f(f"{amount:g} mmol", {"number": amount, "unit": "mmol"}, ev)
    if eq is not None or molpct is not None:
        p["equivalence_reference"] = f(ref, ref, ev)
    return p


def step(sid, ev, pids, temp=None, time=None, atmosphere=None, operation="reaction mixture held under the reported conditions"):
    s = st.new_step(sid, ev)
    s["participant_ids"] = list(pids)
    s["operation"] = f(operation, operation, ev)
    if temp is not None:
        s["temperature"] = f(temp, temp if isinstance(temp, str) else temp, ev)
    if time is not None:
        if isinstance(time, dict):
            q = time.get("raw_text") or f"{time.get('qualifier', '')}{time.get('number', '')} {time.get('unit', '')}".strip()
        else:
            q = time
        s["time"] = f(q, time, ev)
    if atmosphere is not None:
        s["atmosphere"] = f(atmosphere, atmosphere, ev)
    return s


def conflict_temperature(ev):
    a = f("0–5 °C", {"minimum": 0, "maximum": 5, "unit": "°C", "raw_text": "0–5 °C"}, ev)
    b = f("rt", "r.t.", ev)
    return f(["0–5 °C", "rt"], None, ev, status="conflict", derivation="Scheme 3 header prints both temperature expressions; neither is selected.") | {"alternatives": [a, b]}


def outcome(oid, kind, targets, ev, raw, norm, basis=None, scope=None, ytype=None):
    o = st.new_outcome(oid, kind, targets, ev)
    o["value"] = f(raw, norm, ev)
    if basis is not None:
        o["measurement_basis"] = f(basis, basis, ev)
    if scope is not None:
        o["scope"] = f(scope, scope, ev)
    if ytype is not None:
        o["yield_type"] = f(ytype, ytype, ev)
    if kind in ("yield", "product_amount"):
        o["product_role"] = f("reported product", "reported product", ev)
    return o


def classify(r, ev, codes=None, coverage="partial", limitation=None):
    c = r["classification"]
    c["mapped_reaction_smiles"] = f()
    if not codes:
        c["bond_changes"] = []
        c["labels"] = inferred([], ev, "No observed bond-forming product is assigned for this record.") if coverage == "not_applicable" else f()
        c["coverage"] = inferred(coverage, ev, "No complete atom mapping was asserted.")
        c["limitations"] = [limitation or "No complete source-verified atom mapping is available."]
        c["model_votes"] = []
        return
    label_meta = {row["code"]: row["number"] for row in json.loads((ROOT / "skills/chemical-reaction-extraction/references/bond_taxonomy.json").read_text())["categories"]}
    pairs = {"C-C": ("C", "C"), "C-O": ("C", "O"), "O-H": ("O", "H")}
    changes = []
    labels = []
    for i, code in enumerate(codes, 1):
        a, b = pairs[code]
        cid = f"b{i}"
        ch = {"change_id": cid, "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
              "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain" if code == "O-H" else "formed"}
        if code == "O-H":
            ch["hydrogen_count"] = 1
        changes.append(f(None, ch, ev, method="inference", derivation="Class-level net bond change inferred from the source transformation drawing/mechanistic description; atom identities are unresolved."))
        labels.append({"category_number": label_meta[code], "code": code, "supporting_change_ids": [cid]})
    c["bond_changes"] = changes
    c["labels"] = f(None, labels, ev, method="inference", derivation="Taxonomy labels follow the partial class-level bond changes, not a complete atom mapping.")
    c["coverage"] = inferred("partial", ev, "Class-level changes only; exact atoms, multiplicities, proton sources, and coproducts are unresolved.")
    c["limitations"] = [limitation or "Exact atom identities, full atom mapping, hydrogen sources, coproducts, and all bond-order changes remain unresolved."]
    c["model_votes"] = []


def reaction(suffix, ev, record_type="specific_experiment", categories=None, name=None, transform=None):
    rid = f"{DOC}:{suffix}"
    r = st.new_reaction(rid, DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories or [], categories or [], ev)
    r["location"] = f({"pages": sorted({int(x.split(':page')[1].split(':')[0]) for x in ev}), "source_region": suffix},
                      {"pages": sorted({int(x.split(':page')[1].split(':')[0]) for x in ev}), "source_region": suffix}, ev)
    if name:
        r["reaction_name"] = f(name, name, ev)
    if transform:
        r["transformation_type"] = f(transform, transform, ev)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
                   "unresolved": ["No associated SI was supplied; machine record requires expert source review."]}
    r["completeness"] = "partial"
    B["reactions"].append(r)
    return r


# Core reagents/solvents and source-identified reactant families.
simple = [
    ("m_dbu", "DBU", "1,8-diazabicyclo[5.4.0]undec-7-ene", "C1CCC2=NCCCN2CC1"),
    ("m_dabco", "DABCO", "1,4-diazabicyclo[2.2.2]octane", "C1CN2CCN1CC2"),
    ("m_dmap", "DMAP", "4-dimethylaminopyridine", "CN(C)c1ccncc1"),
    ("m_et3n", "Et3N", "triethylamine", "CCN(CC)CC"),
    ("m_pph3", "PPh3", "triphenylphosphine", "P(c1ccccc1)(c1ccccc1)c1ccccc1"),
    ("m_mecn", "acetonitrile", "acetonitrile", "CC#N"),
    ("m_toluene", "toluene", "toluene", "Cc1ccccc1"),
    ("m_dioxane", "1,4-dioxane", "1,4-dioxane", "O1CCOCC1"),
    ("m_dcm", "DCM", "dichloromethane", "ClCCl"),
    ("m_thf", "THF", "tetrahydrofuran", "C1CCOC1"),
    ("m_etoh", "EtOH", "ethanol", "CCO"),
    ("m_benzonitrile", "benzonitrile", "benzonitrile", "N#Cc1ccccc1"),
    ("m_dmso", "DMSO", "dimethyl sulfoxide", "CS(C)=O"),
]
for mid, lab, name, smi in simple:
    if mid == "m_dabco":
        add_molecule(mid, lab, [E[2], E[3]], name=name, smiles=smi, chirality="unspecified", geometry="not_applicable", representation="partial",
                     note="Candidate graph diagnostics expose unassigned potential tetrahedral stereo; the source makes no stereochemical assignment.")
    else:
        add_molecule(mid, lab, [E[2], E[3]], name=name, smiles=smi, chirality="achiral", geometry="not_applicable", representation="encoded")

isatin_defs = {
    "NMe": ("N-methyl isatin", "CN1C(=O)C(=O)c2ccccc21"),
    "NPr": ("N-propyl isatin", "CCCN1C(=O)C(=O)c2ccccc21"),
    "NBn": ("N-benzyl isatin", "O=C1C(=O)N(Cc2ccccc2)c2ccccc21"),
    "NPropargyl": ("N-propargyl isatin", "C#CCN1C(=O)C(=O)c2ccccc21"),
    "NMe_5Cl": ("N-methyl-5-chloroisatin", "CN1C(=O)C(=O)c2cc(Cl)ccc21"),
    "NPr_5Cl": ("N-propyl-5-chloroisatin", "CCCN1C(=O)C(=O)c2cc(Cl)ccc21"),
    "NBn_5Cl": ("N-benzyl-5-chloroisatin", "O=C1C(=O)N(Cc2ccccc2)c2cc(Cl)ccc21"),
    "NMe_5Br": ("N-methyl-5-bromoisatin", "CN1C(=O)C(=O)c2cc(Br)ccc21"),
    "NPr_5Br": ("N-propyl-5-bromoisatin", "CCCN1C(=O)C(=O)c2cc(Br)ccc21"),
    "NBn_5Br": ("N-benzyl-5-bromoisatin", "O=C1C(=O)N(Cc2ccccc2)c2cc(Br)ccc21"),
    "NMe_5Me": ("N-methyl methyl-substituted isatin (Scheme 2/3)", None),
    "NMe_dimethyl": ("N-methyl dimethyl-substituted isatin (Scheme 3)", None),
    "NBn_5Me": ("N-benzyl methyl-substituted isatin (Scheme 2)", None),
    "NMe_5OMe": ("N-methyl-5-methoxyisatin", "COc1ccc2C(=O)C(=O)N(C)c2c1"),
    "NMe_5NO2": ("N-methyl-5-nitroisatin", "CN1C(=O)C(=O)c2cc([N+](=O)[O-])ccc21"),
}
isatin_ids = {}
for key, (name, smi) in isatin_defs.items():
    isatin_ids[key] = add_molecule(f"m_isatin_{key}", name, [E[2], E[3]], name=name, smiles=smi,
                                   chirality="achiral" if smi else "unspecified", geometry="not_applicable" if smi else "unspecified",
                                   note="Ring substitution position/identity follows the compact source drawing; uncertain variants are deliberately unencoded.")


def allenoate_smiles(r2, r3):
    ester = {"Me": "OC", "Et": "OCC", "tBu": "OC(C)(C)C"}[r2]
    left = {"H": "C=", "Me": "CC=", "Et": "CCC=", "Bn": "c1ccccc1CC=", "Ph": "c1ccccc1C=", "tBu": "CC(C)(C)C="}[r3]
    return f"{left}C=CC(=O){ester}"


allenoate_ids = {}
for r2, r3 in sorted({("Et", "Me"), ("Me", "Me"), ("tBu", "Me"), ("Et", "Bn"), ("Et", "H"),
                      ("Me", "Et"), ("Me", "tBu"), ("tBu", "tBu"), ("Me", "Ph"), ("Et", "Ph")}):
    mid = f"m_allenoate_{r2}_{r3}"
    axial = [{"fragment": "C=C=C", "configuration": "unassigned", "basis": "no axial descriptor in source"}]
    ch = "unspecified" if r3 != "H" else "achiral"
    allenoate_ids[(r2, r3)] = add_molecule(mid, f"allenoate 2 (R2={r2}, R3={r3})", [E[2], E[3]],
                                           name=f"{r2} 3-{r3}-allenoate (source substituent notation)",
                                           smiles=allenoate_smiles(r2, r3), chirality=ch, geometry="unspecified",
                                           axial=axial, representation="partial",
                                           note="Candidate connectivity is derived from the displayed allenoate family; allene-axis configuration is not reported.")

# Generic R-group families required for Scheme 1/current-work records.
add_template("m_generic_isatin", "isatin family 1", [E[2]], "isatin core with N-R and optional ring R1", ["R", "R1"], "Source shows an R-group family, not one concrete molecule.")
add_template("m_generic_allenoate", "allenoate family 2", [E[2]], "3-substituted allenoate with ester R2 and gamma R3", ["R2", "R3"], "Source shows an R-group family, not one concrete molecule.")
add_template("m_generic_spiro", "spirofuran oxindole family 3", [E[2]], "spirofuran bis-oxindole framework", ["R", "R1", "R2", "R3"], "Current-work generic product family.")
add_template("m_generic_gamma", "gamma-functionalized allenoate family 4", [E[2], E[3]], "isatin-derived gamma-functionalized allenoate", ["R", "R1", "R2", "R3"], "Current-work generic product family.")
for key, lab in [("jia", "Jia et al. 2011 product family"), ("nair", "Nair et al. 2015 product family"), ("yao", "Yao et al. 2019 product family")]:
    add_template(f"m_lit_{key}", lab, [E[2]], "source-drawn cited product scaffold", ["R", "R1", "R2", "R3"], "Cited prior-work R-group drawing.")

# Product registries.
scope3 = [
    ("3a", "NMe", "Et", "Me", 75), ("3b", "NMe", "Me", "Me", 86), ("3c", "NMe", "tBu", "Me", 53),
    ("3d", "NMe", "Et", "Bn", 58), ("3e", "NMe", "Et", "H", 31), ("3f", "NPr", "Et", "Me", 80),
    ("3g", "NPr", "tBu", "Me", 73), ("3h", "NBn", "Et", "Me", 65), ("3i", "NBn_5Me", "Et", "Me", 70),
    ("3j", "NBn", "tBu", "Me", 60), ("3k", "NBn_5Me", "tBu", "Me", 63),
    ("3l", "NMe_5Cl", "Et", "Me", 70), ("3m", "NMe_5Cl", "tBu", "Me", 71),
    ("3n", "NPr_5Cl", "Et", "Me", 73), ("3o", "NPr_5Cl", "tBu", "Me", 67), ("3p", "NPr_5Cl", "Me", "Me", 80),
    ("3q", "NBn_5Cl", "Et", "Me", 60), ("3r", "NBn_5Cl", "tBu", "Me", 51),
    ("3s", "NMe_5Br", "Et", "Me", 62), ("3t", "NMe_5Br", "tBu", "Me", 56), ("3u", "NMe_5Br", "Et", "Bn", 45),
    ("3v", "NPr_5Br", "Et", "Me", 69), ("3w", "NPr_5Br", "tBu", "Me", 62),
    ("3x", "NBn_5Br", "Et", "Me", 48), ("3y", "NBn_5Br", "tBu", "Me", 40),
    ("3z", "NMe_5Me", "Et", "Me", 81), ("3aa", "NMe_5Me", "tBu", "tBu", 52),
    ("3ab", "NMe_5Me", "Et", "Bn", 61), ("3ac", "NMe_5Me", "Et", "Me", 64),
    ("3ad", "NMe_5Me", "tBu", "Me", 76), ("3ae", "NMe_5Me", "Me", "tBu", 57),
    ("3af", "NMe_5Me", "Me", "Me", 88), ("3ag", "NMe", "tBu", "Me", 71),
    ("3ah", "NPropargyl", "Et", "Me", 72), ("3ai", "NPropargyl", "tBu", "Me", 61),
    ("3aj", "NPropargyl", "Me", "Me", 81),
]
prod3 = {}
for lab, ikey, r2, r3, y in scope3:
    prod3[lab] = add_molecule(f"m_product_{lab}", lab, [E[2]], name=f"spirofuran oxindole {lab}",
                              chirality="partially_specified", geometry="unspecified", representation="not_encoded",
                              atom_centers=[{"status": "relative stereochemistry shown by wedge/dash in Scheme 2; absolute configuration unassigned"}],
                              assignment="Source drawing specifies relative stereochemistry; 3a is described as a single diastereomer and 3ag was confirmed by X-ray, but no absolute configuration is assigned.",
                              note=f"Scheme 2 caption: {lab}; isatin={ikey}, R2={r2}, R3={r3}, yield={y}%.")

scope4 = [
    ("4a", "NMe", "Et", "Me", 87, "1:0.68"), ("4b", "NMe", "Et", "H", 85, "1:0.82"),
    ("4c", "NMe", "Et", "Bn", 70, "1:0.81"), ("4d", "NMe", "Et", "Ph", 68, "1:0.10"),
    ("4e", "NMe", "Me", "Ph", 77, "1:0.43"), ("4f", "NPr", "Et", "Me", 81, "1:0.94"),
    ("4g", "NBn", "Et", "Me", 76, "1:0.91"), ("4h", "NMe_5Cl", "Et", "Me", 79, "1:0.86"),
    ("4i", "NPr_5Cl", "Et", "Me", 72, "1:0.88"), ("4j", "NBn_5Cl", "Et", "Me", 65, "1:0.62"),
    ("4k", "NMe_5Br", "Et", "Me", 72, "1:0.95"), ("4l", "NMe_5Br", "Me", "Ph", 54, "1:0.54"),
    ("4m", "NPr_5Br", "Et", "Me", 67, "1:0.84"), ("4n", "NBn_5Br", "Et", "Me", 59, "1:0.79"),
    ("4o", "NMe_5Me", "Et", "Me", 79, "1:0.88"), ("4p", "NMe_5Me", "Et", "Ph", 65, "1:0.43"),
    ("4q", "NMe_5Me", "Me", "Ph", 71, "1:0.42"), ("4r", "NMe_dimethyl", "Et", "Me", 68, "1:0.95"),
    ("4s", "NPropargyl", "Et", "Me", 54, "1:0.83"),
]
prod4 = {}
for lab, ikey, r2, r3, y, dr in scope4:
    prod4[lab] = add_molecule(f"m_product_{lab}", lab, [E[3]], name=f"gamma-functionalized allenoate oxindole {lab}",
                              chirality="mixture", geometry="unspecified", representation="not_encoded",
                              atom_centers=[{"status": "stereogenic center(s) shown; absolute configurations unassigned"}],
                              axial=[{"fragment": "retained allene", "configuration": "unassigned"}],
                              mixture="Diastereomeric mixture; experiment-specific ratio is stored as a selectivity outcome.",
                              assignment="No absolute stereochemical assignment; 4h connectivity/relative structure was confirmed by X-ray.",
                              note=f"Scheme 3 caption: {lab}; isatin={ikey}, R2={r2}, R3={r3}, yield={y}%, dr={dr}.")

# Scheme 1 cited families and current generic reaction.
for suffix, product, rawname, detail in [
    ("literature_jia_2011", "m_lit_jia", "Jia et al. 2011 cited isatin/allenoate/isocyanide cycloaddition", "PhMe, R4-NC, 100 °C, 18–72 h"),
    ("literature_nair_2015", "m_lit_nair", "Nair et al. 2015 cited phosphine-mediated spiro tetrahydrofuran formation", "PPh3, MeCN, rt, 0.5 h"),
    ("literature_yao_2019", "m_lit_yao", "Yao et al. 2019 cited NHC-catalyzed gamma-functionalization", "NHC, MeCN, 35 °C"),
]:
    r = reaction(suffix, [E[2]], "literature_example", [], rawname, "cited prior-work transformation")
    r["participants"] = [participant("p_isatin", "m_generic_isatin", "reactant", [E[2]]), participant("p_allenoate", "m_generic_allenoate", "reactant", [E[2]]), participant("p_product", product, "product", [E[2]], product_role="cited product family")]
    r["notes"] = [f(detail, detail, [E[2]])]
    classify(r, [E[2]], coverage="unresolved")

r = reaction("generic_divergent_scheme", [E[2], E[3]], "generic_scheme", [], "Base-enabled divergent reactivity of isatins and allenoates", "generic divergent annulation/gamma-functionalization")
r["participants"] = [participant("p_isatin", "m_generic_isatin", "reactant", [E[2], E[3]]), participant("p_allenoate", "m_generic_allenoate", "reactant", [E[2], E[3]]), participant("p_spiro", "m_generic_spiro", "product", [E[2]], product_role="DBU branch"), participant("p_gamma", "m_generic_gamma", "product", [E[2], E[3]], product_role="DABCO branch")]
r["notes"] = [f("DBU gives spirofuran oxindoles; DABCO gives gamma-functionalized allenoates.", "DBU gives spirofuran oxindoles; DABCO gives gamma-functionalized allenoates.", [E[2], E[3]])]
classify(r, [E[2], E[3]], coverage="unresolved")

# Table 1: every row is an independent optimization experiment. Entry 25 also is Scheme 2 product 3a.
table = [
    (1,1.0,"DBU",1.0,"m_mecn","r.t.",2,"yield",18), (2,1.5,"DBU",1.0,"m_mecn","r.t.",2,"yield",27),
    (3,2.0,"DBU",1.0,"m_mecn","r.t.",2,"yield",10), (4,3.0,"DBU",1.0,"m_mecn","r.t.",2,"trace","Trace"),
    (5,0.6,"DBU",1.0,"m_mecn","r.t.",2,"yield",43), (6,0.5,"DBU",1.0,"m_mecn","r.t.",2,"yield",32),
    (7,0.6,"DBU",1.5,"m_mecn","r.t.",2,"yield",49), (8,0.6,"DBU",2.0,"m_mecn","r.t.",2,"yield",56),
    (9,0.6,"DBU",2.5,"m_mecn","r.t.",2,"yield",28), (10,0.6,"DBU",0.5,"m_mecn","r.t.",2,"trace","Trace"),
    (11,0.6,"DMAP",2.0,"m_mecn","r.t.",24,"dash","—"), (12,0.6,"DABCO",2.0,"m_mecn","r.t.",24,"gamma",36),
    (13,0.6,"Et3N",2.0,"m_mecn","r.t.",24,"dash","—"), (14,0.6,"PPh3",2.0,"m_mecn","r.t.",24,"dash","—"),
    (15,0.6,"DBU",2.0,"m_toluene","r.t.",2,"trace","Trace"), (16,0.6,"DBU",2.0,"m_dioxane","r.t.",2,"trace","Trace"),
    (17,0.6,"DBU",2.0,"m_dcm","r.t.",2,"yield",15), (18,0.6,"DBU",2.0,"m_thf","r.t.",2,"yield",7),
    (19,0.6,"DBU",2.0,"m_etoh","r.t.",2,"dash","—"), (20,0.6,"DBU",2.0,"m_benzonitrile","r.t.",2,"yield",22),
    (21,0.6,"DBU",2.0,None,"r.t.",2,"trace","Trace"), (22,0.6,"DBU",2.0,"m_mecn",{"minimum":0,"maximum":5,"unit":"°C","raw_text":"0–5 °C"},2,"yield",24),
    (23,0.6,"DBU",2.0,"m_mecn",{"number":50,"unit":"°C"},2,"yield",31), (24,0.6,"DBU",2.0,"m_mecn","r.t.",0.5,"yield",33),
    (25,0.6,"DBU",2.0,"m_mecn","r.t.",1,"yield",75), (26,0.6,"DBU",2.0,"m_mecn","r.t.",1.5,"yield",66),
]
base_ids = {"DBU":"m_dbu","DABCO":"m_dabco","DMAP":"m_dmap","Et3N":"m_et3n","PPh3":"m_pph3"}
for n, aeq, base, beq, solvent, temp, hours, result, val in table:
    ev = [E[3]] + ([E[2]] if n in (1, 25) else [])
    suffix = f"table1_entry{n}"
    cats = ["optimization"] + (["substrate_scope"] if n == 25 else [])
    r = reaction(suffix, ev, "specific_experiment", cats, "Optimization of DBU-mediated spirofuran oxindole formation", "MBH/cascade annulation")
    pis = participant("p_isatin", isatin_ids["NMe"], "reactant", ev, eq=1.0, amount=0.31, ref=isatin_ids["NMe"])
    pal = participant("p_allenoate", allenoate_ids[("Et","Me")], "reactant", ev, eq=aeq, ref=isatin_ids["NMe"])
    pb = participant("p_base", base_ids[base], "reagent", ev, eq=beq, ref=isatin_ids["NMe"])
    parts = [pis, pal, pb]
    if solvent:
        parts.append(participant("p_solvent", solvent, "solvent", ev))
    if result in ("yield", "trace"):
        parts.append(participant("p_product_3a", prod3["3a"], "product", ev, product_role="desired spirofuran oxindole"))
    if result == "gamma":
        parts += [participant("p_desired_3a", prod3["3a"], "product", ev, product_role="desired spiro product not reported"), participant("p_product_4a", prod4["4a"], "product", ev, product_role="divergent gamma-functionalized product")]
    r["participants"] = parts
    r["reaction_scale"] = f("1: 0.31 mmol", {"reference_molecule_id": isatin_ids["NMe"], "amount": {"number": 0.31, "unit": "mmol"}}, ev)
    if solvent:
        r["solvent_conditions"]["composition"] = f(solvent.replace("m_", ""), [solvent], ev)
    else:
        r["solvent_conditions"]["composition"] = f("neat", [], ev)
    pids = [p["participant_id"] for p in parts if p["stage"] == "reaction"]
    timen = {"number": hours, "unit": "h"}
    r["steps"] = [step("reaction", ev, pids, temp=temp, time=timen, atmosphere="Ar")]
    if result == "yield":
        r["outcomes"] = [outcome("yield_3a", "yield", [prod3["3a"]], ev, f"{val}%", {"number": val, "unit": "%"}, basis="isolated yield", ytype="isolated")]
        classify(r, ev, ["C-C", "C-O", "O-H"])
    elif result == "trace":
        r["outcomes"] = [outcome("trace_3a", "trace", [prod3["3a"]], ev, "Trace", "Trace", basis="isolated yield")]
        classify(r, ev, ["C-C", "C-O", "O-H"])
    elif result == "gamma":
        r["outcomes"] = [outcome("table_dash_3a", "other", [prod3["3a"]], ev, "—", "table dash (—)", scope="desired spiro product"), outcome("yield_4a", "yield", [prod4["4a"]], ev, "36%", {"number":36,"unit":"%"}, basis="isolated yield", ytype="isolated")]
        r["notes"] = [f("Footnote b: DABCO afforded product 4a in 36% yield.", "DABCO afforded product 4a in 36% yield.", ev)]
        classify(r, ev, ["C-C", "O-H"])
    else:
        r["outcomes"] = [outcome("table_dash", "other", [], ev, "—", "table dash (—)", scope="yield column")]
        classify(r, ev, coverage="not_applicable", limitation="The table dash is preserved and is not converted to zero or a claimed no-reaction result.")
    if n == 25:
        r["general_scheme_ids"] = [f"{DOC}:generic_divergent_scheme"]
        r["notes"] = [f("Same experiment appears as Table 1 entry 25 and Scheme 2 product 3a; merged rather than double-counted.", "Same experiment appears as Table 1 entry 25 and Scheme 2 product 3a; merged rather than double-counted.", ev)]

# Remaining Scheme 2 scope experiments (3a is already Table 1 entry 25).
for lab, ikey, r2, r3, y in scope3:
    if lab == "3a":
        continue
    ev = [E[2]]
    r = reaction(f"scheme2_{lab}", ev, "specific_experiment", ["substrate_scope"], "Scope of the reaction for spirofuran oxindoles", "DBU-mediated MBH/cascade annulation")
    parts = [participant("p_isatin", isatin_ids[ikey], "reactant", ev, eq=1.0, ref=isatin_ids[ikey]),
             participant("p_allenoate", allenoate_ids[(r2,r3)], "reactant", ev, eq=0.6, ref=isatin_ids[ikey]),
             participant("p_dbu", "m_dbu", "reagent", ev, eq=2.0, ref=isatin_ids[ikey]),
             participant("p_mecn", "m_mecn", "solvent", ev),
             participant("p_product", prod3[lab], "product", ev, product_role="reported spirofuran oxindole")]
    r["participants"] = parts
    r["solvent_conditions"]["composition"] = f("Acetonitrile", ["m_mecn"], ev)
    r["steps"] = [step("reaction", ev, [p["participant_id"] for p in parts], temp="r.t.", time={"number":1,"unit":"h"}, atmosphere="Ar")]
    r["procedure"]["general_procedure_reference"] = f("optimized conditions stated in page-2 prose", "optimized conditions stated in page-2 prose", ev)
    r["outcomes"] = [outcome(f"yield_{lab}", "yield", [prod3[lab]], ev, f"{y}%", {"number":y,"unit":"%"})]
    r["general_scheme_ids"] = [f"{DOC}:generic_divergent_scheme"]
    r["review"]["unresolved"].append("Scheme 2 does not state the yield measurement basis on the supplied pages.")
    classify(r, ev, ["C-C", "C-O", "O-H"])

# Explicitly reported non-cyclizing scope limitations.
for suffix, ikey, akey, wording in [
    ("limit_5OMe_isatin", "NMe_5OMe", ("Et","Me"), "5-OMe isatin did not furnish the expected cyclized product"),
    ("limit_5NO2_isatin", "NMe_5NO2", ("Et","Me"), "5-NO2 isatin did not furnish the expected cyclized product"),
    ("limit_gamma_phenyl_allenoate", "NMe", ("Et","Ph"), "gamma-phenyl allenoate did not furnish the expected cyclized product"),
]:
    ev=[E[2]]
    r=reaction(suffix,ev,"specific_experiment",["substrate_scope"],"Reported substrate limitation","DBU-mediated annulation attempt")
    parts=[participant("p_isatin",isatin_ids[ikey],"reactant",ev,eq=1.0,ref=isatin_ids[ikey]),participant("p_allenoate",allenoate_ids[akey],"reactant",ev,eq=0.6,ref=isatin_ids[ikey]),participant("p_dbu","m_dbu","reagent",ev,eq=2.0,ref=isatin_ids[ikey]),participant("p_mecn","m_mecn","solvent",ev)]
    r["participants"]=parts;r["solvent_conditions"]["composition"]=f("Acetonitrile",["m_mecn"],ev)
    r["steps"]=[step("reaction",ev,[p["participant_id"] for p in parts],temp="r.t.",time={"number":1,"unit":"h"},atmosphere="Ar")]
    r["outcomes"]=[outcome("desired_product_failure","other",[],ev,wording,wording,scope="expected cyclized product")]
    r["notes"]=[f(writing:=wording,writing,ev)]
    r["review"]["unresolved"].append("Conversion, recovery, byproducts, and detection method are unreported.")
    classify(r,ev,coverage="not_applicable",limitation="No observed cyclized product; expected bond formation is not assigned.")

# Scheme 3 scope, with the printed temperature conflict preserved for every record.
for lab, ikey, r2, r3, y, dr in scope4:
    ev=[E[3]]
    r=reaction(f"scheme3_{lab}",ev,"specific_experiment",["substrate_scope"],"Generality of the reaction for gamma-functionalized allenoate oxindoles","DABCO-mediated gamma-functionalization")
    parts=[participant("p_isatin",isatin_ids[ikey],"reactant",ev,eq=1.0,ref=isatin_ids[ikey]),participant("p_allenoate",allenoate_ids[(r2,r3)],"reactant",ev,eq=1.5,ref=isatin_ids[ikey]),participant("p_dabco","m_dabco","reagent",ev,eq=1.8,ref=isatin_ids[ikey]),participant("p_dmso","m_dmso","solvent",ev),participant("p_product",prod4[lab],"product",ev,product_role="reported gamma-functionalized allenoate")]
    r["participants"]=parts;r["solvent_conditions"]["composition"]=f("DMSO",["m_dmso"],ev)
    s=step("reaction",ev,[p["participant_id"] for p in parts],time={"number":0.5,"unit":"h"},atmosphere="Ar")
    s["temperature"]=conflict_temperature(ev);r["steps"]=[s]
    r["outcomes"]=[outcome(f"yield_{lab}","yield",[prod4[lab]],ev,f"{y}%",{"number":y,"unit":"%"}),outcome(f"dr_{lab}","selectivity",[prod4[lab]],ev,f"dr {dr}",dr,scope={"metric":"dr","reported_order":"1:x","isomer identities":"unassigned"})]
    r["general_scheme_ids"]=[f"{DOC}:generic_divergent_scheme"]
    r["notes"]=[f("Scheme 3 header prints DMSO, 0–5 °C and also Ar, rt; prose states 0–5 °C for 0.5 h.","Scheme 3 header prints DMSO, 0–5 °C and also Ar, rt; prose states 0–5 °C for 0.5 h.",ev)]
    r["review"]["unresolved"] += ["Temperature conflict (0–5 °C versus rt) is preserved without selecting a value.","Yield and dr measurement methods and the identities/order of the diastereomers are unreported."]
    classify(r,ev,["C-C","O-H"])

# Scheme 4 mechanistic controls.  Page 4 continues the interpretation of the
# page-3 drawings, so both source pages are associated without duplicating the
# experiments.
ev=[E[3], E[4]]
r=reaction("scheme4a_intermediate_then_cyclization",ev,"specific_experiment",["mechanistic_study","control"],"Rapid formation of 4a followed by cyclization","sequential gamma-functionalization/cascade cyclization")
parts=[participant("p_isatin",isatin_ids["NMe"],"reactant",ev),participant("p_allenoate",allenoate_ids[("Et","Me")],"reactant",ev),participant("p_dbu","m_dbu","reagent",ev),participant("p_mecn","m_mecn","solvent",ev),participant("p_4a",prod4["4a"],"product",ev,product_role="observed intermediate product"),participant("p_3a",prod3["3a"],"product",ev,product_role="cyclization product")]
r["participants"]=parts;r["solvent_conditions"]["composition"]=f("acetonitrile",["m_mecn"],ev)
r["steps"]=[step("rapid_stage",ev,["p_isatin","p_allenoate","p_dbu","p_mecn","p_4a"],temp="r.t.",time={"number":1,"unit":"min","qualifier":"gt","raw_text":">1 min"},atmosphere="Ar",operation="formation of 4a observed before cyclization"),step("cyclization_stage",ev,["p_4a","p_isatin","p_dbu","p_3a"],temp="r.t.",atmosphere="Ar",operation="subsequent cyclization to 3a; duration unreported")]
r["outcomes"]=[outcome("yield_4a","yield",[prod4["4a"]],ev,"63%",{"number":63,"unit":"%"}),outcome("formation_3a","other",[prod3["3a"]],ev,"cyclization to 3a","cyclization to 3a")]
r["review"]["unresolved"] += ["DBU loading, substrate charges, 4a yield method, and duration/extent of the second stage are unreported."]
classify(r,ev,["C-C","C-O","O-H"])

r=reaction("scheme4b_isolated_4a_cyclization",ev,"specific_experiment",["mechanistic_study","control"],"DBU-mediated cyclization of isolated 4a with isatin 1a","cascade cyclization")
parts=[participant("p_4a",prod4["4a"],"reactant",ev,amount=0.17),participant("p_isatin",isatin_ids["NMe"],"reactant",ev,amount=0.17),participant("p_dbu","m_dbu","reagent",ev,eq=2.0,ref=isatin_ids["NMe"]),participant("p_mecn","m_mecn","solvent",ev),participant("p_3a",prod3["3a"],"product",ev,product_role="reported cyclization product")]
r["participants"]=parts;r["reaction_scale"]=f("4a and 1a: 0.17 mmol each",{"reference_molecule_id":prod4["4a"],"amount":{"number":0.17,"unit":"mmol"}},ev);r["solvent_conditions"]["composition"]=f("acetonitrile",["m_mecn"],ev)
r["steps"]=[step("reaction",ev,[p["participant_id"] for p in parts],temp="r.t.",time={"number":1,"unit":"h"},atmosphere="Ar")]
r["outcomes"]=[outcome("yield_3a","yield",[prod3["3a"]],ev,"78%",{"number":78,"unit":"%"})];r["review"]["unresolved"].append("Yield measurement basis is unreported.");classify(r,ev,["C-C","C-O","O-H"])

r=reaction("scheme4c_dabco_no_cyclization",ev,"specific_experiment",["mechanistic_study","control"],"DABCO reaction giving 4a without further cyclization","gamma-functionalization control")
parts=[participant("p_isatin",isatin_ids["NMe"],"reactant",ev),participant("p_allenoate",allenoate_ids[("Et","Me")],"reactant",ev),participant("p_dabco","m_dabco","reagent",ev),participant("p_dmso","m_dmso","solvent",ev),participant("p_4a",prod4["4a"],"product",ev,product_role="observed gamma-functionalized product"),participant("p_3a",prod3["3a"],"product",ev,product_role="desired cyclization product not formed")]
r["participants"]=parts;r["solvent_conditions"]["composition"]=f("DMSO",["m_dmso"],ev);r["steps"]=[step("reaction",ev,[p["participant_id"] for p in parts],temp={"minimum":0,"maximum":5,"unit":"°C","raw_text":"0–5 °C"},time={"number":24,"unit":"h"},atmosphere="Ar")]
r["outcomes"]=[outcome("yield_4a","yield",[prod4["4a"]],ev,"52%",{"number":52,"unit":"%"}),outcome("no_cyclization_3a","no_reaction",[prod3["3a"]],ev,"cyclization to 3a not observed","cyclization to 3a not observed",scope="further cyclization only")]
r["review"]["unresolved"] += ["DABCO loading, substrate charges, and yield measurement basis are unreported."];classify(r,ev,["C-C","O-H"])

# Proposed mechanism, kept separate from executed controls.
ev=[E[4]]
r=reaction("scheme5_proposed_mechanism",ev,"proposed_mechanism",[],"Plausible mechanistic pathway","author-proposed MBH/cascade mechanism")
r["participants"]=[participant("p_isatin",isatin_ids["NMe"],"reactant",ev),participant("p_allenoate",allenoate_ids[("Et","Me")],"reactant",ev),participant("p_dbu","m_dbu","reagent",ev),participant("p_4a",prod4["4a"],"product",ev,product_role="proposed/observed intermediate product"),participant("p_3a",prod3["3a"],"product",ev,product_role="final product")]
r["notes"]=[f("Author-proposed sequence A→B→C→4a→D→E→F→3a; intermediates A–F are not treated as isolated experiments.","Author-proposed sequence A→B→C→4a→D→E→F→3a; intermediates A–F are not treated as isolated experiments.",ev)]
classify(r,ev,coverage="unresolved",limitation="Proposed mechanism is not a complete executed experiment or atom-mapped net reaction.")

# Non-synthetic context and explicit negative inventory.
B["source_context"] = [
    {"context_id": f"{DOC}:metadata_and_bioactive_context", "document_id": DOC, "kind": "metadata_and_non_synthetic_context", "evidence_ids": [E[1]], "note": "Title/DOI/abstract and Figure 1 bioactive scaffolds; no current synthesis experiment."},
    {"context_id": f"{DOC}:structure_assignment", "document_id": DOC, "kind": "characterization", "evidence_ids": [E[2], E[3]], "note": "3a relative stereochemistry was assigned by NMR/HRMS; 3ag and 4h were confirmed by single-crystal X-ray (CCDC 1979517 and 1979628). External CIFs were not supplied."},
    {"context_id": f"{DOC}:mechanistic_discussion", "document_id": DOC, "kind": "mechanistic_discussion", "evidence_ids": [E[3], E[4]], "note": "Author discussion of DBU/DABCO equilibrium and beta-ammonium enolate stabilization; executed controls are separate reaction records."},
    {"context_id": f"{DOC}:negative_spectral_inventory", "document_id": DOC, "kind": "negative_inventory", "evidence_ids": [E[1], E[2], E[3], E[4]], "note": "No standalone mass-spectrum peak table, molecular-weight table, or fragment-relative-abundance block occurs on pages 1–4."},
]

# Save raw inventory and stereo diagnostics separately from submitted records.
raw_inventory = {
    "document_id": DOC,
    "pages": {
        "1": {"regions": ["title/abstract", "introduction", "Figure 1"], "reaction_decision": "no current concrete experiment"},
        "2": {"regions": ["Scheme 1", "Scheme 2", "optimization prose"], "inventory": {"literature_examples": 3, "generic_schemes": 1, "scheme2_products": 36, "explicit_failed_scope_mentions": 3}},
        "3": {"regions": ["Table 1", "Scheme 3", "Scheme 4", "mechanistic prose"], "inventory": {"table_rows": 26, "scheme3_products": 19, "controls": 3}},
        "4": {"regions": ["Scheme 5", "mechanistic discussion", "references"], "inventory": {"proposed_mechanisms": 1}},
    },
    "deduplication": ["Table 1 entry 1 includes the page-2 initial 18% feasibility report.", "Table 1 entry 25 and Scheme 2 product 3a are one experiment."],
    "limitations": ["No SI was supplied.", "Compact structures are retained unresolved where safe full encoding was not achieved.", "Scheme 3 temperature conflict is preserved."],
}
(RUN / "raw_model/page_extraction.json").write_text(json.dumps(raw_inventory, ensure_ascii=False, indent=2) + "\n")
diagnostics = {}
for m in B["molecules"]:
    smi = st.value(m["smiles"])
    if smi:
        diagnostics[m["molecule_id"]] = st.stereo_inspect(smi)
(RUN / "raw_model/stereo_checks.json").write_text(json.dumps(diagnostics, ensure_ascii=False, indent=2) + "\n")
(RUN / "submitted.json").write_text(json.dumps(B, ensure_ascii=False, indent=2) + "\n")

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 1,
            "state": "attempted",
            "regions": [
                {"region": "title, citation, abstract, and ESI notice", "finding": "Metadata and article summary; ESI is advertised but no SI file was supplied or associated.", "evidence_ids": [E[1]]},
                {"region": "two-column introduction", "finding": "Background and cited literature discussion; no current concrete experiment arrow or table row.", "evidence_ids": [E[1]]},
                {"region": "Figure 1", "finding": "Selected bioactive oxindole/allene scaffolds; non-synthetic context, not experiments from this paper.", "evidence_ids": [E[1]]},
            ],
            "record_ids": [],
            "unresolved": ["No SI was supplied despite the availability notice; CCDC files were not treated as attached sources."],
        },
        {
            "page": 2,
            "state": "attempted",
            "regions": [
                {"region": "Scheme 1", "finding": "Three cited prior-work reaction families plus one current generic divergent scheme; classified separately from current concrete experiments.", "evidence_ids": [E[2]]},
                {"region": "left-column feasibility and optimization prose", "finding": "The 18% feasibility result duplicates Table 1 entry 1; the optimized 75% result duplicates Table 1 entry 25/Scheme 2 product 3a.", "evidence_ids": [E[2]]},
                {"region": "Scheme 2", "finding": "All 36 printed products 3a–3aj inventoried; 3a is cross-page-deduplicated with Table 1 entry 25, and 35 additional scope records were created.", "evidence_ids": [E[2]]},
                {"region": "right-column scope discussion", "finding": "Three explicit failed scope limitations (5-OMe isatin, 5-NO2 isatin, gamma-phenyl allenoate) retained without fabricated products.", "evidence_ids": [E[2]]},
            ],
            "record_ids": [r["reaction_id"] for r in B["reactions"] if E[2] in r["evidence_ids"]],
            "unresolved": ["Scheme 2 does not state yield measurement methods for products other than the deduplicated 3a table experiment.", "Most compact product structures cannot be safely encoded as complete SMILES; source labels/drawings and stereo status are retained."],
        },
        {
            "page": 3,
            "state": "attempted",
            "regions": [
                {"region": "Table 1 and footnotes", "finding": "All 26 optimization rows inventoried separately; dashes and trace values preserved qualitatively; entry 12 separately records 4a at 36% isolated yield.", "evidence_ids": [E[3]]},
                {"region": "Scheme 3", "finding": "All 19 gamma-functionalized allenoate outcomes 4a–4s recorded with printed yields and dr values.", "evidence_ids": [E[3]]},
                {"region": "Scheme 3 header and optimization prose", "finding": "Temperature conflict retained: the scheme prints both 0–5 °C and rt, while prose states 0–5 °C for 0.5 h.", "evidence_ids": [E[3]]},
                {"region": "Scheme 4", "finding": "Three distinct mechanistic controls recorded: rapid 4a formation then cyclization, isolated-4a cyclization, and DABCO formation of 4a without further cyclization.", "evidence_ids": [E[3]]},
            ],
            "record_ids": [r["reaction_id"] for r in B["reactions"] if E[3] in r["evidence_ids"]],
            "unresolved": ["Scheme 3 yield/dr measurement methods and dr component identities/order are unreported.", "Several Scheme 4 reagent charges and yield methods are unreported."],
        },
        {
            "page": 4,
            "state": "attempted",
            "regions": [
                {"region": "Scheme 5", "finding": "One author-proposed mechanism A–F retained as a proposed-mechanism record; intermediates were not expanded into executed experiments.", "evidence_ids": [E[4]]},
                {"region": "left-column continuation of control interpretation", "finding": "Cross-page discussion supports the Scheme 4 controls and DBU/DABCO divergence; linked to existing controls rather than duplicated.", "evidence_ids": [E[4]]},
                {"region": "summary, acknowledgements, conflicts, and references", "finding": "Non-experimental text; no additional current synthesis arrow, optimization row, or procedure.", "evidence_ids": [E[4]]},
            ],
            "record_ids": [r["reaction_id"] for r in B["reactions"] if E[4] in r["evidence_ids"]],
            "unresolved": ["The gamma-functionalized-allenoate mechanism in Scheme S2 is mentioned but unavailable because SI was not supplied.", "No standalone mass-spectrum peak table, molecular-weight table, or fragment-abundance block was found on pages 1–4."],
        },
    ],
}
(RUN / "coverage").mkdir(exist_ok=True)
(RUN / "coverage/pages_0001_0004.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"molecules": len(B["molecules"]), "reactions": len(B["reactions"]), "source_context": len(B["source_context"])}, ensure_ascii=False))
