#!/usr/bin/env python3
"""Source-derived v3.1 extraction for doc_5dd58f28f6bfaa1f.

This file is a machine extraction artifact. It deliberately preserves unknown
structures and unavailable SI details; it is not an expert-reviewed record.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_5dd58f28f6bfaa1f"
RUN = ROOT / "runs/current"
DATASET = json.loads((RUN / "dataset.json").read_text())
E = {p: f"{DOC}:page{p}:image" for p in range(1, 6)}
T = {p: f"{DOC}:page{p}:text" for p in range(1, 6)}


def F(raw=None, norm=None, ev=(), status=None):
    if norm is None and raw is not None:
        norm = raw
    return st.fact(raw, norm, ev, status=status)


def I(norm, ev, why, raw=None):
    return st.fact(raw, norm, ev, method="inference", derivation=why)


def D(norm, ev, why, raw=None):
    return st.fact(raw, norm, ev, method="structure_derivation", derivation=why)


def M():
    return st.fact()


def conflict(raw, alternatives, ev):
    f = st.fact(raw, None, ev, status="conflict")
    f["alternatives"] = alternatives
    return f


bundle = st.empty(DATASET["documents"])
molecules = {}


def add_molecule(label, ev, name=None, smiles=None, kind="unresolved", note=None,
                 chirality="unresolved", geometry="unresolved", representation=None,
                 axial=None, mixture=None, assignment=None):
    mid = f"{DOC}:mol:{label}"
    if mid in molecules:
        return mid
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, ev)
    m["names"]["original"] = F(name or label, name or label, ev)
    m["names"]["english"] = F(name or label, name or label, ev)
    m["names"]["abbreviations"] = F(label, [label], ev) if label != (name or label) else M()
    m["names"]["aliases"] = M()
    m["cas_numbers"] = M()
    m["identity_note"] = F(note, note, ev) if note else M()
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = D(smiles, ev, "Candidate connectivity transcribed from the source drawing; RDKit parsing is not source verification.", raw="source structure drawing")
        m["generic_template"] = None
    else:
        m["representation_type"] = kind
        m["smiles"] = M()
        m["generic_template"] = None
    rep = representation or ("partial" if smiles else "unresolved")
    for key, value in (("chirality_status", chirality), ("geometry_status", geometry), ("representation_status", rep)):
        m["stereochemistry"][key] = I(value, ev, "Stereochemical status assessed from the supplied source drawing; missing labels were not converted to achiral or racemic.", raw="source drawing")
    if axial is not None:
        m["stereochemistry"]["axial"] = F(axial, axial, ev)
    if mixture is not None:
        m["stereochemistry"]["mixture"] = F(mixture, mixture, ev)
    if assignment is not None:
        m["stereochemistry"]["assignment"] = F(assignment, assignment, ev)
    molecules[mid] = m
    return mid


def q(number, unit, **extra):
    return dict(number=number, unit=unit, **extra)


def qraw(value):
    if not isinstance(value, dict):
        return str(value)
    if value.get("raw_text"):
        return value["raw_text"]
    if "number" in value:
        prefix = {"lt":"<", "le":"<=", "gt":">", "ge":">=", "approx":"~"}.get(value.get("qualifier"), "")
        return f"{prefix}{value['number']} {value['unit']}"
    return f"{value.get('minimum')}–{value.get('maximum')} {value.get('unit')}"


def participant(r, pid, mid, role, ev, amount=None, mass=None, volume=None,
                equiv=None, mol_percent=None, ref=None, product_role=None):
    p = st.new_participant(pid, mid, role, ev)
    if product_role:
        p["product_role"] = F(product_role, product_role, ev)
    for key, value in (("amount", amount), ("mass", mass), ("volume", volume),
        ("equivalents", equiv), ("mol_percent", mol_percent)):
        if value is not None:
            p["charges"][key] = F(qraw(value), value, ev)
    if ref is not None:
        p["equivalence_reference"] = F(ref, ref, ev)
    r["participants"].append(p)
    return p


def classify_cc(r, ev, detail="At least one new C-C framework connection is evident from reactant/product drawings; atom-complete mapping is unavailable."):
    change = {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
              "before_order": 0, "after_order": 1, "change_type": "formed"}
    label = {"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1"]}
    r["classification"]["bond_changes"] = [D(change, ev, detail)]
    r["classification"]["labels"] = D([label], ev, detail)
    r["classification"]["coverage"] = I("partial", ev, detail)
    r["classification"]["limitations"] = [
        "No atom-mapped complete structures were supplied or safely reconstructed.",
        "Only the source-evident C-C framework formation is asserted; other net bond-order, heteroatom, and hydrogen changes remain unresolved.",
    ]


def unresolved_class(r, ev, reason, not_applicable=False):
    r["classification"]["labels"] = I([], ev, reason)
    r["classification"]["coverage"] = I("not_applicable" if not_applicable else "unresolved", ev, reason)
    r["classification"]["limitations"] = [reason]


def reaction_base(suffix, ev, record_type, categories, location, name, transformation=None):
    rid = f"{DOC}:{suffix}"
    r = st.new_reaction(rid, DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = F(categories, categories, ev) if categories is not None else M()
    r["location"] = F(location, location, ev)
    r["reaction_name"] = F(name, name, ev) if name else M()
    r["transformation_type"] = I(transformation, ev, "Normalized transformation description based on the source scheme and prose.", raw=name) if transformation else M()
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    return r


def add_step(r, sid, ev, pids, operation=None, temperature=None, time=None):
    s = st.new_step(sid, ev)
    s["participant_ids"] = list(pids)
    if operation is not None:
        s["operation"] = F(operation, operation, ev)
    if temperature is not None:
        s["temperature"] = F(qraw(temperature), temperature, ev)
    if time is not None:
        s["time"] = F(qraw(time), time, ev)
    r["steps"].append(s)
    return s


def yield_outcome(r, oid, mid, value, ev, isolated=True, conflict_value=None):
    o = st.new_outcome(oid, "yield", [mid], ev)
    o["product_role"] = F("desired product", "desired product", ev)
    if conflict_value is None:
        o["value"] = F(qraw(value), value, ev)
    else:
        o["value"] = conflict([value, conflict_value], [
            F(str(value), value, [E[2]]), F(str(conflict_value), conflict_value, [T[3]])
        ], ev)
    if isolated:
        o["yield_type"] = F("isolated yield", "isolated", ev)
        o["measurement_basis"] = F("Isolated yield", "isolated", ev)
    r["outcomes"].append(o)
    return o


def selectivity_outcome(r, oid, mid, raw, metric, ev, method=None, conflict_fact=None):
    o = st.new_outcome(oid, "selectivity", [mid], ev)
    o["product_role"] = F("desired product", "desired product", ev)
    o["value"] = conflict_fact if conflict_fact is not None else F(raw, raw, ev)
    o["measurement_basis"] = F(method, method, ev) if method else M()
    o["scope"] = F({"metric": metric, "reported_order": metric.split("/")},
                   {"metric": metric, "reported_order": metric.split("/")}, ev)
    r["outcomes"].append(o)
    return o


# Core paper molecules and standard-condition agents.
m1a = add_molecule("1a", [E[2]], "nitrone 1a",
    "[O-][N+](c1ccccc1)=C(c1ccccc1)C=Cc1ccccc1",
    note="Connectivity candidate only; the two potential alkene/nitrone geometries are not encoded.",
    chirality="unresolved", geometry="unspecified", representation="partial")

allenoate_smiles = {
    "2a": "CC=C=CC(=O)OCc1ccccc1", "2b": "CC=C=CC(=O)OC",
    "2c": "CC=C=CC(=O)OCCCC", "2d": "CC=C=CC(=O)OC(C)C",
    "2e": "CC=C=CC(=O)OC(C)(C)C", "2f": "CC=C=CC(=O)Oc1ccccc1",
    "2g": "CC=C=CC(=O)OCC=Cc1ccccc1", "2h": "C=C=CC(=O)OCc1ccccc1",
    "2i": "CCCCC=C=CC(=O)OCc1ccccc1", "2j": "c1ccccc1CC=C=CC(=O)OCc1ccccc1",
    "2k": "CC(C)C=C=CC(=O)OCc1ccccc1", "2l": "CC(C)(C)C=C=CC(=O)OCc1ccccc1",
    "2m": "c1ccccc1C=C=CC(=O)OCc1ccccc1",
}
allenoates = {}
for label in [f"2{x}" for x in "abcdefghijklmn"]:
    ev = [E[2]] if label == "2a" else [E[3]]
    if label == "2h":
        ch, ax = "achiral", {"axis": "allene", "status": "not stereogenic because one terminal carbon bears two H in the drawn structure"}
    elif label == "2n":
        ch, ax = "partially_specified", {"axis": "allene", "status": "configuration not assigned; (-)-menthyl ester contains specified source chirality"}
    else:
        ch, ax = "unspecified", {"axis": "allene", "status": "potential axial stereochemistry; configuration not assigned"}
    allenoates[label] = add_molecule(label, ev, f"allenoate {label}", allenoate_smiles.get(label),
        note="Allene axial configuration is not reported; toolkit diagnostics do not establish axial stereochemistry." if label != "2h" else "Unsubstituted allene terminus removes allene axial stereogenicity in the drawn scaffold.",
        chirality=ch, geometry="unspecified", representation="partial" if label in allenoate_smiles else "not_encoded", axial=ax)

common_specs = {
    "NiCl2": ("nickel(II) chloride", "[Cl-].[Cl-].[Ni+2]"),
    "FeCl3": ("iron(III) chloride", "Cl[Fe](Cl)Cl"),
    "ZnCl2": ("zinc chloride", "Cl[Zn]Cl"),
    "L1": ("2,2'-bipyridine (L1)", "n1ccccc1-c1ccccn1"),
    "L2": ("1,10-phenanthroline (L2)", "c1cnc2c(c1)ccc1cccnc12"),
    "L3": ("2,9-dimethyl-1,10-phenanthroline (L3)", "Cc1ccc2c(c1)nc1c(C)cccc1c2"),
    "MgSO4": ("magnesium sulfate", "[Mg+2].[O-]S(=O)(=O)[O-]"),
    "NaOH": ("sodium hydroxide", "[Na+].[OH-]"),
    "toluene": ("toluene", "Cc1ccccc1"), "THF": ("tetrahydrofuran", "C1CCOC1"),
    "DCE": ("1,2-dichloroethane", "ClCCCl"), "DCM": ("dichloromethane", "ClCCl"),
    "MeCN": ("acetonitrile", "CC#N"), "DMSO": ("dimethyl sulfoxide", "CS(C)=O"),
    "iPrOAc": ("isopropyl acetate", "CC(=O)OC(C)C"), "PhH": ("benzene", "c1ccccc1"),
}
common = {}
for label, (name, smi) in common_specs.items():
    ev = [E[1]] if label in ("PhH", "iPrOAc") else [E[2]]
    common[label] = add_molecule(label, ev, name, smi, chirality="achiral", geometry="not_applicable", representation="not_applicable")
common["L4"] = add_molecule("L4", [E[2]], "Pybox ligand L4", note="Full source drawing retained; unique SMILES not assigned.", chirality="unresolved", geometry="not_applicable")

# Scope molecules: unresolved dense structures are kept independent by source label.
nitrones1 = {"1a": m1a}
for x in "bcdefghijklmnopqrstuvwx":
    label = f"1{x}"
    nitrones1[label] = add_molecule(label, [E[2]], f"N-aryl nitrone {label}",
        note="Dense scope structure retained by source label; complete SMILES and alkene geometry require expert structure review.",
        chirality="unresolved", geometry="unspecified", representation="not_encoded")

nitrones4 = {}
for x in "abcdef":
    label = f"4{x}"
    nitrones4[label] = add_molecule(label, [E[3]], f"N-vinyl nitrone {label}",
        note="Dense scope structure retained by source label; complete SMILES and geometry require expert structure review.",
        chirality="unresolved", geometry="unspecified", representation="not_encoded")

products3 = {}
for label in ([f"3{x}a" for x in "abcdefghijklmnopqrstuvwx"] + [f"3a{x}" for x in "bcdefghijklmn"]):
    ev = [E[2]] if label.endswith("a") else [E[3]]
    if label == "3aa":
        geom, mix, assign = "mixture", "Experiment-specific Z/E ratios range from 1:1 to 20:1; the major Z assignment is supported by NOESY.", "Major enamine isomer assigned Z by 2D NOESY; Table 1 entry 9 prose/table disagree on single-isomer versus 20:1."
    elif label == "3an":
        geom, mix, assign = "partially_specified", "1:1 diastereomeric mixture in the reported experiment.", "A single enamine geometry is drawn without an explicit E/Z descriptor; the experiment reports dr = 1:1."
    else:
        geom, mix, assign = "partially_specified", None, "A single enamine geometry is drawn; no product-specific E/Z descriptor or ratio is printed."
    products3[label] = add_molecule(label, ev, f"N-vinylindole {label}",
        note="Product structure is source-identified but not converted to a complete checked SMILES.",
        chirality="mixture" if label == "3an" else "unresolved", geometry=geom, representation="not_encoded",
        mixture=mix, assignment=assign)

products5 = {}
for label in ("5aa", "5ba", "5ca", "5da", "5ea", "5fa", "5ab", "5af", "5ag"):
    products5[label] = add_molecule(label, [E[3]], f"N-vinylpyrrole {label}",
        note="Product structure is source-identified but not converted to a complete checked SMILES.",
        chirality="unresolved", geometry="partially_specified", representation="not_encoded",
        assignment="A single enamine geometry is drawn; no product-specific E/Z ratio is printed.")


# Scheme 1: three prior-work arrows and one current-work generic scheme.
litA_r = add_molecule("litA_benzopyranone_nitrone", [E[1]], "benzopyranone nitrone (Scheme 1A)")
litA_a = add_molecule("litA_allenoate_family", [E[1]], "allenoate family (Scheme 1A)", chirality="unspecified", geometry="unspecified", representation="not_encoded", axial={"axis":"allene","status":"R-substituted family; configuration unspecified"})
litA_p = add_molecule("litA_benzoindolizine_product", [E[1]], "benzoindolizine product (Scheme 1A)")
r = reaction_base("scheme1A_prior_ishar", [E[1]], "literature_example", ["other"], "Scheme 1A, page 1", "metal-free cycloaddition of benzopyranone nitrone with allenoates", "prior-work cascade cycloaddition")
participant(r,"nitrone",litA_r,"reactant",[E[1]]); participant(r,"allenoate",litA_a,"reactant",[E[1]]); participant(r,"benzene",common["PhH"],"solvent",[E[1]]); participant(r,"product",litA_p,"product",[E[1]],product_role="major product")
add_step(r,"reflux",[E[1]],["nitrone","allenoate","benzene","product"],"metal-free reaction","reflux",None)
yield_outcome(r,"yield",litA_p,q(60,"%"),[E[1]],isolated=False)
unresolved_class(r,[E[1]],"Prior-work structures are generic and no atom mapping is available.")
r["review"]["unresolved"]=["Yield measurement basis, time, quantities, and full identities are not reported in Scheme 1A."]
bundle["reactions"].append(r)

litB_r = add_molecule("litB_nitrone_family", [E[1]], "N-aryl alpha,beta-unsaturated nitrone family (Scheme 1B)")
litB_a = add_molecule("litB_unsubstituted_allenoate", [E[1]], "unsubstituted methyl allenoate (Scheme 1B)", chirality="achiral", geometry="unspecified", representation="not_encoded", axial={"axis":"allene","status":"not stereogenic at an unsubstituted terminus"})
for tag, solvent, pname in (("toluene","toluene","dihydrocarbazole major product"),("iproac","iPrOAc","dihydropyridoindole major product")):
    pm = add_molecule(f"litB_{tag}_product", [E[1]], pname)
    r = reaction_base(f"scheme1B_prior_{tag}", [E[1]], "literature_example", ["other"], "Scheme 1B, page 1", "metal-free solvent-controlled cycloaddition", "prior-work solvent-controlled cascade")
    participant(r,"nitrone",litB_r,"reactant",[E[1]]); participant(r,"allenoate",litB_a,"reactant",[E[1]]); participant(r,"solvent",common[solvent],"solvent",[E[1]]); participant(r,"product",pm,"product",[E[1]],product_role="major product")
    add_step(r,"reaction",[E[1]],["nitrone","allenoate","solvent","product"],"metal-free solvent-controlled reaction")
    o=st.new_outcome("major_product","other",[pm],[E[1]]); o["value"]=F("Major product","major product",[E[1]]); o["product_role"]=F("major product","major product",[E[1]]); r["outcomes"].append(o)
    unresolved_class(r,[E[1]],"Prior-work family structures and full conditions are not available in the current paper.")
    r["review"]["unresolved"]=["No yield, time, temperature, quantities, or measurement basis is printed for this Scheme 1B arrow."]
    bundle["reactions"].append(r)

fam1=add_molecule("family_1_nitrones",[E[1]],"alpha,beta-unsaturated nitrone family 1",geometry="unspecified")
fam2=add_molecule("family_2_allenoates",[E[1]],"allenoate family 2",chirality="unspecified",geometry="unspecified",representation="not_encoded",axial={"axis":"allene","status":"substituted family; configuration unspecified"})
fam3=add_molecule("family_3_products",[E[1]],"N-vinylindole family 3",geometry="partially_specified",representation="not_encoded",assignment="A single product geometry is drawn generically without an explicit E/Z descriptor.")
r = reaction_base("scheme1C_current_generic", [E[1]], "generic_scheme", None, "Scheme 1C and graphical abstract, page 1", "Ni(II)-catalyzed cycloaddition of nitrones with allenoates", "generic cascade [3+2] cycloaddition/rearrangement/cyclization")
participant(r,"nitrone_family",fam1,"reactant",[E[1]]); participant(r,"allenoate_family",fam2,"reactant",[E[1]]); participant(r,"product_family",fam3,"product",[E[1]],product_role="desired products")
add_step(r,"generic",[E[1]],["nitrone_family","allenoate_family","product_family"],"Ni(II)-controlled cascade")
classify_cc(r,[E[1]])
r["review"]["unresolved"]=["Generic R-group identities, quantities, and atom-complete bond changes are intentionally unresolved."]
bundle["reactions"].append(r)
generic_id=r["reaction_id"]


# Table 1 optimization rows.
table_rows = {
1:(None,None,"toluene",26,"1:1"),2:(None,None,"THF",32,"2:1"),3:(None,None,"DCE",42,"2:1"),4:(None,None,"DCM",58,"2:1"),
5:(None,None,"MeCN","<5",None),6:(None,None,"DMSO","<5",None),7:("FeCl3",None,"DCM",28,"10:1"),8:("ZnCl2",None,"DCM",40,"10:1"),
9:("NiCl2",None,"DCM",65,"20:1"),10:("NiCl2","L1","DCM",81,"20:1"),11:("NiCl2","L2","DCM",74,"20:1"),
12:("NiCl2","L3","DCM",51,"20:1"),13:("NiCl2","L4","DCM",56,"20:1"),14:("NiCl2","L1","DCM",83,"20:1"),15:("NiCl2","L1","DCM","<5",None)}
for n,(cat,lig,solv,yv,ratio) in table_rows.items():
    r = reaction_base(f"table1_entry{n:02d}",[E[2],T[2]],"specific_experiment",["optimization"],f"Table 1 entry {n}, page 2","Optimization of the reaction conditions","Ni(II)-controlled cascade to N-vinylindole 3aa")
    r["general_scheme_ids"]=[generic_id]
    participant(r,"1a",m1a,"reactant",[E[2]],amount=q(0.2,"mmol")); participant(r,"2a",allenoates["2a"],"reactant",[E[2]],amount=q(0.4,"mmol")); participant(r,"solvent",common[solv],"solvent",[E[2]],volume=q(2.0,"mL")); participant(r,"3aa",products3["3aa"],"product",[E[2]],product_role="desired product")
    if cat: participant(r,"catalyst",common[cat],"catalyst",[E[2]],mol_percent=q(20,"mol%"),ref=m1a)
    if lig: participant(r,"ligand",common[lig],"ligand",[E[2]],mol_percent=q(24,"mol%"),ref=m1a)
    if n==14: participant(r,"MgSO4",common["MgSO4"],"additive",[E[2]],mass=q(200,"mg"))
    if n==15: participant(r,"NaOH",common["NaOH"],"additive",[E[2]],equiv=q(2.0,"equiv"),ref=m1a)
    r["reaction_scale"]=F("1a (0.2 mmol)",{"reference_molecule_id":m1a,"amount":q(0.2,"mmol")},[E[2]])
    r["solvent_conditions"]["composition"]=F(solv,solv,[E[2]]); r["solvent_conditions"]["total_volume"]=F("2.0 mL",q(2.0,"mL"),[E[2]])
    tval=q(24,"h") if n==1 else {"minimum":16,"maximum":24,"unit":"h","raw_text":"16-24 h"}
    add_step(r,"reaction",[E[2],T[2]],[p["participant_id"] for p in r["participants"]],"reaction under Table 1 conditions","rt",tval)
    r["procedure"]["original_text"]=F("Table 1 reaction conditions: 1a (0.2 mmol), 2a (0.4 mmol), cat. (20 mol%), L (24 mol%), solvent (2.0 mL), rt, 16-24 h.","Table 1 footnote a",[E[2]])
    y = q(5,"%",qualifier="lt",raw_text="<5%") if yv=="<5" else q(yv,"%")
    yield_outcome(r,"yield_3aa",products3["3aa"],y,[E[2]],isolated=True)
    if ratio:
        if n==9:
            cf=conflict(["20:1","Z isomer only"],[F("20:1","20:1",[E[2]]),F("Z isomer only","Z only",[T[2]])],[E[2],T[2]])
            selectivity_outcome(r,"ze_3aa",products3["3aa"],ratio,"Z/E",[E[2],T[2]],method="2D NOESY spectra (assignment of Z isomer)",conflict_fact=cf)
            r["review"]["unresolved"].append("Entry 9 Table 1 reports Z/E 20:1, while prose says Z isomer only; both are retained as a conflict.")
        else: selectivity_outcome(r,"ze_3aa",products3["3aa"],ratio,"Z/E",[E[2]])
    if n in (5,6,15):
        unresolved_class(r,[E[2],T[2]],"Product formation was reported as <5% in the table and not observed/completely prohibited in prose; no formed-bond label is asserted.",not_applicable=True)
        r["notes"].append(F("3aa was not observed" if n in (5,6) else "NaOH completely prohibited formation of 3aa",None,[T[2]]))
        r["review"]["unresolved"].append("The table's <5% isolated-yield upper bound and prose no-product wording are both preserved; detection limit is unreported.")
    else: classify_cc(r,[E[2]])
    r["review"]["unresolved"] += ["Atmosphere, pressure, workup and purification are unreported in the supplied main text."]
    bundle["reactions"].append(r)


def add_scope_reaction(suffix, ev, substrate_mid, allenoate_mid, product_mid, product_label, yv, time_range,
                       temp="rt", scale_amount=0.2, scale_mass=None, isolated=True, category="substrate_scope"):
    r = reaction_base(suffix,ev,"specific_experiment",[category],suffix.replace("_"," "),"Ni(II)-catalyzed [3 + 2] cycloaddition of nitrone and allenoate","cascade synthesis of N-vinyl heterocycle")
    r["general_scheme_ids"]=[generic_id]
    a_amt=q(scale_amount,"mmol") if scale_amount is not None else None
    participant(r,"nitrone",substrate_mid,"reactant",ev,amount=a_amt,mass=q(scale_mass,"g") if scale_mass else None)
    participant(r,"allenoate",allenoate_mid,"reactant",ev,amount=q(0.4,"mmol") if scale_amount==0.2 else None)
    participant(r,"NiCl2",common["NiCl2"],"catalyst",ev,mol_percent=q(20,"mol%"),ref=substrate_mid)
    participant(r,"L1",common["L1"],"ligand",ev,mol_percent=q(24,"mol%"),ref=substrate_mid)
    participant(r,"MgSO4",common["MgSO4"],"additive",ev,mass=q(200,"mg"))
    participant(r,"DCM",common["DCM"],"solvent",ev,volume=q(2.0,"mL"))
    participant(r,"product",product_mid,"product",ev,product_role="desired product")
    if scale_amount is not None:
        r["reaction_scale"]=F(f"nitrone ({scale_amount} mmol)",{"reference_molecule_id":substrate_mid,"amount":q(scale_amount,"mmol")},ev)
    r["solvent_conditions"]["composition"]=F("DCM","DCM",ev); r["solvent_conditions"]["total_volume"]=F("2.0 mL",q(2.0,"mL"),ev)
    add_step(r,"reaction",ev,[p["participant_id"] for p in r["participants"]],"reaction under scheme footnote conditions",temp,time_range)
    r["procedure"]["general_procedure_reference"]=F("scheme footnote reaction conditions","scheme footnote reaction conditions",ev)
    yield_outcome(r,f"yield_{product_label}",product_mid,q(yv,"%"),ev,isolated=isolated)
    classify_cc(r,ev)
    r["review"]["unresolved"]=["Atmosphere, pressure, workup, purification and product-specific procedure are unavailable without explicitly associated SI.","Classification is partial because dense structures were not atom mapped."]
    return r


# Scheme 2 scope plus separate scale-up.
scheme2_yields={"3aa":83,"3ba":56,"3ca":59,"3da":74,"3ea":62,"3fa":83,"3ga":65,"3ha":52,"3ia":59,"3ja":47,"3ka":67,"3la":68,"3ma":49,"3na":81,"3oa":35,"3pa":79,"3qa":48,"3ra":19,"3sa":50,"3ta":24,"3ua":70,"3va":73,"3wa":84,"3xa":60}
scope_time={"minimum":16,"maximum":72,"unit":"h","raw_text":"16-72 h"}
for pl,yv in scheme2_yields.items():
    substrate=f"1{pl[1]}" if pl!="3aa" else "1a"
    scope_ev=[E[2],T[3]] if pl in {"3ja","3ka","3sa","3va","3wa","3xa"} else [E[2]]
    r=add_scope_reaction(f"scheme2_scope_{pl}",scope_ev,nitrones1[substrate],allenoates["2a"],products3[pl],pl,yv,scope_time,temp=q(60,"°C") if pl=="3xa" else "rt")
    if pl=="3xa":
        r["outcomes"][0]["value"]=conflict(["60% in Scheme 2","63% in prose"],[F("60%",q(60,"%"),[E[2]]),F("63%",q(63,"%"),[T[3]])],[E[2],T[3]])
        r["review"]["unresolved"].append("Scheme 2 prints 60% for 3xa; page-3 prose states 63%. Neither value is selected.")
    if pl=="3sa": selectivity_outcome(r,"rr_3sa",products3[pl],"1.5/1","rearrangement regioselectivity",[E[2],T[3]])
    if pl=="3va": selectivity_outcome(r,"rr_3va",products3[pl],"single isomer; high rearrangement regioselectivity","rearrangement regioselectivity",[T[3]])
    bundle["reactions"].append(r)

r=add_scope_reaction("scheme2_scaleup_3aa",[E[2],T[2]],m1a,allenoates["2a"],products3["3aa"],"3aa",67,scope_time,scale_amount=3.0,scale_mass=1.0,category="scale_up")
r["review"]["unresolved"].append("Scheme footnote gives 1a 1.0 g/3.0 mmol but does not state whether 2a, MgSO4, DCM and time were rescaled; printed standard-condition values are retained as inherited conditions.")
bundle["reactions"].append(r)

# Scheme 3 allenoate scope.
scheme3_yields={"3ab":79,"3ac":74,"3ad":57,"3ae":31,"3af":58,"3ag":55,"3ah":42,"3ai":69,"3aj":61,"3ak":54,"3al":25,"3am":18,"3an":35}
for pl,yv in scheme3_yields.items():
    al=f"2{pl[2]}"
    r=add_scope_reaction(f"scheme3_scope_{pl}",[E[3]],m1a,allenoates[al],products3[pl],pl,yv,scope_time,temp=q(60,"°C") if pl=="3ag" else "rt")
    if pl=="3an": selectivity_outcome(r,"dr_3an",products3[pl],"1/1","dr",[E[3],T[3]])
    bundle["reactions"].append(r)

# Scheme 4 N-vinyl nitrone scope.
scheme4_map={"5aa":("4a","2a",44),"5ba":("4b","2a",33),"5ca":("4c","2a",75),"5da":("4d","2a",61),"5ea":("4e","2a",49),"5fa":("4f","2a",39),"5ab":("4a","2b",38),"5af":("4a","2f",29),"5ag":("4a","2g",36)}
scope4_time={"minimum":16,"maximum":24,"unit":"h","raw_text":"16-24 h"}
for pl,(nl,al,yv) in scheme4_map.items():
    r=add_scope_reaction(f"scheme4_scope_{pl}",[E[3]],nitrones4[nl],allenoates[al],products5[pl],pl,yv,scope4_time,temp=q(60,"°C") if pl=="5ag" else "rt")
    bundle["reactions"].append(r)

# Partial real mechanistic study described in the main text; detailed experiment is only advertised in unavailable SI.
r=reaction_base("mechanistic_o18_label",[T[3],T[4]],"specific_experiment",["mechanistic_study"],"Mechanistic discussion on pages 3-4","18O-labeled experiment","oxygen-source labeling study")
participant(r,"1a",m1a,"reactant",[T[3]]); participant(r,"2a",allenoates["2a"],"reactant",[T[3]]); participant(r,"3aa",products3["3aa"],"product",[T[3]],product_role="labeled-product observation target")
add_step(r,"labeling_experiment",[T[3],T[4]],["1a","2a","3aa"],"18O-labeled experiment; detailed procedure in unassociated SI")
o=st.new_outcome("oxygen_origin","other",[products3["3aa"]],[T[3]]); o["value"]=F("the hydroxyl group in 3aa might be derived from water","authors state that the labeling result suggests the hydroxyl group might derive from water",[T[3]]); o["product_role"]=F("mechanistic observation target","mechanistic observation target",[T[3]]); r["outcomes"].append(o)
unresolved_class(r,[T[3],T[4]],"Label source, incorporation data, conditions, structures and atom mapping are unavailable in the supplied main text.")
r["review"]["unresolved"]=["The isotope source, enrichment, numerical incorporation, conditions, workup and analytical method are unavailable because SI was not supplied or associated."]
bundle["reactions"].append(r)

# Scheme 5 proposed mechanism (not an experimental proof record).
intermediates={x:add_molecule(f"mechanism_{x}",[E[3]],f"proposed intermediate {x}") for x in "ABCDE"}
complex_mid=add_molecule("NiCl2_L1_complex",[E[3]],"NiCl2-L1 complex")
r=reaction_base("scheme5_proposed_mechanism",[E[3],T[3],T[4]],"proposed_mechanism",None,"Scheme 5 and continuation on page 4","Proposed mechanism","author-proposed cascade mechanism")
participant(r,"1a",m1a,"reactant",[E[3]]); participant(r,"2a",allenoates["2a"],"reactant",[E[3]]); participant(r,"complex",complex_mid,"catalyst",[E[3]])
for x in "ABCDE": participant(r,x,intermediates[x],"unknown",[E[3]])
participant(r,"3aa",products3["3aa"],"product",[E[3]],product_role="proposed final product")
ops=[("cycloaddition","[3+2] cycloaddition to A"),("rearrangement","[3,3]-rearrangement to B"),("retro_mannich","retro-Mannich reaction/enolation to imine C"),("aza_michael","intramolecular aza-Michael (1,4-addition) to D"),("elimination","elimination/HCl trapping to E"),("closure","intramolecular 1,4-addition to 3aa")]
for sid,op in ops: add_step(r,sid,[E[3],T[3],T[4]],["1a","2a","complex","A","B","C","D","E","3aa"],op)
unresolved_class(r,[E[3],T[3],T[4]],"Scheme 5 is an author-proposed mechanism, not an atom-mapped experimental transformation record.")
r["notes"].append(F("Intermediate E: calculated m/z 727.2334; found 727.2336","HRMS trace shown for proposed intermediate E",[E[3]]))
r["review"]["unresolved"]=["Elementary-step atom mapping, transient-species identities and HCl/water mass balance are not established in the supplied pages."]
bundle["reactions"].append(r)

# Non-synthetic observations are retained in source_context, not converted to synthesis records.
bundle["source_context"]=[
    {"context_id":f"{DOC}:context:hrms_trace_E","document_id":DOC,"kind":"mechanistic_mass_spectrometry",
     "evidence_ids":[E[3],T[3]],"raw_value":{"assignment":"proposed intermediate E","calculated_m_z":727.2334,"found_m_z":727.2336},
     "note":"Standalone HRMS trace observation; not a synthetic experiment record."},
    {"context_id":f"{DOC}:context:si_notice","document_id":DOC,"kind":"supplementary_information_notice",
     "evidence_ids":[T[4]],"raw_value":"Detailed mechanistic studies, experimental procedures, characteristic data, and copies of NMR spectra of compounds 3 and 5 (PDF)",
     "note":"Availability notice only; no SI file was supplied or explicitly associated."}
]

bundle["molecules"]=list(molecules.values())
st.write(RUN/"submitted.json",bundle)

# Save graph-only stereo diagnostics separately; they are not proof of source identity.
diagnostics={}
for label,mid in {"1a":m1a,**allenoates}.items():
    smi=st.value(molecules[mid]["smiles"])
    if smi:
        diagnostics[label]=st.stereo_inspect(smi)
st.write(RUN/"raw_model/stereo_diagnostics.json",diagnostics)

summary={"document_id":DOC,"molecules":len(bundle["molecules"]),"reactions":len(bundle["reactions"]),
         "specific_experiments":sum(r["record_type"]=="specific_experiment" for r in bundle["reactions"]),
         "machine_status":"needs_review","semantic_accuracy_verified":False}
st.write(RUN/"raw_model/extraction_summary.json",summary)
print(json.dumps(summary,ensure_ascii=False,indent=2))
