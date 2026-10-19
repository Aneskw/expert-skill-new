#!/usr/bin/env python3
"""Build source-derived v3.1 records for doc_71ecfcbf02cd66de, PDF pages 1-6."""
from pathlib import Path
import json, sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DATASET = ROOT / "runs/current/dataset.json"
OUT = ROOT / "runs/current/submitted.json"
DOC = "doc_71ecfcbf02cd66de"

data = st.read(DATASET)
bundle = st.empty(data["documents"])

def ev(page, image=False):
    return f"{DOC}:page{page}:{'image' if image else 'text'}"

def f(raw=None, norm=None, pages=(), method="original_extraction", derivation=None,
      status=None, alternatives=()):
    evidence = []
    for p in pages:
        evidence.extend([ev(p), ev(p, True)] if p <= 6 else [ev(p)])
    return st.fact(raw, norm, evidence, method=method, derivation=derivation,
                   status=status, alternatives=alternatives)

def set_status(m, chirality="unresolved", geometry="unresolved", representation="not_encoded",
               pages=(1,), assignment=None, mixture=None, axial=None, double_bonds=None):
    stereo = m["stereochemistry"]
    stereo["chirality_status"] = f(chirality, chirality, pages)
    stereo["geometry_status"] = f(geometry, geometry, pages)
    stereo["representation_status"] = f(representation, representation, pages)
    if assignment is not None:
        stereo["assignment"] = f(assignment, assignment, pages)
    if mixture is not None:
        stereo["mixture"] = f(mixture, mixture, pages)
    if axial is not None:
        stereo["axial"] = f(axial, axial, pages)
    if double_bonds is not None:
        stereo["double_bonds"] = f(double_bonds, double_bonds, pages)

def molecule(mid, label, name=None, pages=(2,), smiles=None, chirality="unresolved",
             geometry="unresolved", representation=None, assignment=None, mixture=None,
             axial=None, double_bonds=None, note=None):
    m = st.new_molecule(f"{DOC}:{mid}", DOC, [ev(p) for p in pages])
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, pages)
    m["names"]["original"] = f(name or label, name or label, pages)
    m["names"]["english"] = f(name or label, name or label, pages)
    m["cas_numbers"] = f()
    if smiles:
        diag = st.stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = f(smiles, diag["canonical_isomeric_smiles"], pages,
                         method="structure_derivation",
                         derivation="Machine structure candidate checked with standard stereo; not source verification.")
        m["stereochemistry"]["isomeric_smiles"] = f(smiles, diag["canonical_isomeric_smiles"], pages,
            method="structure_derivation", derivation="Canonical isomeric SMILES from the submitted candidate.")
        m["stereochemistry"]["identity_key"] = f(None, diag["stereo_graph_key"], pages,
            method="structure_derivation", derivation="standard stereo graph key for the machine candidate.")
        rep = representation or ("partial" if any(x["specified"] == "Unspecified" for x in diag["potential_stereo"]) else "encoded")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = f()
        rep = representation or "not_encoded"
    set_status(m, chirality, geometry, rep, pages, assignment, mixture, axial, double_bonds)
    m["identity_note"] = f(note or "Source label/drawing retained; complete stereochemical structure is not machine-verified.",
                            note or "Source label/drawing retained; complete stereochemical structure is not machine-verified.", pages)
    m["external_properties"] = []
    bundle["molecules"].append(m)
    return m["molecule_id"]

mids = {}
def M(key, *args, **kwargs):
    if key not in mids:
        mids[key] = molecule(key, *args, **kwargs)
    return mids[key]

# Exact simple components checked with standard stereo. CAS numbers are intentionally absent.
M("tempo", "TEMPO", "2,2,6,6-tetramethylpiperidin-1-oxyl", (2,3,5,6),
  "CC1(C)CCCC(C)(C)N1[O]", chirality="achiral", geometry="not_applicable",
  assignment="No stereogenic element in the submitted machine graph; source radical drawing retained.")
M("k2co3", "K2CO3", "potassium carbonate", (1,2,3,4,6),
  "[K+].[K+].[O-]C([O-])=O", chirality="achiral", geometry="not_applicable")
M("cs2co3", "Cs2CO3", "cesium carbonate", (2,),
  "[Cs+].[Cs+].[O-]C([O-])=O", chirality="achiral", geometry="not_applicable")
M("net3", "NEt3", "triethylamine", (2,4,5), "CCN(CC)CC", chirality="achiral", geometry="not_applicable")
M("mecn", "MeCN", "acetonitrile", (2,3,5,6), "CC#N", chirality="achiral", geometry="not_applicable")
M("acetone", "acetone", "acetone", (2,4), "CC(C)=O", chirality="achiral", geometry="not_applicable")
M("thf", "THF", "tetrahydrofuran", (5,), "C1CCOC1", chirality="achiral", geometry="not_applicable")
M("dmf", "DMF", "dimethylformamide", (5,), "CN(C)C=O", chirality="achiral", geometry="not_applicable")
M("oxalyl_chloride", "oxalyl chloride", "oxalyl chloride", (5,), "O=C(Cl)C(=O)Cl", chirality="achiral", geometry="not_applicable")
M("nahco3", "NaHCO3", "sodium bicarbonate", (5,), "[Na+].O=C([O-])O", chirality="achiral", geometry="not_applicable")
M("nmethylhydroxylamine_hcl", "N-methylhydroxylamine hydrochloride", pages=(5,), chirality="achiral", geometry="not_applicable")
M("cdi", "1,1'-carbonyldiimidazole", pages=(5,), chirality="achiral", geometry="not_applicable")
M("nah", "NaH", "sodium hydride", (5,), "[Na+].[H-]", chirality="achiral", geometry="not_applicable")
M("dncbenzene", "1-chloro-2,4-dinitrobenzene", pages=(5,), chirality="achiral", geometry="not_applicable")

# Central reaction molecules. Allene axial configuration is not assigned by the source.
for lab in ["1a","1b","1c","1d","1e","1f","1g","1h","1i","1j","1k","1l"]:
    pages = (2,3,6) if lab in ("1a","1b","1c") else (3,)
    M(lab, lab, pages=pages, chirality="unspecified", geometry="not_applicable",
      axial={"element":"allene axis","configuration":"unassigned"},
      note="Allene substrate drawing retained. Axial configuration and any other stereochemical configuration are not reported.")

for lab in ["2a-H","2a","3a","4a","5a","3aa","3ab","3ac","3ad","3b","5b","3c","3d","2e","3f","3g","2h","3i","2j","2k","5k","2l","5l","3ab_crossover"]:
    pages = (2,3,6) if lab in ("2a-H","2a","3a","4a","5a","3b","5b","3c") else (3,)
    if lab.startswith("3") and lab not in ("3i",):
        M(lab, lab, pages=pages, chirality="mixture", geometry="not_applicable",
          mixture={"type":"diastereomeric mixture where a dr is reported; configurations unassigned"})
    elif lab in ("2a","2e","2h","3i","2j","2k","2l"):
        geom = "mixture" if lab == "2e" else ("specified" if lab in ("2a","3i","2j") else "unspecified")
        db = ({"configuration":"E for isolated 2a; Z/E identities are source-dependent"} if lab == "2a" else
              {"configuration":"Z as assigned in the source prose"} if lab in ("3i","2j") else None)
        M(lab, lab, pages=pages, chirality="unspecified", geometry=geom,
          double_bonds=db, note="TMP enol-ether product; alkene geometry is retained only where explicitly assigned by the source.")
    else:
        M(lab, lab, pages=pages, chirality="unspecified", geometry="unspecified")

# TEMPO derivatives used in Table 2.
for key,label in [("ho_tempo","4-hydroxy-TEMPO"),("meo_tempo","4-methoxy-TEMPO"),
                  ("bzo_tempo","4-benzoyloxy-TEMPO"),("oxo_tempo","4-oxo-TEMPO")]:
    M(key,label,pages=(3,5),chirality="unspecified",geometry="not_applicable",
      note="Source radical drawing/name retained; complete structure not machine-verified in this extraction.")

def participant(pid, mid, role, pages, eq=None, amount=None, mass=None, volume=None,
                molpct=None, stage="reaction", product_role=None, reference=None):
    p = st.new_participant(pid, mid, role, [ev(x) for x in pages], stage)
    p["product_role"] = f(product_role, product_role, pages) if product_role else f()
    vals = {"mass":mass,"amount":amount,"volume":volume,"equivalents":eq,"mol_percent":molpct}
    for k,v in vals.items():
        p["charges"][k] = f(str(v) if v is not None else None, v, pages)
    p["equivalence_reference"] = f(reference, reference, pages) if reference else f()
    return p

def labels_for(codes, pages):
    changes=[]; labels=[]
    mapping={"C-N":("C","N",1,5),"C-O":("C","O",1,3),"DB-C-O":("C","O",2,14),"DB-C-C":("C","C",2,13),"N-H":("N","H",1,29)}
    for i,code in enumerate(codes,1):
        a,b,order,num=mapping[code]
        ch={"change_id":f"b{i}","atoms":[{"atom_map":None,"element":a},{"atom_map":None,"element":b}],
            "before_order":0,"after_order":order,"change_type":"formed" if b != "H" else "hydrogen_gain"}
        changes.append(f(None,ch,pages,method="structure_derivation",derivation="Partial net bond interpretation from source drawings; atom mapping unresolved."))
        labels.append({"category_number":num,"code":code,"supporting_change_ids":[f"b{i}"]})
    return changes, labels

def reaction(suffix, pages, rtype="specific_experiment", categories=None, name=None,
             transformation=None, codes=(), unresolved=None):
    rid=f"{DOC}:{suffix}"
    r=st.new_reaction(rid,DOC,[ev(p) for p in pages])
    r["record_type"]=rtype
    r["experiment_categories"]=f(categories or [],categories or [],pages)
    r["location"]=f({"pages":list(pages)}, {"pages":list(pages)}, pages)
    r["reaction_name"]=f(name,name,pages) if name else f()
    r["transformation_type"]=f(transformation,transformation,pages) if transformation else f()
    if codes:
        changes,labs=labels_for(codes,pages)
        r["classification"]["bond_changes"]=changes
        r["classification"]["labels"]=f(None,labs,pages,method="structure_derivation",derivation="Partial multi-label classification from the depicted net transformation.")
        r["classification"]["coverage"]=f(None,"partial",pages,method="inference",derivation="No atom-mapped reaction SMILES; only visually supported net new-bond classes are represented.")
        r["classification"]["limitations"]=["Atom mapping and complete bond-change accounting are unresolved.","Hydrogen transfers and leaving-group bonds are not fully assigned."]
    else:
        r["classification"]["coverage"]=f(None,"unresolved",pages,method="inference",derivation="Insufficient source structure detail for safe classification.")
        r["classification"]["labels"]=f(None,[],pages,method="inference",derivation="No bond-formation label asserted.")
    r["review"]["unresolved"]=unresolved or ["Complete source-level structures and atom mapping require expert review."]
    bundle["reactions"].append(r)
    return r

def step(r, sid, pages, pids, operation, temp=None, time=None, atmosphere=None, stage="reaction"):
    s=st.new_step(sid,[ev(p) for p in pages],stage)
    s["participant_ids"]=list(pids)
    s["operation"]=f(operation,operation,pages)
    s["temperature"]=f(str(temp) if temp is not None else None,temp,pages)
    s["time"]=f(str(time) if time is not None else None,time,pages)
    s["atmosphere"]=f(atmosphere,atmosphere,pages) if atmosphere else f()
    s["pressure"]=f()
    r["steps"].append(s)

def outcome(r, oid, kind, targets, pages, value=None, basis=None, yield_type=None, scope=None, role="product"):
    o=st.new_outcome(oid,kind,targets,[ev(p) for p in pages])
    o["product_role"]=f(role,role,pages)
    o["value"]=f(str(value) if value is not None else None,value,pages)
    o["yield_type"]=f(yield_type,yield_type,pages) if yield_type else f()
    o["measurement_basis"]=f(basis,basis,pages) if basis else f()
    o["scope"]=f(scope,scope,pages) if scope else f()
    r["outcomes"].append(o)

def set_common_photo(r, substrate_key, product_keys, pages, tempo_key="tempo", tempo_eq=3, base_key="k2co3", base_eq=2,
                     solvent_key="mecn", concentration={"number":0.025,"unit":"M"}, scale=None):
    sub=M(substrate_key,substrate_key,pages=pages,chirality="unspecified",geometry="not_applicable")
    r["participants"].append(participant("substrate",sub,"reactant",pages,eq={"number":1,"unit":"equiv"},reference=sub))
    r["participants"].append(participant("tempo",mids[tempo_key],"reagent",pages,eq={"number":tempo_eq,"unit":"equiv"},reference=sub))
    if base_key:
        r["participants"].append(participant("base",mids[base_key],"reagent",pages,eq={"number":base_eq,"unit":"equiv"},reference=sub))
    r["participants"].append(participant("solvent",mids[solvent_key],"solvent",pages))
    for key in product_keys:
        r["participants"].append(participant("product_"+key.replace("-","_"),mids[key],"product",pages,product_role="reported product"))
    r["reaction_scale"]=f(str(scale) if scale else None,{"reference_molecule_id":sub,"amount":scale} if scale else None,pages)
    r["solvent_conditions"]["composition"]=f(solvent_key,solvent_key,pages)
    r["solvent_conditions"]["concentration"]=f(str(concentration),concentration,pages) if concentration else f()
    r["solvent_conditions"]["concentration_basis"]=f("substrate",sub,pages) if concentration else f()
    r["special_conditions"]["light_source"]=f("Kessil blue LED lamp","Kessil blue LED lamp",pages)
    r["special_conditions"]["wavelength"]=f("440 nm",{"number":440,"unit":"nm"},pages)
    r["special_conditions"]["reactor"]=f("25 mL Schlenk tube; lamp approximately 3 cm from vessel",{"vessel":"25 mL Schlenk tube","lamp_distance":"approximately 3 cm"},pages)
    r["procedure"]["general_procedure_reference"]=f("General Procedure for the Visible-Light-Mediated Amidation of Allenes and Optimization Studies", "page 6 general procedure", pages)
    r["procedure"]["addition_order"]=f("substrate and base; evacuate/purge; TEMPO; dry solvent; irradiate",["substrate","base","TEMPO","dry solvent","irradiation"],pages)
    step(r,"irradiation",pages,[p["participant_id"] for p in r["participants"] if p["stage"]=="reaction"],"irradiate until TLC showed complete substrate consumption",None,{"minimum":4,"maximum":5,"unit":"h","raw_text":"usually 4-5 h"},"N2")

# Scheme 1 prior work and the present generic transformation.
for suffix,name in [("lit_montgomery_2016","Montgomery sequential borylcupration of a monosubstituted allene"),
                    ("lit_zhai_2018","Zhai tandem allene trifunctionalization"),
                    ("lit_schomaker_2014","Schomaker aziridination/dihydroxylation/ring-opening trifunctionalization")]:
    r=reaction(suffix,(2,),rtype="literature_example",categories=[],name=name,transformation="prior-work allene trifunctionalization",
               unresolved=["Literature summary only; not a performed experiment of the current paper."])
    r["completeness"]="partial"

generic = reaction("general_allene_trifunctionalization",(1,2),rtype="generic_scheme",name="Allene trifunctionalization via amidyl radical cyclization and TEMPO trapping",transformation="visible-light radical cyclization/TEMPO trapping",codes=("C-N","C-O","DB-C-O"))
generic["participants"]=[participant("generic_substrate",mids["1a"],"reactant",(1,2)),participant("generic_tempo",mids["tempo"],"reagent",(1,2)),participant("generic_product",mids["3a"],"product",(1,2),product_role="desired trifunctionalized product")]

# Table 1: 12 independent optimization rows plus the separately described dark control.
table1={
1:("acetone",3,"k2co3",2,None,{"2a":14,"3a":36,"4a":20}),
2:("acetone",3,None,None,None,{"1a":37,"2a":8,"3a":19}),
3:("mecn",3,"k2co3",2,None,{"2a":4,"3a":44,"4a":14,"5a":14}),
4:("mecn",3,"k2co3",2,"no precautions",{"2a":5,"3a":42,"4a":11,"5a":13}),
5:("mecn",2,"k2co3",2,None,{"2a-H":7,"2a":7,"3a":33,"4a":5,"5a":11}),
6:("mecn",4,"k2co3",2,None,{"3a":42,"4a":14,"5a":7}),
7:("mecn",3,"k2co3",2,{"number":0.05,"unit":"M"},{"1a":16,"2a":11,"3a":29,"4a":18}),
8:("mecn",3,"k2co3",2,{"number":0.01,"unit":"M"},{"2a-H":8,"2a":8,"3a":39,"4a":5,"5a":13}),
9:("mecn",3,"k2co3",0.5,None,{"2a":7,"3a":38,"4a":3,"5a":6}),
10:("mecn",3,"k2co3",4,None,{"2a":8,"3a":36,"4a":18,"5a":5}),
11:("mecn",3,"cs2co3",2,None,{"2a":13,"3a":43,"4a":11}),
12:("mecn",3,"net3",2,None,{"2a-H":5,"2a":24,"3a":9}),
}
for n,(solv,teq,bkey,beq,mod,results) in table1.items():
    products=[k for k in results if k!="1a"]
    r=reaction(f"table1_entry_{n:02d}",(2,),categories=["optimization"],name="Optimization of reaction conditions",transformation="allene amidyl-radical cyclization/TEMPO trapping",codes=("C-N","C-O","DB-C-O"),unresolved=["Blank table cells are unreported, not zero.","Complete structures and all side-product bond changes remain unresolved."])
    conc=mod if isinstance(mod,dict) else {"number":0.025,"unit":"M"}
    set_common_photo(r,"1a",products,(2,),tempo_eq=teq,base_key=bkey,base_eq=beq or 0,solvent_key=solv,concentration=conc)
    r["general_scheme_ids"]=[generic["reaction_id"]]
    if mod=="no precautions": r["notes"].append(f("No precautions taken.","No precautions taken.",(2,)))
    for key,val in results.items():
        kind="recovery" if key=="1a" else "yield"
        outcome(r,f"{kind}_{key}",kind,[mids[key]],(2,),{"number":val,"unit":"%"},"1H NMR with 1,3,5-trimethoxybenzene internal standard",scope={"table":"Table 1","entry":n},role="recovered starting material" if kind=="recovery" else "product")
    if "3a" in results:
        outcome(r,"dr_3a","selectivity",[mids["3a"]],(2,),"1.3:1 to 1:1",None,scope={"metric":"dr","order":"unassigned diastereomers"})

dark=reaction("table1_dark_control",(2,),categories=["control","mechanistic_study"],name="Dark control",transformation="control without irradiation",unresolved=["Control duration and analytical method are not reported."])
dark["participants"]=[participant("substrate",mids["1a"],"reactant",(2,),eq={"number":1,"unit":"equiv"},reference=mids["1a"]),participant("tempo",mids["tempo"],"reagent",(2,),eq={"number":3,"unit":"equiv"},reference=mids["1a"]),participant("base",mids["k2co3"],"reagent",(2,),eq={"number":2,"unit":"equiv"},reference=mids["1a"]),participant("solvent",mids["acetone"],"solvent",(2,))]
step(dark,"dark_hold",(2,),[p["participant_id"] for p in dark["participants"]],"kept in the dark",None,None,"N2")
outcome(dark,"recovered_1a","other",[mids["1a"]],(2,),"only 1a was recovered",None,role="recovered starting material")
dark["classification"]["coverage"]=f(None,"not_applicable",(2,),method="inference",derivation="No conversion was observed; no product bond formation asserted.")

# Table 2 TEMPO derivative screen.
for n,(tk,pk,y,dr,impurity) in enumerate([
    ("ho_tempo","3aa",52,"1:1.4",None),("meo_tempo","3ab",50,"1:1",None),
    ("bzo_tempo","3ac",65,"1:1","approximately 15% inseparable impurity"),("oxo_tempo","3ad",43,"1:1.1",None)],1):
    r=reaction(f"table2_entry_{n}",(3,),categories=["optimization"],name="Impact of TEMPO derivatives on trifunctionalization",transformation="allene trifunctionalization",codes=("C-N","C-O","DB-C-O"),unresolved=["Yield measurement basis is not stated on page 3.","Absolute and relative configurations of the reported diastereomers are unassigned."])
    set_common_photo(r,"1a",[pk],(3,),tempo_key=tk)
    r["general_scheme_ids"]=[generic["reaction_id"]]
    outcome(r,"yield","yield",[mids[pk]],(3,),{"number":y,"unit":"%"},None,scope={"table":"Table 2","entry":n})
    r["outcomes"][-1]["measurement_basis"] = f("Following the general procedure; isolated after column chromatography", "isolated after column chromatography", (8,))
    outcome(r,"dr","selectivity",[mids[pk]],(3,),dr,None,scope={"metric":"dr","order":"unassigned diastereomers"})
    if impurity:r["notes"].append(f(impurity,impurity,(3,)))

# Table 3 scope experiments. The 3a scale-up is a distinct experiment.
scope_rows=[
 ("scope_1a_small","1a",["3a"],{"3a":39},"1:1.6","isolated after column chromatography",{"number":0.104,"unit":"mmol"},(3,6)),
 ("scope_1a_scaleup","1a",["2a","3a"],{"2a":2,"3a":41},"1:1.4","isolated after column chromatography",{"number":1.130,"unit":"mmol"},(3,6)),
 ("scope_1b","1b",["3b","5b"],{"3b":25,"5b":10},"1.2:1","isolated after column chromatography",{"number":0.115,"unit":"mmol"},(3,6)),
 ("scope_1c","1c",["3c"],{"3c":68},"1:1","isolated after column chromatography",{"number":0.110,"unit":"mmol"},(3,6)),
 ("scope_1d","1d",["3d"],{"3d":61},"1:1.8",None,None,(3,)),
 ("scope_1e","1e",["2e"],{"2e":36},"Z:E 10:1",None,None,(3,)),
 ("scope_1f","1f",["3f"],{"3f":54},"1:1.1",None,None,(3,)),
 ("scope_1g","1g",["3g"],{"3g":20},"1.2:1",None,None,(3,)),
 ("scope_1h","1h",["2h"],{"2h":45},None,None,None,(3,)),
 ("scope_1i","1i",["3i"],{"3i":15},None,"crude 1H NMR with 1,3,5-trimethoxybenzene internal standard",None,(3,)),
 ("scope_1j","1j",["2j"],{"2j":40},None,"isolated yield stated in prose",None,(3,)),
 ("scope_1k","1k",["2k","5k"],{"2k":30,"5k":43},None,None,None,(3,)),
 ("scope_1l","1l",["2l","5l"],{"2l":16,"5l":47},None,None,None,(3,)),
]
for suffix,sub,products,yields,sel,basis,scale,pages in scope_rows:
    codes=("C-N","C-O","DB-C-O") if any(p.startswith("3") for p in products) else ("C-N","C-O") if any(p.startswith("2") for p in products) else ("C-N","DB-C-O","DB-C-C")
    unresolved=["Complete structures and atom mapping remain unresolved."]
    if basis is None: unresolved.append("Yield measurement basis is not stated in the page-3 scope figure/prose.")
    r=reaction(suffix,pages,categories=["substrate_scope"] + (["scale_up"] if suffix.endswith("scaleup") else []),name="Effect of allene substitution on trifunctionalization",transformation="visible-light allene radical cyclization/TEMPO trapping",codes=codes,unresolved=unresolved)
    set_common_photo(r,sub,products,pages,scale=scale)
    r["general_scheme_ids"]=[generic["reaction_id"]]
    for pk,y in yields.items():
        outcome(r,"yield_"+pk,"yield",[mids[pk]],pages,{"number":y,"unit":"%"},basis,scope={"table":"Table 3","substrate":sub})
        if basis is None:
            r["outcomes"][-1]["measurement_basis"] = f("Following the general procedure; isolated after column chromatography", "isolated after column chromatography", (7,))
    if sel:
        target=mids[products[0]]
        metric="E/Z" if "E" in sel or "Z" in sel else "dr"
        outcome(r,"selectivity", "selectivity",[target],pages,sel,None,scope={"metric":metric,"order":"as reported; configurations otherwise unassigned"})
    if suffix=="scope_1e":r["notes"].append(f("No detectable 3e.","No detectable 3e.",(3,)))
    if suffix=="scope_1b":r["notes"].append(f("2b-H co-eluted with one diastereomer of 3b; isolated dr/yield differed from crude material.","2b-H co-elution affected isolated yield and dr.",(6,)))
    # Product-specific charges and isolated amounts printed on page 6.
    detailed={
      "scope_1a_small":({"number":33.5,"unit":"mg"}, {"number":0.104,"unit":"mmol"}, [("3a",{"number":12.6,"unit":"mg"},{"number":0.041,"unit":"mmol"})], {"number":4.5,"unit":"h"}),
      "scope_1a_scaleup":(None,{"number":1.130,"unit":"mmol"},[("2a",{"number":6.9,"unit":"mg"},{"number":0.023,"unit":"mmol"}),("3a",{"number":143.2,"unit":"mg"},{"number":0.461,"unit":"mmol"})],{"number":4,"unit":"h"}),
      "scope_1b":({"number":43.5,"unit":"mg"},{"number":0.115,"unit":"mmol"},[("3b",{"number":10.6,"unit":"mg"},{"number":0.029,"unit":"mmol"}),("5b",{"number":2.4,"unit":"mg"},{"number":0.011,"unit":"mmol"})],None),
      "scope_1c":({"number":38.6,"unit":"mg"},{"number":0.110,"unit":"mmol"},[("3c",{"number":25.4,"unit":"mg"},{"number":0.075,"unit":"mmol"})],None),
    }
    if suffix in detailed:
        smass,samt,pvals,ptime=detailed[suffix]
        subp=next(p for p in r["participants"] if p["participant_id"]=="substrate")
        if smass: subp["charges"]["mass"]=f(str(smass),smass,(6,))
        subp["charges"]["amount"]=f(str(samt),samt,(6,))
        if ptime: r["steps"][0]["time"]=f(str(ptime),ptime,(6,))
        for pk,pmass,pamt in pvals:
            outcome(r,"product_mass_"+pk,"product_amount",[mids[pk]],(6,),pmass,"isolated after column chromatography")
            outcome(r,"product_amount_"+pk,"product_amount",[mids[pk]],(6,),pamt,"isolated after column chromatography")

# Proposed mechanism and two actual control experiments on page 4.
mech=reaction("scheme2_proposed_mechanism",(3,4),rtype="proposed_mechanism",name="Proposed mechanisms for desired allene trifunctionalization product and major byproducts",transformation="proposed EDA/amidyl-radical cyclization and TEMPO trapping",unresolved=["Mechanistic arrows and intermediates are proposed, not isolated experiments."])

cross=reaction("control_3a_meotempo_crossover",(4,),categories=["control","mechanistic_study"],name="TEMPO crossover experiment",transformation="crossover test",unresolved=["Amounts, solvent, temperature, time, and analytical method are not reported on page 4."])
cross["participants"]=[participant("reactant_3a",mids["3a"],"reactant",(4,)),participant("meotempo",mids["meo_tempo"],"reagent",(4,)),participant("candidate_3ab",mids["3ab"],"product",(4,),product_role="crossover product not detected")]
step(cross,"treatment",(4,),["reactant_3a","meotempo"],"treated 3a with MeO-TEMPO")
outcome(cross,"no_3ab","no_reaction",[mids["3ab"]],(4,),"no 3ab was noted",None,role="not detected")
cross["classification"]["coverage"]=f(None,"not_applicable",(4,),method="inference",derivation="No crossover product observed; no bond formation asserted.")

reexp=reaction("control_3a_reexposure",(4,),categories=["control","mechanistic_study"],name="Re-exposure of 3a to reaction conditions",transformation="stability/control experiment",unresolved=["Exact duration and analytical method are not reported on page 4."])
reexp["participants"]=[participant("reactant_3a",mids["3a"],"reactant",(4,)),participant("tempo",mids["tempo"],"reagent",(4,)),participant("base",mids["k2co3"],"reagent",(4,)),participant("solvent",mids["mecn"],"solvent",(4,)),participant("candidate_5a",mids["5a"],"product",(4,),product_role="product not detected")]
step(reexp,"reexposure",(4,),[p["participant_id"] for p in reexp["participants"]],"re-exposed 3a to the reaction conditions")
outcome(reexp,"no_5a","no_reaction",[mids["5a"]],(4,),"did not result in formation of 5a",None,role="not detected")
reexp["classification"]["coverage"]=f(None,"not_applicable",(4,),method="inference",derivation="No 5a formation observed; no product bond formation asserted.")

# Starting-material preparation molecules and generic procedures.
for key,label in [("acid_b","corresponding homoallenic carboxylic acid for C-1b"),("acid_d","corresponding homoallenic carboxylic acid for C-1d"),("acid_k","corresponding homoallenic carboxylic acid for C-1k"),("acid_l","corresponding homoallenic carboxylic acid for C-1l"),("alcohol_d1","allenic alcohol D-1")]:
    M(key,label,pages=(5,),chirality="unspecified",geometry="not_applicable",note="Source precursor identity is indicated generically; complete structure is not stated in the text block.")
for key,label in [("c1b","C-1b"),("c1d","C-1d"),("c1k","C-1k"),("c1l","C-1l"),("c1m","C-1m"),("1m","1m")]:
    M(key,label,pages=(5,6) if key in ("c1m","1m") else (5,6),chirality="unspecified",geometry="not_applicable",
      axial={"element":"allene axis","configuration":"unassigned"},note="Allenic precursor; axial configuration is not reported.")

gp1=reaction("general_hydroxamic_acid_oxalyl_chloride",(5,),rtype="generic_scheme",name="General Procedure for the Preparation of N-Methyl Hydroxamic Acids",transformation="carboxylic acid to N-methyl hydroxamic acid",codes=("C-N",))
gp1["participants"]=[participant("acid",mids["acid_b"],"reactant",(5,),eq={"number":1,"unit":"equiv"},reference=mids["acid_b"]),participant("thf",mids["thf"],"solvent",(5,)),participant("dmf",mids["dmf"],"additive",(5,)),participant("oxalyl_chloride",mids["oxalyl_chloride"],"reagent",(5,),eq={"number":1.3,"unit":"equiv"},reference=mids["acid_b"]),participant("nmha_hcl",mids["nmethylhydroxylamine_hcl"],"reagent",(5,),eq={"number":1.1,"unit":"equiv"},reference=mids["acid_b"]),participant("nahco3",mids["nahco3"],"reagent",(5,),eq={"number":2.2,"unit":"equiv"},reference=mids["acid_b"])]
step(gp1,"acid_chloride",(5,),["acid","thf","dmf","oxalyl_chloride"],"dropwise oxalyl chloride addition; stir until complete; concentrate",None,None,"N2")
step(gp1,"hydroxamate_addition",(5,),["nmha_hcl","nahco3"],"cool hydroxylamine/NaHCO3 suspension, then add acid chloride solution dropwise",{"number":0,"unit":"°C"},{"number":15,"unit":"min"},"N2")
step(gp1,"warm_hold",(5,),["acid","nmha_hcl","nahco3"],"warm to room temperature and stir overnight","room temperature","overnight","N2")

gp2=reaction("general_hydroxamic_acid_cdi",(5,),rtype="generic_scheme",name="General Procedure for the Preparation of N-Methyl Hydroxamic Acids Using 1,1'-Carbonyldiimidazole",transformation="carboxylic acid to N-methyl hydroxamic acid",codes=("C-N",))
gp2["participants"]=[participant("acid",mids["acid_l"],"reactant",(5,),eq={"number":1,"unit":"equiv"},reference=mids["acid_l"]),participant("dmf",mids["dmf"],"solvent",(5,),volume={"number":15,"unit":"mL"}),participant("cdi",mids["cdi"],"reagent",(5,),eq={"number":2,"unit":"equiv"},reference=mids["acid_l"]),participant("nmha_hcl",mids["nmethylhydroxylamine_hcl"],"reagent",(5,),eq={"number":3,"unit":"equiv"},reference=mids["acid_l"])]
step(gp2,"activation",(5,),["acid","dmf","cdi"],"add CDI in one portion and stir",None,{"number":30,"unit":"min"})
step(gp2,"reaction",(5,),["acid","cdi","nmha_hcl"],"add N-methylhydroxylamine hydrochloride and stir","room temperature",{"number":12,"unit":"h"})

gp3=reaction("general_nmethylaryloxy_carbamate",(5,),rtype="generic_scheme",name="General Procedure for the Preparation of N-Methylaryloxy Carbamate",transformation="allenic alcohol to N-methyl hydroxycarbamate",codes=("C-O","C-N"))
gp3["participants"]=[participant("alcohol",mids["alcohol_d1"],"reactant",(5,),eq={"number":1,"unit":"equiv"},reference=mids["alcohol_d1"]),participant("thf",mids["thf"],"solvent",(5,)),participant("cdi",mids["cdi"],"reagent",(5,)),participant("nmha_hcl",mids["nmethylhydroxylamine_hcl"],"reagent",(5,),eq={"number":1.2,"unit":"equiv"},reference=mids["alcohol_d1"]),participant("net3",mids["net3"],"reagent",(5,),eq={"number":1.2,"unit":"equiv"},reference=mids["alcohol_d1"])]
step(gp3,"activation",(5,),["alcohol","thf","cdi"],"stir until complete consumption by TLC",None,None,"N2")
step(gp3,"reaction",(5,),["nmha_hcl","net3"],"add N-methylhydroxylamine and Et3N; stir","room temperature",{"number":24,"unit":"h"},"N2")

gp4=reaction("general_aryloxy_amide",(5,6),rtype="generic_scheme",name="General Procedure for the Preparation of N-Methylaryloxy Amides",transformation="O-arylation of N-methyl hydroxamic acids",codes=("C-O",))
gp4["participants"]=[participant("hydroxamic_acid",mids["c1b"],"reactant",(5,6),eq={"number":1,"unit":"equiv"},reference=mids["c1b"]),participant("thf",mids["thf"],"solvent",(5,6)),participant("nah",mids["nah"],"reagent",(5,6),eq={"number":1.1,"unit":"equiv"},reference=mids["c1b"]),participant("aryl_chloride",mids["dncbenzene"],"reagent",(5,6),eq={"number":1.1,"unit":"equiv"},reference=mids["c1b"])]
step(gp4,"cool",(5,6),["hydroxamic_acid","thf"],"cool and stir",{"number":0,"unit":"°C"},{"number":15,"unit":"min"},"N2")
step(gp4,"deprotonation",(5,6),["hydroxamic_acid","nah"],"add NaH and stir",{"number":0,"unit":"°C"},{"number":1,"unit":"h"},"N2")
step(gp4,"arylation",(5,6),["aryl_chloride"],"add aryl chloride in three portions; warm to room temperature and stir overnight","room temperature","overnight","N2")

# Ten concrete precursor preparations reported on pages 5-6.
hydrox_specs=[
 ("prep_c1b","acid_b","c1b",{"number":289.7,"unit":"mg"},{"number":1.371,"unit":"mmol"},{"number":500.0,"unit":"mg"},{"number":2.74,"unit":"mmol"},50,gp1),
 ("prep_c1d","acid_d","c1d",{"number":360,"unit":"mg"},{"number":2.350,"unit":"mmol"},{"number":311.8,"unit":"mg"},{"number":1.580,"unit":"mmol"},67,gp1),
 ("prep_c1k","acid_k","c1k",{"number":180.6,"unit":"mg"},{"number":2.162,"unit":"mmol"},{"number":144.2,"unit":"mg"},{"number":1.663,"unit":"mmol"},51,gp1),
 ("prep_c1l","acid_l","c1l",{"number":400.0,"unit":"mg"},{"number":2.410,"unit":"mmol"},{"number":328.2,"unit":"mg"},{"number":1.681,"unit":"mmol"},70,gp2),
 ("prep_c1m","alcohol_d1","c1m",{"number":614.8,"unit":"mg"},{"number":7.362,"unit":"mmol"},{"number":172,"unit":"mg"},{"number":2.945,"unit":"mmol"},37,gp3),
]
for suffix,src,prod,smass,samt,pmass,pamt,y,gp in hydrox_specs:
    r=reaction(suffix,(5,),categories=["other"],name="Preparation of "+prod,transformation="starting-material preparation",codes=("C-N",),unresolved=["Printed mass/amount values are preserved without stoichiometric correction.","Yield is reported over two steps; intermediate acid chloride or activated intermediate was not isolated."])
    r["general_scheme_ids"]=[gp["reaction_id"]]
    r["participants"]=[participant("reactant",mids[src],"reactant",(5,),mass=smass,amount=samt),participant("product",mids[prod],"product",(5,),product_role="isolated precursor")]
    r["reaction_scale"]=f(str(samt),{"reference_molecule_id":mids[src],"amount":samt},(5,))
    if gp is gp1:
        r["participants"] += [participant("thf",mids["thf"],"solvent",(5,)),participant("dmf",mids["dmf"],"additive",(5,)),participant("oxalyl_chloride",mids["oxalyl_chloride"],"reagent",(5,),eq={"number":1.3,"unit":"equiv"},reference=mids[src]),participant("nmha_hcl",mids["nmethylhydroxylamine_hcl"],"reagent",(5,),eq={"number":1.1,"unit":"equiv"},reference=mids[src]),participant("nahco3",mids["nahco3"],"reagent",(5,),eq={"number":2.2,"unit":"equiv"},reference=mids[src])]
        step(r,"acid_chloride",(5,),["reactant"],"convert corresponding acid to acid chloride with oxalyl chloride/DMF in anhydrous THF; concentrate",None,None,"N2")
        step(r,"hydroxamate_addition",(5,),["reactant"],"add acid chloride solution dropwise to N-methylhydroxylamine/NaHCO3 after the cold pre-hold",{"number":0,"unit":"°C"},{"number":15,"unit":"min"},"N2")
        step(r,"warm_hold",(5,),["reactant"],"warm to room temperature and stir","room temperature","overnight","N2")
    elif gp is gp2:
        r["participants"] += [participant("dmf",mids["dmf"],"solvent",(5,),volume={"number":15,"unit":"mL"}),participant("cdi",mids["cdi"],"reagent",(5,),eq={"number":2,"unit":"equiv"},reference=mids[src]),participant("nmha_hcl",mids["nmethylhydroxylamine_hcl"],"reagent",(5,),eq={"number":3,"unit":"equiv"},reference=mids[src])]
        step(r,"activation",(5,),["reactant"],"activate acid with CDI in DMF",None,{"number":30,"unit":"min"})
        step(r,"reaction",(5,),["reactant"],"add N-methylhydroxylamine hydrochloride and stir","room temperature",{"number":12,"unit":"h"})
    else:
        r["participants"] += [participant("thf",mids["thf"],"solvent",(5,)),participant("cdi",mids["cdi"],"reagent",(5,)),participant("nmha_hcl",mids["nmethylhydroxylamine_hcl"],"reagent",(5,),eq={"number":1.2,"unit":"equiv"},reference=mids[src]),participant("net3",mids["net3"],"reagent",(5,),eq={"number":1.2,"unit":"equiv"},reference=mids[src])]
        step(r,"activation",(5,),["reactant"],"activate allenic alcohol with CDI in THF until complete consumption by TLC",None,None,"N2")
        step(r,"reaction",(5,),["reactant"],"add N-methylhydroxylamine and Et3N and stir","room temperature",{"number":24,"unit":"h"},"N2")
    outcome(r,"product_mass","product_amount",[mids[prod]],(5,),pmass,"isolated after column chromatography")
    outcome(r,"product_amount","product_amount",[mids[prod]],(5,),pamt,"isolated after column chromatography")
    outcome(r,"yield","yield",[mids[prod]],(5,),{"number":y,"unit":"%"},"isolated after column chromatography",yield_type="over two steps",scope={"start":mids[src],"end":mids[prod]})

ary_specs=[
 ("prep_1b","c1b","1b",{"number":282.4,"unit":"mg"},{"number":1.336,"unit":"mmol"},{"number":468.8,"unit":"mg"},{"number":1.242,"unit":"mmol"},93),
 ("prep_1d","c1d","1d",{"number":286.1,"unit":"mg"},{"number":1.450,"unit":"mmol"},{"number":505.5,"unit":"mg"},{"number":1.391,"unit":"mmol"},96),
 ("prep_1k","c1k","1k",{"number":129.7,"unit":"mg"},{"number":0.766,"unit":"mmol"},{"number":165.1,"unit":"mg"},{"number":492.4,"unit":"mmol"},64),
 ("prep_1l","c1l","1l",{"number":200.0,"unit":"mg"},{"number":1.02,"unit":"mmol"},{"number":341.0,"unit":"mg"},{"number":0.944,"unit":"mmol"},92),
 ("prep_1m","c1m","1m",{"number":155.7,"unit":"mg"},{"number":990.6,"unit":"mmol"},{"number":271.5,"unit":"mg"},{"number":0.840,"unit":"mmol"},85),
]
for suffix,src,prod,smass,samt,pmass,pamt,y in ary_specs:
    r=reaction(suffix,(6,),categories=["other"],name="Preparation of "+prod,transformation="O-arylation of N-methyl hydroxamic acid/carbamate",codes=("C-O",),unresolved=["Printed amount values are preserved as printed even where chemically implausible; no silent unit correction was made."])
    r["general_scheme_ids"]=[gp4["reaction_id"]]
    r["participants"]=[participant("reactant",mids[src],"reactant",(6,),mass=smass,amount=samt),participant("nah",mids["nah"],"reagent",(6,),eq={"number":1.1,"unit":"equiv"},reference=mids[src]),participant("aryl_chloride",mids["dncbenzene"],"reagent",(6,),eq={"number":1.1,"unit":"equiv"},reference=mids[src]),participant("thf",mids["thf"],"solvent",(6,)),participant("product",mids[prod],"product",(6,),product_role="isolated allene precursor")]
    step(r,"cool",(6,),["reactant","thf"],"cool and stir",{"number":0,"unit":"°C"},{"number":15,"unit":"min"},"N2")
    step(r,"deprotonation",(6,),["reactant","nah"],"add NaH and stir",{"number":0,"unit":"°C"},{"number":1,"unit":"h"},"N2")
    step(r,"arylation",(6,),["aryl_chloride"],"add in three portions; warm and stir","room temperature","overnight","N2")
    outcome(r,"product_mass","product_amount",[mids[prod]],(6,),pmass,"isolated after column chromatography")
    outcome(r,"product_amount","product_amount",[mids[prod]],(6,),pamt,"isolated after column chromatography")
    outcome(r,"yield","yield",[mids[prod]],(6,),{"number":y,"unit":"%"},"isolated after column chromatography")

# General visible-light procedure as its own generic record.
gphoto=reaction("general_visible_light_amidation",(6,),rtype="generic_scheme",name="General Procedure for the Visible-Light-Mediated Amidation of Allenes and Optimization Studies",transformation="visible-light allene cyclization/TEMPO trapping",codes=("C-N","C-O","DB-C-O"))
set_common_photo(gphoto,"1a",["3a"],(6,))
gphoto["notes"].append(f("After reaction, decant/concentrate; add 1,3,5-trimethoxybenzene for crude 1H NMR; purify by column chromatography.","Workup/analysis sequence as printed.",(6,)))

# Source context: characterization and non-experimental discussion are retained, not promoted to reactions.
bundle["source_context"]=[
 {"context_id":f"{DOC}:page1_background","document_id":DOC,"evidence_ids":[ev(1),ev(1,True)],"kind":"article_background","note":"Abstract/introduction and graphical abstract; present-work generic transformation is linked, not a separate experiment."},
 {"context_id":f"{DOC}:page4_mechanistic_discussion","document_id":DOC,"evidence_ids":[ev(4),ev(4,True)],"kind":"mechanistic_discussion","note":"Author-proposed pathways and reactivity summary; only explicit controls were modeled as experiments."},
 {"context_id":f"{DOC}:page5_characterization","document_id":DOC,"evidence_ids":[ev(5),ev(5,True)],"kind":"characterization","note":"NMR/HRMS peak lists and general analytical methods; not independent syntheses."},
 {"context_id":f"{DOC}:page6_characterization","document_id":DOC,"evidence_ids":[ev(6),ev(6,True)],"kind":"characterization","note":"NMR/HRMS peak lists for precursors/products; not independent syntheses."},
 {"context_id":f"{DOC}:si_notice","document_id":DOC,"evidence_ids":[ev(1)],"kind":"source_availability","note":"Article indicates Supporting Information exists online, but no SI file was explicitly provided or associated in this workspace."}
]

st.write(OUT,bundle)
print(json.dumps({"molecules":len(bundle["molecules"]),"reactions":len(bundle["reactions"]),"out":str(OUT)},ensure_ascii=False))
