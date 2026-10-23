#!/usr/bin/env python3
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
TEMPLATE = RUN / "raw_model/submitted_template_pages_0001_0005.json"
OUT = RUN / "submitted.json"
DOC = "doc_3d06842c48e86b8b"


def ev(page):
    return [f"{DOC}:page{page}:text", f"{DOC}:page{page}:image"]


def fx(raw, normalized, evidence, method="original_extraction", derivation=None,
       status=None):
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation,
                   status=status)


def derived(normalized, evidence, why):
    return fx(None, normalized, evidence, "structure_derivation", why)


def inferred(normalized, evidence, why):
    return fx(None, normalized, evidence, "inference", why)


def missing(evidence=(), state="unreported"):
    return st.fact(evidence=evidence, status=state)


bundle = json.loads(TEMPLATE.read_text())
molecules = {}


def add_molecule(mid, label, page, *, name=None, smiles=None, note=None,
                 chirality="unresolved", geometry="unresolved",
                 representation="unresolved", axial=None, centers=None,
                 double_bonds=None, mixture=None, assignment=None,
                 template=None):
    evidence = ev(page)
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = fx(label, label, evidence)
    m["names"]["original"] = fx(name or label, name or label, evidence)
    m["names"]["english"] = missing(evidence)
    m["names"]["abbreviations"] = fx(label, [label], evidence)
    m["names"]["aliases"] = missing(evidence)
    m["cas_numbers"] = missing(evidence)
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = derived(smiles, evidence,
            f"Candidate connectivity transcribed from the source drawing for {label}; parsing does not verify the drawing.")
        m["generic_template"] = None
    elif template:
        m["representation_type"] = "generic_template"
        m["smiles"] = missing(evidence)
        variables = []
        for var_label, definition in template.get("variables", []):
            variables.append({
                "label": var_label,
                "attachment_points": fx("source drawing", "source-labeled attachment point", evidence),
                "definition": fx(definition, definition, evidence),
                "allowed_values": fx(definition, [definition], evidence),
                "constraints": missing(evidence),
            })
        m["generic_template"] = {
            "fixed_structure": fx(template["fixed"], template["fixed"], evidence),
            "variables": variables,
            "instance_mapping": [],
            "reason": fx("generic/R-group drawing", template.get("reason", "Source uses a scaffold drawing with labeled substituent values; no unique complete SMILES is asserted."), evidence),
        }
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing(evidence)
        m["generic_template"] = None
    s = m["stereochemistry"]
    s["atom_centers"] = fx(centers, centers, evidence) if centers is not None else missing(evidence)
    s["double_bonds"] = fx(double_bonds, double_bonds, evidence) if double_bonds is not None else missing(evidence)
    s["axial"] = fx(axial, axial, evidence) if axial is not None else missing(evidence)
    s["mixture"] = fx(mixture, mixture, evidence) if mixture is not None else missing(evidence)
    s["assignment"] = fx(assignment, assignment, evidence) if assignment is not None else missing(evidence)
    for key, value in (("chirality_status", chirality), ("geometry_status", geometry),
                       ("representation_status", representation)):
        s[key] = inferred(value, evidence,
            f"Status records what is and is not specified in the source drawing/text for {label}; it does not infer racemic or achiral identity from silence.")
    s["isomeric_smiles"] = missing(evidence)
    s["identity_key"] = missing(evidence)
    m["identity_note"] = fx(note or f"Source identity {label}; complete structure not independently verified.",
                            note or f"Source identity {label}; complete structure not independently verified.", evidence)
    m["external_properties"] = []
    molecules[mid] = m
    return mid


def simple(mid, label, smiles, page=2, name=None):
    return add_molecule(mid, label, page, name=name or label, smiles=smiles,
                        chirality="achiral", geometry="not_applicable",
                        representation="not_applicable",
                        note=f"Candidate structure for source reagent/solvent {label}; no CAS was supplied.")


# Shared source structures and reagents.
add_molecule("m_1a", "1a", 2, name="N-(2,3-allenyl)glycine ethyl ester 1a",
             smiles="C=C=CCN(C)CC(=O)OCC", chirality="achiral",
             geometry="not_applicable", representation="not_applicable",
             note="Source drawing shows a terminal allene and no stereogenic allene axis for 1a.")
add_molecule("m_2a", "2a", 2, name="2-(trimethylsilyl)phenyl triflate 2a",
             smiles="C[Si](C)(C)c1ccccc1OS(=O)(=O)C(F)(F)F",
             chirality="achiral", geometry="not_applicable", representation="not_applicable")
add_molecule("m_3a", "3a", 2, name="2-vinylallyamine 3a",
             smiles="CCOC(=O)C(N(C)c1ccccc1)C(=C)C=C",
             chirality="unspecified", geometry="not_applicable", representation="not_encoded",
             centers=[{"site": "alpha carbon", "configuration": "not assigned in source"}],
             note="Candidate connectivity from the table drawing; the newly formed stereogenic center is not assigned.")
add_molecule("m_3ap", "3a′", 2, name="1-amino-1,3-diene 3a′",
             template={"fixed": "isomerized 1-amino-1,3-diene scaffold shown for 3a′",
                       "variables": [("alkene geometry", "Z and E forms, assignments not structurally encoded")],
                       "reason": "The source reports a Z/E mixture and the visual transcription is retained without inventing a single isomeric SMILES."},
             chirality="unspecified", geometry="mixture", representation="not_encoded",
             double_bonds={"configuration": "Z/E mixture", "absolute drawing-to-isomer mapping": "not encoded"},
             mixture={"type": "geometric-isomer mixture", "composition": "experiment-specific"})

simple("m_KF", "KF", "[K+].[F-]", name="potassium fluoride")
simple("m_CsF", "CsF", "[Cs+].[F-]", name="cesium fluoride")
simple("m_ZnF2", "ZnF2", "[Zn+2].[F-].[F-]", name="zinc fluoride")
simple("m_NaF", "NaF", "[Na+].[F-]", name="sodium fluoride")
simple("m_NH4F", "NH4F", "[NH4+].[F-]", name="ammonium fluoride")
simple("m_TBAF", "TBAF", "CCCC[N+](CCCC)(CCCC)CCCC.[F-]", name="tetrabutylammonium fluoride")
simple("m_18c6", "18-crown-6", "O1CCOCCOCCOCCOCC1", name="18-crown-6")
simple("m_MeCN", "MeCN", "CC#N", name="acetonitrile")
simple("m_THF", "THF", "C1CCOC1", name="tetrahydrofuran")
simple("m_DMF", "DMF", "CN(C)C=O", name="N,N-dimethylformamide")
simple("m_dioxane", "1,4-dioxane", "O1CCOCC1")
simple("m_toluene", "toluene", "Cc1ccccc1")
simple("m_DCE", "DCE", "ClCCCl", name="1,2-dichloroethane")
simple("m_DEE", "EtOCH2CH2OEt", "CCOCCOCC", name="1,2-diethoxyethane")


def generic_label(mid, label, page, description, *, stereo=None):
    stereo = stereo or {}
    return add_molecule(mid, label, page, name=description,
        template={"fixed": description, "variables": [("source substituents", "values are those printed for this row/arrow")],
                  "reason": "The complete source drawing is retained as a labeled scaffold; a unique full SMILES was not asserted."},
        chirality=stereo.get("chirality", "unspecified"),
        geometry=stereo.get("geometry", "unspecified"),
        representation=stereo.get("representation", "not_encoded"),
        axial=stereo.get("axial"), centers=stereo.get("centers"),
        double_bonds=stereo.get("double_bonds"), mixture=stereo.get("mixture"),
        assignment=stereo.get("assignment"), note=description)


# Table 2 identities. 1a and 3a are already concrete candidates.
t2_rows = {
  1:("1a","3a",76,None,None), 2:("1b","3b",65,None,None),
  3:("1c","3c",73,None,None), 4:("1d","3d",48,None,None),
  5:("1e","3e",0,None,None), 6:("1f","3f′",49,">99 : 1 Z/E",None),
  7:("1g","3g′",44,">99 : 1 Z/E",None), 8:("1h","3h′",91,"25 : 75 Z/E",None),
  9:("1i","3i",73,None,None), 10:("1j","3j",65,None,None),
  11:("1k","3k",0,None,None), 12:("1l","3l",90,"10 : 90 Z/E",None),
  13:("1m","3m",63,None,None), 14:("1n","3n",69,None,None),
  15:("1o","3o",85,"<1 : 99 Z/E","70% ee"),
  16:("1p","3p",74,"<1 : 99 Z/E","81% ee"),
  17:("1q","3q",39,None,"90% ee"),
}
substrate_desc = {
  "1b":"tert-butyl ester analogue 1b (EWG = CO2tBu)",
  "1c":"benzyl ester analogue 1c (EWG = CO2Bn)",
  "1d":"p-nitrophenyl-substituted allenylamine 1d",
  "1e":"allenylamine amide 1e (EWG = CONH2)",
  "1f":"allenylamine ketone 1f (EWG = COMe)",
  "1g":"allenylamine ketone 1g (EWG = COPh)",
  "1h":"allenylamine nitrile 1h (EWG = CN)",
  "1i":"N-alkyl allenylglycine ester 1i (R1 = n-butyl)",
  "1j":"N-allenyl-substituted allenylglycine ester 1j",
  "1k":"N-phenyl analogue 1k",
  "1l":"substituted axially unspecified allene 1l (R2 = n-C7H15, R3 = H)",
  "1m":"trisubstituted allene 1m (R2 = R3 = Et)",
  "1n":"cyclically substituted allene 1n (R2,R3 = (CH2)5)",
  "1o":"enantioenriched axially chiral allene 1o (R2 = Cy, >99% ee)",
  "1p":"enantioenriched axially chiral allene 1p (R2 = Ph, 98% ee)",
  "1q":"L-proline-derived allenylamine 1q (>99% ee)",
}
for entry,(sub,prod,_,zr,ee) in t2_rows.items():
    if sub != "1a":
        sp = {}
        if sub in ("1o","1p"):
            sp = {"chirality":"partially_specified","geometry":"not_applicable","representation":"not_encoded",
                  "axial":{"axis":"allene","absolute configuration":"not reported in table"},
                  "assignment":{"ee":">99%" if sub=="1o" else "98%","absolute configuration":"unassigned"}}
        elif sub == "1q":
            sp = {"chirality":"partially_specified","geometry":"not_applicable","representation":"not_encoded",
                  "centers":{"source":"L-proline-derived drawing","absolute descriptor":"not printed"},
                  "assignment":{"ee":">99%","absolute configuration":"unassigned"}}
        generic_label(f"m_{sub}", sub, 3, substrate_desc[sub], stereo=sp)
    if prod != "3a":
        ps = {}
        if zr:
            ps.update(geometry="mixture", representation="not_encoded",
                      double_bonds={"reported order":"Z/E","ratio":zr.replace(" Z/E","")},
                      mixture={"type":"Z/E mixture","ratio":zr.replace(" Z/E","")})
        if ee:
            ps.update(chirality="partially_specified", representation="not_encoded",
                      centers={"stereogenic element":"product center shown/implicit; absolute configuration not printed"},
                      assignment={"ee":ee,"absolute configuration":"unassigned"})
        generic_label(f"m_{prod.replace('′','p')}", prod, 3,
                      f"source product {prod} from Table 2 entry {entry}", stereo=ps)

add_molecule("m_side_1e", "2-(methyl(phenyl)amino)acetamide", 5,
             name="2-(methyl(phenyl)amino)acetamide",
             smiles="CN(c1ccccc1)CC(N)=O", chirality="achiral",
             geometry="not_applicable", representation="not_applicable",
             note="Footnote 13 identifies this formal-substitution product from 1e; candidate connectivity follows the printed name.")

# Scheme 3 and Table 3 source identities.
for label,desc in {
 "1r":"N-(2,3-allenyl)-N-allylglycine ester 1r", "3r":"2-vinylallyamine 3r",
 "4a":"homoallylamine 4a", "1s":"N-(2,3-allenyl)-N-propargylglycine ester 1s",
 "3s":"2-vinylallyamine 3s", "5a":"propargyl-migration product 5a",
 "1t":"N-(2,3-allenyl)-N-benzylglycine ester 1t", "3t":"2-vinylallyamine 3t",
 "6a":"benzyl [1,2]-rearrangement product 6a"}.items():
    generic_label(f"m_{label}", label, 3, desc)

for label,desc in {
 "2b":"6-methyl aryne precursor 2b", "2c":"4-methyl aryne precursor 2c",
 "2d":"4,5-dimethyl aryne precursor 2d", "2e":"4,5-dimethoxy aryne precursor 2e",
 "2f":"naphthyne precursor 2f", "2g":"2,3-pyridyne precursor 2g",
 "3u":"3-methyl product 3u", "3v":"4-methyl regioisomer 3v",
 "3u4":"alternate 4-methyl regioisomer labelled 3u in Table 3",
 "3w":"3,4-dimethyl product 3w", "3x":"3,4-dimethoxy product 3x",
 "3y":"naphthyl product 3y", "3z":"pyridyl target 3z"}.items():
    generic_label(f"m_{label}", label, 3, desc)

# Generic Scheme 2 identities.
for label,desc in {
 "g1":"generic tertiary 2,3-allenylamine with R/R′/EWG substituents",
 "g2":"generic aryne", "g3":"generic 2-vinylallyamine rearrangement product",
 "g4":"generic base-isomerized 1-amino-1,3-diene"}.items():
    generic_label(f"m_{label}", label, 2, desc)


def charge(p, key, raw, normalized, evidence, method="original_extraction"):
    p["charges"][key] = fx(raw, normalized, evidence, method=method)


def participant(pid, mid, role, evidence, *, amount=None, volume=None,
                product_role=None, method="original_extraction"):
    p = st.new_participant(pid, mid, role, evidence)
    if amount is not None:
        charge(p, "amount", f"{amount} mmol", {"number": amount, "unit":"mmol"}, evidence, method)
    if volume is not None:
        charge(p, "volume", f"{volume} mL", {"number": volume, "unit":"mL"}, evidence, method)
    if product_role is not None:
        p["product_role"] = fx(product_role, product_role, evidence, method=method)
    return p


def manual_classification(r, evidence, label_codes, why):
    changes=[]
    support={}
    num={"C-C":2,"C-N":5}
    pairs={"C-C":("C","C"),"C-N":("C","N")}
    for i,code in enumerate(label_codes,1):
        change_id=f"b{i}"
        a,b=pairs[code]
        ch={"change_id":change_id,"atoms":[{"atom_map":None,"element":a},{"atom_map":None,"element":b}],
            "before_order":0,"after_order":1,"change_type":"formed"}
        changes.append(derived(ch,evidence,why))
        support.setdefault(code,[]).append(change_id)
    labels=[{"category_number":num[c],"code":c,"supporting_change_ids":support[c]} for c in support]
    r["classification"]["mapped_reaction_smiles"] = missing(evidence)
    r["classification"]["bond_changes"] = changes
    r["classification"]["labels"] = derived(labels,evidence,why)
    r["classification"]["coverage"] = inferred("partial",evidence,
        "Only source-legible net formed-bond classes are represented; full atom mapping, broken bonds and proton transfers remain unresolved.")
    r["classification"]["limitations"] = [
        "No atom-mapped reaction SMILES was asserted because several source structures are generic or stereochemically incomplete.",
        "C-C/C-N labels describe net product connectivity, not an elementary mechanism.",
        "Implicit-H changes and leaving-group coproducts are unresolved."]
    r["classification"]["model_votes"] = []


def no_reaction_classification(r, evidence, why):
    r["classification"]["mapped_reaction_smiles"] = missing(evidence)
    r["classification"]["bond_changes"] = []
    r["classification"]["labels"] = inferred([], evidence, why)
    r["classification"]["coverage"] = inferred("not_applicable", evidence, why)
    r["classification"]["limitations"] = [why]
    r["classification"]["model_votes"] = []


def base_reaction(rid, page, category, location, participants, *, time=3,
                  solvent_mid=None, solvent_name=None, solvent_volume=None,
                  evidence=None, method="original_extraction"):
    evidence = evidence or ev(page)
    r=st.new_reaction(f"{DOC}:{rid}",DOC,evidence)
    r["record_type"]="specific_experiment"
    r["experiment_categories"]=fx(category,[category],evidence,method=method)
    r["location"]=fx(location,{"page":page,"location":location},evidence,method=method)
    r["reaction_name"]=fx("aryne-mediated [2,3]-sigmatropic rearrangement",
                           "aryne-mediated [2,3]-sigmatropic rearrangement",evidence,method=method)
    r["transformation_type"]=inferred("sigmatropic rearrangement",evidence,
        "Normalized from the paper's reaction name; not an independent mechanistic proof.")
    r["participants"]=participants
    if solvent_mid:
        r["solvent_conditions"]["composition"]=fx(solvent_name, [solvent_name], evidence, method=method)
    if solvent_volume is not None:
        r["solvent_conditions"]["total_volume"]=fx(f"{solvent_volume} mL",{"number":solvent_volume,"unit":"mL"},evidence,method=method)
    if any(p["molecule_id"]=="m_1a" and p["charges"]["amount"]["normalized_value"] for p in participants):
        r["reaction_scale"]=fx("1a (0.20 mmol)",{"reference_molecule_id":"m_1a","amount":{"number":0.20,"unit":"mmol"}},evidence,method=method)
    step=st.new_step("reaction",evidence)
    step["participant_ids"]=[p["participant_id"] for p in participants if p["stage"]=="reaction"]
    step["operation"]=fx("combined under the stated reaction conditions","reaction under stated conditions",evidence,method=method)
    step["temperature"]=fx("rt","rt",evidence,method=method)
    step["time"]=fx(f"{time} h",{"number":time,"unit":"h"},evidence,method=method)
    r["steps"]=[step]
    r["procedure"]["original_text"]=fx(location,location,evidence,method=method)
    r["review"]={"status":"needs_review","method":"model","reviewer":None,"checked_fields":[],
                 "unresolved":["No explicitly associated SI was supplied; detailed addition order, workup and purification remain unavailable.",
                               "Full atom mapping and implicit-hydrogen accounting remain unresolved."]}
    r["completeness"]="partial"
    return r


def outcome(oid, kind, targets, raw, norm, evidence, *, yield_type=None,
            basis=None, scope=None, method="original_extraction"):
    o=st.new_outcome(oid,kind,targets,evidence)
    o["value"]=fx(raw,norm,evidence,method=method)
    if yield_type is not None:o["yield_type"]=fx(yield_type,yield_type,evidence,method=method)
    if basis is not None:o["measurement_basis"]=fx(basis,basis,evidence,method=method)
    if scope is not None:o["scope"]=fx(scope,scope,evidence,method=method)
    return o


reactions=[]

# Scheme 2 generic design record (not a concrete experiment).
e=ev(2)
g=st.new_reaction(f"{DOC}:scheme2_generic",DOC,e)
g["record_type"]="generic_scheme"
g["experiment_categories"]=fx("generic design",["other"],e)
g["location"]=fx("Scheme 2, page 2",{"page":2,"scheme":2},e)
g["reaction_name"]=fx("Design of a new [2,3]-sigmatropic rearrangement","[2,3]-sigmatropic rearrangement",e)
g["transformation_type"]=inferred("generic aryne activation / [2,3]-sigmatropic rearrangement / optional base isomerization",e,
    "Normalized from the arrows and caption; the scheme is a proposal/design, not a concrete experiment.")
g["participants"]=[participant("g_sub","m_g1","reactant",e),participant("g_aryne","m_g2","reactant",e),
                   participant("g_prod","m_g3","product",e,product_role="generic rearrangement product"),
                   participant("g_diene","m_g4","product",e,product_role="generic base-isomerized product")]
manual_classification(g,e,["C-N","C-C"],"Net bond classes inferred from the generic Scheme 2 drawings.")
g["notes"]=[fx("Scheme 2 is a design/proposed sequence and reports no isolated experiment.",
                 "generic scheme, not a specific experiment",e)]
g["review"]={"status":"needs_review","method":"model","reviewer":None,"checked_fields":[],
             "unresolved":["Generic R/R′/EWG structures and aryne regioisomers are not uniquely instantiated."]}
g["completeness"]="partial"
reactions.append(g)


# Table 1: every optimization row is an independent experiment.
t1=[
 (1,"m_CsF","CsF",False,"m_MeCN","MeCN",34,"83 : 17",3),
 (2,"m_KF","KF",False,"m_MeCN","MeCN","Trace",None,3),
 (3,"m_ZnF2","ZnF2",False,"m_MeCN","MeCN",0,None,3),
 (4,"m_NaF","NaF",False,"m_MeCN","MeCN",0,None,3),
 (5,"m_NH4F","NH4F",False,"m_MeCN","MeCN","Trace",None,3),
 (6,"m_TBAF","TBAF",False,"m_MeCN","MeCN",82,"2 : 98",3),
 (7,"m_CsF","CsF",True,"m_MeCN","MeCN",65,"30 : 70",3),
 (8,"m_KF","KF",True,"m_MeCN","MeCN",70,"84 : 16",3),
 (9,"m_KF","KF",True,"m_THF","THF",79,"96 : 4",3),
 (10,"m_KF","KF",True,"m_DMF","DMF",40,"<1 : 99",3),
 (11,"m_KF","KF",True,"m_dioxane","1,4-dioxane",57,"97 : 3",3),
 (12,"m_KF","KF",True,"m_toluene","toluene",37,"92 : 8",3),
 (13,"m_KF","KF",True,"m_DCE","DCE",68,"92 : 8",3),
 (14,"m_KF","KF",True,"m_DEE","EtOCH2CH2OEt",70,"97 : 3",3),
 (15,"m_KF","KF",True,"m_THF","THF",76,"84 : 16",20),
]
for entry,fmid,fname,crown,smid,sname,yld,ratio,hours in t1:
    e=ev(2)
    ps=[participant("p_1a","m_1a","reactant",e,amount=.20),
        participant("p_2a","m_2a","reactant",e,amount=.24),
        participant("p_f",fmid,"reagent",e,amount=.40),
        participant("p_solvent",smid,"solvent",e,volume=1.0),
        participant("p_3a","m_3a","product",e,product_role="combined-yield product"),
        participant("p_3ap","m_3ap","product",e,product_role="combined-yield isomerization product")]
    if crown: ps.insert(3,participant("p_crown","m_18c6","additive",e,amount=.40))
    r=base_reaction(f"table1_entry_{entry:02d}",2,"optimization",f"Table 1 entry {entry}",ps,
                    time=hours,solvent_mid=smid,solvent_name=sname,solvent_volume=1.0,evidence=e)
    if entry==6:
        r["notes"].append(fx("A solution of TBAF (1.0 M) in THF was used.",
                             {"TBAF_solution_concentration":{"number":1.0,"unit":"M"},"carrier":"THF","carrier_volume":"unreported"},e))
    val=yld if isinstance(yld,str) else {"number":yld,"unit":"%"}
    r["outcomes"].append(outcome("combined_yield","yield",["m_3a","m_3ap"],
        f"{yld}%" if isinstance(yld,int) else yld,val,e,yield_type="combined isolated yield",basis="isolated"))
    if ratio is not None:
        r["outcomes"].append(outcome("ratio_3a_3ap","selectivity",["m_3a","m_3ap"],ratio,
            {"ratio":ratio,"order":["3a","3a′"]},e,basis="1H NMR analysis of crude products",
            scope="3a:3a′ ratio"))
    else:
        r["outcomes"].append(outcome("ratio_not_reported","other",["m_3a","m_3ap"],"—","—",e,
            scope="Table dash for 3a:3a′ ratio"))
    zerofail=(isinstance(yld,int) and yld==0)
    if zerofail:
        no_reaction_classification(r,e,"Table reports 0% combined yield and no observed product ratio; no product bond class is asserted.")
    else:
        manual_classification(r,e,["C-N","C-C"],"Net C-N arylation and C-C rearrangement connectivity inferred from the Table 1 reactant/product drawings.")
    if entry==6:
        r["outcomes"].append(outcome("z_e_3ap","selectivity",["m_3ap"],"53 : 47 Z/E",
            {"ratio":"53:47","order":["Z","E"]},e,basis="1H NMR analysis of crude products",scope="3a′ Z/E ratio"))
    elif not zerofail:
        r["outcomes"].append(outcome("z_e_3ap","selectivity",["m_3ap"],"nearly 1 : 1 Z/E",
            {"ratio":"nearly 1:1","order":["Z","E"]},e,basis="1H NMR analysis of crude products",scope="3a′ Z/E ratio"))
    reactions.append(r)


def standard_scope_participants(sub_mid,prod_mids,evidence, precursor_mid="m_2a", precursor_amount=.24,
                                method="original_extraction"):
    ps=[participant("p_sub",sub_mid,"reactant",evidence,amount=.20,method=method),
        participant("p_aryne",precursor_mid,"reactant",evidence,amount=precursor_amount,method=method),
        participant("p_KF","m_KF","reagent",evidence,amount=.40,method=method),
        participant("p_crown","m_18c6","additive",evidence,amount=.40,method=method),
        participant("p_THF","m_THF","solvent",evidence,volume=1.0,method=method)]
    for i,(mid,role) in enumerate(prod_mids,1):
        ps.append(participant(f"p_prod{i}",mid,"product",evidence,product_role=role,method=method))
    return ps


# Table 2 substrate scope (standard-scale experiments).
for entry,(sub,prod,yld,zr,ee) in t2_rows.items():
    e=ev(3) + (ev(5) if entry==5 else [])
    submid=f"m_{sub}"
    prodmid=f"m_{prod.replace('′','p')}"
    prod_list=[(prodmid,"reported target product")]
    if entry==5:prod_list.append(("m_side_1e","formal-substitution byproduct"))
    ps=standard_scope_participants(submid,prod_list,e)
    r=base_reaction(f"table2_entry_{entry:02d}",3,"substrate_scope",f"Table 2 entry {entry}",ps,
                    solvent_mid="m_THF",solvent_name="THF",solvent_volume=1.0,evidence=e)
    r["reaction_scale"]=fx(f"{sub} (0.20 mmol)",{"reference_molecule_id":submid,"amount":{"number":.20,"unit":"mmol"}},e)
    r["outcomes"].append(outcome("isolated_yield","yield",[prodmid],f"{yld}%",{"number":yld,"unit":"%"},e,
        yield_type="isolated yield",basis="isolated"))
    if entry==5:
        r["outcomes"].append(outcome("side_product_yield","yield",["m_side_1e"],"30%",{"number":30,"unit":"%"},ev(5),
            yield_type="reported yield",basis=None,scope="formal substitution of the 2,3-allenyl group with phenyl"))
        r["review"]["unresolved"].append("Footnote 13 does not state the 30% side-product yield measurement basis.")
    if zr:
        r["outcomes"].append(outcome("z_e_ratio","selectivity",[prodmid],zr,
            {"ratio":zr.replace(" Z/E",""),"order":["Z","E"]},ev(3),basis="1H NMR analysis of crude products",scope="Z/E ratio"))
    if ee:
        r["outcomes"].append(outcome("product_ee","selectivity",[prodmid],ee,
            {"number":float(ee.replace("% ee","")),"unit":"%","metric":"ee"},ev(3),
            basis="HPLC analysis on a chiral stationary phase",scope="product enantiomeric excess; absolute configuration unassigned"))
    if yld==0 and entry==11:
        no_reaction_classification(r,e,"Table reports 0% desired product and the prose reports a complex mixture; no observed product class is asserted.")
    elif entry==5:
        manual_classification(r,e,["C-N"],"Observed formal substitution product contains the new N-phenyl bond; desired [2,3]-product yield is 0%.")
    else:
        manual_classification(r,e,["C-N","C-C"],"Net C-N arylation and C-C rearrangement connectivity inferred from the Table 2 drawings.")
    reactions.append(r)

# Separate 6 mmol scale-up for 1a -> 3a, Table 2 entry 1 footnote b.
e=ev(3)
ps=[participant("p_sub","m_1a","reactant",e,amount=6.0),
    participant("p_aryne","m_2a","reactant",e),participant("p_KF","m_KF","reagent",e),
    participant("p_crown","m_18c6","additive",e),participant("p_THF","m_THF","solvent",e),
    participant("p_prod","m_3a","product",e,product_role="scale-up product")]
r=base_reaction("table2_entry_01_scaleup",3,"scale_up","Table 2 entry 1, parenthetical 78% and footnote b",ps,
                solvent_mid="m_THF",solvent_name="THF",evidence=e,method="inherited_general_condition")
r["reaction_scale"]=fx("6 mmol scale",{"reference_molecule_id":"m_1a","amount":{"number":6.0,"unit":"mmol"}},e)
r["outcomes"]=[outcome("isolated_yield","yield",["m_3a"],"78%",{"number":78,"unit":"%"},e,
                       yield_type="isolated yield",basis="isolated")]
r["review"]["unresolved"].append("Footnote b gives only the 6 mmol scale; scaled reagent quantities and solvent volume are not reported in the supplied main text and were not calculated.")
manual_classification(r,e,["C-N","C-C"],"Net connectivity inferred from the 1a/2a/3a drawings; scale-up does not change the bond-class assignment.")
reactions.append(r)


# Scheme 3 competing migrations.
scheme3=[("a","m_1r",[("m_3r","2,3-allenyl migration product"),("m_4a","allyl migration product")],76,"43 : 57",["3r","4a"]),
         ("b","m_1s",[("m_3s","2,3-allenyl migration product"),("m_5a","propargyl migration product not observed")],53,"100 : 0",["3s","5a"]),
         ("c","m_1t",[("m_3t","2,3-allenyl migration product"),("m_6a","benzyl [1,2]-migration product")],75,"80 : 20",["3t","6a"])]
for letter,sub,prods,yld,ratio,order in scheme3:
    e=ev(3)
    ps=standard_scope_participants(sub,prods,e,method="inherited_general_condition")
    r=base_reaction(f"scheme3_{letter}",3,"substrate_scope",f"Scheme 3{letter}",ps,
                    solvent_mid="m_THF",solvent_name="THF",solvent_volume=1.0,evidence=e,method="inherited_general_condition")
    r["reaction_scale"]=fx("standard conditions inherited from the page-3 condition set",
        {"reference_molecule_id":sub,"amount":{"number":.20,"unit":"mmol"}},e,method="inherited_general_condition")
    targets=[mid for mid,_ in prods]
    r["outcomes"]=[outcome("combined_yield","yield",targets,f"{yld}%",{"number":yld,"unit":"%"},e,
                            yield_type="combined reported yield",basis=None),
                   outcome("product_ratio","selectivity",targets,ratio,{"ratio":ratio,"order":order},e,
                           basis=None,scope="product ratio")]
    r["review"]["unresolved"].append("Scheme 3 does not state the yield or product-ratio measurement basis in the supplied main text.")
    manual_classification(r,e,["C-N","C-C","C-C"],"Competing net rearrangement branches share N-arylation and form branch-specific C-C bonds; atom positions remain unmapped.")
    reactions.append(r)


# Table 3 aryne scope.
t3=[(1,"m_2b",[("m_3u","main product")],92,None),
    (2,"m_2c",[("m_3v","regioisomer"),("m_3u4","alternate regioisomer")],60,"45 : 55"),
    (3,"m_2d",[("m_3w","main product")],72,None),
    (4,"m_2e",[("m_3x","main product")],76,None),
    (5,"m_2f",[("m_3y","main product")],70,None),
    (6,"m_2g",[("m_3z","desired product not isolated")],0,None)]
for entry,prec,prods,yld,ratio in t3:
    e=ev(3)+(ev(4) if entry==6 else [])
    ps=standard_scope_participants("m_1a",prods,e,precursor_mid=prec)
    r=base_reaction(f"table3_entry_{entry:02d}",3,"substrate_scope",f"Table 3 entry {entry}",ps,
                    solvent_mid="m_THF",solvent_name="THF",solvent_volume=1.0,evidence=e)
    targets=[mid for mid,_ in prods]
    r["outcomes"]=[outcome("isolated_yield","yield",targets,f"{yld}%",{"number":yld,"unit":"%"},ev(3),
                            yield_type="isolated yield",basis="isolated")]
    if ratio:
        r["outcomes"].append(outcome("regioisomer_ratio","selectivity",targets,ratio,
            {"ratio":ratio,"order":["3v","3u (R = 4-Me)"]},ev(3),basis="1H NMR analysis of crude product",scope="regioisomer ratio"))
    if entry==6:
        r["outcomes"].append(outcome("complex_mixture","decomposition",[],"complex mixture","complex mixture",ev(4),scope="desired product not isolated"))
        no_reaction_classification(r,e,"The reaction gave a complex mixture and 0% desired product; no product bond class is asserted.")
    else:
        manual_classification(r,e,["C-N","C-C"],"Net C-N arylation and C-C rearrangement connectivity inferred from the Table 3 drawings.")
    reactions.append(r)


bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
bundle["source_context"] = [
  {"context_id":f"{DOC}:metadata","document_id":DOC,"evidence_ids":ev(1),
   "kind":"bibliographic_metadata","title":"Aryne-mediated [2,3]-sigmatropic rearrangement of tertiary 2,3-allenylamines bearing an electron-withdrawing group at the alpha-position",
   "doi":"10.1039/d1ob00887k","authors":["Lu Han","Xue-Ting Zhang","Dong Xie","Shi-Kai Tian"]},
  {"context_id":f"{DOC}:scheme1","document_id":DOC,"evidence_ids":ev(1),"kind":"prior_work_and_conceptual_comparison",
   "note":"Scheme 1 combines cited prior-work reaction classes with a conceptual 'this work' summary; it is not counted as a concrete experiment."},
  {"context_id":f"{DOC}:mechanistic_discussion","document_id":DOC,"evidence_ids":ev(1)+ev(2),"kind":"author_mechanistic_discussion",
   "note":"Quaternary ammonium ylide generation, [2,3]-rearrangement and base-promoted isomerization are author discussion/design, not separately verified mechanistic experiments."},
  {"context_id":f"{DOC}:back_matter","document_id":DOC,"evidence_ids":ev(4)+ev(5),"kind":"non_experimental_back_matter",
   "note":"Conclusions, conflicts, acknowledgements and references are not current-paper reaction records; footnote 13 is linked to Table 2 entry 5."},
]

OUT.write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
inventory={"document_id":DOC,"pages_reviewed":[1,2,3,4,5],"specific_experiments":42,
           "generic_schemes":1,"reaction_ids":[r["reaction_id"] for r in reactions],
           "notes":["All five rendered page images were visually inspected.",
                    "No SI document was explicitly supplied or associated.",
                    "Raw table rows and scheme arrows were transcribed into the generated submission; all records remain machine needs_review."]}
(RUN/"raw_model/page_inventory.json").write_text(json.dumps(inventory,ensure_ascii=False,indent=2)+"\n")
diagnostics={}
for m in bundle["molecules"]:
    smi=m["smiles"]["normalized_value"]
    if smi:
        diagnostics[m["molecule_id"]]={"input_smiles":smi,"diagnostic":st.stereo_inspect(smi)}
(RUN/"raw_model/stereo_diagnostics.json").write_text(json.dumps(diagnostics,ensure_ascii=False,indent=2)+"\n")
print(json.dumps({"molecules":len(bundle["molecules"]),"reactions":len(reactions),"specific":42,"generic":1},ensure_ascii=False))
