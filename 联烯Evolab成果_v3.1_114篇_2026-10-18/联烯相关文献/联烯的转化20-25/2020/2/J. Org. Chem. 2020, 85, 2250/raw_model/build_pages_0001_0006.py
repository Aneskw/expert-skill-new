#!/usr/bin/env python3
"""Build evidence-linked v3.1 records for source PDF pages 1--6.

This is a source-specific extraction artifact. It does not modify the tools,
schema, or scientific validation rules.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_2edd574f8525f831"
OUT = ROOT / "runs/current"
E = {p: f"{DOC}:page{p}:image" for p in range(1, 7)}
T = {p: f"{DOC}:page{p}:text" for p in range(1, 7)}


def F(raw=None, norm=None, pages=(), method="original_extraction", derivation=None,
      status=None, alternatives=()):
    ev = []
    for p in pages:
        ev.extend([E[p], T[p]])
    return st.fact(raw, norm, ev, method=method, derivation=derivation,
                   status=status, alternatives=alternatives)


def missing(pages=(), status="unreported"):
    return F(pages=pages, status=status)


def q(number, unit, raw=None, qualifier=None):
    x = {"number": number, "unit": unit}
    if qualifier:
        x["qualifier"] = qualifier
    if raw:
        x["raw_text"] = raw
    return x


def srcraw(value):
    """Render a normalized quantity back to a source-like raw expression."""
    if isinstance(value, dict):
        if value.get("raw_text"):
            return value["raw_text"]
        unit = value.get("unit", "")
        if "number" in value:
            mark = {"lt":"<","le":"≤","gt":">","ge":"≥","approx":"~"}.get(value.get("qualifier"), "")
            return f"{mark}{value['number']:g} {unit}".strip()
        if "minimum" in value and "maximum" in value:
            return f"{value['minimum']:g}–{value['maximum']:g} {unit}".strip()
        return json.dumps(value, ensure_ascii=False)
    return str(value)


def molecule(mid, label, name, pages, smiles=None, chirality=None,
             geometry=None, representation=None, assignment=None,
             mixture=None, original=None):
    ev = [z for p in pages for z in (E[p], T[p])]
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, pages)
    m["names"]["original"] = F(original or name, original or name, pages)
    m["names"]["english"] = F(name, name, pages)
    m["cas_numbers"] = missing(pages)
    m["identity_note"] = F(
        "Structure identity remains a machine candidate pending expert comparison with the source drawing.",
        "machine candidate; source drawing not chemically certified", pages,
        method="inference", derivation="Record is generated from the visually inspected drawing/name; validation is not structure certification.")
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = F(None, smiles, pages, method="structure_derivation",
                         derivation="Candidate connectivity/isomeric SMILES derived from the source drawing and compound name; not source-verified.")
        diag = st.stereo_inspect(smiles)
        m["stereochemistry"]["isomeric_smiles"] = F(
            None, diag["canonical_isomeric_smiles"], pages, method="structure_derivation",
            derivation="RDKit canonical isomeric SMILES calculated from the candidate source-derived graph.")
        m["stereochemistry"]["identity_key"] = F(
            None, diag["stereo_graph_key"], pages, method="structure_derivation",
            derivation="RDKit stereo graph key calculated from the candidate graph; not a sample identity certificate.")
        if chirality is None:
            chirality = "unspecified" if any(x.get("specified") == "Unspecified" and str(x.get("type", "")).startswith("Atom_") for x in diag["potential_stereo"]) else "achiral"
        if geometry is None:
            geometry = "unspecified" if any(x.get("specified") == "Unspecified" and "Bond" in str(x.get("type", "")) for x in diag["potential_stereo"]) else "not_applicable"
        representation = representation or ("encoded" if geometry == "specified" else "not_applicable")
    else:
        m["representation_type"] = "unresolved"
        chirality = chirality or "unresolved"
        geometry = geometry or "unresolved"
        representation = representation or "unresolved"
    if geometry == "specified":
        m["stereochemistry"]["double_bonds"] = F(
            "(1E) in source name/drawing", [{"configuration": "E", "source_locator": "product alkene"}], pages)
    if assignment is not None:
        m["stereochemistry"]["assignment"] = F(assignment, assignment, pages)
    if mixture is not None:
        m["stereochemistry"]["mixture"] = F(mixture, mixture, pages)
    for key, val in (("chirality_status", chirality), ("geometry_status", geometry),
                     ("representation_status", representation)):
        m["stereochemistry"][key] = F(None, val, pages, method="inference",
            derivation="Status assigned from the candidate molecular graph together with explicit source stereochemical labels; absence of a label is not treated as racemic.")
    return m


bundle = st.read(OUT / "raw_model/pages_0001_0006_template.json")
mols = {}


def M(mid, label=None, name=None, pages=(2,), smiles=None, **kw):
    if mid not in mols:
        mols[mid] = molecule(mid, label or mid, name or label or mid, pages, smiles, **kw)
    return mid


# Frequently used reagents/materials. Small-molecule SMILES are candidates only.
M("m_eosinY", "EY", "eosin Y", (1, 2, 3, 4, 5, 6))
M("m_rose", "Rose bengal", "Rose bengal", (2,))
M("m_rubpy", "Ru(bpy)3Cl2·6H2O", "tris(2,2'-bipyridine)ruthenium(II) chloride hexahydrate", (2,))
M("m_ptsna", "p-TolSO2Na", "sodium 4-methylbenzenesulfinate", (2, 3, 4, 5, 6), "[Na+].[O-]S(=O)c1ccc(C)cc1")
M("m_phsna", "PhSO2Na", "sodium benzenesulfinate", (3, 4, 6), "[Na+].[O-]S(=O)c1ccccc1")
M("m_brsna", "4-BrC6H4SO2Na", "sodium 4-bromobenzenesulfinate", (3,), "[Na+].[O-]S(=O)c1ccc(Br)cc1")
M("m_nosna", "4-NO2C6H4SO2Na", "sodium 4-nitrobenzenesulfinate", (3,), "[Na+].[O-]S(=O)c1ccc([N+](=O)[O-])cc1")
M("m_thsna", "2-thienylSO2Na", "sodium thiophene-2-sulfinate", (3, 6), "[Na+].[O-]S(=O)c1cccs1")
M("m_mesna", "MeSO2Na", "sodium methanesulfinate", (3, 6), "CS(=O)[O-].[Na+]")
M("m_cfsna", "CF3SO2Na", "sodium trifluoromethanesulfinate (Langlois reagent)", (3,), "[Na+].[O-]S(=O)C(F)(F)F")
M("m_k2co3", "K2CO3", "potassium carbonate", (2, 3, 4, 5, 6), "[K+].[K+].[O-]C([O-])=O")
M("m_tempo", "TEMPO", "2,2,6,6-tetramethylpiperidin-1-oxyl", (4,), "CC1(C)CCCC(C)(C)N1[O]")
for mid, lab, name, smi, pp in [
    ("m_dmf", "DMF", "N,N-dimethylformamide", "CN(C)C=O", (1,2,3,4,5,6)),
    ("m_tfe", "TFE", "2,2,2-trifluoroethanol", "OCC(F)(F)F", (1,2,3,4,5,6)),
    ("m_water", "H2O", "water", "O", (2,)),
    ("m_d2o", "D2O", "deuterium oxide", "[2H]O[2H]", (3,4)),
    ("m_etoh", "EtOH", "ethanol", "CCO", (2,)),
    ("m_chcl3", "CHCl3", "chloroform", "ClC(Cl)Cl", (2,)),
    ("m_hfip", "HFIP", "hexafluoroisopropanol", "OC(C(F)(F)F)C(F)(F)F", (2,)),
    ("m_hoac", "AcOH", "acetic acid", "CC(=O)O", (4,)),
    ("m_o2", "O2", "oxygen", "O=O", (4,)),
    ("m_toluene", "toluene", "toluene", "Cc1ccccc1", (4,)),
    ("m_koh", "KOH", "potassium hydroxide", "[K+].[OH-]", (4,5)),
    ("m_tbab", "TBAB", "tetrabutylammonium bromide", "CCCC[N+](CCCC)(CCCC)CCCC.[Br-]", (4,)),
    ("m_dmso", "DMSO", "dimethyl sulfoxide", "CS(C)=O", (4,5)),
    ("m_propargyl_br", "propargyl bromide", "propargyl bromide", "C#CCBr", (4,5)),
    ("m_paraform", "paraformaldehyde", "paraformaldehyde", None, (5,)),
    ("m_cubr", "CuBr", "copper(I) bromide", "[Cu]Br", (5,)),
    ("m_dipea", "NHiPr2", "diisopropylamine", "CC(C)NC(C)C", (5,)),
    ("m_dioxane", "1,4-dioxane", "1,4-dioxane", "O1CCOCC1", (5,)),
    ("m_et2o", "Et2O", "diethyl ether", "CCOCC", (5,)),
    ("m_etoac", "EtOAc", "ethyl acetate", "CCOC(C)=O", (4,5)),
    ("m_hexane", "hexane", "hexane", "CCCCCC", (4,5)),
    ("m_na2so4", "Na2SO4", "sodium sulfate", "[Na+].[Na+].[O-]S(=O)(=O)[O-]", (5,)),
]: M(mid, lab, name, pp, smi)

# Allene substrates and E-vinyl-sulfone products.
ind = {
    "1a": "n1c(-c2ccccc2)cc2ccccc21",
    "1b": "n1c(-c2ccc(C)cc2)cc2ccccc21",
    "1c": "n1ccc2ccccc21",
    "1d": "n1cccc1",
    "1e": "n1ccc2cc(OC)ccc21",
    "1f": "n1ccc2cc(Br)ccc21",
    "1g": "n1ccc2ccc(Br)cc21",
    "1h": "n1ccc2cccc(C#N)c21",
    "1i": "n1cc(C)c2ccccc21",
}
names1 = {
    "1a":"2-phenyl-1-propadienyl-1H-indole", "1b":"2-(4-methylphenyl)-1-propadienyl-1H-indole",
    "1c":"1-propadienyl-1H-indole", "1d":"1-propadienyl-1H-pyrrole",
    "1e":"5-methoxy-1-propadienyl-1H-indole", "1f":"5-bromo-1-propadienyl-1H-indole",
    "1g":"6-bromo-1-propadienyl-1H-indole", "1h":"4-cyano-1-propadienyl-1H-indole",
    "1i":"3-methyl-1-propadienyl-1H-indole",
}
for lab, frag in ind.items():
    M(f"m_{lab}", lab, names1[lab], (2,3,4,5,6), "C=C=C"+frag)

sulf = {"Ts":"S(=O)(=O)c1ccc(C)cc1", "Ph":"S(=O)(=O)c1ccccc1",
        "Th":"S(=O)(=O)c1cccs1", "Me":"S(=O)(=O)C"}
product_defs = {
    "2a":("1a","Ts"), "2b":("1i","Ts"), "2c":("1f","Ts"), "2d":("1g","Ph"),
    "2e":("1e","Ts"), "2f":("1h","Ts"), "2g":("1g","Th"), "2h":("1c","Ts"),
    "2i":("1c","Ph"), "2j":("1c","Me"), "2k":("1e","Th"), "2m":("1d","Ph"),
    "2l":("1b","Ts"),
}
for lab,(sub,so) in product_defs.items():
    smi = f"C/C(=C\\{ind[sub]}){sulf[so]}"
    M(f"m_{lab}", lab, f"(E)-vinyl sulfone {lab}", (2,3,5,6), smi,
      geometry="specified", representation="encoded", assignment="Source drawing/name specifies the E alkene.")
mols["m_2a"]["stereochemistry"]["assignment"] = F(
    "The structure of 2a was unambiguously determined by single-crystal X-ray diffraction (CCDC 1922457).",
    {"method":"single-crystal X-ray diffraction","ccdc":"1922457","source_claim":"unambiguously determined"},
    (2,))
M("m_2lprime", "2l′", "allylic regioisomer 2l′", (3,),
  f"C=C(C{ind['1b']}){sulf['Ts']}", geometry="unspecified", representation="not_encoded")

aryl_allenes = {
    "3a":("propadienylbenzene","C=C=Cc1ccccc1"),
    "3b":("2-propadienylnaphthalene","C=C=Cc1ccc2ccccc2c1"),
    "3c":("1-methyl-4-propadienylbenzene","C=C=Cc1ccc(C)cc1"),
    "3d":("1-methoxy-4-propadienylbenzene","C=C=Cc1ccc(OC)cc1"),
    "3e":("buta-2,3-dien-1-ylbenzene","C=C=CCc1ccccc1"),
    "3f":("(propadienyloxy)benzene","C=C=COc1ccccc1"),
    "3g":("1-(buta-2,3-dien-1-yloxy)-4-(propan-2-yl)benzene","C=C=COc1ccc(C(C)C)cc1"),
    "3h":("1-propadienyl-4-(trifluoromethyl)benzene","C=C=Cc1ccc(C(F)(F)F)cc1"),
}
for lab,(name,smi) in aryl_allenes.items(): M(f"m_{lab}",lab,name,(3,4,5),smi)

def vinyl(aryl, so):
    return f"C/C(=C\\{aryl}){so}"

prod4 = {
 "4a":vinyl("c1ccccc1",sulf["Ts"]), "4b":vinyl("c1ccccc1",sulf["Ph"]),
 "4c":vinyl("c1ccccc1","S(=O)(=O)c1ccc(Br)cc1"),
 "4d":vinyl("c1ccccc1","S(=O)(=O)c1ccc([N+](=O)[O-])cc1"),
 "4e":vinyl("c1ccccc1",sulf["Th"]), "4f":vinyl("c1ccccc1",sulf["Me"]),
 "4g":vinyl("c1ccc(C)cc1",sulf["Ts"]), "4h":vinyl("c1ccc2ccccc2c1",sulf["Ts"]),
 "4i":vinyl("Oc1ccccc1",sulf["Ts"]), "4l":vinyl("c1ccc(OC)cc1",sulf["Ph"]),
 "4m":vinyl("c1ccc(C(F)(F)F)cc1",sulf["Ph"]),
 "4n":vinyl("c1ccc(C(F)(F)F)cc1","S(=O)(=O)c1ccc(Br)cc1"),
}
for lab,smi in prod4.items():
    M(f"m_{lab}",lab,f"(E)-vinyl sulfone {lab}",(3,4),smi,geometry="specified",representation="encoded",assignment="Source drawing/name specifies the E alkene.")
M("m_4bprime","4b′","regioisomeric vinyl sulfone 4b′",(4,),"C=C(Cc1ccccc1)S(=O)(=O)c1ccccc1",geometry="unspecified",representation="not_encoded")
M("m_indole","indole","1H-indole",(4,),"c1ccc2[nH]ccc2c1")
M("m_ctrl1d","deuterated vinyl sulfone","deuterated vinyl sulfone from 3a",(4,),chirality="unresolved",geometry="specified",representation="not_encoded",assignment="Source drawing places D in the methyl group; exact isotope-aware E SMILES was not asserted.")
M("m_ctrl2d","deuterated vinyl sulfone","deuterated indolyl vinyl sulfone",(4,),chirality="unresolved",geometry="specified",representation="not_encoded",assignment="Source drawing places D in the methyl group; exact isotope-aware E SMILES was not asserted.")

# Precursors for starting-allene syntheses.
M("m_propargyl_1a","2-phenyl-N-propargyl indole","2-phenyl-N-propargyl indole",(4,),"C#CCn1c(-c2ccccc2)cc2ccccc21")
M("m_propargyl_1b","2-(4-methylphenyl)-N-propargyl indole","2-(4-methylphenyl)-N-propargyl indole",(4,),"C#CCn1c(-c2ccc(C)cc2)cc2ccccc21")
substrate_names = {
 "1c":"indole", "1d":"pyrrole", "1e":"5-methoxyindole", "1f":"5-bromoindole",
 "1g":"6-bromoindole", "1h":"4-cyanoindole", "1i":"3-methylindole"}
substrate_smi = {
 "1c":"c1ccc2[nH]ccc2c1", "1d":"c1cc[nH]c1", "1e":"COc1ccc2[nH]ccc2c1",
 "1f":"Brc1ccc2[nH]ccc2c1", "1g":"Brc1cc2cc[nH]c2cc1", "1h":"N#Cc1cccc2[nH]ccc12",
 "1i":"Cc1c[nH]c2ccccc12"}
for lab in substrate_names:
    M(f"m_start_{lab}",substrate_names[lab],substrate_names[lab],(4,5),substrate_smi[lab])
M("m_alkyne_3g","alkyne precursor to 3g","1-ethynylmethyl ether precursor to 3g",(5,))

reactions = []


def base_reaction(rid, pages, rtype="specific_experiment", categories=None,
                  name="", trans="", location=""):
    ev = [z for p in pages for z in (E[p], T[p])]
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = rtype
    r["experiment_categories"] = F(str(categories or []), categories or [], pages)
    r["location"] = F(location, {"pages": list(pages), "region": location}, pages)
    r["reaction_name"] = F(name, name, pages) if name else missing(pages)
    r["transformation_type"] = F(None, trans, pages, method="inference", derivation="Transformation type assigned from the observed source reaction.") if trans else missing(pages)
    r["review"]["unresolved"] = ["Machine extraction requires expert verification of structure drawings, chemical identities, and completeness."]
    r["review"]["status"] = "needs_review"
    r["review"]["method"] = "model"
    return r


def participant(r, pid, mid, role, pages, mass=None, amount=None, volume=None,
                equiv=None, molpct=None, ref=None, stage="reaction"):
    p = st.new_participant(pid, mid, role, [z for pg in pages for z in (E[pg],T[pg])], stage)
    vals = {"mass":mass,"amount":amount,"volume":volume,"equivalents":equiv,"mol_percent":molpct}
    for key,val in vals.items():
        if val is not None:
            p["charges"][key] = F(srcraw(val), val, pages)
    if equiv is not None or molpct is not None:
        p["equivalence_reference"] = F(ref, ref, pages)
    r["participants"].append(p)
    return p


def step(r, sid, pages, operation, participants=(), temp=None, time=None, atmosphere=None, stage="reaction"):
    s = st.new_step(sid, [z for pg in pages for z in (E[pg],T[pg])], stage)
    s["participant_ids"] = list(participants)
    s["operation"] = F(operation, operation, pages)
    if temp is not None: s["temperature"] = F(srcraw(temp), temp, pages)
    if time is not None: s["time"] = F(srcraw(time), time, pages)
    if atmosphere is not None: s["atmosphere"] = F(str(atmosphere), atmosphere, pages)
    r["steps"].append(s)
    return s


def outcome(r, oid, kind, targets, raw, norm, pages, basis=None, yield_type=None, scope=None,
            alternatives=()):
    o = st.new_outcome(oid, kind, targets, [z for pg in pages for z in (E[pg],T[pg])])
    o["value"] = F(raw, norm, pages, alternatives=alternatives,
                   status="conflict" if alternatives else None)
    if basis is not None: o["measurement_basis"] = F(basis, basis, pages)
    if yield_type is not None: o["yield_type"] = F(yield_type, yield_type, pages)
    if scope is not None: o["scope"] = F(str(scope), scope, pages)
    r["outcomes"].append(o)
    return o


def set_solvent(r, pages, composition, ratio=None, total=None, concentration=None, basis=None):
    r["solvent_conditions"]["composition"] = F(srcraw(composition), composition, pages)
    if ratio is not None: r["solvent_conditions"]["mixture_ratio"] = F(srcraw(ratio), ratio, pages)
    if total is not None: r["solvent_conditions"]["total_volume"] = F(srcraw(total), total, pages)
    if concentration is not None: r["solvent_conditions"]["concentration"] = F(srcraw(concentration), concentration, pages)
    if basis is not None: r["solvent_conditions"]["concentration_basis"] = F(srcraw(basis), basis, pages)


def set_light(r, pages, source="green LED", wavelength=None, power=None):
    r["special_conditions"]["light_source"] = F(source, source, pages)
    if wavelength is not None: r["special_conditions"]["wavelength"] = F(srcraw(wavelength), wavelength, pages)
    if power is not None: r["special_conditions"]["power"] = F(srcraw(power), power, pages)


def classify_hydrosulf(r, pages):
    changes = [
      {"change_id":"b1","atoms":[{"atom_map":None,"element":"C"},{"atom_map":None,"element":"S"}],"before_order":0,"after_order":1,"change_type":"formed"},
      {"change_id":"b2","atoms":[{"atom_map":None,"element":"C"},{"atom_map":None,"element":"H"}],"before_order":0,"after_order":1,"change_type":"hydrogen_gain","hydrogen_count":1},
      {"change_id":"b3","atoms":[{"atom_map":None,"element":"C"},{"atom_map":None,"element":"C"}],"before_order":2,"after_order":1,"change_type":"order_decreased"},
    ]
    r["classification"]["bond_changes"] = [F(None,c,pages,method="structure_derivation",derivation="Partial net bond inventory derived from the allene-to-vinyl-sulfone source drawing; atom correspondence is not asserted as complete.") for c in changes]
    labels = [
      {"category_number":4,"code":"C-S","supporting_change_ids":["b1"]},
      {"category_number":1,"code":"C-H","supporting_change_ids":["b2"]},
    ]
    r["classification"]["labels"] = F(None,labels,pages,method="structure_derivation",derivation="Taxonomy labels follow the partial net bond inventory; proton source and complete atom mapping remain unresolved.")
    r["classification"]["coverage"] = F(None,"partial",pages,method="inference",derivation="C-S formation, net C-H gain and one allene bond-order decrease are represented; coproducts and atom mapping are incomplete.")
    r["classification"]["limitations"] = ["No complete atom mapping.","Proton source is not assigned at the experiment-record level.","Unlisted spectators/coproducts are excluded."]


def classify_na(r, pages):
    r["classification"]["labels"] = F(None,[],pages,method="inference",derivation="No observed target transformation; bond formation classification is not applicable.")
    r["classification"]["coverage"] = F(None,"not_applicable",pages,method="inference",derivation="Negative/no-conversion experiment.")


# Page 1 graphical abstract and page 2 prior-art/this-work schemes.
r = base_reaction("graphical_abstract",(1,),"generic_scheme",[],"Photoredox-catalyzed hydrosulfonylation of arylallenes","hydrosulfonylation","page 1 graphical abstract")
r["notes"]=[F("24 examples; 25–82%; R1 = indol-1-yl, pyrrol-1-yl, phenyl, naphthyl; R2 = aryl, alkyl","24 examples; 25–82%",(1,))]
classify_hydrosulf(r,(1,)); reactions.append(r)
for i,desc in enumerate([
 "Lei, 2017: photoredox hydrosulfonylation of terminal alkynes; uncommon Markovnikov selectivity.",
 "Lei, 2016: allene hydrosulfonylation under O2/pyridine, DCE, 45 °C; 2 examples, 23–34%.",
 "Miao, Ren, 2016: AcOH-mediated tunable hydrosulfonylation of activated cyclopropylidene allenes."] ,1):
    r=base_reaction(f"scheme1_literature_{i}",(2,),"literature_example",[],"Recent example of vinyl sulfone synthesis through hydrosulfonylation",location=f"page 2 Scheme 1, prior-work arrow {i}")
    r["notes"]=[F(desc,desc,(2,))];r["classification"]["coverage"]=F(None,"unresolved",(2,),method="inference",derivation="Prior-work overview was not atom-mapped.");reactions.append(r)
r=base_reaction("scheme1_this_work",(2,),"generic_scheme",[],"Photoredox-catalyzed hydrosulfonylation of allenes","hydrosulfonylation","page 2 Scheme 1, this-work arrow (4)")
r["notes"]=[F("Markovnikov selectivity; regio- and (E)-selective; 24 examples","Markovnikov; regio- and E-selective; 24 examples",(2,))]
classify_hydrosulf(r,(2,));reactions.append(r)

# Table 1: one record per optimization row.
table_rows = {
1:("m_eosinY",3,["DMF","H2O"],"3:1",5,"green LED",40,False),
2:("m_eosinY",5,["DMF","H2O"],"H2O 5 equiv",14,"green LED",62,False),
3:("m_eosinY",5,["DMF"],None,20,"green LED",0,False),
4:("m_eosinY",5,["DMF","H2O"],"H2O 5 equiv",14,"no light",0,False),
5:("m_rose",10,["DMF","H2O"],"H2O 5 equiv",1,"green LED",60,False),
6:("m_rubpy",3,["DMF","H2O"],"H2O 5 equiv",10,"blue LED",52,False),
7:("m_eosinY",5,["EtOH"],None,14,"green LED",35,False),
8:("m_eosinY",5,["DMF","CHCl3"],"3:1",6,"green LED",46,False),
9:("m_eosinY",5,["DMF","TFE"],"1:0.3",4,"green LED",61,False),
10:("m_eosinY",5,["DMF","HFIP"],"2:1",5,"green LED",56,False),
11:("m_eosinY",5,["DMF","TFE"],"1:0.3",4,"green LED",63,True),
12:("m_eosinY",5,["TFE"],None,4,"green LED",38,True),
13:(None,None,["DMF","TFE"],"1:0.3",4,"green LED",0,False),
14:(None,None,["DMF","TFE"],"1:0.3",4,"no light",0,False),
}
solvent_mid={"DMF":"m_dmf","H2O":"m_water","EtOH":"m_etoh","CHCl3":"m_chcl3","TFE":"m_tfe","HFIP":"m_hfip"}
for n,(cat,loading,solvs,ratio,hours,light,yld,base) in table_rows.items():
    r=base_reaction(f"table1_entry_{n:02d}",(2,),categories=["optimization"],name="Hydrosulfonylation optimization",trans="hydrosulfonylation",location=f"page 2 Table 1, entry {n}")
    participant(r,"allene", "m_1a","reactant",(2,),amount=q(0.581,"mmol"))
    participant(r,"sulfinate","m_ptsna","reagent",(2,),amount=q(1.742,"mmol"),equiv=q(3,"equiv"),ref="m_1a")
    participant(r,"product","m_2a","product",(2,))
    if cat: participant(r,"photocatalyst",cat,"catalyst",(2,),molpct=q(loading,"mol%"),ref="m_1a")
    if base: participant(r,"base","m_k2co3","additive",(2,),equiv=q(2,"equiv"),ref="m_1a")
    for j,s in enumerate(solvs): participant(r,f"solvent{j+1}",solvent_mid[s],"solvent",(2,))
    set_solvent(r,(2,),solvs,ratio,concentration=q(0.092,"M"),basis="allene 1a")
    step(r,"irradiation",(2,),light,tuple(p["participant_id"] for p in r["participants"]),time=q(hours,"h"),atmosphere="N2")
    if light!="no light":
        if n==6:set_light(r,(2,),light,q(455,"nm"),q(30,"W"))
        else:set_light(r,(2,),light,q(530,"nm"),"one 30 W LED or two 15 W LEDs")
    else:set_light(r,(2,),"no light")
    r["reaction_scale"]=F("allene 1a (0.581 mmol)",{"reference_molecule_id":"m_1a","amount":q(0.581,"mmol")},(2,))
    outcome(r,"yield","yield",["m_2a"],f"{yld}%",q(yld,"%"),(2,),basis="after column chromatography",yield_type="isolated")
    r["procedure"]["original_text"]=F(f"Table 1 entry {n} and footnote a",f"Table 1 entry {n}; general scale 0.581 mmol under N2",(2,))
    if yld==0: classify_na(r,(2,))
    else: classify_hydrosulf(r,(2,))
    reactions.append(r)

# General scope schemes.
for rid,title,region in [("scheme2_general","Photoredox-catalyzed hydrosulfonylation of N-(1,2-propadien-1-yl)indoles","page 3 Scheme 2 generic arrow"),("scheme3_general","Photoredox-catalyzed hydrosulfonylation of arylallenes","page 3 Scheme 3 generic arrow")]:
    r=base_reaction(rid,(3,),"generic_scheme",[],title,"hydrosulfonylation",region);classify_hydrosulf(r,(3,));reactions.append(r)

# Page 3 Scheme 2, merged with exact experimental charges on pages 5--6 where supplied.
scheme2 = {
 "2b":("1i","m_ptsna",43),"2c":("1f","m_ptsna",34),"2d":("1g","m_phsna",32),
 "2e":("1e","m_ptsna",41),"2f":("1h","m_ptsna",41),"2g":("1g","m_thsna",29),
 "2h":("1c","m_ptsna",42),"2i":("1c","m_phsna",30),"2j":("1c","m_mesna",39),
 "2k":("1e","m_thsna",25),"2m":("1d","m_phsna",27),
}
detail = {
 "2b":dict(page=6,am=(.096,.5681),ey=(.020,.0284),sf=(.303,1.7043),k=(.157,1.1362),dmf=4.745,tfe=1.424,py=42,pm=.078),
 "2c":dict(page=6,am=(.137,.5855),ey=(.020,.029),sf=(.242,1.756),k=(.162,1.1708),dmf=4.890,tfe=1.487,py=43,pm=.094),
 "2d":dict(page=6,am=(.137,.5855),ey=(.020,.029),sf=(.288,1.757),k=(.162,1.1708),dmf=4.890,tfe=1.487,py=32,pm=.070),
 "2e":dict(page=6,am=(.092,.4973),ey=(.017,.025),sf=(.266,1.492),k=(.137,.9946),dmf=4.153,tfe=1.246,py=45,pm=.077),
 "2f":dict(page=6,am=(.090,.500),ey=(.017,.025),sf=(.267,1.500),k=(.138,1.00),dmf=4.176,tfe=1.253,py=41,pm=.069),
 "2g":dict(page=6,am=(.137,.5855),ey=(.020,.029),sf=(.299,1.757),k=(.162,1.1708),dmf=4.890,tfe=1.487,py=29,pm=.065),
 "2h":dict(page=6,am=(.072,.5548),ey=(.019,.0277),sf=(.296,1.6644),k=(.153,1.110),dmf=4.634,tfe=1.390,py=42,pm=.072),
 "2i":dict(page=6,am=(.090,.5807),ey=(.020,.0290),sf=(.286,1.7421),k=(.160,1.1614),dmf=4.850,tfe=1.455,py=30,pm=.051),
 "2j":dict(page=6,am=(.090,.5807),ey=(.020,.0290),sf=(.177,1.7421),k=(.160,1.1614),dmf=4.850,tfe=1.455,py=39,pm=.053),
}
for prod,(sub,sf,yld) in scheme2.items():
    pages=(3,detail[prod]["page"]) if prod in detail else (3,)
    r=base_reaction(f"scope_{prod}",pages,categories=["substrate_scope"],name="Photoredox-catalyzed hydrosulfonylation of N-allenyl heteroarenes",trans="hydrosulfonylation",location=f"page 3 Scheme 2, product {prod}"+("; page 6 product paragraph" if prod in detail else ""))
    d=detail.get(prod)
    participant(r,"allene",f"m_{sub}","reactant",pages,mass=q(d["am"][0],"g") if d else None,amount=q(d["am"][1],"mmol") if d else q(.581,"mmol"),equiv=q(1,"equiv"),ref=f"m_{sub}")
    participant(r,"sulfinate",sf,"reagent",pages,mass=q(d["sf"][0],"g") if d else None,amount=q(d["sf"][1],"mmol") if d else q(1.742,"mmol"),equiv=q(3,"equiv"),ref=f"m_{sub}")
    participant(r,"EY","m_eosinY","catalyst",pages,mass=q(d["ey"][0],"g") if d else None,amount=q(d["ey"][1],"mmol") if d else q(.029,"mmol"),molpct=q(5,"mol%"),ref=f"m_{sub}")
    participant(r,"K2CO3","m_k2co3","additive",pages,mass=q(d["k"][0],"g") if d else None,amount=q(d["k"][1],"mmol") if d else q(1.161,"mmol"),equiv=q(2,"equiv"),ref=f"m_{sub}")
    participant(r,"DMF","m_dmf","solvent",pages,volume=q(d["dmf"],"mL") if d else q(4.850,"mL"))
    participant(r,"TFE","m_tfe","solvent",pages,volume=q(d["tfe"],"mL") if d else q(1.455,"mL"))
    participant(r,"product",f"m_{prod}","product",pages)
    set_solvent(r,pages,["DMF","TFE"],"1:0.3",q((d["dmf"]+d["tfe"]) if d else 6.305,"mL"))
    step(r,"irradiation",pages,"irradiate with green LED and stir; TLC tracking",tuple(p["participant_id"] for p in r["participants"]),temp=q(25,"°C"),time={"minimum":3,"maximum":5,"unit":"h","raw_text":"3–5 h"},atmosphere="N2")
    step(r,"purification",pages,"evaporate solvent, silica-gel chromatography; product crystallized/washed where described",("product",),stage="purification")
    set_light(r,pages,"green LED",q(530,"nm"),"one 30 W LED or two 15 W LEDs")
    r["reaction_scale"]=F(f"allene {sub} ({d['am'][1] if d else .581} mmol)",{"reference_molecule_id":f"m_{sub}","amount":q(d["am"][1] if d else .581,"mmol")},pages)
    if d and d["py"]!=yld:
        alts=[F(f"Scheme 2: {yld}%",q(yld,"%"),(3,)),F(f"product paragraph: {d['py']}%",q(d['py'],"%"),(6,))]
        outcome(r,"yield","yield",[f"m_{prod}"],f"conflicting source yields: {yld}% and {d['py']}%",None,pages,basis="isolated after column chromatography",yield_type="isolated",alternatives=alts)
        r["review"]["unresolved"].append(f"Yield conflict for {prod}: Scheme 2 gives {yld}% whereas the page-6 product paragraph gives {d['py']}%; neither value was selected.")
    else:
        outcome(r,"yield","yield",[f"m_{prod}"],f"{yld}%",q(yld,"%"),pages,basis="isolated after column chromatography",yield_type="isolated")
    if d: outcome(r,"product_mass","product_amount",[f"m_{prod}"],f"{d['pm']} g",q(d["pm"],"g"),pages,basis="reported isolated product mass")
    r["procedure"]["general_procedure_reference"]=F("General Procedure for the Hydrosulfonylation of N-Indolylallenes","page 5 general procedure",pages)
    if prod=="2j":
        r["notes"].append(F("Printed page-6 paragraph identifies the 0.090 g / 0.5807 mmol allene charge as 1a, while Scheme 2 and the product identity correspond to unsubstituted allene 1c.","source identity conflict: paragraph says 1a; scheme/product imply 1c",pages))
        r["review"]["unresolved"].append("Reactant label conflict for 2j remains unresolved; participant follows the Scheme 2/product connectivity, while the printed 1a statement is preserved in notes.")
    classify_hydrosulf(r,pages);reactions.append(r)

# 2l/2l' combined outcome.
r=base_reaction("scope_2l_2lprime",(3,),categories=["substrate_scope"],name="Photoredox-catalyzed hydrosulfonylation giving regioisomers 2l and 2l′",trans="hydrosulfonylation",location="page 3 Scheme 2, combined 2l + 2l′ entry")
for args in [("allene","m_1b","reactant"),("sulfinate","m_ptsna","reagent"),("EY","m_eosinY","catalyst"),("K2CO3","m_k2co3","additive"),("DMF","m_dmf","solvent"),("TFE","m_tfe","solvent"),("major","m_2l","product"),("minor","m_2lprime","product")]:participant(r,*args,(3,))
step(r,"irradiation",(3,),"green LED irradiation; TLC tracking",tuple(p["participant_id"] for p in r["participants"]),temp=q(25,"°C"),time={"minimum":3,"maximum":5,"unit":"h"},atmosphere="N2")
outcome(r,"combined_yield","yield",["m_2l","m_2lprime"],"58%",q(58,"%"),(3,),basis="isolated by column chromatography",yield_type="combined isolated",scope={"from":"1b","to":["2l","2l′"]})
outcome(r,"regioisomer_ratio","selectivity",["m_2l","m_2lprime"],"2l:2l′ = 2.45:1",{"ratio":[2.45,1],"order":["m_2l","m_2lprime"]},(3,),scope={"metric":"regioisomer ratio","order":["2l","2l′"]})
classify_hydrosulf(r,(3,));reactions.append(r)

# Scheme 3 arylallene scope.
scope3 = {
 "4a":("3a","m_ptsna",74),"4b":("3a","m_phsna",78),"4c":("3a","m_brsna",74),
 "4d":("3a","m_nosna","trace"),"4e":("3a","m_thsna",48),"4f":("3a","m_mesna",82),
 "4g":("3c","m_ptsna",65),"4h":("3b","m_ptsna",76),"4i":("3f","m_ptsna",17),
 "4j":("3e","m_ptsna",0),"4k":("3g","m_ptsna",0),"4l":("3d","m_phsna",74),
 "4m":("3h","m_phsna",69),"4n":("3h","m_brsna",59),
}
for prod,(sub,sf,res) in scope3.items():
    r=base_reaction(f"scope_{prod}",(3,),categories=["substrate_scope"],name="Photoredox-catalyzed hydrosulfonylation of arylallenes",trans="hydrosulfonylation",location=f"page 3 Scheme 3, entry {prod}")
    participant(r,"allene",f"m_{sub}","reactant",(3,),amount=q(.603,"mmol"),equiv=q(1,"equiv"),ref=f"m_{sub}")
    participant(r,"sulfinate",sf,"reagent",(3,),amount=q(1.810,"mmol"),equiv=q(3,"equiv"),ref=f"m_{sub}")
    participant(r,"EY","m_eosinY","catalyst",(3,),amount=q(.0302,"mmol"),molpct=q(5,"mol%"),ref=f"m_{sub}")
    participant(r,"DMF","m_dmf","solvent",(3,),volume=q(5.0,"mL"));participant(r,"TFE","m_tfe","solvent",(3,),volume=q(1.5,"mL"))
    if prod in prod4:participant(r,"product",f"m_{prod}","product",(3,))
    set_solvent(r,(3,),["DMF","TFE"],"1:0.3",q(6.5,"mL"))
    step(r,"irradiation",(3,),"green LED irradiation; TLC tracking",tuple(p["participant_id"] for p in r["participants"]),temp=q(25,"°C"),time={"minimum":2,"maximum":4,"unit":"h","raw_text":"2–4 h"},atmosphere="N2")
    set_light(r,(3,),"green LED",q(530,"nm"),"one 30 W LED or two 15 W LEDs")
    r["reaction_scale"]=F("allene 3 (0.603 mmol)",{"reference_molecule_id":f"m_{sub}","amount":q(.603,"mmol")},(3,))
    if res=="trace":
        outcome(r,"trace","trace",[f"m_{prod}"],"trace","trace",(3,),basis="source reports trace; measurement method not stated")
        classify_hydrosulf(r,(3,))
    elif res==0:
        outcome(r,"no_product","no_reaction",[],"0%" if prod=="4k" else "complex inseparable mixture / 0% target","0% target" if prod=="4k" else "complex inseparable mixture",(3,),basis="source scheme/text; exact analytical method not reported")
        if prod=="4k": r["notes"].append(F("Scheme 3 reports 4k, 0%.","0%",(3,)))
        else:r["notes"].append(F("Benzylallene gave a complex inseparable mixture.","complex inseparable mixture",(3,)))
        classify_na(r,(3,))
    else:
        raw=f"{res}%"+(" (81% on 2 mmol scale)" if prod=="4b" else "")
        outcome(r,"yield","yield",[f"m_{prod}"],raw,q(res,"%"),(3,),basis="isolated by column chromatography",yield_type="isolated")
        if prod=="4b": outcome(r,"scaleup_yield","yield",["m_4b"],"81% on 2 mmol scale",q(81,"%"),(3,),basis="isolated by column chromatography",yield_type="isolated",scope={"scale":q(2,"mmol")})
        classify_hydrosulf(r,(3,))
    reactions.append(r)

# Page 4 control experiments.
controls=[
 ("control_d2o_3a","m_3a","m_phsna",["m_ctrl1d"],"5 h",[("yield","m_ctrl1d","52%",52),("yield_brsm","m_ctrl1d","84% after recovering starting allene",84)]),
 ("control_d2o_1m","m_1c","m_ptsna",["m_ctrl2d","m_indole"],"5 h",[("yield","m_ctrl2d","32%",32),("byproduct_yield","m_indole","40%",40)]),
 ("control_acetic_acid","m_3a","m_phsna",["m_4b","m_4bprime"],"2 h",[("yield_4b","m_4b","40%",40),("yield_4bprime","m_4bprime","20%",20)]),
 ("control_tempo","m_3a","m_phsna",[],"5 h",[]),
 ("control_oxygen","m_1a","m_ptsna",["m_2a"],"12 h",[("trace","m_2a","trace",None)]),
]
for idx,(rid,sub,sf,products,tm,outs) in enumerate(controls,1):
    r=base_reaction(rid,(4,),categories=["control","mechanistic_study"],name="Control experiment for allene hydrosulfonylation",trans="hydrosulfonylation",location=f"page 4 Scheme 4, equation {idx}")
    participant(r,"allene",sub,"reactant",(4,),amount=q(.431,"mmol"),equiv=q(1,"equiv"),ref=sub)
    participant(r,"sulfinate",sf,"reagent",(4,),amount=q(1.293,"mmol"),equiv=q(3,"equiv"),ref=sub)
    participant(r,"EY","m_eosinY","catalyst",(4,),amount=q(.022,"mmol"),molpct=q(5,"mol%"),ref=sub)
    participant(r,"DMF","m_dmf","solvent",(4,),volume=q(3.6,"mL"))
    if idx in (1,2):participant(r,"D2O","m_d2o","solvent",(4,),volume=q(1.6,"mL"))
    else:participant(r,"TFE","m_tfe","solvent",(4,),volume=q(1.6,"mL"))
    if idx==2:participant(r,"K2CO3","m_k2co3","additive",(4,),equiv=q(2,"equiv"),ref=sub)
    if idx==3:participant(r,"AcOH","m_hoac","additive",(4,))
    if idx==4:participant(r,"TEMPO","m_tempo","additive",(4,),equiv=q(3,"equiv"),ref=sub)
    if idx==5:participant(r,"oxygen","m_o2","reagent",(4,))
    for j,p in enumerate(products):participant(r,f"product{j+1}",p,"product",(4,))
    set_solvent(r,(4,),["DMF","D2O" if idx in (1,2) else ("DMF","TFE")][0] if False else ["DMF","D2O" if idx in (1,2) else "TFE"],total=q(5.2,"mL"))
    step(r,"irradiation",(4,),"green LED irradiation",tuple(p["participant_id"] for p in r["participants"]),temp=q(25,"°C"),time=q(float(tm.split()[0]),"h"),atmosphere="O2" if idx==5 else "N2")
    set_light(r,(4,),"green LED",q(530,"nm"),"one 30 W LED or two 12 W LEDs")
    for oid,target,raw,num in outs:
        if oid=="trace":outcome(r,oid,"trace",[target],raw,"trace",(4,),basis="trace in inseparable mixture with numerous byproducts")
        else:outcome(r,oid,"yield",[target],raw,q(num,"%"),(4,),basis="isolated by column chromatography",yield_type="isolated" if "brsm" not in oid else "yield after recovering starting allene")
    if idx==4:
        outcome(r,"no_reaction","no_reaction",[],"0%","no reaction",(4,),basis="complete quenching in presence of TEMPO");classify_na(r,(4,))
    else:classify_hydrosulf(r,(4,))
    reactions.append(r)

r=base_reaction("scheme5_proposed_mechanism",(4,),"proposed_mechanism",[],"Proposed reaction pathway","radical-polar crossover hydrosulfonylation","page 4 Scheme 5")
r["notes"]=[F("Excited eosin Y is reductively quenched by sulfinate; sulfonyl radical addition gives allyl radical A, reduction gives anion B, protonation gives 2a′/2a; base-catalyzed isomerization is proposed.","author-proposed mechanism; intermediates A/B not isolated",(4,))]
r["classification"]["coverage"]=F(None,"unresolved",(4,),method="inference",derivation="Proposed elementary steps are not recorded as a verified net experiment.");reactions.append(r)

# Starting allene 1a/1b preparations.
for lab,prec,charges in [
 ("1a","m_propargyl_1a",dict(pm=.500,pa=2.2,tbabm=.071,tbaba=.22,kohm=.121,koha=2.2,tol=27)),
 ("1b","m_propargyl_1b",dict(pm=.300,pa=1.2,tbabm=.039,tbaba=.12,kohm=.069,koha=1.2,tol=16)),
]:
    r=base_reaction(f"prep_{lab}",(4,),categories=["other"],name="Synthesis of 2-substituted N-allenylindole",trans="propargyl-to-allene isomerization",location=f"page 4 experimental paragraph for {lab}")
    participant(r,"precursor",prec,"reactant",(4,),mass=q(charges["pm"],"g"),amount=q(charges["pa"],"mmol"),equiv=q(1,"equiv"),ref=prec)
    participant(r,"TBAB","m_tbab","additive",(4,),mass=q(charges["tbabm"],"g"),amount=q(charges["tbaba"],"mmol"),equiv=q(.1,"equiv"),ref=prec)
    participant(r,"KOH","m_koh","reagent",(4,),mass=q(charges["kohm"],"g"),amount=q(charges["koha"],"mmol"),equiv=q(1,"equiv"),ref=prec)
    participant(r,"toluene","m_toluene","solvent",(4,),volume=q(charges["tol"],"mL"));participant(r,"product",f"m_{lab}","product",(4,))
    set_solvent(r,(4,),["toluene"],total=q(charges["tol"],"mL"))
    step(r,"reaction",(4,),"stir mixture",("precursor","TBAB","KOH","toluene"),temp="rt",time=q(24,"h"))
    step(r,"workup",(4,),"filter and evaporate solvent; crude allene used in further reactions",("product",),stage="workup")
    step(r,"analytical_purification",(4,),"thin silica layer neutralized with Et3N; EtOAc–hexane eluent",("product",),stage="purification")
    outcome(r,"yield","yield",[f"m_{lab}"],"considering the reaction yield to be quantitative",q(100,"%"),(4,),basis="crude allene after filtration and solvent evaporation; analytical sample chromatographed",yield_type="assumed quantitative by authors")
    r["classification"]["coverage"]=F(None,"unresolved",(4,),method="inference",derivation="Propargyl/allene hydrogen and bond-order rearrangement was not atom-mapped.")
    reactions.append(r)

# General alkylation/isomerization procedure and individual products 1c--1i.
r=base_reaction("general_procedure_1c_1i",(4,5),"generic_scheme",[],"General procedure for synthesis of 1-propadienyl indoles or pyrroles 1c–i","N-alkylation/propargyl-allene isomerization","page 4 lower right continuing to page 5 upper left")
participant(r,"heteroarene","m_start_1c","reactant",(4,5),amount=q(1.282,"mmol"),equiv=q(1,"equiv"),ref="m_start_1c")
participant(r,"KOH","m_koh","reagent",(4,5),mass=q(.360,"g"),equiv=q(5,"equiv"),ref="m_start_1c")
participant(r,"DMSO","m_dmso","solvent",(4,5),volume=q(2.5,"mL"))
participant(r,"propargyl_bromide","m_propargyl_br","reagent",(4,5),volume=q(.158,"mL"),equiv=q(1.1,"equiv"),ref="m_start_1c")
step(r,"preactivation",(4,),"stir heteroarene with KOH in DMSO",("heteroarene","KOH","DMSO"),temp="rt",time=q(45,"min"))
step(r,"addition",(4,),"add 80% propargyl bromide solution in toluene slowly; exothermic",("propargyl_bromide",),temp={"minimum":28,"maximum":30,"unit":"°C","raw_text":"28–30 °C"})
step(r,"reaction",(4,),"stir resulting brown mixture; TLC",tuple(p["participant_id"] for p in r["participants"]),temp="rt",time={"minimum":45,"maximum":60,"unit":"min","raw_text":"45–60 min"})
step(r,"workup",(5,),"dilute with water; extract with Et2O or EtOAc; wash organic layer with water to neutral aqueous layer; dry Na2SO4; evaporate",stage="workup")
step(r,"analytical_purification",(5,),"rapid thin-layer silica chromatography neutralized with Et3N; EtOAc–hexane 1:20 or 1:10",stage="purification")
r["classification"]["coverage"]=F(None,"unresolved",(4,5),method="inference",derivation="Generic procedure gives mixtures of propargyl and allenyl products; complete net mapping is unresolved.");reactions.append(r)

prep_yields={"1c":(.185,93),"1d":(.705,45),"1e":(.195,82),"1f":(.261,87),"1g":(.264,88),"1h":(.208,92),"1i":(.206,95)}
for lab,(mass,yld) in prep_yields.items():
    record_pages=(5,) if lab=="1d" else (4,5)
    r=base_reaction(f"prep_{lab}",record_pages,categories=["other"],name="Synthesis of 1-propadienyl heteroarene",trans="N-alkylation/propargyl-allene isomerization",location=f"page 5 product paragraph {lab}; page 4–5 general procedure when applicable")
    participant(r,"heteroarene",f"m_start_{lab}","reactant",(5,))
    participant(r,"propargyl_bromide","m_propargyl_br","reagent",(5,))
    participant(r,"product",f"m_{lab}","product",(5,))
    if lab=="1d":
        r["general_scheme_ids"]=[]
        r["procedure"]["original_text"]=F("It was synthesized according to the literature procedure; 1.034 mL (14.925 mmol) pyrrole was used.","literature procedure; pyrrole 1.034 mL (14.925 mmol)",(5,))
        r["participants"][0]["charges"]["volume"]=F("1.034 mL",q(1.034,"mL"),(5,));r["participants"][0]["charges"]["amount"]=F("14.925 mmol",q(14.925,"mmol"),(5,))
        step(r,"reported_literature_procedure",(5,),"synthesized according to cited literature procedure; exact operations, temperature and time are not supplied",tuple(p["participant_id"] for p in r["participants"]))
        r["review"]["unresolved"].append("The cited preparation procedure for 1d is not part of the supplied source; temperature, time, reagent charges and workup remain unreported.")
    else:
        r["general_scheme_ids"]=[f"{DOC}:general_procedure_1c_1i"]
        r["procedure"]["general_procedure_reference"]=F("General Procedure for 1c–i","page 4–5 general procedure",(5,))
        step(r,"gp_preactivation",(4,),"inherit general procedure: stir heteroarene with KOH in DMSO",("heteroarene",),temp="rt",time=q(45,"min"))
        step(r,"gp_addition",(4,),"inherit general procedure: add 80% propargyl bromide solution slowly; exothermic",("propargyl_bromide",),temp={"minimum":28,"maximum":30,"unit":"°C","raw_text":"28–30 °C"})
        step(r,"gp_reaction",(4,),"inherit general procedure: stir resulting mixture; TLC",tuple(p["participant_id"] for p in r["participants"]),temp="rt",time={"minimum":45,"maximum":60,"unit":"min","raw_text":"45–60 min"})
        step(r,"gp_workup",(5,),"inherit general procedure: water dilution, extraction, water washing to neutral, Na2SO4 drying and evaporation",stage="workup")
        step(r,"gp_analytical_purification",(5,),"inherit general procedure: rapid thin-layer silica chromatography neutralized with Et3N for analytical sample",("product",),stage="purification")
    outcome(r,"product_mass","product_amount",[f"m_{lab}"],f"{mass} g",q(mass,"g"),(5,),basis="reported product mass after the stated extraction/evaporation sequence")
    outcome(r,"yield","yield",[f"m_{lab}"],f"{yld}%",q(yld,"%"),(5,),basis="reported mass/percent yield; analytical sample may have been rapidly chromatographed",yield_type="reported")
    r["classification"]["coverage"]=F(None,"unresolved",(5,),method="inference",derivation="N-alkylation/allene formation was not atom-mapped and the procedure can also form propargyl isomer.")
    reactions.append(r)

# Detailed synthesis of allene 3g. Other 3a--f/h literature-procedure statements are source context below.
r=base_reaction("prep_3g",(5,),categories=["other"],name="Synthesis of 1-(buta-2,3-dien-1-yloxy)-4-(propan-2-yl)benzene",trans="allene synthesis",location="page 5 right column, compound 3g")
participant(r,"alkyne","m_alkyne_3g","reactant",(5,),mass=q(.500,"g"),amount=q(2.873,"mmol"),equiv=q(1,"equiv"),ref="m_alkyne_3g")
participant(r,"paraformaldehyde","m_paraform","reagent",(5,),mass=q(.172,"g"),amount=q(5.873,"mmol"),equiv=q(2,"equiv"),ref="m_alkyne_3g")
participant(r,"CuBr","m_cubr","reagent",(5,),mass=q(.205,"g"),amount=q(1.4365,"mmol"),equiv=q(.5,"equiv"),ref="m_alkyne_3g")
participant(r,"diisopropylamine","m_dipea","reagent",(5,),volume=q(.728,"mL"),amount=q(5.1714,"mmol"),equiv=q(1.8,"equiv"),ref="m_alkyne_3g")
participant(r,"dioxane","m_dioxane","solvent",(5,),volume=q(14,"mL"));participant(r,"product","m_3g","product",(5,))
set_solvent(r,(5,),["1,4-dioxane"],total=q(14,"mL"))
step(r,"addition",(5,),"add diisopropylamine in one portion",("diisopropylamine",))
step(r,"reaction",(5,),"reflux bright-green-yellow mixture",tuple(p["participant_id"] for p in r["participants"]),temp="reflux",time=q(2,"h"))
step(r,"workup",(5,),"cool to rt; dilute with water; extract with Et2O; dry Na2SO4; evaporate",stage="workup")
step(r,"purification",(5,),"flash chromatography on silica gel, EtOAc–hexane 1:10",("product",),stage="purification")
outcome(r,"product_mass","product_amount",["m_3g"],"0.303 g",q(.303,"g"),(5,),basis="pure allene after flash chromatography")
outcome(r,"yield","yield",["m_3g"],"56%",q(56,"%"),(5,),basis="pure allene after flash chromatography",yield_type="isolated")
r["classification"]["coverage"]=F(None,"unresolved",(5,),method="inference",derivation="The alkyne/paraformaldehyde-to-allene transformation was not atom-mapped.");reactions.append(r)

# General hydrosulfonylation procedure on page 5.
r=base_reaction("general_procedure_indolyl_hydrosulfonylation",(5,),"generic_scheme",[],"General Procedure for the Hydrosulfonylation of N-Indolylallenes","hydrosulfonylation","page 5 lower right")
participant(r,"allene","m_1c","reactant",(5,),equiv=q(1,"equiv"),ref="m_1c")
participant(r,"EY","m_eosinY","catalyst",(5,),molpct=q(5,"mol%"),ref="m_1c")
participant(r,"sulfinate","m_ptsna","reagent",(5,),equiv=q(3,"equiv"),ref="m_1c")
participant(r,"K2CO3","m_k2co3","additive",(5,),equiv=q(2,"equiv"),ref="m_1c")
participant(r,"DMF","m_dmf","solvent",(5,));participant(r,"TFE","m_tfe","solvent",(5,))
set_solvent(r,(5,),["DMF","TFE"],"1:0.3")
step(r,"irradiation",(5,),"subject mixture under N2 to green light irradiation and stir for appropriate time (TLC)",tuple(p["participant_id"] for p in r["participants"]),atmosphere="N2")
step(r,"purification",(5,),"evaporate solvent; silica chromatography with EtOAc–hexane; crystallize in hexane with a few drops EtOAc and wash",stage="purification")
set_light(r,(5,),"green LED");classify_hydrosulf(r,(5,));reactions.append(r)

# Exact 2a preparation paragraph starts page 5 and yield continues page 6; distinct scale from Table 1.
r=base_reaction("prep_2a",(5,6),categories=["substrate_scope"],name="Hydrosulfonylation preparation of 2a",trans="hydrosulfonylation",location="page 5 lower right continuing page 6 upper left")
participant(r,"allene","m_1a","reactant",(5,6),mass=q(.070,"g"),amount=q(.3043,"mmol"),equiv=q(1,"equiv"),ref="m_1a")
participant(r,"EY","m_eosinY","catalyst",(5,6),mass=q(.0105,"g"),amount=q(.0152,"mmol"),molpct=q(5,"mol%"),ref="m_1a")
participant(r,"sulfinate","m_ptsna","reagent",(5,6),mass=q(.162,"g"),amount=q(.9129,"mmol"),equiv=q(3,"equiv"),ref="m_1a")
participant(r,"K2CO3","m_k2co3","additive",(5,6),mass=q(.084,"g"),amount=q(.6086,"mmol"),equiv=q(2,"equiv"),ref="m_1a")
participant(r,"DMF","m_dmf","solvent",(5,6),volume=q(2.542,"mL"));participant(r,"TFE","m_tfe","solvent",(5,6),volume=q(.763,"mL"));participant(r,"product","m_2a","product",(5,6))
set_solvent(r,(5,6),["DMF","TFE"],"1:0.3",q(3.305,"mL"))
step(r,"irradiation",(5,6),"general N-indolylallene hydrosulfonylation under green LED; exact time not restated in product paragraph",tuple(p["participant_id"] for p in r["participants"]),atmosphere="N2")
step(r,"purification",(5,6),"evaporate, silica chromatography, crystallize and wash according to general procedure",("product",),stage="purification")
outcome(r,"product_mass","product_amount",["m_2a"],"0.074 g",q(.074,"g"),(6,),basis="isolated white crystals")
outcome(r,"yield","yield",["m_2a"],"63%",q(63,"%"),(6,),basis="isolated after chromatography/crystallization",yield_type="isolated")
r["reaction_scale"]=F("1a, 0.3043 mmol",{"reference_molecule_id":"m_1a","amount":q(.3043,"mmol")},(5,))
r["procedure"]["general_procedure_reference"]=F("General Procedure for the Hydrosulfonylation of N-Indolylallenes","page 5 general procedure",(5,6))
r["notes"].append(F("Product heading says phenylsulfonyl, while the charged reagent and structure are p-toluenesulfinate / Ts-substituted.","source naming/reagent inconsistency preserved",(5,6)))
r["review"]["unresolved"].append("The 2a heading omits the p-methyl substituent in 'phenylsulfonyl'; the reagent charge and source structure indicate Ts. Expert resolution required.")
classify_hydrosulf(r,(5,6));reactions.append(r)

# Source context: title, characterization, literature-only starting-material statements.
bundle["source_context"] = [
 {"context_id":f"{DOC}:title","kind":"bibliographic_metadata","value":"Photoredox-Catalyzed Hydrosulfonylation of Arylallenes","evidence_ids":[E[1],T[1]]},
 {"context_id":f"{DOC}:page1_bioactive_sulfones","kind":"non_synthetic_figure","value":"Figure 1 shows representative biologically active sulfones; no experiment in this paper.","evidence_ids":[E[1],T[1]]},
 {"context_id":f"{DOC}:pages4_6_characterization","kind":"characterization","value":"NMR/HRMS/IR/melting-point blocks are characterization evidence attached to compounds, not standalone synthesis records.","evidence_ids":[E[4],T[4],E[5],T[5],E[6],T[6]]},
 {"context_id":f"{DOC}:page5_literature_allenes","kind":"literature_procedure_reference","value":"Compounds 3a–3f and 3h are stated to have been synthesized according to cited literature procedures; the cited procedures were not supplied and were not promoted to fully specified experiments.","evidence_ids":[E[5],T[5]]},
 {"context_id":f"{DOC}:page3_photophysics","kind":"mechanistic_observation","value":"The article reports sulfinate quenching of excited eosin Y, no fluorescence quenching by phenylallene, and quantum yield 0.88; experimental details are deferred to unavailable SI.","evidence_ids":[E[3],T[3]]},
 {"context_id":f"{DOC}:si_unavailable","kind":"availability_note","value":"The main article points to Supporting Information, but no SI file was explicitly provided/associated in this workspace.","evidence_ids":[E[1],T[1]]},
]

bundle["molecules"] = list(mols.values())
bundle["reactions"] = reactions
st.write(OUT / "submitted.json", bundle)

raw = {
 "document_id":DOC,
 "batch_id":"allene-stereo-v31-200-20261017",
 "pages":[1,2,3,4,5,6],
 "semantic_accuracy_verified":False,
 "visual_review":"All six supplied original-resolution page images were opened and inspected.",
 "record_count":len(reactions),
 "molecule_count":len(mols),
 "decisions":[
  "Optimization Table 1 was represented row-by-row (14 experiments).",
  "Scope drawings and matching product paragraphs were merged when they describe the same labeled experiment; differing source yields remain conflicts.",
  "Scheme 4 controls are specific mechanistic-study experiments; Scheme 5 is a proposed mechanism.",
  "Standalone spectra/HRMS/IR blocks are source_context, not synthetic experiments.",
  "Literature-only procedures for 3a–3f/3h remain source_context because the cited procedures are not supplied.",
  "No SI association was inferred.",
 ],
}
st.write(OUT / "raw_model/pages_0001_0006_raw_extraction.json", raw)

stereo_checks = []
for m in bundle["molecules"]:
    smi = st.value(m["smiles"])
    if smi:
        stereo_checks.append({"molecule_id":m["molecule_id"],"smiles":smi,"standard_stereo":st.stereo_inspect(smi)})
st.write(OUT / "raw_model/pages_0001_0006_standard_stereo_checks.json", {
    "source_verified":False,
    "note":"Output of the same graph diagnostic exposed by `reaction_tools.py standard stereo`; source drawing verification remains pending.",
    "checks":stereo_checks,
})

mentions = st.source_mentions(st.read(OUT / "dataset.json"), ["1a","2a","2b","2c","2e","2j","4b","3g"], context=180)
st.write(OUT / "raw_model/pages_0001_0006_source_mentions.json", mentions)

page_regions = {
 1:[
  ("title, abstract and graphical abstract","True title and a generic 24-example hydrosulfonylation overview; recorded as a generic scheme, not one concrete experiment."),
  ("Figure 1 and introduction","Biologically active sulfones and literature discussion are non-synthetic context; no experiment in this paper."),
 ],
 2:[
  ("Scheme 1","Three cited prior-work arrows are literature_example; arrow (4) is this-work generic scheme."),
  ("Table 1","All 14 optimization rows were represented independently, including four 0%/no-reaction rows and footnotes for light, wavelength and K2CO3."),
  ("Results prose below/alongside Table 1","Narrative duplicates Table 1 and reports 2-phenylindole byproduct plus X-ray assignment of 2a; it was not double-counted."),
 ],
 3:[
  ("Scheme 2","Twelve N-allenyl scope outcomes were inventoried, including the 2l/2l′ mixture and 2m scale-up footnote."),
  ("Scheme 3","Fourteen arylallene scope outcomes were inventoried, including trace 4d, zero/failed 4j and 4k, and 4b scale-up yield."),
  ("lower-column mechanistic discussion","Control experiments are only previewed here and mapped to Scheme 4 records on page 4; fluorescence/quantum-yield observations remain source_context because details are deferred to SI."),
 ],
 4:[
  ("Scheme 4","Five control/mechanistic-study experiments were recorded separately, including D2O, AcOH, TEMPO and O2 experiments."),
  ("Scheme 5","Author-proposed radical-polar crossover pathway recorded as proposed_mechanism, not as a verified experiment."),
  ("Experimental section right column","General information is context; detailed preparations of 1a and 1b and the beginning of the 1c–i general procedure were extracted."),
 ],
 5:[
  ("upper left continuation and 1c–1i paragraphs","Cross-page workup/purification sequence and seven reported starting-allene outcomes were extracted without duplicating operations."),
  ("right-column starting arylallenes","Detailed 3g synthesis is a specific experiment; 3a–3f/3h only cite external literature procedures and remain source_context."),
  ("lower right hydrosulfonylation general procedure and 2a start","General N-indolylallene hydrosulfonylation and exact 2a charges were extracted; the 2a yield continues on page 6."),
 ],
 6:[
  ("2a continuation and 2b–2j product paragraphs","Exact charges, product masses and yields were linked to the matching scope records; 2a is a distinct 0.3043 mmol preparation."),
  ("NMR/HRMS/IR/melting-point blocks","Characterization was retained as source_context and not promoted to extra synthesis records."),
  ("source discrepancies","Scheme/paragraph yield conflicts for 2b, 2c and 2e, the printed 2j reactant label conflict, and sulfonyl naming inconsistencies were preserved for review."),
 ],
}
coverage={"document_id":DOC,"semantic_accuracy_verified":False,"pages":[]}
for p in range(1,7):
    ids=[r["reaction_id"] for r in reactions if E[p] in r["evidence_ids"]]
    unresolved=[]
    if p==1:unresolved=["Graphical abstract is generic and does not identify all 24 individual experiments."]
    if p==2:unresolved=["Exact temperature is not reported in Table 1; byproduct quantity for 2-phenylindole is not reported."]
    if p==3:unresolved=["Compact drawings remain machine structure candidates; SI-only photophysical methods are unavailable."]
    if p==4:unresolved=["D2O isotope-aware product structures and exact measurement bases require expert review; mechanism is proposed only."]
    if p==5:unresolved=["Cited procedures for 3a–3f/3h were not supplied; individual reaction times for product 2a are inherited only from the general procedure and not restated."]
    if p==6:unresolved=["Yield conflicts 2b/2c/2e and reactant-label conflict for 2j are intentionally unresolved; product headings inconsistently say phenylsulfonyl where Ts reagent/structure is shown."]
    coverage["pages"].append({"page":p,"state":"attempted","regions":[{"region":a,"finding":b,"evidence_ids":[E[p],T[p]]} for a,b in page_regions[p]],"record_ids":ids,"unresolved":unresolved})
st.write(OUT / "coverage/pages_0001_0006.json", coverage)

tasks = st.read(OUT / "tasks.json")
for task in tasks:
    if task.get("document_id") == DOC and task.get("evidence_id") in {T[p] for p in range(1,7)}:
        task["state"] = "completed"
        task["reason"] = "Original page image visually attempted and mapped in coverage/pages_0001_0006.json; completion is not chemical certification."
st.write(OUT / "tasks.json", tasks)

summary = f"""# Extraction summary — Photoredox-Catalyzed Hydrosulfonylation of Arylallenes

- Source: J. Org. Chem. 2020, 85, 2250–2259; DOI 10.1021/acs.joc.9b02960.
- Visually attempted pages: original PDF pages 1–6 only.
- Formal v3.1 output: {len(mols)} molecules/materials and {len(reactions)} reaction records; every machine record remains `needs_review` and `partial`.
- Table 1: 14 independent optimization entries. Scheme 2/3 scope, Scheme 4 controls, Scheme 5 proposed mechanism, and experimental preparations on pages 4–6 are separately typed.

## Reagent-charge checks

- Common N-indolyl scope: allene 0.581 mmol (1 equiv), sulfinate 1.742 mmol (3 equiv), eosin Y 0.029 mmol (5 mol%), K2CO3 1.161 mmol (2 equiv), DMF 4.850 mL, TFE 1.455 mL; 3–5 h, 25 °C, N2, green LED.
- Common arylallene scope: allene 0.603 mmol, sulfinate 1.810 mmol (3 equiv), eosin Y 0.0302 mmol (5 mol%), DMF 5.0 mL, TFE 1.5 mL; 2–4 h, 25 °C, N2, green LED.
- 2a paragraph: 1a 0.070 g/0.3043 mmol; eosin Y 0.0105 g/0.0152 mmol; sodium p-toluenesulfinate 0.162 g/0.9129 mmol; K2CO3 0.084 g/0.6086 mmol; DMF 2.542 mL; TFE 0.763 mL; product 0.074 g, 63%.
- 1a preparation: propargyl precursor 0.500 g/2.2 mmol, TBAB 0.071 g/0.22 mmol, KOH 0.121 g/2.2 mmol, toluene 27 mL, rt 24 h.
- 1b preparation: propargyl precursor 0.300 g/1.2 mmol, TBAB 0.039 g/0.12 mmol, KOH 0.069 g/1.2 mmol, toluene 16 mL, rt 24 h.
- 1c–i general procedure: heteroarene 1.282 mmol, DMSO 2.5 mL, KOH 0.360 g/5 equiv, then propargyl bromide 0.158 mL (80% in toluene, 1.1 equiv) at 28–30 °C, followed by 45–60 min at rt.
- 3g: alkyne 0.500 g/2.873 mmol, paraformaldehyde 0.172 g/5.873 mmol, CuBr 0.205 g/1.4365 mmol, diisopropylamine 0.728 mL/5.1714 mmol, 1,4-dioxane 14 mL; reflux 2 h; 0.303 g, 56%.

## Preserved conflicts and limitations

- 2b: Scheme 2 gives 43%, product paragraph gives 42%; 2c: 34% vs 43%; 2e: 41% vs 45%. Both values and their page evidence are retained with no selected normalized yield.
- The printed 2j paragraph calls the 0.090 g/0.5807 mmol allene `1a`, while Scheme 2 and product connectivity correspond to `1c`; this remains unresolved.
- Several product headings say “phenylsulfonyl” while the charged reagent/source drawing is p-toluenesulfinate/Ts; the inconsistency is preserved.
- No Supporting Information file was provided or associated. Cited external procedures were not opened. Candidate SMILES and stereo diagnostics are not chemical certification.
"""
(OUT / "extraction_summary.md").write_text(summary)
print(json.dumps({"molecules":len(mols),"reactions":len(reactions)},ensure_ascii=False))
