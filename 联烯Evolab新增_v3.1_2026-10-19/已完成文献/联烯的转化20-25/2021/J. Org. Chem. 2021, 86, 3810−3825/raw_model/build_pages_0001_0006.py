#!/usr/bin/env python3
"""Build evidence-linked v3.1 records for source PDF pages 1--6 only."""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as sr

RUN = ROOT / "runs/current"
DATASET = RUN / "dataset.json"
SUBMITTED = RUN / "submitted.json"
DOC_ID = "doc_0341d857f4f70e0f"
TITLE = (
    "Palladium-Catalyzed Synthesis of 1-Vinyltetrahydro-β-carbolines and "
    "Aza-spiroindolenines: Access to the Syntheses of 1-Vinyl-β-carbolines "
    "and Eudistomins Y1 and Y2"
)


def dump(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


dataset = json.loads(DATASET.read_text())
doc = next(d for d in dataset["documents"] if d["document_id"] == DOC_ID)
doc["title"] = TITLE
DATASET.write_text(json.dumps(dataset, ensure_ascii=False, indent=2) + "\n")
bundle = sr.empty(dataset["documents"])


def ev(page: int) -> list[str]:
    return [f"{DOC_ID}:page{page}:text", f"{DOC_ID}:page{page}:image"]


def f(raw, normalized=None, evidence=(), method="original_extraction", derivation=None,
      status=None):
    if normalized is None and raw is not None:
        normalized = raw
    return sr.fact(raw=raw, normalized=normalized, evidence=evidence, method=method,
                   derivation=derivation, status=status)


molecules: dict[str, dict] = {}
diagnostics: dict[str, dict] = {}


def molecule(key: str, label: str, page: int, name: str | None = None,
             smiles: str | None = None, chirality="unresolved", geometry="unresolved",
             representation="unresolved", assignment=None, mixture=None,
             note: str | None = None) -> str:
    mid = f"{DOC_ID}:mol:{key}"
    evidence = ev(page)
    if mid in molecules:
        return mid
    m = sr.new_molecule(mid, DOC_ID, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    m["names"]["original"] = f(name or label, name or label, evidence)
    m["identity_note"] = f(note, note, evidence) if note else sr.fact(evidence=evidence)
    if smiles:
        diag = sr.stereo_inspect(smiles)
        diagnostics[mid] = diag
        m["representation_type"] = "smiles"
        m["smiles"] = f(
            {"compound_label": label, "source": "structure drawing/name on inspected page"},
            diag["canonical_isomeric_smiles"], evidence, "structure_derivation",
            "Machine structure candidate transcribed from the source label/drawing; parsing and "
            "canonicalization do not establish source correctness."
        )
        m["stereochemistry"]["isomeric_smiles"] = f(
            diag["canonical_isomeric_smiles"], diag["canonical_isomeric_smiles"], evidence,
            "structure_derivation", "RDKit canonical isomeric SMILES from the machine candidate."
        )
        m["stereochemistry"]["identity_key"] = f(
            diag["stereo_graph_key"], diag["stereo_graph_key"], evidence,
            "structure_derivation", "Stereo graph key calculated from the machine candidate."
        )
    stereo_reason = (
        "Status assigned after inspecting the page drawing, labels, wedge/dash bonds, alkene "
        "geometry, and reported stereochemical text; unresolved information was not converted "
        "to achiral or racemic."
    )
    for field, value in (("chirality_status", chirality), ("geometry_status", geometry),
                         ("representation_status", representation)):
        m["stereochemistry"][field] = f(None, value, evidence, "inference", stereo_reason)
    if assignment is not None:
        m["stereochemistry"]["assignment"] = f(assignment, assignment, evidence)
    if mixture is not None:
        m["stereochemistry"]["mixture"] = f(mixture, mixture, evidence)
    molecules[mid] = m
    return mid


# Exact, non-stereogenic small structures are machine candidates only.
SMILES = {
    "7a": "Ic1ccccc1", "7b": "COC(=O)c1ccc(I)cc1", "7c": "Ic1cccc2ccccc12",
    "7d": "Ic1cccs1", "7e": "FC(F)(F)c1ccc(I)cc1", "7f": "Fc1ccc(I)cc1",
    "7g": "Cc1ccc(I)cc1", "7h": "COc1ccc(I)cc1", "7i": "CC(=O)c1ccc(I)cc1",
    "DMF": "CN(C)C=O", "DMSO": "CS(C)=O", "DMA": "CC(=O)N(C)C", "MeCN": "CC#N",
    "toluene": "Cc1ccccc1", "dioxane": "O1CCOCC1", "THF": "C1CCOC1",
    "K2CO3": "[K+].[K+].[O-]C([O-])=O", "Ag2CO3": "[Ag+].[Ag+].[O-]C([O-])=O",
    "KOtBu": "[K+].[O-]C(C)(C)C", "KOH": "[K+].[OH-]", "NaH": "[Na+].[H-]",
    "PPh3": "P(c1ccccc1)(c1ccccc1)c1ccccc1", "CCl4": "ClC(Cl)(Cl)Cl",
    "water": "O", "AlCl3": "Cl[Al](Cl)Cl", "pTsOH": "Cc1ccc(S(=O)(=O)O)cc1.O",
}


def named(key, page, name=None, **kw):
    return molecule(key, key, page, name or key, SMILES.get(key),
                    "achiral" if key in SMILES else kw.pop("chirality", "unresolved"),
                    "not_applicable" if key in SMILES else kw.pop("geometry", "unresolved"),
                    "not_applicable" if key in SMILES else kw.pop("representation", "unresolved"),
                    **kw)


# Recurrent source-labelled substrates/products. Complex drawings remain unresolved candidates.
for key, name in {
    "6a": "allenamide 6a", "6b": "allenamide 6b", "6c": "allenamide 6c",
    "6d": "allenamide 6d", "6e": "allenamide 6e", "6f": "allenamide 6f",
    "6g": "allenamide 6g", "6h": "allenamide 6h", "6i": "allenamide 6i",
}.items():
    ch = "unspecified" if key in {"6e", "6f", "6g"} else "achiral"
    molecule(key, key, 4 if key != "6a" else 3, name, chirality=ch,
             geometry="not_applicable", representation="not_applicable")
for key, name in {
    "1aa": "THBC 1aa", "1a": "THBC 1a", "8": "rearranged product 8",
    **{x: f"THBC {x}" for x in ["1ba","1bb","1ac","1ad","1ab","1ae","1af","1ag",
                                      "1ah","1ca","1ch","1db","1de","1ia","1eb","1fb",
                                      "1gb","1ga","1ha","1hb","1hh"]},
}.items():
    special = key in {"1eb", "1fb", "1gb", "1ga"}
    molecule(key, key, 4 if key != "1aa" else 3, name,
             chirality="mixture" if special else "unspecified",
             geometry="not_applicable", representation="not_encoded" if special else "unresolved",
             assignment={"relative_configuration": "C1/C4 cis", "absolute_configuration": "unreported"}
             if special else None,
             mixture={"metric": "dr", "reported": ">99:1", "major_absolute_configuration": "unreported"}
             if special else None)
for key in ["2a","2b","2c","2d","2e","2f","2g","2h","10a","10b","11a","11b"]:
    molecule(key, key, 5 if key.startswith("2") else 6, key,
             chirality="unresolved", geometry="unresolved", representation="unresolved")
for key in ["9a", "9b", "9i", "9g"]:
    ratio = {"9a":"1.5:1","9b":"9:1","9i":"9:1","9g":"2.2:1"}[key]
    assignment = ({"major_relative_configuration": "3R,11S (relative assignment as stated by authors)",
                   "basis": "NOESY correlations Ha/He; absolute configuration not established"}
                  if key == "9b" else None)
    molecule(key, key, 6, f"aza-spiroindolenine {key}", chirality="mixture",
             geometry="not_applicable", representation="not_encoded", assignment=assignment,
             mixture={"metric":"dr", "reported":ratio, "major_absolute_configuration":"unreported"})
for key, name in {"7a":"iodobenzene", "7b":"methyl 4-iodobenzoate",
                  "7c":"1-iodonaphthalene", "7d":"2-iodothiophene",
                  "7e":"4-iodobenzotrifluoride", "7f":"1-fluoro-4-iodobenzene",
                  "7g":"4-iodotoluene", "7h":"4-iodoanisole",
                  "7i":"4-iodoacetophenone"}.items():
    named(key, 4 if key not in {"7i"} else 6, name)


def participant(r, mid, role, pid, evidence, charge=None, reference=None):
    p = sr.new_participant(pid, mid, role, evidence)
    if charge:
        field, raw, norm = charge
        p["charges"][field] = f(raw, norm, evidence)
        p["equivalence_reference"] = f(reference, reference, evidence) if reference else sr.fact(evidence=evidence)
    r["participants"].append(p)
    return p


def reaction(suffix, page, record_type="specific_experiment", categories=None,
             name=None, transformation=None, extra_evidence=()):
    evidence = ev(page) + list(extra_evidence)
    r = sr.new_reaction(f"{DOC_ID}:{suffix}", DOC_ID, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, evidence) if categories is not None else sr.fact(evidence=evidence)
    r["location"] = f(f"PDF page {page}", {"pdf_page":page}, evidence)
    r["reaction_name"] = f(name, name, evidence) if name else sr.fact(evidence=evidence)
    r["transformation_type"] = f(transformation, transformation, evidence) if transformation else sr.fact(evidence=evidence)
    r["review"]["unresolved"] = ["Structures and chemistry require expert source review; machine records remain needs_review."]
    r["completeness"] = "partial"
    bundle["reactions"].append(r)
    return r


def add_step(r, page, temp=None, time=None, atmosphere=None, operation="stirred/heated under reported conditions"):
    evidence = ev(page)
    s = sr.new_step("reaction", evidence)
    s["participant_ids"] = [p["participant_id"] for p in r["participants"]]
    s["operation"] = f(operation, operation, evidence)
    if temp is not None:
        s["temperature"] = f(temp if isinstance(temp, str) else f"{temp} °C",
                               temp if isinstance(temp, str) else {"number":temp,"unit":"°C"}, evidence)
    if time is not None:
        if isinstance(time, tuple):
            raw = f"{time[0]} {time[1]}"; norm = {"number":time[0],"unit":time[1]}
        else:
            raw = norm = time
        s["time"] = f(raw, norm, evidence)
    if atmosphere:
        s["atmosphere"] = f(atmosphere, atmosphere, evidence)
    r["steps"] = [s]


def add_yield(r, page, mid, number, isolated=True, basis_missing=False, oid="yield"):
    evidence = ev(page)
    o = sr.new_outcome(oid, "yield", [mid], evidence)
    o["product_role"] = f("product", "product", evidence)
    o["value"] = f(f"{number}%", {"number":number,"unit":"%"}, evidence)
    if isolated:
        o["yield_type"] = f("isolated pure product", "isolated", evidence)
        o["measurement_basis"] = f("yield of isolated pure product", "isolated pure product", evidence)
    elif basis_missing:
        r["review"]["unresolved"].append(f"Measurement basis for {oid} is not reported on inspected page {page}.")
    r["outcomes"].append(o)


def add_selectivity(r, page, mid, ratio, method=None, assignment="major/minor configurations unreported"):
    evidence = ev(page)
    o = sr.new_outcome("dr", "selectivity", [mid], evidence)
    o["value"] = f(f"dr {ratio}", {"metric":"dr", "ratio":ratio}, evidence)
    o["measurement_basis"] = f(method, method, evidence) if method else sr.fact(evidence=evidence)
    o["scope"] = f({"metric":"dr", "isomer_order":assignment},
                   {"metric":"dr", "isomer_order":assignment}, evidence)
    r["outcomes"].append(o)
    if not method:
        r["review"]["unresolved"].append("The dr measurement method and absolute identities of major/minor isomers are unreported on pages 1–6.")


def classify(r, page, mode):
    evidence = ev(page)
    modes = {
        "cc": ([{"change_id":"b1","atoms":[{"atom_map":None,"element":"C"},{"atom_map":None,"element":"C"}],
                 "before_order":0,"after_order":1,"change_type":"formed"}], [{"category_number":2,"code":"C-C","supporting_change_ids":["b1"]}]),
        "aromatize": ([
            {"change_id":"b1","atoms":[{"atom_map":None,"element":"C"},{"atom_map":None,"element":"C"}],"before_order":1,"after_order":2,"change_type":"order_increased"},
            {"change_id":"b2","atoms":[{"atom_map":None,"element":"C"},{"atom_map":None,"element":"N"}],"before_order":1,"after_order":2,"change_type":"order_increased"}],
            [{"category_number":13,"code":"DB-C-C","supporting_change_ids":["b1"]},{"category_number":16,"code":"DB-C-N","supporting_change_ids":["b2"]}]),
        "carbonyl": ([{"change_id":"b1","atoms":[{"atom_map":None,"element":"C"},{"atom_map":None,"element":"O"}],"before_order":0,"after_order":2,"change_type":"formed"}],
                       [{"category_number":14,"code":"DB-C-O","supporting_change_ids":["b1"]}]),
        "nh": ([{"change_id":"b1","atoms":[{"atom_map":None,"element":"N"},{"atom_map":None,"element":"H"}],"before_order":0,"after_order":1,"change_type":"hydrogen_gain","hydrogen_count":1}],
                [{"category_number":29,"code":"N-H","supporting_change_ids":["b1"]}]),
    }
    if mode is None:
        r["classification"]["labels"] = f([], [], evidence, "inference", "No formed-bond classification asserted because the arrow is a literature summary, generic scheme, failed experiment, or structurally unresolved branch.")
        r["classification"]["coverage"] = f(None, "unresolved", evidence, "inference", "Insufficient inspected-page structure detail for safe bond mapping.")
        r["classification"]["limitations"] = ["No atom mapping; classification unresolved."]
        return
    changes, labels = modes[mode]
    r["classification"]["bond_changes"] = [f(c, c, evidence, "structure_derivation", "Inferred from inspected reactant/product drawings without atom mapping.") for c in changes]
    r["classification"]["labels"] = f(labels, labels, evidence, "structure_derivation", "Taxonomy labels inferred from the depicted net transformation; mapping remains incomplete.")
    r["classification"]["coverage"] = f(None, "partial", evidence, "inference", "Unmapped source drawings; only clearly depicted formed/promoted bonds are classified.")
    r["classification"]["limitations"] = ["No atom-mapped reaction SMILES; bond changes are partial machine inferences.", "Implicit-H and leaving-group changes not shown in the inspected figure remain unresolved unless explicitly represented."]


def procedure(r, page, text):
    r["procedure"]["original_text"] = f(text, text, ev(page))


# Page 2: prior-work arrows and the authors' generic route. Page 1 is duplicate graphical summary evidence.
for tag, citation in [("a","Shi, 2013; ref 5a"),("b","Zhang, 2017; ref 5b"),
                      ("c","Maestri, 2018; ref 5c"),("d","Haak/Baire, 2018; refs 6c-d"),
                      ("e","Greaney, 2011; ref 6e")]:
    r = reaction(f"scheme1_literature_{tag}", 2, "literature_example", name=f"Scheme 1{tag} prior work")
    procedure(r, 2, citation)
    classify(r, 2, None)
r = reaction("scheme1_present_work", 2, "generic_scheme", name="Present-work generic route",
             transformation="allenamide/aryl iodide cyclocondensation followed by aromatization",
             extra_evidence=ev(1))
procedure(r, 2, "6 + 7 → THBC 1; elimination/oxidation → β-carboline 2")
classify(r, 2, None)


# Table 1, page 3: 14 independent optimization experiments.
table1 = [
    ("Pd(OAc)2","PPh3","K2CO3","DMF",80,(26,"h"),0),
    ("Pd(OAc)2","PPh3","Ag2CO3","DMF",80,(24,"h"),15),
    ("Pd2(dba)3","PPh3","Ag2CO3","DMF",80,(1,"h"),26),
    ("Pd2(dba)3","t-Bu-XPhos","Ag2CO3","DMF",80,(1,"h"),71),
    ("Pd2(dba)3","t-Bu-XantPhos","Ag2CO3","DMF",80,(1,"h"),75),
    ("Pd2(dba)3","t-Bu-XantPhos","Ag2CO3","DMF",70,(3,"h"),70),
    ("Pd2(dba)3","Xantphos","Ag2CO3","DMF",80,(4,"h"),30),
    ("Pd2(dba)3","N-XantPhos","Ag2CO3","DMF",80,(5,"h"),22),
    ("Pd(dba)2","t-Bu-XantPhos","Ag2CO3","DMF",80,(1,"h"),20),
    ("Pd(PPh3)4","t-Bu-XantPhos","Ag2CO3","DMF",80,(10,"h"),40),
    ("Pd2(dba)3","t-Bu-XantPhos","Ag2CO3","DMSO",80,(1,"h"),40),
    ("Pd2(dba)3","t-Bu-XantPhos","Ag2CO3","DMA",80,(2,"h"),29),
    ("Pd2(dba)3","t-Bu-XantPhos","Ag2CO3","MeCN",80,(3,"h"),38),
    ("Pd2(dba)3","t-Bu-XantPhos","Ag2CO3","toluene",80,(1,"h"),35),
]
for i,(cat,lig,base,solv,temp,time,yld) in enumerate(table1,1):
    r = reaction(f"table1_entry_{i:02d}", 3, categories=["optimization"], name=f"Table 1 entry {i}", transformation="Pd-catalyzed cyclocondensation")
    participant(r, molecules[f"{DOC_ID}:mol:6a"]["molecule_id"], "reactant", "p_6a", ev(3), ("amount","0.10 mmol",{"number":0.10,"unit":"mmol"}))
    participant(r, named("7a",4,"iodobenzene"), "reactant", "p_7a", ev(3))
    for key,role,charge in [(cat,"catalyst",("mol_percent","5 mol%",{"number":5,"unit":"mol%"})),
                            (lig,"ligand",("mol_percent","12 mol%",{"number":12,"unit":"mol%"})),
                            (base,"reagent",("amount","0.2 mmol",{"number":0.2,"unit":"mmol"})),
                            (solv,"solvent",None)]:
        mid = named(key,3,key)
        participant(r,mid,role,"p_"+key.replace("(","").replace(")","").replace("-", "_"),ev(3),charge, molecules[f"{DOC_ID}:mol:6a"]["molecule_id"] if charge and charge[0]=="mol_percent" else None)
    participant(r, molecules[f"{DOC_ID}:mol:1aa"]["molecule_id"], "product", "p_1aa", ev(3))
    r["reaction_scale"] = f("6a (0.10 mmol)", {"reference_molecule_id":f"{DOC_ID}:mol:6a","amount":{"number":0.10,"unit":"mmol"}}, ev(3))
    r["solvent_conditions"]["composition"] = f(solv, [solv], ev(3))
    r["solvent_conditions"]["concentration"] = f("0.03 M", {"number":0.03,"unit":"M"}, ev(3))
    r["solvent_conditions"]["concentration_basis"] = f("6a", f"{DOC_ID}:mol:6a", ev(3))
    add_step(r,3,temp,time)
    add_yield(r,3,f"{DOC_ID}:mol:1aa",yld)
    procedure(r,3,"6a (0.10 mmol), Pd catalyst (5 mol%), ligand (12 mol%), base (0.2 mmol), reported solvent (0.03 M), reported temperature/time")
    classify(r,3,"cc")


# Table 2 scope and Scheme 2, page 4.
table2 = [
    ("1aa","6a","7a",75,None),("1ba","6b","7a",57,None),("1bb","6b","7b",62,None),
    ("1ac","6a","7c",60,None),("1ad","6a","7d",55,None),("1ab","6a","7b",95,None),
    ("1ae","6a","7e",78,None),("1af","6a","7f",65,None),("1ag","6a","7g",50,None),
    ("1ah","6a","7h",48,None),("1ca","6c","7a",60,None),("1ch","6c","7h",47,None),
    ("1db","6d","7b",45,None),("1de","6d","7e",37,None),("1ia","6i","7a",76,None),
    ("1eb","6e","7b",64,">99:1"),("1fb","6f","7b",70,">99:1"),
    ("1gb","6g","7b",61,">99:1"),("1ga","6g","7a",75,">99:1"),
    ("1ha","6h","7a",40,None),("1hb","6h","7b",45,None),("1hh","6h","7h",35,None),
]
for prod,sub,iod,yld,dr in table2:
    r = reaction(f"table2_{prod}",4,categories=["substrate_scope"],name=f"Table 2 product {prod}", transformation="Pd-catalyzed cyclocondensation")
    ref=f"{DOC_ID}:mol:{sub}"
    participant(r,ref,"reactant",f"p_{sub}",ev(4),("equivalents","1.0 equiv",{"number":1.0,"unit":"equiv"}),ref)
    participant(r,f"{DOC_ID}:mol:{iod}","reactant",f"p_{iod}",ev(4),("equivalents","1.3 equiv",{"number":1.3,"unit":"equiv"}),ref)
    for key,role,n in [("Pd2(dba)3","catalyst",0.05),("t-Bu-XantPhos","ligand",0.12),("Ag2CO3","reagent",2.0)]:
        participant(r,named(key,4,key),role,f"p_{key}",ev(4),("equivalents",f"{n} equiv",{"number":n,"unit":"equiv"}),ref)
    participant(r,named("DMF",4,"DMF"),"solvent","p_DMF",ev(4))
    participant(r,f"{DOC_ID}:mol:{prod}","product",f"p_{prod}",ev(4))
    r["reaction_scale"] = f(f"{sub} (1.0 equiv)", {"reference_molecule_id":ref,"amount":{"number":1.0,"unit":"equiv"}}, ev(4))
    r["solvent_conditions"]["composition"] = f("DMF",["DMF"],ev(4))
    add_step(r,4,80,(1,"h"),"argon")
    add_yield(r,4,f"{DOC_ID}:mol:{prod}",yld)
    if dr:add_selectivity(r,4,f"{DOC_ID}:mol:{prod}",dr)
    procedure(r,4,"Allenamide 6 (1.0 equiv), iodide 7 (1.3 equiv), Pd2(dba)3 (0.05 equiv), t-Bu-XantPhos (0.12 equiv), Ag2CO3 (2.0 equiv), DMF, argon, 80 °C, 1 h")
    classify(r,4,"cc")
r = reaction("scheme2_6a_to_1a",4,categories=["other"],name="Scheme 2 acid-catalyzed cyclization",transformation="intramolecular cyclization")
participant(r,f"{DOC_ID}:mol:6a","reactant","p_6a",ev(4),("equivalents","1.0 equiv",{"number":1.0,"unit":"equiv"}),f"{DOC_ID}:mol:6a")
participant(r,named("pTsOH",4,"p-TsOH·H2O"),"catalyst","p_pTsOH",ev(4),("mol_percent","10 mol%",{"number":10,"unit":"mol%"}),f"{DOC_ID}:mol:6a")
participant(r,named("THF",4,"dry THF"),"solvent","p_THF",ev(4));participant(r,f"{DOC_ID}:mol:1a","product","p_1a",ev(4))
add_step(r,4,"rt","overnight","argon");add_yield(r,4,f"{DOC_ID}:mol:1a",43);procedure(r,4,"6a (1.0 equiv), p-TsOH·H2O (10 mol%), dry THF, argon, rt, overnight");classify(r,4,"cc")


# Table 3, page 5.
table3 = [
    ("KOtBu","dioxane",110,(24,"h"),0,0,"oxygen balloon",2.0),
    ("KOtBu","dioxane",110,(5,"h"),5,40,"open air",2.0),
    ("KOtBu","dioxane",110,(5,"h"),0,86,"argon",2.0),
    ("KOtBu","DMF","rt",(5,"min"),None,10,"argon",2.0),
    ("KOH","DMSO","rt",(2,"h"),0,20,"argon",2.0),
    ("NaH","DMF","rt",(5,"h"),0,60,"argon",2.0),
    ("DDQ","toluene",100,(6,"h"),80,0,"argon",3.0),
    ("DDQ","toluene",90,(12,"h"),80,0,"argon",3.0),
    ("DDQ","toluene",100,(12,"h"),80,0,"argon",2.0),
    ("DDQ","DMF",100,(12,"h"),5,0,"argon",3.0),
    ("DDQ","DMSO",100,(12,"h"),5,0,"argon",3.0),
    ("DDQ","MeCN",100,(12,"h"),40,0,"argon",3.0),
    ("DDQ","dioxane",100,(12,"h"),35,0,"argon",3.0),
]
for i,(agent,solv,temp,time,y2,y8,atm,eqv) in enumerate(table3,1):
    r=reaction(f"table3_entry_{i:02d}",5,categories=["optimization"],name=f"Table 3 entry {i}",transformation="THBC detosylation/aromatization or rearrangement")
    ref=f"{DOC_ID}:mol:1aa"
    participant(r,ref,"reactant","p_1aa",ev(5),("equivalents","1.0 equiv",{"number":1.0,"unit":"equiv"}),ref)
    participant(r,named(agent,5,agent),"reagent",f"p_{agent}",ev(5),("equivalents",f"{eqv} equiv",{"number":eqv,"unit":"equiv"}),ref)
    participant(r,named(solv,5,solv),"solvent",f"p_{solv}",ev(5))
    if y2 is not None:participant(r,f"{DOC_ID}:mol:2a","product","p_2a",ev(5))
    if y8 is not None:participant(r,f"{DOC_ID}:mol:8","product","p_8",ev(5))
    add_step(r,5,temp,time,atm)
    if y2 is not None:add_yield(r,5,f"{DOC_ID}:mol:2a",y2,oid="yield_2a")
    else:r["review"]["unresolved"].append("Table 3 entry 4 does not report a 2a yield; the blank cell is retained as unreported, not zero.")
    if y8 is not None:add_yield(r,5,f"{DOC_ID}:mol:8",y8,oid="yield_8")
    if i==1:
        o=sr.new_outcome("no_reaction","no_reaction",[ref],ev(5));o["value"]=f("starting compound 1aa remained unchanged","starting compound remained unchanged",ev(5));o["scope"]=f("TLC observation","TLC observation",ev(5));r["outcomes"].append(o)
    if i>=10:
        o=sr.new_outcome("recovered_1aa","recovery",[ref],ev(5));o["value"]=f("10–34% recovered across entries 10–13",{"minimum":10,"maximum":34,"unit":"%","raw_text":"10–34%"},ev(5));o["measurement_basis"]=f("substrate recovered","recovered starting material",ev(5));o["scope"]=f("Shared footnote range for entries 10–13; entry-specific value not assigned","group range; entry-specific attribution unresolved",ev(5));r["outcomes"].append(o);r["review"]["unresolved"].append("Footnote gives 10–34% recovered 1aa collectively for entries 10–13; no entry-specific recovery can be assigned.")
    procedure(r,5,"1aa (1.0 equiv), base (2.0 equiv) or oxidant (3.0 equiv; entry 9 uses 2.0 equiv), reported solvent/temperature/time; argon except entries 1–2")
    classify(r,5,"aromatize" if y2 and y2>0 else None)


# Table 4, page 5.
for prod,hours,yld in [("2a",6,80),("2b",5,79),("2c",5,80),("2d",5,85),
                       ("2e",8,74),("2f",8,63),("2g",7,87),("2h",6,93)]:
    precursor={"2a":"1aa","2b":"1ag","2c":"1af","2d":"1ab","2e":"1ca","2f":"1ch","2g":"1db","2h":"1ga"}[prod]
    r=reaction(f"table4_{prod}",5,categories=["substrate_scope"],name=f"Table 4 product {prod}",transformation="DDQ-mediated aromatization")
    if f"{DOC_ID}:mol:{precursor}" not in molecules:molecule(precursor,precursor,4,precursor,chirality="unspecified")
    ref=f"{DOC_ID}:mol:{precursor}"
    participant(r,ref,"reactant",f"p_{precursor}",ev(5),("equivalents","1.0 equiv",{"number":1.0,"unit":"equiv"}),ref)
    participant(r,named("DDQ",5,"DDQ"),"reagent","p_DDQ",ev(5),("equivalents","3.0 equiv",{"number":3.0,"unit":"equiv"}),ref)
    participant(r,named("toluene",5,"toluene"),"solvent","p_toluene",ev(5));participant(r,f"{DOC_ID}:mol:{prod}","product",f"p_{prod}",ev(5))
    add_step(r,5,100,(hours,"h"));add_yield(r,5,f"{DOC_ID}:mol:{prod}",yld);procedure(r,5,"THBC 1 (1.0 equiv), DDQ (3.0 equiv), toluene, 100 °C, reported time");classify(r,5,"aromatize")


# Page 6: Scheme 3, Scheme 4, Scheme 5.
for prod,iod,yld,dr,hours in [("9a","7a",63,"1.5:1",5),("9b","7b",65,"9:1",2),("9i","7i",61,"9:1",1),("9g","7g",56,"2.2:1",5)]:
    r=reaction(f"scheme3_{prod}",6,categories=["substrate_scope"],name=f"Scheme 3 product {prod}",transformation="Pd-catalyzed aza-spiroindolenine formation")
    ref=f"{DOC_ID}:mol:6h"
    participant(r,ref,"reactant","p_6h",ev(6),("equivalents","1.0 equiv",{"number":1.0,"unit":"equiv"}),ref)
    participant(r,f"{DOC_ID}:mol:{iod}","reactant",f"p_{iod}",ev(6),("equivalents","1.1 equiv",{"number":1.1,"unit":"equiv"}),ref)
    for key,role,n in [("Pd(PPh3)4","catalyst",0.05),("PCy3","ligand",0.12),("Cs2CO3","reagent",4.0)]:participant(r,named(key,6,key),role,f"p_{key}",ev(6),("equivalents",f"{n} equiv",{"number":n,"unit":"equiv"}),ref)
    participant(r,named("DMF",6,"dry DMF"),"solvent","p_DMF",ev(6));participant(r,f"{DOC_ID}:mol:{prod}","product",f"p_{prod}",ev(6))
    add_step(r,6,100,(hours,"h"),"argon");add_yield(r,6,f"{DOC_ID}:mol:{prod}",yld,isolated=False,basis_missing=True);add_selectivity(r,6,f"{DOC_ID}:mol:{prod}",dr,"1H NMR spectroscopy")
    if prod in {"9a","9g"}:r["notes"].append(f("inseparable mixture of diastereomers","inseparable diastereomer mixture",ev(6)))
    procedure(r,6,"6h (1.0 equiv), iodide 7 (1.1 equiv), Pd(PPh3)4 (0.05 equiv), PCy3 (0.12 equiv), Cs2CO3 (4.0 equiv), dry DMF, argon, 100 °C, reported time")
    classify(r,6,"cc")
for prod,sub,yld in [("1ha","9a",63),("1hb","9b",69)]:
    r=reaction(f"scheme4_{prod}",6,categories=["derivatization"],name=f"Scheme 4 {sub} to {prod}",transformation="acid-mediated spiroindolenine-to-THBC rearrangement")
    ref=f"{DOC_ID}:mol:{sub}";participant(r,ref,"reactant",f"p_{sub}",ev(6),("equivalents","1.0 equiv",{"number":1.0,"unit":"equiv"}),ref)
    participant(r,named("pTsOH",6,"p-TsOH·H2O"),"catalyst","p_pTsOH",ev(6),("mol_percent","10 mol%",{"number":10,"unit":"mol%"}),ref);participant(r,named("THF",6,"dry THF"),"solvent","p_THF",ev(6));participant(r,f"{DOC_ID}:mol:{prod}","product",f"p_{prod}",ev(6))
    add_step(r,6,"rt","overnight","argon");add_yield(r,6,f"{DOC_ID}:mol:{prod}",yld);procedure(r,6,"9 (1.0 equiv), p-TsOH·H2O (10 mol%), dry THF, argon, rt, overnight");classify(r,6,"cc")
for sub,prod,yld in [("2e","10a",45),("2f","10b",40)]:
    r=reaction(f"scheme5_{sub}_to_{prod}",6,categories=["derivatization"],name=f"Scheme 5 {sub} to {prod}",transformation="oxidative cleavage to aroyl β-carboline")
    ref=f"{DOC_ID}:mol:{sub}";participant(r,ref,"reactant",f"p_{sub}",ev(6),("equivalents","1.0 equiv",{"number":1.0,"unit":"equiv"}),ref)
    for key,role,field,raw,norm in [("RuCl3","catalyst","mol_percent","5 mol%",{"number":5,"unit":"mol%"}),("NaIO4","reagent","equivalents","6 equiv",{"number":6,"unit":"equiv"})]:participant(r,named(key,6,key),role,f"p_{key}",ev(6),(field,raw,norm),ref)
    for key in ["MeCN","CCl4","water"]:participant(r,named(key,6,key),"solvent",f"p_{key}",ev(6))
    participant(r,f"{DOC_ID}:mol:{prod}","product",f"p_{prod}",ev(6));r["solvent_conditions"]["composition"]=f("CH3CN/CCl4/H2O",["CH3CN","CCl4","H2O"],ev(6));r["solvent_conditions"]["mixture_ratio"]=f("2:1:1",[2,1,1],ev(6));add_step(r,6,"rt",(1,"h"));add_yield(r,6,f"{DOC_ID}:mol:{prod}",yld);procedure(r,6,"2e–f (1.0 equiv), RuCl3 (5 mol%), NaIO4 (6 equiv), CH3CN/CCl4/H2O 2:1:1, rt, 1 h");classify(r,6,"carbonyl")
for sub,prod,yld in [("10a","11a",85),("10b","11b",80)]:
    r=reaction(f"scheme5_{sub}_to_{prod}",6,categories=["derivatization"],name=f"Scheme 5 {sub} to {prod}",transformation="AlCl3-mediated debenzylation")
    ref=f"{DOC_ID}:mol:{sub}";participant(r,ref,"reactant",f"p_{sub}",ev(6),("equivalents","1.0 equiv",{"number":1.0,"unit":"equiv"}),ref);participant(r,named("AlCl3",6,"anhydrous AlCl3"),"reagent","p_AlCl3",ev(6),("equivalents","3.0 equiv",{"number":3.0,"unit":"equiv"}),ref);participant(r,named("toluene",6,"dry toluene"),"solvent","p_toluene",ev(6));participant(r,f"{DOC_ID}:mol:{prod}","product",f"p_{prod}",ev(6));add_step(r,6,"rt",(5,"h"));add_yield(r,6,f"{DOC_ID}:mol:{prod}",yld);procedure(r,6,"10 (1.0 equiv), anhydrous AlCl3 (3.0 equiv), dry toluene, rt, 5 h");classify(r,6,"nh")
r=reaction("scheme5_ref20a_formal_total_synthesis",6,"literature_example",name="Reference 20a conversion from 11b to eudistomins Y1/Y2")
participant(r,f"{DOC_ID}:mol:11b","reactant","p_11b",ev(6));procedure(r,6,"ref. 20a (from 11b); conditions and yields not given on inspected page");classify(r,6,None)


bundle["molecules"] = list(molecules.values())
bundle["source_context"] = [
    {"context_id":"page1_figure1","document_id":DOC_ID,"evidence_ids":ev(1),"type":"nonreaction_structure_context","text":"Figure 1 shows biologically active THBC/β-carboline examples; no synthetic arrow or experiment is reported."},
    {"context_id":"page3_results_not_shown","document_id":DOC_ID,"evidence_ids":ev(3),"type":"unresolved_discussion","text":"Different organic/inorganic bases are said to have been tried, but identities/results are not shown; no experiments fabricated."},
    {"context_id":"page4_unreactive_statement","document_id":DOC_ID,"evidence_ids":ev(4),"type":"unresolved_control","text":"Aryl bromide and an indole-N-sulfonyl allenamide are described as unreactive, but exact substrate identities/conditions are insufficient for separate records."},
    {"context_id":"page6_noesy","document_id":DOC_ID,"evidence_ids":ev(6),"type":"stereochemical_evidence","text":"Figure 2 NOESY correlations for 9b support the authors' relative 3R/11S assignment; cited SI Figure S10/Table S5 was not provided."},
    {"context_id":"page6_mechanism_lead_in","document_id":DOC_ID,"evidence_ids":ev(6),"type":"mechanism_context","text":"Mechanistic prose begins on page 6, but Scheme 6 and continuation are outside the declared page 1–6 coverage."},
]

coverage_records = {p:[] for p in range(1,7)}
for r in bundle["reactions"]:
    for p in range(1,7):
        if any(eid.startswith(f"{DOC_ID}:page{p}:") for eid in r["evidence_ids"]):
            coverage_records[p].append(r["reaction_id"])
coverage = {"document_id":DOC_ID,"semantic_accuracy_verified":False,"pages":[
    {"page":1,"state":"attempted","regions":[
        {"region":"title, abstract, and graphical abstract","finding":"Title metadata recovered; graphical abstract is a duplicate generic summary of the authors' route, not a specific experiment.","evidence_ids":ev(1)},
        {"region":"Figure 1 and introduction","finding":"Biologically active example structures and background only; no reaction record created.","evidence_ids":ev(1)}],"record_ids":coverage_records[1],"unresolved":["Small structures in Figure 1 were not registered because they are contextual examples, not experimental participants on pages 1–6."]},
    {"page":2,"state":"attempted","regions":[
        {"region":"Scheme 1a–e","finding":"Five prior-work arrows classified as literature examples; not counted as experiments of this paper.","evidence_ids":ev(2)},
        {"region":"Scheme 1f and surrounding text","finding":"Present-work route captured as a generic scheme, without assigning scope yields.","evidence_ids":ev(2)}],"record_ids":coverage_records[2],"unresolved":["Graphical structures are generic and conditions for cited work are intentionally not reconstructed from external references."]},
    {"page":3,"state":"attempted","regions":[
        {"region":"Table 1 rows 1–14 and footnote","finding":"Fourteen independent optimization experiments captured with catalysts, ligands, bases, solvents, temperature, time, and isolated yields.","evidence_ids":ev(3)},
        {"region":"results/discussion around Table 1","finding":"Qualitative discussion cross-checked; unspecified 'results not shown' trials retained as source context only.","evidence_ids":ev(3)}],"record_ids":coverage_records[3],"unresolved":["Structures of complex 6a/1aa remain machine-unverified; unspecified trials cannot be enumerated."]},
    {"page":4,"state":"attempted","regions":[
        {"region":"Table 2 substrate scope and footnotes","finding":"Twenty-two independent scope experiments captured; four dr >99:1 results target their products and cis C1/C4 relationship retained.","evidence_ids":ev(4)},
        {"region":"Scheme 2","finding":"Acid-catalyzed 6a to 1a experiment captured with 43% isolated yield.","evidence_ids":ev(4)},
        {"region":"opening Table 3 discussion","finding":"Narrative leading to Table 3 checked; the table itself is recorded from page 5.","evidence_ids":ev(4)}],"record_ids":coverage_records[4],"unresolved":["Absolute configurations and dr measurement method for Table 2 are unreported on pages 1–6; cited SI was not provided.","Exact identities for the qualitative aryl-bromide and indole-N-sulfonyl no-reaction statements are insufficient."]},
    {"page":5,"state":"attempted","regions":[
        {"region":"Table 3 rows 1–13 and footnotes","finding":"Thirteen independent optimization experiments captured, including zero yields, entry-4 blank 2a cell, atmosphere exceptions, and the shared 10–34% recovery range.","evidence_ids":ev(5)},
        {"region":"Table 4 eight entries","finding":"Eight DDQ scope experiments captured with time and isolated yields.","evidence_ids":ev(5)},
        {"region":"spiroindolenine discussion","finding":"Narrative checked; specific Scheme 3 arrows are recorded from page 6.","evidence_ids":ev(5)}],"record_ids":coverage_records[5],"unresolved":["The 10–34% recovered 1aa range is not attributable to individual entries 10–13.","Complex product structures remain machine-unverified."]},
    {"page":6,"state":"attempted","regions":[
        {"region":"Scheme 3 four arrows","finding":"Four aza-spiroindolenine scope experiments captured with times, yields, and product-targeted dr values measured by 1H NMR.","evidence_ids":ev(6)},
        {"region":"Figure 2 NOESY","finding":"Stereochemical evidence retained as source context and linked to molecule 9b; it is not a synthetic experiment.","evidence_ids":ev(6)},
        {"region":"Scheme 4 two arrows","finding":"Two acid-mediated conversions to THBCs captured.","evidence_ids":ev(6)},
        {"region":"Scheme 5 arrows","finding":"Four experiments of this paper captured; the arrow from 11b to eudistomins is separately classified as a literature example.","evidence_ids":ev(6)},
        {"region":"mechanism lead-in text","finding":"Mechanistic discussion begins, but the mechanism scheme/continuation lies beyond page 6 and is not claimed covered.","evidence_ids":ev(6)}],"record_ids":coverage_records[6],"unresolved":["Scheme 3 yield measurement basis is not stated on page 6.","NOESY details cited in SI are unavailable.","Mechanistic Scheme 6 is outside this batch's declared page coverage."]},
]}

# Literal cross-region hints are saved separately and never used to merge experiments.
mentions = sr.source_mentions(dataset, ["1aa", "2a", "6a", "9b", "11b"], 140)
dump(RUN / "raw_model/stereo_diagnostics_pages_0001_0006.json", diagnostics)
dump(RUN / "raw_model/pages_0001_0006_mentions.json", mentions)
dump(RUN / "raw_model/pages_0001_0006_extraction.json", bundle)
dump(SUBMITTED, bundle)
dump(RUN / "coverage/pages_0001_0006.json", coverage)
print(json.dumps({"molecules":len(bundle["molecules"]),"reactions":len(bundle["reactions"]),
                  "specific_experiments":sum(r["record_type"]=="specific_experiment" for r in bundle["reactions"]),
                  "page_counts":{str(p):len(v) for p,v in coverage_records.items()}},ensure_ascii=False,indent=2))
