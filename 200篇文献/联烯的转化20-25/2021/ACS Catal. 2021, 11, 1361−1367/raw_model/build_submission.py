#!/usr/bin/env python3
"""Author the page-1--6 v3.1 extraction from manually inspected source images."""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DOC = "doc_448eb03f81b6e3bb"


def write(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def ev(page: int):
    return [f"{DOC}:page{page}:image", f"{DOC}:page{page}:text"]


def fact(raw=None, normalized=None, evidence=(), method="original_extraction", derivation=None,
         status=None):
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation, status=status)


def derived(normalized, evidence, derivation, raw=None, method="inference"):
    return fact(raw, normalized, evidence, method=method, derivation=derivation)


def qty(number, unit, raw=None, qualifier=None):
    value = {"number": number, "unit": unit}
    if qualifier:
        value["qualifier"] = qualifier
    if raw:
        value["raw_text"] = raw
    return value


def molecule(mid, label, page, *, name=None, smiles=None, chirality="unresolved",
             geometry="unresolved", representation=None, assignment=None, mixture=None,
             note=None):
    evidence = ev(page)
    m = st.new_molecule(f"{DOC}:mol:{mid}", DOC, evidence)
    m["occurrences"][0]["compound_label"] = fact(label, label, evidence)
    m["names"]["original"] = fact(name or label, name or label, evidence)
    m["names"]["english"] = fact(name, name, evidence) if name else fact(evidence=evidence)
    m["names"]["abbreviations"] = fact(evidence=evidence)
    m["names"]["aliases"] = fact(evidence=evidence)
    m["cas_numbers"] = fact(evidence=evidence)
    if smiles:
        m["identity_status"] = "source_identified"
        m["representation_type"] = "smiles"
        m["smiles"] = derived(
            smiles, evidence,
            "Candidate connectivity was transcribed from the displayed structure and checked for RDKit parsing; parsing is not source verification.",
            raw=f"structure drawing for {label}", method="structure_derivation")
        m["stereochemistry"]["isomeric_smiles"] = derived(
            smiles, evidence,
            "No additional stereochemical token was added beyond the candidate source drawing.",
            raw=f"candidate isomeric SMILES for {label}", method="structure_derivation")
    else:
        m["identity_status"] = "source_identified"
        m["representation_type"] = "unresolved"
        m["identity_note"] = fact(
            note or "Source label/drawing identified; full connectivity was not digitized from the main-text figure.",
            note or "Source label/drawing identified; full connectivity was not digitized from the main-text figure.", evidence)
    rep = representation or ("not_applicable" if smiles and chirality == "achiral" else "not_encoded" if chirality in {"specified", "partially_specified", "racemic", "mixture"} else "unresolved")
    for key, value in (("chirality_status", chirality), ("geometry_status", geometry),
                       ("representation_status", rep)):
        m["stereochemistry"][key] = derived(
            value, evidence,
            f"Status assigned conservatively from the source drawing for {label}; absence of a tag was not treated as proof of achirality or racemic composition.",
            raw=f"source drawing/status for {label}")
    if assignment is not None:
        m["stereochemistry"]["assignment"] = fact(assignment, assignment, evidence)
    if mixture is not None:
        m["stereochemistry"]["mixture"] = fact(mixture, mixture, evidence)
    return m


def participant(pid, mid, role, evidence, *, amount=None, volume=None, mol_percent=None,
                reference=None, product_role=None, stage="reaction"):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    if product_role:
        p["product_role"] = fact(product_role, product_role, evidence)
    if amount:
        p["charges"]["amount"] = fact(amount[0], amount[1], evidence)
    if volume:
        p["charges"]["volume"] = fact(volume[0], volume[1], evidence)
    if mol_percent:
        p["charges"]["mol_percent"] = fact(mol_percent[0], mol_percent[1], evidence)
    if reference:
        p["equivalence_reference"] = fact(reference, reference, evidence)
    return p


def outcome(oid, kind, targets, evidence, *, raw=None, value=None, product_role=None,
            yield_type=None, basis=None, scope=None):
    o = st.new_outcome(oid, kind, targets, evidence)
    o["value"] = fact(raw, value, evidence)
    if product_role:
        o["product_role"] = fact(product_role, product_role, evidence)
    if yield_type:
        o["yield_type"] = fact(yield_type, yield_type, evidence)
    if basis:
        o["measurement_basis"] = fact(basis, basis, evidence)
    if scope is not None:
        o["scope"] = fact(scope, scope, evidence)
    return o


def set_classification(r, labels=None, changes=None, coverage="unresolved", limitation=None):
    evidence = r["evidence_ids"]
    r["classification"]["coverage"] = derived(
        coverage, evidence,
        "Net bond classification was inferred from the displayed transformation; no complete atom mapping was asserted.")
    if changes:
        r["classification"]["bond_changes"] = [
            derived(change, evidence, "Partial element-level bond change inferred from the source reaction drawing.", method="structure_derivation")
            for change in changes
        ]
    if labels is not None:
        r["classification"]["labels"] = derived(
            labels, evidence,
            "Taxonomy labels include all and only the new/promoted bonds represented in this partial change list.",
            method="structure_derivation")
    r["classification"]["limitations"] = [limitation or "Full atom mapping and complete coproduct/hydrogen balance are unresolved from the supplied main-text pages."]


def cc2():
    changes = [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
    ]
    labels = [{"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1", "b2"]}]
    return labels, changes


def one_label(code, number, atoms, *, count=1, change_type="formed", before=0, after=1):
    changes = []
    for i in range(count):
        changes.append({"change_id": f"b{i+1}", "atoms": [{"atom_map": None, "element": atoms[0]}, {"atom_map": None, "element": atoms[1]}],
                        "before_order": before, "after_order": after, "change_type": change_type})
    return [{"category_number": number, "code": code, "supporting_change_ids": [c["change_id"] for c in changes]}], changes


def reaction(rid, pages, *, record_type="specific_experiment", categories=None,
             location=None, name=None, transformation=None, unresolved=None):
    evidence = []
    for page in pages:
        evidence.extend(ev(page))
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    if categories is not None:
        r["experiment_categories"] = fact(categories, categories, evidence)
    r["location"] = fact(location or rid, location or rid, evidence)
    if name:
        r["reaction_name"] = fact(name, name, evidence)
    if transformation:
        r["transformation_type"] = derived(
            transformation, evidence, "Transformation type normalized from the article's reaction name and displayed arrow.")
    r["review"]["unresolved"] = unresolved or ["Full structures, atom mapping, and SI-only operational details require expert/source review."]
    return r


def add_step(r, sid, participant_ids, *, operation, temperature=None, time=None):
    evidence = r["evidence_ids"]
    s = st.new_step(sid, evidence)
    s["participant_ids"] = participant_ids
    s["operation"] = fact(operation, operation, evidence)
    if temperature is not None:
        s["temperature"] = fact(temperature[0], temperature[1], evidence)
    if time is not None:
        s["time"] = fact(time[0], time[1], evidence)
    r["steps"].append(s)


def add_solvent(r, name, volume=None):
    r["solvent_conditions"]["composition"] = fact(name, [name], r["evidence_ids"])
    if volume:
        r["solvent_conditions"]["total_volume"] = fact(volume[0], volume[1], r["evidence_ids"])


def add_scale(r, mid, raw, value):
    r["reaction_scale"] = fact(raw, {"reference_molecule_id": mid, "amount": value}, r["evidence_ids"])


dataset = json.loads((RUN / "dataset.json").read_text())
dataset["documents"][0]["title"] = "Dielectrophilic Allenic Ketone-Enabled [4 + 2] Annulation with 3,3′-Bisoxindoles: Enantioselective Creation of Two Contiguous Quaternary Stereogenic Centers"
dataset["documents"][0]["supplementary_material"] = {
    "status": "mentioned_not_provided",
    "note": "The article states that Supporting Information exists, but no SI file was supplied or explicitly associated in this session."
}
write(RUN / "dataset.json", dataset)
bundle = st.empty(dataset["documents"])

# Molecules with candidate connectivity; only 1a/1b are primary source-drawing candidates used to reduce identity ambiguity.
molecules = {}
def add_mol(key, *args, **kwargs):
    m = molecule(key, *args, **kwargs)
    molecules[key] = m
    bundle["molecules"].append(m)
    return m["molecule_id"]


M = {}
M["1a"] = add_mol("1a", "1a", 2, name="β′-acetate allenic ketone 1a", smiles="C=C=C(COC(C)=O)C(=O)c1ccccc1", chirality="achiral", geometry="not_applicable")
M["1b"] = add_mol("1b", "1b", 2, name="β′-acetate allenoate 1b", smiles="C=C=C(COC(C)=O)C(=O)OCc1ccccc1", chirality="achiral", geometry="not_applicable")
simple = {
    "PPh3": ("PPh3", "triphenylphosphine", "P(c1ccccc1)(c1ccccc1)c1ccccc1", 2),
    "toluene": ("toluene", "toluene", "Cc1ccccc1", 2),
    "THF": ("THF", "tetrahydrofuran", "C1CCOC1", 3),
    "CH2Cl2": ("CH2Cl2", "dichloromethane", "ClCCl", 3),
    "CH3CN": ("CH3CN", "acetonitrile", "CC#N", 3),
    "dioxane": ("1,4-dioxane", "1,4-dioxane", "O1CCOCC1", 3),
    "MeOH": ("MeOH", "methanol", "CO", 4),
    "MeI": ("CH3I", "iodomethane", "CI", 4),
    "MeNH2": ("MeNH2", "methylamine", "CN", 4),
    "O3": ("O3", "ozone", "[O-][O+]=O", 4),
}
for key, (label, name, smiles, page) in simple.items():
    M[key] = add_mol(key, label, page, name=name, smiles=smiles, chirality="achiral", geometry="not_applicable")

for key, label, name, page in [
    ("Cs2CO3", "Cs2CO3", "cesium carbonate", 2), ("Na2CO3", "Na2CO3", "sodium carbonate", 3),
    ("K2CO3", "K2CO3", "potassium carbonate", 3), ("Li2CO3", "Li2CO3", "lithium carbonate", 3),
    ("PbOAc4", "Pb(OAc)4", "lead tetraacetate", 4), ("Jones", "Jones reagent", "Jones reagent", 4),
    ("TFA", "TFA", "trifluoroacetic acid", 4), ("LDA", "LDA", "lithium diisopropylamide", 4),
    ("DIBAL", "DIBAL-H", "diisobutylaluminium hydride", 4), ("LiAlH4", "LiAlH4", "lithium aluminium hydride", 4),
    ("RedAl", "Red-Al", "sodium bis(2-methoxyethoxy)aluminium hydride", 4),
]:
    M[key] = add_mol(key, label, page, name=name, chirality="unresolved", geometry="not_applicable", representation="unresolved")

# Source-labelled substrates/products.
for suffix in "abcdef":
    M[f"2{suffix}"] = add_mol(f"2{suffix}", f"2 (row giving 3{suffix})", 2, chirality="unresolved", geometry="unresolved")
    M[f"3{suffix}"] = add_mol(f"3{suffix}", f"3{suffix}", 2, chirality="unresolved", geometry="unresolved")
    M[f"3{suffix}p"] = add_mol(f"3{suffix}p", f"3{suffix}′", 2, chirality="unresolved", geometry="unresolved")

M["4a"] = add_mol("4a", "4a", 3, name="3,3′-bisoxindole 4a", chirality="partially_specified", geometry="not_applicable", assignment="Relative stereochemical drawing shown; detailed assignment not stated on pages 1–6.")
M["4aD"] = add_mol("4aD", "4a-D", 4, name="deuterated 3,3′-bisoxindole 4a-D", chirality="partially_specified", geometry="not_applicable", assignment="Deuterated substrate used for KIE; isotope position requires Eq. 1/SI review.")
M["5arac"] = add_mol("5arac", "5a-rac", 3, chirality="racemic", geometry="not_applicable", mixture="Racemic product as explicitly labelled 5a-rac.")
M["5aprac"] = add_mol("5aprac", "5a′-rac", 3, chirality="racemic", geometry="not_applicable", mixture="Racemic product as explicitly labelled 5a′-rac.")
M["5a"] = add_mol("5a", "5a", 3, chirality="specified", geometry="not_applicable", assignment="Displayed stereoisomer; configurations of scope products other than 5g were assigned by analogy to X-ray-defined 5g.")
for i in range(1, 7):
    M[f"P{i}"] = add_mol(f"P{i}", f"P{i}", 3, name=f"bifunctional phosphine catalyst P{i}", chirality="specified", geometry="not_applicable", assignment="Stereochemical bonds are shown in the catalyst drawing; no isomeric SMILES encoded.")
for suffix in "bcdefghijklm":
    M[f"4{suffix}"] = add_mol(
        f"4{suffix}", f"4{suffix}", 4, chirality="partially_specified",
        geometry="not_applicable",
        assignment="Relative stereochemical drawing is shown for the bisoxindole substrate; absolute descriptors are not stated on pages 1–6.")
for suffix in "nopq":
    M[f"1{suffix}"] = add_mol(f"1{suffix}", f"1{suffix}", 4, name=f"allenic ketone 1{suffix}", chirality="unresolved", geometry="unresolved")
for suffix in "bcdefghijklmnopq":
    assignment = "Absolute configuration determined by X-ray crystallography (CCDC 1970725)." if suffix == "g" else "Displayed stereoisomer; assigned by analogy to X-ray-defined 5g."
    M[f"5{suffix}"] = add_mol(f"5{suffix}", f"5{suffix}", 4, chirality="specified", geometry="not_applicable", assignment=assignment)
for key in ["6", "7", "8", "9", "10", "11"]:
    M[key] = add_mol(key, key, 4, chirality="specified", geometry="not_applicable", assignment="Stereochemical drawing shown in Scheme 1; configuration retained through the sequence unless otherwise stated.")
M["folicanthine"] = add_mol("folicanthine", "(−)-folicanthine", 4, name="(−)-folicanthine", chirality="specified", geometry="not_applicable", assignment="Product identity and negative optical sign are explicit; full isomeric SMILES not encoded.")

# Generic/literature records on pages 1--2.
for letter, label in zip("abcd", ["Lu [3+2] annulation", "Kwon [4+2] annulation", "allenic-ketone [4+2] annulation", "Tong [4+1] annulation"]):
    r = reaction(f"fig1_literature_{letter}", [1], record_type="literature_example", location=f"Figure 1{letter}", name=label,
                 transformation="phosphine-mediated allene annulation", unresolved=["Cited prior-work generic structures were not expanded into current-paper experiments."])
    set_classification(r, coverage="unresolved", limitation="Literature-example generic scheme; no current-paper experiment or complete atom mapping asserted.")
    bundle["reactions"].append(r)

r = reaction("fig2_hypothesis", [2], record_type="generic_scheme", location="Figure 2", name="Design of a novel dielectrophilic C4 synthon",
             transformation="hypothesized phosphine activation and [4+2] annulation", unresolved=["Figure 2 is a hypothesis, not a concrete experiment."])
set_classification(r, coverage="unresolved", limitation="Hypothetical generic transformation; substituent variables and atom mapping unresolved.")
bundle["reactions"].append(r)

r = reaction("fig3_general_41", [2], record_type="generic_scheme", location="Figure 3 header", name="PPh3-catalyzed [4+1] annulation",
             transformation="[4+1] annulation", unresolved=["Generic row; concrete products are represented by separate Figure 3 records."])
labels, changes = cc2(); set_classification(r, labels, changes, "partial")
bundle["reactions"].append(r)

# Figure 3: twelve current-paper comparison experiments.
fig3_yields = {"a": (85, 57), "b": (87, 68), "c": (87, 64), "d": (99, 78), "e": (98, 78), "f": (96, 69)}
for suffix, (yield_a, yield_b) in fig3_yields.items():
    for allene, product_key, y in [("1a", f"3{suffix}", yield_a), ("1b", f"3{suffix}p", yield_b)]:
        rid = f"fig3_{product_key}"
        r = reaction(rid, [2], categories=["other"], location=f"Figure 3, product {product_key.replace('p','′')}", name="PPh3-catalyzed [4+1] annulation", transformation="[4+1] annulation")
        ids = []
        for pid, key, role, amount in [
            ("allene", allene, "reactant", ("0.12 mmol", qty(0.12, "mmol"))),
            ("partner", f"2{suffix}", "reactant", ("0.10 mmol", qty(0.10, "mmol"))),
            ("base", "Cs2CO3", "reagent", ("0.12 mmol", qty(0.12, "mmol"))),
            ("cat", "PPh3", "catalyst", ("0.02 mmol", qty(0.02, "mmol"))),
        ]:
            r["participants"].append(participant(pid, M[key], role, r["evidence_ids"], amount=amount,
                                                 mol_percent=("20 mol%", qty(20, "mol%")) if pid == "cat" else None,
                                                 reference=M[f"2{suffix}"] if pid == "cat" else None)); ids.append(pid)
        r["participants"].append(participant("solvent", M["toluene"], "solvent", r["evidence_ids"], volume=("2 mL", qty(2, "mL")))); ids.append("solvent")
        r["participants"].append(participant("product", M[product_key], "product", r["evidence_ids"], product_role="main product")); ids.append("product")
        add_scale(r, M[f"2{suffix}"], "2 (0.10 mmol)", qty(0.10, "mmol")); add_solvent(r, "toluene", ("2 mL", qty(2, "mL")))
        add_step(r, "reaction", ids, operation="PPh3-catalyzed annulation", temperature=("room temperature", "r.t."), time=("5 h", qty(5, "h")))
        r["outcomes"].append(outcome("yield", "yield", [M[product_key]], r["evidence_ids"], raw=f"{y}%", value=qty(y, "%"), product_role="main product", yield_type="isolated yield", basis="isolated product"))
        r["outcomes"].append(outcome("conversion", "conversion", [M[allene]], r["evidence_ids"], raw="fully consumed", value="fully consumed", basis="source statement; analytical method unreported"))
        r["general_scheme_ids"] = [f"{DOC}:fig3_general_41"]
        labels, changes = cc2(); set_classification(r, labels, changes, "partial")
        bundle["reactions"].append(r)

# Figure 4 feasibility pair.
for allene, product_key, raw_y, y in [("1a", "5arac", ">99% yield", qty(99, "%", raw=">99%", qualifier="gt")), ("1b", "5aprac", "35% yield", qty(35, "%"))]:
    r = reaction(f"fig4_{allene}", [2, 3], categories=["other"], location=f"Figure 4, {allene}", name="Feasibility study of projected [4+2] annulation", transformation="[4+2] annulation")
    specs = [("allene", allene, "reactant", "0.10 mmol", 0.10), ("bisox", "4a", "reactant", "0.12 mmol", 0.12),
             ("base", "Na2CO3", "reagent", "0.12 mmol", 0.12), ("cat", "PPh3", "catalyst", "0.02 mmol", 0.02)]
    ids=[]
    for pid,key,role,raw,n in specs:
        r["participants"].append(participant(pid,M[key],role,r["evidence_ids"],amount=(raw,qty(n,"mmol")),
                                             mol_percent=("20 mol%",qty(20,"mol%")) if pid=="cat" else None,
                                             reference=M[allene] if pid=="cat" else None));ids.append(pid)
    r["participants"].append(participant("solvent",M["toluene"],"solvent",r["evidence_ids"],volume=("2 mL",qty(2,"mL"))));ids.append("solvent")
    r["participants"].append(participant("product",M[product_key],"product",r["evidence_ids"],product_role="main product"));ids.append("product")
    add_scale(r,M[allene],f"{allene} (0.10 mmol)",qty(0.10,"mmol"));add_solvent(r,"toluene",("2 mL",qty(2,"mL")))
    add_step(r,"reaction",ids,operation="PPh3-catalyzed [4+2] annulation",temperature=("room temperature","r.t."),time=("7 h",qty(7,"h")))
    r["outcomes"].append(outcome("yield","yield",[M[product_key]],r["evidence_ids"],raw=raw_y,value=y,product_role="main product",yield_type="isolated yield",basis="isolated product"))
    r["outcomes"].append(outcome("conversion","conversion",[M[allene]],r["evidence_ids"],raw="fully consumed",value="fully consumed",basis="source statement; analytical method unreported"))
    labels,changes=cc2();set_classification(r,labels,changes,"partial")
    bundle["reactions"].append(r)

# Shared generic asymmetric [4+2] scheme for Table 1/Figure 5.
r = reaction("generic_asymmetric_42", [3, 4], record_type="generic_scheme", location="Table 1 header and Figure 5 header", name="Enantioselective [4+2] annulation of allenic ketones with 3,3′-bisoxindoles", transformation="enantioselective [4+2] annulation", unresolved=["R-group scope represented by separate concrete rows."])
labels,changes=cc2();set_classification(r,labels,changes,"partial")
bundle["reactions"].append(r)

# Table 1 optimization rows.
table = [
    (1,"P1","toluene","Na2CO3",4,99,15),(2,"P2","toluene","Na2CO3",4,99,40),
    (3,"P3","toluene","Na2CO3",6,99,50),(4,"P4","toluene","Na2CO3",8,98,77),
    (5,"P5","toluene","Na2CO3",8,98,96),(6,"P6","toluene","Na2CO3",8,98,10),
    (7,"P5","toluene","K2CO3",3,77,96),(8,"P5","toluene","Li2CO3",20,85,96),
    (9,"P5","THF","Na2CO3",8,30,96),(10,"P5","CH2Cl2","Na2CO3",8,90,96),
    (11,"P5","CH3CN","Na2CO3",8,"trace",None),(12,"P5","dioxane","Na2CO3",8,"no reaction",None),
]
for entry,cat,solvent,base,h,y,ee in table:
    r = reaction(f"table1_entry_{entry:02d}",[3],categories=["optimization"],location=f"Table 1 entry {entry}",name="Reaction screening",transformation="enantioselective [4+2] annulation")
    ids=[]
    for pid,key,role,n in [("allene","1a","reactant",0.12),("bisox","4a","reactant",0.10),("base",base,"reagent",0.12),("cat",cat,"catalyst",0.02)]:
        r["participants"].append(participant(pid,M[key],role,r["evidence_ids"],amount=(f"{n:g} mmol",qty(n,"mmol")),
                                             mol_percent=("20 mol%",qty(20,"mol%")) if pid=="cat" else None,
                                             reference=M["4a"] if pid=="cat" else None));ids.append(pid)
    r["participants"].append(participant("solvent",M[solvent],"solvent",r["evidence_ids"],volume=("2 mL",qty(2,"mL"))));ids.append("solvent")
    r["participants"].append(participant("product",M["5a"],"product",r["evidence_ids"],product_role="main product"));ids.append("product")
    add_scale(r,M["4a"],"4a (0.10 mmol)",qty(0.10,"mmol"));add_solvent(r,"1,4-dioxane" if solvent=="dioxane" else solvent,("2 mL",qty(2,"mL")))
    add_step(r,"reaction",ids,operation="Catalyst/base/solvent screening",temperature=("room temperature","r.t."),time=(f"{h} h",qty(h,"h")))
    if isinstance(y,int):
        r["outcomes"].append(outcome("yield","yield",[M["5a"]],r["evidence_ids"],raw=f"{y}%",value=qty(y,"%"),product_role="main product",yield_type="isolated yield",basis="isolated product"))
    elif y=="trace":
        r["outcomes"].append(outcome("trace","trace",[M["5a"]],r["evidence_ids"],raw="trace",value="trace",product_role="main product",basis="table yield column; quantitative method not applicable"))
    else:
        r["outcomes"].append(outcome("no_reaction","no_reaction",[],r["evidence_ids"],raw="No reaction",value=True,basis="Table 1 footnote d"))
    if ee is not None:
        r["outcomes"].append(outcome("ee","selectivity",[M["5a"]],r["evidence_ids"],raw=f"{ee}% ee",value=qty(ee,"%"),product_role="main product",basis="HPLC analysis on a chiral stationary phase",scope={"metric":"ee","major_configuration":"not stated in Table 1"}))
    r["general_scheme_ids"]=[f"{DOC}:generic_asymmetric_42"]
    labels,changes=cc2();set_classification(r,labels,changes,"partial")
    bundle["reactions"].append(r)

# Figure 5 scope. The caption gives no general reaction time; only 5j/5k explicitly use 12 h.
scope_rows = {
    "b":(94,87,None),"c":(96,98,None),"d":(95,96,None),"e":(93,94,None),"f":(94,91,None),
    "g":(91,93,"4:1"),"h":(95,91,"4:1"),"i":(92,87,"4:1"),"j":(94,95,"1.8:1"),
    "k":(84,96,"5.4:1"),"l":(88,95,">20:1"),"m":(81,97,"3.5:1"),
    "n":(95,90,None),"o":(95,93,None),"p":(97,95,None),"q":(98,91,None),
}
for suffix,(y,ee,rr) in scope_rows.items():
    allene = "1a" if suffix in "bcdefghijklm" else f"1{suffix}"
    bisox = f"4{suffix}" if suffix in "bcdefghijklm" else "4a"
    product_key=f"5{suffix}"
    r=reaction(f"fig5_{product_key}",[3,4],categories=["substrate_scope"],location=f"Figure 5, {product_key}",name="Scope of the [4+2] annulation",transformation="enantioselective [4+2] annulation")
    ids=[]
    for pid,key,role,n in [("allene",allene,"reactant",0.10),("bisox",bisox,"reactant",0.12),("base","Na2CO3","reagent",0.12),("cat","P5","catalyst",0.02)]:
        r["participants"].append(participant(pid,M[key],role,r["evidence_ids"],amount=(f"{n:g} mmol",qty(n,"mmol")),
                                             mol_percent=("20 mol%",qty(20,"mol%")) if pid=="cat" else None,
                                             reference=M[allene] if pid=="cat" else None));ids.append(pid)
    r["participants"].append(participant("solvent",M["toluene"],"solvent",r["evidence_ids"],volume=("2 mL",qty(2,"mL"))));ids.append("solvent")
    r["participants"].append(participant("product",M[product_key],"product",r["evidence_ids"],product_role="main product"));ids.append("product")
    add_scale(r,M[allene],f"{allene} (0.10 mmol)",qty(0.10,"mmol"));add_solvent(r,"toluene",("2 mL",qty(2,"mL")))
    add_step(r,"reaction",ids,operation="P5-catalyzed scope reaction",temperature=("room temperature","r.t."),
             time=("12 h",qty(12,"h")) if suffix in "jk" else ("8 h",qty(8,"h")))
    r["outcomes"].append(outcome("yield","yield",[M[product_key]],r["evidence_ids"],raw=f"{y}% yield",value=qty(y,"%"),product_role="main product",yield_type="isolated yield",basis="isolated product"))
    r["outcomes"].append(outcome("ee","selectivity",[M[product_key]],r["evidence_ids"],raw=f"{ee}% ee",value=qty(ee,"%"),product_role="main product",basis="HPLC analysis on a chiral stationary phase",scope={"metric":"ee","major_configuration":"displayed/assigned product; detailed descriptor not stated"}))
    if rr:
        value={"ratio":rr,"qualifier":"gt" if rr.startswith(">") else "exact","major_regioisomer":"displayed 5"+suffix}
        r["outcomes"].append(outcome("rr","selectivity",[M[product_key]],r["evidence_ids"],raw=f"rr = {rr}",value=value,product_role="main product",basis="1H NMR analysis",scope={"metric":"regioisomer ratio","order":"displayed major:minor"}))
    r["general_scheme_ids"]=[f"{DOC}:generic_asymmetric_42"]
    labels,changes=cc2();set_classification(r,labels,changes,"partial")
    bundle["reactions"].append(r)

# Real mechanistic measurements and proposed mechanism are kept distinct.
r=reaction("mechanistic_31P_NMR",[4],categories=["mechanistic_study"],location="Page 4, Mechanism Insight paragraph",name="31P NMR monitoring",transformation="in situ reaction monitoring",unresolved=["Assignments of the two observed phosphorus resonances are not established in the main text."])
for pid,key,role in [("allene","1a","reactant"),("bisox","4a","reactant"),("cat","P5","catalyst")]:r["participants"].append(participant(pid,M[key],role,r["evidence_ids"]))
add_step(r,"monitoring",["allene","bisox","cat"],operation="Monitor phosphorus species during the reaction by 31P NMR")
r["outcomes"].append(outcome("resonances","other",[],r["evidence_ids"],raw="two new resonances at 31.7 and 28.5 ppm",value={"new_resonances_ppm":[31.7,28.5],"P5_ppm":-23.9,"phosphine_oxide_ppm":33.8,"intermediate_I_ppm":22.7},basis="31P NMR"))
set_classification(r,coverage="not_applicable",limitation="Analytical monitoring record; no net isolated transformation classified.")
bundle["reactions"].append(r)

r=reaction("mechanistic_reaction_orders",[4],categories=["mechanistic_study"],location="Page 4, kinetic studies",name="Kinetic reaction-order study",transformation="kinetic measurement",unresolved=["Underlying rate data and detailed protocol are in unavailable SI."])
for pid,key,role in [("allene","1a","reactant"),("bisox","4a","reactant"),("cat","P5","catalyst")]:r["participants"].append(participant(pid,M[key],role,r["evidence_ids"]))
add_step(r,"kinetics",["allene","bisox","cat"],operation="Measure reaction orders")
r["outcomes"].append(outcome("orders","other",[],r["evidence_ids"],raw="first-order in P5, allenic ketone, and bisoxindole",value={"P5":1,"allenic_ketone":1,"bisoxindole":1},basis="kinetic studies; details in unavailable SI"))
set_classification(r,coverage="not_applicable",limitation="Kinetic measurement record; no separate net transformation classified.")
bundle["reactions"].append(r)

r=reaction("mechanistic_KIE",[4],categories=["mechanistic_study"],location="Page 4, eq 1 discussion",name="Deuterium kinetic isotope effect experiment",transformation="kinetic isotope effect measurement",unresolved=["Eq. 1 structures/rate data and deuterium position are not legible/available in the supplied page image; SI details are unavailable."])
for pid,key,role in [("allene","1a","reactant"),("bisox_H","4a","reactant"),("bisox_D","4aD","reactant"),("cat","P5","catalyst"),("base","Na2CO3","reagent"),("solvent","toluene","solvent")]:r["participants"].append(participant(pid,M[key],role,r["evidence_ids"]))
add_step(r,"KIE_comparison",[p["participant_id"] for p in r["participants"]],operation="Compare otherwise identical reactions of 4a and 4a-D",temperature=("otherwise identical optimized conditions; room temperature in surrounding study","r.t."))
r["outcomes"].append(outcome("KIE","other",[],r["evidence_ids"],raw="KIE = 4",value={"KIE":4,"comparison":"4a/4a-D"},basis="relative reaction-rate comparison"))
set_classification(r,coverage="not_applicable",limitation="Comparative kinetic measurement; no distinct product mapping classified.")
bundle["reactions"].append(r)

r=reaction("fig6_proposed_mechanism",[3,4],record_type="proposed_mechanism",location="Figure 6",name="Proposed mechanism",transformation="proposed phosphine-catalytic cycle",unresolved=["Author-proposed intermediates I–V are not experimental reaction records."])
set_classification(r,coverage="not_applicable",limitation="Proposed catalytic cycle, not a net concrete experiment.")
bundle["reactions"].append(r)

# Scheme 1 total synthesis transformations.
def scheme_reaction(rid,start,products,steps,yields,extra_participants,classification):
    r=reaction(rid,[4,5],categories=["derivatization"],location=f"Scheme 1, {start} to {'/'.join(products)}",name="Concise total synthesis of (−)-folicanthine",transformation="multistep derivatization")
    r["participants"].append(participant("start",M[start],"reactant",r["evidence_ids"])); all_ids=["start"]
    for idx,(key,role) in enumerate(extra_participants,1):
        pid=f"agent{idx}";r["participants"].append(participant(pid,M[key],role,r["evidence_ids"]));all_ids.append(pid)
    for idx,key in enumerate(products,1):
        pid=f"product{idx}";r["participants"].append(participant(pid,M[key],"product",r["evidence_ids"],product_role="main product" if idx==1 else "co-product"));all_ids.append(pid)
    for sid,op,temp,time in steps:add_step(r,sid,all_ids,operation=op,temperature=temp,time=time)
    for idx,(key,raw_y,y) in enumerate(yields,1):
        r["outcomes"].append(outcome(f"yield{idx}","yield",[M[key]],r["evidence_ids"],raw=raw_y,value=qty(y,"%"),product_role="main product" if idx==1 else "co-product",yield_type="reported yield",basis="yield basis not stated in Scheme 1"))
    if classification:
        set_classification(r,*classification,coverage="partial")
    else:set_classification(r,coverage="unresolved")
    bundle["reactions"].append(r)

scheme_reaction("scheme1_5a_to_6","5a",["6"],[
    ("ozonolysis","Ozonolysis in CH2Cl2",("−78 °C",qty(-78,"°C")),None),
    ("reductive_workup","Then PPh3 at room temperature",("r.t.","r.t."),("1 h",qty(1,"h"))),
],[("6","85% yield",85)],[("O3","reagent"),("CH2Cl2","solvent"),("PPh3","reagent")],one_label("DB-C-O",14,("C","O"),after=2))
scheme_reaction("scheme1_6_to_7","6",["7"],[("oxidation","Pb(OAc)4 in MeOH",("0 °C",qty(0,"°C")),("2 h",qty(2,"h")))],[("7","94% yield",94)],[("PbOAc4","reagent"),("MeOH","solvent")],one_label("C-O",3,("C","O")))

# 7 -> 8 is explicitly a 3-step overall yield.
labels_cn,changes_cn=one_label("C-N",5,("C","N"),count=2)
labels_co,changes_co=one_label("C-O",3,("C","O"),count=1)
combined_changes=[]
for i,c in enumerate(changes_cn+changes_co,1):c=copy.deepcopy(c);c["change_id"]=f"b{i}";combined_changes.append(c)
combined_labels=[{"category_number":5,"code":"C-N","supporting_change_ids":["b1","b2"]},{"category_number":3,"code":"C-O","supporting_change_ids":["b3"]}]
scheme_reaction("scheme1_7_to_8","7",["8"],[
    ("jones_oxidation","Jones oxidation",("0 °C",qty(0,"°C")),("6 h",qty(6,"h"))),
    ("boc_removal","CH2Cl2/TFA",("0 °C",qty(0,"°C")),("0.5 h",qty(0.5,"h"))),
    ("methylation","THF/Cs2CO3, CH3I",("r.t.","r.t."),("4 h",qty(4,"h"))),
],[("8","73% yield for 3 steps",73)],[("Jones","reagent"),("CH2Cl2","solvent"),("TFA","reagent"),("THF","solvent"),("Cs2CO3","reagent"),("MeI","reagent")],(combined_labels,combined_changes))
scheme_reaction("scheme1_8_to_9","8",["9"],[("aminolysis","MeNH2 in MeOH",("50 °C",qty(50,"°C")),("24 h",qty(24,"h")))],[("9","98% yield",98)],[("MeNH2","reagent"),("MeOH","solvent")],one_label("C-N",5,("C","N"),count=2))
scheme_reaction("scheme1_9_to_10_11","9",["10","11"],[
    ("deprotonation","LDA in THF",("−30 °C",qty(-30,"°C")),("3 h",qty(3,"h"))),
    ("reductive_amination","DIBAL-H",("r.t.","r.t."),("12 h",qty(12,"h"))),
],[("10","16% yield",16),("11","36% yield",36)],[("LDA","reagent"),("THF","solvent"),("DIBAL","reagent")],one_label("C-N",5,("C","N"),count=1))
scheme_reaction("scheme1_10_to_11","10",["11"],[("cyclization","LiAlH4 in THF",("r.t.","r.t."),("5 h",qty(5,"h")))],[("11","47% yield",47)],[("LiAlH4","reagent"),("THF","solvent")],one_label("C-N",5,("C","N"),count=1))
scheme_reaction("scheme1_11_to_folicanthine","11",["folicanthine"],[("lactam_reduction","Red-Al",("r.t.","r.t."),("5 h",qty(5,"h")))],[("folicanthine","75% yield",75)],[("RedAl","reagent")],one_label("C-H",1,("C","H"),change_type="hydrogen_gain"))

# Source context: article metadata, page-level mechanistic interpretation, and unavailable SI statement.
bundle["source_context"] = [
    {"context_id":f"{DOC}:article_metadata","document_id":DOC,"kind":"article_metadata","title":dataset["documents"][0]["title"],"doi":"10.1021/acscatal.0c05225","evidence_ids":ev(1)},
    {"context_id":f"{DOC}:graphical_abstract","document_id":DOC,"kind":"generic_current_work_summary","note":"Graphical abstract repeats the generic [4+2] annulation and downstream folicanthine synthesis; it is not a separate experiment.","evidence_ids":ev(1)},
    {"context_id":f"{DOC}:mechanistic_interpretation","document_id":DOC,"kind":"author_interpretation","note":"The authors infer fast phosphonium formation/addition and likely rate-determining proton transfer; the second intramolecular proton transfer is proposed as crucial.","evidence_ids":ev(4)},
    {"context_id":f"{DOC}:si_statement","document_id":DOC,"kind":"supplementary_information_statement","note":"The article advertises experimental procedures, optimization tables, characterization data, and CIF, but no SI file was supplied or associated in this session.","evidence_ids":ev(5)},
]

write(RUN / "submitted.json", bundle)

raw = {
    "document_id":DOC,"pages":[1,2,3,4,5,6],"semantic_accuracy_verified":False,
    "findings":{
        "figure3_yields":fig3_yields,
        "figure4":[{"allene":"1a","product":"5a-rac","yield":">99%"},{"allene":"1b","product":"5a′-rac","yield":"35%"}],
        "table1":[{"entry":e,"catalyst":c,"solvent":s,"base":b,"time_h":h,"yield":y,"ee":ee} for e,c,s,b,h,y,ee in table],
        "figure5":scope_rows,
        "scheme1":["5a→6 85%","6→7 94%","7→8 73% over 3 steps","8→9 98%","9→10 16% + 11 36%","10→11 47%","11→(−)-folicanthine 75%"],
    },
    "limitations":["No SI was supplied.","Most complex structures are retained by source label/drawing rather than guessed SMILES.","All records are machine-authored and require review."],
}
write(RUN / "raw_model/pages_0001_0006_extraction.json", raw)

by_id={r["reaction_id"]:r for r in bundle["reactions"]}
page_records={1:[],2:[],3:[],4:[],5:[],6:[]}
for rid,r in by_id.items():
    pages={int(e.split(":page")[1].split(":")[0]) for e in r["evidence_ids"] if ":page" in e}
    for p in pages:
        if p in page_records:page_records[p].append(rid)
coverage={"document_id":DOC,"semantic_accuracy_verified":False,"pages":[
    {"page":1,"state":"attempted","regions":[
        {"region":"title/abstract/graphical abstract","finding":"Article metadata and a generic current-work summary; no standalone concrete experiment.","evidence_ids":ev(1)},
        {"region":"Figure 1a–d","finding":"Four cited prior-work generic annulation schemes, classified as literature examples.","evidence_ids":ev(1)}],"record_ids":sorted(page_records[1]),"unresolved":["Generic Figure 1 variables and exact literature conditions were not expanded."]},
    {"page":2,"state":"attempted","regions":[
        {"region":"Figure 2","finding":"Author hypothesis for a dielectrophilic C4 synthon; generic scheme, not a concrete experiment.","evidence_ids":ev(2)},
        {"region":"Figure 3","finding":"Twelve concrete [4+1] comparison experiments (3a–f and 3a′–f′) with isolated yields and shared conditions.","evidence_ids":ev(2)},
        {"region":"right-column feasibility prose","finding":"Narrative duplicates the two Figure 4 feasibility experiments; evidence merged with page 3 records.","evidence_ids":ev(2)}],"record_ids":sorted(page_records[2]),"unresolved":["Six C1 partner structures are identified by row/drawing but not assigned unsupported labels or complete SMILES."]},
    {"page":3,"state":"attempted","regions":[
        {"region":"Figure 4","finding":"Two concrete feasibility experiments using 1a or 1b, producing racemic 5a/5a′.","evidence_ids":ev(3)},
        {"region":"Table 1 entries 1–12","finding":"Twelve distinct optimization experiments transcribed row by row, including trace and no-reaction outcomes.","evidence_ids":ev(3)},
        {"region":"scope/mechanism prose","finding":"Scope identities and stereochemical assignment discussion cross-linked to Figure 5 records; proposed mechanism continues on page 4.","evidence_ids":ev(3)}],"record_ids":sorted(page_records[3]),"unresolved":["Table 1 does not provide absolute configuration descriptors for 5a."]},
    {"page":4,"state":"attempted","regions":[
        {"region":"Figure 5","finding":"Sixteen concrete substrate-scope experiments 5b–5q with yields/ee and rr where reported.","evidence_ids":ev(4)},
        {"region":"Mechanism Insight prose","finding":"Three real mechanistic measurement records (31P NMR, reaction orders, KIE) distinguished from author interpretation.","evidence_ids":ev(4)},
        {"region":"Figure 6","finding":"Proposed catalytic cycle, classified as proposed mechanism rather than experimental proof.","evidence_ids":ev(4)},
        {"region":"Scheme 1","finding":"Seven synthetic transformations in the total synthesis, including a 3-step overall yield and a two-product cyclization.","evidence_ids":ev(4)}],"record_ids":sorted(page_records[4]),"unresolved":["Most complex structures lack machine-verifiable SMILES; Figure 5 non-5g configurations are assignments by analogy.","Mechanistic raw data and SI-only procedures are unavailable."]},
    {"page":5,"state":"attempted","regions":[
        {"region":"left-column Scheme 1 continuation","finding":"Narrative confirms transformations already recorded from Scheme 1; no additional reaction arrow.","evidence_ids":ev(5)},
        {"region":"Associated Content","finding":"SI is stated to exist but was not supplied or associated; no SI facts imported.","evidence_ids":ev(5)},
        {"region":"remaining page","finding":"Conclusions, author information, acknowledgments, and references; no new concrete experiment.","evidence_ids":ev(5)}],"record_ids":sorted(page_records[5]),"unresolved":["Experimental procedures and characterization remain missing because SI was not provided."]},
    {"page":6,"state":"attempted","regions":[
        {"region":"entire page","finding":"References only; cited papers were not converted into experiments performed in this article.","evidence_ids":ev(6)}],"record_ids":[],"unresolved":[]},
]}
write(RUN / "coverage/pages_0001_0006.json", coverage)

tasks=json.loads((RUN/"tasks.json").read_text())
for task in tasks:
    if task["document_id"]==DOC and any(task["evidence_id"]==f"{DOC}:page{p}:text" for p in range(1,7)):
        task["state"]="attempted"
write(RUN/"tasks.json",tasks)

print(json.dumps({"molecules":len(bundle["molecules"]),"reactions":len(bundle["reactions"]),"source_context":len(bundle["source_context"])},ensure_ascii=False))
