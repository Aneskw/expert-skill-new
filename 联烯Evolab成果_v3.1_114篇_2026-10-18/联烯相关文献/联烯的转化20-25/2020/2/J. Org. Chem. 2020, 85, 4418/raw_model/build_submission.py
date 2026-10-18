#!/usr/bin/env python3
"""Build the page-1--6 v3.1 machine submission from the inspected source images."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DOC = "doc_2878b0a078d500b2"
EV = {p: f"{DOC}:page{p}:image" for p in range(1, 7)}
TXT = {p: f"{DOC}:page{p}:text" for p in range(1, 7)}


def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def fx(raw=None, normalized=None, evidence=(), method="original_extraction", derivation=None,
       status=None):
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation, status=status)


def setf(obj, key, raw, normalized, evidence, method="original_extraction", derivation=None,
         status=None):
    obj[key] = fx(raw, normalized, evidence, method, derivation, status)


documents = json.loads((RUN / "dataset.json").read_text())["documents"]
bundle = st.empty(documents)
molecules = {}
stereo_diagnostics = []


def mol(label, evidence, name=None, smiles=None, chirality="unspecified", geometry="unspecified",
        representation=None, note=None, assignment=None, mixture=None):
    mid = f"{DOC}:mol:{label}"
    if mid in molecules:
        m = molecules[mid]
        existing = {x for o in m["occurrences"] for x in o["evidence_ids"]}
        new = [x for x in evidence if x not in existing]
        if new:
            m["occurrences"].append({"document_id": DOC,
                "compound_label": fx(label, label, new), "evidence_ids": new})
        return mid
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = fx(label, label, evidence)
    display = name or label
    m["names"]["original"] = fx(display, display, evidence)
    m["names"]["english"] = fx(display, display, evidence)
    m["names"]["abbreviations"] = fx(label, [label], evidence)
    m["names"]["aliases"] = fx()
    m["cas_numbers"] = fx()
    m["identity_note"] = fx(note, note, evidence) if note else fx()
    m["stereochemistry"]["chirality_status"] = fx(chirality, chirality, evidence)
    m["stereochemistry"]["geometry_status"] = fx(geometry, geometry, evidence)
    rep = representation or ("partial" if smiles else "unresolved")
    m["stereochemistry"]["representation_status"] = fx(rep, rep, evidence)
    if assignment:
        m["stereochemistry"]["assignment"] = fx(assignment, assignment, evidence)
    if mixture:
        m["stereochemistry"]["mixture"] = fx(mixture, mixture, evidence)
    if smiles:
        m["representation_type"] = "smiles"
        deriv = "Connectivity candidate derived from the source name/drawing; parsing is not source verification."
        m["smiles"] = fx(None, smiles, evidence, method="structure_derivation", derivation=deriv)
        diag = st.stereo_inspect(smiles)
        stereo_diagnostics.append({"molecule_id": mid, "input_smiles": smiles, **diag})
        if rep == "encoded":
            iso = diag["canonical_isomeric_smiles"]
            m["stereochemistry"]["isomeric_smiles"] = fx(None, iso, evidence,
                method="structure_derivation", derivation="Canonical isomeric SMILES from standard stereo graph diagnostics.")
            m["stereochemistry"]["identity_key"] = fx(None, diag["stereo_graph_key"], evidence,
                method="structure_derivation", derivation="Stereo graph key from standard stereo diagnostics; not a sample identity.")
        else:
            m["stereochemistry"]["isomeric_smiles"] = fx()
            m["stereochemistry"]["identity_key"] = fx()
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = fx()
    molecules[mid] = m
    return mid


# Named/labelled source molecules. Product structures remain unresolved unless the source name is enough
# for a conservative connectivity candidate. All allene axes and unspecified stereocentres stay unspecified.
S = {}
S["1a"] = mol("1a", [EV[3], EV[4]], "N-(But-3-yn-1-yl)-4-methyl-N-(propa-1,2-dien-1-yl)benzenesulfonamide",
    "C#CCCN(C=C=C)S(=O)(=O)c1ccc(C)cc1")
S["1b"] = mol("1b", [EV[5]], "N-(But-3-yn-1-yl)-N-(propa-1,2-dien-1-yl)benzenesulfonamide",
    "C#CCCN(C=C=C)S(=O)(=O)c1ccccc1")
S["1c"] = mol("1c", [EV[5]], "4-Bromo-N-(but-3-yn-1-yl)-N-(propa-1,2-dien-1-yl)benzenesulfonamide",
    "C#CCCN(C=C=C)S(=O)(=O)c1ccc(Br)cc1")
S["1d"] = mol("1d", [EV[5]], "N-(But-3-yn-1-yl)-4-nitro-N-(propa-1,2-dien-1-yl)benzenesulfonamide",
    "C#CCCN(C=C=C)S(=O)(=O)c1ccc([N+](=O)[O-])cc1")
S["1e"] = mol("1e", [EV[5]], "N-(But-3-yn-1-yl)-N-(propa-1,2-dien-1-yl)methanesulfonamide",
    "C#CCCN(C=C=C)S(C)(=O)=O")
S["1f"] = mol("1f", [EV[5]], "4-Methyl-N-(2-methylbut-3-yn-1-yl)-N-(propa-1,2-dien-1-yl)benzenesulfonamide",
    "C#CC(C)CN(C=C=C)S(=O)(=O)c1ccc(C)cc1")
S["1g"] = mol("1g", [EV[5]], "4-Methyl-N-(1-phenylbut-3-yn-1-yl)-N-(propa-1,2-dien-1-yl)benzenesulfonamide",
    "C#CCC(c1ccccc1)N(C=C=C)S(=O)(=O)c1ccc(C)cc1")
S["1h"] = mol("1h", [EV[5]], "N-(2-Ethynylcyclohexyl)-4-methyl-N-(propa-1,2-dien-1-yl)benzenesulfonamide",
    "C#CC1CCCCC1N(C=C=C)S(=O)(=O)c1ccc(C)cc1", chirality="mixture", representation="partial",
    mixture="mixture of two inseparable diastereoisomers; dr approximately 4:1")
S["1i"] = mol("1i", [EV[6]], "N-(But-3-yn-1-yl)-N-(buta-1,2-dien-1-yl)-4-methylbenzenesulfonamide",
    "C#CCCN(C=C=CC)S(=O)(=O)c1ccc(C)cc1", representation="partial",
    note="The substituted allene may carry an axis; no absolute axial assignment is reported in the inspected pages.")
S["4"] = mol("4", [EV[6]], "4-Methyl-N-(pent-4-yn-1-yl)-N-(propa-1,2-dien-1-yl)benzenesulfonamide",
    "C#CCCCN(C=C=C)S(=O)(=O)c1ccc(C)cc1")
S["6"] = mol("6", [EV[6]], "N-(Hex-5-yn-1-yl)-4-methyl-N-(propa-1,2-dien-1-yl)benzenesulfonamide",
    "C#CCCCCN(C=C=C)S(=O)(=O)c1ccc(C)cc1")
S["2a"] = mol("2a", [EV[3], EV[4], EV[5], EV[6]], "methyl phenyldiazoacetate",
    "COC(=O)C(=[N+]=[N-])c1ccccc1", chirality="not_applicable", geometry="not_applicable", representation="partial")

for lab in [f"2{x}" for x in "bcdefghijklmnopqrst"]:
    S[lab] = mol(lab, [EV[4]], lab, note="Structure is drawn in Scheme 2 through variable substituents; no complete SMILES was promoted.")
for lab in ["2_acceptor_H", "2_acceptor_acceptor"]:
    S[lab] = mol(lab, [EV[3]], lab.replace("_", " "), note="Only a diazo-substrate class is reported; exact identity is unreported.")
for lab in [f"3a{x}" for x in "abcdefghijklmnopqrst"]:
    S[lab] = mol(lab, [EV[4]], lab, note="Product drawing/label in Scheme 2; concrete SMILES not independently verified.")
for lab in ["3ba", "3ca", "3da", "3ea", "3fa", "3ga", "3ha", "3ia"]:
    kwargs = {}
    if lab in ("3fa", "3ga"):
        kwargs = {"chirality":"mixture", "mixture":"source reports an inseparable mixture and a printed dr value", "representation":"unresolved"}
    if lab in ("3ha", "3ia"):
        kwargs = {"chirality":"racemic", "assignment":"racemic product; relative stereochemistry established by X-ray according to the text", "representation":"not_encoded"}
    S[lab] = mol(lab, [EV[5]], ("rac-" if lab in ("3ha", "3ia") else "") + lab, **kwargs)
S["3_acetyl"] = mol("3_acetyl", [EV[3], EV[5]], "expected product from the R=acetyl allenyne", note="Expected product only; not observed.")
S["1_acetyl"] = mol("1_acetyl", [EV[3], EV[5]], "1,6-allenyne with R=acetyl", note="Exact compound label is not printed in the inspected scheme.")
for lab in [f"5{x}" for x in "abcdefghi"]:
    S[lab] = mol(lab, [EV[6]], lab, note="Product drawing/label in Scheme 4; alkene geometry is not textually assigned on the inspected page.")
S["7"] = mol("7", [EV[6]], "cross-coupling bisallene 7", geometry="unspecified", representation="unresolved")

common_names = {
    "CuPPh3Br":"Cu(PPh3)3Br", "CuCl":"CuCl", "CuBr":"CuBr", "CuI":"CuI",
    "CuPhenPPh3Br":"Cu(1,10-Phen)(PPh3)Br", "CuOTf":"CuOTf", "Cutfacac2":"Cu(tfacac)2",
    "IPrCuOTf":"IPrCuOTf", "CH2Cl2":"CH2Cl2", "DCE":"1,2-dichloroethane (DCE)",
    "MeCN":"acetonitrile (MeCN)", "THF":"tetrahydrofuran (THF)", "Et3N":"Et3N",
    "Cs2CO3":"Cs2CO3", "BHT":"BHT", "TEMPO":"TEMPO", "K2CO3":"K2CO3",
    "acetone":"acetone", "propargyl_bromide":"propargyl bromide", "1_bromo_2_butyne":"1-bromo-2-butyne",
    "tBuOK":"tBuOK", "N2":"N2", "Celite":"Celite", "hexanes":"hexanes", "EtOAc":"EtOAc"
}
for k, n in common_names.items():
    S[k] = mol(k, [EV[3], EV[4], EV[5], EV[6]], n, chirality="not_applicable", geometry="not_applicable")


def participant(pid, mid, role, evidence, stage="reaction", amount=None, equivalents=None,
                mol_percent=None, reference=None, mass=None, volume=None):
    p = st.new_participant(pid, mid, role, evidence, stage)
    if amount is not None: p["charges"]["amount"] = fx(str(amount), amount, evidence)
    if mass is not None: p["charges"]["mass"] = fx(str(mass), mass, evidence)
    if volume is not None: p["charges"]["volume"] = fx(str(volume), volume, evidence)
    if equivalents is not None: p["charges"]["equivalents"] = fx(str(equivalents), equivalents, evidence)
    if mol_percent is not None: p["charges"]["mol_percent"] = fx(str(mol_percent), mol_percent, evidence)
    if reference:
        p["equivalence_reference"] = fx(reference, reference, evidence)
    return p


def step(step_id, stage, pids, evidence, operation, temperature=None, time=None, atmosphere=None):
    s = st.new_step(step_id, evidence, stage)
    s["participant_ids"] = list(pids)
    s["operation"] = fx(operation, operation, evidence)
    if temperature is not None: s["temperature"] = fx(str(temperature), temperature, evidence)
    if time is not None: s["time"] = fx(str(time), time, evidence)
    if atmosphere is not None: s["atmosphere"] = fx(str(atmosphere), atmosphere, evidence)
    return s


def classification(evidence, codes=(), count_by_code=None, coverage="partial", limitation=None):
    count_by_code = count_by_code or {}
    changes = []
    labels = []
    number = {x["code"]: x["number"] for x in json.loads((ROOT / "skills/chemical-reaction-extraction/references/bond_taxonomy.json").read_text())["categories"]}
    elems = {"C-C":("C","C",1), "C-N":("C","N",1), "DB-C-C":("C","C",2)}
    for code in codes:
        support = []
        for _ in range(count_by_code.get(code, 1)):
            cid = f"b{len(changes)+1}"
            a, b, order = elems[code]
            ch = {"change_id":cid, "atoms":[{"atom_map":None,"element":a},{"atom_map":None,"element":b}],
                  "before_order":0 if not code.startswith("DB-") else 1, "after_order":order,
                  "change_type":"formed" if order == 1 else "order_increased"}
            changes.append(fx(None, ch, evidence, method="structure_derivation",
                derivation="Partial net bond inventory inferred from the inspected source transformation drawing; atom mapping is unresolved."))
            support.append(cid)
        labels.append({"category_number":number[code], "code":code, "supporting_change_ids":support})
    return {"standard_id":st.STANDARD, "mapped_reaction_smiles":fx(), "bond_changes":changes,
        "labels":fx(None, labels, evidence, method="structure_derivation",
            derivation="Labels follow the 20261016 taxonomy and the partial source-drawing bond inventory."),
        "coverage":fx(None, coverage, evidence, method="inference",
            derivation="No complete source-verified atom mapping is available; classification is partial."),
        "limitations":[limitation or "Partial connectivity-level classification; implicit-H and full atom correspondence remain unresolved."],
        "model_votes":[]}


def reaction(rid, evidence, record_type="specific_experiment", categories=("other",), name=None,
             transform=None, location=None, class_codes=(), class_counts=None, class_coverage="partial"):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = fx(list(categories), list(categories), evidence)
    r["location"] = fx(location or rid, {"pages":sorted({int(x.split(":page")[1].split(":")[0]) for x in evidence}), "region":location or rid}, evidence)
    if name: r["reaction_name"] = fx(name, name, evidence)
    if transform:
        r["transformation_type"] = fx(None, transform, evidence, method="inference", derivation="Concise transformation description inferred from source drawing/text.")
    r["classification"] = classification(evidence, class_codes, class_counts, class_coverage)
    r["review"] = {"status":"needs_review", "method":"model", "reviewer":None, "checked_fields":[],
        "unresolved":["Machine record; chemical identity and structures require expert review."]}
    r["completeness"] = "partial"
    return r


def outcome(r, oid, kind, targets, evidence, raw, normalized, basis=None, yield_type=None, scope=None):
    o = st.new_outcome(oid, kind, targets, evidence)
    o["value"] = fx(raw, normalized, evidence)
    if basis is not None: o["measurement_basis"] = fx(basis, basis, evidence)
    if yield_type is not None: o["yield_type"] = fx(yield_type, yield_type, evidence)
    if scope is not None: o["scope"] = fx(scope, scope, evidence)
    r["outcomes"].append(o)


def common_cascade(r, substrate, diazo, product, evidence, solvent="CH2Cl2", temperature="r.t.", catalyst="CuPPh3Br",
                   scale=True, syringe=True, additives=()):
    sub = S[substrate]; dz = S[diazo]; prod = S.get(product)
    r["participants"] = [
        participant("substrate", sub, "reactant", evidence, amount={"number":0.2,"unit":"mmol"} if scale else None,
                    equivalents={"number":1,"unit":"equiv"} if scale else None, reference=sub if scale else None),
        participant("diazo", dz, "reactant", evidence, amount={"number":0.4,"unit":"mmol"} if scale else None,
                    equivalents={"number":2,"unit":"equiv"} if scale else None, reference=sub if scale else None),
    ]
    if catalyst:
        r["participants"].append(participant("catalyst", S[catalyst], "catalyst", evidence,
            mol_percent={"number":10,"unit":"mol%"}, reference=sub))
    r["participants"].append(participant("solvent", S[solvent], "solvent", evidence))
    for i,(lab,amt,role) in enumerate(additives):
        r["participants"].append(participant(f"additive{i+1}", S[lab], role, evidence,
            amount=amt if isinstance(amt,dict) and amt.get("unit") == "mmol" else None,
            equivalents=amt if isinstance(amt,dict) and amt.get("unit") == "equiv" else None,
            reference=sub))
    if prod:
        r["participants"].append(participant("product", prod, "product", evidence))
    r["reaction_scale"] = fx("0.2 mmol allenyne", {"reference_molecule_id":sub,"amount":{"number":0.2,"unit":"mmol"}}, evidence) if scale else fx()
    r["solvent_conditions"]["composition"] = fx(solvent, [solvent], evidence)
    if scale:
        r["solvent_conditions"]["total_volume"] = fx("1.0 mL + 1.0 mL", {"number":2.0,"unit":"mL"}, evidence,
            method="inference", derivation="Sum of the two explicitly reported solvent portions; not an independently measured total volume.")
    pids = [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"]
    if syringe:
        r["steps"] = [
            step("addition", "reaction", pids, evidence, "inject diazo solution into allenyne/catalyst solution via syringe pump", temperature, {"number":2,"unit":"h"}),
            step("hold", "reaction", pids, evidence, "continue reaction", temperature, {"number":12,"unit":"h"}),
        ]
        r["procedure"]["addition_order"] = fx("diazo solution injected into allenyne/catalyst solution", "diazo solution -> allenyne/catalyst solution", evidence)
    else:
        r["steps"] = [step("one_pot_charge", "reaction", pids, evidence, "one-pot reaction without syringe pump", temperature),
                      step("hold", "reaction", pids, evidence, "continue reaction", temperature, {"number":12,"unit":"h"})]
        r["procedure"]["modifications"] = fx("without a syringe pump", "no syringe pump", evidence)


reactions = []

# Scheme 1: eight cited prior-work arrows, each retained as a literature example; structures are generic.
lit_specs = [
    ("scheme1a_Au", "Scheme 1a left [Au] arrow", "[Au]-promoted cyclization"),
    ("scheme1a_I", "Scheme 1a right I+ arrow", "iodine electrophilic cyclization"),
    ("scheme1b_Pd", "Scheme 1b left [Pd] arrow", "[2+2] cycloaddition of bisallene"),
    ("scheme1b_heat", "Scheme 1b right heat arrow", "thermal cycloaddition of bisallene"),
    ("scheme1c_Au", "Scheme 1c [Au] arrow to bisallene", "in situ bisallene formation"),
    ("scheme1c_cyclize", "Scheme 1c bisallene-to-product arrow", "cyclization of in situ bisallene"),
    ("scheme1d_cross_coupling", "Scheme 1d [Cu] arrow to bisallene", "previously reported copper cross-coupling"),
    ("scheme1d_cyclize", "Scheme 1d bisallene-to-product arrow", "previously reported intramolecular [2+2] cycloaddition"),
]
for rid, loc, trans in lit_specs:
    ev = [EV[2]]
    a = mol(rid+"_reactant", ev, rid+" generic reactant", note="Generic R-group drawing from cited prior work.")
    b = mol(rid+"_product", ev, rid+" generic product", note="Generic R-group drawing from cited prior work.")
    r = reaction(rid, ev, "literature_example", ("other",), "cited prior-work transformation", trans, loc, (), None, "unresolved")
    r["participants"] = [participant("reactant", a, "reactant", ev), participant("product", b, "product", ev)]
    r["steps"] = [step("scheme_arrow", "reaction", ["reactant","product"], ev, loc)]
    r["review"]["unresolved"].append("Cited literature example only; conditions and exact variable structures were not imported as this paper's experiments.")
    reactions.append(r)

# Generic current-work arrow, shown graphically on pages 1 and 2.
generic_in = mol("generic_1n_allenyne", [EV[1],EV[2]], "generic 1,n-allenyne (x=0 or 1)", note="R-group template drawing; left unresolved rather than expanded.")
generic_dz = mol("generic_aryl_diazoacetate", [EV[1],EV[2]], "generic aryl diazoacetate")
generic_out = mol("generic_azabicyclo_product", [EV[1],EV[2]], "generic 3-azabicyclo[m.2.0] product")
r = reaction("scheme1e_generic", [EV[1],EV[2]], "generic_scheme", ("other",),
    "copper-catalyzed cascade reaction of diazo compounds with 1,n-allenynes",
    "tandem cross-coupling/[2+2] cycloaddition", "page 1 graphical abstract / Scheme 1e", ("C-C",), {"C-C":3})
r["participants"] = [participant("allenyne",generic_in,"reactant",r["evidence_ids"]), participant("diazo",generic_dz,"reactant",r["evidence_ids"]),
    participant("copper",S["CuPPh3Br"],"catalyst",r["evidence_ids"]), participant("product",generic_out,"product",r["evidence_ids"])]
r["steps"] = [step("scheme_arrow","reaction",["allenyne","diazo","copper","product"],r["evidence_ids"],"generic copper-catalyzed cascade")]
reactions.append(r)

# Table 1 optimization rows.
opts = {
  1:(None,"CH2Cl2",None,"no_reaction",None,True,()), 2:("CuCl","CH2Cl2",45,"yield",None,True,()),
  3:("CuBr","CH2Cl2",27,"yield",None,True,()), 4:("CuI","CH2Cl2",None,"no_reaction",None,True,()),
  5:("CuPPh3Br","CH2Cl2",83,"yield",None,True,()), 6:("CuPhenPPh3Br","CH2Cl2",None,"no_reaction",None,True,()),
  7:("CuOTf","CH2Cl2",None,"no_reaction",None,True,()), 8:("Cutfacac2","CH2Cl2",None,"no_reaction",None,True,()),
  9:("IPrCuOTf","CH2Cl2",None,"no_reaction",None,True,()), 10:("CuPPh3Br","DCE",32,"yield",None,True,()),
  11:("CuPPh3Br","MeCN",None,"no_reaction",None,True,()), 12:("CuPPh3Br","THF","<10","yield",None,True,()),
  13:("CuPPh3Br","CH2Cl2",None,"no_reaction",{"number":0,"unit":"°C"},True,()),
  14:("CuPPh3Br","CH2Cl2",46,"yield",{"number":60,"unit":"°C"},True,()),
  15:("CuPPh3Br","CH2Cl2","<10","yield",None,False,()),
  16:("CuPPh3Br","CH2Cl2",37,"yield",None,True,(("Et3N",{"number":0.2,"unit":"mmol"},"additive"),)),
  17:("CuPPh3Br","CH2Cl2",None,"no_reaction",None,True,(("Cs2CO3",{"number":0.2,"unit":"mmol"},"additive"),)),
}
for n,(cat,solv,val,kind,temp,syr,adds) in opts.items():
    ev=[EV[3]]; r=reaction(f"table1_entry_{n:02d}",ev,"specific_experiment",("optimization",),
        "optimization of the copper-catalyzed cascade","cascade cross-coupling/[2+2] cycloaddition",f"Table 1 entry {n}",
        ("C-C",) if kind=="yield" else (), {"C-C":3} if kind=="yield" else None,
        "partial" if kind=="yield" else "not_applicable")
    common_cascade(r,"1a","2a","3aa",ev,solv,temp or "r.t.",cat,True,syr,adds)
    if kind=="yield":
        q={"number":val,"unit":"%"} if isinstance(val,(int,float)) else {"number":10,"unit":"%","qualifier":"lt","raw_text":"<10%"}
        outcome(r,"yield","yield",[S["3aa"]],ev,f"{val}%" if val != "<10" else "<10%",q)
    else:
        outcome(r,"no_reaction","no_reaction",[],ev,"no reaction / no catalytic activity",True)
    r["review"]["unresolved"].append("Yield measurement basis is not stated in Table 1; blank cells are represented as no reaction only where the page-2/3 prose explicitly says so.")
    reactions.append(r)

# Scheme 2 scope and its explicit 5 mmol repeat.
scope2 = {"3aa":83,"3ab":81,"3ac":70,"3ad":65,"3ae":87,"3af":85,"3ag":78,"3ah":91,"3ai":71,"3aj":68,
          "3ak":66,"3al":61,"3am":69,"3an":53,"3ao":23,"3ap":60,"3aq":63,"3ar":58,"3as":69,"3at":0}
for i,(prod,y) in enumerate(scope2.items()):
    diazo = "2" + chr(ord("a")+i)
    ev=[EV[4]]; r=reaction(f"scheme2_{prod}",ev,"specific_experiment",("substrate_scope",),
        "scope of diazo substrates","cascade cross-coupling/[2+2] cycloaddition",f"Scheme 2 {prod}",
        ("C-C",) if y else (), {"C-C":3} if y else None, "partial" if y else "not_applicable")
    common_cascade(r,"1a",diazo,prod if y else None,ev)
    if y:
        outcome(r,"yield","yield",[S[prod]],ev,f"{y}%",{"number":y,"unit":"%"})
        r["review"]["unresolved"].append("Scheme 2 does not state the yield measurement basis on the inspected page.")
    else:
        outcome(r,"no_reaction","no_reaction",[],ev,"0%",True)
    reactions.append(r)

r=reaction("scheme2_3aa_scaleup",[EV[3],EV[4]],"specific_experiment",("scale_up",),"gram-scale synthesis of 3aa",
    "cascade cross-coupling/[2+2] cycloaddition","Scheme 2 parenthetical / page-3 prose",("C-C",),{"C-C":3})
common_cascade(r,"1a","2a","3aa",r["evidence_ids"],scale=False)
r["reaction_scale"] = fx("5.0 mmol scale",{"reference_molecule_id":S["1a"],"amount":{"number":5.0,"unit":"mmol"}},r["evidence_ids"])
outcome(r,"yield","yield",[S["3aa"]],r["evidence_ids"],"72%",{"number":72,"unit":"%"})
r["review"]["unresolved"] += ["Only the 5.0 mmol scale and 72% yield are available within pages 1-6; detailed scale-up solvent volumes are not imported from later pages.","Yield basis is not stated in the inspected region."]
reactions.append(r)

# Two unsuccessful diazo classes explicitly discussed on page 3.
for lab in ("2_acceptor_H","2_acceptor_acceptor"):
    ev=[EV[3]]; r=reaction(f"scope_failure_{lab}",ev,"specific_experiment",("substrate_scope","control"),
        "unsuitable diazo-substrate class","attempted cascade reaction",f"page 3: {lab} not suitable",(),None,"not_applicable")
    common_cascade(r,"1a",lab,None,ev)
    outcome(r,"no_reaction","no_reaction",[],ev,"not suitable for this reaction",True)
    r["review"]["unresolved"].append("Exact diazo compound identity and analytical detection method are unreported.")
    reactions.append(r)

# GP-A generic record and eleven concrete product preparations.
gp_ev=[EV[4],EV[5],EV[6]]
gp_r=reaction("GP_A_generic",gp_ev,"generic_scheme",("other",),"General Procedure for the Synthesis of 1,n-Allenyne Substrates 1, 4, and 6",
    "N-propargylation followed by base-promoted propargyl-to-allene isomerization","page 4 GP-A",("C-N","DB-C-C"))
gp_prec=mol("GP_A_aminoalkyne",gp_ev,"corresponding aminoalkyne",note="Generic precursor varies by product.")
gp_prod=mol("GP_A_allenyne",gp_ev,"generic 1,n-allenyne product",note="Generic product class.")
gp_r["participants"]=[participant("aminoalkyne",gp_prec,"reactant",gp_ev),participant("propargyl_bromide",S["propargyl_bromide"],"reactant",gp_ev),
    participant("K2CO3",S["K2CO3"],"reagent",gp_ev),participant("acetone",S["acetone"],"solvent",gp_ev),participant("tBuOK",S["tBuOK"],"reagent",gp_ev),
    participant("THF",S["THF"],"solvent",gp_ev),participant("product",gp_prod,"product",gp_ev)]
gp_r["steps"]=[step("alkylation","reaction",["aminoalkyne","propargyl_bromide","K2CO3","acetone"],gp_ev,"add propargyl bromide in one portion; stir until TLC completion",{"number":60,"unit":"°C"},"until completed (monitored by TLC)","N2 balloon"),
    step("isomerization_addition","reaction",["tBuOK","THF"],gp_ev,"cool crude propargylated product and add tBuOK in portions",{"number":0,"unit":"°C"}),
    step("isomerization_hold","reaction",["tBuOK","THF","product"],gp_ev,"warm to room temperature and stir until TLC completion","r.t.","until completed (monitored by TLC)"),
    step("purification","purification",["product"],gp_ev,"flash column chromatography on silica gel; hexanes/EtOAc eluent")]
gp_r["procedure"]["original_text"] = fx("GP-A on page 4", "two-stage GP-A with filtration/concentration between stages", gp_ev)
reactions.append(gp_r)

prep_data={"1a":(4,60,1568.0,10,"propargyl_bromide"),"1b":(4,71,878.0,5,"propargyl_bromide"),"1c":(5,74,1207.0,5,"propargyl_bromide"),
 "1d":(5,38,555.0,5,"propargyl_bromide"),"1e":(5,66,611.2,5,"propargyl_bromide"),"1f":(5,69,950.0,5,"propargyl_bromide"),
 "1g":(5,57,961.6,5,"propargyl_bromide"),"1h":(5,86,1356.0,5,"propargyl_bromide"),"1i":(6,47,647.2,5,"1_bromo_2_butyne"),
 "4":(6,44,1212.0,10,"propargyl_bromide"),"6":(6,23,665.6,10,"propargyl_bromide")}
for lab,(pg,y,mg,scale,bromide) in prep_data.items():
    ev=[EV[pg],EV[4]] if pg!=4 else [EV[4]]
    prec=mol(f"aminoalkyne_for_{lab}",ev,f"corresponding aminoalkyne precursor for {lab}",note="Identity is referenced to literature and not fully named on the inspected page.")
    r=reaction(f"prep_{lab}",ev,"specific_experiment",("other",),f"preparation of allenyne {lab}",
        "N-alkylation/propargyl-to-allene isomerization",f"individual paragraph for {lab} following GP-A",("C-N","DB-C-C"))
    r["general_scheme_ids"]=[gp_r["reaction_id"]]
    r["participants"]=[participant("aminoalkyne",prec,"reactant",ev,equivalents={"number":1,"unit":"equiv"},reference=prec),
      participant("bromide",S[bromide],"reactant",ev,equivalents={"number":2.5,"unit":"equiv"},reference=prec),
      participant("K2CO3",S["K2CO3"],"reagent",ev,equivalents={"number":5,"unit":"equiv"},reference=prec),
      participant("acetone",S["acetone"],"solvent",ev),participant("tBuOK",S["tBuOK"],"reagent",ev,mol_percent={"number":30,"unit":"mol%"},reference=prec),
      participant("THF",S["THF"],"solvent",ev),participant("N2",S["N2"],"additive",ev),
      participant("Celite",S["Celite"],"workup_agent",ev,"workup"),participant("CH2Cl2",S["CH2Cl2"],"workup_agent",ev,"workup"),
      participant("hexanes",S["hexanes"],"purification_agent",ev,"purification"),participant("EtOAc",S["EtOAc"],"purification_agent",ev,"purification"),
      participant("product",S[lab],"product",ev)]
    r["reaction_scale"] = fx(f"{scale} mmol scale",{"reference_molecule_id":prec,"amount":{"number":scale,"unit":"mmol"}},ev)
    r["solvent_conditions"]["composition"] = fx("acetone (0.5 M), then THF (0.33 M)",["acetone","THF"],ev)
    r["steps"]=[step("alkylation","reaction",["aminoalkyne","bromide","K2CO3","acetone"],ev,"add bromide in one portion and stir until TLC completion",{"number":60,"unit":"°C"},"until completed (monitored by TLC)","N2 balloon"),
      step("first_workup","workup",["Celite"],ev,"cool to room temperature; filter through Celite, wash with acetone, concentrate"),
      step("base_addition","reaction",["tBuOK","THF"],ev,"dissolve crude in THF, cool, add tBuOK in portions",{"number":0,"unit":"°C"}),
      step("isomerization","reaction",["tBuOK","THF","product"],ev,"warm to room temperature and stir until TLC completion","r.t.","until completed (monitored by TLC)"),
      step("second_workup","workup",["Celite","CH2Cl2"],ev,"filter through Celite, wash with CH2Cl2, concentrate"),
      step("purification","purification",["hexanes","EtOAc","product"],ev,"flash column chromatography on silica gel (hexanes/EtOAc)")]
    r["procedure"]["general_procedure_reference"] = fx("following GP-A",gp_r["reaction_id"],ev)
    outcome(r,"yield","yield",[S[lab]],ev,f"{y}%",{"number":y,"unit":"%"},basis="isolated after flash column chromatography",yield_type="isolated")
    outcome(r,"product_mass","product_amount",[S[lab]],ev,f"{mg} mg",{"number":mg,"unit":"mg"},basis="reported isolated product amount")
    if lab=="1h":
        outcome(r,"dr","selectivity",[S[lab]],ev,"dr ~ 4:1","~4:1",basis="1H NMR analysis",scope={"metric":"dr","order":"unassigned diastereomers"})
    r["review"]["unresolved"].append("Reaction times and exact chromatography ratio are unreported in GP-A/individual paragraph.")
    reactions.append(r)

# Scheme 3 successful substrate scope plus the explicit acetyl failure.
scope3={"3ba":80,"3ca":76,"3da":71,"3ea":73,"3fa":82,"3ga":73,"3ha":35,"3ia":74}
sub_for={"3ba":"1b","3ca":"1c","3da":"1d","3ea":"1e","3fa":"1f","3ga":"1g","3ha":"1h","3ia":"1i"}
for prod,y in scope3.items():
    ev=[EV[5]]; r=reaction(f"scheme3_{prod}",ev,"specific_experiment",("substrate_scope",),"scope of 1,6-allenyne substrates",
        "cascade cross-coupling/[2+2] cycloaddition",f"Scheme 3 {prod}",("C-C",),{"C-C":3})
    common_cascade(r,sub_for[prod],"2a",prod,ev)
    basis="flash chromatography followed by recrystallization" if prod=="3ha" else None
    outcome(r,"yield","yield",[S[prod]],ev,f"{y}%",{"number":y,"unit":"%"},basis=basis,yield_type="isolated" if basis else None)
    if prod in ("3fa","3ga"):
        ratio="~1.5:1" if prod=="3fa" else "~2:1"
        outcome(r,"dr","selectivity",[S[prod]],ev,f"dr {ratio}",ratio,basis="1H NMR analysis",scope={"metric":"dr","order":"unassigned components"})
        r["review"]["unresolved"].append("Text calls the products inseparable regioisomers while the scheme labels the ratio as dr; the printed terminology is retained without repair.")
    if prod in ("3ha","3ia"):
        r["notes"].append(fx("one diastereomer was obtained", "one diastereomer obtained", ev))
    if basis is None: r["review"]["unresolved"].append("Yield measurement basis is not stated in Scheme 3 for this entry.")
    reactions.append(r)

ev=[EV[3],EV[5]]; r=reaction("scheme3_acetyl_failure",ev,"specific_experiment",("substrate_scope","control"),"R=acetyl allenyne scope failure",
    "attempted cascade reaction", "page-3 discussion of Scheme 3",(),None,"not_applicable")
common_cascade(r,"1_acetyl","2a",None,ev)
outcome(r,"no_reaction","no_reaction",[],ev,"no desired product was obtained",True)
r["review"]["unresolved"].append("Exact substrate label and detection method are unreported.")
reactions.append(r)

# Scheme 4 scope.
scope4={"5a":61,"5b":65,"5c":68,"5d":60,"5e":58,"5f":55,"5g":52,"5h":38,"5i":47}
for i,(prod,y) in enumerate(scope4.items()):
    diazo="2"+chr(ord("a")+i); ev=[EV[6]]
    r=reaction(f"scheme4_{prod}",ev,"specific_experiment",("substrate_scope",),"scope of diazo substrates with 1,7-allenyne",
        "cascade cross-coupling/[2+2] cycloaddition",f"Scheme 4 {prod}",("C-C",),{"C-C":3})
    common_cascade(r,"4",diazo,prod,ev)
    outcome(r,"yield","yield",[S[prod]],ev,f"{y}%",{"number":y,"unit":"%"})
    r["review"]["unresolved"] += ["Yield measurement basis is not stated in Scheme 4.","Alkene E/Z assignments are not textually stated in the inspected scheme and remain unspecified."]
    reactions.append(r)

# Scheme 5 controls.
for additive,y in (("BHT",80),("TEMPO",84)):
    ev=[EV[6]]; r=reaction(f"scheme5_{additive}",ev,"specific_experiment",("control","mechanistic_study"),
        "radical inhibition experiment","cascade reaction in the presence of radical inhibitor",f"Scheme 5a {additive}",("C-C",),{"C-C":3})
    common_cascade(r,"1a","2a","3aa",ev,additives=((additive,{"number":1,"unit":"equiv"},"additive"),))
    outcome(r,"yield","yield",[S["3aa"]],ev,f"{y}%",{"number":y,"unit":"%"})
    r["review"]["unresolved"].append("Yield basis is not stated in Scheme 5a; no radical-trapping detection limit is reported.")
    reactions.append(r)

ev=[EV[6]]; r=reaction("scheme5_bisallene7",ev,"specific_experiment",("control","mechanistic_study"),
    "1,8-allenyne control experiment","cross-coupling to bisallene without cycloaddition","Scheme 5b",("C-C",),{"C-C":1})
common_cascade(r,"6","2a","7",ev)
outcome(r,"yield","yield",[S["7"]],ev,"56%",{"number":56,"unit":"%"})
r["review"]["unresolved"].append("Yield basis is not stated in Scheme 5b on the inspected page.")
reactions.append(r)

# Page-3 prose-only failed derivatizations and proposed mechanism description.
for rid,label in (("failed_reduction_3aa","reduction of a double bond"),("failed_thermal_rearrangement_3aa","thermal rearrangement")):
    ev=[EV[3]]; r=reaction(rid,ev,"specific_experiment",("derivatization","control"),label,"unsuccessful product transformation",f"page 3: {label}",(),None,"not_applicable")
    r["participants"]=[participant("3aa",S["3aa"],"reactant",ev)]
    r["steps"]=[step("attempt","reaction",["3aa"],ev,label+" was attempted; experimental conditions are unreported")]
    outcome(r,"no_reaction","no_reaction",[],ev,"all failed",True)
    r["review"]["unresolved"] += ["No reagents, temperature, duration, or analytical method are reported for this failed effort."]
    reactions.append(r)

r=reaction("proposed_mechanism_text",[EV[3]],"proposed_mechanism",("other",),"plausible mechanism",
    "copper-carbene formation, migratory insertion, protonation to bisallene, cupracycle formation, reductive elimination","page 3 mechanism discussion referring to Scheme 6",("C-C",),{"C-C":3})
r["participants"]=[participant("1a",S["1a"],"reactant",[EV[3]]),participant("2a",S["2a"],"reactant",[EV[3]]),participant("3aa",S["3aa"],"product",[EV[3]])]
r["notes"].append(fx("a plausible mechanism has been proposed", "author-proposed, not experimentally established", [EV[3]]))
r["review"]["unresolved"].append("Scheme 6 itself is on page 7 and is outside this pass; only the page-3 prose was recorded.")
reactions.append(r)

bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
bundle["source_context"] = [
 {"context_id":f"{DOC}:context:page1_natural_products","document_id":DOC,"evidence_ids":[EV[1]],"kind":"nonreaction_figure","note":"Figure 1 selected natural products; not converted to synthesis experiments."},
 {"context_id":f"{DOC}:context:pages4_6_characterization","document_id":DOC,"evidence_ids":[EV[4],EV[5],EV[6]],"kind":"characterization","note":"NMR/HRMS peak lists and melting points retained as characterization context, not separate reactions."},
 {"context_id":f"{DOC}:context:page4_safety","document_id":DOC,"evidence_ids":[EV[4]],"kind":"safety_note","note":"Diazo compounds: fume cupboard and blast shield."},
 {"context_id":f"{DOC}:context:page4_catalyst_diazo_literature","document_id":DOC,"evidence_ids":[EV[4]],"kind":"literature_preparation_reference","note":"Cu(PPh3)3Br and aryldiazoacetates 2a-t are stated to follow cited literature procedures; no new preparation experiment was invented."},
]

dump(RUN / "submitted.json", bundle)
dump(RUN / "raw_model/stereo_diagnostics.json", stereo_diagnostics)
dump(RUN / "raw_model/visual_transcription_pages_0001_0006.json", {
 "document_id":DOC,"pages":[1,2,3,4,5,6],"source":"manual visual inspection of supplied original-detail page images",
 "record_count":len(reactions),"molecule_count":len(molecules),
 "inventory":{"scheme1_prior_arrows":8,"scheme1_current_generic":1,"table1_rows":17,"scheme2_scope":20,
  "scheme2_scaleup":1,"diazo_class_failures":2,"GP_A_generic":1,"GP_A_concrete_preparations":11,
  "scheme3_successes":8,"scheme3_failure":1,"scheme4_scope":9,"scheme5_controls":3,
  "prose_failed_derivatizations":2,"proposed_mechanism_text":1},
 "limitations":["No SI was supplied or associated.","Product SMILES were not promoted when drawings/names were not independently structure-verified.","Pages after 6 were not claimed as visually covered."]})

rid = lambda x: f"{DOC}:{x}"
coverage = {"document_id":DOC,"semantic_accuracy_verified":False,"pages":[
 {"page":1,"state":"attempted","regions":[
   {"region":"graphical abstract at upper right","finding":"Current-work generic copper-catalyzed allenyne/diazo cascade; duplicate evidence for Scheme 1e, not a separate experiment.","evidence_ids":[EV[1]]},
   {"region":"Figure 1 and introduction","finding":"Selected natural products and background discussion; nonreaction context, no synthesis experiment extracted.","evidence_ids":[EV[1]]}],
  "record_ids":[rid("scheme1e_generic")],"unresolved":["Generic R groups and product structure are not expanded to concrete SMILES."]},
 {"page":2,"state":"attempted","regions":[
   {"region":"Scheme 1a-d","finding":"Eight arrows from cited prior work, recorded as literature examples rather than this paper's experiments.","evidence_ids":[EV[2]]},
   {"region":"Scheme 1e","finding":"One current-report generic cascade arrow; merged with page-1 graphical-abstract evidence.","evidence_ids":[EV[2]]},
   {"region":"Results and Discussion lower half","finding":"Optimization narrative identifies no-reaction/no-activity entries and explains solvent, temperature, and syringe-pump effects; linked to page-3 Table 1 records.","evidence_ids":[EV[2]]}],
  "record_ids":[rid(x) for x,_,_ in lit_specs]+[rid("scheme1e_generic")],
  "unresolved":["Cited examples use generic structures and incomplete conditions; no external papers were imported."]},
 {"page":3,"state":"attempted","regions":[
   {"region":"Table 1 including footnotes a-f","finding":"All 17 optimization rows recorded independently; footnote additives and modified temperatures/no-pump condition are structured.","evidence_ids":[EV[3]]},
   {"region":"left-column scope discussion","finding":"5 mmol 3aa repeat, two unsuitable diazo classes, and R=acetyl failure identified; exact identities remain unresolved where not printed.","evidence_ids":[EV[3]]},
   {"region":"right-column controls/mechanism discussion","finding":"Two failed 3aa derivatizations and the author-proposed mechanism prose recorded separately; Scheme 6 drawing is outside this pass.","evidence_ids":[EV[3]]}],
  "record_ids":[rid(f"table1_entry_{n:02d}") for n in range(1,18)]+[rid("scheme2_3aa_scaleup"),rid("scope_failure_2_acceptor_H"),rid("scope_failure_2_acceptor_acceptor"),rid("scheme3_acetyl_failure"),rid("failed_reduction_3aa"),rid("failed_thermal_rearrangement_3aa"),rid("proposed_mechanism_text")],
  "unresolved":["Table 1 and prose do not report yield measurement methods.","Failed derivatization reagents, conditions, and analytical methods are unreported.","Scheme 6 image is on page 7 and was not visually covered."]},
 {"page":4,"state":"attempted","regions":[
   {"region":"Scheme 2 and footnotes","finding":"Twenty 3aa-3at scope arrows/outcomes plus the distinct 5 mmol 3aa repeat; thiophene example is an explicit 0% experiment.","evidence_ids":[EV[4]]},
   {"region":"Experimental Section GP-A","finding":"Two-stage general allenyne synthesis procedure with ordered reaction/workup/purification operations.","evidence_ids":[EV[4]]},
   {"region":"individual paragraphs 1a-1b","finding":"Two concrete substrate preparations with yield, mass, and scale; NMR/HRMS retained as context only.","evidence_ids":[EV[4]]},
   {"region":"Safety/catalyst/diazo-substrate paragraphs","finding":"Safety note and cited literature preparation references; no unsupported catalyst/diazo syntheses invented.","evidence_ids":[EV[4]]}],
  "record_ids":[rid(f"scheme2_3a{x}") for x in "abcdefghijklmnopqrst"]+[rid("scheme2_3aa_scaleup"),rid("GP_A_generic"),rid("prep_1a"),rid("prep_1b")],
  "unresolved":["Scheme 2 yield measurement basis is not stated.","Most product structures remain unverified drawings/labels with no promoted SMILES.","GP-A completion times and chromatography ratio are unreported."]},
 {"page":5,"state":"attempted","regions":[
   {"region":"Scheme 3 and footnotes a-d","finding":"Eight successful 1,6-allenyne scope entries, one prose-linked R=acetyl failure, printed dr values, racemic labels, X-ray notes, and the special recrystallization basis for 3ha.","evidence_ids":[EV[5]]},
   {"region":"individual paragraphs 1c-1h","finding":"Six concrete allenyne preparations; 1h retained as a two-diastereomer mixture with dr approximately 4:1.","evidence_ids":[EV[5]]},
   {"region":"NMR/HRMS blocks","finding":"Characterization only; associated with preparations but not converted into separate reaction records.","evidence_ids":[EV[5]]}],
  "record_ids":[rid("GP_A_generic")]+[rid(f"prep_1{x}") for x in "cdefgh"]+[rid(f"scheme3_{x}") for x in ("3ba","3ca","3da","3ea","3fa","3ga","3ha","3ia")]+[rid("scheme3_acetyl_failure")],
  "unresolved":["Scheme 3 calls 3fa/3ga inseparable regioisomers while printing dr; terminology is preserved as a source inconsistency.","Unspecified stereocentres/axes are not converted to achiral or racemic.","Yield basis remains unreported except 3ha."]},
 {"page":6,"state":"attempted","regions":[
   {"region":"Scheme 4 and footnote","finding":"Nine 1,7-allenyne scope entries recorded independently with common syringe-pump conditions.","evidence_ids":[EV[6]]},
   {"region":"Scheme 5a","finding":"BHT and TEMPO radical-inhibitor controls recorded as two mechanistic-study experiments with 1.0 equiv additives and product yields.","evidence_ids":[EV[6]]},
   {"region":"Scheme 5b","finding":"1,8-Allenyne control gives bisallene 7 in 56%, recorded as a distinct mechanistic-study experiment.","evidence_ids":[EV[6]]},
   {"region":"individual paragraphs 1i, 4, and 6","finding":"Three concrete allenyne preparations with yield, product mass, and scale; spectra retained as context.","evidence_ids":[EV[6]]}],
  "record_ids":[rid("GP_A_generic"),rid("prep_1i"),rid("prep_4"),rid("prep_6")]+[rid(f"scheme4_5{x}") for x in "abcdefghi"]+[rid("scheme5_BHT"),rid("scheme5_TEMPO"),rid("scheme5_bisallene7")],
  "unresolved":["Scheme 4/5 yield measurement basis is unreported on this page.","Scheme 4 alkene geometry is not textually assigned in the inspected region.","Substituted allene 1i has no reported axial assignment; representation remains partial."]}
]}
dump(RUN / "coverage/pages_0001_0006.json", coverage)
tasks_path = RUN / "tasks.json"
if tasks_path.exists():
    tasks = json.loads(tasks_path.read_text())
    for task in tasks:
        page = int(task["evidence_id"].split(":page",1)[1].split(":",1)[0])
        if 1 <= page <= 6:
            task["state"] = "attempted"
            task["reason"] = "Original page image visually inspected in bounded pages 1-6 extraction; see coverage/pages_0001_0006.json."
    dump(tasks_path, tasks)
print(json.dumps({"molecules":len(molecules),"reactions":len(reactions),"stereo_diagnostics":len(stereo_diagnostics)},ensure_ascii=False))
