#!/usr/bin/env python3
"""Build the evidence-linked v3.1 bundle from the visually reviewed article pages."""
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DOC = "doc_5ee9ec3ae7e5b7c5"
EV = {p: [f"{DOC}:page{p}:image", f"{DOC}:page{p}:text"] for p in range(1, 7)}


def f(raw=None, value=None, evidence=(), method="original_extraction", derivation=None,
      status=None):
    return st.fact(raw, value, evidence, method=method, derivation=derivation, status=status)


def missing(evidence=()):
    return f(evidence=evidence)


def mid(label):
    return f"{DOC}:mol:{label}"


bundle = st.read(RUN / "submitted.json")
bundle["molecules"] = []
bundle["reactions"] = []
bundle["source_context"] = []
molecules = {}
stereo_diags = {}


def add_molecule(label, evidence, *, name=None, smiles=None, stereo_kind="unspecified",
                 geometry="unspecified", note=None, template=None):
    m = st.new_molecule(mid(label), DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    display = name or label
    m["names"]["original"] = f(display, display, evidence)
    m["names"]["english"] = f(display, display, evidence)
    m["cas_numbers"] = missing(evidence)
    m["identity_note"] = f(
        note or "Identity retained from the source label; no CAS number was reported.",
        note or "Identity retained from the source label; no CAS number was reported.",
        evidence,
    )
    if template is not None:
        m["representation_type"] = "generic_template"
        m["generic_template"] = template
        m["smiles"] = missing(evidence)
        rep = "not_encoded"
    elif smiles:
        diag = st.stereo_inspect(smiles)
        stereo_diags[label] = diag
        canonical = diag["canonical_isomeric_smiles"]
        m["representation_type"] = "smiles"
        m["smiles"] = f(
            {"source_label": label, "candidate_transcription": smiles}, canonical, evidence,
            method="structure_derivation",
            derivation="Candidate connectivity transcribed from the displayed structure; RDKit canonicalization checks graph syntax only, not source identity.",
        )
        m["stereochemistry"]["isomeric_smiles"] = f(
            smiles, canonical, evidence, method="structure_derivation",
            derivation="Canonical isomeric SMILES calculated from the candidate source transcription.",
        )
        m["stereochemistry"]["identity_key"] = f(
            None, diag["stereo_graph_key"], evidence, method="structure_derivation",
            derivation="stereo_graph_key calculated by standard stereo; it is not a sample identity.",
        )
        if any(p["specified"] == "Unspecified" for p in diag["potential_stereo"]):
            rep = "partial"
            stereo_kind = "unspecified"
            geometry = "unspecified"
        else:
            rep = "encoded"
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing(evidence)
        rep = "not_encoded"
    stereo_raw = {"source_label": label, "observation": "no wedge/dash or absolute/relative stereochemical descriptor was assigned in the displayed structure"}
    m["stereochemistry"]["chirality_status"] = f(stereo_raw, stereo_kind, evidence)
    m["stereochemistry"]["geometry_status"] = f(
        {"source_label": label, "observation": "double-bond geometry was not explicitly assigned unless noted"},
        geometry, evidence,
    )
    m["stereochemistry"]["representation_status"] = f(
        {"candidate_smiles_present": bool(smiles)}, rep, evidence,
        method="structure_derivation" if smiles else "inference",
        derivation="Compared the displayed stereochemical annotations with the stored candidate representation." if smiles else "No complete stereochemical structure was encoded from the dense source drawing.",
    )
    if smiles and stereo_kind == "achiral":
        m["stereochemistry"]["assignment"] = f(
            None, {"basis": "candidate graph has no assigned/potential tetrahedral or allene axis; terminal allene end is CH2 where applicable"}, evidence,
            method="structure_derivation",
            derivation="standard stereo graph diagnostics plus source scaffold inspection; axial support is limited and remains a review item.",
        )
    molecules[label] = m
    bundle["molecules"].append(m)
    return m


def generic_template(label, evidence, variables):
    return {
        "fixed_structure": f(
            {"figure": "Scheme 1", "label": label},
            {"description": "Generic allenyne/cyclized framework as drawn; atom-complete structure not encoded."},
            evidence,
        ),
        "variables": [
            {
                "label": v,
                "attachment_points": f({"variable": v, "location": "source drawing"}, {"location": "source drawing"}, evidence),
                "definition": f(v, v, evidence),
                "allowed_values": missing(evidence),
                "constraints": missing(evidence),
            }
            for v in variables
        ],
        "instance_mapping": [],
        "reason": f(
            "The source uses variable groups (R/Ar/E/n); a wildcard SMILES would falsely imply a concrete molecule.",
            "generic source structure with unresolved variables", evidence,
        ),
    }


# Candidate source structures. The allene terminus is CH2, so these candidates do not encode axial chirality.
tosyl = "S(=O)(=O)c1ccc(C)cc1"
def allenyne(alkyne_end, allene_sub="C"):
    return f"C=C=C({allene_sub})CN(CC#C{alkyne_end}){tosyl}"

aryl = {
    "1a": "c1ccccc1", "1b": "c1ccccc1C", "1c": "c1cccc(C)c1", "1d": "c1ccc(C)cc1",
    "1e": "c1ccc(C(C)(C)C)cc1", "1f": "c1ccc(OC)cc1", "1g": "c1ccc(F)cc1",
    "1h": "c1ccc(Br)cc1", "1i": "c1ccc(C(F)(F)F)cc1", "1j": "c1ccc(C#N)cc1",
    "1k": "c1ccc(C(=O)OC)cc1", "1l": "c1cccc2ccccc12", "1m": "c1ccsc1",
}
for label, end in aryl.items():
    add_molecule(label, EV[2] + EV[3], smiles=allenyne(end), stereo_kind="achiral", geometry="not_applicable",
                 note="Candidate allenyne connectivity transcribed from Scheme 3; source assigns no stereochemistry. Substituent position follows the displayed drawing and remains machine reviewable.")
add_molecule("1n", EV[2] + EV[3], note="Trisubstituted-allene substrate drawing retained by label; exact atom-complete transcription is unresolved.")
add_molecule("1o", EV[2] + EV[3], smiles=allenyne("C"), stereo_kind="achiral", geometry="not_applicable",
             note="Alkyl-substituted alkyne candidate transcribed from the no-reaction entry.")

for label in [f"2{x}" for x in "abcdefghijklmno"]:
    add_molecule(label, EV[2] + EV[3], note="Dense cyclized product drawing retained by source label; complete stereochemical SMILES was not guessed. The new all-carbon quaternary center is shown without wedge/dash assignment.")

for label, end in {
    "3a": "c1ccccc1", "3b": "c1cccc(CC)c1", "3c": "c1cccc(OC)c1",
    "3d": "c1cccc(F)c1", "3e": "c1cccc(C(=O)OC)c1", "3f": "C",
}.items():
    add_molecule(label, EV[3] + EV[4], smiles=allenyne(end, "c1ccccc1"), stereo_kind="achiral", geometry="not_applicable",
                 note="Candidate phenyl-substituted allenyne connectivity transcribed from Scheme 4; the terminal allene end is CH2 and no source stereochemistry is assigned.")
for label in [f"4{x}" for x in "abcdef"]:
    add_molecule(label, EV[3] + EV[4], note="3,5-Dihydroazepine product drawing retained by label; alkene geometry and any conformational/relative stereochemistry are not explicitly assigned.")
for label in ["5", "6", "7", "8"]:
    add_molecule(label, EV[3], note="Scheme 5 derivatization product retained by label; complete stereochemical SMILES was not guessed from the compact drawing.")

# Reagents/solvents and prior-work condition components.
simple = {
    "BCl3": ("boron trichloride", "B(Cl)(Cl)Cl"), "CH2Cl2": ("dichloromethane", "ClCCl"),
    "pinacol": ("pinacol", "CC(C)(O)C(C)(C)O"), "Et3N": ("triethylamine", "CCN(CC)CC"),
    "THF": ("tetrahydrofuran", "C1CCOC1"), "H2O": ("water", "O"),
    "PhI": ("iodobenzene", "Ic1ccccc1"), "Cs2CO3": ("cesium carbonate", "[Cs+].[Cs+].[O-]C([O-])=O"),
    "MeOH": ("methanol", "CO"), "p-I-toluene": ("1-iodo-4-methylbenzene", "Cc1ccc(I)cc1"),
    "Pd(OAc)2": ("palladium(II) acetate", "CC(=O)[O-].CC(=O)[O-].[Pd+2]"),
    "BQ": ("1,4-benzoquinone", "O=C1C=CC(=O)C=C1"), "DCE": ("1,2-dichloroethane", "ClCCCl"),
    "toluene": ("toluene", "Cc1ccccc1"), "B2Pin2": ("bis(pinacolato)diboron", "CC1(C)OB(B2OC(C)(C)C(C)(C)O2)OC1(C)C"),
}
for label, (name, smi) in simple.items():
    add_molecule(label, EV[2] + EV[3], name=name, smiles=smi, stereo_kind="achiral", geometry="not_applicable")
for label, name in {
    "NaBO3·4H2O": "sodium peroxoborate tetrahydrate", "Pd(PPh3)4": "tetrakis(triphenylphosphine)palladium(0)",
    "KHF2": "potassium bifluoride", "LiOAc·2H2O": "lithium acetate dihydrate", "BF3·Et2O": "boron trifluoride diethyl etherate",
}.items():
    add_molecule(label, EV[2] + EV[3], name=name, note="Named reagent is explicit, but a concrete salt/solvate SMILES was not required or guessed.")

# Generic prior/current-work structures and proposed intermediates.
for label, vars_ in {
    "prior-a-substrate": ["R", "E", "n"], "prior-a-product-top": ["R", "E"],
    "prior-a-product-middle": ["R", "E"], "prior-a-product-bottom": ["R", "E"],
    "prior-b-substrate": ["R", "E"], "prior-b-product-1": ["R", "E"], "prior-b-product-2": ["R", "E"],
    "current-generic-substrate-Me": ["Ar"], "current-generic-product-Me": ["Ar"],
    "current-generic-substrate-Ph": ["Ar"], "current-generic-product-Ph": ["Ar"],
}.items():
    add_molecule(label, EV[1] + EV[2], template=generic_template(label, EV[1] + EV[2], vars_),
                 note="Generic Scheme 1 identity; variables are retained rather than expanded into a false concrete SMILES.")
for label in list("ABCDEF"):
    add_molecule(f"intermediate-{label}", EV[4], note="Proposed mechanistic intermediate; identity is source-proposed, not an isolated experimental product.")


def participant(r, label, role, evidence, *, stage="reaction", equiv=None, molpct=None,
                reference=None, product_role=None):
    pid = f"p{len(r['participants'])+1}_{label.replace('·','_').replace('(','').replace(')','').replace(' ','_')}"
    p = st.new_participant(pid, mid(label), role, evidence, stage=stage)
    if product_role:
        p["product_role"] = f(product_role, product_role, evidence)
    if equiv is not None:
        p["charges"]["equivalents"] = f(f"{equiv} equiv", {"number": equiv, "unit": "equiv"}, evidence)
        p["equivalence_reference"] = f(reference, mid(reference), evidence)
    if molpct is not None:
        p["charges"]["mol_percent"] = f(f"{molpct} mol%", {"number": molpct, "unit": "mol%"}, evidence)
        p["equivalence_reference"] = f(reference, mid(reference), evidence)
    r["participants"].append(p)
    return pid


def set_common(r, evidence, location, record_type, categories, name, transformation):
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, evidence)
    r["location"] = f(location, location, evidence)
    r["reaction_name"] = f(name, name, evidence)
    r["transformation_type"] = f(
        transformation, transformation, evidence, method="inference",
        derivation="Transformation class inferred from the displayed reactant/product relationship and article terminology.",
    )
    r["review"]["unresolved"] = []


def add_step(r, step_id, evidence, pids, operation, temperature=None, time=None, stage="reaction"):
    s = st.new_step(step_id, evidence, stage=stage)
    s["participant_ids"] = list(pids)
    s["operation"] = f(operation, operation, evidence)
    if temperature is not None:
        s["temperature"] = f(temperature, temperature, evidence)
    if time is not None:
        s["time"] = f(time if isinstance(time, str) else f"{time['number']} {time['unit']}", time, evidence)
    r["steps"].append(s)


def set_classification(r, evidence, pairs, coverage="partial", limitations=None):
    changes = []
    labels = {}
    for i, (a, b, before, after, kind) in enumerate(pairs, 1):
        cid = f"b{i}"
        change = {"change_id": cid, "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
                  "before_order": before, "after_order": after, "change_type": kind}
        changes.append(change)
        if kind in ("formed", "hydrogen_gain", "order_increased"):
            lab = st.bond_label(a, b, after)
            labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(cid)
    r["classification"]["bond_changes"] = [f(None, c, evidence, method="structure_derivation",
        derivation="Partial net bond change inferred from the displayed structures; atom-complete mapping was not asserted.") for c in changes]
    r["classification"]["labels"] = f(None, list(labels.values()), evidence, method="structure_derivation",
        derivation="Labels calculated from the explicitly retained partial bond changes under bond-formation-20261016.")
    r["classification"]["coverage"] = f(None, coverage, evidence, method="inference",
        derivation="Coverage is partial because complete atom mapping and hydrogen/source balance are unavailable from the supplied main text.")
    r["classification"]["limitations"] = limitations or [
        "No atom-complete mapped reaction is asserted.",
        "Only source-evident net formed bonds are classified; unresolved bond-order changes and coproducts remain outside the mapping.",
    ]
    r["classification"]["model_votes"] = []


def no_classification(r, evidence):
    r["classification"]["bond_changes"] = []
    r["classification"]["labels"] = f([], [], evidence, method="not_applicable")
    r["classification"]["coverage"] = f("no desired reaction observed", "not_applicable", evidence, method="not_applicable")
    r["classification"]["limitations"] = ["No desired product was observed; no intended bond-formation label was invented."]
    r["classification"]["model_votes"] = []


def add_outcome(r, evidence, kind, target, raw, value, *, basis=None, scope=None):
    o = st.new_outcome(f"o{len(r['outcomes'])+1}", kind, [mid(target)] if target else [], evidence)
    o["value"] = f(raw, value, evidence)
    if target:
        o["product_role"] = f("desired product", "main_product", evidence)
    if kind == "yield":
        o["yield_type"] = f("reported yield", "reported", evidence)
    if basis:
        o["measurement_basis"] = f(basis, basis, evidence)
    if scope:
        o["scope"] = f(scope, scope, evidence)
    r["outcomes"].append(o)


def new_reaction(suffix, evidence):
    r = st.new_reaction(f"{DOC}:{suffix}", DOC, evidence)
    bundle["reactions"].append(r)
    return r


# Scheme 1: four depicted prior-work arrows.
prior_specs = [
    ("scheme1a_n0", "prior-a-substrate", ["prior-a-product-top"], "THF", None, None),
    ("scheme1a_n1_lioac", "prior-a-substrate", ["prior-a-product-middle"], "DCE", "LiOAc·2H2O", 20),
    ("scheme1a_n1_bf3", "prior-a-substrate", ["prior-a-product-bottom"], "THF", "BF3·Et2O", 20),
]
for suffix, sub, prods, solvent, additive, loading in prior_specs:
    r = new_reaction(suffix, EV[2])
    set_common(r, EV[2], {"page": 2, "figure": "Scheme 1a", "arrow": suffix}, "literature_example", [],
               "Pd-catalyzed oxidative borylative carbocyclization (cited work)", "carbocyclization/borylation")
    ps = participant(r, sub, "reactant", EV[2]); participant(r, "B2Pin2", "reagent", EV[2])
    pp = [participant(r, x, "product", EV[2], product_role="depicted product") for x in prods]
    pcat = participant(r, "Pd(OAc)2", "catalyst", EV[2], molpct=2, reference=sub)
    pbq = participant(r, "BQ", "reagent", EV[2], equiv=1.1, reference=sub)
    psol = participant(r, solvent, "solvent", EV[2])
    pids = [ps, pcat, pbq, psol] + pp
    if additive:
        pids.append(participant(r, additive, "additive", EV[2], molpct=loading, reference=sub))
    r["solvent_conditions"]["composition"] = f(solvent, [mid(solvent)], EV[2])
    add_step(r, "reaction", EV[2], pids, "single depicted reaction arrow", {"number": 50, "unit": "°C"})
    set_classification(r, EV[2], [("C", "C", 0, 1, "formed"), ("C", "B", 0, 1, "formed")])
    r["review"]["unresolved"] = ["Cited-work scheme is generic; time, substrate charge, yield, product stereochemistry and complete mapping are not reported here."]

r = new_reaction("scheme1b", EV[2])
set_common(r, EV[2], {"page": 2, "figure": "Scheme 1b"}, "literature_example", [],
           "Pd-catalyzed formal hydroborylative carbocyclization (cited work)", "carbocyclization/borylation")
ps = participant(r, "prior-b-substrate", "reactant", EV[2])
pb = participant(r, "B2Pin2", "reagent", EV[2], equiv=1.2, reference="prior-b-substrate")
pc = participant(r, "Pd(OAc)2", "catalyst", EV[2], molpct=5, reference="prior-b-substrate")
pm = participant(r, "MeOH", "reagent", EV[2], equiv=1.0, reference="prior-b-substrate")
pt = participant(r, "toluene", "solvent", EV[2])
pp1 = participant(r, "prior-b-product-1", "product", EV[2], product_role="depicted product")
pp2 = participant(r, "prior-b-product-2", "product", EV[2], product_role="depicted product")
r["solvent_conditions"]["composition"] = f("toluene", [mid("toluene")], EV[2])
add_step(r, "reaction", EV[2], [ps, pb, pc, pm, pt, pp1, pp2], "single depicted reaction arrow", {"number": 50, "unit": "°C"})
set_classification(r, EV[2], [("C", "C", 0, 1, "formed"), ("C", "B", 0, 1, "formed")])
r["review"]["unresolved"] = ["Cited-work scheme is generic and depicts two products without yields; time and complete mapping are unreported."]

# Scheme 1c/current graphical abstract: two generic branches, evidence-merged with page 1.
for tag, sub, prod in [("Me", "current-generic-substrate-Me", "current-generic-product-Me"),
                       ("Ph", "current-generic-substrate-Ph", "current-generic-product-Ph")]:
    evidence = EV[1] + EV[2]
    r = new_reaction(f"scheme1c_{tag.lower()}", evidence)
    set_common(r, evidence, {"pages": [1, 2], "figures": ["graphical abstract", "Scheme 1c"], "branch": f"R={tag}"},
               "generic_scheme", [], "Borylative cyclization of 1,6-allenynes driven by BCl3", "carbocyclization/borylation")
    ps = participant(r, sub, "reactant", evidence); pb = participant(r, "BCl3", "reagent", evidence)
    ppn = participant(r, "pinacol", "reagent", evidence); pp = participant(r, prod, "product", evidence, product_role="generic product")
    add_step(r, "generic_arrow", evidence, [ps, pb, ppn, pp], "BCl3 then pinacol")
    set_classification(r, evidence, [("C", "C", 0, 1, "formed"), ("C", "B", 0, 1, "formed"), ("C", "Cl", 0, 1, "formed")])
    r["review"]["unresolved"] = ["Generic scheme only; Ar identity, charges, time, workup, yield and atom-complete mapping are not specified."]

# Schemes 2/3: 15 scope experiments. Scheme 2 and Scheme 3 2a are one experiment.
scope3 = [
    ("a", 90, 1.5, "isolated"), ("b", 80, 1.5, None), ("c", 85, 1.5, None),
    ("d", 87, 1.5, None), ("e", 82, 1.5, None), ("f", 26, 1.5, None),
    ("g", 85, 1.5, None), ("h", 84, 1.5, None), ("i", 89, 1.5, None),
    ("j", 78, 1.5, "isolated"), ("k", 87, 2.0, None), ("l", 82, 1.5, None),
    ("m", 81, 1.5, None), ("n", 81, 1.5, None), ("o", None, 1.5, None),
]
for code, yield_value, bcl3_eq, basis in scope3:
    evidence = EV[2] + (EV[3] if code in "ijklmno" else [])
    r = new_reaction(f"scope_2{code}", evidence)
    set_common(r, evidence, {"page": 2, "figure": "Scheme 3", "entry": f"2{code}"}, "specific_experiment", ["substrate_scope"],
               "Borylative cyclization of 1,6-allenyne", "intramolecular carbocyclization/borylation")
    ps = participant(r, f"1{code}", "reactant", evidence)
    pb = participant(r, "BCl3", "reagent", evidence, equiv=bcl3_eq, reference=f"1{code}")
    pd = participant(r, "CH2Cl2", "solvent", evidence)
    pi = participant(r, "pinacol", "reagent", evidence, equiv=1.1, reference=f"1{code}")
    pe = participant(r, "Et3N", "reagent", evidence, equiv=15, reference=f"1{code}")
    product_pid = None
    if yield_value is not None:
        product_pid = participant(r, f"2{code}", "product", evidence, product_role="desired product")
    r["solvent_conditions"]["composition"] = f("CH2Cl2", [mid("CH2Cl2")], evidence)
    t1 = {"number": 10, "unit": "min"} if code == "a" else None
    add_step(r, "bcl3_cyclization", evidence, [ps, pb, pd], "treat allenyne with BCl3 in CH2Cl2", "r.t.", t1)
    add_step(r, "pinacol_trapping", evidence, [pi, pe] + ([product_pid] if product_pid else []), "treat with pinacol and Et3N", "r.t.")
    if yield_value is None:
        add_outcome(r, evidence, "no_reaction", None, "0%; no reaction / no product observed", "no desired product observed")
        no_classification(r, evidence)
        r["review"]["unresolved"] = ["No desired product was observed; detection limit, substrate recovery/conversion and decomposition products are unreported."]
    else:
        add_outcome(r, evidence, "yield", f"2{code}", f"{yield_value}%", {"number": yield_value, "unit": "%"}, basis=basis)
        set_classification(r, evidence, [("C", "C", 0, 1, "formed"), ("C", "B", 0, 1, "formed"), ("C", "Cl", 0, 1, "formed")])
        unresolved = ["Substrate charge, solvent volume, stage-2 duration, workup and purification are unavailable without the unprovided SI."]
        if not basis: unresolved.append("Yield measurement basis is not stated for this entry in the supplied main text.")
        if code == "n": unresolved.append("Complete substrate/product structures and stereochemical representation remain unresolved from the compact drawing.")
        r["review"]["unresolved"] = unresolved
    if code == "a":
        r["location"] = f({"page": 2, "figures": ["Scheme 2", "Scheme 3"], "entry": "2a"},
                          {"page": 2, "figures": ["Scheme 2", "Scheme 3"], "entry": "2a"}, evidence)
        r["procedure"]["original_text"] = f(
            "1a and BCl3 in CH2Cl2 at ambient temperature for 10 min; treatment with pinacol and Et3N gave 2a in 90% isolated yield.",
            "1a and BCl3 in CH2Cl2 at ambient temperature for 10 min; then pinacol and Et3N", evidence)

# Scheme 4: six phenyl-allene scope entries.
scope4 = [("a", 87, 1.5), ("b", 83, 1.5), ("c", 79, 1.5), ("d", 76, 1.5), ("e", 79, 2.0), ("f", None, 1.5)]
for code, yield_value, bcl3_eq in scope4:
    evidence = EV[3]
    r = new_reaction(f"scope_4{code}", evidence)
    set_common(r, evidence, {"page": 3, "figure": "Scheme 4", "entry": f"4{code}"}, "specific_experiment", ["substrate_scope"],
               "Borylative cyclization of phenyl-substituted 1,6-allenyne", "intramolecular carbocyclization/borylation")
    ps = participant(r, f"3{code}", "reactant", evidence)
    pb = participant(r, "BCl3", "reagent", evidence, equiv=bcl3_eq, reference=f"3{code}")
    pd = participant(r, "CH2Cl2", "solvent", evidence)
    pi = participant(r, "pinacol", "reagent", evidence, equiv=1.1, reference=f"3{code}")
    pe = participant(r, "Et3N", "reagent", evidence, equiv=15, reference=f"3{code}")
    product_pid = participant(r, f"4{code}", "product", evidence, product_role="desired product") if yield_value is not None else None
    r["solvent_conditions"]["composition"] = f("CH2Cl2", [mid("CH2Cl2")], evidence)
    add_step(r, "bcl3_cyclization", evidence, [ps, pb, pd], "treat phenyl-substituted allenyne with BCl3 in CH2Cl2", "r.t.")
    add_step(r, "pinacol_trapping", evidence, [pi, pe] + ([product_pid] if product_pid else []), "treat with pinacol and Et3N", "r.t.")
    if yield_value is None:
        add_outcome(r, evidence, "no_reaction", None, "0%; not detected; substrate disappeared", "desired product not detected")
        no_classification(r, evidence)
        r["review"]["unresolved"] = ["The substrate disappeared, but conversion, decomposition products, detection method and detection limit are unreported."]
    else:
        add_outcome(r, evidence, "yield", f"4{code}", f"{yield_value}%", {"number": yield_value, "unit": "%"})
        set_classification(r, evidence, [("C", "C", 0, 1, "formed"), ("C", "B", 0, 1, "formed"), ("C", "Cl", 0, 1, "formed")])
        r["review"]["unresolved"] = ["Reaction time, scale, yield basis, workup, purification and complete product stereochemical structure are unreported in the supplied main text."]

# Scheme 5 derivatizations.
derivs = [
    ("2a_to_5", "2a", "5", ["NaBO3·4H2O", "THF", "H2O"], 89, [("C", "O", 1, 2, "order_increased")], "oxidation to ketone"),
    ("2a_to_6", "2a", "6", ["PhI", "Pd(PPh3)4", "Cs2CO3", "THF"], 82, [("C", "C", 0, 1, "formed")], "Suzuki coupling"),
    ("2a_to_7", "2a", "7", ["KHF2", "MeOH"], 81, [("B", "F", 0, 1, "formed"), ("B", "F", 0, 1, "formed"), ("B", "F", 0, 1, "formed")], "conversion to organotrifluoroborate"),
    ("4a_to_8", "4a", "8", ["p-I-toluene", "Pd(PPh3)4", "Cs2CO3", "THF"], 83, [("C", "C", 0, 1, "formed")], "Suzuki coupling"),
]
for suffix, sub, prod, agents, yld, changes, trans in derivs:
    evidence = EV[3]
    r = new_reaction(suffix, evidence)
    set_common(r, evidence, {"page": 3, "figure": "Scheme 5", "arrow": suffix}, "specific_experiment", ["derivatization"],
               "Derivatization of borylative cyclization product", trans)
    ps = participant(r, sub, "reactant", evidence)
    pids = [ps]
    for a in agents:
        role = "solvent" if a in {"THF", "H2O", "MeOH"} else ("catalyst" if a == "Pd(PPh3)4" else "reagent")
        pids.append(participant(r, a, role, evidence))
    pp = participant(r, prod, "product", evidence, product_role="desired product"); pids.append(pp)
    solvents = [a for a in agents if a in {"THF", "H2O", "MeOH"}]
    if solvents: r["solvent_conditions"]["composition"] = f(", ".join(solvents), [mid(a) for a in solvents], evidence)
    add_step(r, "reaction", evidence, pids, "single Scheme 5 reaction arrow")
    add_outcome(r, evidence, "yield", prod, f"{yld}%", {"number": yld, "unit": "%"})
    set_classification(r, evidence, changes)
    r["review"]["unresolved"] = ["Reagent doses, scale, temperature, time, workup, purification and yield measurement basis are not reported in the supplied main text.",
                                    "Complete reactant/product structures and atom mapping remain unresolved from the compact drawings."]

# Scheme 6 proposed mechanisms: one record per substrate pathway, not experiments.
for tag, sub, prod, ints in [("1a", "1a", "2a", ["intermediate-A", "intermediate-B", "intermediate-C"]),
                             ("3a", "3a", "4a", ["intermediate-D", "intermediate-E", "intermediate-F"])]:
    evidence = EV[3] + EV[4]
    r = new_reaction(f"mechanism_{tag}", evidence)
    set_common(r, evidence, {"page": 4, "figure": "Scheme 6", "pathway": tag}, "proposed_mechanism", [],
               "Proposed mechanism", "proposed zwitterionic carbocyclization/borylation pathway")
    pids = [participant(r, sub, "reactant", evidence), participant(r, "BCl3", "reagent", evidence)]
    pids += [participant(r, x, "unknown", evidence) for x in ints]
    pids += [participant(r, "pinacol", "reagent", evidence), participant(r, "Et3N", "reagent", evidence),
             participant(r, prod, "product", evidence, product_role="proposed net product")]
    add_step(r, "proposed_sequence", evidence, pids, "author-proposed sequence of Lewis-acid activation, cyclization and pinacol trapping")
    set_classification(r, evidence, [("C", "C", 0, 1, "formed"), ("C", "B", 0, 1, "formed"), ("C", "Cl", 0, 1, "formed")],
                       limitations=["Author-proposed mechanism, not a distinct performed experiment.", "Intermediates A–F are proposed and not reported as isolated.", "No atom-complete mapping is asserted."])
    r["review"]["unresolved"] = ["Elementary-step atom mapping, charges and stereochemical course were not independently verified; this record preserves the proposed mechanism only."]

# Source context, including non-reaction pages and explicitly unprovided SI.
bundle["source_context"] = [
    {"context_id": f"{DOC}:context:abstract", "document_id": DOC, "kind": "article_summary", "evidence_ids": EV[1],
     "text": "Abstract and graphical abstract describe the metal-free BCl3-driven borylative cyclization; the graphic duplicates the two current-work generic branches."},
    {"context_id": f"{DOC}:context:dft", "document_id": DOC, "kind": "mechanistic_discussion", "evidence_ids": EV[3] + EV[4],
     "text": "DFT discussion attributes five- versus seven-membered cyclization to conjugation, steric hindrance and carbocation stability; this is not a separate synthetic experiment."},
    {"context_id": f"{DOC}:context:si", "document_id": DOC, "kind": "supplement_notice", "evidence_ids": EV[4],
     "text": "The article advertises SI with procedures, characterization, NMR and X-ray data, but no SI file was explicitly supplied or associated in this session."},
    {"context_id": f"{DOC}:context:ccdc", "document_id": DOC, "kind": "crystallography", "evidence_ids": EV[2] + EV[3] + EV[4],
     "text": "CCDC 2095006 (2d) and 2095011 (4a) are cited; crystallographic depictions are supporting identity evidence, not additional syntheses."},
    {"context_id": f"{DOC}:context:references", "document_id": DOC, "kind": "bibliography", "evidence_ids": EV[5] + EV[6],
     "text": "Pages 5–6 contain references only; bibliography entries were not promoted to current-paper experiments."},
]

st.write(RUN / "submitted.json", bundle)

# Raw visual transcription and decision table retained separately from the formal bundle.
raw = {
    "document_id": DOC,
    "title": "Borylative Cyclization of 1,6-Allenynes Driven by BCl3",
    "visual_pages_reviewed": [1, 2, 3, 4, 5, 6],
    "specific_experiments": {
        "scheme_3": [{"substrate": f"1{x}", "product": f"2{x}", "yield_percent": y, "BCl3_equiv": eq, "basis": basis}
                     for x, y, eq, basis in scope3],
        "scheme_4": [{"substrate": f"3{x}", "product": f"4{x}" if y is not None else None, "yield_percent": y, "BCl3_equiv": eq}
                     for x, y, eq in scope4],
        "scheme_5": [{"reaction": x[0], "substrate": x[1], "product": x[2], "yield_percent": x[4]} for x in derivs],
    },
    "deduplication": ["Scheme 2 1a→2a and Scheme 3 entry 2a are the same experiment and share one record.",
                       "Page-1 graphical abstract and Scheme 1c are duplicate presentations of the same two generic branches."],
    "non_specific_records": {"prior_work_arrows": 4, "current_generic_branches": 2, "proposed_mechanism_pathways": 2},
    "limitations": ["No SI was supplied.", "Dense product drawings were not converted into guessed complete stereochemical SMILES.",
                    "Machine records remain needs_review; visual inspection is documented separately and does not certify chemical accuracy."],
}
st.write(RUN / "raw_model/extracted_rows.json", raw)
st.write(RUN / "raw_model/stereo_diagnostics.json", {
    "document_id": DOC,
    "source_verified": False,
    "diagnostics": stereo_diags,
    "decision": "Any candidate with toolkit-reported unspecified potential stereo is stored as partial/unspecified; no racemate or axial assignment is invented."
})

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title/abstract", "finding": "Article identity, abstract and current-work summary; no specific experimental row.", "evidence_ids": EV[1]},
            {"region": "graphical abstract", "finding": "Two generic current-work branches (R=Me and R=Ph), duplicate of Scheme 1c and linked to the same generic records.", "evidence_ids": EV[1]},
            {"region": "two-column introduction", "finding": "Background and cited work discussion; no additional current-paper experiment.", "evidence_ids": EV[1]},
        ], "record_ids": [f"{DOC}:scheme1c_me", f"{DOC}:scheme1c_ph"],
         "unresolved": ["Graphical abstract omits charges, time and yield; details are not inferred."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Scheme 1a", "finding": "Three distinct cited-work arrows with n/additive/solvent differences; stored as literature examples.", "evidence_ids": EV[2]},
            {"region": "Scheme 1b", "finding": "One cited-work experiment depicting two products; stored as one literature-example record.", "evidence_ids": EV[2]},
            {"region": "Scheme 1c", "finding": "Two current-work generic branches, evidence-merged with page-1 graphic.", "evidence_ids": EV[2]},
            {"region": "Scheme 2 and adjoining prose", "finding": "Model experiment 1a→2a, 10 min first stage, 90% isolated yield; merged with Scheme 3 entry 2a.", "evidence_ids": EV[2]},
            {"region": "Scheme 3", "finding": "Entries 2a–2o: 14 positive-yield scope experiments and one no-reaction entry; common reagent equivalents and 2k exception transcribed.", "evidence_ids": EV[2]},
        ], "record_ids": [f"{DOC}:scheme1a_n0", f"{DOC}:scheme1a_n1_lioac", f"{DOC}:scheme1a_n1_bf3", f"{DOC}:scheme1b",
                          f"{DOC}:scheme1c_me", f"{DOC}:scheme1c_ph"] + [f"{DOC}:scope_2{x}" for x in "abcdefghijklmno"],
         "unresolved": ["Most Scheme 3 yield measurement bases and reaction times are not stated.", "Exact structure/stereochemistry of compact 1n/2n drawing remains unresolved."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "continued Scheme 3 discussion", "finding": "Textual yields/functional-group discussion and no-product statement associated with existing Scheme 3 records.", "evidence_ids": EV[3]},
            {"region": "Scheme 4", "finding": "Entries 4a–4f: five positive-yield scope experiments and one substrate-disappearance/no-desired-product result; common equivalents and 4e exception transcribed.", "evidence_ids": EV[3]},
            {"region": "Scheme 5", "finding": "Four independent derivatization arrows 2a→5, 2a→6, 2a→7 and 4a→8 with printed reagents/solvents/yields.", "evidence_ids": EV[3]},
            {"region": "DFT/mechanistic prose", "finding": "Author interpretation and pointer to proposed Scheme 6; not a separate performed experiment.", "evidence_ids": EV[3]},
        ], "record_ids": [f"{DOC}:scope_2{x}" for x in "ijklmno"] + [f"{DOC}:scope_4{x}" for x in "abcdef"] +
                         [f"{DOC}:{x}" for x in ["2a_to_5", "2a_to_6", "2a_to_7", "4a_to_8", "mechanism_1a", "mechanism_3a"]],
         "unresolved": ["Scheme 4 times and most yield bases are unreported.", "Scheme 5 doses, temperatures, times, workups, purification and yield bases are unreported.",
                        "Dense product stereochemical structures remain partial/unencoded."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "Scheme 6 upper pathway", "finding": "Proposed 1a→A→B→C→2a mechanism; recorded as proposed mechanism, not an experiment.", "evidence_ids": EV[4]},
            {"region": "Scheme 6 lower pathway", "finding": "Proposed 3a→D→E→F→4a mechanism; recorded separately as proposed mechanism.", "evidence_ids": EV[4]},
            {"region": "summary and Associated Content", "finding": "Article summary, SI availability notice and CCDC codes; no new synthetic experiment.", "evidence_ids": EV[4]},
            {"region": "author/acknowledgment column", "finding": "Administrative content and beginning of references; no reaction record.", "evidence_ids": EV[4]},
        ], "record_ids": [f"{DOC}:mechanism_1a", f"{DOC}:mechanism_3a"],
         "unresolved": ["Proposed intermediate structures are not treated as isolated compounds or experimental observations.", "SI is advertised but not supplied/associated."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "references 1–9", "finding": "Bibliography only; cited papers were not reconstructed as experiments from citations.", "evidence_ids": EV[5]},
        ], "record_ids": [], "unresolved": ["No reaction arrows, experimental tables or independent characterization peak tables occur on this page."]},
        {"page": 6, "state": "attempted", "regions": [
            {"region": "references 9–10 continuation", "finding": "Bibliography continuation only; no current-paper experiment.", "evidence_ids": EV[6]},
            {"region": "remaining blank page area", "finding": "No additional content or reaction region.", "evidence_ids": [EV[6][0]]},
        ], "record_ids": [], "unresolved": ["No reaction arrows, experimental tables or independent characterization peak tables occur on this page."]},
    ],
}
(RUN / "coverage").mkdir(parents=True, exist_ok=True)
st.write(RUN / "coverage/pages_0001_0006.json", coverage)

# Mark only the six bounded page tasks attempted.
tasks = st.read(RUN / "tasks.json")
for task in tasks:
    if task.get("document_id") == DOC and any(task.get("evidence_id") == f"{DOC}:page{p}:text" for p in range(1, 7)):
        task["state"] = "attempted"
        task["reason"] = "Original rendered page visually inspected; findings and unresolved items are mapped in coverage/pages_0001_0006.json."
st.write(RUN / "tasks.json", tasks)

print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
                  "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in bundle["reactions"])}, ensure_ascii=False))
