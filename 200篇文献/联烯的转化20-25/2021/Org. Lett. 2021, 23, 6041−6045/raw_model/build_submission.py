#!/usr/bin/env python3
"""Build the machine v3.1 submission from the five visually inspected main-text pages."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DATASET = json.loads((RUN / "dataset.json").read_text())
DOC = "doc_bdc226a45f665dd1"
P = {n: [f"{DOC}:page{n}:image", f"{DOC}:page{n}:text"] for n in range(1, 6)}


def F(raw=None, norm=None, ev=None, method="original_extraction", derivation=None,
      status=None, verification=None, alternatives=()):
    if norm is None and raw is not None:
        norm = raw
    return st.fact(raw, norm, ev or [], method=method, derivation=derivation,
                   status=status, verification=verification, alternatives=alternatives)


def inferred(value, ev, why):
    return F(None, value, ev, method="inference", derivation=why)


bundle = st.empty(DATASET["documents"])
molecules = {}


def molecule(mid, label, ev, *, name=None, smiles=None, chirality="unresolved",
             geometry="unresolved", representation=None, stereo_note=None):
    if mid in molecules:
        return mid
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, ev)
    m["names"]["original"] = F(name or label, name or label, ev)
    m["names"]["english"] = F(name or label, name or label, ev)
    m["names"]["abbreviations"] = F()
    m["names"]["aliases"] = F()
    m["cas_numbers"] = F()
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = F(smiles, smiles, ev, method="structure_derivation",
                         derivation="Candidate molecular graph transcribed from the displayed structure or explicit chemical name; parsing is not source verification.")
        d = st.stereo_inspect(smiles)
        rep = representation or ("partial" if any(x["specified"] == "Unspecified" for x in d["potential_stereo"]) else "encoded")
        m["stereochemistry"]["isomeric_smiles"] = F(None, d["canonical_isomeric_smiles"], ev,
            method="structure_derivation", derivation="Computed with standard stereo diagnostics from the candidate SMILES.")
        m["stereochemistry"]["identity_key"] = F(None, d["stereo_graph_key"], ev,
            method="structure_derivation", derivation="Stereo graph key computed from candidate SMILES; it is not a sample identity.")
        m["stereochemistry"]["representation_status"] = inferred(rep, ev,
            "Compared the candidate graph with potential stereo elements reported by the standard stereo diagnostic; axial/allene support remains limited.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = F()
        m["stereochemistry"]["representation_status"] = inferred(representation or "not_encoded", ev,
            "The source drawing/label was retained, but an exact atom- and stereo-complete graph was not safely encoded from the main text.")
    m["generic_template"] = None
    m["stereochemistry"]["chirality_status"] = inferred(chirality, ev,
        "Status assigned from the displayed graph and labels; absence of a printed descriptor is not treated as racemic.")
    m["stereochemistry"]["geometry_status"] = inferred(geometry, ev,
        "Status assigned from the displayed alkene/allene geometry and explicit ratio annotations.")
    if stereo_note:
        m["stereochemistry"]["assignment"] = F(stereo_note, stereo_note, ev)
    m["identity_note"] = F(
        "Machine identity from the main article only; SI was not provided and exact dense structures remain candidates." if smiles else
        "Main-text label/drawing retained; exact graph and complete stereochemical encoding require expert review or the unprovided SI.",
        None if False else ("Machine identity from the main article only; SI was not provided and exact dense structures remain candidates." if smiles else
        "Main-text label/drawing retained; exact graph and complete stereochemical encoding require expert review or the unprovided SI."), ev)
    m["external_properties"] = []
    molecules[mid] = m
    return mid


# Concrete identities whose graphs are sufficiently simple to retain as candidates.
molecule("mol_1a_prime", "1a′", P[2], name="α-methyl phenylallene 1a′",
         smiles="C=C=C(C)c1ccccc1", chirality="achiral", geometry="not_applicable")
molecule("mol_1a_dprime", "1a″", P[2], name="allenylmethyl carbonate 1a″",
         smiles="COC(=O)OCC(=C=C)c1ccccc1", chirality="achiral", geometry="not_applicable")
molecule("mol_1a", "1a", P[1] + P[2] + P[3] + P[4], name="β-allenyl trimethylsilane 1a",
         smiles="C=C=C(C[Si](C)(C)C)c1ccccc1", chirality="achiral", geometry="not_applicable")
molecule("mol_2a", "2a", P[2] + P[3] + P[4], name="diethyl methylmalonate",
         smiles="CCOC(=O)C(C)C(=O)OCC", chirality="achiral", geometry="not_applicable")

simple = {
    "cubr2": ("CuBr2", "copper(II) bromide", "[Cu+2].[Br-].[Br-]"),
    "cucl2": ("CuCl2", "copper(II) chloride", "[Cu+2].[Cl-].[Cl-]"),
    "cuoac2": ("Cu(OAc)2", "copper(II) acetate", "[Cu+2].CC(=O)[O-].CC(=O)[O-]"),
    "cubr": ("CuBr", "copper(I) bromide", "[Cu+].[Br-]"),
    "cui": ("CuI", "copper(I) iodide", "[Cu+].[I-]"),
    "cucn": ("CuCN", "copper(I) cyanide", "[Cu+].[C-]#N"),
    "dtbp": ("DTBP", "di-tert-butyl peroxide", "CC(C)(C)OOC(C)(C)C"),
    "dmso": ("DMSO", "dimethyl sulfoxide", "CS(C)=O"),
    "pcy3": ("PCy3", "tricyclohexylphosphine", "C1CCC(CC1)P(C1CCCCC1)C1CCCCC1"),
    "bpy": ("2,2′-bipyridyl", "2,2′-bipyridine", "n1ccccc1-c1ccccn1"),
    "cy2nh": ("Cy2NH", "dicyclohexylamine", "C1CCC(CC1)NC1CCCCC1"),
    "n2": ("N2", "nitrogen", "N#N"),
    "toluene": ("toluene", "toluene", "Cc1ccccc1"),
    "npmi": ("N-phenylmaleimide", "N-phenylmaleimide", "O=C1C=CC(=O)N1c1ccccc1"),
    "maleic": ("maleic anhydride", "maleic anhydride", "O=C1OC(=O)C=C1"),
    "tcne": ("tetracyanoethylene", "tetracyanoethylene", "N#CC(=C(C#N)C#N)C#N"),
    "lialh4": ("LiAlH4", "lithium aluminium hydride", "[Li+].[AlH4-]"),
    "et2o": ("Et2O", "diethyl ether", "CCOCC"),
    "mecn": ("CH3CN", "acetonitrile", "CC#N"),
    "tempo": ("TEMPO", "2,2,6,6-tetramethylpiperidine N-oxide", "CC1(C)CCCC(C)(C)N1[O]"),
    "bht": ("BHT", "2,6-di-tert-butyl-4-methylphenol", "Cc1cc(C(C)(C)C)c(O)c(C(C)(C)C)c1"),
}
for key, (lab, name, smi) in simple.items():
    pages = P[2]
    if key in {"toluene", "npmi", "maleic", "tcne", "lialh4", "et2o", "mecn"}: pages = P[3]
    if key in {"tempo", "bht"}: pages = P[3] + P[4]
    molecule("mol_" + key, lab, pages, name=name, smiles=smi,
             chirality="unresolved" if key == "lialh4" else "achiral", geometry="not_applicable")
molecule("mol_bf3et2o", "BF3·Et2O", P[3], name="boron trifluoride diethyl etherate",
         chirality="not_applicable", geometry="not_applicable")

# Source-labelled substrates/products remain separate identities. Dense exact graphs are not guessed.
molecule("mol_3aa", "3aa", P[1] + P[2] + P[3] + P[4], name="branched 1,3-diene 3aa",
         chirality="unresolved", geometry="specified", representation="not_encoded",
         stereo_note="A single diene geometry is drawn, but no E/Z descriptor is printed in the main text.")
for letter in "bcdefghijklmnopqr":
    molecule(f"mol_1{letter}", f"1{letter}", P[3], name=f"β-allenyl silane 1{letter}",
             chirality="unresolved", geometry="not_applicable", representation="not_encoded")
    g = "mixture" if letter in "opq" else ("unresolved" if letter == "r" else "specified")
    c = "unresolved"
    note = "E/Z mixture explicitly reported for the product." if letter in "opq" else "Displayed product geometry retained without assigning an unprinted E/Z descriptor."
    molecule(f"mol_3{letter}a", f"3{letter}a", P[3], name=f"branched 1,3-diene 3{letter}a",
             chirality=c, geometry=g, representation="not_encoded", stereo_note=note)
for letter in "bcdefghijklmnopq":
    molecule(f"mol_2{letter}", f"2{letter}", P[3], name=f"tertiary C–H substrate 2{letter}",
             chirality="unresolved", geometry="not_applicable", representation="not_encoded")
    molecule(f"mol_3a{letter}", f"3a{letter}", P[3], name=f"branched 1,3-diene 3a{letter}",
             chirality="unresolved", geometry="specified", representation="not_encoded",
             stereo_note="Displayed product geometry retained without assigning an unprinted E/Z descriptor.")
for n in range(4, 9):
    molecule(f"mol_{n}", str(n), P[3], name=f"derivatization product {n}",
             chirality="unresolved", geometry="unresolved", representation="not_encoded")
molecule("mol_9", "9", P[4], name="TEMPO radical-trapping adduct 9",
         chirality="unresolved", geometry="not_applicable", representation="not_encoded")
for lab in ("A", "B", "C"):
    molecule(f"mol_int_{lab}", lab, P[4], name=f"proposed intermediate {lab}",
             chirality="unresolved", geometry="unresolved", representation="not_encoded")
molecule("mol_tmsotbu", "TMSOtBu", P[4], name="tert-butoxytrimethylsilane",
         smiles="CC(C)(C)O[Si](C)(C)C", chirality="achiral", geometry="not_applicable")


def participant(pid, mid, role, ev, *, amount=None, equiv=None, molpct=None, volume=None,
                ref=None, stage="reaction", product_role=None, method="original_extraction"):
    p = st.new_participant(pid, mid, role, ev, stage=stage)
    if product_role is not None:
        p["product_role"] = F(product_role, product_role, ev)
    values = {"amount": amount, "equivalents": equiv, "mol_percent": molpct, "volume": volume}
    for key, val in values.items():
        if val is not None:
            raw = val.get("raw_text") if isinstance(val, dict) and val.get("raw_text") else val
            p["charges"][key] = F(raw, val, ev, method=method)
    if any(v is not None for v in (equiv, molpct)):
        p["equivalence_reference"] = F(ref, ref, ev, method=method)
    return p


def step(sid, ev, pids, operation, temp=None, time=None, atmosphere=None,
         *, temp_method="original_extraction", time_method="original_extraction"):
    s = st.new_step(sid, ev)
    s["participant_ids"] = list(pids)
    s["operation"] = F(operation, operation, ev)
    if temp is not None:
        s["temperature"] = F(temp.get("raw_text", temp) if isinstance(temp, dict) else temp, temp, ev, method=temp_method)
    if time is not None:
        s["time"] = F(time.get("raw_text", time) if isinstance(time, dict) else time, time, ev, method=time_method)
    if atmosphere is not None:
        s["atmosphere"] = F(atmosphere, atmosphere, ev)
    return s


def classification(ev, labels=(), *, not_applicable=False, limitations=None):
    if not_applicable:
        return {"standard_id": st.STANDARD, "mapped_reaction_smiles": F(), "bond_changes": [],
                "labels": inferred([], ev, "No observed transformation; no bond-formation label assigned."),
                "coverage": inferred("not_applicable", ev, "Negative/no-reaction control."),
                "limitations": limitations or ["No observed transformation was classified."], "model_votes": []}
    changes = []
    outlabels = []
    for i, code in enumerate(labels, 1):
        order = 2 if code.startswith("DB-") else 3 if code.startswith("TB-") else 1
        pair = code.removeprefix("DB-").removeprefix("TB-").split("-")
        if pair == ["O", "H"]: atoms = [{"atom_map": None, "element": "O"}, {"atom_map": None, "element": "H"}]
        elif pair == ["C", "H"]: atoms = [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}]
        else: atoms = [{"atom_map": None, "element": pair[0]}, {"atom_map": None, "element": pair[1]}]
        change = {"change_id": f"b{i}", "atoms": atoms, "before_order": 0 if order == 1 else 1,
                  "after_order": order, "change_type": "formed" if order == 1 else "order_increased"}
        changes.append(F(None, change, ev, method="structure_derivation",
                         derivation="Category-level change inferred from the displayed reactant/product drawings; no atom-complete mapping asserted."))
        row = st.bond_label(pair[0], pair[1], order)
        row["supporting_change_ids"] = [f"b{i}"]
        outlabels.append(row)
    return {"standard_id": st.STANDARD, "mapped_reaction_smiles": F(), "bond_changes": changes,
            "labels": F(None, outlabels, ev, method="structure_derivation",
                        derivation="All asserted labels correspond to displayed new connections or promoted bonds; unchanged product bonds are excluded."),
            "coverage": inferred("partial", ev, "Category-level source-drawing comparison; no atom-complete mapping."),
            "limitations": limitations or ["Dense source drawings were not atom-mapped.", "Hydrogen losses and coproduct atom balances are incomplete."],
            "model_votes": []}


def base_reaction(suffix, ev, category, location, name, transformation, labels=("C-C", "DB-C-C"),
                  record_type="specific_experiment", no_change=False):
    r = st.new_reaction(f"{DOC}:{suffix}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = F(category, category, ev)
    r["location"] = F(location, location, ev)
    r["reaction_name"] = F(name, name, ev)
    r["transformation_type"] = inferred(transformation, ev, "Transformation type assigned from the displayed net structures and authors' description.")
    r["classification"] = classification(ev, labels, not_applicable=no_change)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
                   "unresolved": ["Machine extraction from the main article; SI procedures and characterization were not provided."]}
    r["completeness"] = "partial"
    return r


def set_solvent(r, mid, ev, raw=None, ratio=None, volume=None):
    r["solvent_conditions"]["composition"] = F(raw or mid, [mid], ev)
    if ratio is not None: r["solvent_conditions"]["mixture_ratio"] = F(ratio, ratio, ev)
    if volume is not None: r["solvent_conditions"]["total_volume"] = F(volume.get("raw_text", volume), volume, ev)


def outcome(oid, kind, targets, ev, value, *, yield_type=None, basis=None, scope=None):
    o = st.new_outcome(oid, kind, targets, ev)
    o["product_role"] = F("main product" if kind in {"yield", "product_amount"} else kind,
                           "main product" if kind in {"yield", "product_amount"} else kind, ev)
    raw = value.get("raw_text") if isinstance(value, dict) and value.get("raw_text") else value
    o["value"] = F(raw, value, ev)
    if yield_type is not None: o["yield_type"] = F(yield_type, yield_type, ev)
    if basis is not None: o["measurement_basis"] = F(basis, basis, ev)
    if scope is not None: o["scope"] = F(scope, scope, ev)
    return o


reactions = []

# Scheme 1c current-work generic transformation; 1a/1b are contextual prior art only.
r = base_reaction("scheme1c_generic", P[1], ["other"], "PDF p1, Scheme 1c and graphical abstract",
                  "Cu-catalyzed dienylation of tertiary C(sp3)–H bonds with β-allenyl silanes",
                  "oxidative cross-nucleophile coupling to branched 1,3-dienes", record_type="generic_scheme")
r["participants"] = [participant("gs_allene", "mol_1a", "reactant", P[1]), participant("gs_ch", "mol_2a", "reactant", P[1]),
                     participant("gs_product", "mol_3aa", "product", P[1], product_role="representative product")]
r["notes"] = [F("Scheme 1a and 1b summarize cited prior work/selectivity concerns and are not counted as current specific experiments.",
                  "Scheme 1a and 1b summarize cited prior work/selectivity concerns and are not counted as current specific experiments.", P[1])]
reactions.append(r)
GENERAL = r["reaction_id"]

# Scheme 2 discovery experiments. The 1a/44% run is identical to Table 1 entry 1 and is represented there once.
for suffix, sub, product, val, conversion in [
    ("discovery_1a_prime", "mol_1a_prime", "mol_3aa", 20, None),
    ("discovery_1a_dprime", "mol_1a_dprime", None, None, 0),
]:
    r = base_reaction(suffix, P[2], ["optimization"], "PDF p2, Scheme 2",
        "reaction discovery with α-substituted allene", "allene/malonate oxidative coupling")
    r["general_scheme_ids"] = [GENERAL]
    r["participants"] = [participant("allene", sub, "reactant", P[2], amount={"number": 0.15, "unit": "mmol"}),
                         participant("malonate", "mol_2a", "reactant", P[2], equiv={"number": 1.5, "unit": "equiv"}, ref=sub),
                         participant("cu", "mol_cubr2", "catalyst", P[2], molpct={"number": 20, "unit": "mol%"}, ref=sub),
                         participant("oxidant", "mol_dtbp", "reagent", P[2], equiv={"number": 2, "unit": "equiv"}, ref=sub),
                         participant("solvent", "mol_dmso", "solvent", P[2])]
    if product: r["participants"].append(participant("product", product, "product", P[2], product_role="coupling product"))
    r["steps"] = [step("reaction", P[2], [p["participant_id"] for p in r["participants"]], "heat under the displayed conditions",
                       {"number": 70, "unit": "°C"}, {"number": 12, "unit": "h"})]
    set_solvent(r, "mol_dmso", P[2])
    if val is not None:
        r["outcomes"] = [outcome("yield", "yield", [product], P[2], {"number": val, "unit": "%"},
                                  yield_type="reported yield", basis="isolated yield")]
    else:
        r["outcomes"] = [outcome("conversion", "conversion", [sub], P[2], {"number": conversion, "unit": "%"},
                                  basis="not reported", scope="substrate conversion")]
        r["classification"] = classification(P[2], not_applicable=True)
    reactions.append(r)

# Table 1: every row is an independent optimization experiment.
table = [
    (1,"cubr2",20,70,44,None),(2,"cucl2",20,70,18,None),(3,"cuoac2",20,70,"trace",None),
    (4,"cubr",20,70,28,None),(5,"cui",20,70,42,None),(6,"cucn",20,70,"trace",None),
    (7,"cubr2",20,70,38,"pcy3"),(8,"cubr2",20,70,37,"bpy"),(9,"cubr2",20,80,19,None),
    (10,"cubr2",20,60,47,None),(11,"cubr2",20,50,13,None),(12,"cubr2",10,60,39,None),
    (13,"cubr2",40,60,54,None),(14,"cubr2",50,60,52,None),(15,"cubr2",40,60,31,"cy2nh"),
    (16,"cubr2",40,60,68,None),(17,"cubr2",40,60,70,None),
]
for entry, copper, loading, temp, y, additive in table:
    ev = P[2] + (P[3] if entry == 17 else [])
    r = base_reaction(f"table1_entry{entry:02d}", ev, ["optimization"], f"PDF p2, Table 1 entry {entry}",
                      "Cu-catalyzed tertiary alkylation of β-allenyl silane", "oxidative coupling to 3aa")
    r["general_scheme_ids"] = [GENERAL]
    limiting = "mol_2a" if entry >= 16 else "mol_1a"
    a1 = {"number": 0.21 if entry >= 16 else 0.15, "unit": "mmol"}
    a2 = {"number": 0.15 if entry >= 16 else 0.225, "unit": "mmol"}
    oxidant_eq = 3 if entry == 17 else 2
    r["participants"] = [participant("allene", "mol_1a", "reactant", ev, amount=a1),
                         participant("malonate", "mol_2a", "reactant", ev, amount=a2),
                         participant("copper", "mol_"+copper, "catalyst", ev, molpct={"number": loading, "unit": "mol%"}, ref=limiting),
                         participant("oxidant", "mol_dtbp", "reagent", ev,
                                     amount=None if entry == 17 else {"number": 0.3, "unit": "mmol"},
                                     equiv={"number": oxidant_eq, "unit": "equiv"}, ref=limiting),
                         participant("solvent", "mol_dmso", "solvent", ev, volume={"number": 1.5, "unit": "mL"}),
                         participant("atmosphere", "mol_n2", "reagent", ev),
                         participant("product", "mol_3aa", "product", ev, product_role="main product")]
    if additive:
        dose = {"number": 0.3, "unit": "mmol"} if additive == "cy2nh" else None
        mp = None if additive == "cy2nh" else {"number": 20, "unit": "mol%"}
        r["participants"].append(participant("additive", "mol_"+additive, "additive", ev, amount=dose, molpct=mp, ref=limiting if mp else None))
    r["reaction_scale"] = F(f"{a1['number']} mmol 1a; {a2['number']} mmol 2a",
                            {"reference_molecule_id": limiting, "amount": a2 if entry >= 16 else a1}, ev)
    r["steps"] = [step("reaction", ev, [p["participant_id"] for p in r["participants"]], "conduct reaction in DMSO under N2",
                       {"number": temp, "unit": "°C"}, {"number": 12, "unit": "h"}, "N2")]
    set_solvent(r, "mol_dmso", ev, volume={"number": 1.5, "unit": "mL"})
    if y == "trace":
        r["outcomes"] = [outcome("trace_3aa", "trace", ["mol_3aa"], ev, "trace", basis="isolated yield (Table 1 footnote b)")]
    else:
        r["outcomes"] = [outcome("yield_3aa", "yield", ["mol_3aa"], ev, {"number": y, "unit": "%"},
                                  yield_type="isolated yield", basis="isolated yield")]
    if entry == 1:
        r["location"] = F("PDF p2, Scheme 2 (1a run) and Table 1 entry 1", "PDF p2, Scheme 2 (1a run) and Table 1 entry 1", ev)
        r["notes"] = [F("The 44% 1a discovery experiment is repeated as Table 1 entry 1 and is stored once.",
                          "The 44% 1a discovery experiment is repeated as Table 1 entry 1 and is stored once.", ev)]
    reactions.append(r)

# Scheme 3 scope (entry 17 supplies the standard 3aa/70% experiment); the 1 mmol run is distinct.
r = base_reaction("scheme3_3aa_scaleup", P[2] + P[3], ["substrate_scope", "scale_up"],
                  "PDF p2 prose and p3 Scheme 3 footnote a", "scale-up preparation of 3aa",
                  "oxidative coupling to 3aa")
r["general_scheme_ids"] = [GENERAL]
r["participants"] = [participant("allene", "mol_1a", "reactant", P[2]+P[3], amount={"number": 1.4, "unit": "mmol"}),
                     participant("malonate", "mol_2a", "reactant", P[2]+P[3], amount={"number": 1, "unit": "mmol"}),
                     participant("copper", "mol_cubr2", "catalyst", P[3], molpct={"number":40,"unit":"mol%"}, ref="mol_2a"),
                     participant("oxidant", "mol_dtbp", "reagent", P[3], equiv={"number":3,"unit":"equiv"}, ref="mol_2a"),
                     participant("solvent", "mol_dmso", "solvent", P[3]), participant("product", "mol_3aa", "product", P[3], product_role="main product")]
r["reaction_scale"] = F("On 1 mmol scale", {"reference_molecule_id":"mol_2a","amount":{"number":1,"unit":"mmol"}}, P[2]+P[3])
r["steps"] = [step("reaction", P[3], [p["participant_id"] for p in r["participants"]], "standard Scheme 3 conditions",
                   {"number":60,"unit":"°C"},{"number":12,"unit":"h"})]
set_solvent(r,"mol_dmso",P[3])
r["outcomes"] = [outcome("yield","yield",["mol_3aa"],P[2]+P[3],{"number":71,"unit":"%"},yield_type="reported yield")]
r["review"]["unresolved"].append("The main text does not state the measurement basis of the 1 mmol 71% yield.")
reactions.append(r)

scope3 = {"b":71,"c":65,"d":63,"e":64,"f":57,"g":56,"h":38,"i":64,"j":56,"k":61,"l":34,
          "m":70,"n":50,"o":67,"p":63,"q":55,"r":"n.d."}
ratios = {"o":"9.1:1","p":"10:1","q":"10:1"}
for letter, y in scope3.items():
    product = f"mol_3{letter}a"
    r = base_reaction(f"scheme3_3{letter}a", P[3], ["substrate_scope"], f"PDF p3, Scheme 3, 3{letter}a",
                      "variation of β-allenyl silanes", f"oxidative coupling to 3{letter}a",
                      no_change=(letter=="r"))
    r["general_scheme_ids"] = [GENERAL]
    temp = 80 if letter in "kl" else 60
    r["participants"] = [participant("allene",f"mol_1{letter}","reactant",P[3],equiv={"number":1.4,"unit":"equiv"},ref="mol_2a",method="inherited_general_condition"),
                         participant("malonate","mol_2a","reactant",P[3],equiv={"number":1,"unit":"equiv"},ref="mol_2a",method="inherited_general_condition"),
                         participant("copper","mol_cubr2","catalyst",P[3],molpct={"number":40,"unit":"mol%"},ref="mol_2a"),
                         participant("oxidant","mol_dtbp","reagent",P[3],equiv={"number":3,"unit":"equiv"},ref="mol_2a"),
                         participant("solvent","mol_dmso","solvent",P[3]), participant("product",product,"product",P[3],product_role="desired product")]
    r["steps"]=[step("reaction",P[3],[p["participant_id"] for p in r["participants"]],"Scheme 3 conditions",
                         {"number":temp,"unit":"°C"},{"number":12,"unit":"h"})]
    set_solvent(r,"mol_dmso",P[3])
    if y == "n.d.":
        r["outcomes"]=[outcome("not_detected","no_reaction",[product],P[3],"not detected",scope="desired 3ra")]
    else:
        r["outcomes"]=[outcome("yield","yield",[product],P[3],{"number":y,"unit":"%"},yield_type="reported yield")]
        r["review"]["unresolved"].append("Scheme 3 does not state the yield measurement basis for this scope entry.")
    if letter in ratios:
        r["outcomes"].append(outcome("ez_ratio","selectivity",[product],P[3],ratios[letter],basis="1H NMR analysis of the crude product",
                                             scope={"metric":"E/Z ratio","order":"E/Z","product":f"3{letter}a"}))
    reactions.append(r)

# Scheme 4 scope.
scope4 = {"b":(64,60),"c":(66,70),"d":(62,60),"e":(62,60),"f":(64,60),"g":(59,60),
          "h":(66,80),"i":(66,60),"j":(70,60),"k":(38,60),"l":(62,70),"m":(61,50),
          "n":(67,60),"o":(75,60),"p":(61,60),"q":(71,60)}
for letter,(y,temp) in scope4.items():
    sub, product = f"mol_2{letter}", f"mol_3a{letter}"
    r=base_reaction(f"scheme4_3a{letter}",P[3],["substrate_scope"],f"PDF p3, Scheme 4, 3a{letter}",
                    "variation of tertiary C(sp3)–H bonds",f"oxidative coupling to 3a{letter}")
    r["general_scheme_ids"]=[GENERAL]
    r["participants"]=[participant("allene","mol_1a","reactant",P[3],equiv={"number":1.4,"unit":"equiv"},ref=sub,method="inherited_general_condition"),
                       participant("ch_substrate",sub,"reactant",P[3],equiv={"number":1,"unit":"equiv"},ref=sub,method="inherited_general_condition"),
                       participant("copper","mol_cubr2","catalyst",P[3],molpct={"number":40,"unit":"mol%"},ref=sub),
                       participant("oxidant","mol_dtbp","reagent",P[3],equiv={"number":3,"unit":"equiv"},ref=sub),
                       participant("solvent","mol_dmso","solvent",P[3]),participant("product",product,"product",P[3],product_role="desired product")]
    r["steps"]=[step("reaction",P[3],[p["participant_id"] for p in r["participants"]],"Scheme 4 conditions",
                       {"number":temp,"unit":"°C"},{"number":12,"unit":"h"})]
    set_solvent(r,"mol_dmso",P[3])
    r["outcomes"]=[outcome("yield","yield",[product],P[3],{"number":y,"unit":"%"},yield_type="reported yield")]
    r["review"]["unresolved"].append("Scheme 4 does not state the yield measurement basis for this scope entry.")
    reactions.append(r)

# Scheme 5 derivatization arrows.
derivs=[
    ("4","mol_npmi","mol_toluene",170,24,81,("C-C","DB-C-C"),"Diels–Alder cycloaddition"),
    ("5","mol_maleic","mol_toluene",170,36,72,("C-C","DB-C-C"),"Diels–Alder cycloaddition"),
    ("6","mol_tcne","mol_toluene",160,12,71,("C-C","DB-C-C"),"Diels–Alder cycloaddition"),
    ("8","mol_bf3et2o","mol_mecn",100,10,74,("C-C",),"Lewis acid mediated intramolecular cyclization"),
]
for prod,reagent,solvent,temp,hours,y,labels,name in derivs:
    r=base_reaction(f"scheme5_product{prod}",P[3],["derivatization"],f"PDF p3, Scheme 5, product {prod}",name,name,labels=labels)
    r["participants"]=[participant("diene","mol_3aa","reactant",P[3]),participant("reagent",reagent,"reagent",P[3]),
                       participant("solvent",solvent,"solvent",P[3]),participant("product",f"mol_{prod}","product",P[3],product_role="derivatization product")]
    r["steps"]=[step("reaction",P[3],[p["participant_id"] for p in r["participants"]],name,
                       {"number":temp,"unit":"°C"},{"number":hours,"unit":"h"})]
    set_solvent(r,solvent,P[3])
    r["outcomes"]=[outcome("yield","yield",[f"mol_{prod}"],P[3],{"number":y,"unit":"%"},yield_type="reported yield")]
    r["review"]["unresolved"].append("Scheme 5 does not state equivalents, scale, workup, purification, or yield measurement basis.")
    reactions.append(r)

# LiAlH4 arrow includes the explicitly shown 0 °C -> rt program.
r=base_reaction("scheme5_product7",P[3],["derivatization"],"PDF p3, Scheme 5, product 7",
                "LiAlH4 reduction of the 1,3-diester moiety","ester reduction to 1,3-diol",labels=("C-H","O-H"))
r["participants"]=[participant("diene","mol_3aa","reactant",P[3]),participant("reagent","mol_lialh4","reagent",P[3]),
                   participant("solvent","mol_et2o","solvent",P[3]),participant("product","mol_7","product",P[3],product_role="1,3-diol")]
r["steps"]=[step("low_temperature","P3" if False else P[3],["diene","reagent","solvent"],"combine/reduce at the displayed initial temperature",
                       {"number":0,"unit":"°C"},None),
            step("warm_to_rt",P[3],["diene","reagent","solvent","product"],"warm from 0 °C to room temperature",
                       "r.t.",{"number":2,"unit":"h"})]
set_solvent(r,"mol_et2o",P[3])
r["outcomes"]=[outcome("yield","yield",["mol_7"],P[3],{"number":91,"unit":"%"},yield_type="reported yield")]
r["review"]["unresolved"] += ["The source prints '0 °C to rt, 2 h' but does not separately report addition time or ramp duration.",
                                "Scheme 5 does not state equivalents, scale, workup, purification, or yield measurement basis."]
reactions.append(r)

# Mechanistic controls reported in prose and Scheme 6a.
for suffix, omitted in [("control_no_copper","copper salt"),("control_no_dtbp","DTBP oxidant")]:
    r=base_reaction(suffix,P[3],["control","mechanistic_study"],"PDF p3, control-experiment prose",
                    f"standard reaction in the absence of {omitted}","negative control",labels=(),no_change=True)
    r["participants"]=[participant("allene","mol_1a","reactant",P[3]),participant("malonate","mol_2a","reactant",P[3]),
                       participant("solvent","mol_dmso","solvent",P[3])]
    if omitted != "copper salt": r["participants"].append(participant("copper","mol_cubr2","catalyst",P[3],molpct={"number":40,"unit":"mol%"},ref="mol_2a",method="inherited_general_condition"))
    if omitted != "DTBP oxidant": r["participants"].append(participant("oxidant","mol_dtbp","reagent",P[3],equiv={"number":3,"unit":"equiv"},ref="mol_2a",method="inherited_general_condition"))
    r["steps"]=[step("reaction",P[3],[p["participant_id"] for p in r["participants"]],f"standard conditions omitting {omitted}",
                       {"number":60,"unit":"°C"},{"number":12,"unit":"h"})]
    set_solvent(r,"mol_dmso",P[3])
    r["outcomes"]=[outcome("no_alkylated_product","no_reaction",[],P[3],"not detected",scope="alkylated product 3aa")]
    r["review"]["unresolved"].append("Detection method and limit for 'no alkylated product was detected' are not reported in the main text.")
    reactions.append(r)

# Split TEMPO and BHT because the added agents differ; the 86% structure 9 is the TEMPO adduct.
for reagent,equiv in [("tempo",3),("bht",1)]:
    ev=P[3]+P[4]
    r=base_reaction(f"scheme6_{reagent}",ev,["control","mechanistic_study"],f"PDF p3 prose and p4 Scheme 6a, {reagent.upper()} run",
                    f"radical trapping/inhibition with {reagent.upper()}","radical trapping control",
                    labels=("O-C",) if False else (("C-O",) if reagent=="tempo" else ()), no_change=(reagent=="bht"))
    r["participants"]=[participant("allene","mol_1a","reactant",ev,equiv={"number":1.4,"unit":"equiv"},ref="mol_2a"),
                       participant("malonate","mol_2a","reactant",ev,amount={"number":0.15,"unit":"mmol"}),
                       participant("trap",f"mol_{reagent}","additive",ev,equiv={"number":equiv,"unit":"equiv"},ref="mol_2a"),
                       participant("copper","mol_cubr2","catalyst",ev,molpct={"number":40,"unit":"mol%"},ref="mol_2a",method="inherited_general_condition"),
                       participant("oxidant","mol_dtbp","reagent",ev,equiv={"number":3,"unit":"equiv"},ref="mol_2a",method="inherited_general_condition"),
                       participant("solvent","mol_dmso","solvent",ev),participant("coupling_product","mol_3aa","product",ev,product_role="suppressed coupling product")]
    if reagent=="tempo": r["participants"].append(participant("adduct","mol_9","product",ev,product_role="radical-trapping adduct"))
    r["steps"]=[step("reaction",ev,[p["participant_id"] for p in r["participants"]],"apply standard conditions with radical trap/inhibitor",
                       {"number":60,"unit":"°C"},{"number":12,"unit":"h"})]
    set_solvent(r,"mol_dmso",ev)
    r["outcomes"]=[outcome("trace_3aa","trace",["mol_3aa"],ev,"trace",scope="coupling product 3aa")]
    if reagent=="tempo": r["outcomes"].append(outcome("yield_9","yield",["mol_9"],ev,{"number":86,"unit":"%"},yield_type="reported yield"))
    r["review"]["unresolved"].append("Scheme 6a combines 'TEMPO or BHT' above one arrow; adduct 9 is assigned only to the TEMPO run because the drawn product contains the TEMPO unit, while the BHT-specific product outcome is not stated.")
    reactions.append(r)

# Proposed pathway is not an experimental proof record.
r=base_reaction("scheme6b_proposed_mechanism",P[3]+P[4],["mechanistic_study"],"PDF p3–4, Scheme 6b and mechanism prose",
                "plausible reaction pathway","proposed radical addition/oxidation/β-Si elimination",
                record_type="proposed_mechanism")
r["participants"]=[participant("oxidant","mol_dtbp","reagent",P[4]),participant("malonate","mol_2a","reactant",P[4]),
                   participant("A","mol_int_A","unknown",P[4]),participant("allene","mol_1a","reactant",P[4]),
                   participant("B","mol_int_B","unknown",P[4]),participant("C","mol_int_C","unknown",P[4]),
                   participant("product","mol_3aa","product",P[4],product_role="proposed final product"),
                   participant("side_product","mol_tmsotbu","product",P[4],product_role="side product detected by GC-MS")]
r["steps"]=[step("hat",P[3]+P[4],["oxidant","malonate","A"],"thermal DTBP homolysis and proposed H-atom abstraction"),
            step("radical_addition",P[4],["A","allene","B"],"proposed radical addition to the central allene carbon"),
            step("set",P[4],["B","C"],"proposed Cu(II)-mediated single-electron oxidation"),
            step("beta_si_elimination",P[4],["C","product","side_product"],"proposed β-Si elimination")]
r["notes"]=[F("TMSOtBu side-product formation was stated to be validated by GC-MS; the supporting data were not provided.",
                "TMSOtBu side-product formation was stated to be validated by GC-MS; the supporting data were not provided.",P[4])]
r["review"]["unresolved"].append("Arrow-pushing intermediates A–C are proposed, not experimentally established; GC-MS evidence for TMSOtBu is only stated and not shown in the supplied PDF.")
reactions.append(r)

bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id":f"{DOC}:scheme1_prior_art","document_id":DOC,"evidence_ids":P[1],"kind":"literature_context",
     "note":"Scheme 1a–b summarizes cited prior C–H dienylation and radical-allene selectivity concerns; not current experiments."},
    {"context_id":f"{DOC}:graphical_abstract","document_id":DOC,"evidence_ids":P[1],"kind":"duplicate_summary",
     "note":"Graphical abstract repeats the current generic transformation and is not an independent experiment."},
    {"context_id":f"{DOC}:si_notice","document_id":DOC,"evidence_ids":P[4],"kind":"availability_notice",
     "note":"The article states that experimental procedures, characterization data, and NMR spectra are in SI; no SI file was supplied or explicitly associated."},
    {"context_id":f"{DOC}:gcms_statement","document_id":DOC,"evidence_ids":P[4],"kind":"characterization_statement",
     "note":"TMSOtBu was reportedly detected by GC-MS; no peak table or spectrum is present in the supplied pages, so no separate synthetic experiment was created."},
    {"context_id":f"{DOC}:references_page5","document_id":DOC,"evidence_ids":P[5],"kind":"references",
     "note":"Page 5 contains references only and no current-work reaction experiment."},
]

# Update page tasks to attempted without claiming semantic verification.
tasks=json.loads((RUN/"tasks.json").read_text())
for t in tasks:
    t["state"]="attempted"
    t["visual_evidence_id"]=f"{DOC}:page{t['evidence_id'].split(':page')[1].split(':')[0]}:image"
(RUN/"tasks.json").write_text(json.dumps(tasks,ensure_ascii=False,indent=2)+"\n")

inventory={"document_id":DOC,"pages_reviewed":[1,2,3,4,5],"molecule_count":len(bundle["molecules"]),
           "reaction_count":len(bundle["reactions"]),"reaction_ids":[r["reaction_id"] for r in bundle["reactions"]],
           "notes":["Raw machine inventory; formal import remains needs_review.","No SI associated."]}
(RUN/"raw_model/inventory.json").write_text(json.dumps(inventory,ensure_ascii=False,indent=2)+"\n")
(RUN/"submitted.json").write_text(json.dumps(bundle,ensure_ascii=False,indent=2)+"\n")
print(json.dumps({"molecules":len(bundle["molecules"]),"reactions":len(bundle["reactions"])},ensure_ascii=False))
