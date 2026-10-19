#!/usr/bin/env python3
"""Source transcription for doc_26447b9b5c1dd197, PDF pages 1--6.

This is a session-local machine extraction.  It deliberately leaves uncertain
structures and unavailable SI details unresolved; the standard merge will
quarantine every record as needs_review.
"""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_26447b9b5c1dd197"
RUN = ROOT / "runs/current"
DATASET = st.read(RUN / "dataset.json")
BUNDLE = st.empty(DATASET["documents"])
P = {i: f"{DOC}:page{i}:image" for i in range(1, 7)}
T = {i: f"{DOC}:page{i}:text" for i in range(1, 8)}


def f(raw=None, norm=None, ev=(), method="original_extraction", derivation=None,
      status=None, verification=None, alternatives=()):
    if norm is None and raw is not None:
        norm = copy.deepcopy(raw)
    return st.fact(raw, norm, ev, method=method, derivation=derivation,
                   status=status, verification=verification,
                   alternatives=alternatives)


def missing(ev=()):
    return st.fact(evidence=ev)


def derived(norm, ev, why):
    return f(None, norm, ev, method="structure_derivation", derivation=why)


def inferred(norm, ev, why):
    return f(None, norm, ev, method="inference", derivation=why)


MOLECULES = {}


def molecule(label, name=None, ev=(P[2],), smiles=None, stereo=None):
    """Create/reuse a paper-namespaced molecule with explicit stereo states."""
    mid = f"{DOC}:mol:{label.replace(' ', '_').replace('/', '_')}"
    if mid in MOLECULES:
        occ = MOLECULES[mid]["occurrences"][0]
        for eid in ev:
            if eid not in occ["evidence_ids"]:
                occ["evidence_ids"].append(eid)
        return mid
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    m["names"]["original"] = f(name or label, name or label, ev)
    m["names"]["english"] = f(name, name, ev) if name else missing(ev)
    m["names"]["abbreviations"] = f([label], [label], ev)
    m["names"]["aliases"] = missing(ev)
    m["cas_numbers"] = missing(ev)
    m["identity_note"] = f(
        "Source label/drawing retained; CAS was not reported in the supplied main-text pages.",
        "Source label/drawing retained; structure remains a machine candidate unless a SMILES is given.", ev)
    if smiles:
        diag = st.stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = derived(smiles, ev, "Candidate graph transcribed from the displayed structure; parsing does not verify identity.")
        m["stereochemistry"]["identity_key"] = derived(
            diag["stereo_graph_key"], ev, "Computed by standard stereo from the candidate isomeric graph.")
        m["stereochemistry"]["isomeric_smiles"] = derived(
            diag["canonical_isomeric_smiles"], ev, "Canonical isomeric SMILES computed from the candidate graph.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing(ev)
        m["generic_template"] = None
        m["stereochemistry"]["identity_key"] = missing(ev)
        m["stereochemistry"]["isomeric_smiles"] = missing(ev)
    states = stereo or {
        "chirality_status": "unresolved",
        "geometry_status": "unresolved",
        "representation_status": "unresolved",
    }
    for key, value in states.items():
        if key in ("chirality_status", "geometry_status", "representation_status"):
            m["stereochemistry"][key] = inferred(
                value, ev, "Conservative status assigned after inspection of the source drawing; no unreported stereoisomer was invented.")
        elif key in m["stereochemistry"]:
            m["stereochemistry"][key] = f(value, value, ev)
    MOLECULES[mid] = m
    return mid


SIMPLE = {
    "1a": ("N-benzyl-N-tosyl allenamide 1a", "C=C=CN(Cc1ccccc1)S(=O)(=O)c1ccc(C)cc1",
           {"chirality_status": "achiral", "geometry_status": "not_applicable", "representation_status": "not_applicable",
            "axial": {"allene": "terminal CH2; no stereogenic allene axis in candidate graph"}}),
    "2a": ("bromobenzoyloxylation N,O-aminal 2a", "C=C(Br)C(OC(=O)c1ccccc1)N(Cc1ccccc1)S(=O)(=O)c1ccc(C)cc1",
           {"chirality_status": "unspecified", "geometry_status": "not_applicable", "representation_status": "not_encoded",
            "atom_centers": [{"site": "aminal carbon", "configuration": "unspecified in source"}]}),
    "3a": ("phenylacetylene", "C#Cc1ccccc1",
           {"chirality_status": "achiral", "geometry_status": "not_applicable", "representation_status": "not_applicable"}),
    "4aa": ("(Z)-skipped bromo enyne 4aa",
            "N(/C=C(\\Br)CC#Cc1ccccc1)(Cc1ccccc1)S(=O)(=O)c1ccc(C)cc1",
            {"chirality_status": "achiral", "geometry_status": "specified", "representation_status": "encoded",
             "double_bonds": [{"site": "enamide C=C", "configuration": "Z", "assignment": "source states exclusive Z"}],
             "assignment": "Z assignment is stated in the main text; assignment details are deferred to unavailable SI."}),
}
for label, (name, smi, stereo) in SIMPLE.items():
    molecule(label, name, (P[2],), smi, stereo)


def unresolved_product(label, ev, name=None, geometry=False, axial=False, natural=False):
    stereo = {
        "chirality_status": "specified" if natural else ("unspecified" if axial else "unresolved"),
        "geometry_status": "specified" if geometry else "unresolved",
        "representation_status": "not_encoded",
    }
    if geometry:
        stereo["double_bonds"] = [{"site": "displayed enamide C=C", "configuration": "Z"}]
    if axial:
        stereo["axial"] = [{"site": "displayed allene axis", "configuration": "unspecified"}]
    if natural:
        stereo["atom_centers"] = "Stereochemical bonds are drawn in the natural-product fragment; individual centers were not transcribed."
    return molecule(label, name, ev, stereo=stereo)


def common(label, name=None, ev=(P[2],), smiles=None):
    stereo = {"chirality_status": "achiral", "geometry_status": "not_applicable",
              "representation_status": "not_applicable"} if smiles else None
    return molecule(label, name or label, ev, smiles, stereo)


# Common condition components.  Candidate SMILES are supplied only for simple,
# unambiguous structures; complexes and source-defined screens remain unresolved.
for args in [
    ("NBS", "N-bromosuccinimide", "O=C1CCC(=O)N1Br"),
    ("NCS", "N-chlorosuccinimide", "O=C1CCC(=O)N1Cl"),
    ("NIS", "N-iodosuccinimide", "O=C1CCC(=O)N1I"),
    ("BzOH", "benzoic acid", "O=C(O)c1ccccc1"),
    ("MeOH", "methanol", "CO"),
    ("AcOH", "acetic acid", "CC(=O)O"),
    ("CH2Cl2", "dichloromethane", "ClCCl"),
    ("DIPEA", "Hünig's base (DIPEA)", "CCN(C(C)C)C(C)C"),
    ("Et3N", "triethylamine", "CCN(CC)CC"),
    ("DBU", "DBU", None),
    ("i-Pr2NH", "diisopropylamine", "CC(C)NC(C)C"),
    ("BF3.OEt2", "boron trifluoride diethyl etherate", None),
    ("ZnBr2", "zinc bromide", "[Br-].[Br-].[Zn+2]"),
    ("Li-3a", "lithium phenylacetylide", "[Li+].[C-]#Cc1ccccc1"),
    ("THF", "tetrahydrofuran", "C1CCOC1"),
    ("KOtBu", "potassium tert-butoxide", "[K+].CC(C)(C)[O-]"),
    ("H2", "hydrogen", "[H][H]"),
    ("quinoline", "quinoline", "c1ccc2ncccc2c1"),
    ("toluene", "toluene", "Cc1ccccc1"),
    ("K2CO3", "potassium carbonate", "[K+].[K+].[O-]C([O-])=O"),
    ("H2O", "water", "O"),
    ("1,4-dioxane", "1,4-dioxane", "O1CCOCC1"),
]:
    common(args[0], args[1], (P[2], P[4]), args[2])


def participant(pid, mid, role, ev, product_role=None):
    p = st.new_participant(pid, mid, role, ev)
    if product_role:
        p["product_role"] = f(product_role, product_role, ev)
    return p


def set_charge(p, field, raw, norm, ev, ref_mid=None):
    p["charges"][field] = f(raw, norm, ev)
    if field in ("equivalents", "mol_percent"):
        p["equivalence_reference"] = f(ref_mid, ref_mid, ev)


def reaction(rid, ev, record_type="specific_experiment", categories=("other",),
             location=None, name=None, transform=None):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = f(list(categories), list(categories), ev)
    r["location"] = f(location or rid, location or rid, ev)
    r["reaction_name"] = f(name, name, ev) if name else missing(ev)
    r["transformation_type"] = inferred(transform, ev, "Transformation class inferred from the displayed reactant/product connectivity.") if transform else missing(ev)
    r["review"]["unresolved"] = [
        "Machine transcription; chemical identity and every field require expert source review.",
        "No SI was supplied or associated; procedures, purification and characterization remain incomplete.",
    ]
    return r


def add_step(r, sid, ev, pids, operation=None, temp=None, time=None, atmosphere=None, pressure=None, stage="reaction"):
    s = st.new_step(sid, ev, stage)
    s["participant_ids"] = list(pids)
    if operation is not None:
        s["operation"] = f(operation, operation, ev)
    if temp is not None:
        raw, norm = temp if isinstance(temp, tuple) else (temp, temp)
        s["temperature"] = f(raw, norm, ev)
    if time is not None:
        raw, norm = time if isinstance(time, tuple) else (time, time)
        s["time"] = f(raw, norm, ev)
    if atmosphere is not None:
        s["atmosphere"] = f(atmosphere, atmosphere, ev)
    if pressure is not None:
        raw, norm = pressure if isinstance(pressure, tuple) else (pressure, pressure)
        s["pressure"] = f(raw, norm, ev)
    r["steps"].append(s)
    return s


def outcome(r, oid, kind, targets, ev, raw=None, norm=None, basis=None, product_role="product", scope=None, yield_type=None):
    o = st.new_outcome(oid, kind, targets, ev)
    o["product_role"] = f(product_role, product_role, ev)
    o["value"] = f(raw, norm, ev) if raw is not None or norm is not None else missing(ev)
    o["measurement_basis"] = f(basis, basis, ev) if basis else missing(ev)
    o["scope"] = f(scope, scope, ev) if scope is not None else missing(ev)
    o["yield_type"] = f(yield_type, yield_type, ev) if yield_type else missing(ev)
    r["outcomes"].append(o)
    return o


def classify(r, ev, changes=None, not_applicable=False, limitations=()):
    c = r["classification"]
    c["mapped_reaction_smiles"] = missing(ev)
    c["bond_changes"] = []
    labels = {}
    if not_applicable:
        c["labels"] = inferred([], ev, "No observed target transformation; no formation label assigned.")
        c["coverage"] = inferred("not_applicable", ev, "No observed target transformation.")
        c["limitations"] = list(limitations) or ["Negative/failed experiment; intended product bonds were not assigned."]
        return
    for i, (a, b, before, after, ctype, note) in enumerate(changes or (), 1):
        change = {"change_id": f"b{i}", "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
                  "before_order": before, "after_order": after, "change_type": ctype, "site_note": note}
        c["bond_changes"].append(derived(change, ev, "Partial bond change inferred from the displayed net connectivity; exact atom mapping was not asserted."))
        if ctype in ("formed", "order_increased", "hydrogen_gain"):
            lab = st.bond_label(a, b, after)
            labels.setdefault(lab["code"], dict(lab, supporting_change_ids=[]))["supporting_change_ids"].append(f"b{i}")
    c["labels"] = derived(list(labels.values()), ev, "Taxonomy labels calculated from the listed partial bond changes.")
    c["coverage"] = inferred("partial", ev, "Only visually secure net changes were encoded; no complete atom mapping or coproduct balance.")
    c["limitations"] = list(limitations) or ["No complete atom mapping; hydrogen/coproduct balance and exact changed-atom identities remain unresolved."]


CC = ("C", "C", 0, 1, "formed", "new substrate-to-alkynyl carbon connection")
CBr = ("C", "Br", 0, 1, "formed", "new carbon-bromine connection")
CCl = ("C", "Cl", 0, 1, "formed", "new carbon-chlorine connection")
CI = ("C", "I", 0, 1, "formed", "new carbon-iodine connection")
CO = ("C", "O", 0, 1, "formed", "new carbon-oxygen connection")
DBCC = ("C", "C", 1, 2, "order_increased", "new/promoted C=C bond in product")
DBCN = ("C", "N", 1, 2, "order_increased", "iminium C=N bond formation")
TBCC = ("C", "C", 2, 3, "order_increased", "new/promoted C#C bond")
CH = ("C", "H", 0, 1, "hydrogen_gain", "net C-H formation")
BREAK_CO = ("C", "O", 1, 0, "broken", "C-O leaving-group bond cleavage")
BREAK_CBR = ("C", "Br", 1, 0, "broken", "C-Br bond cleavage")
BREAK_CH = ("C", "H", 1, 0, "hydrogen_loss", "net C-H cleavage")
DOWN_CC_2_1 = ("C", "C", 2, 1, "order_decreased", "C=C bond-order decrease")
DOWN_CC_3_2 = ("C", "C", 3, 2, "order_decreased", "C#C bond-order decrease to C=C")
DOWN_CC_3_1 = ("C", "C", 3, 1, "order_decreased", "C#C bond-order decrease to C-C")
ALKYNYL = (CC, DBCC, BREAK_CO, DOWN_CC_2_1)


def add_standard_second_step(r, substrate_mid, alkyne_mid, product_mids, ev, scale=0.1, alkyne_amount=0.13,
                             time_h=4, omit=(), modifications=()):
    psub = participant("substrate", substrate_mid, "reactant", ev)
    set_charge(psub, "amount", f"{scale} mmol", {"number": scale, "unit": "mmol"}, ev)
    r["participants"].append(psub)
    palk = participant("alkyne", alkyne_mid, "reactant", ev)
    set_charge(palk, "amount", f"{alkyne_amount} mmol", {"number": alkyne_amount, "unit": "mmol"}, ev)
    r["participants"].append(palk)
    for pid, mid in zip([f"product{i+1}" for i in range(len(product_mids))], product_mids):
        r["participants"].append(participant(pid, mid, "product", ev, "major" if len(product_mids) == 1 else "reported regioisomer"))
    defaults = [("znbr2", common("ZnBr2"), "catalyst", "mol_percent", "10 mol%", {"number": 10, "unit": "mol%"}),
                ("bf3", common("BF3.OEt2"), "reagent", "equivalents", "1.1 equiv", {"number": 1.1, "unit": "equiv"}),
                ("base", common("DIPEA"), "reagent", "equivalents", "2.5 equiv", {"number": 2.5, "unit": "equiv"})]
    for pid, mid, role, field, raw, norm in defaults:
        if pid in omit:
            continue
        p = participant(pid, mid, role, ev)
        set_charge(p, field, raw, norm, ev, substrate_mid)
        r["participants"].append(p)
    solv = participant("solvent", common("CH2Cl2"), "solvent", ev)
    r["participants"].append(solv)
    r["reaction_scale"] = f(f"{scale} mmol of isolated N,O-aminal", {"reference_molecule_id": substrate_mid, "amount": {"number": scale, "unit": "mmol"}}, ev)
    r["solvent_conditions"]["composition"] = f("CH2Cl2", [common("CH2Cl2")], ev)
    r["solvent_conditions"]["concentration"] = f("0.1 M", {"number": 0.1, "unit": "M"}, ev)
    r["solvent_conditions"]["concentration_basis"] = f("substrate N,O-aminal", substrate_mid, ev)
    pids = [p["participant_id"] for p in r["participants"] if p["role"] != "product"]
    add_step(r, "reaction", ev, pids, operation="standard second-step alkynylation", temp="rt",
             time=(f"{time_h} h", {"number": time_h, "unit": "h"}))
    r["procedure"]["general_procedure_reference"] = f("Scheme 2 standard condition", "Scheme 2 standard condition", ev)
    if modifications:
        r["procedure"]["modifications"] = f(list(modifications), list(modifications), ev)


# Scheme 1C: current-work generic transformation.  Scheme 1A/B are retained
# below as prior-work/source context, not current-paper experiments.
r = reaction("scheme1c_generic_haloalkynylation", (P[2],), "generic_scheme", ("other",),
             "PDF page 2, Scheme 1C", "Synthesis of functionalizable skipped halo enynes", "formal haloalkynylation")
for pid, mid, role in [("allene", molecule("allenamide_generic", "generic EWG-substituted allenamide", (P[2],)), "reactant"),
                       ("alkyne", molecule("terminal_alkyne_generic", "generic terminal alkyne", (P[2],)), "reactant"),
                       ("product", molecule("halo_enyne_generic", "generic functionalizable skipped halo enyne", (P[2],)), "product")]:
    r["participants"].append(participant(pid, mid, role, (P[2],), "generic product" if role == "product" else None))
classify(r, (P[2],), (CC, CBr), limitations=("Generic scheme covers X = Cl, Br, I and EWG/alkyne variables; no exact atom map.",))
BUNDLE["reactions"].append(r)


# Scheme 2A: first-stage N,O-aminal syntheses.
first_stage = [
    ("2a", "NBS", "BzOH", 89),
    ("2a-i", "NBS", "p-MeO-BzOH", 73),
    ("2a-ii", "NBS", "m-NO2-BzOH", 84),
    ("2a-iii", "NBS", "AcOH", 68),
    ("2a-iv", "NCS", "MeOH", 76),
    ("2a-v", "NBS", "MeOH", 89),
    ("2a-vi", "NIS", "MeOH", 93),
]
for prod, nxs, roh, y in first_stage:
    candidate = {
        "2a-v": "C=C(Br)C(OC)N(Cc1ccccc1)S(=O)(=O)c1ccc(C)cc1",
    }.get(prod)
    prod_mid = molecule(prod, f"N,O-aminal {prod}", (P[2],), smiles=candidate, stereo={
        "chirality_status": "unspecified", "geometry_status": "not_applicable", "representation_status": "not_encoded",
        "atom_centers": [{"site": "new aminal carbon", "configuration": "unspecified"}]}) if prod != "2a" else common("2a")
    acid_mid = common(roh, roh, (P[2],))
    r = reaction(f"scheme2_first_step_{prod}", (P[2],), categories=("optimization",),
                 location=f"PDF page 2, Scheme 2A, product {prod}", name="Synthesis of N,O-aminal",
                 transform="halocarboxylation/haloalkoxylation of allenamide")
    for pid, mid, role in [("allene", common("1a"), "reactant"), ("nxs", common(nxs), "reagent"),
                           ("roh", acid_mid, "reagent"), ("solvent", common("CH2Cl2"), "solvent"),
                           ("product", prod_mid, "product")]:
        p = participant(pid, mid, role, (P[2],), "product" if role == "product" else None)
        if pid == "allene": set_charge(p, "amount", "0.2 mmol", {"number": 0.2, "unit": "mmol"}, (P[2],))
        if pid in ("nxs", "roh"): set_charge(p, "equivalents", "1.1 equiv", {"number": 1.1, "unit": "equiv"}, (P[2],), common("1a"))
        r["participants"].append(p)
    r["reaction_scale"] = f("1a (0.2 mmol)", {"reference_molecule_id": common("1a"), "amount": {"number": 0.2, "unit": "mmol"}}, (P[2],))
    r["solvent_conditions"]["composition"] = f("CH2Cl2", [common("CH2Cl2")], (P[2],))
    r["solvent_conditions"]["concentration"] = f("0.1 M", {"number": 0.1, "unit": "M"}, (P[2],))
    r["solvent_conditions"]["concentration_basis"] = f("1a", common("1a"), (P[2],))
    add_step(r, "halofunctionalization", (P[2],), ["allene", "nxs", "roh", "solvent"],
             "combine under Scheme 2 first-step conditions", ("0 °C", {"number": 0, "unit": "°C"}),
             ("5 min", {"number": 5, "unit": "min"}))
    outcome(r, "yield", "yield", [prod_mid], (P[2],), f"{y}%", {"number": y, "unit": "%"},
            basis=None, yield_type="first-step yield")
    r["review"]["unresolved"].append("Scheme 2A does not state the measurement/isolation basis for this first-step yield.")
    classify(r, (P[2],), (CBr if nxs == "NBS" else CCl if nxs == "NCS" else CI, CO, DOWN_CC_2_1))
    BUNDLE["reactions"].append(r)


def second_step_record(suffix, substrate_label, product_label, conv=None, yields=(), modifications=(), grouped=False,
                       expected=True, alkyne_label="3a"):
    ev = (P[2],)
    sub = common(substrate_label)
    product = common(product_label) if product_label == "4aa" else unresolved_product(product_label, ev, geometry=True)
    r = reaction(f"scheme2_table_{suffix}", ev, categories=("optimization",),
                 location=f"PDF page 2, Scheme 2B, entry {suffix}", name="Optimization of the second step",
                 transform="alkynylation of N,O-aminal")
    add_standard_second_step(r, sub, common(alkyne_label), [product], ev, modifications=modifications)
    if grouped:
        r["review"]["unresolved"].append("This source row aggregates several mutually exclusive screening runs; the main text does not expose individual values, so it is retained as one row-level record.")
    if conv is not None:
        outcome(r, "conversion", "conversion", [sub], ev, conv[0], conv[1], basis="reported table conversion", product_role="reactant")
    for oid, raw, norm, basis in yields:
        outcome(r, oid, "yield", [product], ev, raw, norm, basis=basis, yield_type="yield of product from isolated N,O-aminal")
    if expected:
        classify(r, ev, ALKYNYL)
    else:
        classify(r, ev, not_applicable=True)
    BUNDLE["reactions"].append(r)
    return r


# Entry 1 and 2 have prose evidence in addition to table values.
r = second_step_record("1", "2a", "4aa", ("100%", {"number": 100, "unit": "%"}),
                       (("yield_nmr", "96%", {"number": 96, "unit": "%"}, "NMR using CH2Br2 internal standard"),
                        ("yield_isolated", "[94%]", {"number": 94, "unit": "%"}, "isolated")))
outcome(r, "selectivity_z", "selectivity", [common("4aa")], (P[2],), "exclusive Z-selectivity", "exclusive Z-selectivity",
        basis="stereochemistry assignment details deferred to Supporting Information",
        scope={"metric": "alkene geometry", "major_isomer": "Z"})
r["notes"].append(f("Table reports 96% NMR and 94% isolated yield; both are retained with their own measurement basis.",
                    "96% NMR and 94% isolated are distinct measurements, not a chosen-value conflict.", (P[2],)))

r = second_step_record("2", "2a-iii", "4aa", ("100%", {"number": 100, "unit": "%"}),
                       (("yield_table_nmr", "73%", {"number": 73, "unit": "%"}, "NMR using CH2Br2 internal standard"),
                        ("yield_prose_isolated", "isolated in 73% yield", {"number": 73, "unit": "%"}, "isolated")),
                       ("N,O-aminal 2a-iii used instead",))
r["notes"].append(f("The table footnote gives an NMR basis while the prose calls the same 73% value isolated; both bases are retained separately.",
                    "same numeric value, distinct stated measurement bases", (P[2],)))
for suffix, sub, prod, y in [("3_2a-i", "2a-i", "4aa", 87), ("3_2a-ii", "2a-ii", "4aa", 93),
                             ("4_2a-iv", "2a-iv", "5a", 73), ("4_2a-v", "2a-v", "4aa", 81),
                             ("4_2a-vi", "2a-vi", "5b", 67)]:
    second_step_record(suffix, sub, prod, ("100%", {"number": 100, "unit": "%"}),
                       (("yield_nmr", f"{y}%", {"number": y, "unit": "%"}, "NMR using CH2Br2 internal standard"),),
                       (f"{sub} used instead of 2a",))

# Aggregated screening rows 5--10.  Alternatives are represented by a single
# explicitly named screen entity so they are not falsely encoded as co-added.
screen_rows = [
    ("5", "Zn-salt screen", "ZnCl2, ZnI2 or Zn(OTf)2 instead of ZnBr2", "5–47%", {"minimum": 5, "maximum": 47, "unit": "%", "raw_text": "5–47%"}, "up to 43% (maximum from ZnI2)", {"number": 43, "unit": "%", "qualifier": "le", "raw_text": "up to 43%"}),
    ("6", "Cu-salt screen", "Cu(I)Cl/BrI/TfC or Cu(OTf)2 instead of ZnBr2", "0–100%", {"minimum": 0, "maximum": 100, "unit": "%", "raw_text": "0–100%"}, "up to 81% (maximum from CuBr)", {"number": 81, "unit": "%", "qualifier": "le", "raw_text": "up to 81%"}),
    ("7", "Lewis/Brønsted acid screen", "BF3·OEt2 replaced by TMSOTf, TfOH or p-TSA", "0–62%", {"minimum": 0, "maximum": 62, "unit": "%", "raw_text": "0–62%"}, "up to 16% (maximum from TfOH)", {"number": 16, "unit": "%", "qualifier": "le", "raw_text": "up to 16%"}),
    ("8", "base screen", "Hünig's base replaced by Et3N, DBU or i-Pr2NH", "32–45%", {"minimum": 32, "maximum": 45, "unit": "%", "raw_text": "32–45%"}, "up to 43% (maximum from DBU)", {"number": 43, "unit": "%", "qualifier": "le", "raw_text": "up to 43%"}),
    ("9", "solvent screen", "other solvents instead of CH2Cl2", "0–100%", {"minimum": 0, "maximum": 100, "unit": "%", "raw_text": "0–100%"}, "up to 67% (maximum from acetonitrile)", {"number": 67, "unit": "%", "qualifier": "le", "raw_text": "up to 67%"}),
    ("10", "BF3 loading screen", "lowering equivalent of BF3·OEt2 from 1.1 to 0.2", "10–58%", {"minimum": 10, "maximum": 58, "unit": "%", "raw_text": "10–58%"}, "up to 54% (maximum from 0.75 equiv)", {"number": 54, "unit": "%", "qualifier": "le", "raw_text": "up to 54%"}),
]
for entry, label, mod, craw, cnorm, yraw, ynorm in screen_rows:
    r = second_step_record(entry, "2a", "4aa", (craw, cnorm),
                           (("yield_nmr", yraw, ynorm, "NMR using CH2Br2 internal standard"),),
                           (mod,), grouped=True)
    screen_mid = molecule(label, label, (P[2],))
    r["participants"].append(participant("screen_alternatives", screen_mid, "unknown", (P[2],)))
    if entry == "10":
        bf = next(p for p in r["participants"] if p["participant_id"] == "bf3")
        set_charge(bf, "equivalents", "0.2–1.1 equiv screen; maximum at 0.75 equiv",
                   {"minimum": 0.2, "maximum": 1.1, "unit": "equiv", "raw_text": "0.2–1.1 equiv screen"}, (P[2],), common("2a"))

r = second_step_record("11", "2a", "4aa", ("100%", {"number": 100, "unit": "%"}),
                       (("yield_nmr", "0%", {"number": 0, "unit": "%"}, "NMR using CH2Br2 internal standard"),),
                       ("Li-acetylide of 3a used instead of ZnBr2",), expected=False)
r["participants"].append(participant("li_acetylide", common("Li-3a"), "reagent", (P[2],)))
outcome(r, "unidentified_mixture", "decomposition", [], (P[2],), "unidentified mixture of compounds", "unidentified mixture of compounds", basis=None)

# Entry 12 is an overall sequence from 1a.  Replace the isolated 2a participant
# with 1a and retain the basic-workup qualification.
r = second_step_record("12", "1a", "4aa", ("100%", {"number": 100, "unit": "%"}),
                       (("yield_nmr", "67%", {"number": 67, "unit": "%"}, "NMR using CH2Br2 internal standard"),
                        ("yield_isolated", "[66%]", {"number": 66, "unit": "%"}, "isolated")),
                       ("from 1a (after basic work-up of 2a)",))
r["review"]["unresolved"].append("Basic-workup reagent and operation details are not reported in the supplied main text.")
r["notes"].append(f("Entry 12 reports 67% NMR and 66% isolated yield after basic workup; neither value replaces the other.",
                    "67% NMR and 66% isolated retained as distinct measurements", (P[2],)))
classify(r, (P[2],), (CBr, CC), limitations=("Overall sequence from 1a; transient C-O formation is not classified as a net product bond.",))

# Prose-only sequential one-pot experiment without successful imide removal.
r = reaction("scheme2_onepot_33pct", (P[2],), categories=("optimization",),
             location="PDF page 2, prose below Scheme 2", name="Sequential one-pot haloalkynylation", transform="formal bromoalkynylation")
for pid, mid, role in [("allene", common("1a"), "reactant"), ("alkyne", common("3a"), "reactant"), ("product", common("4aa"), "product")]:
    r["participants"].append(participant(pid, mid, role, (P[2],), "product" if role == "product" else None))
outcome(r, "yield", "yield", [common("4aa")], (P[2],), "33%", {"number": 33, "unit": "%"}, basis=None, yield_type="overall one-pot yield")
r["review"]["unresolved"] += ["Exact charges and measurement basis are not reported in the supplied prose.", "Authors attribute the poor result to residual imide."]
add_step(r, "sequential_onepot", (P[2],), ["allene", "alkyne"],
         "sequential one-pot reaction; exact stage conditions are not separately reported in the prose", None, None)
classify(r, (P[2],), (CBr, CC))
BUNDLE["reactions"].append(r)

for omit, suffix in [("znbr2", "13_without_ZnBr2"), ("base", "13_without_Hunigs_base")]:
    r = second_step_record(suffix, "2a", "4aa", ("100%", {"number": 100, "unit": "%"}),
                           (("yield_nmr", "0%", {"number": 0, "unit": "%"}, "NMR using CH2Br2 internal standard"),),
                           (f"without {omit}",), expected=False)
    r["participants"] = [p for p in r["participants"] if p["participant_id"] != omit]
    r["steps"][0]["participant_ids"] = [pid for pid in r["steps"][0]["participant_ids"] if pid != omit]


# Scheme 3 substrate scope.  Alkyne and product drawings are retained as
# source-identified unresolved structures except the representative 3a/4aa.
scope3 = [
    ("3a", "4aa", 94, 4), ("3b", "4ab", 92, 4), ("3c", "4ac", 89, 4), ("3d", "4ad", 90, 4),
    ("3e", "4ae", 91, 4), ("3f", "4af", 95, 4), ("3g", "4ag", 93, 4), ("3h", "4ah", 86, 4),
    ("3i", "4ai", 85, 4), ("3j", "4aj", 88, 4), ("3k", "4ak", 90, 4), ("3l", "4al", 92, 4),
    ("3m", "4am", 86, 4), ("3n", "4an", 83, 4), ("3o", "4ao", 76, 12), ("3p", "4ap", 61, 12),
    ("3q", "4aq", 90, 4), ("3r", "4ar", 89, 4), ("3s", "4as", 87, 4),
]
for alk, prod, y, hours in scope3:
    alk_mid = common("3a") if alk == "3a" else molecule(alk, f"alkyne {alk}", (P[2],))
    prod_mid = common("4aa") if prod == "4aa" else unresolved_product(prod, (P[2],), geometry=True)
    r = reaction(f"scheme3_scope_{prod}", (P[2],), categories=("substrate_scope",),
                 location=f"PDF page 2, Scheme 3, {prod}", name="Formal bromoalkynylation of 1a", transform="alkynylation of N,O-aminal")
    add_standard_second_step(r, common("2a"), alk_mid, [prod_mid], (P[2],), scale=0.2, alkyne_amount=0.26, time_h=hours)
    outcome(r, "yield", "yield", [prod_mid], (P[2],), f"{y}%", {"number": y, "unit": "%"}, basis="isolated", yield_type="isolated yield")
    classify(r, (P[2],), ALKYNYL)
    BUNDLE["reactions"].append(r)

# Separate 5 mmol 4aa scale-up, distinct from the 0.2 mmol scope point.
r = reaction("scheme3_scaleup_4aa_5mmol", (P[2], P[3]), categories=("scale_up",),
             location="PDF pages 2–3, Scheme 3 footnote c and prose", name="Scale-up bromoalkynylation to 4aa", transform="alkynylation of N,O-aminal")
add_standard_second_step(r, common("2a"), common("3a"), [common("4aa")], (P[2], P[3]), scale=5.0, alkyne_amount=6.5, time_h=4)
outcome(r, "yield", "yield", [common("4aa")], (P[2], P[3]), "96%", {"number": 96, "unit": "%"}, basis="isolated", yield_type="isolated yield")
r["notes"].append(f("without chromatographic purification", "without chromatographic purification", (P[3],)))
classify(r, (P[2], P[3]), ALKYNYL)
BUNDLE["reactions"].append(r)


# Scheme 4: allene-substituent/halogen scope, regioisomer controls and failures.
scheme4_main = [("2b", "4ba", 73), ("2c", "4ca", 90), ("2d", "4da", 91), ("2e", "4ea", 91)]
for sub, prod, y in scheme4_main:
    sub_mid = molecule(sub, f"N,O-aminal {sub}", (P[3],))
    prod_mid = unresolved_product(prod, (P[3],), geometry=True)
    r = reaction(f"scheme4_{prod}", (P[3],), categories=("substrate_scope",),
                 location=f"PDF page 3, Scheme 4A, {prod}", name="Formal bromoalkynylation with phenylacetylene", transform="alkynylation of N,O-aminal")
    add_standard_second_step(r, sub_mid, common("3a"), [prod_mid], (P[3],))
    outcome(r, "yield", "yield", [prod_mid], (P[3],), f"{y}%", {"number": y, "unit": "%"}, basis=None)
    outcome(r, "rr", "selectivity", [prod_mid], (P[3],), ">49:1", ">49:1", basis="1H NMR", scope={"metric": "r.r.", "order": "γ:α"})
    outcome(r, "z_e", "selectivity", [prod_mid], (P[3],), ">49:1", ">49:1", basis="1H NMR", scope={"metric": "Z:E", "order": "Z:E"})
    if prod == "4ba":
        r["notes"].append(f("Alkynylation without isolation of unstable bromobenzoyloxylation product; reported yield is overall.",
                            "overall yield without isolation of unstable intermediate", (P[3],)))
    r["review"]["unresolved"].append("Scheme 4 does not state a yield measurement basis for this value.")
    classify(r, (P[3],), ALKYNYL)
    BUNDLE["reactions"].append(r)

for x, y, rr, hal in [("Cl", 87, "55:45", "C-Cl"), ("Br", 84, "36:64", "C-Br"), ("I", 27, "<1:49", "C-I")]:
    sub_mid = molecule(f"2f-{x}", f"carbazole-derived N,O-aminal, X={x}", (P[3],))
    gamma = unresolved_product(f"gamma-4fa-{x}", (P[3],), geometry=True)
    alpha = unresolved_product(f"alpha-4fa-{x}", (P[3],), geometry=False)
    r = reaction(f"scheme4_1f_{x}", (P[3], P[4]), categories=("substrate_scope",),
                 location=f"PDF page 3, Scheme 4B, 1f X={x}", name="Haloalkynylation of carbazole allenamine", transform="regiodivergent alkynylation")
    add_standard_second_step(r, sub_mid, common("3a"), [gamma, alpha], (P[3], P[4]))
    yield_targets = [alpha] if x == "I" else [gamma, alpha]
    outcome(r, "yield", "yield", yield_targets, (P[3],), f"{y}%", {"number": y, "unit": "%"}, basis="NMR" if x == "I" else None,
            product_role="α-regioisomer" if x == "I" else "regioisomer mixture")
    outcome(r, "rr", "selectivity", [gamma, alpha], (P[3],), rr, rr, basis="1H NMR", scope={"metric": "r.r.", "order": "γ:α"})
    if x != "I":
        outcome(r, "z_e", "selectivity", [gamma], (P[3],), ">49:1", ">49:1", basis="1H NMR", scope={"metric": "Z:E", "order": "Z:E"})
    r["review"]["unresolved"].append("Regioisomer structures are retained separately, but no definitive SMILES or individual isolated amounts are available in the supplied pages.")
    classify(r, (P[3],), ALKYNYL)
    BUNDLE["reactions"].append(r)

sub_mid = molecule("2g", "allenyl-ether-derived N,O-aminal", (P[3],))
gamma = unresolved_product("gamma-4ga", (P[3],), geometry=True)
alpha = unresolved_product("alpha-4ga", (P[3],))
r = reaction("scheme4_1g_Br", (P[3], P[4]), categories=("substrate_scope",),
             location="PDF page 3, Scheme 4B, 1g", name="Bromoalkynylation of allenyl ether", transform="regiodivergent alkynylation")
add_standard_second_step(r, sub_mid, common("3a"), [gamma, alpha], (P[3], P[4]))
outcome(r, "yield", "yield", [gamma, alpha], (P[3],), "71%", {"number": 71, "unit": "%"}, basis=None, product_role="regioisomer mixture")
outcome(r, "rr", "selectivity", [gamma, alpha], (P[3],), "46:54", "46:54", basis="1H NMR", scope={"metric": "r.r.", "order": "γ:α"})
outcome(r, "z_e", "selectivity", [gamma], (P[3],), ">49:1", ">49:1", basis="1H NMR", scope={"metric": "Z:E", "order": "Z:E"})
classify(r, (P[3],), ALKYNYL)
BUNDLE["reactions"].append(r)

for sub, first_ok in [("1h", True), ("1i", False), ("1j", False), ("1k", False)]:
    sub_mid = molecule(sub, f"allene substrate {sub}", (P[3],), stereo={"chirality_status": "unresolved", "geometry_status": "not_applicable", "representation_status": "not_encoded"})
    r = reaction(f"scheme4_failed_{sub}", (P[3], P[4]), categories=("control",),
                 location=f"PDF page 3, Scheme 4C, {sub}", name="Substrate limitation test", transform=None)
    for pid, mid, role in [("allene", sub_mid, "reactant"), ("nbs", common("NBS"), "reagent"), ("bzoh", common("BzOH"), "reagent"), ("solvent", common("CH2Cl2"), "solvent")]:
        r["participants"].append(participant(pid, mid, role, (P[3], P[4])))
    if first_ok:
        r["participants"].append(participant("alkyne", common("3a"), "reactant", (P[3], P[4])))
        add_step(r, "halofunctionalization", (P[3],), ["allene", "nbs", "bzoh", "solvent"],
                 "formation of N,O-aminal under Scheme 2 first-step conditions",
                 ("0 °C", {"number": 0, "unit": "°C"}), ("5 min", {"number": 5, "unit": "min"}))
        add_step(r, "attempted_alkynylation", (P[3], P[4]), ["allene", "alkyne"],
                 "attempted alkynylation under the referenced standard condition", "rt",
                 ("4 h", {"number": 4, "unit": "h"}))
        outcome(r, "no_alkynylation", "no_reaction", [], (P[3],), "N,O-aminal was successfully formed; alkynylation did not furnish a reported product", True, basis=None)
    else:
        add_step(r, "failed_halofunctionalization", (P[3], P[4]), ["allene", "nbs", "bzoh", "solvent"],
                 "attempted bromobenzoyloxylation under Scheme 2 first-step conditions",
                 ("0 °C", {"number": 0, "unit": "°C"}), ("5 min", {"number": 5, "unit": "min"}))
        outcome(r, "failed_first_step", "no_reaction", [], (P[3],), "bromobenzoyloxylation failed; alkynylation was unable to perform", True, basis=None)
    classify(r, (P[3],), not_applicable=True)
    BUNDLE["reactions"].append(r)

for sub, prod, y, hours, halchange in [("2a-iv", "5a", 91, 4, CCl), ("2a-vi", "5b", 78, 12, CI)]:
    prod_mid = unresolved_product(prod, (P[3],), geometry=True)
    r = reaction(f"scheme4_other_halogen_{prod}", (P[3], P[4]), categories=("substrate_scope",),
                 location=f"PDF page 3, Scheme 4D, {prod}", name="Formal haloalkynylation", transform="alkynylation of N,O-aminal")
    add_standard_second_step(r, common(sub), common("3a"), [prod_mid], (P[3], P[4]), time_h=hours)
    outcome(r, "yield", "yield", [prod_mid], (P[3],), f"{y}%", {"number": y, "unit": "%"}, basis=None)
    r["review"]["unresolved"].append("Yield measurement basis is not reported in the supplied main text.")
    classify(r, (P[3],), ALKYNYL)
    BUNDLE["reactions"].append(r)

one_l = molecule("1l", "gamma-substituted allenamide 1l", (P[3],), stereo={
    "chirality_status": "unspecified", "geometry_status": "not_applicable", "representation_status": "not_encoded",
    "axial": [{"site": "substituted allene axis", "configuration": "not assigned"}]})
six = unresolved_product("6", (P[3],), geometry=False)
r = reaction("scheme4_1l_elimination_6", (P[3], P[4]), categories=("substrate_scope",),
             location="PDF page 3, Scheme 4E", name="Vinylogous elimination", transform="vinylogous elimination")
for pid, mid, role in [("allene", one_l, "reactant"), ("nbs", common("NBS"), "reagent"), ("bzoh", common("BzOH"), "reagent"),
                       ("alkyne", common("3a"), "reactant"), ("product", six, "product")]:
    r["participants"].append(participant(pid, mid, role, (P[3], P[4]), "elimination product" if role == "product" else None))
outcome(r, "yield", "yield", [six], (P[3],), "42%", {"number": 42, "unit": "%"}, basis=None)
r["notes"].append(f("Reaction performed after basic workup of N,O-aminal.", "basic workup before alkynylation", (P[3],)))
r["review"]["unresolved"] += ["Basic-workup details and yield measurement basis are not reported.", "Phenylacetylene is charged under the standard sequence but is not incorporated into displayed product 6."]
add_step(r, "halofunctionalization", (P[3],), ["allene", "nbs", "bzoh"],
         "form N,O-aminal under the referenced first-step conditions",
         ("0 °C", {"number": 0, "unit": "°C"}), ("5 min", {"number": 5, "unit": "min"}))
add_step(r, "basic_workup", (P[3],), ["allene"], "basic workup; reagent and duration unreported", None, None, stage="workup")
add_step(r, "alkynylation_elimination", (P[3], P[4]), ["allene", "alkyne"],
         "apply referenced standard alkynylation condition; elimination product 6 observed", "rt",
         ("4 h", {"number": 4, "unit": "h"}))
classify(r, (P[3],), (CBr,), limitations=("Only the visually secure new C-Br bond was classified; elimination bond-order/hydrogen changes remain unresolved.",))
BUNDLE["reactions"].append(r)


# Figure 2 mechanistic NMR experiments: actual treatments are experimental
# records; peak lists/HRMS values remain source_context below.
iminium = molecule("A", "alpha,beta-unsaturated N-tosyliminium intermediate A", (P[1], P[3]), stereo={
    "chirality_status": "achiral", "geometry_status": "unspecified", "representation_status": "not_encoded"})
r = reaction("figure2_nmr_2a_to_A", (P[3], T[4]), categories=("mechanistic_study",),
             location="PDF page 3 Figure 2 and page 4 prose", name="NMR generation of N-tosyliminium intermediate A", transform="Lewis-acid ionization")
for pid, mid, role in [("substrate", common("2a"), "reactant"), ("bf3", common("BF3.OEt2"), "reagent"),
                       ("cdcl3", common("CDCl3", "deuterated chloroform", (P[3],)), "solvent"), ("intermediate", iminium, "product")]:
    p = participant(pid, mid, role, (P[3], T[4]), "observed intermediate" if role == "product" else None)
    if pid == "bf3": set_charge(p, "equivalents", "1.1 equiv", {"number": 1.1, "unit": "equiv"}, (P[3], T[4]), common("2a"))
    r["participants"].append(p)
r["solvent_conditions"]["composition"] = f("CDCl3", [common("CDCl3")], (P[3],))
r["solvent_conditions"]["concentration"] = f("0.3 M", {"number": 0.3, "unit": "M"}, (P[3],))
r["special_conditions"]["reactor"] = f("sealed NMR tube", "sealed NMR tube", (P[3],))
add_step(r, "nmr_treatment", (P[3], T[4]), ["substrate", "bf3", "cdcl3"], "treat in sealed NMR tube", "rt", ("60 min", {"number": 60, "unit": "min"}))
outcome(r, "intermediate_observed", "other", [iminium], (P[3], T[4]), "new 1H NMR peak at 9.25 ppm; vinylic protons shift to 6.91 ppm", "intermediate A observed by NMR", basis="1H/13C NMR and HSQC")
classify(r, (P[3],), (DBCN, BREAK_CO), limitations=("Exact ionic/counterion mapping and remaining bond-order changes are unresolved.",))
BUNDLE["reactions"].append(r)

two_av = common("2a-v")
r = reaction("figure2_nmr_2a-v_complex", (P[3], T[4]), categories=("mechanistic_study",),
             location="PDF page 4 prose discussing Figure 2", name="Analogous BF3 treatment of methanol-derived N,O-aminal 2a-v", transform="Lewis-acid treatment")
for pid, mid, role in [("substrate", two_av, "reactant"), ("bf3", common("BF3.OEt2"), "reagent")]:
    p = participant(pid, mid, role, (P[3], T[4]))
    if pid == "bf3": set_charge(p, "equivalents", "analogous experiment; exact amount not separately restated", None, (T[4],), two_av)
    r["participants"].append(p)
add_step(r, "nmr_treatment", (T[4],), ["substrate", "bf3"], "analogous NMR experiment", "rt", None)
outcome(r, "complex_spectrum", "other", [], (T[4],), "complex spectrum", "complex spectrum", basis="NMR")
r["review"]["unresolved"].append("The analogous experiment's exact concentration, duration and charge are not separately reported in prose; they were not inherited as exact facts.")
classify(r, (T[4],), not_applicable=True, limitations=("Complex spectrum; no single observed product connectivity assigned.",))
BUNDLE["reactions"].append(r)


# Scheme 5A natural-product-bound alkyne partners.
for alk, prod, y, source in [("3t", "4at", 73, "(+)-menthol"), ("3u", "4au", 67, "cholesterol")]:
    alk_mid = molecule(alk, f"propargylic ether derived from {source}", (P[4],), stereo={
        "chirality_status": "specified", "geometry_status": "not_applicable", "representation_status": "not_encoded",
        "atom_centers": f"source drawing retains the stereochemical bonds of {source}; centers not individually transcribed"})
    prod_mid = unresolved_product(prod, (P[4],), geometry=True, natural=True, name=f"bromo enynamide {prod} derived from {source}")
    r = reaction(f"scheme5_scope_{prod}", (P[4], P[5]), categories=("substrate_scope",),
                 location=f"PDF page 4, Scheme 5A, {prod}", name="Natural-product-bound bromoalkynylation", transform="alkynylation of N,O-aminal")
    add_standard_second_step(r, common("2a"), alk_mid, [prod_mid], (P[4], P[5]))
    outcome(r, "yield", "yield", [prod_mid], (P[4],), f"{y}% (from {source})", {"number": y, "unit": "%"}, basis=None)
    r["review"]["unresolved"].append("Scheme 5A does not state the yield measurement basis; natural-product stereocenters were not converted to candidate SMILES.")
    classify(r, (P[4],), ALKYNYL)
    BUNDLE["reactions"].append(r)


def scheme5_reaction(label, product_name, y, reagents, solvent_label, temp, time, changes, notes=(), geometry=False, axial=False):
    ev = (P[4], P[5])
    product_mid = unresolved_product(label, ev, name=product_name, geometry=geometry, axial=axial)
    r = reaction(f"scheme5_{label}", ev, categories=("derivatization",),
                 location=f"PDF page 4, Scheme 5B, product {label}", name="Synthetic transformation of 4aa", transform=product_name)
    r["participants"].append(participant("substrate", common("4aa"), "reactant", ev))
    for i, spec in enumerate(reagents, 1):
        lab, role, field, raw, norm, ref = spec
        mid = common(lab, lab, ev)
        p = participant(f"reagent{i}", mid, role, ev)
        if field:
            set_charge(p, field, raw, norm, ev, ref or common("4aa"))
        r["participants"].append(p)
    solv_mid = common(solvent_label, solvent_label, ev)
    r["participants"].append(participant("solvent", solv_mid, "solvent", ev))
    r["participants"].append(participant("product", product_mid, "product", ev, "product"))
    r["reaction_scale"] = f("0.1 mmol of 4aa", {"reference_molecule_id": common("4aa"), "amount": {"number": 0.1, "unit": "mmol"}}, ev)
    r["solvent_conditions"]["composition"] = f(solvent_label, [solv_mid], ev)
    rawt, normt = temp
    rawtime, normtime = time
    add_step(r, "reaction", ev, [p["participant_id"] for p in r["participants"] if p["role"] != "product"],
             "perform under Scheme 5 condition", (rawt, normt), (rawtime, normtime))
    outcome(r, "yield", "yield", [product_mid], ev, f"{y}%", {"number": y, "unit": "%"}, basis=None)
    for note in notes:
        r["notes"].append(f(note, note, ev))
    r["review"]["unresolved"].append("Scheme 5 does not state the yield measurement/isolation basis or workup/purification details.")
    classify(r, ev, changes)
    BUNDLE["reactions"].append(r)
    return r


r7 = scheme5_reaction("7", "allene enamide 7", 86,
    [("Et3N", "reagent", "equivalents", "1.5 equiv", {"number": 1.5, "unit": "equiv"}, None)],
    "THF", ("rt", "rt"), ("16 h", {"number": 16, "unit": "h"}), (DBCC, CH, DOWN_CC_3_2, BREAK_CH), axial=True)

# Product 8 has an explicit 0 °C -> rt program.  Replace the single generated
# step with ordered initial and warming/overall-duration steps without inventing
# a stage-specific time allocation.
r8 = scheme5_reaction("8", "skipped diynamide 8", 73,
    [("KOtBu", "reagent", "equivalents", "1.1 equiv", {"number": 1.1, "unit": "equiv"}, None)],
    "THF", ("0 °C to rt", "rt"), ("2 h", {"number": 2, "unit": "h"}), (TBCC, BREAK_CBR, BREAK_CH))
r8["steps"] = []
add_step(r8, "initial_0C", (P[4],), ["substrate", "reagent1", "solvent"], "KOtBu treatment at initial temperature",
         ("0 °C", {"number": 0, "unit": "°C"}), None)
add_step(r8, "warm_to_rt_overall", (P[4],), ["substrate", "reagent1", "solvent"], "warm from 0 °C to rt; duration is reported only for the overall program",
         "rt", ("2 h overall for 0 °C to rt program", {"number": 2, "unit": "h"}))
r8["review"]["unresolved"].append("The source does not allocate the 2 h between the 0 °C stage and warming/room-temperature stage, and gives no ramp time.")

r9 = scheme5_reaction("9", "cyclobutene amide 9", 67,
    [("Et3N", "reagent", "equivalents", "1.5 equiv", {"number": 1.5, "unit": "equiv"}, None)],
    "THF", ("80 °C", {"number": 80, "unit": "°C"}), ("24 h", {"number": 24, "unit": "h"}), (CC,),
    notes=("one-pot base-induced isomerization and 4pi-electrocyclization", "structure confirmed by X-ray crystallography, CCDC 2019531"))

r10 = scheme5_reaction("10", "full hydrogenation product 10", 91,
    [("Lindlar catalyst", "catalyst", "mass", "cal. 5 mg", {"number": 5, "unit": "mg", "qualifier": "approx", "raw_text": "cal. 5 mg"}, None),
     ("H2", "reagent", None, None, None, None)],
    "toluene", ("rt", "rt"), ("16 h", {"number": 16, "unit": "h"}), (CH, CH, CH, CH, DOWN_CC_3_1), geometry=True)
r10["steps"][0]["atmosphere"] = f("H2 balloon", "H2", (P[4],))
r10["steps"][0]["pressure"] = f("1 atm", "1 atm", (P[4],))

r11 = scheme5_reaction("11", "semihydrogenation product 11", 93,
    [("Lindlar catalyst", "catalyst", "mass", "cal. 5 mg", {"number": 5, "unit": "mg", "qualifier": "approx", "raw_text": "cal. 5 mg"}, None),
     ("quinoline", "additive", "mol_percent", "20 mol%", {"number": 20, "unit": "mol%"}, None),
     ("H2", "reagent", None, None, None, None)],
    "toluene", ("rt", "rt"), ("12 h", {"number": 12, "unit": "h"}), (CH, CH, DOWN_CC_3_2), geometry=True)
r11["steps"][0]["atmosphere"] = f("H2 balloon", "H2", (P[4],))
r11["steps"][0]["pressure"] = f("1 atm", "1 atm", (P[4],))
r11["solvent_conditions"]["concentration"] = f("0.3 M", {"number": 0.3, "unit": "M"}, (P[4],))
r11["solvent_conditions"]["concentration_basis"] = f("4aa", common("4aa"), (P[4],))

r12 = scheme5_reaction("12", "allene-appended conjugated ynenamide 12", 73,
    [("3-ethynyl anisole", "reactant", "equivalents", "1.3 equiv", {"number": 1.3, "unit": "equiv"}, None),
     ("Pd(PPh3)2Cl2", "catalyst", "mol_percent", "5 mol%", {"number": 5, "unit": "mol%"}, None),
     ("CuI", "catalyst", "mol_percent", "10 mol%", {"number": 10, "unit": "mol%"}, None),
     ("Et3N", "reagent", None, None, None, None)],
    "Et3N/THF (1:2)", ("rt", "rt"), ("12 h", {"number": 12, "unit": "h"}), (CC, DBCC, CH, BREAK_CBR, DOWN_CC_3_2, BREAK_CH), axial=True,
    notes=("one-pot two-step transformation",))
r12["participants"] = [p for p in r12["participants"] if p["participant_id"] != "solvent"]
r12["participants"].append(participant("solvent_thf", common("THF"), "solvent", (P[4], P[5])))
r12["steps"][0]["participant_ids"] = ["solvent_thf" if pid == "solvent" else pid for pid in r12["steps"][0]["participant_ids"]]
r12["solvent_conditions"]["composition"] = f("Et3N/THF", [common("Et3N"), common("THF")], (P[4],))
r12["solvent_conditions"]["mixture_ratio"] = f("Et3N/THF (1:2)", {"components": ["Et3N", "THF"], "ratio": "1:2"}, (P[4],))

r13 = scheme5_reaction("13", "arylated skipped ynenamide 13", 71,
    [("3-methoxyphenylboronic acid", "reactant", "equivalents", "1.5 equiv", {"number": 1.5, "unit": "equiv"}, None),
     ("Pd(PPh3)2Cl2", "catalyst", "mol_percent", "5 mol%", {"number": 5, "unit": "mol%"}, None),
     ("K2CO3", "reagent", "equivalents", "3 equiv", {"number": 3, "unit": "equiv"}, None)],
    "1,4-dioxane/H2O (2:1)", ("80 °C", {"number": 80, "unit": "°C"}), ("24 h", {"number": 24, "unit": "h"}), (CC, BREAK_CBR), geometry=True,
    notes=("stereospecific Suzuki coupling",))
r13["participants"] = [p for p in r13["participants"] if p["participant_id"] != "solvent"]
r13["participants"].append(participant("solvent_dioxane", common("1,4-dioxane"), "solvent", (P[4], P[5])))
r13["participants"].append(participant("solvent_water", common("H2O"), "solvent", (P[4], P[5])))
r13["steps"][0]["participant_ids"] = [pid for pid in r13["steps"][0]["participant_ids"] if pid != "solvent"] + ["solvent_dioxane", "solvent_water"]
r13["solvent_conditions"]["composition"] = f("1,4-dioxane/H2O", [common("1,4-dioxane"), common("H2O")], (P[4],))
r13["solvent_conditions"]["mixture_ratio"] = f("1,4-dioxane/H2O (2:1)", {"components": ["1,4-dioxane", "H2O"], "ratio": "2:1"}, (P[4],))


# Non-synthetic observations, prior work and bibliography-only regions.
BUNDLE["source_context"] = [
    {"context_id": f"{DOC}:context:bibliography", "document_id": DOC, "kind": "bibliographic_metadata",
     "title": "Characterization and Utilization of the Elusive alpha,beta-Unsaturated N-Tosyliminium: the Synthesis of Highly Functionalizable Skipped Halo Enynes",
     "authors": ["Tapas R. Pradhan", "Mukti Paudel", "Jordan L. Harper", "Paul Ha-Yeon Cheong", "Jin Kyoon Park"],
     "journal": "Organic Letters", "year": 2021, "volume": 23, "pages": "1427-1433", "doi": "10.1021/acs.orglett.1c00103",
     "evidence_ids": [P[1], T[1]]},
    {"context_id": f"{DOC}:context:graphical_abstract", "document_id": DOC, "kind": "current_work_summary",
     "note": "Graphical abstract summarizes the current formal haloalkynylation and downstream products; it is not counted as an additional experiment.", "evidence_ids": [P[1]]},
    {"context_id": f"{DOC}:context:scheme1_prior_work", "document_id": DOC, "kind": "prior_work",
     "note": "Scheme 1A contains general/cited iminium chemistry and Scheme 1B the authors' previous alkynylations; these are not current-paper experiments.", "evidence_ids": [P[2]]},
    {"context_id": f"{DOC}:context:figure1_computation", "document_id": DOC, "kind": "computation",
     "note": "Figure 1 reports computed free energies of chlorinated intermediates; no synthetic experiment record created.", "evidence_ids": [P[3], T[4]]},
    {"context_id": f"{DOC}:context:figure2_spectra", "document_id": DOC, "kind": "characterization",
     "note": "Figure 2 NMR/HSQC/HRMS spectra support intermediate A. Peak positions and m/z values are characterization context, not separate synthesis reactions.",
     "reported_values": {"1H_NMR": ["9.25 ppm", "5.73 and 5.91 ppm to 6.91 ppm"], "ESI_HRMS_mz": ["378.0158", "380.0162"]},
     "evidence_ids": [P[3], T[4]]},
    {"context_id": f"{DOC}:context:figure3_computation", "document_id": DOC, "kind": "computation",
     "note": "Figure 3 gives computed 1,2/1,4-addition transition states; proposed/computed pathways are not observed synthetic experiments.", "evidence_ids": [P[3], T[4]]},
    {"context_id": f"{DOC}:context:si_notice", "document_id": DOC, "kind": "source_availability",
     "note": "Page 5 states that SI contains procedures, characterization and spectra, but no SI file was supplied or explicitly associated in this session.", "evidence_ids": [P[5], T[5]]},
    {"context_id": f"{DOC}:context:ccdc", "document_id": DOC, "kind": "crystallography_notice",
     "note": "CCDC 2019531 is cited for product 9; no CIF was supplied or associated.", "evidence_ids": [P[4], P[5]]},
    {"context_id": f"{DOC}:context:page5_text", "document_id": DOC, "kind": "discussion_and_metadata",
     "note": "Page 5 prose repeats Scheme 5 transformations and then gives conclusions, SI/CCDC notices, authors and references; no additional reaction arrow was found.", "evidence_ids": [P[5]]},
    {"context_id": f"{DOC}:context:page6_references", "document_id": DOC, "kind": "bibliography_only",
     "note": "Page 6 contains cited literature references only; reaction wording in titles was not treated as current-paper experiments.", "evidence_ids": [P[6]]},
]

BUNDLE["molecules"] = list(MOLECULES.values())

# Raw inventory/diagnostics are deliberately separate from the authoritative
# submitted bundle.
inventory = {
    "document_id": DOC,
    "scope_pages": [1, 2, 3, 4, 5, 6],
    "molecule_count": len(BUNDLE["molecules"]),
    "reaction_count": len(BUNDLE["reactions"]),
    "reaction_ids": [r["reaction_id"] for r in BUNDLE["reactions"]],
    "source_context_ids": [c["context_id"] for c in BUNDLE["source_context"]],
    "chemical_accuracy_verified": False,
}
st.write(RUN / "raw_model/extraction_inventory.json", inventory)
diagnostics = {}
for label in ("1a", "2a", "2a-v", "3a", "4aa"):
    mid = f"{DOC}:mol:{label}"
    diagnostics[label] = st.stereo_inspect(st.value(MOLECULES[mid]["smiles"]))
st.write(RUN / "raw_model/stereo_diagnostics.json", diagnostics)
mentions = st.source_mentions(DATASET, ["2a", "3a", "4aa", "5b", "isolated yield", "Supporting Information"])
st.write(RUN / "raw_model/source_mentions.json", mentions)
st.write(RUN / "submitted.json", BUNDLE)
print(json.dumps({"molecules": len(BUNDLE["molecules"]), "reactions": len(BUNDLE["reactions"]),
                  "source_context": len(BUNDLE["source_context"])}, ensure_ascii=False, indent=2))
