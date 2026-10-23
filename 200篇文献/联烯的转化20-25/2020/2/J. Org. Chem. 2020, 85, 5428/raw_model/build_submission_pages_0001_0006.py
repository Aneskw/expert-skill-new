#!/usr/bin/env python3
"""Build the bounded page-1--6 v3.1 extraction for doc_a413938d5770a9e8."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
import standard_records as st

DOC = "doc_a413938d5770a9e8"
P = {n: f"{DOC}:page{n}:text" for n in range(1, 7)}
PI = {n: f"{DOC}:page{n}:image" for n in range(1, 7)}
bundle = json.loads((RUN / "submitted.json").read_text())
bundle["documents"] = json.loads((RUN / "dataset.json").read_text())["documents"]
bundle.update(molecules=[], reactions=[], source_context=[])
molecules = {}


def F(raw=None, norm=None, ev=(), method="original_extraction", derivation=None, status=None):
    return st.fact(raw, norm, ev, method=method, derivation=derivation, status=status)


def mol(mid, label, ev, name=None, smiles=None, stereo=None, note=None):
    ev = list(dict.fromkeys(ev))
    m = st.new_molecule(mid, DOC, ev)
    m["occurrences"][0]["compound_label"] = F(label, label, ev)
    m["names"]["original"] = F(name or label, name or label, ev)
    m["names"]["english"] = F(name or label, name or label, ev)
    m["identity_note"] = F(note or "Source label/name retained; no external identity lookup.",
                           note or "Source label/name retained; no external identity lookup.", ev)
    if smiles:
        m["identity_status"] = "source_identified"
        m["representation_type"] = "smiles"
        m["smiles"] = F(label, smiles, ev, method="structure_derivation",
                         derivation="Connectivity candidate transcribed from the drawing/systematic name; parsing is not source verification.")
        for key, value, why in (
            ("chirality_status", "unresolved", "Missing source assignment is not achiral/racemic."),
            ("geometry_status", "unresolved", "Only source-explicit geometry is promoted."),
            ("representation_status", "partial", "Connectivity candidate does not claim unreported stereo.")):
            m["stereochemistry"][key] = F(value, value, ev, method="inference", derivation=why)
    if stereo:
        for key, value in stereo.items():
            if key in {"chirality_status", "geometry_status", "representation_status"}:
                raw, norm, why = value
                m["stereochemistry"][key] = F(raw, norm, ev, method="inference", derivation=why)
            else:
                m["stereochemistry"][key] = F(value, value, ev)
    molecules[mid] = m
    bundle["molecules"].append(m)
    return mid


def part(pid, mid, role, ev, stage="reaction", **charges):
    p = st.new_participant(pid, mid, role, ev, stage=stage)
    if role == "product":
        p["product_role"] = F("product", "product", ev)
    for key in ("mass", "amount", "volume", "equivalents", "mol_percent"):
        if charges.get(key) is not None:
            raw, norm = charges[key]
            p["charges"][key] = F(raw, norm, ev)
    ref = charges.get("reference")
    if ref:
        p["equivalence_reference"] = F("relative to " + ref, ref, ev)
    return p


def step(sid, stage, pids, ev, operation, temperature=None, time=None, atmosphere=None):
    s = st.new_step(sid, ev, stage=stage)
    s["participant_ids"] = list(pids)
    s["operation"] = F(operation, operation, ev)
    if temperature:
        s["temperature"] = F(temperature[0], temperature[1], ev)
    if time:
        s["time"] = F(time[0], time[1], ev)
    if atmosphere:
        s["atmosphere"] = F(atmosphere, atmosphere, ev)
    return s


def out(oid, kind, targets, ev, raw, norm, yield_type=None, basis=None, scope=None):
    o = st.new_outcome(oid, kind, targets, ev)
    o["value"] = F(raw, norm, ev)
    if yield_type:
        o["yield_type"] = F(yield_type, yield_type, ev)
    if basis:
        o["measurement_basis"] = F(basis, basis, ev)
    if scope is not None:
        o["scope"] = F(scope, scope, ev)
    if kind in {"yield", "product_amount"}:
        o["product_role"] = F("product", "product", ev)
    return o


def classify(r, ev, mode):
    c = r["classification"]
    c["model_votes"] = []
    if mode == "none":
        c["coverage"] = F("negative/no observed transformation", "not_applicable", ev, method="inference",
                          derivation="No intended-product bond formation is classified for a negative control.")
        c["labels"] = F([], [], ev, method="inference", derivation="No observed transformation.")
        return
    maps = {
        "cascade": [("b1", "formed", 0, 1, "C", "C"), ("b2", "formed", 0, 1, "C", "C"),
                    ("b3", "order_increased", 1, 2, "C", "C")],
        "cc": [("b1", "formed", 0, 1, "C", "C")],
        "co": [("b1", "formed", 0, 1, "C", "O")],
        "co2": [("b1", "formed", 0, 1, "C", "O"), ("b2", "formed", 0, 1, "C", "O")],
        "ch": [("b1", "hydrogen_gain", 0, 1, "C", "H")],
        "cb": [("b1", "formed", 0, 1, "C", "B")],
    }
    changes = [{"change_id": x[0], "change_type": x[1], "before_order": x[2], "after_order": x[3],
                "atoms": [{"atom_map": None, "element": x[4]}, {"atom_map": None, "element": x[5]}]}
               for x in maps[mode]]
    labels = {}
    for ch in changes:
        label = st.bond_label(ch["atoms"][0]["element"], ch["atoms"][1]["element"], ch["after_order"])
        labels.setdefault(label["code"], {**label, "supporting_change_ids": []})["supporting_change_ids"].append(ch["change_id"])
    why = "Partial net-bond inventory inferred from the displayed arrow/mechanism; no complete atom mapping is asserted."
    c["bond_changes"] = [F(None, ch, ev, method="structure_derivation", derivation=why) for ch in changes]
    c["labels"] = F(None, list(labels.values()), ev, method="structure_derivation", derivation=why)
    c["coverage"] = F("partial", "partial", ev, method="inference", derivation=why)
    c["limitations"] = ["No source-verified atom mapping.", "Hydrogen/coproduct balance may be incomplete."]


def rxn(suffix, ev, record_type="specific_experiment", categories=(), name=None, location=None, mode="cascade"):
    ev = list(dict.fromkeys(ev))
    r = st.new_reaction(f"{DOC}:{suffix}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = F(list(categories), list(categories), ev)
    r["location"] = F(location or suffix, location or suffix, ev)
    if name:
        r["reaction_name"] = F(name, name, ev)
        r["transformation_type"] = F(name, name, ev)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
                   "unresolved": ["Machine extraction; chemistry requires expert review."]}
    r["completeness"] = "partial"
    classify(r, ev, mode)
    bundle["reactions"].append(r)
    return r


def solvent(r, mid, ev, volume=None):
    r["solvent_conditions"]["composition"] = F(molecules[mid]["names"]["original"]["raw_value"], [mid], ev)
    if volume:
        r["solvent_conditions"]["total_volume"] = F(volume[0], volume[1], ev)


def q(number, unit):
    return {"number": number, "unit": unit}


# Molecule registry: conservative candidates for named allenyne starting materials;
# compact scope/product drawings remain unresolved.
core_smiles = {
    "1a": "COC(=O)C(CC#C)(C=C=C(C)C)C(=O)OC",
    "1b": "CCOC(=O)C(CC#C)(C=C=C(C)C)C(=O)OCC",
    "1c": "COC(=O)C(CC#C)(C=C=C1CCCCC1)C(=O)OC",
    "1d": "COC(=O)C(CC#C)(C=C=C(C)C(C)(C)C)C(=O)OC",
    "1ax": "OCC(CO)(CC#C)C=C=C(C)C", "1cx": "OCC(CO)(CC#C)C=C=C1CCCCC1",
    "1dx": "OCC(CO)(CC#C)C=C=C(C)C(C)(C)C",
    "1e": "C#CCC(COC)(COC)C=C=C(C)C", "1f": "C#CCC(COC)(COC)C=C=C1CCCCC1",
    "1g": "C#CCC(COC)(COC)C=C=C(C)C(C)(C)C",
    "1h": "C#CCC1(C=C=C(C)C)COC(C)(C)OC1", "1i": "C#CCC1(C=C=C(C)C(C)(C)C)COC(C)(C)OC1",
    "1j": "C#CCC(COCc1ccccc1)(COCc1ccccc1)C=C=C(C)C",
}
for label, smi in core_smiles.items():
    pages = [P[2], P[3]] if label in {"1a", "1b", "1c", "1d"} else [P[2], P[5], P[6]]
    mol("m_" + label, label, pages, smiles=smi)
for label in ("1d", "1g", "1i"):
    m = molecules["m_" + label]
    ev = m["occurrences"][0]["evidence_ids"]
    m["stereochemistry"]["axial"] = F("potentially stereogenic allene axis; no assignment",
        [{"type": "allene_axis", "assignment": "unspecified"}], ev)
    m["stereochemistry"]["chirality_status"] = F("allene axial chirality unassigned", "unspecified", ev,
        method="inference", derivation="Unlike substituents occur at both allene termini, but no R_a/S_a or P/M assignment is reported.")

mol("m_2a", "2a", [P[1], P[2], P[3], P[4], P[6]], note="Enallene drawing retained; complete independently checked SMILES not established.")
for label in ["2" + x for x in "klmnopqrstuv"]:
    mol("m_" + label, label, [P[4]], note="Table 3 enallene row identity; full graph unresolved.")
mol("m_2a_d6", "2a-d6", [P[4]], note="KIE isotopologue; isotope placement not encoded as complete SMILES.")

product_stereo = {
    "double_bonds": [{"source_assignment": "middle double bond Z", "scope": "cross-conjugated polyene"}],
    "geometry_status": ("middle double bond Z; remaining geometry not completely encoded", "partially_specified",
                        "The Results text explicitly assigns the middle double bond as Z."),
    "representation_status": ("source geometry retained textually", "not_encoded",
                              "Complex product drawing was not independently transcribed to isomeric SMILES."),
}
for label in ["3" + x for x in "abcdefghijklmnopqrstuv"]:
    pages = [P[2], P[3]] if label[-1] in "abcdefghij" else [P[4]]
    if label in {"3a", "3b", "3c"}:
        pages.append(P[6])
    stereo = dict(product_stereo)
    if label == "3u":
        stereo["mixture"] = {"type": "E/Z mixture", "ratio": None}
        stereo["geometry_status"] = ("E/Z mixture", "mixture", "Table 3 marks 3u as E/Z without a readable ratio.")
    mol("m_" + label, label, pages, stereo=stereo, note="Cross-conjugated product drawing retained; full graph unresolved.")
mol("m_3a_d5", "3a-d5", [P[4]], stereo=product_stereo, note="KIE product isotopologue; full isotope-resolved graph unresolved.")
for label in ("4a", "4b", "5"):
    mol("m_" + label, label, [P[2], P[3]], note="Source drawing retained; full graph unresolved.")
for label in ("int_1", "int_2", "int_3", "int_4"):
    mol("m_" + label, label.replace("_", "-").title(), [P[4], P[5]], note="Author-proposed intermediate, not an isolated compound.")

common = {
    "pd_oac2": ("Pd(OAc)2", "palladium(II) acetate", "CC(=O)[O-].CC(=O)[O-].[Pd+2]"),
    "bq": ("BQ", "1,4-benzoquinone", "O=C1C=CC(=O)C=C1"),
    "na2co3": ("Na2CO3", "sodium carbonate", "[Na+].[Na+].[O-]C([O-])=O"),
    "co_salophen": ("Co(salophen)", "cobalt salophen complex", None), "o2": ("O2", "molecular oxygen", "O=O"),
    "anisole": ("anisole", "anisole", "COc1ccccc1"), "et3n": ("Et3N", "triethylamine", "CCN(CC)CC"),
    "acoh": ("AcOH", "acetic acid", "CC(=O)O"), "pdtfa2": ("Pd(TFA)2", "palladium(II) trifluoroacetate", None),
    "pdpph3cl2": ("Pd(PPh3)2Cl2", "bis(triphenylphosphine)palladium(II) dichloride", None),
    "pd_mecn_cl2": ("Pd(MeCN)2Cl2", "bis(acetonitrile)palladium(II) dichloride", None),
    "pd_acac2": ("Pd(acac)2", "palladium(II) acetylacetonate", None),
    "nah": ("NaH", "sodium hydride", "[Na+].[H-]"), "lialh4": ("LiAlH4", "lithium aluminium hydride", "[Li+].[AlH4-]"),
    "mei": ("MeI", "methyl iodide", "CI"), "bnbr": ("BnBr", "benzyl bromide", "BrCc1ccccc1"),
    "dmp": ("2,2-dimethoxypropane", "2,2-dimethoxypropane", "COC(C)(C)OC"),
    "ptsoh": ("p-TsOH·H2O", "p-toluenesulfonic acid monohydrate", None),
    "propargyl_malonate": ("propargyl malonate", "propargyl malonate substrate (ester groups generic in Scheme 5)", None),
    "bromoallene": ("bromoallene", "bromoallene", None),
    "bromoallene_d": ("1-bromo-3,4,4-trimethyl-1,2-pentadiene", "1-bromo-3,4,4-trimethyl-1,2-pentadiene", None),
    "water": ("H2O", "water", "O"), "nabhco3": ("NaHCO3", "sodium bicarbonate", "[Na+].O=C([O-])O"),
    "naso4": ("Na2SO4", "sodium sulfate", "[Na+].[Na+].[O-]S(=O)(=O)[O-]"),
    "brine": ("brine", "brine", None), "silica": ("silica gel", "silica gel", None),
    "b2pin2": ("B2pin2", "bis(pinacolato)diboron", None), "cuiprcl": ("Cu(iPr)Cl", "Scheme 3A copper complex", None),
    "pd_dpephos": ("Pd(DPEphos)Cl2", "Scheme 3A palladium complex", None),
    "nmm": ("N-methylmaleimide", "N-methylmaleimide", "CN1C(=O)C=CC1=O"),
    "mecn": ("CH3CN", "acetonitrile", "CC#N"), "dce": ("DCE", "1,2-dichloroethane", "ClCCCl"),
    "dioxane": ("dioxane", "1,4-dioxane", "O1CCOCC1"), "dmso": ("DMSO", "dimethyl sulfoxide", "CS(C)=O"),
    "thf": ("THF", "tetrahydrofuran", "C1CCOC1"), "meoh": ("MeOH", "methanol", "CO"),
    "acetone": ("acetone", "acetone", "CC(C)=O"), "toluene": ("toluene", "toluene", "Cc1ccccc1"),
    "et2o": ("Et2O", "diethyl ether", "CCOCC"), "dmf": ("DMF", "N,N-dimethylformamide", "CN(C)C=O"),
    "pentane": ("pentane", "pentane", "CCCCC"), "etoac": ("ethyl acetate", "ethyl acetate", "CCOC(C)=O"),
}
for key, (label, name, smi) in common.items():
    mol("m_" + key, label, [P[2], P[3], P[4], P[5], P[6]], name=name, smiles=smi)
mol("m_diethyl_propargyl_malonate", "diethyl propargylmalonate", [P[5], PI[5]],
    name="diethyl propargylmalonate", smiles="CCOC(=O)C(CC#C)C(=O)OCC")
mol("m_dimethyl_propargyl_malonate", "dimethyl propargylmalonate", [P[5], PI[5]],
    name="dimethyl propargylmalonate", smiles="COC(=O)C(CC#C)C(=O)OC",
    note="Identity follows the dimethyl 1d product and same-manner statement; exact charge is not repeated.")

for mid, label, ev in [
    ("m_generic_allenyne", "generic allenyne X", [P[1], P[2]]), ("m_generic_enallene", "generic enallene Y", [P[1], P[2]]),
    ("m_generic_polyene", "generic cross-conjugated polyene Z", [P[1], P[2]]),
    ("m_prior_alkyne", "alkyne in ref. 13", [P[2]]), ("m_prior_bisboronate", "bisboronate intermediate in ref. 13", [P[2]]),
    ("m_prior_allenoate", "2nd allenoate in ref. 13", [P[2]]), ("m_prior_dendralene", "[4]/[5]dendralene in ref. 13", [P[2]]),
    ("m_scheme2_dendralene", "dendralene in Scheme 2", [P[1]]), ("m_scheme2_tris", "tris adducts in Scheme 2", [P[1]]),
    ("m_scheme2_tetracycle", "tetracyclic [5.6.6.6]fenestrane", [P[1]])]:
    mol(mid, label, ev, note="Cited prior-work/generic identity; not a new compound of this paper.")


# Page 1/2 prior-work arrows and current generic transformation.
for suffix, product in (("literature_scheme2_tris_adduct", "m_scheme2_tris"), ("literature_scheme2_fenestrane", "m_scheme2_tetracycle")):
    r = rxn(suffix, [P[1], PI[1]], "literature_example", name="Cited dendralene cycloaddition", location="Scheme 2 arrow", mode="cc")
    r["participants"] = [part("d", "m_scheme2_dendralene", "reactant", r["evidence_ids"]), part("nmm", "m_nmm", "reactant", r["evidence_ids"]), part("p", product, "product", r["evidence_ids"])]
    r["notes"].append(F("ref. 12", "Cited prior work; not a current-paper experiment.", r["evidence_ids"]))

r = rxn("literature_scheme3a_borylation", [P[2], PI[2]], "literature_example", name="Borylation in ref. 13", location="Scheme 3A first arrow", mode="cb")
r["participants"] = [part("s", "m_prior_alkyne", "reactant", r["evidence_ids"]), part("b", "m_b2pin2", "reagent", r["evidence_ids"], equivalents=("2 equiv", q(2, "equiv")), reference="m_prior_alkyne"), part("cu", "m_cuiprcl", "catalyst", r["evidence_ids"]), part("p", "m_prior_bisboronate", "product", r["evidence_ids"])]
r["steps"] = [step("heat", "reaction", [x["participant_id"] for x in r["participants"]], r["evidence_ids"], "heat", ("75 °C", q(75, "°C")))]
r = rxn("literature_scheme3a_suzuki", [P[2], PI[2]], "literature_example", name="Suzuki-Miyaura reaction in ref. 13", location="Scheme 3A second arrow", mode="cc")
r["participants"] = [part("s1", "m_prior_bisboronate", "reactant", r["evidence_ids"]), part("s2", "m_prior_allenoate", "reactant", r["evidence_ids"]), part("pd", "m_pd_dpephos", "catalyst", r["evidence_ids"]), part("p", "m_prior_dendralene", "product", r["evidence_ids"])]

r = rxn("generic_current_cascade", [P[1], PI[1], P[2], PI[2]], "generic_scheme", name="Palladium-catalyzed oxidative C-C bond forming cascade", location="graphical abstract and Scheme 3B")
r["participants"] = [part("x", "m_generic_allenyne", "reactant", r["evidence_ids"]), part("y", "m_generic_enallene", "reactant", r["evidence_ids"]), part("z", "m_generic_polyene", "product", r["evidence_ids"]), part("pd", "m_pd_oac2", "catalyst", r["evidence_ids"])]


def make_method(letter):
    ev = [P[2], PI[2], P[3], PI[3], P[4], PI[4]]
    r = rxn("general_method_" + letter.lower(), ev, "generic_scheme", name="Method " + letter, location="Scheme 4/Table 2/Table 3 footnotes")
    ps = [part("s1", "m_1a", "reactant", ev), part("s2", "m_2a", "reactant", ev, equivalents=("2 equiv", q(2, "equiv")), reference="m_1a"),
          part("pd", "m_pd_oac2", "catalyst", ev, mol_percent=("5 mol%", q(5, "mol%")), reference="m_1a"),
          part("base", "m_na2co3", "additive", ev, mol_percent=("20 mol%", q(20, "mol%")), reference="m_1a"), part("solv", "m_mecn", "solvent", ev)]
    if letter == "A":
        ps.append(part("bq", "m_bq", "reagent", ev, equivalents=("2 equiv", q(2, "equiv")), reference="m_1a"))
    else:
        ps += [part("bq", "m_bq", "reagent", ev, mol_percent=("20 mol%", q(20, "mol%")), reference="m_1a"), part("co", "m_co_salophen", "catalyst", ev, mol_percent=("5 mol%", q(5, "mol%")), reference="m_1a"), part("o2", "m_o2", "reagent", ev)]
    r["participants"] = ps
    r["steps"] = [step("reaction", "reaction", [x["participant_id"] for x in ps], ev, "stir in CH3CN", ("80 °C", q(80, "°C")), ("24 h", q(24, "h")), "O2 balloon" if letter == "B" else None)]
    solvent(r, "m_mecn", ev)
    return r["reaction_id"]


METHOD_A, METHOD_B = make_method("A"), make_method("B")

scheme5 = [
    ("generic_scheme5_alkylation", ["m_propargyl_malonate", "m_bromoallene"], ["m_1a", "m_1b", "m_1c", "m_1d"], "NaH alkylation to 1a-d", "cc"),
    ("generic_scheme5_reduction", ["m_1a", "m_1c", "m_1d", "m_lialh4"], ["m_1ax", "m_1cx", "m_1dx"], "LiAlH4 reduction to diols", "ch"),
    ("generic_scheme5_etherification", ["m_1ax", "m_1cx", "m_1dx", "m_nah"], ["m_1e", "m_1f", "m_1g", "m_1j"], "NaH/RX etherification", "co"),
    ("generic_scheme5_acetalization", ["m_1ax", "m_1dx", "m_dmp", "m_ptsoh"], ["m_1h", "m_1i"], "acetonide formation", "co2")]
for suffix, inputs, outputs, name, mode in scheme5:
    r = rxn(suffix, [P[2], PI[2]], "generic_scheme", name=name, location="Scheme 5 arrow", mode=mode)
    r["participants"] = [part("in" + str(i), mid, "reactant" if i < 3 else "reagent", r["evidence_ids"]) for i, mid in enumerate(inputs)] + [part("out" + str(i), mid, "product", r["evidence_ids"]) for i, mid in enumerate(outputs)]


# Table 1: 18 independent optimization experiments.
table1 = [
    (1, "m_pd_oac2", "m_dce", 64, 6, None), (2, "m_pd_oac2", "m_dioxane", 51, 6, None),
    (3, "m_pd_oac2", "m_dmso", 55, None, None), (4, "m_pd_oac2", "m_thf", 43, None, None),
    (5, "m_pd_oac2", "m_meoh", 50, 4, None), (6, "m_pd_oac2", "m_acetone", 33, None, None),
    (7, "m_pd_oac2", "m_toluene", 65, 7, None), (8, "m_pd_oac2", "m_mecn", 75, None, None),
    (9, "m_pdtfa2", "m_mecn", None, None, (94, 91)), (10, "m_pdpph3cl2", "m_mecn", None, None, (70, 71)),
    (11, "m_pd_mecn_cl2", "m_mecn", None, None, (71, 47)), (12, "m_pd_acac2", "m_mecn", 60, 6, None),
    (13, "m_pd_oac2", "m_mecn", 80, None, None), (14, "m_pd_oac2", "m_mecn", 83, None, None),
    (15, "m_pd_oac2", "m_mecn", 88, None, None), (16, "m_pd_oac2", "m_mecn", 94, None, None),
    (17, "m_pd_oac2", "m_mecn", 91, None, None), (18, "m_pd_oac2", "m_mecn", 90, None, None)]
for entry, catalyst, solv, y3, y5, rec in table1:
    ev = [P[2], PI[2], P[3]]
    r = rxn(f"table1_entry_{entry:02d}", ev, categories=["optimization"], name="Optimization of oxidative cascade",
            location=f"Table 1 entry {entry}", mode="none" if entry in {9, 10, 11} else "cascade")
    r["general_scheme_ids"] = [METHOD_A]
    e2 = 1.5 if entry == 13 else 2.0 if entry >= 14 else 1.1
    pdload = 1 if entry == 14 else 5
    ps = [part("p1a", "m_1a", "reactant", ev, amount=("0.1 mmol", q(0.1, "mmol")), equivalents=("1.0 equiv", q(1, "equiv")), reference="m_1a"),
          part("p2a", "m_2a", "reactant", ev, amount=("0.11 mmol", q(0.11, "mmol")) if entry <= 12 else None, equivalents=(f"{e2} equiv", q(e2, "equiv")), reference="m_1a"),
          part("pd", catalyst, "catalyst", ev, mol_percent=(f"{pdload} mol%", q(pdload, "mol%")), reference="m_1a"),
          part("bq", "m_bq", "reagent", ev, equivalents=(f"{e2} equiv", q(e2, "equiv")), reference="m_1a"),
          part("solv", solv, "solvent", ev, volume=("1 mL", q(1, "mL")))]
    if entry in {15, 16}:
        ps.append(part("base", "m_na2co3", "additive", ev, mol_percent=("20 mol%", q(20, "mol%")), reference="m_1a"))
    if entry == 17:
        ps.append(part("add", "m_et3n", "additive", ev, mol_percent=("20 mol%", q(20, "mol%")), reference="m_1a"))
    if entry == 18:
        ps.append(part("add", "m_acoh", "additive", ev, mol_percent=("20 mol%", q(20, "mol%")), reference="m_1a"))
    if y3 is not None:
        ps.append(part("prod3", "m_3a", "product", ev))
    if y5 is not None:
        ps.append(part("prod5", "m_5", "product", ev))
    r["participants"] = ps
    temp = 60 if entry in {4, 5, 6} else 80
    r["steps"] = [step("reaction", "reaction", [x["participant_id"] for x in ps], ev, "stir in indicated solvent",
                       (f"{temp} °C", q(temp, "°C")), ("24 h", q(24, "h")) if entry in {16, 17, 18} else None)]
    solvent(r, solv, ev, ("1 mL", q(1, "mL")))
    r["reaction_scale"] = F("1a, 0.1 mmol", {"reference_molecule_id": "m_1a", "amount": q(0.1, "mmol")}, ev)
    if y3 is not None:
        r["outcomes"].append(out("yield3a", "yield", ["m_3a"], ev, f"{y3}%", q(y3, "%"), "NMR yield", "1H NMR with anisole internal standard"))
    if y5 is not None:
        r["outcomes"].append(out("yield5", "yield", ["m_5"], ev, f"{y5}%", q(y5, "%"), "NMR yield", "1H NMR with anisole internal standard"))
    if rec:
        r["outcomes"] += [out("recovery1a", "recovery", ["m_1a"], ev, f"{rec[0]}%", q(rec[0], "%"), basis="1H NMR with anisole internal standard"),
                          out("recovery2a", "recovery", ["m_2a"], ev, f"{rec[1]}%", q(rec[1], "%"), basis="1H NMR with anisole internal standard"),
                          out("no_rxn", "no_reaction", [], ev, "not suitable; starting materials recovered", True, basis="1H NMR with anisole internal standard")]
    r["notes"].append(F("Table dashes were not converted to zero.", "Dash/blank outcomes remain missing, not zero.", ev))


def scope_record(substrate, enallene, product, method, yield_pct, table, extra_ev=()):
    ev = ([P[3], PI[3]] if table == 2 else [P[4], PI[4]]) + list(extra_ev)
    r = rxn(f"table{table}_{product[2:]}_method_{method.lower()}", ev, categories=["substrate_scope"],
            name="Oxidative coupling-carbocyclization", location=f"Table {table}, {product[2:]}, Method {method}")
    r["general_scheme_ids"] = [METHOD_A if method == "A" else METHOD_B]
    ps = [part("substrate", substrate, "reactant", ev, equivalents=("1.0 equiv", q(1, "equiv")), reference=substrate),
          part("enallene", enallene, "reactant", ev, equivalents=("2.0 equiv", q(2, "equiv")), reference=substrate),
          part("product", product, "product", ev),
          part("pd", "m_pd_oac2", "catalyst", ev, mol_percent=("5 mol%", q(5, "mol%")), reference=substrate),
          part("base", "m_na2co3", "additive", ev, mol_percent=("20 mol%", q(20, "mol%")), reference=substrate),
          part("solv", "m_mecn", "solvent", ev)]
    if method == "A":
        ps.append(part("bq", "m_bq", "reagent", ev, equivalents=("2 equiv", q(2, "equiv")), reference=substrate))
    else:
        ps += [part("bq", "m_bq", "reagent", ev, mol_percent=("20 mol%", q(20, "mol%")), reference=substrate),
               part("co", "m_co_salophen", "catalyst", ev, mol_percent=("5 mol%", q(5, "mol%")), reference=substrate), part("o2", "m_o2", "reagent", ev)]
    r["participants"] = ps
    r["reaction_scale"] = F("substrate 1, 1.0 equiv", {"reference_molecule_id": substrate, "amount": None, "basis": "1.0 equiv"}, ev)
    r["steps"] = [step("reaction", "reaction", [x["participant_id"] for x in ps], ev, "stir in CH3CN", ("80 °C", q(80, "°C")), ("24 h", q(24, "h")), "O2 balloon" if method == "B" else None),
                  step("purification", "purification", [], ev, "column chromatography")]
    solvent(r, "m_mecn", ev)
    r["outcomes"] = [out("yield", "yield", [product], ev, f"{yield_pct}%", q(yield_pct, "%"), "isolated yield", "isolated yield after column chromatography")]
    return r


for label, y in {"3a": 90, "3b": 89, "3c": 83, "3d": 82, "3e": 80, "3f": 77, "3g": 81, "3h": 75, "3i": 34, "3j": 83}.items():
    scope_record("m_1" + label[-1], "m_2a", "m_" + label, "A", y, 2, [P[2], PI[2], P[6], PI[6]] if label == "3a" else [])
for label, y in {"3a": 85, "3d": 78}.items():
    scope_record("m_1" + label[-1], "m_2a", "m_" + label, "B", y, 2, [P[2], PI[2], P[6], PI[6]] if label == "3a" else [])
for label, y in {"3k": 85, "3l": 72, "3m": 78, "3n": 75, "3o": 61, "3p": 76, "3q": 57, "3r": 71, "3s": 80, "3t": 80, "3u": 60, "3v": 63}.items():
    scope_record("m_1a", "m_2" + label[-1], "m_" + label, "A", y, 3)
for label, y in {"3k": 87, "3m": 80, "3n": 77}.items():
    scope_record("m_1a", "m_2" + label[-1], "m_" + label, "B", y, 3)

byid = {r["reaction_id"]: r for r in bundle["reactions"]}
u = byid[f"{DOC}:table3_3u_method_a"]
u["outcomes"].append(out("geometry", "selectivity", ["m_3u"], u["evidence_ids"], "E/Z", "E/Z mixture; ratio unreadable/unreported",
                         scope={"metric": "E/Z", "order": "E:Z", "ratio": None}))
u["review"]["unresolved"].append("Table 3 marks 3u as E/Z but no numerical ratio is readable in the supplied image.")


# Page-6 detailed 3a procedures enrich the already-created Table 2 records.
def enrich_3a(method):
    r = byid[f"{DOC}:table2_3a_method_{method.lower()}"]
    charges = {
        "substrate": ((47.2, "mg"), (0.2, "mmol")), "enallene": ((77, "mg"), (0.4, "mmol")),
        "pd": ((2.2, "mg"), (0.01, "mmol")), "base": ((4.2, "mg"), (0.04, "mmol")),
        "bq": ((43.2 if method == "A" else 4.4, "mg"), (0.4 if method == "A" else 0.04, "mmol")),
        "co": ((3.7, "mg"), (0.01, "mmol"))}
    for p in r["participants"]:
        if p["participant_id"] in charges:
            a, b = charges[p["participant_id"]]
            p["charges"]["mass"] = F(f"{a[0]} {a[1]}", q(*a), [P[6]])
            p["charges"]["amount"] = F(f"{b[0]} {b[1]}", q(*b), [P[6]])
        if p["participant_id"] == "solv":
            p["charges"]["volume"] = F("2.0 mL", q(2.0, "mL"), [P[6]])
    r["solvent_conditions"]["total_volume"] = F("2.0 mL", q(2.0, "mL"), [P[6]])
    r["steps"] += [step("concentrate", "workup", [], [P[6]], "concentrate in vacuo"),
                   step("short_column", "purification", [], [P[6]], "short silica-gel column; ethyl acetate/pentane 5:95")]
    mass = 77 if method == "A" else 73
    r["outcomes"].append(out("product_mass", "product_amount", ["m_3a"], [P[6]], f"{mass} mg", q(mass, "mg"), basis="after short column chromatography"))
    r["procedure"]["original_text"] = F(f"Representative Procedure: Synthesis of 3a (Method {method})", f"Detailed Method {method} procedure", [P[6]])


enrich_3a("A")
enrich_3a("B")

r = rxn("scaleup_3a_method_a", [P[6], PI[6]], categories=["scale_up"], name="Method A scale-up synthesis of 3a", location="page 6 scale-up")
r["general_scheme_ids"] = [METHOD_A]
r["participants"] = [part("s", "m_1a", "reactant", r["evidence_ids"], mass=("0.307 g", q(0.307, "g")), amount=("1.3 mmol", q(1.3, "mmol"))),
                     part("e", "m_2a", "reactant", r["evidence_ids"]), part("p", "m_3a", "product", r["evidence_ids"])]
r["reaction_scale"] = F("1a, 1.3 mmol", {"reference_molecule_id": "m_1a", "amount": q(1.3, "mmol")}, r["evidence_ids"])
r["steps"] = [step("method_a", "reaction", ["s", "e", "p"], r["evidence_ids"], "Method A; scaled absolute coproportional charges not stated"),
              step("column", "purification", [], r["evidence_ids"], "short silica-gel column; ethyl acetate/pentane 5:95")]
r["outcomes"] = [out("mass", "product_amount", ["m_3a"], r["evidence_ids"], "0.45 g", q(0.45, "g"), basis="after short column chromatography"),
                 out("yield", "yield", ["m_3a"], r["evidence_ids"], "81%", q(81, "%"), "isolated yield", "after short column chromatography")]
r["review"]["unresolved"].append("Absolute scaled charges other than 1a/product are not stated.")


# Scheme 6: the 2a positive row is duplicate evidence for 3a Method A; 4a/4b are separate controls.
for label in ("4a", "4b"):
    ev = [P[3], PI[3]]
    r = rxn("scheme6_control_" + label, ev, categories=["control"], name="Comparative experiment without pendant olefin",
            location="Scheme 6 " + label, mode="none")
    r["general_scheme_ids"] = [METHOD_A]
    r["participants"] = [part("s1", "m_1a", "reactant", ev, equivalents=("1.0 equiv", q(1, "equiv")), reference="m_1a"),
                         part("s2", "m_" + label, "reactant", ev, equivalents=("2.0 equiv", q(2, "equiv")), reference="m_1a"),
                         part("pd", "m_pd_oac2", "catalyst", ev, mol_percent=("5 mol%", q(5, "mol%")), reference="m_1a"),
                         part("bq", "m_bq", "reagent", ev, equivalents=("2.0 equiv", q(2, "equiv")), reference="m_1a"),
                         part("base", "m_na2co3", "additive", ev, equivalents=("0.2 equiv", q(0.2, "equiv")), reference="m_1a"),
                         part("solv", "m_mecn", "solvent", ev)]
    r["steps"] = [step("reaction", "reaction", [x["participant_id"] for x in r["participants"]], ev, "Method A comparative experiment", ("80 °C", q(80, "°C")), ("24 h", q(24, "h")))]
    solvent(r, "m_mecn", ev)
    r["outcomes"] = [out("no_3", "no_reaction", [], ev, "no detectable amounts of 3; complex mixture", True)]
    r["notes"].append(F("Complex mixture", "Complex mixture; no product 3 detected.", ev))


# Equations 1-3 are three actual mechanistic experiments.
ev = [P[4], PI[4]]
r = rxn("kie_eq1_competition", ev, categories=["mechanistic_study"], name="Intermolecular competition KIE experiment", location="equation 1")
r["participants"] = [part("a", "m_1a", "reactant", ev, equivalents=("1.0 equiv", q(1, "equiv")), reference="m_1a"),
                     part("h", "m_2a", "reactant", ev, equivalents=("1.0 equiv", q(1, "equiv")), reference="m_1a"),
                     part("d", "m_2a_d6", "reactant", ev, equivalents=("1.0 equiv", q(1, "equiv")), reference="m_1a"),
                     part("ph", "m_3a", "product", ev), part("pd5", "m_3a_d5", "product", ev),
                     part("cat", "m_pd_oac2", "catalyst", ev, mol_percent=("5 mol%", q(5, "mol%")), reference="m_1a"),
                     part("bq", "m_bq", "reagent", ev, equivalents=("2 equiv", q(2, "equiv")), reference="m_1a"),
                     part("base", "m_na2co3", "additive", ev, equivalents=("0.2 equiv", q(0.2, "equiv")), reference="m_1a"), part("solv", "m_mecn", "solvent", ev)]
r["steps"] = [step("reaction", "reaction", [x["participant_id"] for x in r["participants"]], ev, "intermolecular competition", ("75 °C", q(75, "°C")), ("1 h", q(1, "h")))]
r["outcomes"] = [out("ratio", "selectivity", ["m_3a", "m_3a_d5"], ev, "3a/3a-d5 = 1.8:1", "1.8:1", scope={"metric": "product ratio", "order": ["m_3a", "m_3a_d5"]}),
                 out("kie", "other", [], ev, "kH/kD = 3.5", 3.5, scope="competitive KIE"),
                 out("conversion", "conversion", ["m_1a"], ev, "Conv = 72%", q(72, "%"))]
r["review"]["unresolved"].append("The figure does not explicitly assign Conv = 72% to a named substrate; provisionally targeted to limiting 1a.")

for number, enallene, product in ((2, "m_2a", "m_3a"), (3, "m_2a_d6", "m_3a_d5")):
    r = rxn(f"kie_eq{number}_parallel", ev, categories=["mechanistic_study"], name="Parallel initial-rate KIE experiment", location=f"equation {number}")
    r["participants"] = [part("a", "m_1a", "reactant", ev, equivalents=("0.5 equiv", q(0.5, "equiv")), reference=enallene),
                         part("e", enallene, "reactant", ev, equivalents=("1.0 equiv", q(1, "equiv")), reference=enallene), part("p", product, "product", ev),
                         part("cat", "m_pd_oac2", "catalyst", ev, mol_percent=("5 mol%", q(5, "mol%")), reference=enallene),
                         part("bq", "m_bq", "reagent", ev, equivalents=("1 equiv", q(1, "equiv")), reference=enallene),
                         part("base", "m_na2co3", "additive", ev, equivalents=("0.2 equiv", q(0.2, "equiv")), reference=enallene), part("solv", "m_mecn", "solvent", ev)]
    r["steps"] = [step("reaction", "reaction", [x["participant_id"] for x in r["participants"]], ev, "parallel initial-rate experiment", ("75 °C", q(75, "°C")), ("30 min", q(30, "min")))]
    r["outcomes"] = [out("yield", "yield", [product], ev, "<33% yield", {"number": 33, "unit": "%", "qualifier": "lt", "raw_text": "<33%"}),
                     out("kie", "other", [], ev, "kH/kD = 3.4 from eqs 2 and 3", 3.4, scope="paired initial-rate KIE")]
    r["review"]["unresolved"].append("Measurement basis for <33% yield and full initial-rate calculation are unavailable/not reported on this page.")

r = rxn("scheme7_proposed_mechanism", [P[4], PI[4], P[5], PI[5]], "proposed_mechanism", name="Proposed mechanism for formation of 3", location="Scheme 7")
r["participants"] = [part("one", "m_1a", "reactant", r["evidence_ids"]), part("two", "m_2a", "reactant", r["evidence_ids"]), part("three", "m_3a", "product", r["evidence_ids"])] + [part(x, "m_" + x, "unknown", r["evidence_ids"]) for x in ("int_1", "int_2", "int_3", "int_4")]
r["notes"].append(F("C-H cleavage, alkyne insertion, allene insertion/carbocyclization, beta-hydride elimination", "Author-proposed sequence; intermediates were not isolated.", r["evidence_ids"]))


def prep(suffix, start, product, ev, name, mode, mass, y, ref=None, purification=None):
    r = rxn("prep_" + suffix, ev, categories=["other"], name=name, location=name, mode=mode)
    r["participants"] = [part("start", start, "reactant", ev), part("product", product, "product", ev)]
    r["outcomes"] = [out("mass", "product_amount", [product], ev, mass[0], mass[1], basis="after column chromatography"),
                     out("yield", "yield", [product], ev, f"{y}%", q(y, "%"), "isolated yield", "after column chromatography")]
    if ref:
        r["procedure"]["general_procedure_reference"] = F(ref, ref, ev)
        r["review"]["unresolved"].append("Exact charges are not repeated in the same-manner paragraph.")
    if purification:
        r["steps"].append(step("column", "purification", [], ev, purification))
    return r


ev5 = [P[5], PI[5]]
r = prep("1b", "m_diethyl_propargyl_malonate", "m_1b", ev5, "Synthesis of 1b", "cc", ("0.71 g", q(0.71, "g")), 31)
r["participants"] += [part("nah", "m_nah", "reagent", ev5, mass=("0.456 g", q(0.456, "g")), amount=("11.4 mmol", q(11.4, "mmol"))),
                      part("bromo", "m_bromoallene", "reactant", ev5, mass=("2.6 g", q(2.6, "g")), amount=("17.6 mmol", q(17.6, "mmol"))),
                      part("thf", "m_thf", "solvent", ev5),
                      part("water", "m_water", "workup_agent", ev5, "workup", volume=("10 mL", q(10, "mL"))),
                      part("ether", "m_et2o", "workup_agent", ev5, "workup"), part("dry", "m_naso4", "workup_agent", ev5, "workup"),
                      part("silica", "m_silica", "purification_agent", ev5, "purification"), part("pentane", "m_pentane", "purification_agent", ev5, "purification"), part("etoac", "m_etoac", "purification_agent", ev5, "purification")]
r["participants"][0]["charges"]["mass"] = F("2.0 g", q(2.0, "g"), ev5)
r["participants"][0]["charges"]["amount"] = F("8.83 mmol", q(8.83, "mmol"), ev5)
r["reaction_scale"] = F("diethyl propargylmalonate, 8.83 mmol", {"reference_molecule_id": "m_diethyl_propargyl_malonate", "amount": q(8.83, "mmol")}, ev5)
r["steps"] = [step("add", "reaction", ["start", "nah", "thf"], ev5, "add malonate solution to NaH at 0 °C", ("0 °C", q(0, "°C"))),
              step("hold", "reaction", ["start", "nah"], ev5, "stir after addition", ("room temperature", "room temperature"), ("20 min", q(20, "min"))),
              step("bromo", "reaction", ["bromo"], ev5, "add bromoallene solution", ("room temperature", "room temperature")),
              step("reflux", "reaction", ["start", "nah", "bromo"], ev5, "reflux", ("reflux", "reflux"), ("20 h", q(20, "h"))),
              step("workup", "workup", ["water", "ether", "dry"], ev5, "cool, evaporate most solvent, dilute, quench, extract, dry and evaporate"),
              step("column", "purification", ["silica", "pentane", "etoac"], ev5, "silica gel column; pentane/ethyl acetate 30:1")]

r = prep("1d", "m_dimethyl_propargyl_malonate", "m_1d", ev5, "Synthesis of 1d", "cc", ("0.79 g", q(0.79, "g")), 32, "prepared in the same manner as 1b", "silica gel column; pentane/ethyl acetate 30:1")
r["participants"].append(part("bromo", "m_bromoallene_d", "reactant", ev5))
r["steps"].insert(0, step("same", "reaction", ["start", "bromo"], ev5, "prepared in same manner as 1b; exact charges not repeated"))

r = prep("1ax", "m_1a", "m_1ax", ev5, "Synthesis of 1ax", "ch", ("203 mg", q(203, "mg")), 79)
r["participants"] += [part("lah", "m_lialh4", "reagent", ev5, mass=("171 mg", q(171, "mg")), amount=("4.5 mmol", q(4.5, "mmol"))),
                      part("ether", "m_et2o", "solvent", ev5), part("water", "m_water", "workup_agent", ev5, "workup", volume=("5 mL", q(5, "mL"))),
                      part("dry", "m_naso4", "workup_agent", ev5, "workup"), part("silica", "m_silica", "purification_agent", ev5, "purification")]
r["participants"][0]["charges"]["mass"] = F("0.34 g", q(0.34, "g"), ev5)
r["participants"][0]["charges"]["amount"] = F("1.5 mmol", q(1.5, "mmol"), ev5)
r["reaction_scale"] = F("1a, 1.5 mmol", {"reference_molecule_id": "m_1a", "amount": q(1.5, "mmol")}, ev5)
r["steps"] = [step("add", "reaction", ["start", "lah", "ether"], ev5, "add 1a solution to LiAlH4 at 0 °C", ("0 °C", q(0, "°C"))),
              step("hold", "reaction", ["start", "lah"], ev5, "stir", ("room temperature", "room temperature"), ("2 h", q(2, "h"))),
              step("workup", "workup", ["water", "ether", "dry"], ev5, "carefully quench, extract, wash, dry and evaporate"),
              step("column", "purification", ["silica"], ev5, "silica gel column; pentane/ethyl acetate 1:3")]
for label, start, mass, y, ratio in (("1cx", "m_1c", 268, 81, "1:2.5"), ("1dx", "m_1d", 323, 97, "1:3")):
    r = prep(label, start, "m_" + label, ev5, "Synthesis of " + label, "ch", (f"{mass} mg", q(mass, "mg")), y, "prepared in the same manner as 1ax", f"silica gel column; pentane/ethyl acetate {ratio}")
    r["participants"].append(part("lah", "m_lialh4", "reagent", ev5))
    r["steps"].insert(0, step("same", "reaction", ["start", "lah"], ev5, "prepared in same manner as 1ax; exact charges not repeated"))

ev56 = [P[5], PI[5], P[6], PI[6]]
r = prep("1e", "m_1ax", "m_1e", ev56, "Synthesis of 1e", "co", ("150 mg", q(150, "mg")), 81)
r["participants"] += [part("nah", "m_nah", "reagent", ev56, mass=("142 mg", q(142, "mg")), amount=("3.6 mmol", q(3.6, "mmol"))),
                      part("mei", "m_mei", "reagent", ev56, volume=("0.45 mL", q(0.45, "mL")), amount=("7.1 mmol", q(7.1, "mmol"))),
                      part("thf", "m_thf", "solvent", ev56, volume=("10 mL", q(10, "mL"))),
                      part("water", "m_water", "workup_agent", ev56, "workup"), part("ether", "m_et2o", "workup_agent", ev56, "workup"),
                      part("brine", "m_brine", "workup_agent", ev56, "workup", volume=("20 mL", q(20, "mL"))), part("dry", "m_naso4", "workup_agent", ev56, "workup"), part("silica", "m_silica", "purification_agent", ev56, "purification")]
r["participants"][0]["charges"]["mass"] = F("160 mg", q(160, "mg"), ev56)
r["participants"][0]["charges"]["amount"] = F("0.89 mmol", q(0.89, "mmol"), ev56)
r["reaction_scale"] = F("1ax, 0.89 mmol", {"reference_molecule_id": "m_1ax", "amount": q(0.89, "mmol")}, ev56)
r["steps"] = [step("deprotonate", "reaction", ["start", "nah", "thf"], ev56, "add NaH and stir", ("0 °C", q(0, "°C")), ("30 min", q(30, "min"))),
              step("alkylate", "reaction", ["start", "mei"], ev56, "add MeI, warm to room temperature and stir", ("room temperature", "room temperature"), ("1 h", q(1, "h"))),
              step("workup", "workup", ["water", "ether", "brine", "dry"], ev56, "quench, extract, wash, dry and evaporate"),
              step("column", "purification", ["silica"], ev56, "flash silica gel column; pentane/ethyl acetate 98:2")]

ev6 = [P[6], PI[6]]
for label, start, mass, y in (("1f", "m_1cx", 209, 95), ("1g", "m_1dx", 189, 85)):
    r = prep(label, start, "m_" + label, ev6, "Synthesis of " + label, "co", (f"{mass} mg", q(mass, "mg")), y, "prepared in the same manner as 1e", "flash silica gel column; pentane/ethyl acetate 98:2")
    r["participants"] += [part("nah", "m_nah", "reagent", ev6), part("mei", "m_mei", "reagent", ev6)]
    r["steps"].insert(0, step("same", "reaction", ["start", "nah", "mei"], ev6, "prepared in same manner as 1e; exact charges not repeated"))

r = prep("1h", "m_1ax", "m_1h", ev6, "Synthesis of 1h", "co2", ("111 mg", q(111, "mg")), 84)
r["participants"] += [part("dmp", "m_dmp", "reagent", ev6, mass=("72 mg", q(72, "mg")), amount=("0.69 mmol", q(0.69, "mmol"))),
                      part("acid", "m_ptsoh", "catalyst", ev6, mass=("5.7 mg", q(5.7, "mg")), amount=("0.03 mmol", q(0.03, "mmol"))),
                      part("acetone", "m_acetone", "solvent", ev6, volume=("5 mL", q(5, "mL"))),
                      part("bicarb", "m_nabhco3", "workup_agent", ev6, "workup"), part("ether", "m_et2o", "workup_agent", ev6, "workup"),
                      part("brine", "m_brine", "workup_agent", ev6, "workup"), part("dry", "m_naso4", "workup_agent", ev6, "workup")]
r["participants"][0]["charges"]["mass"] = F("108 mg", q(108, "mg"), ev6)
r["participants"][0]["charges"]["amount"] = F("0.60 mmol", q(0.60, "mmol"), ev6)
r["steps"] = [step("reaction", "reaction", ["start", "dmp", "acid", "acetone"], ev6, "stir", ("room temperature", "room temperature"), ("3 h", q(3, "h"))),
              step("workup", "workup", ["bicarb", "ether", "brine", "dry"], ev6, "quench with saturated NaHCO3, evaporate acetone, extract, wash, dry and evaporate"),
              step("column", "purification", [], ev6, "flash silica gel column; pentane/ethyl acetate 97:3")]

r = prep("1i", "m_1dx", "m_1i", ev6, "Synthesis of 1i", "co2", ("128 mg", q(128, "mg")), 81, "prepared in the same manner as 1h", "flash silica gel column; pentane/ethyl acetate 97:3")
r["participants"] += [part("dmp", "m_dmp", "reagent", ev6), part("acid", "m_ptsoh", "catalyst", ev6)]
r["steps"].insert(0, step("same", "reaction", ["start", "dmp", "acid"], ev6, "prepared in same manner as 1h; exact charges not repeated"))

r = prep("1j", "m_1ax", "m_1j", ev6, "Synthesis of 1j", "co", ("221 mg", q(221, "mg")), 79)
r["participants"] += [part("nah", "m_nah", "reagent", ev6, mass=("124 mg", q(124, "mg")), amount=("3.1 mmol", q(3.1, "mmol"))),
                      part("bnbr", "m_bnbr", "reagent", ev6, volume=("170 uL", q(170, "uL")), amount=("1.4 mmol", q(1.4, "mmol"))),
                      part("dmf", "m_dmf", "solvent", ev6), part("ether", "m_et2o", "workup_agent", ev6, "workup"),
                      part("water", "m_water", "workup_agent", ev6, "workup"), part("dry", "m_naso4", "workup_agent", ev6, "workup")]
r["participants"][0]["charges"]["mass"] = F("140 mg", q(140, "mg"), ev6)
r["participants"][0]["charges"]["amount"] = F("0.78 mmol", q(0.78, "mmol"), ev6)
r["steps"] = [step("add", "reaction", ["start", "nah", "dmf"], ev6, "add 1ax solution to NaH", ("0 °C", q(0, "°C"))),
              step("hold", "reaction", ["start", "nah"], ev6, "stir", ("0 °C", q(0, "°C")), ("1 h", q(1, "h"))),
              step("benzylate", "reaction", ["bnbr"], ev6, "add BnBr at 0 °C, then stir at room temperature", ("room temperature", "room temperature"), ("overnight", "overnight")),
              step("workup", "workup", ["ether", "water", "dry"], ev6, "dilute, quench, extract, wash, dry and evaporate"),
              step("column", "purification", [], ev6, "silica gel column; pentane/diethyl ether 100:1")]


bundle["source_context"] = [
    {"context_id": f"{DOC}:page1_scheme1", "kind": "nonreaction_context", "evidence_ids": [P[1], PI[1]],
     "description": "Scheme 1 is a panel of example polyene structures without reaction arrows."},
    {"context_id": f"{DOC}:page5_general_information", "kind": "general_information", "evidence_ids": [P[5], PI[5]],
     "description": "General solvent, TLC, chromatography, NMR and HRMS information; not an individual synthesis."},
    {"context_id": f"{DOC}:pages5_6_characterization", "kind": "characterization", "evidence_ids": [P[5], PI[5], P[6], PI[6]],
     "description": "NMR/HRMS blocks for starting materials and products are characterization context, not separate reaction records."},
    {"context_id": f"{DOC}:si_notice", "kind": "supplement_notice", "evidence_ids": [P[2], P[4]],
     "description": "The article refers to Supporting Information, but no SI file was supplied or explicitly associated."},
]

records = [r["reaction_id"] for r in bundle["reactions"]]


def match(text):
    return [rid for rid in records if text in rid]


coverage = {"document_id": DOC, "semantic_accuracy_verified": False, "pages": [
    {"page": 1, "state": "attempted", "regions": [
        {"region": "title/abstract and graphical abstract", "finding": "True title and generic current-work oxidative allenyne/enallene cascade; merged with Scheme 3B.", "evidence_ids": [P[1], PI[1]]},
        {"region": "Scheme 1", "finding": "Polyene examples without reaction arrows; retained as nonreaction context.", "evidence_ids": [P[1], PI[1]]},
        {"region": "Scheme 2", "finding": "Two arrows from cited ref. 12; recorded as literature examples, not current experiments.", "evidence_ids": [P[1], PI[1]]}],
     "record_ids": [f"{DOC}:generic_current_cascade", f"{DOC}:literature_scheme2_tris_adduct", f"{DOC}:literature_scheme2_fenestrane"],
     "unresolved": ["Generic graphical structures were not converted to source-verified SMILES."]},
    {"page": 2, "state": "attempted", "regions": [
        {"region": "Scheme 3A", "finding": "Two cited prior-work arrows recorded separately.", "evidence_ids": [P[2], PI[2]]},
        {"region": "Scheme 3B/Scheme 4", "finding": "Current generic transformation and Method A/B general conditions.", "evidence_ids": [P[2], PI[2]]},
        {"region": "Scheme 5", "finding": "Four distinct generic starting-material preparation arrows.", "evidence_ids": [P[2], PI[2]]},
        {"region": "Table 1", "finding": "All 18 optimization rows recorded independently; dashes not converted to zero.", "evidence_ids": [P[2], PI[2]]}],
     "record_ids": [f"{DOC}:literature_scheme3a_borylation", f"{DOC}:literature_scheme3a_suzuki", f"{DOC}:generic_current_cascade", METHOD_A, METHOD_B] + match("generic_scheme5_") + match("table1_entry_"),
     "unresolved": ["Table 1 entries 1-15 have no explicit reaction time except where later footnotes state 24 h."]},
    {"page": 3, "state": "attempted", "regions": [
        {"region": "Table 1 discussion", "finding": "Narrative cross-check; no duplicate optimization experiments created.", "evidence_ids": [P[3], PI[3]]},
        {"region": "Scheme 6", "finding": "2a positive row is duplicate evidence for Table 2 3a Method A; 4a/4b are separate negative controls.", "evidence_ids": [P[3], PI[3]]},
        {"region": "Table 2", "finding": "Ten Method A and two distinct Method B scope experiments recorded.", "evidence_ids": [P[3], PI[3]]}],
     "record_ids": match("table1_entry_") + match("scheme6_control_") + match("table2_"),
     "unresolved": ["Complete source-verified product SMILES and atom mapping were not established from compact drawings."]},
    {"page": 4, "state": "attempted", "regions": [
        {"region": "Table 3", "finding": "Twelve Method A and three distinct Method B scope experiments; 3u marked E/Z without readable ratio.", "evidence_ids": [P[4], PI[4]]},
        {"region": "Equations 1-3", "finding": "Three mechanistic KIE experiments: one competition and two parallel initial-rate reactions.", "evidence_ids": [P[4], PI[4]]},
        {"region": "Scheme 7", "finding": "Author-proposed mechanism separated from experimental KIE records.", "evidence_ids": [P[4], PI[4]]}],
     "record_ids": match("table3_") + match("kie_eq") + [f"{DOC}:scheme7_proposed_mechanism"],
     "unresolved": ["3u E/Z ratio unreadable/unreported; full initial-rate details reside in unavailable SI."]},
    {"page": 5, "state": "attempted", "regions": [
        {"region": "Scheme 7 continuation/conclusion", "finding": "Mechanism prose merged with proposed-mechanism record; conclusion is narrative.", "evidence_ids": [P[5], PI[5]]},
        {"region": "general information", "finding": "Analytical/purification conventions retained as source context.", "evidence_ids": [P[5], PI[5]]},
        {"region": "preparations 1b, 1d, 1ax, 1cx, 1dx and start of 1e", "finding": "Concrete preparations recorded individually; same-manner exact charges remain missing.", "evidence_ids": [P[5], PI[5]]},
        {"region": "NMR/HRMS blocks", "finding": "Characterization context, not synthesis records.", "evidence_ids": [P[5], PI[5]]}],
     "record_ids": [f"{DOC}:scheme7_proposed_mechanism"] + [f"{DOC}:prep_{x}" for x in ("1b", "1d", "1ax", "1cx", "1dx", "1e")],
     "unresolved": ["Exact charges for same-manner preparations are not repeated."]},
    {"page": 6, "state": "attempted", "regions": [
        {"region": "completion of 1e and preparations 1f-1j", "finding": "1e cross-page procedure merged; 1f/1g/1i same-manner and detailed 1h/1j records created.", "evidence_ids": [P[6], PI[6]]},
        {"region": "Method A/B representative synthesis of 3a", "finding": "Detailed procedures merged with Table 2 3a A/B; 1.3 mmol run is a separate scale-up.", "evidence_ids": [P[6], PI[6]]},
        {"region": "3a-3c characterization", "finding": "Spectral/HRMS blocks retained as characterization context.", "evidence_ids": [P[6], PI[6]]}],
     "record_ids": [f"{DOC}:prep_{x}" for x in ("1e", "1f", "1g", "1h", "1i", "1j")] + [f"{DOC}:table2_3a_method_a", f"{DOC}:table2_3a_method_b", f"{DOC}:scaleup_3a_method_a"],
     "unresolved": ["Scale-up sentence omits absolute charges other than 1a and product."]},
]}

diagnostics = {}
for label, smi in core_smiles.items():
    try:
        diagnostics[label] = st.stereo_inspect(smi)
    except ValueError as exc:
        diagnostics[label] = {"error": str(exc), "source_verified": False}

inventory = {"document_id": DOC, "bounded_pages": [1, 2, 3, 4, 5, 6], "semantic_accuracy_verified": False,
             "title": "Synthesis of Cross-Conjugated Polyenes via Palladium-Catalyzed Oxidative C-C Bond Forming Cascade Reactions of Allenes",
             "reaction_count": len(bundle["reactions"]), "molecule_count": len(bundle["molecules"]), "record_ids": records,
             "notes": ["All six supplied original page images were opened at original detail.", "No SI was supplied/associated.",
                       "Pre-existing workflow content was for another DOI and was not used as evidence."]}

(RUN / "coverage").mkdir(parents=True, exist_ok=True)
(RUN / "raw_model").mkdir(parents=True, exist_ok=True)
(RUN / "coverage/pages_0001_0006.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/visual_inventory_pages_0001_0006.json").write_text(json.dumps(inventory, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/stereo_diagnostics.json").write_text(json.dumps(diagnostics, ensure_ascii=False, indent=2) + "\n")
(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]), "coverage_pages": len(coverage["pages"])}, ensure_ascii=False))
