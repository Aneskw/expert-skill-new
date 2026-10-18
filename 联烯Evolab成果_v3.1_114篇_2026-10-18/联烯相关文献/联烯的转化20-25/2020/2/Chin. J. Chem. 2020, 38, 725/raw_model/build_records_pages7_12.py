#!/usr/bin/env python3
"""Append source-derived v3.1 records for PDF pages 7--12.

The page images were inspected in the current session.  This remains a machine
extraction: every reaction stays needs_review, drawing-defined complex
structures are not promoted to verified SMILES, and unavailable SI details are
left missing.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_bd193f10bf2e3fd4"
RUN = ROOT / "runs/current"
BUNDLE = st.read(RUN / "standard_records.json")
TEMPLATE = st.read(RUN / "raw_model" / "template_pages_0007_0012.json")
# The prepared template is authoritative for document/evidence identity.  The
# prior submitted bundle covered pages 1--6 and therefore did not yet carry the
# page-image evidence entries for pages 7--12.
BUNDLE["documents"] = TEMPLATE["documents"]
# This page-bounded builder is rerunnable after a merge: retain the prior
# pages-1--6 base and regenerate the pages-7--12 increment deterministically.
BUNDLE["reactions"] = BUNDLE["reactions"][:52]
BUNDLE["molecules"] = BUNDLE["molecules"][:102]
BUNDLE["source_context"] = [
    item for item in BUNDLE.get("source_context", [])
    if not any(f":page{page}:" in evidence
               for page in range(7, 13)
               for evidence in item.get("evidence_ids", []))
]
MOLS = {m["molecule_id"]: m for m in BUNDLE["molecules"]}
RXNS = {r["reaction_id"]: r for r in BUNDLE["reactions"]}


def ev(page: int) -> list[str]:
    return [f"{DOC}:page{page}:text", f"{DOC}:page{page}:image"]


def evs(*pages: int) -> list[str]:
    return list(dict.fromkeys(e for p in pages for e in ev(p)))


def Q(number, unit, *, qualifier=None, raw_text=None):
    value = {"number": number, "unit": unit}
    if qualifier:
        value["qualifier"] = qualifier
    if raw_text:
        value["raw_text"] = raw_text
    return value


def QR(minimum, maximum, unit, raw_text):
    return {"minimum": minimum, "maximum": maximum, "unit": unit, "raw_text": raw_text}


def add_molecule(mid, label, pages, name, *, chirality="unresolved",
                 geometry="unresolved", representation="not_encoded",
                 assignment=None, mixture=None, note=None, smiles=None):
    if mid in MOLS:
        return mid
    evidence = evs(*pages)
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"] = [{
        "document_id": DOC,
        "compound_label": st.fact(label, label, evidence),
        "evidence_ids": evidence,
    }]
    m["names"]["original"] = st.fact(name, name, evidence)
    m["names"]["english"] = st.fact(name, name, evidence)
    m["names"]["abbreviations"] = st.fact()
    m["names"]["aliases"] = st.fact()
    m["cas_numbers"] = st.fact()
    m["generic_template"] = None
    if smiles:
        diag = st.stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = st.fact(
            f"candidate connectivity transcribed from source label/drawing {label}",
            smiles, evidence, method="structure_derivation",
            derivation="Candidate graph was transcribed from the supplied page; parsing is not source verification.")
        m["stereochemistry"]["isomeric_smiles"] = st.fact(
            None, diag["canonical_isomeric_smiles"], evidence,
            method="structure_derivation",
            derivation="Canonical isomeric SMILES calculated from the unreviewed candidate graph.")
        m["stereochemistry"]["identity_key"] = st.fact(
            None, diag["stereo_graph_key"], evidence,
            method="structure_derivation",
            derivation="RDKit stereo graph key; not a sample identity or drawing verification.")
        representation = "encoded"
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = st.fact()
        m["stereochemistry"]["isomeric_smiles"] = st.fact()
        m["stereochemistry"]["identity_key"] = st.fact()
    m["stereochemistry"]["atom_centers"] = (
        st.fact(f"source drawing specifies stereochemical identity for {label}",
                f"source drawing specifies stereochemical identity for {label}", evidence)
        if chirality in {"specified", "partially_specified", "racemic", "mixture"}
        else st.fact())
    m["stereochemistry"]["double_bonds"] = (
        st.fact(f"source drawing/label specifies alkene geometry for {label}",
                f"source drawing/label specifies alkene geometry for {label}", evidence)
        if geometry in {"specified", "partially_specified", "mixture"}
        else st.fact())
    m["stereochemistry"]["axial"] = st.fact()
    m["stereochemistry"]["mixture"] = (
        st.fact(mixture, mixture, evidence) if mixture else st.fact())
    m["stereochemistry"]["assignment"] = (
        st.fact(assignment, assignment, evidence) if assignment else st.fact())
    for key, value in (("chirality_status", chirality),
                       ("geometry_status", geometry),
                       ("representation_status", representation)):
        derived = value in {"achiral", "not_applicable"}
        m["stereochemistry"][key] = st.fact(
            value, value, evidence,
            method="structure_derivation" if derived else "original_extraction",
            derivation=("Status inferred from the complete visible drawing; machine record remains unverified."
                        if derived else None))
    m["identity_note"] = st.fact(
        note or "Drawing-defined identity retained by source label; exact stereochemical SMILES was not safely encoded.",
        note or "Drawing-defined identity retained by source label; exact stereochemical SMILES was not safely encoded.",
        evidence)
    m["external_properties"] = []
    MOLS[mid] = m
    return mid


# Drawing-defined intermediates and geometrical-isomer identities.
complex_molecules = [
    ("cmp34b", "34b", (7,), "BT-sulfone/BT-thioether linker 34b", "unresolved", "not_applicable"),
    ("cmp34c", "34c", (7,), "PT-sulfone/PT-thioether linker 34c", "unresolved", "not_applicable"),
    ("cmp35b_E", "35b (E major)", (7,), "E-alkene PT-thioether 35b", "specified", "specified"),
    ("cmp35b_Z", "35b (Z minor)", (7,), "Z-alkene PT-thioether 35b", "specified", "specified"),
    ("cmp37a", "37a", (7, 8, 9, 10), "PT-sulfone 37a (R = TMS)", "specified", "specified"),
    ("cmp37b", "37b", (7,), "PT-sulfone 37b (R = H)", "specified", "specified"),
    ("cmp38a_E", "38a (C14-C15 E major)", (7,), "ethyl ester polyene 38a, E-major C14-C15 isomer", "specified", "specified"),
    ("cmp38a_Z", "38a (C14-C15 Z minor)", (7,), "ethyl ester polyene 38a, Z-minor C14-C15 isomer", "specified", "specified"),
    ("cmp38b_E", "38b (C14-C15 E major)", (7,), "methyl ester polyene 38b, E-major C14-C15 isomer", "specified", "specified"),
    ("cmp38b_Z", "38b (C14-C15 Z minor)", (7,), "methyl ester polyene 38b, Z-minor C14-C15 isomer", "specified", "specified"),
    ("cmp39", "39", (7,), "intended seco-acid 39", "specified", "specified"),
    ("cmp40", "40", (8,), "dienoic acid 40", "specified", "specified"),
    ("cmp41a", "41a", (8, 9, 10), "seco-acid aldehyde 41a", "specified", "specified"),
    ("cmp41b", "41b", (8, 9, 10), "9-BBN borinate of seco-acid aldehyde 41b", "specified", "specified"),
    ("cmp42_E", "42 (C14-C15 E major)", (8, 9, 10), "seco-acid 42, E-major C14-C15 isomer", "specified", "specified"),
    ("cmp42_Z", "42 (C14-C15 Z minor)", (8, 9, 10), "seco-acid 42, Z-minor C14-C15 isomer", "specified", "specified"),
    ("cmp43", "43", (8,), "3-(benzothiazol-2-ylthio)propan-1-ol 43", "achiral", "not_applicable"),
    ("cmp8a", "8a", (8, 9), "BT/PT double-headed sulfone 8a", "achiral", "not_applicable"),
    ("cmp8b", "8b", (8, 9), "bis-BT double-headed sulfone 8b", "achiral", "not_applicable"),
    ("cmp36a_Z", "36a (C17-C18 Z minor)", (9,), "BT-sulfone 36a, Z-minor C17-C18 isomer", "specified", "specified"),
    ("cmp5_lasonolide_a", "5", (9,), "(–)-lasonolide A", "specified", "specified"),
    ("mech_A", "intermediate A", (8,), "proposed bissulfone anion intermediate A", "unresolved", "not_applicable"),
    ("mech_B", "intermediate B", (8,), "proposed spiro degradation intermediate B", "unresolved", "not_applicable"),
    ("tmpl_rcho_scheme18", "RCHO", (8,), "generic aldehyde RCHO in Scheme 18", "unresolved", "not_applicable"),
    ("mech_syn_IA", "syn-IA", (12,), "proposed organoborane intermediate syn-IA", "specified", "not_applicable"),
    ("mech_syn_IB", "syn-IB", (12,), "proposed isomerized alkene intermediate syn-IB", "specified", "specified"),
]
for mid, label, pages, name, chir, geom in complex_molecules:
    assignment = None
    mixture = None
    if mid == "cmp40":
        assignment = "Authors report single-crystal analysis established the tetrahydropyran stereogenic center and three alkene geometries."
    if mid.endswith("_Z") or mid.endswith("_E"):
        assignment = "Geometrical-isomer identity follows the source E/Z outcome order; exact isomeric SMILES was not encoded."
    add_molecule(mid, label, pages, name, chirality=chir, geometry=geom,
                 assignment=assignment, mixture=mixture)

# Simple named participants.  Only safe, non-stereogenic graphs are encoded.
simple = [
    ("solv_dme", "DME", (7, 9, 10), "1,2-dimethoxyethane", "COCCOC"),
    ("solv_dmf", "DMF", (8, 9), "N,N-dimethylformamide", "CN(C)C=O"),
    ("reag_propanediol", "1,3-diol", (8,), "propane-1,3-diol", "OCCCO"),
    ("reag_hcl", "HCl", (8, 9, 10), "hydrochloric acid", "Cl"),
    ("reag_lioh", "LiOH", (7,), "lithium hydroxide", "[Li+].[OH-]"),
    ("reag_tbsoh", "TBAOH", (7,), "tetrabutylammonium hydroxide", "CCCC[N+](CCCC)(CCCC)CCCC.[OH-]"),
    ("reag_pyridine", "pyridine", (9,), "pyridine", "n1ccccc1"),
    ("reag_trichlorobenzoyl_chloride", "2,4,6-trichlorobenzoyl chloride", (9,), "2,4,6-trichlorobenzoyl chloride", "O=C(Cl)c1c(Cl)cc(Cl)cc1Cl"),
]
for mid, label, pages, name, smiles in simple:
    add_molecule(mid, label, pages, name, chirality="achiral",
                 geometry="not_applicable", representation="encoded", smiles=smiles,
                 note="Simple source-named participant; candidate graph retained for diagnostics, not source verification.")

unresolved_agents = [
    ("reag_hmpa", "HMPA", (7, 9, 10), "hexamethylphosphoramide"),
    ("reag_dmpu", "DMPU", (9,), "N,N′-dimethylpropyleneurea"),
    ("reag_bu3sn2o", "(Bu3Sn)2O", (7,), "bis(tributyltin) oxide"),
    ("reag_mcpba", "MCPBA", (8,), "m-chloroperoxybenzoic acid"),
    ("reag_dead", "DEAD", (8,), "diethyl azodicarboxylate"),
    ("reag_diad", "DIAD", (8,), "diisopropyl azodicarboxylate"),
    ("reag_pph3", "PPh3", (8,), "triphenylphosphine"),
    ("reag_btsh", "BTSH", (8,), "2-mercaptobenzothiazole"),
    ("reag_ptsh", "PTSH", (8,), "5-mercapto-1-phenyl-1H-tetrazole"),
    ("reag_hfpy", "HF·py", (9,), "hydrogen fluoride–pyridine"),
    ("reag_yamaguchi", "Yamaguchi conditions", (9,), "Yamaguchi macrolactonization conditions"),
]
for mid, label, pages, name in unresolved_agents:
    add_molecule(mid, label, pages, name, chirality="unresolved",
                 geometry="not_applicable", representation="not_encoded",
                 note="Source-named agent retained without inventing a unique formulation or stereochemical graph.")


def new_reaction(suffix, pages, *, record_type="specific_experiment",
                 categories=("other",), name=None, transformation=None, location=None):
    rid = f"{DOC}:{suffix}"
    if rid in RXNS:
        raise ValueError(f"duplicate reaction id {rid}")
    evidence = evs(*pages)
    r = st.new_reaction(rid, DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = st.fact(list(categories), list(categories), evidence)
    loc = location or {"pages": list(pages)}
    r["location"] = st.fact(str(loc), loc, evidence)
    r["reaction_name"] = st.fact(name, name, evidence) if name else st.fact()
    r["transformation_type"] = (
        st.fact(None, transformation, evidence, method="inference",
                derivation="Transformation description inferred from the source arrow and adjacent prose.")
        if transformation else st.fact())
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    RXNS[rid] = r
    return r


def add_participant(r, pid, mid, role, pages, *, stage="reaction", mass=None,
                    amount=None, volume=None, equiv=None, mol_percent=None,
                    ref=None, product_role=None):
    p = st.new_participant(pid, mid, role, evs(*pages), stage=stage)
    if product_role:
        p["product_role"] = st.fact(product_role, product_role, evs(*pages))
    for key, value in (("mass", mass), ("amount", amount), ("volume", volume),
                       ("equivalents", equiv), ("mol_percent", mol_percent)):
        if value is not None:
            raw = value.get("raw_text", str(value)) if isinstance(value, dict) else str(value)
            p["charges"][key] = st.fact(raw, value, evs(*pages))
    if equiv is not None or mol_percent is not None:
        p["equivalence_reference"] = st.fact(ref, ref, evs(*pages))
    r["participants"].append(p)


def add_step(r, sid, pages, pids, operation, *, temperature=None, time=None,
             stage="reaction"):
    s = st.new_step(sid, evs(*pages), stage=stage)
    s["participant_ids"] = pids
    s["operation"] = st.fact(operation, operation, evs(*pages))
    if temperature is not None:
        s["temperature"] = st.fact(str(temperature), temperature, evs(*pages))
    if time is not None:
        s["time"] = st.fact(str(time), time, evs(*pages))
    r["steps"].append(s)


def set_solvent(r, pages, composition):
    r["solvent_conditions"]["composition"] = st.fact(
        str(composition), composition, evs(*pages))


def set_scale(r, pages, mid, amount):
    r["reaction_scale"] = st.fact(
        f"{mid}: {amount}", {"reference_molecule_id": mid, "amount": amount}, evs(*pages))


def add_outcome(r, oid, kind, targets, pages, value, *, raw=None,
                yield_type=None, basis=None, scope=None, alternatives=()):
    o = st.new_outcome(oid, kind, targets, evs(*pages))
    if alternatives:
        o["value"] = st.fact(raw, None, evs(*pages), status="conflict", alternatives=alternatives)
    else:
        o["value"] = st.fact(raw if raw is not None else str(value), value, evs(*pages))
    if yield_type:
        o["yield_type"] = st.fact(yield_type, yield_type, evs(*pages))
    if basis:
        o["measurement_basis"] = st.fact(basis, basis, evs(*pages))
    if scope:
        o["scope"] = st.fact(str(scope), scope, evs(*pages))
    r["outcomes"].append(o)


def note(r, text, pages):
    r["notes"].append(st.fact(text, text, evs(*pages)))


LABELS = {
    "C-H": (1, "C", "H", 1, "hydrogen_gain"),
    "C-C": (2, "C", "C", 1, "formed"),
    "C-O": (3, "C", "O", 1, "formed"),
    "C-S": (4, "C", "S", 1, "formed"),
    "DB-C-C": (13, "C", "C", 2, "formed"),
    "DB-C-O": (14, "C", "O", 2, "order_increased"),
    "DB-O-S": (18, "O", "S", 2, "formed"),
    "O-H": (24, "O", "H", 1, "hydrogen_gain"),
}


def classify(r, labels, pages, limitation=None):
    changes, normalized = [], []
    for i, code in enumerate(labels, 1):
        number, a, b, order, kind = LABELS[code]
        change = {"change_id": f"b{i}", "atoms": [{"atom_map": None, "element": a},
                  {"atom_map": None, "element": b}],
                  "before_order": 0 if kind in {"formed", "hydrogen_gain"} else 1,
                  "after_order": order, "change_type": kind}
        if kind == "hydrogen_gain":
            change["hydrogen_count"] = 1
        changes.append(st.fact(None, change, evs(*pages), method="structure_derivation",
                               derivation="Partial net bond-change class inferred from the source drawing; atoms are not fully mapped."))
        normalized.append({"category_number": number, "code": code,
                           "supporting_change_ids": [f"b{i}"]})
    r["classification"]["bond_changes"] = changes
    r["classification"]["labels"] = st.fact(
        None, normalized, evs(*pages), method="structure_derivation",
        derivation="Taxonomy labels inferred from the visible net transformation; no complete atom mapping claimed.")
    r["classification"]["coverage"] = st.fact(
        None, "partial", evs(*pages), method="inference",
        derivation="Complex drawing-defined structures and atom correspondence remain unresolved.")
    r["classification"]["limitations"] = [
        limitation or "Partial source-drawing classification; full atom mapping and all implicit-H changes are unresolved."]


def classify_na(r, pages, reason):
    r["classification"]["labels"] = st.fact(None, [], evs(*pages), method="inference", derivation=reason)
    r["classification"]["coverage"] = st.fact(None, "not_applicable", evs(*pages), method="inference", derivation=reason)
    r["classification"]["limitations"] = [reason]


def classify_unresolved(r, pages, reason):
    r["classification"]["coverage"] = st.fact(None, "unresolved", evs(*pages), method="inference", derivation=reason)
    r["classification"]["limitations"] = [reason]


# Scheme 12: two first Julia variants and oxidation of the PT thioether.
for suffix, linker, products, temp, tm, yld, scheme_ratio in [
    ("34b", "cmp34b", ("cmp35_E", "cmp35_Z"), QR(-70, -30, "°C", "–70 to –30 °C"), Q(6, "h"), 63, "7/1"),
    ("34c", "cmp34c", ("cmp35b_E", "cmp35b_Z"), Q(-70, "°C"), Q(3, "h"), 82, ">20/1"),
]:
    r = new_reaction(f"scheme12_21_{suffix}_to_35", (7,), name="Julia/Julia–Kocienski olefination",
                     transformation="acetal conversion followed by aldehyde–sulfone olefination",
                     location={"page": 7, "region": f"Scheme 12, {suffix} branch"})
    for args in [("s12_21", "cmp21", "reactant"), ("s12_link", linker, "reactant"),
                 ("s12_e", products[0], "product"), ("s12_z", products[1], "product"),
                 ("s12_tms", "reag_tmsotf", "reagent"), ("s12_lut", "reag_26lutidine", "reagent"),
                 ("s12_li", "reag_lihmds", "reagent"), ("s12_dcm", "solv_dcm", "solvent"),
                 ("s12_dme", "solv_dme", "solvent"), ("s12_hmpa", "reag_hmpa", "additive")]:
        add_participant(r, *args, (7,))
    set_solvent(r, (7,), ["CH2Cl2", "DME/HMPA"])
    add_step(r, "acetal_to_aldehyde", (7,), ["s12_21", "s12_tms", "s12_lut", "s12_dcm"],
             "TMSOTf/2,6-lutidine treatment in CH2Cl2")
    add_step(r, "julia", (7,), ["s12_link", "s12_li", "s12_dme", "s12_hmpa", "s12_e", "s12_z"],
             "LiHMDS-promoted Julia olefination in DME/HMPA", temperature=temp, time=tm)
    add_outcome(r, f"yield_{suffix}", "yield", list(products), (7,), Q(yld, "%"),
                raw=f"{yld}%", scope={"span": f"21 and {suffix} to products"})
    if suffix == "34c":
        alternatives = [
            st.fact("Scheme 12: E/Z >20/1", {"ratio": ">20/1"}, evs(7)),
            st.fact("prose: E/Z 11/1", {"ratio": "11/1"}, evs(7)),
        ]
        add_outcome(r, "ez_34c", "selectivity", list(products), (7,), None,
                    raw="Scheme 12: E/Z >20/1; prose: E/Z 11/1",
                    scope={"metric": "E/Z", "order": list(products)}, alternatives=alternatives)
        note(r, "Source conflict retained: Scheme 12 reports E/Z >20/1 for 34c, whereas the adjacent prose reports E/Z 11/1.", (7,))
        r["review"]["unresolved"].append("34c E/Z conflicts between Scheme 12 (>20/1) and prose (11/1).")
    else:
        add_outcome(r, "ez_34b", "selectivity", list(products), (7,), {"ratio": scheme_ratio},
                    raw="E/Z 7/1", scope={"metric": "E/Z", "order": list(products)})
    r["review"]["unresolved"] += ["Step-1 time/temperature and reagent charges are unreported.",
                                     "Yield and E/Z measurement bases are unreported on page 7."]
    classify(r, ["DB-C-C"], (7,))

r = new_reaction("scheme12_35b_to_37a_37b", (7,), name="thioether oxidation",
                 transformation="PT-thioether oxidation to PT-sulfones",
                 location={"page": 7, "region": "Scheme 12 lower oxidation"})
for args in [("s12b_35", "cmp35b_E", "reactant"), ("s12b_37a", "cmp37a", "product"),
             ("s12b_37b", "cmp37b", "product"), ("s12b_mo", "reag_ammonium_molybdate", "reagent"),
             ("s12b_h2o2", "reag_h2o2", "reagent")]:
    add_participant(r, *args, (7,))
add_step(r, "oxidation", (7,), ["s12b_35", "s12b_mo", "s12b_h2o2", "s12b_37a", "s12b_37b"],
         "ammonium heptamolybdate/H2O2 oxidation; neutral work-up favors 37a")
add_outcome(r, "yield37a", "yield", ["cmp37a"], (7,), Q(61, "%"), raw="37a, 61%")
add_outcome(r, "yield37b", "yield", ["cmp37b"], (7,), Q(22, "%"), raw="37b, 22%")
note(r, "The prose says 37a or 37b can be isolated under optimized work-up conditions; whether these yields derive from one batch or distinct work-ups is not explicit.", (7,))
r["review"]["unresolved"] = ["Relationship between the 37a and 37b isolations, exact work-up, charges, solvent, temperature, time and yield bases are unreported."]
classify(r, ["DB-O-S"], (7,))

# Scheme 13, two substrate/product variants.
for suffix, substrate, products, yld, ratio in [
    ("38a", "cmp32a", ("cmp38a_E", "cmp38a_Z"), 89, "4.6/1"),
    ("38b", "cmp32b", ("cmp38b_E", "cmp38b_Z"), 73, "4.5/1"),
]:
    r = new_reaction(f"scheme13_37a_32_to_{suffix}", (7,), name="second Julia–Kocienski olefination",
                     transformation="aldehyde–sulfone olefination to polyene",
                     location={"page": 7, "region": f"Scheme 13, {suffix}"})
    for args in [("s13_37", "cmp37a", "reactant"), ("s13_32", substrate, "reactant"),
                 ("s13_e", products[0], "product"), ("s13_z", products[1], "product"),
                 ("s13_k", "reag_khmds", "reagent"), ("s13_thf", "solv_thf", "solvent")]:
        add_participant(r, *args, (7,))
    set_solvent(r, (7,), "THF")
    add_step(r, "base_pretreatment", (7,), ["s13_37", "s13_k", "s13_thf"],
             "treat 37a with KHMDS in THF", temperature=Q(-78, "°C"))
    add_step(r, "aldehyde_addition", (7,), ["s13_32", "s13_e", "s13_z"],
             "immediately react with aldehyde 32a or 32b", temperature=Q(-78, "°C"))
    add_outcome(r, f"yield_{suffix}", "yield", list(products), (7,), Q(yld, "%"), raw=f"{yld}%")
    add_outcome(r, f"ez_{suffix}", "selectivity", list(products), (7,), {"ratio": ratio},
                raw=f"E/Z = {ratio} for C14/C15", scope={"metric": "C14-C15 E/Z", "order": list(products)})
    r["review"]["unresolved"] = ["Charges, reaction time, work-up, purification and measurement bases are unreported on page 7."]
    classify(r, ["DB-C-C"], (7,))

# Scheme 14 is a grouped failed screen, not a fabricated set of row-resolved experiments.
r = new_reaction("scheme14_hydrolysis_screen_38_to_39", (7,), record_type="unresolved",
                 categories=("optimization",), name="failed ester hydrolysis screen",
                 transformation="attempted dienoic ester hydrolysis",
                 location={"page": 7, "region": "Scheme 14"})
for args in [("s14_38a", "cmp38a_E", "reactant"), ("s14_38b", "cmp38b_E", "reactant"),
             ("s14_lioh", "reag_lioh", "reagent"), ("s14_naoh", "reag_naoh", "reagent"),
             ("s14_tbaoh", "reag_tbsoh", "reagent"), ("s14_sn", "reag_bu3sn2o", "reagent")]:
    add_participant(r, *args, (7,))
add_step(r, "alternative_hydrolysis_trials", (7,), ["s14_38a", "s14_38b", "s14_lioh", "s14_naoh", "s14_tbaoh", "s14_sn"],
         "separate unspecified hydrolysis trials using LiOH, NaOH, TBAOH, (Bu3Sn)2O, etc.; the page does not resolve individual conditions")
add_outcome(r, "no_desired_39", "no_reaction", [], (7,), True, raw="hydrolysis to 39 problematic; desired ester not hydrolyzed")
note(r, "This arrow summarizes multiple unspecified trials; listed reagents are alternatives, not a single reaction mixture.", (7,))
r["review"]["unresolved"] = ["Individual hydrolysis trials, solvents, charges, times and outcomes are not enumerated; undesired side-chain hydrolysis product is not assigned a source label."]
classify_na(r, (7,), "Grouped failed screen; no desired 39 formation was reported and no expected-product label is assigned.")

# Scheme 15 contains cited literature chemistry and is kept as context rather than current experiments.
BUNDLE["source_context"].extend([
    {"context_id": f"{DOC}:page7_scheme15_kouster", "evidence_ids": evs(7),
     "kind": "literature_example",
     "description": "Scheme 15 equation 1 reproduces Köster (1993) borane–carboxylic-acid borinate chemistry; it is cited prior work, not a current-paper experiment."},
    {"context_id": f"{DOC}:page7_scheme15_brown", "evidence_ids": evs(7),
     "kind": "literature_example",
     "description": "Scheme 15 equation 2 reproduces Brown (1992) dicyclohexylchloroborane/Et3N enolboration chemistry; it is cited prior work, not a current-paper experiment."},
])

# Scheme 16: telescoped preparation of 40, oxidation, failed direct Julia trials,
# base-tolerance control, and the two protected-acid Julia variants.
r = new_reaction("scheme16_10_to_40", (8,), name="telescoped reduction/HWE/hydrolysis",
                 transformation="diester reduction, HWE chain elongation and ester hydrolysis",
                 location={"page": 8, "region": "Scheme 16, 10 to 40 and adjacent prose"})
for args in [("s16_10", "cmp10", "reactant"), ("s16_30", "cmp30", "reagent"),
             ("s16_9b", "cmp9b", "reactant"), ("s16_40", "cmp40", "product"),
             ("s16_dibal", "reag_dibalh", "reagent"), ("s16_nh4", "reag_nh4cl", "workup_agent"),
             ("s16_li", "reag_lihmds", "reagent"), ("s16_naoh", "reag_naoh", "reagent"),
             ("s16_thf", "solv_thf", "solvent"), ("s16_meoh", "solv_meoh", "solvent")]:
    stage = "workup" if args[1] == "reag_nh4cl" else "reaction"
    add_participant(r, *args, (8,), stage=stage)
set_solvent(r, (8,), ["THF", "MeOH"])
add_step(r, "reduction", (8,), ["s16_10", "s16_dibal", "s16_30"], "DIBAL-H reduction of 10 to aldehyde 30")
add_step(r, "quench", (8,), ["s16_nh4"], "then NH4Cl", stage="workup")
add_step(r, "hwe", (8,), ["s16_30", "s16_9b", "s16_li", "s16_thf"], "LiHMDS-promoted HWE reaction with 9b in THF")
add_step(r, "hydrolysis", (8,), ["s16_naoh", "s16_meoh", "s16_40"], "then NaOH/MeOH hydrolysis")
add_outcome(r, "yield40_overall", "yield", ["cmp40"], (8,), Q(67, "%"), raw="67% (2 steps)",
            scope={"span": "10 to 40", "reported_steps": 2})
add_outcome(r, "yield40_from30", "yield", ["cmp40"], (8,), Q(83, "%"), raw="prose: 30 to 40 in 83% yield",
            scope={"span": "30 to 40", "basis": "prose stage yield"})
add_outcome(r, "ez40", "selectivity", ["cmp40"], (8,), {"ratio": "11/1"}, raw="E/Z = 11/1",
            scope={"metric": "E/Z", "isomer order not assigned": True})
note(r, "The 67% scheme outcome spans 10→40, whereas the 83% prose outcome spans aldehyde 30→40; both are retained as distinct yield scopes, not forced into a conflict.", (8,))
r["review"]["unresolved"] = ["Exact charges, temperatures, times, work-up, purification and measurement bases are unreported; the E/Z isomer order is not explicitly assigned."]
classify(r, ["C-H", "DB-C-C", "O-H"], (8,))

r = new_reaction("scheme16_40_to_41a", (8,), name="Dess–Martin oxidation",
                 transformation="allylic alcohol oxidation to aldehyde",
                 location={"page": 8, "region": "Scheme 16, 40 to 41a"})
for args in [("s16b_40", "cmp40", "reactant"), ("s16b_41", "cmp41a", "product"),
             ("s16b_dmp", "reag_dmp_oxidant", "reagent")]:
    add_participant(r, *args, (8,))
add_step(r, "oxidation", (8,), ["s16b_40", "s16b_dmp", "s16b_41"], "DMP oxidation")
add_outcome(r, "yield41_good", "yield", ["cmp41a"], (8,), "good yield", raw="good yield")
r["review"]["unresolved"] = ["Numerical yield, yield basis, charges, solvent, temperature, time, work-up and purification are unreported."]
classify(r, ["DB-C-O"], (8,))

for sulfone in ("cmp36a", "cmp36b"):
    tag = "36a" if sulfone == "cmp36a" else "36b"
    r = new_reaction(f"scheme16_direct_julia_41a_{tag}_no_reaction", (8,), categories=("control",),
                     name="direct Julia olefination of unprotected acid", transformation="attempted Julia olefination",
                     location={"page": 8, "region": f"Scheme 16, direct 41a + {tag} trial"})
    add_participant(r, "direct_41", "cmp41a", "reactant", (8,))
    add_participant(r, "direct_sulfone", sulfone, "reactant", (8,))
    add_step(r, "attempted_julia", (8,), ["direct_41", "direct_sulfone"],
             f"direct Julia reaction of free acid 41a with {tag}; detailed conditions not reported")
    add_outcome(r, "no_desired_product", "no_reaction", [], (8,), True, raw="did not produce any desired coupling product")
    r["review"]["unresolved"] = ["Base, solvent, temperature, time, charges and detection method are unreported."]
    classify_na(r, (8,), "No desired coupling product was observed; no expected-product bond label assigned.")

r = new_reaction("scheme16_9bbn_base_tolerance_control", (8,), categories=("control", "mechanistic_study"),
                 name="9-BBN protected-acid base-tolerance control",
                 transformation="reversible borinate protection under Julia-basic conditions",
                 location={"page": 8, "region": "Scheme 16, 41a/41b recovery control"})
for args in [("ctrl_41a", "cmp41a", "reactant"), ("ctrl_41b", "cmp41b", "reagent"),
             ("ctrl_9bbn", "reag_9bbn", "reagent"), ("ctrl_k", "reag_khmds", "reagent"),
             ("ctrl_hcl", "reag_hcl", "workup_agent")]:
    eq = Q(1.0, "equiv") if args[1] == "reag_9bbn" else QR(4, 6, "equiv", "4–6 equiv") if args[1] == "reag_khmds" else None
    stage = "workup" if args[1] == "reag_hcl" else "reaction"
    add_participant(r, *args, (8,), equiv=eq, ref="cmp41a" if eq else None, stage=stage)
add_step(r, "borinate_pretreatment", (8,), ["ctrl_41a", "ctrl_9bbn", "ctrl_41b"], "pretreat 41a with 9-BBN (1.0 equiv)")
add_step(r, "base_addition", (8,), ["ctrl_41b", "ctrl_k"], "mix borinate 41b with KHMDS (4–6 equiv)", temperature=Q(-45, "°C"))
add_step(r, "warm", (8,), ["ctrl_41b", "ctrl_k"], "warm from –45 °C to room temperature and hold for 30 min; warming-ramp duration not separately reported", time=Q(30, "min"))
add_step(r, "acid_quench", (8,), ["ctrl_hcl"], "quench with aqueous HCl", stage="workup")
add_outcome(r, "recovery41a", "recovery", ["cmp41a"], (8,), Q(90, "%", qualifier="gt", raw_text=">90%"), raw=">90% yield")
add_outcome(r, "epimerization", "other", ["cmp41a"], (8,), "without appreciable epimerization", raw="without appreciable epimerization")
r["review"]["unresolved"] = ["Recovery measurement basis, borinate pretreatment temperature/time, warming profile, work-up and purification are unreported."]
classify_na(r, (8,), "Starting acid 41a was recovered after reversible protection; no net product-forming bond label assigned.")

def julia_41_record(suffix, sulfone, ratio, yld, pages=(8,), scale=None):
    r = new_reaction(suffix, pages, categories=("scale_up",) if scale else ("other",),
                     name="Julia–Kocienski olefination with in-situ borinate protection",
                     transformation="aldehyde–sulfone olefination under reversible acid protection",
                     location={"pages": list(pages), "region": "Scheme 16 condition A/B" if pages == (8,) else "Scheme 16 condition A and Scheme 19 scale-up"})
    for args in [("j41", "cmp41a", "reactant"), ("j41b", "cmp41b", "reagent"),
                 ("jsulf", sulfone, "reactant"), ("j42e", "cmp42_E", "product"),
                 ("j42z", "cmp42_Z", "product"), ("j9bbn", "reag_9bbn", "reagent"),
                 ("jli", "reag_lihmds", "reagent"), ("jhcl", "reag_hcl", "workup_agent")]:
        eq = Q(1.0, "equiv") if args[1] == "reag_9bbn" else None
        stage = "workup" if args[1] == "reag_hcl" else "reaction"
        add_participant(r, *args, pages, equiv=eq, ref="cmp41a" if eq else None, stage=stage)
    if scale:
        set_scale(r, pages, "cmp41a", scale)
    add_step(r, "borinate", pages, ["j41", "j9bbn", "j41b"], "generate borinate 41b in situ with 9-BBN")
    add_step(r, "sulfone_deprotonation", pages, ["jsulf", "jli"], "deprotonate sulfone with LiHMDS")
    add_step(r, "coupling", pages, ["j41b", "jsulf", "j42e", "j42z"], "react deprotonated sulfone with in-situ borinate aldehyde 41b")
    add_step(r, "acid_workup", pages, ["jhcl"], "aqueous HCl work-up", stage="workup")
    add_outcome(r, "yield42", "yield", ["cmp42_E", "cmp42_Z"], pages, Q(yld, "%"),
                yield_type="isolated", basis="isolated" if suffix.endswith("scaleup") or yld == 65 else None, raw=f"{yld}%")
    add_outcome(r, "ez42", "selectivity", ["cmp42_E", "cmp42_Z"], pages, {"ratio": ratio},
                raw=f"E/Z = {ratio} for C14-C15", scope={"metric": "C14-C15 E/Z", "order": ["cmp42_E", "cmp42_Z"]})
    r["review"]["unresolved"] = ["Exact charges beyond 9-BBN, solvent, temperatures, times and E/Z measurement method are unreported in Scheme 16."]
    classify(r, ["DB-C-C"], pages)
    return r

# Condition A is the same 65%/12:1 preparative experiment repeated in Scheme 19;
# evidence is merged rather than counted twice.
julia_41_record("scheme16_19_conditionA_scaleup", "cmp36a", "12/1", 65, pages=(8, 9), scale=Q(200, "mg"))
julia_41_record("scheme16_conditionB", "cmp36b", "5.4/1", 43, pages=(8,))

# Scheme 17, three preparative arrows.
r = new_reaction("scheme17_propanediol_to_43", (8,), name="Mitsunobu thioether formation",
                 transformation="alcohol substitution by BTSH",
                 location={"page": 8, "region": "Scheme 17, 1,3-diol to 43"})
for args in [("s17_diol", "reag_propanediol", "reactant"), ("s17_43", "cmp43", "product"),
             ("s17_dead", "reag_dead", "reagent"), ("s17_pph3", "reag_pph3", "reagent"),
             ("s17_btsh", "reag_btsh", "reactant")]:
    add_participant(r, *args, (8,))
add_step(r, "mitsunobu", (8,), ["s17_diol", "s17_dead", "s17_pph3", "s17_btsh", "s17_43"], "DEAD/PPh3/BTSH Mitsunobu reaction")
add_outcome(r, "yield43", "yield", ["cmp43"], (8,), Q(83, "%"), raw="83%")
r["review"]["unresolved"] = ["Stoichiometry, solvent, temperature, time, work-up, purification and yield basis are unreported."]
classify(r, ["C-S"], (8,))

for product, azo, thiol, yld in [("cmp8b", "reag_dead", "reag_btsh", 76), ("cmp8a", "reag_diad", "reag_ptsh", 85)]:
    tag = "8b" if product == "cmp8b" else "8a"
    r = new_reaction(f"scheme17_43_to_{tag}", (8,), name="Mitsunobu substitution/oxidation",
                     transformation="second thioether installation followed by sulfone oxidation",
                     location={"page": 8, "region": f"Scheme 17, 43 to {tag}"})
    for args in [("s17b_43", "cmp43", "reactant"), ("s17b_prod", product, "product"),
                 ("s17b_azo", azo, "reagent"), ("s17b_pph3", "reag_pph3", "reagent"),
                 ("s17b_thiol", thiol, "reactant"), ("s17b_mcpba", "reag_mcpba", "reagent")]:
        add_participant(r, *args, (8,))
    add_step(r, "mitsunobu", (8,), ["s17b_43", "s17b_azo", "s17b_pph3", "s17b_thiol"], "Mitsunobu substitution")
    add_step(r, "oxidation", (8,), ["s17b_mcpba", "s17b_prod"], "then MCPBA oxidation")
    add_outcome(r, f"yield{tag}", "yield", [product], (8,), Q(yld, "%"), raw=f"{yld}% (2 steps)",
                scope={"span": f"43 to {tag}", "reported_steps": 2})
    r["review"]["unresolved"] = ["Stoichiometry, solvent, temperature, time, work-up, purification and yield basis are unreported."]
    classify(r, ["C-S", "DB-O-S"], (8,))

r = new_reaction("scheme18_proposed_bissulfone_degradation", (8,), record_type="proposed_mechanism",
                 name="proposed degradation pathway of bissulfone 8a",
                 transformation="base-induced spiro degradation competing with Julia olefination",
                 location={"page": 8, "region": "Scheme 18"})
for args in [("m18_rcho", "tmpl_rcho_scheme18", "reactant"),
             ("m18_8a", "cmp8a", "reactant"), ("m18_a", "mech_A", "reagent"),
             ("m18_b", "mech_B", "product"), ("m18_36a", "cmp36a", "product"),
             ("m18_k", "reag_khmds", "reagent")]:
    eq = (Q(1, "equiv") if args[1] == "tmpl_rcho_scheme18" else
          Q(3, "equiv") if args[1] == "cmp8a" else
          Q(2.5, "equiv") if args[1] == "reag_khmds" else None)
    add_participant(r, *args, (8,), equiv=eq,
                    ref="tmpl_rcho_scheme18" if eq else None)
add_step(r, "deprotonation", (8,), ["m18_8a", "m18_k", "m18_a"], "proposed KHMDS deprotonation of 8a to intermediate A")
add_step(r, "intended_aldehyde_trapping", (8,), ["m18_a", "m18_rcho", "m18_36a"], "desired/proposed trapping of intermediate A with generic aldehyde RCHO to give 36a")
add_step(r, "spiro_degradation", (8,), ["m18_a", "m18_b"], "proposed fast intramolecular formation of spiro intermediate B followed by slow decomposition")
note(r, "This is an author-proposed competition between desired aldehyde trapping and spiro degradation; RCHO is generic and no observed yield is reported.", (8,))
r["review"]["unresolved"] = ["Generic R group, exact atom correspondence and experimental observation status of intermediates A/B are unresolved."]
classify_unresolved(r, (8,), "Proposed generic mechanism with unresolved R group and no safe atom mapping.")

# Table 4: every row is a distinct optimization experiment.
table4 = [
    (1, "reag_khmds", 2.5, "cmp8a", "DME/HMPA", -70, "0", None, None),
    (2, "reag_khmds", 2.5, "cmp8a", "DME/HMPA", -70, ">90", "8.7/1", 73),
    (3, "reag_khmds", 2.5, "cmp8a", "DME", -70, ">90", "4.9/1", None),
    (4, "reag_khmds", 2.5, "cmp8a", "THF/HMPA", -78, ">90", "5.0/1", None),
    (5, "reag_khmds", 2.5, "cmp8a", "DME/18-c-6", -70, ">90", "5.4/1", None),
    (6, "reag_khmds", 2.5, "cmp8a", "DME/HMPA/18-c-6", -70, ">90", "7.3/1", None),
    (7, "reag_khmds", 2.5, "cmp8a", "DMF/HMPA", -60, ">90", ">20/1", 83),
    (8, "reag_khmds", 3.0, "cmp8b", "THF/HMPA", -78, "44", "2.1/1", None),
    (9, "reag_lihmds", 3.0, "cmp8b", "THF/HMPA", -78, "39", "3.9/1", 30),
    (10, "reag_lihmds", 6.0, "cmp8b", "THF/HMPA", -78, "<10", None, None),
    (11, "reag_lihmds", 3.0, "cmp8b", "DME/HMPA", -70, "28", None, None),
    (12, "reag_lihmds", 3.0, "cmp8b", "DMF/DMPU", -70, "<5", None, None),
]


def solvent_participants(r, prefix, text, pages):
    mapping = {"DME": "solv_dme", "THF": "solv_thf", "DMF": "solv_dmf",
               "HMPA": "reag_hmpa", "18-c-6": "reag_18c6", "DMPU": "reag_dmpu"}
    mids = []
    for i, token in enumerate(text.split("/"), 1):
        mid = mapping[token]
        role = "solvent" if token in {"DME", "THF", "DMF"} else "additive"
        pid = f"{prefix}_{i}"
        add_participant(r, pid, mid, role, pages)
        mids.append(pid)
    return mids


for entry, base, beq, sulfone, solvent, temp, conv, ratio, yld in table4:
    r = new_reaction(f"table4_entry{entry}", (9,), categories=("optimization",),
                     name="sequential Julia olefination optimization, step 1",
                     transformation="bissulfone–aldehyde Julia olefination",
                     location={"page": 9, "region": f"Table 4 entry {entry}"})
    add_participant(r, "t4_33", "cmp33", "reactant", (9,), amount=Q(0.05, "mmol"), equiv=Q(1, "equiv"), ref="cmp33")
    add_participant(r, "t4_sulf", sulfone, "reactant", (9,), equiv=Q(3.0, "equiv"), ref="cmp33")
    add_participant(r, "t4_base", base, "reagent", (9,), equiv=Q(beq, "equiv"), ref="cmp33")
    set_scale(r, (9,), "cmp33", Q(0.05, "mmol"))
    spids = solvent_participants(r, "t4_solv", solvent, (9,))
    if entry != 1:
        add_participant(r, "t4_e", "cmp36a", "product", (9,), product_role="C17-C18 E isomer")
        if ratio:
            add_participant(r, "t4_z", "cmp36a_Z", "product", (9,), product_role="C17-C18 Z isomer")
    set_solvent(r, (9,), solvent)
    if entry == 1:
        add_step(r, "premix_sulfone_base", (9,), ["t4_sulf", "t4_base"] + spids,
                 "premix bissulfone and base; unlike the other entries")
        add_step(r, "aldehyde_addition", (9,), ["t4_33"], "add aldehyde 33 after the premixed bissulfone and base", temperature=Q(temp, "°C"))
    else:
        add_step(r, "premix", (9,), ["t4_33", "t4_sulf"] + spids, "premix aldehyde 33 and bissulfone")
        add_step(r, "base_addition", (9,), ["t4_base"] + (["t4_e", "t4_z"] if ratio else ["t4_e"]), "add base after premixing", temperature=Q(temp, "°C"))
    if conv == "0":
        add_outcome(r, "conversion", "conversion", ["cmp33"], (9,), Q(0, "%"), basis="1H NMR (400 MHz, CDCl3)", raw="0%")
        add_outcome(r, "no_reaction", "no_reaction", [], (9,), True, raw="conversion 0%")
        classify_na(r, (9,), "Entry 1 reported 0% conversion; no observed product bond label assigned.")
    else:
        if conv.startswith(">"):
            cv = Q(float(conv[1:]), "%", qualifier="gt", raw_text=conv + "%")
        elif conv.startswith("<"):
            cv = Q(float(conv[1:]), "%", qualifier="lt", raw_text=conv + "%")
        else:
            cv = Q(float(conv), "%")
        add_outcome(r, "conversion", "conversion", ["cmp33"], (9,), cv, basis="1H NMR (400 MHz, CDCl3)", raw=conv + "%")
        if ratio:
            add_outcome(r, "ez", "selectivity", ["cmp36a", "cmp36a_Z"], (9,), {"ratio": ratio},
                        basis="1H NMR (400 MHz, CDCl3)", raw=f"E/Z {ratio}",
                        scope={"metric": "C17-C18 E/Z", "order": ["cmp36a", "cmp36a_Z"]})
        if yld is not None:
            targets = ["cmp36a", "cmp36a_Z"] if ratio else ["cmp36a"]
            add_outcome(r, "yield36a", "yield", targets, (9,), Q(yld, "%"), yield_type="isolated", basis="isolated", raw=f"{yld}%")
        if entry == 3:
            note(r, "Footnote e states that isolation/yield determination was not performed.", (9,))
        classify(r, ["DB-C-C"], (9,))
    r["review"]["unresolved"] += ["Reaction time, work-up and purification are not reported for the Table 4 row."]

# Scheme 19 final macrocyclization/desilylation.
r = new_reaction("scheme19_42_to_lasonolideA", (9,), name="Yamaguchi macrolactonization/global desilylation",
                 transformation="macrolactonization followed by global desilylation",
                 location={"page": 9, "region": "Scheme 19, 42 to (–)-lasonolide A"})
for args in [("s19_42", "cmp42_E", "reactant"), ("s19_5", "cmp5_lasonolide_a", "product"),
             ("s19_yam", "reag_yamaguchi", "reagent"),
             ("s19_tcb", "reag_trichlorobenzoyl_chloride", "reagent"),
             ("s19_hf", "reag_hfpy", "reagent")]:
    add_participant(r, *args, (9,))
add_step(r, "macrolactonization", (9,), ["s19_42", "s19_yam", "s19_tcb"], "Yamaguchi esterification using 2,4,6-trichlorobenzoyl chloride")
add_step(r, "desilylation", (9,), ["s19_hf", "s19_5"], "immediate global desilylation with HF·pyridine")
add_outcome(r, "yield5", "yield", ["cmp5_lasonolide_a"], (9,), Q(58, "%"), yield_type="isolated", basis="isolated", raw="58%")
note(r, "The 58% outcome spans macrolactonization and immediate global desilylation; it is not assigned to desilylation alone.", (9,))
r["review"]["unresolved"] = ["Exact Yamaguchi co-reagents beyond the named acid chloride, charges, solvents, temperatures, times, work-up and purification are unreported on page 9."]
classify(r, ["C-O", "O-H"], (9,))

# Table 5: every row is a distinct step-2 optimization experiment.
table5 = [
    (1, "cmp36a", "reag_lihmds", 2.0, "DME", -75, 6.5, ">95", "12/1", 63),
    (2, "cmp36a", "reag_lihmds", 2.0, "DME/HMPA", -70, 3, ">95", "6.7/1", None),
    (3, "cmp37a", "reag_khmds", 2.0, "THF", -78, 2, ">95", "5.4/1", 43),
    (4, "cmp37a", "reag_khmds", 2.0, "DME", -70, 2, ">95", "3.7/1", None),
    (5, "cmp37a", "reag_lihmds", 2.0, "DME", -70, 3, ">95", "5.1/1", None),
]
for entry, sulfone, base, beq, solvent, temp, tm, conv, ratio, yld in table5:
    r = new_reaction(f"table5_entry{entry}", (10,), categories=("optimization",),
                     name="sequential Julia olefination optimization, step 2",
                     transformation="sulfone–borinate aldehyde Julia olefination",
                     location={"page": 10, "region": f"Table 5 entry {entry}"})
    add_participant(r, "t5_sulf", sulfone, "reactant", (10,), amount=Q(0.05, "mmol"), equiv=Q(1, "equiv"), ref=sulfone)
    add_participant(r, "t5_41b", "cmp41b", "reactant", (10,), equiv=Q(1.4, "equiv"), ref=sulfone)
    add_participant(r, "t5_base", base, "reagent", (10,), equiv=Q(beq, "equiv"), ref=sulfone)
    add_participant(r, "t5_e", "cmp42_E", "product", (10,), product_role="C14-C15 E isomer")
    add_participant(r, "t5_z", "cmp42_Z", "product", (10,), product_role="C14-C15 Z isomer")
    add_participant(r, "t5_hcl", "reag_hcl", "workup_agent", (10,), stage="workup")
    set_scale(r, (10,), sulfone, Q(0.05, "mmol"))
    spids = solvent_participants(r, "t5_solv", solvent, (10,))
    set_solvent(r, (10,), solvent)
    add_step(r, "julia", (10,), ["t5_sulf", "t5_41b", "t5_base"] + spids + ["t5_e", "t5_z"],
             "second Julia olefination", temperature=Q(temp, "°C"), time=Q(tm, "h"))
    add_step(r, "acid_workup", (10,), ["t5_hcl"], "work up with aqueous 1 N HCl after the second Julia reaction", stage="workup")
    cv = Q(float(conv[1:]), "%", qualifier="gt", raw_text=conv + "%")
    add_outcome(r, "conversion", "conversion", [sulfone], (10,), cv,
                basis="1H NMR (400 MHz, CDCl3)", raw=conv + "%")
    add_outcome(r, "ez", "selectivity", ["cmp42_E", "cmp42_Z"], (10,), {"ratio": ratio},
                basis="1H NMR (400 MHz, CDCl3)", raw=f"E/Z {ratio}",
                scope={"metric": "C14-C15 E/Z", "order": ["cmp42_E", "cmp42_Z"]})
    if yld is not None:
        add_outcome(r, "yield42", "yield", ["cmp42_E", "cmp42_Z"], (10,), Q(yld, "%"),
                    yield_type="isolated", basis="isolated", raw=f"{yld}%")
    if entry == 2:
        note(r, "Footnote e states that isolation/yield determination was not performed.", (10,))
    if sulfone == "cmp37a":
        note(r, "Footnote d states that 37a was prepared by the two-step sequence in Scheme 12; those preparation steps are not duplicated here.", (10,))
    r["review"]["unresolved"] = ["Addition order, work-up details beyond 1 N HCl, and purification are unreported for this Table 5 row."]
    classify(r, ["DB-C-C"], (10,))

# Page 12 footnote [31] contains an author-proposed side-product pathway.
r = new_reaction("footnote31_proposed_23c_pathway", (12,), record_type="proposed_mechanism",
                 name="proposed reversible-hydroboration pathway to syn-23c",
                 transformation="hydroboration, beta-boron elimination, rehydroboration and oxidation",
                 location={"page": 12, "region": "reference footnote [31] mechanism drawing"})
for args in [("m31_22", "cmp22_rac_Z", "reactant"), ("m31_ia", "mech_syn_IA", "reagent"),
             ("m31_ib", "mech_syn_IB", "reagent"), ("m31_23c", "cmp23c_rac", "product"),
             ("m31_thex", "reag_thexbh2", "reagent"), ("m31_h2o2", "reag_h2o2", "reagent"),
             ("m31_naoh", "reag_naoh", "reagent")]:
    add_participant(r, *args, (12,))
add_step(r, "hydroboration", (12,), ["m31_22", "m31_thex", "m31_ia"], "proposed hydroboration of (Z)-diol-22 to syn-IA")
add_step(r, "beta_boron_elimination", (12,), ["m31_ia", "m31_ib"], "proposed β-boron elimination to isomerized syn-IB")
add_step(r, "rehydroboration_oxidation", (12,), ["m31_ib", "m31_thex", "m31_h2o2", "m31_naoh", "m31_23c"], "rapid rehydroboration followed by H2O2/NaOH oxidation to syn-23c")
note(r, "The authors state that syn-IB was not isolable; this drawing is a proposed explanation for side product 23c, not a separate observed preparation.", (12,))
r["review"]["unresolved"] = ["Intermediate structures are proposed and not isolated; exact atom mapping and boron species stoichiometry are unresolved."]
classify(r, ["C-H", "C-O"], (12,), "Proposed side-product pathway only; product bond classes are not experimental atom mapping.")

# Non-reaction page regions and explicit SI limitation.
BUNDLE["source_context"].extend([
    {"context_id": f"{DOC}:page10_experimental_si_notice", "evidence_ids": evs(10),
     "kind": "source_availability",
     "description": "The article says procedures and characterization are in Supporting Materials online, but no SI was explicitly supplied or associated in this workspace; no SI facts were imported."},
    {"context_id": f"{DOC}:page10_conclusions", "evidence_ids": evs(10),
     "kind": "discussion",
     "description": "Conclusions summarize the synthesis and do not add a separately executable experiment beyond Table 5 and Scheme 19."},
    {"context_id": f"{DOC}:page11_references", "evidence_ids": evs(11),
     "kind": "references",
     "description": "Page 11 contains bibliography entries only; cited-paper experiments were not imported as current-paper experiments."},
    {"context_id": f"{DOC}:page12_references", "evidence_ids": evs(12),
     "kind": "references",
     "description": "Page 12 is bibliography/footnotes; only the authors' explicit proposed syn-IA/syn-IB pathway in footnote [31] was represented as a proposed mechanism."},
    {"context_id": f"{DOC}:pages7_12_characterization_triage", "evidence_ids": evs(7, 8, 9, 10, 11, 12),
     "kind": "characterization_triage",
     "description": "No standalone mass-spectrum peak table, molecular-weight table, or fragment-relative-abundance table appears on pages 7–12."},
])

BUNDLE["molecules"] = list(MOLS.values())
BUNDLE["reactions"] = list(RXNS.values())
st.write(RUN / "submitted.json", BUNDLE)

# Coverage is intentionally bounded to this session's six pages.
records_by_page = {
    7: [rid for rid in RXNS if any(f":page7:" in e for e in RXNS[rid]["evidence_ids"]) and rid not in {r["reaction_id"] for r in BUNDLE["reactions"][:52]}],
    8: [rid for rid in RXNS if any(f":page8:" in e for e in RXNS[rid]["evidence_ids"]) and rid not in {r["reaction_id"] for r in BUNDLE["reactions"][:52]}],
    9: [rid for rid in RXNS if any(f":page9:" in e for e in RXNS[rid]["evidence_ids"]) and rid not in {r["reaction_id"] for r in BUNDLE["reactions"][:52]}],
    10: [rid for rid in RXNS if any(f":page10:" in e for e in RXNS[rid]["evidence_ids"]) and rid not in {r["reaction_id"] for r in BUNDLE["reactions"][:52]}],
    11: [],
    12: [f"{DOC}:footnote31_proposed_23c_pathway"],
}
regions = {
    7: [
        {"region": "Scheme 12 and adjacent prose", "finding": "Two Julia variants and the subsequent thioether oxidation were inventoried; the 34c E/Z conflict was retained.", "evidence_ids": evs(7)},
        {"region": "Scheme 13", "finding": "Ethyl- and methyl-ester Julia couplings were recorded as separate experiments.", "evidence_ids": evs(7)},
        {"region": "Scheme 14", "finding": "Grouped unsuccessful hydrolysis screen retained as unresolved because individual trials are not enumerated.", "evidence_ids": evs(7)},
        {"region": "Scheme 15", "finding": "Köster and Brown equations are cited literature examples, retained in source_context rather than current experiments.", "evidence_ids": evs(7)},
    ],
    8: [
        {"region": "Scheme 16 and prose", "finding": "Preparation of 40/41a, two failed direct Julia trials, a base-tolerance control, and protected-acid Julia conditions A/B were recorded.", "evidence_ids": evs(8)},
        {"region": "Scheme 17", "finding": "Three preparative arrows for 43, 8a and 8b were recorded independently.", "evidence_ids": evs(8)},
        {"region": "Scheme 18", "finding": "Author-proposed bissulfone degradation pathway recorded as proposed mechanism, not an observed synthesis.", "evidence_ids": evs(8)},
    ],
    9: [
        {"region": "Table 4 entries 1–12 and footnotes", "finding": "All twelve step-1 optimization rows were recorded independently with row-specific order, base, solvent, temperature and outcomes.", "evidence_ids": evs(9)},
        {"region": "Scheme 19 and prose", "finding": "The 200-mg condition-A coupling was merged with its Scheme 16 duplicate evidence; final macrolactonization/desilylation was recorded separately.", "evidence_ids": evs(9)},
        {"region": "Conclusions", "finding": "Narrative summary only; no extra experiment created.", "evidence_ids": evs(9)},
    ],
    10: [
        {"region": "Table 5 entries 1–5 and footnotes", "finding": "All five step-2 optimization rows were recorded independently, including NMR basis, isolation footnotes and 1 N HCl work-up.", "evidence_ids": evs(10)},
        {"region": "Experimental/Supporting Information notice", "finding": "The article points to online Supporting Materials, but none was explicitly provided or associated; no SI facts were inferred.", "evidence_ids": evs(10)},
        {"region": "Conclusions, acknowledgements and references", "finding": "No additional current synthesis experiment on this page.", "evidence_ids": evs(10)},
    ],
    11: [
        {"region": "References [8]–[28]", "finding": "Bibliography only; no current-paper reaction arrows, tables or experimental procedures.", "evidence_ids": evs(11)},
    ],
    12: [
        {"region": "Reference footnote [31] text and mechanism drawing", "finding": "Author-proposed reversible-hydroboration/isomerization pathway to 23c recorded as proposed mechanism; syn-IB explicitly not isolable.", "evidence_ids": evs(12)},
        {"region": "References [29]–[47] and manuscript dates", "finding": "Bibliography/metadata only; cited experiments were not imported as current-paper experiments.", "evidence_ids": evs(12)},
    ],
}
unresolved = {
    7: ["Scheme 14 does not enumerate individual hydrolysis trials or identify the undesired hydrolyzed side-chain product.", "Scheme 12 34c E/Z differs between scheme (>20/1) and prose (11/1)."],
    8: ["Most exact charges, solvents, work-ups and measurement methods require unavailable SI.", "Scheme 18 intermediates are proposed and generic."],
    9: ["Table 4 reaction times and work-up/purification are absent; Scheme 19 exact Yamaguchi conditions require unavailable SI."],
    10: ["Table 5 addition order and purification are absent; SI was mentioned but not provided/associated."],
    11: ["No reaction content to resolve on this reference-only page."],
    12: ["The footnote [31] pathway is proposed; intermediate syn-IB was not isolated and exact mapping is unresolved."],
}
coverage = {"document_id": DOC, "semantic_accuracy_verified": False, "pages": []}
for page in range(7, 13):
    coverage["pages"].append({"page": page, "state": "attempted", "regions": regions[page],
                              "record_ids": records_by_page[page], "unresolved": unresolved[page]})
st.write(RUN / "coverage/pages_0007_0012.json", coverage)
st.write(RUN / "raw_model/page_inventory_pages_0007_0012.json",
         {"document_id": DOC, "pages": regions, "record_ids_by_page": records_by_page})
print(json.dumps({"total_molecules": len(MOLS), "total_reactions": len(RXNS),
                  "new_reactions": len(RXNS) - 52,
                  "records_by_page": {str(k): len(v) for k, v in records_by_page.items()}},
                 ensure_ascii=False, indent=2))
