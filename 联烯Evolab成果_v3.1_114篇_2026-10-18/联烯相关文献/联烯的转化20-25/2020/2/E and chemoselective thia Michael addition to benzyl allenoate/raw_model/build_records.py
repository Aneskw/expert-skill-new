#!/usr/bin/env python3
"""Build page-1--8 machine records for doc_7ba97fbb78e354a5.

This script is an extraction artifact, not a modification of the bundled tools.
All facts remain machine-generated and need review.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DOC = "doc_7ba97fbb78e354a5"
EV = {p: f"{DOC}:page{p}:text" for p in range(1, 9)}
IMG = {p: f"{DOC}:page{p}:image" for p in range(1, 9)}


def sf(raw, normalized=None, evidence=(), method="original_extraction", derivation=None,
       status=None, alternatives=()):
    if normalized is None and raw is not None and status != "conflict":
        normalized = raw
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation,
                   status=status, alternatives=alternatives)


def df(raw, normalized, evidence, derivation):
    return sf(raw, normalized, evidence, method="structure_derivation", derivation=derivation)


def inf(raw, normalized, evidence, derivation):
    return sf(raw, normalized, evidence, method="inference", derivation=derivation)


def molecule(mid, label, name, evidence, smiles=None, geometry=None, aliases=None):
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = sf(label, label, evidence)
    m["names"]["original"] = sf(name, name, evidence)
    m["names"]["english"] = sf(name, name, evidence)
    if aliases:
        m["names"]["aliases"] = sf(aliases, aliases, evidence)
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = df({"source": "name and/or structure drawing", "label": label}, smiles,
                         evidence, "Candidate SMILES manually transcribed from the cited source name/drawing; parsing is not source verification.")
        diag = st.stereo_inspect(smiles)
        if geometry in ("E", "Z"):
            m["stereochemistry"]["double_bonds"] = sf(
                f"({geometry}) product geometry", [{"configuration": geometry, "location": "but-2-enoate C=C"}], evidence)
            m["stereochemistry"]["geometry_status"] = sf(f"({geometry})", "specified", evidence)
            m["stereochemistry"]["chirality_status"] = inf(
                "No stereogenic center is shown in the source product", "achiral", evidence,
                "Candidate graph has no tetrahedral/axial stereogenic element; the E alkene is tracked separately.")
            m["stereochemistry"]["representation_status"] = inf(
                f"{geometry} geometry encoded in candidate SMILES", "encoded", evidence,
                f"Directional SMILES bonds encode the {geometry} alkene identity.")
            m["stereochemistry"]["isomeric_smiles"] = df(
                {"source": f"({geometry}) name/drawing or E/Z-ratio counterpart"}, smiles, evidence,
                f"Isomeric SMILES candidate encodes the {geometry} geometry.")
        else:
            m["stereochemistry"]["chirality_status"] = inf(
                "No stereogenic element shown", "achiral", evidence,
                "Candidate graph and source drawing show no stereogenic center or stereogenic axis.")
            m["stereochemistry"]["geometry_status"] = inf(
                "No applicable E/Z or ring-relative assignment", "not_applicable", evidence,
                "No stereogenic alkene or ring-relative geometry applies to this candidate identity.")
            m["stereochemistry"]["representation_status"] = inf(
                "No stereochemical encoding applicable", "not_applicable", evidence,
                "The candidate identity has no applicable stereochemical element to encode.")
        m["stereochemistry"]["identity_key"] = df(
            {"tool": "standard stereo", "source_verified": False}, diag["stereo_graph_key"], evidence,
            "RDKit stereo graph key calculated from the machine candidate SMILES; this is not sample identity verification.")
        m["identity_note"] = inf(
            "Machine structure candidate", "Source name/drawing recognized; exact candidate remains needs_review.", evidence,
            "The record is machine-extracted and structure parsing does not establish correctness.")
    else:
        m["representation_type"] = "unresolved"
        for key in ("chirality_status", "geometry_status", "representation_status"):
            m["stereochemistry"][key] = sf(
                "Source background drawing/name not safely encoded", "unresolved", evidence)
        m["identity_note"] = sf(
            "Exact structure not transcribed from background schematic", "Exact structure intentionally unresolved.", evidence)
    return m


def generic_molecule(mid, label, name, evidence, product=False, product_geometry="E"):
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = sf(label, label, evidence)
    m["names"]["original"] = sf(name, name, evidence)
    m["names"]["english"] = sf(name, name, evidence)
    m["representation_type"] = "generic_template"
    fixed = f"({product_geometry})-BnO2C-CH=C(CH3)-S-R" if product else "R-SH"
    m["generic_template"] = {
        "fixed_structure": sf(fixed, {"description": fixed}, evidence),
        "variables": [{
            "label": "R",
            "attachment_points": sf("R-S", "R is attached to sulfur", evidence),
            "definition": sf("RSH / SR in the source schemes", "aryl, heteroaryl, or alkyl substituent", evidence),
            "allowed_values": sf("1a-1l", ["aryl", "heteroaryl", "alkyl"], evidence),
            "constraints": sf("Scope examples 1a-1l", "Values are restricted to source-shown examples; no unreported substituents inferred.", evidence),
        }],
        "instance_mapping": [],
        "reason": sf("Generic RSH/SR drawing", "Source uses an unresolved R variable.", evidence),
    }
    m["stereochemistry"]["chirality_status"] = sf("Generic R group", "unresolved", evidence)
    if product:
        m["stereochemistry"]["double_bonds"] = sf(">20:1 E/Z generic product relationship", [{"configuration": product_geometry, "location": "enoate C=C"}], evidence)
        m["stereochemistry"]["geometry_status"] = sf(f"{product_geometry} product identity", "specified", evidence)
        m["stereochemistry"]["representation_status"] = sf("Generic template cannot encode every instance", "not_encoded", evidence)
    else:
        m["stereochemistry"]["geometry_status"] = sf("Generic RSH", "unresolved", evidence)
        m["stereochemistry"]["representation_status"] = sf("Generic template", "not_encoded", evidence)
    m["identity_note"] = sf("Generic source template", "Not a unique concrete molecule.", evidence)
    return m


def participant(pid, mid, role, evidence, stage="reaction", product_role=None):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    if product_role:
        p["product_role"] = sf(product_role, product_role, evidence)
    return p


def charge(p, key, raw, normalized, evidence, reference=None):
    p["charges"][key] = sf(raw, normalized, evidence)
    if reference:
        p["equivalence_reference"] = sf(raw, reference, evidence)


def step(step_id, evidence, stage="reaction", pids=(), operation=None,
         temperature=None, time=None, atmosphere=None):
    s = st.new_step(step_id, evidence, stage=stage)
    s["participant_ids"] = list(pids)
    if operation is not None:
        s["operation"] = sf(operation, operation, evidence)
    if temperature is not None:
        s["temperature"] = sf(temperature[0], temperature[1], evidence)
    if time is not None:
        s["time"] = sf(time[0], time[1], evidence)
    if atmosphere is not None:
        s["atmosphere"] = sf(atmosphere, atmosphere, evidence)
    return s


def outcome(oid, kind, target, raw, normalized, evidence, product_role=None,
            yield_type=None, basis=None, scope=None, value_fact=None):
    targets = target if isinstance(target, list) else [target] if target else []
    o = st.new_outcome(oid, kind, targets, evidence)
    o["value"] = value_fact or sf(raw, normalized, evidence)
    if product_role:
        o["product_role"] = sf(product_role, product_role, evidence)
    if yield_type:
        o["yield_type"] = sf(yield_type, yield_type, evidence)
    if basis:
        o["measurement_basis"] = sf(basis, basis, evidence)
    if scope:
        o["scope"] = sf(scope, scope, evidence)
    return o


def current_classification(r, evidence):
    changes = [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "S"}],
         "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}],
         "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
        {"change_id": "b3", "atoms": [{"atom_map": None, "element": "S"}, {"atom_map": None, "element": "H"}],
         "before_order": 1, "after_order": 0, "change_type": "hydrogen_loss", "hydrogen_count": 1},
        {"change_id": "b4", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
         "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
    ]
    labels = [
        {"category_number": 4, "code": "C-S", "supporting_change_ids": ["b1"]},
        {"category_number": 1, "code": "C-H", "supporting_change_ids": ["b2"]},
    ]
    deriv = "Net thia-Michael transformation inferred from the source allenoate/thiol and thioenoate drawings; no complete atom mapping is asserted."
    r["classification"]["bond_changes"] = [df({"source": "reaction drawing"}, c, evidence, deriv) for c in changes]
    r["classification"]["labels"] = df({"source": "reaction drawing"}, labels, evidence, deriv)
    r["classification"]["coverage"] = inf("Partial source-drawing classification", "partial", evidence, deriv)
    r["classification"]["limitations"] = [
        "No complete atom mapping was asserted.",
        "Hydrogen transfer is a net change; an elementary proton-transfer mechanism is not claimed.",
        "Spectator and work-up species are outside the bond-change classification.",
    ]


def base_reaction(rid, evidence, record_type="specific_experiment", categories=None,
                  location=None, name=None):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    if categories is not None:
        r["experiment_categories"] = sf(categories, categories, evidence)
    if location:
        r["location"] = sf(location, location, evidence)
    if name:
        r["reaction_name"] = sf(name, name, evidence)
    r["completeness"] = "partial"
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    return r


def set_common_scope(r, thiol_mid, product_mid, z_product_mid, evidence, tea_molpct=20,
                     thiol_mass=None, thiol_volume=None):
    pa = participant("allenoate", "mol_benzyl_allenoate", "reactant", evidence)
    charge(pa, "volume", "123 uL", {"number": 123, "unit": "uL"}, [EV[5]])
    charge(pa, "amount", "0.75 mmol", {"number": 0.75, "unit": "mmol"}, [EV[5]])
    charge(pa, "equivalents", "150 mol%", {"number": 1.5, "unit": "equiv"}, evidence, thiol_mid)
    pt = participant("thiol", thiol_mid, "reactant", evidence)
    charge(pt, "amount", "0.5 mmol", {"number": 0.5, "unit": "mmol"}, [EV[5]])
    charge(pt, "equivalents", "100 mol%", {"number": 1.0, "unit": "equiv"}, [EV[5]], thiol_mid)
    if thiol_mass:
        charge(pt, "mass", thiol_mass, {"number": float(thiol_mass.split()[0]), "unit": "mg"}, evidence)
    if thiol_volume:
        charge(pt, "volume", thiol_volume, {"number": float(thiol_volume.split()[0]), "unit": "uL"}, evidence)
    pe = participant("triethylamine", "mol_triethylamine", "catalyst", evidence)
    charge(pe, "mol_percent", f"{tea_molpct} mol%", {"number": tea_molpct, "unit": "mol%"}, evidence, thiol_mid)
    if tea_molpct == 20:
        charge(pe, "volume", "19.2 uL", {"number": 19.2, "unit": "uL"}, [EV[5]])
        charge(pe, "amount", "0.1 mmol", {"number": 0.1, "unit": "mmol"}, [EV[5]])
    ps = participant("thf", "mol_thf", "solvent", evidence)
    charge(ps, "volume", "1 mL", {"number": 1, "unit": "mL"}, [EV[5]])
    pp = participant("product", product_mid, "product", evidence, product_role="main_product")
    pz = participant("z_product", z_product_mid, "product", evidence, product_role="minor_geometric_isomer")
    psg = participant("silica", "mol_silica_gel", "purification_agent", [EV[5]], stage="purification")
    r["participants"] = [pa, pt, pe, ps, pp, pz, psg]
    r["reaction_scale"] = sf("nucleophile 0.5 mmol (100 mol%)",
                             {"reference_molecule_id": thiol_mid, "amount": {"number": 0.5, "unit": "mmol"}}, [EV[5]])
    r["solvent_conditions"]["composition"] = sf("freshly distilled THF", ["THF"], [EV[5]])
    r["solvent_conditions"]["total_volume"] = sf("THF (1 mL)", {"number": 1, "unit": "mL"}, [EV[5]])
    r["solvent_conditions"]["concentration"] = sf("THF (1 mL, 0.5 M)", {"number": 0.5, "unit": "M"}, [EV[5]])
    r["solvent_conditions"]["concentration_basis"] = sf("nucleophile 0.5 mmol in 1 mL THF", thiol_mid, [EV[5]])
    pids = ["allenoate", "thiol", "triethylamine", "thf", "product", "z_product"]
    r["steps"] = [
        step("charge_and_purge", [EV[5]], pids=pids,
             operation="Charge oven-dried sealed tube in the reported order; add freshly distilled THF and purge with nitrogen.",
             atmosphere="nitrogen"),
        step("heat_and_stir", evidence, pids=pids,
             operation="Place in an oil bath and stir.",
             temperature=("70 °C", {"number": 70, "unit": "°C"}),
             time=("24 h", {"number": 24, "unit": "h"}), atmosphere="nitrogen"),
        step("evaporate_and_chromatograph", [EV[5]], stage="purification", pids=["silica", "product", "z_product"],
             operation="Evaporate solvent over dry silica gel and purify by flash column chromatography."),
    ]
    r["special_conditions"]["reactor"] = sf("13 x 100 mm oven-dried sealed test tube", "sealed test tube", [EV[5]])
    r["procedure"]["original_text"] = sf(
        "Oven-dried sealed tube was charged with benzyl buta-2,3-dienoate (123 uL, 0.75 mmol, 150 mol%), distilled triethylamine (19.2 uL, 0.1 mmol, 20 mol%) and nucleophile 1a-s (0.5 mmol, 100 mol%). Freshly distilled THF (1 mL, 0.5 M) was added and purged with nitrogen. The mixture was placed in an oil bath at 70 °C and was allowed to stir for 24 h. Solvent was evaporated over dry silica gel and purification by flash column chromatography gave the desired product 2a-s.",
        "Section 4.2 general procedure", [EV[5]])
    r["procedure"]["general_procedure_reference"] = sf("Section 4.2", "Section 4.2", [EV[5]])
    r["procedure"]["addition_order"] = sf(
        "allenoate, triethylamine, nucleophile; then THF; purge with nitrogen",
        ["benzyl allenoate", "triethylamine", "thiol", "THF", "nitrogen purge"], [EV[5]])
    r["procedure"]["pretreatment"] = sf(
        "oven-dried sealed tube; freshly distilled THF; distilled triethylamine",
        ["oven-dried reactor", "freshly distilled THF", "distilled triethylamine"], [EV[5]])
    if tea_molpct != 20:
        r["procedure"]["modifications"] = sf(
            f"Reaction with {tea_molpct} mol% triethylamine", {"triethylamine_mol_percent": tea_molpct}, evidence)
        r["notes"].append(inf(
            "General procedure lists 20 mol% triethylamine, while the product-specific/table condition is 150 mol%.",
            "The 150 mol% condition is treated as an explicit modification; the general-procedure 19.2 uL/0.1 mmol charge is not inherited.", evidence,
            "The Scheme 2 footnote and Table 1 entry 6 explicitly identify 150 mol% for 2b."))


def scope_reaction(code, thiol_mid, product_mid, y, specific_page=None, tea=20,
                   thiol_mass=None, thiol_volume=None, conflict_yield=False):
    evidence = [EV[4], IMG[4], EV[5]]
    if specific_page:
        pages = [specific_page] if isinstance(specific_page, int) else list(specific_page)
        for page in pages:
            evidence += [EV[page], IMG[page]]
    if code == "2a": evidence += [EV[2], EV[3], IMG[3]]
    if code == "2b": evidence += [EV[2], EV[3], IMG[3]]
    if code in ("2c", "2d", "2e", "2f"): evidence += [EV[3]]
    evidence = list(dict.fromkeys(evidence))
    r = base_reaction(f"scope_{code}", evidence, categories=["substrate_scope", "optimization"] if code in ("2a", "2b") else ["substrate_scope"],
                      location={"pages": sorted({int(x.split(":page")[1].split(":")[0]) for x in evidence}), "scheme": "Scheme 2", "compound": code},
                      name="thia-Michael addition")
    r["general_scheme_ids"] = [f"{DOC}:generic_thia_michael"]
    r["transformation_type"] = inf("allenoate + thiol -> E beta-thioenoate", "thia-Michael addition", evidence,
                                   "Transformation type inferred from the source reaction name and structures.")
    z_product_mid = product_mid + "_z"
    set_common_scope(r, thiol_mid, product_mid, z_product_mid, evidence, tea, thiol_mass, thiol_volume)
    if conflict_yield:
        alt1 = sf("99%", {"number": 99, "unit": "%"}, [EV[4], IMG[4]])
        alt2 = sf("quantitative yield", "quantitative", [EV[5], EV[6]])
        vf = sf(["99%", "quantitative yield"], None, [EV[4], IMG[4], EV[5], EV[6]], status="conflict", alternatives=[alt1, alt2])
        r["outcomes"].append(outcome("yield", "yield", product_mid, None, None, evidence,
                                     product_role="main_product", yield_type="isolated", basis="isolated yield", value_fact=vf))
        r["review"]["unresolved"].append("Product 2h is shown as 99% in Scheme 2 but described as quantitative in prose/Section 4.10; both are retained without choosing.")
    else:
        r["outcomes"].append(outcome("yield", "yield", product_mid, f"{y}%", {"number": y, "unit": "%"}, evidence,
                                     product_role="main_product", yield_type="isolated", basis="isolated yield"))
    ratio_basis = "NMR analysis of the crude reaction mixture" if code == "2b" else None
    r["outcomes"].append(outcome(
        "ez_selectivity", "selectivity", [product_mid, z_product_mid], ">20:1 E/Z",
        {"ratio": [20, 1], "qualifier": "gt", "order": ["E", "Z"], "raw_text": ">20:1"}, evidence,
        product_role="main_product", basis=ratio_basis,
        scope={"metric": "E/Z", "order": [product_mid, z_product_mid], "Z_identity_note": "ratio counterpart; not reported as isolated"}))
    if not ratio_basis:
        r["review"]["unresolved"].append("The source reports >20:1 E/Z but does not state a measurement method for this scope entry.")
    r["review"]["unresolved"].extend([
        "No explicitly associated SI file is available; detailed chromatography parameters and complete reproducibility details remain unavailable.",
        "Candidate structures and partial bond changes require expert review.",
    ])
    current_classification(r, evidence)
    return r


def optimization_reaction(entry, thiol_mid, product_mid, tea, y, minor=False):
    evidence = [EV[3], IMG[3]]
    if entry in (1, 3): evidence.append(EV[2])
    if minor: evidence += [EV[4], IMG[4]]
    r = base_reaction(f"table1_entry_{entry}", evidence, categories=["optimization"],
                      location={"page": 3, "table": "Table 1", "entry": entry}, name="thia-Michael addition")
    r["general_scheme_ids"] = [f"{DOC}:generic_thia_michael"]
    r["transformation_type"] = inf("allenoate + thiol -> E beta-thioenoate", "thia-Michael addition", evidence,
                                   "Transformation type inferred from the source name and scheme.")
    pa = participant("allenoate", "mol_benzyl_allenoate", "reactant", evidence)
    charge(pa, "equivalents", "150 mol%", {"number": 1.5, "unit": "equiv"}, evidence, thiol_mid)
    pt = participant("thiol", thiol_mid, "reactant", evidence)
    charge(pt, "equivalents", "RSH basis", {"number": 1.0, "unit": "equiv"}, evidence, thiol_mid)
    parts = [pa, pt]
    pids = ["allenoate", "thiol"]
    if tea is not None:
        pe = participant("triethylamine", "mol_triethylamine", "catalyst", evidence)
        charge(pe, "mol_percent", f"{tea} mol%", {"number": tea, "unit": "mol%"}, evidence, thiol_mid)
        parts.append(pe); pids.append("triethylamine")
    ps = participant("thf", "mol_thf", "solvent", evidence)
    parts.append(ps); pids.append("thf")
    pp = participant("product", product_mid, "product", evidence, product_role="main_product")
    parts.append(pp); pids.append("product")
    z_product_mid = product_mid + "_z"
    parts.append(participant("z_product", z_product_mid, "product", evidence, product_role="minor_geometric_isomer")); pids.append("z_product")
    if minor:
        parts.append(participant("minor_3b", "mol_3b", "product", evidence, product_role="minor_product"))
        pids.append("minor_3b")
    r["participants"] = parts
    r["solvent_conditions"]["composition"] = sf("THF", ["THF"], evidence)
    r["solvent_conditions"]["concentration"] = sf("THF (0.5 M)", {"number": 0.5, "unit": "M"}, evidence)
    r["steps"] = [step("reaction", evidence, pids=pids, operation="Stir under the Table 1 conditions.",
                       temperature=("70 °C", {"number": 70, "unit": "°C"}),
                       time=("24 h", {"number": 24, "unit": "h"}))]
    r["outcomes"] = [
        outcome("yield", "yield", product_mid, f"{y}%", {"number": y, "unit": "%"}, evidence,
                product_role="main_product", yield_type="isolated", basis="isolated"),
        outcome("ez_selectivity", "selectivity", [product_mid, z_product_mid], ">20:1 E/Z",
                {"ratio": [20, 1], "qualifier": "gt", "order": ["E", "Z"], "raw_text": ">20:1"}, evidence,
                product_role="main_product", basis="NMR analysis of crude reaction mixture" if product_mid == "mol_2b" else None,
                scope={"metric": "E/Z", "order": [product_mid, z_product_mid], "Z_identity_note": "ratio counterpart; not reported as isolated"}),
    ]
    if minor:
        r["outcomes"].append(outcome(
            "minor_3b_observation", "trace", "mol_3b", "minor amount of physically different crystals",
            "minor amount; not isolated in quantity", evidence, product_role="minor_product",
            basis="single-crystal X-ray analysis", scope="Observed only in the no-triethylamine experiment; isolation and spectroscopic characterization unsuccessful."))
        r["notes"].append(sf(
            "there was low conversion", "low conversion (no numeric value reported)", [EV[3]]))
    r["review"]["unresolved"] = [
        "Optimization-table scale, complete charge details, work-up and purification are not reported in the covered pages.",
        "Candidate structures and partial bond changes require expert review.",
    ]
    if product_mid != "mol_2b":
        r["review"]["unresolved"].append("The >20:1 E/Z measurement method is not stated for this entry.")
    current_classification(r, evidence)
    return r


def literature_record(rid, evidence, reactants, products, outcomes, note, additive=None):
    r = base_reaction(rid, evidence, record_type="literature_example", categories=["other"],
                      location={"page": 3, "scheme": "Scheme 1", "context": "Previous work"},
                      name="previous-work thia-Michael addition")
    r["transformation_type"] = inf("Previous work schematic", "thia-Michael addition", evidence,
                                   "Classified from the source heading; not a current-paper experiment.")
    parts = []
    for i, mid in enumerate(reactants): parts.append(participant(f"reactant_{i+1}", mid, "reactant", evidence))
    for i, mid in enumerate(products): parts.append(participant(f"product_{i+1}", mid, "product", evidence, product_role="reported_product"))
    if additive:
        p = participant("additive", additive[0], "catalyst", evidence)
        if additive[1] is not None:
            charge(p, "mol_percent", f"{additive[1]} mol%", {"number": additive[1], "unit": "mol%"}, evidence, reactants[0])
        parts.append(p)
    r["participants"] = parts
    r["outcomes"] = []
    for i, (mid, raw, norm, extra) in enumerate(outcomes):
        r["outcomes"].append(outcome(f"reported_outcome_{i+1}", extra.get("kind", "yield"), mid, raw, norm, evidence,
                                     product_role="reported_product", yield_type=extra.get("yield_type"),
                                     basis=extra.get("basis"), scope=extra.get("scope")))
    r["notes"].append(sf(note, note, evidence))
    r["review"]["unresolved"] = [
        "This is explicitly previous work, not an experiment executed in the current paper.",
        "The background schematic omits complete procedure and measurement-basis details; exact product structures remain unresolved where not safely encoded.",
    ]
    return r


def main():
    data = json.loads((RUN / "dataset.json").read_text())
    # The submission skeleton is produced by `standard init`; this script only
    # fills its molecule/reaction/source-context arrays from page evidence.
    bundle = json.loads((RUN / "raw_model/init_pages_0007_0008.json").read_text())
    # Concrete current-paper molecules.
    mols = []
    mols.append(molecule("mol_benzyl_allenoate", "benzyl buta-2,3-dienoate", "benzyl buta-2,3-dienoate", [EV[2], EV[3], IMG[3], EV[4], IMG[4], EV[5]], "O=C(OCc1ccccc1)C=C=C"))
    mols.append(molecule("mol_triethylamine", "Et3N", "triethylamine", [EV[2], EV[3], IMG[3], EV[4], IMG[4], EV[5]], "CCN(CC)CC"))
    mols.append(molecule("mol_thf", "THF", "tetrahydrofuran", [EV[3], IMG[3], EV[4], IMG[4], EV[5]], "C1CCOC1"))
    mols.append(molecule("mol_silica_gel", "silica gel", "silica gel", [EV[5]], None))
    thiols = {
        "1a": ("mol_1a", "thiophenol", "Sc1ccccc1", [EV[2], EV[3], IMG[3], EV[4], IMG[4], EV[5]]),
        "1b": ("mol_1b", "2-aminothiophenol", "Nc1ccccc1S", [EV[2], EV[3], IMG[3], EV[4], IMG[4], EV[5]]),
        "1c": ("mol_1c", "2-amino-4-chlorothiophenol", "Sc1c(N)cc(Cl)cc1", [EV[3], EV[4], IMG[4], EV[6]]),
        "1d": ("mol_1d", "2-bromothiophenol", "Sc1ccccc1Br", [EV[3], EV[4], IMG[4], EV[6]]),
        "1e": ("mol_1e", "2,6-dichlorothiophenol", "Sc1c(Cl)cccc1Cl", [EV[3], EV[4], IMG[4], EV[6]]),
        "1f": ("mol_1f", "2-naphthalenethiol", "Sc1ccc2ccccc2c1", [EV[3], EV[4], IMG[4], EV[6]]),
        "1g": ("mol_1g", "5-chlorobenzo[d]thiazole-2-thiol", "Sc1nc2ccc(Cl)cc2s1", [EV[4], IMG[4], EV[5], EV[6]]),
        "1h": ("mol_1h", "2-pyridinethiol", "Sc1ccccn1", [EV[4], IMG[4], EV[5], EV[6]]),
        "1i": ("mol_1i", "5-(4-nitrophenyl)-1,3,4-oxadiazole-2-thiol", "Sc1nnc(-c2ccc([N+](=O)[O-])cc2)o1", [EV[4], IMG[4], EV[5], EV[6]]),
        "1j": ("mol_1j", "2-mercaptoethanol", "OCCS", [EV[4], IMG[4], EV[5], EV[6], EV[7], IMG[7]]),
        "1k": ("mol_1k", "tert-butylthiol", "CC(C)(C)S", [EV[4], IMG[4], EV[5], EV[7], IMG[7]]),
        "1l": ("mol_1l", "1-dodecanethiol", "CCCCCCCCCCCCS", [EV[4], IMG[4], EV[5], EV[7], IMG[7]]),
    }
    for label, (mid, name, smi, ev) in thiols.items(): mols.append(molecule(mid, label, name, ev, smi))
    product_specs = {
        "2a": ("mol_2a", "benzyl (E)-3-(phenylthio)but-2-enoate", "O=C(OCc1ccccc1)/C=C(/Sc1ccccc1)C", [EV[2], EV[3], IMG[3], EV[4], IMG[4], EV[5]]),
        "2b": ("mol_2b", "benzyl (E)-3-(2-aminophenylthio)but-2-enoate", "O=C(OCc1ccccc1)/C=C(/Sc1ccccc1N)C", [EV[2], EV[3], IMG[3], EV[4], IMG[4], EV[5]]),
        "2c": ("mol_2c", "benzyl (E)-3-(2-amino-4-chlorophenylthio)but-2-enoate", "O=C(OCc1ccccc1)/C=C(/Sc1c(N)cc(Cl)cc1)C", [EV[3], EV[4], IMG[4], EV[6]]),
        "2d": ("mol_2d", "benzyl (E)-3-(2-bromophenylthio)but-2-enoate", "O=C(OCc1ccccc1)/C=C(/Sc1ccccc1Br)C", [EV[3], EV[4], IMG[4], EV[5], EV[6]]),
        "2e": ("mol_2e", "benzyl (E)-3-(2,6-dichlorophenylthio)but-2-enoate", "O=C(OCc1ccccc1)/C=C(/Sc1c(Cl)cccc1Cl)C", [EV[3], EV[4], IMG[4], EV[5], EV[6]]),
        "2f": ("mol_2f", "benzyl (E)-3-(naphthalen-2-ylthio)but-2-enoate", "O=C(OCc1ccccc1)/C=C(/Sc1ccc2ccccc2c1)C", [EV[3], EV[4], IMG[4], EV[5], EV[6]]),
        "2g": ("mol_2g", "benzyl (E)-3-(5-chlorobenzo[d]thiazol-2-ylthio)but-2-enoate", "O=C(OCc1ccccc1)/C=C(/Sc1nc2ccc(Cl)cc2s1)C", [EV[4], IMG[4], EV[5], EV[6]]),
        "2h": ("mol_2h", "benzyl (E)-3-(pyridin-2-ylthio)but-2-enoate", "O=C(OCc1ccccc1)/C=C(/Sc1ccccn1)C", [EV[4], IMG[4], EV[5], EV[6]]),
        "2i": ("mol_2i", "benzyl (E)-3-(5-(4-nitrophenyl)-1,3,4-oxadiazol-2-ylthio)but-2-enoate", "O=C(OCc1ccccc1)/C=C(/Sc1nnc(-c2ccc([N+](=O)[O-])cc2)o1)C", [EV[4], IMG[4], EV[5], EV[6]]),
        "2j": ("mol_2j", "benzyl (E)-3-(2-hydroxyethylthio)but-2-enoate", "O=C(OCc1ccccc1)/C=C(/SCCO)C", [EV[4], IMG[4], EV[5], EV[6], EV[7], IMG[7]]),
        "2k": ("mol_2k", "benzyl (E)-3-(tert-butylthio)but-2-enoate", "O=C(OCc1ccccc1)/C=C(/SC(C)(C)C)C", [EV[4], IMG[4], EV[5], EV[7], IMG[7]]),
        "2l": ("mol_2l", "benzyl (E)-3-(dodecylthio)but-2-enoate", "O=C(OCc1ccccc1)/C=C(/SCCCCCCCCCCCC)C", [EV[4], IMG[4], EV[5], EV[7], IMG[7]]),
    }
    for label, (mid, name, smi, ev) in product_specs.items():
        mols.append(molecule(mid, label, name, ev, smi, geometry="E"))
        z_smi = smi.replace("C=C(/S", "C=C(\\S")
        zmid = mid + "_z"
        zm = molecule(zmid, f"Z-{label} ratio counterpart", f"Z geometric isomer corresponding to {label}", ev, z_smi, geometry="Z")
        zm["identity_status"] = "unresolved"
        zm["occurrences"][0]["compound_label"] = inf(
            ">20:1 E/Z", f"Z-{label} ratio counterpart", ev,
            "A Z denominator is explicitly implied by the reported E/Z metric, but it is not reported as isolated or characterized.")
        zm["names"]["original"] = inf(
            ">20:1 E/Z", f"Z geometric isomer corresponding to {label}", ev,
            "Machine identity created only to keep the source E and Z ratio targets distinct.")
        zm["names"]["english"] = zm["names"]["original"]
        zm["identity_note"] = inf(
            "Z isomer in >20:1 E/Z ratio", "Ratio-counterpart identity; not reported as isolated, characterized, or detected above a stated limit.", ev,
            "The Z structure is the constitutional counterpart of the named E product with alkene geometry inverted; this requires expert review.")
        mols.append(zm)
    mols.append(molecule("mol_3b", "3b", "unconjugated thia-Michael product 3b", [EV[3], EV[4], IMG[4]], "C=C(Sc1ccccc1N)CC(=O)OCc1ccccc1"))
    mols += [generic_molecule("mol_generic_rsh", "RSH", "generic thiol RSH", [EV[2], EV[3], IMG[3], EV[4], IMG[4]]),
             generic_molecule("mol_generic_product", "E beta-thioenoate", "generic E beta-thioenoate", [EV[2], EV[3], IMG[3], EV[4], IMG[4]], product=True, product_geometry="E"),
             generic_molecule("mol_generic_product_z", "Z beta-thioenoate ratio counterpart", "generic Z beta-thioenoate ratio counterpart", [EV[2], EV[3], IMG[3], EV[4], IMG[4]], product=True, product_geometry="Z")]
    # Background-only identities from Scheme 1 are deliberately unresolved where drawings/names are insufficient.
    mols += [
        molecule("mol_prior_eth_allenoate", "ethyl 2-methylbuta-2,3-dienoate", "ethyl 2-methylbuta-2,3-dienoate", [EV[3], IMG[3]], "CCOC(=O)C(C)=C=C"),
        molecule("mol_prior_4cl_thiol", "4-ClC6H4SH", "4-chlorothiophenol", [EV[3], IMG[3]], "Sc1ccc(Cl)cc1"),
        molecule("mol_prior_conjugated", "prior conjugated product", "conjugated thia-Michael product shown in Scheme 1", [EV[3], IMG[3]], None),
        molecule("mol_prior_nonconjugated", "prior non-conjugated product", "non-conjugated thia-Michael product shown in Scheme 1", [EV[3], IMG[3]], None),
        molecule("mol_pbu3", "PBu3", "tributylphosphine", [EV[3], IMG[3]], "CCCCP(CCCC)CCCC"),
        molecule("mol_methyl_2_butynoate", "methyl 2-butynoate", "methyl but-2-ynoate", [EV[2], EV[3], IMG[3]], "COC(=O)C#CC"),
        molecule("mol_phsna", "PhSNa", "sodium thiophenolate", [EV[3], IMG[3]], "[Na+].[S-]c1ccccc1"),
        molecule("mol_prior_alkynoate_z", "Ref 6b product", "Z-major beta-phenylthioenoate product mixture", [EV[3], IMG[3]], None),
        molecule("mol_prior_alkynoate_ez", "Ref 6c product", "1:1 E/Z beta-phenylthioenoate product mixture", [EV[2], EV[3], IMG[3]], None),
    ]
    bundle["molecules"] = mols

    reactions = []
    # Generic current-paper transformation (graphical abstract + Scheme 1 + Scheme 2).
    ge = [EV[2], IMG[2], EV[3], IMG[3], EV[4], IMG[4]]
    g = base_reaction("generic_thia_michael", ge, record_type="generic_scheme", categories=["other"],
                      location={"pages": [2, 3, 4], "regions": ["graphical abstract", "Scheme 1 this work", "Scheme 2 header"]},
                      name="E- and chemoselective thia-Michael addition")
    g["transformation_type"] = inf("benzyl allenoate + RSH -> E beta-thioenoate", "thia-Michael addition", ge,
                                   "Generic transformation inferred from the source reaction drawings and title.")
    g["participants"] = [participant("allenoate", "mol_benzyl_allenoate", "reactant", ge),
                         participant("rsh", "mol_generic_rsh", "reactant", ge),
                         participant("triethylamine", "mol_triethylamine", "catalyst", ge),
                         participant("thf", "mol_thf", "solvent", ge),
                         participant("product", "mol_generic_product", "product", ge, product_role="generic_main_product"),
                         participant("z_product", "mol_generic_product_z", "product", ge, product_role="generic_minor_geometric_isomer")]
    charge(g["participants"][2], "mol_percent", "20 mol%", {"number": 20, "unit": "mol%"}, ge, "mol_generic_rsh")
    g["solvent_conditions"]["composition"] = sf("THF", ["THF"], ge)
    g["solvent_conditions"]["concentration"] = sf("0.5 M", {"number": 0.5, "unit": "M"}, ge)
    g["steps"] = [step("generic_reaction", ge, pids=["allenoate", "rsh", "triethylamine", "thf", "product", "z_product"],
                       temperature=("70 °C", {"number": 70, "unit": "°C"}), time=("24 h", {"number": 24, "unit": "h"}))]
    g["outcomes"] = [outcome("generic_ez", "selectivity", ["mol_generic_product", "mol_generic_product_z"], ">20:1 E/Z",
                              {"ratio": [20, 1], "qualifier": "gt", "order": ["E", "Z"], "raw_text": ">20:1"}, ge,
                              product_role="generic_product_pair", scope={"metric": "E/Z", "order": ["mol_generic_product", "mol_generic_product_z"], "examples": 12})]
    g["review"]["unresolved"] = ["Generic R does not identify one concrete experiment; individual examples are separate records.",
                                   "The generic E/Z measurement basis is not reported."]
    current_classification(g, ge)
    reactions.append(g)

    # Scope records; 2a and 2b also represent Table 1 entries 2 and 6, respectively.
    scope = [
        ("2a", "mol_1a", "mol_2a", 72, 5, 20, None, "51 uL", False),
        ("2b", "mol_1b", "mol_2b", 92, 5, 150, "62.5 mg", None, False),
        ("2c", "mol_1c", "mol_2c", 93, 6, 20, "100.8 mg", None, False),
        ("2d", "mol_1d", "mol_2d", 84, 6, 20, "94.5 mg", None, False),
        ("2e", "mol_1e", "mol_2e", 82, 6, 20, "89.5 mg", None, False),
        ("2f", "mol_1f", "mol_2f", 68, 6, 20, "80.1 mg", None, False),
        ("2g", "mol_1g", "mol_2g", 83, 6, 20, "100.87 mg", None, False),
        ("2h", "mol_1h", "mol_2h", 99, 6, 20, "55.58 mg", None, True),
        ("2i", "mol_1i", "mol_2i", 75, 6, 20, "111.60 mg", None, False),
        ("2j", "mol_1j", "mol_2j", 92, (6, 7), 20, "39 mg", None, False),
        ("2k", "mol_1k", "mol_2k", 68, 7, 20, "45 mg", None, False),
        ("2l", "mol_1l", "mol_2l", 82, 7, 20, "101.2 mg", None, False),
    ]
    for args in scope:
        reactions.append(scope_reaction(args[0], args[1], args[2], args[3], specific_page=args[4], tea=args[5],
                                        thiol_mass=args[6], thiol_volume=args[7], conflict_yield=args[8]))
    reactions += [
        optimization_reaction(1, "mol_1a", "mol_2a", 150, 60),
        optimization_reaction(3, "mol_1b", "mol_2b", 20, 68),
        optimization_reaction(4, "mol_1b", "mol_2b", 10, 60),
        optimization_reaction(5, "mol_1b", "mol_2b", None, 26, minor=True),
    ]
    # Explicit previous-work rows/arrows from Scheme 1.
    le = [EV[3], IMG[3]]
    reactions += [
        literature_record("literature_ref7_no_additive", le,
                          ["mol_prior_eth_allenoate", "mol_prior_4cl_thiol"],
                          ["mol_prior_conjugated", "mol_prior_nonconjugated"],
                          [("mol_prior_conjugated", "30%", {"number": 30, "unit": "%"}, {"yield_type": None}),
                           ("mol_prior_nonconjugated", "14%", {"number": 14, "unit": "%"}, {"yield_type": None})],
                          "Previous work, Ref. 7; no additive."),
        literature_record("literature_ref7_pbu3", le,
                          ["mol_prior_eth_allenoate", "mol_prior_4cl_thiol"],
                          ["mol_prior_nonconjugated"],
                          [("mol_prior_nonconjugated", "76%", {"number": 76, "unit": "%"}, {"yield_type": None})],
                          "Previous work, Ref. 7; PBu3 gives only the non-conjugated product in the source discussion.",
                          additive=("mol_pbu3", 20)),
        literature_record("literature_ref6b_alkynoate", le,
                          ["mol_methyl_2_butynoate", "mol_phsna"], ["mol_prior_alkynoate_z"],
                          [("mol_prior_alkynoate_z", "77%, Z major", {"number": 77, "unit": "%"}, {"scope": {"geometry": "Z major"}})],
                          "Previous work, Ref. 6b; Z-major product."),
        literature_record("literature_ref6c_alkynoate", [EV[2], EV[3], IMG[3]],
                          ["mol_methyl_2_butynoate", "mol_phsna"], ["mol_prior_alkynoate_ez"],
                          [("mol_prior_alkynoate_ez", "85%", {"number": 85, "unit": "%"}, {}),
                           ("mol_prior_alkynoate_ez", "E/Z 1:1", {"ratio": [1, 1], "order": ["E", "Z"]}, {"kind": "selectivity", "scope": {"metric": "E/Z"}})],
                          "Previous work, Ref. 6c; the prose specifies sodium hydroxide whereas the compact Scheme 1 arrow is labelled PhSNa, so complete conditions are unresolved."),
    ]
    # Current substrate preparation is mentioned but cannot be reconstructed without the unprovided SI.
    ue = [EV[5], IMG[5]]
    u = base_reaction("substrate_benzyl_allenoate_preparation_unresolved", ue, record_type="unresolved", categories=["other"],
                      location={"page": 5, "section": "4.1"}, name="Synthesis of substrate benzyl buta-2,3-dienoate")
    u["participants"] = [participant("product", "mol_benzyl_allenoate", "product", ue, product_role="reported_substrate_product")]
    u["procedure"]["original_text"] = sf(
        "The benzyl buta-2,3-dienoate was synthesize in a slightly modified procedure to the previously reported method. See the Supplemental Materials for the detailed experimental procedure.",
        "Current-paper preparation stated; detailed procedure unavailable.", ue)
    u["review"]["unresolved"] = [
        "Reactants, modifications, conditions, work-up, purification and outcome are unavailable because no SI file was explicitly supplied or associated.",
        "This record is unresolved and is not treated as a reproducible specific experiment.",
    ]
    reactions.append(u)
    bundle["reactions"] = reactions
    bundle["source_context"] = [
        {"context_id": f"{DOC}:bibliography", "kind": "bibliographic_metadata", "evidence_ids": [EV[1], IMG[1], EV[2], IMG[2]],
         "title": "E- and chemoselective thia-Michael addition to benzyl allenoate", "doi": "10.1080/10426507.2020.1799365"},
        {"context_id": f"{DOC}:figure1_xray", "kind": "structural_characterization", "evidence_ids": [EV[3], EV[4], IMG[4]],
         "finding": "Single-crystal X-ray structures support the C-S bond and E geometry of 2b and identify minor unconjugated product 3b in the no-triethylamine experiment."},
        {"context_id": f"{DOC}:characterization_pages5_7", "kind": "characterization", "evidence_ids": [EV[5], IMG[5], EV[6], IMG[6], EV[7], IMG[7]],
         "finding": "NMR/HRMS/IR and physical descriptions for 2a-2l are characterization attached to the corresponding syntheses, not standalone reactions."},
        {"context_id": f"{DOC}:si_boundary", "kind": "source_boundary", "evidence_ids": [EV[1], EV[5], EV[7], IMG[7]],
         "finding": "The article links/mentions Supplemental Materials, but no SI file is explicitly present or associated in source_mapping.json; no SI-only detail was used."},
        {"context_id": f"{DOC}:references_pages7_8", "kind": "bibliographic_context", "evidence_ids": [EV[7], IMG[7], EV[8], IMG[8]],
         "finding": "The reference list begins on page 7 and occupies all of page 8. It contains no current-paper reaction arrow, table row, or experimental procedure."},
    ]

    # Provenance artifacts.
    stereo = {}
    for m in mols:
        smi = st.value(m["smiles"])
        if smi:
            stereo[m["molecule_id"]] = st.stereo_inspect(smi)
    (RUN / "raw_model/stereo_checks.json").write_text(json.dumps(stereo, ensure_ascii=False, indent=2) + "\n")
    page_map = {str(p): [r["reaction_id"] for r in reactions if any(e.startswith(f"{DOC}:page{p}:") for e in r["evidence_ids"])] for p in range(1, 9)}
    (RUN / "raw_model/record_ids_by_page.json").write_text(json.dumps(page_map, ensure_ascii=False, indent=2) + "\n")
    inventory = {
        "1": ["publisher cover, title, DOI, SI link; no page-1 reaction record"],
        "2": ["graphical abstract", "introduction prior-work prose", "results prose for Table 1 entries 1-3"],
        "3": ["Scheme 1 previous work/current generic scheme", "Table 1 entries 1-6", "2b/3b discussion", "scope prose 2c-2f"],
        "4": ["Figure 1 X-ray structures 2b/3b", "Scheme 2 scope 2a-2l with footnotes"],
        "5": ["scope prose 2d-2l", "general experimental", "unresolved substrate synthesis 4.1", "general procedure 4.2", "specific 2a/2b"],
        "6": ["specific syntheses/characterization 2c-2j (2j continues on page 7)"],
        "7": ["2j HRMS/IR continuation", "specific syntheses/characterization 2k and 2l", "supplementary-data statement", "disclosure/funding", "references 1-5"],
        "8": ["references 5-9 only; no current-paper reaction arrow, table, or procedure"],
    }
    (RUN / "raw_model/page_inventory.json").write_text(json.dumps(inventory, ensure_ascii=False, indent=2) + "\n")
    (RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"molecules": len(mols), "reactions": len(reactions), "page_map": page_map}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
