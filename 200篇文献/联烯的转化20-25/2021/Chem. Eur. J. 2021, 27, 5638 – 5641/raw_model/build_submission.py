#!/usr/bin/env python3
"""Build the source-derived v3.1 bundle for Chem. Eur. J. 2021, 27, 5638–5641.

This builder reads only the current prepared dataset/template and the decisions
transcribed from direct inspection of PDF pages 1–4. It does not read legacy
experiment records or a prior submitted bundle.
"""
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DOC = "doc_18eac4838543527a"
E = {p: [f"{DOC}:page{p}:image", f"{DOC}:page{p}:text"] for p in range(1, 5)}


def F(raw=None, normalized=None, ev=(), method="original_extraction", derivation=None,
      status=None):
    return st.fact(raw, normalized, ev, method=method, derivation=derivation,
                   status=status)


def missing(ev=(), status="unreported"):
    return F(evidence=ev, status=status) if False else st.fact(
        None, None, ev, status=status, verification="missing"
    )


def derived(normalized, ev, derivation, raw=None):
    return F(raw, normalized, ev, method="inference", derivation=derivation)


def occurrence(label, ev):
    return {
        "document_id": DOC,
        "compound_label": F(label, label, ev),
        "evidence_ids": list(ev),
    }


def unresolved_molecule(mid, label, ev, family="unresolved", extra_occurrences=()):
    """Concrete source drawing retained without unsafe structure transcription."""
    m = st.new_molecule(mid, DOC, ev)
    m["occurrences"] = [occurrence(label, ev)] + [occurrence(l, e) for l, e in extra_occurrences]
    m["identity_status"] = "source_identified"
    m["names"]["original"] = F(label, label, ev)
    m["names"]["english"] = missing(ev)
    m["names"]["abbreviations"] = missing(ev)
    m["names"]["aliases"] = missing(ev)
    m["cas_numbers"] = missing(ev)
    m["representation_type"] = "unresolved"
    m["smiles"] = missing(ev)
    m["generic_template"] = None
    if family == "vinylallene":
        m["stereochemistry"]["axial"] = derived(
            {"type": "allene axis", "assignment": "unspecified"}, ev,
            "The source drawing contains a substituted allene but gives no axial descriptor; exact substitution and configurational stability were not encoded."
        )
        chiral = "unspecified"
        geom = "unspecified"
        note = "Vinylallene drawing was visually located; no R_a/S_a, P/M, E/Z, wedge/dash, racemate, or enantiopurity statement is supplied."
    elif family in {"tropone", "cycloadduct"}:
        chiral = "unspecified"
        geom = "unspecified"
        note = "Product drawing has unsaturated and potentially stereochemical features but no explicit stereochemical assignment; no configuration was inferred."
    else:
        chiral = "unresolved"
        geom = "unresolved"
        note = "Compact source drawing was not transcribed to a complete stereochemical graph; stereochemistry remains unresolved."
    for key, value in (("chirality_status", chiral), ("geometry_status", geom),
                       ("representation_status", "not_encoded")):
        m["stereochemistry"][key] = derived(value, ev, note)
    for key in ("atom_centers", "double_bonds", "mixture", "assignment",
                "isomeric_smiles", "identity_key"):
        if key != "axial" or family != "vinylallene":
            m["stereochemistry"][key] = missing(ev)
    m["identity_note"] = F(note, note, ev)
    m["external_properties"] = []
    return m


def smiles_molecule(mid, label, name, smiles, ev, diagnostics):
    diag = st.stereo_inspect(smiles)
    diagnostics[mid] = diag
    m = st.new_molecule(mid, DOC, ev)
    m["occurrences"] = [occurrence(label, ev)]
    m["identity_status"] = "source_identified"
    m["names"]["original"] = F(label, label, ev)
    m["names"]["english"] = F(name, name, ev)
    m["names"]["abbreviations"] = F(label, [label], ev) if label != name else missing(ev)
    m["names"]["aliases"] = missing(ev)
    m["cas_numbers"] = missing(ev)
    m["representation_type"] = "smiles"
    basis = f"Candidate graph derived from the explicit source chemical name/abbreviation; RDKit canonicalized it to {diag['canonical_isomeric_smiles']}. Parsing is not source verification."
    m["smiles"] = F(label, diag["canonical_isomeric_smiles"], ev,
                    method="structure_derivation", derivation=basis)
    m["generic_template"] = None
    m["stereochemistry"]["atom_centers"] = derived(diag["atom_centers"], ev, basis)
    m["stereochemistry"]["double_bonds"] = derived(diag["double_bonds"], ev, basis)
    m["stereochemistry"]["axial"] = missing(ev)
    m["stereochemistry"]["mixture"] = missing(ev)
    m["stereochemistry"]["assignment"] = missing(ev)
    m["stereochemistry"]["chirality_status"] = derived("achiral", ev, basis)
    m["stereochemistry"]["geometry_status"] = derived("not_applicable", ev, basis)
    m["stereochemistry"]["representation_status"] = derived("encoded", ev, basis)
    m["stereochemistry"]["isomeric_smiles"] = F(
        label, diag["canonical_isomeric_smiles"], ev,
        method="structure_derivation", derivation=basis
    )
    m["stereochemistry"]["identity_key"] = derived(diag["stereo_graph_key"], ev, basis)
    m["identity_note"] = F(
        f"Source identifies {name}; candidate SMILES is derived and not visually structure-verified.",
        f"Source identifies {name}; candidate SMILES is derived and not visually structure-verified.", ev
    )
    m["external_properties"] = []
    return m


def vinylallene_smiles_molecule(mid, label, smiles, ev, diagnostics):
    """Candidate connectivity for a legible vinylallene; source leaves its axis unassigned."""
    diag = st.stereo_inspect(smiles)
    diagnostics[mid] = diag
    m = st.new_molecule(mid, DOC, ev)
    m["occurrences"] = [occurrence(label, ev)]
    m["identity_status"] = "source_identified"
    m["names"]["original"] = F(label, label, ev)
    m["names"]["english"] = missing(ev)
    m["names"]["abbreviations"] = missing(ev)
    m["names"]["aliases"] = missing(ev)
    m["cas_numbers"] = missing(ev)
    m["representation_type"] = "smiles"
    basis = ("Candidate connectivity transcribed from the source vinylallene drawing and its printed R assignment; "
             f"RDKit canonicalized it to {diag['canonical_isomeric_smiles']}. The source gives no allene-axis assignment, "
             "and toolkit graph diagnostics do not certify axial stereochemistry.")
    m["smiles"] = F(label, diag["canonical_isomeric_smiles"], ev,
                    method="structure_derivation", derivation=basis)
    m["generic_template"] = None
    m["stereochemistry"]["atom_centers"] = derived(diag["atom_centers"], ev, basis)
    m["stereochemistry"]["double_bonds"] = derived(diag["double_bonds"], ev, basis)
    m["stereochemistry"]["axial"] = derived({"type": "allene axis", "assignment": "unspecified"}, ev, basis)
    m["stereochemistry"]["mixture"] = missing(ev)
    m["stereochemistry"]["assignment"] = F("no axial descriptor shown", "allene-axis configuration unassigned", ev)
    m["stereochemistry"]["chirality_status"] = derived("unspecified", ev, basis)
    m["stereochemistry"]["geometry_status"] = derived("unspecified", ev, basis)
    m["stereochemistry"]["representation_status"] = derived("partial", ev, basis)
    m["stereochemistry"]["isomeric_smiles"] = F(label, diag["canonical_isomeric_smiles"], ev,
                                                  method="structure_derivation", derivation=basis)
    m["stereochemistry"]["identity_key"] = derived(diag["stereo_graph_key"], ev, basis)
    m["identity_note"] = F("source drawing plus printed R assignment",
                           "Candidate connectivity only; allene axial configuration is not encoded or assigned.", ev)
    m["external_properties"] = []
    return m


def generic_molecule(mid, label, ev, kind):
    m = st.new_molecule(mid, DOC, ev)
    m["occurrences"] = [occurrence(label, ev)]
    m["identity_status"] = "source_identified"
    m["names"]["original"] = F(label, label, ev)
    m["names"]["english"] = missing(ev)
    m["names"]["abbreviations"] = missing(ev)
    m["names"]["aliases"] = missing(ev)
    m["cas_numbers"] = missing(ev)
    m["representation_type"] = "generic_template"
    m["smiles"] = missing(ev)
    if kind == "substrate":
        fixed = {"description": "vinylallene framework shown in Table 2", "features": ["C=C=C", "terminal alkene"]}
        variables = ["R1", "R2", "R3"]
    else:
        fixed = {"description": "seven-membered tropone [4+2+1] cycloadduct framework shown in Table 2"}
        variables = ["R1", "R2", "R3"]
    m["generic_template"] = {
        "fixed_structure": F(label, fixed, ev),
        "variables": [{
            "label": v,
            "attachment_points": F(f"{v} position in the source drawing", f"{v} position in the source drawing", ev),
            "definition": F(f"Variable substituent {v}", f"Variable substituent {v}", ev),
            "allowed_values": missing(ev),
            "constraints": missing(ev),
        } for v in variables],
        "instance_mapping": [],
        "reason": F("R-group structure in Table 2 header", "R-group source drawing is not a single concrete molecule.", ev),
    }
    for key, val in (("chirality_status", "unresolved"), ("geometry_status", "unresolved"),
                     ("representation_status", "not_encoded")):
        m["stereochemistry"][key] = derived(val, ev, "R-group template cannot establish a unique stereochemical graph.")
    for key in ("atom_centers", "double_bonds", "axial", "mixture", "assignment",
                "isomeric_smiles", "identity_key"):
        m["stereochemistry"][key] = missing(ev)
    m["identity_note"] = F("R-group structure", "Generic source template; no wildcard concrete SMILES asserted.", ev)
    m["external_properties"] = []
    return m


def participant(pid, mid, role, ev, stage="reaction", charges=None, ref=None,
                product_role=None):
    p = st.new_participant(pid, mid, role, ev, stage)
    p["product_role"] = (F(product_role, product_role, ev) if product_role else missing(ev))
    charges = charges or {}
    for key in p["charges"]:
        value = charges.get(key)
        if value is None:
            p["charges"][key] = missing(ev)
        else:
            raw = value.pop("raw", None) if isinstance(value, dict) and "raw" in value else None
            p["charges"][key] = F(raw or value, value, ev)
    p["equivalence_reference"] = F(ref, ref, ev) if ref else missing(ev)
    return p


def step(step_id, pids, ev, temperature=None, time=None, atmosphere=None,
         pressure=None, operation=None):
    s = st.new_step(step_id, ev)
    s["participant_ids"] = list(pids)
    s["operation"] = F(operation, operation, ev) if operation else missing(ev)
    s["temperature"] = F(temperature, temperature, ev) if temperature is not None else missing(ev)
    s["time"] = F(time, time, ev) if time is not None else missing(ev)
    s["atmosphere"] = F(atmosphere, atmosphere, ev) if atmosphere else missing(ev)
    s["pressure"] = F(pressure, pressure, ev) if pressure is not None else missing(ev)
    return s


def outcome(oid, kind, targets, ev, value=None, raw=None, yield_type=None,
            basis=None, role="main product", scope=None):
    o = st.new_outcome(oid, kind, targets, ev)
    o["product_role"] = F(role, role, ev) if role else missing(ev)
    o["value"] = F(raw if raw is not None else value, value, ev) if value is not None else missing(ev)
    o["yield_type"] = F(yield_type, yield_type, ev) if yield_type else missing(ev)
    o["measurement_basis"] = F(basis, basis, ev) if basis else missing(ev)
    o["scope"] = F(scope, scope, ev) if scope else missing(ev)
    return o


def classification_unresolved(ev, reason="Exact atom mapping and net bond changes are unresolved.", not_applicable=False):
    coverage = "not_applicable" if not_applicable else "unresolved"
    return {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": missing(ev),
        "bond_changes": [],
        "labels": derived([], ev, reason),
        "coverage": derived(coverage, ev, reason),
        "limitations": [reason],
        "model_votes": [],
    }


def classification_cc(ev, formed=3, promoted=False, label_reason=None):
    reason = label_reason or (
        "Class-level changes inferred from the source-depicted ring-forming transformation; exact atom identities and full mapping are unavailable."
    )
    changes = []
    for i in range(formed):
        ch = {"change_id": f"b{i+1}", "atoms": [{"element": "C"}, {"element": "C"}],
              "before_order": 0, "after_order": 1, "change_type": "formed"}
        changes.append(F(None, ch, ev, method="structure_derivation", derivation=reason))
    labels = [{"category_number": 2, "code": "C-C",
               "supporting_change_ids": [f"b{i+1}" for i in range(formed)]}]
    if promoted:
        cid = f"b{formed+1}"
        ch = {"change_id": cid, "atoms": [{"element": "C"}, {"element": "C"}],
              "before_order": 1, "after_order": 2, "change_type": "order_increased"}
        changes.append(F(None, ch, ev, method="structure_derivation", derivation=reason))
        labels.append({"category_number": 13, "code": "DB-C-C", "supporting_change_ids": [cid]})
    return {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": missing(ev),
        "bond_changes": changes,
        "labels": F(None, labels, ev, method="structure_derivation", derivation=reason),
        "coverage": derived("partial", ev, reason),
        "limitations": [reason, "CO C–O bond-order change, allene/alkene bond-order changes, implicit hydrogens, coproducts, and exact atom correspondence are not fully mapped."],
        "model_votes": [],
    }


def base_reaction(rid, ev, record_type="specific_experiment", categories=(), location="",
                  reaction_name=None, transformation=None):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = F(list(categories), list(categories), ev) if categories else missing(ev)
    r["location"] = F(location, {"description": location}, ev)
    r["reaction_name"] = F(reaction_name, reaction_name, ev) if reaction_name else missing(ev)
    r["transformation_type"] = derived(transformation, ev, "Normalized from the source reaction name/depiction.") if transformation else missing(ev)
    r["general_scheme_ids"] = []
    r["participants"] = []
    r["reaction_scale"] = missing(ev)
    for key in r["solvent_conditions"]:
        r["solvent_conditions"][key] = missing(ev)
    r["steps"] = []
    for key in r["special_conditions"]:
        r["special_conditions"][key] = missing(ev)
    for key in r["procedure"]:
        r["procedure"][key] = missing(ev)
    r["outcomes"] = []
    r["classification"] = classification_unresolved(ev)
    r["notes"] = []
    r["review"] = {
        "status": "needs_review", "method": "model", "reviewer": None,
        "checked_fields": [],
        "unresolved": ["Supporting Information was not supplied; full charges, addition order, workup, purification, and characterization are unavailable."],
    }
    r["completeness"] = "partial"
    return r


def set_solvent(r, solvent_mid, raw_name, concentration=None, basis_mid=None, ev=()):
    r["solvent_conditions"]["composition"] = F(raw_name, [solvent_mid], ev)
    r["solvent_conditions"]["mixture_ratio"] = missing(ev)
    r["solvent_conditions"]["total_volume"] = missing(ev)
    r["solvent_conditions"]["concentration"] = F(concentration, concentration, ev) if concentration else missing(ev)
    r["solvent_conditions"]["concentration_basis"] = F(basis_mid, basis_mid, ev) if basis_mid else missing(ev)


def add_common_participants(r, substrate, product, catalyst, solvent, ev, scale=None,
                            catalyst_loading=5, include_product=True):
    pids = []
    sub_charge = {"amount": scale} if scale else {}
    r["participants"].append(participant("substrate", substrate, "reactant", ev, charges=sub_charge))
    pids.append("substrate")
    r["participants"].append(participant("co", "mol_CO", "reagent", ev))
    pids.append("co")
    r["participants"].append(participant("catalyst", catalyst, "catalyst", ev,
                                         charges={"mol_percent": {"number": catalyst_loading, "unit": "mol%"}},
                                         ref=substrate))
    pids.append("catalyst")
    r["participants"].append(participant("solvent", solvent, "solvent", ev))
    pids.append("solvent")
    if include_product and product:
        r["participants"].append(participant("product", product, "product", ev, product_role="main product"))
        pids.append("product")
    return pids


def make_table1(entry, catalyst, loading, pressure, solvent_mid, solvent_name, conc,
                temp, yld, product=True, ev=()):
    rid = f"table1_entry{entry}_nmr"
    r = base_reaction(rid, ev, categories=["optimization"],
                      location=f"Table 1, entry {entry}",
                      reaction_name="[4+2+1] reaction",
                      transformation="rhodium-catalyzed three-component [4+2+1] cycloaddition")
    scale = {"number": 0.1, "unit": "mmol"}
    pids = add_common_participants(r, "mol_1d", "mol_2d", catalyst, solvent_mid, ev,
                                   scale=scale, catalyst_loading=loading,
                                   include_product=product)
    if yld != "c.m.":
        # Explicitly keep the NMR internal standard at the characterization stage,
        # not as a reaction-time reagent. Its amount is not reported.
        r["participants"].append(participant("nmr_internal_standard", "mol_nmr_internal_standard",
                                             "reagent", ev, stage="characterization"))
    r["reaction_scale"] = F("All reactions were performed on 0.1 mmol scale.",
                            {"reference_molecule_id": "mol_1d", "amount": scale}, ev)
    set_solvent(r, solvent_mid, solvent_name, {"number": conc, "unit": "M"}, "mol_1d", ev)
    temperature = temp if isinstance(temp, str) else {"number": temp, "unit": "°C"}
    # The bundled quantity vocabulary has no pressure unit; preserve the explicit
    # source expression rather than inventing or changing the validator.
    pressure_obj = f"{pressure} atm"
    r["steps"] = [step("reaction", pids, ev, temperature=temperature,
                       time={"number": 12, "unit": "h"}, atmosphere="CO",
                       pressure=pressure_obj,
                       operation="React under the Table 1 conditions for 12 h.")]
    r["procedure"]["original_text"] = F(
        f"Table 1 entry {entry}: catalyst {loading} mol%, P(CO) {pressure} atm, {solvent_name}, {conc} M, {temp}, yield {yld}.",
        f"Table 1 entry {entry}", ev)
    if yld == "c.m.":
        r["outcomes"] = [outcome("complex_mixture", "other", [], ev, value="complex mixture", raw="c.m.", role=None,
                                 scope="observed reaction mixture")]
        r["classification"] = classification_unresolved(ev, "Complex mixture has no identified net product or bond map.")
        r["review"]["unresolved"].append("The complex-mixture composition, conversion, and recovery were not reported.")
    else:
        r["outcomes"] = [outcome("nmr_yield", "yield", ["mol_2d"], ev,
                                 value={"number": yld, "unit": "%"}, raw=f"{yld}",
                                 yield_type="1H NMR yield",
                                 basis="1H-NMR yield using 1,3,5-trimethoxybenzene as an internal standard")]
        r["classification"] = classification_cc(ev)
    r["notes"].append(F("All Table 1 runs use substrate 1d; yields outside parentheses are 1H NMR yields.",
                        "All Table 1 runs use substrate 1d; yields outside parentheses are 1H NMR yields.", ev))
    r["review"]["unresolved"].append("The stoichiometric charge of CO and the detailed operational procedure are not reported in the main text.")
    return r


def make_scope(entry, sub, prod, yld, temp=0, loading=5, warm=False, merged_opt=False, ev=()):
    rid = f"table2_entry{entry}"
    categories = ["substrate_scope"] + (["optimization"] if merged_opt else [])
    r = base_reaction(rid, ev, categories=categories,
                      location=f"Table 2, entry {entry}" + ("; also Table 1 entry 10 isolated 0.3 mmol result" if merged_opt else ""),
                      reaction_name="three-component [4+2+1] reaction",
                      transformation="rhodium-catalyzed three-component [4+2+1] cycloaddition")
    scale = {"number": 0.3, "unit": "mmol"}
    pids = add_common_participants(r, sub, prod, "mol_rh_co_cl_dimer", "mol_DCM", ev,
                                   scale=scale, catalyst_loading=loading)
    r["reaction_scale"] = F("The reactions were carried out on 0.3 mmol scale.",
                            {"reference_molecule_id": sub, "amount": scale}, ev)
    set_solvent(r, "mol_DCM", "DCM", {"number": 0.1, "unit": "M"}, sub, ev)
    if warm:
        r["steps"] = [
            step("hold_0C", pids, ev, temperature={"number": 0, "unit": "°C"},
                 time={"number": 12, "unit": "h"}, atmosphere="CO",
                 pressure="1 atm",
                 operation="Stir at 0 °C; starting material remained by TLC."),
            step("warm_and_hold_RT", pids, ev, temperature="RT",
                 time={"number": 24, "unit": "h"}, atmosphere="CO",
                 pressure="1 atm",
                 operation="Warm to room temperature and stir for 24 h."),
        ]
        r["review"]["unresolved"].append("The warming/ramp duration from 0 °C to RT was not reported.")
    else:
        temperature = "RT" if temp == "RT" else {"number": temp, "unit": "°C"}
        r["steps"] = [step("reaction", pids, ev, temperature=temperature,
                           time={"number": 12, "unit": "h"}, atmosphere="CO",
                           pressure="1 atm",
                           operation="React under the Table 2 conditions.")]
    r["procedure"]["original_text"] = F(
        f"Table 2 entry {entry}: 0.3 mmol scale, [Rh(CO)2Cl]2 {loading} mol%, CO 1 atm, DCM 0.1 M, product yield {yld}%.",
        f"Table 2 entry {entry}", ev)
    r["outcomes"] = [outcome("isolated_yield", "yield", [prod], ev,
                             value={"number": yld, "unit": "%"}, raw=f"{yld}",
                             yield_type="isolated yield",
                             basis="isolated yield, average of two runs")]
    r["classification"] = classification_cc(ev)
    r["notes"].append(F("Substrates and products were not bench stable when stored for a long time.",
                        "Substrates and products were not bench stable when stored for a long time.", ev))
    if merged_opt:
        r["notes"].append(F(
            "The 70% isolated 0.3 mmol result in Table 1 entry 10 is treated as the same reported result as Table 2 entry 1; the separate 0.1 mmol/81% NMR run remains its own record.",
            "Cross-table duplicate merged; the 0.1 mmol NMR run is separate.", ev))
    return r


def main():
    dataset = json.loads((RUN / "dataset.json").read_text())
    doc = next(d for d in dataset["documents"] if d["document_id"] == DOC)
    doc["title"] = "Rhodium(I)-Catalyzed Three-Component [4+2+1] Cycloaddition of Two Vinylallenes and CO"
    (RUN / "dataset.json").write_text(json.dumps(dataset, ensure_ascii=False, indent=2) + "\n")
    bundle = st.empty(dataset["documents"])
    diagnostics = {}

    # Source-labelled molecules. Candidate connectivity is transcribed only where
    # the compact drawing and printed R assignment are unambiguous enough.
    molecules = []
    vinylallene_candidates = {
        "1d": "c1ccccc1CCC=C=CC=C",
        "1e": "CCCCCCC=C=CC=C",
        "1f": "c1ccccc1CC=C=CC=C",
        "1g": "C1CCCCC1CC=C=CC=C",
        "1h": "CC(=O)OCC=C=CC=C",
        "1i": "CC(C)(C)[Si](C)(C)OCC=C=CC=C",
        "1j": "C1CCCCC1C=C=CC=C",
        "1k": "c1ccccc1C=C=CC=C",
        "1m": "c1ccccc1CCC=C=CC(=C)C",
        "1n": "CCCCCC=C=CC(=C)C",
    }
    for label in [f"1{x}" for x in "abcdefghijklmnopqr"]:
        page = 1 if label in {"1a", "1b", "1c"} else 2 if label == "1d" else 3
        extra = []
        if label in {"1a", "1b", "1c", "1d"}:
            extra = [(label, E[3])]
        if label in vinylallene_candidates:
            m = vinylallene_smiles_molecule(f"mol_{label}", label, vinylallene_candidates[label], E[page], diagnostics)
            m["occurrences"] += [occurrence(l, e) for l, e in extra]
            molecules.append(m)
        else:
            molecules.append(unresolved_molecule(f"mol_{label}", label, E[page], "vinylallene", extra))
    for label in [f"2{x}" for x in "abcdefghijklmno"]:
        page = 1 if label in {"2a", "2b", "2c"} else 2 if label == "2d" else 3
        extra = [(label, E[3])] if label == "2d" else []
        molecules.append(unresolved_molecule(f"mol_{label}", label, E[page], "tropone", extra))
    molecules += [
        smiles_molecule("mol_3", "3", "4-phenylbut-1-yne", "c1ccccc1CCC#C", E[3], diagnostics),
        unresolved_molecule("mol_4", "4", E[3], "cycloadduct"),
        unresolved_molecule("mol_5", "5", E[3], "cycloadduct"),
        unresolved_molecule("mol_rh_co_cl_dimer", "[Rh(CO)2Cl]2", E[2], "unresolved", [("[Rh(CO)2Cl]2", E[3])]),
        unresolved_molecule("mol_rh_cod2_sbf6", "Rh(COD)2SbF6", E[2], "unresolved"),
        generic_molecule("mol_table2_vinylallene_generic", "Table 2 generic vinylallene", E[3], "substrate"),
        generic_molecule("mol_table2_tropone_generic", "Table 2 generic tropone product", E[3], "product"),
    ]
    molecules += [
        smiles_molecule("mol_CO", "CO", "carbon monoxide", "[C-]#[O+]", E[1] + E[2] + E[3], diagnostics),
        smiles_molecule("mol_DCE", "DCE", "1,2-dichloroethane", "ClCCCl", E[1] + E[2] + E[3], diagnostics),
        smiles_molecule("mol_DCM", "DCM", "dichloromethane", "ClCCl", E[2] + E[3], diagnostics),
        smiles_molecule("mol_toluene", "toluene", "toluene", "Cc1ccccc1", E[2], diagnostics),
        smiles_molecule("mol_MeCN", "MeCN", "acetonitrile", "CC#N", E[2], diagnostics),
        smiles_molecule("mol_HFIP", "HFIP", "hexafluoroisopropanol", "OC(C(F)(F)F)C(F)(F)F", E[2], diagnostics),
        smiles_molecule("mol_dioxane", "dioxane", "1,4-dioxane", "O1CCOCC1", E[2], diagnostics),
        smiles_molecule("mol_nmr_internal_standard", "1,3,5-trimethoxybenzene", "1,3,5-trimethoxybenzene", "COc1cc(OC)cc(OC)c1", E[2], diagnostics),
    ]
    bundle["molecules"] = molecules

    reactions = []
    # Scheme 1a: seven cited generic arrows, retained as literature examples.
    prior_types = ["[4+2] (diene or alkyne)", "[4+2]", "[4+1]", "[4+4]",
                   "[4+4+1]", "[3+2]/[2+3]", "[4+2] (dieneophile)"]
    for i, typ in enumerate(prior_types, 1):
        r = base_reaction(f"scheme1a_prior_arrow{i}", E[1], record_type="literature_example",
                          location=f"Scheme 1a, prior-work arrow {i}", reaction_name=typ,
                          transformation="cited vinylallene cycloaddition")
        r["classification"] = classification_unresolved(E[1], "Generic cited prior-work drawing lacks a complete atom-mapped concrete reaction.")
        r["notes"] = [F("Cited prior work summarized in Scheme 1a; not an experiment executed in this paper.",
                        "Cited prior work; not a current-paper experiment.", E[1])]
        r["review"]["unresolved"] = ["Generic substituents and individual cited example identities are not resolved from the overview figure."]
        reactions.append(r)
    # Scheme 1b has two conceptual arrows for the same cited prior intramolecular sequence.
    for i, name in enumerate(["direct intramolecular [4+2+1] cycloaddition", "1,3-acyloxy migration to vinylallene intermediate"], 1):
        r = base_reaction(f"scheme1b_prior_arrow{i}", E[1], record_type="literature_example",
                          location=f"Scheme 1b, prior-work arrow {i}", reaction_name=name,
                          transformation="cited intramolecular vinylallene sequence")
        r["classification"] = classification_unresolved(E[1], "Cited generic scheme is not a concrete experiment and is not atom mapped.")
        r["notes"] = [F("Cited prior work [1q], not an experiment executed in this paper.",
                        "Cited prior work; not a current-paper experiment.", E[1])]
        r["review"]["unresolved"] = ["Generic R groups and experimental details belong to cited work and are not supplied here."]
        reactions.append(r)

    # Scheme 1c preliminary current-paper reactions 1a–1c.
    for lab, yld in (("a", 50), ("b", 21), ("c", 42)):
        ev = E[1] + E[2]
        r = base_reaction(f"scheme1c_{lab}", ev, categories=["substrate_scope"],
                          location=f"Scheme 1c, substrate 1{lab}", reaction_name="three-component [4+2+1] cycloaddition",
                          transformation="rhodium-catalyzed three-component [4+2+1] cycloaddition")
        pids = add_common_participants(r, f"mol_1{lab}", f"mol_2{lab}", "mol_rh_co_cl_dimer", "mol_DCE", ev)
        set_solvent(r, "mol_DCE", "DCE", None, f"mol_1{lab}", ev)
        r["steps"] = [step("reaction", pids, ev, temperature={"number": 40, "unit": "°C"},
                           time={"number": 12, "unit": "h"}, atmosphere="CO",
                           pressure="1 atm", operation="React under Scheme 1c conditions.")]
        r["participants"][2]["charges"]["mol_percent"] = F("5 mol%", {"number": 5, "unit": "mol%"}, ev)
        r["procedure"]["original_text"] = F("1 atm CO; 5 mol% [Rh(CO)2Cl]2; DCE; 40 °C; 12 h", "Scheme 1c conditions", ev)
        r["outcomes"] = [outcome("reported_yield", "yield", [f"mol_2{lab}"], ev,
                                 value={"number": yld, "unit": "%"}, raw=f"{yld}%",
                                 basis="isolated product (prose states products 2a–2c were isolated)")]
        r["classification"] = classification_cc(ev)
        r["review"]["unresolved"].append("Reaction scale, substrate charge, concentration, workup, purification, and exact yield protocol are not given in the main text.")
        reactions.append(r)

    # Proposed mechanism in Scheme 2.
    mech = base_reaction("scheme2_proposed_mechanism", E[2], record_type="proposed_mechanism",
                         location="Scheme 2", reaction_name="proposed mechanism of [4+2+1] reaction of two vinylallenes and CO",
                         transformation="author-proposed oxidative cyclization/allene insertion/CO insertion/reductive elimination sequence")
    mech["participants"] = [participant("vinylallene", "mol_1d", "reactant", E[2]), participant("co", "mol_CO", "reagent", E[2])]
    mech["procedure"]["original_text"] = F(
        "Oxidative cyclization to I, allene insertion to II, CO insertion to III, then reductive elimination.",
        "Author-proposed sequence; not an executed experimental procedure.", E[2])
    mech["classification"] = classification_unresolved(E[2], "Proposed elementary-step sequence is not a single experimentally mapped net reaction.")
    mech["review"]["unresolved"] = ["Intermediates I–III are proposed, not experimentally established in the supplied pages; exact atom mapping is absent."]
    reactions.append(mech)

    table1_rows = [
        (1, "mol_rh_co_cl_dimer", 5, 1, "mol_DCE", "DCE", 0.1, 60, 49, True),
        (2, "mol_rh_co_cl_dimer", 5, 1, "mol_DCE", "DCE", 0.02, 60, 44, True),
        (3, "mol_rh_co_cl_dimer", 5, 1, "mol_toluene", "toluene", 0.1, 60, 12, True),
        (4, "mol_rh_co_cl_dimer", 5, 1, "mol_MeCN", "MeCN", 0.1, 60, 60, True),
        (5, "mol_rh_co_cl_dimer", 5, 1, "mol_HFIP", "HFIP", 0.1, 60, 27, True),
        (6, "mol_rh_co_cl_dimer", 5, 1, "mol_dioxane", "dioxane", 0.1, 60, 76, True),
        (7, "mol_rh_co_cl_dimer", 5, 1, "mol_DCM", "DCM", 0.1, "RT", 78, True),
        (8, "mol_rh_co_cl_dimer", 5, 0.2, "mol_DCM", "DCM", 0.1, "RT", 62, True),
        (9, "mol_rh_co_cl_dimer", 2, 1, "mol_DCM", "DCM", 0.1, "RT", 76, True),
        (10, "mol_rh_co_cl_dimer", 5, 1, "mol_DCM", "DCM", 0.1, 0, 81, True),
        (11, "mol_rh_cod2_sbf6", 10, 1, "mol_DCM", "DCM", 0.1, "RT", "c.m.", False),
    ]
    for row in table1_rows:
        reactions.append(make_table1(*row, ev=E[2]))

    # Parenthetical 66% is a second, 0.3 mmol isolated experiment in Table 1 entry 9.
    r = base_reaction("table1_entry9_isolated_scale", E[2], categories=["optimization", "scale_up"],
                      location="Table 1, entry 9, parenthetical isolated 0.3 mmol result",
                      reaction_name="[4+2+1] reaction",
                      transformation="rhodium-catalyzed three-component [4+2+1] cycloaddition")
    scale = {"number": 0.3, "unit": "mmol"}
    pids = add_common_participants(r, "mol_1d", "mol_2d", "mol_rh_co_cl_dimer", "mol_DCM", E[2], scale=scale, catalyst_loading=2)
    r["reaction_scale"] = F("isolated product based on a 0.3 mmol scale reaction", {"reference_molecule_id": "mol_1d", "amount": scale}, E[2])
    set_solvent(r, "mol_DCM", "DCM", {"number": 0.1, "unit": "M"}, "mol_1d", E[2])
    r["steps"] = [step("reaction", pids, E[2], temperature="RT", time={"number": 12, "unit": "h"}, atmosphere="CO", pressure="1 atm")]
    r["outcomes"] = [outcome("isolated_yield", "yield", ["mol_2d"], E[2], value={"number": 66, "unit": "%"}, raw="(66)", yield_type="isolated yield", basis="isolated product based on a 0.3 mmol scale reaction")]
    r["classification"] = classification_cc(E[2])
    r["review"]["unresolved"].append("The main text does not provide the workup or purification for the isolated 0.3 mmol run.")
    reactions.append(r)

    # Table 2 generic R-group arrow.
    gen = base_reaction("table2_generic_scheme", E[3], record_type="generic_scheme", location="Table 2 header",
                        reaction_name="scope of the three-component [4+2+1] reaction",
                        transformation="generic rhodium-catalyzed [4+2+1] cycloaddition")
    gen["participants"] = [participant("generic_substrate", "mol_table2_vinylallene_generic", "reactant", E[3]),
                           participant("co", "mol_CO", "reagent", E[3]),
                           participant("generic_product", "mol_table2_tropone_generic", "product", E[3], product_role="main product")]
    gen["classification"] = classification_cc(E[3])
    gen["review"]["unresolved"] = ["R-group generic scheme is not one concrete experiment; instance structures are retained separately without unsafe SMILES."]
    reactions.append(gen)

    scope_rows = [
        (1, "mol_1d", "mol_2d", 70, 0, 5, False, True),
        (2, "mol_1e", "mol_2e", 74, 0, 5, False, False),
        (3, "mol_1f", "mol_2f", 52, 0, 5, False, False),
        (4, "mol_1g", "mol_2g", 58, 0, 5, False, False),
        (5, "mol_1h", "mol_2h", 55, 0, 5, False, False),
        (6, "mol_1i", "mol_2i", 51, 0, 5, False, False),
        (7, "mol_1j", "mol_2j", 51, 0, 5, False, False),
        (8, "mol_1k", "mol_2k", 72, "RT", 10, False, False),
        (9, "mol_1l", "mol_2l", 56, 0, 5, True, False),
        (10, "mol_1m", "mol_2m", 72, 0, 5, False, False),
        (11, "mol_1n", "mol_2n", 68, 0, 5, False, False),
        (12, "mol_1o", "mol_2o", 33, 0, 5, True, False),
    ]
    for row in scope_rows:
        r = make_scope(*row, ev=E[3] + (E[2] if row[0] == 1 else []))
        r["general_scheme_ids"] = [f"{DOC}:table2_generic_scheme"]
        reactions.append(r)

    # Text reports 1a–1c under optimized DCM conditions at 40 °C, exact yields only as "similar".
    for lab in "abc":
        r = base_reaction(f"optimized_40C_1{lab}", E[3], categories=["substrate_scope"],
                          location=f"Page 3 prose: 1{lab} under optimized conditions at 40 °C",
                          reaction_name="three-component [4+2+1] reaction",
                          transformation="rhodium-catalyzed three-component [4+2+1] cycloaddition")
        pids = add_common_participants(r, f"mol_1{lab}", f"mol_2{lab}", "mol_rh_co_cl_dimer", "mol_DCM", E[3])
        set_solvent(r, "mol_DCM", "DCM", {"number": 0.1, "unit": "M"}, f"mol_1{lab}", E[3])
        r["steps"] = [step("reaction", pids, E[3], temperature={"number": 40, "unit": "°C"},
                           time={"number": 12, "unit": "h"}, atmosphere="CO", pressure="1 atm",
                           operation="Use optimized conditions but raise temperature to 40 °C.")]
        r["outcomes"] = [outcome("qualitative_yield", "other", [f"mol_2{lab}"], E[3],
                                 value="similar yield to Scheme 1c", raw="gave similar yields to those obtained in Scheme 1c",
                                 role="main product", scope=f"product 2{lab}; exact percentage not reported")]
        r["classification"] = classification_cc(E[3])
        r["review"]["unresolved"].append("Exact yield, scale, and measurement basis for the 40 °C rerun are not reported; Scheme 1c percentages were not copied onto this distinct condition.")
        reactions.append(r)

    # Scheme 3 Pauson–Khand product.
    r = base_reaction("scheme3_product4", E[3], categories=["other"], location="Scheme 3",
                      reaction_name="[2+2+1] Pauson–Khand reaction",
                      transformation="rhodium-catalyzed three-component [2+2+1] cycloaddition")
    r["participants"] = [
        participant("vinylallene", "mol_1h", "reactant", E[3], charges={"equivalents": {"number": 1, "unit": "equiv"}}, ref="mol_1h"),
        participant("alkyne", "mol_3", "reactant", E[3], charges={"equivalents": {"number": 2, "unit": "equiv"}}, ref="mol_1h"),
        participant("co", "mol_CO", "reagent", E[3]),
        participant("catalyst", "mol_rh_co_cl_dimer", "catalyst", E[3], charges={"mol_percent": {"number": 5, "unit": "mol%"}}, ref="mol_1h"),
        participant("solvent", "mol_DCE", "solvent", E[3]),
        participant("product", "mol_4", "product", E[3], product_role="Pauson–Khand product"),
    ]
    set_solvent(r, "mol_DCE", "DCE", None, "mol_1h", E[3])
    r["steps"] = [step("syringe_pump_reaction", [p["participant_id"] for p in r["participants"]], E[3],
                       temperature={"number": 80, "unit": "°C"}, time={"number": 15, "unit": "h"},
                       atmosphere="CO", pressure="1 atm",
                       operation="Scheme labels 1h delivery with a syringe pump; exact addition rate/duration is not reported.")]
    r["procedure"]["original_text"] = F("1h (1 eq., syringe pump) + 3 (2 eq.) + CO; [Rh(CO)2Cl]2 5 mol%; DCE, 80 °C, 15 h", "Scheme 3 conditions", E[3])
    r["outcomes"] = [outcome("yield", "yield", ["mol_4"], E[3], value={"number": 23, "unit": "%"}, raw="23%", basis=None)]
    r["classification"] = classification_cc(E[3])
    r["review"]["unresolved"] += ["Yield measurement basis is not reported.", "Syringe-pump addition rate, addition duration, scale, concentration, workup, and purification are not reported."]
    reactions.append(r)

    # Scheme 4: one successful [4+2] cyclization and two complex-mixture outcomes.
    for lab, result in (("p", 81), ("q", "complex mixture"), ("r", "complex mixture")):
        r = base_reaction(f"scheme4_1{lab}", E[3], categories=["other"], location=f"Scheme 4, substrate 1{lab}",
                          reaction_name="attempted [4+2+1] reaction",
                          transformation="intramolecular [4+2] cycloaddition" if lab == "p" else "unresolved competing reaction")
        pids = add_common_participants(r, f"mol_1{lab}", "mol_5" if lab == "p" else None,
                                       "mol_rh_co_cl_dimer", "mol_DCE", E[3], include_product=(lab == "p"))
        set_solvent(r, "mol_DCE", "DCE", None, f"mol_1{lab}", E[3])
        r["steps"] = [step("reaction", pids, E[3], temperature={"number": 60, "unit": "°C"},
                           time={"number": 24, "unit": "h"}, atmosphere="CO", pressure="1 atm")]
        r["procedure"]["original_text"] = F("[Rh(CO)2Cl]2 (5 mol%), CO (1 atm), DCE, 60 °C, 24 h", "Scheme 4 conditions", E[3])
        if lab == "p":
            r["outcomes"] = [outcome("yield", "yield", ["mol_5"], E[3], value={"number": 81, "unit": "%"}, raw="81%")]
            r["classification"] = classification_cc(E[3], formed=2, promoted=True,
                                                      label_reason="Two C–C connections and one promoted C=C are inferred at class level from the source-depicted intramolecular [4+2] cycloadduct; exact atom mapping is unavailable.")
            r["review"]["unresolved"].append("Yield measurement basis, scale, concentration, workup, and purification are not reported.")
        else:
            r["outcomes"] = [outcome("complex_mixture", "other", [], E[3], value="complex mixture", raw="Complex mixture", role=None, scope="unidentified reaction mixture")]
            r["classification"] = classification_unresolved(E[3], "Complex mixture has no identified product or classifiable net transformation.")
            r["review"]["unresolved"].append("No product identity, yield, conversion, recovery, or mixture composition is reported.")
        reactions.append(r)

    bundle["reactions"] = reactions
    bundle["source_context"] = [
        {
            "context_id": f"{DOC}:bibliographic_metadata",
            "kind": "bibliographic_metadata",
            "title": doc["title"],
            "doi": "10.1002/chem.202005443",
            "evidence_ids": E[1] + E[4],
        },
        {
            "context_id": f"{DOC}:replacement_component_failures",
            "kind": "non_specific_experimental_observation",
            "text": "The prose states that allene, alkyne, and isoprene could not be incorporated into the [4+2+1] cycloadducts, but identifies no complete experiments on the supplied pages and sends details to unavailable SI; no concrete records were fabricated.",
            "evidence_ids": E[3],
        },
        {
            "context_id": f"{DOC}:regioselectivity_and_xray",
            "kind": "discussion_and_characterization",
            "text": "The authors report high regioselectivity with no other isomeric carbocyclic ketones detected and say the structure of 2l was confirmed by X-ray analysis of a derivative. The derivative identity/data are not supplied here.",
            "evidence_ids": E[3],
        },
        {
            "context_id": f"{DOC}:spectral_inventory",
            "kind": "coverage_note",
            "text": "Pages 1–4 contain no standalone mass-spectrum peak table, molecular-weight table, or fragment-relative-abundance table. No characterization-only synthesis record was created.",
            "evidence_ids": E[1] + E[2] + E[3] + E[4],
        },
        {
            "context_id": f"{DOC}:page4_summary",
            "kind": "summary_and_references",
            "text": "Page 4 contains the article summary, acknowledgements, conflict statement, keywords, references, and publication dates; no current-paper experimental arrow or table row occurs.",
            "evidence_ids": E[4],
        },
    ]

    st.write(RUN / "submitted.json", bundle)
    st.write(RUN / "raw_model/stereo_checks.json", diagnostics)
    inventory = {
        "document_id": DOC,
        "title": doc["title"],
        "source_is_presentation_conversion": False,
        "si_association": "unestablished",
        "pages": [
            {"page": 1, "visual_evidence_id": E[1][0], "findings": ["Scheme 1a: seven cited prior-work arrows", "Scheme 1b: two cited prior-work arrows", "Scheme 1c: three current experiments 1a–1c"]},
            {"page": 2, "visual_evidence_id": E[2][0], "findings": ["Scheme 2 proposed mechanism", "Table 1: eleven optimization rows", "Table 1 entries 9 and 10 additionally report separate 0.3 mmol isolated results"]},
            {"page": 3, "visual_evidence_id": E[3][0], "findings": ["Table 2: twelve scope rows", "Three prose-reported 1a–1c reruns at 40 °C", "Scheme 3: one Pauson–Khand experiment", "Scheme 4: three attempted experiments"]},
            {"page": 4, "visual_evidence_id": E[4][0], "findings": ["Summary and references only; no new experimental arrow/table row"]},
        ],
        "decisions": [
            "Table 1 entry 10 parenthetical 70% at 0.3 mmol and Table 2 entry 1 are merged as one reported result; the Table 1 0.1 mmol/81% NMR run remains distinct.",
            "Table 1 entry 9 parenthetical 66% at 0.3 mmol is a separate concrete experiment from its 0.1 mmol/76% NMR row.",
            "Scheme 1 cited work is retained as literature_example, Scheme 2 as proposed_mechanism, and Table 2 header as generic_scheme.",
            "No SI was inferred or associated.",
        ],
    }
    st.write(RUN / "raw_model/page_extraction.json", inventory)

    # Coverage record IDs are deliberately explicit and page bounded.
    p1_ids = [f"{DOC}:scheme1a_prior_arrow{i}" for i in range(1, 8)] + [f"{DOC}:scheme1b_prior_arrow{i}" for i in range(1, 3)] + [f"{DOC}:scheme1c_{x}" for x in "abc"]
    p2_ids = [f"{DOC}:scheme2_proposed_mechanism"] + [f"{DOC}:table1_entry{i}_nmr" for i in range(1, 12)] + [f"{DOC}:table1_entry9_isolated_scale", f"{DOC}:table2_entry1"]
    p3_ids = [f"{DOC}:table2_generic_scheme"] + [f"{DOC}:table2_entry{i}" for i in range(1, 13)] + [f"{DOC}:optimized_40C_1{x}" for x in "abc"] + [f"{DOC}:scheme3_product4"] + [f"{DOC}:scheme4_1{x}" for x in "pqr"]
    coverage = {
        "document_id": DOC,
        "semantic_accuracy_verified": False,
        "pages": [
            {"page": 1, "state": "attempted", "regions": [
                {"region": "Scheme 1a, upper-right prior-work overview", "finding": "Seven generic arrows summarize cited vinylallene cycloadditions; classified as literature examples, not this paper's experiments.", "evidence_ids": E[1]},
                {"region": "Scheme 1b, middle-right", "finding": "Two-arrow cited intramolecular [4+2+1]/1,3-acyloxy-migration sequence; classified as prior work.", "evidence_ids": E[1]},
                {"region": "Scheme 1c, lower-right", "finding": "Three current 1a–1c experiments under 1 atm CO, 5 mol% catalyst, DCE, 40 °C, 12 h, giving 2a–2c in 50%, 21%, and 42%.", "evidence_ids": E[1]},
            ], "record_ids": p1_ids, "unresolved": ["Scheme 1a/b generic substituents and cited-work experimental details are not resolved.", "Scheme 1c scale, concentration, detailed operation, and exact yield protocol are absent."]},
            {"page": 2, "state": "attempted", "regions": [
                {"region": "Scheme 2, left column", "finding": "Author-proposed oxidative cyclization, allene insertion, CO insertion and reductive elimination; not an executed experiment.", "evidence_ids": E[2]},
                {"region": "Table 1, upper-right", "finding": "All 11 optimization rows recorded separately; parenthetical isolated yields at 0.3 mmol are separated from the 0.1 mmol NMR experiments.", "evidence_ids": E[2]},
                {"region": "Table 1 footnotes and adjacent prose", "finding": "12 h header, 0.1 mmol scale, NMR/internal-standard basis, and parenthetical isolated 0.3 mmol basis were retained.", "evidence_ids": E[2]},
            ], "record_ids": p2_ids, "unresolved": ["SI procedure, exact CO charge, workup and purification are unavailable.", "Table 1 entry 11 complex-mixture composition/conversion/recovery are unreported."]},
            {"page": 3, "state": "attempted", "regions": [
                {"region": "Table 2 and footnotes, left column", "finding": "Generic R-group arrow plus all 12 substrate-scope entries recorded; entries 9 and 12 retain 0 °C/12 h then RT/24 h as ordered steps; entry 8 retains RT and 10 mol% catalyst.", "evidence_ids": E[3]},
                {"region": "Top-right prose", "finding": "Distinct 1a–1c reruns under optimized conditions at 40 °C recorded with qualitative 'similar yields' only; Scheme 1c percentages were not copied.", "evidence_ids": E[3]},
                {"region": "Scheme 3, middle-right", "finding": "One 1h/3/CO [2+2+1] Pauson–Khand experiment gives product 4 in 23%; syringe-pump rate is unreadable/unreported.", "evidence_ids": E[3]},
                {"region": "Scheme 4, lower-right", "finding": "1p gives [4+2] product 5 in 81%; 1q and 1r give complex mixtures under the same shown conditions.", "evidence_ids": E[3]},
            ], "record_ids": p3_ids, "unresolved": ["Table 2/Schemes 3–4 lack SI operation, workup, and purification.", "Scheme 3 and Scheme 4 yield measurement bases are not stated.", "Exact structures/atom maps and full stereochemical graphs were not safely transcribed."]},
            {"page": 4, "state": "attempted", "regions": [
                {"region": "Left-column summary/acknowledgements", "finding": "Summary of the Rh(I)-catalyzed [4+2+1] reaction; no new experiment or arrow.", "evidence_ids": E[4]},
                {"region": "References and publication history", "finding": "Citations and manuscript dates only; cited works not reclassified as experiments on this page.", "evidence_ids": E[4]},
            ], "record_ids": [], "unresolved": ["No additional current-paper experimental details are present on this page."]},
        ],
    }
    st.write(RUN / "coverage/pages_0001_0004.json", coverage)


if __name__ == "__main__":
    main()
