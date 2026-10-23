#!/usr/bin/env python3
"""Build the evidence-linked v3.1 extraction for DOI 10.1021/acs.orglett.1c02321.

The source is the five-page main article only.  No SI is associated.  Dense drawn
scope structures are deliberately retained as unresolved molecular identities
rather than guessed SMILES.
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
DOC = "doc_a306c90fe6187bb2"
E = {p: f"{DOC}:page{p}:image" for p in range(1, 6)}
T = {p: f"{DOC}:page{p}:text" for p in range(1, 6)}


def f(raw=None, normalized=None, ev=(), method="original_extraction", derivation=None,
      status=None, alternatives=()):
    return st.fact(raw, normalized, ev, method=method, derivation=derivation,
                   status=status, alternatives=alternatives)


def inferred(normalized, ev, derivation):
    return f(None, normalized, ev, method="inference", derivation=derivation)


def occurrence(label, ev):
    return {"document_id": DOC, "compound_label": f(label, label, ev),
            "evidence_ids": list(ev)}


def set_stereo(m, ev, chirality, geometry, representation, note, *, atom=None,
               double=None, axial=None, assignment=None, mixture=None,
               isomeric=None, identity_key=None):
    deriv = note
    values = {
        "atom_centers": atom,
        "double_bonds": double,
        "axial": axial,
        "mixture": mixture,
        "assignment": assignment,
        "chirality_status": chirality,
        "geometry_status": geometry,
        "representation_status": representation,
        "isomeric_smiles": isomeric,
        "identity_key": identity_key,
    }
    for key, val in values.items():
        if val is None:
            m["stereochemistry"][key] = f()
        else:
            m["stereochemistry"][key] = f(None, val, ev, method="structure_derivation",
                                            derivation=deriv)


def unresolved_molecule(mid, label, ev, kind="neutral", note=None):
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"] = [occurrence(label, ev)]
    m["names"]["original"] = f(label, label, ev)
    m["names"]["english"] = f()
    m["names"]["abbreviations"] = f()
    m["names"]["aliases"] = f()
    m["cas_numbers"] = f()
    m["representation_type"] = "unresolved"
    m["smiles"] = f()
    m["generic_template"] = None
    if kind == "product":
        set_stereo(
            m, ev, "partially_specified", "unresolved", "not_encoded",
            "The source drawing depicts relative stereochemical bonds on the bicyclic product, "
            "but no absolute assignment or safely checked complete isomeric SMILES is available.",
            atom={"source": "drawn wedge/dash relative stereochemistry", "absolute_assignment": None},
            assignment={"type": "relative drawing", "absolute_configuration": "unassigned"},
        )
    elif kind == "allene_substrate":
        set_stereo(
            m, ev, "unresolved", "unresolved", "not_encoded",
            "An allene is drawn in the source; axial stereogenicity and any other stereochemistry "
            "were not assigned and are not inferred from absence of labels.",
            axial={"motif": "allene", "assignment": "unreported"},
        )
    else:
        set_stereo(
            m, ev, "unresolved", "unresolved", "unresolved",
            "The source identity/drawing was retained without assuming achirality, racemic state, "
            "E/Z geometry, or a complete graph.",
        )
    m["identity_note"] = inferred(
        note or "Dense or generic source drawing retained by label; complete SMILES was not guessed.",
        ev,
        "The main-text figure identifies the compound label, but the absent SI and dense drawing "
        "prevent a conservative, atom-complete structure transcription in this run.",
    )
    m["external_properties"] = []
    return m


def concrete_molecule(mid, label, name, smiles, ev):
    diag = st.stereo_inspect(smiles)
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"] = [occurrence(label, ev)]
    m["names"]["original"] = f(label, label, ev)
    m["names"]["english"] = f(name, name, ev)
    m["names"]["abbreviations"] = f()
    m["names"]["aliases"] = f()
    m["cas_numbers"] = f()
    m["representation_type"] = "smiles"
    m["smiles"] = f(label, diag["canonical_isomeric_smiles"], ev,
                    method="structure_derivation",
                    derivation="Candidate connectivity derived from the explicit chemical name in the source; RDKit parsing is not source verification.")
    m["generic_template"] = None
    potential = diag["potential_stereo"]
    chirality = "achiral" if not any(x["type"].startswith("Atom_") for x in potential) else "unresolved"
    geometry = "not_applicable" if not any(x["type"].startswith("Bond_") for x in potential) else "unresolved"
    set_stereo(
        m, ev, chirality, geometry, "not_applicable" if chirality == "achiral" and geometry == "not_applicable" else "partial",
        "RDKit graph diagnostics on a name-derived candidate; the source did not report stereochemical assignments.",
        atom=diag["atom_centers"] if diag["atom_centers"] else None,
        double=diag["double_bonds"] if diag["double_bonds"] else None,
        isomeric=diag["canonical_isomeric_smiles"], identity_key=diag["stereo_graph_key"],
    )
    m["identity_note"] = inferred(
        "Name-derived candidate structure; parsing succeeded but the structure remains machine-derived.", ev,
        "The explicit source name/abbreviation was expanded to a candidate molecular graph and checked with standard stereo.",
    )
    m["external_properties"] = []
    return m, diag


def template_molecule(mid, label, ev, variables, fixed):
    m = unresolved_molecule(mid, label, ev, note="Generic R-group scaffold, not a concrete molecule.")
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": f({"figure_label": label}, fixed, ev),
        "variables": [
            {"label": var, "attachment_points": f("source drawing", "source-defined attachment", ev),
             "definition": f(var, "substituent variable", ev),
             "allowed_values": f(), "constraints": f()}
            for var in variables
        ],
        "instance_mapping": [],
        "reason": f("R-group structure in scheme", "The figure defines a family rather than one concrete molecular identity.", ev),
    }
    m["stereochemistry"]["representation_status"] = inferred(
        "not_encoded", ev, "A generic R-group drawing cannot be represented as a complete concrete isomeric SMILES."
    )
    return m


def amount_fact(raw, num, unit, ev, method="original_extraction"):
    return f(raw, {"number": num, "unit": unit}, ev, method=method)


def participant(pid, mid, role, ev, *, amount=None, volume=None, equiv=None,
                molpct=None, ref=None, stage="reaction", product_role=None):
    p = st.new_participant(pid, mid, role, ev, stage=stage)
    if product_role is not None:
        p["product_role"] = f(product_role, product_role, ev)
    if amount is not None:
        p["charges"]["amount"] = amount_fact(amount[0], amount[1], amount[2], ev,
                                                  amount[3] if len(amount) > 3 else "original_extraction")
    if volume is not None:
        p["charges"]["volume"] = amount_fact(volume[0], volume[1], volume[2], ev,
                                                  volume[3] if len(volume) > 3 else "original_extraction")
    if equiv is not None:
        p["charges"]["equivalents"] = amount_fact(equiv[0], equiv[1], "equiv", ev,
                                                       equiv[2] if len(equiv) > 2 else "original_extraction")
    if molpct is not None:
        p["charges"]["mol_percent"] = amount_fact(molpct[0], molpct[1], "mol%", ev,
                                                      molpct[2] if len(molpct) > 2 else "original_extraction")
    if ref is not None:
        p["equivalence_reference"] = f(ref, ref, ev,
                                           method="inference",
                                           derivation="The explicitly reported equivalents/loading are relative to the substrate named in the same conditions.")
    return p


def c_c_classification(r, ev, explanation="Two net C-C connections are depicted for the intramolecular [5+2] cascade."):
    changes = []
    for idx in (1, 2):
        ch = {"change_id": f"cc{idx}",
              "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
              "before_order": 0, "after_order": 1, "change_type": "formed"}
        changes.append(f(None, ch, ev, method="structure_derivation", derivation=explanation))
    labels = [{"category_number": 2, "code": "C-C", "supporting_change_ids": ["cc1", "cc2"]}]
    r["classification"]["mapped_reaction_smiles"] = f()
    r["classification"]["bond_changes"] = changes
    r["classification"]["labels"] = f(None, labels, ev, method="structure_derivation", derivation=explanation)
    r["classification"]["coverage"] = inferred(
        "partial", ev,
        "Category-level net connectivity was read from the scheme; atom mapping, bond-order changes, migrating bond identity, and hydrogen provenance remain unresolved."
    )
    r["classification"]["limitations"] = [
        "No atom-complete mapping is asserted.",
        "Pinacol-rearrangement bond migration and pi-bond-order changes are not atom-located.",
        "Implicit-hydrogen changes are unresolved and no proton source is assigned.",
    ]
    r["classification"]["model_votes"] = []


def base_reaction(rid, ev, record_type="specific_experiment", categories=("other",),
                  location="", name="electrochemical ODI-[5+2] cycloaddition/pinacol rearrangement cascade"):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = f(list(categories), list(categories), ev)
    r["location"] = f(location, {"description": location}, ev)
    r["reaction_name"] = f(name, name, ev)
    r["transformation_type"] = inferred(
        "oxidative dearomatization-induced [5+2] cycloaddition/pinacol rearrangement",
        ev, "Normalized from the reaction name and arrows printed in the article."
    )
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": [
                       "No explicitly associated SI was supplied; detailed operations, workup, purification, characterization, and many structure assignments remain unavailable.",
                       "Dense primary structures and absolute stereochemical assignments require expert review.",
                       "Bond classification is category-level and partial rather than atom-complete.",
                   ]}
    r["completeness"] = "partial"
    return r


def set_electro_conditions(r, ev, substrate_mid, substrate_pid, substrate_amount=0.20,
                           current=5, solvent_mid="mol_mecn", solvent_name="MeCN",
                           solvent_volume=10, electrolyte_mid="mol_nbu4nbf4",
                           electrolyte_name="n-Bu4NBF4", additive_mid="mol_onp",
                           additive_name="o-nitrophenol", additive_amount=0.60,
                           additive_equiv=None, time=4, atmosphere="argon", scale_method="original_extraction"):
    r["reaction_scale"] = f(
        f"{substrate_amount:g} mmol {substrate_mid}",
        {"reference_molecule_id": substrate_mid, "amount": {"number": substrate_amount, "unit": "mmol"}},
        ev, method=scale_method,
    )
    r["solvent_conditions"]["composition"] = f(
        f"{solvent_name}; {electrolyte_name} (0.1 M)",
        [{"molecule_id": solvent_mid}, {"molecule_id": electrolyte_mid, "concentration": {"number": 0.1, "unit": "M"}}],
        ev, method=scale_method,
    )
    r["solvent_conditions"]["total_volume"] = f(
        f"{solvent_volume:g} mL", {"number": solvent_volume, "unit": "mL"}, ev, method=scale_method
    )
    r["solvent_conditions"]["concentration"] = f()
    r["solvent_conditions"]["mixture_ratio"] = f()
    r["solvent_conditions"]["concentration_basis"] = f()
    s = st.new_step("electrolysis", ev)
    s["participant_ids"] = [p["participant_id"] for p in r["participants"]]
    s["operation"] = f("constant-current electrolysis in an undivided cell", "constant-current electrolysis", ev, method=scale_method)
    s["temperature"] = f("rt", "rt", ev, method=scale_method)
    s["time"] = f(f"{time:g} h", {"number": time, "unit": "h"}, ev, method=scale_method)
    s["atmosphere"] = f(atmosphere, atmosphere, ev, method=scale_method)
    s["pressure"] = f()
    r["steps"] = [s]
    r["special_conditions"]["current"] = f(f"{current:g} mA", {"number": current, "unit": "mA"}, ev, method=scale_method)
    r["special_conditions"]["electrodes"] = f("graphite anode; graphite cathode", {"anode": "graphite", "cathode": "graphite"}, ev, method=scale_method)
    r["special_conditions"]["reactor"] = f("undivided cell", "undivided cell", ev, method=scale_method)
    r["procedure"]["original_text"] = f(
        "undivided cell, graphite anode, graphite cathode, constant current; substrate, additive, electrolyte, argon, solvent, rt",
        "Conditions transcribed from the applicable table/scheme footnote.", ev, method=scale_method,
    )


def yield_outcome(oid, product_mid, raw, number, ev, *, basis="isolated", scope=None):
    o = st.new_outcome(oid, "yield", [product_mid], ev)
    o["product_role"] = f("product", "main_product", ev)
    o["value"] = f(raw, {"number": number, "unit": "%"}, ev)
    o["yield_type"] = f("isolated yield" if basis == "isolated" else basis, basis, ev)
    o["measurement_basis"] = f("Isolated yields" if basis == "isolated" else None,
                                   basis if basis == "isolated" else None, ev if basis == "isolated" else ())
    o["scope"] = f(scope, scope, ev) if scope is not None else f()
    return o


def dr_outcome(oid, product_mid, raw, ev):
    o = st.new_outcome(oid, "selectivity", [product_mid], ev)
    o["product_role"] = f("product", "main_product", ev)
    o["value"] = f(raw, raw, ev)
    o["yield_type"] = f()
    o["measurement_basis"] = f()
    o["scope"] = f("dr", {"metric": "dr", "ratio_order": "major:minor; stereoisomer identities unassigned"}, ev)
    return o


def add_note(r, text, ev, method="inference", derivation=None):
    if method == "inference":
        r["notes"].append(f(None, text, ev, method=method,
                            derivation=derivation or "Extraction/review decision based on the cited figure and footnote."))
    else:
        r["notes"].append(f(text, text, ev, method=method))


def current_scope_record(family, sub_label, prod_label, yield_value, dr, page, generic_id):
    ev = [E[page]]
    sub_mid, prod_mid = f"mol_{sub_label}", f"mol_{prod_label}"
    rid = f"scheme_{family}_{prod_label}_electro"
    r = base_reaction(rid, ev, categories=("substrate_scope",),
                      location=f"Scheme {family}, product {prod_label}, electrochemical value")
    r["general_scheme_ids"] = [f"{DOC}:{generic_id}"]
    r["participants"] = [
        participant("substrate", sub_mid, "reactant", ev, amount=("0.20 mmol", 0.20, "mmol")),
        participant("product", prod_mid, "product", ev, product_role="main product"),
        participant("additive", "mol_onp", "additive", ev, amount=("0.60 mmol", 0.60, "mmol")),
        participant("electrolyte", "mol_nbu4nbf4", "reagent", ev),
        participant("solvent", "mol_mecn", "solvent", ev, volume=("10 mL", 10, "mL")),
    ]
    set_electro_conditions(r, ev, sub_mid, "substrate")
    r["outcomes"] = [yield_outcome("yield", prod_mid, f"{yield_value}%", yield_value, ev,
                                           scope={"method": "electrochemical", "figure_product": prod_label})]
    if dr:
        r["outcomes"].append(dr_outcome("dr", prod_mid, dr, ev))
    add_note(r, "The electrochemical value and parenthetical PIFA/HFIP comparator printed together in the figure were split into distinct experiment records.", ev)
    c_c_classification(r, ev)
    return r


def pifa_scope_record(family, sub_label, prod_label, yield_value, dr, page, generic_id, minutes):
    ev = [E[page]]
    sub_mid, prod_mid = f"mol_{sub_label}", f"mol_{prod_label}"
    rid = f"scheme_{family}_{prod_label}_pifa_comparator"
    r = base_reaction(rid, ev, categories=("substrate_scope", "control"),
                      location=f"Scheme {family}, product {prod_label}, parenthetical hypervalent-iodine comparison")
    r["general_scheme_ids"] = [f"{DOC}:{generic_id}"]
    r["participants"] = [
        participant("substrate", sub_mid, "reactant", ev,
                    amount=("0.20 mmol (scope charge inherited for the paired comparison)", 0.20, "mmol", "inherited_general_condition")),
        participant("product", prod_mid, "product", ev, product_role="main product"),
        participant("pifa", "mol_pifa", "reagent", ev, amount=("0.22 mmol", 0.22, "mmol")),
        participant("solvent", "mol_hfip", "solvent", ev, volume=("10 mL", 10, "mL")),
    ]
    r["reaction_scale"] = f("0.20 mmol paired scope substrate", {"reference_molecule_id": sub_mid, "amount": {"number": 0.20, "unit": "mmol"}}, ev, method="inherited_general_condition")
    r["solvent_conditions"]["composition"] = f("HFIP", [{"molecule_id": "mol_hfip"}], ev)
    r["solvent_conditions"]["total_volume"] = f("10 mL", {"number": 10, "unit": "mL"}, ev)
    s = st.new_step("hypervalent_iodine_oxidation", ev)
    s["participant_ids"] = [p["participant_id"] for p in r["participants"]]
    s["operation"] = f("PIFA oxidation in HFIP", "hypervalent iodine oxidation", ev)
    s["temperature"] = f("0 °C", {"number": 0, "unit": "°C"}, ev)
    s["time"] = f(f"{minutes} min", {"number": minutes, "unit": "min"}, ev)
    s["atmosphere"] = f()
    s["pressure"] = f()
    r["steps"] = [s]
    r["procedure"]["original_text"] = f(
        f"PIFA (0.22 mmol), HFIP (10 mL), 0 °C, {minutes} min",
        "Parenthetical comparator conditions", ev,
    )
    r["outcomes"] = [yield_outcome("yield", prod_mid, f"{yield_value}% (parenthetical)", yield_value, ev,
                                           scope={"method": "PIFA/HFIP comparator", "figure_product": prod_label})]
    if dr:
        r["outcomes"].append(dr_outcome("dr", prod_mid, dr, ev))
    add_note(r, "This is the parenthetical hypervalent-iodine comparison, not the electrochemical experiment reported beside it.", ev)
    c_c_classification(r, ev)
    return r


def main():
    bundle = json.loads((RUN / "submitted.json").read_text())
    molecules = []
    diagnostics = {}

    # Generic scheme families.
    molecules += [
        template_molecule("mol_generic_4", "generic vinylphenol 4", [E[2], E[3]], ["R1", "R2"], "source-drawn vinylphenol family 4"),
        template_molecule("mol_generic_5", "generic bicyclo[3.2.1]octane 5", [E[2], E[3]], ["R1", "R2"], "source-drawn product family 5"),
        template_molecule("mol_generic_7", "generic ethynylphenol 7", [E[3]], ["R1", "R2", "R3"], "source-drawn ethynylphenol family 7"),
        template_molecule("mol_generic_8", "generic bicyclo[3.2.1]octane 8", [E[3]], ["R1", "R2", "R3"], "source-drawn product family 8"),
        template_molecule("mol_generic_9", "generic allenylphenol 9", [E[3], E[4]], ["R1", "R2"], "source-drawn allenylphenol family 9"),
        template_molecule("mol_generic_10", "generic bicyclo[3.2.1]octane 10", [E[3], E[4]], ["R1", "R2"], "source-drawn product family 10"),
    ]

    # Scope identities.  Labels are retained independently; no structures are merged by scaffold/name.
    for prefix, letters, page, kind in [("4", "abcdefghijklmnopqrstuvw", 3, "neutral"),
                                         ("5", "abcdefghijklmnopqrstuvw", 3, "product"),
                                         ("7", "abcdefghijklmnopqrs", 3, "neutral"),
                                         ("8", "abcdefghijklmnopqrs", 3, "product"),
                                         ("9", "abcdefghijklmnop", 4, "allene_substrate"),
                                         ("10", "abcdefghijklmnop", 4, "product")]:
        for letter in letters:
            label = prefix + letter
            ev = [E[2], E[3]] if label in ("4a", "5a") else [E[page]]
            molecules.append(unresolved_molecule(f"mol_{label}", label, ev, kind=kind))

    # Referenced Yamamura example and its depicted intermediate.
    for label, kind in [("1", "neutral"), ("2", "neutral"), ("3", "product")]:
        molecules.append(unresolved_molecule(f"mol_lit_{label}", f"Yamamura compound {label}", [E[2]], kind=kind,
                                             note="Literature-example structure in Scheme 1b; not a current-paper scope compound."))

    common = [
        ("mol_onp", "o-nitrophenol", "2-nitrophenol", "O=[N+]([O-])c1ccccc1O", [E[2], E[3], E[4]]),
        ("mol_mecn", "MeCN", "acetonitrile", "CC#N", [E[2], E[3], E[4]]),
        ("mol_dce", "DCE", "1,2-dichloroethane", "ClCCCl", [E[2]]),
        ("mol_hfip", "HFIP", "1,1,1,3,3,3-hexafluoro-2-propanol", "OC(C(F)(F)F)C(F)(F)F", [E[2], E[3], E[4]]),
        ("mol_meoh", "MeOH", "methanol", "CO", [E[2]]),
        ("mol_hcl", "HCl", "hydrogen chloride", "Cl", [E[2]]),
        ("mol_tfa", "TFA", "trifluoroacetic acid", "O=C(O)C(F)(F)F", [E[2], E[3]]),
        ("mol_ac2o", "Ac2O", "acetic anhydride", "CC(=O)OC(C)=O", [E[2]]),
        ("mol_dhp", "DHP", "3,4-dihydro-2H-pyran", "C1=CCCOC1", [E[3]]),
        ("mol_ch2cl2", "CH2Cl2", "dichloromethane", "ClCCl", [E[3]]),
        ("mol_na2co3", "Na2CO3", "sodium carbonate", "[Na+].[Na+].[O-]C([O-])=O", [E[3]]),
        ("mol_na3po4", "Na3PO4", "trisodium phosphate", "[Na+].[Na+].[Na+].[O-]P([O-])([O-])=O", [E[3]]),
    ]
    for args in common:
        m, diag = concrete_molecule(*args)
        molecules.append(m)
        diagnostics[m["molecule_id"]] = diag

    for mid, label, ev in [
        ("mol_nbu4nbf4", "n-Bu4NBF4", [E[1], E[2], E[3], E[4]]),
        ("mol_nbu4npf6", "n-Bu4NPF6", [E[2]]),
        ("mol_pifa", "PIFA", [E[2], E[3], E[4]]),
        ("mol_ptsoh_h2o", "p-TsOH·H2O", [E[3]]),
    ]:
        molecules.append(unresolved_molecule(mid, label, ev, note="Named reagent retained without a forced complete ionic/hydrate representation."))

    reactions = []

    # Current generic schemes for the three substrate classes.
    for family, sub, prod, page, minutes in [
        ("vinyl", "mol_generic_4", "mol_generic_5", 3, 4),
        ("ethynyl", "mol_generic_7", "mol_generic_8", 3, 4),
        ("allenyl", "mol_generic_9", "mol_generic_10", 4, 4),
    ]:
        rid = f"generic_{family}_electrochemical"
        generic_ev = [E[1], E[2], E[3]] if family == "vinyl" else [E[1], E[page]]
        r = base_reaction(rid, generic_ev, record_type="generic_scheme", categories=("other",),
                          location=f"Scheme {'2' if family=='vinyl' else '3' if family=='ethynyl' else '4'} general arrow")
        r["participants"] = [
            participant("substrate_family", sub, "reactant", [E[page]], amount=("0.20 mmol", 0.20, "mmol")),
            participant("product_family", prod, "product", [E[page]], product_role="product family"),
            participant("additive", "mol_onp", "additive", [E[page]], amount=("0.60 mmol", 0.60, "mmol")),
            participant("electrolyte", "mol_nbu4nbf4", "reagent", [E[page]]),
            participant("solvent", "mol_mecn", "solvent", [E[page]], volume=("10 mL", 10, "mL")),
        ]
        set_electro_conditions(r, [E[page]], sub, "substrate_family")
        r["outcomes"] = []
        if family == "vinyl":
            conflict_value = f(
                ["41–95% (abstract)", "41–96% (Scheme 1c)"], None, [E[1], E[2]],
                status="conflict",
                alternatives=[
                    f("41–95%", {"minimum": 41, "maximum": 95, "unit": "%", "raw_text": "41–95%"}, [E[1]]),
                    f("41–96%", {"minimum": 41, "maximum": 96, "unit": "%", "raw_text": "41–96%"}, [E[2]]),
                ],
            )
            o = st.new_outcome("reported_method_yield_range", "yield", [prod], [E[1], E[2]])
            o["product_role"] = f("product family", "product_family", [E[1], E[2]])
            o["value"] = conflict_value
            o["yield_type"] = f()
            o["measurement_basis"] = f()
            o["scope"] = f("method-level range", {"scope": "method summary rather than one experiment"}, [E[1], E[2]])
            r["outcomes"].append(o)
            r["review"]["unresolved"].append("The abstract reports 41–95% whereas Scheme 1c prints 41–96%; both values are retained without choosing one.")
        c_c_classification(r, [E[page]])
        reactions.append(r)

    # Yamamura cited experiment.
    r = base_reaction("literature_yamamura_1991", [E[2]], record_type="literature_example",
                      categories=("other",), location="Scheme 1b, cited Yamamura work (1991)",
                      name="anodic oxidation of vinylphenol")
    r["participants"] = [
        participant("substrate", "mol_lit_1", "reactant", [E[2]]),
        participant("intermediate", "mol_lit_2", "unknown", [E[2]]),
        participant("product", "mol_lit_3", "product", [E[2]], product_role="literature product"),
        participant("medium", "mol_ac2o", "solvent", [E[2]]),
    ]
    r["solvent_conditions"]["composition"] = f("Ac2O", [{"molecule_id": "mol_ac2o"}], [E[2]])
    s = st.new_step("anodic_oxidation", [E[2]])
    s["participant_ids"] = [p["participant_id"] for p in r["participants"]]
    s["operation"] = f("CCE at 2.5 mA", "constant-current electrolysis", [E[2]])
    r["steps"] = [s]
    r["special_conditions"]["current"] = f("2.5 mA", {"number": 2.5, "unit": "mA"}, [E[2]])
    r["outcomes"] = [yield_outcome("yield", "mol_lit_3", "59%", 59, [E[2]], basis="unreported",
                                           scope={"source": "cited literature example"})]
    c_c_classification(r, [E[2]])
    reactions.append(r)

    # Prior PIFA generic method shown as previous work.
    r = base_reaction("literature_previous_pifa_generic", [E[2]], record_type="literature_example",
                      categories=("other",), location="Scheme 1c, authors' previous work",
                      name="iodine(III)-mediated ODI-[5+2] cascade")
    r["participants"] = [
        participant("substrate_family", "mol_generic_4", "reactant", [E[2]]),
        participant("product_family", "mol_generic_5", "product", [E[2]], product_role="product family"),
        participant("pifa", "mol_pifa", "reagent", [E[2]]),
        participant("hfip", "mol_hfip", "solvent", [E[2]]),
    ]
    r["solvent_conditions"]["composition"] = f("HFIP", [{"molecule_id": "mol_hfip"}], [E[2]])
    s = st.new_step("pifa_oxidation", [E[2]])
    s["participant_ids"] = [p["participant_id"] for p in r["participants"]]
    s["operation"] = f("PIFA, HFIP", "hypervalent iodine oxidation", [E[2]])
    s["temperature"] = f("0 °C", {"number": 0, "unit": "°C"}, [E[2]])
    r["steps"] = [s]
    o = st.new_outcome("yield_range", "yield", ["mol_generic_5"], [E[2]])
    o["product_role"] = f("product family", "main_product", [E[2]])
    o["value"] = f("35–85% yields", {"minimum": 35, "maximum": 85, "unit": "%", "raw_text": "35–85% yields"}, [E[2]])
    o["yield_type"] = f()
    o["measurement_basis"] = f()
    o["scope"] = f("our previous work", {"source": "authors' cited previous work"}, [E[2]])
    r["outcomes"] = [o]
    c_c_classification(r, [E[2]])
    reactions.append(r)

    # Proposed cascade mechanism, not a performed experiment.
    r = base_reaction("proposed_odi_pinacol_mechanism", [E[2]], record_type="proposed_mechanism",
                      categories=("mechanistic_study",), location="Scheme 1c, proposed ionic intermediates and pinacol rearrangement")
    r["general_scheme_ids"] = [f"{DOC}:generic_vinyl_electrochemical"]
    r["participants"] = [
        participant("substrate_family", "mol_generic_4", "reactant", [E[2]]),
        participant("product_family", "mol_generic_5", "product", [E[2]], product_role="product family"),
    ]
    add_note(r, "The phenoxy-cation/dearomatized intermediate and pinacol rearrangement arrows are author-proposed intermediates, not separate performed experiments.", [E[2]])
    c_c_classification(r, [E[2]])
    reactions.append(r)

    # Table 1 optimization rows.
    deviations = {
        1: ("none", 90), 2: ("DCE as the solvent", 54), 3: ("HFIP as the solvent", 25),
        4: ("MeOH as the solvent", 15), 5: ("n-Bu4NPF6 as the electrolyte", 85),
        6: ("2 mA", 48), 7: ("10 mA", 72), 8: ("HCl (2 M, aqueous) as the additive", 25),
        9: ("TFA as the additive", 84), 10: ("2 equiv of o-nitrophenol", 70),
        11: ("4a (15 mmol); 50 h", 89),
    }
    solvent_map = {2: ("mol_dce", "DCE"), 3: ("mol_hfip", "HFIP"), 4: ("mol_meoh", "MeOH")}
    for entry, (dev, y) in deviations.items():
        ev = [E[2]]
        r = base_reaction(f"table1_entry_{entry:02d}", ev, categories=("optimization",) + (("scale_up",) if entry == 11 else ()),
                          location=f"Table 1 entry {entry}")
        r["general_scheme_ids"] = [f"{DOC}:generic_vinyl_electrochemical"]
        sub_amount = 15.0 if entry == 11 else 0.20
        solvent_mid, solvent_name = solvent_map.get(entry, ("mol_mecn", "MeCN"))
        electrolyte_mid = "mol_nbu4npf6" if entry == 5 else "mol_nbu4nbf4"
        electrolyte_name = "n-Bu4NPF6" if entry == 5 else "n-Bu4NBF4"
        current = 2 if entry == 6 else 10 if entry == 7 else 5
        additive_mid = "mol_hcl" if entry == 8 else "mol_tfa" if entry == 9 else "mol_onp"
        additive_name = "HCl (2 M, aqueous)" if entry == 8 else "TFA" if entry == 9 else "o-nitrophenol"
        additive_amount = None if entry in (8, 9, 10, 11) else 0.60
        add_equiv = 2 if entry == 10 else 3 if entry not in (8, 9) else None
        r["participants"] = [
            participant("substrate", "mol_4a", "reactant", ev, amount=(f"{sub_amount:g} mmol", sub_amount, "mmol")),
            participant("product", "mol_5a", "product", ev, product_role="main product"),
            participant("additive", additive_mid, "additive", ev,
                        amount=("0.60 mmol", additive_amount, "mmol") if additive_amount is not None else None,
                        equiv=(f"{add_equiv:g} equiv", add_equiv) if add_equiv is not None else None,
                        ref="mol_4a" if add_equiv is not None else None),
            participant("electrolyte", electrolyte_mid, "reagent", ev),
            participant("solvent", solvent_mid, "solvent", ev, volume=("10 mL", 10, "mL", "inherited_general_condition") if entry in (2,3,4) else ("10 mL", 10, "mL")),
        ]
        set_electro_conditions(r, ev, "mol_4a", "substrate", substrate_amount=sub_amount,
                               current=current, solvent_mid=solvent_mid, solvent_name=solvent_name,
                               electrolyte_mid=electrolyte_mid, electrolyte_name=electrolyte_name,
                               additive_mid=additive_mid, additive_name=additive_name,
                               additive_amount=additive_amount or 0, additive_equiv=add_equiv,
                               time=50 if entry == 11 else 4,
                               scale_method="inherited_general_condition" if entry != 1 else "original_extraction")
        r["procedure"]["modifications"] = f(dev, dev, ev)
        if entry == 8:
            add_note(r, "The source identifies aqueous 2 M HCl as the replacement additive; its added volume/amount is unreported.", ev)
        if entry == 9:
            add_note(r, "The source identifies TFA as the replacement additive; its added amount is unreported.", ev)
        if entry == 11:
            add_note(r, "The scale-up row changes the substrate charge to 15 mmol and the footnote changes the duration to 50 h; other standard-condition quantities are not arithmetically scaled here.", ev)
        r["outcomes"] = [yield_outcome("yield", "mol_5a", f"{y}%", y, ev, scope={"table_entry": entry})]
        c_c_classification(r, ev)
        reactions.append(r)

    # Table 1 entry-1 PIFA comparison is an independent experiment.
    r = pifa_scope_record("1", "4a", "5a", 78, None, 2, "generic_vinyl_electrochemical", 20)
    r["reaction_id"] = f"{DOC}:table1_entry_01_pifa_comparator"
    r["location"] = f("Table 1 entry 1 parenthetical yield and footnote c", {"table": 1, "entry": 1, "comparison": "PIFA"}, [E[2]])
    r["experiment_categories"] = f(["optimization", "control"], ["optimization", "control"], [E[2]])
    reactions.append(r)

    scheme2 = {
        "b": (93, 80, None, None), "c": (92, 85, None, None), "d": (88, 66, None, None),
        "e": (75, 64, None, None), "f": (78, 70, None, None), "g": (89, 60, None, None),
        "h": (94, 74, None, None), "i": (93, 70, None, None), "j": (86, 67, None, None),
        "k": (70, 15, ">20:1 dr", ">20:1 dr"), "l": (76, 35, ">20:1 dr", ">20:1 dr"),
        "m": (96, 82, ">20:1 dr", ">20:1 dr"), "n": (96, 76, "10:1 dr", "10:1 dr"),
        "o": (88, 70, ">20:1 dr", ">20:1 dr"), "p": (88, 84, None, None),
        "q": (85, 70, "15:1 dr", "15:1 dr"), "r": (81, 53, None, None),
        "s": (86, 60, "1.5:1 dr", "1.5:1 dr"), "t": (75, 63, None, None),
        "u": (45, 37, None, None), "v": (71, 58, None, None), "w": (84, 60, None, None),
    }
    scheme3 = {
        "a": (88, 65, None, None), "b": (85, 45, None, None), "c": (84, 38, None, None),
        "d": (79, 34, None, None), "e": (90, 42, None, None), "f": (87, 37, None, None),
        "g": (74, 13, ">20:1 dr", ">20:1 dr"), "h": (70, 34, ">20:1 dr", ">20:1 dr"),
        "i": (90, 50, ">20:1 dr", ">20:1 dr"), "j": (91, 51, "1.8:1 dr", "1.8:1 dr"),
        "k": (87, 41, None, None), "l": (91, 68, None, None), "m": (76, 45, None, None),
        "n": (80, 39, None, None), "o": (95, 48, None, None), "p": (80, 0, None, None),
        "q": (75, 0, None, None), "r": (72, 0, None, None), "s": (81, 0, None, None),
    }
    scheme4 = {
        "a": (92, 25, None, None), "b": (91, 16, None, None), "c": (89, 20, None, None),
        "d": (92, 23, None, None), "e": (74, 10, None, None), "f": (71, 8, None, None),
        "g": (85, 12, None, None), "h": (91, 18, None, None), "i": (87, 20, None, None),
        "j": (76, 11, None, None), "k": (85, 21, None, None), "l": (85, 7, ">20:1 dr", ">20:1 dr"),
        "m": (90, 18, ">20:1 dr", ">20:1 dr"), "n": (87, 17, "4.6:1 dr", "4.6:1 dr"),
        "o": (86, 13, None, None), "p": (41, 6, None, None),
    }
    transcribed = {"scheme2": scheme2, "scheme3": scheme3, "scheme4": scheme4}

    for family, data, sub_prefix, prod_prefix, page, gid, minutes in [
        ("2", scheme2, "4", "5", 3, "generic_vinyl_electrochemical", 20),
        ("3", scheme3, "7", "8", 3, "generic_ethynyl_electrochemical", 40),
        ("4", scheme4, "9", "10", 4, "generic_allenyl_electrochemical", 20),
    ]:
        for letter, (electro_y, pifa_y, electro_dr, pifa_dr) in data.items():
            reactions.append(current_scope_record(family, sub_prefix + letter, prod_prefix + letter,
                                                  electro_y, electro_dr, page, gid))
            if family == "2" and letter == "v":
                ev = [E[3]]
                r = base_reaction("scheme_2_5v_previous_multistep", ev,
                                  categories=("substrate_scope", "control"),
                                  location="Scheme 2, 5v, previous conditions (58%, 2 steps)")
                r["general_scheme_ids"] = [f"{DOC}:{gid}"]
                r["participants"] = [
                    participant("substrate", "mol_4v", "reactant", ev),
                    participant("product", "mol_5v", "product", ev, product_role="main product"),
                    participant("dhp", "mol_dhp", "reagent", ev),
                    participant("ptsoh", "mol_ptsoh_h2o", "reagent", ev),
                    participant("ch2cl2", "mol_ch2cl2", "solvent", ev),
                    participant("pifa", "mol_pifa", "reagent", ev),
                    participant("na2co3", "mol_na2co3", "additive", ev),
                    participant("na3po4", "mol_na3po4", "additive", ev),
                    participant("hfip", "mol_hfip", "solvent", ev),
                    participant("tfa", "mol_tfa", "reagent", ev),
                ]
                for sid, op, temp, pids in [
                    ("protection", "DHP, p-TsOH·H2O, CH2Cl2", "rt", ["substrate", "dhp", "ptsoh", "ch2cl2"]),
                    ("oxidation", "PIFA, Na2CO3, Na3PO4, HFIP", {"number": 0, "unit": "°C"}, ["pifa", "na2co3", "na3po4", "hfip"]),
                    ("tfa_treatment", "then TFA", None, ["tfa", "product"]),
                ]:
                    s = st.new_step(sid, ev)
                    s["participant_ids"] = pids
                    s["operation"] = f(op, op, ev)
                    if temp == "rt": s["temperature"] = f("rt", "rt", ev)
                    elif isinstance(temp, dict): s["temperature"] = f("0 °C", temp, ev)
                    r["steps"].append(s)
                r["outcomes"] = [yield_outcome("overall_yield", "mol_5v", "58% (2 steps)", 58, ev,
                                                       scope={"span": "4v to 5v", "reported_as": "2 steps"})]
                r["procedure"]["original_text"] = f("1. DHP, p-TsOH·H2O, CH2Cl2, rt; 2. PIFA, Na2CO3, Na3PO4, HFIP, 0 °C; then TFA", "previous conditions printed for 5v", ev)
                add_note(r, "The 58% value is retained as an overall two-step yield from 4v to 5v, not assigned to the final TFA operation alone.", ev)
                r["review"]["unresolved"].append("Reagent quantities and all step durations for the previous 5v sequence are unreported in the supplied main text.")
                c_c_classification(r, ev)
                reactions.append(r)
            else:
                reactions.append(pifa_scope_record(family, sub_prefix + letter, prod_prefix + letter,
                                                   pifa_y, pifa_dr, page, gid, minutes))

    # Source context: non-reaction examples and prose-only observations.
    bundle["source_context"] = [
        {"context_id": f"{DOC}:scheme1_natural_products", "document_id": DOC,
         "kind": "structural_examples", "evidence_ids": [E[2]],
         "content": "Scheme 1a depicts oridonin, platensimycin, grayanotoxin I, kopsifoline A, aconitine, and nominine as representative natural products; these are not experiments performed in this paper."},
        {"context_id": f"{DOC}:scope_discussion", "document_id": DOC,
         "kind": "author_discussion", "evidence_ids": [E[2], E[3], E[4]],
         "content": "Author discussion includes decomposition of a catechol-type substrate, double-bond migration for 7k/8k, and exclusive Diels–Alder products 8p–s; no standalone spectral peak table appears in the supplied five pages."},
        {"context_id": f"{DOC}:si_notice", "document_id": DOC,
         "kind": "availability_notice", "evidence_ids": [E[4]],
         "content": "The article advertises Supporting Information, but no SI file was supplied or explicitly associated with this workspace."},
    ]

    bundle["molecules"] = molecules
    bundle["reactions"] = reactions
    (RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
    (RUN / "raw_model/stereo_diagnostics.json").write_text(json.dumps(diagnostics, ensure_ascii=False, indent=2) + "\n")
    (RUN / "raw_model/scope_transcription.json").write_text(json.dumps(transcribed, ensure_ascii=False, indent=2) + "\n")
    mention_terms = ["4a", "5a", "4v", "5v", "7k", "8k", "8p", "8s", "9a", "10a"]
    mentions = st.source_mentions(json.loads((RUN / "dataset.json").read_text()), mention_terms)
    (RUN / "raw_model/source_mentions.json").write_text(json.dumps(mentions, ensure_ascii=False, indent=2) + "\n")
    by_page = {p: [r["reaction_id"] for r in reactions if E[p] in r["evidence_ids"]] for p in range(1, 6)}
    coverage = {
        "document_id": DOC,
        "semantic_accuracy_verified": False,
        "pages": [
            {"page": 1, "state": "attempted", "regions": [
                {"region": "title, abstract, and graphical abstract", "finding": "Current electrochemical ODI-[5+2]/pinacol cascade overview for vinyl-, ethynyl-, and allenylphenols; the abstract's 41–95% range differs from Scheme 1c's 41–96% range and is retained as a conflict.", "evidence_ids": [E[1], T[1]]},
                {"region": "two-column introduction and lower-right optimization lead-in", "finding": "Background/cited work plus the start of the current standard-condition description; no additional discrete current experiment beyond records continued on page 2.", "evidence_ids": [E[1], T[1]]},
            ], "record_ids": by_page[1], "unresolved": ["Graphical-abstract structures are generic and were not treated as concrete molecular identities."]},
            {"page": 2, "state": "attempted", "regions": [
                {"region": "Scheme 1a", "finding": "Representative natural-product structures only; retained as source context, not current-paper experiments.", "evidence_ids": [E[2], T[2]]},
                {"region": "Scheme 1b", "finding": "One cited Yamamura constant-current experiment (1→3, 59%) plus depicted intermediate 2; classified as literature_example.", "evidence_ids": [E[2], T[2]]},
                {"region": "Scheme 1c", "finding": "Authors' previous PIFA method, current electrochemical generic transformation, and author-proposed oxidative-dearomatization/pinacol mechanism kept as distinct record types.", "evidence_ids": [E[2], T[2]]},
                {"region": "Table 1 and footnotes", "finding": "All 11 electrochemical optimization rows plus the independent entry-1 PIFA/HFIP comparison were recorded; entry 11 retains 15 mmol and 50 h.", "evidence_ids": [E[2], T[2]]},
                {"region": "right-column scope discussion", "finding": "Narrative results for 5b–w, including 5v/5w previous-method comparisons and the catechol decomposition observation; scope arrows are on page 3.", "evidence_ids": [E[2], T[2]]},
            ], "record_ids": by_page[2], "unresolved": ["No SI procedures were supplied; workup and purification details are missing.", "The catechol-type decomposition statement lacks a compound label, amount, and quantified outcome, so it remains source context rather than a fabricated experiment."]},
            {"page": 3, "state": "attempted", "regions": [
                {"region": "Scheme 2 and footnotes", "finding": "Twenty-two vinylphenol electrochemical scope experiments (5b–w) and twenty-two distinct previous/PIFA comparison experiments; 5v previous conditions are a separate two-step overall-yield record.", "evidence_ids": [E[3], T[3]]},
                {"region": "Scheme 3 and footnotes", "finding": "Nineteen ethynylphenol electrochemical scope experiments (8a–s) and nineteen distinct PIFA comparisons; 8k migration and 8p–s Diels–Alder outcome class are preserved in context/labels.", "evidence_ids": [E[3], T[3]]},
                {"region": "left-column discussion", "finding": "Author interpretation of ethynyl and allenyl scope; not converted into extra duplicate experiments.", "evidence_ids": [E[3], T[3]]},
            ], "record_ids": by_page[3], "unresolved": ["Dense scope drawings, absolute configurations, exact atom mappings, and dr measurement methods remain unresolved without SI."]},
            {"page": 4, "state": "attempted", "regions": [
                {"region": "Scheme 4 and footnotes", "finding": "Sixteen allenylphenol electrochemical scope experiments (10a–p) and sixteen distinct PIFA comparisons, with all printed yields and dr values retained as product-targeted outcomes.", "evidence_ids": [E[4], T[4]]},
                {"region": "conclusion and Associated Content", "finding": "Conclusion only; SI availability notice is recorded, but no SI association is inferred.", "evidence_ids": [E[4], T[4]]},
                {"region": "author information and references start", "finding": "Bibliographic material; no additional current-work experiment.", "evidence_ids": [E[4], T[4]]},
            ], "record_ids": by_page[4], "unresolved": ["Allene axial stereogenicity and product absolute configurations are unassigned in the supplied main-text figures; representations remain not_encoded/unresolved as applicable."]},
            {"page": 5, "state": "attempted", "regions": [
                {"region": "two-column references", "finding": "References only; cited papers were not converted into additional current-paper experiments.", "evidence_ids": [E[5], T[5]]},
            ], "record_ids": by_page[5], "unresolved": []},
        ],
    }
    (RUN / "coverage").mkdir(exist_ok=True)
    (RUN / "coverage/pages_0001_0005.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")
    tasks = json.loads((RUN / "tasks.json").read_text())
    for task in tasks:
        if task.get("document_id") == DOC and 1 <= int(task["evidence_id"].split(":page", 1)[1].split(":", 1)[0]) <= 5:
            task["state"] = "attempted"
    (RUN / "tasks.json").write_text(json.dumps(tasks, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"molecules": len(molecules), "reactions": len(reactions),
                      "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in reactions),
                      "stereo_diagnostics": len(diagnostics)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
