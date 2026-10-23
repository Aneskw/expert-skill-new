#!/usr/bin/env python3
"""Source-grounded transcription for doc_d83675624ef94ad8, PDF pages 1-6.

This is a session-local extraction artifact. It does not alter the extraction
tools, schemas, taxonomy, tests, or source assets.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_d83675624ef94ad8"
RUN = ROOT / "runs/current"


def E(page: int):
    return [f"{DOC}:page{page}:image", f"{DOC}:page{page}:text"]


def f(raw=None, normalized=None, ev=(), method="original_extraction", derivation=None,
      status=None):
    return st.fact(raw, normalized, ev, method=method, derivation=derivation,
                   status=status)


def inferred(normalized, ev, derivation):
    return f(None, normalized, ev, method="inference", derivation=derivation)


def set_stereo(m, ev, chirality, geometry, representation, *, assignment=None,
               axial=None, mixture=None, double_bonds=None, note=None):
    stereo = m["stereochemistry"]
    stereo["chirality_status"] = inferred(chirality, ev, note or "Assigned conservatively from the inspected source drawing; no unreported absolute configuration was added.")
    stereo["geometry_status"] = inferred(geometry, ev, note or "Assigned conservatively from the inspected source drawing.")
    stereo["representation_status"] = inferred(representation, ev, "Reports whether the source stereochemical information is encoded in a concrete isomeric SMILES.")
    if assignment is not None:
        stereo["assignment"] = f(assignment, assignment, ev)
    if axial is not None:
        stereo["axial"] = f(axial, axial, ev)
    if mixture is not None:
        stereo["mixture"] = f(mixture, mixture, ev)
    if double_bonds is not None:
        stereo["double_bonds"] = f(double_bonds, double_bonds, ev)


def source_molecule(mid, label, name, page, *, stereo_kind="unresolved", note=None):
    ev = E(page)
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    m["names"]["original"] = f(name, name, ev)
    m["names"]["english"] = f(name, name, ev)
    m["identity_note"] = f(note, note, ev) if note else f(
        "Source label and drawing inspected; complete connectivity was not encoded as SMILES.",
        "Source label and drawing inspected; complete connectivity was not encoded as SMILES.", ev)
    if stereo_kind == "substrate_allene":
        set_stereo(
            m, ev, "unspecified", "not_applicable", "not_encoded",
            axial="The source family drawing has a terminal CH2 allene; no axial descriptor is assigned. Carbon-centered configuration is not specified.",
            note="The allenol drawing may contain a stereogenic alcohol-bearing carbon, but no absolute/relative assignment is printed; unspecified is not treated as racemic.")
    elif stereo_kind == "substrate_alkene":
        set_stereo(
            m, ev, "unspecified", "not_applicable", "not_encoded",
            double_bonds="The source family drawing shows a terminal alkene; E/Z is not applicable. Carbon-centered configuration is not specified.",
            note="No absolute/relative configuration is printed for the allylic alcohol substrate; unspecified is not treated as racemic.")
    elif stereo_kind == "cis_product":
        set_stereo(
            m, ev, "partially_specified", "not_applicable", "not_encoded",
            assignment="cis-fused ring-junction relative stereochemistry as drawn; absolute configuration not assigned",
            note="The source drawing/prose specifies relative cis fusion but not absolute R/S; the complex structure is retained without a potentially wrong stereochemical SMILES.")
    elif stereo_kind == "drawn_product":
        set_stereo(
            m, ev, "partially_specified", "unresolved", "not_encoded",
            assignment="relative stereochemical bonds are shown in the source drawing; absolute configuration is not assigned",
            note="The derivatization drawing contains stereochemical bonds, but full center-by-center and alkene assignments were not safely encodable from the supplied main-text image alone.")
    else:
        set_stereo(
            m, ev, "unresolved", "unresolved", "unresolved",
            note="The source identifies the species, but the supplied representation is generic, intermediate-like, or insufficient for a complete stereochemical graph.")
    return m


def smiles_molecule(mid, name, smiles, page):
    ev = E(page)
    diag = st.stereo_inspect(smiles)
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["representation_type"] = "smiles"
    m["occurrences"][0]["compound_label"] = f(name, name, ev)
    m["names"]["original"] = f(name, name, ev)
    m["names"]["english"] = f(name, name, ev)
    m["smiles"] = f(name, diag["canonical_isomeric_smiles"], ev,
                    method="structure_derivation",
                    derivation="Candidate SMILES transcribed from the explicit reagent/solvent identity and checked by RDKit; parsing does not certify source identity.")
    potential = diag["potential_stereo"]
    atom_potential = [x for x in potential if str(x["type"]).startswith("Atom_")]
    bond_potential = [x for x in potential if str(x["type"]).startswith("Bond_")]
    chirality = "unspecified" if atom_potential else "achiral"
    geometry = "unspecified" if bond_potential else "not_applicable"
    set_stereo(m, ev, chirality, geometry,
               "partial" if potential else "not_applicable",
               note="RDKit graph diagnostics were run on the candidate SMILES; absent toolkit flags do not establish source stereochemistry.")
    m["stereochemistry"]["isomeric_smiles"] = f(
        name, diag["canonical_isomeric_smiles"], ev,
        method="structure_derivation",
        derivation="Canonical isomeric SMILES generated by standard stereo diagnostics from the candidate graph.")
    m["stereochemistry"]["identity_key"] = f(
        None, diag["stereo_graph_key"], ev, method="structure_derivation",
        derivation="Stereo graph key generated by standard stereo diagnostics; it is not a sample identity.")
    m["stereochemistry"]["atom_centers"] = f(
        None, diag["atom_centers"], ev, method="structure_derivation",
        derivation="RDKit potential/assigned atom-center diagnostics from the candidate SMILES.")
    m["stereochemistry"]["double_bonds"] = f(
        None, diag["double_bonds"], ev, method="structure_derivation",
        derivation="RDKit double-bond diagnostics from the candidate SMILES.")
    m["identity_note"] = f(
        "Candidate graph for an explicitly named reagent or solvent; source drawing identity still needs expert review.",
        "Candidate graph for an explicitly named reagent or solvent; source drawing identity still needs expert review.", ev)
    return m, diag


def participant(pid, mid, role, ev, *, stage="reaction", mol_percent=None,
                reference=None, product_role=None):
    p = st.new_participant(pid, mid, role, ev, stage=stage)
    if product_role:
        p["product_role"] = f(product_role, product_role, ev)
    if mol_percent is not None:
        p["charges"]["mol_percent"] = f(
            f"{mol_percent} mol%", {"number": mol_percent, "unit": "mol%"}, ev)
        p["equivalence_reference"] = f(reference, reference, ev)
    return p


def manual_classification(r, labels, ev, detail):
    changes = []
    label_items = {}
    for i, code in enumerate(labels, 1):
        if code == "C-C":
            atoms, before, after, kind = [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], 0, 1, "formed"
        elif code == "C-N":
            atoms, before, after, kind = [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "N"}], 0, 1, "formed"
        elif code == "C-Br":
            atoms, before, after, kind = [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "Br"}], 0, 1, "formed"
        elif code == "C-S":
            atoms, before, after, kind = [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "S"}], 0, 1, "formed"
        elif code == "N-H":
            atoms, before, after, kind = [{"atom_map": None, "element": "N"}, {"atom_map": None, "element": "H"}], 0, 1, "hydrogen_gain"
        elif code == "C-H":
            atoms, before, after, kind = [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}], 0, 1, "hydrogen_gain"
        elif code == "DB-C-C":
            atoms, before, after, kind = [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], 1, 2, "order_increased"
        else:
            raise ValueError(code)
        change = {"change_id": f"b{i}", "atoms": atoms, "before_order": before,
                  "after_order": after, "change_type": kind}
        changes.append(change)
        tax = st.bond_label(atoms[0]["element"], atoms[1]["element"], after)
        label_items.setdefault(tax["code"], {**tax, "supporting_change_ids": []})["supporting_change_ids"].append(change["change_id"])
    r["classification"]["bond_changes"] = [
        f(None, c, ev, method="structure_derivation",
          derivation=detail + " Exact atom identities/maps and coproduct/hydrogen balance remain unresolved.")
        for c in changes
    ]
    r["classification"]["labels"] = f(
        None, list(label_items.values()), ev, method="structure_derivation",
        derivation=detail + " Labels cover only visually supported formed/promoted bond classes.")
    r["classification"]["coverage"] = inferred(
        "partial", ev,
        "Source drawings support bond-class assignment, but exact atom mapping, implicit-H balance, leaving groups, and coproducts are incomplete.")
    r["classification"]["limitations"] = [
        "No complete atom-mapped reaction was asserted.",
        "Implicit hydrogen and omitted coproduct balance remain unresolved.",
        "Classification is a structure-derived partial net-change description, not an elementary mechanism.",
    ]


def base_reaction(rid, page, categories, location, name, transformation, *,
                  record_type="specific_experiment", evidence=None):
    ev = evidence or E(page)
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = inferred(categories, ev, "Record category inferred from the source region (scope, derivatization, optimization, or proposed mechanism).")
    r["location"] = f(location, {"pdf_page": page, "region": location}, ev)
    r["reaction_name"] = f(name, name, ev)
    r["transformation_type"] = inferred(transformation, ev, "Transformation type summarized from the source arrow and surrounding prose.")
    r["review"] = {
        "status": "needs_review", "method": "model", "reviewer": None,
        "checked_fields": [],
        "unresolved": [
            "Machine transcription remains needs_review.",
            "Associated SI was mentioned by the paper but was not supplied; missing charges, scale, workup, purification, and measurement details were not inferred.",
            "Complex source structures and exact atom mapping require expert review.",
        ],
    }
    return r


def set_solvent(r, text, ev, *, ratio=None):
    r["solvent_conditions"]["composition"] = f(text, text, ev)
    if ratio:
        r["solvent_conditions"]["mixture_ratio"] = f(ratio, ratio, ev)


def add_step(r, sid, ev, pids, *, operation, temperature=None, time=None,
             atmosphere=None, pressure=None, stage="reaction"):
    s = st.new_step(sid, ev, stage=stage)
    s["participant_ids"] = list(pids)
    s["operation"] = f(operation, operation, ev)
    if temperature is not None:
        raw, norm = temperature
        s["temperature"] = f(raw, norm, ev)
    if time is not None:
        raw, norm = time
        s["time"] = f(raw, norm, ev)
    if atmosphere is not None:
        s["atmosphere"] = f(atmosphere, atmosphere, ev)
    if pressure is not None:
        raw, norm = pressure
        s["pressure"] = f(raw, norm, ev)
    r["steps"].append(s)
    return s


def add_yield(r, oid, product_mid, raw, value, ev, *, isolated=False,
              scope="single-step product yield"):
    o = st.new_outcome(oid, "yield", [product_mid], ev)
    o["product_role"] = f("main product", "main product", ev)
    o["value"] = f(raw, {"number": value, "unit": "%"}, ev)
    if isolated:
        o["yield_type"] = f("isolated yield", "isolated", ev)
        o["measurement_basis"] = f("isolated yield", "isolated", ev)
    o["scope"] = f(scope, scope, ev)
    r["outcomes"].append(o)
    return o


def add_selectivity(r, oid, product_mid, raw, normalized, ev, scope,
                    measurement=None):
    o = st.new_outcome(oid, "selectivity", [product_mid], ev)
    o["product_role"] = f("main product stereochemical/conformational result", "main product stereochemical/conformational result", ev)
    o["value"] = f(raw, normalized, ev)
    o["scope"] = f(scope, scope, ev)
    if measurement:
        o["measurement_basis"] = f(measurement, measurement, ev)
    r["outcomes"].append(o)
    return o


def main():
    bundle = json.loads((RUN / "submitted.json").read_text())
    molecules = []
    stereo_diagnostics = {}

    # Explicitly named simple reagents/solvents. Complex drawings remain unresolved.
    simple = {
        "m_dce": ("1,2-dichloroethane (DCE)", "ClCCCl", 2),
        "m_mecn": ("acetonitrile", "CC#N", 2),
        "m_dbu": ("DBU", "C1CCC2=NCCCN2CC1", 3),
        "m_et2o": ("diethyl ether", "CCOCC", 3),
        "m_nbs": ("N-bromosuccinimide (NBS)", "O=C1N(Br)C(=O)CC1", 3),
        "m_thf": ("tetrahydrofuran (THF)", "C1CCOC1", 3),
        "m_h2": ("hydrogen", "[H][H]", 3),
        "m_meoh": ("methanol", "CO", 3),
        "m_lialh4": ("lithium aluminium hydride", "[Li+].[AlH4-]", 3),
        "m_tfa": ("trifluoroacetic acid (TFA)", "O=C(O)C(F)(F)F", 3),
        "m_dcm": ("dichloromethane (DCM)", "ClCCl", 3),
        "m_tbuamine": ("tert-butylamine", "CC(C)(C)N", 3),
        "m_phsh": ("thiophenol (PhSH)", "Sc1ccccc1", 3),
        "m_ptolboronic": ("4-methylphenylboronic acid", "Cc1ccc(B(O)O)cc1", 3),
        "m_k2co3": ("potassium carbonate", "[K+].[K+].[O-]C([O-])=O", 3),
        "m_dioxane": ("1,4-dioxane", "O1CCOCC1", 3),
        "m_water": ("water", "O", 3),
        "m_diene": ("2,3-dimethylbuta-1,3-diene", "C=C(C)C(=C)C", 3),
        "m_toluene": ("toluene", "Cc1ccccc1", 3),
        "m_chcl3": ("chloroform", "ClC(Cl)Cl", 3),
        "m_tf2cch2": ("1,1-bis(triflyl)ethene (Tf2C=CH2)", "C=C(S(=O)(=O)C(F)(F)F)S(=O)(=O)C(F)(F)F", 2),
    }
    for mid, (name, smi, page) in simple.items():
        m, diag = smiles_molecule(mid, name, smi, page)
        molecules.append(m)
        stereo_diagnostics[mid] = diag

    for mid, label, name, page in [
        ("m_1", "1", "2-(2-fluoropyridinium-1-yl)-1,1-bis(triflyl)ethan-1-ide", 2),
        ("m_d1", "[D]-1", "deuterated Yanai reagent [D]-1", 2),
        ("m_pdc", "Pd/C", "palladium on carbon", 3),
        ("m_pdpph3_4", "Pd(PPh3)4", "tetrakis(triphenylphosphine)palladium(0)", 3),
        ("m_silica", "silica gel", "silica gel chromatography stationary phase", 3),
        ("m_hexanes", "hexanes", "hexanes (composition not specified)", 3),
        ("m_generic_A", "A", "generic anilide-derived allenol/alkenol A", 2),
        ("m_generic_product", "generic 3/7", "generic gem-bis(triflyl) cyclopenta[b]indoline product", 2),
        ("m_generic_2_or_6", "2,6", "generic allenol or alkenol substrate family 2/6", 3),
        ("m_generic_3_or_7", "3,7", "generic gem-bis(triflyl)indoline family 3/7", 3),
        ("m_int1", "INT-1", "proposed zwitterionic intermediate INT-1", 3),
        ("m_int2", "INT-2", "proposed bicyclic intermediate INT-2", 3),
        ("m_int3", "INT-3", "proposed oxonium intermediate INT-3", 3),
        ("m_int4", "INT-4", "proposed 2H-indol-1-ium intermediate INT-4", 3),
    ]:
        molecules.append(source_molecule(mid, label, name, page))

    # Scope substrates and products.
    for label in "abcdefghijklmno":
        molecules.append(source_molecule(f"m_2{label}", f"2{label}", f"anilide-derived allenol 2{label}", 2,
                                         stereo_kind="substrate_allene"))
        mix = None
        m = source_molecule(f"m_3{label}", f"3{label}", f"gem-bis(triflyl) tricyclic indoline 3{label}", 2,
                            stereo_kind="cis_product")
        if label == "d":
            m["stereochemistry"]["mixture"] = f("minor NMR signal set 8%, assigned by authors to rotamer rather than trans isomer", {"minor_rotamer_percent": 8}, E(2))
        if label == "o":
            m["stereochemistry"]["mixture"] = f("minor NMR signal set 12%, assigned by authors to rotamer rather than trans isomer", {"minor_rotamer_percent": 12}, E(2))
        molecules.append(m)
    molecules.append(source_molecule("m_d3a", "[D]-3a", "deuterated gem-bis(triflyl)indoline [D]-3a", 2,
                                     stereo_kind="cis_product"))
    for label in "abcdefg":
        molecules.append(source_molecule(f"m_6{label}", f"6{label}", f"anilide-derived allylic alcohol 6{label}", 3,
                                         stereo_kind="substrate_alkene"))
        molecules.append(source_molecule(f"m_7{label}", f"7{label}", f"gem-bis(triflyl) tricyclic indoline 7{label}", 3,
                                         stereo_kind="cis_product"))
    for number in range(8, 20):
        kind = "drawn_product" if number in (9, 12, 14, 15, 16, 17) else "unresolved"
        molecules.append(source_molecule(f"m_{number}", str(number), f"derivatization product {number}", 3,
                                         stereo_kind=kind))

    bundle["molecules"] = molecules
    reactions = []

    # Scheme 1c: current-work generic design. Prior-work panels 1a/1b remain source context.
    r = base_reaction("scheme1c_current_design", 2, ["other"], "Scheme 1c, current work",
                      "Current design for cyclopenta[b]indolines",
                      "bis(triflyl)ethylation–azacyclization–carbocyclization cascade",
                      record_type="generic_scheme")
    ev = E(2)
    r["participants"] = [
        participant("p_A", "m_generic_A", "reactant", ev),
        participant("p_tf2cch2", "m_tf2cch2", "reagent", ev),
        participant("p_generic_product", "m_generic_product", "product", ev, product_role="generic target product"),
    ]
    manual_classification(r, ["C-C", "C-N", "C-C"], ev,
                          "Scheme 1c explicitly depicts sequential C-C/C-N/C-C bond formation.")
    reactions.append(r)

    # Initial unexpected result in MeCN at room temperature (distinct from optimized 62% experiment).
    r = base_reaction("preliminary_2a_to_3a_35pct", 2, ["optimization"],
                      "page 2 prose below Scheme 2", "Reaction of carbamate-derived substrate 2a",
                      "cascade formation of tricyclic indoline 3a")
    ev = E(2)
    r["participants"] = [
        participant("p_2a", "m_2a", "reactant", ev), participant("p_1", "m_1", "reagent", ev),
        participant("p_mecn", "m_mecn", "solvent", ev), participant("p_3a", "m_3a", "product", ev, product_role="major product"),
    ]
    set_solvent(r, "acetonitrile", ev)
    add_step(r, "reaction", ev, [p["participant_id"] for p in r["participants"]],
             operation="2a reacted with Yanai reagent 1 in acetonitrile",
             temperature=("room temperature", "r.t."))
    add_yield(r, "yield_3a", "m_3a", "35% yield", 35, ev)
    r["notes"].append(f("3a was the major product along with several unidentified compounds.",
                        "3a was the major product along with several unidentified compounds.", ev))
    r["procedure"]["original_text"] = f(
        "reaction of aniline-derived allenol 2a with reagent 1 in acetonitrile at room temperature gave 3a as major product (35%) along with unidentified compounds",
        "reaction of aniline-derived allenol 2a with reagent 1 in acetonitrile at room temperature gave 3a as major product (35%) along with unidentified compounds", ev)
    manual_classification(r, ["C-C", "C-N", "C-C"], ev,
                          "The paper summarizes the cascade as formation of two C-C bonds and one C-N bond.")
    reactions.append(r)

    # Figure 2 / Scheme 2 scope. 3a is recorded once with both locations on the same page.
    fig2 = {
        "a": (62, 2.0), "b": (61, 1.0), "c": (66, 1.0), "d": (55, 1.0),
        "e": (61, 2.0), "f": (47, 1.0), "g": (57, 1.0), "h": (77, 1.0),
        "i": (68, 0.5), "j": (54, 1.0), "k": (49, 1.0), "l": (56, 0.5),
        "m": (40, 1.0), "n": (65, 1.0), "o": (41, 1.0),
    }
    for label, (yield_value, hours) in fig2.items():
        categories = ["substrate_scope"] + (["optimization"] if label == "a" else [])
        location = "Scheme 2 and Figure 2 (3a)" if label == "a" else f"Figure 2, product 3{label}"
        r = base_reaction(f"scope_2{label}_to_3{label}", 2, categories, location,
                          "Synthesis of bis(triflyl)-containing tricyclic indolines 3",
                          "metal-free cascade cyclization")
        ev = E(2)
        r["participants"] = [
            participant(f"p_2{label}", f"m_2{label}", "reactant", ev),
            participant("p_1", "m_1", "reagent", ev), participant("p_dce", "m_dce", "solvent", ev),
            participant(f"p_3{label}", f"m_3{label}", "product", ev, product_role="main product"),
        ]
        set_solvent(r, "1,2-dichloroethane (DCE)", ev)
        r["special_conditions"]["reactor"] = f("in a sealed tube", "sealed tube", ev)
        add_step(r, "reaction", ev, [p["participant_id"] for p in r["participants"]],
                 operation="heat substrate and reagent 1 in DCE in a sealed tube",
                 temperature=("130 °C", {"number": 130, "unit": "°C"}),
                 time=(("30 min", {"number": 30, "unit": "min"}) if hours == 0.5 else
                       (f"{hours:g} h", {"number": hours, "unit": "h"})))
        add_yield(r, f"yield_3{label}", f"m_3{label}", f"{yield_value}%", yield_value, ev, isolated=True)
        if label == "a":
            add_selectivity(r, "selectivity_3a", "m_3a", "single cis isomer",
                            "single cis isomer", ev,
                            "relative ring-junction stereochemistry; cis major/only isomer; absolute configuration unassigned")
            r["notes"].append(f("Use of bench-grade solvents did not affect the reaction outcome.",
                                "Use of bench-grade solvents did not affect the reaction outcome.", ev))
        if label in ("e", "f", "g"):
            add_selectivity(r, f"selectivity_3{label}", f"m_3{label}", "total selectivity",
                            "total selectivity", ev,
                            "author's qualitative product selectivity; absolute configurations unassigned")
        if label == "h":
            add_selectivity(r, "selectivity_3h", "m_3h", "single isomer",
                            "single isomer", ev,
                            "single product isomer; absolute configuration unassigned")
        if label in ("d", "o"):
            pct = 8 if label == "d" else 12
            add_selectivity(r, f"rotamer_3{label}", f"m_3{label}",
                            f"minor signals {pct}% as estimated by 1H NMR",
                            {"minor_rotamer_percent": pct}, ev,
                            "minor NMR signal set attributed by authors to rotamer, not trans stereoisomer",
                            measurement="1H NMR estimate")
        r["procedure"]["original_text"] = f(
            "DCE, 130 °C, in a sealed tube", "DCE, 130 °C, in a sealed tube", ev)
        manual_classification(r, ["C-C", "C-N", "C-C"], ev,
                              "The paper summarizes the observed cascade as formation of two C-C bonds and one C-N bond.")
        reactions.append(r)

    r = base_reaction("scope_2a_d1_to_d3a", 2, ["substrate_scope"],
                      "Figure 2, [D]-3a", "Synthesis of deuterated gem-bis(triflyl)indoline [D]-3a",
                      "deuterated metal-free cascade cyclization")
    ev = E(2)
    r["participants"] = [
        participant("p_2a", "m_2a", "reactant", ev), participant("p_d1", "m_d1", "reagent", ev),
        participant("p_dce", "m_dce", "solvent", ev), participant("p_d3a", "m_d3a", "product", ev, product_role="main deuterated product"),
    ]
    set_solvent(r, "1,2-dichloroethane (DCE)", ev)
    r["special_conditions"]["reactor"] = f("in a sealed tube", "sealed tube", ev)
    add_step(r, "reaction", ev, [p["participant_id"] for p in r["participants"]],
             operation="heat 2a with [D]-1 in DCE in a sealed tube",
             temperature=("130 °C", {"number": 130, "unit": "°C"}),
             time=("2 h", {"number": 2, "unit": "h"}))
    add_yield(r, "yield_d3a", "m_d3a", "59%", 59, ev, isolated=True)
    manual_classification(r, ["C-C", "C-N", "C-C"], ev,
                          "The deuterated example follows the same source-depicted C-C/C-N/C-C cascade.")
    reactions.append(r)

    # Figure 3 allylic-alcohol scope.
    fig3 = {"a": (64, 2.0), "b": (63, 1.0), "c": (34, 3.0), "d": (86, 2.0),
            "e": (52, 1.0), "f": (49, 2.5), "g": (67, 2.0)}
    for label, (yield_value, hours) in fig3.items():
        r = base_reaction(f"scope_6{label}_to_7{label}", 3, ["substrate_scope"],
                          f"Figure 3, product 7{label}",
                          "Synthesis of bis(triflyl)-containing tricyclic indolines 7",
                          "metal-free cascade cyclization of allylic alcohol")
        ev = E(3)
        r["participants"] = [
            participant(f"p_6{label}", f"m_6{label}", "reactant", ev), participant("p_1", "m_1", "reagent", ev),
            participant("p_dce", "m_dce", "solvent", ev), participant(f"p_7{label}", f"m_7{label}", "product", ev, product_role="main product"),
        ]
        set_solvent(r, "1,2-dichloroethane (DCE)", ev)
        r["special_conditions"]["reactor"] = f("in a sealed tube", "sealed tube", ev)
        add_step(r, "reaction", ev, [p["participant_id"] for p in r["participants"]],
                 operation="heat substrate and reagent 1 in DCE in a sealed tube",
                 temperature=("130 °C", {"number": 130, "unit": "°C"}),
                 time=(f"{hours:g} h", {"number": hours, "unit": "h"}))
        add_yield(r, f"yield_7{label}", f"m_7{label}", f"{yield_value}%", yield_value, ev, isolated=True)
        r["procedure"]["original_text"] = f("DCE, 130 °C, in a sealed tube", "DCE, 130 °C, in a sealed tube", ev)
        manual_classification(r, ["C-C", "C-N", "C-C"], ev,
                              "The 7-series follows the source-depicted C-C/C-N/C-C cascade.")
        reactions.append(r)

    # Scheme 3 derivatizations. Each printed arrow/starting-material alternative is independent.
    deriv_specs = []

    def deriv(rid, source, product, participants_spec, solvent, steps, yield_value,
              labels, transformation, *, notes=(), selectivity=None, mixture_ratio=None,
              reactor=None, purification=()):
        ev = E(3)
        r = base_reaction(rid, 3, ["derivatization"], f"Scheme 3, {rid}",
                          "Synthetic transformation of gem-bis(triflyl)indoline",
                          transformation)
        plist = [participant(f"p_{source}", f"m_{source}", "reactant", ev)]
        for pid, mid, role, kwargs in participants_spec:
            plist.append(participant(pid, mid, role, ev, **kwargs))
        plist.append(participant(f"p_{product}", f"m_{product}", "product", ev, product_role="main product"))
        r["participants"] = plist
        if solvent:
            set_solvent(r, solvent, ev, ratio=mixture_ratio)
        if reactor:
            r["special_conditions"]["reactor"] = f(reactor, reactor, ev)
        all_pids = [p["participant_id"] for p in plist if p["stage"] == "reaction"]
        for spec in steps:
            add_step(r, spec["sid"], ev, spec.get("pids", all_pids), operation=spec["operation"],
                     temperature=spec.get("temperature"), time=spec.get("time"),
                     atmosphere=spec.get("atmosphere"), pressure=spec.get("pressure"),
                     stage=spec.get("stage", "reaction"))
        for spec in purification:
            add_step(r, spec["sid"], ev, spec["pids"], operation=spec["operation"], stage=spec["stage"])
        add_yield(r, f"yield_{product}", f"m_{product}", f"{yield_value}%", yield_value, ev)
        if selectivity:
            add_selectivity(r, f"selectivity_{product}", f"m_{product}", selectivity[0], selectivity[1], ev, selectivity[2], selectivity[3] if len(selectivity) > 3 else None)
        for text in notes:
            r["notes"].append(f(text, text, ev))
        r["procedure"]["original_text"] = f(
            "; ".join(spec["operation"] for spec in steps),
            "; ".join(spec["operation"] for spec in steps), ev)
        manual_classification(r, labels, ev, f"Partial net-change classes derived from the inspected Scheme 3 arrow for {source} to {product}.")
        reactions.append(r)

    deriv("scheme3a_3a_to_8", "3a", "8",
          [("p_dbu", "m_dbu", "reagent", {}), ("p_et2o", "m_et2o", "solvent", {})],
          "diethyl ether", [{"sid": "reaction", "operation": "treat 3a with DBU in Et2O", "temperature": ("0 °C", {"number": 0, "unit": "°C"}), "time": ("5 min", {"number": 5, "unit": "min"})}],
          91, ["DB-C-C"], "base-induced elimination/ring opening to dienyl triflone")
    deriv("scheme3a_3a_to_9", "3a", "9",
          [("p_nbs", "m_nbs", "reagent", {}), ("p_thf", "m_thf", "solvent", {})],
          "THF", [{"sid": "reaction", "operation": "treat 3a with NBS in THF", "temperature": ("40 °C", {"number": 40, "unit": "°C"}), "time": ("12 h", {"number": 12, "unit": "h"})}],
          97, ["C-Br"], "bromination")
    deriv("scheme3b_3a_to_10", "3a", "10",
          [("p_h2", "m_h2", "reagent", {}), ("p_pdc", "m_pdc", "catalyst", {}), ("p_meoh", "m_meoh", "solvent", {}),
           ("p_hexanes", "m_hexanes", "workup_agent", {"stage": "workup"}), ("p_chcl3", "m_chcl3", "workup_agent", {"stage": "workup"})],
          "methanol", [{"sid": "reaction", "operation": "hydrogenate 3a over Pd/C in MeOH", "temperature": ("rt", "r.t."), "time": ("30 min", {"number": 30, "unit": "min"}), "atmosphere": "H2", "pressure": ("1 atm", "1 atm")}],
          75, ["C-H"], "hydrogenative ring opening",
          purification=[{"sid": "workup", "stage": "workup", "pids": ["p_hexanes", "p_chcl3"], "operation": "washing-up with hexanes:CHCl3 (9:1)"}])
    deriv("scheme3b_3a_to_11", "3a", "11",
          [("p_lialh4", "m_lialh4", "reagent", {}), ("p_et2o", "m_et2o", "solvent", {}),
           ("p_silica", "m_silica", "purification_agent", {"stage": "purification"})],
          "diethyl ether", [
              {"sid": "addition_0c", "operation": "add/treat with LiAlH4 in Et2O at 0 °C", "temperature": ("0 °C", {"number": 0, "unit": "°C"})},
              {"sid": "warm_to_rt", "operation": "warm from 0 °C to room temperature; printed 30 min is retained for this stage", "temperature": ("rt", "r.t."), "time": ("30 min", {"number": 30, "unit": "min"})},
          ], 63, ["C-H"], "LiAlH4 reduction/ring opening",
          notes=["The source prints 0 °C→rt, 30 min; separate addition and warming durations are not reported."],
          purification=[{"sid": "purification", "stage": "purification", "pids": ["p_silica"], "operation": "silica gel chromatography"}])
    for src, y in [("3c", 87), ("3d", 79)]:
        deriv(f"scheme3c_{src}_to_12", src, "12",
              [("p_tfa", "m_tfa", "reagent", {}), ("p_dcm", "m_dcm", "solvent", {})],
              "DCM", [{"sid": "reaction", "operation": f"treat {src} with TFA in DCM", "temperature": ("40 °C", {"number": 40, "unit": "°C"}), "time": ("12 h", {"number": 12, "unit": "h"})}],
              y, ["N-H"], "acidic N-deprotection")
    deriv("scheme3c_3b_to_13", "3b", "13",
          [("p_h2", "m_h2", "reagent", {}), ("p_pdc", "m_pdc", "catalyst", {}), ("p_meoh", "m_meoh", "solvent", {}),
           ("p_silica", "m_silica", "purification_agent", {"stage": "purification"})],
          "methanol", [{"sid": "reaction", "operation": "hydrogenate 3b over Pd/C in MeOH", "temperature": ("55 °C", {"number": 55, "unit": "°C"}), "time": ("1 h", {"number": 1, "unit": "h"}), "atmosphere": "H2", "pressure": ("1 atm", "1 atm")}],
          55, ["C-H", "N-H"], "hydrogenative reduction/deprotection and ring opening",
          purification=[{"sid": "purification", "stage": "purification", "pids": ["p_silica"], "operation": "silica gel chromatography"}])
    deriv("scheme3d_3a_to_14", "3a", "14",
          [("p_tbuamine", "m_tbuamine", "reagent", {}), ("p_et2o", "m_et2o", "solvent", {})],
          "diethyl ether", [
              {"sid": "addition_0c", "operation": "combine 3a and tert-butylamine in Et2O at 0 °C", "temperature": ("0 °C", {"number": 0, "unit": "°C"})},
              {"sid": "warm_to_rt", "operation": "warm from 0 °C to room temperature; printed 1.5 h retained for this stage", "temperature": ("rt", "r.t."), "time": ("1.5 h", {"number": 1.5, "unit": "h"})},
          ], 84, ["C-N", "DB-C-C"], "amine addition following elimination/ring opening",
          notes=["The source prints 0 °C→rt, 1.5 h; separate addition and warming durations are not reported."],
          selectivity=("dr = 80:20", {"major": 80, "minor": 20, "metric": "dr"}, "diastereomeric ratio; major/minor absolute configurations unassigned"))
    deriv("scheme3d_3a_to_15", "3a", "15",
          [("p_dbu", "m_dbu", "reagent", {}), ("p_phsh", "m_phsh", "reagent", {}), ("p_et2o", "m_et2o", "solvent", {})],
          "diethyl ether", [
              {"sid": "elimination", "operation": "treat 3a with DBU in Et2O", "temperature": ("0 °C", {"number": 0, "unit": "°C"}), "time": ("5 min", {"number": 5, "unit": "min"})},
              {"sid": "thiol_addition", "operation": "add PhSH in Et2O", "temperature": ("0 °C", {"number": 0, "unit": "°C"}), "time": ("30 min", {"number": 30, "unit": "min"})},
          ], 63, ["C-S", "DB-C-C"], "thiol addition following elimination/ring opening")
    deriv("scheme3e_3j_to_16", "3j", "16",
          [("p_ptolboronic", "m_ptolboronic", "reagent", {}),
           ("p_pdpph3_4", "m_pdpph3_4", "catalyst", {"mol_percent": 5, "reference": "m_3j"}),
           ("p_k2co3", "m_k2co3", "reagent", {}), ("p_dioxane", "m_dioxane", "solvent", {}), ("p_water", "m_water", "solvent", {})],
          "1,4-dioxane:H2O", [{"sid": "reaction", "operation": "Suzuki–Miyaura coupling with ArB(OH)2, Pd(PPh3)4 and K2CO3", "temperature": ("60 °C", {"number": 60, "unit": "°C"}), "time": ("1 h", {"number": 1, "unit": "h"})}],
          77, ["C-C", "C-H"], "Suzuki–Miyaura cross-coupling accompanied by detriflylation",
          mixture_ratio="1,4-dioxane:H2O = 2:1",
          notes=["Ar = 4-MeC6H4; the prose states that detriflylation also occurred."])
    deriv("scheme3e_3a_to_17", "3a", "17",
          [("p_dbu", "m_dbu", "reagent", {}), ("p_et2o", "m_et2o", "solvent", {}),
           ("p_diene", "m_diene", "reagent", {}), ("p_toluene", "m_toluene", "solvent", {})],
          "Et2O then toluene", [
              {"sid": "elimination", "operation": "treat 3a with DBU in Et2O to form dienyl triflone 8 in situ", "temperature": ("0 °C", {"number": 0, "unit": "°C"}), "time": ("5 min", {"number": 5, "unit": "min"})},
              {"sid": "diels_alder", "operation": "react with 2,3-dimethylbuta-1,3-diene in toluene in a sealed tube", "temperature": ("120 °C", {"number": 120, "unit": "°C"}), "time": ("48 h", {"number": 48, "unit": "h"})},
          ], 45, ["DB-C-C", "C-C", "C-C"], "elimination followed by Diels–Alder cycloaddition",
          reactor="sealed tube",
          notes=["The displayed 45% is retained as the overall 3a-to-17 sequence yield, not assigned to the second step alone."])
    deriv("scheme3f_7a_to_18", "7a", "18",
          [("p_k2co3", "m_k2co3", "reagent", {}), ("p_dioxane", "m_dioxane", "solvent", {}), ("p_water", "m_water", "solvent", {})],
          "1,4-dioxane:water", [{"sid": "reaction", "operation": "treat 7a with K2CO3 in 1,4-dioxane:water", "temperature": ("70 °C", {"number": 70, "unit": "°C"}), "time": ("3 h", {"number": 3, "unit": "h"})}],
          86, ["DB-C-C"], "elimination to alkenyl triflone", mixture_ratio="1,4-dioxane:water = 2:1")
    deriv("scheme3f_7a_to_19", "7a", "19",
          [("p_lialh4", "m_lialh4", "reagent", {}), ("p_et2o", "m_et2o", "solvent", {}),
           ("p_silica", "m_silica", "purification_agent", {"stage": "purification"})],
          "diethyl ether", [
              {"sid": "addition_0c", "operation": "add/treat 7a with LiAlH4 in Et2O at 0 °C", "temperature": ("0 °C", {"number": 0, "unit": "°C"})},
              {"sid": "warm_to_rt", "operation": "warm from 0 °C to room temperature; printed 40 min retained for this stage", "temperature": ("rt", "r.t."), "time": ("40 min", {"number": 40, "unit": "min"})},
          ], 80, ["C-H"], "LiAlH4 reduction/ring opening",
          notes=["The source prints 0 °C→rt, 40 min; separate addition and warming durations are not reported."],
          purification=[{"sid": "purification", "stage": "purification", "pids": ["p_silica"], "operation": "silica gel chromatography"}])

    # Scheme 4 is a proposed mechanism, not a verified experiment.
    r = base_reaction("scheme4_proposed_mechanism", 3, ["other"], "Scheme 4",
                      "Tentative pathway for the formation of gem-bis(triflyl)indolines 3 and 7",
                      "proposed electrophilic addition/azacyclization/dehydration/carbocyclization sequence",
                      record_type="proposed_mechanism")
    ev = E(3) + E(4)
    r["evidence_ids"] = ev
    r["participants"] = [
        participant("p_2_or_6", "m_generic_2_or_6", "reactant", ev), participant("p_tf2cch2", "m_tf2cch2", "reagent", ev),
        participant("p_int1", "m_int1", "reagent", ev), participant("p_int2", "m_int2", "reagent", ev),
        participant("p_int3", "m_int3", "reagent", ev), participant("p_int4", "m_int4", "reagent", ev),
        participant("p_3_or_7", "m_generic_3_or_7", "product", ev, product_role="proposed final product family"),
    ]
    for sid, operation in [
        ("proposed_1", "electrophilic attack of Tf2C=CH2 at the beta carbon to form INT-1"),
        ("proposed_2", "amide-nitrogen cyclization to INT-2"),
        ("proposed_3", "proton transfer and protonation to INT-3"),
        ("proposed_4", "dehydration to INT-4"),
        ("proposed_5", "intramolecular ionic carbocyclization to 3/7"),
    ]:
        add_step(r, sid, ev, [p["participant_id"] for p in r["participants"]], operation=operation)
    r["notes"].append(f(
        "The paper labels Scheme 4 a tentative pathway; DFT gives a 23.1 kcal mol−1 first-step barrier and a 1.5 kcal mol−1 following C−N-forming barrier.",
        "The paper labels Scheme 4 a tentative pathway; DFT gives a 23.1 kcal mol−1 first-step barrier and a 1.5 kcal mol−1 following C−N-forming barrier.", ev))
    manual_classification(r, ["C-C", "C-N", "C-C"], ev,
                          "The proposed final cascade contains two C-C formations and one C-N formation, but this record is explicitly mechanistic rather than experimental.")
    r["review"]["unresolved"].append("Proposed intermediates and elementary steps are author hypotheses, not independently verified experiments.")
    reactions.append(r)

    bundle["reactions"] = reactions
    bundle["source_context"] = [
        {"context_id": f"{DOC}:article_metadata", "document_id": DOC, "kind": "article_metadata",
         "evidence_ids": E(1), "summary": "Article title, authors, journal citation, abstract, received/published dates and DOI."},
        {"context_id": f"{DOC}:graphical_abstract", "document_id": DOC, "kind": "graphical_summary",
         "evidence_ids": E(1), "summary": "Graphical abstract summarizes 23 examples and up to 86% yield; it is not counted as an additional experiment."},
        {"context_id": f"{DOC}:figure1", "document_id": DOC, "kind": "background",
         "evidence_ids": E(1), "summary": "Natural-product/bioactive indoline structures are background examples, not experiments in this paper."},
        {"context_id": f"{DOC}:scheme1_prior", "document_id": DOC, "kind": "prior_literature",
         "evidence_ids": E(2), "summary": "Scheme 1a-b depicts cited prior electrophilic dearomatization and cycloaddition work; excluded from current-paper experiments."},
        {"context_id": f"{DOC}:scheme2_alternatives", "document_id": DOC, "kind": "unobserved_candidates",
         "evidence_ids": E(2), "summary": "Compounds 4 and 5 are depicted as products not obtained in place of major 3a; no synthetic records were fabricated for them."},
        {"context_id": f"{DOC}:xray", "document_id": DOC, "kind": "characterization",
         "evidence_ids": E(2) + E(3) + E(6), "summary": "X-ray structures/CCDC entries: 3a 2056816, 3e 2056817, 7d 2056818, 7g 2056819; characterization is linked context, not additional synthesis."},
        {"context_id": f"{DOC}:dft", "document_id": DOC, "kind": "computational_mechanistic_support",
         "evidence_ids": E(4) + E(6), "summary": "DFT discussion supports the tentative pathway and reports barriers; it is not a preparative experiment."},
        {"context_id": f"{DOC}:si_notice", "document_id": DOC, "kind": "supporting_information_notice",
         "evidence_ids": E(4), "summary": "The paper states SI exists, but no SI file was supplied or explicitly associated in this workspace."},
        {"context_id": f"{DOC}:bibliography", "document_id": DOC, "kind": "bibliography",
         "evidence_ids": E(4) + E(5) + E(6), "summary": "References and explanatory footnotes were inspected as context; cited syntheses are not current-paper experiments."},
    ]

    # Raw session artifacts and formal submitted bundle.
    (RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
    (RUN / "raw_model/stereo_diagnostics.json").write_text(json.dumps(stereo_diagnostics, ensure_ascii=False, indent=2) + "\n")
    inventory = {
        "document_id": DOC,
        "pages_inspected": [1, 2, 3, 4, 5, 6],
        "molecule_count": len(molecules),
        "reaction_count": len(reactions),
        "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in reactions),
        "generic_schemes": sum(r["record_type"] == "generic_scheme" for r in reactions),
        "proposed_mechanisms": sum(r["record_type"] == "proposed_mechanism" for r in reactions),
        "reaction_ids": [r["reaction_id"] for r in reactions],
        "note": "Inventory is a machine transcription aid, not source-reviewed chemical truth.",
    }
    (RUN / "raw_model/reaction_inventory.json").write_text(json.dumps(inventory, ensure_ascii=False, indent=2) + "\n")

    p2_ids = [r["reaction_id"] for r in reactions if any(e.startswith(f"{DOC}:page2:") for e in r["evidence_ids"])]
    p3_ids = [r["reaction_id"] for r in reactions if any(e.startswith(f"{DOC}:page3:") for e in r["evidence_ids"])]
    coverage = {
        "document_id": DOC,
        "semantic_accuracy_verified": False,
        "pages": [
            {"page": 1, "state": "attempted", "regions": [
                {"region": "title/abstract/graphical abstract", "finding": "Article metadata and a graphical summary of the current cascade; summary is not an additional experiment.", "evidence_ids": E(1)},
                {"region": "Figure 1 and introduction", "finding": "Bioactive/natural products and cited prior approaches are background, not current-paper experiments.", "evidence_ids": E(1)}],
             "record_ids": [], "unresolved": ["Graphical-abstract generic structures are summaries and were not used to invent compound-level structures or charges."]},
            {"page": 2, "state": "attempted", "regions": [
                {"region": "Scheme 1a-b", "finding": "Cited prior-work reaction concepts; retained as source context, not current-paper records.", "evidence_ids": E(2)},
                {"region": "Scheme 1c", "finding": "Current-work generic C-C/C-N/C-C cascade design; one generic_scheme record.", "evidence_ids": E(2)},
                {"region": "Scheme 2 and adjacent prose", "finding": "One preliminary 2a/1 experiment in MeCN at room temperature (3a, 35%) and the optimized DCE/130 °C sealed-tube experiment (3a, 62%, 2 h); compounds 4/5 are unobserved alternatives.", "evidence_ids": E(2)},
                {"region": "Figure 2", "finding": "Sixteen scope points: 3a-3o plus [D]-3a; 3a is merged with duplicate Scheme 2 evidence on the same page.", "evidence_ids": E(2)},
                {"region": "X-ray images/prose", "finding": "3a and 3e characterization and relative-stereochemistry context; not extra syntheses.", "evidence_ids": E(2)}],
             "record_ids": p2_ids, "unresolved": ["No SI was supplied; substrate/reagent charges, scale, workup and purification are unreported in the supplied main-text page.", "Complex 2/3 structures and absolute configurations remain unencoded; 3d/3o minor NMR signals are retained as author-assigned rotamers."]},
            {"page": 3, "state": "attempted", "regions": [
                {"region": "Figure 3", "finding": "Seven allylic-alcohol scope points 6a-g to 7a-g, one record per displayed product.", "evidence_ids": E(3)},
                {"region": "Scheme 3a-f", "finding": "Thirteen distinct derivatization records, including separate 3c/3d deprotections to 12 and the overall two-stage 3a-to-17 sequence.", "evidence_ids": E(3)},
                {"region": "Scheme 3 prose", "finding": "Author descriptions, product 17 stereoselectivity claim and synthetic-utility discussion mapped to records/notes; missing methods preserved.", "evidence_ids": E(3)},
                {"region": "Scheme 4", "finding": "Tentative proposed mechanism recorded separately, not as an experimental synthesis.", "evidence_ids": E(3)},
                {"region": "X-ray images", "finding": "7d/7g characterization context, not additional reactions.", "evidence_ids": E(3)}],
             "record_ids": p3_ids, "unresolved": ["Scheme 3 yield measurement bases and most charges/scales are not reported on the supplied page.", "Complex product connectivity/absolute stereochemistry and complete atom maps remain unresolved."]},
            {"page": 4, "state": "attempted", "regions": [
                {"region": "mechanistic discussion and conclusion", "finding": "Continuation of tentative Scheme 4 discussion and DFT barriers; no new preparative experiment.", "evidence_ids": E(4)},
                {"region": "Associated Content/Accession Codes", "finding": "SI and CCDC availability notices only; SI/CIF files were not supplied or associated.", "evidence_ids": E(4)},
                {"region": "author information/references", "finding": "Metadata and cited literature, not current-paper experiments.", "evidence_ids": E(4)}],
             "record_ids": [f"{DOC}:scheme4_proposed_mechanism"], "unresolved": ["DFT details said to be in unavailable SI were not reconstructed."]},
            {"page": 5, "state": "attempted", "regions": [
                {"region": "references 4-11", "finding": "Bibliography-only page; cited reactions are prior literature and were not imported as current experiments.", "evidence_ids": E(5)}],
             "record_ids": [], "unresolved": []},
            {"page": 6, "state": "attempted", "regions": [
                {"region": "references 12-16", "finding": "Bibliography and crystallographic footnote context; no new current-paper reaction arrow.", "evidence_ids": E(6)},
                {"region": "footnotes 17-20", "finding": "Rotamer/NOE, speculative reduction mechanism and DFT/mechanistic remarks; contextual support only, not independent synthetic experiments.", "evidence_ids": E(6)}],
             "record_ids": [f"{DOC}:scheme4_proposed_mechanism"], "unresolved": ["NOE and Figure S1 details are referenced but unavailable because SI was not supplied."]},
        ],
    }
    (RUN / "coverage/pages_0001_0006.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")

    print(json.dumps(inventory, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
