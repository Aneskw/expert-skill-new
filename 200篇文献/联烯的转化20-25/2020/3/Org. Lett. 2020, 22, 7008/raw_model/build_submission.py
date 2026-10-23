#!/usr/bin/env python3
"""Build the source-derived v3.1 submission for DOI 10.1021/acs.orglett.0c02558.

This generator is intentionally document-specific.  It preserves unresolved dense
structures and source conflicts instead of inventing complete structures or SI data.
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
DOC = "doc_6809ee7b65b00f0c"
E = {p: f"{DOC}:page{p}:image" for p in range(1, 6)}
T = {p: f"{DOC}:page{p}:text" for p in range(1, 6)}
TITLE = "Phosphine-Catalyzed Stereoselective Tandem Annulation Reaction for the Synthesis of Chromeno[4,3-b]pyrroles"
_AUTO = object()


def f(raw=None, normalized=_AUTO, ev=(), method="original_extraction", derivation=None,
      status=None, verification=None, alternatives=()):
    if normalized is _AUTO:
        normalized = raw
    return st.fact(raw, normalized, ev, method=method, derivation=derivation,
                   status=status, verification=verification, alternatives=alternatives)


def inferred(normalized, ev, derivation, raw=None):
    return f(raw, normalized, ev, method="inference", derivation=derivation)


def derived(normalized, ev, derivation, raw=None):
    return f(raw, normalized, ev, method="structure_derivation", derivation=derivation)


def conflict(alternatives):
    return f([a[0] for a in alternatives], None,
             sorted({eid for _, _, ids in alternatives for eid in ids}),
             status="conflict", alternatives=[f(raw, norm, ids) for raw, norm, ids in alternatives])


def q(number, unit, qualifier=None, raw_text=None):
    value = {"number": number, "unit": unit}
    if qualifier:
        value["qualifier"] = qualifier
    if raw_text:
        value["raw_text"] = raw_text
    return value


bundle = json.loads((RUN / "submitted.json").read_text())
bundle["schema_version"] = "3.1"
bundle["tool_version"] = "3.1"
bundle["molecules"] = []
bundle["reactions"] = []
bundle["source_context"] = []
bundle["documents"][0]["title"] = TITLE

molecules = {}


def base_molecule(mid, label, page, name=None):
    ev = [E[page]]
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    m["names"]["original"] = f(name or f"compound {label}", name or f"compound {label}", ev)
    m["names"]["english"] = f()
    m["names"]["abbreviations"] = f(label, [label], ev)
    m["names"]["aliases"] = f()
    m["cas_numbers"] = f()
    m["identity_note"] = f("Identity retained from the source label/drawing; CAS not reported.",
                             "Identity retained from source label/drawing; CAS unreported.", ev)
    m["external_properties"] = []
    molecules[mid] = m
    return m


def set_stereo(m, page, chirality, geometry, representation, atom=None, double=None,
               axial=None, assignment=None, mixture=None, derivation=None):
    ev = [E[page]]
    method = "structure_derivation" if derivation else "original_extraction"
    kwargs = {"method": method, "derivation": derivation}
    m["stereochemistry"]["atom_centers"] = f(atom, atom, ev, **kwargs) if atom is not None else f()
    m["stereochemistry"]["double_bonds"] = f(double, double, ev, **kwargs) if double is not None else f()
    m["stereochemistry"]["axial"] = f(axial, axial, ev, **kwargs) if axial is not None else f()
    m["stereochemistry"]["assignment"] = f(assignment, assignment, ev, **kwargs) if assignment is not None else f()
    m["stereochemistry"]["mixture"] = f(mixture, mixture, ev, **kwargs) if mixture is not None else f()
    m["stereochemistry"]["chirality_status"] = f(chirality, chirality, ev, **kwargs)
    m["stereochemistry"]["geometry_status"] = f(geometry, geometry, ev, **kwargs)
    m["stereochemistry"]["representation_status"] = f(representation, representation, ev, **kwargs)
    m["stereochemistry"]["isomeric_smiles"] = f()
    m["stereochemistry"]["identity_key"] = f()


def add_smiles(mid, label, page, smiles, name, chirality="achiral", geometry="not_applicable",
               representation="not_applicable", stereo_raw=None):
    m = base_molecule(mid, label, page, name)
    m["representation_type"] = "smiles"
    m["smiles"] = derived(smiles, [E[page]],
                           "Candidate SMILES derived from the source label/drawing; parsing does not verify source identity.",
                           raw=stereo_raw or f"source label/drawing: {label}")
    set_stereo(m, page, chirality, geometry, representation,
               assignment="Source stereochemical state retained; no unreported configuration was invented.",
               derivation="Status assessed from the displayed structure/name; candidate graph remains needs_review.")
    if representation == "encoded":
        m["stereochemistry"]["isomeric_smiles"] = derived(
            smiles, [E[page]], "Candidate isomeric SMILES encodes the source-stated geometry.", raw=stereo_raw or label)
    return m


def add_unresolved(mid, label, page, name=None, kind="ordinary"):
    m = base_molecule(mid, label, page, name)
    m["representation_type"] = "unresolved"
    m["smiles"] = f()
    if kind == "product":
        set_stereo(
            m, page, "partially_specified", "partially_specified", "not_encoded",
            atom="The product drawing depicts three contiguous stereogenic centers with wedge/dash relative stereochemistry.",
            double="An exocyclic alkene is depicted; the article discusses Z/E selectivity, but absolute atom-level encoding is unresolved.",
            assignment="Relative configuration is depicted; absolute configuration and major enantiomer identity are not assigned in the supplied main text.")
    elif kind == "allenoate":
        set_stereo(
            m, page, "unspecified", "not_applicable", "not_encoded",
            axial="Allene axis is present in the source drawing; axial configuration is not specified.",
            assignment="Potential central/axial stereochemistry is not assigned in the supplied source.")
    elif kind == "imine":
        set_stereo(m, page, "unspecified", "unspecified", "not_encoded",
                   double="Imine geometry is not assigned for this source-labeled scope substrate.",
                   assignment="No absolute or relative stereochemical assignment is supplied.")
    else:
        set_stereo(m, page, "unresolved", "unresolved", "unresolved",
                   assignment="Dense or abbreviated source structure retained without unsafe stereochemical encoding.")
    return m


def add_template(mid, label, page, fixed, variables):
    m = base_molecule(mid, label, page, f"generic compound {label}")
    m["representation_type"] = "generic_template"
    m["smiles"] = f()
    ev = [E[page]]
    m["generic_template"] = {
        "fixed_structure": f(fixed, {"description": fixed}, ev),
        "variables": [{
            "label": v,
            "attachment_points": inferred("source-depicted attachment point", ev, "Attachment follows the generic source drawing.", raw=v),
            "definition": f(v, "variable substituent", ev),
            "allowed_values": f(),
            "constraints": f(),
        } for v in variables],
        "instance_mapping": [],
        "reason": f("Generic R-group structure in the source.", "Generic R-group source structure cannot be a concrete SMILES.", ev),
    }
    set_stereo(m, page, "unresolved", "unresolved", "not_encoded",
               assignment="Generic variables prevent a unique stereochemical identity.")
    return m


# Generic source structures.
add_template("m_generic_1", "1", 1, "aldimine malonate bearing an ortho-phenol and variable aryl substituents", ["R1"])
add_template("m_generic_2", "2", 1, "delta-acetoxy allenoate with ester and variable R2/R3 substituents", ["R2", "R3"])
add_template("m_generic_3", "3", 1, "chromeno[4,3-b]pyrrole annulation product with variable R1/R2/R3 groups", ["R1", "R2", "R3"])

# Concrete model substrates and controls. Candidate structures are deliberately limited to simple, legible cases.
add_smiles("m_1a", "1a", 2, "CCOC(=O)C(N=Cc1ccccc1O)C(=O)OCC",
           "diethyl 2-[(2-hydroxybenzylidene)amino]malonate",
           geometry="unspecified", representation="partial")
m2a = add_smiles("m_2a", "2a", 2, "CCOC(=O)C=C=CC(OC(C)=O)c1ccccc1",
                  "ethyl delta-acetoxy phenyl-substituted allenoate",
                  chirality="unspecified", representation="partial")
m2a["stereochemistry"]["axial"] = f(
    "Allene axis depicted without an axial descriptor.",
    [{"axis": "C=C=C", "configuration": "unspecified"}], [E[2]])
m2a["stereochemistry"]["atom_centers"] = inferred(
    [{"site": "acetoxy/phenyl-bearing carbon", "configuration": "unspecified"}], [E[2]],
    "The displayed carbon has four different substituent classes; configuration is not assigned.", raw="source drawing 2a")
add_smiles("m_4", "4", 3, "CCOC(=O)C/N=C/c1ccccc1O",
           "ethyl (E)-2-[(2-hydroxybenzylidene)amino]acetate",
           geometry="specified", representation="encoded", stereo_raw="ethyl (E)-2-[(2-hydroxybenzylidene)amino]acetate")
add_smiles("m_5", "5", 3, "CCOC(=O)C(/N=C/c1ccccc1)C(=O)OCC",
           "diethyl (E)-2-(benzylideneamino)malonate",
           geometry="specified", representation="encoded", stereo_raw="diethyl (E)-2-(benzylideneamino)malonate")

# Scope substrate labels whose exact complete SMILES/stereochemistry are not safely reconstructed.
for label in [f"1{x}" for x in "bcdefghijk"]:
    add_unresolved(f"m_{label}", label, 2, kind="imine")
for label in [f"2{x}" for x in "bcdefghijklmn"]:
    add_unresolved(f"m_{label}", label, 2, kind="allenoate")

# Product labels across racemic and asymmetric scopes.
product_labels = {
    "3aa", "3ab", "3ac", "3ad", "3ae", "3af", "3ag", "3ah", "3ai", "3aj", "3ak", "3al", "3am",
    "3ba", "3ca", "3da", "3ea", "3fa", "3ga", "3ha", "3ia", "3ja", "3ka", "3bm", "3cm", "3dm"
}
for label in sorted(product_labels):
    add_unresolved(f"m_{label}", label, 2 if label not in {"3am", "3bm", "3cm", "3dm"} else 3, kind="product")

# Scheme 7 drawn starting structure with R=H is separated because prose calls it 3aa while the drawing bears CO2CHPh2.
add_unresolved("m_scheme7_start_RH", "Scheme 7 starting material (R=H, CO2CHPh2 drawing)", 3, kind="product")
add_unresolved("m_6", "6", 3, kind="product")
add_unresolved("m_7", "7", 3, kind="product")

# Reagents, catalysts and solvents. These candidates are graph diagnostics, not external identity certification.
simple = [
    ("m_PPh3", "PPh3", 2, "P(c1ccccc1)(c1ccccc1)c1ccccc1", "triphenylphosphine"),
    ("m_PEtPh2", "PEtPh2", 2, "CCP(c1ccccc1)c1ccccc1", "ethyldiphenylphosphine"),
    ("m_PEt2Ph", "PEt2Ph", 2, "CCP(CC)c1ccccc1", "diethylphenylphosphine"),
    ("m_PBu3", "PBu3", 2, "CCCCP(CCCC)CCCC", "tributylphosphine"),
    ("m_K2CO3", "K2CO3", 2, "[K+].[K+].[O-]C([O-])=O", "potassium carbonate"),
    ("m_PhMe", "PhMe", 2, "Cc1ccccc1", "toluene"),
    ("m_THF", "THF", 2, "C1CCOC1", "tetrahydrofuran"),
    ("m_CHCl3", "CHCl3", 2, "ClC(Cl)Cl", "chloroform"),
    ("m_CH2Cl2", "CH2Cl2", 2, "ClCCl", "dichloromethane"),
    ("m_CH3CN", "CH3CN", 2, "CC#N", "acetonitrile"),
    ("m_Et2O", "Et2O", 2, "CCOCC", "diethyl ether"),
    ("m_acetone", "acetone", 2, "CC(=O)C", "acetone"),
    ("m_NaBH4", "NaBH4", 3, "[BH4-].[Na+]", "sodium borohydride"),
    ("m_PhB_OH2", "PhB(OH)2", 3, "OB(O)c1ccccc1", "phenylboronic acid"),
    ("m_Pd_OAc2", "Pd(OAc)2", 3, "CC(=O)[O-].CC(=O)[O-].[Pd+2]", "palladium(II) acetate"),
]
for args in simple:
    add_smiles(*args)
add_unresolved("m_C10", "C10", 3, "chiral phosphine catalyst C10")
add_unresolved("m_XPhos", "XPhos", 3, "XPhos ligand")


def participant(pid, mid, role, ev, stage="reaction", amount=None, mass=None, volume=None,
                equiv=None, molpct=None, reference=None, product_role=None):
    p = st.new_participant(pid, mid, role, ev, stage)
    if product_role:
        p["product_role"] = f(product_role, product_role, ev)
    for key, value in (("amount", amount), ("mass", mass), ("volume", volume),
                       ("equivalents", equiv), ("mol_percent", molpct)):
        if value is not None:
            p["charges"][key] = f(value.get("raw_text", f"{value['number']} {value['unit']}"), value, ev)
    if reference:
        p["equivalence_reference"] = inferred(reference, ev,
            "Equivalence/loading denominator is the displayed 1.0 equiv limiting aldimine substrate.", raw=reference)
    return p


def step(step_id, ev, pids, operation=None, temperature=None, time=None):
    s = st.new_step(step_id, ev)
    s["participant_ids"] = list(pids)
    if operation is not None:
        s["operation"] = f(operation, operation, ev)
    if temperature is not None:
        s["temperature"] = f(temperature if isinstance(temperature, str) else temperature.get("raw_text"), temperature, ev)
    if time is not None:
        s["time"] = f(time if isinstance(time, str) else time.get("raw_text"), time, ev)
    return s


def outcome(oid, kind, target, ev, value=None, basis=None, yield_type=None, scope=None, value_fact=None):
    o = st.new_outcome(oid, kind, target, ev)
    o["product_role"] = f("main product", "main_product", ev) if target else f()
    o["value"] = value_fact if value_fact is not None else f(value if isinstance(value, str) else value.get("raw_text") if isinstance(value, dict) else value, value, ev)
    if basis is not None:
        o["measurement_basis"] = f(basis, basis, ev)
    if yield_type is not None:
        o["yield_type"] = f(yield_type, yield_type, ev)
    if scope is not None:
        o["scope"] = f(scope, scope, ev)
    return o


def classification(ev, mode="annulation"):
    if mode == "none":
        return {
            "standard_id": st.STANDARD, "mapped_reaction_smiles": f(), "bond_changes": [], "labels": f([], [], ev),
            "coverage": inferred("not_applicable", ev, "No target-product transformation was observed."),
            "limitations": ["No observed product; intended transformation is not classified."], "model_votes": []}
    if mode == "unresolved":
        return {
            "standard_id": st.STANDARD, "mapped_reaction_smiles": f(), "bond_changes": [], "labels": f(),
            "coverage": inferred("unresolved", ev, "Generic prior-work drawing was not atom-mapped."),
            "limitations": ["Generic literature example; exact bond changes were not asserted."], "model_votes": []}
    if mode == "annulation":
        changes = [
            {"change_id": "b1", "atoms": [{"element": "C"}, {"element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b2", "atoms": [{"element": "C"}, {"element": "N"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b3", "atoms": [{"element": "C"}, {"element": "O"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b4", "atoms": [{"element": "N"}, {"element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
            {"change_id": "b5", "atoms": [{"element": "O"}, {"element": "H"}], "before_order": 1, "after_order": 0, "change_type": "hydrogen_loss", "hydrogen_count": 1},
        ]
        labels = [
            {"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1"]},
            {"category_number": 5, "code": "C-N", "supporting_change_ids": ["b2"]},
            {"category_number": 3, "code": "C-O", "supporting_change_ids": ["b3"]},
            {"category_number": 29, "code": "N-H", "supporting_change_ids": ["b4"]},
        ]
        deriv = "Partial net changes derived from the highlighted product-forming bonds and N-H/O-H difference in the source drawings; exact atom mapping and allene/acetate accounting remain unresolved."
    elif mode == "reduction":
        changes = [
            {"change_id": "b1", "atoms": [{"element": "C"}, {"element": "O"}], "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
            {"change_id": "b2", "atoms": [{"element": "C"}, {"element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 2},
            {"change_id": "b3", "atoms": [{"element": "O"}, {"element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
            {"change_id": "b4", "atoms": [{"element": "C"}, {"element": "O"}], "before_order": 1, "after_order": 0, "change_type": "broken"},
        ]
        labels = [
            {"category_number": 1, "code": "C-H", "supporting_change_ids": ["b2"]},
            {"category_number": 24, "code": "O-H", "supporting_change_ids": ["b3"]},
        ]
        deriv = "Net ester-to-alcohol changes derived from Scheme 7; hydrogen source and leaving alcohol are not fully mapped."
    else:  # Suzuki coupling
        changes = [
            {"change_id": "b1", "atoms": [{"element": "C"}, {"element": "Cl"}], "before_order": 1, "after_order": 0, "change_type": "broken"},
            {"change_id": "b2", "atoms": [{"element": "C"}, {"element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        ]
        labels = [{"category_number": 2, "code": "C-C", "supporting_change_ids": ["b2"]}]
        deriv = "Net aryl C-Cl replacement by aryl C-C coupling derived from Scheme 7; boron-containing coproducts are omitted."
    return {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": f(),
        "bond_changes": [derived(ch, ev, deriv, raw="source reaction drawing") for ch in changes],
        "labels": derived(labels, ev, deriv, raw="source reaction drawing"),
        "coverage": inferred("partial", ev, deriv, raw="partial source-derived classification"),
        "limitations": [deriv, "No complete atom-mapped reaction SMILES was asserted."],
        "model_votes": [],
    }


def reaction(rid, page, record_type, categories, location, name, transformation, participants=(),
             steps=(), outcomes=(), scheme_ids=(), class_mode="unresolved", evidence=None, unresolved=None, notes=()):
    ev = list(evidence or [E[page]])
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, ev) if categories is not None else f()
    r["location"] = f(location, {"page": page, "region": location}, ev)
    r["general_scheme_ids"] = [f"{DOC}:{x}" for x in scheme_ids]
    r["reaction_name"] = f(name, name, ev) if name else f()
    r["transformation_type"] = inferred(transformation, ev, "Normalized transformation type from the source reaction title/drawing.", raw=name or location) if transformation else f()
    r["participants"] = list(participants)
    scale_participant = next((p for p in r["participants"]
                              if p["role"] == "reactant" and p["charges"]["amount"]["normalized_value"] is not None), None)
    if scale_participant:
        r["reaction_scale"] = inferred(
            {"reference_molecule_id": scale_participant["molecule_id"],
             "amount": scale_participant["charges"]["amount"]["normalized_value"]},
            ev, "Reaction scale normalized from the displayed limiting-reactant amount.",
            raw={"reference": scale_participant["molecule_id"],
                 "amount": scale_participant["charges"]["amount"]["raw_value"]})
    r["steps"] = list(steps)
    r["outcomes"] = list(outcomes)
    r["classification"] = classification(ev, class_mode)
    r["notes"] = list(notes)
    r["review"] = {
        "status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
        "unresolved": list(unresolved or ["Complete structures, exact atom correspondence, and omitted SI procedure details require expert review."])}
    r["completeness"] = "partial"
    bundle["reactions"].append(r)
    return r


def set_solvent(r, label, volume, ev, basis_mid="m_1a"):
    r["solvent_conditions"]["composition"] = f(label, [label], ev)
    if volume is not None:
        r["solvent_conditions"]["total_volume"] = f(f"{volume} mL", q(volume, "mL", raw_text=f"{volume} mL"), ev)
    if volume and (value_amount := next((p["charges"]["amount"]["normalized_value"] for p in r["participants"] if p["molecule_id"] == basis_mid), None)):
        if isinstance(value_amount, dict) and value_amount.get("unit") == "mmol":
            conc = value_amount["number"] / volume
            r["solvent_conditions"]["concentration"] = inferred(
                q(conc, "M", raw_text=f"derived {conc:g} M"), ev,
                "Computed from displayed limiting-substrate mmol divided by displayed total solvent volume.", raw="displayed scale/volume")
            r["solvent_conditions"]["concentration_basis"] = inferred(
                basis_mid, ev, "Limiting aldimine substrate used as concentration basis.", raw=basis_mid)


def annulation_participants(ald, allene, product, catalyst, solvent, ev, catalyst_molpct,
                            ald_amount=.1, allene_amount=.12, allene_equiv=None,
                            base_amount=None, base_equiv=1.1, ald_mass=None, allene_mass=None,
                            catalyst_amount=None):
    parts = [
        participant("p_aldimine", f"m_{ald}", "reactant", ev, amount=q(ald_amount, "mmol", raw_text=f"{ald_amount} mmol"),
                    mass=q(ald_mass, "mg", raw_text=f"{ald_mass} mg") if ald_mass else None,
                    equiv=q(1.0, "equiv", raw_text="1.0 equiv"), reference=f"m_{ald}"),
        participant("p_allenoate", f"m_{allene}", "reactant", ev,
                    amount=q(allene_amount, "mmol", raw_text=f"{allene_amount} mmol"),
                    mass=q(allene_mass, "mg", raw_text=f"{allene_mass} mg") if allene_mass else None,
                    equiv=q(allene_equiv, "equiv", raw_text=f"{allene_equiv} equiv") if allene_equiv else None,
                    reference=f"m_{ald}" if allene_equiv else None),
        participant("p_base", "m_K2CO3", "reagent", ev,
                    amount=q(base_amount, "mmol", raw_text=f"{base_amount} mmol") if base_amount else None,
                    equiv=q(base_equiv, "equiv", raw_text=f"{base_equiv} equiv"), reference=f"m_{ald}"),
        participant("p_catalyst", catalyst, "catalyst", ev,
                    amount=q(catalyst_amount, "mmol", raw_text=f"{catalyst_amount} mmol") if catalyst_amount else None,
                    molpct=q(catalyst_molpct, "mol%", raw_text=f"{catalyst_molpct} mol%"), reference=f"m_{ald}"),
        participant("p_solvent", solvent, "solvent", ev),
    ]
    if product:
        parts.append(participant("p_product", f"m_{product}", "product", ev, product_role="main product"))
    return parts


# Page 1: current generic reaction plus two prior-work examples. The graphical abstract and Scheme 1 current-work arrow are one design, not duplicate experiments.
ev1 = [E[1]]
reaction("scheme1_current", 1, "generic_scheme", None, "graphical abstract and Scheme 1, This work",
         "phosphine-catalyzed stereoselective tandem annulation", "tandem annulation",
         participants=[participant("p1", "m_generic_1", "reactant", ev1), participant("p2", "m_generic_2", "reactant", ev1), participant("p3", "m_generic_3", "product", ev1, product_role="main product")],
         class_mode="annulation", unresolved=["Generic R-group structures and exact atom mapping remain unresolved."],
         notes=[f("The graphical abstract and Scheme 1 'This work' depict the same generic design and are not counted twice.",
                  "duplicate source depiction", ev1)])
reaction("scheme1_geib_2005", 1, "literature_example", None, "Scheme 1, Previous work: Geib (2005)",
         "double intramolecular cyclization", None, class_mode="unresolved",
         unresolved=["Prior literature example; not an experiment newly reported in this paper."])
reaction("scheme1_tong_prior", 1, "literature_example", None, "Scheme 1, Tong's work",
         "1,3-dipolar cycloaddition with delta-acetoxy allenoates", None, class_mode="unresolved",
         unresolved=["Prior literature example; generic substituents and exact literature conditions are not imported as current experiments."])

# Page 2 generic scope arrows.
ev2 = [E[2]]
reaction("scheme2_generic", 2, "generic_scheme", None, "Scheme 2 generic arrow", "substrate scope of allenoates 2",
         "tandem annulation",
         participants=[participant("p1", "m_1a", "reactant", ev2), participant("p2", "m_generic_2", "reactant", ev2), participant("p3", "m_generic_3", "product", ev2, product_role="main product")],
         class_mode="annulation", scheme_ids=["scheme1_current"])
reaction("scheme3_generic", 2, "generic_scheme", None, "Scheme 3 generic arrow", "substrate scope of aldimine esters 1",
         "tandem annulation",
         participants=[participant("p1", "m_generic_1", "reactant", ev2), participant("p2", "m_2a", "reactant", ev2), participant("p3", "m_generic_3", "product", ev2, product_role="main product")],
         class_mode="annulation", scheme_ids=["scheme1_current"])

# Table 1, each row is independent. Entry 11 receives the duplicate Scheme 3 3aa evidence and its explicit 4 h caption.
table = [
    (1, "m_PPh3", "m_PhMe", "rt", "55%", q(55, "%", raw_text="55%")),
    (2, "m_PEtPh2", "m_PhMe", "rt", "26%", q(26, "%", raw_text="26%")),
    (3, "m_PEt2Ph", "m_PhMe", "rt", "<10%", q(10, "%", "lt", "<10%")),
    (4, "m_PBu3", "m_PhMe", "rt", "trace", "trace"),
    (5, "m_PPh3", "m_THF", "rt", "31%", q(31, "%", raw_text="31%")),
    (6, "m_PPh3", "m_CHCl3", "rt", "43%", q(43, "%", raw_text="43%")),
    (7, "m_PPh3", "m_CH2Cl2", "rt", "20%", q(20, "%", raw_text="20%")),
    (8, "m_PPh3", "m_CH3CN", "rt", "27%", q(27, "%", raw_text="27%")),
    (9, "m_PPh3", "m_Et2O", "rt", "<10%", q(10, "%", "lt", "<10%")),
    (10, "m_PPh3", "m_acetone", "rt", "trace", "trace"),
    (11, "m_PPh3", "m_PhMe", q(40, "°C", raw_text="40 °C"), "83%", q(83, "%", raw_text="83%")),
    (12, "m_PPh3", "m_PhMe", q(60, "°C", raw_text="60 °C"), "62%", q(62, "%", raw_text="62%")),
    (13, "m_PPh3", "m_PhMe", q(80, "°C", raw_text="80 °C"), "42%", q(42, "%", raw_text="42%")),
]
for entry, cat, solv, temp, raw_yield, yval in table:
    parts = annulation_participants("1a", "2a", "3aa", cat, solv, ev2, 20,
                                    allene_equiv=1.2, base_amount=.11, ald_mass=28,
                                    allene_mass=31, catalyst_amount=.02)
    r = reaction(f"table1_entry{entry:02d}", 2, "specific_experiment", ["optimization"],
                 f"Table 1 entry {entry}", "optimization of the tandem annulation", "tandem annulation",
                 participants=parts,
                 steps=[step("addition_and_reaction", ev2, [p["participant_id"] for p in parts],
                             "A solution of 1a in 1 mL solvent was slowly added to 2a, catalyst, and K2CO3 in 1 mL solvent.",
                             temp, q(4, "h", raw_text="4 h") if entry == 11 else None)],
                 outcomes=[outcome("yield", "yield", ["m_3aa"], ev2, yval, "isolated", "isolated",
                                   value_fact=conflict([("Table 1: trace", "trace", ev2),
                                                        ("prose: PBu3 could not generate 3aa", "not generated", ev2)]) if entry == 4 else None)],
                 scheme_ids=["scheme1_current"], class_mode="annulation",
                 unresolved=["Reaction time is unreported in Table 1." if entry != 11 else "Workup and purification are unavailable because SI was not supplied."])
    set_solvent(r, molecules[solv]["occurrences"][0]["compound_label"]["normalized_value"], 2.0, ev2)
    r["solvent_conditions"]["total_volume"] = inferred(
        q(2.0, "mL", raw_text="derived 2.0 mL"), ev2,
        "Table 1 footnote reports 1 mL for the 2a/catalyst/base solution plus 1 mL for the 1a solution; 2.0 mL is their arithmetic sum, not an author-reported total.",
        raw="1 mL + 1 mL")
    r["procedure"]["addition_order"] = f(
        "To 2a, catalyst, and K2CO3 in solvent (1 mL), a solution of 1a in solvent (1 mL) was slowly added.",
        ["2a/catalyst/K2CO3 solution", "slow addition of 1a solution"], ev2)
    if entry == 11:
        r["notes"].append(f("Scheme 3 repeats 3aa at 83% under the same displayed conditions and supplies the 4 h scope-caption time.",
                            "duplicate evidence merged", ev2))
    if entry == 4:
        r["review"]["unresolved"].append("Table 1 reports trace for PBu3, while the prose says PBu3 could not generate 3aa; both statements are retained without selecting one.")
        r["notes"].append(conflict([("Table 1 entry 4: trace", "trace", ev2),
                                     ("prose: PBu3 could not generate 3aa", "not generated", ev2)]))

# Scheme 2 scope; 3an is a genuine no-product experiment with no fabricated product participant.
scope2 = [("2b", "3ab", 61), ("2c", "3ac", 75), ("2d", "3ad", 79), ("2e", "3ae", 81),
          ("2f", "3af", 75), ("2g", "3ag", 73), ("2h", "3ah", 68), ("2i", "3ai", 56),
          ("2j", "3aj", 53), ("2k", "3ak", 58), ("2l", "3al", 62), ("2n", None, 0)]
for allene, prod, y in scope2:
    parts = annulation_participants("1a", allene, prod, "m_PPh3", "m_PhMe", ev2, 20, base_equiv=1.1)
    pids = [p["participant_id"] for p in parts]
    outs = ([outcome("yield", "yield", [f"m_{prod}"], ev2, q(y, "%", raw_text=f"{y}%"))]
            if prod else [outcome("no_product", "no_reaction", [], ev2, "0% target product", scope="target compound not obtained")])
    if prod in {"3ab", "3ac", "3ad", "3ae", "3af", "3ag", "3ah"}:
        outs.append(outcome("diastereoselectivity", "selectivity", [f"m_{prod}"], ev2,
                            "single diastereomer", scope={"metric": "diastereoselectivity", "major_identity": "unassigned"}))
    r = reaction(f"scheme2_{prod or '3an_no_product'}", 2, "specific_experiment",
                 ["substrate_scope"] + (["control"] if not prod else []), f"Scheme 2 {prod or '3an'}",
                 "substrate scope of allenoates 2", "tandem annulation", parts,
                 [step("reaction", ev2, pids, "reaction under Scheme 2 conditions", q(40, "°C", raw_text="40 °C"), q(4, "h", raw_text="4 h"))],
                 outs, ["scheme2_generic"], "annulation" if prod else "none",
                 unresolved=["Yield measurement basis is not reported in Scheme 2; workup/purification require unavailable SI."] if prod else ["Detection limit/method for the reported 0% target product is unreported."])
    set_solvent(r, "PhMe", 2.0, ev2)

# Scheme 3 additional scope. 3aa is already Table 1 entry 11 and is not duplicated.
scope3 = [("1b", "3ba", 68), ("1c", "3ca", 75), ("1d", "3da", 73), ("1e", "3ea", 78),
          ("1f", "3fa", 70), ("1g", "3ga", 63), ("1h", "3ha", 59), ("1i", "3ia", 57),
          ("1j", "3ja", 30), ("1k", "3ka", 56)]
for ald, prod, y in scope3:
    parts = annulation_participants(ald, "2a", prod, "m_PPh3", "m_PhMe", ev2, 20, base_equiv=1.1)
    pids = [p["participant_id"] for p in parts]
    r = reaction(f"scheme3_{prod}", 2, "specific_experiment", ["substrate_scope"], f"Scheme 3 {prod}",
                 "substrate scope of aldimine esters 1", "tandem annulation", parts,
                 [step("reaction", ev2, pids, "reaction under Scheme 3 conditions", q(40, "°C", raw_text="40 °C"), q(4, "h", raw_text="4 h"))],
                 [outcome("yield", "yield", [f"m_{prod}"], ev2, q(y, "%", raw_text=f"{y}%"))],
                 ["scheme3_generic"], "annulation",
                 unresolved=["Yield measurement basis and individual product procedure are not reported in the supplied main text."])
    set_solvent(r, "PhMe", 2.0, ev2, basis_mid=f"m_{ald}")

# Page 3 controls: prose explicitly says both gave no product under standard conditions.
ev3 = [E[3]]
for label in ("4", "5"):
    parts = annulation_participants(label, "2a", None, "m_PPh3", "m_PhMe", ev3, 20,
                                    base_equiv=1.1, allene_equiv=1.2)
    r = reaction(f"scheme4_control_{label}", 3, "specific_experiment", ["control"],
                 f"Scheme 4 control with {label}", "control experiment", "attempted tandem annulation", parts,
                 [step("reaction", ev3, [p["participant_id"] for p in parts],
                       "treated with 2a under the standard condition", q(40, "°C", raw_text="40 °C"), q(4, "h", raw_text="standard 4 h"))],
                 [outcome("no_product", "no_reaction", [], ev3, "No product", scope="target annulation product not obtained")],
                 ["scheme1_current"], "none",
                 unresolved=["No detection method/limit is reported; the 4 h time is inherited from the explicitly referenced standard conditions."])
    set_solvent(r, "PhMe", 2.0, ev3, basis_mid=f"m_{label}")
    inherited_derivation = "Scheme 4 says 'under our standard condition'; the 0.10/0.12 mmol substrate charges, 1.2 equiv allenoate, 2.0 mL solvent, and 4 h time are inherited from the displayed standard scope conditions."
    for pid in ("p_aldimine", "p_allenoate"):
        p = next(x for x in r["participants"] if x["participant_id"] == pid)
        for key in ("amount", "equivalents"):
            fact = p["charges"][key]
            if fact["normalized_value"] is not None:
                fact["acquisition_method"] = "inherited_general_condition"
                fact["derivation"] = inherited_derivation
    r["steps"][0]["time"] = f("standard condition (4 h)", q(4, "h", raw_text="standard 4 h"), ev3,
                                     method="inherited_general_condition", derivation=inherited_derivation)
    r["solvent_conditions"]["total_volume"] = f("standard condition (2.0 mL)", q(2.0, "mL", raw_text="standard 2.0 mL"), ev3,
                                                     method="inherited_general_condition", derivation=inherited_derivation)
    r["procedure"]["general_procedure_reference"] = f("under our standard condition", "Scheme 2/3 standard conditions", ev3)

# Proposed mechanism is a non-experimental reaction record.
reaction("scheme5_proposed_mechanism", 3, "proposed_mechanism", None, "Scheme 5",
         "proposed mechanism", "author-proposed tandem annulation mechanism",
         participants=[participant("p1", "m_generic_1", "reactant", ev3), participant("p2", "m_generic_2", "reactant", ev3), participant("p3", "m_generic_3", "product", ev3, product_role="main product")],
         class_mode="annulation",
         unresolved=["Intermediates I–VI and H-shift/elimination sequence are author proposals, not isolated or independently demonstrated experiments."],
         notes=[f("The authors suggest addition–elimination, C-C attack, intramolecular addition, H-shift, nucleophilic attack, and elimination.",
                  "author-proposed interpretation", ev3)])

# Scheme 6 generic arrow and eleven specific asymmetric experiments.
reaction("scheme6_generic", 3, "generic_scheme", None, "Scheme 6 generic arrow",
         "preliminary asymmetric phosphine-catalyzed annulation", "asymmetric tandem annulation",
         participants=[participant("p1", "m_generic_1", "reactant", ev3), participant("p2", "m_generic_2", "reactant", ev3), participant("pc", "m_C10", "catalyst", ev3), participant("p3", "m_generic_3", "product", ev3, product_role="main product")],
         class_mode="annulation", scheme_ids=["scheme1_current"],
         unresolved=["C10 structure is depicted but not converted to an unsafe complete stereochemical SMILES."])
asym = [
    ("1a", "2a", "3aa", 79, 86), ("1a", "2c", "3ac", 66, 81), ("1a", "2d", "3ad", 62, 84),
    ("1a", "2f", "3af", 62, 80), ("1a", "2h", "3ah", 57, 82), ("1a", "2l", "3al", 52, 78),
    ("1a", "2m", "3am", 76, 91), ("1b", "2m", "3bm", 71, 80), ("1c", "2m", "3cm", 68, 72),
    ("1d", "2m", "3dm", 56, 84), ("1d", "2a", "3da", 73, 59),
]
for ald, allene, prod, y, ee in asym:
    parts = annulation_participants(ald, allene, prod, "m_C10", "m_CHCl3", ev3, 10, base_equiv=1.1)
    pids = [p["participant_id"] for p in parts]
    if prod == "3aa":
        yfact = conflict([("Scheme 6: 79%", q(79, "%", raw_text="79%"), ev3),
                          ("prose: 75%", q(75, "%", raw_text="75%"), ev3)])
        eefact = conflict([("Scheme 6: 86% ee", q(86, "%", raw_text="86% ee"), ev3),
                           ("prose: 87% ee", q(87, "%", raw_text="87% ee"), ev3)])
    else:
        yfact = f(f"{y}%", q(y, "%", raw_text=f"{y}%"), ev3)
        eefact = f(f"{ee}% ee", q(ee, "%", raw_text=f"{ee}% ee"), ev3)
    outs = [
        outcome("yield", "yield", [f"m_{prod}"], ev3, value_fact=yfact),
        outcome("dr", "selectivity", [f"m_{prod}"], ev3, {"ratio": "20:1", "qualifier": "gt", "raw_text": ">20:1 dr"},
                scope={"metric": "dr", "major_identity": "unassigned", "minor_identity": "unassigned"}),
        outcome("Z_E", "selectivity", [f"m_{prod}"], ev3, {"ratio": "25:1", "qualifier": "gt", "raw_text": ">25:1 Z/E"},
                scope={"metric": "Z/E", "order": ["Z", "E"]}),
        outcome("ee", "selectivity", [f"m_{prod}"], ev3, value_fact=eefact,
                scope={"metric": "ee", "major_enantiomer_absolute_configuration": "unassigned"}),
    ]
    r = reaction(f"scheme6_{prod}", 3, "specific_experiment", ["substrate_scope"], f"Scheme 6 {prod}",
                 "preliminary asymmetric phosphine-catalyzed annulation", "asymmetric tandem annulation", parts,
                 [step("reaction", ev3, pids, "reaction under Scheme 6 conditions", q(0, "°C", raw_text="0 °C"), q(12, "h", raw_text="12 h"))],
                 outs, ["scheme6_generic"], "annulation",
                 unresolved=["Yield and selectivity measurement methods are unreported in the supplied main text; absolute configuration is unassigned."] +
                            (["Scheme 6 and prose conflict for both 3aa yield (79 vs 75%) and ee (86 vs 87%); no value selected."] if prod == "3aa" else []))
    set_solvent(r, "CHCl3", 2.0, ev3, basis_mid=f"m_{ald}")

# Scheme 7 derivatizations.
parts6 = [participant("p_start", "m_scheme7_start_RH", "reactant", ev3), participant("p_NaBH4", "m_NaBH4", "reagent", ev3),
          participant("p_THF", "m_THF", "solvent", ev3), participant("p_product", "m_6", "product", ev3, product_role="main product")]
r6 = reaction("scheme7_product6", 3, "specific_experiment", ["derivatization"], "Scheme 7 upper arrow",
              "NaBH4 reduction", "ester reduction", parts6,
              [step("reduction", ev3, [p["participant_id"] for p in parts6], "treatment with NaBH4 in THF")],
              [outcome("yield", "yield", ["m_6"], ev3, q(80, "%", raw_text="80%"))],
              class_mode="reduction",
              unresolved=["Scheme 7 drawing shows an R=H, CO2CHPh2 starting structure, while prose calls it 3aa; identity conflict remains unresolved.",
                          "Reagent amount, temperature, time, yield basis, workup and purification are unreported in the supplied main text."],
              notes=[conflict([("Scheme 7 drawing: R=H and CO2CHPh2", "m_scheme7_start_RH", ev3),
                               ("prose: treatment of 3aa", "m_3aa", ev3)])])
set_solvent(r6, "THF", None, ev3, basis_mid="m_scheme7_start_RH")

parts7 = [participant("p_start", "m_3bm", "reactant", ev3), participant("p_boronic", "m_PhB_OH2", "reagent", ev3),
          participant("p_Pd", "m_Pd_OAc2", "catalyst", ev3), participant("p_XPhos", "m_XPhos", "ligand", ev3),
          participant("p_base", "m_K2CO3", "reagent", ev3), participant("p_tol", "m_PhMe", "solvent", ev3),
          participant("p_product", "m_7", "product", ev3, product_role="main product")]
r7 = reaction("scheme7_product7", 3, "specific_experiment", ["derivatization"], "Scheme 7 lower arrow",
              "coupling reaction", "Suzuki-type aryl coupling", parts7,
              [step("coupling", ev3, [p["participant_id"] for p in parts7], "PhB(OH)2, Pd(OAc)2, XPhos, K2CO3 in toluene",
                    q(110, "°C", raw_text="110 °C"), q(12, "h", raw_text="12 h"))],
              [outcome("yield", "yield", ["m_7"], ev3, q(71, "%", raw_text="71%"))],
              class_mode="coupling",
              unresolved=["Reagent loadings, yield basis, workup and purification are unreported in the supplied main text."])
set_solvent(r7, "PhMe", None, ev3, basis_mid="m_3bm")


# Non-synthetic source context, including characterization and unprovided SI.
bundle["source_context"] = [
    {"context_id": f"{DOC}:figure1_bioactive_context", "document_id": DOC, "type": "biological_context",
     "evidence_ids": [E[1]], "note": "Figure 1 shows representative bioactive chromenopyrrole structures; these are not syntheses performed in this paper."},
    {"context_id": f"{DOC}:xray_3aa", "document_id": DOC, "type": "characterization",
     "evidence_ids": [E[2]], "note": "The structure of 3aa is stated to be confirmed by X-ray crystallographic analysis; no additional synthesis record created."},
    {"context_id": f"{DOC}:si_notice", "document_id": DOC, "type": "supplementary_information_notice",
     "evidence_ids": [E[4]], "note": "The article advertises SI with procedures and spectra, but no SI file was supplied or explicitly associated in this workspace."},
    {"context_id": f"{DOC}:ccdc", "document_id": DOC, "type": "crystallographic_accession",
     "evidence_ids": [E[4], E[5]], "note": "CCDC 1976344 and 2020535 are accession context, not independent synthetic experiments."},
    {"context_id": f"{DOC}:references_pages4_5", "document_id": DOC, "type": "bibliography",
     "evidence_ids": [E[4], E[5]], "note": "Bibliographic references were reviewed as context and not imported as current-paper experiments."},
]

bundle["molecules"] = list(molecules.values())
st.write(RUN / "submitted.json", bundle)

# Raw source inventory and literal cross-region hints.
inventory = {
    "document_id": DOC,
    "title": TITLE,
    "pages": {
        "1": {"figures": ["graphical abstract", "Figure 1", "Scheme 1"], "current_records": [f"{DOC}:scheme1_current"],
              "other": [f"{DOC}:scheme1_geib_2005", f"{DOC}:scheme1_tong_prior", "Figure 1 biological context"]},
        "2": {"tables": ["Table 1 entries 1-13"], "schemes": ["Scheme 2", "Scheme 3"],
              "notes": ["Scheme 3 3aa duplicates Table 1 entry 11"]},
        "3": {"schemes": ["Scheme 4 controls", "Scheme 5 proposed mechanism", "Scheme 6 asymmetric scope", "Scheme 7 transformations"],
              "notes": ["3aa Scheme/prose numerical conflict", "Scheme 7 product 6 starting-material identity conflict"]},
        "4": {"other": ["conclusion", "SI availability notice", "CCDC accessions", "references"]},
        "5": {"other": ["references and CCDC note continuation"]},
    },
    "limitations": ["No SI supplied.", "Dense structures retained without unsafe complete SMILES.", "All records remain needs_review."],
}
st.write(RUN / "raw_model/pages_0001_0005_extraction.json", inventory)
terms = ["1a", "2a", "3aa", "3am", "3bm", "ethyl ( E)-2-", "NaBH", "C10"]
st.write(RUN / "raw_model/source_mentions.json", st.source_mentions(json.loads((RUN / "dataset.json").read_text()), terms, 320))

checks = {}
for mid, m in molecules.items():
    smi = m["smiles"]["normalized_value"]
    if smi:
        try:
            checks[mid] = st.stereo_inspect(smi)
        except Exception as exc:
            checks[mid] = {"error": str(exc)}
st.write(RUN / "raw_model/stereo_checks.json", checks)

# Coverage contract: exactly the requested five pages, all attempted, with real page evidence.
page2_ids = [f"{DOC}:scheme2_generic", f"{DOC}:scheme3_generic"] + [f"{DOC}:table1_entry{i:02d}" for i in range(1, 14)]
page2_ids += [f"{DOC}:scheme2_{p or '3an_no_product'}" for _, p, _ in scope2]
page2_ids += [f"{DOC}:scheme3_{p}" for _, p, _ in scope3]
page3_ids = [f"{DOC}:scheme4_control_4", f"{DOC}:scheme4_control_5", f"{DOC}:scheme5_proposed_mechanism", f"{DOC}:scheme6_generic"]
page3_ids += [f"{DOC}:scheme6_{p}" for _, _, p, _, _ in asym]
page3_ids += [f"{DOC}:scheme7_product6", f"{DOC}:scheme7_product7"]
coverage = {
    "document_id": DOC, "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "graphical abstract and Scheme 1", "finding": "Current generic annulation plus two explicitly labeled prior-work examples; duplicate current-work depiction merged.", "evidence_ids": [E[1]]},
            {"region": "Figure 1 and surrounding prose", "finding": "Bioactive-scaffold context, not a synthesis performed in this paper.", "evidence_ids": [E[1]]}],
         "record_ids": [f"{DOC}:scheme1_current", f"{DOC}:scheme1_geib_2005", f"{DOC}:scheme1_tong_prior"],
         "unresolved": ["Generic R groups and exact atom mapping are unresolved."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Table 1 rows 1-13 and footnotes", "finding": "Thirteen independent optimization entries; entry 11 is duplicated by Scheme 3 product 3aa and merged.", "evidence_ids": [E[2]]},
            {"region": "Scheme 2 arrow, products 3ab-3al and 3an, caption", "finding": "Twelve allenoate-scope outcomes, including a 0% no-product experiment for 3an.", "evidence_ids": [E[2]]},
            {"region": "Scheme 3 arrow, products 3aa-3ka, caption", "finding": "Ten new aldimine-scope outcomes; 3aa is duplicate evidence for Table 1 entry 11.", "evidence_ids": [E[2]]}],
         "record_ids": page2_ids,
         "unresolved": ["Scope yield measurement basis is not reported; dense product structures require expert review."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Scheme 4 and prose", "finding": "Two true no-product control experiments under standard conditions.", "evidence_ids": [E[3]]},
            {"region": "Scheme 5 and mechanism prose", "finding": "Author-proposed mechanism, not an experimentally verified sequence.", "evidence_ids": [E[3]]},
            {"region": "Scheme 6 and asymmetric-annulation prose", "finding": "Eleven specific asymmetric experiments; Scheme/prose conflict for 3aa yield and ee retained.", "evidence_ids": [E[3]]},
            {"region": "Scheme 7 and transformation prose", "finding": "Two derivatizations; the product-6 starting-material identity conflicts between drawing and prose.", "evidence_ids": [E[3]]}],
         "record_ids": page3_ids,
         "unresolved": ["No SI procedures; selectivity methods and absolute configurations are unreported; two source conflicts remain open."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "conclusion and Associated Content", "finding": "No new reaction arrow/table; SI is advertised but was not supplied.", "evidence_ids": [E[4]]},
            {"region": "CCDC/accession and references", "finding": "Characterization/accession and bibliography only, not new experiments.", "evidence_ids": [E[4]]}],
         "record_ids": [], "unresolved": ["CCDC files and SI were not supplied and were not inferred as associated source documents."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "full printed page", "finding": "Bibliography and CCDC note continuation only; no current-paper reaction experiment.", "evidence_ids": [E[5]]}],
         "record_ids": [], "unresolved": []},
    ]}
st.write(RUN / "coverage/pages_0001_0005.json", coverage)

tasks = json.loads((RUN / "tasks.json").read_text())
for task in tasks:
    task["state"] = "attempted"
    task["reason"] = "Original page image and extracted text were inspected; see coverage/pages_0001_0005.json."
st.write(RUN / "tasks.json", tasks)

print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
                  "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in bundle["reactions"])}, ensure_ascii=False))
