#!/usr/bin/env python3
"""Build the source-derived v3.1 submission for Chem. Sci. 2020, 11, 2401.

This is retained as the raw machine extraction recipe.  It deliberately leaves
image-only structures unresolved rather than manufacturing SMILES.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_c6189dacd18a4bec"
E = {p: f"{DOC}:page{p}:image" for p in range(1, 7)}
ET = {p: f"{DOC}:page{p}:text" for p in range(1, 7)}
TITLE = ("Enantioselective nickel-catalyzed arylative and alkenylative "
         "intramolecular 1,2-allylations of tethered allene–ketones")


def f(raw=None, normalized=None, ev=(), method="original_extraction", derivation=None,
      status=None):
    return st.fact(raw, normalized, ev, method=method, derivation=derivation,
                   status=status)


def inferred(value, ev, derivation):
    return f(None, value, ev, method="inference", derivation=derivation)


def conflict(raw_values, alternatives):
    return {
        "raw_value": raw_values,
        "normalized_value": None,
        "acquisition_method": "original_extraction",
        "evidence_ids": list(dict.fromkeys(eid for alt in alternatives for eid in alt[2])),
        "data_status": "conflict",
        "verification": "unreviewed",
        "derivation": None,
        "alternatives": [f(raw, normalized, ev) for raw, normalized, ev in alternatives],
    }


def na(ev=()):
    return f(None, None, ev, method="not_applicable", status="not_applicable")


def q(number, unit, qualifier=None, raw_text=None):
    out = {"number": number, "unit": unit}
    if qualifier:
        out["qualifier"] = qualifier
    if raw_text:
        out["raw_text"] = raw_text
    return out


def mol(mid, label, page, *, stereo="unresolved", assignment=None,
        mixture=None, note=None):
    ev = [E[page]]
    m = st.new_molecule(f"{DOC}:mol:{mid}", DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    m["names"]["original"] = f(label, label, ev)
    m["names"]["english"] = f(label, label, ev)
    m["names"]["abbreviations"] = f()
    m["names"]["aliases"] = f()
    m["cas_numbers"] = f()
    m["representation_type"] = "unresolved"
    m["smiles"] = f()
    m["generic_template"] = None
    m["stereochemistry"]["atom_centers"] = f()
    m["stereochemistry"]["double_bonds"] = f()
    m["stereochemistry"]["axial"] = f()
    m["stereochemistry"]["mixture"] = f(mixture, mixture, ev) if mixture else f()
    if assignment:
        assign_ev = ev + ([E[6]] if page != 6 else [])
        m["stereochemistry"]["assignment"] = f(assignment, assignment, assign_ev)
    else:
        m["stereochemistry"]["assignment"] = f()
    m["stereochemistry"]["chirality_status"] = inferred(
        stereo, ev, "Status follows explicit wedge/dash, racemate, or mixture information in the source drawing/text; no unshown configuration was added.")
    m["stereochemistry"]["geometry_status"] = inferred(
        "unspecified", ev, "No E/Z descriptor is printed for this source label; this is not a claim of absence of alkene stereochemistry.")
    m["stereochemistry"]["representation_status"] = inferred(
        "not_encoded", ev, "The source structure was inspected visually but no source-verified isomeric SMILES is available from the supplied main text.")
    m["stereochemistry"]["isomeric_smiles"] = f()
    m["stereochemistry"]["identity_key"] = f()
    identity_note = note or "Identity is retained by the paper label and source drawing; exact machine-readable structure remains unresolved without supplied SI/source molfile."
    m["identity_note"] = f(identity_note, identity_note, ev)
    m["external_properties"] = []
    return m


def generic_mol(mid, label, page, fixed, variables):
    ev = [E[page]]
    m = mol(mid, label, page, note="Generic R-group structure shown in the paper; not expanded into a concrete molecule.")
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": f(fixed, fixed, ev),
        "variables": [],
        "instance_mapping": [],
        "reason": f("Generic structure with R/X/n variables in the source figure.",
                    "Generic structure with unresolved variables.", ev),
    }
    for label_var, definition in variables:
        m["generic_template"]["variables"].append({
            "label": label_var,
            "attachment_points": f(),
            "definition": f(definition, definition, ev),
            "allowed_values": f(),
            "constraints": f(),
        })
    return m


def participant(pid, mid, role, ev, *, amount=None, equiv=None, molpct=None,
                ref=None, stage="reaction", product_role=None):
    p = st.new_participant(pid, f"{DOC}:mol:{mid}", role, ev, stage=stage)
    if product_role:
        p["product_role"] = f(product_role, product_role, ev)
    if amount is not None:
        p["charges"]["amount"] = f(amount[0], amount[1], ev)
    if equiv is not None:
        p["charges"]["equivalents"] = f(equiv[0], equiv[1], ev)
        p["equivalence_reference"] = f(ref, f"{DOC}:mol:{ref}", ev)
    if molpct is not None:
        p["charges"]["mol_percent"] = f(molpct[0], molpct[1], ev)
        p["equivalence_reference"] = f(ref, f"{DOC}:mol:{ref}", ev)
    return p


def add_step(r, ev, pids, temp=None, time=None, operation="reaction"):
    s = st.new_step(f"{r['reaction_id']}:step1", ev)
    s["participant_ids"] = list(pids)
    s["operation"] = f(operation, operation, ev)
    if temp:
        s["temperature"] = f(temp[0], temp[1], ev)
    if time:
        s["time"] = f(time[0], time[1], ev)
    r["steps"].append(s)


def outcome(r, oid, kind, target_mid, raw, normalized, ev, *, basis=None,
            yield_type=None, scope=None):
    o = st.new_outcome(f"{r['reaction_id']}:{oid}", kind,
                       [f"{DOC}:mol:{target_mid}"] if target_mid else [], ev)
    o["value"] = f(raw, normalized, ev)
    if basis:
        o["measurement_basis"] = f(basis, basis, ev)
    if yield_type:
        o["yield_type"] = f(yield_type, yield_type, ev)
    if scope:
        o["scope"] = f(scope, scope, ev)
    if kind in ("yield", "product_amount"):
        o["product_role"] = f("main product", "main_product", ev)
    r["outcomes"].append(o)


def arylative_classification(r, ev):
    changes = [
        {"change_id": "c1", "atoms": [{"element": "C", "source_role": "allene"}, {"element": "C", "source_role": "aryl_or_alkenyl"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "c2", "atoms": [{"element": "C", "source_role": "allyl"}, {"element": "C", "source_role": "ketone_carbonyl"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "c3", "atoms": [{"element": "O", "source_role": "former_carbonyl_oxygen"}, {"element": "H", "source_role": "proton_source_unresolved"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
    ]
    r["classification"]["bond_changes"] = [f(None, c, ev, method="structure_derivation", derivation="Net bond change inferred from the source reactant/product drawings; atom mapping is not claimed.") for c in changes]
    labels = [
        {"category_number": 2, "code": "C-C", "supporting_change_ids": ["c1", "c2"]},
        {"category_number": 24, "code": "O-H", "supporting_change_ids": ["c3"]},
    ]
    r["classification"]["labels"] = f(None, labels, ev, method="structure_derivation", derivation="Two net C-C formations and alcohol O-H formation are visible across the generic/source transformation; detailed mapping and proton source remain unresolved.")
    r["classification"]["coverage"] = inferred("partial", ev, "Image-derived net transformation without complete atom mapping or coproduct accounting.")
    r["classification"]["limitations"] = ["No complete atom mapping was asserted.", "Proton source and boron-containing coproducts are not identified in the supplied main text."]


def not_applicable_classification(r, ev, reason):
    r["classification"]["bond_changes"] = []
    r["classification"]["labels"] = inferred([], ev, reason)
    r["classification"]["coverage"] = inferred("not_applicable", ev, reason)
    r["classification"]["limitations"] = [reason]


def base_reaction(suffix, ev, record_type="specific_experiment", categories=None,
                  reaction_name="nickel-catalyzed arylative intramolecular 1,2-allylation"):
    r = st.new_reaction(f"{DOC}:{suffix}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories or ["other"], categories or ["other"], ev)
    r["location"] = f({"source": suffix}, {"source": suffix}, ev)
    r["reaction_name"] = f(reaction_name, reaction_name, ev)
    r["transformation_type"] = inferred(
        "carbonickelation of an allene followed by intramolecular 1,2-allylation of a ketone",
        ev, "Normalized from the title, abstract, and source reaction drawing.")
    r["review"]["unresolved"] = [
        "No supporting-information file was explicitly supplied; detailed operations, workup, purification, and full characterization remain unavailable.",
        "Image-only structures are retained by source labels and have not been promoted to source-verified SMILES.",
        "Net bond classification is partial because no complete atom mapping or proton-source accounting is available.",
    ]
    return r


def set_solvent(r, ev, names, ratio=None):
    r["solvent_conditions"]["composition"] = f(names, names, ev)
    if ratio:
        r["solvent_conditions"]["mixture_ratio"] = f(ratio, ratio, ev)


def add_standard_conditions(r, ev, substrate_mid, reagent_mid, product_mid,
                            *, scale=None, solvent="TFE", ligand="L3",
                            ligand_loading=5, ni_loading=5, time_h=24,
                            reagent_equiv=1.5):
    amount = (f"{scale} mmol", q(scale, "mmol")) if scale is not None else None
    r["participants"] = [
        participant("substrate", substrate_mid, "reactant", ev, amount=amount),
        participant("coupling_partner", reagent_mid, "reactant", ev,
                    equiv=(f"{reagent_equiv} equiv.", q(reagent_equiv, "equiv")), ref=substrate_mid),
        participant("nickel", "ni_acetate", "catalyst", ev,
                    molpct=(f"{ni_loading} mol%", q(ni_loading, "mol%")), ref=substrate_mid),
        participant("ligand", ligand.lower(), "ligand", ev,
                    molpct=(f"{ligand_loading} mol%", q(ligand_loading, "mol%")), ref=substrate_mid),
        participant("solvent", solvent.lower().replace("/", "_"), "solvent", ev),
        participant("product", product_mid, "product", ev, product_role="main product"),
    ]
    if scale is not None:
        r["reaction_scale"] = f(f"{scale} mmol of substrate", {"reference_molecule_id": f"{DOC}:mol:{substrate_mid}", "amount": q(scale, "mmol")}, ev)
    set_solvent(r, ev, [solvent])
    add_step(r, ev, [p["participant_id"] for p in r["participants"]],
             ("80 °C", q(80, "°C")), (f"{time_h} h", q(time_h, "h")))
    arylative_classification(r, ev)


def product_assignment(label):
    if label in {"2a", "5a", "5d"}:
        return "Relative and absolute configuration determined by X-ray crystallography (paper note 17); machine-readable stereochemistry not encoded here."
    return "Relative and absolute configuration assigned by analogy in paper note 17; machine-readable stereochemistry not encoded here."


def build():
    dataset_path = ROOT / "runs/current/dataset.json"
    data = json.loads(dataset_path.read_text())
    data["documents"][0]["title"] = TITLE
    st.write(dataset_path, data)
    bundle = st.empty(data["documents"])

    molecules = {}
    def add(m):
        molecules[m["molecule_id"].split(":mol:", 1)[1]] = m
    # Generic schemes and previous-work examples.
    add(generic_mol("litA_substrate", "Scheme 1A generic allene substrate (ref. 8a)", 1,
                    "cyclohexadienone-tethered allene", [("R", "source variable"), ("X", "NTs or O")]))
    add(generic_mol("litA_product", "Scheme 1A generic product (ref. 8a)", 1,
                    "arylative 1,4-allylation product", [("R", "source variable"), ("Ar", "aryl")]))
    add(generic_mol("litB_substrate", "Scheme 1B generic activated allene (ref. 8b)", 1,
                    "activated allene bearing R groups", [("R1", "source variable"), ("R2", "source variable"), ("R3", "source variable")]))
    add(generic_mol("litB_product", "Scheme 1B generic annulation product (ref. 8b)", 1,
                    "benzannulated product", [("R1", "source variable"), ("R2", "source variable"), ("R3", "source variable")]))
    add(generic_mol("generic_1", "tethered allene–ketone 1 (generic)", 1,
                    "allene tethered to alpha-ketoamide", [("R1", "ketone substituent"), ("R2", "nitrogen substituent"), ("Ar", "arylboronic-acid-derived aryl")]))
    add(generic_mol("generic_2", "cyclized tertiary-alcohol product 2 (generic)", 1,
                    "pyrrolidin-2-one with tertiary alcohol and alkenyl substituent", [("R1", "ketone substituent"), ("R2", "nitrogen substituent"), ("Ar", "aryl/heteroaryl")]))

    # Common reagents/solvents/ligands retained by names; CAS and SMILES are not inferred.
    common = [
        ("phenylboronic_acid", "PhB(OH)2 / phenylboronic acid", 2),
        ("ni_acetate", "Ni(OAc)2·4H2O", 2), ("l1", "L1", 2),
        ("l2", "L2 / (S)-i-Pr-PHOX", 2), ("l3", "L3 / (S)-t-Bu-PHOX", 2),
        ("l4", "L4 / (R)-Ph-PHOX", 2), ("mecn", "MeCN", 2),
        ("1,4-dioxane", "1,4-dioxane", 2),
        ("mecn_1,4-dioxane", "MeCN/1,4-dioxane", 2),
        ("tfe", "TFE / 2,2,2-trifluoroethanol", 2),
        ("vinyl_bf3k", "potassium vinyltrifluoroborate", 3),
        ("can", "CAN", 4), ("water", "H2O", 4),
        ("naio4", "NaIO4", 4), ("osmate", "K2[OsO2(OH)4]", 4),
        ("thf", "THF", 4),
    ]
    for mid, label, page in common:
        add(mol(mid, label, page))

    # Specific substrates and products.
    for label in ["1a", "1b", "1c", "1d", "1e", "1f", "1g"]:
        add(mol(label, f"tethered allene–ketone {label}", 2 if label in {"1a", "1g"} else 3,
                note="Specific source drawing/label retained; allene axial and other stereochemical status are not inferred from absent descriptors."))
    products2 = ["2a", "2b", "2c", "2d", "2e", "2f", "2g", "2h", "2i", "2j", "2k", "2l", "2m", "2n", "2o", "2p", "2q", "2r", "3"]
    for label in products2:
        page = 2 if label in {"2a", "2r"} else 3
        mix = "7:1 mixture of inseparable diastereomers" if label == "2q" else None
        stereo = "mixture" if label == "2q" else "specified"
        add(mol(label, f"product {label}", page, stereo=stereo,
                assignment=product_assignment("2a" if label == "2a" else label), mixture=mix))
    add(mol("ent_2a", "ent-2a", 4, stereo="specified",
            assignment="Minor enantiomer ent-2a shown in the proposed stereochemical model."))
    add(mol("rac_2a", "racemic 2a", 2, stereo="racemic", mixture="racemic product"))
    add(mol("unspecified_2a", "2a (enantiomeric composition unreported)", 2, stereo="unspecified"))

    # Scope-specific aryl/heteroaryl reagents are not assigned names beyond visible substituents.
    ar_reagents = {
        "2g": "para-acetoxyphenylboronic acid", "2h": "para-chlorophenylboronic acid",
        "2i": "para-vinylphenylboronic acid", "2j": "meta-fluorophenylboronic acid",
        "2k": "meta-cyanophenylboronic acid", "2l": "ortho-methylphenylboronic acid",
        "2m": "ortho-fluorophenylboronic acid", "2n": "2-aminophenylboronic acid pinacol ester",
        "2o": "bromo-methylphenylboronic acid", "2p": "dichlorophenylboronic acid",
        "2q": "furylboronic acid",
    }
    for label, name in ar_reagents.items():
        add(mol(f"ar_{label}", name, 3))

    # Table 3 substrates/products and reagents.
    for label in ["5a", "5b", "5c", "5d", "5e", "5f", "5g", "5h", "5i", "5j", "5k"]:
        mix = None
        add(mol(label, f"product {label}", 3, stereo="specified",
                assignment=product_assignment(label), mixture=mix))
        substrate_label = "4d" if label == "5d" else f"4a–4i member leading to {label} (exact label not printed beside product)"
        add(mol(f"sub_{label}", substrate_label, 3,
                note="The table shows the generic 4a–4i family but does not print an unambiguous individual substrate label beside this product."))
        add(mol(f"ar_{label}", f"arylboronic acid used for {label}", 3,
                note="Identity retained from the product substituent/source drawing; no external-name expansion was used."))
    add(mol("sample_5d_4d", "7.7:1 inseparable mixture of 5d and starting allene 4d", 3,
            stereo="mixture", mixture="7.7:1 mixture of product 5d and starting allene 4d",
            note="This sample-composition identity is separate from molecular product 5d; the source does not report a separate recovery percentage."))

    # Failed controls and derivatization product.
    add(mol("4j", "substrate 4j", 4))
    add(mol("4k", "substrate 4k", 4))
    add(mol("12", "product 12", 4, stereo="specified",
            assignment="Configuration depicted as retained from 2a; no separate absolute assignment is stated."))

    bundle["molecules"] = list(molecules.values())
    reactions = []

    # Scheme 1 previous work A/B.
    for code, sub, prod, raw_name in [
        ("scheme1A_literature", "litA_substrate", "litA_product", "Enantioselective Ni-catalyzed arylative intramolecular 1,4-allylations (ref. 8a)"),
        ("scheme1B_literature", "litB_substrate", "litB_product", "Ni-catalyzed annulations of activated allenes (ref. 8b)"),
    ]:
        r = base_reaction(code, [E[1]], record_type="literature_example", reaction_name=raw_name)
        r["participants"] = [participant("substrate", sub, "reactant", [E[1]]), participant("product", prod, "product", [E[1]], product_role="literature product")]
        r["review"]["unresolved"] = ["Cited prior work is preserved as a literature example, not a current-paper experiment.", "Generic structures and conditions were not expanded beyond the source scheme."]
        r["classification"]["coverage"] = inferred("unresolved", [E[1]], "Prior-work generic scheme not atom-mapped in this extraction.")
        reactions.append(r)

    # Scheme 1C current generic transformation.
    r = base_reaction("scheme1C_generic", [E[1]], record_type="generic_scheme")
    r["participants"] = [participant("substrate", "generic_1", "reactant", [E[1]]), participant("coupling_partner", "phenylboronic_acid", "reactant", [E[1]]), participant("product", "generic_2", "product", [E[1]], product_role="generic product")]
    arylative_classification(r, [E[1]])
    reactions.append(r)

    # Optimization table 1, one row per experiment.
    opt_rows = [
        (1, None, "mecn_1,4-dioxane", "MeCN/1,4-dioxane", "3:2", 32, None, "rac_2a"),
        (2, "L1", "mecn_1,4-dioxane", "MeCN/1,4-dioxane", "3:2", 77, None, "unspecified_2a"),
        (3, "L2", "mecn_1,4-dioxane", "MeCN/1,4-dioxane", "3:2", 55, 90, "2a"),
        (4, "L3", "mecn_1,4-dioxane", "MeCN/1,4-dioxane", "3:2", 78, 98, "2a"),
        (5, "L4", "mecn_1,4-dioxane", "MeCN/1,4-dioxane", "3:2", 61, -91, "ent_2a"),
        (6, "L3", "mecn", "MeCN", None, 78, 98, "2a"),
        (7, "L3", "1,4-dioxane", "1,4-dioxane", None, 59, 98, "2a"),
        (8, "L3", "tfe", "TFE", None, ">99", 96, "2a"),
    ]
    for entry, ligand, solvent_mid, solvent_name, ratio, yld, ee, product_mid in opt_rows:
        ev = [E[2]]
        r = base_reaction(f"table1_entry{entry}", ev, categories=["optimization"])
        parts = [
            participant("substrate", "1a", "reactant", ev, amount=("0.10 mmol", q(0.10, "mmol"))),
            participant("boronic_acid", "phenylboronic_acid", "reactant", ev, equiv=("1.5 equiv.", q(1.5, "equiv")), ref="1a"),
            participant("nickel", "ni_acetate", "catalyst", ev, molpct=("5 mol%", q(5, "mol%")), ref="1a"),
            participant("solvent", solvent_mid, "solvent", ev),
            participant("product", product_mid, "product", ev, product_role="main product"),
        ]
        if ligand:
            parts.insert(3, participant("ligand", ligand.lower(), "ligand", ev, molpct=("5 mol%", q(5, "mol%")), ref="1a"))
        r["participants"] = parts
        r["reaction_scale"] = f("0.10 mmol of 1a", {"reference_molecule_id": f"{DOC}:mol:1a", "amount": q(0.10, "mmol")}, ev)
        set_solvent(r, ev, [solvent_name], ratio)
        add_step(r, ev, [p["participant_id"] for p in parts], ("80 °C", q(80, "°C")), ("24 h", q(24, "h")))
        yval = q(99, "%", "gt", ">99%") if yld == ">99" else q(yld, "%")
        outcome(r, "yield", "yield", product_mid, f"{yld}%", yval, ev,
                basis="1H NMR analysis of the crude reaction using 1,3,5-dimethoxybenzene as internal standard", yield_type="NMR yield")
        if ee is not None:
            raw_ee = f"{ee}%" if ee >= 0 else f"{ee}% (major product was the enantiomer of 2a)"
            outcome(r, "ee", "selectivity", product_mid, raw_ee, q(abs(ee), "%"), ev,
                    basis="HPLC analysis on a chiral stationary phase", scope="enantiomeric excess; negative sign in entry 5 denotes the opposite major enantiomer")
        if entry == 1:
            outcome(r, "racemic", "selectivity", product_mid, "racemic 2a", "racemic", ev,
                    scope="enantiomeric composition")
            outcome(r, "dr", "selectivity", product_mid, "one observable diastereomer", "single observable diastereomer", ev, scope="diastereoselectivity")
        else:
            outcome(r, "dr", "selectivity", product_mid, ">19:1 dr", {"ratio": ">19:1"}, ev, scope="diastereomeric ratio shown in the Table 1 reaction graphic")
        arylative_classification(r, ev)
        reactions.append(r)

    # Equation (1): unprotected substrate 1g.
    r = base_reaction("eq1_2r", [E[2]], categories=["substrate_scope"])
    add_standard_conditions(r, [E[2]], "1g", "phenylboronic_acid", "2r", scale=None)
    outcome(r, "yield", "yield", "2r", "90%", q(90, "%"), [E[2]])
    outcome(r, "dr", "selectivity", "2r", ">19:1 dr", {"ratio": ">19:1"}, [E[2]], scope="diastereomeric ratio")
    outcome(r, "ee", "selectivity", "2r", "99% ee", q(99, "%"), [E[2]], basis="HPLC analysis on a chiral stationary phase", scope="enantiomeric excess")
    reactions.append(r)

    # Table 2, including alternative MeCN condition for 2f as a separate experiment.
    table2 = [
        ("2a", "1a", "phenylboronic_acid", "TFE", "99", "96", None),
        ("2b", "1b", "phenylboronic_acid", "TFE", "99", "97", None),
        ("2c", "1c", "phenylboronic_acid", "TFE", "99", "98", None),
        ("2d", "1d", "phenylboronic_acid", "TFE", ">99", ">99", None),
        ("2e", "1e", "phenylboronic_acid", "TFE", ">99", ">99", None),
        ("2f", "1f", "phenylboronic_acid", "TFE", "99", "87", None),
        ("2f", "1f", "phenylboronic_acid", "MeCN", "65", ">99", "b"),
        ("2g", "1a", "ar_2g", "TFE", ">99", "97", None),
        ("2h", "1a", "ar_2h", "TFE", "99", "98", None),
        ("2i", "1a", "ar_2i", "TFE", "71", "99", None),
        ("2j", "1a", "ar_2j", "TFE", ">99", "98", None),
        ("2k", "1a", "ar_2k", "TFE", "90", "99", None),
        ("2l", "1a", "ar_2l", "TFE", "69", "84", None),
        ("2m", "1a", "ar_2m", "TFE", "68", "97", None),
        ("2n", "1a", "ar_2n", "TFE", "48", "56", "c"),
        ("2o", "1a", "ar_2o", "TFE", ">99", "98", None),
        ("2p", "1a", "ar_2p", "TFE", ">99", "98", None),
        ("2q", "1a", "ar_2q", "TFE", "91", "91", "d"),
    ]
    seen = {}
    for product, substrate, reagent, solvent, yraw, eeraw, footnote in table2:
        seen[product] = seen.get(product, 0) + 1
        suffix = f"table2_{product}" + ("_mecn" if seen[product] > 1 else "")
        ev = [E[3]]
        r = base_reaction(suffix, ev, categories=["substrate_scope"])
        add_standard_conditions(r, ev, substrate, reagent, product, scale=0.30, solvent=solvent)
        yn = q(99, "%", "gt", ">99%") if yraw == ">99" else q(int(yraw), "%")
        en = q(99, "%", "gt", ">99% ee") if eeraw == ">99" else q(int(eeraw), "%")
        outcome(r, "yield", "yield", product, f"{yraw}%", yn, ev, basis="isolated product", yield_type="isolated yield")
        outcome(r, "ee", "selectivity", product, f"{eeraw}% ee", en, ev, basis="HPLC analysis on a chiral stationary phase", scope="enantiomeric excess")
        dr_raw = "7:1 dr; isolated as a 7:1 mixture of inseparable diastereomers" if product == "2q" else ">19:1 dr"
        dr_norm = {"ratio": "7:1"} if product == "2q" else {"ratio": ">19:1"}
        outcome(r, "dr", "selectivity", product, dr_raw, dr_norm, ev, basis="1H NMR analysis of the crude reaction mixture", scope="diastereomeric ratio")
        if footnote == "c":
            r["notes"].append(f("Conducted using the pinacol boronate instead of the boronic acid.", "pinacol boronate used", ev))
        if product == "2q":
            r["evidence_ids"] = [E[2], E[3]]
            r["outcomes"][0]["value"] = conflict(
                ["Table 2: 91%", "page-2 prose: 93%"],
                [("91%", q(91, "%"), [E[3]]), ("93%", q(93, "%"), [E[2]])],
            )
            r["outcomes"][1]["value"] = conflict(
                ["Table 2: 91% ee", "page-2 prose: 98% ee"],
                [("91% ee", q(91, "%"), [E[3]]), ("98% ee", q(98, "%"), [E[2]])],
            )
            r["notes"].append(f(
                "Source conflict: Table 2 gives 2q as 91% yield and 91% ee, whereas the preceding prose gives 93% yield and 98% ee.",
                "table/prose conflict retained without choosing a value", [E[2], E[3]]))
            r["review"]["unresolved"].append("Resolve the 2q table/prose yield and ee conflict against an authoritative correction or SI; neither value was selected here.")
        reactions.append(r)

    # Equation (2): alkenylative example, text begins on p2 and scheme is p3.
    r = base_reaction("eq2_product3", [E[2], E[3]], categories=["substrate_scope"],
                      reaction_name="nickel-catalyzed alkenylative intramolecular 1,2-allylation")
    add_standard_conditions(r, [E[3]], "1a", "vinyl_bf3k", "3", scale=None, solvent="TFE")
    r["evidence_ids"] = [E[2], E[3]]
    outcome(r, "yield", "yield", "3", "65%", q(65, "%"), [E[3]])
    outcome(r, "dr", "selectivity", "3", ">19:1 dr", {"ratio": ">19:1"}, [E[3]], scope="diastereomeric ratio")
    outcome(r, "ee", "selectivity", "3", "93% ee", q(93, "%"), [E[3]], basis="HPLC analysis on a chiral stationary phase", scope="enantiomeric excess")
    reactions.append(r)

    # Table 3, plus alternative MeCN condition for 5c.
    table3 = [
        ("5a", "TFE", 92, 84, 24, 5, None),
        ("5b", "TFE", 96, 83, 24, 5, None),
        ("5c", "TFE", 62, 75, 24, 5, None),
        ("5c", "MeCN", 54, 80, 24, 5, "d"),
        ("5d", "TFE", 73, 99, 24, 5, "b"),
        ("5e", "TFE", 91, 90, 24, 5, None),
        ("5f", "MeCN", 33, 90, 24, 5, "d"),
        ("5g", "MeCN", 47, 93, 24, 5, "d"),
        ("5h", "MeCN", 34, 89, 24, 5, "d"),
        ("5i", "TFE", 50, 99, 48, 10, "c,e"),
        ("5j", "MeCN", 67, 85, 24, 10, "d,e"),
        ("5k", "MeCN", 76, 76, 48, 10, "c,d,e"),
    ]
    seen = {}
    for product, solvent, yld, ee, hours, loading, footnotes in table3:
        seen[product] = seen.get(product, 0) + 1
        suffix = f"table3_{product}" + ("_mecn" if seen[product] > 1 else "")
        ev = [E[3]]
        r = base_reaction(suffix, ev, categories=["substrate_scope"])
        add_standard_conditions(r, ev, f"sub_{product}", f"ar_{product}", product,
                                scale=None, solvent=solvent, ligand_loading=loading,
                                ni_loading=loading, time_h=hours)
        outcome(r, "yield", "yield", product, f"{yld}%", q(yld, "%"), ev,
                basis="isolated product", yield_type="isolated yield")
        outcome(r, "ee", "selectivity", product, f"{ee}% ee", q(ee, "%"), ev,
                basis="HPLC analysis on a chiral stationary phase", scope="enantiomeric excess")
        outcome(r, "dr", "selectivity", product, ">19:1 dr", {"ratio": ">19:1"}, ev,
                basis="1H NMR analysis of the crude reaction mixture", scope="diastereomeric ratio")
        if product == "5d":
            r["participants"].append(participant("isolated_mixture", "sample_5d_4d", "product", ev, product_role="inseparable isolated sample"))
            outcome(r, "sample_composition", "other", "sample_5d_4d",
                    "7.7:1 inseparable mixture of 5d and starting allene 4d",
                    {"ratio": "7.7:1", "order": ["5d", "4d"]}, ev,
                    scope="isolated sample composition")
            r["notes"].append(f("A 7.7:1 inseparable mixture of 5d and starting allene 4d was obtained; the yield of 5d was adjusted accordingly.",
                                "yield adjusted for 7.7:1 product/starting-allene mixture", ev))
            r["review"]["unresolved"].append("The source reports an adjusted isolated yield for 5d but does not provide a separate recovered-starting-material percentage.")
        reactions.append(r)

    # Figure 2: two real no-reaction controls.
    for label in ["4j", "4k"]:
        ev = [E[4]]
        r = base_reaction(f"fig2_{label}_no_reaction", ev, categories=["control"])
        r["participants"] = [participant("substrate", label, "reactant", ev)]
        outcome(r, "no_reaction", "no_reaction", label, "unreacted starting material returned", "no reaction; starting material returned", ev,
                scope="author-reported failed arylative cyclization")
        not_applicable_classification(r, ev, "No observed transformation; no bond-formation label is assigned.")
        r["review"]["unresolved"] = ["Exact charges, conditions, and recovery amount are not printed in Figure 2/main-text prose; no SI was supplied."]
        reactions.append(r)

    # Scheme 2: proposed mechanism, not an experiment.
    r = base_reaction("scheme2_proposed_mechanism", [E[4]], record_type="proposed_mechanism",
                      categories=["other"], reaction_name="possible catalytic cycle")
    r["participants"] = []
    r["transformation_type"] = f("Possible catalytic cycle involving species 6–11 and competing 10a/10b conformations.",
                                  "proposed nickel catalytic cycle", [E[4]])
    r["classification"]["coverage"] = inferred("unresolved", [E[4]], "A proposed elementary-step mechanism is not treated as an observed net experiment.")
    r["review"]["unresolved"] = ["Author-proposed mechanism; five-coordinate alternatives are explicitly not excluded."]
    reactions.append(r)

    # Scheme 3 derivatizations.
    ev = [E[4]]
    r = base_reaction("scheme3_deprotection_2a_to_2r", ev, categories=["derivatization"],
                      reaction_name="CAN removal of the 4-methoxyphenyl protecting group")
    r["participants"] = [
        participant("substrate", "2a", "reactant", ev),
        participant("can", "can", "reagent", ev, equiv=("2.5 equiv.", q(2.5, "equiv")), ref="2a"),
        participant("mecn", "mecn", "solvent", ev), participant("water", "water", "solvent", ev),
        participant("product", "2r", "product", ev, product_role="main product"),
    ]
    set_solvent(r, ev, ["MeCN", "H2O"], "1:1")
    add_step(r, ev, [p["participant_id"] for p in r["participants"]], ("−10 °C", q(-10, "°C")), ("30 min", q(30, "min")), "CAN deprotection")
    outcome(r, "yield", "yield", "2r", "98%", q(98, "%"), ev)
    changes = [
        {"change_id": "d1", "atoms": [{"element": "C"}, {"element": "N"}], "before_order": 1, "after_order": 0, "change_type": "broken"},
        {"change_id": "d2", "atoms": [{"element": "N"}, {"element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
    ]
    r["classification"]["bond_changes"] = [f(None, c, ev, method="structure_derivation", derivation="Net change inferred from removal of N-PMP and formation of N-H in source drawings.") for c in changes]
    r["classification"]["labels"] = f(None, [{"category_number": 29, "code": "N-H", "supporting_change_ids": ["d2"]}], ev, method="structure_derivation", derivation="Only newly formed N-H is classified; broken C-N is retained as a bond change but not a formation label.")
    r["classification"]["coverage"] = inferred("partial", ev, "Oxidation coproducts and atom mapping are unavailable.")
    reactions.append(r)

    r = base_reaction("scheme3_oxidative_cleavage_2a_to_12", ev, categories=["derivatization"],
                      reaction_name="oxidative cleavage of the alkene of 2a")
    r["participants"] = [
        participant("substrate", "2a", "reactant", ev),
        participant("periodate", "naio4", "reagent", ev, equiv=("6.0 equiv.", q(6.0, "equiv")), ref="2a"),
        participant("osmate", "osmate", "reagent", ev, molpct=("0.10 mol%", q(0.10, "mol%")), ref="2a"),
        participant("thf", "thf", "solvent", ev), participant("water", "water", "solvent", ev),
        participant("product", "12", "product", ev, product_role="main product"),
    ]
    set_solvent(r, ev, ["THF", "H2O"], "5:2")
    pids = [p["participant_id"] for p in r["participants"]]
    add_step(r, ev, pids, ("0 °C", q(0, "°C")), None,
             "initiate oxidative cleavage at 0 °C")
    warm = st.new_step(f"{r['reaction_id']}:step2", ev)
    warm["participant_ids"] = pids
    warm["operation"] = f("0 °C to RT", "warm from 0 °C to room temperature", ev)
    warm["temperature"] = f("RT", "room temperature", ev)
    warm["time"] = f("18 h", q(18, "h"), ev)
    r["steps"].append(warm)
    outcome(r, "yield", "yield", "12", "91%", q(91, "%"), ev)
    changes = [
        {"change_id": "o1", "atoms": [{"element": "C"}, {"element": "O"}], "before_order": 0, "after_order": 2, "change_type": "formed"},
        {"change_id": "o2", "atoms": [{"element": "C"}, {"element": "C"}], "before_order": 1, "after_order": 0, "change_type": "broken"},
    ]
    r["classification"]["bond_changes"] = [f(None, c, ev, method="structure_derivation", derivation="Net alkene cleavage/carbonyl formation inferred from Scheme 3 drawings; omitted coproduct prevents complete mapping.") for c in changes]
    r["classification"]["labels"] = f(None, [{"category_number": 14, "code": "DB-C-O", "supporting_change_ids": ["o1"]}], ev, method="structure_derivation", derivation="New carbonyl C=O bond is visible; C-C cleavage is retained but not a formation label.")
    r["classification"]["coverage"] = inferred("partial", ev, "Only the depicted main product is mapped conceptually; coproduct identity is omitted.")
    reactions.append(r)

    bundle["reactions"] = reactions
    bundle["source_context"] = [
        {"context_id": f"{DOC}:context:figure1", "kind": "background_structures", "evidence_ids": [E[1]], "raw_value": "Fig. 1 shows pramanicin, norsecurinamine A, cytochalasin Z10, and a herbicide as motivation; these are not experiments from this paper."},
        {"context_id": f"{DOC}:context:esi_notice", "kind": "supplement_notice", "evidence_ids": [E[1]], "raw_value": "Electronic supplementary information is stated to be available, but no SI file was supplied in this session; no association was inferred."},
        {"context_id": f"{DOC}:context:mechanistic_discussion", "kind": "author_interpretation", "evidence_ids": [E[4]], "raw_value": "The authors favor 10a over 10b on steric grounds and explicitly state that alternative five-coordinate nickel complexes cannot be excluded."},
        {"context_id": f"{DOC}:context:stereo_assignment", "kind": "stereochemical_assignment_note", "evidence_ids": [E[6]], "raw_value": "Relative and absolute configurations of 2a, 5a, and 5d were determined by X-ray crystallography; remaining products were assigned by analogy (CCDC 1945864–1945866)."},
        {"context_id": f"{DOC}:context:bibliography", "kind": "bibliography", "evidence_ids": [E[5], E[6]], "raw_value": "Pages 5–6 are acknowledgements and cited references; no new current-paper experimental arrows or tables appear there."},
    ]
    return bundle


if __name__ == "__main__":
    bundle = build()
    st.write(ROOT / "runs/current/submitted.json", bundle)
    inventory = {
        "document_id": DOC,
        "title": TITLE,
        "pages_visually_inspected": [1, 2, 3, 4, 5, 6],
        "molecule_count": len(bundle["molecules"]),
        "reaction_count": len(bundle["reactions"]),
        "reaction_ids": [r["reaction_id"] for r in bundle["reactions"]],
        "notes": [
            "All six rendered source pages were inspected at original detail.",
            "Structures that could not be encoded confidently from the main-text image remain unresolved.",
            "No SI association was inferred.",
        ],
    }
    st.write(ROOT / "runs/current/raw_model/extraction_inventory.json", inventory)
    ids_by_page = {
        p: [r["reaction_id"] for r in bundle["reactions"] if E[p] in r["evidence_ids"]]
        for p in range(1, 7)
    }
    coverage = {
        "document_id": DOC,
        "semantic_accuracy_verified": False,
        "pages": [
            {
                "page": 1,
                "state": "attempted",
                "regions": [
                    {"region": "title, abstract, and bibliographic header", "finding": "Current-paper identity and overall nickel-catalyzed allene–ketone cyclization claim; no concrete experiment here.", "evidence_ids": [E[1], ET[1]]},
                    {"region": "Fig. 1", "finding": "Background natural-product/herbicide structures, not current-paper experiments; retained as source_context.", "evidence_ids": [E[1]]},
                    {"region": "Scheme 1A and 1B", "finding": "Two cited prior-work reaction arrows; recorded as literature_example, not current-paper experiments.", "evidence_ids": [E[1]]},
                    {"region": "Scheme 1C", "finding": "Generic current-work transformation and proposed allylnickel intermediates; one generic_scheme record.", "evidence_ids": [E[1]]},
                    {"region": "ESI footnote", "finding": "ESI availability is stated, but no SI file was supplied or associated in this session.", "evidence_ids": [E[1], ET[1]]}
                ],
                "record_ids": ids_by_page[1],
                "unresolved": ["Generic R groups and Scheme 1 structures are not concrete molecule identities.", "The stated ESI was not supplied."]
            },
            {
                "page": 2,
                "state": "attempted",
                "regions": [
                    {"region": "Table 1 reaction graphic, ligand panel, rows 1–8, and footnotes", "finding": "Eight separate optimization experiments with 0.10 mmol 1a, row-specific ligand/solvent/yield/ee, NMR yield basis, and HPLC ee basis.", "evidence_ids": [E[2], ET[2]]},
                    {"region": "right-column Table 2 discussion", "finding": "Scope narrative cross-checked against page-3 structures, including 2f solvent alternative, 2q 7:1 dr, and 2n pinacol boronate exception.", "evidence_ids": [E[2], ET[2]]},
                    {"region": "equation (1)", "finding": "Specific 1g to 2r experiment: 90% yield, >19:1 dr, 99% ee.", "evidence_ids": [E[2]]},
                    {"region": "bottom-right equation (2) lead-in", "finding": "Text identifies 1a plus potassium vinyltrifluoroborate to product 3; scheme and conditions continue on page 3.", "evidence_ids": [E[2], ET[2]]}
                ],
                "record_ids": ids_by_page[2],
                "unresolved": ["Main-text page does not provide detailed operations/workup for optimization or equation (1).", "Equation (1) yield measurement basis is not stated on this page."]
            },
            {
                "page": 3,
                "state": "attempted",
                "regions": [
                    {"region": "Table 2 general reaction and product grid 2a–2q", "finding": "Seventeen products plus a distinct MeCN experiment for 2f: 18 substrate-scope records; isolated yields, HPLC ee, and crude-mixture NMR dr retained.", "evidence_ids": [E[3], ET[3]]},
                    {"region": "Table 2 footnotes a–d", "finding": "0.30 mmol scale, isolated-yield basis, analytical methods, MeCN exception, pinacol-boronate exception, and 2q inseparable 7:1 mixture applied only where indicated.", "evidence_ids": [E[3]]},
                    {"region": "equation (2)", "finding": "Specific alkenylative experiment affording 3 in 65% yield, >19:1 dr, 93% ee.", "evidence_ids": [E[3]]},
                    {"region": "Table 3 general reaction and product grid 5a–5k", "finding": "Eleven products plus a distinct MeCN experiment for 5c: 12 records with product-specific yield/ee and footnoted time, solvent, and catalyst-loading variations.", "evidence_ids": [E[3], ET[3]]},
                    {"region": "Table 3 footnotes a–e", "finding": "Isolated-yield and analytical bases retained; 5d adjusted-yield/product-starting-material mixture noted; 48 h, MeCN, and 10 mol% exceptions structured.", "evidence_ids": [E[3]]}
                ],
                "record_ids": ids_by_page[3],
                "unresolved": ["Table 3 does not print an unambiguous individual 4a–4i label beside most products; reaction-specific substrate placeholders preserve that uncertainty.", "Exact machine-readable structures were not promoted from page images."]
            },
            {
                "page": 4,
                "state": "attempted",
                "regions": [
                    {"region": "Fig. 2 structures 4j and 4k plus adjacent prose", "finding": "Two real failed/no-reaction controls returning unreacted starting materials; no products or expected-bond labels invented.", "evidence_ids": [E[4], ET[4]]},
                    {"region": "Scheme 2 and mechanistic discussion", "finding": "Author-proposed catalytic cycle and stereochemical model, recorded once as proposed_mechanism; alternative five-coordinate complexes explicitly remain possible.", "evidence_ids": [E[4], ET[4]]},
                    {"region": "Scheme 3 upper arrow", "finding": "CAN deprotection of 2a to 2r: 2.5 equiv CAN, MeCN/H2O 1:1, −10 °C, 30 min, 98% yield.", "evidence_ids": [E[4], ET[4]]},
                    {"region": "Scheme 3 lower arrow", "finding": "Oxidative alkene cleavage of 2a to 12: NaIO4 6.0 equiv, K2[OsO2(OH)4] 0.10 mol%, THF/H2O 5:2, 0 °C to RT, 18 h, 91% yield.", "evidence_ids": [E[4], ET[4]]}
                ],
                "record_ids": ids_by_page[4],
                "unresolved": ["Figure 2 does not state exact charges, conditions, or recovered-starting-material percentages.", "Scheme 3 yield measurement bases and detailed workups are unreported in the supplied main text."]
            },
            {
                "page": 5,
                "state": "attempted",
                "regions": [
                    {"region": "Acknowledgements and references 1–15", "finding": "Bibliography/publication context only; cited reactions are not new experiments from this paper and no new reaction arrow/table is present.", "evidence_ids": [E[5], ET[5]]}
                ],
                "record_ids": ids_by_page[5],
                "unresolved": []
            },
            {
                "page": 6,
                "state": "attempted",
                "regions": [
                    {"region": "references 15–22", "finding": "Bibliography and author notes; no new reaction experiment or table.", "evidence_ids": [E[6], ET[6]]},
                    {"region": "note 17", "finding": "Stereochemical provenance: 2a, 5a, and 5d by X-ray; remaining products assigned by analogy. Retained as molecule assignment evidence and source_context.", "evidence_ids": [E[6], ET[6]]},
                    {"region": "notes 20–21", "finding": "Author mechanistic qualifications concerning directing groups and associative/dissociative ligand substitution; no separate experiment.", "evidence_ids": [E[6], ET[6]]}
                ],
                "record_ids": ids_by_page[6],
                "unresolved": ["Crystallographic CIF/ESI files were not supplied, so absolute structures were not independently encoded or checked."]
            }
        ]
    }
    st.write(ROOT / "runs/current/coverage/pages_0001_0006.json", coverage)
    print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"])}, ensure_ascii=False))
