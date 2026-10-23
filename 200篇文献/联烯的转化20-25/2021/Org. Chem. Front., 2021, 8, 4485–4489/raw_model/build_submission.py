#!/usr/bin/env python3
"""Build the source-derived v3.1 submission for doc_3fe9d12641be7de4.

This session-local builder only serializes facts manually read from PDF pages 1--5.
It does not modify the extraction tools, schemas, or validation rules.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st  # noqa: E402


DOC = "doc_3fe9d12641be7de4"
RUN = ROOT / "runs/current"
EV = {
    p: [f"{DOC}:page{p}:image", f"{DOC}:page{p}:text"] for p in range(1, 6)
}


def extracted(raw, normalized, ev):
    return st.fact(raw, normalized, ev)


def inferred(normalized, ev, derivation, raw=None, method="inference"):
    return st.fact(raw, normalized, ev, method=method, derivation=derivation)


def missing(ev=()):
    return st.fact(evidence=ev)


def not_applicable(ev=()):
    return st.fact(evidence=ev, method="not_applicable", status="not_applicable")


def quantity(number, unit, raw_text=None, qualifier=None):
    value = {"number": number, "unit": unit}
    if raw_text is not None:
        value["raw_text"] = raw_text
    if qualifier is not None:
        value["qualifier"] = qualifier
    return value


bundle = json.loads((RUN / "submitted.json").read_text())
bundle["molecules"] = []
bundle["reactions"] = []
bundle["source_context"] = []
molecules = {}


def add_molecule(
    mid,
    label,
    ev,
    *,
    name=None,
    smiles=None,
    chirality="unresolved",
    geometry="unresolved",
    representation="unresolved",
    atom_centers=None,
    axial=None,
    assignment=None,
    mixture=None,
    double_bonds=None,
    identity_note=None,
):
    if mid in molecules:
        return molecules[mid]
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = extracted(label, label, ev)
    m["occurrences"][0]["evidence_ids"] = list(ev)
    m["names"]["original"] = extracted(name or label, name or label, ev)
    m["names"]["english"] = extracted(name, name, ev) if name else missing(ev)
    m["names"]["abbreviations"] = missing(ev)
    m["names"]["aliases"] = missing(ev)
    m["cas_numbers"] = missing(ev)
    m["external_properties"] = []
    m["identity_note"] = (
        extracted(identity_note, identity_note, ev) if identity_note else missing(ev)
    )
    if smiles:
        diag = st.stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = inferred(
            smiles,
            ev,
            "Candidate graph manually transcribed from the labelled source drawing/name; RDKit parsing and stereo diagnostics do not establish source-structure correctness.",
            raw={"source_structure_label": label},
            method="structure_derivation",
        )
        m["stereochemistry"]["isomeric_smiles"] = inferred(
            diag["canonical_isomeric_smiles"],
            ev,
            "Canonical isomeric SMILES calculated from the candidate registry SMILES with the standard stereo diagnostic.",
            method="structure_derivation",
        )
        m["stereochemistry"]["identity_key"] = inferred(
            diag["stereo_graph_key"],
            ev,
            "Stereo graph key calculated from the candidate registry SMILES; it is not a sample identity.",
            method="structure_derivation",
        )
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing(ev)
    stereo = m["stereochemistry"]
    stereo["atom_centers"] = (
        extracted(atom_centers, atom_centers, ev) if atom_centers is not None else missing(ev)
    )
    stereo["double_bonds"] = (
        extracted(double_bonds, double_bonds, ev) if double_bonds is not None else missing(ev)
    )
    stereo["axial"] = extracted(axial, axial, ev) if axial is not None else missing(ev)
    stereo["mixture"] = extracted(mixture, mixture, ev) if mixture is not None else missing(ev)
    stereo["assignment"] = (
        extracted(assignment, assignment, ev) if assignment is not None else missing(ev)
    )
    for key, value in (
        ("chirality_status", chirality),
        ("geometry_status", geometry),
        ("representation_status", representation),
    ):
        stereo[key] = inferred(
            value,
            ev,
            "Status assigned from the displayed molecular graph and explicit stereochemical annotations; absence of an assignment was not converted to racemic or achiral.",
        )
    molecules[mid] = m
    bundle["molecules"].append(m)
    return m


def template_fact(raw, normalized, ev):
    return extracted(raw, normalized, ev)


def add_generic_molecule(mid, label, ev, fixed, variables, reason, product=False):
    m = add_molecule(
        mid,
        label,
        ev,
        name=label,
        chirality="unspecified" if product else "unresolved",
        geometry="unresolved",
        representation="not_encoded",
        atom_centers=(
            [{"description": "spiro quaternary stereogenic center drawn without R/S assignment"}]
            if product
            else None
        ),
        identity_note=reason,
    )
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": template_fact(fixed, fixed, ev),
        "variables": [
            {
                "label": v,
                "attachment_points": template_fact(f"Source drawing attachment for {v}", f"Source drawing attachment for {v}", ev),
                "definition": template_fact(desc, desc, ev),
                "allowed_values": missing(ev),
                "constraints": missing(ev),
            }
            for v, desc in variables
        ],
        "instance_mapping": [],
        "reason": template_fact(reason, reason, ev),
    }
    return m


# Core substrates/products and the three current-work generic structures.
add_molecule(
    f"{DOC}:mol:1a",
    "1a",
    EV[2] + EV[3],
    name="allenic ketone 1a",
    smiles="C=C=C(COC(C)=O)C(=O)c1ccccc1",
    chirality="achiral",
    geometry="not_applicable",
    representation="not_applicable",
    axial=[{"allene": "terminal CH2 makes axial chirality not applicable"}],
)
add_molecule(
    f"{DOC}:mol:1a_prime",
    "1a′",
    EV[2],
    name="allenoate 1a′",
    smiles="C=C=C(COC(C)=O)C(=O)OCc1ccccc1",
    chirality="achiral",
    geometry="not_applicable",
    representation="not_applicable",
    axial=[{"allene": "terminal CH2 makes axial chirality not applicable"}],
)
add_molecule(
    f"{DOC}:mol:2a",
    "2a",
    EV[2] + EV[3],
    name="N-Boc oxindole 2a",
    smiles="CC(C)(C)OC(=O)N1C(=O)Cc2ccccc21",
    chirality="achiral",
    geometry="not_applicable",
    representation="not_applicable",
)
add_molecule(
    f"{DOC}:mol:3a",
    "3a",
    EV[2] + EV[3],
    name="enantioenriched spirocyclopentene oxindole 3a",
    chirality="unspecified",
    geometry="not_applicable",
    representation="not_encoded",
    atom_centers=[{"center": "spiro quaternary carbon", "configuration": "not assigned in main text"}],
    assignment="ee is reported per experiment, but the major enantiomer is not assigned R/S in pages 1–5",
    identity_note="Paper-local product identity; exact stereochemical SMILES was not safely established from the main-text raster.",
)
add_molecule(
    f"{DOC}:mol:3a_rac",
    "3a-rac",
    EV[2],
    name="racemic spirocyclopentene oxindole 3a-rac",
    chirality="racemic",
    geometry="not_applicable",
    representation="not_encoded",
    atom_centers=[{"center": "spiro quaternary carbon", "configuration": "racemic label"}],
    mixture="racemic",
)
add_molecule(
    f"{DOC}:mol:3a_prime_rac",
    "3a′-rac",
    EV[2],
    name="racemic allenoate-derived annulation product 3a′-rac",
    chirality="racemic",
    geometry="not_applicable",
    representation="not_encoded",
    atom_centers=[{"center": "spiro quaternary carbon", "configuration": "racemic label"}],
    mixture="racemic",
)
add_generic_molecule(
    f"{DOC}:mol:generic_allenic_ketone_1",
    "generic allenic ketone 1",
    EV[2] + EV[3],
    "Drawn allenic-ketone C4 synthon bearing OAc and a ketone substituent",
    [("R1", "substituent on the ketone-side group as drawn")],
    "The source uses an R-group drawing, so no concrete SMILES is asserted.",
)
add_generic_molecule(
    f"{DOC}:mol:generic_oxindole_2",
    "generic oxindole 2",
    EV[2] + EV[3],
    "Drawn oxindole core",
    [("R2", "benzene-ring substituent as drawn"), ("PG", "nitrogen protecting group")],
    "The source uses R2 and PG variables, so no concrete SMILES is asserted.",
)
add_generic_molecule(
    f"{DOC}:mol:generic_product_3",
    "generic spirocyclopentene oxindole 3",
    EV[2] + EV[3],
    "Drawn 3-spirocyclopentene-2-oxindole product core",
    [("R1", "ketone-side substituent"), ("R2", "oxindole-ring substituent"), ("PG", "nitrogen protecting group")],
    "The source uses an R-group product drawing; the spiro stereocenter is shown but not assigned R/S.",
    product=True,
)


# Simple named agents/solvents with candidate graphs; complex catalysts remain unresolved.
simple = [
    ("PPh3", "triphenylphosphine", "P(c1ccccc1)(c1ccccc1)c1ccccc1", EV[2] + EV[3], "achiral"),
    ("Na2CO3", "sodium carbonate", "[Na+].[Na+].[O-]C([O-])=O", EV[2] + EV[3], "not_applicable"),
    ("CH2Cl2", "dichloromethane", "ClCCl", EV[3], "achiral"),
    ("THF", "tetrahydrofuran", "C1CCOC1", EV[3], "achiral"),
    ("toluene", "toluene", "Cc1ccccc1", EV[2] + EV[3], "achiral"),
    ("dioxane_1_4", "1,4-dioxane", "O1CCOCC1", EV[3], "achiral"),
    ("TBME", "tert-butyl methyl ether (TBME)", "COC(C)(C)C", EV[2] + EV[3], "achiral"),
    ("ozone", "ozone", "[O-][O+]=O", EV[3], "not_applicable"),
    ("methanol", "methanol", "CO", EV[3], "achiral"),
    ("TFA", "trifluoroacetic acid", "O=C(O)C(F)(F)F", EV[3], "achiral"),
    ("methylamine", "methylamine", "CN", EV[3], "achiral"),
    ("NaH", "sodium hydride", "[Na+].[H-]", EV[3], "not_applicable"),
    ("prenyl_bromide", "1-bromo-3-methylbut-2-ene", "BrCC=C(C)C", EV[3], "achiral"),
]
for key, name, smi, ev, chiral in simple:
    add_molecule(
        f"{DOC}:mol:{key}",
        key,
        ev,
        name=name,
        smiles=smi,
        chirality=chiral,
        geometry="not_applicable",
        representation="not_applicable",
    )
add_molecule(f"{DOC}:mol:ether", "Ether", EV[3], name="Ether (identity not expanded in source)")
for catalyst in [f"P{i}" for i in range(1, 9)]:
    add_molecule(
        f"{DOC}:mol:{catalyst}",
        catalyst,
        EV[3],
        name=f"chiral phosphine catalyst {catalyst}",
        chirality="specified",
        geometry="not_applicable",
        representation="not_encoded",
        atom_centers="stereogenic centers are depicted by wedge/dash bonds in Table 1 catalyst drawings",
        assignment="source drawing specifies relative/absolute stereochemical depiction without a safely encoded full SMILES here",
    )


# Failure-control oxindoles.
for key, label, name, smi in [
    ("oxindole_NH", "unprotected oxindole", "unprotected oxindole", "O=C1Cc2ccccc2N1"),
    ("oxindole_NMe", "N-methyl oxindole", "N-methyl oxindole", "CN1C(=O)Cc2ccccc21"),
    ("oxindole_NBn", "N-benzyl oxindole", "N-benzyl oxindole", "O=C1Cc2ccccc2N1Cc1ccccc1"),
]:
    add_molecule(
        f"{DOC}:mol:{key}", label, EV[2], name=name, smiles=smi,
        chirality="achiral", geometry="not_applicable", representation="not_applicable",
    )


# Scope-specific identities. One common fully drawn reactant is retained per row;
# the varied partner and product use source-local identities instead of guessed SMILES.
scope = {
    "3b": (94, 93, "oxindole precursor displayed for 3b (fluoro-substituted)"),
    "3c": (90, 91, "oxindole precursor displayed for 3c (chloro-substituted)"),
    "3d": (85, 90, "oxindole precursor displayed for 3d (bromo-substituted)"),
    "3e": (80, 92, "oxindole precursor displayed for 3e (methyl-substituted)"),
    "3f": (95, 91, "oxindole precursor displayed for 3f (methoxy-substituted)"),
    "3g": (95, 94, "oxindole precursor displayed for 3g (chloro-substituted regioisomer)"),
    "3h": (85, 95, "oxindole precursor displayed for 3h (bromo-substituted regioisomer)"),
    "3i": (80, 81, "7-substituted oxindole precursor displayed for 3i"),
    "3j": (83, 76, "7-substituted oxindole precursor displayed for 3j"),
    "3k": (85, 67, "allenic ketone precursor displayed for 3k (ortho-methyl aryl)"),
    "3l": (95, 92, "allenic ketone precursor displayed for 3l (meta-methyl aryl)"),
    "3m": (95, 93, "allenic ketone precursor displayed for 3m (para-methyl aryl)"),
    "3n": (96, 94, "cyclohexyl-containing allenic ketone precursor displayed for 3n"),
    "3o": (94, 92, "allenic ketone precursor displayed for 3o (para-fluoro aryl)"),
    "3p": (95, 91, "allenic ketone precursor displayed for 3p (para-chloro aryl)"),
    "3q": (93, 91, "allenic ketone precursor displayed for 3q (para-bromo aryl)"),
}
for label, (_, _, partner_name) in scope.items():
    partner_kind = "oxindole" if label in {f"3{x}" for x in "bcdefghij"} else "allenic_ketone"
    add_molecule(
        f"{DOC}:mol:scope_{partner_kind}_{label}",
        partner_name,
        EV[2] + EV[3],
        name=partner_name,
        chirality="achiral",
        geometry="not_applicable",
        representation="not_encoded",
        identity_note="Full source drawing was inspected, but a definitive concrete SMILES was not asserted from the raster.",
    )
    add_molecule(
        f"{DOC}:mol:{label}",
        label,
        EV[2] + EV[3],
        name=f"spirocyclopentene oxindole {label}",
        chirality="unspecified",
        geometry="not_applicable",
        representation="not_encoded",
        atom_centers=[{"center": "spiro quaternary carbon", "configuration": "not assigned in main text"}],
        assignment="ee reported for the sample; major absolute configuration not printed in pages 1–5",
        identity_note="Paper-local product structure read from Scheme 2; exact stereochemical SMILES remains unencoded.",
    )


# Scheme 3 compounds and less safely encodable reagents.
for label in ["4", "5", "6", "7", "8"]:
    add_molecule(
        f"{DOC}:mol:{label}",
        label,
        EV[3],
        name=f"total-synthesis intermediate {label}",
        chirality="partially_specified",
        geometry="not_applicable",
        representation="not_encoded",
        atom_centers="stereochemical wedge/dash depiction is retained in the Scheme 3 drawing",
        assignment="drawn stereochemical relationship retained; no R/S nomenclature printed on pages 1–5",
        identity_note="Exact isomeric SMILES was not safely established from the main-text raster.",
    )
add_molecule(
    f"{DOC}:mol:9",
    "9",
    EV[3],
    name="(+)-debromoflustramine B",
    chirality="specified",
    geometry="not_applicable",
    representation="not_encoded",
    atom_centers="stereochemical wedge/dash depiction shown for compound 9",
    assignment="absolute configuration reported as confirmed by comparison of optical rotation with literature; numerical rotation not printed on pages 1–5",
    identity_note="Natural product 9; exact isomeric SMILES and numerical optical rotation require the associated supporting information or expert encoding.",
)
for key, name in [
    ("Pb_OAc4", "lead tetraacetate"),
    ("n_BuLi", "n-butyllithium"),
    ("iPrPPh3I", "i-PrPPh3I (label as printed; exact salt notation unresolved)"),
    ("AlH3_EtMe2N", "AlH3·EtMe2N"),
]:
    add_molecule(f"{DOC}:mol:{key}", key, EV[3], name=name)


def set_charge(p, key, raw, normalized, ev, method="original_extraction", derivation=None):
    if method == "original_extraction":
        p["charges"][key] = extracted(raw, normalized, ev)
    else:
        p["charges"][key] = inferred(normalized, ev, derivation, raw=raw, method=method)


def participant(pid, mid, role, ev, *, charges=None, reference=None, product_role=None):
    p = st.new_participant(pid, mid, role, ev)
    if product_role:
        p["product_role"] = extracted(product_role, product_role, ev)
    elif role != "product":
        p["product_role"] = not_applicable(ev)
    if charges:
        for key, value in charges.items():
            raw, normalized = value
            set_charge(p, key, raw, normalized, ev)
    if reference:
        p["equivalence_reference"] = extracted(reference, reference, ev)
    return p


def base_reaction(rid, ev, record_type, categories, location, name, transformation):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = extracted(categories, categories, ev)
    r["location"] = extracted(location, location, ev)
    r["reaction_name"] = extracted(name, name, ev)
    r["transformation_type"] = inferred(
        transformation,
        ev,
        "Transformation type normalized from the paper's reaction name and displayed net transformation.",
    )
    r["review"] = {
        "status": "needs_review",
        "method": "model",
        "reviewer": None,
        "checked_fields": [],
        "unresolved": [],
    }
    r["completeness"] = "partial"
    return r


def set_solvent(r, name, volume_raw, volume_value, ev):
    r["solvent_conditions"]["composition"] = extracted(name, [name], ev)
    r["solvent_conditions"]["mixture_ratio"] = missing(ev)
    r["solvent_conditions"]["total_volume"] = (
        extracted(volume_raw, volume_value, ev) if volume_raw else missing(ev)
    )
    r["solvent_conditions"]["concentration"] = missing(ev)
    r["solvent_conditions"]["concentration_basis"] = missing(ev)


def make_step(step_id, ev, participant_ids, operation, temperature=None, time=None):
    s = st.new_step(step_id, ev)
    s["participant_ids"] = list(participant_ids)
    s["operation"] = extracted(operation, operation, ev)
    s["temperature"] = (
        extracted(temperature[0], temperature[1], ev) if temperature else missing(ev)
    )
    s["time"] = extracted(time[0], time[1], ev) if time else missing(ev)
    s["atmosphere"] = missing(ev)
    s["pressure"] = missing(ev)
    return s


def yield_outcome(oid, target, raw, value, ev, *, isolated=False, scope=None):
    o = st.new_outcome(oid, "yield", [target], ev)
    o["product_role"] = extracted("product", "main_product", ev)
    o["value"] = extracted(raw, value, ev)
    o["yield_type"] = extracted("isolated product", "isolated", ev) if isolated else missing(ev)
    o["measurement_basis"] = extracted("isolated", "isolated", ev) if isolated else missing(ev)
    o["scope"] = extracted(scope, scope, ev) if scope else missing(ev)
    return o


def ee_outcome(oid, target, ee, ev):
    o = st.new_outcome(oid, "selectivity", [target], ev)
    o["product_role"] = extracted("product", "main_product", ev)
    if ee < 0:
        normalized = {"metric": "ee", "signed_value": ee, "unit_symbol": "%", "absolute_configuration": "unassigned"}
    else:
        normalized = quantity(ee, "%")
    o["value"] = extracted(f"{ee}% ee" if ee >= 0 else f"−{abs(ee)}% ee", normalized, ev)
    o["yield_type"] = not_applicable(ev)
    o["measurement_basis"] = extracted(
        "Determined by HPLC analysis on a chiral stationary phase",
        "chiral-stationary-phase HPLC",
        ev,
    )
    o["scope"] = extracted(
        "ee; major absolute configuration not assigned",
        {"metric": "ee", "major_enantiomer": "unassigned", "minor_enantiomer": "unassigned"},
        ev,
    )
    return o


def no_reaction_outcome(oid, raw, ev):
    o = st.new_outcome(oid, "no_reaction", [], ev)
    o["product_role"] = not_applicable(ev)
    o["value"] = extracted(raw, raw, ev)
    o["yield_type"] = not_applicable(ev)
    o["measurement_basis"] = missing(ev)
    o["scope"] = extracted("desired [4 + 1] annulation product", "desired product", ev)
    return o


def classification_no_reaction(r, ev):
    c = r["classification"]
    c["mapped_reaction_smiles"] = missing(ev)
    c["bond_changes"] = []
    c["labels"] = inferred([], ev, "No observed product was reported, so intended bonds are not classified.")
    c["coverage"] = inferred("not_applicable", ev, "No observed transformation was reported.")
    c["limitations"] = ["No intended-product bond labels are asserted for a no-reaction control."]
    c["model_votes"] = []


def classification_changes(r, ev, changes, limitations):
    c = r["classification"]
    c["mapped_reaction_smiles"] = missing(ev)
    c["bond_changes"] = [
        inferred(ch, ev, "Partial bond-change assignment from the displayed net transformation; exact atom mapping was not asserted.")
        for ch in changes
    ]
    labels = {}
    for ch in changes:
        if ch["change_type"] not in {"formed", "order_increased", "hydrogen_gain"}:
            continue
        code = st.bond_label(ch["atoms"][0]["element"], ch["atoms"][1]["element"], ch["after_order"])
        labels.setdefault(code["code"], {**code, "supporting_change_ids": []})["supporting_change_ids"].append(ch["change_id"])
    c["labels"] = inferred(list(labels.values()), ev, "Labels generated from the explicitly retained partial bond changes under bond-formation-20261016.")
    c["coverage"] = inferred("partial", ev, "Only unambiguous net bonds visible in the scheme are represented; no complete atom map is claimed.")
    c["limitations"] = list(limitations)
    c["model_votes"] = []


def cc_annulation_changes():
    return [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed", "source_description": "first ring-closing C-C connection"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed", "source_description": "second ring-closing C-C connection at the oxindole spiro carbon"},
    ]


def add_reaction(r):
    bundle["reactions"].append(r)


# Current-work generic Scheme 1.
r = base_reaction(
    "scheme1_current_generic", EV[2], "generic_scheme", ["other"],
    {"page": 2, "region": "Scheme 1, This work"},
    "phosphine-catalyzed enantioselective [4 + 1] annulation",
    "[4 + 1] annulation forming a spirocyclopentene oxindole",
)
r["participants"] = [
    participant("p1", f"{DOC}:mol:generic_allenic_ketone_1", "reactant", EV[2]),
    participant("p2", f"{DOC}:mol:generic_oxindole_2", "reactant", EV[2]),
    participant("p3", f"{DOC}:mol:generic_product_3", "product", EV[2], product_role="product"),
]
r["reaction_scale"] = missing(EV[2])
classification_changes(r, EV[2], cc_annulation_changes(), ["Generic R-group scheme; exact structures, hydrogens, and atom mapping are unresolved."])
r["review"]["unresolved"] = ["Generic R1/R2/PG variables are not expanded to concrete molecules."]
add_reaction(r)


def comparison_reaction(rid, allene_mid, product_mid, product_label, yield_raw, yield_value):
    ev = EV[2]
    r = base_reaction(
        rid, ev, "specific_experiment", ["control"],
        {"page": 2, "region": f"equation {1 if rid.endswith('eq1') else 2}"},
        "phosphine-catalyzed [4 + 1] annulation reactivity comparison",
        "[4 + 1] annulation forming a racemic spirocyclopentene oxindole",
    )
    r["general_scheme_ids"] = [f"{DOC}:scheme1_current_generic"]
    r["participants"] = [
        participant("allene", allene_mid, "reactant", ev, charges={"amount": ("0.12 mmol", quantity(0.12, "mmol"))}),
        participant("oxindole", f"{DOC}:mol:2a", "reactant", ev, charges={"amount": ("0.10 mmol", quantity(0.10, "mmol"))}),
        participant("PPh3", f"{DOC}:mol:PPh3", "catalyst", ev, charges={"mol_percent": ("20 mol%", quantity(20, "mol%"))}, reference=f"{DOC}:mol:2a"),
        participant("base", f"{DOC}:mol:Na2CO3", "reagent", ev),
        participant("solvent", f"{DOC}:mol:toluene", "solvent", ev, charges={"volume": ("2 mL", quantity(2, "mL"))}),
        participant("product", product_mid, "product", ev, product_role="product"),
    ]
    r["reaction_scale"] = extracted(
        "2a (0.10 mmol)",
        {"reference_molecule_id": f"{DOC}:mol:2a", "amount": quantity(0.10, "mmol")},
        ev,
    )
    set_solvent(r, "toluene", "2 mL", quantity(2, "mL"), ev)
    r["steps"] = [make_step("reaction", ev, ["allene", "oxindole", "PPh3", "base", "solvent"], "react at room temperature", ("r.t.", "r.t."), ("2 h", quantity(2, "h")))]
    r["procedure"]["original_text"] = extracted(
        f"{product_label}: PPh3 (20 mol%), Na2CO3, toluene (2 mL), r.t., 2 h; {yield_raw}",
        f"PPh3/Na2CO3/toluene comparison reaction to {product_label}",
        ev,
    )
    r["outcomes"] = [yield_outcome("yield", product_mid, yield_raw, yield_value, ev)]
    classification_changes(r, ev, cc_annulation_changes(), ["Product graph is not atom-mapped; classification is limited to the two annulation C-C closures."])
    r["review"]["unresolved"] = ["Yield measurement basis is not stated in the equation/main-text prose.", "Workup, purification, atmosphere, and exact base charge are unavailable without associated SI."]
    return r


add_reaction(comparison_reaction(
    "reactivity_eq1", f"{DOC}:mol:1a", f"{DOC}:mol:3a_rac", "3a-rac",
    ">99% yield", quantity(99, "%", ">99%", "gt"),
))
add_reaction(comparison_reaction(
    "reactivity_eq2", f"{DOC}:mol:1a_prime", f"{DOC}:mol:3a_prime_rac", "3a′-rac",
    "65% yield", quantity(65, "%"),
))


def standard_annulation_participants(allene_mid, oxindole_mid, product_mid, ev, *, product=True):
    items = [
        participant("allene", allene_mid, "reactant", ev, charges={"amount": ("0.12 mmol", quantity(0.12, "mmol"))}),
        participant("oxindole", oxindole_mid, "reactant", ev, charges={"amount": ("0.10 mmol", quantity(0.10, "mmol"))}),
        participant("P8", f"{DOC}:mol:P8", "catalyst", ev, charges={"amount": ("0.02 mmol", quantity(0.02, "mmol")), "mol_percent": ("20 mol%", quantity(20, "mol%"))}, reference=oxindole_mid),
        participant("base", f"{DOC}:mol:Na2CO3", "reagent", ev, charges={"amount": ("0.10 mmol", quantity(0.10, "mmol"))}),
        participant("solvent", f"{DOC}:mol:TBME", "solvent", ev, charges={"volume": ("2 mL", quantity(2, "mL"))}),
    ]
    if product:
        items.append(participant("product", product_mid, "product", ev, product_role="product"))
    return items


# Allenoate under optimized conditions: distinct from equation (2).
ev = EV[2] + EV[3]
r = base_reaction(
    "allenoate_optimized_no_reaction", ev, "specific_experiment", ["control"],
    {"page": 2, "region": "prose below Table 1 discussion; optimized conditions inherited from Table 1 entry 16"},
    "optimized-condition allenoate control",
    "no observed [4 + 1] annulation",
)
r["general_scheme_ids"] = [f"{DOC}:scheme1_current_generic"]
r["participants"] = standard_annulation_participants(f"{DOC}:mol:1a_prime", f"{DOC}:mol:2a", None, ev, product=False)
r["reaction_scale"] = extracted("2a (0.10 mmol)", {"reference_molecule_id": f"{DOC}:mol:2a", "amount": quantity(0.10, "mmol")}, ev)
set_solvent(r, "TBME", "2 mL", quantity(2, "mL"), ev)
r["steps"] = [make_step("reaction", ev, ["allene", "oxindole", "P8", "base", "solvent"], "apply optimized annulation conditions", ("−35 °C", quantity(-35, "°C")), ("12 h", quantity(12, "h")))]
r["procedure"]["general_procedure_reference"] = extracted("optimized reaction conditions", "Table 1 entry 16 conditions", ev)
r["outcomes"] = [no_reaction_outcome("no_reaction", "no reaction was observed", EV[2])]
classification_no_reaction(r, ev)
r["review"]["unresolved"] = ["Detection method and limit for the no-reaction observation are not reported.", "Conditions are explicitly described as optimized conditions; charges are inherited from Table 1 entry 16 rather than reprinted beside this control."]
add_reaction(r)


# Three N-protection failure controls under the Scheme 2 standard conditions.
for suffix, mid, label in [
    ("unprotected", f"{DOC}:mol:oxindole_NH", "unprotected oxindole"),
    ("n_methyl", f"{DOC}:mol:oxindole_NMe", "N-methyl oxindole"),
    ("n_benzyl", f"{DOC}:mol:oxindole_NBn", "N-benzyl oxindole"),
]:
    ev = EV[2] + EV[3]
    r = base_reaction(
        f"scope_failure_{suffix}", ev, "specific_experiment", ["substrate_scope", "control"],
        {"page": 2, "region": "substrate-scope discussion; standard conditions from Scheme 2 caption on page 3"},
        f"failed [4 + 1] annulation with {label}",
        "no desired [4 + 1] annulation product",
    )
    r["general_scheme_ids"] = [f"{DOC}:scheme1_current_generic"]
    r["participants"] = standard_annulation_participants(f"{DOC}:mol:1a", mid, None, ev, product=False)
    r["reaction_scale"] = extracted(label + " (0.10 mmol)", {"reference_molecule_id": mid, "amount": quantity(0.10, "mmol")}, ev)
    set_solvent(r, "TBME", "2 mL", quantity(2, "mL"), ev)
    r["steps"] = [make_step("reaction", ev, ["allene", "oxindole", "P8", "base", "solvent"], "perform Scheme 2 standard annulation", ("−35 °C", quantity(-35, "°C")), ("12 h", quantity(12, "h")))]
    r["procedure"]["general_procedure_reference"] = extracted("Scheme 2 standard conditions", "Scheme 2 standard conditions", ev)
    r["outcomes"] = [no_reaction_outcome("no_desired_product", "no desired [4 + 1] annulation products were observed", EV[2])]
    classification_no_reaction(r, ev)
    r["review"]["unresolved"] = ["No product identity, detection method, conversion, or recovered starting material is reported.", "The prose names the failed substrate class; standard charges/conditions are inherited from the Scheme 2 caption."]
    add_reaction(r)


# Table 1: nineteen distinct optimization rows.
table_rows = [
    (1, "P1", "CH2Cl2", 2, "r.t.", 99, -9),
    (2, "P2", "CH2Cl2", 2, "r.t.", 99, -11),
    (3, "P3", "CH2Cl2", 2, "r.t.", 99, -64),
    (4, "P4", "CH2Cl2", 2, "r.t.", 98, -53),
    (5, "P5", "CH2Cl2", 2, "r.t.", 98, 74),
    (6, "P6", "CH2Cl2", 2, "r.t.", 98, 71),
    (7, "P7", "CH2Cl2", 2, "r.t.", 98, 66),
    (8, "P8", "CH2Cl2", 2, "r.t.", 98, 75),
    (9, "P8", "THF", 2, "r.t.", 80, 84),
    (10, "P8", "toluene", 2, "r.t.", 99, 85),
    (11, "P8", "dioxane_1_4", 2, "r.t.", 85, 70),
    (12, "P8", "ether", 2, "r.t.", 98, 87),
    (13, "P8", "TBME", 2, "r.t.", 98, 88),
    (14, "P8", "TBME", 4, 0, 98, 88),
    (15, "P8", "TBME", 6, -25, 95, 90),
    (16, "P8", "TBME", 12, -35, 95, 92),
    (17, "P8", "TBME", 24, -40, 93, 92),
    (18, "P8", "TBME", 12, -35, 87, 92),
    (19, "P8", "TBME", 12, -35, 71, 92),
]
for entry, cat, solvent, hours, temp, yld, ee in table_rows:
    ev = EV[3]
    r = base_reaction(
        f"table1_entry_{entry:02d}", ev, "specific_experiment", ["optimization"] + (["scale_up"] if entry == 19 else []),
        {"page": 3, "region": "Table 1", "entry": entry},
        "phosphine-catalyzed enantioselective [4 + 1] annulation",
        "[4 + 1] annulation forming 3a",
    )
    r["general_scheme_ids"] = [f"{DOC}:scheme1_current_generic"]
    large = entry == 19
    allene_amt, ox_amt, base_amt, cat_amt, volume = (
        (4.8, 4.0, 4.0, 0.8, 20.0) if large else (0.12, 0.10, 0.10, 0.02, 2.0)
    )
    catalyst_charges = {"mol_percent": (("10 mol%" if entry == 18 else "20 mol%"), quantity(10 if entry == 18 else 20, "mol%"))}
    if entry != 18:
        catalyst_charges["amount"] = (f"{cat_amt:g} mmol", quantity(cat_amt, "mmol"))
    r["participants"] = [
        participant("1a", f"{DOC}:mol:1a", "reactant", ev, charges={"amount": (f"{allene_amt:g} mmol", quantity(allene_amt, "mmol"))}),
        participant("2a", f"{DOC}:mol:2a", "reactant", ev, charges={"amount": (f"{ox_amt:g} mmol", quantity(ox_amt, "mmol"))}),
        participant("catalyst", f"{DOC}:mol:{cat}", "catalyst", ev, charges=catalyst_charges, reference=f"{DOC}:mol:2a"),
        participant("base", f"{DOC}:mol:Na2CO3", "reagent", ev, charges={"amount": (f"{base_amt:g} mmol", quantity(base_amt, "mmol"))}),
        participant("solvent", f"{DOC}:mol:{solvent}", "solvent", ev, charges={"volume": (f"{volume:g} mL", quantity(volume, "mL"))}),
        participant("3a", f"{DOC}:mol:3a", "product", ev, product_role="product"),
    ]
    r["reaction_scale"] = extracted(f"2a ({ox_amt:g} mmol)", {"reference_molecule_id": f"{DOC}:mol:2a", "amount": quantity(ox_amt, "mmol")}, ev)
    solvent_name = {"CH2Cl2": "CH2Cl2", "THF": "THF", "toluene": "Toluene", "dioxane_1_4": "1,4-Dioxane", "ether": "Ether", "TBME": "TBME"}[solvent]
    set_solvent(r, solvent_name, f"{volume:g} mL", quantity(volume, "mL"), ev)
    temp_fact = ("r.t.", "r.t.") if temp == "r.t." else (f"{temp} °C", quantity(temp, "°C"))
    r["steps"] = [make_step("reaction", ev, ["1a", "2a", "catalyst", "base", "solvent"], "reaction-screening run", temp_fact, (f"{hours} h", quantity(hours, "h")))]
    r["procedure"]["original_text"] = extracted(f"Table 1 entry {entry}: {cat}, {solvent_name}, {hours} h, {temp if temp == 'r.t.' else str(temp)+' °C'}", f"Table 1 entry {entry}", ev)
    if entry == 18:
        r["procedure"]["modifications"] = extracted("Catalyst loading decreased to 10 mol%", "P8 loading 10 mol%", ev)
    if entry == 19:
        r["procedure"]["modifications"] = extracted("1a (4.8 mmol), 2a (4 mmol), Na2CO3 (4 mmol), catalyst (0.8 mmol), TBME (20 mL)", "gram-scale charges", ev)
    r["outcomes"] = [
        yield_outcome("isolated_yield", f"{DOC}:mol:3a", f"{yld}%", quantity(yld, "%"), ev, isolated=True),
        ee_outcome("ee", f"{DOC}:mol:3a", ee, ev),
    ]
    classification_changes(r, ev, cc_annulation_changes(), ["Exact product atom map and hydrogen transfers are unresolved; only the two new annulation C-C connections are classified."])
    r["review"]["unresolved"] = ["Workup, purification details, atmosphere, and absolute configuration of the major enantiomer require associated SI or expert review."]
    if entry == 18:
        r["review"]["unresolved"].append("The general 0.02 mmol catalyst footnote is superseded by the explicit 10 mol% entry-18 exception; an absolute catalyst amount is not invented.")
    add_reaction(r)


# Scheme 2: sixteen substrate-scope products.
for label, (yld, ee, partner_name) in scope.items():
    ev = EV[2] + EV[3]
    variant_is_oxindole = label in {f"3{x}" for x in "bcdefghij"}
    allene_mid = f"{DOC}:mol:1a" if variant_is_oxindole else f"{DOC}:mol:scope_allenic_ketone_{label}"
    oxindole_mid = f"{DOC}:mol:scope_oxindole_{label}" if variant_is_oxindole else f"{DOC}:mol:2a"
    product_mid = f"{DOC}:mol:{label}"
    r = base_reaction(
        f"scheme2_{label}", ev, "specific_experiment", ["substrate_scope"],
        {"page": 3, "region": "Scheme 2", "product": label},
        "phosphine-catalyzed enantioselective [4 + 1] annulation",
        f"[4 + 1] annulation forming {label}",
    )
    r["general_scheme_ids"] = [f"{DOC}:scheme1_current_generic"]
    r["participants"] = standard_annulation_participants(allene_mid, oxindole_mid, product_mid, ev)
    r["reaction_scale"] = extracted("2 (0.1 mmol)", {"reference_molecule_id": oxindole_mid, "amount": quantity(0.10, "mmol")}, ev)
    set_solvent(r, "TBME", "2 mL", quantity(2, "mL"), ev)
    r["steps"] = [make_step("reaction", ev, ["allene", "oxindole", "P8", "base", "solvent"], "Scheme 2 standard annulation", ("−35 °C", quantity(-35, "°C")), ("12 h", quantity(12, "h")))]
    r["procedure"]["general_procedure_reference"] = extracted("Scheme 2 caption conditions", "Scheme 2 standard conditions", EV[3])
    r["outcomes"] = [
        yield_outcome("isolated_yield", product_mid, f"{yld}% yield", quantity(yld, "%"), EV[3], isolated=True),
        ee_outcome("ee", product_mid, ee, EV[3]),
    ]
    classification_changes(r, ev, cc_annulation_changes(), ["Exact source structures are retained by labels/drawings, but atom mapping and complete hydrogen changes are unresolved."])
    r["review"]["unresolved"] = ["Exact definitive SMILES/atom mapping and major-enantiomer absolute configuration are not established from the main-text raster.", "Workup and purification details require explicitly associated SI."]
    add_reaction(r)


def add_scheme3_reaction(
    rid, start, product, name, participants_extra, steps, yld, classification_specs, limitations, unresolved, solvent=None, outcome_scope=None
):
    ev = EV[3]
    r = base_reaction(
        rid, ev, "specific_experiment", ["derivatization"],
        {"page": 3, "region": "Scheme 3", "arrow": f"{start} to {product}"},
        name,
        name,
    )
    r["participants"] = [
        participant("start", f"{DOC}:mol:{start}", "reactant", ev),
        *participants_extra,
        participant("product", f"{DOC}:mol:{product}", "product", ev, product_role="product"),
    ]
    r["reaction_scale"] = missing(ev)
    if solvent:
        set_solvent(r, solvent[0], solvent[1], solvent[2], ev)
    r["steps"] = steps
    r["procedure"]["original_text"] = extracted(name, name, ev)
    r["outcomes"] = [yield_outcome("yield", f"{DOC}:mol:{product}", f"{yld}% yield", quantity(yld, "%"), ev, scope=outcome_scope)]
    if classification_specs is None:
        r["classification"]["mapped_reaction_smiles"] = missing(ev)
        r["classification"]["bond_changes"] = []
        r["classification"]["labels"] = inferred([], ev, "No bond labels are asserted where the exact net graph change is unresolved from the raster.")
        r["classification"]["coverage"] = inferred("unresolved", ev, "The exact net bond changes could not be assigned safely from pages 1–5.")
        r["classification"]["limitations"] = limitations
        r["classification"]["model_votes"] = []
    else:
        classification_changes(r, ev, classification_specs, limitations)
    r["review"]["unresolved"] = unresolved + ["Yield measurement/isolation basis, exact scale, atmosphere, workup, and purification are not reported in the main text."]
    add_reaction(r)


# 3a -> 4: ordered ozonolysis and PPh3 workup.
add_scheme3_reaction(
    "scheme3_3a_to_4", "3a", "4", "ozonolysis of 3a followed by PPh3 treatment",
    [
        participant("ozone", f"{DOC}:mol:ozone", "reagent", EV[3]),
        participant("PPh3", f"{DOC}:mol:PPh3", "reagent", EV[3]),
        participant("solvent", f"{DOC}:mol:CH2Cl2", "solvent", EV[3]),
    ],
    [
        make_step("ozonolysis", EV[3], ["start", "ozone", "solvent"], "treat with O3 in CH2Cl2", ("−78 °C", quantity(-78, "°C"))),
        make_step("PPh3_treatment", EV[3], ["PPh3"], "then PPh3", None, None),
    ],
    90,
    [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}], "before_order": 0, "after_order": 2, "change_type": "formed", "source_description": "one carbonyl formed by alkene ozonolysis"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}], "before_order": 0, "after_order": 2, "change_type": "formed", "source_description": "second carbonyl formed by alkene ozonolysis"},
        {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 2, "after_order": 0, "change_type": "broken", "source_description": "cyclopentene C=C cleaved"},
    ],
    ["Carbonyl formation/alkene cleavage is visible, but exact atom maps and coproduct oxygen accounting are unresolved."],
    ["Duration and temperature of the PPh3 stage are unreported."],
    solvent=("CH2Cl2", None, None),
)


# 4 -> 5: exact oxidative rearrangement graph is not asserted.
add_scheme3_reaction(
    "scheme3_4_to_5", "4", "5", "lead tetraacetate oxidation of 4 in methanol",
    [
        participant("PbOAc4", f"{DOC}:mol:Pb_OAc4", "reagent", EV[3]),
        participant("solvent", f"{DOC}:mol:methanol", "solvent", EV[3]),
    ],
    [make_step("oxidation", EV[3], ["start", "PbOAc4", "solvent"], "Pb(OAc)4 in MeOH", None, ("1 h", quantity(1, "h")))],
    87, None,
    ["The raster supports oxidation to ester 5, but exact atom origins/bond changes are not safely mapped."],
    ["Reaction temperature is not printed beside the arrow."],
    solvent=("MeOH", None, None),
)


# 5 -> 6: two-stage deprotection/Wittig sequence; 70% is explicitly retained as two-step overall.
add_scheme3_reaction(
    "scheme3_5_to_6", "5", "6", "N-Boc cleavage followed by Wittig olefination",
    [
        participant("TFA", f"{DOC}:mol:TFA", "reagent", EV[3]),
        participant("CH2Cl2", f"{DOC}:mol:CH2Cl2", "solvent", EV[3]),
        participant("nBuLi", f"{DOC}:mol:n_BuLi", "reagent", EV[3]),
        participant("iPrPPh3I", f"{DOC}:mol:iPrPPh3I", "reagent", EV[3]),
    ],
    [
        make_step("deprotection", EV[3], ["start", "TFA", "CH2Cl2"], "TFA:CH2Cl2 = 1:10", ("0 °C", quantity(0, "°C"))),
        make_step("Wittig", EV[3], ["nBuLi", "iPrPPh3I"], "n-BuLi and i-PrPPh3I; warm from 0 °C to r.t.", ("0 °C to r.t.", {"number": 0, "unit": "°C", "end": "r.t.", "transition_time": None})),
    ],
    70,
    [{"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 2, "change_type": "formed", "source_description": "Wittig alkene formed from the aldehyde"}],
    ["Only the unambiguous new Wittig C=C is classified; Boc deprotection and other broken bonds are not atom-mapped."],
    ["Stage durations and the warm-up duration are unreported.", "The printed reagent label i-PrPPh3I is retained without silently repairing its salt notation."],
    solvent=("TFA:CH2Cl2", None, None), outcome_scope="two-step overall yield from 5 to 6",
)


# 6 -> 7 amidation.
add_scheme3_reaction(
    "scheme3_6_to_7", "6", "7", "amidation with methylamine",
    [
        participant("methylamine", f"{DOC}:mol:methylamine", "reagent", EV[3]),
        participant("solvent", f"{DOC}:mol:methanol", "solvent", EV[3]),
    ],
    [make_step("amidation", EV[3], ["start", "methylamine", "solvent"], "CH3NH2 in CH3OH", ("100 °C", quantity(100, "°C")), ("24 h", quantity(24, "h")))],
    95,
    [{"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "N"}], "before_order": 0, "after_order": 1, "change_type": "formed", "source_description": "methylamide C-N bond"}],
    ["C-N formation is classified; exact leaving-group and hydrogen transfers are unresolved."],
    [], solvent=("CH3OH", None, None),
)


# 7 -> 8 N-prenylation.
add_scheme3_reaction(
    "scheme3_7_to_8", "7", "8", "N-alkylation with 1-bromo-3-methylbut-2-ene",
    [
        participant("NaH", f"{DOC}:mol:NaH", "reagent", EV[3]),
        participant("prenyl_bromide", f"{DOC}:mol:prenyl_bromide", "reagent", EV[3]),
    ],
    [make_step("N_alkylation", EV[3], ["start", "NaH", "prenyl_bromide"], "NaH and 1-bromo-3-methylbut-2-ene; warm from 0 °C to r.t.", ("0 °C to r.t.", {"number": 0, "unit": "°C", "end": "r.t.", "transition_time": None}))],
    88,
    [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "N"}], "before_order": 0, "after_order": 1, "change_type": "formed", "source_description": "N-prenyl C-N bond"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "Br"}], "before_order": 1, "after_order": 0, "change_type": "broken", "source_description": "C-Br bond of prenyl bromide broken"},
    ],
    ["The new C-N bond and C-Br cleavage are represented without a complete atom map."],
    ["Solvent, reaction time, and warm-up duration are not printed beside the arrow."],
)


# 8 -> 9 reductive cyclization/reduction. C-H gains are retained qualitatively
# without claiming a complete hydrogen count.
add_scheme3_reaction(
    "scheme3_8_to_9", "8", "9", "AlH3·EtMe2N-induced reductive cyclization and further reduction",
    [participant("AlH3_complex", f"{DOC}:mol:AlH3_EtMe2N", "reagent", EV[3])],
    [
        make_step("reductive_cyclization", EV[3], ["start", "AlH3_complex"], "AlH3·EtMe2N", ("−15 °C", quantity(-15, "°C"))),
        make_step("further_reduction", EV[3], ["AlH3_complex"], "AlH3·EtMe2N", ("r.t.", "r.t.")),
    ],
    85,
    [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "N"}], "before_order": 0, "after_order": 1, "change_type": "formed", "source_description": "new ring-closing C-N connection in reductive cyclization"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": "unresolved", "source_description": "net C-H gain(s) from carbonyl reduction; exact count unresolved"},
    ],
    ["The source states reductive cyclization then further reduction; C-N and net C-H formation are retained, but exact atom mapping/hydrogen count and the unshown tricyclic lactam intermediate are unresolved."],
    ["Stage durations, reagent loadings, solvent, and the identity/yield of the unshown tricyclic lactam intermediate are not reported."],
)


# Non-reaction source context, bibliography metadata, and route-level results.
bundle["source_context"] = [
    {
        "context_id": f"{DOC}:context:bibliography",
        "kind": "bibliographic_metadata",
        "title": "Phosphine-catalyzed enantioselective [4 + 1] annulation of oxindoles with allenic ketones for the construction of spirocyclopentene oxindoles",
        "authors": ["Xiaodong Tang", "Huanzhen Ni", "Yixin Lu"],
        "journal": "Organic Chemistry Frontiers",
        "year": 2021,
        "volume": 8,
        "pages": "4485–4489",
        "doi": "10.1039/d1qo00669j",
        "evidence_ids": EV[1],
    },
    {
        "context_id": f"{DOC}:context:esi_notice",
        "kind": "supplementary_information_notice",
        "finding": "Page 1 states that ESI is available, but no SI file was supplied or explicitly associated in this workspace; no SI facts were imported.",
        "evidence_ids": EV[1],
    },
    {
        "context_id": f"{DOC}:context:figure1",
        "kind": "background_figure",
        "finding": "Figure 1 shows natural products and MI-5 as background structural examples, not experiments performed in this paper.",
        "evidence_ids": EV[1],
    },
    {
        "context_id": f"{DOC}:context:scheme1_prior_work",
        "kind": "literature_context",
        "finding": "The three upper Scheme 1 transformations are prior-work examples attributed to Marinetti (2010), Lu/Barbas (2011), and Chen (2016); they were not imported as current-paper experiments.",
        "evidence_ids": EV[2],
    },
    {
        "context_id": f"{DOC}:context:route_overall_yield",
        "kind": "route_level_outcome",
        "finding": "The sequence from 3a to (+)-debromoflustramine B (9) is reported in 38.9% overall yield; this was not reassigned to the final 8→9 arrow.",
        "value": {"number": 38.9, "unit": "%"},
        "start_molecule_id": f"{DOC}:mol:3a",
        "end_molecule_id": f"{DOC}:mol:9",
        "evidence_ids": EV[3],
    },
    {
        "context_id": f"{DOC}:context:optical_rotation",
        "kind": "characterization_assignment",
        "finding": "The absolute configuration of synthetic (+)-debromoflustramine B was reported as confirmed by comparison of optical rotation with literature; no numerical rotation is printed on pages 1–5.",
        "target_molecule_id": f"{DOC}:mol:9",
        "evidence_ids": EV[3],
    },
    {
        "context_id": f"{DOC}:context:references_pages4_5",
        "kind": "bibliography_only",
        "finding": "Pages 4–5 contain conclusion/administrative text and references; reaction names in cited titles were not treated as experiments of this paper.",
        "evidence_ids": EV[4] + EV[5],
    },
]


# Raw session artifacts: literal mention hints, inventory, stereo diagnostics, and coverage.
dataset = json.loads((RUN / "dataset.json").read_text())
mentions = st.source_mentions(dataset, ["1a", "2a", "3a", "3b", "3q", "debromoflustramine"], 180)
(RUN / "raw_model/mentions.json").write_text(json.dumps(mentions, ensure_ascii=False, indent=2) + "\n")

stereo_diagnostics = {}
for m in bundle["molecules"]:
    smi = st.value(m["smiles"])
    if smi:
        stereo_diagnostics[m["molecule_id"]] = st.stereo_inspect(smi)
(RUN / "raw_model/stereo_diagnostics.json").write_text(json.dumps(stereo_diagnostics, ensure_ascii=False, indent=2) + "\n")

inventory = {
    "document_id": DOC,
    "pages_inspected": [1, 2, 3, 4, 5],
    "planned_and_serialized_records": {
        "generic_current_work_scheme": 1,
        "page2_comparison_experiments": 2,
        "page2_no_reaction_controls": 4,
        "table1_optimization_rows": 19,
        "scheme2_scope_rows": 16,
        "scheme3_arrows": 6,
        "total": len(bundle["reactions"]),
    },
    "non_experiment_context": [
        "Figure 1 background structures",
        "three cited prior-work examples in Scheme 1",
        "ESI availability notice without supplied SI",
        "38.9% route-level overall yield",
        "optical-rotation configuration statement",
        "references on pages 4–5",
    ],
    "source_accuracy_verified": False,
}
(RUN / "raw_model/extraction_inventory.json").write_text(json.dumps(inventory, ensure_ascii=False, indent=2) + "\n")

page2_ids = [
    f"{DOC}:scheme1_current_generic",
    f"{DOC}:reactivity_eq1",
    f"{DOC}:reactivity_eq2",
    f"{DOC}:allenoate_optimized_no_reaction",
    f"{DOC}:scope_failure_unprotected",
    f"{DOC}:scope_failure_n_methyl",
    f"{DOC}:scope_failure_n_benzyl",
]
page3_ids = [f"{DOC}:table1_entry_{i:02d}" for i in range(1, 20)] + [
    f"{DOC}:scheme2_{label}" for label in scope
] + [
    f"{DOC}:scheme3_3a_to_4",
    f"{DOC}:scheme3_4_to_5",
    f"{DOC}:scheme3_5_to_6",
    f"{DOC}:scheme3_6_to_7",
    f"{DOC}:scheme3_7_to_8",
    f"{DOC}:scheme3_8_to_9",
]
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 1,
            "state": "attempted",
            "regions": [
                {"region": "title, abstract, and bibliographic block", "finding": "Current paper metadata and high-level [4 + 1] annulation claim; no specific experiment conditions.", "evidence_ids": EV[1]},
                {"region": "Figure 1", "finding": "Background natural products/MI-5 structures; not experiments of this paper.", "evidence_ids": EV[1]},
                {"region": "ESI footnote", "finding": "ESI stated available, but no SI file was supplied or explicitly associated.", "evidence_ids": EV[1]},
            ],
            "record_ids": [],
            "unresolved": ["No SI procedures or characterization were available for association."],
        },
        {
            "page": 2,
            "state": "attempted",
            "regions": [
                {"region": "Scheme 1 upper three transformations", "finding": "Prior-work [3 + 2] examples; retained as literature context, not current-paper experiments.", "evidence_ids": EV[2]},
                {"region": "Scheme 1 This work", "finding": "One generic current-work [4 + 1] annulation scheme.", "evidence_ids": EV[2]},
                {"region": "equations (1) and (2)", "finding": "Two distinct PPh3/Na2CO3/toluene comparison experiments yielding 3a-rac (>99%) and 3a′-rac (65%).", "evidence_ids": EV[2]},
                {"region": "optimization and scope prose", "finding": "Distinct optimized-condition allenoate no-reaction control and three N-protection substrate failures; Table 1/Scheme 2 conditions are on page 3.", "evidence_ids": EV[2]},
            ],
            "record_ids": page2_ids,
            "unresolved": ["Comparison-reaction yield measurement basis and no-reaction detection methods are not reported.", "Exact structures of generic R-group species are not concrete."],
        },
        {
            "page": 3,
            "state": "attempted",
            "regions": [
                {"region": "Table 1 including catalyst drawings and footnotes", "finding": "Nineteen distinct optimization experiments; entry 18 has 10 mol% catalyst exception and entry 19 is a 4 mmol scale run.", "evidence_ids": EV[3]},
                {"region": "Scheme 2 and caption", "finding": "Sixteen distinct scope experiments 3b–3q with isolated yields and chiral-HPLC ee under common conditions.", "evidence_ids": EV[3]},
                {"region": "Scheme 3", "finding": "Six distinct synthetic arrows from 3a to 9; 5→6 and 8→9 are multistage arrows, and 38.9% is route-level overall yield rather than the last-step yield.", "evidence_ids": EV[3]},
                {"region": "prose below Scheme 3", "finding": "Absolute configuration of 9 attributed to optical-rotation comparison; characterization context, not an extra synthesis.", "evidence_ids": EV[3]},
            ],
            "record_ids": page3_ids,
            "unresolved": ["Scheme 3 scales, most stage durations, atmospheres, workups, purification, and yield measurement bases are absent from the main text.", "Exact isomeric SMILES/atom maps for complex scope and route compounds are not safely established from the raster."],
        },
        {
            "page": 4,
            "state": "attempted",
            "regions": [
                {"region": "conclusion, conflicts, acknowledgements", "finding": "Narrative/administrative content; no new experimental arrow or table row.", "evidence_ids": EV[4]},
                {"region": "references 1–5", "finding": "Cited literature only; reaction titles were not imported as current-paper experiments.", "evidence_ids": EV[4]},
            ],
            "record_ids": [],
            "unresolved": [],
        },
        {
            "page": 5,
            "state": "attempted",
            "regions": [
                {"region": "references 5–15", "finding": "Bibliography only; no current-paper experiment.", "evidence_ids": EV[5]},
            ],
            "record_ids": [],
            "unresolved": [],
        },
    ],
}
(RUN / "coverage/pages_0001_0005.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")

(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]), "source_context": len(bundle["source_context"])}, ensure_ascii=False))
