#!/usr/bin/env python3
"""Build the source-derived v3.1 machine submission for doc_885b8c2ada399389.

This is an extraction artifact, not a tool or schema modification. Complex structures
remain unresolved when the supplied main-text drawings do not support safe encoding.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st  # noqa: E402

DOC = "doc_885b8c2ada399389"
OUT = ROOT / "runs/current"
EV = {p: f"{DOC}:page{p}:image" for p in range(1, 7)}
TXT = {p: f"{DOC}:page{p}:text" for p in range(1, 7)}
SENTINEL = object()


def F(raw=None, normalized=SENTINEL, evidence=(), method="original_extraction", derivation=None,
      status=None):
    if normalized is SENTINEL:
        normalized = raw
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation, status=status)


def derived(normalized, evidence, derivation):
    return F(None, normalized, evidence, method="inference", derivation=derivation)


bundle = json.loads((OUT / "submitted.json").read_text())
molecules = {}


def molecule(mid, label, evidence, name=None, smiles=None, chirality="unresolved",
             geometry="unresolved", representation="unresolved", atom_centers=None,
             double_bonds=None, axial=None, assignment=None, mixture=None, identity_note=None):
    m = st.new_molecule(f"{DOC}:mol:{mid}", DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, evidence)
    m["names"]["original"] = F(name or label, name or label, evidence)
    m["names"]["english"] = F(name or label, name or label, evidence)
    m["names"]["abbreviations"] = F(label, [label], evidence)
    m["identity_note"] = F(identity_note, identity_note, evidence) if identity_note else F()
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = F(
            name or label, smiles, evidence, method="structure_derivation",
            derivation="Candidate SMILES derived from the explicitly named/simple drawn species; parsing is not source-structure certification.",
        )
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = F()
    stereo = m["stereochemistry"]
    stereo["chirality_status"] = derived(chirality, evidence, "Status follows the printed stereochemical label/drawing; absence of a label is retained as unspecified or unresolved, never racemic.")
    stereo["geometry_status"] = derived(geometry, evidence, "Status follows printed E/Z labels and the main-text drawing.")
    stereo["representation_status"] = derived(representation, evidence, "Complex source structures were not encoded as definitive isomeric SMILES.")
    if atom_centers is not None:
        stereo["atom_centers"] = F(atom_centers, atom_centers, evidence)
    if double_bonds is not None:
        stereo["double_bonds"] = F(double_bonds, double_bonds, evidence)
    if axial is not None:
        stereo["axial"] = F(axial, axial, evidence)
    if assignment is not None:
        stereo["assignment"] = F(assignment, assignment, evidence)
    if mixture is not None:
        stereo["mixture"] = F(mixture, mixture, evidence)
    molecules[mid] = m
    return m["molecule_id"]


def generic_molecule(mid, label, evidence, fixed, variable, chirality, geometry):
    molecule(mid, label, evidence, name=label, chirality=chirality, geometry=geometry,
             representation="not_encoded", identity_note="Generic R-group structure from Scheme 1c.")
    m = molecules[mid]
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": F("Scheme 1c fixed skeleton", fixed, evidence),
        "variables": [{
            "label": "R",
            "attachment_points": F("R position shown in Scheme 1c", [variable], evidence),
            "definition": F("R = aryl, heteroaryl, alkenyl, and 1°, 2°, 3° alkyl group",
                            ["aryl", "heteroaryl", "alkenyl", "primary alkyl", "secondary alkyl", "tertiary alkyl"], evidence),
            "allowed_values": F("R = aryl, heteroaryl, alkenyl, and 1°, 2°, 3° alkyl group",
                                ["aryl", "heteroaryl", "alkenyl", "1° alkyl", "2° alkyl", "3° alkyl"], evidence),
            "constraints": F(),
        }],
        "instance_mapping": [],
        "reason": F("Scheme 1c is an R-group scheme, not one concrete molecule.",
                    "A generic template is required because the source leaves R variable.", evidence),
    }
    if geometry == "specified":
        m["stereochemistry"]["double_bonds"] = F(
            "E-selective", [{"assignment": "E", "site": "conjugated enone in Scheme 1c"}], evidence)
    return m["molecule_id"]


# Generic transformation and the complete Table 3 family.
generic_molecule("generic_1", "generic 2,3-allenol 1", [EV[1]],
                 {"description": "2,3-allenol skeleton drawn in Scheme 1c", "variable": "R"}, "R on carbinol carbon",
                 "unspecified", "not_applicable")
generic_molecule("generic_E2", "generic (E)-2", [EV[1]],
                 {"description": "(2,5-dihydrofuran-3-yl)but-2(E)-en-1-one skeleton", "variable": "R"}, "R substituents shown on both fragments",
                 "unspecified", "specified")

scope = [
    ("a", "Ph", 12.5, 62), ("b", "4-FC6H4", 11, 57), ("c", "2-ClC6H4", 7.5, 53),
    ("d", "3-ClC6H4", 10.5, 54), ("e", "3-BrC6H4", 11, 53), ("f", "4-BrC6H4", 12, 52),
    ("g", "4-MeOC6H4", 11, 57), ("h", "benzodioxol-5-yl", 12.5, 44),
    ("i", "2-naphthyl", 14, 52), ("j", "2-furyl", 12, 43), ("k", "2-thienyl", 11, 49),
    ("l", "cinnamyl", 13, 41), ("m", "n-C6H13", 18, 64), ("n", "n-C9H19", 13, 62),
    ("o", "n-C11H23", 12, 63), ("p", "i-Bu", 9, 51), ("q", "Bn", 10, 54),
    ("r", "(CH2)2Ph", 11, 62), ("s", "(CH2)5Cl", 12, 62), ("t", "(CH2)5OH", 12, 56),
    ("u", "Cy", 12.5, 60), ("v", "t-Bu", 10, 50),
]
for letter, substituent, _, _ in scope:
    input_label = f"(E)-1{letter}" if letter == "l" else f"1{letter}"
    input_geometry = "specified" if letter == "l" else "not_applicable"
    input_db = [{"assignment": "E", "site": "cinnamyl double bond"}] if letter == "l" else None
    molecule(f"1{letter}", input_label, [EV[2]], f"2,3-allenol {input_label}; R={substituent}",
             chirality="unspecified", geometry=input_geometry, representation="not_encoded", double_bonds=input_db,
             identity_note="Paper-local structure retained without guessed concrete SMILES.")
    product_label = f"(E,E,E)-2{letter}" if letter == "l" else f"(E)-2{letter}"
    db = [{"assignment": "E", "site": "conjugated enone"}]
    if letter == "l":
        db = [{"assignment": "E,E,E", "site": "three printed acyclic alkene assignments"}]
    molecule(f"E2{letter}", product_label, [EV[2]], f"{product_label}; R={substituent}",
             chirality="unspecified", geometry="specified", representation="not_encoded", double_bonds=db,
             identity_note="Printed E geometry retained; stereogenic ring center is not assigned in Table 3.")

# Simple reagents/solvents explicitly present in the main text.
simple = {
    "rhcat": ("[Cp*RhCl2]2", None), "cuoac": ("Cu(OAc)2·H2O", None),
    "methanol": ("MeOH", "CO"), "acetonitrile": ("CH3CN", "CC#N"), "oxygen": ("O2", "O=O"),
    "dce": ("DCE", "ClCCCl"), "dmf": ("DMF", "CN(C)C=O"), "dma": ("DMA", "CC(=O)N(C)C"),
    "dmso": ("DMSO", "CS(C)=O"), "dioxane": ("dioxane", "O1CCOCC1"),
    "water": ("H2O", "O"), "hoac": ("HOAc", "CC(=O)O"), "net3": ("NEt3", "CCN(CC)CC"),
    "mesitylene": ("mesitylene", "Cc1cc(C)cc(C)c1"),
    "toluene": ("toluene", "Cc1ccccc1"), "phenylacetylene": ("phenylacetylene", "C#Cc1ccccc1"),
    "pd_sono": ("Pd(PPh3)2Cl2", None), "ipr2nh": ("i-Pr2NH", "CC(C)NC(C)C"), "cui": ("CuI", "[Cu]I"),
    "thf": ("THF", "C1CCOC1"), "boronic": ("(4-methoxyphenyl)boronic acid", "COc1ccc(B(O)O)cc1"),
    "pd_suzuki": ("Pd(PPh3)4", None), "k3po4": ("K3PO4", None), "petasis": ("Cp2TiMe2", None),
    "cd3od": ("CD3OD", "[2H]C([2H])([2H])O"),
}
for mid, (name, smi) in simple.items():
    pages = [EV[1], EV[2], EV[3], EV[4]]
    molecule(mid, name, pages, name=name, smiles=smi, chirality="achiral" if smi else "not_applicable",
             geometry="not_applicable", representation="not_applicable",
             identity_note="Simple named/drawn reaction component; no CAS was supplied.")

# Chiral Table 4 compounds are distinct registry identities from unspecified Table 3 compounds.
chiral_rows = [
    ("S1a", "(S)-1a", "S", 98, "SE2a", "(S,E)-2a", 56, 98, 12.5),
    ("S1b", "(S)-1b", "S", ">99", "SE2b", "(S,E)-2b", 61, ">99", 12),
    ("S1e", "(S)-1e", "S", 95, "SE2e", "(S,E)-2e", 48, 93, 12),
    ("S1j", "(S)-1j", "S", 96, "SE2j", "(S,E)-2j", 42, 95, 12),
    ("R1m", "(R)-1m", "R", 98, "RE2m", "(R,E)-2m", 61, 98, 12),
    ("S1n", "(S)-1n", "S", 96, "SE2n", "(S,E)-2n", 67, 98, 14.5),
]
for imid, ilabel, config, _, pmid, plabel, _, _, _ in chiral_rows:
    molecule(imid, ilabel, [EV[3]], ilabel, chirality="specified", geometry="not_applicable",
             representation="not_encoded", atom_centers=[{"assignment": config, "site": "carbinol stereocenter"}],
             assignment={"absolute_configuration": config, "printed_label": ilabel})
    molecule(pmid, plabel, [EV[3]], plabel, chirality="specified", geometry="specified",
             representation="not_encoded", atom_centers=[{"assignment": config, "site": "dihydrofuran stereocenter"}],
             double_bonds=[{"assignment": "E", "site": "enone"}],
             assignment={"absolute_configuration": config, "alkene_geometry": "E", "printed_label": plabel})

# Cross-coupling, scale-up/derivatization, byproducts, and mechanistic-study species.
extra_molecules = [
    ("1x", "1x", "tertiary 2,3-allenol 1x", "unspecified", "not_applicable"),
    ("E2xe", "(E)-2xe", "cross-coupling product (E)-2xe", "unspecified", "specified"),
    ("E8f", "(E)-8f", "byproduct (E)-8f", "unspecified", "specified"),
    ("piperonal", "piperonyl aldehyde", "piperonyl aldehyde", "achiral", "not_applicable"),
    ("other_byproducts", "mixture of other byproducts", "mixture of other byproducts", "mixture", "unresolved"),
    ("S1f", "(S)-1f", "(S)-1f", "specified", "not_applicable"),
    ("SE2f", "(S,E)-2f", "(S,E)-2f", "specified", "specified"),
    ("SE3f", "(S,E)-3f", "conjugated triene (S,E)-3f", "specified", "specified"),
    ("SE4f", "(S,E)-4f", "Sonogashira product (S,E)-4f", "specified", "specified"),
    ("SE5f", "(S,E)-5f", "mono-Suzuki product (S,E)-5f", "specified", "specified"),
    ("SE6f", "(S,E)-6f", "double-Suzuki product (S,E)-6f", "specified", "specified"),
    ("1w", "1w", "2-heptyl-2,3-butadienol 1w", "unspecified", "not_applicable"),
    ("E2w", "(E)-2w", "expected dimeric product (E)-2w", "unspecified", "specified"),
    ("7w", "7w", "3-heptyl-2,5-dihydrofuran 7w", "unspecified", "not_applicable"),
    ("8e", "8e", "allenyl ketone 8e", "not_applicable", "unspecified"),
    ("9a", "9a", "alcohol 9a (dr = 1:1)", "mixture", "unspecified"),
    ("Da1a", "[Dᵃ]-1a", "deuterated allenol [Dᵃ]-1a", "unspecified", "not_applicable"),
    ("Db1m", "[Dᵇ]-1m", "deuterated allenol [Dᵇ]-1m", "unspecified", "not_applicable"),
    ("DE2m", "[D]-(E)-2m", "deuterated product [D]-(E)-2m", "unspecified", "specified"),
]
for mid, label, name, chi, geo in extra_molecules:
    ev = [EV[2]] if mid in {"E8f", "piperonal", "other_byproducts"} else [EV[3]] if mid in {"1x", "E2xe", "S1f", "SE2f", "SE3f", "SE4f", "SE5f", "SE6f"} else [EV[4]]
    atom = None
    assignment = None
    if mid in {"S1f", "SE2f", "SE3f", "SE4f", "SE5f", "SE6f"}:
        atom = [{"assignment": "S", "site": "printed stereocenter"}]
        assignment = {"absolute_configuration": "S", "printed_label": label}
    db = [{"assignment": "E", "site": "printed acyclic alkene"}] if geo == "specified" else None
    mix = {"description": "dr = 1:1"} if mid == "9a" else ({"description": "unknown-composition byproduct mixture"} if mid == "other_byproducts" else None)
    molecule(mid, label, ev, name, chirality=chi, geometry=geo, representation="not_encoded",
             atom_centers=atom, double_bonds=db, assignment=assignment, mixture=mix,
             identity_note="Paper-local graphical identity; exact SMILES was not guessed from the main-text raster.")


def mid(key):
    return molecules[key]["molecule_id"]


def set_charge(p, key, raw, normalized, ev):
    p["charges"][key] = F(raw, normalized, ev)


def add_part(r, key, role, ev, *, stage="reaction", amount=None, mass=None, volume=None,
             equiv=None, molpct=None, ref=None, product_role=None, pid=None):
    pid = pid or f"p{len(r['participants'])+1}_{key}"
    p = st.new_participant(pid, mid(key), role, ev, stage=stage)
    if product_role:
        p["product_role"] = F(product_role, product_role, ev)
    for field, val in (("amount", amount), ("mass", mass), ("volume", volume),
                       ("equivalents", equiv), ("mol_percent", molpct)):
        if val is not None:
            set_charge(p, field, val[0], val[1], ev)
    if ref is not None:
        p["equivalence_reference"] = F(ref, mid(ref), ev)
    r["participants"].append(p)
    return pid


def reaction(rid, ev, record_type="specific_experiment", categories=("other",), location=None,
             name="Rh-catalyzed oxidative homo-coupling cyclization of 2,3-allenols",
             transformation="oxidative dimerization/cyclization"):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = F(list(categories), list(categories), ev)
    r["location"] = F(location or rid, location or rid, ev)
    r["reaction_name"] = F(name, name, ev)
    r["transformation_type"] = F(None, transformation, ev, method="inference",
                                  derivation="Transformation class inferred from the source title, scheme, and depicted net change.")
    r["review"]["unresolved"] = [
        "No explicitly associated SI was supplied; workup, purification, and unreported procedural detail remain unknown.",
        "Complex primary structures and atom mapping require expert encoding from the source drawings.",
    ]
    return r


def set_scale(r, key, raw, number, ev):
    r["reaction_scale"] = F(raw, {"reference_molecule_id": mid(key), "amount": {"number": number, "unit": "mmol"}}, ev)


def set_solvents(r, keys, ev, ratio=None, total=None, concentration=None, basis=None):
    r["solvent_conditions"]["composition"] = F([molecules[k]["names"]["original"]["raw_value"] for k in keys], [mid(k) for k in keys], ev)
    if ratio is not None:
        r["solvent_conditions"]["mixture_ratio"] = F(ratio, ratio, ev)
    if total is not None:
        r["solvent_conditions"]["total_volume"] = F(total[0], total[1], ev)
    if concentration is not None:
        r["solvent_conditions"]["concentration"] = F(concentration[0], concentration[1], ev)
    if basis is not None:
        r["solvent_conditions"]["concentration_basis"] = F(basis, basis, ev)


def add_step(r, ev, operation, temp=None, time=None, atmosphere=None, pids=None, stage="reaction"):
    s = st.new_step(f"step_{len(r['steps'])+1}", ev, stage=stage)
    s["participant_ids"] = list(pids if pids is not None else [p["participant_id"] for p in r["participants"] if p["stage"] == stage and p["role"] != "product"])
    s["operation"] = F(operation, operation, ev)
    if temp is not None:
        s["temperature"] = F(temp[0], temp[1], ev)
    if time is not None:
        s["time"] = F(time[0], time[1], ev)
    if atmosphere is not None:
        s["atmosphere"] = F(atmosphere, atmosphere, ev)
    r["steps"].append(s)
    return s


def add_outcome(r, oid, kind, targets, ev, raw, normalized, *, yield_type=None, basis=None, scope=None,
                basis_reported=True):
    o = st.new_outcome(oid, kind, [mid(t) for t in targets], ev)
    o["value"] = F(raw, normalized, ev)
    if yield_type is not None:
        o["yield_type"] = F(yield_type, yield_type, ev)
    if basis is not None:
        o["measurement_basis"] = F(basis, basis, ev)
    elif basis_reported is False:
        o["measurement_basis"] = F()
    if scope is not None:
        o["scope"] = F(scope, scope, ev)
    o["product_role"] = F("target product", "target product", ev) if kind in {"yield", "product_amount"} else F()
    r["outcomes"].append(o)
    return o


TAX = {row["code"]: row["number"] for row in json.loads((ROOT / "skills/chemical-reaction-extraction/references/bond_taxonomy.json").read_text())["categories"]}


def set_class(r, codes, ev, rationale):
    changes = []
    supports = {}
    for i, code in enumerate(codes, 1):
        cid = f"b{i}"
        if code == "C-C": atoms, before, after, typ = [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], 0, 1, "formed"
        elif code == "C-O": atoms, before, after, typ = [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}], 0, 1, "formed"
        elif code == "DB-C-O": atoms, before, after, typ = [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}], 1, 2, "order_increased"
        elif code == "DB-C-C": atoms, before, after, typ = [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], 0, 2, "formed"
        else:
            raise ValueError(code)
        change = {"change_id": cid, "atoms": atoms, "before_order": before, "after_order": after, "change_type": typ}
        changes.append(F(None, change, ev, method="structure_derivation", derivation=rationale))
        supports.setdefault(code, []).append(cid)
    labels = [{"category_number": TAX[c], "code": c, "supporting_change_ids": ids} for c, ids in supports.items()]
    r["classification"]["bond_changes"] = changes
    r["classification"]["labels"] = F(None, labels, ev, method="inference", derivation=rationale)
    r["classification"]["coverage"] = F(None, "partial", ev, method="inference", derivation=rationale)
    r["classification"]["limitations"] = [
        "Partial graphical net-change classification only; no atom-mapped reaction is asserted.",
        "Hydrogen gains/losses and all bond-order changes cannot be completed safely from the main-text raster alone.",
    ]


def set_not_applicable(r, ev, reason):
    r["classification"]["coverage"] = F(None, "not_applicable", ev, method="inference", derivation=reason)
    r["classification"]["labels"] = F(None, [], ev, method="inference", derivation=reason)
    r["classification"]["limitations"] = [reason]


def standard_dimer_participants(r, substrate, product, ev, substrate_amount=None, catalyst=5, cu_equiv=2.5,
                                solvents=("acetonitrile", "methanol"), solvent_volumes=None, oxygen=True,
                                catalyst_amount=None, cu_amount=None):
    amount = (f"{substrate_amount} mmol", {"number": substrate_amount, "unit": "mmol"}) if substrate_amount is not None else None
    s = add_part(r, substrate, "reactant", ev, amount=amount, equiv=("1.0 equiv.", {"number": 1.0, "unit": "equiv"}), ref=substrate)
    add_part(r, product, "product", ev, product_role="target product")
    cat_amount = (f"{catalyst_amount} mmol", {"number": catalyst_amount, "unit": "mmol"}) if catalyst_amount is not None else None
    copper_amount = (f"{cu_amount} mmol", {"number": cu_amount, "unit": "mmol"}) if cu_amount is not None else None
    add_part(r, "rhcat", "catalyst", ev, amount=cat_amount,
             molpct=(f"{catalyst} mol%", {"number": catalyst, "unit": "mol%"}), ref=substrate)
    add_part(r, "cuoac", "reagent", ev, amount=copper_amount,
             equiv=(f"{cu_equiv} equiv.", {"number": cu_equiv, "unit": "equiv"}), ref=substrate)
    for idx, key in enumerate(solvents):
        vol = None
        if solvent_volumes and key in solvent_volumes:
            val = solvent_volumes[key]
            vol = (f"{val[0]} {val[1]}", {"number": val[0], "unit": val[1]})
        add_part(r, key, "solvent", ev, volume=vol)
    if oxygen:
        add_part(r, "oxygen", "reagent", ev)
    return s


reactions = []

# Scheme 1c: generic current-work transformation; panels 1a-b are source_context below.
r = reaction("scheme1c_generic", [EV[1]], record_type="generic_scheme", categories=("other",), location="Scheme 1c",
             transformation="Rh(III)-catalyzed oxidative homo-coupling cyclization")
standard_dimer_participants(r, "generic_1", "generic_E2", [EV[1]], catalyst=5, cu_equiv=2.5)
set_solvents(r, ["acetonitrile", "methanol"], [EV[1]], ratio="20/1")
add_step(r, [EV[1]], "generic reaction arrow", temp=("50 °C", {"number": 50, "unit": "°C"}), atmosphere="O2 (balloon)")
set_class(r, ["C-C", "C-O", "DB-C-O"], [EV[1]], "Partial bond changes inferred from the generic substrate/product drawings and oxidative coupling/cyclization description.")
reactions.append(r)

# Table 1, eight solvent-screen experiments.
table1 = [
    (1, "methanol", 14, 0, 0, False), (2, "methanol", 36, 15, 0, True),
    (3, "dce", 72, 15, "Trace", True), (4, "dmf", 72, 0, 33, True),
    (5, "dma", 72, 0, 25, True), (6, "dmso", 72, 2, 53, True),
    (7, "dioxane", 72, 7, 18, True), (8, "acetonitrile", 72, 20, 23, True),
]
for entry, solvent, hours, yld, rec, with_cu in table1:
    ev = [EV[2], TXT[1]] if entry in (1, 2, 8) else [EV[2]]
    r = reaction(f"table1_entry{entry:02d}", ev, categories=("optimization",), location=f"Table 1 entry {entry}")
    s = add_part(r, "1a", "reactant", ev, amount=("0.3 mmol", {"number": 0.3, "unit": "mmol"}), ref=None)
    r["participants"][-1]["charges"]["equivalents"] = F("1.0 equiv.", {"number": 1.0, "unit": "equiv"}, ev)
    r["participants"][-1]["equivalence_reference"] = F("1a", mid("1a"), ev)
    add_part(r, "E2a", "product", ev, product_role="target product")
    add_part(r, "rhcat", "catalyst", ev, molpct=("2.5 mol%", {"number": 2.5, "unit": "mol%"}), ref="1a")
    if with_cu:
        add_part(r, "cuoac", "reagent", ev, amount=("0.6 mmol", {"number": 0.6, "unit": "mmol"}),
                 equiv=("2.0 equiv.", {"number": 2.0, "unit": "equiv"}), ref="1a")
    add_part(r, solvent, "solvent", ev, volume=("1.5 mL", {"number": 1.5, "unit": "mL"}))
    add_part(r, "mesitylene", "additive", ev, stage="characterization")
    set_scale(r, "1a", "1a (0.3 mmol)", 0.3, ev)
    set_solvents(r, [solvent], ev, total=("1.5 mL", {"number": 1.5, "unit": "mL"}),
                 concentration=("c = 0.2 M", {"number": 0.2, "unit": "M"}), basis="1a")
    add_step(r, ev, "reaction under air balloon", temp=("50 °C", {"number": 50, "unit": "°C"}),
             time=(f"{hours} h", {"number": hours, "unit": "h"}), atmosphere="air (balloon)")
    add_outcome(r, "yield_E2a", "yield", ["E2a"], ev, f"{yld}%", {"number": yld, "unit": "%"},
                yield_type="1H NMR yield", basis="1H NMR analysis of crude products with mesitylene as internal standard")
    rec_norm = rec if isinstance(rec, str) else {"number": rec, "unit": "%"}
    rec_raw = rec if isinstance(rec, str) else f"{rec}%"
    add_outcome(r, "recovery_1a", "recovery", ["1a"], ev, rec_raw, rec_norm,
                basis="1H NMR analysis of crude products with mesitylene as internal standard")
    if entry == 1:
        add_outcome(r, "decomposition_1a", "decomposition", ["1a"], ev, "decomposition of 1a", True, basis_reported=False)
        r["notes"].append(F("A1, B1, and D1 were not detected; 1a decomposed.",
                            "Expected A1/B1/D1 not detected; decomposition observed.", ev))
    if yld > 0:
        set_class(r, ["C-C", "C-O", "DB-C-O"], ev, "Partial graphical classification of the observed dimerization product (E)-2a.")
    else:
        set_not_applicable(r, ev, "No (E)-2a formation was observed in this optimization row; no intended-product bond labels are assigned.")
    reactions.append(r)

# Table 2 optimization under O2. Entry 16 is also Table 3 entry 1 and is represented once.
table2 = [
    (1, 2.5, ["acetonitrile"], None, 60, 30, 9, 2.0, 0.3, None, None),
    (2, 4, ["acetonitrile"], None, 43, 38, 4, 2.0, 0.3, None, None),
    (3, 4, ["acetonitrile", "water"], "20/1", 36, 20, 3, 2.0, 0.3, None, None),
    (4, 4, ["acetonitrile", "hoac"], "20/1", 36, 40, 0, 2.0, 0.3, None, None),
    (5, 4, ["acetonitrile", "net3"], "20/1", 36, 0, 22, 2.0, 0.3, None, None),
    (6, 4, ["acetonitrile", "methanol"], "20/1", 36, 47, 8, 2.0, 0.3, None, None),
    (7, 4, ["acetonitrile", "methanol"], "10/1", 41, 40, 14, 2.0, 0.3, None, None),
    (8, 4, ["acetonitrile", "methanol"], "50/1", 41, 42, 9, 2.0, 0.3, None, None),
    (9, 5, ["acetonitrile", "methanol"], "20/1", 43, 49, 4, 2.0, 0.3, None, None),
    (10, 6, ["acetonitrile", "methanol"], "20/1", 43, 47, 0, 2.0, 0.3, None, None),
    (11, 7, ["acetonitrile", "methanol"], "20/1", 43, 38, 0, 2.0, 0.3, None, None),
    (12, 5, ["acetonitrile", "methanol"], "20/1", 48, 43, 3, 1.5, 0.3, None, None),
    (13, 5, ["acetonitrile", "methanol"], "20/1", 38, 49, 0, 2.5, 0.3, None, None),
    (14, 5, ["acetonitrile", "methanol"], "20/1", 60, 48, 2, 2.5, 1.0, 46, 0.2),
    (15, 5, ["acetonitrile", "methanol"], "20/1", 38, 58, 0, 2.5, 1.0, 58, 0.48),
    (16, 5, ["acetonitrile", "methanol"], "20/1", 12.5, 62, 0, 2.5, 1.0, 62, 0.95),
    (17, 5, ["acetonitrile", "methanol"], "20/1", 8, 60, 0, 2.5, 1.0, None, 1.9),
]
for entry, x, solvents, ratio, hours, yld, rec, cu_eq, scale, isolated, conc in table2:
    ev = [EV[2], TXT[1]] if entry <= 8 else [EV[2]]
    cats = ["optimization"]
    location = f"Table 2 entry {entry}"
    if entry == 16:
        cats.append("substrate_scope")
        location += "; duplicated as Table 3 entry 1"
    r = reaction(f"table2_entry{entry:02d}", ev, categories=tuple(cats), location=location)
    add_part(r, "1a", "reactant", ev, amount=(f"{scale} mmol", {"number": scale, "unit": "mmol"}),
             equiv=("1.0 equiv.", {"number": 1.0, "unit": "equiv"}), ref="1a")
    add_part(r, "E2a", "product", ev, product_role="target product")
    add_part(r, "rhcat", "catalyst", ev, molpct=(f"{x} mol%", {"number": x, "unit": "mol%"}), ref="1a")
    cu_amount = ("0.6 mmol", {"number": 0.6, "unit": "mmol"}) if entry <= 11 else (("2.5 mmol", {"number": 2.5, "unit": "mmol"}) if entry == 16 else None)
    add_part(r, "cuoac", "reagent", ev, amount=cu_amount,
             equiv=(f"{cu_eq} equiv.", {"number": cu_eq, "unit": "equiv"}), ref="1a")
    for key in solvents:
        solvent_volume = None
        if entry == 16 and key == "acetonitrile":
            solvent_volume = ("1.0 mL", {"number": 1.0, "unit": "mL"})
        elif entry == 16 and key == "methanol":
            solvent_volume = ("50 uL", {"number": 50, "unit": "uL"})
        add_part(r, key, "solvent", ev, volume=solvent_volume)
    add_part(r, "oxygen", "reagent", ev)
    add_part(r, "mesitylene", "additive", ev, stage="characterization")
    if entry == 16:
        r["participants"][2]["charges"]["amount"] = F("0.05 mmol", {"number": 0.05, "unit": "mmol"}, ev)
    set_scale(r, "1a", f"1a ({scale} mmol)", scale, ev)
    concentration = conc if conc is not None else 0.2
    total = ("1.5 mL", {"number": 1.5, "unit": "mL"}) if scale == 0.3 else None
    set_solvents(r, solvents, ev, ratio=ratio, total=total,
                 concentration=(f"{concentration} M", {"number": concentration, "unit": "M"}), basis="1a")
    add_step(r, ev, "reaction under O2 balloon", temp=("50 °C", {"number": 50, "unit": "°C"}),
             time=(f"{hours} h", {"number": hours, "unit": "h"}), atmosphere="O2 (balloon)")
    add_outcome(r, "nmr_yield_E2a", "yield", ["E2a"], ev, f"{yld}%", {"number": yld, "unit": "%"},
                yield_type="1H NMR yield", basis="1H NMR analysis of crude products with mesitylene as internal standard")
    if isolated is not None:
        add_outcome(r, "isolated_yield_E2a", "yield", ["E2a"], ev, f"({isolated})%", {"number": isolated, "unit": "%"},
                    yield_type="isolated yield", basis="isolated")
    add_outcome(r, "recovery_1a", "recovery", ["1a"], ev, f"{rec}%", {"number": rec, "unit": "%"},
                basis="1H NMR analysis of crude products with mesitylene as internal standard")
    if entry == 16:
        r["notes"].append(F("Table 2 entry 16 and Table 3 entry 1 both report 1a, 12.5 h, 62% isolated yield under the selected 0.95 M conditions.",
                            "Repeated source depiction treated as one experiment, not two.", ev))
    if yld > 0:
        set_class(r, ["C-C", "C-O", "DB-C-O"], ev, "Partial graphical classification of observed (E)-2a formation.")
    else:
        set_not_applicable(r, ev, "No (E)-2a formation observed; intended bond formation is not classified.")
    reactions.append(r)

# Table 3 entries 2-22; entry 1 is the already represented Table 2 entry 16.
for index, (letter, substituent, hours, yld) in enumerate(scope[1:], 2):
    ev = [EV[2]]
    scale = 10.0 if letter == "f" else 5.82 if letter == "n" else 1.0
    r = reaction(f"table3_entry{index:02d}", ev, categories=("substrate_scope",) + (("scale_up",) if letter in {"f", "n"} else ()),
                 location=f"Table 3 entry {index}")
    volumes = {"acetonitrile": (1.0, "mL"), "methanol": (50, "uL")} if scale == 1.0 else None
    standard_dimer_participants(r, f"1{letter}", f"E2{letter}", ev, substrate_amount=scale, catalyst=5, cu_equiv=2.5,
                                solvent_volumes=volumes,
                                catalyst_amount=0.05 if scale == 1.0 else None,
                                cu_amount=2.5 if scale == 1.0 else None)
    set_scale(r, f"1{letter}", f"1{letter} ({scale} mmol)", scale, ev)
    set_solvents(r, ["acetonitrile", "methanol"], ev, ratio="20/1")
    if scale == 1.0:
        r["solvent_conditions"]["concentration"] = F(
            None, {"number": 0.95, "unit": "M"}, ev, method="inference",
            derivation="Calculated as 1.0 mmol divided by the printed 1.0 mL CH3CN + 50 uL MeOH nominal volumes; not an author-reported measured total volume.")
        r["solvent_conditions"]["concentration_basis"] = F(f"1{letter}", mid(f"1{letter}"), ev)
    add_step(r, ev, "standard scope reaction", temp=("50 °C", {"number": 50, "unit": "°C"}),
             time=(f"{hours} h", {"number": hours, "unit": "h"}), atmosphere="O2 (balloon)")
    add_outcome(r, f"isolated_yield_E2{letter}", "yield", [f"E2{letter}"], ev, f"{yld}%", {"number": yld, "unit": "%"},
                yield_type="isolated yield", basis="isolated")
    if letter == "f":
        add_part(r, "E8f", "product", ev, product_role="byproduct")
        add_outcome(r, "byproduct_E8f", "other", ["E8f"], ev, "2 mol% byproduct (E)-8f", {"number": 2, "unit": "mol%"},
                    scope="reported byproduct amount; exact measurement basis unreported", basis_reported=False)
    if letter == "h":
        add_part(r, "piperonal", "product", ev, mass=("3.1 mg", {"number": 3.1, "unit": "mg"}), product_role="byproduct")
        add_part(r, "other_byproducts", "product", ev, mass=("8.1 mg", {"number": 8.1, "unit": "mg"}), product_role="byproduct mixture")
        add_outcome(r, "yield_piperonal", "yield", ["piperonal"], ev, "piperonyl aldehyde (3.1 mg, 2%)", {"number": 2, "unit": "%"},
                    yield_type="reported byproduct yield", basis="unreported")
        add_outcome(r, "amount_other_byproducts", "product_amount", ["other_byproducts"], ev, "mixture of other byproducts (8.1 mg)", {"number": 8.1, "unit": "mg"},
                    basis_reported=False)
        r["review"]["unresolved"].append("The 8.1 mg byproduct mixture has no reported composition.")
    set_class(r, ["C-C", "C-O", "DB-C-O"], ev, "Partial graphical classification of the Table 3 oxidative dimerization product.")
    reactions.append(r)

# Table 4: six optically active starting materials and product-specific yield/ee outcomes.
for entry, row in enumerate(chiral_rows, 1):
    imid, ilabel, config, iee, pmid, plabel, yld, pee, hours = row
    ev = [EV[3]]
    r = reaction(f"table4_entry{entry:02d}", ev, categories=("substrate_scope",), location=f"Table 4 entry {entry}")
    volumes = {"acetonitrile": (0.5, "mL"), "methanol": (25, "uL")}
    standard_dimer_participants(r, imid, pmid, ev, substrate_amount=0.5, catalyst=5, cu_equiv=2.5,
                                solvent_volumes=volumes, catalyst_amount=0.025, cu_amount=1.25)
    set_scale(r, imid, f"1 (0.5 mmol), specifically {ilabel}", 0.5, ev)
    set_solvents(r, ["acetonitrile", "methanol"], ev, ratio="20/1")
    add_step(r, ev, "optically active substrate scope reaction", temp=("50 °C", {"number": 50, "unit": "°C"}),
             time=(f"{hours} h", {"number": hours, "unit": "h"}), atmosphere="O2 (balloon)")
    add_outcome(r, f"starting_ee_{imid}", "selectivity", [imid], ev, f"{iee}% ee" if isinstance(iee, int) else f"{iee}% ee",
                {"number": int(str(iee).replace(">", "")), "unit": "%", **({"qualifier": "gt", "raw_text": ">99% ee"} if isinstance(iee, str) else {})},
                scope={"metric": "ee", "molecule": ilabel, "configuration": config}, basis_reported=False)
    add_outcome(r, f"isolated_yield_{pmid}", "yield", [pmid], ev, f"{yld}%", {"number": yld, "unit": "%"},
                yield_type="isolated yield", basis="isolated")
    add_outcome(r, f"product_ee_{pmid}", "selectivity", [pmid], ev, f"{pee}% ee" if isinstance(pee, int) else f"{pee}% ee",
                {"number": int(str(pee).replace(">", "")), "unit": "%", **({"qualifier": "gt", "raw_text": ">99% ee"} if isinstance(pee, str) else {})},
                scope={"metric": "ee", "molecule": plabel, "configuration": config}, basis_reported=False)
    r["review"]["unresolved"].append("The main text reports ee values but does not state the analytical method; measurement basis remains unreported.")
    set_class(r, ["C-C", "C-O", "DB-C-O"], ev, "Partial graphical classification; the printed stereochemistry is retained separately from bond formation.")
    reactions.append(r)

# Equation (1): cross-coupling of 1x and 1e, with two products.
r = reaction("eq1_cross_coupling", [EV[3]], categories=("substrate_scope",), location="Equation (1)",
             name="cross-coupling reaction of tertiary 2,3-allenol 1x and secondary 2,3-allenol 1e",
             transformation="cross-coupling/cyclization plus competing homo-coupling")
add_part(r, "1x", "reactant", [EV[3]], amount=("0.5 mmol", {"number": 0.5, "unit": "mmol"}), ref=None)
add_part(r, "1e", "reactant", [EV[3]], equiv=("1.0 equiv.", {"number": 1.0, "unit": "equiv"}), ref="1x")
add_part(r, "E2xe", "product", [EV[3]], product_role="cross-coupling product")
add_part(r, "E2e", "product", [EV[3]], product_role="homo-coupling product")
add_part(r, "rhcat", "catalyst", [EV[3]], molpct=("5 mol%", {"number": 5, "unit": "mol%"}), ref="1x")
add_part(r, "cuoac", "reagent", [EV[3]], equiv=("2.5 equiv.", {"number": 2.5, "unit": "equiv"}), ref="1x")
for key in ("acetonitrile", "methanol"):
    add_part(r, key, "solvent", [EV[3]])
add_part(r, "oxygen", "reagent", [EV[3]])
set_scale(r, "1x", "1x (0.5 mmol)", 0.5, [EV[3]])
set_solvents(r, ["acetonitrile", "methanol"], [EV[3]], ratio="20/1")
add_step(r, [EV[3]], "cross-coupling reaction", temp=("50 °C", {"number": 50, "unit": "°C"}),
         time=("11.5 h", {"number": 11.5, "unit": "h"}), atmosphere="O2 (balloon)")
add_outcome(r, "yield_E2xe", "yield", ["E2xe"], [EV[3]], "49%", {"number": 49, "unit": "%"}, basis_reported=False)
add_outcome(r, "yield_E2e", "yield", ["E2e"], [EV[3]], "17%", {"number": 17, "unit": "%"}, basis_reported=False)
r["review"]["unresolved"].append("The equation does not state whether the 49% and 17% values are isolated or analytically measured.")
set_class(r, ["C-C", "C-O", "DB-C-O"], [EV[3]], "Partial graphical classification of the two observed coupling/cyclization products.")
reactions.append(r)

# Scheme 2 gram-scale preparation of (S,E)-2f.
r = reaction("scheme2_scaleup_SE2f", [EV[3]], categories=("scale_up",), location="Scheme 2, standard-conditions gram-scale arrow")
standard_dimer_participants(r, "S1f", "SE2f", [EV[3]], substrate_amount=10.0, catalyst=5, cu_equiv=2.5)
add_part(r, "E8f", "product", [EV[3]], product_role="byproduct")
set_scale(r, "S1f", "(S)-1f, 10 mmol", 10.0, [EV[3]])
set_solvents(r, ["acetonitrile", "methanol"], [EV[3]], ratio="20/1")
add_step(r, [EV[3]], "standard conditions", temp=("50 °C (inherited from standard conditions)", {"number": 50, "unit": "°C"}),
         atmosphere="O2 (balloon; inherited from standard conditions)")
add_outcome(r, "yield_SE2f", "yield", ["SE2f"], [EV[3]], "50%", {"number": 50, "unit": "%"}, basis_reported=False)
add_outcome(r, "amount_SE2f", "product_amount", ["SE2f"], [EV[3]], "1.124 g", {"number": 1.124, "unit": "g"}, basis_reported=False)
add_outcome(r, "ee_SE2f", "selectivity", ["SE2f"], [EV[3]], ">99% ee", {"number": 99, "unit": "%", "qualifier": "gt", "raw_text": ">99% ee"},
            scope={"metric": "ee", "configuration": "S"}, basis_reported=False)
add_outcome(r, "byproduct_E8f", "other", ["E8f"], [EV[3]], "1%", {"number": 1, "unit": "%"}, scope="byproduct abundance; basis unreported", basis_reported=False)
r["review"]["unresolved"].extend(["The standard-conditions arrow does not print a reaction time.", "Yield and ee measurement methods are unreported in the main-text scheme."])
set_class(r, ["C-C", "C-O", "DB-C-O"], [EV[3]], "Partial graphical classification of the gram-scale oxidative dimerization.")
reactions.append(r)

# Scheme 2 derivatizations (i)-(iv), all explicitly under nitrogen.
deriv_specs = [
    ("i", "SE3f", 68, 65, 15, [("petasis", "reagent", "2.0 equiv.", 2.0)], ["toluene"], ["DB-C-C"], None),
    ("ii", "SE4f", 78, 40, 15, [("phenylacetylene", "reactant", "1.3 equiv.", 1.3), ("pd_sono", "catalyst", "3 mol%", 3),
                                            ("ipr2nh", "reagent", "3 equiv.", 3), ("cui", "catalyst", "6 mol%", 6)], ["thf"], ["C-C"], None),
    ("iii", "SE5f", 61, 85, 13, [("boronic", "reactant", "0.91 equiv.", 0.91), ("pd_suzuki", "catalyst", "3 mol%", 3),
                                              ("k3po4", "reagent", "3 equiv.", 3)], ["dioxane"], ["C-C"], "12/1 over (S,E)-6f"),
    ("iv", "SE6f", 65, 85, 15, [("boronic", "reactant", "3 equiv.", 3), ("pd_suzuki", "catalyst", "5 mol%", 5),
                                             ("k3po4", "reagent", "6 equiv.", 6)], ["dioxane"], ["C-C", "C-C"], None),
]
for roman, product, yld, temp, hours, agents, solvents, codes, ratio_note in deriv_specs:
    r = reaction(f"scheme2_deriv_{roman}", [EV[3]], categories=("derivatization",), location=f"Scheme 2 transformation ({roman})",
                 name=f"derivatization of (S,E)-2f, transformation ({roman})", transformation="product derivatization")
    add_part(r, "SE2f", "reactant", [EV[3]], equiv=("1.0 equiv. reference", {"number": 1.0, "unit": "equiv"}), ref="SE2f")
    add_part(r, product, "product", [EV[3]], product_role="target product")
    if roman == "iii":
        add_part(r, "SE6f", "product", [EV[3]], product_role="minor double-coupling product")
    for key, role, raw, number in agents:
        if "mol%" in raw:
            add_part(r, key, role, [EV[3]], molpct=(raw, {"number": number, "unit": "mol%"}), ref="SE2f")
        else:
            add_part(r, key, role, [EV[3]], equiv=(raw, {"number": number, "unit": "equiv"}), ref="SE2f")
    for key in solvents:
        add_part(r, key, "solvent", [EV[3]])
    set_solvents(r, solvents, [EV[3]])
    if roman == "i":
        r["solvent_conditions"]["concentration"] = F("0.1 M in toluene", {"number": 0.1, "unit": "M"}, [EV[3]])
        r["solvent_conditions"]["concentration_basis"] = F("Cp2TiMe2 solution", mid("petasis"), [EV[3]])
    add_step(r, [EV[3]], f"Scheme 2 transformation ({roman})", temp=(f"{temp} °C", {"number": temp, "unit": "°C"}),
             time=(f"{hours} h", {"number": hours, "unit": "h"}), atmosphere="nitrogen")
    add_outcome(r, f"yield_{product}", "yield", [product], [EV[3]], f"{yld}%", {"number": yld, "unit": "%"}, basis_reported=False)
    if roman in {"i", "ii", "iii"}:
        add_outcome(r, f"ee_{product}", "selectivity", [product], [EV[3]], ">99% ee", {"number": 99, "unit": "%", "qualifier": "gt", "raw_text": ">99% ee"},
                    scope={"metric": "ee", "configuration": "S"}, basis_reported=False)
    if ratio_note:
        add_outcome(r, "selectivity_5f_6f", "selectivity", ["SE5f", "SE6f"], [EV[3]], "12/1", "12/1",
                    scope={"metric": "product ratio", "order": ["(S,E)-5f", "(S,E)-6f"]}, basis_reported=False)
    r["review"]["unresolved"].append("The main-text Scheme 2 does not state the yield/ee analytical method.")
    set_class(r, codes, [EV[3]], f"Partial net bond formation inferred from Scheme 2 transformation ({roman}); no atom mapping asserted.")
    reactions.append(r)

# Scheme 3 mechanistic experiments.
r = reaction("scheme3_eq1_1w", [EV[4], TXT[3]], categories=("mechanistic_study", "control"), location="Scheme 3 equation (1)",
             name="cycloisomerization control with 1w", transformation="cycloisomerization")
standard_dimer_participants(r, "1w", "E2w", [EV[4]], substrate_amount=1.0, catalyst=5, cu_equiv=2.5)
add_part(r, "7w", "product", [EV[4]], product_role="observed cycloisomerization product")
set_scale(r, "1w", "1w, 1.0 mmol", 1.0, [EV[4]])
set_solvents(r, ["acetonitrile", "methanol"], [EV[4]], ratio="20/1")
add_step(r, [EV[4]], "standard-conditions mechanistic control", temp=("50 °C", {"number": 50, "unit": "°C"}),
         time=("12 h", {"number": 12, "unit": "h"}), atmosphere="O2 (balloon)")
add_outcome(r, "yield_E2w", "yield", ["E2w"], [EV[4]], "0%", {"number": 0, "unit": "%"}, basis_reported=False)
add_outcome(r, "yield_7w", "yield", ["7w"], [EV[4], TXT[3]], "14%", {"number": 14, "unit": "%"}, yield_type="isolated yield", basis="isolated")
set_class(r, ["C-O"], [EV[4]], "The observed 2,5-dihydrofuran 7w drawing supports a new ring C-O bond; full mapping remains unresolved.")
reactions.append(r)

r = reaction("scheme3_eq2_allenyl_ketone", [EV[4], TXT[3]], categories=("mechanistic_study", "control"), location="Scheme 3 equation (2)",
             name="control excluding allenyl ketone 8e as intermediate", transformation="no observed cross-coupling")
add_part(r, "1e", "reactant", [EV[4]])
add_part(r, "8e", "reactant", [EV[4]])
add_part(r, "E2e", "product", [EV[4]], product_role="intended product")
for key, role in (("rhcat", "catalyst"), ("cuoac", "reagent"), ("acetonitrile", "solvent"), ("methanol", "solvent"), ("oxygen", "reagent")):
    if key == "rhcat": add_part(r, key, role, [EV[4]], molpct=("5 mol%", {"number": 5, "unit": "mol%"}), ref="1e")
    elif key == "cuoac": add_part(r, key, role, [EV[4]], equiv=("2.5 equiv.", {"number": 2.5, "unit": "equiv"}), ref="1e")
    else: add_part(r, key, role, [EV[4]])
set_solvents(r, ["acetonitrile", "methanol"], [EV[4]], ratio="20/1")
add_step(r, [EV[4]], "mechanistic control", temp=("50 °C", {"number": 50, "unit": "°C"}),
         time=("12 h", {"number": 12, "unit": "h"}), atmosphere="O2 (balloon)")
add_outcome(r, "yield_E2e", "yield", ["E2e"], [EV[4]], "0% yield", {"number": 0, "unit": "%"}, basis_reported=False)
add_outcome(r, "recovery_1e", "recovery", ["1e"], [EV[4]], "no recovery of 1e", {"number": 0, "unit": "%"}, basis_reported=False)
add_outcome(r, "recovery_8e", "recovery", ["8e"], [EV[4]], "no recovery of 8e", {"number": 0, "unit": "%"}, basis_reported=False)
set_not_applicable(r, [EV[4]], "No desired product or other defined transformation was observed.")
reactions.append(r)

r = reaction("scheme3_eq3_alcohol9a", [EV[4], TXT[3]], categories=("mechanistic_study", "control"), location="Scheme 3 equation (3)",
             name="control excluding oxidation of alcohol 9a", transformation="no observed oxidation to (E)-2a")
add_part(r, "9a", "reactant", [EV[4]])
add_part(r, "E2a", "product", [EV[4]], product_role="intended product")
add_part(r, "rhcat", "catalyst", [EV[4]], molpct=("3 mol%", {"number": 3, "unit": "mol%"}), ref="9a")
add_part(r, "cuoac", "reagent", [EV[4]], equiv=("1.5 equiv.", {"number": 1.5, "unit": "equiv"}), ref="9a")
for key in ("acetonitrile", "methanol"):
    add_part(r, key, "solvent", [EV[4]])
add_part(r, "oxygen", "reagent", [EV[4]])
set_solvents(r, ["acetonitrile", "methanol"], [EV[4]], ratio="20/1")
add_step(r, [EV[4]], "mechanistic control", temp=("50 °C", {"number": 50, "unit": "°C"}),
         time=("20 h", {"number": 20, "unit": "h"}), atmosphere="O2 (balloon)")
add_outcome(r, "yield_E2a", "yield", ["E2a"], [EV[4]], "0% yield", {"number": 0, "unit": "%"}, basis_reported=False)
add_outcome(r, "recovery_9a", "recovery", ["9a"], [EV[4]], "87% recovery of 9a", {"number": 87, "unit": "%"}, basis_reported=False)
set_not_applicable(r, [EV[4]], "No (E)-2a formation was observed from alcohol 9a.")
reactions.append(r)

r = reaction("scheme3_eq4_Da1a", [EV[4], TXT[3]], categories=("mechanistic_study",), location="Scheme 3 equation (4)",
             name="deuterium-labeling experiment with [Dᵃ]-1a", transformation="oxidative dimerization without retained deuterium")
standard_dimer_participants(r, "Da1a", "E2a", [EV[4]], catalyst=5, cu_equiv=2.5, solvents=("acetonitrile", "cd3od"))
set_solvents(r, ["acetonitrile", "cd3od"], [EV[4]], ratio="20/1")
add_step(r, [EV[4]], "deuterium-labeling experiment", temp=("50 °C", {"number": 50, "unit": "°C"}),
         time=("12.5 h", {"number": 12.5, "unit": "h"}))
add_outcome(r, "starting_D_Da1a", "other", ["Da1a"], [EV[4]], "D = 78% D", {"number": 78, "unit": "%"}, scope="starting-material deuterium incorporation", basis_reported=False)
add_outcome(r, "yield_E2a", "yield", ["E2a"], [EV[4]], "44%", {"number": 44, "unit": "%"}, basis_reported=False)
add_outcome(r, "product_D_E2a", "other", ["E2a"], [EV[4], TXT[3]], "without any deuterium incorporation", {"number": 0, "unit": "%"},
            scope="product deuterium incorporation", basis_reported=False)
set_class(r, ["C-C", "C-O", "DB-C-O"], [EV[4]], "Partial graphical classification of observed (E)-2a; isotopic atom mapping remains unresolved.")
reactions.append(r)

r = reaction("scheme3_eq5_Db1m", [EV[4], TXT[4]], categories=("mechanistic_study",), location="Scheme 3 equation (5)",
             name="deuterium-labeling experiment with [Dᵇ]-1m", transformation="oxidative dimerization with 1,4-D migration")
standard_dimer_participants(r, "Db1m", "DE2m", [EV[4]], catalyst=5, cu_equiv=2.5)
set_solvents(r, ["acetonitrile", "methanol"], [EV[4]], ratio="20/1")
add_step(r, [EV[4]], "deuterium-labeling experiment", temp=("50 °C", {"number": 50, "unit": "°C"}),
         time=("13 h", {"number": 13, "unit": "h"}), atmosphere="O2 (balloon)")
add_outcome(r, "starting_D_Db1m", "other", ["Db1m"], [EV[4]], "D > 99% D", {"number": 99, "unit": "%", "qualifier": "gt", "raw_text": ">99% D"},
            scope="starting-material deuterium incorporation", basis_reported=False)
add_outcome(r, "yield_DE2m", "yield", ["DE2m"], [EV[4]], "61%", {"number": 61, "unit": "%"}, basis_reported=False)
add_outcome(r, "product_D_DE2m", "other", ["DE2m"], [EV[4], TXT[4]], "D = 99% D at newly formed deuterated group", {"number": 99, "unit": "%"},
            scope="site-specific product deuterium incorporation shown in Scheme 3", basis_reported=False)
set_class(r, ["C-C", "C-O", "DB-C-O"], [EV[4]], "Partial graphical classification; deuterium migration is retained as an outcome rather than a complete atom map.")
reactions.append(r)

# Scheme 4 is an author-proposed catalytic mechanism, not a directly observed experiment.
r = reaction("scheme4_proposed_mechanism", [EV[4]], record_type="proposed_mechanism", categories=("mechanistic_study",), location="Scheme 4",
             name="a possible mechanism", transformation="proposed oxyrhodation/insertion/1,4-Rh migration/beta-H elimination sequence")
add_part(r, "Db1m", "reactant", [EV[4]])
add_part(r, "DE2m", "product", [EV[4]], product_role="depicted final product")
add_part(r, "rhcat", "catalyst", [EV[4]])
add_part(r, "cuoac", "reagent", [EV[4]])
add_part(r, "oxygen", "reagent", [EV[4]])
r["procedure"]["original_text"] = F("Scheme 4 A possible mechanism.", "Author-proposed catalytic cycle", [EV[4]])
r["notes"].append(F("IN-1, IN-2/IN-2′, and IN-3 are proposed rhodium intermediates; the arrows are mechanistic proposals, not separate experiments.",
                    "Proposed intermediates are not treated as experimentally isolated products.", [EV[4]]))
r["classification"]["coverage"] = F(None, "unresolved", [EV[4]], method="inference", derivation="A proposed elementary mechanism is not classified as an observed net experiment without reviewed atom mapping.")
r["classification"]["labels"] = F(None, None, [EV[4]], method="inference", derivation="No observed net bond classification asserted for the proposed mechanism.")
r["classification"]["limitations"] = ["Proposed mechanism only; not direct experimental verification."]
r["review"]["unresolved"].append("Intermediate structures IN-1/IN-2/IN-2′/IN-3 are not encoded as isolated molecular identities.")
reactions.append(r)

bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": f"{DOC}:scheme1_prior_work", "document_id": DOC, "evidence_ids": [EV[1]],
     "kind": "literature_context", "note": "Scheme 1a-b summarizes cited Pd(II) and Au(III) prior work; it is not imported as current-paper experiments."},
    {"context_id": f"{DOC}:fig1_xray", "document_id": DOC, "evidence_ids": [EV[3]],
     "kind": "characterization", "note": "Fig. 1 is an ORTEP/X-ray representation of (E)-2k (CCDC 1877067), not a separate synthesis."},
    {"context_id": f"{DOC}:mechanistic_discussion", "document_id": DOC, "evidence_ids": [TXT[4], EV[4]],
     "kind": "mechanistic_discussion", "note": "The prose interprets Scheme 3 and proposes Scheme 4; author interpretation is kept separate from observed control records."},
    {"context_id": f"{DOC}:references_pages5_6", "document_id": DOC, "evidence_ids": [EV[5], EV[6]],
     "kind": "bibliography", "note": "Pages 5-6 contain references only; reactions named in cited titles are not current-paper experiments."},
]

(OUT / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
inventory = {
    "document_id": DOC,
    "machine_records_remain_needs_review": True,
    "molecule_count": len(bundle["molecules"]),
    "reaction_count": len(bundle["reactions"]),
    "record_types": {kind: sum(r["record_type"] == kind for r in bundle["reactions"]) for kind in sorted({r["record_type"] for r in bundle["reactions"]})},
    "specific_experiment_count": sum(r["record_type"] == "specific_experiment" for r in bundle["reactions"]),
    "notes": [
        "Table 2 entry 16 and Table 3 entry 1 are represented once with duplicate-source note.",
        "No SI file was associated; main-text-only missing details remain unreported.",
        "Complex graphical structures were not converted into guessed definitive SMILES.",
    ],
}
(OUT / "raw_model/inventory.json").write_text(json.dumps(inventory, ensure_ascii=False, indent=2) + "\n")
stereo_diagnostics = {}
for key, m in molecules.items():
    smiles = m["smiles"]["normalized_value"]
    if smiles:
        stereo_diagnostics[m["molecule_id"]] = st.stereo_inspect(smiles)
(OUT / "raw_model/stereo_diagnostics.json").write_text(
    json.dumps(stereo_diagnostics, ensure_ascii=False, indent=2) + "\n")
print(json.dumps(inventory, ensure_ascii=False, indent=2))
