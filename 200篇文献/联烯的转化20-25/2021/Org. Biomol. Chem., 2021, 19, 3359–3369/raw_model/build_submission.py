#!/usr/bin/env python3
"""Build the source-derived v3.1 submission for doc_ce3af83f9dcdab7f.

The data below are transcribed from visual inspection of PDF pages 1--11.
Pages 7--11 were reviewed in the current bounded session; earlier page evidence
is retained only to merge cross-page occurrences of the same experiments.
No record is promoted beyond needs_review by this script.
"""
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_ce3af83f9dcdab7f"
RUN = ROOT / "runs/current"
E = {p: f"{DOC}:page{p}:text" for p in range(1, 12)}
EI = {p: f"{DOC}:page{p}:image" for p in range(1, 12)}

bundle = json.loads((RUN / "raw_model/submitted_template_pages_0007_0011.json").read_text())
bundle["molecules"] = []
bundle["reactions"] = []
bundle["source_context"] = []

molecules = {}


def fact(raw=None, normalized=None, pages=(), method="original_extraction", status=None, derivation=None):
    ev = []
    for p in pages:
        ev.append(E[p])
        if p in EI:
            ev.append(EI[p])
    return st.fact(raw, normalized, ev, method=method, status=status, derivation=derivation)


def page_ev(*pages):
    out = []
    for p in pages:
        out.append(E[p])
        if p in EI:
            out.append(EI[p])
    return list(dict.fromkeys(out))


def missing(status="unreported"):
    return st.fact(status=status, method="not_applicable" if status == "not_applicable" else "original_extraction")


def inferred(normalized, pages, derivation):
    return fact(None, normalized, pages, method="inference", derivation=derivation)


def derived(normalized, pages, derivation):
    return fact(None, normalized, pages, method="structure_derivation", derivation=derivation)


def conflict(alternatives, pages):
    """alternatives: [(raw, normalized, pages, method)]"""
    f = fact([x[0] for x in alternatives], None, pages, status="conflict")
    f["alternatives"] = [fact(raw, norm, ps, method=method) for raw, norm, ps, method in alternatives]
    return f


def add_molecule(mid, label, pages, name=None, smiles=None, stereo="achiral", note=None):
    if mid in molecules:
        return molecules[mid]
    ev = page_ev(*pages)
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"] = [{
        "document_id": DOC,
        "compound_label": fact(label, label, pages),
        "evidence_ids": ev,
    }]
    display = name or label
    m["names"]["original"] = fact(display, display, pages)
    m["names"]["english"] = fact(display, display, pages)
    m["cas_numbers"] = missing()
    m["generic_template"] = None
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = derived(
            smiles,
            pages,
            "Connectivity candidate derived from the source drawing/name; RDKit parsing is a graph diagnostic and does not verify the source structure.",
        )
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing()
    if stereo == "unspecified":
        m["stereochemistry"]["chirality_status"] = fact(
            {"source_drawing": "stereogenic sulfoximine/other center shown without an R/S or wedge/dash assignment"},
            "unspecified", pages,
        )
        m["stereochemistry"]["geometry_status"] = fact(
            {"source_drawing": "no E/Z or cis/trans assignment shown"}, "not_applicable", pages,
        )
        m["stereochemistry"]["representation_status"] = fact(
            {"candidate_smiles": "connectivity only; source configuration not assigned"}, "not_encoded", pages,
        )
        m["stereochemistry"]["assignment"] = fact(
            {"assignment": "not stated in pages inspected"},
            {"absolute_configuration": None, "relative_configuration": None}, pages,
        )
    elif stereo == "unresolved":
        m["stereochemistry"]["chirality_status"] = fact(
            {"source_drawing": "configuration cannot be safely resolved at available scale"}, "unresolved", pages,
        )
        m["stereochemistry"]["geometry_status"] = fact(
            {"source_drawing": "geometry cannot be safely resolved at available scale"}, "unresolved", pages,
        )
        m["stereochemistry"]["representation_status"] = fact(
            {"structure": "label/substituent-defined only"}, "unresolved", pages,
        )
    else:
        m["stereochemistry"]["chirality_status"] = derived(
            "achiral", pages, "No tetrahedral/axial stereogenic element is present in the derived connectivity candidate."
        )
        m["stereochemistry"]["geometry_status"] = derived(
            "not_applicable", pages, "No source-assigned E/Z or ring cis/trans relationship applies to this connectivity candidate."
        )
        m["stereochemistry"]["representation_status"] = derived(
            "not_applicable", pages, "No stereochemical encoding is required for the derived connectivity candidate."
        )
    m["stereochemistry"]["atom_centers"] = missing()
    m["stereochemistry"]["double_bonds"] = missing()
    m["stereochemistry"]["axial"] = missing()
    m["stereochemistry"]["mixture"] = missing()
    m["stereochemistry"]["isomeric_smiles"] = missing()
    m["stereochemistry"]["identity_key"] = missing()
    if smiles and stereo == "achiral":
        diag = st.stereo_inspect(smiles)
        m["stereochemistry"]["atom_centers"] = derived(
            diag["atom_centers"], pages, "standard stereo graph diagnostic; empty center list supports the achiral candidate only at graph level."
        )
        m["stereochemistry"]["isomeric_smiles"] = derived(
            diag["canonical_isomeric_smiles"], pages, "Canonical isomeric SMILES generated by standard stereo from the source-derived connectivity candidate."
        )
        m["stereochemistry"]["identity_key"] = derived(
            diag["stereo_graph_key"], pages, "stereo_graph_key generated by standard stereo; not a sample identity."
        )
    m["identity_note"] = fact(note, note, pages) if note else missing()
    m["external_properties"] = []
    molecules[mid] = m
    return m


def ensure_agent(mid, label, pages=(3,), name=None):
    return add_molecule(mid, label, pages, name=name, stereo="unresolved")


def participant(pid, mid, role, pages, stage="reaction", mass=None, amount=None, volume=None,
                equivalents=None, mol_percent=None, reference=None, product_role=None):
    p = st.new_participant(pid, mid, role, page_ev(*pages), stage=stage)
    if product_role:
        p["product_role"] = fact(product_role, product_role, pages)
    if mass is not None:
        p["charges"]["mass"] = fact(str(mass), mass, pages)
    if amount is not None:
        p["charges"]["amount"] = fact(str(amount), amount, pages)
    if volume is not None:
        p["charges"]["volume"] = fact(str(volume), volume, pages)
    if equivalents is not None:
        p["charges"]["equivalents"] = fact(str(equivalents), equivalents, pages)
    if mol_percent is not None:
        p["charges"]["mol_percent"] = fact(str(mol_percent), mol_percent, pages)
    if reference is not None:
        p["equivalence_reference"] = fact(reference, reference, pages)
    return p


def set_charge_conflict(p, key, alternatives, pages):
    p["charges"][key] = conflict(alternatives, pages)


def step(step_id, pages, participant_ids, operation=None, temperature=None, time=None,
         atmosphere=None, pressure=None, stage="reaction"):
    s = st.new_step(step_id, page_ev(*pages), stage=stage)
    s["participant_ids"] = list(participant_ids)
    if operation is not None:
        s["operation"] = fact(operation, operation, pages)
    if temperature is not None:
        s["temperature"] = fact(str(temperature), temperature, pages)
    if time is not None:
        s["time"] = fact(str(time), time, pages)
    if atmosphere is not None:
        s["atmosphere"] = fact(atmosphere, atmosphere, pages)
    if pressure is not None:
        s["pressure"] = fact(str(pressure), pressure, pages)
    return s


def set_step_conflict(s, key, alternatives, pages):
    s[key] = conflict(alternatives, pages)


def outcome(oid, kind, targets, pages, value=None, yield_type=None, basis=None, scope=None):
    o = st.new_outcome(oid, kind, targets, page_ev(*pages))
    if value is not None:
        o["value"] = fact(str(value), value, pages)
    if yield_type is not None:
        o["yield_type"] = fact(yield_type, yield_type, pages)
    if basis is not None:
        raw, normalized, method, derivation = basis
        o["measurement_basis"] = fact(raw, normalized, pages, method=method, derivation=derivation)
    if scope is not None:
        o["scope"] = fact(scope, scope, pages)
    return o


def classification(labels, pages, derivation, coverage="partial"):
    c = {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": missing(),
        "bond_changes": [],
        "labels": missing(),
        "coverage": inferred(coverage, pages, derivation),
        "limitations": [
            "Bond identities are inferred from the source-visible net transformation without atom mapping.",
            "Hydrogen partners and omitted coproducts are not fully mapped.",
        ],
        "model_votes": [],
    }
    if coverage == "unresolved":
        return c
    changes = []
    grouped = {}
    for i, code in enumerate(labels, 1):
        if code == "C-H":
            atoms = [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}]
            kind = "hydrogen_gain"
        else:
            a, b = code.split("-")[-2:]
            atoms = [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}]
            kind = "formed"
        ch = {"change_id": f"b{i}", "atoms": atoms, "before_order": 0, "after_order": 1, "change_type": kind}
        if kind == "hydrogen_gain":
            ch["hydrogen_count"] = 1
        changes.append(ch)
        lab = st.bond_label(atoms[0]["element"], atoms[1]["element"], 1)
        grouped.setdefault(lab["code"], dict(lab, supporting_change_ids=[]))["supporting_change_ids"].append(ch["change_id"])
    c["bond_changes"] = [derived(ch, pages, derivation) for ch in changes]
    c["labels"] = derived(list(grouped.values()), pages, derivation)
    return c


def no_reaction_classification(pages):
    return {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": missing(),
        "bond_changes": [],
        "labels": inferred([], pages, "No observed product transformation is reported; no expected bond label is assigned."),
        "coverage": inferred("not_applicable", pages, "No observed product transformation is reported."),
        "limitations": ["Negative control; expected chemistry is not classified as observed chemistry."],
        "model_votes": [],
    }


def base_reaction(rid, pages, record_type, categories, location, reaction_name, transformation,
                  unresolved=None):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, page_ev(*pages))
    r["record_type"] = record_type
    r["experiment_categories"] = fact(categories, categories, pages) if categories is not None else missing("not_applicable")
    r["location"] = fact(location, location, pages)
    r["reaction_name"] = fact(reaction_name, reaction_name, pages) if reaction_name else missing()
    r["transformation_type"] = inferred(
        transformation, pages, "Transformation description inferred from the source reaction drawing and prose."
    ) if transformation else missing()
    r["review"] = {
        "status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
        "unresolved": list(unresolved or ["Drawing-derived identities and bond changes require expert source review."]),
    }
    r["completeness"] = "partial"
    return r


def set_solvents(r, pages, composition, volumes=None, ratio=None):
    r["solvent_conditions"]["composition"] = fact(composition, composition, pages)
    if volumes is not None:
        r["solvent_conditions"]["total_volume"] = fact(volumes, volumes, pages)
    if ratio is not None:
        r["solvent_conditions"]["mixture_ratio"] = fact(ratio, ratio, pages)


# Core molecules and source-visible candidates.
add_molecule("m_1a", "1a", (3, 5), "N-(2-iodophenyl)-4-methylbenzenesulfonamide",
             "Cc1ccc(S(=O)(=O)Nc2ccccc2I)cc1", "achiral")
add_molecule("m_2a", "2a", (3, 5), "3-bromoprop-1-yne", "C#CCBr", "achiral")
add_molecule("m_3a", "3a", (2, 3, 4, 5), "methyl phenyl NH-sulfoximine",
             "CS(=O)(=N)c1ccccc1", "unspecified",
             "The sulfur stereogenic center is not assigned in the inspected drawings; no racemic composition is inferred.")
add_molecule("m_4a", "4a", (2, 3, 4, 5),
             "N-(methyl(oxo)(phenyl)-λ6-sulfanylidene)-2-(1-tosyl-1H-indol-3-yl)acetamide",
             "CS(=O)(=NC(=O)Cc1cn(S(=O)(=O)c2ccc(C)cc2)c2ccccc12)c1ccccc1", "unspecified")
add_molecule("m_8", "8", (4, 5, 9, 10),
             "methyl(phenyl)((2-(1-tosyl-1H-indol-3-yl)ethyl)imino)-λ6-sulfanone",
             "CS(=O)(=NCCc1cn(S(=O)(=O)c2ccc(C)cc2)c2ccccc12)c1ccccc1", "unspecified")
add_molecule("m_acetanilide_direct", "N-(2-iodophenyl)acetamide", (4,),
             "N-(2-iodophenyl)acetamide", "CC(=O)Nc1ccccc1I", "achiral")
add_molecule("m_5", "5", (4,), "N-(2-iodophenyl)-N-(prop-2-yn-1-yl)acetamide",
             "CC(=O)N(CC#C)c1ccccc1I", "achiral")
add_molecule("m_6a", "6a", (4, 5, 9),
             "2-(1-acetyl-1H-indol-3-yl)-N-(methyl(oxo)(phenyl)-λ6-sulfanylidene)acetamide",
             "CS(=O)(=NC(=O)Cc1cn(C(C)=O)c2ccccc12)c1ccccc1", "unspecified")
add_molecule("m_7a", "7a", (4, 5, 9), "Ts-protected N-allenyl-2-iodoaniline 7a",
             "C=C=CN(c1ccccc1I)S(=O)(=O)c1ccc(C)cc1", "unresolved",
             "Connectivity candidate only; the source does not assign an allene axis and the terminal allene drawing requires expert review.")
add_molecule("m_CO", "CO", (2, 3, 4, 5), "carbon monoxide", "[C-]#[O+]", "achiral")

# Common agents/solvents are source-identified by label; no unsupported structure is asserted.
for mid, label in [
    ("m_K2CO3", "K2CO3"), ("m_acetone", "acetone"), ("m_tBuOK", "tBuOK"),
    ("m_THF", "THF"), ("m_PdOAc2", "Pd(OAc)2"), ("m_PPh3", "PPh3"),
    ("m_DBU", "DBU"), ("m_toluene", "toluene"), ("m_pyridine", "pyridine (Py)"),
    ("m_TEA", "triethylamine (TEA)"), ("m_Cs2CO3", "Cs2CO3"), ("m_K3PO4", "K3PO4"),
    ("m_DABCO", "DABCO"), ("m_PdMeCNCl2", "Pd(CH3CN)2Cl2"),
    ("m_PdTFA2", "Pd(OCOCF3)2"), ("m_Pddba2", "Pd(dba)2"), ("m_PdCl2", "PdCl2"),
    ("m_PTol3", "P(4-MeC6H4)3"), ("m_PAnis3", "P(4-MeOC6H4)3"),
    ("m_PFPh3", "P(4-FC6H4)3"), ("m_PCy3", "PCy3"), ("m_DPPM", "DPPM"),
    ("m_SPhos", "SPhos"), ("m_dioxane", "1,4-dioxane"), ("m_DMF", "DMF"),
    ("m_MeCN", "CH3CN"), ("m_MnCO5Br", "Mn(CO)5Br"),
    ("m_silica", "silica gel"), ("m_EtOAc", "ethyl acetate"),
    ("m_petroleum_ether", "petroleum ether"), ("m_BH3Me2S", "BH3·Me2S"),
    ("m_CH2Cl2", "CH2Cl2"), ("m_MeOH", "MeOH"), ("m_water", "H2O"),
    ("m_MgSO4", "MgSO4"),
]:
    ensure_agent(mid, label, (3, 4, 5))

# Generic/source-family molecule labels for Schemes 1 and 2.
for mid, label in [("m_scheme1_substrate", "N-aryl-N-allenamide (Scheme 1)"),
                   ("m_scheme1_product", "ketenimination product (Scheme 1)"),
                   ("m_family_1", "ortho-iodoaniline family 1"),
                   ("m_family_2", "propargyl bromide family 2"),
                   ("m_family_3", "NH-sulfoximine family 3"),
                   ("m_family_4", "β-indolyl sulfoximidoyl amide family 4")]:
    add_molecule(mid, label, (2,), stereo="unresolved")

# Scheme 1: explicitly prior work, not a current-paper experiment.
r = base_reaction("scheme1_prior_work", (2,), "literature_example", None,
                  {"page": 2, "region": "Scheme 1"},
                  "Ketenimination reaction of N-aryl-N-allenamides with isocyanides",
                  "prior-work ketenimination",
                  ["Generic cited literature transformation; no concrete current-paper experiment is asserted."])
r["participants"] = [participant("lit_s", "m_scheme1_substrate", "reactant", (2,)),
                     participant("lit_p", "m_scheme1_product", "product", (2,), product_role="depicted product family")]
r["classification"] = classification([], (2,), "Generic prior-work scheme is not atom-mapped.", coverage="unresolved")
bundle["reactions"].append(r)

# Scheme 2 eqn (1): proposed mechanism only.
r = base_reaction("scheme2_proposed_mechanism", (2,), "proposed_mechanism", None,
                  {"page": 2, "region": "Scheme 2, eqn (1), intermediates A--F"},
                  "Proposed synthetic mechanism", "proposed Pd annulation/carbonyl insertion/nucleophilic addition",
                  ["Arrow-pushing proposal; intermediates A--F are not isolated experiments.",
                   "Exact atom mapping and elementary-step bond accounting remain unresolved."])
r["classification"] = classification(["C-N", "C-C", "C-C", "C-N"], (2,),
                                     "The paper states that four bonds are formed overall; element-pair labels are inferred from the depicted proposal.")
bundle["reactions"].append(r)

# Scheme 2 eqn (2): current-paper generic scheme.
r = base_reaction("scheme2_generic", (2,), "generic_scheme", None,
                  {"page": 2, "region": "Scheme 2, eqn (2)"},
                  "Pd-catalyzed four-component one-pot reaction", "indole annulation/carbonyl insertion/C-N bond formation",
                  ["Generic R-group scheme; concrete instances are recorded separately.",
                   "R-group structures and atom mapping require expert review."])
r["participants"] = [participant("g1", "m_family_1", "reactant", (2,)),
                     participant("g2", "m_family_2", "reactant", (2,)),
                     participant("gco", "m_CO", "reagent", (2,)),
                     participant("g3", "m_family_3", "reactant", (2,)),
                     participant("g4", "m_family_4", "product", (2,), product_role="desired product family")]
r["classification"] = classification(["C-N", "C-C", "C-C", "C-N"], (2,),
                                     "Net labels inferred from the four-component source drawing; no atom mapping is asserted.")
bundle["reactions"].append(r)


table_rows = [
    (1, "m_PdOAc2", "m_PPh3", "m_pyridine", "m_toluene", {"number": 0, "unit": "%"}),
    (2, "m_PdOAc2", "m_PPh3", "m_TEA", "m_toluene", "trace"),
    (3, "m_PdOAc2", "m_PPh3", "m_K2CO3", "m_toluene", {"number": 27, "unit": "%"}),
    (4, "m_PdMeCNCl2", "m_PPh3", "m_K2CO3", "m_toluene", {"number": 20, "unit": "%"}),
    (5, "m_PdTFA2", "m_PPh3", "m_K2CO3", "m_toluene", {"number": 19, "unit": "%"}),
    (6, "m_Pddba2", "m_PPh3", "m_K2CO3", "m_toluene", {"number": 11, "unit": "%"}),
    (7, "m_PdCl2", "m_PPh3", "m_K2CO3", "m_toluene", {"number": 12, "unit": "%"}),
    (8, "m_PdOAc2", "m_PPh3", "m_Cs2CO3", "m_toluene", {"number": 34, "unit": "%"}),
    (9, "m_PdOAc2", "m_PPh3", "m_K3PO4", "m_toluene", {"number": 29, "unit": "%"}),
    (10, "m_PdOAc2", "m_PPh3", "m_tBuOK", "m_toluene", {"number": 40, "unit": "%"}),
    (11, "m_PdOAc2", "m_PPh3", "m_tBuOK", "m_THF", {"number": 26, "unit": "%"}),
    (12, "m_PdOAc2", "m_PPh3", "m_DBU", "m_toluene", {"number": 94, "unit": "%"}),
    (13, "m_PdOAc2", "m_PPh3", "m_DABCO", "m_toluene", {"number": 10, "unit": "%"}),
    (14, "m_PdOAc2", "m_PTol3", "m_DBU", "m_toluene", {"number": 90, "unit": "%"}),
    (15, "m_PdOAc2", "m_PAnis3", "m_DBU", "m_toluene", {"number": 86, "unit": "%"}),
    (16, "m_PdOAc2", "m_PFPh3", "m_DBU", "m_toluene", {"number": 82, "unit": "%"}),
    (17, "m_PdOAc2", "m_PCy3", "m_DBU", "m_toluene", "trace"),
    (18, "m_PdOAc2", "m_DPPM", "m_DBU", "m_toluene", {"number": 84, "unit": "%"}),
    (19, "m_PdOAc2", "m_SPhos", "m_DBU", "m_toluene", {"number": 75, "unit": "%"}),
    (20, "m_PdOAc2", "m_PPh3", "m_DBU", "m_THF", {"number": 80, "unit": "%"}),
    (21, "m_PdOAc2", "m_PPh3", "m_DBU", "m_dioxane", {"number": 87, "unit": "%"}),
    (22, "m_PdOAc2", "m_PPh3", "m_DBU", "m_DMF", {"number": 82, "unit": "%"}),
    (23, "m_PdOAc2", "m_PPh3", "m_DBU", "m_MeCN", {"number": 92, "unit": "%"}),
    (24, "m_PdOAc2", "m_PPh3", "m_DBU", "mixed", {"number": 88, "unit": "%"}),
    (25, "m_PdOAc2", "m_PPh3", "m_DBU", "m_toluene", {"number": 66, "unit": "%"}),
]


def add_common_table_participants(r, row, pd, ligand, base, solvent):
    pages = (3, 4, 5) if row == 12 else (3,)
    parts = [
        participant("p_1a", "m_1a", "reactant", pages, amount={"number": 0.44, "unit": "mmol"}),
        participant("p_2a", "m_2a", "reactant", pages, amount={"number": 0.44, "unit": "mmol"}),
        participant("p_k2co3", "m_K2CO3", "reagent", pages, amount={"number": 0.44, "unit": "mmol"}),
        participant("p_acetone", "m_acetone", "solvent", pages, volume={"number": 1.0, "unit": "mL"}),
        participant("p_tbuok", "m_tBuOK", "reagent", pages, amount={"number": 0.56, "unit": "mmol"}),
        participant("p_thf", "m_THF", "solvent", pages, volume={"number": 1.5, "unit": "mL"}),
        participant("p_3a", "m_3a", "reactant", pages, amount={"number": 0.40, "unit": "mmol"}),
        participant("p_pd", pd, "catalyst", pages, mol_percent={"number": 5, "unit": "mol%"}, reference="m_3a"),
        participant("p_ligand", ligand, "ligand", pages, mol_percent={"number": 20, "unit": "mol%"}, reference="m_3a"),
        participant("p_base", base, "reagent", pages, amount={"number": 0.80, "unit": "mmol"}),
        participant("p_product", "m_4a", "product", pages, product_role="desired product"),
    ]
    if row == 12:
        set_charge_conflict(parts[1], "amount", [
            ("0.44 mmol (Table 1 footnote / Scheme 3 caption)", {"number": 0.44, "unit": "mmol"}, (3, 4), "original_extraction"),
            ("45 μL, 1.2 mmol (experimental procedure)", {"number": 1.2, "unit": "mmol"}, (5,), "original_extraction"),
        ], (3, 4, 5))
        parts[1]["charges"]["volume"] = fact("45 μL", {"number": 45, "unit": "uL"}, (5,))
        set_charge_conflict(parts[5], "volume", [
            ("1.5 mL (Table 1 footnote / Scheme 3 caption)", {"number": 1.5, "unit": "mL"}, (3, 4), "original_extraction"),
            ("1.0 mL (experimental procedure)", {"number": 1.0, "unit": "mL"}, (5,), "original_extraction"),
        ], (3, 4, 5))
    if row == 24:
        parts[5]["charges"]["volume"] = fact("THF (1.0 mL) in step II was not removed", {"number": 1.0, "unit": "mL"}, (3,))
        parts.append(participant("p_toluene", "m_toluene", "solvent", (3,), volume={"number": 1.0, "unit": "mL"}))
    else:
        smid = solvent
        parts.append(participant("p_solvent", smid, "solvent", pages, volume={"number": 2.0, "unit": "mL"}))
    if row == 25:
        parts.append(participant("p_mn", "m_MnCO5Br", "reagent", (3,), equivalents={"number": 1.0, "unit": "equiv"}, reference="m_3a"))
    else:
        parts.append(participant("p_co", "m_CO", "reagent", pages))
    if row == 12:
        parts += [
            participant("p_silica", "m_silica", "workup_agent", (5,), stage="workup"),
            participant("p_etoac", "m_EtOAc", "workup_agent", (5,), stage="workup", volume={"number": 15.0, "unit": "mL"}),
            participant("p_pe", "m_petroleum_ether", "purification_agent", (5,), stage="purification"),
        ]
    r["participants"] = parts


for row, pd, ligand, base, solvent, yld in table_rows:
    pages = (3, 4, 5) if row == 12 else (2, 3) if row in (1, 2, 3) else (3,)
    cats = ["optimization", "substrate_scope"] if row == 12 else ["optimization"]
    loc = {"page": 3, "region": f"Table 1 entry {row}"}
    if row == 12:
        loc["additional_regions"] = ["Scheme 3 product 4a, page 4", "General procedure/product 4a, page 5"]
    r = base_reaction(f"table1_entry_{row:02d}", pages, "specific_experiment", cats, loc,
                      "Pd-catalyzed multi-step one-pot cascade", "four-component indole annulation/carbonylation",
                      ["Drawing-derived structures and bond changes require expert review.",
                       "Table reaction drawing and footnote report different durations for steps I--III."])
    if row == 12:
        r["review"]["unresolved"] += [
            "Propargyl bromide amount conflicts: 0.44 mmol in Table 1/Scheme 3 versus 45 μL, 1.2 mmol in the experimental procedure.",
            "THF volume conflicts: 1.5 mL in Table 1/Scheme 3 versus 1.0 mL in the experimental procedure.",
            "Scheme 3 caption cites Table 1 entry 8 although its listed DBU conditions correspond to entry 12.",
        ]
    if row == 24:
        r["review"]["unresolved"].append("Entry 24 table ratio THF:toluene = 1:2 conflicts with footnote b volumes of 1.0 mL THF plus 1.0 mL toluene.")
    add_common_table_participants(r, row, pd, ligand, base, solvent)
    r["reaction_scale"] = fact("NH-sulfoximine 3a (0.40 mmol)", {"reference_molecule_id": "m_3a", "amount": {"number": 0.40, "unit": "mmol"}}, pages)
    if row == 24:
        set_solvents(r, (3,), ["acetone", "THF", "toluene"], {"step I": "1.0 mL acetone", "step II retained": "1.0 mL THF", "step III added": "1.0 mL toluene"})
        r["solvent_conditions"]["mixture_ratio"] = conflict([
            ("THF:toluene = 1:2 (table row)", "1:2", (3,), "original_extraction"),
            ("1.0 mL THF retained + 1.0 mL toluene added (footnote b)", "1:1 by stated volumes", (3,), "original_extraction"),
        ], (3,))
    else:
        set_solvents(r, pages, ["acetone", "THF", molecules[solvent]["names"]["original"]["normalized_value"]],
                     {"step I": "1.0 mL", "step II": "1.5 mL", "step III": "2.0 mL"})
    s1 = step("step_I_heat", pages, ["p_1a", "p_2a", "p_k2co3", "p_acetone"],
              "combine 1a, 2a and K2CO3 in acetone; heat", {"number": 60, "unit": "°C"})
    set_step_conflict(s1, "time", [
        ("1 h (reaction drawing)", {"number": 1, "unit": "h"}, (3,), "original_extraction"),
        ("2 h (Table 1 footnote)", {"number": 2, "unit": "h"}, (3,), "original_extraction"),
    ] + ([("1 h (experimental procedure)", {"number": 1, "unit": "h"}, (5,), "original_extraction")] if row == 12 else []), pages)
    s2a = step("step_II_cool_and_add", pages, ["p_tbuok", "p_thf"],
               "add tBuOK in THF at 0 °C", {"number": 0, "unit": "°C"})
    s2b = step("step_II_warm_hold", pages, ["p_tbuok", "p_thf"],
               "allow to warm from 0 °C to room temperature and stir",
               {"number": 0, "unit": "°C", "end": "r.t.", "raw_text": "0 °C to r.t."})
    set_step_conflict(s2b, "time", [
        ("1 h (reaction drawing)", {"number": 1, "unit": "h"}, (3,), "original_extraction"),
        ("2 h (Table 1 footnote)", {"number": 2, "unit": "h"}, (3,), "original_extraction"),
    ] + ([("additional 1 h (experimental procedure)", {"number": 1, "unit": "h"}, (5,), "original_extraction")] if row == 12 else []), pages)
    atm = "N2" if row == 25 else "CO balloon"
    press = None if row == 25 else "1 atm"
    s3 = step("step_III_carbonylation", pages,
              [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"],
              "add catalyst, ligand, base, NH-sulfoximine and step-III solvent; stir", "r.t.", atmosphere=atm, pressure=press)
    set_step_conflict(s3, "time", [
        ("2–6 h (reaction drawing)", {"minimum": 2, "maximum": 6, "unit": "h"}, (3,), "original_extraction"),
        ("2–3 h (Table 1 footnote)", {"minimum": 2, "maximum": 3, "unit": "h"}, (3,), "original_extraction"),
    ] + ([("2–3 h (experimental procedure)", {"minimum": 2, "maximum": 3, "unit": "h"}, (5,), "original_extraction")] if row == 12 else []), pages)
    r["steps"] = [s1, s2a, s2b, s3]
    if row == 12:
        r["steps"] += [
            step("workup_filter", (5,), ["p_silica", "p_etoac"], "filter through a short pad of silica gel and wash with ethyl acetate", stage="workup"),
            step("purification_column", (5,), ["p_pe", "p_etoac", "p_product"], "silica gel column chromatography; petroleum ether/ethyl acetate = 10:1–3:1 (v/v)", stage="purification"),
        ]
    basis = None
    if row in (3, 12, 25):
        basis = ("isolated" if row != 12 else "purified by silica gel column chromatography to afford 4a", "isolated", "original_extraction", None)
    kind = "trace" if yld == "trace" else "yield"
    r["outcomes"] = [outcome(f"yield_{row}", kind, ["m_4a"], pages, yld, "reported yield" if kind == "yield" else None, basis,
                              "Table 1 product 4a")]
    if row == 12:
        r["outcomes"].append(outcome("amount_4a", "product_amount", ["m_4a"], (5,), {"number": 175.3, "unit": "mg"},
                                            basis=("purified product mass", "isolated mass", "original_extraction", None), scope="product 4a"))
        r["procedure"]["original_text"] = fact(
            "General procedure for the synthesis of β-indolyl sulfoximidoyl amide 4a, Steps I–III and chromatography",
            "Three-stage one-pot procedure followed by filtration and column chromatography", (5,))
        r["procedure"]["addition_order"] = fact(
            "Step III: Pd(OAc)2, PPh3, DBU, NH-sulfoximine 3a, and dry toluene were added sequentially; purged with CO balloon 3 times",
            ["Pd(OAc)2", "PPh3", "DBU", "3a", "dry toluene", "CO purge ×3"], (5,))
    else:
        r["procedure"]["original_text"] = fact("Table 1 footnote a conditions", "Table 1 three-stage optimization protocol", (3,))
    r["classification"] = classification(["C-N", "C-C", "C-C", "C-N"], pages,
                                         "Four new bonds are stated in the paper; C–C/C–N identities are inferred from the source-visible net cascade.")
    r["general_scheme_ids"] = [f"{DOC}:scheme2_generic"]
    r["notes"] = [fact(
        "The reaction drawing shows 1 h / 1 h / 2–6 h, whereas Table 1 footnote a shows 2 h / 2 h / 2–3 h.",
        "Source duration conflict retained without choosing a preferred value.", (3,))]
    if row == 12:
        r["notes"] += [fact(
            "Scheme 3 caption cites Table 1, entry 8 while listing Pd(OAc)2/PPh3/DBU conditions and 4a 94%, which match entry 12.",
            "Unresolved source cross-reference conflict.", (3, 4)),
            fact("2a: Table/Scheme 0.44 mmol versus experimental procedure 45 μL, 1.2 mmol; THF 1.5 mL versus 1.0 mL.",
                 "Charge/volume conflicts retained in structured fields.", (3, 4, 5))]
    bundle["reactions"].append(r)


scope = [
    ("4b", 70, "4-chloro indole product", "1"), ("4c", 93, "5-methyl indole product", "1"),
    ("4d", 75, "5-methoxy indole product", "1"), ("4e", 79, "5-fluoro indole product", "1"),
    ("4f", 83, "5-chloro indole product", "1"), ("4g", 88, "6-methyl indole product", "1"),
    ("4h", 83, "6-fluoro indole product", "1"), ("4i", 88, "6-chloro indole product", "1"),
    ("4j", 70, "6-bromo indole product", "1"), ("4k", 82, "5-trifluoromethyl indole product", "1"),
    ("4l", 46, "R2=H, R3=Ph product", "2"), ("4m", 38, "R2=Me, R3=Ph product", "2"),
    ("4n", 90, "S-(2-methylphenyl), S-methyl product", "3"),
    ("4o", 72, "S-(2-chlorophenyl), S-methyl product", "3"),
    ("4p", 86, "S-(2-bromophenyl), S-methyl product", "3"),
    ("4q", 93, "S-(4-methoxyphenyl), S-methyl product", "3"),
    ("4r", 89, "S-(4-chlorophenyl), S-methyl product", "3"),
    ("4s", 68, "S-ethyl, S-phenyl product", "3"), ("4t", 63, "S-cyclopropyl, S-phenyl product", "3"),
    ("4u", 75, "S-benzyl, S-phenyl product", "3"), ("4v", 92, "S-methyl, S-naphthyl product", "3"),
    ("4w", 87, "S,S-diphenyl product", "3"), ("4x", 88, "S-phenyl, S-p-tolyl product", "3"),
    ("4y", 79, "S,S-bis(4-chlorophenyl) product", "3"),
    ("4z", 74, "S-(4-benzoylphenyl), S-p-tolyl product", "3"),
]
page6_masses = {"4b": 140.2, "4c": 178.6, "4d": 148.8, "4e": 153.0, "4f": 166.0,
                "4g": 169.0, "4h": 160.7, "4i": 176.0, "4j": 152.3}
page7_masses = {"4k": 175.2, "4l": 99.7, "4m": 84.5, "4n": 172.8,
                "4o": 144.0, "4p": 187.1, "4q": 184.5, "4r": 178.0}
page8_masses = {"4s": 130.6, "4t": 124.0, "4u": 162.6, "4v": 199.0,
                "4w": 183.8, "4x": 190.8, "4y": 188.3, "4z": 191.2}
scope_mass_pages = {**{k: (6, v) for k, v in page6_masses.items()},
                    **{k: (7, v) for k, v in page7_masses.items()},
                    **{k: (8, v) for k, v in page8_masses.items()}}
scope_source_names = {
    "4k": "N-(Methyl(oxo)(phenyl)-λ6-sulfanylidene)-2-(1-tosyl-5-(trifluoromethyl)-1H-indol-3-yl)acetamide",
    "4l": "N-(Methyl(oxo)(phenyl)-λ6-sulfanylidene)-2-phenyl-2-(1-tosyl-1H-indol-3-yl)acetamide",
    "4m": "N-(Methyl(oxo)(phenyl)-λ6-sulfanylidene)-2-(2-methyl-1-tosyl-1H-indol-3-yl)-2-phenylacetamide",
    "4n": "N-(Methyl(oxo)(o-tolyl)-λ6-sulfanylidene)-2-(1-tosyl-1H-indol-3-yl)acetamide",
    "4o": "N-((2-Chlorophenyl)(methyl)(oxo)-λ6-sulfanylidene)-2-(1-tosyl-1H-indol-3-yl)acetamide",
    "4p": "N-((2-Bromophenyl)(methyl)(oxo)-λ6-sulfanylidene)-2-(1-tosyl-1H-indol-3-yl)acetamide",
    "4q": "N-((4-Methoxyphenyl)(methyl)(oxo)-λ6-sulfanylidene)-2-(1-tosyl-1H-indol-3-yl)acetamide",
    "4r": "N-((4-Chlorophenyl)(methyl)(oxo)-λ6-sulfanylidene)-2-(1-tosyl-1H-indol-3-yl)acetamide",
    "4s": "N-(Ethyl(oxo)(phenyl)-λ6-sulfanylidene)-2-(1-tosyl-1H-indol-3-yl)acetamide",
    "4t": "N-(Cyclopropyl(oxo)(phenyl)-λ6-sulfanylidene)-2-(1-tosyl-1H-indol-3-yl)acetamide",
    "4u": "N-(Benzyl(oxo)(phenyl)-λ6-sulfanylidene)-2-(1-tosyl-1H-indol-3-yl)acetamide",
    "4v": "N-(Methyl(naphthalen-2-yl)(oxo)-λ6-sulfanylidene)-2-(1-tosyl-1H-indol-3-yl)acetamide",
    "4w": "N-(Oxodiphenyl-λ6-sulfanylidene)-2-(1-tosyl-1H-indol-3-yl)acetamide",
    "4x": "N-(Oxo(phenyl)(p-tolyl)-λ6-sulfanylidene)-2-(1-tosyl-1H-indol-3-yl)acetamide",
    "4y": "N-(Bis(4-chlorophenyl)(oxo)-λ6-sulfanylidene)-2-(1-tosyl-1H-indol-3-yl)acetamide",
    "4z": "N-((4-Benzoylphenyl)(oxo)(p-tolyl)-λ6-sulfanylidene)-2-(1-tosyl-1H-indol-3-yl)acetamide",
}

for label, y, desc, variable in scope:
    product_mid = f"m_{label}"
    product_pages = (4, scope_mass_pages[label][0]) if label in scope_mass_pages else (4,)
    # The 4j characterization begins on page 6 and its final lines continue at
    # the top of page 7; retain both pages without creating a second reaction.
    if label == "4j":
        product_pages = (4, 6, 7)
    add_molecule(product_mid, label, product_pages, scope_source_names.get(label, desc), stereo="unspecified",
                 note="Product drawing contains a sulfoximine sulfur center without source-assigned absolute configuration; no racemate is inferred.")
    if variable == "1":
        var_mid = f"m_1_for_{label}"
        add_molecule(var_mid, f"1 ({desc})", (4,), f"2-iodoaniline 1 variant giving {label}", stereo="unresolved")
        mids = [var_mid, "m_2a", "m_3a"]
    elif variable == "2":
        var_mid = f"m_2_for_{label}"
        add_molecule(var_mid, f"2 ({desc})", (4,), f"propargyl bromide 2 variant giving {label}", stereo="unresolved")
        mids = ["m_1a", var_mid, "m_3a"]
    else:
        var_mid = f"m_3_for_{label}"
        add_molecule(var_mid, f"3 ({desc})", (4,), f"NH-sulfoximine 3 variant giving {label}", stereo="unspecified")
        mids = ["m_1a", "m_2a", var_mid]
    pages = tuple(dict.fromkeys(product_pages + ((5,) if label != "4z" else ())))
    r = base_reaction(f"scheme3_{label}", pages, "specific_experiment", ["substrate_scope"],
                      {"page": 4, "region": f"Scheme 3 product {label}"},
                      "Pd-catalyzed multi-step one-pot cascade", "four-component indole annulation/carbonylation",
                      ["Drawing-derived substituent identity and bond changes require expert review.",
                       "Scheme drawing reports 1 h for steps I/II while the caption reports 2 h.",
                       "Scheme 3 caption cites Table 1 entry 8 although listed DBU conditions correspond to entry 12."])
    r["participants"] = [
        participant("p_sub1", mids[0], "reactant", pages, amount={"number": 0.44, "unit": "mmol"}),
        participant("p_sub2", mids[1], "reactant", pages, amount={"number": 0.44, "unit": "mmol"}),
        participant("p_k2co3", "m_K2CO3", "reagent", pages, amount={"number": 0.44, "unit": "mmol"}),
        participant("p_acetone", "m_acetone", "solvent", pages, volume={"number": 1.0, "unit": "mL"}),
        participant("p_tbuok", "m_tBuOK", "reagent", pages, amount={"number": 0.56, "unit": "mmol"}),
        participant("p_thf", "m_THF", "solvent", pages, volume={"number": 1.5, "unit": "mL"}),
        participant("p_sub3", mids[2], "reactant", pages, amount={"number": 0.40, "unit": "mmol"}),
        participant("p_co", "m_CO", "reagent", pages),
        participant("p_pd", "m_PdOAc2", "catalyst", pages, mol_percent={"number": 5, "unit": "mol%"}, reference=mids[2]),
        participant("p_ligand", "m_PPh3", "ligand", pages, mol_percent={"number": 20, "unit": "mol%"}, reference=mids[2]),
        participant("p_dbu", "m_DBU", "reagent", pages, amount={"number": 0.80, "unit": "mmol"}),
        participant("p_toluene", "m_toluene", "solvent", pages, volume={"number": 2.0, "unit": "mL"}),
        participant("p_product", product_mid, "product", pages, product_role="desired product"),
    ]
    if label != "4z":
        r["participants"] += [
            participant("p_silica", "m_silica", "workup_agent", (5,), stage="workup"),
            participant("p_etoac", "m_EtOAc", "workup_agent", (5,), stage="workup", volume={"number": 15.0, "unit": "mL"}),
            participant("p_pe", "m_petroleum_ether", "purification_agent", (5,), stage="purification"),
        ]
    r["reaction_scale"] = fact(f"NH-sulfoximine 3 ({mids[2]}) 0.40 mmol", {"reference_molecule_id": mids[2], "amount": {"number": 0.40, "unit": "mmol"}}, pages)
    set_solvents(r, pages, ["acetone", "THF", "toluene"], {"step I": "1.0 mL", "step II": "1.5 mL", "step III": "2.0 mL"})
    s1 = step("step_I_heat", pages, ["p_sub1", "p_sub2", "p_k2co3", "p_acetone"], "heat in acetone", {"number": 60, "unit": "°C"})
    set_step_conflict(s1, "time", [
        ("1 h (Scheme 3 drawing)", {"number": 1, "unit": "h"}, (4,), "original_extraction"),
        ("2 h (Scheme 3 caption)", {"number": 2, "unit": "h"}, (4,), "original_extraction"),
    ], pages)
    s2a = step("step_II_cool_add", pages, ["p_tbuok", "p_thf"], "add tBuOK in THF at 0 °C", {"number": 0, "unit": "°C"})
    s2b = step("step_II_warm_hold", pages, ["p_tbuok", "p_thf"], "warm from 0 °C to room temperature",
               {"number": 0, "unit": "°C", "end": "r.t.", "raw_text": "0 °C to r.t."})
    set_step_conflict(s2b, "time", [
        ("1 h (Scheme 3 drawing)", {"number": 1, "unit": "h"}, (4,), "original_extraction"),
        ("2 h (Scheme 3 caption)", {"number": 2, "unit": "h"}, (4,), "original_extraction"),
    ], pages)
    s3 = step("step_III", pages, [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"],
              "carbonylation/annulation stage", "r.t.", {"minimum": 2, "maximum": 6, "unit": "h"}, "CO", "1 atm")
    r["steps"] = [s1, s2a, s2b, s3]
    basis = None
    if label != "4z":
        basis = ("Products 4b to 4y were prepared through a similar procedure to 4a; 4a was purified by silica gel chromatography",
                 "isolated after silica gel chromatography (inherited general procedure)", "inherited_general_condition", None)
        r["steps"] += [
            step("workup_filter", (5,), ["p_silica", "p_etoac"], "filter through silica pad and wash with ethyl acetate (inherited general procedure)", stage="workup"),
            step("purification_column", (5,), ["p_pe", "p_etoac", "p_product"], "silica gel column chromatography, petroleum ether/ethyl acetate 10:1–3:1 (inherited general procedure)", stage="purification"),
        ]
    r["outcomes"] = [outcome(f"yield_{label}", "yield", [product_mid], pages, {"number": y, "unit": "%"}, "reported yield", basis, f"Scheme 3 product {label}")]
    if label in scope_mass_pages:
        mass_page, mass_value = scope_mass_pages[label]
        r["outcomes"].append(outcome(f"amount_{label}", "product_amount", [product_mid], (mass_page,), {"number": mass_value, "unit": "mg"},
                                            basis=("Yield: reported mass and percentage", "reported product mass", "original_extraction", None), scope=f"characterization heading for {label}"))
    r["procedure"]["original_text"] = fact("Scheme 3 reaction drawing and caption", "Three-stage one-pot scope procedure", (4,))
    r["procedure"]["general_procedure_reference"] = fact("The other products 4b to 4y were prepared through a similar procedure to that of product 4a", "General procedure for 4a", (5,)) if label != "4z" else missing()
    r["classification"] = classification(["C-N", "C-C", "C-C", "C-N"], pages,
                                         "Net element-pair labels inferred from the source-visible cascade; no atom mapping is asserted.")
    r["general_scheme_ids"] = [f"{DOC}:scheme2_generic"]
    r["notes"] = [fact("Scheme 3 caption says Table 1, entry 8, but lists DBU/PPh3 conditions corresponding to entry 12.",
                       "Unresolved source cross-reference conflict retained.", (4,))]
    bundle["reactions"].append(r)


# Concrete failed direct Ac-protected four-component trial described in page-4 prose.
r = base_reaction("ac_protected_direct_no_reaction", (4,), "specific_experiment", ["control"],
                  {"page": 4, "region": "prose immediately above Scheme 4"},
                  "Ac-protected four-component trial", "attempted one-pot cascade",
                  ["Exact charges, catalyst loading, time, and solvent are not stated for this failed trial.",
                   "No desired product is created in the record."])
r["participants"] = [participant("p_ac", "m_acetanilide_direct", "reactant", (4,)),
                     participant("p_2a", "m_2a", "reactant", (4,)), participant("p_co", "m_CO", "reagent", (4,)),
                     participant("p_3a", "m_3a", "reactant", (4,))]
r["steps"] = [step("attempt", (4,), ["p_ac", "p_2a", "p_co", "p_3a"], "reaction did not work under the Pd-catalyzed conditions", atmosphere="CO", pressure="1 atm")]
r["outcomes"] = [outcome("no_reaction", "no_reaction", [], (4,), True, scope="no desired product reported")]
r["classification"] = no_reaction_classification((4,))
bundle["reactions"].append(r)

# Scheme 4, distinct propargyl amide route to 6a.
r = base_reaction("scheme4_5_to_6a", (4,), "specific_experiment", ["substrate_scope"],
                  {"page": 4, "region": "Scheme 4"},
                  "Pd-catalyzed three-component tandem annulation", "propargyl amide annulation/carbonylation",
                  ["Yield measurement basis is not stated in the inspected source region.",
                   "Drawing-derived bond changes require expert review."])
r["participants"] = [participant("p_5", "m_5", "reactant", (4,)), participant("p_3a", "m_3a", "reactant", (4,)),
                     participant("p_co", "m_CO", "reagent", (4,)),
                     participant("p_pd", "m_PdOAc2", "catalyst", (4,), mol_percent={"number": 5, "unit": "mol%"}, reference="m_3a"),
                     participant("p_lig", "m_PPh3", "ligand", (4,), mol_percent={"number": 10, "unit": "mol%"}, reference="m_3a"),
                     participant("p_dbu", "m_DBU", "reagent", (4,), equivalents={"number": 2.0, "unit": "equiv"}, reference="m_3a"),
                     participant("p_tol", "m_toluene", "solvent", (4,), volume={"number": 2.0, "unit": "mL"}),
                     participant("p_6a", "m_6a", "product", (4,), product_role="desired product")]
set_solvents(r, (4,), ["toluene"], {"reaction": "2.0 mL"})
r["steps"] = [step("reaction", (4,), [p["participant_id"] for p in r["participants"]], "stir under CO", "r.t.", {"number": 12, "unit": "h"}, "CO", "1 atm")]
r["outcomes"] = [outcome("yield_6a_42", "yield", ["m_6a"], (4,), {"number": 42, "unit": "%"}, "reported yield", None, "Scheme 4 product 6a")]
r["classification"] = classification(["C-C", "C-C", "C-N"], (4,), "Net bond types inferred from the source-visible transformation of preformed propargyl amide 5.")
bundle["reactions"].append(r)


# Scheme 5 variants: separate records for each product/trace substrate.
add_molecule("m_7_ac", "7 (R=Ac)", (5, 8, 9), "Ac-protected N-allenyl-2-iodoaniline 7",
             "CC(=O)N(C=C=C)c1ccccc1I", stereo="unresolved")
add_molecule("m_7_boc", "7 (R=Boc)", (5, 8, 9), "Boc-protected N-allenyl-2-iodoaniline 7",
             "CC(C)(C)OC(=O)N(C=C=C)c1ccccc1I", stereo="unresolved")
scheme5_product_names = {
    "6b": "tert-Butyl-3-(2-((methyl(oxo)(phenyl)-λ6-sulfanylidene)amino)-2-oxoethyl)-1H-indole-1-carboxylate",
    "6c": "tert-Butyl-3-(2-((methyl(naphthalen-2-yl)(oxo)-λ6-sulfanylidene)amino)-2-oxoethyl)-1H-indole-1-carboxylate",
    "6d": "tert-Butyl-3-(2-((methyl(oxo)(o-tolyl)-λ6-sulfanylidene)amino)-2-oxoethyl)-1H-indole-1-carboxylate",
    "6e": "tert-Butyl-3-(2-(((4-chlorophenyl)(methyl)(oxo)-λ6-sulfanylidene)amino)-2-oxoethyl)-1H-indole-1-carboxylate",
}
for suf, desc in scheme5_product_names.items():
    add_molecule(f"m_{suf}", suf, (5, 9), desc, stereo="unspecified")
for suf, desc in [("6c", "S-naphthyl NH-sulfoximine 3"), ("6d", "S-o-tolyl NH-sulfoximine 3"),
                  ("6e", "S-(4-chlorophenyl) NH-sulfoximine 3")]:
    add_molecule(f"m_3_{suf}", f"3 ({desc})", (5,), desc, stereo="unspecified")
scheme5_rows = [
    ("6a", "m_7_ac", "m_3a", 84, "yield"), ("6b", "m_7_boc", "m_3a", 90, "yield"),
    ("6c", "m_7_boc", "m_3_6c", 87, "yield"), ("6d", "m_7_boc", "m_3_6d", 93, "yield"),
    ("6e", "m_7_boc", "m_3_6e", 81, "yield"),
]
scheme5_masses = {"6a": 119.0, "6b": 148.4, "6c": 160.8, "6d": 158.5, "6e": 144.5}
for rgroup in ("H", "Me", "Bn"):
    sm = f"m_7_{rgroup.lower()}"
    pm = f"m_6f_{rgroup.lower()}"
    add_molecule(sm, f"7 (R={rgroup})", (5,), f"N-allenyl-2-iodoaniline 7, R={rgroup}", stereo="unresolved")
    add_molecule(pm, f"6f (R={rgroup})", (5,), f"desired product 6f, R={rgroup}", stereo="unspecified")
    scheme5_rows.append((f"6f_{rgroup}", sm, "m_3a", "trace", "trace"))

for label, start_mid, sul_mid, y, kind in scheme5_rows:
    product_mid = "m_6a" if label == "6a" else f"m_{label.lower()}"
    if product_mid not in molecules:
        product_mid = f"m_6f_{label.split('_')[-1].lower()}"
    record_pages = (4, 5, 8, 9) if kind == "yield" else (4, 5)
    unresolved = [
        "Scheme 5 arrow states PPh3 10 mol% whereas the caption and page-9 general procedure state 20 mol%.",
        "Scheme 5 reports 2–6 h whereas the pages 8–9 general procedure reports 2–3 h.",
        "The pages 8–9 general-procedure heading/body calls the starting material N-allenyl-2-iodoaniline 5, while Scheme 5 and the surrounding prose identify it as N-allenylaniline 7.",
        "Product stereochemistry at sulfur is not assigned.",
    ]
    if kind == "trace":
        unresolved.append("The trace result has no stated measurement basis or isolated product characterization.")
    r = base_reaction(f"scheme5_{label.lower()}", record_pages, "specific_experiment", ["substrate_scope"],
                      {"page": 5, "region": f"Scheme 5 product {label}"},
                      "Pd-catalyzed annulation of N-allenylaniline 7", "allenamide annulation/carbonylation",
                      unresolved)
    procedure_pages = (5, 8, 9) if kind == "yield" else (5,)
    parts = [participant("p_7", start_mid, "reactant", procedure_pages, amount={"number": 0.48, "unit": "mmol"}),
             participant("p_3", sul_mid, "reactant", procedure_pages, amount={"number": 0.40, "unit": "mmol"}),
             participant("p_co", "m_CO", "reagent", procedure_pages),
             participant("p_pd", "m_PdOAc2", "catalyst", procedure_pages, mol_percent={"number": 5, "unit": "mol%"}, reference=sul_mid),
             participant("p_lig", "m_PPh3", "ligand", procedure_pages, mol_percent={"number": 10, "unit": "mol%"}, reference=sul_mid),
             participant("p_dbu", "m_DBU", "reagent", procedure_pages, amount={"number": 0.80, "unit": "mmol"}),
             participant("p_tol", "m_toluene", "solvent", procedure_pages, volume={"number": 2.0, "unit": "mL"}),
             participant("p_product", product_mid, "product", (5, 9) if kind == "yield" else (5,), product_role="desired product")]
    if kind == "yield":
        parts += [
            participant("p_silica", "m_silica", "workup_agent", (8, 9), stage="workup"),
            participant("p_etoac", "m_EtOAc", "workup_agent", (8, 9), stage="workup", volume={"number": 15.0, "unit": "mL"}),
            participant("p_pe", "m_petroleum_ether", "purification_agent", (8, 9), stage="purification"),
        ]
    set_charge_conflict(parts[4], "mol_percent", [
        ("10 mol% (Scheme 5 arrow)", {"number": 10, "unit": "mol%"}, (5,), "original_extraction"),
        ("20 mol% (Scheme 5 caption and pages 8–9 general procedure)", {"number": 20, "unit": "mol%"}, (5, 8, 9), "original_extraction"),
    ], record_pages)
    r["participants"] = parts
    r["reaction_scale"] = fact("NH-sulfoximine 3 (0.40 mmol)", {"reference_molecule_id": sul_mid, "amount": {"number": 0.40, "unit": "mmol"}}, procedure_pages)
    set_solvents(r, procedure_pages, ["toluene"], {"reaction": "2.0 mL"})
    reaction_step = step("reaction", procedure_pages, [p["participant_id"] for p in parts if p["stage"] == "reaction"], "purge with CO three times and stir", "r.t.", {"minimum": 2, "maximum": 6, "unit": "h"}, "CO", "1 atm")
    if kind == "yield":
        set_step_conflict(reaction_step, "time", [
            ("2–6 h (Scheme 5)", {"minimum": 2, "maximum": 6, "unit": "h"}, (5,), "original_extraction"),
            ("2–3 h (pages 8–9 general procedure)", {"minimum": 2, "maximum": 3, "unit": "h"}, (8, 9), "original_extraction"),
        ], procedure_pages)
    r["steps"] = [reaction_step]
    if kind == "yield":
        r["steps"] += [
            step("workup_filter", (8, 9), ["p_silica", "p_etoac"], "filter through a short silica pad and wash with ethyl acetate", stage="workup"),
            step("purification_column", (8, 9), ["p_pe", "p_etoac", "p_product"], "silica gel column chromatography, petroleum ether/ethyl acetate 10:1–3:1 (v/v)", stage="purification"),
        ]
    value = "trace" if kind == "trace" else {"number": y, "unit": "%"}
    basis = ("filtered and purified by silica gel column chromatography to afford compound 6",
             "isolated after silica gel chromatography", "inherited_general_condition", None) if kind == "yield" else None
    r["outcomes"] = [outcome(f"outcome_{label}", kind, [product_mid], (4, 5, 8, 9) if kind == "yield" else (4, 5), value,
                              "reported yield" if kind == "yield" else None, basis, f"Scheme 5 product {label}")]
    if kind == "yield":
        r["outcomes"].append(outcome(f"amount_{label}", "product_amount", [product_mid], (9,),
                                             {"number": scheme5_masses[label], "unit": "mg"},
                                             basis=("Yield: reported mass and percentage", "reported isolated product mass", "original_extraction", None),
                                             scope=f"characterization heading for {label}"))
    r["procedure"]["original_text"] = fact("Scheme 5 reaction drawing/caption and pages 8–9 general procedure", "CO-purged room-temperature annulation followed by silica workup/purification" if kind == "yield" else "CO-purged room-temperature annulation", procedure_pages)
    r["classification"] = classification(["C-C", "C-C", "C-N"], (5,), "Net element-pair labels inferred from the source-visible transformation of preformed N-allenylaniline 7.")
    r["notes"] = [
        fact("PPh3 10 mol% in Scheme 5 arrow; 20 mol% in Scheme 5 caption/pages 8–9 procedure.", "Ligand-loading conflict retained.", record_pages),
        fact("Scheme 5: N-allenylaniline 7; pages 8–9 general procedure: N-allenyl-2-iodoaniline 5.", "Source compound-number inconsistency retained; the record follows the Scheme 5 transformation and does not silently rename the source text.", record_pages),
        fact("Scheme 5: 2–6 h; pages 8–9 general procedure: 2–3 h.", "Reaction-time conflict retained.", record_pages),
    ]
    bundle["reactions"].append(r)


# Scheme 6 eqn (1): gram-scale experiment; page 9 is visually verified in this session.
r = base_reaction("scheme6_gram_scale_4a", (4, 5, 9), "specific_experiment", ["scale_up"],
                  {"page": 5, "region": "Scheme 6 eqn (1)", "additional_region": "page 9 gram-scale procedure"},
                  "Gram-scale synthesis of 4a", "allenamide annulation/carbonylation",
                  ["Scheme 6 shows PPh3 10 mol% and toluene 2.0 mL; page-9 extracted text states PPh3 20 mol% and toluene 10 mL.",
                   "Scheme 6 and the page-9 procedure are both retained as evidence; neither conflicting value is preferred."])
parts = [participant("p_7a", "m_7a", "reactant", (5, 9), mass={"number": 987.0, "unit": "mg"}, amount={"number": 2.4, "unit": "mmol"}),
         participant("p_3a", "m_3a", "reactant", (5, 9), mass={"number": 310, "unit": "mg"}, amount={"number": 2.0, "unit": "mmol"}),
         participant("p_co", "m_CO", "reagent", (5, 9)),
         participant("p_pd", "m_PdOAc2", "catalyst", (5, 9), mol_percent={"number": 5, "unit": "mol%"}, reference="m_3a"),
         participant("p_lig", "m_PPh3", "ligand", (5, 9), mol_percent={"number": 10, "unit": "mol%"}, reference="m_3a"),
         participant("p_dbu", "m_DBU", "reagent", (5, 9), equivalents={"number": 2.0, "unit": "equiv"}, reference="m_3a"),
         participant("p_tol", "m_toluene", "solvent", (5, 9), volume={"number": 2.0, "unit": "mL"}),
         participant("p_4a", "m_4a", "product", (5, 9), product_role="desired product")]
set_charge_conflict(parts[4], "mol_percent", [
    ("10 mol% (Scheme 6)", {"number": 10, "unit": "mol%"}, (5,), "original_extraction"),
    ("0.4 mmol, 20 mol% (page-9 extracted procedure)", {"number": 20, "unit": "mol%"}, (9,), "original_extraction"),
], (5, 9))
set_charge_conflict(parts[6], "volume", [
    ("2.0 mL (Scheme 6)", {"number": 2.0, "unit": "mL"}, (5,), "original_extraction"),
    ("10 mL (page-9 extracted procedure)", {"number": 10, "unit": "mL"}, (9,), "original_extraction"),
], (5, 9))
r["participants"] = parts
r["reaction_scale"] = fact("3a, 310 mg, 2.0 mmol", {"reference_molecule_id": "m_3a", "amount": {"number": 2.0, "unit": "mmol"}}, (5, 9))
r["solvent_conditions"]["composition"] = fact("toluene", ["toluene"], (5, 9))
r["solvent_conditions"]["total_volume"] = conflict([
    ("2.0 mL (Scheme 6)", {"number": 2.0, "unit": "mL"}, (5,), "original_extraction"),
    ("10 mL (page-9 extracted procedure)", {"number": 10, "unit": "mL"}, (9,), "original_extraction"),
], (5, 9))
r["steps"] = [step("reaction", (5, 9), [p["participant_id"] for p in parts], "purge with CO three times and stir", "r.t.", {"number": 6, "unit": "h"}, "CO balloon", "1 atm")]
r["outcomes"] = [outcome("yield_gram_4a", "yield", ["m_4a"], (4, 5, 9), {"number": 86, "unit": "%"}, "reported yield",
                                 ("purified by silica gel column chromatography to afford 4a, 801 mg, 86% yield", "isolated after silica gel chromatography", "original_extraction", None),
                                 "gram-scale 7a to 4a"),
                 outcome("amount_gram_4a", "product_amount", ["m_4a"], (5, 9), {"number": 801, "unit": "mg"},
                         basis=("afford 4a, 801 mg", "isolated product mass", "original_extraction", None), scope="gram-scale product 4a")]
r["classification"] = classification(["C-C", "C-C", "C-N"], (5,), "Net bond types inferred from preformed allenamide 7a to product 4a.")
r["notes"] = [fact("Scheme 6: PPh3 10 mol%, toluene 2.0 mL; page-9 extracted procedure: PPh3 20 mol%, toluene 10 mL.", "Conflicts retained without choosing.", (5, 9))]
bundle["reactions"].append(r)


# Scheme 6 eqn (2): reduction of 4a to 8; pages 9--10 are visually verified.
r = base_reaction("scheme6_reduction_4a_to_8", (4, 5, 9, 10), "specific_experiment", ["derivatization"],
                  {"page": 5, "region": "Scheme 6 eqn (2)", "additional_regions": ["page 9 procedure and product heading", "page 10 characterization continuation"]},
                  "Reduction of 4a to tryptamine-containing sulfoximine 8", "amide carbonyl reduction",
                  ["The page-9 procedure calls the purified material 'product 7' although Scheme 6 and the characterization heading identify compound 8.",
                   "Hydrogen atom sources are not completely mapped."])
r["participants"] = [participant("p_4a", "m_4a", "reactant", (5, 9), amount={"number": 0.3, "unit": "mmol"}),
                     participant("p_bh3", "m_BH3Me2S", "reagent", (5, 9), equivalents={"number": 3.0, "unit": "equiv"}, reference="m_4a"),
                     participant("p_dcm", "m_CH2Cl2", "solvent", (5, 9), volume={"number": 3.0, "unit": "mL"}),
                     participant("p_meoh", "m_MeOH", "workup_agent", (9,), stage="workup"),
                     participant("p_water", "m_water", "workup_agent", (9,), stage="workup"),
                     participant("p_mgso4", "m_MgSO4", "workup_agent", (9,), stage="workup"),
                     participant("p_silica", "m_silica", "purification_agent", (9,), stage="purification"),
                     participant("p_8", "m_8", "product", (5, 10), product_role="desired product")]
r["reaction_scale"] = fact("4a (0.3 mmol)", {"reference_molecule_id": "m_4a", "amount": {"number": 0.3, "unit": "mmol"}}, (9,))
set_solvents(r, (5, 9), ["CH2Cl2"], {"reaction": "3.0 mL"})
s_warm = step("warm_hold", (5, 9), ["p_4a", "p_bh3", "p_dcm"], "remove ice bath, warm to room temperature, and stir",
              {"number": 0, "unit": "°C", "end": "r.t.", "raw_text": "0 °C to r.t."}, atmosphere="N2")
set_step_conflict(s_warm, "time", [
    ("12 h (Scheme 6 arrow)", {"number": 12, "unit": "h"}, (5,), "original_extraction"),
    ("overnight (page-9 extracted procedure)", "overnight", (9,), "original_extraction"),
], (5, 9))
r["steps"] = [
    step("cold_addition", (9,), ["p_4a", "p_bh3", "p_dcm"], "add borane complex dropwise to ice-cold solution", {"number": 0, "unit": "°C"}, atmosphere="N2"),
    step("cold_hold", (9,), ["p_4a", "p_bh3", "p_dcm"], "stir at 0 °C", {"number": 0, "unit": "°C"}, {"number": 3, "unit": "h"}, "N2"),
    s_warm,
    step("quench_extract", (9,), ["p_meoh", "p_water", "p_dcm"], "carefully quench with MeOH and H2O; extract with CH2Cl2", stage="workup"),
    step("dry_concentrate", (9,), ["p_mgso4"], "dry combined organic extracts over MgSO4 and evaporate solvent", stage="workup"),
    step("purification_column", (9,), ["p_silica", "p_8"], "purify the crude residue through silica gel chromatography", stage="purification"),
]
r["outcomes"] = [outcome("yield_8", "yield", ["m_8"], (4, 5, 9, 10), {"number": 65, "unit": "%"}, "reported yield",
                          ("purified through silica gel chromatography to yield the pure product ... in 65% yield", "isolated after silica gel chromatography", "original_extraction", None), "4a to 8"),
                 outcome("amount_8", "product_amount", ["m_8"], (9,), {"number": 117.5, "unit": "mg"},
                         basis=("Yield: 117.5 mg, 65%", "reported isolated product mass", "original_extraction", None), scope="product 8")]
r["procedure"]["original_text"] = fact("Scheme 6 eqn (2) and General procedure for the preparation of compound 8", "cold borane addition, 0 °C hold, warm overnight, quench/extract/dry/purify", (5, 9))
r["classification"] = classification(["C-H", "C-H"], (5, 9), "Net carbonyl reduction adds two C–H bonds at the carbonyl carbon; carbonyl oxygen removal and all proton sources remain only partially mapped.")
r["notes"] = [fact("the pure product 7 as a white solid in 65% yield", "Source compound-number inconsistency: procedure says 7; Scheme 6 and characterization identify product 8.", (9,))]
bundle["reactions"].append(r)


bundle["molecules"] = list(molecules.values())
bundle["source_context"] = [
    {"context_id": f"{DOC}:page1_title_abstract", "document_id": DOC, "evidence_ids": page_ev(1),
     "kind": "bibliographic_and_abstract", "note": "Title, authors, DOI and abstract; no concrete experiment on page 1."},
    {"context_id": f"{DOC}:page1_figure1", "document_id": DOC, "evidence_ids": page_ev(1),
     "kind": "non_synthetic_context", "note": "Figure 1 presents active sulfoximine/tryptamine analogues and target scaffold A; it is not a reaction experiment."},
    {"context_id": f"{DOC}:esi_notice", "document_id": DOC, "evidence_ids": page_ev(1),
     "kind": "supplement_notice", "note": "The article says ESI is available, but no SI file was supplied or explicitly associated in this workspace."},
    {"context_id": f"{DOC}:page6_characterization", "document_id": DOC, "evidence_ids": page_ev(6),
     "kind": "characterization", "note": "Product headings, masses, yields, melting points, NMR/IR/HRMS for 4b–4j. Spectral peak lists are context for the corresponding syntheses, not separate reactions."},
    {"context_id": f"{DOC}:page5_characterization_4a", "document_id": DOC, "evidence_ids": page_ev(5),
     "kind": "characterization", "note": "Product 4a mass/yield and NMR/IR/HRMS; not a separate synthesis record."},
    {"context_id": f"{DOC}:page7_characterization", "document_id": DOC, "evidence_ids": page_ev(7),
     "kind": "characterization", "note": "Continuation of 4j and product headings, masses, yields, melting points, NMR/IR/HRMS for 4k–4r. Spectral peak lists are source context for existing Scheme 3 records, not separate reactions."},
    {"context_id": f"{DOC}:page8_characterization_and_gp", "document_id": DOC, "evidence_ids": page_ev(8),
     "kind": "characterization_and_general_procedure", "note": "Product headings, masses, yields, melting points, NMR/IR/HRMS for 4s–4z and the beginning of the general procedure for products 6. The procedure continues on page 9 and is linked to Scheme 5 records."},
    {"context_id": f"{DOC}:page9_procedures_and_characterization", "document_id": DOC, "evidence_ids": page_ev(9),
     "kind": "procedure_and_characterization", "note": "Completion of the products-6 general procedure, characterization headings for 6a–6e, the gram-scale 4a procedure, and the 4a-to-8 reduction procedure. These support existing Scheme 5/6 records rather than duplicate experiments."},
    {"context_id": f"{DOC}:page10_characterization_and_references", "document_id": DOC, "evidence_ids": page_ev(10),
     "kind": "characterization_and_references", "note": "Completion of product 8 characterization followed by conflicts, acknowledgments, and references; no additional synthesis arrow is present."},
    {"context_id": f"{DOC}:page11_references", "document_id": DOC, "evidence_ids": page_ev(11),
     "kind": "references", "note": "References only, including CCDC/ESI and precursor-literature citations; they are not current-paper experiments and do not establish an SI association in this workspace."},
]

(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")

raw = {
    "document_id": DOC,
    "pages_visually_inspected": list(range(1, 12)),
    "bounded_session_pages_visually_inspected": [7, 8, 9, 10, 11],
    "table1": [{"entry": row[0], "yield": row[-1]} for row in table_rows],
    "scheme3": [{"product": x[0], "yield_percent": x[1], "description": x[2]} for x in scope],
    "scheme4": ["direct Ac-protected four-component trial: did not work", "5 + 3a + CO -> 6a: 42%"],
    "scheme5": [{"product": x[0], "result": x[3]} for x in scheme5_rows],
    "scheme6": ["gram-scale 7a -> 4a: 801 mg, 86%", "4a -> 8: 65%"],
    "source_conflicts": [
        "Table/Scheme drawings show step I/II 1 h while captions/footnotes show 2 h.",
        "Table 1 entry 12/Scheme 3 versus page-5 procedure: 2a 0.44 mmol versus 45 μL, 1.2 mmol; THF 1.5 versus 1.0 mL.",
        "Scheme 3 cites Table 1 entry 8 but lists entry-12 DBU conditions.",
        "Scheme 5 PPh3 is 10 mol% in arrow and 20 mol% in caption/page-9 text.",
        "Scheme 5 gives 2–6 h while the pages 8–9 general procedure gives 2–3 h.",
        "Scheme 5 calls the reactant 7 while the pages 8–9 general procedure calls it N-allenyl-2-iodoaniline 5.",
        "Scheme 6 gram-scale arrow versus page-9 text: PPh3 10 versus 20 mol%; toluene 2 versus 10 mL.",
        "The reduction procedure says pure product 7, while Scheme 6 and the characterization heading identify product 8.",
    ],
    "note": "Raw visual transcription; bounded pages 7–11 were viewed at original resolution in this session. Chemistry and drawing-derived identities remain needs_review.",
}
(RUN / "raw_model/visual_transcription.json").write_text(json.dumps(raw, ensure_ascii=False, indent=2) + "\n")
stereo_candidates = {m["molecule_id"]: m["smiles"]["normalized_value"] for m in bundle["molecules"]
                     if m["representation_type"] == "smiles" and m["smiles"]["normalized_value"]}
stereo_diagnostics = {mid: st.stereo_inspect(smi) for mid, smi in stereo_candidates.items()}
(RUN / "raw_model/stereo_diagnostics.json").write_text(json.dumps(stereo_diagnostics, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"])}, ensure_ascii=False))
