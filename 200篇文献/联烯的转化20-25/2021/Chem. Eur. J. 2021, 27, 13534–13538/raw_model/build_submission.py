#!/usr/bin/env python3
"""Build a source-bounded v3.1 submission for doc_08515448f9baa759.

This is a machine extraction artifact.  It intentionally retains unresolved
structures, missing SI-only details, and partial bond-change coverage.
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
DOC = "doc_08515448f9baa759"
bundle = st.read(RUN / "submitted.json")
bundle["documents"] = st.read(RUN / "dataset.json")["documents"]


def ev(page: int, image: bool = True) -> str:
    return f"{DOC}:page{page}:{'image' if image else 'text'}"


def f(raw=None, normalized=None, pages=(), method="original_extraction", derivation=None,
      status=None):
    evidence = [ev(p) for p in pages]
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation,
                   status=status)


def inferred(normalized, pages, derivation, raw=None):
    return f(raw, normalized, pages, method="inference", derivation=derivation)


def derived(normalized, pages, derivation, raw=None):
    return f(raw, normalized, pages, method="structure_derivation", derivation=derivation)


molecules = {}
stereo_checks = []


def ensure_molecule(mid, label, pages, name=None, smiles=None, geometry="unresolved",
                    chirality="unresolved", representation=None, note=None):
    pages = tuple(sorted(set(pages)))
    if mid in molecules:
        occ = molecules[mid]["occurrences"][0]
        ids = set(occ["evidence_ids"]) | {ev(p) for p in pages}
        occ["evidence_ids"] = sorted(ids)
        occ["compound_label"]["evidence_ids"] = sorted(ids)
        return molecules[mid]
    eids = [ev(p) for p in pages]
    m = st.new_molecule(mid, DOC, eids)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, pages)
    m["names"]["original"] = f(name or label, name or label, pages)
    m["names"]["english"] = f(name, name, pages) if name else st.fact(evidence=eids)
    m["cas_numbers"] = st.fact(evidence=eids)
    m["identity_note"] = inferred(
        note or "Identity is source-label bounded; no SI was supplied and the structure drawing was not promoted to a source-verified registry structure.",
        pages,
        "The main-text scheme label and drawing identify the source occurrence, but machine extraction remains needs_review.",
        raw=label,
    )
    m["stereochemistry"]["chirality_status"] = inferred(
        chirality, pages,
        "Status records whether the supplied main-text drawing explicitly specifies stereochemistry; absence of a wedge/name is not treated as achiral or racemic.",
        raw="stereochemistry inspected in cited scheme",
    )
    m["stereochemistry"]["geometry_status"] = inferred(
        geometry, pages,
        "Status records explicit/unspecified alkene or allene geometry in the supplied main-text drawing.",
        raw="geometry inspected in cited scheme",
    )
    if smiles:
        diag = st.stereo_inspect(smiles)
        stereo_checks.append({"molecule_id": mid, "compound_label": label, "input_smiles": smiles, **diag})
        m["representation_type"] = "smiles"
        m["smiles"] = derived(smiles, pages,
            "Candidate SMILES generated from the explicitly named simple reagent/solvent and checked with standard stereo; parsing does not verify source identity.",
            raw={"compound_label": label, "source": "scheme/prose name"})
        m["stereochemistry"]["atom_centers"] = derived(diag["atom_centers"], pages,
            "standard stereo graph diagnostics on candidate SMILES")
        m["stereochemistry"]["double_bonds"] = derived(diag["double_bonds"], pages,
            "standard stereo graph diagnostics on candidate SMILES")
        m["stereochemistry"]["chirality_status"] = inferred(
            "achiral" if not diag["atom_centers"] else "unresolved", pages,
            "Graph diagnostics found no potential atom stereocenter for this simple candidate." if not diag["atom_centers"] else "Potential centers require source assignment.")
        m["stereochemistry"]["geometry_status"] = inferred(
            "not_applicable", pages,
            "No configurational E/Z element was identified in the candidate graph.")
        m["stereochemistry"]["representation_status"] = inferred(
            "not_applicable", pages,
            "No stereogenic element requires encoding in this simple candidate graph.")
        m["stereochemistry"]["isomeric_smiles"] = derived(
            diag["canonical_isomeric_smiles"], pages,
            "Canonical isomeric SMILES returned by standard stereo for the candidate graph.")
        m["stereochemistry"]["identity_key"] = derived(
            diag["stereo_graph_key"], pages,
            "Stereo graph key returned by standard stereo; it is not a sample identity.")
        m["identity_note"] = inferred(
            "Simple-molecule candidate structure derived from the explicit source name; still not source-reviewed after machine merge.",
            pages,
            "Name-to-structure derivation was limited to unambiguous common solvents/reagents.", raw=label)
    else:
        m["representation_type"] = representation or "unresolved"
        m["smiles"] = st.fact(evidence=eids)
        m["stereochemistry"]["representation_status"] = inferred(
            "unresolved", pages,
            "No complete isomeric SMILES is asserted for this source-bounded machine record.")
    molecules[mid] = m
    return m


def ensure_generic(mid, label, pages, fixed, variables):
    m = ensure_molecule(mid, label, pages, name=label, representation="generic_template",
                        geometry="unresolved", chirality="unresolved",
                        note="Generic R/Z-bearing scheme class; it is not a single concrete molecule.")
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": f(fixed, {"description": fixed}, pages),
        "variables": [],
        "instance_mapping": [],
        "reason": f("R/Z variables shown in the scheme", "The source drawing is generic and must not be converted to one concrete SMILES.", pages),
    }
    for var in variables:
        m["generic_template"]["variables"].append({
            "label": var,
            "attachment_points": f(f"scheme-labelled attachment for {var}", f"scheme-labelled attachment for {var}", pages),
            "definition": f(f"{var} is varied in the source scheme", f"{var} is varied in the source scheme", pages),
            "allowed_values": st.fact(evidence=[ev(p) for p in pages]),
            "constraints": st.fact(evidence=[ev(p) for p in pages]),
        })
    return m


def part(pid, mid, role, pages, product_role=None, equiv=None, equiv_ref=None,
         mol_percent=None, stage="reaction"):
    p = st.new_participant(pid, mid, role, [ev(x) for x in pages], stage=stage)
    if product_role:
        p["product_role"] = f(product_role, product_role, pages)
    if equiv is not None:
        p["charges"]["equivalents"] = f(f"{equiv} equiv", {"number": equiv, "unit": "equiv"}, pages)
        p["equivalence_reference"] = inferred(equiv_ref, pages, "Equivalence denominator is the explicitly labelled substrate reference.")
    if mol_percent is not None:
        p["charges"]["mol_percent"] = f(f"{mol_percent} mol%", {"number": mol_percent, "unit": "mol%"}, pages)
        p["equivalence_reference"] = inferred(equiv_ref, pages, "Catalyst loading denominator is the reaction substrate.")
    return p


TAX = {
    "C-H": 1, "C-C": 2, "C-N": 5, "DB-C-C": 13, "DB-C-O": 14,
}


def classification(labels, pages, unresolved=False, not_applicable=False, extra_limitations=()):
    c = st.new_reaction(f"{DOC}:_tmp", DOC, [ev(p) for p in pages])["classification"]
    c["mapped_reaction_smiles"] = st.fact(evidence=[ev(p) for p in pages])
    c["bond_changes"] = []
    c["model_votes"] = []
    c["limitations"] = [
        "No atom-mapped reaction SMILES is asserted because the supplied main text does not provide machine-verifiable atom identities.",
        "Implicit-H and omitted coproduct changes may be incomplete; labels describe only visually supported net formed/promoted bonds.",
        *extra_limitations,
    ]
    if not_applicable:
        c["labels"] = inferred([], pages, "No observed product transformation; the taxonomy is not applied.", raw="no reaction/complex mixture")
        c["coverage"] = inferred("not_applicable", pages, "No single observed product graph is available for bond-formation classification.")
        return c
    if unresolved:
        c["labels"] = inferred([], pages, "The net bond changes could not be assigned safely from the supplied figure.")
        c["coverage"] = inferred("unresolved", pages, "Source drawing or overall one-pot sequence is insufficient for safe net mapping.")
        return c
    change_ids = []
    counts = {}
    for i, code in enumerate(labels, 1):
        counts[code] = counts.get(code, 0) + 1
        cid = f"b{i}"
        if code == "C-C":
            obj = {"change_id": cid, "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
                   "before_order": 0, "after_order": 1, "change_type": "formed"}
        elif code == "C-N":
            obj = {"change_id": cid, "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "N"}],
                   "before_order": 0, "after_order": 1, "change_type": "formed"}
        elif code == "DB-C-C":
            obj = {"change_id": cid, "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
                   "before_order": 1, "after_order": 2, "change_type": "order_increased"}
        elif code == "DB-C-O":
            obj = {"change_id": cid, "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}],
                   "before_order": 0, "after_order": 2, "change_type": "formed"}
        elif code == "C-H":
            obj = {"change_id": cid, "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}],
                   "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1}
        else:
            raise ValueError(code)
        c["bond_changes"].append(derived(obj, pages,
            "Net reactant/product drawing comparison; atom identities remain unresolved."))
        change_ids.append((code, cid))
    ordered = []
    for code in dict.fromkeys(labels):
        ordered.append({"category_number": TAX[code], "code": code,
                        "supporting_change_ids": [cid for ccode, cid in change_ids if ccode == code]})
    c["labels"] = derived(ordered, pages,
        "Taxonomy labels calculated from the explicitly listed partial net bond changes.")
    c["coverage"] = inferred("partial", pages,
        "Only the net formed/promoted bonds that can be read safely from the main-text drawings are included.")
    return c


def step(step_id, pages, operation, temperature=None, time=None, participant_ids=()):
    s = st.new_step(step_id, [ev(p) for p in pages])
    s["participant_ids"] = list(participant_ids)
    s["operation"] = f(operation, operation, pages)
    if temperature is not None:
        raw, norm = temperature
        s["temperature"] = f(raw, norm, pages)
    if time is not None:
        raw, norm = time
        s["time"] = f(raw, norm, pages)
    return s


def outcome(outcome_id, kind, targets, pages, raw=None, normalized=None,
            isolated=False, scope=None, basis_raw=None):
    o = st.new_outcome(outcome_id, kind, targets, [ev(p) for p in pages])
    o["value"] = f(raw, normalized, pages)
    if targets:
        o["product_role"] = f("product shown for this entry", "main_product" if kind in ("yield", "product_amount") else "reported_target", pages)
    if isolated:
        o["yield_type"] = f("isolated yield", "isolated", pages)
        o["measurement_basis"] = f("[isolated yield, t]", "isolated", pages)
    elif basis_raw:
        o["measurement_basis"] = f(basis_raw, basis_raw, pages)
    if scope is not None:
        o["scope"] = f(scope, scope, pages)
    return o


reactions = []


def reaction(rid, pages, record_type, categories, location, name, transform,
             participants, steps, outcomes, classif, general=(), notes=(), unresolved=()):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, [ev(p) for p in pages])
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, pages)
    r["location"] = f(location, location, pages)
    r["general_scheme_ids"] = list(general)
    r["reaction_name"] = f(name, name, pages) if name else st.fact(evidence=[ev(p) for p in pages])
    r["transformation_type"] = inferred(transform, pages, "Transformation type assigned from the source caption and net structure drawings.")
    r["participants"] = participants
    r["steps"] = steps
    r["outcomes"] = outcomes
    r["classification"] = classif
    if steps:
        r["procedure"]["original_text"] = f("; ".join(str(s["operation"]["raw_value"]) for s in steps),
                                             "; ".join(str(s["operation"]["normalized_value"]) for s in steps), pages)
    r["notes"] = [f(n, n, pages) for n in notes]
    r["review"]["unresolved"] = list(unresolved) or [
        "No explicitly associated SI was supplied; charges, workup, purification, and characterization details absent from the main-text figure remain unreported.",
        "Structures and atom mapping remain machine candidates/unresolved; chemistry correctness is not certified.",
    ]
    r["completeness"] = "partial"
    reactions.append(r)
    return r


# Common, unambiguous simple substances.  These are the only concrete SMILES
# promoted; all paper-specific structures remain source-labelled/unresolved.
simple = [
    ("m_etoac", "EtOAc", (1, 2, 3, 4), "ethyl acetate", "CCOC(C)=O"),
    ("m_water", "H2O", (1, 2, 3, 4), "water", "O"),
    ("m_ethanol", "EtOH", (2,), "ethanol", "CCO"),
    ("m_mecn", "MeCN", (2,), "acetonitrile", "CC#N"),
    ("m_toluene", "toluene", (2,), "toluene", "Cc1ccccc1"),
    ("m_dce", "DCE", (2,), "1,2-dichloroethane", "ClCCCl"),
    ("m_acetone", "acetone", (2,), "acetone", "CC(C)=O"),
    ("m_tempo", "TEMPO", (2,), "2,2,6,6-tetramethylpiperidine 1-oxyl", "CC1(C)CCCC(C)(C)N1[O]"),
    ("m_tbuamine", "t-BuNH2", (3, 4), "tert-butylamine", "CC(C)(C)N"),
    ("m_tfa", "TFA", (4,), "trifluoroacetic acid", "O=C(O)C(F)(F)F"),
    ("m_dcm", "DCM", (4,), "dichloromethane", "ClCCl"),
    ("m_h2", "H2", (4,), "hydrogen", "[H][H]"),
    ("m_ether", "Et2O", (4,), "diethyl ether", "CCOCC"),
    ("m_ozone", "O3", (4,), "ozone", "[O-][O+]=O"),
    ("m_iproh", "i-PrOH", (4,), "2-propanol", "CC(O)C"),
]
for args in simple:
    ensure_molecule(*args)

# Source-labelled paper compounds.
for label, pages in {
    "1": (1, 2, 3, 4), "[D2]-1": (3,),
    **{f"2{x}": (2,) for x in "abcdefghijklmnopq"},
    **{f"3{x}": (2,) for x in "abcdefghijklmnopq"},
    "3p-Me": (2,), "4a": (3,), "4b": (3,), "4c": (3,), "4d": (3,),
    "3r": (3,), "3s": (3,), "5": (3,), "6": (3,),
    "7a": (3,), "7h": (3,), "7j": (3,), "7n": (3,),
    **{f"8{x}": (3,) for x in "abcdefgh"},
    **{f"9{x}": (3,) for x in "abcdefgh"},
    "[D2]-9b": (3,), "10c": (4,), "7fa": (4,), "11h": (4,), "12c": (4,),
    "13a": (4,), "14a": (4,), "15b": (4,), "16b": (4,),
}.items():
    mid = ("m_" + label.replace("[", "").replace("]", "").replace("-", "_")).lower()
    geom = "mixture" if label in ("2p", "3p") else "specified" if label == "10c" else "unspecified"
    chir = "unspecified"
    ensure_molecule(mid, label, pages, name=label, geometry=geom, chirality=chir)

# Explicit stereo/mixture statements retained without inventing isomeric SMILES.
for label in ("2p", "3p"):
    m = molecules["m_" + label]
    m["stereochemistry"]["mixture"] = f("58:42 ratio Z:E; same ratio in precursor 2p",
        {"metric": "Z:E", "ratio": "58:42"}, (2,))
    m["stereochemistry"]["geometry_status"] = inferred("mixture", (2,),
        "Scheme 3 footnote explicitly reports a Z:E mixture.")
molecules["m_10c"]["stereochemistry"]["double_bonds"] = f(
    "(E)-4-styryl-3-(triflyl)quinoline 10c", [{"group": "styryl", "configuration": "E"}], (3, 4))
molecules["m_10c"]["stereochemistry"]["geometry_status"] = inferred(
    "specified", (3, 4), "The prose explicitly names product 10c as (E).")
molecules["m_10c"]["stereochemistry"]["representation_status"] = inferred(
    "not_encoded", (4,), "The explicit E assignment is retained textually because no complete isomeric SMILES is asserted.")

# Other reagents/catalysts are source-labelled and deliberately unresolved.
for mid, label, pages, name in [
    ("m_1", "1", (1, 2, 3, 4), "Yanai's reagent / 2-fluoropyridinium betaine 1"),
    ("m_d2_1", "[D2]-1", (3,), "deuterated Yanai's reagent [D2]-1"),
    ("m_p_tol_boronic", "ArB(OH)2, Ar = 4-MeC6H4", (4,), "4-methylphenylboronic acid"),
    ("m_k2co3", "K2CO3", (4,), "potassium carbonate"),
    ("m_pdpph3", "Pd(PPh3)4", (4,), "tetrakis(triphenylphosphine)palladium"),
    ("m_pdc", "Pd/C", (4,), "palladium on carbon"),
    ("m_lialh4", "LiAlH4", (4,), "lithium aluminium hydride"),
    ("m_silica", "silica gel", (4,), "silica gel"),
    ("m_ipraucl", "IPrAuCl", (4,), "IPrAuCl"),
    ("m_agotf", "AgOTf", (4,), "silver triflate"),
].__iter__():
    ensure_molecule(mid, label, pages, name=name)

# Generic scheme classes.
ensure_generic("m_g_2", "2 / generic alkenol", (1, 2, 4),
               "1-(2-aminophenyl)prop-2-en-1-ol class shown with N-protecting group and ring/alkene substituents", ["Ar/R1", "R2", "Z"])
ensure_generic("m_g_3", "3 / generic bis(triflyl) THQ", (1, 2, 4),
               "bis(triflyl)tetrahydroquinoline class bearing the retained unsaturated substituent", ["Ar/R1", "R2", "Z"])
ensure_generic("m_g_4", "4a-c / generic allylic alcohol", (3,),
               "isomeric allylic alcohol class shown in Scheme 4a", ["R1", "R2", "R3"])
ensure_generic("m_g_7", "7 / generic triflylated DHQ", (3,),
               "triflylated 1,2-dihydroquinoline class shown in Scheme 5", ["R1", "R2", "Z"])
ensure_generic("m_g_8", "8 / generic allenol", (3, 4),
               "1-(2-aminophenyl)buta-2,3-dien-1-ol class shown in Scheme 6", ["R1", "Z"])
ensure_generic("m_g_9", "9 / generic allenyl bis(triflyl) THQ", (3, 4),
               "allenyl bis(triflyl)tetrahydroquinoline class shown in Scheme 6", ["R1", "Z"])
ensure_generic("m_g_alkynol", "functionalized anilide bearing an alkynol", (3,),
               "(2-aminophenyl)methanol class in which the alkene is replaced by an alkyne", ["substituents"])
ensure_generic("m_int_i", "INT-I", (4,), "aza-ortho-quinone methide proposed intermediate", ["Ar", "R", "Z"])
ensure_generic("m_int_ii", "INT-II", (4,), "zwitterionic proposed intermediate", ["Ar", "R", "Z"])
ensure_molecule("m_A", "A", (4,), name="bis((trifluoromethyl)sulfonyl)methane A")

# Generic overview and source scheme records.
g1 = reaction("scheme1d_general", (1,), "generic_scheme", ["other"],
    {"page": 1, "scheme": "1d"}, "Incidental discovery", "formal [4+2] heterocyclization",
    [part("g1_s", "m_g_2", "reactant", (1,)), part("g1_r", "m_1", "reagent", (1,)),
     part("g1_sol", "m_etoac", "solvent", (1,)), part("g1_p", "m_g_3", "product", (1,), "main_product")],
    [step("g1_reaction", (1,), "1 as source of Tf2C=CH2; EtOAc, rt; without catalyst; -H2O, -2F-pyridine (volatile)",
          ("rt", "rt"), None, ("g1_s", "g1_r", "g1_sol", "g1_p"))], [],
    classification(["C-C", "C-N"], (1,)),
    notes=("Generic overview only; no entry-specific yield is assigned.",))

g3 = reaction("scheme3_general", (2,), "generic_scheme", ["substrate_scope"],
    {"page": 2, "scheme": "3"}, "Controlled synthesis of bis(triflyl)ated tetrahydroquinolines 3a-q",
    "formal [4+2] heterocyclization",
    [part("g3_s", "m_g_2", "reactant", (2,)), part("g3_r", "m_1", "reagent", (2,)),
     part("g3_sol", "m_etoac", "solvent", (2,)), part("g3_p", "m_g_3", "product", (2,), "main_product")],
    [step("g3_reaction", (2,), "1; EtOAc, rt; -2F-pyridine, -H2O", ("rt", "rt"), None,
          ("g3_s", "g3_r", "g3_sol", "g3_p"))], [], classification(["C-C", "C-N"], (2,)),
    general=(g1["reaction_id"],), notes=("Entry-specific yields and times are recorded separately.",))

g4 = reaction("scheme4a_general", (2, 3), "generic_scheme", ["substrate_scope"],
    {"page": 3, "scheme": "4a"}, "Synthesis from isomeric allylic alcohols 4a-c", "formal [4+2] heterocyclization",
    [part("g4_s", "m_g_4", "reactant", (3,)), part("g4_r", "m_1", "reagent", (3,)),
     part("g4_sol", "m_etoac", "solvent", (3,)), part("g4_p", "m_g_3", "product", (3,), "main_product")],
    [step("g4_reaction", (3,), "1; EtOAc, rt; -2F-pyridine, -H2O", ("rt", "rt"), None,
          ("g4_s", "g4_r", "g4_sol", "g4_p"))], [], classification(["C-C", "C-N"], (3,)),
    general=(g1["reaction_id"],))

g5 = reaction("scheme5_general", (2, 3), "generic_scheme", ["substrate_scope", "derivatization"],
    {"page": 3, "scheme": "5"}, "Controlled synthesis of triflylated 1,2-dihydroquinolines 7",
    "one-pot heterocyclization and hydrodetriflylation",
    [part("g5_s", "m_g_2", "reactant", (3,)), part("g5_r", "m_1", "reagent", (3,)),
     part("g5_sol", "m_etoac", "solvent", (3,)), part("g5_b", "m_tbuamine", "reagent", (3,)),
     part("g5_p", "m_g_7", "product", (3,), "main_product")],
    [step("g5_1", (3,), "1, EtOAc, rt", ("rt", "rt"), None, ("g5_s", "g5_r", "g5_sol")),
     step("g5_2", (3,), "then t-BuNH2, 0 °C; -2F-pyridine, -H2O, -TfH", ("0 °C", {"number": 0, "unit": "°C"}), None,
          ("g5_b", "g5_p"))], [], classification(["C-C", "C-N", "DB-C-C"], (3,)),
    general=(g1["reaction_id"],), notes=("The figure reports entry-associated total time but does not allocate time between the two stages.",))

g6 = reaction("scheme6_general", (3,), "generic_scheme", ["substrate_scope"],
    {"page": 3, "scheme": "6"}, "Controlled synthesis of bis(triflyl)ated tetrahydroquinolines 9a-h",
    "formal [4+2] heterocyclization of allenols",
    [part("g6_s", "m_g_8", "reactant", (3,)), part("g6_r", "m_1", "reagent", (3,)),
     part("g6_sol", "m_etoac", "solvent", (3,)), part("g6_p", "m_g_9", "product", (3,), "main_product")],
    [step("g6_reaction", (3,), "1 or [D2]-1; EtOAc, 90 °C, sealed tube; -2F-pyridine, -H2O",
          ("90 °C", {"number": 90, "unit": "°C"}), None, ("g6_s", "g6_r", "g6_sol", "g6_p"))], [],
    classification(["C-C", "C-N"], (3,)), general=(g1["reaction_id"],),
    notes=("Entry-specific reagent isotope, yield and time are recorded separately.",))

# Scheme 2 optimization: one row per condition.
opt_rows = [
    ("water_90c_24h", "m_water", "H2O, 90 °C, 24 h", ("90 °C", {"number": 90, "unit": "°C"}), ("24 h", {"number": 24, "unit": "h"}), 0, 1.0, True),
    ("ethanol_rt_10min", "m_ethanol", "EtOH, rt, 10 min", ("rt", "rt"), ("10 min", {"number": 10, "unit": "min"}), 23, 1.0, False),
    ("mecn_rt_30min", "m_mecn", "MeCN, rt, 30 min", ("rt", "rt"), ("30 min", {"number": 30, "unit": "min"}), 72, 1.0, False),
    ("toluene_40c_1h", "m_toluene", "toluene, 40 °C, 1 h", ("40 °C", {"number": 40, "unit": "°C"}), ("1 h", {"number": 1, "unit": "h"}), 69, 1.0, False),
    ("dce_rt_30min", "m_dce", "DCE, rt, 30 min", ("rt", "rt"), ("30 min", {"number": 30, "unit": "min"}), 76, 1.0, False),
    ("etoac_rt_30min", "m_etoac", "EtOAc, rt, 30 min", ("rt", "rt"), ("30 min", {"number": 30, "unit": "min"}), 74, 1.0, False),
    ("acetone_rt_10min_2equiv", "m_acetone", "acetone, rt, 10 min; 1 (2.0 equiv)", ("rt", "rt"), ("10 min", {"number": 10, "unit": "min"}), 70, 2.0, False),
]
for key, solvent, cond, temp, tm, yld, reagent_eq, no_rxn in opt_rows:
    ps = [part("substrate", "m_2a", "reactant", (2,), equiv=1.0, equiv_ref="m_2a"),
          part("reagent1", "m_1", "reagent", (2,), equiv=reagent_eq, equiv_ref="m_2a"),
          part("solvent", solvent, "solvent", (2,))]
    outs = []
    if no_rxn:
        outs.append(outcome("no_reaction", "no_reaction", [], (2,), "3a (0%); reaction did not proceed in hot water", True,
                            scope={"scheme": 2, "row": key}))
        cls = classification([], (2,), not_applicable=True)
    else:
        ps.append(part("product", "m_3a", "product", (2,), "main_product"))
        outs.append(outcome("yield_3a", "yield", ["m_3a"], (2,), f"3a ({yld}%)",
                            {"number": yld, "unit": "%"}, scope={"scheme": 2, "row": key}))
        cls = classification(["C-C", "C-N"], (2,))
    reaction(f"scheme2_{key}", (2,), "specific_experiment", ["optimization"],
        {"page": 2, "scheme": 2, "row": key}, "Impact of solvent and temperature on preparation of 3a",
        "formal [4+2] heterocyclization", ps,
        [step("reaction", (2,), cond, temp, tm, tuple(p["participant_id"] for p in ps))], outs, cls,
        general=(g1["reaction_id"],),
        unresolved=("Yield measurement basis is not stated for Scheme 2; it remains unreported.",
                    "No SI was supplied; charges beyond displayed equivalents and workup/purification remain unreported."))

# TEMPO mechanistic control.
reaction("tempo_control", (2,), "specific_experiment", ["control", "mechanistic_study"],
    {"page": 2, "prose": "TEMPO control"}, "Reaction in the presence of TEMPO", "radical-trap control of heterocyclization",
    [part("substrate", "m_2a", "reactant", (2,)), part("reagent1", "m_1", "reagent", (2,)),
     part("tempo", "m_tempo", "additive", (2,)), part("product", "m_3a", "product", (2,), "main_product")],
    [step("reaction", (2,), "conducting the reaction in the presence of TEMPO did not trim the yield", None, None,
          ("substrate", "reagent1", "tempo", "product"))],
    [outcome("yield_effect", "other", ["m_3a"], (2,), "did not trim the yield", "yield not decreased")],
    classification(["C-C", "C-N"], (2,)), general=(g1["reaction_id"],),
    unresolved=("TEMPO quantity, exact conditions, numerical yield and measurement basis are not reported in the supplied main text.",))

# Scheme 3 scope entries.
scope3 = [
    ("a", 74, "30 min", False), ("b", 67, "30 min", False), ("c", 90, "30 min", False),
    ("d", 58, "30 min", False), ("e", 78, "30 min", False), ("f", 81, "30 min", False),
    ("g", 83, "30 min", False), ("h", 68, "30 min", False), ("i", 78, "30 min", False),
    ("j", "quantitative", "4 h", False), ("k", 81, "30 min", False), ("l", "quantitative", "4 h", False),
    ("m", 75, "30 min", False), ("n", 61, "2 h", True), ("o", 48, "2 h", True),
    ("p", 60, "30 min", False), ("q", 84, "30 min", False),
]
for letter, yld, tmraw, sealed90 in scope3:
    sm, pm = f"m_2{letter}", f"m_3{letter}"
    ps = [part("substrate", sm, "reactant", (2,)), part("reagent1", "m_1", "reagent", (2,)),
          part("solvent", "m_etoac", "solvent", (2,)), part("product", pm, "product", (2,), "main_product")]
    tnum, tunit = tmraw.split()
    temp = ("90 °C", {"number": 90, "unit": "°C"}) if sealed90 else ("rt", "rt")
    op = "1; EtOAc, 90 °C, sealed tube" if sealed90 else "1; EtOAc, rt"
    ynorm = yld if isinstance(yld, str) else {"number": yld, "unit": "%"}
    yraw = f"3{letter} ({yld if isinstance(yld, str) else str(yld)+'%'}, {tmraw})"
    outs = [outcome(f"yield_3{letter}", "yield", [pm], (2,), yraw, ynorm, isolated=True,
                    scope={"scheme": 3, "entry": f"3{letter}"})]
    if letter == "p":
        sel = outcome("selectivity_3p", "selectivity", [pm], (2,),
                      "58:42 ratio Z:E; same ratio in precursor 2p", {"ratio": "58:42", "order": "Z:E"},
                      scope={"metric": "Z:E", "order": ["Z", "E"], "precursor_ratio": "same"})
        outs.append(sel)
    reaction(f"scheme3_2{letter}_to_3{letter}", (2,), "specific_experiment", ["substrate_scope"],
        {"page": 2, "scheme": 3, "entry": f"3{letter}"}, "Controlled synthesis of bis(triflyl)ated tetrahydroquinolines",
        "formal [4+2] heterocyclization", ps,
        [step("reaction", (2,), op, temp, (tmraw, {"number": int(tnum), "unit": "min" if tunit == "min" else "h"}),
              tuple(p["participant_id"] for p in ps))], outs, classification(["C-C", "C-N"], (2,)),
        general=(g3["reaction_id"],),
        unresolved=("No SI was supplied; displayed scope entry supports product, isolated yield, time and indicated temperature only.",
                    "Complete structures/atom mapping and stereochemical assignments other than the explicit 3p Z:E ratio remain unresolved."))

# Nonproductive precursor 3p-Me.
ensure_molecule("m_3p_me", "3p-Me", (2,), name="non-productive precursor 3p-Me", geometry="unspecified", chirality="unspecified")
reaction("scheme3_3p_me_nonproductive", (2,), "specific_experiment", ["substrate_scope", "control"],
    {"page": 2, "scheme": 3, "entry": "3p-Me"}, "Protocol applied to 2,3-disubstituted alkenol 3p-Me",
    "failed annulation / complex mixture",
    [part("substrate", "m_3p_me", "reactant", (2,)), part("reagent1", "m_1", "reagent", (2,)),
     part("solvent", "m_etoac", "solvent", (2,))],
    [step("reaction", (2,), "protocol applied to 3p-Me; complicated mixture formed", ("rt", "rt"), None,
          ("substrate", "reagent1", "solvent"))],
    [outcome("complex_mixture", "decomposition", [], (2,), "complicated mixture; desired azacycle could not be attained", "complex mixture")],
    classification([], (2,), not_applicable=True), general=(g3["reaction_id"],),
    unresolved=("Exact time, mixture composition, conversion and measurement method are unreported.",))

# Scheme 4a scope and 4b alternative product.
for src, prod, yld, tmraw in [("4a", "3a", 72, "30 min"), ("4b", "3r", 57, "2 h"),
                              ("4c", "3s", 88, "30 min"), ("4d", "3c", 74, "30 min")]:
    tnum, tunit = tmraw.split()
    ps = [part("substrate", f"m_{src}", "reactant", (3,)), part("reagent1", "m_1", "reagent", (3,)),
          part("solvent", "m_etoac", "solvent", (3,)), part("product", f"m_{prod}", "product", (3,), "main_product")]
    reaction(f"scheme4_{src}_to_{prod}", (2, 3), "specific_experiment", ["substrate_scope"],
        {"page": 3, "scheme": "4a", "entry": prod}, "Synthesis from isomeric allylic alcohol",
        "formal [4+2] heterocyclization", ps,
        [step("reaction", (3,), "1; EtOAc, rt; -2F-pyridine, -H2O", ("rt", "rt"),
              (tmraw, {"number": int(tnum), "unit": "min" if tunit == "min" else "h"}),
              tuple(p["participant_id"] for p in ps))],
        [outcome(f"yield_{prod}", "yield", [f"m_{prod}"], (3,), f"{prod} ({yld}%, {tmraw})",
                 {"number": yld, "unit": "%"}, scope={"scheme": "4a", "entry": prod})],
        classification(["C-C", "C-N"], (3,)), general=(g4["reaction_id"],),
        unresolved=("Yield measurement basis and SI-level procedure are unreported in the supplied main text.",))

ps = [part("substrate", "m_5", "reactant", (3,)), part("reagent1", "m_1", "reagent", (3,)),
      part("solvent", "m_etoac", "solvent", (3,)), part("product", "m_6", "product", (3,), "exclusive product")]
reaction("scheme4_5_to_6", (2, 3), "specific_experiment", ["substrate_scope"],
    {"page": 3, "scheme": "4b"}, "Reaction of tertiary allylic alcohol 5", "bis(triflyl)ethylation/cyclization",
    ps, [step("reaction", (3,), "1; EtOAc, 60 °C, 30 min; -2F-pyridine, -H2O",
              ("60 °C", {"number": 60, "unit": "°C"}), ("30 min", {"number": 30, "unit": "min"}),
              tuple(p["participant_id"] for p in ps))],
    [outcome("yield_6", "yield", ["m_6"], (3,), "6 (88%)", {"number": 88, "unit": "%"},
             scope={"scheme": "4b", "product_statement": "exclusive product"})],
    classification(["C-C", "C-N"], (3,)), general=(g1["reaction_id"],),
    unresolved=("Yield measurement basis and complete net atom mapping are unreported.",))

# Scheme 5 one-pot entries.
for letter, yld, tmraw in [("a", 91, "35 min"), ("h", 95, "35 min"), ("j", 87, "4 h"), ("n", 91, "2 h")]:
    sm, pm = f"m_2{letter}", f"m_7{letter}"
    ps = [part("substrate", sm, "reactant", (3,)), part("reagent1", "m_1", "reagent", (3,)),
          part("solvent", "m_etoac", "solvent", (3,)), part("base", "m_tbuamine", "reagent", (3,)),
          part("product", pm, "product", (3,), "main_product")]
    tnum, tunit = tmraw.split()
    reaction(f"scheme5_2{letter}_to_7{letter}", (2, 3), "specific_experiment", ["substrate_scope", "derivatization"],
        {"page": 3, "scheme": 5, "entry": f"7{letter}"}, "Controlled synthesis of triflylated 1,2-dihydroquinoline",
        "one-pot heterocyclization and hydrodetriflylation", ps,
        [step("annulation", (3,), "1, EtOAc, rt", ("rt", "rt"), None, ("substrate", "reagent1", "solvent")),
         step("base_treatment", (3,), "then t-BuNH2, 0 °C; reported product-associated total time shown with entry",
              ("0 °C", {"number": 0, "unit": "°C"}),
              (tmraw, {"number": int(tnum), "unit": "min" if tunit == "min" else "h"}), ("base", "product"))],
        [outcome(f"yield_7{letter}", "yield", [pm], (3,), f"7{letter} ({yld}%, {tmraw})",
                 {"number": yld, "unit": "%"}, isolated=True, scope={"scheme": 5, "entry": f"7{letter}", "overall": True})],
        classification(["C-C", "C-N", "DB-C-C"], (3,)), general=(g5["reaction_id"],),
        unresolved=("The displayed time is associated with the overall entry; allocation between rt annulation and 0 °C base treatment is not stated.",
                    "No SI was supplied; complete charges, workup and purification remain unreported."))

# Alkynol negative-scope observation.
reaction("alkynol_complex_mixture", (3,), "specific_experiment", ["control", "substrate_scope"],
    {"page": 3, "prose": "alkynol instead of alkenol"}, "Reaction of functionalized anilides bearing alkynol substituents",
    "failed/off-target reaction yielding complex mixtures",
    [part("substrate", "m_g_alkynol", "reactant", (3,)), part("reagent1", "m_1", "reagent", (3,))],
    [step("reaction", (3,), "functionalized anilides with alkynol instead of alkenol substituents gave complex mixtures", None, None,
          ("substrate", "reagent1"))],
    [outcome("complex_mixture", "decomposition", [], (3,), "complex mixtures of products", "complex mixture")],
    classification([], (3,), not_applicable=True),
    unresolved=("Specific substrate labels, charges, temperature, time, conversion and mixture composition are not reported in the supplied main text.",))

# Scheme 6 allenol scope.
for letter, yld, tmraw in [("a", 82, "30 min"), ("b", 74, "30 min"), ("c", 64, "40 min"),
                           ("d", 52, "30 min"), ("e", 54, "30 min"), ("f", 63, "1 h"),
                           ("g", 69, "30 min"), ("h", 68, "30 min")]:
    sm, pm = f"m_8{letter}", f"m_9{letter}"
    ps = [part("substrate", sm, "reactant", (3,)), part("reagent1", "m_1", "reagent", (3,)),
          part("solvent", "m_etoac", "solvent", (3,)), part("product", pm, "product", (3,), "main_product")]
    tnum, tunit = tmraw.split()
    reaction(f"scheme6_8{letter}_to_9{letter}", (3,), "specific_experiment", ["substrate_scope"],
        {"page": 3, "scheme": 6, "entry": f"9{letter}"}, "Controlled synthesis of bis(triflyl)ated tetrahydroquinoline from allenol",
        "formal [4+2] heterocyclization of allenol", ps,
        [step("reaction", (3,), "1; EtOAc, 90 °C, sealed tube; -2F-pyridine, -H2O",
              ("90 °C", {"number": 90, "unit": "°C"}),
              (tmraw, {"number": int(tnum), "unit": "min" if tunit == "min" else "h"}), tuple(p["participant_id"] for p in ps))],
        [outcome(f"yield_9{letter}", "yield", [pm], (3,), f"9{letter} ({yld}%, {tmraw})",
                 {"number": yld, "unit": "%"}, isolated=True, scope={"scheme": 6, "entry": f"9{letter}"})],
        classification(["C-C", "C-N"], (3,)), general=(g6["reaction_id"],),
        unresolved=("Allene axial stereochemistry is not assigned in the source; terminal allene drawings are retained without inventing an axial descriptor.",
                    "No SI was supplied; complete structures/charges/workup/purification remain unreported."))

ps = [part("substrate", "m_8b", "reactant", (3,)), part("reagent1", "m_d2_1", "reagent", (3,)),
      part("solvent", "m_etoac", "solvent", (3,)), part("product", "m_d2_9b", "product", (3,), "main_product")]
reaction("scheme6_d2_8b_to_d2_9b", (3,), "specific_experiment", ["substrate_scope", "mechanistic_study"],
    {"page": 3, "scheme": 6, "entry": "[D2]-9b"}, "Reaction of 8b with deuterated Yanai's reagent [D2]-1",
    "isotope-labelled formal [4+2] heterocyclization", ps,
    [step("reaction", (3,), "[D2]-1; EtOAc, 90 °C, sealed tube, 30 min; -2F-pyridine, -H2O",
          ("90 °C", {"number": 90, "unit": "°C"}), ("30 min", {"number": 30, "unit": "min"}),
          tuple(p["participant_id"] for p in ps))],
    [outcome("yield_d2_9b", "yield", ["m_d2_9b"], (3,), "[D2]-9b (77%, 30 min)",
             {"number": 77, "unit": "%"}, isolated=True, scope={"scheme": 6, "isotope": "D2"})],
    classification(["C-C", "C-N"], (3,), extra_limitations=("C-D bonds already present in [D2]-1 are not classified as newly formed C-H bonds.",)),
    general=(g6["reaction_id"],),
    unresolved=("Exact deuterium placement is retained from the drawing but no complete isotope-aware SMILES/atom mapping is asserted.",))

# Scheme 7 product elaborations.  Each arrow/overall sequence is a separate record.
def scheme7_record(rid, substrate, product, participants_extra, steps_list, yld, transform, labels, unresolved_cls=False):
    ps = [part("substrate", substrate, "reactant", (4,)), *participants_extra,
          part("product", product, "product", (4,), "main_product")]
    return reaction(rid, (3, 4), "specific_experiment", ["derivatization"],
        {"page": 4, "scheme": 7, "entry": product.replace("m_", "")}, "Product elaboration", transform, ps,
        steps_list,
        [outcome("yield", "yield", [product], (4,), f"{product.replace('m_', '').replace('_', '-')} ({yld}%)",
                 {"number": yld, "unit": "%"}, scope={"scheme": 7, "overall_sequence": True})],
        classification(labels, (4,), unresolved=unresolved_cls), general=(g1["reaction_id"],),
        unresolved=("Scheme 7 reports an overall one-pot sequence from the indicated precursor; yield basis and stage-specific material balances are not stated.",
                    "No SI was supplied; complete charges, workup, purification and atom mapping remain unreported."))

scheme7_record("scheme7_2c_to_10c", "m_2c", "m_10c",
    [part("reagent1", "m_1", "reagent", (4,)), part("solvent1", "m_etoac", "solvent", (4,)),
     part("base", "m_tbuamine", "reagent", (4,)), part("acid", "m_tfa", "reagent", (4,)), part("solvent2", "m_dcm", "solvent", (4,))],
    [step("annulation", (4,), "1, EtOAc, rt, 30 min", ("rt", "rt"), ("30 min", {"number": 30, "unit": "min"}), ("substrate", "reagent1", "solvent1")),
     step("base", (4,), "t-BuNH2, EtOAc, 0 °C, 5 min", ("0 °C", {"number": 0, "unit": "°C"}), ("5 min", {"number": 5, "unit": "min"}), ("base", "solvent1")),
     step("acid", (4,), "TFA, DCM, 40 °C, 3 h", ("40 °C", {"number": 40, "unit": "°C"}), ("3 h", {"number": 3, "unit": "h"}), ("acid", "solvent2", "product"))],
    67, "one-pot annulation, detrifylation and aromatization", ["C-C", "C-N", "DB-C-C"])

scheme7_record("scheme7_2f_to_7fa", "m_2f", "m_7fa",
    [part("reagent1", "m_1", "reagent", (4,)), part("solvent", "m_etoac", "solvent", (4,)),
     part("boronic", "m_p_tol_boronic", "reagent", (4,)), part("base", "m_k2co3", "reagent", (4,)),
     part("catalyst", "m_pdpph3", "catalyst", (4,), mol_percent=5, equiv_ref="m_2f"), part("water", "m_water", "solvent", (4,))],
    [step("annulation", (4,), "1, EtOAc, rt, 30 min", ("rt", "rt"), ("30 min", {"number": 30, "unit": "min"}), ("substrate", "reagent1", "solvent")),
     step("coupling", (4,), "ArB(OH)2, K2CO3, 5 mol% Pd(PPh3)4, 1,4-dioxane:H2O (2:1), 60 °C, 2 h",
          ("60 °C", {"number": 60, "unit": "°C"}), ("2 h", {"number": 2, "unit": "h"}),
          ("boronic", "base", "catalyst", "water", "product"))],
    51, "one-pot annulation/Suzuki-Miyaura coupling with TfH release", ["C-C", "C-N", "DB-C-C", "C-C"])

scheme7_record("scheme7_2h_to_11h", "m_2h", "m_11h",
    [part("reagent1", "m_1", "reagent", (4,)), part("solvent", "m_etoac", "solvent", (4,)),
     part("base", "m_tbuamine", "reagent", (4,)), part("hydrogen", "m_h2", "reagent", (4,)),
     part("catalyst", "m_pdc", "catalyst", (4,))],
    [step("annulation", (4,), "1, EtOAc, rt, 30 min", ("rt", "rt"), ("30 min", {"number": 30, "unit": "min"}), ("substrate", "reagent1", "solvent")),
     step("base", (4,), "t-BuNH2, EtOAc, 0 °C, 5 min", ("0 °C", {"number": 0, "unit": "°C"}), ("5 min", {"number": 5, "unit": "min"}), ("base", "solvent")),
     step("hydrogenation", (4,), "H2 (1 atm), Pd/C, EtOAc, rt, 72 h", ("rt", "rt"), ("72 h", {"number": 72, "unit": "h"}), ("hydrogen", "catalyst", "product"))],
    59, "one-pot annulation/detrifylation/hydrogenation", ["C-C", "C-N", "DB-C-C", "C-H"])

scheme7_record("scheme7_2c_to_12c", "m_2c", "m_12c",
    [part("reagent1", "m_1", "reagent", (4,)), part("solvent", "m_etoac", "solvent", (4,)),
     part("base", "m_tbuamine", "reagent", (4,)), part("hydrogen", "m_h2", "reagent", (4,)),
     part("catalyst", "m_pdc", "catalyst", (4,))],
    [step("annulation", (4,), "1, EtOAc, rt, 30 min", ("rt", "rt"), ("30 min", {"number": 30, "unit": "min"}), ("substrate", "reagent1", "solvent")),
     step("base", (4,), "t-BuNH2, EtOAc, 0 °C, 5 min", ("0 °C", {"number": 0, "unit": "°C"}), ("5 min", {"number": 5, "unit": "min"}), ("base", "solvent")),
     step("hydrogenation", (4,), "H2 (1 atm), Pd/C, EtOAc, rt, 2 h", ("rt", "rt"), ("2 h", {"number": 2, "unit": "h"}), ("hydrogen", "catalyst", "product"))],
    78, "one-pot annulation/detrifylation/hydrogenation", ["C-C", "C-N", "DB-C-C", "C-H"])

scheme7_record("scheme7_2a_to_13a", "m_2a", "m_13a",
    [part("reagent1", "m_1", "reagent", (4,)), part("solvent1", "m_etoac", "solvent", (4,)),
     part("reductant", "m_lialh4", "reagent", (4,)), part("solvent2", "m_ether", "solvent", (4,)),
     part("silica", "m_silica", "purification_agent", (4,), stage="purification")],
    [step("annulation", (4,), "1, EtOAc, rt, 30 min", ("rt", "rt"), ("30 min", {"number": 30, "unit": "min"}), ("substrate", "reagent1", "solvent1")),
     step("reduction", (4,), "LiAlH4, Et2O, 0 °C→rt, 3 h", ("0 °C→rt", "warming from zero degrees Celsius to room temperature"), ("3 h", {"number": 3, "unit": "h"}), ("reductant", "solvent2")),
     step("purification", (4,), "silica gel chromatography", None, None, ("silica", "product"))],
    65, "one-pot annulation followed by reductive ring opening", [], unresolved_cls=True)

scheme7_record("scheme7_2p_to_14a", "m_2p", "m_14a",
    [part("reagent1", "m_1", "reagent", (4,)), part("solvent1", "m_etoac", "solvent", (4,)),
     part("ozone", "m_ozone", "reagent", (4,)), part("solvent2", "m_dcm", "solvent", (4,))],
    [step("annulation", (4,), "1, EtOAc, rt, 30 min", ("rt", "rt"), ("30 min", {"number": 30, "unit": "min"}), ("substrate", "reagent1", "solvent1")),
     step("ozonolysis", (4,), "O3, DCM, -78 °C, 2 h", ("-78 °C", {"number": -78, "unit": "°C"}), ("2 h", {"number": 2, "unit": "h"}), ("ozone", "solvent2", "product"))],
    55, "one-pot annulation and ozonolysis to enal", ["C-C", "C-N", "DB-C-O"])

scheme7_record("scheme7_8b_to_15b", "m_8b", "m_15b",
    [part("reagent1", "m_1", "reagent", (4,)), part("solvent", "m_etoac", "solvent", (4,)),
     part("base", "m_tbuamine", "reagent", (4,))],
    [step("annulation", (4,), "1, EtOAc, 90 °C, sealed tube, 30 min", ("90 °C", {"number": 90, "unit": "°C"}), ("30 min", {"number": 30, "unit": "min"}), ("substrate", "reagent1", "solvent")),
     step("base", (4,), "t-BuNH2, EtOAc, 0 °C, 5 min", ("0 °C", {"number": 0, "unit": "°C"}), ("5 min", {"number": 5, "unit": "min"}), ("base", "solvent", "product"))],
    77, "one-pot allenol annulation and detrifylation", ["C-C", "C-N", "DB-C-C"])

scheme7_record("scheme7_15b_to_16b", "m_15b", "m_16b",
    [part("catalyst1", "m_ipraucl", "catalyst", (4,), mol_percent=5, equiv_ref="m_15b"),
     part("catalyst2", "m_agotf", "catalyst", (4,), mol_percent=5, equiv_ref="m_15b"),
     part("solvent", "m_iproh", "solvent", (4,))],
    [step("hydration", (4,), "5 mol% IPrAuCl, 5 mol% AgOTf, i-PrOH, 80 °C, 24 h",
          ("80 °C", {"number": 80, "unit": "°C"}), ("24 h", {"number": 24, "unit": "h"}),
          ("substrate", "catalyst1", "catalyst2", "solvent", "product"))],
    65, "Au-catalysed hydration of allenyl/trienyl substrate to ketone", ["DB-C-O", "C-H"])

# Proposed mechanism is not a performed experiment.
reaction("scheme8_proposed_mechanism", (3, 4), "proposed_mechanism", ["other"],
    {"page": 4, "scheme": 8}, "Rationalization for generation of bis(triflyl)ated tetrahydroquinolines 3 and 9",
    "proposed aza-ortho-quinone-methide/zwitterion pathway",
    [part("precursor", "m_g_2", "reactant", (4,)), part("acid", "m_A", "reagent", (4,)),
     part("int_i", "m_int_i", "unknown", (4,)), part("int_ii", "m_int_ii", "unknown", (4,)),
     part("product", "m_g_3", "product", (4,), "proposed product class")],
    [step("proposal", (4,), "dehydration to INT-I; nucleophilic attack on in situ Tf2C=CH2 to INT-II; ring closure to 3/9", None, None,
          ("precursor", "acid", "int_i", "int_ii", "product"))], [],
    classification(["C-C", "C-N"], (4,), extra_limitations=("Arrow-pushing proposal is not experimental proof of elementary steps.",)),
    notes=("Authors use 'should' for the proposed intermediates and pathway; this record is not a specific experiment.",),
    unresolved=("Intermediate structures and atom identities are proposed, not experimentally established in the supplied main text.",))

bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": f"{DOC}:page1:prior_work", "document_id": DOC, "kind": "literature_examples",
     "evidence_ids": [ev(1)], "content": "Scheme 1a-c depicts cited Povarov reactions, partial quinoline reduction and metal-catalysed cyclizations; these are prior-work examples, not experiments of this paper."},
    {"context_id": f"{DOC}:page1:title", "document_id": DOC, "kind": "bibliographic_metadata",
     "evidence_ids": [ev(1)], "title": "A Convenient Formal [4 + 2] Heterocylization Route to Bis(triflyl)tetrahydroquinolines",
     "doi": "10.1002/chem.202102254"},
    {"context_id": f"{DOC}:page5:si_note", "document_id": DOC, "kind": "availability_note",
     "evidence_ids": [ev(5)], "content": "Reference [12] points readers to Supporting Information for TLC pictures, but no SI file was explicitly supplied or associated in this session."},
    {"context_id": f"{DOC}:page5:references", "document_id": DOC, "kind": "non_experimental_page",
     "evidence_ids": [ev(5)], "content": "Page 5 contains references [10]-[13] and manuscript dates; no new reaction arrow or experimental entry is present."},
]

st.write(RUN / "submitted.json", bundle)

inventory = {
    "document_id": DOC,
    "pages_inspected": [1, 2, 3, 4, 5],
    "machine_extraction": True,
    "semantic_accuracy_verified": False,
    "reaction_count": len(reactions),
    "molecule_count": len(molecules),
    "breakdown": {
        "generic_scheme": sum(r["record_type"] == "generic_scheme" for r in reactions),
        "specific_experiment": sum(r["record_type"] == "specific_experiment" for r in reactions),
        "proposed_mechanism": sum(r["record_type"] == "proposed_mechanism" for r in reactions),
    },
    "notes": [
        "Every page image was viewed; scheme crops were generated with source-lineage receipts.",
        "Paper-specific structures remain unresolved rather than guessing complete/isomeric SMILES.",
        "Only explicitly named simple solvents/reagents receive candidate SMILES and standard-stereo diagnostics.",
        "No supporting-information association was inferred.",
    ],
}
st.write(RUN / "raw_model/extraction_inventory.json", inventory)
st.write(RUN / "raw_model/stereo_checks.json", {
    "tool": "standard_records.stereo_inspect (same implementation as reaction_tools.py standard stereo)",
    "source_verified": False,
    "checks": stereo_checks,
})

by_page = {
    page: [r["reaction_id"] for r in reactions if ev(page) in r["evidence_ids"]]
    for page in range(1, 6)
}
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 1,
            "state": "attempted",
            "regions": [
                {"region": "title, abstract and bibliographic header", "finding": "Article title, authors and DOI; no entry-specific experimental procedure.", "evidence_ids": [ev(1)]},
                {"region": "Scheme 1a-c", "finding": "Cited prior-work Povarov, quinoline reduction and metal-catalysed cyclization examples; retained as source_context, not current-paper experiments.", "evidence_ids": [ev(1)]},
                {"region": "Scheme 1d", "finding": "Generic incidental-discovery formal [4+2] heterocyclization in EtOAc at rt without catalyst; one generic-scheme record.", "evidence_ids": [ev(1)]},
            ],
            "record_ids": by_page[1],
            "unresolved": ["Scheme 1d is a generic overview; compound-specific yield, time, charges and workup are not supplied on this page."],
        },
        {
            "page": 2,
            "state": "attempted",
            "regions": [
                {"region": "left-column optimization discussion and TEMPO sentence", "finding": "Seven Scheme 2 optimization rows plus one qualitative TEMPO radical-trap control.", "evidence_ids": [ev(2)]},
                {"region": "Scheme 2", "finding": "One condition per row: hot water no reaction; EtOH, MeCN, toluene, DCE, EtOAc and acetone entries with displayed equivalents, temperature, time and yield.", "evidence_ids": [ev(2)]},
                {"region": "Scheme 3", "finding": "Generic scope and 17 separate 2a-q to 3a-q entries; isolated yields/times, 90 °C sealed-tube footnote for 3n/3o, and 58:42 Z:E footnote for 3p.", "evidence_ids": [ev(2)]},
                {"region": "Scheme 3 lower failed precursor and page-bottom prose", "finding": "3p-Me gives a complicated mixture; prose also introduces Scheme 4/5 experiments continued visually on page 3.", "evidence_ids": [ev(2)]},
            ],
            "record_ids": by_page[2],
            "unresolved": ["Scheme 2 yield measurement basis is not stated.", "TEMPO amount and numerical result are not stated.", "SI-level charges, workup, purification and characterization were not available."],
        },
        {
            "page": 3,
            "state": "attempted",
            "regions": [
                {"region": "Scheme 4a-b", "finding": "Four allylic-alcohol scope entries and the distinct tertiary-alcohol 5 to exclusive product 6 transformation.", "evidence_ids": [ev(3)]},
                {"region": "Scheme 5", "finding": "Four separate one-pot heterocyclization/base-treatment entries 2a,h,j,n to 7a,h,j,n; overall entry times cannot be allocated between stages.", "evidence_ids": [ev(3)]},
                {"region": "lower-left prose", "finding": "Alkynol-substituted anilides gave complex mixtures; retained as a failed/negative scope experiment without invented product.", "evidence_ids": [ev(3)]},
                {"region": "Scheme 6", "finding": "Eight allenol scope entries 8a-h to 9a-h and one deuterated [D2]-1/[D2]-9b experiment, each with displayed isolated yield and time.", "evidence_ids": [ev(3)]},
                {"region": "right-column product-elaboration and mechanism lead-in", "finding": "Prose identifies Scheme 7 one-pot elaborations and introduces the proposed Scheme 8 pathway; records merge this text with page-4 figures.", "evidence_ids": [ev(3)]},
            ],
            "record_ids": by_page[3],
            "unresolved": ["Allene axial stereochemistry is not assigned in the source.", "Alkynol control lacks concrete substrate labels and conditions.", "Scheme 5 stage-specific time allocation is not reported."],
        },
        {
            "page": 4,
            "state": "attempted",
            "regions": [
                {"region": "Scheme 7a-d", "finding": "Eight product-elaboration arrows/sequences recorded separately, including Suzuki coupling, hydrogenations, reductive ring opening, ozonolysis, detrifylation and Au-catalysed hydration.", "evidence_ids": [ev(4)]},
                {"region": "Scheme 8", "finding": "Author-proposed aza-ortho-quinone-methide/zwitterion pathway; recorded as proposed_mechanism, not experimental proof.", "evidence_ids": [ev(4)]},
                {"region": "summary and acknowledgements", "finding": "Narrative summary states EtOAc/ambient, metal/additive/irradiation-free main protocol; no additional entry-specific experiment.", "evidence_ids": [ev(4)]},
            ],
            "record_ids": by_page[4],
            "unresolved": ["Scheme 7 yields do not state measurement basis.", "Several one-pot overall sequences omit intermediate material balances and complete charges.", "Scheme 8 intermediates are proposed, not observed."],
        },
        {
            "page": 5,
            "state": "attempted",
            "regions": [
                {"region": "references [10]-[13] and manuscript dates", "finding": "No new reaction arrow, table row or experimental procedure; reference [12] only points to SI TLC pictures.", "evidence_ids": [ev(5)]},
            ],
            "record_ids": by_page[5],
            "unresolved": ["Supporting Information is mentioned by the paper but was not explicitly supplied or associated in this session."],
        },
    ],
}
st.write(RUN / "coverage/pages_0001_0005.json", coverage)
print(json.dumps(inventory, ensure_ascii=False, indent=2))
