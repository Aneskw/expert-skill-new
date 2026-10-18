#!/usr/bin/env python3
"""Build the source-specific v3.1 bundle for Chem. Commun. 2021, 57, 5394–5397.

The script consumes only the clean standard-init template and facts transcribed from
PDF pages 1–4. It deliberately does not read previous experiment records.
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
DOC = "doc_a996ed1ac1a2b8f0"
EV = {p: [f"{DOC}:page{p}:text", f"{DOC}:page{p}:image"] for p in range(1, 5)}


def write(path, obj):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n")


def src(raw, normalized=None, ev=(), method="original_extraction"):
    if normalized is None:
        normalized = raw
    return st.fact(raw, normalized, ev, method=method)


def missing(ev=()):
    return st.fact(evidence=ev)


def inferred(normalized, ev, derivation):
    return st.fact(None, normalized, ev, method="inference", derivation=derivation)


def derived(normalized, ev, derivation):
    return st.fact(None, normalized, ev, method="structure_derivation", derivation=derivation)


def conflict(alternatives, ev):
    f = st.fact(alternatives[0][0], None, ev)
    f["data_status"] = "conflict"
    f["alternatives"] = [src(raw, norm, aev) for raw, norm, aev in alternatives]
    return f


bundle = json.loads((RUN / "raw_model/standard_init.json").read_text())
bundle["documents"][0]["title"] = "Engaging yne-allenones in tunable catalytic silane-mediated conjugate transfer reductions"
bundle["documents"][0]["doi"] = "10.1039/d1cc01422f"
molecules = {}
stereo_checks = {}


def set_status(m, key, value, ev, reason):
    m["stereochemistry"][key] = inferred(value, ev, reason)


def molecule(mid, label, name, ev, smiles=None, *, chirality="unspecified", geometry="unspecified", rep=None,
             note=None, assignment=None, axial=None, mixture=None):
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = src(label, label, ev)
    m["names"]["original"] = src(name, name, ev)
    m["names"]["english"] = src(name, name, ev)
    m["cas_numbers"] = missing(ev)
    m["identity_note"] = src(note, note, ev) if note else missing(ev)
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = derived(smiles, ev, "Candidate graph transcribed from the source drawing/name; parsing does not verify identity.")
        diag = st.stereo_inspect(smiles)
        stereo_checks[mid] = diag
        m["stereochemistry"]["isomeric_smiles"] = derived(diag["canonical_isomeric_smiles"], ev, "RDKit canonicalization of the source-derived candidate graph.")
        m["stereochemistry"]["identity_key"] = derived(diag["stereo_graph_key"], ev, "Stereo graph key calculated from the candidate isomeric SMILES.")
        rep = rep or "partial"
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing(ev)
        rep = rep or "not_encoded"
    set_status(m, "chirality_status", chirality, ev, "Status preserves source assignment or explicit lack of assignment; it is not inferred from missing wedges alone.")
    set_status(m, "geometry_status", geometry, ev, "Status preserves source geometry information and unresolved elements.")
    set_status(m, "representation_status", rep, ev, "Representation status reflects whether the source stereochemistry is encoded in the candidate graph.")
    if assignment is not None:
        m["stereochemistry"]["assignment"] = src(assignment, assignment, ev)
    if axial is not None:
        m["stereochemistry"]["axial"] = inferred(axial, ev, "Local allene-axis assessment from the source-derived candidate graph; RDKit axial support is limited.")
    if mixture is not None:
        m["stereochemistry"]["mixture"] = src(mixture, mixture, ev)
    molecules[mid] = m
    return mid


def generic_molecule(mid, label, name, fixed, variables, ev):
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = src(label, label, ev)
    m["names"]["original"] = src(name, name, ev)
    m["names"]["english"] = src(name, name, ev)
    m["cas_numbers"] = missing(ev)
    m["representation_type"] = "generic_template"
    m["smiles"] = missing(ev)
    m["generic_template"] = {
        "fixed_structure": src(fixed, fixed, ev),
        "variables": [{
            "label": v,
            "attachment_points": src(f"source-drawn attachment for {v}", f"source-drawn attachment for {v}", ev),
            "definition": src(f"{v} as printed in Scheme 1", f"{v} as printed in Scheme 1", ev),
            "allowed_values": missing(ev),
            "constraints": missing(ev),
        } for v in variables],
        "instance_mapping": [],
        "reason": src("R-group family in Scheme 1", "R-group family in Scheme 1", ev),
    }
    for k, v in (("chirality_status", "unresolved"), ("geometry_status", "unresolved"), ("representation_status", "not_encoded")):
        set_status(m, k, v, ev, "Generic R-group drawing does not identify one fully specified stereochemical molecule.")
    m["identity_note"] = src("Generic family; not a concrete wildcard SMILES.", "Generic family; not a concrete wildcard SMILES.", ev)
    molecules[mid] = m
    return mid


# Generic source families.
generic_molecule(f"{DOC}:generic_yne_allenone_1", "1", "generic yne-allenone 1", "benzene-tethered yne-allenone with R1/R2 variables", ["R1", "R2"], EV[2])
generic_molecule(f"{DOC}:generic_silane_2", "2", "generic hydrosilane 2", "[Si]-H hydrosilane family", ["[Si]"], EV[2])
generic_molecule(f"{DOC}:generic_product_3", "3", "generic dihydrocyclobuta[a]naphthalen-4-ol 3", "Fe-path cyclobutarene alcohol with R1/R2 variables", ["R1", "R2"], EV[2])
generic_molecule(f"{DOC}:generic_product_4", "4", "generic cyclobuta[a]naphthalen-4(2H)-one 4", "Cu-path cyclobutarenone with R1/R2 variables", ["R1", "R2"], EV[3])
generic_molecule(f"{DOC}:prior_key_intermediate", "key intermediate", "prior-work cyclobutene intermediate family", "intramolecular [2+2] cycloadduct with R1", ["R1"], EV[2])


# Simple reagents/solvents. Candidate SMILES are name-derived and remain unverified.
simple = [
    ("phsih3", "2a", "phenylsilane", "[SiH3]c1ccccc1", EV[1] + EV[2]),
    ("feacac3", "Fe(acac)3", "iron(III) acetylacetonate", None, EV[1] + EV[2]),
    ("ruacac3", "Ru(acac)3", "ruthenium(III) acetylacetonate", None, EV[2]),
    ("meoh", "CH3OH", "methanol", "CO", EV[2]),
    ("thf", "THF", "tetrahydrofuran", "C1CCOC1", EV[2]),
    ("etoh", "EtOH", "ethanol", "CCO", EV[1] + EV[2] + EV[3] + EV[4]),
    ("et3sih", "Et3SiH", "triethylsilane", "CC[SiH](CC)CC", EV[2] + EV[3]),
    ("ph3sih", "Ph3SiH", "triphenylsilane", "[SiH](c1ccccc1)(c1ccccc1)c1ccccc1", EV[2] + EV[3]),
    ("pmhs", "PMHS", "polymethylhydrosiloxane", None, EV[2]),
    ("medimethoxysilane", "2b", "methyldimethoxysilane", "C[SiH](OC)OC", EV[3] + EV[4]),
    ("cucn", "CuCN", "copper(I) cyanide", None, EV[3] + EV[4]),
    ("dppf", "DPPF", "1,1'-bis(diphenylphosphino)ferrocene", None, EV[3] + EV[4]),
    ("naome", "NaOMe", "sodium methoxide", "C[O-].[Na+]", EV[3] + EV[4]),
    ("dce", "DCE", "1,2-dichloroethane", "ClCCCl", EV[3]),
    ("cui", "CuI", "copper(I) iodide", None, EV[3]),
    ("cubr", "CuBr", "copper(I) bromide", None, EV[3]),
    ("cuoac2", "Cu(OAc)2", "copper(II) acetate", None, EV[3]),
    ("tempo", "TEMPO", "2,2,6,6-tetramethyl-1-piperidinyloxy", "CC1(C)CCCC(C)(C)N1[O]", EV[3]),
]
for key, label, name, smi, ev in simple:
    molecule(f"{DOC}:{key}", label, name, ev, smi,
             chirality="unspecified" if key in {"phsih3", "et3sih", "ph3sih", "medimethoxysilane"} else "achiral" if smi else "unresolved",
             geometry="not_applicable" if smi else "unresolved", rep="partial" if smi else "not_encoded",
             note="Simple reagent/solvent candidate; no CAS number was supplied.")

# The cobalt salt identity conflicts between Table 1 and prose.
cobalt = molecule(f"{DOC}:coacac_conflict", "Co(acac)n", "cobalt acetylacetonate (source-conflicted oxidation state)", EV[2], None,
                  chirality="unresolved", geometry="unresolved", rep="not_encoded")
molecules[cobalt]["names"]["original"] = conflict([
    ("Co(acac)3 (Table 1 entry 6)", "Co(acac)3", EV[2]),
    ("Co(acac)2 (adjacent prose)", "Co(acac)2", EV[2]),
], EV[2])


# Concrete yne-allenones. Internal-ring variants are conservatively unencoded when position is compact/ambiguous.
base_smi = "C=C=CC(=O)c1ccccc1C#C{}"
aryl = {
    "ph": "c1ccccc1", "ptol": "c1ccc(C)cc1", "mtol": "c1cccc(C)c1", "pcl": "c1ccc(Cl)cc1",
    "ocl": "c1ccccc1Cl", "mcl": "c1cccc(Cl)c1", "pbr": "c1ccc(Br)cc1", "pcf3": "c1ccc(C(F)(F)F)cc1",
    "pco2me": "c1ccc(C(=O)OC)cc1", "thienyl": "c1ccsc1", "pet": "c1ccc(CC)cc1",
}
fe_specs = [
    ("1a", "ph"), ("1b", "ptol"), ("1c", "mtol"), ("1d", "pcl"), ("1e", "ocl"), ("1f", "pbr"),
    ("1g", "pcf3"), ("1h", "pco2me"), ("1i", "thienyl"),
]
substrates = {}
for label, code in fe_specs:
    mid = f"{DOC}:fe_{label}"
    molecule(mid, label, f"yne-allenone {label}", EV[2], base_smi.format(aryl[code]), chirality="achiral",
             geometry="unspecified", rep="partial", axial=[{"axis": "terminal-CH2 allene", "status": "not stereogenic on candidate local graph"}],
             note="Source-derived candidate graph; terminal CH2 allene has identical substituents at one end.")
    substrates[("fe", label)] = mid
for label in ["1j", "1k", "1l", "1m", "1n", "1o", "1p"]:
    mid = f"{DOC}:fe_{label}"
    molecule(mid, label, f"internally substituted yne-allenone {label}", EV[2], None, chirality="unspecified",
             geometry="unspecified", rep="not_encoded", axial=[{"axis": "allene", "status": "unassigned; terminal CH2 is visible"}],
             note="Compact internal-arene substitution drawing was not forced into a complete SMILES.")
    substrates[("fe", label)] = mid

# Cu-scope label series is separately namespaced because the paper reuses 1b–1m for a different enumeration.
cu_defs = [("1a", substrates[("fe", "1a")]), ("1b", substrates[("fe", "1b")]), ("1c", substrates[("fe", "1c")])]
for label, code in [("1d", "pet"), ("1e", "pcl"), ("1f", "mcl"), ("1g", "ocl"), ("1h", "thienyl")]:
    mid = f"{DOC}:cu_{label}"
    molecule(mid, label, f"Cu-scope yne-allenone {label}", EV[3], base_smi.format(aryl[code]), chirality="achiral",
             geometry="unspecified", rep="partial", axial=[{"axis": "terminal-CH2 allene", "status": "not stereogenic on candidate local graph"}])
    cu_defs.append((label, mid))
for label in ["1i", "1j", "1k", "1l", "1m"]:
    mid = f"{DOC}:cu_{label}"
    molecule(mid, label, f"Cu-scope internally substituted yne-allenone {label}", EV[3], None, chirality="unspecified",
             geometry="unspecified", rep="not_encoded", axial=[{"axis": "allene", "status": "unassigned; terminal CH2 is visible"}])
    cu_defs.append((label, mid))
cu_substrates = dict(cu_defs)

# Scope products. No source stereochemical assignment is invented.
products3, products4 = {}, {}
for c in "abcdefghijklmnop":
    label = f"3{c}"
    products3[label] = molecule(f"{DOC}:{label}", label, f"dihydrocyclobuta[a]naphthalen-4-ol {label}", EV[1] + EV[2] if label == "3a" else EV[2],
                                None, chirality="unspecified", geometry="unspecified", rep="not_encoded",
                                note="No R/S, dr, ee, or racemic assignment is reported on pages 1–4.")
for c in "abcdefghijklm":
    label = f"4{c}"
    products4[label] = molecule(f"{DOC}:{label}", label, f"cyclobuta[a]naphthalen-4(2H)-one {label}", EV[3],
                                None, chirality="unspecified", geometry="unspecified", rep="not_encoded",
                                note="No R/S, dr, ee, or racemic assignment is reported on pages 1–4.")


reactions = []


def participant(pid, mid, role, ev, *, amount=None, mass=None, volume=None, equiv=None, molpct=None, ref=None, stage="reaction", method="original_extraction"):
    p = st.new_participant(pid, mid, role, ev, stage=stage)
    p["product_role"] = src("reported product", "reported product", ev) if role == "product" else missing(ev)
    vals = {"amount": amount, "mass": mass, "volume": volume, "equivalents": equiv, "mol_percent": molpct}
    for key, val in vals.items():
        if val is not None:
            raw = val.get("raw_text") or f"{val.get('number')} {val.get('unit')}"
            p["charges"][key] = src(raw, val, ev, method=method)
    if equiv is not None or molpct is not None:
        p["equivalence_reference"] = src(ref, ref, ev, method=method)
    return p


def set_reaction_common(r, *, categories, location, name, transform, ev, scheme_ids=()):
    r["record_type"] = "specific_experiment"
    r["experiment_categories"] = src(categories, categories, ev)
    r["location"] = src(location, location, ev)
    r["general_scheme_ids"] = list(scheme_ids)
    r["reaction_name"] = src(name, name, ev)
    r["transformation_type"] = src(transform, transform, ev)
    r["review"]["unresolved"] = ["SI was advertised but not supplied; full experimental details remain unavailable."]


def set_solvent(r, name, volume, ev, method="original_extraction"):
    r["solvent_conditions"]["composition"] = src(name, name, ev, method=method)
    if volume is not None:
        r["solvent_conditions"]["total_volume"] = src(f"{volume['number']} {volume['unit']}", volume, ev, method=method)


def add_step(r, sid, pids, ev, temp=None, time=None, atmosphere=None, operation="reaction mixture held under the stated conditions", method="original_extraction"):
    s = st.new_step(sid, ev)
    s["participant_ids"] = pids
    s["operation"] = src(operation, operation, ev, method=method)
    if temp is not None:
        s["temperature"] = src(temp if isinstance(temp, str) else f"{temp['number']} {temp['unit']}", temp, ev, method=method)
    if time is not None:
        s["time"] = src(time if isinstance(time, str) else f"{time['number']} {time['unit']}", time, ev, method=method)
    if atmosphere is not None:
        s["atmosphere"] = src(atmosphere, atmosphere, ev, method=method)
    r["steps"].append(s)


def outcome(r, oid, kind, targets, value, ev, *, basis=None, yield_type=None, scope=None):
    o = st.new_outcome(oid, kind, targets, ev)
    raw = value if isinstance(value, str) else value.get("raw_text") or f"{value.get('number')} {value.get('unit')}"
    o["value"] = src(raw, value, ev)
    if basis is not None:
        o["measurement_basis"] = src(basis, basis, ev)
    if yield_type is not None:
        o["yield_type"] = src(yield_type, yield_type, ev)
    if scope is not None:
        o["scope"] = src(scope, scope, ev)
    r["outcomes"].append(o)


def class_partial(r, mode, ev):
    if mode == "fe":
        specs = [("C", "C", "formed"), ("C", "C", "formed"), ("C", "H", "hydrogen_gain"), ("O", "H", "hydrogen_gain")]
    else:
        specs = [("C", "C", "formed"), ("C", "C", "formed"), ("C", "H", "hydrogen_gain"), ("C", "H", "hydrogen_gain")]
    changes, labels = [], {}
    for i, (a, b, typ) in enumerate(specs, 1):
        change = {"change_id": f"b{i}", "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
                  "before_order": 0, "after_order": 1, "change_type": typ}
        changes.append(change)
        lab = st.bond_label(a, b, 1)
        labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(f"b{i}")
    r["classification"]["bond_changes"] = [derived(c, ev, "Class-level net change inferred from the source-drawn [2+2]/conjugate-reduction transformation; exact atom identities are unresolved.") for c in changes]
    r["classification"]["labels"] = derived(list(labels.values()), ev, "Taxonomy labels cover all inferred formed bonds in the partial class-level change set.")
    r["classification"]["coverage"] = inferred("partial", ev, "No safe complete atom mapping was asserted.")
    r["classification"]["limitations"] = ["Exact atom identities and full atom mapping are unresolved.", "Hydrogen donors/acceptors and coproducts are not fully mapped."]


def class_na(r, ev):
    r["classification"]["labels"] = inferred([], ev, "No observed desired transformation; intended bond formation is not classified.")
    r["classification"]["coverage"] = inferred("not_applicable", ev, "No observed desired transformation.")
    r["classification"]["limitations"] = ["Negative/no-reaction result; expected-product bonds were not assigned."]


def make_generic_records():
    lit_id = f"{DOC}:scheme1_prior_work"
    r = st.new_reaction(lit_id, DOC, EV[2]); r["record_type"] = "literature_example"
    r["experiment_categories"] = inferred([], EV[2], "Cited prior-work reaction family, not a current experiment.")
    r["location"] = src("Scheme 1a, previous works", {"page": 2, "region": "Scheme 1a"}, EV[2])
    r["reaction_name"] = src("previous [2+2] cycloaddition chemistry of yne-allenones", "previous [2+2] cycloaddition chemistry of yne-allenones", EV[2])
    r["transformation_type"] = src("cited yne-allenone cycloaddition/reaction family", "cited yne-allenone cycloaddition/reaction family", EV[2])
    r["participants"] = [participant("prior_1", f"{DOC}:generic_yne_allenone_1", "reactant", EV[2]), participant("prior_A", f"{DOC}:prior_key_intermediate", "product", EV[2])]
    r["classification"]["coverage"] = inferred("unresolved", EV[2], "Grouped cited family is not atom-mapped.")
    r["review"]["unresolved"] = ["Grouped cited work is retained as literature context, not a current-paper experiment."]
    reactions.append(r)

    gid = f"{DOC}:scheme1_current_generic"
    r = st.new_reaction(gid, DOC, EV[1] + EV[2]); r["record_type"] = "generic_scheme"
    r["experiment_categories"] = inferred([], EV[1] + EV[2], "Generic R-group scheme, not a specific experiment.")
    r["location"] = src("Scheme 1b-c and page-1 prose", {"pages": [1, 2], "region": "Scheme 1b-c"}, EV[1] + EV[2])
    r["reaction_name"] = src("tunable catalytic [2+2] cycloaddition/silane-mediated conjugate transfer reductions", "tunable catalytic [2+2] cycloaddition/silane-mediated conjugate transfer reductions", EV[1] + EV[2])
    r["transformation_type"] = src("Fe-mediated 1,6- or Cu-mediated 1,4-conjugate transfer reduction", "chemodivergent cycloaddition/conjugate reduction", EV[1] + EV[2])
    r["participants"] = [participant("gen_1", f"{DOC}:generic_yne_allenone_1", "reactant", EV[2]), participant("gen_2", f"{DOC}:generic_silane_2", "reagent", EV[2]), participant("gen_3", f"{DOC}:generic_product_3", "product", EV[2]), participant("gen_4", f"{DOC}:generic_product_4", "product", EV[2])]
    r["classification"]["coverage"] = inferred("unresolved", EV[1] + EV[2], "Generic family contains two product pathways and no atom mapping.")
    r["review"]["unresolved"] = ["R-group family; instances are represented by specific table/scope records."]
    reactions.append(r)
    return gid


GEN = make_generic_records()


def fe_participants(substrate, product, ev, *, silane=f"{DOC}:phsih3", metal=f"{DOC}:feacac3", solvent=f"{DOC}:etoh",
                    sub_amt=0.2, sil_amt=0.2, metal_amt=0.4, solvent_vol=2.5, inherited=False):
    method = "inherited_general_condition" if inherited else "original_extraction"
    parts = [participant("substrate", substrate, "reactant", ev, amount={"number": sub_amt, "unit": "mmol"}, ref=None, method=method),
             participant("silane", silane, "reagent", ev, amount={"number": sil_amt, "unit": "mmol"}, method=method),
             participant("metal", metal, "catalyst", ev, amount={"number": metal_amt, "unit": "mmol"}, method=method),
             participant("solvent", solvent, "solvent", ev, volume={"number": solvent_vol, "unit": "mL"}, method=method),
             participant("product", product, "product", ev)]
    return parts


table1 = [
    (1, "None", 81), (2, "1.5 equiv. of Fe(acac)3", 64), (3, "1.0 equiv. of Fe(acac)3", 38),
    (4, "2.5 equiv. of Fe(acac)3", 67), (5, "Ru(acac)3 instead of Fe(acac)3", 62),
    (6, "Co(acac)3 instead of Fe(acac)3", 59), (7, "CH3OH instead of EtOH", 58),
    (8, "THF instead of EtOH", 42), (9, "At 50 °C", 72), (10, "At 80 °C", 60),
    (11, "Et3SiH instead of PhSiH3", 45), (12, "Ph3SiH instead of PhSiH3", 50),
    (13, "PMHS instead of PhSiH3", 26),
]
for entry, variation, y in table1:
    ev = EV[1] + EV[2] if entry == 1 else EV[2]
    r = st.new_reaction(f"{DOC}:table1_entry{entry}", DOC, ev)
    set_reaction_common(r, categories=["optimization"] + (["substrate_scope"] if entry == 1 else []),
                        location=f"Table 1 entry {entry}" + ("; Scheme 2 product 3a; page-1 prose" if entry == 1 else ""),
                        name="Fe-mediated [2+2] cycloaddition/1,6-conjugate reduction", transform="cycloaddition and 1,6-conjugate transfer reduction", ev=ev, scheme_ids=[GEN])
    sub, prod = substrates[("fe", "1a")], products3["3a"]
    sil, metal, solv = f"{DOC}:phsih3", f"{DOC}:feacac3", f"{DOC}:etoh"
    sil_amt, metal_amt, solv_vol, temp = 0.2, 0.4, 2.5, 60
    if entry == 2: metal_amt = 0.3
    if entry == 3: metal_amt = 0.2
    if entry == 4: metal_amt = 0.5
    if entry == 5: metal = f"{DOC}:ruacac3"
    if entry == 6: metal = f"{DOC}:coacac_conflict"
    if entry == 7: solv = f"{DOC}:meoh"
    if entry == 8: solv = f"{DOC}:thf"
    if entry == 9: temp = 50
    if entry == 10: temp = 80
    if entry == 11: sil = f"{DOC}:et3sih"
    if entry == 12: sil = f"{DOC}:ph3sih"
    if entry == 13: sil = f"{DOC}:pmhs"
    r["participants"] = fe_participants(sub, prod, ev, silane=sil, metal=metal, solvent=solv, sil_amt=sil_amt, metal_amt=metal_amt, solvent_vol=solv_vol)
    r["reaction_scale"] = src("1a (0.2 mmol)", {"reference_molecule_id": sub, "amount": {"number": 0.2, "unit": "mmol"}}, EV[2])
    set_solvent(r, molecules[solv]["names"]["original"]["normalized_value"] or "source-conflicted solvent", {"number": solv_vol, "unit": "mL"}, EV[2])
    add_step(r, "reaction", [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"], ev,
             temp={"number": temp, "unit": "°C"}, time={"number": 1.5, "unit": "h"}, atmosphere="Ar")
    outcome(r, "yield_3a", "yield", [prod], {"number": y, "unit": "%"}, ev,
            basis="isolated yield based on 1a", yield_type="isolated yield")
    r["notes"].append(src(f"Variation from standard conditions: {variation}", variation, EV[2]))
    if entry == 1:
        r["notes"].append(inferred("Table 1 entry 1 and Scheme 2 product 3a are the same 0.2 mmol/81% standard experiment and were merged.", ev, "Conditions, substrate/product labels, and yield match across the two regions."))
    if entry == 6:
        r["notes"].append(src("Table 1 prints Co(acac)3, while adjacent prose says Co(acac)2.", "source conflict retained", EV[2]))
        r["review"]["unresolved"].append("Cobalt acetylacetonate identity conflicts: Table 1 Co(acac)3 versus prose Co(acac)2.")
    class_partial(r, "fe", ev)
    reactions.append(r)


fe_scope = [
    ("3b", "1b", 71), ("3c", "1c", 59), ("3d", "1d", 58), ("3e", "1e", 68), ("3f", "1f", 50),
    ("3g", "1g", 55), ("3h", "1h", 60), ("3i", "1i", 88), ("3j", "1j", 81), ("3k", "1k", 75),
    ("3l", "1l", 45), ("3m", "1m", 47), ("3n", "1n", 47), ("3o", "1o", 82), ("3p", "1p", 48),
]
for plabel, slabel, y in fe_scope:
    ev = EV[2]; sub, prod = substrates[("fe", slabel)], products3[plabel]
    r = st.new_reaction(f"{DOC}:scheme2_{plabel}", DOC, ev)
    set_reaction_common(r, categories=["substrate_scope"], location=f"Scheme 2 product {plabel}",
                        name="Fe-mediated [2+2] cycloaddition/1,6-conjugate reduction", transform="cycloaddition and 1,6-conjugate transfer reduction", ev=ev, scheme_ids=[GEN])
    r["participants"] = fe_participants(sub, prod, ev, inherited=True)
    r["reaction_scale"] = src("standard conditions referenced; 1 at 0.2 mmol", {"reference_molecule_id": sub, "amount": {"number": 0.2, "unit": "mmol"}}, ev, method="inherited_general_condition")
    set_solvent(r, "EtOH", {"number": 2.5, "unit": "mL"}, ev, method="inherited_general_condition")
    add_step(r, "reaction", [p["participant_id"] for p in r["participants"]], ev, temp={"number": 60, "unit": "°C"}, time={"number": 1.5, "unit": "h"}, atmosphere="Ar", method="inherited_general_condition")
    outcome(r, f"yield_{plabel}", "yield", [prod], {"number": y, "unit": "%"}, ev, yield_type="reported scheme yield")
    r["review"]["unresolved"].append("Scheme 2 does not state the yield measurement basis; Table 1's isolated-yield footnote was not propagated.")
    class_partial(r, "fe", ev); reactions.append(r)

# Scale-up.
sub, prod = substrates[("fe", "1a")], products3["3a"]
r = st.new_reaction(f"{DOC}:scheme2_scale_up_3a", DOC, EV[2])
set_reaction_common(r, categories=["scale_up"], location="Scheme 2b amplification reaction", name="Fe-mediated scale-up reduction",
                    transform="scale-up [2+2]/1,6-conjugate reduction", ev=EV[2], scheme_ids=[GEN])
r["participants"] = [
    participant("substrate", sub, "reactant", EV[2], amount={"number": 4.0, "unit": "mmol"}),
    participant("silane", f"{DOC}:phsih3", "reagent", EV[2]),
    participant("metal", f"{DOC}:feacac3", "catalyst", EV[2]),
    participant("solvent", f"{DOC}:etoh", "solvent", EV[2]),
    participant("product", prod, "product", EV[2]),
]
r["participants"][-1]["charges"]["mass"] = src("0.69 g", {"number": 0.69, "unit": "g"}, EV[2])
r["reaction_scale"] = src("1a (4.0 mmol)", {"reference_molecule_id": sub, "amount": {"number": 4.0, "unit": "mmol"}}, EV[2])
set_solvent(r, "EtOH", None, EV[2], method="inherited_general_condition")
add_step(r, "reaction", [p["participant_id"] for p in r["participants"]], EV[2], temp={"number": 60, "unit": "°C"}, time={"number": 1.5, "unit": "h"}, atmosphere="Ar", method="inherited_general_condition")
outcome(r, "yield_3a", "yield", [prod], {"number": 70, "unit": "%"}, EV[2], yield_type="reported scheme yield")
outcome(r, "amount_3a", "product_amount", [prod], {"number": 0.69, "unit": "g"}, EV[2])
r["review"]["unresolved"] += ["Scale-up yield measurement basis is unreported.", "EtOH volume is not reported for the 4.0 mmol experiment; the 2.5 mL small-scale volume was not copied."]
class_partial(r, "fe", EV[2]); reactions.append(r)


def cu_participants(substrate, product, ev, *, silane=f"{DOC}:medimethoxysilane", copper=f"{DOC}:cucn", inherited=False):
    method = "inherited_general_condition" if inherited else "original_extraction"
    return [
        participant("substrate", substrate, "reactant", ev, amount={"number": 0.2, "unit": "mmol"}, method=method),
        participant("silane", silane, "reagent", ev, equiv={"number": 1.5, "unit": "equiv"}, ref=substrate, method=method),
        participant("copper", copper, "catalyst", ev, molpct={"number": 10, "unit": "mol%"}, ref=substrate, method=method),
        participant("ligand", f"{DOC}:dppf", "ligand", ev, molpct={"number": 10, "unit": "mol%"}, ref=substrate, method=method),
        participant("base", f"{DOC}:naome", "reagent", ev, equiv={"number": 1.5, "unit": "equiv"}, ref=substrate, method=method),
        participant("solvent", f"{DOC}:dce", "solvent", ev, volume={"number": 2.5, "unit": "mL"}, method=method),
        participant("product", product, "product", ev),
    ]


table2 = [
    (1, "None", "yield", {"number": 90, "unit": "%"}),
    (2, "CuI instead of CuCN", "trace", "Trace"),
    (3, "CuBr instead of CuCN", "yield", {"number": 38, "unit": "%"}),
    (4, "Cu(OAc)2 instead of CuCN", "trace", "Trace"),
    (5, "Ph3SiH instead of Me(MeO)2SiH", "yield", {"number": 48, "unit": "%"}),
    (6, "Et3SiH instead of Me(MeO)2SiH", "no_reaction", "NR"),
    (7, "PhSiH3 instead of Me(MeO)2SiH", "no_reaction", "NR"),
]
for entry, variation, kind, val in table2:
    ev = EV[3]; sub, prod = substrates[("fe", "1a")], products4["4a"]
    copper, silane = f"{DOC}:cucn", f"{DOC}:medimethoxysilane"
    if entry == 2: copper = f"{DOC}:cui"
    if entry == 3: copper = f"{DOC}:cubr"
    if entry == 4: copper = f"{DOC}:cuoac2"
    if entry == 5: silane = f"{DOC}:ph3sih"
    if entry == 6: silane = f"{DOC}:et3sih"
    if entry == 7: silane = f"{DOC}:phsih3"
    r = st.new_reaction(f"{DOC}:table2_entry{entry}", DOC, ev)
    set_reaction_common(r, categories=["optimization"] + (["substrate_scope"] if entry == 1 else []),
                        location=f"Table 2 entry {entry}" + ("; Scheme 3 product 4a" if entry == 1 else ""),
                        name="Cu-mediated [2+2] cycloaddition/1,4-conjugate reduction", transform="cycloaddition and 1,4-conjugate transfer reduction", ev=ev, scheme_ids=[GEN])
    r["participants"] = cu_participants(sub, prod, ev, silane=silane, copper=copper)
    r["reaction_scale"] = src("1a (0.2 mmol)", {"reference_molecule_id": sub, "amount": {"number": 0.2, "unit": "mmol"}}, ev)
    set_solvent(r, "DCE", {"number": 2.5, "unit": "mL"}, ev)
    add_step(r, "reaction", [p["participant_id"] for p in r["participants"]], ev, temp={"number": 50, "unit": "°C"}, time={"number": 3, "unit": "h"}, atmosphere="Ar")
    basis = "isolated yield based on 1a" if kind == "yield" else "Table 2 yield column; footnote b states isolated yield based on 1a"
    outcome(r, f"result_4a", kind, [prod] if kind != "no_reaction" else [], val, ev, basis=basis if kind != "no_reaction" else None,
            yield_type="isolated yield" if kind == "yield" else None)
    r["notes"].append(src(f"Variation from standard conditions: {variation}", variation, ev))
    if entry == 1:
        r["notes"].append(inferred("Table 2 entry 1 and Scheme 3 product 4a are the same 0.2 mmol/90% standard experiment and were merged.", ev, "Conditions, labels, and yield match."))
    if kind == "no_reaction": class_na(r, ev)
    else: class_partial(r, "cu", ev)
    reactions.append(r)


cu_scope = [
    ("4b", "1b", 85), ("4c", "1c", 81), ("4d", "1d", 87), ("4e", "1e", 74), ("4f", "1f", 73),
    ("4g", "1g", 70), ("4h", "1h", 83), ("4i", "1i", 85), ("4j", "1j", 78), ("4k", "1k", 73),
    ("4l", "1l", 76), ("4m", "1m", 79),
]
for plabel, slabel, y in cu_scope:
    ev = EV[3]; sub, prod = cu_substrates[slabel], products4[plabel]
    r = st.new_reaction(f"{DOC}:scheme3_{plabel}", DOC, ev)
    set_reaction_common(r, categories=["substrate_scope"], location=f"Scheme 3 product {plabel}",
                        name="Cu-mediated [2+2] cycloaddition/1,4-conjugate reduction", transform="cycloaddition and 1,4-conjugate transfer reduction", ev=ev, scheme_ids=[GEN])
    r["participants"] = cu_participants(sub, prod, ev, inherited=True)
    r["reaction_scale"] = src("standard conditions referenced; 1 at 0.2 mmol", {"reference_molecule_id": sub, "amount": {"number": 0.2, "unit": "mmol"}}, ev, method="inherited_general_condition")
    set_solvent(r, "DCE", {"number": 2.5, "unit": "mL"}, ev, method="inherited_general_condition")
    add_step(r, "reaction", [p["participant_id"] for p in r["participants"]], ev, temp={"number": 50, "unit": "°C"}, time={"number": 3, "unit": "h"}, atmosphere="Ar", method="inherited_general_condition")
    outcome(r, f"yield_{plabel}", "yield", [prod], {"number": y, "unit": "%"}, ev, yield_type="reported scheme yield")
    r["review"]["unresolved"].append("Scheme 3 does not state the yield measurement basis; Table 2's isolated-yield footnote was not propagated.")
    class_partial(r, "cu", ev); reactions.append(r)


# TEMPO controls are separate experiments.
for mode, val, kind in [("fe", "trace", "trace"), ("cu", {"number": 79, "unit": "%"}, "yield")]:
    ev = EV[3]; sub = substrates[("fe", "1a")]; prod = products3["3a"] if mode == "fe" else products4["4a"]
    rid = f"{DOC}:scheme4_tempo_{mode}"
    r = st.new_reaction(rid, DOC, ev)
    set_reaction_common(r, categories=["control", "mechanistic_study"], location=f"Scheme 4a TEMPO control ({mode.upper()} system)",
                        name=f"TEMPO control in the {mode.upper()}-mediated reduction", transform="radical-scavenger control", ev=ev, scheme_ids=[GEN])
    if mode == "fe":
        r["participants"] = fe_participants(sub, prod, ev, inherited=True)
        temp, tm, solv = 60, 1.5, "EtOH"
    else:
        r["participants"] = cu_participants(sub, prod, ev, inherited=True)
        temp, tm, solv = 50, 3, "DCE"
    r["participants"].insert(-1, participant("tempo", f"{DOC}:tempo", "additive", ev, equiv={"number": 2.0, "unit": "equiv"}, ref=sub))
    set_solvent(r, solv, None, ev, method="inherited_general_condition")
    add_step(r, "reaction", [p["participant_id"] for p in r["participants"]], ev, temp={"number": temp, "unit": "°C"}, time={"number": tm, "unit": "h"}, atmosphere="Ar", method="inherited_general_condition")
    outcome(r, f"result_{prod.rsplit(':',1)[-1]}", kind, [prod], val, ev, yield_type="reported scheme result")
    if mode == "fe":
        r["notes"].append(src("the progress of the reaction was completely suppressed", "author interpretation: radical process may be involved", ev))
    else:
        r["notes"].append(src("4a in 79% yield; process did not involve formation of a hydrogen radical", "author interpretation: hydrogen radical not involved", ev))
        r["review"]["unresolved"].append("The 79% yield measurement basis is unreported.")
    class_partial(r, mode, ev); reactions.append(r)


# Proposed mechanisms, not executed intermediate syntheses.
for suffix, page, target, title, transform in [
    ("scheme4b_fe_mechanism", 3, f"{DOC}:generic_product_3", "Plausible reaction pathway for forming 3", "author-proposed Fe-H radical 1,6-reduction cycle"),
    ("scheme5_cu_mechanism", 4, f"{DOC}:generic_product_4", "Plausible mechanism for the formation of 4", "author-proposed Cu-H 1,4-reduction cycle"),
]:
    ev = EV[page]; r = st.new_reaction(f"{DOC}:{suffix}", DOC, ev); r["record_type"] = "proposed_mechanism"
    r["experiment_categories"] = inferred([], ev, "Author-proposed mechanism, not a concrete experiment.")
    r["location"] = src(title, {"page": page, "region": title}, ev)
    r["reaction_name"] = src(title, title, ev); r["transformation_type"] = src(transform, transform, ev)
    r["participants"] = [participant("substrate", f"{DOC}:generic_yne_allenone_1", "reactant", ev), participant("product", target, "product", ev)]
    r["classification"]["coverage"] = inferred("unresolved", ev, "Mechanistic intermediates are not atom-mapped or treated as executed experiments.")
    r["review"]["unresolved"] = ["Plausible mechanism only; intermediates A–F were not isolated experiments."]
    reactions.append(r)


bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": f"{DOC}:metadata_and_esi", "document_id": DOC, "kind": "metadata",
     "evidence_ids": EV[1], "description": "Title, DOI, dates, authors, and an ESI availability notice. No SI file was supplied or associated."},
    {"context_id": f"{DOC}:scope_structure_context", "document_id": DOC, "kind": "structure_drawings",
     "evidence_ids": EV[2] + EV[3], "description": "Scheme 2/3 compact structures support labels and substitution patterns; many complex products remain deliberately unencoded."},
    {"context_id": f"{DOC}:mechanistic_discussion", "document_id": DOC, "kind": "mechanistic_context",
     "evidence_ids": EV[3] + EV[4], "description": "Author discussion of Fe-H radical and Cu-H pathways; executed TEMPO controls are reaction records, proposed cycles are mechanism records."},
    {"context_id": f"{DOC}:negative_spectral_inventory", "document_id": DOC, "kind": "coverage_note",
     "evidence_ids": EV[1] + EV[2] + EV[3] + EV[4], "description": "No standalone mass-spectrum peak table, molecular-weight table, or fragment-relative-abundance block occurs on pages 1–4."},
]

write(RUN / "submitted.json", bundle)
write(RUN / "raw_model/stereo_checks.json", stereo_checks)

# Raw page inventory is deliberately separate from submitted machine records.
page_inventory = {
    "document_id": DOC,
    "pages": [
        {"page": 1, "image_evidence": EV[1][1], "regions": [
            {"region": "title/abstract/ESI notice", "finding": "metadata and unassociated ESI notice"},
            {"region": "right-column current-work prose", "finding": "generic Fe/Cu pathways and standard 3a experiment beginning"},
        ]},
        {"page": 2, "image_evidence": EV[2][1], "regions": [
            {"region": "Scheme 1", "finding": "one prior-work family and one two-branch current generic scheme"},
            {"region": "Table 1", "finding": "13 independent rows; row 1 duplicates Scheme 2 3a; row 6 cobalt identity conflict"},
            {"region": "Scheme 2a", "finding": "16 products 3a–3p; 3a merged with Table 1 row 1"},
            {"region": "Scheme 2b", "finding": "one 4.0 mmol scale-up experiment"},
        ]},
        {"page": 3, "image_evidence": EV[3][1], "regions": [
            {"region": "Table 2", "finding": "7 independent rows; row 1 duplicates Scheme 3 4a"},
            {"region": "Scheme 3", "finding": "13 products 4a–4m; 4a merged with Table 2 row 1"},
            {"region": "Scheme 4a", "finding": "two independent TEMPO controls"},
            {"region": "Scheme 4b", "finding": "one proposed Fe mechanism"},
        ]},
        {"page": 4, "image_evidence": EV[4][1], "regions": [
            {"region": "Scheme 5/top-left prose", "finding": "one proposed Cu mechanism and explanatory prose"},
            {"region": "remainder", "finding": "conclusion, acknowledgements, conflicts, references; no new experiment"},
        ]},
    ],
    "deduplication": [
        "Table 1 entry 1 / Scheme 2 product 3a / page-1 prose are one experiment.",
        "Table 2 entry 1 / Scheme 3 product 4a are one experiment.",
    ],
    "limitations": ["No SI supplied.", "Compact products are not forced into complete SMILES.", "Scope/control yield methods are unreported unless a table footnote explicitly states them."],
}
write(RUN / "raw_model/page_extraction.json", page_inventory)

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 1,
            "state": "attempted",
            "regions": [
                {"region": "title, abstract, authors, DOI, and ESI notice", "finding": "Metadata and abstract; ESI is advertised but no SI file was supplied or associated.", "evidence_ids": EV[1]},
                {"region": "right-column current-work discussion and first standard-result paragraph", "finding": "Current Fe/Cu chemodivergent transformation and the beginning of the 1a→3a standard experiment (81%); linked to the generic scheme and Table 1 entry 1 rather than duplicated.", "evidence_ids": EV[1]},
                {"region": "introduction and biological/method context", "finding": "Background and cited methods; no additional current-paper concrete reaction arrow.", "evidence_ids": EV[1]},
            ],
            "record_ids": [f"{DOC}:scheme1_current_generic", f"{DOC}:table1_entry1"],
            "unresolved": ["No supplied SI; experimental workup and purification details are unavailable."]
        },
        {
            "page": 2,
            "state": "attempted",
            "regions": [
                {"region": "Scheme 1", "finding": "One grouped cited prior-work family and one current two-branch generic Fe/Cu scheme.", "evidence_ids": EV[2]},
                {"region": "Table 1 entries 1–13 and footnotes", "finding": "Thirteen independent Fe-system optimization experiments; entry 1 duplicates Scheme 2 product 3a. Table/prose conflict retained for Co(acac)3 versus Co(acac)2 in entry 6.", "evidence_ids": EV[2]},
                {"region": "Scheme 2a substrate scope", "finding": "Sixteen printed products 3a–3p; 3a merged with Table 1 entry 1, leaving fifteen additional scope records.", "evidence_ids": EV[2]},
                {"region": "Scheme 2b amplification reaction", "finding": "One 4.0 mmol 1a scale-up giving 3a in 70% (0.69 g).", "evidence_ids": EV[2]},
            ],
            "record_ids": ([f"{DOC}:scheme1_prior_work", f"{DOC}:scheme1_current_generic"] +
                           [f"{DOC}:table1_entry{i}" for i in range(1, 14)] +
                           [f"{DOC}:scheme2_3{c}" for c in "bcdefghijklmnop"] +
                           [f"{DOC}:scheme2_scale_up_3a"]),
            "unresolved": ["Scheme 2 yields and the scale-up yield lack an explicit measurement basis.", "Compact product structures and several internal-arene substrate positions remain unencoded."]
        },
        {
            "page": 3,
            "state": "attempted",
            "regions": [
                {"region": "Table 2 entries 1–7 and footnotes", "finding": "Seven independent Cu-system optimization experiments, including two Trace and two NR results; entry 1 duplicates Scheme 3 product 4a.", "evidence_ids": EV[3]},
                {"region": "Scheme 3 substrate scope", "finding": "Thirteen printed products 4a–4m; 4a merged with Table 2 entry 1, leaving twelve additional scope records.", "evidence_ids": EV[3]},
                {"region": "Scheme 4a control experiment", "finding": "Two distinct TEMPO experiments: Fe system gives trace 3a, Cu system gives 4a in 79%.", "evidence_ids": EV[3]},
                {"region": "Scheme 4b and adjacent prose", "finding": "Author-proposed Fe-H radical cycle for 3; mechanism record, not executed intermediate experiments.", "evidence_ids": EV[3]},
            ],
            "record_ids": ([f"{DOC}:table2_entry{i}" for i in range(1, 8)] +
                           [f"{DOC}:scheme3_4{c}" for c in "bcdefghijklm"] +
                           [f"{DOC}:scheme4_tempo_fe", f"{DOC}:scheme4_tempo_cu", f"{DOC}:scheme4b_fe_mechanism"]),
            "unresolved": ["Scheme 3 and Scheme 4a yield measurement bases are not stated.", "No absolute/relative stereochemical assignments, ee, er, or dr are reported for products on this page."]
        },
        {
            "page": 4,
            "state": "attempted",
            "regions": [
                {"region": "Scheme 5 and mechanism discussion", "finding": "Author-proposed Cu-H cycle for formation of 4; no executed intermediate synthesis was created.", "evidence_ids": EV[4]},
                {"region": "conclusion, acknowledgements, conflict statement, and references", "finding": "No additional current-paper synthesis experiment, optimization row, or reaction arrow.", "evidence_ids": EV[4]},
                {"region": "whole-page spectral inventory", "finding": "No standalone mass-spectrum peak table, molecular-weight table, or fragment-relative-abundance block.", "evidence_ids": EV[4]},
            ],
            "record_ids": [f"{DOC}:scheme5_cu_mechanism"],
            "unresolved": ["Mechanistic intermediates A and D–F are proposed and not independently established on the supplied pages."]
        },
    ]
}
write(RUN / "coverage/pages_0001_0004.json", coverage)

print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]), "stereo_checks": len(stereo_checks)}, ensure_ascii=False))
