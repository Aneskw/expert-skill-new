#!/usr/bin/env python3
"""Build the page-1--5 evidence-linked v3.1 submission for the current paper.

This is a source-extraction artifact, not a reusable parser.  It intentionally
leaves structures and SI-only details unresolved when the five-page main text
does not establish them.
"""
from __future__ import annotations

import json
import sys
from copy import deepcopy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st  # noqa: E402


RUN = ROOT / "runs/current"
DATASET = RUN / "dataset.json"
SUBMITTED = RUN / "submitted.json"
DOC = "doc_776e16cc231b4049"
TITLE = "Regioselective Intramolecular Allene Amidation Enabled by an EDA Complex"
E = {p: f"{DOC}:page{p}:image" for p in range(1, 6)}
T = {p: f"{DOC}:page{p}:text" for p in range(1, 6)}


def F(raw=None, norm=None, ev=(), method="original_extraction", status=None, derivation=None):
    return st.fact(raw, norm, ev, method=method, status=status, derivation=derivation)


def inferred(norm, ev, derivation, raw=None):
    return F(raw, norm, ev, method="inference", derivation=derivation)


def derived(norm, ev, derivation, raw=None):
    return F(raw, norm, ev, method="structure_derivation", derivation=derivation)


def conflict(alternatives, ev):
    f = F([a[0] for a in alternatives], None, ev)
    f["data_status"] = "conflict"
    f["alternatives"] = [F(raw, norm, aev) for raw, norm, aev in alternatives]
    return f


data = st.read(DATASET)
mistyped_reactions = {
    f"{DOC}:table3_trapping_2aaa", f"{DOC}:table3_trapping_2aab",
    f"{DOC}:table3_trapping_2aac", f"{DOC}:table3_trapping_2aad",
}
mistyped_molecules = {
    f"{DOC}:mol:2aaa_mixture", f"{DOC}:mol:2aab_mixture",
    f"{DOC}:mol:2aac_mixture", f"{DOC}:mol:2aad_mixture",
}
# Remove only the four typo-bearing records produced by the immediately prior
# run of this same script.  This is not a general deletion or migration path.
data["experiments"] = [x for x in data["experiments"] if x.get("record_id") not in mistyped_reactions]
if "standard" in data:
    data["standard"]["reactions"] = [x for x in data["standard"]["reactions"] if x.get("reaction_id") not in mistyped_reactions]
    data["standard"]["molecules"] = [x for x in data["standard"]["molecules"] if x.get("molecule_id") not in mistyped_molecules]
if data["experiments"]:
    assert "standard" in data and all(
        x.get("document_id") == DOC and x.get("status") == "needs_review" for x in data["experiments"]
    ), "Refuse to use or overwrite unrelated or reviewed experiment records"
data["documents"][0]["title"] = TITLE
st.write(DATASET, data)
bundle = st.read(SUBMITTED)
bundle["documents"] = deepcopy(data["documents"])
bundle["source_context"] = []
molecules = {}
reactions = []
page_records = {p: [] for p in range(1, 6)}


def _set_stereo(m, ev, smiles=None, kind="ordinary"):
    if smiles:
        diag = st.stereo_inspect(smiles)
        m["stereochemistry"]["isomeric_smiles"] = derived(
            diag["canonical_isomeric_smiles"], ev,
            "RDKit graph diagnostic on the candidate SMILES; this does not verify the source drawing."
        )
        m["stereochemistry"]["identity_key"] = derived(
            diag["stereo_graph_key"], ev,
            "SHA-256 stereo graph key calculated from the candidate isomeric SMILES."
        )
    if kind == "allene_unspecified":
        m["stereochemistry"]["axial"] = F(
            "allene drawing has no R_a/S_a or P/M assignment", [{"type": "allene", "configuration": "unspecified"}], ev
        )
        m["stereochemistry"]["chirality_status"] = inferred(
            "unspecified", ev, "The substituted allene can require axial-stereochemistry review, but the source gives no assignment."
        )
        m["stereochemistry"]["geometry_status"] = inferred("not_applicable", ev, "Allene axial stereochemistry is not E/Z geometry.")
        m["stereochemistry"]["representation_status"] = inferred(
            "not_encoded", ev, "Candidate SMILES does not certify allene axial configuration."
        )
    elif kind == "ez_mixture":
        m["stereochemistry"]["double_bonds"] = F(
            "product is reported as E/Z (or Z/E) isomers; absolute structures are not separately assigned",
            [{"configuration": "E/Z mixture", "individual_isomer_identity": "unassigned"}], ev
        )
        m["stereochemistry"]["mixture"] = F(
            "E/Z mixture; experiment-specific ratios are stored on reactions", "geometric-isomer mixture", ev
        )
        m["stereochemistry"]["chirality_status"] = inferred(
            "unspecified", ev, "The cyclic product drawing does not assign any possible stereogenic center."
        )
        m["stereochemistry"]["geometry_status"] = F("E/Z mixture", "mixture", ev)
        m["stereochemistry"]["representation_status"] = inferred(
            "partial", ev, "Connectivity may be represented, but the reported E/Z mixture is not a single isomeric SMILES."
        )
    else:
        m["stereochemistry"]["chirality_status"] = inferred(
            "unresolved", ev, "Stereochemistry was reviewed at the source label/drawing level; no positive assignment is asserted."
        )
        m["stereochemistry"]["geometry_status"] = inferred(
            "unresolved", ev, "No E/Z or ring cis/trans assignment is asserted for this registry entry."
        )
        m["stereochemistry"]["representation_status"] = inferred(
            "unresolved" if not smiles else "partial", ev,
            "The machine record does not claim that all source stereochemistry is fully encoded."
        )


def mol(label, ev, name=None, smiles=None, stereo="ordinary"):
    mid = f"{DOC}:mol:{label.replace(' ', '_')}"
    ev = list(dict.fromkeys(ev))
    if mid in molecules:
        m = molecules[mid]
        occ = m["occurrences"][0]
        occ["evidence_ids"] = list(dict.fromkeys(occ["evidence_ids"] + ev))
        occ["compound_label"]["evidence_ids"] = list(dict.fromkeys(occ["compound_label"]["evidence_ids"] + ev))
        return mid
    m = st.new_molecule(mid, DOC, ev)
    m["occurrences"][0]["compound_label"] = F(label, label, ev)
    m["names"]["original"] = F(name or label, name or label, ev)
    m["identity_note"] = F(
        "Identity is limited to the main-text label/drawing; no SI or external database was supplied.",
        "main-text identity; machine record needs review", ev
    )
    if smiles:
        m["identity_status"] = "source_identified"
        m["representation_type"] = "smiles"
        m["smiles"] = derived(
            smiles, ev,
            "Candidate connectivity was transcribed from the displayed structure/name and checked for SMILES parsing; source identity remains needs_review.",
            raw={"source": "displayed structure or explicit chemical name", "candidate_smiles": smiles},
        )
    else:
        m["identity_status"] = "unresolved"
        m["representation_type"] = "unresolved"
        m["identity_note"] = F(
            "Main-text label/drawing retained; a complete stereochemically reliable SMILES was not established from pages 1-5.",
            "unresolved main-text identity", ev
        )
    _set_stereo(m, ev, smiles, stereo)
    molecules[mid] = m
    return mid


def generic_mol(label, ev, fixed, variables):
    mid = f"{DOC}:mol:{label.replace(' ', '_')}"
    if mid in molecules:
        return mid
    m = st.new_molecule(mid, DOC, ev)
    m["occurrences"][0]["compound_label"] = F(label, label, ev)
    m["names"]["original"] = F(label, label, ev)
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": F(fixed, fixed, ev),
        "variables": [
            {
                "label": v,
                "attachment_points": F(f"source drawing attachment for {v}", f"source drawing attachment for {v}", ev),
                "definition": F(None, None, ev),
                "allowed_values": F(None, None, ev),
                "constraints": F(None, None, ev),
            }
            for v in variables
        ],
        "instance_mapping": [],
        "reason": F("R-group/general-scheme structure", "Generic drawing is not one concrete molecule.", ev),
    }
    m["identity_note"] = F("Generic structure in the source figure.", "generic source structure", ev)
    _set_stereo(m, ev)
    molecules[mid] = m
    return mid


def reaction(rid, ev, record_type="specific_experiment", categories=(), page=None):
    ev = list(dict.fromkeys(ev))
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = F(list(categories), list(categories), ev)
    r["location"] = F(
        {"pages": sorted({int(x.split(":page")[1].split(":")[0]) for x in ev}), "source_region": rid},
        {"pages": sorted({int(x.split(":page")[1].split(":")[0]) for x in ev}), "source_region": rid}, ev
    )
    r["review"] = {
        "status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
        "unresolved": [
            "No supporting information was explicitly supplied; SI-only quantities, procedures, workup, purification, and characterization remain unreported.",
            "Machine extraction and visual inspection do not establish chemical correctness or experimental reproducibility."
        ],
    }
    reactions.append(r)
    if page:
        page_records[page].append(r["reaction_id"])
    return r


def participant(r, pid, mid, role, ev, eq=None, molpct=None, product_role=None, eq_ref=None):
    p = st.new_participant(pid, mid, role, ev)
    if product_role is not None:
        p["product_role"] = F(product_role, product_role, ev)
    if eq is not None:
        p["charges"]["equivalents"] = F(f"{eq} equiv", {"number": eq, "unit": "equiv"}, ev)
        p["equivalence_reference"] = F(eq_ref, eq_ref, ev)
    if molpct is not None:
        p["charges"]["mol_percent"] = F(f"{molpct} mol%", {"number": molpct, "unit": "mol%"}, ev)
        p["equivalence_reference"] = F(eq_ref, eq_ref, ev)
    r["participants"].append(p)
    return p


def participant_conflict_eq(r, pid, mid, role, ev, alts, eq_ref):
    p = participant(r, pid, mid, role, ev)
    p["charges"]["equivalents"] = conflict(alts, ev)
    p["equivalence_reference"] = F(eq_ref, eq_ref, ev)
    return p


def solvent(r, pid, mid, ev, name, concentration=None):
    participant(r, pid, mid, "solvent", ev)
    r["solvent_conditions"]["composition"] = F(name, [name], ev)
    if concentration is not None:
        r["solvent_conditions"]["concentration"] = F(
            f"{concentration} M", {"number": concentration, "unit": "M"}, ev
        )
        r["solvent_conditions"]["concentration_basis"] = F(
            "table/scheme reaction concentration", "reaction mixture concentration (basis compound not explicitly stated)", ev
        )


def step(r, ev, temp=None, time=None, operation="irradiate"):
    s = st.new_step(f"step_{len(r['steps'])+1}", ev)
    s["participant_ids"] = [p["participant_id"] for p in r["participants"]]
    s["operation"] = F(operation, operation, ev)
    if temp is not None:
        if isinstance(temp, dict):
            raw, norm = "-20 °C or r.t.", temp
        elif isinstance(temp, str):
            raw, norm = temp, temp
        else:
            raw, norm = f"{temp} °C", {"number": temp, "unit": "°C"}
        s["temperature"] = F(raw, norm, ev)
    if time is not None:
        n, u = time
        s["time"] = F(f"{n} {u}", {"number": n, "unit": u}, ev)
    r["steps"].append(s)


def light(r, ev, source, wavelength=None, power=None):
    r["special_conditions"]["light_source"] = F(source, source, ev)
    if wavelength is not None:
        r["special_conditions"]["wavelength"] = F(f"{wavelength} nm", {"number": wavelength, "unit": "nm"}, ev)
    if power is not None:
        r["special_conditions"]["power"] = F(f"{power} W", {"number": power, "unit": "W"}, ev)


def outcome(r, oid, kind, targets, ev, raw, norm, yield_type=None, basis=None, scope=None):
    o = st.new_outcome(oid, kind, targets, ev)
    o["value"] = F(raw, norm, ev)
    if yield_type:
        o["yield_type"] = F(yield_type, yield_type, ev)
    if basis:
        o["measurement_basis"] = F(basis, basis, ev)
    if scope:
        o["scope"] = F(scope, scope, ev)
    if kind in ("yield", "product_amount"):
        o["product_role"] = F("reported product", "reported product", ev)
    r["outcomes"].append(o)
    return o


def manual_classification(r, ev, specs, coverage="partial", limitation=None):
    changes = []
    labels = {}
    for i, (kind, a, b, before, after) in enumerate(specs, 1):
        change = {
            "change_id": f"b{i}", "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
            "before_order": before, "after_order": after, "change_type": kind,
        }
        changes.append(change)
        if kind in ("formed", "order_increased", "hydrogen_gain"):
            lab = st.bond_label(a, b, after)
            labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(change["change_id"])
    deriv = "Net bond changes inferred from the displayed substrate/product connectivity; atom mapping is unresolved."
    r["classification"]["bond_changes"] = [derived(c, ev, deriv) for c in changes]
    r["classification"]["labels"] = derived(list(labels.values()), ev, deriv)
    r["classification"]["coverage"] = inferred(coverage, ev, "Classification is partial because full atom mapping and coproduct identities are unavailable.")
    r["classification"]["limitations"] = [
        limitation or "No complete atom mapping; implicit-H and leaving-group provenance are not fully established."
    ]


def no_reaction_classification(r, ev):
    r["classification"]["bond_changes"] = []
    r["classification"]["labels"] = inferred([], ev, "No observed target conversion; expected-product bonds were not assigned.")
    r["classification"]["coverage"] = inferred("not_applicable", ev, "No observed reaction to classify.")
    r["classification"]["limitations"] = ["Negative control; no expected product bond labels inferred."]


def unresolved_classification(r, ev, why):
    r["classification"]["coverage"] = inferred("unresolved", ev, why)
    r["classification"]["limitations"] = [why]


def note(r, text, ev):
    r["notes"].append(F(text, text, ev))


# Core structures and named reagents used repeatedly.
m1a = mol("1a", [E[1], E[2], E[3], E[4], E[5]], "aryloxyamide allene 1a",
          "CC=C=CCC(=O)N(C)Oc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1", "allene_unspecified")
m2a = mol("2a_EZ_mixture", [E[1], E[2], E[3], E[4], E[5]], "gamma-lactam 2a (E/Z mixture)",
          "CN1C(=O)CCC1C=CC", "ez_mixture")
m3a = mol("3a", [E[2], E[4], E[5]], "H-atom abstraction product 3a", "CC=C=CCC(=O)NC", "allene_unspecified")
k2co3 = mol("K2CO3", [E[1], E[2], E[3], E[4], E[5]], "potassium carbonate", "[K+].[K+].[O-]C([O-])=O")
chd = mol("1,4-CHD", [E[1], E[2], E[3], E[4]], "1,4-cyclohexadiene", "C1=CCCC=C1")
acetone = mol("acetone", [E[1], E[2], E[3], E[4]], "acetone", "CC(=O)C")
me_cn = mol("MeCN", [E[2], E[3]], "acetonitrile", "CC#N")
acetone_d6 = mol("acetone-d6", [E[2]], "acetone-d6", "[2H]C([2H])([2H])C(=O)C([2H])([2H])[2H]")
eosin = mol("Eosin_Y", [E[1], E[2]], "Eosin Y")
dha = mol("9,10-DHA", [E[2]], "9,10-dihydroanthracene", "c1ccc2c(c1)Cc1ccccc1C2")
hantzsch = mol("Hantzsch_ester", [E[2]], "Hantzsch ester", "CCOC(=O)C1=C(C)NC(C)=C(C(=O)OCC)1")


# Page 1: general/prior-work arrows and this-work overview.
gsub = generic_mol("Scheme1A_generic_substrate", [E[1]], "amide tethered to alkene/alkyne with N-RPG", ["R", "X", "RPG"])
gprod = generic_mol("Scheme1A_generic_product", [E[1]], "cyclized N-heterocycle radical/product", ["R", "X"])
r = reaction("scheme1A_general_prior", [E[1]], "generic_scheme", (), 1)
r["reaction_name"] = F("General Scheme for addition of amidyl radicals to alkenes and alkynes", "addition of amidyl radicals to alkenes and alkynes", [E[1]])
participant(r, "generic_substrate", gsub, "reactant", [E[1]])
participant(r, "generic_product", gprod, "product", [E[1]], product_role="generic cyclized product")
step(r, [E[1]], operation="initiator or proton-coupled electron transfer under hν/Δ")
unresolved_classification(r, [E[1]], "Generic R-group scheme has no single atom-mapped transformation.")

ksub = generic_mol("Knowles2015_substrate", [E[1]], "Knowles 2015 generic N-radical cyclization substrate", ["Ar", "R1", "R2", "X"])
kprod = generic_mol("Knowles2015_product", [E[1]], "Knowles 2015 generic cyclized product", ["Ar", "R1", "R2", "X"])
kir = mol("Knowles_Ir_catalyst", [E[1]], "Ir(dF(CH3)ppy)2(bpy)PF6")
kph = mol("PhSH", [E[1]], "thiophenol", "Sc1ccccc1")
kadd = mol("NBu4OP(O)(OBu)2", [E[1]], "NBu4OP(O)(OBu)2")
dcm = mol("CH2Cl2", [E[1]], "dichloromethane", "ClCCl")
r = reaction("scheme1A_knowles_2015", [E[1]], "literature_example", (), 1)
r["reaction_name"] = F("Knowles and co-workers (2015)", "prior-work photocatalytic cyclization", [E[1]])
participant(r, "substrate", ksub, "reactant", [E[1]])
participant(r, "product", kprod, "product", [E[1]], product_role="literature-example product")
participant(r, "ir", kir, "catalyst", [E[1]], molpct=2, eq_ref=ksub)
participant(r, "additive", kadd, "additive", [E[1]], molpct=20, eq_ref=ksub)
participant(r, "phsh", kph, "reagent", [E[1]], molpct=10, eq_ref=ksub)
solvent(r, "dcm", dcm, [E[1]], "CH2Cl2", 0.3)
step(r, [E[1]], temp="r.t.", operation="irradiate with blue LEDs")
outcome(r, "yield_range", "yield", [kprod], [E[1]], "68-94%", {"minimum": 68, "maximum": 94, "unit": "%"})
unresolved_classification(r, [E[1]], "Prior-work generic structures contain R groups; not classified as this paper's experiment.")

lsub = generic_mol("Leonori2016_substrate", [E[1]], "Leonori 2016 aryloxyamide substrate", ["R1", "R2", "R3", "X"])
lprod = generic_mol("Leonori2016_product", [E[1]], "Leonori 2016 cyclized product", ["R1", "R2", "R3", "X"])
r = reaction("scheme1A_leonori_2016", [E[1]], "literature_example", (), 1)
r["reaction_name"] = F("Leonori and co-workers (2016)", "prior-work Eosin Y photocatalytic cyclization", [E[1]])
participant(r, "substrate", lsub, "reactant", [E[1]])
participant(r, "product", lprod, "product", [E[1]], product_role="literature-example product")
participant(r, "eosin", eosin, "catalyst", [E[1]], molpct=2, eq_ref=lsub)
participant(r, "chd", chd, "reagent", [E[1]], eq=2, eq_ref=lsub)
participant(r, "k2co3", k2co3, "additive", [E[1]], eq=2, eq_ref=lsub)
solvent(r, "acetone", acetone, [E[1]], "acetone")
r["solvent_conditions"]["concentration"] = F("0.05-0.1 M", {"minimum": 0.05, "maximum": 0.1, "unit": "M"}, [E[1]])
step(r, [E[1]], temp="r.t.", operation="irradiate with green LEDs")
outcome(r, "yield_range", "yield", [lprod], [E[1]], "29-79%", {"minimum": 29, "maximum": 79, "unit": "%"})
unresolved_classification(r, [E[1]], "Prior-work generic structures contain R groups; not classified as this paper's experiment.")

this_sub = generic_mol("Scheme1B_substrate", [E[1]], "this-work generic allene substrate", ["R1", "R2", "R3", "R4", "X"])
this_prod = generic_mol("Scheme1B_product", [E[1]], "this-work generic trapped cyclized product", ["R1", "R2", "R3", "R4", "X", "Y"])
trap_y = generic_mol("trapping_agent_Y", [E[1]], "trapping agent Y", ["Y"])
r = reaction("scheme1B_this_work_general", [E[1]], "generic_scheme", (), 1)
r["reaction_name"] = F("Regioselective addition of amidyl radicals to allenes", "intramolecular allene amidation", [E[1]])
participant(r, "substrate", this_sub, "reactant", [E[1]])
participant(r, "trap", trap_y, "reagent", [E[1]])
participant(r, "k2co3", k2co3, "additive", [E[1]], eq=2, eq_ref=this_sub)
participant(r, "product", this_prod, "product", [E[1]], product_role="generic trapped cyclization product")
solvent(r, "acetone", acetone, [E[1]], "acetone")
step(r, [E[1]], temp="variable cold or room temperature", operation="irradiate with blue LEDs")
r["steps"][-1]["temperature"] = F("-20 °C or r.t.", "variable cold or room temperature", [E[1]])
outcome(r, "yield_range", "yield", [this_prod], [E[1]], "35-85%", {"minimum": 35, "maximum": 85, "unit": "%"})
outcome(r, "ez_range", "selectivity", [this_prod], [E[1]], "E:Z 1.1:1 to 19:1", {"metric": "E:Z", "range": ["1.1:1", "19:1"]}, scope="17 examples; isomer identities unassigned")
unresolved_classification(r, [E[1]], "Generic trapping agent Y gives multiple possible bond types; no single classification asserted.")


# Page 2 Table 1: every optimization row is a separate experiment.
opt_rows = [
    (1, "Green LED", None, 13, "acetone", 0.1, "1,4-CHD", 21, 41, "2.1:1", "rt"),
    (2, "Green LED", None, 13, "acetone", 0.1, "1,4-CHD", 21, 18, "2.2:1", "rt"),
    (3, "440 nm", 440, 40, "acetone", 0.1, "1,4-CHD", 5, 76, "2.2:1", "rt"),
    (4, "440 nm", 440, 40, "acetone", 0.025, "1,4-CHD", 1.5, 78, "2.1:1", "rt"),
    (5, "440 nm", 440, 40, "MeCN", 0.025, "1,4-CHD", 3, 71, "2.1:1", "rt"),
    (6, "440 nm", 440, 40, "acetone", 0.025, "9,10-DHA", 4, 71, "2.1:1", "rt"),
    (7, "440 nm", 440, 40, "acetone", 0.025, "Hantzsch ester", 1, 74, "2.3:1", "rt"),
    (8, "no light", None, None, "acetone", 0.025, "1,4-CHD", 1.5, 0, None, "rt"),
    (9, "440 nm", 440, 40, "acetone", 0.025, "1,4-CHD", 24, 0, None, "rt"),
    (10, "440 nm", 440, 40, "acetone-d6", 0.025, "conflict", 24, 49, "2:1", "rt"),
    (11, "440 nm", 440, 40, "acetone", 0.025, "conflict", 16, 0, None, "rt"),
    (12, "440 nm", 440, 40, "acetone", 0.025, "1,4-CHD", 18, 85, "3.2:1", -20),
    (13, "440 nm", 440, 40, "acetone", 0.025, "1,4-CHD", 92, 83, "3.6:1", -40),
]
donors = {"1,4-CHD": chd, "9,10-DHA": dha, "Hantzsch ester": hantzsch}
solvents = {"acetone": acetone, "MeCN": me_cn, "acetone-d6": acetone_d6}
for ent, src, wl, pw, solv_name, conc, donor_name, hours, yld, ez, temp in opt_rows:
    ev = [E[2], T[2]]
    r = reaction(f"table1_entry_{ent:02d}", ev, "specific_experiment", ("optimization", "control") if ent in (8, 9, 10, 11) else ("optimization",), 2)
    r["reaction_name"] = F("Optimization of amidyl radical addition to allenes", "intramolecular amidyl-radical allene cyclization", ev)
    r["transformation_type"] = inferred("5-exo-trig allene amidation", ev, "Named from the author's discussion of Table 1.")
    participant(r, "substrate_1a", m1a, "reactant", ev)
    participant(r, "product_2a", m2a, "product", ev, product_role="target cyclization product")
    participant(r, "product_3a", m3a, "product", ev, product_role="competing H-atom abstraction product")
    if ent == 9:
        note(r, "Table 1 footnote [c]: no K2CO3.", ev)
    elif ent == 11:
        participant_conflict_eq(
            r, "k2co3", k2co3, "additive", ev,
            [("Table 1 general condition: 2 equiv K2CO3", {"number": 2, "unit": "equiv"}, [E[2]]),
             ("prose: removing both 1,4-CHD and K2CO3", {"number": 0, "unit": "equiv"}, [T[2]])], m1a
        )
        r["review"]["unresolved"].append("Entry 11 conflicts: Table 1 general header implies K2CO3, while prose says K2CO3 was removed.")
    else:
        participant(r, "k2co3", k2co3, "additive", ev, eq=2, eq_ref=m1a)
    if donor_name == "conflict":
        p = participant_conflict_eq(
            r, "h_donor", chd, "reagent", ev,
            [("Table 1 H-donor column: 1,4-CHD", {"number": 2, "unit": "equiv"}, [E[2]]),
             ("prose: absence/removal of 1,4-CHD", {"number": 0, "unit": "equiv"}, [T[2]])], m1a
        )
        r["review"]["unresolved"].append(f"Entry {ent} conflicts on whether 1,4-CHD was present; both table and prose retained.")
    else:
        participant(r, "h_donor", donors[donor_name], "reagent", ev, eq=2, eq_ref=m1a)
    if ent == 1:
        participant(r, "eosin_y", eosin, "catalyst", ev, molpct=2, eq_ref=m1a)
    solvent(r, "solvent", solvents[solv_name], ev, solv_name, conc)
    step(r, ev, temp=temp, time=(hours, "h"), operation="irradiate" if src != "no light" else "hold without light")
    light(r, ev, src, wl, pw)
    r["procedure"]["original_text"] = F({"Table 1 entry": ent}, {"Table 1 entry": ent}, ev)
    outcome(r, "yield_2a", "yield", [m2a], ev, f"{yld}%", {"number": yld, "unit": "%"},
            yield_type="1H NMR yield of E/Z isomers", basis="1H NMR with 1,3,5-trimethoxybenzene as internal standard")
    if ez:
        outcome(r, "selectivity_ez", "selectivity", [m2a], ev, f"E:Z {ez}", {"metric": "E:Z", "ratio": ez},
                scope="2a E/Z isomers; absolute isomer assignments not supplied")
    else:
        outcome(r, "no_target_reaction", "no_reaction", [], ev, "0% 2a; table dash for E:Z", True)
    if ent == 10:
        outcome(r, "deuterium_incorporation", "other", [m2a], ev, "3:1 D:H ratio (footnote [e]); prose states 37% D-incorporation",
                {"table_D:H": "3:1", "prose_D_incorporation": {"number": 37, "unit": "%"}}, scope="acetone-d6 experiment")
    if ent in (12, 13):
        outcome(r, "product_3a_not_observed", "other", [m3a], ev, "no H-atom abstraction product 3a observed", "not observed")
    if yld == 0:
        no_reaction_classification(r, ev)
    else:
        manual_classification(r, ev, [("formed", "C", "N", 0, 1), ("hydrogen_gain", "C", "H", 0, 1), ("broken", "N", "O", 1, 0), ("order_decreased", "C", "C", 2, 1)])


# Page 3 Table 2: 13 substrate-scope arrows.
scope_rows = [
    ("a", 69, "3:1", -20, "isolated yield"), ("b", 65, None, -20, None),
    ("c", 83, "2.4:1", -20, None), ("d", 60, "1.1:1", -20, None),
    ("e", 55, "9:1", -20, None), ("f", 35, "19:1", "r.t.", None),
    ("g", 43, "1.6:1", -20, None), ("h", 52, "1.5:1", "r.t.", None),
    ("i", 76, "3:1", -20, None), ("j", 39, "2.2:1", "r.t.", None),
    ("k", 64, "2.7:1", "r.t.", "1H NMR yield using 1,3,5-trimethoxybenzene as internal standard"),
    ("l", 35, "2.7:1", "r.t.", None), ("m", 40, None, "r.t.", None),
]
for letter, yld, ez, temp, basis in scope_rows:
    ev = [E[3], T[3]] + ([T[2]] if letter in ("a", "b", "c") else [])
    sub = m1a if letter == "a" else mol(f"1{letter}", ev, f"allene substrate 1{letter}")
    prod = m2a if letter == "a" else mol(f"2{letter}_mixture", ev, f"scope product 2{letter} (geometry as shown/reported)", None, "ez_mixture" if ez else "ordinary")
    r = reaction(f"table2_scope_1{letter}_to_2{letter}", ev, "specific_experiment", ("substrate_scope",), 3)
    if letter in ("a", "b", "c"):
        page_records[2].append(r["reaction_id"])
    r["reaction_name"] = F("Scope of intramolecular addition of N-radicals to allenes", "intramolecular amidyl-radical allene cyclization", ev)
    participant(r, "substrate", sub, "reactant", ev)
    participant(r, "product", prod, "product", ev, product_role="scope cyclization product")
    participant(r, "chd", chd, "reagent", ev, eq=2, eq_ref=sub)
    participant(r, "k2co3", k2co3, "additive", ev, eq=2, eq_ref=sub)
    solvent(r, "acetone", acetone, ev, "acetone", 0.025)
    step(r, ev, temp=temp, operation="irradiate; reaction time not reported in main-text table")
    light(r, ev, "Kessil LED", 440, 40)
    outcome(r, "yield", "yield", [prod], ev, f"{yld}%", {"number": yld, "unit": "%"},
            yield_type="reported yield", basis=basis)
    if ez:
        outcome(r, "selectivity_ez", "selectivity", [prod], ev, f"E:Z {ez}", {"metric": "E:Z", "ratio": ez},
                scope=f"2{letter}; individual E and Z structures/configurations not assigned")
    if letter == "m":
        manual_classification(r, ev, [("formed", "C", "N", 0, 1), ("order_increased", "C", "C", 1, 3), ("broken", "C", "Br", 1, 0)],
                              limitation="2m alkyne formation is attributed by the authors to radical bromine elimination; full mapping is unavailable.")
    else:
        manual_classification(r, ev, [("formed", "C", "N", 0, 1), ("hydrogen_gain", "C", "H", 0, 1), ("broken", "N", "O", 1, 0), ("order_decreased", "C", "C", 2, 1)])
    if sub != m1a or prod != m2a:
        r["review"]["unresolved"].append(f"Complete stereochemically reliable structures for 1{letter}/2{letter} were not established from the compact Table 2 drawing; labels and outcomes are retained.")


# Text-reported scope failures/side reactions on pages 2-3.
tri_sub = mol("1-Me-3-iPr-trisubstituted_allene", [T[2]], "1-Me-3-iPr-trisubstituted allene (unlabeled scope substrate)")
r = reaction("text_scope_trisubstituted_no_cyclization", [T[2], E[2]], "specific_experiment", ("substrate_scope", "control"), 2)
participant(r, "substrate", tri_sub, "reactant", [T[2], E[2]])
r["reaction_name"] = F("irradiation under the reaction conditions", "attempted allene amidation", [T[2]])
step(r, [T[2], E[2]], operation="irradiate under the reaction conditions; details deferred to unavailable SI")
outcome(r, "no_cyclized_product", "no_reaction", [], [T[2]], "gave no cyclized product", True)
no_reaction_classification(r, [T[2], E[2]])

n0_sub = mol("n0_XCH2_substrate", [T[3], E[3]], "substrate 1 with n=0 and X=CH2")
n0_prod = mol("n0_H_abstraction_product", [T[3], E[3]], "unidentified H-abstraction product from n=0, X=CH2 substrate")
r = reaction("text_scope_n0_H_abstraction", [T[3], E[3]], "specific_experiment", ("substrate_scope", "control"), 3)
participant(r, "substrate", n0_sub, "reactant", [T[3], E[3]])
participant(r, "product", n0_prod, "product", [T[3], E[3]], product_role="H-abstraction product")
r["reaction_name"] = F("only H-abstraction was noted", "H-atom abstraction without productive cyclization", [T[3]])
step(r, [T[3], E[3]], operation="reaction conditions not detailed in main text; see unavailable SI")
outcome(r, "only_h_abstraction", "other", [n0_prod], [T[3]], "only H-abstraction was noted", "H-abstraction only")
manual_classification(r, [T[3], E[3]], [("formed", "N", "H", 0, 1), ("broken", "N", "O", 1, 0)], limitation="Product structure and hydrogen source are not supplied in the main text.")


# Page 3 Table 3: radical trapping experiments.
ncs = mol("NCS", [E[3]], "N-chlorosuccinimide", "O=C1CCC(=O)N1Cl")
pome3 = mol("P(OMe)3", [E[3]], "trimethyl phosphite", "COP(OC)OC")
brmal = mol("diethyl_bromomalonate", [E[3]], "diethyl bromomalonate", "CCOC(=O)C(Br)C(=O)OCC")
seph = mol("N-phenylseleno_phthalimide", [E[3]], "N-(phenylseleno)phthalimide", "O=C1N([Se]c2ccccc2)C(=O)c2ccccc21")
trap_rows = [
    ("aa", ncs, "NCS", 10, acetone, "acetone", "r.t.", 71, "Z:E", "2.4:1", "Cl"),
    ("ab", pome3, "P(OMe)3", 10, me_cn, "MeCN", "r.t.", 86, "Z:E", "1.3:1", "P"),
    ("ac", brmal, "diethyl bromomalonate", 6, acetone, "acetone", "r.t.", 60, "Z:E", "2.3:1", "Br"),
    ("ad", seph, "N-(phenylseleno)phthalimide", 10, acetone, "acetone", -20, None, None, None, "Se"),
]
for suffix, trap, trap_name, eqv, smid, sname, temp, yld, metric, ratio, atom in trap_rows:
    ev = [E[3], T[3]]
    product_code = f"2{suffix}"
    prod = mol(f"{product_code}_mixture", ev, f"trapped cyclization product {product_code} (geometric-isomer mixture)", None, "ez_mixture")
    r = reaction(f"table3_trapping_{product_code}", ev, "specific_experiment", ("derivatization",), 3)
    r["reaction_name"] = F("Intermolecular trapping of the vinyl radical intermediates", "allene cyclization and radical trapping", ev)
    participant(r, "substrate", m1a, "reactant", ev)
    participant(r, "trap", trap, "reagent", ev, eq=eqv, eq_ref=m1a)
    participant(r, "k2co3", k2co3, "additive", ev, eq=2, eq_ref=m1a)
    participant(r, "product", prod, "product", ev, product_role="trapped alkene product")
    solvent(r, "solvent", smid, ev, sname, 0.025)
    step(r, ev, temp=temp, operation="irradiate; reaction time not reported in Table 3")
    light(r, ev, "440 nm irradiation", 440)
    if suffix != "ad":
        outcome(r, "yield", "yield", [prod], ev, f"{yld}%", {"number": yld, "unit": "%"}, yield_type="reported yield")
        outcome(r, "selectivity", "selectivity", [prod], ev, f"{metric} {ratio}", {"metric": metric, "ratio": ratio},
                scope=f"{product_code}; individual geometric isomers unassigned")
    else:
        o = st.new_outcome("yield_conflict", "yield", [prod], ev)
        o["value"] = conflict(
            [("Table 3: 42%", {"number": 42, "unit": "%"}, [E[3]]),
             ("prose: 69% yield", {"number": 69, "unit": "%"}, [T[3]])], ev
        )
        o["product_role"] = F("2ad trapped product", "reported product", ev)
        r["outcomes"].append(o)
        o = st.new_outcome("selectivity_conflict", "selectivity", [prod], ev)
        o["value"] = conflict(
            [("Table 3: Z:E 3.3:1", {"metric": "Z:E", "ratio": "3.3:1"}, [E[3]]),
             ("prose: E:Z ratio 3.3:1", {"metric": "E:Z", "ratio": "3.3:1"}, [T[3]])], ev
        )
        o["scope"] = F("2ad; metric/order conflicts between table and prose", "unresolved E/Z order", ev)
        r["outcomes"].append(o)
        r["review"]["unresolved"].append("2ad conflict: Table 3 gives 42% and Z:E 3.3:1; prose gives 69% and E:Z 3.3:1.")
        note(r, "Source conflict retained without selecting either Table 3 or prose values for 2ad.", ev)
    specs = [("formed", "C", "N", 0, 1), ("broken", "N", "O", 1, 0), ("order_decreased", "C", "C", 2, 1), ("formed", "C", atom, 0, 1)]
    manual_classification(r, ev, specs)


# Page 4 Table 4: salt effects, all rows independent.
salt_defs = [
    (1, "KCl", "[K+].[Cl-]", 25, "h", 62, "2.1:1", 16, 2),
    (2, "KHCO3", "[K+].O=C([O-])O", 9, "h", 76, "2.2:1", 0, 2),
    (3, "K2SO4", "[K+].[K+].[O-]S(=O)(=O)[O-]", 24, "h", 0, None, 64, 2),
    (4, "KOtBu", "CC(C)(C)[O-].[K+]", 1.2, "h", 23, "2.2:1", 0, 2),
    (5, "Cs2CO3", "[Cs+].[Cs+].[O-]C([O-])=O", 50, "min", 72, "2.3:1", 0, 2),
    (6, "Li2CO3", "[Li+].[Li+].[O-]C([O-])=O", 72, "h", 60, "2.1:1", 10, 2),
    (7, "Na2CO3", "[Na+].[Na+].[O-]C([O-])=O", 5, "h", 32, "2.3:1", 9, 2),
    (8, "TBAOAc", "CCCC[N+](CCCC)(CCCC)CCCC.CC(=O)[O-]", 2, "h", 65, "2.2:1", 0, 2),
    (9, "K2CO3", "[K+].[K+].[O-]C([O-])=O", 5.6, "h", 83, "2.3:1", 0, 0.2),
]
for ent, salt_name, smi, tnum, tunit, yld, ez, recovery, salt_eq in salt_defs:
    ev = [E[4], T[4]]
    salt = k2co3 if salt_name == "K2CO3" else mol(salt_name, ev, salt_name, smi)
    r = reaction(f"table4_salt_entry_{ent:02d}", ev, "specific_experiment", ("optimization", "mechanistic_study", "control") if ent == 3 else ("optimization", "mechanistic_study"), 4)
    r["reaction_name"] = F("Investigation into the effect of the salt in the reaction", "salt-effect amidyl-radical cyclization", ev)
    participant(r, "substrate_1a", m1a, "reactant", ev)
    participant(r, "product_2a", m2a, "product", ev, product_role="target cyclization product")
    participant(r, "product_3a", m3a, "product", ev, product_role="H-atom abstraction product reported as not observed")
    participant(r, "salt", salt, "additive", ev, eq=salt_eq, eq_ref=m1a)
    participant(r, "chd", chd, "reagent", ev, eq=2, eq_ref=m1a)
    solvent(r, "acetone", acetone, ev, "acetone", 0.025)
    step(r, ev, temp="r.t.", time=(tnum, tunit), operation="irradiate")
    light(r, ev, "440 nm irradiation", 440)
    outcome(r, "yield_2a", "yield", [m2a], ev, f"{yld}%", {"number": yld, "unit": "%"},
            yield_type="1H NMR yield", basis="1H NMR with 1,3,5-trimethoxybenzene as internal standard")
    if ez:
        outcome(r, "selectivity_ez", "selectivity", [m2a], ev, f"E:Z {ez}", {"metric": "E:Z", "ratio": ez},
                scope="2a E/Z isomers; selectivity measurement method not explicitly stated")
    else:
        outcome(r, "no_cyclized_product", "no_reaction", [], ev, "0% 2a; table dash for E:Z", True)
    outcome(r, "recovery_1a", "recovery", [m1a], ev, f"1a {recovery}%", {"number": recovery, "unit": "%"},
            basis="1H NMR with 1,3,5-trimethoxybenzene as internal standard")
    outcome(r, "product_3a_not_observed", "other", [m3a], ev, "no 3a observed", "not observed")
    if yld == 0:
        no_reaction_classification(r, ev)
    else:
        manual_classification(r, ev, [("formed", "C", "N", 0, 1), ("hydrogen_gain", "C", "H", 0, 1), ("broken", "N", "O", 1, 0), ("order_decreased", "C", "C", 2, 1)])


# Page 4 Scheme 2: TEMPO radical trapping.
tempo = mol("TEMPO", [E[4]], "TEMPO", "CC1(C)CCCC(C)(C)N1[O]")
tempo_prod = mol("2a-T", [E[4]], "TEMPO-trapped product 2a-T")
r = reaction("scheme2_TEMPO_trapping", [E[4], T[4]], "specific_experiment", ("mechanistic_study",), 4)
r["reaction_name"] = F("Employing TEMPO as a radical trapping agent", "TEMPO radical-trapping experiment", [E[4]])
participant(r, "substrate", m1a, "reactant", [E[4]])
participant(r, "k2co3", k2co3, "additive", [E[4]], eq=2, eq_ref=m1a)
participant(r, "tempo", tempo, "reagent", [E[4]], eq=2, eq_ref=m1a)
participant(r, "product_2a", m2a, "product", [E[4]], product_role="cyclization product")
participant(r, "product_trap", tempo_prod, "product", [E[4]], product_role="TEMPO-trapped product")
solvent(r, "acetone", acetone, [E[4]], "acetone", 0.025)
step(r, [E[4]], temp="r.t.", time=(1.5, "h"), operation="irradiate")
light(r, [E[4]], "440 nm irradiation", 440)
outcome(r, "recovery_1a", "recovery", [m1a], [E[4]], "1a 14%", {"number": 14, "unit": "%"}, basis="1H NMR with 1,3,5-(MeO)3benzene internal standard")
outcome(r, "yield_2a", "yield", [m2a], [E[4]], "2a 11%", {"number": 11, "unit": "%"}, yield_type="1H NMR yield", basis="1H NMR with 1,3,5-(MeO)3benzene internal standard")
outcome(r, "yield_2aT", "yield", [tempo_prod], [E[4]], "2a-T 16%", {"number": 16, "unit": "%"}, yield_type="1H NMR yield", basis="1H NMR with 1,3,5-(MeO)3benzene internal standard")
manual_classification(r, [E[4]], [("formed", "C", "N", 0, 1), ("formed", "C", "O", 0, 1), ("broken", "N", "O", 1, 0), ("order_decreased", "C", "C", 2, 1)],
                      limitation="2a-T structure is retained by label; C-O trapping connectivity is inferred from the displayed TEMPO adduct.")


# Page 4 Scheme 3: N-O cleavage without an allene.
m6 = mol("6", [E[4]], "amidyl radical precursor 6", "CC(=O)N(C)Oc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1")
m7 = mol("7", [E[4]], "H-atom abstraction product 7", "CC(=O)NC")
r = reaction("scheme3_no_allene_cleavage", [E[4], T[4]], "specific_experiment", ("mechanistic_study", "control"), 4)
r["reaction_name"] = F("Cleavage of the N-O bond in the absence of an allene", "amidyl precursor N-O cleavage/H-atom abstraction", [E[4]])
participant(r, "substrate_6", m6, "reactant", [E[4]])
participant(r, "product_7", m7, "product", [E[4]], product_role="H-atom abstraction product")
participant(r, "chd", chd, "reagent", [E[4]], eq=2, eq_ref=m6)
participant(r, "k2co3", k2co3, "additive", [E[4]], eq=2, eq_ref=m6)
solvent(r, "acetone", acetone, [E[4]], "acetone", 0.025)
step(r, [E[4]], temp="r.t.", time=(1.5, "h"), operation="irradiate")
light(r, [E[4]], "440 nm irradiation", 440)
outcome(r, "yield_7", "yield", [m7], [E[4]], "88% 1H NMR yield", {"number": 88, "unit": "%"}, yield_type="1H NMR yield", basis="1H NMR; internal standard not stated for Scheme 3")
manual_classification(r, [E[4]], [("formed", "N", "H", 0, 1), ("broken", "N", "O", 1, 0)])


# Page 5 Scheme 4: proposed mechanism, not a performed experiment.
r = reaction("scheme4_proposed_mechanism", [E[5], T[5]], "proposed_mechanism", (), 5)
r["reaction_name"] = F("Mechanistic proposal", "proposed EDA/SET amidyl-radical cyclization mechanism", [E[5]])
participant(r, "substrate_1a", m1a, "reactant", [E[5]])
participant(r, "product_2a", m2a, "product", [E[5]], product_role="proposed cyclization product")
participant(r, "product_3a", m3a, "product", [E[5]], product_role="proposed H-abstraction product")
step(r, [E[5]], operation="author-proposed EDA complex formation, SET, N-O fragmentation, cyclization, and H-atom trapping")
no_reaction_classification(r, [E[5]])
r["classification"]["coverage"] = inferred("not_applicable", [E[5]], "Proposed elementary mechanism is not classified as one observed net experiment.")
note(r, "Scheme 4 is an author-proposed mechanism; intermediates A-F are not treated as isolated synthetic products.", [E[5]])


# Non-synthetic observations and unenumerated qualitative screens belong in source_context.
bundle["source_context"] = [
    {
        "context_id": f"{DOC}:metadata",
        "kind": "bibliographic_metadata",
        "title": TITLE,
        "authors": ["Lu Liu", "Robert M. Ward", "Jennifer M. Schomaker"],
        "doi": "10.1002/chem.202002533",
        "evidence_ids": [E[1], T[1]],
    },
    {
        "context_id": f"{DOC}:figure1_uvvis",
        "kind": "non_synthetic_observation",
        "description": "UV-vis spectra of K2CO3, 1a, and 1a/K2CO3 in acetone; retained as mechanistic source context, not a synthesis reaction.",
        "evidence_ids": [E[4], T[4]],
    },
    {
        "context_id": f"{DOC}:unresolved_other_traps",
        "kind": "qualitative_screen_summary",
        "description": "Prose says Michael acceptors, O2, and various heteroatom donors gave only 2a or decomposition; individual agents/conditions/outcomes are not enumerated and were not fabricated as separate experiments.",
        "evidence_ids": [T[3], E[3]],
    },
    {
        "context_id": f"{DOC}:lights_on_off",
        "kind": "mechanistic_observation",
        "description": "A lights on/off experiment is mentioned, but its numerical sequence is only in unavailable SI; no specific experiment record was invented.",
        "evidence_ids": [T[4]],
    },
]

bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
st.write(SUBMITTED, bundle)


raw = {
    "document_id": DOC,
    "title": TITLE,
    "pages_visually_inspected": [1, 2, 3, 4, 5],
    "record_count": len(reactions),
    "table_counts": {"Table 1": 13, "Table 2": 13, "Table 3": 4, "Table 4": 9},
    "conflicts": [
        "Table 3 2ad: 42%, Z:E 3.3:1 versus prose 69%, E:Z 3.3:1.",
        "Table 1 entries 10-11: H-donor column/general header versus prose absence/removal of 1,4-CHD and K2CO3."
    ],
    "limitations": [
        "No SI was supplied or associated.",
        "Compact scope structures not confidently transcribed remain unresolved.",
        "Machine records remain needs_review."
    ],
}
st.write(RUN / "raw_model/extraction_summary.json", raw)


coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 1, "state": "attempted",
            "regions": [
                {"region": "title, abstract, and bibliographic header", "finding": "paper identity and stated transformation; no specific experiment", "evidence_ids": [E[1], T[1]]},
                {"region": "Scheme 1A top", "finding": "generic previous-studies scheme, recorded as generic_scheme", "evidence_ids": [E[1]]},
                {"region": "Scheme 1A Knowles and Leonori arrows", "finding": "two cited prior-work literature examples, not experiments of this paper", "evidence_ids": [E[1]]},
                {"region": "Scheme 1B", "finding": "this-work generic scheme and reported aggregate scope", "evidence_ids": [E[1]]},
            ],
            "record_ids": page_records[1],
            "unresolved": ["Generic R groups prevent concrete molecule identities and complete bond mapping."]
        },
        {
            "page": 2, "state": "attempted",
            "regions": [
                {"region": "Table 1 rows 1-13 and footnotes", "finding": "13 optimization/control experiments recorded separately", "evidence_ids": [E[2], T[2]]},
                {"region": "prose surrounding Table 1", "finding": "control interpretations, acetone-d6 D incorporation, and low-temperature observations linked to rows", "evidence_ids": [T[2]]},
                {"region": "bottom-right scope introduction", "finding": "scope 1a-1c mentions linked to Table 2 records; trisubstituted allene no-cyclization control recorded", "evidence_ids": [E[2], T[2]]},
            ],
            "record_ids": list(dict.fromkeys(page_records[2])),
            "unresolved": ["Entries 10-11 conflict between table and prose on 1,4-CHD/K2CO3 presence.", "Wavelength screen and full solvent/H-donor screens are only in unavailable SI."]
        },
        {
            "page": 3, "state": "attempted",
            "regions": [
                {"region": "Table 2 arrows 1a-m to 2a-m and footnotes", "finding": "13 substrate-scope experiments recorded separately", "evidence_ids": [E[3], T[3]]},
                {"region": "scope prose", "finding": "n=0 H-abstraction-only experiment and qualitative discussion retained", "evidence_ids": [T[3]]},
                {"region": "Table 3 four trapping arrows", "finding": "four radical-trapping experiments recorded separately", "evidence_ids": [E[3], T[3]]},
                {"region": "other trapping-agent prose", "finding": "unnamed/unenumerated screens retained in source_context, not fabricated as experiments", "evidence_ids": [T[3]]},
            ],
            "record_ids": page_records[3],
            "unresolved": ["Compact Table 2 structures remain unresolved where a reliable stereochemical SMILES could not be established.", "2ad table/prose yield and E/Z-order conflict retained.", "Most scope yield measurement bases are not stated; only 2a isolated and 2k NMR are explicit."]
        },
        {
            "page": 4, "state": "attempted",
            "regions": [
                {"region": "Table 4 rows 1-9 and footnotes", "finding": "nine salt-effect mechanistic optimization experiments recorded separately", "evidence_ids": [E[4], T[4]]},
                {"region": "Scheme 2", "finding": "TEMPO trapping experiment with recovery and two product yields", "evidence_ids": [E[4], T[4]]},
                {"region": "Scheme 3", "finding": "N-O cleavage/H-abstraction control without an allene", "evidence_ids": [E[4], T[4]]},
                {"region": "Figure 1", "finding": "UV-vis observation stored as source_context, not synthesis", "evidence_ids": [E[4], T[4]]},
            ],
            "record_ids": page_records[4],
            "unresolved": ["UV-vis numerical traces were not digitized.", "Lights on/off numerical data are in unavailable SI."]
        },
        {
            "page": 5, "state": "attempted",
            "regions": [
                {"region": "Scheme 4", "finding": "author-proposed EDA/SET mechanism recorded as proposed_mechanism, not a performed experiment", "evidence_ids": [E[5], T[5]]},
                {"region": "conclusion, acknowledgements, references", "finding": "no additional specific synthetic experiment", "evidence_ids": [E[5], T[5]]},
            ],
            "record_ids": page_records[5],
            "unresolved": ["Intermediates A-F are proposed and not treated as isolated molecules or validated elementary steps."]
        },
    ],
}
(RUN / "coverage").mkdir(exist_ok=True)
st.write(RUN / "coverage/pages_0001_0005.json", coverage)


workflow = f"""# Workflow log

- Batch: `allene-stereo-v31-200-20261017`
- Document: `{DOC}`; DOI `10.1002/chem.202002533`
- Title: {TITLE}
- Actual visual coverage: original PDF pages 1-5, all five rendered page images inspected.
- Existing experiments checked before extraction: 0; no old experiment record was read as a new result.
- SI association: unestablished; no SI file was supplied, so SI-only procedures/details were not inferred.
- Formal extraction: {len(reactions)} v3.1 reactions, {len(molecules)} molecule registry entries, plus 4 source-context observations.
- Row/arrow accounting: Scheme 1 has 4 records; Table 1 has 13 rows; Table 2 has 13 arrows; Table 3 has 4 arrows; two prose-only scope/control observations are separate; Table 4 has 9 rows; Schemes 2 and 3 have one mechanistic experiment each; Scheme 4 is a proposed mechanism.
- Stereochemistry: each molecule carries explicit chirality/geometry/representation status. 1a/3a retain unspecified allene axial status; reported product E/Z mixtures retain mixture/partial representation and experiment-specific selectivity. No E or Z absolute identity was guessed.
- Structure handling: candidate SMILES were used only where the main-text drawing/name was sufficiently clear; compact Table 2 identities and several trapped products remain unresolved rather than guessed. All machine identities remain needs_review.
- Classification: partial net changes were inferred from displayed connectivity; no automatic atom mapping or model votes were invented. Failed controls use not_applicable classification.
- Conflict retained: Table 3 shows 2ad as 42%, Z:E 3.3:1, while prose states 69%, E:Z 3.3:1. Table 1 entries 10-11 conflict between displayed general/H-donor conditions and prose removal/absence statements.
- Measurement bases: Table 1 and Table 4 yields are 1H NMR with 1,3,5-trimethoxybenzene internal standard; 2a scope is explicitly isolated in prose; 2k is explicitly NMR by Table 2 footnote; other scope/trapping bases remain unreported.
- Non-synthetic material: Figure 1 UV-vis and the unquantified lights on/off mention are source_context, not synthetic experiments. Unenumerated failed trapping agents are not expanded into fabricated records.
- Audit disposition: no reaction_identity or outcome_target warnings remain. structure_gap warnings for compact Table 2 structures are intentionally unresolved rather than guessed. measurement_basis warnings are retained where the main text does not state a yield method; null/missing fields are not replaced by assumptions.
- Validation is schema/consistency checking only; semantic accuracy remains false and records remain needs_review.
"""
(RUN / "workflow.md").write_text(workflow)

print(json.dumps({"reactions": len(reactions), "molecules": len(molecules), "pages": [1, 2, 3, 4, 5]}, ensure_ascii=False))
