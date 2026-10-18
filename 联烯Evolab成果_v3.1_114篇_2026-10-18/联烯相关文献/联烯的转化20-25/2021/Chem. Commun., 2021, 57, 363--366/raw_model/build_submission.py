#!/usr/bin/env python3
"""Build the source-specific v3.1 bundle for doc_1fb7cef55cdd0f91."""

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills" / "chemical-reaction-extraction" / "scripts"
sys.path.insert(0, str(SCRIPTS))

import standard_records as st

DOC = "doc_1fb7cef55cdd0f91"
RUN = ROOT / "runs" / "current"
E = {p: [f"{DOC}:page{p}:text", f"{DOC}:page{p}:image"] for p in range(1, 5)}


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def F(raw=None, normalized=None, ev=(), method="original_extraction", status=None,
      derivation=None):
    return st.fact(raw, normalized, ev, method=method, status=status,
                   derivation=derivation)


def qraw(value):
    if isinstance(value, dict) and "number" in value and "unit" in value:
        return value.get("raw_text") or f"{value['number']} {value['unit']}"
    return value


def set_fact(container, key, raw, normalized, ev, method="original_extraction",
             status=None, derivation=None):
    container[key] = F(raw, normalized, ev, method, status, derivation)


def generic_molecule(mid, label, name, fixed, variables, ev):
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["representation_type"] = "generic_template"
    m["occurrences"][0]["compound_label"] = F(label, label, ev)
    m["names"]["original"] = F(name, name, ev)
    m["names"]["english"] = F(name, name, ev)
    m["cas_numbers"] = F()
    var_objs = []
    for variable in variables:
        var_objs.append({
            "label": variable,
            "attachment_points": F(f"Source drawing attachment for {variable}",
                                   f"source-defined attachment for {variable}", ev),
            "definition": F(f"{variable} as drawn in the generic scheme", None, ev,
                            status="unreadable"),
            "allowed_values": F(),
            "constraints": F(),
        })
    m["generic_template"] = {
        "fixed_structure": F(fixed, {"description": fixed}, ev),
        "variables": var_objs,
        "instance_mapping": [],
        "reason": F("Generic/R-group source drawing", "The source drawing contains variable substituents and is not a unique molecule.", ev),
    }
    set_fact(m["stereochemistry"], "chirality_status", "generic structure; configurations vary",
             "unresolved", ev)
    set_fact(m["stereochemistry"], "geometry_status", "generic structure; geometry varies",
             "unresolved", ev)
    set_fact(m["stereochemistry"], "representation_status", "generic template",
             "not_encoded", ev)
    m["identity_note"] = F("Generic source structure", "Variable structure retained without wildcard SMILES.", ev)
    return m


def molecule(mid, label, name, ev, smiles=None, chirality="unspecified",
             geometry="unresolved", representation=None, assignment=None,
             mixture=None, double_bonds=None, axial=None, note=None):
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, ev)
    m["names"]["original"] = F(name, name, ev)
    m["names"]["english"] = F(name, name, ev)
    m["cas_numbers"] = F()
    if smiles:
        diag = st.stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = F(
            f"Candidate graph transcribed/assigned for {label}", smiles, ev,
            method="structure_derivation",
            derivation="Candidate connectivity derived from the source drawing/name; RDKit parsing does not verify identity.",
        )
        m["stereochemistry"]["isomeric_smiles"] = F(
            None, diag["canonical_isomeric_smiles"], ev, method="structure_derivation",
            derivation="RDKit canonical isomeric SMILES for the candidate graph; source identity remains unverified.")
        m["stereochemistry"]["identity_key"] = F(
            None, diag["stereo_graph_key"], ev, method="structure_derivation",
            derivation="RDKit stereo graph key for the candidate graph only.")
        representation = representation or "partial"
    else:
        m["representation_type"] = "unresolved"
        representation = representation or "not_encoded"
    set_fact(m["stereochemistry"], "chirality_status", chirality, chirality, ev)
    set_fact(m["stereochemistry"], "geometry_status", geometry, geometry, ev)
    set_fact(m["stereochemistry"], "representation_status", representation,
             representation, ev)
    if assignment is not None:
        m["stereochemistry"]["assignment"] = F(assignment, assignment, ev)
    if mixture is not None:
        m["stereochemistry"]["mixture"] = F(mixture, mixture, ev)
    if double_bonds is not None:
        m["stereochemistry"]["double_bonds"] = F(double_bonds, double_bonds, ev)
    if axial is not None:
        m["stereochemistry"]["axial"] = F(axial, axial, ev)
    m["identity_note"] = F(note, note, ev) if note else F()
    return m


def label_for(molecules, mid, label, ev):
    m = molecules[mid]
    if not any(o["document_id"] == DOC and o["compound_label"]["normalized_value"] == label
               for o in m["occurrences"]):
        m["occurrences"].append({"document_id": DOC,
                                 "compound_label": F(label, label, ev),
                                 "evidence_ids": list(ev)})


def participant(pid, mid, role, ev, stage="reaction", amount=None, mass=None,
                volume=None, equiv=None, mol_percent=None, reference=None,
                product_role=None):
    p = st.new_participant(pid, mid, role, ev, stage)
    if product_role is not None:
        p["product_role"] = F(product_role, product_role, ev)
    values = {"amount": amount, "mass": mass, "volume": volume,
              "equivalents": equiv, "mol_percent": mol_percent}
    raw_names = {"amount": "amount", "mass": "mass", "volume": "volume",
                 "equivalents": "equivalents", "mol_percent": "mol%"}
    for key, value in values.items():
        if value is not None:
            p["charges"][key] = F(qraw(value), value, ev)
    if reference is not None:
        p["equivalence_reference"] = F(reference, reference, ev)
    return p


def step(step_id, ev, participants, time=None, temperature=None, operation=None,
         atmosphere=None, stage="reaction"):
    s = st.new_step(step_id, ev, stage)
    s["participant_ids"] = list(participants)
    if operation is not None:
        s["operation"] = F(operation, operation, ev)
    if time is not None:
        s["time"] = F(qraw(time), time, ev)
    if temperature is not None:
        s["temperature"] = F(qraw(temperature), temperature, ev)
    if atmosphere is not None:
        s["atmosphere"] = F(atmosphere, atmosphere, ev)
    return s


def outcome(oid, kind, targets, value, ev, role=None, yield_type=None,
            basis=None, scope=None):
    o = st.new_outcome(oid, kind, targets, ev)
    raw = qraw(value)
    o["value"] = F(raw, value, ev)
    if role is not None:
        o["product_role"] = F(role, role, ev)
    if yield_type is not None:
        o["yield_type"] = F(yield_type, yield_type, ev)
    if basis is not None:
        o["measurement_basis"] = F(basis, basis, ev)
    if scope is not None:
        o["scope"] = F(scope, scope, ev)
    return o


def classification(labels, ev, extra_changes=()):
    changes = []
    label_rows = []
    element_pairs = {
        "C-C": ("C", "C", 1), "C-N": ("C", "N", 1),
        "C-S": ("C", "S", 1), "O-H": ("O", "H", 1),
    }
    for code in labels:
        a, b, order = element_pairs[code]
        cid = f"b{len(changes)+1}"
        changes.append({"change_id": cid,
                        "atoms": [{"atom_map": None, "element": a},
                                  {"atom_map": None, "element": b}],
                        "before_order": 0, "after_order": order,
                        "change_type": "hydrogen_gain" if b == "H" else "formed",
                        **({"hydrogen_count": 1} if b == "H" else {})})
        lab = st.bond_label(a, b, order)
        label_rows.append({**lab, "supporting_change_ids": [cid]})
    for ch in extra_changes:
        changes.append(ch)
    derivation = "Partial class-level bond changes inferred from the source transformation drawing/name; exact atoms and complete mapping are unresolved."
    return {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": F(),
        "bond_changes": [F(None, ch, ev, "inference", derivation=derivation) for ch in changes],
        "labels": F(None, label_rows, ev, "inference", derivation=derivation),
        "coverage": F(None, "partial" if labels else "unresolved", ev, "inference",
                      derivation=derivation),
        "limitations": ["Exact atom identities and complete atom mapping are unresolved.",
                        "Hydrogen/proton source and coproducts are not established."],
        "model_votes": [],
    }


def reaction(rid, ev, record_type, categories, name, transformation, labels=(),
             location=None, general=(), unresolved=()):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = F(categories, categories, ev)
    r["location"] = F(location or rid, location or rid, ev)
    r["general_scheme_ids"] = list(general)
    r["reaction_name"] = F(name, name, ev)
    r["transformation_type"] = F(None, transformation, ev, "inference",
                                  derivation="Normalized transformation class inferred from the source wording/drawing.")
    r["classification"] = classification(labels, ev)
    r["review"]["unresolved"] = list(unresolved) or [
        "Supporting Information was not supplied; detailed procedure and full structure verification remain unavailable."]
    return r


def set_scale(r, raw, reference_mid, amount, ev):
    r["reaction_scale"] = F(raw, {"reference_molecule_id": reference_mid,
                                  "amount": amount}, ev)


def set_solvent(r, mid, volume, ev):
    r["solvent_conditions"]["composition"] = F(mid, [mid], ev)
    if volume is not None:
        r["solvent_conditions"]["total_volume"] = F(qraw(volume), volume, ev)


def add_note(r, text, ev, method="original_extraction", derivation=None):
    r["notes"].append(F(text, text, ev, method, derivation=derivation))


dataset = json.loads((RUN / "dataset.json").read_text())
bundle = st.empty(dataset["documents"])
mols = {}

# Generic precedent/current-scheme structures.
mols["lit_lu_reactants"] = generic_molecule("lit_lu_reactants", "Scheme 1a reactants",
    "imine and 2,3-allenoate generic reactants", "Imine plus allenoate as drawn", ["Ar", "EWG"], E[1])
mols["lit_lu_product"] = generic_molecule("lit_lu_product", "Scheme 1a product",
    "chiral pyrroline generic product", "Pyrroline bearing EWG/Ar substituents", ["Ar", "EWG"], E[1])
mols["lit_gong_reactants"] = generic_molecule("lit_gong_reactants", "Scheme 1b reactants",
    "2,3-allenoate and azomethine ylide generic reactants", "Allenoate plus azomethine ylide", ["R", "Ar"], E[1])
mols["lit_gong_product"] = generic_molecule("lit_gong_product", "Scheme 1b product",
    "chiral pyrrolidine generic product", "Substituted pyrrolidine product", ["R", "Ar"], E[1])
mols["generic_1"] = generic_molecule("generic_1", "1", "4-isothiocyanato pyrazolone family 1",
    "4-isothiocyanato pyrazolone core", ["R", "R1"], E[1] + E[2])
mols["generic_2"] = generic_molecule("generic_2", "2", "allenyl ketone family 2",
    "Allenyl ketone core", ["R2"], E[1] + E[2])
mols["generic_4"] = generic_molecule("generic_4", "4", "spiro[pyrrole-pyrazolone] family 4",
    "Spiro[pyrrole-pyrazolone] bearing ketone side chains", ["R", "R1", "R2"], E[1] + E[2])

# Concrete reactants and common agents.
mols["1a"] = molecule("1a", "1a", "4-isothiocyanato-1,3-diphenylpyrazolone 1a", E[1] + E[2],
    smiles="O=C1N(c2ccccc2)N=C(c2ccccc2)C1N=C=S", chirality="unspecified",
    geometry="unresolved", representation="partial",
    assignment="Candidate graph contains an unassigned tetrahedral C4 center; source gives no configuration.",
    note="Candidate connectivity only; tautomerism and exact source drawing require expert review.")
for label in list("bcdefghijk"):
    mid = f"1{label}"
    mols[mid] = molecule(mid, mid, f"4-isothiocyanato pyrazolone {mid}", E[2],
                         chirality="unspecified", geometry="unresolved", representation="not_encoded",
                         note="Compact source drawing retained without an asserted complete SMILES.")
mols["2a"] = molecule("2a", "2a", "phenyl allenone 2a", E[1] + E[2] + E[3],
    smiles="C=C=CC(=O)c1ccccc1", chirality="achiral", geometry="not_applicable",
    representation="partial",
    axial="The terminal CH2 end gives two identical substituents, so the depicted allene axis is not stereogenic on the candidate graph.",
    assignment="Achiral status is a structure-derived statement for the candidate graph, not a source stereochemical assignment.")
for label in list("bcdefghij"):
    mid = f"2{label}"
    mols[mid] = molecule(mid, mid, f"allenyl ketone {mid}", E[2] + E[3],
                         chirality="unspecified", geometry="unresolved", representation="not_encoded",
                         note="Source label/drawing retained; complete substituent graph was not asserted.")
mols["3aa"] = molecule("3aa", "3aa", "spiro[pyrrole-pyrazolone] intermediate 3aa", E[2],
                        chirality="unspecified", geometry="unspecified", representation="not_encoded")

simple = {
    "Q1": ("Q1", "quinine", None), "Q2": ("Q2", "quinidine", None),
    "Q3": ("Q3", "quinine-derived thiourea catalyst Q3", None),
    "Q4": ("Q4", "quinine-derived squaramide catalyst Q4", None),
    "DABCO": ("DABCO", "1,4-diazabicyclo[2.2.2]octane", "C1CN2CCN1CC2"),
    "Et3N": ("Et3N", "triethylamine", "CCN(CC)CC"),
    "DCM": ("DCM", "dichloromethane", "ClCCl"),
    "CHCl3": ("CHCl3", "chloroform", "ClC(Cl)Cl"),
    "THF": ("THF", "tetrahydrofuran", "C1CCOC1"),
    "toluene": ("toluene", "toluene", "Cc1ccccc1"),
    "MeCN": ("CH3CN", "acetonitrile", "CC#N"),
    "DCE": ("DCE", "1,2-dichloroethane", "ClCCCl"),
}
for mid, (label, name, smi) in simple.items():
    mols[mid] = molecule(mid, label, name, E[2] + E[3] + E[4], smiles=smi,
                         chirality="achiral" if smi and mid not in {"Q1", "Q2", "Q3", "Q4"} else "specified" if mid in {"Q1", "Q2", "Q3", "Q4"} else "unresolved",
                         geometry="not_applicable" if smi else "unresolved",
                         representation="encoded" if smi else "not_encoded",
                         assignment=("Source-drawn chiral catalyst identity; complete stereochemical graph not encoded."
                                     if mid in {"Q1", "Q2", "Q3", "Q4"} else
                                     "Candidate graph has no potential stereogenic elements in standard graph inspection." if smi else None))
mols["DABCO"]["stereochemistry"]["chirality_status"] = F(
    "Toolkit flags an unspecified potential element in the cage graph; no source stereochemical claim.",
    "unresolved", E[3])
mols["DABCO"]["stereochemistry"]["representation_status"] = F(
    "Candidate connectivity encoded but stereochemical toolkit diagnostic unresolved.",
    "partial", E[3])

# Product-4 identities. Q2 and catalyst-free samples are kept distinct from the Q4-family sample.
scope4 = ["4aa", "4ba", "4ca", "4da", "4ea", "4fa", "4ga", "4ha", "4ia", "4ja", "4ka",
          "4ab", "4ac", "4ad", "4ae", "4af", "4ag", "4ah", "4ai", "4aj"]
for mid in scope4:
    is_af = mid == "4af"
    mols[mid] = molecule(mid, mid, f"spiro[pyrrole-pyrazolone] {mid}", E[2] + (E[3] if mid[1] == "a" and mid != "4aa" else []),
        chirality="unspecified", geometry="specified" if is_af else "unspecified",
        representation="not_encoded",
        double_bonds=({"side_chain_alkene": "Z", "assignment_method": "NOE study; details in unavailable ESI"} if is_af else None),
        assignment=("Side-chain olefin geometry Z was assigned by NOE; absolute spiro configuration is not stated on the supplied pages." if is_af else None),
        note="Experiment-specific ee is stored on reaction outcomes; missing absolute configuration is not converted to a racemate.")
mols["4aa_antipode"] = molecule("4aa_antipode", "4aa (antipode favored)",
    "4aa sample favoring the antipode", E[2], chirality="unspecified", geometry="unspecified",
    representation="not_encoded", mixture="Entry 2 reports 11% ee favoring the antipode of the Q1/Q4 product.")
mols["rac_4aa"] = molecule("rac_4aa", "rac-4aa", "racemic 4aa", E[2] + E[4],
    chirality="racemic", geometry="unspecified", representation="not_encoded",
    mixture="Racemic sample; Table 1 entry 12 reports 0% ee.")
for mid in ["rac_4ba", "rac_4ha", "rac_4ia", "rac_4ja", "rac_4af"]:
    label = mid.removeprefix("rac_")
    mols[mid] = molecule(mid, f"(±)-{label}", f"racemic {label}", E[4],
        chirality="racemic", geometry="specified" if label == "4af" else "unspecified",
        representation="not_encoded", mixture="Racemic starting material in Scheme 5.",
        double_bonds=({"side_chain_alkene": "Z", "assignment_basis": "4af NOE statement"} if label == "4af" else None))

# Product-5 identities, including the distinct diastereomer/enantiomer series.
five_labels = ["5aa", "5prime_aa", "5ba", "5prime_ba", "5da", "5ea", "5ga",
               "5ha", "5prime_ha", "5ia", "5prime_ia", "5ja", "5prime_ja",
               "5ka", "5ah", "5af", "5prime_af", "5doubleprime_af"]
display5 = {
    "5prime_aa": "5′aa", "5prime_ba": "5′ba", "5prime_ha": "5′ha",
    "5prime_ia": "5′ia", "5prime_ja": "5′ja", "5prime_af": "5′af",
    "5doubleprime_af": "5′′af",
}
for mid in five_labels:
    label = display5.get(mid, mid)
    xray = mid in {"5af", "5doubleprime_af"}
    assignment = None
    if mid == "5af":
        assignment = "Absolute configuration determined by X-ray crystallographic analysis (CCDC 2034086); R/S descriptors are not printed on the supplied pages."
    elif mid == "5doubleprime_af":
        assignment = "Absolute configuration determined by X-ray crystallographic analysis (CCDC 2040382); this is the enantiomer of 5′af."
    mols[mid] = molecule(mid, label, f"spiro[thiopyranopyrrole-pyrazolone] {label}", E[3] + (E[4] if "prime" in mid or mid in {"5aa", "5ba", "5ha", "5ia", "5ja", "5af"} else []),
        chirality="partially_specified" if xray else "unspecified", geometry="unresolved",
        representation="not_encoded", assignment=assignment,
        note="Compact polycyclic structure retained without a complete stereochemical SMILES; ee/dr are experiment-specific outcomes.")
for product in ["5aa", "5da", "5ea", "5ga", "5ha", "5ja", "5ka", "5ah"]:
    mid = f"mix_{product}"
    mols[mid] = molecule(mid, f"{product} total diastereomer mixture",
        f"total isolated diastereomeric mixture reported for {product}", E[3],
        chirality="mixture", geometry="unresolved", representation="not_encoded",
        mixture="Scheme 4 caption reports isolated yield of the total; the experiment-specific dr is stored on the reaction.",
        note="Mixture identity separates total isolated yield from the major product whose ee is reported.")

reactions = []

# Scheme 1 precedents and current generic scheme.
r = reaction("scheme1a_lu_precedent", E[1], "literature_example", [],
             "Lu phosphine-catalyzed [3+2] cyclization", "[3+2] annulation", ["C-C", "C-N"], "Scheme 1a")
r["participants"] = [participant("lit1_reactants", "lit_lu_reactants", "reactant", E[1]),
                     participant("lit1_product", "lit_lu_product", "product", E[1], product_role="reported precedent product")]
r["steps"] = [step("literature_arrow", E[1], ["lit1_reactants", "lit1_product"], operation="Chiral phosphine-catalyzed cyclization")]
reactions.append(r)
r = reaction("scheme1b_gong_precedent", E[1], "literature_example", [],
             "chiral phosphoric acid-catalyzed [3+2] cycloaddition", "[3+2] cycloaddition", ["C-C"], "Scheme 1b")
r["participants"] = [participant("lit2_reactants", "lit_gong_reactants", "reactant", E[1]),
                     participant("lit2_product", "lit_gong_product", "product", E[1], product_role="reported precedent product")]
r["steps"] = [step("literature_arrow", E[1], ["lit2_reactants", "lit2_product"], operation="Bisphosphoric acid-catalyzed cycloaddition")]
reactions.append(r)
GEN = f"{DOC}:scheme1c_current_generic"
r = reaction("scheme1c_current_generic", E[1] + E[2], "generic_scheme", [],
             "Michael/cyclization sequence for chiral pyrrole derivatives", "asymmetric annulation", ["C-C", "C-N"], "Scheme 1c / Scheme 2 generic header")
r["participants"] = [participant("g1", "generic_1", "reactant", E[1] + E[2]),
                     participant("g2", "generic_2", "reactant", E[1] + E[2]),
                     participant("g4", "generic_4", "product", E[1] + E[2], product_role="generic product")]
r["steps"] = [step("generic_arrow", E[1], ["g1", "g2", "g4"], operation="Chiral-catalyst Michael/cyclization sequence")]
reactions.append(r)

# Table 1 optimization.
table1 = [
    (1, "Q1", "DCM", 10, 0.5, "trace", 59, 20, 31, "4aa"),
    (2, "Q2", "DCM", 10, 0.5, "trace", 57, 19, 11, "4aa_antipode"),
    (3, "Q3", "DCM", 10, 0.5, "trace", 59, 18, 52, "4aa"),
    (4, "Q4", "DCM", 10, 0.5, "trace", 59, 15, 91, "4aa"),
    (5, "Q4", "CHCl3", 10, 0.5, "trace", 60, 19, 96, "4aa"),
    (6, "Q4", "THF", 10, 0.5, "trace", 66, 18, 94, "4aa"),
    (7, "Q4", "toluene", 10, 0.5, "trace", 66, 18, 99, "4aa"),
    (8, "Q4", "MeCN", 10, 0.5, "trace", 57, 20, 54, "4aa"),
    (9, "Q4", "DCE", 10, 0.5, "trace", 60, 15, 93, "4aa"),
    (10, "Q4", "toluene", 5, 0.5, "trace", 66, 15, 97, "4aa"),
    (11, "Q4", "toluene", 2.5, 0.5, "trace", 70, 13, 82, "4aa"),
    (12, None, "toluene", None, 6, "trace", 80, "trace", 0, "rac_4aa"),
]
for entry, cat, solv, loading, hours, y3, y4, y5, ee, product4 in table1:
    ev = E[2]
    r = reaction(f"table1_entry{entry}", ev, "specific_experiment", ["optimization"],
                 "cascade spiroannulation/sulfa-Michael addition optimization",
                 "sequential annulation/sulfa-Michael addition", ["C-C", "C-N", "C-S"],
                 f"Table 1 entry {entry}", [GEN])
    ps = [participant("p1a", "1a", "reactant", ev, equiv={"number": 1, "unit": "equiv"}, reference="1a"),
          participant("p2a", "2a", "reactant", ev, equiv={"number": 2.5, "unit": "equiv"}, reference="1a"),
          participant("p3aa", "3aa", "product", ev, product_role="trace annulation intermediate"),
          participant("p4aa", product4, "product", ev, product_role="major cascade product"),
          participant("p5aa", "5aa", "product", ev, product_role="aldol cyclization product"),
          participant("psolv", solv, "solvent", ev, volume={"number": 1, "unit": "mL"})]
    if cat:
        ps.append(participant("pcat", cat, "catalyst", ev,
                              mol_percent={"number": loading, "unit": "mol%"}, reference="1a"))
    r["participants"] = ps
    set_scale(r, "0.1 mmol scale with 1a", "1a", {"number": 0.1, "unit": "mmol"}, ev)
    set_solvent(r, solv, {"number": 1, "unit": "mL"}, ev)
    r["steps"] = [step("reaction", ev, [p["participant_id"] for p in ps if p["stage"] == "reaction"],
                       time={"number": hours, "unit": "h"}, temperature="r.t.",
                       operation="Reaction under the Table 1 conditions")]
    r["procedure"]["original_text"] = F("Table 1 footnote: 1a (1.0 equiv.), 2a (2.5 equiv.), catalyst in solvent; isolated yields.",
                                           "Table 1 experiment", ev)
    r["outcomes"] = [
        outcome("trace_3aa", "trace", ["3aa"], y3, ev, "trace intermediate", basis="isolated-yield table column"),
        outcome("yield_4aa", "yield", [product4], {"number": y4, "unit": "%"}, ev,
                "cascade product", "isolated yield", "isolated"),
        outcome("yield_5aa" if y5 != "trace" else "trace_5aa", "yield" if y5 != "trace" else "trace", ["5aa"],
                {"number": y5, "unit": "%"} if y5 != "trace" else "trace", ev,
                "aldol cyclization product", "isolated yield" if y5 != "trace" else None, "isolated"),
        outcome("ee_4aa", "selectivity", [product4], {"number": ee, "unit": "%"}, ev,
                basis="chiral HPLC", scope={"metric": "ee", "major_enantiomer": "antipode" if entry == 2 else "Q1/Q4-series or racemate"}),
    ]
    if entry == 2:
        r["outcomes"][-1]["value"]["raw_value"] = "−11%"
        add_note(r, "The table prints −11% ee; prose states that the antipode of 4aa is favored. The normalized magnitude is 11% with antipode direction retained in scope.", ev)
    reactions.append(r)

# Scheme 2 substrate scope, including the distinct 0.2 mmol 4aa experiment.
scope2 = [
    ("4aa", "1a", "2a", 66, 99), ("4ba", "1b", "2a", 73, 95),
    ("4ca", "1c", "2a", 73, 98), ("4da", "1d", "2a", 73, 99),
    ("4ea", "1e", "2a", 70, 79), ("4fa", "1f", "2a", 67, 99),
    ("4ga", "1g", "2a", 69, 77), ("4ha", "1h", "2a", 40, 97),
    ("4ia", "1i", "2a", 81, 92), ("4ja", "1j", "2a", 70, 89),
    ("4ka", "1k", "2a", 60, 95), ("4ab", "1a", "2b", 67, 98),
    ("4ac", "1a", "2c", 60, 96), ("4ad", "1a", "2d", 50, 92),
    ("4ae", "1a", "2e", 62, 92), ("4af", "1a", "2f", 71, 99),
    ("4ag", "1a", "2g", 76, 92), ("4ah", "1a", "2h", 70, 89),
    ("4ai", "1a", "2i", 73, 95), ("4aj", "1a", "2j", 57, 75),
]
scheme2_ids = []
for product, m1, m2, yld, ee in scope2:
    ev = E[2] + (E[3] if product in {f"4a{x}" for x in "bcdefghij"} else [])
    rid = f"scheme2_{product}"
    scheme2_ids.append(f"{DOC}:{rid}")
    r = reaction(rid, ev, "specific_experiment", ["substrate_scope"],
                 "asymmetric sequential annulation/sulfa-Michael addition",
                 "sequential annulation/sulfa-Michael addition", ["C-C", "C-N", "C-S"],
                 f"Scheme 2 product {product}", [GEN])
    r["participants"] = [
        participant("p1", m1, "reactant", ev, equiv={"number": 1, "unit": "equiv"}, reference=m1),
        participant("p2", m2, "reactant", ev, equiv={"number": 2.5, "unit": "equiv"}, reference=m1),
        participant("pQ4", "Q4", "catalyst", ev, mol_percent={"number": 10, "unit": "mol%"}, reference=m1),
        participant("ptol", "toluene", "solvent", ev, volume={"number": 2, "unit": "mL"}),
        participant("pprod", product, "product", ev, product_role="isolated spiro[pyrrole-pyrazolone]"),
    ]
    set_scale(r, "0.2 mmol scale with substrate 1", m1, {"number": 0.2, "unit": "mmol"}, ev)
    set_solvent(r, "toluene", {"number": 2, "unit": "mL"}, ev)
    r["steps"] = [step("reaction", ev, [p["participant_id"] for p in r["participants"]],
                       time={"number": 0.5, "unit": "h"}, temperature="r.t.",
                       operation="Q4-catalyzed Scheme 2 reaction")]
    r["outcomes"] = [
        outcome("yield", "yield", [product], {"number": yld, "unit": "%"}, ev,
                "desired product", "isolated yield", "isolated"),
        outcome("ee", "selectivity", [product], {"number": ee, "unit": "%"}, ev,
                basis="chiral HPLC", scope={"metric": "ee", "absolute_configuration": "unassigned"}),
    ]
    if product == "4af":
        add_note(r, "The side-chain olefin geometry was determined to be Z by NOE study; detailed NOE data are in unavailable ESI.", E[3])
    if product == "4aa":
        add_note(r, "This Scheme 2 experiment is not merged with Table 1 entry 7 because its reported scale/solvent volume are 0.2 mmol/2.0 mL rather than 0.1 mmol/1.0 mL.", ev,
                 method="inference", derivation="Comparison of the Table 1 and Scheme 2 captions.")
    reactions.append(r)

# Gram-scale and prolonged one-pot experiment.
r = reaction("scheme3_gram_scale_4aa", E[3], "specific_experiment", ["scale_up"],
             "gram-scale cascade spiroannulation/sulfa-Michael addition",
             "sequential annulation/sulfa-Michael addition", ["C-C", "C-N", "C-S"],
             "Scheme 3", [GEN], ["Solvent volume, atmosphere, workup, purification, and yield/ee measurement methods are unreported on pages 1–4."])
r["participants"] = [participant("p1a", "1a", "reactant", E[3], amount={"number": 2.5, "unit": "mmol"}),
                     participant("p2a", "2a", "reactant", E[3], amount={"number": 6.3, "unit": "mmol"}),
                     participant("pQ4", "Q4", "catalyst", E[3], amount={"number": 0.25, "unit": "mmol"}),
                     participant("ptol", "toluene", "solvent", E[3]),
                     participant("p4aa", "4aa", "product", E[3], mass={"number": 0.97, "unit": "g"}, product_role="gram-scale product")]
set_scale(r, "1a: 2.5 mmol", "1a", {"number": 2.5, "unit": "mmol"}, E[3])
set_solvent(r, "toluene", None, E[3])
r["steps"] = [step("reaction", E[3], [p["participant_id"] for p in r["participants"]],
                   time={"number": 0.5, "unit": "h"}, temperature="r.t.", operation="Gram-scale reaction")]
r["outcomes"] = [outcome("product_amount", "product_amount", ["4aa"], {"number": 0.97, "unit": "g"}, E[3], "product"),
                 outcome("yield", "yield", ["4aa"], {"number": 66, "unit": "%"}, E[3], "product"),
                 outcome("ee", "selectivity", ["4aa"], {"number": 97, "unit": "%"}, E[3], scope={"metric": "ee", "absolute_configuration": "unassigned"})]
reactions.append(r)

r = reaction("prolonged_onepot_5aa", E[2] + E[3], "specific_experiment", ["other"],
             "prolonged sequential annulation/aldol process to 5aa",
             "sequential annulation/sulfa-Michael/vinylogous aldol", ["C-C", "C-N", "C-S", "O-H"],
             "Page 3 prose", [GEN], ["The prose gives only 36 h and 68% yield; exact scale, doses, workup, purification, and measurement basis are unreported."])
r["participants"] = [participant("p1a", "1a", "reactant", E[2] + E[3]),
                     participant("p2a", "2a", "reactant", E[2] + E[3]),
                     participant("p5aa", "5aa", "product", E[3], product_role="prolonged-time polyheterocycle")]
r["steps"] = [step("reaction", E[3], ["p1a", "p2a", "p5aa"], time={"number": 36, "unit": "h"},
                   operation="Prolonged reaction")]
r["outcomes"] = [outcome("yield", "yield", ["5aa"], {"number": 68, "unit": "%"}, E[3], "product")]
add_note(r, "The authors presume 5aa derives from 4aa through an intramolecular vinylogous aldol reaction.", E[3])
reactions.append(r)

# Table 2 optimization and Scheme 4 scope. Entry 3 and Scheme 4 5af are one experiment.
table2 = [
    (1, "DABCO", 24, "1:3", 14, 31, 99, 99),
    (2, "Et3N", 24, "1:1.4", 27, 33, 99, 99),
    (3, "Q4", 12, ">20:1", 85, None, 99, None),
]
aldol_ids = []
for entry, cat, hours, dr, y_main, y_minor, ee_main, ee_minor in table2:
    rid = f"table2_entry{entry}" + ("_scheme4_5af" if entry == 3 else "")
    aldol_ids.append(f"{DOC}:{rid}")
    cats = ["optimization"] + (["substrate_scope"] if entry == 3 else [])
    r = reaction(rid, E[3], "specific_experiment", cats,
                 "asymmetric vinylogous aldol cyclization of 4af",
                 "intramolecular vinylogous aldol cyclization", ["C-C", "O-H"],
                 f"Table 2 entry {entry}" + (" / Scheme 4 product 5af" if entry == 3 else ""),
                 unresolved=["Complete stereochemical structures and atom mapping are unresolved; SI procedure was not supplied."])
    r["participants"] = [participant("p4af", "4af", "reactant", E[3]),
                         participant("pcat", cat, "catalyst", E[3], mol_percent={"number": 10, "unit": "mol%"}, reference="4af"),
                         participant("pDCM", "DCM", "solvent", E[3], volume={"number": 1, "unit": "mL"}),
                         participant("p5af", "5af", "product", E[3], product_role="aldol product 5af"),
                         participant("p5paf", "5prime_af", "product", E[3], product_role="diastereomer 5′af")]
    set_scale(r, "0.1 mmol scale with 4af", "4af", {"number": 0.1, "unit": "mmol"}, E[3])
    set_solvent(r, "DCM", {"number": 1, "unit": "mL"}, E[3])
    r["steps"] = [step("reaction", E[3], [p["participant_id"] for p in r["participants"]],
                       time={"number": hours, "unit": "h"}, temperature="r.t.", operation="Catalyzed aldol cyclization")]
    r["outcomes"] = [outcome("yield_5af", "yield", ["5af"], {"number": y_main, "unit": "%"}, E[3],
                                     "5af", "isolated yield", "isolated"),
                     outcome("ee_5af", "selectivity", ["5af"], {"number": ee_main, "unit": "%"}, E[3],
                             basis="chiral HPLC", scope={"metric": "ee", "configuration": "5af"}),
                     outcome("dr", "selectivity", ["5af", "5prime_af"], dr, E[3],
                             basis="1H NMR of crude reaction mixture", scope={"metric": "dr", "order": ["5af", "5′af"]})]
    if y_minor is not None:
        r["outcomes"] += [outcome("yield_5prime_af", "yield", ["5prime_af"], {"number": y_minor, "unit": "%"}, E[3],
                                  "5′af", "isolated yield", "isolated"),
                          outcome("ee_5prime_af", "selectivity", ["5prime_af"], {"number": ee_minor, "unit": "%"}, E[3],
                                  basis="chiral HPLC", scope={"metric": "ee", "configuration": "5′af"})]
    else:
        r["outcomes"].append(outcome("minor_product_table_dash", "other", ["5prime_af"], "table dash (—)", E[3],
                                             scope={"column": "yield/ee of 5′af", "meaning": "not numerically reported"}))
    reactions.append(r)

scheme4 = [
    ("5aa", "4aa", 89, 93, ">20:1"), ("5da", "4da", 90, 91, "7.3:1"),
    ("5ea", "4ea", 83, 99, "4:1"), ("5ga", "4ga", 81, 99, "6:1"),
    ("5ha", "4ha", 88, 95, ">20:1"), ("5ja", "4ja", 81, 99, "4:1"),
    ("5ka", "4ka", 90, 95, ">20:1"), ("5ah", "4ah", 93, 87, "4.5:1"),
]
for product, substrate, yld, ee, dr in scheme4:
    rid = f"scheme4_{product}"
    aldol_ids.append(f"{DOC}:{rid}")
    r = reaction(rid, E[3], "specific_experiment", ["substrate_scope"],
                 "asymmetric vinylogous aldol cyclization",
                 "intramolecular vinylogous aldol cyclization", ["C-C", "O-H"],
                 f"Scheme 4 product {product}", unresolved=["The minor diastereomer identity and complete stereochemical graph are not established on the supplied pages."])
    mix = f"mix_{product}"
    r["participants"] = [participant("psub", substrate, "reactant", E[3]),
                         participant("pQ4", "Q4", "catalyst", E[3], mol_percent={"number": 10, "unit": "mol%"}, reference=substrate),
                         participant("pDCM", "DCM", "solvent", E[3], volume={"number": 1, "unit": "mL"}),
                         participant("pmix", mix, "product", E[3], product_role="total isolated diastereomer mixture"),
                         participant("pprod", product, "product", E[3], product_role="major aldol product used for ee reporting")]
    set_scale(r, f"0.1 mmol scale with {substrate}", substrate, {"number": 0.1, "unit": "mmol"}, E[3])
    set_solvent(r, "DCM", {"number": 1, "unit": "mL"}, E[3])
    r["steps"] = [step("reaction", E[3], [p["participant_id"] for p in r["participants"]],
                       time={"number": 12, "unit": "h"}, temperature="r.t.", operation="Q4-catalyzed aldol cyclization")]
    r["outcomes"] = [outcome("total_isolated_yield", "yield", [mix], {"number": yld, "unit": "%"}, E[3],
                                     "reported product/total", "isolated yield of the total", "isolated yield of the total"),
                     outcome("ee", "selectivity", [product], {"number": ee, "unit": "%"}, E[3],
                             basis="chiral HPLC", scope={"metric": "ee", "configuration": "unassigned except 5af X-ray"}),
                     outcome("dr", "selectivity", [mix], dr, E[3], basis="1H NMR of crude reaction mixture",
                             scope={"metric": "dr", "major": product, "minor_identity": "unreported"})]
    reactions.append(r)

# Scheme 5 parallel resolution; each drawn product pair is one reaction experiment.
scheme5 = [
    ("aa", "rac_4aa", "5aa", "5prime_aa", 43, 99, 44, 89),
    ("ba", "rac_4ba", "5ba", "5prime_ba", 41, 99, 37, 99),
    ("ha", "rac_4ha", "5ha", "5prime_ha", 47, 97, 47, 91),
    ("ia", "rac_4ia", "5ia", "5prime_ia", 43, 99, 48, 87),
    ("ja", "rac_4ja", "5ja", "5prime_ja", 41, 97, 46, 89),
    ("af", "rac_4af", "5af", "5doubleprime_af", 45, 99, 43, 91),
]
scheme5_ids = []
for code, substrate, p1, p2, y1, ee1, y2, ee2 in scheme5:
    rid = f"scheme5_parallel_resolution_{code}"
    scheme5_ids.append(f"{DOC}:{rid}")
    ev = E[4] + (E[3] if code == "af" else [])
    r = reaction(rid, ev, "specific_experiment", ["substrate_scope"],
                 "parallel resolution by asymmetric vinylogous aldol cyclization",
                 "intramolecular vinylogous aldol cyclization/parallel resolution", ["C-C", "O-H"],
                 f"Scheme 5 pair {code}", unresolved=["Yield/ee methods are reported, but full structures, dr, workup, purification, and absolute configurations are not fully provided on pages 1–4."])
    r["participants"] = [participant("psub", substrate, "reactant", ev),
                         participant("pQ4", "Q4", "catalyst", ev, mol_percent={"number": 10, "unit": "mol%"}, reference=substrate),
                         participant("pchl", "CHCl3", "solvent", ev, volume={"number": 1, "unit": "mL"}),
                         participant("pp1", p1, "product", ev, product_role="first resolved product"),
                         participant("pp2", p2, "product", ev, product_role="paired resolved product")]
    set_scale(r, f"0.1 mmol scale with (±)-4{code}", substrate, {"number": 0.1, "unit": "mmol"}, ev)
    set_solvent(r, "CHCl3", {"number": 1, "unit": "mL"}, ev)
    r["steps"] = [step("reaction", ev, [p["participant_id"] for p in r["participants"]],
                       time={"number": 12, "unit": "h"}, temperature="r.t.", operation="Q4-catalyzed parallel resolution")]
    r["outcomes"] = [outcome("yield_first", "yield", [p1], {"number": y1, "unit": "%"}, ev, "first product", "isolated yield", "isolated"),
                     outcome("ee_first", "selectivity", [p1], {"number": ee1, "unit": "%"}, ev, basis="chiral HPLC", scope={"metric": "ee", "product": p1}),
                     outcome("yield_second", "yield", [p2], {"number": y2, "unit": "%"}, ev, "paired product", "isolated yield", "isolated"),
                     outcome("ee_second", "selectivity", [p2], {"number": ee2, "unit": "%"}, ev, basis="chiral HPLC", scope={"metric": "ee", "product": p2})]
    if code == "af":
        add_note(r, "5′′af is reported as the enantiomer of 5′af and its absolute configuration was determined by X-ray crystallographic analysis.", E[3] + E[4])
    reactions.append(r)

bundle["molecules"] = list(mols.values())
bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": f"{DOC}:metadata", "document_id": DOC, "kind": "bibliographic_metadata",
     "evidence_ids": E[1], "title": dataset["documents"][0]["title"], "doi": "10.1039/d0cc07113g"},
    {"context_id": f"{DOC}:esi_ccdc_notice", "document_id": DOC, "kind": "external_availability_notice",
     "evidence_ids": E[1] + E[3], "text": "ESI is advertised; CCDC 2034086 (5af) and 2040382 (5′′af) are cited. No SI/CIF was supplied or associated."},
    {"context_id": f"{DOC}:stereo_assignments", "document_id": DOC, "kind": "structure_assignment_context",
     "evidence_ids": E[3], "text": "4af side-chain Z geometry was assigned by NOE; 5af and 5′′af absolute configurations were assigned by X-ray. Detailed supporting data are unavailable."},
    {"context_id": f"{DOC}:negative_spectral_inventory", "document_id": DOC, "kind": "negative_inventory",
     "evidence_ids": E[1] + E[2] + E[3] + E[4], "text": "No standalone mass-spectrum peak table, molecular-weight table, or fragment-relative-abundance block occurs on pages 1–4."},
]

write(RUN / "submitted.json", bundle)

# Raw extraction/page inventory kept separate from formal records.
page_inventory = {
    "document_id": DOC,
    "pages": [
        {"page": 1, "inventory": ["title/abstract/metadata", "Scheme 1a literature precedent", "Scheme 1b literature precedent", "Scheme 1c current generic scheme", "ESI/CCDC availability notice"],
         "decisions": ["No concrete current experiment is reported on page 1.", "External SI/CIF files are not associated."]},
        {"page": 2, "inventory": ["Table 1 entries 1–12", "Scheme 2 products 4aa–4ka and 4ab–4aj", "prose product distribution and catalyst/solvent discussion"],
         "decisions": ["All 12 table rows are separate.", "All 20 Scheme 2 outcomes are separate; Scheme 2 4aa is not Table 1 entry 7 because scale/volume differ."]},
        {"page": 3, "inventory": ["Scheme 2 allenone-scope discussion", "Scheme 3 scale-up", "36 h/68% 5aa prose experiment", "Table 2 entries 1–3", "Scheme 4 products 5aa/5da/5ea/5ga/5ha/5ja/5ka/5ah/5af", "parallel-resolution prose"],
         "decisions": ["Table 2 entry 3 and Scheme 4 5af are merged.", "Scheme 5 5af/5′′af prose and page-4 drawing are one experiment."]},
        {"page": 4, "inventory": ["Scheme 5 six parallel-resolution arrows", "conclusion", "references"],
         "decisions": ["Six product-pair arrows are six experiments.", "References are not imported as current-paper experiments."]},
    ],
    "limitations": ["Compact structures were retained without unsafe complete SMILES.", "SI procedures and characterization were not supplied."],
}
write(RUN / "raw_model" / "page_extraction.json", page_inventory)

page2_ids = [f"{DOC}:table1_entry{i}" for i in range(1, 13)] + scheme2_ids
page3_ids = ([f"{DOC}:scheme2_4a{x}" for x in "bcdefghij"] +
             [f"{DOC}:scheme3_gram_scale_4aa", f"{DOC}:prolonged_onepot_5aa"] + aldol_ids +
             [f"{DOC}:scheme5_parallel_resolution_af"])
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title, abstract, author/DOI and ESI notice", "finding": "Bibliographic context; ESI/CCDC are advertised but no associated files were supplied.", "evidence_ids": E[1]},
            {"region": "Scheme 1a–b", "finding": "Two cited prior-work reaction families, retained as literature examples.", "evidence_ids": E[1]},
            {"region": "Scheme 1c", "finding": "Current generic Michael/cyclization strategy; not a concrete experiment.", "evidence_ids": E[1]},
        ], "record_ids": [f"{DOC}:scheme1a_lu_precedent", f"{DOC}:scheme1b_gong_precedent", GEN],
           "unresolved": ["Detailed ESI procedures and external CCDC structures are unavailable."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "left-column prose above Table 1", "finding": "Current reaction/product-distribution and optimization discussion linked to Table 1 records.", "evidence_ids": E[2]},
            {"region": "Table 1 including reaction drawings, catalyst structures, entries 1–12 and footnotes", "finding": "Twelve separate optimization experiments with isolated-yield columns and chiral-HPLC ee.", "evidence_ids": E[2]},
            {"region": "Scheme 2 including header, twenty product drawings, caption and right-column discussion", "finding": "Twenty distinct 0.2 mmol substrate-scope experiments; 4aa remains distinct from Table 1 entry 7 because scale and solvent volume differ.", "evidence_ids": E[2]},
        ], "record_ids": page2_ids,
           "unresolved": ["Complete compact-product SMILES and absolute configurations are not established; no SI procedure is available."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "upper-left Scheme 2 continuation prose", "finding": "Allenyl-ketone scope discussion; 4af side-chain Z geometry is assigned by NOE, with details deferred to unavailable ESI.", "evidence_ids": E[3]},
            {"region": "Scheme 3 and scale-up prose", "finding": "One gram-scale 4aa experiment (2.5/6.3/0.25 mmol; 0.97 g, 66%, 97% ee).", "evidence_ids": E[3]},
            {"region": "middle-left 5aa prose", "finding": "One 36 h prolonged one-pot experiment gives 5aa in 68%; exact doses and measurement basis are not reported.", "evidence_ids": E[3]},
            {"region": "Table 2 entries 1–3", "finding": "Three aldol-cyclization optimization experiments; entry 3 duplicates Scheme 4 product 5af and is merged.", "evidence_ids": E[3]},
            {"region": "Scheme 4 and caption", "finding": "Nine scope outcomes; eight new records plus the merged Table 2 entry 3/5af record.", "evidence_ids": E[3]},
            {"region": "bottom-right parallel-resolution prose", "finding": "5af/5′′af values and X-ray assignment linked to Scheme 5's matching arrow.", "evidence_ids": E[3]},
        ], "record_ids": page3_ids,
           "unresolved": ["Scheme-level procedures/workups and most absolute configurations are unavailable; full structures remain unencoded where unsafe."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "Scheme 5 and caption", "finding": "Six distinct parallel-resolution experiments, each yielding a pair of enantioenriched products under common conditions.", "evidence_ids": E[4]},
            {"region": "conclusion", "finding": "Summary of current annulation/aldol chemistry; no additional concrete experiment.", "evidence_ids": E[4]},
            {"region": "references", "finding": "Bibliography only; no current-paper synthesis or standalone spectral table.", "evidence_ids": E[4]},
        ], "record_ids": scheme5_ids,
           "unresolved": ["Scheme 5 does not report dr, workup/purification, or complete absolute configurations on page 4."]},
    ],
}
write(RUN / "coverage" / "pages_0001_0004.json", coverage)

stereo_checks = {}
for mid, m in mols.items():
    smi = m["smiles"]["normalized_value"]
    if smi:
        stereo_checks[mid] = st.stereo_inspect(smi)
write(RUN / "raw_model" / "stereo_checks.json", {
    "source_verified": False,
    "note": "Standard graph diagnostics only; source identity and axial/relative stereochemistry still require expert review.",
    "molecules": stereo_checks,
})

print(json.dumps({"molecules": len(bundle["molecules"]),
                  "reactions": len(bundle["reactions"]),
                  "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in bundle["reactions"]),
                  "source_context": len(bundle["source_context"])}, ensure_ascii=False))
