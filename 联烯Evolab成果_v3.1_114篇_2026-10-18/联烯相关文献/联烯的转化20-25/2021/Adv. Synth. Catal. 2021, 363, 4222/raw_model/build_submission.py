#!/usr/bin/env python3
"""Build the source-specific v3.1 bundle for doc_cb03cf42d5a5f24b.

The script only encodes facts manually transcribed from PDF pages 1--6.
It does not read legacy experiment records or infer an SI association.
"""
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DOC = "doc_cb03cf42d5a5f24b"


def ev(page, kind="image"):
    return f"{DOC}:page{page}:{kind}"


def F(raw=None, normalized=None, evidence=(), method="original_extraction",
      status=None, verification=None, derivation=None, alternatives=()):
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   verification=verification, derivation=derivation,
                   alternatives=alternatives)


def Q(raw, number, unit, evidence, qualifier=None, method="original_extraction"):
    value = {"number": number, "unit": unit}
    if qualifier:
        value["qualifier"] = qualifier
        value["raw_text"] = raw
    return F(raw, value, evidence, method=method)


def missing(evidence=()):
    return F(evidence=evidence)


def NA(evidence=(), reason="Not applicable for this source item."):
    return F(evidence=evidence, method="not_applicable", status="not_applicable",
             verification="missing", derivation=reason)


def inferred(normalized, evidence, derivation, raw=None):
    return F(raw, normalized, evidence, method="inference", verification="inferred",
             derivation=derivation)


def conflict(raw, alternatives, evidence):
    return F(raw, None, evidence, status="conflict", alternatives=alternatives)


template = st.read(RUN / "raw_model/standard_init.json")
bundle = template
bundle["molecules"] = []
bundle["reactions"] = []
bundle["source_context"] = []
molecules = {}
reactions = []


def add_molecule(mid, label, evidence, name=None, smiles=None, stereo="unresolved",
                 assignment=None, mixture=None, note=None):
    evidence = list(dict.fromkeys(evidence))
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, evidence)
    m["names"]["original"] = F(name or label, name or label, evidence)
    m["names"]["english"] = F(name or label, name or label, evidence)
    m["cas_numbers"] = missing()
    m["identity_note"] = F(note, note, evidence) if note else missing()
    if smiles:
        diag = st.stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = F({"source_drawing": label}, diag["canonical_isomeric_smiles"], evidence,
                         method="structure_derivation", verification="inferred",
                         derivation="Candidate graph transcribed from the source drawing; RDKit parsing/canonicalization does not verify identity.")
        m["stereochemistry"]["isomeric_smiles"] = F(None, diag["canonical_isomeric_smiles"], evidence,
                                                      method="structure_derivation", verification="inferred",
                                                      derivation="Canonical isomeric SMILES generated from the candidate graph; unassigned source stereo remains unassigned.")
        m["stereochemistry"]["identity_key"] = F(None, diag["stereo_graph_key"], evidence,
                                                   method="structure_derivation", verification="inferred",
                                                   derivation="RDKit stereo graph key for the candidate SMILES only.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing()
    if stereo == "allene_unassigned":
        m["stereochemistry"]["axial"] = F("allene drawn without R_a/S_a or P/M descriptor",
                                             [{"type": "allene", "configuration": "unspecified"}], evidence)
        m["stereochemistry"]["double_bonds"] = inferred(
            [{"type": "cumulated_allene", "configuration": "not_E_Z_assigned"}], evidence,
            "The source drawing contains a cumulated allene; E/Z labels are not used to replace its unassigned axial configuration.")
        m["stereochemistry"]["chirality_status"] = F("unassigned allene axis", "unspecified", evidence)
        m["stereochemistry"]["geometry_status"] = inferred("not_applicable", evidence,
                                                              "The relevant stereogenic element is an allene axis, not an alkene E/Z assignment.")
        m["stereochemistry"]["representation_status"] = inferred("partial", evidence,
                                                                    "Connectivity is candidate-only and the unassigned allene axis is not encoded by the toolkit.")
    elif stereo == "hydrazone_unassigned":
        m["stereochemistry"]["double_bonds"] = F("C=N drawn without E/Z descriptor",
                                                    [{"type": "hydrazone_C=N", "configuration": "unspecified"}], evidence)
        m["stereochemistry"]["chirality_status"] = inferred("unspecified", evidence,
                                                               "Complete stereogenic-element assessment is not established from the compact source drawing.")
        m["stereochemistry"]["geometry_status"] = F("no E/Z descriptor", "unspecified", evidence)
        m["stereochemistry"]["representation_status"] = inferred("partial" if smiles else "not_encoded", evidence,
                                                                    "Hydrazone geometry is unassigned in the source.")
    elif stereo == "dr_mixture":
        m["stereochemistry"]["mixture"] = F(mixture or "trans/cis diastereomeric mixture",
                                              mixture or "trans/cis diastereomeric mixture", evidence)
        m["stereochemistry"]["assignment"] = F("relative trans/cis labels; no absolute configuration",
                                                 {"relative": ["trans", "cis"], "absolute": "unreported"}, evidence)
        m["stereochemistry"]["chirality_status"] = F("diastereomeric mixture", "mixture", evidence)
        m["stereochemistry"]["geometry_status"] = inferred("not_applicable", evidence,
                                                              "The reported trans/cis descriptor concerns cyclopropane relative stereochemistry, not alkene E/Z.")
        m["stereochemistry"]["representation_status"] = inferred("not_encoded", evidence,
                                                                    "Individual relative stereoisomers are not encoded in a complete graph.")
    elif stereo == "relative_trans" or stereo == "relative_cis":
        rel = "trans" if stereo.endswith("trans") else "cis"
        m["stereochemistry"]["assignment"] = F(rel, {"relative_configuration": rel, "absolute": "unreported"}, evidence)
        m["stereochemistry"]["chirality_status"] = F(rel, "partially_specified", evidence)
        m["stereochemistry"]["geometry_status"] = inferred("not_applicable", evidence,
                                                              "Cyclopropane relative configuration is not alkene geometry.")
        m["stereochemistry"]["representation_status"] = inferred("not_encoded", evidence,
                                                                    "Relative configuration is retained textually because no complete isomeric graph was asserted.")
    elif stereo == "stereocenter_unassigned":
        m["stereochemistry"]["atom_centers"] = F("stereogenic center present; configuration unreported",
                                                   [{"configuration": "unspecified"}], evidence)
        m["stereochemistry"]["chirality_status"] = F("configuration unreported", "unspecified", evidence)
        m["stereochemistry"]["geometry_status"] = inferred("not_applicable", evidence,
                                                              "No alkene E/Z or other geometric descriptor is reported.")
        m["stereochemistry"]["representation_status"] = inferred("not_encoded", evidence,
                                                                    "The source does not assign an absolute configuration.")
    else:
        m["stereochemistry"]["chirality_status"] = inferred("unresolved", evidence,
                                                               "Center/axis/ring stereochemistry was considered but cannot be established from the supplied compact representation.")
        m["stereochemistry"]["geometry_status"] = inferred("unresolved", evidence,
                                                              "E/Z and ring-relative geometry cannot be established for this source item.")
        m["stereochemistry"]["representation_status"] = inferred("unresolved" if not smiles else "partial", evidence,
                                                                    "No source-verified complete stereochemical representation was asserted.")
    if assignment and stereo not in ("dr_mixture", "relative_trans", "relative_cis"):
        m["stereochemistry"]["assignment"] = F(assignment, assignment, evidence)
    molecules[mid] = m
    bundle["molecules"].append(m)
    return mid


def add_template(mid, label, evidence, fixed, variables, reason):
    add_molecule(mid, label, evidence, name=label, stereo="unresolved", note=reason)
    m = molecules[mid]
    m["representation_type"] = "generic_template"
    m["smiles"] = missing()
    m["generic_template"] = {
        "fixed_structure": F(fixed, fixed, evidence),
        "variables": [],
        "instance_mapping": [],
        "reason": F(reason, reason, evidence),
    }
    for label_, definition in variables:
        m["generic_template"]["variables"].append({
            "label": label_,
            "attachment_points": F("source drawing attachment", ["source-defined"], evidence),
            "definition": F(definition, definition, evidence),
            "allowed_values": missing(),
            "constraints": missing(),
        })
    return mid


def ensure_named(mid, label, evidence, name=None):
    if mid not in molecules:
        add_molecule(mid, label, evidence, name=name or label)
    return mid


# Core generic structures and source-labelled molecules.
add_template(f"{DOC}:generic_allenylphosphine_oxide_1", "allenylphosphine oxide 1",
             [ev(2), ev(4)], "ArO-CH2-C(P(O)Ph2)=C=C(H)-R",
             [("Ar", "aryl group"), ("R", "alkyl or aryl substituent")],
             "The paper depicts an R-group family; it is not a single concrete molecule.")
add_template(f"{DOC}:generic_tosylhydrazone_2", "N-tosylhydrazone 2", [ev(2), ev(4)],
             "R1R2C=N-NHTs", [("R1", "carbon substituent"), ("R2", "carbon substituent")],
             "The paper depicts an R-group family with variable R1 and R2.")
add_template(f"{DOC}:generic_product_3", "alkynylcyclopropane 3", [ev(2), ev(4)],
             "relative-trans/cis phosphinyl alkynylcyclopropane family",
             [("R", "alkynyl substituent"), ("R1", "cyclopropane substituent"), ("R2", "cyclopropane substituent")],
             "The source scheme is generic and includes unencoded relative stereochemistry.")
add_template(f"{DOC}:lit_rearrangement_substrate", "prior-work allenylphosphine oxide", [ev(2)],
             "allenylphosphine oxide bearing ArO and variable R/R'", [("R", "source-defined substituent"), ("R'", "source-defined substituent")],
             "Cited prior-work drawing is generic.")
add_template(f"{DOC}:lit_rearrangement_product", "prior-work rearrangement product", [ev(2)],
             "phosphinyl enyne product family", [("R", "source-defined substituent"), ("R'", "source-defined substituent")],
             "Cited prior-work product is generic.")
add_template(f"{DOC}:lit_cyclopropane_product", "prior-work cyclopropane product", [ev(2)],
             "phosphinyl cyclopropane family", [("R", "source-defined substituent"), ("R1", "source-defined substituent"), ("R2", "source-defined substituent")],
             "Cited prior-work product is generic.")

one_ids = {}
one_smiles = "Cc1cccc(C)c1OCC(P(=O)(c2ccccc2)c2ccccc2)=C=Cc1ccccc1"
one_subs = {
    "a": "Ph", "b": "4-OMe-Ph", "c": "4-Me-Ph", "d": "4-CF3-Ph",
    "e": "4-F-Ph", "f": "4-Cl-Ph", "g": "3-Me-Ph", "h": "2-Me-Ph",
    "i": "2-Naph", "j": "Et", "k": "Cy",
}
for suffix, sub in one_subs.items():
    mid = f"{DOC}:mol_1{suffix}"
    one_ids[suffix] = mid
    add_molecule(mid, f"1{suffix}", [ev(3), ev(4), ev(5)] if suffix == "a" else [ev(4)],
                 name=f"allenylphosphine oxide 1{suffix} (R={sub})",
                 smiles=one_smiles if suffix == "a" else None,
                 stereo="allene_unassigned",
                 note="The source gives no axial descriptor; missing assignment is not treated as achiral or racemic.")

two_ids = {}
two_defs = {
    "a": "R1=Me, R2=PMP", "b": "R1=Me, R2=Ph", "c": "R1=Me, R2=4-F-Ph",
    "d": "R1=Me, R2=4-Cl-Ph", "e": "R1=Me, R2=4-Br-Ph", "f": "R1=Me, R2=4-Me-Ph",
    "g": "R1=Me, R2=4-C6H4-Ph", "h": "R1=Me, R2=2-Naph", "i": "R1=Bn, R2=Ph",
    "j": "R1=Et, R2=Ph", "k": "R1=Et, R2=PMP", "l": "R1=Me, R2=3-Me-Ph",
    "m": "R1=R2=Ph", "n": "R1=R2=4-F-Ph", "o": "R1=R2=4-Cl-Ph",
    "p": "R1=R2=Bn", "q": "R1=R2=Et", "r": "indan-derived cyclic hydrazone",
    "s": "tetrahydronaphthalene-derived cyclic hydrazone",
}
two_smiles = "COc1ccc(C(C)=NNS(=O)(=O)c2ccc(C)cc2)cc1"
for suffix, desc in two_defs.items():
    mid = f"{DOC}:mol_2{suffix}"
    two_ids[suffix] = mid
    add_molecule(mid, f"2{suffix}", [ev(3), ev(4), ev(5)] if suffix == "a" else [ev(4)],
                 name=f"N-tosylhydrazone 2{suffix} ({desc})",
                 smiles=two_smiles if suffix == "a" else None,
                 stereo="hydrazone_unassigned",
                 note="Hydrazone C=N geometry is not assigned in the supplied pages.")

prod_ids = {}
a_products = ["aa", "ba", "ca", "da", "ea", "fa", "ga", "ha", "ia", "ja", "ka"]
b_products = ["ab", "ac", "ad", "ae", "af", "ag", "ah", "ai", "aj", "ak", "al", "am", "an", "ao", "ap", "aq", "ar", "as"]
dr_products = set(a_products + b_products[:11])
for code in a_products + b_products:
    mid = f"{DOC}:mol_3{code}"
    prod_ids[code] = mid
    if code in dr_products:
        stereo = "dr_mixture"
        mix = "reported trans/cis mixture; ratio is experiment-specific"
    elif code in {"ar", "as"}:
        stereo = "relative_trans"
        mix = None
    else:
        stereo = "stereocenter_unassigned"
        mix = None
    add_molecule(mid, f"3{code}" if code not in {"ar", "as"} else f"trans-3{code}",
                 [ev(3), ev(4), ev(5)] if code == "aa" else [ev(4)],
                 name=f"alkynylcyclopropane 3{code}", stereo=stereo, mixture=mix,
                 note="Complete graph was not asserted from the compact scope drawing; relative/absolute configuration limitations are retained.")

trans3aa = add_molecule(f"{DOC}:mol_trans_3aa", "trans-3aa", [ev(3), ev(4), ev(5)],
                        name="trans-3aa", stereo="relative_trans")
cis3aa = add_molecule(f"{DOC}:mol_cis_3aa", "cis-3aa", [ev(3), ev(4), ev(5)],
                      name="cis-3aa", stereo="relative_cis")
four_a = add_molecule(f"{DOC}:mol_4a", "4a", [ev(4), ev(5)],
                      name="conjugated enyne intermediate 4a", stereo="unresolved")
for label, desc in [("5a", "tosyl-substituted allene 5a"), ("6a", "isopropyl-substituted allene 6a"),
                    ("7a", "p-tolyl-substituted allene 7a")]:
    add_molecule(f"{DOC}:mol_{label}", label, [ev(4), ev(5)], name=desc,
                 stereo="allene_unassigned")
trans5aa = add_molecule(f"{DOC}:mol_trans_5aa", "trans-5aa", [ev(4)],
                        name="trans cyclopropane product 5aa", stereo="relative_trans")

# Named reaction, work-up, purification, and characterization species.
named = {
    "K2CO3": "potassium carbonate", "tBuOK": "potassium tert-butoxide",
    "KHCO3": "potassium bicarbonate", "Na2CO3": "sodium carbonate", "DBU": "DBU",
    "DMAP": "DMAP", "dioxane": "1,4-dioxane", "toluene": "toluene", "THF": "THF",
    "nitromethane": "nitromethane", "water": "water", "DMF": "DMF", "DMSO": "DMSO",
    "ethanol": "ethanol", "DCE": "DCE", "A1": "1,10-phenanthroline",
    "A2": "bathophenanthroline", "A3": "3,4,7,8-tetramethyl-1,10-phenanthroline",
    "A4": "4,7-dimethoxy-1,10-phenanthroline", "A5": "bathocuproine",
    "A6": "1,10-phenanthroline-5,6-dione", "TEPA": "triethyl phosphonoacetate",
    "DCM": "dichloromethane", "Na2SO4": "anhydrous sodium sulfate",
    "petroleum_ether": "petroleum ether", "ethyl_acetate": "ethyl acetate",
    "silica": "silica gel",
}
named_ids = {}
for key, name in named.items():
    evidence = [ev(3)] if key not in {"DCM", "Na2SO4", "petroleum_ether", "ethyl_acetate", "silica"} else [ev(5)]
    if key in {"K2CO3", "dioxane", "A1", "water", "DCE"}:
        evidence = [ev(3), ev(4), ev(5)]
    mid = f"{DOC}:mol_{key}"
    named_ids[key] = mid
    add_molecule(mid, key, evidence, name=name)


def reaction(rid, evidence, record_type, categories, location, name, transformation, unresolved=None):
    evidence = list(dict.fromkeys(evidence))
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = F(categories, categories, evidence) if categories is not None else missing()
    r["location"] = F(location, location, evidence)
    r["reaction_name"] = F(name, name, evidence) if name else missing()
    r["transformation_type"] = F(transformation, transformation, evidence) if transformation else missing()
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": list(unresolved or [])}
    r["completeness"] = "partial"
    reactions.append(r)
    bundle["reactions"].append(r)
    return r


def participant(r, pid, mid, role, evidence, stage="reaction", amount=None, mass=None,
                volume=None, equiv=None, molpct=None, reference=None, product_role=None):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    if product_role:
        p["product_role"] = F(product_role, product_role, evidence)
    if amount is not None:
        p["charges"]["amount"] = amount if isinstance(amount, dict) and "raw_value" in amount else Q(str(amount), amount, "mmol", evidence)
    if mass is not None:
        p["charges"]["mass"] = mass if isinstance(mass, dict) and "raw_value" in mass else mass
    if volume is not None:
        p["charges"]["volume"] = volume if isinstance(volume, dict) and "raw_value" in volume else volume
    if equiv is not None:
        p["charges"]["equivalents"] = equiv if isinstance(equiv, dict) and "raw_value" in equiv else Q(str(equiv), equiv, "equiv", evidence)
    if molpct is not None:
        p["charges"]["mol_percent"] = molpct if isinstance(molpct, dict) and "raw_value" in molpct else Q(str(molpct), molpct, "mol%", evidence)
    if reference:
        p["equivalence_reference"] = inferred(reference, evidence,
                                                "The source expresses equivalents/loading relative to the principal substrate.",
                                                raw=molecules[reference]["occurrences"][0]["compound_label"]["raw_value"])
    r["participants"].append(p)
    return p


def step(r, sid, stage, pids, evidence, operation=None, temperature=None, time=None,
         atmosphere=None, pressure=None):
    s = st.new_step(sid, evidence, stage=stage)
    s["participant_ids"] = list(pids)
    if operation is not None:
        s["operation"] = F(operation, operation, evidence)
    if temperature is not None:
        s["temperature"] = temperature if isinstance(temperature, dict) and "raw_value" in temperature else temperature
    if time is not None:
        s["time"] = time if isinstance(time, dict) and "raw_value" in time else time
    if atmosphere is not None:
        s["atmosphere"] = F(atmosphere, atmosphere, evidence)
    if pressure is not None:
        s["pressure"] = F(pressure, pressure, evidence)
    r["steps"].append(s)
    return s


def outcome(r, oid, kind, targets, evidence, value=None, yield_type=None,
            basis=None, scope=None, product_role=None):
    o = st.new_outcome(oid, kind, targets, evidence)
    if value is not None:
        o["value"] = value if isinstance(value, dict) and "raw_value" in value else F(value, value, evidence)
    if yield_type is not None:
        o["yield_type"] = F(yield_type, yield_type, evidence)
    if basis is not None:
        o["measurement_basis"] = F(basis, basis, evidence)
    if scope is not None:
        o["scope"] = F(scope, scope, evidence)
    if product_role is not None:
        o["product_role"] = F(product_role, product_role, evidence)
    r["outcomes"].append(o)
    return o


def classify(r, mode, evidence):
    evidence = list(dict.fromkeys(evidence))
    c = r["classification"]
    c["mapped_reaction_smiles"] = missing()
    c["model_votes"] = []
    if mode == "cascade":
        changes = [
            {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 2, "after_order": 3, "change_type": "order_increased"},
        ]
        labels = [
            {"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1", "b2"]},
            {"category_number": 21, "code": "TB-C-C", "supporting_change_ids": ["b3"]},
        ]
        deriv = "Partial class-level difference inferred from the source's cyclopropanation/[3+2] description and allene-to-alkyne product drawing; atom identities are not mapped."
        c["bond_changes"] = [inferred(x, evidence, deriv) for x in changes]
        c["labels"] = inferred(labels, evidence, deriv)
        c["coverage"] = inferred("partial", evidence, deriv)
        c["limitations"] = ["Exact atom identities, N/Ts-containing coproducts, all hydrogen changes, and complete atom mapping are unresolved.",
                              "Hydrogen changes were not guessed from an incomplete overall equation."]
    elif mode == "rearrangement":
        changes = [
            {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 1, "after_order": 2, "change_type": "order_increased"},
            {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 2, "after_order": 3, "change_type": "order_increased"},
        ]
        labels = [{"category_number": 13, "code": "DB-C-C", "supporting_change_ids": ["b1"]},
                  {"category_number": 21, "code": "TB-C-C", "supporting_change_ids": ["b2"]}]
        deriv = "Partial bond-order changes inferred from the source's allene-to-conjugated-enyne drawing; exact atom mapping and hydrogen loss are unresolved."
        c["bond_changes"] = [inferred(x, evidence, deriv) for x in changes]
        c["labels"] = inferred(labels, evidence, deriv)
        c["coverage"] = inferred("partial", evidence, deriv)
        c["limitations"] = ["Candidate structures are incomplete; exact atom identities and net hydrogen accounting are unresolved."]
    elif mode == "none":
        deriv = "No observed chemical transformation is classified; intended transformations are not assigned to negative controls."
        c["bond_changes"] = []
        c["labels"] = inferred([], evidence, deriv)
        c["coverage"] = inferred("not_applicable", evidence, deriv)
        c["limitations"] = ["Negative/no-detection experiment; no product bond changes were asserted."]
    else:
        c["bond_changes"] = []
        c["labels"] = missing()
        c["coverage"] = inferred("unresolved", evidence, "The compact cited/generic drawing was not assigned a complete bond-difference graph.")
        c["limitations"] = ["Complete structures and atom mapping are unavailable on the supplied pages."]


# Scheme 1 prior work and current generic transformation.
lit_a = reaction("scheme1_prior_rearrangement", [ev(2)], "literature_example", [],
                 {"page": 2, "scheme": "1B", "subpanel": "a"},
                 "Selective rearrangement of allenylphosphine oxides", "photocatalytic allene rearrangement",
                 ["Cited prior work; no concrete charge, yield, or complete structure is reported here."])
participant(lit_a, "prior_a_start", f"{DOC}:lit_rearrangement_substrate", "reactant", [ev(2)])
participant(lit_a, "prior_a_product", f"{DOC}:lit_rearrangement_product", "product", [ev(2)], product_role="reported prior-work product family")
classify(lit_a, "unknown", [ev(2)])

lit_b = reaction("scheme1_prior_cyclopropanation", [ev(2)], "literature_example", [],
                 {"page": 2, "scheme": "1B", "subpanel": "b"},
                 "Intermolecular cyclopropanation of allenylphosphine oxides", "cyclopropanation",
                 ["Cited prior work; no concrete charge, yield, or complete structure is reported here."])
participant(lit_b, "prior_b_start", f"{DOC}:lit_rearrangement_substrate", "reactant", [ev(2)])
participant(lit_b, "prior_b_hydrazone", f"{DOC}:generic_tosylhydrazone_2", "reactant", [ev(2)])
participant(lit_b, "prior_b_product", f"{DOC}:lit_cyclopropane_product", "product", [ev(2)], product_role="reported prior-work product family")
classify(lit_b, "unknown", [ev(2)])

generic = reaction("scheme1_current_generic", [ev(2), ev(4)], "generic_scheme", [],
                   {"pages": [2, 4], "scheme": "1C / Scheme 2 generic header"},
                   "Cascade enyne rearrangement and cyclopropanation", "allene-to-alkynylcyclopropane cascade",
                   ["R-group scheme is not a concrete experiment; exact atom mapping and hydrogen balance are unresolved."])
participant(generic, "generic_1", f"{DOC}:generic_allenylphosphine_oxide_1", "reactant", [ev(2), ev(4)])
participant(generic, "generic_2", f"{DOC}:generic_tosylhydrazone_2", "reactant", [ev(2), ev(4)])
participant(generic, "generic_3", f"{DOC}:generic_product_3", "product", [ev(2), ev(4)], product_role="generic trans/cis products")
participant(generic, "generic_A1", named_ids["A1"], "additive", [ev(2), ev(4)], molpct=20, reference=f"{DOC}:generic_allenylphosphine_oxide_1")
participant(generic, "generic_K2CO3", named_ids["K2CO3"], "reagent", [ev(2), ev(4)], equiv=3, reference=f"{DOC}:generic_allenylphosphine_oxide_1")
classify(generic, "cascade", [ev(2), ev(4)])
GENERIC_ID = generic["reaction_id"]


# Table 1: 24 independent optimization entries.
table_rows = [
    (1, None, "K2CO3", "dioxane", "4.9:1", 65),
    (2, None, "tBuOK", "dioxane", "3.7:1", 33),
    (3, None, "KHCO3", "dioxane", "2.0:1", 61),
    (4, None, "Na2CO3", "dioxane", "3.0:1", 60),
    (5, None, "DBU", "dioxane", "2.9:1", 40),
    (6, None, "DMAP", "dioxane", "3.3:1", 33),
    (7, None, "K2CO3", "toluene", "3.4:1", 42),
    (8, None, "K2CO3", "THF", "2.0:1", 38),
    (9, None, "K2CO3", "nitromethane", None, 0),
    (10, None, "K2CO3", "water", None, "<5"),
    (11, None, "K2CO3", "DMF", "3.3:1", 42),
    (12, None, "K2CO3", "DMSO", "3.0:1", 20),
    (13, None, "K2CO3", "ethanol", "3.6:1", 30),
    (14, None, "K2CO3", "DCE", "3.0:1", 13),
    (15, ("A1", 20), "K2CO3", "dioxane", "3.0:1", 91),
    (16, ("A2", 20), "K2CO3", "dioxane", "3.2:1", 70),
    (17, ("A3", 20), "K2CO3", "dioxane", "2.7:1", 64),
    (18, ("A4", 20), "K2CO3", "dioxane", "2.7:1", 67),
    (19, ("A5", 20), "K2CO3", "dioxane", "2.7:1", 70),
    (20, ("A6", 20), "K2CO3", "dioxane", "2.9:1", 80),
    (21, ("A1", 10), "K2CO3", "dioxane", "2.9:1", 80),
    (22, ("A1", 50), "K2CO3", "dioxane", "2.8:1", 72),
    (23, ("A1", 100), "K2CO3", "dioxane", "2.9:1", 70),
    (24, ("A1", 20), None, "dioxane", None, 0),
]
table_ids = []
for entry, additive, base, solvent, dr, yld in table_rows:
    rid = f"table1_entry{entry}"
    r = reaction(rid, [ev(3)], "specific_experiment", ["optimization"],
                 {"page": 3, "table": 1, "entry": entry},
                 "Optimization of the cascade enyne rearrangement and cyclopropanation",
                 "allene-to-alkynylcyclopropane cascade",
                 ["SI was advertised but not supplied; product structure and full mass balance remain partial."])
    r["general_scheme_ids"] = [GENERIC_ID]
    p1 = participant(r, "p_1a", one_ids["a"], "reactant", [ev(3)], amount=0.3, equiv=1.0, reference=one_ids["a"])
    p2 = participant(r, "p_2a", two_ids["a"], "reactant", [ev(3)], amount=0.6, equiv=2.0, reference=one_ids["a"])
    pp = participant(r, "p_3aa", prod_ids["aa"], "product", [ev(3)], product_role="desired trans/cis product mixture")
    if additive:
        participant(r, f"p_{additive[0]}", named_ids[additive[0]], "additive", [ev(3)],
                    molpct=additive[1], reference=one_ids["a"])
    if base:
        participant(r, f"p_{base}", named_ids[base], "reagent", [ev(3)], equiv=3.0, reference=one_ids["a"])
    participant(r, f"p_{solvent}", named_ids[solvent], "solvent", [ev(3)],
                volume=Q("3 mL", 3, "mL", [ev(3)]))
    if dr:
        participant(r, "p_TEPA", named_ids["TEPA"], "unknown", [ev(3)], stage="characterization")
    r["reaction_scale"] = F("1a, 0.3 mmol", {"reference_molecule_id": one_ids["a"], "amount": {"number": 0.3, "unit": "mmol"}}, [ev(3)])
    r["solvent_conditions"]["composition"] = F(named[solvent], [named_ids[solvent]], [ev(3)])
    r["solvent_conditions"]["total_volume"] = Q("3 mL", 3, "mL", [ev(3)])
    reaction_pids = [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"]
    step(r, "reaction", "reaction", reaction_pids, [ev(3)], operation="combine and heat",
         temperature=Q("110 °C", 110, "°C", [ev(3)]), time=Q("18 h", 18, "h", [ev(3)]), atmosphere="N2")
    if yld == "<5":
        yfact = Q("<5%", 5, "%", [ev(3)], qualifier="lt")
    else:
        yfact = Q(f"{yld}%", yld, "%", [ev(3)])
    outcome(r, "yield_3aa", "yield", [prod_ids["aa"]], [ev(3)], yfact,
            yield_type="isolated", basis="isolated yield", product_role="desired product mixture")
    if dr:
        outcome(r, "dr_3aa", "selectivity", [prod_ids["aa"]], [ev(3)],
                F(dr, {"ratio": dr, "order": ["trans", "cis"]}, [ev(3)]),
                basis="31P NMR spectroscopy of the products with TEPA as internal standard",
                scope={"metric": "dr", "order": ["trans", "cis"], "major": "trans"},
                product_role="trans/cis product mixture")
    classify(r, "none" if yld == 0 else "cascade", [ev(3)])
    table_ids.append(r["reaction_id"])


# Scheme 2 substrate scope, with the page-4/page-5 amount conflict retained.
scope_a = [
    ("aa", "a", 91, "3:1"), ("ba", "b", 77, "8.6:1"), ("ca", "c", 86, "1.9:1"),
    ("da", "d", 58, "1.8:1"), ("ea", "e", 90, "1.8:1"), ("fa", "f", 78, "1.5:1"),
    ("ga", "g", 79, "2.2:1"), ("ha", "h", 83, "1.3:1"), ("ia", "i", 80, "1.1:1"),
    ("ja", "j", 57, "3:1"), ("ka", "k", 76, "2.6:1"),
]
scope_b = [
    ("ab", "b", 73, "1.5:1"), ("ac", "c", 87, "2.3:1"), ("ad", "d", 71, "3.4:1"),
    ("ae", "e", 51, "3.6:1"), ("af", "f", 73, "3.3:1"), ("ag", "g", 61, "3.7:1"),
    ("ah", "h", 64, "3:1"), ("ai", "i", 58, "1.76:1"), ("aj", "j", 73, "1.1:1"),
    ("ak", "k", 64, "1.4:1"), ("al", "l", 59, "4.4:1"), ("am", "m", 58, None),
    ("an", "n", 65, None), ("ao", "o", 31, None), ("ap", "p", 26, None),
    ("aq", "q", 20, None), ("ar", "r", 32, "single trans diastereoisomer"),
    ("as", "s", 36, "single trans diastereoisomer"),
]
scope_ids = []


def scheme_amount(label, direct, general, evidence):
    return conflict([f"{direct} mmol (Scheme 2)", f"{general} mmol (General Procedure)"], [
        Q(f"{direct} mmol", direct, "mmol", [ev(4)]),
        Q(f"{general} mmol", general, "mmol", [ev(5)], method="inherited_general_condition"),
    ], evidence)


def add_scope(code, one_suffix, two_suffix, yld, dr):
    evidence = [ev(4), ev(5)]
    rid = f"scheme2_3{code}"
    unresolved = [
        "Scheme 2 uses 0.2/0.4 mmol while the supplied General Procedure uses 0.3/0.6 mmol; both are retained without choosing one.",
        "Complete product structure and exact atom mapping remain unresolved; SI was not supplied.",
    ]
    if dr:
        unresolved.append("The Scheme 2 dr/single-diastereoisomer measurement method is unreported on pages 1-6.")
    r = reaction(rid, evidence, "specific_experiment", ["substrate_scope"],
                 {"page": 4, "scheme": 2, "product": f"3{code}", "procedure_page": 5},
                 "Cascade enyne rearrangement and cyclopropanation",
                 "allene-to-alkynylcyclopropane cascade", unresolved)
    r["general_scheme_ids"] = [GENERIC_ID]
    p1 = participant(r, f"p_1{one_suffix}", one_ids[one_suffix], "reactant", evidence,
                     amount=scheme_amount(f"1{one_suffix}", 0.2, 0.3, evidence), equiv=1.0,
                     reference=one_ids[one_suffix])
    p2 = participant(r, f"p_2{two_suffix}", two_ids[two_suffix], "reactant", evidence,
                     amount=scheme_amount(f"2{two_suffix}", 0.4, 0.6, evidence), equiv=2.0,
                     reference=one_ids[one_suffix])
    pp = participant(r, f"p_3{code}", prod_ids[code], "product", evidence,
                     product_role="desired alkynylcyclopropane")
    pphen = participant(r, "p_A1", named_ids["A1"], "additive", evidence,
                        mass=Q("11.9 mg", 11.9, "mg", [ev(5)], method="inherited_general_condition"),
                        molpct=20, reference=one_ids[one_suffix])
    pk = participant(r, "p_K2CO3", named_ids["K2CO3"], "reagent", evidence,
                     mass=Q("124 mg", 124, "mg", [ev(5)], method="inherited_general_condition"),
                     amount=Q("0.9 mmol", 0.9, "mmol", [ev(5)], method="inherited_general_condition"),
                     equiv=3.0, reference=one_ids[one_suffix])
    pdiox = participant(r, "p_dioxane", named_ids["dioxane"], "solvent", evidence,
                        volume=Q("3 mL", 3, "mL", [ev(5)], method="inherited_general_condition"))
    for pid, key, role, stage_ in [
        ("p_DCM", "DCM", "workup_agent", "workup"), ("p_water", "water", "workup_agent", "workup"),
        ("p_Na2SO4", "Na2SO4", "workup_agent", "workup"),
        ("p_petroleum_ether", "petroleum_ether", "purification_agent", "purification"),
        ("p_ethyl_acetate", "ethyl_acetate", "purification_agent", "purification"),
        ("p_DCM_purification", "DCM", "purification_agent", "purification"),
        ("p_silica", "silica", "purification_agent", "purification"),
    ]:
        participant(r, pid, named_ids[key], role, [ev(5)], stage=stage_)
    alt_scale1 = F(f"1{one_suffix}, 0.2 mmol", {"reference_molecule_id": one_ids[one_suffix], "amount": {"number": 0.2, "unit": "mmol"}}, [ev(4)])
    alt_scale2 = F(f"1, 0.3 mmol", {"reference_molecule_id": one_ids[one_suffix], "amount": {"number": 0.3, "unit": "mmol"}}, [ev(5)], method="inherited_general_condition")
    r["reaction_scale"] = conflict([alt_scale1["raw_value"], alt_scale2["raw_value"]], [alt_scale1, alt_scale2], evidence)
    r["solvent_conditions"]["composition"] = F("degassed 1,4-dioxane", [named_ids["dioxane"]], [ev(5)], method="inherited_general_condition")
    r["solvent_conditions"]["total_volume"] = Q("3 mL", 3, "mL", [ev(5)], method="inherited_general_condition")
    r["special_conditions"]["reactor"] = F("10 mL two-necked flask equipped with condenser and stirring bar",
                                                {"type": "two-necked flask", "volume": {"number": 10, "unit": "mL"}, "condenser": True},
                                                [ev(5)], method="inherited_general_condition")
    rxn_pids = [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"]
    step(r, "reaction", "reaction", rxn_pids, evidence,
         operation="charge flask under nitrogen and heat in an oil bath until complete consumption of 1 monitored by TLC",
         temperature=Q("110 °C", 110, "°C", evidence), time=Q("18 h", 18, "h", evidence), atmosphere="N2")
    step(r, "workup_remove_volatiles", "workup", ["p_DCM", "p_water"], [ev(5)],
         operation="cool to room temperature and remove all volatiles under vacuum",
         temperature=F("room temperature", "room temperature", [ev(5)]))
    step(r, "workup_extract", "workup", ["p_DCM", "p_water"], [ev(5)],
         operation="extract with DCM and H2O")
    step(r, "workup_dry", "workup", ["p_Na2SO4"], [ev(5)],
         operation="combine organic layers, dry over anhydrous Na2SO4, and evaporate in vacuo")
    step(r, "purification", "purification", ["p_petroleum_ether", "p_ethyl_acetate", "p_DCM_purification", "p_silica"], [ev(5)],
         operation="flash chromatography on silica gel using petroleum ether/ethyl acetate/dichloromethane")
    r["procedure"]["original_text"] = F("General Procedure for compounds (3), page 5", "page-5 General Procedure", [ev(5)], method="inherited_general_condition")
    r["procedure"]["general_procedure_reference"] = F("General Procedure for compounds (3)", "page-5 General Procedure", [ev(5)], method="inherited_general_condition")
    r["procedure"]["addition_order"] = F("flask charged with 1, 2, 1,10-Phen, K2CO3 and degassed 1,4-dioxane",
                                             [f"1{one_suffix}", f"2{two_suffix}", "1,10-Phen", "K2CO3", "degassed 1,4-dioxane"],
                                             [ev(5)], method="inherited_general_condition")
    outcome(r, f"yield_3{code}", "yield", [prod_ids[code]], [ev(4)], Q(f"{yld}%", yld, "%", [ev(4)]),
            yield_type="isolated", basis="isolated yield by column chromatography",
            product_role="desired alkynylcyclopropane")
    if dr and dr.startswith("single"):
        outcome(r, f"selectivity_3{code}", "selectivity", [prod_ids[code]], [ev(4)],
                F(dr, {"metric": "diastereomer identity", "result": "single trans diastereoisomer"}, [ev(4)]),
                scope={"metric": "diastereoselectivity", "major": "trans", "minor": "not reported"},
                product_role="trans product")
    elif dr:
        outcome(r, f"dr_3{code}", "selectivity", [prod_ids[code]], [ev(4)],
                F(dr, {"ratio": dr, "order": ["trans", "cis"]}, [ev(4)]),
                scope={"metric": "dr", "order": ["trans", "cis"], "major": "trans"},
                product_role="trans/cis product mixture")
    classify(r, "cascade", [ev(4), ev(5)])
    r["notes"].append(F("Source conflict: Scheme 2 uses 1/2 = 0.2/0.4 mmol; the General Procedure uses 0.3/0.6 mmol.",
                          "Both source alternatives retained; no preferred amount selected.", evidence))
    scope_ids.append(r["reaction_id"])


for code, one_suffix, yld, dr in scope_a:
    add_scope(code, one_suffix, "a", yld, dr)
for code, two_suffix, yld, dr in scope_b:
    add_scope(code, "a", two_suffix, yld, dr)


# Scheme 3 scale-up.
scale = reaction("scheme3_scale_up", [ev(4)], "specific_experiment", ["scale_up"],
                 {"page": 4, "scheme": 3}, "Gram-scale cascade enyne rearrangement and cyclopropanation",
                 "allene-to-alkynylcyclopropane cascade",
                 ["The scale-up yield/selectivity measurement methods and detailed work-up are not reported on pages 1-6."])
scale["general_scheme_ids"] = [GENERIC_ID]
participant(scale, "p_1a", one_ids["a"], "reactant", [ev(4)], amount=2.5,
            mass=Q("1.10 g", 1.10, "g", [ev(4)]), equiv=1.0, reference=one_ids["a"])
participant(scale, "p_2a", two_ids["a"], "reactant", [ev(4)], amount=5.0,
            mass=Q("1.59 g", 1.59, "g", [ev(4)]), equiv=2.0, reference=one_ids["a"])
participant(scale, "p_A1", named_ids["A1"], "additive", [ev(3), ev(4)], molpct=20, reference=one_ids["a"])
participant(scale, "p_K2CO3", named_ids["K2CO3"], "reagent", [ev(3), ev(4)], equiv=3.0, reference=one_ids["a"])
participant(scale, "p_dioxane", named_ids["dioxane"], "solvent", [ev(3), ev(4)])
participant(scale, "p_trans3aa", trans3aa, "product", [ev(4)], product_role="major trans product")
participant(scale, "p_cis3aa", cis3aa, "product", [ev(4)], product_role="minor cis product")
scale["reaction_scale"] = F("1a, 2.5 mmol", {"reference_molecule_id": one_ids["a"], "amount": {"number": 2.5, "unit": "mmol"}}, [ev(4)])
scale["solvent_conditions"]["composition"] = F("standard condition: 1,4-dioxane", [named_ids["dioxane"]], [ev(3), ev(4)], method="inherited_general_condition")
step(scale, "reaction", "reaction", [p["participant_id"] for p in scale["participants"] if p["stage"] == "reaction"], [ev(3), ev(4)],
     operation="conduct under the standard condition", temperature=Q("110 °C", 110, "°C", [ev(3)], method="inherited_general_condition"),
     time=Q("18 h", 18, "h", [ev(3)], method="inherited_general_condition"), atmosphere="N2")
outcome(scale, "yield_trans3aa", "yield", [trans3aa], [ev(4)], Q("62%", 62, "%", [ev(4)]), product_role="trans-3aa")
outcome(scale, "mass_trans3aa", "product_amount", [trans3aa], [ev(4)], Q("715 mg", 715, "mg", [ev(4)]), product_role="trans-3aa")
outcome(scale, "yield_cis3aa", "yield", [cis3aa], [ev(4)], Q("21%", 21, "%", [ev(4)]), product_role="cis-3aa")
outcome(scale, "mass_cis3aa", "product_amount", [cis3aa], [ev(4)], Q("246 mg", 246, "mg", [ev(4)]), product_role="cis-3aa")
outcome(scale, "total_yield", "yield", [trans3aa, cis3aa], [ev(4)], Q("83%", 83, "%", [ev(4)]), product_role="combined trans/cis products")
outcome(scale, "dr", "selectivity", [trans3aa, cis3aa], [ev(4)],
        F("2.9:1", {"ratio": "2.9:1", "order": [trans3aa, cis3aa]}, [ev(4)]),
        scope={"metric": "dr", "order": ["trans", "cis"], "major": "trans"})
classify(scale, "cascade", [ev(4)])


# Scheme 4 controls.
ctrl_ids = []
ctrl1 = reaction("scheme4_eq_a_K2CO3", [ev(4)], "specific_experiment", ["control", "mechanistic_study"],
                 {"page": 4, "scheme": 4, "equation": "a", "branch": "K2CO3"},
                 "Base-promoted formation of conjugated enyne 4a", "allene-to-enyne rearrangement",
                 ["Atmosphere, solvent volume, detailed isolation, and complete structures are unreported on pages 1-6."])
participant(ctrl1, "p_1a", one_ids["a"], "reactant", [ev(4)], amount=0.2, reference=one_ids["a"])
participant(ctrl1, "p_K2CO3", named_ids["K2CO3"], "reagent", [ev(4)], equiv=1.0, reference=one_ids["a"])
participant(ctrl1, "p_dioxane", named_ids["dioxane"], "solvent", [ev(4)])
participant(ctrl1, "p_4a", four_a, "product", [ev(4)], product_role="isolated conjugated enyne")
ctrl1["reaction_scale"] = F("1a, 0.2 mmol", {"reference_molecule_id": one_ids["a"], "amount": {"number": 0.2, "unit": "mmol"}}, [ev(4)])
ctrl1["solvent_conditions"]["composition"] = F("1,4-dioxane", [named_ids["dioxane"]], [ev(4)])
step(ctrl1, "reaction", "reaction", [p["participant_id"] for p in ctrl1["participants"]], [ev(4)],
     operation="treat 1a with K2CO3", temperature=Q("110 °C", 110, "°C", [ev(4)]), time=Q("2 h", 2, "h", [ev(4)]))
outcome(ctrl1, "yield_4a", "yield", [four_a], [ev(4)], Q("85%", 85, "%", [ev(4)]),
        yield_type="isolated", basis="isolated yield", product_role="conjugated enyne intermediate")
classify(ctrl1, "rearrangement", [ev(4)])
ctrl_ids.append(ctrl1["reaction_id"])

ctrl2 = reaction("scheme4_eq_a_phenanthroline_only", [ev(4)], "specific_experiment", ["control", "mechanistic_study"],
                 {"page": 4, "scheme": 4, "equation": "a", "branch": "1,10-phenanthroline"},
                 "Phenanthroline-only control", "no observed reaction",
                 ["N.R. is qualitative; recovery, detection method, and product identity are unreported."])
participant(ctrl2, "p_1a", one_ids["a"], "reactant", [ev(4)], amount=0.2, reference=one_ids["a"])
participant(ctrl2, "p_A1", named_ids["A1"], "additive", [ev(4)], molpct=20, reference=one_ids["a"])
participant(ctrl2, "p_dioxane", named_ids["dioxane"], "solvent", [ev(4)])
ctrl2["reaction_scale"] = F("1a, 0.2 mmol", {"reference_molecule_id": one_ids["a"], "amount": {"number": 0.2, "unit": "mmol"}}, [ev(4)])
ctrl2["solvent_conditions"]["composition"] = F("1,4-dioxane", [named_ids["dioxane"]], [ev(4)])
step(ctrl2, "reaction", "reaction", [p["participant_id"] for p in ctrl2["participants"]], [ev(4)],
     operation="treat 1a with 1,10-phenanthroline only", temperature=Q("110 °C", 110, "°C", [ev(4)]), time=Q("18 h", 18, "h", [ev(4)]))
outcome(ctrl2, "no_reaction", "no_reaction", [], [ev(4)], F("N.R.", "no reaction reported", [ev(4)]))
classify(ctrl2, "none", [ev(4)])
ctrl_ids.append(ctrl2["reaction_id"])

ctrl3 = reaction("scheme4_eq_b_4a_plus_2a", [ev(4)], "specific_experiment", ["control", "mechanistic_study"],
                 {"page": 4, "scheme": 4, "equation": "b"},
                 "Reaction of proposed intermediate 4a with 2a", "cyclopropanation of conjugated enyne",
                 ["The dr measurement method and complete atom mapping are unreported."])
participant(ctrl3, "p_4a", four_a, "reactant", [ev(4)], amount=0.1, equiv=1.0, reference=four_a)
participant(ctrl3, "p_2a", two_ids["a"], "reactant", [ev(4)], amount=0.2, equiv=2.0, reference=four_a)
participant(ctrl3, "p_A1", named_ids["A1"], "additive", [ev(4)], molpct=20, reference=four_a)
participant(ctrl3, "p_K2CO3", named_ids["K2CO3"], "reagent", [ev(4)], equiv=2.0, reference=four_a)
participant(ctrl3, "p_dioxane", named_ids["dioxane"], "solvent", [ev(4)])
participant(ctrl3, "p_3aa", prod_ids["aa"], "product", [ev(4)], product_role="trans/cis product mixture")
ctrl3["reaction_scale"] = F("4a, 0.1 mmol", {"reference_molecule_id": four_a, "amount": {"number": 0.1, "unit": "mmol"}}, [ev(4)])
ctrl3["solvent_conditions"]["composition"] = F("1,4-dioxane", [named_ids["dioxane"]], [ev(4)])
step(ctrl3, "reaction", "reaction", [p["participant_id"] for p in ctrl3["participants"]], [ev(4)],
     operation="react 4a with 2a", temperature=Q("110 °C", 110, "°C", [ev(4)]), time=Q("18 h", 18, "h", [ev(4)]))
outcome(ctrl3, "yield_3aa", "yield", [prod_ids["aa"]], [ev(4)], Q("86%", 86, "%", [ev(4)]),
        yield_type="isolated", basis="isolated yield", product_role="desired product mixture")
outcome(ctrl3, "dr_3aa", "selectivity", [prod_ids["aa"]], [ev(4)],
        F("2.9:1", {"ratio": "2.9:1", "order": ["trans", "cis"]}, [ev(4)]),
        scope={"metric": "dr", "order": ["trans", "cis"], "major": "trans"})
classify(ctrl3, "cascade", [ev(4)])
ctrl_ids.append(ctrl3["reaction_id"])


def eqc_control(label, product_mid=None, yld=None):
    evidence = [ev(4), ev(5)] if label in {"6a", "7a"} else [ev(4)]
    r = reaction(f"scheme4_eq_c_{label}", evidence, "specific_experiment", ["control", "mechanistic_study"],
                 {"page": 4, "scheme": 4, "equation": "c", "substrate": label},
                 f"Cyclization control with {label}", "cyclopropanation control",
                 ["Standard-condition quantities are inherited from Table 1 entry 15; the Scheme 4 arrow does not repeat them.",
                  "Complete substrate/product structures and exact atom mapping are unresolved."])
    start = f"{DOC}:mol_{label}"
    participant(r, f"p_{label}", start, "reactant", evidence, amount=Q("0.3 mmol (standard condition)", 0.3, "mmol", [ev(3)], method="inherited_general_condition"), equiv=1.0, reference=start)
    participant(r, "p_2a", two_ids["a"], "reactant", evidence, amount=Q("0.6 mmol (standard condition)", 0.6, "mmol", [ev(3)], method="inherited_general_condition"), equiv=2.0, reference=start)
    participant(r, "p_A1", named_ids["A1"], "additive", evidence, molpct=20, reference=start)
    participant(r, "p_K2CO3", named_ids["K2CO3"], "reagent", evidence, equiv=3.0, reference=start)
    participant(r, "p_dioxane", named_ids["dioxane"], "solvent", evidence, volume=Q("3 mL (standard condition)", 3, "mL", [ev(3)], method="inherited_general_condition"))
    if product_mid:
        participant(r, "p_product", product_mid, "product", evidence, product_role="reported trans product")
    r["reaction_scale"] = F(f"{label}, 0.3 mmol (standard condition)", {"reference_molecule_id": start, "amount": {"number": 0.3, "unit": "mmol"}}, [ev(3)], method="inherited_general_condition")
    r["solvent_conditions"]["composition"] = F("1,4-dioxane (standard condition)", [named_ids["dioxane"]], [ev(3)], method="inherited_general_condition")
    r["solvent_conditions"]["total_volume"] = Q("3 mL (standard condition)", 3, "mL", [ev(3)], method="inherited_general_condition")
    step(r, "reaction", "reaction", [p["participant_id"] for p in r["participants"]], [ev(3), ev(4)],
         operation="apply standard condition", temperature=Q("110 °C", 110, "°C", [ev(3)], method="inherited_general_condition"),
         time=Q("18 h", 18, "h", [ev(3)], method="inherited_general_condition"), atmosphere="N2")
    if product_mid:
        outcome(r, "yield_product", "yield", [product_mid], [ev(4)], Q(f"{yld}%", yld, "%", [ev(4)]),
                yield_type="isolated", basis="isolated yield", product_role="trans product")
        classify(r, "cascade", evidence)
    else:
        outcome(r, "not_detected", "no_reaction", [], [ev(4)], F("N.D.", "proposed corresponding product not detected", [ev(4)]))
        outcome(r, "rapid_decomposition", "decomposition", [start], [ev(5)], F("decomposed rapidly", True, [ev(5)]))
        classify(r, "none", evidence)
    ctrl_ids.append(r["reaction_id"])
    return r


eqc_control("5a", trans5aa, 25)
eqc_control("6a")
eqc_control("7a")


# Scheme 5 proposed mechanism.
mechanism = reaction("scheme5_proposed_mechanism", [ev(5)], "proposed_mechanism", [],
                     {"page": 5, "scheme": 5}, "Proposed mechanism for the cascade enyne rearrangement and cyclopropanation",
                     "proposed SET/diazo/[3+2]/N2-extrusion pathway",
                     ["This is an author-proposed mechanism, not a set of executed synthetic experiments.",
                      "Intermediates I-IV and A-E are not treated as isolated compounds; complete atom mapping and hydrogen/coproduct balance are unresolved."])
participant(mechanism, "p_1a", one_ids["a"], "reactant", [ev(5)])
participant(mechanism, "p_2a", two_ids["a"], "reactant", [ev(5)])
participant(mechanism, "p_4a", four_a, "unknown", [ev(5)])
participant(mechanism, "p_3aa", prod_ids["aa"], "product", [ev(5)], product_role="proposed overall trans/cis products")
mechanism["notes"].append(F("ΔE = E(D)-E(E) = -0.91 kJ/mol; cis-3aa is 8.70 kJ/mol higher than trans-3aa",
                              {"D_minus_E": {"number": -0.91, "unit": "kJ/mol"}, "cis_minus_trans": {"number": 8.70, "unit": "kJ/mol"}}, [ev(5)]))
classify(mechanism, "cascade", [ev(5)])


# Non-synthetic and negative inventory context.
bundle["source_context"] = [
    {"context_id": f"{DOC}:context_page1_figure1", "document_id": DOC,
     "evidence_ids": [ev(1)], "kind": "biological_context",
     "description": "Figure 1 shows medicinally relevant alkynylcyclopropane scaffolds; it is not a current-paper synthesis experiment."},
    {"context_id": f"{DOC}:context_table1_additive_structures", "document_id": DOC,
     "evidence_ids": [ev(3)], "kind": "additive_identity_context",
     "description": "Structures and R-group definitions for additives A1-A6 support Table 1 identities."},
    {"context_id": f"{DOC}:context_xray_3ba", "document_id": DOC,
     "evidence_ids": [ev(2), ev(4)], "kind": "structure_assignment",
     "description": "The trans relative stereochemistry of 3ba is assigned from X-ray crystallography, CCDC 2034008; external CIF data were not supplied."},
    {"context_id": f"{DOC}:context_pages1_6_negative_spectral_inventory", "document_id": DOC,
     "evidence_ids": [ev(i) for i in range(1, 7)], "kind": "negative_inventory",
     "description": "No standalone mass-spectrum peak table, molecular-weight table, or fragment-relative-abundance block occurs on pages 1-6; references on page 6 were not converted to experiments."},
]


# Coverage with exactly the requested six PDF pages.
page2_ids = [lit_a["reaction_id"], lit_b["reaction_id"], generic["reaction_id"]]
page3_ids = table_ids
page4_ids = scope_ids + [scale["reaction_id"]] + ctrl_ids
page5_ids = scope_ids + [rid for rid in ctrl_ids if rid.endswith(("6a", "7a"))] + [mechanism["reaction_id"]]
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title/abstract/keywords/SI notice", "finding": "Paper metadata and online-SI notice; no SI file was supplied or associated.", "evidence_ids": [ev(1)]},
            {"region": "Figure 1 and introduction", "finding": "Medicinal alkynylcyclopropane examples and narrative context; not a current synthesis experiment.", "evidence_ids": [ev(1)]}],
         "record_ids": [], "unresolved": ["No reaction arrow or concrete current-paper experiment on this page."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Scheme 1A", "finding": "Conceptual literature-method map; no individual experiment was created from the schematic taxonomy.", "evidence_ids": [ev(2)]},
            {"region": "Scheme 1B panels a-b", "finding": "Two cited prior-work reaction families recorded as literature_example records.", "evidence_ids": [ev(2)]},
            {"region": "Scheme 1C and Results prose", "finding": "Current R-group generic transformation; optimization and scope are introduced.", "evidence_ids": [ev(2)]}],
         "record_ids": page2_ids, "unresolved": ["Scheme 1 prior-work structures are generic and no experimental charges/yields are printed."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Table 1 entries 1-24", "finding": "All 24 optimization rows recorded independently, including 0%, <5%, additive-loading footnotes, no-base control, isolated-yield basis, and dr method.", "evidence_ids": [ev(3)]},
            {"region": "Table 1 additive structures A1-A6", "finding": "Identity/context only; structures support participant labels rather than separate syntheses.", "evidence_ids": [ev(3)]},
            {"region": "bottom discussion", "finding": "Scope narrative and yields for symmetric/spiro cases; linked to Scheme 2 records on page 4.", "evidence_ids": [ev(3)]}],
         "record_ids": page3_ids, "unresolved": ["Compact product 3aa graph and complete atom mapping remain unencoded."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "Scheme 2A", "finding": "Eleven allenylphosphine-oxide scope outcomes 3aa-3ka recorded separately.", "evidence_ids": [ev(4)]},
            {"region": "Scheme 2B", "finding": "Eighteen N-tosylhydrazone scope outcomes 3ab-3as recorded separately; trans-3ar/3as retained as single diastereoisomers.", "evidence_ids": [ev(4)]},
            {"region": "Scheme 3", "finding": "One 2.5 mmol scale-up with separately reported trans/cis masses and yields, total yield, and dr.", "evidence_ids": [ev(4)]},
            {"region": "Scheme 4 equations a-c", "finding": "Six concrete controls recorded: two equation-a branches, one 4a+2a test, and three substituted-allene tests.", "evidence_ids": [ev(4)]}],
         "record_ids": page4_ids, "unresolved": ["Scope/general-procedure 0.2/0.4 versus 0.3/0.6 mmol conflict retained; several compact structures and Scheme 2 dr methods remain unresolved."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "Scheme 5 and mechanism discussion", "finding": "One proposed_mechanism record; intermediates were not expanded into executed experiments.", "evidence_ids": [ev(5)]},
            {"region": "top-left continuation", "finding": "6a and 7a are reported to decompose rapidly without detected proposed products; linked to page-4 controls.", "evidence_ids": [ev(5)]},
            {"region": "General Procedure", "finding": "Reactor, all charges, N2/110 C/18 h, TLC endpoint, workup, drying, evaporation, and chromatography associated with Scheme 2 product records; amount conflict retained.", "evidence_ids": [ev(5)]}],
         "record_ids": page5_ids, "unresolved": ["The General Procedure does not identify a single product label; association to compounds 3 is retained as inherited general condition, not a source-reviewed claim."]},
        {"page": 6, "state": "attempted", "regions": [
            {"region": "References [1]-[8]", "finding": "Bibliographic references only; no current-paper synthesis arrow, procedure, or table row.", "evidence_ids": [ev(6)]}],
         "record_ids": [], "unresolved": ["Page 7 and later were not part of the requested coverage and are not claimed as reviewed."]},
    ],
}

raw_inventory = {
    "document_id": DOC,
    "pages_attempted": [1, 2, 3, 4, 5, 6],
    "source_images": [{"page": i, "evidence_id": ev(i), "asset": f"assets/{DOC}/page{i}.jpg"} for i in range(1, 7)],
    "inventory": {
        "page1": {"figures": ["Figure 1 medicinal context"], "current_experiments": 0},
        "page2": {"schemes": ["Scheme 1A conceptual map", "Scheme 1B two prior-work families", "Scheme 1C current generic scheme"]},
        "page3": {"tables": [{"table": 1, "entries": list(range(1, 25))}], "current_experiments": 24},
        "page4": {"scheme2_products": [f"3{x}" for x in a_products + b_products], "scheme3": ["scale-up"], "scheme4_controls": ["eq a/base", "eq a/phen", "eq b", "eq c/5a", "eq c/6a", "eq c/7a"]},
        "page5": {"schemes": ["Scheme 5 proposed mechanism"], "procedure": "General Procedure for compounds 3"},
        "page6": {"content": "references only", "current_experiments": 0},
    },
    "deduplication": [
        "Table 1 entry 15 and Scheme 2 product 3aa were not merged because their printed substrate amounts differ (0.3/0.6 vs 0.2/0.4 mmol).",
        "Scheme 2 rows and the page-5 General Procedure were linked as repeated/general evidence, not counted as duplicate experiments.",
    ],
    "si_association": "unestablished; no SI file supplied",
}

st.write(RUN / "submitted.json", bundle)
st.write(RUN / "coverage/pages_0001_0006.json", coverage)
st.write(RUN / "raw_model/page_extraction.json", raw_inventory)
st.write(RUN / "raw_model/stereo_checks.json", {
    "1a": st.stereo_inspect(one_smiles),
    "2a": st.stereo_inspect(two_smiles),
    "limitations": "Graph diagnostics only; source identity and unassigned allene/hydrazone stereochemistry remain needs_review.",
})
st.write(RUN / "raw_model/verification_summary.json", {
    "expected_reactions": 64,
    "actual_reactions": len(bundle["reactions"]),
    "expected_concrete_experiments": 60,
    "scope_amount_conflict_records": scope_ids,
    "markdown_tokens_to_check": ["11.9", "124", "0.9", "3 mL", "110", "18 h", "DCM", "H2O", "Na2SO4", "petroleum ether", "ethyl acetate", "dichloromethane", "Source conflict"],
})
print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
                  "concrete_experiments": sum(r["record_type"] == "specific_experiment" for r in bundle["reactions"])}, ensure_ascii=False))
