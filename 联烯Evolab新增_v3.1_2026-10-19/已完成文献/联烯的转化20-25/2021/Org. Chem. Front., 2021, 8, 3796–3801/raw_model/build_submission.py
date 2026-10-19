#!/usr/bin/env python3
"""Build the source-derived v3.1 bundle for doc_bee5eb068d068a29.

This is a machine-authored extraction artifact. It intentionally preserves
unencoded graphical structures and source conflicts for expert review.
"""
from pathlib import Path
import copy
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DATASET = json.loads((RUN / "dataset.json").read_text())
DOC = "doc_bee5eb068d068a29"
E = {p: [f"{DOC}:page{p}:text", f"{DOC}:page{p}:image"] for p in range(1, 7)}


def F(raw=None, norm=None, ev=(), method="original_extraction", derivation=None,
      status=None, verification=None, alternatives=()):
    return st.fact(raw, norm, ev, method=method, derivation=derivation,
                   status=status, verification=verification, alternatives=alternatives)


def inferred(norm, ev, why):
    return F(None, norm, ev, method="inference", derivation=why)


def missing(ev=()):
    return F(ev=ev)


def source_molecule(mid, label, ev, kind="named", name=None):
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, ev)
    m["names"]["original"] = F(name or label, name or label, ev)
    m["names"]["english"] = F(name or label, name or label, ev)
    m["identity_note"] = F(
        "Identity is retained from the paper-local label/name; no CAS number is reported.",
        "Paper-local identity only; external identity not asserted.", ev)
    if kind == "substrate":
        m["stereochemistry"]["axial"] = F(
            "The drawn ene-vinylidenecyclopropane contains a VDCP/allene axis; no R_a/S_a or P/M descriptor is printed.",
            {"axis": "vinylidenecyclopropane/allene", "configuration": "unspecified"}, ev)
        m["stereochemistry"]["chirality_status"] = F("no configuration printed", "unspecified", ev)
        m["stereochemistry"]["geometry_status"] = F("alkene/allene geometry not assigned", "unspecified", ev)
        m["stereochemistry"]["representation_status"] = F("graphical structure not encoded as SMILES", "not_encoded", ev)
    elif kind == "product":
        m["stereochemistry"]["assignment"] = F(
            "Product drawing contains wedge/dashed stereobonds, but no absolute R/S descriptor is printed.",
            {"type": "depicted relative stereochemistry", "absolute_assignment": "unreported"}, ev)
        m["stereochemistry"]["chirality_status"] = F("relative stereochemistry depicted; absolute configuration not named", "partially_specified", ev)
        m["stereochemistry"]["geometry_status"] = F("alkene geometry not assigned in the source drawing", "unspecified", ev)
        m["stereochemistry"]["representation_status"] = F("depicted stereochemistry not encoded", "not_encoded", ev)
    else:
        m["stereochemistry"]["chirality_status"] = F("stereochemical status not established from name alone", "unresolved", ev)
        m["stereochemistry"]["geometry_status"] = F("geometric stereochemistry not established", "unresolved", ev)
        m["stereochemistry"]["representation_status"] = F("no stereochemical structure encoding asserted", "unresolved", ev)
    m["review_note"] = "graphical/unencoded" if kind in {"substrate", "product"} else "named/unencoded"
    m.pop("review_note")
    return m


def smiles_molecule(mid, label, name, smiles, ev):
    m = source_molecule(mid, label, ev, name=name)
    diag = st.stereo_inspect(smiles)
    why = "Candidate SMILES derived from the explicit source name/drawing; RDKit graph diagnostics do not certify source identity."
    m["representation_type"] = "smiles"
    m["smiles"] = F(label, smiles, ev, method="structure_derivation", derivation=why)
    m["stereochemistry"]["atom_centers"] = F(None, diag["atom_centers"], ev, method="structure_derivation", derivation=why)
    m["stereochemistry"]["double_bonds"] = F(None, diag["double_bonds"], ev, method="structure_derivation", derivation=why)
    m["stereochemistry"]["axial"] = F(None, [], ev, method="structure_derivation", derivation=why)
    m["stereochemistry"]["mixture"] = F("single named reagent/solvent", None, ev, status="unreadable") if False else missing(ev)
    m["stereochemistry"]["assignment"] = F(None, {"tool_diagnostic": "no specified stereogenic elements"}, ev, method="structure_derivation", derivation=why)
    m["stereochemistry"]["chirality_status"] = F(None, "achiral", ev, method="structure_derivation", derivation=why)
    m["stereochemistry"]["geometry_status"] = F(None, "not_applicable", ev, method="structure_derivation", derivation=why)
    m["stereochemistry"]["representation_status"] = F(None, "not_applicable", ev, method="structure_derivation", derivation=why)
    m["stereochemistry"]["isomeric_smiles"] = F(None, diag["canonical_isomeric_smiles"], ev, method="structure_derivation", derivation=why)
    m["stereochemistry"]["identity_key"] = F(None, diag["stereo_graph_key"], ev, method="structure_derivation", derivation=why)
    return m, diag


def template_molecule(mid, label, fixed, variables, ev, kind):
    m = source_molecule(mid, label, ev, kind=kind, name=label)
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": F(fixed, {"description": fixed}, ev),
        "variables": [{
            "label": v,
            "attachment_points": F(f"source-labelled attachment point {v}", {"label": v}, ev),
            "definition": F(desc, desc, ev),
            "allowed_values": F("values shown in Schemes 2/3", "source-scope values only", ev),
            "constraints": F("correlated variables follow the individual source drawings", "do not freely recombine", ev),
        } for v, desc in variables],
        "instance_mapping": [],
        "reason": F("Generic R-group drawing", "Source presents a variable scaffold; it is not a concrete SMILES.", ev),
    }
    return m


def participant(pid, mid, role, ev, stage="reaction", eq=None, molpct=None,
                amount=None, volume=None, ref=None, product_role=None):
    p = st.new_participant(pid, mid, role, ev, stage=stage)
    if product_role:
        p["product_role"] = F(product_role, product_role, ev)
    for key, val in (("equivalents", eq), ("mol_percent", molpct),
                     ("amount", amount), ("volume", volume)):
        if val is not None:
            p["charges"][key] = F(val.get("raw_text", str(val)), val, ev)
    if ref:
        p["equivalence_reference"] = F(ref, ref, ev)
    return p


def conflict_molpct(ev):
    a = F("0.5 mol% on the Scheme 2 reaction arrow and in the 2a prose", {"number": 0.5, "unit": "mol%"}, ev)
    b = F("5 mol% in the Scheme 2 caption/Table 1", {"number": 5, "unit": "mol%"}, ev)
    return F(["0.5 mol%", "5 mol%"], None, ev, status="conflict", alternatives=[a, b])


def reaction(rid, ev, rtype, categories, location, name, transform):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = rtype
    r["experiment_categories"] = F(categories, categories, ev)
    r["location"] = F(location, {"source_location": location}, ev)
    r["reaction_name"] = F(name, name, ev)
    r["transformation_type"] = F(transform, transform, ev)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    return r


def set_standard_conditions(r, ev, substrate_mid, substrate_pid, reagent_mid,
                            catalyst_mid="mol_fac_ir", catalyst_loading=None,
                            solvent_mid="mol_dioxane", atmosphere="Ar",
                            time=True, scale=True):
    r["participants"].append(participant("p_substrate", substrate_mid, "reactant", ev,
                                         eq={"number": 1, "unit": "equiv", "raw_text": "1.0 equiv."},
                                         amount={"number": 0.2, "unit": "mmol", "raw_text": "0.2 mmol"},
                                         ref=substrate_mid))
    r["participants"].append(participant("p_rf", reagent_mid, "reagent", ev,
                                         eq={"number": 1.2, "unit": "equiv", "raw_text": "1.2 equiv."},
                                         amount={"number": 0.24, "unit": "mmol", "raw_text": "0.24 mmol"},
                                         ref=substrate_mid))
    if catalyst_mid:
        cat = participant("p_pc", catalyst_mid, "catalyst", ev, ref=substrate_mid)
        if catalyst_loading == "conflict":
            cat["charges"]["mol_percent"] = conflict_molpct(ev)
        else:
            value = 5 if catalyst_loading is None else catalyst_loading
            cat["charges"]["mol_percent"] = F(f"{value} mol%", {"number": value, "unit": "mol%"}, ev)
        r["participants"].append(cat)
    if solvent_mid:
        r["participants"].append(participant("p_solvent", solvent_mid, "solvent", ev,
                                             volume={"number": 10, "unit": "mL", "raw_text": "10 mL"}))
        r["solvent_conditions"]["composition"] = F("dioxane", ["dioxane"], ev)
        r["solvent_conditions"]["total_volume"] = F("10 mL", {"number": 10, "unit": "mL"}, ev)
    if scale:
        r["reaction_scale"] = F("1 (0.2 mmol)", {"reference_molecule_id": substrate_mid,
                                                  "amount": {"number": 0.2, "unit": "mmol"}}, ev)
    step = st.new_step("reaction", ev)
    step["participant_ids"] = [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"]
    step["operation"] = F("irradiation with a 30 W blue LED", "visible-light irradiation", ev)
    step["temperature"] = F("ambient temperature", "ambient temperature", ev)
    if time:
        step["time"] = F("12 hours", {"number": 12, "unit": "h"}, ev)
    step["atmosphere"] = F(atmosphere, atmosphere, ev)
    r["steps"] = [step]
    r["special_conditions"]["light_source"] = F("blue LED", "blue LED", ev)
    r["special_conditions"]["power"] = F("30 W", {"number": 30, "unit": "W"}, ev)
    r["procedure"]["original_text"] = F("standard conditions shown in the source scheme/table", "source scheme/table conditions", ev)


def outcome(oid, kind, targets, ev, raw, norm, basis=None, yield_type=None,
            product_role="desired product", scope=None):
    o = st.new_outcome(oid, kind, targets, ev)
    o["value"] = F(raw, norm, ev)
    if product_role is not None:
        o["product_role"] = F(product_role, product_role, ev)
    if yield_type:
        o["yield_type"] = F(yield_type, yield_type, ev)
    if basis:
        o["measurement_basis"] = F(basis, basis, ev)
    if scope:
        o["scope"] = F(scope, scope, ev)
    return o


def set_classification(r, ev, mode="cascade"):
    c = r["classification"]
    c["mapped_reaction_smiles"] = missing(ev)
    c["model_votes"] = []
    if mode == "none":
        c["bond_changes"] = []
        c["labels"] = inferred([], ev, "No observed synthetic transformation; no formed-bond labels assigned.")
        c["coverage"] = inferred("not_applicable", ev, "No observed synthetic transformation.")
        c["limitations"] = ["Negative/control observation; intended product bonds were not classified."]
        return
    if mode == "unresolved":
        c["bond_changes"] = []
        c["labels"] = inferred([], ev, "Experiment does not itself provide a net preparative transformation to classify.")
        c["coverage"] = inferred("unresolved", ev, "Net bond inventory is not assigned for this measurement-only experiment.")
        c["limitations"] = ["Measurement-only mechanistic experiment; no atom-mapped net reaction supplied."]
        return
    if mode == "tempo":
        changes = [{"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}],
                    "before_order": 0, "after_order": 1, "change_type": "formed"}]
        labels = [{"category_number": 3, "code": "C-O", "supporting_change_ids": ["b1"]}]
        why = "Partial structure derivation from the TEMPO-perfluoroalkyl adduct drawing; no atom map is asserted."
    else:
        changes = [
            {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "I"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b4", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 2, "after_order": 3, "change_type": "order_increased"},
        ]
        labels = [
            {"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1", "b2"]},
            {"category_number": 10, "code": "C-I", "supporting_change_ids": ["b3"]},
            {"category_number": 21, "code": "TB-C-C", "supporting_change_ids": ["b4"]},
        ]
        why = "Partial net-bond interpretation of the source reaction drawing: fluoroalkyl addition/cyclization, iodination, and alkyne formation; exact atoms are not mapped."
    c["bond_changes"] = [F(None, x, ev, method="structure_derivation", derivation=why) for x in changes]
    c["labels"] = F(None, labels, ev, method="structure_derivation", derivation=why)
    c["coverage"] = inferred("partial", ev, "Graphical structures are not fully encoded or atom mapped; only clearly depicted net bond classes are retained.")
    c["limitations"] = ["No atom-mapped reaction SMILES.", "Complex graphical structures remain unencoded.",
                          "Partial labels do not claim a complete elementary mechanism or complete hydrogen inventory."]


bundle = st.empty(copy.deepcopy(DATASET["documents"]))
stereo_diagnostics = {}

# Paper-local substrates and products.
for letter in "abcdefghijklmnopqr":
    ev = E[2] + E[3] if letter == "a" else E[3]
    bundle["molecules"].append(source_molecule(f"mol_1{letter}", f"1{letter}", ev, kind="substrate",
                                                name=f"ene-vinylidenecyclopropane 1{letter}"))
    bundle["molecules"].append(source_molecule(f"mol_2{letter}", f"2{letter}", ev, kind="product",
                                                name=f"fluorinated cyclization product 2{letter}"))
for letter in "abcdefgh":
    ev = E[3] + E[4] if letter == "h" else E[3]
    bundle["molecules"].append(source_molecule(f"mol_3{letter}", f"3{letter}", ev, kind="product",
                                                name=f"perfluoroalkylated cyclization product 3{letter}"))
bundle["molecules"].append(source_molecule("mol_4", "4", E[3] + E[4], kind="named",
                                            name="TEMPO-incorporated perfluoroalkyl radical adduct 4"))

# Generic source drawings used by overview/scope schemes.
bundle["molecules"].append(template_molecule(
    "mol_generic_1", "generic ene-VDCP 1", "ene-vinylidenecyclopropane scaffold with tether X and variables n, R1, R2",
    [("X", "NTs, NBs, O, or carbon tether as drawn"), ("n", "1 or 2"), ("R1", "H or aryl"), ("R2", "source scope substituent")], E[1] + E[3], "substrate"))
bundle["molecules"].append(template_molecule(
    "mol_generic_2", "generic product 2", "iododifluoroacetate-substituted cyclic alkynyl product scaffold",
    [("X", "inherited from substrate"), ("n", "1 or 2"), ("R1", "inherited"), ("R2", "inherited")], E[1] + E[3], "product"))
bundle["molecules"].append(template_molecule(
    "mol_generic_3", "generic product 3", "iodoperfluorobutyl-substituted cyclic alkynyl product scaffold",
    [("X", "inherited from substrate"), ("n", "1 or 2"), ("R1", "inherited"), ("R2", "inherited")], E[3] + E[4], "product"))

# Concrete small molecules for which a candidate SMILES and standard stereo diagnostic are practical.
small = [
    ("mol_icf2co2et", "ICF2CO2Et", "ethyl iododifluoroacetate", "CCOC(=O)C(F)(F)I", E[1] + E[2] + E[3] + E[4]),
    ("mol_brcf2co2et", "BrCF2CO2Et", "ethyl bromodifluoroacetate", "CCOC(=O)C(F)(F)Br", E[2]),
    ("mol_ic4f9", "ICF2CF2CF2CF3", "1-iodoperfluorobutane", "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)I", E[3] + E[4]),
    ("mol_brc6f13", "BrCF2(CF2)4CF3", "1-bromoperfluorohexane", "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)Br", E[2]),
    ("mol_dioxane", "dioxane", "1,4-dioxane", "O1CCOCC1", E[2] + E[3]),
    ("mol_mecn", "MeCN", "acetonitrile", "CC#N", E[2]),
    ("mol_dcm", "DCM", "dichloromethane", "ClCCl", E[2]),
    ("mol_dmf", "DMF", "N,N-dimethylformamide", "CN(C)C=O", E[2]),
    ("mol_fluorobenzene", "fluorobenzene", "fluorobenzene", "Fc1ccccc1", E[2]),
]
for args in small:
    m, diag = smiles_molecule(*args)
    bundle["molecules"].append(m)
    stereo_diagnostics[args[0]] = diag

for mid, label, name, ev in [
    ("mol_fac_ir", "fac-Ir(ppy)3", "fac-tris(2-phenylpyridine)iridium(III)", E[1] + E[2] + E[3] + E[4]),
    ("mol_pc2", "Ir[dF(CF3)(ppy)]2[dCF3(bpy)]PF6", "Ir[dF(CF3)(ppy)]2[dCF3(bpy)]PF6", E[2]),
    ("mol_pc3", "Ir[dF(CF3)(ppy)]2(dtbbpy)PF6", "Ir[dF(CF3)(ppy)]2(dtbbpy)PF6", E[2]),
    ("mol_pc4", "Ir(ppy)2(dtbbpy)PF6", "Ir(ppy)2(dtbbpy)PF6", E[2]),
    ("mol_pc5", "Ir(dFCF3ppy)2(bpy)PF6", "Ir(dFCF3ppy)2(bpy)PF6", E[2]),
    ("mol_tempo", "TEMPO", "2,2,6,6-tetramethyl-1-piperidinyloxy", E[3] + E[4]),
]:
    bundle["molecules"].append(source_molecule(mid, label, ev, name=name))

# Generic difluoroacetate scheme (Scheme 1c + Scheme 2 overview).
r = reaction("generic_difluoroacetate_cascade", E[1] + E[3], "generic_scheme", ["other"],
             "Scheme 1c and Scheme 2 generic arrow", "fluoroalkylation and cycloisomerization cascade",
             "tandem ring opening/iodofluoroalkylation/cycloisomerization")
r["participants"] = [participant("p_substrate", "mol_generic_1", "reactant", E[1] + E[3]),
                     participant("p_rf", "mol_icf2co2et", "reagent", E[1] + E[3]),
                     participant("p_product", "mol_generic_2", "product", E[1] + E[3], product_role="desired product")]
cat = participant("p_pc", "mol_fac_ir", "catalyst", E[1] + E[3], ref="mol_generic_1")
cat["charges"]["mol_percent"] = conflict_molpct(E[3])
r["participants"].append(cat)
r["notes"] = [F("Scheme 2 arrow: 0.5 mol%; Scheme 2 caption: 5 mol%.",
                  "Photocatalyst loading conflict retained without choosing a value.", E[3])]
r["review"]["unresolved"] = ["Scheme 2 gives conflicting fac-Ir(ppy)3 loadings (0.5 mol% arrow versus 5 mol% caption).",
                               "Generic R-group structures are not concrete molecules."]
set_classification(r, E[1] + E[3])
bundle["reactions"].append(r)

# Generic perfluorobutyl scope scheme.
r = reaction("generic_perfluorobutyl_cascade", E[3], "generic_scheme", ["other"], "Scheme 3 generic arrow",
             "perfluoroalkylation and cycloisomerization cascade", "iodoperfluoroalkylation/cycloisomerization")
r["participants"] = [participant("p_substrate", "mol_generic_1", "reactant", E[3]),
                     participant("p_rf", "mol_ic4f9", "reagent", E[3]),
                     participant("p_pc", "mol_fac_ir", "catalyst", E[3], molpct={"number": 5, "unit": "mol%", "raw_text": "5 mol%"}, ref="mol_generic_1"),
                     participant("p_product", "mol_generic_3", "product", E[3], product_role="desired product")]
r["review"]["unresolved"] = ["Generic R-group structures are not concrete molecules."]
set_classification(r, E[3])
bundle["reactions"].append(r)

# Initial catalyst-free screening statement (Table S1 is not supplied).
r = reaction("initial_catalyst_free_80_nmr", E[2], "specific_experiment", ["optimization"],
             "page 2 prose; details referred to unavailable Table S1", "catalyst-free visible-light fluoroalkylation/cyclization",
             "fluoroalkylation/cycloisomerization")
r["participants"] = [participant("p_substrate", "mol_1a", "reactant", E[2]),
                     participant("p_rf", "mol_icf2co2et", "reagent", E[2]),
                     participant("p_product", "mol_2a", "product", E[2], product_role="desired product")]
step = st.new_step("reaction", E[2]); step["participant_ids"] = ["p_substrate", "p_rf", "p_product"]
step["operation"] = F("blue LED irradiation without adding a photocatalyst", "visible-light irradiation; no photocatalyst", E[2])
step["temperature"] = F("room temperature", "room temperature", E[2]); step["atmosphere"] = F("argon atmosphere", "Ar", E[2])
r["steps"] = [step]
r["special_conditions"]["light_source"] = F("blue LED", "blue LED", E[2])
r["outcomes"] = [outcome("o_nmr_yield", "yield", ["mol_2a"], E[2], "80% NMR yield", {"number": 80, "unit": "%"}, basis="NMR", yield_type="NMR yield")]
r["review"]["unresolved"] = ["Solvent, lamp wattage, reaction time, exact charges, and NMR internal standard are in unavailable Table S1/ESI.",
                               "Complex mixture prevented isolation of pure 2a; composition is not reported."]
set_classification(r, E[2])
bundle["reactions"].append(r)

# Table 1 optimization rows; row 1 is the same 2a experiment repeated in Scheme 2.
rows = {
    1: ("mol_icf2co2et", "mol_fac_ir", "mol_dioxane", "Ar", "blue", 92),
    2: ("mol_icf2co2et", "mol_pc2", "mol_dioxane", "Ar", "blue", 60),
    3: ("mol_icf2co2et", "mol_pc3", "mol_dioxane", "Ar", "blue", 85),
    4: ("mol_icf2co2et", "mol_pc4", "mol_dioxane", "Ar", "blue", 35),
    5: ("mol_icf2co2et", "mol_pc5", "mol_dioxane", "Ar", "blue", 51),
    6: ("mol_icf2co2et", "mol_fac_ir", "mol_mecn", "Ar", "blue", 48),
    7: ("mol_icf2co2et", "mol_fac_ir", "mol_dcm", "Ar", "blue", 54),
    8: ("mol_icf2co2et", "mol_fac_ir", "mol_dmf", "Ar", "blue", "Trace"),
    9: ("mol_brcf2co2et", "mol_fac_ir", "mol_dioxane", "Ar", "blue", 65),
    10: ("mol_brc6f13", "mol_fac_ir", "mol_dioxane", "Ar", "blue", "Nd"),
    11: ("mol_icf2co2et", "mol_fac_ir", "mol_dioxane", "Ar", "none", "Nd"),
    12: ("mol_icf2co2et", "mol_fac_ir", "mol_dioxane", "Air", "blue", "Nd"),
}
table_ids = []
for n, (rf, pc, solv, atm, light, result) in rows.items():
    rid = "table1_entry01_scope_2a" if n == 1 else f"table1_entry{n:02d}"
    ev = E[2] + E[3] if n == 1 else E[2]
    r = reaction(rid, ev, "specific_experiment", ["optimization"] + (["substrate_scope"] if n == 1 else []),
                 f"Table 1 entry {n}" + ("; Scheme 2 product 2a" if n == 1 else ""),
                 "visible-light fluoroalkylation/cycloisomerization", "fluoroalkylation/cycloisomerization")
    r["general_scheme_ids"] = [f"{DOC}:generic_difluoroacetate_cascade"]
    product_observed = result != "Nd"
    set_standard_conditions(r, ev, "mol_1a", "p_substrate", rf, catalyst_mid=pc,
                            catalyst_loading="conflict" if n == 1 else 5,
                            solvent_mid=solv, atmosphere=atm)
    if light == "none":
        r["steps"][0]["operation"] = F("w/o light", "no light", E[2])
        r["special_conditions"]["light_source"] = F("w/o light", "none", E[2])
        r["special_conditions"]["power"] = missing(E[2])
    if product_observed:
        r["participants"].append(participant("p_product", "mol_2a", "product", ev, product_role="desired product"))
    r["participants"].append(participant("p_internal_standard", "mol_fluorobenzene", "reagent", E[2], stage="characterization"))
    if isinstance(result, int):
        r["outcomes"] = [outcome("o_nmr_yield", "yield", ["mol_2a"], E[2], f"{result}%", {"number": result, "unit": "%"},
                                  basis="19F NMR using fluorobenzene as internal standard", yield_type="19F NMR yield")]
        if n == 1:
            r["outcomes"].append(outcome("o_isolated_yield", "yield", ["mol_2a"], ev, "(83) isolated", {"number": 83, "unit": "%"},
                                                 basis="isolated", yield_type="isolated yield"))
            r["outcomes"].append(outcome("o_dr", "selectivity", ["mol_2a"], E[2] + E[3], "single diastereoisomer", "single diastereoisomer",
                                                 scope={"metric": "diastereoselectivity", "configuration": "absolute assignment unreported"}))
    elif result == "Trace":
        r["outcomes"] = [outcome("o_trace", "trace", ["mol_2a"], E[2], "Trace", "trace",
                                  basis="19F NMR using fluorobenzene as internal standard")]
    else:
        r["outcomes"] = [outcome("o_no_product", "no_reaction", [], E[2], "Nd", "not detected",
                                  basis="19F NMR using fluorobenzene as internal standard", product_role=None)]
    if n == 1:
        r["notes"] = [F("0.5 mol% in prose/Scheme 2 arrow; 5 mol% in Table 1/Scheme 2 caption.",
                          "Photocatalyst-loading conflict retained.", ev)]
        r["review"]["unresolved"].append("fac-Ir(ppy)3 loading conflicts: 0.5 mol% in prose/Scheme 2 arrow versus 5 mol% in Table 1/Scheme 2 caption.")
    if n in (9, 10):
        r["notes"].append(F("The prose says the alternative bromide reacts 'with 2a', whereas Table 1 depicts substrate 1a and reports the 2a column.",
                              "Source wording discrepancy retained; Table 1 substrate identity used for this row.", E[2]))
        r["review"]["unresolved"].append("Page-2 prose names 2a as reaction partner for entries 9/10, conflicting with the Table 1 substrate drawing 1a.")
    r["review"]["unresolved"].append("Workup and purification details are unavailable because no SI was provided.")
    set_classification(r, ev, "cascade" if product_observed else "none")
    bundle["reactions"].append(r); table_ids.append(r["reaction_id"])

# Scheme 2 scope entries 2b-r (2a is merged above).
scope2_yields = dict(b=78, c=73, d=70, e=67, f=60, g=57, h=80, i=85, j=73,
                     k=50, l=54, m=43, n=76, o=57, p=40, q="Trace", r="Trace")
scope2_ids = [f"{DOC}:table1_entry01_scope_2a"]
for letter, y in scope2_yields.items():
    ev = E[3]
    rid = f"scheme2_scope_2{letter}"
    r = reaction(rid, ev, "specific_experiment", ["substrate_scope"], f"Scheme 2, 1{letter} to 2{letter}",
                 "visible-light fluoroalkylation/cycloisomerization", "fluoroalkylation/cycloisomerization")
    r["general_scheme_ids"] = [f"{DOC}:generic_difluoroacetate_cascade"]
    set_standard_conditions(r, ev, f"mol_1{letter}", "p_substrate", "mol_icf2co2et", catalyst_loading="conflict")
    r["participants"].append(participant("p_product", f"mol_2{letter}", "product", ev, product_role="desired product"))
    if isinstance(y, int):
        r["outcomes"] = [outcome("o_isolated_yield", "yield", [f"mol_2{letter}"], ev, f"{y}%", {"number": y, "unit": "%"},
                                  basis="isolated", yield_type="isolated yield")]
        if letter in "bcdefg":
            r["outcomes"].append(outcome("o_dr", "selectivity", [f"mol_2{letter}"], ev, "single diastereoisomer", "single diastereoisomer",
                                                 scope={"metric": "diastereoselectivity", "configuration": "absolute assignment unreported"}))
    else:
        r["outcomes"] = [outcome("o_trace", "trace", [f"mol_2{letter}"], ev, "trace", "trace", basis="isolated-product scope; trace method not separately reported")]
    r["notes"] = [F("Scheme 2 arrow: fac-Ir(ppy)3 (0.5 mol%); caption: fac-Ir(ppy)3 (5 mol%).",
                      "Photocatalyst-loading conflict retained.", ev)]
    r["review"]["unresolved"] = ["fac-Ir(ppy)3 loading conflicts within Scheme 2 (0.5 mol% arrow versus 5 mol% caption).",
                                   "Exact complex structure and absolute stereochemical assignment remain unencoded.",
                                   "Workup/purification are unavailable without SI."]
    set_classification(r, ev)
    bundle["reactions"].append(r); scope2_ids.append(r["reaction_id"])

# Scheme 3 perfluorobutyl scope.
scope3 = {"a": ("a", 70), "b": ("n", 65), "c": ("h", 92), "d": ("i", 93),
          "e": ("j", 82), "f": ("k", 48), "g": ("l", 45), "h": ("q", 23)}
scope3_ids = []
for product_letter, (substrate_letter, y) in scope3.items():
    ev = E[3]
    rid = f"scheme3_scope_3{product_letter}"
    r = reaction(rid, ev, "specific_experiment", ["substrate_scope"],
                 f"Scheme 3, structure corresponding to 1{substrate_letter} to 3{product_letter}",
                 "visible-light perfluoroalkylation/cycloisomerization", "perfluoroalkylation/cycloisomerization")
    r["general_scheme_ids"] = [f"{DOC}:generic_perfluorobutyl_cascade"]
    set_standard_conditions(r, ev, f"mol_1{substrate_letter}", "p_substrate", "mol_ic4f9", catalyst_loading=5)
    r["participants"].append(participant("p_product", f"mol_3{product_letter}", "product", ev, product_role="desired product"))
    r["outcomes"] = [outcome("o_isolated_yield", "yield", [f"mol_3{product_letter}"], ev, f"{y}%", {"number": y, "unit": "%"},
                              basis="isolated", yield_type="isolated yield")]
    r["review"]["unresolved"] = ["Substrate label is assigned by matching the Scheme 3 drawing to Scheme 2; Scheme 3 does not print the individual 1x label.",
                                   "Exact complex structure and absolute stereochemical assignment remain unencoded.",
                                   "Workup/purification are unavailable without SI."]
    set_classification(r, ev)
    bundle["reactions"].append(r); scope3_ids.append(r["reaction_id"])

# Fluorescence quenching experiment.
r = reaction("fluorescence_quenching", E[3] + E[4], "specific_experiment", ["mechanistic_study"],
             "Fig. 1 and page 4 prose", "fluorescence quenching experiment", "photophysical quenching measurement")
r["participants"] = [participant("p_pc", "mol_fac_ir", "catalyst", E[3] + E[4]),
                     participant("p_quencher", "mol_icf2co2et", "reagent", E[3] + E[4])]
s = st.new_step("fluorescence_measurement", E[3] + E[4], stage="characterization")
s["participant_ids"] = ["p_pc", "p_quencher"]
s["operation"] = F("Fluorescence quenching experiments; Q = 0, 0.2, 0.4, 0.6, 0.8 mM", "fluorescence quenching series", E[3])
r["steps"] = [s]
r["outcomes"] = [outcome("o_quenching", "other", [], E[3] + E[4],
                          "ICF2CO2Et was an effective quencher for the excited-state photocatalyst",
                          {"observation": "effective quenching", "quencher_concentrations_mM": [0, 0.2, 0.4, 0.6, 0.8]},
                          basis="fluorescence intensity/Stern-Volmer plot", product_role=None)]
r["review"]["unresolved"] = ["Excitation wavelength, solution composition, photocatalyst concentration, and full protocol are only in unavailable ESI."]
set_classification(r, E[3] + E[4], "unresolved")
bundle["reactions"].append(r)

# TEMPO radical trapping experiment.
r = reaction("tempo_radical_trap", E[3] + E[4], "specific_experiment", ["control", "mechanistic_study"],
             "Scheme 4 and page 4 prose", "radical trapping study", "TEMPO trapping of perfluoroalkyl radical")
set_standard_conditions(r, E[3] + E[4], "mol_1h", "p_substrate", "mol_ic4f9", catalyst_loading=5)
r["participants"].append(participant("p_tempo", "mol_tempo", "reagent", E[3] + E[4],
                                     eq={"number": 2, "unit": "equiv", "raw_text": "2.0 equiv."}, ref="mol_1h"))
r["participants"].append(participant("p_expected", "mol_3h", "product", E[3] + E[4], product_role="desired product (not formed)"))
r["participants"].append(participant("p_adduct", "mol_4", "product", E[3] + E[4], product_role="radical-trap adduct"))
r["outcomes"] = [outcome("o_3h", "yield", ["mol_3h"], E[3], "3h, 0%", {"number": 0, "unit": "%"}, yield_type="reported yield"),
                 outcome("o_4", "yield", ["mol_4"], E[3] + E[4], "4, 20%; isolated in 20% yield", {"number": 20, "unit": "%"},
                         basis="isolated", yield_type="isolated yield", product_role="TEMPO adduct")]
r["review"]["unresolved"] = ["Measurement basis for the printed 0% of 3h is not separately reported.",
                               "Exact structure of adduct 4 remains graphical/unencoded."]
set_classification(r, E[3] + E[4], "tempo")
bundle["reactions"].append(r)

# Light on/off experiment, distinct from the Scheme 3 preparative 3h run.
r = reaction("light_on_off_1h_to_3h", E[4], "specific_experiment", ["control", "mechanistic_study"],
             "Scheme 5", "light on/off experiment", "perfluoroalkylation/cycloisomerization under intermittent irradiation")
set_standard_conditions(r, E[4], "mol_1h", "p_substrate", "mol_ic4f9", catalyst_loading=5)
r["participants"].append(participant("p_product", "mol_3h", "product", E[4], product_role="desired product"))
r["steps"][0]["operation"] = F("alternating light on/off intervals shown in Scheme 5", "intermittent blue-light irradiation", E[4])
r["steps"][0]["time"] = F("time-course plotted from 0 to 140 min", {"minimum": 0, "maximum": 140, "unit": "min", "raw_text": "0–140 min"}, E[4])
r["outcomes"] = [outcome("o_dark_progress", "other", ["mol_3h"], E[4],
                          "the reaction also proceeded when the light was turned off",
                          "product formation continued during dark intervals", basis="yield-time plot",
                          product_role="desired product", scope={"metric": "light on/off time course"})]
r["review"]["unresolved"] = ["Individual time-point yields are graphical and were not digitized as exact values.",
                               "Full on/off schedule and analytical method are in unavailable ESI."]
set_classification(r, E[4])
bundle["reactions"].append(r)

# Quantum-yield measurement; exact substrate identity is not stated in supplied main-text passage.
r = reaction("quantum_yield_phi23", E[4], "specific_experiment", ["mechanistic_study"],
             "page 4 prose; details referred to ESI", "quantum yield measurement", "photochemical quantum-yield measurement")
r["participants"] = [participant("p_substrate", "mol_generic_1", "reactant", E[4]),
                     participant("p_rf", "mol_icf2co2et", "reagent", E[4]),
                     participant("p_pc", "mol_fac_ir", "catalyst", E[4])]
s = st.new_step("quantum_yield_measurement", E[4]); s["participant_ids"] = ["p_substrate", "p_rf", "p_pc"]
s["operation"] = F("quantum yield was measured", "quantum-yield measurement", E[4]); r["steps"] = [s]
r["outcomes"] = [outcome("o_phi", "other", [], E[4], "Φ = 23", {"metric": "quantum yield", "symbol": "Phi", "number": 23},
                          basis="method/details in unavailable ESI", product_role=None)]
r["review"]["unresolved"] = ["Exact substrate/product identity and quantum-yield protocol are not stated in the supplied main-text passage; ESI was not provided."]
set_classification(r, E[4], "unresolved")
bundle["reactions"].append(r)

# Proposed arrow-pushing mechanism, not an observed experiment.
r = reaction("scheme6_proposed_mechanism", E[4], "proposed_mechanism", ["other"], "Scheme 6",
             "proposed mechanism", "SET/radical-chain fluoroalkylation, cyclization, ring opening, and iodine transfer")
r["participants"] = [participant("p_substrate", "mol_generic_1", "reactant", E[4]),
                     participant("p_rf", "mol_icf2co2et", "reagent", E[4]),
                     participant("p_pc", "mol_fac_ir", "catalyst", E[4]),
                     participant("p_product", "mol_generic_2", "product", E[4], product_role="proposed final product 2")]
r["notes"] = [F("Scheme 6 caption says 'production of 3a', while the drawing/text use ICF2CO2Et and final product 2.",
                  "Caption/content target discrepancy retained without correction.", E[4]),
              F("Paths I and II are proposed by the authors, not directly observed elementary steps.",
                  "proposed mechanism only", E[4])]
r["review"]["unresolved"] = ["Scheme 6 caption identifies 3a, but scheme contents and prose identify difluoroacetate product 2.",
                               "Intermediates A-D and elementary atom transfers are proposed, not independently verified here."]
set_classification(r, E[4])
bundle["reactions"].append(r)

bundle["source_context"] = [
    {"context_id": f"{DOC}:scheme1_prior_work", "kind": "prior_work_summary",
     "evidence_ids": E[1], "note": "Scheme 1a-b summarizes cited prior work and is not imported as this paper's experiments."},
    {"context_id": f"{DOC}:xray_2k", "kind": "structure_characterization",
     "evidence_ids": E[2] + E[3], "note": "The structure of 2k is stated to be determined by X-ray crystallography (CCDC 2047670); no separate synthesis record is created."},
    {"context_id": f"{DOC}:mechanistic_discussion", "kind": "author_interpretation",
     "evidence_ids": E[4], "note": "The authors infer a radical-chain contribution from light on/off behavior and Phi = 23."},
    {"context_id": f"{DOC}:bibliography_pages_5_6", "kind": "bibliography_only",
     "evidence_ids": E[5] + E[6], "note": "Pages 5-6 are references only; cited-paper reactions are not current-paper experiments."},
    {"context_id": f"{DOC}:characterization_inventory", "kind": "coverage_note",
     "evidence_ids": E[1] + E[2] + E[3] + E[4] + E[5] + E[6],
     "note": "No standalone mass-spectrum peak table, molecular-weight table, or fragment-relative-abundance dataset occurs on pages 1-6."},
]

# Persist raw intermediate inventories separately from the formal submission.
raw_inventory = {
    "document_id": DOC,
    "source_title": DATASET["documents"][0]["title"],
    "page_scope": [1, 2, 3, 4, 5, 6],
    "table1_record_ids": table_ids,
    "scheme2_record_ids": scope2_ids,
    "scheme3_record_ids": scope3_ids,
    "mechanistic_experiment_ids": [f"{DOC}:fluorescence_quenching", f"{DOC}:tempo_radical_trap",
                                     f"{DOC}:light_on_off_1h_to_3h", f"{DOC}:quantum_yield_phi23"],
    "proposed_mechanism_ids": [f"{DOC}:scheme6_proposed_mechanism"],
    "visual_findings": {
        "page1": "Title/abstract and Scheme 1; panels a-b are prior work, panel c is the current generic transformation.",
        "page2": "Table 1 entries 1-12 plus prose describing catalyst-free screening and source wording/loading discrepancies.",
        "page3": "Scheme 2 entries 2a-r, Scheme 3 entries 3a-h, Fig. 1 fluorescence quenching, and Scheme 4 TEMPO trap.",
        "page4": "Scheme 5 light on/off, quantum-yield statement, Scheme 6 proposed mechanism, and continuation of TEMPO/quenching prose.",
        "page5": "Bibliography only.", "page6": "Bibliography only."
    },
    "semantic_accuracy_verified": False,
}
st.write(RUN / "raw_model/source_inventory.json", raw_inventory)
st.write(RUN / "raw_model/stereo_diagnostics.json", stereo_diagnostics)
st.write(RUN / "raw_model/source_mentions.json",
         st.source_mentions(DATASET, ["1a", "2a", "1h", "3h", "2k"], context=220))
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title, abstract, and ESI footnote", "finding": "Current-paper identity and summary; ESI is said to exist but no SI file was supplied or associated.", "evidence_ids": E[1]},
            {"region": "Scheme 1a-b", "finding": "Prior-work reaction classes; retained as source context, not current-paper experiments.", "evidence_ids": E[1]},
            {"region": "Scheme 1c", "finding": "Current-work generic ene-VDCP fluoroalkylation/cycloisomerization arrow.", "evidence_ids": E[1]},
        ], "record_ids": [f"{DOC}:generic_difluoroacetate_cascade"],
         "unresolved": ["Scheme 1c is generic and does not enumerate concrete R-group identities or conditions.", "ESI association remains unestablished."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "left-column screening prose", "finding": "Catalyst-free 80% NMR screening result and optimized 2a prose; Table S1/ESI details unavailable.", "evidence_ids": E[2]},
            {"region": "Table 1 reaction arrow and entries 1-12", "finding": "Twelve independent optimization/control rows; row 1 also reports 83% isolated yield and is the same 2a preparation repeated in Scheme 2.", "evidence_ids": E[2]},
            {"region": "Table 1 footnotes", "finding": "1a 0.2 mmol/1.0 equiv, fluoroalkyl reagent 0.24 mmol/1.2 equiv, PC 5 mol%, solvent 10 mL, ambient temperature, 30 W blue LED, 12 h; 19F NMR/fluorobenzene basis and isolated-yield footnote.", "evidence_ids": E[2]},
            {"region": "lower-right scope discussion", "finding": "Narrative corroboration for Scheme 2 yields/no-product outcomes and X-ray assignment of 2k; no extra experiments counted.", "evidence_ids": E[2]},
        ], "record_ids": [f"{DOC}:initial_catalyst_free_80_nmr"] + table_ids,
         "unresolved": ["Prose/Scheme 2 report 0.5 mol% fac-Ir(ppy)3 while Table 1/Scheme 2 caption report 5 mol%; conflict retained.", "Entries 9/10 prose says reaction with 2a while the table depicts substrate 1a.", "No SI procedures, workup, or purification were supplied."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Scheme 2 generic arrow, products 2a-r, caption", "finding": "Eighteen difluoroacetate scope entries; 2a is merged with Table 1 entry 1, 2q/2r are trace outcomes.", "evidence_ids": E[3]},
            {"region": "Scheme 3 generic arrow, products 3a-h, caption", "finding": "Eight isolated-yield perfluorobutyl scope experiments.", "evidence_ids": E[3]},
            {"region": "Fig. 1", "finding": "Fluorescence-quenching mechanistic experiment with 0-0.8 mM quencher series and Stern-Volmer plot.", "evidence_ids": E[3]},
            {"region": "Scheme 4", "finding": "TEMPO (2.0 equiv) radical-trap control: 3h 0%, adduct 4 20%; prose continues on page 4.", "evidence_ids": E[3]},
            {"region": "X-ray inset for 2k", "finding": "Structure-characterization evidence only; no additional synthetic experiment.", "evidence_ids": E[3]},
        ], "record_ids": [f"{DOC}:generic_difluoroacetate_cascade", f"{DOC}:generic_perfluorobutyl_cascade"] + scope2_ids + scope3_ids + [f"{DOC}:fluorescence_quenching", f"{DOC}:tempo_radical_trap"],
         "unresolved": ["Scheme 2 internally conflicts on catalyst loading (0.5 versus 5 mol%).", "Complex substrate/product drawings and absolute stereochemical identities are not encoded as exact SMILES.", "Scheme 3 individual substrate labels are assigned by graphical correspondence and remain review items."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "Scheme 5 and associated prose", "finding": "Distinct light on/off mechanistic time-course using 1h; product formation continues during dark intervals.", "evidence_ids": E[4]},
            {"region": "TEMPO/quenching prose continuation", "finding": "Confirms complete inhibition by TEMPO, isolated 20% adduct 4, and effective quenching by ICF2CO2Et.", "evidence_ids": E[4]},
            {"region": "quantum-yield statement", "finding": "Measured Phi = 23; exact protocol and reaction identity referred to unavailable ESI.", "evidence_ids": E[4]},
            {"region": "Scheme 6 and mechanism prose", "finding": "Proposed SET/radical-chain mechanism, not a directly observed experiment; caption says 3a while contents/text show ICF2CO2Et/product 2.", "evidence_ids": E[4]},
            {"region": "summary and start of references", "finding": "Author summary and bibliographic text; no additional experiment.", "evidence_ids": E[4]},
        ], "record_ids": [f"{DOC}:fluorescence_quenching", f"{DOC}:tempo_radical_trap", f"{DOC}:light_on_off_1h_to_3h", f"{DOC}:quantum_yield_phi23", f"{DOC}:scheme6_proposed_mechanism"],
         "unresolved": ["Light on/off graph points were not digitized as exact yields.", "Quantum-yield experiment identity/protocol unavailable without ESI.", "Scheme 6 caption/content target discrepancy remains unresolved."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "both columns", "finding": "References 1-12 only; cited reactions are prior literature and were not imported as current-paper experiments.", "evidence_ids": E[5]},
        ], "record_ids": [], "unresolved": ["Cited-paper experimental details cannot be extracted from bibliography entries alone."]},
        {"page": 6, "state": "attempted", "regions": [
            {"region": "top of both columns", "finding": "Continuation of references 12-15; remainder blank apart from footer.", "evidence_ids": E[6]},
        ], "record_ids": [], "unresolved": ["Cited-paper experimental details cannot be extracted from bibliography entries alone."]},
    ],
}
st.write(RUN / "coverage/pages_0001_0006.json", coverage)
st.write(RUN / "submitted.json", bundle)
print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
                  "source_context": len(bundle["source_context"])}, ensure_ascii=False))
