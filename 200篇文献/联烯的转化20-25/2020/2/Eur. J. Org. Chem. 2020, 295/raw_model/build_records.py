#!/usr/bin/env python3
"""Build page-1--6 machine candidates for doc_88f54768a8114d5d.

This is a source-extraction artifact, not a source-review certificate.  It uses
the prepared evidence registry without adding or altering evidence.
"""
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st
import reaction_tools as rt

DOC = "doc_88f54768a8114d5d"
RUN = ROOT / "runs/current"
DATASET = json.loads((RUN / "dataset.json").read_text())
BUNDLE = json.loads((RUN / "submitted.json").read_text())
BUNDLE["molecules"] = []
BUNDLE["reactions"] = []
BUNDLE["source_context"] = []
MOLECULES = {}
REACTIONS = []
STEREO_CHECKS = {}


def ev(page):
    return [f"{DOC}:page{page}:text", f"{DOC}:page{page}:image"] if page <= 6 else [f"{DOC}:page{page}:text"]


def f(raw=None, normalized=None, evidence=(), method="original_extraction", status=None,
      derivation=None, verification=None, alternatives=()):
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   derivation=derivation, verification=verification, alternatives=alternatives)


def missing(status="unreported"):
    return f(status=status, method="not_applicable" if status == "not_applicable" else "original_extraction")


def direct(raw, normalized=None, evidence=()):
    return f(raw, raw if normalized is None else normalized, evidence)


def inferred(normalized, evidence, derivation, raw=None):
    return f(raw, normalized, evidence, method="inference", derivation=derivation)


def derived(raw, normalized, evidence, derivation):
    return f(raw, normalized, evidence, method="structure_derivation", derivation=derivation)


def conflict(raw, evidence, alternatives):
    return f(raw, None, evidence, status="conflict", alternatives=alternatives)


def q(number, unit, raw_text=None):
    out = {"number": number, "unit": unit}
    if raw_text:
        out["raw_text"] = raw_text
    return out


def mol(mid, label, page, name=None, smiles=None, chirality="unresolved", geometry="unresolved",
        representation=None, stereo_assignment=None, identity_note=None):
    evidence = ev(page)
    if mid in MOLECULES:
        occ = MOLECULES[mid]["occurrences"][0]
        for eid in evidence:
            if eid not in occ["evidence_ids"]:
                occ["evidence_ids"].append(eid)
            if eid not in occ["compound_label"]["evidence_ids"]:
                occ["compound_label"]["evidence_ids"].append(eid)
        return mid
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = direct(label, label, evidence)
    m["names"]["original"] = direct(name or label, name or label, evidence)
    m["names"]["english"] = direct(name or label, name or label, evidence)
    m["names"]["abbreviations"] = missing()
    m["names"]["aliases"] = missing()
    m["cas_numbers"] = missing()
    m["external_properties"] = []
    if smiles:
        diag = st.stereo_inspect(smiles)
        STEREO_CHECKS[mid] = diag
        canonical = diag["canonical_isomeric_smiles"]
        m["representation_type"] = "smiles"
        m["smiles"] = derived({"source": "named compound and displayed structure", "label": label}, canonical,
                              evidence, "Candidate SMILES transcribed from the displayed/name-defined connectivity; parsing is not source verification.")
        m["stereochemistry"]["isomeric_smiles"] = derived(label, canonical, evidence,
                                                            "Canonical isomeric SMILES generated from the candidate graph.")
        m["stereochemistry"]["identity_key"] = derived(label, diag["stereo_graph_key"], evidence,
                                                         "RDKit stereo graph key from the candidate SMILES.")
        m["stereochemistry"]["chirality_status"] = inferred(chirality, evidence,
            "Status assigned from the displayed/named structure; absence of an assignment was not treated as racemic.")
        m["stereochemistry"]["geometry_status"] = inferred(geometry, evidence,
            "Status assigned from the displayed/named structure and toolkit diagnostic.")
        rep = representation or ("partial" if any(p["specified"] == "Unspecified" for p in diag["potential_stereo"]) else "encoded")
        m["stereochemistry"]["representation_status"] = inferred(rep, evidence,
            "Candidate graph was inspected with standard stereo; partial retains toolkit-unassigned potential stereo.")
        m["stereochemistry"]["atom_centers"] = derived(label, diag["atom_centers"], evidence,
                                                         "Toolkit diagnostic on candidate SMILES.")
        m["stereochemistry"]["double_bonds"] = derived(label, diag["double_bonds"], evidence,
                                                         "Toolkit diagnostic on candidate SMILES.")
        m["stereochemistry"]["axial"] = inferred("not_applicable" if "C=C=C" not in canonical else "unassigned terminal allene", evidence,
                                                   "The displayed allenes have a terminal CH2 end; no axial descriptor is reported.")
        m["stereochemistry"]["mixture"] = missing("not_applicable")
        if stereo_assignment:
            m["stereochemistry"]["assignment"] = direct(stereo_assignment, stereo_assignment, evidence)
        else:
            m["stereochemistry"]["assignment"] = missing()
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing()
        for key, val in (("chirality_status", chirality), ("geometry_status", geometry),
                         ("representation_status", representation or "unresolved")):
            m["stereochemistry"][key] = inferred(val, evidence,
                "Source label/drawing was retained, but a complete stereo-safe candidate graph was not asserted.")
        m["stereochemistry"]["assignment"] = (direct(stereo_assignment, stereo_assignment, evidence)
                                                 if stereo_assignment else missing())
        m["stereochemistry"]["atom_centers"] = missing()
        m["stereochemistry"]["double_bonds"] = missing()
        m["stereochemistry"]["axial"] = missing()
        m["stereochemistry"]["mixture"] = missing()
        m["stereochemistry"]["isomeric_smiles"] = missing()
        m["stereochemistry"]["identity_key"] = missing()
    m["identity_note"] = direct(identity_note, identity_note, evidence) if identity_note else missing()
    MOLECULES[mid] = m
    return mid


def agent(name, page, mid=None):
    mid = mid or "agent_" + "".join(c.lower() if c.isalnum() else "_" for c in name).strip("_")
    return mol(mid, name, page, name=name, chirality="unresolved", geometry="unresolved", representation="unresolved")


def participant(pid, mid, role, evidence, stage="reaction", mass=None, amount=None,
                volume=None, equivalents=None, mol_percent=None, reference=None, product_role=None):
    p = st.new_participant(pid, mid, role, evidence, stage)
    p["product_role"] = direct(product_role, product_role, evidence) if product_role else missing("not_applicable")
    vals = {"mass": mass, "amount": amount, "volume": volume, "equivalents": equivalents, "mol_percent": mol_percent}
    for key, val in vals.items():
        p["charges"][key] = direct(val.get("raw", val), val.get("value", val), evidence) if val is not None else missing()
    p["equivalence_reference"] = direct(reference, reference, evidence) if reference else missing()
    return p


def step(step_id, stage, pids, evidence, operation=None, temperature=None, time=None,
         atmosphere=None, pressure=None, inherited=False):
    s = st.new_step(step_id, evidence, stage)
    s["participant_ids"] = list(pids)
    method = "inherited_general_condition" if inherited else "original_extraction"
    def sf(val):
        if val is None:
            return missing()
        raw = val.get("raw", val) if isinstance(val, dict) else val
        norm = val.get("value", val) if isinstance(val, dict) else val
        return f(raw, norm, evidence, method=method)
    s["operation"] = sf(operation)
    s["temperature"] = sf(temperature)
    s["time"] = sf(time)
    s["atmosphere"] = sf(atmosphere)
    s["pressure"] = sf(pressure)
    return s


def outcome(oid, kind, targets, evidence, value=None, yield_type=None, basis=None, scope=None,
            product_role=None, value_status=None):
    o = st.new_outcome(oid, kind, targets, evidence)
    if value_status == "unreadable":
        o["value"] = f(value, None, evidence, status="unreadable")
    else:
        if value is None:
            o["value"] = missing()
        elif isinstance(value, dict) and ("raw" in value or "value" in value):
            o["value"] = direct(value.get("raw", value), value.get("value", value), evidence)
        else:
            o["value"] = direct(value, value, evidence)
    o["yield_type"] = direct(yield_type, yield_type, evidence) if yield_type else missing()
    o["measurement_basis"] = direct(basis, basis, evidence) if basis else missing()
    o["scope"] = direct(scope, scope, evidence) if scope else missing()
    o["product_role"] = direct(product_role, product_role, evidence) if product_role else missing("not_applicable")
    return o


def classification(changes, evidence, limitations=()):
    c = {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": missing(),
        "bond_changes": [],
        "labels": None,
        "coverage": inferred("partial" if changes else "unresolved", evidence,
                             "Candidate bond changes were inferred from the displayed net structures; no complete atom mapping was asserted."),
        "limitations": list(limitations) or ["No atom-mapped source structures; classification is a partial machine inference."],
        "model_votes": [],
    }
    labels = {}
    for i, (a, b, before, after, kind) in enumerate(changes, 1):
        ch = {"change_id": f"b{i}", "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
              "before_order": before, "after_order": after, "change_type": kind}
        c["bond_changes"].append(derived("displayed net transformation", ch, evidence,
                                         "Abstract bond change inferred from the displayed reactant/product structures."))
        if kind in ("formed", "hydrogen_gain", "order_increased"):
            lab = st.bond_label(a, b, after)
            labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(f"b{i}")
    c["labels"] = derived("displayed net transformation", list(labels.values()), evidence,
                           "Labels generated from the partial inferred bond-change list using bond-formation-20261016.")
    return c


def reaction(rid, page, record_type="specific_experiment", categories=None, location=None,
             name=None, transformation=None, participants=(), steps=(), outcomes=(), changes=(),
             general=(), notes=(), unresolved=(), extra_evidence=()):
    evidence = []
    for eid in ev(page) + list(extra_evidence):
        if eid not in evidence:
            evidence.append(eid)
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = direct(categories or ["other"], categories or ["other"], evidence)
    r["location"] = direct(location or f"PDF page {page}", {"page": page, "location": location or "source region"}, evidence)
    r["general_scheme_ids"] = list(general)
    r["reaction_name"] = direct(name, name, evidence) if name else missing()
    r["transformation_type"] = inferred(transformation, evidence,
        "Transformation label inferred from the displayed source reaction.") if transformation else missing()
    r["participants"] = list(participants)
    r["steps"] = list(steps)
    r["outcomes"] = list(outcomes)
    r["classification"] = classification(changes, evidence)
    r["reaction_scale"] = missing()
    r["solvent_conditions"] = {k: missing() for k in ("composition", "mixture_ratio", "total_volume", "concentration", "concentration_basis")}
    r["special_conditions"] = {k: missing("not_applicable") for k in ("light_source", "wavelength", "power", "current", "voltage", "electrodes", "flow_parameters", "reactor")}
    r["procedure"] = {k: missing() for k in ("original_text", "general_procedure_reference", "modifications", "addition_order", "pretreatment")}
    r["notes"] = [direct(n, n, evidence) for n in notes]
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": list(unresolved) or ["Machine extraction requires expert chemistry review."]}
    r["completeness"] = "partial"
    REACTIONS.append(r)
    return r


def set_solvent(r, raw, normalized, evidence, volume=None, ratio=None, conflict_fact=None):
    r["solvent_conditions"]["composition"] = conflict_fact or direct(raw, normalized, evidence)
    if volume is not None:
        r["solvent_conditions"]["total_volume"] = direct(volume, volume, evidence)
    if ratio is not None:
        r["solvent_conditions"]["mixture_ratio"] = direct(ratio, ratio, evidence)


def set_procedure(r, evidence, text=None, gp=None, addition=None, pretreatment=None, modifications=None):
    if text:
        r["procedure"]["original_text"] = direct(text, text, evidence)
    if gp:
        r["procedure"]["general_procedure_reference"] = direct(gp, gp, evidence)
    if addition:
        r["procedure"]["addition_order"] = direct(addition, addition, evidence)
    if pretreatment:
        r["procedure"]["pretreatment"] = direct(pretreatment, pretreatment, evidence)
    if modifications:
        r["procedure"]["modifications"] = direct(modifications, modifications, evidence)


# Concrete source compounds with conservative candidate structures for Scheme 2 only.
mol("m_1a", "1a", 2, "p-methoxybenzylamine", "COc1ccc(CN)cc1", "achiral", "not_applicable")
mol("m_1b", "1b", 2, "benzylamine", "NCc1ccccc1", "achiral", "not_applicable")
mol("m_2a", "2a", 6, "tert-Butyl 4-Methoxybenzyl(prop-2-yn-1-yl)carbamate", "CC(C)(C)OC(=O)N(CC#C)Cc1ccc(OC)cc1", "achiral", "not_applicable")
mol("m_2b", "2b", 6, "tert-Butyl Benzyl(prop-2-yn-1-yl)carbamate", "CC(C)(C)OC(=O)N(CC#C)Cc1ccccc1", "achiral", "not_applicable")
mol("m_3a", "3a", 6, "tert-Butyl Buta-2,3-dien-1-yl(4-methoxybenzyl)carbamate", "CC(C)(C)OC(=O)N(CC=C=C)Cc1ccc(OC)cc1", "achiral", "not_applicable", "partial")
mol("m_3b", "3b", 6, "tert-Butyl Benzyl(buta-2,3-dien-1-yl)carbamate", "CC(C)(C)OC(=O)N(CC=C=C)Cc1ccccc1", "achiral", "not_applicable", "partial")
mol("m_4a", "4a", 6, "N-(4-Methoxybenzyl)buta-2,3-dien-1-amine", "C=C=CCNCc1ccc(OC)cc1", "achiral", "not_applicable", "partial")
mol("m_4b", "4b", 6, "N-Benzylbuta-2,3-dien-1-amine", "C=C=CCNCc1ccccc1", "achiral", "not_applicable", "partial")
mol("m_5a", "5a", 6, "(S)-tert-Butyl 2-(Buta-2,3-dien-1-yl(4-methoxybenzyl)carbamoyl)pyrrolidine-1-carboxylate", "CC(C)(C)OC(=O)N1CCC[C@H]1C(=O)N(CC=C=C)Cc1ccc(OC)cc1", "specified", "not_applicable", "partial", "(S) stated in the compound name")
mol("m_5b", "5b", 6, "(S)-tert-Butyl 2-(Benzyl(buta-2,3-dien-1-yl)carbamoyl)pyrrolidine-1-carboxylate", "CC(C)(C)OC(=O)N1CCC[C@H]1C(=O)N(CC=C=C)Cc1ccccc1", "specified", "not_applicable", "partial", "(S) stated in the compound name")
mol("m_6a", "6a", 6, "(S)-N-(Buta-2,3-dien-1-yl)-N-(4-methoxybenzyl)pyrrolidine-2-carboxamide", "N1CCC[C@H]1C(=O)N(CC=C=C)Cc1ccc(OC)cc1", "specified", "not_applicable", "partial", "(S) stated in the compound name")
mol("m_6b", "6b", 6, "(S)-N-Benzyl-N-(buta-2,3-dien-1-yl)pyrrolidine-2-carboxamide", "N1CCC[C@H]1C(=O)N(CC=C=C)Cc1ccccc1", "specified", "not_applicable", "partial", "(S) stated in the compound heading; procedure continues on page 7")

# Products and route intermediates whose stereo-safe structures are retained as source labels/drawings.
mol("m_7a", "7a", 3, "cis bicyclic product 7a", chirality="partially_specified", geometry="not_applicable", representation="not_encoded", stereo_assignment="Relative configuration depicted; 7b, not 7a, is explicitly assigned trans in the page-2 prose.")
mol("m_7b", "7b", 3, "trans bicyclic product 7b", chirality="partially_specified", geometry="not_applicable", representation="not_encoded", stereo_assignment="Major product explicitly established as trans by NOESY (page 2).")
for mid, label, page, name in [
    ("m_7c","7c",3,"PMB-protected trans bicyclic product 7c"),("m_8","8",3,"pyrrole 8"),
    ("m_9","9",3,"trifluoroacetate 9"),("m_10","10",3,"dibromopyrrole trifluoroacetate 10"),
    ("m_11","11",3,"dibromopyrrole alcohol 11"),("m_longamide_b","Longamide B",3,"Longamide B"),
    ("m_stylisine_d","Stylisine D",3,"Stylisine D"),("m_12ba","12ba",5,"diastereomer 12ba"),
    ("m_13","13",5,"oxidized pyrrole 13"),("m_14","14",5,"hydroboration/deprotection product 14"),
    ("m_15","15",5,"dibrominated trifluoroacetate 15"),("m_16","16",5,"dibromopyrrole alcohol 16"),
    ("m_17","17",5,"longamide derivative 17")]:
    mol(mid, label, page, name, chirality="partially_specified", geometry="unresolved", representation="not_encoded",
        stereo_assignment="Stereochemical bonds are depicted in the source; absolute/relative assignment is not fully encoded in a candidate graph.")

# Frequently used agents.
for name, page in [
    ("Boc2O",2),("NaOH",2),("THF",2),("H2O",2),("propargyl bromide",2),("NaH",2),("DMF",2),
    ("paraformaldehyde",2),("diisopropylamine",2),("CuBr",2),("dioxane",2),("TFA",2),("CH2Cl2",2),
    ("Boc-L-proline",2),("EDCI",2),("HOBt",2),("DIPEA",2),("silica gel",6),("ethyl acetate",6),
    ("saturated NaHCO3",6),("brine",6),("anhydrous Na2SO4",6),("AgNO3",3),("CaCO3",3),
    ("acetone",3),("MnO2",3),("9-BBN",3),("H2O2",3),("anisole",3),("NBS",3),
    ("TEMPO",3),("PhI(OAc)2",3),("NaHCO3",3),("MeCN",3),("NH4OH",3),("Pd(OAc)2",4),
    ("PPh3",4),("XantPhos",4),("K2CO3",4),("3-chloro-4-methyliodobenzene",4),
    ("1,2-dichloro-4-iodobenzene",4),("1-iodo-3,5-dimethoxybenzene",4),("1-fluoro-4-iodobenzene",4),
    ("2-iodonaphthalene",4),("AuCl3",3),("molecular sieves 4A",3),("Ph3PAuCl",3),("Cu(OTf)2",3),
    ("TiCl4",3),("K2PtCl4",3),("Ag2O",3),("Ag2CO3",3),("silver acetate",3),("AgOTf",3)]:
    agent(name, page)


# Generic schemes and proposed strategies.
g1 = reaction("scheme1b_generic_annulation", 2, "generic_scheme", ["other"], "Scheme 1b",
              "Annulations of pyrrolidine ring", "generic pyrrolidine annulation",
              notes=["Background generic transformation; no conditions or current-paper yield are reported."],
              unresolved=["Substituents X and n are generic; this is not counted as a concrete experiment."])
g2 = reaction("scheme1c_synthetic_strategy", 2, "generic_scheme", ["other"], "Scheme 1c",
              "Strategy for pyrroloketopiperazine synthesis", "retrosynthetic strategy",
              notes=["Retrosynthetic arrows and proposed route, not an executed experiment on this page."],
              unresolved=["No experimental quantities or yield are attached to the strategy drawing."])
g_table1 = reaction("table1_generic_ag_annulation", 3, "generic_scheme", ["optimization"], "Table 1 header",
                    "Metal-catalysed cyclisation of allene-amine", "intramolecular allene amination",
                    notes=["Generic header for the 20 separately recorded optimization rows."],
                    unresolved=["Header alone does not identify a single execution."])
g4 = reaction("scheme4_generic_access", 4, "generic_scheme", ["other"], "Scheme 4",
              "Access to alpha-substituted longamide derivatives", "retrosynthetic strategy",
              notes=["Conceptual access route; no execution conditions or yield are shown."],
              unresolved=["Ar is generic and no concrete experiment is attached to Scheme 4."])
g_t2 = reaction("table2_generic_pd_cyclisation", 4, "generic_scheme", ["optimization"], "Table 2 header",
                "Pd-catalysed cyclisation of proline derivative 6a", "Pd-catalysed allene cyclisation",
                notes=["Generic header; entries a-c are separate executions."], unresolved=["Not a single execution."])
g_t3 = reaction("table3_generic_pd_scope", 4, "generic_scheme", ["substrate_scope"], "Table 3 header",
                "Variation of the iodide in the Pd-catalysed cyclisation", "Pd-catalysed allene arylative cyclisation",
                notes=["Generic header; entries a-e are separate executions."], unresolved=["PG and Ar are generic in the header."])


# Scheme 2 / page-6 concrete preparations.  Each product is a separate execution.
def route2(rid, start, product, product_label, yield_pct, page6_data, steps_data, reagents, changes, notes=(), unresolved=()):
    evidence = ev(2) + ev(6)
    ps = [participant("substrate", start, "reactant", evidence,
                      mass=page6_data.get("start_mass"), amount=page6_data.get("start_amount"))]
    for i, spec in enumerate(reagents):
        mid = agent(spec[0], spec[1] if len(spec) > 1 else 2)
        kwargs = spec[2] if len(spec) > 2 else {}
        ps.append(participant(f"reagent_{i+1}", mid, kwargs.pop("role", "reagent"), evidence, **kwargs))
    ps.append(participant("product", product, "product", evidence, mass=page6_data.get("product_mass"),
                          amount=page6_data.get("product_amount"), product_role="main product"))
    ss = [step(f"s{i+1}", d.get("stage", "reaction"), d.get("pids", [p["participant_id"] for p in ps]), evidence,
               d.get("operation"), d.get("temperature"), d.get("time"), d.get("atmosphere"), d.get("pressure"), d.get("inherited", False))
          for i, d in enumerate(steps_data)]
    outs = []
    if yield_pct is not None:
        outs.append(outcome("yield", "yield", [product], evidence, q(yield_pct, "%"),
                            "reported product yield", None, {"start": start, "end": product}, "main product"))
    r = reaction(rid, 6, "specific_experiment", ["other"], f"Scheme 2; Experimental Section product {product_label}",
                 None, None, ps, ss, outs, changes, notes=notes,
                 unresolved=list(unresolved) + ["Yield measurement basis is not reported; literature-referenced operations remain incomplete."],
                 extra_evidence=ev(2))
    r["reaction_scale"] = direct(page6_data.get("start_amount"),
        {"reference_molecule_id": start, "amount": page6_data.get("start_amount")}, evidence) if page6_data.get("start_amount") else missing()
    return r


r_2a = route2("scheme2_1a_to_2a", "m_1a", "m_2a", "2a", 79,
    {"start_mass":q(4.20,"g"),"start_amount":q(30.66,"mmol"),"product_mass":q(6.67,"g"),"product_amount":q(24.25,"mmol")},
    [{"operation":"Boc protection (literature procedure)","temperature":{"raw":"0 °C to r.t.","value":{"minimum":0,"maximum":25,"unit":"°C","raw_text":"0 °C to r.t."}}},
     {"operation":"Propargylation (literature procedure)","temperature":{"raw":"0 °C to r.t.","value":{"minimum":0,"maximum":25,"unit":"°C","raw_text":"0 °C to r.t."}}}],
    [("Boc2O",2,{}),("NaOH",2,{}),("THF",2,{"role":"solvent"}),("H2O",2,{"role":"solvent"}),("propargyl bromide",2,{}),("NaH",2,{}),("DMF",2,{"role":"solvent"})],
    [("C","N",0,1,"formed"),("C","N",0,1,"formed")], notes=["79% is the reported two-operation yield from p-methoxybenzylamine to 2a."])
r_2b = route2("scheme2_1b_to_2b", "m_1b", "m_2b", "2b", 72,
    {"start_mass":q(3.00,"g"),"start_amount":q(28.04,"mmol"),"product_mass":q(4.95,"g"),"product_amount":q(20.19,"mmol")},
    [{"operation":"Prepared following the procedure for 2a","temperature":{"raw":"0 °C to r.t.","value":{"minimum":0,"maximum":25,"unit":"°C","raw_text":"0 °C to r.t."}},"inherited":True}],
    [("Boc2O",2,{}),("NaOH",2,{}),("THF",2,{"role":"solvent"}),("H2O",2,{"role":"solvent"}),("propargyl bromide",2,{}),("NaH",2,{}),("DMF",2,{"role":"solvent"})],
    [("C","N",0,1,"formed"),("C","N",0,1,"formed")])
r_3a = route2("scheme2_2a_to_3a", "m_2a", "m_3a", "3a", 72,
    {"start_mass":q(5.00,"g"),"start_amount":q(18.18,"mmol"),"product_mass":q(3.79,"g"),"product_amount":q(13.11,"mmol")},
    [{"operation":"Crabbé allene homologation (literature procedure)","temperature":q(110,"°C")}],
    [("paraformaldehyde",2,{}),("diisopropylamine",2,{}),("CuBr",2,{}),("dioxane",2,{"role":"solvent"})],
    [("C","C",0,1,"formed"),("C","C",2,1,"order_decreased"),("C","C",1,2,"order_increased")])
r_3b = route2("scheme2_2b_to_3b", "m_2b", "m_3b", "3b", 69,
    {"start_mass":q(4.95,"g"),"start_amount":q(20.19,"mmol"),"product_mass":q(3.61,"g"),"product_amount":q(13.93,"mmol")},
    [{"operation":"Prepared following the procedure for 3a","temperature":q(110,"°C"),"inherited":True}],
    [("paraformaldehyde",2,{}),("diisopropylamine",2,{}),("CuBr",2,{}),("dioxane",2,{"role":"solvent"})],
    [("C","C",0,1,"formed"),("C","C",2,1,"order_decreased"),("C","C",1,2,"order_increased")])

# Detailed deprotections/couplings on page 6.
e26 = ev(2)+ev(6)
r_4a = route2("scheme2_3a_to_4a", "m_3a", "m_4a", "4a", 96,
    {"start_mass":q(3.41,"g"),"start_amount":q(11.80,"mmol"),"product_mass":q(2.14,"g"),"product_amount":q(11.32,"mmol")},
    [{"operation":"TFA added to stirred solution; stir","temperature":"room temperature","time":q(12,"h")},
     {"stage":"workup","operation":"evaporate; redissolve in ethyl acetate; wash 3 × saturated NaHCO3 then brine; dry Na2SO4; evaporate"}],
    [("CH2Cl2",6,{"role":"solvent","volume":q(25,"mL")}),("TFA",6,{"volume":q(9.0,"mL"),"amount":q(118,"mmol")}),
     ("ethyl acetate",6,{"role":"workup_agent","stage":"workup","volume":q(30,"mL")}),("saturated NaHCO3",6,{"role":"workup_agent","stage":"workup","volume":{"raw":"3 × 30 mL","value":{"number":90,"unit":"mL","raw_text":"3 × 30 mL"}}}),("brine",6,{"role":"workup_agent","stage":"workup","volume":q(30,"mL")}),("anhydrous Na2SO4",6,{"role":"workup_agent","stage":"workup"})],
    [("N","H",0,1,"hydrogen_gain")])
set_solvent(r_4a,"dichloromethane (25 mL)",["CH2Cl2"],e26,q(25,"mL"))
set_procedure(r_4a,e26,"TFA was added to a stirred solution of 3a in dichloromethane; after 12 h at room temperature the mixture was evaporated and worked up as reported.", addition="TFA added to stirred substrate solution")

r_4b = route2("scheme2_3b_to_4b", "m_3b", "m_4b", "4b", 99,
    {"start_mass":q(923,"mg"),"start_amount":q(3.56,"mmol"),"product_mass":q(561,"mg"),"product_amount":q(3.53,"mmol")},
    [{"operation":"Prepared following procedure for 4a","temperature":"room temperature","time":q(12,"h"),"inherited":True},
     {"stage":"workup","operation":"Work-up inherited from 4a; product-specific solvent volumes not reported","inherited":True}],
    [("CH2Cl2",6,{"role":"solvent"}),("TFA",6,{"volume":q(2.73,"mL"),"amount":q(35.6,"mmol")}),
     ("ethyl acetate",6,{"role":"workup_agent","stage":"workup"}),("saturated NaHCO3",6,{"role":"workup_agent","stage":"workup"}),("brine",6,{"role":"workup_agent","stage":"workup"}),("anhydrous Na2SO4",6,{"role":"workup_agent","stage":"workup"})],
    [("N","H",0,1,"hydrogen_gain")], notes=["The procedure is inherited from 4a; exact product-specific work-up volumes are not supplied."])
set_solvent(r_4b,"dichloromethane (procedure for 4a)",["CH2Cl2"],e26)
set_procedure(r_4b,e26,gp="procedure for 4a")

def coupling_record(rid,start,product,label,yld,data,inherited=False):
    evidence=e26
    specs=[("Boc-L-proline",{"mass":data["pro_mass"],"amount":data["pro_amount"]}),
           ("EDCI",{"mass":data["edci_mass"],"amount":data["edci_amount"]}),
           ("HOBt",{"mass":data["hobt_mass"],"amount":data["hobt_amount"]}),
           ("DIPEA",{"volume":data["dipea_volume"],"amount":data["dipea_amount"]}),
           ("CH2Cl2",{"role":"solvent","volume":data.get("dcm_volume")}),
           ("silica gel",{"role":"purification_agent","stage":"purification"})]
    ps=[participant("substrate",start,"reactant",evidence,mass=data["start_mass"],amount=data["start_amount"])]
    for i,(name,kw) in enumerate(specs):
        kw={k:v for k,v in kw.items() if v is not None}; role=kw.pop("role","reagent"); stage=kw.pop("stage","reaction")
        ps.append(participant(f"agent_{i+1}",agent(name,6),role,evidence,stage=stage,**kw))
    ps.append(participant("product",product,"product",evidence,mass=data["product_mass"],amount=data["product_amount"],product_role="main product"))
    ss=[step("s1","reaction",[p["participant_id"] for p in ps if p["stage"]=="reaction"],evidence,
             "cool coupling mixture then add amine substrate",q(0,"°C"),None,"nitrogen" if not inherited else None,inherited=inherited),
        step("s2","reaction",[p["participant_id"] for p in ps if p["stage"]=="reaction"],evidence,
             "stir", "room temperature", q(12,"h"), inherited=inherited),
        step("s3","purification",[p["participant_id"] for p in ps if p["stage"]=="purification"],evidence,
             "concentrate under reduced pressure and purify by silica-gel chromatography", inherited=inherited)]
    r=reaction(rid,6,"specific_experiment",["other"],f"Scheme 2; Experimental Section product {label}",None,"amide coupling",ps,ss,
               [outcome("yield","yield",[product],evidence,q(yld,"%"),"reported product yield",None,{"start":start,"end":product},"main product")],
               [("C","N",0,1,"formed")],extra_evidence=ev(2),
               notes=["Reported reagent masses and mmol values are retained verbatim even where their arithmetic appears inconsistent."],
               unresolved=["Yield measurement basis is not reported.","EDCI/HOBt mass-to-amount arithmetic in the 5a paragraph requires expert/source review; values were not repaired."] if not inherited else ["Yield measurement basis is not reported.","Exact inherited CH2Cl2 volume and atmosphere for 5b are not product-specifically reported."])
    set_solvent(r,"dichloromethane",["CH2Cl2"],evidence,data.get("dcm_volume"))
    set_procedure(r,evidence,gp="procedure for 5a" if inherited else None,addition="coupling reagents to Boc-L-proline, cool, then add amine substrate")
    return r

r_5a=coupling_record("scheme2_4a_to_5a","m_4a","m_5a","5a",88,
    {"start_mass":q(1.70,"g"),"start_amount":q(9.00,"mmol"),"pro_mass":q(2.13,"g"),"pro_amount":q(9.89,"mmol"),
     "edci_mass":q(2.16,"g"),"edci_amount":q(1.35,"mmol"),"hobt_mass":q(1.80,"g"),"hobt_amount":q(1.35,"mmol"),
     "dipea_volume":q(3.12,"mL"),"dipea_amount":q(19.78,"mmol"),"dcm_volume":q(30,"mL"),
     "product_mass":q(3.36,"g"),"product_amount":q(8.70,"mmol")})
r_5b=coupling_record("scheme2_4b_to_5b","m_4b","m_5b","5b",76,
    {"start_mass":q(470,"mg"),"start_amount":q(2.96,"mmol"),"pro_mass":q(700,"mg"),"pro_amount":q(3.25,"mmol"),
     "edci_mass":q(682,"mg"),"edci_amount":q(3.55,"mmol"),"hobt_mass":q(480,"mg"),"hobt_amount":q(3.55,"mmol"),
     "dipea_volume":q(1.03,"mL"),"dipea_amount":q(5.92,"mmol"),"dcm_volume":None,
     "product_mass":q(800,"mg"),"product_amount":q(2.25,"mmol")},True)

r_6a=route2("scheme2_5a_to_6a","m_5a","m_6a","6a",99,
    {"start_mass":q(3.00,"g"),"start_amount":q(7.77,"mmol"),"product_mass":q(2.20,"g"),"product_amount":q(7.69,"mmol")},
    [{"operation":"TFA added to stirred solution; stir","temperature":"room temperature","time":q(12,"h")},
     {"stage":"workup","operation":"evaporate; redissolve in ethyl acetate; wash 3 × saturated NaHCO3 and once with brine; dry Na2SO4; evaporate"}],
    [("CH2Cl2",6,{"role":"solvent","volume":q(20,"mL")}),("TFA",6,{"volume":q(5.90,"mL"),"amount":q(77.7,"mmol")}),
     ("ethyl acetate",6,{"role":"workup_agent","stage":"workup","volume":q(30,"mL")}),("saturated NaHCO3",6,{"role":"workup_agent","stage":"workup","volume":{"raw":"3 × 30 mL","value":{"number":90,"unit":"mL","raw_text":"3 × 30 mL"}}}),("brine",6,{"role":"workup_agent","stage":"workup","volume":q(30,"mL")}),("anhydrous Na2SO4",6,{"role":"workup_agent","stage":"workup"})],
    [("N","H",0,1,"hydrogen_gain")])
set_solvent(r_6a,"dichloromethane (20 mL)",["CH2Cl2"],e26,q(20,"mL"))
set_procedure(r_6a,e26,"TFA was added to 5a in dichloromethane; after 12 h at room temperature the mixture was evaporated and worked up as reported.",addition="TFA added to stirred substrate solution")

# 6b heading is visible on page 6, but its procedure continues outside the bounded range.
r_6b=route2("scheme2_5b_to_6b","m_5b","m_6b","6b",None,{},
    [{"operation":"Prepared following the procedure for 6a; continuation is on page 7","temperature":"room temperature","inherited":True}],
    [("CH2Cl2",6,{"role":"solvent"}),("TFA",6,{})],[("N","H",0,1,"hydrogen_gain")],
    notes=["The page-6 heading and first line establish the experiment, but the product-specific charges and yield continue beyond page 6."],
    unresolved=["Product-specific quantities and yield are outside the bounded page-1--6 extraction and were not imported."])
set_solvent(r_6b,"dichloromethane (inherited procedure)",["CH2Cl2"],e26)
set_procedure(r_6b,e26,gp="procedure for 6a; text continues on page 7")


# Table 1: one record per optimization row.
table1 = [
 ("a",["AuCl3"],5,"DCM","5 d",48,"39:61"),("b",["AuCl3","CaCO3","molecular sieves 4A"],5,"DCM + 5 % ACN","5 d",55,"21:79"),
 ("c",["Ph3PAuCl"],5,"DCM","5 d","/","/"),("d",["Cu(OTf)2"],5,"DCM","5 d","Trace","/"),
 ("e",["Cu(OTf)2"],5,"Acetone","5 d",57,"21:79"),("f",["Cu(OTf)2","CaCO3","molecular sieves 4A"],5,"Acetone","5 d","Trace","/"),
 ("g",["TiCl4"],5,"DCM","5 d","/","/"),("h",["K2PtCl4"],20,"Acetone","5 d","/","/"),
 ("i",["Ag2O"],20,"Acetone","18 h",28,"29:71"),("j",["Ag2CO3"],20,"Acetone","18 h","Trace","/"),
 ("k",["AgNO3"],20,"Acetone","18 h",19,"31:69"),("l",["silver acetate"],20,"Acetone","18 h",21,"30:70"),
 ("m",["AgOTf"],20,"Acetone","3 d","/","/"),("n",["AgOTf"],20,"DCM","18 h",38,"26:74"),
 ("o",["AgNO3","K2CO3"],20,"Acetone + H2O","18 h",25,"33:67"),("p",["AgNO3","CaCO3"],20,"Acetone + H2O","18 h",21,"31:69"),
 ("q",["AgNO3","CaCO3"],40,"Acetone + H2O","18 h",55,"14:86"),("r",["AgNO3","CaCO3"],80,"Acetone + H2O","18 h",83,"17:83"),
 ("s",["AgOTf","CaCO3"],20,"Acetone + H2O","18 h",50,"19:81"),("t",["AgOTf","CaCO3"],40,"Acetone + H2O","18 h",74,"23:77")]
for row,names,loading,solv,t,yld,ratio in table1:
    evidence=ev(3)
    ps=[participant("substrate","m_6b","reactant",evidence)]
    primary=agent(names[0],3)
    ps.append(participant("metal",primary,"catalyst",evidence,mol_percent=q(loading,"mol%"),reference="m_6b"))
    for i,n in enumerate(names[1:],1):ps.append(participant(f"additive_{i}",agent(n,3),"additive",evidence))
    solvent_names=[x.strip() for x in solv.replace(" + 5 % ACN"," + ACN").split("+")]
    for i,n in enumerate(solvent_names):ps.append(participant(f"solvent_{i+1}",agent(n.strip(),3),"solvent",evidence))
    has_products=yld != "/"
    if has_products:
        ps += [participant("product_7a","m_7a","product",evidence,product_role="diastereomer 7a"),participant("product_7b","m_7b","product",evidence,product_role="diastereomer 7b")]
    outs=[]
    if isinstance(yld,(int,float)):
        basis="isolated" if row=="r" else "NMR"
        outs.append(outcome("yield","yield",["m_7a","m_7b"],evidence,q(yld,"%"),"combined yield",basis,{"table_entry":row},"combined diastereomers"))
    elif yld=="Trace":
        outs.append(outcome("trace","trace",["m_7a","m_7b"],evidence,"Trace",None,"NMR",{"table_entry":row},"combined diastereomers"))
    if ratio!="/":
        a,b=map(int,ratio.split(":"))
        outs.append(outcome("ratio","selectivity",["m_7a","m_7b"],evidence,{"raw":ratio,"value":{"metric":"diastereomer_ratio","order":["7a","7b"],"ratio":[a,b]}},None,None,{"order":["m_7a","m_7b"]},"diastereomer ratio"))
    r=reaction(f"table1_{row}",3,"specific_experiment",["optimization"],f"Table 1 entry {row}","Optimization of metal-catalysed cyclisation of allene-amine","intramolecular allene amination",ps,
               [step("reaction","reaction",[p["participant_id"] for p in ps if p["stage"]=="reaction"],evidence,"stir",None,{"raw":t,"value":q(float(t.split()[0]),"d" if t.endswith("d") else "h")})],outs,
               [("C","N",0,1,"formed"),("C","H",0,1,"hydrogen_gain"),("C","C",2,1,"order_decreased")] if has_products else [],
               general=[g_table1["reaction_id"]],notes=[f"Table yield marker: {yld}; ratio marker: {ratio}."],
               unresolved=["Absolute substrate scale, atmosphere, pressure, work-up and purification are not reported in Table 1."] + ([] if row=="r" else ["Ratio measurement method is not reported; the NMR column heading applies to yield only."]))
    set_solvent(r,solv,solvent_names,evidence,ratio="5% ACN" if "5 % ACN" in solv else None)


# Scheme 3 concrete route (source-level conditions; product-specific quantities are outside pages 1-6).
scheme3_specs=[
 ("scheme3_6a_to_7c","m_6a","m_7c",69,["AgNO3","CaCO3","acetone","H2O"],"room temperature",None,[("C","N",0,1,"formed"),("C","H",0,1,"hydrogen_gain"),("C","C",2,1,"order_decreased")],"Ag-catalysed cyclisation"),
 ("scheme3_7c_to_8","m_7c","m_8",57,["MnO2","THF"],q(85,"°C"),None,[],"oxidation"),
 ("scheme3_8_to_9","m_8","m_9",55,["9-BBN","THF","NaOH","H2O2","anisole","TFA"],q(60,"°C"),None,[("C","O",0,1,"formed"),("C","H",0,1,"hydrogen_gain"),("N","H",0,1,"hydrogen_gain")],"two-stage hydroboration/oxidation and PMB deprotection"),
 ("scheme3_9_to_10","m_9","m_10",85,["NBS","THF"],"room temperature",None,[("C","Br",0,1,"formed"),("C","Br",0,1,"formed")],"bromination"),
 ("scheme3_10_to_11","m_10","m_11",99,["NaOH","THF","H2O"],"room temperature",None,[("O","H",0,1,"hydrogen_gain")],"trifluoroacetate hydrolysis"),
 ("scheme3_11_to_longamide_b","m_11","m_longamide_b",83,["TEMPO","PhI(OAc)2","NaHCO3","MeCN","H2O"],"room temperature",None,[("C","O",1,2,"order_increased"),("C","O",0,1,"formed")],"alcohol oxidation"),
 ("scheme3_longamide_b_to_stylisine_d","m_longamide_b","m_stylisine_d",81,["NH4OH","EDCI","HOBt","DIPEA","THF"],"room temperature",None,[("C","N",0,1,"formed")],"amidation")]
for rid,start,prod,yld,names,temp,tim,changes,name in scheme3_specs:
    evidence=ev(3);ps=[participant("substrate",start,"reactant",evidence)]
    for i,n in enumerate(names):ps.append(participant(f"agent_{i+1}",agent(n,3),"solvent" if n in ("THF","acetone","H2O","MeCN") else "reagent",evidence))
    ps.append(participant("product",prod,"product",evidence,product_role="main product"))
    ss=[step("reaction","reaction",[p["participant_id"] for p in ps],evidence,"conditions shown in Scheme 3",temp,tim)]
    if rid=="scheme3_8_to_9":
        ss=[step("hydroboration_oxidation","reaction",[p["participant_id"] for p in ps],evidence,"(i) 9-BBN in THF, then NaOH/H2O2"),
            step("deprotection","reaction",[p["participant_id"] for p in ps],evidence,"(ii) anisole/TFA",q(60,"°C"))]
    r=reaction(rid,3,"specific_experiment",["derivatization"],"Scheme 3",name,name,ps,ss,
               [outcome("yield","yield",[prod],evidence,q(yld,"%"),"reported scheme yield",None,{"start":start,"end":prod,"overall_two_steps":rid=="scheme3_8_to_9"},"main product")],changes,
               unresolved=["Absolute scale, reaction time, atmosphere, work-up, purification and yield basis are not reported in the bounded Scheme 3 evidence.","Complex structures/mapping require expert review."])
    solvents=[n for n in names if n in ("THF","acetone","H2O","MeCN")]
    if solvents:set_solvent(r,"/".join(solvents),solvents,evidence,ratio="1:1" if set(solvents)=={"acetone","H2O"} or set(solvents)=={"MeCN","H2O"} else None)


# Table 2 optimization with explicit source conflicts retained.
for row,ligand,solvent,yld in [("a","PPh3","DMF",20),("b","XantPhos","DMF",38),("c","XantPhos","MeCN",58)]:
    evidence=ev(4)
    ar=agent("3-chloro-4-methyliodobenzene",4)
    ps=[participant("allene","m_6a","reactant",evidence,amount=q(0.1,"mmol")),
        participant("aryl_iodide",ar,"reactant",evidence,amount=q(0.1,"mmol")),
        participant("pd",agent("Pd(OAc)2",4),"catalyst",evidence,amount=q(0.01,"mmol")),
        participant("ligand",agent(ligand,4),"ligand",evidence,amount=q(0.02 if ligand=="PPh3" else 0.01,"mmol")),
        participant("base",agent("K2CO3",4),"reagent",evidence,amount=q(0.23,"mmol")),
        participant("solvent",agent(solvent,4),"solvent",evidence,volume=q(3,"mL")),
        participant("product_a",agent("12aa",4,"m_12aa"),"product",evidence,product_role="diastereomer 12aa"),
        participant("product_b","m_12ba","product",evidence,product_role="diastereomer 12ba")]
    eq_conf=conflict({"prose":"ArI 1.2 equivalents","footnote":"allene 0.1 mmol; ArI 0.1 mmol"},evidence,[
        direct("1.2 equivalents",q(1.2,"equiv"),evidence),
        inferred(q(1.0,"equiv"),evidence,"0.1 mmol ArI divided by 0.1 mmol allene in the table footnote.","0.1 mmol / 0.1 mmol")])
    ps[1]["charges"]["equivalents"]=eq_conf;ps[1]["equivalence_reference"]=direct("allene 6a","m_6a",evidence)
    r=reaction(f"table2_{row}",4,"specific_experiment",["optimization"],f"Table 2 entry {row}","Pd-catalysed cyclisation of proline derivative 6a","Pd-catalysed allene arylative cyclisation",ps,
               [step("reaction","reaction",[p["participant_id"] for p in ps],evidence,"reflux", "reflux", q(18,"h"))],
               [outcome("yield","yield",["m_12aa","m_12ba"],evidence,q(yld,"%"),"combined yield","isolated",{"table_entry":row},"combined diastereomers")],
               [("C","C",0,1,"formed"),("C","N",0,1,"formed")],general=[g_t2["reaction_id"]],
               notes=["Source conflict retained: prose says 1.2 equivalents ArI, while the footnote lists 0.1 mmol each for allene and ArI."] + (["Source conflict retained: row solvent is DMF, while the shared footnote prints MeCN (3 mL)."] if solvent=="DMF" else []),
               unresolved=["ArI loading conflict requires expert interpretation."] + (["DMF row versus MeCN shared-footnote conflict requires expert interpretation."] if solvent=="DMF" else []))
    if solvent=="DMF":
        cf=conflict({"row":solvent,"footnote":"MeCN (3 mL)"},evidence,[direct(solvent,solvent,evidence),direct("MeCN (3 mL)","MeCN",evidence)])
        set_solvent(r,solvent,solvent,evidence,q(3,"mL"),conflict_fact=cf)
    else:set_solvent(r,"MeCN (3 mL)",["MeCN"],evidence,q(3,"mL"))


# Table 3: five concrete aryl iodides; each row has two distinct diastereomer products.
scope_rows=[("a","m_6a","3-chloro-4-methyliodobenzene",58,(55,45),"12aa","12ba"),
            ("b","m_6a","1,2-dichloro-4-iodobenzene",59,(56,44),"12ab","12bb"),
            ("c","m_6a","1-iodo-3,5-dimethoxybenzene",53,(56,44),"12ac","12bc"),
            ("d","m_6a","1-fluoro-4-iodobenzene",68,(57,43),"12ad","12bd"),
            ("e","m_6b","2-iodonaphthalene",75,(57,43),"12ae","12be")]
for row,allene,aryl_name,yld,ratio,la,lb in scope_rows:
    evidence=ev(4);ma=agent(la,4,f"m_{la}");mb="m_12ba" if lb=="12ba" else agent(lb,4,f"m_{lb}")
    ar=agent(aryl_name,4)
    ps=[participant("allene",allene,"reactant",evidence,amount=q(0.1,"mmol")),participant("aryl_iodide",ar,"reactant",evidence,amount=q(0.1,"mmol")),
        participant("pd",agent("Pd(OAc)2",4),"catalyst",evidence,amount=q(0.01,"mmol")),participant("ligand",agent("XantPhos",4),"ligand",evidence,amount=q(0.01,"mmol")),
        participant("base",agent("K2CO3",4),"reagent",evidence,amount=q(0.23,"mmol")),participant("solvent",agent("MeCN",4),"solvent",evidence,volume=q(3,"mL")),
        participant("product_a",ma,"product",evidence,product_role=f"diastereomer {la}"),participant("product_b",mb,"product",evidence,product_role=f"diastereomer {lb}")]
    outs=[outcome("yield","yield",[ma,mb],evidence,q(yld,"%"),"combined yield","isolated",{"table_entry":row},"combined diastereomers"),
          outcome("ratio","selectivity",[ma,mb],evidence,{"raw":f"{ratio[0]}:{ratio[1]}","value":{"metric":"diastereomer_ratio","order":[la,lb],"ratio":list(ratio)}},None,"1H NMR",{"order":[ma,mb]},"diastereomer ratio")]
    if row=="a":
        outs.append(outcome("yield_12ba","yield",["m_12ba"],ev(5),q(26,"%"),"product-specific yield",None,{"scheme":"Scheme 5 step a"},"diastereomer 12ba"))
    r=reaction(f"table3_{row}",4,"specific_experiment",["substrate_scope"],f"Table 3 entry {row}","Variation of the iodide in the Pd-catalysed cyclisation","Pd-catalysed allene arylative cyclisation",ps,
               [step("reaction","reaction",[p["participant_id"] for p in ps],evidence,"reflux","reflux",q(18,"h"))],outs,
               [("C","C",0,1,"formed"),("C","N",0,1,"formed")],general=[g_t3["reaction_id"]],
               notes=["Table yield is isolated combined yield; diastereomer ratio was established by 1H NMR."] + (["Scheme 5 separately reports 26% for the 12ba branch; it is retained as a product-specific outcome, not a duplicate experiment."] if row=="a" else []),
               unresolved=["Absolute configuration of the two table-header diastereomer labels requires expert review; no stereo-safe product SMILES asserted."] + (["Measurement basis for the product-specific 26% 12ba yield is not stated."] if row=="a" else []),
               extra_evidence=ev(5) if row=="a" else ())
    r["reaction_scale"]=direct("allene 0.1 mmol",{"reference_molecule_id":allene,"amount":q(0.1,"mmol")},evidence)
    set_solvent(r,"MeCN (3 mL)",["MeCN"],evidence,q(3,"mL"))


# Scheme 5 downstream derivatizations. Step a is already merged into Table 3a.
scheme5_specs=[
 ("scheme5_12ba_to_13","m_12ba","m_13",56,["MnO2","THF"],q(85,"°C"),[("C","C",1,2,"order_increased")],"oxidation"),
 ("scheme5_13_to_14","m_13","m_14",32,["9-BBN","THF","NaOH","H2O2","anisole","TFA"],q(60,"°C"),[("C","O",0,1,"formed"),("C","H",0,1,"hydrogen_gain"),("N","H",0,1,"hydrogen_gain")],"two-stage hydroboration/oxidation and PMB deprotection"),
 ("scheme5_14_to_15","m_14","m_15",75,["NBS","THF"],"room temperature",[("C","Br",0,1,"formed"),("C","Br",0,1,"formed")],"bromination"),
 ("scheme5_15_to_16","m_15","m_16",99,["NaOH","THF","H2O"],"room temperature",[("O","H",0,1,"hydrogen_gain")],"trifluoroacetate hydrolysis"),
 ("scheme5_16_to_17","m_16","m_17",83,["TEMPO","PhI(OAc)2","NaHCO3","MeCN","H2O"],"room temperature",[("C","O",1,2,"order_increased"),("C","O",0,1,"formed")],"alcohol oxidation")]
for rid,start,prod,yld,names,temp,changes,name in scheme5_specs:
    evidence=ev(5);ps=[participant("substrate",start,"reactant",evidence)]
    for i,n in enumerate(names):ps.append(participant(f"agent_{i+1}",agent(n,5),"solvent" if n in ("THF","H2O","MeCN") else "reagent",evidence))
    ps.append(participant("product",prod,"product",evidence,product_role="main product"))
    ss=[step("reaction","reaction",[p["participant_id"] for p in ps],evidence,"conditions shown in Scheme 5",temp)]
    if rid=="scheme5_13_to_14":
        ss=[step("hydroboration_oxidation","reaction",[p["participant_id"] for p in ps],evidence,"(i) 9-BBN in THF, then NaOH/H2O2"),step("deprotection","reaction",[p["participant_id"] for p in ps],evidence,"(ii) anisole/TFA",q(60,"°C"))]
    r=reaction(rid,5,"specific_experiment",["derivatization"],"Scheme 5",name,name,ps,ss,
               [outcome("yield","yield",[prod],evidence,q(yld,"%"),"reported scheme yield",None,{"start":start,"end":prod,"overall_two_steps":rid=="scheme5_13_to_14"},"main product")],changes,
               unresolved=["Absolute scale, time, atmosphere, work-up, purification and yield basis are not reported in the bounded Scheme 5 evidence.","Complex stereochemistry and atom mapping require expert review."])
    solvents=[n for n in names if n in ("THF","H2O","MeCN")]
    if solvents:set_solvent(r,"/".join(solvents),solvents,evidence,ratio="1:1" if set(solvents)=={"MeCN","H2O"} else None)


# Figure 1 and Figure 2 are proposed transition-state models, not mechanistic experiments.
reaction("figure1_pd_transition_state_model",5,"proposed_mechanism",["mechanistic_study"],"Figure 1","Transition state for the Pd-catalysed cyclisation reaction","proposed transition-state model",
         notes=["Two proposed Re/Si-face transition-state modes lead to 12a and 12b; no independent control, isotope or kinetic experiment is shown."],
         unresolved=["Author model only; not experimentally demonstrated on page 5."])
reaction("figure2_ag_transition_state_model",5,"proposed_mechanism",["mechanistic_study"],"Figure 2","Transition state for the Ag-catalysed cyclization reaction","proposed transition-state model",
         notes=["Two proposed allene orientations lead to 7a/7b; the model rationalizes trans-7b as major."],
         unresolved=["Author model only; not experimentally demonstrated on page 5."])


# Source contexts: discussion/bioactivity/characterization do not become syntheses.
BUNDLE["source_context"] = [
 {"context_id":f"{DOC}:page1_article_context","document_id":DOC,"evidence_ids":ev(1),"kind":"article_metadata_and_background",
  "summary":"Article title/authors/DOI, abstract, background bioactivity and conceptual discussion; no concrete reaction execution on page 1."},
 {"context_id":f"{DOC}:page2_scheme1_natural_products","document_id":DOC,"evidence_ids":ev(2),"kind":"structure_context",
  "summary":"Scheme 1a displays selected natural products; these structures are context, not current-paper syntheses on this page."},
 {"context_id":f"{DOC}:page5_stereochemical_discussion","document_id":DOC,"evidence_ids":ev(5),"kind":"mechanistic_discussion",
  "summary":"Author discussion attributes Ag/Pd diastereoselectivity differences to proposed steric interactions; Figures 1-2 are models."},
 {"context_id":f"{DOC}:page6_characterization","document_id":DOC,"evidence_ids":ev(6),"kind":"characterization",
  "summary":"NMR/IR/HRMS/optical-rotation/TLC data for 2a-6a accompany real preparations but do not create additional synthesis records; 5b is described as a rotamer mixture."},
]


# Fill procedure/solvent defaults for records not already specialized.
for r in REACTIONS:
    if r["record_type"] == "specific_experiment" and r["procedure"]["original_text"]["normalized_value"] is None:
        r["procedure"]["original_text"] = direct(r["location"]["raw_value"], r["location"]["raw_value"], r["evidence_ids"])

BUNDLE["molecules"] = list(MOLECULES.values())
BUNDLE["reactions"] = REACTIONS
(RUN / "submitted.json").write_text(json.dumps(BUNDLE, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/stereo_checks.json").write_text(json.dumps(STEREO_CHECKS, ensure_ascii=False, indent=2) + "\n")

inventory = {
 "document_id": DOC,
 "bounded_pages": [1,2,3,4,5,6],
 "records": [{"reaction_id":r["reaction_id"],"page":r["location"]["normalized_value"].get("page"),"record_type":r["record_type"],"location":r["location"]["raw_value"]} for r in REACTIONS],
 "page_findings": {
  "1":["article metadata, abstract, introduction; no concrete reaction execution"],
  "2":["Scheme 1 background/generic strategy","Scheme 2 five arrows expanded into compound-specific preparations","optimization discussion"],
  "3":["Table 1 rows a-t","Scheme 3 arrows a-g","stereochemical/product discussion"],
  "4":["Scheme 4 generic strategy","Table 2 rows a-c","Table 3 rows a-e"],
  "5":["Scheme 5; step a merged with Table 3a, steps b-f separate","Figures 1-2 proposed transition-state models","conclusion"],
  "6":["Experimental procedures/characterization for 2a-6a and start of 6b","spectral data retained as context, not extra reactions"]}}
(RUN / "raw_model/page_visual_inventory.json").write_text(json.dumps(inventory, ensure_ascii=False, indent=2) + "\n")

mentions = st.source_mentions(DATASET, ["6a","6b","7a","7b","12ba","Longamide B","17"], 180)
(RUN / "raw_model/source_mentions.json").write_text(json.dumps(mentions, ensure_ascii=False, indent=2) + "\n")

by_page={i:[] for i in range(1,7)}
for r in REACTIONS:
    pages={int(e.split(":page")[1].split(":")[0]) for e in r["evidence_ids"] if ":page" in e}
    for p in pages & set(by_page):by_page[p].append(r["reaction_id"])
coverage={"document_id":DOC,"semantic_accuracy_verified":False,"pages":[]}
regions={
 1:[{"region":"title, abstract and Introduction columns","finding":"Article identity, background bioactivity and conceptual discussion; no concrete synthesis execution on this page.","evidence_ids":ev(1)}],
 2:[{"region":"Scheme 1a-c","finding":"Natural-product context plus generic annulation/retrosynthetic strategy; not current concrete experiments.","evidence_ids":ev(2)},
    {"region":"Scheme 2 and caption","finding":"Five route arrows with a-e conditions; compound-specific 1a/b→2a/b→3a/b→4a/b→5a/b→6a/b preparations mapped to records and page-6 prose where available.","evidence_ids":ev(2)},
    {"region":"Results and Discussion text","finding":"Optimization narrative and entry-r 17:83/trans-7b context; linked to Table 1 records on page 3.","evidence_ids":ev(2)}],
 3:[{"region":"Table 1 including footnote","finding":"Twenty optimization executions a-t; NMR-yield column, isolated-yield footnote for r, and 7a/7b ratios captured separately.","evidence_ids":ev(3)},
    {"region":"Scheme 3 and caption","finding":"Seven concrete synthetic arrows a-g captured; 8→9 retained as a two-stage overall yield.","evidence_ids":ev(3)},
    {"region":"left-column discussion","finding":"Deprotection, bromination, hydrolysis and natural-product outcome discussion; duplicate evidence merged with Scheme 3 records.","evidence_ids":ev(3)}],
 4:[{"region":"Scheme 4","finding":"Generic access strategy, not a concrete execution.","evidence_ids":ev(4)},
    {"region":"Table 2 including prose and footnote","finding":"Three optimization executions; ArI-loading and DMF/MeCN source conflicts retained.","evidence_ids":ev(4)},
    {"region":"Table 3 including footnote","finding":"Five aryl-iodide scope executions with combined isolated yields and 1H-NMR diastereomer ratios.","evidence_ids":ev(4)}],
 5:[{"region":"Scheme 5 and caption","finding":"Table-3a/12ba duplicate evidence merged; five downstream derivatizations b-f recorded separately.","evidence_ids":ev(5)},
    {"region":"Figures 1 and 2 with discussion","finding":"Author-proposed Pd and Ag transition-state models; classified as proposed mechanisms, not mechanistic experiments.","evidence_ids":ev(5)},
    {"region":"Conclusions","finding":"Narrative summary only; no additional experiment.","evidence_ids":ev(5)}],
 6:[{"region":"Experimental Section, left column","finding":"General information and preparations 2a, 2b, 3a, 3b, 4a, 4b; preparation facts merged with Scheme 2 records.","evidence_ids":ev(6)},
    {"region":"Experimental Section, right column","finding":"Preparations 5a, 5b, 6a and start of 6b; reagent charges/work-up retained; 6b continuation lies outside bounded pages.","evidence_ids":ev(6)},
    {"region":"NMR/IR/HRMS/TLC/optical rotation blocks","finding":"Characterization linked as source context, not additional synthesis records.","evidence_ids":ev(6)}]}
unresolved={
 1:["No concrete reaction execution is reported on this page."],
 2:["Scheme 1 structures/strategy lack execution conditions; Scheme 2 literature-referenced procedures omit full operations."],
 3:["Table 1 absolute scale/work-up/purification and ratio method are mostly unreported; Scheme 3 lacks scale/time/work-up in the bounded evidence."],
 4:["Table 2 contains source conflicts; product stereo-safe SMILES and full atom mapping remain unresolved."],
 5:["Scheme 5 omits scale/time/work-up; Figures 1-2 are untested models; product 17 stereochemistry becomes a mixture after oxidation according to prose."],
 6:["5a reports mass/mmol combinations that appear arithmetically inconsistent and were preserved verbatim; 6b procedure continues on page 7; SI was not supplied."]}
for p in range(1,7):
    coverage["pages"].append({"page":p,"state":"attempted","regions":regions[p],"record_ids":sorted(set(by_page[p])),"unresolved":unresolved[p]})
(RUN / "coverage/pages_0001_0006.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")

# If the formal merge has already occurred, preserve the read-only audit result
# beside the other raw/provenance artifacts.  Before merge this remains an audit
# of the prepared empty dataset and is therefore intentionally skipped.
current_dataset = json.loads((RUN / "dataset.json").read_text())
if current_dataset.get("experiments"):
    (RUN / "raw_model/audit.json").write_text(json.dumps(rt.dataset_audit(current_dataset), ensure_ascii=False, indent=2) + "\n")

print(json.dumps({"molecules":len(MOLECULES),"reactions":len(REACTIONS),"specific":sum(r["record_type"]=="specific_experiment" for r in REACTIONS),"generic":sum(r["record_type"]=="generic_scheme" for r in REACTIONS),"mechanisms":sum(r["record_type"]=="proposed_mechanism" for r in REACTIONS)}, ensure_ascii=False))
