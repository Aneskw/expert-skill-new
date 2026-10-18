#!/usr/bin/env python3
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_8fa26764c07f8ee4"
RUN = ROOT / "runs/current"
P = {n: [f"{DOC}:page{n}:text", f"{DOC}:page{n}:image"] for n in range(1, 5)}


def sf(raw=None, norm=None, ev=(), method="original_extraction", derivation=None,
       status=None):
    return st.fact(raw, norm, ev, method=method, derivation=derivation, status=status)


def missing(ev=()):
    return sf(ev=ev)


def molecule(mid, label, name, ev, smiles=None, *, chirality="unspecified",
             geometry="unspecified", representation=None, axial=None,
             mixture=None, assignment=None, note=None):
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = sf(label, label, ev)
    m["names"]["original"] = sf(name, name, ev)
    m["names"]["english"] = sf(name, name, ev)
    m["names"]["abbreviations"] = missing(ev)
    m["names"]["aliases"] = missing(ev)
    m["cas_numbers"] = missing(ev)
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = sf(smiles, smiles, ev, method="structure_derivation",
                         derivation="Candidate graph transcribed from the displayed structure/name; parsing is not source identity verification.")
        diag = st.stereo_inspect(smiles)
        m["stereochemistry"]["isomeric_smiles"] = sf(
            None, diag["canonical_isomeric_smiles"], ev, method="structure_derivation",
            derivation="Canonical isomeric SMILES generated from the candidate graph by the standard stereo routine.")
        m["stereochemistry"]["identity_key"] = sf(
            None, diag["stereo_graph_key"], ev, method="structure_derivation",
            derivation="Stereo graph key generated from the candidate graph; it is not a sample identity certificate.")
        if representation is None:
            representation = "partial" if "=C=" in smiles else "encoded"
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing(ev)
        m["stereochemistry"]["isomeric_smiles"] = missing(ev)
        m["stereochemistry"]["identity_key"] = missing(ev)
        representation = representation or "not_encoded"
    m["stereochemistry"]["atom_centers"] = missing(ev)
    m["stereochemistry"]["double_bonds"] = missing(ev)
    m["stereochemistry"]["axial"] = (
        sf(axial, axial, ev, method="structure_derivation",
           derivation="Allene-axis assessment is based only on the depicted local substitution pattern; source gives no axial descriptor.")
        if axial is not None else missing(ev)
    )
    m["stereochemistry"]["mixture"] = sf(mixture, mixture, ev) if mixture is not None else missing(ev)
    m["stereochemistry"]["assignment"] = sf(assignment, assignment, ev) if assignment is not None else missing(ev)
    m["stereochemistry"]["chirality_status"] = sf(chirality, chirality, ev)
    m["stereochemistry"]["geometry_status"] = sf(geometry, geometry, ev)
    m["stereochemistry"]["representation_status"] = sf(representation, representation, ev)
    m["identity_note"] = sf(note, note, ev) if note else missing(ev)
    m["external_properties"] = []
    return m


def template(mid, label, name, ev, fixed, variables):
    m = molecule(mid, label, name, ev, representation="not_encoded")
    m["representation_type"] = "generic_template"
    m["smiles"] = missing(ev)
    m["generic_template"] = {
        "fixed_structure": sf(fixed, fixed, ev),
        "variables": [{
            "label": v,
            "attachment_points": sf(f"source-drawn attachment for {v}", f"source-drawn attachment for {v}", ev),
            "definition": sf("defined by the R-group labels in the source drawing", "defined by the R-group labels in the source drawing", ev),
            "allowed_values": missing(ev),
            "constraints": missing(ev),
        } for v in variables],
        "instance_mapping": [],
        "reason": sf("Generic R-group drawing; no single concrete molecule is asserted.",
                     "Generic R-group drawing; no single concrete molecule is asserted.", ev),
    }
    return m


def participant(pid, mid, role, ev, *, equiv=None, stage="reaction", product_role=None):
    p = st.new_participant(pid, mid, role, ev, stage=stage)
    if product_role:
        p["product_role"] = sf(product_role, product_role, ev)
    if equiv is not None:
        p["charges"]["equivalents"] = sf(f"{equiv} equiv", {"number": equiv, "unit": "equiv"}, ev)
        p["equivalence_reference"] = sf("borane reagent", "m_bcf3" if mid != "m_mebcf2" else "m_mebcf2", ev)
    return p


def set_equiv(p, equiv, ref, ev):
    p["charges"]["equivalents"] = sf(f"{equiv} equiv", {"number": equiv, "unit": "equiv"}, ev)
    p["equivalence_reference"] = sf(ref, ref, ev)


def step(step_id, pids, ev, temp=None, time=None, operation="reaction"):
    s = st.new_step(step_id, ev)
    s["participant_ids"] = list(pids)
    s["operation"] = sf(operation, operation, ev)
    if temp is not None:
        tnorm = temp if isinstance(temp, dict) else temp
        traw = f"{temp['number']} {temp['unit']}" if isinstance(temp, dict) and "number" in temp else temp
        s["temperature"] = sf(traw, tnorm, ev)
    if time is not None:
        timeraw = f"{time['number']} {time['unit']}" if isinstance(time, dict) and "number" in time else time
        s["time"] = sf(timeraw, time, ev)
    return s


def outcome(oid, kind, targets, ev, raw, norm, *, yield_type=None, basis=None, scope=None):
    o = st.new_outcome(oid, kind, targets, ev)
    o["value"] = sf(raw, norm, ev)
    if yield_type:
        o["yield_type"] = sf(yield_type, yield_type, ev)
    if basis:
        o["measurement_basis"] = sf(basis, basis, ev)
    if scope:
        o["scope"] = sf(scope, scope, ev)
    return o


def classify(r, ev, changes, limitation):
    facts = []
    by_code = {}
    for i, (a, b, before, after, ctype) in enumerate(changes, 1):
        cid = f"b{i}"
        ch = {"change_id": cid, "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
              "before_order": before, "after_order": after, "change_type": ctype}
        facts.append(sf(None, ch, ev, method="inference",
                        derivation="Class-level bond change inferred from the source transformation drawing/mechanistic description; exact atom mapping is unresolved."))
        if ctype in ("formed", "hydrogen_gain", "order_increased"):
            lab = st.bond_label(a, b, after)
            by_code.setdefault(lab["code"], dict(lab, supporting_change_ids=[]))["supporting_change_ids"].append(cid)
    r["classification"]["bond_changes"] = facts
    r["classification"]["labels"] = sf(None, list(by_code.values()), ev, method="inference",
        derivation="Taxonomy labels are inferred from the class-level bond changes, not from a complete atom mapping.")
    r["classification"]["coverage"] = sf(None, "partial", ev, method="inference",
        derivation="The compact drawings and absent SI do not support complete atom mapping or all hydrogen/coproduct changes.")
    r["classification"]["limitations"] = [limitation, "Exact atom identities, full bond-order changes, hydrogen transfer, and coproduct accounting are unresolved."]
    r["classification"]["model_votes"] = []


def no_classification(r, ev, coverage="unresolved", reason="No executed transformation is atom-mapped."):
    r["classification"]["bond_changes"] = []
    r["classification"]["labels"] = sf(None, [], ev, method="inference", derivation=reason)
    r["classification"]["coverage"] = sf(None, coverage, ev, method="inference", derivation=reason)
    r["classification"]["limitations"] = [reason]
    r["classification"]["model_votes"] = []


def reaction(suffix, ev, rtype, categories, location, name, transform):
    rid = f"{DOC}:{suffix}"
    r = st.new_reaction(rid, DOC, ev)
    r["record_type"] = rtype
    r["experiment_categories"] = sf(categories, categories, ev)
    r["location"] = sf(location, location, ev)
    r["reaction_name"] = sf(name, name, ev)
    r["transformation_type"] = sf(None, transform, ev, method="inference",
        derivation="Normalized transformation description inferred from the source wording and drawing.")
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    return r


docs = json.loads((RUN / "submitted.json").read_text())["documents"]
bundle = st.empty(docs)

# Molecule registry. CAS identifiers are intentionally left unreported.
mols = [
    template("m_lit_alkyne", "R-C≡C-R'", "generic alkyne", P[1], "R-C≡C-R'", ["R", "R'"]),
    template("m_lit_carboboration", "Erker 2010 product", "generic 1,1-carboboration product", P[1], "alkenyl borane with transferred R''", ["R", "R'", "R''"]),
    template("m_arylallene_generic", "Ar-allene", "generic arylallene", P[1] + P[2], "Ar-CH=C=CH2", ["Ar", "R"]),
    template("m_indene_generic", "pentafluorophenyl indene", "generic pentafluorophenyl-substituted indene", P[1] + P[2], "R-substituted indene bearing C6F5", ["R"]),
    molecule("m_1", "1", "phenylallene", P[1] + P[2] + P[3], "C=C=Cc1ccccc1",
             axial={"allene": "not stereogenic on candidate local graph because one terminus is CH2; no source axial descriptor"},
             note="Candidate structure from the displayed phenylallene drawing."),
    molecule("m_bcf3", "B(C6F5)3", "tris(pentafluorophenyl)borane", P[1] + P[2] + P[3]),
    molecule("m_dcmd2", "DCM-d2", "dichloromethane-d2", P[1] + P[2], "[2H]C([2H])(Cl)Cl", chirality="not_applicable", geometry="not_applicable"),
    molecule("m_dcm", "DCM", "dichloromethane", P[2] + P[3], "ClCCl", chirality="not_applicable", geometry="not_applicable"),
    molecule("m_2", "2", "pentafluorophenyl-substituted indene 2", P[1] + P[2] + P[3],
             chirality="unspecified", geometry="unspecified", assignment={"method": "SCXRD", "CCDC": "2059168"},
             note="Compact fused-ring product drawing is source identified; no complete candidate SMILES is asserted."),
    molecule("m_3", "3", "allylborane side product 3", P[1] + P[2], chirality="unspecified", geometry="unspecified"),
    molecule("m_vinylborane_minor", "minor vinylborane", "minor vinylborane side product", P[1], chirality="unspecified", geometry="unspecified"),
    molecule("m_dimer_minor", "2-boryl-hexa-1,5-diene", "Piers' borane-induced allene dimerization product", P[1], chirality="unspecified", geometry="unspecified"),
    molecule("m_tmb", "1,3,5-trimethoxybenzene", "1,3,5-trimethoxybenzene", P[2], "COc1cc(OC)cc(OC)c1", chirality="not_applicable", geometry="not_applicable"),
    molecule("m_norbornene", "norbornene", "norbornene", P[1] + P[2], "C1=CC2CCC1C2", chirality="unspecified", geometry="not_applicable", representation="partial"),
    molecule("m_naph_allene", "1-naphthylallene", "1-naphthylallene", P[1] + P[2], "C=C=Cc1cccc2ccccc12",
             axial={"allene": "not stereogenic on candidate local graph because one terminus is CH2; no source axial descriptor"}),
    molecule("m_4mix", "4", "benzindene 4 positional-isomer mixture", P[1] + P[2], chirality="unspecified", geometry="mixture",
             mixture={"type": "double-bond positional isomers", "ratio": "2:1"}),
    molecule("m_ptol_allene", "p-tolylallene", "p-tolylallene", P[1] + P[2], "C=C=Cc1ccc(C)cc1",
             axial={"allene": "not stereogenic on candidate local graph because one terminus is CH2; no source axial descriptor"}),
    molecule("m_5mix", "5", "indene 5 regioisomer mixture", P[1] + P[2], chirality="unspecified", geometry="mixture",
             mixture={"type": "regioisomers", "ratio": "2:1"}),
    molecule("m_pipr_allene", "p-(isopropyl)phenylallene", "p-(isopropyl)phenylallene", P[1] + P[2], "C=C=Cc1ccc(C(C)C)cc1",
             axial={"allene": "not stereogenic on candidate local graph because one terminus is CH2; no source axial descriptor"}),
    molecule("m_6mix", "6", "indene 6 regioisomer mixture", P[1] + P[2], chirality="unspecified", geometry="mixture",
             mixture={"type": "regioisomers", "ratio": "2:1"}),
    molecule("m_pcl_allene", "p-(chloro)phenylallene", "p-(chloro)phenylallene", P[2], "C=C=Cc1ccc(Cl)cc1",
             axial={"allene": "not stereogenic on candidate local graph because one terminus is CH2; no source axial descriptor"}),
    molecule("m_7mix", "7", "indene 7 regioisomer mixture", P[2], chirality="unspecified", geometry="mixture",
             mixture={"type": "regioisomers", "ratio": "4:1"}),
    molecule("m_pf_allene", "p-(fluoro)phenylallene", "p-(fluoro)phenylallene", P[2], "C=C=Cc1ccc(F)cc1",
             axial={"allene": "not stereogenic on candidate local graph because one terminus is CH2; no source axial descriptor"}),
    molecule("m_8mix", "8", "indene 8 regioisomer mixture", P[2], chirality="unspecified", geometry="mixture",
             mixture={"type": "regioisomers", "ratio": "4:1"}),
    molecule("m_mebcf2", "MeB(C6F5)2", "methylbis(pentafluorophenyl)borane", P[2]),
    molecule("m_9mix", "9", "methylindene 9 regioisomer mixture", P[2], chirality="unspecified", geometry="mixture",
             mixture={"type": "regioisomers", "ratio": "8:1"}),
    molecule("m_hbcf2", "HB(C6F5)2", "Piers' borane", P[1] + P[2] + P[3]),
]
bundle["molecules"] = mols

reactions = []

# Cited precedent and the current generic transformation.
r = reaction("literature_erker_2010", P[1], "literature_example", ["other"], "page 1, Scheme 1, Erker 2010",
             "1,1-carboboration of alkynes", "generic alkyne 1,1-carboboration")
r["participants"] = [participant("lit_s", "m_lit_alkyne", "reactant", P[1]), participant("lit_p", "m_lit_carboboration", "product", P[1], product_role="reported precedent product")]
classify(r, P[1], [("C", "C", 0, 1, "formed"), ("C", "B", 0, 1, "formed")], "Generic cited precedent; substituents and exact atom mapping are unresolved.")
r["review"]["unresolved"] = ["Generic literature example only; experimental conditions and yields are not supplied on page 1."]
reactions.append(r)

r = reaction("generic_arylallene_indene", P[1] + P[2], "generic_scheme", [], "pages 1-2, Schemes 1 and 3 header",
             "B(C6F5)3-induced indene formation", "arylallene cyclization, aryl transfer, and retro-hydroboration")
r["participants"] = [participant("g_s", "m_arylallene_generic", "reactant", P[1] + P[2]), participant("g_b", "m_bcf3", "reagent", P[1] + P[2]), participant("g_p", "m_indene_generic", "product", P[1] + P[2], product_role="generic desired product")]
classify(r, P[1] + P[2], [("C", "C", 0, 1, "formed"), ("C", "C", 0, 1, "formed")], "Generic R-group scheme; only the two class-level C-C formations are asserted.")
r["review"]["unresolved"] = ["Generic scheme is not a concrete experiment; R-group structures and exact atom mapping remain unresolved."]
reactions.append(r)

# Initial qNMR experiment, including reported side products.
r = reaction("scheme2_initial_qnmr", P[1] + P[2], "specific_experiment", ["other"], "pages 1-2, Scheme 2 and adjacent prose",
             "initial B(C6F5)3-mediated indene formation", "indene formation with competing allene hydroboration/dimerization")
r["general_scheme_ids"] = [f"{DOC}:generic_arylallene_indene"]
parts = [participant("s1", "m_1", "reactant", P[1] + P[2]), participant("b1", "m_bcf3", "reagent", P[1] + P[2]),
         participant("sol1", "m_dcmd2", "solvent", P[1] + P[2]), participant("p2", "m_2", "product", P[1] + P[2], product_role="desired product"),
         participant("p3", "m_3", "product", P[1] + P[2], product_role="side product"),
         participant("pv", "m_vinylborane_minor", "product", P[1], product_role="minor side product"),
         participant("pd", "m_dimer_minor", "product", P[1], product_role="minor side product"),
         participant("std", "m_tmb", "reagent", P[2], stage="characterization")]
r["participants"] = parts
r["solvent_conditions"]["composition"] = sf("DCM-d2", ["m_dcmd2"], P[1] + P[2])
r["steps"] = [step("reaction_rt_45min", [p["participant_id"] for p in parts[:3]], P[1] + P[2], temp="r.t.", time={"number": 45, "unit": "min"})]
r["outcomes"] = [
    outcome("y2", "yield", ["m_2"], P[1] + P[2], "48%", {"number": 48, "unit": "%"}, yield_type="NMR yield", basis="qNMR with 1,3,5-trimethoxybenzene as internal standard"),
    outcome("y3", "yield", ["m_3"], P[2], "22%", {"number": 22, "unit": "%"}, yield_type="NMR yield", basis="qNMR with 1,3,5-trimethoxybenzene as internal standard"),
    outcome("minor_v", "trace", ["m_vinylborane_minor"], P[1], "minor amounts", "minor amounts", scope="side product from hydroboration of the internal double bond"),
    outcome("minor_d", "trace", ["m_dimer_minor"], P[1], "minor amounts", "minor amounts", scope="Piers' borane-induced allene dimerization product"),
]
classify(r, P[1] + P[2], [("C", "C", 0, 1, "formed"), ("C", "B", 0, 1, "formed"), ("C", "H", 0, 1, "hydrogen_gain")], "Multiple observed channels are present; labels are a conservative union of main cyclization/aryl transfer and the drawn allylborane hydroboration channel.")
r["review"]["unresolved"] = ["Absolute charges/equivalents, workup, and isolation details are not reported in the supplied pages.", "Minor side-product quantities and measurement methods are not reported; SI is unavailable."]
reactions.append(r)


def scope_record(suffix, substrate, product, yield_value, ev, *, temp="r.t.", allene_equiv=2,
                 norbornene=False, ratio=None, note=None):
    r = reaction(suffix, ev, "specific_experiment", ["substrate_scope"], "page 2, Scheme 3",
                 "B(C6F5)3-mediated formation of indene", "arylallene cyclization, pentafluorophenyl transfer, and retro-hydroboration")
    r["general_scheme_ids"] = [f"{DOC}:generic_arylallene_indene"]
    ps = [participant("substrate", substrate, "reactant", ev), participant("borane", "m_bcf3", "reagent", ev),
          participant("solvent", "m_dcm", "solvent", ev), participant("product", product, "product", ev, product_role="desired product")]
    set_equiv(ps[0], allene_equiv, "m_bcf3", ev)
    set_equiv(ps[1], 1, "m_bcf3", ev)
    if norbornene:
        nb = participant("norbornene", "m_norbornene", "reagent", ev)
        set_equiv(nb, 1, "m_bcf3", ev)
        ps.insert(2, nb)
    r["participants"] = ps
    r["reaction_scale"] = sf("yield reported with respect to B(C6F5)3, the limiting reagent", {"reference_molecule_id": "m_bcf3"}, ev)
    r["solvent_conditions"]["composition"] = sf("DCM", ["m_dcm"], ev)
    r["steps"] = [step("reaction_45min", [p["participant_id"] for p in ps if p["role"] != "product"], ev,
                        temp=temp, time={"number": 45, "unit": "min"})]
    r["outcomes"] = [outcome("yield", "yield", [product], ev, f"{yield_value}%", {"number": yield_value, "unit": "%"},
                               yield_type="isolated yield", basis="product isolated by column chromatography")]
    if ratio:
        r["outcomes"].append(outcome("ratio", "selectivity", [product], ev, ratio, ratio,
                                     scope="reported positional/regioisomer ratio in the displayed order"))
    classify(r, ev, [("C", "C", 0, 1, "formed"), ("C", "C", 0, 1, "formed")], "Class-level ring-closing C-C bond formation and pentafluorophenyl transfer are inferred from the source; exact atoms are unresolved.")
    r["review"]["unresolved"] = ["Absolute quantities, workup details, and complete atom mapping are not reported in the supplied pages."]
    if ratio:
        r["review"]["unresolved"].append("The ratio measurement method and absolute identities/order of mixture components are not reported.")
    if note:
        r["notes"].append(sf(note, note, ev))
    return r


reactions += [
    scope_record("scheme3_2_68", "m_1", "m_2", 68, P[1] + P[2]),
    scope_record("scheme3_2_norbornene_58", "m_1", "m_2", 58, P[1] + P[2], allene_equiv=1, norbornene=True,
                 note="Footnote a: one equivalent phenylallene and one equivalent norbornene; this is distinct from the 68% two-equivalent experiment."),
    scope_record("scheme3_4", "m_naph_allene", "m_4mix", 65, P[1] + P[2], ratio="2:1"),
    scope_record("scheme3_5_rt", "m_ptol_allene", "m_5mix", 46, P[1] + P[2], ratio="2:1"),
    scope_record("scheme3_5_0c", "m_ptol_allene", "m_5mix", 53, P[2], temp={"number": 0, "unit": "°C"}, ratio="2:1", note="Footnote b: reaction temperature 0 °C."),
    scope_record("scheme3_6_rt", "m_pipr_allene", "m_6mix", 56, P[1] + P[2], ratio="2:1"),
    scope_record("scheme3_6_0c", "m_pipr_allene", "m_6mix", 62, P[2], temp={"number": 0, "unit": "°C"}, ratio="2:1", note="Footnote b: reaction temperature 0 °C."),
    scope_record("scheme3_7_rt", "m_pcl_allene", "m_7mix", 37, P[2], ratio="4:1"),
    scope_record("scheme3_7_60c", "m_pcl_allene", "m_7mix", 45, P[2], temp={"number": 60, "unit": "°C"}, ratio="4:1", note="Footnote c: oil-bath temperature 60 °C."),
    scope_record("scheme3_8_rt", "m_pf_allene", "m_8mix", 47, P[2], ratio="4:1"),
    scope_record("scheme3_8_60c", "m_pf_allene", "m_8mix", 58, P[2], temp={"number": 60, "unit": "°C"}, ratio="4:1", note="Footnote c: oil-bath temperature 60 °C."),
]

# Hammett study.
r = reaction("figure1_hammett", P[2] + P[3], "specific_experiment", ["mechanistic_study"], "page 2, Figure 1 and pages 2-3 prose",
             "Hammett plot analysis", "kinetic substituent-effect study of arylallene cyclization")
r["general_scheme_ids"] = [f"{DOC}:generic_arylallene_indene"]
r["participants"] = [participant("allene_me", "m_ptol_allene", "reactant", P[2]), participant("allene_h", "m_1", "reactant", P[2]), participant("allene_f", "m_pf_allene", "reactant", P[2]), participant("allene_cl", "m_pcl_allene", "reactant", P[2]), participant("borane", "m_bcf3", "reagent", P[2]), participant("solvent", "m_dcm", "solvent", P[2]), participant("products", "m_indene_generic", "product", P[2], product_role="major isomers shown generically")]
for pp in r["participants"][:4]:
    set_equiv(pp, 2, "m_bcf3", P[2])
r["solvent_conditions"]["composition"] = sf("DCM", ["m_dcm"], P[2])
r["steps"] = [step("kinetic_runs", ["allene_me", "allene_h", "allene_f", "allene_cl", "borane", "solvent"], P[2], temp={"number": 0, "unit": "°C"}, operation="Hammett kinetic measurements for separate R = Me, H, F, Cl substrate runs")]
r["outcomes"] = [outcome("rho", "other", [], P[2] + P[3], "ρ = −3.9 ± 0.7", {"value": -3.9, "uncertainty": 0.7}, scope="slope of log(kX/kH) versus σp")]
classify(r, P[2] + P[3], [("C", "C", 0, 1, "formed"), ("C", "C", 0, 1, "formed")], "The rate study monitors the same indene-forming transformation; individual rate constants and exact atom mapping are not printed.")
r["review"]["unresolved"] = ["Individual kinetic protocols, concentrations, rate constants, uncertainties per point, and product yields are not reported in the supplied pages."]
reactions.append(r)

# Methyl transfer reaction.
r = reaction("scheme4_methylindene9", P[2], "specific_experiment", ["substrate_scope"], "page 2, Scheme 4",
             "formation of methylindene 9", "borane-induced cyclization with methyl transfer")
r["participants"] = [participant("substrate", "m_pipr_allene", "reactant", P[2]), participant("borane", "m_mebcf2", "reagent", P[2]), participant("solvent", "m_dcm", "solvent", P[2]), participant("product", "m_9mix", "product", P[2], product_role="regioisomer mixture")]
set_equiv(r["participants"][0], 2, "m_mebcf2", P[2]); set_equiv(r["participants"][1], 1, "m_mebcf2", P[2])
r["solvent_conditions"]["composition"] = sf("DCM", ["m_dcm"], P[2])
r["steps"] = [step("reaction_3d", ["substrate", "borane", "solvent"], P[2], temp={"number": 60, "unit": "°C"}, time={"number": 3, "unit": "d"})]
r["outcomes"] = [outcome("yield", "yield", ["m_9mix"], P[2], "30%", {"number": 30, "unit": "%"}),
                 outcome("ratio", "selectivity", ["m_9mix"], P[2], "8:1", "8:1", scope="reported regioisomer ratio in displayed order")]
classify(r, P[2], [("C", "C", 0, 1, "formed"), ("C", "C", 0, 1, "formed")], "Class-level ring closure and methyl transfer each form a C-C bond; other changes and atom identities are unresolved.")
r["review"]["unresolved"] = ["Yield measurement basis, absolute quantities, workup, purification, and ratio method are not reported."]
reactions.append(r)

# Proposed mechanism is distinct from executed experiments.
r = reaction("scheme5_proposed_mechanism", P[2] + P[3], "proposed_mechanism", [], "pages 2-3, Figure 2 and Scheme 5",
             "proposed mechanism for pentafluorophenyl transfer", "zwitterionic addition, Friedel-Crafts cyclization, aryl transfer, and retro-hydroboration")
r["participants"] = [participant("s", "m_1", "reactant", P[2] + P[3]), participant("b", "m_bcf3", "reagent", P[2] + P[3]), participant("p", "m_2", "product", P[2] + P[3], product_role="proposed overall product"), participant("hb", "m_hbcf2", "product", P[2] + P[3], product_role="proposed coproduct")]
classify(r, P[2] + P[3], [("C", "C", 0, 1, "formed"), ("C", "C", 0, 1, "formed"), ("B", "H", 0, 1, "formed")], "Author-proposed mechanism, not an independently executed synthesis; labels describe the proposed overall connectivity changes only.")
r["notes"] = [sf("INT-1 formation; intramolecular Friedel-Crafts cyclization via TS-2; intermolecular H shift; barrierless C6F5 transfer; retro-hydroboration liberating HB(C6F5)2.",
                 "INT-1 formation; intramolecular Friedel-Crafts cyclization via TS-2; intermolecular H shift; barrierless C6F5 transfer; retro-hydroboration liberating HB(C6F5)2.", P[2] + P[3])]
r["review"]["unresolved"] = ["Proposed intermediates INT-1 to INT-4 and transition states are not isolated experiments; exact atom mapping and proton-transfer partners remain unresolved."]
reactions.append(r)

bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": f"{DOC}:metadata", "document_id": DOC, "evidence_ids": P[1], "kind": "bibliographic_metadata",
     "value": {"title": "Indene formation upon borane-induced cyclization of arylallenes, 1,1-carboboration, and retro-hydroboration", "doi": "10.1039/d1cc01750k"}},
    {"context_id": f"{DOC}:structure_assignment", "document_id": DOC, "evidence_ids": P[1] + P[2], "kind": "characterization",
     "value": "Product 2 was fully characterized; SCXRD/CCDC 2059168 supports the assignment. External CIF was not supplied."},
    {"context_id": f"{DOC}:dft", "document_id": DOC, "evidence_ids": P[2] + P[3], "kind": "computational_mechanistic_context",
     "value": "DFT energy profile at revDSD-PBEP86-D4/def2-QZVPP//PBEh-3c with SMD(DCM); TS-1 is proposed rate determining."},
    {"context_id": f"{DOC}:negative_inventory", "document_id": DOC, "evidence_ids": P[1] + P[2] + P[3] + P[4], "kind": "coverage_note",
     "value": "No standalone mass-spectrum peak table, molecular-weight table, or fragment-relative-abundance block was converted into a synthesis experiment on pages 1-4; page 4 is references only."},
]

(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")

stereo = {}
for m in mols:
    smi = st.value(m["smiles"])
    if smi:
        stereo[m["molecule_id"]] = st.stereo_inspect(smi)
(RUN / "raw_model/stereo_checks.json").write_text(json.dumps(stereo, ensure_ascii=False, indent=2) + "\n")

page_inventory = {
    "document_id": DOC,
    "pages": [
        {"page": 1, "inventory": ["title/abstract/ESI notice", "Scheme 1 cited precedent and current generic scheme", "prose initial qNMR experiment and side products", "scope prose through products 4-6"], "visual_evidence_id": P[1][1]},
        {"page": 2, "inventory": ["Scheme 2 qNMR experiment", "Scheme 3 eleven condition-specific scope experiments", "Figure 1 Hammett study", "Scheme 4 methyl transfer", "DFT/mechanism prose"], "visual_evidence_id": P[2][1]},
        {"page": 3, "inventory": ["Figure 2 DFT energy profile", "Scheme 5 proposed mechanism", "summary and references"], "visual_evidence_id": P[3][1]},
        {"page": 4, "inventory": ["references only; no current-paper experiment"], "visual_evidence_id": P[4][1]},
    ],
    "deduplication": ["Scheme 2 and page-1 initial qNMR prose are one experiment.", "Scheme 3 parenthetical yields are separate condition variants under footnotes a-c, not alternative values in one experiment."],
}
(RUN / "raw_model/page_extraction.json").write_text(json.dumps(page_inventory, ensure_ascii=False, indent=2) + "\n")

ids = lambda *s: [f"{DOC}:{x}" for x in s]
coverage = {"document_id": DOC, "semantic_accuracy_verified": False, "pages": [
    {"page": 1, "state": "attempted", "regions": [
        {"region": "title, abstract, and introductory prose", "finding": "Current transformation, initial qNMR experiment, side products, scope summary, DOI, and unprovided ESI notice.", "evidence_ids": P[1]},
        {"region": "Scheme 1", "finding": "One cited generic Erker 1,1-carboboration family and one current generic arylallene-to-indene scheme.", "evidence_ids": P[1]}],
     "record_ids": ids("literature_erker_2010", "generic_arylallene_indene", "scheme2_initial_qnmr", "scheme3_2_68", "scheme3_2_norbornene_58", "scheme3_4", "scheme3_5_rt", "scheme3_6_rt"),
     "unresolved": ["SI is advertised but not supplied; preparative details and minor-product data are unavailable."]},
    {"page": 2, "state": "attempted", "regions": [
        {"region": "Scheme 2", "finding": "Phenylallene/B(C6F5)3 in DCM-d2 at room temperature for 45 min: 2 (48%) and 3 (22%) by qNMR.", "evidence_ids": P[2]},
        {"region": "Scheme 3 and footnotes a-c", "finding": "Eleven condition-specific isolated-yield records after splitting norbornene, 0 °C, and 60 °C variants; regioisomer ratios retained.", "evidence_ids": P[2]},
        {"region": "Figure 1", "finding": "Hammett mechanistic study for R = Me, H, F, Cl at 0 °C; rho = -3.9 ± 0.7.", "evidence_ids": P[2]},
        {"region": "Scheme 4", "finding": "p-(isopropyl)phenylallene + MeB(C6F5)2 in DCM, 60 °C, 3 d gives 9 in 30%, 8:1.", "evidence_ids": P[2]},
        {"region": "mechanism/DFT prose", "finding": "Computationally supported proposed mechanism; not expanded into executed intermediate experiments.", "evidence_ids": P[2]}],
     "record_ids": ids("scheme2_initial_qnmr", "scheme3_2_68", "scheme3_2_norbornene_58", "scheme3_4", "scheme3_5_rt", "scheme3_5_0c", "scheme3_6_rt", "scheme3_6_0c", "scheme3_7_rt", "scheme3_7_60c", "scheme3_8_rt", "scheme3_8_60c", "figure1_hammett", "scheme4_methylindene9", "scheme5_proposed_mechanism"),
     "unresolved": ["Scope absolute charges and ratio measurement methods are not reported; Scheme 4 yield method is not reported."]},
    {"page": 3, "state": "attempted", "regions": [
        {"region": "Figure 2", "finding": "DFT free-energy profile and optimized TS-1; computational context, not a synthesis experiment.", "evidence_ids": P[3]},
        {"region": "Scheme 5", "finding": "Author-proposed INT-2 to INT-4 sequence and pentafluorophenyl transfer; one proposed-mechanism record.", "evidence_ids": P[3]},
        {"region": "summary, acknowledgements, and references", "finding": "No additional current concrete experiment.", "evidence_ids": P[3]}],
     "record_ids": ids("figure1_hammett", "scheme5_proposed_mechanism"),
     "unresolved": ["Proposed intermediates are not isolated; exact atom mapping and proton-transfer partner are unresolved."]},
    {"page": 4, "state": "attempted", "regions": [
        {"region": "two-column reference list", "finding": "References 6-15 only; no current-paper reaction arrow, table row, procedure, or synthesis experiment.", "evidence_ids": P[4]}],
     "record_ids": [], "unresolved": ["External cited papers and crystallographic files were not supplied and were not imported."]},
]}
(RUN / "coverage").mkdir(exist_ok=True)
(RUN / "coverage/pages_0001_0004.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")

print(json.dumps({"molecules": len(mols), "reactions": len(reactions), "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in reactions)}, ensure_ascii=False))
