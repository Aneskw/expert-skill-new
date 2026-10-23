#!/usr/bin/env python3
"""Build the source-derived v3.1 submission for PDF pages 1--6."""
from collections import OrderedDict
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_35ea965f4db4cdd1"
RUN = ROOT / "runs/current"
TEMPLATE = RUN / "raw_model/submitted_template_pages_0001_0006.json"
bundle = st.read(TEMPLATE)
diagnostics = {}


def ev(page):
    return [f"{DOC}:page{page}:image", f"{DOC}:page{page}:text"]


def direct(raw, normalized, evidence):
    return st.fact(raw, normalized, evidence)


def inferred(normalized, evidence, derivation):
    return st.fact(None, normalized, evidence, method="inference", derivation=derivation)


def derived(normalized, evidence, derivation):
    return st.fact(None, normalized, evidence, method="structure_derivation", derivation=derivation)


def q(number, unit, raw=None):
    value = {"number": number, "unit": unit}
    if raw is not None:
        value["raw_text"] = raw
    return value


def molecule(mid, label, evidence, name=None, smiles=None, chirality="unresolved",
             geometry="unresolved", representation=None, assignment=None,
             atom_centers=None, mixture=None, axial=None, note=None):
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = direct(label, label, evidence)
    m["names"]["original"] = direct(name or label, name or label, evidence)
    if name:
        m["names"]["english"] = direct(name, name, evidence)
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = derived(smiles, evidence, "Candidate SMILES transcribed from the source label/drawing; toolkit parsing does not verify source identity.")
        diag = st.stereo_inspect(smiles)
        diagnostics[mid] = diag
        m["stereochemistry"]["identity_key"] = derived(diag["stereo_graph_key"], evidence, "Computed by standard stereo from the candidate isomeric graph.")
        if representation is None:
            representation = "not_applicable" if chirality == "achiral" and geometry == "not_applicable" else "partial"
        if representation == "encoded":
            m["stereochemistry"]["isomeric_smiles"] = derived(diag["canonical_isomeric_smiles"], evidence, "Canonical isomeric SMILES computed from the candidate source-derived graph.")
    if representation is None:
        representation = "not_encoded" if chirality not in ("unresolved", "achiral", "not_applicable") else "unresolved"
    m["stereochemistry"]["chirality_status"] = inferred(chirality, evidence, "Status assigned conservatively from the supplied drawing/text; absence of a label is not treated as racemic.")
    m["stereochemistry"]["geometry_status"] = inferred(geometry, evidence, "E/Z and ring stereochemistry status assessed from the supplied drawing/text.")
    m["stereochemistry"]["representation_status"] = inferred(representation, evidence, "Compares source stereochemical information with the candidate molecular representation.")
    if assignment is not None:
        m["stereochemistry"]["assignment"] = direct(assignment, assignment, evidence)
    if atom_centers is not None:
        m["stereochemistry"]["atom_centers"] = direct(atom_centers, atom_centers, evidence)
    if mixture is not None:
        m["stereochemistry"]["mixture"] = direct(mixture, mixture, evidence)
    if axial is not None:
        m["stereochemistry"]["axial"] = direct(axial, axial, evidence)
    if note is not None:
        m["identity_note"] = inferred(note, evidence, "Machine identity note explaining the conservative representation decision.")
    bundle["molecules"].append(m)
    return mid


def template_molecule(mid, label, evidence, fixed, variables, note):
    molecule(mid, label, evidence, note=note)
    m = bundle["molecules"][-1]
    m["representation_type"] = "generic_template"
    m["smiles"] = st.fact(evidence=evidence)
    m["generic_template"] = {
        "fixed_structure": derived(fixed, evidence, "Generic scaffold description derived from the source drawing; R groups are not expanded."),
        "variables": [{
            "label": v,
            "attachment_points": derived(f"source drawing attachment for {v}", evidence, "Attachment is retained only at the generic drawing level."),
            "definition": inferred("R-group as drawn; exact identity varies", evidence, "Variable definition inferred from the generic source drawing."),
            "allowed_values": st.fact(evidence=evidence),
            "constraints": st.fact(evidence=evidence),
        } for v in variables],
        "instance_mapping": [],
        "reason": inferred(note, evidence, "Machine explanation for using a generic template instead of a concrete SMILES."),
    }
    return mid


def reaction(suffix, evidence, record_type="specific_experiment", categories=None,
             location=None, reaction_name=None, transformation=None, unresolved=None):
    rid = f"{DOC}:{suffix}"
    r = st.new_reaction(rid, DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = inferred(categories, evidence, "Record category assigned from the source region and report context.") if categories is not None else st.fact(evidence=evidence)
    r["location"] = inferred(location or suffix, evidence, "Machine locator for the source region.")
    if reaction_name:
        r["reaction_name"] = direct(reaction_name, reaction_name, evidence)
    if transformation:
        r["transformation_type"] = inferred(transformation, evidence, "Transformation class inferred from the source reaction arrow and accompanying prose.")
    r["review"] = {
        "status": "needs_review", "method": "model", "reviewer": None,
        "checked_fields": [],
        "unresolved": unresolved or ["Machine extraction; all structures, stereochemical assignments, and field associations require expert review."],
    }
    bundle["reactions"].append(r)
    return r


def participant(r, pid, mid, role, evidence, amount=None, equiv=None, mol_percent=None,
                volume=None, mass=None, reference=None, stage="reaction", product_role=None):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    if amount:
        p["charges"]["amount"] = direct(amount[0], amount[1], evidence)
    if equiv:
        p["charges"]["equivalents"] = direct(equiv[0], equiv[1], evidence)
    if mol_percent:
        p["charges"]["mol_percent"] = direct(mol_percent[0], mol_percent[1], evidence)
    if volume:
        p["charges"]["volume"] = direct(volume[0], volume[1], evidence)
    if mass:
        p["charges"]["mass"] = direct(mass[0], mass[1], evidence)
    if reference:
        p["equivalence_reference"] = inferred(reference, evidence, "Equivalent/mol% denominator derived from the caption's 1.0-equivalent reference substrate.")
    if product_role:
        p["product_role"] = inferred(product_role, evidence, "Role assigned from the reaction arrow and outcome label.")
    r["participants"].append(p)
    return pid


def step(r, sid, evidence, pids, operation=None, temperature=None, time=None,
         atmosphere=None, pressure=None, stage="reaction"):
    s = st.new_step(sid, evidence, stage=stage)
    s["participant_ids"] = list(pids)
    if operation:
        s["operation"] = direct(operation, operation, evidence)
    if temperature:
        s["temperature"] = direct(temperature[0], temperature[1], evidence)
    if time:
        s["time"] = direct(time[0], time[1], evidence)
    if atmosphere:
        s["atmosphere"] = direct(atmosphere, atmosphere, evidence)
    if pressure:
        s["pressure"] = direct(pressure, pressure, evidence)
    r["steps"].append(s)
    return s


def outcome(r, oid, kind, targets, evidence, raw=None, normalized=None,
            yield_type=None, basis=None, scope=None, product_role=None):
    o = st.new_outcome(oid, kind, targets, evidence)
    if raw is not None or normalized is not None:
        o["value"] = direct(raw, normalized, evidence)
    if yield_type:
        o["yield_type"] = direct(yield_type, yield_type, evidence)
    if basis:
        o["measurement_basis"] = direct(basis, basis, evidence)
    if scope is not None:
        o["scope"] = inferred(scope, evidence, "Structured scope/metric interpretation of the source outcome.")
    if product_role:
        o["product_role"] = inferred(product_role, evidence, "Outcome role assigned from the source reaction arrow.")
    r["outcomes"].append(o)
    return o


def note(r, text, evidence, method="original_extraction"):
    if method == "original_extraction":
        r["notes"].append(direct(text, text, evidence))
    else:
        r["notes"].append(inferred(text, evidence, text))


def set_classification(r, evidence, bonds=None, coverage="partial", limitation=None):
    c = r["classification"]
    if bonds is None:
        c["coverage"] = inferred("unresolved", evidence, "The source drawing was not mapped with sufficient confidence for a bond-change assertion.")
        c["limitations"] = [limitation or "Exact atom correspondence is unresolved."]
        return
    if not bonds:
        c["labels"] = inferred([], evidence, "No observed target transformation is classified for this negative/failed experiment.")
        c["coverage"] = inferred("not_applicable", evidence, "No applicable observed bond formation was assigned.")
        c["limitations"] = [limitation or "No target bond formation observed."]
        return
    changes = []
    labels = OrderedDict()
    for idx, (a, b, change_type, before, after) in enumerate(bonds, 1):
        cid = f"b{idx}"
        change = {"change_id": cid, "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
                  "before_order": before, "after_order": after, "change_type": change_type}
        changes.append(derived(change, evidence, "Partial drawing-derived bond change; exact atom maps and full atom balance are unresolved."))
        if change_type in ("formed", "hydrogen_gain", "order_increased"):
            lab = st.bond_label(a, b, after)
            labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(cid)
    c["bond_changes"] = changes
    c["labels"] = derived(list(labels.values()), evidence, "Taxonomy labels calculated from the listed partial source-drawing bond changes.")
    c["coverage"] = inferred(coverage, evidence, "Only source-visible net formed bonds are represented; exact mapping, bond-order shifts, hydrogen balance, and coproducts remain unresolved.")
    c["limitations"] = [limitation or "Partial source-drawing classification without a complete atom-mapped reaction."]


def common_carbonylation_participants(r, evidence, substrate1, substrate2, product,
                                      scale1=None, scale2=None, pd_loading=("5 mol%", q(5, "mol%")),
                                      ligand_loading=("12 mol%", q(12, "mol%")), additive="m_agbf4",
                                      additive_equiv=("1.0 equiv.", q(1.0, "equiv")), product_role="desired product"):
    pids = []
    pids.append(participant(r, "p_substrate1", substrate1, "reactant", evidence,
                            amount=scale1, equiv=("1.0 equiv.", q(1.0, "equiv")) if scale1 else None,
                            reference=substrate1 if scale1 else None))
    pids.append(participant(r, "p_substrate2", substrate2, "reactant", evidence,
                            amount=scale2, equiv=(f"{scale2[1]['number']/scale1[1]['number']:.1f} equiv.", q(scale2[1]["number"]/scale1[1]["number"], "equiv")) if scale1 and scale2 else None,
                            reference=substrate1 if scale1 and scale2 else None))
    pids.append(participant(r, "p_pd", "m_pd2dba3", "catalyst", evidence, mol_percent=pd_loading, reference=substrate1))
    pids.append(participant(r, "p_ligand", "m_btfm_garphos", "ligand", evidence, mol_percent=ligand_loading, reference=substrate1))
    pids.append(participant(r, "p_dipea", "m_dipea", "reagent", evidence, equiv=("2.0 equiv.", q(2.0, "equiv")), reference=substrate1))
    if additive:
        pids.append(participant(r, "p_additive", additive, "additive", evidence, equiv=additive_equiv, reference=substrate1))
    pids.append(participant(r, "p_co", "m_co", "reagent", evidence))
    pids.append(participant(r, "p_dmac", "m_dmac", "solvent", evidence))
    pids.append(participant(r, "p_product", product, "product", evidence, product_role=product_role))
    return pids


def set_solvent(r, evidence, raw="DMAc (1.0 mL)", volume=1.0):
    r["solvent_conditions"]["composition"] = direct("DMAc", ["DMAc"], evidence)
    r["solvent_conditions"]["total_volume"] = direct(raw, q(volume, "mL", raw), evidence)


def er_outcome(r, oid, product, ratio, evidence, basis="chiral HPLC", scope_note=None):
    scope = {"metric": "enantiomeric ratio", "order": "major:minor", "configuration_order": "unassigned"}
    if scope_note:
        scope.update(scope_note)
    outcome(r, oid, "selectivity", [product], evidence, ratio, ratio, basis=basis, scope=scope)


P1, P2, P3, P4 = ev(1), ev(2), ev(3), ev(4)

# Core molecules and named conditions.
molecule("m_1a", "1a", P1 + P2 + P3 + P4, "2-iodoaniline", "Nc1ccccc1I", "achiral", "not_applicable", "not_applicable")
molecule("m_2a", "2a", P1 + P2 + P3 + P4, "3-phenyl-1,2-butadiene", "C=C=C(C)c1ccccc1", "achiral", "not_applicable", "not_applicable",
         note="Terminal CH2 makes this allene non-axially chiral; toolkit reports an unspecified cumulated-bond stereo possibility, so the source drawing remains authoritative.")
molecule("m_3a", "(S)-3a", P1 + P2 + P3 + P4, "(S)-3a", "C=C1N[C@@](C)(c2ccccc2)C(=O)c2ccccc21",
         "specified", "not_applicable", "encoded", assignment="absolute configuration (S), assigned by X-ray crystallographic analysis",
         atom_centers=[{"source_label": "aza-quaternary center", "configuration": "S"}],
         note="Candidate connectivity and encoded S center; source identity still requires expert structure review.")
molecule("m_3a_racemic", "racemic 3a", P2, "racemic 3a", "C=C1NC(C)(c2ccccc2)C(=O)c2ccccc21",
         "racemic", "not_applicable", "not_encoded", mixture="racemic product", note="Distinct sample identity from enantioenriched (S)-major 3a.")
molecule("m_pd2dba3", "Pd2(dba)3", P2 + P3 + P4, "tris(dibenzylideneacetone)dipalladium(0)")
molecule("m_btfm_garphos", "(R)-BTFM-GarPhos", P2 + P3 + P4, "(R)-BTFM-GarPhos", chirality="specified", geometry="not_applicable", representation="not_encoded", axial="(R)-configured axially chiral ligand")
molecule("m_r_binap", "(R)-BINAP", P2, "(R)-BINAP", chirality="specified", geometry="not_applicable", representation="not_encoded", axial="(R)-configured axially chiral ligand")
molecule("m_dipea", "iPr2NEt", P2 + P3 + P4, "N,N-diisopropylethylamine", "CCN(C(C)C)C(C)C", "achiral", "not_applicable", "not_applicable")
molecule("m_agbf4", "AgBF4", P2 + P3 + P4, "silver tetrafluoroborate", "[Ag+].[B-](F)(F)(F)F", "achiral", "not_applicable", "not_applicable")
molecule("m_agotf", "AgOTf", P2 + P3, "silver trifluoromethanesulfonate", "[Ag+].[O-]S(=O)(=O)C(F)(F)F", "achiral", "not_applicable", "not_applicable")
molecule("m_dmac", "DMAc", P2 + P3 + P4, "N,N-dimethylacetamide", "CC(=O)N(C)C", "achiral", "not_applicable", "not_applicable")
molecule("m_co", "CO", P1 + P2 + P3 + P4, "carbon monoxide", "[C-]#[O+]", "achiral", "not_applicable", "not_applicable")
molecule("m_decane", "decane", P3, "decane", "CCCCCCCCCC", "achiral", "not_applicable", "not_applicable")
template_molecule("m_petroleum_ether", "petroleum ether", P3 + P4, "petroleum ether hydrocarbon mixture", ["boiling_range/composition"], "Named solvent mixture; exact composition is not reported.")
molecule("m_dcm", "dichloromethane", P3 + P4, "dichloromethane", "ClCCl", "achiral", "not_applicable", "not_applicable")

# Scheme 1: four prior-work arrows and one generic arrow for this work.
for suffix, loc, trans in [
    ("scheme1B_aerobic_oxidation_ref7", "Scheme 1B, aerobic oxidation, ref. 7", "aerobic C-N bond formation"),
    ("scheme1B_radical_bicyclization_ref8", "Scheme 1B, intramolecular radical bicyclization, ref. 8", "intramolecular radical bicyclization"),
    ("scheme1B_cycloaddition_ref9", "Scheme 1B, [4+2] cycloaddition, ref. 9", "[4+2] cycloaddition"),
    ("scheme1C_prior_cyclocarbonylation_ref12", "Scheme 1C, Pd-catalyzed non-asymmetric cyclocarbonylation, ref. 12", "non-asymmetric cyclocarbonylation"),
]:
    rin = template_molecule(f"m_{suffix}_input", f"{loc} input(s)", P1, "generic structures exactly as drawn in Scheme 1", ["R1", "R2", "R3", "R4", "Ar", "X"], "Prior-work generic structure; variables are not expanded.")
    rout = template_molecule(f"m_{suffix}_product", f"{loc} product", P1, "generic product exactly as drawn in Scheme 1", ["R1", "R2", "R3", "R4", "Ar", "X"], "Prior-work generic product; variables are not expanded.")
    r = reaction(suffix, P1, "literature_example", None, loc, transformation=trans,
                 unresolved=["Prior-work illustration only; detailed experimental conditions and exact instances are not supplied on this page."])
    participant(r, "p_input", rin, "reactant", P1)
    participant(r, "p_product", rout, "product", P1, product_role="illustrated literature product")
    step(r, "illustrated_arrow", P1, ["p_input", "p_product"], operation=loc)
    set_classification(r, P1, None, limitation="Generic prior-work drawing with variable groups; no complete atom mapping asserted.")

generic_1 = template_molecule("m_scheme1d_generic_iodoaniline", "generic o-iodoaniline", P1, "o-iodoaniline scaffold", ["R"], "Scheme 1D generic substrate.")
generic_2 = template_molecule("m_scheme1d_generic_allene", "generic aryl allene", P1, "aryl-substituted allene scaffold", ["Ar"], "Scheme 1D generic substrate.")
generic_3 = template_molecule("m_scheme1d_generic_product", "generic dihydroquinolinone product", P1, "methylene-2,3-dihydroquinolin-4-one scaffold", ["R", "Ar"], "Scheme 1D generic product.")
r = reaction("scheme1D_this_work_general", P1, "generic_scheme", None, "Scheme 1D", "carbonylative construction of chiral aza-quaternary carbon center via C-N bond formation", "asymmetric carbonylative Heck annulation", ["Generic summary; reported maxima are not assigned to a single substrate pair."])
for pid, mid, role in [("p_1", generic_1, "reactant"), ("p_2", generic_2, "reactant"), ("p_co", "m_co", "reagent"), ("p_prod", generic_3, "product")]:
    participant(r, pid, mid, role, P1, product_role="generic product" if role == "product" else None)
step(r, "general_arrow", P1, [p["participant_id"] for p in r["participants"]], operation="[Pd]/L*, CO balloon, 80 °C, 6 h",
     temperature=("80 °C", q(80, "°C")), time=("6 h", q(6, "h")), atmosphere="CO balloon")
outcome(r, "general_yield_range", "other", [generic_3], P1, "up to 91% yield", "up to 91% yield", scope={"type": "generic maximum; not a single experiment"})
outcome(r, "general_er", "selectivity", [generic_3], P1, "95:5 e.r. (up to 99:1 e.r. after crystallization)", "95:5 e.r. (up to 99:1 e.r. after crystallization)", scope={"metric": "generic enantiomeric-ratio maxima"})
set_classification(r, P1, [("C", "C", "formed", 0, 1), ("C", "C", "formed", 0, 1), ("C", "N", "formed", 0, 1)])

# Page 2 prose-only partial experiments and controls.
for suffix, ligand, yld, er in [
    ("ligand_screen_r_binap", "m_r_binap", 88, None),
    ("ligand_screen_btfm_garphos", "m_btfm_garphos", None, "86:14"),
]:
    r = reaction(suffix, P2, categories=["optimization"], location="page 2 ligand-screen discussion",
                 reaction_name="palladium-catalysed enantioselective carbonylative Heck reaction",
                 transformation="asymmetric carbonylative Heck annulation",
                 unresolved=["The detailed SI Table S1 row, substrate charges, catalyst/ligand loadings, and any unquoted outcome are unavailable."])
    for pid, mid, role in [("p_1a", "m_1a", "reactant"), ("p_2a", "m_2a", "reactant"), ("p_pd", "m_pd2dba3", "catalyst"), ("p_ligand", ligand, "ligand"), ("p_base", "m_dipea", "reagent"), ("p_co", "m_co", "reagent"), ("p_dmac", "m_dmac", "solvent"), ("p_3a", "m_3a", "product")]:
        participant(r, pid, mid, role, P2, product_role="desired product" if role == "product" else None)
    step(r, "screening_reaction", P2, [p["participant_id"] for p in r["participants"]], operation="ligand screening in DMAc under CO",
         temperature=("80 °C", q(80, "°C")), time=("24 h", q(24, "h")), atmosphere="CO atmosphere")
    if yld is not None:
        outcome(r, "reported_yield", "yield", ["m_3a"], P2, f"{yld}%", q(yld, "%"), product_role="desired product")
    if er is not None:
        er_outcome(r, "reported_er", "m_3a", er, P2, basis=None)
    set_classification(r, P2, [("C", "C", "formed", 0, 1), ("C", "C", "formed", 0, 1), ("C", "N", "formed", 0, 1)])

molecule("m_2_bromoaniline", "2-bromoaniline", P2, "2-bromoaniline", "Nc1ccccc1Br", "achiral", "not_applicable", "not_applicable")
molecule("m_2_aminophenyl_triflate", "2-amino phenyltriflate", P2, "2-aminophenyl triflate", "Nc1ccccc1OS(=O)(=O)C(F)(F)F", "achiral", "not_applicable", "not_applicable")
for suffix, sub in [("control_2_bromoaniline", "m_2_bromoaniline"), ("control_2_aminophenyl_triflate", "m_2_aminophenyl_triflate")]:
    r = reaction(suffix, P2, categories=["control"], location="page 2 negative-substrate discussion", transformation="attempted carbonylative annulation",
                 unresolved=["The detailed Table S7 conditions and detection method are not supplied; no expected product is invented."])
    participant(r, "p_substrate", sub, "reactant", P2)
    participant(r, "p_allene", "m_2a", "reactant", P2)
    participant(r, "p_co", "m_co", "reagent", P2)
    step(r, "reported_control", P2, ["p_substrate", "p_allene", "p_co"], operation="o-iodoaniline was replaced by the indicated less reactive substrate", atmosphere="CO atmosphere")
    outcome(r, "no_conversion", "conversion", [sub], P2, "no conversion", "no conversion", scope={"target": "input substrate"})
    set_classification(r, P2, [], limitation="Source explicitly reports no conversion; no expected-product bond labels are assigned.")

molecule("m_pdi2", "PdI2", P2, "palladium(II) iodide", "[Pd+2].[I-].[I-]", "achiral", "not_applicable", "not_applicable")
molecule("m_pd_c", "Pd/C", P2, "palladium on carbon", note="Named heterogeneous catalyst; composition/loading is not given.")
for suffix, catalyst in [("ligand_free_pd2dba3", "m_pd2dba3"), ("ligand_free_pdi2", "m_pdi2"), ("ligand_free_pd_c", "m_pd_c")]:
    r = reaction(suffix, P2, categories=["control", "optimization"], location="page 2 ligand-free catalyst discussion", transformation="ligand-free carbonylative annulation",
                 unresolved=["Detailed Table S6 quantities, temperature/time, quantitative conversion value, and trace-product measurement basis are unavailable."])
    for pid, mid, role in [("p_1a", "m_1a", "reactant"), ("p_2a", "m_2a", "reactant"), ("p_cat", catalyst, "catalyst"), ("p_co", "m_co", "reagent"), ("p_3a", "m_3a_racemic", "product")]:
        participant(r, pid, mid, role, P2, product_role="minor racemic product" if role == "product" else None)
    step(r, "ligand_free_control", P2, [p["participant_id"] for p in r["participants"]], operation="reaction without phosphine ligand", atmosphere="CO atmosphere")
    outcome(r, "substrate_conversion", "conversion", ["m_1a"], P2, "quantitatively consumed", "quantitatively consumed", scope={"target": "2-iodoaniline"})
    outcome(r, "racemic_3a_trace", "trace", ["m_3a_racemic"], P2, "only a little racemic product 3a", "only a little racemic product 3a", scope={"qualitative": True})
    set_classification(r, P2, [("C", "C", "formed", 0, 1), ("C", "C", "formed", 0, 1), ("C", "N", "formed", 0, 1)], limitation="Labels describe the observed minor racemic 3a formation, not the uncharacterized consumption pathways.")

# Table 1 entries. Entry 4 also represents Scheme 2 product 3a and prose on pages 2--3.
table_rows = [
    (1, None, "one portion", 0.4, 34, "90:10", None),
    (2, None, "syringe 6 h", 0.3, 50, "89:11", None),
    (3, "m_agbf4", "one portion", 0.4, 29, "92:8", None),
    (4, "m_agbf4", "syringe 6 h", 0.3, 96, "89:11", {"isolated": 88, "recrystallized": 66, "recryst_er": "99:1"}),
    (5, "m_agbf4", "syringe 2 h + hold 4 h", 0.3, 90, "88:12", None),
    (6, "m_agotf", "syringe 6 h", 0.3, 95, "85:15", None),
]
for n, additive, mode, mmol2, yld, er, extra in table_rows:
    record_evidence = P2 + P3 if n == 4 else P3
    evidence = P3
    r = reaction(f"table1_entry_{n}", record_evidence, categories=["optimization"], location=f"Table 1 entry {n}",
                 reaction_name="palladium-catalysed enantioselective carbonylative Heck reaction", transformation="asymmetric carbonylative Heck annulation",
                 unresolved=["No SI procedure is associated; workup/purification details beyond the stated recrystallization are unavailable."])
    pids = common_carbonylation_participants(r, evidence, "m_1a", "m_2a", "m_3a",
        scale1=("0.2 mmol", q(0.2, "mmol")), scale2=(f"{mmol2} mmol", q(mmol2, "mmol")), additive=additive)
    # Table footnote gives absolute catalyst/ligand/base/additive charges.
    for p in r["participants"]:
        if p["participant_id"] == "p_pd": p["charges"]["amount"] = direct("0.01 mmol", q(0.01, "mmol"), evidence)
        if p["participant_id"] == "p_ligand": p["charges"]["amount"] = direct("0.024 mmol", q(0.024, "mmol"), evidence)
        if p["participant_id"] == "p_dipea": p["charges"]["amount"] = direct("0.4 mmol", q(0.4, "mmol"), evidence)
        if p["participant_id"] == "p_additive": p["charges"]["amount"] = direct("0.2 mmol", q(0.2, "mmol"), evidence)
        if p["participant_id"] == "p_dmac": p["charges"]["volume"] = direct("DMAc (1.0 mL)", q(1.0, "mL"), evidence)
    set_solvent(r, evidence)
    r["reaction_scale"] = direct("1a (0.2 mmol)", {"reference_molecule_id": "m_1a", "amount": q(0.2, "mmol")}, evidence)
    r["procedure"]["addition_order"] = direct(mode, mode, evidence)
    if mode == "one portion":
        step(r, "reaction", evidence, pids, operation="one-portion addition; react for 6 h", temperature=("80 °C", q(80, "°C")), time=("6 h", q(6, "h")), atmosphere="CO balloon")
    elif mode == "syringe 2 h + hold 4 h":
        step(r, "syringe_addition", evidence, pids, operation="1a in 0.5 mL DMAc fed into 0.5 mL reaction solution by programmable syringe pump", temperature=("80 °C", q(80, "°C")), time=("over 2 h", q(2, "h")), atmosphere="CO", pressure="atmospheric CO pressure")
        step(r, "post_addition_hold", evidence, pids, operation="then further reacted", temperature=("80 °C", q(80, "°C")), time=("4 h", q(4, "h")), atmosphere="CO", pressure="atmospheric CO pressure")
    else:
        step(r, "syringe_addition_reaction", evidence, pids, operation="1a in 0.5 mL DMAc fed into 0.5 mL reaction solution by programmable syringe pump", temperature=("80 °C", q(80, "°C")), time=("over 6 h", q(6, "h")), atmosphere="CO", pressure="atmospheric CO pressure")
    outcome(r, "gc_yield", "yield", ["m_3a"], evidence, f"{yld}%", q(yld, "%"), basis="GC using decane as internal standard", scope={"stage": "reaction mixture"}, product_role="desired product")
    er_outcome(r, "er", "m_3a", er, evidence)
    if extra:
        outcome(r, "isolated_yield", "yield", ["m_3a"], evidence, f"{extra['isolated']}%", q(extra["isolated"], "%"), yield_type="isolated", basis="isolated yield", scope={"stage": "before recrystallization"}, product_role="desired product")
        participant(r, "p_petroleum_ether", "m_petroleum_ether", "purification_agent", evidence, stage="purification")
        participant(r, "p_dcm_pur", "m_dcm", "purification_agent", evidence, stage="purification")
        step(r, "recrystallization", evidence, ["p_product", "p_petroleum_ether", "p_dcm_pur"], operation="recrystallization with petroleum ether/dichloromethane", stage="purification")
        outcome(r, "yield_after_recrystallization", "yield", ["m_3a"], evidence, f"{extra['recrystallized']}%", q(extra["recrystallized"], "%"), yield_type="after recrystallization", basis="isolated after recrystallization", scope={"starting_point": "isolated 3a", "end_point": "recrystallized 3a"}, product_role="recrystallized product")
        er_outcome(r, "er_after_recrystallization", "m_3a", extra["recryst_er"], evidence, scope_note={"stage": "after recrystallization", "configuration_order": "S-major:minor"})
        r["general_scheme_ids"] = [f"{DOC}:scheme2_general"]
        note(r, "Page 2 discusses the syringe-pump/AgBF4 benefit; the numerical charges and outcomes are taken from Table 1 and Scheme 2 on page 3.", P2, method="inference")
    set_classification(r, evidence, [("C", "C", "formed", 0, 1), ("C", "C", "formed", 0, 1), ("C", "N", "formed", 0, 1)])

# Scheme 2 generic arrow and scope entries 3b--3u; 3a is Table 1 entry 4.
sg1 = template_molecule("m_scheme2_generic_1", "1", P3, "substituted o-iodoaniline scaffold", ["R1", "R2"], "Scheme 2 generic substrate family.")
sg2 = template_molecule("m_scheme2_generic_2", "2", P3, "substituted allene scaffold", ["Ar"], "Scheme 2 generic allene family.")
sg3 = template_molecule("m_scheme2_generic_3", "3", P3, "substituted methylene-dihydroquinolinone scaffold", ["R1", "R2", "Ar"], "Scheme 2 generic product family.")
r = reaction("scheme2_general", P3, "generic_scheme", None, "Scheme 2 general arrow", "substrate scope of o-iodoanilines and allenes", "asymmetric carbonylative Heck annulation", ["Generic substrate-family arrow; individual products are separate records."])
for pid, mid, role in [("p_1", sg1, "reactant"), ("p_2", sg2, "reactant"), ("p_co", "m_co", "reagent"), ("p_3", sg3, "product")]:
    participant(r, pid, mid, role, P3, product_role="generic product" if role == "product" else None)
step(r, "generic_scope_arrow", P3, [p["participant_id"] for p in r["participants"]], operation="Pd2(dba)3/(R)-BTFM-GarPhos, iPr2NEt, AgBF4, CO balloon, DMAc", temperature=("80 °C", q(80, "°C")), time=("6 h", q(6, "h")), atmosphere="CO balloon")
set_classification(r, P3, [("C", "C", "formed", 0, 1), ("C", "C", "formed", 0, 1), ("C", "N", "formed", 0, 1)])

scope_rows = [
    ("3b", 90, "80:20", None), ("3c", 87, "88:12", (43, "99:1")),
    ("3d", 80, "87:13", None), ("3e", 81, "88:12", None),
    ("3f", 85, "90:10", (62, "99.5:0.5")), ("3g", 85, "86:14", None),
    ("3h", 85, "82:18", None), ("3i", 47, "88:12", None),
    ("3j", 38, "90:10", None), ("3k", 82, "82:18", None),
    ("3l", 71, "73:27", None), ("3m", 52, "88:12", None),
    ("3n", 62, "91:9", None), ("3o", 89, "81:19", None),
    ("3p", 91, "79:21", None), ("3q", 74, "95:5", None),
    ("3r", 80, "79:21", None), ("3s", 15, "85:15", None),
]
for label, yld, ratio, recryst in scope_rows:
    product = molecule(f"m_{label}", label, P3, label, chirality="unspecified", geometry="not_applicable", representation="not_encoded",
                       atom_centers=[{"source": "solid-wedge product drawing", "absolute_configuration": "not assigned"}],
                       note="Source drawing reviewed; exact connectivity/substituent position and absolute configuration are not encoded as SMILES pending expert structure review.")
    if label <= "3i":
        sub1 = molecule(f"m_1_for_{label}", f"o-iodoaniline substrate for {label}", P3, chirality="unresolved", geometry="not_applicable", representation="not_encoded", note="Distinct source-drawn ring-substituted o-iodoaniline; exact SMILES not asserted.")
        sub2 = "m_2a"
    elif label <= "3r":
        sub1 = "m_1a"
        sub2 = molecule(f"m_2_for_{label}", f"allene substrate for {label}", P3, chirality="unresolved", geometry="unresolved", representation="not_encoded", note="Distinct source-drawn substituted allene; exact SMILES and axial-stereo applicability are unresolved.")
    else:
        sub1 = molecule(f"m_1_for_{label}", f"N-substituted o-iodoaniline substrate for {label}", P3,
                        name="N-methyl-2-iodoaniline", smiles="CNc1ccccc1I", chirality="achiral",
                        geometry="not_applicable", representation="not_applicable",
                        note="Candidate source-derived connectivity for the N-methyl substrate used for 3s.")
        sub2 = "m_2a"
    r = reaction(f"scheme2_{label}", P3, categories=["substrate_scope"], location=f"Scheme 2 product {label}", reaction_name="carbonylative Heck reaction between o-iodoanilines and allenes", transformation="asymmetric carbonylative Heck annulation",
                 unresolved=["No associated SI procedure/characterization is supplied; exact source-drawn connectivity and absolute configuration require expert review."])
    pids = common_carbonylation_participants(r, P3, sub1, sub2, product,
        scale1=("0.2 mmol", q(0.2, "mmol")), scale2=("0.3 mmol", q(0.3, "mmol")))
    for p in r["participants"]:
        if p["participant_id"] == "p_pd": p["charges"]["amount"] = direct("0.01 mmol", q(0.01, "mmol"), P3)
        if p["participant_id"] == "p_ligand": p["charges"]["amount"] = direct("0.024 mmol", q(0.024, "mmol"), P3)
        if p["participant_id"] == "p_dipea": p["charges"]["amount"] = direct("0.4 mmol", q(0.4, "mmol"), P3)
        if p["participant_id"] == "p_additive": p["charges"]["amount"] = direct("0.2 mmol", q(0.2, "mmol"), P3)
        if p["participant_id"] == "p_dmac": p["charges"]["volume"] = direct("DMAc (1.0 mL)", q(1.0, "mL"), P3)
    set_solvent(r, P3)
    r["reaction_scale"] = direct("1 (0.2 mmol)", {"reference_molecule_id": sub1, "amount": q(0.2, "mmol")}, P3)
    r["procedure"]["addition_order"] = direct("1 was dissolved in 0.5 mL DMAc and fed into the reaction solution (0.5 mL DMAc) with a programmable syringe pump over 6 h", "substrate 1 syringe-fed over 6 h; solvent split 0.5 mL + 0.5 mL", P3)
    step(r, "syringe_addition_reaction", P3, pids, operation="substrate 1 in 0.5 mL DMAc fed into 0.5 mL reaction solution by programmable syringe pump", temperature=("80 °C", q(80, "°C")), time=("over 6 h", q(6, "h")), atmosphere="CO balloon")
    outcome(r, "isolated_yield", "yield", [product], P3, f"{yld}% yield", q(yld, "%"), yield_type="isolated", basis="isolated yield", product_role="desired product")
    er_outcome(r, "er", product, ratio, P3)
    if recryst:
        participant(r, "p_petroleum_ether", "m_petroleum_ether", "purification_agent", P3, stage="purification")
        participant(r, "p_dcm_pur", "m_dcm", "purification_agent", P3, stage="purification")
        step(r, "recrystallization", P3, ["p_product", "p_petroleum_ether", "p_dcm_pur"], operation="recrystallization with petroleum ether/dichloromethane", stage="purification")
        outcome(r, "yield_after_recrystallization", "yield", [product], P3, f"{recryst[0]}%", q(recryst[0], "%"), yield_type="after recrystallization", basis="isolated after recrystallization", scope={"starting_point": f"isolated {label}", "end_point": f"recrystallized {label}"}, product_role="recrystallized product")
        er_outcome(r, "er_after_recrystallization", product, recryst[1], P3, scope_note={"stage": "after recrystallization"})
    r["general_scheme_ids"] = [f"{DOC}:scheme2_general"]
    set_classification(r, P3, [("C", "C", "formed", 0, 1), ("C", "C", "formed", 0, 1), ("C", "N", "formed", 0, 1)])

# Failed Scheme 2 entries 3t and 3u.
molecule("m_3t", "3t", P3, "desired 3t", chirality="unspecified", geometry="not_applicable", representation="not_encoded", note="Desired product shown but reported n.d.")
molecule("m_1_for_3t", "N-acetyl o-iodoaniline", P3, "N-acetyl o-iodoaniline", "CC(=O)Nc1ccccc1I", "achiral", "not_applicable", "not_applicable")
molecule("m_benzoxazinone", "2-methyl-4H-3,1-benzoxazin-4-one", P3, "2-methyl-4H-3,1-benzoxazin-4-one", "CC1=NC2=CC=CC=C2OC1=O", "achiral", "not_applicable", "not_applicable")
r = reaction("scheme2_3t", P3, categories=["substrate_scope"], location="Scheme 2 desired product 3t", transformation="competing intramolecular cyclocarbonylation",
             unresolved=["Desired 3t detection method and the alternative-product yield/measurement basis are unreported; exact bond mapping is unresolved."])
pids = common_carbonylation_participants(r, P3, "m_1_for_3t", "m_2a", "m_3t", scale1=("0.2 mmol", q(0.2, "mmol")), scale2=("0.3 mmol", q(0.3, "mmol")))
for p in r["participants"]:
    if p["participant_id"] == "p_pd": p["charges"]["amount"] = direct("0.01 mmol", q(0.01, "mmol"), P3)
    if p["participant_id"] == "p_ligand": p["charges"]["amount"] = direct("0.024 mmol", q(0.024, "mmol"), P3)
    if p["participant_id"] == "p_dipea": p["charges"]["amount"] = direct("0.4 mmol", q(0.4, "mmol"), P3)
    if p["participant_id"] == "p_additive": p["charges"]["amount"] = direct("0.2 mmol", q(0.2, "mmol"), P3)
    if p["participant_id"] == "p_dmac": p["charges"]["volume"] = direct("DMAc (1.0 mL)", q(1.0, "mL"), P3)
participant(r, "p_alt_product", "m_benzoxazinone", "product", P3, product_role="alternative product")
set_solvent(r, P3)
r["reaction_scale"] = direct("1 (0.2 mmol)", {"reference_molecule_id": "m_1_for_3t", "amount": q(0.2, "mmol")}, P3)
r["procedure"]["addition_order"] = direct("1 was dissolved in 0.5 mL DMAc and fed into the reaction solution (0.5 mL DMAc) with a programmable syringe pump over 6 h", "substrate 1 syringe-fed over 6 h; solvent split 0.5 mL + 0.5 mL", P3)
step(r, "scope_attempt", P3, [p["participant_id"] for p in r["participants"]], operation="Scheme 2 standard syringe-pump conditions", temperature=("80 °C", q(80, "°C")), time=("6 h", q(6, "h")), atmosphere="CO balloon")
outcome(r, "desired_product_not_detected", "no_reaction", ["m_3t"], P3, "n.d.", "not detected", scope={"target": "desired 3t"})
outcome(r, "alternative_product_formed", "other", ["m_benzoxazinone"], P3, "2-methyl-4H-3,1-benzoxazin-4-one instead", "formed; amount unreported", scope={"product": "alternative"})
r["general_scheme_ids"] = [f"{DOC}:scheme2_general"]
set_classification(r, P3, None, limitation="Alternative cyclocarbonylation product is identified, but exact source-to-product atom correspondence is not mapped.")

molecule("m_3u", "3u", P3, "desired 3u", chirality="unspecified", geometry="not_applicable", representation="not_encoded", note="Desired product shown but reported n.d.")
molecule("m_1_for_3u", "N-tosyl o-iodoaniline", P3, "N-tosyl o-iodoaniline", "Cc1ccc(S(=O)(=O)Nc2ccccc2I)cc1", "achiral", "not_applicable", "not_applicable")
molecule("m_3u_side_products", "inseparable side products from N-tosyl o-iodoaniline", P3, note="Unknown-composition side-product mixture; no unique structure is invented.")
r = reaction("scheme2_3u", P3, categories=["substrate_scope"], location="Scheme 2 desired product 3u", transformation="failed target annulation with side-product formation",
             unresolved=["Desired-product detection method and identities/amounts of inseparable side products are unreported."])
pids = common_carbonylation_participants(r, P3, "m_1_for_3u", "m_2a", "m_3u", scale1=("0.2 mmol", q(0.2, "mmol")), scale2=("0.3 mmol", q(0.3, "mmol")))
for p in r["participants"]:
    if p["participant_id"] == "p_pd": p["charges"]["amount"] = direct("0.01 mmol", q(0.01, "mmol"), P3)
    if p["participant_id"] == "p_ligand": p["charges"]["amount"] = direct("0.024 mmol", q(0.024, "mmol"), P3)
    if p["participant_id"] == "p_dipea": p["charges"]["amount"] = direct("0.4 mmol", q(0.4, "mmol"), P3)
    if p["participant_id"] == "p_additive": p["charges"]["amount"] = direct("0.2 mmol", q(0.2, "mmol"), P3)
    if p["participant_id"] == "p_dmac": p["charges"]["volume"] = direct("DMAc (1.0 mL)", q(1.0, "mL"), P3)
participant(r, "p_side_products", "m_3u_side_products", "product", P3, product_role="inseparable side products")
set_solvent(r, P3)
r["reaction_scale"] = direct("1 (0.2 mmol)", {"reference_molecule_id": "m_1_for_3u", "amount": q(0.2, "mmol")}, P3)
r["procedure"]["addition_order"] = direct("1 was dissolved in 0.5 mL DMAc and fed into the reaction solution (0.5 mL DMAc) with a programmable syringe pump over 6 h", "substrate 1 syringe-fed over 6 h; solvent split 0.5 mL + 0.5 mL", P3)
step(r, "scope_attempt", P3, [p["participant_id"] for p in r["participants"]], operation="Scheme 2 standard syringe-pump conditions", temperature=("80 °C", q(80, "°C")), time=("6 h", q(6, "h")), atmosphere="CO balloon")
outcome(r, "desired_product_not_detected", "no_reaction", ["m_3u"], P3, "n.d.", "not detected", scope={"target": "desired 3u"})
outcome(r, "side_products", "other", ["m_3u_side_products"], P3, "a series of inseparable side products", "inseparable side products", scope={"composition": "unreported"})
r["general_scheme_ids"] = [f"{DOC}:scheme2_general"]
set_classification(r, P3, None, limitation="Unidentified side-product mixture prevents meaningful bond classification.")

# Scheme 3 scale-up and three downstream transformations.
r = reaction("scheme3_scale_up", P3 + P4, categories=["scale_up"], location="Scheme 3 scale-up reaction", reaction_name="scale-up carbonylative Heck reaction", transformation="asymmetric carbonylative Heck annulation",
             unresolved=["The 86% and 60% yield measurement bases are not explicitly stated in Scheme 3; the detailed workup is unavailable without SI."])
pids = common_carbonylation_participants(r, P3 + P4, "m_1a", "m_2a", "m_3a", scale1=("1.0 mmol", q(1.0, "mmol")), scale2=("1.2 mmol", q(1.2, "mmol")), pd_loading=("2.5 mol%", q(2.5, "mol%")), ligand_loading=("7.5 mol%", q(7.5, "mol%")))
set_solvent(r, P4, raw="DMAc; volume not reported in Scheme 3", volume=0)
r["solvent_conditions"]["total_volume"] = st.fact(evidence=P4)
r["reaction_scale"] = direct("1a (1.0 mmol)", {"reference_molecule_id": "m_1a", "amount": q(1.0, "mmol")}, P4)
step(r, "scale_up_reaction", P4, pids, operation="scale-up under Scheme 3 conditions", temperature=("80 °C", q(80, "°C")), time=("6 h", q(6, "h")), atmosphere="CO balloon")
outcome(r, "reported_yield", "yield", ["m_3a"], P4, "86%", q(86, "%"), product_role="3a before recrystallization")
er_outcome(r, "er", "m_3a", "89:11", P4, basis=None)
participant(r, "p_petroleum_ether", "m_petroleum_ether", "purification_agent", P3 + P4, stage="purification")
participant(r, "p_dcm_pur", "m_dcm", "purification_agent", P3 + P4, stage="purification")
step(r, "recrystallization", P3 + P4, ["p_product", "p_petroleum_ether", "p_dcm_pur"], operation="recrystallization with petroleum ether/dichloromethane", stage="purification")
outcome(r, "yield_after_recrystallization", "yield", ["m_3a"], P3 + P4, "60%", q(60, "%"), yield_type="after recrystallization", scope={"starting_point": "1a/2a scale-up", "end_point": "recrystallized 3a"}, product_role="recrystallized product")
er_outcome(r, "er_after_recrystallization", "m_3a", "99:1", P3 + P4, basis=None, scope_note={"stage": "after recrystallization", "configuration_order": "S-major:minor"})
set_classification(r, P4, [("C", "C", "formed", 0, 1), ("C", "C", "formed", 0, 1), ("C", "N", "formed", 0, 1)])

for mid, label, name, smi in [
    ("m_ethyl_diazoacetate", "ethyl 2-diazoacetate", "ethyl 2-diazoacetate", "CCOC(=O)[CH-][N+]#N"),
    ("m_nitrosobenzene", "nitrosobenzene", "nitrosobenzene", "O=Nc1ccccc1"),
    ("m_koac", "KOAc", "potassium acetate", "[K+].CC(=O)[O-]"),
    ("m_dce", "DCE", "1,2-dichloroethane", "ClCCCl"),
    ("m_p_methylthiophenol", "para-methylthiophenol", "4-methylbenzenethiol", "Cc1ccc(S)cc1"),
    ("m_et3n", "Et3N", "triethylamine", "CCN(CC)CC"),
    ("m_toluene", "toluene", "toluene", "Cc1ccccc1"),
    ("m_indole", "indole", "indole", "c1ccc2[nH]ccc2c1"),
    ("m_diphenyl_phosphate", "(PhO)2P(O)OH", "diphenyl phosphate", "O=P(O)(Oc1ccccc1)Oc1ccccc1"),
]:
    molecule(mid, label, P4, name, smi, "achiral" if smi else "unresolved", "not_applicable" if smi else "unresolved", "not_applicable" if smi else "unresolved")

for mid, label, centers in [("m_4a", "4a", 3), ("m_5a", "5a", 2), ("m_6a", "6a", 2)]:
    molecule(mid, label, P4, label, chirality="partially_specified", geometry="not_applicable", representation="not_encoded",
             atom_centers={"count_reported_in_text": centers, "absolute_configuration": "not fully established"},
             assignment="Source drawings show relative stereochemistry; 4a configurations were discussed as DFT predictions because X-ray assignment failed." if label == "4a" else "Parent 3a stereocenter is retained, while the new center is not assigned by name on the supplied page.",
             note="No full stereochemical SMILES asserted; er/dr are reaction outcomes, not identity completion.")

deriv_specs = [
    ("scheme3_cycloaddition_4a", "m_4a", [("m_ethyl_diazoacetate", "reagent", "2.0 equiv."), ("m_nitrosobenzene", "reagent", "2.0 equiv."), ("m_koac", "catalyst", "10 mol%"), ("m_dce", "solvent", None)], "DCE", "r.t.", 24, 90, "99:1", "95:5", [("C", "C", "formed", 0, 1), ("C", "O", "formed", 0, 1)], "three-component 1,3-dipolar cycloaddition"),
    ("scheme3_michael_5a", "m_5a", [("m_p_methylthiophenol", "reactant", "2.0 equiv."), ("m_et3n", "reagent", "2.0 equiv."), ("m_toluene", "solvent", None)], "toluene", "-20 °C", 24, 96, "99:1", "4.5:1", [("C", "S", "formed", 0, 1), ("C", "H", "hydrogen_gain", 0, 1)], "sulfa-Michael addition"),
    ("scheme3_michael_6a", "m_6a", [("m_indole", "reactant", "2.0 equiv."), ("m_diphenyl_phosphate", "reagent", "2.0 equiv."), ("m_toluene", "solvent", None)], "toluene", "r.t.", 24, 96, "99:1", "4.7:1", [("C", "C", "formed", 0, 1), ("C", "H", "hydrogen_gain", 0, 1)], "indole Michael addition"),
]
for suffix, product, agents, solvent_name, temp_raw, hours, yld, er, dr, bonds, trans in deriv_specs:
    r = reaction(suffix, P4, categories=["derivatization"], location=f"Scheme 3 product {product.split('_')[-1]}", transformation=trans,
                 unresolved=["Yield and er/dr measurement bases are not stated on the supplied page; no SI is associated."])
    participant(r, "p_3a", "m_3a", "reactant", P4)
    pids = ["p_3a"]
    for i, (mid, role, dose) in enumerate(agents, 1):
        pid = f"p_agent_{i}"
        if dose and dose.endswith("equiv."):
            val = float(dose.split()[0]); participant(r, pid, mid, role, P4, equiv=(dose, q(val, "equiv")), reference="m_3a")
        elif dose and dose.endswith("mol%"):
            val = float(dose.split()[0]); participant(r, pid, mid, role, P4, mol_percent=(dose, q(val, "mol%")), reference="m_3a")
        else:
            participant(r, pid, mid, role, P4)
        pids.append(pid)
    participant(r, "p_product", product, "product", P4, product_role="derivatization product")
    pids.append("p_product")
    r["solvent_conditions"]["composition"] = direct(solvent_name, [solvent_name], P4)
    temp_norm = "r.t." if temp_raw == "r.t." else q(-20, "°C")
    step(r, "derivatization_reaction", P4, pids, operation=trans, temperature=(temp_raw, temp_norm), time=(f"24 h", q(hours, "h")))
    outcome(r, "reported_yield", "yield", [product], P4, f"{yld}%", q(yld, "%"), product_role="derivatization product")
    er_outcome(r, "er", product, er, P4, basis=None)
    outcome(r, "dr", "selectivity", [product], P4, f"{dr} dr", dr, scope={"metric": "diastereomeric ratio", "order": "major:minor", "configuration_order": "unassigned"})
    set_classification(r, P4, bonds)

# Non-synthetic observations and page inventory context.
bundle["source_context"] = [
    {"context_id": f"{DOC}:article_metadata", "document_id": DOC, "evidence_ids": P1,
     "title": "Constructing chiral aza-quaternary carbon centers by enantioselective carbonylative Heck reaction of o-iodoanilines with allenes",
     "doi": "10.1039/d0qo01486a", "note": "True source title from page 1."},
    {"context_id": f"{DOC}:scheme1A_natural_products", "document_id": DOC, "evidence_ids": P1,
     "note": "Selected natural products/pharmaceuticals and bioactivity values are context, not synthetic experiments in this paper."},
    {"context_id": f"{DOC}:page2_aggregate_screening", "document_id": DOC, "evidence_ids": P2,
     "note": "Aggregate prose reports testing 29 ligands and further solvent/Pd:P/CO-pressure screening; missing SI rows are not fabricated as separate experiments."},
    {"context_id": f"{DOC}:figure1_dft", "document_id": DOC, "evidence_ids": P4,
     "note": "Figure 1 shows two DFT transition states for 4a formation. It is computational mechanistic context, not a synthetic experiment or observed structure assignment."},
    {"context_id": f"{DOC}:pages5_6_references", "document_id": DOC, "evidence_ids": ev(5) + ev(6),
     "note": "Pages 5-6 contain acknowledgements/references only; no current-paper reaction arrows, tables, or standalone characterization peak lists were found."},
]

st.write(RUN / "submitted.json", bundle)
st.write(RUN / "raw_model/stereo_diagnostics.json", diagnostics)
inventory = {
    "document_id": DOC,
    "reviewed_pages": [1, 2, 3, 4, 5, 6],
    "reaction_count": len(bundle["reactions"]),
    "molecule_count": len(bundle["molecules"]),
    "reaction_ids": [r["reaction_id"] for r in bundle["reactions"]],
    "notes": ["Table 1 entry 4 and Scheme 2 product 3a are one experiment with multiple evidence regions.",
              "Pages 5-6 contain only acknowledgements/references.",
              "Page 7 was used only from extracted text as necessary context and is not claimed as visually covered."],
}
st.write(RUN / "raw_model/page_inventory.json", inventory)
page_regions = {
    1: [
        {"region": "title, abstract, and article metadata", "finding": "Current-paper title/DOI and high-level yield/er claims; the abstract is context rather than a separate experiment.", "evidence_ids": ev(1)},
        {"region": "Scheme 1A", "finding": "Natural products/pharmaceuticals and bioactivity values; non-synthetic context for this paper.", "evidence_ids": ev(1)},
        {"region": "Scheme 1B", "finding": "Three prior-work reaction arrows (aerobic oxidation, radical bicyclization, [4+2] cycloaddition), recorded as literature examples.", "evidence_ids": ev(1)},
        {"region": "Scheme 1C", "finding": "Prior Pd-catalysed non-asymmetric cyclocarbonylation generic arrow, recorded as a literature example.", "evidence_ids": ev(1)},
        {"region": "Scheme 1D", "finding": "Generic scheme for this work with CO balloon, 80 °C, 6 h and generic maximum yield/er claims.", "evidence_ids": ev(1)},
        {"region": "ESI footnote", "finding": "Paper states ESI availability, but no SI file was supplied/associated in this workspace; no SI data were imported.", "evidence_ids": ev(1)},
    ],
    2: [
        {"region": "left-column model-reaction and ligand-screen discussion", "finding": "Model substrates 1a/2a and screening at 80 °C for 24 h; explicit partial results for (R)-BINAP and (R)-BTFM-GarPhos are recorded, while the aggregate 29-ligand inventory remains context.", "evidence_ids": ev(2)},
        {"region": "right-column optimization discussion", "finding": "Explicit no-conversion controls for 2-bromoaniline and 2-aminophenyl triflate, plus three ligand-free catalyst cases giving quantitative 1a consumption and a little racemic 3a.", "evidence_ids": ev(2)},
        {"region": "right-column Table 1 discussion", "finding": "Syringe-pump/additive comparisons refer to Table 1 entries on page 3; entry 4 is cross-linked rather than duplicated.", "evidence_ids": ev(2)},
    ],
    3: [
        {"region": "Table 1 reaction scheme and rows 1-6", "finding": "Six distinct optimization experiments with additive, feeding mode, GC yield and er; entry 4 also has isolated/recrystallized yields and enriched er.", "evidence_ids": ev(3)},
        {"region": "Table 1 footnotes a-d", "finding": "Common charges, GC/decane and chiral-HPLC methods, syringe-pump solvent split, atmospheric CO, isolated-yield and recrystallization notes, and entry-5 2 h addition plus 4 h hold were structured.", "evidence_ids": ev(3)},
        {"region": "Scheme 2 general arrow/caption", "finding": "Generic scope conditions and isolated-yield/chiral-HPLC bases; represented as a generic scheme.", "evidence_ids": ev(3)},
        {"region": "Scheme 2 products 3a-3s", "finding": "3a is the same experiment as Table 1 entry 4; 3b-3s are separate scope experiments with product-specific isolated yields and er, including post-recrystallization results for 3c and 3f.", "evidence_ids": ev(3)},
        {"region": "Scheme 2 products 3t and 3u with lower-right prose", "finding": "Desired products reported n.d.; 3t instead forms 2-methyl-4H-3,1-benzoxazin-4-one, while 3u gives inseparable side products. Missing amounts/detection bases are preserved.", "evidence_ids": ev(3)},
        {"region": "bottom-right scale-up lead-in", "finding": "Text introduces the distinct 1.0 mmol scale-up shown on page 4; evidence is merged into one scale-up record.", "evidence_ids": ev(3)},
    ],
    4: [
        {"region": "Scheme 3 left-hand scale-up arrow", "finding": "One 1.0 mmol 1a / 1.2 mmol 2a scale-up experiment at 80 °C for 6 h, with 86%/89:11 er and 60%/99:1 after recrystallization.", "evidence_ids": ev(4)},
        {"region": "Scheme 3 upper derivatization arrow", "finding": "Sulfa-Michael addition of para-methylthiophenol to 3a at -20 °C for 24 h affords 5a with 96%, 99:1 er, 4.5:1 dr.", "evidence_ids": ev(4)},
        {"region": "Scheme 3 right derivatization arrow", "finding": "Three-component 1,3-dipolar cycloaddition in DCE at room temperature for 24 h affords 4a with 90%, 99:1 er, 95:5 dr.", "evidence_ids": ev(4)},
        {"region": "Scheme 3 lower derivatization arrow", "finding": "Indole Michael addition in toluene at room temperature for 24 h affords 6a with 96%, 99:1 er, 4.7:1 dr.", "evidence_ids": ev(4)},
        {"region": "Figure 1 and DFT discussion", "finding": "Two calculated transition states and predicted configurations for 4a; retained as computational source context, not a synthetic experiment or experimental absolute-configuration proof.", "evidence_ids": ev(4)},
    ],
    5: [{"region": "full page", "finding": "Acknowledgements and references 1-12 only; no current-paper experimental arrow/table or standalone spectrum peak list.", "evidence_ids": ev(5)}],
    6: [{"region": "full page", "finding": "References 13-21 only; cited literature is bibliography rather than new experiments in the current paper.", "evidence_ids": ev(6)}],
}
coverage_pages = []
for page in range(1, 7):
    page_evidence = set(ev(page))
    ids = [rxn["reaction_id"] for rxn in bundle["reactions"] if page_evidence.intersection(rxn["evidence_ids"])]
    unresolved = []
    if page in (1, 2, 3, 4):
        unresolved.append("No explicitly supplied/associated SI; compound-specific procedures, characterization, and omitted table details remain unavailable.")
    if page == 2:
        unresolved.append("Aggregate screening statements do not expose every SI row or all quantitative outcomes.")
    if page == 3:
        unresolved.append("Most scope structures and absolute configurations remain unencoded/unassigned pending expert structure review.")
    if page == 4:
        unresolved.append("Scheme 3 yield and er/dr measurement bases are not stated on the supplied page; 4a non-parent configurations are computational predictions.")
    coverage_pages.append({"page": page, "state": "attempted", "regions": page_regions[page], "record_ids": ids, "unresolved": unresolved})
st.write(RUN / "coverage/pages_0001_0006.json", {"document_id": DOC, "semantic_accuracy_verified": False, "pages": coverage_pages})
print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"])}, ensure_ascii=False))
