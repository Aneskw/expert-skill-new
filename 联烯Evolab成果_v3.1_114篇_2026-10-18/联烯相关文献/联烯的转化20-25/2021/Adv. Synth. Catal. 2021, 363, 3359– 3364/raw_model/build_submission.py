#!/usr/bin/env python3
"""Source-specific v3.1 builder for Zheng et al., Adv. Synth. Catal. 2021, 363, 3359–3364."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DATASET = st.read(RUN / "dataset.json")
DOC = "doc_25b2677dbb7fd376"
EV = {p: f"{DOC}:page{p}:image" for p in range(1, 7)}
EVT = {p: f"{DOC}:page{p}:text" for p in range(1, 7)}


def f(raw=None, normalized=None, pages=(), method="original_extraction", derivation=None,
      status=None, verification=None, alternatives=()):
    evidence = []
    for p in pages:
        evidence.extend([EV[p], EVT[p]])
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation,
                   status=status, verification=verification, alternatives=alternatives)


def missing(pages=(), status="unreported"):
    return f(pages=pages, status=status)


def q(number, unit, raw=None):
    value = {"number": number, "unit": unit}
    if raw is not None:
        value["raw_text"] = raw
    return value


def mol(mid, label, name=None, pages=(2,), smiles=None, stereo="unresolved",
        geometry="unresolved", representation="unresolved", note=None, axial=None,
        assignment=None, mixture=None):
    evidence = [x for p in pages for x in (EV[p], EVT[p])]
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, pages)
    m["names"]["original"] = f(name or label, name or label, pages)
    m["names"]["english"] = f(name, name, pages) if name else missing(pages)
    m["names"]["abbreviations"] = missing(pages)
    m["names"]["aliases"] = missing(pages)
    m["cas_numbers"] = missing(pages)
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = f(
            {"source": f"structure drawing associated with {label}", "representation": "candidate SMILES"},
            smiles, pages, method="structure_derivation",
            derivation="Candidate connectivity transcribed from the source drawing/name; parsing does not establish identity.",
        )
        try:
            diag = st.stereo_inspect(smiles)
            m["stereochemistry"]["isomeric_smiles"] = f(
                "RDKit canonical isomeric form of candidate connectivity",
                diag["canonical_isomeric_smiles"], pages, method="structure_derivation",
                derivation="Computed from the candidate SMILES; source stereochemical assignments remain separate.",
            )
            m["stereochemistry"]["identity_key"] = f(
                None, diag["stereo_graph_key"], pages, method="structure_derivation",
                derivation="stereo_graph_key computed from the candidate graph; not a sample identity.",
            )
        except Exception:
            pass
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing(pages)
    for key, val, why in (
        ("chirality_status", stereo, "Source drawing/name was reviewed for specified center/axial stereochemistry."),
        ("geometry_status", geometry, "Source drawing/name was reviewed for specified E/Z or cis/trans geometry."),
        ("representation_status", representation, "Status reflects whether the candidate representation preserves source stereo information."),
    ):
        m["stereochemistry"][key] = f(
            f"source stereo status for {label}", val, pages, method="inference", derivation=why
        )
    m["stereochemistry"]["atom_centers"] = missing(pages)
    m["stereochemistry"]["double_bonds"] = missing(pages)
    m["stereochemistry"]["axial"] = (
        f(axial, axial, pages) if axial is not None else missing(pages)
    )
    m["stereochemistry"]["mixture"] = (
        f(mixture, mixture, pages) if mixture is not None else missing(pages)
    )
    m["stereochemistry"]["assignment"] = (
        f(assignment, assignment, pages) if assignment is not None else missing(pages)
    )
    m["identity_note"] = f(note, note, pages) if note else missing(pages)
    m["external_properties"] = []
    return m


def generic_mol(mid, label, fixed, variables, pages, reason):
    m = mol(mid, label, pages=pages, stereo="unresolved", geometry="unresolved",
            representation="not_encoded", note=reason)
    m["representation_type"] = "generic_template"
    m["smiles"] = missing(pages)
    m["generic_template"] = {
        "fixed_structure": f(fixed, fixed, pages),
        "variables": [{
            "label": v,
            "attachment_points": f("source drawing attachment point", "source drawing attachment point", pages),
            "definition": f(desc, desc, pages),
            "allowed_values": missing(pages),
            "constraints": missing(pages),
        } for v, desc in variables],
        "instance_mapping": [],
        "reason": f(reason, reason, pages),
    }
    return m


MOLECULES = []
def addm(m):
    MOLECULES.append(m)
    return m["molecule_id"]


# Generic drawings and the repeatedly used concrete core molecules.
addm(generic_mol("m_allenoic_generic", "1 (generic 2,3-allenoic acid)",
                 "aryl-substituted 2,3-allenoic acid drawing", [("R1", "aryl-ring substituent"),
                 ("R2", "substituent at the distal allene carbon"), ("R3", "alpha substituent")],
                 (2, 3), "R-group drawing; a wildcard SMILES would not be a concrete molecule."))
addm(generic_mol("m_oxime_generic", "2 (generic cycloketone oxime ester)",
                 "O-(4-(trifluoromethyl)benzoyl) cycloketone oxime ester drawing",
                 [("R4", "cycloketone substituent"), ("R5", "cycloketone substituent")],
                 (2, 3), "R-group drawing and variable ring substitution."))
addm(generic_mol("m_butenolide_generic", "3/4 (generic cyanoalkylsulfonylated butenolide)",
                 "cyanoalkylsulfonylated furan-2(5H)-one drawing",
                 [("R1", "aryl substituent"), ("R2", "C5 substituent"), ("R3", "alkene substituent"),
                  ("R4", "cyanoalkyl-chain substituent"), ("R5", "cyanoalkyl-chain substituent")],
                 (2, 3), "Product-family template with multiple source variables."))

addm(mol("m_1a", "1a", "2-methyl-4-phenyl-2,3-allenoic acid", (2, 3, 4),
         "CC(=C=Cc1ccccc1)C(=O)O", stereo="unspecified", geometry="unspecified",
         representation="partial",
         axial={"axis": "C=C=C", "assignment": "not specified", "note": "Both allene termini bear unlike substituents; the source gives no axial descriptor."},
         note="Candidate allene graph; RDKit flags cumulene double-bond stereo as unspecified and does not establish axial configuration."))
addm(mol("m_2a", "2a", "cyclobutanone O-(4-(trifluoromethyl)benzoyl) oxime", (2, 3, 4),
         "O=C(ON=C1CCC1)c1ccc(C(F)(F)F)cc1", stereo="unspecified", geometry="unspecified",
         representation="partial", note="Oxime geometry is not assigned in the supplied pages."))
addm(mol("m_3a", "3a", "cyanoalkylsulfonylated butenolide 3a", (2, 3, 4),
         "CC1=C(S(=O)(=O)CCC#N)C(c2ccccc2)OC1=O", stereo="unspecified",
         geometry="not_applicable", representation="partial",
         assignment={"tetrahedral_center": "C5 configuration not specified"},
         note="Candidate graph contains an unassigned tetrahedral center; no racemate or absolute configuration is reported."))

# Scope identities. Exact substitution names are not guessed where the figure is the only source.
for label in list("bcdefghijklmno"):
    addm(mol(f"m_1{label}", f"1{label}", f"2,3-allenoic acid substrate leading to 3{label}", (3,),
             stereo="unresolved", geometry="unresolved", representation="unresolved",
             note="Structure is visually distinguished in Scheme 2, but a complete candidate graph is withheld pending expert structure transcription."))
    addm(mol(f"m_3{label}", f"3{label}", f"cyanoalkylsulfonylated butenolide 3{label}", (3,),
             stereo="unspecified", geometry="unresolved", representation="unresolved",
             assignment={"source": "no absolute/relative stereochemical descriptor printed"},
             note="Product drawing and label retained; exact candidate SMILES is unresolved."))

for label in list("abcdefghijklm"):
    addm(mol(f"m_2_for_4{label}", f"2 (substrate leading to 4{label})",
             f"substituted cyclobutanone oxime ester leading to 4{label}", (3,),
             stereo="unresolved", geometry="unresolved", representation="unresolved",
             note="Scheme 3 identifies the substrate by its paired product drawing; an individual 2-series label is not printed."))
    has_dr = label in "abcdefgijkm"
    addm(mol(f"m_4{label}", f"4{label}", f"cyanoalkylsulfonylated butenolide 4{label}", (3,),
             stereo="mixture" if has_dr else "unspecified", geometry="unresolved", representation="not_encoded",
             mixture={"type": "diastereomeric mixture", "ratio": "experiment-specific"} if has_dr else None,
             note="Scheme 3 product drawing retained; reported dr, where present, is stored on the experiment."))
addm(mol("m_2n", "2n", "unsuccessful cycloketone oxime substrate 2n", (3,),
         stereo="unresolved", geometry="unresolved", representation="unresolved"))
addm(mol("m_2o", "2o", "unsuccessful Boc-protected azetidinone oxime substrate 2o", (3,),
         stereo="unresolved", geometry="unresolved", representation="unresolved"))

# Reagents, catalysts, solvents, workup agents, controls and derived products.
simple = [
    ("m_k2s2o5", "K2S2O5", "potassium metabisulfite", (2, 3, 4), "[K+].[K+].[O-]S(=O)OS(=O)[O-]"),
    ("m_na2s2o5", "Na2S2O5", "sodium metabisulfite", (2,), "[Na+].[Na+].[O-]S(=O)OS(=O)[O-]"),
    ("m_fecl2", "FeCl2", "iron(II) chloride", (2, 3, 4), "[Cl-].[Cl-].[Fe+2]"),
    ("m_dce", "DCE", "1,2-dichloroethane", (2, 3, 4), "ClCCCl"),
    ("m_mecn", "CH3CN", "acetonitrile", (2,), "CC#N"),
    ("m_thf", "THF", "tetrahydrofuran", (2,), "C1CCOC1"),
    ("m_dioxane", "1,4-dioxane", "1,4-dioxane", (2,), "C1COCCO1"),
    ("m_dcm", "dichloromethane", "dichloromethane", (4,), "ClCCl"),
    ("m_na2co3", "sat. aq. Na2CO3", "saturated aqueous sodium carbonate", (4,), None),
    ("m_na2so4", "Na2SO4", "sodium sulfate", (4,), "[Na+].[Na+].[O-]S(=O)(=O)[O-]"),
    ("m_nbuoh", "n-BuOH", "1-butanol", (3,), "CCCCO"),
    ("m_zncl2", "ZnCl2", "zinc chloride", (3,), "[Cl-].[Cl-].[Zn+2]"),
    ("m_h2so4", "H2SO4", "sulfuric acid", (3,), "OS(=O)(=O)O"),
    ("m_tbuoac", "t-BuOAc", "tert-butyl acetate", (3,), "CC(=O)OC(C)(C)C"),
    ("m_temppo", "TEMPO", "2,2,6,6-tetramethylpiperidin-1-oxyl", (3, 4), "CC1(C)CCCC(C)(C)N1[O]"),
    ("m_bht", "BHT", "2,6-di-tert-butyl-4-methylphenol", (3, 4), "Cc1cc(C(C)(C)C)c(O)c(C(C)(C)C)c1"),
    ("m_dpe", "1,1-diphenylethene", "1,1-diphenylethene", (3, 4), "C=C(c1ccccc1)c1ccccc1"),
]
for mid, label, name, pages, smi in simple:
    resonance_partial = mid in {"m_k2s2o5", "m_na2s2o5"}
    addm(mol(mid, label, name, pages, smi, stereo="achiral" if smi else "unresolved",
             geometry="unspecified" if resonance_partial else "not_applicable" if smi else "unresolved",
             representation="partial" if resonance_partial else "encoded" if smi else "unresolved",
             assignment={"basis": "candidate graph has no unassigned potential stereocenter in standard stereo diagnostics"} if smi else None))

for mid, label, name in [
    ("m_cuoac2", "Cu(OAc)2", "copper(II) acetate"), ("m_cuoac", "CuOAc", "copper(I) acetate"),
    ("m_niacac2", "Ni(acac)2", "nickel(II) acetylacetonate"),
    ("m_cuacn4bf4", "Cu(CH3CN)4BF4", "tetrakis(acetonitrile)copper(I) tetrafluoroborate"),
    ("m_feotf2", "Fe(OTf)2", "iron(II) triflate"), ("m_feotf3", "Fe(OTf)3", "iron(III) triflate"),
    ("m_fecl3", "FeCl3", "iron(III) chloride"), ("m_feoac2", "Fe(OAc)2", "iron(II) acetate"),
    ("m_dabso", "DABSO", "DABCO-bis(sulfur dioxide) adduct"),
]:
    addm(mol(mid, label, name, (2,), stereo="unresolved", geometry="not_applicable",
             representation="unresolved", note="Identity named in Table 1; no source SMILES is asserted."))

for mid, label, name, pages in [
    ("m_5a", "5a", "butanolysis product 5a", (3,)), ("m_5b", "5b", "N-tert-butyl amide 5b", (3,)),
    ("m_6", "6", "TEMPO-trapped cyanoalkyl radical product 6", (3, 4)),
    ("m_7", "7", "DPE trapping product 7", (3, 4)), ("m_8", "8", "DPE trapping product 8", (3, 4)),
    ("m_9", "9", "ethyl 2-methyl-4-phenyl-2,3-allenoate", (3, 4)),
    ("m_10", "10", "preformed butenolide control substrate 10", (3, 4)),
]:
    smi = "CCOC(=O)C(C)=C=Cc1ccccc1" if mid == "m_9" else None
    addm(mol(mid, label, name, pages, smi, stereo="unspecified" if mid == "m_9" else "unresolved",
             geometry="unspecified" if mid == "m_9" else "unresolved",
             representation="partial" if mid == "m_9" else "unresolved",
             axial={"axis": "C=C=C", "assignment": "not specified"} if mid == "m_9" else None))
for label in "ABCDE":
    addm(mol(f"m_mech_{label}", f"mechanistic intermediate {label}",
             f"proposed intermediate {label}", (4,), stereo="unresolved", geometry="unresolved",
             representation="unresolved", note="Proposed intermediate, not an isolated experiment product."))


def participant(pid, mid, role, pages, stage="reaction", amount=None, mass=None, volume=None,
                equiv=None, molpct=None, ref=None, product_role=None, inherited=False):
    evidence = [x for p in pages for x in (EV[p], EVT[p])]
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    method = "inherited_general_condition" if inherited else "original_extraction"
    if amount is not None:
        p["charges"]["amount"] = f(str(amount), amount, pages, method=method)
    if mass is not None:
        p["charges"]["mass"] = f(str(mass), mass, pages, method=method)
    if volume is not None:
        p["charges"]["volume"] = f(str(volume), volume, pages, method=method)
    if equiv is not None:
        p["charges"]["equivalents"] = f(str(equiv), equiv, pages, method=method)
    if molpct is not None:
        p["charges"]["mol_percent"] = f(str(molpct), molpct, pages, method=method)
    if ref:
        p["equivalence_reference"] = f(ref, ref, pages, method=method)
    if product_role:
        p["product_role"] = f(product_role, product_role, pages)
    return p


def manual_class(r, changes, pages, coverage="partial", limits=None):
    facts = []
    labels = {}
    for ch in changes:
        cid, a, b, before, after, kind = ch
        obj = {"change_id": cid, "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
               "before_order": before, "after_order": after, "change_type": kind}
        facts.append(f(None, obj, pages, method="inference",
                       derivation="Class-level net change inferred from the drawn reactant/product family; exact atom mapping is unresolved."))
        if kind in ("formed", "hydrogen_gain", "order_increased"):
            lab = st.bond_label(a, b, after)
            labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(cid)
    r["classification"]["bond_changes"] = facts
    r["classification"]["labels"] = f(None, list(labels.values()), pages, method="inference",
                                           derivation="Taxonomy labels computed from the retained class-level bond changes.")
    r["classification"]["coverage"] = f(None, coverage, pages, method="inference",
                                             derivation="Exact atom mapping and complete coproduct accounting are unavailable.")
    r["classification"]["limitations"] = limits or [
        "Class-level changes only; exact atom identities and full atom balance are unresolved.",
        "No atom-mapped reaction SMILES is asserted.",
    ]
    r["classification"]["model_votes"] = []


MAIN_CHANGES = [
    ("c_s", "C", "S", 0, 1, "formed"),
    ("c_o", "C", "O", 0, 1, "formed"),
    ("cn_triple", "C", "N", 2, 3, "order_increased"),
    ("ring_cc_break", "C", "C", 1, 0, "broken"),
    ("no_break", "N", "O", 1, 0, "broken"),
]


def base_reaction(rid, pages, record_type="specific_experiment", categories=("other",),
                  name=None, transformation=None, location=None):
    evidence = [x for p in pages for x in (EV[p], EVT[p])]
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = f(list(categories), list(categories), pages) if record_type == "specific_experiment" else missing(pages, "not_applicable")
    r["location"] = f(location or f"PDF page(s) {','.join(map(str, pages))}",
                        {"pages": list(pages), "region": location or rid}, pages)
    r["reaction_name"] = f(name, name, pages) if name else missing(pages)
    r["transformation_type"] = f(None, transformation, pages, method="inference",
                                     derivation="Transformation class inferred from source title/caption and structures.") if transformation else missing(pages)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": [
                       "Supporting information was advertised but not supplied or associated.",
                       "Machine extraction remains needs_review; schema validation is not chemistry verification.",
                   ]}
    r["completeness"] = "partial"
    return r


def set_conditions(r, substrate_mid, oxime_mid, product_mid, pages, *, catalyst_mid="m_fecl2",
                   sulfur_mid="m_k2s2o5", solvent_mid="m_dce", temp=60, hours=24,
                   catalyst_pct=10, scope=True, yield_value=None, yield_raw=None, isolated=True,
                   include_workup=False):
    method_inherited = scope
    r["participants"] = [
        participant("substrate", substrate_mid, "reactant", pages, amount=q(0.2, "mmol"), equiv=q(1, "equiv"), ref=substrate_mid, inherited=method_inherited),
        participant("oxime", oxime_mid, "reactant", pages, amount=q(0.3, "mmol"), equiv=q(1.5, "equiv"), ref=substrate_mid, inherited=method_inherited),
        participant("sulfur", sulfur_mid, "reagent", pages, amount=q(0.4, "mmol") if scope else None,
                    equiv=q(2, "equiv"), ref=substrate_mid, inherited=method_inherited),
        participant("catalyst", catalyst_mid, "catalyst", pages, amount=q(0.02, "mmol") if scope else None,
                    molpct=q(catalyst_pct, "mol%"), ref=substrate_mid, inherited=method_inherited),
        participant("solvent", solvent_mid, "solvent", pages, volume=q(2, "mL"), inherited=method_inherited),
    ]
    if product_mid:
        r["participants"].append(participant("product", product_mid, "product", pages, product_role="desired product"))
    r["reaction_scale"] = f("0.2 mmol of allenoic acid", {"reference_molecule_id": substrate_mid, "amount": q(0.2, "mmol")}, pages,
                                method="inherited_general_condition" if scope else "original_extraction")
    r["solvent_conditions"]["composition"] = f("DCE or table-listed replacement", [solvent_mid], pages,
                                                    method="inherited_general_condition" if scope else "original_extraction")
    r["solvent_conditions"]["total_volume"] = f("2 mL", q(2, "mL"), pages,
                                                    method="inherited_general_condition" if scope else "original_extraction")
    step = st.new_step("reaction", r["evidence_ids"], "reaction")
    step["participant_ids"] = [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"]
    step["operation"] = f("heated", "heated", pages, method="inherited_general_condition" if scope else "original_extraction")
    step["temperature"] = f(f"{temp} °C", q(temp, "°C"), pages, method="inherited_general_condition" if scope else "original_extraction")
    step["time"] = f(f"{hours} h", q(hours, "h"), pages, method="inherited_general_condition" if scope else "original_extraction")
    step["atmosphere"] = f("Ar atmosphere", "Ar", pages, method="inherited_general_condition" if scope else "original_extraction")
    r["steps"] = [step]
    r["solvent_conditions"]["mixture_ratio"] = missing(pages)
    r["solvent_conditions"]["concentration"] = missing(pages)
    r["solvent_conditions"]["concentration_basis"] = missing(pages)
    r["special_conditions"]["reactor"] = f("sealed tube (10 mL)", {"type": "sealed tube", "volume": q(10, "mL")}, pages,
                                                method="inherited_general_condition") if scope else missing(pages)
    r["procedure"]["general_procedure_reference"] = f("General Procedure for the Synthesis of Cyanoalkylsulfonylated Butenolides", "page 4 general procedure", pages,
                                                           method="inherited_general_condition") if scope else missing(pages)
    r["procedure"]["addition_order"] = f("A mixture of allenoic acid, oxime, K2S2O5, FeCl2 and DCE was heated.",
                                              [substrate_mid, oxime_mid, sulfur_mid, catalyst_mid, solvent_mid], pages,
                                              method="inherited_general_condition" if scope else "original_extraction")
    if include_workup:
        for pid, mid, role, stage, vol in [
            ("dcm_workup", "m_dcm", "workup_agent", "workup", q(15, "mL")),
            ("carbonate_wash", "m_na2co3", "workup_agent", "workup", q(10, "mL")),
            ("drying_agent", "m_na2so4", "workup_agent", "workup", None),
        ]:
            r["participants"].append(participant(pid, mid, role, pages, stage=stage, volume=vol, inherited=True))
        steps = [
            ("cool_dilute", "workup", ["dcm_workup"], "cooled to room temperature and diluted with dichloromethane (5 × 3 mL)"),
            ("wash", "workup", ["carbonate_wash"], "combined organic layers washed with saturated aqueous Na2CO3 (10 mL)"),
            ("dry_concentrate", "workup", ["drying_agent"], "dried over Na2SO4 and concentrated in vacuo"),
            ("column", "purification", [], "purified by column chromatography"),
        ]
        for sid, stage, pids, text in steps:
            s = st.new_step(sid, r["evidence_ids"], stage)
            s["participant_ids"] = pids
            s["operation"] = f(text, text, pages, method="inherited_general_condition")
            r["steps"].append(s)
        r["procedure"]["original_text"] = f(
            "A mixture ... was heated at 60 °C for 24 h in a sealed tube (10 mL) ... diluted with dichloromethane (5 × 3 mL) ... washed with sat. aq. Na2CO3 (10 mL), dried (Na2SO4), and concentrated in vacuo. The residue was purified by column chromatography.",
            "page 4 general procedure", pages, method="inherited_general_condition")
    if yield_value is not None:
        out = st.new_outcome("yield", "yield", [product_mid], r["evidence_ids"])
        out["product_role"] = f("desired product", "desired product", pages)
        out["value"] = f(yield_raw or f"{yield_value}%", q(yield_value, "%"), pages)
        out["yield_type"] = f("isolated yield", "isolated", pages) if isolated else missing(pages)
        out["measurement_basis"] = f("Isolated yield", "isolated", pages) if isolated else missing(pages)
        out["scope"] = f("single reaction outcome", {"start": substrate_mid, "end": product_mid}, pages)
        r["outcomes"].append(out)


REACTIONS = []
def addr(r):
    REACTIONS.append(r)
    return r


# Scheme 1: two cited reaction families and the current generic family.
for rid, loc, title, trans in [
    ("scheme1a_prior_pd", "Scheme 1a", "Palladium-catalyzed coupling cyclization of 2,3-dienoic acids", "cited coupling cyclization"),
    ("scheme1b_prior_cu", "Scheme 1b", "Copper-catalyzed free-radical cyclization of 2,3-dienoic acids", "cited free-radical cyclization"),
]:
    r = base_reaction(rid, (2,), "literature_example", (), title, trans, loc)
    r["participants"] = [participant("generic_allenoic", "m_allenoic_generic", "reactant", (2,)),
                         participant("generic_product", "m_butenolide_generic", "product", (2,), product_role="generic cited product")]
    r["notes"] = [f("Cited prior-work family; no concrete charge or yield is given in this paper.",
                     "Cited prior-work family; no concrete charge or yield is given in this paper.", (2,))]
    manual_class(r, [], (2,), "unresolved", ["Generic cited scheme; exact partners and atom mapping are not supplied here."])
    addr(r)

r = base_reaction("scheme1_current_generic", (2, 3), "generic_scheme", (),
                  "Iron-catalyzed three-component cyanoalkylsulfonylation", "radical cascade cyclization", "Scheme 1, this work; Schemes 2–3 generic headers")
r["participants"] = [participant("allenoic", "m_allenoic_generic", "reactant", (2, 3)),
                     participant("oxime", "m_oxime_generic", "reactant", (2, 3)),
                     participant("sulfur", "m_k2s2o5", "reagent", (2, 3)),
                     participant("product", "m_butenolide_generic", "product", (2, 3), product_role="generic product family")]
manual_class(r, MAIN_CHANGES, (2, 3))
addr(r)


# Table 1: every row independently; entry 18 also represents Scheme 2 product 3a.
rows = [
    (1, "m_cuoac2", "m_dce", 60, 18, "m_k2s2o5", 18),
    (2, "m_cuoac", "m_dce", 60, 18, "m_k2s2o5", "trace"),
    (3, "m_fecl2", "m_dce", 60, 18, "m_k2s2o5", 65),
    (4, "m_niacac2", "m_dce", 60, 18, "m_k2s2o5", 0),
    (5, "m_cuacn4bf4", "m_dce", 60, 18, "m_k2s2o5", 20),
    (6, "m_feotf2", "m_dce", 60, 18, "m_k2s2o5", 35),
    (7, "m_feotf3", "m_dce", 60, 18, "m_k2s2o5", 50),
    (8, "m_fecl3", "m_dce", 60, 18, "m_k2s2o5", 54),
    (9, "m_feoac2", "m_dce", 60, 18, "m_k2s2o5", 45),
    (10, "m_fecl2", "m_mecn", 60, 18, "m_k2s2o5", 25),
    (11, "m_fecl2", "m_thf", 60, 18, "m_k2s2o5", "trace"),
    (12, "m_fecl2", "m_dioxane", 60, 18, "m_k2s2o5", 8),
    (13, "m_fecl2", "m_dce", 60, 18, "m_na2s2o5", 42),
    (14, "m_fecl2", "m_dce", 60, 18, "m_dabso", 10),
    (15, "m_fecl2", "m_dce", 50, 18, "m_k2s2o5", 42),
    (16, "m_fecl2", "m_dce", 70, 18, "m_k2s2o5", 60),
    (17, "m_fecl2", "m_dce", 60, 12, "m_k2s2o5", 55),
    (18, "m_fecl2", "m_dce", 60, 24, "m_k2s2o5", 75),
]
for num, cat, solv, temp, hrs, sulfur, result in rows:
    pages = (2, 3, 4) if num == 18 else (2,)
    r = base_reaction(f"table1_entry{num}", pages, "specific_experiment", ("optimization",),
                      "Optimization of cyanoalkylsulfonylation", "radical cascade cyclization",
                      f"Table 1 entry {num}" + ("; duplicate Scheme 2 3a outcome merged" if num == 18 else ""))
    r["general_scheme_ids"] = [f"{DOC}:scheme1_current_generic"]
    set_conditions(r, "m_1a", "m_2a", "m_3a", pages, catalyst_mid=cat, sulfur_mid=sulfur,
                   solvent_mid=solv, temp=temp, hours=hrs, scope=False,
                   yield_value=result if isinstance(result, int) else None, isolated=True,
                   include_workup=(num == 18))
    if num == 18:
        next(p for p in r["participants"] if p["participant_id"] == "sulfur")["charges"]["amount"] = f(
            "0.4 mmol", q(0.4, "mmol"), (4,), method="inherited_general_condition")
        next(p for p in r["participants"] if p["participant_id"] == "catalyst")["charges"]["amount"] = f(
            "0.02 mmol", q(0.02, "mmol"), (4,), method="inherited_general_condition")
    if result == "trace":
        out = st.new_outcome("trace_3a", "trace", ["m_3a"], r["evidence_ids"])
        out["product_role"] = f("desired product", "desired product", pages)
        out["value"] = f("trace", "trace", pages)
        out["measurement_basis"] = f("Yield column footnote: isolated yield", "isolated", pages)
        r["outcomes"].append(out)
    if result == 0:
        r["review"]["unresolved"].append("A 0% desired-product yield is reported; conversion, recovery, and byproducts are unreported.")
        manual_class(r, [], pages, "not_applicable", ["No observed product transformation is classified from a 0% desired-product result."])
    elif result == "trace":
        manual_class(r, MAIN_CHANGES, pages)
    else:
        manual_class(r, MAIN_CHANGES, pages)
    if num == 18:
        r["notes"].append(f("Table 1 entry 18 and Scheme 2 compound 3a show the same charges, conditions, and 75% isolated yield; treated as one experiment with multiple evidence regions.",
                            "Table 1 entry 18 and Scheme 2 3a merged as duplicate reporting.", pages, method="inference",
                            derivation="Exact equality of substrate identities, charges, catalyst, solvent, temperature, time, atmosphere and yield across the two source regions."))
    addr(r)


# Scheme 2 additional substrate-scope rows (3a is merged into Table 1 entry 18).
scheme2 = {"b": 78, "c": 75, "d": 67, "e": 41, "f": 45, "g": 59, "h": 37,
           "i": 22, "j": 40, "k": 45, "l": 70, "m": 75, "n": 61, "o": 56}
for label, y in scheme2.items():
    r = base_reaction(f"scheme2_3{label}", (3, 4), "specific_experiment", ("substrate_scope",),
                      "Cyanoalkylsulfonylation of 2,3-allenoic acids", "radical cascade cyclization", f"Scheme 2, product 3{label}")
    r["general_scheme_ids"] = [f"{DOC}:scheme1_current_generic"]
    set_conditions(r, f"m_1{label}", "m_2a", f"m_3{label}", (3, 4), yield_value=y, isolated=True, include_workup=True)
    manual_class(r, MAIN_CHANGES, (3, 4))
    if label == "c":
        r["notes"].append(f("The crystal structure of 3c was confirmed by X-ray crystallographic analysis; CCDC 2073119.",
                            {"assignment": "X-ray", "ccdc": "2073119"}, (3, 6)))
    addr(r)


# Scheme 3: thirteen successful products and two unsuccessful substrates.
scheme3 = {
    "a": (36, "1.5:1"), "b": (38, "1.2:1"), "c": (60, "1.1:1"), "d": (31, "1.4:1"),
    "e": (43, "1.1:1"), "f": (37, "1.9:1"), "g": (40, "1.1:1"), "h": (31, None),
    "i": (25, "1:1"), "j": (39, "1:1"), "k": (36, "1.1:1"), "l": (51, None), "m": (19, "1.1:1"),
}
for label, (y, dr) in scheme3.items():
    r = base_reaction(f"scheme3_4{label}", (3, 4), "specific_experiment", ("substrate_scope",),
                      "Cyanoalkylsulfonylation using substituted cyclobutanone oxime esters", "radical cascade cyclization", f"Scheme 3, product 4{label}")
    r["general_scheme_ids"] = [f"{DOC}:scheme1_current_generic"]
    set_conditions(r, "m_1a", f"m_2_for_4{label}", f"m_4{label}", (3, 4), yield_value=y, isolated=True, include_workup=True)
    if dr:
        out = st.new_outcome("dr", "selectivity", [f"m_4{label}"], r["evidence_ids"])
        out["product_role"] = f("diastereomeric product mixture", "diastereomeric product mixture", (3,))
        out["value"] = f(f"dr = {dr}", dr, (3,))
        out["measurement_basis"] = missing((3,))
        out["scope"] = f("dr", {"metric": "dr", "order": "major/minor identities unassigned", "ratio": dr}, (3,))
        r["outcomes"].append(out)
        r["review"]["unresolved"].append("The dr is printed, but diastereomer identities/order assignment and measurement method are not reported.")
    manual_class(r, MAIN_CHANGES, (3, 4))
    addr(r)

for label in ("n", "o"):
    r = base_reaction(f"scheme3_2{label}_failed", (3, 4), "specific_experiment", ("substrate_scope",),
                      "Attempted cyanoalkylsulfonylation of an unsuccessful cycloketone oxime", "attempted radical cascade cyclization", f"Scheme 3 unsuccessful substrate 2{label}")
    r["general_scheme_ids"] = [f"{DOC}:scheme1_current_generic"]
    set_conditions(r, "m_1a", f"m_2{label}", None, (3, 4), yield_value=None, include_workup=False)
    out = st.new_outcome("desired_product_not_formed", "no_reaction", [], r["evidence_ids"])
    out["value"] = f("failed to form the desired product", True, (3,))
    out["scope"] = f("desired product formation", "desired product not identified", (3,))
    r["outcomes"] = [out]
    r["review"]["unresolved"].append("The source says the desired product failed to form; conversion, recovery, byproducts, and detection method are unreported.")
    manual_class(r, [], (3, 4), "not_applicable", ["No identified transformation product; intended chemistry is not classified."])
    addr(r)


# Scheme 4 derivatizations.
for rid, product, reagents, temp, y, trans, changes in [
    ("scheme4_5a", "m_5a", [("nbuoh", "m_nbuoh", "reagent"), ("zncl2", "m_zncl2", "reagent")], 120, 70,
     "lactone ring opening and butanolysis", [("new_ester_co", "C", "O", 0, 1, "formed"), ("ring_co_break", "C", "O", 1, 0, "broken")]),
    ("scheme4_5b", "m_5b", [("h2so4", "m_h2so4", "reagent"), ("tbuoac", "m_tbuoac", "reagent")], 70, 76,
     "Ritter reaction of nitrile to N-tert-butyl amide", [("amide_cn", "C", "N", 0, 1, "formed"),
      ("amide_co", "C", "O", 0, 2, "formed"), ("amide_nh", "N", "H", 0, 1, "hydrogen_gain"),
      ("nitrile_decrease", "C", "N", 3, 1, "order_decreased")]),
]:
    r = base_reaction(rid, (3,), "specific_experiment", ("derivatization",), trans, trans, f"Scheme 4, product {product[-2:]}")
    r["participants"] = [participant("start", "m_3a", "reactant", (3,)), participant("product", product, "product", (3,), product_role="desired product")]
    for pid, mid, role in reagents:
        r["participants"].append(participant(pid, mid, role, (3,)))
    step = st.new_step("reaction", r["evidence_ids"], "reaction")
    step["participant_ids"] = [p["participant_id"] for p in r["participants"]]
    step["operation"] = f(trans, trans, (3,))
    step["temperature"] = f(f"{temp} °C", q(temp, "°C"), (3,))
    step["time"] = missing((3,))
    r["steps"] = [step]
    out = st.new_outcome("yield", "yield", [product], r["evidence_ids"])
    out["product_role"] = f("desired product", "desired product", (3,))
    out["value"] = f(f"{y}%", q(y, "%"), (3,))
    out["yield_type"] = missing((3,))
    out["measurement_basis"] = missing((3,))
    out["scope"] = f("single Scheme 4 arrow", {"start": "m_3a", "end": product}, (3,))
    r["outcomes"] = [out]
    r["review"]["unresolved"].extend(["Reaction time, reagent quantities, workup, purification, and yield measurement basis are not reported on the supplied page."])
    manual_class(r, changes, (3,))
    addr(r)


# Scheme 5 controls: TEMPO and BHT are alternatives and are retained separately.
for trap, trap_mid in [("TEMPO", "m_temppo"), ("BHT", "m_bht")]:
    r = base_reaction(f"scheme5_{trap.lower()}_control", (3, 4), "specific_experiment", ("control", "mechanistic_study"),
                      f"Radical scavenger control with {trap}", "radical-trapping control", f"Scheme 5 experiment 1 ({trap} alternative)")
    set_conditions(r, "m_1a", "m_2a", "m_3a", (3, 4), scope=False, yield_value=0, isolated=False)
    r["participants"].append(participant("trap", trap_mid, "additive", (3, 4), equiv=q(2, "equiv"), ref="m_1a"))
    if trap == "TEMPO":
        r["participants"].append(participant("trap_product", "m_6", "product", (3, 4), product_role="radical-trapping product"))
        out = st.new_outcome("compound6_detected", "other", ["m_6"], r["evidence_ids"])
        out["product_role"] = f("radical-trapping product", "radical-trapping product", (4,))
        out["value"] = f("detected by GC-MS", "detected", (4,))
        out["measurement_basis"] = f("GC-MS", "GC-MS", (4,))
        r["outcomes"].append(out)
        manual_class(r, [("trap_co", "C", "O", 0, 1, "formed"), ("cn_triple", "C", "N", 2, 3, "order_increased")], (3, 4))
    else:
        manual_class(r, [], (3, 4), "not_applicable", ["Only suppression of desired product is assigned to this BHT alternative; no BHT adduct is identified."])
    r["review"]["unresolved"].append("The combined Scheme 5 alternative and prose do not fully delimit whether GC-MS detection of 6 applies beyond the TEMPO condition; 6 is assigned only to TEMPO because its drawn structure is a TEMPO adduct.")
    addr(r)

r = base_reaction("scheme5_dpe_control", (3, 4), "specific_experiment", ("control", "mechanistic_study"),
                  "1,1-Diphenylethene radical-trapping control", "radical-trapping control", "Scheme 5 experiment 2")
set_conditions(r, "m_1a", "m_2a", "m_3a", (3, 4), scope=False, yield_value=0, isolated=False)
r["participants"].append(participant("dpe", "m_dpe", "additive", (3, 4), equiv=q(3, "equiv"), ref="m_1a"))
for label, y in [("7", 34), ("8", 13)]:
    r["participants"].append(participant(f"product_{label}", f"m_{label}", "product", (3, 4), product_role="radical-trapping product"))
    o = st.new_outcome(f"yield_{label}", "yield", [f"m_{label}"], r["evidence_ids"])
    o["product_role"] = f("radical-trapping product", "radical-trapping product", (3, 4))
    o["value"] = f(f"{y}%", q(y, "%"), (3, 4))
    o["yield_type"] = f("isolated", "isolated", (3, 4))
    o["measurement_basis"] = f("could be isolated", "isolated", (3, 4))
    o["scope"] = f("DPE control", {"product": label}, (3, 4))
    r["outcomes"].append(o)
manual_class(r, [("dpe_cc", "C", "C", 0, 1, "formed"), ("cn_triple", "C", "N", 2, 3, "order_increased")], (3, 4))
addr(r)

for rid, start_mid, loc, statement in [
    ("scheme5_allenoate9_control", "m_9", "Scheme 5 experiment 3", "2,3-allenoate 9 gave 3a in 0% yield"),
    ("scheme5_butenolide10_control", "m_10", "Scheme 5 experiment 4", "10 and 2a did not yield 3a under standard conditions"),
]:
    r = base_reaction(rid, (3, 4), "specific_experiment", ("control", "mechanistic_study"), statement, "negative mechanistic control", loc)
    set_conditions(r, start_mid, "m_2a", "m_3a", (3, 4), scope=False, yield_value=0, isolated=False)
    r["review"]["unresolved"].append("The desired product is reported as 0%; conversion, recovery, byproducts, and analytical basis are unreported.")
    manual_class(r, [], (3, 4), "not_applicable", ["No observed transformation is classified for the negative control."])
    addr(r)


# Scheme 6 is a proposed mechanism, not a set of executed experiments.
r = base_reaction("scheme6_proposed_mechanism", (4,), "proposed_mechanism", (),
                  "Proposed free-radical cascade pathway", "proposed SET/radical addition/oxidation/cyclization", "Scheme 6")
r["participants"] = [
    participant("oxime", "m_2a", "reactant", (4,)), participant("allenoic", "m_1a", "reactant", (4,)),
    participant("sulfur", "m_k2s2o5", "reagent", (4,)), participant("catalyst", "m_fecl2", "catalyst", (4,)),
    participant("product", "m_3a", "product", (4,), product_role="proposed final product"),
] + [participant(f"intermediate_{x}", f"m_mech_{x}", "unknown", (4,)) for x in "ABCDE"]
r["notes"] = [f(
    "Author-proposed sequence: Fe(II) SET and N–O cleavage to iminyl radical A; C–C cleavage to cyanoalkyl radical B; SO2 capture to sulfonyl radical C; addition to 1a to radical D; Fe(III) oxidation to allylic cation E; intramolecular cyclization to 3a.",
    ["A", "B", "C", "D", "E", "3a"], (4,))]
manual_class(r, MAIN_CHANGES, (4,), "partial", ["Proposed mechanism only; intermediates A–E are not treated as executed experiments.",
                                                "Exact atom mapping and proton transfers are unresolved."])
addr(r)


BUNDLE = st.empty(DATASET["documents"])
BUNDLE["molecules"] = MOLECULES
BUNDLE["reactions"] = REACTIONS
BUNDLE["source_context"] = [
    {"context_id": f"{DOC}:figure1", "document_id": DOC, "kind": "biological_context",
     "evidence_ids": [EV[1], EVT[1]], "description": "Figure 1 shows biologically active furan-2(5H)-ones; it is not a synthesis experiment."},
    {"context_id": f"{DOC}:xray_3c", "document_id": DOC, "kind": "structure_assignment",
     "evidence_ids": [EV[3], EVT[3], EV[6], EVT[6]], "description": "Product 3c is assigned by X-ray crystallography; CCDC 2073119."},
    {"context_id": f"{DOC}:si_status", "document_id": DOC, "kind": "source_limitation",
     "evidence_ids": [EV[1], EVT[1]], "description": "The article advertises supporting information online, but no SI file was supplied or explicitly associated in this workspace."},
    {"context_id": f"{DOC}:spectra_inventory", "document_id": DOC, "kind": "negative_inventory",
     "evidence_ids": [EV[p] for p in range(1, 7)], "description": "No standalone mass-spectrum peak table, molecular-weight table, or fragment-relative-abundance block occurs on pages 1–6. GC-MS detection of compound 6 is retained as a control outcome."},
]

errors = st.validate(BUNDLE, RUN)
if errors:
    raise SystemExit(json.dumps(errors, ensure_ascii=False, indent=2))
st.write(RUN / "submitted.json", BUNDLE)

raw = {
    "document_id": DOC,
    "title": DATASET["documents"][0]["title"],
    "page_inventory": {
        "1": ["title/metadata/abstract", "Figure 1 biological context", "introduction; no current-paper reaction arrow"],
        "2": ["Scheme 1a cited Pd family", "Scheme 1b cited Cu family", "Scheme 1 current generic family", "Table 1 entries 1–18"],
        "3": ["Scheme 2 products 3a–3o", "Scheme 3 products 4a–4m and failed 2n/2o", "Scheme 4 arrows to 5a and 5b", "control-experiment prose"],
        "4": ["Scheme 5 controls: TEMPO/BHT alternatives, DPE, allenoate 9, butenolide 10", "Scheme 6 proposed mechanism", "general procedure"],
        "5": ["references only"], "6": ["references and CCDC 2073119 notice only"],
    },
    "deduplication": {"table1_entry18_and_scheme2_3a": "merged into one reaction because identity, charges, conditions, atmosphere and 75% isolated yield match"},
    "counts": {"molecules": len(MOLECULES), "reactions": len(REACTIONS)},
    "raw_status": "model-generated source interpretation; needs_review",
}
st.write(RUN / "raw_model/page_extraction.json", raw)

checks = {}
for mid in ("m_1a", "m_2a", "m_3a", "m_9", "m_temppo", "m_dpe"):
    m = next(x for x in MOLECULES if x["molecule_id"] == mid)
    smi = st.value(m["smiles"])
    checks[mid] = st.stereo_inspect(smi) if smi else {"source_verified": False, "note": "No candidate SMILES."}
st.write(RUN / "raw_model/stereo_checks.json", checks)
print(json.dumps({"molecules": len(MOLECULES), "reactions": len(REACTIONS)}, ensure_ascii=False))
