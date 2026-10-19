#!/usr/bin/env python3
"""Build the source-derived v3.1 submission for JACS 2021, 143, 4935-4941.

This file is a reproducible transcription helper, not an alternate schema or a
source-review certificate.  It writes only run artifacts and keeps every
reaction in needs_review/partial state.
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
DOC = "doc_2c8b9ffd03c82435"
EV = {n: f"{DOC}:page{n}:image" for n in range(1, 8)}
EVT = {n: f"{DOC}:page{n}:text" for n in range(1, 8)}
TITLE = "CuH-Catalyzed Regio- and Enantioselective Hydrocarboxylation of Allenes: Toward Carboxylic Acids with Acyclic Quaternary Centers"


def f(raw=None, normalized=None, evidence=(), method="original_extraction", derivation=None,
      status=None, verification=None, alternatives=()):
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation,
                   status=status, verification=verification, alternatives=alternatives)


def missing(evidence=()):
    return f(evidence=evidence)


def inferred(normalized, evidence, derivation):
    return f(None, normalized, evidence, method="inference", derivation=derivation)


def qty(number, unit, qualifier=None, raw_text=None):
    value = {"number": number, "unit": unit}
    if qualifier:
        value["qualifier"] = qualifier
    if raw_text:
        value["raw_text"] = raw_text
    return value


def range_qty(minimum, maximum, unit, raw_text):
    return {"minimum": minimum, "maximum": maximum, "unit": unit, "raw_text": raw_text}


bundle = json.loads((RUN / "submitted.json").read_text())
# The prior page-1-to-6 submission predates registration of the rendered page-7
# image.  `standard init` now discovers this source asset, so carry the same
# evidence object into the cumulative bundle before attaching page-7 facts.
document = next(doc for doc in bundle["documents"] if doc["document_id"] == DOC)
if not any(item["evidence_id"] == EV[7] for item in document["evidence"]):
    document["evidence"].append({
        "evidence_id": EV[7],
        "document_id": DOC,
        "page": 7,
        "kind": "page_image",
        "text": None,
        "asset": f"assets/{DOC}/page7.jpg",
        "bbox": None,
        "bbox_units": None,
        "provenance": "rendered_source_pdf",
    })
bundle["molecules"] = []
bundle["reactions"] = []
bundle["source_context"] = []
molecules = {}


def add_molecule(key, label, page, smiles=None, name=None, stereo="achiral",
                 assignment=None, representation_status=None, identity_note=None):
    if key in molecules:
        return key
    ev = [EV[page]]
    m = st.new_molecule(key, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    m["names"]["original"] = f(name or label, name or label, ev)
    m["names"]["english"] = f(name, name, ev) if name else missing(ev)
    m["names"]["abbreviations"] = missing(ev)
    m["names"]["aliases"] = missing(ev)
    m["cas_numbers"] = missing(ev)
    has_unspecified = False
    if smiles:
        diag = st.stereo_inspect(smiles)
        has_unspecified = any(item["specified"] == "Unspecified" for item in diag["potential_stereo"])
        canonical = diag["canonical_isomeric_smiles"]
        m["representation_type"] = "smiles"
        m["smiles"] = f({"source": "drawn structure", "label": label}, canonical, ev,
                         method="structure_derivation",
                         derivation="Candidate connectivity was transcribed from the page drawing and canonicalized with RDKit; parsing does not certify source identity.")
        m["stereochemistry"]["identity_key"] = f(None, diag["stereo_graph_key"], ev,
                                                     method="structure_derivation",
                                                     derivation="RDKit stereo graph key from the candidate registry SMILES.")
        if representation_status is None:
            representation_status = "partial" if has_unspecified else ("encoded" if stereo == "achiral" else "partial")
        if representation_status == "encoded":
            m["stereochemistry"]["isomeric_smiles"] = f(None, canonical, ev,
                                                           method="structure_derivation",
                                                           derivation="Canonical isomeric SMILES for a candidate graph without unresolved source stereochemistry.")
        else:
            m["stereochemistry"]["isomeric_smiles"] = missing(ev)
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing(ev)
        m["stereochemistry"]["identity_key"] = missing(ev)
        m["stereochemistry"]["isomeric_smiles"] = missing(ev)
        representation_status = representation_status or "unresolved"
    if stereo == "achiral":
        m["stereochemistry"]["chirality_status"] = inferred(
            "achiral", ev, "The transcribed candidate graph has no tetrahedral or axial stereogenic element; this is a structure-derived statement, not an inference from missing wedges.")
        m["stereochemistry"]["atom_centers"] = inferred([], ev, "No tetrahedral stereocenter in the candidate graph.")
        m["stereochemistry"]["axial"] = inferred([], ev, "Terminal allene or nonallene candidate has no configurationally stable allene axis.")
    elif stereo == "specified":
        m["stereochemistry"]["chirality_status"] = f("solid wedge in source drawing", "specified", ev)
        m["stereochemistry"]["atom_centers"] = f(
            "solid wedge at product stereocenter; no printed R/S descriptor",
            [{"center": "drawn product stereocenter", "configuration": "specified by wedge; absolute descriptor unassigned"}], ev)
        m["stereochemistry"]["axial"] = inferred([], ev, "The product contains no allene axis after hydrocarboxylation.")
    elif stereo == "named_S":
        m["stereochemistry"]["chirality_status"] = f("(S)-Indobufen", "specified", ev)
        m["stereochemistry"]["atom_centers"] = f("(S)-Indobufen", [{"center": "alpha carbon", "configuration": "S"}], ev)
        m["stereochemistry"]["assignment"] = f("(S)-Indobufen (10)", {"alpha_center": "S"}, ev)
        m["stereochemistry"]["axial"] = inferred([], ev, "No axial stereogenic element in the product connectivity.")
    else:
        m["stereochemistry"]["chirality_status"] = inferred("unresolved", ev, "Source label or drawing was not safely expanded to a complete stereo graph.")
    if assignment and stereo != "named_S":
        m["stereochemistry"]["assignment"] = f(assignment, assignment, ev)
    elif stereo == "achiral" and has_unspecified:
        m["stereochemistry"]["assignment"] = inferred(
            {"assessment": "no configurational assignment asserted", "reason": "candidate contains toolkit-flagged unspecified stereo"},
            ev,
            "The candidate is kept partial because standard stereo diagnostics flag an unspecified potential element; source connectivity and molecular symmetry/terminal substitution were considered, but no @/@@ or E/Z encoding is asserted.")
    elif stereo != "named_S":
        m["stereochemistry"]["assignment"] = missing(ev)
    m["stereochemistry"]["double_bonds"] = inferred([], ev, "No source-assigned E/Z alkene is present; terminal vinyl/allene geometry is not an E/Z assignment.")
    m["stereochemistry"]["mixture"] = missing(ev)
    m["stereochemistry"]["geometry_status"] = inferred("not_applicable", ev, "The depicted terminal vinyl/allene connectivity has no source-applicable E/Z assignment.")
    m["stereochemistry"]["representation_status"] = inferred(
        representation_status, ev,
        "Status compares the candidate SMILES with the source stereo drawing; wedge-defined products without safe @/@@ assignment remain partial or not_encoded.")
    m["identity_note"] = f(identity_note, identity_note, ev) if identity_note else f(
        "Structure candidate remains a machine transcription pending expert review.",
        "Structure candidate remains a machine transcription pending expert review.", ev)
    m["external_properties"] = []
    molecules[key] = m
    bundle["molecules"].append(m)
    return key


def add_unresolved(key, label, page, name=None, stereo="unresolved", assignment=None):
    return add_molecule(key, label, page, None, name, stereo=stereo,
                        assignment=assignment, representation_status="not_encoded" if stereo == "specified" else "unresolved",
                        identity_note="The source label/drawing was retained without forcing a complete SMILES.")


AD = "OC12CC3CC(CC(C3)C1)C2"

# Common reagents and clear identities.
add_molecule("m_1a", "1a", 2, "C=C=C(C)c1ccccc1", "1-phenyl-1-methylallene")
add_molecule("m_2a", "2a", 2, "O=C(Oc1ccccc1)Oc1ccccc1", "diphenyl carbonate")
add_molecule("m_2b", "2b", 2, "CC(C)(C)OC(=O)OC(=O)OC(C)(C)C", "di-tert-butyl dicarbonate")
add_molecule("m_2c", "2c", 2, "COC(=O)Cl", "methyl chloroformate")
add_molecule("m_2d", "2d", 2, "O=C(F)" + AD, "1-adamantyl fluoroformate")
add_molecule("m_cu_oac2", "Cu(OAc)2", 2, "CC(=O)[O-].CC(=O)[O-].[Cu+2]", "copper(II) acetate")
add_molecule("m_dmms", "(MeO)2MeSiH", 2, "CO[SiH](C)OC", "dimethoxy(methyl)silane")
add_molecule("m_thf", "THF", 2, "C1CCOC1", "tetrahydrofuran")
add_molecule("m_dme", "DME", 3, "COCCOC", "1,2-dimethoxyethane")
add_molecule("m_meoh", "MeOH", 3, "CO", "methanol")
add_molecule("m_tfa", "TFA", 3, "O=C(O)C(F)(F)F", "trifluoroacetic acid")
add_molecule("m_nh4f", "NH4F", 3, "[NH4+].[F-]", "ammonium fluoride")
add_unresolved("m_L1", "L1", 2, "(R,R)-Ph-BPE", stereo="specified", assignment="(R,R)-Ph-BPE")
add_unresolved("m_L2", "L2", 2, "(R)-DTBM-SEGPHOS", stereo="specified", assignment="(R)-DTBM-SEGPHOS")
add_unresolved("m_L3", "L3", 2, "SL-J011-1", stereo="specified", assignment="chiral ferrocene ligand SL-J011-1 as drawn")
add_molecule("m_liome", "LiOMe", 2, "C[O-].[Li+]", "lithium methoxide")
add_molecule("m_csobz", "CsOBz", 2, "O=C([O-])c1ccccc1.[Cs+]", "cesium benzoate")

# Explicit unresolved ligand identity for the Table 1 entry-5 source conflict.
conf = add_unresolved("m_table1_entry5_ligand_conflict", "L1/L2 conflict", 2)
molecules[conf]["names"]["original"] = f(
    ["Table 1 entry 5: L1", "page-2 prose: L2"], None, [EV[2]], status="conflict",
    alternatives=[f("Table 1 entry 5: L1", "L1", [EV[2]]), f("prose: (R)-DTBM-SEGPHOS (L2), Table 1 entry 5", "L2", [EV[2]])])

# Scope substrates and products.  Dense drawings remain unresolved rather than guessed.
scope = {
    "a": ("C=C=C(C)c1ccccc1", "C=CC(C)(c1ccccc1)C(=O)O", "acid"),
    "b": ("C=C=C(C)c1ccc(C(F)(F)F)cc1", "C=CC(C)(c1ccc(C(F)(F)F)cc1)C(=O)O", "acid"),
    "c": ("C=C=C(C)c1ccc(OC)cc1", "C=CC(C)(c1ccc(OC)cc1)C(=O)O", "acid"),
    "d": ("C=C=C(C)c1cccc(Br)c1", "C=CC(C)(c1cccc(Br)c1)C(=O)" + AD, "ester"),
    "e": ("C=C=C(C)c1ccccc1F", "C=CC(C)(c1ccccc1F)C(=O)" + AD, "ester"),
    "f": ("C=C=C(C)c1ccc2OCOc2c1", "C=CC(C)(c1ccc2OCOc2c1)C(=O)" + AD, "ester"),
    "g": (None, None, "ester"),
    "h": (None, None, "ester"),
    "i": (None, None, "acid"),
    "j": ("C=C=C(CC)c1ccccc1", "C=CC(CC)(c1ccccc1)C(=O)O", "acid"),
    "k": (None, None, "acid"),
    "l": (None, None, "acid"),
    "m": ("C=C=C(CCO[Si](C(C)C)(C(C)C)C(C)C)c1ccccc1", "C=CC(CCO[Si](C(C)C)(C(C)C)C(C)C)(c1ccccc1)C(=O)" + AD, "ester"),
    "n": ("C=C=C(C)C1CCCCC1", "C=CC(C)(C1CCCCC1)C(=O)" + AD, "ester"),
    "o": ("C=C=Cc1ccccc1", "C=CC(c1ccccc1)C(=O)" + AD, "ester"),
    "p": ("C=C=Cc1ccc(SC)cc1", "C=CC(c1ccc(SC)cc1)C(=O)" + AD, "ester"),
    "q": ("C=C=CC1CCCCC1", "C=CC(C1CCCCC1)C(=O)" + AD, "ester"),
}
for letter, (smi_in, smi_out, _) in scope.items():
    if letter == "a":
        # Reuse 1a but add the page-3 occurrence.
        molecules["m_1a"]["occurrences"].append({"document_id": DOC, "compound_label": f("1a", "1a", [EV[3]]), "evidence_ids": [EV[3]]})
    else:
        (add_molecule if smi_in else add_unresolved)(f"m_1{letter}", f"1{letter}", 3, smi_in) if smi_in else add_unresolved(f"m_1{letter}", f"1{letter}", 3)
    if smi_out:
        add_molecule(f"m_3{letter}", f"3{letter}", 3, smi_out, stereo="specified", representation_status="partial")
    else:
        add_unresolved(f"m_3{letter}", f"3{letter}", 3, stereo="specified", assignment="solid wedge in the source product drawing; absolute descriptor unassigned")

# Page-7 footnote 101 identifies an isolated ester intermediate in the 3l
# sequence, but neither names nor redraws it.  Keep it as a separate unresolved
# molecule instead of guessing a structure or stereochemical encoding.
m3l_ester = add_unresolved(
    "m_3l_ester_intermediate",
    "the corresponding ester",
    7,
    name="corresponding ester intermediate of 3l",
)
molecules[m3l_ester]["occurrences"][0]["compound_label"] = f(
    "the corresponding ester", "corresponding ester intermediate of 3l", [EV[7]]
)
molecules[m3l_ester]["names"]["original"] = f(
    "the corresponding ester", "corresponding ester intermediate of 3l", [EV[7]]
)
molecules[m3l_ester]["names"]["english"] = inferred(
    "corresponding ester intermediate of 3l",
    [EV[7]],
    "Descriptive English registry name derived from footnote 101; the source does not provide a systematic name or structure.",
)
for stereo_field in ("atom_centers", "double_bonds", "axial", "mixture", "assignment"):
    molecules[m3l_ester]["stereochemistry"][stereo_field] = f(
        "the corresponding ester; structure not depicted",
        "unresolved",
        [EV[7]],
    )
molecules[m3l_ester]["stereochemistry"]["chirality_status"] = f(
    "the corresponding ester; structure not depicted", "unresolved", [EV[7]]
)
molecules[m3l_ester]["stereochemistry"]["geometry_status"] = f(
    "the corresponding ester; structure not depicted", "unresolved", [EV[7]]
)
molecules[m3l_ester]["stereochemistry"]["representation_status"] = inferred(
    "unresolved",
    [EV[7]],
    "No structure drawing, systematic name, or stereochemical descriptor is provided for the isolated ester intermediate.",
)
molecules[m3l_ester]["stereochemistry"]["isomeric_smiles"] = missing([EV[7]])
molecules[m3l_ester]["stereochemistry"]["identity_key"] = missing([EV[7]])

# Optimization products depend on electrophile; 2d product is distinct from isolated acid 3a.
add_molecule("m_opt_prod_2a", "3 (phenyl ester from 2a)", 2, "C=CC(C)(c1ccccc1)C(=O)Oc1ccccc1", stereo="specified", representation_status="partial")
add_molecule("m_opt_prod_2b", "desired tert-butyl ester from 2b", 2, "C=CC(C)(c1ccccc1)C(=O)OC(C)(C)C", stereo="specified", representation_status="partial")
add_molecule("m_opt_prod_2c", "desired methyl ester from 2c", 2, "C=CC(C)(c1ccccc1)C(=O)OC", stereo="specified", representation_status="partial")
add_molecule("m_opt_prod_2d", "3 (1-adamantyl ester from 2d)", 2, "C=CC(C)(c1ccccc1)C(=O)" + AD, stereo="specified", representation_status="partial")

# Scheme 1 identities.
add_molecule("m_6", "6", 3, "C=CC(C)(N)c1ccccc1", stereo="specified", representation_status="partial")
add_unresolved("m_reagent7", "7", 3, name="electrophilic nitrogen reagent 7")
add_molecule("m_8", "8", 3, "NCCC(C)(c1cccc(Br)c1)C(=O)" + AD, stereo="specified", representation_status="partial")
add_molecule("m_dp_pazide", "(PhO)2P(O)N3", 3, "O=P(N=[N+]=[N-])(Oc1ccccc1)Oc1ccccc1", "diphenyl phosphoryl azide")
add_molecule("m_et3n", "Et3N", 3, "CCN(CC)CC", "triethylamine")
add_molecule("m_toluene", "toluene", 3, "Cc1ccccc1", "toluene")
add_molecule("m_hcl", "5 M HCl", 3, "Cl", "hydrochloric acid")
add_molecule("m_nh2oh_hcl", "NH2OH·HCl", 3, "NO.Cl", "hydroxylamine hydrochloride")
add_molecule("m_pto2", "PtO2", 3, "O=[Pt]=O", "platinum dioxide")
add_molecule("m_h2", "H2", 3, "[H][H]", "hydrogen")
add_molecule("m_dcm", "DCM", 3, "ClCCl", "dichloromethane")

indole = "N2Cc3ccccc3C2=O"
add_molecule("m_1r", "1r", 3, f"C=C=Cc1ccc({indole})cc1", stereo="achiral")
add_molecule("m_3r", "3r", 3, f"C=CC(c1ccc({indole})cc1)C(=O){AD}", stereo="specified", representation_status="partial")
add_molecule("m_9", "9", 3, f"CCC(c1ccc({indole})cc1)C(=O){AD}", stereo="specified", representation_status="partial")
add_molecule("m_10", "(S)-Indobufen (10)", 3, f"CCC(c1ccc({indole})cc1)C(=O)O", stereo="named_S", representation_status="not_encoded")


def set_charge(p, kind, raw, normalized, ev, reference=None, method="original_extraction"):
    p["charges"][kind] = f(raw, normalized, ev, method=method,
                             derivation="Inherited from the table-wide condition and applied to this row." if method == "inherited_general_condition" else None)
    if kind in ("equivalents", "mol_percent"):
        p["equivalence_reference"] = f(reference, reference, ev, method=method,
                                          derivation="The table identifies the limiting/reference reagent." if method == "inherited_general_condition" else None)


def participant(r, pid, mid, role, ev, stage="reaction", charge=None, reference=None, inherited=False):
    p = st.new_participant(pid, mid, role, ev, stage=stage)
    if charge:
        kind, raw, normalized = charge
        set_charge(p, kind, raw, normalized, ev, reference,
                   "inherited_general_condition" if inherited else "original_extraction")
    r["participants"].append(p)
    return p


def add_step(r, sid, stage, pids, ev, operation=None, temperature=None, time=None,
             atmosphere=None, pressure=None, inherited=False):
    s = st.new_step(sid, ev, stage=stage)
    s["participant_ids"] = list(pids)
    method = "inherited_general_condition" if inherited else "original_extraction"
    deriv = "Applied from the table-wide condition to this entry." if inherited else None
    if operation is not None:
        s["operation"] = f(operation, operation, ev, method=method, derivation=deriv)
    if temperature is not None:
        raw, norm = temperature
        s["temperature"] = f(raw, norm, ev, method=method, derivation=deriv)
    if time is not None:
        raw, norm = time
        s["time"] = f(raw, norm, ev, method=method, derivation=deriv)
    if atmosphere is not None:
        s["atmosphere"] = f(atmosphere, atmosphere, ev, method=method, derivation=deriv)
    if pressure is not None:
        raw, norm = pressure
        s["pressure"] = f(raw, norm, ev, method=method, derivation=deriv)
    r["steps"].append(s)
    return s


def manual_classification(r, changes, ev, limitation):
    labels = {}
    facts = []
    for idx, change in enumerate(changes, 1):
        c = dict(change)
        c["change_id"] = f"b{idx}"
        facts.append(f(None, c, ev, method="structure_derivation",
                       derivation="Net bond change inferred from the source reactant/product drawings; no complete atom mapping is claimed."))
        if c["change_type"] in ("formed", "hydrogen_gain", "order_increased"):
            lab = st.bond_label(c["atoms"][0]["element"], c["atoms"][1]["element"], c["after_order"])
            labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(c["change_id"])
    r["classification"]["mapped_reaction_smiles"] = missing(ev)
    r["classification"]["bond_changes"] = facts
    r["classification"]["labels"] = f(None, list(labels.values()), ev, method="structure_derivation",
                                            derivation="Taxonomy labels are calculated from the manually transcribed partial net bond changes.")
    r["classification"]["coverage"] = inferred("partial", ev, "Dense drawings and omitted coproducts prevent complete atom mapping.")
    r["classification"]["limitations"] = [limitation, "No complete atom mapping or elementary-mechanism correspondence is asserted."]
    r["classification"]["model_votes"] = []


C_C = {"atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"}
C_H = {"atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1}
O_H = {"atoms": [{"atom_map": None, "element": "O"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1}
C_N = {"atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "N"}], "before_order": 0, "after_order": 1, "change_type": "formed"}
N_H = {"atoms": [{"atom_map": None, "element": "N"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1}
C_C_REDUCE = {"atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 2, "after_order": 1, "change_type": "order_decreased"}
C_O_BREAK = {"atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}], "before_order": 1, "after_order": 0, "change_type": "broken"}
C_C_BREAK = {"atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 1, "after_order": 0, "change_type": "broken"}


def reaction(rid, page, record_type, categories, location, name, transformation, evidence_extra=()):
    ev = [EV[page], *evidence_extra]
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, ev)
    r["location"] = f(location, {"page": page, "region": location}, ev)
    r["reaction_name"] = f(name, name, ev)
    r["transformation_type"] = inferred(transformation, ev, "Transformation type assigned from the depicted net reactant/product change.")
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
                   "unresolved": ["No Supporting Information was supplied; detailed procedure, exact workup timing, and characterization remain unavailable."]}
    r["completeness"] = "partial"
    bundle["reactions"].append(r)
    return r, ev


def outcome(r, oid, kind, targets, ev, raw, normalized, product_role=None,
            yield_type=None, basis=None, scope=None, inherited=False):
    o = st.new_outcome(oid, kind, targets, ev)
    method = "inherited_general_condition" if inherited else "original_extraction"
    deriv = "Applied from the table-wide footnote to this entry." if inherited else None
    o["value"] = f(raw, normalized, ev, method=method, derivation=deriv)
    if product_role:
        o["product_role"] = f(product_role, product_role, ev)
    if yield_type:
        o["yield_type"] = f(yield_type, yield_type, ev, method=method, derivation=deriv)
    if basis:
        o["measurement_basis"] = f(basis, basis, ev, method=method, derivation=deriv)
    if scope:
        o["scope"] = f(scope, scope, ev)
    r["outcomes"].append(o)
    return o


# Page 1: two prior-strategy illustrations and one current generic scheme.
for rid, location, name, note in [
    ("figure1b_prior_allylic_alkylation", "Figure 1B, upper-left prior-strategy arrow", "prior allylic alkylation strategy", "Prior-work class uses geometrically pure substrate and superstoichiometric organometallic reagent."),
    ("figure1b_prior_alpha_functionalization", "Figure 1B, lower-left prior-strategy arrow", "prior alpha-functionalization strategy", "Prior-work class typically requires a directing group or electron-withdrawing group."),
]:
    r, ev = reaction(rid, 1, "literature_example", ["other"], location, name, "generic C-C bond formation")
    r["notes"] = [f(note, note, ev)]
    manual_classification(r, [C_C], ev, "Generic R-group literature illustration; no concrete experiment, product yield, or complete structures are asserted.")

r, ev = reaction("figure1b_this_work_generic", 1, "generic_scheme", ["other"], "Figure 1B, blue 'this work' arrow", "CuH-catalyzed hydrocarboxylation of allenes", "branched allene hydrocarboxylation")
r["notes"] = [f("Exclusive branched-selectivity; capable of synthesizing both alpha-quaternary and tertiary acids.", "Exclusive branched-selectivity; capable of synthesizing both alpha-quaternary and tertiary acids.", ev)]
manual_classification(r, [C_C, C_H], ev, "Generic R-group scheme; ester hydrolysis and product form vary by workup.")


# Table 1 optimization entries.
table1 = [
    (1, "m_L1", "m_2a", "m_opt_prod_2a", 40, 2.0, 42, "10:90", None),
    (2, "m_L1", "m_2b", "m_opt_prod_2b", 40, 2.0, "<5", None, None),
    (3, "m_L1", "m_2c", "m_opt_prod_2c", 25, 1.5, "<5", None, "Either LiOMe (1.1 equiv) or CsOBz (1.1 equiv) was used as an additive."),
    (4, "m_L1", "m_2d", "m_opt_prod_2d", 25, 2.0, 83, "13:87", None),
    (5, conf, "m_2d", "m_opt_prod_2d", 25, 1.2, 92, "99:1", "Ligand conflict: Table 1 cell reads L1; prose identifies L2 for entry 5."),
    (6, "m_L1", "m_2d", "m_opt_prod_2d", 25, 1.0, 77, "96:4", "Entry 6 reverses the limiting-reagent convention; the exact charges are structured in participants."),
]
for entry, lig, electrophile, product, temp, allene_eq, yld, er, note in table1:
    r, ev = reaction(f"table1_entry_{entry}", 2, "specific_experiment", ["optimization"], f"Table 1, entry {entry}", "CuH-catalyzed hydrocarboxylation of allene", "branched allene hydroesterification")
    participant(r, "allene", "m_1a", "reactant", ev, charge=("equivalents", f"{allene_eq} equiv", qty(allene_eq, "equiv")), reference=electrophile, inherited=True)
    e_eq = 1.2 if entry == 6 else 1.0
    participant(r, "electrophile", electrophile, "reactant", ev, charge=("equivalents", f"{e_eq} equiv", qty(e_eq, "equiv")), reference=("m_1a" if entry == 6 else electrophile), inherited=True)
    if entry != 6:
        r["reaction_scale"] = f("0.10 mmol 2", {"reference_molecule_id": electrophile, "amount": qty(0.10, "mmol")}, ev, method="inherited_general_condition", derivation="Table 1 general condition applies to this row.")
    else:
        r["reaction_scale"] = f("0.10 mmol limiting-reagent scale; entry-6 reference changes to 1a", {"reference_molecule_id": "m_1a", "amount": qty(0.10, "mmol")}, ev,
                                  status="pending_review")
        r["review"]["unresolved"].append("Entry 6 changes the 1a/2d equivalence reference; the printed 0.10 mmol header is retained as a limiting-reagent-scale interpretation.")
    participant(r, "copper", "m_cu_oac2", "catalyst", ev, charge=("mol_percent", "5.0 mol%", qty(5.0, "mol%")), reference=("m_1a" if entry == 6 else electrophile), inherited=True)
    participant(r, "ligand", lig, "ligand", ev, charge=("mol_percent", "5.5 mol%", qty(5.5, "mol%")), reference=("m_1a" if entry == 6 else electrophile), inherited=True)
    participant(r, "silane", "m_dmms", "reagent", ev, charge=("equivalents", "3.0 equiv", qty(3.0, "equiv")), reference=("m_1a" if entry == 6 else electrophile), inherited=True)
    participant(r, "solvent", "m_thf", "solvent", ev)
    if entry == 3:
        participant(r, "additive_liome", "m_liome", "additive", ev, charge=("equivalents", "1.1 equiv", qty(1.1, "equiv")), reference=electrophile)
        participant(r, "additive_csobz", "m_csobz", "additive", ev, charge=("equivalents", "1.1 equiv", qty(1.1, "equiv")), reference=electrophile)
        r["review"]["unresolved"].append("Entry 3 reports either LiOMe or CsOBz, not simultaneous use; the table does not split them into separate yield rows.")
    participant(r, "product", product, "product", ev)
    r["solvent_conditions"]["composition"] = f("THF", ["THF"], ev, method="inherited_general_condition", derivation="Table 1 general condition.")
    r["solvent_conditions"]["concentration"] = f("0.5 M", qty(0.5, "M"), ev, method="inherited_general_condition", derivation="Table 1 general condition.")
    r["solvent_conditions"]["concentration_basis"] = f("2 (or limiting reagent for entry 6)", electrophile if entry != 6 else "m_1a", ev, method="inherited_general_condition", derivation="Table heading and entry-6 footnote.")
    add_step(r, "reaction", "reaction", [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction" and p["role"] != "product"], ev,
             operation="CuH-catalyzed hydroesterification", temperature=(f"{temp} °C", qty(temp, "°C")), time=("24 h", qty(24, "h")), inherited=True)
    yval = qty(5, "%", "lt", "<5%") if yld == "<5" else qty(yld, "%")
    outcome(r, "yield", "yield", [product], ev, f"{yld}%" if yld != "<5" else "<5%", yval, product_role="desired branched product",
            yield_type="crude reaction mixture yield", basis="1H NMR spectroscopy using 1,3,5-trimethoxybenzene as internal standard", inherited=True)
    if er:
        er_norm = {"metric": "er", "ratio": [float(x) for x in er.split(":")], "order": "as printed; enantiomer identities not assigned"}
        o = outcome(r, "er", "selectivity", [product], ev, f"{er} er", er_norm, basis="SFC analysis", scope="enantiomeric ratio", inherited=True)
        if entry == 1:
            o["value"] = f(["Table 1: 10:90 er", "prose: 90:10 er"], None, ev, status="conflict",
                             alternatives=[f("Table 1: 10:90 er", {"metric": "er", "ratio": [10, 90], "order": "table order"}, ev),
                                           f("prose: 90:10 er", {"metric": "er", "ratio": [90, 10], "order": "prose order"}, ev)])
            r["review"]["unresolved"].append("Entry 1 er order conflicts between Table 1 (10:90) and prose (90:10).")
    if note:
        r["notes"].append(f(note, note, ev))
    if entry == 5:
        r["review"]["unresolved"].append("Ligand identity conflicts between the Table 1 cell (L1) and the page-2 prose (L2).")
    manual_classification(r, [C_C, C_H], ev, "Branched hydroesterification inferred from the depicted allene-to-alpha-vinyl ester transformation.")
    r["procedure"]["original_text"] = f("Conditions: 0.10 mmol 2 (1.0 equiv), 1a (2.0 equiv), copper(II) acetate (5.0 mol%), ligand (5.5 mol%), dimethoxy(methyl)silane (3.0 equiv) in THF (0.5 M), temp, 24 h; row footnotes modify listed entries.",
                                                   "Table 1 conditions plus row-specific footnote", ev)


# Table 2 substrate scope.
scope_rows = {
    "a": (89, "99:1", []), "b": (93, "96:4", []), "c": (55, "97.5:2.5", ["b"]),
    "d": (70, "99:1", ["c"]), "e": (74, "96:4", ["b"]), "f": (77, "98:2", ["b"]),
    "g": (88, "98.5:1.5", []), "h": (70, "96.5:3.5", []), "i": (76, "97.5:2.5", ["d"]),
    "j": (87, "98:2", ["b"]), "k": (59, "99:1", ["b"]), "l": (78, "99:1", ["b", "e"]),
    "m": (75, "98.5:1.5", ["b"]), "n": (70, "94.5:5.5", ["d"]),
    "o": (70, "93:7", ["d", "f"]), "p": (72, "93:7", ["d", "f"]), "q": (80, "93:7", ["d", "f"]),
}
for letter, (yld, er, flags) in scope_rows.items():
    r, ev = reaction(f"table2_3{letter}", 3, "specific_experiment", ["substrate_scope"], f"Table 2, product 3{letter}", "CuH-catalyzed hydrocarboxylation of allenes", "branched allene hydrocarboxylation")
    substrate = f"m_1{letter}"
    product = f"m_3{letter}"
    allene_eq = 1.1 if "e" in flags else 1.2
    ligand = "m_L3" if "d" in flags else "m_L2"
    participant(r, "electrophile", "m_2d", "reactant", ev, charge=("equivalents", "1.0 equiv", qty(1.0, "equiv")), reference="m_2d", inherited=True)
    participant(r, "allene", substrate, "reactant", ev, charge=("equivalents", f"{allene_eq} equiv", qty(allene_eq, "equiv")), reference="m_2d", inherited=True)
    participant(r, "copper", "m_cu_oac2", "catalyst", ev, charge=("mol_percent", "5.0 mol%", qty(5.0, "mol%")), reference="m_2d", inherited=True)
    participant(r, "ligand", ligand, "ligand", ev, charge=("mol_percent", "5.5 mol%", qty(5.5, "mol%")), reference="m_2d", inherited=True)
    participant(r, "silane", "m_dmms", "reagent", ev, charge=("equivalents", "3.0 equiv", qty(3.0, "equiv")), reference="m_2d", inherited=True)
    solvent_mid = "m_dme" if "f" in flags else "m_thf"
    participant(r, "solvent", solvent_mid, "solvent", ev)
    participant(r, "product", product, "product", ev)
    if letter == "l":
        participant(r, "ester_intermediate", "m_3l_ester_intermediate", "product", [EV[3], EV[7]], stage="reaction")
    workup_a = scope[letter][2] == "acid"
    participant(r, "workup_nh4f", "m_nh4f", "workup_agent", ev, stage="workup")
    participant(r, "workup_meoh", "m_meoh", "workup_agent", ev, stage="workup")
    if workup_a:
        participant(r, "hydrolysis_tfa", "m_tfa", "workup_agent", ev, stage="workup")
    r["reaction_scale"] = f("0.50 mmol 2d", {"reference_molecule_id": "m_2d", "amount": qty(0.50, "mmol")}, ev, method="inherited_general_condition", derivation="Table 2 general conditions.")
    r["solvent_conditions"]["composition"] = f("DME" if "f" in flags else "THF", ["DME" if "f" in flags else "THF"], ev,
                                                    method="inherited_general_condition", derivation="Table 2 general conditions and footnote f where applicable.")
    if "f" in flags:
        r["solvent_conditions"]["total_volume"] = f("DME, 1.0 mL", qty(1.0, "mL"), ev)
    else:
        r["solvent_conditions"]["concentration"] = f("0.5 M", qty(0.5, "M"), ev, method="inherited_general_condition", derivation="Table 2 general condition.")
        r["solvent_conditions"]["concentration_basis"] = f("2d", "m_2d", ev, method="inherited_general_condition", derivation="Table 2 limiting reagent.")
    temp = 0 if "f" in flags else 40 if "b" in flags else 30 if "c" in flags else "rt"
    temp_norm = qty(temp, "°C") if isinstance(temp, int) else "rt"
    add_step(r, "reaction", "reaction", ["electrophile", "allene", "copper", "ligand", "silane", "solvent"], ev,
             operation="CuH-catalyzed hydrocarboxylation", temperature=(f"{temp} °C" if isinstance(temp, int) else "rt", temp_norm),
             time=("6-24 h", range_qty(6, 24, "h", "6-24 h")), inherited=True)
    workup_target = "ester_intermediate" if letter == "l" else "product"
    add_step(r, "nh4f_meoh_workup", "workup", ["workup_nh4f", "workup_meoh", workup_target], ev, operation="NH4F/MeOH workup")
    if letter == "l":
        page7_ev = [EV[7], EVT[7]]
        r["evidence_ids"].extend(page7_ev)
        add_step(
            r,
            "isolate_ester_intermediate",
            "purification",
            ["ester_intermediate"],
            page7_ev,
            operation="first isolating the corresponding ester in order to facilitate its purification",
        )
        r["notes"].append(f(
            "The only exception was the hydrolysis to give carboxylic acid 3l, which was carried out after first isolating the corresponding ester in order to facilitate its purification.",
            "For 3l, the corresponding ester was isolated before hydrolysis to the carboxylic acid.",
            page7_ev,
        ))
        r["review"]["unresolved"].append(
            "Footnote 101 states that the corresponding ester was isolated before hydrolysis, but does not report its structure, isolation method, amount, yield, time, or temperature."
        )
    if workup_a:
        hydrolysis_participants = ["hydrolysis_tfa", "product"]
        hydrolysis_ev = ev
        if letter == "l":
            hydrolysis_participants.insert(0, "ester_intermediate")
            hydrolysis_ev = [EV[3], EV[7], EVT[7]]
        add_step(r, "tfa_hydrolysis", "workup", hydrolysis_participants, hydrolysis_ev, operation="hydrolysis using TFA")
    outcome(r, "yield", "yield", [product], ev, f"{yld}% yield", qty(yld, "%"), product_role="isolated product",
            yield_type="average isolated yield of two runs", basis="average isolated yield of two runs", inherited=True)
    outcome(r, "er", "selectivity", [product], ev, f"{er} er", {"metric": "er", "ratio": [float(x) for x in er.split(":")], "order": "as printed; absolute enantiomer order not assigned"},
            scope="enantiomeric ratio; method not reported in main-text scheme")
    changes = [C_C, C_H] + ([O_H, C_O_BREAK] if workup_a else [])
    manual_classification(r, changes, ev, "Overall isolated-product record includes TFA ester hydrolysis for workup-A acids; workup-B entries stop at adamantyl ester.")
    r["procedure"]["original_text"] = f("Conditions: 0.50 mmol 2d (1.0 equiv), 1 (1.2 equiv), copper(II) acetate (5.0 mol%), L2 (5.5 mol%), dimethoxy(methyl)silane (3.0 equiv) in THF (0.5 M); workup A or B; row footnotes apply.",
                                                   "Table 2 general conditions plus product footnotes", ev)
    if not molecules[substrate]["smiles"]["normalized_value"] or not molecules[product]["smiles"]["normalized_value"]:
        r["review"]["unresolved"].append("Dense substrate/product drawing was retained by label without a complete SMILES candidate.")


# Scheme 1A: one multistage Curtius experiment.
r, ev = reaction("scheme1a_3a_to_6", 3, "specific_experiment", ["derivatization"], "Scheme 1A", "Curtius Rearrangement - alpha-Tertiary Amine", "Curtius rearrangement and hydrolysis")
for pid, mid, role in [("substrate", "m_3a", "reactant"), ("dp_pa", "m_dp_pazide", "reagent"), ("et3n", "m_et3n", "reagent"), ("toluene", "m_toluene", "solvent"), ("hcl", "m_hcl", "reagent"), ("product", "m_6", "product")]:
    participant(r, pid, mid, role, ev)
add_step(r, "curtius", "reaction", ["substrate", "dp_pa", "et3n", "toluene"], ev, operation="Curtius rearrangement", temperature=("110 °C", qty(110, "°C")), time=("2 h", qty(2, "h")))
add_step(r, "acid_hydrolysis", "reaction", ["hcl", "product"], ev, operation="5 M HCl treatment", temperature=("100 °C", qty(100, "°C")), time=("2 h", qty(2, "h")))
r["solvent_conditions"]["composition"] = f("toluene", ["toluene"], ev)
outcome(r, "yield", "yield", ["m_6"], ev, "75% yield", qty(75, "%"), product_role="product 6", yield_type="yield; isolation basis not printed")
outcome(r, "es", "selectivity", ["m_6"], ev, "100% es", {"metric": "enantiospecificity", "number": 100, "unit": "%", "configuration_order": "stereoretentive; absolute product descriptor not printed"}, scope="enantiospecificity; method not reported")
manual_classification(r, [C_N, N_H, N_H, C_C_BREAK], ev, "Partial Curtius net model: carboxyl carbon loss and amine formation; omitted nitrogen-containing coproducts prevent complete balance.")

# Scheme 1B hydroamination.
r, ev = reaction("scheme1b_3d_to_8", 3, "specific_experiment", ["derivatization"], "Scheme 1B", "Hydroamination - gamma-Amino Ester", "CuH-catalyzed alkene hydroamination")
for pid, mid, role in [("substrate", "m_3d", "reactant"), ("n_source", "m_reagent7", "reactant"), ("copper", "m_cu_oac2", "catalyst"), ("ligand", "m_L2", "ligand"), ("silane", "m_dmms", "reagent"), ("thf", "m_thf", "solvent"), ("nh2oh_hcl", "m_nh2oh_hcl", "reagent"), ("product", "m_8", "product")]:
    p = participant(r, pid, mid, role, ev)
    if pid == "copper": set_charge(p, "mol_percent", "5.0 mol%", qty(5.0, "mol%"), ev, "m_3d")
    if pid == "ligand": set_charge(p, "mol_percent", "5.5 mol%", qty(5.5, "mol%"), ev, "m_3d")
    if pid == "silane": set_charge(p, "equivalents", "3.0 equiv", qty(3.0, "equiv"), ev, "m_3d")
r["solvent_conditions"]["composition"] = f("THF", ["THF"], ev)
add_step(r, "hydroamination", "reaction", ["substrate", "n_source", "copper", "ligand", "silane", "thf"], ev, operation="CuH-catalyzed hydroamination", temperature=("40 °C", qty(40, "°C")), time=("24 h", qty(24, "h")))
add_step(r, "nh2oh_hcl_treatment", "reaction", ["nh2oh_hcl", "product"], ev, operation="then NH2OH·HCl; time and temperature not printed")
outcome(r, "yield", "yield", ["m_8"], ev, "53% yield", qty(53, "%"), product_role="product 8", yield_type="yield; isolation basis not printed")
outcome(r, "es", "selectivity", ["m_8"], ev, "100% es", {"metric": "enantiospecificity", "number": 100, "unit": "%", "configuration_order": "not assigned"}, scope="enantiospecificity; method not reported")
manual_classification(r, [C_N, C_H, C_C_REDUCE], ev, "Partial net hydroamination of the product-3d vinyl group.")


def indobufen_route(tag, one_r_eq, two_d_eq, final_yield):
    # Arrow 1: hydrocarboxylation.
    r1, ev = reaction(f"scheme1c_{tag}_1r_to_3r", 3, "specific_experiment", ["derivatization"], f"Scheme 1C, route {tag}, arrow 1", "CuH-catalyzed hydrocarboxylation toward (S)-Indobufen", "branched allene hydroesterification")
    participant(r1, "allene", "m_1r", "reactant", ev, charge=("equivalents", f"{one_r_eq} equiv", qty(one_r_eq, "equiv")), reference=("m_1r" if one_r_eq == 1.0 else "m_2d"))
    participant(r1, "electrophile", "m_2d", "reactant", ev, charge=("equivalents", f"{two_d_eq} equiv", qty(two_d_eq, "equiv")), reference=("m_1r" if one_r_eq == 1.0 else "m_2d"))
    for pid, mid, role, val in [("copper", "m_cu_oac2", "catalyst", 5.0), ("ligand", "m_L3", "ligand", 5.5)]:
        p = participant(r1, pid, mid, role, ev)
        set_charge(p, "mol_percent", f"{val} mol%", qty(val, "mol%"), ev, "m_1r" if one_r_eq == 1.0 else "m_2d")
    p = participant(r1, "silane", "m_dmms", "reagent", ev)
    set_charge(p, "equivalents", "3.0 equiv", qty(3.0, "equiv"), ev, "m_1r" if one_r_eq == 1.0 else "m_2d")
    for pid, mid, role, stage in [("dme", "m_dme", "solvent", "reaction"), ("nh4f", "m_nh4f", "workup_agent", "workup"), ("intermediate", "m_3r", "product", "reaction")]:
        participant(r1, pid, mid, role, ev, stage=stage)
    r1["solvent_conditions"]["composition"] = f("DME", ["DME"], ev)
    add_step(r1, "hydrocarboxylation", "reaction", ["allene", "electrophile", "copper", "ligand", "silane", "dme"], ev, operation="CuH-catalyzed hydrocarboxylation", temperature=("0 °C", qty(0, "°C")), time=("6 h", qty(6, "h")))
    add_step(r1, "nh4f_workup", "workup", ["nh4f", "intermediate"], ev, operation="NH4F workup")
    outcome(r1, "not_isolated", "other", ["m_3r"], ev, "3r, not isolated", "not isolated", product_role="intermediate 3r")
    manual_classification(r1, [C_C, C_H], ev, "Branched hydroesterification; intermediate 3r was not isolated.")

    # Arrow 2: hydrogenation.
    r2, ev = reaction(f"scheme1c_{tag}_3r_to_9", 3, "specific_experiment", ["derivatization"], f"Scheme 1C, route {tag}, arrow 2", "hydrogenation toward (S)-Indobufen", "alkene hydrogenation")
    for pid, mid, role in [("intermediate", "m_3r", "reactant"), ("pto2", "m_pto2", "catalyst"), ("hydrogen", "m_h2", "reagent"), ("meoh", "m_meoh", "solvent"), ("dcm", "m_dcm", "solvent"), ("product", "m_9", "product")]:
        participant(r2, pid, mid, role, ev)
    r2["solvent_conditions"]["composition"] = f("MeOH, DCM", ["MeOH", "DCM"], ev)
    add_step(r2, "hydrogenation", "reaction", ["intermediate", "pto2", "hydrogen", "meoh", "dcm"], ev, operation="PtO2-catalyzed hydrogenation", temperature=("rt", "rt"), time=("20 h", qty(20, "h")), pressure=("1 atm H2", "1 atm H2"))
    outcome(r2, "not_isolated", "other", ["m_9"], ev, "9, not isolated", "not isolated", product_role="intermediate 9")
    manual_classification(r2, [C_H, C_H, C_C_REDUCE], ev, "Net hydrogenation includes two carbon-hydrogen gains and C=C bond-order decrease; no new C-C sigma bond is classified.")

    # Arrow 3: ester hydrolysis.  Yield is explicitly overall from 1r to 10.
    r3, ev = reaction(f"scheme1c_{tag}_9_to_10", 3, "specific_experiment", ["derivatization"], f"Scheme 1C, route {tag}, arrow 3 and overall outcome", "hydrolysis to (S)-Indobufen", "adamantyl ester hydrolysis")
    for pid, mid, role in [("intermediate", "m_9", "reactant"), ("tfa", "m_tfa", "reagent"), ("dcm", "m_dcm", "solvent"), ("product", "m_10", "product")]:
        participant(r3, pid, mid, role, ev)
    r3["solvent_conditions"]["composition"] = f("DCM", ["DCM"], ev)
    add_step(r3, "hydrolysis", "reaction", ["intermediate", "tfa", "dcm"], ev, operation="TFA hydrolysis", temperature=("rt", "rt"), time=("4 h", qty(4, "h")))
    outcome(r3, "overall_yield", "yield", ["m_10"], ev, f"{final_yield}% overall yield", qty(final_yield, "%"), product_role="(S)-Indobufen (10)",
            yield_type="overall yield", scope=f"overall route {tag}: 1r to 10; not a last-step-only yield")
    outcome(r3, "er", "selectivity", ["m_10"], ev, "92:8 er", {"metric": "er", "ratio": [92, 8], "order": "as printed for (S)-Indobufen; minor enantiomer not separately named"}, scope="final product 10; method not reported in main text")
    r3["notes"].append(f(f"The {final_yield}% value is printed as overall yield for the multistep sequence beginning with 1r.", f"The {final_yield}% value is printed as overall yield for the multistep sequence beginning with 1r.", ev))
    manual_classification(r3, [O_H, C_O_BREAK], ev, "Ester hydrolysis classification only; the reported yield spans the entire 1r-to-10 route.")


indobufen_route("b", 1.0, 1.2, 76)
indobufen_route("c", 1.2, 1.0, 85)

# Page 4: proposed mechanism, not an experiment.
r, ev = reaction("figure2_proposed_mechanism", 4, "proposed_mechanism", ["other"], "Figure 2", "Proposed mechanism for the CuH-catalyzed hydrocarboxylation of allenes", "proposed catalytic cycle")
r["notes"] = [
    f("Allene hydrocupration gives rapidly equilibrating allylcopper species B and C; B reacts through proposed six-membered transition state D; beta-fluoride elimination gives product 3 and CuF; sigma-bond metathesis regenerates CuH.",
      "Allene hydrocupration; regioselective enantioselective addition; beta-fluoride elimination; sigma-bond metathesis.", ev),
    f("Authors state that stereochemical details of the enantio-determining transition state cannot be definitively pinpointed.",
      "stereochemical details not definitively established", ev),
]
manual_classification(r, [C_C, C_H], ev, "Overall catalytic-cycle transformation only; elementary intermediates are proposed by analogy and not entered as verified experiments.")

# Nonreaction/source-context items and six-page coverage.
bundle["source_context"] = [
    {"context_id": f"{DOC}:page3_crystal_structure_3a", "document_id": DOC, "kind": "crystallographic_structure_illustration", "evidence_ids": [EV[3]], "value": "Crystal Structure of 3a; CCDC 2050451 is cited on page 4. Not a separate synthesis experiment."},
    {"context_id": f"{DOC}:page4_si_notice", "document_id": DOC, "kind": "supplementary_information_notice", "evidence_ids": [EV[4]], "value": "Supporting Information is advertised, but no SI PDF was supplied or associated in this session."},
    {"context_id": f"{DOC}:pages5_6_references", "document_id": DOC, "kind": "bibliography", "evidence_ids": [EV[5], EV[6]], "value": "Pages 5-6 contain references only; cited papers were not imported as current-paper experiments."},
    {"context_id": f"{DOC}:page7_references", "document_id": DOC, "kind": "bibliography", "evidence_ids": [EV[7], EVT[7]], "value": "Page 7 contains references 91-100 and 102-112. These cited works were not imported as experiments of the current paper; footnote 101 is instead linked to the existing 3l substrate-scope record."},
]

page_records = {n: [] for n in range(1, 8)}
for rxn in bundle["reactions"]:
    for n in range(1, 8):
        if EV[n] in rxn["evidence_ids"]:
            page_records[n].append(rxn["reaction_id"])

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title and abstract", "finding": "Article identity and current-work summary; no concrete experiment row.", "evidence_ids": [EV[1]]},
            {"region": "Figure 1A", "finding": "Representative bioactive compounds; contextual structures, not current synthesis experiments.", "evidence_ids": [EV[1]]},
            {"region": "Figure 1B", "finding": "Two prior-strategy literature illustrations and one current generic hydrocarboxylation scheme inventoried separately.", "evidence_ids": [EV[1]]},
        ], "record_ids": page_records[1], "unresolved": ["Generic R groups prevent concrete molecular identities and complete mapping."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Table 1 reaction scheme, six rows, and footnotes", "finding": "Six independent optimization experiments transcribed with row-specific ratios, additives, yields, and er values.", "evidence_ids": [EV[2]]},
            {"region": "Table 1 ligand and electrophile drawings", "finding": "L1-L3 and 2a-2d inventoried; dense chiral ligands retained without forced SMILES.", "evidence_ids": [EV[2]]},
            {"region": "page-2 prose", "finding": "Corroborates entries and creates explicit Table-1/prose conflicts for entry-1 er order and entry-5 ligand identity.", "evidence_ids": [EV[2], EVT[2]]},
        ], "record_ids": page_records[2], "unresolved": ["Entry 3 reports either LiOMe or CsOBz in one table row; separate yields are not provided.", "SI-only screening details are unavailable."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Table 2 substrate scope and footnotes", "finding": "All 17 product panels 3a-3q transcribed as independent substrate-scope experiments.", "evidence_ids": [EV[3]]},
            {"region": "Table 2 crystal structure inset", "finding": "3a crystallographic illustration retained as source context, not a synthesis record.", "evidence_ids": [EV[3]]},
            {"region": "Scheme 1A", "finding": "One two-stage Curtius derivatization experiment.", "evidence_ids": [EV[3]]},
            {"region": "Scheme 1B", "finding": "One hydroamination derivatization experiment with subsequent NH2OH·HCl treatment.", "evidence_ids": [EV[3]]},
            {"region": "Scheme 1C and footnotes b/c", "finding": "Three arrows inventoried for each of two distinct 1r/2d charge-ratio routes; overall yields remain attached to the full 1r-to-10 span.", "evidence_ids": [EV[3]]},
        ], "record_ids": page_records[3], "unresolved": ["Experimental details and selectivity measurement methods are deferred to unavailable SI.", "Several dense scope structures remain label/drawing-based."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "Figure 2 catalytic cycle", "finding": "Author-proposed mechanism recorded as proposed_mechanism, not an experiment.", "evidence_ids": [EV[4]]},
            {"region": "conclusion and Supporting Information notice", "finding": "Narrative summary and SI availability notice; no additional concrete experiment arrow.", "evidence_ids": [EV[4]]},
            {"region": "author information, acknowledgments, beginning of references", "finding": "No synthesis experiment or analytical peak table.", "evidence_ids": [EV[4]]},
        ], "record_ids": page_records[4], "unresolved": ["SI is advertised but not supplied; no SI association was inferred.", "Proposed transition-state stereochemistry is explicitly not definitive."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "references 12-51", "finding": "Bibliography only; cited literature was not imported as this paper's experiments.", "evidence_ids": [EV[5]]},
        ], "record_ids": page_records[5], "unresolved": []},
        {"page": 6, "state": "attempted", "regions": [
            {"region": "references 51-90", "finding": "Bibliography only; no current reaction arrow, table row, or standalone characterization table.", "evidence_ids": [EV[6]]},
        ], "record_ids": page_records[6], "unresolved": []},
    ],
}

coverage_page7 = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 7,
            "state": "attempted",
            "regions": [
                {
                    "region": "left column, references 91-100 and 102-109",
                    "finding": "Bibliography entries only; cited literature was not imported as experiments of this paper. No reaction arrow, experimental table row, or standalone characterization peak table is present.",
                    "evidence_ids": [EV[7]],
                },
                {
                    "region": "left column, footnote 101 between references 100 and 102",
                    "finding": "Current-paper procedural clarification for product 3l: the corresponding ester was first isolated before hydrolysis to carboxylic acid 3l to facilitate purification. Linked to the existing Table 2 3l record rather than counted as a new experiment.",
                    "evidence_ids": [EV[7], EVT[7]],
                },
                {
                    "region": "right column, references 110-112 and page footer",
                    "finding": "Bibliography continuation and article footer only; no new current-paper experiment or analytical result.",
                    "evidence_ids": [EV[7]],
                },
            ],
            "record_ids": page_records[7],
            "unresolved": [
                "Footnote 101 does not provide the ester intermediate structure, isolation/purification method, amount, yield, time, or temperature.",
                "No Supporting Information was explicitly supplied or associated, so the missing 3l procedural details remain unavailable.",
            ],
        }
    ],
}

raw = {
    "document_id": DOC,
    "title": TITLE,
    "pages_visually_inspected": [1, 2, 3, 4, 5, 6],
    "table1_rows": [{"entry": e, "ligand": l, "electrophile": x, "temperature_C": t, "allene_equiv": ae, "yield_percent": y, "er": er, "note": n} for e, l, x, _, t, ae, y, er, n in table1],
    "table2_rows": [{"product": f"3{k}", "yield_percent": v[0], "er": v[1], "footnotes": v[2]} for k, v in scope_rows.items()],
    "scheme1": {"A": "3a to 6; 75% yield, 100% es", "B": "3d plus 7 to 8; 53% yield, 100% es", "C": "two charge-ratio routes through nonisolated 3r and 9 to (S)-Indobufen 10; 76% or 85% overall yield, both 92:8 er"},
    "limitations": ["Raw visual transcription is a machine artifact and remains needs_review.", "No SI was associated.", "Pages outside 1-6 were used only as extracted-text context and are not claimed as visually covered."],
}

raw_page7 = {
    "document_id": DOC,
    "title": TITLE,
    "pages_visually_inspected": [7],
    "page_findings": [
        {
            "page": 7,
            "regions": [
                "references 91-100 and 102-112",
                "footnote 101 concerning the 3l ester-isolation/hydrolysis order",
                "article footer",
            ],
            "new_reaction_records": [],
            "existing_reaction_records_updated": [f"{DOC}:table2_3l"],
            "nonreaction_context": "Bibliography continuation; cited papers are not experiments of this paper.",
        }
    ],
    "limitations": [
        "Machine visual transcription remains needs_review.",
        "The corresponding ester intermediate in footnote 101 is not structurally depicted or named on page 7.",
        "No SI was associated.",
        "Only original PDF page 7 is claimed as visually covered in this batch segment.",
    ],
}

st.write(RUN / "submitted.json", bundle)
st.write(RUN / "coverage/pages_0001_0006.json", coverage)
st.write(RUN / "coverage/pages_0007_0007.json", coverage_page7)
st.write(RUN / "raw_model/visual_transcription.json", raw)
st.write(RUN / "raw_model/page7_visual_transcription.json", raw_page7)
st.write(RUN / "raw_model/source_mentions.json", st.source_mentions(json.loads((RUN / "dataset.json").read_text()), ["1a", "2d", "3a", "3d", "1r", "Indobufen"], context=180))
diagnostics = {}
for mid, mol in molecules.items():
    smi = mol["smiles"]["normalized_value"]
    if smi:
        diagnostics[mid] = st.stereo_inspect(smi)
st.write(RUN / "raw_model/stereo_diagnostics.json", diagnostics)
print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]), "page_record_counts": {str(k): len(v) for k, v in page_records.items()}}, ensure_ascii=False, indent=2))
