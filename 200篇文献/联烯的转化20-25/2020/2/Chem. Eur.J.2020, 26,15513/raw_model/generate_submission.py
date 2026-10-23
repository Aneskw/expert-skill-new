#!/usr/bin/env python3
"""Build the evidence-linked v3.1 submission for doc_2e89295460c922b2.

This is a source-derived construction script.  It deliberately keeps structures
unresolved where the six supplied main-text pages do not support a safe complete
SMILES, and it never treats schema validation as chemical verification.
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_2e89295460c922b2"
RUN = ROOT / "runs/current"
DATASET = json.loads((RUN / "dataset.json").read_text())
EV = {p: f"{DOC}:page{p}:text" for p in range(1, 7)}
IMG = {p: f"{DOC}:page{p}:image" for p in range(1, 7)}


def expand_evidence(evidence):
    """Keep text and the actually inspected rendered page together.

    Text is useful for searchable wording, while structures, wedges, tables and
    plot marks are source-supported by the rendered PDF page.  Do not create a
    synthetic crop evidence ID: crops retain receipts, but the dataset's real
    evidence object is the full page image.
    """
    expanded = []
    for eid in evidence:
        expanded.append(eid)
        if eid.endswith(":text"):
            image_eid = eid[:-5] + ":image"
            if image_eid not in expanded:
                expanded.append(image_eid)
    return list(dict.fromkeys(expanded))


def q(number, unit, **extra):
    return {"number": number, "unit": unit, **extra}


def direct(raw, normalized, evidence):
    return st.fact(raw, normalized, expand_evidence(evidence))


def derived(raw, normalized, evidence, derivation, method="structure_derivation"):
    return st.fact(raw, normalized, expand_evidence(evidence), method=method, derivation=derivation)


def inherited(raw, normalized, evidence, derivation):
    return st.fact(
        raw,
        normalized,
        expand_evidence(evidence),
        method="inherited_general_condition",
        derivation=derivation,
    )


def missing(evidence=(), status="unreported"):
    return st.fact(evidence=expand_evidence(evidence), status=status)


def unreadable(raw, evidence):
    return st.fact(raw, None, expand_evidence(evidence), status="unreadable")


bundle = json.loads((RUN / "raw_model/standard_init_template.json").read_text())
molecules = {}
stereo_diagnostics = {}


def add_molecule(
    mid,
    label,
    evidence,
    *,
    name=None,
    smiles=None,
    chirality="unresolved",
    geometry="unresolved",
    representation="unresolved",
    atom_centers=None,
    double_bonds=None,
    axial=None,
    mixture=None,
    assignment=None,
    note=None,
):
    evidence = expand_evidence(evidence)
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = direct(label, label, evidence)
    m["names"]["original"] = direct(name or label, name or label, evidence)
    m["names"]["english"] = direct(name or label, name or label, evidence)
    m["names"]["abbreviations"] = missing(evidence)
    m["names"]["aliases"] = missing(evidence)
    m["cas_numbers"] = missing(evidence)
    m["external_properties"] = []
    if smiles:
        diag = st.stereo_inspect(smiles)
        stereo_diagnostics[mid] = diag
        canonical = diag["canonical_isomeric_smiles"]
        m["representation_type"] = "smiles"
        m["smiles"] = derived(
            {"source_label": label, "source_representation": "drawn structure or explicit chemical name"},
            canonical,
            evidence,
            "Candidate SMILES was transcribed from the supplied source drawing/name and canonicalized with RDKit; parsing does not verify the chemical identity.",
        )
        m["generic_template"] = None
        m["stereochemistry"]["atom_centers"] = derived(
            "standard stereo graph diagnostic", diag["atom_centers"], evidence,
            "Computed from the candidate SMILES with standard stereo; the source identity remains needs_review."
        )
        m["stereochemistry"]["double_bonds"] = derived(
            "standard stereo graph diagnostic", diag["double_bonds"], evidence,
            "Computed from the candidate SMILES with standard stereo; ordinary ring/aromatic double bonds are not E/Z assignments."
        )
        m["stereochemistry"]["isomeric_smiles"] = derived(
            "standard stereo canonical isomeric SMILES", canonical, evidence,
            "Generated from the candidate SMILES by standard stereo."
        )
        m["stereochemistry"]["identity_key"] = derived(
            "standard stereo graph key", diag["stereo_graph_key"], evidence,
            "Generated from the canonical isomeric molecular graph; this is not a sample identity."
        )
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing(evidence)
        m["generic_template"] = None
        m["stereochemistry"]["atom_centers"] = (
            derived("source drawing", atom_centers, evidence, "Stereochemical content was read from the source drawing; a complete atom-indexed structure is unavailable.")
            if atom_centers is not None else missing(evidence)
        )
        m["stereochemistry"]["double_bonds"] = (
            derived("source drawing", double_bonds, evidence, "Alkene geometry status was assessed from the source drawing; no unsupported E/Z descriptor was assigned.")
            if double_bonds is not None else missing(evidence)
        )
        m["stereochemistry"]["isomeric_smiles"] = missing(evidence)
        m["stereochemistry"]["identity_key"] = missing(evidence)
    m["stereochemistry"]["axial"] = (
        derived("source allene drawing", axial, evidence, "Allene axial stereochemistry was assessed from the drawn substitution pattern and explicit labels; absent labels were not converted to racemic or achiral.")
        if axial is not None else missing(evidence)
    )
    m["stereochemistry"]["mixture"] = (
        direct("source mixture notation", mixture, evidence) if mixture is not None else missing(evidence)
    )
    m["stereochemistry"]["assignment"] = (
        derived("source stereochemical drawing/notation", assignment, evidence, "Assignment preserves only the explicit wedge/dash or mixture notation visible in the source.")
        if assignment is not None else missing(evidence)
    )
    for key, value in (
        ("chirality_status", chirality),
        ("geometry_status", geometry),
        ("representation_status", representation),
    ):
        m["stereochemistry"][key] = derived(
            {"source_label": label, "stereochemical_review": "supplied drawing/name"}, value, evidence,
            "Status assigned conservatively after review of the supplied drawing and standard stereo diagnostics when a candidate SMILES exists.",
            method="inference",
        )
    default_note = (
        "Identity is tied to the paper label; the candidate SMILES is a needs-review transcription and parsing is not structure verification."
        if smiles else
        "Identity is tied to the paper label. A missing SMILES records unresolved graphical structure rather than absence of structure."
    )
    m["identity_note"] = direct(
        note or default_note,
        note or default_note,
        evidence,
    )
    molecules[mid] = m
    bundle["molecules"].append(m)
    return mid


def add_template(mid, label, evidence, fixed, variables, note):
    add_molecule(mid, label, evidence, name=label, chirality="unresolved", geometry="unresolved", representation="not_encoded", note=note)
    m = molecules[mid]
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": derived("generic structure in Scheme 2", fixed, evidence, "Only the invariant connectivity and labels in the generic source drawing are retained."),
        "variables": [
            {
                "label": v,
                "attachment_points": derived(f"attachment point for {v}", info["attachment"], evidence, "Attachment position follows the generic source drawing."),
                "definition": direct(info["definition"], info["definition"], evidence),
                "allowed_values": direct(info.get("allowed_raw"), info.get("allowed"), evidence) if info.get("allowed") is not None else missing(evidence),
                "constraints": direct(info.get("constraints"), info.get("constraints"), evidence) if info.get("constraints") is not None else missing(evidence),
            }
            for v, info in variables.items()
        ],
        "instance_mapping": [],
        "reason": direct(note, note, evidence),
    }
    return mid


# Generic current-work transformation (Scheme 2b).
add_template(
    "m_generic_bisallene_1", "generic bisallene 1", [EV[2]],
    {"motif": "olefin-assisted O-linked bisallene", "specified_labels": ["R", "R'"], "complete_structure": False},
    {
        "R": {"attachment": "left allene substituent", "definition": "R as drawn in Scheme 2b"},
        "Rprime": {"attachment": "right allene substituent", "definition": "R' as drawn in Scheme 2b"},
    },
    "Scheme 2b is a generic R-group drawing; no wildcard SMILES is asserted.",
)
add_template(
    "m_generic_cycle_3", "generic seven-membered product 3", [EV[2]],
    {"motif": "O-containing seven-membered carbocyclization product", "specified_labels": ["R", "R'", "X"], "complete_structure": False},
    {
        "R": {"attachment": "ring substituent", "definition": "R as drawn in Scheme 2b"},
        "Rprime": {"attachment": "ring substituent", "definition": "R' as drawn in Scheme 2b"},
        "X": {"attachment": "newly installed group", "definition": "X in Scheme 2b", "allowed_raw": "X = Bpin, Ph or CO2Me", "allowed": ["Bpin", "Ph", "CO2Me"]},
    },
    "Scheme 2b is a generic product family; alternatives are not collapsed to one concrete molecule.",
)

# Common reagents and solvents.  Each candidate SMILES is checked below.
add_molecule("m_b2pin2", "2a", [EV[2], EV[3], EV[4]], name="bis(pinacolato)diboron (B2pin2)", smiles="CC1(C)OB(B2OC(C)(C)C(C)(C)O2)OC1(C)C", chirality="achiral", geometry="not_applicable", representation="not_applicable")
add_molecule("m_pd_oac2", "Pd(OAc)2", [EV[2], EV[3], EV[4]], name="palladium(II) acetate", smiles="CC(=O)[O-].CC(=O)[O-].[Pd+2]", chirality="achiral", geometry="not_applicable", representation="not_applicable")
add_molecule("m_o2", "O2", [EV[1], EV[2], EV[3], EV[4], EV[5]], name="molecular oxygen", smiles="O=O", chirality="achiral", geometry="not_applicable", representation="not_applicable")
add_molecule("m_acetone", "acetone", [EV[2], EV[3], EV[4]], name="acetone", smiles="CC(C)=O", chirality="achiral", geometry="not_applicable", representation="not_applicable")
add_molecule("m_bq", "BQ", [EV[2], EV[4]], name="p-benzoquinone", smiles="O=C1C=CC(=O)C=C1", chirality="achiral", geometry="not_applicable", representation="not_applicable")
add_molecule("m_hq", "HQ", [EV[1], EV[2], EV[3], EV[4]], name="hydroquinone", smiles="Oc1ccc(O)cc1", chirality="achiral", geometry="not_applicable", representation="not_applicable")
add_molecule("m_phboronic", "PhB(OH)2", [EV[4]], name="phenylboronic acid", smiles="OB(O)c1ccccc1", chirality="achiral", geometry="not_applicable", representation="not_applicable")
add_molecule("m_dmso", "DMSO", [EV[4]], name="dimethyl sulfoxide", smiles="CS(C)=O", chirality="achiral", geometry="not_applicable", representation="not_applicable")
add_molecule("m_meoh", "MeOH", [EV[4]], name="methanol", smiles="CO", chirality="achiral", geometry="not_applicable", representation="not_applicable")
add_molecule("m_co", "CO", [EV[4]], name="carbon monoxide", smiles="[C-]#[O+]", chirality="achiral", geometry="not_applicable", representation="not_applicable")
add_molecule("m_dce", "DCE", [EV[4]], name="1,2-dichloroethane", smiles="ClCCCl", chirality="achiral", geometry="not_applicable", representation="not_applicable")
add_molecule("m_thf", "THF", [EV[4]], name="tetrahydrofuran", smiles="C1CCOC1", chirality="achiral", geometry="not_applicable", representation="not_applicable")
add_molecule("m_dmbq", "DMBQ", [EV[6]], name="2,6-dimethyl-p-benzoquinone", smiles="CC1=CC(=O)C(C)=CC1=O", chirality="achiral", geometry="not_applicable", representation="not_applicable")
add_molecule("m_anisole", "anisole", [EV[2], EV[3]], name="anisole (NMR internal standard)", smiles="COc1ccccc1", chirality="achiral", geometry="not_applicable", representation="not_applicable", note="Anisole is a characterization-stage internal standard, not asserted as a reaction-time reagent.")

# Catalyst synthesis inputs and products.
add_molecule("m_salicylaldehyde_hq", "salicylaldehyde-hydroquinone", [EV[2]], name="salicylaldehyde-hydroquinone ligand precursor", chirality="unresolved", geometry="unresolved", representation="not_encoded")
add_molecule("m_opda", "o-phenylenediamine", [EV[2]], name="o-phenylenediamine", smiles="Nc1ccccc1N", chirality="achiral", geometry="not_applicable", representation="not_applicable")
add_molecule("m_4cl_opda", "4-chloro-o-phenylenediamine", [EV[2]], name="4-chloro-o-phenylenediamine", smiles="Nc1ccc(Cl)cc1N", chirality="achiral", geometry="not_applicable", representation="not_applicable")
add_molecule("m_4ome_opda", "4-methoxy-o-phenylenediamine", [EV[2]], name="4-methoxy-o-phenylenediamine", smiles="COc1ccc(N)c(N)c1", chirality="achiral", geometry="not_applicable", representation="not_applicable")
add_molecule("m_en", "ethylenediamine", [EV[2]], name="ethylenediamine", smiles="NCCN", chirality="achiral", geometry="not_applicable", representation="not_applicable")
add_molecule("m_trans_dach_rac", "(±)-trans-1,2-diaminocyclohexane", [EV[2]], name="racemic trans-1,2-diaminocyclohexane", chirality="racemic", geometry="not_applicable", representation="not_encoded", atom_centers=[{"source": "opposite wedge/dash at C1/C2", "absolute_configuration": "not assigned"}], mixture={"type": "racemate", "ratio": "not reported"}, assignment={"relative": "trans", "absolute": "not assigned", "source_mark": "(±)"})
add_molecule("m_tripod_diamine", "N-(3-aminopropyl)-N-methylpropane-1,3-diamine", [EV[2]], name="N-(3-aminopropyl)-N-methylpropane-1,3-diamine", smiles="CN(CCCN)CCCN", chirality="achiral", geometry="not_applicable", representation="not_applicable")
for metal, smi in {
    "co": "CC(=O)[O-].CC(=O)[O-].[Co+2]",
    "fe": "CC(=O)[O-].CC(=O)[O-].[Fe+2]",
    "ni": "CC(=O)[O-].CC(=O)[O-].[Ni+2]",
    "cu": "CC(=O)[O-].CC(=O)[O-].[Cu+2]",
}.items():
    add_molecule(f"m_{metal}_oac2", f"{metal.capitalize()}(OAc)2", [EV[2]], name=f"{metal.capitalize()}(OAc)2", smiles=smi, chirality="achiral", geometry="not_applicable", representation="not_applicable")
for i in range(1, 10):
    kwargs = {}
    if i == 8:
        kwargs = dict(chirality="racemic", geometry="not_applicable", representation="not_encoded", mixture={"type": "racemic catalyst from (±)-trans diamine", "ratio": "not reported"}, assignment={"relative": "trans ligand backbone", "absolute": "not assigned"})
    else:
        kwargs = dict(chirality="unresolved", geometry="not_applicable", representation="not_encoded")
    add_molecule(f"m_cat{i}", f"Cat. {i}", [EV[2], EV[3]] if i == 1 else [EV[2], EV[3]], name="Co(salophen)-HQ" if i == 1 else f"hybrid catalyst Cat. {i}", **kwargs)
add_molecule("m_co_salophen", "Co(salophen)", [EV[2], EV[4]], name="cobalt salophen complex", chirality="unresolved", geometry="not_applicable", representation="not_encoded")

# Main-text labelled substrates and products.  Full SMILES are intentionally not guessed.
for letter in "abcdefghijklmno":
    ch = "specified" if letter == "m" else "unspecified"
    assignment = {"embedded_fragment": "estrone-derived stereocenters shown by wedge/dash", "absolute_configuration": "not independently assigned"} if letter == "m" else None
    add_molecule(
        f"m_1{letter}", f"1{letter}", [EV[3]], name=f"bisallene 1{letter}",
        chirality="unspecified", geometry="not_applicable", representation="not_encoded",
        axial=[{"feature": "allene unit(s)", "configuration": "not assigned in source"}],
        note="Specific source drawing reviewed; complete atom-level transcription is deferred because the main-text figure is insufficient for safe stereochemically complete SMILES.",
    )
    add_molecule(
        f"m_3{letter}", f"3{letter}", [EV[3]], name=f"seven-membered borylation product 3{letter}",
        chirality=ch, geometry="unspecified", representation="not_encoded",
        double_bonds=[{"feature": "drawn alkene(s)", "configuration": "no E/Z label in supplied source"}],
        assignment=assignment,
        note="Specific product drawing reviewed; terminal/exocyclic alkene geometry is not promoted to an E/Z assignment, and no complete SMILES is asserted.",
    )
# Extend occurrences for labels repeated elsewhere.
for mid in ("m_1a", "m_3a"):
    extra = [EV[2], EV[3], EV[4]]
    molecules[mid]["occurrences"][0]["evidence_ids"] = extra
    molecules[mid]["occurrences"][0]["compound_label"]["evidence_ids"] = extra
for mid in ("m_1g", "m_3g"):
    extra = [EV[3], EV[4]]
    molecules[mid]["occurrences"][0]["evidence_ids"] = extra
    molecules[mid]["occurrences"][0]["compound_label"]["evidence_ids"] = extra
add_molecule("m_1s", "1s", [EV[4]], name="sterically hindered bisallene 1s", chirality="unspecified", geometry="not_applicable", representation="not_encoded", axial=[{"feature": "allene unit(s)", "configuration": "not assigned"}])
add_molecule("m_1aprime", "1a'", [EV[4]], name="ethyl analogue 1a' lacking the pending vinyl group", chirality="unspecified", geometry="not_applicable", representation="not_encoded", axial=[{"feature": "allene unit(s)", "configuration": "not assigned"}])
add_molecule("m_nts_analogue", "N analogue of 3g precursor", [EV[3], EV[4]], name="NTs analogue attempted for an N-containing seven-membered heterocycle", chirality="unresolved", geometry="unresolved", representation="not_encoded")
add_molecule("m_b2neop2", "bis(neopentyl glycolato)diboron", [EV[3], EV[4]], name="bis(neopentyl glycolato)diboron", chirality="achiral", geometry="not_applicable", representation="not_encoded")
add_molecule("m_b2hex2", "bis(hexylene glycolato)diboron", [EV[3], EV[4]], name="bis(hexylene glycolato)diboron", chirality="unresolved", geometry="not_applicable", representation="not_encoded", note="The exact hexylene-glycolato isomer is not named in the supplied pages; no unique structure is invented.")
for label, page, name in [("4", 4, "arylated carbocyclization product 4"), ("5", 4, "methoxycarbonylated carbocyclization product 5"), ("6", 4, "beta-H-elimination carbocyclization product 6")]:
    add_molecule(f"m_{label}", label, [EV[page]], name=name, chirality="unspecified", geometry="unspecified", representation="not_encoded", double_bonds=[{"feature": "drawn alkene(s)", "configuration": "not explicitly assigned"}])


def participant(r, pid, mid, role, evidence, *, equiv=None, molpct=None, amount=None, stage="reaction", eqref=None):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    if equiv is not None:
        p["charges"]["equivalents"] = direct(f"{equiv} equiv", q(equiv, "equiv"), evidence)
    if molpct is not None:
        p["charges"]["mol_percent"] = direct(f"{molpct} mol%", q(molpct, "mol%"), evidence)
    if amount is not None:
        number, unit = amount
        p["charges"]["amount"] = direct(f"{number} {unit}", q(number, unit), evidence)
    if equiv is not None or molpct is not None:
        p["equivalence_reference"] = direct(eqref, eqref, evidence)
    r["participants"].append(p)
    return p


def base_reaction(rid, evidence, record_type, categories, location, reaction_name, transformation):
    evidence = expand_evidence(evidence)
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = direct(categories, categories, evidence)
    r["location"] = direct(location, location, evidence)
    r["reaction_name"] = direct(reaction_name, reaction_name, evidence)
    r["transformation_type"] = derived(
        reaction_name, transformation, evidence,
        "Normalized transformation description derived from the source reaction title/caption without claiming an elementary mechanism.",
        method="inference",
    )
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    bundle["reactions"].append(r)
    return r


def set_procedure(r, raw, evidence, general_ref=None):
    r["procedure"]["original_text"] = direct(raw, raw, evidence)
    if general_ref:
        r["procedure"]["general_procedure_reference"] = direct(general_ref, general_ref, evidence)


def set_solvent(r, mids, evidence, concentration=None, basis=None):
    r["solvent_conditions"]["composition"] = direct(mids, mids, evidence)
    if concentration is not None:
        r["solvent_conditions"]["concentration"] = direct(f"{concentration} M", q(concentration, "M"), evidence)
        r["solvent_conditions"]["concentration_basis"] = direct(basis, basis, evidence)


def add_reaction_step(r, evidence, participant_ids, *, operation, temp=None, time=None, atmosphere=None, pressure=None, sid="reaction"):
    s = st.new_step(sid, evidence)
    s["participant_ids"] = participant_ids
    s["operation"] = direct(operation, operation, evidence)
    if temp is not None:
        s["temperature"] = direct(f"{temp} °C", q(temp, "°C"), evidence)
    if time is not None:
        s["time"] = direct(f"{time} h", q(time, "h"), evidence)
    if atmosphere is not None:
        s["atmosphere"] = direct(atmosphere, atmosphere, evidence)
    if pressure is not None:
        s["pressure"] = direct(pressure, pressure, evidence)
    r["steps"].append(s)
    return s


def add_yield(r, oid, product_mid, percent, evidence, basis=None, yield_type=None):
    o = st.new_outcome(oid, "yield", [product_mid], evidence)
    o["value"] = direct(f"{percent}%", q(percent, "%"), evidence)
    if yield_type:
        o["yield_type"] = direct(yield_type, yield_type, evidence)
    if basis:
        o["measurement_basis"] = direct(basis, basis, evidence)
    o["product_role"] = direct("desired product", "desired product", evidence)
    r["outcomes"].append(o)
    return o


def add_no_reaction(r, oid, evidence, raw, basis=None):
    o = st.new_outcome(oid, "no_reaction", [], evidence)
    o["value"] = direct(raw, True, evidence)
    o["scope"] = direct(raw, raw, evidence)
    if basis:
        o["measurement_basis"] = direct(basis, basis, evidence)
    r["outcomes"].append(o)


def add_recovery(r, oid, target_mid, value, evidence, basis=None):
    o = st.new_outcome(oid, "recovery", [target_mid], evidence)
    if isinstance(value, tuple):
        number, qualifier, raw = value
        o["value"] = direct(raw, q(number, "%", qualifier=qualifier, raw_text=raw), evidence)
    else:
        o["value"] = direct(f"{value}%", q(value, "%"), evidence)
    if basis:
        o["measurement_basis"] = direct(basis, basis, evidence)
    o["scope"] = direct("recovered starting material", "recovered starting material", evidence)
    r["outcomes"].append(o)


def set_classification(r, evidence, changes, limitation):
    normalized = []
    labels = {}
    for i, (a, b, before, after, kind, source) in enumerate(changes, 1):
        cid = f"b{i}"
        ch = {"change_id": cid, "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}], "before_order": before, "after_order": after, "change_type": kind}
        normalized.append(ch)
        if kind in ("formed", "order_increased", "hydrogen_gain"):
            lab = st.bond_label(a, b, after)
            labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(cid)
        r["classification"]["bond_changes"].append(derived(source, ch, evidence, "Net bond change inferred from the highlighted/drawn reactant-to-product transformation; atom mapping is incomplete."))
    r["classification"]["labels"] = derived(
        "formed/highlighted bonds in source drawing", list(labels.values()), evidence,
        "Taxonomy labels are calculated from the explicitly retained net formed/promoted bonds only.",
    )
    r["classification"]["coverage"] = derived(
        "source structure comparison", "partial", evidence,
        "Only unambiguous highlighted/new connections are classified; allene bond-order shifts, H loss/gain, and complete atom mapping remain unresolved.",
        method="inference",
    )
    r["classification"]["limitations"] = [limitation, "No complete atom mapping is asserted.", "model_votes intentionally left empty."]


def set_no_transform_classification(r, evidence, reason, coverage="not_applicable"):
    r["classification"]["labels"] = derived(reason, [], evidence, "No observed product transformation supports a formed-bond label.", method="inference")
    r["classification"]["coverage"] = derived(reason, coverage, evidence, "The source reports no reaction/product, or the record is non-preparative; expected products are not classified.", method="inference")
    r["classification"]["limitations"] = [reason, "model_votes intentionally left empty."]


def set_unresolved_classification(r, evidence, reason):
    r["classification"]["coverage"] = derived(reason, "unresolved", evidence, "The supplied source does not support a complete net bond inventory for this record.", method="inference")
    r["classification"]["limitations"] = [reason, "No bond labels asserted.", "model_votes intentionally left empty."]


# Scheme 2b generic transformation.
r = base_reaction("scheme2b_generic", [EV[2]], "generic_scheme", ["other"], {"page": 2, "region": "Scheme 2b"}, "palladium-catalyzed oxidative carbocyclization", "generic oxidative carbocyclization family")
participant(r, "generic_substrate", "m_generic_bisallene_1", "reactant", [EV[2]])
participant(r, "generic_product", "m_generic_cycle_3", "product", [EV[2]])
r["notes"].append(direct("X = Bpin, Ph or CO2Me; B2pin2, PhB(OH)2 or CO/MeOH shown as trapping systems.", "Alternative trapping systems are represented by later specific records.", [EV[2]]))
set_unresolved_classification(r, [EV[2]], "Generic alternatives and R groups prevent one complete mapped net transformation.")
r["review"]["unresolved"] = ["Generic R/R'/X structures are retained as templates; no single concrete reaction is implied."]


# Scheme 3 catalyst preparations: one record per product/yield row.
cat_rows = [
    (1, "m_opda", "m_co_oac2", 85), (2, "m_opda", "m_fe_oac2", 85),
    (3, "m_opda", "m_ni_oac2", 90), (4, "m_opda", "m_cu_oac2", 89),
    (5, "m_4cl_opda", "m_co_oac2", 95), (6, "m_4ome_opda", "m_co_oac2", 79),
    (7, "m_en", "m_co_oac2", 78), (8, "m_trans_dach_rac", "m_co_oac2", 91),
    (9, "m_tripod_diamine", "m_co_oac2", 80),
]
for idx, diamine, metal, yld in cat_rows:
    ev = [EV[2]]
    r = base_reaction(f"scheme3_cat{idx}", ev, "specific_experiment", ["other"], {"page": 2, "region": f"Scheme 3 Cat. {idx} row"}, "synthesis of Schiff base-hydroquinone hybrid catalyst", "two-stage condensation and metal coordination")
    participant(r, "ligand_precursor", "m_salicylaldehyde_hq", "reactant", ev, equiv=2, eqref=diamine)
    participant(r, "diamine", diamine, "reactant", ev)
    participant(r, "metal_acetate", metal, "reagent", ev)
    participant(r, "catalyst_product", f"m_cat{idx}", "product", ev)
    s1 = add_reaction_step(r, ev, ["ligand_precursor", "diamine"], operation="condensation of two equivalents of salicylaldehyde-hydroquinone with the indicated diamine", sid="condensation")
    s2 = add_reaction_step(r, ev, ["metal_acetate", "catalyst_product"], operation="coordination of the resulting unpurified salen-hydroquinone ligand with M(OAc)2", sid="coordination")
    add_yield(r, "catalyst_yield", f"m_cat{idx}", yld, ev)
    set_procedure(r, "Scheme 3: 1. diamine; 2. M(OAc)2; resulting ligand used directly without further purification.", ev)
    set_unresolved_classification(r, ev, "The intermediate ligand composition and metal coordination representation are not fully specified on the supplied page.")
    r["review"]["unresolved"] = ["Temperature, time, solvent, metal-acetate charge, workup, purification and yield measurement basis are not reported on the supplied main-text page.", "Complete catalyst structure/coordination and full net bond inventory require SI review."]


def build_borylation(
    rid,
    evidence,
    location,
    substrate,
    product,
    *,
    categories,
    diboron="m_b2pin2",
    diboron_equiv=1.3,
    pd=5,
    etm=(),
    yld=None,
    yield_basis=None,
    yield_type=None,
    scale=None,
    b2_amount=None,
    recovery=None,
    no_reaction_raw=None,
    cat1_loading=None,
    general_ref=None,
):
    ev = list(dict.fromkeys(evidence))
    r = base_reaction(rid, ev, "specific_experiment", categories, location, "aerobic borylative carbocyclization", "oxidative borylative carbocyclization to a seven-membered ring")
    participant(r, "substrate", substrate, "reactant", ev, equiv=1.0, amount=scale, eqref=substrate)
    if product is not None:
        participant(r, "product", product, "product", ev)
    participant(r, "diboron", diboron, "reagent", ev, equiv=diboron_equiv, amount=b2_amount, eqref=substrate)
    participant(r, "pd_catalyst", "m_pd_oac2", "catalyst", ev, molpct=pd, eqref=substrate)
    for j, (mid, loading, role) in enumerate(etm, 1):
        participant(r, f"etm{j}", mid, role, ev, molpct=loading, eqref=substrate)
    participant(r, "oxygen", "m_o2", "reagent", ev)
    participant(r, "solvent", "m_acetone", "solvent", ev)
    if scale:
        r["reaction_scale"] = direct(f"{scale[0]} {scale[1]} of substrate", {"reference_molecule_id": substrate, "amount": q(scale[0], scale[1])}, ev)
    set_solvent(r, ["m_acetone"], ev, 0.1, substrate)
    add_reaction_step(r, ev, [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"], operation="stir under the stated aerobic carbocyclization conditions", temp=25, time=30, atmosphere="O2", pressure="1 atm")
    set_procedure(r, "bisallene, diboron reagent, Pd(OAc)2 and specified ETM(s) in 0.1 M acetone under O2 (1 atm) at 25 °C for 30 h", ev, general_ref=general_ref)
    if yld is not None and product is not None:
        add_yield(r, "product_yield", product, yld, ev, basis=yield_basis, yield_type=yield_type)
    if no_reaction_raw:
        add_no_reaction(r, "no_reaction", ev, no_reaction_raw, basis=yield_basis)
    if recovery is not None:
        add_recovery(r, "starting_material_recovery", substrate, recovery, ev)
    if product is not None and (yld is None or yld > 0):
        set_classification(r, ev, [("C", "C", 0, 1, "formed", "new ring-closing C-C bond highlighted in source"), ("C", "B", 0, 1, "formed", "new C-B bond to BX group shown in product")], "Hydrogen transfer and all allene-to-alkene bond-order changes are not atom mapped.")
    else:
        set_no_transform_classification(r, ev, no_reaction_raw or "0% product yield; no observed transformation classified")
    if not yield_basis and yld is not None:
        r["review"]["unresolved"].append("Yield measurement basis is not reported for this record on the supplied pages.")
    r["review"]["unresolved"].append("Complete reactant/product SMILES and atom mapping remain unresolved from the main-text figure; record is partial.")
    return r


# Table 1 rows.  Entry 5 also cites Figure 1 because it is the same condition/yield.
table_rows = [
    (1, [], None, None, 93, "no desired borylation product observed (0% yield)"),
    (2, [("m_bq", 20, "additive")], "m_3a", 15, None, None),
    (3, [("m_co_salophen", 10, "catalyst"), ("m_bq", 20, "additive")], "m_3a", 48, None, None),
    (4, [("m_co_salophen", 10, "catalyst"), ("m_hq", 20, "additive")], "m_3a", 35, None, None),
    (5, [("m_cat1", 10, "catalyst")], "m_3a", 78, (1, "lt", "<1%"), None),
    (6, [("m_bq", 100, "reagent")], "m_3a", 78, (1, "lt", "<1%"), None),
]
for entry, etms, prod, yld, recov, no_rxn in table_rows:
    ev = [EV[2], EV[3]] if entry == 5 else [EV[2]]
    build_borylation(
        f"table1_entry{entry}", ev,
        {"page": 2, "region": f"Table 1 entry {entry}", "duplicate_evidence": "Figure 1 Cat.1" if entry == 5 else None},
        "m_1a", prod, categories=["optimization"], etm=etms, yld=yld,
        yield_basis="1H NMR using anisole as internal standard",
        scale=(0.1, "mmol"), b2_amount=(0.13, "mmol"), recovery=recov,
        no_reaction_raw=no_rxn, general_ref="Table 1 footnote a",
    )


# Figure 1 catalyst screen excluding Cat.1, already represented by Table 1 entry 5.
fig1 = {2: 23, 3: 0, 4: 5, 5: 77, 6: 74, 7: 58, 8: 65, 9: 63}
for idx, yld in fig1.items():
    build_borylation(
        f"figure1_cat{idx}", [EV[3]], {"page": 3, "region": f"Figure 1 Cat. {idx} bar"},
        "m_1a", None if yld == 0 else "m_3a", categories=["optimization"],
        etm=[(f"m_cat{idx}", 10, "catalyst")], yld=None if yld == 0 else yld,
        yield_basis="1H NMR using anisole as internal standard", scale=(0.1, "mmol"),
        b2_amount=(0.13, "mmol"), no_reaction_raw="0% yield of 3a" if yld == 0 else None,
        general_ref="Figure 1 caption",
    )


# Scheme 4 scope: 15 distinct substrate/product experiments.
scope_yields = dict(zip("abcdefghijklmno", [75, 78, 60, 30, 59, 65, 69, 52, 66, 72, 70, 75, 78, 40, 52]))
for letter, yld in scope_yields.items():
    diboron = "m_b2neop2" if letter == "n" else "m_b2hex2" if letter == "o" else "m_b2pin2"
    loading = 15 if letter == "f" else 10
    r = build_borylation(
        f"scheme4_1{letter}_to_3{letter}", [EV[3]], {"page": 3, "region": f"Scheme 4 product 3{letter}"},
        f"m_1{letter}", f"m_3{letter}", categories=["substrate_scope"], diboron=diboron,
        etm=[("m_cat1", loading, "catalyst")], yld=yld, yield_basis="isolated yield", yield_type="isolated",
        general_ref="Scheme 4 caption and footnote a" if letter == "f" else "Scheme 4 caption",
    )
    if letter == "f":
        r["notes"].append(direct("[a] 15 mol% Cat. 1 was added.", "Cat. 1 loading exception is structured on the participant.", [EV[3]]))


# Limitations and control experiments described around Scheme 4/Eq. 1.
r = base_reaction("nts_analogue_unsuccessful", [EV[3], EV[4]], "specific_experiment", ["control"], {"pages": [3, 4], "region": "paragraph spanning pages 3–4"}, "attempted synthesis of N-containing seven-membered heterocycle", "unsuccessful oxidative carbocyclization")
participant(r, "substrate", "m_nts_analogue", "reactant", [EV[3], EV[4]])
add_reaction_step(r, [EV[3], EV[4]], ["substrate"], operation="attempted reaction; exact procedure and conditions are referred to unavailable Supporting Information")
add_no_reaction(r, "no_product", [EV[3], EV[4]], "attempts were unsuccessful")
set_procedure(r, "Attempts using the nitrogen analogue (O in the ring replaced by NTs) were unsuccessful; details referred to Supporting Information.", [EV[3], EV[4]])
set_no_transform_classification(r, [EV[3], EV[4]], "No product was reported; intended bonds are not classified.")
r["review"]["unresolved"] = ["Exact substrate structure, charges, conditions, number of attempts, detection method and SI procedure are unavailable."]

r = build_borylation(
    "scheme4_1s_no_borylation", [EV[3], EV[4]], {"pages": [3, 4], "region": "Scheme 4 caption and continuation paragraph above Eq. (1)"},
    "m_1s", None, categories=["control"], etm=[("m_cat1", 10, "catalyst")],
    no_reaction_raw="the seven-membered borylation product from 1s was not observed",
    general_ref="Scheme 4 standard conditions inherited by the immediately continuing scope limitation",
)
r["notes"].append(direct("Authors attribute the result to steric hindrance of substrate 1s.", "author interpretation: steric hindrance", [EV[4]]))
# The negative result is a continuation of the Scheme 4 scope paragraph.  Its
# conditions are inherited from the caption on page 3, not directly restated
# beside the page-4 sentence.
inherit_reason_1s = "Inherited from the Scheme 4 caption because the 1s statement is the immediately continuing scope limitation; the page-4 sentence directly supports only substrate identity and the no-product observation."
for p in r["participants"]:
    if p["participant_id"] != "substrate":
        for field in ("equivalents", "mol_percent", "amount"):
            fact = p["charges"][field]
            if fact["normalized_value"] is not None:
                p["charges"][field] = inherited(fact["raw_value"], fact["normalized_value"], [EV[3], EV[4]], inherit_reason_1s)
for field in ("composition", "concentration", "concentration_basis"):
    fact = r["solvent_conditions"][field]
    if fact["normalized_value"] is not None:
        r["solvent_conditions"][field] = inherited(fact["raw_value"], fact["normalized_value"], [EV[3], EV[4]], inherit_reason_1s)
for field in ("temperature", "time", "atmosphere", "pressure"):
    fact = r["steps"][0][field]
    r["steps"][0][field] = inherited(fact["raw_value"], fact["normalized_value"], [EV[3], EV[4]], inherit_reason_1s)
r["review"]["unresolved"].append("The page-4 sentence does not restate charges or conditions; their association is an explicit general-condition inheritance from Scheme 4, not a direct page-4 extraction.")

r = build_borylation(
    "eq1_1aprime_control", [EV[3], EV[4]], {"pages": [3, 4], "region": "Eq. (1), with standard-condition inheritance from Scheme 4"},
    "m_1aprime", None, categories=["control", "mechanistic_study"], etm=[("m_cat1", 10, "catalyst")],
    no_reaction_raw="3a', 0%; attempted reaction did not give any product",
    general_ref="Eq. (1) standard conditions",
)
r["notes"].append(direct("The vinyl group in 1a was replaced by an ethyl group in 1a'.", "control tests the necessity of the pending olefin", [EV[4]]))
# Eq. (1) explicitly prints Pd/Cat.1, acetone, O2, r.t. and 30 h.  The
# concentration, O2 pressure and reagent/substrate equivalents are supplied
# only by the phrase "standard conditions" and the Scheme 4 caption.
r["steps"][0]["temperature"] = direct("r.t.", "r.t.", [EV[4]])
inherit_reason_eq1 = "Inherited from the Scheme 4 caption through the page-4 statement that 1a' was reacted under standard conditions; Eq. (1) itself omits this numeric detail."
for pid, fields in {
    "substrate": ("equivalents",),
    "diboron": ("equivalents",),
}.items():
    p = next(x for x in r["participants"] if x["participant_id"] == pid)
    for field in fields:
        fact = p["charges"][field]
        p["charges"][field] = inherited(fact["raw_value"], fact["normalized_value"], [EV[3], EV[4]], inherit_reason_eq1)
for field in ("concentration", "concentration_basis"):
    fact = r["solvent_conditions"][field]
    r["solvent_conditions"][field] = inherited(fact["raw_value"], fact["normalized_value"], [EV[3], EV[4]], inherit_reason_eq1)
pressure_fact = r["steps"][0]["pressure"]
r["steps"][0]["pressure"] = inherited(pressure_fact["raw_value"], pressure_fact["normalized_value"], [EV[3], EV[4]], inherit_reason_eq1)
r["review"]["unresolved"].append("Eq. (1) directly shows room temperature rather than a numeric 25 °C; inherited numeric details are limited to those omitted by the equation but supplied by Scheme 4 standard conditions.")


def scheme5_reaction(rid, substrate, product, evidence, yld, participants_spec, solvent, concentration, atmosphere, pressure, changes, name, transformation):
    r = base_reaction(rid, evidence, "specific_experiment", ["other"], {"page": 4, "region": f"Scheme 5 {rid[-1]}"}, name, transformation)
    participant(r, "substrate", substrate, "reactant", evidence, equiv=1.0, eqref=substrate)
    participant(r, "product", product, "product", evidence)
    for pid, mid, role, equiv, molpct in participants_spec:
        participant(r, pid, mid, role, evidence, equiv=equiv, molpct=molpct, eqref=substrate)
    participant(r, "solvent", solvent, "solvent", evidence)
    set_solvent(r, [solvent], evidence, concentration, substrate)
    add_reaction_step(r, evidence, [p["participant_id"] for p in r["participants"]], operation=name, temp=25, time=30, atmosphere=atmosphere, pressure=pressure)
    set_procedure(r, f"Scheme 5 caption conditions for {rid[-1]}", evidence)
    add_yield(r, "product_yield", product, yld, evidence)
    set_classification(r, evidence, changes, "Only the clearly drawn new ring/termination bonds are classified; complete atom mapping and H balance are unresolved.")
    r["review"]["unresolved"] = ["Yield measurement basis is not stated in the supplied Scheme 5 caption.", "Workup, purification, detailed charges and complete structures require SI review."]
    return r


scheme5_reaction(
    "scheme5a", "m_1a", "m_4", [EV[4]], 48,
    [("phenylboronic", "m_phboronic", "reagent", 1.3, None), ("pd", "m_pd_oac2", "catalyst", None, 5), ("cat1", "m_cat1", "catalyst", None, 10), ("oxygen", "m_o2", "reagent", None, None)],
    "m_acetone", 0.1, "O2", "1 atm",
    [("C", "C", 0, 1, "formed", "new ring-closing C-C bond highlighted in yellow"), ("C", "C", 0, 1, "formed", "new aryl C-C bond highlighted in red")],
    "arylative carbocyclization", "oxidative carbocyclization-arylation",
)
scheme5_reaction(
    "scheme5b", "m_1a", "m_5", [EV[4]], 78,
    [("pd", "m_pd_oac2", "catalyst", None, 5), ("dmso", "m_dmso", "additive", None, 20), ("bq", "m_bq", "reagent", 1.5, None), ("meoh", "m_meoh", "reagent", 5.0, None), ("co", "m_co", "reagent", None, None)],
    "m_dce", 0.1, "CO", "1 atm",
    [("C", "C", 0, 1, "formed", "new ring-closing C-C bond highlighted in yellow"), ("C", "C", 0, 1, "formed", "new ring-to-carbonyl-carbon C-C bond shown for CO2Me installation")],
    "carbonylative carbocyclization", "oxidative carbocyclization-methoxycarbonylation",
)
scheme5_reaction(
    "scheme5c", "m_1s", "m_6", [EV[4]], 79,
    [("pd", "m_pd_oac2", "catalyst", None, 5), ("cat1", "m_cat1", "catalyst", None, 10), ("oxygen", "m_o2", "reagent", None, None)],
    "m_thf", 0.1, "O2", "1 atm",
    [("C", "C", 0, 1, "formed", "new ring-closing C-C bond highlighted in yellow"), ("C", "C", 1, 2, "order_increased", "new alkene from beta-H elimination shown in product")],
    "beta-H-eliminative carbocyclization", "intramolecular oxidative coupling followed by beta-H elimination",
)


# Figure 2 contains six real reaction-profile experiments.  The graph is retained
# without inventing exact point values; the caption explicitly defers details to SI.
profile_specs = [
    ("fig2a_cat1", "m_1g", "m_3g", "m_cat1", "Co(salophen)-HQ (Cat. 1)", "borylation"),
    ("fig2a_bq", "m_1g", "m_3g", "m_bq", "Co(salophen) and BQ", "borylation"),
    ("fig2a_hq", "m_1g", "m_3g", "m_hq", "Co(salophen) and HQ", "borylation"),
    ("fig2b_cat1", "m_1s", "m_6", "m_cat1", "Co(salophen)-HQ (Cat. 1)", "elimination"),
    ("fig2b_bq", "m_1s", "m_6", "m_bq", "Co(salophen) and BQ", "elimination"),
    ("fig2b_hq", "m_1s", "m_6", "m_hq", "Co(salophen) and HQ", "elimination"),
]
for rid, substrate, product, etm_mid, curve, kind in profile_specs:
    panel = "a" if rid.startswith("fig2a") else "b"
    r = base_reaction(rid, [EV[4]], "specific_experiment", ["mechanistic_study"], {"page": 4, "region": f"Figure 2{panel}, {curve} curve"}, "reaction profile with different ETMs", "kinetic reaction-profile experiment")
    participant(r, "substrate", substrate, "reactant", [EV[4]])
    participant(r, "product", product, "product", [EV[4]])
    if "bq" in rid or "hq" in rid:
        participant(r, "co_salophen", "m_co_salophen", "catalyst", [EV[4]])
    participant(r, "etm", etm_mid, "additive" if etm_mid in ("m_bq", "m_hq") else "catalyst", [EV[4]])
    add_reaction_step(r, [EV[4]], [p["participant_id"] for p in r["participants"]], operation="time-course sampling for the plotted reaction profile; detailed conditions referred to Supporting Information")
    o = st.new_outcome("reaction_profile", "other", [product], [EV[4]])
    o["value"] = unreadable(f"Figure 2{panel} plots yield (%) versus time (h) for {curve}", [EV[4]])
    o["scope"] = direct({"x_axis": "Time (h)", "y_axis": f"Yield of {product.replace('m_', '')} (%)", "curve": curve}, {"x_axis": "Time (h)", "y_axis": f"Yield of {product.replace('m_', '')} (%)", "curve": curve}, [EV[4]])
    r["outcomes"].append(o)
    set_procedure(r, "Figure 2 caption: for details, please see the supporting information.", [EV[4]])
    if kind == "borylation":
        set_classification(r, [EV[4]], [("C", "C", 0, 1, "formed", "corresponding 1g to 3g ring closure"), ("C", "B", 0, 1, "formed", "corresponding 1g to 3g borylation")], "Kinetic curve identity is clear, but individual point values and complete mapping are unavailable.")
    else:
        set_classification(r, [EV[4]], [("C", "C", 0, 1, "formed", "corresponding 1s to 6 ring closure"), ("C", "C", 1, 2, "order_increased", "corresponding beta-H-elimination alkene formation")], "Kinetic curve identity is clear, but individual point values and complete mapping are unavailable.")
    r["review"]["unresolved"] = ["Exact plotted datapoints, replicate/error information, analytical method, charges and detailed reaction conditions require the unavailable SI."]


# Proposed mechanism is not converted into an experimental proof.
r = base_reaction("scheme6_proposed_mechanism", [EV[4], EV[5]], "proposed_mechanism", ["other"], {"pages": [4, 5], "region": "mechanistic discussion and Scheme 6"}, "proposed catalytic cycle", "proposed PdII/Pd0 oxidative carbocyclization cycle")
participant(r, "generic_substrate", "m_generic_bisallene_1", "reactant", [EV[4], EV[5]])
participant(r, "generic_product", "m_generic_cycle_3", "product", [EV[4], EV[5]])
r["notes"].append(direct("Int-I coordination; allenic C(sp3)-H cleavage to Int-II; distal allene coordination to Int-III; insertion to seven-membered pi-allyl Pd Int-IV; trapping; aerobic Pd reoxidation.", "Author-proposed sequence, not a directly observed preparative experiment.", [EV[4], EV[5]]))
set_no_transform_classification(r, [EV[4], EV[5]], "Proposed elementary cycle is non-preparative and is not treated as an observed net experiment.")
r["review"]["unresolved"] = ["Intermediates Int-I–Int-IV are proposed, not established as isolated or directly observed species on supplied pages."]


# Footnote 21 reports a current-method DMBQ comparison without a numerical result.
r = base_reaction("dmbq_less_efficient_observation", [EV[6]], "specific_experiment", ["optimization"], {"page": 6, "region": "reference footnote 21"}, "DMBQ comparison in the borylative carbocyclization", "optimization observation")
participant(r, "substrate", "m_1a", "reactant", [EV[6]])
participant(r, "product", "m_3a", "product", [EV[6]])
participant(r, "dmbq", "m_dmbq", "additive", [EV[6]])
add_reaction_step(r, [EV[6]], ["substrate", "product", "dmbq"], operation="DMBQ tested in this reaction; exact conditions not stated in footnote")
o = st.new_outcome("relative_efficiency", "other", ["m_3a"], [EV[6]])
o["value"] = direct("DMBQ was less efficient in this reaction", "less efficient than the reported preferred mediator; numerical value unreported", [EV[6]])
o["scope"] = direct("qualitative comparison", "qualitative comparison", [EV[6]])
r["outcomes"].append(o)
set_procedure(r, "Footnote 21 provides only the qualitative statement and literature context.", [EV[6]])
set_classification(r, [EV[6]], [("C", "C", 0, 1, "formed", "same stated borylative carbocyclization family"), ("C", "B", 0, 1, "formed", "same stated borylative carbocyclization family")], "Specific yield and full experimental conditions are unreported.")
r["review"]["unresolved"] = ["Numerical yield, DMBQ loading, other conditions and analytical basis are not reported in the supplied footnote."]


# Non-experiment context and explicit exclusions.
bundle["source_context"] = [
    {"context_id": f"{DOC}:metadata", "document_id": DOC, "kind": "bibliographic_metadata", "title": "Efficient Palladium-Catalyzed Aerobic Oxidative Carbocyclization to Seven-Membered Heterocycles", "authors": ["Jie Liu", "Jan-E. Bäckvall"], "doi": "10.1002/chem.202004265", "evidence_ids": [EV[1]]},
    {"context_id": f"{DOC}:scheme1_prior_strategies", "document_id": DOC, "kind": "prior_work_context", "note": "Scheme 1 summarizes literature strategies and is not counted as a present-paper experiment.", "evidence_ids": [EV[1]]},
    {"context_id": f"{DOC}:figure2_unreadable_points", "document_id": DOC, "kind": "kinetic_plot_context", "note": "Figure 2 contains real time-course curves; six experiment records retain curve identities, while exact point values remain unreadable/unreported without SI.", "evidence_ids": [EV[4]]},
    {"context_id": f"{DOC}:pages5_6_bibliography", "document_id": DOC, "kind": "bibliography_and_summary", "note": "Most of pages 5–6 are summary/references. Cited literature is not imported as present-paper experiments.", "evidence_ids": [EV[5], EV[6]]},
    {"context_id": f"{DOC}:no_characterization_tables", "document_id": DOC, "kind": "coverage_exclusion", "note": "No standalone MS/NMR peak table or fragment-abundance table occurs on supplied pages 1–6; therefore none was misclassified as a synthetic reaction.", "evidence_ids": [EV[1], EV[2], EV[3], EV[4], EV[5], EV[6]]},
]


# Raw inventory, required page coverage, and task states.
page_records = defaultdict(list)
for r in bundle["reactions"]:
    for eid in r["evidence_ids"]:
        page = int(eid.split(":page", 1)[1].split(":", 1)[0])
        page_records[page].append(r["reaction_id"])

regions = {
    1: [
        {"region": "title, abstract and main text", "finding": "Bibliographic metadata and present-work overview; no row-level preparative experiment.", "evidence_ids": [EV[1]]},
        {"region": "Scheme 1", "finding": "Prior-work strategy overview, excluded from present-paper experiment count.", "evidence_ids": [EV[1]]},
    ],
    2: [
        {"region": "Scheme 2", "finding": "Prior-work generic concept (2a) and current-work generic transformation (2b); only 2b recorded as a generic scheme.", "evidence_ids": [EV[2]]},
        {"region": "Scheme 3 and catalyst-yield table", "finding": "Nine catalyst-preparation rows, Cat.1–Cat.9, recorded separately.", "evidence_ids": [EV[2]]},
        {"region": "Table 1 including footnotes", "finding": "Six optimization entries; entry 1 includes 93% recovered 1a and entries 5–6 include <1% recovered 1a.", "evidence_ids": [EV[2]]},
    ],
    3: [
        {"region": "Figure 1 and caption", "finding": "Nine catalyst-screen bars; Cat.1 duplicates Table 1 entry 5 and is merged by evidence, Cat.2–Cat.9 are eight records.", "evidence_ids": [EV[3]]},
        {"region": "Scheme 4 and footnote", "finding": "Fifteen isolated-yield scope arrows 1a–1o to 3a–3o; 3f uses 15 mol% Cat.1. The caption is also the explicitly marked inherited-condition source for the immediately continuing 1s limitation and the page-4 control described as using standard conditions.", "evidence_ids": [EV[3]]},
        {"region": "scope discussion at page bottom", "finding": "N-containing analogue attempt continues onto page 4 and is retained as one unsuccessful control with cross-page evidence.", "evidence_ids": [EV[3]]},
    ],
    4: [
        {"region": "top-left continuation and Eq. (1)", "finding": "NTs analogue unsuccessful; 1s borylation not observed; 1a' control gave 0% product.", "evidence_ids": [EV[4]]},
        {"region": "Scheme 5", "finding": "Three separate arylative, carbonylative and beta-H-eliminative carbocyclization experiments.", "evidence_ids": [EV[4]]},
        {"region": "Figure 2", "finding": "Six kinetic reaction-profile curves retained as six mechanistic-study experiments; exact datapoints unavailable without SI.", "evidence_ids": [EV[4]]},
        {"region": "mechanistic discussion", "finding": "Author-proposed cycle discussion linked to Scheme 6; not treated as experimental proof.", "evidence_ids": [EV[4]]},
    ],
    5: [
        {"region": "Scheme 6", "finding": "Proposed catalytic cycle, recorded as proposed_mechanism rather than a preparative experiment.", "evidence_ids": [EV[5]]},
        {"region": "summary, acknowledgements and references 1–20", "finding": "Narrative/bibliographic context only; cited work not imported as current experiments.", "evidence_ids": [EV[5]]},
    ],
    6: [
        {"region": "reference footnote 21", "finding": "Qualitative present-reaction observation that DMBQ was less efficient; numerical yield and conditions unreported.", "evidence_ids": [EV[6]]},
        {"region": "references 20–23 and manuscript dates", "finding": "Bibliographic metadata only; no additional present-paper reaction arrows or characterization tables.", "evidence_ids": [EV[6]]},
    ],
}
unresolved_by_page = {
    1: ["No experiment-level conditions on this page; Scheme 1 is prior-work context."],
    2: ["Scheme 3 detailed procedures and yield bases require SI.", "Table 1 does not state the measurement basis for starting-material recovery."],
    3: ["Complete scope structures and stereochemical assignments are not safely encodable from the compact main-text drawings alone."],
    4: ["Figure 2 exact datapoints and detailed conditions are deferred to unavailable SI.", "Scheme 5 yield measurement basis is not stated."],
    5: ["Proposed intermediates are not direct experimental observations."],
    6: ["DMBQ comparison is qualitative only; no numerical result or loading is supplied."],
}
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": p, "state": "attempted", "regions": regions[p], "record_ids": list(dict.fromkeys(page_records[p])), "unresolved": unresolved_by_page[p]}
        for p in range(1, 7)
    ],
}

# Coverage must point to the real rendered page evidence as well as searchable
# page text because every required PDF page was visually inspected.
for page_entry in coverage["pages"]:
    for region in page_entry["regions"]:
        region["evidence_ids"] = expand_evidence(region["evidence_ids"])

for context in bundle["source_context"]:
    context["evidence_ids"] = expand_evidence(context["evidence_ids"])

raw_inventory = {
    "document_id": DOC,
    "source_pages_viewed": [1, 2, 3, 4, 5, 6],
    "source_images_viewed": [f"runs/current/assets/{DOC}/page{p}.jpg" for p in range(1, 7)],
    "reaction_count": len(bundle["reactions"]),
    "molecule_count": len(bundle["molecules"]),
    "records": [{"reaction_id": r["reaction_id"], "record_type": r["record_type"], "evidence_ids": r["evidence_ids"]} for r in bundle["reactions"]],
    "stereo_diagnostics_file": "stereo_diagnostics.json",
    "limitations": ["Raw model inventory is not a source-review certification.", "No SI file was supplied or associated."],
}

st.write(RUN / "submitted.json", bundle)
st.write(RUN / "raw_model/model_extraction.json", raw_inventory)
st.write(RUN / "raw_model/stereo_diagnostics.json", stereo_diagnostics)
st.write(RUN / "coverage/pages_0001_0006.json", coverage)

tasks_path = RUN / "tasks.json"
tasks = json.loads(tasks_path.read_text())
for task in tasks:
    task["state"] = "attempted"
    page = int(task["evidence_id"].split(":page", 1)[1].split(":", 1)[0])
    task["visual_evidence_id"] = IMG[page]
    task["reason"] = "The rendered original PDF page and extracted page text were inspected; experiment mapping and unresolved decisions are recorded in coverage/pages_0001_0006.json and runs/workflow.md."
st.write(tasks_path, tasks)

print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]), "stereo_checked_smiles": len(stereo_diagnostics)}, ensure_ascii=False))
