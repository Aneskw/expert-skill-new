#!/usr/bin/env python3
"""Build page-bounded v3.1 machine records for doc_f020a06cbd6cce19.

The source images, not this script, are the authority.  All structures are
machine candidates and every imported record remains needs_review.
"""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_f020a06cbd6cce19"
RUN = ROOT / "runs/current"
DATASET = json.loads((RUN / "dataset.json").read_text())
BUNDLE = json.loads((RUN / "submitted.json").read_text())
BUNDLE["documents"] = copy.deepcopy(DATASET["documents"])
BUNDLE["molecules"] = []
BUNDLE["reactions"] = []
BUNDLE["source_context"] = []

P = {n: f"{DOC}:page{n}:image" for n in range(1, 6)}
PT = {n: f"{DOC}:page{n}:text" for n in range(1, 6)}
MOLECULES: dict[str, dict] = {}
REACTIONS: list[dict] = []
STEREO_CHECKS: dict[str, dict] = {}


def F(raw=None, norm=None, ev=(), method="original_extraction", derivation=None,
      status=None):
    return st.fact(raw, norm, ev, method=method, derivation=derivation,
                   status=status)


def missing(method="original_extraction"):
    return st.fact(method=method)


def molecule(mid, label, ev, smiles=None, name=None, geometry=None,
             representation=None, note=None):
    if mid in MOLECULES:
        return mid
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, ev)
    m["names"]["original"] = F(name or label, name or label, ev)
    m["names"]["english"] = F(name or label, name or label, ev)
    m["names"]["abbreviations"] = missing()
    m["names"]["aliases"] = missing()
    m["cas_numbers"] = missing()
    m["external_properties"] = []
    m["identity_note"] = F(
        note or "Machine identity candidate from the labelled source drawing; expert structure review remains required.",
        note or "Machine identity candidate from the labelled source drawing; expert structure review remains required.",
        ev,
    )
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = F(
            {"source": "labelled structure drawing", "compound_label": label},
            smiles, ev, "structure_derivation",
            "Connectivity and any slash/backslash geometry were transcribed from the labelled source drawing; parsing does not verify the drawing assignment.",
        )
        diag = st.stereo_inspect(smiles)
        STEREO_CHECKS[mid] = diag
        m["stereochemistry"]["atom_centers"] = F(
            "No tetrahedral R/S label was reported in the inspected main-text drawing.",
            "unspecified", ev)
        m["stereochemistry"]["axial"] = F(
            "No axial assignment was reported in the inspected main-text drawing.",
            "unspecified", ev)
        m["stereochemistry"]["mixture"] = missing()
        if geometry in ("E", "Z"):
            m["stereochemistry"]["double_bonds"] = F(
                f"{label} is drawn/reported as the {geometry}-isomer.",
                [{"configuration": geometry, "location": "product alkene"}], ev)
            m["stereochemistry"]["assignment"] = F(
                f"Relative alkene geometry reported as {geometry}.",
                {"alkene": geometry, "assignment_type": "relative"}, ev)
            m["stereochemistry"]["geometry_status"] = F(
                geometry, "specified", ev)
            has_other_unassigned = any(p["specified"] == "Unspecified" for p in diag["potential_stereo"])
            representation_status = "partial" if has_other_unassigned else "encoded"
            m["stereochemistry"]["representation_status"] = F(
                "product-alkene stereobond encoded; other potential stereochemistry remains unassigned" if has_other_unassigned else "stereobond candidate encoded in SMILES",
                representation_status, ev,
                "structure_derivation",
                "Directional SMILES encodes the source-reported alkene geometry; any additional toolkit-flagged potential element remains unassigned, and toolkit diagnostics are graph checks only.")
            m["stereochemistry"]["isomeric_smiles"] = F(
                {"source_geometry": geometry}, smiles, ev, "structure_derivation",
                "Isomeric SMILES generated from the source-labelled E/Z drawing.")
        else:
            m["stereochemistry"]["double_bonds"] = F(
                "No E/Z assignment was reported for this identity in the inspected source region.",
                "unspecified", ev)
            m["stereochemistry"]["assignment"] = missing()
            m["stereochemistry"]["geometry_status"] = F(
                "not assigned in source", "unresolved", ev)
            m["stereochemistry"]["representation_status"] = F(
                "stereochemical representation not asserted", "unresolved", ev)
            m["stereochemistry"]["isomeric_smiles"] = missing()
        m["stereochemistry"]["chirality_status"] = F(
            "No center/axis configuration was assigned in the inspected source region.",
            "unspecified", ev)
        m["stereochemistry"]["identity_key"] = F(
            None, diag["stereo_graph_key"], ev, "structure_derivation",
            "Computed by standard stereo from the candidate isomeric graph; not a sample identity.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing()
        m["generic_template"] = None
        m["stereochemistry"]["atom_centers"] = F(
            "Atom-centered stereochemistry could not be safely encoded from the inspected main-text miniature.",
            "unresolved", ev)
        m["stereochemistry"]["double_bonds"] = F(
            f"Source geometry statement: {geometry}" if geometry else "No complete geometry assignment safely encoded.",
            ([{"configuration": geometry}] if geometry in ("E", "Z") else "unresolved"), ev)
        m["stereochemistry"]["axial"] = F(
            "Axial stereochemistry not assigned in source.", "unspecified", ev)
        m["stereochemistry"]["mixture"] = missing()
        m["stereochemistry"]["assignment"] = F(
            f"Source reports {geometry}-geometry; structure encoding unresolved." if geometry else None,
            ({"alkene": geometry} if geometry else None), ev) if geometry else missing()
        m["stereochemistry"]["chirality_status"] = F(
            "No center/axis assignment reported.", "unspecified", ev)
        m["stereochemistry"]["geometry_status"] = F(
            geometry or "unresolved", "specified" if geometry else "unresolved", ev)
        m["stereochemistry"]["representation_status"] = F(
            "Source geometry is not encoded" if geometry else "Structure unresolved",
            "not_encoded" if geometry else "unresolved", ev)
        m["stereochemistry"]["isomeric_smiles"] = missing()
        m["stereochemistry"]["identity_key"] = missing()
    MOLECULES[mid] = m
    return mid


def template_molecule(mid, label, ev, fixed, variables, note):
    if mid in MOLECULES:
        return mid
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, ev)
    m["names"]["original"] = F(label, label, ev)
    m["names"]["english"] = F(label, label, ev)
    m["names"]["abbreviations"] = missing()
    m["names"]["aliases"] = missing()
    m["cas_numbers"] = missing()
    m["representation_type"] = "generic_template"
    m["smiles"] = missing()
    m["generic_template"] = {
        "fixed_structure": F(fixed, fixed, ev),
        "variables": [{
            "label": v,
            "attachment_points": F(f"Source-labelled attachment point for {v}", f"source attachment for {v}", ev),
            "definition": F(f"Variable {v} in source drawing", f"variable {v}", ev),
            "allowed_values": missing(),
            "constraints": missing(),
        } for v in variables],
        "instance_mapping": [],
        "reason": F(note, note, ev),
    }
    for key, val in (("atom_centers", "unresolved"), ("double_bonds", "unresolved"),
                     ("axial", "unresolved"), ("mixture", "unresolved"),
                     ("assignment", "unresolved")):
        m["stereochemistry"][key] = F("Generic source drawing", val, ev)
    m["stereochemistry"]["chirality_status"] = F("Generic source drawing", "unresolved", ev)
    m["stereochemistry"]["geometry_status"] = F("Generic source drawing", "unresolved", ev)
    m["stereochemistry"]["representation_status"] = F("Generic source drawing", "not_encoded", ev)
    m["stereochemistry"]["isomeric_smiles"] = missing()
    m["stereochemistry"]["identity_key"] = missing()
    m["identity_note"] = F(note, note, ev)
    m["external_properties"] = []
    MOLECULES[mid] = m
    return mid


def participant(pid, mid, role, ev, stage="reaction", eq=None, eq_ref=None,
                volume=None, mol_percent=None, product_role=None):
    p = st.new_participant(pid, mid, role, ev, stage)
    if product_role:
        p["product_role"] = F(product_role, product_role, ev)
    if eq is not None:
        p["charges"]["equivalents"] = F(
            f"{eq:g} equiv.", {"number": eq, "unit": "equiv"}, ev)
        p["equivalence_reference"] = F(
            f"relative to {eq_ref}", eq_ref, ev)
    if volume is not None:
        raw, norm = volume
        p["charges"]["volume"] = F(raw, norm, ev)
    if mol_percent is not None:
        p["charges"]["mol_percent"] = F(
            f"{mol_percent:g} mol%", {"number": mol_percent, "unit": "mol%"}, ev)
        p["equivalence_reference"] = F(f"relative to {eq_ref}", eq_ref, ev)
    return p


def base_reaction(rid, ev, record_type="specific_experiment", categories=None,
                  location=None, name=None, transformation=None):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = F(categories or [], categories or [], ev)
    r["location"] = F(location or rid, {"description": location or rid}, ev)
    r["reaction_name"] = F(name, name, ev) if name else missing()
    r["transformation_type"] = F(
        None, transformation, ev, "inference",
        "Transformation label inferred from the source reactant/product drawings; not an author-named reaction." if transformation else None) if transformation else missing()
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    return r


def set_class(r, bond, ev):
    element = {"thioether": "S", "ether": "O", "azide": "N"}.get(bond)
    if not element:
        r["classification"]["coverage"] = F(
            None, "unresolved", ev, "inference",
            "Net atom mapping is unresolved for this transformation; no bond labels asserted.")
        r["classification"]["limitations"] = [
            "No complete atom mapping was asserted from the main-text miniature.",
            "Classification remains unresolved and requires expert review.",
        ]
        return
    changes = [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": element}],
         "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}],
         "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
        {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
         "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
    ]
    if bond in ("thioether", "ether"):
        changes.append({"change_id": "b4", "atoms": [{"atom_map": None, "element": element}, {"atom_map": None, "element": "H"}],
                        "before_order": 1, "after_order": 0, "change_type": "hydrogen_loss", "hydrogen_count": 1})
    labels = []
    for c in changes:
        if c["change_type"] in ("formed", "hydrogen_gain", "order_increased"):
            lab = st.bond_label(c["atoms"][0]["element"], c["atoms"][1]["element"], c["after_order"])
            lab["supporting_change_ids"] = [c["change_id"]]
            labels.append(lab)
    deriv = "Partial net bond changes inferred from the generic source drawings; atom identities and proton source are not completely mapped."
    r["classification"]["bond_changes"] = [F(None, c, ev, "structure_derivation", deriv) for c in changes]
    r["classification"]["labels"] = F(None, labels, ev, "structure_derivation", deriv)
    r["classification"]["coverage"] = F(None, "partial", ev, "inference", deriv)
    r["classification"]["limitations"] = [
        "Partial graph comparison only; no complete mapped reaction SMILES asserted.",
        "C-H is a net gain and the proton source is not atom-mapped.",
        "The decreased C-C bond order is retained but is not classified as a newly formed C-C sigma bond.",
    ]


def add_outcome(r, oid, kind, targets, raw, norm, ev, basis_raw=None,
                basis_norm=None, scope=None, product_role="principal_product"):
    o = st.new_outcome(oid, kind, targets, ev)
    o["product_role"] = F(product_role, product_role, ev)
    o["value"] = F(raw, norm, ev)
    o["yield_type"] = missing()
    if basis_raw is not None:
        o["measurement_basis"] = F(basis_raw, basis_norm or basis_raw, ev)
    else:
        o["measurement_basis"] = missing()
    o["scope"] = F(scope, scope, ev) if scope else missing()
    r["outcomes"].append(o)


def reaction_step(sid, ev, pids, temp, time, atmosphere, operation):
    s = st.new_step(sid, ev, "reaction")
    s["participant_ids"] = list(pids)
    s["operation"] = F(operation, operation, ev)
    if isinstance(temp, (int, float)):
        s["temperature"] = F(f"{temp:g} °C", {"number": temp, "unit": "°C"}, ev)
    elif temp:
        s["temperature"] = F(temp, temp, ev)
    if isinstance(time, tuple):
        lo, hi = time
        s["time"] = F(f"{lo:g}–{hi:g} h", {"minimum": lo, "maximum": hi, "unit": "h", "raw_text": f"{lo:g}–{hi:g} h"}, ev)
    elif isinstance(time, (int, float)):
        s["time"] = F(f"{time:g} h", {"number": time, "unit": "h"}, ev)
    elif time:
        s["time"] = F(time, time, ev)
    if atmosphere:
        s["atmosphere"] = F(atmosphere, atmosphere, ev)
    return s


def add_workup(r, ev, prefix):
    mids = {
        "water": molecule("water", "water", [P[4]], "O", "water"),
        "etoac": molecule("ethyl_acetate", "ethyl acetate", [P[4]], "CCOC(C)=O", "ethyl acetate"),
        "brine": molecule("brine", "brine solution", [P[4]], None, "brine solution",
                          note="Named aqueous mixture; composition/concentration not specified, so no unique SMILES asserted."),
        "na2so4": molecule("sodium_sulfate", "anhydrous sodium sulfate", [P[4]], "[Na+].[Na+].[O-]S(=O)(=O)[O-]", "anhydrous sodium sulfate"),
        "silica": molecule("silica_gel", "silica gel", [P[4]], None, "silica gel",
                           note="Named heterogeneous purification material; no unique molecular SMILES asserted."),
    }
    wp = [
        participant(prefix+"_water", mids["water"], "workup_agent", ev, "workup"),
        participant(prefix+"_etoac", mids["etoac"], "workup_agent", ev, "workup",
                    volume=("3 × 20 mL", {"number": 20, "unit": "mL", "replicates": 3, "raw_text": "3 × 20 mL"})),
        participant(prefix+"_brine", mids["brine"], "workup_agent", ev, "workup"),
        participant(prefix+"_na2so4", mids["na2so4"], "workup_agent", ev, "workup"),
        participant(prefix+"_silica", mids["silica"], "purification_agent", ev, "purification"),
    ]
    r["participants"].extend(wp)
    ids = {x["molecule_id"]: x["participant_id"] for x in wp}
    s1 = st.new_step(prefix+"_quench_extract", ev, "workup")
    s1["participant_ids"] = [ids[mids["water"]], ids[mids["etoac"]]]
    s1["operation"] = F("quenched with water and extracted with ethyl acetate (3 × 20 mL)",
                         "quench with water; extract three times with 20 mL ethyl acetate", ev)
    s2 = st.new_step(prefix+"_wash_dry", ev, "workup")
    s2["participant_ids"] = [ids[mids["brine"]], ids[mids["na2so4"]]]
    s2["operation"] = F("organic part was washed with brine solution and dried with anhydrous sodium sulfate",
                         "wash with brine; dry over anhydrous sodium sulfate", ev)
    s3 = st.new_step(prefix+"_evaporate", ev, "workup")
    s3["participant_ids"] = []
    s3["operation"] = F("solvent was evaporated under reduced pressure", "evaporate under reduced pressure", ev)
    s4 = st.new_step(prefix+"_chromatography", ev, "purification")
    s4["participant_ids"] = [ids[mids["silica"]]]
    s4["operation"] = F("purified by column chromatography using silica gel in a specified solvent combination",
                         "silica-gel column chromatography; eluent not reported in main text", ev)
    r["steps"].extend([s1, s2, s3, s4])


# Shared source identities.
template_molecule("generic_allenic_ketone_1", "allenic ketone 1", [P[1], P[2], P[3]],
                  "R1-C(=O)-allene-R2", ["R1", "R2"],
                  "Generic R-group drawing in Schemes/Tables 2–5; not a concrete SMILES.")
template_molecule("generic_beta_vinyl_product", "E beta-vinyl alpha,beta-unsaturated ketone", [P[1], P[2], P[3]],
                  "R1-C(=O)-C=C(XR)-CH2R2", ["R1", "R2", "R", "X"],
                  "Generic E-product family with X=S/O or N3; variable groups are not expanded as one molecule.")
template_molecule("generic_benzodioxin_5", "benzodioxin 5", [P[2], P[3]],
                  "substituted benzodioxin ketone", ["R", "R1"],
                  "Generic benzodioxin family in Table 5; substitution variables prevent a concrete SMILES.")
template_molecule("generic_ortho_hydroxybenzaldehyde", "substituted ortho-hydroxybenzaldehyde", [P[2], P[3]],
                  "ortho-HO-C6H3(R1)-CHO", ["R1"],
                  "Generic source reactant for the benzodioxin cascade.")

common = {
    "li2co3": molecule("lithium_carbonate", "Li2CO3", [P[1], P[2], P[3], P[4]], "[Li+].[Li+].[O-]C([O-])=O", "lithium carbonate"),
    "thf": molecule("thf", "THF", [P[1], P[2], P[3], P[4]], "C1CCOC1", "tetrahydrofuran"),
    "dmf": molecule("dmf", "DMF", [P[2], P[3], P[4]], "CN(C)C=O", "N,N-dimethylformamide"),
    "oxone": molecule("oxone", "Oxone", [P[2], P[3], P[4]], None, "Oxone",
                      note="Commercial named oxidant mixture; composition is not uniquely specified in the article, so no concrete SMILES asserted."),
    "tmsn3": molecule("tms_azide", "TMS-N3", [P[1], P[3]], "C[Si](C)(C)N=[N+]=[N-]", "trimethylsilyl azide"),
    "toluene": molecule("toluene", "toluene", [P[4]], "Cc1ccccc1", "toluene"),
}


def generic_scope(rid, kind, ev, page_location, general_name):
    r = base_reaction(rid, ev, "generic_scheme", ["other"], page_location,
                      general_name, f"nucleophilic addition to an allenic ketone ({kind})")
    r["participants"] = [
        participant("substrate", "generic_allenic_ketone_1", "reactant", ev),
        participant("product", "generic_beta_vinyl_product", "product", ev, product_role="principal_product"),
    ]
    r["steps"] = [reaction_step("generic_reaction", ev, ["substrate", "product"], 35, (4, 7), "open-air", "stir under the conditions printed over the generic arrow")]
    set_class(r, kind, ev)
    r["review"]["unresolved"] = ["Generic R groups are not concrete molecular identities."]
    REACTIONS.append(r)
    return r["reaction_id"]


GEN_THIO = generic_scope("generic_thioether", "thioether", [P[1], P[2]], "Scheme 2 and Table 2 header", "regioselective (E) beta-vinyl thioether synthesis")
GEN_ETHER = generic_scope("generic_ether", "ether", [P[1], P[3]], "Scheme 2 and Table 3 header", "regioselective (E) beta-vinyl ether synthesis")
GEN_AZIDE = generic_scope("generic_azide", "azide", [P[1], P[3]], "Scheme 2 and Table 4 header", "regioselective (E) beta-vinyl azide synthesis")

# Table 5 generic arrow and proposed cascade are distinct from concrete rows.
gen_benzo = base_reaction("generic_benzodioxin", [P[2], P[3]], "generic_scheme", ["other"],
                          "Scheme 3 'This Work' and Table 5 header", "benzodioxin synthesis", "one-pot benzodioxin annulation")
gen_benzo["participants"] = [
    participant("allene", "generic_allenic_ketone_1", "reactant", [P[2], P[3]]),
    participant("aldehyde", "generic_ortho_hydroxybenzaldehyde", "reactant", [P[2], P[3]]),
    participant("product", "generic_benzodioxin_5", "product", [P[2], P[3]], product_role="principal_product"),
]
gen_benzo["steps"] = [reaction_step("annulation", [P[2], P[3]], ["allene", "aldehyde", "product"], 25, 12, "open-air", "one-pot annulation under Oxone/Li2CO3 in DMF")]
set_class(gen_benzo, None, [P[2], P[3]])
gen_benzo["review"]["unresolved"] = ["Generic R/R1 identities and complete cascade atom mapping are unresolved."]
REACTIONS.append(gen_benzo)
GEN_BENZO = gen_benzo["reaction_id"]

mech = base_reaction("proposed_benzodioxin_cascade", [P[3]], "proposed_mechanism", ["other"],
                     "Table 5 lower dashed pathway", "proposed cyclization sequence", "proposed addition/oxidation/cyclization cascade")
mech["participants"] = [
    participant("allene", "generic_allenic_ketone_1", "reactant", [P[3]]),
    participant("aldehyde", "generic_ortho_hydroxybenzaldehyde", "reactant", [P[3]]),
    participant("product", "generic_benzodioxin_5", "product", [P[3]], product_role="proposed_product"),
]
set_class(mech, None, [P[3]])
mech["review"]["unresolved"] = ["Dashed pathway is an author-proposed sequence, not a separately demonstrated mechanistic experiment."]
REACTIONS.append(mech)

# Scheme 1: one prior-work experiment produces four regio/geometric products.
lit_allene = template_molecule("scheme1_allenic_ketone", "Sugita allenic ketone", [P[1]], "R-C(=O)-allene-R1", ["R", "R1"], "Generic prior-work reactant in Scheme 1.")
lit_thiol = template_molecule("scheme1_thiol", "RSH", [P[1]], "R-SH", ["R"], "Generic thiol in cited Sugita work.")
for x in ("I", "II", "III", "IV"):
    template_molecule(f"scheme1_product_{x}", f"Sugita product {x}", [P[1]], f"drawn vinyl sulfide isomer {x}", ["R", "R1", "R2"], "One of four generic regio/geometric products drawn for cited work.")
lit = base_reaction("literature_sugita_scheme1", [P[1]], "literature_example", ["other"], "Scheme 1", "Sugita report", "thiol addition to allenic ketone")
lit["participants"] = [participant("allene", lit_allene, "reactant", [P[1]]), participant("thiol", lit_thiol, "reactant", [P[1]])]
for x in ("I", "II", "III", "IV"):
    lit["participants"].append(participant(f"product_{x}", f"scheme1_product_{x}", "product", [P[1]], product_role="reported_isomer"))
lit["steps"] = [reaction_step("prior_reaction", [P[1]], [p["participant_id"] for p in lit["participants"]], None, None, None, "RSH, Et3N in deuterated solvent in an NMR tube")]
set_class(lit, None, [P[1]])
lit["review"]["unresolved"] = ["Cited prior work: substrate identities, amounts, temperature, time and yields are not supplied in this article."]
REACTIONS.append(lit)

# Scheme 3 has seven visibly distinct cited strategy arrows.  Preserve each as a
# literature example without importing conditions absent from the main-text graphic.
scheme3_rows = [
    ("ref12", "vinyl acetate; mercuric acetate; sulfuric acid"),
    ("ref11", "K3PO4·3H2O; DME; 6 h; 100 °C; substituted DCM"),
    ("ref10", "KOH; Cu(OAc)2; L; NaHCO3; DMF; 120 °C"),
    ("ref9", "iodine; morpholine"),
    ("ref16", "allenic ketone; Oxone; Li2CO3; DMF; open-air; r.t.; 12 h"),
    ("refs14_15", "acetone/DMAP/SOCl2 or ketones/TFA/trifluoroacetic anhydride"),
    ("ref13", "aldehydes; DABCO"),
]
for code, cond in scheme3_rows:
    precursor = template_molecule(f"scheme3_precursor_{code}", f"Scheme 3 precursor {code}", [P[2]],
                                  "substituted aromatic precursor as drawn", ["R"], "Cited generic precursor; miniature does not establish a concrete identity.")
    r = base_reaction(f"literature_scheme3_{code}", [P[2]], "literature_example", ["other"],
                      f"Scheme 3, {code}", f"cited benzodioxin strategy {code}", "benzodioxin synthesis")
    r["participants"] = [participant("precursor", precursor, "reactant", [P[2]]), participant("product", "generic_benzodioxin_5", "product", [P[2]], product_role="reported_product_family")]
    r["steps"] = [reaction_step("literature_arrow", [P[2]], ["precursor", "product"], None, None, None, cond)]
    set_class(r, None, [P[2]])
    r["review"]["unresolved"] = ["Cited literature arrow only; amounts, yield basis and full procedure are not supplied in this article."]
    REACTIONS.append(r)


# Concrete molecule builders used by Tables 1–4.
ACYL = {
    "benzyl": "c1ccccc1CC(=O)",
    "phenethyl": "c1ccccc1CCC(=O)",
    "phenyl": "O=C(c1ccccc1)",
    "p_tolyl": "O=C(c1ccc(C)cc1)",
    "p_anisyl": "O=C(c1ccc(OC)cc1)",
    "o_anisyl": "O=C(c1ccccc1OC)",
    "o_chlorophenyl": "O=C(c1ccccc1Cl)",
    "p_fluorophenyl": "O=C(c1ccc(F)cc1)",
    "p_bromophenyl": "O=C(c1ccc(Br)cc1)",
    "dichlorophenyl": "O=C(c1c(Cl)cccc1Cl)",
    "o_tolyl": "O=C(c1ccccc1C)",
    "m_bromophenyl": "O=C(c1cccc(Br)c1)",
    "naphthyl": "O=C(c1cccc2ccccc12)",
    "thienyl": "O=C(c1ccsc1)",
    "aliphatic": "CCCCCCC(=O)",
    "allyloxyphenyl": "O=C(c1ccccc1OCC=C)",
    "benzyloxyphenyl": "O=C(c1ccccc1OCc2ccccc2)",
    "methoxy_methylphenyl": "O=C(c1cc(C)ccc1OC)",
}


def substrate(key, ev, r2=None):
    suffix = "_r2ph" if r2 == "Ph" else ""
    mid = f"allene_{key}{suffix}"
    if mid in MOLECULES:
        return mid
    smi = ACYL[key] + "C=C=C" + ("c1ccccc1" if r2 == "Ph" else "")
    note = "Candidate allenic ketone reconstructed from the product R1/R2 drawing; no axial configuration was assigned by the source."
    return molecule(mid, f"allenic ketone precursor ({key}{', R2=Ph' if r2 else ''})", ev, smi, note=note)


def e_product(mid, label, acyl, xfrag, ev, r2=None, geometry="E"):
    slash = "/" if geometry == "E" else "\\"
    tail = "Cc1ccccc1" if r2 == "Ph" else "C"
    smi = ACYL[acyl] + "/C=C(" + slash + xfrag + ")" + tail
    return molecule(mid, label, ev, smi, geometry=geometry)


def simple(mid, label, ev, smiles=None, name=None, note=None):
    return molecule(mid, label, ev, smiles, name, note=note)


BASE_SMILES = {
    "K2CO3": "[K+].[K+].[O-]C([O-])=O", "Na2CO3": "[Na+].[Na+].[O-]C([O-])=O",
    "Cs2CO3": "[Cs+].[Cs+].[O-]C([O-])=O", "Li2CO3": "[Li+].[Li+].[O-]C([O-])=O",
    "Et3N": "CCN(CC)CC", "Pyridine": "n1ccccc1",
}
SOLV_SMILES = {"CH3CN": "CC#N", "THF": "C1CCOC1", "1,4-Dioxane": "O1CCOCC1",
               "Benzene": "c1ccccc1", "CHCl3": "ClC(Cl)Cl", "DCM": "ClCCl",
               "DCE": "ClCCCl", "CDCl3": "[2H]C(Cl)(Cl)Cl"}
for name, smi in BASE_SMILES.items():
    simple("base_" + name.lower().replace(",", "").replace(" ", "_"), name, [P[2]], smi, name)
for name, smi in SOLV_SMILES.items():
    simple("solvent_" + name.lower().replace(",", "").replace(" ", "_").replace("-", "_"), name, [P[2]], smi, name)

sub_1f = substrate("benzyl", [P[1], P[2]])
thiocresol = simple("p_thiocresol", "p-thiocresol 2", [P[1], P[2]], "Cc1ccc(S)cc1", "4-methylbenzenethiol")
prod_2f_e = e_product("product_2f_E", "2f (E-isomer)", "benzyl", "Sc1ccc(C)cc1", [P[1], P[2], P[3]], geometry="E")
prod_2f_z = e_product("product_2f_Z", "2f′ (Z-isomer)", "benzyl", "Sc1ccc(C)cc1", [P[2]], geometry="Z")

opt_rows = [
    (1,"CH3CN","K2CO3","r.t.",60,"60:40"),(2,"THF","K2CO3","r.t.",70,"65:35"),
    (3,"THF","Na2CO3","r.t.",75,"65:35"),(4,"THF","Na2CO3",40,80,"85:15"),
    (5,"THF","Cs2CO3",40,70,"85:15"),(6,"THF","Li2CO3",40,85,"90:10"),
    (7,"THF","Li2CO3",35,96,"100:0"),(8,"THF","Li2CO3","r.t.",96,"55:45"),
    (9,"1,4-Dioxane","Li2CO3",35,85,"100:0"),(10,"Benzene","Li2CO3",35,70,"65:35"),
    (11,"CHCl3","Li2CO3",35,78,"100:0"),(12,"DCM","Li2CO3",35,75,"100:0"),
    (13,"DCE","Li2CO3",35,80,"100:0"),(14,"THF","Et3N",35,75,"83:17"),
    (15,"THF","Pyridine",35,70,"65:35"),(16,"CDCl3","Et3N",35,72,"80:20"),
    (17,"CDCl3","Li2CO3",35,80,"100:0"),
]
for n, solv, base, temp, yld, ratio in opt_rows:
    ev=[P[2]]
    r=base_reaction(f"table1_entry_{n:02d}",ev,"specific_experiment",["optimization"],f"Table 1 entry {n}","optimization for formation of 2f and 2f′","thiol addition to allenic ketone")
    base_mid="base_"+base.lower().replace(",","").replace(" ","_")
    solv_mid="solvent_"+solv.lower().replace(",","").replace(" ","_").replace("-","_")
    r["participants"]=[
        participant("substrate",sub_1f,"reactant",ev,eq=1.0,eq_ref=sub_1f),
        participant("thiol",thiocresol,"reactant",ev,eq=1.1,eq_ref=sub_1f),
        participant("base",base_mid,"reagent",ev,eq=1.1,eq_ref=sub_1f),
        participant("solvent",solv_mid,"solvent",ev,volume=("5 mL",{"number":5,"unit":"mL"})),
        participant("product_E",prod_2f_e,"product",ev,product_role="E-isomer"),
        participant("product_Z",prod_2f_z,"product",ev,product_role="Z-isomer"),
    ]
    r["solvent_conditions"]["composition"]=F(solv,solv,ev)
    r["solvent_conditions"]["total_volume"]=F("5 mL",{"number":5,"unit":"mL"},ev)
    r["steps"]=[reaction_step("reaction",ev,[p["participant_id"] for p in r["participants"][:4]],temp,5,"air atmosphere","combine reagents in dry solvent and hold under Table 1 conditions")]
    add_outcome(r,"yield","yield",[prod_2f_e,prod_2f_z],f"{yld}%",{"number":yld,"unit":"%"},ev,scope="combined 2f/2f′ product yield")
    add_outcome(r,"ez_ratio","selectivity",[prod_2f_e,prod_2f_z],ratio,{"metric":"E/Z","order":[prod_2f_e,prod_2f_z],"ratio":ratio},ev,"E/Z is determined by 1H-NMR","1H NMR",scope="E/Z ratio in order 2f:2f′",product_role="isomer_distribution")
    set_class(r,"thioether",ev)
    r["review"]["unresolved"]=["Yield measurement basis is not stated in Table 1.","Absolute reaction scale and work-up are not reported for this optimization row."]
    REACTIONS.append(r)


THIO = [
 ("2a",96,"p_anisyl","SCCC",None,"1-propanethiol"),("2b",90,"p_tolyl","SCCO",None,"2-mercaptoethanol"),
 ("2c",94,"o_anisyl","SCc1ccccc1",None,"benzyl mercaptan"),("2d",88,"p_anisyl","SCCO",None,"2-mercaptoethanol"),
 ("2e",95,"allyloxyphenyl","SC1CCCCC1",None,"cyclohexanethiol"),("2f",96,"benzyl","Sc1ccc(C)cc1",None,"p-thiocresol"),
 ("2g",93,"benzyl","Sc1ccccn1",None,"pyridine-2-thiol"),("2h",93,"phenethyl","Sc1ccccc1",None,"thiophenol"),
 ("2i",90,"benzyl","Sc1ccc(OC)cc1",None,"4-methoxythiophenol"),("2j",90,"benzyl","Sc1ccc(Br)cc1",None,"4-bromothiophenol"),
 ("2k",70,"phenethyl","Sc1ccc(N)cc1",None,"4-aminothiophenol"),("2l",88,"benzyl","Sc1ccc(F)cc1",None,"4-fluorothiophenol"),
 ("2m",94,"p_tolyl","Sc1ccc(C)cc1",None,"p-thiocresol"),("2n",92,"p_tolyl","Sc1ccc(Br)cc1",None,"4-bromothiophenol"),
 ("2o",90,"dichlorophenyl","Sc1ccc(F)cc1",None,"4-fluorothiophenol"),("2p",88,"o_anisyl","Sc1ccc(OC)cc1",None,"4-methoxythiophenol"),
 ("2q",90,"dichlorophenyl","Sc1ccc(Cl)cc1",None,"4-chlorothiophenol"),("2r",92,"p_bromophenyl","Sc1ccc(C)cc1",None,"p-thiocresol"),
 ("2s",92,"o_anisyl","Sc1ccccn1",None,"pyridine-2-thiol"),("2t",89,"dichlorophenyl","Sc1ccc(Br)cc1",None,"4-bromothiophenol"),
 ("2u",88,"p_tolyl","Sc1ccc(C)cc1","Ph","p-thiocresol"),("2v",90,"aliphatic","Sc1ccc(C)cc1",None,"p-thiocresol"),
 ("2w",85,"aliphatic","SCCO",None,"2-mercaptoethanol"),("2x",86,"benzyl","Sc1ccccc1N",None,"2-aminothiophenol"),
]

def reagent_from_x(label, xfrag, ev, name, kind):
    mid=f"nucleophile_{label}"
    if kind=="thioether":
        smi = {"SCCC":"CCCS","SCCO":"OCCS","SCc1ccccc1":"c1ccccc1CS","SC1CCCCC1":"SC1CCCCC1"}.get(xfrag,xfrag)
    elif kind=="ether":
        smi = {"OCC":"CCO","OC(C)C":"CC(C)O","OC":"CO"}.get(xfrag,xfrag)
    else:
        smi=None
    return simple(mid,name,ev,smi,name)


def scope_record(label,yld,acyl,xfrag,r2,nuc_name,kind,general_id,page=2):
    ev=[P[page],P[4]]
    sub=substrate(acyl,[P[page]],r2)
    prod = prod_2f_e if label=="2f" else e_product(f"product_{label}_E",f"{label} (E-isomer)",acyl,xfrag,[P[page]],r2,"E")
    nuc=reagent_from_x(label,xfrag,[P[page]],nuc_name,kind)
    rid=f"table{2 if kind=='thioether' else 3}_{label}"
    r=base_reaction(rid,ev,"specific_experiment",["substrate_scope"],f"Table {2 if kind=='thioether' else 3}, {label}",f"(E) beta-vinyl {kind} synthesis",f"nucleophilic addition to allenic ketone ({kind})")
    r["general_scheme_ids"]=[general_id]
    psub=participant("substrate",sub,"reactant",ev,eq=1.0,eq_ref=sub)
    pnuc=participant("nucleophile",nuc,"reactant",ev,eq=1.1,eq_ref=sub)
    pbase=participant("base",common["li2co3"],"reagent",ev,eq=1.1,eq_ref=sub)
    psolv=participant("solvent",common["thf"],"solvent",ev,volume=("5 mL",{"number":5,"unit":"mL"}))
    pprod=participant("product",prod,"product",ev,product_role="principal E-isomer")
    r["participants"]=[psub,pnuc,pbase,psolv,pprod]
    r["solvent_conditions"]["composition"]=F("THF","THF",ev)
    r["solvent_conditions"]["total_volume"]=F("5 mL",{"number":5,"unit":"mL"},[P[4]],"inherited_general_condition")
    r["steps"]=[reaction_step("reaction",ev,["substrate","nucleophile","base","solvent"],35,(4,7),"open atmosphere","General Procedure (A): allenic ketone in THF, then Li2CO3 and corresponding nucleophile; hold at 35 °C")]
    add_workup(r,[P[4]],label)
    r["procedure"]["general_procedure_reference"]=F("General Procedure (A)","General Procedure (A)",[P[4]],"inherited_general_condition")
    r["procedure"]["addition_order"]=F("allenic ketone in THF followed by Li2CO3 and corresponding nucleophile","substrate in THF; add Li2CO3; add nucleophile",[P[4]],"inherited_general_condition")
    r["procedure"]["pretreatment"]=F("oven dried single neck round-bottomed flask (10 mL capacity)","oven-dried 10 mL single-neck round-bottom flask",[P[4]],"inherited_general_condition")
    add_outcome(r,"yield","yield",[prod],f"{yld}%",{"number":yld,"unit":"%"},[P[page]],scope=f"Table scope product {label}")
    add_outcome(r,"geometry","selectivity",[prod],"E-isomer",{"metric":"alkene geometry","major_isomer":"E","reported_value":"E-isomer"},[P[page]],scope=f"geometry of {label}",product_role="reported E-isomer")
    set_class(r,kind,[P[page]])
    r["review"]["unresolved"]=["Yield measurement basis is not stated in the inspected main-text table.","Absolute substrate amount and product-specific chromatography eluent are not reported in the supplied main text; SI is unassociated."]
    REACTIONS.append(r)

for row in THIO:
    scope_record(*row,"thioether",GEN_THIO,2)

ETHER = [
 ("3a",96,"p_tolyl","OCC",None,"ethanol"),("3b",90,"o_chlorophenyl","OCC",None,"ethanol"),
 ("3c",88,"o_chlorophenyl","OC(C)C",None,"isopropanol"),("3d",92,"o_chlorophenyl","OC",None,"methanol"),
 ("3e",94,"phenethyl","OC",None,"methanol"),("3f",92,"p_tolyl","Oc1cccc(C=O)c1",None,"hydroxybenzaldehyde isomer drawn for 3f"),
 ("3g",92,"p_tolyl","Oc1ccc([N+](=O)[O-])cc1",None,"4-nitrophenol"),("3h",94,"p_tolyl","Oc1ccccc1[N+](=O)[O-]",None,"2-nitrophenol"),
 ("3i",90,"p_tolyl","Oc1ccccn1",None,"2-hydroxypyridine"),("3j",94,"p_tolyl","Oc1c(OC)cccc1C=O",None,"methoxy-hydroxybenzaldehyde drawn for 3j"),
 ("3k",95,"p_fluorophenyl","Oc1ccccc1C=O",None,"2-hydroxybenzaldehyde"),("3l",91,"p_fluorophenyl","Oc1ccccc1C=O",None,"hydroxybenzaldehyde isomer drawn for 3l"),
 ("3m",92,"p_tolyl","Oc1ccccc1C=O",None,"2-hydroxybenzaldehyde"),("3n",90,"p_tolyl","Oc1ccc(C=O)cc1",None,"4-hydroxybenzaldehyde"),
 ("3o",86,"aliphatic","OC",None,"methanol"),
]
for row in ETHER:
    scope_record(*row,"ether",GEN_ETHER,3)

AZIDE = [
 ("4a",86,"o_tolyl"),("4b",88,"m_bromophenyl"),("4c",84,"p_fluorophenyl"),("4d",85,"methoxy_methylphenyl"),
 ("4e",83,"allyloxyphenyl"),("4f",82,"dichlorophenyl"),("4g",90,"naphthyl"),("4h",85,"allyloxyphenyl"),
 ("4i",88,"benzyloxyphenyl"),("4j",83,"thienyl"),("4k",75,"benzyl"),("4l",78,"phenethyl"),
]
for label,yld,acyl in AZIDE:
    ev=[P[3],P[4]]
    sub=substrate(acyl,[P[3]])
    prod=e_product(f"product_{label}_E",f"{label} (E-isomer)",acyl,"N=[N+]=[N-]",[P[3]],None,"E")
    r=base_reaction(f"table4_{label}",ev,"specific_experiment",["substrate_scope"],f"Table 4, {label}","(E) beta-vinyl azide synthesis","azide addition to allenic ketone")
    r["general_scheme_ids"]=[GEN_AZIDE]
    r["participants"]=[participant("substrate",sub,"reactant",ev,eq=1.0,eq_ref=sub),participant("tmsn3",common["tmsn3"],"reagent",ev,eq=1.1,eq_ref=sub),participant("solvent",common["thf"],"solvent",ev,volume=("5 mL",{"number":5,"unit":"mL"})),participant("product",prod,"product",ev,product_role="principal E-isomer")]
    r["solvent_conditions"]["composition"]=F("THF","THF",ev)
    r["solvent_conditions"]["total_volume"]=F("5 mL",{"number":5,"unit":"mL"},[P[4]],"inherited_general_condition")
    r["steps"]=[reaction_step("reaction",ev,["substrate","tmsn3","solvent"],35,(4,7),"open-air","General Procedure (A), except no external base for azide; TMS-N3 in THF")]
    add_workup(r,[P[4]],label)
    r["procedure"]["general_procedure_reference"]=F("General Procedure (A), Li2CO3 excepted for azide","General Procedure (A), no Li2CO3 for azide",[P[3],P[4]],"inherited_general_condition")
    add_outcome(r,"yield","yield",[prod],f"{yld}%",{"number":yld,"unit":"%"},[P[3]],scope=f"Table scope product {label}")
    add_outcome(r,"geometry","selectivity",[prod],"E-isomer",{"metric":"alkene geometry","major_isomer":"E","reported_value":"E-isomer"},[P[3]],scope=f"geometry of {label}",product_role="reported E-isomer")
    set_class(r,"azide",[P[3]])
    r["review"]["unresolved"]=["Yield measurement basis is not stated in the inspected main-text table.","Absolute substrate amount and product-specific chromatography eluent are not reported in the supplied main text; SI is unassociated."]
    REACTIONS.append(r)


# Table 5 concrete one-pot experiments.  Product drawings are retained as
# unresolved identities rather than forcing unsafe benzodioxin SMILES.
BENZO = [
 ("5a",94,"p_fluorophenyl","salicylaldehyde"),
 ("5b",91,"p_tolyl","salicylaldehyde"),
 ("5c",92,"p_tolyl","methoxy-salicylaldehyde"),
 ("5d",90,"p_anisyl","hydroxy-naphthaldehyde"),
 ("5e",93,"naphthyl","methyl-salicylaldehyde"),
]
for label,yld,acyl,aldname in BENZO:
    ev=[P[3],P[4]]
    sub=substrate(acyl,[P[3]])
    ald=simple(f"aldehyde_{label}",aldname,[P[3]],None,aldname,note="Labelled aromatic aldehyde drawing is legible as a substitution pattern, but a concrete SMILES is withheld pending expert positional review.")
    prod=molecule(f"product_{label}",f"{label} benzodioxin",[P[3]],None,f"benzodioxin {label}",note="Concrete benzodioxin drawing retained as image evidence; machine did not safely encode the fused-ring/substitution connectivity.")
    r=base_reaction(f"table5_{label}",ev,"specific_experiment",["substrate_scope"],f"Table 5, {label}","benzodioxin skeleton synthesis","one-pot benzodioxin annulation")
    r["general_scheme_ids"]=[GEN_BENZO]
    r["participants"]=[participant("substrate",sub,"reactant",ev,eq=1.0,eq_ref=sub),participant("aldehyde",ald,"reactant",ev,eq=1.1,eq_ref=sub),participant("li2co3",common["li2co3"],"reagent",ev,eq=1.2,eq_ref=sub),participant("oxone",common["oxone"],"reagent",ev,eq=1.2,eq_ref=sub),participant("solvent",common["dmf"],"solvent",ev,volume=("5 mL",{"number":5,"unit":"mL"})),participant("product",prod,"product",ev,product_role="principal_product")]
    r["solvent_conditions"]["composition"]=F("DMF","DMF",ev)
    r["solvent_conditions"]["total_volume"]=F("5 mL",{"number":5,"unit":"mL"},[P[4]],"inherited_general_condition")
    r["steps"]=[reaction_step("reaction",ev,["substrate","aldehyde","li2co3","oxone","solvent"],25,12,"open atmosphere","General Procedure (B): combine allenic ketone in DMF, Li2CO3, substituted ortho-hydroxybenzaldehyde and Oxone; hold at room temperature")]
    add_workup(r,[P[4]],label)
    r["procedure"]["general_procedure_reference"]=F("General Procedure (B)","General Procedure (B)",[P[4]],"inherited_general_condition")
    r["procedure"]["addition_order"]=F("allenic ketone in DMF followed by Li2CO3, substituted ortho-hydroxy benzaldehyde and oxone","substrate in DMF; add Li2CO3; add aldehyde; add Oxone",[P[4]],"inherited_general_condition")
    r["procedure"]["pretreatment"]=F("oven dried single neck round-bottomed flask (10 mL capacity)","oven-dried 10 mL single-neck round-bottom flask",[P[4]],"inherited_general_condition")
    add_outcome(r,"yield","yield",[prod],f"{yld}%",{"number":yld,"unit":"%"},[P[3]],scope=f"Table 5 one-pot product {label}")
    set_class(r,None,[P[3]])
    r["review"]["unresolved"]=["Yield measurement basis is not stated in Table 5.","Benzodioxin product and aldehyde substitution positions need expert structure review; no concrete SMILES asserted.","Absolute substrate amount and product-specific chromatography eluent are not reported; SI is unassociated."]
    REACTIONS.append(r)


def oxidation_product(mid,label,start_label,level,aryl,ev):
    acyl="p_tolyl"
    x = f"S(=O){aryl}" if level==1 else f"S(=O)(=O){aryl}"
    return e_product(mid,label,acyl,x,ev,None,"E")

prod_6a=oxidation_product("product_6a_E","6a (E sulfoxide)","2m",1,"c1ccc(C)cc1",[P[4]])
prod_6b=oxidation_product("product_6b_E","6b (E sulfoxide)","2n",1,"c1ccc(Br)cc1",[P[4]])
prod_7=oxidation_product("product_7_E","7 (E sulfone)","2m",2,"c1ccc(C)cc1",[P[4]])
prod_2m=MOLECULES["product_2m_E"]["molecule_id"]
prod_2n=MOLECULES["product_2n_E"]["molecule_id"]
for label,start,prod,eq,yld in [("6a",prod_2m,prod_6a,1.1,90),("6b",prod_2n,prod_6b,1.1,85),("7",prod_2m,prod_7,2.2,90)]:
    ev=[P[4]]
    r=base_reaction(f"scheme4_{label}",ev,"specific_experiment",["derivatization"],f"Scheme 4, {label}","post-functionalization with Oxone","sulfide oxidation")
    r["participants"]=[participant("substrate",start,"reactant",ev,eq=1.0,eq_ref=start),participant("oxone",common["oxone"],"reagent",ev,eq=eq,eq_ref=start),participant("solvent",common["dmf"],"solvent",ev),participant("product",prod,"product",ev,product_role="principal_product")]
    r["solvent_conditions"]["composition"]=F("DMF","DMF",ev)
    r["steps"]=[reaction_step("oxidation",ev,["substrate","oxone","solvent"],"r.t.",12,None,"treat beta-vinyl thioether with Oxone in DMF")]
    add_outcome(r,"yield","yield",[prod],f"{yld}%",{"number":yld,"unit":"%"},ev,scope=f"Scheme 4 product {label}")
    changes=[]
    count=1 if label in ("6a","6b") else 2
    for i in range(count):
        changes.append({"change_id":f"b{i+1}","atoms":[{"atom_map":None,"element":"S"},{"atom_map":None,"element":"O"}],"before_order":0,"after_order":2,"change_type":"formed"})
    labels=[{"category_number":18,"code":"DB-O-S","supporting_change_ids":[c["change_id"] for c in changes]}]
    deriv="S=O formation inferred from the Scheme 4 sulfide/sulfoxide/sulfone drawings; oxygen atoms are not fully mapped to Oxone."
    r["classification"]["bond_changes"]=[F(None,c,ev,"structure_derivation",deriv) for c in changes]
    r["classification"]["labels"]=F(None,labels,ev,"structure_derivation",deriv)
    r["classification"]["coverage"]=F(None,"partial",ev,"inference",deriv)
    r["classification"]["limitations"]=["Oxone oxygen provenance and coproducts are not mapped."]
    r["review"]["unresolved"]=["Yield measurement basis, reaction scale, atmosphere and work-up are not reported for this Scheme 4 arrow."]
    REACTIONS.append(r)

# 3m -> 5b is a separate 94% post-functionalization, not the 91% one-pot Table 5 entry.
prod_3m=MOLECULES["product_3m_E"]["molecule_id"]
prod_5b=MOLECULES["product_5b"]["molecule_id"]
r=base_reaction("scheme4_3m_to_5b",[P[4]],"specific_experiment",["derivatization"],"Scheme 4, 3m to 5b","post-functionalization of beta-vinyl ether 3m","oxidative cyclization to benzodioxin")
r["participants"]=[participant("substrate",prod_3m,"reactant",[P[4]],eq=1.0,eq_ref=prod_3m),participant("oxone",common["oxone"],"reagent",[P[4]],eq=1.2,eq_ref=prod_3m),participant("solvent",common["dmf"],"solvent",[P[4]]),participant("product",prod_5b,"product",[P[4]],product_role="principal_product")]
r["solvent_conditions"]["composition"]=F("DMF","DMF",[P[4]])
r["steps"]=[reaction_step("cyclization",[P[4]],["substrate","oxone","solvent"],"r.t.",12,None,"treat 3m with Oxone in DMF")]
add_outcome(r,"yield","yield",[prod_5b],"94%",{"number":94,"unit":"%"},[P[4]],scope="Scheme 4 post-functionalization 3m to 5b")
set_class(r,None,[P[4]])
r["review"]["unresolved"]=["Yield measurement basis, reaction scale, atmosphere, work-up and exact net atom mapping are not reported/resolved."]
REACTIONS.append(r)

# Azide thermolysis gives two concrete oxazole examples.
for label,start_label,aryl,yld in [("8a","4a","c1ccccc1C",75),("8b","4c","c1ccc(F)cc1",74)]:
    start=MOLECULES[f"product_{start_label}_E"]["molecule_id"]
    prod=simple(f"product_{label}",f"{label} oxazole",[P[4]],f"Cc1nc(-{aryl})oc1",f"2-methyl-5-aryl oxazole {label}",note="Candidate oxazole connectivity transcribed from Scheme 4; expert review remains required.")
    r=base_reaction(f"scheme4_{label}",[P[4]],"specific_experiment",["derivatization"],f"Scheme 4, {label}","thermolysis of beta-vinyl azide","2,5-disubstituted oxazole synthesis")
    r["participants"]=[participant("substrate",start,"reactant",[P[4]]),participant("solvent",common["toluene"],"solvent",[P[4]]),participant("product",prod,"product",[P[4]],product_role="principal_product")]
    r["solvent_conditions"]["composition"]=F("toluene","toluene",[P[4]])
    r["steps"]=[reaction_step("thermolysis",[P[4]],["substrate","solvent"],120,None,"N2","heat beta-vinyl azide in toluene via a drawn 2H-azirine intermediate")]
    add_outcome(r,"yield","yield",[prod],f"{yld}%",{"number":yld,"unit":"%"},[P[4]],scope=f"Scheme 4 product {label}")
    set_class(r,None,[P[4]])
    r["review"]["unresolved"]=["Reaction duration, yield measurement basis, scale, work-up and complete atom mapping are not reported/resolved."]
    REACTIONS.append(r)


# Non-synthetic observations remain source context.
BUNDLE["source_context"] = [
    {"context_id": f"{DOC}:figure1_biological_structures", "document_id": DOC,
     "evidence_ids": [P[1]], "kind": "background_structures",
     "note": "Figure 1 shows biologically important benzodioxin-containing structures; it is not a synthesis experiment in this paper."},
    {"context_id": f"{DOC}:figure2_noe", "document_id": DOC,
     "evidence_ids": [P[3]], "kind": "stereochemical_characterization",
     "targets": ["product_2f_E", "product_2s_E", "product_3d_E", "product_4j_E", "product_4l_E"],
     "note": "Multiple nOe contacts are drawn for 2f, 2s, 3d, 4j and 4l to support relative product geometry; Figure 2 is characterization evidence, not five additional syntheses."},
    {"context_id": f"{DOC}:page5_references", "document_id": DOC,
     "evidence_ids": [P[5]], "kind": "references_only",
     "note": "Page 5 contains bibliography continuation and received date; no reaction arrow, table row, or experimental procedure is present."},
]

BUNDLE["molecules"] = list(MOLECULES.values())
BUNDLE["reactions"] = REACTIONS
(RUN / "submitted.json").write_text(json.dumps(BUNDLE, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/stereo_checks.json").write_text(json.dumps(STEREO_CHECKS, ensure_ascii=False, indent=2) + "\n")

page_ids = {n: [] for n in range(1, 6)}
for r in REACTIONS:
    for n in range(1, 6):
        if P[n] in r["evidence_ids"]:
            page_ids[n].append(r["reaction_id"])

inventory = {
    "document_id": DOC,
    "pages": {
        "1": {"findings": ["Scheme 1 cited prior-work arrow with four products", "Scheme 2 generic thioether/ether/azide arrows", "Figure 1 background structures", "optimization discussion"], "record_ids": page_ids[1]},
        "2": {"findings": ["Scheme 3 seven cited strategy arrows plus current generic route", "Table 1 entries 1-17", "Table 2 entries 2a-2x"], "record_ids": page_ids[2]},
        "3": {"findings": ["Table 3 entries 3a-3o", "Table 4 entries 4a-4l", "Figure 2 nOe characterization", "Table 5 entries 5a-5e and proposed cascade"], "record_ids": page_ids[3]},
        "4": {"findings": ["Scheme 4 six concrete post-functionalization examples", "General Procedures A and B"], "record_ids": page_ids[4]},
        "5": {"findings": ["references only; no reaction record"], "record_ids": page_ids[5]},
    },
    "counts": {"molecules": len(MOLECULES), "reactions": len(REACTIONS)},
    "limitations": ["No SI was explicitly supplied or associated.", "Main-text yield bases are often not stated.", "All structures and classifications are machine candidates requiring expert review."],
}
(RUN / "raw_model/page_visual_inventory.json").write_text(json.dumps(inventory, ensure_ascii=False, indent=2) + "\n")

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "upper title and abstract", "finding": "Article title/DOI and current-work summary; no separate experiment arrow.", "evidence_ids": [P[1]]},
            {"region": "lower-left Scheme 1", "finding": "One cited Sugita prior-work reaction arrow leading to four drawn isomers; recorded only as literature_example.", "evidence_ids": [P[1]]},
            {"region": "upper-right Scheme 2", "finding": "Three current generic transformations: thioether, ether and azide additions; no concrete scope row counted here.", "evidence_ids": [P[1]]},
            {"region": "middle-right Figure 1", "finding": "Background benzodioxin-containing structures; non-synthetic observation retained in source_context.", "evidence_ids": [P[1]]},
            {"region": "bottom-right optimization discussion", "finding": "Narrative identifies 1f/p-thiocresol and optimum entry 7; concrete rows are enumerated on page 2.", "evidence_ids": [P[1]]},
        ], "record_ids": page_ids[1], "unresolved": ["Supporting information is mentioned but no SI file was explicitly supplied or associated."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "upper-left Scheme 3", "finding": "Seven distinct cited benzodioxin strategy arrows and the current generic route; prior work kept separate from current experiments.", "evidence_ids": [P[2]]},
            {"region": "middle-left Table 1", "finding": "All 17 optimization rows transcribed independently, including solvent/base/temperature/yield and 2f:2f′ ratio; E/Z method is 1H NMR.", "evidence_ids": [P[2]]},
            {"region": "right Table 2", "finding": "All 24 concrete thioether scope products 2a–2x transcribed as separate experiments under the printed generic conditions.", "evidence_ids": [P[2]]},
        ], "record_ids": page_ids[2], "unresolved": ["Table 1 and Table 2 yield measurement bases are not stated.", "Several miniature substituent positions/chain lengths remain machine structure candidates requiring expert review."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "left upper Table 3", "finding": "All 15 ether scope products 3a–3o transcribed independently.", "evidence_ids": [P[3]]},
            {"region": "left lower Table 4", "finding": "All 12 azide scope products 4a–4l transcribed independently, with no external base inherited.", "evidence_ids": [P[3]]},
            {"region": "right upper Figure 2", "finding": "nOe characterization for 2f, 2s, 3d, 4j and 4l retained in source_context and not counted as synthesis experiments.", "evidence_ids": [P[3]]},
            {"region": "right lower Table 5", "finding": "Generic one-pot annulation, proposed dashed cascade, and five concrete benzodioxin products 5a–5e transcribed; proposed pathway is not treated as a demonstrated mechanism experiment.", "evidence_ids": [P[3]]},
        ], "record_ids": page_ids[3], "unresolved": ["Tables 3–5 do not state yield measurement bases.", "Benzodioxin and several phenoxy substitution patterns were not forced into concrete SMILES where positional encoding was unsafe."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "left upper Scheme 4", "finding": "Six concrete post-functionalization examples recorded: 6a, 6b, 7, a distinct 3m→5b execution, 8a and 8b.", "evidence_ids": [P[4]]},
            {"region": "left lower General Procedure A", "finding": "Shared addition order, THF volume, open atmosphere, 35 °C/4–7 h, work-up and silica chromatography inherited by Tables 2–4; Li2CO3 explicitly excepted for azides.", "evidence_ids": [P[4]]},
            {"region": "right upper General Procedure B", "finding": "Shared DMF volume, reagent equivalents, open atmosphere, 25 °C/12 h, work-up and purification inherited by Table 5.", "evidence_ids": [P[4]]},
        ], "record_ids": page_ids[4], "unresolved": ["Scheme 4 omits yield measurement bases, absolute scales and most work-up details.", "General procedures give no absolute substrate amount or product-specific chromatography eluent."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "both columns", "finding": "Bibliography continuation and received date only; no reaction arrow, table row, control, mechanism experiment or procedure.", "evidence_ids": [P[5]]},
        ], "record_ids": page_ids[5], "unresolved": []},
    ],
}
(RUN / "coverage").mkdir(parents=True, exist_ok=True)
(RUN / "coverage/pages_0001_0005.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")
print(json.dumps(inventory["counts"], ensure_ascii=False))
