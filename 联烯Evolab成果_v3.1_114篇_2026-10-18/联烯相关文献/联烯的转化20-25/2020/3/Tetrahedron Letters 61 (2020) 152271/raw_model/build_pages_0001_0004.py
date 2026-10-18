#!/usr/bin/env python3
"""Source-derived v3.1 extraction for Nakao et al., PDF pages 1--4.

This script is a transparent extraction artifact.  It uses only the prepared
document registry and the four visually inspected page images/text evidence.
"""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills" / "chemical-reaction-extraction" / "scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs" / "current"
DATASET = json.loads((RUN / "dataset.json").read_text())
DOC = "doc_5d57d6fc4002ccdc"
P1 = f"{DOC}:page1:text"
P2 = f"{DOC}:page2:text"
P3 = f"{DOC}:page3:text"
P4 = f"{DOC}:page4:text"
IMG1 = f"{DOC}:page1:image"
IMG2 = f"{DOC}:page2:image"
IMG3 = f"{DOC}:page3:image"
IMG4 = f"{DOC}:page4:image"


def F(raw=None, normalized=None, ev=(), method="original_extraction", derivation=None,
      status=None):
    return st.fact(raw, normalized, ev, method=method, derivation=derivation,
                   status=status)


def D(normalized, ev, derivation, raw=None):
    return F(raw, normalized, ev, method="structure_derivation", derivation=derivation)


def I(normalized, ev, derivation, raw=None):
    return F(raw, normalized, ev, method="inference", derivation=derivation)


def q(number, unit, raw_text=None, qualifier=None):
    out = {"number": number, "unit": unit}
    if raw_text is not None:
        out["raw_text"] = raw_text
    if qualifier is not None:
        out["qualifier"] = qualifier
    return out


def mid(label):
    return f"{DOC}:mol:{label}"


def rid(label):
    return f"{DOC}:{label}"


bundle = st.empty(DATASET["documents"])
molecules = {}


def set_status(m, field, value, ev, why):
    m["stereochemistry"][field] = I(value, ev, why)


def concrete(label, name, smiles, ev, stereo="simple"):
    molecule_id = mid(label)
    m = st.new_molecule(molecule_id, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"] = [{
        "document_id": DOC,
        "compound_label": F(label, label, ev),
        "evidence_ids": list(ev),
    }]
    m["names"]["original"] = F(name, name, ev)
    m["names"]["english"] = F(name, name, ev)
    m["representation_type"] = "smiles"
    m["smiles"] = D(
        smiles, ev,
        "Connectivity SMILES was transcribed from the displayed structure or explicit chemical name; "
        "RDKit parsing is only a graph check and not independent source verification.",
    )
    m["generic_template"] = None
    m["identity_note"] = F(
        {"source_label": label, "structure_source": "displayed structure and/or explicit name"},
        "Source-identified identity; machine structure remains needs review.", ev,
    )
    diag = st.stereo_inspect(smiles)
    s = m["stereochemistry"]
    s["mixture"] = F()
    if stereo == "allene_unspecified":
        s["axial"] = I(
            [{"axis": "cumulated C=C=C", "assignment": "unspecified"}], ev,
            "The displayed allene has different substituents at both terminal carbons, so axial "
            "stereogenicity is possible; the source gives no R_a/S_a or P/M assignment.",
        )
        s["atom_centers"] = I([], ev, "No tetrahedral stereocenter is asserted for the displayed substrate.")
        s["double_bonds"] = I([], ev, "Ordinary E/Z nomenclature is not assigned to the cumulated allene axis.")
        s["assignment"] = F(
            {"source_marker": f"displayed compound {label}", "stereo_label": "none"},
            "Axial configuration unassigned by the source.", ev,
        )
        set_status(m, "chirality_status", "unspecified", ev,
                   "Potential allene axial chirality is present but unassigned.")
        set_status(m, "geometry_status", "not_applicable", ev,
                   "The relevant element is an allene axis, not an ordinary E/Z alkene.")
        set_status(m, "representation_status", "not_encoded", ev,
                   "The connectivity SMILES does not encode the unassigned allene axis.")
    elif stereo == "allene_symmetric":
        s["axial"] = I(
            [{"axis": "cumulated C=C=C", "assignment": "not stereogenic because one terminal carbon bears two phenyl groups"}],
            ev, "The two substituents on one terminal allene carbon are both phenyl.")
        s["atom_centers"] = I([], ev, "No tetrahedral stereocenter is present in the candidate graph.")
        s["double_bonds"] = I([], ev, "No ordinary E/Z element is assigned to the cumulated allene.")
        s["assignment"] = I("no stereogenic allene axis", ev,
                              "Identical phenyl substituents remove axial stereogenicity.")
        set_status(m, "chirality_status", "achiral", ev,
                   "Identical phenyl substituents remove the possible allene axis.")
        set_status(m, "geometry_status", "not_applicable", ev,
                   "The cumulated system is not assigned ordinary E/Z geometry.")
        set_status(m, "representation_status", "partial", ev,
                   "Connectivity is encoded; toolkit diagnostics report cumulated double bonds but do not certify axial stereochemistry.")
    elif stereo == "center_unspecified":
        center = next((x for x in diag["atom_centers"] if x["CIP"] == "?"), None)
        s["atom_centers"] = I(
            [{"atom_index": center["atom_index"] if center else None,
              "assignment": "unspecified",
              "description": "carbon attached to Ar, R, H, and thiophene"}], ev,
            "The displayed product has four different substituents at this carbon, but no wedge/dash or R/S label.",
        )
        s["double_bonds"] = I([], ev, "No independently assigned E/Z alkene is present in this aromatic thiophene product.")
        s["axial"] = I([], ev, "No axial stereogenic element is asserted for this product.")
        s["assignment"] = F(
            {"source_marker": f"displayed compound {label}", "wedge_or_RS": "not shown"},
            "Tetrahedral configuration unspecified.", ev,
        )
        set_status(m, "chirality_status", "unspecified", ev,
                   "A potential stereocenter exists but the source does not assign it.")
        set_status(m, "geometry_status", "not_applicable", ev,
                   "No ordinary E/Z element is present.")
        set_status(m, "representation_status", "partial", ev,
                   "Connectivity is encoded while the unassigned tetrahedral configuration is intentionally omitted.")
    elif stereo == "center_symmetric":
        s["atom_centers"] = I([], ev, "The candidate carbon bears two identical phenyl substituents and is not stereogenic.")
        s["double_bonds"] = I([], ev, "No independently assigned E/Z alkene is present.")
        s["axial"] = I([], ev, "No axial stereogenic element is present.")
        s["assignment"] = I("no stereogenic center", ev,
                              "The candidate carbon has two identical phenyl substituents.")
        set_status(m, "chirality_status", "achiral", ev,
                   "The candidate graph has no potential stereocenter.")
        set_status(m, "geometry_status", "not_applicable", ev,
                   "No E/Z element is present.")
        set_status(m, "representation_status", "encoded", ev,
                   "The complete non-stereogenic graph is represented by the SMILES.")
        s["isomeric_smiles"] = D(diag["canonical_isomeric_smiles"], ev,
                                   "Generated by standard stereo graph diagnostics from the candidate SMILES.")
        s["identity_key"] = D(diag["stereo_graph_key"], ev,
                                "Generated by standard stereo graph diagnostics; not a sample identity.")
    elif stereo == "alkene_relative":
        s["atom_centers"] = I([], ev, "No tetrahedral stereocenter is asserted for 5a.")
        s["axial"] = I([], ev, "No axial stereogenic element is present in 5a.")
        s["double_bonds"] = F(
            {"source": "drawn 5a and statement that geometry was determined by 1H-1H NOEs",
             "relative_drawing": "Ph and CH2CO2Me are drawn on the same side; Me and SCH2CO2Me are drawn on the same side"},
            [{"bond": "central alkene of 5a", "assignment": "relative geometry shown; E/Z not explicitly named"}], ev,
        )
        s["assignment"] = F(
            "The geometry of 5a was determined on the basis of the 1H-1H nuclear Overhauser effects (NOEs).",
            {"method": "1H-1H NOE", "absolute_EZ_label": "not stated"}, ev,
        )
        set_status(m, "chirality_status", "not_applicable", ev,
                   "No tetrahedral or axial chirality is asserted for 5a.")
        set_status(m, "geometry_status", "specified", ev,
                   "The source drawing and NOE statement specify relative alkene geometry.")
        set_status(m, "representation_status", "not_encoded", ev,
                   "The candidate SMILES deliberately omits E/Z because the source does not name E/Z and the relative drawing was not converted to CIP.")
    else:
        s["atom_centers"] = I([], ev, "Standard stereo diagnostics found no tetrahedral stereocenter in the candidate graph.")
        s["double_bonds"] = I([], ev, "Standard stereo diagnostics found no relevant E/Z element.")
        s["axial"] = I([], ev, "No axial stereogenic element is present in the candidate graph.")
        s["assignment"] = I("no stereogenic element", ev,
                              "Standard stereo graph diagnostics found no potential stereogenic element.")
        set_status(m, "chirality_status", "achiral", ev,
                   "No potential stereogenic element was found in the candidate graph.")
        set_status(m, "geometry_status", "not_applicable", ev,
                   "No relevant E/Z element is present.")
        set_status(m, "representation_status", "encoded", ev,
                   "The complete non-stereogenic candidate graph is encoded.")
        s["isomeric_smiles"] = D(diag["canonical_isomeric_smiles"], ev,
                                   "Generated by standard stereo graph diagnostics from the candidate SMILES.")
        s["identity_key"] = D(diag["stereo_graph_key"], ev,
                                "Generated by standard stereo graph diagnostics; not a sample identity.")
    molecules[molecule_id] = m
    return m


def generic(label, name, ev, fixed, variables, instances, stereo_note):
    molecule_id = mid(label)
    m = st.new_molecule(molecule_id, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"] = [{"document_id": DOC,
                          "compound_label": F(label, label, ev),
                          "evidence_ids": list(ev)}]
    m["names"]["original"] = F(name, name, ev)
    m["names"]["english"] = F(name, name, ev)
    m["representation_type"] = "generic_template"
    m["smiles"] = F()
    var_rows = []
    for vlabel, definition, values in variables:
        var_rows.append({
            "label": vlabel,
            "attachment_points": F(f"Scheme attachment for {vlabel}", f"scheme-defined attachment for {vlabel}", ev),
            "definition": F(definition, definition, ev),
            "allowed_values": F(values, values, ev),
            "constraints": F(),
        })
    mappings = []
    for target, assignments in instances:
        mappings.append({
            "molecule_id": mid(target),
            "compound_label": F(target, target, ev),
            "assignments": F(assignments, assignments, ev),
        })
    m["generic_template"] = {
        "fixed_structure": F(fixed, fixed, ev),
        "variables": var_rows,
        "instance_mapping": mappings,
        "reason": F("Scheme uses Ar and R variables.", "Generic source drawing with Ar/R variables; no unique concrete SMILES.", ev),
    }
    s = m["stereochemistry"]
    s["atom_centers"] = I(None, ev, stereo_note) if False else F()
    s["double_bonds"] = F()
    s["axial"] = I(stereo_note, ev, stereo_note)
    s["mixture"] = F()
    s["assignment"] = I("depends on Ar/R instance", ev, stereo_note)
    set_status(m, "chirality_status", "unresolved", ev, stereo_note)
    set_status(m, "geometry_status", "unresolved", ev, stereo_note)
    set_status(m, "representation_status", "not_encoded", ev,
               "A generic Ar/R template cannot be represented as a unique concrete stereochemical SMILES.")
    m["identity_note"] = F(fixed, "Generic source structure; variables retained without expansion.", ev)
    molecules[molecule_id] = m
    return m


# Concrete substrates and products.
ar = {
    "a": ("Ph", "Me", "c1ccccc1", "C"),
    "b": ("4-MeOC6H4", "Me", "c1ccc(OC)cc1", "C"),
    "c": ("4-MeC6H4", "Me", "c1ccc(C)cc1", "C"),
    "d": ("4-ClC6H4", "Me", "c1ccc(Cl)cc1", "C"),
    "e": ("4-NO2C6H4", "Me", "c1ccc([N+](=O)[O-])cc1", "C"),
    "f": ("Ph", "Et", "c1ccccc1", "CC"),
    "g": ("Ph", "Ph", "c1ccccc1", "c1ccccc1"),
}
substrate_smiles = {
    "a": "COC(=O)C=C=C(C)c1ccccc1",
    "b": "COC(=O)C=C=C(C)c1ccc(OC)cc1",
    "c": "COC(=O)C=C=C(C)c1ccc(C)cc1",
    "d": "COC(=O)C=C=C(C)c1ccc(Cl)cc1",
    "e": "COC(=O)C=C=C(C)c1ccc([N+](=O)[O-])cc1",
    "f": "COC(=O)C=C=C(CC)c1ccccc1",
    "g": "COC(=O)C=C=C(c1ccccc1)c1ccccc1",
}
side = {
    "a": "C(C)c2ccccc2",
    "b": "C(C)c2ccc(OC)cc2",
    "c": "C(C)c2ccc(C)cc2",
    "d": "C(C)c2ccc(Cl)cc2",
    "e": "C(C)c2ccc([N+](=O)[O-])cc2",
    "f": "C(CC)c2ccccc2",
    "g": "C(c2ccccc2)c2ccccc2",
}
for suffix in ar:
    ev = [P2, IMG2] if suffix not in ("a",) else [P2, P3, P4, IMG2, IMG3, IMG4]
    concrete(f"1{suffix}", f"allenyl ester 1{suffix}", substrate_smiles[suffix], ev,
             "allene_symmetric" if suffix == "g" else "allene_unspecified")
    concrete(f"3{suffix}", f"2,3,5-trisubstituted thiophene 3{suffix}",
             f"COC(=O)c1sc({side[suffix]})cc1O", ev,
             "center_symmetric" if suffix == "g" else "center_unspecified")
    concrete(f"4{suffix}", f"2,3,4-trisubstituted thiophene 4{suffix}",
             f"COC(=O)c1c({side[suffix]})scc1O", ev,
             "center_symmetric" if suffix == "g" else "center_unspecified")

concrete("2", "methyl mercaptoacetate", "COC(=O)CS", [P2, P3, P4, IMG2, IMG3, IMG4])
concrete("5a", "thia-Michael adduct 5a", "COC(=O)CC(=C(C)c1ccccc1)SCC(=O)OC",
         [P2, P3, IMG2, IMG3], "alkene_relative")

# Named agents, solvents, workup and purification materials.
agent_defs = {
    "DBU": ("1,8-diazabicyclo[5.4.0]undec-7-ene (DBU)", "C1CCC2=NCCCN2CC1", [P1, P2, P3, P4, IMG1, IMG2, IMG3, IMG4]),
    "Et3N": ("triethylamine", "CCN(CC)CC", [P1, P2, P3, P4, IMG1, IMG2, IMG3, IMG4]),
    "LiCl": ("lithium chloride", "[Li+].[Cl-]", [P1, P2, IMG1, IMG2]),
    "MgBr2": ("magnesium bromide", "[Mg+2].[Br-].[Br-]", [P1, P2, P3, P4, IMG1, IMG2, IMG3, IMG4]),
    "THF": ("tetrahydrofuran (THF)", "C1CCOC1", [P2, P3, P4, IMG2, IMG3, IMG4]),
    "CH2Cl2": ("dichloromethane", "ClCCl", [P2, IMG2]),
    "dioxane": ("1,4-dioxane", "O1CCOCC1", [P2, IMG2]),
    "HCl": ("hydrochloric acid", "Cl", [P4, IMG4]),
    "CHCl3": ("chloroform", "ClC(Cl)Cl", [P4, IMG4]),
    "MgSO4": ("anhydrous magnesium sulfate", "[Mg+2].[O-]S(=O)(=O)[O-]", [P4, IMG4]),
    "silica": ("Silica Gel PSQ 60B", "O=[Si]=O", [P4, IMG4]),
    "hexane": ("n-hexane", "CCCCCC", [P4, IMG4]),
    "EtOAc": ("ethyl acetate (AcOEt)", "CCOC(C)=O", [P4, IMG4]),
}
for label, (name, smi, ev) in agent_defs.items():
    concrete(label, name, smi, ev)

# Generic compounds in Scheme 1 and intermediates in Scheme 4.
instances_1 = [(f"1{s}", {"Ar": ar[s][0], "R": ar[s][1]}) for s in ar]
instances_3 = [(f"3{s}", {"Ar": ar[s][0], "R": ar[s][1]}) for s in ar]
instances_4 = [(f"4{s}", {"Ar": ar[s][0], "R": ar[s][1]}) for s in ar]
variables = [("Ar", "aryl substituent shown as Ar", [x[0] for x in ar.values()]),
             ("R", "substituent shown as R", [x[1] for x in ar.values()])]
generic("g1", "allenyl esters 1 (generic)", [P2, IMG2],
        "Ar,R-substituted allenyl methyl ester shown in Scheme 1", variables, instances_1,
        "Allene axial chirality depends on whether both terminal substituent pairs differ; no generic assignment is possible.")
generic("g3", "2,3,5-trisubstituted thiophenes 3 (generic)", [P2, P3, IMG2, IMG3],
        "Scheme-1 2,3,5-trisubstituted hydroxythiophene bearing Ar/R substituent and CO2Me", variables, instances_3,
        "Tetrahedral chirality at the Ar/R-bearing carbon depends on whether Ar and R differ and is not generically assigned.")
generic("g4", "2,3,4-trisubstituted thiophenes 4 (generic)", [P2, P3, IMG2, IMG3],
        "Scheme-1 2,3,4-trisubstituted hydroxythiophene bearing Ar/R substituent and CO2Me", variables, instances_4,
        "Tetrahedral chirality at the Ar/R-bearing carbon depends on whether Ar and R differ and is not generically assigned.")
generic("g5", "thia-Michael adducts 5 (generic)", [P2, P3, IMG2, IMG3],
        "Scheme-1 thia-Michael adduct retaining two methyl ester groups", variables,
        [("5a", {"Ar": "Ph", "R": "Me"})],
        "Alkene geometry may vary; only 5a has a source NOE assignment in the supplied main text.")
for label, name in [("A", "kinetic enolate A"), ("B", "thermodynamic enolate B"),
                    ("C", "magnesium enolate C")]:
    generic(label, name, [P3, IMG3], f"Proposed intermediate {label} drawn in Scheme 4",
            variables, [], "Proposed generic enolate; stereochemical identity is not established experimentally here.")

bundle["molecules"] = list(molecules.values())


def charge(raw, norm, ev):
    return F(raw, norm, ev)


def participant(pid, molecule, role, ev, charges=None, ref=None, stage="reaction", product_role=None):
    p = st.new_participant(pid, mid(molecule), role, ev, stage=stage)
    for key, spec in (charges or {}).items():
        p["charges"][key] = charge(spec[0], spec[1], ev)
    if ref is not None:
        p["equivalence_reference"] = F(f"relative to {ref}", mid(ref), ev)
    if product_role is not None:
        p["product_role"] = F(product_role, product_role, ev)
    return p


def step(step_id, ev, stage, pids, operation, temperature=None, time=None,
         atmosphere=None, pressure=None):
    s = st.new_step(step_id, ev, stage=stage)
    s["participant_ids"] = list(pids)
    s["operation"] = F(operation, operation, ev)
    if temperature is not None:
        s["temperature"] = F(temperature[0], temperature[1], ev)
    if time is not None:
        s["time"] = F(time[0], time[1], ev)
    if atmosphere is not None:
        s["atmosphere"] = F(atmosphere, atmosphere, ev)
    if pressure is not None:
        s["pressure"] = F(pressure, pressure, ev)
    return s


def outcome(oid, kind, targets, ev, raw, norm, product_role=None, yield_type=None,
            basis=None, scope=None):
    o = st.new_outcome(oid, kind, [mid(x) for x in targets], ev)
    o["value"] = F(raw, norm, ev)
    if product_role is not None:
        o["product_role"] = F(product_role, product_role, ev)
    if yield_type is not None:
        o["yield_type"] = F(yield_type, yield_type, ev)
    if basis is not None:
        o["measurement_basis"] = F(basis, basis, ev)
    if scope is not None:
        o["scope"] = F(scope, scope, ev)
    return o


def base_reaction(name, ev, record_type, categories, location, reaction_name,
                  transformation, unresolved):
    r = st.new_reaction(rid(name), DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = F(categories, categories, ev)
    r["location"] = F(location, location, ev)
    r["reaction_name"] = F(reaction_name, reaction_name, ev)
    r["transformation_type"] = I(transformation, ev,
        "Transformation description is a machine normalization of the source reaction name/drawing.")
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": list(unresolved)}
    r["completeness"] = "partial"
    return r


def labels_for(changes):
    labels = {}
    for ch in changes:
        if ch["change_type"] in ("formed", "order_increased", "hydrogen_gain"):
            lab = st.bond_label(ch["atoms"][0]["element"], ch["atoms"][1]["element"], ch["after_order"])
            labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(ch["change_id"])
    return list(labels.values())


def classify_manual(r, changes, ev, limitation):
    c = r["classification"]
    c["mapped_reaction_smiles"] = F()
    c["bond_changes"] = [D(ch, ev,
        "Manual graph-difference candidate from the source structures; atom mapping is not asserted complete.") for ch in changes]
    c["labels"] = D(labels_for(changes), ev,
        "Taxonomy labels were calculated from the explicitly listed candidate bond changes.")
    c["coverage"] = I("partial", ev,
        "Only source-supported net graph changes are asserted; complete atom mapping, proton source, and coproducts remain unresolved.")
    c["limitations"] = [limitation,
        "No atom-mapped reaction SMILES is claimed.",
        "Hydrogen gains/losses are net changes; proton donor/acceptor identity is not established."]
    c["model_votes"] = []


def classify_none(r, ev):
    c = r["classification"]
    c["mapped_reaction_smiles"] = F()
    c["bond_changes"] = []
    c["labels"] = I([], ev, "No observed transformation; no bond-formation label is assigned.")
    c["coverage"] = I("not_applicable", ev, "The source reports no reaction under these conditions.")
    c["limitations"] = ["Negative control; intended products are not used to infer a classification."]
    c["model_votes"] = []


def classify_unresolved(r, ev, note):
    c = r["classification"]
    c["mapped_reaction_smiles"] = F()
    c["bond_changes"] = []
    c["labels"] = F()
    c["coverage"] = I("unresolved", ev, note)
    c["limitations"] = [note]
    c["model_votes"] = []


def ch(cid, a, b, before, after, kind, count=None, branch=None):
    out = {"change_id": cid, "atoms": [{"atom_map": None, "element": a},
                                          {"atom_map": None, "element": b}],
           "before_order": before, "after_order": after, "change_type": kind}
    if count is not None:
        out["hydrogen_count"] = count
    if branch is not None:
        out["branch"] = branch
    return out


changes3 = [
    ch("p3_cs", "C", "S", 0, 1, "formed"),
    ch("p3_cc", "C", "C", 0, 2, "formed"),
    ch("p3_ch", "C", "H", 0, 1, "hydrogen_gain", 1),
    ch("p3_oh", "O", "H", 0, 1, "hydrogen_gain", 1),
    ch("p3_ccdec", "C", "C", 2, 1, "order_decreased"),
    ch("p3_codec", "C", "O", 2, 1, "order_decreased"),
    ch("p3_sh", "S", "H", 1, 0, "hydrogen_loss", 1),
]
changes4 = [
    ch("p4_cs", "C", "S", 0, 1, "formed"),
    ch("p4_cc", "C", "C", 0, 1, "formed"),
    ch("p4_ccinc", "C", "C", 1, 2, "order_increased"),
    ch("p4_ch", "C", "H", 0, 1, "hydrogen_gain", 1),
    ch("p4_oh", "O", "H", 0, 1, "hydrogen_gain", 1),
    ch("p4_ccdec", "C", "C", 2, 1, "order_decreased"),
    ch("p4_codec", "C", "O", 2, 1, "order_decreased"),
    ch("p4_sh", "S", "H", 1, 0, "hydrogen_loss", 1),
]
changes5 = [
    ch("p5_cs", "C", "S", 0, 1, "formed"),
    ch("p5_ch", "C", "H", 0, 1, "hydrogen_gain", 1),
    ch("p5_ccdec", "C", "C", 2, 1, "order_decreased"),
    ch("p5_sh", "S", "H", 1, 0, "hydrogen_loss", 1),
]
changes_cyclize3 = [
    ch("cyc_cc_new", "C", "C", 0, 2, "formed"),
    ch("cyc_cc_promote", "C", "C", 1, 2, "order_increased"),
    ch("cyc_ch_gain", "C", "H", 0, 1, "hydrogen_gain", 1),
    ch("cyc_oh", "O", "H", 0, 1, "hydrogen_gain", 1),
    ch("cyc_ch_loss", "C", "H", 1, 0, "hydrogen_loss", 1),
    ch("cyc_ccdec", "C", "C", 2, 1, "order_decreased"),
    ch("cyc_codec", "C", "O", 2, 1, "order_decreased"),
]
changes_cyclize4 = [
    ch("cyc4_cc_new", "C", "C", 0, 1, "formed"),
    ch("cyc4_cc_promote1", "C", "C", 1, 2, "order_increased"),
    ch("cyc4_cc_promote2", "C", "C", 1, 2, "order_increased"),
    ch("cyc4_ch_gain", "C", "H", 0, 1, "hydrogen_gain", 1),
    ch("cyc4_oh", "O", "H", 0, 1, "hydrogen_gain", 1),
    ch("cyc4_ch_loss", "C", "H", 1, 0, "hydrogen_loss", 2),
    ch("cyc4_ccdec", "C", "C", 2, 1, "order_decreased"),
    ch("cyc4_codec", "C", "O", 2, 1, "order_decreased"),
]


reactions = []

# Scheme 1 generic transformation.
r = base_reaction("scheme1_general", [P2, IMG2], "generic_scheme", [],
                  {"page": 2, "region": "Scheme 1"},
                  "Tandem thia-Michael/Dieckmann condensation",
                  "thia-Michael addition followed by regioselective Dieckmann condensation",
                  ["Generic Ar/R structure is not a concrete experiment.",
                   "Branch-specific atom mapping is unresolved."])
r["participants"] = [
    participant("g1", "g1", "reactant", [P2, IMG2]),
    participant("g2", "2", "reactant", [P2, IMG2]),
    participant("g5", "g5", "product", [P2, IMG2], product_role="intermediate"),
    participant("g3", "g3", "product", [P2, IMG2], product_role="2,3,5-regioisomer"),
    participant("g4", "g4", "product", [P2, IMG2], product_role="2,3,4-regioisomer"),
]
r["procedure"]["original_text"] = F(
    "Scheme 1: allenyl esters 1 + methyl mercaptoacetate (2) -> thia-Michael adducts -> thiophenes 3 or 4",
    "Generic two-stage reaction map; no single set of experimental conditions.", [P2, IMG2])
classify_manual(r, changes3 + changes4, [P2, IMG2],
                "Union of the two displayed regioisomer branches; individual branch mapping remains partial.")
reactions.append(r)


table1 = [
    (1, ["Et3N"], 1.0, "THF", 0, None, None),
    (2, ["DBU"], 1.0, "THF", 99, (97, 3), "400 MHz 1H NMR (CDCl3)"),
    (3, ["DBU"], 1.0, "CH2Cl2", 81, (93, 7), "400 MHz 1H NMR (CDCl3)"),
    (4, ["DBU"], 1.0, "dioxane", 70, (97, 3), "500 MHz 1H NMR (CDCl3)"),
    (5, ["LiCl", "Et3N"], 1.2, "THF", 77, (0, 100), "500 MHz 1H NMR (CDCl3)"),
    (6, ["LiCl", "DBU"], 1.2, "THF", 91, (31, 69), "500 MHz 1H NMR (CDCl3)"),
    (7, ["MgBr2", "Et3N"], 1.2, "THF", 85, (0, 100), "500 MHz 1H NMR (CDCl3)"),
    (8, ["MgBr2", "DBU"], 1.2, "THF", 80, (0, 100), "500 MHz 1H NMR (CDCl3)"),
]


def detailed_participants(entry, ev, x, solvent, agents):
    if entry == 2:
        ps = [
            participant("substrate", "1a", "reactant", ev,
                        {"mass": ("53.8 mg", q(53.8, "mg")), "amount": ("0.286 mmol", q(0.286, "mmol")),
                         "equivalents": ("reference substrate", q(1, "equiv"))}, ref="1a"),
            participant("thiol", "2", "reactant", ev,
                        {"volume": ("26.0 uL", q(26.0, "uL")), "amount": ("0.286 mmol", q(0.286, "mmol")),
                         "equivalents": ("1 mol eq", q(1, "equiv"))}, ref="1a"),
            participant("DBU", "DBU", "reagent", ev,
                        {"volume": ("42.7 uL", q(42.7, "uL")), "amount": ("0.286 mmol", q(0.286, "mmol")),
                         "equivalents": ("1 mol eq", q(1, "equiv"))}, ref="1a"),
            participant("solvent", "THF", "solvent", ev, {"volume": ("2 mL", q(2, "mL"))}),
        ]
    elif entry == 7:
        ps = [
            participant("substrate", "1a", "reactant", ev,
                        {"mass": ("44.0 mg", q(44.0, "mg")), "amount": ("0.234 mmol", q(0.234, "mmol")),
                         "equivalents": ("reference substrate", q(1, "equiv"))}, ref="1a"),
            participant("thiol", "2", "reactant", ev,
                        {"volume": ("25.5 uL", q(25.5, "uL")), "amount": ("0.280 mmol", q(0.280, "mmol")),
                         "equivalents": ("1.2 mol eq", q(1.2, "equiv"))}, ref="1a"),
            participant("MgBr2", "MgBr2", "reagent", ev,
                        {"mass": ("51.6 mg", q(51.6, "mg")), "amount": ("0.280 mmol", q(0.280, "mmol")),
                         "equivalents": ("1.2 mol eq", q(1.2, "equiv"))}, ref="1a"),
            participant("Et3N", "Et3N", "reagent", ev,
                        {"volume": ("39.2 uL", q(39.2, "uL")), "amount": ("0.280 mmol", q(0.280, "mmol")),
                         "equivalents": ("1.2 mol eq", q(1.2, "equiv"))}, ref="1a"),
            participant("solvent", "THF", "solvent", ev, {"volume": ("2 mL", q(2, "mL"))}),
        ]
    else:
        ps = [
            participant("substrate", "1a", "reactant", ev),
            participant("thiol", "2", "reactant", ev,
                        {"equivalents": (f"{x:g} mol eq", q(x, "equiv"))}, ref="1a"),
        ]
        for i, a in enumerate(agents):
            role = "additive" if a == "LiCl" else "reagent"
            ps.append(participant(f"agent{i+1}", a, role, ev,
                                  {"equivalents": (f"{x:g} mol eq", q(x, "equiv"))}, ref="1a"))
        ps.append(participant("solvent", solvent, "solvent", ev))
    ps += [participant("product3", "3a", "product", ev, product_role="2,3,5-regioisomer"),
           participant("product4", "4a", "product", ev, product_role="2,3,4-regioisomer")]
    if entry in (2, 7):
        ps += [participant("acid", "HCl", "workup_agent", [P4, IMG4],
                           {"volume": ("3 mL of 1 N HCl aq", q(3, "mL"))}, stage="workup"),
               participant("extractant", "CHCl3", "workup_agent", [P4, IMG4],
                           {"volume": ("20 mL x 3", q(20, "mL"))}, stage="workup"),
               participant("drying", "MgSO4", "workup_agent", [P4, IMG4], stage="workup"),
               participant("stationary", "silica", "purification_agent", [P4, IMG4], stage="purification"),
               participant("eluent1", "hexane", "purification_agent", [P4, IMG4], stage="purification"),
               participant("eluent2", "EtOAc", "purification_agent", [P4, IMG4], stage="purification")]
    return ps


for entry, agents, x, solvent, yld, ratio, ratio_basis in table1:
    ev = [P1, P2, IMG1, IMG2] + ([P4, IMG4] if entry in (2, 7) else [])
    unresolved = ["No explicitly associated SI was supplied."]
    if entry not in (2, 7):
        unresolved.append("Only table-level conditions are available; charge masses/volumes, atmosphere, workup and purification are unreported.")
    r = base_reaction(f"table1_entry{entry}", ev, "specific_experiment", ["optimization"],
                      {"page": 2, "region": "Table 1", "entry": entry},
                      "Tandem thia-Michael/Dieckmann condensation",
                      "regioselective cyclization to trisubstituted thiophene(s)", unresolved)
    r["general_scheme_ids"] = [rid("scheme1_general")]
    r["participants"] = detailed_participants(entry, ev, x, solvent, agents)
    input_pids = [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction" and p["role"] != "product"]
    r["steps"] = [step("reaction", ev, "reaction", input_pids,
                       "stir allenyl ester 1a, methyl mercaptoacetate and listed agent(s) in solvent",
                       ("rt", "r.t."), ("3 h", q(3, "h")), "argon" if entry in (2, 7) else None)]
    r["solvent_conditions"]["composition"] = F(solvent if entry not in (2, 7) else "anhydrous THF",
                                                 [mid(solvent)], ev)
    if entry in (2, 7):
        r["solvent_conditions"]["total_volume"] = F("2 mL", q(2, "mL"), [P4, IMG4])
        scale = 0.286 if entry == 2 else 0.234
        r["reaction_scale"] = F(f"1a, {scale:.3f} mmol",
                                 {"reference_molecule_id": mid("1a"), "amount": q(scale, "mmol")}, [P4, IMG4])
        ratio_eluent = "8:1" if entry == 2 else "9:1"
        r["steps"] += [
            step("acid_treatment", [P4, IMG4], "workup", ["acid"], "treat with 1 N aqueous HCl"),
            step("extraction", [P4, IMG4], "workup", ["extractant"], "extract with CHCl3 three times"),
            step("dry_filter_concentrate", [P4, IMG4], "workup", ["drying"],
                 "dry over anhydrous MgSO4, filter, and concentrate in vacuo"),
            step("chromatography", [P4, IMG4], "purification", ["stationary", "eluent1", "eluent2"],
                 f"column chromatography on Silica Gel PSQ 60B with n-hexane-AcOEt {ratio_eluent}"),
        ]
        if entry == 2:
            r["steps"].append(step("analytical_ptlc", [P4, IMG4], "purification", ["product3"],
                                   "further PTLC purification to obtain an analytical sample of 3a"))
            paragraph = "Typical procedure [18]: 1a (53.8 mg, 0.286 mmol), 2 (26.0 uL, 0.286 mmol), DBU (42.7 uL, 0.286 mmol), anhydrous THF (2 mL), room temperature under argon, 3 h; HCl workup, CHCl3 extraction, MgSO4 drying, and silica chromatography (8:1)."
        else:
            paragraph = "Typical procedure [23]: MgBr2 (51.6 mg, 0.280 mmol), 1a (44.0 mg, 0.234 mmol), 2 (25.5 uL, 0.280 mmol), Et3N (39.2 uL, 0.280 mmol), anhydrous THF (2 mL), room temperature under argon, 3 h; HCl workup, CHCl3 extraction, MgSO4 drying, and silica chromatography (9:1)."
        r["procedure"]["original_text"] = F(paragraph, paragraph, [P4, IMG4])
        r["procedure"]["pretreatment"] = F("anhydrous THF; under argon", ["anhydrous THF", "argon"], [P4, IMG4])
        if entry == 2:
            r["procedure"]["addition_order"] = F("DBU was added to a solution of 1a and 2 in anhydrous THF.",
                                                   ["solution: 1a + 2 in THF", "add DBU"], [P4, IMG4])
        else:
            r["procedure"]["addition_order"] = F("Et3N was added to a solution of MgBr2, 1a and 2 in anhydrous THF.",
                                                   ["solution: MgBr2 + 1a + 2 in THF", "add Et3N"], [P4, IMG4])
    else:
        r["procedure"]["original_text"] = F(
            {"table": "Table 1", "entry": entry, "agents": agents, "X_mol_eq": x,
             "solvent": solvent, "temperature": "rt", "time": "3 h"},
            "Table-level optimization entry.", [P2, IMG2])
    if entry == 1:
        r["outcomes"] = [
            outcome("yield", "yield", ["3a", "4a"], ev, "0", q(0, "%"),
                    "intended thiophene products", "isolated", "isolated yield",
                    {"span": "1a + 2 to 3a + 4a"}),
            outcome("no_reaction", "no_reaction", [], ev, "no thia-Michael addition occurred", True,
                    scope={"observation": "no thia-Michael addition"}),
        ]
        classify_none(r, ev)
    else:
        yield_targets = ["3a"] if ratio == (100, 0) else ["4a"] if ratio == (0, 100) else ["3a", "4a"]
        yield_raw = f"{yld}%"
        if entry == 2:
            yield_raw = "Table 1: 99%; procedure [18]: 74.0 mg, 99%, 3a/4a = 97:3"
        elif entry == 7:
            yield_raw = "Table 1: 85%; procedure [23]: 4a 52.3 mg, 85%"
        r["outcomes"] = [outcome("yield", "yield", yield_targets, ev, yield_raw, q(yld, "%"),
                                  "reported product(s)", "isolated", "isolated yield",
                                  {"span": "1a + 2 to reported thiophene product(s)"})]
        if entry == 2:
            r["outcomes"].append(outcome("product_amount", "product_amount", ["3a", "4a"], [P4, IMG4],
                                                 "74.0 mg mixture", q(74.0, "mg"), "inseparable mixture"))
        if entry == 7:
            r["outcomes"].append(outcome("product_amount", "product_amount", ["4a"], [P4, IMG4],
                                                 "52.3 mg", q(52.3, "mg"), "sole product"))
        r["outcomes"].append(outcome(
            "regioselectivity", "selectivity", ["3a", "4a"], ev,
            f"3a/4a = {ratio[0]}:{ratio[1]}",
            {"metric": "3a/4a molar ratio", "ordered_molecule_ids": [mid("3a"), mid("4a")],
             "ratio": list(ratio)}, "regioisomer ratio", basis=ratio_basis,
            scope={"metric": "3a/4a", "order": [mid("3a"), mid("4a")]},
        ))
        observed = changes3 if ratio == (100, 0) else changes4 if ratio == (0, 100) else changes3 + changes4
        classify_manual(r, observed, ev,
                        "Candidate net changes cover the observed regioisomer branch(es); aromatic tautomer/bond localization and atom mapping remain partial.")
    reactions.append(r)


# Table 2 substrate scope under DBU.
table2 = {
    "b": (93, (95, 5), "500 MHz 1H NMR (CDCl3)"),
    "c": (91, (98, 2), "400 MHz 1H NMR (CDCl3)"),
    "d": (92, (98, 2), "400 MHz 1H NMR (CDCl3)"),
    "e": (95, (97, 3), "500 MHz 1H NMR (CDCl3)"),
    "f": (92, (98, 2), "400 MHz 1H NMR (CDCl3)"),
    "g": (79, (100, 0), "500 MHz 1H NMR (CDCl3)"),
}
for entry, suffix in enumerate("bcdefg", 1):
    yld, ratio, basis = table2[suffix]
    ev = ([P1, P2, IMG1, IMG2] if suffix in "bcde" else [P2, IMG2])
    r = base_reaction(f"table2_entry{entry}", ev, "specific_experiment", ["substrate_scope"],
                      {"page": 2, "region": "Table 2", "entry": entry},
                      "Tandem thia-Michael/Dieckmann condensation in the presence of DBU",
                      "regioselective cyclization to 2,3,5-trisubstituted thiophene as major product",
                      ["Only table-level conditions are supplied; charge masses/volumes, atmosphere, workup and purification are unreported.",
                       "No explicitly associated SI was supplied."])
    r["general_scheme_ids"] = [rid("scheme1_general")]
    r["participants"] = [
        participant("substrate", f"1{suffix}", "reactant", ev),
        participant("thiol", "2", "reactant", ev, {"equivalents": ("1 mol eq", q(1, "equiv"))}, ref=f"1{suffix}"),
        participant("DBU", "DBU", "reagent", ev, {"equivalents": ("1 mol eq", q(1, "equiv"))}, ref=f"1{suffix}"),
        participant("solvent", "THF", "solvent", ev),
        participant("product3", f"3{suffix}", "product", ev, product_role="major 2,3,5-regioisomer"),
        participant("product4", f"4{suffix}", "product", ev, product_role="minor 2,3,4-regioisomer"),
    ]
    r["steps"] = [step("reaction", ev, "reaction", ["substrate", "thiol", "DBU", "solvent"],
                       "stir allenyl ester, methyl mercaptoacetate and DBU in THF",
                       ("rt", "r.t."), ("3 h", q(3, "h")))]
    r["solvent_conditions"]["composition"] = F("THF", [mid("THF")], ev)
    r["procedure"]["original_text"] = F(
        {"table": "Table 2", "entry": entry, "substrate": f"1{suffix}", "DBU": "1 mol eq",
         "methyl_mercaptoacetate": "1 mol eq", "solvent": "THF", "temperature": "rt", "time": "3 h"},
        "Table-level substrate-scope entry.", [P2, IMG2])
    yield_targets = [f"3{suffix}"] if ratio == (100, 0) else [f"3{suffix}", f"4{suffix}"]
    r["outcomes"] = [
        outcome("yield", "yield", yield_targets, ev,
                f"{yld}% ({' + '.join(yield_targets)})", q(yld, "%"),
                "reported thiophene product(s)", "isolated", "isolated yield",
                {"span": f"1{suffix} + 2 to reported thiophene product(s)"}),
        outcome("regioselectivity", "selectivity", [f"3{suffix}", f"4{suffix}"], ev,
                f"3{suffix}/4{suffix} = {ratio[0]}:{ratio[1]}",
                {"metric": f"3{suffix}/4{suffix} molar ratio",
                 "ordered_molecule_ids": [mid(f"3{suffix}"), mid(f"4{suffix}")], "ratio": list(ratio)},
                "regioisomer ratio", basis=basis,
                scope={"metric": f"3{suffix}/4{suffix}",
                       "order": [mid(f"3{suffix}"), mid(f"4{suffix}")]})]
    classify_manual(r, changes3 if ratio == (100, 0) else changes3 + changes4, ev,
                    "Candidate net changes cover the observed regioisomer branch(es); atom mapping remains partial.")
    reactions.append(r)


# Table 3 substrate scope under MgBr2/Et3N.
table3 = {"b": 80, "c": 92, "d": 80, "e": 89, "f": 86, "g": 77}
for entry, suffix in enumerate("bcdefg", 1):
    yld = table3[suffix]
    ev = [P2, P3, IMG2, IMG3]
    r = base_reaction(f"table3_entry{entry}", ev, "specific_experiment", ["substrate_scope"],
                      {"page": 3, "region": "Table 3", "entry": entry},
                      "Tandem thia-Michael/Dieckmann condensation in the presence of MgBr2/Et3N",
                      "regioselective cyclization to 2,3,4-trisubstituted thiophene",
                      ["Only table-level conditions are supplied; charge masses/volumes, atmosphere, workup and purification are unreported.",
                       "No explicitly associated SI was supplied."])
    r["general_scheme_ids"] = [rid("scheme1_general")]
    r["participants"] = [
        participant("substrate", f"1{suffix}", "reactant", ev),
        participant("thiol", "2", "reactant", ev, {"equivalents": ("1.2 mol eq", q(1.2, "equiv"))}, ref=f"1{suffix}"),
        participant("MgBr2", "MgBr2", "reagent", ev, {"equivalents": ("1.2 mol eq", q(1.2, "equiv"))}, ref=f"1{suffix}"),
        participant("Et3N", "Et3N", "reagent", ev, {"equivalents": ("1.2 mol eq", q(1.2, "equiv"))}, ref=f"1{suffix}"),
        participant("solvent", "THF", "solvent", ev),
        participant("product", f"4{suffix}", "product", ev, product_role="sole 2,3,4-regioisomer"),
    ]
    r["steps"] = [step("reaction", ev, "reaction", ["substrate", "thiol", "MgBr2", "Et3N", "solvent"],
                       "stir allenyl ester, methyl mercaptoacetate, MgBr2 and Et3N in THF",
                       ("rt", "r.t."), ("3 h", q(3, "h")))]
    r["solvent_conditions"]["composition"] = F("THF", [mid("THF")], ev)
    r["procedure"]["original_text"] = F(
        {"table": "Table 3", "entry": entry, "substrate": f"1{suffix}",
         "MgBr2": "1.2 mol eq", "Et3N": "1.2 mol eq", "methyl_mercaptoacetate": "1.2 mol eq",
         "solvent": "THF", "temperature": "rt", "time": "3 h"},
        "Table-level substrate-scope entry.", [P3, IMG3])
    r["outcomes"] = [outcome("yield", "yield", [f"4{suffix}"], ev,
                                      f"{yld}% (4{suffix})", q(yld, "%"),
                                      "sole 2,3,4-regioisomer", "isolated", "isolated yield",
                                      {"span": f"1{suffix} + 2 to 4{suffix}"})]
    classify_manual(r, changes4, ev,
                    "Candidate net changes describe the sole displayed 2,3,4-regioisomer branch; atom mapping remains partial.")
    reactions.append(r)


# Scheme 2: low-temperature mechanistic experiment.
ev = [P2, P3, IMG2, IMG3]
r = base_reaction("scheme2_low_temperature", ev, "specific_experiment", ["mechanistic_study"],
                  {"page": 3, "region": "Scheme 2"},
                  "Reaction of allenyl ester 1a in the presence of DBU at -78 °C",
                  "low-temperature thia-Michael addition with minor tandem cyclization",
                  ["Yield measurement basis for 4a and 5a is not reported in the supplied main text.",
                   "Detection method/limit for '3a was not observed' is unreported.",
                   "No explicitly associated SI was supplied."])
r["general_scheme_ids"] = [rid("scheme1_general")]
r["participants"] = [
    participant("substrate", "1a", "reactant", ev),
    participant("thiol", "2", "reactant", ev, {"equivalents": ("1 mol eq", q(1, "equiv"))}, ref="1a"),
    participant("DBU", "DBU", "reagent", ev, {"equivalents": ("1 mol eq", q(1, "equiv"))}, ref="1a"),
    participant("solvent", "THF", "solvent", ev),
    participant("product4", "4a", "product", ev, product_role="minor cyclized product"),
    participant("product5", "5a", "product", ev, product_role="major thia-Michael adduct"),
    participant("unobserved3", "3a", "product", ev, product_role="not observed"),
]
r["steps"] = [step("reaction", ev, "reaction", ["substrate", "thiol", "DBU", "solvent"],
                   "stir 1a, 2 and DBU in THF at low temperature",
                   ("-78 °C", q(-78, "°C")), ("3 h", q(3, "h")))]
r["solvent_conditions"]["composition"] = F("THF", [mid("THF")], ev)
r["procedure"]["original_text"] = F(
    "Scheme 2: 1a + 2 (1 mol eq), DBU (1 mol eq), THF, -78 °C, 3 h -> 4a 6% + 5a 91%; 3a was not observed.",
    "Low-temperature mechanistic experiment.", ev)
r["outcomes"] = [
    outcome("yield_4a", "yield", ["4a"], ev, "4a 6%", q(6, "%"), "minor product",
            scope={"product": mid("4a")}),
    outcome("yield_5a", "yield", ["5a"], ev, "5a 91%", q(91, "%"), "major product",
            scope={"product": mid("5a")}),
    outcome("not_observed_3a", "other", ["3a"], ev, "3a was not observed", "not observed",
            "unobserved product", scope={"product": mid("3a")}),
]
scheme2_changes = [dict(x, change_id="s2_" + x["change_id"], branch="4a") for x in changes4]
scheme2_changes += [dict(x, change_id="s2_" + x["change_id"], branch="5a") for x in changes5]
classify_manual(r, scheme2_changes, ev,
                "Union of the observed 4a and 5a branches; branch-resolved complete atom mapping is unavailable.")
reactions.append(r)


# Scheme 3: two distinct experiments from the common isolated intermediate 5a.
ev = [P2, P3, IMG2, IMG3]
r = base_reaction("scheme3_dbu_cyclization", ev, "specific_experiment", ["mechanistic_study"],
                  {"page": 3, "region": "Scheme 3", "branch": "DBU"},
                  "Cyclization of thia-Michael adduct 5a with DBU",
                  "base-promoted cyclization to 2,3,5-trisubstituted thiophene",
                  ["Yield and regioisomer-ratio measurement bases are unreported in the supplied main text.",
                   "No explicitly associated SI was supplied."])
r["general_scheme_ids"] = [rid("scheme1_general")]
r["participants"] = [
    participant("substrate", "5a", "reactant", ev),
    participant("DBU", "DBU", "reagent", ev, {"equivalents": ("1 mol eq", q(1, "equiv"))}, ref="5a"),
    participant("solvent", "THF", "solvent", ev),
    participant("product3", "3a", "product", ev, product_role="major product"),
    participant("product4", "4a", "product", ev, product_role="minor product"),
]
r["steps"] = [step("reaction", ev, "reaction", ["substrate", "DBU", "solvent"],
                   "stir 5a with DBU in THF", ("rt", "r.t."), ("3 h", q(3, "h")))]
r["solvent_conditions"]["composition"] = F("THF", [mid("THF")], ev)
r["procedure"]["original_text"] = F(
    "Scheme 3 upper branch: 5a, DBU (1 mol eq), THF, rt, 3 h -> 99% (3a/4a = >99:<1).",
    "DBU-promoted cyclization experiment.", ev)
r["outcomes"] = [
    outcome("yield", "yield", ["3a", "4a"], ev, "99%", q(99, "%"), "product mixture",
            scope={"span": "5a to 3a + 4a"}),
    outcome("regioselectivity", "selectivity", ["3a", "4a"], ev, "3a/4a = >99:<1",
            {"metric": "3a/4a ratio", "ordered_molecule_ids": [mid("3a"), mid("4a")],
             "first": {"number": 99, "qualifier": "gt"}, "second": {"number": 1, "qualifier": "lt"}},
            "regioisomer ratio", scope={"metric": "3a/4a", "order": [mid("3a"), mid("4a")]}),
]
cyc_union = [dict(x, change_id="s3_dbu_" + x["change_id"]) for x in changes_cyclize3]
cyc_union += [dict(x, change_id="s3_minor_" + x["change_id"], branch="minor 4a") for x in changes_cyclize4]
classify_manual(r, cyc_union, ev,
                "Major 3a cyclization and source-reported <1 part 4a are represented as a union; mapping remains partial.")
reactions.append(r)

r = base_reaction("scheme3_mgbr2_cyclization", ev, "specific_experiment", ["mechanistic_study"],
                  {"page": 3, "region": "Scheme 3", "branch": "MgBr2/Et3N"},
                  "Cyclization of thia-Michael adduct 5a with MgBr2/Et3N",
                  "Lewis-acid/base-promoted partial cyclization to 2,3,5-trisubstituted thiophene",
                  ["Yield and recovery measurement bases are unreported in the supplied main text.",
                   "No explicitly associated SI was supplied."])
r["general_scheme_ids"] = [rid("scheme1_general")]
r["participants"] = [
    participant("substrate", "5a", "reactant", ev),
    participant("MgBr2", "MgBr2", "reagent", ev, {"equivalents": ("1.2 mol eq", q(1.2, "equiv"))}, ref="5a"),
    participant("Et3N", "Et3N", "reagent", ev, {"equivalents": ("1.2 mol eq", q(1.2, "equiv"))}, ref="5a"),
    participant("solvent", "THF", "solvent", ev),
    participant("product", "3a", "product", ev, product_role="product"),
]
r["steps"] = [step("reaction", ev, "reaction", ["substrate", "MgBr2", "Et3N", "solvent"],
                   "stir 5a with MgBr2 and Et3N in THF", ("rt", "r.t."), ("3 h", q(3, "h")))]
r["solvent_conditions"]["composition"] = F("THF", [mid("THF")], ev)
r["procedure"]["original_text"] = F(
    "Scheme 3 lower branch: 5a, MgBr2 (1.2 mol eq), Et3N (1.2 mol eq), THF, rt, 3 h -> 3a 15%; recovery of 5a 77%.",
    "MgBr2/Et3N cyclization experiment.", ev)
r["outcomes"] = [
    outcome("yield", "yield", ["3a"], ev, "3a 15%", q(15, "%"), "product",
            scope={"span": "5a to 3a"}),
    outcome("recovery", "recovery", ["5a"], ev, "Recovery of 5a 77%", q(77, "%"),
            "recovered starting material", scope={"starting_material": mid("5a")}),
]
classify_manual(r, changes_cyclize3, ev,
                "Candidate net changes describe the observed 5a-to-3a branch; recovered 5a is not treated as a product yield.")
reactions.append(r)


# Scheme 4 is an author-proposed pathway, not an observed experiment.
ev = [P3, IMG3]
r = base_reaction("scheme4_proposed_mechanism", ev, "proposed_mechanism", [],
                  {"page": 3, "region": "Scheme 4"},
                  "Plausible reaction pathways for regioselective thiophene formation",
                  "proposed enolate-controlled regioselective tandem pathway",
                  ["The mechanism is explicitly described as plausible and not fully understood.",
                   "Intermediates A, B and C are proposed, not isolated experimental records."])
r["general_scheme_ids"] = [rid("scheme1_general")]
r["participants"] = [
    participant("substrate", "g1", "reactant", ev), participant("thiol", "2", "reactant", ev),
    participant("DBU", "DBU", "reagent", ev), participant("MgBr2", "MgBr2", "reagent", ev),
    participant("Et3N", "Et3N", "reagent", ev), participant("A", "A", "product", ev, product_role="proposed intermediate"),
    participant("B", "B", "product", ev, product_role="proposed intermediate"),
    participant("C", "C", "product", ev, product_role="proposed intermediate"),
    participant("adduct", "g5", "product", ev, product_role="proposed/observed class of adduct"),
    participant("product3", "g3", "product", ev, product_role="DBU branch product class"),
    participant("product4", "g4", "product", ev, product_role="MgBr2/Et3N branch product class"),
]
r["procedure"]["original_text"] = F(
    "Scheme 4: plausible pathways via kinetic enolate A, thia-Michael adduct 5, thermodynamic enolate B, or magnesium enolate C.",
    "Author-proposed mechanism; not a concrete experiment.", ev)
r["notes"] = [F(
    "Although the reaction mechanism underlying the tandem thia-Michael/Dieckmann condensation is not fully understood, a plausible pathway is shown in Scheme 4.",
    "Author explicitly qualifies the pathway as plausible.", ev)]
classify_unresolved(r, ev, "Proposed mechanism has no single observed atom-mapped experiment; no taxonomy labels are asserted.")
reactions.append(r)

bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": f"{DOC}:bibliography", "document_id": DOC, "evidence_ids": [P1, IMG1],
     "kind": "bibliographic_metadata",
     "title": "Tandem thia-Michael/Dieckmann condensation of allenyl esters for the regioselective synthesis of trisubstituted thiophenes",
     "authors": ["Michiyasu Nakao", "Munehisa Toguchi", "Yuki Shimabukuro", "Shigeki Sano"],
     "journal": "Tetrahedron Letters", "volume": 61, "year": 2020, "article": "152271",
     "doi": "10.1016/j.tetlet.2020.152271"},
    {"context_id": f"{DOC}:spectroscopy_3a", "document_id": DOC, "evidence_ids": [P4, IMG4],
     "kind": "characterization", "molecule_id": mid("3a"),
     "note": "Page 4 reports color, 1H/13C NMR, IR and HRMS for 3a; retained as source context, not a synthetic experiment."},
    {"context_id": f"{DOC}:spectroscopy_4a", "document_id": DOC, "evidence_ids": [P4, IMG4],
     "kind": "characterization", "molecule_id": mid("4a"),
     "note": "Page 4 reports color, 1H/13C NMR, IR, HRMS and elemental analysis for 4a; retained as source context, not a synthetic experiment."},
    {"context_id": f"{DOC}:si_statement", "document_id": DOC, "evidence_ids": [P4, IMG4],
     "kind": "source_boundary",
     "note": "The main text states that supplementary data are online, but no SI file was explicitly supplied or associated in this workspace; SI was not used."},
    {"context_id": f"{DOC}:references", "document_id": DOC, "evidence_ids": [P4, IMG4],
     "kind": "references",
     "note": "Cited literature entries are not experiments of the current paper and were not converted into reaction records."},
]

raw = RUN / "raw_model" / "pages_0001_0004_extraction.json"
submitted = RUN / "submitted.json"
st.write(raw, bundle)
st.write(submitted, bundle)
print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
                  "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in bundle["reactions"]),
                  "raw": str(raw), "submitted": str(submitted)}, ensure_ascii=False, indent=2))
