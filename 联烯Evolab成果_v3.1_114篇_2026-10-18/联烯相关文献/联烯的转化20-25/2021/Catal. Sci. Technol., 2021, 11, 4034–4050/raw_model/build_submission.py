#!/usr/bin/env python3
"""Source-specific v3.1 builder for doc_837adcd3d5d8cb7b, PDF pages 1-6.

This script contains only source-derived candidates and explicit extraction decisions.
It does not alter the schema, tools, or quality rules.
"""
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
SKILL = ROOT / "skills" / "chemical-reaction-extraction"
sys.path.insert(0, str(SKILL / "scripts"))
from standard_records import (  # noqa: E402
    fact, new_molecule, new_reaction, new_participant, new_step,
    new_outcome, stereo_inspect, source_mentions, write, bond_label,
)

DOC = "doc_837adcd3d5d8cb7b"
OUT = ROOT / "runs" / "current"
RAW = OUT / "raw_model"


def et(page):
    return f"{DOC}:page{page}:text"


def ei(page):
    return f"{DOC}:page{page}:image"


def ev(*pages, text=True, image=True):
    result = []
    for page in pages:
        if text:
            result.append(et(page))
        if image and page <= 6:
            result.append(ei(page))
    return result


def f(raw=None, norm=None, evidence=(), method="original_extraction",
      status=None, derivation=None):
    verification = "inferred" if method in ("inference", "structure_derivation") and (raw is not None or norm is not None) else None
    return fact(raw=raw, normalized=norm, evidence=evidence, method=method,
                status=status, verification=verification, derivation=derivation)


def q(number, unit, qualifier=None, raw_text=None):
    value = {"number": number, "unit": unit}
    if qualifier:
        value["qualifier"] = qualifier
    if raw_text:
        value["raw_text"] = raw_text
    return value


def missing(evidence=()):
    return fact(evidence=evidence)


def set_status(m, key, value, evidence, raw=None, derivation=None, method="inference"):
    m["stereochemistry"][key] = f(raw if raw is not None else value, value, evidence, method=method,
                                     derivation=derivation or "Status retained conservatively from the source drawing/name; absence of a descriptor is not treated as achiral or racemic.")


def molecule(slug, name, label, evidence, smiles=None, abbrev=None,
             chirality="unresolved", geometry="unresolved", rep="not_encoded",
             axial_raw=None, identity_note=None):
    mid = f"{DOC}:mol_{slug}"
    m = new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    m["names"]["original"] = f(name, name, evidence)
    m["names"]["english"] = f(name, name, evidence)
    if abbrev:
        m["names"]["abbreviations"] = f(abbrev, [abbrev], evidence)
    m["cas_numbers"] = fact(evidence=evidence, status="unreported")
    m["identity_note"] = f(identity_note, identity_note, evidence) if identity_note else fact(evidence=evidence, status="unreported")
    if smiles:
        diag = stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = f(smiles, diag["canonical_isomeric_smiles"], evidence,
                          method="structure_derivation",
                          derivation="Candidate graph transcribed from the source name/formula or a simple source drawing; parsing is not independent identity verification.")
        m["stereochemistry"]["isomeric_smiles"] = f(smiles, diag["canonical_isomeric_smiles"], evidence,
                                                       method="structure_derivation",
                                                       derivation="Canonical isomeric SMILES returned by the standard stereo graph routine.")
        m["stereochemistry"]["identity_key"] = f(diag["stereo_graph_key"], diag["stereo_graph_key"], evidence,
                                                  method="structure_derivation",
                                                  derivation="Stereo graph key computed by the standard stereo routine; it is not a sample identity.")
        m["stereochemistry"]["atom_centers"] = f(diag["atom_centers"], diag["atom_centers"], evidence,
                                                   method="structure_derivation", derivation="Standard stereo graph diagnostics.")
        m["stereochemistry"]["double_bonds"] = f(diag["double_bonds"], diag["double_bonds"], evidence,
                                                   method="structure_derivation", derivation="Standard stereo graph diagnostics.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = fact(evidence=evidence, status="unreported")
    if axial_raw:
        m["stereochemistry"]["axial"] = f(axial_raw, axial_raw, evidence)
    set_status(m, "chirality_status", chirality, evidence)
    set_status(m, "geometry_status", geometry, evidence)
    set_status(m, "representation_status", rep, evidence)
    return m


M = {}


def add_mol(*args, **kwargs):
    m = molecule(*args, **kwargs)
    M[m["molecule_id"].split(":mol_")[-1]] = m
    return m


# Registry: exact organozirconium/borazane graphs are deliberately not guessed.
add_mol("complex1", "Cp2Zr(Cl)(mu-Me3SiC3SiMe3)Zr(Cl)Cp2", "1", ev(2,3,4,5,6), abbrev="1",
        chirality="unspecified", geometry="unresolved", rep="not_encoded",
        axial_raw="allenediide-bridged dinuclear complex; no axial descriptor reported",
        identity_note="Dinuclear zirconocene allenediide-bridged precatalyst; complete organometallic connectivity is retained by source drawing rather than forced SMILES.")
add_mol("complex2", "dinuclear zirconocene dimethyl complex", "2", ev(2,5,6), abbrev="2",
        chirality="unspecified", geometry="unresolved", rep="not_encoded")
add_mol("complex4", "Cp2Zr(NMe2)(mu-Me3SiC3SiMe3)Zr(NMe2)Cp2", "4", ev(4,5,13), abbrev="4",
        chirality="unspecified", geometry="unresolved", rep="not_encoded")
add_mol("complex5", "dinuclear zirconocene borate complex", "5", ev(4,5,14), abbrev="5",
        chirality="unspecified", geometry="unresolved", rep="not_encoded")
add_mol("complex6", "singly methylated zirconocene complex", "6", ev(5,6,14), abbrev="6",
        chirality="unspecified", geometry="mixture", rep="not_encoded",
        identity_note="The solid contains three Cl/Me orientation isomers; no absolute stereochemical descriptor is assigned.")
add_mol("me_li", "methyl lithium", "MeLi", ev(2,3,5,6,15), smiles="[Li+].[CH3-]", abbrev="MeLi",
        chirality="achiral", geometry="not_applicable", rep="encoded")
add_mol("methylamine_borane", "methylamine borane", "H3B·NMeH2", ev(2,3), smiles="[BH3-][NH2+]C",
        chirality="achiral", geometry="not_applicable", rep="encoded")
add_mol("polyaminoborane", "poly(methylaminoborane)", "[HBNHMe]n", ev(2,3),
        chirality="unresolved", geometry="unresolved", rep="not_encoded")
add_mol("dimethylamine_borane", "dimethylamine borane", "H3B·NMe2H", ev(1,2,3,4,5,6,15),
        smiles="[BH3-][NH+](C)C", abbrev="DMAB", chirality="achiral", geometry="not_applicable", rep="encoded")
add_mol("hydrogen", "hydrogen", "H2", ev(1,2,3,4,6,15), smiles="[H][H]",
        chirality="achiral", geometry="not_applicable", rep="encoded")
add_mol("methane", "methane", "CH4", ev(2), smiles="C", chirality="achiral", geometry="not_applicable", rep="encoded")
add_mol("cyclic_diborazane", "cyclic diborazane", "(H2BNMe2)2", ev(1,3,4),
        chirality="unspecified", geometry="unresolved", rep="not_encoded")
add_mol("linear_diborazane", "linear diborazane", "H3B·NMe2BH2·NMe2H", ev(1,3,4),
        chirality="unspecified", geometry="unresolved", rep="not_encoded")
add_mol("diaminoborane", "diaminoborane", "HB(NMe2)2", ev(3,4),
        chirality="unspecified", geometry="unresolved", rep="not_encoded")
add_mol("aminoborane", "dimethylaminoborane", "H2B=NMe2", ev(1,3,4),
        chirality="unspecified", geometry="unspecified", rep="not_encoded")
add_mol("mu_aminodiborane", "mu-aminodiborane", "BH2(mu-Me2N)(mu-H)BH2", ev(4),
        chirality="unspecified", geometry="unresolved", rep="not_encoded")
add_mol("dimethylamine", "dimethylamine", "NMe2H", ev(4), smiles="CNC", abbrev="NMe2H",
        chirality="achiral", geometry="not_applicable", rep="encoded")
add_mol("bh3_thf", "borane-tetrahydrofuran complex", "BH3·THF", ev(4), abbrev="BH3·THF",
        chirality="unresolved", geometry="unresolved", rep="not_encoded")
add_mol("li_nme2", "lithium dimethylamide", "LiNMe2", ev(4,13), smiles="[Li+].[N-](C)C", abbrev="LiNMe2",
        chirality="achiral", geometry="not_applicable", rep="encoded")
add_mol("nabh4", "sodium borohydride", "NaBH4", ev(4,5,14), smiles="[Na+].[BH4-]", abbrev="NaBH4",
        chirality="achiral", geometry="not_applicable", rep="encoded")
add_mol("licl", "lithium chloride", "LiCl", ev(4,5), smiles="[Li+].[Cl-]",
        chirality="achiral", geometry="not_applicable", rep="encoded")
add_mol("nacl", "sodium chloride", "NaCl", ev(4,14), smiles="[Na+].[Cl-]",
        chirality="achiral", geometry="not_applicable", rep="encoded")
add_mol("d3b_nme2h", "B-deuterated dimethylamine borane", "D3B·NMe2H", ev(4),
        smiles="[2H][B-]([2H])([2H])[NH+](C)C", chirality="achiral", geometry="not_applicable", rep="encoded")
add_mol("h3b_nme2d", "N-deuterated dimethylamine borane", "H3B·NMe2D", ev(4),
        smiles="[BH3-][N+]([2H])(C)C", chirality="achiral", geometry="not_applicable", rep="encoded")
add_mol("d3b_nme2d", "B,N-deuterated dimethylamine borane", "D3B·NMe2D", ev(4),
        smiles="[2H][B-]([2H])([2H])[N+]([2H])(C)C", chirality="achiral", geometry="not_applicable", rep="encoded")
add_mol("li_amidoborane", "lithium dimethylamidoborane", "Li[NMe2BH3]", ev(6),
        smiles="[Li+].[BH3-][N](C)C", abbrev="Li[NMe2BH3]", chirality="achiral", geometry="not_applicable", rep="encoded")
add_mol("toluene", "toluene", "toluene", ev(3,4,5,6,15), smiles="Cc1ccccc1",
        chirality="achiral", geometry="not_applicable", rep="encoded")
add_mol("toluene_d8", "toluene-d8", "toluene-d8", ev(3,5),
        chirality="achiral", geometry="not_applicable", rep="not_encoded")
add_mol("benzene", "benzene", "benzene", ev(5,6,13,14), smiles="c1ccccc1",
        chirality="achiral", geometry="not_applicable", rep="encoded")
add_mol("diethyl_ether", "diethyl ether", "Et2O", ev(4,6,13,14,15), smiles="CCOCC", abbrev="Et2O",
        chirality="achiral", geometry="not_applicable", rep="encoded")


def pid(rid, index):
    return rid.split(":")[-1] + f"_p{index}"


def participant(rid, index, slug, role, evidence, amounts=None, reference=None, product_role=None, stage="reaction"):
    p = new_participant(pid(rid, index), M[slug]["molecule_id"], role, evidence, stage=stage)
    if product_role:
        p["product_role"] = f(product_role, product_role, evidence)
    for key, value in (amounts or {}).items():
        if isinstance(value, dict) and "_fact" in value:
            p["charges"][key] = value["_fact"]
        else:
            p["charges"][key] = f(value, value, evidence)
    if reference:
        p["equivalence_reference"] = f(reference, M[reference]["molecule_id"], evidence)
    return p


def step(rid, index, evidence, participants, operation, temperature=None, time=None, atmosphere=None, stage="reaction"):
    s = new_step(rid.split(":")[-1] + f"_s{index}", evidence, stage=stage)
    s["participant_ids"] = participants
    s["operation"] = f(operation, operation, evidence)
    if temperature is not None:
        s["temperature"] = f(temperature, temperature, evidence)
    if time is not None:
        s["time"] = f(time, time, evidence)
    if atmosphere is not None:
        s["atmosphere"] = f(atmosphere, atmosphere, evidence)
    return s


def outcome(rid, index, kind, targets, evidence, value, basis=None, role=None, scope=None, yield_type=None):
    o = new_outcome(rid.split(":")[-1] + f"_o{index}", kind, [M[t]["molecule_id"] for t in targets], evidence)
    o["value"] = f(value, value, evidence)
    if basis is not None:
        o["measurement_basis"] = f(basis, basis, evidence)
    if role is not None:
        o["product_role"] = f(role, role, evidence)
    if scope is not None:
        o["scope"] = f(scope, scope, evidence)
    if yield_type is not None:
        o["yield_type"] = f(yield_type, yield_type, evidence)
    return o


def set_classification(r, evidence, changes=None, coverage="unresolved", limitation=None):
    c = r["classification"]
    if coverage == "not_applicable":
        c["coverage"] = f("not applicable", "not_applicable", evidence, method="inference", derivation="No observed chemical transformation is classified.")
        c["labels"] = f([], [], evidence, method="inference", derivation="Negative/control or physical-characterization record; no intended bond formation is asserted.")
        c["bond_changes"] = []
        c["limitations"] = [limitation or "No observed transformation to classify."]
        return
    if not changes:
        c["coverage"] = f("unresolved", "unresolved", evidence, method="inference", derivation="Source pages do not support safe atom-level mapping.")
        c["labels"] = fact(evidence=evidence, status="unreported")
        c["bond_changes"] = []
        c["limitations"] = [limitation or "Exact atom mapping and bond changes are unresolved."]
        return
    labels = {}
    c["bond_changes"] = []
    for idx, (a, b, before, after, change_type, note) in enumerate(changes, 1):
        cid = f"b{idx}"
        ch = {"change_id": cid, "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
              "before_order": before, "after_order": after, "change_type": change_type,
              "source_note": note}
        c["bond_changes"].append(f(note, ch, evidence, method="inference",
                                           derivation="Class-level net change inferred from the displayed overall transformation; atom identities are not mapped."))
        if change_type in ("formed", "hydrogen_gain", "order_increased"):
            lab = bond_label(a, b, after)
            labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(cid)
    c["labels"] = f("class-level inferred labels", list(labels.values()), evidence, method="inference",
                     derivation="All and only formed/promoted class-level changes above are indexed.")
    c["coverage"] = f("partial", "partial", evidence, method="inference",
                       derivation="Connectivity is class-level only; no complete atom mapping is asserted.")
    c["limitations"] = [limitation or "Exact atom identities, all hydrogen transfers, coproducts, and complete atom mapping remain unresolved."]


R = []


def reaction(slug, record_type, evidence, categories, name, transformation, parts, steps, outcomes,
             general=(), class_changes=None, class_coverage=None, notes=(), unresolved=()):
    rid = f"{DOC}:{slug}"
    r = new_reaction(rid, DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, evidence)
    r["location"] = f({"pages": sorted({int(x.split(":page")[1].split(":")[0]) for x in evidence if ":page" in x}), "source": slug},
                        {"pages": sorted({int(x.split(":page")[1].split(":")[0]) for x in evidence if ":page" in x}), "source": slug}, evidence)
    r["general_scheme_ids"] = [f"{DOC}:{g}" for g in general]
    r["reaction_name"] = f(name, name, evidence)
    r["transformation_type"] = f(transformation, transformation, evidence)
    r["participants"] = parts
    r["steps"] = steps
    r["outcomes"] = outcomes
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": list(unresolved) or ["Machine extraction requires expert review."]}
    r["notes"] = [f(n, n, evidence) for n in notes]
    if class_coverage:
        set_classification(r, evidence, class_changes, class_coverage)
    else:
        set_classification(r, evidence, class_changes)
    R.append(r)
    return r


DEHYDRO = [("H", "H", 0, 1, "formed", "Net H-H formation is reported by H2 evolution."),
           ("B", "N", 0, 1, "formed", "Net B-N coupling is represented by the observed diborazane products.")]
ZR_C_ONE = [("C", "Zr", 0, 1, "formed", "One Zr-C bond is formed on methylation of one chlorido site."),
            ("Cl", "Zr", 1, 0, "broken", "One Zr-Cl bond is displaced with LiCl formation.")]
ZR_C_TWO = ZR_C_ONE + [("C", "Zr", 0, 1, "formed", "A second Zr-C bond is formed on methylation."),
                       ("Cl", "Zr", 1, 0, "broken", "The second Zr-Cl bond is displaced.")]


# Page 2: cited prior work and current generic scheme.
e2 = ev(2,3)
rid = f"{DOC}:literature_dehydropolymerisation"
parts = [participant(rid, 1, "methylamine_borane", "reactant", e2), participant(rid, 2, "complex1", "catalyst", e2),
         participant(rid, 3, "me_li", "reagent", e2), participant(rid, 4, "polyaminoborane", "product", e2, product_role="main product"),
         participant(rid, 5, "hydrogen", "product", e2, product_role="gas coproduct")]
reaction("literature_dehydropolymerisation", "literature_example", e2, [], "previous dehydropolymerisation of methylamine borane",
         "dehydropolymerisation", parts, [], [], class_changes=DEHYDRO,
         notes=["Scheme 1a is explicitly labelled previous work; proposed intermediates are not current-paper experiments."],
         unresolved=["Prior-work conditions and product molecular-weight distribution are not reported on pages 1-6."])

rid = f"{DOC}:scheme1_current_dehydrocoupling"
parts = [participant(rid, 1, "dimethylamine_borane", "reactant", e2), participant(rid, 2, "complex1", "catalyst", e2),
         participant(rid, 3, "me_li", "reagent", e2), participant(rid, 4, "cyclic_diborazane", "product", e2, product_role="main product"),
         participant(rid, 5, "linear_diborazane", "product", e2, product_role="minor product"),
         participant(rid, 6, "aminoborane", "product", e2, product_role="minor product"),
         participant(rid, 7, "diaminoborane", "product", e2, product_role="minor product"),
         participant(rid, 8, "hydrogen", "product", e2, product_role="gas coproduct")]
reaction("scheme1_current_dehydrocoupling", "generic_scheme", e2, [], "catalytic dehydrocoupling of dimethylamine borane",
         "amine-borane dehydrocoupling", parts, [], [], class_changes=DEHYDRO,
         unresolved=["Scheme 1b is a transformation summary, not one uniquely quantified experiment."])

e56 = ev(5,6)
rid = f"{DOC}:scheme_activation_sequence"
parts = [participant(rid, 1, "complex1", "reactant", e56), participant(rid, 2, "me_li", "reagent", e56),
         participant(rid, 3, "complex6", "product", e56, product_role="intermediate"), participant(rid, 4, "complex2", "product", e56, product_role="final methylated complex"),
         participant(rid, 5, "licl", "product", e56, product_role="salt coproduct")]
reaction("scheme_activation_sequence", "generic_scheme", e56, [], "stepwise methylation of complex 1", "precatalyst activation",
         parts, [], [], class_changes=ZR_C_TWO,
         unresolved=["The scheme summarizes a concentration-dependent sequence; it is not a single kinetic experiment."])


def catalytic_parts(rid, evidence, catalyst_slug="complex1", catalyst_loading=2.5, substrate_slug="dimethylamine_borane",
                    me_li=True, extra=()):
    parts = [participant(rid, 1, substrate_slug, "reactant", evidence,
                         amounts={"mass": q(78.4, "mg"), "amount": q(1.33, "mmol")})]
    parts.append(participant(rid, 2, catalyst_slug, "catalyst", evidence,
                             amounts={"mol_percent": q(catalyst_loading, "mol%")}, reference=substrate_slug))
    idx = 3
    if me_li:
        parts.append(participant(rid, idx, "me_li", "reagent", evidence,
                                 amounts={"equivalents": q(2, "equiv")}, reference=catalyst_slug)); idx += 1
    for slug, role, amounts, reference in extra:
        parts.append(participant(rid, idx, slug, role, evidence, amounts=amounts, reference=reference)); idx += 1
    parts += [participant(rid, idx, "toluene", "solvent", evidence, amounts={"volume": q(4, "mL")}),
              participant(rid, idx+1, "hydrogen", "product", evidence, product_role="measured gas product")]
    return parts


# Page 3 Figure 1: six catalyst-loading experiments. General procedure text on page 15 supplies explicit charge/volume/Ar setup.
loading_info = {
    "0_63": (0.63, "similar volumetric trace to the 1.25 mol% run; much slower than 2.0/2.5 mol%"),
    "1_25": (1.25, "similar volumetric trace to the 0.63 mol% run; much slower than 2.0/2.5 mol%"),
    "2_00": (2.00, "among the faster catalyst-loading runs"),
    "2_50": (2.50, "approximately 0.85 equivalents of H2; no full substrate conversion"),
    "3_75": (3.75, "much slower than the 2.0/2.5 mol% runs"),
    "5_00": (5.00, "much slower than the 2.0/2.5 mol% runs"),
}
for key, (loading, observation) in loading_info.items():
    evidence = ev(3) + [et(15)]
    slug = f"catalyst_loading_{key}"
    rid = f"{DOC}:{slug}"
    parts = catalytic_parts(rid, evidence, catalyst_loading=loading)
    pids = [p["participant_id"] for p in parts if p["stage"] == "reaction"]
    st = [step(rid, 1, evidence, pids, "add MeLi and monitor gas evolution", q(25, "°C"), None, "Ar")]
    if key == "2_50":
        outs = [outcome(rid, 1, "product_amount", ["hydrogen"], evidence, q(0.85, "equiv", raw_text="n(H2)/n(H3B·NMe2H) = 0.85"),
                        basis="volumetric gas evolution relative to initial H3B·NMe2H", role="gas product")]
    else:
        outs = [outcome(rid, 1, "other", ["hydrogen"], evidence, observation,
                        basis="volumetric gas-evolution trace", role="gas product")]
    r = reaction(slug, "specific_experiment", evidence, ["optimization", "mechanistic_study"],
                 "catalytic dehydrocoupling catalyst-loading series", "amine-borane dehydrocoupling",
                 parts, st, outs, general=("scheme1_current_dehydrocoupling",), class_changes=DEHYDRO,
                 notes=["Page 15 general procedure gives 1.33 mmol substrate in 4 mL toluene, consistent with page 3's 0.33 M; page 15 was used only as text context and is not claimed as visually covered."],
                 unresolved=["Exact endpoint and time are not textually reported for this individual curve." if key != "2_50" else "The source reports gas equivalents, not a product-isolation yield."])
    r["reaction_scale"] = f("H3B·NMe2H 1.33 mmol", {"reference_molecule_id": M["dimethylamine_borane"]["molecule_id"], "amount": q(1.33, "mmol")}, evidence)
    r["solvent_conditions"]["composition"] = f("toluene", [M["toluene"]["molecule_id"]], evidence)
    r["solvent_conditions"]["total_volume"] = f("4 mL", q(4, "mL"), evidence)
    r["solvent_conditions"]["concentration"] = f("0.33 M", q(0.33, "M"), evidence)
    r["solvent_conditions"]["concentration_basis"] = f("[H3B·NMe2H]", M["dimethylamine_borane"]["molecule_id"], evidence)


# NMR product-distribution run, distinct from the gas-buret loading series because solvent/loading are explicitly different.
e3 = ev(3,4)
rid = f"{DOC}:nmr_product_distribution"
parts = [participant(rid, 1, "dimethylamine_borane", "reactant", e3),
         participant(rid, 2, "complex1", "catalyst", e3, amounts={"mol_percent": q(2.5, "mol%")}, reference="dimethylamine_borane"),
         participant(rid, 3, "me_li", "reagent", e3, amounts={"mol_percent": q(5.2, "mol%")}, reference="dimethylamine_borane"),
         participant(rid, 4, "toluene_d8", "solvent", e3),
         participant(rid, 5, "cyclic_diborazane", "product", e3, product_role="main product"),
         participant(rid, 6, "linear_diborazane", "product", e3, product_role="minor product"),
         participant(rid, 7, "diaminoborane", "product", e3, product_role="minor product"),
         participant(rid, 8, "mu_aminodiborane", "product", e3, product_role="minor product"),
         participant(rid, 9, "hydrogen", "product", e3, product_role="gas coproduct")]
st = [step(rid, 1, e3, [p["participant_id"] for p in parts], "dehydrocoupling followed by 11B NMR analysis", q(297, "K"))]
outs = [outcome(rid, 1, "other", ["cyclic_diborazane"], e3, "main product", basis="11B NMR (96.7 MHz)", role="main product"),
        outcome(rid, 2, "other", ["linear_diborazane", "diaminoborane", "mu_aminodiborane"], e3, "minor components", basis="11B NMR", role="minor products"),
        outcome(rid, 3, "recovery", ["dimethylamine_borane"], e3, "considerable residual starting material", basis="11B NMR", role="residual reactant")]
reaction("nmr_product_distribution", "specific_experiment", e3, ["mechanistic_study"], "11B NMR product-distribution experiment",
         "amine-borane dehydrocoupling", parts, st, outs, general=("scheme1_current_dehydrocoupling",), class_changes=DEHYDRO,
         unresolved=["Individual product yields and the quantity of residual substrate are not reported."])


# Page 4 isotope experiments.
for slug, substrate, observation in [
    ("isotope_d3b_nme2h", "d3b_nme2h", "virtually identical H2 trace with slightly lower initial activity; B-H/D activation may affect slow catalyst formation"),
    ("isotope_h3b_nme2d", "h3b_nme2d", "similar onset followed by significantly decreased activity at higher conversion with linear H2 evolution"),
    ("isotope_d3b_nme2d", "d3b_nme2d", "reduced initial rate and significantly decreased activity at higher conversion"),
]:
    evidence = ev(4)
    rid = f"{DOC}:{slug}"
    parts = [participant(rid, 1, substrate, "reactant", evidence), participant(rid, 2, "complex1", "catalyst", evidence),
             participant(rid, 3, "me_li", "reagent", evidence), participant(rid, 4, "hydrogen", "product", evidence, product_role="gas/isotopic gas outcome")]
    st = [step(rid, 1, evidence, [p["participant_id"] for p in parts], "isotopomer dehydrocoupling and volumetric monitoring")]
    outs = [outcome(rid, 1, "other", ["hydrogen"], evidence, observation, basis="volumetric gas-evolution trace")]
    reaction(slug, "specific_experiment", evidence, ["mechanistic_study"], "isotopomer dehydrocoupling experiment",
             "isotope-effect probe", parts, st, outs, general=("scheme1_current_dehydrocoupling",), class_changes=DEHYDRO,
             notes=["The authors explicitly refrain from assigning a kinetic isotope effect because the system changes kinetic regime."],
             unresolved=["Exact catalyst loading, isotopic gas composition, and numerical initial rate are in unavailable SI figures S5/S6."])


# Page 4 Figure 2 catalyst poisoning and isolated-complex catalytic tests; page 15 supplies the common apparatus/scale.
fig2_defs = [
    ("nme2h_poison", "complex1", True, [("dimethylamine", "additive", {"mol_percent": q(5.0, "mol%")}, "dimethylamine_borane")],
     q(0.4, "equiv", qualifier="approx", raw_text="approximately 0.4 equivalents after 10 hours"), "NMe2H inhibition", False),
    ("bh3_poison", "complex1", True, [("bh3_thf", "additive", {"mol_percent": q(5.0, "mol%")}, "dimethylamine_borane")],
     "almost complete loss of activity", "BH3·THF inhibition", False),
    ("catalysis_complex4", "complex4", False, [], q(1.0, "equiv", qualifier="approx", raw_text="approximately 1.0 equivalent of H2"), "full conversion after a pronounced induction period", True),
    ("catalysis_complex5", "complex5", False, [], "no activity", "inactive borate-complex control", False),
]
for slug, catalyst, use_meli, extra, result, note, full in fig2_defs:
    evidence = ev(4,5) + [et(15)]
    rid = f"{DOC}:{slug}"
    parts = catalytic_parts(rid, evidence, catalyst_slug=catalyst, catalyst_loading=2.5, me_li=use_meli, extra=extra)
    st = [step(rid, 1, evidence, [p["participant_id"] for p in parts], "catalytic dehydrocoupling with volumetric monitoring", q(25, "°C"), q(10, "h"), "Ar")]
    if slug == "catalysis_complex5":
        outs = [outcome(rid, 1, "no_reaction", ["hydrogen"], evidence, result, basis="volumetric gas-evolution trace")]
        coverage = "not_applicable"; changes = None
    elif slug == "bh3_poison":
        outs = [outcome(rid, 1, "other", ["hydrogen"], evidence, result, basis="volumetric gas-evolution trace")]
        coverage = None; changes = DEHYDRO
    else:
        outs = [outcome(rid, 1, "product_amount", ["hydrogen"], evidence, result, basis="volumetric gas evolution relative to initial H3B·NMe2H")]
        if full:
            outs.append(outcome(rid, 2, "conversion", ["dimethylamine_borane"], evidence, q(100, "%"), basis="inferred by authors from approximately 1.0 equivalent H2 formation"))
        coverage = None; changes = DEHYDRO
    reaction(slug, "specific_experiment", evidence, ["control", "mechanistic_study"], note, "amine-borane dehydrocoupling",
             parts, st, outs, general=("scheme1_current_dehydrocoupling",), class_changes=changes, class_coverage=coverage,
             notes=["Figure 2 caption gives 2.5 mol% precatalyst, toluene, 25 °C."],
             unresolved=["Individual borazane yields are not reported for this run." if slug != "catalysis_complex5" else "No product was invented for the inactive control."])


# Scheme 2 syntheses, enriched only with the matching experimental-section text as necessary context (pages 13-14 not claimed as covered).
e4_13 = ev(4) + [et(13)]
rid = f"{DOC}:synthesis_complex4"
parts = [participant(rid, 1, "complex1", "reactant", e4_13, amounts={"mass": q(400, "mg"), "amount": q(0.57, "mmol")}),
         participant(rid, 2, "li_nme2", "reagent", e4_13, amounts={"mass": q(58.6, "mg"), "amount": q(1.15, "mmol")}),
         participant(rid, 3, "diethyl_ether", "solvent", e4_13, amounts={"volume": q(15, "mL")}),
         participant(rid, 4, "benzene", "solvent", e4_13), participant(rid, 5, "complex4", "product", e4_13, product_role="isolated product"),
         participant(rid, 6, "licl", "product", e4_13, product_role="salt coproduct")]
pids = [p["participant_id"] for p in parts]
st = [step(rid, 1, e4_13, pids[:3], "mix in Et2O and stir", None, "overnight"),
      step(rid, 2, e4_13, pids, "remove Et2O, dissolve residue in benzene, filter, wash with benzene 3 x 2 mL, concentrate", stage="workup"),
      step(rid, 3, e4_13, pids, "crystallize from pentane", stage="purification")]
outs = [outcome(rid, 1, "product_amount", ["complex4"], e4_13, q(326, "mg"), basis="isolated crystals"),
        outcome(rid, 2, "yield", ["complex4"], e4_13, q(82, "%"), basis="isolated crystals", yield_type="isolated")]
reaction("synthesis_complex4", "specific_experiment", e4_13, ["other"], "synthesis of complex 4", "salt metathesis / amide substitution",
         parts, st, outs, class_changes=[("N", "Zr", 0, 1, "formed", "Two Zr-NMe2 linkages are formed; represented as one class-level change."),
                                                ("Cl", "Zr", 1, 0, "broken", "Zr-Cl linkages are displaced with LiCl formation.")],
         notes=["Page 13 text was used as necessary procedure context; page 13 is not included in this session's visual coverage."],
         unresolved=["Reaction temperature for the overnight Et2O stir is unreported."])

e4_14 = ev(4,5) + [et(14)]
rid = f"{DOC}:synthesis_complex5"
parts = [participant(rid, 1, "complex1", "reactant", e4_14, amounts={"mass": q(200, "mg"), "amount": q(0.26, "mmol")}),
         participant(rid, 2, "nabh4", "reagent", e4_14, amounts={"mass": q(20, "mg"), "amount": q(0.53, "mmol")}),
         participant(rid, 3, "diethyl_ether", "solvent", e4_14, amounts={"volume": q(15, "mL")}),
         participant(rid, 4, "benzene", "solvent", e4_14), participant(rid, 5, "complex5", "product", e4_14, product_role="isolated product"),
         participant(rid, 6, "nacl", "product", e4_14, product_role="salt coproduct")]
pids = [p["participant_id"] for p in parts]
st = [step(rid, 1, e4_14, pids[:3], "stir in Et2O", q(40, "°C"), q(1, "week") if False else "one week"),
      step(rid, 2, e4_14, pids, "cannula-filter, remove volatiles, dissolve foam in benzene, warm and slowly cool overnight", None, "overnight", stage="workup"),
      step(rid, 3, e4_14, pids, "collect orange crystals", stage="purification")]
outs = [outcome(rid, 1, "product_amount", ["complex5"], e4_14, q(152, "mg"), basis="isolated crystals"),
        outcome(rid, 2, "yield", ["complex5"], e4_14, q(76, "%"), basis="isolated crystals", yield_type="isolated")]
reaction("synthesis_complex5", "specific_experiment", e4_14, ["other"], "synthesis of complex 5", "salt metathesis / borohydride substitution",
         parts, st, outs, class_changes=None,
         notes=["The source states that complex 5 could not be synthesized in high purity for melting point and CHN analysis; page 14 text is context, not claimed visual coverage."],
         unresolved=["Exact Zr-B/H connectivity changes are not safely classifiable from the compact drawing; product purity is limited."])


# Page 5 activation experiments.
e5 = ev(5)
rid = f"{DOC}:activation_uvvis_2eq"
parts = [participant(rid, 1, "complex1", "reactant", e5), participant(rid, 2, "me_li", "reagent", e5, amounts={"equivalents": q(2, "equiv")}, reference="complex1"),
         participant(rid, 3, "toluene", "solvent", e5), participant(rid, 4, "complex6", "product", e5, product_role="assigned intermediate"),
         participant(rid, 5, "complex2", "product", e5, product_role="methylation product"), participant(rid, 6, "licl", "product", e5, product_role="salt coproduct")]
st = [step(rid, 1, e5, [p["participant_id"] for p in parts], "add two equivalents MeLi and monitor by operando UV-vis", q(25, "°C"), q(60, "min"))]
outs = [outcome(rid, 1, "other", ["complex6", "complex2"], e5, "bands evolve at 352 and 424 nm with minimum at 335 nm; more than two UV-vis-active species", basis="UV-vis spectroscopy and singular-value analysis")]
r = reaction("activation_uvvis_2eq", "specific_experiment", e5, ["mechanistic_study"], "UV-vis monitoring of precatalyst activation",
             "stepwise methylation", parts, st, outs, general=("scheme_activation_sequence",), class_changes=ZR_C_TWO,
             unresolved=["Assignment of the second UV-vis pattern to complex 6 assumes complex 2 at full conversion; multiple absorbing/aggregated species prevent quantitative kinetics."])
r["solvent_conditions"]["composition"] = f("toluene", [M["toluene"]["molecule_id"]], e5)
r["solvent_conditions"]["concentration"] = f("[1] = 3.1 x 10^-4 M; [MeLi] = 6.3 x 10^-4 M", {"complex1": q(3.1e-4, "M"), "MeLi": q(6.3e-4, "M")}, e5)

rid = f"{DOC}:activation_nmr_2eq"
parts = [participant(rid, 1, "complex1", "reactant", e5), participant(rid, 2, "me_li", "reagent", e5),
         participant(rid, 3, "benzene", "solvent", e5), participant(rid, 4, "complex2", "product", e5, product_role="dominant product")]
st = [step(rid, 1, e5, [p["participant_id"] for p in parts], "mix and analyze by 1H NMR", "room temperature")]
outs = [outcome(rid, 1, "other", ["complex2"], e5, "rapid/exclusive formation; Cp resonances at 5.90 and 5.92 ppm", basis="1H NMR")]
r = reaction("activation_nmr_2eq", "specific_experiment", e5, ["mechanistic_study"], "concentrated NMR activation with MeLi",
             "methylation to complex 2", parts, st, outs, general=("scheme_activation_sequence",), class_changes=ZR_C_TWO,
             unresolved=["Caption concentration ratio [MeLi]/[1] = 2.27 rather than exactly 2; no reaction time is stated."])
r["solvent_conditions"]["composition"] = f("benzene", [M["benzene"]["molecule_id"]], e5)
r["solvent_conditions"]["concentration"] = f("[1] = 6.6 x 10^-2 M; [MeLi] = 1.5 x 10^-1 M", {"complex1": q(6.6e-2, "M"), "MeLi": q(1.5e-1, "M")}, e5)

rid = f"{DOC}:activation_nmr_1eq"
parts = [participant(rid, 1, "complex1", "reactant", e5), participant(rid, 2, "me_li", "reagent", e5, amounts={"equivalents": q(1, "equiv")}, reference="complex1"),
         participant(rid, 3, "benzene", "solvent", e5), participant(rid, 4, "complex6", "product", e5, product_role="main species"),
         participant(rid, 5, "complex2", "product", e5, product_role="minor species")]
st = [step(rid, 1, e5, [p["participant_id"] for p in parts], "add one equivalent MeLi and analyze by 1H NMR", "room temperature")]
outs = [outcome(rid, 1, "other", ["complex6"], e5, "main species; Cp resonances 6.05, 6.05, 6.00, and 5.94 ppm", basis="1H NMR")]
r = reaction("activation_nmr_1eq", "specific_experiment", e5, ["mechanistic_study"], "one-equivalent NMR activation",
             "formation of monomethyl complex 6", parts, st, outs, general=("scheme_activation_sequence",), class_changes=ZR_C_ONE,
             unresolved=["Exact reaction time and product ratio are not reported."])
r["solvent_conditions"]["composition"] = f("benzene", [M["benzene"]["molecule_id"]], e5)
r["solvent_conditions"]["concentration"] = f("[1] = [MeLi] = 3.3 x 10^-2 M", {"complex1": q(3.3e-2, "M"), "MeLi": q(3.3e-2, "M")}, e5)

rid = f"{DOC}:meli_uvvis_dilution"
parts = [participant(rid, 1, "me_li", "reactant", e5), participant(rid, 2, "toluene", "solvent", e5)]
st = [step(rid, 1, e5, [p["participant_id"] for p in parts], "stepwise dilute MeLi in toluene and monitor by UV-vis")]
outs = [outcome(rid, 1, "other", ["me_li"], e5, "nonlinear extinction versus concentration; at least two UV-vis-active species", basis="UV-vis spectroscopy and singular-value decomposition")]
reaction("meli_uvvis_dilution", "specific_experiment", e5, ["mechanistic_study", "control"], "MeLi aggregation/dilution control",
         "physical speciation study", parts, st, outs, class_coverage="not_applicable",
         unresolved=["Individual dilution concentrations and aggregation identities are in unavailable SI figures S80/S81."])


# Page 6 isolation/crystallographic composition of 6; detailed page-14 procedure is context.
e6_14 = ev(5,6) + [et(14)]
rid = f"{DOC}:isolation_complex6"
source_meli_volume = fact(raw="0.27 μL, 1.6 M in Et2O", normalized=None, evidence=[et(14)], status="conflict", verification="unreviewed", derivation=None,
                          alternatives=[f("0.27 μL as printed", q(0.27, "uL"), [et(14)]),
                                        f("0.43 mmol at 1.6 M implies about 0.27 mL", q(0.27, "mL"), [et(14)], method="inference", derivation="Amount/concentration arithmetic; not chosen over the printed unit.")])
parts = [participant(rid, 1, "complex1", "reactant", e6_14, amounts={"mass": q(300, "mg"), "amount": q(0.43, "mmol")}),
         participant(rid, 2, "me_li", "reagent", e6_14, amounts={"amount": q(0.43, "mmol"), "volume": {"_fact": source_meli_volume}, "equivalents": q(1, "equiv")}, reference="complex1"),
         participant(rid, 3, "benzene", "solvent", e6_14, amounts={"volume": q(7, "mL")}), participant(rid, 4, "diethyl_ether", "solvent", e6_14),
         participant(rid, 5, "complex6", "product", e6_14, product_role="isolated mixed crystals"), participant(rid, 6, "complex2", "product", e6_14, product_role="co-crystallized minor component")]
pids = [p["participant_id"] for p in parts]
st = [step(rid, 1, e6_14, pids[:3], "add MeLi slowly in benzene and stir", q(40, "°C"), "overnight"),
      step(rid, 2, e6_14, pids, "filter, remove benzene, dissolve in Et2O and concentrate", stage="workup"),
      step(rid, 3, e6_14, pids, "crystallize mixed 6/2 material from Et2O", stage="purification")]
outs = [outcome(rid, 1, "product_amount", ["complex6"], e6_14, q(122, "mg"), basis="isolated orange crystals"),
        outcome(rid, 2, "yield", ["complex6"], e6_14, q(40, "%"), basis="isolated orange crystals", yield_type="isolated"),
        outcome(rid, 3, "selectivity", ["complex6", "complex2"], e6_14, "89:11", basis="combined crystallography and 1H NMR", scope={"metric": "crystal composition", "order": [M["complex6"]["molecule_id"], M["complex2"]["molecule_id"]]}),
        outcome(rid, 4, "selectivity", ["complex6"], ev(6), "48:34:18", basis="crystallographic occupancy with calculated thermodynamic comparison", scope={"metric": "orientation-isomer distribution", "order": ["inCl-inMe", "inCl-outMe", "outCl-inMe"]})]
reaction("isolation_complex6", "specific_experiment", e6_14, ["mechanistic_study", "other"], "synthesis and crystallization of complex 6",
         "single methylation of complex 1", parts, st, outs, general=("scheme_activation_sequence",), class_changes=ZR_C_ONE,
         notes=["The source prints MeLi volume as 0.27 μL while also stating 0.43 mmol at 1.6 M; both the printed value and arithmetic alternative are preserved with no chosen correction.",
                "Page 14 procedure text is necessary context and is not claimed as visually covered."],
         unresolved=["Printed MeLi volume/amount/concentration conflict remains unresolved.", "The isolated crystals contain 11% complex 2 and three orientation isomers of 6."])


# Page 6 Figure 5 alternative activation/catalyst protocols. Parent 1/MeLi curve is already catalyst_loading_2_50.
alt_defs = [
    ("catalysis_complex6_meli", "complex6", True, None, "similar volumetric profile to 1/MeLi", "complex 6/MeLi activation"),
    ("catalysis_complex2", "complex2", False, None, "full conversion but considerably slower than 1/MeLi", "isolated complex 2 as precatalyst"),
    ("catalysis_in_situ_complex2", "complex1", True, q(1, "h"), "much less H2 than direct 1/MeLi", "one-hour preactivation to in situ complex 2"),
    ("catalysis_li_amidoborane", "complex1", False, None, "similar conversion to 1/MeLi with a significant induction period", "Li[NMe2BH3] cocatalyst"),
]
for slug, catalyst, use_meli, pretime, result, title in alt_defs:
    evidence = ev(6) + [et(15)]
    rid = f"{DOC}:{slug}"
    if slug == "catalysis_li_amidoborane":
        extra = [("li_amidoborane", "reagent", {"mol_percent": q(5.0, "mol%")}, "dimethylamine_borane")]
    else:
        extra = []
    parts = catalytic_parts(rid, evidence, catalyst_slug=catalyst, catalyst_loading=2.5, me_li=use_meli, extra=extra)
    pids = [p["participant_id"] for p in parts]
    st = []
    if pretime:
        st.append(step(rid, 1, evidence, pids, "preactivate complex 1 with MeLi before substrate addition", q(25, "°C"), pretime))
        st.append(step(rid, 2, evidence, pids, "add substrate and monitor H2 evolution", q(25, "°C"), None, "Ar"))
    else:
        st.append(step(rid, 1, evidence, pids, "combine catalyst system and substrate; monitor H2 evolution", q(25, "°C"), None, "Ar"))
    outs = [outcome(rid, 1, "other", ["hydrogen"], evidence, result, basis="volumetric gas-evolution trace")]
    reaction(slug, "specific_experiment", evidence, ["mechanistic_study", "control"], title, "amine-borane dehydrocoupling",
             parts, st, outs, general=("scheme1_current_dehydrocoupling",), class_changes=DEHYDRO,
             unresolved=["Exact endpoint/time is read qualitatively from the source discussion; individual borazane yields are unreported."])


# Add page-6 evidence to the already-counted direct 1/MeLi baseline record rather than duplicate it.
baseline = next(r for r in R if r["reaction_id"].endswith("catalyst_loading_2_50"))
for eid in ev(6):
    if eid not in baseline["evidence_ids"]:
        baseline["evidence_ids"].append(eid)
baseline["notes"].append(f("Figure 5 reuses the direct 1/MeLi baseline for comparison; no duplicate experiment record was created.",
                            "Figure 5 reuses the direct 1/MeLi baseline for comparison; no duplicate experiment record was created.", ev(6)))


bundle = json.loads((RAW / "standard_init.json").read_text())
bundle["molecules"] = list(M.values())
bundle["reactions"] = R
bundle["source_context"] = [
    {"context_id": f"{DOC}:context_metadata", "document_id": DOC, "evidence_ids": ev(1,2),
     "kind": "bibliographic_metadata", "title": "Mechanistic insights into dehydrocoupling of amine boranes using dinuclear zirconocene complexes",
     "doi": "10.1039/d1cy00531f", "note": "The paper advertises ESI, but no SI file was supplied or associated."},
    {"context_id": f"{DOC}:context_spectroscopy", "document_id": DOC, "evidence_ids": ev(3,4,5),
     "kind": "spectroscopy_and_kinetics", "note": "11B NMR peak assignments, UV-vis traces, NMR resonances, and singular-value discussion support reaction records but are not separate syntheses."},
    {"context_id": f"{DOC}:context_xrd", "document_id": DOC, "evidence_ids": ev(6),
     "kind": "structure_assignment", "note": "Figure 4 X-ray structure and crystal occupancies support complex 6 identity/composition; the structure plot is not a synthesis experiment."},
    {"context_id": f"{DOC}:context_negative_spectral_inventory", "document_id": DOC, "evidence_ids": ev(1,2,3,4,5,6),
     "kind": "negative_inventory", "note": "No standalone mass-spectrum peak table, molecular-weight table, or fragment-relative-abundance block on pages 1-6 was converted into a synthesis record."},
    {"context_id": f"{DOC}:context_cover_art", "document_id": DOC, "evidence_ids": ev(1),
     "kind": "cover_art", "note": "The graphical cover summary is illustrative and duplicates the paper's reaction/speciation story; it is not counted as an experiment."},
]
write(OUT / "submitted.json", bundle)

# Raw visual inventory and deduplication decisions.
page_records = {page: [] for page in range(1, 7)}
for r in R:
    for page in range(1, 7):
        if et(page) in r["evidence_ids"] or ei(page) in r["evidence_ids"]:
            page_records[page].append(r["reaction_id"])

inventory = {
    "document_id": DOC,
    "pages_visually_inspected": [1, 2, 3, 4, 5, 6],
    "source_images": {str(p): f"assets/{DOC}/page{p}.jpg" for p in range(1, 7)},
    "page_inventory": [
        {"page": 1, "items": ["cover artwork and graphical reaction/speciation summary", "no independent experimental table or arrow"], "record_ids": page_records[1]},
        {"page": 2, "items": ["title/abstract/DOI", "Scheme 1a cited prior work", "Scheme 1b current generic transformation", "ESI availability notice only"], "record_ids": page_records[2]},
        {"page": 3, "items": ["Figure 1 six catalyst-loading curves", "Figure 1 11B NMR product-distribution run", "Scheme 1 repeated generic transformation"], "record_ids": page_records[3]},
        {"page": 4, "items": ["three deuterated-substrate mechanistic experiments", "Figure 2 baseline plus four distinct poison/complex comparisons", "Scheme 2 two synthesis arrows"], "record_ids": page_records[4]},
        {"page": 5, "items": ["Scheme/figure 3 stepwise methylation", "UV-vis 2-equiv activation", "two NMR activation experiments", "MeLi dilution/speciation control", "complex 5 inactivity continuation"], "record_ids": page_records[5]},
        {"page": 6, "items": ["complex 6 isolation/XRD composition", "Figure 5 direct baseline duplicate plus four alternative catalyst/activation protocols"], "record_ids": page_records[6]},
    ],
    "deduplication": [
        "Figure 5 direct 1/MeLi curve was merged into catalyst_loading_2_50.",
        "Scheme 1b is a generic scheme and was not counted as a concrete catalytic run.",
        "The page-1 cover graphic was retained as source context only.",
        "NMR/UV-vis/MS/XRD observations were not converted into standalone synthesis records unless they described a real activation/control experiment.",
    ],
    "limitations": [
        "Only pages 1-6 were visually inspected and claimed as covered.",
        "Pages 13-15 extracted text was used only as necessary matching procedure context; those pages were not visually covered.",
        "ESI was advertised but not supplied and was not associated.",
    ],
}
write(RAW / "page_extraction.json", inventory)

coverage_regions = {
    1: [{"region": "full-page cover artwork and caption", "finding": "Illustrative dehydrocoupling/speciation summary; no independent experiment counted.", "evidence_ids": [ei(1), et(1)]}],
    2: [{"region": "title, abstract, DOI and ESI notice", "finding": "Bibliographic metadata and unavailable-ESI notice.", "evidence_ids": [ei(2), et(2)]},
        {"region": "Scheme 1a", "finding": "Explicitly cited previous dehydropolymerisation; stored as literature_example.", "evidence_ids": [ei(2), et(2)]},
        {"region": "Scheme 1b", "finding": "Current generic dehydrocoupling transformation; not a uniquely quantified run.", "evidence_ids": [ei(2), et(2)]}],
    3: [{"region": "Figure 1 top", "finding": "Six catalyst-loading gas-evolution experiments inventoried separately.", "evidence_ids": [ei(3), et(3)]},
        {"region": "Figure 1 bottom and adjacent prose", "finding": "Distinct 11B NMR product-distribution experiment with main/minor/residual species.", "evidence_ids": [ei(3), et(3)]},
        {"region": "Scheme 1 and Results prose", "finding": "Repeated generic scheme and current protocol context; no duplicate record.", "evidence_ids": [ei(3), et(3)]}],
    4: [{"region": "left-column isotope discussion", "finding": "Three distinct deuterated-substrate mechanistic experiments; exact SI rates unavailable.", "evidence_ids": [ei(4), et(4)]},
        {"region": "Figure 2 and poison/catalyst prose", "finding": "Baseline plus NMe2H, BH3·THF, complex 4 and complex 5 comparisons; baseline deduplicated.", "evidence_ids": [ei(4), et(4)]},
        {"region": "Scheme 2", "finding": "Two concrete synthesis arrows for complexes 4 and 5.", "evidence_ids": [ei(4), et(4)]}],
    5: [{"region": "Figure 3 top and activation prose", "finding": "Two-equivalent UV-vis activation plus stepwise 1→6→2 scheme.", "evidence_ids": [ei(5), et(5)]},
        {"region": "Figure 3 bottom", "finding": "One-equivalent and concentrated two-equivalent NMR activation experiments.", "evidence_ids": [ei(5), et(5)]},
        {"region": "left-column MeLi dilution discussion", "finding": "Physical MeLi aggregation/speciation control; no chemical product invented.", "evidence_ids": [ei(5), et(5)]}],
    6: [{"region": "left column and Figure 4", "finding": "Complex 6 isolation and mixed-crystal/XRD composition; spectroscopy/structure not duplicated as synthesis.", "evidence_ids": [ei(6), et(6)]},
        {"region": "Figure 5 top", "finding": "Direct baseline reused; complex 6/MeLi, isolated 2, and one-hour in situ 2 protocols recorded separately.", "evidence_ids": [ei(6), et(6)]},
        {"region": "Figure 5 bottom", "finding": "MeLi baseline reused and Li[NMe2BH3] cocatalyst experiment recorded separately.", "evidence_ids": [ei(6), et(6)]}],
}
coverage = {"document_id": DOC, "semantic_accuracy_verified": False, "pages": []}
for page in range(1, 7):
    unresolved = []
    if page == 1:
        unresolved.append("Cover artwork is not an experimental protocol.")
    if page == 2:
        unresolved.append("Advertised ESI was not supplied; proposed intermediates are not treated as executed experiments.")
    if page == 3:
        unresolved.append("Exact numerical endpoints for five of six catalyst-loading curves are not stated in prose.")
    if page == 4:
        unresolved.append("Isotope-study numerical rates and full protocols are in unavailable SI figures S5/S6.")
    if page == 5:
        unresolved.append("MeLi aggregation species and complete dilution series are unresolved; operando UV-vis is not kinetically interpretable.")
    if page == 6:
        unresolved.append("Several Figure 5 curves have qualitative rather than textually quantified outcomes; no values were digitized from plots.")
    coverage["pages"].append({"page": page, "state": "attempted", "regions": coverage_regions[page],
                              "record_ids": page_records[page], "unresolved": unresolved})
write(OUT / "coverage" / "pages_0001_0006.json", coverage)

# Save graph diagnostics separately as raw machine output.
stereo = {}
for slug, m in M.items():
    smi = m["smiles"]["normalized_value"]
    if smi:
        stereo[m["molecule_id"]] = stereo_inspect(smi)
write(RAW / "stereo_checks.json", stereo)

prepared_dataset = json.loads((OUT / "dataset.json").read_text())
mentions = source_mentions(prepared_dataset,
                           ["H3B·NMe2H", "Li[NMe2BH3]", "complex 4", "complex 5", "complex 6", "two equivalents of MeLi"],
                           context=320)
write(RAW / "source_mentions.json", mentions)

summary = {"document_id": DOC, "molecules": len(M), "reactions": len(R),
           "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in R),
           "generic_schemes": sum(r["record_type"] == "generic_scheme" for r in R),
           "literature_examples": sum(r["record_type"] == "literature_example" for r in R),
           "all_needs_review": all(r["review"]["status"] == "needs_review" for r in R),
           "pages_covered": [1, 2, 3, 4, 5, 6]}
write(RAW / "verification_summary.json", summary)
print(json.dumps(summary, ensure_ascii=False, indent=2))
