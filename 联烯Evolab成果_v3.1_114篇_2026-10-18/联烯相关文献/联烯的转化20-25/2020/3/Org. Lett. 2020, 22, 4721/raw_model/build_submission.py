#!/usr/bin/env python3
"""Build the source-derived v3.1 submission for DOI 10.1021/acs.orglett.0c01523.

This is a reproducible machine transcription of the five supplied main-text pages.
It deliberately leaves unavailable SI procedures and unsafe dense product graphs unresolved.
"""
from pathlib import Path
import copy
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_d9d234c1ff998484"
E = {p: f"{DOC}:page{p}:image" for p in range(1, 6)}
T = {p: f"{DOC}:page{p}:text" for p in range(1, 6)}


def f(raw=None, normalized=None, pages=(), method="original_extraction", derivation=None,
      status=None):
    evidence = []
    for page in pages:
        evidence.extend([E[page], T[page]])
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation,
                   status=status)


def inferred(raw, normalized, pages, why):
    return f(raw, normalized, pages, method="inference", derivation=why)


def derived(raw, normalized, pages, why):
    return f(raw, normalized, pages, method="structure_derivation", derivation=why)


bundle = st.read(RUN / "submitted.json")
bundle["molecules"] = []
bundle["reactions"] = []
bundle["source_context"] = []
stereo_checks = {}


def molecule(label, pages, name=None, smiles=None, chirality="unresolved",
             geometry="unresolved", representation=None, stereo_raw=None,
             identity_note=None):
    mid = f"{DOC}:mol:{label}"
    ev = [x for p in pages for x in (E[p], T[p])]
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, pages)
    m["names"]["original"] = f(name or label, name or label, pages)
    if name:
        m["names"]["english"] = f(name, name, pages)
    m["cas_numbers"] = f()
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = derived(
            f"source drawing for {label}", smiles, pages,
            "Candidate connectivity transcribed from the supplied source drawing; parsing does not verify source identity."
        )
        diag = st.stereo_inspect(smiles)
        stereo_checks[mid] = diag
        m["stereochemistry"]["isomeric_smiles"] = derived(
            smiles, diag["canonical_isomeric_smiles"], pages,
            "RDKit canonical isomeric SMILES for the candidate graph; allene axial configuration is not certified."
        )
        m["stereochemistry"]["identity_key"] = derived(
            None, diag["stereo_graph_key"], pages,
            "RDKit stereo graph key for comparison only; not a sample identity."
        )
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = f()
    if representation is None:
        representation = "partial" if chirality in ("racemic", "specified", "partially_specified", "unspecified") else "unresolved"
    status_why = "Status follows the source drawing/text and does not infer an unreported stereoisomer."
    m["stereochemistry"]["chirality_status"] = (
        f(stereo_raw or chirality, chirality, pages)
        if chirality in ("racemic", "specified", "partially_specified", "unspecified", "mixture")
        else inferred(stereo_raw or chirality, chirality, pages, status_why)
    )
    m["stereochemistry"]["geometry_status"] = inferred(
        geometry, geometry, pages,
        "Geometry status is assigned from the depicted functional groups; no unreported E/Z descriptor is introduced."
    )
    m["stereochemistry"]["representation_status"] = inferred(
        representation, representation, pages,
        "Representation status compares the candidate/unresolved graph with the stereochemical information visible in the source."
    )
    if chirality == "racemic":
        m["stereochemistry"]["axial"] = f(
            "racemic γ-substituted allene", [{"axis": "C=C=C", "assignment": "racemic; individual R_a/S_a not assigned"}], pages
        )
        m["stereochemistry"]["mixture"] = f(
            "racemic", {"type": "racemate", "ratio": "not reported"}, pages
        )
        m["stereochemistry"]["assignment"] = f(
            "racemic substrate; no individual axial descriptor", "axial configuration unassigned", pages
        )
    elif chirality == "unspecified" and "allene" in (name or "").casefold():
        m["stereochemistry"]["axial"] = f(
            "allene axis shown without R_a/S_a", [{"axis": "C=C=C", "assignment": "unspecified"}], pages
        )
        m["stereochemistry"]["assignment"] = f(
            "configuration not specified", "unassigned", pages
        )
    elif chirality in ("specified", "partially_specified"):
        m["stereochemistry"]["atom_centers"] = f(
            stereo_raw or "wedge/dash stereochemistry in source drawing",
            [{"assignment": "source-depicted stereochemistry; exact R/S descriptors not safely reconstructed"}], pages
        )
        m["stereochemistry"]["assignment"] = f(
            stereo_raw or "source-depicted stereochemistry",
            "specified in source but not encoded in a complete SMILES", pages
        )
    elif chirality == "achiral":
        m["stereochemistry"]["assignment"] = inferred(
            "no stereogenic element in candidate graph", "achiral", pages,
            "The candidate graph/drawing has no stereogenic center or substituted allene axis."
        )
    m["identity_note"] = f(identity_note, identity_note, pages) if identity_note else f()
    bundle["molecules"].append(m)
    return mid


# Azaheptafulvenes and catalysts are retained by labels because the dense drawings and
# sulfonyl/catalyst stereochemistry should not be promoted to unsafe complete SMILES.
M = {}
for label, name, pages in [
    ("1a", "N-nosyl azaheptafulvene", (1, 2)),
    ("1b", "N-tosyl azaheptafulvene", (2,)),
    ("1c", "N-(4-methoxyphenylsulfonyl) azaheptafulvene", (2,)),
    ("1d", "N-mesyl azaheptafulvene", (2,)),
]:
    M[label] = molecule(label, pages, name=name, chirality="achiral", geometry="not_applicable",
                        representation="unresolved", identity_note="Complete connectivity retained by source label/drawing; no CAS supplied.")

allene_smiles = {
    "2a": "CCOC(=O)C=C=C",
    "2b": "CC(C)(C)OC(=O)C=C=C",
    "2c": "O=C(OC1=CC=CC=C1)C=C=C",
    "2d": "O=C(Sc1ccccc1)C=C=C",
    "2e": "O=C(c1ccccc1)C=C=C",
    "2f": "O=C(N(c1ccccc1)c1ccccc1)C=C=C",
    "2g": "O=C(N(c1ccccc1)c1ccccc1)C=C=Cc1ccccc1",
    "2h": "O=C(N(c1ccccc1)c1ccccc1)C=C=Cc1ccc(Cl)cc1",
    "2i": "O=C(N(c1ccccc1)c1ccccc1)C=C=Cc1ccc(F)cc1",
    "2j": "O=C(N(c1ccccc1)c1ccccc1)C=C=Cc1ccc(OC)cc1",
    "2k": "O=C(N(c1ccccc1)c1ccccc1)C=C=Cc1ccc2c(c1)OCO2",
    "2l": "O=C(N(c1ccccc1)c1ccccc1)C=C=Cc1ccccc1Cl",
    "2m": "O=C(N(c1ccccc1)c1ccccc1)C=C=Cc1cccc(Br)c1",
    "2n": "O=C(N(c1ccccc1)c1ccccc1)C=C=Cc1cc(OC)c(OC)c(OC)c1",
    "2o": "O=C(N(c1ccccc1)c1ccccc1)C=C=Cc1ccc2ccccc2c1",
    "2p": "O=C(N(c1ccccc1)c1ccccc1)C=C=CC",
    "2q": "O=C(N(c1ccccc1)c1ccccc1)C=C=CCC",
    "2r": "O=C(N(c1ccccc1)c1ccccc1)C=C=CC(C)C",
}
allene_names = {
    "2a": "ethyl 2,3-butadienoate", "2b": "tert-butyl 2,3-butadienoate",
    "2c": "phenyl 2,3-butadienoate", "2d": "phenyl 2,3-butadienethioate",
    "2e": "allene phenyl ketone 2e", "2f": "N,N-diphenyl 2,3-butadienamide",
}
for label, smi in allene_smiles.items():
    substituted = label not in {"2a", "2b", "2c", "2d", "2e", "2f"}
    M[label] = molecule(
        label, (2, 3) if substituted else (1, 2),
        name=allene_names.get(label, f"γ-substituted allenic amide {label}"), smiles=smi,
        chirality="racemic" if substituted else "achiral", geometry="not_applicable",
        representation="partial" if substituted else "not_applicable",
        stereo_raw="racemic γ-substituted allenic amide" if substituted else "terminal allene candidate graph",
        identity_note=("Candidate regioisomer transcribed from Scheme 2; dense aryl substitution and allene-axis encoding require expert review."
                       if substituted else "Candidate connectivity from Table 1/Scheme 2; structure remains machine-derived.")
    )

for label in ["3a", "3b", "3c", "3d", "3e", "3f", "3g"]:
    M[label] = molecule(
        label, (2,), name=f"chiral amino-acid-derived amidophosphine {label}",
        chirality="specified", geometry="not_applicable", representation="not_encoded",
        stereo_raw="wedge/dash stereochemistry shown in Table 1 catalyst drawing",
        identity_note="Dense catalyst graph and exact absolute descriptors are not reconstructed without the SI."
    )

M["toluene"] = molecule("toluene", (2,), name="toluene", smiles="Cc1ccccc1",
                         chirality="achiral", geometry="not_applicable", representation="not_applicable")
M["nBu3P"] = molecule("nBu3P", (5,), name="tri-n-butylphosphine", smiles="CCCCP(CCCC)CCCC",
                       chirality="achiral", geometry="not_applicable", representation="not_applicable")
M["ethyl_hexa"] = molecule(
    "ethyl_hexa", (5,), name="ethyl hexa-2,3-dienoate", smiles="CCOC(=O)C=C=CCC",
    chirality="unspecified", geometry="not_applicable", representation="partial",
    stereo_raw="γ-substituted allenoate; axial configuration not reported",
    identity_note="Footnote 12 substrate; no racemic composition or R_a/S_a assignment is reported."
)

# Products: the main-text drawings contain stereochemical wedges, but complete graphs
# and exact R/S descriptors are intentionally not reconstructed from the dense artwork.
for i, label in enumerate([f"4{x}" for x in "abcdefghijklmnopqrstu"]):
    pages = (2, 3) if label in {"4j", "4k", "4s", "4t", "4u"} else (2,)
    assignment = "X-ray absolute/relative structure" if label in {"4g", "4k", "4s"} else "source-depicted/assigned stereochemical outcome"
    M[label] = molecule(label, pages, name=f"[8+2] cycloadduct {label}", chirality="specified",
                        geometry="not_applicable", representation="not_encoded", stereo_raw=assignment,
                        identity_note="Dense azaazulenoid product retained by source label and drawing; complete stereochemical SMILES unavailable from supplied main text.")
for label in ["5a", "5b", "5c"]:
    assignment = "X-ray absolute configuration; 90% ee" if label == "5a" else "source-depicted [8+4] cycloadduct"
    M[label] = molecule(label, (2,), name=f"[8+4] cycloadduct {label}", chirality="specified",
                        geometry="not_applicable", representation="not_encoded", stereo_raw=assignment,
                        identity_note="Product retained by source label/drawing; no complete stereochemical SMILES asserted.")
M["5_ethyl_hexa"] = molecule(
    "5_ethyl_hexa", (5,), name="minor [8+4] cycloaddition product of type 5",
    chirality="unresolved", geometry="unresolved", representation="unresolved",
    identity_note="Footnote 12 reports only a product of type 5 detected by NMR; exact structure and stereochemistry are not printed."
)


def base_reaction(suffix, pages, record_type="specific_experiment", categories=None,
                  location=None):
    ev = [x for p in pages for x in (E[p], T[p])]
    r = st.new_reaction(f"{DOC}:{suffix}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories or [], categories or [], pages)
    r["location"] = f(location or suffix, {"pages": list(pages), "region": location or suffix}, pages)
    r["reaction_name"] = f("enantioselective higher order [8+2] cycloaddition", "enantioselective [8+2] cycloaddition", pages)
    r["transformation_type"] = f("higher order cycloaddition", "[8+2] annulation", pages)
    r["review"] = {
        "status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
        "unresolved": ["No explicitly associated SI was supplied; detailed procedure, workup, purification, and characterization remain unavailable."]
    }
    r["completeness"] = "partial"
    return r


def participant(r, pid, mid_key, role, pages, amount=None, mol_percent=None):
    p = st.new_participant(pid, M[mid_key], role, [x for page in pages for x in (E[page], T[page])])
    if amount is not None:
        p["charges"]["amount"] = f(f"{amount} mmol", {"number": amount, "unit": "mmol"}, pages)
    if mol_percent is not None:
        p["charges"]["mol_percent"] = f(f"{mol_percent} mol %", {"number": mol_percent, "unit": "mol%"}, pages)
        p["equivalence_reference"] = inferred(
            None, r["participants"][0]["molecule_id"] if r["participants"] else M[mid_key], pages,
            "Catalyst-loading denominator is inferred as the displayed limiting azaheptafulvene substrate; the caption does not name the denominator."
        )
    r["participants"].append(p)
    return p


def set_conditions(r, pages, scale, temperature, time=None, catalyst="3g", cat_loading=10,
                   substrate1=None, substrate2=None, products=(), include_solvent=True):
    if substrate1:
        participant(r, "aza", substrate1, "reactant", pages, amount=scale)
    if substrate2:
        participant(r, "allene", substrate2, "reactant", pages, amount=scale)
    if catalyst:
        participant(r, "cat", catalyst, "catalyst", pages, mol_percent=cat_loading)
    if include_solvent:
        participant(r, "solvent", "toluene", "solvent", pages)
    for n, prod in enumerate(products, 1):
        participant(r, f"product{n}", prod, "product", pages)
    r["reaction_scale"] = f(
        f"{scale} mmol of 1 and 2", {"reference_molecule_id": M[substrate1], "amount": {"number": scale, "unit": "mmol"}}, pages
    ) if substrate1 and scale is not None else f()
    if include_solvent:
        r["solvent_conditions"]["composition"] = f("toluene", [M["toluene"]], pages)
        r["solvent_conditions"]["concentration"] = f("0.05 M", {"number": 0.05, "unit": "M"}, pages)
        r["solvent_conditions"]["concentration_basis"] = inferred(
            None, M[substrate1] if substrate1 else None, pages,
            "The caption gives 0.05 M but does not explicitly name its concentration basis; limiting azaheptafulvene is used as a provisional reference."
        ) if substrate1 else f()
    step = st.new_step("reaction", [x for page in pages for x in (E[page], T[page])])
    step["participant_ids"] = [p["participant_id"] for p in r["participants"] if p["role"] != "product"]
    step["operation"] = f("reaction under the displayed conditions", "combine and react", pages)
    if temperature is not None:
        raw_t = f"{temperature} °C" if isinstance(temperature, (int, float)) else temperature
        norm_t = {"number": temperature, "unit": "°C"} if isinstance(temperature, (int, float)) else temperature
        step["temperature"] = f(raw_t, norm_t, pages)
    if time is not None:
        if isinstance(time, tuple):
            raw_time = f"{time[0]}–{time[1]} h"
            norm_time = {"minimum": time[0], "maximum": time[1], "unit": "h", "raw_text": raw_time}
            step["time"] = f(raw_time, norm_time, pages, method="inherited_general_condition")
        else:
            step["time"] = f(f"{time} h", {"number": time, "unit": "h"}, pages)
    r["steps"] = [step]
    return r


def add_outcome(r, oid, kind, targets, pages, raw, normalized, product_role=None,
                yield_type=None, basis=None, scope=None):
    o = st.new_outcome(oid, kind, [M[t] for t in targets], [x for p in pages for x in (E[p], T[p])])
    o["value"] = f(raw, normalized, pages)
    if product_role:
        o["product_role"] = f(product_role, product_role, pages)
    if yield_type:
        o["yield_type"] = f(yield_type, yield_type, pages)
    if basis:
        o["measurement_basis"] = f(basis, basis, pages)
    if scope:
        o["scope"] = f(scope, scope, pages)
    r["outcomes"].append(o)
    return o


def classify_cycloaddition(r, pages, note="Two new C–C sigma connections visible in the cycloadduct drawing."):
    changes = [
        {"change_id": "b1", "atoms": [{"element": "C", "source": "azaheptafulvene terminus"}, {"element": "C", "source": "allene component"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b2", "atoms": [{"element": "C", "source": "azaheptafulvene terminus"}, {"element": "C", "source": "allene component"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
    ]
    r["classification"]["bond_changes"] = [derived(None, c, pages, note) for c in changes]
    r["classification"]["labels"] = derived(
        None, [{"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1", "b2"]}], pages,
        "Partial source-drawing classification: both ring-closing connections are C–C bonds."
    )
    r["classification"]["coverage"] = inferred(
        None, "partial", pages,
        "Exact atom mapping, pi-bond migration, and proton transfer are unresolved; only visually secure new C–C connections are asserted."
    )
    r["classification"]["limitations"] = [
        "No complete atom mapping is claimed.",
        "Implicit-H relocation and all bond-order changes require expert review of complete structures.",
        "Labels are derived from source drawings, not author-reported classification."
    ]


def classify_not_applicable(r, pages, why):
    r["classification"]["bond_changes"] = []
    r["classification"]["labels"] = inferred(None, [], pages, why)
    r["classification"]["coverage"] = inferred(None, "not_applicable", pages, why)
    r["classification"]["limitations"] = [why]


# Scheme 1: three cited prior-work arrows and one generic current-work arrow.
literature_findings = {
    "scheme1_eq1_literature": "Previous work: phosphine-catalyzed enantioselective [3+2] cycloaddition, γ-substituted allene as C3 component (Ref. 5).",
    "scheme1_eq2_literature": "Previous work: γ,δ-reactivity of an enolizable γ-substituted allenoate (Ref. 6).",
    "scheme1_eq3_literature": "Previous work: γ,δ-reactivity of a δ-acetoxy allenoate with a leaving group (Ref. 7).",
}
for suffix, text in literature_findings.items():
    r = base_reaction(suffix, (1,), "literature_example", [], "Scheme 1 previous-work arrow")
    r["reaction_name"] = f(text, text, (1,))
    r["transformation_type"] = f("cited prior work", "literature example", (1,))
    r["notes"] = [f(text, text, (1,))]
    r["review"]["unresolved"] = ["Generic cited drawing only; no current-paper experiment or conditions are asserted."]
    bundle["reactions"].append(r)

r = base_reaction("scheme1_eq4_generic", (1,), "generic_scheme", [], "Scheme 1 eq 4, this work")
r["notes"] = [f("This work: enantioselective [8+2] cycloaddition; γ-substituted allene acts as C2 component via β,γ-reactivity.",
                  "generic current-study transformation", (1,))]
r["review"]["unresolved"] = ["Generic scheme; substituent variables and exact atom mapping are unresolved."]
bundle["reactions"].append(r)


# Table 1, one experiment per row.
table = [
    (1, "1a", "2a", "3a", "4a", "5a", "0.3:1", 29, 28, 90),
    (2, "1a", "2b", "3a", "4b", "5b", "0.2:1", 47, 54, 87),
    (3, "1a", "2c", "3a", "4c", "5c", "2.4:1", 37, 18, 10),
    (4, "1a", "2d", "3a", "4d", None, ">20:1", 22, 16, None),
    (5, "1a", "2e", "3a", "4e", None, ">20:1", 72, 22, None),
    (6, "1a", "2f", "3a", "4f", None, ">20:1", 88, 96, None),
    (7, "1b", "2f", "3a", "4g", None, ">20:1", 96, 95, None),
    (8, "1b", "2f", "3b", "4g", None, ">20:1", 98, 87, None),
    (9, "1b", "2f", "3c", "4g", None, ">20:1", 99, 92, None),
    (10, "1b", "2f", "3d", "4g", None, ">20:1", 70, 74, None),
    (11, "1b", "2f", "3e", "4g", None, ">20:1", 92, 80, None),
    (12, "1b", "2f", "3f", "4g", None, ">20:1", 99, 96, None),
    (13, "1b", "2f", "3g", "4g", None, ">20:1", 98, 99, None),
]
for entry, aza, allene, cat, p4, p5, ratio, yld, ee4, ee5 in table:
    products = [p4] + ([p5] if p5 else [])
    r = base_reaction(f"table1_entry{entry:02d}", (1, 2) if entry == 1 else (2,),
                      "specific_experiment", ["optimization"], f"Table 1 entry {entry}")
    set_conditions(r, (2,), 0.05, 25, catalyst=cat, substrate1=aza, substrate2=allene, products=products)
    y_targets = products if entry <= 3 else [p4]
    add_outcome(r, "yield", "yield", y_targets, (2,), f"{yld}%" + (" combined" if entry <= 3 else ""),
                {"number": yld, "unit": "%"}, "combined cycloadducts" if entry <= 3 else "main [8+2] product",
                "combined yield for 4 and 5" if entry <= 3 else "isolated pure product",
                None if entry <= 3 else "isolated pure product")
    ratio_norm = {"order": [p4, p5 or "5"], "ratio": ratio}
    add_outcome(r, "ratio_4_5", "selectivity", products, (2,), ratio, ratio_norm,
                basis="1H NMR analysis of crude reaction mixture", scope="4/5 product ratio")
    add_outcome(r, "ee_4", "selectivity", [p4], (2,), f"{ee4}% ee", {"number": ee4, "unit": "%"},
                basis="HPLC on a chiral stationary phase", scope=f"enantiomeric excess of {p4}")
    if ee5 is not None:
        add_outcome(r, "ee_5", "selectivity", [p5], (2,), f"{ee5}% ee", {"number": ee5, "unit": "%"},
                    scope=f"enantiomeric excess of {p5}")
        r["review"]["unresolved"].append(f"The analytical method for the footnoted {p5} ee is not separately stated; Table footnote d explicitly applies HPLC to product 4.")
    if entry == 1:
        r["notes"].append(f("The structure and absolute configuration of 5a (90% ee) were determined by X-ray diffraction analysis.",
                            "5a absolute configuration established by X-ray", (2,)))
    if entry <= 3:
        r["review"]["unresolved"].append("The combined 4+5 yield measurement/isolation basis is not explicitly stated by footnote e.")
    r["review"]["unresolved"].append("Reaction time, atmosphere, workup, and purification are not reported in the supplied main-text table.")
    classify_cycloaddition(r, (2,))
    bundle["reactions"].append(r)


r = base_reaction("scheme2_generic", (2,), "generic_scheme", [], "Scheme 2 generic scope arrow")
r["notes"] = [f("Reactions of azaheptafulvenes 1a–d and allenic amides 2f–r to cycloadducts 4f–u.",
                  "generic substrate-scope transformation", (2,))]
r["review"]["unresolved"] = ["Generic variable structure; concrete rows are represented separately."]
bundle["reactions"].append(r)


scope = [
    ("4f", "1a", "2f", 25, (2, 5), 92, ">99"),
    ("4g", "1b", "2f", 25, (2, 5), 98, "99"),
    ("4h", "1c", "2f", 25, (2, 5), 98, "97"),
    ("4i", "1d", "2f", 25, (2, 5), 98, "93"),
    ("4j", "1b", "2g", 50, (18, 36), 87, ">99"),
    ("4k", "1b", "2h", 50, (18, 36), 71, ">99"),
    ("4l", "1b", "2i", 50, (18, 36), 76, ">99"),
    ("4m", "1b", "2j", 50, (18, 36), 72, ">99"),
    ("4n", "1b", "2k", 50, (18, 36), 93, ">99"),
    ("4o", "1b", "2l", 50, (18, 36), 71, "99"),
    ("4p", "1b", "2m", 50, (18, 36), 86, "99"),
    ("4q", "1b", "2n", 50, (18, 36), 55, "99"),
    ("4r", "1b", "2o", 50, (18, 36), 72, ">99"),
    ("4s", "1b", "2p", 50, None, 69, "96"),
    ("4t", "1b", "2q", 50, None, 61, "98"),
    ("4u", "1b", "2r", 50, 96, "<5", "n.d."),
]
for product, aza, allene, temp, time, yld, ee in scope:
    pages = (2, 3) if product in {"4j", "4k", "4s", "4t", "4u"} else (2,)
    r = base_reaction(f"scheme2_{product}", pages, "specific_experiment", ["substrate_scope"], f"Scheme 2 {product}")
    set_conditions(r, pages, 0.1, temp, time=time, substrate1=aza, substrate2=allene, products=(product,))
    if yld == "<5":
        yval = {"number": 5, "unit": "%", "qualifier": "lt", "raw_text": "<5%"}
    else:
        yval = {"number": yld, "unit": "%"}
    add_outcome(r, "yield", "yield", [product], (2,), f"{yld}%", yval, "[8+2] cycloadduct",
                "isolated pure product", "isolated pure product")
    if ee == "n.d.":
        add_outcome(r, "ee", "selectivity", [product], (2,), "n.d.", "not determined",
                    scope=f"enantiomeric excess of {product}")
    else:
        eval_ = {"number": int(ee.replace(">", "")), "unit": "%"}
        if ee.startswith(">"):
            eval_["qualifier"] = "gt"
            eval_["raw_text"] = f"{ee}% ee"
        add_outcome(r, "ee", "selectivity", [product], (2,), f"{ee}% ee", eval_,
                    basis="HPLC on a chiral stationary phase", scope=f"enantiomeric excess of {product}")
    r["review"]["unresolved"].append("Exact per-row reaction time is unavailable where only a group-level range is stated; detailed procedure is in unsupplied SI.")
    if product == "4u":
        r["notes"].append(f("Scheme 2: <5%, n.d.; prose: did not afford the corresponding cycloadduct after 96 h at 50 or 80 °C.",
                            "Scheme/prose descriptions retained together", (2, 3)))
        r["review"]["unresolved"].append("The relationship between the scheme's <5% entry and prose 'did not afford' wording is not resolved; no value is rewritten.")
    classify_cycloaddition(r, pages)
    bundle["reactions"].append(r)


# Parenthetical scale-up for 4j is a separate experiment, not a second yield on the 0.1 mmol row.
r = base_reaction("scheme2_4j_scaleup", (2, 3), "specific_experiment",
                  ["substrate_scope", "scale_up"], "Scheme 2 4j parenthetical scale-up")
# The parenthetical result and prose identify 4j (normally from 2g), while caption
# footnote c literally prints "1b and 2f". Retain both source alternatives without
# silently repairing the apparent substrate-label conflict.
set_conditions(r, (2, 3), 1.0, 50, substrate1="1b", substrate2="2g", products=("4j",))
allene_2g = next(p for p in r["participants"] if p["participant_id"] == "allene")
allene_2g["role"] = "unknown"
allene_2g["charges"]["amount"] = f()
participant(r, "caption_allene", "2f", "unknown", (2,))
r["reaction_scale"] = f(["4j scale-up at 1.0 mmol", "caption: 1b and 2f"], None, (2, 3), status="conflict")
r["reaction_scale"]["alternatives"] = [
    inferred("4j evaluated at a higher 1.0 mmol scale", {"reference_molecule_id": M["2g"], "amount": {"number": 1.0, "unit": "mmol"}}, (3,),
             "4j is the product mapped to allenic amide 2g in Scheme 2."),
    f("Reaction carried out on 1.0 mmol scale of 1b and 2f", {"reference_molecule_id": M["2f"], "amount": {"number": 1.0, "unit": "mmol"}}, (2,)),
]
add_outcome(r, "yield", "yield", ["4j"], (2,), "73%", {"number": 73, "unit": "%"},
            "[8+2] cycloadduct", "isolated pure product", "isolated pure product")
add_outcome(r, "ee", "selectivity", ["4j"], (2,), "99% ee", {"number": 99, "unit": "%"},
            basis="HPLC on a chiral stationary phase", scope="enantiomeric excess of 4j")
r["notes"].append(f("Scale-up substrate-label conflict: product placement/prose identify 4j (mapped to 2g), but Scheme 2 footnote c prints 1b and 2f.",
                    "Both source alternatives retained; no correction chosen.", (2, 3)))
r["review"]["unresolved"].extend([
    "Scale-up reaction time is not individually printed in the main-text scheme.",
    "Scale-up substrate identity conflicts between the 4j placement/prose and caption footnote c (2f); both alternatives are retained."
])
classify_cycloaddition(r, (2, 3))
bundle["reactions"].append(r)


# The prose explicitly reports a second failed temperature for bulky 2r.
r = base_reaction("control_2r_80C", (3,), "specific_experiment", ["control"],
                  "Page 3 bulky γ-isopropyl allene at 80 °C")
set_conditions(r, (3,), 0.1, 80, time=96, substrate1="1b", substrate2="2r", products=())
add_outcome(r, "no_cycloadduct", "no_reaction", [], (3,),
            "did not afford the corresponding cycloadduct after 96 h at 80 °C", True,
            scope="corresponding [8+2] cycloadduct not obtained")
classify_not_applicable(r, (3,), "No observed cycloadduct at 80 °C; no expected bond formation is assigned.")
r["review"]["unresolved"].append("No conversion, recovery, decomposition basis, or analytical detection method is reported.")
bundle["reactions"].append(r)


# Figure 1 is an author-proposed stereochemical model, not an experimental reaction.
r = base_reaction("figure1_stereochemical_model", (3,), "proposed_mechanism", [], "Figure 1 stereochemical model")
r["reaction_name"] = f("tentative stereochemical model", "proposed stereochemical model", (3,))
r["transformation_type"] = f("allyl phosphonium ylide arrangement with catalyst-substrate H-bond", "proposed transition arrangement", (3,))
r["notes"] = [f("An H-bond is proposed between the catalyst amide group and the azaheptafulvene sulfonyl group.",
                  "author-proposed H-bond model", (3,))]
r["review"]["unresolved"] = ["Tentative model only; not treated as experimental proof or an atom-mapped reaction."]
bundle["reactions"].append(r)


# Footnote 12: two real control experiments.
r = base_reaction("footnote12_3g_decomposition", (2, 5), "specific_experiment", ["control"],
                  "Footnote 12, ethyl hexa-2,3-dienoate under optimized conditions")
set_conditions(r, (2, 5), 0.1, 25, catalyst="3g", substrate1="1a", substrate2="ethyl_hexa", products=())
add_outcome(r, "decomposition", "decomposition", ["ethyl_hexa"], (5,), "only decomposition products were obtained",
            "decomposition only", scope="ethyl hexa-2,3-dienoate reaction outcome")
classify_not_applicable(r, (5,), "Only decomposition products were reported; no defined product connectivity is available.")
r["review"]["unresolved"].append("Footnote 12 does not state time, conversion, recovery, or a measurement method; scale/solvent/temperature/loading are inherited from the explicitly referenced optimized conditions.")
bundle["reactions"].append(r)

r = base_reaction("footnote12_nBu3P_trace5", (5,), "specific_experiment", ["control"],
                  "Footnote 12, nBu3P catalysis")
# Only substrate identity and catalyst replacement are directly safe; other conditions remain missing.
participant(r, "aza", "1a", "reactant", (5,))
participant(r, "allene", "ethyl_hexa", "reactant", (5,))
participant(r, "cat", "nBu3P", "catalyst", (5,))
participant(r, "product", "5_ethyl_hexa", "product", (5,))
step = st.new_step("reaction", [E[5], T[5]])
step["participant_ids"] = ["aza", "allene", "cat"]
step["operation"] = f("reaction under nBu3P catalysis", "react under nBu3P catalysis", (5,))
r["steps"] = [step]
add_outcome(r, "trace_product", "trace", ["5_ethyl_hexa"], (5,),
            "minor amounts of [8+4] cycloaddition product of type 5", "minor amounts",
            product_role="minor [8+4] product", basis="detected in the reaction mixture by NMR",
            scope="[8+4] cycloaddition product of type 5")
classify_cycloaddition(r, (5,), "Product type 5 denotes an annulation with two new C–C ring-closing connections; exact atom mapping is unavailable.")
r["review"]["unresolved"].extend([
    "nBu3P loading, solvent, concentration, temperature, time, conversion, and exact product structure are unreported in the supplied footnote.",
    "'Minor amounts' is retained qualitatively and is not converted to a numerical yield."
])
bundle["reactions"].append(r)


bundle["source_context"] = [
    {"context_id": f"{DOC}:context:abstract_graphical", "document_id": DOC, "pages": [1],
     "evidence_ids": [E[1], T[1]], "kind": "overview",
     "content": "Abstract and graphical abstract summarize racemic γ-substituted allenes, chiral amidophosphine catalysis, and high regio-/diastereo-/enantioselectivity; not duplicated as experiments."},
    {"context_id": f"{DOC}:context:group_times", "document_id": DOC, "pages": [2, 3],
     "evidence_ids": [E[2], T[2], E[3], T[3]], "kind": "scope_context",
     "content": "4f–i are described collectively as completing after 2–5 h; aryl-amide products 4j–r are described collectively as 18–36 h at 50 °C. Exact per-row times are not supplied."},
    {"context_id": f"{DOC}:context:xray_assignments", "document_id": DOC, "pages": [2, 3],
     "evidence_ids": [E[2], T[2], E[3], T[3]], "kind": "stereochemical_assignment",
     "content": "X-ray diffraction confirms 4g and 5a; 4k and 4s are also determined by X-ray. Other 4 products are assigned assuming identical stereochemical outcome."},
    {"context_id": f"{DOC}:context:figure1", "document_id": DOC, "pages": [3],
     "evidence_ids": [E[3], T[3]], "kind": "proposed_model",
     "content": "Figure 1 is a tentative H-bonded stereochemical model, not a measured synthetic experiment."},
    {"context_id": f"{DOC}:context:si_notice", "document_id": DOC, "pages": [3],
     "evidence_ids": [E[3], T[3]], "kind": "availability_notice",
     "content": "The paper states that SI contains complete procedures and characterization, but no SI file was explicitly supplied or associated in this session."},
    {"context_id": f"{DOC}:context:bibliography", "document_id": DOC, "pages": [4, 5],
     "evidence_ids": [E[4], T[4], E[5], T[5]], "kind": "bibliography",
     "content": "Page 4 and most of page 5 are references; footnote 12 on page 5 is separately represented by two control experiments."},
]


st.write(RUN / "submitted.json", bundle)
st.write(RUN / "raw_model/stereo_checks.json", stereo_checks)

page_records = {p: [] for p in range(1, 6)}
for reaction in bundle["reactions"]:
    pages = reaction["location"]["normalized_value"].get("pages", [])
    for page in pages:
        if page in page_records:
            page_records[page].append(reaction["reaction_id"])

raw = {
    "document_id": DOC,
    "pages_attempted": [1, 2, 3, 4, 5],
    "machine_records_remain_needs_review": True,
    "molecule_count": len(bundle["molecules"]),
    "reaction_count": len(bundle["reactions"]),
    "page_record_ids": page_records,
    "limitations": [
        "No SI was supplied or associated.",
        "Dense products and catalysts are retained by labels/drawings rather than unsafe complete SMILES.",
        "Candidate allene SMILES and bond classifications are machine-derived and require expert review.",
        "Schema validation does not establish chemical accuracy or semantic completeness."
    ],
    "bundle": copy.deepcopy(bundle),
}
st.write(RUN / "raw_model/pages_0001_0005_extraction.json", raw)

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title/abstract/graphical abstract", "finding": "Article metadata and overview only; no standalone experiment duplicated.", "evidence_ids": [E[1], T[1]]},
            {"region": "Scheme 1 eqs 1–3", "finding": "Three cited prior-work generic arrows, recorded as literature examples.", "evidence_ids": [E[1], T[1]]},
            {"region": "Scheme 1 eq 4", "finding": "Generic current-work [8+2] transformation, not a concrete experiment.", "evidence_ids": [E[1], T[1]]},
            {"region": "bottom prose crossing to page 2", "finding": "Table 1 entry 1 setup and products 4a/5a; concrete values are on page 2 and the record cites both pages.", "evidence_ids": [E[1], T[1], E[2], T[2]]}
        ], "record_ids": page_records[1], "unresolved": ["No SI procedure is supplied; Scheme 1 prior-work substituent details are generic."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Table 1 structure header and catalyst panel", "finding": "Azaheptafulvenes 1a/b, allenes 2a–f, catalysts 3a–g, and product branches 4/5 visually inspected.", "evidence_ids": [E[2], T[2]]},
            {"region": "Table 1 entries 1–13 and footnotes a–h", "finding": "Thirteen independent optimization experiments with yields, 4/5 ratios, ee values, and measurement footnotes.", "evidence_ids": [E[2], T[2]]},
            {"region": "Scheme 2 generic arrow and products 4f–u", "finding": "Sixteen standard scope rows plus one parenthetical 1.0 mmol scale-up; 4u retains <5%/n.d. without forcing a numerical zero.", "evidence_ids": [E[2], T[2]]},
            {"region": "Scheme 2 footnotes and surrounding prose", "finding": "0.1 mmol, 10 mol% 3g, toluene 0.05 M, 25/50 °C, isolated-yield basis, chiral-HPLC ee basis, group time ranges, and scale-up caveat captured.", "evidence_ids": [E[2], T[2]]}
        ], "record_ids": page_records[2], "unresolved": ["Exact per-row times and SI procedures are unavailable; dense product/catalyst complete structures need expert review."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "scope discussion", "finding": "4j scale-up, 4s/4t success, 2r failure after 96 h at 50 or 80 °C, and X-ray/analogy stereochemical assignments.", "evidence_ids": [E[3], T[3]]},
            {"region": "Figure 1 and mechanism discussion", "finding": "Tentative H-bonded stereochemical model recorded as proposed mechanism, not experimental proof.", "evidence_ids": [E[3], T[3]]},
            {"region": "Associated Content", "finding": "SI availability stated by the article, but no SI file supplied; association remains unestablished.", "evidence_ids": [E[3], T[3]]},
            {"region": "conclusion/authors/references start", "finding": "Context and bibliography; no additional concrete experiment.", "evidence_ids": [E[3], T[3]]}
        ], "record_ids": page_records[3], "unresolved": ["The prose 'did not afford' and Scheme 2 '<5%' descriptions for 4u are both retained without choosing an interpretation."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "two-column reference list", "finding": "Bibliography only; cited reactions are not imported as current-paper experiments.", "evidence_ids": [E[4], T[4]]}
        ], "record_ids": page_records[4], "unresolved": []},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "reference 11 continuation", "finding": "Bibliographic continuation only.", "evidence_ids": [E[5], T[5]]},
            {"region": "footnote 12", "finding": "Two controls: optimized conditions give decomposition only; nBu3P gives minor type-5 [8+4] product detected by NMR.", "evidence_ids": [E[5], T[5]]},
            {"region": "reference 13", "finding": "Cited computational literature only, not a current-paper experiment.", "evidence_ids": [E[5], T[5]]}
        ], "record_ids": page_records[5], "unresolved": ["Footnote 12 does not report nBu3P loading or full conditions, and the minor product's exact structure is absent."]}
    ]
}
st.write(RUN / "coverage/pages_0001_0005.json", coverage)

tasks = json.loads((RUN / "tasks.json").read_text())
for task in tasks:
    if task.get("document_id") == DOC and task.get("evidence_id") in T.values():
        task["state"] = "attempted"
        task["reason"] = "Original page image and extracted text inspected; see coverage/pages_0001_0005.json."
st.write(RUN / "tasks.json", tasks)

# Literal mention search is a cross-region review aid only; it does not merge experiments.
mentions = st.source_mentions(json.loads((RUN / "dataset.json").read_text()),
                              ["1a", "1b", "2f", "2r", "4g", "4j", "4u", "5a"], context=220)
st.write(RUN / "raw_model/source_mentions.json", mentions)
import reaction_tools as legacy_tools
current_dataset = json.loads((RUN / "dataset.json").read_text())
if "standard" in current_dataset:
    st.write(RUN / "raw_model/audit.json", legacy_tools.dataset_audit(current_dataset))

print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
                  "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in bundle["reactions"])},
                 ensure_ascii=False, indent=2))
