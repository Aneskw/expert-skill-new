#!/usr/bin/env python3
"""Build the evidence-linked v3.1 submission for Titov et al., Synlett 2020.

This is a deterministic transcription of the page 1--5 visual/source review.
It deliberately leaves complex structures unencoded when a safe isomeric
SMILES was not established from the supplied paper images.
"""

from __future__ import annotations

import copy
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st  # noqa: E402


DOC = "doc_776101937d954ae0"
RUN = ROOT / "runs/current"
E = {n: f"{DOC}:page{n}:text" for n in range(1, 6)}
I = {n: f"{DOC}:page{n}:image" for n in range(1, 6)}


def ev(*pages: int) -> list[str]:
    values: list[str] = []
    for page in pages:
        values.extend([E[page], I[page]])
    return values


def fact(raw=None, normalized=None, evidence=(), method="original_extraction", derivation=None,
         status=None):
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation, status=status)


def derived(raw, normalized, evidence, derivation):
    return fact(raw, normalized, evidence, method="structure_derivation", derivation=derivation)


def inferred(raw, normalized, evidence, derivation):
    return fact(raw, normalized, evidence, method="inference", derivation=derivation)


def label_fact(label: str, evidence: list[str]):
    return fact(label, label, evidence)


def set_stereo(m, evidence, *, chirality="unresolved", geometry="unresolved",
               representation="not_encoded", axial=None, assignment=None, mixture=None,
               double_bonds=None, note=None):
    s = m["stereochemistry"]
    s["chirality_status"] = inferred(
        chirality, chirality, evidence,
        "Status records only what the supplied drawing/text specifies; missing wedges or names are not treated as achiral or racemic.",
    )
    s["geometry_status"] = inferred(
        geometry, geometry, evidence,
        "E/Z and ring geometry were reviewed separately from central/axial chirality; no unreported assignment was invented.",
    )
    s["representation_status"] = inferred(
        representation, representation, evidence,
        "No complete isomeric SMILES was asserted for this source drawing unless separately provided.",
    )
    if axial is not None:
        s["axial"] = fact(axial, axial, evidence)
    if assignment is not None:
        s["assignment"] = fact(assignment, assignment, evidence)
    if mixture is not None:
        s["mixture"] = fact(mixture, mixture, evidence)
    if double_bonds is not None:
        s["double_bonds"] = fact(double_bonds, double_bonds, evidence)
    if note:
        m["identity_note"] = fact(note, note, evidence)


def unresolved_molecule(mid: str, label: str, evidence: list[str], *, name=None,
                        chirality="unresolved", geometry="unresolved",
                        representation="not_encoded", axial=None, assignment=None,
                        mixture=None, double_bonds=None, note=None):
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["representation_type"] = "unresolved"
    m["occurrences"][0]["compound_label"] = label_fact(label, evidence)
    m["names"]["original"] = fact(name or label, name or label, evidence)
    set_stereo(
        m, evidence, chirality=chirality, geometry=geometry,
        representation=representation, axial=axial, assignment=assignment,
        mixture=mixture, double_bonds=double_bonds,
        note=note or "Source identity is retained, but no safe complete SMILES was assigned from the dense drawing.",
    )
    return m


def smiles_molecule(mid: str, label: str, smiles: str, evidence: list[str], *, name=None):
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["representation_type"] = "smiles"
    m["occurrences"][0]["compound_label"] = label_fact(label, evidence)
    m["names"]["original"] = fact(name or label, name or label, evidence)
    diag = st.stereo_inspect(smiles)
    derivation = (
        f"Candidate connectivity for the explicit source name/abbreviation {label!r}; "
        "RDKit graph diagnostics check syntax and potential stereo but do not certify source identity."
    )
    m["smiles"] = derived(label, smiles, evidence, derivation)
    m["stereochemistry"]["chirality_status"] = inferred(
        label, "achiral", evidence,
        "The candidate graph has no tetrahedral or axial stereogenic element; this is a structure-derived status.",
    )
    m["stereochemistry"]["geometry_status"] = inferred(
        label, "not_applicable", evidence,
        "No source-relevant E/Z or ring cis/trans assignment is present in the candidate graph.",
    )
    m["stereochemistry"]["representation_status"] = inferred(
        smiles, "encoded", evidence,
        "The candidate graph contains no unresolved stereogenic element according to the saved graph diagnostic.",
    )
    m["stereochemistry"]["isomeric_smiles"] = derived(
        smiles, diag["canonical_isomeric_smiles"], evidence, derivation
    )
    m["stereochemistry"]["identity_key"] = derived(
        diag["canonical_isomeric_smiles"], diag["stereo_graph_key"], evidence,
        "SHA-256 stereo graph key calculated by the bundled standard stereo helper; not a sample identity.",
    )
    m["identity_note"] = fact(
        label,
        "Explicitly named simple participant; candidate SMILES is structure-derived and remains machine-review status.",
        evidence,
    )
    stereo_checks[mid] = diag
    return m


def template_molecule(mid: str, label: str, evidence: list[str], *, fixed, variables,
                      chirality="unresolved", geometry="unresolved", axial=None,
                      assignment=None, note=None):
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["representation_type"] = "generic_template"
    m["occurrences"][0]["compound_label"] = label_fact(label, evidence)
    m["names"]["original"] = fact(label, label, evidence)
    var_rows = []
    for var in variables:
        var_rows.append({
            "label": var["label"],
            "attachment_points": fact(var.get("attachment"), var.get("attachment"), evidence),
            "definition": fact(var.get("definition"), var.get("definition"), evidence),
            "allowed_values": fact(var.get("allowed"), var.get("allowed"), evidence),
            "constraints": fact(var.get("constraints"), var.get("constraints"), evidence),
        })
    m["generic_template"] = {
        "fixed_structure": fact(fixed, fixed, evidence),
        "variables": var_rows,
        "instance_mapping": [],
        "reason": fact(
            "Generic R-group/ring-size drawing in the source",
            "Retained as a template because the source does not define one concrete molecule.",
            evidence,
        ),
    }
    set_stereo(
        m, evidence, chirality=chirality, geometry=geometry,
        representation="not_encoded", axial=axial, assignment=assignment,
        note=note or "Generic source template; a wildcard SMILES was not used as a concrete identity.",
    )
    return m


def add_name_occurrence(m, label, evidence):
    existing = set(m["occurrences"][0]["evidence_ids"])
    for item in evidence:
        if item not in existing:
            m["occurrences"][0]["evidence_ids"].append(item)


def participant(pid, mid, role, evidence, *, stage="reaction", amount=None, volume=None,
                mass=None, product_role=None):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    if amount is not None:
        p["charges"]["amount"] = fact(amount["raw"], amount["value"], evidence)
    if volume is not None:
        p["charges"]["volume"] = fact(volume["raw"], volume["value"], evidence)
    if mass is not None:
        p["charges"]["mass"] = fact(mass["raw"], mass["value"], evidence)
    if product_role is not None:
        p["product_role"] = fact(product_role, product_role, evidence)
    return p


def step(step_id, stage, pids, evidence, *, operation=None, temperature=None, time=None):
    s = st.new_step(step_id, evidence, stage=stage)
    s["participant_ids"] = list(pids)
    if operation is not None:
        s["operation"] = fact(operation, operation, evidence)
    if temperature is not None:
        s["temperature"] = fact(temperature["raw"], temperature["value"], evidence)
    if time is not None:
        s["time"] = fact(time["raw"], time["value"], evidence)
    return s


def outcome(oid, kind, targets, evidence, *, raw, value, product_role=None,
            yield_type=None, basis=None, scope=None):
    o = st.new_outcome(oid, kind, targets, evidence)
    o["value"] = fact(raw, value, evidence)
    if product_role is not None:
        o["product_role"] = fact(product_role, product_role, evidence)
    if yield_type is not None:
        o["yield_type"] = fact(yield_type, yield_type, evidence)
    if basis is not None:
        o["measurement_basis"] = fact(basis, basis, evidence)
    if scope is not None:
        o["scope"] = fact(scope, scope, evidence)
    return o


def set_location(r, raw, evidence, normalized):
    r["location"] = fact(raw, normalized, evidence)


def set_identity(r, categories, evidence, reaction_name, transformation):
    r["experiment_categories"] = fact(categories, categories, evidence)
    r["reaction_name"] = fact(reaction_name, reaction_name, evidence)
    r["transformation_type"] = inferred(
        transformation, transformation, evidence,
        "Transformation description inferred from the source arrow and accompanying prose.",
    )


def set_classification(r, evidence, changes, *, coverage="partial", limitations=()):
    normalized_changes = []
    labels = {}
    for index, change in enumerate(changes, 1):
        item = copy.deepcopy(change)
        item.setdefault("change_id", f"b{index}")
        normalized_changes.append(item)
        if item["change_type"] in ("formed", "hydrogen_gain", "order_increased"):
            lab = st.bond_label(item["atoms"][0]["element"], item["atoms"][1]["element"], item["after_order"])
            labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(item["change_id"])
    derivation = (
        "Partial source-drawing comparison only; no complete atom mapping was asserted. "
        "The listed changes capture only unambiguous net bond classes visible in the reaction scheme/mechanistic description."
    )
    r["classification"]["bond_changes"] = [
        derived("source-visible net bond change", item, evidence, derivation)
        for item in normalized_changes
    ]
    r["classification"]["labels"] = derived(
        "formed/promoted bonds from the listed partial changes", list(labels.values()), evidence, derivation
    )
    r["classification"]["coverage"] = inferred(coverage, coverage, evidence, derivation)
    r["classification"]["limitations"] = list(limitations) or [
        "Complex source drawings were not assigned a complete atom mapping.",
        "Only source-visible net bond classes are listed; omitted bond-order shifts or hydrogen partners may remain.",
    ]
    r["classification"]["model_votes"] = []


def unresolved_classification(r, evidence, reason):
    r["classification"]["labels"] = inferred([], [], evidence, reason)
    r["classification"]["coverage"] = inferred("unresolved", "unresolved", evidence, reason)
    r["classification"]["limitations"] = [reason]
    r["classification"]["model_votes"] = []


def mark_review(r, unresolved):
    r["review"] = {
        "status": "needs_review",
        "method": "model",
        "reviewer": None,
        "checked_fields": [],
        "unresolved": list(unresolved),
    }
    r["completeness"] = "partial"


def add_common_microwave(r, evidence, label, scale_mmol, time_h_or_min, tlc_ratio,
                         *, procedure, dry_word="dry"):
    sub_mid = f"{DOC}:mol_{label}"
    p_sub = participant(
        "substrate", sub_mid, "reactant", evidence,
        amount={"raw": f"{scale_mmol} mmol", "value": {"number": scale_mmol, "unit": "mmol"}},
    )
    p_tol = participant(
        "toluene", f"{DOC}:mol_toluene", "solvent", evidence,
        volume={"raw": "6.0 mL", "value": {"number": 6.0, "unit": "mL"}},
    )
    r["participants"].extend([p_sub, p_tol])
    r["reaction_scale"] = fact(
        f"{label} ({scale_mmol} mmol)",
        {"reference_molecule_id": sub_mid, "amount": {"number": scale_mmol, "unit": "mmol"}},
        evidence,
    )
    r["solvent_conditions"]["composition"] = fact(
        f"{dry_word} toluene", [{"molecule_id": f"{DOC}:mol_toluene", "descriptor": dry_word}], evidence
    )
    r["solvent_conditions"]["total_volume"] = fact(
        "6.0 mL", {"number": 6.0, "unit": "mL"}, evidence
    )
    if isinstance(time_h_or_min, tuple):
        time_obj = {"raw": time_h_or_min[0], "value": time_h_or_min[1]}
    else:
        time_obj = {"raw": f"{time_h_or_min} h", "value": {"number": time_h_or_min, "unit": "h"}}
    r["steps"].append(step(
        "reaction", "reaction", ["substrate", "toluene"], evidence,
        operation="Microwave irradiation with stirring; progress monitored by TLC",
        temperature={"raw": "150 °C", "value": {"number": 150, "unit": "°C"}},
        time=time_obj,
    ))
    r["special_conditions"]["reactor"] = fact("microwave reactor", "microwave reactor", evidence)
    r["procedure"]["original_text"] = fact(procedure, procedure, evidence)
    r["procedure"]["general_procedure_reference"] = fact(
        procedure.split(".", 1)[0], procedure.split(".", 1)[0], evidence
    )
    r["procedure"]["addition_order"] = fact(
        f"A solution of the appropriate substrate in {dry_word} toluene was placed in the microwave reactor.",
        ["dissolve substrate in toluene", "place solution in microwave reactor"],
        evidence,
    )
    r["procedure"]["pretreatment"] = fact(
        f"{dry_word} toluene", f"toluene specified as {dry_word}", evidence
    )
    r["participants"].extend([
        participant("tlc_etoac", f"{DOC}:mol_ethyl_acetate", "unknown", evidence, stage="characterization"),
        participant("tlc_hexane", f"{DOC}:mol_hexane", "unknown", evidence, stage="characterization"),
    ])
    r["steps"].append(step(
        "tlc_monitoring", "characterization", ["tlc_etoac", "tlc_hexane"], evidence,
        operation=f"TLC monitoring with EtOAc–hexane {tlc_ratio}",
    ))


stereo_checks = {}
bundle = json.loads((RUN / "submitted.json").read_text())
bundle["molecules"] = []
bundle["reactions"] = []
bundle["source_context"] = []


# Explicit simple participants. Candidate SMILES are checked and saved below.
simple = [
    ("mol_toluene", "toluene", "Cc1ccccc1", "toluene", ev(1, 2, 4, 5)),
    ("mol_ethyl_acetate", "EtOAc", "CCOC(=O)C", "ethyl acetate", ev(4, 5)),
    ("mol_hexane", "hexane", "CCCCCC", "n-hexane", ev(4, 5)),
    ("mol_diethyl_ether", "Et2O", "CCOCC", "diethyl ether", ev(4, 5)),
    ("mol_methyl_propiolate", "methyl propiolate", "COC(=O)C#C", "methyl propiolate", ev(1, 2)),
    ("mol_trifluoroethanol", "CF3CH2OH", "OCC(F)(F)F", "2,2,2-trifluoroethanol", ev(1, 2)),
    ("mol_dbu", "DBU", "C1CCC2=NCCCN2CC1", "1,8-diazabicyclo[5.4.0]undec-7-ene", ev(1, 2)),
    ("mol_thf", "THF", "C1CCOC1", "tetrahydrofuran", ev(1)),
    ("mol_acetonitrile", "MeCN", "CC#N", "acetonitrile", ev(1)),
    ("mol_benzene", "benzene", "c1ccccc1", "benzene", ev(2)),
    ("mol_chloroacetyl_fluoride", "ClCH2COF", "O=C(F)CCl", "chloroacetyl fluoride", ev(1)),
    ("mol_dichloromethane", "CH2Cl2", "ClCCl", "dichloromethane", ev(1)),
    ("mol_potassium_carbonate", "K2CO3", "[K+].[K+].[O-]C([O-])=O", "potassium carbonate", ev(1)),
    ("mol_cyclohexadiene", "1,2-cyclohexadiene", "C1CCC=C=C1", "1,2-cyclohexadiene", ev(2)),
]
for suffix, label, smiles, name, evidence in simple:
    bundle["molecules"].append(smiles_molecule(f"{DOC}:{suffix}", label, smiles, evidence, name=name))

# RDKit's graph diagnostic does not certify axial/allene stereochemistry.  The
# prior-work 1,2-cyclohexadiene drawing has no axial assignment, so its
# connectivity candidate must not be promoted to an achiral, fully encoded
# stereoisomer merely because the toolkit reports no tetrahedral center.
cyclohexadiene = next(m for m in bundle["molecules"] if m["molecule_id"] == f"{DOC}:mol_cyclohexadiene")
cyclohexadiene["stereochemistry"]["chirality_status"] = inferred(
    "No axial configuration is shown for 1,2-cyclohexadiene",
    "unspecified", ev(2),
    "The source allene has no R_a/S_a or P/M assignment; toolkit axial support is insufficient to call it achiral.",
)
cyclohexadiene["stereochemistry"]["axial"] = fact(
    "1,2-cyclohexadiene allene axis; configuration not specified",
    {"feature": "allene axis", "configuration": "unspecified"}, ev(2),
)
cyclohexadiene["stereochemistry"]["representation_status"] = inferred(
    "Connectivity SMILES does not encode the allene axis", "partial", ev(2),
    "C1CCC=C=C1 retains connectivity only; no axial stereochemical claim is made.",
)
cyclohexadiene["stereochemistry"]["isomeric_smiles"] = fact()
cyclohexadiene["stereochemistry"]["identity_key"] = fact()
cyclohexadiene["identity_note"] = fact(
    "1,2-cyclohexadiene source drawing",
    "Connectivity candidate only; axial configuration is unreported and not encoded.", ev(2),
)

# Named but composition/structure-unsafely encoded auxiliaries.
bundle["molecules"].extend([
    unresolved_molecule(f"{DOC}:mol_silica_gel", "silica gel", ev(4), name="silica gel",
                        chirality="not_applicable", geometry="not_applicable",
                        note="Chromatography stationary phase; composition is not represented as one molecular SMILES."),
    unresolved_molecule(f"{DOC}:mol_ether_unspecified", "ether", ev(2), name="ether (identity unspecified)",
                        chirality="unresolved", geometry="not_applicable",
                        note="Scheme 2 prints only 'ether'; its exact identity is not inferred."),
    unresolved_molecule(f"{DOC}:mol_trimethylaluminum", "Me3Al", ev(1), name="trimethylaluminum",
                        chirality="unresolved", geometry="not_applicable",
                        note="Reagent identity is explicit, but aggregation/solution form is not reported and no concrete registry SMILES is asserted."),
])

# Current-paper compounds 4a--e, 7a--e and 8a--e.
subs = {
    "a": ("OMe", "Me"), "b": ("OMe", "i-Pr"), "c": ("OMe", "Bn"),
    "d": ("OMe", "Ph"), "e": ("H", "Me"),
}
name4 = {
    "b": "Methyl 8-isopropyl-10,11-dimethoxy-3-methyl-6-phenyl-3-azabenzo[d]cyclodeca-4,6,7-triene-5-carboxylate",
    "d": "Methyl 10,11-dimethoxy-3-methyl-6,8-diphenyl-3-azabenzo[d]cyclodeca-4,6,7-triene-5-carboxylate",
}
name7a = "Methyl 6,7-dimethoxy-2,4b-dimethyl-4a-phenyl-1,2,4a,4b,8b,8c-hexahydro-2-azabenzo[a]cyclopropa[cd]azulene-4-carboxylate"
name8a = "Methyl 5,6-dimethoxy-2,11-dimethyl-14-phenyl-11-azatetracyclo[7.5.0.0^2,12.0^3,8]tetradeca-3,5,7,13-tetraene-13-carboxylate"
for letter, (r_value, r1_value) in subs.items():
    extra_pages = (4,) if letter in ("b", "d") else ()
    evidence4 = ev(1, 2, *extra_pages)
    axial = {
        "feature": "cyclic allene C=C=C shown in the source drawing",
        "assignment": "no R_a/S_a or P/M assignment reported",
        "substituents": {"R": r_value, "R1": r1_value},
    }
    bundle["molecules"].append(unresolved_molecule(
        f"{DOC}:mol_4{letter}", f"4{letter}", evidence4,
        name=name4.get(letter, f"azacyclic allene 4{letter} (R={r_value}, R1={r1_value})"),
        chirality="unspecified", geometry="unspecified", axial=axial,
        double_bonds="Source drawing contains cumulene/enamine unsaturation; no E/Z designation is reported.",
        note="Concrete source compound with table-defined R/R1 values; cyclic-allene axial configuration is not reported and no complete isomeric SMILES is assigned.",
    ))
    prod_ev = ev(2, 3) if letter == "a" else ev(2)
    assignment = {
        "source_statement": "thermolysis products 7 and 8 were formed as single diastereoisomers",
        "absolute_configuration": "not assigned",
        "xray": "7a and 8a confirmed by single-crystal X-ray" if letter == "a" else "not reported on supplied pages",
    }
    for family, full_name in (("7", name7a if letter == "a" else None), ("8", name8a if letter == "a" else None)):
        pages = list(prod_ev)
        if letter == "a":
            pages += ev(4, 5)
        bundle["molecules"].append(unresolved_molecule(
            f"{DOC}:mol_{family}{letter}", f"{family}{letter}", pages,
            name=full_name or (f"2-azabenzo[a]cyclopropa[cd]azulene {family}{letter}" if family == "7" else f"(epiminomethano)cyclopenta[a]indene {family}{letter}"),
            chirality="partially_specified", geometry="unspecified",
            assignment=assignment, mixture="single diastereoisomer; enantiomeric composition not reported",
            double_bonds="No E/Z assignment is named on the supplied pages.",
            note="Source reports a single diastereoisomer but provides no absolute R/S assignment; dense polycyclic stereochemistry is retained textually and not encoded as SMILES.",
        ))

# Generic templates for the current thermolysis sequence and literature arrows.
template_ev = ev(1, 2, 3)
vars_4 = [
    {"label": "R", "attachment": "two aromatic-ring positions", "definition": "aromatic substituent", "allowed": ["H", "OMe"]},
    {"label": "R1", "attachment": "allene terminus/ring carbon shown in Scheme/Table", "definition": "substituent", "allowed": ["Me", "i-Pr", "Bn", "Ph"]},
]
for family, fixed, chirality in [
    ("4", "azacyclic allene/enamine framework shown for 4", "unspecified"),
    ("7", "2-azabenzo[a]cyclopropa[cd]azulene framework shown for 7", "partially_specified"),
    ("8", "(epiminomethano)cyclopenta[a]indene framework shown for 8", "partially_specified"),
]:
    bundle["molecules"].append(template_molecule(
        f"{DOC}:mol_generic_{family}", f"generic {family}", template_ev,
        fixed=fixed, variables=vars_4, chirality=chirality, geometry="unspecified",
        axial={"feature": "allene axis", "assignment": "unassigned"} if family == "4" else None,
        assignment={"source_statement": "products 7 and 8 are formed as single diastereoisomers", "absolute_configuration": "unreported"} if family in ("7", "8") else None,
    ))

literature_templates = [
    ("lit_salt1", "6-ethynyl-1-carboxymethyl-1,2,5,6-tetrahydropyridinium salt template", ev(1), "pyridinium salt drawing", [{"label": "R", "definition": "source substituent", "allowed": None}]),
    ("lit_allene1", "azacyclononatriene 1 template", ev(1), "azacyclononatriene framework 1", [{"label": "R", "definition": "source substituent", "allowed": None}]),
    ("lit_salt2", "isothiochromenium triflate template", ev(1), "1-alkynyl isothiochromenium triflate drawing", [{"label": "R", "definition": "alkynyl substituent", "allowed": ["Me", "Bu", "Ph"]}]),
    ("lit_allene2", "sulfur-containing cyclic allene 2 template", ev(1, 2), "sulfur-containing cyclic allene framework 2", [{"label": "R", "definition": "substituent", "allowed": ["Me", "Bu", "Ph"]}]),
    ("lit_alkyne3", "2-(phenylethynyl)-N-methylpiperidine/azepine template", ev(1), "cyclic propargyl amine drawing", [{"label": "n", "definition": "ring-size variable", "allowed": [1, 2]}]),
    ("lit_lactam3", "allenyl lactam 3 template", ev(1), "3-chloro-1-aza-deca/undeca-4,5-diene-2-one framework", [{"label": "n", "definition": "ring-size variable", "allowed": [1, 2]}]),
    ("lit_isoquinoline", "1-substituted 1-(phenylethynyl)-1,2,3,4-tetrahydroisoquinoline template", ev(1, 2), "isoquinoline precursor drawing", vars_4),
    ("lit_electron_alkyne", "electron-deficient alkyne template", ev(1), "alkyne bearing X", [{"label": "X", "definition": "electron-withdrawing group", "allowed": ["Ac", "CO2Me"]}]),
    ("lit_5a", "thermal dimer 5a", ev(2), "path-a carbocyclic product drawing", []),
    ("lit_5b", "low-temperature dimer 5b", ev(2), "path-b polycyclic product drawing", []),
    ("lit_diene6", "1,3-diene 6 template", ev(2), "sulfur-containing 1,3-diene framework 6", [{"label": "R", "definition": "substituent", "allowed": ["Me", "Bu", "Ph"]}]),
]
for suffix, label, evidence, fixed, variables in literature_templates:
    bundle["molecules"].append(template_molecule(
        f"{DOC}:mol_{suffix}", label, evidence, fixed=fixed, variables=variables,
        chirality="unresolved", geometry="unspecified",
    ))

# Concrete but unencoded precursors for the two novel allenes.
for letter in ("b", "d"):
    bundle["molecules"].append(unresolved_molecule(
        f"{DOC}:mol_precursor_4{letter}", f"corresponding 1-(phenylethynyl)isoquinoline for 4{letter}", ev(1, 2),
        name=f"corresponding 1-(phenylethynyl)isoquinoline precursor to 4{letter}",
        chirality="unspecified", geometry="unspecified",
        note="The prose identifies the precursor class but does not print a full name or safe line notation for this specific member.",
    ))


def literature_reaction(suffix, pages, reactant_mid, product_mid, *, condition_parts,
                        temp=None, yield_raw=None, yield_value=None, name, transformation,
                        changes=()):
    evidence = ev(*pages)
    rid = f"{DOC}:{suffix}"
    r = st.new_reaction(rid, DOC, evidence)
    r["record_type"] = "literature_example"
    set_location(r, name, evidence, {"pages": list(pages), "source_type": "cited prior-work scheme"})
    set_identity(r, ["other"], evidence, name, transformation)
    r["participants"] = [
        participant("reactant", f"{DOC}:{reactant_mid}", "reactant", evidence),
        participant("product", f"{DOC}:{product_mid}", "product", evidence, product_role="reported literature product"),
    ]
    for index, (mid, role) in enumerate(condition_parts, 1):
        r["participants"].append(participant(f"agent_{index}", f"{DOC}:{mid}", role, evidence))
    pids = [p["participant_id"] for p in r["participants"]]
    r["steps"] = [step(
        "literature_arrow", "reaction", pids, evidence,
        operation="Conditions printed on the cited prior-work reaction arrow",
        temperature=temp,
    )]
    if yield_raw is not None:
        r["outcomes"] = [outcome(
            "reported_yield", "yield", [f"{DOC}:{product_mid}"], evidence,
            raw=yield_raw, value=yield_value, product_role="reported literature product",
            yield_type="reported yield in cited scheme", basis=None,
            scope="Prior publication result reproduced in this paper's contextual scheme.",
        )]
    if changes:
        set_classification(r, evidence, list(changes))
    else:
        unresolved_classification(r, evidence, "Prior-work drawing is retained, but a complete safe atom correspondence was not established from this contextual scheme.")
    mark_review(r, [
        "This is a contextual literature example, not a newly reported experiment in the present paper.",
        "Exact scale, time, workup, purification and measurement basis are not supplied on these pages.",
    ])
    bundle["reactions"].append(r)


literature_reaction(
    "scheme1_literature_allene1", (1,), "mol_lit_salt1", "mol_lit_allene1",
    condition_parts=[("mol_dbu", "reagent"), ("mol_thf", "solvent")],
    temp={"raw": "r.t.", "value": "r.t."},
    yield_raw="30–40%", yield_value={"minimum": 30, "maximum": 40, "unit": "%", "raw_text": "30–40%"},
    name="Sashida and co-workers' synthesis of allene 1", transformation="formation of a heterocyclic allene",
)
literature_reaction(
    "scheme1_literature_allene2", (1,), "mol_lit_salt2", "mol_lit_allene2",
    condition_parts=[("mol_dbu", "reagent"), ("mol_acetonitrile", "solvent")],
    yield_raw="60–78%", yield_value={"minimum": 60, "maximum": 78, "unit": "%", "raw_text": "60–78%"},
    name="Sashida and co-workers' synthesis of allene 2", transformation="formation of a sulfur-containing cyclic allene",
)
literature_reaction(
    "scheme1_literature_lactam3", (1,), "mol_lit_alkyne3", "mol_lit_lactam3",
    condition_parts=[("mol_chloroacetyl_fluoride", "reagent"), ("mol_trimethylaluminum", "reagent"),
                     ("mol_dichloromethane", "solvent"), ("mol_potassium_carbonate", "reagent")],
    temp={"raw": "0 °C to r.t.", "value": {"number": 0, "unit": "°C", "raw_text": "0 °C to r.t.", "transition_to": "r.t."}},
    yield_raw="70–77%", yield_value={"minimum": 70, "maximum": 77, "unit": "%", "raw_text": "70–77%"},
    name="Perscheid and co-workers' synthesis of medium-sized ring allenyl lactams 3", transformation="[3,3]-sigmatropic rearrangement to an allenyl lactam",
)
literature_reaction(
    "scheme1_previous_allene4", (1, 2), "mol_lit_isoquinoline", "mol_generic_4",
    condition_parts=[("mol_lit_electron_alkyne", "reactant"), ("mol_trifluoroethanol", "solvent")],
    temp={"raw": "+7 °C", "value": {"number": 7, "unit": "°C"}},
    yield_raw="80–90%", yield_value={"minimum": 80, "maximum": 90, "unit": "%", "raw_text": "80–90%"},
    name="Authors' previous synthesis of benzazecines 4 containing an allene moiety", transformation="formation of an azacyclic allene",
)
cc_change = ({"atoms": [{"element": "C", "source_locator": "first allene-derived carbon"}, {"element": "C", "source_locator": "second allene-derived carbon"}], "before_order": 0, "after_order": 1, "change_type": "formed"},)
literature_reaction(
    "scheme2_literature_path_a", (2,), "mol_cyclohexadiene", "mol_lit_5a",
    condition_parts=[("mol_ether_unspecified", "solvent")],
    temp={"raw": "reflux (path a)", "value": "reflux"},
    name="Mahlokozera and co-workers' path-a reaction of 1,2-cyclohexadiene", transformation="thermal dimerization/cycloaddition", changes=cc_change,
)
literature_reaction(
    "scheme2_literature_path_b", (2,), "mol_cyclohexadiene", "mol_lit_5b",
    condition_parts=[("mol_ether_unspecified", "solvent")],
    temp={"raw": "−80 °C (path b)", "value": {"number": -80, "unit": "°C"}},
    name="Mahlokozera and co-workers' path-b reaction of 1,2-cyclohexadiene", transformation="low-temperature dimerization/cycloaddition", changes=cc_change,
)
literature_reaction(
    "scheme2_literature_allene2_to_6", (2,), "mol_lit_allene2", "mol_lit_diene6",
    condition_parts=[("mol_dbu", "reagent"), ("mol_benzene", "solvent")],
    temp={"raw": "reflux", "value": "reflux"},
    yield_raw="45–70%", yield_value={"minimum": 45, "maximum": 70, "unit": "%", "raw_text": "45–70%"},
    name="Sashida and co-workers' allene isomerisation 2 to 6", transformation="thermal/base-promoted isomerization of a cyclic allene to a 1,3-diene",
)

# Generic current-paper thermolysis sequence.
generic_id = f"{DOC}:generic_thermolysis_4_to_7_8"
g_ev = ev(1, 2, 3)
g = st.new_reaction(generic_id, DOC, g_ev)
g["record_type"] = "generic_scheme"
set_location(g, "Graphical abstract, Table 1 generic arrow and Scheme 3", g_ev,
             {"pages": [1, 2, 3], "regions": ["graphical abstract", "Table 1", "Scheme 3"]})
set_identity(g, ["other"], g_ev, "Thermal transformation of azacyclic allenes 4 into cyclopropanes 7 and indenes 8",
             "microwave thermolysis/rearrangement")
g["participants"] = [
    participant("generic_4", f"{DOC}:mol_generic_4", "reactant", g_ev),
    participant("generic_7", f"{DOC}:mol_generic_7", "product", g_ev, product_role="primary thermolysis product/intermediate"),
    participant("generic_8", f"{DOC}:mol_generic_8", "product", g_ev, product_role="final thermolysis product"),
    participant("toluene", f"{DOC}:mol_toluene", "solvent", g_ev),
]
g["solvent_conditions"]["composition"] = fact("toluene", [{"molecule_id": f"{DOC}:mol_toluene"}], g_ev)
g["steps"] = [
    step("primary_thermolysis", "reaction", ["generic_4", "generic_7", "generic_8", "toluene"], g_ev,
         operation="Microwave-assisted thermolysis giving 7 and 8",
         temperature={"raw": "150 °C", "value": {"number": 150, "unit": "°C"}},
         time={"raw": "30 min", "value": {"number": 30, "unit": "min"}}),
    step("prolonged_thermolysis", "reaction", ["generic_4", "generic_7", "generic_8", "toluene"], g_ev,
         operation="Continued microwave heating to consume 4 and 7 and form 8",
         temperature={"raw": "150 °C", "value": {"number": 150, "unit": "°C"}},
         time={"raw": "2–7.5 hours (graphical abstract)", "value": {"minimum": 2, "maximum": 7.5, "unit": "h", "raw_text": "2–7.5 hours"}}),
]
g["special_conditions"]["reactor"] = fact("MW / microwave-assisted synthesis", "microwave reactor", g_ev)
g["notes"] = [fact(
    "Cyclopropanes 7 were established as primary thermolysis products/intermediates; further heating gave indenes 8.",
    "Author interpretation: 7 is an intermediate en route to 8.", g_ev
)]
set_classification(g, g_ev, [
    {"atoms": [{"element": "C", "source_locator": "allene-derived carbon"}, {"element": "C", "source_locator": "ring carbon"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
    {"atoms": [{"element": "C", "source_locator": "migrating-hydrogen acceptor"}, {"element": "H", "source_locator": "benzyl hydrogen"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain"},
])
mark_review(g, [
    "Generic R/R1 structures are templates; no complete atom mapping or stereochemical encoding is asserted.",
    "The graphical-abstract range 2–7.5 h summarizes routes; substrate-specific completion times are recorded separately from Table 1.",
])
bundle["reactions"].append(g)

# Two novel allene preparations reported by product yield/characterization, with incomplete procedures.
for letter, mass_mg in (("b", 438), ("d", 473)):
    evidence = ev(1, 2, 4)
    rid = f"{DOC}:synthesis_4{letter}"
    r = st.new_reaction(rid, DOC, evidence)
    r["record_type"] = "specific_experiment"
    set_location(r, f"Scheme 1/prose and Note {13 if letter == 'b' else 14}", evidence,
                 {"pages": [1, 2, 4], "compound": f"4{letter}"})
    set_identity(r, ["other"], evidence, f"Synthesis of novel azacyclic allene 4{letter}", "formation of an azacyclic allene")
    r["general_scheme_ids"] = []
    r["participants"] = [
        participant("isoquinoline_precursor", f"{DOC}:mol_precursor_4{letter}", "reactant", evidence),
        participant("methyl_propiolate", f"{DOC}:mol_methyl_propiolate", "reactant", evidence),
        participant("trifluoroethanol", f"{DOC}:mol_trifluoroethanol", "solvent", evidence),
        participant("product", f"{DOC}:mol_4{letter}", "product", evidence, product_role="novel azacyclic allene"),
    ]
    r["solvent_conditions"]["composition"] = fact("trifluoroethanol", [{"molecule_id": f"{DOC}:mol_trifluoroethanol"}], evidence)
    r["steps"] = [step(
        "reaction", "reaction", [p["participant_id"] for p in r["participants"]], evidence,
        operation="Reaction of the corresponding 1-(phenylethynyl)isoquinoline with methyl propiolate in trifluoroethanol",
        temperature={"raw": "7 °C", "value": {"number": 7, "unit": "°C"}},
    )]
    r["outcomes"] = [
        outcome("yield", "yield", [f"{DOC}:mol_4{letter}"], evidence, raw="92%",
                value={"number": 92, "unit": "%"}, product_role="novel azacyclic allene",
                yield_type="reported product yield", basis="yield reported with isolated colorless crystals; exact measurement/isolation method not stated",
                scope=f"formation of 4{letter}"),
        outcome("product_amount", "product_amount", [f"{DOC}:mol_4{letter}"], evidence,
                raw=f"{mass_mg} mg", value={"number": mass_mg, "unit": "mg"},
                product_role="novel azacyclic allene", basis="reported mass of colorless crystals",
                scope=f"formation of 4{letter}"),
    ]
    r["procedure"]["general_procedure_reference"] = fact(
        "Scheme 1 / prose cites the method used for previously reported analogues; no detailed current-paper procedure is printed.",
        "incomplete source procedure", evidence
    )
    unresolved_classification(r, evidence, "The current pages establish reactant classes and product identity but do not support a safe complete atom mapping for the ring-expansion/allene-forming process.")
    mark_review(r, [
        "Exact substrate name/structure, charges, reaction time, atmosphere, workup, purification and yield measurement basis are not reported on the supplied pages.",
        "Scheme 1 is captioned as previous work, while prose calls 4b/4d novel; only source-explicit current compound yields are assigned to these partial records.",
        "Cyclic-allene axial configuration is not assigned.",
    ])
    bundle["reactions"].append(r)


table = {
    "a": {"short7": 59, "short8": 1, "long_time": 3, "long8": 50, "tlc": "1:1", "eluent7": "1:3", "eluent8": "1:2"},
    "b": {"short7": 5, "short8": 2, "long_time": 4, "long8": 58, "tlc": "1:2", "eluent7": "1:5", "eluent8": "1:5"},
    "c": {"short7": 9, "short8": 16, "long_time": 5, "long8": 50, "tlc": "1:1", "eluent7": "1:3", "eluent8": "1:1"},
    "d": {"short7": "trace", "short8": 30, "long_time": 3, "long8": 81, "tlc": "1:2", "eluent7": None, "eluent8": None},
    "e": {"short7": 26, "short8": 2, "long_time": 7.5, "long8": 50, "tlc": "1:2", "eluent7": "1:12", "eluent8": "1:10"},
}


procedure_a_text = (
    "General Procedure A. A solution of the appropriate allene derivative (0.68 mmol) in dry toluene (6.0 mL) "
    "was placed in a microwave reactor and heated at 150 °C for 30 min while monitored by TLC; solvent was "
    "evaporated in vacuo and products were isolated by chromatography on silica gel."
)
procedure_b_text = (
    "General Procedure B. A solution of the appropriate allene derivative (0.16 mmol) in dry toluene (6.0 mL) "
    "was irradiated in a microwave reactor at 150 °C for 3–7.5 h with stirring while monitored by TLC; solvent "
    "was evaporated in vacuo and products 8a–e were crystallized from Et2O."
)


thermolysis_changes = [
    {"atoms": [{"element": "C", "source_locator": "allene-derived carbon"}, {"element": "C", "source_locator": "new ring junction"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
    {"atoms": [{"element": "C", "source_locator": "hydrogen-migration acceptor"}, {"element": "H", "source_locator": "benzyl hydrogen"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain"},
]


for letter, row in table.items():
    # Procedure A: one actual experiment per Table 1 row, with both products as outcomes.
    pages = [2, 4]
    if letter == "a":
        pages += [3, 5]
    evidence = ev(*pages)
    rid = f"{DOC}:procedureA_4{letter}"
    r = st.new_reaction(rid, DOC, evidence)
    r["record_type"] = "specific_experiment"
    set_location(r, f"Table 1 row 4{letter}; General Procedure A", evidence,
                 {"pages": pages, "table": "Table 1", "row": f"4{letter}", "procedure": "A"})
    set_identity(r, ["substrate_scope"], evidence, "Thermolysis of allene 4 to primary cyclopropane 7 and indene 8",
                 "short-time microwave thermolysis/rearrangement")
    r["general_scheme_ids"] = [generic_id]
    add_common_microwave(
        r, evidence, f"4{letter}", 0.68,
        ("30 min", {"number": 30, "unit": "min"}), row["tlc"], procedure=procedure_a_text,
    )
    r["participants"].extend([
        participant("product_7", f"{DOC}:mol_7{letter}", "product", evidence, product_role="primary thermolysis product"),
        participant("product_8", f"{DOC}:mol_8{letter}", "product", evidence, product_role="final thermolysis product formed during short run"),
        participant("silica", f"{DOC}:mol_silica_gel", "purification_agent", evidence, stage="purification"),
        participant("purif_etoac", f"{DOC}:mol_ethyl_acetate", "purification_agent", evidence, stage="purification"),
        participant("purif_hexane", f"{DOC}:mol_hexane", "purification_agent", evidence, stage="purification"),
    ])
    r["steps"][0]["participant_ids"].extend(["product_7", "product_8"])
    r["steps"].append(step("evaporation", "workup", ["substrate", "toluene", "product_7", "product_8"], evidence,
                           operation="Solvent evaporated in vacuo"))
    eluent_desc = (
        f"silica-gel chromatography; EtOAc–hexane {row['eluent7']} for 7{letter}; EtOAc–hexane {row['eluent8']} for 8{letter}"
        if row["eluent7"] and row["eluent8"] else
        f"silica-gel chromatography; specific eluent not assigned in Procedure A for {'7d and 8d' if letter == 'd' else 'this product pair'}"
    )
    r["steps"].append(step("chromatography", "purification", ["product_7", "product_8", "silica", "purif_etoac", "purif_hexane"], evidence,
                           operation=eluent_desc))
    short7 = row["short7"]
    if short7 == "trace":
        o7 = outcome("trace_7", "trace", [f"{DOC}:mol_7{letter}"], evidence,
                     raw="trace", value="trace", product_role="primary thermolysis product",
                     basis="Table 1 qualitative result; no quantitative detection or isolation basis reported",
                     scope="30 min after the start of the experiment")
    else:
        o7 = outcome("yield_7", "yield", [f"{DOC}:mol_7{letter}"], evidence,
                     raw=f"{short7}%", value={"number": short7, "unit": "%"},
                     product_role="primary thermolysis product", yield_type="isolated yield",
                     basis="Yield of isolated product after chromatography (Table 1 footnote a)",
                     scope="30 min after the start of the experiment (Table 1 footnote b)")
    r["outcomes"] = [
        o7,
        outcome("yield_8", "yield", [f"{DOC}:mol_8{letter}"], evidence,
                raw=f"{row['short8']}%", value={"number": row["short8"], "unit": "%"},
                product_role="final thermolysis product formed during short run", yield_type="isolated yield",
                basis="Yield of isolated product after chromatography (Table 1 footnote a)",
                scope="30 min after the start of the experiment (Table 1 footnote b)"),
    ]
    if letter == "a":
        r["outcomes"].extend([
            outcome("amount_7a", "product_amount", [f"{DOC}:mol_7a"], evidence,
                    raw="162 mg", value={"number": 162, "unit": "mg"}, product_role="primary thermolysis product",
                    basis="isolated colorless crystals (Procedure A)", scope="7a from Procedure A"),
            outcome("amount_8a", "product_amount", [f"{DOC}:mol_8a"], evidence,
                    raw="3 mg", value={"number": 3, "unit": "mg"}, product_role="final thermolysis product formed during short run",
                    basis="isolated colorless crystals (Procedure A)", scope="8a from Procedure A"),
        ])
    r["notes"] = [fact(
        "The thermolysis products were formed as single diastereoisomers (7 and 8).",
        "Single diastereoisomer reported for each thermolysis product; absolute configuration and enantiomeric composition unreported.",
        ev(2),
    )]
    if letter == "d":
        r["notes"].append(fact(
            "Procedure A names isolated products 7a–c, 7e, and 8a–e; Table 1 reports 7d only as trace.",
            "7d is retained as a trace outcome, not an isolated-yield outcome.", evidence
        ))
    set_classification(r, evidence, thermolysis_changes)
    unresolved = [
        "Complete structures/atom mapping and absolute stereochemistry are not encoded; all complex molecule identities remain source-label records.",
        "Atmosphere and pressure are not reported on the supplied pages.",
        "Supporting Information is mentioned but was not supplied or associated.",
    ]
    if letter == "d":
        unresolved.append("Procedure A does not assign a specific chromatography eluent to trace 7d or isolated 8d.")
    mark_review(r, unresolved)
    bundle["reactions"].append(r)

    # Procedure B: same table row, distinct scale/time/condition and final product only.
    evidence = ev(2, 4, *( [3, 5] if letter == "a" else [] ))
    rid = f"{DOC}:procedureB_4{letter}"
    r = st.new_reaction(rid, DOC, evidence)
    r["record_type"] = "specific_experiment"
    set_location(r, f"Table 1 row 4{letter}; General Procedure B", evidence,
                 {"pages": [2, 4] + ([3, 5] if letter == "a" else []), "table": "Table 1", "row": f"4{letter}", "procedure": "B"})
    set_identity(r, ["substrate_scope"], evidence, "Complete thermolysis of allene 4 to indene 8",
                 "prolonged microwave thermolysis/rearrangement")
    r["general_scheme_ids"] = [generic_id]
    add_common_microwave(r, evidence, f"4{letter}", 0.16, row["long_time"], row["tlc"], procedure=procedure_b_text)
    r["participants"].extend([
        participant("product_8", f"{DOC}:mol_8{letter}", "product", evidence, product_role="final thermolysis product"),
        participant("ether_crystallization", f"{DOC}:mol_diethyl_ether", "purification_agent", evidence, stage="purification"),
    ])
    r["steps"][0]["participant_ids"].append("product_8")
    r["steps"].extend([
        step("evaporation", "workup", ["substrate", "toluene", "product_8"], evidence,
             operation="Solvent evaporated in vacuo"),
        step("crystallization", "purification", ["product_8", "ether_crystallization"], evidence,
             operation="Product crystallized from Et2O"),
    ])
    r["outcomes"] = [outcome(
        "yield_8", "yield", [f"{DOC}:mol_8{letter}"], evidence,
        raw=f"{row['long8']}%", value={"number": row["long8"], "unit": "%"},
        product_role="final thermolysis product", yield_type="isolated yield",
        basis="Table 1 footnote a says isolated product after chromatography; Procedure B instead states crystallization from Et2O",
        scope=f"after completion of thermolysis; {row['long_time']} h for 4{letter} and 7{letter} to disappear by TLC",
    )]
    if letter == "a":
        r["outcomes"].append(outcome(
            "amount_8a", "product_amount", [f"{DOC}:mol_8a"], evidence,
            raw="32 mg", value={"number": 32, "unit": "mg"}, product_role="final thermolysis product",
            basis="isolated colorless crystals (Procedure B)", scope="8a from Procedure B",
        ))
    r["notes"] = [
        fact("Table 1 footnote a: yield of isolated product after chromatography.",
             "Table footnote assigns chromatography to reported yields.", ev(2)),
        fact("General Procedure B: products 8a–e were crystallized from Et2O.",
             "Procedure B assigns crystallization from Et2O, not chromatography.", ev(4)),
        fact("The thermolysis products were formed as single diastereoisomers (7 and 8).",
             "Single diastereoisomer reported; absolute configuration/enantiomeric composition unreported.", ev(2)),
    ]
    set_classification(r, evidence, thermolysis_changes)
    mark_review(r, [
        "Table 1 footnote a describes isolated product after chromatography, whereas General Procedure B says crystallization from Et2O; both source statements are retained without silently choosing one.",
        "Complete structures/atom mapping and absolute stereochemistry are not encoded.",
        "Atmosphere and pressure are unreported; Supporting Information is mentioned but not supplied or associated.",
    ])
    bundle["reactions"].append(r)

# Procedure C: isolated 7a to 8a.
c_ev = ev(3, 5)
c = st.new_reaction(f"{DOC}:procedureC_7a_to_8a", DOC, c_ev)
c["record_type"] = "specific_experiment"
set_location(c, "Page 3 prose and General Procedure C (page 5)", c_ev,
             {"pages": [3, 5], "procedure": "C"})
set_identity(c, ["mechanistic_study"], c_ev, "Thermal rearrangement of isolated 7a to 8a",
             "vinylcyclopropane-to-cyclopentene rearrangement")
add_common_microwave(c, c_ev, "7a", 0.24, 2, "1:1", procedure=(
    "General Procedure C. A solution of 7a (0.24 mmol) in absolute toluene (6.0 mL) was irradiated in a microwave "
    "reactor at 150 °C for 2 h with stirring while monitored by TLC; solvent was evaporated in vacuo and 8a was "
    "crystallized from Et2O."
), dry_word="absolute")
c["participants"].extend([
    participant("product_8a", f"{DOC}:mol_8a", "product", c_ev, product_role="rearrangement product"),
    participant("ether_crystallization", f"{DOC}:mol_diethyl_ether", "purification_agent", c_ev, stage="purification"),
])
c["steps"][0]["participant_ids"].append("product_8a")
c["steps"].extend([
    step("evaporation", "workup", ["substrate", "toluene", "product_8a"], c_ev,
         operation="Solvent evaporated in vacuo"),
    step("crystallization", "purification", ["product_8a", "ether_crystallization"], c_ev,
         operation="8a crystallized from Et2O"),
])
c["outcomes"] = [
    outcome("yield", "yield", [f"{DOC}:mol_8a"], c_ev, raw="50%",
            value={"number": 50, "unit": "%"}, product_role="rearrangement product",
            yield_type="isolated yield", basis="product crystallized from Et2O",
            scope="isolated 7a to 8a under Procedure C"),
    outcome("product_amount", "product_amount", [f"{DOC}:mol_8a"], c_ev, raw="49 mg",
            value={"number": 49, "unit": "mg"}, product_role="rearrangement product",
            basis="product crystallized from Et2O", scope="isolated 7a to 8a under Procedure C"),
]
c["notes"] = [fact(
    "The authors use this experiment to establish cyclopropane 7a as an intermediate leading to 8a.",
    "Mechanistic-study interpretation: isolated-intermediate conversion supports 7a→8a.", c_ev
)]
set_classification(c, c_ev, [
    {"atoms": [{"element": "C", "source_locator": "cyclopropane carbon"}, {"element": "C", "source_locator": "new cyclopentene bond partner"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
    {"atoms": [{"element": "C", "source_locator": "cleaved cyclopropane carbon"}, {"element": "C", "source_locator": "cleaved ring partner"}], "before_order": 1, "after_order": 0, "change_type": "broken"},
    {"atoms": [{"element": "C", "source_locator": "new cyclopentene alkene carbon"}, {"element": "C", "source_locator": "new cyclopentene alkene carbon"}], "before_order": 1, "after_order": 2, "change_type": "order_increased"},
    {"atoms": [{"element": "C", "source_locator": "former vinyl bond carbon"}, {"element": "C", "source_locator": "former vinyl bond carbon"}], "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
])
mark_review(c, [
    "Complete atom mapping and absolute stereochemistry are not encoded; classification is a partial net rearrangement description.",
    "Atmosphere and pressure are not reported; Supporting Information is mentioned but was not supplied or associated.",
])
bundle["reactions"].append(c)

# Proposed radical mechanism is not a performed experiment.
m_ev = ev(3)
m = st.new_reaction(f"{DOC}:scheme3_proposed_biradical_mechanism", DOC, m_ev)
m["record_type"] = "proposed_mechanism"
set_location(m, "Scheme 3", m_ev, {"page": 3, "scheme": "Scheme 3"})
set_identity(m, ["mechanistic_study"], m_ev, "Proposed mechanism via biradical intermediates",
             "author-proposed radical rearrangement pathway")
m["participants"] = [
    participant("generic_4", f"{DOC}:mol_generic_4", "reactant", m_ev),
    participant("generic_7", f"{DOC}:mol_generic_7", "product", m_ev, product_role="proposed primary product/intermediate"),
    participant("generic_8", f"{DOC}:mol_generic_8", "product", m_ev, product_role="proposed final product"),
]
m["steps"] = [step(
    "proposed_pathway", "reaction", ["generic_4", "generic_7", "generic_8"], m_ev,
    operation="Author-proposed sequence 4→A→B→C⇌7 and C⇌D→8 via biradicals and C–H insertion",
)]
m["notes"] = [fact(
    "Because rearrangement occurred upon heating to 150 °C, a radical mechanism can be presumed; the authors assume biradical intermediates A–D.",
    "Author proposal, not direct experimental proof of intermediates A–D.", m_ev
)]
m["classification"]["labels"] = inferred([], [], m_ev, "A proposed elementary pathway is not classified as one performed net experiment.")
m["classification"]["coverage"] = inferred("not_applicable", "not_applicable", m_ev, "No independent performed reaction record is represented by the proposed mechanism drawing.")
m["classification"]["limitations"] = ["Biradical intermediates A–D are proposed, not experimentally isolated or directly observed on the supplied pages."]
m["classification"]["model_votes"] = []
mark_review(m, ["Proposed biradical intermediates are author interpretation, not experimentally verified species."])
bundle["reactions"].append(m)

# Non-synthetic observations and metadata. These do not create reaction records.
bundle["source_context"] = [
    {
        "context_id": f"{DOC}:bibliographic_metadata",
        "kind": "bibliographic_metadata",
        "evidence_ids": ev(1),
        "title": "Unusual Transformations of Cyclic Allenes with an Enamine Moiety into Complex Frameworks",
        "journal": "Synlett",
        "year": 2020,
        "volume": 31,
        "doi": "10.1055/s-0039-1691586",
    },
    {
        "context_id": f"{DOC}:thermolysis_stereochemistry_nmr",
        "kind": "stereochemistry_and_characterization",
        "evidence_ids": ev(2, 3),
        "finding": "Products 7 and 8 are reported as single diastereoisomers; diagnostic 1H NMR ranges are discussed. No absolute R/S, allene-axis, ee/er, or E/Z assignment is reported.",
    },
    {
        "context_id": f"{DOC}:xray_7a_8a",
        "kind": "xray_characterization",
        "evidence_ids": ev(3, 5),
        "finding": "Single-crystal X-ray analyses confirm the structures of 7a and 8a; CCDC 1830405 and 1579823 are cited. These observations are characterization, not separate syntheses.",
    },
    {
        "context_id": f"{DOC}:spectral_characterization_4b_4d_7a_8a",
        "kind": "spectral_characterization",
        "evidence_ids": ev(4, 5),
        "finding": "IR, NMR, LC-MS and elemental-analysis data for 4b, 4d, 7a and 8a are retained as characterization context and are not separate reaction records.",
    },
    {
        "context_id": f"{DOC}:si_availability_unassociated",
        "kind": "supplementary_material_statement",
        "evidence_ids": ev(4),
        "finding": "The article states that Supporting Information and Primary Data are available online, but no SI file was supplied or explicitly associated in this workspace; experimental gaps remain unfilled.",
    },
    {
        "context_id": f"{DOC}:page5_references",
        "kind": "bibliography",
        "evidence_ids": ev(5),
        "finding": "The remainder of page 5 contains references after Procedure C; citation entries are not reconstructed as experiments.",
    },
]

# Raw page inventory and exact coverage contract.
page_records = {
    1: [
        f"{DOC}:scheme1_literature_allene1", f"{DOC}:scheme1_literature_allene2",
        f"{DOC}:scheme1_literature_lactam3", f"{DOC}:scheme1_previous_allene4",
        generic_id, f"{DOC}:synthesis_4b", f"{DOC}:synthesis_4d",
    ],
    2: [
        f"{DOC}:scheme1_previous_allene4", f"{DOC}:scheme2_literature_path_a",
        f"{DOC}:scheme2_literature_path_b", f"{DOC}:scheme2_literature_allene2_to_6",
        generic_id, f"{DOC}:synthesis_4b", f"{DOC}:synthesis_4d",
    ] + [f"{DOC}:procedureA_4{x}" for x in subs] + [f"{DOC}:procedureB_4{x}" for x in subs],
    3: [generic_id, f"{DOC}:procedureA_4a", f"{DOC}:procedureB_4a",
        f"{DOC}:procedureC_7a_to_8a", f"{DOC}:scheme3_proposed_biradical_mechanism"],
    4: [f"{DOC}:synthesis_4b", f"{DOC}:synthesis_4d"] +
       [f"{DOC}:procedureA_4{x}" for x in subs] + [f"{DOC}:procedureB_4{x}" for x in subs],
    5: [f"{DOC}:procedureA_4a", f"{DOC}:procedureB_4a", f"{DOC}:procedureC_7a_to_8a"],
}

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 1,
            "state": "attempted",
            "regions": [
                {"region": "title, graphical abstract and abstract", "finding": "Bibliographic metadata and a duplicate summary of the current thermolysis sequence; no extra experiment beyond Table 1/Procedures A–B.", "evidence_ids": ev(1)},
                {"region": "Scheme 1", "finding": "Four contextual prior-work arrows: syntheses of 1, 2, 3 and the authors' earlier generic route to 4; classified as literature examples, not current experiments.", "evidence_ids": ev(1)},
                {"region": "lower prose", "finding": "Introduction/prior-work discussion; cross-page continuation only, with no independent current experiment on this page.", "evidence_ids": ev(1)},
            ],
            "record_ids": page_records[1],
            "unresolved": ["Generic Scheme 1 structures contain R/n/X variables and are retained as templates.", "The title graphic is a summary and does not specify substrate-specific long times."],
        },
        {
            "page": 2,
            "state": "attempted",
            "regions": [
                {"region": "Scheme 2", "finding": "Three cited prior-work transformations: cyclohexadiene paths a/b and allene 2 isomerization to diene 6; not current experiments.", "evidence_ids": ev(2)},
                {"region": "Table 1", "finding": "Five short Procedure-A experiments with two product outcomes each, plus five distinct completion/Procedure-B experiments; all five rows and footnotes inventoried.", "evidence_ids": ev(2)},
                {"region": "prose surrounding Table 1", "finding": "Identifies 4b/4d as novel, describes TLC evolution 4/7→8, and reports products 7/8 as single diastereoisomers; no absolute configurations assigned.", "evidence_ids": ev(2)},
            ],
            "record_ids": page_records[2],
            "unresolved": ["Dense structures 4a–e, 7a–e and 8a–e are source-labelled but not safely encoded as complete isomeric SMILES.", "Exact enantiomeric composition and absolute stereochemistry are unreported."],
        },
        {
            "page": 3,
            "state": "attempted",
            "regions": [
                {"region": "Figures 1 and 2", "finding": "Single-crystal X-ray renderings of 7a and 8a; characterization evidence, not separate syntheses.", "evidence_ids": ev(3)},
                {"region": "Scheme 3", "finding": "Author-proposed biradical pathway A–D, recorded as proposed mechanism rather than experimental proof.", "evidence_ids": ev(3)},
                {"region": "right-column prose", "finding": "Reports the performed isolated-intermediate experiment 7a→8a in 50% yield and explains its mechanistic interpretation.", "evidence_ids": ev(3)},
            ],
            "record_ids": page_records[3],
            "unresolved": ["Biradicals A–D are proposed and not directly observed.", "X-ray pictures confirm connectivity/relative arrangement but no absolute R/S naming is supplied on the page."],
        },
        {
            "page": 4,
            "state": "attempted",
            "regions": [
                {"region": "Notes 13–14", "finding": "Product yield/amount and characterization for novel allenes 4b and 4d; spectral peak lists retained as source context, not reactions.", "evidence_ids": ev(4)},
                {"region": "General Procedure A", "finding": "0.68 mmol substrate, dry toluene 6.0 mL, microwave 150 °C/30 min, TLC systems, evaporation, silica chromatography and product-specific eluents transcribed.", "evidence_ids": ev(4)},
                {"region": "General Procedure B", "finding": "0.16 mmol substrate, dry toluene 6.0 mL, microwave 150 °C/3–7.5 h, TLC, evaporation and Et2O crystallization transcribed.", "evidence_ids": ev(4)},
                {"region": "Note 17 and 8a characterization start", "finding": "7a mass/yield and 8a Procedure-A/Procedure-B amounts/yields tied to their respective experiments; NMR/MS/analysis remain characterization context.", "evidence_ids": ev(4)},
            ],
            "record_ids": page_records[4],
            "unresolved": ["No detailed preparation procedure for 4b/4d is printed.", "Procedure A does not assign an eluent to 8d; 7d is trace.", "Table footnote says chromatography for yields while Procedure B says crystallization from Et2O."],
        },
        {
            "page": 5,
            "state": "attempted",
            "regions": [
                {"region": "upper-left continuation", "finding": "Continuation of 8a characterization and CCDC identifiers; no new synthesis record.", "evidence_ids": ev(5)},
                {"region": "General Procedure C", "finding": "Distinct mechanistic-study experiment: 7a 0.24 mmol in absolute toluene 6.0 mL, microwave 150 °C/2 h, TLC, evaporation, Et2O crystallization, 8a 49 mg (50%).", "evidence_ids": ev(5)},
                {"region": "remaining references", "finding": "Bibliography only; cited papers are not reconstructed as experiments.", "evidence_ids": ev(5)},
            ],
            "record_ids": page_records[5],
            "unresolved": ["Procedure C atmosphere and pressure are unreported.", "No supplied SI is associated despite the article's online-SI statement."],
        },
    ],
}

raw_inventory = {
    "batch_id": "allene-stereo-v31-200-20261017",
    "document_id": DOC,
    "title": "Unusual Transformations of Cyclic Allenes with an Enamine Moiety into Complex Frameworks",
    "visual_pages_reviewed": [1, 2, 3, 4, 5],
    "reaction_record_count": len(bundle["reactions"]),
    "molecule_count": len(bundle["molecules"]),
    "record_ids_by_page": page_records,
    "decisions": {
        "specific_current_experiments": 13,
        "literature_examples": 7,
        "generic_schemes": 1,
        "proposed_mechanisms": 1,
        "short_table_rows": 5,
        "completion_table_rows": 5,
        "novel_allene_preparations": 2,
        "isolated_intermediate_experiments": 1,
        "si_associated": False,
        "spectral_peak_tables_as_reactions": False,
    },
}

st.write(RUN / "submitted.json", bundle)
st.write(RUN / "raw_model/pages_0001_0005_extraction.json", raw_inventory)
st.write(RUN / "raw_model/stereo_checks.json", stereo_checks)
st.write(RUN / "coverage/pages_0001_0005.json", coverage)
print(json.dumps({
    "molecules": len(bundle["molecules"]),
    "reactions": len(bundle["reactions"]),
    "source_context": len(bundle["source_context"]),
    "stereo_checks": len(stereo_checks),
}, ensure_ascii=False, indent=2))
