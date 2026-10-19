#!/usr/bin/env python3
"""Source transcription for DOI 10.1021/acs.orglett.1c01234, PDF pages 1-5.

This file is a machine-authored extraction aid.  It deliberately leaves dense
source-drawn structures unresolved and all imported records as needs_review.
"""
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_175faee9a12c1b16"
RUN = ROOT / "runs/current"


def ev(page, image=True):
    return [f"{DOC}:page{page}:{'image' if image else 'text'}"]


def ev_both(page):
    return [f"{DOC}:page{page}:image", f"{DOC}:page{page}:text"]


def f(raw=None, norm=None, evidence=(), method="original_extraction", status=None,
      derivation=None):
    return st.fact(raw, norm, evidence, method=method, status=status,
                   derivation=derivation)


def inferred(raw, norm, evidence, derivation):
    return f(raw, norm, evidence, method="inference", derivation=derivation)


def derived(raw, norm, evidence, derivation):
    return f(raw, norm, evidence, method="structure_derivation", derivation=derivation)


def not_applicable(evidence):
    return f(evidence=evidence, method="not_applicable", status="not_applicable")


def mid(label):
    return f"{DOC}:mol_{label}"


bundle = json.loads((RUN / "submitted.json").read_text())
molecules = []
stereo_diagnostics = {}


def add_unresolved(label, name, evidence, chirality="unresolved", geometry="unresolved",
                   assignment=None, axial=None, double_bonds=None, note=None):
    m = st.new_molecule(mid(label), DOC, evidence)
    m["identity_status"] = "source_identified"
    m["representation_type"] = "unresolved"
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    m["names"]["original"] = f(name, name, evidence)
    m["names"]["english"] = f(name, name, evidence)
    m["cas_numbers"] = f(evidence=evidence)
    m["smiles"] = f(evidence=evidence)
    m["generic_template"] = None
    if assignment is not None:
        m["stereochemistry"]["assignment"] = f(assignment, assignment, evidence)
        m["stereochemistry"]["atom_centers"] = f(
            "source-drawn stereochemical bonds; exact atom indices unavailable",
            [{"assignment": assignment, "atom_indices": None}], evidence)
    else:
        m["stereochemistry"]["assignment"] = f(evidence=evidence)
        m["stereochemistry"]["atom_centers"] = (
            derived("source drawing has no tetrahedral stereocenter", [], evidence,
                    "Visual comparison of the complete source drawing; no stereogenic tetrahedral center is present.")
            if chirality == "achiral" else f(evidence=evidence))
    if axial is not None:
        m["stereochemistry"]["axial"] = derived(
            axial, axial, evidence, "Allene/biaryl axial status derived from the source drawing and explicit ligand name.")
    elif chirality == "achiral" and "allen" in name.lower():
        axis = {"axis": "terminal allene C=C=C", "status": "not stereogenic",
                "reason": "the terminal CH2 end bears two identical hydrogen substituents"}
        m["stereochemistry"]["axial"] = derived(
            "terminal allene in the source drawing", axis, evidence,
            "A terminal CH2=C=C unit cannot carry configurational allene axial chirality at that end.")
    else:
        m["stereochemistry"]["axial"] = f(evidence=evidence)
    m["stereochemistry"]["double_bonds"] = (
        derived("source-drawn terminal allene/alkene; E/Z is not applicable", double_bonds or [], evidence,
                "Terminal allene/alkene geometry does not define an E/Z pair.")
        if geometry == "not_applicable" else f(evidence=evidence))
    m["stereochemistry"]["mixture"] = f(evidence=evidence)
    status_method = "structure_derivation" if chirality == "achiral" else "original_extraction"
    status_deriv = "Source structure inspected for tetrahedral and axial stereogenicity." if status_method == "structure_derivation" else None
    m["stereochemistry"]["chirality_status"] = f(
        assignment or ("source drawing" if chirality != "unresolved" else "not resolved"),
        chirality, evidence, method=status_method, derivation=status_deriv)
    m["stereochemistry"]["geometry_status"] = (
        derived("terminal unsaturation", geometry, evidence,
                "Only terminal alkene/allene unsaturation is shown, so E/Z is not applicable.")
        if geometry == "not_applicable" else f("not resolved", geometry, evidence))
    m["stereochemistry"]["representation_status"] = inferred(
        "structure retained as source drawing", "not_encoded", evidence,
        "The dense source structure and its stereochemical bonds were not converted to a complete isomeric SMILES.")
    m["stereochemistry"]["isomeric_smiles"] = f(evidence=evidence)
    m["stereochemistry"]["identity_key"] = f(evidence=evidence)
    m["identity_note"] = f(
        note or "Identity is tied to the compound label and source drawing; complete stereochemical SMILES is unreported.",
        note or "Identity is tied to the compound label and source drawing; complete stereochemical SMILES is unreported.",
        evidence)
    m["external_properties"] = []
    molecules.append(m)
    return m


def add_smiles(label, name, smiles, evidence):
    diag = st.stereo_inspect(smiles)
    stereo_diagnostics[mid(label)] = {"input_smiles": smiles, **diag}
    m = st.new_molecule(mid(label), DOC, evidence)
    m["identity_status"] = "source_identified"
    m["representation_type"] = "smiles"
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    m["names"]["original"] = f(name, name, evidence)
    m["names"]["english"] = f(name, name, evidence)
    m["cas_numbers"] = f(evidence=evidence)
    m["smiles"] = derived(
        {"source_label": label, "source_name": name}, smiles, evidence,
        "Candidate SMILES transcribed from an unambiguous small named reagent/solvent and checked by the bundled standard stereo graph diagnostic; source identity is still needs_review.")
    m["generic_template"] = None
    atom_potential = [x for x in diag["potential_stereo"] if str(x["type"]).startswith("Atom_")]
    bond_potential = [x for x in diag["potential_stereo"] if str(x["type"]).startswith("Bond_")]
    chirality = "unresolved" if atom_potential else "achiral"
    geometry = "unresolved" if bond_potential else "not_applicable"
    m["stereochemistry"]["atom_centers"] = derived(
        "standard stereo diagnostic", diag["atom_centers"], evidence,
        "RDKit graph diagnostic; absence of a toolkit center is combined with visual review only for these simple named molecules.")
    m["stereochemistry"]["double_bonds"] = derived(
        "standard stereo diagnostic", diag["double_bonds"], evidence,
        "RDKit graph diagnostic; this does not infer unreported E/Z configuration.")
    m["stereochemistry"]["axial"] = f(evidence=evidence)
    m["stereochemistry"]["mixture"] = f(evidence=evidence)
    m["stereochemistry"]["assignment"] = f(evidence=evidence)
    m["stereochemistry"]["chirality_status"] = derived(
        "simple named structure and stereo diagnostic", chirality, evidence,
        "No potential tetrahedral stereocenter is present in the candidate graph." if chirality == "achiral" else "Potential stereochemistry remains unresolved.")
    m["stereochemistry"]["geometry_status"] = derived(
        "standard stereo diagnostic", geometry, evidence,
        "No potential E/Z bond is present in the candidate graph." if geometry == "not_applicable" else "Potential bond geometry remains unresolved.")
    m["stereochemistry"]["representation_status"] = derived(
        "candidate contains no applicable stereochemical encoding", "not_applicable" if chirality == "achiral" and geometry == "not_applicable" else "partial",
        evidence, "Representation status follows the candidate graph diagnostic, not a claim of source verification.")
    m["stereochemistry"]["isomeric_smiles"] = derived(
        smiles, diag["canonical_isomeric_smiles"], evidence,
        "Canonical isomeric SMILES generated by the bundled standard stereo diagnostic.")
    m["stereochemistry"]["identity_key"] = derived(
        "standard stereo graph key", diag["stereo_graph_key"], evidence,
        "SHA-256 graph key from canonical isomeric SMILES; not a sample identity.")
    m["identity_note"] = f(
        "Candidate structure for a simple named participant; parsing/graph diagnostics do not establish source correctness.",
        "Candidate structure for a simple named participant; parsing/graph diagnostics do not establish source correctness.", evidence)
    m["external_properties"] = []
    molecules.append(m)
    return m


def add_template(label, name, evidence, fixed, variables):
    m = st.new_molecule(mid(label), DOC, evidence)
    m["identity_status"] = "source_identified"
    m["representation_type"] = "generic_template"
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    m["names"]["original"] = f(name, name, evidence)
    m["names"]["english"] = f(name, name, evidence)
    m["cas_numbers"] = f(evidence=evidence)
    m["smiles"] = f(evidence=evidence)
    template_vars = []
    for item in variables:
        template_vars.append({
            "label": item,
            "attachment_points": f(f"source attachment of {item}", {"label": item, "attachment": "as drawn"}, evidence),
            "definition": f(item, item, evidence),
            "allowed_values": f(evidence=evidence),
            "constraints": f(evidence=evidence),
        })
    m["generic_template"] = {
        "fixed_structure": f("source generic drawing", fixed, evidence),
        "variables": template_vars,
        "instance_mapping": [],
        "reason": f("generic R-group structure in Scheme 1", "Source is a generic family; no concrete SMILES was invented.", evidence),
    }
    for key in ("atom_centers", "double_bonds", "axial", "mixture", "assignment"):
        m["stereochemistry"][key] = f(evidence=evidence)
    m["stereochemistry"]["chirality_status"] = f("generic family", "unresolved", evidence)
    m["stereochemistry"]["geometry_status"] = f("generic family", "unresolved", evidence)
    m["stereochemistry"]["representation_status"] = inferred(
        "generic source drawing", "not_encoded", evidence,
        "Variable source family cannot be represented as one complete isomeric SMILES.")
    m["stereochemistry"]["isomeric_smiles"] = f(evidence=evidence)
    m["stereochemistry"]["identity_key"] = f(evidence=evidence)
    m["identity_note"] = f("generic Scheme 1 family", "generic Scheme 1 family", evidence)
    m["external_properties"] = []
    molecules.append(m)
    return m


# Generic previous/current-work families on page 2.
add_template("litA_sub", "previous-work malonate-tethered 3-allenylindole family A", ev(2),
             {"scaffold": "malonate-tethered 3-allenylindole", "product_class": "tetrahydrocarbazole"}, ["R"])
add_template("litA_prod", "previous-work tetrahydrocarbazole family A", ev(2),
             {"scaffold": "functionalized tetrahydrocarbazole"}, ["R"])
add_template("litB_sub", "previous-work malonate-tethered 3-allenylindole family B", ev(2),
             {"scaffold": "malonate-tethered 3-allenylindole"}, ["R"])
add_template("litB_prod", "previous-work five-membered spiroindoline family B", ev(2),
             {"scaffold": "five-membered spiroindoline"}, ["R"])
add_template("fam1", "N-allenyltryptamine family 1", ev(2),
             {"scaffold": "3-tethered N-allenyltryptamine", "terminal_allene": True}, ["R", "NPG"])
add_template("fam2", "six-membered vinylic spiroindolenine family 2", ev(2),
             {"scaffold": "six-membered aza-spiroindolenine", "terminal_alkene": True}, ["R", "NPG"])
add_template("fam3", "malonate-tethered 3-allenylindole family 3", ev(2),
             {"scaffold": "malonate-tethered 3-allenylindole", "terminal_allene": True}, ["R"])
add_template("fam4", "malonate-tethered vinylic spiroindolenine family 4", ev(2),
             {"scaffold": "six-membered carbocyclic spiroindolenine", "terminal_alkene": True}, ["R"])

# Source-labelled substrates and products.
one_names = {
    "1a": "tosyl-protected N-allenyltryptamine 1a",
    "1b": "p-nosyl-protected N-allenyltryptamine 1b",
    "1c": "o-nosyl-protected N-allenyltryptamine 1c",
    "1d": "benzenesulfonyl-protected N-allenyltryptamine 1d",
    "1e": "MBS-protected N-allenyltryptamine 1e",
    "1f": "Cbz-protected N-allenyltryptamine 1f",
    "1g": "5-methoxy p-nosyl N-allenyltryptamine 1g",
    "1h": "5-OTBS p-nosyl N-allenyltryptamine 1h",
    "1i": "tryptophan-derived methyl ester N-allenyl substrate 1i",
}
for label, name in one_names.items():
    add_unresolved(label, name, ev(2),
                   chirality="partially_specified" if label == "1i" else "achiral",
                   geometry="not_applicable",
                   assignment="tryptophan-derived stereocenter present; exact descriptor not printed on this page" if label == "1i" else None,
                   note="Complete source drawing inspected; no complete SMILES was asserted because dense connectivity/protecting-group expansion and stereochemistry require SI-level confirmation.")

two_names = {
    "2a": "(S,S)-spiroindolenine 2a (N-Ts)",
    "2b": "spiroindolenine 2b (N-p-Ns)",
    "2c": "spiroindolenine 2c (N-o-Ns)",
    "2d": "spiroindolenine 2d (N-SO2Ph)",
    "2e": "spiroindolenine 2e (N-MBS)",
    "2f": "spiroindolenine 2f (N-Cbz)",
    "2g": "5-methoxy spiroindolenine 2g (N-p-Ns)",
    "2h": "5-OTBS spiroindolenine 2h (N-p-Ns)",
    "2i": "tryptophan-derived methyl ester spiroindolenine 2i",
}
for label, name in two_names.items():
    if label == "2a":
        assignment = "absolute configuration (S,S), determined by X-ray crystal structure"
        status = "specified"
    elif label == "2i":
        assignment = "major product stereochemistry shown; corresponding centers assigned by analogy to 2a; four-diastereomer ratio reported"
        status = "partially_specified"
    else:
        assignment = "absolute configuration assigned by analogy to X-ray structure of 2a"
        status = "specified"
    add_unresolved(label, name, ev_both(2), chirality=status, geometry="not_applicable",
                   assignment=assignment,
                   note="Six-membered source-drawn spiroindolenine; stereochemistry is retained textually and is not encoded as SMILES.")

three_names = {
    "3a": "unsubstituted malonate-tethered 3-allenylindole 3a",
    "3b": "5-bromo malonate-tethered 3-allenylindole 3b",
    "3c": "5-iodo malonate-tethered 3-allenylindole 3c",
    "3d": "source-defined methoxy malonate-tethered 3-allenylindole 3d",
    "3e": "source-defined methoxy positional isomer 3e",
    "3f": "source-defined methyl malonate-tethered 3-allenylindole 3f",
}
for label, name in three_names.items():
    add_unresolved(label, name, ev(3), chirality="achiral", geometry="not_applicable",
                   note="Exact substitution position for 3d-3f is retained by label/source drawing and not guessed from the dense figure.")

four_names = {
    "4a": "unsubstituted malonate-tethered vinylic spiroindolenine 4a",
    "4b": "5-bromo malonate-tethered vinylic spiroindolenine 4b",
    "4c": "5-iodo malonate-tethered vinylic spiroindolenine 4c",
    "4d": "source-defined methoxy spiroindolenine 4d",
    "4e": "source-defined methoxy positional isomer spiroindolenine 4e",
    "4f": "source-defined methyl spiroindolenine 4f",
}
for label, name in four_names.items():
    add_unresolved(label, name, ev_both(3), chirality="specified", geometry="not_applicable",
                   assignment="absolute configuration assigned by analogy to 2a; relative stereochemical bonds shown",
                   note="The (rac) labels on the 4b/4c ORTEP images describe centrosymmetric crystallographic samples and are not used to overwrite the enantioenriched scope-product identity.")

for label, name, assignment in [
    ("5", "selectively imine-reduced spiroindoline 5", "major diastereomer drawn; configuration retained from 2b, exact descriptors not printed"),
    ("6a", "fully hydrogenated p-ABS spiroindoline 6a", "major diastereomer drawn; exact descriptors not printed"),
    ("6b", "fully hydrogenated MBS spiroindoline 6b", "major diastereomer drawn; exact descriptors not printed"),
    ("7", "spirooxindole 7", "major diastereomer drawn; configuration retained from 2b, exact descriptors not printed"),
    ("8", "deprotected spiroindolenine 8", "major diastereomer drawn; configuration retained from 2b, exact descriptors not printed"),
]:
    add_unresolved(label, name, ev_both(4), chirality="partially_specified", geometry="not_applicable",
                   assignment=assignment,
                   note="No loss of optical purity was stated collectively for Scheme 4 transformations; that statement is not converted into an unreported absolute descriptor.")

# Named catalysts, ligands, reagents, additives, and solvents.  Only simple,
# unambiguous species receive candidate SMILES.
add_unresolved("Rh_cod_Cl_dimer", "[Rh(cod)Cl]2", ev_both(2), note="Source catalyst formula; organometallic connectivity not encoded.")
add_unresolved("R_DTBM_Segphos", "(R)-DTBM-Segphos", ev_both(2), chirality="specified",
               axial={"axis": "biaryl", "configuration": "R"},
               note="Commercial chiral ligand drawing/name retained; full organophosphorus SMILES not asserted.")
add_unresolved("R_BINAP", "(R)-BINAP", ev_both(2), chirality="specified",
               axial={"axis": "binaphthyl", "configuration": "R"},
               note="Axial configuration is explicit in the ligand name and drawing; full SMILES not asserted.")
add_smiles("PPTS", "pyridinium p-toluenesulfonate", "Cc1ccc(S(=O)(=O)[O-])cc1.c1cc[nH+]cc1", ev_both(2))
add_smiles("DCM", "dichloromethane", "ClCCl", ev_both(2))
add_unresolved("HEH", "Hantzsch ester (HEH)", ev(2), note="Source drawing retained; no candidate SMILES needed for current-work experiments.")

small = [
    ("benzaldehyde", "benzaldehyde", "O=Cc1ccccc1", ev_both(3)),
    ("benzoic_acid", "benzoic acid", "O=C(O)c1ccccc1", ev_both(3)),
    ("ethanol", "ethanol", "CCO", ev_both(3)),
    ("phenol", "phenol", "Oc1ccccc1", ev_both(3)),
    ("terminal_allene_additive", "penta-3,4-dien-1-ylbenzene", "C=C=CCCc1ccccc1", ev_both(3)),
    ("internal_alkyne_additive", "prop-1-yn-1-ylbenzene", "CC#Cc1ccccc1", ev_both(3)),
    ("TMS_diene", "2-TMS-buta-1,3-diene", "C=C([Si](C)(C)C)C=C", ev_both(3)),
    ("water", "water", "O", ev_both(3)),
    ("MeOH", "methanol", "CO", ev_both(4)),
    ("H2", "hydrogen", "[H][H]", ev_both(4)),
    ("2methyl2butene", "2-methyl-2-butene", "CC=C(C)C", ev_both(4)),
    ("THF", "tetrahydrofuran", "C1CCOC1", ev_both(4)),
    ("p_thiocresole", "p-thiocresole", "Cc1ccc(S)cc1", ev_both(4)),
    ("K2CO3", "potassium carbonate", "[K+].[K+].[O-]C([O-])=O", ev_both(4)),
    ("DMF", "N,N-dimethylformamide", "CN(C)C=O", ev_both(4)),
    ("TsOH", "p-toluenesulfonic acid", "Cc1ccc(S(=O)(=O)O)cc1", ev_both(4)),
    ("HCl", "hydrogen chloride", "Cl", ev_both(4)),
]
for args in small:
    add_smiles(*args)
add_unresolved("4_bromopyrazole", "4-bromopyrazole", ev_both(3), note="Positional heteroaromatic structure is retained by source name; no candidate SMILES asserted.")
add_unresolved("NaBH3CN", "sodium cyanoborohydride", ev_both(4), note="Ionic reagent formula retained; no candidate SMILES asserted.")
add_unresolved("Pd_C", "Pd/C (5 wt %, 10 mol % Pd)", ev_both(4), note="Heterogeneous catalyst/material; not a single molecular SMILES.")
add_unresolved("NaClO2", "sodium chlorite", ev_both(4), note="Ionic reagent formula retained; no candidate SMILES asserted.")
add_unresolved("NaH2PO4_aq", "aqueous NaH2PO4 (1 M)", ev_both(4), note="Aqueous reagent solution retained as named; concentration is recorded at reaction level.")

bundle["molecules"] = molecules
reactions = []


def new_reaction(suffix, page, record_type, categories, location, name, transformation, evidence=None):
    evidence = evidence or ev_both(page)
    r = st.new_reaction(f"{DOC}:{suffix}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, evidence)
    r["location"] = f(location, {"page": page, "region": location}, evidence)
    r["reaction_name"] = f(name, name, evidence) if name else f(evidence=evidence)
    r["transformation_type"] = inferred(
        transformation, transformation, evidence,
        "Transformation class assigned from the source arrow/product class; it is not a verbatim reaction name.") if transformation else f(evidence=evidence)
    r["review"] = {
        "status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
        "unresolved": [
            "Supporting Information was announced but not supplied or explicitly associated.",
            "Workup and purification are unreported in the supplied main-text pages.",
            "Dense source-drawn primary structures are not represented by complete atom-mapped isomeric SMILES.",
        ],
    }
    r["completeness"] = "partial"
    reactions.append(r)
    return r


def participant(r, pid, mol_label, role, evidence, stage="reaction", product_role=None):
    p = st.new_participant(pid, mid(mol_label), role, evidence, stage=stage)
    if product_role:
        p["product_role"] = f(product_role, product_role, evidence)
    r["participants"].append(p)
    return p


def set_charge(p, key, raw, norm, evidence, reference=None):
    p["charges"][key] = f(raw, norm, evidence)
    if reference is not None:
        p["equivalence_reference"] = f(reference, mid(reference), evidence)


def solvent(r, label, evidence, raw_name=None):
    p = participant(r, f"p_{label}", label, "solvent", evidence)
    r["solvent_conditions"]["composition"] = f(raw_name or label, [mid(label)], evidence)
    return p


def reaction_step(r, step_id, evidence, operation, temperature=None, time=None,
                  atmosphere=None, pressure=None):
    s = st.new_step(step_id, evidence)
    s["participant_ids"] = [p["participant_id"] for p in r["participants"]]
    s["operation"] = f(operation, operation, evidence)
    if temperature is not None:
        s["temperature"] = f(temperature[0], temperature[1], evidence)
    if time is not None:
        s["time"] = f(time[0], time[1], evidence)
    if atmosphere is not None:
        s["atmosphere"] = f(atmosphere, atmosphere, evidence)
    if pressure is not None:
        s["pressure"] = f(pressure[0], pressure[1], evidence)
    r["steps"].append(s)
    return s


def outcome(r, oid, kind, target, raw, norm, evidence, basis=None, yield_type=None, scope=None):
    targets = [mid(target)] if isinstance(target, str) else [mid(x) for x in target]
    o = st.new_outcome(oid, kind, targets, evidence)
    o["product_role"] = f("reported product/result", "reported product/result", evidence)
    o["value"] = f(raw, norm, evidence)
    if yield_type is not None:
        o["yield_type"] = f(yield_type, yield_type, evidence)
    if basis is not None:
        o["measurement_basis"] = f(basis, basis, evidence)
    if scope is not None:
        o["scope"] = f(scope, scope, evidence)
    r["outcomes"].append(o)
    return o


def classification(r, specs, evidence, coverage="partial", limitation=None):
    changes = []
    labels = {}
    for idx, spec in enumerate(specs, 1):
        a, b, before, after, kind, text = spec
        change = {
            "change_id": f"b{idx}",
            "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
            "before_order": before, "after_order": after, "change_type": kind,
        }
        if kind.startswith("hydrogen_"):
            change["hydrogen_count"] = 1
        r["classification"]["bond_changes"].append(derived(
            text, change, evidence,
            "Category-level net change inferred from the source reactant/product drawings; atom identities are not mapped."))
        if kind in ("formed", "order_increased", "hydrogen_gain"):
            lab = st.bond_label(a, b, after)
            labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(change["change_id"])
    r["classification"]["labels"] = derived(
        "category-level source-structure comparison", list(labels.values()), evidence,
        "Labels are generated from the explicitly listed partial bond changes under bond-formation-20261016.")
    r["classification"]["coverage"] = inferred(
        "source drawings compared without atom mapping", coverage, evidence,
        "Complete atom mapping and hydrogen provenance are unavailable; coverage is therefore not complete.")
    r["classification"]["limitations"] = [
        limitation or "Partial category-level classification from source drawings; no complete atom mapping or stoichiometric hydrogen mapping.",
        "Model votes were not requested or fabricated.",
    ]
    r["classification"]["model_votes"] = []


def no_observed_classification(r, evidence, coverage="not_applicable"):
    r["classification"]["bond_changes"] = []
    r["classification"]["labels"] = inferred(
        "no observed product transformation", [], evidence,
        "No bond-formation category is assigned without an observed transformation.")
    r["classification"]["coverage"] = inferred(
        "no observed transformation", coverage, evidence,
        "Classification is not applicable to a no-product control." if coverage == "not_applicable" else "The item is a proposed mechanism, not an observed mapped reaction.")
    r["classification"]["limitations"] = ["No observed product bond changes are asserted."]
    r["classification"]["model_votes"] = []


MAIN_CHANGES = [
    ("C", "C", 0, 1, "formed", "new intramolecular C-C bond closes the six-membered spirocycle"),
    ("C", "H", 0, 1, "hydrogen_gain", "net C-H gain in allene hydrofunctionalization; proton source/atom mapping unresolved"),
]


def add_standard_system(r, substrate, product, ligand_label, temp, hours, evidence,
                        substrate_amount=None, substrate_mass=None):
    ps = participant(r, "p_substrate", substrate, "reactant", evidence)
    pp = participant(r, "p_product", product, "product", evidence, product_role="major_product")
    pc = participant(r, "p_rh", "Rh_cod_Cl_dimer", "catalyst", evidence)
    pl = participant(r, "p_ligand", ligand_label, "ligand", evidence)
    pa = participant(r, "p_ppts", "PPTS", "additive", evidence)
    solvent(r, "DCM", evidence, "DCM")
    set_charge(pc, "mol_percent", "2.5 mol%", {"number": 2.5, "unit": "mol%"}, evidence, substrate)
    set_charge(pl, "mol_percent", "7.5 mol%", {"number": 7.5, "unit": "mol%"}, evidence, substrate)
    set_charge(pa, "mol_percent", "50 mol%", {"number": 50, "unit": "mol%"}, evidence, substrate)
    if substrate_amount is not None:
        set_charge(ps, "amount", substrate_amount[0], substrate_amount[1], evidence)
        r["reaction_scale"] = f(substrate_amount[0], {"reference_molecule_id": mid(substrate), "amount": substrate_amount[1]}, evidence)
    if substrate_mass is not None:
        set_charge(ps, "mass", substrate_mass[0], substrate_mass[1], evidence)
        r["reaction_scale"] = f(substrate_mass[0], {"reference_molecule_id": mid(substrate), "mass": substrate_mass[1]}, evidence)
    reaction_step(r, "reaction", evidence, "rhodium-catalyzed intramolecular cyclization",
                  (f"{temp} °C", {"number": temp, "unit": "°C"}),
                  (f"{hours} h", {"number": hours, "unit": "h"}))
    classification(r, MAIN_CHANGES, evidence)
    return ps, pp


def add_yield_er_dr(r, product, y, er=None, dr=None, evidence=None,
                    yield_basis="isolated", er_basis="chiral HPLC", dr_basis="1H NMR analysis of the crude product"):
    evidence = evidence or r["evidence_ids"]
    outcome(r, "yield", "yield", product, f"{y}%", {"number": y, "unit": "%"}, evidence,
            basis=yield_basis, yield_type=yield_basis,
            scope={"from": "source substrate label", "to": product})
    if er is not None:
        outcome(r, "er", "selectivity", product, er,
                {"metric": "er", "ratio": er, "order": "reported major:minor enantiomer"}, evidence,
                basis=er_basis, scope={"metric": "er", "isomer_order": "reported major:minor; exact configurations follow source assignment"})
    if dr is not None:
        outcome(r, "dr", "selectivity", product, dr,
                {"metric": "dr", "ratio": dr, "order": "reported major:minor diastereomer(s)"}, evidence,
                basis=dr_basis, scope={"metric": "dr", "isomer_order": "drawn major product first; minor configuration(s) unassigned"})


# Scheme 1 previous-work literature examples.
for suffix, sub, prod, ppts, extra, outcomes in [
    ("scheme1A_previous_work", "litA_sub", "litA_prod", 50, None,
     [("yield_range", "yield", "76-99%", {"minimum": 76, "maximum": 99, "unit": "%"}),
      ("ee_range", "selectivity", "87-98% ee", {"metric": "ee", "minimum": 87, "maximum": 98, "unit": "%"})]),
    ("scheme1B_previous_work", "litB_sub", "litB_prod", 10, "HEH",
     [("yield_range", "yield", "71-95%", {"minimum": 71, "maximum": 95, "unit": "%"}),
      ("dr_range", "selectivity", ">20-8:1 dr", {"metric": "dr", "reported_range": ">20:1 to 8:1"}),
      ("ee_range", "selectivity", "19-98% ee", {"metric": "ee", "minimum": 19, "maximum": 98, "unit": "%"})]),
]:
    e = ev(2)
    r = new_reaction(suffix, 2, "literature_example", [], "Scheme 1 previous work", "previous work by the Breit group", "cited prior-work allene cyclization", e)
    participant(r, "p_substrate", sub, "reactant", e)
    participant(r, "p_product", prod, "product", e, product_role="reported product family")
    pc = participant(r, "p_rh", "Rh_cod_Cl_dimer", "catalyst", e)
    pl = participant(r, "p_ligand", "R_DTBM_Segphos", "ligand", e)
    pp = participant(r, "p_ppts", "PPTS", "additive", e)
    set_charge(pc, "mol_percent", "2.5 mol%", {"number": 2.5, "unit": "mol%"}, e, sub)
    set_charge(pl, "mol_percent", "7.5 mol%", {"number": 7.5, "unit": "mol%"}, e, sub)
    set_charge(pp, "mol_percent", f"{ppts} mol%", {"number": ppts, "unit": "mol%"}, e, sub)
    if extra:
        participant(r, "p_heh", extra, "reagent", e)
    solvent(r, "DCM", e, "DCM")
    reaction_step(r, "reaction", e, "reported previous-work cyclization", ("40 °C", {"number": 40, "unit": "°C"}), ("16 h", {"number": 16, "unit": "h"}))
    for oid, kind, raw, norm in outcomes:
        outcome(r, oid, kind, prod, raw, norm, e, scope={"reported_examples": 18})
    classification(r, MAIN_CHANGES, e)
    r["review"]["unresolved"].append("This is a cited prior-work range, not a newly performed experiment in the current paper.")


# Current-work generic arrows (Scheme 1 C and D).
for suffix, sub, prod, ligand, temp, tag in [
    ("scheme1C_generic_malonate", "fam3", "fam4", "R_DTBM_Segphos", 40, "C"),
    ("scheme1D_generic_tryptamine", "fam1", "fam2", "R_BINAP", 60, "D"),
]:
    e = ev(2)
    r = new_reaction(suffix, 2, "generic_scheme", [], f"Scheme 1 this work, arrow {tag}", "rhodium-catalyzed stereoselective cyclization", "generic intramolecular allene hydroarylation", e)
    add_standard_system(r, sub, prod, ligand, temp, 16, e)
    r["review"]["unresolved"].append("The adjacent 54-98%/dr/ee/14-example summary spans the combined current-work families and is not duplicated as a family-specific outcome.")


# Text-described initial experiment on 1a.
e = ev_both(2)
r = new_reaction("optimization_initial_1a", 2, "specific_experiment", ["optimization"],
                 "page 2 prose below Scheme 1", "cyclization of tosyl-protected model substrate 1a", "intramolecular allene hydroarylation", e)
add_standard_system(r, "1a", "2a", "R_DTBM_Segphos", 40, 24, e)
outcome(r, "yield", "yield", "2a", "70%", {"number": 70, "unit": "%"}, e,
        scope={"from": "1a", "to": "2a"})
outcome(r, "er", "selectivity", "2a", "82:18 er", {"metric": "er", "ratio": [82, 18]}, e,
        scope={"metric": "er", "order": "reported major:minor; exact configuration not restated for this optimization"})
outcome(r, "dr", "selectivity", "2a", "2:1 dr", {"metric": "dr", "ratio": [2, 1]}, e,
        scope={"metric": "dr", "order": "reported major:minor"})
r["procedure"]["general_procedure_reference"] = f("synthetic protocol from the tetrahydrocarbazole synthesis (Scheme 1A)", "Scheme 1A conditions inherited except the explicitly reported 24 h duration", e)
r["review"]["unresolved"].append("The 70% yield basis and er/dr measurement methods are not stated for this prose-only optimization result.")


# Scheme 2 scope.
scheme2 = [
    ("a", 90, [96, 4], [4, 1]), ("b", 89, [96.5, 3.5], [5, 1]),
    ("c", 92, [95, 5], [4, 1]), ("d", 97, [96.5, 3.5], [4, 1]),
    ("e", 82, [98, 2], [5, 1]), ("f", 71, [90, 10], [3, 1]),
    ("g", 85, [96, 4], [6, 1]), ("h", 87, [97.5, 2.5], [5, 1]),
    ("i", 95, None, [10, 1, 0, 0]),
]
for letter, y, eratio, dratio in scheme2:
    e = ev(2)
    r = new_reaction(f"scheme2_2{letter}", 2, "specific_experiment", ["substrate_scope"],
                     f"Scheme 2 entry 2{letter}", "rhodium-catalyzed cyclization of N-allenyltryptamine", "intramolecular allene hydroarylation", e)
    r["general_scheme_ids"] = [f"{DOC}:scheme1D_generic_tryptamine"]
    add_standard_system(r, f"1{letter}", f"2{letter}", "R_BINAP", 60, 16, e)
    add_yield_er_dr(r, f"2{letter}", y, eratio, dratio, e)


# Text-described lower-temperature 3a experiment.
e = ev_both(2)
r = new_reaction("optimization_3a_30C", 2, "specific_experiment", ["optimization"],
                 "page 2 prose after Scheme 2", "lower-temperature cyclization of model substrate 3a", "intramolecular allene hydroarylation", e)
r["general_scheme_ids"] = [f"{DOC}:scheme1C_generic_malonate"]
add_standard_system(r, "3a", "4a", "R_DTBM_Segphos", 30, 16, e)
add_yield_er_dr(r, "4a", 70, [95, 5], [13, 1], e)
r["review"]["unresolved"].append("The 16 h duration is inherited from the adjacent displayed catalytic protocol; the prose explicitly changes only temperature to 30 °C.")


# Scheme 3 scope.
scheme3 = [
    ("a", 99, [92, 8], [14, 1]), ("b", 63, [85, 15], [4, 1]),
    ("c", 76, [74, 26], [4, 1]), ("d", 90, [93, 7], [">20", 1]),
    ("e", 98, [95, 5], [16, 1]), ("f", 54, [94, 6], [9, 1]),
]
for letter, y, eratio, dratio in scheme3:
    e = ev(3)
    r = new_reaction(f"scheme3_4{letter}", 3, "specific_experiment", ["substrate_scope"],
                     f"Scheme 3 entry 4{letter}", "rhodium-catalyzed cyclization of malonate-tethered 3-allenylindole", "intramolecular allene hydroarylation", e)
    r["general_scheme_ids"] = [f"{DOC}:scheme1C_generic_malonate"]
    add_standard_system(r, f"3{letter}", f"4{letter}", "R_DTBM_Segphos", 40, 16, e)
    add_yield_er_dr(r, f"4{letter}", y, eratio, dratio, e)


# Table 1 robustness controls.
robustness = [
    (1, "benzaldehyde", 88, [7.5, 1], [98.5, 1.5]),
    (2, "benzoic_acid", 83, [5, 1], [93, 7]),
    (3, "4_bromopyrazole", 62, [6, 1], [75, 25]),
    (4, "ethanol", 83, [7.5, 1], [96.5, 3.5]),
    (5, "phenol", 84, [6.5, 1], [97, 3]),
    (6, "terminal_allene_additive", 29, [6.5, 1], [97, 3]),
    (7, "internal_alkyne_additive", 83, [6, 1], [97, 3]),
    (8, "TMS_diene", 71, [5, 1], [95, 5]),
    (9, None, 84, [9, 1], [95.5, 4.5]),
    (10, "water", 81, [6, 1], [96, 4]),
]
for entry, additive, y, dratio, eratio in robustness:
    e = ev(3)
    r = new_reaction(f"table1_entry{entry}", 3, "specific_experiment", ["control"],
                     f"Table 1 entry {entry}", "robustness screening of 1b cyclization", "intramolecular allene hydroarylation", e)
    add_standard_system(r, "1b", "2b", "R_BINAP", 60, 16, e,
                        substrate_amount=("0.2 mmol", {"number": 0.2, "unit": "mmol"}))
    if additive:
        pa = participant(r, "p_robustness_additive", additive, "additive", e)
        set_charge(pa, "equivalents", "1.0 equiv additive", {"number": 1.0, "unit": "equiv"}, e, "1b")
        r["steps"][0]["participant_ids"].append(pa["participant_id"])
    else:
        r["procedure"]["modifications"] = f("non flame-dried flask", "non-flame-dried flask", e)
    add_yield_er_dr(r, "2b", y, eratio, dratio, e,
                    yield_basis="isolated", er_basis="chiral HPLC",
                    dr_basis="1H NMR analysis of the isolated compound")
    if entry == 7:
        outcome(r, "alkyne_recovery", "recovery", "internal_alkyne_additive", "70% unreacted alkyne",
                {"number": 70, "unit": "%"}, e, basis="crude 1H NMR",
                scope={"entity": "prop-1-yn-1-ylbenzene additive", "state": "unreacted"})


# Proposed mechanism and transition-state rationale (not an observed experiment).
e = [f"{DOC}:page2:text", f"{DOC}:page3:image", f"{DOC}:page3:text"]
r = new_reaction("proposed_mechanism", 3, "proposed_mechanism", [],
                 "pages 2-3 mechanistic prose and Figure 2", "postulated rhodium-hydride/allyl-rhodium cyclization pathway", "author-proposed mechanism", e)
r["procedure"]["original_text"] = f(
    "formation of a rhodium(III) hydride with PPTS, hydrometalation of the allene, and nucleophilic attack by indole C3; chairlike TS proposed for the six-membered product",
    "author-proposed pathway and transition-state rationale", e)
no_observed_classification(r, e, coverage="unresolved")
r["review"]["unresolved"].append("The pathway and transition states are author proposals, not atom-mapped or experimentally observed elementary steps.")


# Scheme 4 gram-scale experiment.
e = ev(4)
r = new_reaction("scheme4_gram_scale_2b", 4, "specific_experiment", ["scale_up"],
                 "Scheme 4 gram-scale vertical arrow", "gram-scale rhodium-catalyzed cyclization of 1b", "intramolecular allene hydroarylation", e)
r["general_scheme_ids"] = [f"{DOC}:scheme1D_generic_tryptamine"]
add_standard_system(r, "1b", "2b", "R_BINAP", 60, 16, e,
                    substrate_mass=("1 g", {"number": 1, "unit": "g"}))
add_yield_er_dr(r, "2b", 87, [95.5, 4.5], [5, 1], e,
                yield_basis="isolated", er_basis="chiral HPLC",
                dr_basis="1H NMR analysis of the isolated product")


def transformation_record(suffix, substrate, product, operation, evidence, solvent_label,
                          temp_raw, temp_norm, time_raw, time_norm, changes):
    r = new_reaction(suffix, 4, "specific_experiment", ["derivatization"],
                     f"Scheme 4 arrow {suffix.split('_')[-1]}", operation, operation, evidence)
    participant(r, "p_substrate", substrate, "reactant", evidence)
    participant(r, "p_product", product, "product", evidence, product_role="major_product")
    solvent(r, solvent_label, evidence, solvent_label)
    reaction_step(r, "reaction", evidence, operation,
                  (temp_raw, temp_norm) if temp_raw is not None else None,
                  (time_raw, time_norm) if time_raw is not None else None)
    classification(r, changes, evidence)
    return r


# (a) selective imine reduction.
e = ev(4)
r = transformation_record("scheme4_a", "2b", "5", "selective reduction of the imine moiety", e, "MeOH",
                          None, None, "16 h", {"number": 16, "unit": "h"},
                          [("C", "N", 2, 1, "order_decreased", "C=N bond reduced to C-N"),
                           ("C", "H", 0, 1, "hydrogen_gain", "net C-H gain at the reduced imine carbon"),
                           ("N", "H", 0, 1, "hydrogen_gain", "net N-H gain in the reduced product")])
pr = participant(r, "p_NaBH3CN", "NaBH3CN", "reagent", e)
set_charge(pr, "equivalents", "1.0 equiv", {"number": 1.0, "unit": "equiv"}, e, "2b")
r["steps"][0]["participant_ids"].append(pr["participant_id"])
r["steps"][0]["temperature"] = f("0 °C to rt", None, e, status="unreadable")
r["steps"][0]["operation"] = f("NaBH3CN, MeOH, 0 °C to rt, 16 h", "reaction starts at 0 °C and warms to room temperature; stage durations unreported", e)
add_yield_er_dr(r, "5", 93, [97, 3], [6, 1], e,
                yield_basis="isolated", er_basis="chiral HPLC",
                dr_basis="1H NMR analysis of the isolated product")
r["review"]["unresolved"].append("The source gives an overall 0 °C-to-rt/16 h program but does not allocate addition, warming, or hold times; no stage duration was invented.")


# (b) hydrogenation of 2b and 2e are separate experiments.
for product, substrate, y, eratio, dratio in [
    ("6a", "2b", 94, [97.5, 2.5], [6, 1]),
    ("6b", "2e", 98, [98.5, 1.5], [8, 1]),
]:
    r = transformation_record(f"scheme4_b_{product}", substrate, product,
                              "Pd/C hydrogenation of imine and terminal alkene", e, "MeOH",
                              "rt", "rt", "16 h", {"number": 16, "unit": "h"},
                              [("C", "N", 2, 1, "order_decreased", "imine C=N bond reduced"),
                               ("C", "C", 2, 1, "order_decreased", "terminal C=C bond reduced"),
                               ("C", "H", 0, 1, "hydrogen_gain", "net C-H gain during hydrogenation"),
                               ("N", "H", 0, 1, "hydrogen_gain", "net N-H gain during hydrogenation")])
    pcat = participant(r, "p_PdC", "Pd_C", "catalyst", e)
    set_charge(pcat, "mol_percent", "10 mol % Pd (Pd/C, 5 wt %)", {"number": 10, "unit": "mol%"}, e, substrate)
    ph2 = participant(r, "p_H2", "H2", "reagent", e)
    r["steps"][0]["participant_ids"] += [pcat["participant_id"], ph2["participant_id"]]
    r["steps"][0]["atmosphere"] = f("H2 (via balloon)", "hydrogen atmosphere via balloon", e)
    r["steps"][0]["pressure"] = f("1 atm", "1 atm", e)
    r["notes"].append(f("H2 (1 atm, via balloon)", "Pressure is preserved as source text because atm is not a supported structured quantity unit in the fixed schema/tooling.", e))
    add_yield_er_dr(r, product, y, eratio, dratio, e,
                    yield_basis="isolated", er_basis="chiral HPLC",
                    dr_basis="1H NMR analysis of the isolated product")
    if product == "6a":
        r["notes"].append(f("ABS = 4-aminobenzenesulfonyl", "p-nosyl nitro group is also reduced to the corresponding amine", e))
    else:
        r["notes"].append(f("6b was synthesized from 2e; R = MBS", "MBS substrate 2e avoids p-nosyl nitro-group reduction", e))


# (c) Pinnick oxidation.
r = transformation_record("scheme4_c", "2b", "7", "Pinnick oxidation to spirooxindole", e, "THF",
                          "rt", "rt", "7 h", {"number": 7, "unit": "h"},
                          [("C", "O", 0, 2, "formed", "new oxindole carbonyl C=O bond"),
                           ("C", "N", 2, 1, "order_decreased", "indolenine C=N becomes lactam C-N")])
for pid, label, raw, norm in [
    ("p_scavenger", "2methyl2butene", "10 equiv", {"number": 10, "unit": "equiv"}),
    ("p_NaClO2", "NaClO2", "5.0 equiv", {"number": 5.0, "unit": "equiv"}),
    ("p_NaH2PO4", "NaH2PO4_aq", "1.5 equiv (aq. 1 M)", {"number": 1.5, "unit": "equiv"}),
]:
    p = participant(r, pid, label, "reagent" if label != "NaH2PO4_aq" else "additive", e)
    set_charge(p, "equivalents", raw, norm, e, "2b")
    r["steps"][0]["participant_ids"].append(p["participant_id"])
r["notes"].append(f("NaH2PO4 (aq. 1 M)", "NaH2PO4 was supplied as a 1 M aqueous solution; its solution volume was not reported", e))
add_yield_er_dr(r, "7", 98, [97.5, 2.5], [6, 1], e,
                yield_basis="isolated", er_basis="chiral HPLC",
                dr_basis="1H NMR analysis of the isolated product")


# (d) p-nosyl deprotection.
r = transformation_record("scheme4_d", "2b", "8", "p-nosyl deprotection", e, "DMF",
                          "50 °C", {"number": 50, "unit": "°C"}, "16 h", {"number": 16, "unit": "h"},
                          [("N", "S", 1, 0, "broken", "protecting-group N-S bond cleaved"),
                           ("N", "H", 0, 1, "hydrogen_gain", "N-H bond present in deprotected product")])
for pid, label, raw, norm in [
    ("p_thiocresole", "p_thiocresole", "1.5 equiv", {"number": 1.5, "unit": "equiv"}),
    ("p_K2CO3", "K2CO3", "3.0 equiv", {"number": 3.0, "unit": "equiv"}),
]:
    p = participant(r, pid, label, "reagent", e)
    set_charge(p, "equivalents", raw, norm, e, "2b")
    r["steps"][0]["participant_ids"].append(p["participant_id"])
add_yield_er_dr(r, "8", 69, [97.5, 2.5], [5, 1], e,
                yield_basis="isolated", er_basis="chiral HPLC after conversion to derivative 2a",
                dr_basis="1H NMR analysis of the isolated product")
r["notes"].append(f("er determined after conversion to derivative 2a", "analytical derivatization conditions are unreported and therefore are not fabricated as a separate preparative reaction", e))


# (e) and (f) are distinct no-product controls.
for suffix, reagent_label, solvent_label, temp_raw, temp_norm, time_h, extra in [
    ("scheme4_e", "TsOH", "THF", "rt to 60 °C", None, 4, "TsOH (30 mol %)"),
    ("scheme4_f", "HCl", "THF", "rt", "rt", 16, "HCl(g) in THF (4.5 M)"),
]:
    r = new_reaction(suffix, 4, "specific_experiment", ["control"],
                     f"Scheme 4 arrow {suffix[-1]}", "attempted acid-catalyzed migration", "no observed migration product", e)
    participant(r, "p_substrate", "2b", "reactant", e)
    preag = participant(r, "p_acid", reagent_label, "reagent", e)
    solvent(r, solvent_label, e, "THF")
    if reagent_label == "TsOH":
        set_charge(preag, "mol_percent", "30 mol %", {"number": 30, "unit": "mol%"}, e, "2b")
    reaction_step(r, "reaction", e, f"{extra}; no product",
                  (temp_raw, temp_norm) if temp_norm is not None else None,
                  (f"{time_h} h", {"number": time_h, "unit": "h"}))
    if temp_norm is None:
        r["steps"][0]["temperature"] = f(temp_raw, None, e, status="unreadable")
        r["review"]["unresolved"].append("The source gives rt-to-60 °C over an overall 4 h but does not allocate a ramp or hold duration.")
    if reagent_label == "HCl":
        r["solvent_conditions"]["composition"] = f("HCl (g) in THF (4.5 M)", [mid("HCl"), mid("THF")], e)
        r["notes"].append(f("HCl (g) in THF (4.5 M)", "HCl solution concentration reported; charge/volume unreported", e))
    outcome(r, "no_product", "no_reaction", [], "no product", "no product", e,
            basis="visual/isolated-product reporting basis not stated", scope={"attempted_target": "annulated 7-membered heterocycle"})
    if reagent_label == "TsOH":
        outcome(r, "substrate_recovery", "recovery", "2b", "96% of 2b isolated", {"number": 96, "unit": "%"}, e,
                basis="isolated", scope={"entity": "unchanged 2b"})
    no_observed_classification(r, e)


bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": f"{DOC}:context_page1", "type": "article_context", "evidence_ids": ev_both(1),
     "note": "Title, abstract, duplicate graphical summary, introduction, and Figure 1 natural products; no current-paper experimental arrow counted from this page."},
    {"context_id": f"{DOC}:context_scheme2_xray", "type": "stereochemical_assignment", "evidence_ids": ev_both(2),
     "note": "2a absolute configuration is (S,S) by X-ray; other 2-series compounds are assigned by analogy. Relative configuration was first compared with literature 1H NMR data."},
    {"context_id": f"{DOC}:context_scheme3_crystals", "type": "crystallography", "evidence_ids": ev_both(3),
     "note": "Centrosymmetric (rac) ORTEP structures are shown for 4b and 4c as support for relative diastereochemistry; this is not a separate synthesis or a basis to relabel enantioenriched scope products as racemic."},
    {"context_id": f"{DOC}:context_figure2", "type": "mechanistic_discussion", "evidence_ids": ev_both(3),
     "note": "Figure 2 compares proposed chairlike and envelope transition-state geometries; it is an author hypothesis, not an observed experiment."},
    {"context_id": f"{DOC}:context_scheme4_optical_purity", "type": "collective_observation", "evidence_ids": ev_both(3),
     "note": "The prose states no loss of optical purity in the transformations; individual reported er values remain attached to their specific Scheme 4 products."},
    {"context_id": f"{DOC}:context_si_notice", "type": "supplement_notice", "evidence_ids": ev_both(4),
     "note": "The paper announces experimental details and analytical data in Supporting Information, but no SI file was supplied or explicitly associated in this session."},
    {"context_id": f"{DOC}:context_page5", "type": "bibliography", "evidence_ids": ev_both(5),
     "note": "Page 5 contains references only; no new current-paper experiment is extracted."},
]

st.write(RUN / "submitted.json", bundle)
st.write(RUN / "raw_model/stereo_diagnostics.json", stereo_diagnostics)
st.write(RUN / "raw_model/transcription.json", {
    "document_id": DOC,
    "pages_inspected": [1, 2, 3, 4, 5],
    "reaction_count": len(reactions),
    "molecule_count": len(molecules),
    "scheme2": scheme2,
    "scheme3": scheme3,
    "table1": robustness,
    "notes": [
        "Machine transcription; all formal records remain needs_review.",
        "No SI was associated.",
        "Dense source-drawn structures were not forced into guessed complete SMILES.",
    ],
})
print(json.dumps({"molecules": len(molecules), "reactions": len(reactions),
                  "smiles_diagnostics": len(stereo_diagnostics)}, ensure_ascii=False))
