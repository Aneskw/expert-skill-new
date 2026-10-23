#!/usr/bin/env python3
"""Build the evidence-linked v3.1 submission for Li, Jandl, and Bach (2020).

This is a source transcription helper, not an independent chemistry validator.
Dense product structures remain unresolved rather than being forced into unsafe
SMILES.  The generated bundle is intentionally machine-review status.
"""

from __future__ import annotations

import copy
import json
import sys
from pathlib import Path


WORKSPACE = Path(__file__).resolve().parents[3]
SCRIPTS = WORKSPACE / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))

import standard_records as st  # noqa: E402


RUN = WORKSPACE / "runs/current"
DATASET_PATH = RUN / "dataset.json"
DOC = "doc_349d3090350bfe00"
EV = {page: f"{DOC}:page{page}:image" for page in range(1, 6)}
TXT = {page: f"{DOC}:page{page}:text" for page in range(1, 6)}
TITLE = "Visible-Light-Mediated Enantioselective Photoreactions of 3-Alkylquinolones with 4-O-Tethered Alkenes and Allenes"


def f(raw=None, normalized=None, evidence=(), method="original_extraction", derivation=None,
      status=None):
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation, status=status)


def missing(evidence=()):
    return f(evidence=evidence)


def derived(normalized, evidence, derivation):
    return f(None, normalized, evidence, method="structure_derivation", derivation=derivation)


def inferred(normalized, evidence, derivation):
    return f(None, normalized, evidence, method="inference", derivation=derivation)


def q(number, unit, raw_text=None, qualifier=None):
    value = {"number": number, "unit": unit}
    if raw_text is not None:
        value["raw_text"] = raw_text
    if qualifier is not None:
        value["qualifier"] = qualifier
    return value


def range_q(minimum, maximum, unit, raw_text):
    return {"minimum": minimum, "maximum": maximum, "unit": unit, "raw_text": raw_text}


def set_status_fact(container, key, value, evidence, rationale, source=False):
    if source:
        container[key] = f(value, value, evidence)
    else:
        container[key] = inferred(value, evidence, rationale)


def base_molecule(mid, label, evidence, name=None):
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    m["names"]["original"] = f(name or label, name or label, evidence)
    m["names"]["english"] = missing(evidence)
    m["names"]["abbreviations"] = missing(evidence)
    m["names"]["aliases"] = missing(evidence)
    m["cas_numbers"] = missing(evidence)
    m["external_properties"] = []
    return m


STEREO_CHECKS = {}


def smiles_molecule(mid, label, smiles, evidence, name, identity_note):
    m = base_molecule(mid, label, evidence, name)
    m["representation_type"] = "smiles"
    m["smiles"] = f({"source": "structure drawing", "candidate_smiles": smiles}, smiles, evidence)
    m["generic_template"] = None
    diag = st.stereo_inspect(smiles)
    STEREO_CHECKS[mid] = diag
    m["stereochemistry"]["atom_centers"] = derived(
        diag["atom_centers"], evidence,
        "RDKit graph diagnostic on the transcribed candidate SMILES; not source verification."
    )
    m["stereochemistry"]["double_bonds"] = derived(
        diag["double_bonds"], evidence,
        "RDKit graph diagnostic on the transcribed candidate SMILES; aromatic bonds excluded from E/Z claims."
    )
    m["stereochemistry"]["axial"] = inferred(
        None, evidence,
        "Allene axial assignment cannot be certified by the toolkit; source drawing was inspected separately."
    )
    m["stereochemistry"]["mixture"] = missing(evidence)
    m["stereochemistry"]["assignment"] = inferred(
        {"candidate_graph": "no source stereodescriptor is printed for this input structure"}, evidence,
        "The candidate graph has no encoded tetrahedral or E/Z descriptor; absence of source descriptors is not treated as racemic."
    )
    set_status_fact(m["stereochemistry"], "chirality_status", "achiral", evidence,
                    "Candidate graph and source drawing have no stereogenic center; allenes here have a terminal CH2 or duplicate substituents.")
    set_status_fact(m["stereochemistry"], "geometry_status", "not_applicable", evidence,
                    "No source-relevant E/Z element is present in the candidate graph.")
    set_status_fact(m["stereochemistry"], "representation_status", "not_applicable", evidence,
                    "No stereochemical encoding is required for the transcribed candidate graph.")
    m["stereochemistry"]["isomeric_smiles"] = derived(
        diag["canonical_isomeric_smiles"], evidence,
        "Canonical isomeric SMILES generated from the candidate graph; this does not verify transcription."
    )
    m["stereochemistry"]["identity_key"] = derived(
        diag["stereo_graph_key"], evidence,
        "Stereo graph key generated from the candidate SMILES for within-bundle comparison only."
    )
    m["identity_note"] = f(identity_note, identity_note, evidence)
    return m


def unresolved_molecule(mid, label, evidence, note, chirality="unresolved", geometry="unresolved",
                        representation="unresolved", atom_centers=None, axial=None, assignment=None,
                        mixture=None, name=None):
    m = base_molecule(mid, label, evidence, name)
    m["representation_type"] = "unresolved"
    m["smiles"] = missing(evidence)
    m["generic_template"] = None
    m["identity_note"] = f(note, note, evidence)
    m["stereochemistry"]["atom_centers"] = (
        f(atom_centers, atom_centers, evidence) if atom_centers is not None else missing(evidence)
    )
    m["stereochemistry"]["double_bonds"] = missing(evidence)
    m["stereochemistry"]["axial"] = f(axial, axial, evidence) if axial is not None else missing(evidence)
    m["stereochemistry"]["mixture"] = f(mixture, mixture, evidence) if mixture is not None else missing(evidence)
    m["stereochemistry"]["assignment"] = (
        f(assignment, assignment, evidence) if assignment is not None else missing(evidence)
    )
    set_status_fact(m["stereochemistry"], "chirality_status", chirality, evidence,
                    "Status retained from the source drawing/label without inventing R/S descriptors.", source=True)
    set_status_fact(m["stereochemistry"], "geometry_status", geometry, evidence,
                    "Status retained from the source drawing/label without inventing E/Z descriptors.", source=True)
    set_status_fact(m["stereochemistry"], "representation_status", representation, evidence,
                    "Dense source structure was not forced into an unverified stereochemical SMILES.", source=True)
    m["stereochemistry"]["isomeric_smiles"] = missing(evidence)
    m["stereochemistry"]["identity_key"] = missing(evidence)
    return m


def template_molecule(mid, label, evidence, fixed, variables, reason, chirality="unresolved",
                      geometry="unresolved", representation="not_encoded"):
    m = base_molecule(mid, label, evidence)
    m["representation_type"] = "generic_template"
    m["smiles"] = missing(evidence)
    template_vars = []
    for variable in variables:
        template_vars.append({
            "label": variable["label"],
            "attachment_points": f(variable.get("attachment"), variable.get("attachment"), evidence),
            "definition": f(variable.get("definition"), variable.get("definition"), evidence),
            "allowed_values": f(variable.get("allowed"), variable.get("allowed"), evidence),
            "constraints": f(variable.get("constraints"), variable.get("constraints"), evidence),
        })
    m["generic_template"] = {
        "fixed_structure": f(fixed, fixed, evidence),
        "variables": template_vars,
        "instance_mapping": [],
        "reason": f(reason, reason, evidence),
    }
    m["identity_note"] = f(reason, reason, evidence)
    m["stereochemistry"]["atom_centers"] = missing(evidence)
    m["stereochemistry"]["double_bonds"] = missing(evidence)
    m["stereochemistry"]["axial"] = missing(evidence)
    m["stereochemistry"]["mixture"] = missing(evidence)
    m["stereochemistry"]["assignment"] = (
        f(
            "Generic product drawing uses wedge/dash bonds to specify relative stereochemistry; exact R/S descriptors vary by instance and are not printed.",
            {"source": "wedge/dash bonds in generic product drawing", "R/S": "not printed"},
            evidence,
        )
        if chirality == "specified" else missing(evidence)
    )
    set_status_fact(m["stereochemistry"], "chirality_status", chirality, evidence, "Generic R-group structure.", source=True)
    set_status_fact(m["stereochemistry"], "geometry_status", geometry, evidence, "Generic R-group structure.", source=True)
    set_status_fact(m["stereochemistry"], "representation_status", representation, evidence, "Generic R-group structure.", source=True)
    m["stereochemistry"]["isomeric_smiles"] = missing(evidence)
    m["stereochemistry"]["identity_key"] = missing(evidence)
    return m


def base_reaction(suffix, evidence, record_type="specific_experiment", categories=None,
                  location=None, reaction_name=None, transformation=None):
    rid = f"{DOC}:{suffix}"
    r = st.new_reaction(rid, DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = (
        f(categories, categories, evidence) if categories is not None else missing(evidence)
    )
    r["location"] = f(location, location, evidence)
    r["reaction_name"] = f(reaction_name, reaction_name, evidence) if reaction_name else missing(evidence)
    r["transformation_type"] = (
        inferred(transformation, evidence, "Normalized from the source reaction title/drawing.")
        if transformation else missing(evidence)
    )
    r["review"] = {
        "status": "needs_review",
        "method": "model",
        "reviewer": None,
        "checked_fields": [],
        "unresolved": [],
    }
    r["completeness"] = "partial"
    return r


def participant(pid, mid, role, evidence, mol_percent=None, reference=None, stage="reaction",
                product_role=None):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    if product_role is not None:
        p["product_role"] = f(product_role, product_role, evidence)
    if mol_percent is not None:
        raw = f"{mol_percent:g} mol %"
        p["charges"]["mol_percent"] = f(raw, q(mol_percent, "mol%", raw), evidence)
        p["equivalence_reference"] = f(reference, reference, evidence)
    return p


def reaction_step(step_id, participants, evidence, temperature=None, time=None,
                  operation="visible-light irradiation"):
    s = st.new_step(step_id, evidence)
    s["participant_ids"] = participants
    s["operation"] = f(operation, operation, evidence)
    if temperature is not None:
        raw, value = temperature
        s["temperature"] = f(raw, value, evidence)
    if time is not None:
        raw, value, method = time
        s["time"] = f(raw, value, evidence, method=method)
    return s


def solvent_fields(r, composition, evidence, concentration=None, concentration_basis=None):
    r["solvent_conditions"]["composition"] = f(composition, composition, evidence)
    if concentration is not None:
        raw, value = concentration
        r["solvent_conditions"]["concentration"] = f(raw, value, evidence)
    if concentration_basis is not None:
        r["solvent_conditions"]["concentration_basis"] = f(
            concentration_basis, concentration_basis, evidence
        )


def light_fields(r, evidence, source="hν", wavelength=420):
    r["special_conditions"]["light_source"] = f(source, source, evidence)
    r["special_conditions"]["wavelength"] = f(
        f"λ = {wavelength} nm", q(wavelength, "nm", f"λ = {wavelength} nm"), evidence
    )


def product_outcome(oid, kind, target, raw, value, evidence, role="main product",
                    yield_type=None, basis=None, scope=None):
    o = st.new_outcome(oid, kind, [target] if isinstance(target, str) else target, evidence)
    o["product_role"] = f(role, role, evidence)
    o["value"] = f(raw, value, evidence)
    if yield_type is not None:
        o["yield_type"] = f(yield_type, yield_type, evidence)
    if basis is not None:
        o["measurement_basis"] = f(basis, basis, evidence)
    if scope is not None:
        o["scope"] = f(scope, scope, evidence)
    return o


def yield_and_ee(prefix, product_mid, yield_value, ee_value, evidence, yield_basis=None,
                 role="main product", ee_note="major enantiomer corresponds to the drawn product"):
    return [
        product_outcome(
            f"{prefix}_yield", "yield", product_mid, f"{yield_value}%", q(yield_value, "%", f"{yield_value}%"),
            evidence, role=role, basis=yield_basis,
            scope={"span": "substrate to named product", "product": product_mid},
        ),
        product_outcome(
            f"{prefix}_ee", "selectivity", product_mid, f"{ee_value}% ee", q(ee_value, "%", f"{ee_value}% ee"),
            evidence, role=role, scope={"metric": "ee", "assignment": ee_note},
        ),
    ]


def set_partial_cycloaddition_classification(r, evidence, formed_bonds=2, note=None):
    changes = []
    supporting = []
    for index in range(1, formed_bonds + 1):
        cid = f"b{index}"
        supporting.append(cid)
        change = {
            "change_id": cid,
            "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
            "before_order": 0,
            "after_order": 1,
            "change_type": "formed",
        }
        changes.append(derived(
            change, evidence,
            "The displayed intramolecular [2+2] product has a cyclobutane made by two new C-C sigma connections; exact atom maps are unresolved."
        ))
    label = [{"category_number": 2, "code": "C-C", "supporting_change_ids": supporting}]
    r["classification"]["bond_changes"] = changes
    r["classification"]["labels"] = derived(
        label, evidence,
        "Taxonomy label follows the source-drawn new cyclobutane connectivity; no complete atom mapping is asserted."
    )
    r["classification"]["coverage"] = inferred(
        "partial", evidence,
        "Only the unambiguous new C-C connections are classified; bond-order decreases, hydrogens, and atom maps are unresolved."
    )
    r["classification"]["limitations"] = [
        "Exact atom mapping and all bond-order/hydrogen changes are unresolved.",
        note or "Classification is structure-derived from the compact source drawing and remains pending review.",
    ]


def set_unresolved_classification(r, evidence, reason, not_applicable=False):
    r["classification"]["bond_changes"] = []
    r["classification"]["labels"] = inferred(
        [], evidence, "No bond-formation label is asserted for this non-observed or unresolved transformation."
    )
    r["classification"]["coverage"] = inferred(
        "not_applicable" if not_applicable else "unresolved", evidence, reason
    )
    r["classification"]["limitations"] = [reason]


def add_common_photoreaction(r, substrate_mid, product_mids, evidence, catalyst_loading,
                             solvent_mids, solvent_text, temperature, time,
                             wavelength=420, light_source="hν", operation="visible-light irradiation"):
    substrate_pid = f"p_{substrate_mid}"
    r["participants"].append(participant(substrate_pid, substrate_mid, "reactant", evidence))
    for index, product_mid in enumerate(product_mids, 1):
        r["participants"].append(participant(
            f"p_{product_mid}", product_mid, "product", evidence,
            product_role="main product" if len(product_mids) == 1 else f"product branch {index}"
        ))
    r["participants"].append(participant(
        "p_catalyst6", "mol_catalyst_6", "catalyst", evidence,
        mol_percent=catalyst_loading, reference=substrate_mid
    ))
    for solvent_mid in solvent_mids:
        r["participants"].append(participant(f"p_{solvent_mid}", solvent_mid, "solvent", evidence))
    step_pids = [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"]
    r["steps"] = [reaction_step(
        "irradiation", step_pids, evidence, temperature=temperature, time=time,
        operation=operation
    )]
    solvent_fields(r, solvent_text, evidence)
    light_fields(r, evidence, source=light_source, wavelength=wavelength)


def main():
    data = json.loads(DATASET_PATH.read_text())
    if any(r.get("status") == "reviewed" for r in data.get("experiments", [])):
        raise RuntimeError("Refusing to rebuild when source-reviewed records are present.")
    doc = next(d for d in data["documents"] if d["document_id"] == DOC)
    doc["title"] = TITLE
    doc["supplementary_material"] = {
        "status": "not_provided",
        "note": "The article states that Supporting Information is available, but no SI file was explicitly supplied or associated in this session."
    }
    DATASET_PATH.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n")
    bundle = st.empty(copy.deepcopy(data["documents"]))

    # Concrete candidate structures are limited to clearly drawn inputs and simple media.
    bundle["molecules"].extend([
        smiles_molecule("mol_5a", "5a", "O=C1NC2=CC=CC=C2C(OCCC=C)=C1", [EV[2]],
                        "4-(but-3-enyloxy)quinolin-2(1H)-one",
                        "Candidate SMILES transcribed from Scheme 1; no stereochemical descriptor is printed."),
        smiles_molecule("mol_5b", "5b", "CC1=C(OCCC=C)c2ccccc2NC1=O", [EV[2]],
                        "4-(but-3-enyloxy)-3-methylquinolin-2(1H)-one",
                        "Candidate SMILES transcribed from Scheme 1 and the accompanying prose."),
        smiles_molecule("mol_phcf3", "PhCF3", "FC(F)(F)c1ccccc1", [EV[2], EV[3], EV[4]],
                        "trifluorotoluene", "Solvent identity expanded from the printed PhCF3 abbreviation."),
        smiles_molecule("mol_hfx", "HFX", "FC(F)(F)c1cccc(C(F)(F)F)c1", [EV[3]],
                        "hexafluoro-m-xylene", "Solvent identity is explicitly expanded in the Scheme 4 title."),
        smiles_molecule("mol_8", "8", "CC1=C(OCCC(F)=C(F)F)c2ccccc2NC1=O", [EV[3]],
                        "4-(3,4,4-trifluorobut-3-enyloxy)-3-methylquinolin-2(1H)-one",
                        "Candidate SMILES transcribed from Scheme 4 and its printed name."),
        smiles_molecule("mol_10a", "10a", "CC1=C(OCCC=C=C)c2ccccc2NC1=O", [EV[3]],
                        "allene substrate 10a", "Candidate connectivity transcribed from Scheme 5 (R=H, X=H, R1=Me)."),
        smiles_molecule("mol_10b", "10b", "CC1=C(OCCC(C)=C=C)c2ccccc2NC1=O", [EV[3]],
                        "allene substrate 10b", "Candidate connectivity transcribed from Scheme 5 (R=Me, X=H, R1=Me)."),
        smiles_molecule("mol_10c", "10c", "CC1=C(OCCC(C)=C=C)c2cc(Cl)ccc2NC1=O", [EV[3]],
                        "allene substrate 10c", "Candidate connectivity transcribed from Scheme 5 (R=Me, X=Cl, R1=Me); benzo-ring position follows the drawing."),
        smiles_molecule("mol_10d", "10d", "CCC1=C(OCCC(C)=C=C)c2ccccc2NC1=O", [EV[3]],
                        "allene substrate 10d", "Candidate connectivity transcribed from Scheme 5 (R=Me, X=H, R1=Et)."),
        smiles_molecule("mol_12a", "12a", "CC1=C(OCC=C=C)c2cc(Cl)ccc2NC1=O", [EV[4]],
                        "6-chloro-4-(buta-2,3-dienyloxy)-3-methylquinolin-2(1H)-one",
                        "Candidate SMILES transcribed from the upper arrow of Scheme 6; benzo-ring position follows the drawing."),
        smiles_molecule("mol_15", "15", "CC1=C(OCC=C=C(C)C)c2ccccc2NC1=O", [EV[4]],
                        "terminally dimethyl-substituted allene 15",
                        "Candidate SMILES transcribed from the lower arrow of Scheme 6."),
        smiles_molecule("mol_h2so4", "H2SO4", "OS(=O)(=O)O", [EV[3], EV[4]],
                        "sulfuric acid", "Acid identity is printed in the page-3 ring-opening discussion."),
    ])

    bundle["molecules"].append(unresolved_molecule(
        "mol_catalyst_6", "6", [EV[2], EV[3], EV[4]],
        "Chiral thioxanthone sensitizer 6 with a lactam hydrogen-bonding motif; dense polycyclic structure retained from the source drawing.",
        chirality="specified", geometry="not_applicable", representation="not_encoded",
        atom_centers={"source": "wedge/dash stereochemistry in catalyst drawing", "absolute_descriptors": "not printed"},
        assignment="Chiral catalyst drawing specifies relative/absolute spatial arrangement, but no R/S labels are transcribed."
    ))

    # Scheme 2 substrates 5c-5m are defined by the generic drawing and their product panels;
    # no unsafe concrete SMILES are supplied from the compact product-only variants.
    for label in ["5c", "5d", "5e", "5f", "5g", "5h", "5i", "5j", "5k", "5l", "5m"]:
        bundle["molecules"].append(unresolved_molecule(
            f"mol_{label}", label, [EV[2]],
            f"Quinolone substrate {label}, defined by the generic substrate 5 and the corresponding product panel in Scheme 2; compact variable assignment not forced into SMILES.",
            chirality="unspecified", geometry="unspecified", representation="unresolved"
        ))
    bundle["molecules"].append(unresolved_molecule(
        "mol_5n", "5n", [EV[3]],
        "Fluoro-substituted 3-methylquinolone 5n bearing a methyl-substituted terminal alkene tether, as drawn in Scheme 3.",
        chirality="unspecified", geometry="not_applicable", representation="unresolved"
    ))

    assignment_analogy = (
        "The drawn major stereoisomer is specified by wedge/dash bonds. Product 7j was assigned by anomalous single-crystal X-ray diffraction; "
        "the absolute configurations of the other Scheme 2 products were assigned by analogy."
    )
    for label in ["7b", "7c", "7d", "7e", "7f", "7g", "7h", "7i", "7j", "7k", "7l", "7m"]:
        bundle["molecules"].append(unresolved_molecule(
            f"mol_{label}", label, [EV[2]],
            f"Fused cyclobutane product {label} as drawn in Scheme 2; dense ring stereochemistry retained textually rather than encoded in SMILES.",
            chirality="specified", geometry="not_applicable", representation="not_encoded",
            atom_centers={"source": "multiple wedge/dash-defined ring stereocenters", "R/S": "not transcribed"},
            assignment=assignment_analogy
        ))
    bundle["molecules"].append(unresolved_molecule(
        "mol_7n", "7n", [EV[3]],
        "Fluoro-substituted fused cyclobutane product 7n as drawn in Scheme 3.",
        chirality="specified", geometry="not_applicable", representation="not_encoded",
        atom_centers={"source": "multiple wedge/dash-defined ring stereocenters", "R/S": "not printed"},
        assignment="Major stereoisomer shown with wedge/dash bonds; absolute R/S labels are not printed in Scheme 3."
    ))
    bundle["molecules"].append(unresolved_molecule(
        "mol_9", "9", [EV[3]],
        "Trifluorinated fused cyclobutane 9 as drawn in Scheme 4.",
        chirality="specified", geometry="not_applicable", representation="not_encoded",
        atom_centers={"source": "wedge/dash-defined cyclobutane stereocenters", "R/S": "not printed"},
        assignment="Major enantiomer is represented by the drawn structure; absolute descriptors are not printed."
    ))
    for label in ["11a", "11b", "11c", "11d"]:
        bundle["molecules"].append(unresolved_molecule(
            f"mol_{label}", label, [EV[3]],
            f"Allene [2+2] cycloadduct {label} with an exocyclic methylene group, as drawn in Scheme 5.",
            chirality="specified", geometry="not_applicable", representation="not_encoded",
            atom_centers={"source": "wedge/dash-defined fused-ring stereocenters", "R/S": "not printed"},
            assignment="A single regioisomer and the drawn major enantiomer are reported; absolute R/S descriptors are not printed."
        ))
    bundle["molecules"].extend([
        unresolved_molecule(
            "mol_13a", "13a", [EV[4]],
            "6-Chloro 4,5,5a-trihydrocyclobuta-2H-furan product 13a as drawn in Scheme 6.",
            chirality="specified", geometry="not_applicable", representation="not_encoded",
            atom_centers={"source": "wedge/dash-defined fused-ring stereocenters", "R/S": "not printed"},
            assignment="Drawn major enantiomer; absolute descriptors are not printed."
        ),
        unresolved_molecule(
            "mol_16", "16", [EV[4]],
            "Cyclobutane product 16 from substrate 15 as drawn in Scheme 6.",
            chirality="specified", geometry="not_applicable", representation="not_encoded",
            atom_centers={"source": "wedge/dash-defined fused-ring stereocenters", "R/S": "not printed"},
            assignment="Drawn major enantiomer; absolute descriptors are not printed."
        ),
        unresolved_molecule(
            "mol_18", "18", [EV[4]],
            "Spirocyclic diene product 18 from the competing pathway in Scheme 6.",
            chirality="specified", geometry="unspecified", representation="not_encoded",
            atom_centers={"source": "wedge-defined spirocyclic stereocenter", "R/S": "not printed"},
            assignment="Drawn major enantiomer; 92% ee is reported, but absolute descriptors and alkene E/Z labels are not printed."
        ),
        unresolved_molecule(
            "mol_17", "17", [EV[3], EV[4]],
            "Proposed 1,4-diradical intermediate 17 shown in brackets in Scheme 6; not an isolated product.",
            chirality="unresolved", geometry="unresolved", representation="unresolved",
            assignment="Proposed intermediate only; stereochemical drawing is mechanistic, not an isolated-sample assignment."
        ),
    ])

    # Generic structures retain their R/X variables rather than masquerading as molecules.
    bundle["molecules"].extend([
        template_molecule(
            "mol_5_generic", "5", [EV[2]],
            "3-alkyl-4-alkenyloxyquinolin-2(1H)-one scaffold with variables R, R1, and X",
            [
                {"label": "R", "attachment": "quinolone C3", "definition": "3-alkyl substituent", "allowed": ["Me", "Et"]},
                {"label": "R1", "attachment": "tethered alkene", "definition": "alkene substitution/cyclic tether variable", "allowed": "values shown by products 7b-7m"},
                {"label": "X", "attachment": "benzo ring", "definition": "ring substituent", "allowed": ["H", "Me", "Cl", "CN", "OMe", "F"]},
            ], "Scheme 2 uses an R-group drawing; it is not one concrete molecule.",
            chirality="unspecified", geometry="unspecified"
        ),
        template_molecule(
            "mol_7_generic", "7", [EV[2]],
            "Fused cyclobutane products of generic substrates 5 with wedge/dash-defined relative stereochemistry",
            [
                {"label": "R", "attachment": "ring junction/C3-derived center", "definition": "alkyl substituent", "allowed": ["Me", "Et"]},
                {"label": "R1", "attachment": "cyclobutane substituent", "definition": "alkene-derived substituent", "allowed": "values shown for 7b-7m"},
                {"label": "X", "attachment": "benzo ring", "definition": "ring substituent", "allowed": ["H", "Me", "Cl", "CN", "OMe", "F"]},
            ], "Scheme 2 product is a generic R-group drawing and cannot be a concrete SMILES.",
            chirality="specified", geometry="not_applicable"
        ),
        template_molecule(
            "mol_10_generic", "10", [EV[3]],
            "3-alkyl-4-(penta-3,4-dienyloxy)quinolone scaffold with variables R, R1, and X",
            [
                {"label": "R", "attachment": "allene terminal carbon", "definition": "allene substituent", "allowed": ["H", "Me"]},
                {"label": "R1", "attachment": "quinolone C3", "definition": "3-alkyl group", "allowed": ["Me", "Et"]},
                {"label": "X", "attachment": "benzo ring", "definition": "ring substituent", "allowed": ["H", "Cl"]},
            ], "Scheme 5 uses a generic variable drawing; concrete instances are registered separately.",
            chirality="unresolved", geometry="not_applicable"
        ),
        template_molecule(
            "mol_11_generic", "11", [EV[3]],
            "Fused cyclobutane products 11 bearing an exocyclic methylene and variable R/R1/X groups",
            [
                {"label": "R", "attachment": "cyclobutane substituent", "definition": "allene-derived substituent", "allowed": ["H", "Me"]},
                {"label": "R1", "attachment": "ring stereocenter", "definition": "3-alkyl group", "allowed": ["Me", "Et"]},
                {"label": "X", "attachment": "benzo ring", "definition": "ring substituent", "allowed": ["H", "Cl"]},
            ], "Scheme 5 product is a generic variable drawing with non-encoded wedge/dash stereochemistry.",
            chirality="specified", geometry="not_applicable"
        ),
        template_molecule(
            "mol_13_generic", "13", [EV[3], EV[4]],
            "Strained cyclobuta-dihydrofuran intermediate/product family 13 with benzo-ring variable X",
            [{"label": "X", "attachment": "benzo ring", "definition": "ring substituent", "allowed": "not enumerated in the main text"}],
            "Scheme 6 ring-opening panel is generic in X and represents an unstable product family.",
            chirality="specified", geometry="not_applicable"
        ),
        template_molecule(
            "mol_rac14_generic", "rac-14", [EV[3], EV[4]],
            "Racemic furan-containing azocinone ring-opening product family 14 with variable X",
            [{"label": "X", "attachment": "benzo ring", "definition": "ring substituent", "allowed": "X=H product reported isolated; other values not specified"}],
            "Generic ring-opening product; the main text states products 14 are racemic.",
            chirality="racemic", geometry="unspecified"
        ),
    ])

    # Prior literature structures in Figure 1 are grouped as a single literature context record.
    for label, stereo, note in [
        ("rac1", "racemic", "Prior-literature cyclobutane rac-1; source product structure only."),
        ("rac2", "racemic", "Prior-literature cyclobutane rac-2; source product structure only."),
        ("rac3", "racemic", "Prior-literature cyclobutane rac-3; source product structure only."),
        ("4", "specified", "Prior-literature enantioenriched dihydroquinolone 4; wedge/dash structure shown."),
    ]:
        bundle["molecules"].append(unresolved_molecule(
            f"mol_lit_{label}", f"rac-{label[-1]}" if label.startswith("rac") else label, [EV[1]], note,
            chirality=stereo, geometry="not_applicable", representation="not_encoded",
            atom_centers={"source": "Figure 1 stereochemical drawing", "absolute_descriptors": "not fully printed"},
            mixture="racemic" if stereo == "racemic" else None,
            assignment="Prior-literature assignment; not a current-paper experiment."
        ))

    # Literature context record, no current-paper reaction inferred from the product montage.
    lit = base_reaction(
        "figure1_prior_literature", [EV[1]], record_type="literature_example",
        location={"page": 1, "figure": "Figure 1"},
        reaction_name="Prior literature examples of α-alkyl-substituted enone photocycloadditions",
        transformation="literature context; intramolecular [2+2] photocycloaddition products"
    )
    lit["participants"] = [
        participant(f"p_lit_{label}", f"mol_lit_{label}", "product", [EV[1]], product_role="literature product")
        for label in ["rac1", "rac2", "rac3", "4"]
    ]
    lit["outcomes"] = yield_and_ee(
        "lit4", "mol_lit_4", 99, 92, [EV[1]],
        ee_note="prior-literature enantioenriched product 4; current paper does not establish the analytical method"
    )
    lit["notes"] = [f(
        "Figure 1 shows rac-1, rac-2, rac-3, and enantioenriched 4 from cited prior work; compound 4 followed dynamic salt formation, photocycloaddition, and base treatment.",
        "Figure 1 is prior literature context, not a new experiment in this paper.", [EV[1]]
    )]
    set_unresolved_classification(lit, [EV[1]], "Only product structures and narrative context are shown; reactants and atom mapping are not supplied.")
    lit["review"]["unresolved"] = ["Prior-work experimental quantities, conditions, and measurement methods are not reported on page 1."]
    bundle["reactions"].append(lit)

    # Scheme 1: control 5a and ambient-temperature 5b experiment.
    control = base_reaction(
        "scheme1_5a_no_reaction", [EV[2], TXT[1]], categories=["control"],
        location={"page": 2, "scheme": "Scheme 1", "row": "5a, R=H"},
        reaction_name="Sensitized intramolecular [2+2] photocycloaddition",
        transformation="no observed reaction under attempted photocycloaddition conditions"
    )
    add_common_photoreaction(
        control, "mol_5a", [], [EV[2]], 10, ["mol_phcf3"], "PhCF3",
        ("25 °C", q(25, "°C", "25 °C")),
        ("t = 1 h", q(1, "h", "t = 1 h"), "original_extraction"),
        light_source="fluorescent lamps / hν"
    )
    no_rxn = st.new_outcome("no_reaction", "no_reaction", ["mol_5a"], [EV[2]])
    no_rxn["value"] = f("— (no reaction)", "no reaction", [EV[2]])
    no_rxn["scope"] = f("substrate 5a", {"target": "mol_5a"}, [EV[2]])
    control["outcomes"] = [no_rxn]
    set_unresolved_classification(control, [EV[2]], "Negative control: no observed product, so intended bond formation is not classified.", not_applicable=True)
    control["review"]["unresolved"] = ["Conversion and no-reaction detection method are not reported; no product is invented.", "Scale, atmosphere, workup, and purification are unavailable without SI."]
    bundle["reactions"].append(control)

    ambient = base_reaction(
        "scheme1_5b_ambient", [EV[2]], categories=["optimization"],
        location={"page": 2, "scheme": "Scheme 1", "row": "5b, R=Me"},
        reaction_name="Sensitized intramolecular [2+2] photocycloaddition",
        transformation="enantioselective intramolecular [2+2] photocycloaddition"
    )
    add_common_photoreaction(
        ambient, "mol_5b", ["mol_7b"], [EV[2]], 10, ["mol_phcf3"], "trifluorotoluene (PhCF3)",
        ("25 °C", q(25, "°C", "25 °C")),
        ("t = 1 h", q(1, "h", "t = 1 h"), "original_extraction"),
        light_source="fluorescent lamps / hν"
    )
    ambient["outcomes"] = yield_and_ee("ambient_7b", "mol_7b", 98, 72, [EV[2]])
    set_partial_cycloaddition_classification(ambient, [EV[2]])
    ambient["review"]["unresolved"] = ["Yield and ee measurement methods are unreported in the supplied main text.", "Scale, concentration, atmosphere, workup, and purification are unavailable without SI."]
    bundle["reactions"].append(ambient)

    # Scheme 2 generic reaction and twelve separate scope entries.
    scheme2 = base_reaction(
        "scheme2_general", [EV[2]], record_type="generic_scheme",
        location={"page": 2, "scheme": "Scheme 2", "region": "generic arrow"},
        reaction_name="Intramolecular enantioselective [2+2] photocycloaddition of 3-alkyl-4-alkenyloxyquinolones 5 to cyclobutanes 7",
        transformation="generic enantioselective intramolecular [2+2] photocycloaddition"
    )
    scheme2["participants"] = [
        participant("p_5_generic", "mol_5_generic", "reactant", [EV[2]]),
        participant("p_7_generic", "mol_7_generic", "product", [EV[2]], product_role="generic major product"),
        participant("p_catalyst6", "mol_catalyst_6", "catalyst", [EV[2]], mol_percent=10, reference="mol_5_generic"),
        participant("p_phcf3", "mol_phcf3", "solvent", [EV[2]]),
    ]
    scheme2["steps"] = [reaction_step(
        "irradiation", [p["participant_id"] for p in scheme2["participants"]], [EV[2]],
        temperature=("−25 °C", q(-25, "°C", "−25 °C")),
        time=("t = 1 h", q(1, "h", "t = 1 h"), "original_extraction")
    )]
    solvent_fields(scheme2, "PhCF3", [EV[2]])
    light_fields(scheme2, [EV[2]])
    set_partial_cycloaddition_classification(scheme2, [EV[2]])
    scheme2["review"]["unresolved"] = ["Generic R/X variables are not a concrete reaction instance."]
    bundle["reactions"].append(scheme2)

    scope2 = [
        ("b", 97, 88), ("c", 99, 99), ("d", 98, 96), ("e", 99, 98),
        ("f", 93, 96), ("g", 99, 93), ("h", 94, 98), ("i", 98, 94),
        ("j", 99, 93), ("k", 94, 96), ("l", 99, 98), ("m", 97, 93),
    ]
    for letter, yld, ee in scope2:
        evidence = [EV[2]]
        r = base_reaction(
            f"scheme2_5{letter}_to_7{letter}", evidence, categories=["substrate_scope"],
            location={"page": 2, "scheme": "Scheme 2", "entry": f"7{letter}"},
            reaction_name="Intramolecular enantioselective [2+2] photocycloaddition",
            transformation="enantioselective intramolecular [2+2] photocycloaddition"
        )
        time_value = 1.5 if letter == "k" else 1
        time_raw = "t = 1.5 h (footnote a)" if letter == "k" else "t = 1 h"
        add_common_photoreaction(
            r, f"mol_5{letter}", [f"mol_7{letter}"], evidence, 10,
            ["mol_phcf3"], "PhCF3",
            ("−25 °C", q(-25, "°C", "−25 °C")),
            (time_raw, q(time_value, "h", time_raw), "inherited_general_condition")
        )
        if letter == "k":
            solvent_fields(r, "PhCF3", evidence, concentration=("c = 1.5 mM", q(1.5, "mM", "c = 1.5 mM")))
            r["review"]["unresolved"].append("The 1.5 mM footnote does not explicitly name the concentration basis.")
        r["general_scheme_ids"] = [f"{DOC}:scheme2_general"]
        ee_note = (
            "absolute configuration determined by anomalous single-crystal X-ray diffraction"
            if letter == "j" else
            "drawn major enantiomer; absolute configuration assigned by analogy to 7j"
        )
        r["outcomes"] = yield_and_ee(f"7{letter}", f"mol_7{letter}", yld, ee, evidence, ee_note=ee_note)
        set_partial_cycloaddition_classification(r, evidence)
        r["review"]["unresolved"].extend([
            "Yield and ee measurement methods are unreported in the supplied main text.",
            "Exact scale, atmosphere, workup, and purification are unavailable without SI.",
        ])
        if letter not in ("b",):
            r["review"]["unresolved"].append("Compact substrate/product variant is retained by label and source drawing; no unsafe product SMILES was invented.")
        bundle["reactions"].append(r)

    # Scheme 3 catalyst-loading series: three independent experiments.
    for loading, yld, ee, tag in [(10, 99, 99, "10mol"), (1.0, 99, 95, "1mol"), (0.5, 99, 93, "0p5mol")]:
        r = base_reaction(
            f"scheme3_5n_{tag}", [EV[3]], categories=["optimization"],
            location={"page": 3, "scheme": "Scheme 3", "catalyst_loading": f"{loading:g} mol%"},
            reaction_name="Enantioselective [2+2] photocycloaddition of 5n to 7n",
            transformation="enantioselective intramolecular [2+2] photocycloaddition"
        )
        add_common_photoreaction(
            r, "mol_5n", ["mol_7n"], [EV[3]], loading, ["mol_phcf3"], "PhCF3",
            ("−25 °C", q(-25, "°C", "−25 °C")), None
        )
        r["general_scheme_ids"] = [f"{DOC}:scheme2_general"]
        r["outcomes"] = yield_and_ee(f"7n_{tag}", "mol_7n", yld, ee, [EV[3]])
        set_partial_cycloaddition_classification(r, [EV[3]])
        r["review"]["unresolved"] = [
            "Irradiation time is not printed in Scheme 3 and is not silently inherited from another scheme.",
            "Yield and ee measurement methods are unreported in the supplied main text.",
            "Scale, concentration, atmosphere, workup, and purification are unavailable without SI.",
            "Exact substrate/product SMILES and absolute descriptors remain unresolved from the compact drawing.",
        ]
        bundle["reactions"].append(r)

    # Scheme 4 fluorinated alkene.
    r4 = base_reaction(
        "scheme4_8_to_9", [EV[3]], categories=["substrate_scope"],
        location={"page": 3, "scheme": "Scheme 4"},
        reaction_name="Intramolecular [2+2] photocycloaddition of 4-(3,4,4-trifluorobut-3-enyloxy)-3-methylquinolone",
        transformation="enantioselective intramolecular [2+2] photocycloaddition"
    )
    add_common_photoreaction(
        r4, "mol_8", ["mol_9"], [EV[3]], 10, ["mol_hfx", "mol_phcf3"], "HFX/PhCF3, 2/1",
        ("−65 °C", q(-65, "°C", "−65 °C")),
        ("t = 2 h", q(2, "h", "t = 2 h"), "original_extraction")
    )
    r4["solvent_conditions"]["mixture_ratio"] = f("HFX/PhCF3, 2/1", {"components": ["HFX", "PhCF3"], "ratio": [2, 1]}, [EV[3]])
    r4["outcomes"] = yield_and_ee("product9", "mol_9", 72, 81, [EV[3]])
    set_partial_cycloaddition_classification(r4, [EV[3]])
    r4["review"]["unresolved"] = ["Yield and ee measurement methods are unreported.", "Scale, atmosphere, workup, and purification are unavailable without SI.", "Product stereochemical SMILES and absolute descriptors remain unresolved."]
    bundle["reactions"].append(r4)

    # Scheme 5 generic and four allene-scope experiments.
    scheme5 = base_reaction(
        "scheme5_general", [EV[3]], record_type="generic_scheme",
        location={"page": 3, "scheme": "Scheme 5", "region": "generic arrow"},
        reaction_name="Enantioselective intramolecular [2+2] photocycloaddition of 3-alkyl-4-(penta-3,4-dienyloxy)quinolones 10",
        transformation="generic enantioselective intramolecular allene [2+2] photocycloaddition"
    )
    scheme5["participants"] = [
        participant("p_10_generic", "mol_10_generic", "reactant", [EV[3]]),
        participant("p_11_generic", "mol_11_generic", "product", [EV[3]], product_role="generic single regioisomer"),
        participant("p_catalyst6", "mol_catalyst_6", "catalyst", [EV[3]], mol_percent=10, reference="mol_10_generic"),
        participant("p_phcf3", "mol_phcf3", "solvent", [EV[3]]),
    ]
    scheme5["steps"] = [reaction_step(
        "irradiation", [p["participant_id"] for p in scheme5["participants"]], [EV[3]],
        temperature=("−25 °C", q(-25, "°C", "−25 °C")),
        time=("t = 1.5–2.5 h", range_q(1.5, 2.5, "h", "t = 1.5–2.5 h"), "original_extraction")
    )]
    solvent_fields(scheme5, "PhCF3", [EV[3]])
    light_fields(scheme5, [EV[3]])
    set_partial_cycloaddition_classification(scheme5, [EV[3]])
    scheme5["review"]["unresolved"] = ["The 1.5–2.5 h range is not assigned to individual entries in the main-text scheme."]
    bundle["reactions"].append(scheme5)

    scope5 = [("a", 85, 91), ("b", 98, 96), ("c", 97, 91), ("d", 90, 92)]
    for letter, yld, ee in scope5:
        r = base_reaction(
            f"scheme5_10{letter}_to_11{letter}", [EV[3]], categories=["substrate_scope"],
            location={"page": 3, "scheme": "Scheme 5", "entry": f"10{letter}/11{letter}"},
            reaction_name="Enantioselective intramolecular [2+2] photocycloaddition of an allene-tethered quinolone",
            transformation="enantioselective intramolecular allene [2+2] photocycloaddition"
        )
        add_common_photoreaction(
            r, f"mol_10{letter}", [f"mol_11{letter}"], [EV[3]], 10, ["mol_phcf3"], "PhCF3",
            ("−25 °C", q(-25, "°C", "−25 °C")),
            ("t = 1.5–2.5 h", range_q(1.5, 2.5, "h", "t = 1.5–2.5 h"), "inherited_general_condition")
        )
        r["general_scheme_ids"] = [f"{DOC}:scheme5_general"]
        r["outcomes"] = yield_and_ee(f"11{letter}", f"mol_11{letter}", yld, ee, [EV[3]])
        set_partial_cycloaddition_classification(r, [EV[3]])
        r["notes"] = [f("a single regioisomer was obtained with an exocyclic double bond", "single regioisomer with exocyclic methylene", [EV[3]])]
        r["review"]["unresolved"] = [
            "Individual irradiation time is not assigned; only the Scheme 5 range 1.5–2.5 h is retained.",
            "Yield and ee measurement methods are unreported.",
            "Scale, concentration, atmosphere, workup, and purification are unavailable without SI.",
            "Product absolute descriptors and stereochemical SMILES remain unresolved.",
        ]
        bundle["reactions"].append(r)

    # Scheme 6 upper [2+2] experiment.
    r12 = base_reaction(
        "scheme6_12a_to_13a", [EV[4]], categories=["substrate_scope"],
        location={"page": 4, "scheme": "Scheme 6", "arrow": "upper"},
        reaction_name="Enantioselective [2+2] photocycloaddition of 4-(buta-2,3-dienyloxy)-3-methylquinolone 12a",
        transformation="enantioselective intramolecular allene [2+2] photocycloaddition"
    )
    add_common_photoreaction(
        r12, "mol_12a", ["mol_13a"], [EV[4]], 10, ["mol_phcf3"], "PhCF3",
        ("−25 °C", q(-25, "°C", "−25 °C")),
        ("t = 1.5 h", q(1.5, "h", "t = 1.5 h"), "original_extraction")
    )
    r12["outcomes"] = yield_and_ee("13a", "mol_13a", 85, 88, [EV[4]])
    set_partial_cycloaddition_classification(r12, [EV[4]])
    r12["review"]["unresolved"] = ["Yield and ee measurement methods are unreported.", "Scale, concentration, atmosphere, workup, and purification are unavailable without SI.", "Product stereochemical SMILES and absolute descriptors remain unresolved."]
    bundle["reactions"].append(r12)

    # Generic acid-promoted ring opening; no specific yield/time is printed.
    opening = base_reaction(
        "scheme6_13_to_rac14_ring_opening", [EV[3], EV[4]], record_type="generic_scheme",
        location={"pages": [3, 4], "scheme": "Scheme 6", "arrow": "middle ring opening / tautomerization"},
        reaction_name="Acid-promoted ring opening and tautomerization of cyclobuta-dihydrofuran products 13",
        transformation="cyclobutane fragmentation/ring opening and furan formation"
    )
    opening["participants"] = [
        participant("p_13", "mol_13_generic", "reactant", [EV[4]]),
        participant("p_rac14", "mol_rac14_generic", "product", [EV[4]], product_role="racemic ring-opening product family"),
        participant("p_h2so4", "mol_h2so4", "reagent", [EV[3], EV[4]], mol_percent=0.5, reference="mol_13_generic"),
    ]
    opening["steps"] = [reaction_step(
        "ring_opening", [p["participant_id"] for p in opening["participants"]], [EV[3], EV[4]],
        operation="ring opening followed by tautomerization"
    )]
    racemic = product_outcome(
        "rac14_racemic", "selectivity", "mol_rac14_generic", "products 14 were racemic", "racemic",
        [EV[3], EV[4]], role="ring-opening product family",
        scope={"metric": "enantiomeric composition", "assignment": "racemic; ratio not numerically reported"}
    )
    opening["outcomes"] = [racemic]
    opening["notes"] = [
        f("Traces of acid (e.g., 0.5 mol % of H2SO4) led to ring opening", "0.5 mol% H2SO4 is an example acid loading stated in the prose", [EV[3]]),
        f("One of the fragmentation products (X = H) was isolated and fully characterized (see the SI)", "X=H product isolated; yield and procedure absent from supplied main text", [EV[3]]),
    ]
    set_unresolved_classification(opening, [EV[3], EV[4]], "Generic X-substituted structures and incomplete atom correspondence prevent a defensible complete bond-change classification.")
    opening["review"]["unresolved"] = ["The exact experiment identity, acid exposure time/temperature, yield, and workup are not reported in the supplied main text.", "The SI is mentioned but was not supplied or associated.", "The X=H isolated product is not assigned a numeric yield in the main text."]
    bundle["reactions"].append(opening)

    # Scheme 6 lower competing-pathway experiment with two products.
    competing = base_reaction(
        "scheme6_15_to_16_and_18", [EV[3], EV[4]], categories=["substrate_scope"],
        location={"pages": [3, 4], "scheme": "Scheme 6", "arrow": "lower"},
        reaction_name="Competing enantioselective photoreactions of terminally dimethyl-substituted allene 15",
        transformation="intramolecular [2+2] cycloaddition competing with C-C bond formation/hydrogen abstraction"
    )
    add_common_photoreaction(
        competing, "mol_15", ["mol_16", "mol_18"], [EV[4]], 10, ["mol_phcf3"], "PhCF3",
        ("−25 °C", q(-25, "°C", "−25 °C")),
        ("t = 1 h", q(1, "h", "t = 1 h"), "original_extraction")
    )
    competing["outcomes"] = (
        yield_and_ee("16", "mol_16", 50, 98, [EV[4]], role="cyclobutane product 16") +
        yield_and_ee("18", "mol_18", 47, 92, [EV[3], EV[4]], yield_basis="isolated", role="spirocyclic diene product 18")
    )
    competing["notes"] = [f(
        "It seems likely that [18] stems from the same 1,4-diradical precursor 17 as cyclobutane 16",
        "Authors propose a shared 1,4-diradical 17 and intramolecular hydrogen abstraction for the 18 branch.", [EV[3], EV[4]]
    )]
    # Three branch-specific new C-C connections are represented only as partial, non-mapped changes:
    # two for 16 and one for 18. No new DB label is asserted without a full atom correspondence.
    set_partial_cycloaddition_classification(
        competing, [EV[3], EV[4]], formed_bonds=3,
        note="The three abstract C-C changes cover two bonds in product 16 and one in product 18; branch-specific atom identities and other bond-order/H changes remain unresolved."
    )
    competing["review"]["unresolved"] = [
        "The 50% yield basis for 16 and both ee measurement methods are unreported; only the 47% yield of 18 is explicitly described as isolated.",
        "Branch-specific atom mapping, hydrogen transfer, and diene bond-order changes are unresolved; only unambiguous C-C formation is classified.",
        "Scale, concentration, atmosphere, workup, and purification are unavailable without SI.",
        "Product stereochemical SMILES and absolute descriptors remain unresolved.",
    ]
    bundle["reactions"].append(competing)

    mechanism = base_reaction(
        "scheme6_proposed_diradical_17", [EV[3], EV[4]], record_type="proposed_mechanism",
        location={"pages": [3, 4], "scheme": "Scheme 6", "region": "intermediate 17 and curved arrows"},
        reaction_name="Proposed common 1,4-diradical pathway to 16 and 18",
        transformation="author-proposed energy-transfer, C-C bond formation, ISC/ring closure or hydrogen abstraction"
    )
    mechanism["participants"] = [
        participant("p_15", "mol_15", "reactant", [EV[3], EV[4]]),
        participant("p_17", "mol_17", "unknown", [EV[3], EV[4]]),
        participant("p_16", "mol_16", "product", [EV[4]], product_role="proposed ring-closure branch product"),
        participant("p_18", "mol_18", "product", [EV[4]], product_role="proposed hydrogen-abstraction branch product"),
    ]
    mechanism["notes"] = [
        f("The quinolone triplet undergoes C-C bond formation ... to a 1,4-diradical, which after ISC generates the product", "author-proposed general photochemical pathway", [EV[3]]),
        f("The additional methyl groups open an alternative reaction pathway to intermediate 17, allowing for intramolecular hydrogen abstraction", "author-proposed competing pathway to 18", [EV[3], EV[4]]),
    ]
    set_unresolved_classification(mechanism, [EV[3], EV[4]], "Proposed intermediate and elementary steps are not treated as an experimentally observed standalone transformation.", not_applicable=True)
    mechanism["review"]["unresolved"] = ["Intermediate 17 is proposed, not isolated or directly observed; no kinetic, trapping, or isotopic experiment is reported on pages 1–5."]
    bundle["reactions"].append(mechanism)

    bundle["source_context"] = [
        {
            "context_id": f"{DOC}:publication_metadata",
            "document_id": DOC,
            "kind": "publication_metadata",
            "title": TITLE,
            "authors": ["Xinyao Li", "Christian Jandl", "Thorsten Bach"],
            "journal": "Organic Letters",
            "year": 2020,
            "volume": 22,
            "pages": "3618-3622",
            "doi": "10.1021/acs.orglett.0c01065",
            "evidence_ids": [EV[1], TXT[1]],
        },
        {
            "context_id": f"{DOC}:figure2_photophysics",
            "document_id": DOC,
            "kind": "non_synthetic_observation",
            "location": {"page": 2, "figure": "Figure 2"},
            "finding": "Fluorescence/phosphorescence spectra and triplet energies: 5a ET=298 kJ/mol, 5b ET=273 kJ/mol at 77 K in ethanol; not a synthesis experiment.",
            "evidence_ids": [EV[2], TXT[2]],
        },
        {
            "context_id": f"{DOC}:figure3_xray_binding",
            "document_id": DOC,
            "kind": "crystallography_and_mechanistic_context",
            "location": {"page": 2, "figure": "Figure 3"},
            "finding": "Absolute configuration of 7j determined by anomalous X-ray diffraction; other product assignments are by analogy. Putative 1:1 complex 6·5j is proposed.",
            "evidence_ids": [EV[2], TXT[2]],
        },
        {
            "context_id": f"{DOC}:si_notice",
            "document_id": DOC,
            "kind": "source_availability",
            "location": {"page": 4, "section": "Associated Content"},
            "finding": "The article states that SI contains synthetic procedures and characterization, but no SI file was explicitly provided in this workspace; no association was inferred.",
            "evidence_ids": [EV[4], TXT[4]],
        },
        {
            "context_id": f"{DOC}:ccdc",
            "document_id": DOC,
            "kind": "crystallographic_accession",
            "finding": "CCDC 1988524 is printed for supplementary crystallographic data.",
            "evidence_ids": [EV[4], TXT[4]],
        },
        {
            "context_id": f"{DOC}:page5_bibliography",
            "document_id": DOC,
            "kind": "bibliography_only",
            "finding": "Page 5 contains references only; no current-paper reaction arrow, table, or characterization dataset is present.",
            "evidence_ids": [EV[5], TXT[5]],
        },
    ]

    # Coverage contract and raw page decisions are deliberately separate from the formal bundle.
    page_records = {
        1: [f"{DOC}:figure1_prior_literature"],
        2: [r["reaction_id"] for r in bundle["reactions"] if "scheme1_" in r["reaction_id"] or "scheme2_" in r["reaction_id"]],
        3: [r["reaction_id"] for r in bundle["reactions"] if any(x in r["reaction_id"] for x in ("scheme3_", "scheme4_", "scheme5_", "scheme6_13", "scheme6_15", "scheme6_proposed"))],
        4: [r["reaction_id"] for r in bundle["reactions"] if r["reaction_id"].startswith(f"{DOC}:scheme6_")],
        5: [],
    }
    coverage = {
        "document_id": DOC,
        "semantic_accuracy_verified": False,
        "pages": [
            {
                "page": 1,
                "state": "attempted",
                "regions": [
                    {"region": "title, abstract, and graphical abstract", "finding": "Current-paper overview and generic visible-light [2+2] transformation; summary values are not treated as standalone experiments.", "evidence_ids": [EV[1], TXT[1]]},
                    {"region": "Figure 1 and adjacent prose", "finding": "rac-1, rac-2, rac-3, and enantioenriched 4 are cited prior-literature products, recorded as literature context rather than current experiments.", "evidence_ids": [EV[1], TXT[1]]},
                    {"region": "right-column discussion ending at Scheme 1 reference", "finding": "Current study framing and the beginning of the 5a failed-attempt discussion; experimental drawing continues on page 2 and is recorded there.", "evidence_ids": [EV[1], TXT[1]]},
                ],
                "record_ids": page_records[1],
                "unresolved": ["Graphical-abstract structures are generic and do not provide row-specific conditions.", "Prior-literature Figure 1 procedures and measurement methods are not reproduced."],
            },
            {
                "page": 2,
                "state": "attempted",
                "regions": [
                    {"region": "Scheme 1", "finding": "Two separate experiments: 5a no reaction and 5b to 7b at 25 °C (98% yield, 72% ee).", "evidence_ids": [EV[2], TXT[2]]},
                    {"region": "Figure 2 and photophysical prose", "finding": "Spectroscopy/triplet-energy observations, preserved as source_context rather than synthesis records.", "evidence_ids": [EV[2], TXT[2]]},
                    {"region": "Scheme 2 generic arrow, twelve product panels 7b-7m, and footnote a", "finding": "One generic scheme plus twelve distinct substrate-scope experiments; 7k alone uses c=1.5 mM and 1.5 h.", "evidence_ids": [EV[2], TXT[2]]},
                    {"region": "Figure 3", "finding": "7j absolute configuration by anomalous X-ray diffraction and proposed 6·5j binding model; crystallography and mechanistic context only.", "evidence_ids": [EV[2], TXT[2]]},
                ],
                "record_ids": page_records[2],
                "unresolved": ["Scheme 2 yield/ee analytical methods, reaction scale, and workup are absent.", "Dense product stereochemical SMILES are not safely recoverable from the compact main-text drawings."],
            },
            {
                "page": 3,
                "state": "attempted",
                "regions": [
                    {"region": "mechanistic prose at upper left", "finding": "Author-proposed triplet energy transfer, 1,4-diradical formation, and ISC; not an observed mechanistic experiment.", "evidence_ids": [EV[3], TXT[3]]},
                    {"region": "Scheme 3", "finding": "Three catalyst-loading experiments (10, 1.0, and 0.5 mol% 6), all 99% yield with 99, 95, and 93% ee, respectively.", "evidence_ids": [EV[3], TXT[3]]},
                    {"region": "Scheme 4", "finding": "8 to 9 under 420 nm light, 10 mol% 6, −65 °C, 2 h, HFX/PhCF3 2/1; 72% yield and 81% ee.", "evidence_ids": [EV[3], TXT[3]]},
                    {"region": "Scheme 5", "finding": "One generic allene [2+2] scheme and four distinct scope entries 10a-d to 11a-d with yields/ee; only a common 1.5–2.5 h range is printed.", "evidence_ids": [EV[3], TXT[3]]},
                    {"region": "lower-right Scheme 6 discussion", "finding": "Narrative evidence for 12/13 regioselectivity, acid-promoted 13 to rac-14 ring opening, and competing products 16/18 via proposed 17; drawings continue on page 4.", "evidence_ids": [EV[3], TXT[3]]},
                ],
                "record_ids": page_records[3],
                "unresolved": ["Scheme 3 irradiation times are not printed.", "Scheme 5 row-specific times are not assigned.", "Ring-opening experiment details and yield are deferred to unavailable SI."],
            },
            {
                "page": 4,
                "state": "attempted",
                "regions": [
                    {"region": "Scheme 6 upper arrow", "finding": "Specific 12a to 13a photocycloaddition: 420 nm, 10 mol% 6, 1.5 h, −25 °C, PhCF3, 85% yield, 88% ee.", "evidence_ids": [EV[4], TXT[4]]},
                    {"region": "Scheme 6 middle arrow", "finding": "Generic 13 to rac-14 ring opening/tautomerization; no numeric yield or time/temperature printed.", "evidence_ids": [EV[4], TXT[4]]},
                    {"region": "Scheme 6 lower arrow and intermediate 17", "finding": "One experiment on 15 gives 16 (50%, 98% ee) and 18 (47%, 92% ee); 17 is a proposed common diradical, not an isolated experiment.", "evidence_ids": [EV[4], TXT[4]]},
                    {"region": "Associated Content and CCDC notice", "finding": "SI availability and CCDC 1988524 are source context; no SI file was supplied or associated.", "evidence_ids": [EV[4], TXT[4]]},
                    {"region": "references beginning in right column", "finding": "Bibliography, no additional current-paper experiment.", "evidence_ids": [EV[4], TXT[4]]},
                ],
                "record_ids": page_records[4],
                "unresolved": ["Exact ring-opening procedure and X=H product yield are unavailable without SI.", "The 50% yield basis for 16 and ee analytical methods are not reported.", "Product absolute descriptors are not printed for 13a, 16, or 18."],
            },
            {
                "page": 5,
                "state": "attempted",
                "regions": [
                    {"region": "both columns", "finding": "References 3–19 only; no current-paper synthesis experiment, table, reaction arrow, or standalone spectral peak table.", "evidence_ids": [EV[5], TXT[5]]},
                ],
                "record_ids": page_records[5],
                "unresolved": [],
            },
        ],
    }
    (RUN / "coverage/pages_0001_0005.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")
    raw = {
        "document_id": DOC,
        "title": TITLE,
        "pages_attempted": [1, 2, 3, 4, 5],
        "specific_experiments": [r["reaction_id"] for r in bundle["reactions"] if r["record_type"] == "specific_experiment"],
        "non_specific_records": [r["reaction_id"] for r in bundle["reactions"] if r["record_type"] != "specific_experiment"],
        "decisions": [
            "Scheme 1 rows are separate because one is a no-reaction control and the other a successful ambient-temperature experiment.",
            "Each Scheme 2 product panel and each Scheme 3 catalyst loading is a separate experiment.",
            "Scheme 5 entries retain the printed 1.5–2.5 h range because row-specific times are absent.",
            "The two products 16 and 18 belong to one experiment on 15, not two experiments.",
            "Proposed intermediate 17 is a proposed_mechanism record, not mechanistic proof.",
            "Figure 2 spectroscopy and page 5 bibliography are source_context, not synthesis records.",
            "No SI was associated; missing procedures and analytical methods remain missing.",
            "Dense product structures retain source labels and stereochemical descriptions without unsafe SMILES.",
        ],
        "chemical_accuracy_verified": False,
    }
    (RUN / "raw_model/pages_0001_0005_extraction.json").write_text(json.dumps(raw, ensure_ascii=False, indent=2) + "\n")
    (RUN / "raw_model/stereo_checks.json").write_text(json.dumps(STEREO_CHECKS, ensure_ascii=False, indent=2) + "\n")
    (RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")

    tasks_path = RUN / "tasks.json"
    tasks = json.loads(tasks_path.read_text())
    for task in tasks:
        task["state"] = "attempted"
        task["reason"] = "Page image and extracted text were inspected; detailed region findings are in coverage/pages_0001_0005.json."
    tasks_path.write_text(json.dumps(tasks, ensure_ascii=False, indent=2) + "\n")

    print(json.dumps({
        "molecules": len(bundle["molecules"]),
        "reactions": len(bundle["reactions"]),
        "specific_experiments": len(raw["specific_experiments"]),
        "non_specific_records": len(raw["non_specific_records"]),
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
