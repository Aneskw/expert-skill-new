#!/usr/bin/env python3
"""Build the evidence-linked v3.1 extraction for Chem. Eur. J. 2021, 27, 12335–12340.

This is a source-specific construction script. It writes generated records only; it
does not modify the extraction tools, schemas, or quality rules.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_c7aae70d07836214"
RUN = ROOT / "runs/current"


def ev(page: int) -> list[str]:
    return [f"{DOC}:page{page}:text", f"{DOC}:page{page}:image"]


def evs(*pages: int) -> list[str]:
    return [item for page in pages for item in ev(page)]


def extracted(raw, normalized, evidence):
    return st.fact(raw, normalized, evidence, method="original_extraction")


def derived(raw, normalized, evidence, why):
    return st.fact(raw, normalized, evidence, method="structure_derivation", derivation=why)


def inferred(raw, normalized, evidence, why):
    return st.fact(raw, normalized, evidence, method="inference", derivation=why)


def inherited(raw, normalized, evidence, why):
    return st.fact(raw, normalized, evidence, method="inherited_general_condition", derivation=why)


def missing(evidence=()):
    return st.fact(evidence=evidence)


def put_fact(container, key, raw, normalized, evidence, method="original_extraction", why=None):
    if method == "original_extraction":
        container[key] = extracted(raw, normalized, evidence)
    elif method == "structure_derivation":
        container[key] = derived(raw, normalized, evidence, why)
    elif method == "inference":
        container[key] = inferred(raw, normalized, evidence, why)
    elif method == "inherited_general_condition":
        container[key] = inherited(raw, normalized, evidence, why)
    else:
        raise ValueError(method)


def concrete_molecule(mid, label, name, smiles, evidence, *, note=None,
                      chirality="achiral", geometry="not_applicable"):
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["representation_type"] = "smiles"
    m["occurrences"][0]["compound_label"] = extracted(label, label, evidence)
    m["names"]["original"] = extracted(name, name, evidence)
    m["names"]["english"] = extracted(name, name, evidence)
    m["names"]["abbreviations"] = extracted(label, [label], evidence)
    m["smiles"] = derived(
        {"source": "structure/name in cited page", "label": label}, smiles, evidence,
        "Candidate SMILES transcribed from the source structure/name; parsing is not identity verification."
    )
    diag = st.stereo_inspect(smiles)
    m["stereochemistry"]["atom_centers"] = derived(
        {"source": "source structure plus graph inspection"}, diag["atom_centers"], evidence,
        "RDKit potential/assigned atom-center inspection of the candidate source-derived graph."
    )
    m["stereochemistry"]["double_bonds"] = derived(
        {"source": "source structure plus graph inspection"}, diag["double_bonds"], evidence,
        "RDKit double-bond stereo inspection of the candidate source-derived graph."
    )
    m["stereochemistry"]["axial"] = derived(
        {"source": "source allene/reagent structure"}, [], evidence,
        "No assignable axial stereogenic element in this candidate graph; monosubstituted allenes have two identical H substituents at one terminus."
    )
    m["stereochemistry"]["mixture"] = extracted("no stereochemical mixture specified", None, evidence)
    m["stereochemistry"]["assignment"] = derived(
        {"source": "source drawing/name"}, "no stereochemical assignment required for candidate graph", evidence,
        "Graph inspection found no specified stereogenic element for this molecule."
    )
    m["stereochemistry"]["chirality_status"] = derived(
        {"source": "source drawing/name"}, chirality, evidence,
        "Status follows the complete candidate graph, not absence of stereochemical marks alone."
    )
    m["stereochemistry"]["geometry_status"] = derived(
        {"source": "source drawing/name"}, geometry, evidence,
        "No independently assignable E/Z geometry applies to the source graph unless stated otherwise."
    )
    m["stereochemistry"]["representation_status"] = derived(
        {"source": "candidate graph"}, "not_applicable", evidence,
        "Candidate graph has no stereochemical element requiring an isomeric tag."
    )
    m["stereochemistry"]["isomeric_smiles"] = derived(
        {"source": "candidate graph"}, diag["canonical_isomeric_smiles"], evidence,
        "Canonical isomeric SMILES calculated from the same candidate graph."
    )
    m["stereochemistry"]["identity_key"] = derived(
        {"source": "candidate graph"}, diag["stereo_graph_key"], evidence,
        "Graph key calculated from canonical isomeric SMILES; it is not a sample identifier."
    )
    m["identity_note"] = inferred(
        note or "Source-derived identity; no CAS number or external-database identity was supplied.",
        note or "Source-derived identity; no CAS number or external-database identity was supplied.",
        evidence, "Identity note records the limits of the source-derived candidate."
    )
    return m


def template_molecule(mid, label, name, evidence, fixed, *, reason,
                      chirality="unresolved", geometry="unresolved",
                      representation="partial", assignment=None, double_bonds=None,
                      mixture=None, variables=None):
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["representation_type"] = "generic_template"
    m["occurrences"][0]["compound_label"] = extracted(label, label, evidence)
    m["names"]["original"] = extracted(name, name, evidence)
    m["names"]["english"] = extracted(name, name, evidence)
    m["names"]["abbreviations"] = extracted(label, [label], evidence)
    m["generic_template"] = {
        "fixed_structure": extracted({"source_drawing": fixed}, fixed, evidence),
        "variables": variables or [],
        "instance_mapping": [],
        "reason": inferred(reason, reason, evidence, "A template is required because the source does not define one finite, fully stereospecified molecule."),
    }
    m["stereochemistry"]["atom_centers"] = extracted("not assigned in source", None, evidence)
    m["stereochemistry"]["double_bonds"] = extracted(
        double_bonds if double_bonds is not None else "not assigned in source",
        double_bonds, evidence
    )
    m["stereochemistry"]["axial"] = extracted("not assigned in source", None, evidence)
    m["stereochemistry"]["mixture"] = extracted(mixture, mixture, evidence) if mixture else missing(evidence)
    m["stereochemistry"]["assignment"] = extracted(assignment, assignment, evidence) if assignment else missing(evidence)
    m["stereochemistry"]["chirality_status"] = inferred(
        chirality, chirality, evidence, "Status preserves unresolved polymer/product stereochemistry rather than treating missing marks as achiral."
    )
    m["stereochemistry"]["geometry_status"] = inferred(
        geometry, geometry, evidence, "Status follows the explicit or missing E/Z/cis/trans information in the source."
    )
    m["stereochemistry"]["representation_status"] = inferred(
        representation, representation, evidence, "Template does not encode every stereochemical and repeat-unit variable."
    )
    m["identity_note"] = inferred(reason, reason, evidence, "Identity limitation is carried forward from the source representation.")
    return m


def polymer_template(mid, label, aryl, regio, evidence):
    variables = [{
        "label": "n",
        "attachment_points": derived("polymer repeat brackets", "repeat count", evidence, "Repeat count is denoted n in the source drawing."),
        "definition": extracted("n", "degree of polymerization variable", evidence),
        "allowed_values": extracted("not fixed", None, evidence),
        "constraints": extracted("polymer sample has a molecular-weight distribution", "positive integer; distributed sample", evidence),
    }]
    return template_molecule(
        mid, label, f"{regio}-substituted poly({aryl}allene)", evidence,
        {"polymer_repeat_unit": f"source-drawn {regio}-substituted polyaryallene repeat", "aryl_substituent": aryl},
        reason="Polymer length, tacticity, conformer-block distribution, and residual-alkene geometry are not fully specified.",
        chirality="unresolved", geometry="unspecified", representation="partial",
        double_bonds="Residual main-chain C=C is drawn without an E/Z assignment.", variables=variables
    )


def add_charge(p, key, raw, normalized, evidence, reference=None, method="original_extraction", why=None):
    put_fact(p["charges"], key, raw, normalized, evidence, method, why)
    if key in ("equivalents", "mol_percent"):
        if reference is None:
            raise ValueError("loading/equivalents requires reference")
        p["equivalence_reference"] = (
            extracted(raw, reference, evidence) if method == "original_extraction"
            else inherited(raw, reference, evidence, why or "Reference inherited with condition.")
        )


def participant(r, pid, mid, role, evidence, stage="reaction", product_role=None):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    if product_role:
        p["product_role"] = extracted(product_role, product_role, evidence)
    r["participants"].append(p)
    return p


def add_step(r, sid, evidence, pids, *, operation=None, temperature=None, time=None,
             stage="reaction", method="original_extraction", why=None):
    s = st.new_step(sid, evidence, stage=stage)
    s["participant_ids"] = list(pids)
    if operation is not None:
        put_fact(s, "operation", operation, operation, evidence, method, why)
    if temperature is not None:
        raw, norm = temperature
        put_fact(s, "temperature", raw, norm, evidence, method, why)
    if time is not None:
        raw, norm = time
        put_fact(s, "time", raw, norm, evidence, method, why)
    r["steps"].append(s)
    return s


def outcome(r, oid, kind, targets, evidence, raw, normalized, *, role=None,
            yield_type=None, basis=None, scope=None, method="original_extraction", why=None):
    o = st.new_outcome(oid, kind, targets, evidence)
    if method == "original_extraction":
        o["value"] = extracted(raw, normalized, evidence)
    elif method == "inference":
        o["value"] = inferred(raw, normalized, evidence, why)
    elif method == "structure_derivation":
        o["value"] = derived(raw, normalized, evidence, why)
    if role is not None:
        o["product_role"] = extracted(role, role, evidence)
    if yield_type is not None:
        o["yield_type"] = extracted(yield_type, yield_type, evidence)
    if basis is not None:
        o["measurement_basis"] = extracted(basis, basis, evidence)
    if scope is not None:
        o["scope"] = extracted(scope, scope, evidence)
    r["outcomes"].append(o)
    return o


def set_classification(r, evidence, mode="unresolved"):
    c = r["classification"]
    c["mapped_reaction_smiles"] = missing(evidence)
    c["model_votes"] = []
    if mode == "polymerization":
        change = {"change_id": "c1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
                  "before_order": 0, "after_order": 1, "change_type": "formed"}
        c["bond_changes"] = [derived("polymer repeat drawing", change, evidence,
                                      "At least one new inter-monomer C–C sigma connection is explicit in the repeat-unit drawing; full repeat mapping is unresolved.")]
        c["labels"] = derived("new inter-monomer C–C bond", [{"category_number": 2, "code": "C-C", "supporting_change_ids": ["c1"]}], evidence,
                              "Taxonomy label follows the single source-supported formed connection retained in the partial mapping.")
        c["coverage"] = inferred("partial repeat-unit analysis", "partial", evidence,
                                 "Polymer chain ends, full atom correspondence, bond-order changes, and tacticity are not fully mapped.")
        c["limitations"] = ["Partial repeat-unit classification only; no complete atom-mapped polymer reaction is asserted."]
    elif mode == "cycloaddition":
        changes = []
        for i in (1, 2):
            changes.append({"change_id": f"c{i}", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
                            "before_order": 0, "after_order": 1, "change_type": "formed"})
        c["bond_changes"] = [derived("[2+2] cycloaddition product drawing", x, evidence,
                                      "Ring closure creates two new intermolecular C–C sigma bonds; atom mapping is not asserted.") for x in changes]
        c["labels"] = derived("two formed C–C bonds", [{"category_number": 2, "code": "C-C", "supporting_change_ids": ["c1", "c2"]}], evidence,
                              "Both source-supported new ring connections fall in taxonomy category 2.")
        c["coverage"] = inferred("partial cycloaddition graph analysis", "partial", evidence,
                                 "New sigma connections are classified, but full atom mapping and all bond-order changes are unresolved.")
        c["limitations"] = ["No complete atom mapping; exocyclic alkene geometries are not individually assigned."]
    elif mode == "not_applicable":
        c["bond_changes"] = []
        c["labels"] = inferred("no observed net transformation/proposed mechanism only", [], evidence,
                               "No observed product bond formation is classified for this record.")
        c["coverage"] = inferred("not applicable", "not_applicable", evidence,
                                 "Negative control, analytical observation, or proposed mechanism is not an observed net product mapping.")
        c["limitations"] = ["No bond-formation label inferred from an unobserved or non-net transformation."]
    else:
        c["bond_changes"] = []
        c["labels"] = missing(evidence)
        c["coverage"] = inferred("unresolved", "unresolved", evidence,
                                 "Source does not provide enough structure/atom correspondence for a defensible net classification.")
        c["limitations"] = ["Classification unresolved; no label asserted."]


def new_reaction(suffix, evidence, *, record_type="specific_experiment", categories=None,
                 location=None, name=None, transformation=None, classification="unresolved"):
    rid = f"{DOC}:{suffix}"
    r = st.new_reaction(rid, DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = extracted(categories or [], categories or [], evidence)
    r["location"] = extracted(location or suffix, location or suffix, evidence)
    if name:
        r["reaction_name"] = extracted(name, name, evidence)
    if transformation:
        r["transformation_type"] = inferred(transformation, transformation, evidence,
                                               "Transformation type normalized from the source reaction name/drawing.")
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    set_classification(r, evidence, classification)
    return r


def set_standard_solvent(r, evidence, volume=600, method="original_extraction", why=None, solvent="m_c6d6"):
    raw_name = "C6D6" if solvent == "m_c6d6" else "C6H6"
    p = participant(r, f"{r['reaction_id'].split(':')[-1]}_solvent", solvent, "solvent", evidence)
    if volume is not None:
        add_charge(p, "volume", f"{volume} μL {raw_name}", {"number": volume, "unit": "uL"}, evidence, method=method, why=why)
    put_fact(r["solvent_conditions"], "composition", raw_name, [raw_name], evidence, method, why)
    if volume is not None:
        put_fact(r["solvent_conditions"], "total_volume", f"{volume} μL", {"number": volume, "unit": "uL"}, evidence, method, why)
    return p


def add_standard_inputs(r, evidence, *, monomer="m_1", monomer_amount=True, fe=1, hydride="m_nh3bh3",
                        hydride_loading=5, solvent=True, solvent_volume=600, method="original_extraction", why=None):
    tag = r["reaction_id"].split(":")[-1]
    pm = participant(r, f"{tag}_monomer", monomer, "reactant", evidence)
    if monomer_amount:
        add_charge(pm, "amount", "0.5 mmol allene", {"number": 0.5, "unit": "mmol"}, evidence, method=method, why=why)
    if fe is not None:
        pf = participant(r, f"{tag}_fe", "m_fe", "catalyst", evidence)
        add_charge(pf, "mol_percent", f"{fe} mol% [Fe]", {"number": fe, "unit": "mol%"}, evidence, reference=monomer, method=method, why=why)
    if hydride is not None:
        ph = participant(r, f"{tag}_hydride", hydride, "reagent", evidence)
        if isinstance(hydride_loading, tuple) and hydride_loading[1] == "equiv":
            add_charge(ph, "equivalents", f"{hydride_loading[0]} equiv", {"number": hydride_loading[0], "unit": "equiv"}, evidence, reference=monomer, method=method, why=why)
        elif hydride_loading is not None:
            add_charge(ph, "mol_percent", f"{hydride_loading} mol% hydride", {"number": hydride_loading, "unit": "mol%"}, evidence, reference=monomer, method=method, why=why)
    if solvent:
        set_standard_solvent(r, evidence, solvent_volume, method=method, why=why)
    return pm


def add_polymer_products(r, evidence, p23="m_p1_23", p12="m_p1_12"):
    tag = r["reaction_id"].split(":")[-1]
    participant(r, f"{tag}_p23", p23, "product", evidence, product_role="major (2,3)-substituted polymer")
    participant(r, f"{tag}_p12", p12, "product", evidence, product_role="minor (1,2)-substituted polymer")


def add_polymer_metrics(r, evidence, monomer, p23, p12, *, conversion, ratio=None, mn=None, dispersity=None, dp=None, tg=None):
    outcome(r, "conversion", "conversion", [monomer], evidence, f"{conversion}%", {"number": conversion, "unit": "%"},
            basis="1H NMR spectroscopy", scope="allene conversion")
    if ratio is not None:
        outcome(r, "regioisomer_ratio", "selectivity", [p23, p12], evidence, f"{ratio[0]}:{ratio[1]}",
                {"ratio": list(ratio), "order": [p23, p12]}, scope=f"{p23}:{p12}")
    if mn is not None:
        raw, norm = mn if isinstance(mn, tuple) else (f"{mn} kDa", {"number": mn, "reported_unit": "kDa"})
        if isinstance(norm, dict) and norm.get("unit") == "kDa":
            norm = {**norm, "reported_unit": norm["unit"]}
            del norm["unit"]
        outcome(r, "mn", "other", [p23, p12], evidence, raw, norm, scope="number-average molecular mass")
    if dispersity is not None:
        outcome(r, "dispersity", "other", [p23, p12], evidence, str(dispersity), dispersity, scope="Ð")
    if dp is not None:
        outcome(r, "degree_of_polymerization", "other", [p23, p12], evidence, str(dp), dp, scope="DP")
    if tg is not None:
        outcome(r, "glass_transition_temperature", "other", [p23, p12], evidence, f"{tg} °C",
                {"number": tg, "unit": "°C"}, basis="DSC", scope="Tg")


def table_reaction(suffix, evidence, *, hydride, conversion, ratio=None, mn=None, dispersity=None, dp=None,
                   fe=1, hydride_loading=5, temp="RT", hours=16, operation="polymerization", note=None):
    r = new_reaction(suffix, evidence, categories=["optimization"], location=suffix,
                     name="phenylallene polymerization", transformation="allene polymerization", classification="polymerization")
    add_standard_inputs(r, evidence, fe=fe, hydride=hydride, hydride_loading=hydride_loading)
    if conversion > 0:
        add_polymer_products(r, evidence)
    add_step(r, "reaction", evidence, [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"],
             operation=operation, temperature=(temp, temp), time=(f"{hours} h", {"number": hours, "unit": "h"}))
    if conversion > 0:
        add_polymer_metrics(r, evidence, "m_1", "m_p1_23", "m_p1_12", conversion=conversion,
                            ratio=ratio, mn=mn, dispersity=dispersity, dp=dp)
    else:
        outcome(r, "conversion", "conversion", ["m_1"], evidence, f"{conversion}%", {"number": conversion, "unit": "%"},
                basis="1H NMR spectroscopy", scope="allene conversion")
        outcome(r, "no_polymerization", "no_reaction", ["m_1"], evidence, "0% conversion", True, scope="polymerization")
        set_classification(r, evidence, "not_applicable")
    r["procedure"]["original_text"] = extracted(
        "Conditions: 600 μL C6D6, 0.5 mmol 1, 1 mol% [Fe], 5 mol% hydride, 16 h, RT; row-specific exceptions apply.",
        "Table 1 conditions with row-specific exception", evidence)
    if note:
        r["notes"].append(inferred(note, note, evidence, "Note records a source-specific limitation or exception."))
    return r


def main():
    dataset = json.loads((RUN / "dataset.json").read_text())
    bundle = st.empty(dataset["documents"])
    molecules = []

    # Source-drawn allene monomers (Table 3 legend fixes the substitution pattern).
    monomers = [
        ("m_1", "1", "phenylallene", "C=C=Cc1ccccc1", evs(1, 2, 3, 5)),
        ("m_2", "2", "1-(prop-1,2-dien-1-yl)-4-methylbenzene (p-methylphenylallene)", "C=C=Cc1ccc(C)cc1", ev(4)),
        ("m_3", "3", "1-(prop-1,2-dien-1-yl)-3-methylbenzene (m-methylphenylallene)", "C=C=Cc1cccc(C)c1", ev(4)),
        ("m_4", "4", "1-(prop-1,2-dien-1-yl)-2-methylbenzene (o-methylphenylallene)", "C=C=Cc1ccccc1C", ev(4)),
        ("m_5", "5", "1-methoxy-4-(prop-1,2-dien-1-yl)benzene (p-methoxyphenylallene)", "C=C=Cc1ccc(OC)cc1", ev(4)),
        ("m_6", "6", "1-fluoro-4-(prop-1,2-dien-1-yl)benzene (p-fluorophenylallene)", "C=C=Cc1ccc(F)cc1", ev(4)),
        ("m_7", "7", "1-chloro-4-(prop-1,2-dien-1-yl)benzene (p-chlorophenylallene)", "C=C=Cc1ccc(Cl)cc1", ev(4)),
    ]
    for row in monomers:
        molecules.append(concrete_molecule(*row, note="Monosubstituted terminal allene; source gives no axial assignment and the CH2 terminus prevents allene axial chirality."))

    # Concrete small molecules used as inputs/solvents.
    small = [
        ("m_hbpin", "HBpin", "pinacolborane", "[H]B1OC(C)(C)C(C)(C)O1", ev(2)),
        ("m_nh3bh3", "H3N·BH3", "ammonia borane", "[BH3-][NH3+]", evs(1, 2, 3, 4, 5)),
        ("m_me2hnbh3", "Me2HN·BH3", "dimethylamine borane", "[BH3-][NH+](C)C", ev(2)),
        ("m_hsioet3", "HSi(OEt)3", "triethoxysilane", "[H][Si](OCC)(OCC)OCC", ev(2)),
        ("m_c6d6", "C6D6", "benzene-d6", "[2H]c1c([2H])c([2H])c([2H])c([2H])c1[2H]", evs(2, 3, 4, 5)),
        ("m_c6h6", "C6H6", "benzene", "c1ccccc1", ev(5)),
        ("m_meoh", "MeOH", "methanol", "CO", evs(3, 4)),
        ("m_h2", "H2", "hydrogen", "[H][H]", ev(3)),
        ("m_pme3", "PMe3", "trimethylphosphine", "CP(C)C", ev(3)),
        ("m_cyclohexylallene", "cyclohexylallene", "cyclohexylallene", "C=C=CC1CCCCC1", ev(3)),
        ("m_methoxyallene", "methoxyallene", "methoxyallene", "COC=C=C", ev(3)),
        ("m_1d2", "1-d2", "dideuterated phenylallene 1-d2", "[2H]C([2H])=C=Cc1ccccc1", ev(3)),
    ]
    for row in small:
        molecules.append(concrete_molecule(*row))

    # Complexes and source labels that cannot be safely flattened to one concrete SMILES.
    molecules += [
        template_molecule("m_fe", "[Fe]", "iron(II) β-diketiminate pre-catalyst [Fe]", evs(1, 2, 3, 4, 5),
                          "source-drawn Fe(II) β-diketiminate complex bearing aryl/isopropyl and trimethylsilyl substituents",
                          reason="Organometallic coordination and exact ligand graph are retained from the drawing rather than guessed as SMILES."),
        template_molecule("m_feh_dimer", "([Fe]H)2", "bridged iron hydride dimer ([Fe]H)2", evs(3, 4),
                          "source-drawn hydride-bridged dimer of the β-diketiminate iron complex",
                          reason="Bridging hydrides and organometallic bonding are not encoded as a concrete source-verified SMILES."),
        template_molecule("m_feI", "[Fe](I)", "Fe(I)-η6-arene species [Fe](I)", evs(4, 5),
                          "source-drawn Fe(I)-η6-arene analogue",
                          reason="Organometallic hapticity and complete structure are not safely representable from the main-text drawing."),
        template_molecule("m_me2dnbh3", "Me2DN·BH3", "N-deuterated dimethylamine borane", evs(3, 5),
                          "dimethylamine borane with N-D label",
                          reason="Isotopologue identity is explicit but source-safe charged/coordinate SMILES is not asserted."),
        template_molecule("m_me2hnbd3", "Me2HN·BD3", "B-deuterated dimethylamine borane", evs(3, 4, 5),
                          "dimethylamine borane with BD3 label",
                          reason="Isotopologue identity is explicit but source-safe charged/coordinate SMILES is not asserted."),
        template_molecule("m_radical_clock", "radical clock", "radical clock (identity not stated in main text)", ev(3),
                          "unidentified radical clock reagent",
                          reason="Main text reports a radical-clock experiment but not the clock structure."),
        template_molecule("m_hydroborated", "hydroborated product", "phenylallene hydroboration product (regiochemistry unresolved)", ev(2),
                          "hydroborated phenylallene product; exact connectivity not shown in the table",
                          reason="Table footnote reports 7% hydroborated product without a source structure or regioisomer assignment."),
    ]

    # Polymer regioisomers; P1 appears on pages 1–5 and P2–P7 on page 4.
    molecules.append(polymer_template("m_p1_23", "P1_2,3", "phenyl", "2,3", evs(1, 2, 3, 4, 5)))
    molecules.append(polymer_template("m_p1_12", "P1_1,2", "phenyl", "1,2", evs(2, 3, 5)))
    molecules.append(template_molecule("m_p1_unresolved", "P1 (regioisomer unresolved)", "polyphenylallene P1, regioisomer distribution unreported", ev(2),
                                       "polyphenylallene formed without reported regioisomer assignment",
                                       reason="Table entry reports conversion but dashes for regioisomer ratio and polymer metrics.",
                                       chirality="unresolved", geometry="unspecified", representation="partial",
                                       double_bonds="Residual polymer alkene geometry not assigned."))
    aryl_names = {2: "p-methylphenyl", 3: "m-methylphenyl", 4: "o-methylphenyl", 5: "p-methoxyphenyl", 6: "p-fluorophenyl", 7: "p-chlorophenyl"}
    for i in range(2, 8):
        molecules.append(polymer_template(f"m_p{i}_23", f"P{i}_2,3", aryl_names[i], "2,3", ev(4)))
        molecules.append(polymer_template(f"m_p{i}_12", f"P{i}_1,2", aryl_names[i], "1,2", ev(4)))
    molecules.append(template_molecule("m_p1d2", "deuterium-retaining P1", "polyphenylallene from 1-d2 with retained deuterium", ev(3),
                                       "P1 polymer repeat derived from 1-d2 with source-reported D retention",
                                       reason="Polymer length, regiosequence, tacticity, and full isotope positions are not represented by a finite SMILES.",
                                       chirality="unresolved", geometry="unspecified", representation="partial",
                                       double_bonds="Residual polymer alkene geometry not assigned."))

    # Generic literature structures from Scheme 1.
    molecules.append(template_molecule("m_arylallene_generic", "Ar-allene", "generic arylallene", ev(1),
                                       "Ar-substituted terminal allene",
                                       reason="Scheme 1 summarizes several monomers using variable Ar.",
                                       chirality="unresolved", geometry="not_applicable", representation="partial",
                                       variables=[{"label": "Ar", "attachment_points": extracted("Ar–allene", "terminal allene carbon", ev(1)),
                                                   "definition": extracted("Ar", "aryl group", ev(1)),
                                                   "allowed_values": extracted("multiple monomers", None, ev(1)),
                                                   "constraints": extracted("not enumerated in this paper", None, ev(1))}]))
    molecules.append(template_molecule("m_lit_polymer_generic", "literature polyaryallene", "generic literature polyaryallene", ev(1),
                                       "source-drawn polyaryallene repeat unit with variable Ar and n",
                                       reason="Scheme 1 is a literature summary and does not enumerate structures in this paper.",
                                       chirality="unresolved", geometry="unspecified", representation="partial"))

    # Explicit cyclobutane diastereomers; exocyclic alkene geometry remains unresolved.
    for code, assignment in (("cis", "cis-1,3 relationship on cyclobutane"), ("trans", "trans-1,3 relationship on cyclobutane")):
        molecules.append(template_molecule(
            f"m_cyclobutane_{code}", f"{code}-cyclobutane dimer",
            f"{code}-1,3-di((Z/E)-benzylidene)cyclobutane", ev(5),
            "1,3-dibenzylidenecyclobutane connectivity with source-specified ring relationship",
            reason="The sample is reported as cis/trans and Z/E-containing isomers; individual exocyclic geometries are not assigned.",
            chirality="unresolved", geometry="partially_specified", representation="partial",
            assignment=assignment, double_bonds="Exocyclic benzylidene units are described as (Z/E) without individual assignments."
        ))

    bundle["molecules"] = molecules
    reactions = []

    # Scheme 1: two cited literature arrows plus the generic reaction reported here.
    for suffix, catalyst, condition, metric in [
        ("scheme1_endo_literature", "[(π-allyl)Ni(OCOCF3)]2", "toluene or THF, 0 °C to RT", "allene/[Ni]=80; Mn=12.6–13 kDa; PDI=1.1"),
        ("scheme1_cui_literature", "(flu-py)In(CH2TMS)2/[Ph3C][B(C6F5)4]/AliBu3", "toluene, RT", "allene/[Gd]=1000 as printed; Mn=22–40 kDa; PDI=1.6–2.1"),
    ]:
        r = new_reaction(suffix, ev(1), record_type="literature_example", categories=["other"], location="Scheme 1",
                         name="previously reported phenylallene/arylallene polymerization", transformation="allene polymerization",
                         classification="polymerization")
        participant(r, f"{suffix}_substrate", "m_arylallene_generic", "reactant", ev(1))
        participant(r, f"{suffix}_product", "m_lit_polymer_generic", "product", ev(1), product_role="literature polymer")
        temperature_norm = None if "0 °C to RT" in condition else condition
        add_step(r, "reaction", ev(1), [p["participant_id"] for p in r["participants"]], operation=f"polymerization with {catalyst}; {condition}", temperature=(condition, temperature_norm))
        outcome(r, "reported_metrics", "other", ["m_lit_polymer_generic"], ev(1), metric, metric, scope="Scheme 1 literature summary")
        r["review"]["unresolved"].append("Cited literature example only; full experimental details and exact monomer identities were not re-extracted from the cited papers.")
        reactions.append(r)
    r = new_reaction("scheme1_this_work_generic", ev(1), record_type="generic_scheme", categories=["other"], location="Scheme 1c",
                     name="iron-catalyzed arylallene polymerization", transformation="allene polymerization", classification="polymerization")
    participant(r, "scheme1_substrate", "m_arylallene_generic", "reactant", ev(1))
    participant(r, "scheme1_fe", "m_fe", "catalyst", ev(1))
    participant(r, "scheme1_hydride", "m_nh3bh3", "reagent", ev(1))
    participant(r, "scheme1_product", "m_lit_polymer_generic", "product", ev(1), product_role="generic polyaryallene")
    add_step(r, "reaction", ev(1), [p["participant_id"] for p in r["participants"]], temperature=("RT", "RT"), operation="polymerization in C6D6")
    r["notes"].append(extracted("7 different monomers; allene/[Fe]=100; Mn=53 to 229 kDa; PDI=1.1 to 2.1", "Scheme 1 summary", ev(1)))
    r["review"]["unresolved"].append("Generic scheme only; individual experiments are recorded separately from Tables 1–3.")
    reactions.append(r)

    # Table 1 (entry 4 is merged with Table 2 entry 1 below).
    reactions.append(table_reaction("table1_entry1_no_hydride", ev(2), hydride=None, conversion=3, fe=1,
                                    operation="polymerization attempt without hydride",
                                    note="Regioisomer ratio and polymer metrics are table dashes; product identity is retained as unresolved polyphenylallene."))
    r = table_reaction("table1_entry2_hbpin_1eq", ev(2), hydride="m_hbpin", hydride_loading=(1, "equiv"),
                       conversion=82, ratio=(8, 1), mn=53, dispersity=1.12, dp=460)
    p = next(p for p in r["participants"] if p["molecule_id"] == "m_hbpin")
    add_charge(p, "amount", "0.5 mmol HBpin", {"number": 0.5, "unit": "mmol"}, ev(2))
    participant(r, "table1_entry2_hydroborated", "m_hydroborated", "product", ev(2), product_role="hydroborated side product")
    outcome(r, "hydroborated_product", "yield", ["m_hydroborated"], ev(2), "7% conversion to hydroborated product observed",
            {"number": 7, "unit": "%"}, role="side product", basis="1H NMR spectroscopy", scope="phenylallene-to-hydroborated-product")
    r["review"]["unresolved"].append("Hydroborated-product connectivity/regiochemistry is not shown in the main-text table.")
    reactions.append(r)
    reactions.append(table_reaction("table1_entry3_hbpin_5mol", ev(2), hydride="m_hbpin", conversion=91,
                                    ratio=(5, 1), mn=152, dispersity=1.36, dp=1310))
    standard = table_reaction("standard_conditions_p1", evs(2, 3), hydride="m_nh3bh3", conversion=93,
                              ratio=(7, 1), mn=189, dispersity=1.19, dp=1630)
    standard["location"] = extracted("Table 1 entry 4; repeated as Table 2 entry 1 and optimized-conditions prose", ["Table 1 entry 4", "Table 2 entry 1"], evs(2, 3))
    standard["evidence_ids"] = evs(2, 3)
    standard["procedure"]["original_text"] = extracted(
        "[Fe] (1 mol%), H3N·BH3 (5 mol%), C6D6 (0.6 mL), RT, 16 h",
        "optimized conditions", evs(2, 3))
    reactions.append(standard)
    reactions.append(table_reaction("table1_entry5_me2hnbh3", ev(2), hydride="m_me2hnbh3", conversion=89,
                                    ratio=(5, 1), mn=10.5, dispersity=6.25, dp=91))
    reactions.append(table_reaction("table1_entry6_hsioet3", ev(2), hydride="m_hsioet3", conversion=92,
                                    ratio=(6, 1), mn=62, dispersity=1.18, dp=530))
    reactions.append(table_reaction("table1_entry7_no_fe", ev(2), hydride="m_nh3bh3", conversion=0, fe=None,
                                    temp="60 °C", hours=72, operation="control without [Fe]"))

    # Table 2 deviations from standard conditions (entry 1 already merged above).
    table2 = [
        ("table2_entry2_no_stirring", "no stirring", 91, (83.9, 1.38), 600, 1, 5, "RT", 16),
        ("table2_entry3_100ul_solvent", "100 μL solvent", 92, ((">250 kDa (footnote; table cell is dash)", {"number": 250, "unit": "kDa", "qualifier": "gt"}), None), 100, 1, 5, "RT", 16),
        ("table2_entry4_no_solvent", "no solvent", 71, (12.6, 1.75), None, 1, 5, "RT", 16),
        ("table2_entry5_80c", "heated at 80 °C", 96, (132, 1.85), 600, 1, 5, "80 °C", 16),
        ("table2_entry6_quench_5min", "quenched after 5 minutes", 19, (10.5, 1.33), 600, 1, 5, "RT", 5/60),
        ("table2_entry7_20mol_fe", "20 mol% catalyst loading", 96, ((">250 kDa (footnote; table cell is dash)", {"number": 250, "unit": "kDa", "qualifier": "gt"}), None), 600, 20, 5, "RT", 16),
    ]
    for suffix, deviation, conv, metrics, solvent_volume, fe_load, hyd_load, temp, hours in table2:
        r = new_reaction(suffix, ev(3), categories=["optimization"], location=f"Table 2, {suffix.split('_entry')[-1]}",
                         name="phenylallene polymerization", transformation="allene polymerization", classification="polymerization")
        add_standard_inputs(r, ev(3), fe=fe_load, hydride="m_nh3bh3", hydride_loading=hyd_load,
                            solvent=solvent_volume is not None, solvent_volume=solvent_volume)
        add_polymer_products(r, ev(3))
        if suffix == "table2_entry6_quench_5min":
            p = participant(r, "table2_entry6_meoh", "m_meoh", "workup_agent", ev(3), stage="workup")
            add_step(r, "reaction", ev(3), [x["participant_id"] for x in r["participants"] if x["stage"] == "reaction"],
                     operation="polymerization before methanol quench", temperature=(temp, temp), time=("5 min", {"number": 5, "unit": "min"}))
            add_step(r, "quench", ev(3), [p["participant_id"]], operation="quenched with methanol", stage="workup")
        else:
            add_step(r, "reaction", ev(3), [x["participant_id"] for x in r["participants"]], operation=deviation,
                     temperature=(temp, temp), time=(f"{hours} h", {"number": hours, "unit": "h"}))
        add_polymer_metrics(r, ev(3), "m_1", "m_p1_23", "m_p1_12", conversion=conv,
                            mn=metrics[0], dispersity=metrics[1])
        r["procedure"]["modifications"] = extracted(deviation, deviation, ev(3))
        if suffix == "table2_entry4_no_solvent":
            mn = next(o for o in r["outcomes"] if o["outcome_id"] == "mn")
            a = extracted("Table cell: 12.6 kDa", {"number": 12.6, "unit": "kDa"}, ev(3))
            b = extracted("Footnote [a]: polymer length was greater than 250 kDa", {"number": 250, "unit": "kDa", "qualifier": "gt"}, ev(3))
            mn["value"] = st.fact([a["raw_value"], b["raw_value"]], None, ev(3), status="conflict", alternatives=[a, b])
            r["notes"].append(extracted("Table 2 entry 4 gives Mn 12.6 kDa while footnote [a] includes entry 4 among conditions >250 kDa.",
                                                "source-internal Mn conflict retained", ev(3)))
            r["review"]["unresolved"].append("Resolve Table 2 entry 4 Mn conflict (12.6 kDa table cell versus >250 kDa footnote) from SI or author clarification; no SI was provided.")
        if suffix == "table2_entry5_80c":
            participant(r, "table2_entry5_cis_dimer", "m_cyclobutane_cis", "product", ev(3), product_role="cyclobutane dimer component")
            participant(r, "table2_entry5_trans_dimer", "m_cyclobutane_trans", "product", ev(3), product_role="cyclobutane dimer component")
            outcome(r, "polymer_fraction", "yield", ["m_p1_23", "m_p1_12"], ev(3), "90% polymer", {"number": 90, "unit": "%"}, role="polymer products", scope="Table 2 footnote [b]")
            outcome(r, "dimer_fraction", "yield", ["m_cyclobutane_cis", "m_cyclobutane_trans"], ev(3), "6% dimer", {"number": 6, "unit": "%"}, role="dimer side product", scope="Table 2 footnote [b]")
            r["review"]["unresolved"].append("Measurement basis for the 90% polymer/6% dimer split is not separately stated; only overall conversion is explicitly assigned to 1H NMR.")
        reactions.append(r)

    # Table 3 substrate scope.
    scope_rows = {
        2: (94, (8, 1), 229, 1.63, 1760, 86.4, 1),
        3: (96, (7, 1), 72, 1.17, 550, 63.0, 1),
        4: (91, (5, 1), 185, 1.62, 1420, 75.4, 1),
        5: (98, (11, 1), 174, 2.05, 1190, 81.6, 1),
        6: (77, (13, 1), 53, 1.12, 400, 57.0, 5),
        7: (74, (6, 1), 166, 1.67, 1110, 97.1, 5),
    }
    for i, (conv, ratio, mn, disp, dp, tg, fe_load) in scope_rows.items():
        suffix = f"table3_entry{i-1}_monomer{i}"
        r = new_reaction(suffix, ev(4), categories=["substrate_scope"], location=f"Table 3 entry {i-1}",
                         name="functionalized arylallene polymerization", transformation="allene polymerization", classification="polymerization")
        add_standard_inputs(r, ev(4), monomer=f"m_{i}", fe=fe_load, hydride="m_nh3bh3", hydride_loading=5)
        add_polymer_products(r, ev(4), f"m_p{i}_23", f"m_p{i}_12")
        add_step(r, "reaction", ev(4), [p["participant_id"] for p in r["participants"]], operation="polymerization",
                 temperature=("RT", "RT"), time=("16 h", {"number": 16, "unit": "h"}))
        add_polymer_metrics(r, ev(4), f"m_{i}", f"m_p{i}_23", f"m_p{i}_12", conversion=conv, ratio=ratio,
                            mn=mn, dispersity=disp, dp=dp, tg=tg)
        r["procedure"]["original_text"] = extracted(
            "600 μL C6D6, 0.5 mmol allene, 1 mol% [Fe], 5 mol% H3N·BH3, 16 h, RT; entries 5 and 6 use 5 mol% [Fe].",
            "Table 3 conditions with footnote [b] exception", ev(4))
        reactions.append(r)

    # Non-aryl allene controls stated in prose; exact charges/conditions are not restated.
    for mid, label in (("m_cyclohexylallene", "cyclohexylallene"), ("m_methoxyallene", "methoxyallene")):
        r = new_reaction(f"control_no_polymerization_{label}", ev(3), categories=["control"], location="page 3 prose",
                         name="non-aryl allene polymerization test", transformation="allene polymerization", classification="not_applicable")
        participant(r, f"{label}_substrate", mid, "reactant", ev(3))
        add_step(r, "polymerization_test", ev(3), [f"{label}_substrate"], operation="polymerization test of non-aryl allene")
        outcome(r, "no_polymerization", "no_reaction", [mid], ev(3), "We do not observe any polymerization", True, scope=label)
        r["review"]["unresolved"].append("Exact charges and conditions for this prose-only non-aryl control are not restated in the main text; optimized conditions were not silently inherited.")
        reactions.append(r)

    # Mechanistic experiments in Scheme 2 and surrounding prose.
    r = new_reaction("mechanistic_kinetic_catalyst_order", ev(3), categories=["mechanistic_study"], location="page 3 prose",
                     name="kinetic analysis", transformation="kinetic study", classification="not_applicable")
    participant(r, "kinetic_monomer", "m_1", "reactant", ev(3))
    participant(r, "kinetic_fe", "m_fe", "catalyst", ev(3))
    add_step(r, "kinetic_analysis", ev(3), ["kinetic_monomer", "kinetic_fe"], operation="kinetic analysis of the polymerization")
    outcome(r, "catalyst_order", "other", ["m_fe"], ev(3), "first order with respect to the catalyst", 1, scope="kinetic order in [Fe]")
    r["review"]["unresolved"].append("Kinetic concentrations, temperature, rate data, and fitting method are not reported in the main text.")
    reactions.append(r)

    r = new_reaction("mechanistic_pme3_poisoning", ev(3), categories=["mechanistic_study", "control"], location="page 3 prose",
                     name="PMe3 poisoning experiment", transformation="catalyst poisoning study", classification="unresolved")
    participant(r, "pme3_monomer", "m_1", "reactant", ev(3))
    participant(r, "pme3_fe", "m_fe", "catalyst", ev(3))
    participant(r, "pme3_additive", "m_pme3", "additive", ev(3))
    add_step(r, "poisoning_test", ev(3), ["pme3_monomer", "pme3_fe", "pme3_additive"], operation="PMe3 poisoning experiment on the polymerization")
    outcome(r, "interpretation", "other", [], ev(3), "reaction is not nanoparticle mediated", "authors' interpretation", scope="PMe3 poisoning experiment")
    r["review"]["unresolved"].append("PMe3 loading, reaction conditions, and quantitative outcome are not reported in the main text.")
    reactions.append(r)

    r = new_reaction("mechanistic_radical_clock", ev(3), categories=["mechanistic_study", "control"], location="page 3 prose",
                     name="radical clock experiment", transformation="radical-clock study", classification="unresolved")
    participant(r, "clock_monomer", "m_1", "reactant", ev(3))
    participant(r, "clock_reagent", "m_radical_clock", "additive", ev(3))
    add_step(r, "radical_clock_test", ev(3), ["clock_monomer", "clock_reagent"], operation="addition of a radical clock during polymerization")
    outcome(r, "rate_effect", "other", [], ev(3), "does not significantly alter rate", "no significant rate alteration", scope="radical clock addition")
    r["review"]["unresolved"].append("Radical-clock identity, dose, conditions, and numerical rate comparison are not stated in the main text.")
    reactions.append(r)

    r = new_reaction("scheme2a_deuterated_monomer", ev(3), categories=["mechanistic_study"], location="Scheme 2a and page 3 prose",
                     name="polymerization of deuterium-labelled phenylallene", transformation="allene polymerization", classification="polymerization")
    add_standard_inputs(r, ev(3), monomer="m_1d2", method="inherited_general_condition", why="Scheme 2a explicitly says standard conditions; the Table 2 standard conditions are on the same page.")
    participant(r, "scheme2a_product", "m_p1d2", "product", ev(3), product_role="deuterium-retaining polymer")
    add_step(r, "reaction", ev(3), [p["participant_id"] for p in r["participants"]], operation="polymerization under standard conditions",
             temperature=("standard conditions: RT", "RT"), time=("standard conditions: 16 h", {"number": 16, "unit": "h"}),
             method="inherited_general_condition", why="Inherited because Scheme 2a explicitly invokes standard conditions.")
    outcome(r, "isotope_retention", "other", ["m_p1d2"], ev(3), "clean isotopic retention; no intramolecular scrambling",
            "deuterium retained without detectable intramolecular scrambling", basis="1H and 2H NMR spectroscopy", scope="polymer from 1-d2")
    r["notes"].append(extracted("3.46 ppm methylene integral matches residual protons; 6.24 ppm α proton integrates to 1H; 2H NMR peaks only at 3.70–3.45 ppm.",
                                        "reported isotope-retention observations", ev(3)))
    reactions.append(r)

    for hydride, code in (("m_me2dnbh3", "me2dnbh3"), ("m_me2hnbd3", "me2hnbd3")):
        r = new_reaction(f"scheme2b_{code}", ev(3), categories=["mechanistic_study"], location="Scheme 2b and page 3 prose",
                         name="polymerization with deuterated hydride source", transformation="allene polymerization", classification="polymerization")
        add_standard_inputs(r, ev(3), hydride=hydride, hydride_loading=5, method="inherited_general_condition",
                            why="Scheme 2b states standard conditions with the named deuterated hydride source.")
        add_polymer_products(r, ev(3))
        add_step(r, "reaction", ev(3), [p["participant_id"] for p in r["participants"]], operation="polymerization followed by early quench",
                 temperature=("standard conditions: RT", "RT"), method="inherited_general_condition",
                 why="Temperature inherited from explicitly invoked standard conditions; quench time is not reported.")
        outcome(r, "deuterium_incorporation", "other", ["m_p1_23", "m_p1_12"], ev(3), "no apparent D-incorporation in the polymer",
                "no apparent deuterium incorporation", basis="NMR (specific method not separately stated)", scope="low-Mn polymer after quench")
        r["review"]["unresolved"].append("Early-quench time and workup identity are not given for Scheme 2b in the main text.")
        reactions.append(r)

    r = new_reaction("scheme2c_feh_dimer", evs(3, 4), categories=["mechanistic_study", "control"], location="Scheme 2c and page 4 prose",
                     name="([Fe]H)2 catalyst test", transformation="allene polymerization", classification="polymerization")
    pm = participant(r, "feh_monomer", "m_1", "reactant", evs(3, 4))
    pc = participant(r, "feh_catalyst", "m_feh_dimer", "catalyst", evs(3, 4))
    add_charge(pc, "mol_percent", "2.5 mol% ([Fe]H)2", {"number": 2.5, "unit": "mol%"}, ev(3), reference="m_1")
    set_standard_solvent(r, ev(3), 600)
    participant(r, "feh_p1", "m_p1_unresolved", "product", evs(3, 4), product_role="low-conversion polymer")
    add_step(r, "reaction", evs(3, 4), [p["participant_id"] for p in r["participants"]], operation="catalyst test without H3N·BH3",
             temperature=("RT", "RT"), time=("16 h", {"number": 16, "unit": "h"}))
    outcome(r, "conversion", "conversion", ["m_1"], evs(3, 4), "5% conversion", {"number": 5, "unit": "%"}, scope="phenylallene conversion")
    r["review"]["unresolved"].append("Measurement basis for the 5% conversion is not separately stated in the prose.")
    reactions.append(r)

    r = new_reaction("scheme2c_hydrogen_control", ev(3), categories=["mechanistic_study", "control"], location="Scheme 2c",
                     name="[Fe] under hydrogen atmosphere", transformation="allene polymerization attempt", classification="not_applicable")
    participant(r, "h2ctrl_monomer", "m_1", "reactant", ev(3))
    p = participant(r, "h2ctrl_fe", "m_fe", "catalyst", ev(3)); add_charge(p, "mol_percent", "1 mol% [Fe]", {"number": 1, "unit": "mol%"}, ev(3), reference="m_1")
    p = participant(r, "h2ctrl_h2", "m_h2", "reagent", ev(3)); add_charge(p, "equivalents", "1 atm H2", None, ev(3), reference="m_1")
    # Pressure belongs to the step; the participant equivalence entry above is removed because 1 atm is not an equivalent.
    p["charges"]["equivalents"] = missing(ev(3)); p["equivalence_reference"] = missing(ev(3))
    set_standard_solvent(r, ev(3), 600)
    s = add_step(r, "reaction", ev(3), [p["participant_id"] for p in r["participants"]], operation="reaction under hydrogen",
                 temperature=("RT", "RT"), time=("16 h", {"number": 16, "unit": "h"}))
    s["atmosphere"] = extracted("H2", "hydrogen", ev(3)); s["pressure"] = extracted("1 atm H2", {"number": 1, "reported_unit": "atm"}, ev(3))
    outcome(r, "no_reaction", "no_reaction", ["m_1"], ev(3), "No reaction in absence of hydride source", True, scope="Scheme 2c")
    reactions.append(r)

    # Scheme 3 proposed pathway and the separately observed Fe(I) experiment.
    r = new_reaction("scheme3_proposed_mechanism", evs(4, 5), record_type="proposed_mechanism", categories=["mechanistic_study"],
                     location="Scheme 3", name="postulated catalyst activation and polymerization process",
                     transformation="proposed Fe(I)/Fe(III) chain-growth mechanism", classification="not_applicable")
    participant(r, "scheme3_fe", "m_fe", "catalyst", ev(4))
    participant(r, "scheme3_monomer", "m_arylallene_generic", "reactant", ev(4))
    participant(r, "scheme3_meoh", "m_meoh", "reagent", ev(4))
    r["notes"].append(extracted("reduction to Fe(I) species I; terminal C–H oxidative addition to II; insertion into Fe–C bond to III; MeOH termination/quench",
                                        "authors' proposed sequence, not experimentally proven elementary steps", evs(4, 5)))
    r["review"]["unresolved"].append("Proposed mechanism only; intermediates I–III are not recorded as independently isolated experiments.")
    reactions.append(r)

    r = new_reaction("fei_arene_catalyst_experiment", ev(5), categories=["mechanistic_study", "control"], location="page 5 prose adjacent to Scheme 3",
                     name="polymerization catalyzed by isolated analogous Fe(I)-η6-arene species", transformation="allene polymerization", classification="polymerization")
    participant(r, "fei_monomer", "m_1", "reactant", ev(5))
    participant(r, "fei_cat", "m_feI", "catalyst", ev(5))
    add_polymer_products(r, ev(5))
    add_step(r, "reaction", ev(5), [p["participant_id"] for p in r["participants"]], operation="polymerization with Fe(I)-η6-arene catalyst without hydride source")
    outcome(r, "conversion", "conversion", ["m_1"], ev(5), "90% conversion to P1", {"number": 90, "unit": "%"}, scope="phenylallene conversion")
    outcome(r, "ratio", "selectivity", ["m_p1_23", "m_p1_12"], ev(5), "6:1 P1_2,3:P1_1,2", {"ratio": [6, 1], "order": ["m_p1_23", "m_p1_12"]}, scope="polymer regioisomer ratio")
    outcome(r, "mn", "other", ["m_p1_23", "m_p1_12"], ev(5), "Mn=190 kDa", {"number": 190, "reported_unit": "kDa"}, scope="number-average molecular mass")
    outcome(r, "dispersity", "other", ["m_p1_23", "m_p1_12"], ev(5), "Ð=1.51", 1.51, scope="dispersity")
    r["notes"].append(extracted("without a hydride source present", "no hydride co-catalyst", ev(5)))
    r["review"]["unresolved"].append("Fe(I) catalyst loading, solvent, temperature, time, and measurement bases are not restated in the main text.")
    reactions.append(r)

    # Scheme 4 isolated cyclobutane experiment.
    r = new_reaction("scheme4_cyclobutane", ev(5), categories=["optimization"], location="Scheme 4",
                     name="[2+2] cycloaddition of phenylallene", transformation="[2+2] allene dimerization", classification="cycloaddition")
    pm = participant(r, "scheme4_monomer", "m_1", "reactant", ev(5))
    pf = participant(r, "scheme4_fe", "m_fe", "catalyst", ev(5)); add_charge(pf, "mol_percent", "5 mol% [Fe]", {"number": 5, "unit": "mol%"}, ev(5), reference="m_1")
    ph = participant(r, "scheme4_hydride", "m_nh3bh3", "reagent", ev(5)); add_charge(ph, "mol_percent", "10 mol% H3N·BH3", {"number": 10, "unit": "mol%"}, ev(5), reference="m_1")
    set_standard_solvent(r, ev(5), None)
    participant(r, "scheme4_cis", "m_cyclobutane_cis", "product", ev(5), product_role="cis cyclobutane isomer")
    participant(r, "scheme4_trans", "m_cyclobutane_trans", "product", ev(5), product_role="trans cyclobutane isomer")
    add_step(r, "reaction", ev(5), [p["participant_id"] for p in r["participants"]], operation="[2+2] cycloaddition",
             temperature=("80 °C", {"number": 80, "unit": "°C"}), time=("16 h", {"number": 16, "unit": "h"}))
    outcome(r, "isolated_yield", "yield", ["m_cyclobutane_cis", "m_cyclobutane_trans"], ev(5), "27% isolated yield",
            {"number": 27, "unit": "%"}, role="combined cyclobutane isomers", yield_type="isolated yield", basis="isolated", scope="combined cis/trans product")
    outcome(r, "cis_trans_ratio", "selectivity", ["m_cyclobutane_cis", "m_cyclobutane_trans"], ev(5), "2:1 cis:trans",
            {"ratio": [2, 1], "order": ["m_cyclobutane_cis", "m_cyclobutane_trans"]}, scope="cis:trans")
    r["review"]["unresolved"].append("Scheme 4 does not state phenylallene charge or C6D6 volume; exocyclic Z/E components are not individually assigned.")
    reactions.append(r)

    r = new_reaction("cyclobutane_dilute_syringe_attempt", ev(5), categories=["optimization"], location="page 5 prose",
                     name="controlled dilute addition for cyclobutane formation", transformation="[2+2] allene dimerization", classification="cycloaddition")
    participant(r, "dilute_monomer", "m_1", "reactant", ev(5))
    pf = participant(r, "dilute_fe", "m_fe", "catalyst", ev(5)); add_charge(pf, "mol_percent", "1 mol% [Fe]", {"number": 1, "unit": "mol%"}, ev(5), reference="m_1")
    ph = participant(r, "dilute_hydride", "m_nh3bh3", "reagent", ev(5)); add_charge(ph, "mol_percent", "5 mol% H3N·BH3", {"number": 5, "unit": "mol%"}, ev(5), reference="m_1")
    ps = participant(r, "dilute_feed_solvent", "m_c6h6", "solvent", ev(5))
    add_charge(ps, "volume", "670 μL C6H6 in the monomer feed", {"number": 670, "unit": "uL"}, ev(5))
    ps2 = participant(r, "dilute_receiving_solvent", "m_c6h6", "solvent", ev(5))
    add_charge(ps2, "volume", "6 mL C6H6 in the receiving flask", {"number": 6, "unit": "mL"}, ev(5))
    participant(r, "dilute_cis", "m_cyclobutane_cis", "product", ev(5), product_role="cyclobutane target component")
    participant(r, "dilute_trans", "m_cyclobutane_trans", "product", ev(5), product_role="cyclobutane target component")
    add_polymer_products(r, ev(5))
    put_fact(r["solvent_conditions"], "composition", "C6H6", ["C6H6"], ev(5))
    put_fact(r["solvent_conditions"], "total_volume", "670 μL monomer solution added to 6 mL receiving solution",
             {"monomer_solution": {"number": 670, "unit": "uL"}, "receiving_solution": {"number": 6, "unit": "mL"}}, ev(5))
    r["special_conditions"]["flow_parameters"] = extracted("20 μL/min syringe-pump addition", {"flow_rate": {"number": 20, "unit": "uL/min"}}, ev(5))
    add_step(r, "addition", ev(5), ["dilute_monomer", "dilute_feed_solvent", "dilute_receiving_solvent"], operation="controlled dilute syringe-pump addition")
    outcome(r, "yield_effect", "other", ["m_cyclobutane_cis", "m_cyclobutane_trans"], ev(5), "failed to give a substantial increase in yield",
            "no substantial yield increase", scope="cyclobutane product")
    outcome(r, "polymer_dominance", "other", ["m_p1_23", "m_p1_12"], ev(5), "polymerization still dominates", "polymerization dominant", scope="competing polymerization")
    r["review"]["unresolved"].append("Phenylallene amount, temperature, total reaction time, and numerical product yields are not reported for the syringe-pump example.")
    reactions.append(r)

    for hydride, code in (("m_me2dnbh3", "me2dnbh3"), ("m_me2hnbd3", "me2hnbd3")):
        r = new_reaction(f"cyclobutane_deuterated_{code}", ev(5), categories=["mechanistic_study"], location="page 5 prose after Scheme 4",
                         name="[2+2] cycloaddition with deuterated hydride source", transformation="[2+2] allene dimerization", classification="cycloaddition")
        participant(r, f"{code}_monomer", "m_1", "reactant", ev(5))
        participant(r, f"{code}_fe", "m_fe", "catalyst", ev(5))
        participant(r, f"{code}_hydride", hydride, "reagent", ev(5))
        participant(r, f"{code}_cis", "m_cyclobutane_cis", "product", ev(5), product_role="cis cyclobutane component")
        participant(r, f"{code}_trans", "m_cyclobutane_trans", "product", ev(5), product_role="trans cyclobutane component")
        add_step(r, "reaction", ev(5), [p["participant_id"] for p in r["participants"]], operation="undertaking this [2+2] cycloaddition with deuterated hydride source")
        outcome(r, "deuterium_incorporation", "other", ["m_cyclobutane_cis", "m_cyclobutane_trans"], ev(5), "does not lead to 2H incorporation in the product",
                "no deuterium incorporation observed", scope="cyclobutane product")
        r["review"]["unresolved"].append("Exact loadings, temperature, time, yield, selectivity, and analytical method are not restated for this isotopic cycloaddition experiment.")
        reactions.append(r)

    bundle["reactions"] = reactions
    bundle["source_context"] = [
        {"context_id": f"{DOC}:figure1_dft", "document_id": DOC, "evidence_ids": ev(2),
         "kind": "computational_conformer_analysis", "summary": "Seven-repeat-unit DFT models: 0.0, +2.7, +6.5, and +8.9 kcal mol−1; not separate synthetic experiments."},
        {"context_id": f"{DOC}:p1_nmr_conformers", "document_id": DOC, "evidence_ids": ev(2),
         "kind": "characterization", "summary": "P1 main-chain proton environments 1.7:1.1:1 and NOE/NOESY discussion; no extra synthesis record created."},
        {"context_id": f"{DOC}:p1_maldi", "document_id": DOC, "evidence_ids": ev(3),
         "kind": "characterization", "summary": "MALDI-ToF repeat units 116 Da and proteo/allenyl end groups; no CH2SiMe3, BH2, or NH2 fragments incorporated."},
        {"context_id": f"{DOC}:p1_dsc_waxd", "document_id": DOC, "evidence_ids": ev(3),
         "kind": "characterization", "summary": "P1 Tg 64.8 °C by DSC; no melting transition; WAXD indicates amorphous material."},
        {"context_id": f"{DOC}:precatalyst_residual", "document_id": DOC, "evidence_ids": ev(5),
         "kind": "mechanistic_observation", "summary": "In situ wide-sweep NMR reports significant unreacted pre-catalyst after full substrate conversion; no separate synthesis record."},
        {"context_id": f"{DOC}:page6_references", "document_id": DOC, "evidence_ids": ev(6),
         "kind": "bibliography", "summary": "Page 6 contains references and manuscript dates; no reaction arrow, table row, or new experiment."},
    ]

    # Raw ledger is kept separate from the formal submission.
    raw = {
        "document_id": DOC,
        "title": dataset["documents"][0]["title"],
        "pages_visually_inspected": [1, 2, 3, 4, 5, 6],
        "detected_record_ids": [r["reaction_id"] for r in reactions],
        "molecule_ids": [m["molecule_id"] for m in molecules],
        "notes": [
            "Table 1 entry 4 and Table 2 entry 1 were treated as one repeated standard-condition experiment.",
            "Table 2 entry 4 Mn conflict retained: 12.6 kDa table cell versus >250 kDa footnote.",
            "No SI file was associated; all missing SI-only details remain unresolved.",
            "All formal records are machine submissions and must remain needs_review after merge.",
        ],
    }
    (RUN / "raw_model/extraction_ledger.json").write_text(json.dumps(raw, ensure_ascii=False, indent=2) + "\n")
    (RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")


if __name__ == "__main__":
    main()
