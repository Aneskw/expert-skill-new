#!/usr/bin/env python3
"""Build the page-1--6 v3.1 submission from manually inspected source images."""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st


DOC = "doc_fb3fead4a044bff7"
RUN = ROOT / "runs/current"
BLANK = RUN / "raw_model/blank_template.json"


def ev(page: int, kind: str = "image") -> str:
    return f"{DOC}:page{page}:{kind}"


def q(number, unit, **extra):
    value = {"number": number, "unit": unit}
    value.update(extra)
    return value


def f(raw=None, normalized=None, evidence=(), method="original_extraction", status=None,
      derivation=None):
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   derivation=derivation)


def missing(status="unreported"):
    return st.fact(status=status, method="not_applicable" if status == "not_applicable" else "original_extraction")


bundle = json.loads(BLANK.read_text())
molecules = {}


def add_molecule(mid, label, page, name=None, smiles=None, chirality="unresolved",
                 geometry="unresolved", representation="unresolved", assignment=None,
                 mixture=None, note=None, extra_evidence=()):
    if mid in molecules:
        return mid
    evidence = [ev(page, "image"), *extra_evidence]
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    m["names"]["original"] = f(name or label, name or label, evidence)
    m["names"]["english"] = f(name or label, name or label, evidence)
    m["cas_numbers"] = missing()
    if smiles:
        diag = st.stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = f(f"source drawing for {label}", smiles, evidence,
                         method="structure_derivation",
                         derivation="Candidate SMILES transcribed from the inspected source drawing; parsing does not certify identity.")
        if representation == "encoded":
            m["stereochemistry"]["isomeric_smiles"] = f(
                "tool canonicalization", diag["canonical_isomeric_smiles"], evidence,
                method="structure_derivation",
                derivation="Computed with standard stereo from the candidate source-derived SMILES.")
            m["stereochemistry"]["identity_key"] = f(
                "standard stereo graph key", diag["stereo_graph_key"], evidence,
                method="structure_derivation",
                derivation="Computed with standard stereo; this is not a sample identity.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing()
    m["stereochemistry"]["chirality_status"] = f(
        chirality, chirality, evidence, method="inference" if chirality in {"achiral", "unresolved"} else "original_extraction",
        derivation="Status assigned after inspection of the source structure; unresolved is retained where a full graph/configuration was not safely encoded." if chirality in {"achiral", "unresolved"} else None)
    m["stereochemistry"]["geometry_status"] = f(
        geometry, geometry, evidence, method="inference" if geometry in {"not_applicable", "unresolved"} else "original_extraction",
        derivation="Status assigned after inspection; terminal allenes are not treated as E/Z elements." if geometry == "not_applicable" else "Geometry could not be fully assessed from the encoded representation." if geometry == "unresolved" else None)
    m["stereochemistry"]["representation_status"] = f(
        representation, representation, evidence, method="inference",
        derivation="Compares the source stereochemical information with the stored candidate representation.")
    if assignment is not None:
        m["stereochemistry"]["assignment"] = f(assignment, assignment, evidence)
    if mixture is not None:
        m["stereochemistry"]["mixture"] = f(mixture, mixture, evidence)
    if note:
        m["identity_note"] = f(note, note, evidence)
    molecules[mid] = m
    return mid


def ensure_named(label, page=3, name=None, **kwargs):
    mid = "m_" + label.lower().replace(" ", "_").replace("/", "_").replace("·", "_").replace("-", "_")
    return add_molecule(mid, label, page, name=name or label, **kwargs)


# Source-derived candidate structures for the six allenyl carbonate substrates.
carbonate_tail = "COC(=O)OCc1ccccc1"
allene_smiles = {
    "1a": "C=C=C(CCCC)" + carbonate_tail,
    "1b": "C=C=C(C)" + carbonate_tail,
    "1c": "C=C=C(CCCCCCCC)" + carbonate_tail,
    "1d": "C=C=C([Si](C)(C)C)" + carbonate_tail,
    "1e": "C=C=C(c1ccccc1)" + carbonate_tail,
    "1f": "C=C=C" + carbonate_tail,
}
for label, smiles in allene_smiles.items():
    add_molecule(f"m_{label}", label, 6 if label != "1a" else 3,
                 name=f"benzyl substituted buta-2,3-dienyl carbonate {label}",
                 smiles=smiles, chirality="achiral", geometry="not_applicable",
                 representation="partial",
                 assignment="Terminal CH2 allene shown; no axial stereogenicity assigned.")

pyrazolone_smiles = "CC1=NN(c2ccccc2)C(=O)C1Cc1ccccc1"
add_molecule("m_2a", "2a", 3, "racemic pyrazol-5-one 2a", smiles=pyrazolone_smiles,
             chirality="racemic", geometry="unresolved", representation="partial",
             mixture="racemic pyrazol-5-one", note="Candidate connectivity is source-derived; enantiomers are not separately encoded.")

pyrazolone_rows = {
    "2b": "Me/o-MeC6H5CH2/Ph", "2c": "Me/m-MeC6H5CH2/Ph", "2d": "Me/p-MeC6H5CH2/Ph",
    "2e": "Me/p-ClC6H5CH2/Ph", "2f": "Me/p-BrC6H5CH2/Ph", "2g": "Me/p-FC6H5CH2/Ph",
    "2h": "Me/p-NCC6H5CH2/Ph", "2i": "Me/p-CF3C6H5CH2/Ph", "2j": "Me/m-MeOC6H5CH2/Ph",
    "2k": "Me/m-O2NC6H5CH2/Ph", "2l": "Me/Et/Ph", "2m": "Me/allyl/Ph",
    "2n": "Me/alpha-NaphthylCH2/Ph", "2o": "Me/CH2CH2OH/Ph", "2p": "Me/CH2CH2OTBS/Ph",
    "2t": "Me/Ph/Ph", "2q": "Et/Bn/Ph", "2r": "iPr/Bn/Ph", "2s": "Ph/Bn/Ph",
    "2u": "Me/Bn/Cy", "2v": "Me/Bn/Bn",
}
for label, substituents in pyrazolone_rows.items():
    add_molecule(f"m_{label}", label, 6, f"racemic pyrazol-5-one {label} ({substituents})",
                 chirality="racemic", geometry="unresolved", representation="not_encoded",
                 mixture="racemic pyrazol-5-one",
                 note="Complete candidate SMILES was withheld because the source table uses substituent shorthand and the stereochemical graph was not independently checked.")


def product3(label, config, page=6, xray=False):
    suffix = f"_{config}" if label == "3am" else ""
    mid = f"m_{label}{suffix}"
    assignment = f"({config})-{label}"
    if xray:
        assignment += "; absolute configuration determined by single-crystal X-ray diffraction"
    return add_molecule(mid, f"({config})-{label}", page, f"C-allenylation product ({config})-{label}",
                        chirality="specified", geometry="not_applicable", representation="not_encoded",
                        assignment=assignment,
                        note="The stereogenic pyrazolone carbon is source-labelled, but a full isomeric SMILES was not asserted.",
                        extra_evidence=(ev(4, "text"),) if xray else ())


m3aa = product3("3aa", "S", page=3)
m4aa = add_molecule("m_4aa", "4aa", 3, "N-allenylation side product 4aa",
                    chirality="unspecified", geometry="not_applicable", representation="not_encoded",
                    assignment="N-allenylation constitutional isomer; no stereochemical descriptor assigned.")

# Common catalysts, ligands, solvents, and characterization identities.
common = {}
for label, page in [
    ("Pd2(dba)3·CHCl3", 3), ("CHCl3", 3), ("THF", 3), ("CH3CN", 3), ("toluene", 3),
    ("DMF", 3), ("CH3OH", 3), ("AcOH", 5), ("mesitylene", 3),
]:
    common[label] = ensure_named(label, page)
for i in range(1, 12):
    common[f"L{i}"] = add_molecule(f"m_L{i}", f"L{i}", 3, f"chiral ligand L{i}",
                                    chirality="partially_specified", geometry="unresolved", representation="not_encoded",
                                    assignment="Chiral ligand drawing includes stereochemical marks; complete configuration/graph not encoded.")
for label in ["DACH-ZYC-Phos-C1", "DACH-ZYC-Phos-C2", "DACH-ZYC-Phos-C3", "DADPE-ZYC-Phos-C2"]:
    common[label] = add_molecule("m_" + label.replace("-", "_"), label, 4 if label.endswith("C1") else 5,
                                 label, chirality="specified", geometry="unresolved", representation="not_encoded",
                                 assignment="Ligand drawing derives from (1R,2R)-DACH or (1R,2R)-DADPE as labelled in Fig. 3.")
common["ent-DACH-ZYC-Phos-C1"] = add_molecule(
    "m_ent_DACH_ZYC_Phos_C1", "ent-DACH-ZYC-Phos-C1", 6, "enantiomeric DACH-ZYC-Phos-C1",
    chirality="specified", geometry="unresolved", representation="not_encoded",
    assignment="Enantiomeric ligand explicitly labelled ent-DACH-ZYC-Phos-C1.")


def set_charge(participant, key, raw, normalized, evidence, reference=None, method="original_extraction"):
    participant["charges"][key] = f(raw, normalized, evidence, method=method,
                                     derivation="Inherited from the figure/table general conditions." if method == "inherited_general_condition" else None)
    if key in {"equivalents", "mol_percent"}:
        participant["equivalence_reference"] = f(reference, reference, evidence,
                                                  method="inference",
                                                  derivation="The table defines loadings relative to allene substrate 1.")


def participant(pid, mid, role, evidence, charges=None, product_role=None):
    p = st.new_participant(pid, mid, role, evidence)
    if product_role:
        p["product_role"] = f(product_role, product_role, evidence)
    for key, data in (charges or {}).items():
        raw, normalized, reference, method = data
        set_charge(p, key, raw, normalized, evidence, reference, method)
    return p


def set_classification(r, labels, evidence, coverage="partial", limitation=None):
    if coverage in {"unresolved", "not_applicable"}:
        r["classification"]["coverage"] = f(coverage, coverage, evidence, method="inference",
                                               derivation="No observed bond changes are asserted for this record." if coverage == "not_applicable" else "Insufficient structural detail for a defensible bond-change inventory.")
        r["classification"]["labels"] = f(None if coverage == "unresolved" else [], None if coverage == "unresolved" else [], evidence,
                                             method="inference",
                                             derivation="No taxonomy labels asserted." if coverage == "unresolved" else "No product formation was observed.")
        r["classification"]["bond_changes"] = []
        r["classification"]["limitations"] = [limitation or "No complete atom mapping was available."]
        return
    changes = []
    taxonomy = {"C-C": 2, "C-O": 3, "C-N": 5, "C-P": 6}
    atoms_for = {"C-C": ("C", "C"), "C-O": ("C", "O"), "C-N": ("C", "N"), "C-P": ("C", "P")}
    label_rows = []
    for index, code in enumerate(labels, 1):
        change_id = f"b{index}"
        a, b = atoms_for[code]
        change = {"change_id": change_id, "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
                  "before_order": 0, "after_order": 1, "change_type": "formed"}
        changes.append(f("source-drawing bond change", change, evidence, method="structure_derivation",
                         derivation="Conservative bond-type assignment from the inspected reaction drawing; atoms are not fully mapped."))
        existing = next((x for x in label_rows if x["code"] == code), None)
        if existing:
            existing["supporting_change_ids"].append(change_id)
        else:
            label_rows.append({"category_number": taxonomy[code], "code": code, "supporting_change_ids": [change_id]})
    r["classification"]["bond_changes"] = changes
    r["classification"]["labels"] = f("source-drawing bond types", label_rows, evidence,
                                         method="structure_derivation",
                                         derivation="Labels include only conservatively recognized newly formed bonds; no complete mapping is claimed.")
    r["classification"]["coverage"] = f("partial", "partial", evidence, method="inference",
                                           derivation="Primary formed-bond types only; broken bonds and implicit-H changes remain unresolved.")
    r["classification"]["limitations"] = [limitation or "No reviewed atom mapping; implicit-H and leaving-group changes are unresolved."]


def make_reaction(rid, page, record_type, categories, name, transformation, participants,
                  temperature=None, time=None, solvent=None, solvent_volume=None,
                  concentration=None, scale=None, outcomes=None, labels=None,
                  evidence_extra=(), unresolved=(), procedure=None, location=None,
                  classification_coverage="partial", notes=()):
    evidence = [ev(page, "image"), ev(page, "text"), *evidence_extra]
    evidence = list(dict.fromkeys(evidence))
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, evidence)
    r["location"] = f(location or f"PDF page {page}", {"page": page, "region": location or "source page"}, evidence)
    r["reaction_name"] = f(name, name, evidence)
    r["transformation_type"] = f(transformation, transformation, evidence, method="inference",
                                    derivation="Normalized transformation description from source terminology and drawing.")
    r["participants"] = participants
    if scale:
        r["reaction_scale"] = f(scale[0], {"reference_molecule_id": scale[1], "amount": scale[2]}, evidence)
    if solvent:
        r["solvent_conditions"]["composition"] = f(solvent, [solvent], evidence)
    if solvent_volume:
        raw, normalized, status = solvent_volume
        r["solvent_conditions"]["total_volume"] = f(raw, normalized, evidence, status=status)
    if concentration:
        raw, normalized, method = concentration
        r["solvent_conditions"]["concentration"] = f(raw, normalized, evidence, method=method,
                                                        derivation="Inherited from the explicitly stated optimal/general condition." if method == "inherited_general_condition" else None)
        r["solvent_conditions"]["concentration_basis"] = f("based on allene substrate 1", "molecule 1", evidence,
                                                              method="inference", derivation="The stated scale and reaction concentration refer to allene 1.")
    step = st.new_step("reaction", evidence)
    step["participant_ids"] = [p["participant_id"] for p in participants if p["stage"] == "reaction"]
    step["operation"] = f(procedure or name, procedure or name, evidence)
    if temperature is not None:
        step["temperature"] = f(temperature[0], temperature[1], evidence, method=temperature[2] if len(temperature) > 2 else "original_extraction",
                                  derivation="Inherited from the table general condition." if len(temperature) > 2 and temperature[2] == "inherited_general_condition" else None)
    if time is not None:
        step["time"] = f(time[0], time[1], evidence, method=time[2] if len(time) > 2 else "original_extraction",
                          derivation="Inherited from the figure/table general condition." if len(time) > 2 and time[2] == "inherited_general_condition" else None)
    r["steps"] = [step]
    r["procedure"]["original_text"] = f(procedure, procedure, evidence) if procedure else missing()
    r["outcomes"] = outcomes or []
    set_classification(r, labels or [], evidence, coverage=classification_coverage)
    r["notes"] = [f(text, text, evidence) for text in notes]
    r["general_scheme_ids"] = []
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
                   "unresolved": list(unresolved) or ["Machine extraction requires expert chemistry review."]}
    r["completeness"] = "partial"
    bundle["reactions"].append(r)
    return r


def outcome(oid, kind, target, raw, normalized, evidence, basis_raw=None, basis_norm=None,
            product_role="main product", yield_type=None, scope=None, status=None):
    o = st.new_outcome(oid, kind, [target] if target else [], evidence)
    o["product_role"] = f(product_role, product_role, evidence)
    o["value"] = f(raw, normalized, evidence, status=status)
    if yield_type:
        o["yield_type"] = f(yield_type, yield_type, evidence)
    if basis_raw is not None:
        o["measurement_basis"] = f(basis_raw, basis_norm or basis_raw, evidence)
    if scope is not None:
        o["scope"] = f(scope, scope, evidence)
    return o


# Page 2: two cited literature examples and the paper's generic current-work scheme.
bg_allene = add_molecule("m_bg_allene", "generic 2,3-allenyl acetate", 2, chirality="unspecified", geometry="unresolved", representation="not_encoded")
bg_pyrazolone = add_molecule("m_bg_pyrazolone", "generic pyrazol-5-one", 2, chirality="unspecified", geometry="unresolved", representation="not_encoded")
bg_spiro1 = add_molecule("m_bg_spiro1", "literature spiro product (Pd route)", 2, chirality="partially_specified", geometry="unresolved", representation="not_encoded", assignment="Source drawing shows stereochemistry; full graph not encoded.")
bg_spiro2 = add_molecule("m_bg_spiro2", "literature spiro product (phosphine route)", 2, chirality="partially_specified", geometry="unresolved", representation="not_encoded", assignment="Source drawing shows stereochemistry; full graph not encoded.")
for rid, product, catalyst_name, base_name in [
    ("fig1b_literature_pd", bg_spiro1, "[Pd(allyl)Cl]2/Xylyl-PhanePhos", "K2CO3"),
    ("fig1b_literature_phosphine", bg_spiro2, "chiral phosphine", "Cs2CO3"),
]:
    e = [ev(2, "image"), ev(2, "text")]
    parts = [participant("allene", bg_allene, "reactant", e), participant("pyrazolone", bg_pyrazolone, "reactant", e),
             participant("product", product, "product", e, product_role="literature product"),
             participant("catalyst", ensure_named(catalyst_name, 2), "catalyst", e),
             participant("base", ensure_named(base_name, 2), "reagent", e)]
    make_reaction(rid, 2, "literature_example", [], "asymmetric spirocyclization of pyrazolones with 2,3-allenyl acetates",
                  "literature spirocyclization", parts, labels=[], classification_coverage="unresolved",
                  unresolved=["Cited prior work; not a current-paper experiment.", "Exact structures and bond inventory not encoded."],
                  location="Fig. 1b")

generic_allene = add_molecule("m_generic_1", "allene 1 generic family", 2, chirality="unspecified", geometry="not_applicable", representation="not_encoded")
generic_pyrazolone = add_molecule("m_generic_2", "pyrazolone 2 generic family", 2, chirality="racemic", geometry="unresolved", representation="not_encoded", mixture="racemic pyrazolone family")
generic_product = add_molecule("m_generic_3", "C-allenylation product 3 generic family", 2, chirality="specified", geometry="not_applicable", representation="not_encoded", assignment="Enantioenriched product family; individual R/S assignments vary by entry.")
generic_pd = ensure_named("[Pd] generic catalyst", 2)
generic_stretchable_ligand = add_molecule("m_generic_stretchable_ligand", "a stretchable chiral ligand", 2,
                                          chirality="specified", geometry="unresolved", representation="not_encoded",
                                          assignment="Generic chiral ligand class; individual ligand structure is not fixed in Fig. 1c.")
make_reaction("fig1c_generic_allenylation", 2, "generic_scheme", [], "asymmetric allenylation of pyrazolones",
              "C-allenylation", [participant("allene1", generic_allene, "reactant", [ev(2, "image")]),
              participant("pyrazolone2", generic_pyrazolone, "reactant", [ev(2, "image")]),
              participant("product3", generic_product, "product", [ev(2, "image")], product_role="generic main product"),
              participant("pd", generic_pd, "catalyst", [ev(2, "image")]),
              participant("ligand", generic_stretchable_ligand, "ligand", [ev(2, "image")])], labels=["C-C"],
              unresolved=["Generic R-group scheme; no single concrete experiment or product structure."], location="Fig. 1c")
generic_scheme_id = f"{DOC}:fig1c_generic_allenylation"


def optimization_record(rid, page, ligand, ligand_loading, solvent, volume, temp, hours,
                        y3, ee3, y4, recovery, scale_mmol=0.2, concentration=None,
                        extra_evidence=(), note=None, volume_unit_known=True):
    evidence = [ev(page, "image"), ev(page, "text"), *extra_evidence]
    parts = [
        participant("allene_1a", "m_1a", "reactant", evidence,
                    {"amount": (f"1a ({scale_mmol} mmol)", q(scale_mmol, "mmol"), None, "original_extraction")}),
        participant("pyrazolone_2a", "m_2a", "reactant", evidence,
                    {"equivalents": ("2a (1.2 equiv)", q(1.2, "equiv"), "m_1a", "original_extraction")}),
        participant("pd_catalyst", common["Pd2(dba)3·CHCl3"], "catalyst", evidence,
                    {"mol_percent": ("2.5 mol%", q(2.5, "mol%"), "m_1a", "original_extraction")}),
        participant("ligand", common[ligand], "ligand", evidence,
                    {"mol_percent": (f"{ligand_loading} mol%", q(ligand_loading, "mol%"), "m_1a", "original_extraction")}) if ligand_loading is not None else participant("ligand", common[ligand], "ligand", evidence),
        participant("solvent", common[solvent], "solvent", evidence,
                    {"volume": (f"{volume} mL", q(volume, "mL"), None, "original_extraction")}) if volume_unit_known else participant("solvent", common[solvent], "solvent", evidence),
        participant("product_3aa", m3aa, "product", evidence, product_role="C-allenylation product"),
        participant("product_4aa", m4aa, "product", evidence, product_role="N-allenylation side product"),
    ]
    outs = []
    if y3 == "trace":
        outs.append(outcome("3aa_trace", "trace", m3aa, "trace", "trace", evidence, product_role="C-allenylation product"))
    else:
        outs.append(outcome("3aa_yield", "yield", m3aa, f"{y3}%", q(y3, "%"), evidence,
                            "1H NMR analysis of crude using mesitylene as internal standard",
                            "1H NMR of crude; mesitylene internal standard", "C-allenylation product", "NMR yield"))
    if ee3 is not None:
        if ee3 == "n.d.":
            ee_norm = "not determined"
        elif isinstance(ee3, (int, float)) and ee3 < 0:
            ee_norm = {"number": abs(ee3), "unit": "%", "direction": "negative sign retained relative to (S)-3aa", "raw_text": f"{ee3}% ee"}
        else:
            ee_norm = q(ee3, "%")
        outs.append(outcome("3aa_ee", "selectivity", m3aa, f"{ee3}% ee" if ee3 != "n.d." else "n.d.", ee_norm, evidence,
                            "chiral HPLC of isolated (S)-3aa" if ee3 != "n.d." else None,
                            "chiral HPLC of isolated product" if ee3 != "n.d." else None,
                            "C-allenylation product", scope={"metric": "ee", "major_isomer": "(S)-3aa"}))
    if y4 == "trace":
        outs.append(outcome("4aa_trace", "trace", m4aa, "trace", "trace", evidence, product_role="N-allenylation side product"))
    else:
        outs.append(outcome("4aa_yield", "yield", m4aa, f"{y4}%", q(y4, "%"), evidence,
                            "1H NMR analysis of crude using mesitylene as internal standard",
                            "1H NMR of crude; mesitylene internal standard", "N-allenylation side product", "NMR yield"))
    outs.append(outcome("1a_recovery", "recovery", "m_1a", f"{recovery}%", q(recovery, "%"), evidence,
                        "1H NMR analysis of crude using mesitylene as internal standard",
                        "1H NMR of crude; mesitylene internal standard", "recovered starting material"))
    observed_labels = []
    if y3 == "trace" or (isinstance(y3, (int, float)) and y3 > 0):
        observed_labels.append("C-C")
    if y4 == "trace" or (isinstance(y4, (int, float)) and y4 > 0):
        observed_labels.append("C-N")
    volume_data = (f"{solvent} ({volume})", q(volume, "mL"), "pending_review") if volume_unit_known else (f"{solvent} ({volume})", None, "unreadable")
    unresolved_items = (["Full atom mapping, leaving-group changes, and implicit-H changes are unresolved."]
                        if observed_labels else ["No C- or N-allenylation product was observed; unaccounted mass/decomposition is not specified."])
    if ligand_loading is None:
        unresolved_items.append("The Table 1 footnote gives ligand loadings only for entries 1-13; entries 14-15 remain unreported.")
    if not volume_unit_known:
        unresolved_items.append("The Table 1 solvent header is '(v)' and does not explicitly state the unit of the parenthetical value.")
    rx = make_reaction(rid, page, "specific_experiment", ["optimization"], "asymmetric allenylation of pyrazolone 2a",
                       "competitive C- and N-allenylation", parts,
                       temperature=(f"{temp} °C", q(temp, "°C")), time=(f"{hours} h", q(hours, "h")),
                       solvent=solvent, solvent_volume=volume_data,
                       concentration=concentration, scale=(f"1a ({scale_mmol} mmol)", "m_1a", q(scale_mmol, "mmol")),
                       outcomes=outs, labels=observed_labels,
                       classification_coverage="partial" if observed_labels else "not_applicable",
                       evidence_extra=extra_evidence,
                       unresolved=unresolved_items,
                       location="Fig. 2 ligand screen" if rid.startswith("fig2_ligand") else "Fig. 2 optimization table",
                       notes=[note] if note else [])
    rx["general_scheme_ids"] = [generic_scheme_id]
    return rx


# Page 3 ligand screen L1--L11.
ligand_rows = [
    ("L1", 33, 4, 42, 0), ("L2", 0, "n.d.", 0, 85), ("L3", 0, "n.d.", 0, 95),
    ("L4", 0, "n.d.", 0, 85), ("L5", 5, "n.d.", 0, 85), ("L6", 0, "n.d.", 0, 91),
    ("L7", 0, "n.d.", 0, 93), ("L8", 36, -11, 32, 21), ("L9", 0, "n.d.", 0, 90),
    ("L10", 73, 50, 9, 13), ("L11", 84, 75, 5, 6),
]
for ligand, y3, ee3, y4, rec in ligand_rows:
    loading = 12 if ligand in {"L1", "L2", "L3", "L4"} else 6
    hours = 4 if ligand == "L11" else 12
    optimization_record(f"fig2_ligand_{ligand.lower()}", 3, ligand, loading, "CHCl3", 3, 30, hours,
                        y3, ee3, y4, rec, extra_evidence=(ev(2, "text"),) if ligand in {"L1", "L8", "L9", "L10", "L11"} else ())

fig2_table = [
    (1, 30, "THF", 3, 79, 82, 2, 18), (2, 30, "CH3CN", 3, 0, "n.d.", 0, 89),
    (3, 30, "toluene", 3, 85, 84, 2, 10), (4, 30, "DMF", 3, 0, "n.d.", 0, 85),
    (5, 30, "CH3OH", 3, 0, "n.d.", 0, 88), (6, 40, "toluene", 3, 89, 82, 2, 5),
    (7, 50, "toluene", 3, 93, 79, 3, 4), (8, 60, "toluene", 3, 94, 78, 2, 3),
    (9, 30, "toluene", 5, 84, 84, 0, 14), (10, 30, "toluene", 10, 75, 86, 0, 25),
    (11, 60, "toluene", 10, 92, 85, 0, 0),
]
for entry, temp, solvent, volume, y3, ee3, y4, rec in fig2_table:
    optimization_record(f"fig2_entry_{entry:02d}", 3, "L11", 5 if entry == 11 else 6, solvent, volume,
                        temp, 3 if entry == 11 else 12, y3, ee3, y4, rec,
                        concentration=("0.02 M", q(0.02, "M"), "original_extraction") if entry == 11 else None,
                        extra_evidence=(ev(2, "text"),))


# Page 4 Fig. 3 ligand-synthesis arrows.
for label in ["Ph2PH", "KOH", "DMSO", "NHS", "EDCI", "DCM", "DACH", "DADPE", "DMAP",
              "ethyl chloroalkanoate n=2", "ethyl chloroalkanoate n=3", "chloroacetic acid",
              "NHS ester C2", "NHS ester C3", "diphenylphosphinoacetic acid"]:
    common.setdefault(label, ensure_named(label, 4))


def ligand_synthesis(rid, reactant_labels, product_mid, conditions, labels, note):
    evidence = [ev(4, "image"), ev(4, "text")]
    parts = []
    reference_mid = common.get(reactant_labels[0][0]) or ensure_named(reactant_labels[0][0], 4)
    for index, (label, role, equiv) in enumerate(reactant_labels):
        mid = common.get(label) or ensure_named(label, 4)
        charges = {}
        if equiv is not None:
            charges["equivalents"] = (f"{equiv} equiv", q(equiv, "equiv"), reference_mid, "original_extraction")
        parts.append(participant(f"r{index+1}", mid, role, evidence, charges))
    parts.append(participant("product", product_mid, "product", evidence, product_role="ligand or ligand precursor"))
    rx = make_reaction(rid, 4, "specific_experiment", ["other"], "synthesis of ZYC-Phos ligand or precursor",
                       "ligand synthesis", parts, labels=labels, procedure=conditions,
                       location="Fig. 3b", unresolved=[note, "Reaction time, yield, workup, and purification are not reported on page 4."],
                       notes=[note])
    # Fig. 3 contains sequential thermal stages; keep both with paired temperature and no invented time.
    if "55 °C" in conditions and "rt" in conditions:
        rx["steps"] = []
        for sid, temp_raw, temp_norm in [("stage_1_55C", "55 °C", q(55, "°C")), ("stage_2_rt", "rt", "rt")]:
            s = st.new_step(sid, evidence)
            s["participant_ids"] = [p["participant_id"] for p in parts]
            s["operation"] = f(conditions, conditions, evidence)
            s["temperature"] = f(temp_raw, temp_norm, evidence)
            rx["steps"].append(s)
    return rx


for n in (2, 3):
    precursor = common[f"NHS ester C{n}"]
    ligand_synthesis(f"fig3_nhs_ester_c{n}", [("Ph2PH", "reactant", None), (f"ethyl chloroalkanoate n={n}", "reactant", 1.1),
                     ("KOH", "reagent", 2.6), ("DMSO", "solvent", None), ("NHS", "reagent", 1.1),
                     ("EDCI", "reagent", 2.0), ("DCM", "solvent", None)], precursor,
                     "1) KOH (2.6 equiv), DMSO, 55 °C; 2) NHS (1.1 equiv), EDCI (2.0 equiv), DCM, rt",
                     ["C-P", "C-O"], "Two-stage arrow; individual stage yields are not reported.")

for rid, precursor, ligand_label, diamine in [
    ("fig3_dach_c2", "NHS ester C2", "DACH-ZYC-Phos-C2", "DACH"),
    ("fig3_dach_c3", "NHS ester C3", "DACH-ZYC-Phos-C3", "DACH"),
    ("fig3_dadpe_c2", "NHS ester C2", "DADPE-ZYC-Phos-C2", "DADPE"),
]:
    ligand_synthesis(rid, [(precursor, "reactant", None), (diamine, "reactant", 0.45), ("DCM", "solvent", None)],
                     common[ligand_label], f"{diamine} (0.45 equiv), DCM, rt", ["C-N", "C-N"],
                     "Two amide bonds are shown in the symmetric product; no complete atom mapping is claimed.")

ligand_synthesis("fig3_dppaa_c1", [("Ph2PH", "reactant", None), ("chloroacetic acid", "reactant", 1.1),
                 ("KOH", "reagent", 3.5), ("DMSO", "solvent", None)], common["diphenylphosphinoacetic acid"],
                 "KOH (3.5 equiv), DMSO, 55 °C", ["C-P"], "Precursor-forming arrow for DACH-ZYC-Phos-C1.")
ligand_synthesis("fig3_dach_c1", [("diphenylphosphinoacetic acid", "reactant", None), ("DACH", "reactant", 0.48),
                 ("EDCI", "reagent", 1.05), ("DMAP", "reagent", 1.0), ("DCM", "solvent", None)],
                 common["DACH-ZYC-Phos-C1"], "DACH (0.48 equiv), EDCI (1.05 equiv), DMAP (1.0 equiv), DCM, rt",
                 ["C-N", "C-N"], "Two amide bonds are shown in the symmetric product; no complete atom mapping is claimed.")


# Page 5 Table 1 optimization entries.
table1 = [
    (1, "DACH-ZYC-Phos-C2", "CHCl3", 3, 14, 85, "trace", 67, 30, 20),
    (2, "DACH-ZYC-Phos-C2", "CHCl3", 3, 75, 87, 4, 6, 60, 12),
    (3, "DACH-ZYC-Phos-C2", "CHCl3", 3, 80, 86, 4, 0, 80, 12),
    (4, "DACH-ZYC-Phos-C2", "CHCl3", 1.5, 67, 82, 7, 4, 60, 14),
    (5, "DACH-ZYC-Phos-C2", "CHCl3", 5, 78, 87, 3, 4, 60, 14),
    (6, "DACH-ZYC-Phos-C2", "CHCl3", 10, 70, 88, 2, 12, 60, 16),
    (7, "DACH-ZYC-Phos-C2", "CHCl3", 10, 76, 89, 2, 17, 60, 12),
    (8, "DADPE-ZYC-Phos-C2", "CHCl3", 10, 79, 79, 5, 8, 60, 12),
    (9, "DACH-ZYC-Phos-C3", "CHCl3", 10, 84, 76, 7, 0, 60, 12.5),
    (10, "DACH-ZYC-Phos-C1", "CHCl3", 10, 79, 90, 2, 10, 60, 12),
    (11, "DACH-ZYC-Phos-C1", "toluene", 10, 100, 95, 0, 0, 60, 12),
    (12, "DACH-ZYC-Phos-C1", "THF", 10, 92, 83, 3, 0, 60, 12),
    (13, "DACH-ZYC-Phos-C1", "CH3CN", 10, 23, "n.d.", 4, 60, 60, 12),
    (14, "DACH-ZYC-Phos-C2", "toluene", 10, 87, 96, 0, 11, 60, 12),
    (15, "DACH-ZYC-Phos-C3", "toluene", 10, 90, 87, 0, 0, 60, 12),
]
for entry, ligand, solvent, vol, y3, ee3, y4, rec, temp, hours in table1:
    load = 6 if entry <= 6 else 5 if entry <= 13 else None
    concentration = None
    if entry in {6, 11}:
        concentration = ("0.02 M", q(0.02, "M"), "original_extraction")
    optimization_record(f"table1_entry_{entry:02d}", 5, ligand, load, solvent, vol, temp, hours,
                        y3, ee3, y4, rec, scale_mmol=0.5 if entry == 11 else 0.2,
                        concentration=concentration, extra_evidence=(ev(2, "text"), ev(4, "text")),
                        note="Table 1 header uses solvent (v); the parenthetical value is preserved without assuming a volume unit." if entry != 11 else None,
                        volume_unit_known=False)


# Page 5 prose: downstream transformations, complex preparations, and no-reaction controls.
m3am_s = product3("3am", "S", page=6)
m3da = product3("3da", "S", page=6)
m5am = add_molecule("m_5am", "(3aS,5S)-5am", 5, "tricyclic Pauson-Khand product (3aS,5S)-5am",
                    chirality="specified", geometry="unresolved", representation="not_encoded", assignment="(3aS,5S)-5am")
m6da = add_molecule("m_6da", "(S,Z)-6da", 5, "seven-membered heterocycle (S,Z)-6da",
                    chirality="specified", geometry="specified", representation="not_encoded", assignment="(S,Z)-6da")
m7da = add_molecule("m_7da", "(S)-7da", 5, "2-propynyl product (S)-7da",
                    chirality="specified", geometry="not_applicable", representation="not_encoded", assignment="(S)-7da")
ophenol = ensure_named("o-(propen-2-yl)phenol", 5)
tic = ensure_named("TiCl4", 5)
for rid, reactants, product_mid, rxn_name, outs, note in [
    ("prose_3am_to_5am", [(m3am_s, "reactant")], m5am, "allenic Pauson-Khand annulation",
     [outcome("yield", "yield", m5am, "55%", q(55, "%"), [ev(5, "text")]),
      outcome("ee", "selectivity", m5am, "99% ee", q(99, "%"), [ev(5, "text")], scope={"metric": "ee", "major_isomer": "(3aS,5S)-5am"})],
     "Catalyst and operational conditions are deferred to Fig. 4/SI and are not present on page 5."),
    ("prose_3da_to_6da", [(m3da, "reactant"), (ophenol, "reactant")], m6da, "annulation with o-(propen-2-yl)phenol",
     [outcome("yield", "yield", m6da, "72%", q(72, "%"), [ev(5, "text")]),
      outcome("ee", "selectivity", m6da, "95% ee", q(95, "%"), [ev(5, "text")], scope={"metric": "ee", "major_isomer": "(S,Z)-6da"}),
      outcome("z_e", "selectivity", m6da, ">20/1 Z/E", {"ratio": [20, 1], "qualifier": "gt", "order": "Z/E"}, [ev(5, "text")], scope={"metric": "Z/E", "order": "Z/E"})],
     "Catalyst and operational conditions are not stated in the page-5 prose."),
    ("prose_3da_to_7da", [(m3da, "reactant"), (tic, "reagent")], m7da, "allene-to-2-propynyl isomerization with TiCl4", [],
     "Page 5 reports the transformation but no yield, dose, temperature, or time."),
]:
    e = [ev(5, "image"), ev(5, "text")]
    parts = [participant(f"r{i+1}", mid, role, e) for i, (mid, role) in enumerate(reactants)]
    parts.append(participant("product", product_mid, "product", e, product_role="main product"))
    make_reaction(rid, 5, "specific_experiment", ["derivatization"], rxn_name, rxn_name, parts,
                  outcomes=outs, classification_coverage="unresolved", labels=[],
                  location="page 5 synthetic-applications prose", unresolved=[note, "Complete structures/bond changes require the cited later figure or SI."], notes=[note])

pdcl2 = ensure_named("PdCl2", 5)
pdoac2 = ensure_named("Pd(OAc)2", 5)
base_excess = ensure_named("unspecified excess bases", 5)
complex3 = add_molecule("m_complex_III", "complex III", 5, "Pd(II)-DACH-ZYC-Phos-C1 complex III",
                        chirality="partially_specified", geometry="unresolved", representation="not_encoded",
                        assignment="Chiral DACH ligand complex; exact coordination graph not encoded.")
complex4 = add_molecule("m_complex_IV", "complex IV", 5, "Pd(II)-DACH-Phenyl Trost ligand complex IV",
                        chirality="partially_specified", geometry="unresolved", representation="not_encoded",
                        assignment="Chiral DACH ligand complex; exact coordination graph not encoded.")
for rid, metal, ligand_mid, solvent, product_mid in [
    ("prose_complex_III_prep", pdcl2, common["DACH-ZYC-Phos-C1"], "toluene", complex3),
    ("prose_complex_IV_prep", pdoac2, common["L11"], "THF", complex4),
]:
    e = [ev(5, "image"), ev(5, "text")]
    parts = [participant("metal", metal, "reactant", e), participant("ligand", ligand_mid, "reactant", e),
             participant("bases", base_excess, "reagent", e), participant("solvent", common[solvent], "solvent", e),
             participant("complex", product_mid, "product", e, product_role="isolated crystalline complex"),
             participant("recrystallization_chcl3", common["CHCl3"], "purification_agent", e),
             participant("recrystallization_hexane", ensure_named("n-hexane", 5), "purification_agent", e)]
    make_reaction(rid, 5, "specific_experiment", ["mechanistic_study"], "preparation of Pd(II)-ligand complex",
                  "coordination-complex preparation", parts, solvent=solvent, labels=[], classification_coverage="unresolved",
                  procedure="reaction with excess bases followed by recrystallization from CHCl3/n-hexane",
                  location="page 5 complex-preparation prose",
                  unresolved=["Base identity/amount, time, temperature, yield, and coordination-bond classification are not reported on page 5."])

for suffix, acetic in [("standard", False), ("acoh", True)]:
    e = [ev(5, "image"), ev(5, "text")]
    parts = [participant("allene", "m_1a", "reactant", e), participant("pyrazolone", "m_2a", "reactant", e),
             participant("complex", complex3, "catalyst", e), participant("solvent", common["toluene"], "solvent", e)]
    if acetic:
        parts.append(participant("acoh", common["AcOH"], "additive", e))
    no_rxn = outcome("no_reaction", "no_reaction", None, "failed to catalyze", True, e, product_role="no observed allenylation product")
    make_reaction(f"prose_complex_III_control_{suffix}", 5, "specific_experiment", ["control", "mechanistic_study"],
                  "allenylation control with isolated complex III", "no observed reaction", parts,
                  temperature=("60 °C (inherited standard conditions)", q(60, "°C"), "inherited_general_condition"),
                  solvent="toluene", concentration=("0.02 M (inherited standard conditions)", q(0.02, "M"), "inherited_general_condition"),
                  outcomes=[no_rxn], labels=[], classification_coverage="not_applicable",
                  procedure="under the standard conditions" + (" with AcOH as protic additive" if acetic else ""),
                  location="page 5 mechanistic-control prose",
                  unresolved=[("AcOH amount and reaction time are not stated on page 5; no product or detection limit is reported."
                               if acetic else "Reaction time, product identity/detection limit, and any recovered material are not stated on page 5.")])


# Page 6 Table 2 substrate scope.
scope_rows = [
    (1, "1a", "2a", "3aa", "S", 96, 95), (2, "1b", "2a", "3ba", "S", 91, 93),
    (3, "1c", "2a", "3ca", "S", 90, 93), (4, "1d", "2a", "3da", "S", 96, 94),
    (5, "1e", "2a", "3ea", "S", ">99", 97), (6, "1f", "2a", "3fa", "R", 86, 78),
    (7, "1a", "2b", "3ab", "R", 97, 95), (8, "1a", "2c", "3ac", "S", 96, 92),
    (9, "1a", "2d", "3ad", "S", 92, 95), (10, "1a", "2e", "3ae", "S", 97, 95),
    (11, "1a", "2f", "3af", "S", 96, 94), (12, "1a", "2g", "3ag", "S", 94, 94),
    (13, "1a", "2h", "3ah", "S", 97, 94), (14, "1a", "2i", "3ai", "S", 96, 94),
    (15, "1a", "2j", "3aj", "S", 95, 94), (16, "1a", "2k", "3ak", "S", 86, 98),
    (17, "1a", "2l", "3al", "S", 87, 94), (18, "1a", "2m", "3am", "S", 96, 93),
    (19, "1a", "2m", "3am", "R", 95, 94), (20, "1a", "2n", "3an", "R", 95, 97),
    (21, "1a", "2o", "3ao", "R", 94, 80), (22, "1a", "2p", "3ap", "R", 91, 93),
    (23, "1a", "2t", "3at", "R", 62, 36), (24, "1a", "2q", "3aq", "S", ">99", 96),
    (25, "1a", "2r", "3ar", "S", 91, 97), (26, "1a", "2s", "3as", "S", 94, 95),
    (27, "1b", "2u", "3bu", "S", 84, 95), (28, "1b", "2v", "3bv", "S", 82, 90),
]
for entry, allene, pyrazolone, product_label, config, yield_value, ee_value in scope_rows:
    product_mid = product3(product_label, config, page=6, xray=(product_label == "3ak"))
    evidence = [ev(6, "image"), ev(6, "text"), ev(4, "text")]
    catalyst_loading = 3.0 if entry in {23, 27, 28} else 2.5
    ligand_label = "ent-DACH-ZYC-Phos-C1" if entry == 19 else "DACH-ZYC-Phos-C1"
    ligand_loading = 6.0 if entry in {23, 27, 28} else 5.0
    pyraz_eq = 1.4 if entry == 17 else 1.2
    temp = 80 if entry in {23, 27, 28} else 60
    hours = 2 if entry == 6 else 13 if entry == 17 else 11 if entry in {18, 19} else 12
    parts = [
        participant("allene", f"m_{allene}", "reactant", evidence,
                    {"amount": (f"{allene} (0.5 mmol)", q(0.5, "mmol"), None, "original_extraction")}),
        participant("pyrazolone", f"m_{pyrazolone}", "reactant", evidence,
                    {"equivalents": (f"{pyrazolone} ({pyraz_eq} equiv)", q(pyraz_eq, "equiv"), f"m_{allene}", "original_extraction")}),
        participant("pd_catalyst", common["Pd2(dba)3·CHCl3"], "catalyst", evidence,
                    {"mol_percent": (f"{catalyst_loading} mol%", q(catalyst_loading, "mol%"), f"m_{allene}", "original_extraction")}),
        participant("ligand", common[ligand_label], "ligand", evidence,
                    {"mol_percent": (f"{ligand_loading} mol%", q(ligand_loading, "mol%"), f"m_{allene}", "original_extraction")}),
        participant("solvent", common["toluene"], "solvent", evidence,
                    {"volume": ("25 mL", q(25, "mL"), None, "original_extraction")}),
        participant("product", product_mid, "product", evidence, product_role="C-allenylation product"),
    ]
    y_norm = q(99, "%", qualifier="gt", raw_text=">99%") if yield_value == ">99" else q(yield_value, "%")
    outs = [
        outcome("yield", "yield", product_mid, f"{yield_value}%", y_norm, evidence,
                product_role="C-allenylation product", yield_type="reported yield"),
        outcome("ee", "selectivity", product_mid, f"{ee_value}% ee", q(ee_value, "%"), evidence,
                "chiral HPLC", "chiral HPLC", "C-allenylation product",
                scope={"metric": "ee", "major_isomer": f"({config})-{product_label}"}),
    ]
    if entry == 23:
        outs.append(outcome("1a_recovery", "recovery", "m_1a", "9%", q(9, "%"), evidence,
                            "1H NMR analysis of crude using mesitylene as internal standard",
                            "1H NMR of crude; mesitylene internal standard", "recovered starting material"))
    rx = make_reaction(f"table2_entry_{entry:02d}", 6, "specific_experiment", ["substrate_scope"],
                       "asymmetric allenylation of pyrazol-5-one", "C-allenylation", parts,
                       temperature=(f"{temp} °C", q(temp, "°C")), time=(f"{hours} h", q(hours, "h")),
                       solvent="toluene", solvent_volume=("25 mL", q(25, "mL"), "pending_review"),
                       concentration=("0.02 M", q(0.02, "M"), "original_extraction"),
                       scale=(f"{allene} (0.5 mmol)", f"m_{allene}", q(0.5, "mmol")),
                       outcomes=outs, labels=["C-C"], location=f"Table 2 entry {entry}",
                       unresolved=["Yield measurement basis is not stated on page 6; it is not assumed to be isolated.",
                                   "Full product SMILES, atom mapping, leaving-group changes, and implicit-H changes remain unresolved."],
                       notes=["Entry-specific exceptions from Table 2 footnotes were applied; all other conditions come from the table general condition."])
    rx["general_scheme_ids"] = [generic_scheme_id]


# Non-synthetic context from the inspected pages.
bundle["source_context"] = [
    {"context_id": f"{DOC}:page1_metadata", "document_id": DOC, "kind": "bibliographic_and_abstract",
     "evidence_ids": [ev(1, "image"), ev(1, "text")],
     "content": "Article title/authors/DOI and abstract; no concrete experiment on page 1."},
    {"context_id": f"{DOC}:fig1a_bioactive_examples", "document_id": DOC, "kind": "background",
     "evidence_ids": [ev(2, "image"), ev(2, "text")],
     "content": "Bioactive pyrazolone examples A-C are contextual structures, not current-paper syntheses."},
    {"context_id": f"{DOC}:fig3a_ligand_design", "document_id": DOC, "kind": "design_concept",
     "evidence_ids": [ev(4, "image"), ev(4, "text")],
     "content": "Conceptual comparison of Feng, Trost, and stretchable ZYC-Phos ligand pockets."},
    {"context_id": f"{DOC}:page5_mechanism_dft", "document_id": DOC, "kind": "mechanistic_and_computational_discussion",
     "evidence_ids": [ev(5, "image"), ev(5, "text"), ev(6, "image"), ev(6, "text")],
     "content": "Proposed mechanism, X-ray geometry discussion, and DFT transition-state energetics; not preparative experiments except separately recorded complex preparations and controls."},
]

bundle["molecules"] = list(molecules.values())

# Raw visual inventory and exact page coverage are generated from the same manually curated record IDs.
by_page = {page: [] for page in range(1, 7)}
for reaction in bundle["reactions"]:
    page = reaction["location"]["normalized_value"]["page"]
    if page in by_page:
        by_page[page].append(reaction["reaction_id"])

inventory = {
    "document_id": DOC,
    "source_review_mode": "manual visual inspection of original-detail rendered page images",
    "pages": {
        "1": ["title/author/DOI block", "abstract", "affiliations; no reaction scheme"],
        "2": ["two-column introduction/optimization narrative", "Fig. 1a bioactive context", "Fig. 1b two cited spirocyclization routes", "Fig. 1c generic current-work allenylation"],
        "3": ["Fig. 2 reaction drawing", "L1-L11 ligand screen", "optimization table entries 1-11", "caption and footnotes"],
        "4": ["Fig. 3a ligand design concept", "Fig. 3b seven ligand/precursor synthesis arrows", "scope narrative"],
        "5": ["Table 1 entries 1-15 and footnotes", "three synthetic-application prose experiments", "two Pd-complex preparations", "two failed complex-III controls", "mechanistic/DFT prose"],
        "6": ["Table 2 entries 1-28 and footnotes", "DFT/non-covalent-interaction prose"],
    },
    "reaction_record_ids_by_page": by_page,
    "semantic_accuracy_verified": False,
}
(RUN / "raw_model/visual_inventory.json").write_text(json.dumps(inventory, ensure_ascii=False, indent=2) + "\n")

stereo_diagnostics = {label: st.stereo_inspect(smiles) for label, smiles in allene_smiles.items()}
stereo_diagnostics["2a"] = st.stereo_inspect(pyrazolone_smiles)
(RUN / "raw_model/stereo_diagnostics.json").write_text(json.dumps(stereo_diagnostics, ensure_ascii=False, indent=2) + "\n")

coverage_regions = {
    1: [{"region": "full page: title, abstract, affiliations", "finding": "Bibliographic/abstract content only; no concrete experiment or reaction arrow.", "evidence_ids": [ev(1, "image"), ev(1, "text")]}],
    2: [
        {"region": "upper two-column prose", "finding": "Optimization narrative and current-work motivation; numeric experiments are linked to their Fig. 2/Table 1 records rather than duplicated.", "evidence_ids": [ev(2, "image"), ev(2, "text")]},
        {"region": "Fig. 1a", "finding": "Bioactive pyrazolone background structures; non-synthetic source context.", "evidence_ids": [ev(2, "image")]},
        {"region": "Fig. 1b", "finding": "Two cited prior-work spirocyclization routes, recorded as literature examples, not current-paper experiments.", "evidence_ids": [ev(2, "image"), ev(2, "text")]},
        {"region": "Fig. 1c", "finding": "Generic asymmetric C-allenylation scheme for the present study.", "evidence_ids": [ev(2, "image")]},
    ],
    3: [
        {"region": "Fig. 2 upper/middle ligand structures and results", "finding": "Eleven separate L1-L11 optimization experiments, including zero-yield and competing N-allenylation results.", "evidence_ids": [ev(3, "image"), ev(3, "text")]},
        {"region": "Fig. 2 lower table entries 1-11", "finding": "Eleven separate solvent/temperature/concentration optimization experiments.", "evidence_ids": [ev(3, "image"), ev(3, "text")]},
        {"region": "Fig. 2 caption/footnotes", "finding": "General loadings, times, NMR yield/recovery basis, and chiral-HPLC ee basis mapped to the corresponding rows.", "evidence_ids": [ev(3, "image"), ev(3, "text")]},
    ],
    4: [
        {"region": "Fig. 3a", "finding": "Ligand-design concept only; stored as source context.", "evidence_ids": [ev(4, "image"), ev(4, "text")]},
        {"region": "Fig. 3b", "finding": "Seven source-visible ligand/precursor synthesis instances across the n=2/n=3 and C1 routes; yields and times are unreported.", "evidence_ids": [ev(4, "image"), ev(4, "text")]},
        {"region": "lower-page scope prose", "finding": "Narrative cross-reference for Table 2 products and stereochemical assignments; merged as evidence with page-6 scope records.", "evidence_ids": [ev(4, "image"), ev(4, "text")]},
    ],
    5: [
        {"region": "Table 1", "finding": "Fifteen independent optimization entries with product yield, ee, N-allenylation, and starting-material recovery retained separately.", "evidence_ids": [ev(5, "image"), ev(5, "text")]},
        {"region": "left-column synthetic applications", "finding": "Three prose-reported derivatizations (5am, 6da, 7da) recorded; detailed later-figure/SI conditions remain unresolved.", "evidence_ids": [ev(5, "image"), ev(5, "text")]},
        {"region": "left-column complex preparation/control discussion", "finding": "Two Pd-complex preparations and two no-reaction complex-III controls recorded as mechanistic-study experiments.", "evidence_ids": [ev(5, "image"), ev(5, "text")]},
        {"region": "right-column DFT discussion", "finding": "Computational/mechanistic context, not a synthetic experiment.", "evidence_ids": [ev(5, "image"), ev(5, "text")]},
    ],
    6: [
        {"region": "Table 2 entries 1-28", "finding": "Twenty-eight independent substrate-scope experiments, including opposite enantiomer 3am as a distinct molecule and footnote-specific condition changes.", "evidence_ids": [ev(6, "image"), ev(6, "text")]},
        {"region": "Table 2 footnotes", "finding": "General conditions, chiral-HPLC ee method, time/load/temperature exceptions, and entry-23 recovery mapped to records.", "evidence_ids": [ev(6, "image"), ev(6, "text")]},
        {"region": "lower-page DFT/NCI prose", "finding": "Computational transition-state analysis only; no new preparative reaction.", "evidence_ids": [ev(6, "image"), ev(6, "text")]},
    ],
}
coverage = {"document_id": DOC, "semantic_accuracy_verified": False, "pages": []}
page_unresolved = {
    1: ["Abstract-level claims are not converted into concrete experiments."],
    2: ["Literature-example exact structures and conditions are only partially visible; no SI association is inferred."],
    3: ["Full structures of chiral ligands and product atom mapping remain unencoded."],
    4: ["Ligand-synthesis yields, times, workup, purification, and several precursor structures are not reported on this page."],
    5: ["Table solvent parenthetical values use header '(v)' without an explicit unit; entries 14-15 ligand loading is not stated by the footnote ranges.", "Derivatization and complex-preparation details deferred to later figures/SI remain missing."],
    6: ["Scope yield measurement basis is not stated on page 6; complete product graphs/SMILES and atom mapping remain unresolved."],
}
for page in range(1, 7):
    coverage["pages"].append({"page": page, "state": "attempted", "regions": coverage_regions[page],
                              "record_ids": by_page[page], "unresolved": page_unresolved[page]})
(RUN / "coverage/pages_0001_0006.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")

# Literal mention search is navigation-only and preserved separately.
dataset = json.loads((RUN / "dataset.json").read_text())
mentions = st.source_mentions(dataset, ["3aa", "4aa", "3am", "3ak", "DACH-ZYC-Phos-C1"], context=220)
(RUN / "raw_model/source_mentions.json").write_text(json.dumps(mentions, ensure_ascii=False, indent=2) + "\n")

st.write(RUN / "submitted.json", bundle)
print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
                  "by_page": {str(k): len(v) for k, v in by_page.items()}}, ensure_ascii=False, indent=2))
