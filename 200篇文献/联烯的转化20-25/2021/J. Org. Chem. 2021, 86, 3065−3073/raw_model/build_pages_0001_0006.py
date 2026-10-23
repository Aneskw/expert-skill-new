#!/usr/bin/env python3
"""Build source-derived v3.1 records for original PDF pages 1--6.

The script is deterministic and uses only the prepared dataset plus facts manually
transcribed from the six inspected page images.  It does not certify chemistry.
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

RUN = ROOT / "runs/current"
DOC = "doc_64f4a881f331038c"
E = {p: f"{DOC}:page{p}:text" for p in range(1, 10)}
I = {p: f"{DOC}:page{p}:image" for p in range(1, 7)}
MOLS: dict[str, dict] = {}
REACTIONS: list[dict] = []


def dump(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


_UNSET = object()


def uniq_ev(evidence):
    return tuple(dict.fromkeys(evidence))


def f(raw=None, normalized=_UNSET, evidence=(), method="original_extraction",
      status=None, derivation=None, alternatives=()):
    evidence = uniq_ev(evidence)
    if normalized is _UNSET:
        normalized = raw
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   derivation=derivation, alternatives=alternatives)


def missing(evidence=(), status="unreported"):
    return f(evidence=evidence, method="original_extraction", status=status)


def not_applicable(evidence=()):
    return f(evidence=evidence, method="not_applicable", status="not_applicable")


def derived(raw, normalized, evidence, derivation):
    return f(raw, normalized, evidence, method="structure_derivation",
             derivation=derivation)


def q(number, unit, raw_text=None, **extra):
    out = {"number": number, "unit": unit}
    if raw_text is not None:
        out["raw_text"] = raw_text
    out.update(extra)
    return out


def ensure_status_facts(m, ev, chirality="unresolved", geometry="unresolved",
                        representation="unresolved", source_raw=None):
    source_raw = source_raw or "stereochemistry reviewed from the cited name/drawing"
    for key, val in (("chirality_status", chirality),
                     ("geometry_status", geometry),
                     ("representation_status", representation)):
        m["stereochemistry"][key] = f(source_raw, val, ev)


def add_molecule(mid, label, name, ev, smiles=None, chirality="unresolved",
                 geometry="unresolved", representation=None, atom_centers=None,
                 double_bonds=None, axial=None, mixture=None, assignment=None,
                 note=None, aliases=None):
    ev = uniq_ev(ev)
    if mid in MOLS:
        return mid
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    m["names"]["original"] = f(name, name, ev)
    m["names"]["english"] = f(name, name, ev)
    m["names"]["abbreviations"] = f(label, [label], ev)
    m["names"]["aliases"] = f(aliases, aliases, ev) if aliases else missing(ev)
    m["cas_numbers"] = missing(ev)
    m["identity_note"] = f(note, note, ev) if note else missing(ev)
    m["external_properties"] = []
    m["stereochemistry"]["atom_centers"] = (
        f(atom_centers, atom_centers, ev) if atom_centers is not None else missing(ev)
    )
    m["stereochemistry"]["double_bonds"] = (
        f(double_bonds, double_bonds, ev) if double_bonds is not None else missing(ev)
    )
    m["stereochemistry"]["axial"] = (
        f(axial, axial, ev) if axial is not None else missing(ev)
    )
    m["stereochemistry"]["mixture"] = (
        f(mixture, mixture, ev) if mixture is not None else missing(ev)
    )
    m["stereochemistry"]["assignment"] = (
        f(assignment, assignment, ev) if assignment is not None else missing(ev)
    )
    if smiles:
        m["representation_type"] = "smiles"
        reason = ("Candidate graph transcribed from the cited structure/name; RDKit parsing and "
                  "stereo diagnostics do not constitute source verification.")
        m["smiles"] = derived(name, smiles, ev, reason)
        diag = st.stereo_inspect(smiles)
        m["stereochemistry"]["isomeric_smiles"] = derived(
            name, diag["canonical_isomeric_smiles"], ev, reason)
        m["stereochemistry"]["identity_key"] = derived(
            name, diag["stereo_graph_key"], ev,
            "Hash of the candidate canonical isomeric graph; not a sample identity.")
        if representation is None:
            representation = "encoded" if (chirality in {"specified", "partially_specified"}
                                              or geometry in {"specified", "partially_specified"}) else "not_applicable"
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing(ev)
        m["stereochemistry"]["isomeric_smiles"] = missing(ev)
        m["stereochemistry"]["identity_key"] = missing(ev)
        representation = representation or "not_encoded"
    ensure_status_facts(m, ev, chirality, geometry, representation)
    MOLS[mid] = m
    return mid


def add_template(mid, label, name, ev, fixed_structure, variables,
                 chirality="unresolved", geometry="unresolved", note=None):
    ev = uniq_ev(ev)
    if mid in MOLS:
        return mid
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    m["names"]["original"] = f(name, name, ev)
    m["names"]["english"] = f(name, name, ev)
    m["names"]["abbreviations"] = f(label, [label], ev)
    m["names"]["aliases"] = missing(ev)
    m["cas_numbers"] = missing(ev)
    m["representation_type"] = "generic_template"
    m["smiles"] = missing(ev)
    m["generic_template"] = {
        "fixed_structure": f(fixed_structure, fixed_structure, ev),
        "variables": [
            {"label": v, "attachment_points": f(f"drawing attachment for {v}", {"label": v}, ev),
             "definition": f(f"R-group {v} in source drawing", f"source-defined variable {v}", ev),
             "allowed_values": missing(ev), "constraints": missing(ev)}
            for v in variables
        ],
        "instance_mapping": [],
        "reason": f("Generic R-group drawing", "Source structure contains unexpanded R groups; no concrete SMILES asserted.", ev),
    }
    m["stereochemistry"]["atom_centers"] = missing(ev)
    m["stereochemistry"]["double_bonds"] = missing(ev)
    m["stereochemistry"]["axial"] = missing(ev)
    m["stereochemistry"]["mixture"] = missing(ev)
    generic_assignment = note or ("Generic drawing indicates stereochemical content but does not fix one concrete configuration."
                                  if chirality in {"specified", "partially_specified"} else None)
    m["stereochemistry"]["assignment"] = (
        f(generic_assignment, generic_assignment, ev) if generic_assignment else missing(ev)
    )
    ensure_status_facts(m, ev, chirality, geometry, "not_encoded",
                        "generic R-group drawing")
    m["stereochemistry"]["isomeric_smiles"] = missing(ev)
    m["stereochemistry"]["identity_key"] = missing(ev)
    m["identity_note"] = f(note, note, ev) if note else missing(ev)
    m["external_properties"] = []
    MOLS[mid] = m
    return mid


def add_common(mid, label, name, smiles, ev=(E[2],)):
    return add_molecule(mid, label, name, ev, smiles=smiles, chirality="achiral",
                        geometry="not_applicable", representation="not_applicable",
                        assignment={"basis": "candidate molecular graph has no assigned stereogenic element"},
                        note="Common reagent/solvent identity is a source-linked candidate; no CAS lookup performed.")


def participant(pid, mid, role, ev, stage="reaction", mass=None, amount=None,
                volume=None, equivalents=None, mol_percent=None, reference=None,
                product_role=None):
    ev = uniq_ev(ev)
    p = st.new_participant(pid, mid, role, ev, stage=stage)
    if product_role is not None:
        p["product_role"] = f(product_role, product_role, ev)
    if mass is not None:
        p["charges"]["mass"] = f(str(mass), mass, ev)
    if amount is not None:
        p["charges"]["amount"] = f(str(amount), amount, ev)
    if volume is not None:
        p["charges"]["volume"] = f(str(volume), volume, ev)
    if equivalents is not None:
        p["charges"]["equivalents"] = f(str(equivalents), equivalents, ev)
    if mol_percent is not None:
        p["charges"]["mol_percent"] = f(str(mol_percent), mol_percent, ev)
    if reference is not None:
        p["equivalence_reference"] = f(reference, reference, ev)
    return p


def step(step_id, stage, pids, ev, operation, temperature=None, time=None,
         atmosphere=None, pressure=None):
    ev = uniq_ev(ev)
    s = st.new_step(step_id, ev, stage=stage)
    s["participant_ids"] = pids
    s["operation"] = f(operation, operation, ev)
    s["temperature"] = f(str(temperature), temperature, ev) if temperature is not None else missing(ev)
    s["time"] = f(str(time), time, ev) if time is not None else missing(ev)
    s["atmosphere"] = f(atmosphere, atmosphere, ev) if atmosphere is not None else missing(ev)
    s["pressure"] = f(pressure, pressure, ev) if pressure is not None else missing(ev)
    return s


def outcome(oid, kind, targets, ev, value=None, yield_type=None,
            basis=None, scope=None, product_role="main product"):
    ev = uniq_ev(ev)
    o = st.new_outcome(oid, kind, targets, ev)
    o["product_role"] = f(product_role, product_role, ev) if targets else not_applicable(ev)
    o["value"] = f(str(value), value, ev) if value is not None else missing(ev)
    o["yield_type"] = f(yield_type, yield_type, ev) if yield_type else missing(ev)
    o["measurement_basis"] = f(basis, basis, ev) if basis else missing(ev)
    o["scope"] = f(scope, scope, ev) if scope is not None else missing(ev)
    return o


def labels_for(changes):
    labels = {}
    for c in changes:
        if c["change_type"] in {"formed", "order_increased", "hydrogen_gain"}:
            lab = st.bond_label(c["atoms"][0]["element"], c["atoms"][1]["element"], c["after_order"])
            labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(c["change_id"])
    return list(labels.values())


def classification(ev, mode="final"):
    ev = uniq_ev(ev)
    if mode == "none":
        return {"standard_id": st.STANDARD, "mapped_reaction_smiles": missing(ev),
                "bond_changes": [],
                "labels": derived("no observed transformation", [], ev, "No bond-formation label is assigned to a zero-yield/no-reaction experiment."),
                "coverage": f("no observed transformation", "not_applicable", ev, method="inference",
                              derivation="The source reports 0% yield/no reaction."),
                "limitations": ["No product graph or intended transformation is classified."], "model_votes": []}
    if mode == "intermediate":
        changes = [
            {"change_id": "b1", "atoms": [{"element": "C"}, {"element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b2", "atoms": [{"element": "C"}, {"element": "N"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b3", "atoms": [{"element": "C"}, {"element": "N"}], "before_order": 1, "after_order": 0, "change_type": "broken"},
            {"change_id": "b4", "atoms": [{"element": "C"}, {"element": "C"}], "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
        ]
    elif mode == "cyclization":
        changes = [
            {"change_id": "b1", "atoms": [{"element": "C"}, {"element": "N"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b2", "atoms": [{"element": "C"}, {"element": "N"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b3", "atoms": [{"element": "C"}, {"element": "Br"}], "before_order": 1, "after_order": 0, "change_type": "broken"},
            {"change_id": "b4", "atoms": [{"element": "C"}, {"element": "Br"}], "before_order": 1, "after_order": 0, "change_type": "broken"},
        ]
    elif mode == "wittig":
        changes = [
            {"change_id": "b1", "atoms": [{"element": "C"}, {"element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b2", "atoms": [{"element": "C"}, {"element": "C"}], "before_order": 1, "after_order": 2, "change_type": "order_increased"},
        ]
    elif mode == "amination":
        changes = [
            {"change_id": "b1", "atoms": [{"element": "C"}, {"element": "N"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b2", "atoms": [{"element": "C"}, {"element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
            {"change_id": "b3", "atoms": [{"element": "C"}, {"element": "O"}], "before_order": 2, "after_order": 0, "change_type": "broken"},
        ]
    elif mode == "methylation_ether":
        changes = [
            {"change_id": "b1", "atoms": [{"element": "C"}, {"element": "N"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b2", "atoms": [{"element": "C"}, {"element": "O"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b3", "atoms": [{"element": "C"}, {"element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 2},
            {"change_id": "b4", "atoms": [{"element": "C"}, {"element": "O"}], "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
        ]
    else:
        changes = [
            {"change_id": "b1", "atoms": [{"element": "C"}, {"element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b2", "atoms": [{"element": "C"}, {"element": "O"}], "before_order": 0, "after_order": 2, "change_type": "formed"},
            {"change_id": "b3", "atoms": [{"element": "C"}, {"element": "N"}], "before_order": 1, "after_order": 0, "change_type": "broken"},
            {"change_id": "b4", "atoms": [{"element": "C"}, {"element": "C"}], "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
        ]
    deriv = ("Partial net-bond interpretation from the source reaction drawing. Atom mapping, "
             "hydrogen/coproduct balance, and elementary mechanism are not asserted.")
    return {"standard_id": st.STANDARD, "mapped_reaction_smiles": missing(ev),
            "bond_changes": [derived("source reaction drawing", c, ev, deriv) for c in changes],
            "labels": derived("source reaction drawing", labels_for(changes), ev, deriv),
            "coverage": f("partial structural interpretation", "partial", ev, method="inference", derivation=deriv),
            "limitations": [deriv, "No three-model votes were performed."], "model_votes": []}


def reaction(rid, ev, record_type="specific_experiment", categories=None,
             location=None, name=None, transformation=None, class_mode="final"):
    ev = uniq_ev(ev)
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories or ["other"], categories or ["other"], ev)
    r["location"] = f(location or rid, location or rid, ev)
    r["reaction_name"] = f(name, name, ev) if name else missing(ev)
    r["transformation_type"] = (f(transformation, transformation, ev, method="inference",
                                    derivation="Transformation label summarized from source reaction name/drawing.")
                                   if transformation else missing(ev))
    r["classification"] = classification(ev, class_mode)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [],
                   "unresolved": ["Machine record; source structure transcription and partial bond classification require expert review."]}
    r["completeness"] = "partial"
    return r


def set_procedure(r, ev, original, general=None, modifications=None, addition_order=None,
                  pretreatment=None):
    r["procedure"]["original_text"] = f(original, original, ev)
    r["procedure"]["general_procedure_reference"] = f(general, general, ev) if general else missing(ev)
    r["procedure"]["modifications"] = f(modifications, modifications, ev) if modifications else missing(ev)
    r["procedure"]["addition_order"] = f(addition_order, addition_order, ev) if addition_order else missing(ev)
    r["procedure"]["pretreatment"] = f(pretreatment, pretreatment, ev) if pretreatment else missing(ev)


# Common source molecules and materials.
common = [
    ("m_meoh", "MeOH", "methanol", "CO"), ("m_dcm", "CH2Cl2", "dichloromethane", "ClCCl"),
    ("m_water", "H2O", "water", "O"), ("m_chcl3", "CHCl3", "chloroform", "ClC(Cl)Cl"),
    ("m_dce", "DCE", "1,2-dichloroethane", "ClCCCl"), ("m_toluene", "PhMe", "toluene", "Cc1ccccc1"),
    ("m_hexane", "hexane", "hexane", "CCCCCC"), ("m_etoac", "EtOAc", "ethyl acetate", "CCOC(C)=O"),
    ("m_thf", "THF", "tetrahydrofuran", "C1CCOC1"), ("m_dioxane", "dioxane", "1,4-dioxane", "O1CCOCC1"),
    ("m_acetone", "acetone", "acetone", "CC(C)=O"), ("m_mecn", "MeCN", "acetonitrile", "CC#N"),
    ("m_dmf", "DMF", "N,N-dimethylformamide", "CN(C)C=O"), ("m_dmso", "DMSO", "dimethyl sulfoxide", "CS(C)=O"),
    ("m_etoh", "EtOH", "ethanol", "CCO"), ("m_silica", "silica gel", "silica gel", None),
    ("m_petroleum_ether", "petroleum ether", "petroleum ether (composition not specified)", None),
    ("m_na2so4", "Na2SO4", "anhydrous sodium sulfate", "[Na+].[Na+].[O-]S(=O)(=O)[O-]"),
    ("m_brine", "brine", "brine (composition/concentration not specified)", None),
    ("m_ether", "ether", "diethyl ether", "CCOCC"),
    ("m_celite", "Celite", "Celite filter aid", None),
]
for mid, label, name, smi in common:
    if smi:
        add_common(mid, label, name, smi, (E[2], E[3], E[4], E[5]))
    else:
        add_molecule(mid, label, name, (E[2], E[3], E[4], E[5]),
                     chirality="not_applicable", geometry="not_applicable",
                     representation="not_applicable",
                     assignment={"reason": "composition or macromolecular material not represented as a unique small-molecule SMILES"},
                     note="Named material with non-unique or macromolecular composition; no concrete SMILES asserted.")


# Core substrates and products.
add_molecule("m_rac1a", "rac-1a", "racemic (E)-1-(4-phenylbut-3-en-2-yl)pyrrolidine",
             (E[1], E[2], I[2]), "CC(N1CCCC1)/C=C/c1ccccc1", "racemic", "specified", "partial",
             atom_centers=[{"center": "allylic carbon", "configuration": "R/S racemate"}],
             double_bonds=[{"bond": "allylic C=C", "configuration": "E"}],
             mixture={"type": "racemate", "ratio": "not reported"},
             assignment={"source": "rac-1a label and E drawing/name"})
add_molecule("m_1a", "1a", "(S,E)-1-(4-phenylbut-3-en-2-yl)pyrrolidine",
             (E[2], E[3], E[5], I[2], I[3], I[5]), "C[C@H](N1CCCC1)/C=C/c1ccccc1", "specified", "specified",
             atom_centers=[{"center": "allylic carbon", "configuration": "S"}],
             double_bonds=[{"bond": "allylic C=C", "configuration": "E"}],
             assignment={"source": "compound heading and Table 3; 97% ee"})
add_molecule("m_2a", "2a", "1-phenylbuta-2,3-dien-1-one",
             (E[1], E[2], E[5], I[1], I[2], I[5]), "C=C=CC(=O)c1ccccc1", "achiral", "not_applicable", "not_applicable",
             axial={"allene": "terminal allene; no axial configuration specified or possible for CH2 terminus"},
             assignment={"basis": "terminal allenone candidate graph"},
             note="Prepared according to literature; no SI was supplied in this session.")
add_molecule("m_rac3a", "rac-3a", "(2Z,6E)-1,5-diphenyl-3-(pyrrolidin-1-yl)octa-2,6-dien-1-one",
             (E[1], E[5], I[1], I[5]),
             "O=C(c1ccccc1)/C=C(\\N1CCCC1)CC(c1ccccc1)/C=C/C",
             "racemic", "specified", "partial",
             atom_centers=[{"center": "C5", "configuration": "R/S, not resolved"}],
             double_bonds=[{"bond": "C2=C3", "configuration": "Z"}, {"bond": "C6=C7", "configuration": "E"}],
             mixture={"type": "racemic", "ratio": "not reported"},
             assignment={"source": "rac-3a heading; 2Z,6E name"})
add_molecule("m_rac4a", "rac-4a", "racemic (E)-1,5-diphenyloct-6-ene-1,3-dione",
             (E[2], I[2]), "O=C(c1ccccc1)CC(=O)CC(c1ccccc1)/C=C/C", "racemic", "specified", "partial",
             atom_centers=[{"center": "C5", "configuration": "R/S racemate"}],
             double_bonds=[{"bond": "C6=C7", "configuration": "E"}],
             mixture={"type": "racemate", "ratio": "not reported"},
             assignment={"source": "rac-4a label and product drawing"},
             note="Named diketone representation; the paper later states products 4 are predominantly enol in NMR spectra.")


allyl = {
    "1a": ("(S,E)-1-(4-phenylbut-3-en-2-yl)pyrrolidine", "C[C@H](N1CCCC1)/C=C/c1ccccc1", "S", ["E"]),
    "1b": ("(S,E)-1-(1-phenylpent-1-en-3-yl)pyrrolidine", "CC[C@H](N1CCCC1)/C=C/c1ccccc1", "S", ["E"]),
    "1c": ("(S,E)-1-(1,5-diphenylpent-1-en-3-yl)pyrrolidine", "c1ccccc1CC[C@H](N2CCCC2)/C=C/c1ccccc1", "S", ["E"]),
    "1d": ("(S,E)-1-(4-methyl-1-phenylpent-1-en-3-yl)pyrrolidine", "CC(C)[C@H](N1CCCC1)/C=C/c1ccccc1", "S", ["E"]),
    "1e": ("(R,E)-1-(1,3-diphenylallyl)pyrrolidine", "N1(CCCC1)[C@@H](c2ccccc2)/C=C/c3ccccc3", "R", ["E"]),
    "1f": ("(S,E)-1-(4-(2-chlorophenyl)but-3-en-2-yl)pyrrolidine", "C[C@H](N1CCCC1)/C=C/c1ccccc1Cl", "S", ["E"]),
    "1g": ("(S,E)-1-(1-(4-methoxyphenyl)pent-1-en-3-yl)pyrrolidine", "CC[C@H](N1CCCC1)/C=C/c1ccc(OC)cc1", "S", ["E"]),
    "1h": ("1-((S,3E,5E)-6-phenylhexa-3,5-dien-2-yl)pyrrolidine", "C[C@H](N1CCCC1)/C=C/C=C/c1ccccc1", "S", ["E", "E"]),
    "1i": ("(S,E)-1-(4-cyclohexylbut-3-en-2-yl)pyrrolidine", "C[C@H](N1CCCC1)/C=C/C2CCCCC2", "S", ["E"]),
    "1j": ("(S,E)-6-(benzyloxy)-N,N,2,5-tetramethylhex-4-en-3-amine", "CC(C)[C@H](N(C)C)/C=C(C)\\COCc1ccccc1", "S", ["E"]),
    "1k": ("(S,E)-1-(3-methyl-4-phenylbut-3-en-2-yl)pyrrolidine", "C[C@H](N1CCCC1)/C(=C/c2ccccc2)C", "S", ["E"]),
    "1l": ("(1S,5S)-N,N,2-trimethyl-5-(prop-1-en-2-yl)cyclohex-2-en-1-amine", "CN(C)[C@@H]1C(C)=CC[C@H](C(=C)C)C1", "1S,5S", []),
}
for label, (name, smi, center, geoms) in allyl.items():
    if label == "1a":
        continue
    page = 3 if label in {"1b", "1c", "1d"} else 4 if label in {"1e", "1f", "1g", "1h", "1i", "1j", "1k"} else 5
    add_molecule(f"m_{label}", label, name, (E[2], E[page], I[page]), smi, "specified",
                 "specified" if geoms else "not_applicable",
                 atom_centers=[{"center": "source-named stereocenter(s)", "configuration": center}],
                 double_bonds=[{"configuration": g} for g in geoms] if geoms else None,
                 assignment={"source": "compound heading/Table 3", "configuration": center})


rac_variants = {
    "1ab": ("racemic N-cyclohexyl? no: piperidinyl allylamine rac-1ab", "CC(N1CCCCC1)/C=C/c1ccccc1"),
    "1ac": ("racemic morpholinyl allylamine rac-1ac", "CC(N1CCOCC1)/C=C/c1ccccc1"),
    "1ad": ("racemic N,N-dimethyl allylamine rac-1ad", "CC(N(C)C)/C=C/c1ccccc1"),
    "1ae": ("racemic N-benzyl-N-methyl allylamine rac-1ae", "CC(N(C)Cc1ccccc1)/C=C/c1ccccc1"),
    "1af": ("racemic N-methyl-N-phenyl allylamine rac-1af", "CC(N(C)c1ccccc1)/C=C/c1ccccc1"),
    "1ag": ("racemic secondary N-methyl allylamine rac-1ag", "CC(NC)/C=C/c1ccccc1"),
    "1ah": ("racemic primary allylamine rac-1ah", "CC(N)/C=C/c1ccccc1"),
}
for label, (name, smi) in rac_variants.items():
    add_molecule(f"m_rac{label}", f"rac-{label}", name, (E[2], I[2]), smi,
                 "racemic", "specified", "partial",
                 atom_centers=[{"center": "allylic carbon", "configuration": "R/S racemate"}],
                 double_bonds=[{"configuration": "E"}],
                 mixture={"type": "racemate", "ratio": "not reported"},
                 assignment={"source": "rac label and Table 2 drawing"})


allenones = {
    "2b": ("1-(4-methoxyphenyl)buta-2,3-dien-1-one", "C=C=CC(=O)c1ccc(OC)cc1", "not_applicable"),
    "2c": ("1-(4-fluorophenyl)buta-2,3-dien-1-one", "C=C=CC(=O)c1ccc(F)cc1", "not_applicable"),
    "2d": ("1-(4-bromophenyl)buta-2,3-dien-1-one", "C=C=CC(=O)c1ccc(Br)cc1", "not_applicable"),
    "2e": ("1-(2-chlorophenyl)buta-2,3-dien-1-one", "C=C=CC(=O)c1ccccc1Cl", "not_applicable"),
    "2f": ("1-mesitylbuta-2,3-dien-1-one", "C=C=CC(=O)c1c(C)cc(C)cc1C", "not_applicable"),
    "2g": ("1-(furan-2-yl)buta-2,3-dien-1-one", "C=C=CC(=O)c1ccco1", "not_applicable"),
    "2h": ("1-(thiophen-2-yl)buta-2,3-dien-1-one", "C=C=CC(=O)c1cccs1", "not_applicable"),
    "2i": ("(E)-1-phenylocta? styryl-substituted allenone 2i", "C=C=CC(=O)/C=C/c1ccccc1", "specified"),
    "2j": ("1-benzylbuta-2,3-dien-1-one", "C=C=CC(=O)Cc1ccccc1", "not_applicable"),
    "2k": ("2-methyl-1-phenylbuta-2,3-dien-1-one", "C=C=C(C)C(=O)c1ccccc1", "not_applicable"),
    "2l": ("4-methyl-1-phenylpenta-2,3-dien-1-one", "CC(C)=C=CC(=O)c1ccccc1", "unresolved"),
}
for label, (name, smi, geom) in allenones.items():
    add_molecule(f"m_{label}", label, name, (E[3], I[3]), smi,
                 "achiral" if label != "2l" else "unresolved", geom,
                 "encoded" if geom == "specified" else "not_applicable" if geom == "not_applicable" else "unresolved",
                 double_bonds=[{"bond": "styryl C=C", "configuration": "E"}] if label == "2i" else None,
                 axial={"allene": "terminal or source-unspecified substituted allene; no axial assignment reported"},
                 assignment={"basis": "candidate graph from Table 4/prose; allene axial stereochemistry not reported"})

# RDKit exposes an unspecified potential cumulene element in 2i; the source fixes
# the styryl E bond but not an allene axis, so the representation is intentionally partial.
MOLS["m_2i"]["stereochemistry"]["representation_status"] = f(
    "styryl E encoded; allene axial assignment not reported", "partial", (E[3], I[3]))


product_data = {
    "4a": ("(S,E)-1,5-diphenyloct-6-ene-1,3-dione", "O=C(c1ccccc1)CC(=O)C[C@H](c1ccccc1)/C=C/C", "S", ["E"]),
    "4b": ("(S,E)-1,5-diphenylnon-6-ene-1,3-dione", "O=C(c1ccccc1)CC(=O)C[C@H](c1ccccc1)/C=C/CC", "S", ["E"]),
    "4c": ("(S,E)-1,5,9-triphenylnon-6-ene-1,3-dione", "O=C(c1ccccc1)CC(=O)C[C@H](c1ccccc1)/C=C/CCc1ccccc1", "S", ["E"]),
    "4d": ("(S,E)-8-methyl-1,5-diphenylnon-6-ene-1,3-dione", "O=C(c1ccccc1)CC(=O)C[C@H](c1ccccc1)/C=C/C(C)C", "S", ["E"]),
    "4e": ("(S,E)-1,5,7-triphenylhept-6-ene-1,3-dione", "O=C(c1ccccc1)CC(=O)C[C@H](c1ccccc1)/C=C/c1ccccc1", "S", ["E"]),
    "4f": ("(S,E)-5-(2-chlorophenyl)-1-phenyloct-6-ene-1,3-dione", "O=C(c1ccccc1)CC(=O)C[C@H](c1ccccc1Cl)/C=C/C", "S", ["E"]),
    "4g": ("(S,E)-5-(4-methoxyphenyl)-1-phenylnon-6-ene-1,3-dione", "O=C(c1ccccc1)CC(=O)C[C@H](c1ccc(OC)cc1)/C=C/CC", "S", ["E"]),
    "4h": ("(S,E)-1-phenyl-5-((E)-styryl)oct-6-ene-1,3-dione", "O=C(c1ccccc1)CC(=O)C[C@H](/C=C/c1ccccc1)/C=C/C", "S", ["E", "E"]),
    "4i": ("(R,E)-5-cyclohexyl-1-phenyloct-6-ene-1,3-dione", "O=C(c1ccccc1)CC(=O)C[C@H](C1CCCCC1)/C=C/C", "R", ["E"]),
    "4j": ("(S,E)-5-((benzyloxy)methyl)-5,8-dimethyl-1-phenylnon-6-ene-1,3-dione", "O=C(c1ccccc1)CC(=O)C[C@](C)(COCc1ccccc1)/C=C/C(C)C", "S", ["E"]),
    "4k": ("(S,E)-6-methyl-1,5-diphenyloct-6-ene-1,3-dione", "O=C(c1ccccc1)CC(=O)C[C@H](c1ccccc1)/C(=C/C)C", "S", ["E"]),
    "4l": ("4-((1S,5R)-2-methyl-5-(prop-1-en-2-yl)cyclohex-2-en-1-yl)-1-phenylbutane-1,3-dione", "O=C(c1ccccc1)CC(=O)C[C@H]1C(C)=CC[C@H](C(=C)C)C1", "1S,5R", []),
    "4m": ("(S,E)-1-(4-methoxyphenyl)-5-phenyloct-6-ene-1,3-dione", "O=C(c1ccc(OC)cc1)CC(=O)C[C@H](c1ccccc1)/C=C/C", "S", ["E"]),
    "4n": ("(S,E)-1-(4-fluorophenyl)-5-phenyloct-6-ene-1,3-dione", "O=C(c1ccc(F)cc1)CC(=O)C[C@H](c1ccccc1)/C=C/C", "S", ["E"]),
    "4o": ("(S,E)-1-(4-bromophenyl)-5-phenyloct-6-ene-1,3-dione", "O=C(c1ccc(Br)cc1)CC(=O)C[C@H](c1ccccc1)/C=C/C", "S", ["E"]),
    "4p": ("(S,E)-1-(2-chlorophenyl)-5-phenyloct-6-ene-1,3-dione", "O=C(c1ccccc1Cl)CC(=O)C[C@H](c1ccccc1)/C=C/C", "S", ["E"]),
    "4q": ("(S,E)-1-mesityl-5-phenyloct-6-ene-1,3-dione", "O=C(c1c(C)cc(C)cc1C)CC(=O)C[C@H](c1ccccc1)/C=C/C", "S", ["E"]),
    "4r": ("(S,E)-1-(furan-2-yl)-5-phenyloct-6-ene-1,3-dione", "O=C(c1ccco1)CC(=O)C[C@H](c1ccccc1)/C=C/C", "S", ["E"]),
    "4s": ("(S,E)-5-phenyl-1-(thiophen-2-yl)oct-6-ene-1,3-dione", "O=C(c1cccs1)CC(=O)C[C@H](c1ccccc1)/C=C/C", "S", ["E"]),
    "4t": ("(S,1E,8E)-1,7-diphenyldeca-1,8-diene-3,5-dione", "c1ccccc1/C=C/C(=O)CC(=O)C[C@H](c1ccccc1)/C=C/C", "S", ["E", "E"]),
    "4u": ("(S,E)-1,6-diphenylnon-7-ene-2,4-dione", "c1ccccc1CC(=O)CC(=O)C[C@H](c1ccccc1)/C=C/C", "S", ["E"]),
}
for label, (name, smi, center, geoms) in product_data.items():
    if label <= "4e":
        ev = (E[2], E[5], I[2], I[5])
    elif label <= "4o":
        ev = (E[2] if label <= "4l" else E[3], E[6], I[2] if label <= "4l" else I[3], I[6])
    else:
        ev = (E[3], E[7], I[3])
    assignment = {"configuration": center,
                  "basis": "4o assigned S by X-ray crystallography; other Table 3/4 products assigned by analogy"}
    add_molecule(f"m_{label}", label, name, ev, smi, "specified",
                 "specified" if geoms else "not_applicable",
                 atom_centers=[{"center": "source-named stereocenter(s)", "configuration": center}],
                 double_bonds=[{"configuration": g} for g in geoms] if geoms else None,
                 assignment=assignment,
                 note="SMILES represents the named diketone tautomer; the paper states products 4 exist predominantly in enol form by NMR.")


# Generic drawings on page 1.
add_template("m_generic_allylamine", "generic allylamine 1", "enantioenriched alpha-chiral allylamine template",
             (E[1], I[1]), "substituted allylic tertiary amine with R1-R6", ["R1", "R2", "R3", "R4", "R5", "R6"],
             "partially_specified", "unresolved", "Alpha chirality is shown conceptually; absolute configuration depends on the instance.")
add_template("m_generic_allenone", "generic allenone 2", "allenone template", (E[1], I[1]),
             "R-C(=O)-C=C=C", ["R"], "unresolved", "not_applicable")
add_template("m_generic_diketone", "generic product 4", "delta-chiral beta-diketone template", (E[1], I[1]),
             "R-substituted beta-diketone bearing transferred allyl substituents", ["R", "R3", "R4", "R5", "R6"],
             "partially_specified", "specified", "Generic product depicts chirality transfer and E alkene; absolute configuration varies.")
add_template("m_prev_allylamine", "Scheme 1a allylamine", "proline-derived alpha-chiral allylamine template", (E[1], I[1]),
             "proline-derived cyclic allylamine", ["R1", "R2"], "partially_specified", "unresolved")
add_template("m_prev_acyl", "R3CH2COX", "generic acyl electrophile R3CH2COX", (E[1], I[1]),
             "R3-CH2-C(O)-X", ["R3", "X"], "unresolved", "not_applicable")
add_template("m_prev_lactam", "Scheme 1a product", "ring-expanded lactam template", (E[1], I[1]),
             "substituted unsaturated ring-expanded lactam", ["R1", "R2", "R3"], "partially_specified", "unresolved")


# Page 1 generic/current/literature/mechanism records.
r = reaction("scheme1a_previous_work", (E[1], I[1]), "literature_example", ["other"],
             "Scheme 1a, previous work", "zwitterionic aza-Claisen ring expansion", "literature aza-Claisen/ring expansion", "final")
r["participants"] = [participant("prev_am", "m_prev_allylamine", "reactant", r["evidence_ids"]),
                     participant("prev_acyl", "m_prev_acyl", "reactant", r["evidence_ids"]),
                     participant("prev_prod", "m_prev_lactam", "product", r["evidence_ids"], product_role="literature product")]
set_procedure(r, r["evidence_ids"], "Scheme 1a: R3CH2COX, base, additive, followed by [3,3] rearrangement.")
r["review"]["unresolved"].append("Literature-example quantities, yield, identities of base/additive, and X are not reported in the figure.")
REACTIONS.append(r)

r = reaction("generic_aza_claisen_hydrolysis", (E[1], I[1]), "generic_scheme", ["other"],
             "graphical abstract and Scheme 1b", "asymmetric zwitterionic aza-Claisen rearrangement", "aza-Claisen rearrangement followed by hydrolysis", "final")
r["participants"] = [participant("gen_amine", "m_generic_allylamine", "reactant", r["evidence_ids"]),
                     participant("gen_allenone", "m_generic_allenone", "reactant", r["evidence_ids"]),
                     participant("gen_water", "m_water", "reagent", r["evidence_ids"], stage="workup"),
                     participant("gen_silica", "m_silica", "workup_agent", r["evidence_ids"], stage="workup"),
                     participant("gen_product", "m_generic_diketone", "product", r["evidence_ids"], product_role="generic main product")]
r["steps"] = [step("generic_rearrangement", "reaction", ["gen_amine", "gen_allenone"], r["evidence_ids"], "catalyst-free rearrangement", "r.t."),
              step("generic_hydrolysis", "workup", ["gen_water", "gen_silica", "gen_product"], r["evidence_ids"], "silica gel/water hydrolysis", "r.t.")]
set_procedure(r, r["evidence_ids"], "1) catalyst-free; 2) silica gel, water (Scheme 1b); graphical abstract: MeOH, rt then H2O/silica gel/CH2Cl2, rt.")
r["review"]["unresolved"].append("The generic drawings contain R groups and do not define one concrete experiment.")
REACTIONS.append(r)

r = reaction("scheme1b_proposed_mechanism", (E[1], I[1]), "proposed_mechanism", ["other"],
             "Scheme 1b lower pathway", "proposed 1,3-chirality-transfer mechanism", "nucleophilic addition, [3,3] rearrangement, hydrolysis", "final")
r["participants"] = [participant("mech_am", "m_generic_allylamine", "reactant", r["evidence_ids"]),
                     participant("mech_allenone", "m_generic_allenone", "reactant", r["evidence_ids"]),
                     participant("mech_water", "m_water", "reagent", r["evidence_ids"]),
                     participant("mech_product", "m_generic_diketone", "product", r["evidence_ids"], product_role="mechanistic endpoint")]
r["steps"] = [step("addition", "reaction", ["mech_am", "mech_allenone"], r["evidence_ids"], "proposed nucleophilic addition"),
              step("sigmatropic", "reaction", ["mech_am", "mech_allenone"], r["evidence_ids"], "proposed chair-like [3,3] rearrangement"),
              step("hydrolysis", "reaction", ["mech_water", "mech_product"], r["evidence_ids"], "proposed hydrolysis with amine release")]
set_procedure(r, r["evidence_ids"], "Arrow-pushing pathway in Scheme 1b; this is a proposal, not a reported experimental procedure.")
r["review"]["unresolved"].append("No mechanistic intermediate was isolated or quantified on page 1.")
REACTIONS.append(r)


# Table 1, fourteen independent solvent-screening experiments.
solvents = [
    (1, "m_dcm", 85), (2, "m_chcl3", 82), (3, "m_dce", 50), (4, "m_toluene", 68),
    (5, "m_hexane", 45), (6, "m_etoac", 61), (7, "m_thf", 70), (8, "m_dioxane", 60),
    (9, "m_acetone", 94), (10, "m_mecn", 65), (11, "m_dmf", 92), (12, "m_dmso", 89),
    (13, "m_meoh", 98), (14, "m_etoh", 83),
]
table1_ids = []
for entry, solvent_mid, yld in solvents:
    ev = (E[2], I[2], E[5], I[5]) + ((E[1], I[1]) if entry == 1 else ())
    rid = f"table1_solvent_entry_{entry:02d}"
    r = reaction(rid, ev, categories=["optimization"], location=f"Table 1 entry {entry}",
                 name="aza-Claisen rearrangement solvent screening", transformation="allylamine + allenone to beta-amino enone", class_mode="intermediate")
    r["participants"] = [participant("rac1a", "m_rac1a", "reactant", ev, amount=q(0.20, "mmol")),
                         participant("a2", "m_2a", "reactant", ev, amount=q(0.24, "mmol")),
                         participant("solvent", solvent_mid, "solvent", ev, volume=q(1.0, "mL")),
                         participant("rac3a", "m_rac3a", "product", ev, product_role="crude beta-amino enone")]
    r["reaction_scale"] = f("rac-1a (0.20 mmol)", {"reference_molecule_id": "m_rac1a", "amount": q(0.20, "mmol")}, ev)
    r["solvent_conditions"]["composition"] = f(MOLS[solvent_mid]["names"]["original"]["normalized_value"], [solvent_mid], ev)
    r["solvent_conditions"]["total_volume"] = f("1.0 mL", q(1.0, "mL"), ev)
    r["steps"] = [step("screen_reaction", "reaction", ["rac1a", "a2", "solvent", "rac3a"], ev,
                       "stirred", "r.t.", q(12, "h"))]
    r["outcomes"] = [outcome("nmr_yield", "yield", ["m_rac3a"], ev, q(yld, "%"), "NMR yield",
                             "1H NMR spectroscopic analysis of crude product", {"table": 1, "entry": entry})]
    set_procedure(r, ev, "rac-1a and 2a were stirred in the listed solvent at room temperature for 12 h, then concentrated under reduced pressure.",
                  "General Procedure for the Screening of Solvents (Table 1)",
                  addition_order="2a added to a stirred solution of rac-1a in solvent")
    REACTIONS.append(r); table1_ids.append(r["reaction_id"])


# Shared final-reaction procedure for Tables 2--4.
def add_final_participants(r, amine_mid, allene_mid, product_mid, ev, product=True):
    parts = [participant("amine", amine_mid, "reactant", ev, amount=q(0.20, "mmol")),
             participant("allenone", allene_mid, "reactant", ev, amount=q(0.24, "mmol")),
             participant("meoh", "m_meoh", "solvent", ev, volume=q(1.0, "mL")),
             participant("dcm", "m_dcm", "solvent", ev, stage="workup", volume=q(5.0, "mL")),
             participant("water", "m_water", "workup_agent", ev, stage="workup"),
             participant("silica", "m_silica", "workup_agent", ev, stage="workup", mass=q(500, "mg")),
             participant("petroleum", "m_petroleum_ether", "purification_agent", ev, stage="purification"),
             participant("etoac", "m_etoac", "purification_agent", ev, stage="purification")]
    if product:
        parts.append(participant("product", product_mid, "product", ev, product_role="main product"))
    r["participants"] = parts
    r["reaction_scale"] = f("allylamine (0.20 mmol)", {"reference_molecule_id": amine_mid, "amount": q(0.20, "mmol")}, ev)
    r["solvent_conditions"]["composition"] = f("methanol; then dichloromethane containing 1.0% (v/v) water", ["m_meoh", "m_dcm", "m_water"], ev)
    r["solvent_conditions"]["mixture_ratio"] = f("CH2Cl2 containing 1.0% (v/v) water", {"water": q(1.0, "%"), "basis": "v/v"}, ev)
    r["solvent_conditions"]["total_volume"] = f("MeOH 1.0 mL; CH2Cl2/water 5.0 mL", {"reaction": q(1.0, "mL"), "workup_stage": q(5.0, "mL")}, ev)
    r["steps"] = [
        step("rearrangement", "reaction", ["amine", "allenone", "meoh"], ev, "stirred", "r.t.", q(12, "h")),
        step("remove_methanol", "workup", ["amine", "allenone"], ev, "concentrated under reduced pressure to remove methanol", "r.t."),
        step("silica_hydrolysis", "workup", ["dcm", "water", "silica"] + (["product"] if product else []), ev, "residue diluted with wet dichloromethane; silica gel added and mixture stirred", "r.t.", q(12, "h")),
        step("filter_wash", "workup", ["dcm", "silica"] + (["product"] if product else []), ev, "filtered; silica gel washed with dichloromethane (3 × 10 mL); combined filtrates concentrated"),
        step("chromatography", "purification", ["petroleum", "etoac"] + (["product"] if product else []), ev, "silica gel chromatography, petroleum ether/ethyl acetate = 10:1"),
    ]
    set_procedure(r, ev, "General procedure: allylamine (0.20 mmol) and allenone (0.24 mmol) in MeOH (1.0 mL), rt 12 h; remove MeOH; wet CH2Cl2/silica, rt 12 h; filter, wash, concentrate, chromatography.",
                  "General Procedure for the Aza-Claisen Rearrangement between Allylamines and Allenones (Tables 2-4)",
                  addition_order="allenone added to a stirred solution of allylamine in methanol")


# Table 2 N-substituent survey.
table2 = [(1, "m_rac1a", 95), (2, "m_rac1ab", 85), (3, "m_rac1ac", 82),
          (4, "m_rac1ad", 90), (5, "m_rac1ae", 93), (6, "m_rac1af", 34),
          (7, "m_rac1ag", 0), (8, "m_rac1ah", 0)]
table2_ids = []
for entry, amine_mid, yld in table2:
    ev = (E[2], I[2], E[5], I[5])
    r = reaction(f"table2_n_substituent_entry_{entry:02d}", ev, categories=["optimization"],
                 location=f"Table 2 entry {entry}", name="aza-Claisen rearrangement/hydrolysis N-substituent survey",
                 transformation="allylamine + allenone to beta-diketone" if yld else "no observed aza-Claisen rearrangement",
                 class_mode="final" if yld else "none")
    add_final_participants(r, amine_mid, "m_2a", "m_rac4a", ev, product=True)
    r["outcomes"] = [outcome("isolated_yield", "yield", ["m_rac4a"], ev, q(yld, "%"), "isolated yield", "isolated", {"table": 2, "entry": entry})]
    if yld == 0:
        r["outcomes"].append(outcome("no_reaction", "no_reaction", [], ev, True, basis=None,
                                     scope={"source_statement": "aza-Claisen rearrangement was not observed"}, product_role=None))
        r["review"]["unresolved"].append("The source reports 0% isolated yield/no rearrangement but no conversion or recovery measurement.")
    REACTIONS.append(r); table2_ids.append(r["reaction_id"])


# Table 3 enantioenriched allylamine scope.
table3 = [
    (1, "1a", "4a", 95, 97, "ee", 97), (2, "1b", "4b", 83, 97, "ee", 97),
    (3, "1c", "4c", 90, 96, "ee", 96), (4, "1d", "4d", 90, 98, "ee", 98),
    (5, "1e", "4e", 78, 96, "ee", 96), (6, "1f", "4f", 90, 99, "ee", 99),
    (7, "1g", "4g", 90, 80, "ee", 80), (8, "1h", "4h", 89, 99, "ee", 99),
    (9, "1i", "4i", 88, 97, "ee", 97), (10, "1j", "4j", 87, ">99", "ee", 95),
    (11, "1k", "4k", 85, 99, "ee", 99), (12, "1l", "4l", 91, "95:5", "dr", "93:7"),
]
table3_ids = []
for entry, am, prod, yld, input_sel, metric, output_sel in table3:
    ev = (E[2], I[2], E[3], I[3], E[5], I[5])
    if entry >= 6:
        ev += (E[6], I[6])
    r = reaction(f"table3_allylamine_scope_entry_{entry:02d}", ev, categories=["substrate_scope"],
                 location=f"Table 3 entry {entry}", name="asymmetric aza-Claisen rearrangement followed by hydrolysis",
                 transformation="enantioenriched allylamine + allenone to delta-chiral beta-diketone", class_mode="final")
    add_final_participants(r, f"m_{am}", "m_2a", f"m_{prod}", ev)
    r["general_scheme_ids"] = [f"{DOC}:generic_aza_claisen_hydrolysis"]
    r["outcomes"] = [outcome("isolated_yield", "yield", [f"m_{prod}"], ev, q(yld, "%"), "isolated yield", "isolated", {"table": 3, "entry": entry}),
                     outcome("substrate_selectivity", "selectivity", [f"m_{am}"], ev,
                             q(input_sel, "%") if isinstance(input_sel, (int, float)) else input_sel,
                             basis="HPLC on a chiral stationary phase" if metric == "ee" else "1H NMR spectroscopic analysis",
                             scope={"metric": metric, "role": "input allylamine"}, product_role="reactant sample"),
                     outcome("product_selectivity", "selectivity", [f"m_{prod}"], ev,
                             q(output_sel, "%") if isinstance(output_sel, (int, float)) else output_sel,
                             basis="HPLC on a chiral stationary phase" if metric == "ee" else "1H NMR spectroscopic analysis of crude product",
                             scope={"metric": metric, "role": "product", "isomer_order": "major:minor where a ratio is reported"})]
    if entry <= 11:
        r["outcomes"].append(outcome("alkene_geometry", "selectivity", [f"m_{prod}"], (E[3], I[3]),
                                             "exclusive E", basis="source statement; method not separately reported",
                                             scope={"metric": "E/Z", "major_isomer": "E"}))
    if prod == "4o":
        r["notes"].append(f("Absolute configuration of 4o was assigned S by X-ray crystallography (CCDC 1489729).", "4o: S by X-ray; CCDC 1489729", (E[3], I[3])))
    else:
        r["notes"].append(f("Absolute configuration assigned by analogy to 4o.", "assignment by analogy", (E[3], I[3])))
    REACTIONS.append(r); table3_ids.append(r["reaction_id"])


# Table 4 allenone scope.
table4 = [(1, "2b", "4m", 88), (2, "2c", "4n", 79), (3, "2d", "4o", 89),
          (4, "2e", "4p", 80), (5, "2f", "4q", 61), (6, "2g", "4r", 83),
          (7, "2h", "4s", 83), (8, "2i", "4t", 75), (9, "2j", "4u", 82)]
table4_ids = []
for entry, allene, prod, yld in table4:
    ev = (E[3], I[3], E[5], I[5])
    if entry <= 3:
        ev += (E[6], I[6])
    else:
        ev += (E[7],)
    r = reaction(f"table4_allenone_scope_entry_{entry:02d}", ev, categories=["substrate_scope"],
                 location=f"Table 4 entry {entry}", name="asymmetric aza-Claisen rearrangement followed by hydrolysis",
                 transformation="enantioenriched allylamine + allenone to delta-chiral beta-diketone", class_mode="final")
    add_final_participants(r, "m_1a", f"m_{allene}", f"m_{prod}", ev)
    r["general_scheme_ids"] = [f"{DOC}:generic_aza_claisen_hydrolysis"]
    r["outcomes"] = [outcome("isolated_yield", "yield", [f"m_{prod}"], ev, q(yld, "%"), "isolated yield", "isolated", {"table": 4, "entry": entry}),
                     outcome("product_ee", "selectivity", [f"m_{prod}"], ev, q(97, "%"),
                             basis="HPLC on a chiral stationary phase", scope={"metric": "ee"}),
                     outcome("alkene_geometry", "selectivity", [f"m_{prod}"], ev, "exclusive E",
                             basis="source statement; method not separately reported", scope={"metric": "E/Z", "major_isomer": "E"})]
    r["notes"].append(f("4o absolute configuration assigned S by X-ray; other Table 4 products assigned by analogy.",
                         "4o X-ray; others by analogy", (E[3], I[3])))
    REACTIONS.append(r); table4_ids.append(r["reaction_id"])


# Two explicitly reported failed bulky-allenone reactions, separate from Table 4.
negative_ids = []
for label in ("2k", "2l"):
    ev = (E[3], I[3], E[5], I[5])
    r = reaction(f"bulky_allenone_{label}_no_reaction", ev, categories=["substrate_scope", "control"],
                 location=f"page 3 prose, {label}", name="aza-Claisen rearrangement under standard conditions",
                 transformation="no reaction with bulky allenone", class_mode="none")
    add_final_participants(r, "m_1a", f"m_{label}", "m_4a", ev, product=False)
    r["outcomes"] = [outcome("no_reaction", "no_reaction", [], ev, True,
                             scope={"allenone": label, "source": "reaction did not occur"}, product_role=None)]
    r["review"]["unresolved"].append("No numerical conversion, recovery, detection limit, or product identity is reported for this failed reaction.")
    REACTIONS.append(r); negative_ids.append(r["reaction_id"])


# Preparation of allylamines 1a-i and 1k: one generic scheme plus ten product-specific experiments.
add_template("m_S1_generic", "S1", "chiral primary allylamine S1 template", (E[3], I[3]),
             "R1/R2/R3-substituted chiral primary allylamine", ["R1", "R2", "R3"],
             "partially_specified", "partially_specified", "Prepared according to cited literature; individual precursor quantities/identities are not separately listed.")
add_common("m_dibromobutane", "Br(CH2)4Br", "1,4-dibromobutane", "BrCCCCBr", (E[3], I[3]))
add_common("m_nahco3", "NaHCO3", "sodium bicarbonate", "[Na+].OC([O-])=O", (E[3], I[3]))
prep_generic = reaction("generic_preparation_allylamines_1a_i_1k", (E[3], I[3]), "generic_scheme", ["other"],
                        "page 3 preparation scheme", "cyclization of chiral primary allylamines with 1,4-dibromobutane",
                        "N-alkylation/ring closure to pyrrolidines", "cyclization")
prep_generic["participants"] = [participant("S1", "m_S1_generic", "reactant", prep_generic["evidence_ids"]),
                                participant("dibromo", "m_dibromobutane", "reactant", prep_generic["evidence_ids"]),
                                participant("base", "m_nahco3", "reagent", prep_generic["evidence_ids"]),
                                participant("solvent", "m_toluene", "solvent", prep_generic["evidence_ids"])]
prep_generic["steps"] = [step("generic_cyclization", "reaction", ["S1", "dibromo", "base", "solvent"], prep_generic["evidence_ids"], "heated", q(115, "°C"), q(23, "h"), "nitrogen")]
set_procedure(prep_generic, prep_generic["evidence_ids"], "S1, 1,4-dibromobutane and NaHCO3 in toluene; nitrogen, 115 °C, 23 h; cool and chromatograph.")
REACTIONS.append(prep_generic)

prep_products = [("1a", 154, 76, 3), ("1b", 148, 69, 3), ("1c", 120, 41, 3), ("1d", 118, 51, 3),
                 ("1e", 105, 40, 4), ("1f", 102, 43, 4), ("1g", 99.0, 40, 4), ("1h", 120, 53, 4),
                 ("1i", 60.0, 29, 4), ("1k", 92.0, 43, 4)]
prep_ids = []
for label, mass_mg, yld, page in prep_products:
    ev = (E[3], I[3], E[page], I[page])
    precursor_mid = f"m_S1_for_{label}"
    add_molecule(precursor_mid, f"S1 precursor to {label}", f"corresponding chiral primary allylamine S1 precursor to {label}", ev,
                 chirality="specified", geometry="specified", representation="not_encoded",
                 atom_centers=[{"configuration": allyl[label][2]}],
                 double_bonds=[{"configuration": g} for g in allyl[label][3]],
                 assignment={"basis": "corresponding precursor relationship in the general preparation; full identity not individually printed"},
                 note="Specific precursor structure is implied by the corresponding product but not independently named; no SMILES asserted.")
    r = reaction(f"preparation_{label}", ev, categories=["other"], location=f"Preparation of Allylamines 1a-i and 1k; {label} paragraph",
                 name="preparation of tertiary allylamine", transformation="double N-alkylation/ring closure", class_mode="cyclization")
    r["general_scheme_ids"] = [prep_generic["reaction_id"]]
    r["participants"] = [participant("S1", precursor_mid, "reactant", ev, amount=q(1.0, "mmol")),
                         participant("dibromo", "m_dibromobutane", "reactant", ev, mass=q(238, "mg"), amount=q(1.1, "mmol")),
                         participant("base", "m_nahco3", "reagent", ev, mass=q(185, "mg"), amount=q(2.2, "mmol")),
                         participant("toluene", "m_toluene", "solvent", ev, volume=q(1.0, "mL")),
                         participant("dcm", "m_dcm", "purification_agent", ev, stage="purification"),
                         participant("meoh", "m_meoh", "purification_agent", ev, stage="purification"),
                         participant("product", f"m_{label}", "product", ev, product_role="isolated tertiary allylamine")]
    r["reaction_scale"] = f("S1 (1.0 mmol)", {"reference_molecule_id": precursor_mid, "amount": q(1.0, "mmol")}, ev)
    r["solvent_conditions"]["composition"] = f("toluene", ["m_toluene"], ev)
    r["solvent_conditions"]["total_volume"] = f("1.0 mL", q(1.0, "mL"), ev)
    r["steps"] = [step("cyclization", "reaction", ["S1", "dibromo", "base", "toluene"], ev, "heated", q(115, "°C"), q(23, "h"), "nitrogen"),
                  step("cool", "workup", ["product"], ev, "cooled to room temperature", "r.t."),
                  step("purify", "purification", ["dcm", "meoh", "product"], ev, "silica gel chromatography, dichloromethane/methanol 15:1")]
    r["outcomes"] = [outcome("product_mass", "product_amount", [f"m_{label}"], ev, q(mass_mg, "mg"), basis="isolated material", scope={"product": label}),
                     outcome("isolated_yield", "yield", [f"m_{label}"], ev, q(yld, "%"), "isolated yield", "isolated", {"product": label})]
    set_procedure(r, ev, "S1 (1.0 mmol), 1,4-dibromobutane (238 mg, 1.1 mmol), NaHCO3 (185 mg, 2.2 mmol), toluene (1.0 mL), N2, 115 °C, 23 h; cool; chromatography DCM/MeOH 15:1.",
                  "Preparation of Allylamines 1a-i and 1k")
    REACTIONS.append(r); prep_ids.append(r["reaction_id"])


# Page 4: preparation of 1j, retained as two experimentally isolated spans.
for mid, label, name, smi in [
    ("m_S2", "S2", "amide S2 (Weinreb amide shown in Scheme)", None),
    ("m_S3", "S3", "phosphorus ylide S3", None),
    ("m_S4", "S4", "methyl (S,E)-4-((tert-butoxycarbonyl)amino)-2,5-dimethylhex-2-enoate", "COC(=O)/C(C)=C/[C@H](NC(=O)OC(C)(C)C)C(C)C"),
    ("m_S5", "S5", "methyl (S,E)-4-(dimethylamino)-2,5-dimethylhex-2-enoate", "COC(=O)/C(C)=C/[C@H](N(C)C)C(C)C"),
]:
    if smi:
        add_molecule(mid, label, name, (E[4], I[4]), smi, "specified", "specified",
                     atom_centers=[{"configuration": "S"}], double_bonds=[{"configuration": "E"}],
                     assignment={"source": "scheme/compound heading"})
    else:
        add_molecule(mid, label, name, (E[4], I[4]), chirality="partially_specified", geometry="unresolved",
                     representation="not_encoded", assignment={"source": "structure drawing on page 4"},
                     note="Structure is present as a drawing, but a concrete candidate SMILES was not asserted in this extraction.")
add_molecule("m_crude_aldehyde", "crude aldehyde", "crude aldehyde from S2", (E[4], I[4]),
             chirality="specified", geometry="not_applicable", representation="not_encoded",
             atom_centers=[{"configuration": "S; retained precursor center"}],
             assignment={"source": "used directly without purification"}, note="Crude intermediate; amount and yield unreported.")
add_molecule("m_crude_alcohol", "crude alcohol", "crude alcohol from S5", (E[4], I[4]),
             chirality="specified", geometry="specified", representation="not_encoded",
             atom_centers=[{"configuration": "S"}], double_bonds=[{"configuration": "E"}],
             assignment={"source": "used without purification"}, note="Crude intermediate; amount and yield unreported.")
reagents_1j = [
    ("m_lialh4", "LiAlH4", "lithium aluminium hydride", "[Li+].[AlH4-]"),
    ("m_khso4", "KHSO4", "potassium bisulfate", "[K+].OS(=O)(=O)[O-]"),
    ("m_tfa", "CF3CO2H", "trifluoroacetic acid", "O=C(O)C(F)(F)F"),
    ("m_naoh", "NaOH", "sodium hydroxide", "[Na+].[OH-]"),
    ("m_formaldehyde", "HCHO", "formaldehyde", "C=O"),
    ("m_nabh3cn", "NaBH3CN", "sodium cyanoborohydride", "[Na+].[BH3-]C#N"),
    ("m_acetic_acid", "AcOH", "acetic acid", "CC(=O)O"),
    ("m_benzyl_bromide", "BnBr", "benzyl bromide", "BrCc1ccccc1"),
]
for args in reagents_1j:
    add_common(*args, ev=(E[4], I[4]))

r = reaction("preparation_1j_s2_to_s4", (E[4], I[4]), categories=["other"],
             location="Preparation of Allylamine 1j, S2 to S4", name="amide reduction/Wittig olefination",
             transformation="S2 through crude aldehyde to S4", class_mode="wittig")
r["participants"] = [participant("S2", "m_S2", "reactant", r["evidence_ids"], mass=q(1.30, "g"), amount=q(5.0, "mmol")),
                     participant("lialh4", "m_lialh4", "reagent", r["evidence_ids"], volume=q(0.63, "mL"), amount=q(2.5, "mmol")),
                     participant("thf", "m_thf", "solvent", r["evidence_ids"], volume=q(25, "mL")),
                     participant("aldehyde", "m_crude_aldehyde", "product", r["evidence_ids"], product_role="crude intermediate"),
                     participant("S3", "m_S3", "reactant", r["evidence_ids"], mass=q(1.92, "g"), amount=q(5.5, "mmol")),
                     participant("dcm", "m_dcm", "solvent", r["evidence_ids"], volume=q(20, "mL")),
                     participant("S4", "m_S4", "product", r["evidence_ids"], product_role="isolated intermediate"),
                     participant("khso4", "m_khso4", "workup_agent", r["evidence_ids"], stage="workup", volume=q(15, "mL")),
                     participant("ether", "m_ether", "workup_agent", r["evidence_ids"], stage="workup", volume=q(40, "mL")),
                     participant("nahco3", "m_nahco3", "workup_agent", r["evidence_ids"], stage="workup"),
                     participant("brine", "m_brine", "workup_agent", r["evidence_ids"], stage="workup"),
                     participant("na2so4", "m_na2so4", "workup_agent", r["evidence_ids"], stage="workup"),
                     participant("petroleum", "m_petroleum_ether", "purification_agent", r["evidence_ids"], stage="purification"),
                     participant("etoac", "m_etoac", "purification_agent", r["evidence_ids"], stage="purification")]
r["reaction_scale"] = f("S2 (5.0 mmol)", {"reference_molecule_id": "m_S2", "amount": q(5.0, "mmol")}, r["evidence_ids"])
r["solvent_conditions"]["composition"] = f("dry THF then dry dichloromethane", ["m_thf", "m_dcm"], r["evidence_ids"])
r["steps"] = [step("lialh4_addition", "reaction", ["S2", "lialh4", "thf"], r["evidence_ids"], "LiAlH4 solution added dropwise", q(0, "°C")),
              step("reduction_hold", "reaction", ["S2", "lialh4", "thf"], r["evidence_ids"], "stirred", q(0, "°C"), q(30, "min"), "nitrogen"),
              step("cool_and_quench", "workup", ["khso4", "ether", "aldehyde"], r["evidence_ids"], "cooled to -15 °C; saturated KHSO4 added carefully; diluted with ether and stirred vigorously", q(-15, "°C"), q(20, "min")),
              step("aldehyde_workup", "workup", ["brine", "na2so4", "aldehyde"], r["evidence_ids"], "organic layer washed, dried, filtered, concentrated; crude aldehyde used directly"),
              step("wittig", "reaction", ["aldehyde", "S3", "dcm"], r["evidence_ids"], "S3 added in one portion; stirred", "r.t.", q(12, "h"), "nitrogen"),
              step("wittig_workup", "workup", ["khso4", "nahco3", "brine", "na2so4", "S4"], r["evidence_ids"], "washed successively with 10% KHSO4, saturated NaHCO3 and brine; dried and concentrated"),
              step("s4_chromatography", "purification", ["petroleum", "etoac", "S4"], r["evidence_ids"], "column chromatography, petroleum ether/ethyl acetate 5:1")]
r["outcomes"] = [outcome("s4_mass", "product_amount", ["m_S4"], r["evidence_ids"], q(1.20, "g"), basis="isolated", scope={"product": "S4"}),
                 outcome("s4_yield", "yield", ["m_S4"], r["evidence_ids"], q(88, "%"), "isolated yield", "isolated", {"span": "S2 to S4"})]
set_procedure(r, r["evidence_ids"], "S2 was reduced with LiAlH4 to a crude aldehyde, followed without purification by Wittig reaction with S3 and chromatography to S4.")
REACTIONS.append(r)

r = reaction("preparation_1j_s4_to_1j", (E[4], I[4]), categories=["other"],
             location="Preparation of Allylamine 1j, four steps from S4", name="Boc deprotection/reductive methylation/ester reduction/benzylation",
             transformation="S4 through S5 and crude alcohol to 1j", class_mode="methylation_ether")
r["participants"] = [participant("S4", "m_S4", "reactant", r["evidence_ids"], mass=q(1.09, "g"), amount=q(4.0, "mmol")),
                     participant("tfa", "m_tfa", "reagent", r["evidence_ids"], volume=q(2.0, "mL")),
                     participant("dcm", "m_dcm", "solvent", r["evidence_ids"], volume=q(8.0, "mL")),
                     participant("naoh", "m_naoh", "workup_agent", r["evidence_ids"], stage="workup", volume=q(4.0, "mL")),
                     participant("mecn", "m_mecn", "solvent", r["evidence_ids"], volume=q(4.0, "mL")),
                     participant("hcho", "m_formaldehyde", "reagent", r["evidence_ids"], volume=q(1.0, "mL")),
                     participant("nabh3cn", "m_nabh3cn", "reagent", r["evidence_ids"], mass=q(628, "mg"), amount=q(10, "mmol")),
                     participant("acoh", "m_acetic_acid", "additive", r["evidence_ids"]),
                     participant("S5", "m_S5", "product", r["evidence_ids"], product_role="crude intermediate"),
                     participant("lialh4", "m_lialh4", "reagent", r["evidence_ids"], mass=q(159, "mg"), amount=q(4.2, "mmol")),
                     participant("thf", "m_thf", "solvent", r["evidence_ids"], volume=q(20, "mL")),
                     participant("alcohol", "m_crude_alcohol", "product", r["evidence_ids"], product_role="crude intermediate"),
                     participant("bnbr", "m_benzyl_bromide", "reactant", r["evidence_ids"], mass=q(684, "mg"), amount=q(4.0, "mmol")),
                     participant("water", "m_water", "solvent", r["evidence_ids"], volume=q(5.0, "mL")),
                     participant("khso4", "m_khso4", "workup_agent", r["evidence_ids"], stage="workup", volume=q(5.0, "mL")),
                     participant("ether", "m_ether", "workup_agent", r["evidence_ids"], stage="workup"),
                     participant("etoac", "m_etoac", "workup_agent", r["evidence_ids"], stage="workup"),
                     participant("brine", "m_brine", "workup_agent", r["evidence_ids"], stage="workup"),
                     participant("na2so4", "m_na2so4", "workup_agent", r["evidence_ids"], stage="workup"),
                     participant("petroleum", "m_petroleum_ether", "purification_agent", r["evidence_ids"], stage="purification"),
                     participant("product", "m_1j", "product", r["evidence_ids"], product_role="isolated allylamine")]
r["reaction_scale"] = f("S4 (4.0 mmol)", {"reference_molecule_id": "m_S4", "amount": q(4.0, "mmol")}, r["evidence_ids"])
r["steps"] = [step("tfa_addition", "reaction", ["S4", "tfa", "dcm"], r["evidence_ids"], "TFA added slowly", q(0, "°C")),
              step("deprotection", "reaction", ["S4", "tfa", "dcm"], r["evidence_ids"], "stirred", "r.t.", q(30, "min")),
              step("deprotection_workup", "workup", ["naoh"], r["evidence_ids"], "solvent removed; 2.0 M NaOH added; extracted with DCM 3 × 5 mL; dried and concentrated"),
              step("reductive_methylation", "reaction", ["mecn", "hcho", "nabh3cn", "acoh", "S5"], r["evidence_ids"], "formaldehyde then NaBH3CN added; acetic acid used to keep pH near neutral; stirred", "r.t.", q(2, "h")),
              step("s5_workup", "workup", ["naoh", "etoac", "brine", "na2so4", "S5"], r["evidence_ids"], "concentrated; 1 M NaOH added; extracted, washed, dried, concentrated to crude S5"),
              step("ester_reduction", "reaction", ["S5", "lialh4", "thf", "alcohol"], r["evidence_ids"], "LiAlH4 added; stirred", "r.t.", q(4, "h")),
              step("alcohol_workup", "workup", ["khso4", "ether", "na2so4", "alcohol"], r["evidence_ids"], "saturated KHSO4 added under vigorous stirring; phases separated; aqueous phase extracted; dried and evaporated"),
              step("benzylation", "reaction", ["alcohol", "bnbr", "naoh", "water", "product"], r["evidence_ids"], "crude alcohol and benzyl bromide stirred vigorously in 50% aqueous NaOH", "r.t.", q(6, "h")),
              step("final_workup", "workup", ["dcm", "na2so4", "product"], r["evidence_ids"], "aqueous phase extracted with DCM 3 × 10 mL; dried and concentrated"),
              step("final_chromatography", "purification", ["petroleum", "etoac", "product"], r["evidence_ids"], "silica gel chromatography, petroleum ether/ethyl acetate 2:1")]
r["outcomes"] = [outcome("product_mass", "product_amount", ["m_1j"], r["evidence_ids"], q(123, "mg"), basis="isolated", scope={"product": "1j"}),
                 outcome("overall_yield", "yield", ["m_1j"], r["evidence_ids"], q(12, "%"), "four-step isolated yield from S4", "isolated", {"start": "S4", "end": "1j", "steps": 4})]
set_procedure(r, r["evidence_ids"], "Four-step sequence from S4: TFA deprotection, formaldehyde/NaBH3CN methylation to crude S5, LiAlH4 ester reduction to crude alcohol, benzyl bromide/aqueous NaOH benzylation, chromatography.")
r["review"]["unresolved"].append("Amounts of crude amine S5 and crude alcohol are not reported; no quantitative intermediate yield was inferred.")
REACTIONS.append(r)


# Page 5: preparation of 1l as two real arrow experiments.
add_common("m_carvone", "(S)-(+)-carvone", "(S)-(+)-carvone", "CC1=CC[C@@H](CC1=O)C(=C)C", (E[5], I[5]))
add_common("m_ti_oipr4", "Ti(OiPr)4", "titanium(IV) isopropoxide", "CC(C)O[Ti](OC(C)C)(OC(C)C)OC(C)C", (E[5], I[5]))
add_common("m_methylamine", "MeNH2", "methylamine", "CN", (E[5], I[5]))
add_common("m_nabh4", "NaBH4", "sodium borohydride", "[Na+].[BH4-]", (E[5], I[5]))
add_molecule("m_S6", "S6", "carvone-derived secondary amine S6", (E[5], I[5]),
             chirality="specified", geometry="not_applicable", representation="not_encoded",
             atom_centers=[{"configuration": "source drawing retains carvone-derived stereochemistry"}],
             assignment={"source": "S6 structure drawing"}, note="Crude intermediate used directly; no isolated yield or concrete SMILES asserted.")
r = reaction("preparation_1l_carvone_to_s6", (E[5], I[5]), categories=["other"],
             location="Preparation of Allylamine 1l, first arrow", name="reductive amination of (S)-(+)-carvone",
             transformation="carvone to crude secondary amine S6", class_mode="amination")
r["participants"] = [participant("carvone", "m_carvone", "reactant", r["evidence_ids"], mass=q(300, "mg"), amount=q(2.0, "mmol")),
                     participant("ti", "m_ti_oipr4", "reagent", r["evidence_ids"], volume=q(1.2, "mL"), mass=q(1.14, "g"), amount=q(4.0, "mmol")),
                     participant("methylamine", "m_methylamine", "reactant", r["evidence_ids"], volume=q(5.0, "mL"), amount=q(10, "mmol")),
                     participant("meoh", "m_meoh", "solvent", r["evidence_ids"], volume=q(10, "mL")),
                     participant("nabh4", "m_nabh4", "reagent", r["evidence_ids"], mass=q(113, "mg"), amount=q(3.0, "mmol")),
                     participant("water", "m_water", "workup_agent", r["evidence_ids"], stage="workup", volume=q(3, "mL")),
                     participant("celite", "m_celite", "workup_agent", r["evidence_ids"], stage="workup"),
                     participant("etoac", "m_etoac", "workup_agent", r["evidence_ids"], stage="workup"),
                     participant("na2so4", "m_na2so4", "workup_agent", r["evidence_ids"], stage="workup"),
                     participant("S6", "m_S6", "product", r["evidence_ids"], product_role="crude intermediate")]
r["reaction_scale"] = f("(S)-(+)-carvone (2.0 mmol)", {"reference_molecule_id": "m_carvone", "amount": q(2.0, "mmol")}, r["evidence_ids"])
r["solvent_conditions"]["composition"] = f("absolute methanol plus methylamine solution in methanol", ["m_meoh"], r["evidence_ids"])
r["steps"] = [step("imine_formation", "reaction", ["carvone", "ti", "methylamine", "meoh"], r["evidence_ids"], "stirred", "r.t.", q(12, "h"), "nitrogen"),
              step("borohydride_reduction", "reaction", ["nabh4", "S6"], r["evidence_ids"], "NaBH4 added; mixture stirred", q(-78, "°C"), q(2, "h")),
              step("quench", "workup", ["water", "S6"], r["evidence_ids"], "water added; stirring maintained", "r.t.", q(20, "min")),
              step("filter_dry", "workup", ["celite", "meoh", "etoac", "na2so4", "S6"], r["evidence_ids"], "precipitate filtered through Celite and washed; combined organics dried and concentrated; crude S6 used directly")]
set_procedure(r, r["evidence_ids"], "(S)-(+)-carvone, Ti(OiPr)4 and methylamine in absolute MeOH, N2, rt 12 h; NaBH4 at -78 °C for 2 h; aqueous workup to crude S6.")
r["review"]["unresolved"].append("Crude S6 mass, amount and yield are unreported.")
REACTIONS.append(r)

r = reaction("preparation_1l_s6_to_1l", (E[5], I[5]), categories=["other"],
             location="Preparation of Allylamine 1l, second arrow", name="reductive methylation of crude S6",
             transformation="secondary amine S6 to tertiary allylamine 1l", class_mode="amination")
r["participants"] = [participant("S6", "m_S6", "reactant", r["evidence_ids"]),
                     participant("mecn", "m_mecn", "solvent", r["evidence_ids"], volume=q(3.0, "mL")),
                     participant("hcho", "m_formaldehyde", "reactant", r["evidence_ids"], volume=q(0.40, "mL")),
                     participant("nabh3cn", "m_nabh3cn", "reagent", r["evidence_ids"], mass=q(81.7, "mg"), amount=q(1.3, "mmol")),
                     participant("acoh", "m_acetic_acid", "additive", r["evidence_ids"]),
                     participant("naoh", "m_naoh", "workup_agent", r["evidence_ids"], stage="workup", volume=q(10, "mL")),
                     participant("etoac", "m_etoac", "workup_agent", r["evidence_ids"], stage="workup"),
                     participant("brine", "m_brine", "workup_agent", r["evidence_ids"], stage="workup", volume=q(20, "mL")),
                     participant("na2so4", "m_na2so4", "workup_agent", r["evidence_ids"], stage="workup"),
                     participant("dcm", "m_dcm", "purification_agent", r["evidence_ids"], stage="purification"),
                     participant("meoh_pur", "m_meoh", "purification_agent", r["evidence_ids"], stage="purification"),
                     participant("product", "m_1l", "product", r["evidence_ids"], product_role="isolated allylamine")]
r["steps"] = [step("reductive_methylation", "reaction", ["S6", "mecn", "hcho", "nabh3cn", "acoh", "product"], r["evidence_ids"], "formaldehyde then NaBH3CN added; a few drops of acetic acid maintained near-neutral pH; stirred", "r.t.", q(2, "h")),
              step("workup", "workup", ["naoh", "etoac", "brine", "na2so4", "product"], r["evidence_ids"], "concentrated; 2 N NaOH added; extracted with EtOAc 3 × 15 mL; washed, dried and concentrated"),
              step("chromatography", "purification", ["dcm", "meoh_pur", "product"], r["evidence_ids"], "column chromatography, dichloromethane/methanol 15:1")]
r["outcomes"] = [outcome("product_mass", "product_amount", ["m_1l"], r["evidence_ids"], q(115, "mg"), basis="isolated", scope={"product": "1l"}),
                 outcome("reported_yield", "yield", ["m_1l"], r["evidence_ids"], q(64, "%"), "reported isolated yield", "isolated", {"source_span": "procedure reports 64%; basis amount of crude S6 is unreported"}),
                 outcome("product_dr", "selectivity", ["m_1l"], r["evidence_ids"], "95:5", basis="1H NMR spectroscopic analysis", scope={"metric": "dr", "isomer_order": "major:minor"})]
set_procedure(r, r["evidence_ids"], "Crude S6 in MeCN was treated with 30% aqueous formaldehyde and NaBH3CN at near-neutral pH, rt 2 h; workup and DCM/MeOH chromatography gave 1l.")
r["review"]["unresolved"].append("The amount of crude S6 is unreported, so the denominator of the stated 64% yield cannot be reconstructed from this page.")
r["notes"].append(derived("115 mg 1l; formula C12H21N from HRMS; source reports 64% after starting from 2.0 mmol carvone in prior step",
                           {"observation": "115 mg corresponds to about 0.64 mmol by the reported formula; this is about 32% of 2.0 mmol carvone, not 64%",
                            "interpretation": "possible unstated denominator or source arithmetic inconsistency; no correction selected"},
                           r["evidence_ids"], "Calculation uses reported mass and molecular formula only to flag, not repair, the source values."))
REACTIONS.append(r)


# Source context: characterization is attached to synthesis records, not duplicated as reactions.
source_context = [
    {"context_id": f"{DOC}:page3_general_information", "document_id": DOC,
     "kind": "source_context", "evidence_ids": [E[3], I[3]],
     "description": "General NMR/HRMS/HPLC/optical-rotation methods and literature-preparation statements; not a synthesis experiment by itself."},
    {"context_id": f"{DOC}:pages3_4_allylamine_characterization", "document_id": DOC,
     "kind": "characterization", "evidence_ids": [E[3], I[3], E[4], I[4]],
     "description": "Product-specific NMR, HRMS, optical rotation and physical-property paragraphs for 1a-i and 1k; linked to preparation records, not separate reactions."},
    {"context_id": f"{DOC}:pages5_6_product_characterization", "document_id": DOC,
     "kind": "characterization", "evidence_ids": [E[5], I[5], E[6], I[6]],
     "description": "NMR, HRMS, HPLC retention data, optical rotation and melting points for 4a-o; not separate synthetic experiments."},
    {"context_id": f"{DOC}:products4_tautomer_statement", "document_id": DOC,
     "kind": "source_context", "evidence_ids": [E[5], I[5]],
     "description": "Authors state all products 4 exist predominantly in the enol form according to NMR spectra; named diketone structures are retained as registry candidates with this limitation."},
]


template = json.loads((RUN / "submitted.json").read_text())
bundle = copy.deepcopy(template)
bundle["molecules"] = list(MOLS.values())
bundle["reactions"] = REACTIONS
bundle["source_context"] = source_context

raw = copy.deepcopy(bundle)
raw["raw_model_metadata"] = {
    "batch_id": "allene-stereo-v31-200-20261017",
    "bounded_pages": [1, 2, 3, 4, 5, 6],
    "semantic_accuracy_verified": False,
    "note": "Machine-authored candidates after direct visual inspection; all formal records remain needs_review after import.",
}

all_ids = [r["reaction_id"] for r in REACTIONS]
page3_prep = [x for x in prep_ids if x.endswith(tuple(f"preparation_{v}" for v in ("1a", "1b", "1c", "1d")))]
page4_prep = [x for x in prep_ids if x not in page3_prep]
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title/abstract and graphical abstract", "finding": "Current catalyst-free aza-Claisen/hydrolysis generic reaction; repeated as a generic scheme, not a concrete 21-example aggregate experiment.", "evidence_ids": [E[1], I[1]]},
            {"region": "Scheme 1a", "finding": "Prior literature ring-expansion example with generic R groups; classified as literature_example.", "evidence_ids": [E[1], I[1]]},
            {"region": "Scheme 1b", "finding": "This-work generic reaction plus proposed addition/[3,3]/hydrolysis mechanism; separated into generic_scheme and proposed_mechanism.", "evidence_ids": [E[1], I[1]]},
            {"region": "lower-right prose", "finding": "Concrete rac-1a/2a dichloromethane experiment gives rac-3a in 85% NMR yield; merged with Table 1 entry 1.", "evidence_ids": [E[1], I[1], E[2], I[2]]}],
         "record_ids": [f"{DOC}:scheme1a_previous_work", f"{DOC}:generic_aza_claisen_hydrolysis", f"{DOC}:scheme1b_proposed_mechanism", table1_ids[0]],
         "unresolved": ["Generic Scheme 1 R groups and previous-work base/additive identities are not expanded."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Table 1 rows 1-14 and footnote", "finding": "Fourteen distinct solvent-screening experiments; crude-product 1H NMR yield basis retained.", "evidence_ids": [E[2], I[2]]},
            {"region": "Table 2 rows 1-8 and footnote", "finding": "Eight N-substituent experiments, including two true 0%/no-rearrangement entries; isolated-yield basis retained.", "evidence_ids": [E[2], I[2]]},
            {"region": "Table 3 rows 1-12 and footnote", "finding": "Twelve enantioenriched allylamine scope experiments; input/product ee or dr, isolated yields and methods retained.", "evidence_ids": [E[2], I[2]]}],
         "record_ids": table1_ids + table2_ids + table3_ids,
         "unresolved": ["Table 3 gives no numerical E/Z ratios; source wording 'exclusive E' is preserved qualitatively."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Table 4 rows 1-9 and footnote", "finding": "Nine allenone scope experiments; isolated yields, 97% ee, X-ray/analogy assignment and E selectivity retained.", "evidence_ids": [E[3], I[3]]},
            {"region": "Table 4 discussion below", "finding": "Two separate bulky allenones 2k and 2l reported not to react; no conversion/recovery invented.", "evidence_ids": [E[3], I[3]]},
            {"region": "Preparation of Allylamines 1a-i and 1k scheme/procedure", "finding": "Generic cyclization plus product-specific preparations 1a-d on this page; characterization remains source_context.", "evidence_ids": [E[3], I[3]]}],
         "record_ids": table4_ids + negative_ids + [prep_generic["reaction_id"]] + page3_prep,
         "unresolved": ["Allenones were prepared by literature procedures; no supplied SI or precursor procedures are associated."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "left-column product paragraphs", "finding": "Product-specific preparations/characterization for 1e-i and 1k; masses and isolated yields linked to the page-3 general procedure.", "evidence_ids": [E[4], I[4]]},
            {"region": "Preparation of Allylamine 1j scheme", "finding": "Two isolated-span records: S2 through crude aldehyde to S4 (88%), and four-step S4 through S5/crude alcohol to 1j (12% overall from S4).", "evidence_ids": [E[4], I[4]]},
            {"region": "right-column multistep procedure", "finding": "Low-temperature addition/hold, -15 °C quench, room-temperature holds, workups and chromatography retained as ordered stages.", "evidence_ids": [E[4], I[4]]}],
         "record_ids": page4_prep + [f"{DOC}:preparation_1j_s2_to_s4", f"{DOC}:preparation_1j_s4_to_1j"],
         "unresolved": ["S2/S3 drawings and crude intermediate quantities are not converted into fully verified concrete structures."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "Preparation of Allylamine 1l scheme/procedure", "finding": "Carvone-to-crude-S6 and S6-to-1l arrows recorded separately; -78 °C reduction, workup, final 64%/115 mg and 95:5 dr retained.", "evidence_ids": [E[5], I[5]]},
            {"region": "General Procedure for Screening of Solvents", "finding": "Table 1 addition order and concentration/removal details linked to all fourteen existing Table 1 records; no duplicates created.", "evidence_ids": [E[5], I[5]]},
            {"region": "General Procedure for Tables 2-4", "finding": "Reaction, wet-silica hydrolysis, filtration/wash and 10:1 chromatography stages inherited by Tables 2-4 records.", "evidence_ids": [E[5], I[5]]},
            {"region": "4a-4e characterization", "finding": "Mass/yield/ee paragraphs corroborate scope records; spectra and retention times retained as source_context rather than reactions.", "evidence_ids": [E[5], I[5]]}],
         "record_ids": [f"{DOC}:preparation_1l_carvone_to_s6", f"{DOC}:preparation_1l_s6_to_1l"] + table1_ids + table2_ids + table3_ids + table4_ids,
         "unresolved": ["Reported 1l mass/formula and prior 2.0 mmol scale imply about 32%, while text states 64%; denominator may be unstated, so no value was corrected."]},
        {"page": 6, "state": "attempted", "regions": [
            {"region": "4e continuation and 4f-4l characterization", "finding": "Product-specific masses/yields/ee or dr and HPLC/NMR methods linked to Table 3 entries; peak lists are source_context.", "evidence_ids": [E[6], I[6]]},
            {"region": "4m-4o characterization and 4p start", "finding": "Table 4 product data on this page linked to entries 1-3; 4p paragraph begins but continues outside the bounded page range.", "evidence_ids": [E[6], I[6]]}],
         "record_ids": table3_ids[4:] + table4_ids[:3],
         "unresolved": ["The 4p characterization paragraph is incomplete on page 6 and was not treated as fully covered; the Table 4 row remains the bounded source for its experiment."]},
    ],
}

mentions = st.source_mentions(json.loads((RUN / "dataset.json").read_text()),
                              ["rac-1a", "2a", "rac-3a", "4a", "1j", "1l", "4o"], context=180)
stereo_diagnostics = {}
for mid, mol in MOLS.items():
    smi = st.value(mol["smiles"])
    if smi:
        stereo_diagnostics[mid] = st.stereo_inspect(smi)

dump(RUN / "raw_model/pages_0001_0006_extraction.json", raw)
dump(RUN / "raw_model/pages_0001_0006_mentions.json", mentions)
dump(RUN / "raw_model/stereo_diagnostics_pages_0001_0006.json", stereo_diagnostics)
dump(RUN / "coverage/pages_0001_0006.json", coverage)
dump(RUN / "submitted.json", bundle)
print(json.dumps({"molecules": len(MOLS), "reactions": len(REACTIONS),
                  "source_context": len(source_context), "pages": 6}, ensure_ascii=False))
