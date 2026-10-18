#!/usr/bin/env python3
"""Build the source-specific v3.1 bundle for Chem. Commun. 2021, 57, 3599-3602.

This script is intentionally source-specific.  It consumes only the prepared dataset and
the clean standard-init template in this workspace.  It does not read historical records.
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

DOC = "doc_eaeb2de1d7d95d0f"
E = {p: f"{DOC}:page{p}:image" for p in range(1, 5)}
T = {p: f"{DOC}:page{p}:text" for p in range(1, 5)}
INIT = ROOT / "runs/current/raw_model/standard_init.json"
DATASET = ROOT / "runs/current/dataset.json"
OUT = ROOT / "runs/current/submitted.json"
RAW = ROOT / "runs/current/raw_model"


def write(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def q(number, unit, **extra):
    result = {"number": number, "unit": unit}
    result.update(extra)
    return result


def fact(raw=None, normalized=None, evidence=(), method="original_extraction",
         status=None, derivation=None, alternatives=()):
    if normalized is None and raw is not None and not alternatives:
        normalized = raw
    return st.fact(raw=raw, normalized=normalized, evidence=evidence, method=method,
                   status=status, derivation=derivation, alternatives=alternatives)


def inferred(normalized, evidence, derivation, raw=None):
    return fact(raw=raw, normalized=normalized, evidence=evidence, method="inference",
                derivation=derivation)


def missing(evidence=(), status="unreported"):
    return fact(evidence=evidence, status=status)


bundle = json.loads(INIT.read_text())
dataset = json.loads(DATASET.read_text())
assert not dataset.get("experiments") and "standard" not in dataset
assert bundle["documents"] == dataset["documents"]

molecules = {}
reactions = []


def molecule(mid, label, evidence, name=None, smiles=None, note=None,
             chirality="unresolved", geometry="unresolved", representation="not_encoded",
             atom_centers=None, double_bonds=None, axial=None, mixture=None,
             assignment=None):
    ev = list(dict.fromkeys(evidence))
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = fact(label, label, ev)
    original = name or label
    m["names"]["original"] = fact(original, original, ev)
    m["names"]["english"] = fact(original, original, ev)
    m["names"]["abbreviations"] = missing(ev)
    m["names"]["aliases"] = missing(ev)
    m["cas_numbers"] = missing(ev)
    if smiles:
        diag = st.stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = fact(
            raw=f"candidate transcription from source drawing for {label}",
            normalized=diag["canonical_isomeric_smiles"], evidence=ev,
            method="structure_derivation",
            derivation="Candidate graph transcribed from the displayed source structure; RDKit parsing/canonicalization does not verify source identity.")
        m["stereochemistry"]["isomeric_smiles"] = fact(
            raw=diag["canonical_isomeric_smiles"], normalized=diag["canonical_isomeric_smiles"],
            evidence=ev, method="structure_derivation",
            derivation="Canonical isomeric SMILES calculated from the unverified candidate graph.")
        m["stereochemistry"]["identity_key"] = fact(
            raw=diag["stereo_graph_key"], normalized=diag["stereo_graph_key"], evidence=ev,
            method="structure_derivation",
            derivation="Stereo graph key calculated from the unverified candidate graph.")
        m["stereochemistry"]["representation_status"] = inferred(
            representation, ev,
            "Status reflects source stereochemical content and candidate-graph/toolkit limits; parsing is not source verification.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing(ev)
        m["stereochemistry"]["isomeric_smiles"] = missing(ev)
        m["stereochemistry"]["identity_key"] = missing(ev)
        m["stereochemistry"]["representation_status"] = inferred(
            representation, ev,
            "A safe complete stereo-aware graph was not asserted from the compact source drawing or name.")
    for key, value in (("atom_centers", atom_centers), ("double_bonds", double_bonds),
                       ("axial", axial), ("mixture", mixture), ("assignment", assignment)):
        m["stereochemistry"][key] = fact(value, value, ev) if value is not None else missing(ev)
    m["stereochemistry"]["chirality_status"] = inferred(
        chirality, ev, "Source drawing/name was inspected; unassigned stereochemistry is retained as unassigned.")
    m["stereochemistry"]["geometry_status"] = inferred(
        geometry, ev, "Source drawing/name was inspected; unassigned E/Z or other geometry is not invented.")
    m["identity_note"] = fact(note, note, ev) if note else missing(ev)
    m["external_properties"] = []
    molecules[mid] = m
    return mid


def generic_molecule(mid, label, evidence, fixed, variables=(), note=None):
    ev = list(dict.fromkeys(evidence))
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = fact(label, label, ev)
    m["names"]["original"] = fact(label, label, ev)
    m["names"]["english"] = fact(label, label, ev)
    m["names"]["abbreviations"] = missing(ev)
    m["names"]["aliases"] = missing(ev)
    m["cas_numbers"] = missing(ev)
    m["representation_type"] = "generic_template"
    m["smiles"] = missing(ev)
    template = {"fixed_structure": fact(fixed, fixed, ev), "variables": [],
                "instance_mapping": [],
                "reason": fact("Source uses an R-group or generic reaction drawing.",
                               "Generic source drawing; no wildcard SMILES was asserted.", ev)}
    for var in variables:
        template["variables"].append({
            "label": var,
            "attachment_points": fact(f"source-labelled attachment for {var}", f"source-labelled attachment for {var}", ev),
            "definition": fact(f"{var} as drawn/defined in the source", f"{var} as drawn/defined in the source", ev),
            "allowed_values": missing(ev),
            "constraints": missing(ev),
        })
    m["generic_template"] = template
    for key in ("atom_centers", "double_bonds", "axial", "mixture", "assignment",
                "isomeric_smiles", "identity_key"):
        m["stereochemistry"][key] = missing(ev)
    m["stereochemistry"]["chirality_status"] = inferred(
        "unresolved", ev, "Generic R-group structure does not establish a unique stereochemical identity.")
    m["stereochemistry"]["geometry_status"] = inferred(
        "unresolved", ev, "Generic R-group structure does not establish unique alkene/allene geometry.")
    m["stereochemistry"]["representation_status"] = inferred(
        "not_encoded", ev, "Generic template is retained instead of an unsafe wildcard SMILES.")
    m["identity_note"] = fact(note, note, ev) if note else missing(ev)
    m["external_properties"] = []
    molecules[mid] = m
    return mid


def add_simple(mid, name, smiles=None, evidence=(E[2],)):
    return molecule(mid, name, evidence, name=name, smiles=smiles,
                    note="Identity is source-named; candidate graph, when present, is unverified.",
                    chirality="achiral" if smiles else "unresolved",
                    geometry="not_applicable" if smiles else "unresolved",
                    representation="encoded" if smiles else "not_encoded",
                    assignment={"basis": "candidate graph has no represented stereogenic element"} if smiles else None)


def charge(p, key, raw, normalized, evidence, method="original_extraction", derivation=None):
    p["charges"][key] = fact(raw, normalized, evidence, method=method, derivation=derivation)


def participant(pid, mid, role, evidence, stage="reaction", product_role=None):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    if product_role:
        p["product_role"] = fact(product_role, product_role, evidence)
    return p


def reaction(rid, record_type, evidence, categories=None, location=None, name=None,
             transformation=None, general=()):
    ev = list(dict.fromkeys(evidence))
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = fact(categories, categories, ev) if categories is not None else missing(ev)
    r["location"] = fact(location, location, ev) if location is not None else missing(ev)
    r["general_scheme_ids"] = list(general)
    r["reaction_name"] = fact(name, name, ev) if name else missing(ev)
    r["transformation_type"] = fact(transformation, transformation, ev) if transformation else missing(ev)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    return r


def set_reaction_scale(r, ref_mid, amount, evidence):
    r["reaction_scale"] = fact(
        raw={"reference_molecule_id": ref_mid, "amount": amount},
        normalized={"reference_molecule_id": ref_mid, "amount": amount}, evidence=evidence)


def set_solvent(r, mid, name, volume, evidence, method="original_extraction"):
    r["solvent_conditions"]["composition"] = fact(name, [mid], evidence, method=method,
                                                   derivation="Explicitly inherited source condition." if method == "inherited_general_condition" else None)
    if volume is not None:
        r["solvent_conditions"]["total_volume"] = fact(volume, volume, evidence, method=method,
                                                       derivation="Explicitly inherited source condition." if method == "inherited_general_condition" else None)


def add_step(r, sid, evidence, pids, operation=None, temp=None, time=None,
             atmosphere=None, stage="reaction"):
    s = st.new_step(sid, evidence, stage=stage)
    s["participant_ids"] = list(pids)
    if operation is not None:
        s["operation"] = fact(operation, operation, evidence)
    if temp is not None:
        s["temperature"] = fact(temp, temp, evidence)
    if time is not None:
        s["time"] = fact(time, time, evidence)
    if atmosphere is not None:
        s["atmosphere"] = fact(atmosphere, atmosphere, evidence)
    r["steps"].append(s)
    return s


def add_outcome(r, oid, kind, targets, raw, normalized, evidence, yield_type=None,
                basis=None, scope=None):
    o = st.new_outcome(oid, kind, targets, evidence)
    o["value"] = fact(raw, normalized, evidence)
    if yield_type is not None:
        o["yield_type"] = fact(yield_type, yield_type, evidence)
    if basis is not None:
        o["measurement_basis"] = fact(basis, basis, evidence)
    if scope is not None:
        o["scope"] = fact(scope, scope, evidence)
    r["outcomes"].append(o)
    return o


def classify(r, codes=None, evidence=None, coverage="partial", limitation=None):
    ev = list(evidence or r["evidence_ids"])
    if coverage in ("unresolved", "not_applicable"):
        r["classification"]["bond_changes"] = []
        if coverage == "not_applicable":
            r["classification"]["labels"] = inferred([], ev, "No observed transformation; intended-product bonds were not assigned.")
        else:
            r["classification"]["labels"] = missing(ev)
        r["classification"]["coverage"] = inferred(
            coverage, ev, "No source-safe bond-change assignment was made." if coverage == "unresolved" else "No observed chemical transformation.")
        r["classification"]["limitations"] = [limitation or "Bond changes are unresolved or not applicable."]
        r["classification"]["model_votes"] = []
        return
    specs = {"C-C": ("C", "C", 0, 1, "formed"),
             "C-O": ("C", "O", 0, 1, "formed"),
             "C-H": ("C", "H", 0, 1, "hydrogen_gain"),
             "DB-C-C": ("C", "C", 1, 2, "order_increased")}
    changes, labels = [], []
    for idx, code in enumerate(codes or [], 1):
        a, b, before, after, kind = specs[code]
        cid = f"b{idx}"
        ch = {"change_id": cid, "atoms": [{"atom_map": None, "element": a},
                                             {"atom_map": None, "element": b}],
              "before_order": before, "after_order": after, "change_type": kind}
        if kind == "hydrogen_gain":
            ch["hydrogen_count"] = 1
        changes.append(ch)
        lab = st.bond_label(a, b, after)
        lab["supporting_change_ids"] = [cid]
        labels.append(lab)
    derivation = "Class-level bond changes inferred from the displayed overall transformation; exact atoms and full mapping are unresolved."
    r["classification"]["bond_changes"] = [
        fact(raw="displayed overall transformation", normalized=c, evidence=ev,
             method="inference", derivation=derivation) for c in changes]
    r["classification"]["labels"] = fact(
        raw="displayed overall transformation", normalized=labels, evidence=ev,
        method="inference", derivation=derivation)
    r["classification"]["coverage"] = inferred("partial", ev, derivation)
    r["classification"]["limitations"] = [
        limitation or "Exact atom identities, bond multiplicities, hydrogen changes, coproducts, and full atom mapping remain unresolved."]
    r["classification"]["model_votes"] = []


# Core source structures and generic families.
GEN1 = generic_molecule(f"{DOC}:mol:generic_1", "generic cyclopropyl allenyl ketone 1", [E[1], E[2]],
                        "cyclopropylidene allene bearing R1 acyl and R2 substituent", ("R1", "R2"))
GEN2 = generic_molecule(f"{DOC}:mol:generic_2", "generic 4-vinyl-2,3-dioxopyrrolidine 2", [E[1], E[3]],
                        "4-arylidene-2,3-dioxopyrrolidine bearing R3 and N-benzyl R4 substituents", ("R3", "R4"))
GEN3 = generic_molecule(f"{DOC}:mol:generic_3", "generic tricyclic product 3", [E[1], E[2], E[3]],
                        "2-oxopyrrolidine-fused tricyclic oxo-IEDDA product bearing R1-R4", ("R1", "R2", "R3", "R4"))

M1 = {}
M1["1a"] = molecule(
    f"{DOC}:mol:1a", "1a", [E[1], E[2], E[3]],
    name="3-cyclopropylideneprop-2-en-1-one 1a (R1=4-MeC6H4, R2=Ph)",
    smiles="Cc1ccc(C(=O)C(c2ccccc2)=C=C2CC2)cc1",
    note="Candidate connectivity follows the source drawing. No axial descriptor is reported; the cyclopropylidene terminus has symmetry-related ring paths, while toolkit cumulene diagnostics remain limited.",
    chirality="unspecified", geometry="unspecified", representation="partial",
    axial={"fragment": "cyclopropylidene allene", "assignment": "unassigned", "source_descriptor": None},
    double_bonds={"allene": "source gives no axial/E/Z descriptor"})

scheme2_subs = {
    "1b": "R1=H, R2=Ph", "1c": "R1=2-MeC6H4, R2=Ph", "1d": "R1=3-MeC6H4, R2=Ph",
    "1e": "R1=4-MeOC6H4, R2=Ph", "1f": "R1=4-FC6H4, R2=Ph", "1g": "R1=4-ClC6H4, R2=Ph",
    "1h": "R1=4-BrC6H4, R2=Ph", "1i": "R1=n-Pr, R2=Ph", "1j": "R1=Bn, R2=Ph",
    "1k": "R1=1-naphthyl, R2=Ph", "1l": "R1=2-thiophenyl, R2=Ph", "1m": "R1=4-MeC6H4, R2=4-MeC6H4",
    "1n": "R1=4-MeC6H4, R2=2-MeOC6H4", "1o": "R1=4-MeC6H4, R2=3-MeOC6H4",
    "1p": "R1=4-MeC6H4, R2=4-MeOC6H4", "1q": "R1=4-MeC6H4, R2=4-PhC6H4",
}
for label, desc in scheme2_subs.items():
    M1[label] = molecule(f"{DOC}:mol:{label}", label, [E[2]],
                         name=f"cyclopropyl allenyl ketone {label} ({desc})",
                         note=f"Source-defined Scheme 2 substrate; {desc}. Complete stereo-aware SMILES not asserted.",
                         chirality="unspecified", geometry="unspecified", representation="not_encoded",
                         axial={"fragment": "cyclopropylidene allene", "assignment": "unassigned"})

M2 = {}
M2["2a"] = molecule(
    f"{DOC}:mol:2a", "2a", [E[1], E[2], E[3]],
    name="1-benzyl-4-benzylidenepyrrolidine-2,3-dione",
    smiles="O=C1N(Cc2ccccc2)CC(=Cc2ccccc2)C1=O",
    note="Candidate connectivity follows the source name/drawing; exocyclic alkene geometry is not assigned.",
    chirality="unspecified", geometry="unspecified", representation="partial",
    double_bonds={"exocyclic_benzylidene": "configuration unassigned"})
M2["2b"] = molecule(
    f"{DOC}:mol:2b", "2b", [E[3]],
    name="4-(4-methylbenzylidene)-1-benzylpyrrolidine-2,3-dione",
    smiles="O=C1N(Cc2ccccc2)CC(=Cc2ccc(C)cc2)C1=O",
    note="Candidate connectivity follows Scheme 3; exocyclic alkene geometry is not assigned.",
    chirality="unspecified", geometry="unspecified", representation="partial",
    double_bonds={"exocyclic_arylidene": "configuration unassigned"})
scheme3_subs = {
    "2c": "R3=4-MeO", "2d": "R3=4-F", "2e": "R3=2-F", "2f": "R3=4-Cl",
    "2g": "R3=4-Br", "2h": "R3=4-CF3", "2i": "R3=4-NO2", "2j": "R3=1-naphthyl",
    "2k": "R3=2-thiophenyl", "2l": "R3=Ph; R4=4-Cl", "2m": "R3=Ph; R4=4-Br",
    "2n": "R3=Ph; R4=3,5-diCl",
}
for label, desc in scheme3_subs.items():
    M2[label] = molecule(f"{DOC}:mol:{label}", label, [E[3]],
                         name=f"4-vinyl-2,3-dioxopyrrolidine {label} ({desc})",
                         note=f"Source-defined Scheme 3 substrate; {desc}. Exocyclic alkene geometry is unassigned.",
                         chirality="unspecified", geometry="unspecified", representation="not_encoded",
                         double_bonds={"exocyclic_arylidene": "configuration unassigned"})

M3 = {}
scheme2_products = {
    "3aa": "R1=4-Me, R2=Ph", "3ba": "R1=H, R2=Ph", "3ca": "R1=2-Me, R2=Ph",
    "3da": "R1=3-Me, R2=Ph", "3ea": "R1=4-OMe, R2=Ph", "3fa": "R1=4-F, R2=Ph",
    "3ga": "R1=4-Cl, R2=Ph", "3ha": "R1=4-Br, R2=Ph", "3ia": "R1=n-Pr, R2=Ph",
    "3ja": "R1=Bn, R2=Ph", "3ka": "R1=1-naphthyl, R2=Ph", "3la": "R1=2-thiophenyl, R2=Ph",
    "3ma": "R1=4-Me, R2=4-Me", "3na": "R1=4-Me, R2=2-OMe", "3oa": "R1=4-Me, R2=3-OMe",
    "3pa": "R1=4-Me, R2=4-OMe", "3qa": "R1=4-Me, R2=4-Ph",
}
for label, desc in scheme2_products.items():
    ev = [E[1], E[2]] if label == "3aa" else [E[2], E[3]] if label == "3la" else [E[2]]
    assignment = {"relative_stereochemistry": "wedge/dash configuration shown in source drawing",
                  "absolute_configuration": "unassigned"}
    note = f"Source-defined tricyclic product; {desc}."
    if label == "3aa":
        note += " Prose reports a single diastereomer."
    if label == "3la":
        note += " Structure assigned by X-ray diffraction; external CIF was not supplied."
    M3[label] = molecule(f"{DOC}:mol:{label}", label, ev, name=f"tricyclic product {label} ({desc})",
                         note=note, chirality="partially_specified", geometry="unresolved",
                         representation="not_encoded",
                         atom_centers={"relative_configuration": "shown by wedge/dash", "absolute": "unassigned"},
                         assignment=assignment,
                         mixture="single diastereomer" if label == "3aa" else None)

scheme3_products = {
    "3ab": "R3=4-Me", "3ac": "R3=4-OMe", "3ad": "R3=4-F", "3ae": "R3=2-F",
    "3af": "R3=4-Cl", "3ag": "R3=4-Br", "3ah": "R3=4-CF3", "3ai": "R3=4-NO2",
    "3aj": "R3=1-naphthyl", "3ak": "R3=2-thiophenyl", "3al": "R3=Ph; R4=4-Cl",
    "3am": "R3=Ph; R4=4-Br", "3an": "R3=Ph; R4=3,5-diCl",
}
for label, desc in scheme3_products.items():
    is_aj = label == "3aj"
    M3[label] = molecule(
        f"{DOC}:mol:{label}", label, [E[3]], name=f"tricyclic product {label} ({desc})",
        note=(f"Source-defined tricyclic product; {desc}." +
              (" Reported as a pair of rotamers with d.s.=3:1 due to steric bulkiness; atropisomer/rotamer identity remains unresolved." if is_aj else "")),
        chirality="mixture" if is_aj else "partially_specified", geometry="unresolved",
        representation="not_encoded",
        atom_centers={"relative_configuration": "shown by wedge/dash", "absolute": "unassigned"},
        axial={"hindered_rotation": "reported as rotamers; axis/configuration not assigned"} if is_aj else None,
        mixture={"type": "reported pair of rotamers", "ratio": "3:1", "metric": "d.s."} if is_aj else None,
        assignment={"relative_stereochemistry": "wedge/dash configuration shown", "absolute": "unassigned"})

D_MID = molecule(f"{DOC}:mol:intermediate_D", "D", [E[1], E[3]],
                 name="cyclobutane-fused furan intermediate D",
                 note="Intermediate identity is source-assigned by 1H NMR tracing; complete stereo-aware structure was not encoded.",
                 chirality="unspecified", geometry="unspecified", representation="not_encoded")
P4 = molecule(f"{DOC}:mol:4", "4", [E[3]], name="pyrrole-fused tetrahydropyran derivative 4",
              note="Compact Scheme 4 product; safe complete stereo-aware SMILES was not asserted.",
              chirality="partially_specified", geometry="unresolved", representation="not_encoded",
              atom_centers={"relative_configuration": "source drawing retained", "absolute": "unassigned"},
              assignment={"source": "Scheme 4 drawing", "absolute": "unassigned"})

# Named reagents, catalysts, solvents, and characterization standard.
named = {
    "AgSbF6": ("silver hexafluoroantimonate", None), "CuCl": ("copper(I) chloride", None),
    "IPrAuCl": ("IPrAuCl", None), "AgNTf2": ("silver bis(trifluoromethanesulfonyl)imide", None),
    "AgOTf": ("silver triflate", None), "AgBF4": ("silver tetrafluoroborate", None),
    "InOTf3": ("indium(III) triflate", None), "BF3Et2O": ("boron trifluoride diethyl etherate", None),
    "YbOTf3": ("ytterbium(III) triflate", None), "YOTf3": ("yttrium(III) triflate", None),
    "GdOTf3": ("gadolinium(III) triflate", None), "PdCl2": ("palladium(II) chloride", "[Cl-].[Cl-].[Pd+2]"),
    "4A_MS": ("4 Å molecular sieves", None), "DCM": ("dichloromethane", "ClCCl"),
    "DCE": ("1,2-dichloroethane", "ClCCCl"), "toluene": ("toluene", "Cc1ccccc1"),
    "THF": ("tetrahydrofuran", "C1CCOC1"), "MeCN": ("acetonitrile", "CC#N"),
    "Et2O": ("diethyl ether", "CCOCC"), "LiAlH4": ("lithium aluminium hydride", "[Li+].[AlH4-]"),
    "ZnCl2": ("zinc chloride", "[Zn+2].[Cl-].[Cl-]"),
    "internal_standard": ("1,3,5-trimethoxybenzene", "COc1cc(OC)cc(OC)c1"),
}
def named_evidence(key):
    if key in {"PdCl2", "Et2O", "LiAlH4", "ZnCl2"}:
        return [E[3]]
    if key in {"AgSbF6", "DCM", "4A_MS"}:
        return [E[2], E[3]]
    return [E[2]]


N = {key: add_simple(f"{DOC}:mol:{key}", name, smi, named_evidence(key))
     for key, (name, smi) in named.items()}
# RDKit exposes unspecified tetrahedral potential around the ionic aluminium-hydride
# candidate.  Preserve that limitation instead of declaring the representation complete.
molecules[N["LiAlH4"]]["stereochemistry"]["chirality_status"] = inferred(
    "unresolved", [E[3]], "Ionic candidate graph has toolkit-exposed unspecified stereochemical potential; no source assignment is available.")
molecules[N["LiAlH4"]]["stereochemistry"]["representation_status"] = inferred(
    "partial", [E[3]], "Candidate ionic graph is parseable but not a complete source-verified stereochemical representation.")

# Literature examples in Scheme 1, each retained as its own source-classified arrow/family.
lit_specs = [
    ("lit_i_a_13_migration", "Scheme 1 I-a", "metal-catalyzed cyclization/1,3-migration to indoles"),
    ("lit_i_b_12_migration", "Scheme 1 I-b", "metal-catalyzed cyclization/1,2-migration of cyclic amines/amides"),
    ("lit_ii_c_allenyl_ketone", "Scheme 1 II-c", "oxo-cyclization/1,2-migration of allenyl ketones"),
    ("lit_ii_d_vcp_ketone", "Scheme 1 II-d", "oxo-cyclization/1,2-migration of vinylidenecyclopropane ketones"),
    ("lit_iii_oxidative_ring_opening", "Scheme 1 III", "oxidative ring-opening of cyclobutane-fused furan D"),
    ("lit_iii_dimerization", "Scheme 1 III", "self-dimerization of cyclobutane-fused furan D"),
    ("lit_iii_gold_nucleophile", "Scheme 1 III", "gold-mediated ring addition/ring opening and nucleophile trapping"),
    ("lit_iii_gold_4plus3", "Scheme 1 III", "gold-mediated [4+3] cycloaddition"),
    ("lit_iii_gold_8plus4", "Scheme 1 III", "gold-mediated [8+4] cycloaddition"),
]
for rid, loc, desc in lit_specs:
    a = generic_molecule(f"{DOC}:mol:{rid}:reactant", f"generic reactant for {loc}", [E[1]], desc + " reactant family")
    p = generic_molecule(f"{DOC}:mol:{rid}:product", f"generic product for {loc}", [E[1]], desc + " product family")
    r = reaction(rid, "literature_example", [E[1], T[1]], location=loc,
                 name=desc, transformation=desc)
    r["participants"] = [participant("reactant", a, "reactant", [E[1]]),
                         participant("product", p, "product", [E[1]], product_role="reported literature product family")]
    r["notes"].append(fact("Cited prior work, not a current-paper experiment.",
                           "Cited prior work; excluded from current experiment count.", [E[1], T[1]]))
    r["review"]["unresolved"] = ["Generic literature structures and conditions are not expanded from the cited external papers."]
    classify(r, coverage="unresolved", limitation="Cited generic family; exact structures and bond mapping were not imported from external papers.")
    reactions.append(r)

GENERAL_ID = f"{DOC}:generic_domino_iedda"
gr = reaction("generic_domino_iedda", "generic_scheme", [E[1], E[2], E[3], T[4]],
              location="Scheme 1 IV; repeated generic headers in Schemes 2 and 3",
              name="silver-catalyzed domino cyclization-migration/inverse electron-demand oxo-Diels-Alder reaction",
              transformation="one-pot cyclization-migration/oxo-IEDDA")
gr["participants"] = [participant("generic_1", GEN1, "reactant", [E[1]]),
                      participant("generic_2", GEN2, "reactant", [E[1]]),
                      participant("generic_3", GEN3, "product", [E[1]], product_role="generic desired product")]
gr["notes"].append(fact("This work", "Current-paper generic transformation; not a concrete experiment.", [E[1]]))
gr["review"]["unresolved"] = ["R-group identities, exact atom mapping, and generic-scheme conditions are represented by specific experiments separately."]
classify(gr, ["C-C", "C-O"], limitation="Class-level overall cycloaddition/cyclization labels only; exact multiplicity and atom mapping are unresolved.")
reactions.append(gr)


def common_current_participants(r, sub1, sub2, product, catalyst, solvent, ev,
                                a1=q(0.20, "mmol"), a2=q(0.22, "mmol"),
                                cat=q(5, "mol%"), ms=q(100, "mg"), volume=q(2.0, "mL"),
                                inherited=False):
    method = "inherited_general_condition" if inherited else "original_extraction"
    deriv = "Explicitly inherited from the scheme/table common-condition footnote." if inherited else None
    p1 = participant("substrate_1", sub1, "reactant", ev)
    p2 = participant("substrate_2", sub2, "reactant", ev)
    pp = participant("product", product, "product", ev, product_role="desired tricyclic product")
    pc = participant("catalyst", catalyst, "catalyst", ev)
    pm = participant("molecular_sieves", N["4A_MS"], "additive", ev)
    ps = participant("solvent", solvent, "solvent", ev)
    charge(p1, "amount", a1, a1, ev, method, deriv)
    charge(p2, "amount", a2, a2, ev, method, deriv)
    charge(pc, "mol_percent", cat, cat, ev, method, deriv)
    pc["equivalence_reference"] = fact("relative to substrate 1", sub1, ev, method=method, derivation=deriv)
    charge(pm, "mass", ms, ms, ev, method, deriv)
    charge(ps, "volume", volume, volume, ev, method, deriv)
    r["participants"] = [p1, p2, pc, pm, ps, pp]
    set_reaction_scale(r, sub1, a1, ev)
    set_solvent(r, solvent, molecules[solvent]["names"]["original"]["normalized_value"], volume, ev, method)


# Table 1: one row per experiment. Entry 1 also represents the duplicate Scheme 2 3aa result.
table = [
    (1, "AgSbF6", "DCM", 94), (2, "CuCl", "DCM", 88), (3, "IPrAuCl", "DCM", None),
    (4, "AgNTf2", "DCM", 85), (5, "AgOTf", "DCM", 91), (6, "AgBF4", "DCM", 90),
    (7, "InOTf3", "DCM", None), (8, "BF3Et2O", "DCM", None), (9, "YbOTf3", "DCM", None),
    (10, "YOTf3", "DCM", None), (11, "GdOTf3", "DCM", None), (12, "AgSbF6", "DCE", 71),
    (13, "AgSbF6", "toluene", 61), (14, "AgSbF6", "THF", 78), (15, "AgSbF6", "MeCN", 70),
    (16, None, "DCM", None),
]
for entry, catalyst, solvent, yld in table:
    evidence = [E[2], T[2]]
    categories = ["optimization", "substrate_scope"] if entry == 1 else ["optimization"]
    loc = f"Table 1 entry {entry}" + ("; duplicate appearance as Scheme 2 product 3aa" if entry == 1 else "")
    r = reaction(f"table1_entry{entry}", "specific_experiment", evidence, categories, loc,
                 "optimization of Ag-catalyzed domino oxo-IEDDA", "domino cyclization-migration/oxo-IEDDA",
                 [GENERAL_ID])
    p1 = participant("substrate_1a", M1["1a"], "reactant", evidence)
    p2 = participant("substrate_2a", M2["2a"], "reactant", evidence)
    pp = participant("product_3aa", M3["3aa"], "product", evidence, product_role="desired product")
    pm = participant("molecular_sieves", N["4A_MS"], "additive", evidence)
    ps = participant("solvent", N[solvent], "solvent", evidence)
    charge(p1, "amount", q(0.20, "mmol"), q(0.20, "mmol"), evidence)
    charge(p2, "amount", q(0.22, "mmol"), q(0.22, "mmol"), evidence)
    charge(pm, "mass", q(100, "mg"), q(100, "mg"), evidence)
    charge(ps, "volume", q(2.0, "mL"), q(2.0, "mL"), evidence)
    plist = [p1, p2]
    if catalyst:
        pc = participant("catalyst", N[catalyst], "catalyst", evidence)
        charge(pc, "mol_percent", q(5, "mol%"), q(5, "mol%"), evidence)
        pc["equivalence_reference"] = fact("relative to 1a", M1["1a"], evidence)
        plist.append(pc)
    plist += [pm, ps]
    standard = participant("nmr_internal_standard", N["internal_standard"], "reagent", evidence,
                           stage="characterization")
    plist += [standard, pp]
    r["participants"] = plist
    set_reaction_scale(r, M1["1a"], q(0.20, "mmol"), evidence)
    set_solvent(r, N[solvent], molecules[N[solvent]]["names"]["original"]["normalized_value"], q(2.0, "mL"), evidence)
    add_step(r, "reaction", evidence, [p["participant_id"] for p in plist if p["stage"] == "reaction"],
             operation="combine under the Table 1 conditions", temp="r.t.", time=q(10, "min"), atmosphere="Ar")
    if yld is not None:
        add_outcome(r, "nmr_yield", "yield", [M3["3aa"]], f"{yld}%", q(yld, "%"), evidence,
                    yield_type="NMR yield",
                    basis="1H NMR using 1,3,5-trimethoxybenzene as internal standard")
        if entry == 1:
            add_outcome(r, "isolated_yield", "yield", [M3["3aa"]], "(90%) isolated", q(90, "%"), evidence,
                        yield_type="isolated yield", basis="isolated yield")
            r["notes"].append(fact(
                "Table 1 entry 1 and Scheme 2 3aa both report the 90% isolated optimized experiment; Scheme 2 gives 25 °C and a collective 5-60 min range, while Table 1 gives the exact 10 min run.",
                "Duplicate sources merged; exact Table 1 time retained and the compatible Scheme 2 range is not treated as a second experiment.", evidence))
        classify(r, ["C-C", "C-O"])
    else:
        wording = "no reaction occurred" if entry == 16 else "none of 3aa was formed" if entry == 3 else "ineffective in delivering 3aa"
        add_outcome(r, "desired_product_not_formed", "no_reaction", [M3["3aa"]],
                    "table dash (—); " + wording, wording, evidence,
                    scope={"desired_product": "3aa", "table_marker": "—"})
        classify(r, coverage="not_applicable", limitation="No desired transformation was observed; intended-product bonds were not assigned.")
    r["review"]["unresolved"] = ["SI procedure and full workup/purification are unavailable."]
    reactions.append(r)

# Scheme 2 additional substrate-scope entries.  3aa is already merged with Table 1 entry 1.
scheme2_rows = [
    ("3ba", "1b", 92), ("3ca", "1c", 83), ("3da", "1d", 72), ("3ea", "1e", 98),
    ("3fa", "1f", 73), ("3ga", "1g", 83), ("3ha", "1h", 83), ("3ia", "1i", 65),
    ("3ja", "1j", 76), ("3ka", "1k", 87), ("3la", "1l", 50), ("3ma", "1m", 81),
    ("3na", "1n", 53), ("3oa", "1o", 71), ("3pa", "1p", 93), ("3qa", "1q", 45),
]
for product, sub, yld in scheme2_rows:
    ev = [E[2], T[2]] + ([E[3], T[3]] if product == "3la" else [])
    r = reaction(f"scheme2_{product}", "specific_experiment", ev, ["substrate_scope"],
                 f"Scheme 2 product {product}",
                 "Ag-catalyzed domino cyclization-migration/inverse electron-demand oxo-Diels-Alder reaction",
                 "domino cyclization-migration/oxo-IEDDA", [GENERAL_ID])
    common_current_participants(r, M1[sub], M2["2a"], M3[product], N["AgSbF6"], N["DCM"], ev, inherited=True)
    add_step(r, "reaction", ev, [p["participant_id"] for p in r["participants"]],
             operation="react under Scheme 2 common conditions", temp=q(25, "°C"),
             time={"minimum": 5, "maximum": 60, "unit": "min", "raw_text": "5 to 60 minutes"}, atmosphere="Ar")
    add_outcome(r, "isolated_yield", "yield", [M3[product]], f"{yld}%", q(yld, "%"), ev,
                yield_type="isolated yield", basis="isolated yield")
    if product == "3la":
        r["notes"].append(fact("The structure of 3la was determined by X-ray diffraction; CIF is in unavailable ESI.",
                               "X-ray assignment reported; external CIF not supplied.", [E[3], T[3]]))
    r["review"]["unresolved"] = ["The scheme reports only a collective 5-60 min range, not the exact time for this product.",
                                  "Workup and purification details are in unavailable ESI."]
    classify(r, ["C-C", "C-O"])
    reactions.append(r)

# Scheme 3 substrate scope.
scheme3_rows = [
    ("3ab", "2b", 82), ("3ac", "2c", 89), ("3ad", "2d", 74), ("3ae", "2e", 45),
    ("3af", "2f", 60), ("3ag", "2g", 83), ("3ah", "2h", 81), ("3ai", "2i", 26),
    ("3aj", "2j", 70), ("3ak", "2k", 65), ("3al", "2l", 78), ("3am", "2m", 68),
    ("3an", "2n", 74),
]
for product, sub, yld in scheme3_rows:
    ev = [E[3], T[3]]
    r = reaction(f"scheme3_{product}", "specific_experiment", ev, ["substrate_scope"],
                 f"Scheme 3 product {product}",
                 "Ag-catalyzed domino cyclization-migration/inverse electron-demand oxo-Diels-Alder reaction",
                 "domino cyclization-migration/oxo-IEDDA", [GENERAL_ID])
    common_current_participants(r, M1["1a"], M2[sub], M3[product], N["AgSbF6"], N["DCM"], ev, inherited=True)
    add_step(r, "reaction", ev, [p["participant_id"] for p in r["participants"]],
             operation="react under Scheme 3 common conditions", temp=q(25, "°C"),
             time={"minimum": 5, "maximum": 60, "unit": "min", "raw_text": "5 to 60 min"}, atmosphere="Ar")
    add_outcome(r, "isolated_yield", "yield", [M3[product]], f"{yld}%", q(yld, "%"), ev,
                yield_type="isolated yield", basis="isolated yield")
    if product == "3aj":
        add_outcome(r, "rotamer_ds", "selectivity", [M3[product]], "d.s. = 3:1", "3:1", ev,
                    basis=None, scope={"metric": "d.s.", "order": "reported rotamer pair; component identities unassigned"})
    if product == "3ai":
        r["notes"].append(fact("Poor solubility of 2i might lead to the lowest yield.",
                               "Author-proposed explanation: poor solubility of 2i.", ev))
    r["review"]["unresolved"] = ["The scheme reports only a collective 5-60 min range, not the exact time for this product.",
                                  "Workup and purification details are in unavailable ESI."]
    if product == "3aj":
        r["review"]["unresolved"].append("The d.s. measurement method and identities/order of the two reported rotamers are not stated.")
    classify(r, ["C-C", "C-O"])
    reactions.append(r)

# Scheme 4: scale-up and derivatization.
ev3 = [E[3], T[3]]
r = reaction("scheme4_scaleup_3ab", "specific_experiment", ev3, ["scale_up"], "Scheme 4 I",
             "gram-scale synthesis of 3ab", "domino cyclization-migration/oxo-IEDDA", [GENERAL_ID])
p1 = participant("substrate_1a", M1["1a"], "reactant", ev3)
p2 = participant("substrate_2b", M2["2b"], "reactant", ev3)
pc = participant("catalyst", N["AgSbF6"], "catalyst", ev3)
pm = participant("molecular_sieves", N["4A_MS"], "additive", ev3)
ps = participant("solvent", N["DCM"], "solvent", ev3)
pp = participant("product_3ab", M3["3ab"], "product", ev3, product_role="desired product")
charge(p1, "amount", q(3.00, "mmol"), q(3.00, "mmol"), ev3); charge(p1, "mass", q(0.780, "g"), q(0.780, "g"), ev3)
charge(p2, "amount", q(3.15, "mmol"), q(3.15, "mmol"), ev3); charge(p2, "mass", q(0.917, "g"), q(0.917, "g"), ev3)
charge(pc, "mol_percent", q(5, "mol%"), q(5, "mol%"), ev3, "inherited_general_condition", "Inherited from the referenced standard conditions.")
pc["equivalence_reference"] = fact("relative to 1a", M1["1a"], ev3, method="inherited_general_condition", derivation="Inherited from standard conditions.")
r["participants"] = [p1, p2, pc, pm, ps, pp]
set_reaction_scale(r, M1["1a"], q(3.00, "mmol"), ev3)
set_solvent(r, N["DCM"], "dichloromethane", None, ev3, "inherited_general_condition")
add_step(r, "reaction", ev3, [p["participant_id"] for p in r["participants"]],
         operation="react under standard conditions", temp=q(25, "°C"),
         time={"minimum": 5, "maximum": 60, "unit": "min", "raw_text": "standard conditions; scheme range 5-60 min"}, atmosphere="Ar")
add_outcome(r, "product_amount_mmol", "product_amount", [M3["3ab"]], "2.36 mmol", q(2.36, "mmol"), ev3)
add_outcome(r, "product_amount_mass", "product_amount", [M3["3ab"]], "1.299 g", q(1.299, "g"), ev3)
add_outcome(r, "yield", "yield", [M3["3ab"]], "79%", q(79, "%"), ev3)
r["notes"].append(fact("standard conditions", "Standard-condition loading/temperature/atmosphere inherited; absolute solvent and molecular-sieve amounts at scale are unreported.", ev3))
r["review"]["unresolved"] = ["Yield measurement basis is not stated in Scheme 4.",
                              "Exact scale-up solvent volume, molecular-sieve mass, exact time, workup, and purification are unreported."]
classify(r, ["C-C", "C-O"])
reactions.append(r)

r = reaction("scheme4_derivatization_4", "specific_experiment", ev3, ["derivatization"], "Scheme 4 II",
             "cascade reduction and dehydration-aromatization of 3ab", "reduction/dehydration-aromatization")
p0 = participant("substrate_3ab", M3["3ab"], "reactant", ev3)
pl = participant("LiAlH4", N["LiAlH4"], "reagent", ev3)
pz = participant("ZnCl2", N["ZnCl2"], "reagent", ev3)
pe = participant("ether", N["Et2O"], "solvent", ev3)
pp = participant("product_4", P4, "product", ev3, product_role="pyrrole-fused tetrahydropyran derivative")
r["participants"] = [p0, pl, pz, pe, pp]
set_solvent(r, N["Et2O"], "diethyl ether", None, ev3)
add_step(r, "reaction", ev3, [p["participant_id"] for p in r["participants"]],
         operation="treat 3ab with LiAlH4 and ZnCl2 in ether", temp=q(0, "°C"))
add_outcome(r, "yield", "yield", [P4], "50%", q(50, "%"), ev3)
r["review"]["unresolved"] = ["Reagent amounts, time, atmosphere, workup, purification, and yield measurement basis are unreported.",
                              "Complex net bond changes and complete reactant/product structures are not safely encoded."]
classify(r, coverage="unresolved", limitation="Reduction/dehydration-aromatization is source-assigned, but exact net bond changes are not safely mapped.")
reactions.append(r)

# Figure 1 NMR tracing experiment.
r = reaction("figure1_nmr_tracing", "specific_experiment", ev3, ["mechanistic_study"], "Figure 1",
             "1H NMR tracing experiment", "sequential formation and trapping of cyclobutane-fused furan D", [GENERAL_ID])
p1 = participant("substrate_1a", M1["1a"], "reactant", ev3)
pc = participant("silver_catalyst", N["AgSbF6"], "catalyst", ev3)
p2 = participant("substrate_2a", M2["2a"], "reactant", ev3)
pd = participant("intermediate_D", D_MID, "product", ev3, product_role="observed intermediate")
pp = participant("product_3aa", M3["3aa"], "product", ev3, product_role="observed final product")
r["participants"] = [p1, pc, p2, pd, pp]
add_step(r, "initial_spectrum", ev3, ["substrate_1a"], operation="record spectrum of only substrate 1a")
add_step(r, "after_catalyst", ev3, ["substrate_1a", "silver_catalyst", "intermediate_D"],
         operation="add catalyst [Ag] and monitor formation of D", time=q(5, "min", qualifier="approx", raw_text="about 5 min"))
add_step(r, "after_2a", ev3, ["intermediate_D", "substrate_2a", "product_3aa"],
         operation="add substrate 2a; D signals disappear and 3aa signals appear", time=q(5, "min", qualifier="approx", raw_text="about 5 min"))
add_outcome(r, "D_observed", "other", [D_MID], "signals of intermediate D appeared", "D observed", ev3,
            basis="1H NMR spectral comparison", scope={"stage": "after catalyst addition"})
add_outcome(r, "product_observed", "other", [M3["3aa"]], "D signals immediately disappeared and 3aa signals appeared after 2a", "3aa observed", ev3,
            basis="1H NMR tracing", scope={"stage": "after 2a addition"})
r["review"]["unresolved"] = ["The figure does not report exact sample quantities, solvent, temperature, atmosphere, or catalyst identity beyond [Ag] in the tracing labels."]
classify(r, ["C-C", "C-O"], limitation="Observed sequential pathway supports class-level C-C/C-O formation; exact stepwise atom mapping is unresolved.")
reactions.append(r)

# Scheme 5 first arrow: independent preparation/observation of D with PdCl2.
r = reaction("scheme5_prepare_D", "specific_experiment", ev3, ["mechanistic_study", "control"], "Scheme 5, preparation of D",
             "preparation of cyclobutane-fused furan D", "PdCl2-promoted cyclization-migration")
p1 = participant("substrate_1a", M1["1a"], "reactant", ev3)
pc = participant("PdCl2", N["PdCl2"], "catalyst", ev3)
ps = participant("DCM", N["DCM"], "solvent", ev3)
pd = participant("intermediate_D", D_MID, "product", ev3, product_role="observed intermediate")
r["participants"] = [p1, pc, ps, pd]
set_solvent(r, N["DCM"], "dichloromethane", None, ev3)
add_step(r, "reaction", ev3, [p["participant_id"] for p in r["participants"]],
         operation="treat 1a with PdCl2 in DCM", temp="r.t.")
add_outcome(r, "D_identity", "other", [D_MID], "intermediate D determined by 1H NMR spectrum", "D observed", ev3,
            basis="1H NMR spectrum")
r["review"]["unresolved"] = ["PdCl2 loading, substrate scale, time, atmosphere, conversion, isolation, and workup are unreported."]
classify(r, ["C-C", "C-O"], limitation="Class-level cyclization/migration labels only; exact bond multiplicity and atom mapping are unresolved.")
reactions.append(r)

# Scheme 5 two independent D-trapping controls.
for suffix, with_ag, yld in (("with_AgSbF6", True, 63), ("without_AgSbF6", False, 60)):
    r = reaction(f"scheme5_D_to_3aa_{suffix}", "specific_experiment", ev3, ["mechanistic_study", "control"],
                 f"Scheme 5 {'i' if with_ag else 'ii'}",
                 "control conversion of preformed D to 3aa", "oxo-IEDDA of D with 2a")
    pd = participant("intermediate_D", D_MID, "reactant", ev3)
    p2 = participant("substrate_2a", M2["2a"], "reactant", ev3)
    pp = participant("product_3aa", M3["3aa"], "product", ev3, product_role="desired product")
    plist = [pd, p2]
    if with_ag:
        plist.append(participant("AgSbF6", N["AgSbF6"], "catalyst", ev3))
    plist.append(pp)
    r["participants"] = plist
    add_step(r, "reaction", ev3, [p["participant_id"] for p in plist],
             operation="combine preformed D and 2a " + ("with AgSbF6" if with_ag else "without AgSbF6"))
    add_outcome(r, "yield", "yield", [M3["3aa"]], f"{yld}%", q(yld, "%"), ev3)
    r["review"]["unresolved"] = ["Amounts, solvent, temperature, time, atmosphere, workup, purification, and yield measurement basis are unreported."]
    classify(r, ["C-C", "C-O"], limitation="Class-level oxo-IEDDA labels only; exact atoms and mapping are unresolved.")
    reactions.append(r)

# Scheme 6 proposed mechanism, with page-4 prose continuation.
r = reaction("scheme6_proposed_mechanism", "proposed_mechanism", [E[3], T[3], E[4], T[4]],
             location="Scheme 6 and page-4 continuation",
             name="proposed Ag(I)-catalyzed cyclization-migration followed by catalyst-free oxo-IEDDA",
             transformation="proposed A-B-C-D sequence and D+2 to 3")
r["participants"] = [participant("generic_1", GEN1, "reactant", [E[3]]),
                     participant("generic_2", GEN2, "reactant", [E[3]]),
                     participant("generic_3", GEN3, "product", [E[3]], product_role="proposed final product")]
r["notes"].append(fact(
    "Carbonyl oxygen attacks Ag(I)-coordinated allene to A; A isomerizes to silver carbenoid B; ring expansion and metal elimination give D; D undergoes oxo-IEDDA with 2 to give 3.",
    "Author-proposed mechanism; intermediates A-C/D are not treated as executed isolated experiments.", [E[3], E[4], T[4]]))
r["review"]["unresolved"] = ["Mechanistic arrows are a proposal, not complete experimental atom mapping."]
classify(r, ["C-C", "C-O"], limitation="Overall class-level inferred labels only; proposed intermediates are not independently verified experiments.")
reactions.append(r)

# Source context: non-synthetic observations, metadata, and negative inventories.
bundle["source_context"] = [
    {"context_id": f"{DOC}:context:metadata", "document_id": DOC,
     "evidence_ids": [E[1], T[1]], "kind": "document_metadata",
     "content": {"title": bundle["documents"][0]["title"], "doi": "10.1039/d1cc00707f",
                 "si_status": "advertised by paper but not supplied or associated",
                 "ccdc": "2052391 advertised externally; CIF not supplied"}},
    {"context_id": f"{DOC}:context:scheme1_background", "document_id": DOC,
     "evidence_ids": [E[1], T[1]], "kind": "literature_context",
     "content": "Scheme 1 prior-work arrows are retained as literature_example records and are not current-paper experiments."},
    {"context_id": f"{DOC}:context:characterization", "document_id": DOC,
     "evidence_ids": [E[3], T[3]], "kind": "characterization",
     "content": "Figure 1 NMR spectra and the reported X-ray assignment for 3la support mechanistic/identity context; no standalone spectrum peak table was converted into a synthesis."},
    {"context_id": f"{DOC}:context:page4_references", "document_id": DOC,
     "evidence_ids": [E[4], T[4]], "kind": "non_experimental_page",
     "content": "Page 4 contains the mechanism prose continuation, conclusion, acknowledgements, conflicts, and references; it has no new current-paper reaction arrow or table row."},
    {"context_id": f"{DOC}:context:negative_spectral_inventory", "document_id": DOC,
     "evidence_ids": [E[1], E[2], E[3], E[4]], "kind": "negative_inventory",
     "content": "No standalone mass-spectrum peak list, molecular-weight table, or fragment-relative-abundance block occurs on pages 1-4."},
]

bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
write(OUT, bundle)

# Raw extraction inventory and literal mention hints remain separate from normalized records.
page_ids = {p: [] for p in range(1, 5)}
for r in reactions:
    pages = {int(e.split(":page", 1)[1].split(":", 1)[0]) for e in r["evidence_ids"] if ":page" in e}
    for p in pages:
        page_ids[p].append(r["reaction_id"])

inventory = {
    "document_id": DOC,
    "source": "visual inspection of rendered source pages; extracted text used as transcription aid",
    "pages": [
        {"page": 1, "image_evidence_id": E[1], "inventory": [
            "title/abstract/metadata", "Scheme 1: nine prior-work arrows/families", "Scheme 1 IV current generic scheme", "ESI/CCDC availability notice"]},
        {"page": 2, "image_evidence_id": E[2], "inventory": [
            "Table 1 entries 1-16", "Scheme 2 products 3aa-3qa", "3aa duplicate of Table 1 entry 1", "Scheme 2 footnotes including 5-60 min and isolated-yield basis"]},
        {"page": 3, "image_evidence_id": E[3], "inventory": [
            "Scheme 3 products 3ab-3an", "Scheme 4 scale-up and derivatization", "Figure 1 NMR tracing", "Scheme 5 preparation of D plus two controls", "Scheme 6 proposed mechanism"]},
        {"page": 4, "image_evidence_id": E[4], "inventory": [
            "Scheme 6 prose continuation", "conclusion/acknowledgements/conflicts/references", "no new experiment arrow/table row"]},
    ],
    "reaction_ids_by_page": page_ids,
    "counts": {"reactions": len(reactions), "molecules": len(molecules),
               "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in reactions),
               "literature_examples": sum(r["record_type"] == "literature_example" for r in reactions)},
}
write(RAW / "page_extraction.json", inventory)
write(RAW / "source_mentions.json", st.source_mentions(dataset, ["1a", "2a", "3aa", "3ab", "intermediate D"], context=180))

stereo = {}
for mid, m in molecules.items():
    smi = st.value(m["smiles"])
    if smi:
        stereo[mid] = st.stereo_inspect(smi)
write(RAW / "stereo_checks.json", stereo)

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title, abstract, DOI, ESI notice", "finding": "Document metadata; ESI and CCDC are advertised but not supplied.", "evidence_ids": [E[1], T[1]]},
            {"region": "Scheme 1 panels I-III", "finding": "Nine cited prior-work arrow/family records; classified as literature examples, not current experiments.", "evidence_ids": [E[1], T[1]]},
            {"region": "Scheme 1 panel IV", "finding": "Current generic domino cyclization-migration/oxo-IEDDA scheme; no concrete yield assigned here.", "evidence_ids": [E[1]]}],
         "record_ids": page_ids[1],
         "unresolved": ["Generic prior-work structures/conditions require external cited papers and were not expanded.", "Advertised SI/CCDC files were not supplied."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Table 1", "finding": "All 16 optimization rows recorded separately; entry 1 includes 94% NMR and 90% isolated yields; dash/no-product rows preserved.", "evidence_ids": [E[2], T[2]]},
            {"region": "Scheme 2", "finding": "All products 3aa-3qa inventoried; 3aa merged with matching Table 1 entry 1 and 3ba-3qa recorded independently.", "evidence_ids": [E[2], T[2]]},
            {"region": "Scheme 2 footnotes", "finding": "Common 0.20/0.22 mmol charges, 5 mol% AgSbF6, 100 mg 4A MS, 2.0 mL anhydrous DCM, 25 C, Ar, 5-60 min, isolated yields retained.", "evidence_ids": [E[2]]}],
         "record_ids": page_ids[2],
         "unresolved": ["Scheme 2 gives only a collective 5-60 min range, not exact product-specific times.", "Workup/purification are in unavailable ESI."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Scheme 3", "finding": "All 13 outcomes 3ab-3an recorded; 3aj retains reported rotamer d.s.=3:1 with unassigned component identities/method.", "evidence_ids": [E[3], T[3]]},
            {"region": "Scheme 4", "finding": "Gram-scale 3ab experiment and separate LiAlH4/ZnCl2 derivatization to 4 recorded.", "evidence_ids": [E[3], T[3]]},
            {"region": "Figure 1", "finding": "1H NMR tracing retained as a mechanistic experiment with two approximately 5 min observation stages.", "evidence_ids": [E[3], T[3]]},
            {"region": "Scheme 5", "finding": "PdCl2/DCM preparation/observation of D and the separate with/without AgSbF6 63%/60% controls recorded as three experiments.", "evidence_ids": [E[3], T[3]]},
            {"region": "Scheme 6", "finding": "Author-proposed mechanism retained as one proposed_mechanism record, not executed intermediate experiments.", "evidence_ids": [E[3], T[3]]}],
         "record_ids": page_ids[3],
         "unresolved": ["Scheme 4/5 yields lack stated measurement bases where explicitly noted in records.", "Many compact tricyclic structures and exact atom mappings remain unencoded.", "Figure 1 exact sample quantities/solvent/temperature are unreported."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "left-column top", "finding": "Continuation of Scheme 6 mechanism prose; linked to the existing proposed-mechanism record.", "evidence_ids": [E[4], T[4]]},
            {"region": "conclusion and acknowledgements", "finding": "Narrative summary and funding; no new reaction experiment.", "evidence_ids": [E[4], T[4]]},
            {"region": "references", "finding": "Bibliography only; cited papers not imported as current-paper experiments.", "evidence_ids": [E[4], T[4]]}],
         "record_ids": page_ids[4],
         "unresolved": ["No new experiment on page 4; mechanism remains author-proposed rather than experimentally complete."]},
    ],
}
write(ROOT / "runs/current/coverage/pages_0001_0004.json", coverage)

tasks = json.loads((ROOT / "runs/current/tasks.json").read_text())
for task in tasks:
    task["state"] = "attempted"
    task["reason"] = "Original page image inspected; findings and record IDs are in coverage/pages_0001_0004.json."
write(ROOT / "runs/current/tasks.json", tasks)

print(json.dumps(inventory["counts"], ensure_ascii=False))
