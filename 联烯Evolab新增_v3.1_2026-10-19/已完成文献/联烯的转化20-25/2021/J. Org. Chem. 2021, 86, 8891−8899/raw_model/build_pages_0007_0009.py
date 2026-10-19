#!/usr/bin/env python3
"""Build evidence-linked v3.1 additions/revisions for PDF pages 7-9."""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_71ecfcbf02cd66de"
CURRENT = ROOT / "runs/current"
BASE = st.read(CURRENT / "standard_records.json")
INIT = st.read(CURRENT / "raw_model/pages_0007_0009_init.json")

# Start from the current v3 registry so updated cross-page experiments replace,
# rather than duplicate, their page-3 identities. The init output supplies the
# immutable prepared document/evidence registry required by standard merge.
bundle = copy.deepcopy(BASE)
bundle["documents"] = copy.deepcopy(INIT["documents"])
molecules = {m["molecule_id"]: m for m in bundle["molecules"]}
reactions = {r["reaction_id"]: r for r in bundle["reactions"]}


def evid(pages):
    out = []
    for page in pages:
        out.extend([f"{DOC}:page{page}:text", f"{DOC}:page{page}:image"])
    return out


def fact(raw=None, normalized=None, pages=(), method="original_extraction",
         status=None, verification=None, derivation=None, alternatives=()):
    return st.fact(raw, normalized, evid(pages), method=method, status=status,
                   verification=verification, derivation=derivation,
                   alternatives=alternatives)


def missing(pages=(), status="unreported", derivation=None):
    return fact(None, None, pages, status=status, verification="missing",
                derivation=derivation)


def add_unique(seq, values):
    for value in values:
        if value not in seq:
            seq.append(value)


def occurrence_evidence(m, pages):
    evs = [f"{DOC}:page{p}:text" for p in pages]
    occ = next(o for o in m["occurrences"] if o["document_id"] == DOC)
    add_unique(occ["evidence_ids"], evs)
    add_unique(occ["compound_label"]["evidence_ids"], evid(pages))


def set_stereo(m, pages, chirality, geometry, representation,
               atom_centers=None, double_bonds=None, axial=None,
               mixture=None, assignment=None):
    s = m["stereochemistry"]
    s["chirality_status"] = fact(chirality, chirality, pages)
    s["geometry_status"] = fact(geometry, geometry, pages)
    s["representation_status"] = fact(representation, representation, pages)
    s["atom_centers"] = fact(atom_centers, atom_centers, pages) if atom_centers is not None else missing(pages)
    s["double_bonds"] = fact(double_bonds, double_bonds, pages) if double_bonds is not None else missing(pages)
    s["axial"] = fact(axial, axial, pages) if axial is not None else missing(pages)
    s["mixture"] = fact(mixture, mixture, pages) if mixture is not None else missing(pages)
    s["assignment"] = fact(assignment, assignment, pages) if assignment is not None else missing(pages)


def molecule_unresolved(key, label, name, pages, note, chirality="unresolved",
                        geometry="unresolved", representation="not_encoded",
                        double_bonds=None, mixture=None, assignment=None):
    mid = f"{DOC}:{key}"
    m = molecules.get(mid) or st.new_molecule(mid, DOC, [f"{DOC}:page{p}:text" for p in pages])
    m["identity_status"] = "source_identified" if name else "unresolved"
    m["occurrences"][0]["compound_label"] = fact(label, label, pages)
    m["names"]["original"] = fact(name or label, name or label, pages)
    m["names"]["english"] = fact(name or label, name or label, pages)
    m["cas_numbers"] = missing()
    m["representation_type"] = "unresolved"
    m["smiles"] = missing()
    m["generic_template"] = None
    set_stereo(m, pages, chirality, geometry, representation,
               double_bonds=double_bonds, mixture=mixture, assignment=assignment)
    m["stereochemistry"]["isomeric_smiles"] = missing(pages)
    m["stereochemistry"]["identity_key"] = missing(pages)
    m["identity_note"] = fact(note, note, pages)
    m["external_properties"] = []
    molecules[mid] = m
    return mid


def molecule_smiles(key, label, name, pages, smiles, chirality, geometry,
                    representation="partial", atom_centers=None,
                    double_bonds=None, axial=None, mixture=None,
                    assignment=None, note=None):
    mid = f"{DOC}:{key}"
    m = molecules.get(mid) or st.new_molecule(mid, DOC, [f"{DOC}:page{p}:text" for p in pages])
    occurrence_evidence(m, pages)
    diag = st.stereo_inspect(smiles)
    canonical = diag["canonical_isomeric_smiles"]
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = fact(label, label, pages)
    m["names"]["original"] = fact(name, name, pages)
    m["names"]["english"] = fact(name, name, pages)
    m["cas_numbers"] = missing()
    m["representation_type"] = "smiles"
    derivation = (
        "Connectivity candidate derived from the printed English compound name and checked by "
        "standard stereo/RDKit. Parsing and formula consistency do not establish source drawing correctness."
    )
    m["smiles"] = fact(smiles, canonical, pages, method="structure_derivation", derivation=derivation)
    m["generic_template"] = None
    centers = atom_centers if atom_centers is not None else diag["atom_centers"]
    set_stereo(m, pages, chirality, geometry, representation,
               atom_centers=centers, double_bonds=double_bonds,
               axial=axial, mixture=mixture, assignment=assignment)
    m["stereochemistry"]["isomeric_smiles"] = fact(
        smiles, canonical, pages, method="structure_derivation",
        derivation="Canonical isomeric SMILES for the connectivity candidate; unassigned stereo remains unencoded."
    )
    m["stereochemistry"]["identity_key"] = fact(
        None, diag["stereo_graph_key"], pages, method="structure_derivation",
        derivation="standard stereo graph key for this machine connectivity candidate."
    )
    m["identity_note"] = fact(
        note or "Machine connectivity candidate from the printed name; stereochemical identity remains needs_review.",
        note or "Machine connectivity candidate from the printed name; stereochemical identity remains needs_review.",
        pages,
    )
    m["external_properties"] = []
    molecules[mid] = m
    return mid


# Purification agents explicitly printed on pages 7-8.
hexane = molecule_smiles("hexane", "hexane", "hexane", (7, 8), "CCCCCC",
                         "achiral", "not_applicable", "encoded",
                         atom_centers=[], assignment="No stereogenic element in the submitted candidate graph.")
isopropanol = molecule_smiles("isopropanol", "2-propanol", "2-propanol", (7, 8), "CC(O)C",
                             "achiral", "not_applicable", "encoded",
                             atom_centers=[], assignment="No stereogenic element in the submitted candidate graph.")
ethyl_acetate = molecule_smiles("ethyl_acetate", "ethyl acetate", "ethyl acetate", (8,), "CCOC(C)=O",
                               "achiral", "not_applicable", "encoded",
                               atom_centers=[], assignment="No stereogenic element in the submitted candidate graph.")
dcm = molecule_smiles("dcm", "dichloromethane", "dichloromethane", (7,), "ClCCl",
                      "achiral", "not_applicable", "encoded",
                      atom_centers=[], assignment="No stereogenic element in the submitted candidate graph.")
ether = molecule_unresolved(
    "ether_source", "ether", "ether", (7, 8),
    "The source says only 'ether'; composition is not expanded, so no unique solvent structure is asserted.",
)

# Page-7 named products. The structures are candidates from explicit names;
# stereocentres are intentionally unassigned and dr mixtures are not collapsed.
named = {}
named["3d"] = molecule_smiles(
    "3d", "3d", "5-(3,3-Dimethyl-2-((2,2,6,6-tetramethylpiperidin-1-yl)oxy)butanoyl)-1-methylpyrrolidin-2-one",
    (7,), "CN1C(=O)CCC1C(=O)C(ON2C(C)(C)CCCC2(C)C)C(C)(C)C",
    "mixture", "not_applicable", mixture={"type": "diastereomeric mixture", "ratio": "1:1.8", "configurations": "unassigned"})
named["2e"] = molecule_smiles(
    "2e", "2e", "1-Methyl-5-(2-phenyl-1-((2,2,6,6-tetramethylpiperidin-1-yl)oxy)vinyl)pyrrolidin-2-one",
    (7,), "CN1C(=O)CCC1C(ON2C(C)(C)CCCC2(C)C)=C(c3ccccc3)",
    "mixture", "mixture", mixture={"type": "geometric/diastereomeric mixture", "ratio": "10:1"},
    double_bonds={"configuration": "conflict: compound heading says (Z)-2e; spectral subheading says (E)-2e"},
    assignment="The paragraph reports a 10:1 mixture; absolute configuration at the lactam stereocentre is unassigned.",
    note="The source internally conflicts between '(Z)-2e' in the compound heading and '(E)-2e' before the NMR block; no single geometry is chosen.")
named["3f"] = molecule_smiles(
    "3f", "3f", "1-Isopropyl-5-(2-((2,2,6,6-tetramethylpiperidin-1-yl)oxy)propanoyl)pyrrolidin-2-one",
    (7,), "CC(C)N1C(=O)CCC1C(=O)C(C)ON2C(C)(C)CCCC2(C)C",
    "mixture", "not_applicable", mixture={"type": "diastereomeric mixture", "ratio": "1:1.1", "configurations": "unassigned"})
named["3g"] = molecule_smiles(
    "3g", "3g", "1-Benzyl-5-(2-((2,2,6,6-tetramethylpiperidin-1-yl)oxy)propanoyl)pyrrolidin-2-one",
    (7,), "c1ccccc1CN2C(=O)CCC2C(=O)C(C)ON3C(C)(C)CCCC3(C)C",
    "mixture", "not_applicable", mixture={"type": "diastereomeric mixture", "ratio": "1.2:1", "configurations": "unassigned"})
named["5g"] = molecule_smiles(
    "5g", "5g", "5-Acryloyl-1-benzylpyrrolidin-2-one", (7,),
    "c1ccccc1CN2C(=O)CCC2C(=O)C=C", "unspecified", "not_applicable")
named["2h"] = molecule_smiles(
    "2h", "2h", "1-Methyl-5-(1-((2,2,6,6-tetramethylpiperidin-1-yl)oxy)vinyl)pyrrolidin-2-one",
    (7,), "CN1C(=O)CCC1C(=C)ON2C(C)(C)CCCC2(C)C", "unspecified", "not_applicable")

# The 1i paragraph has an irreconcilable label conflict, so no concrete SMILES is
# assigned to 2i/3i/5i from this page alone.
named["2i"] = molecule_unresolved(
    "2i", "(Z)-2i", "(Z)-2i", (3, 7),
    "Table 3 and the paragraph identify a (Z)-2i product and NOE assignment, but the page-7 compound heading labels the named structure 3i; structure identity is unresolved.",
    chirality="unresolved", geometry="specified", representation="not_encoded",
    double_bonds={"configuration": "Z", "assignment_method": "NOE"},
    assignment="An analytical sample and NOE support Z geometry; the compound-number/name association conflicts on page 7.")
named["5i"] = molecule_unresolved(
    "5i", "5i", "5i", (7,),
    "The page states that a 1H NMR yield of 5i was obtained but does not print its structure or numeric yield here.")
occurrence_evidence(molecules[f"{DOC}:3i"], (7,))
molecules[f"{DOC}:3i"]["identity_note"] = fact(
    "Page-7 heading assigns the printed name to 3i, while the sentence distinguishes 3i from (Z)-2i and 5i; identity remains unresolved.",
    "Page-7 heading assigns the printed name to 3i, while the sentence distinguishes 3i from (Z)-2i and 5i; identity remains unresolved.", (7,))

named["2j"] = molecule_smiles(
    "2j", "(Z)-2j", "(Z)-5-(tert-Butyldimethylsilyl)-1-methyl-5-(1-((2,2,6,6-tetramethylpiperidin-1-yl)oxy)prop-1-en-1-yl)pyrrolidin-2-one",
    (7,), "CN1C(=O)CCC1([Si](C)(C)C(C)(C)C)C(ON2C(C)(C)CCCC2(C)C)=CC",
    "unspecified", "specified", double_bonds={"configuration": "Z"},
    assignment="Z geometry is explicitly assigned; configuration at the potential tetrahedral centre is unreported.")
named["2k"] = molecule_smiles(
    "2k", "2k", "1-Methyl-5-(2-methyl-1-((2,2,6,6-tetramethylpiperidin-1-yl)oxy)prop-1-en-1-yl)pyrrolidin-2-one",
    (7,), "CN1C(=O)CCC1C(=C(C)C)ON2C(C)(C)CCCC2(C)C", "unspecified", "not_applicable")
named["5k"] = molecule_smiles(
    "5k", "5k", "5-Methacryloyl-1-methylpyrrolidin-2-one", (7,),
    "CN1C(=O)CCC1C(=O)C(=C)C", "unspecified", "not_applicable")
named["2l"] = molecule_smiles(
    "2l", "2l", "5-(Cyclopentylidene((2,2,6,6-tetramethylpiperidin-1-yl)oxy)methyl)-1-methylpyrrolidin-2-one",
    (7,), "CN1C(=O)CCC1C(=C2CCCC2)ON3C(C)(C)CCCC3(C)C", "unspecified", "not_applicable")
named["5l"] = molecule_smiles(
    "5l", "5l", "5-(Cyclopent-1-ene-1-carbonyl)-1-methylpyrrolidin-2-one", (7, 8),
    "CN1C(=O)CCC1C(=O)C2=CCCC2", "unspecified", "not_applicable")

# Table-2 product names are printed on page 8. All are unresolved dr mixtures.
named["3aa"] = molecule_smiles(
    "3aa", "3aa", "5-(2-((4-Hydroxy-2,2,6,6-tetramethylpiperidin-1-yl)oxy)propanoyl)-1-methylpyrrolidin-2-one",
    (8,), "CN1C(=O)CCC1C(=O)C(C)ON2C(C)(C)CC(O)CC2(C)C",
    "mixture", "not_applicable", mixture={"type": "diastereomeric mixture", "ratio": "1:1.4", "configurations": "unassigned"})
named["3ab"] = molecule_smiles(
    "3ab", "3ab", "5-(2-((4-Methoxy-2,2,6,6-tetramethylpiperidin-1-yl)oxy)propanoyl)-1-methylpyrrolidin-2-one",
    (8,), "CN1C(=O)CCC1C(=O)C(C)ON2C(C)(C)CC(OC)CC2(C)C",
    "mixture", "not_applicable", mixture={"type": "diastereomeric mixture", "ratio": "1:1", "configurations": "unassigned"})
named["3ac"] = molecule_smiles(
    "3ac", "3ac", "2,2,6,6-Tetramethyl-1-((1-(1-methyl-5-oxopyrrolidin-2-yl)-1-oxopropan-2-yl)oxy)piperidin-4-yl benzoate",
    (8,), "CN1C(=O)CCC1C(=O)C(C)ON2C(C)(C)CC(OC(=O)c3ccccc3)CC2(C)C",
    "mixture", "not_applicable", mixture={"type": "diastereomeric mixture", "ratio": "1:1", "configurations": "unassigned"})
named["3ad"] = molecule_smiles(
    "3ad", "3ad", "2,2,6,6-Tetramethyl-1-((1-(1-methyl-5-oxopyrrolidin-2-yl)-1-oxopropan-2-yl)oxy)piperidin-4-one",
    (8,), "CN1C(=O)CCC1C(=O)C(C)ON2C(C)(C)CC(=O)CC2(C)C",
    "mixture", "not_applicable", mixture={"type": "diastereomeric mixture", "ratio": "1:1.1", "configurations": "unassigned"})

# Page-8 new substrates and products. A second provisional 1m identity is used
# because the page-8 product scaffold conflicts with the page-6 compound named 1m.
new_mids = {}
new_mids["1m_page8"] = molecule_unresolved(
    "1m_page8", "1m", "substrate 1m (page-8 occurrence)", (8,),
    "The page-8 4m product scaffold is incompatible with the page-6 named penta-2,3-dien-1-yl carbamate also labelled 1m; the occurrences are not merged pending expert review.")
new_mids["1n"] = molecule_unresolved(
    "1n", "1n", "substrate 1n", (8,), "Only the source label is printed on this page; structure is unreported here.")
new_mids["1o"] = molecule_unresolved(
    "1o", "1o", "substrate 1o", (8,), "Only the source label is printed on this page; structure is unreported here.")
new_mids["1prime"] = molecule_unresolved(
    "1prime", "1′", "substrate 1′", (8,), "Only the source label is printed on this page; structure is unreported here.")
new_mids["4m"] = molecule_smiles(
    "4m", "4m", "N-(((2,2,6,6-Tetramethylpiperidin-1-yl)oxy)methyl)octa-5,6-dienamide",
    (8,), "CC=C=CCCCC(=O)NCON1C(C)(C)CCCC1(C)C", "unspecified", "not_applicable",
    axial={"element": "allene axis", "configuration": "unassigned"},
    assignment="The allene axis is potentially stereogenic; no axial configuration or mixture composition is reported.")
new_mids["4n"] = molecule_smiles(
    "4n", "4n", "Penta-2,3-dien-1-yl (((2,2,6,6-tetramethylpiperidin-1-yl)oxy)methyl)carbamate",
    (8,), "CC=C=CCOC(=O)NCON1C(C)(C)CCCC1(C)C", "unspecified", "not_applicable",
    axial={"element": "allene axis", "configuration": "unassigned"},
    assignment="The allene axis is potentially stereogenic; no axial configuration or mixture composition is reported.")
new_mids["3o"] = molecule_smiles(
    "3o", "3o", "tert-Butyl 2-(2-((2,2,6,6-tetramethylpiperidin-1-yl)oxy)propanoyl)pyrrolidine-1-carboxylate",
    (8,), "CC(C)(C)OC(=O)N1CCCC1C(=O)C(C)ON2C(C)(C)CCCC2(C)C",
    "mixture", "not_applicable", mixture={"type": "diastereomeric mixture", "ratio": "2.14:1", "configurations": "unassigned"})
new_mids["4prime"] = molecule_smiles(
    "4prime", "4′", "N-(((2,2,6,6-Tetramethylpiperidin-1-yl)oxy)methyl)acetamide",
    (8,), "CC(=O)NCON1C(C)(C)CCCC1(C)C", "achiral", "not_applicable", "encoded",
    atom_centers=[], assignment="No stereogenic element in the submitted machine graph; source reports rotamers in NMR.")


def participant(pid, mid, role, pages, stage="reaction", mass=None, amount=None,
                equivalents=None, product_role=None, reference=None):
    p = st.new_participant(pid, mid, role, evid(pages), stage)
    p["product_role"] = fact(product_role, product_role, pages) if product_role else missing()
    p["charges"]["mass"] = fact(str(mass), mass, pages) if mass else missing(pages)
    p["charges"]["amount"] = fact(str(amount), amount, pages) if amount else missing(pages)
    p["charges"]["volume"] = missing(pages)
    p["charges"]["equivalents"] = fact(str(equivalents), equivalents, pages) if equivalents else missing(pages)
    p["charges"]["mol_percent"] = missing(pages)
    p["equivalence_reference"] = fact(reference, reference, pages) if reference else missing()
    return p


def outcome(oid, kind, target, pages, raw, value, basis=None, yield_type=None, scope=None, role="product"):
    o = st.new_outcome(oid, kind, [target], evid(pages))
    o["product_role"] = fact(role, role, pages)
    o["value"] = fact(raw, value, pages)
    o["yield_type"] = fact(yield_type, yield_type, pages) if yield_type else missing()
    o["measurement_basis"] = fact(basis, basis, pages) if basis else missing(pages)
    o["scope"] = fact(scope, scope, pages) if scope else missing()
    return o


def formation_labels(r, pages, codes):
    labels = []
    changes = []
    elements = {"C-O": ("C", "O"), "C-N": ("C", "N"), "DB-C-O": ("C", "O"), "N-H": ("N", "H")}
    orders = {"C-O": 1, "C-N": 1, "DB-C-O": 2, "N-H": 1}
    numbers = {"C-O": 3, "C-N": 5, "DB-C-O": 14, "N-H": 29}
    for idx, code in enumerate(codes, 1):
        a, b = elements[code]
        cid = f"b{idx}"
        change = {"change_id": cid, "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
                  "before_order": 0, "after_order": orders[code],
                  "change_type": "hydrogen_gain" if code == "N-H" else "formed"}
        changes.append(fact(None, change, pages, method="structure_derivation",
                            derivation="Partial net bond interpretation from printed product name/structure; atom mapping is unresolved."))
        labels.append({"category_number": numbers[code], "code": code, "supporting_change_ids": [cid]})
    r["classification"]["bond_changes"] = changes
    r["classification"]["labels"] = fact(None, labels, pages, method="structure_derivation",
                                           derivation="Multi-label classification from only the visually/name-supported new bonds.")
    r["classification"]["coverage"] = fact(None, "partial", pages, method="inference",
                                             derivation="No complete atom mapping; leaving-group and hydrogen-source accounting remain unresolved.")
    r["classification"]["limitations"] = [
        "No atom-mapped reaction SMILES was asserted.",
        "Only source-supported net formed-bond classes are represented; all bond breaking and hydrogen provenance are incomplete.",
    ]
    r["classification"]["model_votes"] = []


def set_amounts(r, page, substrate_mass, substrate_amount):
    sub = next(p for p in r["participants"] if p["participant_id"] == "substrate")
    sub["charges"]["mass"] = fact(str(substrate_mass), substrate_mass, (page,))
    sub["charges"]["amount"] = fact(str(substrate_amount), substrate_amount, (page,))
    r["reaction_scale"] = fact(str(substrate_amount),
        {"reference_molecule_id": sub["molecule_id"], "amount": substrate_amount}, (page,))


def add_purification(r, page, agents, operation):
    for pid, mid in agents:
        if not any(p["participant_id"] == pid for p in r["participants"]):
            r["participants"].append(participant(pid, mid, "purification_agent", (page,), stage="purification"))
    sid = "column_chromatography"
    r["steps"] = [s for s in r["steps"] if s["step_id"] != sid]
    s = st.new_step(sid, evid((page,)), "purification")
    s["participant_ids"] = [pid for pid, _ in agents]
    s["operation"] = fact(operation, operation, (page,))
    s["temperature"] = missing((page,))
    s["time"] = missing((page,))
    s["atmosphere"] = missing()
    s["pressure"] = missing()
    r["steps"].append(s)


def add_amount_outcomes(r, page, product_key, mass, amount):
    mid = f"{DOC}:{product_key}"
    for oid in (f"product_mass_{product_key}", f"product_amount_{product_key}"):
        r["outcomes"] = [o for o in r["outcomes"] if o["outcome_id"] != oid]
    r["outcomes"].append(outcome(f"product_mass_{product_key}", "product_amount", mid, (page,), str(mass), mass,
                                  "isolated after column chromatography", scope={"measure": "mass"}))
    r["outcomes"].append(outcome(f"product_amount_{product_key}", "product_amount", mid, (page,), str(amount), amount,
                                  "isolated after column chromatography", scope={"measure": "amount"}))


def replace_outcome(r, new):
    r["outcomes"] = [o for o in r["outcomes"] if o["outcome_id"] != new["outcome_id"]]
    r["outcomes"].append(new)


def touch_reaction(r, pages):
    add_unique(r["evidence_ids"], [f"{DOC}:page{p}:text" for p in pages])
    old_pages = r["location"].get("normalized_value", {}).get("pages", [])
    all_pages = sorted(set(old_pages) | set(pages))
    r["location"] = fact({"pages": all_pages}, {"pages": all_pages}, tuple(all_pages))
    r["review"]["status"] = "needs_review"
    r["review"]["method"] = "model"
    r["completeness"] = "partial"


# Revise existing Table-3 records with the concrete page-7/8 preparation data.
detail = {
    "scope_1d": (7, {"number": 41.6, "unit": "mg"}, {"number": 0.114, "unit": "mmol"}, [("3d", {"number": 24.8, "unit": "mg"}, {"number": 0.074, "unit": "mmol"})], [("ether", ether), ("hexane", hexane), ("isopropanol", isopropanol)]),
    "scope_1e": (7, {"number": 26.1, "unit": "mg"}, {"number": 0.068, "unit": "mmol"}, [("2e", {"number": 8.7, "unit": "mg"}, {"number": 0.024, "unit": "mmol"})], [("ether", ether), ("hexane", hexane)]),
    "scope_1f": (7, {"number": 17.3, "unit": "mg"}, {"number": 0.0495, "unit": "mmol"}, [("3f", {"number": 9.1, "unit": "mg"}, {"number": 0.075, "unit": "mmol"})], [("dcm", dcm), ("ether", ether)]),
    "scope_1g": (7, {"number": 43.4, "unit": "mg"}, {"number": 0.109, "unit": "mmol"}, [("3g", {"number": 8.3, "unit": "mg"}, {"number": 0.021, "unit": "mmol"}), ("5g", {"number": 0.8, "unit": "mg"}, {"number": 0.003, "unit": "mmol"})], [("ether", ether), ("hexane", hexane), ("dcm", dcm)]),
    "scope_1h": (7, {"number": 0.020, "unit": "g"}, {"number": 0.065, "unit": "mmol"}, [("2h", {"number": 8.3, "unit": "mg"}, {"number": 0.030, "unit": "mmol"})], [("dcm", dcm), ("ether", ether)]),
    "scope_1j": (7, {"number": 37.5, "unit": "mg"}, {"number": 0.086, "unit": "mmol"}, [("2j", {"number": 14.2, "unit": "mg"}, {"number": 0.035, "unit": "mmol"})], [("dcm", dcm), ("ether", ether)]),
    "scope_1k": (7, {"number": 36.5, "unit": "mg"}, {"number": 0.109, "unit": "mmol"}, [("2k", {"number": 9.4, "unit": "mg"}, {"number": 0.030, "unit": "mmol"}), ("5k", {"number": 7.9, "unit": "mg"}, {"number": 0.047, "unit": "mmol"})], [("ether", ether), ("hexane", hexane), ("isopropanol", isopropanol)]),
    "scope_1l": (7, {"number": 36.5, "unit": "mg"}, {"number": 0.101, "unit": "mmol"}, [("2l", {"number": 5.3, "unit": "mg"}, {"number": 0.016, "unit": "mmol"}), ("5l", {"number": 9.1, "unit": "mg"}, {"number": 0.047, "unit": "mmol"})], [("ether", ether), ("hexane", hexane), ("isopropanol", isopropanol)]),
}
for suffix, (page, smass, samt, products, agents) in detail.items():
    r = reactions[f"{DOC}:{suffix}"]
    touch_reaction(r, (page,) if suffix != "scope_1l" else (7, 8))
    set_amounts(r, page, smass, samt)
    for key, pmass, pamt in products:
        mid = f"{DOC}:{key}"
        if not any(p["molecule_id"] == mid for p in r["participants"]):
            r["participants"].append(participant(f"product_{key}", mid, "product", (page,), product_role="reported product/byproduct"))
        add_amount_outcomes(r, page, key, pmass, pamt)
    add_purification(r, page, agents, "purification by column chromatography using the printed solvent sequence")

# Product-specific page-7 results and conflicts.
r = reactions[f"{DOC}:scope_1g"]
replace_outcome(r, outcome("yield_5g", "yield", named["5g"], (7,), "3%", {"number": 3, "unit": "%"},
                           "isolated after column chromatography", "isolated", {"product": "5g"}, "byproduct"))

r = reactions[f"{DOC}:scope_1f"]
r["notes"].append(fact(
    "9.1 mg (54%, 0.075 mmol)",
    "The printed product amount 0.075 mmol exceeds the printed 0.0495 mmol substrate charge and is inconsistent with the stated 54% yield; values retained without correction.",
    (7,)))
add_unique(r["review"]["unresolved"], ["Printed 3f product mass/yield/amount are internally inconsistent; no author correction was made."])

r = reactions[f"{DOC}:scope_1e"]
r["notes"].append(fact(
    {"compound_heading": "(Z)-2e", "spectral_heading": "(E)-2e"}, None, (7,),
    status="conflict", alternatives=[
        fact("(Z)-2e", "Z", (7,)), fact("(E)-2e", "E", (7,))
    ]))
add_unique(r["review"]["unresolved"], ["Page-7 2e geometry is internally conflicting: compound heading Z, spectral subheading E; the 10:1 mixture remains unassigned beyond source wording."])

r = reactions[f"{DOC}:scope_1k"]
old30 = fact("30% yield of the TMP enol ether 2k", {"number": 30, "unit": "%"}, (3,))
new28 = fact("9.4 mg (0.030 mmol, 28%) of 2k", {"number": 28, "unit": "%"}, (7,))
conflict = outcome("yield_2k", "yield", named["2k"], (3, 7),
                   {"discussion": "30%", "specific preparation": "28%"}, None,
                   "isolated after column chromatography", "isolated", {"product": "2k"})
conflict["value"] = fact({"discussion": "30%", "specific preparation": "28%"}, None, (3, 7),
                         status="conflict", alternatives=[old30, new28])
replace_outcome(r, conflict)
r["notes"].append(fact("30% in page-3 discussion; 28% in page-7 specific preparation", None, (3, 7),
                       status="conflict", alternatives=[old30, new28]))
add_unique(r["review"]["unresolved"], ["2k isolated yield conflicts between page-3 discussion (30%) and page-7 preparation (28%)."])

# Correct and expand the existing 1i experiment: Table 3 assigns 15% to 2i, not 3i.
r = reactions[f"{DOC}:scope_1i"]
touch_reaction(r, (7,))
set_amounts(r, 7, {"number": 41.8, "unit": "mg"}, {"number": 0.115, "unit": "mmol"})
r["participants"] = [p for p in r["participants"] if p["participant_id"] != "product_3i"]
for key in ("2i", "3i", "5i"):
    mid = f"{DOC}:{key}"
    if not any(p["molecule_id"] == mid for p in r["participants"]):
        r["participants"].append(participant(f"product_{key}", mid, "product", (7,), product_role="reported product"))
r["outcomes"] = [o for o in r["outcomes"] if o["outcome_id"] != "yield_3i"]
r["outcomes"].append(outcome("yield_2i", "yield", named["2i"], (3, 7), "15%", {"number": 15, "unit": "%"},
                             "crude 1H NMR with 1,3,5-trimethoxybenzene as internal standard", "1H NMR", {"table": "Table 3"}))
for key in ("3i", "5i"):
    r["outcomes"].append(outcome(f"nmr_yield_{key}_unprinted", "other", f"{DOC}:{key}", (7,),
                                 f"1H NMR yield of {key} stated but numeric value not printed on page 7",
                                 f"1H NMR yield reported without a numeric value on this page",
                                 "1H NMR with 1,3,5-trimethoxybenzene as internal standard",
                                 scope={"value_state": "unreported on page 7"}))
r["outcomes"].append(outcome("z_2i_noe", "other", named["2i"], (7,),
                             "An analytical sample of (Z)-2i was obtained and an NOE was conducted to support formation of the Z-isomer.",
                             {"geometry": "Z", "assignment_method": "NOE"}, "NOE", scope={"metric": "geometry assignment"}))
r["notes"].append(fact(
    {"heading": "(Z)-... (3i)", "sentence": "3i, (Z)-2i, and 5i"}, None, (7,), status="conflict",
    alternatives=[fact("compound heading labels named structure 3i", "3i", (7,)),
                  fact("sentence identifies a separate (Z)-2i analytical sample", "2i", (7,))]))
add_unique(r["review"]["unresolved"], [
    "The page-7 compound heading and paragraph conflict over the 2i/3i identity; Table 3 assigns the 15% value to 2i.",
    "Numeric 1H NMR yields of 3i and 5i are not printed on page 7.",
])

# Table-2 rows gain their concrete charges and page-8 purification details.
table2 = {
    "table2_entry_1": ({"number": 27.0, "unit": "mg"}, {"number": 0.084, "unit": "mmol"}, "3aa", {"number": 14.4, "unit": "mg"}, {"number": 0.0441, "unit": "mmol"}),
    "table2_entry_2": ({"number": 35.7, "unit": "mg"}, {"number": 0.111, "unit": "mmol"}, "3ab", {"number": 18.8, "unit": "mg"}, {"number": 0.0552, "unit": "mmol"}),
    "table2_entry_3": ({"number": 33.0, "unit": "mg"}, {"number": 0.102, "unit": "mmol"}, "3ac", {"number": 27.7, "unit": "mg"}, {"number": 0.065, "unit": "mmol"}),
    "table2_entry_4": ({"number": 30.1, "unit": "mg"}, {"number": 0.0937, "unit": "mmol"}, "3ad", {"number": 13.2, "unit": "mg"}, {"number": 0.0407, "unit": "mmol"}),
}
for suffix, (smass, samt, product, pmass, pamt) in table2.items():
    r = reactions[f"{DOC}:{suffix}"]
    touch_reaction(r, (8,))
    set_amounts(r, 8, smass, samt)
    add_amount_outcomes(r, 8, product, pmass, pamt)
    add_purification(r, 8, [("ether", ether), ("hexane", hexane), ("isopropanol", isopropanol)],
                     "purification by column chromatography using ether in hexane and then 2-propanol in ether")

r = reactions[f"{DOC}:table2_entry_2"]
r["notes"].append(fact("18.8 mg ... of 3aa", "The 3ab preparation paragraph calls the isolated material 3aa once; heading/Table 2 identify 3ab.", (8,)))
add_unique(r["review"]["unresolved"], ["Page-8 3ab paragraph contains a product-label typo/conflict ('3aa'); no silent source correction was made."])
for suffix, printed in (("table2_entry_3", "1ac"), ("table2_entry_4", "1ad")):
    r = reactions[f"{DOC}:{suffix}"]
    r["notes"].append(fact(
        {"Table 2 substrate": "1a", "page-8 preparation substrate": printed}, None, (3, 8), status="conflict",
        alternatives=[fact("Table 2: 1a", "1a", (3,)), fact(f"page-8 preparation: {printed}", printed, (8,))]))
    add_unique(r["review"]["unresolved"], [f"Substrate label conflicts between Table 2 (1a) and page-8 preparation ({printed})."])


def new_photo_reaction(suffix, substrate_mid, product_mid, title, page, smass, samt,
                       yvalue, ybasis, pmass=None, pamount=None, dr=None,
                       categories=("substrate_scope",), codes=("C-O", "N-H"),
                       purification=(), original_text=None, extra_unresolved=()):
    rid = f"{DOC}:{suffix}"
    r = st.new_reaction(rid, DOC, evid((page,)))
    r["record_type"] = "specific_experiment"
    r["experiment_categories"] = fact(list(categories), list(categories), (page,))
    r["location"] = fact({"pages": [page]}, {"pages": [page]}, (page,))
    r["general_scheme_ids"] = [f"{DOC}:general_visible_light_amidation"]
    r["reaction_name"] = fact(title, title, (page,))
    r["transformation_type"] = fact("visible-light-mediated reaction under the general procedure",
                                     "visible-light-mediated reaction under the general procedure", (page,))
    r["participants"] = [
        participant("substrate", substrate_mid, "reactant", (page,), mass=smass, amount=samt,
                    equivalents={"number": 1, "unit": "equiv"}, reference=substrate_mid),
        participant("tempo", f"{DOC}:tempo", "reagent", (6, page), equivalents={"number": 3, "unit": "equiv"}, reference=substrate_mid),
        participant("base", f"{DOC}:k2co3", "reagent", (6, page), equivalents={"number": 2, "unit": "equiv"}, reference=substrate_mid),
        participant("solvent", f"{DOC}:mecn", "solvent", (6, page)),
        participant("product", product_mid, "product", (page,), product_role="reported product"),
    ]
    r["reaction_scale"] = fact(str(samt), {"reference_molecule_id": substrate_mid, "amount": samt}, (page,))
    r["solvent_conditions"]["composition"] = fact("dry MeCN", f"{DOC}:mecn", (6, page), method="inherited_general_condition")
    r["solvent_conditions"]["concentration"] = fact("0.025 M", {"number": 0.025, "unit": "M"}, (6, page), method="inherited_general_condition")
    r["solvent_conditions"]["concentration_basis"] = fact("substrate", substrate_mid, (6, page), method="inherited_general_condition")
    r["special_conditions"]["light_source"] = fact("Kessil blue LED lamp", "Kessil blue LED lamp", (6, page), method="inherited_general_condition")
    r["special_conditions"]["wavelength"] = fact("440 nm", {"number": 440, "unit": "nm"}, (6, page), method="inherited_general_condition")
    r["special_conditions"]["reactor"] = fact("flame-dried 25 mL Schlenk tube; lamp approximately 3 cm away",
        {"vessel": "25 mL Schlenk tube", "lamp_distance": "approximately 3 cm"}, (6, page), method="inherited_general_condition")
    r["procedure"]["original_text"] = fact(original_text, original_text, (page,))
    r["procedure"]["general_procedure_reference"] = fact("Following the general procedure", "page-6 general visible-light procedure", (6, page), method="inherited_general_condition")
    r["procedure"]["addition_order"] = fact(
        "substrate and K2CO3; evacuate/purge with N2 three times; TEMPO; dry MeCN; irradiate",
        ["substrate", "K2CO3", "N2 purge x3", "TEMPO", "dry MeCN", "irradiation"],
        (6, page), method="inherited_general_condition")
    s = st.new_step("irradiation", evid((6, page)), "reaction")
    s["participant_ids"] = ["substrate", "tempo", "base", "solvent"]
    s["operation"] = fact("irradiated until TLC showed complete consumption", "irradiate until TLC completion", (6, page), method="inherited_general_condition")
    s["temperature"] = missing((6, page))
    s["time"] = fact("usually 4-5 h", {"minimum": 4, "maximum": 5, "unit": "h", "raw_text": "usually 4-5 h"}, (6, page), method="inherited_general_condition")
    s["atmosphere"] = fact("under nitrogen; evacuated and purged with N2 three times", "N2", (6, page), method="inherited_general_condition")
    s["pressure"] = missing()
    r["steps"] = [s]
    replace_outcome(r, outcome("yield", "yield", product_mid, (page,), f"{yvalue['number']}%", yvalue,
                               ybasis, "isolated" if "isolated" in ybasis else "1H NMR", {"product": product_mid}))
    if pmass:
        r["outcomes"].append(outcome("product_mass", "product_amount", product_mid, (page,), str(pmass), pmass,
                                     "isolated after column chromatography", scope={"measure": "mass"}))
    if pamount:
        r["outcomes"].append(outcome("product_amount", "product_amount", product_mid, (page,), str(pamount), pamount,
                                     "isolated after column chromatography", scope={"measure": "amount"}))
    if dr:
        r["outcomes"].append(outcome("dr", "selectivity", product_mid, (page,), dr, dr,
                                     "1H NMR (ratio reported with the NMR yield)", scope={"metric": "dr", "order": "unassigned diastereomers"}))
    if purification:
        add_purification(r, page, list(purification), "purification by column chromatography using the printed solvent sequence")
    formation_labels(r, (page,), codes)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
                   "unresolved": [
                       "Substrate structure is not printed on page 8 and remains unresolved.",
                       "Reaction temperature and exact irradiation time are not reported for this individual experiment.",
                       "Complete atom mapping, bond breaking and hydrogen provenance remain unresolved.",
                       *extra_unresolved,
                   ]}
    r["completeness"] = "partial"
    reactions[rid] = r


new_photo_reaction(
    "scope_1m_alpha_oxidation", new_mids["1m_page8"], new_mids["4m"],
    "Preparation of N-(((2,2,6,6-tetramethylpiperidin-1-yl)oxy)methyl)octa-5,6-dienamide (4m)",
    8, {"number": 37.2, "unit": "mg"}, {"number": 0.111, "unit": "mmol"},
    {"number": 4, "unit": "%"}, "isolated after column chromatography",
    {"number": 1.5, "unit": "mg"}, {"number": 0.0049, "unit": "mmol"},
    purification=(("ether", ether), ("hexane", hexane)),
    original_text="Following the general procedure, 37.2 mg (0.111 mmol) of 1m was irradiated to give 1.5 mg (0.0049 mmol, 4%) of 4m.",
    extra_unresolved=("The page-8 1m occurrence conflicts structurally with the page-6 named 1m occurrence; a provisional separate molecule ID is retained.",),
)
new_photo_reaction(
    "scope_1n_alpha_oxidation", new_mids["1n"], new_mids["4n"],
    "Preparation of penta-2,3-dien-1-yl (((2,2,6,6-tetramethylpiperidin-1-yl)oxy)methyl)carbamate (4n)",
    8, {"number": 27.5, "unit": "mg"}, {"number": 0.085, "unit": "mmol"},
    {"number": 26, "unit": "%"}, "isolated after column chromatography",
    {"number": 6.5, "unit": "mg"}, {"number": 0.022, "unit": "mmol"},
    purification=(("ether", ether), ("hexane", hexane), ("isopropanol", isopropanol)),
    original_text="Following the general procedure, 27.5 mg (0.085 mmol) of 1n was irradiated to give 6.5 mg (0.022 mmol, 26%) of 4n.",
)
new_photo_reaction(
    "scope_1o", new_mids["1o"], new_mids["3o"],
    "Preparation of tert-butyl 2-(2-((2,2,6,6-tetramethylpiperidin-1-yl)oxy)propanoyl)pyrrolidine-1-carboxylate (3o)",
    8, {"number": 36.3, "unit": "mg"}, {"number": 0.0891, "unit": "mmol"},
    {"number": 22, "unit": "%"}, "1H NMR yield", dr="2.14:1",
    codes=("C-N", "C-O", "DB-C-O"),
    original_text="Following the general procedure, 36.3 mg (0.0891 mmol) of 1o was irradiated under the reaction conditions to result in a 22% 1H NMR yield with a 2.14:1 ratio of diastereomers.",
)
new_photo_reaction(
    "control_1prime_alpha_oxidation", new_mids["1prime"], new_mids["4prime"],
    "Preparation of N-(((2,2,6,6-tetramethylpiperidin-1-yl)oxy)methyl)acetamide (4′)",
    8, {"number": 20.0, "unit": "mg"}, {"number": 0.078, "unit": "mmol"},
    {"number": 65, "unit": "%"}, "isolated after column chromatography",
    {"number": 4.9, "unit": "mg"}, {"number": 0.021, "unit": "mmol"},
    categories=("control", "mechanistic_study"), codes=("C-O", "N-H"),
    purification=(("ethyl_acetate", ethyl_acetate), ("hexane", hexane)),
    original_text="Following the general procedure, 20.0 mg (0.078 mmol) of substrate 1′ was irradiated to result in 4.9 mg (0.021 mmol, 65%) of 4′.",
    extra_unresolved=("The source describes the isolated physical form as 'a light-yellow oil solid'; wording is retained without normalization.",),
)
reactions[f"{DOC}:control_1prime_alpha_oxidation"]["notes"].append(
    fact("as a light-yellow oil solid", "source physical-form wording retained verbatim", (8,)))

# Source context: characterization is associated with synthesis records but is
# not promoted to independent reaction records. Page 9 contains no new experiment.
contexts = [
    {"context_id": f"{DOC}:page7_characterization", "document_id": DOC,
     "evidence_ids": [f"{DOC}:page7:text", f"{DOC}:page7:image"], "kind": "characterization",
     "note": "1H/13C NMR and HRMS blocks for 3d, 2e, 3f, 3g/5g, 2h, 2j, 2k/5k and 2l; not independent syntheses."},
    {"context_id": f"{DOC}:page8_characterization", "document_id": DOC,
     "evidence_ids": [f"{DOC}:page8:text", f"{DOC}:page8:image"], "kind": "characterization",
     "note": "Continuation of 5l and characterization blocks for 4m, 4n, 3aa-3ad and 4′; not independent mass-spectrum experiments."},
    {"context_id": f"{DOC}:page9_back_matter", "document_id": DOC,
     "evidence_ids": [f"{DOC}:page9:text", f"{DOC}:page9:image"], "kind": "back_matter_and_references",
     "note": "Author contributions, notes, acknowledgments and bibliography only; cited reactions are prior literature, not experiments of this paper."},
]
for context in contexts:
    if context not in bundle["source_context"]:
        bundle["source_context"].append(context)

bundle["molecules"] = list(molecules.values())
bundle["reactions"] = list(reactions.values())

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 7,
            "state": "attempted",
            "regions": [
                {"region": "two-column experimental text for 3d, 2e, 3f, 3g/5g and 2h",
                 "finding": "Concrete substrate/product charges, isolated yields, dr or Z:E wording, chromatography solvents, and characterization were checked; characterization peaks were not made into reactions.",
                 "evidence_ids": [f"{DOC}:page7:text", f"{DOC}:page7:image"]},
                {"region": "right-column experimental text for 1i products, 2j, 2k/5k and beginning of 2l/5l",
                 "finding": "Existing Table-3 experiment identities were extended. The 1i paragraph has a 2i/3i label conflict; 2k has 30% versus 28% source conflict.",
                 "evidence_ids": [f"{DOC}:page7:text", f"{DOC}:page7:image"]},
            ],
            "record_ids": [f"{DOC}:{x}" for x in ("scope_1d", "scope_1e", "scope_1f", "scope_1g", "scope_1h", "scope_1i", "scope_1j", "scope_1k", "scope_1l")],
            "unresolved": [
                "1i paragraph heading and prose conflict over 2i/3i; numeric NMR yields of 3i and 5i are absent from this page.",
                "2e heading says Z while its spectral subheading says E.",
                "3f printed product amount is inconsistent with substrate amount and stated yield.",
                "Absolute configurations of dr mixtures are not assigned.",
            ],
        },
        {
            "page": 8,
            "state": "attempted",
            "regions": [
                {"region": "left column continuation of 5l and preparations 4m, 4n, 3o, 3aa, 3ab",
                 "finding": "Four independent page-8 experiments (4m, 4n, 3o and later 4′) were added; 5l and Table-2 records were merged with existing identities rather than duplicated.",
                 "evidence_ids": [f"{DOC}:page8:text", f"{DOC}:page8:image"]},
                {"region": "right column preparations/characterization 3ac, 3ad, 4′ and Associated Content",
                 "finding": "Product amounts/yields/dr and chromatography were captured. SI availability is stated, but no SI file was provided or associated, so no SI records were inferred.",
                 "evidence_ids": [f"{DOC}:page8:text", f"{DOC}:page8:image"]},
            ],
            "record_ids": [
                f"{DOC}:scope_1l", f"{DOC}:scope_1m_alpha_oxidation", f"{DOC}:scope_1n_alpha_oxidation", f"{DOC}:scope_1o",
                f"{DOC}:table2_entry_1", f"{DOC}:table2_entry_2", f"{DOC}:table2_entry_3", f"{DOC}:table2_entry_4",
                f"{DOC}:control_1prime_alpha_oxidation",
            ],
            "unresolved": [
                "Page-8 1m occurrence conflicts structurally with the page-6 compound also labelled 1m; provisional identities remain separate.",
                "Page-8 Table-2 prose contains product/substrate label conflicts (3aa in the 3ab paragraph; 1ac/1ad versus Table-2 1a).",
                "SI is mentioned by the article but was not supplied/linked in this workspace.",
            ],
        },
        {
            "page": 9,
            "state": "attempted",
            "regions": [
                {"region": "author contributions, notes, acknowledgments and references",
                 "finding": "No synthesis experiment, optimization row, reaction arrow, table, or mechanistic control of the present paper is shown; bibliography entries are prior literature only.",
                 "evidence_ids": [f"{DOC}:page9:text", f"{DOC}:page9:image"]},
            ],
            "record_ids": [],
            "unresolved": [],
        },
    ],
}

st.write(CURRENT / "raw_model/pages_0007_0009_extraction.json", bundle)
st.write(CURRENT / "submitted.json", bundle)
st.write(CURRENT / "coverage/pages_0007_0009.json", coverage)
print(json.dumps({
    "molecules": len(bundle["molecules"]),
    "reactions": len(bundle["reactions"]),
    "new_reaction_ids": [
        f"{DOC}:scope_1m_alpha_oxidation", f"{DOC}:scope_1n_alpha_oxidation",
        f"{DOC}:scope_1o", f"{DOC}:control_1prime_alpha_oxidation",
    ],
}, ensure_ascii=False))
