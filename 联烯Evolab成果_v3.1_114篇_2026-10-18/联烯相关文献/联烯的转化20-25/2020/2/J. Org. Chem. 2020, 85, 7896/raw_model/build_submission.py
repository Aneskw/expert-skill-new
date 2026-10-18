#!/usr/bin/env python3
"""Build the source-derived v3.1 bundle for PDF pages 1-6.

This is retained as a raw machine construction artifact. It does not claim
source review or chemical verification.
"""
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_b298156f3445f448"
RUN = ROOT / "runs/current"
E = {p: f"{DOC}:page{p}:text" for p in range(1, 10)}
I = {p: f"{DOC}:page{p}:image" for p in range(1, 7)}


def F(raw=None, normalized=None, evidence=(), method="original_extraction", derivation=None,
      status=None):
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation, status=status)


def setf(container, key, raw, normalized, evidence, method="original_extraction", derivation=None,
         status=None):
    container[key] = F(raw, normalized, evidence, method, derivation, status)


def molecule(mid, label, evidence, name=None, smiles=None, product_stereo=False,
             assignment=None, geometry="E"):
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, evidence)
    m["names"]["original"] = F(name or label, name or label, evidence)
    m["names"]["english"] = F(name or label, name or label, evidence)
    m["identity_note"] = F(
        f"Source label {label}; identity remains a machine extraction candidate.",
        f"Source label {label}; identity remains a machine extraction candidate.", evidence)
    if smiles:
        deriv = "Candidate connectivity transcribed from the source drawing/name; RDKit parsing is not source verification."
        m["representation_type"] = "smiles"
        m["smiles"] = F(label, smiles, evidence, "structure_derivation", deriv)
        diag = st.stereo_inspect(smiles)
        if product_stereo:
            m["stereochemistry"]["atom_centers"] = F(
                "solid/dashed stereobond at the newly formed quaternary center",
                [{"center": "new quaternary carbon", "configuration": assignment or "as drawn"}],
                evidence)
            m["stereochemistry"]["double_bonds"] = F(
                f"{geometry} allylic enol-ether geometry", [{"bond": "allylic C=C", "configuration": geometry}], evidence)
            m["stereochemistry"]["assignment"] = F(
                assignment or "configuration shown by stereobond in the source drawing",
                assignment or "configuration shown by stereobond in the source drawing", evidence)
            m["stereochemistry"]["chirality_status"] = F(
                "source specifies one configuration", "specified", evidence)
            m["stereochemistry"]["geometry_status"] = F(
                f"source specifies {geometry}", "specified", evidence)
            m["stereochemistry"]["representation_status"] = F(
                "E geometry encoded; tetrahedral configuration not encoded in candidate SMILES",
                "partial", evidence, "inference",
                "The candidate SMILES preserves the source E alkene but deliberately leaves the quaternary center untagged pending expert structure mapping.")
        else:
            for key in ("chirality_status", "geometry_status", "representation_status"):
                m["stereochemistry"][key] = F(
                    "not established in bounded extraction", "unresolved", evidence, "inference",
                    "No stereochemical status is inferred from absence of source tags.")
        m["stereochemistry"]["isomeric_smiles"] = F(
            label, diag["canonical_isomeric_smiles"], evidence, "structure_derivation", deriv)
        m["stereochemistry"]["identity_key"] = F(
            None, diag["stereo_graph_key"], evidence, "structure_derivation",
            "Computed from the candidate isomeric molecular graph; not a sample identity.")
    return m


def generic_molecule(mid, label, evidence, variables):
    m = molecule(mid, label, evidence)
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": F(
            f"Source generic drawing for {label}",
            {"description": f"Generic source scaffold {label}; fixed bonds are preserved only as a drawing description."}, evidence),
        "variables": [
            {
                "label": v,
                "attachment_points": F(f"attachment marked {v}", f"source attachment marked {v}", evidence),
                "definition": F(f"source variable {v}", f"source variable {v}", evidence),
                "allowed_values": F(None, None, evidence),
                "constraints": F(None, None, evidence),
            }
            for v in variables
        ],
        "instance_mapping": [],
        "reason": F(
            "R-group/generic drawing", "R-group/generic drawing is not a concrete molecule", evidence),
    }
    return m


M = {}


def add_mol(*args, **kwargs):
    m = molecule(*args, **kwargs)
    M[m["molecule_id"]] = m
    return m["molecule_id"]


def add_generic(*args, **kwargs):
    m = generic_molecule(*args, **kwargs)
    M[m["molecule_id"]] = m
    return m["molecule_id"]


# Generic drawings on page 1 and in Schemes 2/3 are templates, not concrete SMILES.
add_generic(f"{DOC}:mol:generic_naphthol", "generic naphthol", [I[1]], ["R1", "R2"])
add_generic(f"{DOC}:mol:generic_alkoxyallene", "generic alkoxyallene", [I[1]], ["R3", "R4"])
add_generic(f"{DOC}:mol:generic_naphthalenone", "generic dearomatized naphthalenone", [I[1]], ["R1", "R2", "R4"])
add_generic(f"{DOC}:mol:generic_allyl_electrophile", "generic allylic alcohol/derivative", [I[1]], ["Ar", "R3"])
add_generic(f"{DOC}:mol:generic_alkyne", "generic alkyne", [I[1]], ["Ar"])
add_generic(f"{DOC}:mol:generic_allenamide", "generic allenamide", [I[1]], ["N-substituent"])
add_generic(f"{DOC}:mol:scheme2_naphthol", "Scheme 2 naphthol template 1", [I[3]], ["R1", "R2", "R3"])
add_generic(f"{DOC}:mol:scheme2_allene", "Scheme 2 alkoxyallene template 2", [I[3]], ["R4"])
add_generic(f"{DOC}:mol:scheme2_product", "Scheme 2 product template 3", [I[3]], ["R1", "R2", "R3", "R4"])
add_generic(f"{DOC}:mol:scheme3_naphthol", "Scheme 3 alpha-naphthol template 4", [I[4]], ["R1", "R2"])
add_generic(f"{DOC}:mol:scheme3_allene", "Scheme 3 alkoxyallene template 2", [I[4]], ["R3"])
add_generic(f"{DOC}:mol:scheme3_product", "Scheme 3 product template 5", [I[4]], ["R1", "R2", "R3"])

# Concrete recurring materials. CAS numbers are deliberately not populated.
add_mol(f"{DOC}:mol:1a", "1a", [I[2], I[3], I[4]], "1,3-dimethylnaphthalen-2-ol",
        "Cc1cc2ccccc2c(C)c1O")
add_mol(f"{DOC}:mol:2a", "2a", [I[2], I[3], I[4]], "benzyloxyallene", "C=C=COCc1ccccc1")
add_mol(f"{DOC}:mol:MeO-1a", "MeO-1a", [I[4]], "2-methoxy-1,3-dimethylnaphthalene",
        "COc1c(C)cc2ccccc2c1C")

agents = {
    "Pd(dba)2": [I[2]], "Pd2(dba)3": [I[2], I[3], I[4]], "Pd(OAc)2": [I[2]],
    "(R)-L1": [I[2]], "(R)-L2": [I[2]], "(R,R)-L3": [I[2]], "(R,R)-L4": [I[2]],
    "(R,R)-L5": [I[2], I[3], I[4]], "DBU": [I[2]], "Et3N": [I[2]], "TMG": [I[2], I[3], I[4]],
    "DCE": [I[2]], "THF": [I[2]], "PO": [I[2], I[3], I[4]], "NaBH4": [I[4], E[5]],
    "MeOH": [I[4], E[5]], "D2O": [I[4]], "EtOAc": [E[4], E[5]], "hexane": [E[4], E[5]],
    "petroleum ether": [E[5]],
}
for label, ev in agents.items():
    add_mol(f"{DOC}:mol:agent:{label.replace(' ', '_')}", label, ev, label)

# Product candidate SMILES. They encode the source E alkene but intentionally do
# not encode the quaternary tetrahedral center.
sm3 = {
    "3a": "CC1=Cc2ccccc2C(C)(C/C=C/OCc2ccccc2)C1=O",
    "3b": "CC1=Cc2ccccc2C(C)(C/C=C/OCc2ccc(F)cc2)C1=O",
    "3c": "CC1=Cc2ccccc2C(C)(C/C=C/OCc2ccc(Cl)cc2)C1=O",
    "3d": "CC1=Cc2ccccc2C(C)(C/C=C/OCc2ccc(Br)cc2)C1=O",
    "3e": "CC1=Cc2ccccc2C(C)(C/C=C/OCc2ccccc2Br)C1=O",
    "3f": "CC1=Cc2ccccc2C(C)(C/C=C/OCc2ccc(CC)cc2)C1=O",
    "3g": "CC1=Cc2ccccc2C(C)(C/C=C/OCc2ccc(C(C)(C)C)cc2)C1=O",
    "3h": "CC1=Cc2ccccc2C(C)(C/C=C/OCc2ccc(OC)cc2)C1=O",
    "3i": "CC1=Cc2ccccc2C(C)(C/C=C/OCc2ccc(OC(C)C)cc2)C1=O",
    "3j": "CC1=Cc2ccccc2C(C)(C/C=C/OCc2cc(OC)cc(OC)c2)C1=O",
    "3k": "CC1=Cc2ccccc2C(C)(C/C=C/OCc2ccc3ccccc3c2)C1=O",
    "3l": "CC1=Cc2ccccc2C(C)(C/C=C/OC(c2ccccc2)c2ccccc2)C1=O",
    "3m": "CC1=Cc2ccccc2C(CC)(C/C=C/OCc2ccccc2)C1=O",
    "3n": "CC1=Cc2ccccc2C(CCC)(C/C=C/OCc2ccccc2)C1=O",
    "3o": "CC1=Cc2ccccc2C(CCCC)(C/C=C/OCc2ccccc2)C1=O",
    "3p": "CC1=Cc2ccccc2C(CCCCCC)(C/C=C/OCc2ccccc2)C1=O",
    "3q": "CC1=Cc2ccccc2C(CC=C)(C/C=C/OCc2ccccc2)C1=O",
    "3r": "CCC1=Cc2ccccc2C(C)(C/C=C/OCc2ccccc2)C1=O",
    "3s": "C=CCC1=Cc2ccccc2C(C)(C/C=C/OCc2ccccc2)C1=O",
    "3t": "O=C1C(c2ccccc2)=Cc2ccccc2C1(C)C/C=C/OCc1ccccc1",
    "3u": "CC1=Cc2ccc(-c3ccccc3)cc2C(C)(C/C=C/OCc2ccccc2)C1=O",
}
sm5 = {
    "5a": "C=CCC1(C/C=C/OCc2ccccc2)C=Cc2ccccc2C1=O",
    "5b": "C=CCC1(C/C=C/OCc2ccc(CC)cc2)C=Cc2ccccc2C1=O",
    "5c": "C=CCC1(C/C=C/OCc2cc(OC)cc(OC)c2)C=Cc2ccccc2C1=O",
    "5d": "C=CCC1(C/C=C/OCc2ccc(F)cc2)C=Cc2ccccc2C1=O",
    "5e": "C=CCC1(C/C=C/OCc2ccccc2)C=C(OC)c2ccccc2C1=O",
    "5f": "C=CCC1(C/C=C/OCc2ccccc2)C=C(OCC)c2ccccc2C1=O",
}

page_for_3 = {x: 5 for x in ["3a", "3b", "3c", "3d", "3e", "3f", "3g", "3h", "3i"]}
page_for_3.update({x: 6 for x in ["3j", "3k", "3l", "3m", "3n", "3o", "3p", "3q", "3r", "3s"]})
for label, smi in sm3.items():
    ev = [I[3]] + ([E[page_for_3[label]]] if label in page_for_3 else []) + ([E[6]] if label == "3i" else [])
    add_mol(f"{DOC}:mol:{label}", label, ev, f"compound {label}", smi, True,
            "R at the new quaternary center" if label in page_for_3 else "configuration shown by Scheme 2 stereobond")
for label, smi in sm5.items():
    add_mol(f"{DOC}:mol:{label}", label, [I[4]], f"compound {label}", smi, True,
            "configuration shown by Scheme 3 stereobond")

# Concrete but not fully structure-transcribed scope inputs remain unresolved.
for label in [f"2{c}" for c in "bcdefghijkl"] + [f"1{c}" for c in "mnopqrstuvw"] + ["4a", "4e", "4f"]:
    if f"{DOC}:mol:{label}" not in M:
        add_mol(f"{DOC}:mol:{label}", label, [I[3] if label.startswith(("1", "2")) else I[4]], f"scope substrate {label}")

add_mol(f"{DOC}:mol:3ab", "3ab", [I[4], E[5]], "alcohol 3ab")
add_mol(f"{DOC}:mol:3a-D", "3a-D(H) mixture", [I[4]], "partially deuterated 3a mixture")


def pid(label):
    return label.replace("(", "").replace(")", "").replace(",", "").replace(" ", "_").replace("/", "_")


def participant(r, molecule_id, role, evidence, label, amount=None, mass=None, volume=None,
                equiv=None, molpct=None, ref=None, stage="reaction"):
    p = st.new_participant(f"p:{pid(label)}:{len(r['participants'])+1}", molecule_id, role, evidence, stage)
    if amount is not None:
        p["charges"]["amount"] = F(str(amount), amount, evidence)
    if mass is not None:
        p["charges"]["mass"] = F(str(mass), mass, evidence)
    if volume is not None:
        p["charges"]["volume"] = F(str(volume), volume, evidence)
    if equiv is not None:
        p["charges"]["equivalents"] = F(str(equiv), equiv, evidence)
    if molpct is not None:
        p["charges"]["mol_percent"] = F(str(molpct), molpct, evidence)
    if ref is not None:
        p["equivalence_reference"] = F(ref, ref, evidence)
    r["participants"].append(p)
    return p["participant_id"]


def reaction(suffix, evidence, record_type="specific_experiment", categories=None, location=None,
             name="Pd-catalyzed dearomative asymmetric allylic alkylation"):
    rid = f"{DOC}:{suffix}"
    r = st.new_reaction(rid, DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = F(categories or [], categories or [], evidence)
    r["location"] = F(location or suffix, location or suffix, evidence)
    r["reaction_name"] = F(name, name, evidence)
    r["transformation_type"] = F(
        None, "dearomative allylic C-alkylation" if "reduction" not in suffix else "ketone reduction",
        evidence, "inference", "Normalized from the source reaction description; not an author taxonomy.")
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": ["Machine extraction requires expert review; chemical accuracy is not verified."]}
    r["completeness"] = "partial"
    return r


def classify_dearomatization(r, evidence):
    changes = [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
         "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}],
         "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
        {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}],
         "before_order": 1, "after_order": 2, "change_type": "order_increased"},
        {"change_id": "b4", "atoms": [{"atom_map": None, "element": "O"}, {"atom_map": None, "element": "H"}],
         "before_order": 1, "after_order": 0, "change_type": "hydrogen_loss", "hydrogen_count": 1},
        {"change_id": "b5", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
         "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
    ]
    labels = [
        {**st.bond_label("C", "C", 1), "supporting_change_ids": ["b1"]},
        {**st.bond_label("C", "H", 1), "supporting_change_ids": ["b2"]},
        {**st.bond_label("C", "O", 2), "supporting_change_ids": ["b3"]},
    ]
    deriv = "Partial net-change interpretation of the source reactant/product drawings; no complete atom mapping is asserted."
    r["classification"]["bond_changes"] = [F(None, c, evidence, "structure_derivation", deriv) for c in changes]
    r["classification"]["labels"] = F(None, labels, evidence, "structure_derivation", deriv)
    r["classification"]["coverage"] = F(None, "partial", evidence, "inference", deriv)
    r["classification"]["limitations"] = [
        "No complete atom-mapped reaction is asserted.",
        "Aromatic bond reorganization, palladium intermediates, and coproduct identities are not exhaustively mapped.",
        "The C-H label records net hydrogen transfer implicit in the product graph; its elementary source is not independently verified.",
    ]


def classify_none(r, evidence, reason):
    r["classification"]["labels"] = F(None, [], evidence, "inference", reason)
    r["classification"]["coverage"] = F(None, "not_applicable", evidence, "inference", reason)
    r["classification"]["limitations"] = [reason]


def classify_unresolved(r, evidence, reason):
    r["classification"]["labels"] = F(None, None, evidence)
    r["classification"]["coverage"] = F(None, "unresolved", evidence, "inference", reason)
    r["classification"]["limitations"] = [reason]


def classify_reduction(r, evidence):
    changes = [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}],
         "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}],
         "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
        {"change_id": "b3", "atoms": [{"atom_map": None, "element": "O"}, {"atom_map": None, "element": "H"}],
         "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
    ]
    labels = [
        {**st.bond_label("C", "H", 1), "supporting_change_ids": ["b2"]},
        {**st.bond_label("O", "H", 1), "supporting_change_ids": ["b3"]},
    ]
    deriv = "Partial net reduction of the ketone drawing to the alcohol; hydride/proton atom mapping is not asserted."
    r["classification"]["bond_changes"] = [F(None, c, evidence, "structure_derivation", deriv) for c in changes]
    r["classification"]["labels"] = F(None, labels, evidence, "structure_derivation", deriv)
    r["classification"]["coverage"] = F(None, "partial", evidence, "inference", deriv)
    r["classification"]["limitations"] = ["No complete atom mapping; NaBH4/MeOH hydrogen sources are not assigned atom-by-atom."]


def add_step(r, sid, stage, participant_ids, evidence, operation=None, temp=None, time=None,
             atmosphere=None, pressure=None):
    s = st.new_step(sid, evidence, stage)
    s["participant_ids"] = participant_ids
    if operation is not None:
        s["operation"] = F(operation, operation, evidence)
    if temp is not None:
        s["temperature"] = F(str(temp), temp, evidence)
    if time is not None:
        s["time"] = F(str(time), time, evidence)
    if atmosphere is not None:
        s["atmosphere"] = F(atmosphere, atmosphere, evidence)
    if pressure is not None:
        s["pressure"] = F(pressure, pressure, evidence)
    r["steps"].append(s)


def outcome(r, oid, kind, target, evidence, raw, normalized, basis=None, yield_type=None, scope=None,
            method="original_extraction", derivation=None, status=None):
    o = st.new_outcome(oid, kind, [target] if target else [], evidence)
    o["value"] = F(raw, normalized, evidence, method, derivation, status)
    if basis is not None:
        o["measurement_basis"] = F(basis, basis, evidence)
    if yield_type is not None:
        o["yield_type"] = F(yield_type, yield_type, evidence)
    if scope is not None:
        o["scope"] = F(scope, scope, evidence)
    r["outcomes"].append(o)


R = []

# Page 1: three cited prior-work arrows and one current-work generic arrow.
for suffix, electrophile, catalyst, finding in [
    ("scheme1a_prior_allylic_derivatives", "generic_allyl_electrophile", "[Pd]/[Ir]/[Rh]", "previous work: allylic alcohol derivatives"),
    ("scheme1a_prior_alkyne", "generic_alkyne", "[Pd]/oxidant or without oxidant", "previous work: atom-economic alkyne coupling"),
    ("scheme1a_prior_allenamide", "generic_allenamide", "[Au] or chiral Brønsted acid", "previous work: allenamides"),
]:
    r = reaction(suffix, [I[1]], "literature_example", [], "Scheme 1a", finding)
    participant(r, M[f"{DOC}:mol:generic_naphthol"]["molecule_id"], "reactant", [I[1]], "generic naphthol")
    participant(r, M[f"{DOC}:mol:{electrophile}"]["molecule_id"], "reactant", [I[1]], electrophile)
    participant(r, M[f"{DOC}:mol:generic_naphthalenone"]["molecule_id"], "product", [I[1]], "generic product")
    r["notes"].append(F(catalyst, catalyst, [I[1]]))
    classify_unresolved(r, [I[1]], "Generic cited prior work is not classified as a current-paper experiment.")
    R.append(r)

r = reaction("scheme1b_current_generic", [I[1]], "generic_scheme", [], "Scheme 1b", "this work: naphthols with alkoxyallenes")
participant(r, M[f"{DOC}:mol:generic_naphthol"]["molecule_id"], "reactant", [I[1]], "generic naphthol")
participant(r, M[f"{DOC}:mol:generic_alkoxyallene"]["molecule_id"], "reactant", [I[1]], "generic alkoxyallene")
participant(r, M[f"{DOC}:mol:generic_naphthalenone"]["molecule_id"], "product", [I[1]], "generic product")
classify_dearomatization(r, [I[1]])
R.append(r)

# Page 2: every optimization row is a distinct experiment.
table = [
    (1, "(R)-L1", "Pd(dba)2", "DBU", "DCE", 60, 64, "55:45"),
    (2, "(R)-L2", "Pd(dba)2", "DBU", "DCE", 60, 36, "56:44"),
    (3, "(R,R)-L3", "Pd(dba)2", "DBU", "DCE", 60, 58, "78:22"),
    (4, "(R,R)-L4", "Pd(dba)2", "DBU", "DCE", 60, None, None),
    (5, "(R,R)-L5", "Pd(dba)2", "DBU", "DCE", 60, 63, "82:18"),
    (6, "(R,R)-L5", "Pd2(dba)3", "DBU", "DCE", 60, 66, "89:11"),
    (7, "(R,R)-L5", "Pd(OAc)2", "DBU", "DCE", 60, 61, "77:23"),
    (8, "(R,R)-L5", "Pd2(dba)3", "DBU", "THF", 60, 25, "91:9"),
    (9, "(R,R)-L5", "Pd2(dba)3", "Et3N", "THF", 60, 35, "60:40"),
    (10, "(R,R)-L5", "Pd2(dba)3", "TMG", "THF", 60, 34, "91.5:8.5"),
    (11, "(R,R)-L5", "Pd2(dba)3", "TMG", "PO", 60, 81, "92:8"),
    (12, "(R,R)-L5", "Pd2(dba)3", "TMG", "PO", 25, None, None),
]
for entry, ligand, catalyst, base, solvent, temp, yld, er in table:
    ev = [I[2], E[2]]
    r = reaction(f"table1_entry{entry:02d}", ev, categories=["optimization"], location=f"Table 1 entry {entry}")
    a = participant(r, f"{DOC}:mol:1a", "reactant", ev, "1a", amount={"number": 0.2, "unit": "mmol"})
    b = participant(r, f"{DOC}:mol:2a", "reactant", ev, "2a", amount={"number": 0.3, "unit": "mmol"},
                    equiv={"number": 1.5, "unit": "equiv"}, ref=f"{DOC}:mol:1a")
    product = participant(r, f"{DOC}:mol:3a", "product", ev, "3a")
    c = participant(r, f"{DOC}:mol:agent:{catalyst}", "catalyst", ev, catalyst,
                    molpct={"number": 2.5, "unit": "mol%"}, ref=f"{DOC}:mol:1a")
    l = participant(r, f"{DOC}:mol:agent:{ligand}", "ligand", ev, ligand,
                    molpct={"number": 5, "unit": "mol%"}, ref=f"{DOC}:mol:1a")
    q = participant(r, f"{DOC}:mol:agent:{base}", "reagent", ev, base,
                    amount={"number": 0.2, "unit": "mmol"}, equiv={"number": 1, "unit": "equiv"},
                    ref=f"{DOC}:mol:1a")
    sol = participant(r, f"{DOC}:mol:agent:{solvent}", "solvent", ev, solvent,
                      volume={"number": 2, "unit": "mL"})
    r["reaction_scale"] = F("1a (0.2 mmol)", {"reference_molecule_id": f"{DOC}:mol:1a", "amount": {"number": 0.2, "unit": "mmol"}}, ev)
    r["solvent_conditions"]["composition"] = F(solvent, [solvent], ev)
    r["solvent_conditions"]["total_volume"] = F("2 mL", {"number": 2, "unit": "mL"}, ev)
    if entry == 11:
        r["solvent_conditions"]["concentration"] = F("PO, 0.1 M", {"number": 0.1, "unit": "M"}, ev)
        conflict = F(["TMG (table)", "DBU (prose)"], None, ev, status="conflict")
        conflict["alternatives"] = [F("TMG (Table 1 entry 11)", "TMG", [I[2]]), F("DBU (optimization prose)", "DBU", [E[2]])]
        r["notes"].append(conflict)
        r["review"]["unresolved"].append("Entry 11 base conflicts: table reports TMG, prose reports DBU; neither is selected as the corrected value.")
    add_step(r, "reaction", "reaction", [a, b, c, l, q, sol, product], ev, "stir under argon",
             {"number": temp, "unit": "°C"}, {"number": 12, "unit": "h"}, "argon")
    r["procedure"]["original_text"] = F("Table 1 general conditions and row-specific cells", "Table 1 general conditions and row-specific cells", ev)
    if yld is not None:
        outcome(r, "yield", "yield", f"{DOC}:mol:3a", ev, f"{yld}%", {"number": yld, "unit": "%"},
                "isolated after silica gel chromatography", "isolated")
        outcome(r, "er", "selectivity", f"{DOC}:mol:3a", ev, f"{er} e.r.",
                {"metric": "e.r.", "ratio": er, "order": "major:minor"}, "chiral stationary-phase HPLC", scope="enantiomer ratio")
        classify_dearomatization(r, ev)
    else:
        r["notes"].append(F("yield and e.r. cells are blank", "yield and e.r. cells are blank; not interpreted as no reaction", [I[2]]))
        r["review"]["unresolved"].append("Blank Table 1 yield/e.r. cells do not establish no reaction or a zero value.")
        classify_unresolved(r, ev, "The table row has blank outcome cells, so an observed transformation is not asserted.")
    if entry == 11:
        outcome(r, "ez", "selectivity", f"{DOC}:mol:3a", ev, "E/Z > 19:1",
                {"metric": "E/Z", "relation": ">", "ratio": "19:1", "major": "E"}, None, scope="geometric isomer ratio")
    R.append(r)

# Scheme 2 generic arrow.
r = reaction("scheme2_generic", [I[3]], "generic_scheme", [], "Scheme 2 generic arrow")
for mid, role, label in [("scheme2_naphthol", "reactant", "1"), ("scheme2_allene", "reactant", "2"), ("scheme2_product", "product", "3")]:
    participant(r, f"{DOC}:mol:{mid}", role, [I[3]], label)
classify_dearomatization(r, [I[3]])
R.append(r)

yields3 = {
    "3a": (81, "92:8"), "3b": (83, "92:8"), "3c": (92, "91:9"), "3d": (75, "93.5:6.5"),
    "3e": (60, "92.5:7.5"), "3f": (70, "93:7"), "3g": (80, "92.5:7.5"), "3h": (95, "92:8"),
    "3i": (74, "93.5:6.5"), "3j": (85, "93.5:6.5"), "3k": (76, "91:9"), "3l": (82, "93:7"),
    "3m": (76, "93.5:6.5"), "3n": (65, "94.5:5.5"), "3o": (61, "94.5:5.5"), "3p": (80, "94.5:5.5"),
    "3q": (65, "93.5:6.5"), "3r": (76, "93.5:6.5"), "3s": (89, "93:8"), "3t": (75, "92.5:7.5"),
    "3u": (71, "91.5:8.5"),
}
masses3 = {"3a": 51.5, "3b": 55.9, "3c": 64.9, "3d": 59.6, "3e": 47.6, "3f": 48.6,
           "3g": 60.0, "3h": 66.3, "3i": 55.8, "3j": 64.4, "3k": 56.1, "3l": 64.8,
           "3m": 50.6, "3n": 45.1, "3o": 44.0, "3p": 62.2, "3q": 44.9, "3r": 50.6, "3s": 61.4}
hplc3 = {
    "3a": "Chiralcel OD, 95:5 hexanes/i-PrOH, 1 mL/min; tR major 9.41 min, minor 8.72 min",
    "3b": "Chiralcel OD, 99:1 hexanes/i-PrOH, 1 mL/min; tR major 15.86 min, minor 17.91 min",
    "3c": "Chiralcel AD, 99:1 hexanes/i-PrOH, 1 mL/min; tR major 20.25 min, minor 16.70 min",
    "3d": "Chiralcel AD, 95:5 hexanes/i-PrOH, 1 mL/min; tR major 9.60 min, minor 8.41 min",
    "3e": "Chiralcel IC, 90:10 hexanes/i-PrOH, 1 mL/min; tR major 8.10 min, minor 9.52 min",
    "3f": "Chiralcel AD, 99:1 hexanes/i-PrOH, 1 mL/min; tR major 15.40 min, minor 14.13 min",
    "3g": "Chiralcel OD, 99:1 hexanes/i-PrOH, 1 mL/min; tR major 11.15 min, minor 13.09 min",
    "3h": "Chiralcel OD, 99:1 hexanes/i-PrOH, 1 mL/min; tR major 28.29 min, minor 26.00 min",
    "3i": "Chiralcel IC, 95:5 hexanes/i-PrOH, 0.5 mL/min; tR major 42.96 min, minor 45.81 min",
    "3j": "Chiralcel IA, 90:10 hexanes/i-PrOH, 1 mL/min; tR major 7.36 min, minor 6.92 min",
    "3k": "Chiralcel AD, 90:10 hexanes/i-PrOH, 1 mL/min; tR major 8.03 min, minor 6.92 min",
    "3l": "Chiralcel IC, 90:10 hexanes/i-PrOH, 1 mL/min; tR major 8.24 min, minor 9.99 min",
    "3m": "Chiralcel OD, 95:5 hexanes/i-PrOH, 1 mL/min; tR major 8.85 min, minor 8.18 min",
    "3n": "Chiralcel IC, 95:5 hexanes/i-PrOH, 1 mL/min; tR major 10.84 min, minor 11.77 min",
    "3o": "Chiralcel IC, 90:10 hexanes/i-PrOH, 1 mL/min; tR major 9.34 min, minor 10.26 min",
    "3p": "Chiralcel AD, 95:5 hexanes/i-PrOH, 1 mL/min; tR major 7.88 min, minor 7.22 min",
    "3q": "Chiralcel AD, 99:1 hexanes/i-PrOH, 1 mL/min; tR major 9.64 min, minor 9.03 min",
    "3r": "Chiralcel IB, 95:5 hexanes/i-PrOH, 1 mL/min; tR major 5.92 min, minor 5.47 min",
    "3s": "Chiralcel OD, 90:10 hexanes/i-PrOH, 1 mL/min; tR major 8.11 min, minor 7.23 min",
}

for label, (yld, er) in yields3.items():
    letter = label[1:]
    input1 = "1a" if label in [f"3{c}" for c in "abcdefghijkl"] else f"1{letter}"
    input2 = f"2{letter}" if label in [f"3{c}" for c in "bcdefghijkl"] else "2a"
    if input1 not in ["1a"] and f"{DOC}:mol:{input1}" not in M:
        add_mol(f"{DOC}:mol:{input1}", input1, [I[3]], f"scope substrate {input1}")
    ev = [I[3], E[4]] + ([E[page_for_3[label]]] if label in page_for_3 else []) + ([E[6]] if label == "3i" else [])
    r = reaction(f"scheme2_{label}", ev, categories=["substrate_scope"], location=f"Scheme 2, {label}")
    r["general_scheme_ids"] = [f"{DOC}:scheme2_generic"]
    s1 = participant(r, f"{DOC}:mol:{input1}", "reactant", ev, input1, amount={"number": 0.2, "unit": "mmol"})
    s2 = participant(r, f"{DOC}:mol:{input2}", "reactant", ev, input2, amount={"number": 0.3, "unit": "mmol"},
                     equiv={"number": 1.5, "unit": "equiv"}, ref=f"{DOC}:mol:{input1}")
    prod = participant(r, f"{DOC}:mol:{label}", "product", ev, label)
    cat = participant(r, f"{DOC}:mol:agent:Pd2(dba)3", "catalyst", ev, "Pd2(dba)3",
                      molpct={"number": 2.5, "unit": "mol%"}, ref=f"{DOC}:mol:{input1}")
    lig = participant(r, f"{DOC}:mol:agent:(R,R)-L5", "ligand", ev, "(R,R)-L5",
                      molpct={"number": 5, "unit": "mol%"}, ref=f"{DOC}:mol:{input1}")
    base = participant(r, f"{DOC}:mol:agent:TMG", "reagent", ev, "TMG", amount={"number": 0.2, "unit": "mmol"},
                       equiv={"number": 1, "unit": "equiv"}, ref=f"{DOC}:mol:{input1}")
    sol = participant(r, f"{DOC}:mol:agent:PO", "solvent", ev, "PO", volume={"number": 2, "unit": "mL"})
    etoac = participant(r, f"{DOC}:mol:agent:EtOAc", "purification_agent", ev, "EtOAc", stage="purification")
    hx = participant(r, f"{DOC}:mol:agent:hexane", "purification_agent", ev, "hexane", stage="purification")
    r["reaction_scale"] = F(f"{input1} (0.2 mmol)", {"reference_molecule_id": f"{DOC}:mol:{input1}", "amount": {"number": 0.2, "unit": "mmol"}}, ev)
    r["solvent_conditions"]["composition"] = F("propylene oxide (PO)", ["PO"], ev)
    r["solvent_conditions"]["total_volume"] = F("2 mL", {"number": 2, "unit": "mL"}, ev)
    r["solvent_conditions"]["concentration"] = F("0.1 M", {"number": 0.1, "unit": "M"}, [E[2]])
    r["solvent_conditions"]["concentration_basis"] = F("relative to 0.2 mmol naphthol in 2 mL PO", f"{input1}", [E[2]])
    add_step(r, "reaction", "reaction", [s1, s2, cat, lig, base, sol, prod], ev,
             "stir in a 5 mL sealed tube until TLC completion; Scheme 2 says overnight",
             {"number": 60, "unit": "°C"}, "overnight", "noble gas protection")
    add_step(r, "concentrate", "workup", [prod], [E[4]], "concentrate under reduced pressure", pressure="reduced pressure")
    add_step(r, "chromatography", "purification", [prod, etoac, hx], [E[4]], "flash chromatography (EtOAc/hexane)")
    r["procedure"]["original_text"] = F("General Procedure for the Dearomative Asymmetric Allylic Reaction of beta-Naphthols with Alkoxyallenes", "sealed tube; PO; catalyst/ligand/TMG; 60 °C; TLC; concentrate; flash chromatography", [E[4]])
    r["procedure"]["general_procedure_reference"] = F("beta-naphthol general procedure", "beta-naphthol general procedure", [E[4]])
    r["procedure"]["addition_order"] = F("substrates in PO, followed by Pd2(dba)3, TMG, and (R,R)-L5", [input1, input2, "PO", "Pd2(dba)3", "TMG", "(R,R)-L5"], [E[4]])
    outcome(r, "yield", "yield", f"{DOC}:mol:{label}", ev, f"{yld}%", {"number": yld, "unit": "%"},
            "isolated product", "isolated")
    er_basis = hplc3.get(label, "chiral HPLC analysis")
    outcome(r, "er", "selectivity", f"{DOC}:mol:{label}", ev, f"{er} e.r.",
            {"metric": "e.r.", "ratio": er, "order": "major:minor"}, er_basis, scope="enantiomer ratio")
    if label in masses3:
        outcome(r, "mass", "product_amount", f"{DOC}:mol:{label}", [E[page_for_3[label]]], f"{masses3[label]} mg",
                {"number": masses3[label], "unit": "mg"}, "isolated product after chromatography")
    if label == "3s":
        r["notes"].append(F("93:8 e.r.", "source ratio sums to 101 and is preserved without correction", [I[3], E[6]]))
        r["review"]["unresolved"].append("The source reports 3s as 93:8 e.r.; the ratio is preserved exactly despite summing to 101.")
    classify_dearomatization(r, ev)
    R.append(r)

# Scheme 3 generic arrow and six alpha-naphthol scope rows.
r = reaction("scheme3_generic", [I[4]], "generic_scheme", [], "Scheme 3 generic arrow")
for mid, role, label in [("scheme3_naphthol", "reactant", "4"), ("scheme3_allene", "reactant", "2"), ("scheme3_product", "product", "5")]:
    participant(r, f"{DOC}:mol:{mid}", role, [I[4]], label)
classify_dearomatization(r, [I[4]])
R.append(r)

yields5 = {"5a": (66, "91.5:8.5", "4a", "2a"), "5b": (62, "92:8", "4a", "2f"),
           "5c": (71, "92:8", "4a", "2j"), "5d": (65, "91.5:8.5", "4a", "2b"),
           "5e": (79, "94:6", "4e", "2a"), "5f": (79, "93:7", "4f", "2a")}
for label, (yld, er, sub, allene) in yields5.items():
    ev = [I[4], E[4], E[5]]
    r = reaction(f"scheme3_{label}", ev, categories=["substrate_scope"], location=f"Scheme 3, {label}")
    r["general_scheme_ids"] = [f"{DOC}:scheme3_generic"]
    p1 = participant(r, f"{DOC}:mol:{sub}", "reactant", ev, sub, amount={"number": 0.2, "unit": "mmol"})
    p2 = participant(r, f"{DOC}:mol:{allene}", "reactant", ev, allene, amount={"number": 0.3, "unit": "mmol"},
                     equiv={"number": 1.5, "unit": "equiv"}, ref=f"{DOC}:mol:{sub}")
    pp = participant(r, f"{DOC}:mol:{label}", "product", ev, label)
    pc = participant(r, f"{DOC}:mol:agent:Pd2(dba)3", "catalyst", ev, "Pd2(dba)3",
                     molpct={"number": 2.5, "unit": "mol%"}, ref=f"{DOC}:mol:{sub}")
    pl = participant(r, f"{DOC}:mol:agent:(R,R)-L5", "ligand", ev, "(R,R)-L5",
                     molpct={"number": 5, "unit": "mol%"}, ref=f"{DOC}:mol:{sub}")
    pb = participant(r, f"{DOC}:mol:agent:TMG", "reagent", ev, "TMG", amount={"number": 0.2, "unit": "mmol"},
                     equiv={"number": 1, "unit": "equiv"}, ref=f"{DOC}:mol:{sub}")
    ps = participant(r, f"{DOC}:mol:agent:PO", "solvent", ev, "PO", volume={"number": 2, "unit": "mL"})
    pe = participant(r, f"{DOC}:mol:agent:EtOAc", "purification_agent", ev, "EtOAc", stage="purification")
    ph = participant(r, f"{DOC}:mol:agent:hexane", "purification_agent", ev, "hexane", stage="purification")
    r["reaction_scale"] = F("Scheme 3 footnote says 1 (0.2 mmol); arrow/general procedure identify substrate 4",
                            {"reference_molecule_id": f"{DOC}:mol:{sub}", "amount": {"number": 0.2, "unit": "mmol"}}, ev)
    r["solvent_conditions"]["composition"] = F("PO", ["PO"], ev)
    r["solvent_conditions"]["total_volume"] = F("2 mL", {"number": 2, "unit": "mL"}, ev)
    add_step(r, "reaction", "reaction", [p1, p2, pc, pl, pb, ps, pp], ev,
             "stir in a 5 mL sealed tube until TLC completion; Scheme 3 says overnight",
             {"number": 60, "unit": "°C"}, "overnight", "noble gas protection")
    add_step(r, "concentrate", "workup", [pp], [E[5]], "concentrate under reduced pressure", pressure="reduced pressure")
    add_step(r, "chromatography", "purification", [pp, pe, ph], [E[5]], "flash chromatography (EtOAc/hexane)")
    r["procedure"]["original_text"] = F("General Procedure for the Dearomative Asymmetric Allylic Reaction of alpha-Naphthols with Alkoxyallenes", "sealed tube; PO; catalyst/ligand/TMG; 60 °C; TLC; concentrate; flash chromatography", [E[4], E[5]])
    r["procedure"]["general_procedure_reference"] = F("alpha-naphthol general procedure", "alpha-naphthol general procedure", [E[4], E[5]])
    r["notes"].append(F("Scheme 3 footnote labels the 0.2 mmol substrate as 1, while the arrow and general procedure label it 4", "label conflict retained", [I[4], E[5]]))
    r["review"]["unresolved"].append("Scheme 3 footnote calls the substrate 1; arrow/general procedure call it 4. Quantity is linked to the row-specific alpha-naphthol without silently correcting the source text.")
    outcome(r, "yield", "yield", f"{DOC}:mol:{label}", ev, f"{yld}%", {"number": yld, "unit": "%"}, "isolated product", "isolated")
    outcome(r, "er", "selectivity", f"{DOC}:mol:{label}", ev, f"{er} e.r.", {"metric": "e.r.", "ratio": er, "order": "major:minor"}, "chiral HPLC analysis", scope="enantiomer ratio")
    classify_dearomatization(r, ev)
    R.append(r)

# Scheme 4a scale-up, preserving the 1.81 g prose / 1.44 g scheme conflict.
ev = [E[3], I[4]]
r = reaction("scheme4a_scaleup_3a", ev, categories=["scale_up"], location="Scheme 4a and page-3 prose")
p1 = participant(r, f"{DOC}:mol:1a", "reactant", ev, "1a", amount={"number": 5.8, "unit": "mmol"}, mass={"number": 1.0, "unit": "g"})
p2 = participant(r, f"{DOC}:mol:2a", "reactant", ev, "2a")
pp = participant(r, f"{DOC}:mol:3a", "product", ev, "3a")
r["reaction_scale"] = F("1a (1.0 g, 5.8 mmol)", {"reference_molecule_id": f"{DOC}:mol:1a", "amount": {"number": 5.8, "unit": "mmol"}}, ev)
add_step(r, "reaction", "reaction", [p1, p2, pp], ev, "standard conditions")
outcome(r, "yield", "yield", f"{DOC}:mol:3a", ev, "78%", {"number": 78, "unit": "%"})
outcome(r, "er", "selectivity", f"{DOC}:mol:3a", ev, "92:8 e.r.", {"metric": "e.r.", "ratio": "92:8", "order": "major:minor"}, scope="enantiomer ratio")
o = st.new_outcome("product_mass_conflict", "product_amount", [f"{DOC}:mol:3a"], ev)
o["value"] = F(["1.81 g in prose", "1.44 g in Scheme 4a"], None, ev, status="conflict")
o["value"]["alternatives"] = [F("1.81 g of 3a", {"number": 1.81, "unit": "g"}, [E[3]]), F("3a (1.44 g, 78% yield)", {"number": 1.44, "unit": "g"}, [I[4]])]
o["scope"] = F("preparative-scale isolated product mass", "preparative-scale isolated product mass", ev)
r["outcomes"].append(o)
r["notes"].append(F("The prose and Scheme 4a disagree on isolated mass; both values are retained.", "source conflict", ev))
r["review"]["unresolved"].append("Scale-up product mass conflict: 1.81 g in page-3 prose versus 1.44 g in Scheme 4a.")
r["review"]["unresolved"].append("The 78% scale-up yield measurement basis is not explicitly reported in the bounded pages.")
classify_dearomatization(r, ev)
R.append(r)

# Scheme 4b ketone reduction, with the detailed page-5 procedure.
ev = [I[4], E[5]]
r = reaction("scheme4b_reduction_3ab", ev, categories=["derivatization"], location="Scheme 4b and reduction procedure", name="diastereoselective reduction of 3a")
p3a = participant(r, f"{DOC}:mol:3a", "reactant", ev, "3a", amount={"number": 0.1, "unit": "mmol"}, mass={"number": 31.8, "unit": "mg"})
pbh = participant(r, f"{DOC}:mol:agent:NaBH4", "reagent", ev, "NaBH4", amount={"number": 0.12, "unit": "mmol"}, mass={"number": 4.9, "unit": "mg"})
pme = participant(r, f"{DOC}:mol:agent:MeOH", "solvent", ev, "MeOH", volume={"number": 2.0, "unit": "mL"})
prod = participant(r, f"{DOC}:mol:3ab", "product", ev, "3ab")
ppe = participant(r, f"{DOC}:mol:agent:petroleum_ether", "purification_agent", ev, "petroleum ether", stage="purification")
pet = participant(r, f"{DOC}:mol:agent:EtOAc", "purification_agent", ev, "EtOAc", stage="purification")
r["reaction_scale"] = F("3a (31.8 mg, 0.1 mmol, 92:8 e.r.)", {"reference_molecule_id": f"{DOC}:mol:3a", "amount": {"number": 0.1, "unit": "mmol"}}, ev)
r["solvent_conditions"]["composition"] = F("anhydrous MeOH", ["anhydrous MeOH"], ev)
r["solvent_conditions"]["total_volume"] = F("2.0 mL", {"number": 2.0, "unit": "mL"}, ev)
add_step(r, "reaction", "reaction", [p3a, pbh, pme, prod], ev, "add NaBH4 to 3a in anhydrous MeOH; stir until TLC completion", "rt", {"number": 5, "unit": "h"})
add_step(r, "remove_solvent", "workup", [prod], ev, "remove solvents under reduced pressure", pressure="reduced pressure")
add_step(r, "crude_dr", "characterization", [prod], ev, "determine crude diastereomeric ratio by 1H NMR")
add_step(r, "chromatography", "purification", [prod, ppe, pet], ev, "silica gel chromatography, PE/EtOAc = 9/1")
r["procedure"]["original_text"] = F("Procedure for the Synthesis of Alcohol 6a [product called 3ab in the procedure]", "NaBH4 reduction procedure affording 3ab", ev)
r["notes"].append(F("Heading says Alcohol 6a, while the text and Scheme 4b identify product 3ab", "source label inconsistency retained", [E[5], I[4]]))
r["review"]["unresolved"].append("Reduction procedure heading calls the alcohol 6a, while the body and scheme call it 3ab.")
outcome(r, "yield", "yield", f"{DOC}:mol:3ab", ev, "95%", {"number": 95, "unit": "%"}, "isolated after silica gel chromatography", "isolated")
outcome(r, "er", "selectivity", f"{DOC}:mol:3ab", ev, "91.5:8.5 e.r.", {"metric": "e.r.", "ratio": "91.5:8.5", "order": "major:minor"}, None, scope="enantiomer ratio")
outcome(r, "dr", "selectivity", f"{DOC}:mol:3ab", ev, ">95:5 d.r.", {"metric": "d.r.", "relation": ">", "ratio": "95:5", "order": "major:minor"}, "1H NMR of crude reaction mixture", scope="diastereomer ratio")
classify_reduction(r, ev)
R.append(r)

# Scheme 4c: real control without base.
ev = [E[3], I[4]]
r = reaction("scheme4c_control_no_base", ev, categories=["control", "mechanistic_study"], location="Scheme 4c")
p1 = participant(r, f"{DOC}:mol:1a", "reactant", ev, "1a")
p2 = participant(r, f"{DOC}:mol:2a", "reactant", ev, "2a")
pp = participant(r, f"{DOC}:mol:3a", "product", ev, "3a")
add_step(r, "reaction", "reaction", [p1, p2, pp], ev, "standard conditions without base", None, {"number": 12, "unit": "h"})
outcome(r, "yield", "yield", f"{DOC}:mol:3a", ev, "31% yield", {"number": 31, "unit": "%"})
outcome(r, "er", "selectivity", f"{DOC}:mol:3a", ev, "60:40 e.r.", {"metric": "e.r.", "ratio": "60:40", "order": "major:minor"}, scope="enantiomer ratio")
r["review"]["unresolved"].append("The 31% control yield measurement basis is not explicitly reported.")
classify_dearomatization(r, ev)
R.append(r)

# Scheme 4d: true negative control; no product is invented.
ev = [E[3], I[4]]
r = reaction("scheme4d_control_methyl_ether", ev, categories=["control", "mechanistic_study"], location="Scheme 4d")
p1 = participant(r, f"{DOC}:mol:MeO-1a", "reactant", ev, "MeO-1a")
p2 = participant(r, f"{DOC}:mol:2a", "reactant", ev, "2a")
add_step(r, "reaction", "reaction", [p1, p2], ev, "standard conditions")
outcome(r, "no_reaction", "no_reaction", None, ev, "no reaction", "no reaction")
classify_none(r, ev, "Source reports no reaction; intended product bonds are not classified as observed chemistry.")
R.append(r)

# Scheme 4e: deuterium-labeling mechanistic experiment.
ev = [E[3], I[4]]
r = reaction("scheme4e_deuterium_labeling", ev, categories=["mechanistic_study"], location="Scheme 4e")
p1 = participant(r, f"{DOC}:mol:1a", "reactant", ev, "1a")
p2 = participant(r, f"{DOC}:mol:2a", "reactant", ev, "2a")
pd2o = participant(r, f"{DOC}:mol:agent:D2O", "additive", ev, "D2O")
next(p for p in r["participants"] if p["participant_id"] == pd2o)["charges"]["volume"] = F(
    "1 drop", "one drop", ev)
pp = participant(r, f"{DOC}:mol:3a-D", "product", ev, "3a-D(H) mixture")
add_step(r, "reaction", "reaction", [p1, p2, pd2o, pp], ev, "standard conditions with D2O (1 drop)")
outcome(r, "yield", "yield", f"{DOC}:mol:3a-D", ev, "72% yield", {"number": 72, "unit": "%"})
outcome(r, "dh", "selectivity", f"{DOC}:mol:3a-D", ev, "D:H = 61:39", {"metric": "D:H", "ratio": "61:39", "site": "middle carbon of allene-derived fragment"}, None, scope="isotopic incorporation ratio")
r["review"]["unresolved"].append("The 72% labeling-experiment yield basis and the D:H measurement method are not explicitly reported in the bounded pages.")
classify_dearomatization(r, ev)
r["classification"]["limitations"].append("The hydrogen-gain change includes H/D isotopic incorporation; the partial classification does not atom-map isotope sources.")
R.append(r)

# Scheme 5 is an author-proposed cycle, not an observed experiment.
ev = [E[3], I[4]]
r = reaction("scheme5_proposed_cycle", ev, "proposed_mechanism", ["mechanistic_study"], "Scheme 5", "proposed catalytic cycle")
r["notes"].append(F("oxidative addition to naphthol; hydropalladation of allene; nucleophilic addition to pi-allyl; reductive elimination", "author-proposed sequence", ev))
classify_none(r, ev, "Proposed catalytic cycle is not an observed net synthesis experiment.")
R.append(r)

bundle = json.loads((RUN / "submitted.json").read_text())
bundle["molecules"] = list(M.values())
bundle["reactions"] = R
bundle["source_context"] = [
    {"context_id": f"{DOC}:context:title", "kind": "bibliographic_metadata", "evidence_ids": [I[1]],
     "value": "Pd-Catalyzed Dearomative Asymmetric Allylic Alkylation of Naphthols with Alkoxyallenes"},
    {"context_id": f"{DOC}:context:general_information", "kind": "experimental_context", "evidence_ids": [E[4]],
     "value": "Analytical and characterization instrumentation, chromatography practice, and preparation of racemic HPLC references."},
    {"context_id": f"{DOC}:context:characterization_page5", "kind": "characterization", "evidence_ids": [E[5], I[5]],
     "value": "Preparation/characterization paragraphs for 3a-3i; masses, yields, and HPLC methods are merged into the matching Scheme 2 records. Spectra/HRMS are not separate reactions."},
    {"context_id": f"{DOC}:context:characterization_page6", "kind": "characterization", "evidence_ids": [E[6], I[6]],
     "value": "Preparation/characterization paragraphs for 3j-3s; masses, yields, and HPLC methods are merged into the matching Scheme 2 records. Spectra/HRMS are not separate reactions."},
    {"context_id": f"{DOC}:context:absolute_configuration", "kind": "stereochemistry_discussion", "evidence_ids": [E[2]],
     "value": "The authors assign R configuration to 3a by calculated ECD comparison and refer to SI Figure S1; the SI itself was not provided or associated."},
]

st.write(RUN / "submitted.json", bundle)
checks = {}
for mid, m in M.items():
    smi = st.value(m["smiles"])
    if smi:
        checks[mid] = st.stereo_inspect(smi)
st.write(RUN / "raw_model/stereo_checks.json", checks)
dataset = json.loads((RUN / "dataset.json").read_text())
st.write(RUN / "raw_model/source_mentions_key_labels.json",
         st.source_mentions(dataset, ["1a", "2a", "3a", "3ab", "5a"], context=220))
st.write(RUN / "raw_model/extraction_inventory.json", {
    "document_id": DOC,
    "pages_attempted": [1, 2, 3, 4, 5, 6],
    "reaction_count": len(R),
    "molecule_count": len(M),
    "record_groups": {
        "literature_examples": 3,
        "generic_schemes": 3,
        "optimization_entries": 12,
        "beta_scope_entries": 21,
        "alpha_scope_entries": 6,
        "scheme4_experiments": 5,
        "proposed_mechanisms": 1,
    },
    "machine_status": "needs_review",
    "semantic_accuracy_verified": False,
})
print(json.dumps({"reactions": len(R), "molecules": len(M), "stereo_checks": len(checks)}, ensure_ascii=False))
