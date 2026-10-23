#!/usr/bin/env python3
"""Build the source-derived v3.1 extraction for original PDF pages 1--6.

This is intentionally paper-specific.  It records machine candidates and
source facts without promoting them to source-reviewed chemistry.
"""

import copy
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_fdf0107b6851842e"
RUN = ROOT / "runs/current"
DATASET = json.loads((RUN / "dataset.json").read_text())
E = {p: f"{DOC}:page{p}:text" for p in range(1, 7)}
I = {p: f"{DOC}:page{p}:image" for p in range(1, 7)}


def f(raw=None, normalized=None, evidence=(), method="original_extraction",
      status=None, derivation=None, alternatives=()):
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   derivation=derivation, alternatives=alternatives)


def derived(normalized, evidence, why):
    return f(None, normalized, evidence, method="structure_derivation",
             derivation=why)


def inferred(normalized, evidence, why):
    return f(None, normalized, evidence, method="inference", derivation=why)


def set_fact(obj, key, raw, normalized, evidence, **kw):
    obj[key] = f(raw, normalized, evidence, **kw)


bundle = st.empty(DATASET["documents"])
molecules = {}
stereo_diagnostics = {}


def add_molecule(mid, label, evidence, name=None, smiles=None,
                 chirality="unresolved", geometry="unresolved",
                 representation=None, stereo_note=None, identity_note=None):
    if mid in molecules:
        mol = molecules[mid]
        if not any(o["compound_label"]["normalized_value"] == label for o in mol["occurrences"]):
            mol["occurrences"].append({
                "document_id": DOC,
                "compound_label": f(label, label, evidence),
                "evidence_ids": list(evidence),
            })
        return mol
    mol = st.new_molecule(mid, DOC, evidence)
    mol["identity_status"] = "source_identified"
    mol["occurrences"][0]["compound_label"] = f(label, label, evidence)
    original = name or label
    mol["names"]["original"] = f(original, original, evidence)
    mol["names"]["english"] = f(original, original, evidence)
    mol["names"]["abbreviations"] = f()
    mol["names"]["aliases"] = f()
    mol["cas_numbers"] = f()
    if smiles:
        mol["representation_type"] = "smiles"
        mol["smiles"] = derived(smiles, evidence,
            "Candidate connectivity transcribed from the source drawing/name; parsing is not source verification.")
        diag = st.stereo_inspect(smiles)
        stereo_diagnostics[mid] = diag
        rep = representation or ("encoded" if any(p["specified"] == "Specified" for p in diag["potential_stereo"]) else "not_encoded")
        mol["stereochemistry"]["isomeric_smiles"] = derived(
            diag["canonical_isomeric_smiles"], evidence,
            "RDKit canonical isomeric SMILES from the candidate graph; source_verified=false.")
        mol["stereochemistry"]["identity_key"] = derived(
            diag["stereo_graph_key"], evidence,
            "RDKit stereo graph key from candidate connectivity; not a sample identity.")
    else:
        mol["representation_type"] = "unresolved"
        mol["smiles"] = f()
        rep = representation or "unresolved"
    mol["stereochemistry"]["atom_centers"] = f()
    mol["stereochemistry"]["double_bonds"] = f()
    mol["stereochemistry"]["axial"] = f()
    mol["stereochemistry"]["mixture"] = f()
    mol["stereochemistry"]["assignment"] = f()
    if stereo_note:
        for key, value in stereo_note.items():
            mol["stereochemistry"][key] = f(value, value, evidence)
    mol["stereochemistry"]["chirality_status"] = inferred(
        chirality, evidence, "Status assigned from the visible candidate graph and explicit source stereochemical labels; no unreported configuration is supplied.")
    mol["stereochemistry"]["geometry_status"] = inferred(
        geometry, evidence, "Status assigned from the visible candidate graph and explicit source stereochemical labels.")
    mol["stereochemistry"]["representation_status"] = inferred(
        rep, evidence, "Whether the candidate SMILES retains the source-visible stereochemical information.")
    mol["identity_note"] = f(identity_note, identity_note, evidence) if identity_note else f()
    molecules[mid] = mol
    return mol


# Common materials.  SMILES are identity candidates used to avoid treating
# named simple materials as blank structures; they remain unreviewed after merge.
common = {
    "ag_i": ("AgI", "[Ag]I"),
    "ag_sbf6": ("AgSbF6", "[Ag+].[Sb-](F)(F)(F)(F)(F)F"),
    "ag_no3": ("AgNO3", "[Ag+].[O-][N+](=O)[O-]"),
    "ag_oac": ("AgOAc", "CC(=O)[O-].[Ag+]"),
    "ag_bf4": ("AgBF4", "[Ag+].[B-](F)(F)(F)F"),
    "ag_tfa": ("CF3CO2Ag", "O=C([O-])C(F)(F)F.[Ag+]"),
    "pd_pph3_4": ("Pd(PPh3)4", None),
    "rh2_oac4": ("Rh2(OAc)4", None),
    "bf3_oet2": ("BF3.OEt2", "CCOCC.FB(F)F"),
    "ticl4": ("TiCl4", "Cl[Ti](Cl)(Cl)Cl"),
    "et3n": ("Et3N", "CCN(CC)CC"),
    "lewis_acid": ("Lewis acid (identity not specified in the general procedure)", None),
    "nh4cl": ("saturated aqueous NH4Cl", "[NH4+].[Cl-]"),
    "nahco3": ("saturated aqueous NaHCO3", "[Na+].OC([O-])=O"),
    "nacl": ("saturated aqueous NaCl", "[Na+].[Cl-]"),
    "mgso4": ("anhydrous MgSO4", "[Mg+2].[O-]S(=O)(=O)[O-]"),
    "ch3cn": ("CH3CN", "CC#N"),
    "dcm": ("DCM", "ClCCl"),
    "dce": ("DCE", "ClCCCl"),
    "thf": ("THF", "C1CCOC1"),
    "toluene": ("toluene", "Cc1ccccc1"),
    "dioxane": ("1,4-dioxane", "O1CCOCC1"),
    "ea": ("ethyl acetate", "CCOC(C)=O"),
    "pe": ("petroleum ether", None),
    "silica": ("silica gel", None),
    "nitrogen": ("nitrogen", "N#N"),
}
PAGE5_MATERIALS = {"ticl4", "et3n", "lewis_acid", "nh4cl", "nahco3", "nacl",
                   "mgso4", "ea", "pe", "silica"}
MATERIAL_EVIDENCE_PAGE = {"rh2_oac4": 3, "bf3_oet2": 4}
for mid, (name, smi) in common.items():
    material_page = 5 if mid in PAGE5_MATERIALS else MATERIAL_EVIDENCE_PAGE.get(mid, 2)
    material_evidence = [E[material_page]]
    add_molecule(f"{DOC}:mol:{mid}", name, material_evidence, name=name, smiles=smi,
                 chirality=("unresolved" if mid == "ag_sbf6" else "achiral") if smi else "unresolved",
                 geometry="not_applicable" if smi else "unresolved",
                 representation=("partial" if mid == "ag_sbf6" else "encoded") if smi else "unresolved",
                 identity_note="Simple-material structure candidate; source naming is retained."
                 if smi else "Named material is not expanded to a unique structure candidate.")

add_molecule(f"{DOC}:mol:beta_keto_alpha_diazo", "beta-keto alpha-diazo carbonyl compounds",
             [E[5]], name="beta-keto alpha-diazo carbonyl compounds",
             identity_note="Generic starting-material class; exact member depends on compound 1 instance.")
add_molecule(f"{DOC}:mol:alkynyl_enone", "alpha-alkynyl alpha,beta-unsaturated ketones",
             [E[5]], name="alpha-alkynyl alpha,beta-unsaturated ketones",
             identity_note="Generic starting-material class; exact member depends on compound 1 instance.")


R1_BY_CODE = {
    "aa": ("Ph", "c1ccccc1"), "ba": ("p-MeC6H4", "c1ccc(C)cc1"),
    "ca": ("p-MeOC6H4", "c1ccc(OC)cc1"), "da": ("p-BrC6H4", "c1ccc(Br)cc1"),
    "ea": ("p-ClC6H4", "c1ccc(Cl)cc1"), "fa": ("p-FC6H4", "c1ccc(F)cc1"),
    "ga": ("3-thienyl", "c1ccsc1"), "ha": ("n-Bu", "CCCC"),
    "ia": ("t-Bu", "C(C)(C)C"), "ja": ("H", None),
    "ab": ("Ph", "c1ccccc1"), "bb": ("p-MeC6H4", "c1ccc(C)cc1"),
    "cb": ("p-FC6H4", "c1ccc(F)cc1"), "db": ("3-thienyl", "c1ccsc1"),
    "eb": ("n-Bu", "CCCC"), "fb": ("t-Bu", "C(C)(C)C"), "gb": ("H", None),
}


def substrate_smiles(code):
    _, second = code
    r = R1_BY_CODE[code][1]
    head = "CCOC(=O)" if second == "a" else "O=C(c1ccccc1)"
    tail = "C" if r is None else (f"C{r}" if r.startswith("1") else r)
    # C8 terminal allene carbon is written as C(R1) for substituents and C for H.
    end = "C" if r is None else ("C(" + r + ")")
    return head + "C(=[N+]=[N-])C(=O)CC(c1ccccc1)C(C(C)=O)=C=" + end


def product2_smiles(code):
    _, second = code
    r = R1_BY_CODE[code][1]
    head = "CCOC(=O)" if second == "a" else "O=C(c1ccccc1)"
    subst = "" if r is None else "(" + r + ")"
    return head + "C(=[N+]=[N-])C(=O)CC(c1ccccc1)c1c(C)oc" + subst + "c1"


codes = ["aa", "ba", "ca", "da", "ea", "fa", "ga", "ha", "ia", "ja",
         "ab", "bb", "cb", "db", "eb", "fb", "gb"]
compound_names = {
    "aa": "Ethyl 6-acetyl-2-diazo-3-oxo-5,8-diphenylocta-6,7-dienoate",
    "ba": "Ethyl 6-acetyl-2-diazo-3-oxo-5-phenyl-8-(p-tolyl)octa-6,7-dienoate",
    "ca": "Ethyl 6-acetyl-2-diazo-8-(4-methoxyphenyl)-3-oxo-5-phenylocta-6,7-dienoate",
    "da": "Ethyl 6-acetyl-8-(4-bromophenyl)-2-diazo-3-oxo-5-phenylocta-6,7-dienoate",
    "ea": "Ethyl 6-acetyl-8-(4-chlorophenyl)-2-diazo-3-oxo-5-phenylocta-6,7-dienoate",
    "fa": "Ethyl 6-acetyl-2-diazo-8-(4-fluorophenyl)-3-oxo-5-phenylocta-6,7-dienoate",
    "ga": "Ethyl 6-acetyl-2-diazo-3-oxo-5-phenyl-8-(thiophen-3-yl)octa-6,7-dienoate",
    "ha": "Ethyl 6-acetyl-2-diazo-3-oxo-5-phenyldodeca-6,7-dienoate",
    "ia": "Ethyl 6-acetyl-2-diazo-9,9-dimethyl-3-oxo-5-phenyldeca-6,7-dienoate",
    "ja": "Ethyl 6-acetyl-2-diazo-3-oxo-5-phenylocta-6,7-dienoate",
    "ab": "2-Diazo-1,5-diphenyl-6-(2-phenylvinylidene)octane-1,3,7-trione",
    "bb": "2-Diazo-1,5-diphenyl-6-(2-(p-tolyl)vinylidene)octane-1,3,7-trione",
    "cb": "2-Diazo-6-(2-(4-fluorophenyl)vinylidene)-1,5-diphenyloctane-1,3,7-trione",
    "db": "2-Diazo-1,5-diphenyl-6-(2-(thiophen-3-yl)vinylidene)octane-1,3,7-trione",
    "eb": "2-Diazo-6-(hex-1-en-1-ylidene)-1,5-diphenyloctane-1,3,7-trione",
    "fb": "2-Diazo-6-(3,3-dimethylbut-1-en-1-ylidene)-1,5-diphenyloctane-1,3,7-trione",
    "gb": "2-Diazo-1,5-diphenyl-6-vinylideneoctane-1,3,7-trione",
}
for code in codes:
    page = 5 if code in codes[:10] else 6
    axial = [] if code in ("ja", "gb") else [{"axis": "allene C=C=C", "configuration": "not specified in source"}]
    identity_note = "Candidate graph from the printed systematic name; absolute/axial configuration is not reported."
    molecule_evidence = [E[2], E[page]]
    if code == "aa":
        identity_note += " Page 1 prose calls 1aa '6-benzoyl', whereas Table 1 and section 4.2.1 call it '6-acetyl'; the latter drawing/name is used only as an unreviewed candidate and the source conflict is retained."
        molecule_evidence = [E[1], E[2], E[5]]
    add_molecule(f"{DOC}:mol:1{code}", f"1{code}", molecule_evidence,
                 name=compound_names[code], smiles=substrate_smiles(code),
                 chirality="unspecified", geometry="not_applicable",
                 representation="not_encoded",
                 stereo_note={"axial": axial} if axial else None,
                 identity_note=identity_note)
    add_molecule(f"{DOC}:mol:2{code}", f"2{code}", [E[2]],
                 name=f"diazo trisubstituted furan 2{code}", smiles=product2_smiles(code),
                 chirality="unspecified", geometry="not_applicable",
                 representation="not_encoded",
                 identity_note="Candidate graph follows the Table 2 substituent drawing; the benzylic stereocenter is unspecified.")


product3_codes = ["aa", "fa", "ha", "ia", "ja", "ab", "bb", "cb", "eb", "gb"]
product4_codes = ["aa", "ba", "ca", "da", "ea", "fa", "ga", "ha", "ia", "ja",
                  "ab", "bb", "cb", "eb", "gb"]
for code in product3_codes:
    add_molecule(f"{DOC}:mol:3{code}", f"3{code}", [E[3], E[4]],
                 name=f"nonisomerization product 3{code}", chirality="unspecified",
                 geometry="specified", representation="not_encoded",
                 stereo_note={"double_bonds": [{"configuration": "Z", "basis": "product label/drawing"}]},
                 identity_note="Specific connectivity is visible, but a concrete SMILES was withheld to avoid overclaiming ring/exocyclic geometry transcription.")


def product4_smiles(code):
    _, second = code
    r_label = R1_BY_CODE[code][0]
    side = {
        "Ph": "CC(=O)c2ccccc2", "p-MeC6H4": "CC(=O)c2ccc(C)cc2",
        "p-MeOC6H4": "CC(=O)c2ccc(OC)cc2", "p-BrC6H4": "CC(=O)c2ccc(Br)cc2",
        "p-ClC6H4": "CC(=O)c2ccc(Cl)cc2", "p-FC6H4": "CC(=O)c2ccc(F)cc2",
        "3-thienyl": "CC(=O)c2ccsc2", "n-Bu": "CC(=O)CCCC",
        "t-Bu": "CC(=O)C(C)(C)C", "H": "CC=O",
    }[r_label]
    acyl = "C(=O)OCC" if second == "a" else "C(=O)c2ccccc2"
    return f"Cc1c({side})c(-c2ccccc2)cc(O)c1{acyl}"


for code in product4_codes:
    add_molecule(f"{DOC}:mol:4{code}", f"4{code}", [E[3], E[4]],
                 name=f"penta-substituted aromatic compound 4{code}",
                 smiles=product4_smiles(code), chirality="achiral",
                 geometry="not_applicable", representation="encoded",
                 identity_note="Candidate connectivity derived from the Scheme 2/Table 3 drawing and substituent table; no stereoisomer is asserted.")


def participant(pid, mid, role, evidence, stage="reaction", amount=None,
                mass=None, volume=None, equiv=None, molpct=None, ref=None, product_role=None):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    if amount is not None: p["charges"]["amount"] = f(str(amount), amount, evidence)
    if mass is not None: p["charges"]["mass"] = f(str(mass), mass, evidence)
    if volume is not None: p["charges"]["volume"] = f(str(volume), volume, evidence)
    if equiv is not None: p["charges"]["equivalents"] = f(str(equiv), equiv, evidence)
    if molpct is not None: p["charges"]["mol_percent"] = f(str(molpct), molpct, evidence)
    if ref is not None: p["equivalence_reference"] = f(ref, ref, evidence)
    if product_role is not None: p["product_role"] = f(product_role, product_role, evidence)
    return p


def base_reaction(rid, evidence, record_type="specific_experiment", categories=None,
                  location=None, name=None, transformation=None):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, evidence) if categories is not None else f()
    r["location"] = f(location, location, evidence) if location else f()
    r["reaction_name"] = f(name, name, evidence) if name else f()
    r["transformation_type"] = inferred(transformation, evidence,
        "Transformation description inferred from the source reaction arrow and discussion.") if transformation else f()
    r["review"]["unresolved"] = []
    return r


def add_step(r, sid, evidence, pids=(), operation=None, temperature=None, time=None,
             atmosphere=None, stage="reaction"):
    s = st.new_step(sid, evidence, stage=stage)
    s["participant_ids"] = list(pids)
    if operation is not None: s["operation"] = f(operation, operation, evidence)
    if temperature is not None: s["temperature"] = f(str(temperature), temperature, evidence)
    if time is not None: s["time"] = f(str(time), time, evidence)
    if atmosphere is not None: s["atmosphere"] = f(atmosphere, atmosphere, evidence)
    r["steps"].append(s)
    return s


def add_outcome(r, oid, kind, targets, evidence, value, raw=None, basis=None, role="main product", scope=None):
    o = st.new_outcome(oid, kind, targets, evidence)
    o["value"] = f(raw if raw is not None else str(value), value, evidence)
    o["product_role"] = f(role, role, evidence)
    if kind == "yield":
        o["yield_type"] = f("isolated yield", "isolated", evidence)
    if basis is not None:
        o["measurement_basis"] = f(basis, basis, evidence)
    if scope is not None:
        o["scope"] = f(scope, scope, evidence)
    r["outcomes"].append(o)
    return o


def set_classification(r, evidence, changes=None, labels=None, coverage="unresolved", limitations=()):
    r["classification"]["mapped_reaction_smiles"] = f()
    r["classification"]["bond_changes"] = [derived(ch, evidence, "Net bond change inferred from source reactant/product drawings; no complete atom mapping asserted.") for ch in (changes or [])]
    r["classification"]["labels"] = derived(labels or [], evidence,
        "Taxonomy labels follow the listed net formed/promoted bonds only.") if coverage != "unresolved" else f()
    r["classification"]["coverage"] = inferred(coverage, evidence,
        "Classification coverage is limited by unverified candidate structures and absent complete atom mapping.")
    r["classification"]["limitations"] = list(limitations)
    r["classification"]["model_votes"] = []


CHANGE_CO = {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}], "before_order": 0, "after_order": 1, "change_type": "formed"}
LABEL_CO = [{"category_number": 3, "code": "C-O", "supporting_change_ids": ["b1"]}]
CHANGE_CC = {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"}
LABEL_CC = [{"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1"]}]
CHANGE_OH = {"change_id": "b1", "atoms": [{"atom_map": None, "element": "O"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1}
CHANGE_DBCC = {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 1, "after_order": 2, "change_type": "order_increased"}
LABEL_ISO = [
    {"category_number": 24, "code": "O-H", "supporting_change_ids": ["b1"]},
    {"category_number": 13, "code": "DB-C-C", "supporting_change_ids": ["b2"]},
]


reactions = []

# Scheme 1: three cited prior-work arrows and two this-work generic arrows.
r = base_reaction("literature_scheme1a_cyclopropanation", [E[2], I[2]], "literature_example",
                  location={"page": 2, "region": "Scheme 1a"},
                  name="Cyclopropanation reaction (refs 7-9)", transformation="carbene cyclopropanation")
set_classification(r, [E[2]], [CHANGE_CC], LABEL_CC, "partial", ["Generic R groups; cited prior work, not a concrete experiment of this paper."])
reactions.append(r)
r = base_reaction("literature_scheme1b_furan_formation", [E[2], I[2]], "literature_example",
                  location={"page": 2, "region": "Scheme 1b first arrow"},
                  name="C(sp2)-H insertion sequence: furan formation (ref 15)", transformation="Pd-catalyzed furan formation")
set_classification(r, [E[2]], [CHANGE_CO], LABEL_CO, "partial", ["Generic prior-work arrow; exact atom mapping not supplied."])
reactions.append(r)
r = base_reaction("literature_scheme1b_ch_insertion", [E[2], I[2]], "literature_example",
                  location={"page": 2, "region": "Scheme 1b second arrow"},
                  name="C(sp2)-H insertion of carbenoid (ref 15)", transformation="Rh(II)-mediated C-H insertion")
set_classification(r, [E[2]], [CHANGE_CC], LABEL_CC, "partial", ["Generic prior-work arrow; C-H cleavage is not mislabeled as formed C-H."])
reactions.append(r)

r = base_reaction("general_furanization", [E[1], E[2], E[6], I[1], I[2], I[6]], "generic_scheme",
                  categories=["other"], location={"pages": [1, 2, 6], "regions": ["graphical abstract", "Scheme 1c first arrow", "section 4.3"]},
                  name="Synthesis of diazo trisubstituted furans 2", transformation="Ag(I)-catalyzed cycloisomerization of diazo cumulated allenyl ketones")
r["participants"] = [
    participant("generic_1", f"{DOC}:mol:1aa", "reactant", [E[2]]),
    participant("generic_2", f"{DOC}:mol:2aa", "product", [E[2]], product_role="representative product"),
    participant("ag_tfa", f"{DOC}:mol:ag_tfa", "catalyst", [E[2]], molpct={"number": 20, "unit": "mol%"}, ref=f"{DOC}:mol:1aa"),
    participant("dce", f"{DOC}:mol:dce", "solvent", [E[2]]),
]
set_classification(r, [E[2]], [CHANGE_CO], LABEL_CO, "partial", ["Generic scheme and representative structures only."])
reactions.append(r)

r = base_reaction("general_rh_cascade", [E[1], E[2], E[3], I[1], I[2], I[3]], "generic_scheme",
                  categories=["other"], location={"pages": [1, 2, 3], "regions": ["graphical abstract", "Scheme 1c second arrow", "Table 3 header"]},
                  name="Rhodium(II)-mediated cascade cyclopropanation/rearrangement/isomerization", transformation="intramolecular cyclopropanation, ring opening and optional aromatization")
r["participants"] = [participant("generic_2", f"{DOC}:mol:2aa", "reactant", [E[3]]),
                     participant("generic_3", f"{DOC}:mol:3aa", "product", [E[3]], product_role="nonisomerization product"),
                     participant("generic_4", f"{DOC}:mol:4aa", "product", [E[3]], product_role="aromatized product"),
                     participant("rh", f"{DOC}:mol:rh2_oac4", "catalyst", [E[3]], molpct={"number": 10, "unit": "mol%"}, ref=f"{DOC}:mol:2aa"),
                     participant("dcm", f"{DOC}:mol:dcm", "solvent", [E[3]])]
set_classification(r, [E[3]], [CHANGE_CC], LABEL_CC, "partial", ["Product-dependent tautomerization changes are not generalized across all entries."])
reactions.append(r)


# Table 1 optimization rows.
table1 = [
    (1, "ag_i", "ch3cn", 90, 18, None, "NR"), (2, "ag_sbf6", "ch3cn", 90, 20, None, "NR"),
    (3, "ag_no3", "ch3cn", 90, 7, 40, "no light; substrate not fully converted"),
    (4, "ag_oac", "ch3cn", 90, 15, 30, "substrate not fully converted"),
    (5, "ag_bf4", "ch3cn", 90, 20, 40, "substrate not fully converted"),
    (6, "ag_tfa", "ch3cn", 90, 16, 42, "substrate not fully converted"),
    (7, "pd_pph3_4", "thf", 25, 66, 30, "under nitrogen; by-product formed; substrate not fully converted"),
    (8, "ag_tfa", "dcm", 25, 8.5, 86, "substrate fully consumed"),
    (9, "ag_tfa", "dce", 25, 8.5, 88, "substrate fully consumed"),
    (10, "ag_tfa", "thf", 25, 35, 24, "substrate not fully converted; by-product formed"),
    (11, "ag_tfa", "toluene", 25, 31, 58, "substrate not fully converted"),
    (12, "ag_tfa", "dioxane", 25, 30, 61, "substrate not fully converted"),
    (13, "ag_tfa", "dce", 40, 2.5, 80, "substrate fully consumed"),
    (14, "ag_tfa", "dce", 60, 0.5, 82, "substrate fully consumed; prolonged time may transform 2aa to unknown side products"),
    (15, "ag_tfa", "dce", 80, 0.4, 80, "substrate fully consumed; prolonged time may transform 2aa to unknown side products"),
    (16, "ag_tfa", "dce", 60, 5, 75, "substrate fully consumed"),
    (17, "ag_tfa", "dce", 60, 2, 80, "substrate fully consumed"),
    (18, "ag_tfa", "dce", 60, 1.5, 80, "substrate fully consumed"),
]
loadings = {"ag_i": 10, "ag_sbf6": 10, "ag_no3": 10, "ag_oac": 10, "ag_bf4": 10,
            "pd_pph3_4": 10, "ag_tfa": 20}
for entry, cat, solv, temp, hours, yld, note in table1:
    loading = {16: 5, 17: 10, 18: 15}.get(entry, loadings[cat])
    ev = [E[2], I[2]]
    r = base_reaction(f"table1_entry_{entry:02d}", ev, categories=["optimization"],
                      location={"page": 2, "region": f"Table 1 entry {entry}"},
                      name="Reaction condition optimization for synthesis of 2aa", transformation="Ag/Pd-catalyzed furan formation")
    r["general_scheme_ids"] = [f"{DOC}:general_furanization"]
    r["participants"] = [
        participant("substrate", f"{DOC}:mol:1aa", "reactant", ev, amount={"number": 0.2, "unit": "mmol"}, mass={"number": 80, "unit": "mg"}),
        participant("product", f"{DOC}:mol:2aa", "product", ev, product_role="desired product"),
        participant("catalyst", f"{DOC}:mol:{cat}", "catalyst", ev, molpct={"number": loading, "unit": "mol%"}, ref=f"{DOC}:mol:1aa"),
        participant("solvent", f"{DOC}:mol:{solv}", "solvent", ev, volume={"number": 3, "unit": "mL"}),
    ]
    r["reaction_scale"] = f("1aa (0.2 mmol, 80 mg)", {"reference_molecule_id": f"{DOC}:mol:1aa", "amount": {"number": 0.2, "unit": "mmol"}}, ev)
    r["solvent_conditions"]["composition"] = f(common[solv][0], [common[solv][0]], ev)
    r["solvent_conditions"]["total_volume"] = f("3 mL", {"number": 3, "unit": "mL"}, ev)
    add_step(r, "reaction", ev, ["substrate", "catalyst", "solvent"], "stir under the tabulated conditions", {"number": temp, "unit": "°C"}, {"number": hours, "unit": "h"}, "nitrogen" if entry == 7 else None)
    if yld is None:
        add_outcome(r, "no_reaction", "no_reaction", [], ev, "NR", basis="not reported", role="no desired product")
        set_classification(r, ev, [], [], "not_applicable", ["No reaction observed; expected transformation is not classified."])
    else:
        add_outcome(r, "yield_2aa", "yield", [f"{DOC}:mol:2aa"], ev,
                    {"number": yld, "unit": "%"}, raw=f"{yld}%", basis="isolated by column chromatography")
        set_classification(r, ev, [CHANGE_CO], LABEL_CO, "partial", ["Candidate net cyclization; complete atom mapping absent."])
    r["notes"].append(f(note, note, ev))
    r["notes"].append(f("Source naming conflict for 1aa: page 1 prose says 6-benzoyl, while Table 1 and section 4.2.1 say 6-acetyl.",
                         "1aa source naming conflict retained", [E[1], E[2], E[5]]))
    r["review"]["unresolved"].append("1aa is called 6-benzoyl in page 1 prose but 6-acetyl in Table 1/section 4.2.1; candidate structure follows the latter without resolving the source conflict.")
    if "not fully" in note or "by-product" in note or "prolonged" in note:
        r["review"]["unresolved"].append(note)
    reactions.append(r)


# Table 2 scope rows.  2aa additionally carries the exact section 4.3 charge/procedure evidence.
scope2 = {
    "aa": (0.5, 82), "ba": (0.4, 78), "ca": (0.5, 80), "da": (0.6, 81), "ea": (0.4, 80),
    "fa": (0.5, 78), "ga": (0.5, 80), "ha": (0.4, 95), "ia": (0.3, 80), "ja": (1.0, 65),
    "ab": (0.5, 80), "bb": (0.5, 80), "cb": (0.5, 75), "db": (0.5, 79), "eb": (0.4, 82),
    "fb": (0.6, 85), "gb": (0.5, 79),
}
for idx, code in enumerate(codes, 1):
    hours, yld = scope2[code]
    ev = [E[2], I[2]] + ([E[6], I[6]] if code == "aa" else [])
    r = base_reaction(f"table2_scope_2{code}", ev, categories=["substrate_scope"],
                      location={"page": 2, "region": f"Table 2 entry {idx}"},
                      name="Formation of diazo trisubstituted furan", transformation="Ag(I)-catalyzed furan formation")
    r["general_scheme_ids"] = [f"{DOC}:general_furanization"]
    sub_amount = {"number": 0.2, "unit": "mmol"}
    r["participants"] = [
        participant("substrate", f"{DOC}:mol:1{code}", "reactant", ev, amount=sub_amount),
        participant("product", f"{DOC}:mol:2{code}", "product", ev, product_role="main product"),
        participant("catalyst", f"{DOC}:mol:ag_tfa", "catalyst", ev, molpct={"number": 20, "unit": "mol%"}, ref=f"{DOC}:mol:1{code}"),
        participant("solvent", f"{DOC}:mol:dce", "solvent", ev, volume={"number": 3, "unit": "mL"}),
    ]
    r["reaction_scale"] = f("substrate 1 (0.2 mmol)", {"reference_molecule_id": f"{DOC}:mol:1{code}", "amount": sub_amount}, [E[2]])
    r["solvent_conditions"]["composition"] = f("DCE", ["DCE"], ev)
    r["solvent_conditions"]["total_volume"] = f("3 mL (Table 2)", {"number": 3, "unit": "mL"}, [E[2]])
    add_step(r, "reaction", ev, ["substrate", "catalyst", "solvent"], "reaction under Table 2 conditions", {"number": 60, "unit": "°C"}, {"number": hours, "unit": "h"})
    add_outcome(r, f"yield_2{code}", "yield", [f"{DOC}:mol:2{code}"], ev,
                {"number": yld, "unit": "%"}, raw=f"{yld}%", basis="isolated by column chromatography")
    set_classification(r, ev, [CHANGE_CO], LABEL_CO, "partial", ["No complete atom mapping; structure candidates remain unreviewed."])
    if code == "aa":
        # Preserve the printed internal inconsistencies rather than repairing them.
        table_scale = f("Table 2: 1aa (0.2 mmol)", {"reference_molecule_id": f"{DOC}:mol:1aa", "amount": {"number": 0.2, "unit": "mmol"}}, [E[2]])
        procedure_scale = f("Section 4.3: 1aa (0.40 g, 0.1 mmol)", {"reference_molecule_id": f"{DOC}:mol:1aa", "amount": {"number": 0.1, "unit": "mmol"}}, [E[6]])
        r["reaction_scale"] = f(["Table 2: 0.2 mmol", "Section 4.3: 0.40 g, 0.1 mmol"], None,
                                  [E[2], E[6]], status="conflict",
                                  alternatives=[table_scale, procedure_scale])
        r["participants"][0]["charges"]["amount"] = f(
            ["0.2 mmol (Table 2)", "0.1 mmol (section 4.3)"], None, [E[2], E[6]],
            status="conflict",
            alternatives=[f("0.2 mmol", {"number": 0.2, "unit": "mmol"}, [E[2]]),
                          f("0.1 mmol", {"number": 0.1, "unit": "mmol"}, [E[6]])])
        r["participants"][0]["charges"]["mass"] = f("0.40 g", {"number": 0.40, "unit": "g"}, [E[6]])
        r["participants"][2]["charges"]["mass"] = f("6.9 mg", {"number": 6.9, "unit": "mg"}, [E[6]])
        r["participants"][2]["charges"]["amount"] = f("0.02 mmol", {"number": 0.02, "unit": "mmol"}, [E[6]])
        r["participants"].extend([
            participant("dce_substrate_solution", f"{DOC}:mol:dce", "solvent", [E[6]], volume={"number": 3, "unit": "mL"}),
            participant("dce_catalyst_solution", f"{DOC}:mol:dce", "solvent", [E[6]], volume={"number": 5, "unit": "mL"}),
        ])
        r["procedure"]["original_text"] = f("Section 4.3 general procedure: '0.40 g, 0.1 mmol' 1aa in DCE (3 mL) added to CF3CO2Ag (6.9 mg, 0.02 mmol) in DCE (5 mL), 60 °C, 30 min; product 2aa 86.7 mg, 82%.",
                                                   "Section 4.3 procedure retained verbatim in substance", [E[6]])
        r["procedure"]["addition_order"] = f("1aa solution in DCE was syringed dropwise into vigorously stirred CF3CO2Ag solution in DCE at 60 °C", "substrate solution into catalyst solution", [E[6]])
        r["notes"].append(f("Conflict/internal inconsistency: Table 2 uses 0.2 mmol substrate and 3 mL DCE, while section 4.3 prints 0.40 g, 0.1 mmol in 3 mL plus catalyst in 5 mL; 0.40 g is also inconsistent with 1aa molecular mass and the reported 86.7 mg product.",
                            "Printed charge/volume inconsistency preserved; no correction selected.", [E[2], E[6]]))
        r["review"]["unresolved"].append("2aa scale and DCE volume conflict between Table 2 and section 4.3; printed 0.40 g/0.1 mmol is internally inconsistent.")
        add_step(r, "concentrate", [E[6]], ["substrate", "catalyst", "solvent"], "concentrate under reduced pressure", stage="workup")
        add_step(r, "purification", [E[6]], ["product"], "silica-gel column chromatography with a suitable PE/EA mixture", stage="purification")
        add_outcome(r, "product_amount_2aa_general_procedure", "product_amount", [f"{DOC}:mol:2aa"], [E[6]], {"number": 86.7, "unit": "mg"}, raw="86.7 mg", basis="isolated after column chromatography; section 4.3 general-procedure paragraph")
        add_outcome(r, "product_amount_2aa_characterization_paragraph", "product_amount", [f"{DOC}:mol:2aa"], [E[6]], {"number": 66.0, "unit": "mg"}, raw="66.0 mg", basis="isolated material; section 4.3.1 product paragraph")
        r["notes"].append(f("Source conflict on isolated mass of 2aa: section 4.3 general procedure gives 86.7 mg, while section 4.3.1 gives 66.0 mg; both state 82%.",
                            "Both reported product masses retained as separate targeted outcomes.", [E[6]]))
        r["review"]["unresolved"].append("The page 6 general procedure reports 86.7 mg 2aa, but section 4.3.1 reports 66.0 mg; both report 82%, and neither value is selected.")
    reactions.append(r)


# Table 3 Rh-mediated scope.
table3 = {
    "aa": ("3", 77), "ba": ("4", 68), "ca": ("4", 82), "da": ("4", 68), "ea": ("4", 73),
    "fa": ("3", 82), "ga": ("4", 80), "ha": ("3", 84), "ia": ("3", 60), "ja": ("3", 82),
    "ab": ("3", 42), "bb": ("3", 65), "cb": ("3", 70), "db": (None, 0), "eb": ("3", 75),
    "fb": (None, 0), "gb": ("3", 40),
}
for idx, code in enumerate(codes, 1):
    series, yld = table3[code]
    ev = [E[3], I[3]]
    r = base_reaction(f"table3_scope_2{code}", ev, categories=["substrate_scope"],
                      location={"page": 3, "region": f"Table 3 entry {idx}"},
                      name="Rh2(OAc)4-mediated cyclization of diazo trisubstituted furan", transformation="cascade intramolecular cyclopropanation and rearrangement")
    r["general_scheme_ids"] = [f"{DOC}:general_rh_cascade"]
    r["participants"] = [
        participant("substrate", f"{DOC}:mol:2{code}", "reactant", ev, amount={"number": 0.2, "unit": "mmol"}),
        participant("catalyst", f"{DOC}:mol:rh2_oac4", "catalyst", ev, molpct={"number": 10, "unit": "mol%"}, ref=f"{DOC}:mol:2{code}"),
        participant("solvent", f"{DOC}:mol:dcm", "solvent", ev, volume={"number": 3, "unit": "mL"}),
    ]
    if series:
        product_mid = f"{DOC}:mol:{series}{code}"
        r["participants"].append(participant("product", product_mid, "product", ev,
                                              product_role="nonisomerization product" if series == "3" else "aromatized product"))
    r["reaction_scale"] = f("furan derivative 2 (0.2 mmol)", {"reference_molecule_id": f"{DOC}:mol:2{code}", "amount": {"number": 0.2, "unit": "mmol"}}, ev)
    r["solvent_conditions"]["composition"] = f("DCM", ["DCM"], ev)
    r["solvent_conditions"]["total_volume"] = f("3 mL", {"number": 3, "unit": "mL"}, ev)
    add_step(r, "reaction", ev, [p["participant_id"] for p in r["participants"] if p["role"] != "product"], "stir under nitrogen", {"number": 25, "unit": "°C"}, {"number": 30, "unit": "min"}, "nitrogen")
    if series:
        add_outcome(r, "yield", "yield", [product_mid], ev, {"number": yld, "unit": "%"}, raw=f"{yld}%", basis="isolated by column chromatography")
        if series == "3":
            add_outcome(r, "ratio_3_4", "selectivity", [product_mid], ev, "100:0", raw="3/4 = 100:0", basis="not reported", scope={"metric": "product ratio", "order": [f"3{code}", f"4{code}"]})
            set_classification(r, ev, [CHANGE_CC], LABEL_CC, "partial", ["Net C-C ring-forming change inferred; N2 loss and other bond-order changes are not fully mapped."])
        else:
            add_outcome(r, "ratio_3_4", "selectivity", [product_mid], ev, "0:100", raw="3/4 = 0:100", basis="not reported", scope={"metric": "product ratio", "order": [f"3{code}", f"4{code}"]})
            changes = [copy.deepcopy(CHANGE_CC), copy.deepcopy(CHANGE_OH), copy.deepcopy(CHANGE_DBCC)]
            changes[1]["change_id"] = "b2"; changes[2]["change_id"] = "b3"
            labels = [
                {"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1"]},
                {"category_number": 24, "code": "O-H", "supporting_change_ids": ["b2"]},
                {"category_number": 13, "code": "DB-C-C", "supporting_change_ids": ["b3"]},
            ]
            set_classification(r, ev, changes, labels, "partial", ["Aromatization labels are net graph changes, not elementary mechanistic steps."])
    else:
        add_outcome(r, "desired_product_yield", "other", [], ev, {"number": 0, "unit": "%"}, raw="0%; by-product was produced", basis="not reported", role="no desired 3 or 4 isolated")
        r["notes"].append(f("By-product was produced; its identity and amount are not reported.", "By-product was produced", ev))
        r["review"]["unresolved"].append("Unknown by-product identity and measurement basis; no expected product was invented.")
        set_classification(r, ev, [], [], "unresolved", ["Observed by-product connectivity is unknown, so no bond-formation label is asserted."])
    reactions.append(r)


# Scheme 2: BF3-promoted isomerization, ten separate arrows/products.
scheme2 = {
    "aa": (1.0, 78), "fa": (1.0, 80), "ha": (1.5, 86), "ia": (1.0, 65), "ja": (5.0, 0),
    "ab": (1.5, 50), "bb": (2.5, 0), "cb": (4.0, 70), "eb": (2.5, 78), "gb": (12.0, 0),
}
for code, (hours, yld) in scheme2.items():
    ev = [E[4], I[4]]
    r = base_reaction(f"scheme2_isomerization_3{code}", ev, categories=["derivatization"],
                      location={"page": 4, "region": f"Scheme 2 product 4{code}"},
                      name="Transformation of compound 3 into penta-substituted aromatic compound 4", transformation="BF3-promoted isomerization/aromatization")
    r["participants"] = [
        participant("substrate", f"{DOC}:mol:3{code}", "reactant", ev, amount={"number": 0.2, "unit": "mmol"}),
        participant("product", f"{DOC}:mol:4{code}", "product", ev, product_role="aromatized product"),
        participant("bf3", f"{DOC}:mol:bf3_oet2", "reagent", ev, equiv={"number": 1, "unit": "equiv"}, ref=f"{DOC}:mol:3{code}"),
        participant("dcm", f"{DOC}:mol:dcm", "solvent", ev, volume={"number": 2, "unit": "mL"}),
    ]
    r["reaction_scale"] = f("3 (0.2 mmol)", {"reference_molecule_id": f"{DOC}:mol:3{code}", "amount": {"number": 0.2, "unit": "mmol"}}, ev)
    r["solvent_conditions"]["composition"] = f("DCM", ["DCM"], ev)
    r["solvent_conditions"]["total_volume"] = f("2 mL", {"number": 2, "unit": "mL"}, ev)
    add_step(r, "reaction", ev, ["substrate", "bf3", "dcm"], "stir under Scheme 2 conditions", {"number": 25, "unit": "°C"}, {"number": hours, "unit": "h"})
    if yld:
        add_outcome(r, "yield", "yield", [f"{DOC}:mol:4{code}"], ev, {"number": yld, "unit": "%"}, raw=f"{yld}%", basis="isolated by column chromatography")
        set_classification(r, ev, [CHANGE_OH, CHANGE_DBCC], LABEL_ISO, "partial", ["Hydrogen source/acceptor and full aromatic bond-order mapping are unresolved."])
    else:
        add_outcome(r, "desired_product_yield", "other", [], ev, {"number": 0, "unit": "%"}, raw="0%; by-product was formed", basis="not reported", role="desired product not isolated")
        r["notes"].append(f("By-product was formed; identity and amount unreported.", "By-product was formed", ev))
        r["review"]["unresolved"].append("Unknown by-product identity; 4 product is shown only as the intended target at 0%.")
        set_classification(r, ev, [], [], "unresolved", ["Actual by-product structure is unknown."])
    reactions.append(r)


# Scheme 3 one-pot attempt is one concrete experiment with three explicit outcomes.
ev = [E[3], E[4], I[4]]
r = base_reaction("scheme3_one_pot_attempt", ev, categories=["control"],
                  location={"page": 4, "region": "Scheme 3"},
                  name="Attempt for one-pot synthesis of penta-substituted aromatic compound 4", transformation="one-pot furanization/cascade/isomerization attempt")
r["participants"] = [
    participant("substrate", f"{DOC}:mol:1aa", "reactant", ev),
    participant("p2", f"{DOC}:mol:2aa", "product", ev, product_role="only isolated product"),
    participant("p3", f"{DOC}:mol:3aa", "product", ev, product_role="intended downstream product; not formed"),
    participant("p4", f"{DOC}:mol:4aa", "product", ev, product_role="intended final product; not formed"),
    participant("rh", f"{DOC}:mol:rh2_oac4", "catalyst", ev),
    participant("ag", f"{DOC}:mol:ag_tfa", "catalyst", ev),
    participant("bf3", f"{DOC}:mol:bf3_oet2", "reagent", ev),
    participant("dcm", f"{DOC}:mol:dcm", "solvent", ev),
]
r["solvent_conditions"]["composition"] = f("DCM", ["DCM"], ev)
add_step(r, "one_pot_reaction", ev, [p["participant_id"] for p in r["participants"] if p["role"] != "product"], "one-pot reaction with Rh2(OAc)4/CF3CO2Ag/BF3.OEt2", {"number": 25, "unit": "°C"}, {"number": 10, "unit": "h"}, "nitrogen")
add_outcome(r, "yield_2aa", "yield", [f"{DOC}:mol:2aa"], ev, {"number": 50, "unit": "%"}, raw="2aa, 50%", basis="not reported")
add_outcome(r, "yield_3aa", "yield", [f"{DOC}:mol:3aa"], ev, {"number": 0, "unit": "%"}, raw="3aa, 0%", basis="not reported", role="not formed")
add_outcome(r, "yield_4aa", "yield", [f"{DOC}:mol:4aa"], ev, {"number": 0, "unit": "%"}, raw="4aa, 0%", basis="not reported", role="not formed")
r["notes"].append(f("Authors suggest Rh2(OAc)4 deactivation as a possible reason for failure of downstream conversion.", "perhaps because of deactivation of Rh2(OAc)4", [E[3]]))
r["review"]["unresolved"].append("Exact reagent loadings and isolated-yield basis are not printed in Scheme 3.")
set_classification(r, ev, [CHANGE_CO], LABEL_CO, "partial", ["Only observed 1aa-to-2aa transformation is classified; unobserved downstream transformations are excluded."])
reactions.append(r)


# Scheme 4 is author-proposed mechanism, not a measured experiment.
ev = [E[3], E[4], I[4]]
r = base_reaction("scheme4_proposed_mechanism", ev, "proposed_mechanism",
                  location={"page": 4, "region": "Scheme 4"},
                  name="Reasonable mechanism for construction of penta-substituted aromatic compound 4", transformation="proposed Rh-carbene cyclopropanation/ring opening/aromatization")
r["notes"] = [f("IA-IF and IB-IB are proposed intermediates; attempted isolation of IB was unsuccessful.", "A reasonable mechanism is proposed; isolation of IB was unsuccessful", ev)]
r["review"]["unresolved"].append("Proposed intermediates are not isolated compounds and are not treated as experimental products.")
set_classification(r, ev, [], [], "unresolved", ["Mechanistic arrow sequence is not a concrete net experiment."])
reactions.append(r)


# Section 4.2 general procedure and 17 reported preparations of compounds 1.
prep_data = {
    "aa": (5, 0.78, 43), "ba": (5, 0.87, 50), "ca": (5, 0.90, 52), "da": (5, 1.49, 50),
    "ea": (5, 1.50, 53), "fa": (5, 1.71, 58), "ga": (5, 0.86, 45), "ha": (5, 0.64, 48),
    "ia": (5, 0.66, 43), "ja": (5, 0.63, 45), "ab": (6, 0.87, 40), "bb": (6, 1.88, 42),
    "cb": (6, 0.77, 46), "db": (6, 2.20, 50), "eb": (6, 1.95, 47), "fb": (6, 0.46, 44),
    "gb": (6, 0.68, 42),
}
for code in codes:
    page, mass_g, yld = prep_data[code]
    ev = [E[5], I[5], E[page], I[page]] if page != 5 else [E[5], I[5]]
    r = base_reaction(f"section4_2_prepare_1{code}", ev, categories=["other"],
                      location={"page": page, "region": f"Section 4.2.{codes.index(code)+1}"},
                      name="Preparation of diazo cumulated allenyl ketone", transformation="Michael addition giving diazo cumulated allenyl ketone")
    r["participants"] = [
        participant("diazo_carbonyl", f"{DOC}:mol:beta_keto_alpha_diazo", "reactant", ev, amount={"number": 10, "unit": "mmol"}, equiv={"number": 1, "unit": "equiv"}, ref=f"{DOC}:mol:beta_keto_alpha_diazo"),
        participant("alkynyl_enone", f"{DOC}:mol:alkynyl_enone", "reactant", ev, amount={"number": 10, "unit": "mmol"}, equiv={"number": 1, "unit": "equiv"}, ref=f"{DOC}:mol:beta_keto_alpha_diazo"),
        participant("product", f"{DOC}:mol:1{code}", "product", ev, product_role="isolated diazo cumulated allenyl ketone"),
        participant("ticl4", f"{DOC}:mol:ticl4", "reagent", ev, amount={"number": 10, "unit": "mmol"}, volume={"number": 1.1, "unit": "mL"}, ref=f"{DOC}:mol:beta_keto_alpha_diazo"),
        participant("et3n", f"{DOC}:mol:et3n", "reagent", ev, amount={"number": 10, "unit": "mmol"}, mass={"number": 1, "unit": "g"}, ref=f"{DOC}:mol:beta_keto_alpha_diazo"),
        participant("lewis_acid", f"{DOC}:mol:lewis_acid", "reagent", ev, amount={"number": 10, "unit": "mmol"}, equiv={"number": 1, "unit": "equiv"}, ref=f"{DOC}:mol:beta_keto_alpha_diazo"),
        participant("dcm_main", f"{DOC}:mol:dcm", "solvent", ev, volume={"number": 30, "unit": "mL"}),
        participant("dcm_enone", f"{DOC}:mol:dcm", "solvent", ev, volume={"number": 10, "unit": "mL"}),
        participant("nh4cl", f"{DOC}:mol:nh4cl", "workup_agent", ev, stage="workup", volume={"number": 30, "unit": "mL"}),
        participant("nahco3", f"{DOC}:mol:nahco3", "workup_agent", ev, stage="workup", volume={"number": 30, "unit": "mL"}),
        participant("nacl", f"{DOC}:mol:nacl", "workup_agent", ev, stage="workup", volume={"number": 30, "unit": "mL"}),
        participant("mgso4", f"{DOC}:mol:mgso4", "workup_agent", ev, stage="workup"),
        participant("silica", f"{DOC}:mol:silica", "purification_agent", ev, stage="purification"),
        participant("ea", f"{DOC}:mol:ea", "purification_agent", ev, stage="purification"),
        participant("pe", f"{DOC}:mol:pe", "purification_agent", ev, stage="purification"),
        participant("nitrogen", f"{DOC}:mol:nitrogen", "additive", ev),
    ]
    r["reaction_scale"] = f("beta-keto alpha-diazo carbonyl compound (10 mmol)", {"reference_molecule_id": f"{DOC}:mol:beta_keto_alpha_diazo", "amount": {"number": 10, "unit": "mmol"}}, [E[5]])
    r["solvent_conditions"]["composition"] = f("anhydrous CH2Cl2: 30 mL initial plus 10 mL enone solution", ["anhydrous CH2Cl2"], [E[5]])
    r["solvent_conditions"]["total_volume"] = f("30 mL + 10 mL (separate additions)", None, [E[5]], status="unreadable")
    add_step(r, "low_temperature_addition", [E[5]], ["diazo_carbonyl", "ticl4", "et3n", "dcm_main", "nitrogen"], "add TiCl4 and Et3N dropwise by syringe to vigorously stirred diazo-carbonyl solution", {"number": -78, "unit": "°C"}, None, "nitrogen")
    add_step(r, "low_temperature_hold", [E[5]], ["diazo_carbonyl", "ticl4", "et3n"], "stir dark-red mixture", {"number": -78, "unit": "°C"}, {"number": 1, "unit": "h"}, "nitrogen")
    add_step(r, "second_addition_and_hold", [E[5]], ["lewis_acid", "alkynyl_enone", "dcm_enone"], "syringe Lewis acid and alpha-alkynyl enone solution; stir", None, {"number": 10, "unit": "h"}, "nitrogen")
    add_step(r, "quench", [E[5]], ["nh4cl"], "quench with saturated aqueous NH4Cl", stage="workup")
    add_step(r, "wash_and_dry", [E[5]], ["nahco3", "nacl", "mgso4"], "separate organic layer; wash with saturated aqueous NaHCO3 then NaCl; dry over anhydrous MgSO4; concentrate in vacuum", stage="workup")
    add_step(r, "column", [E[5]], ["silica", "ea", "pe", "product"], "silica-gel column chromatography, EA/PE 1:20", stage="purification")
    r["procedure"]["general_procedure_reference"] = f("Section 4.2 general procedure", "Section 4.2", [E[5]])
    r["procedure"]["addition_order"] = f("TiCl4 and Et3N -> 1 h at -78 °C -> Lewis acid and enone solution -> 10 h -> quench/workup", ["TiCl4 and Et3N", "Lewis acid and enone", "quench"], [E[5]])
    add_outcome(r, "product_mass", "product_amount", [f"{DOC}:mol:1{code}"], [E[page]], {"number": mass_g, "unit": "g"}, raw=f"{mass_g} g", basis="isolated material")
    add_outcome(r, "yield", "yield", [f"{DOC}:mol:1{code}"], [E[page]], {"number": yld, "unit": "%"}, raw=f"{yld}%", basis="isolated after column chromatography")
    r["notes"].append(f("Temperature during the second 10 h hold and exact identity of the separately named 'Lewis acid' are not reported in the printed general procedure.",
                        "Second-stage temperature and Lewis-acid identity unreported", [E[5]]))
    r["review"]["unresolved"].extend(["Temperature for the second 10 h reaction stage is unreported.", "The separately stated Lewis acid is not identified; it is not silently equated with TiCl4.", "Generic reactant identities depend on the specific 1-series member and are not expanded."])
    set_classification(r, ev, [CHANGE_CC], LABEL_CC, "partial", ["Michael-addition C-C formation inferred; complete mapping and proton transfers unresolved."])
    reactions.append(r)


# Source context: prose, hypothetical products, characterization, and general lab information.
bundle["source_context"] = [
    {"context_id": f"{DOC}:context:page1_title_abstract", "document_id": DOC, "evidence_ids": [E[1], I[1]], "kind": "bibliographic_and_summary", "note": "True title/authors/abstract and graphical-abstract duplicate."},
    {"context_id": f"{DOC}:context:figure1_hypothetical_products", "document_id": DOC, "evidence_ids": [E[3], I[3]], "kind": "hypothetical_products", "note": "P-I, P-II and P-III are possible products discussed as not detected, not reported synthesis records."},
    {"context_id": f"{DOC}:context:page4_general_analytics", "document_id": DOC, "evidence_ids": [E[4], I[4], E[5], I[5]], "kind": "general_characterization", "note": "Instrumental methods and general laboratory handling are context, not reaction experiments."},
    {"context_id": f"{DOC}:context:pages5_6_characterization", "document_id": DOC, "evidence_ids": [E[5], I[5], E[6], I[6]], "kind": "characterization", "note": "NMR/IR/HRMS/mp/Rf blocks support compound identity but are not separate syntheses."},
]

bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions

# Raw visual inventory and exact coverage contract.
page_records = {p: [] for p in range(1, 7)}
for r in reactions:
    for eid in r["evidence_ids"]:
        for p in range(1, 7):
            if eid in (E[p], I[p]) and r["reaction_id"] not in page_records[p]:
                page_records[p].append(r["reaction_id"])

regions = {
    1: [
        {"region": "title, abstract and graphical abstract", "finding": "bibliographic context plus duplicate generic two-stage this-work scheme; no concrete row", "evidence_ids": [E[1], I[1]]},
        {"region": "Results and Discussion opening prose", "finding": "model substrate and optimization discussion; concrete rows are on page 2", "evidence_ids": [E[1], I[1]]},
    ],
    2: [
        {"region": "Scheme 1a-c", "finding": "three cited prior-work arrows and two this-work generic arrows", "evidence_ids": [E[2], I[2]]},
        {"region": "Table 1 including footnotes", "finding": "18 independent optimization entries, including two NR controls and footnoted incomplete conversion/by-products", "evidence_ids": [E[2], I[2]]},
        {"region": "Table 2 including footnote", "finding": "17 independent furan-formation scope entries", "evidence_ids": [E[2], I[2]]},
    ],
    3: [
        {"region": "Figure 1", "finding": "hypothetical P-I/P-II/P-III; reported as not detected, retained in source_context", "evidence_ids": [E[3], I[3]]},
        {"region": "Table 3 including footnotes", "finding": "17 independent Rh-mediated scope entries; db/fb give unidentified by-products and 0% desired products", "evidence_ids": [E[3], I[3]]},
        {"region": "discussion beneath Table 3", "finding": "product-series interpretation, 4ea X-ray statement, Scheme 2/3 prose and author explanation", "evidence_ids": [E[3], I[3]]},
    ],
    4: [
        {"region": "Scheme 2 including footnotes", "finding": "10 separate BF3-promoted isomerization arrows; seven isolated products and three 0%/by-product cases", "evidence_ids": [E[4], I[4]]},
        {"region": "Scheme 3", "finding": "one concrete one-pot control yielding only 2aa (50%), with 3aa and 4aa at 0%", "evidence_ids": [E[4], I[4]]},
        {"region": "Scheme 4 and adjacent prose", "finding": "author-proposed mechanism; intermediates not isolated and not concrete experiments", "evidence_ids": [E[4], I[4]]},
        {"region": "Experimental section 4.1 start", "finding": "general analytical/laboratory context only", "evidence_ids": [E[4], I[4]]},
    ],
    5: [
        {"region": "Section 4.2 general procedure", "finding": "ordered -78 °C additions/hold, 10 h second stage, workup and column purification; second-stage temperature and Lewis-acid identity missing", "evidence_ids": [E[5], I[5]]},
        {"region": "Sections 4.2.1-4.2.10", "finding": "10 reported preparations of 1aa-1ja; characterization retained as context", "evidence_ids": [E[5], I[5]]},
    ],
    6: [
        {"region": "Sections 4.2.11-4.2.17", "finding": "seven reported preparations of 1ab-1gb; characterization retained as context", "evidence_ids": [E[6], I[6]]},
        {"region": "Section 4.3 general procedure and 4.3.1", "finding": "2aa scope experiment cross-linked to Table 2; printed scale/volume/mass inconsistencies preserved", "evidence_ids": [E[6], I[6]]},
        {"region": "Sections 4.3.2-4.3.3", "finding": "product characterization only within page range; corresponding reactions already represented by Table 2 rows", "evidence_ids": [E[6], I[6]]},
    ],
}
unresolved = {
    1: ["Graphical abstract is a duplicate summary and does not supply exact conditions."],
    2: ["Table 1 non-entry-7 atmosphere is not explicitly tabulated; it is not silently inherited.", "Scope structures are machine candidates and remain unverified."],
    3: ["By-product identities for Table 3 entries 14 and 16 are unreported.", "Table 3 caption says furan derivatives 5 although the scheme/table uses series 2; source wording retained as a likely label error without correction."],
    4: ["Scheme 2 by-product identities for 4ja/4bb/4gb are unreported.", "Scheme 3 exact reagent loadings and yield basis are not printed."],
    5: ["Second-stage temperature and identity of the separately named Lewis acid are unreported.", "No SI was supplied or associated."],
    6: ["Section 4.3 prints 0.40 g, 0.1 mmol 1aa and 86.7 mg product, internally inconsistent; Table 2 instead uses 0.2 mmol and 3 mL DCE.", "Section 4.3.1 separately reports 66.0 mg 2aa at the same 82% yield; both masses are retained without selection.", "Only original pages 1-6 are claimed visually covered; later product paragraphs are not claimed as covered."],
}
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [{"page": p, "state": "attempted", "regions": regions[p],
               "record_ids": page_records[p], "unresolved": unresolved[p]}
              for p in range(1, 7)],
}
raw = {
    "batch_id": "allene-stereo-v31-200-20261017",
    "document_id": DOC,
    "visual_pages_opened_original_detail": [1, 2, 3, 4, 5, 6],
    "record_count": len(reactions),
    "molecule_count": len(molecules),
    "page_regions": regions,
    "decisions": [
        "Every Table 1/2/3 row is a separate record.",
        "Scheme 2 products are separate experiments; Scheme 3 is one experiment with three outcomes.",
        "Scheme 4 is proposed_mechanism, not experimental proof.",
        "NMR/IR/HRMS/mp/Rf blocks are source_context, not synthesis records.",
        "Table 2 2aa and section 4.3 are merged as the same labeled experiment with both evidence strands and unresolved conflicts.",
        "The page 1 6-benzoyl versus Table 1/section 4.2.1 6-acetyl naming conflict for 1aa is preserved; the candidate graph follows the latter drawing/name without adjudication.",
        "No SI association was inferred.",
    ],
}

st.write(RUN / "submitted.json", bundle)
st.write(RUN / "coverage/pages_0001_0006.json", coverage)
st.write(RUN / "raw_model/pages_0001_0006_raw_extraction.json", raw)
st.write(RUN / "raw_model/pages_0001_0006_standard_stereo_checks.json", stereo_diagnostics)
st.write(RUN / "raw_model/pages_0001_0006_source_mentions.json",
         st.source_mentions(DATASET, ["1aa", "2aa", "3aa", "4aa", "1gb"], context=180))
print(json.dumps({"molecules": len(molecules), "reactions": len(reactions),
                  "pages": list(regions), "stereo_checks": len(stereo_diagnostics)}, ensure_ascii=False))
