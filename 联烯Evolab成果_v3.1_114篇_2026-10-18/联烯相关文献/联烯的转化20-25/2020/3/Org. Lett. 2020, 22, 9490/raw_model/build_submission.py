#!/usr/bin/env python3
"""Build the source-derived v3.1 submission for doc_19fbeef4228775ee.

This is a reproducible transcription of the five visually reviewed journal pages.
It deliberately leaves unavailable SI details and uncertain structures unresolved.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DATASET_PATH = RUN / "dataset.json"
DOC = "doc_19fbeef4228775ee"
P = {n: f"{DOC}:page{n}:image" for n in range(1, 6)}
T = {n: f"{DOC}:page{n}:text" for n in range(1, 6)}


def F(raw=None, norm=None, ev=(), method="original_extraction", derivation=None,
      status=None, verification=None, alternatives=()):
    return st.fact(raw, norm, ev, method=method, derivation=derivation,
                   status=status, verification=verification, alternatives=alternatives)


def missing(ev=(), status="unreported"):
    return F(ev=ev, status=status)


def direct(raw, norm=None, ev=()):
    return F(raw, raw if norm is None else norm, ev)


def inferred(norm, ev, why, raw=None):
    return F(raw, norm, ev, method="inference", derivation=why)


def derived(norm, ev, why, raw=None):
    return F(raw, norm, ev, method="structure_derivation", derivation=why)


def qty(number, unit, raw_text=None, qualifier=None):
    result = {"number": number, "unit": unit}
    if raw_text is not None:
        result["raw_text"] = raw_text
    if qualifier is not None:
        result["qualifier"] = qualifier
    return result


def set_occurrence(m, label, ev):
    m["occurrences"] = [{
        "document_id": DOC,
        "compound_label": direct(label, label, ev),
        "evidence_ids": list(ev),
    }]


def set_stereo_from_smiles(m, smiles, ev, note=None):
    diag = st.stereo_inspect(smiles)
    unspecified_atoms = [x for x in diag["potential_stereo"]
                         if x["type"].startswith("Atom_") and x["specified"] == "Unspecified"]
    unspecified_bonds = [x for x in diag["potential_stereo"]
                         if x["type"].startswith("Bond_") and x["specified"] == "Unspecified"]
    stereo = m["stereochemistry"]
    stereo["atom_centers"] = derived(
        [{"atom_index": x["atom_index"], "assignment": x["CIP"]} for x in diag["atom_centers"]],
        ev, "RDKit graph diagnostic on the drawing-derived candidate SMILES; '?' is unassigned, not racemic.")
    stereo["double_bonds"] = derived(
        diag["double_bonds"], ev,
        "RDKit graph diagnostic on the candidate SMILES; terminal/non-stereogenic double bonds do not imply E/Z.")
    stereo["axial"] = derived([], ev,
        "No source axial assignment is shown. Toolkit output does not certify absence of all possible axial stereochemistry.")
    stereo["mixture"] = missing(ev)
    stereo["assignment"] = direct(note, note, ev) if note else missing(ev)
    stereo["chirality_status"] = derived(
        "unspecified" if unspecified_atoms else "achiral", ev,
        "Candidate graph has unassigned potential tetrahedral stereochemistry." if unspecified_atoms
        else "No potential tetrahedral center was found in the candidate graph; source identity remains needs_review.")
    stereo["geometry_status"] = derived(
        "unspecified" if unspecified_bonds else "not_applicable", ev,
        "Candidate graph has an E/Z-capable bond without an assignment." if unspecified_bonds
        else "No E/Z-capable double bond was found in the candidate graph.")
    stereo["representation_status"] = derived(
        "partial" if unspecified_atoms or unspecified_bonds else "encoded", ev,
        "Potential stereochemistry is not encoded." if unspecified_atoms or unspecified_bonds
        else "All toolkit-detectable stereochemistry in this candidate graph is represented.")
    stereo["isomeric_smiles"] = derived(
        diag["canonical_isomeric_smiles"], ev,
        "Canonical isomeric SMILES from the drawing-derived candidate; unspecified centers remain unspecified.")
    stereo["identity_key"] = derived(
        diag["stereo_graph_key"], ev,
        "Stereo graph key generated from the candidate SMILES; not a sample-identity certificate.")


def molecule(mid, label, ev, name=None, smiles=None, identity_note=None):
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    set_occurrence(m, label, ev)
    m["names"]["original"] = direct(name or label, name or label, ev)
    m["names"]["english"] = direct(name or label, name or label, ev)
    m["names"]["abbreviations"] = missing(ev)
    m["names"]["aliases"] = direct([label], [label], ev)
    m["cas_numbers"] = missing(ev)
    m["external_properties"] = []
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = derived(smiles, ev,
            "Candidate SMILES transcribed from the source drawing; parsing does not establish source correctness.")
        m["generic_template"] = None
        set_stereo_from_smiles(m, smiles, ev, identity_note)
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing(ev)
        m["generic_template"] = None
        for key in ("atom_centers", "double_bonds", "axial", "mixture", "assignment",
                    "isomeric_smiles", "identity_key"):
            m["stereochemistry"][key] = missing(ev)
        for key in ("chirality_status", "geometry_status", "representation_status"):
            m["stereochemistry"][key] = inferred("unresolved", ev,
                "The complete source structure/stereochemical representation was not safely established from the supplied pages.")
    m["identity_note"] = direct(identity_note, identity_note, ev) if identity_note else missing(ev)
    return m


def template_molecule(mid, label, ev, fixed, variables, reason):
    m = molecule(mid, label, ev, name=label)
    m["representation_type"] = "generic_template"
    m["smiles"] = missing(ev)
    m["generic_template"] = {
        "fixed_structure": direct(fixed, fixed, ev),
        "variables": [{
            "label": v,
            "attachment_points": direct(f"source drawing attachment for {v}", f"source drawing attachment for {v}", ev),
            "definition": direct("variable substituent as drawn", "variable substituent as drawn", ev),
            "allowed_values": missing(ev),
            "constraints": missing(ev),
        } for v in variables],
        "instance_mapping": [],
        "reason": direct(reason, reason, ev),
    }
    for key in ("chirality_status", "geometry_status", "representation_status"):
        m["stereochemistry"][key] = inferred("unresolved", ev,
            "Generic R-group drawing cannot establish a unique stereochemical identity.")
    return m


def participant(pid, mid, role, ev, stage="reaction", product_role=None):
    p = st.new_participant(pid, mid, role, ev, stage)
    if product_role:
        p["product_role"] = direct(product_role, product_role, ev)
    return p


def set_charge(p, key, raw, normalized, ev, reference=None, method="original_extraction", why=None):
    p["charges"][key] = F(raw, normalized, ev, method=method, derivation=why)
    if key in ("equivalents", "mol_percent"):
        p["equivalence_reference"] = F(reference, reference, ev,
            method="inference", derivation="The limiting allenol is the stated equivalents/loading denominator.")


def reaction(rid, record_type, ev, categories=None, location=None):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = direct(categories, categories, ev) if categories is not None else missing(ev)
    r["location"] = direct(location or rid, location or rid, ev)
    r["review"] = {
        "status": "needs_review", "method": "model", "reviewer": None,
        "checked_fields": [],
        "unresolved": [
            "No Supporting Information PDF was supplied or explicitly associated; detailed procedure, workup, purification, and characterization remain unavailable.",
            "Machine transcription and structure candidates require expert source review.",
        ],
    }
    return r


def add_step(r, sid, ev, participants, operation=None, temperature=None, time=None,
             atmosphere=None, operation_method="original_extraction"):
    s = st.new_step(sid, ev)
    s["participant_ids"] = list(participants)
    if operation is not None:
        s["operation"] = F(operation, operation, ev, method=operation_method,
                           derivation="Condition inherited from the explicitly described optimal reaction." if operation_method == "inherited_general_condition" else None)
    if temperature is not None:
        s["temperature"] = F(temperature, temperature, ev,
            method="inherited_general_condition" if operation_method == "inherited_general_condition" else "original_extraction",
            derivation="Room temperature inherited from the displayed optimal conditions." if operation_method == "inherited_general_condition" else None)
    if time is not None:
        s["time"] = direct(f"{time} h", qty(time, "h", f"{time} h"), ev)
    if atmosphere is not None:
        s["atmosphere"] = F(atmosphere, atmosphere, ev,
            method="inherited_general_condition" if operation_method == "inherited_general_condition" else "original_extraction",
            derivation="Argon atmosphere inherited from the displayed optimal conditions." if operation_method == "inherited_general_condition" else None)
    r["steps"].append(s)
    return s


def outcome(r, oid, kind, targets, raw, norm, ev, role=None, yield_type=None, basis=None, scope=None):
    o = st.new_outcome(oid, kind, targets, ev)
    o["value"] = direct(raw, norm, ev)
    if role:
        o["product_role"] = direct(role, role, ev)
    if yield_type:
        o["yield_type"] = direct(yield_type, yield_type, ev)
    if basis:
        o["measurement_basis"] = direct(basis, basis, ev)
    if scope:
        o["scope"] = direct(scope, scope, ev)
    r["outcomes"].append(o)
    return o


def unresolved_classification(r, ev, coverage="unresolved", limitation=None):
    r["classification"]["mapped_reaction_smiles"] = missing(ev)
    r["classification"]["bond_changes"] = []
    r["classification"]["labels"] = inferred(None if coverage == "unresolved" else [], ev,
        "No observed product transformation is available for classification." if coverage == "not_applicable"
        else "Exact atom correspondence is not safely established from the supplied generic or complex drawing.")
    r["classification"]["coverage"] = inferred(coverage, ev,
        "No observed product transformation." if coverage == "not_applicable"
        else "Structures or atom correspondence remain unresolved.")
    r["classification"]["limitations"] = [limitation or "No complete atom mapping is asserted."]
    r["classification"]["model_votes"] = []


def cascade_classification(r, ev):
    changes = [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "S"}],
         "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
         "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}],
         "before_order": 1, "after_order": 2, "change_type": "order_increased"},
    ]
    why = ("Partial net-change assignment from the source drawings and proposed 1,2-aryl migration: "
           "new C-S connection, migrated aryl C-C connection, and alcohol C-O bond-order increase to aldehyde. "
           "Exact atom mapping, allene bond-order changes, proton balance, and coproducts are unresolved.")
    r["classification"]["mapped_reaction_smiles"] = missing(ev)
    r["classification"]["bond_changes"] = [derived(x, ev, why) for x in changes]
    r["classification"]["labels"] = derived([
        {"category_number": 4, "code": "C-S", "supporting_change_ids": ["b1"]},
        {"category_number": 2, "code": "C-C", "supporting_change_ids": ["b2"]},
        {"category_number": 14, "code": "DB-C-O", "supporting_change_ids": ["b3"]},
    ], ev, why)
    r["classification"]["coverage"] = inferred("partial", ev, why)
    r["classification"]["limitations"] = [
        "No complete atom mapping is asserted.",
        "Allene bond-order changes, hydrogen/proton transfers, counterions, nitrogen loss, and coproducts are not fully balanced.",
    ]
    r["classification"]["model_votes"] = []


dataset = json.loads(DATASET_PATH.read_text())
dataset["documents"][0]["title"] = (
    "Visible-Light-Mediated Ru-Catalyzed Synthesis of 3-(Arylsulfonyl)but-3-enals "
    "via Coupling of α-Allenols with Diazonium Salts and Sulfur Dioxide"
)
bundle = st.empty(dataset["documents"])
molecules = {}


def add_mol(m):
    molecules[m["molecule_id"]] = m
    return m["molecule_id"]


# Generic structures used only where the paper itself shows R-group classes.
GEN_ALLENE = add_mol(template_molecule("mol_gen_allene", "generic allene", [P[2]], "R-substituted allene", ["R1", "R2", "R3"], "Scheme 1a is a generic prior-work drawing."))
GEN_SULFONYL = add_mol(template_molecule("mol_gen_sulfonyl_reagent", "ArSO2Y", [P[2]], "Ar-SO2-Y", ["Ar", "Y"], "Scheme 1a defines Y = H, Br, I, or CN."))
GEN_SULFONE = add_mol(template_molecule("mol_gen_allylic_sulfone", "generic allylic sulfone product", [P[2]], "allene-derived allylic sulfone", ["R1", "R2", "R3", "Ar", "Y"], "Scheme 1a shows alternative generic products."))
GEN_ALLENOL = add_mol(template_molecule("mol_gen_alpha_allenol", "α-allenol 1", [P[1], P[2], P[3]], "R1-CH(OH)-C(R2)=C=CH2", ["R1", "R2"], "Current-work generic α-allenol."))
GEN_DIAZO = add_mol(template_molecule("mol_gen_diazonium", "arenediazonium tetrafluoroborate 2", [P[1], P[2], P[3]], "Ar-N2+ BF4-", ["Ar"], "Current-work generic diazonium salt."))
GEN_PRODUCT = add_mol(template_molecule("mol_gen_butenal", "3-(arylsulfonyl)but-3-enal 3", [P[1], P[2], P[3]], "OHC-C(R1)(R2)-C(=CH2)-SO2Ar", ["R1", "R2", "Ar"], "Current-work generic product."))
GEN_ALLENOIC = add_mol(template_molecule("mol_gen_allenoic_acid", "2,3-allenoic acid", [P[2]], "substituted 2,3-allenoic acid", ["R1", "R2", "R3"], "Scheme 1d prior-work generic substrate."))
GEN_LACTONE = add_mol(template_molecule("mol_gen_sulfonyl_lactone", "4-sulfonylated furan-2(5H)-one", [P[2]], "sulfonylated lactone", ["R1", "R2", "R3", "Ar"], "Scheme 1d prior-work generic product."))

# Named reagents, catalysts, and solvents.
DABSO = add_mol(molecule("mol_dabso", "DABSO", [P[1], P[2], P[3]], "DABCO·(SO2)2"))
RU = add_mol(molecule("mol_ru", "[Ru]", [P[1], P[2], P[3]], "[Ru(bpy)3](PF6)2"))
IR = add_mol(molecule("mol_ir", "[Ir]", [P[2]], "[Ir(ppy)2(dtbbpy)](PF6)2"))
EOSIN = add_mol(molecule("mol_eosin_y", "Eosin Y", [P[2]], "Eosin Y"))
MECN = add_mol(molecule("mol_mecn", "MeCN", [P[1], P[2], P[3]], "acetonitrile", "CC#N"))
DCE = add_mol(molecule("mol_dce", "DCE", [P[2]], "1,2-dichloroethane", "ClCCCl"))
THF = add_mol(molecule("mol_thf", "THF", [P[2]], "tetrahydrofuran", "C1CCOC1"))
MEOH = add_mol(molecule("mol_meoh", "MeOH", [P[2]], "methanol", "CO"))
TEMPO = add_mol(molecule("mol_tempo", "TEMPO", [P[3]], "TEMPO"))
TOSYLHYD = add_mol(molecule("mol_tosylhydrazide", "TsNHNH2", [P[2]], "p-toluenesulfonyl hydrazide"))
AR_IODIDE = add_mol(template_molecule("mol_gen_aryl_iodide", "ArI", [P[2]], "aryl iodide", ["Ar"], "Scheme 1c generic coupling partner."))
SODIUM_METABISULFITE = add_mol(molecule("mol_sodium_metabisulfite", "Na2S2O5", [P[2]], "sodium metabisulfite", "[Na+].[Na+].[O-]S(=O)S(=O)(=O)[O-]"))


def allenol_smiles(aryl):
    return f"OC({aryl})C(C)=C=C"


def product_smiles(aryl1, aryl2):
    return f"O=CC(C)({aryl1})C(=C)S(=O)(=O){aryl2}"


def diazonium_smiles(aryl):
    return f"N#[N+]{aryl}.[B-](F)(F)(F)F"


aryl_scope2 = {
    "a": ("p-BrC6H4", "c1ccc(Br)cc1"),
    "b": ("C6H5", "c1ccccc1"),
    "c": ("p-MeC6H4", "c1ccc(C)cc1"),
    "d": ("p-MeOC6H4", "c1ccc(OC)cc1"),
    "e": ("p-CO2MeC6H4", "c1ccc(C(=O)OC)cc1"),
    "f": ("o-MeOC6H4", "c1ccccc1OC"),
    "g": ("o-MeC6H4", "c1ccccc1C"),
    "h": ("2,6-(MeO)2C6H3", "c1c(OC)cccc1OC"),
    "i": ("3,5-(MeO)2C6H3", "c1cc(OC)cc(OC)c1"),
    "j": ("2-thienyl", "c1cccs1"),
    "k": ("3-thienyl", "c1ccsc1"),
}
scope2_sub = {}
scope2_prod = {}
pbr = "c1ccc(Br)cc1"
for suffix, (desc, aryl) in aryl_scope2.items():
    ev = [P[2]]
    scope2_sub[suffix] = add_mol(molecule(
        f"mol_1{suffix}_scheme2", f"1{suffix}", ev,
        f"α-allenol 1{suffix} ({desc})", allenol_smiles(aryl),
        "The alcohol stereocenter is drawn without an absolute or relative configuration; no racemate is inferred."))
    pname = "2-(4-bromophenyl)-3-[(4-bromophenyl)sulfonyl]-2-methylbut-3-enal" if suffix == "a" else f"3{suffix}"
    scope2_prod[suffix] = add_mol(molecule(
        f"mol_3{suffix}", f"3{suffix}", ev, pname, product_smiles(aryl, pbr),
        "The newly formed tetrahedral carbon is drawn without stereochemical assignment; no racemate is inferred."))

DIAZO_2A = add_mol(molecule("mol_2a", "2a", [P[1], P[2], P[3]],
    "4-bromobenzenediazonium tetrafluoroborate", diazonium_smiles(pbr)))

# Scheme 3 drawing-defined 1c is kept separate because it conflicts with Scheme 2's p-Me 1c.
pome = "c1ccc(OC)cc1"
SUB_1C_S3 = add_mol(molecule("mol_1c_scheme3_pome", "1c (Scheme 3 drawing)", [P[3]],
    "α-allenol drawn with p-MeOC6H4 in Scheme 3", allenol_smiles(pome),
    "Source-label conflict: Scheme 3 calls this p-MeO drawing 1c, whereas Scheme 2 lists 1c with R = Me. Identities are not merged."))

scope3_specs = {
    "l": ("p-OMe", "c1ccc(OC)cc1"),
    "m": ("p-Cl", "c1ccc(Cl)cc1"),
    "n": ("p-F", "c1ccc(F)cc1"),
    "o": ("p-CO2Et", "c1ccc(C(=O)OCC)cc1"),
    "p": ("p-NO2", "c1ccc([N+](=O)[O-])cc1"),
    "q": ("o-Br", "c1ccccc1Br"),
    "r": ("H", "c1ccccc1"),
}
scope3_diazo = {}
scope3_prod = {}
for suffix, (desc, aryl2) in scope3_specs.items():
    scope3_diazo[suffix] = add_mol(molecule(
        f"mol_diazo_3{suffix}", f"2b-i member for 3{suffix} (R={desc})", [P[3]],
        f"arenediazonium tetrafluoroborate for 3{suffix} ({desc})", diazonium_smiles(aryl2)))
    scope3_prod[suffix] = add_mol(molecule(
        f"mol_3{suffix}", f"3{suffix}", [P[3]], f"3{suffix}", product_smiles(pome, aryl2),
        "The Scheme 3 product drawing has an unassigned tetrahedral stereogenic center; no racemate is inferred."))

# Dense diaryl substrates/products and Supplementary-Scheme products remain source-identified but structurally unresolved.
complex_ids = {}
for label in ["1a-Ar(for 6a)", "1a-Ar(for 6b)", "1b-Ar(for 6c)", "1b-Ar(for 6d)", "1b-Ar(for 6e)",
              "6a", "6b", "6c", "6d", "6e", "7a", "4a", "5a", "5b"]:
    mid = "mol_" + label.lower().replace(" ", "_").replace("(", "_").replace(")", "").replace("-", "_")
    ev = [P[2]] if label in ("4a", "5a", "5b") else [P[3]]
    complex_ids[label] = add_mol(molecule(mid, label, ev, label,
        identity_note="Complete stereochemical SMILES was not safely established from the supplied main-text drawing/prose."))

# Preserve the one explicit geometric assignment without inventing its unavailable structure.
m5b = molecules[complex_ids["5b"]]
m5b["stereochemistry"]["double_bonds"] = direct("(E,E)-5b", [{"configuration": "E"}, {"configuration": "E"}], [P[2]])
m5b["stereochemistry"]["geometry_status"] = direct("(E,E)-5b", "specified", [P[2]])
m5b["stereochemistry"]["representation_status"] = inferred("not_encoded", [P[2]],
    "The source assigns E,E geometry, but the complete structure is unavailable on the supplied main-text pages.")
m5b["stereochemistry"]["assignment"] = direct(
    "the structure of (E,E)-5b was assigned by single-crystal X-ray diffraction analysis",
    {"geometry": "E,E", "method": "single-crystal X-ray diffraction"}, [P[2]])

bundle["molecules"] = list(molecules.values())
reactions = []


def add_reaction(r):
    reactions.append(r)
    return r["reaction_id"]


# Prior-work drawings and current generic schemes.
def simple_generic(rid, record_type, location, reactants, products, ev, name):
    r = reaction(rid, record_type, ev, None, location)
    r["reaction_name"] = direct(name, name, ev)
    for i, mid in enumerate(reactants):
        r["participants"].append(participant(f"{rid}_r{i+1}", mid, "reactant", ev))
    for i, mid in enumerate(products):
        r["participants"].append(participant(f"{rid}_p{i+1}", mid, "product", ev, product_role="depicted product"))
    unresolved_classification(r, ev, "unresolved", "Generic prior-work drawing has alternative products or unresolved atom correspondence.")
    return add_reaction(r)


LIT_A = simple_generic("scheme1a_prior_work", "literature_example", "page 2 Scheme 1a",
    [GEN_ALLENE, GEN_SULFONYL], [GEN_SULFONE], [P[2]], "previous syntheses of sulfones from allenes")
LIT_B = simple_generic("scheme1b_prior_work", "literature_example", "page 2 Scheme 1b",
    [GEN_ALLENE, TOSYLHYD], [GEN_SULFONE], [P[2]], "hydrosulfonylation of allenes with sulfonyl hydrazides")
LIT_C = simple_generic("scheme1c_prior_work", "literature_example", "page 2 Scheme 1c",
    [GEN_ALLENOL, TOSYLHYD, AR_IODIDE], [GEN_SULFONE], [P[2]], "palladium-catalyzed three-component reaction")
LIT_D = simple_generic("scheme1d_prior_work", "literature_example", "page 2 Scheme 1d",
    [GEN_ALLENOIC, GEN_DIAZO, SODIUM_METABISULFITE], [GEN_LACTONE], [P[2]], "copper-catalyzed oxysulfonylation of allenoic acids")


def current_generic(rid, location, ev):
    r = reaction(rid, "generic_scheme", ev, None, location)
    r["reaction_name"] = direct("sulfonylation-rearrangement cascade", "sulfonylation-rearrangement cascade", ev)
    r["transformation_type"] = inferred("visible-light three-component sulfonylation with 1,2-aryl migration", ev,
        "Transformation description combines the source title/caption and product drawing.")
    for pid, mid, role in [("allenol", GEN_ALLENOL, "reactant"), ("diazo", GEN_DIAZO, "reactant"),
                           ("dabso", DABSO, "reagent"), ("ru", RU, "catalyst"),
                           ("mecn", MECN, "solvent"), ("product", GEN_PRODUCT, "product")]:
        r["participants"].append(participant(f"{rid}_{pid}", mid, role, ev, product_role="product class" if role == "product" else None))
    set_charge(r["participants"][2], "equivalents", "2 equiv", qty(2, "equiv", "2 equiv"), ev, GEN_ALLENOL)
    set_charge(r["participants"][3], "mol_percent", "2 mol%", qty(2, "mol%", "2 mol%"), ev, GEN_ALLENOL)
    r["solvent_conditions"]["composition"] = direct("MeCN", [MECN], ev)
    add_step(r, "reaction", ev, [p["participant_id"] for p in r["participants"]],
             "visible-light irradiation", "RT", None, "Ar")
    r["special_conditions"]["light_source"] = direct("visible light / LED lamp", "visible light", ev)
    cascade_classification(r, ev)
    return add_reaction(r)


GEN_CURRENT = current_generic("scheme1e_current_generic", "page 1 graphical abstract and page 2 Scheme 1e", [P[1], P[2]])
GEN_S2 = current_generic("scheme2_generic", "page 2 Scheme 2", [P[2]])
GEN_S3 = current_generic("scheme3_generic", "page 3 Scheme 3", [P[3]])


solvent_ids = {"MeCN": MECN, "DCE": DCE, "THF": THF, "MeOH": MEOH}
table = [
    (1, IR, "MeCN", 16, 48, None),
    (2, RU, "MeCN", 20, 59, "optimal; duplicate depiction as Scheme 2 product 3a"),
    (3, EOSIN, "MeCN", 16, 21, None),
    (4, RU, "DCE", 20, 0, None),
    (5, RU, "THF", 20, 0, None),
    (6, RU, "MeOH", 20, 0, None),
    (7, RU, "MeCN", 18, 31, "DABSO 0.2 mmol"),
    (8, RU, "MeCN", 18, 35, "2a 0.2 mmol"),
    (9, None, "MeCN", 20, 0, "no photocatalyst"),
    (10, RU, "MeCN", 20, 0, "darkness"),
    (11, RU, "MeCN", 10, 9, "daylight on a sunny day"),
    (12, RU, "MeCN", 10, 30, "air instead of argon"),
]


def add_table_entry(n, catalyst, solvent, hours, yld, note):
    ev = [P[2]]
    if n == 2:
        ev = [P[1], P[2]]
    if n in (9, 10):
        ev = [P[2], P[3]]
    rid = f"table1_entry{n:02d}"
    categories = ["optimization"] + (["substrate_scope"] if n == 2 else []) + (["control", "mechanistic_study"] if n in (9, 10) else [])
    r = reaction(rid, "specific_experiment", ev, categories, f"page 2 Table 1 entry {n}" + ("; page 2 Scheme 2 3a" if n == 2 else ""))
    r["general_scheme_ids"] = [GEN_CURRENT] + ([GEN_S2] if n == 2 else [])
    r["reaction_name"] = direct("Sulfonylation-Rearrangement Cascade of Allenol 1a with Arenediazonium Salt 2a", "sulfonylation-rearrangement cascade", [P[2]])
    r["transformation_type"] = inferred("visible-light three-component sulfonylation with 1,2-aryl migration", ev,
        "Source scheme/prose interpretation; not a complete atom mapping.")
    parts = [
        participant(f"{rid}_1a", scope2_sub["a"], "reactant", ev),
        participant(f"{rid}_2a", DIAZO_2A, "reactant", ev),
        participant(f"{rid}_dabso", DABSO, "reagent", ev),
    ]
    if catalyst:
        parts.append(participant(f"{rid}_pc", catalyst, "catalyst", ev))
    parts += [participant(f"{rid}_solvent", solvent_ids[solvent], "solvent", ev),
              participant(f"{rid}_3a", scope2_prod["a"], "product", ev, product_role="intended/main product")]
    r["participants"] = parts
    set_charge(parts[0], "amount", "0.2 mmol", qty(0.2, "mmol", "0.2 mmol"), [P[2]])
    amount_2a = 0.2 if n == 8 else 0.3
    set_charge(parts[1], "amount", f"{amount_2a} mmol", qty(amount_2a, "mmol", f"{amount_2a} mmol"), [P[2]])
    set_charge(parts[1], "equivalents", "1.0 equiv" if n == 8 else "1.5 equiv",
               qty(1.0 if n == 8 else 1.5, "equiv", "1.0 equiv" if n == 8 else "1.5 equiv"), [P[1], P[2]] if n != 8 else [P[2]], scope2_sub["a"])
    amount_dabso = 0.2 if n == 7 else 0.4
    set_charge(parts[2], "amount", f"{amount_dabso} mmol", qty(amount_dabso, "mmol", f"{amount_dabso} mmol"), [P[2]])
    set_charge(parts[2], "equivalents", "1 equiv" if n == 7 else "2 equiv",
               qty(1 if n == 7 else 2, "equiv", "1 equiv" if n == 7 else "2 equiv"), [P[2]], scope2_sub["a"])
    if catalyst:
        set_charge(parts[3], "mol_percent", "2 mol%", qty(2, "mol%", "2 mol%"), [P[2]], scope2_sub["a"])
    r["reaction_scale"] = direct("1a (0.2 mmol)", {"reference_molecule_id": scope2_sub["a"], "amount": qty(0.2, "mmol", "0.2 mmol")}, [P[2]])
    r["solvent_conditions"]["composition"] = direct(solvent, [solvent_ids[solvent]], [P[2]])
    if n == 2:
        r["solvent_conditions"]["concentration"] = inferred(qty(1/30, "M", "30 mL/mmol 1a"), [P[1]],
            "30 mL solvent per mmol 1a corresponds to 0.0333 mol/L.", raw="MeCN (30 mL/mmol 1a)")
        r["solvent_conditions"]["concentration_basis"] = direct("per mmol 1a", scope2_sub["a"], [P[1]])
        r["solvent_conditions"]["total_volume"] = inferred(qty(6, "mL", "calculated from 0.2 mmol × 30 mL/mmol"), [P[1], P[2]],
            "0.2 mmol 1a multiplied by 30 mL/mmol gives 6 mL; this is calculated, not an author-reported total volume.")
    atmosphere = "air" if n == 12 else "Ar"
    op = "reaction carried out in darkness" if n == 10 else "visible-light irradiation"
    if n == 11:
        op = "daylight on a sunny day"
    add_step(r, "reaction", ev, [p["participant_id"] for p in parts], op, "RT", hours, atmosphere)
    light = "darkness" if n == 10 else "daylight on a sunny day" if n == 11 else "visible light"
    r["special_conditions"]["light_source"] = direct(light, light, [P[2]])
    outcome(r, "yield_3a", "yield", [scope2_prod["a"]], f"{yld}%", qty(yld, "%", f"{yld}%"), [P[2]],
            role="main product", yield_type="pure isolated product", basis="isolated; correct analytical and spectral data")
    if n == 2:
        conversion = outcome(r, "conversion_1a", "conversion", [scope2_sub["a"]], "complete conversion", "complete conversion", [P[2]], scope="starting allenol consumption")
        conversion["measurement_basis"] = direct("Reaction progress was followed by TLC until total consumption of the starting allenol.", "TLC", [P[2]])
        r["notes"].append(inferred("Table 1 entry 2 and Scheme 2 product 3a are one experiment represented once with both contexts.", ev,
            "Identical substrate/product, catalyst, solvent, atmosphere, temperature, time, and 59% yield."))
    if note:
        r["notes"].append(direct(note, note, [P[2]]))
    if yld == 0:
        unresolved_classification(r, ev, "not_applicable", "Zero yield/no observed sulfonyl-butenal; no expected-product bond labels are asserted.")
        r["review"]["unresolved"].append("The source reports 0% in the 3a yield column; no product formation is treated as observed.")
    else:
        cascade_classification(r, ev)
    add_reaction(r)


for row in table:
    add_table_entry(*row)


def add_scope(rid, substrate, diazo, products, yields, hours, ev, scheme_id, location, notes=None):
    r = reaction(rid, "specific_experiment", ev, ["substrate_scope"], location)
    r["general_scheme_ids"] = [scheme_id]
    r["reaction_name"] = direct("sulfonylation-rearrangement cascade", "sulfonylation-rearrangement cascade", ev)
    r["transformation_type"] = inferred("visible-light three-component sulfonylation with 1,2-aryl migration", ev,
        "Source scheme/prose interpretation; no complete atom mapping is asserted.")
    parts = [participant(f"{rid}_sub", substrate, "reactant", ev),
             participant(f"{rid}_diazo", diazo, "reactant", ev),
             participant(f"{rid}_dabso", DABSO, "reagent", ev),
             participant(f"{rid}_ru", RU, "catalyst", ev),
             participant(f"{rid}_mecn", MECN, "solvent", ev)]
    for i, (mid, role) in enumerate(products):
        parts.append(participant(f"{rid}_p{i+1}", mid, "product", ev, product_role=role))
    r["participants"] = parts
    set_charge(parts[1], "equivalents", "1.5 equiv", qty(1.5, "equiv", "1.5 equiv"), [P[1], *ev], substrate,
               method="inherited_general_condition", why="Inherited from the explicit optimal-condition prose on page 1.")
    set_charge(parts[2], "equivalents", "2 equiv", qty(2, "equiv", "2 equiv"), ev, substrate)
    set_charge(parts[3], "mol_percent", "2 mol%", qty(2, "mol%", "2 mol%"), ev, substrate)
    r["solvent_conditions"]["composition"] = direct("MeCN", [MECN], ev)
    r["solvent_conditions"]["concentration"] = F("30 mL/mmol allenol", qty(1/30, "M", "30 mL/mmol allenol"), [P[1], *ev],
        method="inherited_general_condition", derivation="Optimal-condition solvent ratio corresponds to 0.0333 M and is inherited; scope scale is unreported.")
    r["solvent_conditions"]["concentration_basis"] = F("per mmol allenol", substrate, [P[1], *ev],
        method="inherited_general_condition", derivation="Inherited optimal-condition basis.")
    add_step(r, "reaction", ev, [p["participant_id"] for p in parts], "visible-light irradiation", "RT", hours, "Ar")
    r["special_conditions"]["light_source"] = F("LED light bulb / visible light", "visible light", [P[1], *ev],
        method="inherited_general_condition", derivation="LED light source explicitly stated for the optimal reaction and inherited by the scope scheme.")
    for i, (mid, role) in enumerate(products):
        y = yields[i]
        outcome(r, f"yield_{i+1}", "yield", [mid], f"{y}%", qty(y, "%", f"{y}%"), ev,
                role=role, yield_type="isolated yield", basis="isolated")
    conversion = outcome(r, "conversion", "conversion", [substrate], "complete conversion", "complete conversion", ev,
                         scope="starting allenol consumption")
    if scheme_id == GEN_S2:
        conversion["measurement_basis"] = direct(
            "Reaction progress was followed by TLC until total consumption of the starting allenol.", "TLC", ev)
    if notes:
        for note in notes:
            r["notes"].append(direct(note, note, ev))
    cascade_classification(r, ev)
    add_reaction(r)


scope2_yields = {"b": (49, 16), "c": (59, 16), "d": (72, 16), "e": (61, 16),
                 "f": (61, 16), "g": (64, 16), "h": (53, 16), "i": (61, 16),
                 "j": (62, 1.5), "k": (60, 6)}
for suffix, (y, h) in scope2_yields.items():
    add_scope(f"scheme2_3{suffix}", scope2_sub[suffix], DIAZO_2A,
              [(scope2_prod[suffix], "main product")], [y], h, [P[2]], GEN_S2,
              f"page 2 Scheme 2 product 3{suffix}")


scope3_yields = {"l": (51, 16), "m": (43, 1), "n": (69, 1), "o": (38, 2),
                 "p": (42, 3), "q": (59, 1), "r": (52, 2)}
for suffix, (y, h) in scope3_yields.items():
    add_scope(f"scheme3_3{suffix}", SUB_1C_S3, scope3_diazo[suffix],
              [(scope3_prod[suffix], "main product")], [y], h, [P[3]], GEN_S3,
              f"page 3 Scheme 3 product 3{suffix}",
              ["Source-label conflict retained: Scheme 3 draws 1c as p-MeO, while Scheme 2 lists 1c as R = Me."])

# Diaryl examples. Exact precursor/product SMILES are left unresolved; 6a and 7a are one experiment.
diaryl_rows = [
    ("6a", "1a-Ar(for 6a)", DIAZO_2A, [("6a", "major product"), ("7a", "minor product")], [28, 10], 24),
    ("6b", "1a-Ar(for 6b)", add_mol(molecule("mol_diazo_mbr", "2i (m-Br drawing)", [P[3]], "3-bromobenzenediazonium tetrafluoroborate", diazonium_smiles("c1cccc(Br)c1"))), [("6b", "main product")], [41], 5),
    ("6c", "1b-Ar(for 6c)", DIAZO_2A, [("6c", "main product")], [52], 16),
    ("6d", "1b-Ar(for 6d)", DIAZO_2A, [("6d", "main product")], [46], 16),
    ("6e", "1b-Ar(for 6e)", scope3_diazo["q"], [("6e", "main product")], [45], 2),
]
# One molecule was added after the initial bundle assignment.
bundle["molecules"] = list(molecules.values())
for label, sublabel, diazo, product_rows, ys, h in diaryl_rows:
    add_scope(f"scheme3_{label}", complex_ids[sublabel], diazo,
              [(complex_ids[p], role) for p, role in product_rows], ys, h, [P[3]], GEN_S3,
              f"page 3 Scheme 3 product {label}" + (" and side product 7a" if label == "6a" else ""),
              ["7a was obtained as a chromatographically separable minor product."] if label == "6a" else None)

# Main-text qualitative homodimerization; SI Scheme S1 is unavailable.
r = reaction("homodimerization_4a", "specific_experiment", [P[2]], ["substrate_scope"], "page 2 prose; Supplementary Scheme S1 not supplied")
r["reaction_name"] = direct("conversion of oxindole-tethered quaternary allenol 4a into homodimeric derivatives 5a and 5b", "homodimerization of allenol 4a", [P[2]])
r["participants"] = [
    participant("dimer_4a", complex_ids["4a"], "reactant", [P[2]]),
    participant("dimer_2a", DIAZO_2A, "reactant", [P[2]]),
    participant("dimer_dabso", DABSO, "reagent", [P[2]]),
    participant("dimer_ru", RU, "catalyst", [P[2]]),
    participant("dimer_mecn", MECN, "solvent", [P[2]]),
    participant("dimer_5a", complex_ids["5a"], "product", [P[2]], product_role="homodimeric derivative"),
    participant("dimer_5b", complex_ids["5b"], "product", [P[2]], product_role="homodimeric derivative; E,E assignment"),
]
set_charge(r["participants"][1], "equivalents", "1.5 equiv (otherwise identical conditions)", qty(1.5, "equiv", "1.5 equiv"), [P[1], P[2]], complex_ids["4a"],
           method="inherited_general_condition", why="Page 2 states otherwise identical conditions; exact SI procedure is unavailable.")
set_charge(r["participants"][2], "equivalents", "2 equiv (otherwise identical conditions)", qty(2, "equiv", "2 equiv"), [P[1], P[2]], complex_ids["4a"],
           method="inherited_general_condition", why="Page 2 states otherwise identical conditions; exact SI procedure is unavailable.")
set_charge(r["participants"][3], "mol_percent", "2 mol% (otherwise identical conditions)", qty(2, "mol%", "2 mol%"), [P[1], P[2]], complex_ids["4a"],
           method="inherited_general_condition", why="Page 2 states otherwise identical conditions; exact SI procedure is unavailable.")
r["solvent_conditions"]["composition"] = F("MeCN (otherwise identical conditions)", [MECN], [P[1], P[2]], method="inherited_general_condition", derivation="Inherited from the stated optimal conditions.")
add_step(r, "reaction", [P[1], P[2]], [p["participant_id"] for p in r["participants"]],
         "otherwise identical visible-light conditions", "RT", None, "Ar", operation_method="inherited_general_condition")
r["special_conditions"]["light_source"] = F("otherwise identical visible-light conditions", "visible light", [P[1], P[2]], method="inherited_general_condition", derivation="Inherited from the explicit wording 'under otherwise identical conditions'.")
outcome(r, "formation_5a_5b", "other", [complex_ids["5a"], complex_ids["5b"]],
        "smoothly converted into homodimeric derivatives 5a and 5b", "formed", [P[2]], scope="qualitative formation; yields unreported in main text")
unresolved_classification(r, [P[2]], "unresolved", "Homodimer structures and atom correspondence are only in unavailable Supplementary Scheme S1.")
r["review"]["unresolved"].append("Whether every optimal-condition reagent/loading applies unchanged needs the unavailable Supplementary Scheme S1; inherited values are explicitly marked.")
add_reaction(r)

# TEMPO control described in prose; exact Supplementary Scheme S2 identities/charges are unavailable.
r = reaction("tempo_control", "specific_experiment", [P[3]], ["control", "mechanistic_study"], "page 3 control-experiment prose; Supplementary Scheme S2 not supplied")
r["reaction_name"] = direct("radical inhibition experiment with TEMPO", "TEMPO inhibition control", [P[3]])
r["participants"] = [
    participant("tempo_allenol", GEN_ALLENOL, "reactant", [P[3]]),
    participant("tempo_diazo", GEN_DIAZO, "reactant", [P[3]]),
    participant("tempo_dabso", DABSO, "reagent", [P[3]]),
    participant("tempo_ru", RU, "catalyst", [P[3]]),
    participant("tempo_additive", TEMPO, "additive", [P[3]]),
]
add_step(r, "reaction", [P[3]], [p["participant_id"] for p in r["participants"]], "TEMPO present in the reaction medium")
outcome(r, "suppression", "no_reaction", [], "construction of the sulfonyl-butenal scaffold was suppressed", "suppressed", [P[3]], scope="qualitative radical-inhibition observation")
unresolved_classification(r, [P[3]], "not_applicable", "No observed product transformation in the TEMPO-inhibited control.")
r["review"]["unresolved"].append("Exact substrate identity, TEMPO loading, time, temperature, and detection basis are only in unavailable Supplementary Scheme S2.")
add_reaction(r)

# Proposed mechanism is not promoted to experimental proof.
r = reaction("scheme4_proposed_mechanism", "proposed_mechanism", [P[3]], None, "page 3 Scheme 4")
r["reaction_name"] = direct("Tentative Mechanistic Explanation for the Multicomponent Reaction between DABSO, Diazonium Salts, and α-Allenols", "tentative radical sulfonylation-rearrangement mechanism", [P[3]])
for pid, mid, role in [("m_allenol", GEN_ALLENOL, "reactant"), ("m_diazo", GEN_DIAZO, "reactant"),
                       ("m_dabso", DABSO, "reagent"), ("m_ru", RU, "catalyst"), ("m_product", GEN_PRODUCT, "product")]:
    r["participants"].append(participant(pid, mid, role, [P[3]], product_role="proposed final product class" if role == "product" else None))
r["procedure"]["original_text"] = direct(
    "arylsulfonyl radical addition; SET; 1,2-aryl migration; aldehyde formation; deprotonation",
    ["arylsulfonyl radical addition", "SET oxidation", "1,2-aryl migration", "aldehyde formation", "deprotonation"], [P[3]])
r["notes"].append(direct("The authors call Scheme 4 a tentative mechanistic explanation.", "author-proposed, not experimentally demonstrated sequence", [P[3]]))
unresolved_classification(r, [P[3]], "unresolved", "Mechanistic intermediates I–VI are proposed; the record is not treated as an observed experiment.")
add_reaction(r)

bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": f"{DOC}:control_discussion", "document_id": DOC, "evidence_ids": [P[3]],
     "kind": "mechanistic_interpretation",
     "text": "Authors interpret no-photocatalyst, darkness, and TEMPO experiments as support for arylsulfonyl radicals; Table 1 entries 9 and 10 hold the first two observed controls, and a separate partial TEMPO record is retained."},
    {"context_id": f"{DOC}:xray_5b", "document_id": DOC, "evidence_ids": [P[2], P[3]],
     "kind": "characterization",
     "text": "(E,E)-5b was assigned by single-crystal X-ray diffraction; CCDC 1900824 is cited. This is characterization context, not another synthesis."},
    {"context_id": f"{DOC}:si_notice", "document_id": DOC, "evidence_ids": [P[3]],
     "kind": "supplementary_information_notice",
     "text": "The article advertises an SI PDF, but no SI file was supplied or explicitly associated in this session."},
    {"context_id": f"{DOC}:bibliography_pages_4_5", "document_id": DOC, "evidence_ids": [P[4], P[5]],
     "kind": "non_experimental_content",
     "text": "Pages 4-5 are acknowledgments/dedication/references and footnotes; cited literature is not imported as experiments newly reported in this paper."},
]

# Raw page inventory and the required coverage contract.
page_records = {
    1: [GEN_CURRENT, f"{DOC}:table1_entry02"],
    2: [LIT_A, LIT_B, LIT_C, LIT_D, GEN_CURRENT, GEN_S2]
       + [f"{DOC}:table1_entry{i:02d}" for i in range(1, 13)]
       + [f"{DOC}:scheme2_3{x}" for x in "bcdefghijk"]
       + [f"{DOC}:homodimerization_4a"],
    3: [GEN_S3] + [f"{DOC}:scheme3_3{x}" for x in "lmnopqr"]
       + [f"{DOC}:scheme3_6{x}" for x in "abcde"]
       + [f"{DOC}:table1_entry09", f"{DOC}:table1_entry10", f"{DOC}:tempo_control", f"{DOC}:scheme4_proposed_mechanism"],
    4: [],
    5: [],
}
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title/abstract and graphical abstract", "finding": "current-work generic three-component α-allenol/diazonium/DABSO coupling; 23 examples and up to 72% are summary context", "evidence_ids": [P[1]]},
            {"region": "right-column optimization prose", "finding": "explicit optimal 1a/2a/DABSO/Ru conditions and 59% 3a result; mapped to the merged Table 1 entry 2/Scheme 2 3a record", "evidence_ids": [P[1], T[1]]},
        ], "record_ids": page_records[1], "unresolved": ["No SI PDF supplied; detailed procedure and workup are unavailable."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Scheme 1a-d", "finding": "four prior-work generic reaction arrows; classified as literature examples, not current experiments", "evidence_ids": [P[2]]},
            {"region": "Scheme 1e", "finding": "current-work generic sulfonylation-rearrangement proposal", "evidence_ids": [P[2]]},
            {"region": "Table 1 rows 1-12 and footnotes a-h", "finding": "twelve independent optimization/control experiments; all row values and footnote modifications transcribed", "evidence_ids": [P[2], T[2]]},
            {"region": "Scheme 2 products 3a-k", "finding": "eleven scope outcomes; 3a is merged with identical Table 1 entry 2 and 3b-k are separate experiments", "evidence_ids": [P[2]]},
            {"region": "lower-right prose", "finding": "qualitative 4a homodimerization to 5a/5b retained as a partial experiment; (E,E)-5b X-ray assignment retained as characterization", "evidence_ids": [P[2], T[2]]},
        ], "record_ids": page_records[2], "unresolved": ["Scheme S1 yields/structures and exact dimerization procedure are unavailable.", "Scheme 2 assigns 1c as R=Me, conflicting with the p-MeO drawing labeled 1c in Scheme 3."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Scheme 3 products 3l-r", "finding": "seven substrate-scope experiments with isolated yields and times", "evidence_ids": [P[3]]},
            {"region": "Scheme 3 products 6a-e and 7a", "finding": "five diaryl scope experiments; 6a major and chromatographically separable 7a minor are one experiment with two product yields", "evidence_ids": [P[3], T[3]]},
            {"region": "control-experiment prose", "finding": "no-photocatalyst and darkness observations cross-linked to Table 1 entries 9/10; distinct TEMPO suppression retained as a partial mechanistic control", "evidence_ids": [P[3], T[3]]},
            {"region": "Scheme 4", "finding": "tentative author-proposed radical/SET/1,2-aryl-migration mechanism, not an observed experimental sequence", "evidence_ids": [P[3]]},
            {"region": "Associated Content/Accession Codes", "finding": "SI availability notice and CCDC 1900824 characterization context; no SI association inferred", "evidence_ids": [P[3]]},
        ], "record_ids": page_records[3], "unresolved": ["Supplementary Schemes S2/S3 are not supplied; exact TEMPO control and 7a rationale remain incomplete.", "Scheme 3's p-MeO drawing labeled 1c conflicts with Scheme 2's 1c R=Me and is kept as a separate source identity."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "acknowledgments, dedication, and references 1-10", "finding": "no experiment newly reported in this paper; bibliography not imported as current reactions", "evidence_ids": [P[4], T[4]]},
        ], "record_ids": [], "unresolved": []},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "references 10-18 and notes", "finding": "bibliography plus caution/rationale footnotes; no new current-paper reaction arrow or experiment", "evidence_ids": [P[5], T[5]]},
        ], "record_ids": [], "unresolved": ["Footnotes refer to unavailable Supplementary Scheme S3 but do not provide its experiment details."]},
    ],
}

raw = {
    "document_id": DOC,
    "pages_reviewed": [1, 2, 3, 4, 5],
    "visual_review": True,
    "chemical_accuracy_verified": False,
    "reaction_ids": [r["reaction_id"] for r in reactions],
    "molecule_ids": [m["molecule_id"] for m in bundle["molecules"]],
    "page_inventory": coverage["pages"],
    "decisions": [
        "Merged Table 1 entry 2 and Scheme 2 3a as the same experiment.",
        "Kept Scheme 3 p-MeO '1c' separate from Scheme 2 p-Me '1c' because the supplied source is internally inconsistent.",
        "Stored 6a and 7a as major/minor outcomes of one reaction.",
        "Did not import cited references, spectra, or CCDC context as synthesis experiments.",
        "Did not infer or fetch Supporting Information association.",
    ],
}

RUN.joinpath("coverage").mkdir(parents=True, exist_ok=True)
st.write(RUN / "submitted.json", bundle)
st.write(RUN / "raw_model/pages_0001_0005_extraction.json", raw)
st.write(RUN / "coverage/pages_0001_0005.json", coverage)
st.write(DATASET_PATH, dataset)

# Mark only the five mechanically prepared page tasks as attempted; this is not semantic verification.
tasks_path = RUN / "tasks.json"
tasks = json.loads(tasks_path.read_text())
for task in tasks:
    task["state"] = "attempted"
    task["reason"] = "Original page image and extracted text inspected; see coverage/pages_0001_0005.json."
st.write(tasks_path, tasks)

stereo_checks = {}
for m in bundle["molecules"]:
    smi = st.value(m["smiles"])
    if smi:
        stereo_checks[m["molecule_id"]] = st.stereo_inspect(smi)
st.write(RUN / "raw_model/stereo_checks.json", stereo_checks)

print(json.dumps({
    "molecules": len(bundle["molecules"]),
    "reactions": len(reactions),
    "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in reactions),
    "generic_schemes": sum(r["record_type"] == "generic_scheme" for r in reactions),
    "literature_examples": sum(r["record_type"] == "literature_example" for r in reactions),
    "proposed_mechanisms": sum(r["record_type"] == "proposed_mechanism" for r in reactions),
}, ensure_ascii=False, indent=2))
