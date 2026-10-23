#!/usr/bin/env python3
"""Build the page-1--6 machine submission from visually inspected source pages."""
from pathlib import Path
import copy
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
TOOLS = ROOT / "skills" / "chemical-reaction-extraction" / "scripts"
sys.path.insert(0, str(TOOLS))
from standard_records import (  # noqa: E402
    STANDARD,
    bond_label,
    empty,
    fact,
    new_molecule,
    new_outcome,
    new_participant,
    new_reaction,
    new_step,
    stereo_inspect,
    write,
)
from reaction_tools import dataset_audit  # noqa: E402

DOC = "doc_dfa009e2c8318d59"
RUN = ROOT / "runs" / "current"
DATASET_PATH = RUN / "dataset.json"
SUBMITTED_PATH = RUN / "submitted.json"


def ev(page, image=True, text=True):
    ids = []
    if image:
        ids.append(f"{DOC}:page{page}:image")
    if text:
        ids.append(f"{DOC}:page{page}:text")
    return ids


def src(raw, normalized=None, evidence=(), method="original_extraction", derivation=None):
    return fact(
        raw=raw,
        normalized=raw if normalized is None else normalized,
        evidence=evidence,
        method=method,
        derivation=derivation,
    )


def inferred(normalized, evidence, derivation, raw=None):
    return fact(
        raw=raw,
        normalized=normalized,
        evidence=evidence,
        method="inference",
        derivation=derivation,
    )


def derived(normalized, evidence, derivation, raw=None):
    return fact(
        raw=raw,
        normalized=normalized,
        evidence=evidence,
        method="structure_derivation",
        derivation=derivation,
    )


def missing(evidence=(), status="unreported"):
    return fact(evidence=evidence, status=status)


def not_applicable(evidence=()):
    return fact(evidence=evidence, method="not_applicable", status="not_applicable")


def quantity(number, unit, raw_text, qualifier=None):
    value = {"number": number, "unit": unit, "raw_text": raw_text}
    if qualifier:
        value["qualifier"] = qualifier
    return value


def qrange(minimum, maximum, unit, raw_text):
    return {"minimum": minimum, "maximum": maximum, "unit": unit, "raw_text": raw_text}


def make_molecule(mid, label, name, evidence, smiles=None, chirality="unresolved",
                  geometry="unresolved", representation=None, note=None, axial=None,
                  atom_centers=None):
    m = new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = src(label, label, evidence)
    m["names"]["original"] = src(name, name, evidence)
    m["names"]["english"] = src(name, name, evidence)
    m["names"]["abbreviations"] = missing(evidence)
    m["names"]["aliases"] = missing(evidence)
    m["cas_numbers"] = missing(evidence)
    m["external_properties"] = []
    st = m["stereochemistry"]
    if smiles:
        diag = stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = derived(
            diag["canonical_isomeric_smiles"], evidence,
            "Candidate connectivity was derived from the source name/drawing and normalized with RDKit; parsing is not source verification.",
            raw={"source_label": label, "candidate_smiles": smiles},
        )
        m["generic_template"] = None
        st["atom_centers"] = derived(diag["atom_centers"], evidence, "RDKit graph diagnostic on the candidate SMILES.", raw=label)
        st["double_bonds"] = derived(diag["double_bonds"], evidence, "RDKit graph diagnostic; a listed double bond is not necessarily E/Z-stereogenic.", raw=label)
        st["axial"] = derived(axial or [], evidence, "Source drawing and graph were checked for an explicitly assigned axis; toolkit axial support is limited.", raw=label)
        st["mixture"] = src("No stereoisomeric mixture is reported for this named reagent/solvent.", None, evidence)
        st["assignment"] = derived("No configurational assignment is applicable to the candidate graph.", evidence, "No potential specified stereoelement was returned for the candidate graph.", raw=label)
        st["chirality_status"] = derived(chirality, evidence, "Status is based on the candidate graph and the source drawing; it is not inferred merely from absent wedge bonds.", raw=label)
        st["geometry_status"] = derived(geometry, evidence, "No E/Z-capable double bond is present in the candidate graph unless separately stated.", raw=label)
        has_unspecified = any(p["specified"] == "Unspecified" for p in diag["potential_stereo"])
        effective_representation = representation or ("partial" if has_unspecified else "encoded")
        st["representation_status"] = derived(effective_representation, evidence, "The registered candidate SMILES is marked partial whenever the graph diagnostic reports an unspecified potential stereoelement.", raw=label)
        st["isomeric_smiles"] = derived(diag["canonical_isomeric_smiles"], evidence, "Canonical isomeric SMILES generated from the candidate graph.", raw=smiles)
        st["identity_key"] = derived(diag["stereo_graph_key"], evidence, "RDKit stereo graph key; not a sample identity.", raw=smiles)
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing(evidence)
        m["generic_template"] = None
        st["atom_centers"] = derived(atom_centers, evidence, "Potential centers were assessed from the visible source drawing; no absolute assignment is claimed.", raw=label) if atom_centers is not None else missing(evidence)
        st["double_bonds"] = missing(evidence)
        st["axial"] = derived(axial, evidence, "Allene/axial status was assessed from the visible source drawing.", raw=label) if axial is not None else missing(evidence)
        st["mixture"] = missing(evidence)
        st["assignment"] = src("No wedge/dash, R/S, Ra/Sa, E/Z, or cis/trans assignment is shown for this occurrence.", "unspecified", evidence)
        st["chirality_status"] = inferred(chirality, evidence, "The source does not specify configuration; unspecified is not treated as achiral or racemic.", raw=label)
        st["geometry_status"] = inferred(geometry, evidence, "No source assignment is available for geometric stereochemistry.", raw=label)
        st["representation_status"] = inferred(representation or "not_encoded", evidence, "A complete reliable SMILES was not asserted for the source drawing.", raw=label)
        st["isomeric_smiles"] = missing(evidence)
        st["identity_key"] = missing(evidence)
    m["identity_note"] = src(note, note, evidence) if note else missing(evidence)
    return m


def participant(pid, mid, role, evidence, product_role=None, stage="reaction"):
    p = new_participant(pid, mid, role, evidence, stage=stage)
    p["product_role"] = src(product_role, product_role, evidence) if product_role else not_applicable(evidence)
    return p


def add_equiv(p, number, raw_text, reference_mid, evidence, approximate=False):
    p["charges"]["equivalents"] = src(
        raw_text,
        quantity(number, "equiv", raw_text, "approx" if approximate else None),
        evidence,
    )
    p["equivalence_reference"] = inferred(
        reference_mid,
        evidence,
        "The reported equivalent is associated with the stoichiometric allene/adduct sequence; the page does not explicitly restate the denominator.",
        raw="equivalent basis not explicitly restated",
    )


def add_mol_percent(p, number, raw_text, reference_mid, evidence, approximate=False):
    p["charges"]["mol_percent"] = src(
        raw_text,
        quantity(number, "mol%", raw_text, "approx" if approximate else None),
        evidence,
    )
    p["equivalence_reference"] = inferred(
        reference_mid,
        evidence,
        "CuBr loading is reported for cyclization of the 1-aza-1,3,4-triene; the precursor originates stoichiometrically from allene 1.",
        raw="loading denominator not explicitly restated in Scheme 1",
    )


def step(step_id, evidence, pids, operation, temperature=None, time=None, stage="reaction",
         atmosphere=None):
    s = new_step(step_id, evidence, stage=stage)
    s["participant_ids"] = list(pids)
    s["operation"] = src(operation, operation, evidence)
    s["temperature"] = src(temperature[0], temperature[1], evidence) if temperature else missing(evidence)
    s["time"] = src(time[0], time[1], evidence) if time else missing(evidence)
    s["atmosphere"] = src(atmosphere, atmosphere, evidence) if atmosphere else missing(evidence)
    s["pressure"] = missing(evidence)
    return s


def base_reaction(suffix, evidence, record_type, categories, location, name, transformation,
                  unresolved):
    r = new_reaction(f"{DOC}:{suffix}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = src(categories, categories, evidence)
    r["location"] = src(location, location, evidence)
    r["reaction_name"] = src(name, name, evidence)
    r["transformation_type"] = inferred(transformation, evidence, "Transformation description was inferred from the displayed starting material/product connectivity and caption.", raw=name)
    r["review"] = {
        "status": "needs_review",
        "method": "model",
        "reviewer": None,
        "checked_fields": [],
        "unresolved": list(unresolved),
    }
    r["completeness"] = "partial"
    return r


def set_solvent(r, raw, normalized, evidence):
    r["solvent_conditions"]["composition"] = src(raw, normalized, evidence)
    r["solvent_conditions"]["mixture_ratio"] = missing(evidence)
    r["solvent_conditions"]["total_volume"] = missing(evidence)
    r["solvent_conditions"]["concentration"] = missing(evidence)
    r["solvent_conditions"]["concentration_basis"] = missing(evidence)


def set_procedure(r, raw, normalized, evidence, general=None, modifications=None, addition=None):
    r["procedure"]["original_text"] = src(raw, normalized, evidence)
    r["procedure"]["general_procedure_reference"] = src(general, general, evidence) if general else missing(evidence)
    r["procedure"]["modifications"] = src(modifications, modifications, evidence) if modifications else missing(evidence)
    r["procedure"]["addition_order"] = src(addition, addition, evidence) if addition else missing(evidence)
    r["procedure"]["pretreatment"] = missing(evidence)
    for key in r["special_conditions"]:
        r["special_conditions"][key] = not_applicable(evidence)


def outcome(oid, kind, targets, raw_value, normalized, evidence, product_role=None,
            yield_type=None, basis=None, scope=None):
    o = new_outcome(oid, kind, targets, evidence)
    o["product_role"] = src(product_role, product_role, evidence) if product_role else not_applicable(evidence)
    o["value"] = src(raw_value, normalized, evidence)
    o["yield_type"] = src(yield_type, yield_type, evidence) if yield_type else not_applicable(evidence)
    o["measurement_basis"] = src(basis, basis, evidence) if basis else missing(evidence)
    o["scope"] = src(scope, scope, evidence) if scope else missing(evidence)
    return o


def change(change_id, elements, before, after, kind, evidence, description, hydrogen_count=None):
    value = {
        "change_id": change_id,
        "atoms": [{"atom_map": None, "element": elements[0]}, {"atom_map": None, "element": elements[1]}],
        "before_order": before,
        "after_order": after,
        "change_type": kind,
        "description": description,
    }
    if hydrogen_count is not None:
        value["hydrogen_count"] = hydrogen_count
    return derived(value, evidence, "Manual net-connectivity comparison of source drawings; atom identities are not fully mapped.", raw=description)


def set_classification(r, changes, evidence, coverage="partial", limitations=()):
    r["classification"]["standard_id"] = STANDARD
    r["classification"]["mapped_reaction_smiles"] = missing(evidence)
    r["classification"]["bond_changes"] = changes
    labels = {}
    for c in changes:
        value = c["normalized_value"]
        if value["change_type"] in ("formed", "order_increased", "hydrogen_gain"):
            label = bond_label(value["atoms"][0]["element"], value["atoms"][1]["element"], value["after_order"])
            labels.setdefault(label["code"], {**label, "supporting_change_ids": []})["supporting_change_ids"].append(value["change_id"])
    r["classification"]["labels"] = derived(list(labels.values()), evidence, "Labels are the taxonomy projection of the manually listed bond changes.", raw="manual bond-change review")
    r["classification"]["coverage"] = inferred(coverage, evidence, "Coverage is partial because a complete atom mapping, coproduct accounting, and implicit-H provenance are unavailable.", raw="source drawings visually inspected")
    r["classification"]["limitations"] = list(limitations)
    r["classification"]["model_votes"] = []


def set_no_classification(r, evidence, limitation, coverage="not_applicable"):
    r["classification"]["mapped_reaction_smiles"] = missing(evidence)
    r["classification"]["bond_changes"] = []
    r["classification"]["labels"] = derived([], evidence, "No observed transformation is classified." if coverage == "not_applicable" else "No reliable net bond mapping was asserted.", raw=limitation)
    r["classification"]["coverage"] = inferred(coverage, evidence, "No observed product connectivity is available for formed-bond classification.", raw=limitation)
    r["classification"]["limitations"] = [limitation]
    r["classification"]["model_votes"] = []


dataset = json.loads(DATASET_PATH.read_text())
doc = dataset["documents"][0]
doc["title"] = "Synthesis of pyrrole-ferrocene ensembles and their rearrangement into 2-(ferrocenylmethyl)-1,2-dihydro-3H-pyrrol-3-ones"
doc["doi"] = "10.1016/j.jorganchem.2020.121651"
write(DATASET_PATH, dataset)

bundle = empty(dataset["documents"])
molecules = {}


def register(*args, **kwargs):
    m = make_molecule(*args, **kwargs)
    molecules[m["molecule_id"]] = m
    return m["molecule_id"]


# Reagents and solvents with source-derived candidate connectivity.
register("m_allene1", "1", "(ferrocenylmethoxy)allene", ev(2), chirality="achiral", geometry="not_applicable", representation="not_encoded",
         axial=[{"fragment": "C=C=C", "assignment": "not stereogenic", "basis": "one terminal carbon is CH2"}],
         note="The terminal allene is not axially stereogenic, but a complete organometallic SMILES is not asserted.")
small = [
    ("m_nBuLi", "n-BuLi", "n-butyllithium", "[Li]CCCC", 2),
    ("m_thf", "THF", "tetrahydrofuran", "C1CCOC1", 2),
    ("m_hexane", "hexane", "hexane", "CCCCCC", 2),
    ("m_water", "H2O", "water", "O", 2),
    ("m_me_ncs", "MeNCS", "methyl isothiocyanate", "CN=C=S", 2),
    ("m_ipr_ncs", "i-PrNCS", "isopropyl isothiocyanate", "CC(C)N=C=S", 2),
    ("m_voe_ncs", "CH2=CHOCH2CH2NCS", "2-(vinyloxy)ethyl isothiocyanate", "C=COCCN=C=S", 2),
    ("m_cpent_ncs", "cyclo-C5H9NCS", "cyclopentyl isothiocyanate", "S=C=NC1CCCC1", 2),
    ("m_chept_ncs", "cyclo-C7H13NCS", "cycloheptyl isothiocyanate", "S=C=NC1CCCCCC1", 2),
    ("m_ph_ncs", "PhNCS", "phenyl isothiocyanate", "S=C=Nc1ccccc1", 2),
    ("m_mei", "MeI", "methyl iodide", "CI", 2),
    ("m_eti", "EtI", "ethyl iodide", "CCI", 3),
    ("m_propargyl_br", "HC≡CCH2Br", "propargyl bromide", "C#CCBr", 2),
    ("m_benzyl_br", "BnBr", "benzyl bromide", "BrCc1ccccc1", 2),
    ("m_methyl_bromoacetate", "MeOC(O)CH2Br", "methyl bromoacetate", "COC(=O)CBr", 2),
    ("m_tbutyl_bromoacetate", "t-BuOC(O)CH2Br", "tert-butyl bromoacetate", "CC(C)(C)OC(=O)CBr", 2),
    ("m_bromoacetonitrile", "BrCH2CN", "bromoacetonitrile", "N#CCBr", 2),
    ("m_cubr", "CuBr", "copper(I) bromide", "[Cu+].[Br-]", 2),
    ("m_cdcl3", "CDCl3", "deuterated chloroform", "[2H]C(Cl)(Cl)Cl", 4),
    ("m_toluene", "toluene", "toluene", "Cc1ccccc1", 5),
    ("m_allyloxyallene", "13a", "(allyloxy)allene", "C=C=COCC=C", 4),
    ("m_benzyloxyallene", "13b", "(benzyloxy)allene", "C=C=COCc1ccccc1", 4),
    ("m_tbuok", "t-BuOK", "potassium tert-butoxide", "[K+].CC(C)(C)[O-]", 6),
    ("m_dmso", "DMSO", "dimethyl sulfoxide", "CS(C)=O", 6),
]
for mid, label, name, smiles, page in small:
    register(mid, label, name, ev(page), smiles=smiles, chirality="achiral", geometry="not_applicable")
register("m_cux", "CuX", "copper(I) halide (CuI or CuBr)", ev(2), chirality="not_applicable", geometry="not_applicable", representation="not_encoded",
         note="The page-2 prose permits CuI or CuBr and the page-3 generic arrow states CuX; row-specific identity is not assigned from these pages.")

product2_names = {
    "a": "3-(ferrocenylmethoxy)-1-methyl-2-(methylsulfanyl)-1H-pyrrole",
    "b": "3-(ferrocenylmethoxy)-1-isopropyl-2-(methylsulfanyl)-1H-pyrrole",
    "c": "3-(ferrocenylmethoxy)-2-(methylsulfanyl)-1-[2-(vinyloxy)ethyl]-1H-pyrrole",
    "d": "1-cyclopentyl-3-(ferrocenylmethoxy)-2-(methylsulfanyl)-1H-pyrrole",
    "e": "1-cycloheptyl-3-(ferrocenylmethoxy)-2-(methylsulfanyl)-1H-pyrrole",
    "f": "3-(ferrocenylmethoxy)-2-(methylsulfanyl)-1-phenyl-1H-pyrrole",
    "g": "2-(ethylsulfanyl)-3-(ferrocenylmethoxy)-1-[2-(vinyloxy)ethyl]-1H-pyrrole",
    "h": "3-(ferrocenylmethoxy)-1-methyl-2-(prop-2-ynylsulfanyl)-1H-pyrrole",
    "i": "2-(benzylsulfanyl)-3-(ferrocenylmethoxy)-1-methyl-1H-pyrrole",
    "j": "methyl {[3-(ferrocenylmethoxy)-1-methyl-1H-pyrrol-2-yl]sulfanyl}acetate",
    "k": "tert-butyl {[3-(ferrocenylmethoxy)-1-methyl-1H-pyrrol-2-yl]sulfanyl}acetate",
    "l": "{[3-(ferrocenylmethoxy)-1-methyl-1H-pyrrol-2-yl]sulfanyl}acetonitrile",
}
for key, name in product2_names.items():
    register(f"m_2{key}", f"2{key}", name, ev(3), chirality="unspecified", geometry="not_applicable", representation="not_encoded",
             note="The ferrocenyl product drawing was inventoried, but no complete machine SMILES is asserted.")

for key in "abcdefghijk":
    register(f"m_7{key}", f"7{key}", f"2-(ferrocenylmethyl)-2-sulfanyl-1,2-dihydro-3H-pyrrol-3-one 7{key}", ev(4),
             chirality="unspecified", geometry="not_applicable", representation="not_encoded",
             atom_centers=[{"center": "C2 of pyrrol-3-one", "assignment": "unspecified", "source_marking": "no wedge/dash"}],
             note="C2 is potentially stereogenic; the table supplies no absolute or relative assignment.")

register("m_12a", "12a", "3-(allyloxy)-2-(methylsulfanyl)-1-phenyl-1H-pyrrole", ev(4), chirality="unspecified", geometry="not_applicable", representation="not_encoded")
register("m_12b", "12b", "3-(benzyloxy)-2-(methylsulfanyl)-1-phenyl-1H-pyrrole", ev(4), chirality="unspecified", geometry="not_applicable", representation="not_encoded")
register("m_17", "17", "3-(ferrocenylmethoxy)-2-(methylsulfanyl)-1-vinyl-1H-pyrrole", ev(6), chirality="unspecified", geometry="not_applicable", representation="not_encoded")
register("m_18", "18", "2-(ferrocenylmethyl)-2-(methylsulfanyl)-1-vinyl-1,2-dihydro-3H-pyrrol-3-one", ev(6), chirality="unspecified", geometry="not_applicable", representation="not_encoded",
         atom_centers=[{"center": "C2 of pyrrol-3-one", "assignment": "unspecified", "source_marking": "no wedge/dash"}])

bundle["molecules"] = list(molecules.values())
reactions = []

# Generic routes shown on pages 2--3.
r = base_reaction("scheme_table1_onepot", ev(2) + ev(3), "generic_scheme", ["other"],
                  "Page 3, Table 1 reaction arrow; page 2 conditions", "one-pot construction of pyrroles 2",
                  "lithiation/addition/S-alkylation/Cu(I)-catalyzed cyclization",
                  ["R1, R2, catalyst halide, and exact charges vary by row; the generic arrow is not a concrete experiment."])
r["participants"] = [participant("p_allene", "m_allene1", "reactant", ev(2) + ev(3)), participant("p_nBuLi", "m_nBuLi", "reagent", ev(2) + ev(3)), participant("p_thf", "m_thf", "solvent", ev(2) + ev(3)), participant("p_hexane", "m_hexane", "solvent", ev(2) + ev(3)), participant("p_cux", "m_cux", "catalyst", ev(2) + ev(3))]
set_solvent(r, "THF/hexane", ["THF", "hexane"], ev(2) + ev(3))
r["steps"] = [
    step("s1_low_temperature_lithiation", ev(2), ["p_allene", "p_nBuLi", "p_thf", "p_hexane"], "lithiation", ("−100 °C", quantity(-100, "°C", "−100 °C"))),
    step("s2_lithiation_hold", ev(2), ["p_allene", "p_nBuLi"], "hold after lithiation", ("−65 to −40 °C", qrange(-65, -40, "°C", "−65 to −40 °C")), ("5–15 min", qrange(5, 15, "min", "5–15 min"))),
    step("s3_S_alkylation", ev(2), ["p_allene"], "S-alkylation after isothiocyanate addition", ("−80 to −5 °C", qrange(-80, -5, "°C", "−80 to −10 to −5 °C"))),
    step("s4_cyclization", ev(2), ["p_cux"], "Cu(I)-catalyzed intramolecular cyclization", ("30–50 °C", qrange(30, 50, "°C", "30–50 °C")), ("15–50 min", qrange(15, 50, "min", "15–50 min"))),
]
set_procedure(r, {"figure": "Table 1 generic arrow", "sequence": ["n-BuLi", "R1NCS", "R2X", "CuX"]}, "one-pot four-stage sequence", ev(2) + ev(3), addition="n-BuLi, isothiocyanate, alkylating agent, then CuX")
set_classification(r, [change("b1", ("C", "S"), 0, 1, "formed", ev(3), "alkylating-agent carbon–sulfur bond"), change("b2", ("C", "N"), 0, 1, "formed", ev(3), "ring-closing carbon–nitrogen bond")], ev(3), limitations=["R-group atoms and implicit-H changes are not completely mapped."])
reactions.append(r)

r = base_reaction("scheme1_two_step_switch", ev(2), "generic_scheme", ["other"], "Page 2, Scheme 1 lower route", "two-step switching of triene cyclization from thiophene to pyrrole", "isolation of 1-aza-1,3,4-triene followed by CuBr cyclization", ["The scheme is generic in R1 and ester R3; the equivalent and catalyst denominators are not explicitly restated."])
p_allene = participant("p_allene", "m_allene1", "reactant", ev(2)); p_water = participant("p_water", "m_water", "additive", ev(2)); add_equiv(p_water, 1, "1 equivalent", "m_allene1", ev(2)); p_cubr = participant("p_cubr", "m_cubr", "catalyst", ev(2)); add_mol_percent(p_cubr, 9, "~9.0 mol%", "m_allene1", ev(2), approximate=True)
r["participants"] = [p_allene, p_water, p_cubr, participant("p_thf", "m_thf", "solvent", ev(2))]
set_solvent(r, "THF", ["THF"], ev(2))
r["steps"] = [step("s1_water_addition", ev(2), ["p_allene", "p_water", "p_thf"], "add water before S-alkylation", ("−95 to −80 °C", qrange(-95, -80, "°C", "−95 to −80 °C"))), step("s2_low_temperature_alkylation", ev(2), ["p_allene", "p_water"], "add alkyl bromoacetate", ("−80 to −30 °C", qrange(-80, -30, "°C", "−80 to −30 °C")), ("10 min", quantity(10, "min", "10 min"))), step("s3_isolate_triene", ev(2), ["p_allene"], "isolate 1-aza-1,3,4-triene 5", stage="workup"), step("s4_cubr_cyclization", ev(2), ["p_allene", "p_cubr", "p_thf"], "cyclize isolated triene with CuBr", ("r.t.", "r.t."), ("0.5 h", quantity(0.5, "h", "0.5 h")))]
set_procedure(r, {"figure": "Scheme 1 lower route"}, "water-assisted low-temperature alkylation, isolation, then CuBr cyclization", ev(2), modifications="one equivalent water and lower alkylation temperature suppress thiophene formation", addition="water, alkyl bromoacetate, isolate triene, then CuBr")
set_classification(r, [change("b1", ("C", "S"), 0, 1, "formed", ev(2), "alkylating-agent carbon–sulfur bond"), change("b2", ("C", "N"), 0, 1, "formed", ev(2), "ring-closing carbon–nitrogen bond")], ev(2), limitations=["Generic R groups and the isolated intermediate are not completely mapped."])
reactions.append(r)

# Table 1: one concrete substrate-scope record per row.
table1 = {
    "a": ("m_me_ncs", "m_mei", 55, "onepot"), "b": ("m_ipr_ncs", "m_mei", 44, "onepot"),
    "c": ("m_voe_ncs", "m_mei", 67, "onepot"), "d": ("m_cpent_ncs", "m_mei", 22, "onepot"),
    "e": ("m_chept_ncs", "m_mei", 37, "onepot"), "f": ("m_ph_ncs", "m_mei", 41, "onepot"),
    "g": ("m_voe_ncs", "m_eti", 45, "onepot"), "h": ("m_me_ncs", "m_propargyl_br", 18, "onepot"),
    "i": ("m_me_ncs", "m_benzyl_br", 50, "onepot"), "j": ("m_me_ncs", "m_methyl_bromoacetate", 31, "twostep"),
    "k": ("m_me_ncs", "m_tbutyl_bromoacetate", 20, "twostep"), "l": ("m_me_ncs", "m_bromoacetonitrile", 3, "twostep"),
}
ncs_conditions = {
    "m_me_ncs": (("−50 to −40 °C", qrange(-50, -40, "°C", "−50 to −40 °C")), ("5–10 min", qrange(5, 10, "min", "5–10 min"))),
    "m_voe_ncs": (("−45 to −40 °C", qrange(-45, -40, "°C", "−45 to −40 °C")), ("10 min", quantity(10, "min", "10 min"))),
    "m_ph_ncs": (("ca. −65 °C", quantity(-65, "°C", "ca. −65 °C", "approx")), ("15 min", quantity(15, "min", "15 min"))),
    "m_ipr_ncs": (("ca. −30 °C", quantity(-30, "°C", "ca. −30 °C", "approx")), ("20 min", quantity(20, "min", "20 min"))),
    "m_cpent_ncs": (("ca. −30 °C", quantity(-30, "°C", "ca. −30 °C", "approx")), ("20 min", quantity(20, "min", "20 min"))),
    "m_chept_ncs": (("ca. −30 °C", quantity(-30, "°C", "ca. −30 °C", "approx")), ("20 min", quantity(20, "min", "20 min"))),
}
for key, (ncs, alkylator, yld, route) in table1.items():
    evidence = ev(2) + ev(3) + ev(6)
    r = base_reaction(f"table1_2{key}", evidence, "specific_experiment", ["substrate_scope"], f"Page 3, Table 1, row 2{key}", "Synthesis of 3-(ferrocenylmethoxy)-2-sulfanyl-1H-pyrrole", "multicomponent pyrrole construction", ["Exact charges, volumes, concentration, and row-specific CuI/CuBr identity are not stated on pages 1–6.", "Complete ferrocene-containing structures and atom mapping remain unresolved."])
    r["general_scheme_ids"] = [f"{DOC}:scheme_table1_onepot" if route == "onepot" else f"{DOC}:scheme1_two_step_switch"]
    parts = [participant("p_allene", "m_allene1", "reactant", evidence), participant("p_ncs", ncs, "reactant", evidence), participant("p_alkylator", alkylator, "reagent", evidence), participant("p_nBuLi", "m_nBuLi", "reagent", evidence), participant("p_thf", "m_thf", "solvent", evidence), participant("p_hexane", "m_hexane", "solvent", evidence), participant("p_product", f"m_2{key}", "product", evidence, "main product")]
    if route == "onepot":
        parts.append(participant("p_cux", "m_cux", "catalyst", evidence))
    else:
        p_w = participant("p_water", "m_water", "additive", evidence); add_equiv(p_w, 1, "1 equivalent", "m_allene1", ev(2)); parts.append(p_w)
        p_c = participant("p_cubr", "m_cubr", "catalyst", evidence); add_mol_percent(p_c, 9, "~9.0 mol%", "m_allene1", ev(2), approximate=True); parts.append(p_c)
    r["participants"] = parts
    set_solvent(r, "THF/hexane; THF is the reaction solvent and hexane is the n-BuLi carrier", ["THF", "hexane"], evidence)
    temp, tm = ncs_conditions[ncs]
    r["steps"] = [
        step("s1_low_temperature_lithiation", ev(2), ["p_allene", "p_nBuLi", "p_thf", "p_hexane"], "add n-BuLi to allene 1", ("−100 °C", quantity(-100, "°C", "−100 °C")), atmosphere="argon"),
        step("s2_lithiation_hold", ev(2), ["p_allene", "p_nBuLi"], "hold lithiated allene", ("−65 to −40 °C", qrange(-65, -40, "°C", "−65 to −40 °C")), ("5–15 min", qrange(5, 15, "min", "5–15 min")), atmosphere="argon"),
        step("s3_isothiocyanate_addition", ev(2), ["p_allene", "p_ncs"], "add isothiocyanate to lithioallene", temp, tm, atmosphere="argon"),
    ]
    if route == "onepot":
        r["steps"] += [step("s4_S_alkylation", ev(2), ["p_allene", "p_ncs", "p_alkylator"], "S-alkylation", ("−80 to −5 °C", qrange(-80, -5, "°C", "−80 to −10 to −5 °C")), atmosphere="argon"), step("s5_cu_cyclization", ev(2), ["p_cux", "p_product"], "Cu(I)-catalyzed intramolecular cyclization", ("30–50 °C", qrange(30, 50, "°C", "30–50 °C")), ("15–50 min", qrange(15, 50, "min", "15–50 min")), atmosphere="argon"), step("s6_purification", ev(2), ["p_product"], "column chromatography on neutral alumina", stage="purification")]
    else:
        r["steps"] += [step("s4_water_addition", ev(2), ["p_water", "p_allene"], "add water before S-alkylation", ("−95 to −80 °C", qrange(-95, -80, "°C", "−95 to −80 °C")), atmosphere="argon"), step("s5_low_temperature_alkylation", ev(2), ["p_alkylator", "p_allene"], "S-alkylation under modified conditions", ("−80 to −30 °C", qrange(-80, -30, "°C", "−80 to −30 °C")), ("10 min", quantity(10, "min", "10 min")), atmosphere="argon"), step("s6_isolate_triene", ev(2), ["p_allene"], "isolate the 1-aza-1,3,4-triene", stage="workup"), step("s7_cubr_cyclization", ev(2), ["p_cubr", "p_product", "p_thf"], "CuBr-catalyzed cyclization of isolated triene", ("r.t.", "r.t."), ("0.5 h", quantity(0.5, "h", "0.5 h")), atmosphere="argon"), step("s8_purification", ev(2), ["p_product"], "column chromatography on neutral alumina", stage="purification")]
    r["reaction_scale"] = missing(evidence)
    set_procedure(r, {"table": "Table 1", "row": f"2{key}", "footnote": "two preparative steps according to Scheme 1" if route == "twostep" else None}, "general lithiation–isothiocyanate addition–S-alkylation–Cu(I) cyclization sequence", evidence, general="Table 1 generic sequence and page-2 prose", modifications="two preparative steps with water and isolated triene" if route == "twostep" else None, addition="n-BuLi, isothiocyanate, alkylating agent, then Cu(I) catalyst")
    basis = "isolated purified product; exact purification basis is not specified on pages 1–6" if key == "l" else "isolated product after neutral-alumina chromatography"
    r["outcomes"] = [outcome("o_yield", "yield", [f"m_2{key}"], f"{yld}%", quantity(yld, "%", f"{yld}%"), ev(3) + ev(6), "main product", "isolated yield", basis, {"span": "allene 1/isothiocyanate to pyrrole 2", "table_row": f"2{key}"})]
    set_classification(r, [change("b1", ("C", "S"), 0, 1, "formed", ev(3), "S-alkyl bond in product"), change("b2", ("C", "N"), 0, 1, "formed", ev(3), "ring-closing C–N bond in pyrrole")], ev(3), limitations=["Complete atom-to-atom mapping and implicit-H redistribution are unresolved."])
    reactions.append(r)

# Scheme 3 synthesis of the two comparator pyrroles and their negative rearrangement controls.
for key, start in (("a", "m_allyloxyallene"), ("b", "m_benzyloxyallene")):
    evidence = ev(4)
    r = base_reaction(f"scheme3_synthesis_12{key}", evidence, "specific_experiment", ["other"], f"Page 4, Scheme 3, route to 12{key}", "Synthesis of comparator pyrrole 12", "lithiation/addition/S-alkylation/CuBr cyclization", ["Scheme 3 does not report charges, time, temperature, or yield; later experimental pages are outside the bounded coverage."])
    r["general_scheme_ids"] = [f"{DOC}:scheme_table1_onepot"]
    r["participants"] = [participant("p_allene", start, "reactant", evidence), participant("p_phncs", "m_ph_ncs", "reactant", evidence), participant("p_mei", "m_mei", "reagent", evidence), participant("p_nBuLi", "m_nBuLi", "reagent", evidence), participant("p_cubr", "m_cubr", "catalyst", evidence), participant("p_thf", "m_thf", "solvent", evidence), participant("p_hexane", "m_hexane", "solvent", evidence), participant("p_product", f"m_12{key}", "product", evidence, "main comparator product")]
    set_solvent(r, "THF/hexane", ["THF", "hexane"], evidence)
    r["steps"] = [step("s1_sequence", evidence, [p["participant_id"] for p in r["participants"]], "n-BuLi, PhNCS, MeI, then CuBr; detailed conditions not shown")]
    set_procedure(r, {"figure": "Scheme 3", "reagents": ["n-BuLi", "PhNCS", "MeI", "CuBr"]}, "displayed comparator-pyrrole synthesis", evidence, addition="n-BuLi, PhNCS, MeI, CuBr")
    set_classification(r, [change("b1", ("C", "S"), 0, 1, "formed", evidence, "methyl–sulfur bond"), change("b2", ("C", "N"), 0, 1, "formed", evidence, "ring-closing carbon–nitrogen bond")], evidence, limitations=["Complete atom mapping and implicit-H changes are unresolved."])
    reactions.append(r)

    r = base_reaction(f"scheme3_no_rearrangement_12{key}", evidence, "specific_experiment", ["control", "mechanistic_study"], f"Page 4, Scheme 3, CDCl3 control for 12{key}", "No rearrangement of comparator pyrrole 12", "negative rearrangement control", ["The quantity recovered and analytical detection limit are not reported; expected products are not entered as observed products."])
    r["participants"] = [participant("p_start", f"m_12{key}", "reactant", evidence), participant("p_cdcl3", "m_cdcl3", "solvent", evidence)]
    set_solvent(r, "CDCl3", ["CDCl3"], evidence)
    r["steps"] = [step("s1_storage", evidence, ["p_start", "p_cdcl3"], "store solution and assess rearrangement", ("r.t.", "r.t."), ("~20 days", quantity(20, "d", "~20 days", "approx")))]
    set_procedure(r, {"figure": "Scheme 3", "crossed_out_product": f"14{key}"}, "CDCl3 room-temperature no-rearrangement control", evidence)
    r["outcomes"] = [outcome("o_no_reaction", "no_reaction", [f"m_12{key}"], "did not rearrange; starting material recovered", "no rearrangement observed; starting material recovered", evidence, scope={"expected_product": f"14{key}", "duration": "~20 days"})]
    set_no_classification(r, evidence, "No observed transformation; expected rearrangement products were not identified.")
    reactions.append(r)

# Table 2: one monitoring experiment per starting pyrrole, preserving all time points.
table2 = {
    "a": [("12 d", 12, "d", 82), ("19 d", 19, "d", 96)],
    "b": [("12 d", 12, "d", 17)],
    "c": [("20 d", 20, "d", 77), ("55 d", 55, "d", 100)],
    "d": [("13 d", 13, "d", 2), ("28 d", 28, "d", 19)],
    "e": [("21 d", 21, "d", 21)],
    "f": [("2 d", 2, "d", 24), ("8 d", 8, "d", 96), ("19 d", 19, "d", 100)],
    "g": [("during purification on neutral Al2O3", None, None, 12)],
    "h": [("15 d", 15, "d", 42)],
    "i": [("16 d", 16, "d", 100)],
    "j": [("10 d", 10, "d", 16), ("4.5 months", 4.5, "months", 90)],
    "k": [("4 months", 4, "months", 100)],
    "l": [("14 d", 14, "d", 0)],
}
for key, points in table2.items():
    evidence = ev(4) + (ev(5) if key in "acefgijl" else [])
    unresolved = ["The concentration and sample amount are not reported on pages 1–6.", "Ferrocene-containing structures lack a complete checked atom mapping."]
    if key in "jk":
        unresolved.append("Footnote c states that the actual time required for the reported 90/100% conversion is undefined.")
    if key == "l":
        unresolved.append("No product 7l is entered because the source reports 0% conversion/no tendency to rearrange.")
    r = base_reaction(f"table2_rearrangement_2{key}", evidence, "specific_experiment", ["mechanistic_study"], f"Page 4, Table 2, row 2{key}→7{key}", "Rearrangement of pyrrole 2 into pyrrol-3-one 7", "formal [1,3]-O-to-C rearrangement", unresolved)
    r["participants"] = [participant("p_start", f"m_2{key}", "reactant", evidence), participant("p_cdcl3", "m_cdcl3", "solvent", evidence)]
    if key != "l":
        r["participants"].append(participant("p_product", f"m_7{key}", "product", evidence, "rearrangement product"))
    set_solvent(r, "CDCl3", ["CDCl3"], evidence)
    r["steps"] = []
    r["outcomes"] = []
    for idx, (raw_time, number, unit, conv) in enumerate(points, 1):
        pids = [p["participant_id"] for p in r["participants"]]
        if number is None:
            r["steps"].append(step(f"s{idx}_purification", ev(4), pids, "conversion occurred during purification on neutral Al2O3", stage="purification"))
        elif unit == "d":
            r["steps"].append(step(f"s{idx}_monitor", ev(4), pids, "store CDCl3 solution and record 1H NMR spectrum", ("r.t.", "r.t."), (raw_time, quantity(number, "d", raw_time))))
        else:
            month_step = step(f"s{idx}_monitor", ev(4), pids, "store CDCl3 solution and record 1H NMR spectrum", ("r.t.", "r.t."))
            month_step["time"] = fact(raw=raw_time, normalized=None, evidence=ev(4), method="original_extraction", status="pending_review", verification="unreviewed", derivation="The exact source month unit is preserved without inventing a day conversion; months are not a supported normalized quantity unit.")
            r["steps"].append(month_step)
        r["outcomes"].append(outcome(f"o_conversion_{idx}", "conversion", [f"m_2{key}"], f"{conv}%", quantity(conv, "%", f"{conv}%"), ev(4), basis="1H NMR", scope={"starting_material": f"2{key}", "product": None if key == "l" else f"7{key}", "time": raw_time}))
    if key == "l":
        r["outcomes"].append(outcome("o_no_reaction", "no_reaction", ["m_2l"], "did not show a tendency to rearrange", "no rearrangement observed", ev(4) + ev(5), scope={"duration": "14 d", "solvent": "CDCl3"}))
    set_procedure(r, {"table": "Table 2", "row": f"2{key}", "footnotes": ["conversion by 1H NMR"] + (["actual time undefined"] if key in "jk" else [])}, "CDCl3 storage/monitoring experiment", evidence)
    if key == "l":
        set_no_classification(r, evidence, "0% conversion/no observed rearrangement.")
    else:
        set_classification(r, [change("b1", ("C", "O"), 1, 0, "broken", ev(4), "ferrocenylmethyl C–O bond cleavage"), change("b2", ("C", "C"), 0, 1, "formed", ev(4), "ferrocenylmethyl–pyrrole C2 bond"), change("b3", ("C", "O"), 1, 2, "order_increased", ev(4), "pyrrole C3–O promoted to carbonyl")], ev(4), limitations=["Formal net changes only; the paper presents alternative acid-mediated and sigmatropic mechanisms."])
    reactions.append(r)

# Page 5 proposed mechanisms are records, not experimental proof.
for suffix, location, name, mechanism in (
    ("scheme5_proposed_sigmatropic", "Page 5, Scheme 5", "Plausible formal C[1,3] shift mechanism", "tandem O[1,3]-sigmatropic and [3,3]-sigmatropic shifts via A–C"),
    ("scheme6_proposed_ionic", "Page 5, Scheme 6", "Plausible CDCl3 rearrangement mechanism", "C–O heterolysis to ferrocenylmethyl cation D and pyrrolyloxide/carbanion E/F"),
):
    r = base_reaction(suffix, ev(5), "proposed_mechanism", ["other"], location, name, mechanism, ["This is an author-proposed pathway, not a separately observed experiment; intermediate structures are not registered as isolated molecules."])
    set_solvent(r, "CDCl3" if "ionic" in suffix else "not specified for the generic proposal", ["CDCl3"] if "ionic" in suffix else None, ev(5))
    set_procedure(r, {"figure": location}, mechanism, ev(5))
    set_no_classification(r, ev(5), "Proposed elementary pathway; no complete mapped net reaction is asserted.", coverage="unresolved")
    reactions.append(r)

# Thermal rearrangement of 2f in toluene.
r = base_reaction("thermal_rearrangement_2f", ev(5), "specific_experiment", ["mechanistic_study"], "Page 5, thermal rearrangement paragraph", "Thermally induced rearrangement of 2f to 7f", "formal [1,3]-O-to-C rearrangement under heating", ["Sample amount, toluene volume, concentration, and isolation at the monitoring points are not reported on page 5."])
r["participants"] = [participant("p_start", "m_2f", "reactant", ev(5)), participant("p_toluene", "m_toluene", "solvent", ev(5)), participant("p_product", "m_7f", "product", ev(5), "rearrangement product")]
set_solvent(r, "toluene", ["toluene"], ev(5))
r["steps"] = [step("s1_3h", ev(5), ["p_start", "p_toluene", "p_product"], "heat and monitor by 1H NMR", ("98–100 °C", qrange(98, 100, "°C", "98–100 °C")), ("3 h", quantity(3, "h", "3 h"))), step("s2_8h", ev(5), ["p_start", "p_toluene", "p_product"], "continue heating and monitor by 1H NMR", ("98–100 °C", qrange(98, 100, "°C", "98–100 °C")), ("8 h", quantity(8, "h", "8 h")))]
r["outcomes"] = [outcome("o_conversion_3h", "conversion", ["m_2f"], "45%", quantity(45, "%", "45%"), ev(5), basis="1H NMR", scope={"time": "3 h", "product": "7f"}), outcome("o_conversion_8h", "conversion", ["m_2f"], "96%", quantity(96, "%", "96%"), ev(5), basis="1H NMR", scope={"time": "8 h", "product": "7f"})]
set_procedure(r, "a solution of pyrrole 2f in toluene was heated at 98–100 °C and monitored by 1H NMR", "thermal time-course monitoring", ev(5))
set_classification(r, [change("b1", ("C", "O"), 1, 0, "broken", ev(5), "ferrocenylmethyl C–O bond cleavage"), change("b2", ("C", "C"), 0, 1, "formed", ev(5), "ferrocenylmethyl–pyrrole C2 bond"), change("b3", ("C", "O"), 1, 2, "order_increased", ev(5), "formation of pyrrol-3-one carbonyl")], ev(5), limitations=["Formal net changes; a complete atom mapping is not supplied."])
reactions.append(r)

# Page 6 Scheme 7 derivatization.
r = base_reaction("scheme7_2c_to_17", ev(6), "specific_experiment", ["derivatization"], "Page 6, Scheme 7 and section 2.3", "Transformation of 2c into N-vinylpyrrole 17", "base-induced cleavage/elimination to an N-vinylpyrrole", ["Scheme 7 depicts loss as CH2=CHOK whereas the prose calls it elimination of vinyl alcohol; the coproduct identity is therefore not asserted.", "Charges and DMSO volume are not reported on page 6."])
r["participants"] = [participant("p_start", "m_2c", "reactant", ev(6)), participant("p_tbuok", "m_tbuok", "reagent", ev(6)), participant("p_dmso", "m_dmso", "solvent", ev(6)), participant("p_product", "m_17", "product", ev(6), "main product")]
set_solvent(r, "DMSO", ["DMSO"], ev(6))
r["steps"] = [step("s1_elimination", ev(6), ["p_start", "p_tbuok", "p_dmso", "p_product"], "treat 2c with t-BuOK/DMSO", ("r.t.", "r.t."), ("~30 min", quantity(30, "min", "~30 min", "approx")))]
r["outcomes"] = [outcome("o_yield", "yield", ["m_17"], "96%", quantity(96, "%", "96%"), ev(6), "main product", "isolated yield", "isolated product", {"starting_material": "2c", "product": "17"})]
set_procedure(r, {"figure": "Scheme 7", "reagents": ["t-BuOK", "DMSO"]}, "room-temperature base-induced formation of 17", ev(6))
set_classification(r, [change("b1", ("C", "O"), 1, 0, "broken", ev(6), "cleavage of the side-chain C–O bond"), change("b2", ("C", "C"), 1, 2, "order_increased", ev(6), "formation of the N-vinyl C=C bond"), change("b3", ("C", "H"), 1, 0, "hydrogen_loss", ev(6), "net C–H loss accompanying elimination", hydrogen_count=1)], ev(6), limitations=["Coproduct identity/proton transfer and complete hydrogen accounting remain unresolved."])
reactions.append(r)

# Page 6 Scheme 8 rearrangement of 17.
r = base_reaction("scheme8_17_to_18", ev(6), "specific_experiment", ["mechanistic_study"], "Page 6, Scheme 8", "Rearrangement of N-vinylpyrrole 17 into pyrrol-3-one 18", "formal [1,3]-O-to-C rearrangement", ["The scheme gives 25% but does not explicitly distinguish isolated yield from conversion or state the measurement method.", "Temperature, sample amount, solvent volume, and concentration are unreported in the displayed region."])
r["participants"] = [participant("p_start", "m_17", "reactant", ev(6)), participant("p_cdcl3", "m_cdcl3", "solvent", ev(6)), participant("p_product", "m_18", "product", ev(6), "rearrangement product")]
set_solvent(r, "CDCl3", ["CDCl3"], ev(6))
r["steps"] = [step("s1_rearrangement", ev(6), ["p_start", "p_cdcl3", "p_product"], "store CDCl3 solution and observe rearrangement", None, ("10 days", quantity(10, "d", "10 days")))]
r["outcomes"] = [outcome("o_reported_percentage", "yield", ["m_18"], "25%", quantity(25, "%", "25%"), ev(6), "rearrangement product", "unresolved (scheme percentage)", None, {"starting_material": "17", "product": "18", "reported_as": "product label followed by 25%"})]
set_procedure(r, {"figure": "Scheme 8", "condition": "CDCl3, 10 days"}, "CDCl3 rearrangement of 17", ev(6))
set_classification(r, [change("b1", ("C", "O"), 1, 0, "broken", ev(6), "ferrocenylmethyl C–O bond cleavage"), change("b2", ("C", "C"), 0, 1, "formed", ev(6), "ferrocenylmethyl–pyrrole C2 bond"), change("b3", ("C", "O"), 1, 2, "order_increased", ev(6), "formation of pyrrol-3-one carbonyl")], ev(6), limitations=["Formal net changes only; complete atom mapping is unavailable."])
reactions.append(r)

bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": "ctx_page1_bibliography", "type": "bibliographic_and_abstract", "evidence_ids": ev(1), "summary": "Title, authors, DOI, abstract, and introduction; no concrete experiment row or reaction arrow."},
    {"context_id": "ctx_page2_figure1", "type": "characterization", "evidence_ids": ev(2), "summary": "Figure 1 is the crystal structure of 2f with bond metrics; it is not a synthesis experiment."},
    {"context_id": "ctx_page3_xrd", "type": "characterization", "evidence_ids": ev(3), "summary": "Crystal packing and crystallographic data for 2f; no separate synthesis record."},
    {"context_id": "ctx_page4_scheme2", "type": "literature_and_hypothetical_context", "evidence_ids": ev(4), "summary": "Scheme 2 shows classical Claisen rearrangements and hypothetical products 10/11; these were not observed in this paper."},
    {"context_id": "ctx_page5_scheme4", "type": "negative_hypothetical_context", "evidence_ids": ev(5), "summary": "Expected Claisen products 14a, 15b, and 16b are crossed out/not identified; negative controls are represented by the Scheme 3 control records."},
    {"context_id": "ctx_page5_mass_spectrum", "type": "characterization", "evidence_ids": ev(5), "summary": "m/z 199 [FcCH2]+ is mechanistic support/characterization, not a synthetic experiment."},
    {"context_id": "ctx_page6_general_methods", "type": "general_methods", "evidence_ids": ev(6), "summary": "General solvent drying, argon, TLC, instrumentation, and statement that yields refer to isolated products."},
    {"context_id": "ctx_page6_xray", "type": "characterization", "evidence_ids": ev(6), "summary": "Section 4.2 X-ray crystallography details for 2f; no reaction record."},
]

write(SUBMITTED_PATH, bundle)

inventory = {
    "document_id": DOC,
    "pages_reviewed": [1, 2, 3, 4, 5, 6],
    "visual_review": "original-resolution JPEG page render inspected for every page",
    "pages": {
        "1": ["title/abstract/introduction", "no concrete reaction arrow or table"],
        "2": ["generic one-pot conditions", "Scheme 1 two-path switch", "Figure 1 XRD characterization"],
        "3": ["Table 1 generic arrow", "12 substrate-scope rows", "XRD prose"],
        "4": ["Scheme 2 literature/hypothetical context", "Scheme 3 two syntheses and two controls", "Table 2 twelve monitoring rows"],
        "5": ["Scheme 4 expected but unobserved products", "Schemes 5 and 6 proposed mechanisms", "thermal 2f experiment", "MS characterization"],
        "6": ["Schemes 7 and 8 concrete transformations", "general methods", "X-ray methods"],
    },
}
write(RUN / "raw_model" / "page_inventory_0001_0006.json", inventory)

diagnostics = {}
for mid, m in molecules.items():
    smiles = m["smiles"]["normalized_value"]
    if smiles:
        diagnostics[mid] = stereo_inspect(smiles)
write(RUN / "raw_model" / "stereo_diagnostics.json", diagnostics)

record_ids = [r["reaction_id"] for r in reactions]
def ids(prefixes):
    return [rid for rid in record_ids if any(rid.split(":", 1)[1].startswith(p) for p in prefixes)]

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [{"region": "full page: title, abstract, introduction", "finding": "Bibliographic context and high-level synthesis/rearrangement summary; no concrete experimental row, table, or reaction arrow.", "evidence_ids": ev(1)}], "record_ids": [], "unresolved": ["Abstract-level conditions were not promoted into concrete experiments."]},
        {"page": 2, "state": "attempted", "regions": [{"region": "upper and middle columns: section 2.1 prose", "finding": "Generic lithiation, isothiocyanate-addition, alkylation, and Cu(I)-cyclization condition ranges used by Table 1 records.", "evidence_ids": ev(2)}, {"region": "lower-left Scheme 1", "finding": "One-pot thiophene side route and water-assisted two-step pyrrole route; lower route registered as a generic scheme.", "evidence_ids": ev(2)}, {"region": "lower-right Figure 1", "finding": "X-ray structure/metrics for 2f; characterization only.", "evidence_ids": ev(2)}], "record_ids": ids(["scheme_table1", "scheme1_", "table1_"]), "unresolved": ["Row-specific CuI/CuBr identity and exact charges are not specified on this page.", "Scheme 1 upper thiophene products 6 are generic/prior-condition side outcomes and were not converted into concrete product records."]},
        {"page": 3, "state": "attempted", "regions": [{"region": "upper two-thirds: Table 1 reaction diagram and rows 2a–2l", "finding": "Twelve substrate-scope experiments with isolated yields 55, 44, 67, 22, 37, 41, 45, 18, 50, 31, 20, and 3%; rows 2j–2l carry the two-step footnote.", "evidence_ids": ev(3)}, {"region": "lower prose", "finding": "Crystallographic discussion of 2f; characterization, not a separate synthesis.", "evidence_ids": ev(3)}], "record_ids": ids(["scheme_table1", "table1_"]), "unresolved": ["Complete organometallic structures and atom mapping remain unencoded."]},
        {"page": 4, "state": "attempted", "regions": [{"region": "top-left Scheme 2", "finding": "Classical Claisen literature examples and hypothetical 2→10/11 products; not observed experiments.", "evidence_ids": ev(4)}, {"region": "top-right Scheme 3", "finding": "Two comparator-pyrrole syntheses (12a,12b) and two CDCl3 no-rearrangement controls (~20 d).", "evidence_ids": ev(4)}, {"region": "bottom Table 2", "finding": "Twelve substrate-specific CDCl3 rearrangement monitoring experiments; multiple time points retained within the same sample record; conversions measured by 1H NMR.", "evidence_ids": ev(4)}], "record_ids": ids(["scheme3_", "table2_"]), "unresolved": ["Actual time for late 2j/2k conversions is explicitly undefined by footnote c.", "2g conversion occurs during alumina purification, with no numeric duration."]},
        {"page": 5, "state": "attempted", "regions": [{"region": "top Scheme 4", "finding": "Expected Claisen products are crossed out/not identified; retained as negative context, not product records.", "evidence_ids": ev(5)}, {"region": "middle Schemes 5 and 6", "finding": "Two alternative author-proposed mechanisms, separately classified as proposed mechanisms.", "evidence_ids": ev(5)}, {"region": "lower-right prose", "finding": "Concrete thermal 2f→7f time course in toluene at 98–100 °C: 45% conversion at 3 h and 96% at 8 h by 1H NMR.", "evidence_ids": ev(5)}, {"region": "middle prose: m/z 199", "finding": "Mass-spectral fragment used as mechanistic support; characterization only.", "evidence_ids": ev(5)}], "record_ids": ids(["scheme5_", "scheme6_", "thermal_"]) + [f"{DOC}:table2_rearrangement_2{x}" for x in "acefgijl"], "unresolved": ["The proposed mechanisms are not treated as experimentally proven."]},
        {"page": 6, "state": "attempted", "regions": [{"region": "top Schemes 7 and 8", "finding": "Concrete 2c→17 base-induced derivatization (96%) and 17→18 CDCl3 rearrangement (10 d, 25%).", "evidence_ids": ev(6)}, {"region": "section 4.1 General", "finding": "General inert-atmosphere, solvent, TLC, and isolated-yield conventions used as inherited context where applicable.", "evidence_ids": ev(6)}, {"region": "section 4.2", "finding": "X-ray crystallography method for 2f; characterization only.", "evidence_ids": ev(6)}], "record_ids": ids(["scheme7_", "scheme8_"]) + ids(["table1_"]), "unresolved": ["Scheme 8 does not state whether 25% is isolated yield or conversion and gives no measurement basis.", "Scheme 7 prose/byproduct depiction does not uniquely resolve the eliminated coproduct identity."]},
    ],
}
write(RUN / "coverage" / "pages_0001_0006.json", coverage)

# When rerun after the official merge, preserve the advisory audit as a raw
# machine artifact.  It is deliberately not treated as chemical verification.
merged_dataset = json.loads(DATASET_PATH.read_text())
write(RUN / "raw_model" / "audit_after_merge.json", dataset_audit(merged_dataset))

print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]), "source_context": len(bundle["source_context"])}, ensure_ascii=False))
