#!/usr/bin/env python3
"""Build evidence-linked v3.1 records for original PDF pages 7--10.

This source-specific extraction artifact updates the matching page-3 scope
records with the exact experimental paragraphs on pages 7--8 and adds the
distinct 4b scale-up experiment.  It does not modify schemas or validators.
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

DOC = "doc_2edd574f8525f831"
BATCH = "allene-stereo-v31-200-20261017"
OUT = ROOT / "runs/current"
E = {p: f"{DOC}:page{p}:image" for p in range(1, 11)}
T = {p: f"{DOC}:page{p}:text" for p in range(1, 11)}


def evidence(pages, kinds=("image", "text")):
    ids = []
    for p in pages:
        if "image" in kinds:
            ids.append(E[p])
        if "text" in kinds:
            ids.append(T[p])
    return ids


def F(raw=None, norm=None, pages=(), method="original_extraction", derivation=None,
      status=None, alternatives=(), kinds=("image", "text")):
    return st.fact(raw, norm, evidence(pages, kinds), method=method,
                   derivation=derivation, status=status,
                   alternatives=alternatives)


def missing(pages=(), status="unreported"):
    return F(pages=pages, status=status)


def q(number, unit, raw=None, qualifier=None):
    value = {"number": number, "unit": unit}
    if raw is not None:
        value["raw_text"] = raw
    if qualifier is not None:
        value["qualifier"] = qualifier
    return value


def qrange(low, high, unit, raw):
    return {"minimum": low, "maximum": high, "unit": unit, "raw_text": raw}


def rawq(value):
    if isinstance(value, dict):
        if value.get("raw_text"):
            return value["raw_text"]
        if "number" in value:
            return f"{value['number']:g} {value['unit']}"
        if "minimum" in value:
            return f"{value['minimum']:g}–{value['maximum']:g} {value['unit']}"
    return str(value)


template = st.read(OUT / "raw_model/pages_0007_0010_template.json")
existing = st.read(OUT / "standard_records.json")
existing_molecules = {m["molecule_id"]: m for m in existing["molecules"]}
molecules = {}


def reuse_molecule(mid, label, pages, *, chirality=None, geometry=None,
                   representation=None, assignment=None):
    """Reuse a registry entity, adding this source occurrence and stereo review."""
    m = copy.deepcopy(existing_molecules[mid])
    ev = evidence(pages)
    if not any(o.get("document_id") == DOC and o.get("evidence_ids") == ev
               and st.value(o.get("compound_label", {})) == label
               for o in m["occurrences"]):
        m["occurrences"].append({
            "document_id": DOC,
            "compound_label": F(label, label, pages),
            "evidence_ids": ev,
        })
    smi = st.value(m["smiles"])
    if smi:
        diag = st.stereo_inspect(smi)
        m["stereochemistry"]["isomeric_smiles"] = F(
            None, diag["canonical_isomeric_smiles"], pages,
            method="structure_derivation",
            derivation="Standard stereo graph diagnostic on the retained candidate structure; not source certification.")
        m["stereochemistry"]["identity_key"] = F(
            None, diag["stereo_graph_key"], pages,
            method="structure_derivation",
            derivation="Standard stereo graph key on the retained candidate structure; not source certification.")
    if chirality is not None:
        m["stereochemistry"]["chirality_status"] = F(
            None, chirality, pages, method="inference",
            derivation="Center/axis review of the source drawing and candidate graph; no unreported racemic state was assigned.")
    if geometry is not None:
        m["stereochemistry"]["geometry_status"] = F(
            None, geometry, pages, method="inference",
            derivation="E/Z applicability and explicit source labels were reviewed for this occurrence.")
    if representation is not None:
        m["stereochemistry"]["representation_status"] = F(
            None, representation, pages, method="inference",
            derivation="Representation status reflects only stereochemistry explicitly supported by the source occurrence.")
    if assignment is not None:
        m["stereochemistry"]["assignment"] = F(assignment, assignment, pages)
    molecules[mid] = m
    return mid


def unresolved_candidate(mid, label, name, pages, candidate_smiles, note,
                         geometry="unresolved", assignment=None):
    ev = evidence(pages)
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "unresolved"
    m["occurrences"][0]["compound_label"] = F(label, label, pages)
    m["names"]["original"] = F(name, name, pages)
    m["names"]["english"] = F(name, name, pages)
    m["cas_numbers"] = missing(pages)
    m["representation_type"] = "unresolved"
    alternatives = [
        F(text, smi, pp, method="structure_derivation",
          derivation="Machine candidate retained from the indicated source strand; conflicting strands prevent registry selection.")
        for text, smi, pp in candidate_smiles
    ]
    m["smiles"] = F(note, None, pages, status="conflict", alternatives=alternatives)
    m["identity_note"] = F(note, note, pages)
    m["stereochemistry"]["chirality_status"] = F(
        None, "achiral", pages, method="inference",
        derivation="Both retained connectivity candidates lack a stereogenic center or qualifying allene/axis.")
    m["stereochemistry"]["geometry_status"] = F(
        None, geometry, pages, method="inference",
        derivation="The source depicts/names an E product even though full connectivity is unresolved.")
    m["stereochemistry"]["representation_status"] = F(
        None, "unresolved", pages, method="inference",
        derivation="No single connectivity candidate was selected, so stereochemistry is not asserted as encoded.")
    if assignment:
        m["stereochemistry"]["assignment"] = F(assignment, assignment, pages)
        m["stereochemistry"]["double_bonds"] = F(
            "(1E) in source", [{"configuration": "E", "source_locator": "product alkene"}], pages)
    molecules[mid] = m
    return mid


# Registry entities used by this batch. Terminal allenes have no qualifying
# axial chirality and no E/Z geometry; source omission is not called racemic.
for mid, label, pp in [
    ("m_eosinY", "EY", (7, 8)), ("m_k2co3", "K2CO3", (7,)),
    ("m_dmf", "DMF", (7, 8)), ("m_tfe", "TFE", (7, 8)),
    ("m_ptsna", "p-TolSO2Na", (7, 8)), ("m_phsna", "PhSO2Na", (7, 8)),
    ("m_brsna", "4-BrC6H4SO2Na", (7, 8)), ("m_thsna", "2-thienylSO2Na", (7,)),
    ("m_mesna", "MeSO2Na", (7,)),
]:
    reuse_molecule(mid, label, pp)

for mid, label, pp in [
    ("m_1e", "1e", (7,)), ("m_1b", "1b", (7,)), ("m_1d", "1d", (7,)),
    ("m_3a", "3a", (7, 8)), ("m_3b", "3b", (7,)), ("m_3c", "3c", (7,)),
    ("m_3d", "3d", (8,)), ("m_3f", "3f", (8,)), ("m_3h", "3h", (8,)),
]:
    reuse_molecule(mid, label, pp, chirality="achiral",
                   geometry="not_applicable", representation="not_applicable")

for mid, label, pp in [
    ("m_2k", "2k", (3, 7)), ("m_2l", "2l", (3, 7)), ("m_2m", "2m", (3, 7)),
    ("m_4a", "4a", (3, 7)), ("m_4b", "4b", (3, 7, 8)),
    ("m_4c", "4c", (3, 7)), ("m_4e", "4e", (3, 7)),
    ("m_4f", "4f", (3, 7)), ("m_4g", "4g", (3, 7)),
    ("m_4i", "4i", (3, 8)), ("m_4l", "4l", (3, 8)),
    ("m_4m", "4m", (3, 8)), ("m_4n", "4n", (3, 8)),
]:
    reuse_molecule(mid, label, pp, chirality="achiral", geometry="specified",
                   representation="encoded",
                   assignment="Source structure/name specifies the E product alkene; no center, axis, or ring cis/trans stereochemistry is present.")

reuse_molecule("m_2lprime", "2l′", (3, 7), chirality="achiral",
               geometry="not_applicable", representation="not_applicable",
               assignment="Terminal alkene regioisomer; E/Z is not applicable.")

ts_4h = "C/C(=C\\c1ccc2ccccc2c1)S(=O)(=O)c1ccc(C)cc1"
ph_4h = "C/C(=C\\c1ccc2ccccc2c1)S(=O)(=O)c1ccccc1"
unresolved_candidate(
    "m_4h", "4h", "conflicting source identity for (1E)-naphthyl vinyl sulfone 4h",
    (3, 7, 8),
    [("Scheme 3 structure and page-8 formula/NMR support a p-tolylsulfonyl product", ts_4h, (3, 8)),
     ("Page-7 heading and charged sodium benzenesulfinate support a phenylsulfonyl product", ph_4h, (7,))],
    "4h connectivity conflict: Scheme 3 draws Ts and page 8 reports C20H19O2S plus an extra methyl signal, whereas the page-7 heading says phenyl sulfone and charges sodium benzenesulfinate.",
    geometry="specified", assignment="The product alkene is specified as (1E) in the source heading/drawing; sulfonyl aryl identity remains unresolved.")

unresolved_candidate(
    "m_4h_sulfinate_conflict", "sulfinate for 4h",
    "sodium benzenesulfinate / sodium 4-methylbenzenesulfinate source conflict",
    (3, 7),
    [("Page-7 exact charge names sodium benzenesulfinate", "[Na+].[O-]S(=O)c1ccccc1", (7,)),
     ("Scheme-3 product drawing requires a p-tolylsulfonyl source", "[Na+].[O-]S(=O)c1ccc(C)cc1", (3,))],
    "The exact charge (0.222 g, 1.3554 mmol) is printed as sodium benzenesulfinate, but Scheme 3 draws product 4h as the Ts derivative; no identity is selected.",
    geometry="not_applicable")


def base_reaction(suffix, pages, *, record_type="specific_experiment", categories=(),
                  name, location):
    r = st.new_reaction(f"{DOC}:{suffix}", DOC, evidence(pages))
    r["record_type"] = record_type
    r["experiment_categories"] = F(str(list(categories)), list(categories), pages)
    r["location"] = F(location, {"pages": list(pages), "region": location}, pages)
    r["reaction_name"] = F(name, name, pages)
    r["transformation_type"] = F(
        None, "hydrosulfonylation", pages, method="inference",
        derivation="Transformation type assigned from the source allene-to-vinyl-sulfone reaction.")
    r["review"] = {
        "status": "needs_review", "method": "model", "reviewer": None,
        "checked_fields": [],
        "unresolved": ["Machine extraction requires expert verification of structures, source conflicts, and completeness."],
    }
    r["completeness"] = "partial"
    return r


def participant(r, pid, mid, role, pages, *, mass=None, amount=None, volume=None,
                equiv=None, molpct=None, ref=None, stage="reaction"):
    p = st.new_participant(pid, mid, role, evidence(pages), stage)
    for key, value in (("mass", mass), ("amount", amount), ("volume", volume),
                       ("equivalents", equiv), ("mol_percent", molpct)):
        if value is not None:
            p["charges"][key] = F(rawq(value), value, pages)
    if equiv is not None or molpct is not None:
        p["equivalence_reference"] = F(ref, ref, pages)
    r["participants"].append(p)
    return p


def step(r, sid, pages, operation, participant_ids=(), *, temperature=None,
         time=None, atmosphere=None, stage="reaction", method="original_extraction",
         derivation=None):
    s = st.new_step(sid, evidence(pages), stage)
    s["participant_ids"] = list(participant_ids)
    s["operation"] = F(operation, operation, pages, method=method, derivation=derivation)
    if temperature is not None:
        s["temperature"] = F(rawq(temperature), temperature, pages, method=method, derivation=derivation)
    if time is not None:
        s["time"] = F(rawq(time), time, pages, method=method, derivation=derivation)
    if atmosphere is not None:
        s["atmosphere"] = F(str(atmosphere), atmosphere, pages, method=method, derivation=derivation)
    r["steps"].append(s)
    return s


def set_solvent(r, pages, dmf, tfe):
    r["solvent_conditions"]["composition"] = F("DMF–TFE", ["DMF", "TFE"], pages)
    r["solvent_conditions"]["mixture_ratio"] = F("1:0.3", "1:0.3", pages)
    r["solvent_conditions"]["total_volume"] = F(
        f"{dmf + tfe:g} mL", q(dmf + tfe, "mL"), pages,
        method="inference", derivation="Arithmetic sum of separately reported DMF and TFE volumes.")


def set_light(r, pages, power):
    r["special_conditions"]["light_source"] = F("green LED", "green LED", pages)
    r["special_conditions"]["wavelength"] = F("530 nm", q(530, "nm"), pages)
    r["special_conditions"]["power"] = F(power, power, pages)


def product_outcome(r, pages, target, mass, percent, *, scheme_percent=None,
                    combined=False):
    mass_out = st.new_outcome("product_mass", "product_amount", [target], evidence(pages))
    mass_out["value"] = F(f"{mass:g} g", q(mass, "g"), pages)
    mass_out["measurement_basis"] = F("reported isolated product mass", "reported isolated product mass", pages)
    r["outcomes"].append(mass_out)

    y = st.new_outcome("yield", "yield", [target], evidence(pages))
    y["yield_type"] = F("combined isolated" if combined else "isolated",
                         "combined isolated" if combined else "isolated", pages)
    y["measurement_basis"] = F(
        "isolated after the stated chromatography/purification",
        "isolated after the stated chromatography/purification", pages)
    if scheme_percent is not None and scheme_percent != percent:
        y["value"] = F(
            f"source conflict: scope scheme {scheme_percent}% versus product paragraph {percent}%",
            None, (3,) + tuple(pages), status="conflict",
            alternatives=[
                F(f"scope scheme: {scheme_percent}%", q(scheme_percent, "%"), (3,)),
                F(f"product paragraph: {percent}%", q(percent, "%"), pages),
            ])
    else:
        y["value"] = F(f"{percent}%", q(percent, "%"), pages)
    r["outcomes"].append(y)
    return y


def classify(r, pages):
    changes = [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "S"}],
         "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}],
         "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
        {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
         "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
    ]
    r["classification"]["bond_changes"] = [
        F(None, change, pages, method="structure_derivation",
          derivation="Partial net inventory from the source allene and vinyl-sulfone structures; no complete atom mapping asserted.")
        for change in changes
    ]
    r["classification"]["labels"] = F(
        None,
        [{"category_number": 4, "code": "C-S", "supporting_change_ids": ["b1"]},
         {"category_number": 1, "code": "C-H", "supporting_change_ids": ["b2"]}],
        pages, method="structure_derivation",
        derivation="Labels cover actual new C-S and implicit-H C-H bonds; the unchanged product alkene and bond-order decrease are not called newly formed bonds.")
    r["classification"]["coverage"] = F(
        None, "partial", pages, method="inference",
        derivation="C-S formation, net C-H gain, and one allene bond-order decrease are represented; atom mapping and coproducts remain incomplete.")
    r["classification"]["limitations"] = [
        "No complete atom mapping.", "Proton source is not assigned at the experiment-record level.",
        "Unlisted spectators/coproducts are excluded.",
    ]


def common_conditions(r, family, participant_ids):
    if family == "indolyl":
        step(r, "irradiation", (3, 5),
             "under nitrogen, irradiate with green LED and stir; monitor by TLC",
             participant_ids, temperature=q(25, "°C"),
             time=qrange(3, 5, "h", "3–5 h"), atmosphere="N2",
             method="inherited_general_condition",
             derivation="Inherited from Scheme 2 and the page-5 general procedure; the individual page-7 paragraph does not restate time/temperature.")
        product_ids = tuple(p["participant_id"] for p in r["participants"] if p["role"] == "product")
        step(r, "purification", (5,),
             "evaporate solvent; silica chromatography with EtOAc–hexane; crystallize in hexane with a few drops of EtOAc and wash",
             product_ids, stage="purification", method="inherited_general_condition",
             derivation="Inherited from the page-5 general procedure.")
        set_light(r, (3,), "one 30 W LED or two 15 W 530 nm LEDs")
    else:
        step(r, "irradiation", (3, 7),
             "under nitrogen, irradiate with green LED and stir; monitor by TLC",
             participant_ids, temperature=q(25, "°C"),
             time=qrange(2, 4, "h", "2–4 h"), atmosphere="N2",
             method="inherited_general_condition",
             derivation="Inherited from Scheme 3 and the page-7 general procedure; the individual product paragraph does not restate time/temperature.")
        step(r, "evaporation", (7,), "evaporate solvent under reduced pressure",
             ("product",), stage="workup", method="inherited_general_condition",
             derivation="Inherited from the page-7 arylallene general procedure.")
        step(r, "chromatography", (7,),
             "silica-gel chromatography with EtOAc–hexane mixtures",
             ("product",), stage="purification", method="inherited_general_condition",
             derivation="Inherited from the page-7 arylallene general procedure.")
        set_light(r, (3,), "one 30 W LED or two 15 W 530 nm LEDs")


reactions = []

# A general procedure is retained as a generic record and not counted as an
# additional product experiment.
r = base_reaction(
    "general_procedure_aryl_hydrosulfonylation", (3, 7),
    record_type="generic_scheme", name="General Procedure for the Hydrosulfonylation of Arylallenes",
    location="page 7, General Procedure for the Hydrosulfonylation of Arylallenes; Scheme 3 conditions on page 3")
r["experiment_categories"] = F("[]", [], (3, 7))
r["solvent_conditions"]["composition"] = F("DMF–TFE", ["DMF", "TFE"], (7,))
r["solvent_conditions"]["mixture_ratio"] = F("1:0.3", "1:0.3", (7,))
step(r, "irradiation", (3, 7), "mix allene, 5 mol% eosin Y, and 3 equiv sulfinate in DMF–TFE under nitrogen; irradiate with green LED and stir with TLC tracking",
     temperature=q(25, "°C"), time=qrange(2, 4, "h", "2–4 h"), atmosphere="N2")
step(r, "evaporation", (7,), "evaporate solvent under reduced pressure", stage="workup")
step(r, "chromatography", (7,), "silica-gel chromatography with EtOAc–hexane mixtures", stage="purification")
set_light(r, (3,), "one 30 W LED or two 15 W 530 nm LEDs")
classify(r, (3, 7))
reactions.append(r)


indolyl_details = {
    "2k": dict(sub="1e", sf="m_thsna", sm=(0.100, 0.5405), ey=(0.019, 0.027),
               sfq=(0.275, 1.620), k=(0.149, 0.1081), dmf=4.513, tfe=1.354,
               pm=0.045, py=25, scheme=25),
    "2m": dict(sub="1d", sf="m_phsna", sm=(0.500, 4.760), ey=(0.165, 0.238),
               sfq=(2.342, 14.280), k=(1.314, 9.520), dmf=20.0, tfe=6.0,
               pm=0.410, py=35, scheme=27),
}
for prod, d in indolyl_details.items():
    pages = (3, 5, 7)
    r = base_reaction(
        f"scope_{prod}", pages, categories=("substrate_scope",),
        name="Photoredox-catalyzed hydrosulfonylation of N-allenyl heteroarene",
        location=f"page 3 Scheme 2 product {prod}; page 7 exact experimental paragraph")
    participant(r, "allene", f"m_{d['sub']}", "reactant", (7,), mass=q(d["sm"][0], "g"),
                amount=q(d["sm"][1], "mmol"), equiv=q(1, "equiv"), ref=f"m_{d['sub']}")
    participant(r, "EY", "m_eosinY", "catalyst", (7,), mass=q(d["ey"][0], "g"),
                amount=q(d["ey"][1], "mmol"), molpct=q(5, "mol%"), ref=f"m_{d['sub']}")
    participant(r, "sulfinate", d["sf"], "reagent", (7,), mass=q(d["sfq"][0], "g"),
                amount=q(d["sfq"][1], "mmol"), equiv=q(3, "equiv"), ref=f"m_{d['sub']}")
    participant(r, "K2CO3", "m_k2co3", "additive", (7,), mass=q(d["k"][0], "g"),
                amount=q(d["k"][1], "mmol"), equiv=q(2, "equiv"), ref=f"m_{d['sub']}")
    participant(r, "DMF", "m_dmf", "solvent", (7,), volume=q(d["dmf"], "mL"))
    participant(r, "TFE", "m_tfe", "solvent", (7,), volume=q(d["tfe"], "mL"))
    participant(r, "product", f"m_{prod}", "product", (3, 7))
    set_solvent(r, (7,), d["dmf"], d["tfe"])
    common_conditions(r, "indolyl", tuple(p["participant_id"] for p in r["participants"]))
    r["reaction_scale"] = F(
        f"{d['sub']} ({d['sm'][1]:g} mmol)",
        {"reference_molecule_id": f"m_{d['sub']}", "amount": q(d["sm"][1], "mmol")}, (7,))
    r["procedure"]["general_procedure_reference"] = F(
        "General Procedure for the Hydrosulfonylation of N-Indolylallenes",
        "page 5 general procedure", (5, 7))
    product_outcome(r, (7,), f"m_{prod}", d["pm"], d["py"], scheme_percent=d["scheme"])
    if d["scheme"] != d["py"]:
        r["review"]["unresolved"].append(
            f"Yield conflict for {prod}: Scheme 2 gives {d['scheme']}%, while the page-7 experimental paragraph gives {d['py']}%; neither was selected.")
    if prod == "2k":
        r["notes"].append(F(
            "Printed K2CO3 charge is 0.149 g (0.1081 mmol, 2 equiv); the stated mmol is inconsistent with both mass and two equivalents relative to 0.5405 mmol allene.",
            "K2CO3 source quantity conflict retained without correction", (7,)))
        r["review"]["unresolved"].append(
            "2k K2CO3 charge is internally inconsistent (0.149 g; 0.1081 mmol; 2 equiv); printed values are preserved.")
    if prod == "2m":
        r["notes"].append(F(
            "Scheme 2 labels the drawn product 2m at 27%, the adjacent prose calls it 2n at 27%, and the page-7 experimental heading labels it 2m at 35%.",
            "2m/2n label and 27%/35% yield conflict", (3, 7)))
        r["review"]["unresolved"].append(
            "The paper inconsistently uses 2m/2n for the N-allenylpyrrole product; registry retains drawn/experimental label 2m, with prose label 2n preserved in notes.")
    classify(r, pages)
    reactions.append(r)


# Combined regioisomer experiment 2l + 2l′.
r = base_reaction(
    "scope_2l_2lprime", (3, 5, 7), categories=("substrate_scope",),
    name="Photoredox-catalyzed hydrosulfonylation giving regioisomers 2l and 2l′",
    location="page 3 Scheme 2 combined entry; page 7 exact experimental paragraph")
participant(r, "allene", "m_1b", "reactant", (7,), mass=q(0.087, "g"), amount=q(0.3551, "mmol"), equiv=q(1, "equiv"), ref="m_1b")
participant(r, "EY", "m_eosinY", "catalyst", (7,), mass=q(0.012, "g"), amount=q(0.0178, "mmol"), molpct=q(5, "mol%"), ref="m_1b")
participant(r, "sulfinate", "m_ptsna", "reagent", (7,), mass=q(0.190, "g"), amount=q(1.0653, "mmol"), equiv=q(3, "equiv"), ref="m_1b")
participant(r, "K2CO3", "m_k2co3", "additive", (7,), mass=q(0.098, "g"), amount=q(0.7102, "mmol"), equiv=q(2, "equiv"), ref="m_1b")
participant(r, "DMF", "m_dmf", "solvent", (7,), volume=q(2.966, "mL"))
participant(r, "TFE", "m_tfe", "solvent", (7,), volume=q(0.890, "mL"))
participant(r, "major", "m_2l", "product", (3, 7))
participant(r, "minor", "m_2lprime", "product", (3, 7))
set_solvent(r, (7,), 2.966, 0.890)
common_conditions(r, "indolyl", tuple(p["participant_id"] for p in r["participants"]))
r["reaction_scale"] = F("1b (0.3551 mmol)", {"reference_molecule_id": "m_1b", "amount": q(0.3551, "mmol")}, (7,))
r["procedure"]["general_procedure_reference"] = F("General Procedure for the Hydrosulfonylation of N-Indolylallenes", "page 5 general procedure", (5, 7))
mass_o = st.new_outcome("product_mass", "product_amount", ["m_2l", "m_2lprime"], evidence((7,)))
mass_o["value"] = F("0.075 g", q(0.075, "g"), (7,))
mass_o["measurement_basis"] = F("reported mass of isolated regioisomer mixture", "reported mass of isolated regioisomer mixture", (7,))
r["outcomes"].append(mass_o)
y = st.new_outcome("combined_yield", "yield", ["m_2l", "m_2lprime"], evidence((3, 7)))
y["value"] = F("source conflict: Scheme 2 gives 58%; experimental paragraph gives 53%", None, (3, 7), status="conflict", alternatives=[
    F("Scheme 2: 58%", q(58, "%"), (3,)), F("page-7 product paragraph: 53%", q(53, "%"), (7,))])
y["yield_type"] = F("combined isolated", "combined isolated", (3, 7))
y["measurement_basis"] = F("isolated regioisomer mixture after the stated purification", "isolated regioisomer mixture after the stated purification", (7,))
y["scope"] = F("1b to 2l + 2l′", {"from": "m_1b", "to": ["m_2l", "m_2lprime"]}, (3, 7))
r["outcomes"].append(y)
ratio = st.new_outcome("regioisomer_ratio", "selectivity", ["m_2l", "m_2lprime"], evidence((3,)))
ratio["value"] = F("2l:2l′ = 2.45:1", {"metric": "regioisomer_ratio", "ratio": [2.45, 1], "order": ["m_2l", "m_2lprime"]}, (3,))
ratio["scope"] = F("order: 2l, 2l′", {"metric": "regioisomer ratio", "order": ["m_2l", "m_2lprime"]}, (3,))
# Measurement method is intentionally left null/missing: the scheme reports the ratio but not its assay.
r["outcomes"].append(ratio)
r["review"]["unresolved"].extend([
    "Combined yield conflict remains unresolved: Scheme 2 gives 58%, experimental paragraph gives 53%.",
    "The 2l:2l′ = 2.45:1 ratio is explicit, but the measurement method is not reported.",
])
classify(r, (3, 5, 7))
reactions.append(r)


aryl_details = {
    "4a": dict(sub="3a", sf="m_ptsna", sm=(0.088, 0.759), ey=(0.026, 0.0380), sfq=(0.405, 2.276), dmf=6.34, tfe=1.9, pm=0.141, py=68, scheme=74),
    "4b": dict(sub="3a", sf="m_phsna", sm=(0.090, 0.776), ey=(0.027, 0.0388), sfq=(0.394, 2.328), dmf=6.48, tfe=1.94, pm=0.156, py=78, scheme=78),
    "4c": dict(sub="3a", sf="m_brsna", sm=(0.070, 0.6034), ey=(0.021, 0.03017), sfq=(0.440, 1.810), dmf=5.039, tfe=1.311, pm=0.150, py=74, scheme=74),
    "4e": dict(sub="3a", sf="m_thsna", sm=(0.062, 0.539), ey=(0.019, 0.02695), sfq=(0.275, 1.6164), dmf=4.5, tfe=1.35, pm=0.069, py=48, scheme=48),
    "4f": dict(sub="3a", sf="m_mesna", sm=(0.070, 0.6034), ey=(0.021, 0.03017), sfq=(0.185, 1.81), dmf=5.039, tfe=1.511, pm=0.097, py=82, scheme=82),
    "4g": dict(sub="3c", sf="m_ptsna", sm=(0.100, 0.769), ey=(0.027, 0.03845), sfq=(0.411, 2.307), dmf=6.422, tfe=1.927, pm=0.143, py=65, scheme=65),
    "4h": dict(sub="3b", sf="m_4h_sulfinate_conflict", sm=(0.075, 0.4518), ey=(0.016, 0.02259), sfq=(0.222, 1.3554), dmf=3.77, tfe=1.13, pm=0.111, py=76, scheme=76),
    "4i": dict(sub="3f", sf="m_ptsna", sm=(0.076, 0.5757), ey=(0.020, 0.02879), sfq=(0.307, 1.7271), dmf=4.634, tfe=1.390, pm=0.028, py=17, scheme=17),
    "4l": dict(sub="3d", sf="m_phsna", sm=(0.082, 0.5616), ey=(0.019, 0.02808), sfq=(0.276, 2.307), dmf=4.690, tfe=1.407, pm=0.121, py=74, scheme=74),
    "4m": dict(sub="3h", sf="m_phsna", sm=(0.067, 0.3641), ey=(0.013, 0.0182), sfq=(0.179, 1.0923), dmf=3.050, tfe=0.900, pm=0.082, py=69, scheme=69),
    "4n": dict(sub="3h", sf="m_brsna", sm=(0.067, 0.3641), ey=(0.013, 0.0182), sfq=(0.265, 1.0923), dmf=3.050, tfe=0.900, pm=0.084, py=59, scheme=59),
}

for prod, d in aryl_details.items():
    detail_page = 7 if prod in {"4a", "4b", "4c", "4e", "4f", "4g", "4h"} else 8
    pages = (3, 7, 8) if prod == "4h" else (3, 7) if detail_page == 7 else (3, 7, 8)
    r = base_reaction(
        f"scope_{prod}", pages, categories=("substrate_scope",),
        name="Photoredox-catalyzed hydrosulfonylation of arylallene",
        location=f"page 3 Scheme 3 product {prod}; page {detail_page} exact experimental paragraph")
    participant(r, "allene", f"m_{d['sub']}", "reactant", (detail_page,), mass=q(d["sm"][0], "g"),
                amount=q(d["sm"][1], "mmol"), equiv=q(1, "equiv"), ref=f"m_{d['sub']}")
    participant(r, "EY", "m_eosinY", "catalyst", (detail_page,), mass=q(d["ey"][0], "g"),
                amount=q(d["ey"][1], "mmol"), molpct=q(5, "mol%"), ref=f"m_{d['sub']}")
    participant(r, "sulfinate", d["sf"], "reagent", (detail_page,), mass=q(d["sfq"][0], "g"),
                amount=q(d["sfq"][1], "mmol"), equiv=q(3, "equiv"), ref=f"m_{d['sub']}")
    participant(r, "DMF", "m_dmf", "solvent", (detail_page,), volume=q(d["dmf"], "mL"))
    participant(r, "TFE", "m_tfe", "solvent", (detail_page,), volume=q(d["tfe"], "mL"))
    participant(r, "product", f"m_{prod}", "product", (3, detail_page))
    set_solvent(r, (detail_page,), d["dmf"], d["tfe"])
    common_conditions(r, "aryl", tuple(p["participant_id"] for p in r["participants"]))
    r["reaction_scale"] = F(
        f"{d['sub']} ({d['sm'][1]:g} mmol)",
        {"reference_molecule_id": f"m_{d['sub']}", "amount": q(d["sm"][1], "mmol")}, (detail_page,))
    r["procedure"]["general_procedure_reference"] = F(
        "General Procedure for the Hydrosulfonylation of Arylallenes",
        "page 7 general procedure", (7, detail_page))
    product_outcome(r, (detail_page,), f"m_{prod}", d["pm"], d["py"], scheme_percent=d["scheme"])
    if d["scheme"] != d["py"]:
        r["review"]["unresolved"].append(
            f"Yield conflict for {prod}: Scheme 3 gives {d['scheme']}%, while the product paragraph gives {d['py']}%; neither was selected.")
    if prod == "4h":
        r["notes"].append(F(
            "Scheme 3 draws 4h as the Ts product, whereas the page-7 heading says phenyl sulfone and charges sodium benzenesulfinate; page-8 formula/NMR again support the Ts product.",
            "4h product/reagent identity conflict retained", (3, 7, 8)))
        r["review"]["unresolved"].append(
            "4h product and sulfinate identities remain unresolved; candidate structures are retained on unresolved molecule entities.")
    if prod == "4l":
        r["notes"].append(F(
            "Printed sodium benzenesulfinate charge is 0.276 g (2.307 mmol, 3 equiv); 2.307 mmol is inconsistent with both the mass and three equivalents relative to 0.5616 mmol allene.",
            "4l sulfinate source quantity conflict retained without correction", (8,)))
        r["review"]["unresolved"].append(
            "4l sulfinate charge is internally inconsistent (0.276 g; 2.307 mmol; 3 equiv); printed values are preserved.")
    classify(r, pages)
    reactions.append(r)


# Distinct scale-up record; it is not folded into the 0.776 mmol scope record.
r = base_reaction(
    "scaleup_4b", (3, 8), categories=("scale_up",),
    name="Scale-up hydrosulfonylation preparation of 4b",
    location="page 8 Scale-Up Experiment Procedure; page 3 Scheme 3 footnote b")
participant(r, "allene", "m_3a", "reactant", (8,), mass=q(0.232, "g"), amount=q(2, "mmol"), equiv=q(1, "equiv"), ref="m_3a")
participant(r, "EY", "m_eosinY", "catalyst", (8,), mass=q(0.070, "g"), amount=q(0.1, "mmol"), molpct=q(5, "mol%"), ref="m_3a")
participant(r, "sulfinate", "m_phsna", "reagent", (8,), mass=q(0.984, "g"), amount=q(6, "mmol"), equiv=q(3, "equiv"), ref="m_3a")
participant(r, "DMF", "m_dmf", "solvent", (8,), volume=q(16.660, "mL"))
participant(r, "TFE", "m_tfe", "solvent", (8,), volume=q(5.0, "mL"))
participant(r, "product", "m_4b", "product", (8,))
set_solvent(r, (8,), 16.660, 5.0)
step(r, "irradiation", (8,), "under nitrogen, irradiate with green LED and stir; TLC tracking",
     tuple(p["participant_id"] for p in r["participants"]), time=q(5, "h"), atmosphere="N2")
step(r, "evaporation", (8,), "evaporate solvent under reduced pressure", ("product",), stage="workup")
step(r, "chromatography", (8,), "silica-gel chromatography, EtOAc–hexane 1:5", ("product",), stage="purification")
r["special_conditions"]["light_source"] = F("green LED", "green LED", (8,))
r["reaction_scale"] = F("3a (2 mmol)", {"reference_molecule_id": "m_3a", "amount": q(2, "mmol")}, (8,))
product_outcome(r, (8,), "m_4b", 0.418, 81)
classify(r, (3, 8))
reactions.append(r)


template["molecules"] = list(molecules.values())
template["reactions"] = reactions
template["source_context"] = [
    {"context_id": f"{DOC}:page7_characterization", "kind": "characterization",
     "value": "Page 7 begins with continuation of 2j NMR/HRMS/IR and contains characterization blocks for 2k, 2l/2l′, 2m and 4a–4h; these are not standalone synthesis experiments.",
     "evidence_ids": [E[7], T[7]]},
    {"context_id": f"{DOC}:page8_characterization", "kind": "characterization",
     "value": "Page 8 begins with continuation of 4h characterization and includes characterization blocks for 4i, 4l, 4m and 4n; these are not standalone synthesis experiments.",
     "evidence_ids": [E[8], T[8]]},
    {"context_id": f"{DOC}:page8_si_availability", "kind": "availability_note",
     "value": "The article states that Supporting Information is available online, but no SI file was explicitly provided or associated in this workspace.",
     "evidence_ids": [E[8], T[8]]},
    {"context_id": f"{DOC}:pages9_10_references", "kind": "bibliography",
     "value": "Original PDF pages 9–10 contain reference-list entries only; cited work was not promoted to this paper's experiments.",
     "evidence_ids": [E[9], T[9], E[10], T[10]]},
]

st.write(OUT / "submitted.json", template)

st.write(OUT / "raw_model/pages_0007_0010_raw_extraction.json", {
    "document_id": DOC,
    "batch_id": BATCH,
    "pages": [7, 8, 9, 10],
    "semantic_accuracy_verified": False,
    "visual_review": "Original-resolution page images 7, 8, 9, and 10 were each opened and inspected.",
    "records_submitted": [r["reaction_id"] for r in reactions],
    "molecules_submitted": sorted(molecules),
    "decisions": [
        "Matched scope drawings on page 3 were merged with the exact experimental paragraphs on pages 7–8 by compound label; no old experiment was imported as a new result.",
        "The 4b 2 mmol scale-up is a distinct specific experiment, not an outcome on the 0.776 mmol scope record.",
        "Characterization-only blocks remain source_context.",
        "Pages 9–10 are bibliography only and generate no reaction records.",
        "Yield, label, quantity, and 4h identity conflicts remain explicit and unresolved.",
        "No SI association was inferred.",
    ],
})

stereo_checks = []
for m in template["molecules"]:
    smi = st.value(m["smiles"])
    if smi:
        stereo_checks.append({
            "molecule_id": m["molecule_id"], "smiles": smi,
            "standard_stereo": st.stereo_inspect(smi),
            "declared": {
                "chirality_status": st.value(m["stereochemistry"]["chirality_status"]),
                "geometry_status": st.value(m["stereochemistry"]["geometry_status"]),
                "representation_status": st.value(m["stereochemistry"]["representation_status"]),
            },
        })
st.write(OUT / "raw_model/pages_0007_0010_standard_stereo_checks.json", {
    "source_verified": False,
    "note": "Standard graph diagnostics; source drawing comparison remains needs_review. Unresolved 4h candidates are stored as alternatives and are not selected for stereo diagnostics.",
    "checks": stereo_checks,
})

mentions = st.source_mentions(
    st.read(OUT / "dataset.json"),
    ["2k", "2l", "2m", "4a", "4b", "4h", "4l", "Scale-Up Experiment Procedure"],
    context=220)
st.write(OUT / "raw_model/pages_0007_0010_source_mentions.json", mentions)

page7_ids = [f"{DOC}:{suffix}" for suffix in [
    "general_procedure_aryl_hydrosulfonylation", "scope_2k", "scope_2m",
    "scope_2l_2lprime", "scope_4a", "scope_4b", "scope_4c", "scope_4e",
    "scope_4f", "scope_4g", "scope_4h",
]]
page8_ids = [f"{DOC}:{suffix}" for suffix in [
    "scope_4h", "scope_4i", "scope_4l", "scope_4m", "scope_4n", "scaleup_4b",
]]
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 7, "state": "attempted", "regions": [
            {"region": "upper-left continuation and compounds 2k, 2l/2l′, 2m",
             "finding": "2j characterization continuation is context only; exact synthesis charges/yields for 2k, combined 2l/2l′, and 2m were merged with Scheme 2. Conflicts 2k K2CO3, 2l yield, and 2m/2n label plus yield were retained.",
             "evidence_ids": [E[7], T[7]]},
            {"region": "middle/lower columns, arylallene general procedure and compounds 4a–4h",
             "finding": "General procedure and exact synthesis paragraphs for 4a, 4b, 4c, 4e, 4f, 4g, and 4h were extracted. 4a yield and 4h reagent/product identities conflict with Scheme 3 and remain unresolved.",
             "evidence_ids": [E[7], T[7]]},
            {"region": "NMR/HRMS/IR/melting-point blocks",
             "finding": "Characterization was retained as source_context and was not converted into synthetic experiments.",
             "evidence_ids": [E[7], T[7]]},
        ], "record_ids": page7_ids,
         "unresolved": [
             "2k K2CO3 source quantity is internally inconsistent.",
             "2l/2l′ combined yield is 58% in Scheme 2 versus 53% in the experimental paragraph; ratio method for 2.45:1 is not reported.",
             "The pyrrole product is labeled 2m in scheme/experiment but 2n in prose, with 27% versus 35% yields.",
             "4a yield is 74% in Scheme 3 versus 68% in the experimental paragraph.",
             "4h sulfinate/product identity conflicts across structure, heading, charge, formula and NMR.",
         ]},
        {"page": 8, "state": "attempted", "regions": [
            {"region": "upper-left 4h continuation and compounds 4i, 4l, 4m, 4n",
             "finding": "4h characterization continuation was used only to expose the source identity conflict; exact synthesis paragraphs for 4i, 4l, 4m and 4n were extracted. The 4l sulfinate quantity conflict was retained.",
             "evidence_ids": [E[8], T[8]]},
            {"region": "left-column Scale-Up Experiment Procedure for 4b",
             "finding": "A distinct 2 mmol 4b scale-up experiment was recorded with all charges, 5 h irradiation, EtOAc–hexane 1:5 chromatography, 0.418 g product and 81% isolated yield.",
             "evidence_ids": [E[8], T[8]]},
            {"region": "Associated Content, author/acknowledgment and references",
             "finding": "SI availability was retained as context only because no SI file was provided; author material and references are not experiments.",
             "evidence_ids": [E[8], T[8]]},
        ], "record_ids": page8_ids,
         "unresolved": [
             "4l sodium benzenesulfinate charge is internally inconsistent (0.276 g, 2.307 mmol, 3 equiv).",
             "Supporting Information is mentioned but not provided/associated.",
         ]},
        {"page": 9, "state": "attempted", "regions": [
            {"region": "full page reference list",
             "finding": "Bibliography only; cited publications are prior work and were not promoted to experiments in this paper.",
             "evidence_ids": [E[9], T[9]]},
        ], "record_ids": [], "unresolved": []},
        {"page": 10, "state": "attempted", "regions": [
            {"region": "full page reference list",
             "finding": "Bibliography only; no synthesis arrow, table row, control experiment, or characterization experiment belonging to this paper is present.",
             "evidence_ids": [E[10], T[10]]},
        ], "record_ids": [], "unresolved": []},
    ],
}
st.write(OUT / "coverage/pages_0007_0010.json", coverage)

tasks = st.read(OUT / "tasks.json")
for task in tasks:
    if task.get("document_id") == DOC and task.get("evidence_id") in {T[p] for p in range(7, 11)}:
        task["state"] = "completed"
        task["reason"] = "Original page image visually attempted and mapped in coverage/pages_0007_0010.json; completion is not chemical certification."
st.write(OUT / "tasks.json", tasks)
