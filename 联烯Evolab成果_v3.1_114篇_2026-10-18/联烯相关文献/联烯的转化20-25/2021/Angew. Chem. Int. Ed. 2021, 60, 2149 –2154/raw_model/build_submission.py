#!/usr/bin/env python3
"""Build the source-specific v3.1 bundle for Zhao and Ge, Angew. 2021.

This script encodes only facts visible on PDF pages 1--6. It deliberately leaves
unavailable SI details and uncertain compact structures unresolved.
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
DOC = "doc_31953e04d69caf00"
EV = {p: f"{DOC}:page{p}:image" for p in range(1, 7)}
ET = {p: f"{DOC}:page{p}:text" for p in range(1, 7)}


def f(raw=None, normalized=None, pages=(), method="original_extraction", derivation=None,
      status=None):
    evidence = []
    for page in pages:
        evidence.extend([EV[page], ET[page]])
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation,
                   status=status)


def missing():
    return st.fact()


def q(number, unit, raw=None, qualifier=None):
    value = {"number": number, "unit": unit}
    if qualifier:
        value["qualifier"] = qualifier
    if raw:
        value["raw_text"] = raw
    return value


def derive(raw, normalized, pages, why):
    return f(raw, normalized, pages, method="structure_derivation", derivation=why)


def molecule(mid, label, pages, name=None, smiles=None, kind="ordinary",
             geometry=None, assignment=None):
    m = st.new_molecule(f"{DOC}:mol:{mid}", DOC, [EV[pages[0]], ET[pages[0]]])
    evidence = [EV[p] for p in pages] + [ET[p] for p in pages]
    m["identity_status"] = "source_identified"
    m["occurrences"] = [{"document_id": DOC,
                         "compound_label": st.fact(label, label, evidence),
                         "evidence_ids": evidence}]
    m["names"]["original"] = st.fact(name or label, name or label, evidence)
    m["names"]["english"] = st.fact(name or label, name or label, evidence)
    m["names"]["abbreviations"] = st.fact(label, [label], evidence)
    m["names"]["aliases"] = missing()
    m["cas_numbers"] = missing()
    m["external_properties"] = []
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = st.fact(f"candidate transcribed from the source drawing for {label}", smiles,
                              evidence, method="structure_derivation",
                              derivation="Candidate connectivity transcribed from the displayed structure; parsing is not identity verification.")
        m["generic_template"] = None
        m["stereochemistry"]["isomeric_smiles"] = st.fact(
            f"candidate {smiles}", smiles, evidence, method="structure_derivation",
            derivation="No source stereodescriptor beyond the separately recorded status was added to the candidate graph.")
        try:
            diag = st.stereo_inspect(smiles)
            m["stereochemistry"]["identity_key"] = st.fact(
                "RDKit stereo_graph_key", diag["stereo_graph_key"], evidence,
                method="structure_derivation", derivation="Computed from the candidate SMILES with standard stereo; source_verified=false.")
        except Exception:
            pass
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing()
        m["generic_template"] = None
    if kind == "terminal_allene":
        m["stereochemistry"]["axial"] = derive(
            "R-substituted terminal allene drawing R-CH=C=CH2",
            [{"type": "allene", "configuration": "not stereogenic at the allene axis", "basis": "terminal CH2"}],
            pages, "The depicted terminal allene has two identical hydrogens at one terminus; this does not rule out other stereogenic elements in R.")
        m["stereochemistry"]["chirality_status"] = derive(
            "terminal-allene axis inspected; other centers not assigned", "unspecified", pages,
            "The allene axis itself is not stereogenic, but the complete substituent was not asserted achiral when its structure is unresolved.")
        m["stereochemistry"]["geometry_status"] = derive(
            "cumulene drawing without E/Z descriptor", "not_applicable", pages,
            "Ordinary alkene E/Z labels do not apply to the terminal allene axis.")
        m["stereochemistry"]["representation_status"] = derive(
            "candidate/unencoded allene stereochemistry", "partial" if smiles else "not_encoded", pages,
            "Toolkit allene-axis support is limited and no source axial descriptor is present.")
    elif kind == "ez_product":
        m["stereochemistry"]["chirality_status"] = f(
            "no center/axis assignment reported", "unspecified", pages)
        m["stereochemistry"]["double_bonds"] = f(
            "boryl-functionalized skipped (E,Z)-diene", [{"configuration": "E"}, {"configuration": "Z"}], pages)
        m["stereochemistry"]["assignment"] = f(
            "product drawn and described as (E,Z)", {"geometry": "(E,Z)", "absolute_configuration": "unreported"}, pages)
        m["stereochemistry"]["geometry_status"] = f("(E,Z)", "specified", pages)
        m["stereochemistry"]["representation_status"] = derive(
            "source (E,Z) assignment retained outside SMILES", "not_encoded", pages,
            "A complete source-faithful isomeric graph was not asserted for the compact drawing.")
    elif kind == "deuterated_ez":
        m["stereochemistry"]["chirality_status"] = f("not assigned", "unspecified", pages)
        m["stereochemistry"]["double_bonds"] = f(
            "(E,Z)-diene; D on the (Z)-alkenyl carbon", [{"configuration": "E"}, {"configuration": "Z", "isotope": "D"}], pages)
        m["stereochemistry"]["assignment"] = f(
            "D selectively incorporated onto the (Z)-alkenyl carbon", {"geometry": "(E,Z)", "deuterium_position": "Z-alkenyl carbon"}, pages)
        m["stereochemistry"]["geometry_status"] = f("(E,Z)", "specified", pages)
        m["stereochemistry"]["representation_status"] = derive(
            "source geometry/isotope retained textually", "not_encoded", pages,
            "No complete isotope- and stereo-resolved candidate SMILES was asserted.")
    else:
        m["stereochemistry"]["chirality_status"] = f("not assigned", "unresolved", pages)
        m["stereochemistry"]["geometry_status"] = f(geometry or "not assigned", "specified" if geometry else "unresolved", pages)
        m["stereochemistry"]["representation_status"] = derive(
            "representation reviewed only to the stated level", "partial" if smiles else "unresolved", pages,
            "No unreported stereochemical assignment is inferred.")
        if assignment:
            m["stereochemistry"]["assignment"] = f(assignment, assignment, pages)
    m["identity_note"] = f(
        "Source-identified molecule; CAS not reported; machine record requires review.",
        "Source-identified molecule; CAS not reported; machine record requires review.", pages)
    return m


M = {}


def add_mol(mid, *args, **kwargs):
    M[mid] = molecule(mid, *args, **kwargs)
    return M[mid]["molecule_id"]


def mid(short):
    return M[short]["molecule_id"]


def participant(pid, mshort, role, pages, stage="reaction", amount=None, mass=None,
                volume=None, equiv=None, molpct=None, reference=None, product_role=None):
    p = st.new_participant(pid, mid(mshort), role, [EV[pages[0]], ET[pages[0]]], stage=stage)
    if product_role:
        p["product_role"] = f(product_role, product_role, pages)
    values = {"amount": amount, "mass": mass, "volume": volume,
              "equivalents": equiv, "mol_percent": molpct}
    for key, value in values.items():
        if value is not None:
            p["charges"][key] = f(value.get("raw_text") or f"{value['number']} {value['unit']}", value, pages)
    if reference:
        p["equivalence_reference"] = f(reference, mid(reference), pages)
    return p


def step(step_id, pages, participant_ids=(), operation=None, temperature=None,
         time=None, atmosphere=None, stage="reaction"):
    s = st.new_step(step_id, [EV[pages[0]], ET[pages[0]]], stage=stage)
    s["participant_ids"] = list(participant_ids)
    if operation is not None:
        s["operation"] = f(operation, operation, pages)
    if temperature is not None:
        raw = temperature if isinstance(temperature, str) else temperature.get("raw_text") or f"{temperature.get('number')} {temperature.get('unit')}"
        s["temperature"] = f(raw, temperature, pages)
    if time is not None:
        raw = time if isinstance(time, str) else time.get("raw_text") or f"{time.get('number')} {time.get('unit')}"
        s["time"] = f(raw, time, pages)
    if atmosphere is not None:
        s["atmosphere"] = f(atmosphere, atmosphere, pages)
    return s


def outcome(oid, kind, targets, pages, value, basis=None, yield_type=None, scope=None,
            raw=None):
    o = st.new_outcome(oid, kind, [mid(t) for t in targets], [EV[pages[0]], ET[pages[0]]])
    if raw is None:
        raw = value.get("raw_text") if isinstance(value, dict) else str(value)
        if raw is None:
            raw = f"{value.get('number')} {value.get('unit')}" if isinstance(value, dict) else str(value)
    o["value"] = f(raw, value, pages)
    if basis:
        o["measurement_basis"] = f(basis, basis, pages)
    if yield_type:
        o["yield_type"] = f(yield_type, yield_type, pages)
    if scope:
        o["scope"] = f(scope, scope, pages)
    return o


def classification(r, pages, labels=(), changes=(), coverage="partial", limitations=()):
    r["classification"]["coverage"] = derive(
        "class-level source drawing comparison", coverage, pages,
        "Class-level bond changes were inferred from source reactant/product drawings; no complete atom mapping is asserted.")
    r["classification"]["bond_changes"] = [derive(None, c, pages, "Class-level source drawing comparison without exact atom maps.") for c in changes]
    labrows = []
    for code, number, support in labels:
        labrows.append({"category_number": number, "code": code, "supporting_change_ids": list(support)})
    r["classification"]["labels"] = derive(None, labrows, pages, "Labels include all and only the stated formed/promoted class-level bonds.")
    r["classification"]["limitations"] = list(limitations) or [
        "Exact atom identities, complete hydrogen accounting, coproducts, and full atom mapping remain unresolved."]
    r["classification"]["model_votes"] = []


DIMER_CHANGES = (
    {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
    {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "B"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
    {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
)
DIMER_LABELS = (("C-C", 2, ("b1",)), ("C-B", 11, ("b2",)), ("C-H", 1, ("b3",)))


def reaction(rid, pages, record_type="specific_experiment", categories=(), name=None,
             transformation=None, participants=(), steps=(), outcomes=(), notes=(),
             general=(), unresolved=()):
    evidence = []
    for p in pages:
        evidence.extend([EV[p], ET[p]])
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = f(", ".join(categories), list(categories), pages) if categories else missing()
    r["location"] = f("; ".join(f"PDF page {p}" for p in pages), {"pages": list(pages)}, pages)
    r["reaction_name"] = f(name, name, pages) if name else missing()
    r["transformation_type"] = f(transformation, transformation, pages) if transformation else missing()
    r["participants"] = list(participants)
    r["steps"] = list(steps)
    r["outcomes"] = list(outcomes)
    r["general_scheme_ids"] = [f"{DOC}:{x}" for x in general]
    r["notes"] = [f(n, n, pages, method="inference",
                    derivation="Extractor summary of the linked source drawing/prose; not a verbatim quotation.") for n in notes]
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": list(unresolved)}
    r["completeness"] = "partial"
    return r


# Molecules: key source compounds and all independently numbered substrate/product identities.
add_mol("allene_generic", "1", [1, 2, 4], "generic terminal allene R-CH=C=CH2", kind="terminal_allene")
add_mol("product_generic", "2", [1, 2, 4], "generic boryl-functionalized skipped (E,Z)-diene", kind="ez_product")
scope_smiles = {
    "a": "C=C=Cc1ccc(OC)cc1", "b": "C=C=Cc1ccccc1",
    "c": "C=C=Cc1ccc(C(C)(C)C)cc1", "d": "C=C=Cc1ccc(F)cc1",
    "e": "C=C=Cc1ccc(Br)cc1", "f": "C=C=Cc1ccccc1C",
    "m": "C=C=CCc1ccccc1", "n": "C=C=CC1CCCCC1",
    "o": "C=C=C[Si](C)(C)c1ccccc1", "p": "C=C=CCOCc1ccccc1",
}
product_smiles = {
    "b": "c1ccccc1C=CCC(=Cc2ccccc2)B1OC(C)(C)C(C)(C)O1",
    "p": "c1ccccc1COCC=CCC(=CCOCc2ccccc2)B1OC(C)(C)C(C)(C)O1",
}
for letter in "abcdefghijklmnopqrstuvwxyz":
    add_mol(f"1{letter}", f"1{letter}", [2] + ([3, 4] if letter in "anp" else []),
            f"allene 1{letter}", smiles=scope_smiles.get(letter), kind="terminal_allene")
    add_mol(f"2{letter}", f"2{letter}", [2] + ([3, 4] if letter in "abnp" else []),
            f"boryl-functionalized skipped (E,Z)-diene 2{letter}", smiles=product_smiles.get(letter), kind="ez_product")
for label in "abcd":
    add_mol(f"3{label}", f"3{label}", [3], f"cross-dimerization product 3{label}", kind="ez_product")
for label in ["4", "5", "6", "7", "8", "9", "10"]:
    add_mol(label, label, [3], f"synthetic transformation product {label}")
add_mol("2p_d1", "2p-d1", [4], "deuterated (E,Z)-diene 2p-d1", kind="deuterated_ez")

# Reagents, catalysts, ligands, solvents, and collectively described species.
simple = {
    "HBpin": ("HBpin", "pinacolborane", "[BH]1OC(C)(C)C(C)(C)O1"),
    "DBpin": ("DBpin", "deuterated pinacolborane", None),
    "CrCl2": ("CrCl2", "chromium(II) chloride", "[Cl-].[Cl-].[Cr+2]"),
    "CrCl3": ("CrCl3", "chromium(III) chloride", "[Cl-].[Cl-].[Cl-].[Cr+3]"),
    "L1": ("mesPDI (L1)", "pyridine-2,6-diimine ligand mesPDI (L1)", None),
    "L2": ("L2", "PDI ligand L2", None), "L3": ("L3", "PDI ligand L3", None),
    "L4orL5": ("L4 or L5", "PDI ligand L4 or L5 grouped table variation", None),
    "L6orL7": ("L6 or L7", "PDI ligand L6 or L7 grouped table variation", None),
    "NaBHEt3": ("NaBHEt3", "sodium triethylborohydride", "[Na+].[BH-](CC)(CC)CC"),
    "EtMgBr": ("EtMgBr", "ethylmagnesium bromide", "CC[Mg]Br"),
    "NaH": ("NaH", "sodium hydride", "[H-].[Na+]"), "Mg": ("Mg", "magnesium powder", "[Mg]"),
    "Zn": ("Zn", "zinc powder", "[Zn]"), "MCl2": ("MCl2", "MnCl2, FeCl2, CoCl2, NiCl2, or CuCl2 grouped screen", None),
    "THF": ("THF", "tetrahydrofuran", "C1CCOC1"), "dioxane": ("1,4-dioxane", "1,4-dioxane", "C1COCCO1"),
    "DME": ("dimethoxyethane", "1,2-dimethoxyethane", "COCCOC"), "hexane": ("hexane", "hexane", "CCCCCC"),
    "toluene": ("toluene", "toluene", "Cc1ccccc1"), "NaBO3": ("NaBO3·nH2O", "sodium perborate hydrate", None),
    "mCPBA": ("m-CPBA", "meta-chloroperoxybenzoic acid", "OOC(=O)c1cccc(Cl)c1"),
    "paraformaldehyde": ("paraformaldehyde", "paraformaldehyde", None),
    "CoBr2": ("CoBr2", "cobalt(II) bromide", "[Br-].[Br-].[Co+2]"),
    "PrPDI": ("iPr-PDI (L3)", "isopropyl PDI ligand L3", None),
    "CH2Br2": ("CH2Br2", "dibromomethane", "BrCBr"), "Zn2": ("Zn", "zinc", "[Zn]"),
    "ArBr": ("ArBr", "1-bromo-3,5-bis(trifluoromethyl)benzene", "Brc1cc(C(F)(F)F)cc(C(F)(F)F)c1"),
    "nBuLi": ("n-BuLi", "n-butyllithium", "[Li]CCCC"), "Eschenmoser": ("Me2N=CH2 I-", "N,N-dimethylmethyleneammonium iodide", "C[N+](C)=C.[I-]"),
    "NiCl2": ("NiCl2", "nickel(II) chloride", "[Cl-].[Cl-].[Ni+2]"), "L6": ("L6", "ligand L6", None),
    "EtOSiH": ("(EtO)3SiH", "triethoxysilane", "CCO[SiH](OCC)OCC"), "CsF": ("CsF", "cesium fluoride", "[Cs+].[F-]"),
    "arylI_H": ("aryl iodide (R=H)", "iodobenzene", "Ic1ccccc1"), "arylI_OMe": ("aryl iodide (R=OMe)", "iodoanisole isomer shown in Scheme 2F", None),
    "mesPDICrCl": ("(mesPDI)CrCl", "preformed mesPDI-ligated chromium chloride", None),
}
for short, (label, name, smiles) in simple.items():
    add_mol(short, label, [2, 3, 4], name, smiles=smiles)


R = []

# Cited precedent and current general scheme.
prior = reaction("scheme1_prior_borane_dimerization", (1,), "literature_example", (),
                 "borane-induced dimerization of arylallenes", "stoichiometric arylallene dimerization",
                 participants=(participant("prior_allene", "allene_generic", "reactant", (1,)),
                               participant("prior_product", "product_generic", "product", (1,), product_role="reported prior-work product")),
                 notes=("Cited Erker precedent: (C6F5)2BH, toluene, room temperature, 3 days, five examples, up to 54% yield; product described as an (E,E)-1,5-diene.",),
                 unresolved=("Generic prior-work structures and per-example conditions are not expanded on the supplied page.",))
classification(prior, (1,), coverage="unresolved")
R.append(prior)

generic = reaction("scheme1_current_generic", (1, 2), "generic_scheme", (),
                   "chromium-catalyzed selective dimerization/hydroboration of allenes",
                   "allene dimerization/hydroboration",
                   participants=(participant("g_allene", "allene_generic", "reactant", (1,)),
                                 participant("g_hbpin", "HBpin", "reagent", (1,)),
                                 participant("g_product", "product_generic", "product", (1,), product_role="generic desired product")),
                   notes=("This-work generic scheme couples two allene molecules with HBpin to a boryl-functionalized skipped (E,Z)-diene.",),
                   unresolved=("R-group identities and exact atom mapping are generic.",))
classification(generic, (1, 2), DIMER_LABELS, DIMER_CHANGES)
R.append(generic)


def common_dimer_participants(allene="1a", product="2a", pages=(2,), scope=False, solvent="THF",
                              catalyst="CrCl2", ligand="L1", activator="NaBHEt3"):
    if scope:
        return [
            participant("allene", allene, "reactant", pages, amount=q(1.00, "mmol", "1.00 mmol")),
            participant("hbpin", "HBpin", "reagent", pages, amount=q(0.600, "mmol", "0.600 mmol")),
            participant("cr", catalyst, "catalyst", pages, amount=q(15.0, "umol", "15.0 umol"), molpct=q(3.0, "mol%", "3 mol%"), reference=allene),
            participant("ligand", ligand, "ligand", pages, amount=q(15.0, "umol", "15.0 umol"), molpct=q(3.0, "mol%", "3 mol%"), reference=allene),
            participant("activator", activator, "reagent", pages, amount=q(30.0, "umol", "30.0 umol"), molpct=q(6.0, "mol%", "6 mol%"), reference=allene),
            participant("solvent", solvent, "solvent", pages, volume=q(1.0, "mL", "1 mL")),
            participant("product", product, "product", pages, product_role="desired product"),
        ]
    return [
        participant("allene", allene, "reactant", pages, amount=q(0.200, "mmol", "0.200 mmol")),
        participant("hbpin", "HBpin", "reagent", pages, amount=q(0.120, "mmol", "0.120 mmol")),
        participant("cr", catalyst, "catalyst", pages, amount=q(3.0, "umol", "3.0 umol"), molpct=q(3.0, "mol%", "3 mol%"), reference=allene),
        participant("ligand", ligand, "ligand", pages, amount=q(3.0, "umol", "3.0 umol"), molpct=q(3.0, "mol%", "3 mol%"), reference=allene),
        participant("activator", activator, "reagent", pages, amount=q(6.0, "umol", "6.0 umol"), molpct=q(6.0, "mol%", "6 mol%"), reference=allene),
        participant("solvent", "THF", "solvent", pages, volume=q(0.2, "mL", "0.2 mL")),
        participant("product", product, "product", pages, product_role="desired product"),
    ]


# Table 1: all twelve printed rows; row 11 is split into the three explicitly named activator experiments.
opt = [
    ("1", "none", ">99%", 78, None), ("2", "without CrCl2", None, None, "no reaction"),
    ("3", "without L1", None, None, "no reaction"), ("4", "without NaBHEt3", None, None, "no reaction"),
    ("5", "L2 instead of L1", ">99%", 65, "hydrodimerization byproducts"),
    ("6", "L3 instead of L1", ">99%", 23, "hydrodimerization product was major"),
    ("7", "L4 or L5 instead of L1", "52%", "<5%", "hydroboration occurred predominantly"),
    ("8", "L6 or L7 instead of L1", "50%", "<5%", "NaBHEt3 9.0 umol"),
    ("9", "CrCl3 instead of CrCl2", ">99%", 73, None),
    ("10", "EtMgBr instead of NaBHEt3", ">99%", 64, "some hydrodimerization products"),
    ("11a", "NaH instead of NaBHEt3", "<5%", "<5%", "grouped in printed entry 11"),
    ("11b", "Mg powder instead of NaBHEt3", "<5%", "<5%", "grouped in printed entry 11"),
    ("11c", "Zn powder instead of NaBHEt3", "<5%", "<5%", "grouped in printed entry 11"),
    ("12", "MCl2 (M=Mn, Fe, Co, Ni, Cu) instead of CrCl2", "30-50%", "<5%", "grouped metal screen; mixed hydrogenation/hydroboration/hydrodimerization products"),
]
for entry, variation, conv, yld, note in opt:
    ps = common_dimer_participants()
    # Apply the printed variation without silently keeping the replaced component.
    if entry == "2": ps = [p for p in ps if p["participant_id"] != "cr"]
    elif entry == "3": ps = [p for p in ps if p["participant_id"] != "ligand"]
    elif entry == "4": ps = [p for p in ps if p["participant_id"] != "activator"]
    elif entry == "5": ps = [p if p["participant_id"] != "ligand" else participant("ligand", "L2", "ligand", (2,), amount=q(3.0,"umol"), molpct=q(3.0,"mol%"), reference="1a") for p in ps]
    elif entry == "6": ps = [p if p["participant_id"] != "ligand" else participant("ligand", "L3", "ligand", (2,), amount=q(3.0,"umol"), molpct=q(3.0,"mol%"), reference="1a") for p in ps]
    elif entry == "7": ps = [p if p["participant_id"] != "ligand" else participant("ligand", "L4orL5", "ligand", (2,), amount=q(3.0,"umol"), molpct=q(3.0,"mol%"), reference="1a") for p in ps]
    elif entry == "8":
        ps = [p if p["participant_id"] != "ligand" else participant("ligand", "L6orL7", "ligand", (2,), amount=q(3.0,"umol"), molpct=q(3.0,"mol%"), reference="1a") for p in ps]
        ps = [p if p["participant_id"] != "activator" else participant("activator", "NaBHEt3", "reagent", (2,), amount=q(9.0,"umol"), molpct=q(9.0,"mol%"), reference="1a") for p in ps]
    elif entry == "9":
        ps = [p if p["participant_id"] != "cr" else participant("cr", "CrCl3", "catalyst", (2,), amount=q(3.0,"umol"), molpct=q(3.0,"mol%"), reference="1a") for p in ps]
        ps = [p if p["participant_id"] != "activator" else participant("activator", "NaBHEt3", "reagent", (2,), amount=q(9.0,"umol"), molpct=q(9.0,"mol%"), reference="1a") for p in ps]
    elif entry == "10": ps = [p if p["participant_id"] != "activator" else participant("activator", "EtMgBr", "reagent", (2,), amount=q(6.0,"umol"), molpct=q(6.0,"mol%"), reference="1a") for p in ps]
    elif entry.startswith("11"):
        repl = {"11a":"NaH","11b":"Mg","11c":"Zn"}[entry]
        ps = [p if p["participant_id"] != "activator" else participant("activator", repl, "reagent", (2,), amount=q(6.0,"umol"), reference="1a") for p in ps]
    elif entry == "12": ps = [p if p["participant_id"] != "cr" else participant("cr", "MCl2", "catalyst", (2,), amount=q(3.0,"umol"), molpct=q(3.0,"mol%"), reference="1a") for p in ps]
    outs = []
    if conv:
        if conv == ">99%": cv = q(99, "%", ">99%", "gt")
        elif conv == "<5%": cv = q(5, "%", "<5%", "lt")
        elif conv == "52%": cv = q(52, "%")
        elif conv == "50%": cv = q(50, "%")
        else: cv = {"minimum":30,"maximum":50,"unit":"%","raw_text":"30-50%"}
        outs.append(outcome("conversion", "conversion", ("1a",), (2,), cv, raw=conv))
    if yld is not None:
        yv = q(yld, "%") if isinstance(yld, int) else q(5, "%", "<5%", "lt")
        outs.append(outcome("yield", "yield", ("2a",), (2,), yv, basis="isolated yield", yield_type="isolated", raw=str(yld) if isinstance(yld,int) else yld))
    if note == "no reaction": outs.append(outcome("no_reaction", "no_reaction", (), (2,), True, raw="reaction did not occur"))
    notes = [f"Table 1 entry {entry.replace('11a','11').replace('11b','11').replace('11c','11')}: {variation}."]
    if note: notes.append(note)
    r = reaction(f"table1_entry{entry}", (2,), categories=("optimization",),
                 name="optimization of allene dimerization/hydroboration",
                 transformation="allene dimerization/hydroboration", participants=ps,
                 steps=(step("reaction", (2,), [p["participant_id"] for p in ps if p["stage"]=="reaction"],
                             "combine under the tabulated conditions", "RT", q(8,"h")),),
                 outcomes=outs, notes=notes, general=("scheme1_current_generic",),
                 unresolved=("Conversion measurement method is unreported.", "SI procedure was not supplied.") if conv else ("SI procedure was not supplied.",))
    r["solvent_conditions"]["composition"] = f("THF", [mid("THF")], (2,))
    r["solvent_conditions"]["total_volume"] = f("0.2 mL", q(0.2,"mL"), (2,))
    r["reaction_scale"] = f("1a (0.200 mmol)", {"reference_molecule_id":mid("1a"),"amount":q(0.200,"mmol")}, (2,))
    if note == "no reaction" or (isinstance(yld,str) and yld.startswith("<")):
        classification(r, (2,), coverage="not_applicable" if note == "no reaction" else "unresolved")
    else: classification(r, (2,), DIMER_LABELS, DIMER_CHANGES)
    R.append(r)

# Table 2: 26 independently printed substrate-scope products.
scope_yields = dict(a=76,b=75,c=76,d=80,e=80,f=78,g=65,h=67,i=60,j=72,k=87,l=62,m=72,n=87,o=62,p=75,q=52,r=76,s=75,t=73,u=68,v=60,w=77,x=80,y=66,z=68)
for letter, yld in scope_yields.items():
    temp = q(50,"°C") if letter in "iq" else q(40,"°C") if letter in "tu" else "RT"
    ps = common_dimer_participants(f"1{letter}", f"2{letter}", scope=True)
    r = reaction(f"table2_2{letter}", (2,), categories=("substrate_scope",),
                 name="Cr-catalyzed dimerization/hydroboration of allene",
                 transformation="allene dimerization/hydroboration", participants=ps,
                 steps=(step("reaction", (2,), [p["participant_id"] for p in ps], "react under Table 2 conditions", temp, q(8,"h")),),
                 outcomes=(outcome("yield", "yield", (f"2{letter}",), (2,), q(yld,"%"), basis="isolated yield", yield_type="isolated"),),
                 notes=(f"Table 2 product 2{letter}; source assigns the product family as boryl-containing skipped (E,Z)-dienes.",),
                 general=("scheme1_current_generic",), unresolved=("No associated SI was supplied; workup and purification are unreported on pages 1-6.",))
    r["reaction_scale"] = f("allene (1.00 mmol)", {"reference_molecule_id":mid(f"1{letter}"),"amount":q(1.00,"mmol")}, (2,))
    r["solvent_conditions"]["composition"] = f("THF", [mid("THF")], (2,))
    r["solvent_conditions"]["total_volume"] = f("1 mL", q(1,"mL"), (2,))
    classification(r, (2,), DIMER_LABELS, DIMER_CHANGES)
    R.append(r)

# Separate gram-scale reaction.
gps = common_dimer_participants("1p", "2p", (2,), scope=True)
for p in gps:
    if p["participant_id"] == "allene": p["charges"]["amount"] = f("6.00 mmol", q(6.00,"mmol"), (2,))
    elif p["participant_id"] == "hbpin": p["charges"]["amount"] = f("scaled 0.600:1.00 ratio from Table 2", q(3.60,"mmol"), (2,), method="inference", derivation="Calculated only as an inherited ratio; exact gram-scale charge is not stated.")
gram = reaction("gram_scale_2p", (2,), categories=("scale_up",), name="gram-scale dimerization/hydroboration of 1p",
                transformation="allene dimerization/hydroboration", participants=gps,
                steps=(step("reaction", (2,), [p["participant_id"] for p in gps], "gram-scale reaction under the scope conditions", "RT", q(8,"h")),),
                outcomes=(outcome("amount", "product_amount", ("2p",), (2,), q(1.03,"g"), raw="1.03 g"),
                          outcome("yield", "yield", ("2p",), (2,), q(72,"%"), basis="isolated yield", yield_type="isolated")),
                notes=("The article reports 1p (6.00 mmol), 2p (1.03 g), and 72% isolated yield.",),
                general=("scheme1_current_generic",), unresolved=("Exact scaled HBpin/catalyst/ligand/activator charges are not explicitly stated; inherited ratios require review.",))
gram["reaction_scale"] = f("1p (6.00 mmol)", {"reference_molecule_id":mid("1p"),"amount":q(6.00,"mmol")}, (2,))
gram["solvent_conditions"]["composition"] = f("THF (inherited Table 2 condition)", [mid("THF")], (2,), method="inherited_general_condition")
classification(gram, (2,), DIMER_LABELS, DIMER_CHANGES)
R.append(gram)

# Equation 1 cross-dimerizations.
cross_yields = dict(a=45,b=40,c=51,d=42)
for label, yld in cross_yields.items():
    ps = [participant("allene_R1", "allene_generic", "reactant", (3,), equiv=q(1,"equiv"), reference="allene_generic"),
          participant("allene_R2", "allene_generic", "reactant", (3,), equiv=q(2,"equiv"), reference="allene_generic"),
          participant("hbpin", "HBpin", "reagent", (3,)),
          participant("cr", "CrCl2", "catalyst", (3,), molpct=q(3,"mol%"), reference="allene_generic"),
          participant("ligand", "L1", "ligand", (3,), molpct=q(3,"mol%"), reference="allene_generic"),
          participant("activator", "NaBHEt3", "reagent", (3,), molpct=q(6,"mol%"), reference="allene_generic"),
          participant("solvent", "THF", "solvent", (3,)),
          participant("product", f"3{label}", "product", (3,), product_role="desired cross-dimerization product")]
    r = reaction(f"equation1_3{label}", (3,), categories=("substrate_scope",),
                 name="cross-dimerization/hydroboration of two different allenes", transformation="cross-dimerization/hydroboration",
                 participants=ps, steps=(step("reaction", (3,), [p["participant_id"] for p in ps], "react at the displayed 1:2 allene ratio", "RT", q(8,"h")),),
                 outcomes=(outcome("yield", "yield", (f"3{label}",), (3,), q(yld,"%"), basis="isolated yield", yield_type="isolated"),),
                 notes=("The four cross/homo product ratios are collectively reported as varying from 6:1 to 15:1; no per-product ratio is assigned.",),
                 general=("scheme1_current_generic",), unresolved=("Individual cross:homo ratio and exact R1/R2 source-molecule identities are not fully assigned from the compact figure.",))
    r["solvent_conditions"]["composition"] = f("THF", [mid("THF")], (3,))
    classification(r, (3,), DIMER_LABELS, DIMER_CHANGES)
    R.append(r)


def derivative(rid, start, product, reagents, temperature, duration, yld, labels, changes, notes=()):
    ps = [participant("start", start, "reactant", (3,))]
    for i, (mshort, role, amount) in enumerate(reagents):
        kw = {"equiv": amount} if amount and amount.get("unit") == "equiv" else {"molpct": amount} if amount and amount.get("unit") == "mol%" else {}
        ps.append(participant(f"reagent{i+1}", mshort, role, (3,), reference=start if amount and amount.get("unit") in ("equiv","mol%") else None, **kw))
    ps.append(participant("product", product, "product", (3,), product_role="reported product"))
    r = reaction(rid, (3,), categories=("derivatization",), name="synthetic transformation of boryl-functionalized skipped diene",
                 transformation="product derivatization", participants=ps,
                 steps=(step("reaction", (3,), [p["participant_id"] for p in ps], "react under Scheme 2 conditions", temperature, duration),),
                 outcomes=(outcome("yield", "yield", (product,), (3,), q(yld,"%"), basis="isolated yield", yield_type="isolated"),),
                 notes=notes, unresolved=("Workup and purification are not reported on pages 1-6; SI was not supplied.",))
    classification(r, (3,), labels, changes)
    return r


c_co = ({"change_id":"b1","atoms":[{"atom_map":None,"element":"C"},{"atom_map":None,"element":"O"}],"before_order":0,"after_order":1,"change_type":"formed"},)
R.append(derivative("scheme2A_2p_to_4", "2p", "4", (("NaBO3","reagent",q(4,"equiv")),("THF","solvent",None)), "RT", q(3,"h"), 95,
                    (("C-O",3,("b1",)),), c_co))
epox_changes = (
    {"change_id":"b1","atoms":[{"atom_map":None,"element":"C"},{"atom_map":None,"element":"O"}],"before_order":0,"after_order":1,"change_type":"formed"},
    {"change_id":"b2","atoms":[{"atom_map":None,"element":"C"},{"atom_map":None,"element":"O"}],"before_order":0,"after_order":1,"change_type":"formed"},)
R.append(derivative("scheme2B_4_to_5", "4", "5", (("mCPBA","reagent",q(1.1,"equiv")),), q(0,"°C"), q(3,"h"), 72,
                    (("C-O",3,("b1","b2")),), epox_changes, ("m-CPBA epoxidation occurred selectively at the trisubstituted alkene.",)))
c_cc = ({"change_id":"b1","atoms":[{"atom_map":None,"element":"C"},{"atom_map":None,"element":"C"}],"before_order":0,"after_order":1,"change_type":"formed"},)
R.append(derivative("scheme2C_2p_to_6", "2p", "6", (("paraformaldehyde","reagent",q(5,"equiv")),("toluene","solvent",None),("THF","solvent",None)), q(90,"°C"), q(6,"h"), 92,
                    (("C-C",2,("b1",)),), c_cc))
cyclo = (
    {"change_id":"b1","atoms":[{"atom_map":None,"element":"C"},{"atom_map":None,"element":"C"}],"before_order":0,"after_order":1,"change_type":"formed"},
    {"change_id":"b2","atoms":[{"atom_map":None,"element":"C"},{"atom_map":None,"element":"C"}],"before_order":0,"after_order":1,"change_type":"formed"},)
R.append(derivative("scheme2D_6_to_7", "6", "7", (("CoBr2","catalyst",q(5,"mol%")),("PrPDI","ligand",q(6,"mol%")),("CH2Br2","reagent",q(1.5,"equiv")),("Zn2","reagent",q(2,"equiv")),("THF","solvent",None)), "RT", q(24,"h"), 80,
                    (("C-C",2,("b1","b2")),), cyclo))
R.append(derivative("scheme2E_2p_to_8", "2p", "8", (("ArBr","reagent",q(1.2,"equiv")),("nBuLi","reagent",q(1.2,"equiv")),("THF","solvent",None),("Eschenmoser","reagent",q(1.2,"equiv"))), "-78 °C to RT", q(6,"h"), 86,
                    (("C-C",2,("b1",)),), c_cc, ("Scheme 2E also reports an initial -78 °C, 1 h stage before addition of the ammonium reagent and warming; this is expanded below.",)))
# Replace Scheme 2E's single summary step with ordered stages.
R[-1]["steps"] = [step("organolithium_stage", (3,), [p["participant_id"] for p in R[-1]["participants"][:3]], "combine ArBr and n-BuLi with 2p", q(-78,"°C"), q(1,"h")),
                    step("aminomethylation_stage", (3,), [p["participant_id"] for p in R[-1]["participants"]], "add Me2N=CH2 I- and warm", "-78 °C to RT", q(6,"h"))]
R[-1]["steps"][1]["temperature"] = f("-78 °C to RT", "temperature transition", (3,))
hydroary_changes = (c_cc[0], {"change_id":"b2","atoms":[{"atom_map":None,"element":"C"},{"atom_map":None,"element":"H"}],"before_order":0,"after_order":1,"change_type":"hydrogen_gain","hydrogen_count":1})
hydroary_labels = (("C-C",2,("b1",)),("C-H",1,("b2",)))
for suffix, aryl, prod, yld in [("R_H","arylI_H","9",60),("R_OMe","arylI_OMe","10",68)]:
    R.append(derivative(f"scheme2F_{suffix}", "2b", prod, ((aryl,"reactant",None),("NiCl2","catalyst",q(5,"mol%")),("L6","ligand",q(6,"mol%")),("EtOSiH","reagent",q(2.5,"equiv")),("CsF","reagent",q(2,"equiv")),("THF","solvent",None)), "RT", q(16,"h"), yld,
                        hydroary_labels, hydroary_changes))

# Mechanistic controls and measurements.
def mech_measure(rid, pages, participants, observation, steps=(), unresolved=()):
    r = reaction(rid, pages, categories=("mechanistic_study",), name="mechanistic experiment", transformation="mechanistic observation",
                 participants=participants, steps=steps,
                 outcomes=(outcome("observation", "other", (), pages, observation, raw=str(observation)),),
                 unresolved=unresolved)
    classification(r, pages, coverage="not_applicable")
    return r

R.append(mech_measure("epr_CrCl2_L1_no_signal", (3,),
                      [participant("cr","CrCl2","reactant",(3,),equiv=q(1,"equiv"),reference="L1"), participant("ligand","L1","ligand",(3,),equiv=q(1,"equiv"),reference="L1")],
                      {"EPR_signal":"none observed","molar_ratio":"1:1"},
                      (step("epr_sample_preparation",(3,),["cr","ligand"],"mix CrCl2 and mesPDI in a 1:1 molar ratio",None,None),),
                      ("Solvent, concentration, temperature, duration, and acquisition settings are unreported on pages 1-6.",)))
R.append(mech_measure("epr_CrCl2_L1_NaBHEt3", (3,),
                      [participant("cr","CrCl2","reactant",(3,)),participant("ligand","L1","ligand",(3,)),participant("activator","NaBHEt3","reagent",(3,),equiv=q(2,"equiv"),reference="CrCl2")],
                      {"EPR":"singlet resonance","g_value":2.001},
                      (step("activation_and_epr",(3,),["cr","ligand","activator"],"add 2 equiv NaBHEt3 to the CrCl2/mesPDI mixture and record EPR",None,None),),
                      ("Solvent, concentration, temperature, duration, and acquisition settings are not reported on pages 1-6.",)))
R.append(mech_measure("epr_standard_reaction_1a", (3,), common_dimer_participants(), {"EPR":"singlet resonance","g_value":1.998},
                      (step("reaction_monitoring",(3,),["allene","hbpin","cr","ligand","activator","solvent"],"monitor the catalytic reaction by EPR","RT",None),),
                      ("EPR acquisition settings are unavailable without SI.",)))
cr3ps = common_dimer_participants(catalyst="CrCl3")
R.append(mech_measure("epr_CrCl3_reaction", (3,), cr3ps, {"EPR":"singlet resonance","g_value":1.999},
                      (step("reaction_monitoring",(3,),[p["participant_id"] for p in cr3ps],"monitor the CrCl3/mesPDI reaction by EPR","RT",None),),
                      ("EPR acquisition settings are unavailable without SI.",)))

for suffix, active, yld in [("without_NaBHEt3",False,0),("with_NaBHEt3",True,72)]:
    ps=[participant("allene","1a","reactant",(3,)),participant("hbpin","HBpin","reagent",(3,)),participant("complex","mesPDICrCl","catalyst",(3,),molpct=q(5,"mol%"),reference="1a"),participant("solvent","THF","solvent",(3,)),participant("product","2a","product",(3,),product_role="desired product")]
    if active: ps.insert(3,participant("activator","NaBHEt3","reagent",(3,),equiv=q(1.2,"equiv"),reference="mesPDICrCl"))
    r=reaction(f"equation2_{suffix}",(3,),categories=("control","mechanistic_study"),name="preformed (mesPDI)CrCl control",transformation="allene dimerization/hydroboration",participants=ps,
               steps=(step("reaction",(3,),[p["participant_id"] for p in ps],"react under Equation 2 conditions","RT",q(8,"h")),),
               outcomes=(outcome("yield","yield",("2a",),(3,),q(yld,"%"),raw=f"{yld}%"),),
               unresolved=("Yield measurement basis and reagent charges are unreported.",))
    classification(r,(3,),DIMER_LABELS,DIMER_CHANGES) if active else classification(r,(3,),coverage="not_applicable")
    R.append(r)
R.append(mech_measure("epr_preformed_complex_with_NaBHEt3", (3,),
                      [participant("complex","mesPDICrCl","reactant",(3,)),participant("activator","NaBHEt3","reagent",(3,),equiv=q(1.2,"equiv"),reference="mesPDICrCl")],
                      {"EPR":"singlet resonance","g_value":1.998},
                      (step("activation_and_epr",(3,),["complex","activator"],"combine preformed (mesPDI)CrCl and NaBHEt3 and record EPR",None,None),),
                      ("Solvent, concentration, temperature, duration, and other details are in unavailable Supporting Information Figure S5.",)))

# Equation 3 deuterium labeling.
dps=[participant("allene","1p","reactant",(4,)),participant("dbpin","DBpin","reagent",(4,)),participant("cr","CrCl2","catalyst",(4,),molpct=q(3,"mol%"),reference="1p"),participant("ligand","L1","ligand",(4,),molpct=q(3,"mol%"),reference="1p"),participant("activator","NaBHEt3","reagent",(4,),molpct=q(6,"mol%"),reference="1p"),participant("solvent","THF","solvent",(4,)),participant("product","2p_d1","product",(4,),product_role="deuterated product")]
deut=reaction("equation3_deuterium_labeling",(3,4),categories=("mechanistic_study",),name="deuterium-labeling dimerization/hydroboration",transformation="dimerization/deuteroboration",participants=dps,
              steps=(step("reaction",(4,),[p["participant_id"] for p in dps],"react 1p with DBpin under standard conditions","RT",q(8,"h")),),
              outcomes=(outcome("yield","yield",("2p_d1",),(4,),q(74,"%")),
                        outcome("deuterium_incorporation","selectivity",("2p_d1",),(4,),q(95,"%"),scope={"metric":"D incorporation","position":"Z-alkenyl carbon"})),
              unresolved=("Yield and isotope-measurement methods are unreported; SI is unavailable.",))
classification(deut,(3,4),DIMER_LABELS,DIMER_CHANGES,limitations=("C-D is indexed under the hydrogen-gain class C-H; isotope identity is retained in the outcome.","Exact atom mapping and all other bond-order changes remain unresolved."))
R.append(deut)

# Scheme 3A solvent screen: five actual experiments.
for solvent, yld in [("THF",95),("dioxane",74),("DME",83),("hexane",0),("toluene",0)]:
    ps=common_dimer_participants("1n","2n",(4,),scope=False,solvent=solvent)
    # Scheme 3 reports mol% but not the Table-1 micromole/volume charges; remove inherited absolute amounts.
    for p in ps:
        if p["participant_id"] in ("allene","hbpin","cr","ligand","activator","solvent"):
            p["charges"]["amount"] = missing()
            p["charges"]["volume"] = missing()
    r=reaction(f"scheme3A_solvent_{solvent}",(4,),categories=("optimization","mechanistic_study"),name="solvent-effect experiment",transformation="allene dimerization/hydroboration",participants=ps,
               steps=(step("reaction",(4,),[p["participant_id"] for p in ps],"react in the displayed solvent","RT",q(8,"h")),),
               outcomes=(outcome("yield","yield",("2n",),(4,),q(yld,"%"),basis="GC yield" if yld else None,yield_type="GC" if yld else None),),
               notes=("Prose reports full conversion in coordinating solvents and no reaction in toluene/hexane.",),unresolved=("Allene/HBpin quantities, solvent volume, and GC method are unreported.",))
    classification(r,(4,),DIMER_LABELS,DIMER_CHANGES) if yld else classification(r,(4,),coverage="not_applicable")
    R.append(r)

R.append(mech_measure("scheme3B_epr_toluene",(4,),
                      [participant("allene","1n","reactant",(4,)),participant("hbpin","HBpin","reagent",(4,)),participant("cr","CrCl2","catalyst",(4,),molpct=q(3,"mol%"),reference="1n"),participant("ligand","L1","ligand",(4,),molpct=q(3,"mol%"),reference="1n"),participant("solvent","toluene","solvent",(4,))],
                      {"EPR":"very weak signal","g_value":2.004},
                      (step("reaction_monitoring",(4,),["allene","hbpin","cr","ligand","solvent"],"monitor the reaction in toluene by EPR","RT",None),),
                      ("NaBHEt3 charge, duration, and EPR acquisition settings are not stated in the Scheme 3B caption/prose.",)))
R.append(mech_measure("scheme3C_epr_THF_addition",(4,),
                      [participant("allene","1n","reactant",(4,)),participant("hbpin","HBpin","reagent",(4,)),participant("cr","CrCl2","catalyst",(4,),molpct=q(3,"mol%"),reference="1n"),participant("ligand","L1","ligand",(4,),molpct=q(3,"mol%"),reference="1n"),participant("solvent","toluene","solvent",(4,)),participant("thf","THF","additive",(4,),equiv=q(0.5,"equiv"),reference="1n")],
                      {"EPR":"significant active chromium signal","g_value":1.997},
                      (step("thf_addition_and_epr",(4,),["allene","hbpin","cr","ligand","solvent","thf"],"add 0.5 equiv THF to the toluene reaction mixture and record EPR","RT",None),),
                      ("NaBHEt3 charge, duration, and EPR acquisition settings are not stated in the Scheme 3C caption/prose.",)))

ps=common_dimer_participants("1n","2n",(4,),scope=False,solvent="toluene")
for p in ps:
    if p["participant_id"] in ("allene","hbpin","cr","ligand","activator","solvent"):
        p["charges"]["amount"] = missing(); p["charges"]["volume"] = missing()
ps.insert(-1,participant("thf_additive","THF","additive",(4,),molpct=q(50,"mol%"),reference="1n"))
s3d=reaction("scheme3D_toluene_plus_THF",(4,),categories=("optimization","mechanistic_study"),name="THF-addition solvent-effect experiment",transformation="allene dimerization/hydroboration",participants=ps,
             steps=(step("reaction",(4,),[p["participant_id"] for p in ps],"react in toluene with 50 mol% THF","RT",q(8,"h")),),
             outcomes=(outcome("yield","yield",("2n",),(4,),q(92,"%"),basis="GC yield",yield_type="GC"),),
             unresolved=("Substrate and HBpin charges, total solvent volume, and GC protocol are unreported.",))
classification(s3d,(4,),DIMER_LABELS,DIMER_CHANGES);R.append(s3d)

# Main-text kinetic study (individual plots reside in unavailable SI).
kps=[participant("allene","1n","reactant",(4,)),participant("hbpin","HBpin","reagent",(4,)),participant("cr","CrCl2","catalyst",(4,)),participant("ligand","L1","ligand",(4,)),participant("solvent","THF","solvent",(4,))]
kin=reaction("kinetic_study_1n",(4,5),categories=("mechanistic_study",),name="kinetic study of allene dimerization/hydroboration",transformation="kinetic analysis",participants=kps,
             steps=(step("rate_measurement",(4,),[p["participant_id"] for p in kps],"measure initial rates","RT",None),),
             outcomes=(outcome("orders","other",(),(4,),{"Cr_catalyst":1,"allene_1n":0,"HBpin":0},raw="first order in chromium catalyst and zero order in allene and HBpin"),),
             unresolved=("Concentrations, time-course data, and measurement protocol are only in unavailable SI Figures S8-S11.",))
classification(kin,(4,5),coverage="not_applicable");R.append(kin)

mechanism=reaction("scheme4_proposed_mechanism",(4,5),"proposed_mechanism",(),"proposed Cr(I)-hydride catalytic cycle","proposed mechanism",
                   participants=(participant("allene","allene_generic","reactant",(4,)),participant("hbpin","HBpin","reagent",(4,)),participant("product","product_generic","product",(4,),product_role="generic product")),
                   notes=("Author-proposed cycle I-VI: Cr(I)-hydride formation, THF dissociation, two allene insertions, and reaction with HBpin; THF dissociation is proposed as turnover-limiting.",),
                   unresolved=("Proposed intermediates are not treated as isolated experiments or verified species.",))
classification(mechanism,(4,5),coverage="unresolved");R.append(mechanism)

# Source context: metadata, spectroscopy/characterization, non-synthetic pages, and negative inventory.
source_context = [
    {"context_id":f"{DOC}:context:metadata", "document_id":DOC, "kind":"bibliographic_metadata",
     "evidence_ids":[EV[1],ET[1]], "content":{"title":"Chromium-Catalyzed Selective Dimerization/Hydroboration of Allenes to Access Boryl-Functionalized Skipped (E,Z)-Dienes","authors":["Yinsong Zhao","Shaozhong Ge"],"doi":"10.1002/anie.202012344"}},
    {"context_id":f"{DOC}:context:epr", "document_id":DOC, "kind":"spectroscopy_context",
     "evidence_ids":[EV[3],ET[3],EV[4],ET[4]], "content":{"note":"EPR spectra and author interpretation are linked to concrete mechanistic-study records; spectral plots are not synthetic experiments."}},
    {"context_id":f"{DOC}:context:pages5_6", "document_id":DOC, "kind":"non_reaction_content",
     "evidence_ids":[EV[5],ET[5],EV[6],ET[6]], "content":{"page5":"summary, acknowledgements, conflict statement, keywords, and references","page6":"references and manuscript dates"}},
    {"context_id":f"{DOC}:context:negative_spectral_inventory", "document_id":DOC, "kind":"negative_inventory",
     "evidence_ids":[EV[p] for p in range(1,7)], "content":{"finding":"No standalone mass-spectrum peak table, molecular-weight table, or fragment-relative-abundance block occurs on pages 1-6."}},
]

bundle = copy.deepcopy(json.loads((RUN/"raw_model/standard_init.json").read_text()))
bundle["molecules"] = list(M.values())
bundle["reactions"] = R
bundle["source_context"] = source_context

# Page coverage and raw inventory.
by_page = {p: [] for p in range(1,7)}
for r in R:
    for p in range(1,7):
        if EV[p] in r["evidence_ids"]:
            by_page[p].append(r["reaction_id"])
coverage = {"document_id":DOC,"semantic_accuracy_verified":False,"pages":[
    {"page":1,"state":"attempted","regions":[
        {"region":"title, abstract, introduction, and Scheme 1","finding":"Bibliographic metadata; one cited borane-mediated prior-work scheme; one current generic dimerization/hydroboration scheme; bioactive-molecule drawings are context.","evidence_ids":[EV[1],ET[1]]}],"record_ids":by_page[1],"unresolved":["Online SI is advertised but was not supplied or associated."]},
    {"page":2,"state":"attempted","regions":[
        {"region":"Table 1 and footnotes","finding":"Twelve printed optimization rows inventoried; grouped entry 11 is represented by three named activator records, while grouped entry 12 remains an inseparable five-metal screen.","evidence_ids":[EV[2],ET[2]]},
        {"region":"Table 2 and caption","finding":"All 26 products 2a-2z recorded independently with isolated yields and temperature footnotes; prose-only 1p gram-scale reaction is separate.","evidence_ids":[EV[2],ET[2]]}],"record_ids":by_page[2],"unresolved":["Compact structures for several functionalized substrates/products were not assigned complete SMILES.","Entry 12 does not give individual results for each metal chloride."]},
    {"page":3,"state":"attempted","regions":[
        {"region":"Equation 1","finding":"Four cross-dimerization products 3a-3d recorded independently; collective 6:1-15:1 cross:homo range not assigned to individual products.","evidence_ids":[EV[3],ET[3]]},
        {"region":"Scheme 2","finding":"Seven derivatization experiments recorded: A-E and the two R variants in F.","evidence_ids":[EV[3],ET[3]]},
        {"region":"Figure 1, Equation 2 discussion, and EPR prose","finding":"EPR measurements and the two preformed-complex controls recorded as mechanistic studies; plotted spectra also retained as source context.","evidence_ids":[EV[3],ET[3]]}],"record_ids":by_page[3],"unresolved":["EPR acquisition settings and most workup/purification details are in unavailable SI."]},
    {"page":4,"state":"attempted","regions":[
        {"region":"Equations 2-3","finding":"Preformed-complex controls and DBpin labeling experiment linked across the page break; 74% yield and 95% D at the Z-alkenyl carbon retained separately.","evidence_ids":[EV[4],ET[4]]},
        {"region":"Scheme 3A-D","finding":"Five solvent-screen reactions, two EPR solvent-effect measurements, and one 50 mol% THF rescue reaction recorded independently.","evidence_ids":[EV[4],ET[4]]},
        {"region":"kinetics prose and Scheme 4","finding":"One detail-limited kinetic study and one proposed catalytic cycle; proposed intermediates are not experiments.","evidence_ids":[EV[4],ET[4]]}],"record_ids":by_page[4],"unresolved":["Kinetic plots/protocol are only in unavailable SI Figures S8-S11."]},
    {"page":5,"state":"attempted","regions":[
        {"region":"mechanistic conclusion, summary, acknowledgements, and references","finding":"Continuation of kinetic interpretation and proposed turnover-limiting step; no new synthesis arrow or table row.","evidence_ids":[EV[5],ET[5]]}],"record_ids":by_page[5],"unresolved":[]},
    {"page":6,"state":"attempted","regions":[
        {"region":"references 16-19 and manuscript dates","finding":"Bibliography and publication metadata only; no current-paper reaction arrow, table row, procedure, or standalone spectral peak table.","evidence_ids":[EV[6],ET[6]]}],"record_ids":by_page[6],"unresolved":[]},
]}

raw_inventory = {
    "document_id":DOC,"pages":[1,2,3,4,5,6],"source":"direct visual inspection of rendered original page images",
    "counts":{"molecules":len(M),"reactions":len(R)},
    "decisions":[
        "Table 1 entry 11 was split into NaH, Mg, and Zn because prose calls them reactions; their shared cells remain group-level evidence.",
        "Table 1 entry 12 remains one unresolved grouped screen because individual Mn/Fe/Co/Ni/Cu results cannot be separated.",
        "The Table 2 1p experiment and prose gram-scale 1p experiment are different scales and remain separate.",
        "Scheme 2F R=H and R=OMe are separate experiments.",
        "EPR spectra are mechanistic-study outcomes, not synthetic yields; plotted spectra are also source context.",
        "Proposed Scheme 4 intermediates are not expanded into concrete experiments.",
    ]}

stereo_checks = []
for m in M.values():
    smiles = m["smiles"]["normalized_value"]
    if smiles:
        try:
            stereo_checks.append({"molecule_id":m["molecule_id"],"smiles":smiles,
                                  "standard_stereo":st.stereo_inspect(smiles)})
        except Exception as exc:
            stereo_checks.append({"molecule_id":m["molecule_id"],"smiles":smiles,"error":str(exc)})

(RUN/"submitted.json").write_text(json.dumps(bundle,ensure_ascii=False,indent=2)+"\n")
(RUN/"coverage/pages_0001_0006.json").write_text(json.dumps(coverage,ensure_ascii=False,indent=2)+"\n")
(RUN/"raw_model/page_extraction.json").write_text(json.dumps(raw_inventory,ensure_ascii=False,indent=2)+"\n")
(RUN/"raw_model/stereo_checks.json").write_text(json.dumps(stereo_checks,ensure_ascii=False,indent=2)+"\n")
print(json.dumps({"molecules":len(M),"reactions":len(R),"pages":list(by_page),"records_by_page":{str(k):len(v) for k,v in by_page.items()}},ensure_ascii=False,indent=2))
