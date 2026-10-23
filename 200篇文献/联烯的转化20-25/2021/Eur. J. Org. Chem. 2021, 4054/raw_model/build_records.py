"""Build source-derived v3.1 records for doc_1cf04303b61f633d, PDF pages 1-5."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
import standard_records as st

RUN = Path(__file__).resolve().parents[1]
SUBMITTED = RUN / "submitted.json"
DOC = "doc_1cf04303b61f633d"
E = {n: f"{DOC}:page{n}:image" for n in range(1, 6)}


def F(raw=None, normalized=None, ev=(), method="original_extraction", status=None, derivation=None):
    return st.fact(raw, normalized, ev, method=method, status=status, derivation=derivation)


def derived(normalized, ev, derivation, raw=None):
    return F(raw, normalized, ev, method="structure_derivation", derivation=derivation)


def molecule(mid, label, name, ev, *, smiles=None, chirality="unresolved", geometry="unresolved",
             representation=None, mixture=None, note=None):
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, ev)
    m["names"]["original"] = F(name, name, ev)
    m["names"]["english"] = F(name, name, ev)
    m["names"]["abbreviations"] = F()
    m["names"]["aliases"] = F()
    m["cas_numbers"] = F()
    m["identity_note"] = F(note, note, ev) if note else F()
    if smiles:
        diag = st.stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = F({"source": label, "candidate_transcription": smiles}, smiles, ev)
        m["stereochemistry"]["isomeric_smiles"] = derived(
            diag["canonical_isomeric_smiles"], ev,
            "RDKit canonicalization of the source-derived candidate SMILES; graph parsing is not source verification.",
        )
        m["stereochemistry"]["identity_key"] = derived(
            diag["stereo_graph_key"], ev,
            "RDKit stereo graph key for diagnostics only; it is not a sample identity.",
        )
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = F()
        m["stereochemistry"]["isomeric_smiles"] = F()
        m["stereochemistry"]["identity_key"] = F()
    m["generic_template"] = None
    m["stereochemistry"]["chirality_status"] = derived(
        chirality, ev, "Stereo status assessed from the supplied source drawing; absence of a wedge was not treated as racemic."
    )
    m["stereochemistry"]["geometry_status"] = derived(
        geometry, ev, "Geometry status assessed from the supplied source drawing and substituent identities."
    )
    rep = representation or ("not_applicable" if smiles and chirality == "achiral" and geometry == "not_applicable" else "partial" if smiles else "not_encoded")
    m["stereochemistry"]["representation_status"] = derived(
        rep, ev, "Reports whether source-visible stereochemistry is encoded by the candidate representation."
    )
    m["stereochemistry"]["atom_centers"] = F()
    m["stereochemistry"]["double_bonds"] = F()
    m["stereochemistry"]["axial"] = F()
    m["stereochemistry"]["assignment"] = F()
    m["stereochemistry"]["mixture"] = F(mixture, mixture, ev) if mixture else F()
    return m


def participant(pid, mid, role, ev, *, stage="reaction", amount=None, mass=None, volume=None,
                equiv=None, molpct=None, ref="m_1a", product_role=None):
    p = st.new_participant(pid, mid, role, ev, stage)
    if product_role:
        p["product_role"] = F(product_role, product_role, ev)
    values = {"amount": amount, "mass": mass, "volume": volume, "equivalents": equiv, "mol_percent": molpct}
    for key, val in values.items():
        if val is not None:
            p["charges"][key] = F(val[1], val[0], ev)
    if equiv is not None or molpct is not None:
        p["equivalence_reference"] = F(ref, ref, ev)
    return p


def step(sid, ev, pids, operation, *, temperature=None, time=None, atmosphere=None, stage="reaction"):
    s = st.new_step(sid, ev, stage)
    s["participant_ids"] = list(pids)
    s["operation"] = F(operation, operation, ev)
    if temperature is not None:
        s["temperature"] = F(temperature[1], temperature[0], ev)
    if time is not None:
        s["time"] = F(time[1], time[0], ev)
    if atmosphere is not None:
        s["atmosphere"] = F(atmosphere, atmosphere, ev)
    return s


LABEL_NUM = {"C-H": 1, "C-C": 2, "C-O": 3, "C-N": 5, "C-S": 4, "C-B": 11,
             "C-Si": 12, "DB-C-C": 13, "O-H": 24}


def classification(r, ev, specs=None, coverage="unresolved", limitations=()):
    c = r["classification"]
    c["mapped_reaction_smiles"] = F()
    c["bond_changes"] = []
    labels = []
    if specs:
        supports = {}
        for idx, (a, b, before, after, kind, code) in enumerate(specs, 1):
            cid = f"b{idx}"
            change = {"change_id": cid, "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
                      "before_order": before, "after_order": after, "change_type": kind}
            if kind.startswith("hydrogen_"):
                change["hydrogen_count"] = 1
            c["bond_changes"].append(derived(change, ev, "Net source-drawing bond change; atom mapping is incomplete."))
            if code:
                supports.setdefault(code, []).append(cid)
        labels = [{"category_number": LABEL_NUM[code], "code": code, "supporting_change_ids": ids}
                  for code, ids in supports.items()]
        c["labels"] = derived(labels, ev, "Taxonomy labels derived from the retained net bond changes.")
    else:
        c["labels"] = F() if coverage == "unresolved" else derived([], ev, "No observed transformation to classify.")
    c["coverage"] = derived(coverage, ev, "Classification is limited to source-visible net changes; complete atom mapping is not asserted.")
    c["limitations"] = list(limitations) or ["No complete atom mapping; spectator/coproduct and exact hydrogen-source accounting remain unresolved."]
    c["model_votes"] = []


BORYL = [("C", "B", 0, 1, "formed", "C-B"), ("C", "H", 0, 1, "hydrogen_gain", "C-H")]
SILYL = [("C", "Si", 0, 1, "formed", "C-Si"), ("C", "H", 0, 1, "hydrogen_gain", "C-H")]
OXIDIZED = [("C", "O", 0, 1, "formed", "C-O"), ("C", "H", 0, 1, "hydrogen_gain", "C-H")]


def outcome(oid, kind, targets, ev, value=None, raw=None, *, basis=None, yield_type=None, scope=None):
    o = st.new_outcome(oid, kind, targets, ev)
    if value is not None:
        o["value"] = F(raw, value, ev)
    if basis is not None:
        o["measurement_basis"] = F(basis, basis, ev)
    if yield_type is not None:
        o["yield_type"] = F(yield_type, yield_type, ev)
    if scope is not None:
        o["scope"] = F(scope, scope, ev)
    if targets:
        o["product_role"] = F("reported target", "reported target", ev)
    return o


def reaction(rid, ev, *, rtype="specific_experiment", categories=None, location=None,
             name=None, transform=None, general=True):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = rtype
    r["experiment_categories"] = F(categories, categories, ev) if categories is not None else F()
    r["location"] = F(location, location, ev) if location else F()
    r["reaction_name"] = F(name, name, ev) if name else F()
    if transform:
        r["transformation_type"] = derived(transform, ev, "Transformation description inferred from source reactant/product drawings.")
    r["general_scheme_ids"] = [f"{DOC}:scheme1g_general"] if general and rtype == "specific_experiment" else []
    r["review"]["unresolved"] = []
    return r


bundle = json.loads(SUBMITTED.read_text())

# Registry: complex scope structures are intentionally retained by source label rather than guessed SMILES.
molecules = []
scope_labels = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j"]
molecules.append(molecule("m_1a", "1a", "1,1-difluoroallene 1a", [E[2], E[3]],
                          smiles="CC(CCc1ccccc1)=C=C(F)F", chirality="achiral", geometry="not_applicable"))
for x in scope_labels[1:]:
    chir = "unspecified" if x in ("i", "j") else "achiral"
    molecules.append(molecule(f"m_1{x}", f"1{x}", f"gem-difluoroallene 1{x}", [E[3]],
                              chirality=chir, geometry="not_applicable", representation="not_encoded",
                              note="Source-depicted substrate; full candidate SMILES was not safely transcribed from the dense scope table."))

known = [
    ("m_2a", "2a", "bis(pinacolato)diboron", "CC1(C)OB(B2OC(C)(C)C(C)(C)O2)OC1(C)C", [E[2], E[3]]),
    ("m_3", "3", "dimethylphenylsilyl pinacol boronate (PhMe2Si-Bpin)", "C[Si](C)(c1ccccc1)B1OC(C)(C)C(C)(C)O1", [E[2], E[3]]),
    ("m_CuI", "CuI", "copper(I) iodide", "[Cu]I", [E[2], E[3]]),
    ("m_KOtBu", "KOtBu", "potassium tert-butoxide", "CC(C)(C)[O-].[K+]", [E[2], E[3]]),
    ("m_tBuOH", "tBuOH", "tert-butanol", "CC(C)(C)O", [E[2], E[3]]),
    ("m_THF", "THF", "tetrahydrofuran", "C1CCOC1", [E[2], E[3]]),
    ("m_PPh3", "PPh3", "triphenylphosphine", "P(c1ccccc1)(c1ccccc1)c1ccccc1", [E[2]]),
    ("m_NaH", "NaH", "sodium hydride", "[Na+].[H-]", [E[3]]),
    ("m_DMSO", "DMSO", "dimethyl sulfoxide", "CS(C)=O", [E[3]]),
    ("m_ClCH2Br", "ClCH2Br", "bromochloromethane", "ClCBr", [E[3]]),
    ("m_benzaldehyde", "PhCHO", "benzaldehyde", "O=Cc1ccccc1", [E[3]]),
    ("m_TFAA", "TFAA", "trifluoroacetic anhydride", "O=C(OC(=O)C(F)(F)F)C(F)(F)F", [E[3]]),
    ("m_NaOH", "NaOH(aq.)", "aqueous sodium hydroxide", "[Na+].[OH-]", [E[3]]),
    ("m_H2O2", "H2O2", "hydrogen peroxide", "OO", [E[3]]),
]
for mid, label, name, smi, ev in known:
    molecules.append(molecule(mid, label, name, ev, smiles=smi, chirality="achiral", geometry="not_applicable"))

for mid, label, name, ev in [
    ("m_L1", "L1", "N-heterocyclic carbene precursor L1", [E[2]]),
    ("m_L2", "L2", "bis(oxazoline)-type ligand L2", [E[2]]),
    ("m_L3", "L3", "bisphosphine ligand L3", [E[2]]),
    ("m_BpinBdan", "Bpin-Bdan", "pinacolato/1,8-diaminonaphthalene diboron reagent", [E[2], E[3]]),
    ("m_diboron_2", "Table 3 diboron reagent 2 (middle)", "source-depicted diboron reagent", [E[3]]),
    ("m_diboron_3", "Table 3 diboron reagent 2 (right)", "source-depicted diboron reagent", [E[3]]),
    ("m_nBuLi", "nBuLi", "n-butyllithium", [E[3]]),
]:
    molecules.append(molecule(mid, label, name, ev, chirality="unresolved", geometry="not_applicable",
                              note="Identity retained from the source label/drawing; no complete SMILES asserted."))

# Main products; only 4a/5a are transcribed as candidate graphs.
molecules.append(molecule("m_4a", "4a", "gem-difluoroallylboronate 4a", [E[2], E[3]],
                          smiles="CC(CCc1ccccc1)(B1OC(C)(C)C(C)(C)O1)C=C(F)F",
                          chirality="unspecified", geometry="not_applicable", representation="partial"))
molecules.append(molecule("m_5a", "5a", "gem-difluoroallylsilane 5a", [E[2], E[3]],
                          smiles="CC(CCc1ccccc1)([Si](C)(C)c1ccccc1)C=C(F)F",
                          chirality="unspecified", geometry="not_applicable", representation="partial"))
dr_values = {"4d": "dr = 93:7", "4i": "dr = 57:43", "4j": "dr = 55:45",
             "5i": "dr = 55:45", "5j": "dr = 64:36"}
for label in ["4b", "4c", "4d", "4e", "4f", "4g", "4i", "4j", "6",
              "5b", "5c", "5e", "5f", "5h", "5i", "5j"]:
    mix = dr_values.get(label)
    molecules.append(molecule(f"m_{label}", label, f"source-depicted product {label}", [E[3]],
                              chirality="mixture" if mix else "unspecified", geometry="not_applicable",
                              representation="not_encoded", mixture=mix,
                              note="Full source structure retained in the page image; no guessed complete SMILES."))
for label, chir, geom, mix in [
    ("8", "unspecified", "not_applicable", None),
    ("9", "unspecified", "not_applicable", None),
    ("10", "unspecified", "not_applicable", None),
    ("11", "unspecified", "mixture", "Z/E = 34:66"),
    ("12", "mixture", "mixture", "E/Z > 95:5; dr = 60:40"),
]:
    molecules.append(molecule(f"m_{label}", label, f"source-depicted product {label}", [E[3]],
                              chirality=chir, geometry=geom, representation="not_encoded", mixture=mix,
                              note="Full source structure retained in Scheme 2/Table 3; no unique stereochemical SMILES asserted."))

bundle["molecules"] = molecules
reactions = []

# Page 1: six cited-work panels, plus the current-work general reaction design.
prior = {
    "a": "Prior copper-catalyzed alpha-addition of non-fluorinated allenes",
    "b": "Prior beta-addition/protonation or carbonyl capture of non-fluorinated allenes",
    "c": "Prior beta-addition/electrophile capture of non-fluorinated allenes",
    "d": "Hammond 2006 Grignard addition-elimination of gem-difluoroallenes",
    "e": "Ichikawa 2012 gamma-selective organocopper addition/electrophile capture",
    "f": "Ichikawa 2020 gold-catalyzed alpha/gamma C-N or C-S bond formation",
}
for panel, desc in prior.items():
    r = reaction(f"scheme1{panel}_literature", [E[1]], rtype="literature_example", categories=None,
                 location=f"PDF page 1, Scheme 1{panel}", name=desc, general=False)
    r["notes"] = [F({"figure_panel": f"Scheme 1{panel}", "status": "cited prior work"}, desc, [E[1]])]
    if panel == "f":
        classification(r, [E[1]], [("C", "N", 0, 1, "formed", "C-N"), ("C", "S", 0, 1, "formed", "C-S")], "partial",
                       ["Two alternative literature branches are shown; no single experiment is asserted."])
    elif panel == "d":
        classification(r, [E[1]], [("C", "C", 0, 1, "formed", "C-C")], "partial")
    else:
        classification(r, [E[1]], None, "unresolved")
    r["review"]["unresolved"] = ["Cited prior-work panel is not a present-paper experiment; detailed conditions and identities were not expanded."]
    reactions.append(r)

r = reaction("scheme1g_general", [E[1]], rtype="generic_scheme", categories=None,
             location="PDF page 1, Scheme 1g", name="Copper-catalyzed protoborylation/protosilylation of gem-difluoroallenes", general=False)
classification(r, [E[1]], BORYL + SILYL, "partial", ["The figure shows borylation and silylation as alternative generic branches; atom mapping is not complete."])
r["notes"] = [F("This work: copper catalyzed proto-borylation/silylation of gem-difluoroallene (gamma-addition)",
                  "Current-work generic reaction design", [E[1]])]
r["review"]["unresolved"] = ["Generic R1/R2 structures are not concrete molecules; no experiment-specific yield belongs to this record."]
reactions.append(r)


def add_common_reaction_fields(r, ev, solvent="THF", volume_ml=1.0, procedure=None):
    r["solvent_conditions"]["composition"] = F(solvent, solvent, ev)
    if volume_ml is not None:
        r["solvent_conditions"]["total_volume"] = F(f"{volume_ml:g} mL", {"number": volume_ml, "unit": "mL"}, ev)
    if procedure:
        r["procedure"]["original_text"] = F(procedure, procedure, ev)


# Page 2, Table 1: one independent record per optimization row.
table1 = [
    (1, 10, 2.0, 1.0, 28, 72, 65), (2, 10, 2.0, 1.0, 40, 84, None),
    (3, 10, 2.0, 1.0, 50, 86, None), (4, 10, 2.0, 1.0, 60, 15, None),
    (5, 10, 2.0, 2.0, 50, 87, None), (6, 10, 2.0, 3.0, 50, 85, None),
    (7, 10, 1.5, 2.0, 50, 85, None), (8, 10, 1.2, 2.0, 50, 86, None),
    (9, 5, 1.2, 2.0, 50, 88, 86), (10, 5, 1.2, 2.0, 50, 87, None),
    (11, None, 1.2, 2.0, 50, "N.R.", None), (12, 5, 1.2, 2.0, 50, "N.R.", None),
    (13, 10, 2.0, 1.0, 28, 93, 89), (14, 10, 2.0, 1.0, 28, 54, None),
    (15, 10, 2.0, 1.0, 28, 14, None), (16, 10, 2.0, 1.0, 28, 0, None),
    (17, 10, 1.5, 1.0, 28, 94, None), (18, 5, 1.5, 1.0, 28, 91, 86),
]
table1_proc = "Reaction conditions: 1a (0.2 mmol), 2, Cu catalyst, tBuOH, KOtBu (11 mol%) and solvent (1.0 mL) in argon, temperature."
for entry, cu, req, oh, temp, nmr_yield, isolated in table1:
    silyl = entry >= 13
    ev = [E[2], E[3]] if entry in (9, 18) else [E[2]]
    cats = ["optimization", "substrate_scope"] if entry in (9, 18) else ["optimization"]
    loc = f"PDF page 2, Table 1 entry {entry}"
    if entry == 9:
        loc += "; repeated as Table 2 product 4a on PDF page 3"
    if entry == 18:
        loc += "; repeated as Table 2 product 5a on PDF page 3"
    r = reaction(f"table1_entry_{entry:02d}", ev, categories=cats, location=loc,
                 name="protosilylation" if silyl else "protoborylation",
                 transform="gamma-protosilylation of gem-difluoroallene" if silyl else "gamma-protoborylation of gem-difluoroallene")
    reagent = "m_3" if silyl else "m_2a"
    product = "m_5a" if silyl else "m_4a"
    p = [participant("p_substrate", "m_1a", "reactant", ev, amount=({"number": 0.2, "unit": "mmol"}, "0.2 mmol")),
         participant("p_reagent", reagent, "reagent", ev, equiv=({"number": req, "unit": "equiv"}, f"{req} equiv")),
         participant("p_tBuOH", "m_tBuOH", "additive", ev, equiv=({"number": oh, "unit": "equiv"}, f"{oh} equiv")),
         participant("p_THF", "m_THF", "solvent", ev, volume=({"number": 1.0, "unit": "mL"}, "1.0 mL"))]
    if cu is not None:
        p.append(participant("p_CuI", "m_CuI", "catalyst", ev, molpct=({"number": cu, "unit": "mol%"}, f"{cu} mol%")))
    if entry != 12:
        base_pct = 23 if entry == 14 else 11
        p.append(participant("p_KOtBu", "m_KOtBu", "reagent", ev, molpct=({"number": base_pct, "unit": "mol%"}, f"{base_pct} mol%")))
    ligand = {10: ("m_PPh3", 10), 14: ("m_L1", 12), 15: ("m_L2", 12), 16: ("m_L3", 12)}.get(entry)
    if ligand:
        p.append(participant("p_ligand", ligand[0], "ligand", ev, molpct=({"number": ligand[1], "unit": "mol%"}, f"{ligand[1]} mol%")))
    if nmr_yield != "N.R.":
        p.append(participant("p_product", product, "product", ev, product_role="target product"))
    r["participants"] = p
    r["reaction_scale"] = F("1a (0.2 mmol)", {"reference_molecule_id": "m_1a", "amount": {"number": 0.2, "unit": "mmol"}}, ev)
    add_common_reaction_fields(r, ev, procedure=table1_proc)
    r["steps"] = [step("reaction", ev, [x["participant_id"] for x in p if x["stage"] == "reaction"],
                       "combine the listed reaction components and react", temperature=({"number": temp, "unit": "°C"}, f"{temp} °C"),
                       time=({"number": 12 if silyl else 24, "unit": "h"}, f"{12 if silyl else 24} h"), atmosphere="argon")]
    if nmr_yield == "N.R.":
        r["outcomes"] = [outcome("no_reaction", "no_reaction", [], ev, True, "N.R.")]
        classification(r, ev, None, "not_applicable", ["No observed product; no intended bond labels assigned."])
    else:
        r["outcomes"] = [outcome("nmr_yield", "yield", [product], ev,
                                  {"number": nmr_yield, "unit": "%"}, f"{nmr_yield}%",
                                  basis="1H NMR analysis of crude reaction mixture", yield_type="NMR yield")]
        if isolated is not None:
            r["outcomes"].append(outcome("isolated_yield", "yield", [product], ev,
                                                 {"number": isolated, "unit": "%"}, f"({isolated})",
                                                 basis="isolated", yield_type="isolated yield"))
        classification(r, ev, SILYL if silyl and nmr_yield > 0 else BORYL if not silyl and nmr_yield > 0 else None,
                       "partial" if nmr_yield > 0 else "not_applicable",
                       ["Zero-yield entry has no observed bond formation."] if nmr_yield == 0 else ())
    if entry in (11, 12, 16):
        r["review"]["unresolved"] = ["No observed target product; source gives no detection limit for N.R./0 yield."]
    reactions.append(r)

# Page 3, Table 2: two repeated optimum examples are already linked above; add the remaining 16 experiments.
scope = [
    ("4b", "1b", "m_2a", 96, 36, 5, 50, 2.0, None),
    ("4c", "1c", "m_2a", 95, 24, 5, 50, 2.0, None),
    ("4d", "1d", "m_2a", 63, 72, 10, 50, 2.0, "dr = 93:7"),
    ("4e", "1e", "m_2a", 85, 24, 5, 50, 2.0, None),
    ("4f", "1f", "m_2a", 88, 48, 5, 50, 2.0, None),
    ("4g", "1g", "m_2a", 6, 72, 5, 75, 2.0, None),
    ("4i", "1i", "m_2a", 77, 36, 5, 50, 2.0, "dr = 57:43"),
    ("4j", "1j", "m_2a", 92, 30, 5, 50, 2.0, "dr = 55:45"),
    ("6", "1a", "m_BpinBdan", 27, 24, 5, 50, 2.0, None),
    ("5b", "1b", "m_3", 59, 36, 5, 28, 1.0, None),
    ("5c", "1c", "m_3", 82, 30, 5, 28, 1.0, None),
    ("5e", "1e", "m_3", 46, 30, 5, 28, 1.0, None),
    ("5f", "1f", "m_3", 85, 12, 10, 28, 1.0, None),
    ("5h", "1h", "m_3", 86, 12, 5, 28, 1.0, None),
    ("5i", "1i", "m_3", 87, 12, 5, 28, 1.0, "dr = 55:45"),
    ("5j", "1j", "m_3", 93, 12, 5, 28, 1.0, "dr = 64:36"),
]
table2_proc = "Reaction conditions: 1 (0.2 mmol), B2pin2 (0.24 mmol), CuI (5 mol%), KOtBu (11 mol%), tBuOH (2.0 equiv) and THF (1 mL) in argon, 50 °C for 24 h; product-specific footnotes modify these conditions."
for product, substrate, reagent, yld, hours, cu, temp, oh, dr in scope:
    silyl = product.startswith("5")
    r = reaction(f"table2_{product}", [E[3]], categories=["substrate_scope"],
                 location=f"PDF page 3, Table 2 product {product}",
                 name="protosilylation" if silyl else "protoborylation",
                 transform="gamma-protosilylation of gem-difluoroallene" if silyl else "gamma-protoborylation of gem-difluoroallene")
    product_mid = f"m_{product}"
    req_equiv = 1.5 if silyl else 1.2
    p = [participant("p_substrate", f"m_{substrate}", "reactant", [E[3]], amount=({"number": 0.2, "unit": "mmol"}, "0.2 mmol")),
         participant("p_reagent", reagent, "reagent", [E[3]], equiv=({"number": req_equiv, "unit": "equiv"}, f"{req_equiv} equiv"), ref=f"m_{substrate}"),
         participant("p_CuI", "m_CuI", "catalyst", [E[3]], molpct=({"number": cu, "unit": "mol%"}, f"{cu} mol%"), ref=f"m_{substrate}"),
         participant("p_KOtBu", "m_KOtBu", "reagent", [E[3]], molpct=({"number": 11, "unit": "mol%"}, "11 mol%"), ref=f"m_{substrate}"),
         participant("p_tBuOH", "m_tBuOH", "additive", [E[3]], equiv=({"number": oh, "unit": "equiv"}, f"{oh} equiv"), ref=f"m_{substrate}"),
         participant("p_THF", "m_THF", "solvent", [E[3]], volume=({"number": 1, "unit": "mL"}, "1 mL")),
         participant("p_product", product_mid, "product", [E[3]], product_role="target product")]
    r["participants"] = p
    r["reaction_scale"] = F(f"{substrate} (0.2 mmol)", {"reference_molecule_id": f"m_{substrate}", "amount": {"number": 0.2, "unit": "mmol"}}, [E[3]])
    add_common_reaction_fields(r, [E[3]], procedure=table2_proc)
    r["steps"] = [step("reaction", [E[3]], [x["participant_id"] for x in p], "react under the Table 2 conditions",
                       temperature=({"number": temp, "unit": "°C"}, f"{temp} °C"),
                       time=({"number": hours, "unit": "h"}, f"{hours} h"), atmosphere="argon")]
    r["outcomes"] = [outcome("yield", "yield", [product_mid], [E[3]], {"number": yld, "unit": "%"}, f"{yld}%")]
    if dr:
        a, b = [int(x) for x in dr.split("=")[1].split(":")]
        r["outcomes"].append(outcome("dr", "selectivity", [product_mid], [E[3]],
                                             {"ratio": [a, b], "raw_text": dr}, dr,
                                             basis="1H NMR analysis of crude reaction mixture", scope="diastereomeric ratio; major/minor configurations unassigned"))
    classification(r, [E[3]], SILYL if silyl else BORYL, "partial")
    r["review"]["unresolved"] = ["Yield measurement basis is not stated for this Table 2 entry.",
                                     "Dense scope-table structure retained by label/image; complete SMILES not asserted."]
    reactions.append(r)

# Page 3, Table 3: three independent two-stage experiments ending in alcohol 8.
for idx, (reagent, yld) in enumerate([("m_2a", 71), ("m_diboron_2", 43), ("m_diboron_3", 55)], 1):
    r = reaction(f"table3_diboron_{idx}", [E[3]], categories=["substrate_scope"],
                 location=f"PDF page 3, Table 3 diboron reagent {idx}", name="borylation followed by oxidation",
                 transform="two-stage conversion of gem-difluoroallene to gem-difluoroallylic alcohol")
    p = [participant("p_substrate", "m_1a", "reactant", [E[3]], amount=({"number": 0.2, "unit": "mmol"}, "0.2 mmol")),
         participant("p_diboron", reagent, "reagent", [E[3]], amount=({"number": 0.24, "unit": "mmol"}, "0.24 mmol")),
         participant("p_CuI", "m_CuI", "catalyst", [E[3]], molpct=({"minimum": 5, "maximum": 10, "unit": "mol%", "raw_text": "5-10 mol%"}, "5-10 mol%")),
         participant("p_KOtBu", "m_KOtBu", "reagent", [E[3]], molpct=({"number": 11, "unit": "mol%"}, "11 mol%")),
         participant("p_tBuOH", "m_tBuOH", "additive", [E[3]], equiv=({"number": 2, "unit": "equiv"}, "2.0 equiv")),
         participant("p_THF", "m_THF", "solvent", [E[3]], volume=({"number": 1, "unit": "mL"}, "1 mL")),
         participant("p_NaOH", "m_NaOH", "workup_agent", [E[3]], stage="workup"),
         participant("p_H2O2", "m_H2O2", "workup_agent", [E[3]], stage="workup"),
         participant("p_product", "m_8", "product", [E[3]], product_role="oxidized target product")]
    r["participants"] = p
    r["reaction_scale"] = F("1a (0.2 mmol)", {"reference_molecule_id": "m_1a", "amount": {"number": 0.2, "unit": "mmol"}}, [E[3]])
    add_common_reaction_fields(r, [E[3]], procedure="1 (0.2 mmol), 2 (0.24 mmol), CuI (5-10 mol%), KOtBu (11 mol%), tBuOH (2.0 equiv) and THF (1 mL) in argon, 50 °C for 72 h, then NaOH(aq.), H2O2.")
    r["steps"] = [step("borylation", [E[3]], [x["participant_id"] for x in p if x["stage"] == "reaction"], "form unstable intermediate 7",
                       temperature=({"number": 50, "unit": "°C"}, "50 °C"), time=({"number": 72, "unit": "h"}, "72 h"), atmosphere="argon"),
                  step("oxidation", [E[3]], ["p_NaOH", "p_H2O2", "p_product"], "oxidize with aqueous NaOH and H2O2", temperature=("r.t.", "r.t."), stage="workup")]
    r["outcomes"] = [outcome("yield", "yield", ["m_8"], [E[3]], {"number": yld, "unit": "%"}, f"{yld}%")]
    classification(r, [E[3]], OXIDIZED, "partial", ["Unstable intermediate 7 is not independently isolated; labels describe the net two-stage product 8."])
    r["review"]["unresolved"] = ["Specific CuI loading within the reported 5-10 mol% range is not assigned per diboron reagent.",
                                     "Oxidation time and yield measurement basis are unreported."]
    reactions.append(r)

# Scheme 2 scale-up and derivatizations.
for kind, reagent, product, yld, mass_g, temp, hours, eq, cls in [
    ("boryl", "m_2a", "m_4a", 86, 1.669, 50, 24, 1.2, BORYL),
    ("silyl", "m_3", "m_5a", 90, 1.784, 28, 12, 1.5, SILYL),
]:
    ev = [E[3]]
    r = reaction(f"scheme2_scaleup_{kind}", ev, categories=["scale_up"],
                 location=f"PDF page 3, Scheme 2 eq. {'1' if kind == 'boryl' else '2'}",
                 name=f"gram-scale proto{kind}ation", transform=f"gram-scale gamma-proto{kind}ation")
    p = [participant("p_substrate", "m_1a", "reactant", ev, amount=({"number": 6, "unit": "mmol"}, "6 mmol")),
         participant("p_reagent", reagent, "reagent", ev, equiv=({"number": eq, "unit": "equiv"}, f"{eq} equiv")),
         participant("p_CuI", "m_CuI", "catalyst", ev, molpct=({"number": 5, "unit": "mol%"}, "5 mol%")),
         participant("p_KOtBu", "m_KOtBu", "reagent", ev, molpct=({"number": 11, "unit": "mol%"}, "11 mol%")),
         participant("p_tBuOH", "m_tBuOH", "additive", ev, equiv=({"number": 2 if kind == 'boryl' else 1, "unit": "equiv"}, "2.0 equiv" if kind == 'boryl' else "1.0 equiv")),
         participant("p_THF", "m_THF", "solvent", ev),
         participant("p_product", product, "product", ev, product_role="target product")]
    r["participants"] = p
    r["reaction_scale"] = F("1a, 6 mmol", {"reference_molecule_id": "m_1a", "amount": {"number": 6, "unit": "mmol"}}, ev)
    add_common_reaction_fields(r, ev, volume_ml=None)
    r["steps"] = [step("reaction", ev, [x["participant_id"] for x in p], "gram-scale reaction",
                       temperature=({"number": temp, "unit": "°C"}, f"{temp} °C"), time=({"number": hours, "unit": "h"}, f"{hours} h"))]
    r["outcomes"] = [outcome("yield", "yield", [product], ev, {"number": yld, "unit": "%"}, f"{yld}%"),
                     outcome("product_amount", "product_amount", [product], ev, {"number": mass_g, "unit": "g"}, f"{mass_g} g")]
    classification(r, ev, cls, "partial")
    r["review"]["unresolved"] = ["Scheme 2 does not state solvent volume, atmosphere, workup, purification, or yield measurement basis."]
    reactions.append(r)

r = reaction("scheme2_homologation", [E[3]], categories=["derivatization"], location="PDF page 3, Scheme 2 eq. 3",
             name="homologation with a halomethyl lithium reagent", transform="homologation and competing oxygen migration")
p = [participant("p_4a", "m_4a", "reactant", [E[3]]),
     participant("p_ClCH2Br", "m_ClCH2Br", "reagent", [E[3]], equiv=({"number": 2, "unit": "equiv"}, "2.0 equiv"), ref="m_4a"),
     participant("p_nBuLi", "m_nBuLi", "reagent", [E[3]], equiv=({"number": 1.5, "unit": "equiv"}, "1.5 equiv"), ref="m_4a"),
     participant("p_THF", "m_THF", "solvent", [E[3]]),
     participant("p_9", "m_9", "product", [E[3]], product_role="homoallylboronate product"),
     participant("p_10", "m_10", "product", [E[3]], product_role="oxygen migration product")]
r["participants"] = p
add_common_reaction_fields(r, [E[3]], volume_ml=None)
r["steps"] = [step("low_temperature_reaction", [E[3]], [x["participant_id"] for x in p], "react at -78 °C before warming",
                   temperature=({"number": -78, "unit": "°C"}, "-78 °C"), atmosphere="Ar"),
              step("warm_to_room_temperature", [E[3]], [x["participant_id"] for x in p], "warm from -78 °C to room temperature; source prints 2 h for the program",
                   temperature=("r.t.", "r.t."), time=({"number": 2, "unit": "h"}, "2 h"), atmosphere="Ar")]
r["outcomes"] = [outcome("yield_9", "yield", ["m_9"], [E[3]], {"number": 31, "unit": "%"}, "31%"),
                 outcome("yield_10", "yield", ["m_10"], [E[3]], {"number": 28, "unit": "%"}, "28%")]
classification(r, [E[3]], [("C", "C", 0, 1, "formed", "C-C"), ("C", "O", 0, 1, "formed", "C-O")], "partial",
               ["C-C homologation and C-O migration labels span two competing products; product-specific complete atom maps are unavailable."])
r["review"]["unresolved"] = ["The printed 2 h may describe the combined -78 °C-to-room-temperature program; ramp/hold allocation is not separately reported.",
                                 "Yield measurement bases and product-specific complete structures remain unreported/unencoded."]
reactions.append(r)

r = reaction("scheme2_5a_to_11", [E[3]], categories=["derivatization"], location="PDF page 3, Scheme 2 eq. 4",
             name="allyl protonation process", transform="desilylation/rearrangement to difluoromethyl alkene")
p = [participant("p_5a", "m_5a", "reactant", [E[3]]),
     participant("p_NaH", "m_NaH", "reagent", [E[3]], equiv=({"number": 3, "unit": "equiv"}, "3.0 equiv"), ref="m_5a"),
     participant("p_DMSO", "m_DMSO", "solvent", [E[3]]),
     participant("p_11", "m_11", "product", [E[3]], product_role="E/Z product mixture")]
r["participants"] = p
r["solvent_conditions"]["composition"] = F("DMSO", "DMSO", [E[3]])
r["steps"] = [step("reaction", [E[3]], [x["participant_id"] for x in p], "react with NaH in DMSO",
                   temperature=({"number": 25, "unit": "°C"}, "25 °C"), time=({"number": 12, "unit": "h"}, "12 h"), atmosphere="Ar")]
r["outcomes"] = [outcome("yield", "yield", ["m_11"], [E[3]], {"number": 80, "unit": "%"}, "80%"),
                 outcome("ZE", "selectivity", ["m_11"], [E[3]], {"ratio": [34, 66], "raw_text": "Z/E = 34:66"}, "Z/E = 34:66",
                         scope="Z/E ratio; individual isomers not separately identified")]
classification(r, [E[3]], [("C", "Si", 1, 0, "broken", None), ("C", "H", 0, 1, "hydrogen_gain", "C-H"),
                             ("C", "C", 1, 2, "order_increased", "DB-C-C")], "partial")
r["review"]["unresolved"] = ["Yield and Z/E measurement bases, workup, purification, and absolute product configurations are unreported."]
reactions.append(r)

r = reaction("scheme2_4i_to_12", [E[3]], categories=["derivatization"], location="PDF page 3, Scheme 2 eq. 5",
             name="reaction of gem-difluoroallylboronate 4i with benzaldehyde", transform="allylboration of benzaldehyde")
p = [participant("p_4i", "m_4i", "reactant", [E[3]]),
     participant("p_nBuLi", "m_nBuLi", "reagent", [E[3]], equiv=({"number": 1.1, "unit": "equiv"}, "1.1 equiv"), ref="m_4i"),
     participant("p_TFAA", "m_TFAA", "reagent", [E[3]], equiv=({"number": 1.2, "unit": "equiv"}, "1.2 equiv"), ref="m_4i"),
     participant("p_PhCHO", "m_benzaldehyde", "reactant", [E[3]], equiv=({"number": 1.5, "unit": "equiv"}, "1.5 equiv"), ref="m_4i"),
     participant("p_THF", "m_THF", "solvent", [E[3]]),
     participant("p_12", "m_12", "product", [E[3]], product_role="allylation product mixture")]
r["participants"] = p
add_common_reaction_fields(r, [E[3]], volume_ml=None)
r["steps"] = [step("nBuLi_stage", [E[3]], ["p_4i", "p_nBuLi", "p_THF"], "treat with nBuLi",
                   temperature=({"number": -78, "unit": "°C"}, "-78 °C"), time=({"number": 15, "unit": "min"}, "15 min")),
              step("TFAA_stage", [E[3]], ["p_4i", "p_TFAA", "p_THF"], "add TFAA",
                   temperature=({"number": -78, "unit": "°C"}, "-78 °C"), time=({"number": 30, "unit": "min"}, "30 min")),
              step("benzaldehyde_stage", [E[3]], ["p_4i", "p_PhCHO", "p_THF", "p_12"], "add PhCHO at -78 °C then warm to room temperature",
                   temperature=({"number": -78, "unit": "°C"}, "-78 °C")),
              step("warm_to_room_temperature", [E[3]], ["p_4i", "p_PhCHO", "p_THF", "p_12"], "warm to room temperature",
                   temperature=("r.t.", "r.t."))]
r["outcomes"] = [outcome("yield", "yield", ["m_12"], [E[3]], {"number": 70, "unit": "%"}, "70%"),
                 outcome("EZ", "selectivity", ["m_12"], [E[3]], {"ratio": [95, 5], "qualifier": "gt", "raw_text": "E/Z > 95:5"}, "E/Z > 95:5",
                         scope="E/Z ratio; major configuration reported as E"),
                 outcome("dr", "selectivity", ["m_12"], [E[3]], {"ratio": [60, 40], "raw_text": "dr = 60:40"}, "dr = 60:40",
                         scope="diastereomeric ratio; individual configurations unassigned")]
classification(r, [E[3]], [("C", "C", 0, 1, "formed", "C-C"), ("O", "H", 0, 1, "hydrogen_gain", "O-H"),
                             ("C", "O", 2, 1, "order_decreased", None)], "partial")
r["review"]["unresolved"] = ["Time for the final -78 °C-to-room-temperature stage is unreported.",
                                 "Yield/selectivity measurement bases, workup, purification, dr isomer ordering, and absolute configurations are unreported."]
reactions.append(r)

# Page 4 proposed mechanism; page 5 contains references only.
r = reaction("scheme3_proposed_mechanism", [E[4]], rtype="proposed_mechanism", categories=None,
             location="PDF page 4, Scheme 3", name="Proposed copper boryl/silyl catalytic cycle", general=False)
r["notes"] = [F("Scheme 3. Proposed reaction mechanism.", "Author-proposed mechanism; not an experimental proof record.", [E[4]])]
classification(r, [E[4]], None, "not_applicable", ["Mechanistic proposal is not an observed standalone transformation experiment."])
r["review"]["unresolved"] = ["Catalytic intermediates I-III and regioselectivity rationale are author proposals, not directly verified experiments in the supplied pages."]
reactions.append(r)

bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": f"{DOC}:article_metadata", "document_id": DOC, "evidence_ids": [E[1]],
     "kind": "article_metadata", "title": "Copper Catalyzed Protosilylation/Protoborylation of gem-Difluoroallenes",
     "authors": ["Cui-Cui Shan", "Kai-Yang Dai", "Meng Zhao", "Yun-He Xu"], "doi": "10.1002/ejoc.202100304"},
    {"context_id": f"{DOC}:page5_references", "document_id": DOC, "evidence_ids": [E[5]],
     "kind": "non_experimental_page", "note": "Page 5 contains references 11-22 and manuscript dates; no current-paper experiment was identified."},
]

SUBMITTED.write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")

inventory = {
    "document_id": DOC,
    "pages_viewed": [1, 2, 3, 4, 5],
    "page_findings": {
        "1": "Title/abstract/introduction; Scheme 1a-f cited prior work and Scheme 1g current generic reaction design; SI mentioned but not supplied.",
        "2": "Table 1 has 18 independent optimization rows with footnotes; prose discusses conditions and scope.",
        "3": "Table 2 scope, Table 3 three diboron/oxidation rows, five Scheme 2 arrows, and beginning of author-proposed mechanism discussion.",
        "4": "Scheme 3 proposed mechanism; conclusion, acknowledgements, references. No new synthesis experiment.",
        "5": "References and manuscript dates only. No new synthesis experiment.",
    },
    "limitations": ["No SI file supplied.", "Dense complex structures without a safe full transcription retain source labels and image evidence."],
}
(RUN / "raw_model/page1_5_visual_inventory.json").write_text(json.dumps(inventory, ensure_ascii=False, indent=2) + "\n")

raw_candidates = {
    "document_id": DOC,
    "status": "machine_candidates_needing_review",
    "molecule_ids": [m["molecule_id"] for m in molecules],
    "reaction_ids": [r["reaction_id"] for r in reactions],
    "notes": [
        "This is a raw candidate index, not an independently reviewed chemistry dataset.",
        "The authoritative source remains the supplied PDF/page images; missing SI was not inferred.",
    ],
}
(RUN / "raw_model/extraction_candidates.json").write_text(json.dumps(raw_candidates, ensure_ascii=False, indent=2) + "\n")

print(json.dumps({"molecules": len(molecules), "reactions": len(reactions)}, ensure_ascii=False))
